
"""
`module Lsq`

This sub-modules contains code to convert a database `LsqDB` into
linear LSQ systems, apply suitable weights, solve the resulting systems
and return an IP (e.g. an NBodyIP).

The ordering of observations (entries of Y, row-indices of Psi) is given by
the loop ordering
```
for ct in configtypes
   for dt in datatypes  # observation types
      for dat in db.data_groups[ct]
      # or: db.kron_groups[ct][dt][:,:,...]
```
E.g., if the configtypes are "E", "F", then
```
Y = [E, E, E, E...., F, F, F, F, ...]
```
"""
module Lsq

import JuLIP, IPFitting, StatsBase

using StaticArrays
using JuLIP: AbstractCalculator, Atoms
using JuLIP.Potentials: OneBody
using JuLIP.MLIPs: SumIP
using IPFitting: Dat, LsqDB, weighthook, observations,
                      observation, hasobservation, eval_obs, vec_obs
using IPFitting.Data: configtype
using IPFitting.DB: dbpath, _nconfigs, matrows

using LinearAlgebra: lmul!, Diagonal, qr, qr!, cond, norm, svd
using InteractiveUtils: versioninfo
using LowRankApprox

const Err = IPFitting.Errors

export lsqfit, onb



"""
`collect_observations(db::LsqDB, weights::Dict, Vref) -> Y, W, Icfg

Construct vectors of observations and weights from the lsq system specified
by `db`. The observations do not have the weights applied yet.

Weights are specified as follows:
```
   wh = weighthook(obskey, cfg)
   w = w1 * wh
```
where `w1` is obtained from the `weights` dictionary.
```
weights = Dict( "default" => Dict( "E" => 30.0, "F" => 1.0, "V" => 0.3 ),
                "bulk" => Dict( "E" => 100.0, "F" = 20.0, "V" => 0.0 ),
                ... )
```
To compute the weights, the code first looks for the `configtype` entry. If
it doesn't exist, then it reverts to the `default` entry. The weights obtained
from the `weights` dictionary is then multiplied by a

- The above procedure can be overwritten by an entry in `cfg.info["W"]`;
in this case the re-weighting via the `weighthook` is ignored as well.
- If `configtype(cfg) ∈ weights["ignore"]` then this configuration is skipped.
- If `obskey ∈ weights["ignore"]` then this observation is skipped
"""
function collect_observations(db::LsqDB,
                              weights::Dict,
                              Vref )

   nrows = size(db.Ψ, 1)  # total number of observations if we collect
                          # everything in the database
   Y = zeros(Float64, nrows)
   W = zeros(Float64, nrows)
   Icfg = zeros(Int, nrows)
   for (obskey, dat, icfg) in observations(db)  # obskey ∈ {"E","F",...}; d::Dat
      if configtype(dat) in weights["ignore"] || obskey in weights["ignore"]
         continue
      end

      irows = matrows(dat, obskey)
      ct = configtype(dat)

      # ----- Observation ------
      obs = observation(obskey, dat)
      # If Vref != nothing the it is a calculator and we can subtract
      # its value => this allows us to fit from a reference potential
      if haskey(dat.info, "Vref")
         obs = obs - dat.info["Vref"][obskey]
      elseif Vref != nothing
         obs = obs - vec_obs(obskey, eval_obs(obskey, Vref, dat))
      end
      Y[irows] .= obs
      Icfg[irows] .= icfg
      # ------ Weights --------
      # modify the weights from extra information in the dat structure
      W[irows] .= _get_weights(weights,
                               weighthook(obskey, dat),
                               dat, obskey, obs)
   end
   return Y, W, Icfg
end

"""
see documentation in `collect_observations`
"""
function _get_weights(weights, wh, dat, obskey, o)
   # initialise the weight to 0.0. If this isn't overwritten then it means
   # we will ignore this observation.
   w = 0.0

   # if there is a "W" entry in dat.info then this means all defaults
   # are over-written and we use this weight directly as is
   if haskey(dat.info, "W")
      # now check that this observation type (E, F, V) exists in dat.D["W"]
      # if yes return the corresponding weight, if not we ignore that
      # observation.
      if haskey(dat.info["W"], obskey)
         w = dat.info["W"][obskey]
      end
   else
      # if no "W" dict exists extract the correct weights from
      # the `weights` dictionary
      # first check what the current configtype has a separate entry
      # and otherwise use the default
      ct = configtype(dat)
      if haskey(weights, ct)
         cfgkey = ct
      else
         cfgkey = "default"
      end
      # now check whether there is an entry for the current observation  type
      if haskey(weights[cfgkey], obskey)
         w = weights[cfgkey][obskey] * wh
      end
   end

   # transform the weights into a vector (if necessary) and return
   if length(w) == 1
      return w * ones(length(o))
   elseif length(w) == length(o)
      return w
   end
   error("_get_weights: length(w) is neither 1 nor length(o)?!?!?")
end




function _regularise!(Ψ::Matrix{T}, Y::Vector{T}, basis, regularisers;
                      verbose=false, Ibasis = :) where {T}
   # TODO: incorporate Ibasis again !
   # assemble the regularisers
   Ψreg = Matrix{Float64}(undef, 0, size(Ψ, 2))
   Yreg = Float64[]
   for reg in regularisers
      if reg isa Matrix
         P = reg
         q = zeros(size(P, 1))
      else
         P, q = Matrix(reg, basis; verbose=verbose)
      end
      Ψreg = vcat(Ψreg, P)
      Yreg = vcat(Yreg, q)
   end
   # check they all have the correct size
   return vcat(Ψ, Ψreg), vcat(Y, Yreg)
end



"""
`get_lsq_system(db::LsqDB; kwargs...) -> Ψ, Y`

Assemble the least squares system + rhs (observations). The `kwargs` can be used
to select a subset of the available data or basis set, and to adjust the
weights by config_type and observation type. See `lsqfit` for a list
of keyword arguments. Allowed/required kwargs are `verbose, weights,
E0, Vref, Ibasis, Itrain, regularisers`.
"""
@noinline function get_lsq_system(db::LsqDB;
                                  verbose=true,
                                  Ibasis = :,
                                  Itrain = :,
                                  E0 = nothing,
                                  Vref = OneBody(E0),
                                  weights = nothing,
                                  regularisers = [])
   weights = _fix_weights!(weights)

   # get the observations vector and the weights vector
   # the Vref potential is subtracted from the observations
   Y, W, Icfg = collect_observations(db, weights, Vref)

   # check for NaNs
   any(isnan, Y) && error("NaN detected in observations vector")
   any(isnan, W) && error("NaN detected in weights vector")


   # convert : into a vector or list
   if Itrain == : ; Itrain = 1:length(Y); end

   # get the right slice of the big fat huge enourmous LSQ matrix
   # the columns are just the basis functions
   Icols = Ibasis
   # the rows are those which have (a) a non-zero weight
   # and (b) the configuration index is part of the training set
   Irows = intersect( findall(W .!= 0) |> sort,
                      findall(in.(Icfg, Ref(Itrain))) )

   # we make this a view but make sure to copy it before applying the weights
   Ψ = @view db.Ψ[Irows, Icols]
   Y = Y[Irows]
   W = W[Irows]

   if any(isnan, Ψ)
      @error("discovered NaNs in LSQ system matrix")
      naninds = findall(isnan.(Ψ))
      n = min(length(naninds), 10)
      @show length(naninds)
      @show naninds[1:n]
      error("discovered NaNs in LSQ system matrix")
   end

   # regularise
   if !isempty(regularisers)
      Ψ, Y = _regularise!(Ψ, Y, db.basis, regularisers; verbose=verbose, Ibasis=Ibasis)
      # adjust the weight-vector to the new system size
      append!(W, ones(length(Y)-length(W)))
   else
      Ψ = collect(Ψ)
   end

   # now rescale Y and Ψ according to W => Y_W, Ψ_W; then the two systems
   #   \| Y_W - Ψ_W c \| -> min  and \| W (Y - Ψ*c)^T \| -> min
   # are equivalent
   Y = Diagonal(W) * Y
   lmul!(Diagonal(W), Ψ)

   # this should be it ...
   return Ψ, Y
end

_show_free_mem() =
      @info("Free Memory: ≈ $(round(Sys.free_memory()*1e-9, digits=2)) GB")


"""
`lsqfit(db; kwargs...) -> IP, errs`

Given the pre-computed least-squares system `db` (cf `LsqDB`) setup a least
squares system, compute the solution and return an interatomic potential
`IP` as well as the associated error estimates `errs`.

## Keyword Arguments:

* `E0` (required unless `Vref` is provided) : energy of the `OneBody` term
* `Vref` : a reference potential from which to start the fit.
* `weights` (required) : A dictionary specifying the weights for different
types of configurations, see `?collect_observations` for how this is specified.
* `Ibasis` : indices of basis functions to be used in the fit, default is `:`
* `verbose` : true or false
* `solver` : -experimental, still need to  write the docs for this-
* `regularisers` : a list of regularisers to be added to the lsq functional.
Each regulariser `R` can be of an arbtirary type but this type must implement the
conversion to matrix
```
Matrix(R, basis; verbose={true,false})
```
If `Areg = Matrix(R, basis)` then this corresponds to adding
`|| Areg * x ||²` to the least squares functional.
* `deldb = false` : experimental - if `true` then the lsq matrix in the `db` is
deleted after assembling the weighted lsq system
* `asmerrs = false` : experimental - if `true` then the lsq matrix is used to
assemble the fitting errors.
* `saveqr = nothing` : experimental - if `saveqr` is a Dict then the factors
of the QR factorisation are stored in it for future use outside of this
function

## Return types

* `IP`: whatever `JuLIP.MLIPs.combine` returns
* `fitinfo::Dict`: adds all sort of information about the fit
"""
@noinline function lsqfit(db::LsqDB;
                solver=(:qr, ), verbose=true,
                Ibasis = :,
                Itrain = :,
                Itest = nothing,
                E0 = nothing,
                Vref = OneBody(E0),
                weights = nothing,
                regularisers = [],
                deldb = false,
                asmerrs = false,
                saveqr = nothing,
                kwargs...)
   weights = _fix_weights!(weights)
   Jbasis = ((Ibasis == Colon()) ? (1:length(db.basis)) : Ibasis)

   verbose && @info("assemble lsq system")
   verbose && _show_free_mem()
   Ψ, Y = get_lsq_system(db; verbose=verbose, Vref=Vref,
                             Ibasis=Ibasis, Itrain = Itrain,
                             weights = weights,
                             regularisers = regularisers)
   verbose && _show_free_mem()
   if deldb
      @info("Deleting database - `db` can no longer be saved to disk")
      db.Ψ = Matrix{Float64}(undef, 0,0)
      db.dbpath = ""
   end
   GC.gc()
   verbose && _show_free_mem()

   κ = 0.0
   if (solver[1] == :qr) || (solver == :qr)
      verbose && @info("solve $(size(Ψ)) LSQ system using QR factorisation")
      qrΨ = qr!(Ψ)
      κ = cond(qrΨ.R)
      GC.gc();
      verbose && @info("cond(R) = $(cond(qrΨ.R))")
      c = qrΨ \ Y
      rel_rms = norm(qrΨ.Q * (qrΨ.R * c) - Y) / norm(Y)

      if saveqr isa Dict
         saveqr["Q"] = qrΨ.Q
         saveqr["R"] = qrΨ.R
         saveqr["Y"] = Y
      end

      qrΨ = nothing

   elseif solver[1] == :svd
      verbose && @info("solve $(size(Ψ)) LSQ system using SVD factorisation")
      ndiscard = solver[2]
      F = svd(Ψ)
      c = F.V[:,1:(end-ndiscard)] * (Diagonal(F.S[1:(end-ndiscard)]) \ (F.U' * Y)[1:(end-ndiscard)])
      rel_rms = norm(Ψ * c - Y) / norm(Y)

   elseif solver[1] == :rrqr
      verbose && @info("solve $(size(Ψ)) LSQ system using Rank-Revealing QR factorisation")
      qrΨ = pqrfact(Ψ, rtol=solver[2])
      verbose && @info("cond(R) = $(cond(qrΨ.R))")
      c = qrΨ \ Y
      rel_rms = norm(Ψ * c - Y) / norm(Y)

   else
      error("unknown `solver` in `lsqfit`")
   end

   # delete the lsq system and gc again
   Ψ = nothing
   GC.gc()
   verbose && _show_free_mem()

   if verbose
      @info("Relative RMSE on training set: $rel_rms")
   end

   IP = JuLIP.MLIPs.combine(db.basis, c)
   if (Vref != nothing) && (Vref != OneBody(0.0))
      IP = SumIP(Vref, IP)
   end

   infodict = asm_fitinfo(db, IP, c, Ibasis, weights,
                          Vref, solver, E0, regularisers, verbose,
                          Itrain, Itest, asmerrs)
   infodict["kappa"] = κ
   GC.gc()
   return IP, infodict
end


function asm_fitinfo(db, IP, c, Ibasis, weights,
                     Vref, solver, E0, regularisers, verbose,
                     Itrain = :, Itest = nothing, asmerrs=true)
   if Ibasis isa Colon
      Jbasis = collect(1:length(db.basis))
   else
      Jbasis = Ibasis
   end

   cfgtypes = setdiff(unique(configtype.(db.configs)), weights["ignore"])

   # compute errors TODO: still need to fix this!
   if asmerrs
      verbose && @info("Assemble errors table")
      @warn("new error implementation... redo this part please ")
      errs = Err.lsqerrors(db, c, Jbasis;
               cfgtypes=cfgtypes, Vref=OneBody(E0), Icfg=Itrain)
      if Itest != nothing
         errtest = Err.lsqerrors(db, c, Jbasis;
                  cfgtypes=cfgtypes, Vref=OneBody(E0), Icfg=Itest)
      else
         errtest = Dict()
      end
   end
   # --------------------------------------------------------------------
   # ASSEMBLE INFO DICT
   # --------------------------------------------------------------------
   verbose && @info("Assemble Information about the fit")

   # Julia Version Info
   iob = IOBuffer()
   versioninfo(iob)
   juliainfo = String(take!(iob))

   infodict = Dict("solver" => String(solver[1]),
                   "E0"     => E0,
                   "Ibasis" => Vector{Int}(Jbasis),
                   "c"      => c,
                   "dbpath" => dbpath(db),
                   "weights" => weights,
                   "regularisers"  => Dict.(regularisers),
                   "juliaversion"  => juliainfo,
                   "IPFitting_version" => get_pkg_info("IPFitting"),
                  )

   if asmerrs
      infodict["errors"] = errs
      infodict["errtest"] = errtest
   end
   # --------------------------------------------------------------------
   return infodict
end

import Pkg

function get_pkg_info(pkg::AbstractString)
    pkgs = [Pkg.API.check_package_name(pkg)]
    ctx = Pkg.Types.Context()
    Pkg.API.project_resolve!(ctx.env, pkgs)
    Pkg.API.project_deps_resolve!(ctx.env, pkgs)
    Pkg.API.manifest_resolve!(ctx.env, pkgs)
    Pkg.API.ensure_resolved(ctx.env, pkgs)
    i = Pkg.Display.status(ctx, pkgs)[1]
    return Dict("name" => i.name,
                "uuid" => string(i.uuid),
                "ver" => string(i.new.ver))
end


_fix_weights!(::Nothing) = _fix_weights!(Dict{String, Any}())

function _fix_weights!(weights::Dict)
   if !haskey(weights, "ignore")
      weights = Dict{Any,Any}(weights)
      weights["ignore"] = String[]
   end
   if !haskey(weights, "default")
      weights["default"] = Dict("E" => 1.0, "F" => 1.0, "V" => 1.0)
   end
   return weights
end


end





# @noinline function onb(db::LsqDB;
#                          solver=(:qr, ), verbose=true,
#                          Ibasis = :,
#                          Itrain = :,
#                          E0 = nothing,
#                          Vref = OneBody(E0),
#                          weights = nothing,
#                          regularisers = [],
#                          combineIP = nothing,
#                          kwargs...)
#
#    verbose && @info("assemble lsq system")
#    Ψ, _ = get_lsq_system(db; verbose=verbose, Vref=Vref,
#                              Ibasis=Ibasis, Itrain=Itrain,
#                              weights = weights,
#                              regularisers = regularisers,
#                              kwargs...)
#    @assert solver[1] == :qr
#    verbose && @info("QR-factorize Ψ, size=$(size(Ψ))")
#    qrΨ = qr(Ψ)
#    verbose && @info("cond(R) = $(cond(qrΨ.R))")
#    Rinv = pinv(qrΨ.R)
#    basis = db.basis[Ibasis]
#    onb = []
#    for n = 1:size(Rinv, 2)
#       push!(onb, combineIP(basis, Rinv[:,n]))
#    end
#    return [b for b in onb]
# end
