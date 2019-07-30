
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

using LinearAlgebra: lmul!, Diagonal, qr, cond, norm, svd
using InteractiveUtils: versioninfo
using LowRankApprox

const Err = IPFitting.Errors

export lsqfit, onb



"""
`collect_observations(db::LsqDB, configweights::Dict, obsweights::Dict, Vref) -> Vector{Float64}`

construct a vector of observations from the database `db`, specifically from
`db.data`. Only those configurations are considered whose config_type is
in `keys(configweights)` and only those datatypes are loaded whose
ids are in `keys(dataweights)`.
"""
function collect_observations(db::LsqDB,
                              configweights::Dict,
                              obsweights::Dict, Vref )

   nrows = size(db.Ψ, 1)  # total number of observations if we collect
                          # everything in the database
   Y = zeros(Float64, nrows)
   W = zeros(Float64, nrows)
   Icfg = zeros(Int, nrows)
   for (obskey, dat, icfg) in observations(db)  # obskey ∈ {"E","F",...}; d::Dat
      # check that we want to fit this observation
      if !haskey(obsweights, obskey); continue; end
      if !haskey(configweights, configtype(dat)); continue; end

      irows = matrows(dat, obskey)
      ct = configtype(dat)

      # ----- Observation ------
      obs = observation(obskey, dat)
      # If Vref != nothing the it is a calculator and we can subtract
      # its value => this allows us to fit from a reference potential
      if haskey(dat.info, "Vref")
         obs = obs - dat.info["Vref"][obskey]
      elseif Vref != nothing
         obs = obs - vec_obs(obskey, eval_obs(obskey, Vref, Atoms(dat)))
      end
      Y[irows] .= obs
      Icfg[irows] .= icfg
      # ------ Weights --------
      # modify the weights from extra information in the dat structure
      W[irows] .= _get_weights(configweights[ct],
                               obsweights[obskey],
                               weighthook(obskey, dat),
                               dat, obskey, obs)
   end
   return Y, W, Icfg
end


function _get_weights(ctweight, dtweights_dt, wh, dat, dt, o)
   # if there is a "W" entry in dat.D then this means all defaults
   # are over-written
   if haskey(dat.info, "W")
      # now check that this observation type (E, F, V) exists in dat.D["W"]
      # if yes return the corresponding weight, if not make it zero
      # which means that we will simply ignore this observation
      if haskey(dat.info["W"], dt)
         w = dat.info["W"][dt]
      else
         w = 0.0
      end
   else
      # if no "W" dict exists use the default weights
      w = ctweight * dtweights_dt * wh
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
of keyword arguments. Allowed kwargs are `verbose, configweights, obsweights,
E0, Ibasis`.
"""
@noinline function get_lsq_system(db::LsqDB;
                         solver=(:qr, ), verbose=true,
                         Ibasis = :, Itrain = :,
                         E0 = nothing,
                         Vref = OneBody(E0),
                         configweights = nothing,
                         obsweights = nothing,
                         regularisers = [],
                         kwargs...)

   # get the observations vector and the weights vector
   # the Vref potential is subtracted from the observations
   Y, W, Icfg = collect_observations(db, configweights, obsweights, Vref)

   # check for NaNs
   any(isnan, Y) && @error("NaN detected in observations vector")
   any(isnan, W) && @error("NaN detected in weights vector")


   # convert : into a vector or list
   if Itrain == : ; Itrain = 1:length(Y); end

   # get the right slice of the big fat huge enourmous LSQ matrix
   # the columns are just the basis functions
   Icols = Ibasis
   # the rows are those which have (a) a non-zero weight
   # and (b) the configuration index is part of the training set
   Irows = intersect(findall(W .!= 0) |> sort, findall(in.(Icfg, Ref(Itrain))))

   # we can't keep this a view since we need to multiply by weights
   Ψ = db.Ψ[Irows, Icols]
   Y = Y[Irows]
   W = W[Irows]

   if any(isnan, Ψ)
      @info("discovered NaNs in LSQ system matrix")

      naninds = findall(isnan.(Ψ))
      n = min(length(naninds), 10)
      @show length(naninds)
      @show naninds[1:n]

      @error("discovered NaNs in LSQ system matrix")
   end

   # now rescale Y and Ψ according to W => Y_W, Ψ_W; then the two systems
   #   \| Y_W - Ψ_W c \| -> min  and \| W (Y - Ψ*c)^T \| -> min
   # are equivalent
   @. Y = Y * W
   lmul!(Diagonal(W), Ψ)

   # regularise
   if !isempty(regularisers)
      Ψ, Y = _regularise!(Ψ, Y, db.basis, regularisers; verbose=verbose, Ibasis=Ibasis)
   end

   # this should be it ...
   return Ψ, Y
end


@noinline function onb(db::LsqDB;
                         solver=(:qr, ), verbose=true,
                         Ibasis = :,
                         Itrain = :,
                         E0 = nothing,
                         Vref = OneBody(E0),
                         configweights = nothing,
                         obsweights = nothing,
                         regularisers = [],
                         combineIP = nothing,
                         kwargs...)

   verbose && @info("assemble lsq system")
   Ψ, _ = get_lsq_system(db; verbose=verbose, Vref=Vref,
                             Ibasis=Ibasis, Itrain=Itrain,
                             configweights = configweights,
                             obsweights = obsweights,
                             regularisers = regularisers,
                             kwargs...)
   @assert solver[1] == :qr
   verbose && @info("QR-factorize Ψ, size=$(size(Ψ))")
   qrΨ = qr(Ψ)
   verbose && @info("cond(R) = $(cond(qrΨ.R))")
   Rinv = pinv(qrΨ.R)
   basis = db.basis[Ibasis]
   onb = []
   for n = 1:size(Rinv, 2)
      push!(onb, combineIP(basis, Rinv[:,n]))
   end
   return [b for b in onb]
end


"""
`lsqfit(db; kwargs...) -> IP, errs`

Given the pre-computed least-squares system `db` (cf `LsqDB`) setup a least
squares system, compute the solution and return an interatomic potential
`IP` as well as the associated error estiamates `errs`.

## Keyword Arguments:

* `E0` (required unless `Vref` is provided) : energy of the `OneBody` term
* `Vref` : a reference potential from which to start the fit.
* `configweights` (required) : A dictionary specifying the weights for different
types of configurations, e.g.,
```
configweights = Dict("solid" => 10.0, "liquid" => 0.1)
```
The keys, `["solid", "liquid"]` in the above example, specify which configtypes
are to be fitted - all other configtypes are ignored.
* `obsweights` (required) : a `Dict` specifying how different kinds of
observations should be weighted, e.g.,
```
obsweights = Dict("E" => 100.0, "F" => 1.0, "V" => 0.1)
```
* `Ibasis` : indices of basis functions to be used in the fit, default is `:`
* `verbose` : true or false
* `solver` : -experimental, still need to  write the docs for this-
* `regularisers` : a list of regularisers to be added to the lsq functional.
Each regulariser `R` can be of an arbtirary type but this type must implement the
conversion to matrix
* `combineIP` : this is a required kwarg, it tells `lsqfit` how to
combine `basis` and `coeffs` into an IP by calling `combineIP(basis, coeffs)`;
use `(b, c) -> c` to just return the coefficients.
```
Matrix(R, basis; verbose={true,false})
```
If `Areg = Matrix(R, basis)` then this corresponds to adding
`|| Areg * x ||²` to the least squares functional.

## More on Weights

The `configweights` and `obsweights` dictionaries specify weights as follows:
for an observation `o` from a config `cfg` where the `obsweights is `wo` and
the configweight is `wc` the weight on this observation `o` will be `w = w0*wc`.
This given a diagonal weight matrix `W`. The least squares functional then becomes
|| W (A * x - Y) ||²` where `x` are the unknown coefficients.

## Return types

* `IP`: whatever `combineIP` returns
* `errs::LsqErrors`: stores individual RMSE and MAE for the different
configtypes and datatypes. Use `table_relative(errs)` and `table_absolute(errs)`
to display these as tables and `rmse, mae` to access individual errors.
"""
@noinline function lsqfit(db::LsqDB;
                solver=(:qr, ), verbose=true,
                Ibasis = :,
                Itrain = :,
                Itest = nothing,
                E0 = nothing,
                Vref = OneBody(E0),
                configweights = nothing,
                obsweights = nothing,
                regularisers = [],
                combineIP = nothing,
                kwargs...)

   Jbasis = ((Ibasis == Colon()) ? (1:length(db.basis)) : Ibasis)

   verbose && @info("assemble lsq system")
   Ψ, Y = get_lsq_system(db; verbose=verbose, Vref=Vref, Ibasis=Ibasis, Itrain = Itrain,
                             configweights = configweights,
                             obsweights = obsweights,
                             regularisers = regularisers,
                             kwargs...)

   if (solver[1] == :qr) || (solver == :qr)
      verbose && @info("solve $(size(Ψ)) LSQ system using QR factorisation")
      qrΨ = qr(Ψ)
      verbose && @info("cond(R) = $(cond(qrΨ.R))")
      c = qrΨ \ Y
      qrΨ = nothing; GC.gc();

   elseif solver[1] == :svd
      verbose && @info("solve $(size(Ψ)) LSQ system using SVD factorisation")
      ndiscard = solver[2]
      F = svd(Ψ)
      c = F.V[:,1:(end-ndiscard)] * (Diagonal(F.S[1:(end-ndiscard)]) \ (F.U' * Y)[1:(end-ndiscard)])
   elseif solver[1] == :rrqr
       verbose && @info("solve $(size(Ψ)) LSQ system using Rank-Revealing QR factorisation")
       qrΨ = pqrfact(Ψ, rtol=solver[2])
       verbose && @info("cond(R) = $(cond(qrΨ.R))")
       c = qrΨ \ Y
   else
      error("unknown `solver` in `lsqfit`")
   end

   if verbose
      rel_rms = norm(Ψ * c - Y) / norm(Y)
      @info("Relative RMSE on training set: $rel_rms")
   end

   Ψ = nothing; GC.gc()

   IP = JuLIP.MLIPs.combine(db.basis, c)
   if (Vref != nothing) && (Vref != OneBody(0.0))
      IP = SumIP(Vref, IP)
   end

   infodict = asm_fitinfo(db, IP, c, Ibasis, configweights, obsweights,
                          Vref, solver, E0, regularisers, verbose,
                          Itrain, Itest)

   return IP, infodict
end


function asm_fitinfo(db, IP, c, Ibasis, configweights, obsweights,
                     Vref, solver, E0, regularisers, verbose,
                     Itrain = :, Itest = nothing)
   if Ibasis isa Colon
      Jbasis = collect(1:length(db.basis))
   else
      Jbasis = Ibasis
   end
   # compute errors TODO: still need to fix this!
   verbose && @info("Assemble errors table")
   @warn("new error implementation... redo this part please ")
   errs = Err.lsqerrors(db, c, Jbasis;
            cfgtypes=keys(configweights), Vref=OneBody(E0), Icfg=Itrain)
   if Itest != nothing
      errtest = Err.lsqerrors(db, c, Jbasis;
               cfgtypes=keys(configweights), Vref=OneBody(E0), Icfg=Itest)
   else
      errtest = Dict()
   end

   # --------------------------------------------------------------------
   # ASSEMBLE INFO DICT
   # --------------------------------------------------------------------
   verbose && @info("Assemble Information about the fit")

   # Julia Version Info
   iob = IOBuffer()
   versioninfo(iob)
   juliainfo = String(take!(iob))

   infodict = Dict("errors" => errs,
                   "errtest" => errtest,
                   "solver" => String(solver[1]),
                   "E0"     => E0,
                   "Ibasis" => Vector{Int}(Jbasis),
                   "c"      => c,
                   "dbpath" => dbpath(db),
                   "configweights" => configweights,
                   "confignames"   => keys(configweights),
                   "obsweights"    => obsweights,
                   "regularisers"  => Dict.(regularisers),
                   "juliaversion"  => juliainfo,
                   "IPFitting_version" => get_pkg_info("IPFitting"),
                  )
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



end
