
"""
`module Lsq`

This sub-modules contains code to convert a database `LsqDB` into
linear LSQ systems, apply suitable weights, solve the resulting systems
and return an NBodyIP.

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

import JuLIP, NBodyIPFitting, StatsBase

using StaticArrays
using JuLIP: AbstractCalculator, Atoms
using NBodyIPs: OneBody, NBodyIP
using NBodyIPFitting: Dat, LsqDB, weighthook, observations,
                      observation, hasobservation
using NBodyIPFitting.Data: configtype
using NBodyIPFitting.DataTypes: ENERGY
using NBodyIPFitting.DB: dbpath, _nconfigs, matrows

using LinearAlgebra: lmul!, Diagonal, qr, cond, norm
using InteractiveUtils: versioninfo

const Err = NBodyIPFitting.Errors

export lsqfit, onb



_keys(configweights, obsweights) = collect(keys(configweights)),
                                    collect(keys(obsweights))

function _random_subset(db::LsqDB, configtype::AbstractString, p::Real)
   @assert 0 <= p <= 1
   # number of configurations in db for this configtype
   nconfigs = _nconfigs(db, configtype)
   # number of samples we want
   nsamples = ceil(Int, p * nconfigs)
   # draw a random subset
   return StatsBase.sample(1:nconfigs, nsamples, replace = false, ordered = true)
end



"""
`observations(db::LsqDB, configweights::Dict, dataweights::Dict) -> Vector{Float64}`

construct a vector of observations from the database `db`, specifically from
`db.data`. Only those configurations are considered whose config_type is
in `keys(configweights)` and only those datatypes are loaded whose
ids are in `keys(dataweights)`.
"""
function collect_observations(db::LsqDB,
                              configtypes::Vector{String},
                              obstypes::Vector{String},
                              configweights::Dict,
                              obsweights::Dict, E0 )
                              # , Iconfigs::Dict )  # TODO
   nrows = size(db.Ψ, 1)  # total number of observations if we collect
                          # everything in the database
   Y = zeros(Float64, nrows)
   W = zeros(Float64, nrows)
   for (obskey, dat, _) in observations(db)  # obskey ∈ {"E","F",...}; d::Dat
      if !(obskey in obstypes)
         continue
      end
      irows = matrows(dat, obskey)
      ct = configtype(dat)
      # ----- Observation ------
      obs = observation(obskey, dat)
      # TODO: fix this hack!!!! (put in basis function constraints)
      if obskey == "E"
         obs = copy(obs) .- length(dat) * E0
      end
      Y[irows] .= obs
      # ------ Weights --------
      # modify the weights from extra information in the dat structure
      W[irows] .= _get_weights(configweights[ct], # .weight, TODO
                               obsweights[obskey],
                               weighthook(obskey, dat),
                               dat, obskey, obs)
   end
   return Y, W
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
                      verbose=false) where {T}
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
`struct ConfigFitInfo` : stores `configtype`-dependent
weights and also proportion of configurations to fit to
"""
struct ConfigFitInfo
   weight::Float64
   proportion::Float64
   order::Symbol
end

ConfigFitInfo(cfi::ConfigFitInfo) = deepcopy(cfi)
ConfigFitInfo(cfi::Tuple) = ConfigFitInfo(cfi...)
ConfigFitInfo(w::Real) = ConfigFitInfo(w, 1.0, :ignore)
ConfigFitInfo(w::Real, p::Real) = ConfigFitInfo(w, p, :ignore)
ConfigFitInfo(w::Real, p::Real, o::Any) = ConfigFitInfo(w, p, Symbol(o))



"""
`get_lsq_system(db::LsqDB; kwargs...) -> Ψ, Y`

Assemble the least squares system + rhs (observations). The `kwargs` can be used
to select a subset of the available data or basis set, and to adjust the
weights by config_type and observation type. See `lsqfit` for a list
of keyword arguments. Allowed kwargs are `verbose, configweights, obsweights,
E0, Ibasis`.
"""
@noinline function get_lsq_system(db::LsqDB; verbose = true,
                        configweights = nothing,
                        obsweights = nothing,
                        E0 = nothing,
                        Ibasis = :,
                        regularisers = nothing )
   # we need to be able to call `length` on `Ibasis`
   Jbasis = ((Ibasis == Colon()) ? (1:length(db.basis)) : Ibasis)

   # TODO: put this back in!!!
   # # restrict the configurations that we want
   # Iconfigs = Dict( key => _random_subset(db, key, val.proportion)
   #                  for (key, val) in configweights )
   Iconfigs = nothing

   # # TODO => create some suitable "hooks"
   # E0 = lsq.basis[1]()
   # # and while we're at it, subtract E0 from Y
   # Y[idx] -= E0 * len

   # fix some ordering of the configtypes and obstypes
   # even though this can be inferred from configweights and obsweights we
   # need to be paranoid that the ordering does not change!
   # TODO: this is no longer needed?!?!?
   configtypes, obstypes = _keys(configweights, obsweights)
   # get the observations vector and the weights vector
   Y, W = collect_observations(db, configtypes, obstypes, configweights, obsweights,
                       E0) # , Iconfigs)
   # check for NaNs
   any(isnan, Y) && @error("NaN detected in observations vector")
   any(isnan, W) && @error("NaN detected in weights vector")

   # get the slice of the big fat huge enourmous LSQ matrix
   Icols = Ibasis
   Irows = findall(W .!= 0) |> sort
   # we can't keep this a view since we need to multiply by weights
   Ψ = db.Ψ[Irows, Icols]
   Y = Y[Irows]
   W = W[Irows]
   any(isnan, Ψ) && @error("discovered NaNs in LSQ system matrix")

   # now rescale Y and Ψ according to W => Y_W, Ψ_W; then the two systems
   #   \| Y_W - Ψ_W c \| -> min  and (Y - Ψ*c)^T W (Y - Ψ*x) => MIN
   # are equivalent
   # >>>>>> W .= sqrt.(W) <<<<<<< TODO: which one is it?
   @. Y = Y * W
   lmul!(Diagonal(W), Ψ)

   # regularise
   if regularisers != nothing
      Ψ, Y = _regularise!(Ψ, Y, db.basis[Jbasis], regularisers; verbose=verbose)
   end

   # this should be it ...
   return Ψ, Y
end


@noinline function onb(db::LsqDB;
             solver=(:qr, ), verbose=true, E0 = nothing,
             Ibasis = :, configweights=nothing, obsweights = nothing,
             regularisers = [],
             kwargs...)
   @assert E0 != nothing
   Jbasis = ((Ibasis == Colon()) ? (1:length(db.basis)) : Ibasis)

   verbose && @info("assemble lsq system")
   Ψ, _ = get_lsq_system(db; verbose=verbose, E0=E0, Ibasis=Ibasis,
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
      push!(onb, NBodyIP(basis, Rinv[:,n]))
   end
   return [b for b in onb]
end


"""
`lsqfit(db; kwargs...) -> IP, errs`

Given the pre-computed least-squares system `db` (cf `LsqDB`) setup a least
squares system, compute the solution and return an interatomic potential
`IP` as well as the associated error estiamates `errs`.

## Keyword Arguments:

* `E0` (required) : energy of the `OneBody` term
* `configweights` (required) : A dictionary specifying the weights for different
types of configurations, e.g.,
```
configweights = Dict("solid" => 10.0, "liquid" => 0.1)
```
The keys, `["solid", "liquid"]` in the above example, specify which configtypes
are to be fitted - all other configtypes are ignored.
* `dataweights` (required) : a `Dict` specifying how different kinds of
observations (data) should be weighted, e.g.,
```
dataweights = Dict("E" => 100.0, "F" => 1.0, "V" => 0.1)
```
* `Ibasis` : indices of basis functions to be used in the fit, default is `:`
* `verbose` : true or false
* `solver` : at the moment this should be ignored, only admissible choice
is `:qr`. On request we can try others.

## Return types

* `IP::NBodyIP`: see documentation of `NBodyIPS.jl`
* `errs::LsqErrors`: stores individual RMSE and MAE for the different
configtypes and datatypes. Use `table_relative(errs)` and `table_absolute(errs)`
to display these as tables and `rmse, mae` to access individual errors.
"""
@noinline function lsqfit(db::LsqDB;
                solver=(:qr, ), verbose=true, E0 = nothing,
                Ibasis = :, configweights=nothing, obsweights = nothing,
                regularisers = [],
                kwargs...)
   @assert E0 != nothing
   Jbasis = ((Ibasis == Colon()) ? (1:length(db.basis)) : Ibasis)

   verbose && @info("assemble lsq system")
   Ψ, Y = get_lsq_system(db; verbose=verbose, E0=E0, Ibasis=Ibasis,
                             configweights = configweights,
                             obsweights = obsweights,
                             regularisers = regularisers,
                             kwargs...)

   if (solver[1] == :qr) || (solver == :qr)
      verbose && @info("solve $(size(Ψ)) LSQ system using QR factorisation")
      qrΨ = qr(Ψ)
      verbose && @info("cond(R) = $(cond(qrΨ.R))")
      c = qrΨ \ Y

   elseif solver[1] == :svd
      verbose && @info("solve $(size(Ψ)) LSQ system using SVD factorisation")
      ndiscard = solver[2]
      F = svdfact(Ψ)
      c = F[:V][:,1:(end-ndiscard)] * (Diagonal(F[:S][1:(end-ndiscard)]) \ (F[:U]' * Y)[1:(end-ndiscard)])

   else
      error("unknown `solver` in `lsqfit`")
   end

   if verbose
      rel_rms = norm(Ψ * c - Y) / norm(Y)
      @info("Relative RMSE on training set: $rel_rms")
   end

   # compute errors
   verbose && @info("Assemble errors table")
   @warn("new error implementation... redo this part please ")
   errs = Err.lsqerrors(db, c, Jbasis; cfgtypes=keys(configweights), E0=E0)

   if E0 != 0
      basis = [ OneBody(E0); db.basis[Ibasis] ]
      c = [1.0; c]
   else
      basis = db.basis[Ibasis]
   end

   # --------------------------------------------------------------------
   # ASSEMBLE INFO DICT
   # --------------------------------------------------------------------
   verbose && @info("Assemble Information about the fit")

   # Julia Version Info
   iob = IOBuffer()
   versioninfo(iob)
   juliainfo = String(take!(iob))

   # NBodyIPs and NBodyIPFitting Version Info
   # TODO: put back in
   # nbipinfo, nbipfitinfo = get_git_info()

   # number of configurations for each configtype
   # TODO: put back in
   # numconfigs = Dict{String, Int}()
   # for ct in keys(db.data_groups)
   #    cn = configname(ct)
   #    if !haskey(numconfigs, cn)
   #       numconfigs[cn] = 0
   #    end
   #    numconfigs[cn] += length(db.data_groups[ct])
   # end

   infodict = Dict("errors" => errs,
                   "solver" => String(solver[1]),
                   "E0" => E0,
                   "Ibasis" => Vector{Int}(Jbasis),
                   "dbpath" => dbpath(db),
                   "configweights" => configweights,
                   "confignames" => keys(configweights),
                   "obsweights" => obsweights,
                   "regularisers" => Dict.(regularisers),
                   "juliaversion" => juliainfo,
                   # "NBodyIPs_version" => nbipinfo,
                   # "NBodyIPFitting_version" => nbipfitinfo,
                   # "numconfigs" => numconfigs
                  )
   # --------------------------------------------------------------------


   return NBodyIP(basis, c), infodict
end


# ## TODO: THIS NEEDS A REWRITE !!!!
# function get_git_info()
#    @warn("Storing Package git version info is no longer supported; this needs a rewrite.")
#    # nbipinfo = read(`cat $(Pkg.dir("NBodyIPs")*"/.git/refs/heads/master")`, String)[1:end-1]
#    # nbipfitinfo = read(`cat $(Pkg.dir("NBodyIPFitting")*"/.git/refs/heads/master")`, String)[1:end-1]
#    # nbipinfo = readstring(`git -C $(Pkg.dir("NBodyIPs")) rev-parse HEAD`)
#    # nbipfitinfo = readstring(`git -C $(Pkg.dir("NBodyIPFitting")) rev-parse HEAD`)
#    return "", "" # nbipinfo, nbipfitinfo
# end

end
