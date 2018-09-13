
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

using StaticArrays
using JuLIP: AbstractCalculator, Atoms
using NBodyIPs: OneBody, NBodyIP
using NBodyIPFitting: Dat, LsqDB, data, weighthook
using NBodyIPFitting.Data: observation, hasobservation, configname, configtype
using NBodyIPFitting.DataTypes: ENERGY
import NBodyIPFitting
const Err = NBodyIPFitting.Errors

export lsqfit

_keys(configweights, dataweights) = collect(keys(configweights)),
                                    collect(keys(dataweights))

"""
`observations(db::LsqDB, configweights::Dict, dataweights::Dict) -> Vector{Float64}`

construct a vector of observations from the database `db`, specifically from
`db.data`. Only those configurations are considered whose config_type is
in `keys(configweights)` and only those datatypes are loaded whose
ids are in `keys(dataweights)`.
"""
function observations(db::LsqDB,
                      configtypes::Vector{String},
                      datatypes::Vector{String},
                      ctweights::Dict,
                      dtweights::Dict, E0 )
   Y = Float64[]
   W = Float64[]
   # loop through configuration groups
   for ct in configtypes
      ctweight = ctweights[ct]
      # loop through the datatypes / observation types to be assembled
      for dt in datatypes
         # loop through the observations (x, f(x)) of the current configuration group
         for dat in db.data_groups[ct]
            # TODO: move this `if` after the `for dt` line
            if hasobservation(dat, dt)
               w = weighthook(dt, dat)
               o = observation(dat, dt)
               # TODO: This is a hack => can we replace it with a hook?
               if dt == ENERGY # subtract the 1-body reference energy
                  o = copy(o)
                  o[1] -= length(dat) * E0
               end
               append!(Y, o)
               append!(W, ctweight * dtweights[dt] * (w .* ones(length(o))) )
            end
         end
      end
   end
   return Y, W
end


# TODO: this function here suggests that the ordering we are usig at the moment
#       is perfectly fine, and is in fact the more performant ordering
#       and we should not switch it
function lsq_matrix!(Ψ, db, configtypes, datatypes, Ibasis)
   idx = 0
   for ct in configtypes
      k = db.kron_groups[ct]
      for dt in datatypes
         if haskey(k, dt)
            # block[:, i, j] = <data(i), basis(j)>_{dt}
            block = k[dt]
            nrows = size(block , 1) * size(block , 2)
            Ψ[idx+1:idx+nrows, :] .= reshape(block[:, :, Ibasis], nrows, :)
            idx += nrows
         end
      end
   end
   return nothing
end


function _regularise!(Ψ::Matrix{T}, Y::Vector{T}, basis, regularisers) where {T}
   # assemble the regularisers
   Ψreg = Matrix{Float64}(0, size(Ψ, 2))
   Yreg = Float64[]
   for reg in regularisers
      if reg isa Matrix
         P = reg
      else
         P = Matrix(reg, basis)
      end
      Ψreg = vcat(Ψreg, P)
      Yreg = vcat(Yreg, zeros(size(P, 1)))
   end
   # check they all have the correct size
   return vcat(Ψ, Ψreg), vcat(Y, Yreg)
end

"""
generate a new dictionary of configtype weights where configweights is a
dictionary of configname weights
"""
function _extend_configweights(configweights::Dict, configtypes)
   names = keys(configweights) |> collect
   newweights = Dict{String, Float64}()
   for ct in configtypes
      cn = configname(ct)
      if cn in names
         newweights[ct] = configweights[cn]
      end
   end
   return newweights
end



"""
`get_lsq_system(db::LsqDB; kwargs...) -> Ψ, Y`

Assemble the least squares system + rhs (observations). The `kwargs` can be used
to select a subset of the available data or basis set, and to adjust the
weights by config_type and observation type. See `lsqfit` for a list
of keyword arguments. Allowed kwargs are `verbose, configweights, dataweights,
E0, Ibasis`.
"""
function get_lsq_system(db::LsqDB; verbose = true,
                        configweights = nothing,
                        dataweights = nothing,
                        E0 = nothing,
                        Ibasis = :,
                        regularisers = nothing )
   # we need to be able to call `length` on `Ibasis`
   Jbasis = ((Ibasis == Colon()) ? (1:length(db.basis)) : Ibasis)

   configweights = _extend_configweights(configweights, keys(db.data_groups))

   # # TODO TODO TODO => create some suitable "hooks"
   # E0 = lsq.basis[1]()
   # # and while we're at it, subtract E0 from Y
   # Y[idx] -= E0 * len

   # fix some ordering of the configtypes and datatypes
   # even though this can be inferred from configweights and dataweights we
   # need to by paranoid that the ordering does not change!
   configtypes, datatypes = _keys(configweights, dataweights)
   # get the observations vector and the weights vector
   Y, W = observations(db, configtypes, datatypes, configweights, dataweights, E0)
   # check for NaNs
   any(isnan, Y) && error("NaN detected in observations vector")
   any(isnan, W) && error("NaN detected in weights vector")

   # allocate and assemble the big fat huge enourmous LSQ matrix
   Ψ = zeros(length(Y), length(Jbasis))
   lsq_matrix!(Ψ, db, configtypes, datatypes, Jbasis)
   any(isnan, Ψ) && error("discovered NaNs in LSQ system matrix")

   # remove anything with zero-weight
   Idata = find(W .!== 0.0) |> sort
   Y, W, Ψ = Y[Idata], W[Idata], Ψ[Idata, :]

   # now rescale Y and Ψ according to W => Y_W, Ψ_W; then the two systems
   #   \| Y_W - Ψ_W c \| -> min  and (Y - Ψ*c)^T W (Y - Ψ*x) => MIN
   # are equivalent
   # W .= sqrt.(W)
   Y .*= W
   scale!(W, Ψ)

   # regularise
   if regularisers != nothing
      Ψ, Y = _regularise!(Ψ, Y, db.basis[Jbasis], regularisers)
   end

   # this should be it ...
   return Ψ, Y
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
function lsqfit(db::LsqDB;
                solver=:qr, verbose=true, E0 = nothing,
                Ibasis = :, configweights=nothing, dataweights = nothing,
                regularisers = [],
                kwargs...)
   @assert solver == :qr
   @assert E0 != nothing
   Jbasis = ((Ibasis == Colon()) ? (1:length(db.basis)) : Ibasis)

   verbose && info("assemble lsq system")
   Ψ, Y = get_lsq_system(db; verbose=verbose, E0=E0, Ibasis=Ibasis,
                             configweights = configweights,
                             dataweights = dataweights,
                             regularisers = regularisers
                             kwargs...)

   verbose && info("solve $(size(Ψ)) LSQ system using QR factorisation")
   qrΨ = qrfact(Ψ)
   verbose && @show cond(qrΨ[:R])
   c = qrΨ \ Y

   if verbose
      rel_rms = norm(Ψ * c - Y) / norm(Y)
      verbose && println("Relative RMSE on training set: ", rel_rms)
   end

   # compute errors
   errs = Err.lsqerrors(db, c, Jbasis; confignames=keys(configweights), E0=E0)

   if E0 != 0
      basis = [ OneBody(E0); db.basis[Ibasis] ]
      c = [1.0; c]
   else
      basis = db.basis[Ibasis]
   end

   iob = IOBuffer()
   versioninfo(iob)
   juliainfo = String(iob)

   nbipinfo = readstring(`git -C $(Pkg.dir("NBodyIPs")) rev-parse HEAD`)
   nbipfitinfo = readstring(`git -C $(Pkg.dir("NBodyIPFitting")) rev-parse HEAD`)


   info = Dict("errors" => Dict(errs),
               "solver" => String(solver),
               "E0" => E0,
               "Ibasis" => Vector{Int}(Ibasis),
               "dbpath" => dbpath(db),
               "configweights" => configweights,
               "confignames" => confignames,
               "dataweights" => dataweights,
               "regularisers" => string.(typeof.(regularisers)),
               "juliaversion" => juliainfo,
               "NBodyIPs_version" => nbipinfo,
               "NBodyIPFitting_version" => nbipfitinfo
         )

   return NBodyIP(basis, c), info
end




end
