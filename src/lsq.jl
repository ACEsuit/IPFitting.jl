
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
using NBodyIPFitting: Dat, LsqDB, data
using NBodyIPFitting.Data: observation, hasobservation


export fit_nbodyip

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
                      configtypes::Vector{String}, datatypes::Vector{String},
                      ctweights::Dict, dtweights::Dict)
   Y = Float64[]
   W = Float64[]
   for ct in configtypes
      ctweight = ctweights[ct]
      for dt in datatypes
         for dat in db.data_groups[ct]
            if hasobservation(dat, dt)
               o = observation(dat, dt)
               append!(Y, o)
               append!(W, ctweight * dtweights[dt] * ones(length(o)) )
            end
         end
      end
   end
   return Y, W
end


# TODO: this function here suggests that the ordering we are usig at the moment
#       is perfectly fine, and is in fact the more performance ordering
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



function get_lsq_system(db::LsqDB; verbose = true,
                        configweights = nothing,
                        dataweights = nothing,
                        E0 = nothing,
                        Ibasis = : )

   if Ibasis == Colon()
      Jbasis = 1:length(db.basis)
   else
      Jbasis = Ibasis
   end

   # # reference energy => we assume the first basis function is 1B
   # # TODO TODO TODO => create some suitable "hooks"
   # E0 = lsq.basis[1]()
   # # and while we're at it, subtract E0 from Y
   # Y[idx] -= E0 * len

   # fix some ordering of the configtypes and datatypes
   # even though this can be inferred from configweights and dataweights we
   # need to by paranoid that the ordering does not change!
   configtypes, datatypes = _keys(configweights, dataweights)
   # get the observations vector and the weights vector
   Y, W = observations(db, configtypes, datatypes, configweights, dataweights)
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
   W .= sqrt.(W)
   Y .*= W
   scale!(W, Ψ)

   # TODO
   # # regularise
   # if regulariser != nothing
   #    P = regulariser(lsq.basis[Jbasis])
   #    Ψ, Y = regularise(Ψ, Y, P)
   # end

   # this should be it ...
   return Ψ, Y
end


function fit(db::LsqDB; solver=:qr, verbose=true, E0 = nothing, Ibasis = :, kwargs...)
   @assert solver == :qr
   @assert E0 != nothing

   verbose && info("assemble lsq system")
   Ψ, Y = get_lsq_system(db; verbose=verbose, E0=E0, Ibasis=Ibasis, kwargs...)

   verbose && info("solve $(size(Ψ)) LSQ system using QR factorisation")
   qrΨ = qrfact(Ψ)
   verbose && @show cond(qrΨ[:R])
   c = qrΨ \ Y

   if verbose
      rel_rms = norm(Ψ * c - Y) / norm(Y)
      verbose && println("Relative RMSE on training set: ", rel_rms)
   end

   if E0 != 0
      basis = [ OneBody(E0); db.basis[Ibasis] ]
      c = [1.0; c]
   else
      basis = db.basis[Ibasis]
   end

   return NBodyIP(basis, c)
end


"""
`get_lsq_system(lsq; kwargs...) -> Ψ, Y, Ibasis`

Assemble the least squares system + rhs. The `kwargs` can be used to
select a subset of the available data or basis set, and to adjust the
weights by config_type. For more complex weight adjustments, one
can directly modify the `lsq.data[n].w` coefficients.

## Keyword Arguments:

* weights: A dictionary specifying the
```
weights = Dict("E" => 100.0, "F" => 1.0, "V" => 0.01)
```

* config_weights: a tuple of string, value pairs, e.g.,
```
configweights = Dict("solid" => 10.0, "liquid" => 0.1)
```
this adjusts the weights on individual configurations from these categories
if no weight is provided then the weight provided with the is used.
Note in particular that `config_weights` takes precedence of Dat.w!
If a weight 0.0 is used, then those configurations are removed from the LSQ
system.
"""


# # assemble crude error tables and scatter plots
# include("errors.jl")


end
