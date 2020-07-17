
import Base: ==, convert, vec

import JuLIP
using JuLIP: Atoms, JVec, JMat, AbstractCalculator, mat, vecs,
             numbers, positions, cell, pbc

using JuLIP.MLIPs: IPBasis

import JuLIP: read_dict, write_dict

export LsqDB, Dat, basis, configs, config


"""
`Dat`: store one atomistic configuration and the associated
observations (e.g. E, F, V, ...)
"""
mutable struct Dat
   at::Atoms                         # configuration
   configtype::String                # group identifier
   D::Dict{String, Vector{Float64}}  # list of observations
   rows::Dict{String, Vector{Int}}   # row indices for LSQ system
   info::Dict{String, Any}           # anything else...
end

# should D.info also be compared?
==(d1::Dat, d2::Dat) = (
      (d1.configtype == d2.configtype) && (d1.D == d2.D) &&
      (d1.info == d2.info) &&
      all( f(d1.at) == f(d2.at) for f in (positions, numbers, cell, pbc) )
   )

function Dat(at::Atoms, config_type::AbstractString; kwargs...)
   dat = Dat(at, config_type, Dict{String, Vector{Float64}}(),
             Dict{String, Vector{Int}}(), Dict{String, Any}())
   for (key, val) in kwargs
      str_key = string(key)
      if !(ismissing(val) || (val == nothing))
         # convert the observation to a basic vector and store it
         dat.D[str_key] = vec_obs(str_key, val)
      end
   end
   return dat
end

write_dict(d::Dat) =
   Dict("__id__" => "IPFitting.Dat",
         "at" => write_dict(d.at),
         "configtype" => d.configtype,
         "D" => d.D,
         "rows" => d.rows,
         "info" => d.info)


function Dat(D::Dict)
   at = Atoms(D["at"])
   return Dat(at,
              D["configtype"],
              Dict{String, Any}(D["D"]),
              Dict{String, Any}(D["rows"]),
              Dict{String, Any}(D["info"]))
end

read_dict(::Val{Symbol("IPFitting.Dat")}, D::Dict) = Dat(D)

observation(d::Dat, key::String) = d.D[key]
hasobservation(d::Dat, key::String) = haskey(d.D, key)
observation(key::String, d::Dat) = d.D[key]


# -----------------------------------------------------------------


"""
`mutable struct LsqDB{TD, TB}`

A representation of a least-squares system stored on disk, which can
be extended by adding data, or basis functions
"""

mutable struct LsqDB
   basis::IPBasis
   configs::Vector{Dat}
   Ψ::Matrix{Float64}
   dbpath::String
   QR::Dict
end

basis(db::LsqDB) = db.basis
basis(db::LsqDB, i::Integer) = db.basis[i]
configs(db::LsqDB) = db.configs
config(db::LsqDB, i::Integer) = db.configs[i]


# prototypes

"""
convert some real data, in some generic format, into a vector to be stored
in a `Dat` or Lsq system. E.g.,
```
F = forces(...)::Vector{JVecF}
vec_obs(::Val{:F}, F) = mat(F)[:]
```
or equivalently, `vec_obs("F", F)`
"""
function vec_obs end

"""
convert a Vector{T} to some real (atomistic) data, e.g.,
```
x::Vector{Float64}
devec_obs(::Val{:F}, x) = vecs( resize(x, 3, length(x) ÷ 3) )
```
"""
function devec_obs end

"""
evaluate a specific observation type
"""
function eval_obs end

"""
create special weights for different observations
"""
weighthook(::Val, d::Dat) = 1.0
weighthook(s::String, d::Dat) = weighthook(Val(Symbol(s)), d)

"""
create special weights for different observations, similar to
`weighthooks` but this is used for computing errors rather than
for fitting
"""
err_weighthook(::Val, d::Dat) = 1.0
err_weighthook(s::String, d::Dat) = err_weighthook(Val(Symbol(s)), d)
