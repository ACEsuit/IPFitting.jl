
import Base: ==, convert, vec

import JuLIP
using JuLIP: Atoms, JVec, JMat, JVecs, AbstractCalculator, mat, vecs,
             numbers, positions, cell, pbc

export LsqDB, Dat, basis, configs, config


"""
`Dat`: store one atomistic configuration and the associated
observations (e.g. E, F, V, ...)
"""
mutable struct Dat
   at::Atoms                         # configuration
   configtype::String                # group identifier
   D::Dict{String, Vector{Float64}}  # list of observations
   info::Dict{String, Any}           # anything else...
end

# TODO: should D.info also be compared?
==(d1::Dat, d2::Dat) = (
      (d1.configtype == d2.configtype) && (d1.D == d2.D) &&
      all( f(d1.at) == f(d2.at) for f in (positions, numbers, cell, pbc) )
   )

function Dat(at::Atoms, config_type::AbstractString; kwargs...)
   dat = Dat(at, config_type, Dict{String, Vector{Float64}}(), Dict{String, Any}())
   for (key, val) in kwargs
      str_key = string(key)
      if !(ismissing(val) || (val == nothing))
         # convert the observation to a basic vector and store it
         dat.D[str_key] = vec(str_key, val)
      end
   end
   return dat
end

# TODO: using `JuLIP.Dict(::Atoms...)` and
# enable storage of the complete Atoms object,
# which may contain additional information

Base.Dict(d::Dat) =
   Dict("__id__" => "NBodyIPFitting.Dat",
         "X" => positions(d.at) |> mat,
         "Z" => numbers(d.at),
         "cell" => cell(d.at) |> Matrix,
         "pbc" => Int.([pbc(d.at)...]),
         "configtype" => d.configtype,
         "D" => d.D )


function Dat(D::Dict)
   X = JuLIP._auto_X(D["X"])
   at = Atoms( X = JuLIP._auto_X(D["X"]),
               Z = Vector{Int}(D["Z"]),
               cell = JuLIP._auto_cell(D["cell"]),
               pbc = JuLIP._auto_pbc(D["pbc"]) )
   return Dat(at, D["configtype"], Dict{String, Any}(D["D"]), Dict{String, Any}())
end

convert(::Val{Symbol("NBodyIPFitting.Dat")}, D::Dict) = Dat(D)

# -----------------------------------------------------------------


"""
`mutable struct LsqDB{TD, TB}`

A representation of a least-squares system stored on disk, which can
be extended by adding data, or basis functions
"""
mutable struct LsqDB
   basis::Vector{AbstractCalculator}
   configs::Vector{Dat}
   Ψ::Matrix{Float64}
   dbpath::String
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
vec(::Val{:F}, F) = mat(F)[:]
```
or equivalently, `vec("F", F)`
"""
function vec_obs end

"""
convert a Vector{T} to some real (atomistic) data, e.g.,
```
x::Vector{Float64}
devec(::Val{:F}, x) = vecs( resize(x, 3, length(x) ÷ 3) )
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
