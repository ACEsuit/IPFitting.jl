
import Base: ==, convert, vec

import JuLIP
using JuLIP: Atoms, JVec, JMat, JVecs, AbstractCalculator, mat, vecs,
             numbers, positions, cell, pbc

export LsqDB, Dat


"""
`Dat`: store one simulation data point. If `d::Dat`, to obtain the data, use
```
Atoms(d)
energy(d)
forces(d)
virial(d)
# weight(d)
configtype(d)
length(d)    # number of atoms
```
If information is missing, the relevant function will return `nothing` instead
(TODO for J v1.0: change `nothing` to `missing`)
"""
mutable struct Dat
   at::Atoms
   configtype::String
   D::Dict{String, Vector{Float64}}
   info::Dict{String, Any}
end

# function convert(Dat, d::Any)
#    try
#       return Dat(d.at, d.configtype, d.D, Dict{String, Any}())
#    catch
#       @show d
#       error("convert of `d` unsuccesful; edit `prototypes.jl:convert()` to fix this")
#    end
# end

==(d1::Dat, d2::Dat) = (
      (d1.configtype == d2.configtype) && (d1.D == d2.D) &&
      all( f(d1.at) == f(d2.at) for f in (positions, numbers, cell, pbc) )
   )

function Dat(at::Atoms, config_type::AbstractString; kwargs...)
   dat = Dat(at, config_type, Dict{String, Vector{Float64}}(), Dict{String, Any}())
   for (key, val) in kwargs
      str_key = string(key)
      # interpret `nothing` as `missing`
      if val != nothing
         dat.D[str_key] = vec(Val(Symbol(str_key)), val)
      end
   end
   return dat
end

Base.Dict(d::Dat) =
   Dict("__id__" => "NBodyIPFitting.Dat",
         "X" => positions(d.at) |> mat,
         "Z" => numbers(d.at),
         "cell" => cell(d.at) |> Matrix,
         "pbc" => Int.([pbc(d.at)...]),
         "configtype" => d.configtype,
         "D" => d.D )


_read_cell(C::Matrix) = JuLIP.JMatF(C)
_read_cell(C::Vector{Any}) = JuLIP.JMatF( ([C[1] C[2] C[3]])... )

function Dat(D::Dict)
   at = Atoms( X = JuLIP._read_X(D["X"]),
               Z = Vector{Int}(D["Z"]),
               cell = _read_cell(D["cell"]),
               pbc = Bool.(D["pbc"]) )  # tuple(Bool.(D["pbc"])...)
   return Dat(at, D["configtype"], Dict{String, Any}(D["D"]), Dict{String, Any}())
end

convert(::Val{Symbol("NBodyIPFitting.Dat")}, D::Dict) = Dat(D)

# -----------------------------------------------------------------

const KronGroup = Dict{String, Array{Float64, 3}}
const DataGroup = Vector{Dat}

"""
`mutable struct LsqDB{TD, TB}`

A representation of a least-squares system stored on disk, which can
be extended by adding data, or basis functions
"""
mutable struct LsqDB
   basis::Vector{AbstractCalculator}
   data_groups::Dict{String, DataGroup}
   kron_groups::Dict{String, KronGroup}
   dbpath::String
end

basis(db::LsqDB) = db.basis
basis(db::LsqDB, i::Integer) = db.basis[i]
data(db::LsqDB) = db.data


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
vec(s::AbstractString, args...) = vec(Val(Symbol(s)), args...)

"""
convert a Vector{T} to some real (atomistic) data, e.g.,
```
x::Vector{Float64}
devec(::Val{:F}, x) = vecs( resize(x, 3, length(x) ÷ 3) )
```
"""
function devec end

function evaluate_lsq end

weighthook(::Val, d::Dat) = 1.0
weighthook(s::String, d::Dat) = weighthook(Val(Symbol(s)), d)
