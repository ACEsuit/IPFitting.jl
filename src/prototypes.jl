
import Base: ==, convert

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
config_type(d)
length(d)    # number of atoms
```
If information is missing, the relevant function will return `nothing` instead
(TODO for J v1.0: change `nothing` to `missing`)
"""
mutable struct Dat
   at::Atoms
   config_type::Union{Void, String}
   D::Dict{String, Any}
end

==(d1::Dat, d2::Dat) = (
      (d1.config_type == d2.config_type) && (d1.D == d2.D) &&
      all( f(d1.at) == f(d2.at) for f in (positions, numbers, cell, pbc) )
   )

function Dat(at::Atoms, config_type::AbstractString; kwargs...)
   dat = Dat(at, config_type, Dict{String, Any}())
   for (key, val) in kwargs
      # interpret `nothing` as `missing`
      if val != nothing
         dat.D[string(key)] = val
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
         "config_type" => d.config_type,
         "D" => d.D )

function Dat(D::Dict)
   at = Atoms( X = D["X"] |> vecs,
               Z = D["Z"],
               cell = JMat(D["cell"]),
               pbc = tuple(Bool.(D["pbc"])...) )
   return Dat(at, D["config_type"], Dict{String, Any}(D["D"]))
end

convert(::Val{Symbol("NBodyIPFitting.Dat")}, D::Dict) = Dat(D)


"""
`mutable struct LsqDB{TD, TB}`

A representation of a least-squares system stored on disk, which can
be extended by adding data, or
"""
mutable struct LsqDB
   data::Vector{Dat}
   basis::Vector{AbstractCalculator}
   dirname::String
end


# prototypes

import Base.vec

"""
convert some real data, in some generic format, into a vector to be stored
in a `Dat` or Lsq system. E.g.,
```
F = forces(...)::Vector{JVecF}
vec(::Val{:F}, F) = mat(F)[:]
```
"""
vec

"""
convert a Vector{T} to some real (atomistic) data, e.g.,
```
x::Vector{Float64}
devec(::Val{:F}, x) = vecs( resize(x, 3, length(x) ÷ 3) )
```
"""
function devec end




# """
# `mutable struct LsqSys`: type storing all information to perform a
# LSQ fit for an interatomic potential. To assemble the LsqSys use
# ```
# kron(data, basis)
# LsqSys(data, basis)
# ```
# """
# mutable struct LsqSys{TD}
#    data::Vector{TD}
#    basis::Vector{AbstractCalculator}
#    Iord::Vector{Vector{Int}}     # result of split_basis
#    Ψ::Matrix{Float64}
# end
