
using JuLIP: Atoms, JVec, JMat, JVecs, AbstractCalculator

export LsqDB


"""
`Dat`: store one simulation data point. If `d::Dat`, to obtain the data, use
```
Atoms(d)
energy(d)
forces(d)
virial(d)
weight(d)
config_type(d)
length(d)    # number of atoms
```
If information is missing, the relevant function will return `nothing` instead
(TODO for J v1.0: change `nothing` to `missing`)
"""
mutable struct Dat{T}
   at::Atoms
   E::Union{Void, T}         # energy
   F::Union{Void, JVecs{T}}  # forces
   V::Union{Void, JMat{T}}   # stress
   w                         # weight, could be anything?
   config_type::Union{Void, String}
   D::Dict{String, Any}
end

Base.Dict(d::Dat) =
   Dict("id" => "NBodyIPFitting.Dat",
         "at" => Dict(d.at), "E" => d.E, "F" => mat(d.F), "V" => Matrix(d.V),
         "w" => d.w, "config_type" => d.config_type, "D" => d.D)

function Dat(D::Dict)
   at = Atoms(D["at"])
   E = D["E"]::Union{Void, Float64}
   F = D["F"] == nothing ? nothing : vecs(Matrix{Float64}(D["F"]))
   V = D["V"] == nothing ? nothing : JMat(Matrix{Float64}(D["V"]))
   w =
   return Dat(at, E, F, V, D["w"], D["config_type"], Dict{String, Any}(D["D"]))
end

"""
`mutable struct LsqSys`: type storing all information to perform a
LSQ fit for an interatomic potential. To assemble the LsqSys use
```
kron(data, basis)
LsqSys(data, basis)
```
"""
mutable struct LsqSys{TD}
   data::Vector{TD}
   basis::Vector{AbstractCalculator}
   Iord::Vector{Vector{Int}}     # result of split_basis
   Ψ::Matrix{Float64}
end


"""
`mutable struct LsqDB{TD, TB}`

A representation of a least-squares system stored on disk, which can
be extended by adding data, or
"""
mutable struct LsqDB{TD}
   data::Vector{TD}
   basis::Vector{AbstractCalculator}
   dirname::String
end
