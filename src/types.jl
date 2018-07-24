

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
   S::Union{Void, JMat{T}}   # stress
   w                         # weight, could be anything?
   config_type::Union{Void, String}
   D::Dict{String, Any}
end



"""
`mutable struct LsqSys`: type storing all information to perform a
LSQ fit for an interatomic potential. To assemble the LsqSys use
```
kron(data, basis)
LsqSys(data, basis)
```
"""
mutable struct LsqSys{TD, TB}
   data::Vector{TD}
   basis::Vector{TB}
   Iord::Vector{Vector{Int}}     # result of split_basis
   Ψ::Matrix{Float64}
end


mutable struct LsqDB

end
