

using JuLIP, NeighbourLists, StaticArrays, FunctionWrappers

import JuLIP.Potentials: evaluate, evaluate_d, cutoff, energy, forces
using JuLIP.Potentials: @pot

export NBody, NBodies, NBodyIP


"""
`NBody` : this turns a function f defined on an N-simplex into a
calculator that evaluates energies and gradients of atom configurations
"""
NBody(N::Integer, f, f_d, cutoff::T; wrap=true) where {T} =
      NBodies(N, T[1.0], [f], [f_d], cutoff; wrap=wrap)


@pot struct NBodies{N, T, TF, TG}
   _::Val{N}        # body order
   c::Vector{T}     # coefficients
   f::Vector{TF}    # basis functions
   f_d::Vector{TG}  # gradient of basis functions
   cutoff::T        # cut-off radius
end

"""
`NBodies` : A struct that collects several polynomials of the same
body-order together so that their energies can be assembled in one
single loop. It represents the N-body function
```
sum( c[n] * f[n](at)  for n = 1:length(c) )
```
"""
NBodies

NBodies(N::Integer, c, f, f_d, cutoff::T; wrap=true) where {T} =
   NBodies(Val(N), Val((N*(N-1))÷2), c, f, f_d, cutoff, wrap)

function NBodies(valN::Val{N}, ::Val{DIM}, c, fs, fs_d, cutoff, wrap) where {N, DIM}
   if wrap
      fs = [ FWrap{DIM, Float64}(f) for f in fs ]
      fs_d = [ GWrap{DIM, Float64}(f_d) for f_d in fs_d ]
   end
   return NBodies(valN, c, fs, fs_d, cutoff)
end


bodyorder(b::NBodies{N}) where {N} = N

Base.length(V::NBodies) = length(V.c)

cutoff(V::NBodies) = V.cutoff

# TODO: this evaluate should be optimised to make the loop over basis
#       functions INSIDE the maptosites!
function evaluate(V::NBodies{N, T}, at::Atoms{T}) where {N, T}
   E = 0.0
   temp = zeros(T, length(at))
   nlist = neighbourlist(at, cutoff(V))::PairList
   for n = 1:length(V)
      fill!(temp, 0.0)
      maptosites!(V.f[n], temp, nbodies(N, nlist))   # compute site energies
      E += V.c[n] * sum_kbn(temp)    # E = ∑ c[n] * basis[n](at)
   end
   return E
end

function evaluate_d(V::NBodies{N}, at::Atoms{T}) where {N, T}
   dE = zeros(JVec{T}, length(at))
   temp = zeros(JVec{T}, length(at))
   nlist = neighbourlist(at, cutoff(V))::PairList
   for n = 1:length(V)
      fill!(temp, zero(JVec{T}))
      maptosites_d!(V.f_d[n], temp, nbodies(N, nlist))
      dE .+= V.c[n] .* temp
   end
   return dE
end


"""
`NBodyIP` : wraps `NBodies` into a JuLIP calculator, defining
`energy`, `forces` and `cutoff`.
"""
struct NBodyIP <: AbstractCalculator
   orders::Vector{NBodies}
end

# NBodyIP(args...) = NBodyIP( [args...] )
cutoff(V::NBodyIP) = maximum( cutoff.(V.orders) )
energy(V::NBodyIP, at::Atoms) = sum( Vn(at)  for Vn in V.orders )
forces(V::NBodyIP, at::Atoms) = - sum( (@D Vn(at))  for Vn in V.orders )


"""
```
function NBodies( basis::AbstractVector{TB <: NBodies},
                  coeffs::AbstractVector{T <: AbstractFloat} )
```
convert a basis - i.e. a vector of length-1 `NBodies` (=an `NBody`) into a
single `NBodies term with coefficients.
"""
@noinline function NBodies(basis, coeffs)
   # check that what we are given here is really a meaningful basis with
   # a consistent body-order
   rcut = cutoff(basis[1])
   bo = bodyorder(basis[1])
   @assert all(length(b) == 1 for b in basis)
   @assert all(b.c[1] == 1.0 for b in basis)
   @assert all(cutoff(b) == rcut for b in basis)
   @assert all(bodyorder(b) == bo for b in basis)
   # collect into a single `NBodies`
   return NBodies(bo,
                  coeffs,
                  [b.f[1] for b in basis],
                  [b.f_d[1] for b in basis],
                  rcut)
end



function NBodyIP(basis, coeffs)
   orders = NBodies[]
   bos = bodyorder.(basis)
   for bo = 2:maximum(bos)
      Ibo = find(bos .== bo)
      push!(orders, NBodies(basis[Ibo], coeffs[Ibo]))
   end
   return NBodyIP(orders)
end
