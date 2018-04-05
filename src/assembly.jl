

using JuLIP, NeighbourLists, StaticArrays

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

Base.length(V::NBodies) = length(V.c)

cutoff(V::NBodies) = V.cutoff

function evaluate(V::NBodies{N, T}, at::Atoms{T}) where {N, T}
   E = 0.0
   temp = zeros(T, length(at))
   for n = 1:length(V)
      nlist = neighbourlist(at, cutoff(V))::PairList
      fill!(temp, 0.0)
      maptosites!(V.f[n], temp, nbodies(N, nlist))
      E += c[n] * sum_kbn(temp)
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
      dE .+= c[n] .* temp
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

NBodyIP(args...) = NBodyIP( [args...] )
cutoff(V::NBodyIP) = maximum( cutoff.(V.orders) )
energy(V::NBodyIP, at::Atoms) = sum( Vn(at)  for Vn in V.orders )
forces(V::NBodyIP, at::Atoms) = - sum( (@D Vn(at))  for Vn in V.orders )
