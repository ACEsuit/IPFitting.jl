


import JuLIP.Potentials: evaluate, evaluate_d
using JuLIP.Potentials: @pot

@pot struct NBody{N, T, TF, TG}
   _::Val{N}
   f::TF
   f_d::TG
   cutoff::T
end

NBody(N::Integer, f, f_d, cutoff) = NBody(Val(N), f, f_d, cutoff)
cutoff(V::NBody) = V.cutoff

function evaluate(V::NBody{N}, at::Atoms{T}) where {N, T}
   temp = zeros(T, length(at))
   nlist = neighbourlist(at, cutoff(V))::PairList
   maptosites!(V.f, temp, nbodies(N, nlist))
   return sum_kbn(temp)
end

function evaluate_d(V::NBody{N}, at::Atoms{T}) where {N, T}
   temp = zeros(JVec{T}, length(at))
   nlist = neighbourlist(at, cutoff(V))::PairList
   return maptosites_d!(V.f_d, temp, nbodies(N, nlist))
end
