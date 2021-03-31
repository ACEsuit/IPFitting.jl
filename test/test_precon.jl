

module Preconditioners

using JuLIP, LinearAlgebra

# @doc raw"""
# `AdjacancyPrecon`
# ```math
#    P_{ij} = k(n) \qquad \text{if} \quad
#    r_0(n) \leq r_{ij} \leq r_1(n)
# ```
# """
struct AdjacancyPrecon
   intervals::Vector{NamedTuple{(:r0, :r1, :k), Tuple{Float64, Float64, Float64}}}
   stab::Float64
   infrot::Bool
end

import JuLIP.FIO: read_dict, write_dict



AdjacancyPrecon(intervals; stab = 1e-4, infrot=false) =
   AdjacancyPrecon(intervals, stab, infrot)

_cutoff(P::AdjacancyPrecon) =
   maximum( i.r1 for i in P.intervals )

function (P::AdjacancyPrecon)(R::JVecF)
   r = norm(R)
   if !P.infrot
      for i in P.intervals
         if i.r0 <= r <= i.r1
            return i.k * Matrix(I, (3,3))
         end
      end
   else
      error("infrot not implemented yet")
   end
   return 0.0
end

function (precon::AdjacancyPrecon)(at::Atoms)
   P = zeros(3*length(at), 3*length(at))
   nlist = neighbourlist(at, _cutoff(precon))
   for i = 1:length(at)
      Js, Rs = neigs(nlist, i)
      for (j, Rij) in zip(Js, Rs)
         if j < i; continue; end
         Pij = precon(Rij)
         inds = 3*(i-1) .+ (1:3)
         jnds = 3*(j-1) .+ (1:3)
         for a = 1:3, b = 1:3
            P[inds[a], jnds[b]] -= Pij[a,b]
            P[jnds[b], inds[a]] -= Pij[a,b]
            P[inds[a], inds[b]] += Pij[a,b]
            P[jnds[a], jnds[b]] += Pij[a,b]
         end
      end
   end
   P += precon.stab * I
   return P
end

end

#---
using JuLIP

#---

at = bulk(:W, cubic=true) * 2

adjp = Preconditioners.AdjacancyPrecon( [ (r0 = 2.0, r1 = 3.5, k = 1.0),
                                          (r0 = 3.5, r1 = 4.5, k = 0.1) ] )
P = adjp(at)

P' ≈ P
eigvals(P)

weights = Dict()
weights["precon"] = Dict("default" => adjp)

#---


using LinearAlgebra
A = rand(10,10)
P = A' * A

rtP1 = sqrt(P)
rtP1 * rtP1 ≈ P

rtP2 = cholesky(P, Val(false)).L
pinv(rtP2') * pinv(rtP2) ≈ pinv(P)


B = rand(10, 100)
pinv(rtP2) * B ≈ rtP2 \ B


Matrix(I, (3,3))
