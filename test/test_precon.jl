

module Preconditioners

using JuLIP, LinearAlgebra

@doc raw"""
`AdjacancyPrecon` :

If `(r0, r1, k)` is in the `intervals` list and if
`r0 <= rij <= r1` then the 3 x 3 block
```math
   P_{ij} = k \big( (1-c) \hat{R} \otimes \hat{R} + c I \big),
```
where ``\hat{R}`` is the bond direction and `c = innerstab`.
* Choose `innerstab = 0` to get full infinitesimal rotation-invariance
* Choose `innerstab = 1` to get full resistance against infinitesimal rotations
* In optimisation we often found that `innerstab = 0.1` is a decent choice.

In addition the `stab` parameter adds a global stabilisation shifting
the entire matrix by `stab * I`.
"""
struct AdjacancyPrecon
   intervals::Vector{NamedTuple{(:r0, :r1, :k), Tuple{Float64, Float64, Float64}}}
   stab::Float64
   innerstab::Float64
end

import JuLIP.FIO: read_dict, write_dict


AdjacancyPrecon(intervals; stab = 1e-4, innerstab = 0.1) =
   AdjacancyPrecon(intervals, stab, innerstab)

_cutoff(P::AdjacancyPrecon) =
   maximum( i.r1 for i in P.intervals )

function (P::AdjacancyPrecon)(R::JVecF)
   r = norm(R)
   R̂ = R / r
   Π = (1 - P.innerstab) * R̂ * R̂' + P.innerstab * one(JMat{Float64})
   for i in P.intervals
      if i.r0 <= r <= i.r1
         return i.k * Π
      end
   end
   return 0.0 * Π
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

#---
# CAS - please convert this into a test set:

@info("Test evaluation and some basic properties of AdjacancyPrecon")

using JuLIP, Test

at = bulk(:W, cubic=true) * 2

intervals = [ (r0 = 2.0, r1 = 3.5, k = 1.0),
              (r0 = 3.5, r1 = 4.5, k = 0.1) ]

# the block-Identity preconditoner
adjp1 =  Preconditioners.AdjacancyPrecon(intervals; innerstab = 1.0)
# a variation that accounts for some infinitesimal rotational invariance
# the parameter innerstab can interpolate between the two.
adjp2 =  Preconditioners.AdjacancyPrecon(intervals; innerstab = 0.1)  # <= DEFAULT
adjp3 =  Preconditioners.AdjacancyPrecon(intervals; innerstab = 0.0)

for adjp in [adjp1, adjp2, adjp3]
   P = adjp(at)
   println(@test P' ≈ P)
   sumP = sum(P, dims=1)
   println(@test all(sumP .≈ adjp.stab))
   sumP = sum(P, dims=2)
   println(@test norm(sumP .- adjp.stab, Inf) < 1e-12)
   spec = sort(eigvals(P))
   println(@test all(spec[1:3] .≈ adjp.stab))
end

#---
# # how to use:
# weights = Dict()
# weights["precon"] = Dict("default" => adjp)
