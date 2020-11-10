
module Aux

using LinearAlgebra: norm, pinv
using IPFitting: Dat
using JuLIP

function rdf(cfgs::AbstractVector{<: Dat}, rcut=5.0)
   rs = Float64[]
   for cfg in cfgs
      at = cfg.at
      for (i,j,r) in pairs(at, rcut)
         push!(rs, norm(r))
      end
   end
   return rs
end

function binding_energy(cfgs::AbstractVector{<: Dat}, E0, nneigs)
   Es = Float64[]
   for cfg in cfgs
      Nat = length(cfg.at)
      E = cfg.D["E"][1]
      push!(Es, (E - Nat * E0) / Nat)
   end
   return sum(Es) / length(cfgs) / nneigs
end

function extrap_grade(at, cfgs::AbstractVector{<: Dat}, basis)
   A = zeros(length(cfgs), length(basis))
   for (i,cfg) in enumerate(cfgs)
      A[i,:] = energy(basis, cfg.at)
   end
   C = energy(basis, at)
   return maximum(abs.(C' * pinv(A)))
end

end