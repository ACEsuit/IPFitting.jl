
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

function extrap_grade(test, train::AbstractVector{<: Dat}, basis)
   A = zeros(length(train), length(basis))
   for (i,cfg) in enumerate(train)
      A[i,:] = energy(basis, cfg.at)
   end
   A_inv = pinv(A, rtol=1e-10)
   extrap_grades = zeros(length(test))
   for (i,cfg) in enumerate(test)
      C = energy(basis, cfg.at)
      extrap_grades[i] = maximum(abs.(C' * A_inv))
   end
   return extrap_grades
end

end