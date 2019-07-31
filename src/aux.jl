
module Aux

using LinearAlgebra: norm
using IPFitting: Dat

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

end
