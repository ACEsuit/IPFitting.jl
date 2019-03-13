
using NBodyIPFitting, JuLIP, Random, Test
using NBodyIPFitting: observations, tfor_observations

function generate_data(species, L, rmax, N, calc)
   data = Dat[]
   for n = 1:N
      at = bulk(species; cubic=true, pbc=true) * L
      rattle!(at, rand() * rmax)
      E = energy(calc, at)
      F = forces(calc, at)
      push!(data, Dat(at, "rand"; E = E, F = F))
   end
   return data
end

Random.seed!(1)
r0 = rnn(:Cu)
calc = let r0=r0
   LennardJones(r0=r0) * C2Shift(2.5*r0)
end
data = generate_data(:Cu, 3, 0.25*r0, 20, calc)

for (okey, cfg, _) in observations(data)
   cfg.info["for_$okey"] = hash((hash(cfg.at), hash(okey)))
end

tfor_observations(data,
   (okey, cfg, _, lck) -> begin
      lock(lck)
      cfg.info["tfor_$okey"] = hash((hash(cfg.at), hash(okey)))
      unlock(lck)
   end)

for (okey, cfg, _) in observations(data)
   print(@test(cfg.info["for_$okey"] == cfg.info["tfor_$okey"]))
end
println()
