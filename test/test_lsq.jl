using ACE1, JuLIP, Test, IPFitting
using JuLIP.Potentials: evaluate_d
using IPFitting: Dat, LsqDB, observation, vec_obs
using LinearAlgebra: qr, norm

# generate random data
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

r0 = rnn(:Cu)
calc = let r0=r0
   LennardJones(r0=r0) * C2Shift(2.5*r0)
end
data = generate_data(:Cu, 3, 0.33*r0, 50, calc)
rcut2 = cutoff(calc)

B2 = rpi_basis(species = :Cu, N = 2, maxdeg = 6,
   r0 = r0, rcut = cutoff(calc), rin = 0.5*r0,
   pin = 2)

# generate a database so that we generate the row indices
db = LsqDB("", B2, data)

## generate a LSQ system by hand
atlen = length(data[1])
len = 1 + length(observation(data[1], "F"))
Y_man = zeros(length(data) * len)
Ψ_man = zeros(length(data) * len, length(B2))
for (okey, evalb, w) in zip(["E", "F"], (energy, forces), (1.345/sqrt(atlen), 2.987))
   for d in data
      o = observation(d, okey)
      irows = IPFitting.DB.matrows(d, okey)
      leno = length(o)
      Y_man[irows] = w * o
      #for colidx = 1:length(B2)
      Ψ_man[irows, :] = w * vec_obs(okey, evalb(B2, d.at))
      #end
   end
end

## test that Ψ_man, F_man are equivalent to a system generated from a db

Ψ, Y = IPFitting.Lsq.get_lsq_system( db,
               Vref = OneBody(:Cu => 0.0),
               weights = Dict("default" => Dict("E" => 1.345, "F" => 2.987)),
            )

##
print_tf(@test Y ≈ Y_man)
print_tf(@test Ψ ≈ Ψ_man)
println() 