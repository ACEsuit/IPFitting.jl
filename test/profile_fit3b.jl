
using NBodyIPs, JuLIP, Base.Test, NBodyIPFitting
using JuLIP.Potentials: evaluate_d
using NBodyIPFitting: Dat, LsqDB
using NBodyIPs: BondLengthDesc
const Lsq = NBodyIPFitting.Lsq

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

r0 = rnn(:Si)
calc = StillingerWeber()
data = generate_data(:Si, 2, 0.2*r0, 30, calc)

TRANSFORM = "exp( - 2 * (r/$r0 - 1) )"

rcut2 = 2 * cutoff(calc)
CUTOFF2 = "(:cos, $(rcut2-1), $rcut2)"
D2 = BondLengthDesc(TRANSFORM, CUTOFF2)

rcut3 = 2 * cutoff(calc)
CUTOFF3 = "(:cos, $(rcut3-1), $rcut3)"
D3 = BondLengthDesc(TRANSFORM, CUTOFF3)

info("Profiling pure 3B basis assembly")
deg3 = 10
B = nbpolys(3, D3, deg3)
@show length(B)
LsqDB("", B, data[1:2]; verbose=false)
@time LsqDB("", B, data)

info("Profiling combined 2B+3B basis assembly")
deg2, deg3 = 10, 10
B = [nbpolys(2, D2, deg2); nbpolys(3, D3, deg3)]
@show length(B)
LsqDB("", B, data[1:2], verbose=false)
@time LsqDB("", B, data)
