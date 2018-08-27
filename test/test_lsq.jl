
using NBodyIPs, JuLIP, Base.Test, NBodyIPFitting
using JuLIP.Potentials: evaluate_d
using NBodyIPFitting: Dat, LsqDB
using NBodyIPFitting.Data: observation
using NBodyIPs: BondLengthDesc, nbpolys
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

r0 = rnn(:Cu)
calc = let r0=r0
   LennardJones(r0=r0) * C2Shift(2.5*r0)
end
data = generate_data(:Cu, 3, 0.25*r0, 30, calc)
rcut2 = cutoff(calc)
D2 = BondLengthDesc("exp( - 2 * (r/$r0 - 1) )", (:cos, rcut2-1, rcut2))

B2 = nbpolys(2, D2, 6)

## generate a LSQ system by hand
len = 1 + length(observation(data[1], "F"))
Y = zeros(length(data) * len)
Ψ = zeros(length(data) * len, length(B2))
rowidx = 0
for (okey, evalb) in zip(["E", "F"], (energy, forces))
   for d in data
      o = observation(d, okey)
      leno = length(o)
      Y[(rowidx+1):(rowidx+leno)] = o
      for colidx = 1:length(B2)
         Ψ[(rowidx+1):(rowidx+leno), colidx] = vec(Val(Symbol(okey)), evalb(B2[colidx], d.at))
      end
      rowidx += leno
   end
end
## test that Ψ_man, F_man are equivalent to a system generated from a db
Ψ_man, Y_man = Ψ, Y
db = LsqDB("", B2, data)
Ψ, Y = NBodyIPFitting.Lsq.get_lsq_system( db,
               E0 = 0.0,
               configweights = Dict("rand" => 1.0),
               dataweights = Dict("E" => 1.0, "F" => 1.0)
            )

println(@test Y ≈ Y_man)
println(@test Ψ ≈ Ψ_man)

# TODO: this could be expanded to add multiple basis function types and
#       multiple configuratin types
