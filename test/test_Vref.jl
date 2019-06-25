# TODO: review and rewrite this!

using NBodyIPs, JuLIP, Test, IPFitting, DataFrames
using JuLIP.Potentials: evaluate_d
using IPFitting: Dat, LsqDB
using NBodyIPs: BondLengthDesc, BondAngleDesc, bodyorder
Lsq = IPFitting.Lsq
Err = IPFitting.Errors
using LinearAlgebra: norm

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
data = generate_data(:Si, 2, 0.33*r0, 100, calc)

TRANSFORM = "exp( - 2 * (r/$r0 - 1.5) )"

rcut2 = cutoff(calc)*2.0
CUTOFF2 = PolyCutSym(2, rcut2)
D2 = BondLengthDesc(TRANSFORM, CUTOFF2)
rcut3 = cutoff(calc)*2.0
CUTOFF3 = PolyCut(2, rcut3)
D3 = BondLengthDesc(TRANSFORM, CUTOFF3)

B = [nbpolys(2, D2, 6); nbpolys(3, D3, 6)]

db = LsqDB("", B, data)

Ib2 = findall(b -> bodyorder(b) == 2, db.basis)
Ib3 = findall(b -> bodyorder(b) == 3, db.basis)

IP2, _ = Lsq.lsqfit( db,
                     E0 = 0.0,
                     configweights = Dict("rand" => 1.0),
                     obsweights   = Dict("E" => 30.0, "F" => 1.0),
                     combineIP = NBodyIP,
                     Ibasis = Ib2 )

# IP3a is the standard fit with Vref which is evaluated during Lsq assembly
IP3a, _ = Lsq.lsqfit( db,
                      Vref = IP2,
                      configweights = Dict("rand" => 1.0),
                      obsweights   = Dict("E" => 30.0, "F" => 1.0),
                      combineIP = NBodyIP,
                      Ibasis = Ib3 )


# precompute the observations for Vref in parallel
# (use JULIA_NUM_THREADS to use multi-threading for this!)
IP2f = fast(IP2)
Err.add_fits!(IP2f, db.configs, fitkey = "Vref")

# in the following comment Vref = IP2 is used only to construct the
# final potential but it will ot be evaluated, instead the precomputed
# observations from the `add_fits!` line are used.
IP3b, _ = Lsq.lsqfit( db,
                      Vref = IP2,
                      configweights = Dict("rand" => 1.0),
                      obsweights   = Dict("E" => 30.0, "F" => 1.0),
                      combineIP = NBodyIP,
                      Ibasis = Ib3 )

println(@test IP3a.components[1] == IP3b.components[1])
println(@test IP3a.components[2] === IP3b.components[2])
println(@test IP3a.components[3].c ≈ IP3b.components[3].c)
println(@test !(IP3a.components[3] === IP3b.components[3]))
