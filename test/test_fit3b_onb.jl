
using NBodyIPs, JuLIP, Test, NBodyIPFitting, DataFrames
using JuLIP.Potentials: evaluate_d
using NBodyIPFitting: Dat, LsqDB
using NBodyIPs: BondLengthDesc, BondAngleDesc
using NBodyIPs.Polys: NBPoly

const Lsq = NBodyIPFitting.Lsq
const Err = NBodyIPFitting.Errors

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
data = generate_data(:Si, 2, 0.25*r0, 100, calc)

TRANSFORM = "exp( - 2 * (r/$r0 - 1.5) )"

rcut2 = cutoff(calc)*1.6
CUTOFF2 = "(:cos, $(rcut2-1), $rcut2)"
D2 = BondLengthDesc(TRANSFORM, CUTOFF2)

rcut3 = cutoff(calc)*1.6
CUTOFF3 = "(:cos, $(rcut3-1), $rcut3)"
D3 = BondLengthDesc(TRANSFORM, CUTOFF3)

cw = Dict("rand" => 1.0)
dw   = Dict("E" => 100.0, "F" => 1.0)

## 

# create a monomial database
B = nbpolys(3, D3, 8)
@show length(B)
db = LsqDB("", B, data)

# orthogonalise the basis
Bo = Lsq.onb(db, E0 = 0.0, configweights = cw, dataweights = dw)
# extract the pure 3-body terms
Bo = [ b.components[1] for b in Bo ]

# checkout how "nice" the coefficients are
@show [ norm(b.c) for b in Bo ]

# construct a new database, based on the data-ONB
db_onb = LsqDB("", Bo, data)

# fit with ONB
IP_onb, errs = Lsq.lsqfit( db_onb,
      E0 = 0.0,
      configweights = Dict("rand" => 1.0),
      dataweights   = Dict("E" => 100.0, "F" => 1.0)
   )

# fit the monomial basis
IP, errs = Lsq.lsqfit( db,
      E0 = 0.0,
      configweights = Dict("rand" => 1.0),
      dataweights   = Dict("E" => 100.0, "F" => 1.0)
   )

# check for errors
@show maximum( abs(energy(IP_onb, d.at) - energy(IP, d.at))  for d in data  )
