info("Loading libraries...")

using NBodyIPFitting.Data: observation, hasobservation, configname
using NBodyIPs, JuLIP, Base.Test, NBodyIPFitting, DataFrames, ASE
using JuLIP.Potentials: evaluate_d
using NBodyIPFitting: Dat, LsqDB
using NBodyIPs: BondLengthDesc, BondAngleDesc
using NBodyIPs: bodyorder, degree
using NBodyIPFitting
using NBodyIPFitting.Lsq
using Base.Test

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

rcut2 = cutoff(calc)
CUTOFF2 = "(:cos, $(rcut2-1), $rcut2)"
D2 = BondAngleDesc(TRANSFORM, CUTOFF2)

rcut3 = cutoff(calc)
CUTOFF3 = "(:cos, $(rcut3-1), $rcut3)"
D3 = BondAngleDesc(TRANSFORM, CUTOFF3)
basis = [nbpolys(2, D2, 8); nbpolys(3, D3, 10)]

println("Create a database.")
@show length(basis)
db = LsqDB("", basis, data)

configweights = Dict("rand" => 1.0)
IP, errs = Lsq.lsqfit( db,
                      E0 = 0.0,
                      configweights = configweights,
                      dataweights   = Dict("E" => 100.0, "F" => 1.0) )

res_dict = results_dict(data, IP; pathname = "")

errs1 = lsqerrors(db, [IP.components[1].c ; IP.components[2].c] ,
                  collect(1:length(basis)) ; E0 = 0.0)
errs2 = lsqerrors(res_dict, data; E0 = 0.0)

println("Check consistency of the errors")
for f in fieldnames(errs1)
    D1 = getfield(errs1, f)
    D2 = getfield(errs2, f)
    for cn in keys(configweights)
       for ot in ["E","F","V"]
          if !haskey(D1[cn],ot)
             continue
          end
          @test maximum(abs.(D1[cn][ot] - D2[cn][ot])) < 1e-6
       end
    end
end
