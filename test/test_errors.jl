
using JuLIP, NBodyIPs, IPFitting
using IPFitting: Dat, LsqDB
using NBodyIPs: bodyorder, degree
using Test
using LinearAlgebra: norm



# generate random data
function generate_data(species, L, rmax, N, calc; cn="rand")
   data = Dat[]
   for n = 1:N
      at = bulk(species; cubic=true, pbc=true) * L
      rattle!(at, rand() * rmax)
      E = energy(calc, at)
      F = forces(calc, at)
      V = virial(calc, at)
      push!(data, Dat(at, cn; E = E, F = F, V = V))
   end
   return data
end

##
r0 = rnn(:Si)
calc = StillingerWeber()
data1 = generate_data(:Si, 2, 0.33*r0, 20, calc; cn="rand1")
data2 = generate_data(:Si, 2, 0.1*r0, 20, calc; cn="rand2")
data = [data1; data2]

##
@info("generate a 3B fit to SW")
TRANSFORM = "exp( - 2 * (r/$r0 - 1.5) )"
rcut2 = cutoff(calc)*1.4
rcut3 = cutoff(calc)*1.9
CUTOFF2 = "(:cos, $(rcut2-1), $rcut2)"
CUTOFF3 = "(:cos, $(rcut3-1), $rcut3)"
D2 = BondLengthDesc(TRANSFORM, CUTOFF2)
D3 = BondLengthDesc(TRANSFORM, CUTOFF3)
B = [nbpolys(2, D3, 8); nbpolys(3, D3, 6)]
@show length(B)
db = LsqDB("", B, data)
IP, fitinfo = lsqfit( db,
                   E0 = 0.0,
                   configweights = Dict("rand1" => 1.0, "rand2" => 0.5),
                   obsweights   = Dict("E" => 100.0, "F" => 1.0),
                   combineIP = NBodyIP
                   )
IPf = fast(IP)

##
@info("Checking that the two rmse computations agree")

add_fits!(IPf, data)
errs, errsrel = rmse(data)
rmse_table(errs, errsrel)

olderrs = fitinfo["errors"]
rmse_table(rmse(olderrs)...)

for cn in keys(errs)
   for ot in keys(errs[cn])
      println(@test errs[cn][ot] ≈ olderrs["rmse"][cn][ot])
      println(@test errsrel[cn][ot] ≈ olderrs["relrmse"][cn][ot])
   end
end
