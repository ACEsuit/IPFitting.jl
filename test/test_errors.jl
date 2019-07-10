
using JuLIP, IPFitting, SHIPs, Printf
using IPFitting: Dat, LsqDB
using JuLIP.MLIPs: IPSuperBasis
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
      # V = virial(calc, at)
      push!(data, Dat(at, cn; E = E, F = F)) # , V = V))
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
b3basis(deg) = IPSuperBasis(
      PairBasis(deg, PolyTransform(2, r0), 2, cutoff(calc)),
      SHIPBasis(TotalDegree(deg, 1.5), 2, PolyTransform(3, r0), 2, 0.5*r0, cutoff(calc))
   )
B = b3basis(10)
@show length(B)
db = LsqDB("", B, data)
IP, fitinfo = lsqfit( db,
                   E0 = 0.0,
                   configweights = Dict("rand1" => 1.0, "rand2" => 0.5),
                   obsweights   = Dict("E" => 100.0, "F" => 1.0),
                   verbose=true,
                   solver = (:rrqr, 1e-7) )
# note we are using RRQR here to make sure the fit is well-conditioned!

# IPf = fast(IP)
IPf = IP

##
@info("Checking that the two rmse computations agree")

add_fits!(IPf, data)
errs, errsrel = rmse(data)
rmse_table(errs, errsrel)

olderrs = fitinfo["errors"]
rmse_table(rmse(olderrs)...)

# @printf(" cfgtype obstype | abs errs            |  rel errs \n")
# @printf("                 |    new        old   |    new        old   \n")
for cn in keys(errs)
   for ot in keys(errs[cn])
      # @printf(" %6s %5s    | %.3e %.3e | %.3e %.3e \n",
      #          cn, ot,
      #          errs[cn][ot], olderrs["rmse"][cn][ot],
      #          errsrel[cn][ot], olderrs["relrmse"][cn][ot] )
      println(@test errs[cn][ot] ≈ olderrs["rmse"][cn][ot])
      println(@test errsrel[cn][ot] ≈ olderrs["relrmse"][cn][ot])
   end
end
