
using JuLIP, IPFitting, ACE1, Printf
using IPFitting: Dat, LsqDB
using JuLIP.MLIPs: IPSuperBasis
using JuLIP.Testing: print_tf 
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
b3basis(deg) = rpi_basis(species = :Si, N = 2, maxdeg = 10,
                         r0 = r0, rcut = cutoff(calc), rin = 0.5*r0,
                         pin = 2)
B = b3basis(10)
@show length(B)
db = LsqDB("", B, data)
IP, fitinfo = lsqfit( db,
                      Vref = OneBody(:Si => 0.0),
                      weights = Dict("default" => Dict("E"=>100.0, "F"=>1.0),
                                       "rand2" => Dict("E"=>50.0, "F"=>0.5)),
                      verbose=true,
                      solver = Dict("solver" => :rrqr, "rrqr_tol" => 1e-5),
                      error_table=true )
# note we are using RRQR here to make sure the fit is well-conditioned!

# IPf = fast(IP)

##
@info("Checking that the two rmse computations agree")

add_fits!(IP, data)
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
      print_tf(@test errs[cn][ot] ≈ olderrs["rmse"][cn][ot])
      print_tf(@test errsrel[cn][ot] ≈ olderrs["relrmse"][cn][ot])
   end
end
