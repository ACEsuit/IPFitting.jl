
using JuLIP, Test, IPFitting, DataFrames
using SHIPs
using JuLIP.Potentials: evaluate_d
using IPFitting: Dat, LsqDB
Lsq = IPFitting.Lsq
Err = IPFitting.Errors
import Random
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

Random.seed!(1)
r0 = rnn(:Cu)
z0 = atomic_number(:Cu)
calc = let r0=r0
   LennardJones(r0=r0) * C2Shift(2.5*r0)
end
data = generate_data(:Cu, 3, 0.25*r0, 70, calc)
rcut2 = cutoff(calc)

trans = PolyTransform(1, r0)
fcut = PolyCutoff2s(2, 0.5*r0, rcut2)
pairbasis(deg) = SHIPBasis( SparseSHIP(1, :Cu, deg, 1.0), trans, fcut)

##
degrees = [4, 7, 10, 13, 16, 19]

IPt = nothing

for solve_met in [(:qr,), (:rrqr, 1e-12)]
   global err_eunif = Float64[]
   global err_funif = Float64[]
   global err_erms = Float64[]
   global err_frms = Float64[]

   ptrain = 0.8

   for d in degrees
      rr = range(0.9*r0, cutoff(calc), length=200)
      global err_eunif
      global err_funif
      global err_erms
      global err_frms
      global degrees
      local B2
      local db
      B2 = pairbasis(d)
      @show length(B2)
      db = LsqDB("", B2, data)
      Itrain, Itest = splittraintest(db)
      @test isempty(intersect(Itrain, Itest))
      @test sort(union(Itrain, Itest)) == 1:length(db.configs)

      IP, fitinfo = Lsq.lsqfit(db, E0 = 0.0,
                               Itrain = Itrain,
                               Itest = Itest,
                               configweights = Dict("rand" => 1.0),
                               obsweights   = Dict("E" => 100.0, "F" => 1.0),
                               solver = solve_met )
      @info("done fitting...")
      errs = fitinfo["errors"]
      errs_test = fitinfo["errtest"]
      @info("done assembling errors")
      @info("Training Errors")
      Err.rmse_table(rmse(errs)...)
      @info("Test Errors")
      Err.rmse_table(rmse(errs_test)...)
      V2 = r -> JuLIP.Potentials.evaluate(IP, [r*JVec(1.0,0.0,0.0)], [z0,], z0)
      dV2 = r -> JuLIP.Potentials.evaluate_d(IP, [r*JVec(1.0,0.0,0.0)], [z0,], z0)[1][1]
      ev2 = norm(V2.(rr) - 0.5 * calc.(rr), Inf)
      dev2 = norm(dV2.(rr) - 0.5 * evaluate_d.(Ref(calc), rr), Inf)
      println("   V2 - uniform error = ", ev2, " | ", dev2)
      push!(err_eunif, ev2)
      push!(err_funif, dev2)
      push!(err_erms, errs["rmse"]["rand"]["E"])
      push!(err_frms, errs["rmse"]["rand"]["F"])
   end

   df = DataFrame( :degrees => degrees,
                          :unif_E => err_eunif,
                          :unif_F => err_funif,
                          :rms_E => err_erms,
                          :rms_F => err_frms )

   display(df)
   (@test minimum(err_erms) < 1e-6) |> println
   (@test minimum(err_frms) < 2e-4) |> println
   (@test minimum(err_eunif) < 2e-5) |> println
   (@test minimum(err_funif) < 1e-4) |> println
end
