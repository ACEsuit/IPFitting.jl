
using NBodyIPs, JuLIP, Test, IPFitting, DataFrames
using JuLIP.Potentials: evaluate_d
using IPFitting: Dat, LsqDB
using NBodyIPs: BondLengthDesc
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
calc = let r0=r0
   LennardJones(r0=r0) * C2Shift(2.5*r0)
end
data = generate_data(:Cu, 3, 0.25*r0, 100, calc)
rcut2 = cutoff(calc)
D2 = BondLengthDesc(PolyTransform(1, r0), CosCut(rcut2-1, rcut2))

##
degrees = [4, 6, 8, 10, 12, 14, 16]

for solve_met in [(:qr,), (:svd, 2), (:rrqr, 1e-14)]
   global err_eunif = Float64[]
   global err_funif = Float64[]
   global err_erms = Float64[]
   global err_frms = Float64[]

   ptrain = 0.8

   for d in degrees
      rr = range(0.9*r0, stop=cutoff(calc), length=200)
      global err_eunif
      global err_funif
      global err_erms
      global err_frms
      global degrees
      local B2
      local db
      B2 = nbpolys(2, D2, d)
      @show length(B2)
      db = LsqDB("", B2, data)
      Itrain, Itest = splittraintest(db)
      @show Itrain
      @show Itest
      @show length(Itrain), length(Itest)
      @test isempty(intersect(Itrain, Itest))
      @test sort(union(Itrain, Itest)) == 1:length(db.configs)

      IP, fitinfo = Lsq.lsqfit(db, E0 = 0.0,
                               Itrain = Itrain,
                               Itest = Itest,
                               configweights = Dict("rand" => 1.0),
                               obsweights   = Dict("E" => 100.0, "F" => 1.0),
                               solver = solve_met,
                               combineIP = NBodyIP )
      @info("done fitting...")
      errs = fitinfo["errors"]
      errs_test = fitinfo["errtest"]
      @info("done assembling errors")
      @info("Training Errors")
      Err.rmse_table(rmse(errs)...)
      @info("Test Errors")
      Err.rmse_table(rmse(errs_test)...)
      V2 = IP.components[2]
      ev2 = norm(V2.(rr) - 0.5 * calc.(rr), Inf)
      dev2 = norm(evaluate_d.(Ref(V2), rr) - 0.5 * evaluate_d.(Ref(calc), rr), Inf)
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
   (@test minimum(err_erms) < 1e-4) |> println
   (@test minimum(err_frms) < 1e-3) |> println
   (@test minimum(err_eunif) < 3e-4) |> println
   (@test minimum(err_funif) < 4e-4) |> println
end
