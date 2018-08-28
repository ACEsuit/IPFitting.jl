
using NBodyIPs, JuLIP, Base.Test, NBodyIPFitting, DataFrames
using JuLIP.Potentials: evaluate_d
using NBodyIPFitting: Dat, LsqDB
using NBodyIPs: BondLengthDesc
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

r0 = rnn(:Cu)
calc = let r0=r0
   LennardJones(r0=r0) * C2Shift(2.5*r0)
end
data = generate_data(:Cu, 3, 0.25*r0, 30, calc)
rcut2 = cutoff(calc)
# D2 = BondLengthDesc("exp( - 2 * (r/$r0 - 1) )", (:cos, rcut2-1, rcut2))
D2 = BondLengthDesc("($r0)/r", (:cos, rcut2-1, rcut2))

##
err_eunif = Float64[]
err_funif = Float64[]
err_erms = Float64[]
err_frms = Float64[]
degrees = [4, 6, 8, 10, 12, 14, 16]
rr = linspace(0.9*r0, cutoff(calc), 200)
for d in degrees
   B2 = nbpolys(2, D2, d)
   @show length(B2)
   db = LsqDB("", B2, data)
   IP, errs = Lsq.lsqfit(db, E0 = 0.0,
                             configweights = Dict("rand" => 1.0),
                             dataweights   = Dict("E" => 100.0, "F" => 1.0) )
   Err.table_absolute(errs)
   V2 = IP.components[1]
   ev2 = vecnorm(V2.(rr) - 0.5 * calc.(rr), Inf)
   dev2 = vecnorm(evaluate_d.(V2, rr) - 0.5 * evaluate_d.(calc, rr), Inf)
   println("   V2 - uniform error = ", ev2, " | ", dev2)
   push!(err_eunif, ev2)
   push!(err_funif, dev2)
   push!(err_erms, Err.rmse(errs, "rand", "E"))
   push!(err_frms, Err.rmse(errs, "rand", "F"))
end



df = DataFrame( :degrees => degrees,
                :unif_E => err_eunif,
                :unif_F => err_funif,
                :rms_E => err_erms,
                :rms_F => err_frms )
display(df)

(@test minimum(err_erms) < 1e-4) |> println
(@test minimum(err_frms) < 1e-3) |> println
(@test minimum(err_eunif) < 2e-4) |> println
(@test minimum(err_funif) < 4e-4) |> println
