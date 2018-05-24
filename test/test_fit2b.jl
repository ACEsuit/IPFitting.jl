
using NBodyIPs, JuLIP, Base.Test

using JuLIP.Potentials: evaluate_d
using NBodyIPs.Data: Dat

# generate random data
function generate_data(species, L, rmax, N, calc)
   data = Dat[]
   for n = 1:N
      at = bulk(species; cubic=true, pbc=true) * L
      rattle!(at, rand() * rmax)
      E = energy(calc, at)
      F = forces(calc, at)
      push!(data, Dat(at, E, F, nothing))
   end
   return data
end

r0 = rnn(:Cu)
calc = LennardJones(r0=r0) * C2Shift(2.7*r0)
data = generate_data(:Cu, 3, 0.33, 200, calc)

B1 = [NBody(1.0)]
# TRANSFORM = (@analytic r -> (2.9/r)^3)
TRANSFORM = let r0=r0
   (@analytic r -> exp( - 2 * (r/r0 - 1) ))
end
rcut2 = cutoff(calc)
CUTOFF2 = (:cos, rcut2-1, rcut2)
D2 = Dictionary(TRANSFORM, CUTOFF2)

err_eunif = Float64[]
err_funif = Float64[]
err_erms = Float64[]
err_frms = Float64[]
degrees = [4, 6, 8, 10, 12, 14, 16, 18, 20]
rr = linspace(0.9*r0, cutoff(calc), 200)
for d in degrees
   B2 = [B1; gen_basis(2, D2, d)]
   @show length(B2)
   c = regression(B2, data, nforces = Inf, cstab = 0.0)
   @show norm(c, Inf)
   IP = NBodyIP(B2, c)
   rE, rF, mE, mF = fiterrors(IP, data)
   println("   E-rms, E-mae on testset = ", rE, ", ", mE)
   println("   F-rms, F-mae on testset = ", rF, ", ", mF)
   V2 = IP.orders[2]
   ev2 = vecnorm(V2.(rr) - calc.(rr), Inf)
   dev2 = vecnorm(evaluate_d.(V2, rr) - evaluate_d.(calc, rr), Inf)
   println("   V2 - uniform error = ", ev2, " | ", dev2)
   push!(err_eunif, ev2)
   push!(err_funif, dev2)
   push!(err_erms, rE)
   push!(err_frms, rF)
end

using DataFrames
df = DataFrame( :degrees => degrees,
                :unif_E => err_eunif,
                :unif_F => err_funif,
                :rms_E => err_erms,
                :rms_F => err_frms )
display(df)
