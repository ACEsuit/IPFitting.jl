
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

r0 = rnn(:Si)
calc = StillingerWeber()
data = generate_data(:Si, 2, 0.2*r0, 30, calc)

TRANSFORM = "exp( - 2 * (r/$r0 - 1) )"

rcut2 = 2 * cutoff(calc)
CUTOFF2 = "(:cos, $(rcut2-1), $rcut2)"
D2 = BondLengthDesc(TRANSFORM, CUTOFF2)

rcut3 = 2 * cutoff(calc)
CUTOFF3 = "(:cos, $(rcut3-1), $rcut3)"
D3 = BondLengthDesc(TRANSFORM, CUTOFF3)

##
err_erms = Float64[]
err_frms = Float64[]
degrees = [4, 6, 8, 10, 12]
rr = linspace(0.9*r0, cutoff(calc), 200)
for deg3 in degrees
   # B = [B1; gen_basis(2, D2, deg2); gen_basis(3, D3, deg3)]
   B = [nbpolys(2, D2, 8); nbpolys(3, D3, deg3)]
   @show length(B)
   db = LsqDB("", B, data)
   IP, errs = Lsq.lsqfit( db,
                         E0 = 0.0,
                         configweights = Dict("rand" => 1.0),
                         dataweights   = Dict("E" => 100.0, "F" => 1.0) )
   push!(err_erms, Err.relrmse(errs, "rand", "E"))
   push!(err_frms, Err.relrmse(errs, "rand", "F"))
end

##
df = DataFrame( :degrees => degrees,
                :relrms_E => err_erms,
                :relrms_F => err_frms )
display(df)

(@test minimum(err_erms) < 0.0001) |> println
(@test minimum(err_frms) < 0.01) |> println
