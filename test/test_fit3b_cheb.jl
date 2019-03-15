
using NBodyIPs, JuLIP, Test, IPFitting, DataFrames
using JuLIP.Potentials: evaluate_d
using IPFitting: Dat, LsqDB
using NBodyIPs: BondLengthDesc, BondAngleDesc
using NBodyIPs.Polys: Chebyshev
const Lsq = IPFitting.Lsq
const Err = IPFitting.Errors


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

TRANSFORM = "exp( - 2*(r/$r0 - 1))"

# exp( - 2*(0.9*r0/r0-1))

rcut2 = cutoff(calc)*1.6
CUTOFF2 = "(:cos, $(rcut2-1), $rcut2)"
D2 = BondLengthDesc(TRANSFORM, CUTOFF2)

rcut3 = cutoff(calc)*1.6
CUTOFF3 = "(:cos, $(rcut3-1), $rcut3)"
D3 = BondLengthDesc(TRANSFORM, CUTOFF3)

##
err_erms = Float64[]
err_frms = Float64[]
degrees = [4, 6, 8, 10] # , 10, 12]
rr = range(0.9*r0, stop=cutoff(calc), length=200)
for deg3 in degrees
   # B = [B1; gen_basis(2, D2, deg2); gen_basis(3, D3, deg3)]
   # B = [nbpolys(2, D2, 8, Chebyshev()); nbpolys(3, D3, deg3, Chebyshev())]
   B = nbpolys(3, D3, deg3, Chebyshev())
   @show length(B)
   db = LsqDB("", B, data)
   IP, errs = Lsq.lsqfit( db,
                         E0 = 0.0,
                         configweights = Dict("rand" => 1.0),
                         dataweights   = Dict("E" => 100.0, "F" => 1.0) )
   # push!(err_erms, Err.relrmse(errs, "rand", "E"))
   # push!(err_frms, Err.relrmse(errs, "rand", "F"))
end


##
df = DataFrame( :degrees => degrees,
                :relrms_E => err_erms,
                :relrms_F => err_frms )
display(df)

(@test minimum(err_erms) < 0.0001) |> println
(@test minimum(err_frms) < 0.01) |> println


# BOND LENGTH
# │ Row │ degrees │ rms_E       │ rms_F      │
# ├─────┼─────────┼─────────────┼────────────┤
# │ 1   │ 4       │ 0.00134475  │ 0.118625   │
# │ 2   │ 6       │ 0.000244722 │ 0.0362398  │
# │ 3   │ 8       │ 0.000115287 │ 0.0179257  │
# │ 4   │ 10      │ 5.28124e-5  │ 0.00888845 │

# BOND ANGLE
# │ Row │ degrees │ relrms_E    │ relrms_F   │
# ├─────┼─────────┼─────────────┼────────────┤
# │ 1   │ 4       │ 0.000523011 │ 0.070237   │
# │ 2   │ 6       │ 0.000445928 │ 0.0477666  │
# │ 3   │ 8       │ 0.000100544 │ 0.016773   │
# │ 4   │ 10      │ 3.7269e-5   │ 0.00568809 │
# │ 5   │ 12      │ 2.00043e-5  │ 0.00419453 │
