
using NBodyIPs, JuLIP, Test, NBodyIPFitting, DataFrames
using JuLIP.Potentials: evaluate_d
using NBodyIPFitting: Dat, LsqDB
using NBodyIPs: BondLengthDesc, BondAngleDesc
Lsq = NBodyIPFitting.Lsq
Err = NBodyIPFitting.Errors
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

r0 = rnn(:Si)
calc = StillingerWeber()
data = generate_data(:Si, 2, 0.33*r0, 100, calc)

TRANSFORM = "exp( - 2 * (r/$r0 - 1.5) )"

rcut2 = cutoff(calc)*2.0
CUTOFF2 = PolyCutSym(2, rcut2)
D2 = BondLengthDesc(TRANSFORM, CUTOFF2)

rcut3 = cutoff(calc)*2.0
# CUTOFF3 = CosCut(rcut3-1, rcut3)
CUTOFF3 = PolyCut(2, rcut3)
# CUTOFF3 = PolyCut2s(2, 0.0, rnn(:Si), rcut3)

rcut3 = cutoff(calc)*2.0
D3 = BondLengthDesc(TRANSFORM, CUTOFF3)

##
err_erms = Float64[]
err_frms = Float64[]
degrees = [4, 6, 8, 10] # [4, 6, 8, 10, 12]
rr = range(0.9*r0, stop=cutoff(calc), length=200)
for deg3 in degrees
   # B = [B1; gen_basis(2, D2, deg2); gen_basis(3, D3, deg3)]
   B = [nbpolys(2, D2, 10); nbpolys(3, D3, deg3)]
   # B = nbpolys(3, D3, deg3)
   @show length(B)
   db = LsqDB("", B, data)
   IP, fitinfo = Lsq.lsqfit( db,
                         E0 = 0.0,
                         configweights = Dict("rand" => 1.0),
                         obsweights   = Dict("E" => 100.0, "F" => 1.0),
                         combineIP = NBodyIP )
   push!(err_erms, fitinfo["errors"]["relrmse"]["set"]["E"])
   push!(err_frms, fitinfo["errors"]["relrmse"]["set"]["F"])
end



##
df = DataFrame( :degrees => degrees,
                :relrms_E => err_erms,
                :relrms_F => err_frms )
display(df)

(@test minimum(err_erms) < 0.001) |> println
(@test minimum(err_frms) < 0.05) |> println


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
