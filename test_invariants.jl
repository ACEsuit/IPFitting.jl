using StaticArrays, BenchmarkTools

include("Invariant_generation/NBody_4_deg_10_julia_check.jl")

x = rand(6)
@btime invariants_Q6_check($x)


include("Invariant_generation/NBody_5_deg_6_julia_check.jl")
include("test_sec_5_body.jl")

x = rand(10)
# @btime invariants_Q10_check($x)

(Primary_inv, Sec_Inv, Irr_sec) = invariants_Q10_check(x)

display(Sec_inv[1:5])

Sec = invariants_Q10(x)
@btime invariants_Q10($x)

display(SVector(Irr_sec...)-SVector(Sec...))
