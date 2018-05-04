
using StaticArrays, BenchmarkTools, Combinatorics

# include("Invariant_generation/NBody_4_deg_10_julia_check.jl")
#
# x = rand(6)
# @btime invariants_Q6_check($x)



include("Invariant_generation/NBody_5_deg_6_julia_check.jl")
include("test_sec_5_body.jl")

x = rand(10)
@btime invariants_Q10_check($x)

(Primary_inv, Sec_Inv, Irr_sec) = invariants_Q10_check(x)

Sec = invariants_Q10(x)
@btime invariants_Q10($x)

SVector(Irr_sec[1:10]...)-SVector(Sec[1:10]...)
SVector(Irr_sec[11:20]...)-SVector(Sec[11:20]...)
SVector(Irr_sec[21:27]...)-SVector(Sec[21:27]...)

display(Irr_sec-collect(Sec))

display(maximum(abs.(Irr_sec-collect(Sec))))
