
using StaticArrays, BenchmarkTools, Combinatorics

include("../data/NB_5_deg_6_non_efficient_invariants.jl")
include("invariants_5B.jl")

x = @SVector rand(10)

(Primary_inv, Sec_Inv, Irr_sec) = invariants_Q10_check(x)
(Primary,Sec) = invariants_Q10(x)

# ------------------
# Invariant check vs slow version
# ------------------

# Primary comparison
Primary - Primary_inv
maximum(abs.(SVector(Primary...) - Primary_inv))
#dont match yet...dont know why.

# Secondary comparison
SVector(Sec_Inv...) - Sec
maximum(abs.(SVector(Sec_Inv...) - Sec))

# ------------------
# Timings
# ------------------

@btime invariants_Q10_check($x)
@btime invariants_Q10($x)
