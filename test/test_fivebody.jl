
using NBodyIPs
using JuLIP, Base.Test, StaticArrays, BenchmarkTools, Combinatorics


Ball = NBodyIPs.nbody_onlytuples(4,12)
B3 = NBodyIPs.nbody_tuples(4, 12)
