
using NBodyIPs
using JuLIP, Base.Test, StaticArrays, BenchmarkTools, Combinatorics


Ball = NBodyIPs.nbody_onlytuples(3,12)
B3 = NBodyIPs.nbody_tuples(3, 12)

trueNbody((1,0,0))
