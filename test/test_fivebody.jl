
using NBodyIPs
using JuLIP, Base.Test, StaticArrays, BenchmarkTools, Combinatorics


Ball = NBodyIPs.nbody_onlytuples(5,5)
B3 = NBodyIPs.nbody_tuples(4, 10)


for len in 4:10
   Ball = NBodyIPs.nbody_alltuples(4, len)
   B2 = NBodyIPs.nbody_tuples(2, len)
   B3 = NBodyIPs.nbody_tuples(3, len)
   B4 = NBodyIPs.nbody_tuples(4, len)
   @test length(Ball) == length([B2; B3; B4])
end
