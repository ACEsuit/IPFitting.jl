using ManyBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

ex, f, df = ManyBodyIPs.psym_monomial([2,0,2], dict(:inv2, 5, 3.0))

b = NBody(3, f, df, 3.0)
R = 1.0 + @SVector rand(3)
b.f(R)
f(R)

psym_polys(3, dict(:poly, 5, 3.0)...; simplify = true)


using Combinatorics
dim = 3
deg = 4
for i in 1:deg
	for m = 1:dim, alpha in partitions(i, m)
      @show alpha
   end
end
