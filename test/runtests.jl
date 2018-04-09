using NBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools


@testset "NBodyIPs" begin
   @testset "Polynomials" begin include("test_polynomials.jl") end
   @testset "Assembly" begin include("test_assemble.jl") end
end


# D = dict(:poly, 10)
# ex2, _, _ = NBodyIPs.psym_polys_nbody(2, dict(:poly, 14)...; simplify=false)
# ex3, _, _ = NBodyIPs.psym_polys_nbody(3, D...; simplify=false)
# ex4, _, _ = NBodyIPs.psym_polys_nbody(4, D...; simplify=false)
# ex5, _, _ = NBodyIPs.psym_polys_nbody(5, D...; simplify=false)

# ex, f, df = NBodyIPs.psym_monomial([2,0,2], dict(:inv2, 5, 3.0))
#
# b = NBody(4, f, df, 3.0)
# R = 1.0 + @SVector rand(3)
# b.f(R)
# f(R)
#
# psym_polys(3, dict(:inv2, 5, 3.0)...; simplify = true)
#
# # ord = 5
# # deg = (ord * (ord-1)) ÷ 2
#
# using Combinatorics
# dim = 3
# deg = 5
# for i in 1:deg
# 	for m = 1:dim, alpha in collect(partitions(i, m))
#       if length(alpha) < dim
#          append!(alpha, zeros(Int64, dim - length(alpha)))
#       end
#       @show alpha
#    end
# end
#
#
# using Combinatorics
# dim = 3
# deg = 5
# for i in 1:deg
# 	for m = 1:dim, alpha in partitions(i, m)
#       append!(alpha, zeros(Int64, dim - length(alpha)))
#       @show alpha
#    end
# end
