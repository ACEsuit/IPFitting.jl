using NBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

@testset "NBodyIPs" begin
   # TODO: monomials testset
   # TODO: fast_polys testset
   @testset "Misc" begin include("test_misc.jl") end
   @testset "Iterators" begin include("test_iterators.jl") end
   @testset "Invariants" begin include("test_invariants.jl") end
   # @testset "Polynomials" begin include("test_polynomials.jl") end
   # TODO: fitting testset
   # TODO: Data testset
   # TODO: IO / serialization testset
end
