using NBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

@testset "NBodyIPs" begin
   # TODO: monomials testset
   @testset "Invariants" begin include("test_invariants.jl") end
   # TODO: fitting testset
   # TODO: Data testset
   # TODO: IO / serialization testset
end

@btime NBodyIPs.Polys.gen_tuples(3, 8)
@btime NBodyIPs.Polys.gen_tuples_new(3, 8)
