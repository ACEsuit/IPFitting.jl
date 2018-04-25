using NBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

@testset "NBodyIPs" begin
   # TODO: monomials testset
   @testset "Invariants" begin include("test_invariants.jl") end
   # TODO: fitting testset
end
