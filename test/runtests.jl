using NBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

@testset "NBodyIPs" begin
   @testset "Invariants" begin include("test_invariants.jl") end
end

# TODO  - write test set for fitting
