using NBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

@testset "NBodyIPs" begin
   @testset "Invariants" begin include("test_invariants.jl") end
end

# TODO  - write test set for fitting



NBodyIPs.Polynomials.ftrans_analyse((:poly, -1.234)...)
NBodyIPs.Polynomials.fcut_analyse((:cos, 1.234, 2.345)...)
