using NBodyIPFitting
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

@testset "NBodyIPFitting" begin
   @testset "LsqDB" begin include("test_lsq_db.jl") end
   # @testset "Fitting" begin include("test_fit.jl") end
   # TODO: Data testset
   # TODO: IO / serialization testset
end
