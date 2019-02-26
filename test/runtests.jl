using NBodyIPFitting
using JuLIP, Test, StaticArrays
using BenchmarkTools

@testset "NBodyIPFitting" begin
   @testset "Dat" begin include("test_dat.jl") end
   @testset "LsqDB" begin include("test_lsq_db.jl") end
   @testset "Lsq" begin include("test_lsq.jl") end
   @testset "Fitting" begin include("test_fit.jl") end
   @testset "Errors" begin include("test_errors.jl") end 
end
