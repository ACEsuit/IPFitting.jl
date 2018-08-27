using NBodyIPFitting
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

@testset "NBodyIPFitting" begin
   @testset "FIO" begin include("test_fio.jl") end
   @testset "Dat" begin include("test_dat.jl") end
   @testset "LsqDB" begin include("test_lsq_db.jl") end
   @testset "Lsq" begin include("test_lsq.jl") end
   @testset "Fitting" begin include("test_fit.jl") end
end
