using JuLIP, ManyBodyIPs, StaticArrays

# here is a little code snippet that shows how to use
# the methods defined in `assembly.jl`

rcut = 4.1

# define the two-body terms
V2 = let lj = LennardJones() * C2Shift(rcut), rcut = rcut
   NBody(2, r->lj(r[1]), r -> (@D lj(r[1])), rcut, wrap=false)
end

# define the three-body terms
V3 = let rcut = rcut
   f(r) = prod( (r .- rcut).^2 )
   f_d(r::SVector{3,T}) where T =
      SVector{3,T}( 2 * (r[1]-rcut)   * (r[2]-rcut)^2 * (r[3]-rcut)^2,
                    2 * (r[1]-rcut)^2 * (r[2]-rcut)   * (r[3]-rcut)^2,
                    2 * (r[1]-rcut)^2 * (r[2]-rcut)^2 * (r[3]-rcut)    )
   NBody(3, f, f_d, 4.1, wrap = false)
end

# test on some atoms
at = bulk(:Cu, cubic=true, pbc = false) * 3
try
   @show V2(at)
   @show norm(@D V2(at))
   @show V3(at)
   @show norm(@D V3(at))
   @test true
catch
   @test false
end

# assemble into an interatomic potential
V = NBodyIP([V2, V3])
try
   @show energy(V, at)
   @show norm(forces(V, at))
   @test true
   @test JuLIP.Testing.fdtest(V, at)
catch
   @test false
end
