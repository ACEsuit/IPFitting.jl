

# here is a little code snippet that shows how to use
# the methods defined in `assembly.jl`

rcut = 4.1

# define the two-body terms
V2 = let lj = LennardJones() * C2Shift(rcut), rcut = rcut
   NBody(2, r->lj(r[1]), r -> (@D lj(r[1])), rcut)
end

# define the three-body terms
V3 = let rcut = rcut
   f(r) = prod( (r .- rcut).^2 )
   f_d(r::SVector{3,T}) where T =
      SVector{3,T}( 2 * (r[1]-rcut)   * (r[2]-rcut)^2 * (r[3]-rcut)^2,
                    2 * (r[1]-rcut)^2 * (r[2]-rcut)   * (r[3]-rcut)^2,
                    2 * (r[1]-rcut)^2 * (r[2]-rcut)^2 * (r[3]-rcut)    )
   NBody(3, f, f_d, 4.1)
end

# test on some atoms
at = bulk(:Cu, cubic=true, pbc = false) * 3
@show V2.f([1.0])
@show V2(at)
@show norm(@D V2(at))
@show V3.f([1.0, 1.6, 0.8])
@show V3(at)
@show norm(@D V3(at))

# assemble into an interatomic potential
V = NBodyIP([V2, V3])
@show energy(V, at)
@show norm(forces(V, at))
