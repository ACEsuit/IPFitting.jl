using ManyBodyIPs, JuLIP
using Base.Test
using StaticArrays

f(r) = prod(r)
f_d(r::SVector{3,T}) where T = SVector{3,T}(r[2]*r[3], r[1]*r[3], r[1]*r[2])
V3 = NBody(3, f, f_d, 4.1)
at = bulk(:Cu, cubic=true, pbc = false) * 3
V3(at)
@D V3(at)
