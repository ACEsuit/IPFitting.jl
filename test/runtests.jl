using ManyBodyIPs
using Base.Test
using StaticArrays

f(r) = prod(r)
f_d(r::SVector{3,T}) where T = SVector{3,T}(r[2]*r[3], r[1]*r[3], r[1]*r[2])
