
NBody = NBODY;
Deg = DEGREE;

using Combinatorics, StaticArrays, NBodyIPs

include("perm_svector_generator.jl")
include("misc.jl")

display("include ok")

generate_all_irr_sec(NBody,Deg)
