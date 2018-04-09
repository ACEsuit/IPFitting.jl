
using NBodyIPs
using JuLIP, Base.Test, StaticArrays, BenchmarkTools, Combinatorics

D = dict(:poly, 6)


B = NBodyIPs.polys_fourbody(D...)
