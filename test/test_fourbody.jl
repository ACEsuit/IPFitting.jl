
using NBodyIPs
using JuLIP, Base.Test, StaticArrays, BenchmarkTools, Combinatorics

D = dict(:poly, 6)


B = NBodyIPs.polys_fourbody(D...)

display(B[4])
ex, f, df = NBodyIPs.gen_fun(B[4], D...)
