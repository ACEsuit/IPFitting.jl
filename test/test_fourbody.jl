
using NBodyIPs
using JuLIP, Base.Test, StaticArrays, BenchmarkTools, Combinatorics

D = dict(:poly, 8)
B = NBodyIPs.polys_fourbody(D...)
Bsym = NBodyIPs.psym_polys_nbody(4, D...)

# D = dict(:poly, 10)
# B = NBodyIPs.polys_fourbody(D...)
# Bsym = NBodyIPs.psym_polys_nbody(4, D...)


display(B[4])
ex, f, df = NBodyIPs.gen_fun(B[4], D...)



D = dict(:poly, 4)


B = NBodyIPs.polys_fourbody(6)
