using ManyBodyIPs
using JuLIP, Base.Test, StaticArrays

pex, pf = ManyBodyIPs.psym_monomial([0,1,2])

dict(:poly, 3)


p = psym_polys(3, dict(:poly, 3)...)
p






q = psym_polys_tot(3, dict(:poly, 3)...)
q
