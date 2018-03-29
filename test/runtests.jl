using ManyBodyIPs
using JuLIP, Base.Test, StaticArrays

pex, pf = ManyBodyIPs.psym_monomial([0,1,2])

dict(:poly, 3)


pex, pf = psym_polys(3, dict(:poly, 3)...)






qex, qf = psym_polys_tot(3, dict(:poly, 3)...)
ex = qex[8]
g = ManyBodyIPs.gen_psym_grad(3, ex, :x)

R = @SVector rand(3)
@show R
g(R)
