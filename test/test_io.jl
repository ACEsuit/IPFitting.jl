# possible start for a testset

# data = NBodyIPs.Data.read("~/Dropbox/PIBmat/Ti_DFTB_Data/Ti_N54_T2000.xyz")
#
# D = Dictionary("@analytic r -> (2.0/r)^3",
#                "(:cos, 6.0, 9.0)")
#
# s = serialize(D)
# deserialize(Dictionary, s)
#
# println("generating basis functions")
#
# rcuts = [9.2, 6.2, 4.5]   # 9.2
# TRANSFORM = "@analytic r -> (2.9/r)^3"
# CUTOFF = ["(:cos, $(0.66*rcut), $(rcut))" for rcut in rcuts]
#
# D2 = Dictionary(TRANSFORM, CUTOFF[1])
# D3 = Dictionary(TRANSFORM, CUTOFF[2])
# D4 = Dictionary(TRANSFORM, CUTOFF[3])
#
# B1 = [NBody(1.0)]
# B2 = gen_basis(2, D2, 12)
# B3 = gen_basis(3, D3, 10)
# B4 = gen_basis(4, D4, 6)
# B = [B1; B2; B3; B4]
# c = rand(length(B))
# IP = NBodyIP(B, c)
#
# V = IP.orders[3]
# s = serialize(V)
# V3 = deserialize(NBody, s)
#
# NBodyIPs.IO.write(@__DIR__() * "/test.jld", IP)
# IP = NBodyIPs.IO.read(@__DIR__() * "/test.jld")
