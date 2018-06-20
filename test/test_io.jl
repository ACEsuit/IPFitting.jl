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


# testing the custom JLD serialization

using JuLIP, NBodyIPs, DataFrames, Plots, JLD
using NBodyIPs.Data: Dat

include(homedir() * "/Dropbox/PIBmat/Ti_DFTB_Data/Ti.jl")
data = Ti.read_Ti(; exclude = ["wire", "surface", "omega", "hcp"])

# reference energy
B1 = [ NBody(Ti.get_E0()) ]

# generate long-range 2B basis
r0 = rnn(:Ti)
TLONG = "@analytic r -> exp( - 2.5 * (r/$r0 - 1))"
TSHORT = "@analytic r -> ($r0/r)^8"
CUT2 = "(:cos, 5.5, 7.5)"
B2 = [ gen_basis(2, Dictionary(TLONG, CUT2), 9);
       gen_basis(2, Dictionary(TSHORT, CUT2), 6) ]

# 3B BASIS
TRANS3 = "@analytic r -> exp( - 3.0 * (r/$r0 - 1))"
CUT3 = "(:cos, 5.0, 6.5)"
B3 = gen_basis(3, Dictionary(TRANS3, CUT3), 8)

B = [B1; B2; B3]
@show length(B)

lsq = kron(data[1:10], B)
JLD.save("temp.jld", "d1", lsq)
d1 = JLD.load("temp.jld", "d1")


JLD.@load "temp.jld" lsq






using FileIO
save("ip_temp.jld2", "IP", IP4)
IP4b = load("ip_temp.jld2", "IP")
IP4bf = fast(IP4b)
using StaticArrays
for n in [1,3,6]
    for _ = 1:10
        r = SVector((3.0 + rand(n))...)
        @assert IP4(r) ≈ IP4b(r)
        @assert IP4(r) ≈ IP4bf(r)
    end
end
