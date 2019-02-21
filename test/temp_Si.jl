
using JuLIP, NBodyIPs, NBodyIPFitting, StaticArrays, Plots

include(homedir() * "/Dropbox/PIBmat/W_Data/W.jl")
E0 = W.get_E0()

db = LsqDB(homedir() * "/scratch/nbodyips/W_4BBL_2s_med")

db

dataweights = Dict("E" => 50.0, "F" => 1.0, "V" => 1.0)
configweights = Dict("dia" => (1.0),
                     "surface_110" => 1.0,
                     "surface_001" => 1.0,
                     "surface_111" => 1.0,
                     )

r0 = rnn(:Si)
rcuts = unique(cutoff.(db.basis))
rcut4, rcut3, rcut2 = (sort(rcuts)...,)

reg = [ BAReg(2, 0.5*r0, 0.85*r0, creg = 1.0),
        BAReg(3, 0.7*r0, 0.85*r0, creg = 0.1),
        BAReg(4, 0.7*r0, 0.85*r0, creg = 0.1),
        ]


IPreg, inforeg = lsqfit( db; E0 = E0, dataweights = dataweights,
                           configweights = configweights,
                           regularisers = reg )

V20 = IPreg.components[2]
V21 = IPreg.components[3]
r = range(0.8*r0, stop=rcut2, length=200)
P = plot(r, V20.Vr.(r) );
plot!(P, r, V21.Vr.(r)  );
display(P)
