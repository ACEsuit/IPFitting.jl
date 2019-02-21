
using JuLIP, NBodyIPs, NBodyIPFitting, StaticArrays, Plots

include(homedir() * "/Dropbox/PIBmat/W_Data/W.jl")
E0 = W.get_E0()

db = LsqDB(homedir() * "/scratch/nbodyips/W_4BBL_2s_med")

dataweights = Dict("E" => 50.0, "F" => 1.0, "V" => 1.0)
p = 0.03
configweights = Dict("surface"                => (1.0, p),
                     "gamma_surface"          => (1.0, p),
                     "slice_sample"           => (1.0, p),
                     "gamma_surface_vacancy"  => (1.0, p),
                     "md_bulk"                => (1.0, p),
                     "vacancy"                => (1.0, p),
                     "dislocation_quadrupole" => (1.0, p),
                     "vacancy"                => (1.0, p) )

r0 = rnn(:W)
rcuts = unique(cutoff.(db.basis))
rcut4, rcut3, rcut2 = (sort(rcuts)...,)

reg = [ BLReg(2, 0.5*r0, 0.85*r0, creg = 0.1),
        BLReg(3, 0.7*r0, 0.85*r0, creg = 0.1),
        BLReg(4, 0.7*r0, 0.85*r0, creg = 0.1),
        ]

IP, fitinfo = lsqfit( db; E0 = E0, dataweights = dataweights,
                          configweights = configweights,
                          regularisers = reg
        )

table_absolute(fitinfo["errors"])


##
#   p  ||  E [eV] | F[eV/A] | V[eV/A] ||  E [eV] | F[eV/A] | V[eV/A]
# -----||---------|---------|---------||---------|---------|---------
#  1.0 ||  0.0055 |  0.0607 |  0.0685 ||  0.0050 |  0.0442 |  0.0469
# -----||---------|---------|---------||---------|---------|---------
# 0.001||  0.0206 |  2.4942 |  0.9675 ||  0.0069 |  0.3046 |  0.3325
# 0.002||  0.0018 |  0.2982 |  0.1432 ||  0.0015 |  0.0650 |  0.0877
# 0.003||  0.0014 |  0.0751 |  0.1529 ||  0.0013 |  0.0420 |  0.0859
# 0.005||  0.0017 |  0.0687 |  0.1294 ||  0.0015 |  0.0421 |  0.0820
# 0.01 ||  0.0021 |  0.0611 |  0.0641 ||  0.0019 |  0.0432 |  0.0447
# 0.02 ||  0.0028 |  0.0585 |  0.0850 ||  0.0026 |  0.0416 |  0.0582
# 0.04 ||  0.0045 |  0.0586 |  0.0756 ||  0.0041 |  0.0419 |  0.0500
# 0.1  ||  0.0053 |  0.0597 |  0.0633 ||  0.0048 |  0.0435 |  0.0432
