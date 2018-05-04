using Combinatorics, StaticArrays, NBodyIPs

include("invariants_generator.jl")
include("misc.jl")

# Parameters
NBody = 5;
Deg = 6;
# --------------
NBlengths = Int(NBody*(NBody-1)/2);

# -------------------------------------------
#
# Generate irreducible secondaries
#
# -------------------------------------------
filename = "data/NB_$NBody"*"_deg_$Deg"*"_irr_sec_text.jl";
filenamedata = "data/NB_$NBody""_deg_$Deg""_irr_invariants.jl";
preword = "# Irreducible secondaries for NBody=$NBody"*"and deg=$Deg \n"
pref = "IS"

generate_invariants(filenamedata,filename,NBlengths,Deg,preword,pref)


# -------------------------------------------
#
# Generate primary invariants
#
# -------------------------------------------
filename = "data/NB_$NBody"*"_deg_$Deg"*"_prim_text.jl";
filenamedata = "data/NB_$NBody""_deg_$Deg""_prim_invariants.jl";
preword = "# Primary invariants for NBody=$NBody"*"and deg=$Deg \n"
pref = "P"

generate_invariants(filenamedata,filename,NBlengths,Deg,preword,pref)
