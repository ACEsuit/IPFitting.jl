@info("Loading libraries...")
using JuLIP, ASE
using NBodyIPs
using NBodyIPs: bodyorder, degree
using NBodyIPFitting
using NBodyIPFitting.Lsq
using NBodyIPFitting.Data: observation, hasobservation, configname
using FileIO, Plots, DataFrames
using Distributions
using Test

@info("Loading data...")
include(homedir() * "/Dropbox/PIBmat/W_Data/W.jl")
data = W.loaddb()
@show length(data)

@info("Loading db file...")
dbname = homedir() * "/Research/01_En_cours/Post-doc/nbodyips/W_2BBL_r2_2"
@show dbname
@time db = LsqDB(dbname)
info(db)


num_config = Dict(
   "gamma_surface_vacancy"  =>  750,
   "vacancy"  =>  420,
   "slice_sample"  =>  2000,
   "gamma_surface"  =>  6183,
   "surface"  =>  180,
   "md_bulk"  =>  60,
   "dislocation_quadrupole"  =>  100 )

p_config = Dict(
   "gamma_surface_vacancy"  =>  1.0,
   "vacancy"  =>  1.0,
   "slice_sample"  =>  1.0,
   "gamma_surface"  =>  1.0,
   "surface"  =>  1.0,
   "md_bulk"  =>  1.0,
   "dislocation_quadrupole"  =>  1.0 )


function get_weights(db)
   cfgkeys = confignames(db)
   vals = ones(length(cfgkeys))
   configweights = Dict( k => (v, p_config[k]) for (k, v) in zip(cfgkeys, vals))
   dataweights = Dict( "E" => 1.0, "F" => 1.0, "V" => 0.0 )
   return configweights, dataweights
end

function get_Ibasis(db, degrees)
   Ibasis = Int[]
   for (i, (bo, deg)) in enumerate( zip( bodyorder.(db.basis), degree.(db.basis) ) )
      if degrees[bo] >= deg
         push!(Ibasis, i)
      end
   end
   return Ibasis
end


degrees = Dict(2 => 3, 3=>0, 4=>0,  5=>0)
configweights, dataweights = get_weights(db)
Ibasis = get_Ibasis(db, degrees)

IP, lsqinfo = lsqfit( db; E0 = W.get_E0(),
                      dataweights = dataweights,
                      configweights = configweights,
                      Ibasis = Ibasis,
                      verbose = false )

@time res_dict = results_dict(data, IP; confignames = Colon(), pathname = homedir() * "/Research/01_En_cours/Post-doc/nbodyips/res_dict_W.jld2")

errs1 = lsqerrors(db, [IP.components[2].c ; IP.components[3].c] , Ibasis; E0 = W.get_E0())
table(errs1)

errs2 = lsqerrors(res_dict, data; E0 = W.get_E0())
table(errs2)

for f in fieldnames(errs1)
   D1 = getfield(errs1, f)
   D2 = getfield(errs2, f)
   for cn in keys(configweights)
      for ot in ["E","F","V"]
         if !haskey(D1[cn],ot)
            continue
         end
         @test maximum(abs.(D1[cn][ot] - D2[cn][ot])) < 1e-6
      end
   end
end
