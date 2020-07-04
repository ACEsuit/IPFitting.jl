
module Errors

using JuLIP: Atoms, energy, forces, virial
using IPFitting: LsqDB, Dat, configtype, err_weighthook, eval_obs,
                      observation, hasobservation, observations,
                      vec_obs, tfor_observations, truncate_string
using IPFitting.DB: matrows
using FileIO, Printf, Base.Threads

using Statistics: quantile
using LinearAlgebra: norm
using ProgressMeter
using Random: shuffle!

export add_fits!, add_fits_serial!, rmse, mae, rmse_table, rmse_table_wg, mae_table,
       lsqerrors, splittraintest


# ------------------------- USER INTERFACE FUNCTIONS -------------------------

rmse(data::Vector{Dat}; fitkey="fit") = _fiterrs(data; fitkey=fitkey, p=2)
mae(data::Vector{Dat}; fitkey="fit") = _fiterrs(data; fitkey=fitkey, p=1)

rmse(D::Dict) = D["rmse"], D["relrmse"]
mae(D::Dict) = D["mae"], D["relmae"]

rmse_table(data::AbstractVector{Dat}) = rmse_table(rmse(data)...)
mae_table( data::AbstractVector{Dat}) =  mae_table( mae(data)...)

rmse_table(D::Dict) = rmse_table(rmse(D)...)

rmse_table_wg(D::Dict, weights) = rmse_table_wg(rmse(D)..., weights)

mae_table(D::Dict) =  mae_table( mae(D)...)

rmse_table(errs::Dict, errs_rel::Dict; configtypes=:) =
   _err_table(errs, errs_rel, "RMSE", configtypes)

rmse_table_wg(errs::Dict, errs_rel::Dict, weights; configtypes=:) =
      _err_table_wg(errs, errs_rel, "RMSE WEIGHTED TABLE", weights, configtypes)

mae_table(errs::Dict, errs_rel::Dict) = _err_table(errs, errs_rel, "MAE")




"""
`add_fits(IP, data::Vector{Dat}) -> nothing`

for each `d in data`, compute the observations with the IP (fitted to this
data) and store in `d.info["fit"][...]`

The key `"fit"` can be replaced with any string, through the kwarg
`fitkey`, e.g.,
```julia
data = ...
add_fits!(PIP4Benv, data; fitkey = "PIP4Benv")
add_fits!(GAP, data; fitkey = "GAP2010")
```
"""
function add_fits!(IP, configs::Vector{Dat}; fitkey = "fit")
   # create the nec essary dictionaries
   for cfg in configs
      cfg.info[fitkey] = Dict{String, Vector{Float64}}()
   end
   tfor_observations( configs,
         (n, okey, cfg, lck) -> begin
            obs = vec_obs(okey, eval_obs(okey, IP, cfg))
            lock(lck)
            cfg.info[fitkey][okey] = obs
            unlock(lck)
         end;
         msg="Add Fit info to configs")

   return nothing
end


function add_fits_serial!(IP, configs::Vector{Dat}; fitkey = "fit")
   # create the nec essary dictionaries
   for cfg in configs
      cfg.info[fitkey] = Dict{String, Vector{Float64}}()
   end
   for (okey, cfg, _) in observations(configs)
      obs = vec_obs(okey, eval_obs(okey, IP, cfg))
      cfg.info[fitkey][okey] = obs
   end
   return nothing
end

function _relerr(relerrs, ct, ot)
   if haskey(relerrs[ct], ot)
      if relerrs[ct][ot] isa Number
         return 100*relerrs[ct][ot]
      end
   end
   return NaN
end

function _err(errs, ct, ot)
   if  haskey(errs[ct], ot)
      if errs[ct][ot] isa Number
         return errs[ct][ot]
      end
   end
   return NaN
end

function _err_table_wg(errs, relerrs, title, weights, configtypes=:)
   lentitle = length(title)
   print("┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓\n")
   print("┃ "); printstyled(title; bold=true);
   print(repeat(' ', 70-3-lentitle)); print("┃\n")
   print("┣━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┯━━━━━━━━━━━━━━━━━┯━━━━━━━━━━━━━━━━━┫\n")
   print("┃ config type  ┃     E [meV]     │     F [eV/A]    │     V [meV]     ┃\n")
   print("┠──────────────╂────────┬────────┼────────┬────────┼────────┬────────┨\n")
   s_set = ""
   for ct in keys(errs)  # ct in configtypes
      if (configtypes != :)
         if !(ct in configtypes)
            continue
         end
      end

      if haskey(weights, ct)
         w = weights[ct]
      else
         w = weights["default"]
      end

      s = @sprintf("┃ %12s ┃ %6.2f ┊ %5.3f%% │ %6.3f ┊ %5.2f%% │ %6.1f ┊ %5.2f%% ┃\n",
         truncate_string(ct, 12),
         _err(errs, ct, "E")*1000 * w["E"], _relerr(relerrs, ct, "E") * w["E"],
         _err(errs, ct, "F") * w["F"], _relerr(relerrs, ct, "F") * w["F"],
         _err(errs, ct, "V")*1000 * w["V"], _relerr(relerrs, ct, "V") * w["V"])
      if ct == "set"
         s_set = s
      else
         print(s)
      end
   end
   if s_set != ""
      print("┠──────────────╂────────┼────────┼────────┼────────┼────────┼────────┨\n")
      print(s_set)
   end
      print("┗━━━━━━━━━━━━━━┻━━━━━━━━┷━━━━━━━━┷━━━━━━━━┷━━━━━━━━┷━━━━━━━━┷━━━━━━━━┛\n")
end

function _err_table(errs, relerrs, title, configtypes=:)
   lentitle = length(title)
   print("┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓\n")
   print("┃ "); printstyled(title; bold=true);
   print(repeat(' ', 70-3-lentitle)); print("┃\n")
   print("┣━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┯━━━━━━━━━━━━━━━━━┯━━━━━━━━━━━━━━━━━┫\n")
   print("┃ config type  ┃     E [meV]     │     F [eV/A]    │     V [meV]     ┃\n")
   print("┠──────────────╂────────┬────────┼────────┬────────┼────────┬────────┨\n")
   s_set = ""
   for ct in keys(errs)  # ct in configtypes
      if (configtypes != :)
         if !(ct in configtypes)
            continue
         end
      end

      s = @sprintf("┃ %12s ┃ %6.2f ┊ %5.3f%% │ %6.3f ┊ %5.2f%% │ %6.1f ┊ %5.2f%% ┃\n",
         truncate_string(ct, 12),
         _err(errs, ct, "E")*1000, _relerr(relerrs, ct, "E"),
         _err(errs, ct, "F"), _relerr(relerrs, ct, "F"),
         _err(errs, ct, "V")*1000, _relerr(relerrs, ct, "V"))
      if ct == "set"
         s_set = s
      else
         print(s)
      end
   end
   if s_set != ""
      print("┠──────────────╂────────┼────────┼────────┼────────┼────────┼────────┨\n")
      print(s_set)
   end
      print("┗━━━━━━━━━━━━━━┻━━━━━━━━┷━━━━━━━━┷━━━━━━━━┷━━━━━━━━┷━━━━━━━━┷━━━━━━━━┛\n")
end


# ---------------------------------------------------------------------
# ------------------------- PRIVATE FUNCTIONS -------------------------
# ---------------------------------------------------------------------

"""
initialize an error dictionary
"""
function errdict(configtypes)
   D = Dict{String, Dict{String, Float64}}()
   for ct in configtypes
      D[ct] = Dict{String, Float64}()
   end
   D["set"] = Dict{String, Float64}()
   return D
end



function _fiterrs(configs::Vector{Dat}; fitkey="fit", p=2)
   configtypes = unique(configtype.(configs))
   err = errdict(configtypes)
   relerr = errdict(configtypes)
   nrm = errdict(configtypes)
   len = errdict(configtypes)

   # assemble the errors
   for (okey, dat, _) in observations(configs)
      ct = configtype(dat)
      if !haskey(err[ct], okey)
         len[ct][okey] = 0.0
         nrm[ct][okey] = 0.0
         err[ct][okey] = 0.0
      end
      if !haskey(err["set"], okey)
         len["set"][okey] = 0.0
         nrm["set"][okey] = 0.0
         err["set"][okey] = 0.0
      end
      w = err_weighthook(okey, dat)  # the default weight for this observation
      exval = dat.D[okey] * w
      fitval = dat.info[fitkey][okey] * w
      len[ct][okey] += length(exval)
      len["set"][okey] += length(exval)
      nrm[ct][okey] += norm(exval, p)^p
      nrm["set"][okey] += norm(exval, p)^p
      err[ct][okey] += norm(exval-fitval, p)^p
      err["set"][okey] += norm(exval-fitval, p)^p
   end

   # normalise correctly
   for ct in keys(err)
      for okey in keys(err[ct])
         err[ct][okey] = (err[ct][okey]/len[ct][okey])^(1/p)
         nrm[ct][okey] = (nrm[ct][okey]/len[ct][okey])^(1/p)
         relerr[ct][okey] = err[ct][okey] / nrm[ct][okey]
      end
   end

   return err, relerr
end


# TODO: the next fucntion is to be retired!!

@noinline function lsqerrors(db, c, Ibasis; cfgtypes = Colon(), Vref = nothing,
                                            Icfg=Colon())

   cfgtypes isa Colon ?  cfgtypes = configtypes(db) :
                         cfgtypes = collect(cfgtypes)

   # create the dict for the fit errors
   errs = Dict("rmse" => errdict([cfgtypes; "set"]),
               "relrmse" => errdict([cfgtypes; "set"]),
               "mae"  => errdict([cfgtypes; "set"]),
               "maxe" => errdict([cfgtypes; "set"]),
               "nrm2" => errdict([cfgtypes; "set"]),
               "nrminf" => errdict([cfgtypes; "set"]),
               "nrm1" => errdict([cfgtypes; "set"]) )
   lengths = Dict{String, Dict{String, Int}}()
   for ct in [cfgtypes; "set"]
      lengths[ct] = Dict{String, Int}()
   end

   for (okey, dat, icfg) in observations(db, Icfg)

      ct = configtype(dat)
      if !(ct in cfgtypes)
         continue
      end
      if !haskey(errs["rmse"][ct], okey)
         errs["rmse"][ct][okey] = 0.0
         errs["nrm2"][ct][okey] = 0.0
         errs["mae"][ct][okey] = 0.0
         errs["nrm1"][ct][okey] = 0.0
         errs["maxe"][ct][okey] = 0.0
         errs["nrminf"][ct][okey] = 0.0
         lengths[ct][okey] = 0
      end
      if !haskey(errs["rmse"]["set"], okey)
         errs["rmse"]["set"][okey] = 0.0
         errs["nrm2"]["set"][okey] = 0.0
         errs["mae"]["set"][okey] = 0.0
         errs["nrm1"]["set"][okey] = 0.0
         errs["maxe"]["set"][okey] = 0.0
         errs["nrminf"]["set"][okey] = 0.0
         lengths["set"][okey] = 0
      end
      # assemble the observation vector
      y = Float64[]
      o = observation(dat, okey)
      # correct the observation if we are fitting from a reference potential
      if Vref != nothing
         o -= vec_obs(okey, eval_obs(okey, Vref, dat))
      end
      append!(y, o)
      # get the lsq block
      irows = matrows(dat, okey)
      block = db.Ψ[irows, Ibasis]
      # the error on this particular data-point
      e = block * c - y
      # compute the various errors for this ct/okey combination
      w = err_weighthook(okey, dat)
      errs["rmse"][ct][okey] += norm(w .* e)^2
      errs["nrm2"][ct][okey] += norm(w .* y)^2
      errs["mae"][ct][okey]  += norm(w .* e, 1)
      errs["nrm1"][ct][okey] += norm(w .* y, 1)
      errs["maxe"][ct][okey] = max(errs["maxe"][ct][okey], norm(e, Inf))
      errs["nrminf"][ct][okey] = max(errs["nrminf"][ct][okey], norm(y, Inf))
      lengths[ct][okey] += length(y)
      errs["rmse"]["set"][okey] += norm(w .* e)^2
      errs["nrm2"]["set"][okey] += norm(w .* y)^2
      errs["mae"]["set"][okey]  += norm(w .* e, 1)
      errs["nrm1"]["set"][okey] += norm(w .* y, 1)
      errs["maxe"]["set"][okey] = max(errs["maxe"]["set"][okey], norm(e, Inf))
      errs["nrminf"]["set"][okey] = max(errs["nrminf"]["set"][okey], norm(y, Inf))
      lengths["set"][okey] += length(y)
   end

   for ct in [cfgtypes; "set"]
      for okey in keys(errs["rmse"][ct])
         len = lengths[ct][okey]
         errs["rmse"][ct][okey] = sqrt( errs["rmse"][ct][okey] / len )
         errs["nrm2"][ct][okey] = sqrt( errs["nrm2"][ct][okey] / len )
         errs["mae"][ct][okey] = errs["mae"][ct][okey] / len
         errs["nrm1"][ct][okey] = errs["nrm1"][ct][okey] / len
         errs["relrmse"][ct][okey] = errs["rmse"][ct][okey] / errs["nrm2"][ct][okey]
      end
   end

   return errs
end


splittraintest(db::LsqDB, ptrain=0.8) = splittraintest(db.configs, ptrain)

function splittraintest(cfgs::AbstractVector{Dat}, ptrain=0.8)
   cfgtypes = unique(configtype.(cfgs))
   Itrain = Int[]
   Itest = Int[]
   for ct in cfgtypes
      # find all indices of configurations with the current configtype
      Ict = findall(configtype.(cfgs) .== ct)
      # randomly permute the indices
      shuffle!(Ict)
      # first 80% go into train, remaining 20% into test
      ntrain = ceil(Int, length(Ict) * ptrain)
      append!(Itrain, Ict[1:ntrain])
      append!(Itest, Ict[ntrain+1:end])
   end
   return Itrain, Itest
end


end
