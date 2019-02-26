
module Errors

using JuLIP: Atoms, energy, forces, virial
using NBodyIPFitting: LsqDB, Dat, configtype, weighthook, eval_obs,
                      observation, hasobservation, observations,
                      vec_obs
using NBodyIPFitting.DB: matrows
using ASE: ASEAtoms   # TODO: WE SHOULD AVOID THIS!!!!
using FileIO, Printf
using NBodyIPs: fast

using Statistics: quantile
using LinearAlgebra: norm
using ProgressMeter

export add_fits!, rmse, mae, rmse_table, mae_table,
       lsqerrors

# ------------------------- USER INTERFACE FUNCTIONS -------------------------

rmse(data::Vector{Dat}; fitkey="fit") = _fiterrs(data; fitkey=fitkey, p=2)
# TODO: add back in!
# mae(data::Vector{Dat}; fitkey="fit") = _fiterrs(data; fitkey=fitkey, p=1)
# maxe(data::Vector{Dat}; fitkey="fit") = _fiterrs(data; fitkey=fitkey, p=Inf)

rmse_table(data::AbstractVector{Dat}) = rmse_table(rmse(data)...)
rmse_table(errs::Dict, errs_rel::Dict) = _err_table(errs, errs_rel, "RMSE")

rmse(D::Dict) = D["rmse"], D["relrmse"]


"""
`add_fits(IP, data::Vector{Dat}) -> nothing`

for each `d in data`, compute the observations with the IP (fitted to this
data) and store in `d.info["fit"][...]`

The key `"fit"` can be replaced with any string, through the kwarg
`fitkey`, e.g.,
```julia
data = ...
add_fits(PIP4Benv, data; fitkey = "PIP4Benv")
add_fits(GAP, data; fitkey = "GAP2010")
```
"""
function add_fits!(IP, data::Vector{Dat}; fitkey = "fit")
   @showprogress for d in data   # TODO: multi-threaded
      d.info[fitkey] = Dict{String, Vector{Float64}}()
      for okey in keys(d.D)   # d.D are the observations
         d.info[fitkey][okey] = vec_obs(okey, eval_obs(okey, IP, d.at))
      end
   end
   return nothing
end


function truncate_string(s, n)
   if length(s) <= n
      return s
   end
   n1 = n2 = (n-2) ÷ 2
   if n1 + n2 + 2 < n
      n1 += 1
   end
   return "$(s[1:n1])..$(s[end-n2+1:end])"
end

_relerr(relerrs, ct, ot) = haskey(relerrs[ct], ot) ? 100*relerrs[ct][ot] : NaN
_err(errs, ct, ot) = haskey(errs[ct], ot) ? errs[ct][ot] : NaN

function _err_table(errs, relerrs, title)
   lentitle = length(title)
   print("┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓\n")
   print("┃ "); printstyled(title; bold=true);
   print(repeat(' ', 70-3-lentitle)); print("┃\n")
   print("┣━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┯━━━━━━━━━━━━━━━━━┯━━━━━━━━━━━━━━━━━┫\n")
   print("┃ config type  ┃      E [eV]     │     F [eV/A]    │     V [eV/A]    ┃\n")
   print("┠──────────────╂────────┬────────┼────────┬────────┼────────┬────────┨\n")
   s_set = ""
   for ct in keys(errs)  # ct in configtypes
      s = @sprintf("┃ %12s ┃ %6.4f ┊ %5.2f%% │ %6.3f ┊ %5.2f%% │ %6.3f ┊ %5.2f%% ┃\n",
         truncate_string(ct, 12),
         _err(errs, ct, "E"), _relerr(errs, ct, "E"),
         _err(errs, ct, "F"), _relerr(errs, ct, "F"),
         _err(errs, ct, "V"), _relerr(errs, ct, "V") )
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
      w = weighthook(okey, dat)  # the default weight for this observation
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

@noinline function lsqerrors(db, c, Ibasis; cfgtypes = Colon(), E0 = nothing)

   cfgtypes isa Colon ?  cfgtypes = configtypes(db) :
                         cfgtypes = collect(cfgtypes)
   @assert E0 != nothing

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

   for (okey, dat, _) in observations(db)
      ct = configtype(dat)
      if !(ct in cfgtypes)
         continue
      end
      if !haskey(errs["rmse"][ct], okey)
         errs["rmse"][ct][okey] = 0.0
         errs["nrm2"][ct][okey] = 0.0
         errs["mae"][ct][okey] = 0.0
         errs["nrm1"][ct][okey] = 0.0
         # errs["allerr"][ct][okey] = Float64[]
         errs["maxe"][ct][okey] = 0.0
         errs["nrminf"][ct][okey] = 0.0
         lengths[ct][okey] = 0
      end
      if !haskey(errs["rmse"]["set"], okey)
         errs["rmse"]["set"][okey] = 0.0
         errs["nrm2"]["set"][okey] = 0.0
         errs["mae"]["set"][okey] = 0.0
         errs["nrm1"]["set"][okey] = 0.0
         # errs["allerr"]["set"][okey] = Float64[]
         errs["maxe"]["set"][okey] = 0.0
         errs["nrminf"]["set"][okey] = 0.0
         lengths["set"][okey] = 0
      end
      # assemble the observation vector
      y = Float64[]
      o = observation(dat, okey)
      # TODO: hack again!!!
      if okey == "E"
         o = copy(o)
         o[1] -= length(dat) * E0
      end
      append!(y, o)
      # get the lsq block
      irows = matrows(dat, okey)
      block = db.Ψ[irows, Ibasis]
      # the error on this particular data-point
      e = block * c - y
      # compute the various errors for this ct/okey combination
      w = weighthook(okey, dat)
      errs["rmse"][ct][okey] += w^2 * norm(e)^2
      errs["nrm2"][ct][okey] += w^2 * norm(y)^2
      errs["mae"][ct][okey]  += w * norm(e, 1)
      errs["nrm1"][ct][okey] += w * norm(y, 1)
      # append!(errs["allerr"][ct][okey], w * abs.(e) )
      errs["maxe"][ct][okey] = max(errs["maxe"][ct][okey], norm(e, Inf))
      errs["nrminf"][ct][okey] = max(errs["nrminf"][ct][okey], norm(y, Inf))
      lengths[ct][okey] += length(y)
      errs["rmse"]["set"][okey] += w^2 * norm(e)^2
      errs["nrm2"]["set"][okey] += w^2 * norm(y)^2
      errs["mae"]["set"][okey]  += w * norm(e, 1)
      errs["nrm1"]["set"][okey] += w * norm(y, 1)
      errs["maxe"]["set"][okey] = max(errs["maxe"]["set"][okey], norm(e, Inf))
      errs["nrminf"]["set"][okey] = max(errs["nrminf"]["set"][okey], norm(y, Inf))
      lengths["set"][okey] += length(y)
      # append!(errs["allerr"]["set"][okey], w * abs.(block * c - y) )
      # @show length(errs["allerr"]["set"][okey])
      # @show lengths["set"][okey]
   end

   for ct in [cfgtypes; "set"]
      for okey in keys(errs["rmse"][ct])
         len = lengths[ct][okey]
         errs["rmse"][ct][okey] = sqrt( errs["rmse"][ct][okey] / len )
         errs["nrm2"][ct][okey] = sqrt( errs["nrm2"][ct][okey] / len )
         errs["mae"][ct][okey] = errs["mae"][ct][okey] / len
         errs["nrm1"][ct][okey] = errs["nrm1"][ct][okey] / len
         errs["relrmse"][ct][okey] = errs["rmse"][ct][okey] / errs["nrm2"][ct][okey]
         # errs["allerr"][ct][okey] = quantile(errs["allerr"][ct][okey] ,range(0., stop=1., length=nb_points_cdf))
      end
   end

   # # ------- - scatter E -------
   # push!(scatterE[ct][1], E_data/len)
   # push!(scatterE[ct][2], E_fit/len)
   # append!(scatterF[ct][1], f_data)
   # append!(scatterF[ct][2], f_fit)

   return errs
end




end
