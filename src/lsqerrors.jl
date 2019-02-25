
module Errors

using JuLIP: Atoms, energy, forces, virial
using NBodyIPFitting: LsqDB, Dat, configtype, weighthook, evaluate_lsq
using NBodyIPFitting.Data: observation, hasobservation, configname
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
         d.info[fitkey][okey] = vec(okey, evaluate_lsq(okey, IP, d.at))
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
   print("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
   print("                       $title    \n")
   print("───────────────┰─────────────────┬─────────────────┬─────────────────\n")
   print("  config type  ┃      E [eV]     │     F [eV/A]    │     V [eV/A]   \n")
   print("───────────────╂────────┬────────┼────────┬────────┼────────┬────────\n")
   s_set = ""
   for ct in keys(errs)  # ct in configtypes
      s = @sprintf(" %13s ┃ %6.4f ┊ %5.2f%% │ %6.3f ┊ %5.2f%% │ %6.3f ┊ %5.2f%%\n",
         truncate_string(ct, 13),
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
      print("───────────────╂────────┼────────┼────────┼────────┼────────┼────────\n")
      print(s_set)
   end
      print("━━━━━━━━━━━━━━━┷━━━━━━━━┷━━━━━━━━┷━━━━━━━━┷━━━━━━━━┷━━━━━━━━┷━━━━━━━━\n")
end


# ---------------------------------------------------------------------
# ------------------------- PRIVATE FUNCTIONS -------------------------
# ---------------------------------------------------------------------

"""
initialize an error dictionary
"""
function errdict(confignames)
   D = Dict{String, Dict{String, Float64}}()
   for ct in confignames
      D[ct] = Dict{String, Float64}()
   end
   D["set"] = Dict{String, Float64}()
   return D
end



function _fiterrs(data::Vector{Dat}; fitkey="fit", p=2)
   confignames = unique(configname.(data))
   err = errdict(confignames)
   relerr = errdict(confignames)
   nrm = errdict(confignames)
   len = errdict(confignames)

   # assemble the errors
   for d in data
      cn = configname(d)
      for ot in keys(d.D) # the observation types (E, F, ...)
         if !haskey(err[cn], ot)
            len[cn][ot] = 0.0
            nrm[cn][ot] = 0.0
            err[cn][ot] = 0.0
         end
         if !haskey(err["set"], ot)
            len["set"][ot] = 0.0
            nrm["set"][ot] = 0.0
            err["set"][ot] = 0.0
         end
         w = weighthook(ot, d)  # the default weight for this observation
         exval = d.D[ot] * w
         fitval = d.info[fitkey][ot] * w
         len[cn][ot] += length(exval)
         len["set"][ot] += length(exval)
         nrm[cn][ot] += norm(exval, p)^p
         nrm["set"][ot] += norm(exval, p)^p
         err[cn][ot] += norm(exval-fitval, p)^p
         err["set"][ot] += norm(exval-fitval, p)^p
      end
   end

   # normalise correctly
   for cn in keys(err)
      for ot in keys(err[cn])
         err[cn][ot] = (err[cn][ot]/len[cn][ot])^(1/p)
         nrm[cn][ot] = (nrm[cn][ot]/len[cn][ot])^(1/p)
         relerr[cn][ot] = err[cn][ot] / nrm[cn][ot]
      end
   end

   return err, relerr
end


# TODO: the next fucntion is to be retired!!

@noinline function lsqerrors(db, c, Ibasis; confignames = Colon(), E0 = nothing)
   if confignames isa Colon
      confignames = unique(configname.(collect(keys(db.data_groups))))
   else
      confignames = collect(confignames)
   end
   @assert E0 != nothing

   # create the dict for the fit errors
   errs = Dict("rmse" => errdict([confignames; "set"]),
               "relrmse" => errdict([confignames; "set"]),
               "mae"  => errdict([confignames; "set"]),
               "maxe" => errdict([confignames; "set"]),
               "nrm2" => errdict([confignames; "set"]),
               "nrminf" => errdict([confignames; "set"]),
               "nrm1" => errdict([confignames; "set"]) )
   lengths = Dict{String, Dict{String, Int}}()
   for cn in [confignames; "set"]
      lengths[cn] = Dict{String, Int}()
   end

   for ct in keys(db.data_groups)
      cn = configname(ct)
      if !(cn in confignames)
         continue
      end
      # loop through observation types in the current config type
      for ot in keys(db.kron_groups[ct])
         if !haskey(errs["rmse"][cn], ot)
            errs["rmse"][cn][ot] = 0.0
            errs["nrm2"][cn][ot] = 0.0
            errs["mae"][cn][ot] = 0.0
            errs["nrm1"][cn][ot] = 0.0
            # errs["allerr"][cn][ot] = Float64[]
            errs["maxe"][cn][ot] = 0.0
            errs["nrminf"][cn][ot] = 0.0
            lengths[cn][ot] = 0
         end
         if !haskey(errs["rmse"]["set"], ot)
            errs["rmse"]["set"][ot] = 0.0
            errs["nrm2"]["set"][ot] = 0.0
            errs["mae"]["set"][ot] = 0.0
            errs["nrm1"]["set"][ot] = 0.0
            # errs["allerr"]["set"][ot] = Float64[]
            errs["maxe"]["set"][ot] = 0.0
            errs["nrminf"]["set"][ot] = 0.0
            lengths["set"][ot] = 0
         end
         # assemble the observation vector
         y = Float64[]
         for d in db.data_groups[ct]
            o = observation(d, ot)
            # TODO: hack again!!!
            if ot == "E"
               o = copy(o)
               o[1] -= length(d) * E0
            end
            append!(y, o)
         end
         # get the lsq block
         block = reshape(db.kron_groups[ct][ot][:,:,Ibasis], :, length(Ibasis))
         # the error on this particular data-point
         e = block * c - y
         # compute the various errors for this ct/ot combination
         w = weighthook(ot, db.data_groups[ct][1])
         errs["rmse"][cn][ot] += w^2 * norm(e)^2
         errs["nrm2"][cn][ot] += w^2 * norm(y)^2
         errs["mae"][cn][ot]  += w * norm(e, 1)
         errs["nrm1"][cn][ot] += w * norm(y, 1)
         # append!(errs["allerr"][cn][ot], w * abs.(e) )
         errs["maxe"][cn][ot] = max(errs["maxe"][cn][ot], norm(e, Inf))
         errs["nrminf"][cn][ot] = max(errs["nrminf"][cn][ot], norm(y, Inf))
         lengths[cn][ot] += length(y)
         errs["rmse"]["set"][ot] += w^2 * norm(e)^2
         errs["nrm2"]["set"][ot] += w^2 * norm(y)^2
         errs["mae"]["set"][ot]  += w * norm(e, 1)
         errs["nrm1"]["set"][ot] += w * norm(y, 1)
         errs["maxe"]["set"][ot] = max(errs["maxe"]["set"][ot], norm(e, Inf))
         errs["nrminf"]["set"][ot] = max(errs["nrminf"]["set"][ot], norm(y, Inf))
         lengths["set"][ot] += length(y)
         # append!(errs["allerr"]["set"][ot], w * abs.(block * c - y) )
         # @show length(errs["allerr"]["set"][ot])
         # @show lengths["set"][ot]
      end
   end

   for cn in [confignames; "set"]
      for ot in keys(errs["rmse"][cn])
         len = lengths[cn][ot]
         errs["rmse"][cn][ot] = sqrt( errs["rmse"][cn][ot] / len )
         errs["nrm2"][cn][ot] = sqrt( errs["nrm2"][cn][ot] / len )
         errs["mae"][cn][ot] = errs["mae"][cn][ot] / len
         errs["nrm1"][cn][ot] = errs["nrm1"][cn][ot] / len
         errs["relrmse"][cn][ot] = errs["rmse"][cn][ot] / errs["nrm2"][cn][ot]
         # errs["allerr"][cn][ot] = quantile(errs["allerr"][cn][ot] ,range(0., stop=1., length=nb_points_cdf))
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
