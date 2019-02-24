
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

export add_fits!, rmse, mae, rmse_table, mae_table

# ------------------------- USER INTERFACE FUNCTIONS -------------------------

rmse(data::Vector{Dat}; fitkey="fit") = _fiterrs(data; fitkey=fitkey, p=2)
# mae(data::Vector{Dat}; fitkey="fit") = _fiterrs(data; fitkey=fitkey, p=1)
# maxe(data::Vector{Dat}; fitkey="fit") = _fiterrs(data; fitkey=fitkey, p=Inf)

rmse_table(data::AbstractVector{Dat}) = rmse_table(rmse(data)...)
rmse_table(errs::Dict, errs_rel::Dict) = _err_table(errs, errs_rel, "RMSE")



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


# https://github.com/cortner/NBodyIPFitting.jl/blob/master/src/lsqerrors.jl
         # # assemble the observation vector
         # y = Float64[]
         # for d in db.data_groups[ct]
         #    o = observation(d, ot)
         #    # TODO: hack again!!!
         #    if ot == "E"
         #       o = copy(o)
         #       o[1] -= length(d) * E0
         #    end
         #    append!(y, o)
         # end
         # # get the lsq block
         # block = reshape(db.kron_groups[ct][ot][:,:,Ibasis], :, length(Ibasis))
         # # the error on this particular data-point
         # e = block * c - y
         #

# function _fiterrs(db::LsqDB, c, Ibasis; fitkey="fit", p=2, E0=nothing)
#    @assert E0 != nothing
#    confignames = unique(configname.(db.data_groups))
#    err = errdict(confignames)
#    relerr = errdict(confignames)
#    nrm = errdict(confignames)
#    len = errdict(confignames)
#
#    # assemble the errors
#    for ct in keys(db.data_groups)
#       cn = configname(ct)
#       for ot in keys(db.kron_groups[ct])  # the observation types
#          if !haskey(err[cn], ot)
#             len[cn][ot] = 0.0
#             nrm[cn][ot] = 0.0
#             err[cn][ot] = 0.0
#          end
#          if !haskey(err["set"], ot)
#             len["set"][ot] = 0.0
#             nrm["set"][ot] = 0.0
#             err["set"][ot] = 0.0
#          end
#          # get the lsq block
#          block = reshape(db.kron_groups[ct][ot][:,:,Ibasis], :, length(Ibasis))
#          # the error on this particular data-point
#          e = block * c - y
#          w = weighthook(ot, db.data_groups[ct][1])  # the default weight for this observation
#          exval = d.D[ot] * w
#          fitval = d.info[fitkey][ot] * w
#          len[cn][ot] += length(exval)
#          len["set"][ot] += length(exval)
#          nrm[cn][ot] += norm(exval, p)^p
#          nrm["set"][ot] += norm(exval, p)^p
#          err[cn][ot] += norm(exval-fitval, p)^p
#          err["set"][ot] += norm(exval-fitval, p)^p
#       end
#    end
#
#    # normalise correctly
#    for cn in keys(err)
#       for ot in keys(err[cn])
#          err[cn][ot] = (err[cn][ot]/len[cn][ot])^(1/p)
#          nrm[cn][ot] = (nrm[cn][ot]/len[cn][ot])^(1/p)
#          relerr[cn][ot] = err[cn][ot] / nrm[cn][ot]
#       end
#    end
#
#    return err, relerr
# end



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
   print("───────────────╂─────────────────┼─────────────────┼─────────────────\n")
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
      print("───────────────╂─────────────────┼─────────────────┼─────────────────\n")
      print(s_set)
   end
      print("━━━━━━━━━━━━━━━┷━━━━━━━━━━━━━━━━━┷━━━━━━━━━━━━━━━━━┷━━━━━━━━━━━━━━━━━\n")
end




end
