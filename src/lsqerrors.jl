
module Errors

using JuLIP: Atoms, energy, forces
using NBodyIPFitting: LsqDB, Dat, configtype, weighthook
using NBodyIPFitting.Data: observation, hasobservation, configname

export lsqerrors, table, table_relative, table_absolute, relerr_table, abserr_table
# , scatter_E, scatter_F


struct LsqErrors
   rmse::Dict{String, Dict{String,Float64}}
    mae::Dict{String, Dict{String,Float64}}
   nrm2::Dict{String, Dict{String,Float64}}
   nrm1::Dict{String, Dict{String,Float64}}
   # scatterE::Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}
   # scatterF::Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}
end

Base.Dict(errs::LsqErrors) =
   Dict("rmse" => errs.rmse, "mae" => errs.mae,
        "nrm2" => errs.nrm2, "nrm1" => errs.nrm1)

LsqErrors(errs::Dict) =
   LsqErrors(errs["rmse"], errs["mae"], errs["nrm2"], errs["nrm1"])

rmse(errs::LsqErrors, ct, ot) =
      haskey(errs.rmse[ct], ot) ? errs.rmse[ct][ot] : NaN
relrmse(errs::LsqErrors, ct, ot) =
      haskey(errs.rmse[ct], ot) ? errs.rmse[ct][ot]/errs.nrm2[ct][ot] : NaN
mae(errs::LsqErrors, ct, ot) =
      haskey(errs.mae[ct], ot) ? errs.mae[ct][ot] : NaN
relmae(errs::LsqErrors, ct, ot) =
      haskey(errs.mae[ct], ot) ? errs.mae[ct][ot]/errs.nrm1[ct][ot] : NaN

rmse(errs::Dict, ct, ot) = rmse(LsqErrors(errs))
relrmse(errs::Dict, ct, ot) = relrmse(LsqErrors(errs))
mae(errs::Dict, ct, ot) = mae(LsqErrors(errs))
relmae(errs::Dict, ct, ot) = relmae(LsqErrors(errs))

function errdict(configtypes)
   D = Dict{String, Dict{String, Float64}}()
   for ct in configtypes
      D[ct] = Dict{String, Float64}()
   end
   return D
end

LsqErrors(configtypes) =
         LsqErrors(errdict(configtypes), errdict(configtypes),
                   errdict(configtypes), errdict(configtypes) )

function lsqerrors(db, c, Ibasis; confignames = Colon(), E0 = nothing)
   if confignames isa Colon
      confignames = configname.(collect(keys(db.data_groups)))
   else
      confignames = collect(confignames)
   end
   @assert E0 != nothing

   # create the dict for the fit errors
   errs = LsqErrors([confignames; "set"])
   lengths = Dict{String, Dict{String, Int}}()
   for cn in [confignames; "set"]
      lengths[cn] = Dict{String, Int}()
   end


   # scatterE = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
   # scatterF = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()

   idx = 0
   for ct in keys(db.data_groups)
      cn = configname(ct)
      if !(cn in confignames)
         continue
      end
      # loop through observation types in the current config type
      for ot in keys(db.kron_groups[ct])
         if !haskey(errs.rmse[cn], ot)
            errs.rmse[cn][ot] = 0.0
            errs.nrm2[cn][ot] = 0.0
            errs.mae[cn][ot] = 0.0
            errs.nrm1[cn][ot] = 0.0
            lengths[cn][ot] = 0
            errs.rmse["set"][ot] = 0.0
            errs.nrm2["set"][ot] = 0.0
            errs.mae["set"][ot] = 0.0
            errs.nrm1["set"][ot] = 0.0
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
         # compute the various errors for this ct/ot combination
         w = weighthook(ot, db.data_groups[ct][1])
         errs.rmse[cn][ot] += w^2 * norm(block * c - y)^2
         errs.nrm2[cn][ot] += w^2 * norm(y)^2
         errs.mae[cn][ot]  += w * norm(block * c - y, 1)
         errs.nrm1[cn][ot] += w * norm(y, 1)
         lengths[cn][ot] += length(y)
         errs.rmse["set"][ot] += w^2 * norm(block * c - y)^2
         errs.nrm2["set"][ot] += w^2 * norm(y)^2
         errs.mae["set"][ot]  += w * norm(block * c - y, 1)
         errs.nrm1["set"][ot] += w * norm(y, 1)
         lengths["set"][ot] += length(y)
      end
   end

   for cn in [confignames; "set"]
      for ot in keys(errs.rmse[cn])
         len = lengths[cn][ot]
         errs.rmse[cn][ot] = sqrt( errs.rmse[cn][ot] / len )
         errs.nrm2[cn][ot] = sqrt( errs.nrm2[cn][ot] / len )
         errs.mae[cn][ot] = errs.mae[cn][ot] / len
         errs.nrm1[cn][ot] = errs.nrm1[cn][ot] / len
      end
   end

   # # ------- - scatter E -------
   # push!(scatterE[ct][1], E_data/len)
   # push!(scatterE[ct][2], E_fit/len)
   # append!(scatterF[ct][1], f_data)
   # append!(scatterF[ct][2], f_fit)

   return errs
end


function table(errs; relative=false)
   if relative
      table_relative(errs)
   else
      table_absolute(errs)
   end
end

relerr_table(fit_info) = table_relative(fit_info["errors"])
abserr_table(fit_info) = table_absolute(fit_info["errors"])

table_relative(errs::Dict) = table_relative(LsqErrors(errs))
table_absolute(errs::Dict) = table_absolute(LsqErrors(errs))

function table_relative(errs::LsqErrors)
   print("-----------------------------------------------------------------\n")
   print("               ||           RMSE        ||           MAE      \n")
   print("  config type  || E [%] | F [%] | V [%] || E [%] | F [%] | V [%] \n")
   print("---------------||-------|-------|-------||-------|-------|-------\n")
   s_set = ""
   for ct in keys(errs.rmse)  # ct in configtypes
      lct = min(length(ct), 13)
      s = @sprintf(" %13s || %5.2f | %5.2f | %5.2f || %5.2f | %5.2f | %5.2f \n",
         ct[1:lct],
         100*relrmse(errs, ct, "E"),
         100*relrmse(errs, ct, "F"),
         100*relrmse(errs, ct, "V"),
         100* relmae(errs, ct, "E"),
         100* relmae(errs, ct, "F"),
         100* relmae(errs, ct, "V") )
      if ct == "set"
         s_set = s
      else
         print(s)
      end
   end
   # print("---------------||-------|-------|-------||-------|-------|-------\n")
   # print(s_set)
   print("-----------------------------------------------------------------\n")
end


function table_absolute(errs::LsqErrors)
   print("-----------------------------------------------------------------------------\n")
   print("               ||            RMSE             ||             MAE        \n")
   print("  config type  ||  E [eV] | F[eV/A] | V[eV/A] ||  E [eV] | F[eV/A] | V[eV/A] \n")
   print("---------------||---------|---------|---------||---------|---------|---------\n")
   s_set = ""
   for ct in keys(errs.rmse)  # ct in configtypes

      lct = min(length(ct), 16)
      s = @sprintf(" %16s || %7.4f | %7.4f | %7.4f || %7.4f | %7.4f | %7.4f \n",
                   ct[1:lct],
                   rmse(errs, ct, "E"), rmse(errs, ct, "F"), rmse(errs, ct, "V"),
                    mae(errs, ct, "E"),  mae(errs, ct, "F"),  mae(errs, ct, "V") )
      if ct == "set"
         s_set = s
      else
         print(s)
      end
   end
   # print("---------------||---------|---------|---------||---------|---------|---------\n")
   # print(s_set)
   print("-----------------------------------------------------------------------------\n")
end


# function scatter_E(errs::LsqErrors; s=1)
#    D = deepcopy(errs.scatterE)
#    for (key, val) in D
#       E_data, E_fit = val
#       D[key] = (E_data[1:s:end], E_fit[1:s:end])
#    end
#    return D
# end
#
# function scatter_F(errs::LsqErrors; s=20)
#    D = deepcopy(errs.scatterF)
#    for (key, val) in D
#       F_data, F_fit = val
#       D[key] = (F_data[1:s:end], F_fit[1:s:end])
#    end
#    return D
# end

end
