
module Errors

using JuLIP: Atoms, energy, forces
using NBodyIPFitting: LsqDB, Dat, configtype, weighthook
using NBodyIPFitting.Data: observation

export lsqerrors, table, table_relative, table_absolute
# , scatter_E, scatter_F


struct LsqErrors
   rmse::Dict{String, Dict{String,Float64}}
    mae::Dict{String, Dict{String,Float64}}
   nrm2::Dict{String, Dict{String,Float64}}
   nrm1::Dict{String, Dict{String,Float64}}
   # scatterE::Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}
   # scatterF::Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}
end

rmse(errs::LsqErrors, ct, ot) =
      haskey(errs.rmse[ct], ot) ? errs.rmse[ct][ot] : NaN
relrmse(errs::LsqErrors, ct, ot) =
      haskey(errs.rmse[ct], ot) ? errs.rmse[ct][ot]/errs.nrm2[ct][ot] : NaN
mae(errs::LsqErrors, ct, ot) =
      haskey(errs.mae[ct], ot) ? errs.mae[ct][ot] : NaN
relmae(errs::LsqErrors, ct, ot) =
      haskey(errs.mae[ct], ot) ? errs.mae[ct][ot]/errs.nrm1[ct][ot] : NaN

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

function lsqerrors(db, c, Ibasis; configtypes = Colon(), E0 = nothing)
   if configtypes isa Colon
      configtypes = collect(keys(db.data_groups))
   end
   @assert E0 != nothing

   # create the dict for the fit errors
   errs = LsqErrors(configtypes)
   # scatterE = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
   # scatterF = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()

   idx = 0
   for ct in configtypes
      # loop through observation types in the current config type
      for ot in keys(db.kron_groups[ct])
         # assemble the observation vector
         y = Float64[]
         for d in db.data_groups[ct]
            append!(y, observation(d, ot))
         end
         # get the lsq block
         block = reshape(db.kron_groups[ct][ot][:,:,Ibasis], :, length(Ibasis))
         # compute the various errors for this ct/ot combination
         w = weighthook(ot, db.data_groups[ct][1])
         errs.rmse[ct][ot] = w * norm(block * c - y) / sqrt(length(y))
         errs.nrm2[ct][ot] = w * norm(y) / sqrt(length(y))
         errs.mae[ct][ot]  = w * norm(block * c - y, 1) / length(y)
         errs.nrm1[ct][ot] = w * norm(y, 1) / length(y)
      end
   end

   # # ------- - scatter E -------
   # push!(scatterE[ct][1], E_data/len)
   # push!(scatterE[ct][2], E_fit/len)
   # append!(scatterF[ct][1], f_data)
   # append!(scatterF[ct][2], f_fit)

   return errs
end


function table(errs::LsqErrors; relative=false)
   if relative
      table_relative(errs)
   else
      table_absolute(errs)
   end
end

function table_relative(errs)
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

      lct = min(length(ct), 13)
      s = @sprintf(" %13s || %7.4f | %7.4f | %7.4f || %7.4f | %7.4f | %7.4f \n",
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
