# # Computation of the energy-force-virials for a given atomic position and IP
# function ipcomp(at,ot,IP)
#    if ot == "E"
#       return vec_obs(ot,energy(IP,at))
#    elseif ot == "F"
#       return vec_obs(ot,forces(IP,at))
#    elseif ot == "V"
#       return vec_obs(ot,virial(IP,at))
#    else
#       error("observation type not implemented")
#    end
# end
#

# function results_dict(data, IP; confignames = Colon(), pathname = "")
#    IPf = fast(IP)
#    allconfignames = unique(configname.(data))
#    allconfigtypes = unique(configtype.(data))
#    if confignames isa Colon
#       confignames = allconfignames
#    else
#       confignames = collect(confignames)
#    end
#    confignames = allconfignames
#
#    # create the dict for the results
#    results = Dict{String, Dict{String,Vector{Tuple{Vector{Float64},Vector{Float64},Int64}}}}()
#    for cn in confignames
#       results[cn] = Dict{String,Vector{Tuple{Vector{Float64},Vector{Float64},Int64}}}()
#    end
#
#    for (i,dat) in enumerate(data)
#        # print(".")
#       ct = configtype.(dat)
#       cn = configname.(dat)
#       if !(cn in confignames)
#          continue
#       end
#       for ot in ["E", "F", "V"]
#          #
#          if !hasobservation(dat, ot)
#             continue
#          end
#          if !haskey(results[cn], ot)
#             results[cn][ot] = Vector{Tuple{Vector{Float64},Vector{Float64},Int64}}[]
#          end
#          # Store exact data + approximation + indice in the data
#          push!(results[cn][ot], (observation(dat,ot) , (ipcomp(ASEAtoms(dat.at),ot,IPf)),i) )
#       end
#    end
#    if pathname != ""
#       save(pathname, results)
#    end
#    return results
# end
#
#
# function lsqerrors(res_dict, data; confignames = Colon(), nb_points_cdf = 40, E0 = nothing)
#    @assert E0 != nothing
#    if confignames isa Colon
#       confignames = collect(keys(res_dict))
#    else
#       confignames = collect(confignames)
#    end
#
#    # create the dict for the fit errors
#    errs = LsqErrors([confignames; "set"])
#    lengths = Dict{String, Dict{String, Int}}()
#    for cn in [confignames; "set"]
#       lengths[cn] = Dict{String, Int}()
#    end
#
#    for cn in keys(res_dict)
#       if !(cn in confignames)
#          continue
#       end
#       # loop through observation types in the current config type
#       for ot in ["E", "F", "V"]
#          if !haskey(res_dict[cn], ot)
#             continue
#          end
#          if !haskey(errs.rmse[cn], ot)
#             errs.rmse[cn][ot] = 0.0
#             errs.nrm2[cn][ot] = 0.0
#             errs.mae[cn][ot] = 0.0
#             errs.nrm1[cn][ot] = 0.0
#             errs.allerr[cn][ot] = Float64[]
#             errs.maxe[cn][ot] = 0.0
#             errs.nrminf[cn][ot] = 0.0
#             lengths[cn][ot] = 0
#          end
#          if !haskey(errs.rmse["set"], ot)
#             errs.rmse["set"][ot] = 0.0
#             errs.nrm2["set"][ot] = 0.0
#             errs.mae["set"][ot] = 0.0
#             errs.nrm1["set"][ot] = 0.0
#             errs.allerr["set"][ot] = Float64[]
#             errs.maxe["set"][ot] = 0.0
#             errs.nrminf["set"][ot] = 0.0
#             lengths["set"][ot] = 0
#          end
#          # loop over configurations
#          for d in res_dict[cn][ot]
#             i = d[3]
#             y = d[1]
#             yapprox = d[2]
#             if ot == "E"
#                y -= length(data[i].at) * E0
#                yapprox -= length(data[i].at) * E0
#             end
#             # get the approximate computation
#             # compute the various errors for this data/ot combination
#             w = weighthook(ot, data[i])
#             e = yapprox - y
#             errs.rmse[cn][ot] += w^2 * norm(e)^2
#             errs.nrm2[cn][ot] += w^2 * norm(y)^2
#             errs.mae[cn][ot]  += w * norm(e, 1)
#             errs.nrm1[cn][ot] += w * norm(y, 1)
#             append!(errs.allerr[cn][ot], w * abs.(e) )
#             errs.maxe[cn][ot] = max(errs.maxe[cn][ot], norm(e, Inf))
#             errs.nrminf[cn][ot] = max(errs.nrminf[cn][ot], norm(y, Inf))
#             lengths[cn][ot] += length(y)
#             errs.rmse["set"][ot] += w^2 * norm(e)^2
#             errs.nrm2["set"][ot] += w^2 * norm(y)^2
#             errs.mae["set"][ot]  += w * norm(e, 1)
#             errs.nrm1["set"][ot] += w * norm(y, 1)
#             errs.maxe["set"][ot] = max(errs.maxe["set"][ot], norm(e, Inf))
#             errs.nrminf["set"][ot] = max(errs.nrminf["set"][ot], norm(y, Inf))
#             lengths["set"][ot] += length(y)
#             append!(errs.allerr["set"][ot], w * abs.(e) )
#             # @show length(errs.allerr["set"][ot])
#             # @show lengths["set"][ot]
#          end
#
#       end
#    end
#
#    for cn in [confignames; "set"]
#       for ot in keys(errs.rmse[cn])
#          len = lengths[cn][ot]
#          errs.rmse[cn][ot] = sqrt( errs.rmse[cn][ot] / len )
#          errs.nrm2[cn][ot] = sqrt( errs.nrm2[cn][ot] / len )
#          errs.mae[cn][ot] = errs.mae[cn][ot] / len
#          errs.nrm1[cn][ot] = errs.nrm1[cn][ot] / len
#          errs.allerr[cn][ot] = quantile(errs.allerr[cn][ot] ,range(0., stop=1., length=nb_points_cdf))
#       end
#    end
#
#    # # ------- - scatter E -------
#    # push!(scatterE[ct][1], E_data/len)
#    # push!(scatterE[ct][2], E_fit/len)
#    # append!(scatterF[ct][1], f_data)
#    # append!(scatterF[ct][2], f_fit)
#
#    return errs
#
#
# end
#
# # function cdf_energy_forces(res_dict; nb_points_cdf = 40)
# #    absEerr = Float64[]
# #    absFerr = Float64[]
# #    absVerr = Float64[]
# #    for cn in collect(keys(res_dict))
# #       if haskey(res_dict[cn],"E")
# #          for i in 1:length(res_dict[cn]["E"])
# #             append!(absEerr, abs.(res_dict[cn]["E"][i][1] - res_dict[cn]["E"][i][2]) )
# #          end
# #       end
# #       if haskey(res_dict[cn],"F")
# #          for i in 1:length(res_dict[cn]["F"])
# #             append!(absFerr, abs.(res_dict[cn]["F"][i][1] - res_dict[cn]["F"][i][2]) )
# #          end
# #       end
# #    end
# #    return quantile(absEerr, collect(linspace(0.,1.,nb_points_cdf))), quantile(absFerr, collect(linspace(0.,1.,nb_points_cdf)))
# # end
#
