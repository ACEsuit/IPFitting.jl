using JuLIP
using NBodyIPs.Data: Dat, weight, config_type, energy, forces

export fiterrors, scatter_data, scatter_E, scatter_F


_fiterrsdict() = Dict("E-RMS" => 0.0, "F-RMS" => 0.0, "V-RMS" => 0.0,
                      "E-MAE" => 0.0, "F-MAE" => 0.0, "V-MAE" => 0.0)
_cnterrsdict() = Dict("E" => 0, "F" => 0, "V" => 0)

struct FitErrors
   errs::Dict{String, Dict{String,Float64}}
   nrms::Dict{String, Dict{String,Float64}}
   scatterE::Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}
   scatterF::Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}
end

function fiterrors(lsq, c, Ibasis; include=nothing, exclude=nothing)
   include = analyse_include_exclude(lsq, include, exclude)

   E0 = lsq.basis[1]()
   # create the dict for the fit errors
   errs = Dict( "set" => _fiterrsdict() )
   # count number of configurations
   num = Dict("set" => _cnterrsdict())
   obs = Dict("set" => _fiterrsdict())
   scatterE = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
   scatterF = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()

   for ct in include
      errs[ct] = _fiterrsdict()
      num[ct] = _cnterrsdict()
      obs[ct] = _fiterrsdict()
      scatterE[ct] = (Vector{Float64}(), Vector{Float64}())
      scatterF[ct] = (Vector{Float64}(), Vector{Float64}())
   end

   idx = 0
   for d in lsq.data
      ct = config_type(d)
      E_data, F_data, V_data, len = energy(d), forces(d), virial(d), length(d)
      E_data -= E0 * len

      if !(ct in include)
         idx += 1 # E
         if F_data != nothing
            idx += 3 * len
         end
         if virial(d) != nothing
            idx += length(_IS)
         end
         continue
      end


      # ----- Energy error -----
      E_fit = dot(lsq.Ψ[idx+1, Ibasis], c)
      Erms = (E_data - E_fit)^2 / len^2
      Emae = abs(E_data - E_fit) / len
      errs[ct]["E-RMS"] += Erms
      obs[ct]["E-RMS"] += (E_data / len)^2
      errs[ct]["E-MAE"] += Emae
      obs[ct]["E-MAE"] += abs(E_data / len)
      num[ct]["E"] += 1
      # - - - - - - - - - - - - - - - -
      errs["set"]["E-RMS"] += Erms
      obs["set"]["E-RMS"] += (E_data / len)^2
      errs["set"]["E-MAE"] += Emae
      obs["set"]["E-MAE"] += abs(E_data / len)
      num["set"]["E"] += 1
      idx += 1
      # ------- - scatter E -------
      push!(scatterE[ct][1], E_data/len)
      push!(scatterE[ct][2], E_fit/len)

      # ----- Force error -------
      if (F_data != nothing) && (len > 1)
         f_data = mat(F_data)[:]
         f_fit = lsq.Ψ[(idx+1):(idx+3*len), Ibasis] * c
         Frms = norm(f_data - f_fit)^2
         Fmae = norm(f_data - f_fit, 1)
         errs[ct]["F-RMS"] += Frms
         obs[ct]["F-RMS"] += norm(f_data)^2
         errs[ct]["F-MAE"] += Fmae
         obs[ct]["F-MAE"] += norm(f_data, 1)
         num[ct]["F"] += 3 * len
         # - - - - - - - - - - - - - - - -
         errs["set"]["F-RMS"] += Frms
         obs["set"]["F-RMS"] += norm(f_data)^2
         errs["set"]["F-MAE"] += Fmae
         obs["set"]["F-MAE"] += norm(f_data, 1)
         num["set"]["F"] += 3 * len
         idx += 3 * len
         # -------- scatter F ------
         append!(scatterF[ct][1], f_data)
         append!(scatterF[ct][2], f_fit)
      end

      # -------- stress errors ---------
      if V_data != nothing
         V_fit = lsq.Ψ[(idx+1):(idx+length(_IS)),Ibasis] * c
         V_err = norm(V_fit - V_data[_IS], Inf) / len  # replace with operator norm on matrix
         V_nrm = norm(V_data[_IS], Inf)
         errs[ct]["V-RMS"] += V_err^2
         obs[ct]["V-RMS"] += V_nrm^2
         errs[ct]["V-MAE"] += V_err
         obs[ct]["V-MAE"] += V_nrm
         num[ct]["V"] += 1
         # - - - - - - - - - - - - - - - -
         errs["set"]["V-RMS"] += V_err^2
         obs["set"]["V-RMS"] += V_nrm^2
         errs["set"]["V-MAE"] += V_err
         obs["set"]["V-MAE"] += V_nrm
         num["set"]["V"] += 1
         idx += length(_IS)
      end
   end

   # NORMALISE
   for key in keys(errs)
      nE = num[key]["E"]
      nF = num[key]["F"]
      nV = num[key]["V"]
      errs[key]["E-RMS"] = sqrt(errs[key]["E-RMS"] / nE)
      obs[key]["E-RMS"] = sqrt(obs[key]["E-RMS"] / nE)
      errs[key]["F-RMS"] = sqrt(errs[key]["F-RMS"] / nF)
      obs[key]["F-RMS"] = sqrt(obs[key]["F-RMS"] / nF)
      errs[key]["V-RMS"] = sqrt(errs[key]["V-RMS"] / nV)
      obs[key]["V-RMS"] = sqrt(obs[key]["V-RMS"] / nV)
      errs[key]["E-MAE"] = errs[key]["E-MAE"] / nE
      obs[key]["E-MAE"] = obs[key]["E-MAE"] / nE
      errs[key]["F-MAE"] = errs[key]["F-MAE"] / nF
      obs[key]["F-MAE"] = obs[key]["F-MAE"] / nF
      errs[key]["V-MAE"] = errs[key]["V-MAE"] / nV
      obs[key]["V-MAE"] = obs[key]["V-MAE"] / nV
   end

   return FitErrors(errs, obs, scatterE, scatterF)
end


Base.info(errs::FitErrors; kwargs...) = table(errs; kwargs...)

function table(errs::FitErrors; relative=false)
   if relative
      table_relative(errs)
   else
      table_absolute(errs)
   end
end

function table_relative(errs)
   print("---------------------------------------------------------------\n")
   print("             ||           RMSE        ||           MAE      \n")
   print(" config type || E [%] | F [%] | V [%] || E [%] | F [%] | V [%] \n")
   print("-------------||-------|-------|-------||-------|-------|-------\n")
   s_set = ""
   nrms = errs.nrms
   for (key, D) in errs.errs
      nrm = nrms[key]
      lkey = min(length(key), 11)
      s = @sprintf(" %11s || %5.2f | %5.2f | %5.2f || %5.2f | %5.2f | %5.2f \n",
         key[1:lkey],
         100*D["E-RMS"]/nrm["E-RMS"], 100*D["F-RMS"]/nrm["F-RMS"], 100*D["V-RMS"]/nrm["V-RMS"],
         100*D["E-MAE"]/nrm["E-MAE"], 100*D["F-MAE"]/nrm["F-MAE"], 100*D["V-MAE"]/nrm["V-MAE"])
      if key == "set"
         s_set = s
      else
         print(s)
      end
   end
   print("-------------||-------|-------|-------||-------|-------|-------\n")
   print(s_set)
   print("---------------------------------------------------------------\n")
end


function table_absolute(errs::FitErrors)
   print("---------------------------------------------------------------------------\n")
   print("             ||            RMSE             ||             MAE        \n")
   print(" config type ||  E [eV] | F[eV/A] | V[eV/A] ||  E [eV] | F[eV/A] | V[eV/A] \n")
   print("-------------||---------|---------|---------||---------|---------|---------\n")
   s_set = ""
   for (key, D) in errs.errs
      lkey = min(length(key), 11)
      s = @sprintf(" %11s || %7.4f | %7.4f | %7.4f || %7.4f | %7.4f | %7.4f \n",
                   key[1:lkey],
                   D["E-RMS"], D["F-RMS"], D["V-RMS"],
                   D["E-RMS"], D["E-MAE"], D["V-MAE"] )
      if key == "set"
         s_set = s
      else
         print(s)
      end
   end
   print("-------------||---------|---------|---------||---------|---------|---------\n")
   print(s_set)
   print("---------------------------------------------------------------------------\n")
end




function scatter_E(errs::FitErrors; s=1)
   D = deepcopy(errs.scatterE)
   for (key, val) in D
      E_data, E_fit = val
      D[key] = (E_data[1:s:end], E_fit[1:s:end])
   end
   return D
end

function scatter_F(errs::FitErrors; s=20)
   D = deepcopy(errs.scatterF)
   for (key, val) in D
      F_data, F_fit = val
      D[key] = (F_data[1:s:end], F_fit[1:s:end])
   end
   return D
end
