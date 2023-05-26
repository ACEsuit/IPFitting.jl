
"""
# `module Data`

Provides methods to read files containing simulation data. Primarily this is
intended to load `.xyz` files and convert them to JuLIP-compatible data:
```
data = IPFitting.Data.load_data("mydata.xyz")
```
where `mydata.xyz` contains multiple configurations, will read in those
configurations, then convert them into a `JuLIP.Atoms` object, extract
energies, forces, etc store each configuration in a `Dat` and return them
as a `Vector{Dat}`. This can then be stored using `JLD` or `JLD2`.

If `d::Dat`, then one can call `energy(d), forces(d)`, etc to extract the
loaded information.
"""
module Data

using JuLIP, ProgressMeter, FileIO, Printf, StatsBase, PrettyTables, DataFrames
using IPFitting: Dat, vec_obs, devec_obs, observation, hasobservation
using IPFitting.DataTypes
using StringDistances
import JuLIP: Atoms, energy, forces, virial
import Base: length, Dict
using PyCall

import ASE     # use ASE since this will have already figured out how to
               # load `ase` without problems
using ASE: ASEAtoms
#ase_io = ASE.ase_io

const ase_io = PyNULL()

function __init__()
   copy!(ase_io, pyimport_conda("ase.io", "ase", "rmg"))
end

export configtype, weight, load_data, truncate_string

Atoms(d::Dat) = d.at
length(d::Dat) = length(d.at)
configtype(d::Dat) = d.configtype
energy(d::Dat) = haskey(d.D, ENERGY) ? devec_obs(Val(:E), d.D[ENERGY]) : nothing
forces(d::Dat) = haskey(d.D, FORCES) ? devec_obs(Val(:F), d.D[FORCES]) : nothing
virial(d::Dat) = haskey(d.D, VIRIAL) ? devec_obs(Val(:V), d.D[VIRIAL]) : nothing


function read_energy(atpy, energy_key)
   for key in keys(atpy.info)
      if lowercase(key) == lowercase(energy_key)
         return atpy.info[key]
      end
   end
   return nothing
end

function read_forces(atpy, force_key)
   for key in keys(atpy.arrays)
      if lowercase(key) == lowercase(force_key)
         F = collect(vecs(collect(atpy.arrays[key]')))
         return F
      end
   end
   return nothing
end

function read_virial(atpy, virial_key)
   for key in keys(atpy.info)
      if lowercase(key) == lowercase(virial_key)
         return JMat(atpy.info[key]...)
      end
   end
   # for key in keys(atpy.info)
   #    if compare(lowercase(key), "stress", Levenshtein()) > 0.8
   #       @info("No virial key found: converting stress to virial!")
   #       return -1.0 * JMat(atpy.info[key]...) * atpy.get_volume()
   #    end
   # end
   return nothing
end

function read_configtype(atpy)
   for key in keys(atpy.info)
      if lowercase(key) == "config_type"
         return atpy.info[key]
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

count_scalars(::Nothing) = 0 
count_scalars(::Number) = 1 
count_scalars(x::AbstractArray{<:Number}) = length(x)
count_scalars(x::AbstractArray{<: AbstractArray{<: Number}}) = 
      sum(count_scalars, x)


function keys_info(key_dict)
   s = "keys found: "
   for key in sort(collect(keys(key_dict)))
      n = key_dict[key]
      s *= "\"$(key)\" [$n], "
   end
   return s[1:end-2]
end

read_Atoms(atpy) = Atoms(ASEAtoms( atpy ))

"""
```
function read_xyz(fname; verbose=true, index = ":",
                         include = nothing, exclude = nothing)
```
Loads the atoms objects contained in an xyz file, attempts to read the
DFT data stored inside and returns a `Vector{Dat}`.
"""
function read_xyz(fname; energy_key = "dft_energy", force_key = "dft_force", virial_key = "dft_virial", auto=false, verbose=true, index = ":",
                         include = nothing, exclude = nothing)
   if verbose
      println("Reading in $fname ...")
   end
   at_list = ase_io.read(fname, index=index)
   #at_list = ase_io_read(fname, index=index)
   data = Vector{Dat}(undef, length(at_list))
   idx = 0
   if verbose
      println("Processing data ...")
      # tic()
   end

   info_dict = countmap(collect(Iterators.flatten([collect(keys(at.info)) for at in at_list])))
   array_dict = countmap(collect(Iterators.flatten([collect(keys(at.arrays)) for at in at_list])))

   if auto
      energy_key = findmax([[compare(energy_key, key, Levenshtein()), key] for key in keys(info_dict)])[1][2]
      force_key = findmax([[compare(force_key, key, Levenshtein()), key] for key in keys(array_dict)])[1][2]
      virial_key = findmax([[compare(virial_key, key, Levenshtein()), key] for key in keys(info_dict)])[1][2]
   end

   @info("Keys used: E => \"$(energy_key)\", F => \"$(force_key)\", V => \"$(virial_key)\"")

   if verbose
      for (key,val) in Dict("Info " => info_dict, "Array " => array_dict)
         @info(key * keys_info(val))
      end
   end

   ct_dict = Dict()
   ct_datf = DataFrame()#A = Any[], B = Any[], C = Any[], D = Any[], E = Any[], F = Any[])

   dt = verbose ? 1.0 : 0.0
   @showprogress dt for (i,atpy) in enumerate(at_list)

      # get the config type and decide whether to keep or skip this config
      ct = read_configtype(atpy)
      if ct == nothing
         verbose && @warn("$(idx+1) has no config_type")
         ct = "nothing"
      end
      if exclude != nothing
         if ct in exclude
            continue
         end
      end
      if include != nothing
         if !(ct in include)
            continue
         end
      end

      idx += 1
      at = read_Atoms(atpy)
      E = read_energy(atpy, energy_key)
      F = read_forces(atpy, force_key)
      V = read_virial(atpy, virial_key)

      E_c = count_scalars(E)
      F_c = count_scalars(F)
      V_c = count_scalars(V)

      if ct in names(ct_datf)
         ct_datf[!, Symbol(ct)] .+= [1, length(atpy), E_c, F_c, V_c]
      else
         ct_datf[!, Symbol(ct)] = [1, length(atpy), E_c, F_c, V_c]
      end

      EFV = ""
      (E != nothing) && (EFV *= "E")
      (F != nothing) && (EFV *= "F")
      (V != nothing) && (EFV *= "V")

      data[idx] = Dat( at,
                       ct;
                       E = E, F = F, V = V )

      if (haskey(atpy.info, "energy_sigma")
            || haskey(atpy.info, "force_sigma")
            || haskey(atpy.info, "virial_sigma"))
         data[idx].info["W"] = Dict{String,Float64}()
      end
      if haskey(atpy.info, "energy_sigma")
         data[idx].info["W"]["E"] = 1.0/atpy.info["energy_sigma"]
      end
      if haskey(atpy.info, "force_sigma")
         data[idx].info["W"]["F"] = 1.0/atpy.info["force_sigma"]
      end
      if haskey(atpy.info, "virial_sigma")
         data[idx].info["W"]["V"] = 1.0/atpy.info["virial_sigma"]
      end

   end
   # verbose && toc()

   ct_datf_trans = DataFrame([[names(ct_datf)]; collect.(eachrow(ct_datf))], [:column; Symbol.(axes(ct_datf, 1))]) #apparently only "simple" way to take a dataframe transpose......
   rename!(ct_datf_trans, Symbol.(["config_type", "#cfgs", "#envs", "#E", "#F", "#V"]))

   totals = sum.(eachcol(ct_datf_trans[:, 2:end]))
   missings = [0, 0, abs(totals[1] - totals[3]), abs(3*totals[2] - totals[4]), abs(9*totals[1] - totals[5])]
   push!(ct_datf_trans, vcat(["total", totals]...))
   push!(ct_datf_trans, vcat(["missing", missings]...))

   pretty_table(ct_datf_trans, crop = :horizontal, body_hlines = [length(ct_datf_trans[!, 1])-2])
   #
   # totals = [0,0,0,0,0]
   #
   # print("┏━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┓\n")
   # print("┃ config type  ┃ #cfgs ┃ #envs  ┃      #E      ┃      #F      ┃      #V      ┃\n")
   # print("┠──────────────╂───────╂────────╂──────────────╂──────────────╂──────────────┨\n")
   # for ct in sort(collect(keys(ct_dict)))
   #    s = @sprintf("┃ %12s ┃ %5s ┃ %6s ┃ %6s [%4s]┃ %6s [%4s]┃ %6s [%4s]┃\n",
   #       truncate_string(ct, 12), ct_dict[ct][1], ct_dict[ct][2],
   #          ct_dict[ct][3], abs(ct_dict[ct][3] - ct_dict[ct][1]),
   #          ct_dict[ct][4], abs(ct_dict[ct][4] - 3*ct_dict[ct][2]),
   #          ct_dict[ct][5], abs(ct_dict[ct][5] - 9*ct_dict[ct][1]))
   #    print(s)
   #    for i in 1:5
   #       totals[i] += ct_dict[ct][i]
   #    end
   # end
   #
   # print("┠──────────────╂───────╂────────╂──────────────╂──────────────╂──────────────┨\n")
   # s = @sprintf("┃    totals    ┃ %5s ┃ %6s ┃ %12s ┃ %12s ┃ %12s ┃\n",
   #          totals[1], totals[2], totals[3], totals[4], totals[5])
   # print(s)
   # print("┗━━━━━━━━━━━━━━┻━━━━━━━┻━━━━━━━━┻━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━┛\n")

   return data[1:idx]
end



function load_data(fname; kwargs...)
   if fname[end-3:end] == ".xyz"
      return read_xyz(fname; kwargs...)
   elseif fname[end-3:end] == ".jld" || fname[end-4:end] == ".jld2"
      return load(fname)
   end
   error("unknown file format $(fname[end-3:end])")
end


end
