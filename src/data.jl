
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

using JuLIP, ProgressMeter, FileIO, Printf
using IPFitting: Dat, vec_obs, devec_obs, observation, hasobservation
using IPFitting.DataTypes
import JuLIP: Atoms, energy, forces, virial
import Base: length, Dict

import ASE     # use ASE since this will have already figured out how to
               # load `ase` without problems
using ASE: ASEAtoms
ase_io = ASE.ase_io

export configtype, weight, load_data

Atoms(d::Dat) = d.at
length(d::Dat) = length(d.at)
configtype(d::Dat) = d.configtype
energy(d::Dat) = haskey(d.D, ENERGY) ? devec_obs(Val(:E), d.D[ENERGY]) : nothing
forces(d::Dat) = haskey(d.D, FORCES) ? devec_obs(Val(:F), d.D[FORCES]) : nothing
virial(d::Dat) = haskey(d.D, VIRIAL) ? devec_obs(Val(:V), d.D[VIRIAL]) : nothing


function read_energy(atpy, energy_key)
   for key in keys(atpy.info)
      if lowercase(key) == energy_key
         return atpy.info[key]
      end
   end
   return nothing
end

function read_forces(atpy, force_key)
   for key in keys(atpy.arrays)
      if lowercase(key) == force_key
         return atpy.arrays[key]' |> vecs
      end
   end
   return nothing
end

function read_virial(atpy, virial_key)
   for key in keys(atpy.info)
      if lowercase(key) == virial_key
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

read_Atoms(atpy) = Atoms(ASEAtoms( atpy ))

"""
```
function read_xyz(fname; verbose=true, index = ":",
                         include = nothing, exclude = nothing)
```
Loads the atoms objects contained in an xyz file, attempts to read the
DFT data stored inside and returns a `Vector{Dat}`.
"""
function read_xyz(fname; energy_key = "", force_key = "", virial_key = "", verbose=true, index = ":",
                         include = nothing, exclude = nothing)
   if verbose
      println("Reading in $fname ...")
   end
   at_list = ase_io.read(fname, index=index)
   data = Vector{Dat}(undef, length(at_list))
   idx = 0
   if verbose
      println("Processing data ...")
      # tic()
   end

   #E_f, F_f, V_f count the # of energy/force/virial components found

   E_f = 0
   F_f = 0
   V_f = 0

   config_dict = Dict()

   dt = verbose ? 1.0 : 0.0
   @showprogress dt for (i,atpy) in enumerate(at_list)

      # get the config type and decide whether to keep or skip this config
      config_type = read_configtype(atpy)
      if config_type == nothing
         verbose && @warn("$(idx+1) has no config_type")
         config_type = "nothing"
      end
      if exclude != nothing
         if config_type in exclude
            continue
         end
      end
      if include != nothing
         if !(config_type in include)
            continue
         end
      end

      if config_type in keys(config_dict)
         config_dict[config_type] .+= [1,length(atpy)]
      else
         config_dict[config_type] = [1,length(atpy)]
      end

      idx += 1
      at = read_Atoms(atpy)
      E = read_energy(atpy, energy_key)
      F = read_forces(atpy, force_key)
      V = read_virial(atpy, virial_key)

      if E != nothing
         E_f += length(hcat(E...))
      else
         @warn("Energy not found for config: " * string(i))
      end

      if F != nothing
         F_f += length(hcat(F...))
      else
         @warn("Force not found for config: " * string(i))
      end

      if V != nothing
         V_f += length(hcat(V...))
      else
         @warn("Virial not found for config: " * string(i))
      end

      EFV = ""
      (E != nothing) && (EFV *= "E")
      (F != nothing) && (EFV *= "F")
      (V != nothing) && (EFV *= "V")

      data[idx] = Dat( at,
                       config_type;
                       E = E, F = F, V = V )
   end
   # verbose && toc()
   print("┏━━━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━┓\n")
   print("┃ config type  ┃  # cfgs  ┃  # envs  ┃\n")
   print("┠──────────────╂──────────╂──────────┨\n")
   for config_type in sort(collect(keys(config_dict)))
      s = @sprintf("┃ %12s ┃ %8s ┃ %8s ┃\n", config_type, config_dict[config_type][1], config_dict[config_type][2])
      print(s)
   end
   print("┗━━━━━━━━━━━━━━┻━━━━━━━━━━┻━━━━━━━━━━┛\n")

   @info("Configurations found in .xyz file: " * string(length(at_list)))
   @info("Energies found [key: \"" * energy_key * "\"]: " * string(E_f))
   @info("Forces found [key: \"" * force_key * "\"]: " * string(F_f))
   @info("Virials found [key: \"" * virial_key * "\"]: " * string(V_f))

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
