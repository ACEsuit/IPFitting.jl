
"""
# `module Data`

Provides methods to read files containing simulation data. Primarily this is
intended to load `.xyz` files and convert them to JuLIP/NBodyIPs-compatible data:
```
data = NBodyIPs.Data.load_data("mydata.xyz")
```
where `mydata.xyz` contains multiple configurations, will read in those
configurations, then convert them into a `JuLIP.Atoms` object, extract
energies, forces, etc store each configuration in a `Dat` and return them
as a `Vector{Dat}`. This can then be stored using `JLD` or `JLD2`.

If `d::Dat`, then one can call `energy(d), forces(d)`, etc to extract the
loaded information.
"""
module Data

using JuLIP, ProgressMeter, FileIO
using NBodyIPFitting: Dat, vec, devec
using NBodyIPFitting.DataTypes
import JuLIP: Atoms, energy, forces, virial
import Base: length, Dict

import ASE     # use ASE since this will have already figured out how to
               # load `ase` without problems
ase_io = ASE.ase_io

export configtype, configname, weight, load_data

Atoms(d) = d.at
length(d::Dat) = length(d.at)
configtype(d::Dat) = d.configtype
configname(d::Dat) = configname(configtype(d))
configname(s::String) = match(r"[a-z,A-Z,0-9,_,%]*", s).match
energy(d::Dat) = haskey(d.D, ENERGY) ? devec(Val(:E), d.D[ENERGY]) : nothing
forces(d::Dat) = haskey(d.D, FORCES) ? devec(Val(:F), d.D[FORCES]) : nothing
virial(d::Dat) = haskey(d.D, VIRIAL) ? devec(Val(:V), d.D[VIRIAL]) : nothing
observation(d::Dat, key::String) = d.D[key]
hasobservation(d::Dat, key::String) = haskey(d.D, key)


function read_energy(atpy)
   for key in keys(atpy[:info])
      if lowercase(key) == "dft_energy"
         return atpy[:info][key]
      end
   end
   try
      return atpy[:get_potential_energy]()
   catch
   end
   return nothing
end

function read_forces(atpy)
   for key in keys(atpy[:arrays])
      if lowercase(key) == "dft_force"
         return atpy[:arrays][key]' |> vecs
      end
   end
   try
      return atpy[:get_array]("force")' |> vecs
   catch
   end
   return nothing
end

function read_virial(atpy)
   for key in keys(atpy[:info])
      if lowercase(key) == "dft_virial"
         return JMat(atpy[:info][key]...)
      end
   end
   if haskey(atpy[:info], "virial")
      return JMat(atpy[:info]["virial"]...)
   end
   return nothing
end

function read_configtype(atpy)
   if haskey(atpy[:info], "config_type")
      return atpy[:info]["config_type"]
   elseif haskey(atpy[:info], "configtype")
      return atpy[:info]["configtype"]
   end
   return nothing
end

read_Atoms(atpy) = Atoms( vecs(atpy[:get_positions]()'),
                          vecs(atpy[:get_momenta]()'),
                          atpy[:get_masses](),
                          atpy[:get_atomic_numbers](),
                          atpy[:get_cell](),
                          tuple(atpy[:get_pbc]()...) )

"""
```
function read_xyz(fname; verbose=true, index = ":",
                         include = nothing, exclude = nothing)
```
Loads the atoms objects contained in an xyz file, attempts to read the
DFT data stored inside and returns a `Vector{Dat}`.
"""
function read_xyz(fname; verbose=true, index = ":",
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
   dt = verbose ? 1.0 : 0.0
   @showprogress dt for atpy in at_list

      # get the config type and decide whether to keep or skip this config
      config_type = read_configtype(atpy)
      if config_type == nothing
         verbose && warn("$(idx+1) has no config_type")
         config_type = ""
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

      idx += 1
      at = read_Atoms(atpy)
      E = read_energy(atpy)
      F = read_forces(atpy)
      V = read_virial(atpy)
      EFV = ""
      (E != nothing) && (EFV *= "E")
      (F != nothing) && (EFV *= "F")
      (V != nothing) && (EFV *= "V")

      data[idx] = Dat( at,
                       config_type * ":$(length(at))$(EFV)";
                       E = E, F = F, V = V )
   end
   # verbose && toc()
   return data[1:idx]
end



function load_data(fname; kwargs...)
   if fname[end-3:end] == ".xyz"
      return read_xyz(fname; kwargs...)
   elseif fname[end-3:end] == ".jld" || fname[end-4:end] == ".jld2"
      return load(fname)
   end
   error("NBodyIPs.Data.load_data: unknown file format $(fname[end-3:end])")
end


end
