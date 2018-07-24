
"""
# `module Data`

Provides methods to read files containing simulation data. Primarily this is
intended to load `.xyz` files and convert them to JuLIP-compatible data:
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

using JuLIP, ASE, ProgressMeter, FileIO
import JuLIP: Atoms, energy, forces, virial
import Base: length

using PyCall
@pyimport ase.io as ase_io

export config_type, weight, load_data


"""
`Dat`: store one simulation data point. If `d::Dat`, to obtain the data, use
```
Atoms(d)
energy(d)
forces(d)
virial(d)
weight(d)
config_type(d)
length(d)    # number of atoms
```
If information is missing, the relevant function will return `nothing` instead
(TODO for J v1.0: change `nothing` to `missing`)
"""
mutable struct Dat{T}
   at::Atoms
   E::Union{Void, T}         # energy
   F::Union{Void, JVecs{T}}  # forces
   S::Union{Void, JMat{T}}   # stress
   w                         # weight, could be anything?
   config_type::Union{Void, String}
   D::Dict{String, Any}
end

Atoms(d) = d.at
energy(d::Dat) = d.E
forces(d::Dat) = d.F
virial(d::Dat) = d.S
weight(d::Dat) = d.w
length(d::Dat) = length(d.at)
config_type(d::Dat) = d.config_type

function read_energy(atpy)
   for key in keys(atpy[:info])
      if lowercase(key) == "dft_energy"
         return atpy[:info][key]
      end
   end
   try
      return atpy[:get_potential_energy]()
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
   data = Vector{Dat{Float64}}(length(at_list))
   idx = 0
   if verbose
      println("Processing data ...")
      tic()
   end
   dt = verbose ? 1.0 : 0.0
   @showprogress dt for atpy in at_list

      # get the config type and decide whether to keep or skip this config
      config_type = read_configtype(atpy)
      if config_type == nothing
         warn("$idx has no config_type")
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
      data[idx] = Dat( Atoms(ASEAtoms(atpy)),
                       read_energy(atpy),
                       read_forces(atpy),
                       read_virial(atpy),
                       nothing,
                       config_type,
                       Dict{String,Any}() )
   end
   verbose && toc()
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
