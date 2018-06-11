
"""
# `module Data`

Provides methods to read files containing simulation data. At present, only
`.xyz` is supported, via the ASE interface:
```
data = NBodyIPs.Data.read("mydata.xyz")
```
where `mydata.xyz` contains multiple configurations, will read in those
configurations, then convert them into a `JuLIP.Atoms` object, extract
energies, forces, etc store each configuration in a `Dat` and return them
as a `Vector{Dat}`.

If `d::Dat`, then one can call `energy(d), forces(d)`, etc to extract the
loaded information.
"""
module Data

using JuLIP, ASE, ProgressMeter, Base.Threads
import JuLIP: Atoms, energy, forces, virial
import Base: length

using PyCall
@pyimport ase.io as ase_io

export Dat, config_type, weight


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

if information is missing, the relevant function will return `nothing` instead
(TODO for J v1.0: change `nothing` to `missing`)
"""
struct Dat{T}
   at::Atoms
   E::Union{Void,T}         # energy
   F::Union{Void,JVecs{T}}  # forces
   S::Union{Void,JMat{T}}   # stress
   w::T   # weight
   config_type::String
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
   return ""
end


function read_xyz(fname; verbose=true,
                  exclude = [], index = ":" )
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
   @showprogress dt for n = 1:length(at_list)
      atpy = at_list[n]
      config_type = read_configtype(atpy)
      if config_type == ""
         warn("$idx has no config_type")
      end
      if config_type in exclude
         continue
      end
      E = read_energy(atpy)
      if E == nothing
         warn("$idx has not energy")
      end
      data[n] = Dat( Atoms(ASEAtoms(atpy)),
                     E,
                     read_forces(atpy),
                     read_virial(atpy),
                     1.0,
                     config_type )
   end
   verbose && toc()
   return data
end


function read(fname; kwargs...)
   if fname[end-3:end] == ".xyz"
      return read_xyz(fname; kwargs...)
   end
   error("NBodyIPs.Data.read: unknown file format $(fname[end-3:end])")
end


# --------------- FIX JLD Bug --------------


struct DatSerializer
   at
   E
   F
   S
   w
   config_type
end

import JLD
JLD.writeas(d::Dat) = DatSerializer(d.at, d.E, d.F, d.S, d.w, d.config_type)
JLD.readas(d::DatSerializer) = Dat(d.at, d.E, d.F, d.S, d.w, d.config_type)

end
