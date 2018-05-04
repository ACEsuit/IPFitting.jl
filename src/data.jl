
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

using JuLIP, ASE, ProgressMeter
import JuLIP: Atoms, energy, forces
import Base: length

using PyCall
@pyimport ase.io as ase_io


"""
`Dat`: store one simulation data point. If `d::Dat`, to obtain the data, use
```
Atoms(d)
energy(d)
forces(d)
length(d)    # number of atoms
```

## TODO
```
stress(d)
```
"""
struct Dat{T}
   at::Atoms
   E::T
   F::JVecs{T}
end

Atoms(d) = d.at 
energy(d::Dat) = d.E
forces(d::Dat) = d.F
length(d::Dat) = length(d.at)


function read_xyz(fname; index = ":", verbose=true,
                  dt = verbose ? 0.5 : Inf )
   if verbose
      println("Reading in $fname ...")
   end
   at_list = ase_io.read(fname, index=index)
   data = Dat{Float64}[]
   @showprogress dt "Processing ..." for atpy in at_list
      E = atpy[:get_potential_energy]()
      F = atpy[:get_array]("force")' |> vecs
      at = Atoms(ASEAtoms(atpy))
      push!(data, Dat(at, E, F))
   end
   return data
end



function read(fname; kwargs...)
   if fname[end-3:end] == ".xyz"
      return read_xyz(fname; kwargs...)
   end
   error("NBodyIPs.Data.read: unknown file format $(fname[end-3:end])")
end

end
