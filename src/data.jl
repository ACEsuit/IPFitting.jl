

module Data

using JuLIP, ASE
import JuLIP: energy, forces
import Base: read, length

using PyCall
@pyimport ase.io as ase_io


"""
`Dat`: store one simulation data point
"""
struct Dat{T}
   at::Atoms
   E::T
   F::JVecs{T}
end

energy(d::Dat) = d.E
forces(d::Dat) = d.F
length(d::Dat) = length(d.at)


function read_xyz(fname; index = ":", verbose=true)
   if verbose
      println("Reading in $fname ...")
   end
   at_list = ase_io.read(fname, index=index)
   data = Dat{Float64}[]
   for atpy in at_list
      E = atpy[:get_potential_energy]()
      F = atpy[:get_array]("force")' |> vecs
      at = Atoms(ASEAtoms(atpy))
      push!(data, Dat(at, E, F))
   end
   return data
end



function Base.read(fname; kwargs...)
   if fname[end-3:end] == ".xyz"
      return read_xyz(fname; kwargs...)
   end
   error("NBodyIPs.Data.read: unknown file format $(fname[end-3:end])")
end

end
