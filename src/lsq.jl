
module Lsq

using StaticArrays 
using JuLIP: AbstractCalculator
using NBodyIPFitting: Dat

# components of the stress (up to symmetry)
const _IS = SVector(1,2,3,5,6,9)


"""
Take a basis and split it into individual basis groups.
"""
function split_basis(basis)
   # get the types of the individual basis elements
   tps = typeof.(basis)
   Iord = Vector{Int}[]
   Bord = Any[]
   for tp in unique(tps)
      # find which elements of basis have type `tp`
      I = find( tp .== tps )
      push!(Iord, I)
      push!(Bord, [b for b in basis[I]])
   end
   return Bord, Iord
end


# fill the LSQ system, i.e. evaluate basis at data points
function evallsq(d::Dat, B::AbstractVector{TB}
                 ) where {TB <: AbstractCalculator}
   if !(isleaftype(T))
      return evallsq_split(d, B)
   end

   # -----------------------------------------------
   # TB is a leaf-type so we can use "evaluate_many"
   # -----------------------------------------------
   len = length(d)
   at = Atoms(d)

   D = Dict{String, Any}()
   if energy(d) != nothing
      D["E"] = energy(B, at)
   end
   if (forces(d) != nothing) && len > 1
      D["F"] = forces(B, at)
   end
   if (virial(d) != nothing)
      Vs = virial(Bord[n], at)
      # store a vector rather than a matrix
      D["V"] = [ v[_IS] for v in Vs ]
   end
   # TODO: allow generic types of data

   return D
end

function _cat_(As, Iord)
   TA = eltype(As[1])
   A = zeros(TA, sum(size(AA)[end] for AA in As))
   for i = 1:length(As)
      A[Iord[i]] = As[i]
   end
   return A
end

"""
split the Basis `B` into subsets of identical types and evaluate
those independently (fast).
"""
function evallsq_split(d, B)
   # TB is not a leaf-type so we should split the basis to be able to
   # evaluate_many & co
   Bord, Iord = split_basis(basis)
   D_ord = [ evallsq(d, BB)  for BB in Bord ]
   # each D_ord[i] is a Dict storing the LSQ system components.
   # we assume that all D_ord[i] contain the same keys.
   D = Dict{String,Any}()
   for key in keys(D_ord[1])
      D[key] = _cat_([DD[key] for DD in DD_ord], Iord)
   end
   return D
end



end
