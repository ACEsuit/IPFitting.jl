
module DataTypes

using StaticArrays
using JuLIP: vecs, mat, JVec, energy, forces, virial
import NBodyIPFitting: vec, devec, evaluate_lsq

export ENERGY, FORCES, VIRIAL

# using Voigt convention for vectorising  symmetric 3 x 3 matrix
#  1  6  5
#  6  2  4
#  5  4  3

const _IV = [1,5,9,6,3,2]
const _IVst = SVector(1,5,9,6,3,2)


const ENERGY = "E"
const ValE = Val{:E}
vec(::ValE, E::Real) = [E]
devec(::ValE, x::AbstractVector) = ((@assert length(x) == 1); x[1])
evaluate_lsq(::ValE, B, at) = energy(B, at)

const FORCES = "F"
const ValF = Val{:F}
vec(v::ValF, F::Vector{<:JVec}) = vec(v, mat(F))
vec(::ValF, F::Matrix) = F[:]
devec(::ValF, x::AbstractVector) = vecs(reshape(x, 3, :))
evaluate_lsq(::ValF, B, at) = forces(B, at)

const VIRIAL = "V"
const ValV = Val{:V}
vec(::ValV, v::AbstractVector) = (@assert length(v) == 6; Vector(v))
vec(::ValV, V::AbstractMatrix) = (@assert size(V) == (3,3); V[_IV])
devec(::ValV, x::AbstractVector) =
   SMatrix{3,3}(x[1], x[6], x[5], x[6], x[2], x[4], x[5], x[4], x[3])
   #  1  6  5
   #  6  2  4
   #  5  4  3
evaluate_lsq(::ValV, B, at) = virial(B, at)


end
