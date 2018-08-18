
module DataTypes

using StaticArrays
using JuLIP: vecs, mat, JVec,
import NBodyIPFitting: vec, devec

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

const FORCES = "F"
const ValF = Val{:F}
vec(v::ValF, F::Vector{JVec{T}}) where {T} = vec(v, mat(F))
vec(::ValF, F::Matrix) = F[:]
devec(::ValF, x::AbstractVector) = vecs(reshape(x, 3, :))

const VIRIAL = "V"
const ValV = Val{
vec(::ValV, v::AbstractVector) = (@assert length(v) == 6; Vector(v))
vec(::ValV, V::AbstractMatrix) = (@assert size(V) == (3,3); V[_IV])
devec(::ValV, x::AbstractVector{T}) where {T <: AbstractFloat} =
   SMatrix(x[1], x[6], x[5], x[6], x[2], x[4], x[5], x[4], x[3])
   #  1  6  5
   #  6  2  4
   #  5  4  3


end
