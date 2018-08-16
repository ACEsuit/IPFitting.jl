
module DataTypes

using StaticArrays
using JuLIP: vecs, mat
import NBodyIPFitting: vec, devec

export ENERGY, FORCES, VIRIAL

const _IV = [1,2,3,5,6,9]
const _IVst = SVector(1,2,3,5,6,9)


const ENERGY = "E"
vec(::Val{:E}, E::Real) = [E]
devec(::Val{:E}, x::AbstractVector) = ((@assert length(x) == 1); x[1])

const FORCES = "F"
vec(::Val{:F}, F::Vector{JVec{T}}) where {T} = vec(::Val{F}, mat(F))
vec(::Val{:F}, F::Matrix) = F[:]
devec(::Val{:F}, x::AbstractVector) = vecs(reshape(x, 3, :))

const VIRIAL = "V"
vec(::Val{:V}, v::AbstractVector) = (@assert length(v) == 6; Vector(v))
vec(::Val{:V}, V::AbstractMatrix) = (@assert size(V) == (3,3); V[_IV])
devec(::Val{V}, x::AbstractVector{T}) where {T <: AbstractFloat} =
   SMatrix(x[1], x[2], x[3], x[2], x[4], x[5], x[3], x[5], x[6])
   #  1  2  3
   #  2  4  5
   #  3  5  6


end
