"""

`module Regularisers`

## Old Documentation Copy-Pasted ...

`regularise_2b(B, r0, r1; creg = 1e-2, Nquad = 20)`

construct a regularising stabilising matrix that acts only on 2-body
terms.

* `B` : basis
* `r0, r1` : upper and lower bound over which to integrate
* `creg` : multiplier
* `Nquad` : number of quadrature / sample points

```
   I = ∑_r h | ∑_j c_j ϕ_j''(r) |^2
     = ∑_{i,j} c_i c_j  ∑_r h ϕ_i''(r) ϕ_j''(r)
     = c' * M * c
   M_ij = ∑_r h ϕ_i''(r) ϕ_j''(r) = h * Φ'' * Φ''
   Φ''_ri = ϕ_i''(r)
```


`regularise(N, B, r0, r1; kwargs...)`

* `N::Integer` : body order
* `B::Vector` : basis
"""
module Regularisers

using StaticArrays
using JuLIP: AbstractCalculator
using JuLIP.Potentials: evaluate_d
using NBodyIPs: bodyorder, transform, evaluate_many_ricoords!
using NBodyIPFitting.Tools: @def

import Base: Matrix

export BLRegulariser, BLReg, BARegulariser, BAReg

abstract type AbstractRegulariser end

abstract type NBodyRegulariser{N} end

@def nbregfields begin
   N::Int
   npoints::Int
   creg::T
   r0::T
   r1::T
   sequence::Symbol
   freg::Function
   valN::Val{N}
end

struct BLRegulariser{N, T} <: NBodyRegulariser{N}
   @nbregfields
end

struct BARegulariser{N, T} <: NBodyRegulariser{N}
   @nbregfields
end


const BLReg = BLRegulariser
const BAReg = BARegulariser

BLRegulariser(N, r0, r1;
             npoints = Nquad(Val(N)),
             creg = 0.1,
             sequence = :sobol,
             freg = laplace_regulariser) =
   BLRegulariser(N, npoints, creg, r0, r1, sequence, freg, Val(N))

BARegulariser(N, r0, r1;
             npoints = Nquad(Val(N)),
             creg = 0.1,
             sequence = :sobol,
             freg = laplace_regulariser) =
   BARegulariser(N, npoints, creg, r0, r1, sequence, freg, Val(N))

# ---------------------------------------------------------------------------

include("sobol.jl")

# ---------------------------------------------------------------------------


Nquad(::Val{2}) = 100
Nquad(::Val{3}) = 1000
Nquad(::Val{4}) = 10_000

_bainvt(inv_t, x::StaticVector{1}) =
      inv_t.(x), SVector()
_bainvt(inv_t, x::StaticVector{3}) =
      SVector(inv_t(x[1]), inv_t(x[2])), SVector(x[3])
_bainvt(inv_t, x::StaticVector{6}) =
      SVector(inv_t(x[1]), inv_t(x[2]), inv_t(x[3])),  SVector(x[4], x[5], x[6])

#
# this converts the Regulariser type information to a matrix that can be
# attached to the LSQ problem .
#
function Matrix(reg::NBodyRegulariser{N}, B::Vector{<: AbstractCalculator}
                ) where {N}

   # 2B is a bit simpler than the rest, treat it separately
   if N == 2
      return regularise_2b(B, reg.r0, reg.r1, reg.creg, reg.npoints)
   end

   # TODO: this assumes that all elements of B have the same descriptor
   # get the indices of the N-body basis functions
   Ib = find(bodyorder.(B) .== N)
   D = B[Ib[1]].D
   # assume all have the same dictionary
   # if not, then this regularisation is not valid
   @assert all(b.D == D  for b in B[Ib])

   # scalar inverse transform
   inv_t = x -> inv_transform(x, reg.r0, reg.r1, D)

   # filter
   if reg isa BLRegulariser
      # vectorial inverse transform
      inv_tv = x -> inv_t.(x)
      filter = x -> bl_is_simplex( inv_tv(x) )
      x0 = transform(D, reg.r0) * SVector(ones((N*(N-1))÷2)...)
      x1 = transform(D, reg.r1) * SVector(ones((N*(N-1))÷2)...)
   elseif reg isa BARegulariser
      inv_tv = x -> _bainvt(inv_t, x)
      filter = x -> ba_is_simplex( inv_tv(x)... )
      x0 = vcat( transform(D, reg.r0) * SVector(ones(N-1)...),
                 - SVector(ones( ((N-1)*(N-2))÷2 )...) )
      x1 = vcat( transform(D, reg.r1) * SVector(ones(N-1)...),
                 SVector(ones( ((N-1)*(N-2))÷2 )...) )
   else
      error("Unknown type of reg: `typeof(reg) == $(typeof(reg))`")
   end

   if reg.sequence == :sobol
      # construct a low discrepancy sequence
      X = filtered_sobol(x0, x1, filter; npoints = 10*reg.npoints, nfiltered=reg.npoints)
   elseif reg.sequence == :cart
      error("TODO: implement `sequence == :cart`")
   elseif reg.sequence isa Vector # assume it is a data vector
      # X = _data_seq(npoints, N, D, r0, r1; inv_t=inv_t)
      error("TODO: implement `sequence isa Vector`")
   else
      error("unknown argument `sequence = $sequence`")
   end

   # loop through sobol points and collect the laplacians at each point.
   return reg.creg * assemble_reg_matrix(X, [b for b in B[Ib]], length(B), Ib,
                                         inv_tv, reg.freg)
end



function regularise_2b(B::Vector, r0::Number, r1::Number, creg, Nquad)
   I2 = find(bodyorder.(B) .== 2)
   rr = linspace(r0, r1, Nquad)
   Φ = zeros(Nquad, length(B))
   h = (r1 - r0) / (Nquad-1)
   for (ib, b) in zip(I2, B[I2]), (iq, r) in enumerate(rr)
      Φ[iq, ib] = (evaluate_d(b, r+1e-2)-evaluate_d(b, r-1e-2))/(2e-2) * sqrt(creg * h)
   end
   return creg * Φ
end



function assemble_reg_matrix(X, B, nB, Ib, inv_tv, freg)
   Ψ = zeros(length(X), nB)
   temp = zeros(length(Ib))
   for (ix, x) in enumerate(X)
      Ψ[ix, Ib] = freg(x, B, temp, inv_tv)
   end
   return Ψ
end


function laplace_regulariser(x::SVector{DIM,T}, B::Vector{<: AbstractCalculator},
                             temp::Vector{T}, inv_tv) where {DIM, T}
   if !isleaftype(eltype(B))
      warn("laplace_regulariser: `TB` is not a leaf type")
   end
   h = 1e-2
   r = inv_tv(x)
   evaluate_many_ricoords!(temp, B, r)
   L = 2 * DIM * copy(temp)

   EE = @SMatrix eye(DIM)
   for j = 1:DIM
      rp = inv_tv(x + h * EE[:,j])
      evaluate_many_ricoords!(temp, B, rp)
      L -= temp
      rm = inv_tv(x - h * EE[:,j])
      evaluate_many_ricoords!(temp, B, rm)
      L -= temp
   end

   return L/h^2
end


function energy_regulariser(x::SVector{DIM,T}, B::Vector{<: AbstractCalculator},
                            temp::Vector{T}, inv_tv) where {DIM, T}
   evaluate_many_ricoords!(temp, B, r)
   return - temp
end

end
