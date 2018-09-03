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
using NBodyIPs: bodyorder

import Base: Matrix

export BLRegulariser, BLReg, AbstractRegulariser

abstract type AbstractRegulariser{N} end

struct BLRegulariser{N, T} <: AbstractRegulariser{N}
   N::Int
   npoints::Int
   creg::T
   r0::T
   r1::T
   sequence::Symbol
   valN::Val{N}
end

struct BARegulariser{N, T} <: AbstractRegulariser{N}
   N::Int
   npoints::Int
   creg::T
   r0::T
   r1::T
   sequence::Symbol
   valN::Val{N}
end

const BLReg = BLRegulariser
const BAReg = BARegulariser

BLRegulariser(N, r0, r1;
             npoints = Nquad(Val(N)),
             creg = 0.1,
             sequence = :sobol) =
   BLRegulariser(N, npoints, creg, r0, r1, sequence, Val(N))

BARegulariser(N, r0, r1;
             npoints = Nquad(Val(N)),
             creg = 0.1,
             sequence = :sobol) =
   BARegulariser(N, npoints, creg, r0, r1, sequence, Val(N))

# ---------------------------------------------------------------------------

include("sobol.jl")

# ---------------------------------------------------------------------------

Matrix(reg::BLRegulariser, basis) = _regularise(reg, basis)


function regularise_2b(B::Vector, r0::Number, r1::Number, creg, Nquad)
   I2 = find(bodyorder.(B) .== 2)
   rr = linspace(r0, r1, Nquad)
   Φ = zeros(Nquad, length(B))
   h = (r1 - r0) / (Nquad-1)
   for (ib, b) in zip(I2, B[I2]), (iq, r) in enumerate(rr)
      Φ[iq, ib] = (evaluate_d(b, r+1e-2)-evaluate_d(b, r-1e-2))/(2e-2) * sqrt(creg * h)
   end
   return Φ
end


Nquad(::Val{2}) = 20
Nquad(::Val{3}) = 1000
Nquad(::Val{4}) = 10_000

#
# this converts the Regulariser type information to a matrix that can be
# attached to the LSQ problem .
#
function Matrix(reg::AbstractRegulariser{N}, B::Vector{<: AbstractCalculator})

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

   # inverse transform
   inv_t = x -> inv_transform(x, reg.r0, reg.r1, D.transform)

   # filter
   if reg isa BLRegulariser
      filter = x -> bl_is_simplex( inv_t.(x) )
      x0 = D.transform(reg.r0) * SVector(ones((N*(N-1))÷2)...)
      x1 = D.transform(reg.r1) * SVector(ones((N*(N-1))÷2)...)
   else if reg isa BARegulariser
      filter = x -> ba_is_simplex( inv_t.(SVector(x[1], x[2], x[3])),
                                   SVector(x[4], x[5], x[6]) )
      x0 = vcat( D.transform(reg.r0) * SVector(ones(N)...),
                 - SVector(ones( ((N-1)*(N-2))÷2 )...) )
      x1 = vcat( D.transform(reg.r1) * SVector(ones(N)...),
                 SVector(ones( ((N-1)*(N-2))÷2 )...) )
   else
      error("Unknown type of reg: `typeof(reg) == $(typeof(reg))`")
   end

   if reg.sequence == :sobol
      # construct a low discrepancy sequence
      X = filtered_sobol(x0, x1, filter; npoints = 100*reg.npoints, nfiltered=reg.npoints)
   elseif reg.sequence == :cart
      error("TODO: implement `sequence == :cart`")
   elseif reg.sequence isa Vector # assume it is a data vector
      # X = _data_seq(npoints, N, D, r0, r1; inv_t=inv_t)
      error("TODO: implement `sequence isa Vector`")
   else
      error("unknown argument `sequence = $sequence`")
   end

   # loop through sobol points and collect the laplacians at each point.
   return reg.creg * assemble_reg_matrix(X, [b for b in B[Ib]], nB, Ib, inv_t)
end


function assemble_reg_matrix(X, B, nB, Ib, inv_t)
   Ψ = zeros(length(X), length(B))
   temp = zeros(length(Ib))
   for (ix, x) in enumerate(X)
      Ψ[ix, Ib] = laplace_regulariser(x, B, temp, inv_t)
   end
   return Ψ
end


function laplace_regulariser(x::SVector{DIM,T}, B::Vector{TB},
                             temp::Vector{T},
                             inv_t
                             ) where {DIM, T, TB <: AbstractCalculator}
   if !isleaftype(TB)
      warn("laplace_regulariser: `TB` is not a leaf type")
   end
   h = 1e-2
   r = inv_t.(x)
   evaluate_many!(temp, B, r)
   L = 2 * DIM * copy(temp)

   EE = @SMatrix eye(DIM)
   for j = 1:DIM
      rp = inv_t.(x + h * EE[:,j])
      evaluate_many!(temp, B, rp)
      L -= temp
      rm = inv_t.(x - h * EE[:,j])
      evaluate_many!(temp, B, rm)
      L -= temp
   end

   return L/h^2
end
