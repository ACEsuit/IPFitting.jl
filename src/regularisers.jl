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

using Sobol: SobolSeq, next!

export NBRegulariser, NBReg

struct NBRegulariser{N, T}
   N::Int
   npoints::Int
   creg::T
   r0::T
   r1::T
   sequence::Symbol
   valN::Val{N}
end

const NBReg = NBRegulariser

NBRegulariser(N, r0, r1;
             npoints = Nquad(Val(N)),
             creg = 0.1,
             sequence = :sobol) =
   NBRegulariser(N, npoints, creg, r0, r1, sequence, Val(N))


Matrix(reg::NBRegulariser{N}, basis) =
      _regularise(N, basis, reg.r0, reg.r1, reg.sequence, reg.creg, reg.npoints)


function regularise_2b(B::Vector, r0::Number, r1::Number, creg, Nquad)
   I2 = find(bodyorder.(B) .== 2)
   rr = linspace(r0, r1, Nquad)
   Φ = zeros(Nquad, length(B))
   h = (r1 - r0) / (Nquad-1)
   for (ib, b) in zip(I2, B[I2]), (iq, r) in enumerate(rr)
      Φ[iq, ib] = evaluate_dd(b, r) * sqrt(creg * h)
   end
   return Φ
end

function inv_transform(x::T, r0::T, r1::T, transform)::T where {T}
   TOL = 1e-6
   # Secant Bisection Method (transform is monotone)
   # Roots.jl is too slow (why?!?)
   t0 = transform(r0) - x
   t1 = transform(r1) - x
   r = 0.5*(r0+r1)
   while abs(r1 - r0) > TOL
      r = (r1 * t0 - r0 * t1) / (t0 - t1)
      t = transform(r) - x
      if abs(t) < 1e-7
         return r
      end
      if t*t0 > 0
         r0 = r
         t0 = t
      else
         r1 = r
         t1 = t
      end
   end
   return r
end

function cayley_menger(r::SVector{6, T}) where {T}
   A = SMatrix{5,5,T}(
      #       r12     r13     r14
      0.0,    r[1]^2, r[2]^2, r[3]^2, 1.0,
      # r12           r23     r24
      r[1]^2, 0.0,    r[4]^2, r[5]^2, 1.0,
      # r13   r23             r34
      r[2]^2, r[4]^2, 0.0,    r[6]^2, 1.0,
      # r14   r24     r34
      r[3]^2, r[5]^2, r[6]^2, 0.0,    1.0,
      1.0,    1.0,    1.0,    1.0,    0.0 )
   return det(A)
end

is_simplex(r::SVector{3}) = ((r[1] <= r[2]+r[3]+1e-4) &&
                             (r[2] <= r[3]+r[1]+1e-4) &&
                             (r[3] <= r[1]+r[2]+1e-4))

is_simplex(r::SVector{6}) = is_simplex(r[SVector(1,2,3)]) &&
                            cayley_menger(r) >= -1e-4

"""
* N : body-order
* D : dictionary
* r0, r1 : upper and lower bound on the bond-lengths
"""
function _sobol_seq(npoints, N, D, r0, r1;
                    inv_t = x -> inv_transform(x, r0, r1, D.transform),
                    verbose = false)
   # compute the boundary in transformed coordinates
   x0 = D.transform(r1)
   x1 = D.transform(r0)
   # call the inner sobol function (function barrier)
   return _sobol_inner(Val((N*(N-1)) ÷ 2), npoints, inv_t, x0, x1, verbose)
end

function _sobol_inner(::Val{DIM}, npoints, inv_t, x0::T, x1::T, verbose=false ) where {DIM, T}
   # upper and lower bounds on the transformed variables
   @assert x0 < x1
   # Sobol sequence in the [x0, x1]^dim hypercube
   s = SobolSeq(DIM, x0*ones(DIM), x1*ones(DIM))
   # temporary storage
   t = zero(MVector{DIM, T})
   # output storage
   X = SVector{DIM, T}[]
   # generate points
   failed = 0
   succesful = 0
   while length(X) < npoints
      next!(s, t)
      r = inv_t.(t)::SVector{DIM,T}
      if is_simplex(r)
         push!(X, SVector(t))
         succesful += 1
      else
         failed += 1
      end
   end
   if verbose
      @show failed, succesful
   end
   return X
end

function _cartesian_seq(npoint, N, D, r0, r1; inv_t=inv_t)
   dim = (N*(N-1))÷2
   ndim = ceil(Int, npoint^(1/dim))
   x0 = D.transform(r0)
   x1 = D.transform(r1)
   xx = linspace(x0, x1, ndim) |> collect
   oo = ones(ndim)
   if N == 2
      return xx
   end
   if N == 3
      Xmat = [kron(xx, oo, oo) kron(oo, xx, oo) kron(oo, oo, xx)]'
      X = vecs(Xmat)
      X1 = X[ [ is_simplex(inv_t.(x)) for x in X ] ]
      return X1
   end
   error("`_cartesian_seq` is only implemented for N = 2, 3")
end


Nquad(::Val{2}) = 20
Nquad(::Val{3}) = 1000
Nquad(::Val{4}) = 10_000


function _regularise(N::Integer, B::Vector{<: AbstractCalculator}, r0, r1,
                    sequence, creg, npoints, inv_t = nothing )

   # 2B is a bit simpler than the rest, treat it separately
   if N == 2
      return regularise_2b(B, r0, r1, creg, npoints)
   end

   # get the indices of the N-body basis functions
   Ib = find(bodyorder.(B) .== N)
   D = B[Ib[1]].D
   # assume all have the same dictionary
   # if not, then this regularisation is not valid
   @assert all(b.D == D  for b in B[Ib])

   if inv_t == nothing
      inv_t = x -> inv_transform(x, r0, r1, D.transform)
   end

   if sequence == :sobol
      # construct a low discrepancy sequence
      X = _sobol_seq(npoints, N, D, r0, r1; inv_t=inv_t)
   elseif sequence == :cart
      X = _cartesian_seq(npoints, N, D, r0, r1; inv_t=inv_t)
   elseif sequence isa Vector # assume it is a data vector
      X = _data_seq(npoints, N, D, r0, r1; inv_t=inv_t)
   else
      error("unknown argument `sequence = $sequence`")
   end

   # loop through sobol points and collect the laplacians at each point.
   dim = (N * (N-1)) ÷ 2
   Ψ = zeros(length(X), length(B))
   temp = zeros(length(Ib))
   B_N = [b for b in B[Ib]]
   for (ix, x) in enumerate(X)
      Ψ[ix, Ib] = creg * laplace_regulariser(x, B_N, temp, inv_t)
   end
   return Ψ
end

function laplace_regulariser(x::SVector{DIM,T}, B::Vector{TB},
                             temp::Vector{T},
                             inv_t
                             ) where {DIM, T, TB <: NBody}
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


end
