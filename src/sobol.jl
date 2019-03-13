using Sobol: SobolSeq, next!

"""
* N : body-order
* D : dictionary
* r0, r1 : upper and lower bound on the bond-lengths
"""
function filtered_sobol(x0, x1, filter; npoints=nothing,
                                        nfiltered=nothing)
   @assert npoints isa Integer && nfiltered isa Integer
   x0, x1 = min.(x0, x1), max.(x0, x1) # @assert all(x0 .< x1)
   return _sobol_inner(SVector(x0...), SVector(x1...), npoints, nfiltered, filter)
end

function _sobol_inner(x0::SVector{DIM, T}, x1::SVector{DIM, T},
                      npoints, nfiltered, filter ) where {DIM, T}
   # Sobol sequence in the [x0, x1]^dim hypercube
   s = SobolSeq(DIM, x0, x1)
   # temporary storage
   t = zero(MVector{DIM, T})
   # output storage
   X = SVector{DIM, T}[]
   # generate points
   npts = 0
   nfilt = 0
   while (npts < npoints) && (nfilt < nfiltered)
      next!(s, t)
      npts += 1
      if filter(t)
         push!(X, SVector(t))
         nfilt += 1
      end
   end
   return X
end


# function _cartesian_seq(npoint, N, D, r0, r1; inv_t=inv_t)
#    dim = (N*(N-1))÷2
#    ndim = ceil(Int, npoint^(1/dim))
#    x0 = D.transform(r0)
#    x1 = D.transform(r1)
#    xx = linspace(x0, x1, ndim) |> collect
#    oo = ones(ndim)
#    if N == 2
#       return xx
#    end
#    if N == 3
#       Xmat = [kron(xx, oo, oo) kron(oo, xx, oo) kron(oo, oo, xx)]'
#       X = vecs(Xmat)
#       X1 = X[ [ is_simplex(inv_t.(x)) for x in X ] ]
#       return X1
#    end
#    error("`_cartesian_seq` is only implemented for N = 2, 3")
# end



# ------------------------------------------------------------------
#               FILTERS
# ------------------------------------------------------------------

function cayley_menger(r::SVector{6, T}) where {T}
   A = @SMatrix( T[
      #       r12     r13      r14
      0.0     r[1]^2  r[2]^2   r[3]^2   1.0;
      # r12           r23      r24
      r[1]^2  0.0     r[4]^2   r[5]^2   1.0;
      # r13    r23             r34
      r[2]^2  r[4]^2  0.0      r[6]^2   1.0;
      # r14    r24    r34
      r[3]^2  r[5]^2  r[6]^2   0.0      1.0;
      1.0     1.0     1.0      1.0      0.0 ])
   return det(A)
end

function cayley_menger(r::SVector{10, T}) where {T}
   A = @SMatrix( T[
      0.0     1.0     1.0     1.0      1.0      1.0;
      #                r12     r13      r14     r15
      1.0     0.0     r[1]^2  r[2]^2   r[3]^2   r[4]^2;
      #       r12             r23      r24      r25
      1.0     r[1]^2  0.0     r[5]^2   r[6]^2   r[7]^2;
      #       r13     r23              r34      r35
      1.0     r[2]^2  r[5]^2  0.0      r[8]^2   r[9]^2;
      #       r14     r24     r34               r45
      1.0     r[3]^2  r[6]^2  r[6]^2   0.0      r[10]^2;
      #       r15     r25     r35      r45
      1.0     r[4]^2  r[7]^2  r[6]^2   r[10]^2  1.0 ])
   return det(A)
end


# bond-length filters
bl_is_simplex(r::MVector) = bl_is_simplex(SVector(r))

bl_is_simplex(r::SVector{1}) = true

bl_is_simplex(r::SVector{3}) = (abs(r[2]-r[1]) <= r[3] <= r[2]+r[1])

bl_is_simplex(r::SVector{6}) =
   bl_is_simplex(SVector(r[1], r[2], r[4])) &&   # r12, r13, r23
   bl_is_simplex(SVector(r[1], r[3], r[5])) &&   # r12, r14, r2
   bl_is_simplex(SVector(r[2], r[3], r[6])) &&   # r13, r14, r34
   bl_is_simplex(SVector(r[4], r[5], r[6])) &&    # r23, r24, r34
   (cayley_menger(r) >= 0)

# r3^2 = r1^2 + r2^2 - 2 ψ r1 r2
# cosangle3(r1, r2, r3) = (r1^2 + r2^2 - r3^2) / (2 * r1 * r2)

ba_is_simplex(r::SVector{1}) = true

ba_is_simplex(rθ::SVector{3}) = true
ba_is_simplex(r::SVector{2}, θ::SVector{1}) = true

ba_is_simplex(rψ::SVector{6}) = ba_is_simplex(SVector(rψ[1], rψ[2], rψ[3]),
                                              SVector(rψ[4], rψ[5], rψ[6]) )

ba_is_simplex(r::SVector{3}, ψ::SVector{3}) =
   bl_is_simplex(vcat(SVector(1.0,1.0,1.0), sqrt.(2 - 2*ψ)))


"""
an auxiliary function to help transform between real and transformed
   coordinates where the Sobol sequence is constructed ...
"""
function inv_transform(x::T, r0::T, r1::T, D)::T where {T}
   TOL = 1e-6
   # Secant Bisection Method (transform is monotone)
   # Roots.jl is too slow (why?!?)
   t0 = transform(D, r0) - x
   t1 = transform(D, r1) - x
   r = 0.5*(r0+r1)
   while abs(r1 - r0) > TOL
      r = (r1 * t0 - r0 * t1) / (t0 - t1)
      t = transform(D, r) - x
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



# using BenchmarkTools, StaticArrays
# r = @SVector rand(6)
# cayley_menger(r)
# @btime cayley_menger($r)
# @btime bl_is_simplex($r)
# @btime bl_is_simplex_2($r)
#
# # cnt = 0
# # for n = 1:200
# #    r = @SVector rand(6)
# #    if (bl_is_simplex(r) && cayley_menger(r) < -0.001) ||
# #             (!bl_is_simplex(r) && cayley_menger(r) > 0.001)
# #       r = round.(r, 4)
# #       cnt += 1
# #       @show r
# #       @show cayley_menger(r), bl_is_simplex(r)
# #    end
# # end
#
# ##
# tic()
# for n = 1:100_000
#    # r = vcat((@SVector rand(3)), 2 * ((@SVector rand(3)) - 0.5))
#    r = 1.0 + (@SVector rand(6))
#    if bl_is_simplex(r) != bl_is_simplex_2(r)
#       r = round.(r, 4)
#       @show r
#       @show cayley_menger(r)
#       @show bl_is_simplex(r), bl_is_simplex_2(r)
#    end
# end
# toc()
