
module FastPolys

using StaticArrays

export fpoly, fpoly_d

# two auxiliary functions to make for easier assembly of the code
Base.push!(ex::Vector{Expr}, s::String) = push!(ex, parse(s))
append_str!(ex::Vector{Expr}, s::Vector{String}) = append!(ex, parse.(s))

# ACCUMULATOR FUNCTIONS
# ----------------------

# 1st order

_mon_muladd_(I::NTuple{1, Int}) =
   "m += x1_$(I[1])"

_mon_muladd_d_(I::NTuple{1, Int}) =
   ["dm_$(I[1]) += dx1_$(I[1])",]

# 2nd order

_mon_muladd_(I::NTuple{2, Int}) =
   "m = muladd(x1_$(I[1]), x2_$(I[2]), m)"

_mon_muladd_d_(I::NTuple{2, Int}) =
   ["dm_$(I[1]) = muladd(dx1_$(I[1]),  x2_$(I[2]), dm_$(I[1]))",
    "dm_$(I[2]) = muladd( x1_$(I[1]), dx2_$(I[2]), dm_$(I[2]))"]


# 3rd order

_mon_muladd_(I::NTuple{3, Int}) =
   "m = muladd(x1_$(I[1]), x2_$(I[2]) * x3_$(I[3]), m)"

_mon_muladd_d_(I::NTuple{3, Int}) =
   ["dm_$(I[1]) = muladd(dx1_$(I[1]),  x2_$(I[2]) *  x3_$(I[3]), dm_$(I[1]))",
    "dm_$(I[2]) = muladd( x1_$(I[1]), dx2_$(I[2]) *  x3_$(I[3]), dm_$(I[2]))",
    "dm_$(I[3]) = muladd( x1_$(I[1]),  x2_$(I[2]) * dx3_$(I[3]), dm_$(I[3]))"]


# 4th order

_mon_muladd_(I::NTuple{4, Int}) =
   "m = muladd(x1_$(I[1]) * x2_$(I[2]), x3_$(I[3]) * x4_$(I[4]), m)"

_mon_muladd_d_(I::NTuple{4, Int}) =
   ["a12 = x1_$(I[1]) * x2_$(I[2])",
    "a34 = x3_$(I[3]) * x4_$(I[4])",
    "dm_$(I[1]) = muladd(dx1_$(I[1]) *  x2_$(I[2]), a34, dm_$(I[1]))",
    "dm_$(I[2]) = muladd( x1_$(I[1]) * dx2_$(I[2]), a34, dm_$(I[2]))",
    "dm_$(I[3]) = muladd( a12, dx3_$(I[3]) *  x4_$(I[4]), dm_$(I[3]))",
    "dm_$(I[4]) = muladd( a12,  x3_$(I[3]) * dx4_$(I[4]), dm_$(I[4]))"]


# 5th order

_mon_muladd_(I::NTuple{5, Int}) =
   "m = muladd(x1_$(I[1]) * x2_$(I[2]), x3_$(I[3]) * x4_$(I[4]) * x5_$(I[5]), m)"

_mon_muladd_d_(I::NTuple{5, Int}) =
   ["a12 = x1_$(I[1]) * x2_$(I[2])",
    "a123 = a12 * x3_$(I[3])",
    "a45 = x4_$(I[4]) * x5_$(I[5])",
    "a345 = x3_$(I[3]) * a45",
    "dm_$(I[1]) = muladd(dx1_$(I[1]) *  x2_$(I[2]), a345, dm_$(I[1]))",
    "dm_$(I[2]) = muladd( x1_$(I[1]) * dx2_$(I[2]), a345, dm_$(I[2]))",
    "dm_$(I[3]) = muladd( a12, dx3_$(I[3]) *  a45, dm_$(I[3]))",
    "dm_$(I[4]) = muladd( a123,  dx4_$(I[4]) *  x5_$(I[5]), dm_$(I[4]))",
    "dm_$(I[5]) = muladd( a123,   x4_$(I[4]) * dx5_$(I[5]), dm_$(I[5]))" ]


function _check_input_(D, A)
   @assert length(A) == D
   @assert typeof(A[1]) <: NTuple
   @assert eltype(A[1]) <: Integer
   for n = 2:length(A)
      @assert typeof(A[1]) == typeof(A[n])
   end
end

"""
`function fpoly(X, ::Val{A})`

TODO: write documentation
"""
@generated function fpoly(x::NTuple{D, SVector{N, T}},
                          ::Val{A}) where {D, N, T, A}
   # make some assumptions about what the input is
   _check_input_(D, A)

   # generate the code
   # -----------------
   exprs = Expr[]
   # read all the data from the StaticArrays
   #     for very small monomials, this could be the bottleneck, but for
   #     reasonably large ones it should make no difference
   for d = 1:D, n = 1:N
      push!(exprs, "x$(d)_$(n) = x[$(d)][$(n)]")
   end
   # also initialise the output
   push!(exprs, "m = zero($T)")
   # evaluate the monomial
   for I in zip(A...)
      push!(exprs, _mon_muladd_(I))
   end
   quote
      $(Expr(:meta, :inline))
      @inbounds $(Expr(:block, exprs...))
      return m
   end
end


"""
`function fpoly(X, dX, ::Val{A})`

TODO: write documentation
"""
@generated function fpoly_d(
         x::NTuple{D, SVector{N, T}},
         dx::NTuple{D, SVector{N, T}},
         ::Val{A}
      ) where {D, N, T, A}
   # make some assumptions about what the input is
   _check_input_(D, A)

   # generate the code
   # -----------------
   exprs = Expr[]
   # read all the data from the StaticArrays
   #     for very small monomials, this could be the bottleneck, but for
   #     reasonably large ones it should make no difference
   for d = 1:D, n = 1:N
      push!(exprs, "x$(d)_$(n) = x[$(d)][$(n)]")
      push!(exprs, "dx$(d)_$(n) = dx[$(d)][$(n)]")
   end
   # initialise the output
   push!(exprs, "m = zero($T)")
   for n = 1:N
      push!(exprs, "dm_$(n) = zero($T)")
   end
   # evaluate the monomial and its derivative
   for I in zip(A...)
      append_str!(exprs, _mon_muladd_d_(I))
   end
   # collect the dm variables into an SVector
   coll = "SVector(" * prod("dm_$n, " for n = 1:N) * ")"
   quote
      $(Expr(:meta, :inline))
      @inbounds $(Expr(:block, exprs...))
      $(parse(coll))
   end
end


end
