using BenchmarkTools, StaticArrays

Base.push!(ex::Vector{Expr}, s::String) = push!(ex, parse(s))
append_str!(ex::Vector{Expr}, s::Vector{String}) = append!(ex, parse.(s))

# ACCUMULATOR FUNCTIONS
# ----------------------

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


@generated function fmon(x::NTuple{D, SVector{N, T}},
                         ::Val{A}) where {D, N, T, A}
   # make some assumptions about what the input is
   @assert length(A) == D
   @assert typeof(A[1]) <: NTuple
   @assert eltype(A[1]) <: Integer
   for n = 2:length(A)
      @assert typeof(A[1]) == typeof(A[n])
   end

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
      @inbounds $(Expr(:block, exprs...))
      return m
   end
end


@generated function fmon_d(
         x::NTuple{D, SVector{N, T}},
         dx::NTuple{D, SVector{N, T}},
         ::Val{A}
      ) where {D, N, T, A}
   # make some assumptions about what the input is
   @assert length(A) == D
   @assert typeof(A[1]) <: NTuple
   @assert eltype(A[1]) <: Integer
   for n = 2:length(A)
      @assert typeof(A[1]) == typeof(A[n])
   end

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
      @inbounds $(Expr(:block, exprs...))
      $(parse(coll))
   end
end



@generated function gmon(
               x1::SVector{N, T}, a1::Val{A1},
               x2::SVector{N, T}, a2::Val{A2},
               x3::SVector{N, T}, a3::Val{A3}
            ) where {N, T, A1, A2, A3}
   @assert typeof(A1) == typeof(A2) == typeof(A3) <: NTuple
   exprs = Expr[]
   for n = 1:N
      push!(exprs, parse("x1_$n = x1[$n]"))
      push!(exprs, parse("x2_$n = x2[$n]"))
      push!(exprs, parse("x3_$n = x3[$n]"))
   end
   for (i1, i2, i3) in zip(A1, A2, A3)
      push!(exprs, parse("m = muladd(x1_$i1 * x2_$i2, x3_$i3, m)"))
   end
   quote
      m = zero(eltype(x1))
      @inbounds $(Expr(:block, exprs...))
      return m
   end
end


@generated function gmon(
               x1::SVector{N, T}, a1::Val{A1},
               x2::SVector{N, T}, a2::Val{A2},
               x3::SVector{N, T}, a3::Val{A3},
               x4::SVector{N, T}, a4::Val{A4}
            ) where {N, T, A1, A2, A3, A4}
   exprs = Expr[]
   for n = 1:N
      push!(exprs, parse("x1_$n = x1[$n]"))
      push!(exprs, parse("x2_$n = x2[$n]"))
      push!(exprs, parse("x3_$n = x3[$n]"))
      push!(exprs, parse("x4_$n = x4[$n]"))
   end
   for (i1, i2, i3, i4) in zip(A1, A2, A3, A4)
      push!(exprs, parse("m = muladd(x1_$i1 * x2_$i2, x3_$i3 * x4_$i4, m)"))
   end
   quote
      m = zero(eltype(x1))
      @inbounds $(Expr(:block, exprs...))
      return m
   end
end


@generated function gmon(
               x1::SVector{N, T}, a1::Val{A1},
               x2::SVector{N, T}, a2::Val{A2},
               x3::SVector{N, T}, a3::Val{A3},
               x4::SVector{N, T}, a4::Val{A4},
               x5::SVector{N, T}, a5::Val{A5}
            ) where {N, T, A1, A2, A3, A4, A5}
   exprs = Expr[]
   for n = 1:N
      push!(exprs, parse("x1_$n = x1[$n]"))
      push!(exprs, parse("x2_$n = x2[$n]"))
      push!(exprs, parse("x3_$n = x3[$n]"))
      push!(exprs, parse("x4_$n = x4[$n]"))
      push!(exprs, parse("x5_$n = x5[$n]"))
   end
   for (i1, i2, i3, i4, i5) in zip(A1, A2, A3, A4, A5)
      push!(exprs, parse("m = muladd(x1_$i1 * x2_$i2, x3_$i3 * x4_$i4 * x5_$i5, m)"))
   end
   quote
      m = zero(eltype(x1))
      @inbounds $(Expr(:block, exprs...))
      return m
   end
end


@generated function gmon_d(
         x1::SVector{N, T}, dx1::SVector{N,T}, a1::Val{A1},
         x2::SVector{N, T}, dx2::SVector{N,T}, a2::Val{A2},
         x3::SVector{N, T}, dx3::SVector{N,T}, a3::Val{A3}
      ) where {N, T, A1, A2, A3}
   @assert typeof(A1) == typeof(A2) == typeof(A3) <: NTuple
   exprs = Expr[]
   for n = 1:N
      push!(exprs, parse("x1_$n = x1[$n]"))
      push!(exprs, parse("x2_$n = x2[$n]"))
      push!(exprs, parse("x3_$n = x3[$n]"))
      push!(exprs, parse("dx1_$n = dx1[$n]"))
      push!(exprs, parse("dx2_$n = dx2[$n]"))
      push!(exprs, parse("dx3_$n = dx3[$n]"))
      push!(exprs, parse("dm_$n = zero($T)"))
   end
   for (i1, i2, i3) in zip(A1, A2, A3)
      # m += x1_i1 * x2_i2 * x3_i3
      push!(exprs, parse("dm_$i1 = muladd(dx1_$i1 *  x2_$i2,  x3_$i3, dm_$i1)"))
      push!(exprs, parse("dm_$i2 = muladd( x1_$i1 * dx2_$i2,  x3_$i3, dm_$i2)"))
      push!(exprs, parse("dm_$i3 = muladd( x1_$i1 *  x2_$i2, dx3_$i3, dm_$i3)"))
   end
   coll = "SVector(" * prod("dm_$n, " for n = 1:N) * ")"
   quote
      @inbounds $(Expr(:block, exprs...))
      $(parse(coll))
   end
end

@generated function gmon_d(
         x1::SVector{N, T}, dx1::SVector{N,T}, a1::Val{A1},
         x2::SVector{N, T}, dx2::SVector{N,T}, a2::Val{A2},
         x3::SVector{N, T}, dx3::SVector{N,T}, a3::Val{A3},
         x4::SVector{N, T}, dx4::SVector{N,T}, a4::Val{A4}
      ) where {N, T, A1, A2, A3, A4}
   exprs = Expr[]
   for n = 1:N
      push!(exprs, parse("x1_$n = x1[$n]"))
      push!(exprs, parse("dx1_$n = dx1[$n]"))
      push!(exprs, parse("x2_$n = x2[$n]"))
      push!(exprs, parse("dx2_$n = dx2[$n]"))
      push!(exprs, parse("x3_$n = x3[$n]"))
      push!(exprs, parse("dx3_$n = dx3[$n]"))
      push!(exprs, parse("x4_$n = x4[$n]"))
      push!(exprs, parse("dx4_$n = dx4[$n]"))
      push!(exprs, parse("dm_$n = zero($T)"))
   end
   for (i1, i2, i3, i4) in zip(A1, A2, A3, A4)
      push!(exprs, parse("a12 = x1_$i1 * x2_$i2"))
      push!(exprs, parse("a34 = x3_$i3 * x4_$i4"))
      push!(exprs, parse("dm_$i1 = muladd(dx1_$i1 *  x2_$i2,  a34, dm_$i1)"))
      push!(exprs, parse("dm_$i2 = muladd( x1_$i1 * dx2_$i2,  a34, dm_$i2)"))
      push!(exprs, parse("dm_$i3 = muladd( a12, dx3_$i3 *  x4_$i4, dm_$i3)"))
      push!(exprs, parse("dm_$i4 = muladd( a12,  x3_$i3 * dx4_$i4, dm_$i4)"))
   end
   coll = "SVector(" * prod("dm_$n, " for n = 1:N) * ")"
   quote
      @inbounds $(Expr(:block, exprs...))
      $(parse(coll))
   end
end

@generated function gmon_d(
         x1::SVector{N, T}, dx1::SVector{N,T},
         x2::SVector{N, T}, dx2::SVector{N,T},
         x3::SVector{N, T}, dx3::SVector{N,T},
         x4::SVector{N, T}, dx4::SVector{N,T},
         x5::SVector{N, T}, dx5::SVector{N,T},
         a::Val{A}
      ) where {N, T, A}
   exprs = Expr[]
   for n = 1:N
      push!(exprs, parse("x1_$n = x1[$n]"))
      push!(exprs, parse("dx1_$n = dx1[$n]"))
      push!(exprs, parse("x2_$n = x2[$n]"))
      push!(exprs, parse("dx2_$n = dx2[$n]"))
      push!(exprs, parse("x3_$n = x3[$n]"))
      push!(exprs, parse("dx3_$n = dx3[$n]"))
      push!(exprs, parse("x4_$n = x4[$n]"))
      push!(exprs, parse("dx4_$n = dx4[$n]"))
      push!(exprs, parse("x5_$n = x5[$n]"))
      push!(exprs, parse("dx5_$n = dx5[$n]"))
      push!(exprs, parse("dm_$n = zero($T)"))
   end
   for (i1, i2, i3, i4, i5) in zip(A[1], A[2], A[3], A[4], A[5])
      push!(exprs, parse("a12 = x1_$i1 * x2_$i2"))
      push!(exprs, parse("a123 = a12 * x3_$i3"))
      push!(exprs, parse("a45 = x4_$i4 * x5_$i5"))
      push!(exprs, parse("a345 = x3_$i3 * a45"))
      push!(exprs, parse("dm_$i1 = muladd(dx1_$i1 *  x2_$i2,  a345, dm_$i1)"))
      push!(exprs, parse("dm_$i2 = muladd( x1_$i1 * dx2_$i2,  a345, dm_$i2)"))
      push!(exprs, parse("dm_$i3 = muladd( a12, dx3_$i3 *  a45, dm_$i3)"))
      push!(exprs, parse("dm_$i4 = muladd( a123,  dx4_$i4 * x5_$i5, dm_$i4)"))
      push!(exprs, parse("dm_$i5 = muladd( a123,  x4_$i4 * dx5_$i5, dm_$i5)"))
   end
   coll = "SVector(" * prod("dm_$n, " for n = 1:N) * ")"
   quote
      @inbounds $(Expr(:block, exprs...))
      $(parse(coll))
   end
end



# run benchmark tests

const I1 = @SVector [1,1,1,2,1,1,1,2,2,3,4,3,4,2,3,4,5,5,6,7,]
const I2 = @SVector [2,2,3,3,5,5,6,5,5,6,7,6,7,8,8,9,6,8,8,9,]
const I3 = @SVector [4,3,4,4,7,6,7,9,8,10,10,8,9,9,10,10,7,9,10,10,]
const I4 = @SVector [4,4,4,4,4,4,7,7,7,8,8,8,9,8,9,9,9,9,6,6,]
const I5 = @SVector [8,9,10,6,7,5,8,9,10,9,9,10,10,10,10,10,10,10,7,7]

const II1 = I1.data
const II2 = I2.data
const II3 = I3.data
const II4 = I4.data
const II5 = I5.data

const I12 = Val((II1, II2))
const I123 = Val((II1, II2, II3))
const I1234 = Val((II1, II2, II3, II4))
const I12345 = Val((II1, II2, II3, II4, II5))


const J1 = Val(I1.data)
const J2 = Val(I2.data)
const J3 = Val(I3.data)
const J4 = Val(I4.data)
const J5 = Val(I5.data)

x = @SVector rand(10)
dx = @SVector rand(10)


println("Third-Order:")
print("  fmon: "); @btime fmon($(x,x,x), $I123)
print("  gmon: "); @btime gmon($x, $J1, $x, $J2, $x, $J3)
print("fmon_d: "); @btime fmon_d($(x,x,x), $(x,x,x), $I123)
print("gmon_d: "); @btime gmon_d($x, $dx, $J1, $x, $dx, $J2, $x, $dx, $J3)

println("Fourth-Order:")
print("  fmon: "); @btime fmon($(x, x, x, x), $I1234)
print("  gmon: "); @btime gmon($x, $J1, $x, $J2, $x, $J3, $x, $J4)
print("fmon_d: "); @btime fmon_d($(x, x, x, x), $(x, x, x, x), $I1234)
print("gmon_d: "); @btime gmon_d($x, $dx, $J1, $x, $dx, $J2, $x, $dx, $J3, $x, $dx, $J4)

println("Fifth-Order:")
print("  fmon: "); @btime fmon($(x, x, x, x, x), $I12345)
print("  gmon: "); @btime gmon($x, $J1, $x, $J2, $x, $J3, $x, $J4, $x, $J5)
print("fmon_d: "); @btime fmon_d($(x, x, x, x, x), $(x, x, x, x, x), $I12345)
print("gmon_d: "); @btime gmon_d($x, $dx, $x, $dx, $x, $dx, $x, $dx, $x, $dx,
                                 Val((I1, I2, I3, I4, I5)))
