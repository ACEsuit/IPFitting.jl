using BenchmarkTools, StaticArrays

@inline fmon(x1, x2, I1, I2) = dot(x1[I1], x2[I2])

@inline function fmon(x1, x2, x3, I1, I2, I3)
   return dot(x1[I1], x2[I2] .* x3[I3])
end

@inline function fmon_d(x1, dx1, I1,
                        x2, dx2, I2,
                        x3, dx3, I3, t)
   fill!(t, 0.0)
   for n = 1:length(I1)
      @inbounds t[I1[n]] += dx1[I1[n]]* x2[I2[n]] * x3[I3[n]]
      @inbounds t[I2[n]] += x1[I2[n]] *dx2[I2[n]] * x3[I3[n]]
      @inbounds t[I3[n]] += x1[I2[n]] * x2[I2[n]] *dx3[I3[n]]
   end
   return SVector(t)
end


@inline function fmon(x1, x2, x3, x4, I1, I2, I3, I4)
   return dot(x1[I1] .* x2[I2], x3[I3] .* x4[I4])
end

@inline function fmon(x1, x2, x3, x4, x5, I1, I2, I3, I4, I5)
   return dot(x1[I1] .* x2[I2], x3[I3] .* x4[I4] .* x5[I5])
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
         x4::SVector{N, T}, dx4::SVector{N,T}, a4::Val{A4},
         x5::SVector{N, T}, dx5::SVector{N,T}, a5::Val{A5}
      ) where {N, T, A1, A2, A3, A4, A5}
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
   for (i1, i2, i3, i4, i5) in zip(A1, A2, A3, A4, A5)
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

# @generated function gmon_d(
#          x1::SVector{N, T}, dx1::SVector{N,T}, a1::Val{A1},
#          x2::SVector{N, T}, dx2::SVector{N,T}, a2::Val{A2},
#          x3::SVector{N, T}, dx3::SVector{N,T}, a3::Val{A3},
#          x4::SVector{N, T}, dx4::SVector{N,T}, a4::Val{A4},
#          x5::SVector{N, T}, dx5::SVector{N,T}, a5::Val{A5}
#       ) where {N, T, A1, A2, A3, A4, A5}
#    exprs = Expr[]
#    for n = 1:N
#       push!(exprs, parse("x1_$n = x1[$n]"))
#       push!(exprs, parse("dx1_$n = dx1[$n]"))
#       push!(exprs, parse("x2_$n = x2[$n]"))
#       push!(exprs, parse("dx2_$n = dx2[$n]"))
#       push!(exprs, parse("x3_$n = x3[$n]"))
#       push!(exprs, parse("dx3_$n = dx3[$n]"))
#       push!(exprs, parse("x4_$n = x4[$n]"))
#       push!(exprs, parse("dx4_$n = dx4[$n]"))
#       # push!(exprs, parse("x5_$n = x5[$n]"))
#       # push!(exprs, parse("dx5_$n = dx5[$n]"))
#       push!(exprs, parse("dm_$n = zero($T)"))
#    end
#    for (i1, i2, i3, i4, i5) in zip(A1, A2, A3, A4, A5)
#       # m += x1_i1 * x2_i2 * x3_i3 * x4_i4 * x5_i5
#       # push!(exprs, parse("a12 = x1_$i1 * x2_$i2"))
#       # push!(exprs, parse("a123 = a12 * x3_$i3"))
#       # push!(exprs, parse("a45 = x4_$i4 * x5_$i5"))
#       # push!(exprs, parse("a345 = x3_$i3 * a45"))
#       # push!(exprs, parse("dm_$i1 = muladd(dx1_$i1 *  x2_$i2,  a345, dm_$i1)"))
#       # push!(exprs, parse("dm_$i2 = muladd( x1_$i1 * dx2_$i2,  a345, dm_$i2)"))
#       # push!(exprs, parse("dm_$i3 = muladd( a12, dx3_$i3 *  a45, dm_$i3)"))
#       # push!(exprs, parse("dm_$i4 = muladd( a123,  dx4_$i4 * x5_$i5, dm_$i4)"))
#       # push!(exprs, parse("dm_$i5 = muladd( a123,  x4_$i4 * dx5_$i5, dm_$i5)"))
#       push!(exprs, parse("a = x1_$i1 * x2_$i2"))
#       push!(exprs, parse("b = x3_$i3 * x4_$i4"))
#       push!(exprs, parse("dm_$i1 = muladd(dx1_$i1 *  x2_$i2,  b, dm_$i1)"))
#       push!(exprs, parse("dm_$i2 = muladd( x1_$i1 * dx2_$i2,  b, dm_$i2)"))
#       push!(exprs, parse("dm_$i3 = muladd( a, dx3_$i3 *  x4_$i4, dm_$i3)"))
#       push!(exprs, parse("dm_$i4 = muladd( a,  x3_$i3 * dx4_$i4, dm_$i4)"))
#    end
#    coll = "SVector(" * prod("dm_$n, " for n = 1:N) * ")"
#    quote
#       @inbounds $(Expr(:block, exprs...))
#       $(parse(coll))
#    end
# end

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

@btime gmon_d($x, $dx, $x, $dx, $x, $dx, $x, $dx, $x, $dx,
              $(Val((I1, I2, I3, I4, I5))))


# run benchmark tests

const I1 = @SVector [1,1,1,2,1,1,1,2,2,3,4,3,4,2,3,4,5,5,6,7,]
const I2 = @SVector [2,2,3,3,5,5,6,5,5,6,7,6,7,8,8,9,6,8,8,9,]
const I3 = @SVector [4,3,4,4,7,6,7,9,8,10,10,8,9,9,10,10,7,9,10,10,]
const I4 = @SVector [4,4,4,4,4,4,7,7,7,8,8,8,9,8,9,9,9,9,6,6,]
const I5 = @SVector [8,9,10,6,7,5,8,9,10,9,9,10,10,10,10,10,10,10,7,7]

const J1 = Val(I1.data)
const J2 = Val(I2.data)
const J3 = Val(I3.data)
const J4 = Val(I4.data)
const J5 = Val(I5.data)

x = @SVector rand(10)
dx = @SVector rand(10)

println("Third-Order:")
print("  fmon: "); @btime fmon($x, $x, $x, $I1, $I2, $I3)
print("  gmon: "); @btime gmon($x, $J1, $x, $J2, $x, $J3)
print("gmon_d: "); @btime gmon_d($x, $dx, $J1, $x, $dx, $J2, $x, $dx, $J3)

println("Fourth-Order:")
print("  fmon: "); @btime fmon($x, $x, $x, $x, $I1, $I2, $I3, $I4)
print("  gmon: "); @btime gmon($x, $J1, $x, $J2, $x, $J3, $x, $J4)
print("gmon_d: "); @btime gmon_d($x, $dx, $J1, $x, $dx, $J2, $x, $dx, $J3, $x, $dx, $J4)

println("Fifth-Order:")
print("  fmon: "); @btime fmon($x, $x, $x, $x, $x, $I1, $I2, $I3, $I4, $I5)
print("  gmon: "); @btime gmon($x, $J1, $x, $J2, $x, $J3, $x, $J4, $x, $J5)
print("gmon_d: "); @btime gmon_d($x, $dx, $x, $dx, $x, $dx, $x, $dx, $x, $dx,
                                 Val((I1, I2, I3, I4, I5)))
