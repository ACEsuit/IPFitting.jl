
using BenchmarkTools, StaticArrays, ForwardDiff


function invariants_Q6(x::SVector{6, T}) where {T}
   x2 = x.*x
   x3 = x2.*x
   x4 = x3.*x
   x5 = x4.*x

   I3_1 = @SVector [x[1], x[2], x[3] ]
   I3_2 = @SVector [x[6], x[5], x[4] ]

   I4_1 = @SVector [x[1], x[1], x[2], x[3] ]
   I4_2 = @SVector [x[2], x[4], x[4], x[5] ]
   I4_3 = @SVector [x[3], x[5], x[6], x[6] ]

   PV1_1 = @SVector [x[1], x[1], x[2], x[1], x[1], x[2], x[3], x[2], x[3], x[4], x[4], x[5] ]
   PV1_2 = @SVector [x[2], x[3], x[3], x[4], x[5], x[4], x[5], x[6], x[6], x[5], x[6], x[6] ]

   PV1_1_2 = PV1_1.*PV1_1
   PV1_2_2 = PV1_2.*PV1_2

   PV1_1_3 = PV1_1_2.*PV1_1
   PV1_2_3 = PV1_2_2.*PV1_2

   PV1 = dot(PV1_1_2, PV1_2)+dot(PV1_1, PV1_2_2)
   PV2 = dot(PV1_1_3, PV1_2)+dot(PV1_1, PV1_2_3)
   PV3 = sum(x5)

   return SVector(
      # ---------------------------- primary invariants
      # I1
      sum(x),
      # I2
      dot(I3_1, I3_2),
      # I3
      sum(x2),
      # I4
      dot((I4_1.*I4_2), I4_3),
      # I5
      sum(x3),
      # I6
      sum(x4),
      # sneak in an additional secondary "invariant"
      1.0,
      # I7
      PV1,
      # I8
      PV2,
      # I9
      PV3,
      # I10
      PV1*PV1,
      # I11
      PV2*PV3
   )
end


invariants_d11(r::SVector{6}) = ForwardDiff.jacobian(invariants_Q6, r)

x = @SVector rand(6)
@btime invariants_Q6($x)
@btime invariants_d11($x)


using StaticArrays, BenchmarkTools
using DynamicPolynomials: @polyvar
import StaticPolynomials, ForwardDiff
@polyvar x1 x2 x3 x4 x5 x6

const PI2 = StaticPolynomials.Polynomial(
      x1*x6 + x2*x5 + x3*x4 )

const PI4 = StaticPolynomials.Polynomial(
      x1*x2*x3 + x1*x4*x5 + x2*x4*x6 + x3*x5*x6 )

const A = @SMatrix [0 1 1 1 1 0
                    1 0 1 1 0 1
                    1 1 0 0 1 1
                    1 1 0 0 1 1
                    1 0 1 1 0 1
                    0 1 1 1 1 0]

function invariants_co(x::SVector{6, T}) where {T}
   x2 = x.*x
   x3 = x2.*x
   x4 = x3.*x


   I2 = x[1]*x[6] + x[2]*x[5] + x[3]*x[4]
   I4 = x[1]*x[2]*x[3] + x[1]*x[4]*x[5] + x[2]*x[4]*x[6] + x[3]*x[5]*x[6]

   Ax = A*x
   PV1 = dot(x2, Ax)
   PV2 = dot(x3, Ax)
   PV3 = dot(x4, x)
   I11 = PV1 * PV1
   I12 = PV2 * PV3
   return SVector(
      sum(x), I2, sum(x2), I4, sum(x3), sum(x4), 1.0, PV1, PV2, PV3, I11, I12
   )
end

function invariants_co_d(x::SVector{6, T}) where {T}
   x2 = x.*x
   x3 = x2.*x
   x4 = x3.*x
   o = @SVector ones(6)
   z = @SVector zeros(6)
   ∇I2 = @SVector [x[6], x[5], x[4], x[3], x[2], x[1]]
   #∇I4 = StaticPolynomials.gradient(PI4, x)
   # x[1]*x[2]*x[3] + x[1]*x[4]*x[5] + x[2]*x[4]*x[6] + x[3]*x[5]*x[6]
   ∇I4 = @SVector [ x[2]*x[3]+x[4]*x[5],
                    x[1]*x[3]+x[4]*x[6],
                    x[1]*x[2]+x[5]*x[6],
                    x[1]*x[5]+x[2]*x[6],
                    x[1]*x[4]+x[3]*x[6],
                    x[2]*x[4]+x[3]*x[5] ]
   Ax = A*x
   PV1 = dot(x2, Ax)
   PV2 = dot(x3, Ax)
   PV3 = dot(x4, x)
   ∇PV1 = A * x2 + 2 * (x .* Ax)
   ∇PV2 = A * x3 + 3 * (x2 .* Ax)
   ∇PV3 = 5 * x4

   return hcat(o, ∇I2, 2*x, ∇I4, 3*x2, 4*x3, z, ∇PV1, ∇PV2,
               ∇PV3, 2*PV1*∇PV1, PV3*∇PV2 + PV2*∇PV3)'
end

invariants_co_ad(r::SVector{6}) = ForwardDiff.jacobian(invariants_co, r)

x = @SVector rand(6)
@btime invariants_co($x)
@btime invariants_co_d($x)
@btime invariants_co_ad($x)

invariants_co(x) ≈ invariants_Q6(x)
invariants_co_d(x) ≈ invariants_co_ad(x)

import NBodyIPs

for n = 1:10
   r = 1.0 + (@SVector rand(6))
   I = invariants_co(r)
   @assert invariants_co_d(r) ≈ invariants_co_ad(r)
   for rπ in NBodyIPs.simplex_permutations(r)
      Iπ = invariants_co(SVector{6,Float64}(rπ...))
      @assert I ≈ Iπ
   end
   print(".")
end
println()
