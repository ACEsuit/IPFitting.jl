
using BenchmarkTools, NBodyIPs, StaticArrays

function invariants_Q6(x::SVector{6, T}) where {T}
   x2 = x.*x;
   x3 = x2.*x;
   x4 = x3.*x;
   x5 = x4.*x;

   I3_1 = SVector(x[1], x[2], x[3] )
   I3_2 = @SVector [x[6], x[5], x[4] ]

   I4_1 = @SVector [x[1], x[1], x[2], x[3] ]
   I4_2 = @SVector [x[2], x[4], x[4], x[5] ]
   I4_3 = @SVector [x[3], x[5], x[6], x[6] ]

   PV1_1 = SVector(x[1], x[1], x[2], x[1], x[1], x[2], x[3], x[2], x[3], x[4], x[4], x[5] )
   PV1_2 = @SVector [x[2], x[3], x[3], x[4], x[5], x[4], x[5], x[6], x[6], x[5], x[6], x[6] ]

   PV1_1_2 = PV1_1.*PV1_1
   PV1_2_2 = PV1_2.*PV1_2

   PV1_1_3 = PV1_1_2.*PV1_1
   PV1_2_3 = PV1_2_2.*PV1_2

   PV1 = sum(PV1_1_2.*PV1_2)+sum(PV1_1.*PV1_2_2)
   PV2 = sum(PV1_1_3.*PV1_2)+sum(PV1_1.*PV1_2_3)
   PV3 = sum(x5)


   return SVector{6, T}(
      # ---------------------------- primary invariants
      # I1
      (sum(x)),
      # I2
      (sum(I3_1.*I3_2)),
      # I3
      (sum(x2)),
      # I4
      (sum((I4_1.*I4_2).*I4_3)),
      # I5
      (sum(x3)),
      # I6
      (sum(x4))
   ),
      # ---------------------------- secondary invariants
   SVector{6, T}(
      # sneak in an additional secondary "invariant"
      1.0,
      # I7
      (PV1),
      # I8
      (PV2),
      # I9
      (PV3),
      # I10
      (PV1*PV1),
      # I11
      (PV2*PV3)
   )
end

function invariants_d11(r::SVector{6, T}) where {T}
   J12 = ForwardDiff.jacobian( r_ -> vcat(invariants_Q6(r_)...), r )
   return J12
end

x = @SVector rand(6)
@btime invariants_Q6($x)
@btime invariants_d11($x)

invariants_d(x)

@btime invariants($x)
@btime invariants_d($x)
