function invariants_Q6(x)

   return SVector{6, T}(
      # ---------------------------- primary invariants
      # I1
      (Q[1]),
      # I2
      (Q2[2] + Q2[3] + Q2[4]),
      # I3
      (Q[2] * Q[3] * Q[4]),
      # I4
      (Q2_34 + Q2_24 + Q2_23),
      # I5
      (Q2[5] + Q2[6]),
      # I6
      (Q[6] * (Q2[6] - 3*Q2[5]))
   ),
      # ---------------------------- secondary invariants
   SVector{6, T}(
      # sneak in an additional secondary "invariant"
      1.0,
      # I7
      Q[6] * (2*Q2[2] - Q2[3] - Q2[4]) + rt3 * Q[5] * (Q2[3] - Q2[4]),
      # I8
      (Q2[6] - Q2[5]) * (2*Q2[2] - Q2[3] - Q2[4]) - 2 * rt3 * Q_56 * (Q2[3] - Q2[4]),
      # I9
      Q[6] * (2*Q2_34 - Q2_24 - Q2_23) + rt3 * Q[5] * (Q2_24 - Q2_23),
      # I10
      (Q2[6] - Q2[5])*(2*Q2_34 - Q2_24 - Q2_23) - 2 * rt3 * Q_56 * (Q2_24 - Q2_23),
      # I11
      (Q2[3] - Q2[4]) * (Q2[4] - Q2[2]) * (Q2[2] - Q2[3]) * Q[5] * (3*Q2[6] - Q2[5])
   )
end

using StaticArrays
x = [1,2,3,4,5,6]
Primary_invariants= @SVector [
   x[1] + x[2] + x[3] + x[4] + x[5] + x[6],
   x[1]*x[6] + x[2]*x[5] + x[3]*x[4], x[1]^2 +
   x[2]^2 + x[3]^2 + x[4]^2 + x[5]^2 + x[6]^2, x[1]*x[2]*x[3] + x[1]*x[4]*x[5] + x[2]*x[4]*x[6] + x[3]*x[5]*x[6], x[1]^3 + x[2]^3 + x[3]^3 + x[4]^3 + x[5]^3 + x[6]^3, x[1]^4 + x[2]^4 + x[3]^4 + x[4]^4 + x[5]^4 + x[6]^4
   ]


   using StaticArrays
   x = [1,2,3,4,5,6]

   pv=zeros(3,1)
   v=zeros(6,1)

   ww1 = zeros(4,1)
