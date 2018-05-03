# Irreducible secondaries for NBody=4and deg=10 
PV1_1 = @SVector [x[1], x[1], x[2], x[1], x[1], x[2], x[3], x[2], x[3], x[4], x[4], x[5] ] 
PV1_2 = @SVector [x[2], x[3], x[3], x[4], x[5], x[4], x[5], x[6], x[6], x[5], x[6], x[6] ] 
PV1_1_2 =PV1_1.*PV1_1
PV1_2_2 =PV1_2.*PV1_2
PV1 = dot(PV1_1_2,PV1_2)+dot(PV1_1,PV1_2_2)
PV2_1 = @SVector [x[1], x[1], x[2], x[1], x[1], x[2], x[3], x[2], x[3], x[4], x[4], x[5] ] 
PV2_2 = @SVector [x[2], x[3], x[3], x[4], x[5], x[4], x[5], x[6], x[6], x[5], x[6], x[6] ] 
PV2_1_2 =PV2_1.*PV2_1
PV2_2_2 =PV2_2.*PV2_2
PV2_1_3 =PV2_1_2.*PV2_1
PV2_2_3 =PV2_2_2.*PV2_2
PV2 = dot(PV2_1_3,PV2_2)+dot(PV2_1,PV2_2_3)
PV3_1 = @SVector [x[1], x[2], x[3], x[4], x[5], x[6] ] 
PV3 = sum(x5)
