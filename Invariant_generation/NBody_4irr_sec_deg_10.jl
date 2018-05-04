# Irreducible secondaries for NBody=4and deg=10 
A1=[0 1 1 1 1 0 ; 1 0 1 1 0 1 ; 1 1 0 0 1 1 ; 1 1 0 0 1 1 ; 1 0 1 1 0 1 ; 0 1 1 1 1 0 ]

PV1 = x2'*A1 x1
A2=[0 1 1 1 1 0 ; 1 0 1 1 0 1 ; 1 1 0 0 1 1 ; 1 1 0 0 1 1 ; 1 0 1 1 0 1 ; 0 1 1 1 1 0 ]

PV2 = x3'*A2 x1
PV3 = sum(x5)
