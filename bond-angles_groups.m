//4-body

//Define the field
K := RationalField();
//Define the symmetry group
G := PermutationGroup<6| (2,3,1)(5,6,4)>;
//Define the invariant ring
R := InvariantRing(G, K);
R;

//Molien series
M<t> := MolienSeries(G);
M;

R0 := PrimaryInvariants(R);
R0;
[TotalDegree(R0[i]): i in [1..#R0]];

R1 := SecondaryInvariants(R);
R1;
[TotalDegree(R1[i]): i in [1..#R1]];

F := FundamentalInvariants(R);
F;
[TotalDegree(F[i]): i in [1..#F]];

//For 5-body

//Define the field
K := RationalField();
//Define the symmetry group
G := PermutationGroup<10 | (2,3)(5,6)(10,9), (2,3,4)(5,6,7)(8,10,9), (1,3,4,2)(5,6,10,9)(7,8)>;
//Define the invariant ring
R := InvariantRing(G, K);
R;

//Molien series
M<t> := MolienSeries(G);
M;

R0 := PrimaryInvariants(R);
R0;
[TotalDegree(R0[i]): i in [1..#R0]];

R1 := SecondaryInvariants(R);
R1;
[TotalDegree(R1[i]): i in [1..#R1]];

R2 := IrreducibleSecondaryInvariants(R);
R2;
[TotalDegree(R2[i]): i in [1..#R2]];

F := FundamentalInvariants(R);
F;
[TotalDegree(F[i]): i in [1..#F]];

RR2 := IrreducibleSecondaryInvariants(R);
RR2;

A, Q := Algebra(R);
A;
Q;


// //Primary and secondary invariants
// PSI :=MolInvRngGen([nbody],deg);
//
// S<t>:=PowerSeriesRing(Rationals(),deg);
// Write(outputfile, "Group degree and order:");
// Write(outputfile, Degree(Group(PSI))) ;
// Write(outputfile, #Group(PSI)) ;
// dnpr:=[0:i in [0..Precision(S)-1]];
// for f in PSI`PrimaryInvariants do
//  i:=TotalDegree(f);
//  if i lt Precision(S) then
//   dnpr[i+1]:=dnpr[i+1]+1;
//  end if;
// end for;
// Write(outputfile, "Degrees of Primaries:");
// Write(outputfile, [TotalDegree(f): f in PSI`PrimaryInvariants]);
// Write(outputfile, "Number of primaries at each degree, and sums:");
// Write(outputfile, [dnpr[1+i]:i in [0..Precision(S)-1]]);
// Write(outputfile, [&+[dnpr[1+j]:j in [0..i]]:i in [0..Precision(S)-1]]);
// dnb:=S!PSI`HilbertSeries;
// dnpb:=1/&*[1-t^TotalDegree(f):f in PSI`PrimaryInvariants];
// Write(outputfile, "Dimensions of the Primaries Ring:");
// Write(outputfile, [Coefficient(dnpb,i):i in [0..Precision(S)-1]]);
// Write(outputfile, [&+[Coefficient(dnpb,j):j in [0..i]]:i in [0..Precision(S)-1]]);
// dnsc:=dnb*&*[1-t^TotalDegree(f):f in PSI`PrimaryInvariants];
// Write(outputfile, "Expected Numbers of Secondaries, from Degree 0:");
// Write(outputfile, [Coefficient(dnsc,i):i in [0..Precision(S)-1]]);
// Write(outputfile, [&+[Coefficient(dnsc,j):j in [0..i]]:i in [0..Precision(S)-1]]);
// Write(outputfile, "Molien/Hilbert Series from Degree 0, and sums:");
// Write(outputfile, [Coefficient(dnb,i):i in [0..Precision(S)-1]]);
// Write(outputfile, [&+[Coefficient(dnb,j):j in [0..i]]:i in [0..Precision(S)-1]]);
//
//
// Write(outputfile, "Primary invariants");
// Write(outputfile, PSI`PrimaryInvariants);
//
// Write(outputfile, "Secondary invariants");
// Write(outputfile, PSI`SecondaryInvariants);
//
// Write(outputfile, "Degrees of the secondary invariants");
// Write(outputfile, [TotalDegree(PSI`SecondaryInvariants[i]): i in [1..#PSI`SecondaryInvariants]]);
//
//
// //Fundamental invariants
// F := FundamentalInvariants(R);
// Write(outputfile, "Fundamental invariants");
// Fd := [];
// count := 0;
// for f in F do
//  i:=TotalDegree(f);
//  if i le deg then
//  count := count +1;
//   Fd[count]:=f;
//  end if;
// end for;
// Fd;
// #Fd;
//
// Write(outputfile, Fd);
// Write(outputfile, "Degrees of the fundamental invariants");
// Write(outputfile, [TotalDegree(Fd[i]): i in [1..#Fd]]);
//
//
//

exit;
