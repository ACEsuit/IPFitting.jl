
//Group symmetry - body order ([4] for 4-body, [5] for 5-body, etc.)
nbody := NBODY;

//Degree (up to)
deg := DEGREE;

//Outputfile name
// outputfile := "Nbody_output.txt";
//Logfile name
logfile := "logNbody_output.txt";

//Attach package from Braams
Attach("pack_opt_primaries.m");

//Prescribe outputfile and logfile
// SetOutputFile(outputfile: Overwrite := true);
SetLogFile(logfile: Overwrite := true);

printf " Nbody group where N="*IntegerToString(nbody)*"\n";
printf " Maximal polynomial degree="*IntegerToString(deg)*"\n";

// Write(outputfile, "N-body group where N=" cat IntegerToString(nbody));
// Write(outputfile, "Maximal polynomial degree:" cat IntegerToString(deg));

//Define the field
K := RationalField();
//Define the symmetry group
G := MolSymGen([nbody]);
//Define the invariant ring
R := InvariantRing(G, K);

//Primary and secondary invariants
PSI :=MolInvRngGen([nbody],deg);

// S<t>:=PowerSeriesRing(Rationals(),deg+1);
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


// Write(outputfile, "Primary invariants");
// Write(outputfile, PSI`PrimaryInvariants);

printf " Primary_invariants=\n";
PSI`PrimaryInvariants;

// Write(outputfile, "Secondary invariants");
// Write(outputfile, PSI`SecondaryInvariants);
//
// Write(outputfile, "Degrees of the secondary invariants");
// Write(outputfile, [TotalDegree(PSI`SecondaryInvariants[i]): i in [1..#PSI`SecondaryInvariants]]);


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

// Write(outputfile, Fd);
// Write(outputfile, "Degrees of the fundamental invariants");
// Write(outputfile, [TotalDegree(Fd[i]): i in [1..#Fd]]);

// //Molien series
// M<t> := MolienSeries(G);
// Write(outputfile, "Molien series");
// Write(outputfile, M);
Sec := PSI`SecondaryInvariants;

printf " Nb_secondary_invariants=";
#Sec;

exit;
