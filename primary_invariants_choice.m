//Choice for the primary invariants

//Group symmetry - body order ([4] for 4-body, [5] for 5-body, etc.)
n := 5;

//Degree (up to)
deg := 7;

//Attach package from Braams
Attach("pack_modified.m");

//Define the field
K := RationalField();
//Define the symmetry group
G := MolSymGen([n]);
//Define the invariant ring
R0 := InvariantRing(G, K);
_:=HilbertSeries(R0);
R:=PolynomialRing(R0);

d:=R.1;
e:=&+[R.(Rank(R)-n+1+i):i in [1..n-1]];
f:=R.1+R.6;
g:=R.1+R.2+R.3+R.4;
if n eq 2 then
 Prims0:=[d];
elif n eq 3 then
 Prims0:=[d,e^2,e^3];
elif n eq 4 then
 //Prims0:=[d,e^2,d^2,e^3,d^3,e^4];
 //Prims0:=[d,f^2,d^2,f^3,d^3,f^4];
 //Prims0:=[d,R.1*R.6,d^2,f^3,d^3,f^4];
 Prims0:=[d,R.1*R.6,d^2,R.1*R.2*R.3,d^3,d^4];
 elif n eq 5 then
 //Prims0:=[d,e^2,d^2,e^3,d^3,e^4,d^4,e^5,d^5,d^6];
 //Prims0:=[d,f^2,d^2,f^3,d^3,f^4,d^4,f^5,d^5,d^6];
 //Prims0:=[d,R.1*R.5,d^2,R.1*R.2*R.3,d^3,R.1*R.5*R.8*R.10,d^4,R.1*R.2*R.3*R.4*R.5,d^5,d^6];
 Prims0:=[d,R.1*R.2,d^2,R.1*R.2*R.3,d^3,R.1*R.2*R.3*R.4,d^4,R.1*R.2*R.3*R.4*R.5,d^5,d^6];
 elif n eq 6 then
 // just guessing ...
 Prims0:=[d,e^2,d^2,e^3,d^3,R.1^2*R.2,R.1*R.2*R.3,
    e^4,d^4,R.1^3*R.2,e^5,d^5,R.1^4*R.2,e^6,d^6];
end if;
R0`PrimaryInvariants:=[&+[f^g:g in Group(R0)]:f in Prims0];



R0;
RR1 := SecondaryInvariants(R0);
RR1;
RR2 := IrreducibleSecondaryInvariants(R0);
RR2;

A, Q := Algebra(R0);
A;
Q;


//Primary and secondary invariants
logfile := "logNbody_Braams";
SetLogFile(logfile: Overwrite := true);
Attach("pack_modified.m");
PSI :=MolInvRngGen([5],7);
PSI`PrimaryInvariants;



#PSI`SecondaryInvariants;
#PSI`IrreducibleSecondaries;
PSI`IrreducibleSecondaryInvariants := PSI`IrreducibleSecondaries;

A, Q := Algebra(PSI);
A;
Q;



////////////
 Secs, IrrSecs := MolInvSecsToDegree(5,
   PSI:
   IniPrims:=PSI`PrimaryInvariants,
   Names:=[x1,x2,x3,x4,x5,x6,x7,x8,x9,x10])


intrinsic MolInvSecsToDegree
 (deg::RngIntElt,
  R0::RngInvar:
  IniPrims:=[],
  IniIrrSecs:=[],
  Names:=[]) ->
 [], []
{Compute Secondary Invariants and Irreducible Secondaries to degree deg}
R:=PolynomialRing(CoefficientRing(P),Rank(P),"grevlex")
  where P is PolynomialRing(R0);
if Names ne [] then
 AssignNames(~R,Names);
end if;
RedSecsReps:=[[IntegerRing()|]];
if IsVerbose("MolInv") then
 MolInvF95PrintReps(0,RedSecsReps);
end if;
if IniPrims eq [] then
 Prims:=[R!f:f in R0`PrimaryInvariants];
else
 Prims:=[R!f:f in IniPrims];
end if;
IrrSecs:=[R|] cat [R!f:f in IniIrrSecs];
Secs:=[R!1] cat IrrSecs;
// Note, we don't "F95Print" the IniIrrSecs
for d:=1 to deg do
 //time NewRSR:=MolInvNewRedSecsReps(d,Prims,IrrSecs,RedSecsReps);
 NewRSR:=MolInvNewRedSecsReps(d,Prims,IrrSecs,RedSecsReps);
 if IsVerbose("MolInv") then
  MolInvF95PrintReps(#Secs,NewRSR);
 end if;
 RedSecsReps:=RedSecsReps cat NewRSR;
 Secs:=Secs cat [&*IrrSecs[t]:t in NewRSR];
 //time TryIrrSecs:=Sort([R!f:f in InvariantsOfDegree(R0,d)]);
 TryIrrSecs:=Sort([R!f:f in InvariantsOfDegree(R0,d)]);
 //time l:=HomogeneousModuleTestBasis(Prims,Secs,TryIrrSecs);
 l:=HomogeneousModuleTestBasis(Prims,Secs,TryIrrSecs);
 NewIS:=Reverse(TryIrrSecs[l]);
 if IsVerbose("MolInv") then
  MolInvF95PrintIrrSecs(#IrrSecs,#Secs,NewIS);
 end if;
 IrrSecs:=IrrSecs cat NewIS;
 Secs := Secs cat NewIS;
end for;
return Secs, IrrSecs;
end intrinsic;
