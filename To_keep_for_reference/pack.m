// Sample use:
// Attach("Magma/pack.m");
// _:=MolInvRngGen([3,2,1],6);

declare verbose MolInv,2;
SetVerbose("MolInv",1);
MolInvPrec:=10;

declare attributes RngInvar: IrreducibleSecondaries;

intrinsic GeneratorSequence
 (G::Grp) ->
 []
{ordered sequence of generators of G}
return [G.i:i in [1..NumberOfGenerators(G)]];
end intrinsic;

intrinsic MolSymS
 (n::RngIntElt) ->
 GrpPerm
{The complete symmetric group, by a different name}
return Sym(n);
end intrinsic;

intrinsic MolSymSNames
 (n::RngIntElt) ->
 []
{Print Names for MolSymS}
return ["d(i"*IntegerToString(i)*")":i in [0..n-1]];
end intrinsic;

intrinsic MolSymMS
 (n::RngIntElt,
  r::RngIntElt) ->
 GrpPerm
{multiple representation of the symmetric group}
pairs:=[<i,j>:i in [1..n],j in [1..r]];
return PermutationGroup<#pairs|
  [[Index(pairs,<p[1]^g,p[2]>):p in pairs]:g in GeneratorSequence(Sym(n))]>;
end intrinsic;

intrinsic MolSymMSNames
 (n::RngIntElt,
  r::RngIntElt) ->
 []
{Print Names for MolSymMS}
return ["d(i"*IntegerToString(i)*","*IntegerToString(j)*")":
  i in [0..n-1],j in [0..r-1]];
end intrinsic;

intrinsic MolSymCG
 (n::RngIntElt) ->
 GrpPerm
{Permutation group for the complete graph CG(n)}
pairs:=[{i,j}:i in [1..j-1],j in [1..n]];
return PermutationGroup<#pairs|
  [[Index(pairs,p^g):p in pairs]:g in GeneratorSequence(Sym(n))]>;
end intrinsic;

intrinsic MolSymCGNames
 (n::RngIntElt) ->
 []
{Print Names for MolSymCG}
return ["d(i"*IntegerToString(i)*",i"*IntegerToString(j)*")":
  i in [0..j-1],j in [0..n-1]];
end intrinsic;

intrinsic MolSymNewCG
 (n::RngIntElt) ->
 GrpPerm
{Permutation group for the complete graph CG(n)}
pairs:=[{i,j}:i in [1..j-1],j in [1..n]];
return PermutationGroup<1+n+#pairs|[
  [1] cat [1+i^g:i in [1..n]] cat
    [1+n+Index(pairs,p^g):p in pairs]:
  g in GeneratorSequence(Sym(n))]>;
end intrinsic;

intrinsic MolSymNewCGNames
 (n::RngIntElt) ->
 []
{Print Names for MolSymNewCG}
return
  ["b"] cat ["c(i"*IntegerToString(i)*")":i in [0..n-1]] cat
  ["d(i"*IntegerToString(i)*",i"*IntegerToString(j)*")":
    i in [0..j-1],j in [0..n-1]];
end intrinsic;

intrinsic MolSymBG
 (m::RngIntElt,
  n::RngIntElt) ->
 GrpPerm
{Permutation group for the bipartite graph BG(m,n)}
pairs:=[<i,j>:i in [1..m],j in [1..n]];
return PermutationGroup<#pairs|
  [[Index(pairs,<p[1]^g,p[2]>):p in pairs]:g in GeneratorSequence(Sym(m))] cat
  [[Index(pairs,<p[1],p[2]^g>):p in pairs]:g in GeneratorSequence(Sym(n))]>;
end intrinsic;

intrinsic MolSymBGNames
 (m::RngIntElt,
  n::RngIntElt) ->
 []
{Print Names for MolSymBG}
return ["d(i"*IntegerToString(i)*",j"*IntegerToString(j)*")":
  i in [0..m-1],j in [0..n-1]];
end intrinsic;

intrinsic MolSymGen
 (nki::[RngIntElt]) ->
 GrpPerm
{Permutation group for a Molecule}
pairs:=[{i,j}:i in s[k],j in s[l],k in [1..l],l in [1..#nki]|i lt j]
  where s is [[&+nki[1..k-1]+1..&+nki[1..k]]:k in [1..#nki]];
return PermutationGroup<#pairs|
  [[Index(pairs,p^g):p in pairs]:g in GeneratorSequence(G)]>
  where G is DirectProduct([Sym(k):k in nki]);
end intrinsic;

intrinsic MolSymGenNames
 (nki::[RngIntElt]) ->
 []
{Print Names for a Molecule}
s:=[[&+nki[1..k-1]+1..&+nki[1..k]]:k in [1..#nki]];
lets:=["i","j","k","l","m","n","o","p","q"];
w:=[lets[k]*IntegerToString(i):
  i in [0..nki[k]-1],k in [1..#nki]];
return ["d("*w[i]*","*w[j]*")":
  i in s[k],j in s[l],k in [1..l],l in [1..#nki]|i lt j];
end intrinsic;

intrinsic MolInvMerge
 (L::[[]]) ->
 []
{Merge seqs of L by ascending degree}
P:= []; b:=[0:i in [1..#L]]; c:=[#L[i]:i in [1..#L]];
while b lt c do
// append earliest element of minimal degree
 deg:=Infinity(); l:=0;
 for i:=1 to #L do
  if b[i] lt c[i] then
   if TotalDegree(L[i][b[i]+1]) lt deg then
    l:=i; deg:=TotalDegree(L[i][b[i]+1]);
   end if;
  end if;
 end for;
 Append(~P,L[l][b[l]+1]);
 b[l]+:=1;
end while;
return P;
end intrinsic;

intrinsic MolInvF95PrintDims
 (R0::RngInvar:
  prec:=MolInvPrec)
{Print some basic dimension information}
S<t>:=PowerSeriesRing(Rationals(),prec);
print "Group degree and order:", Degree(Group(R0)), #Group(R0);
dnpr:=[0:i in [0..Precision(S)-1]];
for f in R0`PrimaryInvariants do
 i:=TotalDegree(f);
 if i lt Precision(S) then
  dnpr[i+1]:=dnpr[i+1]+1;
 end if;
end for;
print "Degrees of Primaries:";
print [TotalDegree(f): f in R0`PrimaryInvariants];
print "Number of primaries at each degree, and sums:";
print [dnpr[1+i]:i in [0..Precision(S)-1]];
print [&+[dnpr[1+j]:j in [0..i]]:i in [0..Precision(S)-1]];
dnb:=S!R0`HilbertSeries;
dnpb:=1/&*[1-t^TotalDegree(f):f in R0`PrimaryInvariants];
print "Dimensions of the Primaries Ring:";
print [Coefficient(dnpb,i):i in [0..Precision(S)-1]];
print [&+[Coefficient(dnpb,j):j in [0..i]]:i in [0..Precision(S)-1]];
dnsc:=dnb*&*[1-t^TotalDegree(f):f in R0`PrimaryInvariants];
print "Expected Numbers of Secondaries, from Degree 0:";
print [Coefficient(dnsc,i):i in [0..Precision(S)-1]];
print [&+[Coefficient(dnsc,j):j in [0..i]]:i in [0..Precision(S)-1]];
print "Molien/Hilbert Series from Degree 0, and sums:";
print [Coefficient(dnb,i):i in [0..Precision(S)-1]];
print [&+[Coefficient(dnb,j):j in [0..i]]:i in [0..Precision(S)-1]];
end intrinsic;

intrinsic MolInvF95PrintReps
 (base::RngIntElt,
  Reps::[[]])
{Print RedSecsReps for editing into F95 code}
for i:=0 to #Reps-1 do
 if #Reps[i+1] eq 0 then
  printf " v(%o) = 1\n", base+i;
 else
  printf " v(%o) = pv(%o)", base+i, Reps[i+1][1]-1;
  for j:=2 to #Reps[i+1] do
   printf "*pv(%o)", Reps[i+1][j]-1;
  end for;
  printf "\n";
 end if;
end for;
end intrinsic;

intrinsic MolInvF95PrintIrrSecs
 (base0::RngIntElt,
  base1::RngIntElt,
  IrrSecs::[RngMPolElt])
{Print IrrSecs for editing into F95 code}
for i:=0 to #IrrSecs-1 do
 printf " pv(%o) = SYM %o\n",
   base0+i, LeadingMonomial(IrrSecs[i+1]);
end for;
if 0 lt #IrrSecs then
 printf " v(%o:%o) = pv(%o:%o)\n",
   base1, base1+#IrrSecs-1, base0, base0+#IrrSecs-1;
end if;
end intrinsic;

intrinsic MolInvPrimaries
 (R0::RngInvar) ->
 []
{Primaries computed following Kempers email, May 26, 2003}
K:=CoefficientRing(R0);
M:=GModule(Group(R0),K);
Gl:=MatrixGroup(M);
Pi:=[&+[K!Integers()!c(g)*MatrixAlgebra(K,Degree(Gl))!g: g in Gl]/#Gl:
    c in CharacterTable(Gl)];
Pi:=[p: p in Pi | p ne 0];
L:=[sub<M | Image(p)>: p in Pi];
L:=[M: M in L | Dimension(M) ne 1];
prims:=[**];
time
for N in L do
    prim:=PrimaryInvariants(InvariantRing(MatrixGroup(N)));
    Append(~prims,prim);
end for;
return prims;
end intrinsic;

intrinsic MolInvNewRedSecsReps
 (d::RngIntElt,
  Prims::[RngMPolElt],
  IrrSecs::[RngMPolElt],
  RedSecsReps::[[]]) ->
 []
{New Reducible Secondaries at degree d}
TryReps:=[[i] cat g:i in [1..g[1]],
  g in RedSecsReps[[2..#RedSecsReps]] cat [[i]:i in [1..#IrrSecs]]|
  TotalDegree(IrrSecs[i])+&+[TotalDegree(IrrSecs[i]):i in g] eq d];
l:=HomogeneousModuleTestBasis(Prims,
  IrrSecs cat [&*IrrSecs[t]:t in RedSecsReps],
  [&*IrrSecs[t]:t in TryReps]);
return TryReps[l];
end intrinsic;

intrinsic MolInvSecsToDegree
 (mxd::RngIntElt,
  R0::RngInvar:
  IniPrims:=[],
  IniIrrSecs:=[],
  Names:=[]) ->
 [], []
{Compute Secondary Invariants and Irreducible Secondaries to degree mxd}
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
for d:=1 to mxd do
 time NewRSR:=MolInvNewRedSecsReps(d,Prims,IrrSecs,RedSecsReps);
 if IsVerbose("MolInv") then
  MolInvF95PrintReps(#Secs,NewRSR);
 end if;
 RedSecsReps:=RedSecsReps cat NewRSR;
 Secs:=Secs cat [&*IrrSecs[t]:t in NewRSR];
 time TryIrrSecs:=Sort([R!f:f in InvariantsOfDegree(R0,d)]);
 time l:=HomogeneousModuleTestBasis(Prims,Secs,TryIrrSecs);
 NewIS:=Reverse(TryIrrSecs[l]);
 if IsVerbose("MolInv") then
  MolInvF95PrintIrrSecs(#IrrSecs,#Secs,NewIS);
 end if;
 IrrSecs:=IrrSecs cat NewIS;
 Secs := Secs cat NewIS;
end for;
return Secs, IrrSecs;
end intrinsic;

intrinsic MolInvPrimsRngS
 (n::RngIntElt) ->
 RngInvar
{Invariant Ring for MolSymS}
R0:=InvariantRing(MolSymS(n),Rationals());
_:=HilbertSeries(R0);
R:=PolynomialRing(R0);
// Assigned Primaries
Prims0:=[R.1^k:k in [1..n]];
R0`PrimaryInvariants:=[&+[f^g:g in Group(R0)]:f in Prims0];
return R0;
end intrinsic;

intrinsic MolInvPrimsRngMS
 (n::RngIntElt,
  r::RngIntElt) ->
 RngInvar
{Invariant Ring for MolSymMS}
R0:=InvariantRing(MolSymMS(n,r),Rationals());
_:=HilbertSeries(R0);
R:=PolynomialRing(R0);
// Assigned Primaries
Prims0:=[R.(j*n+1)^k:j in [0..r-1], k in [1..n]];
R0`PrimaryInvariants:=[&+[f^g:g in Group(R0)]:f in Prims0];
return R0;
end intrinsic;

intrinsic MolInvPrimsRngCG
 (n::RngIntElt) ->
 RngInvar
{Invariant Ring for MolSymCG}
R0:=InvariantRing(MolSymCG(n),Rationals());
_:=HilbertSeries(R0);
R:=PolynomialRing(R0);
// Assigned (maybe conjectured) Primaries
d:=R.1;
e:=&+[R.(Rank(R)-n+1+i):i in [1..n-1]];
if n eq 2 then
 Prims0:=[d];
elif n eq 3 then
 Prims0:=[d,e^2,e^3];
elif n eq 4 then
 Prims0:=[d,e^2,d^2,e^3,d^3,e^4];
elif n eq 5 then
 Prims0:=[d,e^2,d^2,e^3,d^3,e^4,d^4,e^5,d^5,d^6];
elif n eq 6 then
// just guessing ...
 Prims0:=[d,e^2,d^2,e^3,d^3,R.1^2*R.2,R.1*R.2*R.3,
   e^4,d^4,R.1^3*R.2,e^5,d^5,R.1^4*R.2,e^6,d^6];
elif n eq 7 then
// just guessing ...
 Prims0:=[d,e^2,d^2,e^3,d^3,R.1^2*R.2,R.1*R.2*R.3,R.1*R.3*R.4,
   e^4,d^4,R.1^3*R.2,e^5,d^5,R.1^4*R.2,e^6,d^6,
   e^7,d^7,R.1^6*R.2,d^10,d^12];
elif n eq 8 then
// totally guessing ...
 Prims0:=[d,e^2,d^2,e^3,d^3,R.1^2*R.2,R.1*R.2*R.3,R.1*R.3*R.4,
   e^4,d^4,R.1^3*R.2,e^5,d^5,R.1^4*R.2,e^6,d^6,
   e^7,d^7,R.1^6*R.2,d^10,d^12,
   d^3,d^12,d^12,d^12,d^12,d^12,d^12];
end if;
R0`PrimaryInvariants:=[&+[f^g:g in Group(R0)]:f in Prims0];
return R0;
end intrinsic;

intrinsic MolInvPrimsRngNewCG
 (n::RngIntElt) ->
 RngInvar
{Invariant Ring for MolSymNewCG}
R0:=InvariantRing(MolSymNewCG(n),Rationals());
_:=HilbertSeries(R0);
R:=PolynomialRing(R0);
// Assigned (maybe conjectured) Primaries
b:=R.1;
c:=R.2;
d:=R.(1+n+1);
e:=&+[R.(Rank(R)-n+1+i):i in [1..n-1]];
if n eq 2 then
 Prims0:=[b,c,d,c^2];
elif n eq 3 then
 Prims0:=[b,c,d,c^2,e^2,c^3,e^3];
elif n eq 4 then
 Prims0:=[b,c,d,c^2,e^2,d^2,c^3,e^3,d^3,c^4,e^4];
elif n eq 5 then
 Prims0:=[b,c,d,c^2,e^2,d^2,c^3,e^3,d^3,c^4,e^4,d^4,c^5,e^5,d^5,d^6];
elif n eq 6 then
// just guessing ...
 Prims0:=[b,c,d,c^2,e^2,d^2,c^3,e^3,d^3,R.1^2*R.2,R.1*R.2*R.3,
   c^4,e^4,d^4,R.1^3*R.2,c^5,e^5,d^5,R.1^4*R.2,c^6,e^6,d^6];
elif n eq 7 then
// just guessing ...
 Prims0:=[b,c,d,c^2,e^2,d^2,c^3,e^3,d^3,R.1^2*R.2,R.1*R.2*R.3,R.1*R.3*R.4,
   c^4,e^4,d^4,R.1^3*R.2,c^5,e^5,d^5,R.1^4*R.2,c^6,e^6,d^6,
   c^7,e^7,d^7,R.1^6*R.2,d^10,d^12];
end if;
R0`PrimaryInvariants:=[&+[f^g:g in Group(R0)]:f in Prims0];
//  Prims0 cat:=[R.4];
//  Prims0 cat:=[R.5+R.6,R.5+R.7,R.6+R.7];
//  Prims0 cat:=[R.6+R.7+R.9,R.6+R.8+R.10,R.7+R.8+R.11,R.9+R.10+R.11]
//  Prims0 cat:=[R.7+R.8+R.10+R.13,R.7+R.9+R.11+R.14,R.8+R.9+R.12+R.15,
//    R.10+R.11+R.12+R.16,R.13+R.14+R.15+R.16]
return R0;
end intrinsic;

intrinsic MolInvPrimsRngBG
 (m::RngIntElt,
  n::RngIntElt) ->
 RngInvar
{Invariant Ring for MolSymBG}
R0:=InvariantRing(MolSymBG(m,n),Rationals());
_:=HilbertSeries(R0);
R:=PolynomialRing(R0);
// Assigned (maybe conjectured) Primaries
d:=R.1;
e:=&+[R.(m*j+1):j in [0..n-1]];
f:=&+[R.(i+1):i in [0..m-1]];
if 2 le m then g:=R.2; end if;
if 2 le n then h:=R.(m+1); end if;
if m eq 1 or n eq 1 then
 Prims0:=[d^k:k in [1..Max([m,n])]];
elif m eq 2 and n eq 2 then
 Prims0:=[d,e^2,f^2,d^2];
elif m eq 3 and n eq 2 then
 Prims0:=[d,e^2,f^2,d^2,e^3,d^6];
elif m eq 2 and n eq 3 then
 Prims0:=[d,e^2,f^2,d^2,f^3,d^6];
elif m eq 4 and n eq 2 then
 Prims0:=[d,e^2,f^2,d^2,e^3,e^4,d^4,d^6];
elif m eq 2 and n eq 4 then
 Prims0:=[d,e^2,f^2,d^2,f^3,f^4,d^4,d^6];
elif m eq 5 and n eq 2 then
 Prims0:=[d,e^2,f^2,d^2,e^3,e^4,d^4,e^5,d^6,d^10];
elif m eq 2 and n eq 5 then
 Prims0:=[d,e^2,f^2,d^2,f^3,f^4,d^4,f^5,d^6,d^10];
elif m eq 6 and n eq 2 then
 Prims0:=[d,e^2,f^2,d^2,e^3,d^3,e^4,d^4,e^5,e^6,d^6,d^10];
elif m eq 2 and n eq 6 then
 Prims0:=[d,e^2,f^2,d^2,f^3,d^3,f^4,d^4,f^5,f^6,d^6,d^10];
elif m eq 7 and n eq 2 then
 Prims0:=[d,e^2,f^2,d^2,e^3,d^3,e^4,d^4,e^5,e^6,d^6,e^7,d^10,d^14];
elif m eq 2 and n eq 7 then
 Prims0:=[d,e^2,f^2,d^2,f^3,d^3,f^4,d^4,f^5,f^6,d^6,f^7,d^10,d^14];
elif m eq 3 and n eq 3 then
 Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,d^4,d^6];
elif m eq 4 and n eq 3 then
// just guessing...
 Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,e^4,d^4,d^3*g,d^6,d^12];
elif m eq 3 and n eq 4 then
// just guessing...
 Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,f^4,d^4,d^3*h,d^6,d^12];
elif m eq 5 and n eq 3 then
// totally guessing...
// Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,e^4,d^4,d^3*g,
//   d^5,d^6,d^10,d^12,d^15];
// try this
 Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,e^4,d^4,d^3*g,
   e^5,d^6,d^10,d^12,d^15];
// The following looks proper based on the Molien series alone, but
// it fails the test; one secondary too many is found at degree 7.
// Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,d^2*g,d^2*h,e^4,d^4,
//   e^5,d^5,d^6,d^60];
elif m eq 3 and n eq 5 then
// totally guessing...
// Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,f^4,d^4,d^3*h,
//   d^5,d^6,d^10,d^12,d^15];
 Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,f^4,d^4,d^3*h,
   f^5,d^6,d^10,d^12,d^15];
elif m eq 6 and n eq 3 then
// totally guessing...
 Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,e^4,d^4,d^3*g,d^3*h,
   e^5,d^5,e^6,d^6,d^10,d^12,d^15];
// The Molien series test would allow lower powers, for example the next
// sequence, but in view of the (5,3) experience I don't trust it.
// Prims0:=[d,d^2,d^2,d^2,d^3,d^3,d^3,d^3,d^3,d^3,d^4,d^4,d^4,
//   d^5,d^5,d^6,d^6,d^60];
elif m eq 3 and n eq 6 then
// totally guessing...
 Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,f^4,d^4,d^3*g,d^3*h,
   f^5,d^5,f^6,d^6,d^10,d^12,d^15];
elif m eq 4 and n eq 4 then
// totally guessing...
 Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,d^2*g,d^2*h,e^4,f^4,d^4,d^3*g,d^3*h,
   d^6,d^12];
elif m eq 5 and n eq 4 then
// totally guessing...
 Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,d^2*g,d^2*h,e^4,f^4,d^4,d^3*g,d^3*h,
   e^5,d^6,d^10,d^12,d^15,d^20];
elif m eq 4 and n eq 5 then
// totally guessing...
 Prims0:=[d,e^2,f^2,d^2,e^3,f^3,d^3,d^2*g,d^2*h,e^4,f^4,d^4,d^3*g,d^3*h,
   f^5,d^6,d^10,d^12,d^15,d^20];
end if;
R0`PrimaryInvariants:=[&+[f^g:g in Group(R0)]:f in Prims0];
return R0;
end intrinsic;

intrinsic MolInvPrimsRngGen
 (nki::[RngIntElt]:
  nkj:=nki) ->
 RngInvar
{Invariant Ring for a general Molecular Group}
G:=MolSymGen(nki);
R0:=InvariantRing(G,Rationals());
_:=HilbertSeries(R0);
R:=PolynomialRing(R0);
// Assigned (maybe conjectured) Primaries based on nkj; nkj should
// be a coarsening of nki
List:=[];
base:=0;
for l:=1 to #nkj do
 for k:=1 to l do
  if k lt l then
   R1:=MolInvPrimsRngBG(nkj[k],nkj[l]);
   r:=Rank(PolynomialRing(R1));
  elif 2 le nkj[l] then
   R1:=MolInvPrimsRngCG(nkj[l]);
   r:=Rank(PolynomialRing(R1));
  else
   r:=0;
  end if;
  if r ne 0 then
   h:=hom<PolynomialRing(R1)->R|[R.i:i in [base+1..base+r]]>;
   Append(~List,h(R1`PrimaryInvariants));
   base+:=r;
  end if;
 end for;
end for;
R0`PrimaryInvariants:=MolInvMerge(List);
return R0;
end intrinsic;

intrinsic MolInvIniIrrSecsGen
 (R::RngMPol,
  nki::[RngIntElt],
  mxd::RngIntElt) ->
 []
{Initial Secondaries for SymGen, from the component reps}
List:=[];
base:=0;
for l:=1 to #nki do
 for k:=1 to l do
  if k lt l then
   R1:=MolInvRngBG(nki[k],nki[l],mxd);
   r:=Rank(PolynomialRing(R1));
  elif 2 le nki[l] then
   R1:=MolInvRngCG(nki[l],mxd);
   r:=Rank(PolynomialRing(R1));
  else
   r:=0;
  end if;
  if r ne 0 then
   h:=hom<PolynomialRing(R1)->R|[R.i:i in [base+1..base+r]]>;
   s:=R1`SecondaryInvariants;
   Append(~List,h(s[2..#s]));
   base+:=r;
  end if;
 end for;
end for;
return MolInvMerge(List);
end intrinsic

intrinsic MolInvRngS
 (n::RngIntElt,
  mxd::RngIntElt) ->
 RngInvar
{Invariant Ring with Secondaries to degree mxd for MolSymS}
// Secondaries for Sym(n) are trivial, we know that.
R0:=MolInvPrimsRngS(n);
if IsVerbose("MolInv") then
 MolInvF95PrintDims(R0);
end if;
R0`SecondaryInvariants,R0`IrreducibleSecondaries:=
  MolInvSecsToDegree(mxd,R0:Names:=MolSymSNames(n));
return R0;
end intrinsic;

intrinsic MolInvRngMS
 (n::RngIntElt,
  r::RngIntElt,
  mxd::RngIntElt) ->
 RngInvar
{Invariant Ring with Secondaries to degree mxd for MolSymMS}
R0:=MolInvPrimsRngMS(n,r);
if IsVerbose("MolInv") then
 MolInvF95PrintDims(R0);
end if;
R0`SecondaryInvariants,R0`IrreducibleSecondaries:=
  MolInvSecsToDegree(mxd,R0:Names:=MolSymMSNames(n,r));
return R0;
end intrinsic;

intrinsic MolInvRngCG
 (n::RngIntElt,
  mxd::RngIntElt) ->
 RngInvar
{Invariant Ring with Secondaries to degree mxd for MolSymCG}
R0:=MolInvPrimsRngCG(n);
if IsVerbose("MolInv") then
 MolInvF95PrintDims(R0);
end if;
R0`SecondaryInvariants,R0`IrreducibleSecondaries:=
  MolInvSecsToDegree(mxd,R0:Names:=MolSymCGNames(n));
return R0;
end intrinsic;

intrinsic MolInvRngNewCG
 (n::RngIntElt,
  mxd::RngIntElt) ->
 RngInvar
{Invariant Ring with Secondaries to degree mxd for MolSymNewCG}
R0:=MolInvPrimsRngNewCG(n);
if IsVerbose("MolInv") then
 MolInvF95PrintDims(R0);
end if;
R0`SecondaryInvariants,R0`IrreducibleSecondaries:=
  MolInvSecsToDegree(mxd,R0:Names:=MolSymNewCGNames(n));
return R0;
end intrinsic;

intrinsic MolInvRngBG
 (m::RngIntElt,
  n::RngIntElt,
  mxd::RngIntElt) ->
 RngInvar
{Invariant Ring with Secondaries to degree mxd for MolSymBG}
R0:=MolInvPrimsRngBG(m,n);
if IsVerbose("MolInv") then
 MolInvF95PrintDims(R0);
end if;
R0`SecondaryInvariants,R0`IrreducibleSecondaries:=
  MolInvSecsToDegree(mxd,R0:Names:=MolSymBGNames(m,n));
return R0;
end intrinsic;

intrinsic MolInvRngGen
 (nki::[RngIntElt],
  mxd::RngIntElt:
  nkj:=nki) ->
 RngInvar
{Invariant Ring with Secondaries to degree mxd for MolSymGen}
R0:=MolInvPrimsRngGen(nki:nkj:=nkj);
if IsVerbose("MolInv") then
 MolInvF95PrintDims(R0);
end if;
R0`SecondaryInvariants,R0`IrreducibleSecondaries:=
  MolInvSecsToDegree(mxd,R0:
  Names:=MolSymGenNames(nki));
return R0;
end intrinsic;

intrinsic MolInvRngGenVec
 (nki::[RngIntElt],
  k::RngIntElt,
  mxd::RngIntElt) ->
 RngInvar
{Vector Covariants with Secondaries to degree mxd for MolSymGen}
R0:=MolInvRngGen(nki,mxd);
nkj:=[nki[i]:i in [1..k-1]] cat [1,nki[k]-1] cat [nki[i]:i in [k+1..#nki]];
G:=MolSymGen(nkj);
R1:=InvariantRing(G,Rationals());
R:=PolynomialRing(R1);
_:=HilbertSeries(R1);
R1`PrimaryInvariants:=R0`PrimaryInvariants;
R1`SecondaryInvariants,R1`IrreducibleSecondaries:=
  MolInvSecsToDegree(mxd,R1:
  IniPrims:=[R!f: f in R0`PrimaryInvariants] cat
    [R!f:f in R0`IrreducibleSecondaries],
  Names:=MolSymGenNames(nki));
return R1;
end intrinsic;
