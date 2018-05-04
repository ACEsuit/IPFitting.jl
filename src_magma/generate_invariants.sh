#!/bin/bash

NBODY=5
NBlengths=$((($NBODY*($NBODY-1))/2))
DEGREE=6

ECHO Nbody order= $NBODY
ECHO Nb of lengths= $NBlengths
ECHO Polynomial degree= $DEGREE

# filename_output="NBody_$NBODY""_deg_$DEGREE""_output.txt"
filename_log="NB_$NBODY""_deg_$DEGREE""_log.txt"
fn_jl_check="NB_$NBODY""_deg_$DEGREE""_non_efficient_invariants.jl"
fn_jl_irr_inv="NB_$NBODY""_deg_$DEGREE""_irr_invariants.jl"
fn_jl_prim_inv="NB_$NBODY""_deg_$DEGREE""_prim_invariants.jl"
fn_jl_sec_rel_inv="NB_$NBODY""_deg_$DEGREE""_relations_invariants.jl"

ECHO Output files:

ECHO $filename_log
ECHO $fn_jl_check
ECHO $fn_jl_irr_inv

cp Nbody_inv_auto_generation.m Nbody_run.m;

sed -i -e "s/DEGREE/$DEGREE/g" Nbody_run.m;
sed -i -e "s/NBODY/$NBODY/g" Nbody_run.m;

scp pack_opt_primaries.m dusson@galois.warwick.ac.uk: ;
scp Nbody_run.m dusson@galois.warwick.ac.uk: ;

ssh dusson@galois.warwick.ac.uk << EOF
magma Nbody_run.m
EOF

rm Nbody_run.m;

scp dusson@galois.warwick.ac.uk:magma_invariants/logNbody_output.txt .;

mv logNbody_output.txt $filename_log

# Generate julia file with function computing primary and secondary invariants (not efficient but hopefully correct)
# Pick lines with primaries, irreducible secondaries and secondaries
cp $filename_log $fn_jl_check
# remove things that are not secondaries or primaries
sed -i '' '/v\[1\]/,$!d' $fn_jl_check
sed -i '' '/Total/d' $fn_jl_check

TEMPVAR=$(grep "Nb_secondary_invariants" $fn_jl_check)
NBsecondaries=${TEMPVAR#*=}
ECHO "Nb of secondaries="$NBsecondaries

sed -i '' '/ Nb_secondary_invariants/d' $fn_jl_check

# replace variables for the primaries
# xi -> x[i]
for a in `seq $NBlengths -1 1`; do
	OLD="x$a" ;
	NEW="x[$a]" ;
	sed -i '' "s/$OLD/$NEW/g" $fn_jl_check
done

# replace variables for the secondaries
# and d(ik,il) -> x[i]
count=$NBlengths
for a in `seq $(($NBODY-2)) -1 0`; do
	for b in `seq $((($NBODY-1))) -1 $(($a+1))`; do
		ECHO $a,$b, $count
		OLD="d(i$a,i$b)" ;
		NEW="x\[$count\]" ;
		sed -i '' "s/$OLD/$NEW/g" $fn_jl_check
		count=$(($count-1))
	done
done

cp $fn_jl_check $fn_jl_irr_inv
cp $fn_jl_check $fn_jl_prim_inv
cp $fn_jl_check $fn_jl_sec_rel_inv

sed -i '' '/SYM/d' $fn_jl_check

echo "v=zeros($NBsecondaries"",1);" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "pv=zeros($NBsecondaries"",1);" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "function invariants_Q$NBlengths""_check(x)" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check


echo "return Primary_invariants, v, pv"  >> $fn_jl_check
echo "" >> $fn_jl_check
echo "end" >> $fn_jl_check
echo "x = rand($NBlengths"")"  >> $fn_jl_check
echo "display(invariants_Q$NBlengths""_check(x))"  >> $fn_jl_check

#TODO: put lines in Primary invariants as a single line (otherwise doesnt work.)
#Shortcut cmd+j then cmd+shift+j
#TODO: check that no lines start by + (otherwise wring invariants are computed)

# Generate a file with only monomials of irreducible secondaries
# ---------------------------------------------------------
sed -i '' '/SYM/!d' $fn_jl_irr_inv
sed -i '' 's/SYM/ /' $fn_jl_irr_inv

# Generate a file with only monomials of primaries
# ---------------------------------------------------------
sed -i '' '/Primary/,$!d' $fn_jl_prim_inv

# Generate a file with relations between irreducible and secondary invariants
# ---------------------------------------------------------
sed -i '' '/ v\[/!d' $fn_jl_sec_rel_inv

mv $filename_log ../data_temp/$filename_log
mv $fn_jl_check ../data_temp/$fn_jl_check
mv $fn_jl_irr_inv ../data_temp/$fn_jl_irr_inv
mv $fn_jl_prim_inv ../data_temp/$fn_jl_prim_inv
mv $fn_jl_sec_rel_inv ../data_temp/$fn_jl_sec_rel_inv

rm NBody_run.m-e

# # ----------------------------------
# cp generate_irr_secondaries.jl irr_sec_run.jl;
#
# sed -i -e "s/DEGREE/$DEGREE/g" irr_sec_run.jl;
# sed -i -e "s/NBODY/$NBODY/g" irr_sec_run.jl;
# #
# /Applications/Julia-0.6.app/Contents/Resources/julia/bin/julia irr_sec_run.jl
#
# rm irr_sec_run.jl
