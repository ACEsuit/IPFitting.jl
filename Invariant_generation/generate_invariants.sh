#!/bin/bash

NBODY=5
NBlengths=$((($NBODY*($NBODY-1))/2))
DEGREE=6

ECHO Nbody order= $NBODY
ECHO Nb of lengths= $NBlengths
ECHO Polynomial degree= $DEGREE

# filename_output="NBody_$NBODY""_deg_$DEGREE""_output.txt"
filename_log="NBody_$NBODY""_deg_$DEGREE""_log.txt"
fn_jl_check="NBody_$NBODY""_deg_$DEGREE""_julia_check.jl"
fn_jl_inv="NBody_$NBODY""_deg_$DEGREE""_invariants.jl"

ECHO Output files:

ECHO $filename_log
ECHO $fn_jl_check
ECHO $fn_jl_inv

cp Nbody_inv_auto_generation.m Nbody_run.m;

sed -i -e "s/DEGREE/$DEGREE/g" NBody_run.m;
sed -i -e "s/NBODY/$NBODY/g" NBody_run.m;

# scp pack_opt_primaries.m dusson@galois.warwick.ac.uk:magma_invariants;
# scp Nbody_run.m dusson@galois.warwick.ac.uk:magma_invariants;
#
# ssh dusson@galois.warwick.ac.uk << EOF
# cd magma_invariants
# magma Nbody_run.m
# EOF
#
# scp dusson@galois.warwick.ac.uk:magma_invariants/logNbody_output.txt .;
#
# mv logNbody_output.txt $filename_log

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

cp $fn_jl_check $fn_jl_inv

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

# ---------------------------------------------------------
sed -i '' '/SYM/!d' $fn_jl_inv
sed -i '' 's/SYM/ /' $fn_jl_inv

# ----------------------------------
cp generate_irr_secondaries.jl irr_sec_run.jl;

sed -i -e "s/DEGREE/$DEGREE/g" irr_sec_run.jl;
sed -i -e "s/NBODY/$NBODY/g" irr_sec_run.jl;
#
/Applications/Julia-0.6.app/Contents/Resources/julia/bin/julia irr_sec_run.jl
