#!/bin/bash

NBODY=4
NBlengths=$((($NBODY*($NBODY-1))/2))
NBsecondaries=6
DEGREE=10

ECHO Nbody order= $NBODY
ECHO Nb of lengths= $NBlengths
ECHO Polynomial degree= $DEGREE

# filename_output="NBody_$NBODY""_deg_$DEGREE""_output.txt"
filename_log="NBody_$NBODY""_deg_$DEGREE""_log.txt"
fn_jl_check="NBody_$NBODY""_deg_$DEGREE""_julia_check.jl"

ECHO Output files:

ECHO $filename_log
ECHO $fn_jl_check

cp Nbody_inv_auto_generation.m Nbody_run.m;

# sed -i -e "s/DEGREE/$DEGREE/g" NBody_run.m;
# sed -i -e "s/NBODY/$NBODY/g" NBody_run.m;
#
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


echo "v=zeros($NBsecondaries"",1);" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "pv=zeros($NBsecondaries"",1);" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check
echo "function invariants_Q$NBlengths""_check(x)" | cat - $fn_jl_check > /tmp/tempfile && mv /tmp/tempfile $fn_jl_check


echo "return Primary_invariants, v"  >> $fn_jl_check
echo "" >> $fn_jl_check
echo "end" >> $fn_jl_check
echo "x = rand($NBlengths"")"  >> $fn_jl_check
echo "display(invariants_Q6_check(x))"  >> $fn_jl_check

# # sed '/Primary_invariants/{n;N;N;N;N;d}' $fn_jl_check
# sed -i '' '/Primary/1/{N;d;}' $fn_jl_check
#
# sed s/day/night/ <old >new
# # sed -i $'s/\r//' $fn_jl_check
# # sed -i '' 's/ /\'$'\n/g' $fn_jl_check
# # sed -i '' '1s/^/<added text> \n/' '\n' $fn_jl_check
# # sed -i '' '1i\' $fn_jl_check
#
# # sed -i '' '1s/^/Secondary_invariants /' $fn_jl_check
