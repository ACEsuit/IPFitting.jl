#!/bin/bash

NBODY=4
DEGREE=10

ECHO Nbody order= $NBODY
ECHO Polynomial degree= $DEGREE

filename_output="NBody_$NBODY""_deg_$DEGREE""_output.txt"
filename_log="NBody_$NBODY""_deg_$DEGREE""_log.txt"

ECHO Output files:

ECHO $filename_output
ECHO $filename_log

cp Nbody_inv_auto_generation.m Nbody_run.m;

sed -i -e "s/DEGREE/$DEGREE/g" NBody_run.m;
sed -i -e "s/NBODY/$NBODY/g" NBody_run.m;

scp pack_opt_primaries.m dusson@galois.warwick.ac.uk:magma_invariants;
scp Nbody_run.m dusson@galois.warwick.ac.uk:magma_invariants;

ssh dusson@galois.warwick.ac.uk << EOF
cd magma_invariants
magma Nbody_run.m
EOF

scp dusson@galois.warwick.ac.uk:magma_invariants/Nbody_output.txt .;
scp dusson@galois.warwick.ac.uk:magma_invariants/logNbody_output.txt .;

mv Nbody_output.txt $filename_output
mv logNbody_output.txt $filename_log


# declare nb of variables
# Nb_variables=14

# for a in `seq $Nb_variables -1 1`; do
#     echo "$a/10 to Exit." ;
# 	OLD="x$a" ;
# 	echo $OLD ;
# 	NEW="x[$a]" ;
# 	echo $NEW
# 	sed -i '.bak' "s/$OLD/$NEW/g" 'Nbody_5_degree_7.txt'
# done
