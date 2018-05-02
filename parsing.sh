#!/bin/bash
# declare nb of variables
Nb_variables=14

for a in `seq $Nb_variables -1 1`; do
    echo "$a/10 to Exit." ;
	OLD="x$a" ; 
	echo $OLD ; 
	NEW="x[$a]" ;
	echo $NEW
	sed -i '.bak' "s/$OLD/$NEW/g" 'Nbody_5_degree_7.txt'
done