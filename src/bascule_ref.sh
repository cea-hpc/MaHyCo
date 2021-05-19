#!/bin/bash
#
# This script changes the reference of testcases writted on a list : list_of_cases_to_change
#
function main {
 local readonly curent_dir=$1
 local readonly test_dir=$curent_dir/NONREGRESSION
 echo "lancement"
 cat list_of_cases_to_change.txt
for cas in $(cat list_of_cases_to_change.txt); do
    local readonly cas_dir=${test_dir}/$cas
    echo CAS=$cas_dir
    if [[ -d ${cas_dir} ]]; then
	echo "lancement du $cas_dir"
        echo $curent_dir
	echo $test_dir/$cas
	cp $cas_dir/Donnees.arc .
	# $curent_dir/build/src/Mahyco Donnees.arc
	#	mpiexec -n 4 $curent_dir/build/src/Mahyco Donnees.arc
		mpiexec -n 8 $curent_dir/build/src/Mahyco Donnees.arc
	#  $curent_dir/build/src/Mahyco -arcane_opt max_iteration 10 Donnees.arc
	#  reprise
	#  $curent_dir/build/src/Mahyco -arcane_opt continue Donnees.arc
	# avec mpiexec -n 4 ou 8
	
        paraview $cas_dir/output/depouillement/ensight.case &
	#paraview output/depouillement/ensight.case &
        #diff output $cas_dir/reference
	echo $cas
	echo "Basculer Yes/No"
	read reponse
	if  [[ "$reponse" == "Yes" ]]; then
	    rm -rf $cas_dir/output
	    mv output $cas_dir/output
	    rm -rf $cas_dir/output/listing*
	fi
    fi
done
}

main $@

