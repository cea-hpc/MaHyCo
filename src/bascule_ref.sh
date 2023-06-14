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
	test1=$(grep nsd Donnees.arc | grep 4)
	echo $test1 > test
	echo " resulat test1"
	echo $test1
    test2=$(grep "2" test)
    echo " resulat test2"
    echo $test2
	if [ "$test1" != "" ]
	then
	    if [ "$test2" != "" ]
	    then
            echo " cas sur 8 procs"
            mpiexec -n 8 $curent_dir/build/src/Mahyco Donnees.arc
            RUN_PROC="cas sur 8 procs"
	    else
            echo " cas sur 4 procs"
            mpiexec -n 4 $curent_dir/build/src/Mahyco Donnees.arc
            RUN_PROC="cas sur 4 procs"
	    fi
	else
        test3=$(grep "2 2 2" Donnees.arc)
        echo " resulat test3 "
        echo $test3
        if [ "$test3" != "" ]
	    then
            echo " cas sur 8 procs"
            mpiexec -n 8 $curent_dir/build/src/Mahyco Donnees.arc
            RUN_PROC="cas sur 8 procs"
		else 
            test4=$(grep "2 2" test)
            if [ "$test4" != "" ]
            then
                echo " cas sur 4 procs"
                mpiexec -n 4 $curent_dir/build/src/Mahyco Donnees.arc
                RUN_PROC="cas sur 4 procs"
            else
                echo " cas sur 1 procs"
                $curent_dir/build/src/Mahyco Donnees.arc
                RUN_PROC="cas sur 1 proc"
            fi
		fi
	fi
	#  $curent_dir/build/src/Mahyco -arcane_opt max_iteration 10 Donnees.arc
	#  reprise
	#  $curent_dir/build/src/Mahyco -arcane_opt continue Donnees.arc
	# avec mpiexec -n 4 ou 8
	
    #paraview $cas_dir/output/depouillement/ensight.case &
	#paraview output/depouillement/ensight.case &
        #diff output $cas_dir/reference
	echo $cas
	echo $RUN_PROC
	echo "Basculer Yes/No"
	read reponse
	if  [[ "$reponse" == "Yes" ]]; then
	    rm -rf $cas_dir/output
	    mv output $cas_dir/output
	    rm -rf $cas_dir/output/listing*
	    rm -rf $cas_dir/output/checkpoint_info.xml
	    rm -rf $cas_dir/output/protection
	    rm -rf $cas_dir/output/courbes
	fi
    fi
done
}

main $@

