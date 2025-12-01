#!/bin/bash
#
# This script changes the reference of testcases writted on a list : list_of_cases_to_change
# Lancement :
# cd mahyco
# un fichier list_of_cases_to_change a été créé par ctest au préalable
# ./NONREGRESSION/bascule_ref.sh PATH_RACINE_MAHYCO PATH_DU_BUILD 
#

# -----------------------------------------------------
# Aide
# -----------------------------------------------------
function helpme {
  echo "Procédure de lancement de la mise à jour des cas de non régression"
  echo "Depuis la racine, lancer :"
  echo "./NONREGRESSION/bascule_ref.sh . build_XXX"
  echo ""
  echo "Ce script suppose que la compilation est faite dans le répertoire build_XXX/"
  echo "et que la non régression a tourné avec ctest pour produire le fichier list_of_cases_to_change"

  echo "Pilotage par variable d'environnement :"
  echo "AFFICHE_DIFF : affiche le diff lorsqu'il y a des écarts"
  echo "OUVRE_PARAVIEW : ouvre paraview pour visualiser les écarts"
  echo "PLOT_TIME_HISTORY_DIFF : trace les sorties bilan de time history (exécution et référence)"
  echo "BASCULE_FORCEE : accepte les changements de résultats sans poser la question, incompatible avec AFFICHE_DIFF et OUVRE_PARAVIEW"
}

# -----------------------------------------------------------------------------
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local readonly mpi_launcher=$3
  local return_code=0

  ${mpi_launcher} -n 1 $1 $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    echo "seq-"$(basename ${data_dir}) >>  $data_dir/../../list_of_pb_exec
    return_code=1
  fi
  return ${return_code}
}
# -----------------------------------------------------------------------------
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_seq_pr {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local readonly mpi_launcher=$3
  local return_code=0

  ${mpi_launcher} -n 1 $1 -arcane_opt max_iteration 10 $data_dir/Donnees.arc
  ${mpi_launcher} -n 1 $1 -arcane_opt continue $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    echo "seq_pr-"$(basename ${data_dir}) >>  $data_dir/../../list_of_pb_exec
    return_code=1
  fi
  return ${return_code}
}
# -----------------------------------------------------------------------------
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_para_4 {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local readonly mpi_launcher=$3
  local return_code=0

  ${mpi_launcher} -n 4 $1 $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    echo "para_4-"$(basename ${data_dir}) >>  $data_dir/../../list_of_pb_exec
    return_code=1
  fi
  return ${return_code}
}
# -----------------------------------------------------------------------------
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_para_8 {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local readonly mpi_launcher=$3
  local return_code=0

  ${mpi_launcher} -n 8 $1 $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    echo "para_8-"$(basename ${data_dir}) >>  $data_dir/../../list_of_pb_exec
    return_code=1
  fi
  return ${return_code}
}
# -----------------------------------------------------------------------------
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_cuda_1 {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local readonly mpi_launcher=$3
  local return_code=0
  
  ${mpi_launcher} -n 1 $1 -A,AcceleratorRuntime=cuda $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    echo "cuda_1-"$(basename ${data_dir}) >>  $data_dir/../../list_of_pb_exec
    return_code=1
  fi
  return ${return_code}
}
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_cuda_4 {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local readonly mpi_launcher=$3
  local return_code=0

  ${mpi_launcher} -n 4 $1 -A,AcceleratorRuntime=cuda $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    echo "cuda_4-"$(basename ${data_dir}) >>  $data_dir/../../list_of_pb_exec
    return_code=1
  fi
  return ${return_code}
}

# -----------------------------------------------------
# Lancement de la procédure
# $1 : chemin vers la racine mahyco
# $2 : chemin vers le dossier de build de Mahyco
#
# Pilotage par variable d'environnement :
# AFFICHE_DIFF : affiche le diff lorsqu'il y a des écarts
# OUVRE_PARAVIEW : ouvre paraview pour visualiser les écarts
# BASCULE_FORCEE : accepte les changements de résultats sans poser la question
# -----------------------------------------------------
function main {
  
  local readonly mahyco_root_dir=$1
  local readonly test_dir=$1/NONREGRESSION
  local readonly mahyco_build_dir=$2

  if [ $mahyco_root_dir == "-h" ] ; then 
    helpme
    return 0
  fi

  # Vérification que le répertoire test_dir existe bien
  # ie que le script est lancé depuis la racine de Mahyco
  if [ ! -d $test_dir ] ; then
    echo "Impossible de trouver le chemin vers le répertoire $test_dir"
    helpme
    return -1
  fi

  # Utile uniquement s'il y a des cas tests à mettre à jour identifiés par ctest
  if [ ! -e list_of_cases_to_change ]; then
    echo "Le fichier list_of_cases_to_change n'existe pas."
    echo "Lancer ctest dans le répertoire build pour créer ce fichier."
    return 0
  fi

  host=`hostname -d`
  if [ "$host" == "c-inti.mg1.ccc.ocre.cea.fr" ]; then
    mpi_launcher="/usr/bin/ccc_mprun -E --exclusive"
  else  
    mpi_launcher=`which mpiexec`
  fi

  exe_path=${mahyco_build_dir}"/src/Mahyco"

  echo "============================================================="
  echo "LANCEMENT DE LA PROCEDURE DE MISE A JOUR DES RESULTATS"
  echo "============================================================="
  echo "- Chemin vers le répertoire des cas de non reg : $test_dir"
  for cas in $(cat list_of_cases_to_change); do
    local readonly type="${cas%%-*}"
    local readonly cas_name="${cas#*-}"
    local readonly cas_dir=${test_dir}/$cas_name
    echo CAS=$cas_dir
    if [[ -d ${cas_dir} ]]; then
      echo "lancement du $cas_dir en mode ${type}"
      
      if [ ${type}  = "para_8" ]
      then
          echo " lancement parallele sur 8 coeurs" 
          launch_computation_para_8 ${exe_path} ${cas_dir} ${mpi_launcher}
      elif [ ${type}  = "para_4" ]
      then
          echo " lancement parallele sur 4 coeurs" 
          launch_computation_para_4 ${exe_path} ${cas_dir} ${mpi_launcher}
      elif [ ${type}  = "seq_pr" ]
      then
          echo " lancement sequentiel protection-reprise" 
          launch_computation_seq_pr ${exe_path} ${cas_dir} ${mpi_launcher}
      elif [ ${type}  = "cuda_1" ]
      then
          echo " lancement sequentiel 1 GPU" 
          launch_computation_cuda_1 ${exe_path} ${cas_dir} ${mpi_launcher}
      elif [ ${type}  = "cuda_4" ]
      then
          echo " lancement parallele 4 GPUs" 
          launch_computation_cuda_4 ${exe_path} ${cas_dir} ${mpi_launcher}
      else
          echo " lancement sequentiel" 
          launch_computation ${exe_path} ${cas_dir} ${mpi_launcher}
      fi    
      if [[ $? -ne 0 ]]; then
        echo "Aborting!"
        exit 1
      fi

      # Visualisation des résultats :
      if [ $OUVRE_PARAVIEW ]; then
        echo "----> Lancement de paraview pour le cas "
        rm _execution _reference # si jamais ils existent déjà, on supprime les liens symboliques
        ln -s $cas_dir/output_${type}/depouillement _reference
        ln -s output/depouillement _execution
        paraview _execution/ensight.case &
        paraview _reference/ensight.case
        rm _execution _reference
      fi

      #if [ $PLOT_TIME_HISTORY_DIFF ] ; then
      #  echo "----> Tracé du time_history output [cette exécution] et $cas_dir/output [référence]"
      #  module load python3
      #  for f in $(ls $cas_dir/output/courbes/gnuplot); do
      #    python3 $mahyco_root_dir/utils/plot_time_history.py output/courbes/gnuplot/$f --reference $cas_dir/output/courbes/gnuplot/$f  &
      #  done
      #fi

      if [ $AFFICHE_DIFF ]; then 
        echo "----> Affichage du diff output [cette exécution] et $cas_dir/output [référence]"
        diff -r output $cas_dir/output_${type}
      fi
      
      echo "----> Changement des résultats sur le cas $cas, exécution ${type}"
      if [ $BASCULE_FORCEE ]; then
        reponse="Yes"  # pour mettre à jour les références de façon systématique
      else
        echo "----> Basculer Yes/No ?"
        read reponse  
      fi

      if  [[ "$reponse" == "Yes" ]]; then
        rm -rf $cas_dir/output_${type}
        # Sauvegarde des résultats
        mv output $cas_dir/output_${type}
        # Suppression des éléments inutiles
        rm -rf $cas_dir/output_${type}/listing*
        rm -rf $cas_dir/output_${type}/checkpoint_info.xml
        rm -rf $cas_dir/output_${type}/protection
        rm -rf $cas_dir/output_${type}/courbes/curves.acv
        rm -rf $cas_dir/output_${type}/courbes/gnuplot/*Time
        rm -rf $cas_dir/output_${type}/courbes/gnuplot/TotalMemory
        rm -rf $cas_dir/output_${type}/courbes/time_history.json
        rm -rf $cas_dir/output_${type}/courbes/time_history.xml
      fi
    fi
  done
}

main $@

