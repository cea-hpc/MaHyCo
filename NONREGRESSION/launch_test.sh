#!/bin/bash
# Launch the test and compares obtained results to expected ones
set -euo pipefail 

# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local return_code=0
  $1 $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    return_code=1
  fi
  return ${return_code}
}
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_seq_pr {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local return_code=0
  $1 -arcane_opt max_iteration 10 $data_dir/Donnees.arc
  $1 -arcane_opt continue $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    return_code=1
  fi
  return ${return_code}
}
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_para_4 {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local return_code=0
  mpiexec -n 4 $1 $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    return_code=1
  fi
  return ${return_code}
}
# This function launch the computation by calling the executable with arguments
# taken from args.txt file
function launch_computation_para_8 {
  local readonly exe_path=$1
  local readonly data_dir=$2
  local return_code=0
  mpiexec -n 8 $1 $data_dir/Donnees.arc
  if [[ $? -ne 0 ]]; then
    echo "A problem occured during test execution."
    return_code=1
  fi
  return ${return_code}
}
# This function compares the results obtained in the output directory with those
# expected in the reference directory
function compare_results {
  local readonly reference_dir=$1
  local return_code=0
  ls -l "$reference_dir/output"
  if ! diff -r output/depouillement "$reference_dir/output/depouillement" > /dev/null 2>&1; then
      echo "Test failure! Output is different from reference"
      echo ${test_dir} >>  list_of_cases_to_change
      diff -r output/depouillement "$reference_dir/output/depouillement" > DIFF.txt
    return_code=1
  fi
  return ${return_code}
}

# Main function. Calls launch_computation and compare_results
# For each test a directory is created inside /tmp
function main {
  local readonly exe_path=$1
  local readonly test_dir=$2
  local readonly type=$3

  echo "Executable path : ${exe_path}"
  echo "Executing test $(basename ${test_dir})"
  echo "Executing test mode ${type}"
  if [ ${type} = "para" ]; then
    echo "Executing test parallele mode ${para}"
  fi
  
  tmp_dir=$(mktemp -d -t mahyco-ci-XXXXXXXXXX)

  if [[ ! -d ${tmp_dir} ]]; then
    echo "Unable to create a temporary directory!"
    echo "Aborting!"
    exit 3
  fi

  cd ${tmp_dir}
  echo "This directory contains the output of the test under ${test_dir}" > README.txt
  pwd
  echo ${type}
  if [ ${type}  = "para_8" ]
  then
      echo " lancement parallele sur 8 coeurs" 
      launch_computation_para_8 ${exe_path} ${test_dir}
  elif [ ${type}  = "para_4" ]
  then
      echo " lancement parallele sur 4 coeurs" 
      launch_computation_para_4 ${exe_path} ${test_dir}
  elif [ ${type}  = "seq_pr" ]
  then
      echo " lancement sequentiel protection-reprise" 
      launch_computation_seq_pr ${exe_path} ${test_dir}    
  else
      echo " lancement sequentiel" 
      launch_computation ${exe_path} ${test_dir}
  fi    
  if [[ $? -ne 0 ]]; then
    echo "Aborting!"
    exit 1
  fi

  compare_results ${test_dir}
  if [[ $? -ne 0 ]]; then
    echo "Aborting!"
    exit 2
  else
    echo "Success!"
  fi
  exit 0
}

main $@
