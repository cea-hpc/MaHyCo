#!/bin/bash

#set -x 

VERBOSE=false
btype="release"
MHC_BUILD_TYPE="Release"
NVTX="TRUE"
CUDA_MODE="FALSE"

host=`hostname -d`
kernel=`uname -s`
archi=`uname -m`

# Print usage help
usage() {
    echo "Usage: $0 <PATH_TO_MHC_SRC> <PATH_TO_ARC_INSTALL> [-build=] [-v] [-h]"
    echo "  <PATH_TO_MHC_SRC>     The absolute path to your Mahyco sources"
    echo "  <PATH_TO_ARC_INSTALL> The root of the absolute path of your Arcane installation"
    echo "  -build=               Build type (release, debug, check), release by default"
    echo "  -nvtx=                Activate NVTX profiling, on by default"
    echo "  -acc=                 Want an installation of Arcane that support GPU accelerators"
    echo "  -v                    Enable verbose mode"
    echo "  -h                    Show this help message"
    exit 1
}

if [[ -z "${1}" || -z "${2}" ]]
then
  usage
fi

MAHYCO_SRC_ROOT=`realpath -q ${1}`
ARCANE_INSTALL_PATH=`realpath -q ${2}`
shift 2


for arg in $*
do
  case ${arg} in
    -build=*)
      btype=$(echo "${arg}" | cut -d= -f2)
      if [[ "${btype}" == "release" ]]
      then
        MHC_BUILD_TYPE="Release"
      elif [[ "${btype}" == "debug" ]]
      then
        MHC_BUILD_TYPE="Debug"
      elif [[ "${btype}" == "check" ]]
      then
        MHC_BUILD_TYPE="Check"
      else
        echo "-build= : ${btype} valeur inconnue. Valeurs autorisées : {release, debug, check}"
        exit 1
      fi
      shift 1
    ;;
    
    -nvtx=*)
      nvtx=$(echo "${arg}" | cut -d= -f2)
      if [[ "${vnvtx}" == "on" ]]
      then
        NVTX="TRUE"
      elif [[ "${vnvtx}" == "off" ]]
      then
        NVTX="FALSE"
      else
        echo "-nvtx= : ${vnvtx} valeur inconnue. Valeurs autorisées : {on, off}"
        exit 1
      fi
      shift 1
    ;;

    -v|--verbose)
      VERBOSE=true
      CMAKE_VERBOSE_MAKEFILE=true
      shift 1
    ;;
    
    -h|--help)
      shift 1
      usage
    ;; 
    
    -acc=*)
      MHC_MODE_SUFFIX="${arg#-acc=}"
      MHC_MODE_SUFFIX_PART="${MHC_MODE_SUFFIX:+_${MHC_MODE_SUFFIX}}"
      CUDA_MODE="TRUE"
      shift 1
    ;;
    
    -acc=)
      echo "-acc= : nécessiste une chaine de caractère : {CUDA, HIP}"
      exit 1
    ;;

    *)
      echo "${arg} : argument inconnu"
      usage
    ;;
  esac
done

#########################
# Cluster configuration 
#########################
if [ "$host" == "c-inti.mg1.ccc.ocre.cea.fr" ]; then
  CCCOS=`/ccc/products/ccc_users_env/bin/ccc_os`
  if [ "${CCCOS}" == "Rhel_9__aarch64" ]; then
    CUDA_ARCHI=90
  elif [ "${CCCOS}" == "Rhel_8__x86_64" ]; then
    CUDA_ARCHI=80
  else
    echo "Architecture non reconnue "
    exit
  fi
  MPI_LAUNCHER="/usr/bin/ccc_mprun"
  BASETMPDIR_NONREG=${CCCSCRATCHDIR}
  source ${MAHYCO_SRC_ROOT}/env_gcc12.3_cuda12.4_mpi4.1.7.sh
#########################
# Laptop configuration 
#########################
else  
  CCCOS=$(echo "${kernel}"__"${archi}")
  CUDA_ARCHI=75
  MPI_LAUNCHER=`which mpiexec`
  BASETMPDIR_NONREG="/tmp"
fi


if [[ ! -f "${ARCANE_INSTALL_PATH}/${CCCOS}/${MHC_BUILD_TYPE}/lib/libarcane_accelerator_cuda_runtime.so" ]]; then
  echo "Lib libarcane_accelerator_cuda_runtime.so couldn't be found in ${ARCANE_INSTALL_PATH}. You will not be able to launch Mahyco with support of Nvidia GPUs."
  echo "Please recompile Arcane with support of Nvidia GPUs."
  exit 1
fi

rm -rf   build_${CCCOS}${MHC_MODE_SUFFIX_PART}/${MHC_BUILD_TYPE}
mkdir -p build_${CCCOS}${MHC_MODE_SUFFIX_PART}/${MHC_BUILD_TYPE}
cd       build_${CCCOS}${MHC_MODE_SUFFIX_PART}/${MHC_BUILD_TYPE}


if [ ${VERBOSE} == true ]; then
  echo "Cmake configuration line : "
  echo " "
  echo "cmake  -DWANT_CUDA=${CUDA_MODE} "
  echo "       -DCMAKE_BUILD_TYPE=${MHC_BUILD_TYPE} "
  echo "       -DWANT_PROF_ACC=${NVTX} "
  echo "       -DArcane_ROOT=${ARCANE_INSTALL_PATH}/${CCCOS}/${MHC_BUILD_TYPE} "
  echo "       -DMPI_LAUNCHER=${MPI_LAUNCHER} "
  echo "       -DBASETMPDIR_NONREG=${BASETMPDIR_NONREG} "
  echo "       -DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCHI} "
  echo "       ${MAHYCO_SRC_ROOT} "$*" "
fi

cmake  -DWANT_CUDA=${CUDA_MODE} \
       -DCMAKE_BUILD_TYPE=${MHC_BUILD_TYPE} \
       -DWANT_PROF_ACC=${NVTX} \
       -DArcane_ROOT=${ARCANE_INSTALL_PATH}/${CCCOS}/${MHC_BUILD_TYPE} \
       -DMPI_LAUNCHER=${MPI_LAUNCHER} \
       -DBASETMPDIR_NONREG=${BASETMPDIR_NONREG} \
       -DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCHI} \
       -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
       ${MAHYCO_SRC_ROOT} "$*"

# Build
if [ "$host" == "c-inti.mg1.ccc.ocre.cea.fr" ]; then
  cmake --build . -- -j 128
else
  cmake --build . -- -j 8
fi


# TEST 
cd src
cp ${MAHYCO_SRC_ROOT}/NONREGRESSION/CAS_BiSodCaseX/Donnees.arc .

if [ "$host" == "c-inti.mg1.ccc.ocre.cea.fr" ]; then
  source ${MAHYCO_SRC_ROOT}/env_gcc12.3_cuda12.4_mpi4.1.7.sh
  ccc_mprun -n 4 -c 72 -p gh200-bxi ../../../bin/wrapper_mgpu.bash ./Mahyco -A,AcceleratorRuntime=cuda Donnees.arc
else
  mpiexec -n 1 ./Mahyco -A,AcceleratorRuntime=cuda Donnees.arc
fi

