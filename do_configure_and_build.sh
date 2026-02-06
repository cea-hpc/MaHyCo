#!/bin/sh

#set -x 

VERBOSE=false
btype="release"
MHC_BUILD_TYPE="Release"
NVTX="TRUE"

# Print usage help
usage() {
    echo "Usage: $0 <PATH_TO_MHC_SRC> <PATH_TO_ARC_INSTALL> [-build=] [-v] [-h]"
    echo "  <PATH_TO_MHC_SRC>     The absolute path to your Mahyco sources"
    echo "  <PATH_TO_ARC_INSTALL> The root of the absolute path of your Arcane installation"
    echo "  -build=               Build type (release, debug, check), release by default"
    echo "  -nvtx=                Activate NVTX profiling, on by default"
    echo "  -v                    Enable verbose mode"
    echo "  -h                    Show this help message"
    exit 1
}

MAHYCO_SRC_ROOT=`realpath -q ${1}`
ARCANE_INSTALL_PATH=`realpath -q ${2}`

if [[ -z "${MAHYCO_SRC_ROOT}" || -z "${ARCANE_INSTALL_PATH}" ]]
then
  usage
fi
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

    -v)
      VERBOSE=true
      CMAKE_VERBOSE_MAKEFILE=true
      shift 1
    ;;
    
    -h)
      shift 1
      usage
    ;;

    *)
      echo "${arg} : argument inconnu"
      usage
    ;;
  esac
done

host=`hostname -d`
kernel=`uname -s`
archi=`uname -m`

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
else  
  CCCOS=$(echo "${kernel}"__"${archi}")
  CUDA_ARCHI=75
  MPI_LAUNCHER=`which mpiexec`
  BASETMPDIR_NONREG="/tmp"
fi

rm -rf   build_${CCCOS}/${MHC_BUILD_TYPE}
mkdir -p build_${CCCOS}/${MHC_BUILD_TYPE}
cd       build_${CCCOS}/${MHC_BUILD_TYPE}


if [ ${VERBOSE} == true ]; then
  echo "Cmake configuration line : "
  echo "
	cmake  -DWANT_CUDA=TRUE \
	       -DCMAKE_BUILD_TYPE=${MHC_BUILD_TYPE} \
	       -DWANT_PROF_ACC=${NVTX} \
	       -DArcane_ROOT=${ARCANE_INSTALL_PATH}/${CCCOS}/${MHC_BUILD_TYPE} \
	       -DMPI_LAUNCHER=${MPI_LAUNCHER} \
	       -DBASETMPDIR_NONREG=${BASETMPDIR_NONREG} \
	       -DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCHI} \
	       ${MAHYCO_SRC_ROOT} "$*"
       "
fi

cmake  -DWANT_CUDA=TRUE \
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
cp ${MAHYCO_SRC_ROOT}/NONREGRESSION/CAS_BiSodCaseX_RemapArcane/Donnees.arc .

if [ "$host" == "c-inti.mg1.ccc.ocre.cea.fr" ]; then
  source ${MAHYCO_SRC_ROOT}/env_gcc12.3_cuda12.4_mpi4.1.7.sh
  ccc_mprun -n 4 -c 72 -p gh200-bxi ../../../bin/wrapper_mgpu.bash ./Mahyco -A,AcceleratorRuntime=cuda Donnees.arc
else
  mpiexec -n 4 ./Mahyco -A,AcceleratorRuntime=cuda Donnees.arc
fi

