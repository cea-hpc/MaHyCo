#!/bin/bash

#set -x

VERBOSE=false
btype="release"
ARC_BUILD_TYPE="Release"
LAUNCH_ARC_CTEST=false

host=`hostname -d`
kernel=`uname -s`
archi=`uname -m`

# Print usage help
usage() {
    echo "Usage: $0 <PATH_TO_ARC_SRC> <PATH_TO_ARC_INSTALL> [-build=] [-suff=] [arc_tests] [-v] [-h]"
    echo "  <PATH_TO_ARC_SRC>     The path to your Arcane sources"
    echo "  <PATH_TO_ARC_INSTALL> The root of the path of your Arcane installation"
    echo "                        The tree structure will be :  PATH_TO_ARC_INSTALL/\$OS/\$BUILD_TYPE "
    echo "  -build                Build type (release, debug, check), release by default"
    echo "  -v                    Enable verbose mode"
    echo "  -h                    Show this help message"
    echo "  -arc_tests            Launch Arcane ctest after building (could be very long), deactivated by default"
    echo "  -acc                  Specify which kind of accelerator you want to use"
    echo "  -suff                 Suffix to be added after the Arcane version number in the <PATH_TO_ARC_INSTALL> variable"
    echo "                        For instance, a hash commit."
    echo ""
    echo "This bash script is designed for Laptop computers and clusters from the CEA."
    echo "If you want to adapt it to a new archicture, please look at usage and adapt anything that is required."
    echo "At first, the hostname of your computer should be added to this script."
    echo "Then, verify your dependecies (most likely, the same you use for building Arcane)."
    exit 1
}

if [[ -z "${1}" || -z "${2}" ]]
then
  usage
fi

ARCANE_SRC_ROOT=`realpath -q ${1}`
ARCANE_INSTALL_ROOT=`realpath -q ${2}`
shift 2

for arg in $*
do
  case ${arg} in
    -build=*)
      btype=$(echo "${arg}" | cut -d= -f2)
      if [[ "${btype}" == "release" ]]
      then
        ARC_BUILD_TYPE="Release"
      elif [[ "${btype}" == "debug" ]]
      then
        ARC_BUILD_TYPE="Debug"
      elif [[ "${btype}" == "check" ]]
      then
        ARC_BUILD_TYPE="Check"
      else
        echo "-build= : ${btype} valeur inconnue. Valeurs autorisées : {release, debug, check}"
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
    
    -arc_tests)
      LAUNCH_ARC_CTEST=true
      shift 1
    ;;
    
    -acc=*)
      ACC_MODE_SUFFIX="${arg#-acc=}"
      ACC_MODE_SUFFIX_PART="${ACC_MODE_SUFFIX:+_${ACC_MODE_SUFFIX}}"
      ARCANE_ACCELERATOR_MODE=${ACC_MODE_SUFFIX}
      shift 1
    ;;
    
    -acc=)
      echo "-acc= : nécessiste une chaine de caractère. Valeurs autorisées : {CUDA, ROCM, HIP, SYCL}"
      exit 1
      shift 1
    ;;
    
    -suff=*)
      SUFFIX="${arg#-suff=}"
      SUFFIX_PART="${SUFFIX:+_${SUFFIX}}"
      shift 1
    ;;
    
    -suff=)
      echo "-suff= : nécessiste une chaine de caractère"
      exit 1
    ;;

    *)
      echo "${arg} : argument inconnu"
      usage
    ;;
  esac
done

ARCANE_VERSION=$(<${ARCANE_SRC_ROOT}/arcane/version)
if [[ -z "${ARCANE_VERSION}" ]]
then
  echo "File arcane/version couldn't be found in ${ARCANE_SRC_ROOT}. Please check where the Arcane project is located."
  exit 1
fi
if [ ${VERBOSE} == true ]; then
  echo "Arcane version is : " ${ARCANE_VERSION}
fi


#########################
# Inti Cluster configuration 
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
  export CCCOS
  CCCARCH=`uname -m`
  PRODUCTS_ROOT_BASE=/ccc/home/cont001/arcaneuser/arcaneuser/products
  PRODUCTS_ROOT=${PRODUCTS_ROOT_BASE}/${CCCOS}
  export PRODUCTS_ROOT_BASE
  export PRODUCTS_ROOT
  echo PRODUCTS_ROOT=${PRODUCTS_ROOT}
  
  DOTNET_PATH=${PRODUCTS_ROOT_BASE}/dotnet-${CCCARCH}/8.0.17
  NINJA_PATH=${PRODUCTS_ROOT_BASE}/rhel8-${CCCARCH}/ninja/1.10.2/bin
  PATH=${DOTNET_PATH}:${NINJA_PATH}:${PATH}
  export PATH
  echo PATH=$PATH
  
  _HWLOC_PATH="/ccc/products/hwloc-2.2.0/system/default"
  if [ "${CCCOS}" == "Rhel_9__aarch64" ]; then
    echo "Configuration for GH200"
    # GH200
    _HWLOC_PATH="/ccc/products/hwloc-2.9.2/system__cuda--12.4/default"
  fi
  
  _TBB_PATH="/ccc/home/cont001/arcaneuser/arcaneuser/products/rhel8-${CCCARCH}/tbb/2021.5"
  _GOOGLETEST_PATH="/ccc/home/cont001/arcaneuser/arcaneuser/products/rhel7-${CCCARCH}/googletest/1.10.0"
  _OTF2_PATH="/ccc/products/otf2-2.2/gcc--8.3.0__openmpi--4.0.1/default"
  _PARMETIS_PATH="/ccc/home/cont001/arcaneuser/arcaneuser/products/${CCCOS}/parmetis/4.0.3-ompi405"
  
  COMMON_CMAKE_PREFIX_PATH="${_HWLOC_PATH};${_TBB_PATH};${_GOOGLETEST_PATH};${_OTF2_PATH};${_PARMETIS_PATH}"
  export COMMON_CMAKE_PREFIX_PATH
  echo "COMMON_CMAKE_PREFIX_PATH=${COMMON_CMAKE_PREFIX_PATH}"
  ARCANE_INSTALL_PREFIX=${ARCANE_INSTALL_ROOT}/arcane${ARCANE_VERSION}${SUFFIX_PART}_gcc123${ACC_MODE_SUFFIX_PART}_mpi417/${CCCOS}/${ARC_BUILD_TYPE}
  module purge
  module load cmake/3.26.4 c++/gcc/12 cuda/12.4 hdf5/1.14.3 mpi/openmpi/4.1.7
  MPI_LAUNCHER="/usr/bin/ccc_mprun"

else  

#########################
# Laptop configuration 
#########################
  CCCOS=$(echo "${kernel}"__"${archi}")
  export CCCOS
  CUDA_ARCHI=75
  ARCANE_INSTALL_PREFIX=${ARCANE_INSTALL_ROOT}/arcane${ARCANE_VERSION}${SUFFIX_PART}${ACC_MODE_SUFFIX_PART}/${CCCOS}/${ARC_BUILD_TYPE}
  MPI_LAUNCHER=`which mpiexec`
fi


if [ ${VERBOSE} == true ]; then
  echo "Arcane will be installed at : " ${ARCANE_INSTALL_PREFIX}
fi

if [ -d ${ARCANE_INSTALL_PREFIX} ]; then
  cd ${ARCANE_INSTALL_PREFIX}
  rm -rf   build
  mkdir -p build
  cd       build
else
  echo "Folder ${ARCANE_INSTALL_PREFIX} does not exists, creating it."
  mkdir -p ${ARCANE_INSTALL_PREFIX}
  cd ${ARCANE_INSTALL_PREFIX}
  rm -rf   build
  mkdir -p build
  cd       build
fi


export CXX=`which c++`
export CC=`which gcc`
export CXX CC


# Si on veut compiler en spécifiant directement le compilateur:
# -DCMAKE_CUDA_COMPILER=/ccc/products/nvhpc-21.7/system/default/Linux_x86_64/21.7/cuda/11.4/bin/nvcc

# On désactive le wrapper C# par défaut car il n'est pas utilisé sur inti.

if [ ${VERBOSE} == true ]; then
  if [ -n "${ARCANE_ACCELERATOR_MODE}"  ]; then
    echo "Cmake configuration line : "
    echo "
         cmake -DCMAKE_INSTALL_PREFIX=${ARCANE_INSTALL_PREFIX} 
  	     -DCMAKE_PREFIX_PATH=${COMMON_CMAKE_PREFIX_PATH};${HYPRE_PREFIX} 
  	     -DARCANE_ACCELERATOR_MODE=${ACC_MODE_SUFFIX} 
  	     -DARCANEFRAMEWORK_BUILD_COMPONENTS=Arcane 
  	     -DCMAKE_DISABLE_FIND_PACKAGE_SWIG=TRUE 
  	     -DARCANE_BUILD_TYPE=${ARC_BUILD_TYPE} 
  	     -DCMAKE_CUDA_ARCHITECTURES=${CUDA_ARCHI} 
  	     -DARCANE_CUSTOM_MPI_DRIVER=${MPI_LAUNCHER} 
         ${ARCANE_SRC_ROOT} "$*"
         "
  else
    echo "Cmake configuration line : "
    echo "
         cmake -DCMAKE_INSTALL_PREFIX=${ARCANE_INSTALL_PREFIX} 
  	     -DCMAKE_PREFIX_PATH=${COMMON_CMAKE_PREFIX_PATH};${HYPRE_PREFIX} 
  	     -DARCANEFRAMEWORK_BUILD_COMPONENTS=Arcane 
  	     -DCMAKE_DISABLE_FIND_PACKAGE_SWIG=TRUE 
  	     -DARCANE_BUILD_TYPE=${ARC_BUILD_TYPE} 
  	     -DARCANE_CUSTOM_MPI_DRIVER=${MPI_LAUNCHER} 
         ${ARCANE_SRC_ROOT} "$*"
         "
  fi
fi


if [ -n "${ARCANE_ACCELERATOR_MODE}"  ]; then

$( cmake -DCMAKE_INSTALL_PREFIX="${ARCANE_INSTALL_PREFIX}" \
      -DCMAKE_PREFIX_PATH="${COMMON_CMAKE_PREFIX_PATH};${HYPRE_PREFIX}" \
      -DARCANE_ACCELERATOR_MODE="${ACC_MODE_SUFFIX}" \
      -DARCANEFRAMEWORK_BUILD_COMPONENTS=Arcane \
      -DCMAKE_DISABLE_FIND_PACKAGE_SWIG=TRUE \
      -DARCANE_BUILD_TYPE="${ARC_BUILD_TYPE}" \
      -DCMAKE_CUDA_ARCHITECTURES="${CUDA_ARCHI}" \
      -DARCANE_CUSTOM_MPI_DRIVER="${MPI_LAUNCHER}" \
"${ARCANE_SRC_ROOT}" "$*")
       
else
 
$( cmake -DCMAKE_INSTALL_PREFIX="${ARCANE_INSTALL_PREFIX}" \
      -DCMAKE_PREFIX_PATH="${COMMON_CMAKE_PREFIX_PATH};${HYPRE_PREFIX}" \
      -DARCANEFRAMEWORK_BUILD_COMPONENTS=Arcane \
      -DCMAKE_DISABLE_FIND_PACKAGE_SWIG=TRUE \
      -DARCANE_BUILD_TYPE="${ARC_BUILD_TYPE}" \
      -DARCANE_CUSTOM_MPI_DRIVER="${MPI_LAUNCHER}" \
"${ARCANE_SRC_ROOT}" "$*" )
       
fi

# Build
if [ "$host" == "c-inti.mg1.ccc.ocre.cea.fr" ]; then
  cmake --build . -- -j 288
else
  cmake --build . -- -j 8
fi

# Install
cmake --build . --target install
# Tests
if [ ${LAUNCH_ARC_CTEST} == true ]; then
  cmake --build . --target test
fi

