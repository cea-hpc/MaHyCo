#!/bin/bash


#tentative pour compiler en full nvhpc mais on link qd meme avec cc
#NVCC_BIN=/opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/
#CXX=${NVCC_BIN}/nvc++
#CC=${NVCC_BIN}/nvcc
#export CXX CC


host=`hostname -d`
if [ "$host" = "c-inti.mg1.ccc.ocre.cea.fr" ]; then
    echo "Passage pour INTI"
    # Build with makefiles in parallel
    cmake .. -DCMAKE_BUILD_TYPE=Release -DArcane_ROOT="/home/spraulm/arcane_install/"
    cmake --build . -- -j16 
else
    echo "Passage pour Ubuntu"
    # Build with makefiles in parallel
    cmake .. -DCMAKE_BUILD_TYPE=Release -DArcane_ROOT="/home/spraulm/arcane_install/"
    cmake --build . -- -j8
fi



# cmake --build . --target test

cp ../NONREGRESSION/CAS_BiSodCaseX/Donnees.arc .
mpiexec -n 1 ./Mahyco Donnees.arc
