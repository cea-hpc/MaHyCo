#!/bin/bash


rm -rf build
mkdir build
cd build

#tentative pour compiler en full nvhpc mais on link qd meme avec cc
#NVCC_BIN=/opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/
#CXX=${NVCC_BIN}/nvc++
#CC=${NVCC_BIN}/nvcc
#export CXX CC

# Build with makefiles in parallel
cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DArcane_ROOT="/home/spraulm/arcane_install/"
ninja

# Build with Ninja (natively parallel)
#cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DArcane_ROOT="/home/meltzb/workspace/arcane/" -DWANT_CUDA=TRUE -DWANT_PROF_ACC=TRUE

cmake --build . --target test

#cd /home/spraulm/dev/mahyco/build/src
#cp ../../NONREGRESSION/CAS_BiSodCaseX/Donnees.arc .
#mpiexec -n 1 ./Mahyco Donnees.arc
