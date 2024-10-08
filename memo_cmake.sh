#!/bin/bash



#tentative pour compiler en full nvhpc mais on link qd meme avec cc
#NVCC_BIN=/opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/
#CXX=${NVCC_BIN}/nvc++
#CC=${NVCC_BIN}/nvcc
#export CXX CC
host=`hostname -d`
if [ "$host" = "c-inti.mg1.ccc.ocre.cea.fr" ]; then
    echo "Passage pour INTI"
    rm -r $OWN_CCCWORKDIR/MaHyCo/build
    mkdir $OWN_CCCWORKDIR/MaHyCo/build
    cd $OWN_CCCWORKDIR/MaHyCo/build
    cmake $HOME/workspace/MaHyCo \
         -DCMAKE_BUILD_TYPE=Release \
         -DArcane_ROOT="/ccc/home/cont001/arcaneuser/arcaneuser/products/Rhel_8__x86_64/arcane/3.12.18.0/release/" \
         #-DArcane_ROOT="$HOME/local_arcane" \
         -DCINETIQUE_SRC="$HOME/workspace/Cinetique_chgt_phase/"
    
    cmake --build . -- -j16 
    cd $OWN_CCCWORKDIR/MaHyCo/build
    
    cp $HOME/workspace/MaHyCo/NONREGRESSION/CAS_BiSodCaseX/Donnees.arc .
    mpiexec -n 1 ./src/Mahyco Donnees.arc
    rm -rf Donnees.arc output
else
    echo "Passage pour Ubuntu"
    rm -rf build
    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release -DArcane_ROOT="$HOME/Workspace/arcane/" -DCINETIQUE_SRC="$HOME/Workspace/Cinetique_chgt_phase/"
    #cmake .. -DCMAKE_BUILD_TYPE=Release -DArcane_ROOT="$HOME/Workspace/arcane/"
    #cmake .. -DCMAKE_BUILD_TYPE=Release -DArcane_ROOT="$HOME/Workspace/arcane_cuda/"
    cmake --build . -- -j16 

    # Build with Ninja (natively parallel)
    #cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DArcane_ROOT="$HOME/workspace/arcane/" -DWANT_CUDA=TRUE -DWANT_PROF_ACC=TRUE

    #cmake --build . --target test
    cp ../NONREGRESSION/CAS_BiSodCaseX/Donnees.arc .
    mpiexec -n 1 ./src/Mahyco Donnees.arc
    rm -rf Donnees.arc output
fi

