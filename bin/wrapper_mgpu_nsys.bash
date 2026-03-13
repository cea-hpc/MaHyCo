#!/bin/bash

#set -x
if [[ -z $SLURM_LOCALID ]]; then
    export LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
else
    export LOCAL_RANK=$SLURM_LOCALID
fi
export CUDA_VISIBLE_DEVICES=$LOCAL_RANK

# l'option --capture-range=cudaProfilerApi de nsys pose pb avec Mahyco, mais pas avec pattern4gpu, voir option nsys profile
#exec nsys profile --trace=cuda,nvtx,osrt,mpi --capture-range=cudaProfilerApi --force-overwrite true -o output/mahyco_${SLURM_LOCALID}  $*
exec nsys profile --trace=cuda,nvtx,osrt,mpi --force-overwrite true -o output/mahyco_${SLURM_LOCALID}  $*
