#!/bin/bash

#set -x
if [[ -z $SLURM_LOCALID ]]; then
    export LOCAL_RANK=$OMPI_COMM_WORLD_LOCAL_RANK
else
    export LOCAL_RANK=$SLURM_LOCALID
fi
export CUDA_VISIBLE_DEVICES=$LOCAL_RANK

exec $*
