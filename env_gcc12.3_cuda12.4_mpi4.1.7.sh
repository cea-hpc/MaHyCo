#!/bin/sh

module purge
module load cmake/3.26.4 c++/gcc/12 cuda/12.4 hdf5/1.14.3 swig mpi/openmpi/4.1.7

export CXX=`which c++`
export CC=`which gcc`
export CXX CC

# en attendant que SLURM_JOB_ID soit correctement configuré sur GH200 si allocation via ccc_mprun -x
if [[ ! -v SLURM_JOB_ID ]]; then
  JOB_ID="${CCCSHMDIR#/dev/shm/SLURM_}"
  export SLURM_JOB_ID=${JOB_ID}
fi
