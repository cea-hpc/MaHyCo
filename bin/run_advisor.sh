#!/bin/bash

##########
# LAPTOP #
##########
source /opt/advisor_2021_609709/advixe-vars.sh

cd /home/meltzb/workspace/MaHyCo/build/src
cp ../../NONREGRESSION/CAS_BiSodCaseX/Donnees.arc .

PROG="./Mahyco Donnees.arc"

PROJ_DIR="../../advisor_Mahyco"

SURV_FLAGS="-stackwalk-mode=offline -static-instruction-mix --track-memory-objects --data-transfer-histogram "
MAP_FLAGS="-stackwalk-mode=offline -static-instruction-mix "
TRIP_FLAGS="-flop -enable-cache-simulation -cache-sources -stacks -trip-counts "

advixe-cl --collect=survey     --project-dir=${PROJ_DIR} ${SURV_FLAGS} -- ${PROG}
advixe-cl --collect=map        --project-dir=${PROJ_DIR} ${MAP_FLAGS}  -- ${PROG}
advixe-cl --collect=tripcounts --project-dir=${PROJ_DIR} ${TRIP_FLAGS} -- ${PROG}


##########
#  INTI  #
##########
#module load advisor/21.4.0
source /ccc/home/cont001/ocre/meltzb/home_work/opt/advisor_2021_609709/advixe-vars.sh

cd /ccc/home/cont001/ocre/meltzb/home_work/cea-hpc/MaHyCo/build/src
cp ../../NONREGRESSION/CAS_BiSodCaseX/Donnees.arc .

PROG="./Mahyco Donnees.arc"

PROJ_DIR="../../advisor_Mahyco_advisor_offload"

SURV_FLAGS="-stackwalk-mode=offline -static-instruction-mix "
TRIP_FLAGS="-flop -enable-cache-simulation -cache-sources -stacks -trip-counts --cache-config=4:8w:32k:64l/4:4w:256k:64l/1:16w:6m:64l "
echo 'TRIP_FLAGS = ' ${TRIP_FLAGS}
#advixe-cl --collect=survey     --project-dir=${PROJ_DIR} ${SURV_FLAGS} -- ${PROG}
advixe-cl --collect=tripcounts --project-dir=${PROJ_DIR} ${TRIP_FLAGS} -- ${PROG}
