#!/bin/bash
# ==============================================================================
# generate_msub.sh
# Generates one MSUB file per Donnees*_nsd_*.arc data file
# TGCC environment (ccc_msub / ccc_mprun)
#
# Usage:
#   ./generate_msub.sh -e <executable> [OPTIONS]
#
# Options:
#   -e <path>          Path to executable (mandatory, resolved to absolute path)
#   -p <partition>     TGCC partition                    (default: gh200-bxi)
#   -a <args>          Extra arguments for the executable (default: none)
#   -n <tasks/node>    MPI tasks per node                (default: 288)
#   -t <seconds>       Wall-clock time limit in seconds  (default: 3600)
#   -w <dir>           Working directory                 (default: current dir)
#   -m                 Multi-thread mode: same .arc as MPI, launch with -n 1 -c <ncores>
#   -h                 Show this help
# ==============================================================================

# --- Default parameters ---
PARTITION="gh200-bxi"
EXECUTABLE=""
EXE_ARGS=""
NTASKS_PER_NODE=288
TIME_LIMIT=3600
WORKDIR=$(pwd)
MULTITHREAD=0
# --------------------------

usage() {
    sed -n '/^# Usage/,/^# ====/p' "$0" | sed 's/^# \{0,3\}//'
    exit 0
}

while getopts ":p:e:a:n:t:w:mh" opt; do
    case $opt in
        p) PARTITION="$OPTARG" ;;
        e) EXECUTABLE="$OPTARG" ;;
        a) EXE_ARGS="$OPTARG" ;;
        n) NTASKS_PER_NODE="$OPTARG" ;;
        t) TIME_LIMIT="$OPTARG" ;;
        w) WORKDIR="$OPTARG" ;;
        m) MULTITHREAD=1 ;;
        h) usage ;;
        :) echo "ERROR: option -$OPTARG requires an argument."; exit 1 ;;
        \?) echo "ERROR: unknown option -$OPTARG"; exit 1 ;;
    esac
done

# --- Validate mandatory -e option ---
if [[ -z "$EXECUTABLE" ]]; then
    echo "ERROR: option -e <path> is mandatory."
    echo "Usage: $0 -e <path/to/executable> [OPTIONS]"
    exit 1
fi

# Resolve to absolute path
EXECUTABLE=$(realpath "$EXECUTABLE" 2>/dev/null)
if [[ $? -ne 0 ]]; then
    echo "ERROR: cannot resolve executable path: $EXECUTABLE"
    exit 1
fi

if [[ ! -x "$EXECUTABLE" ]]; then
    echo "ERROR: executable not found or not executable: $EXECUTABLE"
    exit 1
fi

echo "=== generate_msub.sh ==="
echo "  Partition       : $PARTITION"
echo "  Executable      : $EXECUTABLE"
echo "  Extra args      : ${EXE_ARGS:-(none)}"
echo "  Tasks/node      : $NTASKS_PER_NODE"
echo "  Time limit      : ${TIME_LIMIT}s"
echo "  Working dir     : $WORKDIR"
if [[ $MULTITHREAD -eq 1 ]]; then
    echo "  Mode            : multi-thread (-n 1 -c <ncores>)"
else
    echo "  Mode            : pure MPI"
fi
echo ""

shopt -s nullglob
files=(Donnees*_nsd_*.arc)

if [[ ${#files[@]} -eq 0 ]]; then
    echo "ERROR: no file matching Donnees*_nsd_*.arc found in $(pwd)"
    exit 1
fi

count=0

for arc_file in "${files[@]}"; do
    basename_noext="${arc_file%.arc}"
    nsd_part=$(echo "$basename_noext" | grep -oP '\d+x\d+x\d+')

    if [[ -z "$nsd_part" ]]; then
        echo "  [SKIP] unexpected format: $arc_file"
        continue
    fi

    IFS='x' read -r nx ny nz <<< "$nsd_part"
    nprocs=$(( 10#$nx * 10#$ny * 10#$nz ))
    nprocs_fmt=$(printf "%03d" "$nprocs")
    prefix=$(echo "$basename_noext" | grep -oP '^Donnees\w*(?=_nsd_)')

    # --- Pure MPI ---
    if [[ $MULTITHREAD -eq 0 ]]; then
        if [[ -n "$EXE_ARGS" ]]; then
            run_cmd="${EXECUTABLE} ${arc_file} ${EXE_ARGS}"
        else
            run_cmd="${EXECUTABLE} ${arc_file}"
        fi
        nnodes=$(( (nprocs + NTASKS_PER_NODE - 1) / NTASKS_PER_NODE ))
        job_name="${prefix:-Donnees}_${nsd_part}"
        msub_file="job_${basename_noext}.msub"
        log_out="job_${basename_noext}.out"
        log_err="job_${basename_noext}.err"
        listing="listing_${PARTITION}_n${nprocs_fmt}"

        cat > "$msub_file" <<EOF
#!/bin/bash
#MSUB -r ${job_name}
#MSUB -o ${log_out}
#MSUB -e ${log_err}
#MSUB -q ${PARTITION}
#MSUB -N ${nnodes}
#MSUB -n ${nprocs}
#MSUB -c 1
#MSUB -T ${TIME_LIMIT}
#MSUB -x

echo "Job      : ${job_name}"
echo "File     : ${arc_file}"
echo "NSD      : ${nx}x${ny}x${nz}  =>  ${nprocs} MPI tasks on ${nnodes} node(s)"
echo "Listing  : ${listing}"
echo "Start    : \$(date)"

cd ${WORKDIR}

module load c++/gcc/12 cuda/12.4 hdf5/1.14.3 mpi/openmpi/4.1.7

ccc_mprun -n ${nprocs} -N ${nnodes} ${run_cmd} > ${listing}

echo "End      : \$(date)"
EOF
        echo "  [OK] $msub_file  ($nprocs procs, $nnodes node(s), listing: $listing)"

    # --- Multi-thread ---
    else
        if [[ -n "$EXE_ARGS" ]]; then
            run_cmd="${EXECUTABLE} -A,T=${nprocs} ${arc_file} ${EXE_ARGS}"
        else
            run_cmd="${EXECUTABLE} -A,T=${nprocs} ${arc_file}"
        fi
        job_name="${prefix:-Donnees}_mt_${nsd_part}"
        msub_file="job_${prefix:-Donnees}_mt_${nsd_part}.msub"
        log_out="job_${prefix:-Donnees}_mt_${nsd_part}.out"
        log_err="job_${prefix:-Donnees}_mt_${nsd_part}.err"
        listing="listing_${PARTITION}_k${nprocs_fmt}"

        cat > "$msub_file" <<EOF
#!/bin/bash
#MSUB -r ${job_name}
#MSUB -o ${log_out}
#MSUB -e ${log_err}
#MSUB -q ${PARTITION}
#MSUB -N 1
#MSUB -n 1
#MSUB -c ${nprocs}
#MSUB -T ${TIME_LIMIT}
#MSUB -x

echo "Job      : ${job_name}"
echo "File     : ${arc_file}"
echo "NSD      : ${nx}x${ny}x${nz}  =>  ${nprocs} threads"
echo "Listing  : ${listing}"
echo "Start    : \$(date)"

cd ${WORKDIR}

module load c++/gcc/12 cuda/12.4 hdf5/1.14.3 mpi/openmpi/4.1.7

ccc_mprun -n 1 -c ${nprocs} ${run_cmd} > ${listing}

echo "End      : \$(date)"
EOF
        echo "  [OK] $msub_file  (1 task, $nprocs threads, listing: $listing)"
    fi

    (( count++ ))
done

echo ""
echo ">>> $count MSUB file(s) generated."
