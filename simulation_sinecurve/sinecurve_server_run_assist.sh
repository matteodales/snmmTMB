#!/bin/bash
# ---------------------------------------------------------------------------
# sinecurve_server_run_assist.sh
#
# Launches the sine-curve assist simulation as NJOBS independent Rscript
# processes, one per detached `screen` session.
#
#   ./sinecurve_server_run_assist.sh               # NJOBS jobs, all of them
#   ./sinecurve_server_run_assist.sh 40            # 40 jobs, all of them
#   ./sinecurve_server_run_assist.sh 40 3 8 14     # 40 jobs, launch only 3, 8, 14
#
# Run from this directory. Each job writes its own results/ and logs/ file;
# sinecurve_server_merge_and_plots.R merges the results afterwards.
#
# Linux only: this needs `screen` and `taskset`. On any other platform run
# the jobs directly instead, one at a time or backgrounded:
#   Rscript sinecurve_server_simulation_assist.R <job> <njobs>
#
# Before launching, the assist package must be installed:
#   Rscript -e 'install.packages("assist", repos = "https://cloud.r-project.org")'
# ---------------------------------------------------------------------------

set -u

# --------------------------- configuration ---------------------------------

NJOBS=20                                              # default number of jobs
SCRIPT="sinecurve_server_simulation_assist.R"
PREFIX="sinecurve_assist"                             # screen-session name prefix
CPU_OFFSET=140                                        # job i is pinned to core CPU_OFFSET+i; the three
                                                      # methods use disjoint ranges so they can run side by side
USE_TASKSET=1                                         # 0 to let the scheduler decide
NICE_LEVEL=19
LOG_DIR="logs"
RESULTS_DIR="results"

# ------------------------- command-line overrides --------------------------

JOBS_TO_RUN=""
if [ $# -ge 1 ]; then NJOBS="$1"; shift; fi
if [ $# -ge 1 ]; then JOBS_TO_RUN="$*"; fi
if [ -z "$JOBS_TO_RUN" ]; then JOBS_TO_RUN=$(seq 1 "$NJOBS"); fi
FIRST_JOB=$(echo $JOBS_TO_RUN | cut -d' ' -f1)

# --------------------------- core-range check ------------------------------

if [ "$USE_TASKSET" -eq 1 ]; then
  NCPU=$(nproc)
  MAX_JOB=$(echo $JOBS_TO_RUN | tr ' ' '\n' | sort -n | tail -1)
  MAX_CORE=$((CPU_OFFSET + MAX_JOB))
  if [ "$MAX_CORE" -ge "$NCPU" ]; then
    echo "ERROR: job ${MAX_JOB} would be pinned to core ${MAX_CORE}, but this machine has" >&2
    echo "       only ${NCPU} cores (0-$((NCPU - 1)))." >&2
    echo "       Lower CPU_OFFSET to at most $((NCPU - 1 - MAX_JOB)), or set USE_TASKSET=0." >&2
    exit 1
  fi
fi

mkdir -p "$LOG_DIR" "$RESULTS_DIR"

if [ ! -f "$SCRIPT" ]; then
  echo "ERROR: $SCRIPT not found in $(pwd)" >&2
  exit 1
fi

# ------------------------------- launch ------------------------------------

launched=0
for i in $JOBS_TO_RUN; do

  screen_name="${PREFIX}_${i}"
  log_file="${LOG_DIR}/${PREFIX}_job${i}_of${NJOBS}.log"

  # two copies of a job would write the same results file
  if screen -list | grep -q "[.]${screen_name}[[:space:]]"; then
    echo "SKIP  ${screen_name} - a screen with that name already exists"
    continue
  fi

  if [ "$USE_TASKSET" -eq 1 ]; then
    pin="taskset -c $((CPU_OFFSET + i))"
  else
    pin=""
  fi

  job_cmd="OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 nice -n ${NICE_LEVEL} ${pin} Rscript ${SCRIPT} ${i} ${NJOBS} >> ${log_file} 2>&1"

  screen -dmS "${screen_name}"
  screen -S "${screen_name}" -p 0 -X stuff "${job_cmd}\n"

  echo "START ${screen_name}  ->  ${log_file}"
  launched=$((launched + 1))
done

echo
echo "${launched} job(s) started (NJOBS=${NJOBS}, script=${SCRIPT})."
echo "Watch one:      tail -f ${LOG_DIR}/${PREFIX}_job${FIRST_JOB}_of${NJOBS}.log"
echo "List screens:   screen -ls"
echo "Attach:         screen -r ${PREFIX}_${FIRST_JOB}     (detach again with Ctrl-a d)"
echo "Count done:     grep -c '^DONE' ${LOG_DIR}/${PREFIX}_job*.log"
echo "Kill all:       for s in \$(screen -ls | grep -o '${PREFIX}_[0-9]*'); do screen -X -S \$s quit; done"
