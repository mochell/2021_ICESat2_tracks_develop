#!/bin/bash
# Run one analysis_db stage locally, headless, with a log file.
# usage: analysis_db/run_local.sh <script.py> <ID> <batch_key> <flag>
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PYTHONSTARTUP="$ROOT/config/pythonrc_local.py" MPLBACKEND=Agg
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
cd "$ROOT/analysis_db"
stage=$(basename "$1" .py)
log="$ROOT/logs/${stage}_$2.log"
echo "=== START $(date -Iseconds) $*" | tee "$log"
/opt/anaconda3/envs/2026-icesat2-tracks/bin/python "$@" >> "$log" 2>&1
rc=$?
echo "=== EXIT $rc $(date -Iseconds)" | tee -a "$log"
exit $rc
