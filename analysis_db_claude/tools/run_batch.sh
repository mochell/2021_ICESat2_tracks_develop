#!/bin/bash
# Run the Snakemake pipeline for one batch with the standard resource caps.
# usage: tools/run_batch.sh <batch_key> [extra snakemake args...]
#   tools/run_batch.sh SH_dev_small                    # everything up to C01 + index
#   tools/run_batch.sh SH_dev_small --until B04
#   tools/run_batch.sh SH_dev_small -n        # dry run
# Environment: JOBS (heavy jobs, default 4), CORES (default 16),
#              SL_JOBS (concurrent SlideRule downloads, default 1), TDS_JOBS (concurrent THREDDS prior downloads, default 2)
set -u
HERE="$(cd "$(dirname "$0")/.." && pwd)"
BATCH="$1"; shift
if [ -x /opt/anaconda3/envs/2026-icesat2-tracks/bin/snakemake ]; then
  SM=/opt/anaconda3/envs/2026-icesat2-tracks/bin/snakemake
else
  SM=$HOME/.conda/envs/2026-icesat2-tracks/bin/snakemake
fi
cd "$HERE"
exec "$SM" --snakefile Snakefile --config batch="$BATCH" \
  --cores "${CORES:-16}" --resources heavy="${JOBS:-4}" sliderule="${SL_JOBS:-1}" thredds="${TDS_JOBS:-2}" \
  --keep-going --rerun-incomplete --printshellcmds "$@"
