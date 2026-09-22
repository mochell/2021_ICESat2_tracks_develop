#!/bin/bash
# B00 for one batch: CMR query, RGT start points, tracks.csv, overview map.
# usage: tools/discover.sh <batch_key>
set -u
HERE="$(cd "$(dirname "$0")/.." && pwd)"
if [ -x /opt/anaconda3/envs/2026-icesat2-tracks/bin/python ]; then
  PY=/opt/anaconda3/envs/2026-icesat2-tracks/bin/python
else
  PY=$HOME/.conda/envs/2026-icesat2-tracks/bin/python
fi
cd "$HERE" && MPLBACKEND=Agg "$PY" stages/B00_discover.py "$1"
