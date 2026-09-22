#!/bin/bash
# Serve the plots directory (static gallery pages) on cerberus, bound to localhost only.
# usage: tools/serve_gallery.sh [port]     then from the laptop: tools/tunnel.sh [port]
PORT="${1:-8765}"
HERE="$(cd "$(dirname "$0")/.." && pwd)"
if [ -x /opt/anaconda3/envs/2026-icesat2-tracks/bin/python ]; then
  PY=/opt/anaconda3/envs/2026-icesat2-tracks/bin/python
else
  PY=$HOME/.conda/envs/2026-icesat2-tracks/bin/python
fi
PLOTS=$("$PY" -c "import sys; sys.path.insert(0,'$HERE'); import pipeline_config as c; print(c.mconfig['paths']['plot'])")
tmux kill-session -t gallery 2>/dev/null
tmux new -d -s gallery "cd '$PLOTS' && '$PY' -m http.server $PORT --bind 127.0.0.1"
echo "serving $PLOTS on 127.0.0.1:$PORT (tmux session 'gallery')"
