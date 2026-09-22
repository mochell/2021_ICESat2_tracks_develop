#!/bin/bash
# Forward the cerberus gallery port to the laptop and print the URL.
# usage: tools/tunnel.sh [port] [batch_key]
PORT="${1:-8765}"
BATCH="${2:-}"
ssh -O check cerberus >/dev/null 2>&1 || { echo "no ssh master; run: ssh -fN cerberus"; exit 1; }
ssh -O forward -L "$PORT:localhost:$PORT" cerberus && echo "tunnel up: http://localhost:$PORT/"
[ -n "$BATCH" ] && echo "http://localhost:$PORT/${BATCH%%_*}/$BATCH/index.html"
