#!/bin/bash
# Copy the RGT shapefiles (never in git) to cerberus.
HERE="$(cd "$(dirname "$0")/../.." && pwd)"
ssh -o BatchMode=yes cerberus 'mkdir -p /srv/hades/2021_ICESat2_tracks/groundtracks'
rsync -av --include='IS2_mission_points_*_RGT_all.*' --exclude='*' \
  "$HERE/data/groundtracks/" cerberus:/srv/hades/2021_ICESat2_tracks/groundtracks/
