# analysis_db_claude — batch pipeline for ICESat-2 wave spectra

Rewrite of `analysis_db/` as a batch pipeline: a batch is a lat/lon box × time window; every ATL03
track in it is processed B01 → A02 → B02 → B03 → B04 → B05 → B06 → C01, with one status record
per (stage, track), Snakemake as the scheduler and a static HTML gallery for screening.

```
pipeline_config.py   config pick (config/config.json on cerberus, config_local.json on the laptop),
                     sys.path for modules/, mconfig, M/MT/col, figure style, paths_for(), save_fig()
pipeline_status.py   StageRun context manager -> work/<batch>/status/<stage>/<ID>.json (+ .ok), logs/
pipeline_params.py   batches/<batch>.toml, params/<version>.toml, params_hash, materialize()
pipeline_dag.py      stage order, dependencies, threads/resources, key figure per stage
Snakefile            per-track rules; outputs are the .ok sentinels
gallery.py           plots/<hemis>/<batch>/index.html + per-track / per-stage pages
stages/              B00_discover, B01_sliderule_load, A02_ww3_prior, B02_spectra, B03_plot_spectra,
                     B04_angle, B05_define_angle, B06_correct, C01_collect, C01_index
batches/*.toml       batch definitions (region, time, chunking, selection, sliderule nodes)
params/v1.toml       per-stage processing parameters (hashed into every status record)
tools/               run_batch.sh, discover.sh, serve_gallery.sh, tunnel.sh, sync_groundtracks.sh
```

## Commands

```bash
cd analysis_db_claude
tools/discover.sh SH_dev_small                 # B00: CMR query -> work/<batch>/tracks.csv + overview map
tools/run_batch.sh SH_dev_small -n --reason    # dry run: what would run and why
tools/run_batch.sh SH_dev_small                # everything (B01 per chunk, then per track up to C01, index)
tools/run_batch.sh SH_dev_small --until B04    # stop after a stage
tools/run_batch.sh SH_dev_small -R B04         # rerun a stage (+ downstream) after editing its script
tools/run_batch.sh SH_dev_small -f ../data/work/SH_dev_small/status/B04/SH_20190502_05200310.ok   # one track
python gallery.py SH_dev_small                 # (re)build the HTML pages (also done automatically at the end of a run)
```

Environment: `JOBS` (heavy jobs, default 4) and `CORES` (default 16) for `run_batch.sh`. Each job caps
its BLAS threads (`OMP_NUM_THREADS` from the rule's `threads`).

Single stage by hand (e.g. while debugging):
```bash
MPLBACKEND=Agg python stages/B02_spectra.py SH_20190502_05200310 SH_dev_small
IS2_SKIP_UPSTREAM=1 ... python stages/B04_angle.py <ID> <batch>   # ignore missing upstream status records
```
Interactive: open a stage script, run the cells inside `run_stage` with
`run = StageRun.debug('B02', ID, batch_key); prm = load_params(batch_key)['B02']; P = paths_for(batch_key, ID)`.

## On cerberus

```bash
ssh cerberus
cd ~/2021_ICESat2_tracks_develop && git pull
cd analysis_db_claude
tools/discover.sh SH_box140_155_201905
tmux new -d -s is2 "tools/run_batch.sh SH_box140_155_201905 > ../logs/run_SH_box140_155_201905.log 2>&1"
tools/serve_gallery.sh            # http.server on 127.0.0.1:8765 in tmux session 'gallery'
```
On the laptop: `tools/tunnel.sh 8765 SH_box140_155_201905` then open the printed URL.
The RGT shapefiles live in `/srv/hades/2021_ICESat2_tracks/groundtracks/` (`paths.groundtracks` in
`config/config.json`); copy them with `tools/sync_groundtracks.sh`.

## Status records

`work/<batch>/status/<stage>/<ID>.json`: `status` ∈ success | fail | skip | blocked | running,
`reason`, `error_type/error_msg/traceback_tail`, timing, host, git hash, `params_hash`, `params`,
`upstream`, `outputs`, `figures`, `info` (stage metrics). `<ID>.ok` exists for success/skip/blocked
(Snakemake target); a fail has no `.ok` and is retried on the next run. `logs/<stage>/<ID>.log` has
the stage's stdout/stderr, `<ID>.launcher.log` whatever happened before the stage started.

B01 runs per time chunk (`chunk_days` in the batch file): `status/B01/chunk<n>.json` for the
download job and one record per expected track. A track that failed inside a chunk is retried only
with `-R B01` (re-downloads the chunk).

## Reruns

- Edit `params/v1.toml` → `materialize()` rewrites `work/<batch>/params/<stage>.json` only if that
  section changed → Snakemake reruns that stage and everything downstream.
- Edit a stage script → its mtime triggers the rerun (`-R <stage>` if in doubt).
- Edit a module under `modules/` → `-R <stage>` for the stages that use it.
- A skipped/blocked track keeps its `.ok` and is not retried; delete the `.ok` (or the `status/<stage>/<ID>.*`) to force it.

## Database

`C01_collect` writes `work/<batch>/C01_database/C01_<ID>.nc` (corrected PSD(x, beam, k), angle PDF,
prior scalars, metadata; float32, zlib, ~1 MB) and `C01_index.py` builds `index.csv`/`index.parquet`
over the batch.
