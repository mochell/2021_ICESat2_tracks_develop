# Pipeline monitoring & rerun workflow — recommendation

*Drafted 2026-09-22 for the next coding session. Scope: `analysis_db/` pipeline run on cerberus over 1,000–10,000 ICESat-2 tracks, with repeated reruns after changing scripts/parameters.*

## Goal

1. **Status overview**: per stage, how many tracks succeeded / failed / were skipped / are blocked, and *why* they failed.
2. **Figure review**: pick a stage, figure type or failure class and scroll the corresponding figures across many tracks; click one track to see all its figures + log.
3. **Cheap reruns**: change a script or a parameter → only the affected stage and everything downstream is recomputed.

---

## 1. Findings in the current code

| Issue | Where | Consequence |
|---|---|---|
| Failure = *missing* `*_success.json`, no failure record | all stages (`B02_fail.json` is commented out) | cannot tell "not run" vs "crashed" vs "killed" vs "rejected for physics reasons" |
| Logs overwrite each other | makefile: B04/B05/B06 → `log/B04/$*.txt`; A01b/A01c/A02 → `log/A02/$*.txt` | a B06 run erases the B04 traceback |
| Make only reruns on changed *input files* | makefile; only B06 lists its script as a dependency | parameter / code changes are not picked up |
| Parameters hard-coded in scripts | e.g. SlideRule `params` dict in `B01_SL_load_batch.py` | no record of which parameters produced which output |
| Four different track-ID conventions | see §3 | IDs from different generations do not match |
| **Credential in code** | `B01_SL_load_batch.py`, commented `sliderule.authenticate(...)` line | password is in git history → **rotate it**, move to `.netrc` / env variable |
| Figures: ~40–100 per track, PNG/PDF mixed, nested sub-folders | `plots/SH/<batch>/<track>/…` | up to ~10⁵–10⁶ figures per batch; PDFs are slow to browse |

The plot layout itself is good: the folder is the track ID, the filename prefix is the stage (`B03_`, `B06_`, `PB03_`), and names encode beam (`gt1l`) and segment (`_x12`). A figure browser can parse everything it needs from the path.

---

## 2. Code changes

1. **Snakemake instead of the makefile.** Same per-track pattern rules; adds:
   - rerun on code or parameter change (`--rerun-triggers mtime params code input`), with each script listed as a rule input;
   - one log per stage per track (`log:` directive);
   - `--keep-going`, `-n` (dry run: what would run and why), `--summary`, `--report report.html` (figures grouped by stage/track);
   - a **checkpoint** to handle a track-ID change after the selection stage (if one remains, see §3);
   - `resources: nsidc=1` to cap parallel downloads (if downloads remain).
2. **Parameters in `config/params.yaml`**, one section per stage. Each script reads only its own section. Snakemake then reruns only that stage and what depends on it.
3. **Status record per track × stage.** Wrap each script body in a context manager (`pipeline_status.py`) that writes `work/<batch>/status/<stage>/<track>.json` containing:
   - `status`: `success` | `fail` | `skip` (deliberate physics rejection via `raise SkipTrack("reason")`);
   - error type, message, traceback tail;
   - runtime, host, git hash, hash of the parameters used;
   - list of figures produced.
4. **Always save a PNG next to every PDF** (`fig.savefig(path.png, dpi=100)`).
5. **Figure filenames**: either keep the track ID out of filenames inside the track folder, or always put it in the same position. This makes "figure families" (`B06_decomposition_gt1l_x*`) easy to group.
6. **Rotate the SlideRule credential** (see §1).

---

## 3. Track ID — make it consistent with SlideRule

Conventions currently in the repository:

| Source | Example | Last two digits |
|---|---|---|
| ATL07 granule | `20190219063727_08070201_005_01` | granule region |
| legacy A01b ID | `SH_20190219_08070210` | ATL03 granule region (10/12) |
| SlideRule `sct.create_ID_name` | `SH_YYYYMMDD_RRRRCC00` | "segment", default `00` |
| `data/work/database/batch_test` | `SH_20140511_123401_02_003` | – |

Legacy and SlideRule IDs have the same *shape* but different meaning in the last two digits (and possibly in the date), so the same pass gets different IDs. Adopt the SlideRule form `HEMIS_YYYYMMDD_RRRRCCSS` everywhere, after deciding the open points in §5.

---

## 4. Tools

### Status overview (decided)

- **DuckDB** status table built by `build_status_db.py`: scans the status JSONs and, for old batches, infers status from success markers + last `…Error:` line of the logs. It uses the stage dependency graph so that downstream stages show as **blocked**, not failed.
- Outputs: tables `runs` (track, stage, status, error, figures, …) and `overview` (track × stage); SQLite copy for **Datasette** (optional, faceted browsing of the table).
- Refreshed automatically at the end of every Snakemake run (`onsuccess` / `onerror`).

### Figure browsing (keep both, decide later)

| | **FiftyOne** | **Streamlit app** (`status_app.py`) |
|---|---|---|
| What it is | purpose-built browser for large image collections | ~200-line custom app, our code |
| Workflow | sidebar filters on any field (track, stage, figure family, beam, status), fast grid, grouping, tagging, export tagged list | stage/status/error filters → table → paged gallery of one stage across tracks; click row → all figures + log tail of that track; "flag for rerun" button |
| Setup | `pip install fiftyone`; indexing script (~30 lines) parses fields from paths + joins status; `fo.launch_app(ds, remote=True)`, tunnel port 5151 | `pip install streamlit duckdb pymupdf`; `streamlit run status_app.py --server.headless true`, tunnel port 8501 |
| Strength | scales to 10⁵–10⁶ figures, no UI code to maintain | table-first, status overview + figures in one place, fully adaptable |
| Weakness | PNG only; less tailored | we maintain it; slows down if too many images are shown at once (hence paging) |
| Effort | ~1 day | ~½ day to adapt to the real folder layout |

Both read the same folders and the same status table, so running both costs little.

### Comparing parameter versions (optional, later)

Only needed if v1 and v2 of the same track should be compared side by side. Ranked by setup effort vs. fit (all require logging figures from the scripts, which duplicates storage):

1. **Aim**: `pip install aim`, `aim up`; best UI for grouping by parameters; website offline, maintenance uncertain (repo still active Aug 2026).
2. **Trackio**: lightweight, SQLite, wandb-compatible API; young, untested at this scale.
3. **ClearML**: most complete; docker-compose server (Elasticsearch, MongoDB, Redis).
4. W&B (cloud) and MLflow (no cross-run image gallery): not recommended here.

A cheaper alternative is a "compare versions" tab in the Streamlit app, if outputs are written per parameter version (see §5).

---

## 5. Open decisions

1. **ID semantics**: what do the last two digits encode (region polygon index / sub-segment / always `00`)? Date from granule start or first photon (they differ across midnight)? Include a domain name?
2. **Does the A section survive?** If all data comes via SlideRule (`B01_SL_*`), the NSIDC download and ATL07→ATL03 stages go away and the pipeline starts with a SlideRule query stage that writes the track list.
3. **Skip handling**: a deliberately rejected track either writes its success marker with `status: skip` (downstream won't retry) or stays failed.
4. **Parameter versions**: overwrite in place, or write to `work/<batch>/<params_version>/` and `plots/…/<params_version>/` for side-by-side comparison.
5. **Figure browser**: FiftyOne, Streamlit, or both (see §4).
6. **Where the dashboard runs**: on cerberus behind an SSH tunnel (preferred; data stays there) or on the laptop after rsync.

### Naming decisions to confirm

- this file's name and location (`RECOMMENDATION_pipeline_monitoring.md`, repo root)
- status folder `work/<batch>/status/<stage>/<track>.json` and field names (`status`, `error_type`, `error_msg`, `params_hash`, `runtime_s`, `figures`)
- status values `success` / `fail` / `skip` / `blocked` / `not_run`
- log path `work/<batch>/logs/<stage>/<track>.log`
- `config/params.yaml` section names = stage names (`B01`, `B02`, …)
- parameter version label (`params_version: v1`, …)
- file names `build_status_db.py`, `pipeline_status.py`, `status_app.py`, `pipeline_status.duckdb`

---

## 6. Coding-session plan

1. Rotate the SlideRule credential; remove it from the script.
2. Settle §5 items 1–3 (ID, A section, skip).
3. Add `modules/pipeline_status.py`; wrap one stage (B04, the most failure-prone) and test on a small batch.
4. Run `build_status_db.py` on an existing batch on cerberus (works without step 3); check the counts against what you expect.
5. Move parameters of one stage into `config/params.yaml`.
6. Port the makefile to a `Snakefile`; `snakemake -n` must list exactly the failed + blocked tracks; then a real run with `--keep-going`.
7. Add PNG saving where only PDFs exist (B04, B05).
8. Set up the figure browser(s): adapt `status_app.py` to the real plot layout (figure families) and/or write the FiftyOne indexing script.
9. Roll out status wrapping + params to the remaining stages.

## 7. What has already been tested

Drafts of `Snakefile`, `params.yaml`, `pipeline_status.py`, `build_status_db.py`, `status_app.py` and a signac-dashboard example were written in the chat on 2026-09-22 and tested **on a synthetic batch** (320 tracks, real folder layout and ID format), not on real data:

- `build_status_db.py`: correct counts per stage/status; bug with shared legacy logs found and fixed (dependency graph → *blocked*).
- `status_app.py`: runs without errors in headless test mode (filters, table, gallery).
- `Snakefile`: dry run planned exactly the failed + blocked + skipped jobs (348), nothing already completed.

Assumptions to verify on real data: A01b always writes `A01b_<id>_success.json`; figure filenames in the Snakefile `report()` outputs; placeholder names in `params.yaml`.

---

## Links

- Snakemake: https://snakemake.github.io
- DuckDB: https://duckdb.org · Datasette: https://datasette.io
- FiftyOne: https://github.com/voxel51/fiftyone · live demo: https://try.fiftyone.ai
- Streamlit: https://docs.streamlit.io · row-selection tutorial: https://docs.streamlit.io/develop/tutorials/elements/dataframe-row-selections · similar app (benchmark case browser): https://huggingface.co/spaces/OpenHandsCommunity/evaluation
- Aim: https://github.com/aimhubio/aim · docs https://aimstack.readthedocs.io/en/latest/ · demos https://huggingface.co/aimstack
- Trackio: https://github.com/gradio-app/trackio · https://huggingface.co/docs/trackio
- ClearML: https://github.com/allegroai/clearml · server https://github.com/allegroai/clearml-server
- signac-dashboard: https://github.com/glotzerlab/signac-dashboard
