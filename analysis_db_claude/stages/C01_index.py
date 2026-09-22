"""
C01 index: one row per track over all C01_<ID>.nc files of a batch (+ status of every stage).

    python stages/C01_index.py <batch_key>

Writes work/<batch>/C01_database/index.csv and index.parquet
"""
import sys
import glob
import os
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline_config import np, pd, xr, paths_for
from pipeline_status import iter_status
from pipeline_dag import TRACK_STAGES


def build(batch_key):
    P = paths_for(batch_key)
    d = P.stage_dir('C01_database')
    st = iter_status(batch_key)
    status_wide = st.pivot_table(index='ID', columns='stage', values='status', aggfunc='first') if len(st) else pd.DataFrame()

    rows = []
    for f in sorted(glob.glob(d + 'C01_*.nc')):
        with xr.open_dataset(f) as D:
            a = dict(D.attrs)
            psd = D['gFT_PSD_data'].sel(beam='weig') if 'weig' in D.beam.values else D['gFT_PSD_data'].mean('beam')
            row = {k: a.get(k) for k in ('ID', 'batch_key', 'hemis', 'date', 'rgt', 'cycle', 'segment', 'ascending',
                                          'start_lon', 'start_lat', 'end_lon', 'end_lat', 'start_time',
                                          'best_guess_incident_angle_deg', 'theta_applied',
                                          'prior_hs', 'prior_fp', 'prior_dir', 'prior_spr', 'prior_ice')}
            row.update(n_x=int(D.x.size), n_k=int(D.k.size),
                       x_min_km=float(D.x.min() / 1e3), x_max_km=float(D.x.max() / 1e3),
                       psd_max=float(np.nanmax(psd.values)) if psd.size else np.nan,
                       k_peak=float(D.k.values[np.nanargmax(np.nanmean(psd.values, axis=0))]) if psd.size else np.nan,
                       coverage_mean=float(np.nanmean(D['N_per_stancil_fraction'].values)) if 'N_per_stancil_fraction' in D else np.nan,
                       file=os.path.relpath(f, P.batch_work), size_mb=round(os.path.getsize(f) / 1e6, 2))
        for s in TRACK_STAGES:
            row[f'status_{s}'] = status_wide.loc[row['ID'], s] if len(status_wide) and row['ID'] in status_wide.index and s in status_wide else None
        rows.append(row)
    idx = pd.DataFrame(rows)
    idx.to_csv(d + 'index.csv', index=False)
    try:
        idx.to_parquet(d + 'index.parquet', index=False)
    except Exception as e:
        print('parquet not written:', e)
    print(f'{len(idx)} tracks -> {d}index.csv')
    return idx


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    idx = build(args[0] if args else 'SH_dev_small')
    if len(idx):
        print(idx[['ID', 'n_x', 'best_guess_incident_angle_deg', 'prior_hs', 'size_mb']].to_string())
