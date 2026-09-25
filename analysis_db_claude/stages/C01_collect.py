# %%
"""
C01: collect the results of one track into a compact netCDF for the wave-spectra database.

    python stages/C01_collect.py <ID> <batch_key>

Reads  work/<batch>/B06_corrected_separated/B06_<ID>_gFT_k_corrected.nc   (B06)
       work/<batch>/B04_angle/B05_<ID>_angle_pdf.nc                        (B05)
       work/<batch>/A02_prior/A02_<ID>.h5                                   (A02)
       work/<batch>/A01b_ID/A01b_ID_<ID>.json, tracks.csv                   (B01, B00)
       status/*/<ID>.json                                                    (all stages)
Writes work/<batch>/C01_database/C01_<ID>.nc   (~1 MB, float32 + zlib)
"""
import os
import sys
import json
import datetime as dt
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline_config import mconfig, np, pd, xr, MT, paths_for, cli_args
from pipeline_status import StageRun, SkipTrack, require_upstream, read_status
from pipeline_params import load_params
from pipeline_dag import TRACK_STAGES

STAGE = 'C01'
PRIOR_SCALARS = ['hs', 'fp', 'dir', 'dp', 'spr', 't01', 't02', 'ice', 'lon', 'lat']


def run_stage(ID, batch_key, prm, run):
    P = paths_for(batch_key, ID)
    require_upstream(run, ['B06', 'A02'])     # B05 angle optional
    save_path = P.stage_dir('C01_database')

    # %% load
    Gk = xr.open_dataset(P.stage_dir('B06_corrected_separated') + f'B06_{ID}_gFT_k_corrected.nc')
    angle_file = P.stage_dir('B04_angle') + f'B05_{ID}_angle_pdf.nc'
    b05 = read_status(batch_key, 'B05', ID) or {}
    Ga = xr.open_dataset(angle_file) if (b05.get('status') == 'success' and os.path.exists(angle_file)) else None
    Prior = MT.load_pandas_table_dict('/A02_' + ID, P.stage_dir('A02_prior'))['priors_hindcast']
    with open(P.stage_dir('A01b_ID') + f'A01b_ID_{ID}.json') as f:
        IDj = json.load(f)
    tr = None
    if os.path.exists(P.batch_work + 'tracks.csv'):
        tr = pd.read_csv(P.batch_work + 'tracks.csv')
        tr = tr[tr.ID == ID].iloc[0] if (tr.ID == ID).any() else None

    # %% spectra: selected variables, all beams incl. the weighted mean beam 'weig'
    keep = [v for v in prm['k_variables'] if v in Gk.variables]
    D = Gk[keep].copy()
    for c in ('k_corrected', 'x_corrected', 'k_lim', 'N_per_stancil_fraction', 'L', 'Lpoints'):
        if c in Gk.coords or c in Gk.variables:
            D.coords[c] = Gk[c]
    D = D.transpose('x', 'beam', 'k', missing_dims='ignore')

    # %% angle PDF on its own x axis
    if Ga is not None:
        A = xr.Dataset({'angle_PDF': Ga.weighted_angle_PDF.rename({'x': 'x_angle'}),
                        'angle_PDF_smth': Ga.weighted_angle_PDF_smth.rename({'x': 'x_angle'})})
        A.coords['N_data_angle'] = Ga.N_data.rename({'x': 'x_angle'})
        D = xr.merge([D, A])

    # %% scalars and metadata
    theta = float(Gk.attrs.get('best_guess_incident_angle', np.nan))
    st = {s: read_status(batch_key, s, ID) for s in TRACK_STAGES}
    b06 = (st.get('B06') or {}).get('info', {})
    attrs = {
        'ID': ID, 'batch_key': batch_key, 'hemis': P.hemis,
        'granules': json.dumps(IDj['tracks'].get('ATL03')),
        'rgt': int(tr.rgt) if tr is not None else -1, 'cycle': int(tr.cycle) if tr is not None else -1,
        'segment': int(tr.segment) if tr is not None else -1, 'date': str(tr.date) if tr is not None else ID.split('_')[1],
        'ascending': int(bool(IDj['pars'].get('ascending', IDj['pars'].get('poleward')))),
        'x_reference_m': float(IDj['pars'].get('x_reference_m', np.nan)),
        'start_lon': IDj['pars']['start']['longitude'], 'start_lat': IDj['pars']['start']['latitude'],
        'end_lon': IDj['pars']['end']['longitude'], 'end_lat': IDj['pars']['end']['latitude'],
        'start_time': dt.datetime.fromtimestamp(IDj['pars']['start']['delta_time'], dt.UTC).replace(tzinfo=None).isoformat(),
        'best_guess_incident_angle_rad': theta, 'best_guess_incident_angle_deg': float(np.rad2deg(theta)),
        'theta_applied': int(bool(b06.get('theta_applied', not np.isnan(theta)))),
        'angle_status': b05.get('status', 'not_run'), 'angle_reason': b05.get('reason') or '',
        'L': float(Gk.attrs.get('L', np.nan)), 'Lpoints': int(Gk.attrs.get('Lpoints', 0)),
        'created': dt.datetime.now().isoformat(),
    }
    for k in PRIOR_SCALARS:
        if k in Prior.index:
            attrs['prior_' + k] = float(Prior.loc[k]['mean'])
            attrs['prior_' + k + '_std'] = float(Prior.loc[k]['std'])
    for s, r in st.items():
        if r:
            attrs[f'{s}_status'] = r['status']
            attrs[f'{s}_params_hash'] = r.get('params_hash') or ''
            attrs[f'{s}_git_hash'] = r.get('git_hash') or ''
    D.attrs = {k: v for k, v in attrs.items() if v is not None}

    # %% save
    enc = {}
    for v in D.data_vars:
        e = {'zlib': True, 'complevel': 4}
        if prm.get('float32', True) and np.issubdtype(D[v].dtype, np.floating):
            e['dtype'] = 'float32'
        enc[v] = e
    out = save_path + f'C01_{ID}.nc'
    D.to_netcdf(out, encoding=enc)
    size_mb = os.path.getsize(out) / 1e6
    print(f'saved {out} ({size_mb:.2f} MB)')
    run.info(size_mb=round(size_mb, 2), n_x=int(D.x.size), n_k=int(D.k.size), n_x_angle=int(D.x_angle.size) if 'x_angle' in D.dims else 0,
             theta_deg=attrs['best_guess_incident_angle_deg'], prior_hs=attrs.get('prior_hs'), prior_dir=attrs.get('prior_dir'))


def main(ID, batch_key):
    prm = load_params(batch_key)
    section = dict(prm[STAGE], version=prm['version'])
    with StageRun(STAGE, ID, batch_key, params=section, script=__file__) as run:
        run_stage(ID, batch_key, prm[STAGE], run)
    return run.status


if __name__ == '__main__':
    ID, batch_key = cli_args(sys.argv, default=('SH_20190502_05180312', 'SH_testSLsinglefile2'))
    sys.exit(0 if main(ID, batch_key) != 'fail' else 1)
