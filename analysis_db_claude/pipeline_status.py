"""
Status records per (stage, track): work/<batch>/status/<stage>/<ID>.json (+ <ID>.ok sentinel)
and a log per (stage, track): work/<batch>/logs/<stage>/<ID>.log.

    with StageRun('B02', ID, batch_key, params=prm, script=__file__) as run:
        require_upstream(run, ['B01'])
        ...
        if bad: raise SkipTrack('nan fraction > max', nan_fraction=0.97)
        run.info(n_k=440)

Exit semantics
  normal exit    -> status 'success', json + .ok
  SkipTrack      -> 'skip'    (deliberate rejection), json + .ok, exception swallowed, exit code 0
  BlockedTrack   -> 'blocked' (upstream not usable),   json + .ok, swallowed, exit code 0
  other          -> 'fail', json only (no .ok, so schedulers retry), exception re-raised -> exit 1
"""
import os
import sys
import json
import time
import socket
import getpass
import platform
import traceback
import subprocess
import datetime as dt
from pathlib import Path

from pipeline_config import ROOT, paths_for, MT
from pipeline_dag import STAGE_DIR, DEPS

SCHEMA_VERSION = 1
STATUS_VALUES = ('success', 'fail', 'skip', 'blocked', 'running', 'not_run')


class SkipTrack(Exception):
    """Deliberate rejection of a track by a stage (physics / data quality)."""

    def __init__(self, reason, **info):
        super().__init__(reason)
        self.reason = reason
        self.info = info


class BlockedTrack(SkipTrack):
    """Upstream stage did not succeed for this track."""


# ----------------------------------------------------------------------------- helpers
def _now():
    return dt.datetime.now().replace(microsecond=0)


def _git_state():
    try:
        h = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, capture_output=True,
                           text=True, timeout=5).stdout.strip()
        d = subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=ROOT,
                           capture_output=True, text=True, timeout=5).stdout.strip()
        return h or None, bool(d)
    except Exception:
        return None, None


def status_file(batch_key, stage, ID):
    return Path(paths_for(batch_key).status) / stage / f'{ID}.json'


def ok_file(batch_key, stage, ID):
    return Path(paths_for(batch_key).status) / stage / f'{ID}.ok'


def log_file(batch_key, stage, ID):
    return Path(paths_for(batch_key).logs) / stage / f'{ID}.log'


def read_status(batch_key, stage, ID):
    p = status_file(batch_key, stage, ID)
    if not p.exists():
        return None
    with open(p) as f:
        return json.load(f)


def write_status(rec):
    """write a record dict (used by StageRun and by B01 for the tracks it produces)"""
    p = status_file(rec['batch_key'], rec['stage'], rec['ID'])
    MT.mkdirs_r(str(p.parent))
    tmp = p.with_suffix('.json.tmp')
    with open(tmp, 'w') as f:
        json.dump(rec, f, indent=2, default=str)
    os.replace(tmp, p)
    okp = p.with_suffix('.ok')
    if rec['status'] in ('success', 'skip', 'blocked'):
        okp.touch()
    elif okp.exists():
        okp.unlink()
    return p


def iter_status(batch_key, stages=None):
    """long table with one row per status json: stage, ID, status, reason, error_type, runtime_s, ..."""
    import pandas as pd
    rows = []
    root = Path(paths_for(batch_key).status)
    if not root.exists():
        return pd.DataFrame(columns=['stage', 'ID', 'status'])
    for sdir in sorted(root.iterdir()):
        if not sdir.is_dir() or (stages and sdir.name not in stages):
            continue
        for p in sorted(sdir.glob('*.json')):
            try:
                with open(p) as f:
                    r = json.load(f)
            except Exception as e:
                r = {'stage': sdir.name, 'ID': p.stem, 'status': 'fail', 'error_type': 'BadStatusFile',
                     'error_msg': str(e)}
            rows.append({k: r.get(k) for k in ('stage', 'ID', 'batch_key', 'status', 'reason', 'error_type',
                                                'error_msg', 't_start', 't_end', 'runtime_s', 'host',
                                                'params_hash', 'params_version', 'git_hash', 'log')}
                        | {'info': r.get('info', {}), 'figures': r.get('figures', []),
                           'outputs': r.get('outputs', [])})
    return pd.DataFrame(rows)


def require_upstream(run, stages=None):
    """
    Raise BlockedTrack unless every upstream stage has status 'success' for run.ID.
    Defaults to pipeline_dag.DEPS[run.stage]. Records the upstream statuses in the run.
    """
    stages = DEPS[run.stage] if stages is None else stages
    for s in stages:
        st = read_status(run.batch_key, s, run.ID)
        status = st['status'] if st else 'not_run'
        run.upstream[s] = status
        if os.environ.get('IS2_SKIP_UPSTREAM'):      # debugging on batches without status records
            continue
        if status != 'success':
            reason = (st or {}).get('reason') or (st or {}).get('error_msg') or ''
            raise BlockedTrack(f'{s} {status}' + (f': {reason}' if reason else ''))


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, s):
        for st in self.streams:
            st.write(s)
            st.flush()

    def flush(self):
        for st in self.streams:
            st.flush()

    def isatty(self):
        return False


# ----------------------------------------------------------------------------- StageRun
class StageRun:
    def __init__(self, stage, ID, batch_key, params=None, script=None, tee_log=True, write=True):
        self.stage, self.ID, self.batch_key = stage, ID, batch_key
        self.params = params or {}
        self.script = str(Path(script).resolve().relative_to(ROOT)) if script else None
        self.tee_log = tee_log
        self.write = write            # False: debug mode, nothing is written
        self.status = 'running'
        self.reason = None
        self.info_dict = {}
        self.upstream = {}
        self.outputs = []
        self.figures = []
        self.P = paths_for(batch_key, ID)
        self._log_fh = None
        self._orig = None

    @classmethod
    def debug(cls, stage, ID, batch_key, params=None):
        """for running # %% cells interactively: same object, no files written"""
        r = cls(stage, ID, batch_key, params=params, tee_log=False, write=False)
        r.t_start = _now()
        return r

    # -- context manager
    def __enter__(self):
        self.t_start = _now()
        self._t0 = time.time()
        if self.write and self.tee_log:
            lf = log_file(self.batch_key, self.stage, self.ID)
            MT.mkdirs_r(str(lf.parent))
            self._log_fh = open(lf, 'a')
            self._log_fh.write(f'\n===== {self.stage} {self.ID} {self.batch_key} START {self.t_start.isoformat()} '
                               f'host={socket.gethostname()} pid={os.getpid()}\n')
            self._orig = (sys.stdout, sys.stderr)
            sys.stdout = _Tee(sys.__stdout__, self._log_fh)
            sys.stderr = _Tee(sys.__stderr__, self._log_fh)
        if self.write:
            write_status(self._record('running'))
        return self

    def __exit__(self, exc_type, exc, tb):
        err = {}
        if exc_type is None:
            self.status = 'success'
        elif issubclass(exc_type, BlockedTrack):
            self.status, self.reason = 'blocked', exc.reason
            self.info_dict.update(exc.info)
        elif issubclass(exc_type, SkipTrack):
            self.status, self.reason = 'skip', exc.reason
            self.info_dict.update(exc.info)
        else:
            self.status = 'fail'
            tail = traceback.format_exception(exc_type, exc, tb)
            err = {'error_type': 'Interrupted' if exc_type is KeyboardInterrupt else exc_type.__name__,
                   'error_msg': str(exc)[:500],
                   'traceback_tail': ''.join(tail).splitlines()[-30:]}
            print(''.join(tail), file=sys.stderr)
        self.t_end = _now()
        self.runtime_s = round(time.time() - self._t0, 1)
        if self.write:
            self._detect_products()
            write_status(self._record(self.status, **err))
        msg = f'===== {self.stage} {self.ID} {self.status.upper()}'
        if self.reason:
            msg += f' ({self.reason})'
        print(msg + f' {self.runtime_s}s')
        if self._orig:
            sys.stdout, sys.stderr = self._orig
            self._log_fh.close()
        # swallow skip/blocked, re-raise real errors
        return exc_type is not None and issubclass(exc_type, SkipTrack)

    # -- API for the stage body
    def info(self, **kv):
        self.info_dict.update(kv)

    def add_output(self, path):
        self.outputs.append(str(path))

    # -- internals
    def _detect_products(self):
        t0 = self._t0 - 1
        sdir = STAGE_DIR.get(self.stage)
        if sdir:
            d = Path(self.P.batch_work) / sdir
            if d.exists():
                for p in d.glob(f'*{self.ID}*'):
                    if p.is_file() and p.stat().st_mtime >= t0:
                        self.outputs.append(str(p.relative_to(self.P.batch_work)))
        if self.P.plot_track and Path(self.P.plot_track).exists():
            for p in Path(self.P.plot_track).rglob('*'):
                rel = str(p.relative_to(self.P.plot_track))
                # only this stage's figures (name or sub-folder prefix), written during this run
                mine = p.name.startswith(self.stage) or rel.split('/')[0].startswith(self.stage)
                if mine and p.is_file() and p.suffix in ('.png', '.pdf') and p.stat().st_mtime >= t0:
                    self.figures.append(rel)
        self.outputs = sorted(set(self.outputs))
        self.figures = sorted(set(self.figures))

    def _record(self, status, **extra):
        from pipeline_params import params_hash
        git_hash, git_dirty = _git_state()
        rec = {
            'schema_version': SCHEMA_VERSION,
            'stage': self.stage, 'ID': self.ID, 'batch_key': self.batch_key,
            'status': status, 'reason': self.reason,
            'error_type': None, 'error_msg': None, 'traceback_tail': [],
            't_start': self.t_start.isoformat(),
            't_end': getattr(self, 't_end', None) and self.t_end.isoformat(),
            'runtime_s': getattr(self, 'runtime_s', None),
            'host': socket.gethostname().split('.')[0], 'user': getpass.getuser(), 'pid': os.getpid(),
            'python': platform.python_version(),
            'git_hash': git_hash, 'git_dirty': git_dirty,
            'script': self.script,
            'params_version': self.params.get('version') if isinstance(self.params, dict) else None,
            'params_hash': params_hash(self.params),
            'params': self.params,
            'upstream': self.upstream,
            'outputs': self.outputs, 'figures': self.figures,
            'log': f'logs/{self.stage}/{self.ID}.log',
            'info': self.info_dict,
        }
        rec.update(extra)
        return rec


# ----------------------------------------------------------------------------- self test
if __name__ == '__main__':
    """python pipeline_status.py  -> exercises the four exit paths in a throwaway batch"""
    import shutil
    bk = 'SH_statustest'
    P = paths_for(bk)
    shutil.rmtree(P.batch_work, ignore_errors=True)

    with StageRun('B02', 'SH_00000000_00000000', bk, params={'a': 1}, script=__file__) as r:
        print('hello from success')
    assert r.status == 'success' and ok_file(bk, 'B02', r.ID).exists()

    with StageRun('B03', r.ID, bk, params={}, script=__file__) as r2:
        raise SkipTrack('deliberate', nan_fraction=0.97)
    assert r2.status == 'skip' and ok_file(bk, 'B03', r.ID).exists()
    assert read_status(bk, 'B03', r.ID)['info']['nan_fraction'] == 0.97

    with StageRun('B04', r.ID, bk, params={}, script=__file__) as r3:
        require_upstream(r3, ['B02', 'B03'])       # B03 is skip -> blocked
    assert r3.status == 'blocked' and 'B03 skip' in r3.reason

    try:
        with StageRun('B05', r.ID, bk, params={}, script=__file__) as r4:
            raise ValueError('boom')
    except ValueError:
        pass
    st = read_status(bk, 'B05', r.ID)
    assert r4.status == 'fail' and st['error_type'] == 'ValueError' and not ok_file(bk, 'B05', r.ID).exists()
    assert 'hello from success' in log_file(bk, 'B02', r.ID).read_text()
    print(iter_status(bk)[['stage', 'ID', 'status', 'reason', 'error_type', 'runtime_s']])
    shutil.rmtree(P.batch_work, ignore_errors=True)
    print('status self-test OK')
