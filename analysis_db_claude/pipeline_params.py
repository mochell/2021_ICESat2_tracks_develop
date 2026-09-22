"""
Batch definitions (batches/<batch_key>.toml) and per-stage parameters (params/<version>.toml).

    batch  = load_batch('SH_box140_155_201905')
    prm    = load_params('SH_box140_155_201905')['B02']
    h      = params_hash(prm)
    files  = materialize('SH_box140_155_201905')   # work/<batch>/params/<stage>.json, write-if-changed
"""
import json
import hashlib
import tomllib
from pathlib import Path

from pipeline_config import HERE, paths_for, MT

BATCH_DIR = HERE / 'batches'
PARAMS_DIR = HERE / 'params'


def _read_toml(path):
    with open(path, 'rb') as f:
        return tomllib.load(f)


def load_batch(batch_key):
    """batches/<batch_key>.toml -> dict. Falls back to work/<batch>/batch.json (written by B00)."""
    p = BATCH_DIR / f'{batch_key}.toml'
    if p.exists():
        b = _read_toml(p)
    else:
        pj = Path(paths_for(batch_key).batch_work) / 'batch.json'
        if not pj.exists():
            raise FileNotFoundError(f'no batch definition: {p} or {pj}')
        with open(pj) as f:
            b = json.load(f)
    b['batch'].setdefault('key', batch_key)
    if b['batch']['key'] != batch_key:
        raise ValueError(f"batch key mismatch: file says {b['batch']['key']}, asked for {batch_key}")
    b['batch'].setdefault('hemis', batch_key.split('_')[0])
    b['batch'].setdefault('params_version', 'v1')
    b.setdefault('selection', {})
    b['selection'].setdefault('max_tracks', 0)
    b['selection'].setdefault('include_ids', [])
    b['selection'].setdefault('exclude_ids', [])
    b['selection'].setdefault('require_rgt', True)
    b.setdefault('time', {}).setdefault('chunk_days', 0)
    b.setdefault('sliderule', {})
    b['sliderule'].setdefault('desired_nodes', 1)
    b['sliderule'].setdefault('time_to_live', 90)
    return b


def load_params(batch_key=None, version=None):
    """all stage sections of params/<version>.toml (version taken from the batch if not given)"""
    if version is None:
        version = load_batch(batch_key)['batch']['params_version']
    p = PARAMS_DIR / f'{version}.toml'
    prm = _read_toml(p)
    prm.setdefault('version', version)
    return prm


def params_hash(section):
    return hashlib.sha1(json.dumps(section, sort_keys=True, default=str).encode()).hexdigest()[:10]


def materialize(batch_key):
    """
    Write work/<batch>/params/<stage>.json for every stage section, only if the content changed,
    so the file mtimes are stable and can be used as Snakemake inputs.
    Returns {stage: path}.
    """
    P = paths_for(batch_key)
    MT.mkdirs_r(P.params)
    prm = load_params(batch_key)
    out = {}
    for stage, section in prm.items():
        if not isinstance(section, dict):
            continue
        text = json.dumps({'version': prm['version'], 'params_hash': params_hash(section), **section},
                          indent=2, sort_keys=True, default=str)
        path = Path(P.params) / f'{stage}.json'
        if not path.exists() or path.read_text() != text:
            path.write_text(text)
        out[stage] = str(path)
    return out


if __name__ == '__main__':
    import sys
    key = sys.argv[1] if len(sys.argv) > 1 else 'SH_box140_155_201905'
    b = load_batch(key)
    print(json.dumps(b, indent=2, default=str))
    prm = load_params(key)
    for k, v in prm.items():
        if isinstance(v, dict):
            print(k, params_hash(v))
