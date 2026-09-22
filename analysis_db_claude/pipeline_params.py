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


def _unwrap_lons(pts):
    """make longitudes continuous along the ring (no jumps > 180), so a polygon given in
    -180..180 that crosses the antimeridian becomes a single planar polygon"""
    out = [list(pts[0])]
    for lon, lat in pts[1:]:
        prev = out[-1][0]
        while lon - prev > 180:
            lon -= 360
        while lon - prev < -180:
            lon += 360
        out.append([lon, lat])
    return out


def region_polygon(batch):
    """
    Region of a batch as the dict sct.create_polygons() returns, plus antimeridian handling:
    {'list': [{'lat','lon'}, ... closed ring]   the polygon in its continuous (unwrapped) longitudes
     'parts': [ring, ring, ...]                 the same region cut at +-180 into rings within -180..180
                                                (one entry unless the region crosses the antimeridian)
     'shapely': (Multi)Polygon of the parts     for masking the RGT shapefile (which is in -180..180)
     'lons': [min, max], 'lats': [min, max]     bounds in the unwrapped frame
     'kind': 'box' | 'polygon', 'crosses_antimeridian': bool}

    [region] accepts either
        lat = [lat0, lat1]; lon = [lon0, lon1]                       axis-aligned box
    or  polygon = [[lon, lat], [lon, lat], ...]                       any simple polygon (tilted box,
                                                                      quadrilateral, ...), >= 3 vertices,
                                                                      closing vertex optional.
    Longitudes may be given in -180..180 or continued past +-180 (e.g. -192.4 for 167.6 E) when the
    polygon crosses the antimeridian; both are handled.
    """
    from shapely.geometry import Polygon, box, MultiPolygon
    from shapely.geometry.polygon import orient
    from shapely.ops import unary_union
    reg = batch['region']
    if reg.get('polygon'):
        pts = [(float(p[0]), float(p[1])) for p in reg['polygon']]
        if len(pts) < 3:
            raise ValueError('[region] polygon needs at least 3 [lon, lat] vertices')
        pts = _unwrap_lons(pts)
        pg = Polygon(pts)
        if not pg.is_valid:
            raise ValueError('[region] polygon is not a valid simple polygon')
        kind = 'polygon'
    else:
        lat = sorted(float(v) for v in reg['lat'])
        lon = sorted(float(v) for v in reg['lon'])
        pg = Polygon([(lon[1], lat[1]), (lon[0], lat[1]), (lon[0], lat[0]), (lon[1], lat[0])])
        kind = 'box'
    pg = orient(pg, sign=1.0)                       # counter-clockwise
    ring = list(pg.exterior.coords)                 # closed (first == last)
    lons, lats = [c[0] for c in ring], [c[1] for c in ring]

    # cut at the antimeridian: every piece is shifted back into -180..180
    parts = []
    for shift in (0.0, -360.0, 360.0):
        piece = pg.intersection(box(-180 + shift, -90, 180 + shift, 90))
        if piece.is_empty:
            continue
        geoms = piece.geoms if hasattr(piece, 'geoms') else [piece]
        for g in geoms:
            if g.geom_type != 'Polygon' or g.area == 0:
                continue
            g = orient(Polygon([(x - shift, y) for x, y in g.exterior.coords]), sign=1.0)
            parts.append(g)
    crosses = len(parts) > 1 or min(lons) < -180 or max(lons) > 180
    return {'list': [{'lat': la, 'lon': lo} for lo, la in ring],
            'parts': [[{'lat': la, 'lon': lo} for lo, la in g.exterior.coords] for g in parts],
            'shapely': parts[0] if len(parts) == 1 else unary_union(parts),
            'lons': [min(lons), max(lons)], 'lats': [min(lats), max(lats)],
            'kind': kind, 'crosses_antimeridian': bool(crosses)}


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
