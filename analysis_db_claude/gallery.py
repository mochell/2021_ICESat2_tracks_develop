"""
Static HTML gallery + status report for one batch.

    python gallery.py <batch_key>

Reads  work/<batch>/status/*/*.json, work/<batch>/tracks.csv, work/<batch>/logs/, plots/<hemis>/<batch>/<ID>/**
Writes plots/<hemis>/<batch>/index.html                      track x stage grid, counts, filters
       plots/<hemis>/<batch>/_gallery/track/<ID>.html         all figures + status + log tail of one track
       plots/<hemis>/<batch>/_gallery/stage/<stage>.html      key figure of one stage for all tracks
       plots/<hemis>/<batch>/_gallery/thumbs/<ID>/*.jpg       360 px thumbnails (PIL)
       plots/<hemis>/<batch>/status.csv, status_summary.json

Serve the plots root with `tools/serve_gallery.sh` (cerberus) and open it through `tools/tunnel.sh`.
"""
import os
import sys
import json
import html
import glob
import datetime as dt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_config import pd, paths_for, MT
from pipeline_status import iter_status, log_file
from pipeline_dag import STAGES, TRACK_STAGES, DEPS, KEY_FIGURE, STATUS_ORDER, STATUS_COLOR

THUMB_W = 360
LOG_TAIL = 60
GRID_STAGES = ['B01'] + TRACK_STAGES

CSS = """
body{font-family:-apple-system,Helvetica,Arial,sans-serif;font-size:13px;margin:16px;color:#222;background:#fff;color-scheme:light}
h1{font-size:20px;margin:4px 0} h2{font-size:16px;margin:18px 0 6px}
table{border-collapse:collapse} td,th{padding:3px 6px;border:1px solid #ddd;vertical-align:top}
th{background:#f5f5f5;text-align:left}
.cell{display:inline-block;min-width:54px;text-align:center;padding:2px 4px;border-radius:3px;color:#fff;font-size:11px}
.cell a{color:#fff;text-decoration:none}
.badge{display:inline-block;padding:2px 6px;border-radius:3px;color:#fff;font-size:12px}
.grid{display:flex;flex-wrap:wrap;gap:10px}
:root{--card-w:380px}
.card{border:1px solid #ddd;padding:6px;width:var(--card-w);background:#fafafa}
.sizebar{margin:6px 0;font-size:12px} .sizebar button{margin-left:4px}
.card img{width:100%;display:block}
.small{color:#666;font-size:11px}
pre{background:#f7f7f7;border:1px solid #e5e5e5;padding:6px;font-size:11px;max-height:320px;overflow:auto}
.filters{margin:8px 0;padding:6px;background:#f5f5f5;border:1px solid #ddd}
.filters select,.filters input{margin-right:10px}
a{color:#1a5fb4}
tr.hidden,div.hidden{display:none}
"""

JS_SIZE = """
function setCardSize(px){
  document.documentElement.style.setProperty('--card-w', px+'px');
  document.querySelectorAll('img[data-thumb]').forEach(im=>{ im.src = (px>420 && im.dataset.full) ? im.dataset.full : im.dataset.thumb; });
  try{ localStorage.setItem('gallery_card_w', px); }catch(e){}
  document.querySelectorAll('.sizebar button').forEach(b=>b.style.fontWeight = (+b.dataset.px===px)?'bold':'normal');
}
document.addEventListener('DOMContentLoaded',()=>{ let px=380; try{ px=+localStorage.getItem('gallery_card_w')||380; }catch(e){} setCardSize(px); });
"""
JS_SORT = """
function sortCards(mode){
  const grid=document.querySelector('.grid'); if(!grid) return;
  const cards=Array.from(grid.children);
  const key=c=>mode==='id'?c.dataset.id:(c.dataset.mtime||'');
  cards.sort((a,b)=>{ const ka=key(a), kb=key(b); return mode==='newest' ? (kb>ka?1:kb<ka?-1:0) : (ka>kb?1:ka<kb?-1:0); });
  cards.forEach(c=>grid.appendChild(c));
  try{ localStorage.setItem('gallery_sort', mode); }catch(e){}
  const sel=document.getElementById('f_sort'); if(sel) sel.value=mode;
}
document.addEventListener('DOMContentLoaded',()=>{ let m='id'; try{ m=localStorage.getItem('gallery_sort')||'id'; }catch(e){} sortCards(m); });
"""
SORTBAR = ('<span class="sizebar">order: <select id="f_sort" onchange="sortCards(this.value)">'
           '<option value="id">track ID</option><option value="newest">newest first</option>'
           '<option value="oldest">oldest first</option></select></span>')

SIZEBAR = ('<div class="sizebar">image size: ' + ''.join(
    f'<button data-px="{px}" onclick="setCardSize({px})">{lab}</button>'
    for lab, px in (('S', 240), ('M', 380), ('L', 600), ('XL', 900), ('XXL', 1400))) + '</div>')

JS_FILTER = """
function applyFilter(){
  const st=document.getElementById('f_stage').value, sv=document.getElementById('f_status').value,
        et=document.getElementById('f_err').value, tx=document.getElementById('f_id').value.toLowerCase();
  let n=0;
  document.querySelectorAll('[data-row]').forEach(r=>{
    let show=true;
    if(tx && !r.dataset.id.toLowerCase().includes(tx)) show=false;
    if(show && st && sv){ if((r.dataset['s_'+st.toLowerCase()]||'not_run')!==sv) show=false; }  // data-* names are lowercased by HTML
    else if(show && sv){ let any=false; for(const k in r.dataset){ if(k.startsWith('s_') && r.dataset[k]===sv) any=true;} if(!any) show=false; }
    if(show && et){ if(!(r.dataset.err||'').split('|').includes(et)) show=false; }
    r.classList.toggle('hidden',!show); if(show) n++;
  });
  document.getElementById('f_count').textContent=n+' shown';
}
document.addEventListener('DOMContentLoaded',()=>{document.querySelectorAll('.filters select,.filters input').forEach(e=>e.addEventListener('input',applyFilter));applyFilter();});
"""


def _isnan(v):
    return v is None or (isinstance(v, float) and v != v)


def esc(s):
    return '' if _isnan(s) else html.escape(str(s))


def badge(status, text=None):
    return f'<span class="badge" style="background:{STATUS_COLOR.get(status, "#999")}">{esc(text or status)}</span>'


def thumb(png_abs, thumb_abs):
    """make/refresh a JPEG thumbnail, return True if it exists afterwards"""
    try:
        from PIL import Image
        if os.path.exists(thumb_abs) and os.path.getmtime(thumb_abs) >= os.path.getmtime(png_abs):
            return True
        MT.mkdirs_r(os.path.dirname(thumb_abs))
        im = Image.open(png_abs).convert('RGB')
        im.thumbnail((THUMB_W, THUMB_W * 3))
        im.save(thumb_abs, 'JPEG', quality=80)
        return True
    except Exception as e:
        print('thumb failed', png_abs, e)
        return False


def log_tail(batch_key, stage, ID):
    p = log_file(batch_key, stage, ID)
    if not p.exists():
        return ''
    lines = p.read_text(errors='replace').splitlines()
    return '\n'.join(lines[-LOG_TAIL:])


# ----------------------------------------------------------------------------- data assembly
def collect(batch_key):
    P = paths_for(batch_key)
    tracks = pd.read_csv(P.batch_work + 'tracks.csv') if os.path.exists(P.batch_work + 'tracks.csv') else pd.DataFrame()
    st = iter_status(batch_key)
    recs = {}
    for r in st.to_dict('records'):
        recs[(r['stage'], r['ID'])] = r
    ids = sorted(tracks[tracks.selected].ID.unique().tolist()) if len(tracks) else sorted({k[1] for k in recs if k[0] in TRACK_STAGES})
    batch_jobs = sorted({k[1] for k in recs if k[0] in ('B00', 'B01') and k[1] not in ids})

    rows = []
    for ID in ids:
        row = {'ID': ID}
        errs = set()
        for s in GRID_STAGES:
            r = recs.get((s, ID))
            row[s] = r['status'] if r else 'not_run'
            reason = (r or {}).get('reason')
            emsg = (r or {}).get('error_msg')
            row[s + '_reason'] = ('' if _isnan(reason) else reason) or ('' if _isnan(emsg) else emsg) or ''
            if r and isinstance(r.get('error_type'), str):
                errs.add(r['error_type'])
        row['errors'] = sorted(errs)
        if len(tracks):
            t = tracks[tracks.ID == ID].iloc[0]
            row.update(date=t.date, rgt=int(t.rgt), cycle=int(t.cycle), chunk=int(t.chunk), granule=t.granule)
        rows.append(row)
    return P, tracks, recs, ids, batch_jobs, rows


def summary_counts(rows):
    out = {}
    for s in GRID_STAGES:
        c = {k: 0 for k in STATUS_ORDER}
        for r in rows:
            c[r[s]] = c.get(r[s], 0) + 1
        out[s] = c
    return out


# ----------------------------------------------------------------------------- pages
def figures_of(P, ID):
    """all png figures of a track, grouped by stage prefix, relative to the track dir"""
    d = P.plot_batch + ID + '/'
    figs = sorted(str(Path(p).relative_to(d)) for p in glob.glob(d + '**/*.png', recursive=True))
    groups = {}
    for f in figs:
        stage = os.path.basename(f)[:3]
        if stage not in STAGES:
            stage = f.split('/')[0][:3] if f.split('/')[0][:3] in STAGES else 'other'
        groups.setdefault(stage, []).append(f)
    return groups


def track_page(P, batch_key, ID, recs, row):
    rel_track = '../../' + ID + '/'
    figs = figures_of(P, ID)
    parts = [f'<html><head><meta charset="utf-8"><title>{ID}</title><style>{CSS}</style><script>{JS_SIZE}</script></head><body>',
             f'<a href="../../index.html">&larr; {batch_key}</a>', SIZEBAR,
             f'<h1>{ID}</h1><div class="small">rgt {row.get("rgt")} cycle {row.get("cycle")} date {row.get("date")} '
             f'chunk {row.get("chunk")} granule {esc(row.get("granule"))}</div>',
             '<p>' + ' '.join(f'<a href="#{s}">{badge(row[s], s + ": " + row[s])}</a>' for s in GRID_STAGES) + '</p>']
    for s in GRID_STAGES:
        r = recs.get((s, ID))
        parts.append(f'<h2 id="{s}">{s} {badge(row[s])}</h2>')
        if r:
            parts.append(f'<div class="small">{esc(r.get("t_start"))} &rarr; {esc(r.get("t_end"))}  '
                         f'{esc(r.get("runtime_s"))} s  host {esc(r.get("host"))}  params {esc(r.get("params_hash"))}  '
                         f'git {esc(r.get("git_hash"))}</div>')
            if not _isnan(r.get('reason')) and r.get('reason'):
                parts.append(f'<p><b>reason:</b> {esc(r["reason"])}</p>')
            if isinstance(r.get('error_type'), str):
                parts.append(f'<p><b>{esc(r["error_type"])}:</b> {esc(r.get("error_msg"))}</p>')
            if r.get('info'):
                parts.append('<pre>' + esc(json.dumps(r['info'], indent=1, default=str)) + '</pre>')
        for f in figs.get(s, []):
            th = f'../thumbs/{ID}/' + f.replace('/', '__') + '.jpg'
            ok = thumb(P.plot_batch + ID + '/' + f, P.plot_batch + '_gallery/thumbs/' + ID + '/' + f.replace('/', '__') + '.jpg')
            pdf = f[:-4] + '.pdf'
            pdf_link = f' <a href="{rel_track + pdf}">pdf</a>' if os.path.exists(P.plot_batch + ID + '/' + pdf) else ''
            parts.append(f'<div class="card"><a href="{rel_track + f}">'
                         + (f'<img src="{th}" data-thumb="{th}" data-full="{rel_track + f}" loading="lazy">' if ok else esc(f)) +
                         f'</a><div class="small">{esc(f)}{pdf_link}</div></div>')
        if r:
            tail = log_tail(batch_key, s, ID if s != 'B01' else f'chunk{row.get("chunk", 0)}')
            if tail:
                parts.append('<details><summary class="small">log tail</summary><pre>' + esc(tail) + '</pre></details>')
    parts.append('</body></html>')
    out = P.plot_batch + '_gallery/track/' + ID + '.html'
    MT.mkdirs_r(os.path.dirname(out))
    Path(out).write_text('\n'.join(parts))


def stage_page(P, batch_key, stage, rows, recs):
    pattern = KEY_FIGURE.get(stage)
    parts = [f'<html><head><meta charset="utf-8"><title>{batch_key} {stage}</title><style>{CSS}</style>'
             f'<script>{JS_FILTER}</script><script>{JS_SIZE}</script><script>{JS_SORT}</script></head><body>',
             f'<a href="../../index.html">&larr; {batch_key}</a><h1>{stage} — {esc(pattern or "no key figure")}</h1>', SIZEBAR + SORTBAR,
             '<div class="filters">status <select id="f_status"><option value="">all</option>'
             + ''.join(f'<option>{s}</option>' for s in STATUS_ORDER) + '</select>'
             '<input type="hidden" id="f_stage" value="' + stage + '"><input type="hidden" id="f_err" value="">'
             'ID <input id="f_id" placeholder="filter"> <span id="f_count"></span></div><div class="grid">']
    for row in rows:
        ID = row['ID']
        r = recs.get((stage, ID))
        fig = None
        if pattern:
            hits = sorted(glob.glob(P.plot_batch + ID + '/' + pattern))
            fig = str(Path(hits[0]).relative_to(P.plot_batch + ID + '/')) if hits else None
        img = ''
        if fig:
            th_rel = '../thumbs/' + ID + '/' + fig.replace('/', '__') + '.jpg'
            if thumb(P.plot_batch + ID + '/' + fig, P.plot_batch + '_gallery/thumbs/' + ID + '/' + fig.replace('/', '__') + '.jpg'):
                img = f'<a href="../../{ID}/{fig}"><img src="{th_rel}" data-thumb="{th_rel}" data-full="../../{ID}/{fig}" loading="lazy"></a>'
        note = esc(row.get(stage + '_reason', ''))[:120]
        mtime = esc((r or {}).get('t_end') or (r or {}).get('t_start') or '')
        parts.append(f'<div class="card" data-row="1" data-id="{ID}" data-s_{stage}="{row[stage]}" data-mtime="{mtime}">'
                     f'<div><a href="../track/{ID}.html">{ID}</a> {badge(row[stage])} '
                     f'<span class="small">{esc((r or {}).get("runtime_s") or "")} s · {mtime[:16].replace("T", " ")}</span></div>{img}'
                     f'<div class="small">{note}</div></div>')
    parts.append('</div></body></html>')
    out = P.plot_batch + '_gallery/stage/' + stage + '.html'
    MT.mkdirs_r(os.path.dirname(out))
    Path(out).write_text('\n'.join(parts))


def index_page(P, batch_key, rows, recs, batch_jobs, counts):
    batch = {}
    if os.path.exists(P.batch_work + 'batch.json'):
        batch = json.load(open(P.batch_work + 'batch.json'))
    errs = sorted({e for r in rows for e in r['errors']})
    parts = [f'<html><head><meta charset="utf-8"><title>{batch_key}</title><style>{CSS}</style>'
             f'<script>{JS_FILTER}</script></head><body>',
             f'<h1>{batch_key}</h1><div class="small">{esc(batch.get("batch", {}).get("description", ""))} — '
             f'region lat {esc(batch.get("region", {}).get("lat"))} lon {esc(batch.get("region", {}).get("lon"))} — '
             f'time {esc(batch.get("time", {}).get("t0"))} .. {esc(batch.get("time", {}).get("t1"))} — '
             f'generated {dt.datetime.now():%Y-%m-%d %H:%M}</div>']
    if os.path.exists(P.plot_batch + '_batch/B00_overview.png'):
        parts.append('<p><a href="_batch/B00_overview.png"><img src="_batch/B00_overview.png" style="height:220px"></a></p>')

    # counts table
    parts.append('<h2>Stage summary</h2><table><tr><th>stage</th>' + ''.join(f'<th>{s}</th>' for s in STATUS_ORDER) + '</tr>')
    for s in GRID_STAGES:
        parts.append(f'<tr><th><a href="_gallery/stage/{s}.html">{s}</a></th>'
                     + ''.join(f'<td>{counts[s][k] or ""}</td>' for k in STATUS_ORDER) + '</tr>')
    parts.append('</table>')

    # batch-level jobs
    if batch_jobs:
        parts.append('<h2>Batch jobs</h2><table><tr><th>stage</th><th>job</th><th>status</th><th>runtime</th><th>info</th></tr>')
        for j in batch_jobs:
            for s in ('B00', 'B01'):
                r = recs.get((s, j))
                if r:
                    parts.append(f'<tr><td>{s}</td><td>{esc(j)}</td><td>{badge(r["status"])}</td><td>{esc(r.get("runtime_s"))} s</td>'
                                 f'<td class="small">{esc(json.dumps(r.get("info", {}), default=str))[:300]}</td></tr>')
        parts.append('</table>')

    # track grid
    parts.append('<h2>Tracks</h2><div class="filters">stage <select id="f_stage"><option value="">any</option>'
                 + ''.join(f'<option>{s}</option>' for s in GRID_STAGES) + '</select>'
                 'status <select id="f_status"><option value="">all</option>' + ''.join(f'<option>{s}</option>' for s in STATUS_ORDER) + '</select>'
                 'error <select id="f_err"><option value="">all</option>' + ''.join(f'<option>{esc(e)}</option>' for e in errs) + '</select>'
                 'ID <input id="f_id" placeholder="filter"> <span id="f_count"></span></div>')
    parts.append('<table><tr><th>track</th><th>date</th><th>rgt</th>' + ''.join(f'<th><a href="_gallery/stage/{s}.html">{s}</a></th>' for s in GRID_STAGES) + '</tr>')
    for row in rows:
        data = ' '.join(f'data-s_{s}="{row[s]}"' for s in GRID_STAGES)
        parts.append(f'<tr data-row="1" data-id="{row["ID"]}" data-err="{esc("|".join(row["errors"]))}" {data}>'
                     f'<td><a href="_gallery/track/{row["ID"]}.html">{row["ID"]}</a></td><td>{esc(row.get("date", ""))}</td><td>{esc(row.get("rgt", ""))}</td>')
        for s in GRID_STAGES:
            title = esc(row.get(s + '_reason', ''))
            r = recs.get((s, row['ID']))
            rt = f' {r["runtime_s"]:.0f}s' if r and not _isnan(r.get('runtime_s')) else ''
            parts.append(f'<td><span class="cell" style="background:{STATUS_COLOR[row[s]]}" title="{title}">'
                         f'<a href="_gallery/track/{row["ID"]}.html#{s}">{row[s]}{rt}</a></span></td>')
        parts.append('</tr>')
    parts.append('</table></body></html>')
    Path(P.plot_batch + 'index.html').write_text('\n'.join(parts))


# ----------------------------------------------------------------------------- main
def build(batch_key):
    P, tracks, recs, ids, batch_jobs, rows = collect(batch_key)
    MT.mkdirs_r(P.plot_batch + '_gallery/')
    counts = summary_counts(rows)
    for row in rows:
        track_page(P, batch_key, row['ID'], recs, row)
    for s in GRID_STAGES:
        stage_page(P, batch_key, s, rows, recs)
    index_page(P, batch_key, rows, recs, batch_jobs, counts)

    st = iter_status(batch_key)
    if len(st):
        st.drop(columns=['info', 'figures', 'outputs']).to_csv(P.plot_batch + 'status.csv', index=False)
    with open(P.plot_batch + 'status_summary.json', 'w') as f:
        json.dump({'batch': batch_key, 'generated': dt.datetime.now().isoformat(), 'n_tracks': len(rows), 'counts': counts}, f, indent=2)
    print(f'gallery: {len(rows)} tracks -> {P.plot_batch}index.html')
    print(pd.DataFrame(counts).T[STATUS_ORDER].to_string())
    return P.plot_batch + 'index.html'


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('-')]
    build(args[0] if args else 'SH_dev_small')
