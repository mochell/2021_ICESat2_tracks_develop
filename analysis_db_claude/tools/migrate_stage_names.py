"""
Rename the pipeline stages to the block order (2026-10-10):

    B block (data, spectra):        B00 B01 B02 B03                (unchanged)
    A block (angle, optional):      A02 -> A01a, B04 -> A02, B05 -> A03
    C block (correction, database): B06 -> C01, C01 -> C02
    B01 metadata folder/prefix:     A01b_ID -> B01_ID

One regex pass with a lookup table, so A02 -> A01a and B04 -> A02 never chain. A name is only
replaced as a whole token (not inside 'A02c' or 'ATL03') and never after 'analysis_db/' (comments
pointing at the original scripts keep their names).

    python tools/migrate_stage_names.py <batch_key> [...]   # migrate work/<batch>/ and plots/<hemis>/<batch>/
    python tools/migrate_stage_names.py --dry-run <batch>   # list what would change
    rename_text(s) / rename_name(s) are also used to rename the code itself.

Not idempotent by itself (A02 and C01 are old *and* new names), so every migrated root gets a
marker file MARKER and is refused on a second run. Renames inside a directory go through temporary
names first (B06 -> C01 while C01 -> C02 in the same folder).
"""
import os
import re
import sys
from pathlib import Path

MAP = {'A01b_ID': 'B01_ID', 'A02': 'A01a', 'B04': 'A02', 'B05': 'A03', 'B06': 'C01', 'C01': 'C02'}
PATTERN = re.compile(r'(?<![A-Za-z0-9])(?<!analysis_db/)(A01b_ID|A02|B04|B05|B06|C01)(?![a-z0-9])')
TEXT_SUFFIXES = {'.json', '.csv', '.html', '.txt', '.log', '.ok', ''}
MARKER = '.stage_names_2026-10-10'
TMP = '__stage_migration__'


def rename_text(s):
    return PATTERN.sub(lambda m: MAP[m.group(1)], s)


def rename_name(name):
    return rename_text(name)


def migrate_tree(root, dry=False):
    """rename files/dirs below root (deepest first) and rewrite stage names inside text files"""
    root = Path(root)
    if not root.exists():
        return 0, 0
    if (root / MARKER).exists():
        print(f'  {root}: {(root / MARKER).read_text().strip()} ({MARKER}), skipped')
        return 0, 0
    if not dry:                                          # written first: an interrupted run is never repeated
        (root / MARKER).write_text('stage names migration IN PROGRESS (interrupted runs need manual repair)\n')
    n_text = n_names = 0
    for dirpath, dirnames, filenames in os.walk(root, topdown=False):
        filenames = [f for f in filenames if f != MARKER]
        for f in filenames:
            p = Path(dirpath) / f
            if p.suffix in TEXT_SUFFIXES and p.stat().st_size < 20e6:
                try:
                    s = p.read_text()
                except UnicodeDecodeError:
                    s = None
                if s is not None:
                    t = rename_text(s)
                    if t != s:
                        n_text += 1
                        if not dry:
                            st = p.stat()
                            p.write_text(t)
                            os.utime(p, ns=(st.st_atime_ns, st.st_mtime_ns))   # keep mtimes (Snakemake)
        moves = [(name, rename_name(name)) for name in filenames + dirnames if rename_name(name) != name]
        n_names += len(moves)
        if dry:
            for name, new in moves:
                print(f'  {Path(dirpath) / name} -> {new}')
            continue
        for name, new in moves:                          # phase 1: free every old name
            (Path(dirpath) / name).rename(Path(dirpath) / (TMP + new))
        for name, new in moves:                          # phase 2: final names
            dst = Path(dirpath) / new
            if dst.exists():
                raise FileExistsError(f'{dst} exists, refusing to overwrite')
            (Path(dirpath) / (TMP + new)).rename(dst)
    if not dry:
        (root / MARKER).write_text('stage names migrated: ' + ', '.join(f'{k} -> {v}' for k, v in MAP.items()) + '\n')
    return n_text, n_names


def main(argv):
    dry = '--dry-run' in argv
    batches = [a for a in argv if not a.startswith('-')]
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from pipeline_config import paths_for
    for b in batches:
        P = paths_for(b)
        for root in (P.batch_work, P.plot_batch):
            n_text, n_names = migrate_tree(root, dry)
            print(f'{b}: {root}: {n_text} text files rewritten, {n_names} names renamed{" (dry run)" if dry else ""}')


if __name__ == '__main__':
    main(sys.argv[1:])
