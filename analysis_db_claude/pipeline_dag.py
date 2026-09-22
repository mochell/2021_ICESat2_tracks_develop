"""
Stage graph shared by the Snakefile, the status layer and the gallery.
"""

STAGES = ['B00', 'B01', 'A02', 'B02', 'B03', 'B04', 'B05', 'B06', 'C01']

# upstream stages whose status must be 'success' before a stage runs
DEPS = {
    'B00': [],
    'B01': ['B00'],
    'A02': ['B01'],
    'B02': ['B01'],
    'B03': ['B02'],
    'B04': ['B02', 'A02'],
    'B05': ['B04'],
    'B06': ['B02', 'B05'],
    'C01': ['B06', 'B05', 'A02'],
}

# per-track stages (B00 and B01 run per batch / per chunk)
TRACK_STAGES = ['A02', 'B02', 'B03', 'B04', 'B05', 'B06', 'C01']

SCRIPT = {
    'B00': 'stages/B00_discover.py',
    'B01': 'stages/B01_sliderule_load.py',
    'A02': 'stages/A02_ww3_prior.py',
    'B02': 'stages/B02_spectra.py',
    'B03': 'stages/B03_plot_spectra.py',
    'B04': 'stages/B04_angle.py',
    'B05': 'stages/B05_define_angle.py',
    'B06': 'stages/B06_correct.py',
    'C01': 'stages/C01_collect.py',
}

# work/<batch>/<dir>/ where a stage writes its data products (None: figures only)
STAGE_DIR = {
    'B00': None,
    'B01': 'B01_regrid',
    'A02': 'A02_prior',
    'B02': 'B02_spectra',
    'B03': None,
    'B04': 'B04_angle',
    'B05': 'B04_angle',
    'B06': 'B06_corrected_separated',
    'C01': 'C01_database',
}

# BLAS/OpenMP threads per job
THREADS = {'B02': 4, 'B04': 2, 'B06': 2}
THREADS_DEFAULT = 1

# scheduler resources (network concurrency caps); 'heavy' counts against the cerberus job limit
RESOURCES = {
    'B01': {'sliderule': 1, 'heavy': 1},
    'A02': {'thredds': 1},
    'B02': {'heavy': 1},
    'B04': {'heavy': 1},
    'B06': {'heavy': 1},
}

# one representative figure per stage (relative to plots/<hemis>/<batch>/<ID>/), glob allowed
KEY_FIGURE = {
    'B01': 'B01b_beam_statistics.png',
    'A02': 'A02_hindcast_prior.png',
    'B02': None,
    'B03': 'B03_specs_L*.png',
    'B04': 'B04_marginal_distributions.png',
    'B05': 'B05_dir_ov.png',
    'B06': 'B06_correction/B06_atten_ov_simple.png',
    'C01': None,
}

STATUS_ORDER = ['success', 'skip', 'blocked', 'fail', 'running', 'not_run']
STATUS_COLOR = {
    'success': '#4caf50', 'skip': '#ffb300', 'blocked': '#9e9e9e',
    'fail': '#e53935', 'running': '#1e88e5', 'not_run': '#eeeeee',
}


def topo_order(stages):
    """return the given stages in dependency order"""
    return [s for s in STAGES if s in set(stages)]
