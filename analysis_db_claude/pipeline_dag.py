"""
Stage graph shared by the Snakefile, the status layer and the gallery.
"""

# blocks: B = data and spectra, A = wave angle (optional for C), C = correction and database
STAGES = ['B00', 'B01', 'B02', 'B03', 'A01a', 'A02', 'A03', 'C01', 'C02']

# upstream stages whose status must be 'success' before a stage runs
DEPS = {
    'B00': [],
    'B01': ['B00'],
    'B02': ['B01'],
    'B03': ['B02'],
    'A01a': ['B01'],             # WW3 prior: needs the B01 track position/time only
    'A02': ['B02', 'A01a'],
    'A03': ['A02'],
    'C01': ['B02'],              # angle (A02/A03) is optional, see ORDER_AFTER
    'C02': ['C01', 'A01a'],      # prior scalars from A01a; angle optional
}

# stages that must have *finished* (any status) before a stage runs, without being required to
# succeed: C01 / C02 use the A03 angle when it exists and otherwise write uncorrected products
ORDER_AFTER = {'C01': ['A02', 'A03'], 'C02': ['A03']}

# per-track stages (B00 and B01 run per batch / per chunk)
TRACK_STAGES = ['B02', 'B03', 'A01a', 'A02', 'A03', 'C01', 'C02']

SCRIPT = {
    'B00': 'stages/B00_discover.py',
    'B01': 'stages/B01_sliderule_load.py',
    'A01a': 'stages/A01a_ww3_prior.py',
    'B02': 'stages/B02_spectra.py',
    'B03': 'stages/B03_plot_spectra.py',
    'A02': 'stages/A02_angle.py',
    'A03': 'stages/A03_define_angle.py',
    'C01': 'stages/C01_correct.py',
    'C02': 'stages/C02_collect.py',
}

# work/<batch>/<dir>/ where a stage writes its data products (None: figures only)
STAGE_DIR = {
    'B00': None,
    'B01': 'B01_regrid',
    'A01a': 'A01a_prior',
    'B02': 'B02_spectra',
    'B03': None,
    'A02': 'A02_angle',
    'A03': 'A02_angle',
    'C01': 'C01_corrected_separated',
    'C02': 'C02_database',
}

# BLAS/OpenMP threads per job
THREADS = {'B02': 4, 'A02': 2, 'C01': 2}
THREADS_DEFAULT = 1

# scheduler resources (network concurrency caps); 'heavy' counts against the cerberus job limit
RESOURCES = {
    'B01': {'sliderule': 1, 'heavy': 1},
    'A01a': {'thredds': 1},
    'B02': {'heavy': 1},
    'A02': {'heavy': 1},
    'C01': {'heavy': 1},
}

# one representative figure per stage (relative to plots/<hemis>/<batch>/<ID>/), glob allowed
KEY_FIGURE = {
    'B01': 'B01b_beam_statistics.png',
    'A01a': 'A01a_hindcast_prior.png',
    'B02': None,
    'B03': 'B03_specs_L*.png',
    'A02': 'A02_marginal_distributions.png',
    'A03': 'A03_dir_ov.png',
    'C01': 'C01_correction/C01_atten_ov_simple.png',
    'C02': None,
}

STATUS_ORDER = ['success', 'skip', 'blocked', 'fail', 'running', 'not_run']
STATUS_COLOR = {
    'success': '#4caf50', 'skip': '#ffb300', 'blocked': '#9e9e9e',
    'fail': '#e53935', 'running': '#1e88e5', 'not_run': '#eeeeee',
}


def topo_order(stages):
    """return the given stages in dependency order"""
    return [s for s in STAGES if s in set(stages)]
