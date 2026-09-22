"""
Shared setup for the analysis_db_claude stages.

Replaces the old ``exec(open(os.environ['PYTHONSTARTUP']))`` / ``STARTUP_2021_IceSAT2`` pair:
picks the machine config, puts ``modules/`` on ``sys.path``, imports the plotting helpers and
applies the figure style. Every stage starts with::

    from pipeline_config import mconfig, np, pd, xr, plt, M, MT, col, paths_for, save_fig

Config file selection (first hit wins):
  1. ``$IS2_CONFIG``                      explicit override
  2. ``config/config.json``               cerberus (untracked)
  3. ``config/config_local.json``         laptop
"""
import os
import sys
import json
import socket
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]          # repository root
HERE = Path(__file__).resolve().parent               # analysis_db_claude/


def _config_file():
    if os.environ.get('IS2_CONFIG'):
        return Path(os.environ['IS2_CONFIG'])
    if (ROOT / 'config' / 'config.json').exists():
        return ROOT / 'config' / 'config.json'
    return ROOT / 'config' / 'config_local.json'


CONFIG_FILE = _config_file()
with open(CONFIG_FILE) as f:
    mconfig = json.load(f)
mconfig['paths'].setdefault('groundtracks', str(ROOT / 'data' / 'groundtracks') + '/')
HOST = socket.gethostname().split('.')[0]

# project modules (stages import e.g. `spicke_remover`, `generalized_FT` as top-level names)
for _p in (mconfig['paths']['local_script'], mconfig['paths']['local_script'] + '/ICEsat2_SI_tools/'):
    if _p not in sys.path:
        sys.path.insert(0, _p)
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import matplotlib
matplotlib.use(os.environ.get('MPLBACKEND', 'Agg'))
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import string

xr.set_options(display_style='text')

import m_colormanager_ph3 as M_color
import m_tools_ph3 as MT
import m_general_ph3 as M

col = M_color.color(path=str(ROOT / 'config') + '/', name='color_def')
lstrings = iter([i + ') ' for i in list(string.ascii_lowercase)])
fig_sizes = mconfig['fig_sizes']['AMS']


# ----------------------------------------------------------------------------- figure style
def _apply_style(small, medium):
    plt.rc('font', size=small, serif='Helvetica Neue', weight='normal')
    plt.rc('text', usetex='false')
    plt.rc('axes', titlesize=medium, labelweight='normal')
    plt.rc('axes', labelsize=small, labelweight='normal')
    plt.rc('xtick', labelsize=small)
    plt.rc('ytick', labelsize=small)
    plt.rc('legend', fontsize=small, frameon=False)
    plt.rc('figure', titlesize=medium, titleweight='bold', autolayout=True)


def setup_style():
    _apply_style(8, 10)
    plt.rc('path', simplify=True)
    plt.rcParams['figure.figsize'] = (10, 8)
    plt.rc('xtick.major', size=4, width=1)
    plt.rc('ytick.major', size=3.8, width=1)
    plt.rc('axes', labelsize=10, labelweight='normal')
    plt.rc('axes.spines', top=False, right=False)


def font_for_print():
    _apply_style(6, 8)


def font_for_pres():
    _apply_style(10, 12)


setup_style()


# ----------------------------------------------------------------------------- paths
class Paths:
    """All members are strings ending with '/'. Built by paths_for()."""

    def __init__(self, batch_key, ID=None):
        self.batch_key = batch_key
        self.ID = ID
        self.hemis = batch_key.split('_')[0]
        self.work = _slash(mconfig['paths']['work'])
        self.plot = _slash(mconfig['paths']['plot'])
        self.batch_work = self.work + batch_key + '/'
        self.status = self.batch_work + 'status/'
        self.logs = self.batch_work + 'logs/'
        self.params = self.batch_work + 'params/'
        self.plot_batch = self.plot + self.hemis + '/' + batch_key + '/'
        self.plot_track = self.plot_batch + ID + '/' if ID else None
        self.groundtracks = _slash(mconfig['paths']['groundtracks'])

    def stage_dir(self, name, mkdir=True):
        """work/<batch>/<name>/ e.g. stage_dir('B02_spectra')"""
        p = self.batch_work + name + '/'
        if mkdir:
            MT.mkdirs_r(p)
        return p

    def track_plot_dir(self, sub=None, mkdir=True):
        p = self.plot_track + (sub + '/' if sub else '')
        if mkdir:
            MT.mkdirs_r(p)
        return p

    def __repr__(self):
        return f'Paths(work={self.batch_work}, plots={self.plot_batch})'


def _slash(p):
    return p if p.endswith('/') else p + '/'


def paths_for(batch_key, ID=None):
    return Paths(batch_key, ID)


# ----------------------------------------------------------------------------- figures
def save_fig(F, path, name, pdf=True, png=True, close=True):
    """
    Save a m_general_ph3.figure_axis_xy figure (or a bare matplotlib Figure) as PNG and,
    optionally, PDF. Every figure gets a PNG so the gallery can show it.
    """
    MT.mkdirs_r(path)
    fig = F.fig if hasattr(F, 'fig') else F
    M.remove_empty_axes(fig)
    if png:
        fig.savefig(os.path.join(path, name + '.png'), bbox_inches='tight', format='png', dpi=180)
    if pdf:
        fig.savefig(os.path.join(path, name + '.pdf'), bbox_inches='tight', format='pdf', dpi=300)
    if close:
        plt.close(fig)


# ----------------------------------------------------------------------------- CLI
def cli_args(argv, default=None):
    """
    Return (ID, batch_key) from ``python stage.py <ID> <batch_key>``.
    With no arguments (interactive use) return ``default``.
    """
    args = [a for a in argv[1:] if not a.startswith('-')]
    if len(args) >= 2:
        return args[0], args[1]
    if default is None:
        raise SystemExit(f'usage: {os.path.basename(argv[0])} <ID> <batch_key>')
    print('no arguments given, using defaults:', default)
    return default


if __name__ == '__main__':
    print('config file :', CONFIG_FILE)
    print('host        :', HOST)
    print('work        :', mconfig['paths']['work'])
    print('plot        :', mconfig['paths']['plot'])
    print('groundtracks:', mconfig['paths']['groundtracks'])
    print('modules     :', mconfig['paths']['local_script'])
    print('color theme :', col.__class__.__name__, 'M:', M.__name__, 'MT:', MT.__name__)
    print(paths_for('SH_testSLsinglefile2', 'SH_20190502_05180312'))
