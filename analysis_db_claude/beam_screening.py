"""
Beam usability test, shared by B01 (decides which beams and tracks are kept) and B02 (gFT window).

A beam is usable when, after the photon filter (N_photos > N_photos_min, finite h_mean), it has at
least `min_dense_segments` dense gFT-length segments inside one window of `dense_window_km`. A segment
is dense when at least `segment_fill_min` of its dx slots hold a point; the gFT drops stancils below
the same fill (generalized_FT.wavenumber_spectrogram_gFT.min_fill), so a usable beam yields spectra.
The parameters live in the shared [gft] table of params/<v>.toml (copied into B01 and B02).
"""
import numpy as np


def stancil_length(T_max, L_factor, dx):
    """gFT window: L_factor wavelengths of the longest period T_max (deep water); (points, metres)"""
    k_0 = (2 * np.pi / T_max) ** 2 / 9.81
    Lpoints = int(np.round(2 * np.pi / k_0 / dx) * L_factor)
    return Lpoints, Lpoints * dx


def dense_segments(x, xlims, Lmeters, dx, fill_min, window_m):
    """
    Largest number of dense L-segments inside any window of window_m along the track.
    [xlims] is cut into consecutive segments of Lmeters; a segment is dense when at least fill_min
    of its dx slots hold a point.
    """
    edges = np.arange(xlims[0], xlims[1] + Lmeters, Lmeters)
    if edges.size < 2 or x.size == 0:
        return 0
    counts, _ = np.histogram(x, edges)
    dense = (counts / (Lmeters / dx) >= fill_min).astype(int)
    n_win = max(int(window_m // Lmeters), 1)
    return int(np.convolve(dense, np.ones(n_win, dtype=int), mode='valid').max()) if dense.size >= n_win else int(dense.sum())


def dense_segments_per_beam(Ti, beams, gft, dx):
    """{beam: dense segments in the best window} for per-beam tables with x, h_mean, N_photos"""
    _, Lmeters = stancil_length(gft['T_max'], gft['L_factor'], dx)
    out = {}
    for b in beams:
        d = Ti[b]
        m = (d['N_photos'] > gft['N_photos_min']) & np.isfinite(d['h_mean'])
        x = np.asarray(d['x'][m], dtype=float)
        out[b] = 0 if x.size == 0 else dense_segments(x, (x.min(), x.max()), Lmeters, dx,
                                                      gft['segment_fill_min'], gft['dense_window_km'] * 1e3)
    return out
