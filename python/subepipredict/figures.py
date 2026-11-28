"""Helpers shared by the plotting functions: panel grids, date axes, labels.

Ports of ``subepidemicPanelGrid.m``, ``subepidemicEstimateLine.m``,
``subepidemicParameterLabel.m`` and ``GetAdjustedFontSize.m``; the date ticks
stand in for the ``datestr`` / ``xticklabel_rotate`` blocks of the scripts.
"""

from __future__ import annotations

import math

import numpy as np

from .errors import SubEpiPredictError
from .utils import datenum, datestr, factor

# colors of the sub-epidemic components, in the order of the MATLAB scripts
COMPONENT_COLORS = ["r", "b", "g", "m", "c", "k", "y"]
RANKING_COLORS = ["r", "b", "g", "m", "c", "y"]
GRAY7 = (2 / 3, 2 / 3, 2 / 3)  # gray(10), row 7


def GetAdjustedFontSize():
    """Font size of the figure labels.

    MATLAB scales 24 points by the screen resolution; the figures of the port
    are sized in inches, where 13 points gives the same proportions.
    """
    return 13


def subepidemicPanelGrid(n):
    """Rows and columns of a figure with one panel per model.

    Port of ``subepidemicPanelGrid.m``: the grid built from the prime factors
    of ``n`` whenever that grid holds ``n`` panels (n rows by 1 column for a
    prime, factor 1 by factor 2 for two factors, factor 1 * factor 2 by
    factor 3 for three), otherwise the squarest grid that holds ``n`` (4 by 4
    for 16, 5 by 5 for 24, 6 by 6 for 32). Returns ``(rows, cols)``.
    """
    ok = isinstance(n, (int, float, np.integer, np.floating)) and not isinstance(n, (bool, np.bool_))
    if not (ok and math.isfinite(n) and n >= 1 and n == int(n)):
        raise SubEpiPredictError("SubEpiPredict:PanelGrid", "The number of panels must be a positive integer.")
    n = int(n)
    factors = factor(n)
    if len(factors) == 1:
        rows, cols = factors[0], 1
    elif len(factors) == 3:
        rows, cols = factors[0] * factors[1], factors[2]
    else:
        rows, cols = factors[0], factors[1]
    if rows * cols < n:
        cols = int(math.ceil(math.sqrt(n)))
        rows = int(math.ceil(n / cols))
    return rows, cols


def subepidemicEstimateLine(estimate, ax=None):
    """Vertical dashed line at the fitted estimate (``subepidemicEstimateLine.m``).

    Drawn over the full height of the axes, whose vertical limits are kept.
    Nothing is drawn, and None is returned, when ``estimate`` is not a finite
    real number.
    """
    ok = isinstance(estimate, (int, float, np.integer, np.floating)) and not isinstance(estimate, (bool, np.bool_))
    if not (ok and math.isfinite(estimate)):
        return None
    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()
    yl = ax.get_ylim()
    (h,) = ax.plot([estimate, estimate], yl, "k--", linewidth=2)
    ax.set_ylim(yl)
    return h


def _localNumber(x):
    """Three significant digits; plain notation for 0.001 <= |x| < 1e7."""
    x = float(x)
    if not math.isfinite(x):
        return "%g" % x
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    y = round(x / 10 ** (e - 2)) * 10 ** (e - 2)  # three significant digits
    if abs(y) >= 1e7 or abs(y) < 1e-3:
        return "%.3g" % y
    if abs(y) >= 1000:
        return "%d" % round(y)
    return "%.3g" % y


def subepidemicParameterLabel(name, estimate, lower, upper):
    """Title text for one parameter: ``'<name>=<estimate>(95%CI:<lower>,<upper>)'``.

    Port of ``subepidemicParameterLabel.m``. ``estimate`` is the value in the
    ranked fit; ``lower`` and ``upper`` are the 2.5% and 97.5% quantiles of
    the bootstrap refits. Each number is rounded to three significant digits
    and written without an exponent between 0.001 and 10 million.
    """
    if not isinstance(name, str):
        raise SubEpiPredictError("SubEpiPredict:InvalidLabel", "name must be a character row.")
    for v in (estimate, lower, upper):
        ok = isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_))
        if not ok:
            raise SubEpiPredictError(
                "SubEpiPredict:InvalidLabel", "estimate, lower and upper must be real scalars."
            )
    return "%s=%s(95%%CI:%s,%s)" % (name, _localNumber(estimate), _localNumber(lower), _localNumber(upper))


def mathLabel(label):
    """Parameter label with the index as a subscript for a figure title (``r_1=...`` -> ``$r_1$=...``)."""
    name, sep, rest = label.partition("=")
    if "_" in name:
        name = "$%s$" % name
    elif name == "alpha":
        name = r"$\alpha$"
    return name + sep + rest


def series_dates(datevecfirst1, timelags, DT, npoints):
    """Serial date numbers of the points of a series that starts ``timelags``
    time steps after ``datevecfirst1``."""
    datenumIni = datenum(datevecfirst1) + int(timelags) * int(DT)
    return datenumIni + int(DT) * np.arange(int(npoints))


def date_ticks(ax, datevecfirst1, timelags, DT, npoints, fontsize=9):
    """Date labels of the time axis, as the MATLAB scripts set them: every third
    point for daily data, every second for weekly data, the year for yearly data."""
    npoints = int(npoints)
    if DT == 365:
        years = int(datevecfirst1[0]) + int(timelags) + np.arange(npoints)
        ticks = np.arange(npoints)
        labels = [str(v) for v in years]
    else:
        step = 3 if DT == 1 else 2
        dn = series_dates(datevecfirst1, timelags, DT, npoints)
        ticks = np.arange(0, npoints, step)
        labels = [datestr(dn[k], "mm-dd") for k in ticks]
    # keep the labels readable in long series
    if len(ticks) > 40:
        keep = slice(None, None, int(math.ceil(len(ticks) / 40)))
        ticks, labels = ticks[keep], labels[keep]
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=90, fontsize=fontsize)


def plot_curves(ax, x, Y, color, linewidth=0.8, alpha=1.0):
    """Draw the columns of ``Y`` against ``x`` as one collection (fast for many curves)."""
    from matplotlib.collections import LineCollection

    x = np.asarray(x, dtype=float).ravel()
    Y = np.asarray(Y, dtype=float)
    if Y.ndim == 1:
        Y = Y.reshape(-1, 1)
    if Y.shape[1] == 0:
        return None
    segments = [np.column_stack([x, Y[:, j]]) for j in range(Y.shape[1])]
    lc = LineCollection(segments, colors=color, linewidths=linewidth, alpha=alpha)
    ax.add_collection(lc)
    ax.autoscale_view()
    return lc


def histogram(ax, values, bins=10, color=(0.24, 0.15, 0.66)):
    """``hist(values)``: a histogram of 10 equal bins; values that are equal up
    to rounding go into one bar."""
    values = np.asarray(values, dtype=float).ravel()
    values = values[np.isfinite(values)]
    if values.size == 0:
        return
    lo, hi = float(np.min(values)), float(np.max(values))
    if not (hi - lo > 1e-12 * max(1.0, abs(hi), abs(lo))):
        width = 1.0 if lo == 0 else abs(lo) * 0.1
        ax.bar([lo], [values.size], width=width, color=color, edgecolor="k", linewidth=0.5)
        ax.set_xlim(lo - 5.5 * width, lo + 5.5 * width)
        return
    counts, edges = np.histogram(values, bins=bins, range=(lo, hi))
    ax.bar(edges[:-1], counts, width=np.diff(edges), align="edge", color=color, edgecolor="k", linewidth=0.5)
