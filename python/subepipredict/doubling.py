"""Successive doubling times of cumulative incidence (``getDoublingTimeCurve.m``)."""

from __future__ import annotations

import numpy as np

from .rt import _colon


def getDoublingTimeCurve(curve, DT, print1=0):
    """Successive doubling times of the cumulative curve of an incidence curve.

    Port of ``getDoublingTimeCurve.m``.

    Parameters
    ----------
    curve : incidence at successive time points
    DT : time between points
    print1 : 1 to plot the cumulative curve and the doubling times

    Returns
    -------
    ``(tds, C0data, curve2, doublingtimes)``
        ``tds``: rows ``[time, cumulative count]`` at the start and at each
        doubling. ``C0data``: cumulative count at the first point (None for a
        curve of one point). ``curve2``: cumulative curve at the time points,
        ``[time, count]`` (empty when the count never doubles, as in MATLAB).
        ``doublingtimes``: time taken by each successive doubling (at most
        30); empty when the count does not double within the curve.

    The cumulative curve is interpolated on a grid of step 0.01 and the time
    of each doubling is the last grid time at which the curve is at or below
    twice the previous level.
    """
    dt1 = 0.01
    tds = np.zeros((0, 2))
    C0data = None
    curve2 = np.zeros((0, 2))
    doublingtimes = np.zeros(0)

    curve = np.asarray(curve, dtype=float).ravel()
    if curve.size <= 1:
        return tds, C0data, curve2, doublingtimes

    tb = np.arange(curve.size, dtype=float) * DT
    cumulative1 = np.cumsum(curve)
    cum1 = cumulative1[-1]
    C0data = float(cumulative1[0])

    XX = _colon(tb[0], dt1, tb[-1])
    cumulative1 = np.interp(XX, tb, cumulative1)
    index1 = cumulative1 <= cum1
    XX = XX[index1]
    cumulative1 = cumulative1[index1]
    step = int(round(1.0 / dt1))
    curveOut = np.column_stack([XX, cumulative1])[::step, :]

    if print1:
        import matplotlib.pyplot as plt

        plt.subplot(1, 2, 1)
        plt.plot(XX, cumulative1, "-")
        plt.xlabel("Time", fontsize=16)
        plt.ylabel("Cumulative cases", fontsize=16)

    x = cumulative1
    timevect = XX
    monotone = bool(np.all(np.diff(x) >= 0))
    C0 = x[0]
    rows = [(0.0, float(C0))]
    for _ in range(25):
        C0i1 = 2 * C0
        if monotone:
            last = int(np.searchsorted(x, C0i1, side="right")) - 1
        else:
            hits = np.flatnonzero(x <= C0i1)
            last = int(hits[-1]) if hits.size else -1
        if last >= 0:
            rows.append((float(timevect[last]), float(x[last])))
        C0 = C0i1
    tds = np.unique(np.array(rows, dtype=float), axis=0)  # unique(tds,'rows'): sorted rows
    tds = tds[:-1, :]

    if print1:
        import matplotlib.pyplot as plt

        plt.subplot(1, 2, 2)
        if tds.shape[0] > 0:
            plt.plot(tds[:, 1], np.concatenate(([tds[0, 0]], np.diff(tds[:, 0]))), "ro--")
        plt.ylabel("Doubling time", fontsize=16)
        plt.xlabel("Doubling cumulative incidence", fontsize=16)

    if tds.shape[0] <= 1:  # the count never doubles
        return tds, C0data, curve2, doublingtimes

    maxdoublings = 30
    doublingtimes = np.diff(tds[:, 0])
    doublingtimes = doublingtimes[: min(doublingtimes.size, maxdoublings)]
    return tds, C0data, curveOut, doublingtimes
