"""Small helpers shared by the toolbox: MATLAB-compatible statistics, dates,
the random-number stream and the data files.
"""

from __future__ import annotations

import datetime
import math
import os

import numpy as np

from .errors import SubEpiPredictError

# <============================================================================>
# <============================ Struct container ==============================>
# <============================================================================>


class Struct(dict):
    """A dict whose keys can also be read and written as attributes.

    Stands in for a MATLAB struct: ``s.LB`` and ``s['LB']`` are the same thing.
    """

    __slots__ = ()

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name) from None

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        try:
            del self[name]
        except KeyError:
            raise AttributeError(name) from None

    def copy(self):
        return Struct(self)


# <============================================================================>
# <========================= Random-number stream =============================>
# <============================================================================>

_GLOBAL_RNG = np.random.Generator(np.random.MT19937())


def rng(seed=None):
    """Seed the toolbox's random-number stream, like MATLAB's ``rng(seed)``.

    Every function that draws random numbers (start points of the search,
    bootstrap data, forecast noise) uses this stream unless it is given its
    own ``numpy.random.Generator`` through its ``rng`` argument.

    The stream is a Mersenne Twister, as in MATLAB, but NumPy and MATLAB turn
    its output into normal, Poisson and negative binomial variates with
    different algorithms. A seeded Python run therefore repeats itself, not
    the MATLAB run with the same seed. Returns the generator.
    """
    global _GLOBAL_RNG
    _GLOBAL_RNG = np.random.Generator(np.random.MT19937(seed))
    return _GLOBAL_RNG


def get_rng(rng_=None):
    """The generator to draw from: the one given, or the toolbox stream."""
    if rng_ is None:
        return _GLOBAL_RNG
    if isinstance(rng_, np.random.Generator):
        return rng_
    return np.random.Generator(np.random.MT19937(rng_))


# <============================================================================>
# <====================== MATLAB-compatible statistics ========================>
# <============================================================================>


def quantile(x, p, axis=None):
    """Quantiles with MATLAB's definition (``quantile`` / ``prctile``).

    MATLAB treats the sorted sample ``x(1..n)`` as the ``(k-0.5)/n`` quantiles
    and interpolates linearly between them, using the minimum and maximum
    outside that range. This is NumPy's ``method='hazen'``; NumPy's default
    (``'linear'``) gives slightly different interval limits. NaNs are ignored,
    as in MATLAB.

    ``axis=None`` follows MATLAB for a vector (one value per probability). For
    a matrix pass the axis explicitly: MATLAB's ``quantile(X, p)`` works down
    the columns (``axis=0``), ``quantile(X, p, 2)`` along the rows (``axis=1``).
    """
    import warnings

    x = np.asarray(x, dtype=float)
    scalar_p = np.ndim(p) == 0
    p = np.atleast_1d(np.asarray(p, dtype=float))
    if x.size == 0:
        out = np.full(p.shape, np.nan)
        return float(out[0]) if scalar_p else out
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        out = np.nanquantile(x, p, axis=axis, method="hazen")
    if scalar_p:
        out = out[0]
        if np.ndim(out) == 0:
            return float(out)
    return out


def movmean(y, span):
    """Moving average like MATLAB's ``smoothdata(y, 'movmean', span)``.

    The window is centred on each point. An odd ``span`` covers
    ``(span-1)/2`` points on each side; an even one covers ``span/2`` points
    before the current point and ``span/2 - 1`` after it. Near the ends the
    window is cut at the end of the series and the mean is taken over the
    points that are left. A ``span`` that is not a whole number is treated as
    MATLAB does: the window is ``[t - span/2, t + span/2)`` around each point.
    """
    y = np.asarray(y, dtype=float).ravel()
    n = y.size
    if n == 0 or span <= 1:
        return y.copy()
    span = float(span)
    if span == int(span) and int(span) % 2 == 1:
        back = fwd = (int(span) - 1) // 2
    elif span == int(span):
        back = int(span) // 2
        fwd = int(span) // 2 - 1
    else:
        back = int(math.floor(span / 2))
        fwd = int(math.ceil(span / 2)) - 1
    c = np.concatenate(([0.0], np.cumsum(y)))
    lo = np.maximum(np.arange(n) - back, 0)
    hi = np.minimum(np.arange(n) + fwd, n - 1) + 1
    return (c[hi] - c[lo]) / (hi - lo)


def factor(n):
    """Prime factors of a positive integer in increasing order (MATLAB ``factor``)."""
    n = int(n)
    if n < 4:
        return [n]
    out = []
    f = 2
    while f * f <= n:
        while n % f == 0:
            out.append(f)
            n //= f
        f += 1
    if n > 1:
        out.append(n)
    return out


# <============================================================================>
# <=================================== Text ===================================>
# <============================================================================>


def num2str(x, precision=None):
    """Text of a number as MATLAB's ``num2str`` writes it.

    Whole numbers are written without a decimal point. Other numbers get
    ``max(floor(log10(|x|)) + 5, 5)`` significant digits, which is four
    decimals for numbers above 1 (``num2str(1140.55357)`` is ``'1140.5536'``).
    ``num2str(x, n)`` writes ``n`` significant digits.
    """
    if isinstance(x, (bool, np.bool_)):
        return "1" if x else "0"
    xf = float(x)
    if precision is not None:
        return "%.*g" % (int(precision), xf)
    if math.isnan(xf):
        return "NaN"
    if math.isinf(xf):
        return "Inf" if xf > 0 else "-Inf"
    if xf == int(xf) and abs(xf) < 1e15:
        return str(int(xf))
    digits = min(16, max(int(math.floor(math.log10(abs(xf)))) + 5, 5))
    return "%.*g" % (digits, xf)


def mat2str(v, precision=6):
    """Text of a number or a vector like MATLAB's ``mat2str(v, precision)``."""
    a = np.asarray(v, dtype=float)
    if a.ndim == 0:
        return "%.*g" % (precision, float(a))
    if a.size == 0:
        return "zeros(%d,%d)" % (a.shape[0] if a.ndim > 1 else 0, a.shape[-1] if a.ndim > 1 else 0)
    if a.ndim == 1:
        return "[" + " ".join("%.*g" % (precision, x) for x in a) + "]"
    return "[" + ";".join(" ".join("%.*g" % (precision, x) for x in row) for row in a) + "]"


def num2ordinal(num1):
    """Ordinal text of a positive integer: 1st, 2nd, 3rd, 4th, ..., 11th, 21st.

    Port of ``num2ordinal.m``; used in the figure titles of the ranked models.
    """
    ok = isinstance(num1, (int, float, np.integer, np.floating)) and not isinstance(num1, bool)
    if not (ok and math.isfinite(num1) and num1 >= 1 and num1 == int(num1)):
        raise SubEpiPredictError("SubEpiPredict:InvalidRank", "num2ordinal needs one positive integer.")
    num1 = int(num1)
    lastTwo = num1 % 100
    if 11 <= lastTwo <= 13:
        suffix = "th"  # 11th, 12th, 13th, 111th
    else:
        suffix = {1: "st", 2: "nd", 3: "rd"}.get(num1 % 10, "th")
    return "%d%s" % (num1, suffix)


# <============================================================================>
# <================================== Dates ===================================>
# <============================================================================>

_DATENUM_OFFSET = 366  # MATLAB datenum of a date = its proleptic ordinal + 366


def datenum(*args):
    """MATLAB serial date number.

    ``datenum([y, m, d])``, ``datenum(y, m, d)`` or ``datenum('mm-dd-yyyy')``.
    Day numbers outside the month roll over, as in MATLAB.
    """
    if len(args) == 1 and isinstance(args[0], str):
        text = args[0].strip()
        month, day, year = int(text[0:2]), int(text[3:5]), int(text[6:])
        if year < 100:
            year += 2000
    elif len(args) == 1:
        year, month, day = [int(v) for v in list(args[0])[:3]]
    else:
        year, month, day = [int(v) for v in args[:3]]
    year += (month - 1) // 12
    month = (month - 1) % 12 + 1
    return datetime.date(year, month, 1).toordinal() + (day - 1) + _DATENUM_OFFSET


def datevec(dn):
    """``[year, month, day]`` of a MATLAB serial date number."""
    d = datetime.date.fromordinal(int(round(float(dn))) - _DATENUM_OFFSET)
    return [d.year, d.month, d.day]


def datestr(dn, fmt="mm-dd-yyyy"):
    """Text of a serial date number in one of the formats the toolbox uses:
    ``'mm-dd-yyyy'``, ``'mm-dd-yy'`` and ``'mm-dd'``."""
    y, m, d = datevec(dn)
    if fmt == "mm-dd-yyyy":
        return "%02d-%02d-%04d" % (m, d, y)
    if fmt == "mm-dd-yy":
        return "%02d-%02d-%02d" % (m, d, y % 100)
    if fmt == "mm-dd":
        return "%02d-%02d" % (m, d)
    raise ValueError("Unsupported date format %r." % fmt)


def cadtemporal_of(DT):
    """'daily', 'weekly' or 'yearly' for DT = 1, 7 or 365."""
    try:
        return {1: "daily", 7: "weekly", 365: "yearly"}[int(DT)]
    except (KeyError, TypeError, ValueError):
        raise SubEpiPredictError(
            "SubEpiPredict:InvalidSettings", "DT = %s must be 1 (daily), 7 (weekly) or 365 (yearly)." % (DT,)
        ) from None


# <============================================================================>
# <================================ Data files ================================>
# <============================================================================>

_PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

# Where the example datasets are looked for when a file is not in ./input:
# a folder 'input' inside the package (present only if the data are bundled
# with an installed copy), then the MATLAB code folder of the repository this
# 'python' folder belongs to.
_BUNDLED_INPUT_DIRS = (
    os.path.join(_PACKAGE_DIR, "input"),
    os.path.join(os.path.dirname(os.path.dirname(_PACKAGE_DIR)), "ensemble n-subepidemic code v1.0", "input"),
)


def example_input_dir():
    """Folder with the example datasets of the toolbox, or None when it is not found.

    In a checkout of the repository this is ``ensemble n-subepidemic code v1.0/input``.
    """
    for folder in _BUNDLED_INPUT_DIRS:
        if os.path.isdir(folder):
            return folder
    return None


def find_input_file(name, input_dir="./input"):
    """Full path of a data file, or None: first in ``input_dir``, then among the
    example datasets of the repository."""
    for folder in (input_dir,) + _BUNDLED_INPUT_DIRS:
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            return path
    return None


def resolve_input_file(name, input_dir="./input"):
    """Full path of a data file; an error that names the file when it is missing.

    MATLAB looks only in ``./input``. The fallback to the example datasets of
    the repository lets the bundled examples run from any working folder.
    """
    path = find_input_file(name, input_dir)
    if path is None:
        raise FileNotFoundError(
            'The data file "%s" was not found in %s or among the example datasets of the toolbox.'
            % (name, os.path.abspath(input_dir))
        )
    return path


def load_matrix(path):
    """Load a header-free numeric text file as a two-dimensional array (MATLAB ``load``).

    A file that holds anything but numbers (a header line, a column of dates
    written as text) is an error that names the file and the text.
    """
    import warnings

    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # an empty file is reported by the caller
            data = np.loadtxt(path, ndmin=2)
    except ValueError as problem:
        raise SubEpiPredictError(
            "SubEpiPredict:DataFileNotNumeric",
            "The data file %s is not a header-free numeric matrix with the same number of values in every row "
            "(%s). Remove text such as a header line or a column of dates." % (os.path.basename(path), problem),
        ) from None
    return data


def data_file_name(cumulative1, cadtemporal, caddisease, datatype, cadregion, datestring):
    """Name of an input file: ``[cumulative-]<frequency>-<disease>-<datatype>-<region>-<mm-dd-yyyy>.txt``."""
    name = "%s-%s-%s-%s-%s.txt" % (cadtemporal, caddisease.strip(), datatype.strip(), cadregion.strip(), datestring)
    if cumulative1 == 1:
        name = "cumulative-" + name
    return name
