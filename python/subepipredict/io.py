"""Files of a run: data, output names, MAT files and CSV tables.

The output file names are those of the MATLAB toolbox, so results written by
the Python port are named like results written by MATLAB.

MAT files are written with ``scipy.io.savemat`` (MATLAB v5 format, readable
with MATLAB's ``load``). MATLAB saves its whole workspace in the fit files;
the port saves the variables that the fit-plot, ranking, forecast and Rt
functions read, under the same names (``RMSES``, ``PS``, ``Phatss``,
``curves``, ``data1``, ``P0``, ``AICc_best``, ``bootstrapDiagnostics`` ...).
CSV tables have the column names of MATLAB's ``writetable`` output and numbers
with 15 significant digits.
"""

from __future__ import annotations

import json
import os

import numpy as np
from scipy.io import loadmat, savemat
from scipy.io.matlab import mat_struct

from .errors import SubEpiPredictError
from .utils import (
    Struct,
    cadtemporal_of,
    data_file_name,
    datenum,
    load_matrix,
    num2str,
    resolve_input_file,
)

ALPHAQUANTILES = [0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500,
                  0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900, 0.950, 0.975, 0.990]  # fmt: skip

QUANT_NAMES_RANKED = ["Q_%.3f" % q for q in ALPHAQUANTILES]


# <============================================================================>
# <================================== Data ====================================>
# <============================================================================>


def load_area_series(filename, outbreak1, cumulative1, input_dir="./input"):
    """Incidence series of one area (column ``outbreak1``, from 1) of a data file.

    A cumulative file is differenced: ``[C(1); diff(C)]``. Returns
    ``(series, number of columns of the file)``.
    """
    data = load_matrix(resolve_input_file(filename, input_dir))
    if data.size == 0:
        raise SubEpiPredictError("SubEpiPredict:EmptyData", "The dataset is empty")
    ncols = data.shape[1]
    data1 = data[:, int(outbreak1) - 1].astype(float)
    if cumulative1 == 1:
        data1 = np.concatenate(([data1[0]], np.diff(data1)))  # Incidence curve
    return data1, ncols


def getData(cumulative1, cadtemporal, caddisease, datatype, cadregion, DT, datevecfirst1, datevecend1, date1,
            outbreak1, forecastingperiod, input_dir="./input"):  # fmt: skip
    """Observations of the forecast period from the most recent data file.

    Port of ``getData.m``. The file is the one dated ``datevecend1``. For
    daily and weekly data the result holds the ``forecastingperiod`` values
    that start at the date ``date1`` (``[yyyy, mm, dd]``, the first forecast
    date); for yearly data, the values after the year ``date1[0]``. With
    ``forecastingperiod = 0`` the result is NaN.
    """
    if forecastingperiod == 0:
        return np.nan
    filename1 = data_file_name(
        cumulative1, cadtemporal, caddisease, datatype, cadregion, _datestr(datenum(datevecend1))
    )
    data1, _ = load_area_series(filename1, outbreak1, cumulative1, input_dir)
    forecastingperiod = int(forecastingperiod)
    if DT == 365:
        nyears = int(date1[0]) - int(datevecfirst1[0]) + 1
        if data1.size < nyears + forecastingperiod:
            raise SubEpiPredictError(
                "SubEpiPredict:ForecastTooLong",
                "Forecasting period is too long and cannot be evaluated with the available data.",
            )
        return data1[nyears : nyears + forecastingperiod].copy()
    npoints = len(range(datenum(datevecfirst1), datenum(date1) + 1, int(DT)))
    first = npoints - 1  # zero-based index of the date date1
    if first + forecastingperiod > data1.size:
        raise SubEpiPredictError(
            "SubEpiPredict:ForecastTooLong",
            "Forecasting period is too long and cannot be evaluated with the available data.",
        )
    return data1[first : first + forecastingperiod].copy()


def _datestr(dn):
    from .utils import datestr

    return datestr(dn, "mm-dd-yyyy")


# <============================================================================>
# <============================= Output file names ============================>
# <============================================================================>


class OutputNames:
    """Builds the file names of one area and date (same patterns as MATLAB)."""

    def __init__(self, opts, outbreakx, caddate1, output_dir="./output"):
        self.output_dir = output_dir
        self.npatches_fixed = opts.npatches_fixed
        self.onset_fixed = opts.onset_fixed
        self.smoothfactor1 = opts.smoothfactor1
        self.flag1 = opts.flag1
        self.method1 = opts.method1
        self.dist1 = opts.dist1
        self.calibrationperiod1 = opts.calibrationperiod1
        self.cadtemporal = cadtemporal_of(opts.DT)
        self.caddisease = opts.caddisease
        self.datatype = opts.datatype
        self.cadregion = opts.cadregion
        self.outbreakx = outbreakx
        self.caddate1 = caddate1
        # adjusted data file of one area: <frequency>-<disease>-<datatype>-<region>-state-<area>-<date>
        self.cadfilename2 = "%s-%s-%s-%s-state-%s-%s" % (
            self.cadtemporal, self.caddisease, self.datatype, self.cadregion, num2str(outbreakx), caddate1,
        )  # fmt: skip

    def _path(self, name):
        return os.path.join(self.output_dir, name)

    def _fit_core(self):
        return "ensem-npatchesfixed-%s-onsetfixed-%s-smoothing-%s-%s-flag1-%s-method-%s-dist-%s-calibrationperiod-%s" % (
            num2str(self.npatches_fixed), num2str(self.onset_fixed), num2str(self.smoothfactor1), self.cadfilename2,
            num2str(self.flag1), num2str(self.method1), num2str(self.dist1), num2str(self.calibrationperiod1),
        )  # fmt: skip

    # --- MAT files -----------------------------------------------------------

    def abc_mat(self):
        """Screen: ``ABC-ensem-npatchesfixed-...-calibrationperiod-<c>.mat``."""
        return self._path("ABC-%s.mat" % self._fit_core())

    def ranked_mat(self, rank):
        """Ranked model: ``modifiedLogisticPatch-ensem-...-rank-<k>.mat``."""
        return self._path("modifiedLogisticPatch-%s-rank-%s.mat" % (self._fit_core(), num2str(rank)))

    def forecast_mat(self, forecastingperiod, rank):
        """Temporary forecast curves: ``Forecast-modifiedLogisticPatch-ensem-...-forecastingperiod-<h>-rank-<k>.mat``."""
        return self._path(
            "Forecast-modifiedLogisticPatch-%s-forecastingperiod-%s-rank-%s.mat"
            % (self._fit_core(), num2str(forecastingperiod), num2str(rank))
        )

    def aicc_bests_mat(self, topmodels):
        return self._path(
            "AICc_bests-modifiedLogisticPatch-ensem-npatchesfixed-%s-onsetfixed-%s-smoothing-%s-%s-%s-%s-state-%s"
            "-dateini-%s-flag1-%s-method-%s-dist-%s-calibrationperiod-%s-topmodels-%s.mat"
            % (num2str(self.npatches_fixed), num2str(self.onset_fixed), num2str(self.smoothfactor1), self.caddisease,
               self.datatype, self.cadregion, num2str(self.outbreakx), self.caddate1[:6] + self.caddate1[8:],
               num2str(self.flag1), num2str(self.method1), num2str(self.dist1), num2str(self.calibrationperiod1),
               num2str(topmodels))
        )  # fmt: skip

    # --- CSV files -----------------------------------------------------------

    def csv(self, prefix, horizon=None, weight=None, weight_label="weight_type"):
        """``<prefix>-onsetfixed-..-flag1-..-method-..-dist-..[-horizon-h][-weight_type-w]-<frequency>-<disease>-<datatype>-<region>-area-<area>-<date>.csv``."""
        name = "%s-onsetfixed-%s-flag1-%s-method-%s-dist-%s" % (
            prefix, num2str(self.onset_fixed), num2str(self.flag1), num2str(self.method1), num2str(self.dist1),
        )  # fmt: skip
        if horizon is not None:
            name += "-horizon-%s" % num2str(horizon)
        if weight is not None:
            name += "-%s-%s" % (weight_label, num2str(weight))
        name += "-%s-%s-%s-%s-area-%s-%s.csv" % (
            self.cadtemporal, self.caddisease, self.datatype, self.cadregion, num2str(self.outbreakx), self.caddate1,
        )  # fmt: skip
        return self._path(name)


# <============================================================================>
# <================================ MAT files =================================>
# <============================================================================>


def _matlab_value(v):
    """Make a Python value storable by savemat."""
    if v is None:
        return np.zeros((0, 0))
    if isinstance(v, dict):
        if len(v) == 0:
            return np.zeros((0, 0))
        return {str(k): _matlab_value(x) for k, x in v.items()}
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, str):
        return v
    if isinstance(v, (list, tuple)):
        if len(v) == 0:
            return np.zeros((0, 0))
        if all(isinstance(x, str) for x in v) or any(isinstance(x, (dict, list, tuple, np.ndarray)) for x in v):
            out = np.empty(len(v), dtype=object)
            for i, x in enumerate(v):
                out[i] = _matlab_value(x)
            return out
        try:
            return np.asarray(v, dtype=float)
        except (TypeError, ValueError):
            out = np.empty(len(v), dtype=object)
            for i, x in enumerate(v):
                out[i] = _matlab_value(x)
            return out
    if isinstance(v, np.ndarray) and v.dtype == bool:
        return v.astype(np.uint8)
    return v


def save_mat(path, variables):
    """Write a dict of variables to a MAT file (vectors as rows, like MATLAB)."""
    out = {}
    for key, value in variables.items():
        out[key] = _matlab_value(value)
    savemat(path, out, do_compression=True, oned_as="row", long_field_names=True)


def _from_matlab(v):
    if isinstance(v, mat_struct):
        return Struct({name: _from_matlab(getattr(v, name)) for name in v._fieldnames})
    if hasattr(v, "toarray"):  # sparse matrix saved by MATLAB
        return np.asarray(v.toarray(), dtype=float)
    if isinstance(v, np.ndarray):
        if v.dtype == object:
            if v.size == 1 and isinstance(v.ravel()[0], mat_struct):
                return _from_matlab(v.ravel()[0])
            return [_from_matlab(x) for x in v.ravel()]
        if v.dtype.kind in "US":
            if v.size == 0:
                return ""
            return str(v.ravel()[0]) if v.ndim <= 1 or v.shape[0] == 1 else [str(x) for x in v.ravel()]
    return v


def load_mat(path):
    """Read a MAT file into a :class:`Struct`.

    Numeric arrays come back two-dimensional, as MATLAB stores them; use
    :func:`scalar` and :func:`vector` to unwrap them. Structs become Structs,
    cell arrays lists, char arrays strings.
    """
    raw = loadmat(path, squeeze_me=False, struct_as_record=False, mat_dtype=True)
    return Struct({k: _from_matlab(v) for k, v in raw.items() if not k.startswith("__")})


def scalar(x):
    """The single number held by a 1-by-1 MAT array."""
    return float(np.asarray(x, dtype=float).ravel()[0])


def vector(x):
    """A MAT row or column as a one-dimensional float array."""
    return np.asarray(x, dtype=float).ravel()


def to_json(value):
    """Text of a (nested) settings value, for storage in a MAT file."""

    def convert(v):
        if isinstance(v, dict):
            return {str(k): convert(x) for k, x in v.items()}
        if isinstance(v, np.ndarray):
            return convert(v.tolist())
        if isinstance(v, (list, tuple)):
            return [convert(x) for x in v]
        if isinstance(v, (np.bool_,)):
            return bool(v)
        if isinstance(v, np.integer):
            return int(v)
        if isinstance(v, np.floating):
            v = float(v)
        if isinstance(v, float) and not np.isfinite(v):
            return "NaN" if np.isnan(v) else ("Inf" if v > 0 else "-Inf")
        return v

    return json.dumps(convert(value), sort_keys=True)


def from_json(text):
    """Inverse of :func:`to_json`."""

    def convert(v):
        if isinstance(v, dict):
            return Struct({k: convert(x) for k, x in v.items()})
        if isinstance(v, list):
            return [convert(x) for x in v]
        if v == "NaN":
            return float("nan")
        if v == "Inf":
            return float("inf")
        if v == "-Inf":
            return float("-inf")
        return v

    return convert(json.loads(text))


# <============================================================================>
# <================================ CSV tables ================================>
# <============================================================================>


def _fmt(x):
    x = float(x)
    if np.isnan(x):
        return "NaN"
    if np.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    return "%.15g" % x


def writetable(path, array, names):
    """Write a numeric table with a header row, like ``writetable(array2table(A), file)``."""
    array = np.asarray(array, dtype=float)
    if array.ndim == 1:
        array = array.reshape(1, -1)
    if array.size and array.shape[1] != len(names):
        raise ValueError("The table has %d columns but %d names were given." % (array.shape[1], len(names)))
    with open(path, "w", newline="") as fh:
        fh.write(",".join(names) + "\n")
        for row in array:
            fh.write(",".join(_fmt(v) for v in row) + "\n")


def readtable(path):
    """Read a CSV written by :func:`writetable`. Returns ``(names, array)``."""
    with open(path, "r") as fh:
        names = fh.readline().rstrip("\r\n").split(",")
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        array = np.genfromtxt(path, delimiter=",", skip_header=1, ndmin=2)
    if array.size == 0:
        array = np.zeros((0, len(names)))
    return names, array
