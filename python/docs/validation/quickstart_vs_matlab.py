"""Compare the quick-start fit of the Python port with the MATLAB toolbox.

The MATLAB numbers are those shown in the reference figures of the repository
(docs/images/ranking_profiles.png and docs/images/parameters.png), which were
made with fit revision m1 and random seed 1 (docs/images/README_figures_info.txt):
the onset thresholds of the four ranked models, the AICc as printed with six
significant digits, and the estimates of ranked model 1 with three.

Run from the ``python`` folder (the fit takes several minutes)::

    python docs/validation/quickstart_vs_matlab.py                   # fit in ./output, then compare
    python docs/validation/quickstart_vs_matlab.py --output ./output  # compare a fit that is already there
    python docs/validation/quickstart_vs_matlab.py --parallel

The comparison is printed and written to
docs/validation/quickstart_vs_matlab.txt.
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))
sys.dont_write_bytecode = True

import matplotlib  # noqa: E402

matplotlib.use("Agg")

import numpy as np  # noqa: E402

import subepipredict as sp  # noqa: E402
from subepipredict.figures import _localNumber  # noqa: E402

# read off the MATLAB figures
MATLAB_THRESHOLDS = ["1140.5536", "0", "27361.7865", "28501.8401"]
MATLAB_AICC = "-1.09196e+06"
MATLAB_RANK1 = {"r_1": "0.561", "r_2": "4.63", "p_1": "0.826", "p_2": "0.659", "a_1": "1", "a_2": "1",
                "K_1": "80400", "K_2": "24800"}  # fmt: skip


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="./output")
    parser.add_argument("--parallel", action="store_true")
    args = parser.parse_args()

    area, origin = 52, "05-11-2020"
    opts = sp.options()
    names = sp.OutputNames(opts, area, origin, args.output)
    if not os.path.isfile(names.ranked_mat(opts.topmodelsx)):
        sp.Run_Fit_subepidemicFramework(area, origin, options=opts, fit_settings=dict(useParallel=args.parallel),
                                        output_dir=args.output, plot=False)  # fmt: skip
        sp.shutdown_pool()

    abc = sp.load_mat(names.abc_mat())
    RMSES = np.asarray(abc.RMSES, dtype=float)
    PS = np.asarray(abc.PS, dtype=float)

    from subepipredict.utils import num2str

    lines = []
    lines.append("Quick start (U.S. daily COVID-19 deaths, column 52, data to 05-11-2020; method1 = 3, flag1 = 1,")
    lines.append("up to 2 sub-epidemics, 20 start points, random seed 1): MATLAB reference figures against the Python port.")
    lines.append("")
    lines.append("Ranked models")
    lines.append("  rank  sub-epidemics   C_thr MATLAB   C_thr Python     AICc MATLAB    AICc Python   AICc Python (full)")
    ok = True
    for k in range(4):
        thr = num2str(RMSES[k, 1])
        aicc = "%.6g" % RMSES[k, 2]
        same = thr == MATLAB_THRESHOLDS[k] and aicc == MATLAB_AICC and int(RMSES[k, 0]) == 2
        ok = ok and same
        lines.append("  %4d  %13d   %12s   %12s   %13s  %13s   %.4f   %s" % (
            k + 1, RMSES[k, 0], MATLAB_THRESHOLDS[k], thr, MATLAB_AICC, aicc, RMSES[k, 2], "same" if same else "DIFFERENT"))
    lines.append("")
    lines.append("Estimates of ranked model 1 (three significant digits, as in the titles of the parameter figure)")
    lines.append("  parameter      MATLAB      Python   Python (full)")
    P = PS[0]
    values = {"r_1": P[0], "r_2": P[1], "p_1": P[2], "p_2": P[3], "a_1": P[4], "a_2": P[5], "K_1": P[6], "K_2": P[7]}
    for name, matlab in MATLAB_RANK1.items():
        python = _localNumber(values[name])
        relative = abs(values[name] - float(matlab)) / float(matlab)
        same = python == matlab
        note = "same" if same else ("differs by %.2g%% (last digit shown)" % (100 * relative))
        ok = ok and relative < 2e-3
        lines.append("  %-9s %11s %11s   %-14.8g %s" % (name, matlab, python, values[name], note))
    lines.append("")
    lines.append("Result: %s" % ("the ranking and the estimates agree with the MATLAB figures."
                                 if ok else "differences beyond the last digit shown; see above."))
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(HERE, "quickstart_vs_matlab.txt"), "w") as fh:
        fh.write(text + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
