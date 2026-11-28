"""Compare legacy (count-based) and kernel-dependent r bounds on bundled series.

usage: python r_bounds_by_kernel.py <python_pkg_dir> <input_dir> <out_csv> [single|screen|richards3]

single     one sub-epidemic, every window (about 4 CPU minutes)
screen     full screen with up to two sub-epidemics, 8 windows, logistic and Gompertz (about 40 CPU minutes)
richards3  the three Richards windows discussed in r_bounds_by_kernel.txt; NSTART sets the start points

Every fit: Poisson likelihood (method1 = 1), 10 start points, seed 1. "legacy" is rBoundsByKernel = False.
Summary: python r_bounds_by_kernel_summary.py <out_csv> [-v]
"""
import sys, os, time, warnings, itertools, csv
import numpy as np
from multiprocessing import Pool

PKG, INP, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
MODE = sys.argv[4] if len(sys.argv) > 4 else "single"   # single | screen
QUICK = False
NP = 1 if MODE in ("single", "richards3") else 2
NSTART = int(os.environ.get("NSTART", 10))
sys.path.insert(0, PKG)
import matplotlib; matplotlib.use("Agg")
import subepipredict as sp
from subepipredict.fitting import fittingModifiedLogisticFunctionPatchABC
from subepipredict.model import simulateSubepidemic, unpackSubepidemicParameters
from subepipredict.settings import subepidemicFitSettings
from subepipredict.search import subepidemicSearchSpace


def incidence(fname, col, cumulative):
    d = np.loadtxt(os.path.join(INP, fname), ndmin=2)[:, col - 1].astype(float)
    return np.concatenate(([d[0]], np.diff(d))) if cumulative else d


def jobs():
    J = []
    deaths = "cumulative-daily-coronavirus-deaths-USA-05-09-2022.txt"
    D = np.loadtxt(os.path.join(INP, deaths), ndmin=2)
    cols = [52] + [int(c) + 1 for c in np.argsort(D[-1, :51])[::-1][:3]]  # USA + 3 largest states
    ends = range(150, 790, 80)
    if QUICK:
        cols, ends = cols[:1], [230]
    for c in cols:
        y = np.maximum(incidence(deaths, c, True), 0)
        for e in ends:
            J.append(("covid-deaths col%d" % c, y, e, 90, 14, 7))
    if not QUICK:
        y = np.maximum(incidence("cumulative-weekly-monkeypox-cases-USA-07-20-2023.txt", 1, True), 0)
        for e in (12, 16, 20, 28, 36):
            J.append(("mpox weekly USA", y, e, 20, 4, 1))
        for f, lab in (("daily-Ebola-cases-Uganda-11-12-2022.txt", "ebola Uganda"),
                       ("daily-SARS-cases-Singapore-02-25-2003.txt", "SARS Singapore")):
            y = incidence(f, 1, False)
            for e in (30, 45, 60):
                J.append((lab, y, e, 60, 10, 1))
    y = incidence("yearly-HIV-cases-World_Regions-01-01-2019.txt", 1, False)
    J.append(("HIV yearly region1", y, y.size - 3, 25, 3, 1))
    kernels = (3, 4, 5)
    if MODE == "screen":   # full 2-sub-epidemic screen: a subset of windows, logistic and Gompertz
        keep = {("covid-deaths col52", 150), ("covid-deaths col52", 390), ("covid-deaths col52", 630),
                ("covid-deaths col%d" % cols[1], 310), ("covid-deaths col%d" % cols[2], 470),
                ("mpox weekly USA", 20), ("ebola Uganda", 60), ("SARS Singapore", 45)}
        J = [j for j in J if (j[0], j[2]) in keep]
        kernels = (3, 5)
    if MODE == "richards3":
        J = [j for j in J if (j[0], j[2]) in {("covid-deaths col52", 150), ("covid-deaths col52", 310), ("covid-deaths col52", 390)}]
        kernels = (4,)
    out = []
    for (lab, y, e, cal, h, sm), flag1, mode in itertools.product(J, kernels, ("legacy", "new")):
        out.append((lab, y, e, cal, h, sm, flag1, mode))
    return out


def run(job):
    lab, y, end, cal, h, sm, flag1, mode = job
    warnings.simplefilter("ignore")
    S = subepidemicFitSettings({"rBoundsByKernel": mode == "new", "rngSeed": 1})
    win = y[max(0, end - cal):end].copy()
    nz = np.flatnonzero(win > 0)
    win = win[nz[0]:]
    t0 = time.time()
    data1 = np.column_stack([np.arange(win.size, dtype=float), win])
    rng = np.random.default_rng(1)
    try:
        RMSES, PS, n, thr, P, screen = fittingModifiedLogisticFunctionPatchABC(
            "x", data1, 1, np.arange(1, win.size + 1), 0, flag1, NSTART, method1=1, dist1=1,
            npatches_fixed=NP, onset_fixed=0, smoothfactor1=sm, calibrationperiod1=cal,
            fit_settings=S, save=False, plot=False, verbose=False, rng=rng)
    except Exception as ex:  # record failures
        return dict(series=lab, end=end, flag1=flag1, mode=mode, error=str(ex)[:120])
    n = int(n)
    LB, UB, _, _ = subepidemicSearchSpace(n, screen.data if "data" in screen else win, flag1, 1, S)
    rs, ps, as_, Ks = unpackSubepidemicParameters(P, n)[:4]
    tv = np.asarray(screen.timevect, float)
    tf = np.concatenate([tv, tv[-1] + np.arange(1, h + 1)])
    _, x, totinc, _, _ = simulateSubepidemic(tf, screen.I0, n, 0, thr, flag1, rs, ps, as_, Ks)
    fit = np.asarray(totinc[: tv.size], float)
    fc = np.asarray(totinc[tv.size:], float)
    obs = y[end:end + h]
    k = min(fc.size, obs.size)
    ydat = np.asarray(screen.ydata, float)
    return dict(series=lab, end=end, flag1=flag1, mode=mode, first_count=win[0], first2_mean=win[:2].mean(),
                n=n, thr=float(thr), aicc=float(RMSES[0, 2]), ncand=RMSES.shape[0],
                r=";".join("%.4g" % v for v in rs), K=";".join("%.4g" % v for v in Ks),
                a=";".join("%.3g" % v for v in as_),
                rLB=float(LB[0]), rUB=float(UB[0]),
                r_at_LB=int(np.any(np.isclose(rs, LB[:n], rtol=1e-3))),
                r_at_UB=int(np.any(np.isclose(rs, UB[:n], rtol=1e-3))),
                fit_mae=float(np.mean(np.abs(fit - ydat))) if fit.size == ydat.size else np.nan,
                fc_mae=float(np.mean(np.abs(fc[:k] - obs[:k]))) if k else np.nan,
                cum_fit_vs_data=float(Ks.sum() / max(win.sum(), 1)), secs=round(time.time() - t0, 1))


if __name__ == "__main__":
    J = jobs()
    print(len(J), "fits", flush=True)
    keys = None
    with Pool(int(os.environ.get("NPROC", os.cpu_count()))) as pool, open(OUT, "w", newline="") as fh:
        for i, row in enumerate(pool.imap_unordered(run, J)):
            if keys is None:
                keys = ["series", "end", "flag1", "mode", "first_count", "first2_mean", "n", "thr", "aicc", "ncand",
                        "r", "a", "K", "rLB", "rUB", "r_at_LB", "r_at_UB", "fit_mae", "fc_mae", "cum_fit_vs_data",
                        "secs", "error"]
                w = csv.DictWriter(fh, keys, extrasaction="ignore"); w.writeheader()
            w.writerow(row); fh.flush()
            print(i + 1, row.get("series"), row.get("end"), row.get("flag1"), row.get("mode"),
                  row.get("aicc", row.get("error")), row.get("secs"), flush=True)
