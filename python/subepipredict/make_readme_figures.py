"""Make the six reference figures of the README from the quick-start fit.

Run from the ``python`` folder::

    python -m subepipredict.make_readme_figures                  # fit, then plot
    python -m subepipredict.make_readme_figures --output ./output  # plot a fit that is already there
    python -m subepipredict.make_readme_figures --parallel       # local solves on a process pool

The settings are those of ``docs/images/README_figures_info.txt`` of the
MATLAB toolbox: U.S. daily COVID-19 deaths, column 52, data to 05-11-2020,
quick-start options, random seed 1, forecast horizon 30, weight code 1.
The figures are written to ``docs/images`` of the ``python`` folder.
"""

from __future__ import annotations

import argparse
import os
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--output", default=None, help="folder of the fit files (default: ./output; fitted if empty)")
    parser.add_argument("--images", default=None, help="folder for the figures (default: docs/images next to the package)")
    parser.add_argument("--parallel", action="store_true", help="run the local solves on a process pool")
    parser.add_argument("--dpi", type=int, default=100)
    args = parser.parse_args(argv)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import subepipredict as sp

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    images = args.images or os.path.join(here, "docs", "images")
    output = args.output or "./output"
    os.makedirs(images, exist_ok=True)

    area, origin = 52, "05-11-2020"
    opts = sp.options()
    names = sp.OutputNames(opts, area, origin, output)
    settings = dict(useParallel=args.parallel)
    if not os.path.isfile(names.ranked_mat(opts.topmodelsx)):
        sp.Run_Fit_subepidemicFramework(area, origin, options=opts, fit_settings=settings, output_dir=output, plot=False)
        sp.shutdown_pool()

    def save(number, name):
        path = os.path.join(images, name)
        plt.figure(number).savefig(path, dpi=args.dpi, facecolor="white")
        print("wrote", path)

    ranked = sp.plotFit_subepidemicFramework(area, origin, options=opts, output_dir=output, verbose=False)
    save(301 + ranked.shape[0], "model_fit.png")
    save(301, "parameters.png")

    sp.plotRankings_subepidemicFramework(area, origin, options=opts, output_dir=output, verbose=False)
    save(99, "rankings.png")
    save(100, "ranking_profiles.png")

    sp.plotForecast_subepidemicFramework(area, origin, 30, 1, options=opts, output_dir=output, verbose=False)
    save(10, "ranked_forecasts.png")
    save(500, "ensembles.png")
    plt.close("all")
    return 0


if __name__ == "__main__":
    sys.exit(main())
