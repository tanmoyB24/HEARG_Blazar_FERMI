#!/usr/bin/env python
"""
Generate FAKE results in the same format run_fermipy.py writes, so you can test
analyze_population.py and preview every figure/table layout before real fits finish.

All summaries are flagged "synthetic": true, which makes analyze_population.py stamp
"SYNTHETIC DEMO" on every figure. The numbers are random draws with loosely
class-dependent trends - they are NOT results and must never appear in a talk.

    python tools/make_synthetic_results.py --out results_demo
    python analyze_population.py --results results_demo --out demo_out
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import blazarlib as bl  # noqa: E402


def make_one(rng, name, cls, out: Path, nbins=120):
    fsrq = cls == "FSRQ"
    gamma = rng.normal(2.45 if fsrq else 2.0, 0.12)
    beta = abs(rng.normal(0.10 if fsrq else 0.04, 0.03))
    flux0 = 10 ** rng.uniform(-8.3, -6.6)  # ph/cm2/s
    # SED (E^2 dN/dE, MeV cm^-2 s^-1) from a log-parabola on 8 log bins
    edges = np.logspace(2, np.log10(3e5), 9)
    e_ref = np.sqrt(edges[:-1] * edges[1:])
    eb = 1000.0
    x = e_ref / eb
    shape = x ** (2 - gamma - beta * np.log(x))
    e2 = flux0 * 3e3 * shape / (shape[3]) * 1.0e0 * 1e-1
    ts = np.clip(rng.normal(1, 0.15, 8) * (e2 / e2.max()) * 400 * (flux0 / 1e-7), 0, None)
    err = e2 * np.clip(2.0 / np.sqrt(np.maximum(ts, 1)), 0.05, 1.5)
    np.savez(out / "sed.npz", e_min=edges[:-1], e_max=edges[1:], e_ref=e_ref, e2dnde=e2 * rng.normal(1, 0.05, 8),
             e2dnde_err=err, e2dnde_ul95=e2 + 2 * err, ts=ts)
    # light curve: red-noise-like log-normal variability
    sig = rng.uniform(0.6, 1.0) if fsrq else rng.uniform(0.25, 0.6)
    walk = np.cumsum(rng.normal(0, 0.35, nbins))
    walk = (walk - walk.mean()) / walk.std() * sig
    flux = flux0 * np.exp(walk - sig ** 2 / 2)
    ferr = flux * np.clip(rng.normal(0.12, 0.03, nbins), 0.05, None) * (flux0 / flux) ** 0.2
    tmin = bl.iso_to_mjd("2008-08-04") + np.arange(nbins) * 30.0
    lc_ts = np.clip((flux / ferr) ** 2 * rng.normal(1, 0.1, nbins), 0, None)
    np.savez(out / "lightcurve.npz", tmin_mjd=tmin, tmax_mjd=tmin + 30, flux=flux, flux_err=ferr,
             flux_ul95=flux + 2 * ferr, ts=lc_ts)
    lc = dict(np.load(out / "lightcurve.npz"))
    summ = dict(
        name=name, fgl_name="synthetic", cls=cls, ra=0.0, dec=0.0, cat_class=cls.lower(),
        cat_spectrum_type="LogParabola", ts=float(ts.sum()), npred=1e4, flux=float(flux0),
        flux_err=float(flux0 * 0.03), index=float(gamma), index_err=0.03, beta=float(beta), beta_err=0.02,
        ts_curv=float(max(0, rng.normal(60 if fsrq else 10, 25))), model="LogParabola",
        synthetic=True, status={"synthetic": "yes"},
    )
    summ.update(bl.lightcurve_stats(lc))
    bl.dump_json(summ, out / "summary.json")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="results_demo")
    ap.add_argument("--sources", default="sources.csv")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    srcs = bl.read_sources(args.sources)
    for _, r in srcs.iterrows():
        d = Path(args.out) / bl.slugify(r["name"])
        d.mkdir(parents=True, exist_ok=True)
        make_one(rng, r["name"], r["class"], d)
    print(f"wrote {len(srcs)} synthetic sources to {args.out}")


if __name__ == "__main__":
    main()
