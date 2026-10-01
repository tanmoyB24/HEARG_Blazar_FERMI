#!/usr/bin/env python
"""
Per-source Fermi-LAT analysis with fermipy (generalises HEARG Tutorial 3).

For every source in sources.csv:
  setup -> optimize -> free neighbours -> fit (NEWMINUIT) -> residmap/tsmap
  -> prune weak sources -> PowerLaw vs LogParabola test (TS_curv)
  -> localize -> SED -> light curve -> write ROI

Outputs go to results/<slug>/ :
  summary.json, sed.npz, lightcurve.npz, config.yaml, fermipy ROI files

Usage
-----
    python run_fermipy.py --dry-run                 # only write configs, no fermipy needed
    python run_fermipy.py --only "3C 279"           # one source
    python run_fermipy.py                           # all sources (resumes; skips finished)
    python run_fermipy.py --force --skip-lc         # redo everything except light curves

NOTE: the fermipy calls follow the fermipy documentation but were written without
access to the Fermi Science Tools, so test on ONE source first. Every optional step
is wrapped so a failure is recorded in summary.json["status"] instead of aborting.
"""
from __future__ import annotations

import argparse
import logging
import time
import traceback
from pathlib import Path

import numpy as np
import yaml

import blazarlib as bl

log = logging.getLogger("blazar")


def fnum(x, default=np.nan):
    try:
        return float(x)
    except Exception:
        return default


# ----------------------------------------------------------------------------
# Spectral helpers
# ----------------------------------------------------------------------------
def get_par(gta, name, par):
    """(value, error) of a spectral parameter in physical (scaled) units."""
    p = gta.roi[name].spectral_pars.get(par)
    if p is None:
        return np.nan, np.nan
    sc = fnum(p.get("scale", 1.0), 1.0)
    return fnum(p.get("value")) * sc, abs(fnum(p.get("error", np.nan)) * sc)


def photon_index(gta, name):
    """Photon index Gamma (dN/dE ~ E^-Gamma) for LogParabola (alpha) or PowerLaw (Index)."""
    stype = gta.roi[name]["SpectrumType"]
    if stype == "LogParabola":
        return get_par(gta, name, "alpha")
    if stype == "PowerLaw":
        v, e = get_par(gta, name, "Index")
        return -v, e
    return np.nan, np.nan


def ensure_logparabola(gta, name):
    """Convert a PowerLaw source to a LogParabola with the same normalisation/index."""
    src = gta.roi[name]
    stype = src["SpectrumType"]
    if stype == "LogParabola":
        return True
    if stype != "PowerLaw":
        return False
    pp = src.spectral_pars
    pref, ind, scl = pp["Prefactor"], pp["Index"], pp["Scale"]
    gamma = -fnum(ind["value"]) * fnum(ind["scale"])
    new = {
        "norm": {"value": pref["value"], "scale": pref["scale"], "min": pref["min"], "max": pref["max"], "free": True},
        "alpha": {"value": gamma, "scale": 1.0, "min": 0.0, "max": 5.0, "free": True},
        "beta": {"value": 0.0, "scale": 1.0, "min": -1.0, "max": 5.0, "free": True},
        "Eb": {"value": scl["value"], "scale": scl["scale"], "min": scl["min"], "max": scl["max"], "free": False},
    }
    gta.set_source_spectrum(name, "LogParabola", spectrum_pars=new)
    return True


def curvature_test(gta, name, threshold, optimizer):
    """Nested-model test: LogParabola (beta free) vs PowerLaw (beta = 0).

    TS_curv = 2 * (logL_LP - logL_PL). The source is left in the preferred model.
    """
    out = dict(ts_curv=np.nan, ll_pl=np.nan, ll_lp=np.nan, model="unknown")
    if not ensure_logparabola(gta, name):
        out["model"] = f"skipped ({gta.roi[name]['SpectrumType']})"
        return out
    gta.free_parameter(name, "norm", True)
    gta.free_parameter(name, "alpha", True)

    gta.set_parameter(name, "beta", 0.0)
    gta.free_parameter(name, "beta", False)
    out["ll_pl"] = fnum(gta.fit(optimizer=optimizer)["loglike"])

    gta.set_parameter(name, "beta", 0.05)
    gta.free_parameter(name, "beta", True)
    out["ll_lp"] = fnum(gta.fit(optimizer=optimizer)["loglike"])

    d = out["ll_lp"] - out["ll_pl"]
    if d < -1.0:
        log.warning("%s: logL(LP) < logL(PL) by %.2f - check fit convergence / loglike sign", name, -d)
    out["ts_curv"] = max(0.0, 2.0 * d)

    if out["ts_curv"] >= threshold:
        out["model"] = "LogParabola"
    else:
        gta.set_parameter(name, "beta", 0.0)
        gta.free_parameter(name, "beta", False)
        gta.fit(optimizer=optimizer)
        out["model"] = "PowerLaw"
    return out


def arrays_only(d, keys=None):
    """Keep numeric 1-D arrays from a fermipy result dict so they can go in an npz."""
    out = {}
    for k, v in d.items():
        if keys and k not in keys:
            continue
        try:
            a = np.asarray(v, dtype=float)
        except Exception:
            continue
        if a.ndim == 1 and a.size > 1:
            out[k] = a
    return out


def summarize(gta, name, row):
    src = gta.roi[name]
    gamma, gamma_err = photon_index(gta, name)
    beta, beta_err = get_par(gta, name, "beta")
    return dict(
        name=row["name"], fgl_name=row["fgl_name"], cls=row["class"],
        ra=row["ra"], dec=row["dec"], cat_class=row["cat_class"],
        cat_spectrum_type=row["spectrum_type"],
        ts=fnum(src["ts"]), npred=fnum(src["npred"]),
        flux=fnum(src["flux"]), flux_err=fnum(src["flux_err"]),
        eflux=fnum(src["eflux"]), eflux_err=fnum(src["eflux_err"]),
        index=gamma, index_err=gamma_err, beta=beta, beta_err=beta_err,
        synthetic=False,
    )


# ----------------------------------------------------------------------------
# One source
# ----------------------------------------------------------------------------
def process_source(row, s, args):
    name = row["fgl_name"]
    slug = bl.slugify(row["name"])
    outdir = Path(s["_root"]) / s["paths"]["results_dir"] / slug
    outdir.mkdir(parents=True, exist_ok=True)
    if (outdir / "summary.json").exists() and not args.force:
        log.info("%s: already done (use --force to redo)", row["name"])
        return

    evfile, scfile = bl.find_data(slug, s)
    cfg = bl.build_fermipy_config(row, s, outdir, evfile, scfile)
    cfg_path = outdir / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    if args.dry_run:
        log.info("%s: wrote %s (dry run)", row["name"], cfg_path)
        return

    from fermipy.gtanalysis import GTAnalysis  # imported late so --dry-run works anywhere

    a, lcs = s["analysis"], s["lightcurve"]
    opt = a["optimizer"]
    status, t0 = {}, time.time()

    gta = GTAnalysis(str(cfg_path), logging={"verbosity": 3})
    gta.setup()
    gta.optimize()

    # --- fit: neighbours (norm), diffuse (norm), target (all shape parameters) ---
    gta.free_sources(distance=a["free_radius"], pars="norm")
    gta.free_source("galdiff", pars="norm")
    gta.free_source("isodiff", pars="norm")
    gta.free_source(name)
    fit = gta.fit(optimizer=opt)
    status["fit_success"] = bool(fit.get("fit_success", True))

    # --- model diagnostics ---
    for step, call in (("residmap", lambda: gta.residmap(make_plots=True)),
                       ("tsmap", lambda: gta.tsmap(make_plots=True))):
        try:
            call()
            status[step] = "ok"
        except Exception as exc:
            status[step] = f"failed: {exc}"

    if a.get("find_sources"):
        try:
            gta.find_sources(sqrt_ts_threshold=4.0, min_separation=0.5)
            status["find_sources"] = "ok"
        except Exception as exc:
            status["find_sources"] = f"failed: {exc}"

    try:
        gta.delete_sources(minmax_ts=[-np.inf, a["delete_below_ts"]], exclude=[name])
        gta.fit(optimizer=opt)
        status["prune"] = "ok"
    except Exception as exc:
        status["prune"] = f"failed: {exc}"

    # --- curvature test (leaves target in the preferred model) ---
    curv = dict(ts_curv=np.nan, ll_pl=np.nan, ll_lp=np.nan, model="unknown")
    try:
        curv = curvature_test(gta, name, a["curvature_ts_threshold"], opt)
        status["curvature"] = "ok"
    except Exception as exc:
        status["curvature"] = f"failed: {exc}"
        log.error(traceback.format_exc())

    summary = summarize(gta, name, row)
    summary.update(curv)
    summary["time"] = dict(start=s["time"]["start"], end=s["time"]["end"])

    # --- localisation (does not move the source) ---
    try:
        loc = gta.localize(name, update=False, make_plots=False)
        summary.update(loc_ra=fnum(loc.get("ra")), loc_dec=fnum(loc.get("dec")),
                       loc_offset=fnum(loc.get("pos_offset")), loc_r95=fnum(loc.get("pos_r95")))
        status["localize"] = "ok"
    except Exception as exc:
        status["localize"] = f"failed: {exc}"

    # --- SED ---
    if not args.skip_sed:
        try:
            sed = gta.sed(name, make_plots=False, outfile=f"{slug}_sed.fits")
            np.savez(outdir / "sed.npz", **arrays_only(sed))
            status["sed"] = "ok"
        except Exception as exc:
            status["sed"] = f"failed: {exc}"
            log.error(traceback.format_exc())

    # --- light curve ---
    if not args.skip_lc:
        try:
            lc = gta.lightcurve(name, binsz=lcs["bin_days"] * 86400.0, free_radius=lcs["free_radius"],
                                use_scaled_srcmap=True, multithread=True, nthread=lcs["nthread"],
                                save_bin_data=False, make_plots=False,
                                outfile=f"{slug}_lightcurve.fits")
            np.savez(outdir / "lightcurve.npz", **arrays_only(lc))
            summary.update(bl.lightcurve_stats(lc, lcs["ts_min"], lcs["flare_factor"]))
            status["lightcurve"] = "ok"
        except Exception as exc:
            status["lightcurve"] = f"failed: {exc}"
            log.error(traceback.format_exc())

    try:
        gta.write_roi(f"{slug}_fit", make_plots=False)
    except Exception as exc:
        status["write_roi"] = f"failed: {exc}"

    summary["status"] = status
    summary["runtime_min"] = (time.time() - t0) / 60.0
    bl.dump_json(summary, outdir / "summary.json")
    log.info("%s done in %.1f min  TS=%.0f  Gamma=%.2f  TS_curv=%.1f", row["name"],
             summary["runtime_min"], summary["ts"], summary["index"], summary["ts_curv"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--settings", default="settings.yaml")
    ap.add_argument("--only", help="process just this source name (as in sources.csv)")
    ap.add_argument("--force", action="store_true", help="redo sources that already have a summary.json")
    ap.add_argument("--dry-run", action="store_true", help="write configs only")
    ap.add_argument("--skip-lc", action="store_true")
    ap.add_argument("--skip-sed", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    s = bl.load_settings(args.settings)
    cat = bl.Catalog(s["paths"]["catalog"])
    srcs = bl.read_sources(f'{s["_root"]}/{s["paths"]["sources"]}')

    for _, r in srcs.iterrows():
        if args.only and r["name"] != args.only:
            continue
        hit = cat.lookup(r["fgl_name"], r.get("assoc"))
        if hit is None:
            log.error("%s: %s not in catalog - run verify_sources.py", r["name"], r["fgl_name"])
            continue
        row = {**r.to_dict(), **hit}  # catalog values (name, ra, dec, ...) take precedence
        row["name"], row["class"] = r["name"], r["class"]
        row["fgl_name"] = hit["fgl_name"]
        try:
            process_source(row, s, args)
        except Exception:
            log.error("%s FAILED:\n%s", r["name"], traceback.format_exc())


if __name__ == "__main__":
    main()
