"""
blazarlib.py - shared helpers for the Fermi-LAT blazar project.

Everything that does NOT need fermipy lives here so it can be tested and reused
by the pipeline, the population analysis and the notebooks.
"""
from __future__ import annotations

import json
import math
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import yaml

# ----------------------------------------------------------------------------
# Time conversion
# ----------------------------------------------------------------------------
# Fermi mission elapsed time (MET): seconds since 2001-01-01 00:00:00 UTC.
# MJDREF below is the value written in Fermi FITS headers (TT system).
# Treating a UTC date as TT introduces a ~1 minute error, irrelevant for
# multi-year analysis windows.
MJDREF = 51910.0 + 7.428703703703703e-4
_MJD_UNIX_EPOCH = 40587.0


def iso_to_mjd(s: str) -> float:
    dt = datetime.fromisoformat(str(s).strip())
    return _MJD_UNIX_EPOCH + (dt - datetime(1970, 1, 1)).total_seconds() / 86400.0


def mjd_to_met(mjd: float) -> float:
    return (mjd - MJDREF) * 86400.0


def met_to_mjd(met: float) -> float:
    return met / 86400.0 + MJDREF


def iso_to_met(s: str) -> float:
    return mjd_to_met(iso_to_mjd(s))


# ----------------------------------------------------------------------------
# Settings / names
# ----------------------------------------------------------------------------
def load_settings(path: str | Path = "settings.yaml") -> dict:
    path = Path(path)
    with open(path) as fh:
        s = yaml.safe_load(fh)
    s["_root"] = str(path.resolve().parent)
    return s


def slugify(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_").lower()


def normalize_4fgl(name: str) -> str | None:
    """'4FGLJ2253.9+1609', 'J2253.9+1609' or '4FGL J2253.9+1609' -> '4FGL J2253.9+1609'."""
    m = re.search(r"J\d{4}\.\d[+-]\d{4}[a-z]?", str(name))
    return f"4FGL {m.group(0)}" if m else None


def _norm_assoc(s: str) -> str:
    return re.sub(r"[^A-Z0-9+\-.]", "", str(s).upper())


def read_sources(path: str | Path):
    import pandas as pd

    df = pd.read_csv(path, comment="#", skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    for col in ("name", "fgl_name", "class"):
        if col not in df.columns:
            raise ValueError(f"{path} must have columns: name, fgl_name, class")
    df["class"] = df["class"].str.upper().str.strip()
    return df


# ----------------------------------------------------------------------------
# 4FGL catalog access
# ----------------------------------------------------------------------------
class Catalog:
    """Minimal reader for the Fermi-LAT 4FGL point-source FITS catalog."""

    def __init__(self, path: str | Path):
        from astropy.io import fits

        with fits.open(path) as hdul:
            hdu = hdul["LAT_Point_Source_Catalog"] if "LAT_Point_Source_Catalog" in hdul else hdul[1]
            d = hdu.data
            cols = set(d.columns.names)

            def strcol(c):
                return np.char.strip(np.asarray(d[c]).astype(str)) if c in cols else np.array([""] * len(d))

            def numcol(c):
                return np.asarray(d[c], dtype=float) if c in cols else np.full(len(d), np.nan)

            self.name = strcol("Source_Name")
            self.ra = numcol("RAJ2000")
            self.dec = numcol("DEJ2000")
            self.cls = strcol("CLASS1")
            self.assoc = strcol("ASSOC1")
            self.sptype = strcol("SpectrumType")
            self.signif = numcol("Signif_Avg")
            self.r95 = numcol("Conf_95_SemiMajor")
        self.path = str(path)

    def __len__(self):
        return len(self.name)

    def _row(self, i: int) -> dict:
        return dict(
            fgl_name=str(self.name[i]), ra=float(self.ra[i]), dec=float(self.dec[i]),
            cat_class=str(self.cls[i]), assoc=str(self.assoc[i]),
            spectrum_type=str(self.sptype[i]), signif=float(self.signif[i]),
            r95=float(self.r95[i]),
        )

    def lookup(self, fgl_name: str | None = None, assoc: str | None = None) -> dict | None:
        """Find a source by 4FGL name; fall back to its association (e.g. '3C 454.3')."""
        norm = normalize_4fgl(fgl_name) if fgl_name else None
        if norm:
            idx = np.where(self.name == norm)[0]
            if len(idx):
                return self._row(int(idx[0]))
        if assoc:
            target = _norm_assoc(assoc)
            for i, a in enumerate(self.assoc):
                if a and _norm_assoc(a) == target:
                    return self._row(i)
        return None

    def nearest(self, ra: float, dec: float):
        """Return (row_dict, separation_deg) of the catalog source closest to (ra, dec)."""
        r1, d1 = np.radians(ra), np.radians(dec)
        r2, d2 = np.radians(self.ra), np.radians(self.dec)
        a = np.sin((d2 - d1) / 2) ** 2 + np.cos(d1) * np.cos(d2) * np.sin((r2 - r1) / 2) ** 2
        sep = np.degrees(2 * np.arcsin(np.sqrt(np.clip(a, 0, 1))))
        i = int(np.nanargmin(sep))
        return self._row(i), float(sep[i])


def class_group(code: str) -> str:
    """Collapse 4FGL CLASS1 codes into coarse groups (case-insensitive)."""
    c = str(code).strip().lower()
    if c == "fsrq":
        return "FSRQ"
    if c == "bll":
        return "BL Lac"
    if c == "bcu":
        return "BCU"
    if c in {"agn", "rdg", "css", "nlsy1", "sey", "sbg", "ssrq"}:
        return "Other AGN"
    if c in {"psr", "pwn", "snr", "spp", "hmb", "lmb", "bin", "glc", "sfr", "nov", "gal"}:
        return "Galactic"
    if c in {"", "unk"}:
        return "Unassociated"
    return "Other"


# ----------------------------------------------------------------------------
# Variability statistics
# ----------------------------------------------------------------------------
def fractional_variability(flux, err):
    """Excess-variance fractional variability F_var and its error (Vaughan et al. 2003).

    Returns (fvar, fvar_err). fvar is NaN if there are <5 points or the excess
    variance is not positive (no intrinsic variability detected above noise).
    """
    x = np.asarray(flux, float)
    s = np.asarray(err, float)
    n = len(x)
    if n < 5:
        return np.nan, np.nan
    mean = x.mean()
    s2 = x.var(ddof=1)
    mse = np.mean(s ** 2)
    excess = s2 - mse
    if excess <= 0 or mean <= 0:
        return np.nan, np.nan
    fvar = math.sqrt(excess) / mean
    term1 = math.sqrt(1.0 / (2 * n)) * mse / (mean ** 2 * fvar)
    term2 = math.sqrt(mse / n) / mean
    return fvar, math.hypot(term1, term2)


def constancy_test(flux, err):
    """chi-square test against a constant (weighted-mean) flux -> (chi2, dof, p_value)."""
    from scipy.stats import chi2 as chi2dist

    x = np.asarray(flux, float)
    s = np.asarray(err, float)
    if len(x) < 3:
        return np.nan, 0, np.nan
    w = 1.0 / s ** 2
    mu = np.sum(w * x) / np.sum(w)
    chi2 = float(np.sum(((x - mu) / s) ** 2))
    dof = len(x) - 1
    return chi2, dof, float(chi2dist.sf(chi2, dof))


def flare_stats(flux, factor: float = 2.0):
    """Operational flare statistics on detected bins.

    duty_cycle      = fraction of bins with flux > factor * median flux
    peak_to_median  = max flux / median flux
    """
    x = np.asarray(flux, float)
    if len(x) == 0:
        return np.nan, np.nan
    med = np.median(x)
    return float(np.mean(x > factor * med)), float(x.max() / med)


def lightcurve_stats(lc, ts_min: float = 4.0, flare_factor: float = 2.0) -> dict:
    """Variability summary from a light-curve dict/npz with flux, flux_err, ts."""
    key = "flux" if "flux" in lc else "eflux"
    ekey = key + "_err"
    flux = np.asarray(lc[key], float)
    err = np.asarray(lc[ekey], float)
    ts = np.asarray(lc["ts"], float)
    ok = np.isfinite(flux) & np.isfinite(err) & np.isfinite(ts) & (ts >= ts_min) & (err > 0) & (flux > 0)
    out = dict(lc_flux_key=key, lc_nbins=int(len(flux)), lc_ndet=int(ok.sum()))
    if ok.sum() < 5:
        out.update(fvar=np.nan, fvar_err=np.nan, chi2=np.nan, chi2_dof=0, chi2_p=np.nan,
                   duty_cycle=np.nan, peak_to_median=np.nan, lc_mean=np.nan)
        return out
    f, e = flux[ok], err[ok]
    fvar, fvar_err = fractional_variability(f, e)
    chi2, dof, p = constancy_test(f, e)
    duty, p2m = flare_stats(f, flare_factor)
    out.update(fvar=fvar, fvar_err=fvar_err, chi2=chi2, chi2_dof=dof, chi2_p=p,
               duty_cycle=duty, peak_to_median=p2m, lc_mean=float(f.mean()))
    return out


# ----------------------------------------------------------------------------
# fermipy configuration
# ----------------------------------------------------------------------------
def find_data(src_slug: str, settings: dict):
    """Return (events_list_path, spacecraft_file) for a source, creating events.txt."""
    root = Path(settings["_root"])
    d = root / settings["paths"]["data_dir"] / src_slug
    ph = sorted(d.glob("*_PH*.fits"))
    sc = sorted(d.glob("*_SC*.fits"))
    if not ph or not sc:
        raise FileNotFoundError(
            f"No LAT data in {d}. Expected *_PH*.fits photon files and a *_SC*.fits "
            "spacecraft file (run download_lat_data.py or use the LAT data server)."
        )
    ev = d / "events.txt"
    ev.write_text("\n".join(str(p) for p in ph) + "\n")
    return ev, sc[0]


def build_fermipy_config(src: dict, settings: dict, outdir: Path, evfile: Path, scfile: Path) -> dict:
    t, e, r = settings["time"], settings["energy"], settings["roi"]
    sel, dif = settings["selection"], settings["diffuse"]
    return {
        "logging": {"verbosity": 3},
        "fileio": {"outdir": str(outdir)},
        "data": {"evfile": str(evfile), "scfile": str(scfile)},
        "binning": {"roiwidth": r["roiwidth"], "binsz": r["binsz"], "binsperdec": r["binsperdec"]},
        "selection": {
            "emin": e["emin"], "emax": e["emax"], "zmax": sel["zmax"],
            "evclass": sel["evclass"], "evtype": sel["evtype"],
            "tmin": round(iso_to_met(t["start"])), "tmax": round(iso_to_met(t["end"])),
            "ra": round(src["ra"], 4), "dec": round(src["dec"], 4),
        },
        "gtlike": {"edisp": True, "irfs": settings["irfs"], "edisp_disable": ["isodiff"]},
        "model": {
            "src_roiwidth": r["src_roiwidth"],
            "galdiff": dif["galdiff"], "isodiff": dif["isodiff"],
            "catalogs": [settings["paths"]["catalog"]],
        },
    }


def to_jsonable(o):
    """Recursively convert numpy types / NaN so json.dump works."""
    if isinstance(o, dict):
        return {str(k): to_jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return to_jsonable(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def dump_json(obj, path):
    with open(path, "w") as fh:
        json.dump(to_jsonable(obj), fh, indent=2)
