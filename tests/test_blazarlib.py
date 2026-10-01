"""Run with:  python tests/test_blazarlib.py   (or pytest). No Fermi software needed."""
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import blazarlib as bl  # noqa: E402


def test_met_roundtrip_and_mission_start():
    # Fermi mission start: MJD 54682.65603 (TT) <-> MET 239557417
    assert abs(bl.mjd_to_met(54682.65603) - 239557417) < 2
    assert abs(bl.met_to_mjd(bl.mjd_to_met(55555.5)) - 55555.5) < 1e-9
    # UTC date -> MET is within ~1 minute of the true value, fine for multi-year windows
    assert abs(bl.iso_to_met("2008-08-04T15:43:36") - 239557417) < 120


def test_names():
    assert bl.normalize_4fgl("4FGLJ2253.9+1609") == "4FGL J2253.9+1609"
    assert bl.normalize_4fgl("J1256.1-0547") == "4FGL J1256.1-0547"
    assert bl.normalize_4fgl("3C 279") is None
    assert bl.slugify("4C +21.35") == "4c_21_35"
    assert bl.class_group("fsrq") == "FSRQ" and bl.class_group("BLL") == "BL Lac"
    assert bl.class_group("") == "Unassociated" and bl.class_group("PSR") == "Galactic"


def test_fvar_recovers_injected_variability():
    rng = np.random.default_rng(0)
    n, true_fvar = 4000, 0.5
    mean = 1.0
    x = rng.normal(mean, true_fvar * mean, n)
    err = np.full(n, 0.1)
    obs = x + rng.normal(0, err)
    fvar, ferr = bl.fractional_variability(obs, err)
    assert abs(fvar - true_fvar) < 0.03, fvar
    assert ferr > 0


def test_fvar_constant_source_gives_nan():
    rng = np.random.default_rng(1)
    err = np.full(200, 0.2)
    obs = 1.0 + rng.normal(0, err)
    fvar, _ = bl.fractional_variability(obs, err)
    assert np.isnan(fvar) or fvar < 0.1
    chi2, dof, p = bl.constancy_test(obs, err)
    assert p > 0.01 and dof == 199


def test_flare_stats():
    f = np.array([1, 1, 1, 1, 5, 1, 1, 6, 1, 1.0])
    duty, p2m = bl.flare_stats(f, 2.0)
    assert duty == 0.2 and p2m == 6.0


def test_lightcurve_stats_masks_nondetections():
    rng = np.random.default_rng(2)
    n = 100
    flux = np.exp(rng.normal(0, 0.5, n))
    lc = dict(flux=flux, flux_err=0.05 * flux, ts=np.where(np.arange(n) < 10, 1.0, 50.0))
    st = bl.lightcurve_stats(lc, ts_min=4.0)
    assert st["lc_ndet"] == 90 and st["lc_nbins"] == 100 and st["fvar"] > 0.3


def test_catalog_lookup_and_nearest():
    from astropy.io import fits

    cols = [
        fits.Column(name="Source_Name", format="20A", array=np.array(["4FGL J2253.9+1609", "4FGL J1256.1-0547"])),
        fits.Column(name="RAJ2000", format="E", array=np.array([343.49, 194.04])),
        fits.Column(name="DEJ2000", format="E", array=np.array([16.15, -5.79])),
        fits.Column(name="CLASS1", format="10A", array=np.array(["fsrq", "FSRQ"])),
        fits.Column(name="ASSOC1", format="20A", array=np.array(["3C 454.3", "3C 279"])),
        fits.Column(name="SpectrumType", format="20A", array=np.array(["LogParabola", "LogParabola"])),
        fits.Column(name="Signif_Avg", format="E", array=np.array([300.0, 250.0])),
        fits.Column(name="Conf_95_SemiMajor", format="E", array=np.array([0.01, 0.01])),
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name="LAT_Point_Source_Catalog")
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "cat.fit"
        fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(p)
        cat = bl.Catalog(p)
        assert cat.lookup("4FGL J2253.9+1609")["assoc"] == "3C 454.3"
        assert cat.lookup("4FGL J0000.0+0000", assoc="3c 279")["fgl_name"] == "4FGL J1256.1-0547"
        assert cat.lookup("4FGL J0000.0+0000") is None
        row, sep = cat.nearest(194.1, -5.8)
        assert row["fgl_name"] == "4FGL J1256.1-0547" and sep < 0.2


def test_config_builder_has_required_keys():
    s = {
        "_root": ".", "paths": {"catalog": "cat.fit"}, "time": {"start": "2008-08-04", "end": "2018-08-04"},
        "energy": {"emin": 100, "emax": 300000},
        "roi": {"roiwidth": 12.0, "src_roiwidth": 16.0, "binsz": 0.1, "binsperdec": 8},
        "selection": {"zmax": 90, "evclass": 128, "evtype": 3}, "irfs": "P8R3_SOURCE_V3",
        "diffuse": {"galdiff": "g.fits", "isodiff": "i.txt"},
    }
    cfg = bl.build_fermipy_config({"ra": 10.0, "dec": -5.0}, s, Path("out"), Path("ev.txt"), Path("sc.fits"))
    assert cfg["selection"]["tmax"] > cfg["selection"]["tmin"] > 2e8
    for k in ("data", "binning", "selection", "gtlike", "model"):
        assert k in cfg


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print("ok  ", fn.__name__)
    print(f"\n{len(fns)} tests passed")
