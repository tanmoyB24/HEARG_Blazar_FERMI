#!/usr/bin/env python
"""
Download Fermi-LAT photon + spacecraft files for each source (OPTIONAL helper).

Uses astroquery.fermi, which talks to the public LAT Data Server. This script is
a convenience wrapper and has NOT been tested against the live server by the
project author - if it fails, fall back to the web form:
    https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi
and put the downloaded *_PH*.fits and *_SC*.fits files into data/<source_slug>/
(the slug is the lower-case name with non-alphanumerics replaced by '_', e.g. 3c_454_3).

Query settings to use on the web form: search radius 20 deg (must exceed the ROI
half-diagonal), energy 100-300000 MeV, the same dates as settings.yaml, and tick
"spacecraft data".

    python download_lat_data.py --only "3C 279"
    python download_lat_data.py            # all sources
"""
import argparse
import time
import urllib.request
from pathlib import Path

import blazarlib as bl


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--settings", default="settings.yaml")
    ap.add_argument("--only", help="only this source name")
    ap.add_argument("--radius", type=float, default=20.0, help="search radius in deg")
    args = ap.parse_args()

    try:
        from astroquery.fermi import FermiLAT
    except ImportError:
        raise SystemExit("pip install astroquery  (or use the LAT web form, see --help)")

    s = bl.load_settings(args.settings)
    cat = bl.Catalog(s["paths"]["catalog"])
    srcs = bl.read_sources(f'{s["_root"]}/{s["paths"]["sources"]}')
    t0, t1 = s["time"]["start"], s["time"]["end"]

    for _, r in srcs.iterrows():
        if args.only and r["name"] != args.only:
            continue
        hit = cat.lookup(r["fgl_name"], r.get("assoc"))
        if hit is None:
            print(f"skip {r['name']}: not in catalog")
            continue
        slug = bl.slugify(r["name"])
        out = Path(s["_root"]) / s["paths"]["data_dir"] / slug
        if list(out.glob("*_PH*.fits")):
            print(f"{r['name']}: data already present, skipping")
            continue
        out.mkdir(parents=True, exist_ok=True)
        print(f"{r['name']}: querying LAT data server ...")
        try:
            urls = FermiLAT.query_object(
                f"{hit['ra']},{hit['dec']}",
                coordsystem="J2000",
                searchradius=args.radius,
                obsdates=f"{t0} 00:00:00, {t1} 00:00:00",
                timesystem="Gregorian",
                energyrange_MeV="100, 300000",
                spacecraftdata="on",
            )
        except Exception as exc:  # network / server problems
            print(f"  query failed: {exc}")
            continue
        for u in urls:
            dest = out / u.split("/")[-1]
            print(f"  -> {dest.name}")
            urllib.request.urlretrieve(u, dest)
        time.sleep(2)  # be polite to the server


if __name__ == "__main__":
    main()
