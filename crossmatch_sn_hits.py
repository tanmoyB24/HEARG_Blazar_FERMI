#!/usr/bin/env python
"""
Quantify the blazar contamination of your original supernova search.

Input : a CSV of the gamma-ray excesses your SN analysis flagged, with columns
        ra, dec (deg, J2000) and optionally name, ts.
Output: for each hit, the nearest 4FGL source, its separation and class, plus a
        summary of what fraction of hits are blazars. This turns the "problem"
        (SN searches keep finding blazars) into the motivation slide of the project.

    python crossmatch_sn_hits.py my_sn_hits.csv --radius 0.2

A hit is matched to the nearest 4FGL source if the separation is below
max(--radius, that source's 95% positional error semi-major axis). Check by eye
that the matches make sense in crowded fields.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import blazarlib as bl


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("hits", help="CSV with columns ra, dec [, name, ts]")
    ap.add_argument("--settings", default="settings.yaml")
    ap.add_argument("--radius", type=float, default=0.2, help="minimum match radius in deg")
    ap.add_argument("--out", default="tables")
    args = ap.parse_args()

    s = bl.load_settings(args.settings)
    cat = bl.Catalog(s["paths"]["catalog"])
    hits = pd.read_csv(args.hits)
    if not {"ra", "dec"} <= set(hits.columns):
        raise SystemExit("hits CSV needs ra and dec columns (degrees)")

    rows = []
    for i, h in hits.iterrows():
        near, sep = cat.nearest(h["ra"], h["dec"])
        matched = sep <= max(args.radius, near["r95"] if near["r95"] == near["r95"] else 0.0)
        group = bl.class_group(near["cat_class"]) if matched else "No 4FGL match"
        rows.append(dict(name=h.get("name", f"hit_{i}"), ra=h["ra"], dec=h["dec"],
                         ts=h.get("ts", float("nan")), nearest=near["fgl_name"] if matched else "",
                         sep_deg=round(sep, 3), cat_class=near["cat_class"] if matched else "",
                         assoc=near["assoc"] if matched else "", group=group))
    res = pd.DataFrame(rows)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    res.to_csv(out / "sn_hits_crossmatch.csv", index=False)

    counts = res["group"].value_counts()
    n = len(res)
    print(f"{n} hits\n")
    for g, c in counts.items():
        print(f"  {g:<15} {c:>4}  ({100 * c / n:.0f}%)")
    blz = counts.reindex(["FSRQ", "BL Lac", "BCU"]).fillna(0).sum()
    print(f"\nBlazar-class (FSRQ + BL Lac + BCU): {int(blz)}/{n} = {100 * blz / n:.0f}%")

    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    colors = {"FSRQ": "#c0392b", "BL Lac": "#2471a3", "BCU": "#8e44ad"}
    ax.barh(counts.index[::-1], counts.values[::-1], color=[colors.get(g, "#95a5a6") for g in counts.index[::-1]])
    ax.set_xlabel("Number of SN-search gamma-ray hits")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "sn_hits_classes.png", dpi=220)
    print(f"\nwrote {out}/sn_hits_crossmatch.csv and sn_hits_classes.png")


if __name__ == "__main__":
    main()
