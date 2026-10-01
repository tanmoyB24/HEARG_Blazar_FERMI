#!/usr/bin/env python
"""
Check sources.csv against your 4FGL catalog file BEFORE running anything heavy.

For each row it prints the catalog coordinates, class, association and spectral
model, and warns when the name is missing or the catalog class disagrees with the
class you assigned. Writes sources_resolved.csv with the coordinates.

    python verify_sources.py
"""
import argparse
import sys

import pandas as pd

import blazarlib as bl


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--settings", default="settings.yaml")
    args = ap.parse_args()

    s = bl.load_settings(args.settings)
    cat = bl.Catalog(s["paths"]["catalog"])
    srcs = bl.read_sources(f'{s["_root"]}/{s["paths"]["sources"]}')
    print(f"Catalog: {cat.path}  ({len(cat)} sources)\n")

    rows, problems = [], 0
    for _, r in srcs.iterrows():
        hit = cat.lookup(r["fgl_name"], r.get("assoc"))
        if hit is None:
            print(f"[MISSING] {r['name']:<14} {r['fgl_name']}  -> not in catalog. Fix the name in sources.csv.")
            problems += 1
            continue
        flags = []
        if hit["fgl_name"] != bl.normalize_4fgl(r["fgl_name"]):
            flags.append(f"name resolved via association -> use {hit['fgl_name']}")
        want = {"FSRQ": "fsrq", "BLL": "bll"}.get(r["class"], r["class"].lower())
        if hit["cat_class"].lower() != want:
            flags.append(f"class mismatch: you said {r['class']}, catalog says '{hit['cat_class']}'")
        status = "WARN" if flags else "ok"
        problems += bool(flags)
        print(f"[{status:>4}] {r['name']:<14} {hit['fgl_name']:<20} RA={hit['ra']:8.3f} Dec={hit['dec']:8.3f} "
              f"class={hit['cat_class']:<5} assoc={hit['assoc']:<16} model={hit['spectrum_type']:<14} "
              f"sig={hit['signif']:.0f}")
        for f in flags:
            print(f"        ! {f}")
        rows.append({**r.to_dict(), **hit})

    out = f'{s["_root"]}/sources_resolved.csv'
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"\nWrote {out} ({len(rows)}/{len(srcs)} resolved, {problems} issue(s)).")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
