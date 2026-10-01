#!/usr/bin/env python
"""
Population analysis of the blazar sample: table, FSRQ-vs-BL Lac statistics and
conference-ready figures. Reads results/<slug>/{summary.json,sed.npz,lightcurve.npz}.

    python analyze_population.py                 # -> figures/ and tables/
    python analyze_population.py --results results_demo --out demo_out

If any summary.json is flagged "synthetic": true (see tools/make_synthetic_results.py),
every figure carries a SYNTHETIC DEMO watermark so fake data cannot be shown by accident.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import blazarlib as bl

COLORS = {"FSRQ": "#c0392b", "BLL": "#2471a3"}
LABELS = {"FSRQ": "FSRQ", "BLL": "BL Lac"}
plt.rcParams.update({
    "font.size": 11, "axes.labelsize": 12, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 110, "savefig.dpi": 220, "savefig.bbox": "tight", "legend.frameon": False,
})


# ----------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------
def load_table(results_dir, ts_min=4.0, flare_factor=2.0) -> pd.DataFrame:
    rows = []
    for d in sorted(Path(results_dir).iterdir()):
        f = d / "summary.json"
        if not f.is_file():
            continue
        s = json.loads(f.read_text())
        lcf = d / "lightcurve.npz"
        if lcf.exists():  # recompute from the stored light curve so ts_min/flare_factor can be changed here
            s.update(bl.lightcurve_stats(np.load(lcf), ts_min, flare_factor))
        s["slug"] = d.name
        rows.append(s)
    if not rows:
        raise SystemExit(f"No summary.json files found in {results_dir}. Run run_fermipy.py first.")
    df = pd.DataFrame(rows)
    for c in ("flux", "index", "ts_curv", "fvar", "duty_cycle", "ts", "beta"):
        if c not in df:
            df[c] = np.nan
    return df


def watermark(fig, synthetic):
    if synthetic:
        fig.text(0.5, 0.5, "SYNTHETIC DEMO - NOT REAL DATA", rotation=25, ha="center", va="center",
                 fontsize=28, color="red", alpha=0.18, weight="bold")


def save(fig, out, name, synthetic):
    watermark(fig, synthetic)
    p = Path(out) / "figures" / name
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p)
    plt.close(fig)
    return p


def _legend_handles(ax):
    h, l = ax.get_legend_handles_labels()
    seen = dict(zip(l, h))
    ax.legend(seen.values(), seen.keys())


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
def fig_index_flux(df, out, syn):
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    for c, g in df.groupby("cls"):
        ax.errorbar(g["flux"], g["index"], yerr=g["index_err"], fmt="o", ms=6, capsize=2,
                    color=COLORS.get(c, "gray"), label=LABELS.get(c, c))
        for _, r in g.iterrows():
            ax.annotate(r["name"], (r["flux"], r["index"]), fontsize=6.5, xytext=(3, 3),
                        textcoords="offset points", color=COLORS.get(c, "gray"))
    ax.set_xscale("log")
    ax.set_xlabel(r"Flux 0.1-300 GeV  [ph cm$^{-2}$ s$^{-1}$]")
    ax.set_ylabel(r"Photon index $\Gamma$")
    _legend_handles(ax)
    return save(fig, out, "fig1_index_vs_flux.png", syn)


def _strip(ax, df, col, ylabel, log=False):
    groups = [c for c in ("FSRQ", "BLL") if c in set(df["cls"])]
    for i, c in enumerate(groups):
        v = df.loc[df["cls"] == c, col].dropna().values
        if len(v) == 0:
            continue
        x = i + np.random.default_rng(1).uniform(-0.12, 0.12, len(v))
        ax.scatter(x, v, s=36, color=COLORS[c], alpha=0.85, zorder=3)
        ax.hlines(np.median(v), i - 0.3, i + 0.3, color="k", lw=2, zorder=4)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([LABELS[c] for c in groups])
    ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.set_ylabel(ylabel)
    if log:
        ax.set_yscale("log")


def fig_class_comparison(df, out, syn):
    fig, axs = plt.subplots(1, 4, figsize=(12.5, 3.6))
    _strip(axs[0], df, "index", r"Photon index $\Gamma$")
    _strip(axs[1], df, "fvar", r"Fractional variability $F_{\rm var}$")
    _strip(axs[2], df, "ts_curv", r"$TS_{\rm curv}$", log=False)
    axs[2].axhline(4, ls="--", lw=0.8, color="gray")
    _strip(axs[3], df, "duty_cycle", "Flare duty cycle")
    fig.suptitle("Black bar = median", fontsize=9, y=1.0)
    fig.tight_layout()
    return save(fig, out, "fig2_class_comparison.png", syn)


def fig_fvar_index(df, out, syn):
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    for c, g in df.groupby("cls"):
        ax.errorbar(g["index"], g["fvar"], xerr=g["index_err"], yerr=g["fvar_err"], fmt="o", capsize=2,
                    color=COLORS.get(c, "gray"), label=LABELS.get(c, c))
    ax.set_xlabel(r"Photon index $\Gamma$")
    ax.set_ylabel(r"$F_{\rm var}$")
    _legend_handles(ax)
    return save(fig, out, "fig3_fvar_vs_index.png", syn)


def _grid(n, ncols=4, w=3.3, h=2.5):
    nrows = int(np.ceil(n / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(w * ncols, h * nrows), squeeze=False)
    return fig, axs.ravel()


def fig_sed_gallery(df, results, out, syn, ts_min=4.0):
    rows = [(r, Path(results) / r["slug"] / "sed.npz") for _, r in df.iterrows()]
    rows = [(r, p) for r, p in rows if p.exists()]
    if not rows:
        return None
    fig, axs = _grid(len(rows))
    for ax, (r, p) in zip(axs, rows):
        d = np.load(p)
        e, lo, hi = d["e_ref"], d["e_ref"] - d["e_min"], d["e_max"] - d["e_ref"]
        det = d["ts"] >= ts_min
        col = COLORS.get(r["cls"], "gray")
        ax.errorbar(e[det], d["e2dnde"][det], xerr=[lo[det], hi[det]], yerr=d["e2dnde_err"][det],
                    fmt="o", ms=4, color=col, capsize=0)
        ax.errorbar(e[~det], d["e2dnde_ul95"][~det], xerr=[lo[~det], hi[~det]],
                    yerr=0.35 * d["e2dnde_ul95"][~det], uplims=True, fmt="none", color=col, alpha=0.7)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f'{r["name"]} ({LABELS.get(r["cls"], r["cls"])})', fontsize=9)
    for ax in axs[len(rows):]:
        ax.axis("off")
    for ax in axs[:len(rows)]:
        ax.set_xlabel("E [MeV]", fontsize=8)
    axs[0].set_ylabel(r"E$^2$dN/dE [MeV cm$^{-2}$ s$^{-1}$]", fontsize=8)
    fig.tight_layout()
    return save(fig, out, "fig4_sed_gallery.png", syn)


def fig_lc_gallery(df, results, out, syn, ts_min=4.0):
    rows = [(r, Path(results) / r["slug"] / "lightcurve.npz") for _, r in df.iterrows()]
    rows = [(r, p) for r, p in rows if p.exists()]
    if not rows:
        return None
    fig, axs = _grid(len(rows))
    for ax, (r, p) in zip(axs, rows):
        d = np.load(p)
        k = "flux" if "flux" in d else "eflux"
        tc = 0.5 * (d["tmin_mjd"] + d["tmax_mjd"])
        det = d["ts"] >= ts_min
        col = COLORS.get(r["cls"], "gray")
        ax.errorbar(tc[det], d[k][det], yerr=d[k + "_err"][det], fmt="o", ms=2.5, lw=0.8, color=col)
        ax.plot(tc[~det], d[k + "_ul95"][~det], "v", ms=2.5, color=col, alpha=0.5)
        ax.set_title(f'{r["name"]}  ' + (rf'$F_{{var}}$={r["fvar"]:.2f}' if np.isfinite(r["fvar"]) else ""),
                     fontsize=9)
        ax.set_yscale("log")
    for ax in axs[len(rows):]:
        ax.axis("off")
    for ax in axs[:len(rows)]:
        ax.set_xlabel("MJD", fontsize=8)
    axs[0].set_ylabel("Flux [ph cm$^{-2}$ s$^{-1}$]", fontsize=8)
    fig.tight_layout()
    return save(fig, out, "fig5_lightcurve_gallery.png", syn)


def fig_single_lc(row, results, out, syn, ts_min=4.0, flare_factor=2.0):
    """Large light curve with flare threshold, for a talk's 'example source' slide."""
    p = Path(results) / row["slug"] / "lightcurve.npz"
    if not p.exists():
        return None
    d = np.load(p)
    k = "flux" if "flux" in d else "eflux"
    tc = 0.5 * (d["tmin_mjd"] + d["tmax_mjd"])
    det = d["ts"] >= ts_min
    col = COLORS.get(row["cls"], "gray")
    fig, ax = plt.subplots(figsize=(8.5, 3.6))
    ax.errorbar(tc[det], d[k][det], yerr=d[k + "_err"][det], fmt="o", ms=4, lw=1, color=col, label="detections")
    ax.plot(tc[~det], d[k + "_ul95"][~det], "v", color=col, alpha=0.5, label="95% upper limits")
    ax.axhline(flare_factor * np.median(d[k][det]), ls="--", color="k", lw=0.9,
               label=f"{flare_factor:g}x median (flare threshold)")
    ax.set_yscale("log")
    ax.set_xlabel("MJD")
    ax.set_ylabel(r"Flux [ph cm$^{-2}$ s$^{-1}$]")
    ax.set_title(f'{row["name"]} - Fermi-LAT light curve')
    ax.legend(fontsize=8)
    return save(fig, out, f'example_lc_{row["slug"]}.png', syn)


# ----------------------------------------------------------------------------
# Statistics + tables
# ----------------------------------------------------------------------------
def class_stats(df) -> str:
    lines = ["FSRQ vs BL Lac comparison (two-sided Mann-Whitney U; medians in brackets)", "-" * 72]
    a, b = df[df["cls"] == "FSRQ"], df[df["cls"] == "BLL"]
    lines.append(f"N(FSRQ)={len(a)}  N(BLL)={len(b)}")
    if len(a) < 3 or len(b) < 3:
        lines.append("Too few sources per class for a meaningful test (need >= 3, ideally >= 8).")
    for col, label in (("index", "photon index"), ("fvar", "F_var"), ("ts_curv", "TS_curv"),
                       ("duty_cycle", "flare duty cycle"), ("flux", "flux")):
        x, y = a[col].dropna(), b[col].dropna()
        if len(x) < 3 or len(y) < 3:
            lines.append(f"{label:<18} not enough values (FSRQ {len(x)}, BLL {len(y)})")
            continue
        u, p = stats.mannwhitneyu(x, y, alternative="two-sided")
        lines.append(f"{label:<18} FSRQ [{x.median():.3g}]  BLL [{y.median():.3g}]   U={u:.0f}  p={p:.3g}")
    ok = df[["index", "flux"]].dropna()
    if len(ok) >= 5:
        rho, p = stats.spearmanr(ok["flux"], ok["index"])
        lines.append(f"Spearman rho(index, flux) = {rho:.2f}  (p={p:.3g}, N={len(ok)})")
    ok = df[["index", "fvar"]].dropna()
    if len(ok) >= 5:
        rho, p = stats.spearmanr(ok["index"], ok["fvar"])
        lines.append(f"Spearman rho(index, F_var) = {rho:.2f}  (p={p:.3g}, N={len(ok)})")
    lines.append("")
    lines.append("With small N treat p-values as indicative only; say so on the poster.")
    return "\n".join(lines)


def write_tables(df, out):
    tdir = Path(out) / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    cols = ["name", "cls", "ts", "flux", "index", "index_err", "beta", "ts_curv", "model", "fvar", "fvar_err",
            "chi2_p", "duty_cycle", "peak_to_median", "lc_ndet", "lc_nbins"]
    cols = [c for c in cols if c in df.columns]
    df[cols].to_csv(tdir / "population_table.csv", index=False)

    def f(v, fmt):
        return "--" if (v is None or not np.isfinite(v)) else format(v, fmt)

    tex = [r"\begin{tabular}{llrrrrrr}", r"\hline",
           r"Source & Class & TS & Flux$^a$ & $\Gamma$ & $TS_{\rm curv}$ & $F_{\rm var}$ & Duty \\", r"\hline"]
    for _, r in df.iterrows():
        tex.append(f'{r["name"]} & {LABELS.get(r["cls"], r["cls"])} & {f(r["ts"], ".0f")} & '
                   f'{f(r["flux"] * 1e8, ".1f")} & {f(r["index"], ".2f")} & {f(r["ts_curv"], ".1f")} & '
                   f'{f(r["fvar"], ".2f")} & {f(r["duty_cycle"], ".2f")} \\\\')
    tex += [r"\hline", r"\multicolumn{8}{l}{$^a$ 0.1--300 GeV, units of $10^{-8}$ ph cm$^{-2}$ s$^{-1}$} \\",
            r"\end{tabular}"]
    (tdir / "population_table.tex").write_text("\n".join(tex) + "\n")
    (tdir / "class_statistics.txt").write_text(class_stats(df) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default=".")
    ap.add_argument("--ts-min", type=float, default=4.0, help="min TS for a light-curve bin to count as detected")
    ap.add_argument("--flare-factor", type=float, default=2.0)
    ap.add_argument("--example", help="source name for a large single light-curve figure")
    args = ap.parse_args()

    df = load_table(args.results, args.ts_min, args.flare_factor)
    syn = bool(df.get("synthetic", pd.Series([False])).fillna(False).any())
    made = [fig_index_flux(df, args.out, syn), fig_class_comparison(df, args.out, syn),
            fig_fvar_index(df, args.out, syn), fig_sed_gallery(df, args.results, args.out, syn, args.ts_min),
            fig_lc_gallery(df, args.results, args.out, syn, args.ts_min)]
    if args.example:
        row = df[df["name"] == args.example]
        if len(row):
            made.append(fig_single_lc(row.iloc[0], args.results, args.out, syn, args.ts_min, args.flare_factor))
    write_tables(df, args.out)
    print(f"{len(df)} sources analysed" + ("  [SYNTHETIC DEMO DATA]" if syn else ""))
    for m in made:
        if m:
            print("wrote", m)
    print(f"wrote {args.out}/tables/ (csv, tex, class_statistics.txt)")
    print()
    print((Path(args.out) / "tables" / "class_statistics.txt").read_text())


if __name__ == "__main__":
    main()
