# Gamma-ray variability and spectral properties of Fermi-LAT blazars (FSRQs vs BL Lacs)

A conference-oriented project that reuses the HEARG Tutorial 2 (Fermi tools) and Tutorial 3 (fermipy)
workflow, but for a **sample of blazars** instead of a single source or a supernova search.

**Science question.** Do FSRQs and BL Lacs differ in their GeV spectra (photon index, spectral
curvature) and in their variability (fractional variability, flare duty cycle) over ~10 years of
Fermi-LAT data?

**Motivating story.** Your supernova search kept returning blazars. `crossmatch_sn_hits.py` quantifies that
("X% of our SN-search hits are blazars"), which becomes the opening slide: blazars dominate the
extragalactic gamma-ray sky, so they must be understood before any faint transient can be claimed.

---

## Status: what is tested and what is not

| Part | Status |
|---|---|
| `blazarlib.py` (time conversion, catalog, F_var, flare stats, config builder) | Unit-tested (`python tests/test_blazarlib.py`) |
| `analyze_population.py`, figures, tables, statistics | Tested end-to-end on synthetic data |
| `verify_sources.py`, `crossmatch_sn_hits.py`, `run_fermipy.py --dry-run` | Tested with a fake catalog |
| **`run_fermipy.py` fermipy fitting, the notebook 01 fermipy cells, `download_lat_data.py`** | **Not run against real Fermi data** (no Fermi Science Tools were available when this was written). The calls follow the fermipy documentation and your Tutorial 3, but **run ONE source first** and expect to adjust a keyword or two. |

Every optional step in `run_fermipy.py` is wrapped so a failure is recorded in `results/<slug>/summary.json["status"]`
instead of killing the run.

---

## Folder layout

```
settings.yaml            all knobs: catalog path, dates, ROI, energies, bin size, thresholds
sources.csv              the sample (name, 4FGL name, class)
blazarlib.py             shared helpers
verify_sources.py        step 1: check the sample against YOUR 4FGL file
download_lat_data.py     step 2 (optional): fetch photon + spacecraft files
run_fermipy.py           step 3: the per-source pipeline
analyze_population.py    step 4: table, statistics, figures
crossmatch_sn_hits.py    motivation: blazar fraction in your SN search hits
notebooks/               01 single-source walkthrough, 02 population results
tools/                   synthetic-data generator (for testing/layout only)
tests/                   unit tests
data/<slug>/             LAT files per source   (you provide)
results/<slug>/          pipeline output        (generated)
```

## Quick start

```bash
# 0. environment (or reuse your Tutorial docker image)
conda env create -f environment.yml && conda activate blazar-fermi

# 1. edit settings.yaml: set paths.catalog to the gll_psc_*.fit you used in Tutorial 2
python verify_sources.py              # fix any [MISSING]/[WARN] rows in sources.csv

# 2. data: one folder per source, e.g. data/3c_279/ containing *_PH*.fits and *_SC*.fits
python download_lat_data.py --only "3C 279"     # or use the LAT web form (see script docstring)

# 3. test ONE source end-to-end, ideally in notebook 01 so you can inspect each step
python run_fermipy.py --dry-run --only "3C 279"  # sanity-check the generated config
python run_fermipy.py --only "3C 279"

# 4. when happy, run the sample (resumes automatically; Ctrl-C safe)
python run_fermipy.py
python analyze_population.py --example "3C 279"
```

**Preview everything with fake data first** (figures are watermarked "SYNTHETIC DEMO"):

```bash
python tools/make_synthetic_results.py --out results_demo
python analyze_population.py --results results_demo --out demo_out
```

**Quantify the motivation** (optional but recommended):

```bash
# my_sn_hits.csv: ra,dec[,name,ts] of the gamma-ray sources your SN search flagged
python crossmatch_sn_hits.py my_sn_hits.csv --radius 0.2
```

---

## What the pipeline measures

| Quantity | Definition |
|---|---|
| `ts`, `flux`, `index` | Whole-period fit: test statistic, 0.1-300 GeV photon flux, photon index Gamma (dN/dE ~ E^-Gamma) |
| `ts_curv` | 2 (lnL_LogParabola - lnL_PowerLaw); source left in the preferred model (threshold `curvature_ts_threshold`, default 4) |
| `fvar`, `fvar_err` | Fractional variability from the excess variance of the light curve (Vaughan et al. 2003), detected bins only (TS >= `ts_min`) |
| `chi2_p` | p-value of a constant-flux fit; small = significantly variable |
| `duty_cycle` | Fraction of detected bins with flux > `flare_factor` x median (default 2x). This is an *operational* definition; state it on the poster |
| `peak_to_median` | max flux / median flux |
| `loc_offset`, `loc_r95` | Gamma-ray position vs catalog position and its 95% error |

## Choices you should be ready to defend

* **One catalog, one time range for every source.** Set once in `settings.yaml`.
* **Monthly bins (default)** keep most bright blazars above TS = 4 per bin, so F_var is well defined. Weekly bins
  need many more CPU hours and give more upper limits for BL Lacs.
* **F_var uses detected bins only**, so it is biased low for faint sources. Say so, or restrict the F_var
  comparison to sources with > ~50 detected bins.
* **Sample selection bias.** The sources are chosen because they are bright, so conclusions apply to bright
  Fermi blazars, not "all blazars".
* **Small N.** With 8 vs 8 sources, Mann-Whitney p-values are indicative. Increase the sample if time allows.
* **Systematics** (diffuse model, effective area at the few-percent level) are not propagated.

## Compute budget (honest expectation)

Each source needs its own ROI setup (livetime cube, exposure, source maps) and a light curve that refits every
bin. On a laptop expect from tens of minutes to several hours per source depending on time span and CPUs; the
light curve dominates. Start with 6-8 sources (3-4 per class), shorten `time.end` if needed, and scale up later.
`run_fermipy.py` skips finished sources, so you can run it overnight and resume.

## Suggested timeline (8 weeks)

| Week | Goal |
|---|---|
| 1 | Environment, `verify_sources.py`, data for 1 source, run notebook 01 end-to-end |
| 2 | Debug and freeze `settings.yaml`; check residual/TS maps look clean for 2-3 sources |
| 3-4 | Run the full sample; inspect each source's diagnostics; drop/replace failed sources |
| 5 | `analyze_population.py`; decide the key result and the figures to show |
| 6 | Optional extras (see below); SN-hit cross-match for the motivation slide |
| 7 | Draft abstract/poster/talk (`CONFERENCE.md`); show a colleague |
| 8 | Rehearse, prepare backup slides and answers to expected questions |

## Optional extensions

* Spectral hardening with flux: free the index in each light-curve bin and plot Gamma vs flux.
* Weekly or daily bins for one or two flaring sources (flare case study).
* Larger sample from 4LAC, selected by a clear rule (e.g. top N by `Signif_Avg` per class).
* Power spectral density / structure function of the light curves.
* Add multiwavelength context (e.g. radio flux) from the literature.

## Credits

Workflow adapted from the HEARG (CAM-SUST) Tutorials 2 and 3. Built on `fermipy` (Wood et al. 2017; see
`CONFERENCE.md` for references) and the NASA Fermi Science Tools.
