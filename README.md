# Gamma-ray Variability and Spectral Properties of Fermi-LAT Blazars

A Fermi-LAT analysis pipeline for studying the **gamma-ray spectra and long-term variability of blazars**, with a focus on comparing **FSRQs** and **BL Lacs**.

The workflow uses the **Fermi Science Tools** and **fermipy** to perform spectral and light-curve analyses on a population of sources.

## Scientific Goal

The main question is:

> **Do FSRQs and BL Lacs differ in their GeV spectra and long-term gamma-ray variability?**

The analysis measures:

* Photon index
* Spectral curvature
* Gamma-ray flux
* Fractional variability (`F_var`)
* Variability significance
* Flare duty cycle
* Peak-to-median flux ratio

The default analysis covers approximately **10 years of Fermi-LAT data** using a consistent energy range and time interval.

---

## Workflow

```text
sources.csv
     │
     ▼
verify_sources.py
     │
     ▼
Fermi-LAT data
     │
     ▼
run_fermipy.py
     │
     ├── Spectral analysis
     └── Light-curve analysis
              │
              ▼
       Source results
              │
              ▼
   analyze_population.py
              │
              ▼
     FSRQ vs BL Lac
```

## Repository Structure

```text
.
├── settings.yaml
├── sources.csv
├── blazarlib.py
├── verify_sources.py
├── download_lat_data.py
├── run_fermipy.py
├── analyze_population.py
├── crossmatch_sn_hits.py
│
├── notebooks/
│   ├── 01_single_source_analysis.ipynb
│   └── 02_population_analysis.ipynb
│
├── tools/
│   └── make_synthetic_results.py
├── tests/
│   └── test_blazarlib.py
├── data/
└── results/
```

### Main Components

| File                    | Description                                |
| ----------------------- | ------------------------------------------ |
| `settings.yaml`         | Analysis configuration                     |
| `sources.csv`           | Blazar sample and classifications          |
| `blazarlib.py`          | Shared analysis and statistical utilities  |
| `verify_sources.py`     | Validates sources against the 4FGL catalog |
| `download_lat_data.py`  | Prepares/downloads LAT data                |
| `run_fermipy.py`        | Per-source Fermi-LAT/fermipy analysis      |
| `analyze_population.py` | Population statistics and figures          |
| `crossmatch_sn_hits.py` | Optional source cross-matching             |
| `tools/`                | Synthetic-data and testing utilities       |
| `tests/`                | Unit tests                                 |

---

## Requirements

* Python 3
* Fermi Science Tools
* `fermipy`
* NumPy
* SciPy
* pandas
* Matplotlib
* Astropy
* PyYAML
* pytest

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate blazar-fermi
```

A working **Fermi Science Tools + fermipy** installation is required for analysis of real LAT data.

---

## Usage

### 1. Verify the source sample

```bash
python verify_sources.py
```

### 2. Prepare LAT data

For example:

```bash
python download_lat_data.py --only "3C 279"
```

Each source should have its own directory under:

```text
data/<source_slug>/
```

containing the required LAT photon and spacecraft files.

### 3. Test one source

First check the generated configuration:

```bash
python run_fermipy.py --dry-run --only "3C 279"
```

Then run the analysis:

```bash
python run_fermipy.py --only "3C 279"
```

Once validated, process the full sample:

```bash
python run_fermipy.py
```

### 4. Population analysis

```bash
python analyze_population.py
```

---

## Main Measurements

| Quantity         | Description                                         |
| ---------------- | --------------------------------------------------- |
| `flux`           | Integrated gamma-ray photon flux                    |
| `index`          | Photon index                                        |
| `ts_curv`        | Spectral curvature test statistic                   |
| `fvar`           | Fractional variability                              |
| `chi2_p`         | Constant-flux test probability                      |
| `duty_cycle`     | Fraction of detected bins above the flare threshold |
| `peak_to_median` | Maximum/median flux                                 |

The default light-curve binning is **monthly**.

---

## Important Considerations

* The same catalog, energy range, time interval, and analysis settings should be used for all sources.
* `F_var` is calculated from detected bins and may be biased for faint sources.
* Bright-source selection means the sample may not represent the entire blazar population.
* Monthly bins provide a balance between time resolution and detection significance.
* Run and validate **one source first** before processing the complete sample.

---

## Testing

Run the unit tests with:

```bash
pytest
```

Synthetic data can be used to test the population-analysis workflow without real Fermi data:

```bash
python tools/make_synthetic_results.py --out results_demo
python analyze_population.py --results results_demo --out demo_out
```

Synthetic results are intended only for **software testing**, not scientific conclusions.

---

## References

* Wood et al. (2017) — `fermipy`
* Vaughan et al. (2003) — Fractional variability
* Fermi-LAT 4FGL/4LAC catalogs
* NASA Fermi Science Tools

The workflow is adapted from the HEARG (CAM-SUST) Fermi-LAT analysis tutorials.
