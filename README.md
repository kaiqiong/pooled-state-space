# pooled-state-space

Replication code for

> **Pooled likelihood inference for covariate-dependent persistence in repeated short physiological recovery series**
> Kaiqiong Zhao and Xuan Quynh Vo

This package corresponds to the manuscript and supplement dated 13 October 2026. It contains the model and inference code, the simulation generators, the CGMacros application pipeline, the table and figure generators, a small test suite, and checked aggregate summaries of the published results. Raw participant data and large fit caches are **not** included; they are downloaded or regenerated locally.

Intended repository: `https://github.com/kaiqiong/pooled-state-space`

---

## Contents

| Path | What it holds |
|---|---|
| `R/` | Validated model, inference and mixed-model functions (Kalman filter and smoother, EM with the constrained M-step, per-episode observed information, participant-clustered sandwiches) |
| `scripts/application/` | CGMacros data preparation, fits and retained sensitivity checks |
| `scripts/simulation/` | Simulation generators, working-likelihood targets, summaries and numerical checks |
| `scripts/figures/` | Table and figure export, rendered from numeric inputs |
| `tests/` | Small independent numerical and estimation checks (no data needed) |
| `reference/` | Compact checked aggregate summaries of the published results (not fit caches) |
| `config/` | Tested package versions, frozen generating parameters, file-integrity manifest |
| `docs/` | Stage instructions, output map and validation record |
| `data/README.txt` | How to obtain the source data |
| `results/` | Created locally by the scripts; ignored by Git |

---

## Quick start

**Requirements.** R (tested with 4.3.2) and Python 3.7 or later. The numerical analyses run in R; Python's standard library provides the top-level runner and the data downloader. Only the illustrated raw-window figure needs `numpy` and `matplotlib`.

```bash
# 1. Check dependencies (add --install to install missing R packages from CRAN)
python3 run.py check --install

# 2. Optional: Python plotting packages for the raw-window figure
python3 -m pip install -r requirements.txt

# 3. Run the synthetic checks (about 15 seconds; no data download)
python3 run.py smoke

# 4. Recreate the published tables and figures from the shipped reference summaries
python3 run.py figures --source reference
```

Step 4 is a rendering and checking route, not a new Monte Carlo experiment. It does not include participant records or the raw-window figures. `docs/OUTPUTS.txt` and `docs/OUTPUT_MAP.csv` list exactly what is generated and from which inputs.

Required R packages: `numDeriv`, `statmod`, `nlme`, `ggplot2`, `patchwork`, `dplyr`, `tidyr`. Exact tested versions are recorded in `config/R-packages-tested.csv`; the installer does not pin an older CRAN snapshot. Floating-point and optimizer details can differ across platforms, R versions, BLAS libraries and package versions; use the recorded versions when investigating numerical or graphical differences.

---

## Reproducing the CGMacros application

The application uses CGMacros v1.0.0 from PhysioNet (CC BY-NC-SA 4.0; see `data/README.txt` for access, citation and checksum details). Participant records are downloaded locally and never committed.

```bash
# Download the official v1.0.0 archive (about 627 MB) and extract the 46 CSV files
python3 run.py data
# or use an already downloaded official ZIP (its SHA-256 is verified)
python3 run.py data --archive /path/to/CGMacros_dateshifted365.zip
```

Then run the stages in order from a fresh checkout:

```bash
python3 run.py application preprocess    # 65-minute episodes (1174 windows, 44 participants, 14 readings)
                                          # and 120-minute episodes (1059 windows, 25 readings)
python3 run.py application prepare       # positive fixed-variance fit
python3 run.py application primary       # reduced model at zero observation variance; participant-level intervals
python3 run.py application profile       # observation-variance profile
python3 run.py application reference     # recovery-reference sensitivity, B = 2000
python3 run.py application diagnostics   # retained supporting checks
python3 run.py application tables
```

or all at once:

```bash
python3 run.py application all --cores 4
```

`--data DIR` reuses already extracted inputs; `--out DIR` writes to a separate run directory. The simulation and recomputed-figure entry points expect the default `results/application` outputs, so use the default location for a complete replication. The `PSS_B` override exists only for development runs and must never be used to claim replication of the published intervals.

Details of script dependencies, seeds, output files and bootstrap settings: `docs/APPLICATION.txt`.

---

## Reproducing the simulations

Run `application preprocess` first: the clustered designs use the actual retained covariates and episode counts. The generating parameters are frozen at full precision in `config/generating_parameters.R`, independently of any new application fit. The `numerical` stage also needs `application prepare`. Independent-episode benchmarks need no CGMacros inputs.

```bash
python3 run.py simulation independent --cores 4
python3 run.py simulation targets --cores 4
python3 run.py simulation clustered --cores 4
python3 run.py simulation summaries
python3 run.py simulation robustness --cores 4
python3 run.py simulation numerical
python3 run.py simulation diagnostics --cores 4
python3 run.py simulation caps --cores 4
```

or everything, including the supplementary numerical and cap-sensitivity checks:

```bash
python3 run.py simulation all --cores 4
```

Defaults preserve the published replicate counts, bootstrap sizes, seeds, critical values, fitted-model definitions and working-likelihood targets. The full bootstrap studies are substantial computations (hundreds of thousands of refits); caches are resumable. Parallel code uses fork-based `parallel::mclapply`, so on Windows use `--cores 1`. See `docs/SIMULATIONS.txt` for prerequisites and the costly stages.

After recomputation, render from the new results:

```bash
python3 run.py figures --source recomputed --include-application
```

This mode never substitutes shipped reference summaries for missing regenerated results.

---

## What is verified, and what is not

Read `docs/VALIDATION.txt` before interpreting a run. It separates the checks completed while packaging (source identity with the working analysis code, 132 synthetic numerical checks, regeneration of the application episodes, bounded refits, exact reproduction of seven final summary tables from saved fits) from the full Monte Carlo experiments, which were **not** rerun during packaging. Machine-readable evidence is in `docs/validation/`.

Two distinctions are preserved throughout the summary scripts:

- under random participant heterogeneity the **conditional generating coefficients** differ from the **working-likelihood target** that the fitted model estimates, and
- the **fixed** and **estimated** observation-variance designs have their own targets and their own generating parameters.

One historical entry (the fixed-variance row of supplementary Table S1) has no traceable original script; reference mode labels it as a published-table snapshot and recomputed mode reports only the traceable examples. See `docs/OUTPUTS.txt`.

---

## Integrity check

```bash
python3 run.py check
```

verifies every packaged source and reference file against `config/release_manifest.json`. The manifest identifies this prepared snapshot; intentional later edits will change the checksums. The GitHub Actions workflow (`.github/workflows/validate.yml`) runs the integrity check and the synthetic tests on every push.

---

## Citation

Please cite the accompanying manuscript when using this code (see `CITATION.cff`), and cite the CGMacros dataset when using the application pipeline:

> Gutierrez-Osuna R, Kerr D, Mortazavi B, Das A (2025). CGMacros: a scientific dataset for personalized nutrition and diet monitoring (version 1.0.0). PhysioNet. https://doi.org/10.13026/3z8q-x658

---

## License

The code in this repository is released under the [MIT License](LICENSE). The CGMacros dataset is not part of this repository and is distributed by PhysioNet under its own terms (CC BY-NC-SA 4.0; see `data/README.txt`).
