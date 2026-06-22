# Reproducing the Results

This guide explains how to independently verify all numerical results
presented in the paper using Google Colab (free, no local installation
required).

---

## Quick Start

### Step 1 — Open Google Colab

Go to [colab.research.google.com](https://colab.research.google.com/)

### Step 2 — Upload a notebook

Click **File → Upload notebook** and select any notebook from the
`notebooks/` folder. All dependencies (`numpy`, `scipy`,
`matplotlib`) are pre-installed in Colab.

### Step 3 — Upload the required data

Each notebook loads results from previous notebooks via JSON files
in `notebooks/results/`. Upload the corresponding JSON files to
your Colab session when prompted, or upload the entire `results/`
folder.

### Important: Run in chronological order

The notebooks form a sequential pipeline — each one depends on
results produced by the previous ones:

```
NB01 → NB02 → NB03 → NB04 → NB05 (×6) → NB06 → NB07 → NB08 → NB09
```

If you only want to verify a specific result, the pre-computed
JSON files in `notebooks/results/` allow you to start from any
notebook without re-running earlier ones.

---

## Notebook Verification Guide

Each notebook produces results that can be cross-checked against
the paper. The expected outputs are listed below.

### NB01 — Derivation Framework
**Verifies:** The five EU parameters (λ, b, z_trans, ε_bare, ε_IR)
**Paper sections:** II, III
**Expected outputs:**
- λ = 2/3 (exact)
- b = 19/36 (exact)
- z_trans ≈ 5.986
- ε_bare ≈ 0.0543
- ε_IR ≈ 0.0426

### NB02 — GKI/LKI Phenomenology
**Verifies:** Hubble tension resolution
**Paper sections:** IV.A, IV.B
**Expected outputs:**
- H₀_GKI = 68.89 ± 0.28 km/s/Mpc
- H₀_LKI = 72.71 ± 0.30 km/s/Mpc
- SH0ES tension reduced to 0.52σ

### NB03 — Background Cosmology
**Verifies:** Modified Friedmann equations, H(z), distances
**Paper sections:** II.C, IV
**Key check:** EU and ΛCDM background evolutions diverge below
z_trans ≈ 5.986

### NB04 — Observational Validation
**Verifies:** Consistency with Planck, BAO, SNe Ia
**Paper sections:** IV, V
**Key check:** EU predictions within 1σ of all datasets

### NB05 (×6) — MCMC Analyses
**Verifies:** Bayesian parameter estimation
**Paper section:** V
**Expected outputs:**
- All runs converge at R⁻¹ < 0.01
- UV-derived parameters match posteriors to 0.04σ
- Six configurations: A1, A2, C1, C2, D1, D2

> **Note:** MCMC chains are not included in this repository due to
> file size (~2 GB per run). The notebooks load pre-computed summary
> statistics from `notebooks/results/NB05_*.json`. Full chains are
> available upon request.

### NB06 — Model Selection
**Verifies:** ΔAIC, Δχ², information criteria
**Paper section:** V.E
**Expected outputs:**
- ΔAIC = −6.88 (fair-fight configuration)

### NB07 — N-body Validation
**Verifies:** CDM drain signature in P(k)
**Paper section:** VI
**Expected outputs:**
- Drain accuracy: 0.003% at z = 0
- CDM loss: 3.64% at z = 0

> **Note:** Simulation snapshots are not included due to file size.
> The notebook loads pre-computed power spectra from
> `nbody/input/pk_*.txt`. Full Gadget-4 outputs are archived on
> AWS S3 and available upon request.

### NB08 — S₈ Forensics
**Verifies:** Lensing kernel bias and growth tension resolution
**Paper sections:** VIII, IX
**Expected outputs:**
- S₈(EU) = 0.811
- Ŝ₈(DES-Y3, kernel-corrected) = 0.794
- Tension reduced from ~2.6σ to 1.03σ

### NB09 — Euclid Predictions
**Verifies:** Ten zero-parameter predictions
**Paper section:** X
**Expected outputs:**
- P1–P5: Euclid DR1 C_ℓ suppression (−9% to −16%)
- P6: CPL mirage w₀ = −0.855
- P7–P10: LSST/DESI predictions

---

## Directory Reference

```
notebooks/results/    ← Pre-computed JSON outputs from each notebook
class-eu/             ← CLASS-EU solver outputs (background, Cℓ, P(k))
nbody/input/          ← Gadget-4 power spectrum snapshots (15 redshifts)
paper/figures/        ← Publication figures (PDF) + generation scripts (Python)
```

---

## About the Computational Pipeline

The notebooks in this repository **validate and analyze**
pre-computed results. The heavy computations (MCMC chains,
N-body simulation) were performed on cloud infrastructure
and their outputs are provided as JSON summaries and power
spectrum snapshots.

| Component | Status | Data provided |
|:--|:--|:--|
| Boltzmann solver (CLASS-EU) | ✅ Included | `class-eu/*.dat` |
| MCMC chains (6 runs, Cobaya) | ✅ Included | `notebooks/mcmc-chains/RUN_*/` |
| MCMC summaries | ✅ Included | `notebooks/results/NB05_*.json` |
| N-body simulation (Gadget-4) | ✅ Included | `nbody/input/pk_*.txt` |

Full Gadget-4 simulation snapshots (HDF5) are available upon
request.

The modified CLASS source code (CLASS-EU patches), Cobaya YAML
configuration files, and Gadget-4 CDM drain patches will be made
available in this repository following publication.

---

## Questions?

For questions about reproduction or to request full data,
contact:

**Ricardo Alvim** — ricardoalv1m@pm.me

