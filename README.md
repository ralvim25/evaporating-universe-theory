# The Evaporating Universe I
## A Dissipation Principle — Resolving Cosmic Tensions Without Additional Free Parameters

> **A ghost-condensate interacting dark energy model** — All five parameters derived from the Gribov–Zwanziger effective action in curved spacetime. Zero fitted to data.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20877311.svg)](https://doi.org/10.5281/zenodo.20877311)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

---

## 📊 Key Results

| Result | Value |
|---|---|
| Hubble tension | Resolved from 5.8σ to **0.52σ** |
| Growth tension (S₈) | Resolved from ~2.6σ to **0.58σ** |
| Additional free parameters | **0** (Δk = 0) |
| Model selection (ΔAIC vs ΛCDM) | **−6.88** (strong evidence) |
| DESI DR2 CPL mirage | Predicted w₀ = −0.855, observed −0.838 ± 0.055 (**0.3σ**) |
| Perturbation cancellation theorem | Δσ₈ = 0 (CMB-safe by construction) |
| N-body CDM drain validation | 0.016% agreement with analytic prediction |
| Falsifiable predictions | **10** zero-parameter forecasts (Euclid/LSST/DESI) |

**Companion paper:** [The Evaporating Universe II — The KBC Void and Why Standard Candles Cannot See It](https://github.com/ralvim25/EU_PaperII)

---

## 🔬 Analysis Pipeline

| Step | Notebook | What it does |
|---|---|---|
| 1️⃣ | [`NB01_Derivation_Framework`](notebooks/NB01_Derivation_Framework.ipynb) | Derivation of the five EU parameters from the GZ effective action |
| 2️⃣ | [`NB02_GKI_LKI_Phenomenology`](notebooks/NB02_GKI_LKI_Phenomenology.ipynb) | GKI/LKI phenomenology and Hubble tension resolution |
| 3️⃣ | [`NB03_Background_Cosmology`](notebooks/NB03_Background_Cosmology.ipynb) | Background cosmology: H(z), distances, ages |
| 4️⃣ | [`NB04_Observational_Validation`](notebooks/NB04_Observational_Validation.ipynb) | Observational validation against Planck, BAO, SNe Ia |
| 5️⃣ | [`NB05_MCMC`](notebooks/NB05_Nomenclatura_Runs.md) (×6) | MCMC analyses: configurations A1, A2, C1, C2, D1, D2 |
| 6️⃣ | [`NB06_Model_Selection`](notebooks/NB06_Model_Selection.ipynb) | Model selection: ΔAIC, Δχ², information criteria |
| 7️⃣ | [`NB07_Nbody_Validation`](notebooks/NB07_Nbody_Validation.ipynb) | Gadget-4 N-body simulation validation |
| 8️⃣ | [`NB08_S8_Forensics`](notebooks/NB08_S8_Forensics.ipynb) | S₈ forensic analysis: kernel bias and growth tension |
| 9️⃣ | [`NB09_Euclid_Predictions`](notebooks/NB09_Euclid_Predictions.ipynb) | 10 zero-parameter forecasts (Euclid/LSST/DESI) |

---

## 📁 Repository Structure

```
evaporating-universe-theory/
├── README.md
├── REPRODUCING.md                     ← Detailed reproduction guide
├── LICENSE
├── class-eu/                          ← CLASS-EU Boltzmann solver outputs
│   ├── eu_modeB_background.dat        ← EU Mode B expansion history E(z)
│   ├── eu_modeB_cl.dat                ← EU Mode B power spectra
│   ├── eu_modeB_pk.dat                ← EU Mode B matter P(k)
│   ├── lcdm_background.dat            ← ΛCDM reference background
│   ├── lcdm_cl.dat                    ← ΛCDM reference spectra
│   └── lcdm_pk.dat                    ← ΛCDM reference P(k)
├── notebooks/                         ← Analysis pipeline (NB01–NB09)
│   ├── NB01–NB04 (.ipynb)             ← Theory + background + validation
│   ├── NB05_MCMC_A1/A2/C1/C2/D1/D2   ← Six MCMC configurations
│   ├── NB06–NB09 (.ipynb)             ← Model selection, N-body, S₈, predictions
│   ├── mcmc-chains/                   ← MCMC chain files
│   └── results/                       ← JSON outputs from each notebook
├── nbody/                             ← Gadget-4 N-body simulation
│   ├── input/                         ← Initial conditions
│   └── simulation/                    ← Power spectrum snapshots
└── paper/
    ├── eu_paper1.tex                  ← Manuscript source (REVTeX 4.2)
    ├── eu_paper1.bib                  ← Bibliography
    ├── eu_paper1.pdf                  ← Compiled manuscript
    └── figures/                       ← Publication figures (PDF)
```

---

## ⚙️ Computational Infrastructure

| Component | Tool | Description |
|---|---|---|
| Boltzmann solver | CLASS-EU (modified) | CDM → vacuum drain implementation |
| MCMC sampler | Cobaya + Metropolis-Hastings | 6 run configurations (A1–D2) |
| Weak lensing | CosmoLike | DES-Y3 3×2pt likelihood (462 data points) |
| N-body simulation | Gadget-4 (patched) | 1024³ particles, CDM drain at z = 0 |
| Initial conditions | monofonIC | FFTW-MPI parallel IC generation |

---

## ▶️ How to Reproduce

**Option 1 — Browse on GitHub (no install)**
> Click any notebook above — GitHub renders `.ipynb` files with all code, outputs, and figures.

**Option 2 — Run on Google Colab**

Open any notebook directly in Colab:

| Notebook | Colab |
|---|---|
| NB01 — Derivation | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ralvim25/evaporating-universe-theory/blob/main/notebooks/NB01_Derivation_Framework.ipynb) |
| NB02 — GKI/LKI | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ralvim25/evaporating-universe-theory/blob/main/notebooks/NB02_GKI_LKI_Phenomenology.ipynb) |
| NB03 — Background | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ralvim25/evaporating-universe-theory/blob/main/notebooks/NB03_Background_Cosmology.ipynb) |
| NB04 — Validation | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ralvim25/evaporating-universe-theory/blob/main/notebooks/NB04_Observational_Validation.ipynb) |
| NB06 — Model Selection | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ralvim25/evaporating-universe-theory/blob/main/notebooks/NB06_Model_Selection.ipynb) |
| NB07 — N-body | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ralvim25/evaporating-universe-theory/blob/main/notebooks/NB07_Nbody_Validation.ipynb) |
| NB08 — S₈ Forensics | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ralvim25/evaporating-universe-theory/blob/main/notebooks/NB08_S8_Forensics.ipynb) |
| NB09 — Predictions | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/ralvim25/evaporating-universe-theory/blob/main/notebooks/NB09_Euclid_Predictions.ipynb) |

> **Note:** For MCMC notebooks (NB05_*), see [`REPRODUCING.md`](REPRODUCING.md) for full Cobaya/CosmoLike setup instructions.

---

## 📄 Citation

If you use this code or data, please cite:

```bibtex
@article{alvim2026_eu,
  author  = {Alvim, Ricardo},
  title   = {The Evaporating Universe {I}: A dissipation principle ---
             Resolving cosmic tensions without additional free parameters},
  year    = {2026},
  doi     = {10.5281/zenodo.20877311},
  note    = {Submitted to Physics of the Dark Universe}
}
```

---

## 👤 Author

**Ricardo Alvim** — Independent Researcher, Rio de Janeiro, Brazil
- ✉️ ricardoalv1m@pm.me
- 🔗 [ORCID: 0009-0004-7744-175X](https://orcid.org/0009-0004-7744-175X)

---

## 📜 License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).
