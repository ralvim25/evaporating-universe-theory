# The Evaporating Universe

**Paper I: A Dissipation Principle — Resolving Cosmic Tensions Without Additional Free Parameters**

[![Zenodo](https://img.shields.io/badge/Zenodo-10.5281%2Fzenodo.20790999-blue)](https://doi.org/10.5281/zenodo.20790999)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

---

## Overview

The Evaporating Universe (EU) is an interacting dark energy model in which a Gribov–Zwanziger ghost condensate mediates a continuous transfer of energy from cold dark matter (CDM) to the vacuum. All five model parameters are derived from effective field theory in curved spacetime; none are fitted to data (Δk = 0).

### Key Results

| Result | Value |
|:--|:--|
| Hubble tension | Resolved from 5.8σ to **0.52σ** |
| Growth tension (S₈) | Resolved from ~2.6σ to **0.58σ** |
| Additional free parameters | **0** |
| Model selection (ΔAIC) | **−6.88** (strong evidence) |
| DESI DR2 CPL mirage | Predicted w₀ = −0.855, observed −0.838 ± 0.055 (**0.3σ**) |

---

## Repository Structure

```
evaporating-universe-theory/
├── paper/                     # Manuscript source
│   ├── eu_paper1.tex          # LaTeX source (REVTeX 4.2)
│   ├── eu_paper1.bib          # Bibliography
│   └── figures/               # Publication-quality figures (PDF)
├── notebooks/                 # Analysis notebooks (Jupyter)
│   ├── NB01_Derivation_Framework.ipynb
│   ├── NB02_GKI_LKI_Phenomenology.ipynb
│   ├── NB03_Background_Cosmology.ipynb
│   ├── NB04_Observational_Validation.ipynb
│   ├── NB05_MCMC_*.ipynb      # Six MCMC run analyses
│   ├── NB06_Model_Selection.ipynb
│   ├── NB07_Nbody_Validation.ipynb
│   ├── NB08_S8_Forensics.ipynb
│   ├── NB09_Euclid_Predictions.ipynb
│   └── results/               # JSON outputs from each notebook
├── class-eu/                  # CLASS-EU solver outputs
├── nbody/                     # Gadget-4 N-body simulation data
│   └── input/                 # Power spectrum snapshots
└── figures/
    └── scripts/               # Python scripts for figure generation
```

---

## Notebooks

| Notebook | Description |
|:--|:--|
| **NB01** | Derivation of the five EU parameters from the GZ effective action |
| **NB02** | GKI/LKI phenomenology and Hubble tension resolution |
| **NB03** | Background cosmology: H(z), distances, ages |
| **NB04** | Observational validation against Planck, BAO, SNe Ia |
| **NB05** (×6) | MCMC analyses: configurations A1, A2, C1, C2, D1, D2 |
| **NB06** | Model selection: ΔAIC, Δχ², information criteria |
| **NB07** | Gadget-4 N-body simulation validation |
| **NB08** | S₈ forensic analysis: kernel bias and growth tension |
| **NB09** | Euclid/LSST/DESI predictions (10 zero-parameter forecasts) |

---

## Computational Infrastructure

- **Boltzmann solver:** Modified CLASS (CLASS-EU) with CDM drain implementation
- **MCMC sampler:** Cobaya with Metropolis-Hastings
- **Weak lensing:** CosmoLike (DES-Y3 3×2pt likelihood)
- **N-body:** Gadget-4 with CDM drain patches (1024³ particles)
- **Initial conditions:** monofonIC

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{alvim2026_eu,
  author  = {Alvim, Ricardo},
  title   = {The Evaporating Universe~{I}: A dissipation principle ---
             Resolving cosmic tensions without additional free parameters},
  year    = {2026},
  doi     = {10.5281/zenodo.20790999},
  note    = {Manuscript submitted for publication}
}
```

---

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

---

## Contact

Ricardo Alvim — [ricardoalv1m@pm.me](mailto:ricardoalv1m@pm.me)  
ORCID: [0009-0004-7744-175X](https://orcid.org/0009-0004-7744-175X)
