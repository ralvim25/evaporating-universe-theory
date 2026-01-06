# The Evaporating Universe Theory

**Dark Sector Unification via Dissipative Soliton Phase Transition**

This repository contains the source code and numerical simulation scripts used to generate the figures for the paper *The Evaporating Universe Theory* by Ricardo Alvim.

## 🌌 Overview

The code simulates a Unified Dark Fluid (UDF) that behaves as:
1.  **Superfluid Solitons (Dark Matter)** in high-density regimes.
2.  **Dissipative Dark Radiation (Dark Energy)** in low-density regimes.

It addresses the Hubble Tension, Lithium Problem, and $S_8$ Tension using a modified Boltzmann solver approach.

## 🚀 Installation

### Prerequisites
You will need Python 3.8+ and the following libraries:

```bash
pip install -r requirements.txt
```

Note on CLASS (Boltzmann Solver)
Figures 4 and 5 require the classy Python wrapper for the CLASS Boltzmann solver.

If classy is installed, all figures will be generated.

If classy is not installed, the script will skip Figures 4 and 5 but still generate Figures 1, 2, 3, and 6.

To install CLASS, please refer to the [official CLASS repository](https://github.com/lesgourg/class_public).

📊 Usage
Simply run the main script to generate the plots:

```bash
python evaporating_universe_alvim_theory.py
```
The output images will be saved in the root directory as PNG files (300 DPI).

📜 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

✍️ Citation
If you use this code or model in your research, please cite:

Alvim, R. (2026). The Evaporating Universe Theory: Dark Sector Unification via Dissipative Soliton Phase Transition. [Journal Name/Preprint Link].


