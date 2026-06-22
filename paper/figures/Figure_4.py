"""
Figure 4 — Centro de Gravidade S₈ (Forest Plot)
=================================================
Paper I: The Evaporating Universe
Section: §IX.B (Quantitative resolution)

Horizontal forest plot showing 4 lensing surveys:
  - Published S₈ (blue squares + error bars)
  - EU kernel-corrected Ŝ₈ (red diamonds)
  - Arrows showing the correction direction
  - Planck ΛCDM band (grey) vs EU prediction band (green)

Style: PRD single-column (3.375 inches wide), Computer Modern fonts
RULE:  ZERO text overlap with data points, error bars, or bands.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from pathlib import Path

# ================================================================
# PRD STYLE SETUP
# ================================================================
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['CMU Serif', 'Computer Modern Roman', 'Times New Roman'],
    'mathtext.fontset': 'cm',
    'font.size': 8,
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 7.5,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.width': 0.35,
    'ytick.minor.width': 0.35,
    'xtick.major.size': 3.5,
    'ytick.major.size': 3.5,
    'xtick.minor.size': 2.0,
    'ytick.minor.size': 2.0,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'lines.linewidth': 1.0,
    'figure.dpi': 300,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.03,
})

# ================================================================
# LOAD DATA
# ================================================================
results_dir = Path(r'c:\Users\ricar\Desktop\PAPER 1\Current'
                   r'\Paper_I_A_Dissipation_Principle\01_Notebooks\results')

with open(results_dir / 'NB08_S8_results.json') as f:
    nb08 = json.load(f)

s8_data = nb08['gap1_s8_inference']

# Reference values
S8_Planck     = 0.832
S8_Planck_err = 0.013
S8_EU_true    = 0.8112
S8_EU_err     = 0.008  # from spec

# Survey order (bottom to top)
surveys = ['KiDS\n1000', 'HSC-Y3', 'DES-Y3', 'KiDS\nLegacy']
survey_keys = ['KiDS-1000', 'HSC-Y3', 'DES-Y3', 'KiDS-Legacy']
y_pos = [0, 1, 2, 3]

# Extract per-survey values
S8_pub  = [s8_data[s]['S8_published'] for s in survey_keys]
S8_err  = [s8_data[s]['S8_err'] for s in survey_keys]
S8_corr = [s8_data[s]['S8_LCDM_inferred'] for s in survey_keys]
tensions = [s8_data[s]['tension_sigma'] for s in survey_keys]

# ================================================================
# COLOR PALETTE
# ================================================================
col_planck = '#546E7A'   # grey-blue
col_eu     = '#2E7D32'   # green
col_pub    = '#1565C0'   # blue (published)
col_corr   = '#C62828'   # red (corrected)
col_arrow  = '#757575'   # grey arrows

# ================================================================
# FIGURE
# ================================================================
fig, ax = plt.subplots(figsize=(3.375, 2.4))

# --- Reference bands ---
# Planck ΛCDM band
ax.axvspan(S8_Planck - S8_Planck_err, S8_Planck + S8_Planck_err,
           color=col_planck, alpha=0.10, zorder=0)
ax.axvline(x=S8_Planck, color=col_planck, ls='-', lw=0.6, alpha=0.4, zorder=1)

# EU True Prediction band
ax.axvspan(S8_EU_true - S8_EU_err, S8_EU_true + S8_EU_err,
           color=col_eu, alpha=0.12, zorder=0)
ax.axvline(x=S8_EU_true, color=col_eu, ls='-', lw=0.6, alpha=0.4, zorder=1)

# --- Data points ---
for i, survey in enumerate(surveys):
    y = y_pos[i]

    # Published S₈ (blue square + error bar)
    ax.errorbar(S8_pub[i], y, xerr=S8_err[i],
                fmt='s', color=col_pub, ms=4, lw=0.8,
                capsize=2.5, capthick=0.6, zorder=5,
                markeredgewidth=0.5)

    # EU kernel-corrected Ŝ₈ (red diamond)
    ax.plot(S8_corr[i], y, 'D', color=col_corr, ms=4,
            markeredgewidth=0.5, zorder=5)

    # Thin connector line: published → corrected
    ax.plot([S8_pub[i], S8_corr[i]], [y, y],
            color=col_arrow, lw=0.5, alpha=0.4, zorder=3)

# --- Tension annotations (right side, clear of data) ---
x_sigma = S8_corr[2] + 0.005  # align all on same x (DES-Y3 ref)

# DES-Y3: tension with EU
ax.text(x_sigma, y_pos[2] + 0.15,
        r'$1.03\sigma$',
        fontsize=8, color=col_corr, ha='left', va='bottom')

# KiDS-Legacy: tension with EU
ax.text(x_sigma, y_pos[3] + 0.15,
        r'$1.15\sigma$',
        fontsize=8, color=col_corr, ha='left', va='bottom')

# KiDS-1000
ax.text(x_sigma, y_pos[0] + 0.15,
        r'$1.38\sigma$',
        fontsize=8, color=col_corr, ha='left', va='bottom')

# HSC-Y3
ax.text(x_sigma, y_pos[1] + 0.15,
        r'$0.63\sigma$',
        fontsize=8, color=col_corr, ha='left', va='bottom')

# --- Band labels (top of plot, clear zone) ---
ax.text(S8_Planck, 4.3,
        r'Planck $\Lambda$CDM' + '\n' + r'$0.832 \pm 0.013$',
        fontsize=8, color=col_planck, ha='center', va='center',
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                  edgecolor='none', alpha=0.85))

ax.text(S8_EU_true, -0.7,
        r'EU prediction' + '\n' + r'$0.811 \pm 0.008$',
        fontsize=8, color=col_eu, ha='center', va='center',
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                  edgecolor='none', alpha=0.85))

# ================================================================
# AXIS CONFIGURATION
# ================================================================
ax.set_yticks(y_pos)
ax.set_yticklabels(surveys, fontsize=8.5)
ax.set_xlabel(r'$S_8 = \sigma_8 \sqrt{\Omega_m / 0.3}$')
ax.set_xlim(0.73, 0.86)
ax.set_ylim(-1.2, 5.0)

ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', length=0)

ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))

# --- Legend ---
legend_handles = [
    Line2D([0], [0], color=col_pub, marker='s', ms=4, lw=0,
           markeredgewidth=0.5,
           label=r'$S_8^{\rm pub}$ (published)'),
    Line2D([0], [0], color=col_corr, marker='D', ms=4, lw=0,
           markeredgewidth=0.5,
           label=r'$\hat{S}_8^{\rm EU}$ (kernel-corrected)'),
]
ax.legend(handles=legend_handles, loc='upper left',
          frameon=True, framealpha=0.95,
          edgecolor='#BDBDBD', fancybox=False, borderpad=0.25,
          handletextpad=0.3, fontsize=7.5,
          bbox_to_anchor=(-0.02, 1.02))

# ================================================================
# SAVE
# ================================================================
output_dir = Path(r'c:\Users\ricar\Desktop\PAPER 1\Current'
                  r'\Paper_I_A_Dissipation_Principle\04_Figures')
fig.savefig(output_dir / 'fig4_s8_forest.pdf', format='pdf')
fig.savefig(output_dir / 'fig4_s8_forest.png', format='png')
plt.show()
print("Figure 4 saved: fig4_s8_forest.pdf + fig4_s8_forest.png")
