"""
Figure 5 — Kernel Bias Gradient (Drain Freezing)
=================================================
Paper I: The Evaporating Universe
Section: §X.B (Lensing kernel bias and the CDM drain freezing)

ΔW/W kernel bias (%) vs mean redshift z for Euclid DR1 (10 bins)
and LSST Y1 (5 bins). The monotonic decrease in bias magnitude
with z is a direct imprint of the CDM drain freezing at high z.

Style: PRD single-column (3.375 inches wide), Computer Modern fonts
RULE:  ZERO text overlap with curves, data points, or other text.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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
    'xtick.labelsize': 10.5,
    'ytick.labelsize': 10.5,
    'legend.fontsize': 10,
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

with open(results_dir / 'NB09_Euclid_results.json') as f:
    nb09 = json.load(f)

kb = nb09['ip1_kernel_bias']

# Euclid DR1 (10 bins)
euclid = kb['Euclid DR1']
z_euclid = np.array([b['z_mean'] for b in euclid])
bias_euclid = np.array([b['kernel_bias_pct'] for b in euclid])

# LSST Y1 (5 bins)
lsst = kb['LSST (Rubin)']
z_lsst = np.array([b['z_mean'] for b in lsst])
bias_lsst = np.array([b['kernel_bias_pct'] for b in lsst])

# ================================================================
# SMOOTH ANALYTIC CURVE (PCHIP through combined data)
# ================================================================
from scipy.interpolate import PchipInterpolator

# Combine both surveys for the trend curve
z_all = np.concatenate([z_euclid, z_lsst])
b_all = np.concatenate([bias_euclid, bias_lsst])
sort_idx = np.argsort(z_all)
z_sorted = z_all[sort_idx]
b_sorted = b_all[sort_idx]

# Remove duplicates (average if same z)
z_unique, inv = np.unique(z_sorted, return_inverse=True)
b_unique = np.array([b_sorted[inv == i].mean() for i in range(len(z_unique))])

interp = PchipInterpolator(z_unique, b_unique)
z_fine = np.linspace(0.15, 2.0, 200)
bias_fine = interp(z_fine)

# ================================================================
# COLOR PALETTE
# ================================================================
col_euclid = '#C62828'   # dark red (Euclid)
col_lsst   = '#1565C0'   # blue (LSST)
col_curve  = '#424242'   # dark gray (trend line)

# ================================================================
# FIGURE
# ================================================================
fig, ax = plt.subplots(figsize=(3.375, 2.8))

# --- Smooth trend curve ---
ax.plot(z_fine, bias_fine, color=col_curve, lw=1.0, ls='--',
        alpha=0.5, zorder=2, label='_nolegend_')

# --- Euclid DR1 markers (circles) ---
ax.plot(z_euclid, bias_euclid, 'o', color=col_euclid, ms=5,
        markeredgecolor='white', markeredgewidth=0.6, zorder=5,
        label=r'Euclid DR1 (10 bins)')

# --- LSST Y1 markers (hollow squares) ---
ax.plot(z_lsst, bias_lsst, 's', color=col_lsst, ms=6,
        markerfacecolor='none', markeredgecolor=col_lsst,
        markeredgewidth=1.2, zorder=6,
        label=r'LSST Y1 (5 bins)')

# ================================================================
# ANNOTATIONS — following the curve, no overlap
# ================================================================
# "Active CDM Drain" at low-z (left side)
ax.annotate(r'Active CDM Drain',
            xy=(0.35, -7.94), xytext=(0.55, -8.3),
            fontsize=10, color=col_curve, fontweight='bold',
            ha='left', va='center',
            arrowprops=dict(arrowstyle='->', color=col_curve,
                            lw=0.6, shrinkB=5, connectionstyle='arc3,rad=-0.2'),
            bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                      edgecolor='none', alpha=0.85))

# "Drain Freezing" at high-z (right side)
ax.annotate(r'Drain Freezing',
            xy=(1.75, -4.81), xytext=(1.45, -4.4),
            fontsize=10, color=col_curve, fontweight='bold',
            ha='right', va='center',
            arrowprops=dict(arrowstyle='->', color=col_curve,
                            lw=0.6, shrinkB=5, connectionstyle='arc3,rad=0.2'),
            bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                      edgecolor='none', alpha=0.85))

# Gradient magnitude annotation
ax.text(0.85, -5.0,
        r'$40\%$ decrease',
        fontsize=10, color=col_curve, ha='center', va='center',
        alpha=0.6, style='italic')

# ================================================================
# LEGEND
# ================================================================
ax.legend(loc='center right', frameon=True, framealpha=0.85,
          edgecolor='none', fontsize=8.5, handletextpad=0.4,
          borderpad=0.3, bbox_to_anchor=(1.0, 0.38))

# ================================================================
# AXIS CONFIGURATION
# ================================================================
ax.set_xlabel(r'Mean redshift $\bar{z}$')
ax.set_ylabel(r'Kernel bias $\Delta W_i / W_i$ [%]', labelpad=-1)
ax.set_xlim(0.15, 2.05)
ax.set_ylim(-8.8, -4.0)

ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

# ================================================================
# SAVE
# ================================================================
output_dir = Path(r'c:\Users\ricar\Desktop\PAPER 1\Current'
                  r'\Paper_I_A_Dissipation_Principle\04_Figures')
fig.savefig(output_dir / 'fig5_kernel_gradient.pdf', format='pdf')
fig.savefig(output_dir / 'fig5_kernel_gradient.png', format='png')
plt.show()
print("Figure 5 saved: fig5_kernel_gradient.pdf + fig5_kernel_gradient.png")
