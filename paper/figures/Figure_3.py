"""
Figure 3 — Non-Linear Fingerprint (Ratio Plot)
=================================================
Paper I: The Evaporating Universe
Section: §VIII.C (Power spectrum suppression and anemic halos)

S(k) = P_EU(k) / P_LCDM(k) showing 3 regimes:
  I.   Linear Enhancement  (k < 0.1):  ~1.065
  II.  Anemic Halos         (0.3 < k < 2): dip relative to peak
  III. Baryonic Artifact    (k > 3): DMO divergence

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

with open(results_dir / 'NB07_results.json') as f:
    nb07 = json.load(f)

sr = nb07['suppression_ratio']
CUB = sr['CUB_flat_suppression']  # 1.0383

# Load multi-z S(k) data
redshifts = ['z=0.00', 'z=0.50', 'z=1.00', 'z=2.00']
z_labels  = [r'$z = 0$', r'$z = 0.5$', r'$z = 1$', r'$z = 2$']

data = {}
for zkey in redshifts:
    entry = sr['per_redshift'][zkey]
    k = np.array(entry['k_full'])
    S = np.array(entry['S_full'])
    data[zkey] = (k, S)

# ================================================================
# COLOR PALETTE
# ================================================================
col_z0    = '#C62828'   # z=0 — primary
col_z05   = '#E65100'   # z=0.5
col_z1    = '#2E7D32'   # z=1
col_z2    = '#1565C0'   # z=2
z_colors  = [col_z0, col_z05, col_z1, col_z2]
z_lstyles = ['-', '--', '-.', ':']
col_lcdm  = '#546E7A'
col_wall  = '#B71C1C'

# ================================================================
# FIGURE
# ================================================================
fig, ax = plt.subplots(figsize=(3.375, 2.6))

# --- Smoothing function for noisy N-body data ---
def smooth(y, window=11):
    """Simple moving average to smooth N-body noise."""
    w = np.ones(window) / window
    padded = np.concatenate([y[:window], y, y[-window:]])
    smoothed = np.convolve(padded, w, mode='same')
    return smoothed[window:-window]

# --- Baryonic Artifact zone (hatched, behind everything) ---
ax.fill_betweenx([0.96, 1.25], 2.5, 10,
                 color='#FFCDD2', alpha=0.25, zorder=0,
                 hatch='///', edgecolor='#EF9A9A', linewidth=0.3)

# --- Anemic Halos zone (subtle shading for the quasi-linear regime) ---
ax.fill_betweenx([0.96, 1.25], 0.3, 2.5,
                 color='#E1BEE7', alpha=0.10, zorder=0)

# --- Plot S(k) curves for each redshift ---
for i, zkey in enumerate(redshifts):
    k, S = data[zkey]
    S_smooth = smooth(S, window=11)
    mask = k < 3.5  # Extend into DMO zone to show divergence
    lw = 1.1 if i == 0 else 0.8
    alpha = 1.0 if i == 0 else 0.6
    ax.plot(k[mask], S_smooth[mask],
            color=z_colors[i], ls=z_lstyles[i], lw=lw, alpha=alpha,
            label=z_labels[i], zorder=4 - i)

# --- CUB line ---
ax.axhline(y=CUB, color=col_lcdm, ls='--', lw=0.6, alpha=0.4, zorder=2)

# --- S = 1 reference ---
ax.axhline(y=1.0, color='black', ls='-', lw=0.4, alpha=0.25, zorder=1)

# ================================================================
# ANNOTATIONS — zero overlap, positioned in clean zones
# ================================================================

# CUB label — left side, below the dashed line
ax.text(0.035, CUB - 0.004, r'CUB $= 1.038$',
        fontsize=10, color=col_lcdm, ha='left', va='top')

# I. Linear Enhancement — bottom-left clean zone, below all curves
ax.text(0.12, 0.99, 'Linear Enhancement',
        fontsize=10, color=col_z0, ha='center', va='center',
        style='italic',
        bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                  edgecolor='none', alpha=0.85))

# II. Anemic Halos — inside the shaded regime zone
ax.text(0.85, 1.17, 'Anemic Halos',
        fontsize=10, color='#4A148C', ha='center', va='center',
        style='italic',
        bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                  edgecolor='#CE93D8', alpha=0.85, lw=0.4))

# III. DMO Artifact — at CUB height, inside hatched zone
ax.text(2.95, CUB, 'DMO\nArtifact',
        fontsize=10, color=col_wall, ha='center', va='center',
        style='italic', alpha=0.9,
        bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                  edgecolor='none', alpha=0.85))

# ================================================================
# AXIS CONFIGURATION
# ================================================================
ax.set_xscale('log')
ax.set_xlabel(r'$k\;\;[h\,\mathrm{Mpc}^{-1}]$')
ax.set_ylabel(r'$\mathcal{S}(k) = P_{\rm EU} / P_{\Lambda\rm CDM}$',
              labelpad=-1)
ax.set_xlim(0.03, 3.5)
ax.set_ylim(0.97, 1.20)

ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

# Legend — compact, top-left clean zone
ax.legend(loc='upper left', frameon=True, framealpha=0.95,
          edgecolor='#BDBDBD', fancybox=False, borderpad=0.25,
          handletextpad=0.3, handlelength=1.5, fontsize=10,
          bbox_to_anchor=(0.0, 1.0))

# ================================================================
# SAVE
# ================================================================
output_dir = Path(r'c:\Users\ricar\Desktop\PAPER 1\Current'
                  r'\Paper_I_A_Dissipation_Principle\04_Figures')
fig.savefig(output_dir / 'fig3_power_spectrum.pdf', format='pdf')
fig.savefig(output_dir / 'fig3_power_spectrum.png', format='png')
plt.show()
print("Figure 3 saved: fig3_power_spectrum.pdf + fig3_power_spectrum.png")
