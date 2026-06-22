"""
Figure 6 — A Miragem do DESI (CPL Mirage)
===========================================
Paper I: The Evaporating Universe
Section: §X.C (Distance probes and the equation-of-state mirage)

w₀ apparent as a function of forced Ωm prior, showing that the
DESI "phantom dark energy" signal is a prior artifact.

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

mirage = nb09['cpl_mirage']
scan = mirage['Om_scan']

# Extract scan arrays
Om_scan = np.array([p['Om_fixed'] for p in scan])
w0_scan = np.array([p['w0'] for p in scan])

# Key points (matching §X.C text)
Om_postdrain = 0.288   # EU post-drain value
w0_postdrain = -1.045  # from Om_scan
Om_FitB      = 0.3153  # Fit B: Ωm fixed to Planck CMB prior
w0_FitB      = mirage['fit_B_Om_fixed']['w0_pipe']  # -0.8549

# Build PCHIP interpolator for smooth curve + crossover search
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
interp = PchipInterpolator(Om_scan, w0_scan)

# Find crossover Ωm where w₀ = −1
Om_cross = brentq(lambda om: interp(om) - (-1.0), 0.288, 0.315)
w0_cross = -1.0  # by definition

# DESI DR2 measurement
w0_DESI     = -0.838
w0_DESI_err = 0.055

# ================================================================
# COLOR PALETTE
# ================================================================
col_curve  = '#C62828'   # EU red
col_eu_pt  = '#2E7D32'   # green (EU true)
col_cmb_pt = '#C62828'   # red (CMB artifact)
col_desi   = '#1565C0'   # blue (DESI)
col_cc     = '#212121'   # black (cosmological constant)

# ================================================================
# FIGURE
# ================================================================
fig, ax = plt.subplots(figsize=(3.375, 2.8))

# --- Interpolated smooth curve ---
Om_fine = np.linspace(0.265, 0.335, 200)
w0_fine = interp(Om_fine)

ax.plot(Om_fine, w0_fine, color=col_curve, lw=1.3, zorder=4,
        label=r'EU distance fit')

# --- w₀ = −1 reference line ---
ax.axhline(y=-1.0, color=col_cc, ls='--', lw=0.8, alpha=0.6, zorder=2)
ax.text(0.268, -0.99, r'$w_0 = -1$',
        fontsize=9, color=col_cc, ha='left', va='bottom', alpha=0.7)

# --- Ωm = 0.315 vertical line removed (redundant with CMB Prior label) ---

# --- DESI DR2 band ---
ax.axhspan(w0_DESI - w0_DESI_err, w0_DESI + w0_DESI_err,
           color=col_desi, alpha=0.08, zorder=0)
ax.axhline(y=w0_DESI, color=col_desi, ls='-', lw=0.5, alpha=0.3, zorder=1)

# --- Key data points ---
# EU post-drain (green circle)
ax.plot(Om_postdrain, w0_postdrain, 'o', color=col_eu_pt, ms=8,
        markeredgecolor='white', markeredgewidth=0.8, zorder=6)

# Crossover (green diamond — on the w₀=-1 line)
ax.plot(Om_cross, w0_cross, 'D', color=col_eu_pt, ms=7,
        markeredgecolor='white', markeredgewidth=0.8, zorder=6)

# Fit B: Planck CMB prior (red circle)
ax.plot(Om_FitB, w0_FitB, 'o', color=col_cmb_pt, ms=8,
        markeredgecolor='white', markeredgewidth=0.8, zorder=6)

# --- Point labels — following the curve trajectory, offset right ---
# Compute x-position of curve at each label's y-height
from scipy.optimize import brentq as _brentq

def curve_x_at(w0_target):
    """Find Ωm where the curve equals w0_target."""
    return _brentq(lambda om: interp(om) - w0_target, 0.265, 0.335)

label_offset = 0.007  # offset to the right of the curve

# DESI DR2 (blue, at curve height)
x_desi = curve_x_at(-0.81)
ax.text(0.268, -0.79,
        r'DESI DR2' + '\n' + r'$w_0 = -0.838 \pm 0.055$',
        fontsize=10, color=col_desi, ha='left', va='center',
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                  edgecolor='none', alpha=0.55))

# CMB Prior (red, at curve height)
x_cmb = curve_x_at(w0_FitB)
ax.text(x_cmb + label_offset, w0_FitB,
        r'CMB Prior' + '\n' + r'$\Omega_m = 0.315$',
        fontsize=10, color=col_cmb_pt, ha='left', va='center',
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                  edgecolor='none', alpha=0.55))

# EU Crossover (green, at w₀=-1)
ax.text(Om_cross + label_offset, -1.0,
        r'EU Crossover' + '\n' + r'$\Omega_m \approx 0.297$',
        fontsize=10, color=col_eu_pt, ha='left', va='center',
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                  edgecolor='none', alpha=0.85))

# EU Post-Drain (green, closer to its circle)
ax.text(Om_postdrain + 0.003, w0_postdrain - 0.015,
        r'EU Post-Drain' + '\n' + r'$\Omega_m = 0.288$',
        fontsize=10, color=col_eu_pt, ha='left', va='top',
        fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                  edgecolor='none', alpha=0.85))

# ================================================================
# AXIS CONFIGURATION
# ================================================================
ax.set_xlabel(r'Forced prior $\Omega_m$')
ax.set_ylabel(r'Apparent $w_0$', labelpad=-1)
ax.set_xlim(0.265, 0.350)
ax.set_ylim(-1.15, -0.70)

ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

# ================================================================
# SAVE
# ================================================================
output_dir = Path(r'c:\Users\ricar\Desktop\PAPER 1\Current'
                  r'\Paper_I_A_Dissipation_Principle\04_Figures')
fig.savefig(output_dir / 'fig6_cpl_mirage.pdf', format='pdf')
fig.savefig(output_dir / 'fig6_cpl_mirage.png', format='png')
plt.show()
print("Figure 6 saved: fig6_cpl_mirage.pdf + fig6_cpl_mirage.png")
