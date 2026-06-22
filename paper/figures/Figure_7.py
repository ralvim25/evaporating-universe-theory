"""
Figure 7 — Posterior Comparison: Theory-Fixed vs Free-Parameter
================================================================
Paper I: The Evaporating Universe
Location: Appendix B (Free-parameter runs)

Overlaid corner plot comparing:
  - C2 (theory-fixed, baseline)  — blue
  - A1 (free-parameter, EU-blind) — green

Parameters: H₀(GKI), H₀(LKI), σ₈, S₈, Ωm
Demonstrates posterior degradation when EU parameters are freed.

Style: PRD two-column (7 inches wide), Computer Modern fonts
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
from getdist import loadMCSamples

# ================================================================
# STYLE — PRD Computer Modern
# ================================================================
plt.rcParams.update({
    'text.usetex': True,
    'font.family': 'serif',
    'font.serif': ['Computer Modern Roman'],
    'axes.labelsize': 11,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
})

# ================================================================
# PATHS
# ================================================================
chains_dir = Path(r'c:\Users\ricar\Desktop\PAPER 1\Current'
                  r'\Paper_I_A_Dissipation_Principle'
                  r'\01_Notebooks\results\MCMC_Chains')
output_dir = Path(r'c:\Users\ricar\Desktop\PAPER 1\Current'
                  r'\Paper_I_A_Dissipation_Principle\04_Figures')

# ================================================================
# LOAD CHAINS
# ================================================================
print("Loading C2 chains...")
c2 = loadMCSamples(str(chains_dir / 'RUN_C2_Final' / 'eu_NB05C_V2'),
                   settings={'ignore_rows': 0.3})

print("Loading A1 chains...")
a1 = loadMCSamples(str(chains_dir / 'RUN_A1_Final' / 'eu_NB05A'),
                   settings={'ignore_rows': 0.3})

# ================================================================
# PARAMETER SELECTION — 5 key parameters
# ================================================================
params = ['H0', 'H0_LKI', 'sigma8', 'S8', 'Omega_m']

labels = {
    'H0':      r'$H_0^{\rm GKI}$',
    'H0_LKI':  r'$H_0^{\rm LKI}$',
    'sigma8':  r'$\sigma_8$',
    'S8':      r'$S_8$',
    'Omega_m': r'$\Omega_m$',
}

# Theory values (EFT-derived)
theory_vals = {
    'H0':      68.95,
    'H0_LKI':  72.71,
    'sigma8':  0.827,
    'S8':      0.811,
    'Omega_m': 0.288,
}

n = len(params)

# ================================================================
# COLORS
# ================================================================
col_c2      = '#1565C0'    # Blue for theory-fixed
col_c2_68   = '#64B5F6'    # Blue 68%
col_c2_95   = '#BBDEFB'    # Blue 95%
col_a1      = '#2E7D32'    # Green for free-param
col_a1_68   = '#81C784'    # Green 68%
col_a1_95   = '#C8E6C9'    # Green 95%
col_theory  = '#D32F2F'    # Red for theory values

# Standardized EFT line style — same everywhere
eft_style = dict(color=col_theory, ls=':', lw=0.9, alpha=0.75)

# ================================================================
# HELPER
# ================================================================
def get_1d(samples, param):
    """1D marginalized density, normalized to peak = 1."""
    d = samples.get1DDensity(param)
    return d.x, d.P / d.P.max()

# ================================================================
# BUILD FIGURE — 5×5 triangle, two-column width
# ================================================================
fig, axes = plt.subplots(n, n, figsize=(7.0, 6.8))

# Hide upper triangle
for i in range(n):
    for j in range(n):
        if j > i:
            axes[i, j].set_visible(False)

# ================================================================
# DIAGONAL — 1D posteriors
# ================================================================
for i, p in enumerate(params):
    ax = axes[i, i]

    x_c2, y_c2 = get_1d(c2, p)
    ax.fill_between(x_c2, y_c2, alpha=0.30, color=col_c2_68, zorder=2)
    ax.plot(x_c2, y_c2, color=col_c2, lw=1.3, zorder=3)

    x_a1, y_a1 = get_1d(a1, p)
    ax.fill_between(x_a1, y_a1, alpha=0.25, color=col_a1_68, zorder=1)
    ax.plot(x_a1, y_a1, color=col_a1, lw=1.1, ls='--', zorder=3)

    if p in theory_vals:
        ax.axvline(theory_vals[p], **eft_style, zorder=4)

    ax.set_ylim(0, 1.4)
    ax.set_yticks([])
    ax.tick_params(axis='x', labelsize=7, rotation=45)

    if i < n - 1:
        ax.set_xticklabels([])

# ================================================================
# OFF-DIAGONAL — 2D contours
# ================================================================
for i in range(1, n):
    for j in range(i):
        ax = axes[i, j]
        px, py = params[j], params[i]

        # --- C2 (theory-fixed) ---
        try:
            d2 = c2.get2DDensity(px, py, normalized=True)
            levs = d2.getContourLevels([0.95, 0.68])
            ax.contourf(d2.x, d2.y, d2.P,
                       levels=[0] + list(levs) + [d2.P.max()],
                       colors=['white', col_c2_95, col_c2_68,
                               col_c2_68],
                       alpha=0.50, zorder=1)
            ax.contour(d2.x, d2.y, d2.P, levels=levs,
                      colors=[col_c2, col_c2],
                      linewidths=[0.4, 0.9], zorder=3)
        except Exception as e:
            print(f"  C2 {px}x{py}: {e}")

        # --- A1 (free-parameter) ---
        try:
            d2 = a1.get2DDensity(px, py, normalized=True)
            levs = d2.getContourLevels([0.95, 0.68])
            ax.contourf(d2.x, d2.y, d2.P,
                       levels=[0] + list(levs) + [d2.P.max()],
                       colors=['white', col_a1_95, col_a1_68,
                               col_a1_68],
                       alpha=0.35, zorder=0)
            ax.contour(d2.x, d2.y, d2.P, levels=levs,
                      colors=[col_a1, col_a1],
                      linewidths=[0.4, 0.9], linestyles='--',
                      zorder=2)
        except Exception as e:
            print(f"  A1 {px}x{py}: {e}")

        # Theory crosshairs — SAME style as diagonal
        if px in theory_vals:
            ax.axvline(theory_vals[px], **eft_style, zorder=0)
        if py in theory_vals:
            ax.axhline(theory_vals[py], **eft_style, zorder=0)

        ax.tick_params(axis='both', labelsize=7)
        ax.tick_params(axis='x', rotation=45)

        if i < n - 1:
            ax.set_xticklabels([])
        if j > 0:
            ax.set_yticklabels([])

# ================================================================
# AXIS LABELS
# ================================================================
for j in range(n):
    axes[n - 1, j].set_xlabel(labels[params[j]], fontsize=10)

for i in range(1, n):
    axes[i, 0].set_ylabel(labels[params[i]], fontsize=10)

# ================================================================
# LEGEND — in upper-right empty space
# ================================================================
legend_elements = [
    Line2D([0], [0], color=col_c2, lw=1.5,
           label=r'Theory-fixed'),
    Line2D([0], [0], color=col_a1, lw=1.3, ls='--',
           label=r'Free-parameter'),
    Line2D([0], [0], color=col_theory, lw=0.9, ls=':',
           label=r'EFT prediction'),
]

ax_leg = axes[0, n - 1]
ax_leg.set_visible(True)
ax_leg.axis('off')
ax_leg.legend(handles=legend_elements, loc='center',
              fontsize=8.5, frameon=True, fancybox=False,
              edgecolor='0.75', facecolor='white',
              framealpha=0.95, handlelength=2.0)

# Hide remaining upper triangle cells
for i in range(n):
    for j in range(i + 1, n):
        if axes[i, j].get_visible():
            axes[i, j].axis('off')

plt.tight_layout(h_pad=0.4, w_pad=0.4)
plt.subplots_adjust(hspace=0.05, wspace=0.05)

# ================================================================
# SAVE
# ================================================================
fig.savefig(output_dir / 'fig7_corner_comparison.pdf', format='pdf')
fig.savefig(output_dir / 'fig7_corner_comparison.png', format='png')
print("Figure 7 saved: fig7_corner_comparison.pdf + fig7_corner_comparison.png")
