"""
Figure 1 — The Wall and the Escape (Dual Vertical Panel)
=========================================================
Paper I: The Evaporating Universe
Section: §VII.D (The geometric wall) + §VII.A (Two-scale structure)

Top panel:    Profile likelihood Dchi2(H0) — the "geometric wall"
Bottom panel: Waterfall decomposition Planck -> GKI -> LKI -> SH0ES

Style: PRD single-column (3.375 inches wide), Computer Modern fonts
RULE:  ZERO text overlap with curves, data points, or other text.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

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
# DATA (from NB04 JSONs, code freeze)
# ================================================================
H0_best_lcdm = 68.84422110552764
H0_best_eu   = 69.54773869346734
Dchi2_at_73   = 158.40885302329082

H0_Planck     = 67.36
H0_GKI        = 68.889061
GKI_contrib   = 1.529061
LKI_contrib   = 3.81593424
H0_LKI        = 72.70499524
H0_SH0ES      = 73.17
H0_SH0ES_err  = 0.86

# Gap to SH0ES
gap_to_shoes = H0_SH0ES - H0_LKI  # ~0.465

# ================================================================
# COLOR PALETTE
# ================================================================
col_eu       = '#C62828'
col_lcdm     = '#546E7A'
col_shoes    = '#E65100'
col_gki      = '#2E7D32'
col_lki      = '#1565C0'
col_wall     = '#B71C1C'

# ================================================================
# Profile likelihood parabolas
# ================================================================
sigma_eff = (H0_SH0ES - H0_best_lcdm) / np.sqrt(Dchi2_at_73)
H0_scan = np.linspace(64, 76, 500)
Dchi2_lcdm = ((H0_scan - H0_best_lcdm) / sigma_eff)**2
Dchi2_eu   = ((H0_scan - H0_best_eu) / sigma_eff)**2

# ================================================================
# CREATE FIGURE
# ================================================================
fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(3.375, 4.8),
    gridspec_kw={'height_ratios': [1.1, 1], 'hspace': 0.30},
)

# ================================================================
# PANEL (a): Profile Likelihood
# ================================================================
line_lcdm, = ax1.plot(H0_scan, Dchi2_lcdm, color=col_lcdm, lw=1.2,
                      zorder=3)
line_eu, = ax1.plot(H0_scan, Dchi2_eu, color=col_eu, lw=1.2, ls='--',
                    zorder=3)

# SH0ES vertical
ax1.axvline(x=H0_SH0ES, color=col_shoes, ls='-', lw=0.9, alpha=0.7,
            zorder=2)

# Wall fill (gentle)
wall_mask = H0_scan >= 71.5
ax1.fill_between(H0_scan[wall_mask], 0, Dchi2_lcdm[wall_mask],
                 color=col_wall, alpha=0.05, zorder=1)

# Minimum markers (triangles)
mk_lcdm = ax1.plot(H0_best_lcdm, 0, 'v', color=col_lcdm, ms=5,
                    zorder=6)[0]
mk_eu = ax1.plot(H0_best_eu, 0, 'v', color=col_eu, ms=5, zorder=6)[0]

# --- ANNOTATIONS (zero overlap) ---

# SH0ES label: right edge, clear zone above wall
ax1.text(75.3, 175, 'SH0ES',
         fontsize=8, color=col_shoes, ha='center', va='center',
         fontweight='bold')

# Wall annotation box: far right, mid-height
ax1.text(75.0, 95,
         r'$\Delta\chi^2 = 158.4$' + '\n' + r'$(12.59\sigma)$',
         fontsize=8, color=col_wall, ha='center', va='center',
         bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                   edgecolor=col_wall, alpha=0.9, lw=0.4))

# "Geometric wall" — to the right of the SH0ES vertical line, clear zone
ax1.text(74.8, 40, 'Geometric\nwall',
         fontsize=9, color=col_wall, ha='center', va='center',
         style='italic', alpha=0.8,
         bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                   edgecolor='none', alpha=0.85))

ax1.set_xlabel(r'$H_0\;[\mathrm{km\,s^{-1}\,Mpc^{-1}}]$')
ax1.set_ylabel(r'$\Delta\chi^2$', labelpad=-2)
ax1.set_xlim(64.5, 76.8)
ax1.set_ylim(-8, 205)  # Top slightly above 200 so tick label sits inside frame

# "(a)" label — moved right to avoid sitting on the y-axis spine
ax1.text(0.06, 0.93, '(a)', transform=ax1.transAxes,
         fontsize=10, fontweight='bold', va='top')

# LEGEND with markers: triangles + lines, using H0 notation
legend_handles = [
    Line2D([0], [0], color=col_lcdm, lw=1.2, marker='v', ms=5,
           markerfacecolor=col_lcdm, label=r'$\Lambda$CDM ($68.84$)'),
    Line2D([0], [0], color=col_eu, lw=1.2, ls='--', marker='v', ms=5,
           markerfacecolor=col_eu, label=r'$H_0^{\rm GKI}$ ($69.55$)'),
]
ax1.legend(handles=legend_handles, loc='upper center',
           frameon=True, framealpha=0.95,
           edgecolor='#BDBDBD', fancybox=False, borderpad=0.3,
           handletextpad=0.4, fontsize=8,
           bbox_to_anchor=(0.35, 0.98))

ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

# ================================================================
# PANEL (b): Waterfall — 4 bars
# ================================================================
# y=3: Planck baseline (67.36)
# y=2: +GKI increment (67.36 → 68.89)
# y=1: +LKI increment (68.89 → 72.70)
# y=0: Gap to SH0ES (72.70 → 73.17, labeled 0.51σ)

y_pos = [3, 2, 1, 0]
labels = [
    'Planck\n' + r'$\Lambda$CDM',
    r'$H_0^{\rm GKI}$',
    r'$H_0^{\rm LKI}$',
    'Gap\nSH0ES',
]

bar_h = 0.50

# --- y=3: Planck baseline ---
planck_left = 66.3
ax2.barh(y_pos[0], H0_Planck - planck_left, left=planck_left,
         height=bar_h, color=col_lcdm, alpha=0.2,
         edgecolor=col_lcdm, linewidth=0.7, zorder=3)
ax2.plot([H0_Planck, H0_Planck],
         [y_pos[0] - bar_h/2, y_pos[0] + bar_h/2],
         color=col_lcdm, lw=1.2, zorder=5)
ax2.text(H0_Planck + 0.15, y_pos[0], r'$67.36$',
         fontsize=8, color=col_lcdm, ha='left', va='center')

# --- y=2: +GKI ---
ax2.barh(y_pos[1], GKI_contrib, left=H0_Planck, height=bar_h,
         color=col_gki, alpha=0.75, edgecolor=col_gki, linewidth=0.5,
         zorder=3)
ax2.text(H0_Planck + GKI_contrib / 2, y_pos[1],
         f'+{GKI_contrib:.2f}',
         fontsize=8, color='white', ha='center', va='center',
         fontweight='bold', zorder=4)
ax2.text(H0_GKI + 0.15, y_pos[1], f'{H0_GKI:.2f}',
         fontsize=8, color=col_gki, ha='left', va='center')

# Connector: Planck → GKI
ax2.plot([H0_Planck, H0_Planck],
         [y_pos[0] - bar_h/2, y_pos[1] + bar_h/2],
         ls=':', color='gray', lw=0.5, zorder=2)

# --- y=1: +LKI ---
ax2.barh(y_pos[2], LKI_contrib, left=H0_GKI, height=bar_h,
         color=col_lki, alpha=0.75, edgecolor=col_lki, linewidth=0.5,
         zorder=3)
ax2.text(H0_GKI + LKI_contrib / 2, y_pos[2],
         f'+{LKI_contrib:.2f}',
         fontsize=8, color='white', ha='center', va='center',
         fontweight='bold', zorder=4)
ax2.text(H0_LKI + 0.15, y_pos[2], f'{H0_LKI:.2f}',
         fontsize=8, color=col_lki, ha='left', va='center',
         bbox=dict(boxstyle='round,pad=0.08', facecolor='white',
                   edgecolor='none', alpha=0.85),
         zorder=6)

# Connector: GKI → LKI
ax2.plot([H0_GKI, H0_GKI],
         [y_pos[1] - bar_h/2, y_pos[2] + bar_h/2],
         ls=':', color='gray', lw=0.5, zorder=2)

# --- y=0: Gap to SH0ES (0.51σ) ---
ax2.barh(y_pos[3], gap_to_shoes, left=H0_LKI, height=bar_h,
         color=col_shoes, alpha=0.25, edgecolor=col_shoes, linewidth=0.5,
         zorder=3, hatch='///')
# 0.51σ label: same style as other value labels, to the right of gap bar
ax2.text(H0_SH0ES + 0.15, y_pos[3], r'$0.51\sigma$',
         fontsize=8, color=col_shoes, ha='left', va='center',
         bbox=dict(boxstyle='round,pad=0.08', facecolor='white',
                   edgecolor='none', alpha=0.85),
         zorder=6)

# Connector: LKI → Gap
ax2.plot([H0_LKI, H0_LKI],
         [y_pos[2] - bar_h/2, y_pos[3] + bar_h/2],
         ls=':', color='gray', lw=0.5, zorder=2)

# SH0ES band
ax2.axvspan(H0_SH0ES - H0_SH0ES_err, H0_SH0ES + H0_SH0ES_err,
            color=col_shoes, alpha=0.08, zorder=1)
ax2.axvline(x=H0_SH0ES, color=col_shoes, ls='-', lw=0.7,
            alpha=0.5, zorder=2)

# Y-axis
ax2.set_yticks(y_pos)
ax2.set_yticklabels(labels, fontsize=8.5)
ax2.set_xlabel(r'$H_0\;[\mathrm{km\,s^{-1}\,Mpc^{-1}}]$')
ax2.set_xlim(66.3, 75.5)
ax2.set_ylim(-0.5, 4.3)

ax2.spines['left'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.tick_params(axis='y', length=0)

ax2.text(0.06, 0.95, '(b)', transform=ax2.transAxes,
         fontsize=10, fontweight='bold', va='top')

ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))

fig.align_ylabels([ax1, ax2])

# ================================================================
# SAVE
# ================================================================
output_dir = r'c:\Users\ricar\Desktop\PAPER 1\Current\Paper_I_A_Dissipation_Principle\04_Figures'
fig.savefig(f'{output_dir}\\fig1_wall_escape.pdf', format='pdf')
fig.savefig(f'{output_dir}\\fig1_wall_escape.png', format='png')
plt.show()
print("Figure 1 saved: fig1_wall_escape.pdf + fig1_wall_escape.png")
