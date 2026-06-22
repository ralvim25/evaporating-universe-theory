"""
Figure 2 — The Anatomy of Evaporation (Dual Vertical Panel)
============================================================
Paper I: The Evaporating Universe
Section: §II (Theoretical Framework) + §VIII (N-body Validation)

Top panel:    epsilon(z) logistic coupling vs redshift
Bottom panel: f_cdm(z) CDM survival fraction — N-body snapshots vs analytical ODE

Style: PRD single-column (3.375 inches wide), Computer Modern fonts
RULE:  ZERO text overlap with curves, data points, or other text.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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
# EU PARAMETERS (from NB01, code freeze)
# ================================================================
eps_IR = 0.042638612831103685
z_trans = 5.986116748871208
b = 19.0 / 36.0  # = 0.527...
lam = 2.0 / 3.0

# ================================================================
# ANALYTICAL FUNCTIONS
# ================================================================
def epsilon(z):
    """Logistic coupling amplitude with Heaviside cutoff (Eq. 5 of paper)."""
    return np.where(
        z <= z_trans,
        eps_IR / (1.0 + ((1.0 + z) / (1.0 + z_trans))**(1.0 / b)),
        0.0
    )


def f_cdm_analytic(z_array):
    """
    CDM survival fraction from numerical integration of ODE:
    df/d(ln a) = -lambda * epsilon(z) * f
    """
    z_start = 50.0
    n_steps = 50000
    z_fine = np.linspace(z_start, 0.0, n_steps)

    f = 1.0
    f_values = np.zeros(n_steps)
    f_values[0] = 1.0

    for i in range(1, n_steps):
        z_mid = 0.5 * (z_fine[i-1] + z_fine[i])
        dz = z_fine[i-1] - z_fine[i]  # positive (decreasing z)
        eps_val = epsilon(z_mid)
        dlna = dz / (1.0 + z_mid)
        f = f * np.exp(-lam * eps_val * dlna)
        f_values[i] = f

    return np.interp(z_array, z_fine[::-1], f_values[::-1])


# ================================================================
# N-BODY DATA (from NB07_results.json)
# ================================================================
# Production run (full TreePM, gravity ON)
prod_snapshots = [
    {'z': 49.0,    'f': 1.0},
    {'z': 6.9906,  'f': 1.0},
    {'z': 6.0172,  'f': 1.0},
    {'z': 5.5011,  'f': 0.9758038446674817},
    {'z': 5.0228,  'f': 0.9741807345507449},
    {'z': 4.4952,  'f': 0.9722768905407981},
    {'z': 4.0137,  'f': 0.9704451309065601},
    {'z': 3.5051,  'f': 0.9684193047399072},
    {'z': 2.9867,  'f': 0.9663131376354607},
    {'z': 2.5011,  'f': 0.9643204692821956},
    {'z': 2.0049,  'f': 0.9623038383615866},
    {'z': 1.5015,  'f': 0.9603449866837105},
    {'z': 1.0043,  'f': 0.9585683994588211},
    {'z': 0.4992,  'f': 0.9570062946759096},
    {'z': 0.0,     'f': 0.9557888379620736},
]

# Null test (gravity OFF, drain ON)
null_snapshots = [
    {'z': 49.0,    'f': 1.0},
    {'z': 6.9906,  'f': 1.0},
    {'z': 6.0172,  'f': 1.0},
    {'z': 5.5011,  'f': 0.9759468637539533},
    {'z': 5.0228,  'f': 0.9743216998124733},
    {'z': 4.4952,  'f': 0.9724135676528443},
    {'z': 4.0137,  'f': 0.9705757604844771},
    {'z': 3.5051,  'f': 0.9685497711315473},
    {'z': 2.9867,  'f': 0.9664221493860837},
    {'z': 2.5011,  'f': 0.9644067965851837},
    {'z': 2.0049,  'f': 0.9623721252868815},
    {'z': 1.5015,  'f': 0.9603932275123078},
    {'z': 1.0043,  'f': 0.9585963070077556},
    {'z': 0.4992,  'f': 0.9570146290987912},
    {'z': 0.0,     'f': 0.9557808196478058},
]

prod_z = np.array([s['z'] for s in prod_snapshots])
prod_f = np.array([s['f'] for s in prod_snapshots])
null_z = np.array([s['z'] for s in null_snapshots])
null_f = np.array([s['f'] for s in null_snapshots])

# ================================================================
# ANALYTICAL CURVES
# ================================================================
z_curve = np.linspace(0, 10, 2000)
eps_curve = epsilon(z_curve)
f_curve = f_cdm_analytic(z_curve)

# ================================================================
# COLOR PALETTE
# ================================================================
col_eu       = '#C62828'   # EU primary (deep red)
col_nbody    = '#1B5E20'   # N-body production (deep green)
col_null     = '#E65100'   # Null test (deep orange)
col_lcdm     = '#78909C'   # LCDM passive (blue-gray)
col_ztrans   = '#7B1FA2'   # z_trans (purple)
col_fill     = '#FFCDD2'   # Light fill under epsilon

# ================================================================
# CREATE FIGURE
# ================================================================
fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(3.375, 4.6),
    gridspec_kw={'height_ratios': [1, 1.2], 'hspace': 0.06},
    sharex=True
)

# ================================================================
# PANEL (a): epsilon(z)
# ================================================================
ax1.plot(z_curve, eps_curve * 100, color=col_eu, lw=1.3, zorder=3)
ax1.fill_between(z_curve, 0, eps_curve * 100,
                 color=col_fill, alpha=0.35, zorder=1)

# --- Horizontal dashed at eps_IR (ABOVE the line, clear zone) ---
ax1.axhline(y=eps_IR * 100, color=col_lcdm, ls='--', lw=0.5, alpha=0.6)

# Label for eps_IR: right of (a) label, along the dashed line
ax1.text(2.5, 4.65, r'$\varepsilon_{\rm IR} = 4.26\%$',
         fontsize=9.5, color=col_lcdm, va='center', ha='left')

# --- z_trans vertical ---
ax1.axvline(x=z_trans, color=col_ztrans, ls=':', lw=0.7, alpha=0.7)

# Label for z_trans: top-right clear zone (z>6, y>4)
ax1.text(7.8, 4.5, r'$z_{\rm trans} \approx 5.99$',
         fontsize=9.5, color=col_ztrans, ha='center', va='center',
         bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                   edgecolor=col_ztrans, alpha=0.9, lw=0.4))

# "epsilon = 0 (exact)" in the dead zone right of z_trans
ax1.text(8.3, 1.8, r'$\varepsilon = 0$' + '\n(exact)',
         fontsize=9.5, color=col_lcdm, ha='center', va='center',
         style='italic')

# "Ghost condensate active" — centered in the pink fill region
ax1.text(2.8, 1.8, 'Ghost condensate\nactive',
         fontsize=10, color='#B71C1C', ha='center', va='center',
         style='italic', alpha=0.85)

ax1.set_ylabel(r'Coupling $\varepsilon(z)$ [%]')
ax1.set_ylim(-0.3, 5.5)
ax1.set_xlim(-0.3, 10.3)

ax1.text(0.03, 0.92, '(a)', transform=ax1.transAxes,
         fontsize=10, fontweight='bold', va='top')

ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

# Hide x tick labels on top panel (shared axis)
plt.setp(ax1.get_xticklabels(), visible=False)

# ================================================================
# PANEL (b): f_cdm(z)
# ================================================================

# --- Analytical curve ---
ax2.plot(z_curve, f_curve, color=col_eu, lw=1.2, zorder=3,
         label='Analytical ODE')

# --- Null test markers (plot FIRST, behind production) ---
# Slightly larger, hollow — so they peek out behind production
ax2.scatter(null_z[3:], null_f[3:],  # skip z>z_trans (all f=1 there)
            s=30, marker='D', facecolors='none',
            edgecolors=col_null, linewidths=0.8, zorder=4,
            label='Null test (no gravity)')

# --- Production markers (on top) ---
ax2.scatter(prod_z[3:], prod_f[3:],  # skip z>z_trans
            s=18, marker='o', facecolors=col_nbody,
            edgecolors='none', zorder=5,
            label=r'Production ($N$-body)')

# --- Exact points at f=1 (z > z_trans): show 3 snapshots ---
ax2.scatter(prod_z[:3], prod_f[:3],
            s=18, marker='o', facecolors=col_nbody,
            edgecolors='none', zorder=5)
ax2.scatter(null_z[:3], null_f[:3],
            s=30, marker='D', facecolors='none',
            edgecolors=col_null, linewidths=0.8, zorder=4)

# --- Horizontal dashed at f = 1.0 ---
ax2.axhline(y=1.0, color=col_lcdm, ls='--', lw=0.4, alpha=0.4)

# --- z_trans vertical ---
ax2.axvline(x=z_trans, color=col_ztrans, ls=':', lw=0.7, alpha=0.7)

# --- ANNOTATIONS (all in clear zones, NO overlap with data) ---

# "f_cdm = 1 (exact)": top-right, above the f=1 line, right of z_trans
ax2.text(9.2, 1.005, r'$f_{\rm cdm} \equiv 1$',
         fontsize=9, color=col_lcdm, ha='center', va='bottom',
         style='italic')

# "f_cdm(0) = 0.9558": center-bottom, arrow curves BELOW the curve
ax2.annotate(
    r'$f_{\rm cdm}(0) = 0.9558$' + '\n' + r'$(4.42\%\;{\rm drain})$',
    xy=(0.0, 0.9558), xytext=(3.8, 0.9525),
    fontsize=8.5, color=col_eu,
    arrowprops=dict(arrowstyle='->', color=col_eu, lw=0.6,
                    connectionstyle='arc3,rad=-0.3'),
    ha='center', va='center',
    bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
              edgecolor=col_eu, alpha=0.95, lw=0.4))

# "Gravitational decoupling" annotation: top-center, above curve
# At z~3, the curve is at f~0.968. Put text well above at f~1.006
ax2.text(3.5, 1.011,
         r'Gravitational decoupling:' + '\n' + r'$|\Delta M|_{\rm max} < 0.015\%$',
         fontsize=8.5, color=col_lcdm, ha='center', va='center',
         style='italic')

# --- Legend: upper-right, inside the flat f=1 zone (z>6.5, f~0.99) ---
ax2.legend(loc='center right', frameon=True, framealpha=0.95,
           edgecolor='#BDBDBD', fancybox=False, borderpad=0.4,
           handletextpad=0.4, fontsize=6.5,
           bbox_to_anchor=(0.99, 0.18))

ax2.set_xlabel(r'Redshift $z$')
ax2.set_ylabel(r'CDM survival fraction $f_{\rm cdm}(z)$')
ax2.set_ylim(0.947, 1.018)

ax2.text(0.03, 0.92, '(b)', transform=ax2.transAxes,
         fontsize=10, fontweight='bold', va='top')

ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))

# Align y-axis labels across both panels
fig.align_ylabels([ax1, ax2])

# ================================================================
# SAVE
# ================================================================
output_dir = r'c:\Users\ricar\Desktop\PAPER 1\Current\Paper_I_A_Dissipation_Principle\04_Figures'

fig.savefig(f'{output_dir}\\fig2_anatomy.pdf', format='pdf')
fig.savefig(f'{output_dir}\\fig2_anatomy.png', format='png')

plt.show()
print("Figure 2 saved: fig2_anatomy.pdf + fig2_anatomy.png")
