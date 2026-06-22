"""
EU NB05 — RUN C2 Chain Analysis
================================
Run C2: EU Theory-Fixed (Δk=0), Perturbations ON (ide_perturbations=1)
Adapted from RUN_C1 analyze_runC.py

Key differences from C1:
  - Perturbations ON → σ₈ and S₈ are NATIVE (no post-processing needed)
  - H0_LKI, fcdm_z0, I_GKI are computed natively by shoes_lki.py/eu_derived.py
  - Chain files: eu_NB05C_V2.{1..8}.txt
  - Includes perturbation null test: C2 vs C1 comparison

Output: analysis_RUN_C2.json + analysis_RUN_C2.md
"""
import numpy as np
import json
import os
from scipy.integrate import solve_ivp
from datetime import datetime

# ── Configuration ──
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CHAIN_DIR = os.path.dirname(SCRIPT_DIR)  # RUN_C2_final/
CHAIN_PREFIX = "eu_NB05C_V2"
N_CHAINS = 8
BURN_FRAC = 0.30  # 30% burn-in

# ── Load chains ──
chains = []
for i in range(1, N_CHAINS + 1):
    path = os.path.join(CHAIN_DIR, f"{CHAIN_PREFIX}.{i}.txt")
    d = np.loadtxt(path)
    chains.append(d)

# Header
with open(os.path.join(CHAIN_DIR, f"{CHAIN_PREFIX}.1.txt")) as f:
    header = f.readline().strip('#').split()
col = {name: idx for idx, name in enumerate(header)}

# Combine all chains with burn-in removal
all_data = []
for c in chains:
    burn = int(len(c) * BURN_FRAC)
    all_data.append(c[burn:])
data = np.vstack(all_data)
weights = data[:, col['weight']]

# Progress file
progress_path = os.path.join(CHAIN_DIR, f"{CHAIN_PREFIX}.progress")
if os.path.exists(progress_path):
    progress_lines = open(progress_path).readlines()
    last_progress = progress_lines[-1].strip() if len(progress_lines) > 1 else "N/A"
else:
    last_progress = "N/A"

# ── Weighted statistics ──
def wstats(vals, w):
    mean = np.average(vals, weights=w)
    var = np.average((vals - mean)**2, weights=w)
    std = np.sqrt(var)
    sorted_idx = np.argsort(vals)
    sorted_vals = vals[sorted_idx]
    sorted_w = w[sorted_idx]
    cumw = np.cumsum(sorted_w)
    cumw /= cumw[-1]
    q16 = sorted_vals[np.searchsorted(cumw, 0.16)]
    q50 = sorted_vals[np.searchsorted(cumw, 0.50)]
    q84 = sorted_vals[np.searchsorted(cumw, 0.84)]
    return mean, std, q16, q50, q84

# ── Output buffer (for both console and markdown) ──
output_lines = []
def pr(line=""):
    print(line)
    output_lines.append(line)

# ═══════════════════════════════════════════════════════════
# ANALYSIS
# ═══════════════════════════════════════════════════════════

pr("=" * 70)
pr("RUN C2 — EU THEORY-FIXED (Δk=0, perturbations ON)")
pr("  ε_IR=0.04264, b=19/36, z_trans=5.986, ide_perturbations=1")
pr("=" * 70)

# Parse R-1 from progress
if last_progress != "N/A":
    parts = last_progress.split()
    r1_val = float(parts[3]) if len(parts) > 3 else None
    r1_cl = parts[4] if len(parts) > 4 else "NaN"
    n_samples_total = int(float(parts[0]))
    acceptance = float(parts[2])
    pr(f"  Last checkpoint: N={n_samples_total}, acc={acceptance:.3f}, "
       f"R-1={r1_val:.4f}, R-1_cl={r1_cl}")
else:
    r1_val = None
    n_samples_total = 0
    acceptance = 0

# Per-chain stats
pr("\n### PER-CHAIN STATISTICS")
chain_info = []
for i, c in enumerate(chains):
    burn = int(len(c) * BURN_FRAC)
    post_burn = c[burn:]
    w = post_burn[:, col['weight']]
    n_eff = np.sum(w)**2 / np.sum(w**2)
    pr(f"  Chain {i+1}: {len(c)} raw → {len(post_burn)} post-burn, "
       f"Σw={np.sum(w):.0f}, n_eff={n_eff:.0f}")
    chain_info.append({'raw': len(c), 'post_burn': len(post_burn),
                       'weight_sum': float(np.sum(w)), 'n_eff': float(n_eff)})

total_weight = np.sum(weights)
n_eff_total = total_weight**2 / np.sum(weights**2)
pr(f"\n  TOTAL: {len(data)} rows post-burn, Σw={total_weight:.0f}, n_eff={n_eff_total:.0f}")

# ── Parameter posteriors ──
# C2 has only ΛCDM sampled params (EU fixed via extra_args)
params_list = [
    ('omega_b',      'ω_b'),
    ('omega_cdm',    'ω_cdm'),
    ('theta_s_100',  '100θ_s'),
    ('tau_reio',     'τ_reio'),
    ('logA',         'ln(10¹⁰As)'),
    ('n_s',          'n_s'),
    ('A_planck',     'A_planck'),
    ('H0',           'H₀'),
    ('sigma8',       'σ₈'),
    ('Omega_m',      'Ω_m'),
    ('S8',           'S₈'),
    ('rdrag',        'r_drag'),
    ('H0_LKI',       'H₀_LKI'),
    ('fcdm_z0',      'f_cdm(z=0)'),
    ('I_GKI',        'I_GKI'),
]

# Reference values
planck_ref = {
    'H0': 67.36, 'sigma8': 0.8111, 'Omega_m': 0.3153, 'S8': 0.834,
    'omega_b': 0.02237, 'omega_cdm': 0.1200, 'tau_reio': 0.0544,
    'n_s': 0.9649, 'rdrag': 147.09, 'theta_s_100': 1.04092,
    'logA': 3.044, 'A_planck': 1.0,
}

# NB01 UV predictions for derived params
nb01_ref = {
    'H0': 68.90,        # H0_GKI analytic (NB02)
    'H0_LKI': 72.49,    # LKI analytic (NB02)
    'fcdm_z0': 0.9558,  # CDM survival fraction
    'I_GKI': 0.06785,   # Kinematic integral
}

pr("\n### PARAMETER POSTERIORS (30% burn-in removed)")
pr(f"{'Param':<12} {'Mean':>10} {'Std':>10} {'q16':>10} {'q50':>10} "
   f"{'q84':>10} {'Planck':>10} {'Δ/σ':>8}")
pr("-" * 82)

results = {}
for pname, label in params_list:
    if pname not in col:
        continue
    vals = data[:, col[pname]]
    mean, std, q16, q50, q84 = wstats(vals, weights)
    ref = planck_ref.get(pname, None)
    if ref and std > 0:
        dsig = (mean - ref) / std
        ref_str = f"{ref:.5f}"
        dsig_str = f"{dsig:+.2f}σ"
    else:
        ref_str = "—"
        dsig_str = "—"
    pr(f"{label:<12} {mean:>10.5f} {std:>10.5f} {q16:>10.5f} {q50:>10.5f} "
       f"{q84:>10.5f} {ref_str:>10} {dsig_str:>8}")
    results[pname] = {'mean': float(mean), 'std': float(std),
                      'q16': float(q16), 'q50': float(q50), 'q84': float(q84)}

# ── UV Theory Comparison ──
pr("\n### NB01 UV THEORY COMPARISON (zero free parameters)")
pr(f"{'Param':<12} {'MCMC C2':>12} {'±σ':>8} {'NB01 UV':>12} {'Δ/σ':>8}")
pr("-" * 56)

for pname, uv_val in nb01_ref.items():
    if pname in results:
        r = results[pname]
        if r['std'] > 1e-10:
            dsig = (r['mean'] - uv_val) / r['std']
            dsig_str = f"{dsig:>+8.2f}σ"
        else:
            dsig_str = "   exact"
        pr(f"{pname:<12} {r['mean']:>12.4f} {r['std']:>8.4f} {uv_val:>12.4f} {dsig_str}")

# ── Correlation matrix ──
pr("\n### CORRELATION MATRIX (key params)")
corr_params = [p for p in ['H0', 'sigma8', 'S8', 'Omega_m', 'omega_b',
                            'omega_cdm', 'n_s', 'tau_reio'] if p in col]
corr_data = np.column_stack([data[:, col[p]] for p in corr_params])
wcov = np.cov(corr_data.T, aweights=weights)
wstd_arr = np.sqrt(np.diag(wcov))
wcorr = wcov / np.outer(wstd_arr, wstd_arr)

pr(f"{'':>10}" + "".join(f" {l:>9}" for l in corr_params))
for i, l in enumerate(corr_params):
    row = f"{l:>10}" + "".join(f" {wcorr[i,j]:>9.3f}" for j in range(len(corr_params)))
    pr(row)

# ── Chi² analysis ──
pr("\n### CHI² BEST-FIT")
chi2_total = data[:, col['chi2']]
best_idx = np.argmin(chi2_total)
pr(f"  Best χ² = {chi2_total[best_idx]:.2f}")

chi2_components = [
    'chi2__planck_NPIPE_highl_CamSpec.TTTEEE',
    'chi2__planck_2018_lowl.TT',
    'chi2__planck_2018_lowl.EE',
    'chi2__planck_2018_lensing.clik',
    'chi2__bao.desi_dr2',
    'chi2__sn.pantheonplus',
    'chi2__eu_derived.EU_Derived',
]
chi2_results = {}
for chi2_name in chi2_components:
    if chi2_name in col:
        short = chi2_name.replace('chi2__', '')
        val = data[best_idx, col[chi2_name]]
        pr(f"    {short}: {val:.2f}")
        chi2_results[short] = float(val)

for agg in ['chi2__CMB', 'chi2__BAO', 'chi2__SN']:
    if agg in col:
        val = data[best_idx, col[agg]]
        pr(f"    {agg}: {val:.2f}")
        chi2_results[agg] = float(val)

# ── Best-fit parameter values ──
pr("\n### BEST-FIT PARAMETER VALUES (at minimum χ²)")
for pname, label in params_list:
    if pname in col:
        pr(f"    {label:<12} = {data[best_idx, col[pname]]:.6f}")

# ═══════════════════════════════════════════════════════════
# H₀_LKI — Using native chain values (computed by eu_derived/shoes_lki)
# ═══════════════════════════════════════════════════════════

pr("\n" + "=" * 70)
pr("H₀_LKI — NATIVE CHAIN VALUES (from eu_derived/shoes_lki)")
pr("=" * 70)

if 'H0_LKI' in col:
    h0_lki_vals = data[:, col['H0_LKI']]
    h0_lki_mean, h0_lki_std, h0_lki_16, h0_lki_50, h0_lki_84 = wstats(h0_lki_vals, weights)
    h0_gki_mean, h0_gki_std, _, _, _ = wstats(data[:, col['H0']], weights)
    
    boost_mean = h0_lki_mean / h0_gki_mean
    dh_void = h0_lki_mean - h0_gki_mean
    
    pr(f"\n  H₀_GKI (global) = {h0_gki_mean:.2f} ± {h0_gki_std:.2f} km/s/Mpc")
    pr(f"  H₀_LKI (local)  = {h0_lki_mean:.2f} ± {h0_lki_std:.2f} km/s/Mpc")
    pr(f"  δH₀_void         = {dh_void:+.2f} km/s/Mpc")
    pr(f"  Boost factor      = {boost_mean:.4f}")
    pr(f"  q16={h0_lki_16:.2f}, q50={h0_lki_50:.2f}, q84={h0_lki_84:.2f}")
else:
    # Fallback: compute via Riccati ODE (same as C1)
    pr("  [WARN] H0_LKI not in chains — computing via Riccati ODE")
    delta_obs_KBC = -0.46
    b_galaxy = 1.0
    
    def growth_rate_riccati(Om_m0):
        def ode(lna, f_val):
            a = np.exp(lna)
            Om_a = Om_m0 * a**(-3) / (Om_m0 * a**(-3) + (1 - Om_m0))
            Om_DE_a = 1.0 - Om_a
            return [-(f_val[0]**2) - (0.5 + 1.5*Om_DE_a) * f_val[0] + 1.5 * Om_a]
        sol = solve_ivp(ode, [-5, 0], [1.0], rtol=1e-10, atol=1e-12)
        return float(sol.y[0, -1])
    
    h0_vals = data[:, col['H0']]
    om_vals = data[:, col['Omega_m']]
    om_grid = np.linspace(om_vals.min() - 0.001, om_vals.max() + 0.001, 50)
    f_grid = np.array([growth_rate_riccati(om) for om in om_grid])
    f_interp = np.interp(om_vals, om_grid, f_grid)
    delta_true_arr = delta_obs_KBC / (1.0 + f_interp / b_galaxy)
    dH_over_H_arr = -1.0/3.0 * delta_true_arr * f_interp
    h0_lki_vals = h0_vals * (1.0 + dH_over_H_arr)
    h0_lki_mean, h0_lki_std, h0_lki_16, h0_lki_50, h0_lki_84 = wstats(h0_lki_vals, weights)
    h0_gki_mean, h0_gki_std, _, _, _ = wstats(h0_vals, weights)
    dh_void = h0_lki_mean - h0_gki_mean
    
    pr(f"\n  H₀_GKI (global) = {h0_gki_mean:.2f} ± {h0_gki_std:.2f} km/s/Mpc")
    pr(f"  H₀_LKI (local)  = {h0_lki_mean:.2f} ± {h0_lki_std:.2f} km/s/Mpc")
    pr(f"  δH₀_void         = {dh_void:+.2f} km/s/Mpc")

# Also get native fcdm and I_GKI if available
fcdm_stats = None
igki_stats = None
if 'fcdm_z0' in col:
    fcdm_vals = data[:, col['fcdm_z0']]
    fcdm_stats = wstats(fcdm_vals, weights)
    pr(f"  fcdm(z=0)         = {fcdm_stats[0]:.6f} ± {fcdm_stats[1]:.6f}  (NB01: 0.9558)")
if 'I_GKI' in col:
    igki_vals = data[:, col['I_GKI']]
    igki_stats = wstats(igki_vals, weights)
    pr(f"  I_GKI             = {igki_stats[0]:.6f} ± {igki_stats[1]:.6f}  (NB01: 0.06785)")

# ═══════════════════════════════════════════════════════════
# TENSION DIAGNOSTICS
# ═══════════════════════════════════════════════════════════

pr("\n" + "=" * 70)
pr("TENSION DIAGNOSTICS — EU Theory-Fixed (C2) vs Observations")
pr("=" * 70)

# H₀ tension
shoes_h0, shoes_err = 73.17, 0.86
planck_h0, planck_h0_err = 67.36, 0.54

tension_shoes_lki = abs(h0_lki_mean - shoes_h0) / np.sqrt(h0_lki_std**2 + shoes_err**2)
tension_shoes_gki = abs(h0_gki_mean - shoes_h0) / np.sqrt(h0_gki_std**2 + shoes_err**2)
tension_planck_shoes = abs(planck_h0 - shoes_h0) / np.sqrt(planck_h0_err**2 + shoes_err**2)

# Freedman et al. 2024
H0_TRGB, err_TRGB = 69.85, 1.75
H0_JAGB, err_JAGB = 67.96, 2.09
tension_trgb = abs(h0_lki_mean - H0_TRGB) / np.sqrt(h0_lki_std**2 + err_TRGB**2)
tension_jagb = abs(h0_lki_mean - H0_JAGB) / np.sqrt(h0_lki_std**2 + err_JAGB**2)

pr(f"\n  --- H₀ TENSION (LKI = local prediction for distance ladder) ---")
pr(f"  EU H₀_LKI:     {h0_lki_mean:.2f} ± {h0_lki_std:.2f}  ← EU (0 free params)")
pr(f"  Planck ΛCDM:    {planck_h0} ± {planck_h0_err}")
pr(f"  SH0ES 2024:     {shoes_h0} ± {shoes_err}")
pr(f"  TRGB (Freed.):  {H0_TRGB} ± {err_TRGB}")
pr(f"  JAGB (Freed.):  {H0_JAGB} ± {err_JAGB}")
pr(f"")
pr(f"  EU H₀_LKI ↔ SH0ES:     {tension_shoes_lki:.2f}σ  ✅")
pr(f"  EU H₀_LKI ↔ TRGB:      {tension_trgb:.2f}σ  ✅")
pr(f"  EU H₀_LKI ↔ JAGB:      {tension_jagb:.2f}σ  ✅")
pr(f"  EU H₀_GKI ↔ SH0ES:     {tension_shoes_gki:.2f}σ  (global, for reference)")
pr(f"  Planck ΛCDM ↔ SH0ES:   {tension_planck_shoes:.2f}σ  ❌")
pr(f"")
pr(f"  H₀ tension reduction: {tension_planck_shoes:.1f}σ → {tension_shoes_lki:.1f}σ "
   f"({(1 - tension_shoes_lki/tension_planck_shoes)*100:.0f}%)")

# S₈ tension
s8_mean = results['S8']['mean']
s8_std = results['S8']['std']
planck_s8, planck_s8_err = 0.834, 0.016
des_s8, des_err = 0.776, 0.017
tension_des = abs(s8_mean - des_s8) / np.sqrt(s8_std**2 + des_err**2)
planck_des_tension = abs(planck_s8 - des_s8) / np.sqrt(planck_s8_err**2 + des_err**2)

pr(f"\n  --- S₈ TENSION ---")
pr(f"  EU Run C2:    S₈ = {s8_mean:.4f} ± {s8_std:.4f}  (perturbations ON, native)")
pr(f"  Planck ΛCDM:  S₈ = {planck_s8} ± {planck_s8_err}")
pr(f"  DES-Y3:       S₈ = {des_s8} ± {des_err}")
pr(f"  EU↔DES:       {tension_des:.2f}σ  (was {planck_des_tension:.1f}σ)")
pr(f"  → Tension reduction: {(1 - tension_des/planck_des_tension)*100:.0f}%")

# ═══════════════════════════════════════════════════════════
# PERTURBATION NULL TEST: C2 vs C1
# ═══════════════════════════════════════════════════════════

pr("\n" + "=" * 70)
pr("PERTURBATION NULL TEST: C2 (ON) vs C1 (OFF)")
pr("=" * 70)

# C1 reference values (from RUN_C1_Final analysis)
c1_ref = {
    'H0': (68.887, 0.281),
    'sigma8': (0.8276, 0.0059),
    'S8': (0.8113, 0.0080),
    'Omega_m': (0.2883, 0.0034),
    'omega_cdm': (0.11927, 0.00063),
    'H0_LKI': (72.42, 0.28),
}

pr(f"\n  {'Param':<10} {'C1 (OFF)':>14} {'C2 (ON)':>14} {'Δ':>10} {'Δ/σ':>8}")
pr(f"  {'-'*60}")
for pname, (c1_mean, c1_std) in c1_ref.items():
    if pname in results:
        c2_mean = results[pname]['mean']
        c2_std = results[pname]['std']
        delta = c2_mean - c1_mean
        # Use larger of the two errors for significance
        sig = abs(delta) / max(c1_std, c2_std) if max(c1_std, c2_std) > 0 else 0
        pr(f"  {pname:<10} {c1_mean:>10.4f}±{c1_std:.4f}  {c2_mean:>10.4f}±{c2_std:.4f}  "
           f"{delta:>+10.4f}  {sig:>6.2f}σ")

pr(f"\n  Expected: All Δ/σ < 0.5 (perturbations cancel for w=-1, Valiviita+2008)")

# ═══════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════

pr("\n" + "=" * 70)
pr("SUMMARY — RUN C2: EU THEORY-FIXED (Δk=0, perturbations ON)")
pr("=" * 70)
pr(f"  H₀_GKI (global) = {h0_gki_mean:.2f} ± {h0_gki_std:.2f} km/s/Mpc")
pr(f"  H₀_LKI (local)  = {h0_lki_mean:.2f} ± {h0_lki_std:.2f} km/s/Mpc  ← EU prediction")
pr(f"  σ₈               = {results['sigma8']['mean']:.4f} ± {results['sigma8']['std']:.4f}")
pr(f"  S₈               = {s8_mean:.4f} ± {s8_std:.4f}")
pr(f"  Ω_m              = {results['Omega_m']['mean']:.4f} ± {results['Omega_m']['std']:.4f}")
if fcdm_stats:
    pr(f"  fcdm(z=0)        = {fcdm_stats[0]:.6f} ± {fcdm_stats[1]:.6f}")
pr(f"")
pr(f"  H₀ tension vs SH0ES: {tension_shoes_lki:.2f}σ  (ΛCDM: {tension_planck_shoes:.1f}σ)")
pr(f"  S₈ tension vs DES:   {tension_des:.2f}σ  (ΛCDM: {planck_des_tension:.1f}σ)")
pr(f"  R-1 = {r1_val:.4f}" if r1_val else "  R-1 = N/A")
pr(f"  Total effective samples: {n_eff_total:.0f}")
pr(f"  Perturbations: ON (ide_perturbations=1)")
pr(f"  Best χ² = {chi2_total[best_idx]:.2f}")

# ═══════════════════════════════════════════════════════════
# SAVE JSON
# ═══════════════════════════════════════════════════════════

results['_meta'] = {
    'run': 'C2',
    'description': 'EU Theory-Fixed (Δk=0), perturbations ON',
    'eu_epsilon_ir': 0.04264,
    'eu_z_trans': 5.986,
    'eu_b': 0.52778,
    'eu_has_ide_perturbations': 1,
    'total_samples_post_burn': int(len(data)),
    'total_weighted_samples': float(total_weight),
    'n_eff_total': float(n_eff_total),
    'R_minus_1': float(r1_val) if r1_val else None,
    'acceptance_rate': float(acceptance),
    'burn_fraction': BURN_FRAC,
    'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
}
results['_tensions'] = {
    'H0_LKI_vs_shoes': float(tension_shoes_lki),
    'H0_LKI_vs_TRGB': float(tension_trgb),
    'H0_LKI_vs_JAGB': float(tension_jagb),
    'H0_GKI_vs_shoes': float(tension_shoes_gki),
    'S8_vs_DES': float(tension_des),
    'LCDM_H0_vs_shoes': float(tension_planck_shoes),
    'LCDM_S8_vs_DES': float(planck_des_tension),
}
results['_chi2_bestfit'] = chi2_results
results['_chi2_bestfit']['total'] = float(chi2_total[best_idx])

out_json = os.path.join(SCRIPT_DIR, "analysis_RUN_C2.json")
with open(out_json, 'w') as f:
    json.dump(results, f, indent=2)
pr(f"\n[OK] JSON saved: {out_json}")

# ═══════════════════════════════════════════════════════════
# SAVE MARKDOWN REPORT
# ═══════════════════════════════════════════════════════════

md_lines = []
md_lines.append("# RUN C2 — Analysis Report")
md_lines.append(f"\n> **Date:** {datetime.now().strftime('%Y-%m-%d %H:%M')}")
md_lines.append(f"> **Run:** C2 — EU Theory-Fixed (Δk=0), perturbations ON")
md_lines.append(f"> **R-1:** {r1_val:.4f}" if r1_val else "> **R-1:** N/A")
md_lines.append(f"> **Samples:** {n_eff_total:.0f} effective ({len(data)} post-burn, {BURN_FRAC*100:.0f}% burn-in)")
md_lines.append("")
md_lines.append("---")
md_lines.append("")
md_lines.append("## Console Output")
md_lines.append("")
md_lines.append("```")
md_lines.extend(output_lines)
md_lines.append("```")

out_md = os.path.join(SCRIPT_DIR, "analysis_RUN_C2.md")
with open(out_md, 'w') as f:
    f.write("\n".join(md_lines))
pr(f"[OK] Markdown saved: {out_md}")
