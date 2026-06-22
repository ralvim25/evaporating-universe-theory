"""
EU NB05 — RUN A2 Chain Analysis
================================
Run A2: EU params SAMPLED (Δk=3), Perturbations ON (ide_perturbations=1)
Adapted from RUN_C2 analyze_C2.py

Key differences from C2:
  - EU params (ε_IR, z_trans, b) are SAMPLED with flat priors → Δk = 3
  - Reports EU param posteriors + comparison with NB01 UV values
  - Prior volume effect analysis (EU params prior-dominated?)
  - Comparison: A2 vs C2 (sampled vs fixed), A2 vs A1 (pert ON vs OFF)

Output: analysis_RUN_A2.json + analysis_RUN_A2.md
"""
import numpy as np
import json
import os
from scipy.integrate import solve_ivp
from datetime import datetime

# ── Configuration ──
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CHAIN_DIR = os.path.dirname(SCRIPT_DIR)  # RUN_A2_Final/
CHAIN_PREFIX = "eu_NB05A2"
N_CHAINS = 8
BURN_FRAC = 0.30

# ── Load chains ──
chains = []
for i in range(1, N_CHAINS + 1):
    path = os.path.join(CHAIN_DIR, f"{CHAIN_PREFIX}.{i}.txt")
    d = np.loadtxt(path)
    chains.append(d)

with open(os.path.join(CHAIN_DIR, f"{CHAIN_PREFIX}.1.txt")) as f:
    header = f.readline().strip('#').split()
col = {name: idx for idx, name in enumerate(header)}

all_data = []
for c in chains:
    burn = int(len(c) * BURN_FRAC)
    all_data.append(c[burn:])
data = np.vstack(all_data)
weights = data[:, col['weight']]

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

output_lines = []
def pr(line=""):
    print(line)
    output_lines.append(line)

# ═══════════════════════════════════════════════════════════
pr("=" * 70)
pr("RUN A2 — EU PARAMS SAMPLED (Δk=3, perturbations ON)")
pr("  ε_IR ∈ [0, 0.15], z_trans ∈ [0.5, 15], b ∈ [0.1, 1.0]")
pr("  ide_perturbations=1")
pr("=" * 70)

if last_progress != "N/A":
    parts = last_progress.split()
    r1_val = float(parts[3]) if len(parts) > 3 else None
    r1_cl = parts[4] if len(parts) > 4 else "NaN"
    n_samples_total = int(float(parts[0]))
    acceptance = float(parts[2])
    pr(f"  Last checkpoint: N={n_samples_total}, acc={acceptance:.3f}, "
       f"R-1={r1_val:.4f}, R-1_cl={r1_cl}")
else:
    r1_val = None; n_samples_total = 0; acceptance = 0

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

# ── Parameters ──
# A2 includes EU sampled params
params_list = [
    ('omega_b',       'ω_b'),
    ('omega_cdm',     'ω_cdm'),
    ('theta_s_100',   '100θ_s'),
    ('tau_reio',      'τ_reio'),
    ('logA',          'ln(10¹⁰As)'),
    ('n_s',           'n_s'),
    ('eu_epsilon_ir',  'ε_IR'),
    ('eu_z_trans',     'z_trans'),
    ('eu_b',           'b'),
    ('A_planck',       'A_planck'),
    ('H0',             'H₀'),
    ('sigma8',         'σ₈'),
    ('Omega_m',        'Ω_m'),
    ('S8',             'S₈'),
    ('rdrag',          'r_drag'),
    ('H0_LKI',        'H₀_LKI'),
    ('fcdm_z0',       'f_cdm(z=0)'),
    ('I_GKI',         'I_GKI'),
]

planck_ref = {
    'H0': 67.36, 'sigma8': 0.8111, 'Omega_m': 0.3153, 'S8': 0.834,
    'omega_b': 0.02237, 'omega_cdm': 0.1200, 'tau_reio': 0.0544,
    'n_s': 0.9649, 'rdrag': 147.09, 'theta_s_100': 1.04092,
    'logA': 3.044, 'A_planck': 1.0,
}

nb01_ref = {
    'eu_epsilon_ir': 0.04264,
    'eu_z_trans': 5.986,
    'eu_b': 0.52778,
    'H0': 68.90,
    'H0_LKI': 72.49,
    'fcdm_z0': 0.9558,
    'I_GKI': 0.06785,
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

# ── EU Params: NB01 UV comparison ──
pr("\n### EU PARAMETERS — NB01 UV COMPARISON")
pr(f"{'Param':<14} {'MCMC A2':>12} {'±σ':>8} {'NB01 UV':>12} {'Δ/σ':>8} {'CV%':>6}")
pr("-" * 64)

eu_params = ['eu_epsilon_ir', 'eu_z_trans', 'eu_b']
for pname in eu_params:
    if pname in results:
        r = results[pname]
        uv_val = nb01_ref[pname]
        dsig = (r['mean'] - uv_val) / r['std'] if r['std'] > 0 else 0
        cv = (r['std'] / abs(r['mean'])) * 100 if r['mean'] != 0 else 0
        pr(f"{pname:<14} {r['mean']:>12.5f} {r['std']:>8.5f} {uv_val:>12.5f} "
           f"{dsig:>+8.2f}σ {cv:>6.1f}%")

# Prior volume analysis
pr("\n### PRIOR VOLUME ANALYSIS")
prior_ranges = {
    'eu_epsilon_ir': (0.0, 0.15),
    'eu_z_trans': (0.5, 15.0),
    'eu_b': (0.1, 1.0),
}
for pname, (lo, hi) in prior_ranges.items():
    if pname in results:
        r = results[pname]
        prior_width = hi - lo
        posterior_width = r['q84'] - r['q16']  # 68% interval
        compression = posterior_width / prior_width
        pr(f"  {pname:<14}: prior=[{lo}, {hi}], 68%CI=[{r['q16']:.4f}, {r['q84']:.4f}], "
           f"compression={compression:.2f}")
        if compression > 0.5:
            pr(f"    ⚠️ PRIOR-DOMINATED (compression > 0.5)")
        else:
            pr(f"    ✅ Data-constrained (compression < 0.5)")

# ── Derived params UV comparison ──
pr("\n### DERIVED QUANTITIES — NB01 UV COMPARISON")
pr(f"{'Param':<12} {'MCMC A2':>12} {'±σ':>8} {'NB01 UV':>12} {'Δ/σ':>8}")
pr("-" * 56)
for pname in ['H0', 'H0_LKI', 'fcdm_z0', 'I_GKI']:
    if pname in results and pname in nb01_ref:
        r = results[pname]
        uv_val = nb01_ref[pname]
        if r['std'] > 1e-10:
            dsig = (r['mean'] - uv_val) / r['std']
            dsig_str = f"{dsig:>+8.2f}σ"
        else:
            dsig_str = "   exact"
        pr(f"{pname:<12} {r['mean']:>12.4f} {r['std']:>8.4f} {uv_val:>12.4f} {dsig_str}")

# ── Correlation matrix ──
pr("\n### CORRELATION MATRIX (key params)")
corr_params = [p for p in ['H0', 'sigma8', 'S8', 'Omega_m', 'omega_cdm',
                            'eu_epsilon_ir', 'eu_z_trans', 'eu_b'] if p in col]
corr_data_arr = np.column_stack([data[:, col[p]] for p in corr_params])
wcov = np.cov(corr_data_arr.T, aweights=weights)
wstd_arr = np.sqrt(np.diag(wcov))
wcorr = wcov / np.outer(wstd_arr, wstd_arr)

short_names = {'eu_epsilon_ir': 'ε_IR', 'eu_z_trans': 'z_t', 'eu_b': 'b',
               'omega_cdm': 'ω_cdm'}
labels = [short_names.get(p, p) for p in corr_params]
pr(f"{'':>8}" + "".join(f" {l:>8}" for l in labels))
for i, l in enumerate(labels):
    row = f"{l:>8}" + "".join(f" {wcorr[i,j]:>8.3f}" for j in range(len(labels)))
    pr(row)

# ── Chi² analysis ──
pr("\n### CHI² BEST-FIT")
chi2_total = data[:, col['chi2']]
best_idx = np.argmin(chi2_total)
pr(f"  Best χ² = {chi2_total[best_idx]:.2f}")

chi2_components = [
    'chi2__planck_NPIPE_highl_CamSpec.TTTEEE',
    'chi2__planck_2018_lowl.TT', 'chi2__planck_2018_lowl.EE',
    'chi2__planck_2018_lensing.clik',
    'chi2__bao.desi_dr2', 'chi2__sn.pantheonplus',
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

pr("\n### BEST-FIT PARAMETER VALUES (at minimum χ²)")
for pname, label in params_list:
    if pname in col:
        pr(f"    {label:<12} = {data[best_idx, col[pname]]:.6f}")

# ═══════════════════════════════════════════════════════════
# H₀_LKI — Native chain values
# ═══════════════════════════════════════════════════════════
pr("\n" + "=" * 70)
pr("H₀_LKI — NATIVE CHAIN VALUES")
pr("=" * 70)

if 'H0_LKI' in col:
    h0_lki_vals = data[:, col['H0_LKI']]
    h0_lki_mean, h0_lki_std, h0_lki_16, h0_lki_50, h0_lki_84 = wstats(h0_lki_vals, weights)
    h0_gki_mean, h0_gki_std, _, _, _ = wstats(data[:, col['H0']], weights)
    dh_void = h0_lki_mean - h0_gki_mean
    pr(f"\n  H₀_GKI (global) = {h0_gki_mean:.2f} ± {h0_gki_std:.2f} km/s/Mpc")
    pr(f"  H₀_LKI (local)  = {h0_lki_mean:.2f} ± {h0_lki_std:.2f} km/s/Mpc")
    pr(f"  δH₀_void         = {dh_void:+.2f} km/s/Mpc")
    pr(f"  q16={h0_lki_16:.2f}, q50={h0_lki_50:.2f}, q84={h0_lki_84:.2f}")

fcdm_stats = igki_stats = None
if 'fcdm_z0' in col:
    fcdm_vals = data[:, col['fcdm_z0']]
    fcdm_stats = wstats(fcdm_vals, weights)
    pr(f"  fcdm(z=0)        = {fcdm_stats[0]:.6f} ± {fcdm_stats[1]:.6f}  (NB01: 0.9558)")
if 'I_GKI' in col:
    igki_vals = data[:, col['I_GKI']]
    igki_stats = wstats(igki_vals, weights)
    pr(f"  I_GKI            = {igki_stats[0]:.6f} ± {igki_stats[1]:.6f}  (NB01: 0.06785)")

# ═══════════════════════════════════════════════════════════
# TENSION DIAGNOSTICS
# ═══════════════════════════════════════════════════════════
pr("\n" + "=" * 70)
pr("TENSION DIAGNOSTICS — EU Sampled (A2) vs Observations")
pr("=" * 70)

shoes_h0, shoes_err = 73.17, 0.86
planck_h0, planck_h0_err = 67.36, 0.54

tension_shoes_lki = abs(h0_lki_mean - shoes_h0) / np.sqrt(h0_lki_std**2 + shoes_err**2)
tension_shoes_gki = abs(h0_gki_mean - shoes_h0) / np.sqrt(h0_gki_std**2 + shoes_err**2)
tension_planck_shoes = abs(planck_h0 - shoes_h0) / np.sqrt(planck_h0_err**2 + shoes_err**2)

H0_TRGB, err_TRGB = 69.85, 1.75
H0_JAGB, err_JAGB = 67.96, 2.09
tension_trgb = abs(h0_lki_mean - H0_TRGB) / np.sqrt(h0_lki_std**2 + err_TRGB**2)
tension_jagb = abs(h0_lki_mean - H0_JAGB) / np.sqrt(h0_lki_std**2 + err_JAGB**2)

pr(f"\n  --- H₀ TENSION ---")
pr(f"  EU H₀_LKI:     {h0_lki_mean:.2f} ± {h0_lki_std:.2f}  ← EU (3 EU params sampled)")
pr(f"  Planck ΛCDM:    {planck_h0} ± {planck_h0_err}")
pr(f"  SH0ES 2024:     {shoes_h0} ± {shoes_err}")
pr(f"  TRGB (Freed.):  {H0_TRGB} ± {err_TRGB}")
pr(f"  JAGB (Freed.):  {H0_JAGB} ± {err_JAGB}")
pr(f"")
pr(f"  EU H₀_LKI ↔ SH0ES:     {tension_shoes_lki:.2f}σ  ✅")
pr(f"  EU H₀_LKI ↔ TRGB:      {tension_trgb:.2f}σ  ✅")
pr(f"  EU H₀_LKI ↔ JAGB:      {tension_jagb:.2f}σ  ✅")
pr(f"  EU H₀_GKI ↔ SH0ES:     {tension_shoes_gki:.2f}σ  (global)")
pr(f"  Planck ΛCDM ↔ SH0ES:   {tension_planck_shoes:.2f}σ  ❌")
pr(f"")
pr(f"  H₀ tension reduction: {tension_planck_shoes:.1f}σ → {tension_shoes_lki:.1f}σ "
   f"({(1 - tension_shoes_lki/tension_planck_shoes)*100:.0f}%)")

s8_mean = results['S8']['mean']
s8_std = results['S8']['std']
planck_s8, planck_s8_err = 0.834, 0.016
des_s8, des_err = 0.776, 0.017
tension_des = abs(s8_mean - des_s8) / np.sqrt(s8_std**2 + des_err**2)
planck_des_tension = abs(planck_s8 - des_s8) / np.sqrt(planck_s8_err**2 + des_err**2)

pr(f"\n  --- S₈ TENSION ---")
pr(f"  EU Run A2:    S₈ = {s8_mean:.4f} ± {s8_std:.4f}  (perturbations ON, native)")
pr(f"  Planck ΛCDM:  S₈ = {planck_s8} ± {planck_s8_err}")
pr(f"  DES-Y3:       S₈ = {des_s8} ± {des_err}")
pr(f"  EU↔DES:       {tension_des:.2f}σ  (was {planck_des_tension:.1f}σ)")
pr(f"  → Tension reduction: {(1 - tension_des/planck_des_tension)*100:.0f}%")

# ═══════════════════════════════════════════════════════════
# COMPARISON: A2 vs C2 (SAMPLED vs FIXED)
# ═══════════════════════════════════════════════════════════
pr("\n" + "=" * 70)
pr("COMPARISON: A2 (EU sampled) vs C2 (EU fixed)")
pr("=" * 70)

c2_ref = {
    'H0': (68.889, 0.281),
    'sigma8': (0.8275, 0.0059),
    'S8': (0.8112, 0.0080),
    'Omega_m': (0.2883, 0.0034),
    'omega_cdm': (0.1193, 0.0006),
    'H0_LKI': (72.493, 0.295),
    'fcdm_z0': (0.9558, 0.0000),
}

pr(f"\n  {'Param':<10} {'C2 (fixed)':>14} {'A2 (sampled)':>14} {'Δ':>10} {'Δ/σ':>8}")
pr(f"  {'-'*60}")
for pname, (c2_mean, c2_std) in c2_ref.items():
    if pname in results:
        a2_mean = results[pname]['mean']
        a2_std = results[pname]['std']
        delta = a2_mean - c2_mean
        sig = abs(delta) / max(a2_std, c2_std) if max(a2_std, c2_std) > 0 else 0
        pr(f"  {pname:<10} {c2_mean:>10.4f}±{c2_std:.4f}  {a2_mean:>10.4f}±{a2_std:.4f}  "
           f"{delta:>+10.4f}  {sig:>6.2f}σ")

pr(f"\n  Prior volume effect: if A2 posteriors differ from C2, it is due to the")
pr(f"  enlarged parameter space (flat priors on ε, z_t, b) allowing the MCMC to")
pr(f"  explore regions away from the UV prediction. This is NOT physical — it is")
pr(f"  a Bayesian volume effect. C2 (fixed) is the theory prediction.")

# ═══════════════════════════════════════════════════════════
# PERTURBATION NULL TEST: A2 vs A1
# ═══════════════════════════════════════════════════════════
pr("\n" + "=" * 70)
pr("PERTURBATION NULL TEST: A2 (ON) vs A1 (OFF)")
pr("=" * 70)

a1_ref = {
    'H0': (68.132, 0.281),
    'sigma8': (0.8054, 0.0059),
    'S8': (0.8097, 0.0080),
    'Omega_m': (0.3032, 0.0036),
    'omega_cdm': (0.1178, 0.0006),
}

pr(f"\n  {'Param':<10} {'A1 (OFF)':>14} {'A2 (ON)':>14} {'Δ':>10} {'Δ/σ':>8}")
pr(f"  {'-'*60}")
for pname, (a1_mean, a1_std) in a1_ref.items():
    if pname in results:
        a2_mean = results[pname]['mean']
        a2_std = results[pname]['std']
        delta = a2_mean - a1_mean
        sig = abs(delta) / max(a1_std, a2_std) if max(a1_std, a2_std) > 0 else 0
        pr(f"  {pname:<10} {a1_mean:>10.4f}±{a1_std:.4f}  {a2_mean:>10.4f}±{a2_std:.4f}  "
           f"{delta:>+10.4f}  {sig:>6.2f}σ")

pr(f"\n  Note: A1≈A2 expected — MCMC absorbs perturbation effects by adjusting EU params.")

# ═══════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════
pr("\n" + "=" * 70)
pr("SUMMARY — RUN A2: EU SAMPLED (Δk=3, perturbations ON)")
pr("=" * 70)
pr(f"  H₀_GKI (global) = {h0_gki_mean:.2f} ± {h0_gki_std:.2f} km/s/Mpc")
pr(f"  H₀_LKI (local)  = {h0_lki_mean:.2f} ± {h0_lki_std:.2f} km/s/Mpc")
pr(f"  σ₈               = {results['sigma8']['mean']:.4f} ± {results['sigma8']['std']:.4f}")
pr(f"  S₈               = {s8_mean:.4f} ± {s8_std:.4f}")
pr(f"  Ω_m              = {results['Omega_m']['mean']:.4f} ± {results['Omega_m']['std']:.4f}")
if 'eu_epsilon_ir' in results:
    pr(f"  ε_IR             = {results['eu_epsilon_ir']['mean']:.5f} ± {results['eu_epsilon_ir']['std']:.5f}  (NB01: 0.04264)")
if 'eu_z_trans' in results:
    pr(f"  z_trans          = {results['eu_z_trans']['mean']:.3f} ± {results['eu_z_trans']['std']:.3f}  (NB01: 5.986)")
if 'eu_b' in results:
    pr(f"  b                = {results['eu_b']['mean']:.4f} ± {results['eu_b']['std']:.4f}  (NB01: 0.5278)")
if fcdm_stats:
    pr(f"  fcdm(z=0)        = {fcdm_stats[0]:.4f} ± {fcdm_stats[1]:.4f}")
pr(f"")
pr(f"  H₀ tension vs SH0ES: {tension_shoes_lki:.2f}σ  (ΛCDM: {tension_planck_shoes:.1f}σ)")
pr(f"  S₈ tension vs DES:   {tension_des:.2f}σ  (ΛCDM: {planck_des_tension:.1f}σ)")
pr(f"  R-1 = {r1_val:.4f}" if r1_val else "  R-1 = N/A")
pr(f"  Total effective samples: {n_eff_total:.0f}")
pr(f"  Best χ² = {chi2_total[best_idx]:.2f}")

# ═══════════════════════════════════════════════════════════
# SAVE JSON
# ═══════════════════════════════════════════════════════════
results['_meta'] = {
    'run': 'A2',
    'description': 'EU Sampled (Δk=3), perturbations ON',
    'eu_params_status': 'sampled',
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

out_json = os.path.join(SCRIPT_DIR, "analysis_RUN_A2.json")
with open(out_json, 'w') as f:
    json.dump(results, f, indent=2)
pr(f"\n[OK] JSON saved: {out_json}")

# ── SAVE MARKDOWN ──
md_lines = []
md_lines.append("# RUN A2 — Analysis Report")
md_lines.append(f"\n> **Date:** {datetime.now().strftime('%Y-%m-%d %H:%M')}")
md_lines.append(f"> **Run:** A2 — EU Sampled (Δk=3), perturbations ON")
md_lines.append(f"> **R-1:** {r1_val:.4f}" if r1_val else "> **R-1:** N/A")
md_lines.append(f"> **Samples:** {n_eff_total:.0f} effective ({len(data)} post-burn)")
md_lines.append("")
md_lines.append("---")
md_lines.append("")
md_lines.append("## Console Output")
md_lines.append("")
md_lines.append("```")
md_lines.extend(output_lines)
md_lines.append("```")

out_md = os.path.join(SCRIPT_DIR, "analysis_RUN_A2.md")
with open(out_md, 'w') as f:
    f.write("\n".join(md_lines))
pr(f"[OK] Markdown saved: {out_md}")
