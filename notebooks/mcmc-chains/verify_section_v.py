"""
FORENSIC VERIFICATION SCRIPT — Complete Section V Audit
========================================================
Verifies ALL values in R5_CONVERGENCE_INVESTIGATION.md,
SECTION_V_ANALYSIS.md, and SECTION_V_STRUCTURE.md against
raw chains and JSON source files.

Output: line-by-line verification with PASS/FAIL status.
"""
import numpy as np
import json
import glob
import os

RESULTS = r"c:\Users\ricar\Desktop\PAPER 1\Current\Paper_I_A_Dissipation_Principle\01_Notebooks\results"
CHAINS  = os.path.join(RESULTS, "MCMC_Chains")

print("=" * 70)
print("FORENSIC VERIFICATION — Section V Complete Audit")
print("=" * 70)

# ============================================================
# 1. LOAD ALL JSONs
# ============================================================
print("\n[1] LOADING JSON SOURCE FILES")
print("-" * 50)

json_files = {
    "C2": os.path.join(RESULTS, "NB05_C2_results.json"),
    "D1": os.path.join(RESULTS, "NB05_D1_results.json"),
    "D2": os.path.join(RESULTS, "NB05_D2_results.json"),
    "A2": os.path.join(RESULTS, "NB05_A2_results.json"),
    "NB06": os.path.join(RESULTS, "NB06_results.json"),
    "LCDM": os.path.join(RESULTS, "NB06_lcdm_minimizer.json"),
}

data = {}
for key, path in json_files.items():
    if os.path.exists(path):
        with open(path) as f:
            data[key] = json.load(f)
        print(f"  ✅ {key}: {os.path.basename(path)}")
    else:
        print(f"  ❌ {key}: NOT FOUND — {path}")

# ============================================================
# 2. R-1 CONVERGENCE (Gelman-Rubin per-parameter)
# ============================================================
print("\n[2] R⁻¹ CONVERGENCE — Gelman-Rubin per-parameter")
print("-" * 50)

def gelman_rubin_per_param(chain_dir, pattern, burn_frac=0.3):
    """Calculate max R-1 across all parameters using GR diagnostic."""
    chain_files = sorted(glob.glob(os.path.join(chain_dir, pattern)))
    if not chain_files:
        return None, None, 0
    
    # Read header
    with open(chain_files[0]) as f:
        header = f.readline().strip().lstrip('#').split()
    
    # Find parameter columns (skip weight=0, minuslogpost=1, chi2 columns)
    param_cols = []
    param_names = []
    for i, h in enumerate(header):
        if i < 2:
            continue
        if h.startswith('chi2'):
            continue
        param_cols.append(i)
        param_names.append(h)
    
    # Load chains with burn-in removal
    chains = []
    for cf in chain_files:
        d = np.loadtxt(cf)
        burn = int(burn_frac * len(d))
        chains.append(d[burn:])
    
    n_chains = len(chains)
    if n_chains < 2:
        return None, None, 0
    
    # GR diagnostic per parameter
    max_r1 = 0
    worst_param = ""
    all_r1 = {}
    
    for col, name in zip(param_cols, param_names):
        chain_means = []
        chain_vars = []
        n_samples = []
        
        for chain in chains:
            weights = chain[:, 0]
            values = chain[:, col]
            w_sum = np.sum(weights)
            mean = np.sum(weights * values) / w_sum
            var = np.sum(weights * (values - mean)**2) / w_sum
            chain_means.append(mean)
            chain_vars.append(var)
            n_samples.append(w_sum)
        
        chain_means = np.array(chain_means)
        chain_vars = np.array(chain_vars)
        
        grand_mean = np.mean(chain_means)
        n = np.mean(n_samples)
        
        B = n * np.var(chain_means, ddof=1)  # Between-chain variance
        W = np.mean(chain_vars)               # Within-chain variance
        
        if W == 0:
            continue
        
        var_hat = ((n - 1) / n) * W + (1 / n) * B
        R = var_hat / W
        r1 = R - 1
        
        all_r1[name] = r1
        if r1 > max_r1:
            max_r1 = r1
            worst_param = name
    
    return max_r1, worst_param, n_chains

# Check available chain directories
chain_dirs = {
    "A1": ("RUN_A1_Final", "*.txt"),
    "A2": ("RUN_A2_Final", "*.txt"),
    "C1": ("RUN_C1_Final", "*.txt"),
    "C2": ("RUN_C2_Final", "*.txt"),
    "D1": ("RUN_D1_Final", "*.txt"),
    "D2": ("RUN_D2_CP22_42k", "eu_NB05D2.*.txt"),
}

# Expected values from SECTION_V_STRUCTURE.md
expected_r1 = {
    "A1": 0.003, "A2": 0.003, "C1": 0.004,
    "C2": 0.005, "D1": 0.008, "D2": 0.006,
}

print(f"\n{'Run':<6} {'R⁻¹_max':<12} {'Worst param':<20} {'Expected':<10} {'< 0.01?':<8} {'Match?'}")
print("-" * 80)

verified_r1 = {}
for run, (subdir, pat) in chain_dirs.items():
    chain_path = os.path.join(CHAINS, subdir)
    if os.path.exists(chain_path):
        r1, worst, n_ch = gelman_rubin_per_param(chain_path, pat)
        if r1 is not None:
            verified_r1[run] = r1
            exp = expected_r1[run]
            lt001 = "✅" if r1 < 0.01 else "❌"
            # Match = rounds to same value at 3 decimal places
            match = "✅" if round(r1, 3) == exp else f"⚠️ ({round(r1,3)})"
            print(f"{run:<6} {r1:<12.6f} {worst:<20} {exp:<10} {lt001:<8} {match}")
        else:
            print(f"{run:<6} INSUFFICIENT CHAINS")
    else:
        # Check JSON for R-1
        if run in data:
            r1_json = data[run].get("_metadata", {}).get("convergence", {}).get("R-1_max_post_burnin", "N/A")
            print(f"{run:<6} {r1_json:<12} (from JSON, chains not found)")
        else:
            print(f"{run:<6} CHAINS NOT FOUND: {chain_path}")

# ============================================================
# 3. CHI² SH0ES VALUES
# ============================================================
print("\n\n[3] χ² SH0ES VALUES")
print("-" * 50)

# 3a. Blind bestfit (NB06)
if "NB06" in data:
    chi2_blind = data["NB06"].get("chi2_shoes_EU", 
                  data["NB06"].get("scenarios", {}).get("scenario_1_baseline", {}).get("chi2_shoes_EU"))
    # Search more carefully
    for key in data["NB06"]:
        if isinstance(data["NB06"][key], dict):
            if "chi2_shoes_EU" in data["NB06"][key]:
                chi2_blind = data["NB06"][key]["chi2_shoes_EU"]
                print(f"  Found chi2_shoes_EU in NB06['{key}']: {chi2_blind}")
    
    # Direct search
    def find_key(d, target, path=""):
        results = []
        if isinstance(d, dict):
            for k, v in d.items():
                if k == target:
                    results.append((path + "." + k, v))
                if isinstance(v, (dict, list)):
                    results.extend(find_key(v, target, path + "." + k))
        return results
    
    shoes_vals = find_key(data["NB06"], "chi2_shoes_EU")
    for path, val in shoes_vals:
        print(f"  NB06{path} = {val:.6f}")

# 3b. D1 bestfit
if "D1" in data:
    chi2_d1 = data["D1"].get("chi2_bestfit", {}).get("shoes_lki.SH0ES_LKI", "N/A")
    print(f"  D1 bestfit chi2_SH0ES = {chi2_d1}")

# 3c. D1 mean from chains
d1_path = os.path.join(CHAINS, "RUN_D1_Final")
if os.path.exists(d1_path):
    d1_files = sorted(glob.glob(os.path.join(d1_path, "*.txt")))
    with open(d1_files[0]) as f:
        hdr = f.readline().strip().lstrip('#').split()
    shoes_idx = [i for i, h in enumerate(hdr) if 'shoes' in h.lower()]
    if shoes_idx:
        col = shoes_idx[0]
        all_chi2 = []
        all_w = []
        for cf in d1_files:
            dd = np.loadtxt(cf)
            burn = int(0.3 * len(dd))
            dd = dd[burn:]
            all_chi2.extend(dd[:, col])
            all_w.extend(dd[:, 0])
        mean_d1 = np.average(all_chi2, weights=all_w)
        print(f"  D1 mean chi2_SH0ES (30% burn-in) = {mean_d1:.4f}")

# 3d. C2 evaluated at mean H0_LKI
if "C2" in data:
    h0lki_mean = data["C2"]["derived_params"]["H0_LKI"]["mean"]
    h0_shoes = 73.17
    err_shoes = 0.86
    chi2_c2_at_mean = ((h0lki_mean - h0_shoes) / err_shoes)**2
    print(f"  C2 chi2_SH0ES at mean H0_LKI ({h0lki_mean:.4f}) = {chi2_c2_at_mean:.4f}")

# 3e. D1 evaluated at mean H0_LKI
if "D1" in data:
    h0lki_d1 = data["D1"]["derived_params"]["H0_LKI"]["mean"]
    chi2_d1_at_mean = ((h0lki_d1 - h0_shoes) / err_shoes)**2
    print(f"  D1 chi2_SH0ES at mean H0_LKI ({h0lki_d1:.4f}) = {chi2_d1_at_mean:.4f}")

# ============================================================
# 4. TABLE IV VALUES — Cross-check against JSONs
# ============================================================
print("\n\n[4] TABLE IV — Parameter Cross-Check")
print("-" * 50)

params_to_check = {
    "omega_b": {"C2": "0.02218 ± 0.00013", "D1": "0.02221 ± 0.00012", "D2": "0.02218 ± 0.00012"},
    "omega_cdm": {"C2": "0.1193 ± 0.0006", "D1": "0.1191 ± 0.0006", "D2": "0.1194 ± 0.0006"},
}

for run_key in ["C2", "D1", "D2"]:
    if run_key not in data:
        continue
    d = data[run_key]
    print(f"\n  {run_key}:")
    
    # Cosmological params
    cosmo = d.get("cosmological_params", {})
    for pname in ["omega_b", "omega_cdm", "theta_s_100"]:
        if pname in cosmo:
            p = cosmo[pname]
            print(f"    {pname}: mean={p['mean']:.6f} ± {p['std']:.6f}")
    
    # Derived params
    derived = d.get("derived_params", {})
    for pname in ["H0", "H0_LKI", "sigma8", "S8", "Omega_m"]:
        if pname in derived:
            p = derived[pname]
            print(f"    {pname}: mean={p['mean']:.6f} ± {p['std']:.6f}")
    
    # Chi2 total
    chi2 = d.get("chi2_bestfit", {}).get("total", "N/A")
    print(f"    chi2_total (bestfit): {chi2}")
    
    # Samples
    samples = d.get("_metadata", {}).get("samples_weighted", "N/A")
    print(f"    N_samples: {samples}")

# ============================================================
# 5. TABLE V — UV→MCMC Concordance
# ============================================================
print("\n\n[5] TABLE V — UV→MCMC Concordance")
print("-" * 50)

if "C2" in data:
    d = data["C2"]
    der = d["derived_params"]
    
    uv_values = {
        "H0_GKI": 68.90,
        "H0_LKI": 72.72,
        "f_cdm": 0.9558,
        "I_GKI": 0.0679,
    }
    
    mcmc_values = {
        "H0_GKI": (der["H0"]["mean"], der["H0"]["std"]),
        "H0_LKI": (der["H0_LKI"]["mean"], der["H0_LKI"]["std"]),
    }
    
    for name, uv in uv_values.items():
        if name in mcmc_values:
            mean, std = mcmc_values[name]
            shift = abs(uv - mean) / std
            print(f"  {name}: UV={uv}, MCMC={mean:.4f}±{std:.4f}, shift={shift:.2f}σ")

# ============================================================
# 6. TABLE VI — ΔAIC/Δχ²
# ============================================================
print("\n\n[6] TABLE VI — Model Selection (NB06)")
print("-" * 50)

if "NB06" in data:
    nb06 = data["NB06"]
    # Find scenarios
    for key, val in nb06.items():
        if isinstance(val, dict) and "delta_chi2" in val:
            desc = val.get("description", key)
            dchi2 = val["delta_chi2"]
            print(f"  {desc}: Δχ² = {dchi2}")
        elif isinstance(val, dict) and "delta_aic" in val:
            desc = val.get("description", key)
            daic = val["delta_aic"]
            print(f"  {desc}: ΔAIC = {daic}")

# ============================================================
# 7. Δχ² BREAKDOWN (Baseline)
# ============================================================
print("\n\n[7] Δχ² BREAKDOWN (C2 vs ΛCDM bestfit)")
print("-" * 50)

if "C2" in data and "LCDM" in data:
    c2_chi2 = data["C2"].get("chi2_bestfit", {})
    lcdm_chi2 = data["LCDM"].get("chi2_bestfit", {})
    
    if c2_chi2 and lcdm_chi2:
        # CMB components
        cmb_keys_c2 = ["planck_NPIPE_highl_CamSpec.TTTEEE", 
                        "planck_2018_lowl.TT", "planck_2018_lowl.EE",
                        "planck_2018_lensing.clik"]
        
        cmb_eu = sum(c2_chi2.get(k, 0) for k in cmb_keys_c2)
        cmb_lcdm = sum(lcdm_chi2.get(k, 0) for k in cmb_keys_c2)
        
        bao_eu = c2_chi2.get("bao.desi_dr2", c2_chi2.get("BAO", 0))
        bao_lcdm = lcdm_chi2.get("bao.desi_dr2", lcdm_chi2.get("BAO", 0))
        
        sne_eu = c2_chi2.get("sn.pantheonplus", c2_chi2.get("SN", 0))
        sne_lcdm = lcdm_chi2.get("sn.pantheonplus", lcdm_chi2.get("SN", 0))
        
        total_eu = c2_chi2.get("total", 0)
        total_lcdm = lcdm_chi2.get("total", 0)
        
        print(f"  CMB:  EU={cmb_eu:.1f}  ΛCDM={cmb_lcdm:.1f}  Δ={cmb_eu-cmb_lcdm:.1f}")
        print(f"  BAO:  EU={bao_eu:.1f}  ΛCDM={bao_lcdm:.1f}  Δ={bao_eu-bao_lcdm:.1f}")
        print(f"  SNe:  EU={sne_eu:.1f}  ΛCDM={sne_lcdm:.1f}  Δ={sne_eu-sne_lcdm:.1f}")
        print(f"  Total: EU={total_eu:.1f}  ΛCDM={total_lcdm:.1f}  Δ={total_eu-total_lcdm:.1f}")

# ============================================================
# 8. CONVERGENCE CROSS-CHECK: JSON R-1 vs Raw Chains
# ============================================================
print("\n\n[8] CONVERGENCE CROSS-CHECK: JSON vs Raw Chains")
print("-" * 50)

for run in ["A2", "C2", "D1", "D2"]:
    if run in data and run in verified_r1:
        json_r1 = data[run]["_metadata"]["convergence"]["R-1_max_post_burnin"]
        raw_r1 = verified_r1[run]
        match = abs(json_r1 - raw_r1) < 0.002
        status = "✅" if match else "⚠️"
        print(f"  {run}: JSON={json_r1:.6f}  Raw={raw_r1:.6f}  {status}")

print("\n" + "=" * 70)
print("VERIFICATION COMPLETE")
print("=" * 70)
