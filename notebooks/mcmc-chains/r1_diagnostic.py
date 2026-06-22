"""
Independent R-1 Gelman-Rubin diagnostic for D2 chains.
Tests multiple burn-in fractions to understand the discrepancy
between Cobaya's 0.028 and NB05's 0.006.
"""
import numpy as np
import glob
import os

CHAIN_DIR = r"c:\Users\ricar\Desktop\PAPER 1\Current\Paper_I_A_Dissipation_Principle\01_Notebooks\results\MCMC_Chains\RUN_D2_CP22_42k"

# Load all 8 chains
chain_files = sorted(glob.glob(os.path.join(CHAIN_DIR, "eu_NB05D2.*.txt")))
print(f"Found {len(chain_files)} chain files")

# Read header
with open(chain_files[0]) as f:
    header = f.readline().strip().lstrip('#').split()

# Load raw chains (no burn-in)
chains_full = []
for cf in chain_files:
    data = np.loadtxt(cf)
    chains_full.append(data)
    print(f"  {os.path.basename(cf)}: {len(data)} rows")

total_rows = sum(len(c) for c in chains_full)
print(f"\nTotal samples across all chains: {total_rows}")

# Gelman-Rubin R-1 (weighted)
def gelman_rubin(chains, col_idx, weight_idx=0):
    chain_means, chain_vars, chain_n = [], [], []
    for c in chains:
        w = c[:, weight_idx]
        x = c[:, col_idx]
        n = np.sum(w)
        mean = np.average(x, weights=w)
        var = np.average((x - mean)**2, weights=w)
        chain_means.append(mean)
        chain_vars.append(var)
        chain_n.append(n)
    m = len(chains)
    grand_mean = np.mean(chain_means)
    n_avg = np.mean(chain_n)
    B = n_avg / (m - 1) * sum((mu - grand_mean)**2 for mu in chain_means)
    W = np.mean(chain_vars)
    if W < 1e-30:
        return np.nan
    V = (1 - 1/n_avg) * W + B / n_avg
    return V / W - 1

# Cosmological params to check (skip weights, logpost, chi2)
cosmo_params = ['omega_b', 'omega_cdm', 'theta_s_100', 'tau_reio', 
                'logA', 'n_s', 'A_planck', 'H0', 'sigma8', 'S8', 'Omega_m']

param_indices = {}
for p in cosmo_params:
    if p in header:
        param_indices[p] = header.index(p)

# Test different burn-in fractions
burn_fractions = [0.0, 0.10, 0.20, 0.30, 0.50]

print("\n" + "="*90)
print("R-1 DIAGNOSTIC — varying burn-in fraction")
print("="*90)

for burn_frac in burn_fractions:
    chains_trimmed = []
    for c in chains_full:
        burn = int(burn_frac * len(c))
        chains_trimmed.append(c[burn:])
    
    # Cobaya method: use SECOND HALF of each chain
    # Our method: use chain after burn_frac removal
    
    r1_values = {}
    for pname, pidx in param_indices.items():
        r1 = gelman_rubin(chains_trimmed, pidx)
        r1_values[pname] = r1
    
    r1_max = max(r1_values.values())
    r1_max_param = max(r1_values, key=r1_values.get)
    
    total_kept = sum(len(c) for c in chains_trimmed)
    total_weighted = sum(np.sum(c[:, 0]) for c in chains_trimmed)
    
    print(f"\n--- Burn-in = {burn_frac:.0%} (kept {total_kept}/{total_rows} rows, {total_weighted:.0f} weighted samples) ---")
    print(f"{'Parameter':<18} {'R-1':>10}")
    print("-" * 30)
    for pname in cosmo_params:
        if pname in r1_values:
            r1 = r1_values[pname]
            flag = "OK" if r1 < 0.01 else ("~" if r1 < 0.03 else "X")
            print(f"  {pname:<16} {r1:10.6f}  {flag}")
    print(f"\n  MAX R-1 = {r1_max:.6f} ({r1_max_param})")
    print(f"  Status: {'< 0.01 OK' if r1_max < 0.01 else '< 0.03 ~' if r1_max < 0.03 else '> 0.03 X'}")

# Now replicate EXACTLY what Cobaya does: second half
print("\n" + "="*90)
print("COBAYA METHOD — second half of each chain")
print("="*90)

chains_cobaya = []
for c in chains_full:
    half = len(c) // 2
    chains_cobaya.append(c[half:])

r1_cobaya = {}
for pname, pidx in param_indices.items():
    r1 = gelman_rubin(chains_cobaya, pidx)
    r1_cobaya[pname] = r1

r1_max_cobaya = max(r1_cobaya.values())
r1_max_param_cobaya = max(r1_cobaya, key=r1_cobaya.get)

print(f"{'Parameter':<18} {'R-1':>10}")
print("-" * 30)
for pname in cosmo_params:
    if pname in r1_cobaya:
        r1 = r1_cobaya[pname]
        flag = "OK" if r1 < 0.01 else ("~" if r1 < 0.03 else "X")
        print(f"  {pname:<16} {r1:10.6f}  {flag}")
print(f"\n  MAX R-1 = {r1_max_cobaya:.6f} ({r1_max_param_cobaya})")
print(f"  Cobaya checkpoint value: 0.027986")
print(f"  Our replication: {r1_max_cobaya:.6f}")
print(f"  Match: {'YES' if abs(r1_max_cobaya - 0.027986) < 0.005 else 'NO'}")

# Summary table
print("\n" + "="*90)
print("SUMMARY — R-1_max by burn-in method")
print("="*90)
print(f"{'Method':<30} {'Burn-in':>10} {'R-1_max':>12} {'Worst param':<18}")
print("-" * 72)
for burn_frac in burn_fractions:
    chains_trimmed = [c[int(burn_frac * len(c)):] for c in chains_full]
    r1_all = {p: gelman_rubin(chains_trimmed, idx) for p, idx in param_indices.items()}
    r1_max = max(r1_all.values())
    worst = max(r1_all, key=r1_all.get)
    label = f"Remove first {burn_frac:.0%}"
    print(f"  {label:<28} {burn_frac:>8.0%} {r1_max:12.6f} {worst:<18}")

# Cobaya method
print(f"  {'Cobaya (second half)':<28} {'50%':>8} {r1_max_cobaya:12.6f} {r1_max_param_cobaya:<18}")
print(f"  {'Cobaya checkpoint':<28} {'':>8} {'0.027986':>12}")
