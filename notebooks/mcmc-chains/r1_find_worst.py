"""Quick check: which parameter gives R-1 = 0.028 in Cobaya?"""
import numpy as np, glob, os

CHAIN_DIR = r"c:\Users\ricar\Desktop\PAPER 1\Current\Paper_I_A_Dissipation_Principle\01_Notebooks\results\MCMC_Chains\RUN_D2_CP22_42k"
chain_files = sorted(glob.glob(os.path.join(CHAIN_DIR, "eu_NB05D2.*.txt")))

with open(chain_files[0]) as f:
    header = f.readline().strip().lstrip('#').split()

chains_full = [np.loadtxt(cf) for cf in chain_files]

# Cobaya method: second half
chains_half = [c[len(c)//2:] for c in chains_full]

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

# ALL parameters (including DES nuisance)
skip = {'weight', 'minuslogpost', 'minuslogprior', 'minuslogprior__0', 'chi2'}
print(f"{'Parameter':<45} {'R-1 (50%)':>12} {'R-1 (30%)':>12} {'R-1 (0%)':>12}")
print("-" * 85)

chains_30 = [c[int(0.3*len(c)):] for c in chains_full]

all_r1 = []
for j, name in enumerate(header):
    if name in skip or name.startswith('chi2__'):
        continue
    r1_50 = gelman_rubin(chains_half, j)
    r1_30 = gelman_rubin(chains_30, j)
    r1_00 = gelman_rubin(chains_full, j)
    if np.isnan(r1_50):
        continue
    all_r1.append((name, r1_50, r1_30, r1_00))

# Sort by R-1 (50%) descending
all_r1.sort(key=lambda x: x[1], reverse=True)

for name, r1_50, r1_30, r1_00 in all_r1:
    marker = " <-- WORST" if r1_50 > 0.02 else ""
    print(f"  {name:<43} {r1_50:12.6f} {r1_30:12.6f} {r1_00:12.6f}{marker}")

print(f"\nR-1 max (50% burn): {max(x[1] for x in all_r1):.6f} = {max(all_r1, key=lambda x: x[1])[0]}")
print(f"R-1 max (30% burn): {max(x[2] for x in all_r1):.6f} = {max(all_r1, key=lambda x: x[2])[0]}")
print(f"R-1 max (0% burn):  {max(x[3] for x in all_r1):.6f} = {max(all_r1, key=lambda x: x[3])[0]}")
print(f"Cobaya checkpoint:  0.027986")
