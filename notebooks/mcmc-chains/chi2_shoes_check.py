"""Calculate mean chi2_SH0ES from D2 chains too."""
import numpy as np, glob, os

D2_DIR = r"c:\Users\ricar\Desktop\PAPER 1\Current\Paper_I_A_Dissipation_Principle\01_Notebooks\results\MCMC_Chains\RUN_D2_CP22_42k"
chain_files = sorted(glob.glob(os.path.join(D2_DIR, "eu_NB05D2.*.txt")))

with open(chain_files[0]) as f:
    header = f.readline().strip().lstrip('#').split()

shoes_cols = [i for i, h in enumerate(header) if 'shoes' in h.lower()]
print(f"D2 SH0ES columns: {[(i, header[i]) for i in shoes_cols]}")

if shoes_cols:
    col = shoes_cols[0]
    all_chi2 = []
    all_weights = []
    for cf in chain_files:
        data = np.loadtxt(cf)
        burn = int(0.3 * len(data))
        data = data[burn:]
        all_chi2.extend(data[:, col])
        all_weights.extend(data[:, 0])
    
    mean_chi2 = np.average(all_chi2, weights=all_weights)
    print(f"D2 mean chi2_SH0ES (30% burn-in, weighted): {mean_chi2:.4f}")
    print(f"D2 bestfit chi2_SH0ES (from JSON): 1.3200")

# Also compute C2 mean for shoes (C2 doesnt have shoes in likelihood, 
# but we can compute it from H0_LKI column)
C2_DIR = r"c:\Users\ricar\Desktop\PAPER 1\Current\Paper_I_A_Dissipation_Principle\01_Notebooks\results\MCMC_Chains\RUN_C2_Final"
if os.path.exists(C2_DIR):
    c2_files = sorted(glob.glob(os.path.join(C2_DIR, "*.txt")))
    with open(c2_files[0]) as f:
        c2_header = f.readline().strip().lstrip('#').split()
    
    h0lki_idx = [i for i, h in enumerate(c2_header) if 'H0_LKI' in h]
    print(f"\nC2 H0_LKI column: {[(i, c2_header[i]) for i in h0lki_idx]}")
    
    if h0lki_idx:
        col = h0lki_idx[0]
        all_h0lki = []
        all_w = []
        for cf in c2_files:
            data = np.loadtxt(cf)
            burn = int(0.3 * len(data))
            data = data[burn:]
            all_h0lki.extend(data[:, col])
            all_w.extend(data[:, 0])
        
        h0_shoes = 73.17
        err_shoes = 0.86
        # chi2 = ((H0_LKI - H0_shoes) / err_shoes)^2
        chi2_vals = [((h - h0_shoes) / err_shoes)**2 for h in all_h0lki]
        mean_chi2_c2 = np.average(chi2_vals, weights=all_w)
        mean_h0lki = np.average(all_h0lki, weights=all_w)
        print(f"C2 mean H0_LKI: {mean_h0lki:.4f}")
        print(f"C2 mean chi2_SH0ES (computed from H0_LKI): {mean_chi2_c2:.4f}")
        print(f"C2 bestfit chi2_SH0ES (from NB06): 0.2922")
