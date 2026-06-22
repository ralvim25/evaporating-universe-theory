# RUN C2 — Analysis Report

> **Date:** 2026-05-27 15:11
> **Run:** C2 — EU Theory-Fixed (Δk=0), perturbations ON
> **R-1:** 0.0084
> **Samples:** 33150 effective (58270 post-burn, 30% burn-in)

---

## Console Output

```
======================================================================
RUN C2 — EU THEORY-FIXED (Δk=0, perturbations ON)
  ε_IR=0.04264, b=19/36, z_trans=5.986, ide_perturbations=1
======================================================================
  Last checkpoint: N=83236, acc=0.278, R-1=0.0084, R-1_cl=0.061637

### PER-CHAIN STATISTICS
  Chain 1: 10320 raw → 7224 post-burn, Σw=26303, n_eff=4048
  Chain 2: 10272 raw → 7191 post-burn, Σw=26536, n_eff=3920
  Chain 3: 10492 raw → 7345 post-burn, Σw=26337, n_eff=4204
  Chain 4: 10429 raw → 7301 post-burn, Σw=26331, n_eff=4206
  Chain 5: 10565 raw → 7396 post-burn, Σw=26164, n_eff=4191
  Chain 6: 10403 raw → 7283 post-burn, Σw=25905, n_eff=4183
  Chain 7: 10356 raw → 7250 post-burn, Σw=26346, n_eff=4209
  Chain 8: 10399 raw → 7280 post-burn, Σw=26521, n_eff=4215

  TOTAL: 58270 rows post-burn, Σw=210443, n_eff=33150

### PARAMETER POSTERIORS (30% burn-in removed)
Param              Mean        Std        q16        q50        q84     Planck      Δ/σ
----------------------------------------------------------------------------------
ω_b             0.02218    0.00012    0.02206    0.02219    0.02231    0.02237   -1.49σ
ω_cdm           0.11926    0.00063    0.11863    0.11926    0.11989    0.12000   -1.17σ
100θ_s          1.04180    0.00023    1.04158    1.04180    1.04203    1.04092   +3.85σ
τ_reio          0.05578    0.00694    0.04907    0.05562    0.06255    0.05440   +0.20σ
ln(10¹⁰As)      3.04362    0.01397    3.02995    3.04343    3.05725    3.04400   -0.03σ
n_s             0.96371    0.00345    0.96028    0.96373    0.96717    0.96490   -0.34σ
A_planck        1.00081    0.00249    0.99833    1.00085    1.00327    1.00000   +0.33σ
H₀             68.88906    0.28056   68.60417   68.88620   69.16820   67.36000   +5.45σ
σ₈              0.82749    0.00592    0.82171    0.82738    0.83333    0.81110   +2.77σ
Ω_m             0.28829    0.00339    0.28489    0.28828    0.29171    0.31530   -7.97σ
S₈              0.81117    0.00796    0.80318    0.81119    0.81909    0.83400   -2.87σ
r_drag        147.50070    0.18629  147.31625  147.50049  147.68700  147.09000   +2.20σ
H₀_LKI         72.49288    0.29524   72.19308   72.48986   72.78662          —        —
f_cdm(z=0)      0.95577    0.00000    0.95577    0.95577    0.95577          —        —
I_GKI           0.06786    0.00000    0.06786    0.06786    0.06786          —        —

### NB01 UV THEORY COMPARISON (zero free parameters)
Param             MCMC C2       ±σ      NB01 UV      Δ/σ
--------------------------------------------------------
H0                68.8891   0.2806      68.9000    -0.04σ
H0_LKI            72.4929   0.2952      72.4900    +0.01σ
fcdm_z0            0.9558   0.0000       0.9558    exact
I_GKI              0.0679   0.0000       0.0678    exact

### CORRELATION MATRIX (key params)
                  H0    sigma8        S8   Omega_m   omega_b omega_cdm       n_s  tau_reio
        H0     1.000    -0.087    -0.647    -0.974     0.615    -0.880     0.441     0.279
    sigma8    -0.087     1.000     0.804     0.124     0.039     0.163     0.185     0.830
        S8    -0.647     0.804     1.000     0.690    -0.241     0.696    -0.133     0.441
   Omega_m    -0.974     0.124     0.690     1.000    -0.450     0.963    -0.447    -0.275
   omega_b     0.615     0.039    -0.241    -0.450     1.000    -0.273     0.237     0.179
 omega_cdm    -0.880     0.163     0.696     0.963    -0.273     1.000    -0.430    -0.256
       n_s     0.441     0.185    -0.133    -0.447     0.237    -0.430     1.000     0.239
  tau_reio     0.279     0.830     0.441    -0.275     0.179    -0.256     0.239     1.000

### CHI² BEST-FIT
  Best χ² = 12394.77
    planck_NPIPE_highl_CamSpec.TTTEEE: 10544.96
    planck_2018_lowl.TT: 22.92
    planck_2018_lowl.EE: 396.06
    planck_2018_lensing.clik: 9.44
    bao.desi_dr2: 12.54
    sn.pantheonplus: 1408.86
    eu_derived.EU_Derived: -0.00
    chi2__CMB: 10973.37
    chi2__BAO: 12.54
    chi2__SN: 1408.86

### BEST-FIT PARAMETER VALUES (at minimum χ²)
    ω_b          = 0.022246
    ω_cdm        = 0.119019
    100θ_s       = 1.041670
    τ_reio       = 0.054304
    ln(10¹⁰As)   = 3.049748
    n_s          = 0.964615
    A_planck     = 1.005994
    H₀           = 68.982666
    σ₈           = 0.829047
    Ω_m          = 0.287129
    S₈           = 0.811068
    r_drag       = 147.495370
    H₀_LKI       = 72.591377
    f_cdm(z=0)   = 0.955765
    I_GKI        = 0.067864

======================================================================
H₀_LKI — NATIVE CHAIN VALUES (from eu_derived/shoes_lki)
======================================================================

  H₀_GKI (global) = 68.89 ± 0.28 km/s/Mpc
  H₀_LKI (local)  = 72.49 ± 0.30 km/s/Mpc
  δH₀_void         = +3.60 km/s/Mpc
  Boost factor      = 1.0523
  q16=72.19, q50=72.49, q84=72.79
  fcdm(z=0)         = 0.955765 ± 0.000000  (NB01: 0.9558)
  I_GKI             = 0.067864 ± 0.000000  (NB01: 0.06785)

======================================================================
TENSION DIAGNOSTICS — EU Theory-Fixed (C2) vs Observations
======================================================================

  --- H₀ TENSION (LKI = local prediction for distance ladder) ---
  EU H₀_LKI:     72.49 ± 0.30  ← EU (0 free params)
  Planck ΛCDM:    67.36 ± 0.54
  SH0ES 2024:     73.17 ± 0.86
  TRGB (Freed.):  69.85 ± 1.75
  JAGB (Freed.):  67.96 ± 2.09

  EU H₀_LKI ↔ SH0ES:     0.74σ  ✅
  EU H₀_LKI ↔ TRGB:      1.49σ  ✅
  EU H₀_LKI ↔ JAGB:      2.15σ  ✅
  EU H₀_GKI ↔ SH0ES:     4.73σ  (global, for reference)
  Planck ΛCDM ↔ SH0ES:   5.72σ  ❌

  H₀ tension reduction: 5.7σ → 0.7σ (87%)

  --- S₈ TENSION ---
  EU Run C2:    S₈ = 0.8112 ± 0.0080  (perturbations ON, native)
  Planck ΛCDM:  S₈ = 0.834 ± 0.016
  DES-Y3:       S₈ = 0.776 ± 0.017
  EU↔DES:       1.87σ  (was 2.5σ)
  → Tension reduction: 25%

======================================================================
PERTURBATION NULL TEST: C2 (ON) vs C1 (OFF)
======================================================================

  Param            C1 (OFF)        C2 (ON)          Δ      Δ/σ
  ------------------------------------------------------------
  H0            68.8870±0.2810     68.8891±0.2806     +0.0021    0.01σ
  sigma8         0.8276±0.0059      0.8275±0.0059     -0.0001    0.02σ
  S8             0.8113±0.0080      0.8112±0.0080     -0.0001    0.02σ
  Omega_m        0.2883±0.0034      0.2883±0.0034     -0.0000    0.00σ
  omega_cdm      0.1193±0.0006      0.1193±0.0006     -0.0000    0.02σ
  H0_LKI        72.4200±0.2800     72.4929±0.2952     +0.0729    0.25σ

  Expected: All Δ/σ < 0.5 (perturbations cancel for w=-1, Valiviita+2008)

======================================================================
SUMMARY — RUN C2: EU THEORY-FIXED (Δk=0, perturbations ON)
======================================================================
  H₀_GKI (global) = 68.89 ± 0.28 km/s/Mpc
  H₀_LKI (local)  = 72.49 ± 0.30 km/s/Mpc  ← EU prediction
  σ₈               = 0.8275 ± 0.0059
  S₈               = 0.8112 ± 0.0080
  Ω_m              = 0.2883 ± 0.0034
  fcdm(z=0)        = 0.955765 ± 0.000000

  H₀ tension vs SH0ES: 0.74σ  (ΛCDM: 5.7σ)
  S₈ tension vs DES:   1.87σ  (ΛCDM: 2.5σ)
  R-1 = 0.0084
  Total effective samples: 33150
  Perturbations: ON (ide_perturbations=1)
  Best χ² = 12394.77

[OK] JSON saved: /Users/alvim/Documents/Paper_I/Paper_1_A_Constituent_Law_of_Cosmic_Evolution/02.1_New_Notebooks/NB05 - MCMC COBAYA/eu_aws_vf1/Outputs/RUN_C2_final/analyze/analysis_RUN_C2.json
```