# RUN A2 — Analysis Report

> **Date:** 2026-05-27 15:11
> **Run:** A2 — EU Sampled (Δk=3), perturbations ON
> **R-1:** 0.0091
> **Samples:** 33888 effective (59885 post-burn)

---

## Console Output

```
======================================================================
RUN A2 — EU PARAMS SAMPLED (Δk=3, perturbations ON)
  ε_IR ∈ [0, 0.15], z_trans ∈ [0.5, 15], b ∈ [0.1, 1.0]
  ide_perturbations=1
======================================================================
  Last checkpoint: N=85545, acc=0.276, R-1=0.0091, R-1_cl=0.060954

### PER-CHAIN STATISTICS
  Chain 1: 10643 raw → 7451 post-burn, Σw=27309, n_eff=4216
  Chain 2: 10741 raw → 7519 post-burn, Σw=27347, n_eff=4244
  Chain 3: 10709 raw → 7497 post-burn, Σw=27130, n_eff=4199
  Chain 4: 10637 raw → 7446 post-burn, Σw=26783, n_eff=4242
  Chain 5: 10560 raw → 7392 post-burn, Σw=27218, n_eff=4177
  Chain 6: 10798 raw → 7559 post-burn, Σw=26670, n_eff=4325
  Chain 7: 10791 raw → 7554 post-burn, Σw=27075, n_eff=4268
  Chain 8: 10666 raw → 7467 post-burn, Σw=27255, n_eff=4227

  TOTAL: 59885 rows post-burn, Σw=216787, n_eff=33888

### PARAMETER POSTERIORS (30% burn-in removed)
Param              Mean        Std        q16        q50        q84     Planck      Δ/σ
----------------------------------------------------------------------------------
ω_b             0.02230    0.00012    0.02218    0.02231    0.02243    0.02237   -0.54σ
ω_cdm           0.11778    0.00062    0.11717    0.11778    0.11839    0.12000   -3.60σ
100θ_s          1.04191    0.00023    1.04168    1.04191    1.04214    1.04092   +4.34σ
τ_reio          0.05866    0.00710    0.05171    0.05845    0.06555    0.05440   +0.60σ
ln(10¹⁰As)      3.04656    0.01429    3.03256    3.04625    3.06048    3.04400   +0.18σ
n_s             0.96763    0.00336    0.96432    0.96766    0.97092    0.96490   +0.81σ
ε_IR            0.07513    0.04329    0.02420    0.07519    0.12627          —        —
z_trans         7.73207    4.18949    2.83691    7.72285   12.69888          —        —
b               0.55126    0.26060    0.24450    0.55114    0.85882          —        —
A_planck        1.00101    0.00248    0.99855    1.00099    1.00347    1.00000   +0.41σ
H₀             68.12572    0.27881   67.84850   68.12530   68.40371   67.36000   +2.75σ
σ₈              0.80562    0.00586    0.79985    0.80549    0.81138    0.81110   -0.94σ
Ω_m             0.30322    0.00358    0.29965    0.30321    0.30677    0.31530   -3.37σ
S₈              0.80992    0.00798    0.80200    0.80984    0.81779    0.83400   -3.02σ
r_drag        147.76173    0.18561  147.57895  147.76051  147.94872  147.09000   +3.62σ
H₀_LKI         71.68960    0.29339   71.39787   71.68916   71.98213          —        —
f_cdm(z=0)      0.92104    0.05330    0.86066    0.92852    0.97880          —        —
I_GKI           0.12594    0.08822    0.03214    0.11125    0.22508          —        —

### EU PARAMETERS — NB01 UV COMPARISON
Param               MCMC A2       ±σ      NB01 UV      Δ/σ    CV%
----------------------------------------------------------------
eu_epsilon_ir       0.07513  0.04329      0.04264    +0.75σ   57.6%
eu_z_trans          7.73207  4.18949      5.98600    +0.42σ   54.2%
eu_b                0.55126  0.26060      0.52778    +0.09σ   47.3%

### PRIOR VOLUME ANALYSIS
  eu_epsilon_ir : prior=[0.0, 0.15], 68%CI=[0.0242, 0.1263], compression=0.68
    ⚠️ PRIOR-DOMINATED (compression > 0.5)
  eu_z_trans    : prior=[0.5, 15.0], 68%CI=[2.8369, 12.6989], compression=0.68
    ⚠️ PRIOR-DOMINATED (compression > 0.5)
  eu_b          : prior=[0.1, 1.0], 68%CI=[0.2445, 0.8588], compression=0.68
    ⚠️ PRIOR-DOMINATED (compression > 0.5)

### DERIVED QUANTITIES — NB01 UV COMPARISON
Param             MCMC A2       ±σ      NB01 UV      Δ/σ
--------------------------------------------------------
H0                68.1257   0.2788      68.9000    -2.78σ
H0_LKI            71.6896   0.2934      72.4900    -2.73σ
fcdm_z0            0.9210   0.0533       0.9558    -0.65σ
I_GKI              0.1259   0.0882       0.0678    +0.66σ

### CORRELATION MATRIX (key params)
               H0   sigma8       S8  Omega_m    ω_cdm     ε_IR      z_t        b
      H0    1.000   -0.078   -0.642   -0.975   -0.878    0.003   -0.005   -0.010
  sigma8   -0.078    1.000    0.803    0.108    0.140    0.001   -0.003   -0.000
      S8   -0.642    0.803    1.000    0.680    0.680   -0.002    0.001    0.006
 Omega_m   -0.975    0.108    0.680    1.000    0.961   -0.005    0.005    0.010
   ω_cdm   -0.878    0.140    0.680    0.961    1.000   -0.006    0.005    0.010
    ε_IR    0.003    0.001   -0.002   -0.005   -0.006    1.000    0.002    0.006
     z_t   -0.005   -0.003    0.001    0.005    0.005    0.002    1.000   -0.004
       b   -0.010   -0.000    0.006    0.010    0.010    0.006   -0.004    1.000

### CHI² BEST-FIT
  Best χ² = 12393.41
    planck_NPIPE_highl_CamSpec.TTTEEE: 10546.51
    planck_2018_lowl.TT: 23.15
    planck_2018_lowl.EE: 396.39
    planck_2018_lensing.clik: 8.92
    bao.desi_dr2: 13.06
    sn.pantheonplus: 1405.38
    eu_derived.EU_Derived: -0.00
    chi2__CMB: 10974.97
    chi2__BAO: 13.06
    chi2__SN: 1405.38

### BEST-FIT PARAMETER VALUES (at minimum χ²)
    ω_b          = 0.022249
    ω_cdm        = 0.117777
    100θ_s       = 1.041920
    τ_reio       = 0.056220
    ln(10¹⁰As)   = 3.053966
    n_s          = 0.965311
    ε_IR         = 0.047107
    z_trans      = 13.962423
    b            = 0.699917
    A_planck     = 1.007204
    H₀           = 68.078393
    σ₈           = 0.808118
    Ω_m          = 0.303492
    S₈           = 0.812807
    r_drag       = 147.824010
    H₀_LKI       = 71.639799
    f_cdm(z=0)   = 0.932211
    I_GKI        = 0.105295

======================================================================
H₀_LKI — NATIVE CHAIN VALUES
======================================================================

  H₀_GKI (global) = 68.13 ± 0.28 km/s/Mpc
  H₀_LKI (local)  = 71.69 ± 0.29 km/s/Mpc
  δH₀_void         = +3.56 km/s/Mpc
  q16=71.40, q50=71.69, q84=71.98
  fcdm(z=0)        = 0.921040 ± 0.053300  (NB01: 0.9558)
  I_GKI            = 0.125943 ± 0.088221  (NB01: 0.06785)

======================================================================
TENSION DIAGNOSTICS — EU Sampled (A2) vs Observations
======================================================================

  --- H₀ TENSION ---
  EU H₀_LKI:     71.69 ± 0.29  ← EU (3 EU params sampled)
  Planck ΛCDM:    67.36 ± 0.54
  SH0ES 2024:     73.17 ± 0.86
  TRGB (Freed.):  69.85 ± 1.75
  JAGB (Freed.):  67.96 ± 2.09

  EU H₀_LKI ↔ SH0ES:     1.63σ  ✅
  EU H₀_LKI ↔ TRGB:      1.04σ  ✅
  EU H₀_LKI ↔ JAGB:      1.77σ  ✅
  EU H₀_GKI ↔ SH0ES:     5.58σ  (global)
  Planck ΛCDM ↔ SH0ES:   5.72σ  ❌

  H₀ tension reduction: 5.7σ → 1.6σ (72%)

  --- S₈ TENSION ---
  EU Run A2:    S₈ = 0.8099 ± 0.0080  (perturbations ON, native)
  Planck ΛCDM:  S₈ = 0.834 ± 0.016
  DES-Y3:       S₈ = 0.776 ± 0.017
  EU↔DES:       1.81σ  (was 2.5σ)
  → Tension reduction: 27%

======================================================================
COMPARISON: A2 (EU sampled) vs C2 (EU fixed)
======================================================================

  Param          C2 (fixed)   A2 (sampled)          Δ      Δ/σ
  ------------------------------------------------------------
  H0            68.8890±0.2810     68.1257±0.2788     -0.7633    2.72σ
  sigma8         0.8275±0.0059      0.8056±0.0059     -0.0219    3.71σ
  S8             0.8112±0.0080      0.8099±0.0080     -0.0013    0.16σ
  Omega_m        0.2883±0.0034      0.3032±0.0036     +0.0149    4.16σ
  omega_cdm      0.1193±0.0006      0.1178±0.0006     -0.0015    2.47σ
  H0_LKI        72.4930±0.2950     71.6896±0.2934     -0.8034    2.72σ
  fcdm_z0        0.9558±0.0000      0.9210±0.0533     -0.0348    0.65σ

  Prior volume effect: if A2 posteriors differ from C2, it is due to the
  enlarged parameter space (flat priors on ε, z_t, b) allowing the MCMC to
  explore regions away from the UV prediction. This is NOT physical — it is
  a Bayesian volume effect. C2 (fixed) is the theory prediction.

======================================================================
PERTURBATION NULL TEST: A2 (ON) vs A1 (OFF)
======================================================================

  Param            A1 (OFF)        A2 (ON)          Δ      Δ/σ
  ------------------------------------------------------------
  H0            68.1320±0.2810     68.1257±0.2788     -0.0063    0.02σ
  sigma8         0.8054±0.0059      0.8056±0.0059     +0.0002    0.04σ
  S8             0.8097±0.0080      0.8099±0.0080     +0.0002    0.03σ
  Omega_m        0.3032±0.0036      0.3032±0.0036     +0.0000    0.01σ
  omega_cdm      0.1178±0.0006      0.1178±0.0006     -0.0000    0.03σ

  Note: A1≈A2 expected — MCMC absorbs perturbation effects by adjusting EU params.

======================================================================
SUMMARY — RUN A2: EU SAMPLED (Δk=3, perturbations ON)
======================================================================
  H₀_GKI (global) = 68.13 ± 0.28 km/s/Mpc
  H₀_LKI (local)  = 71.69 ± 0.29 km/s/Mpc
  σ₈               = 0.8056 ± 0.0059
  S₈               = 0.8099 ± 0.0080
  Ω_m              = 0.3032 ± 0.0036
  ε_IR             = 0.07513 ± 0.04329  (NB01: 0.04264)
  z_trans          = 7.732 ± 4.189  (NB01: 5.986)
  b                = 0.5513 ± 0.2606  (NB01: 0.5278)
  fcdm(z=0)        = 0.9210 ± 0.0533

  H₀ tension vs SH0ES: 1.63σ  (ΛCDM: 5.7σ)
  S₈ tension vs DES:   1.81σ  (ΛCDM: 2.5σ)
  R-1 = 0.0091
  Total effective samples: 33888
  Best χ² = 12393.41

[OK] JSON saved: /Users/alvim/Documents/Paper_I/Paper_1_A_Constituent_Law_of_Cosmic_Evolution/02.1_New_Notebooks/NB05 - MCMC COBAYA/eu_aws_vf1/Outputs/RUN_A2_Final/analyze/analysis_RUN_A2.json
```