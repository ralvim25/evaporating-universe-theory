# MCMC Dataset Separation: Technical Notes for Referees

## Overview

This document describes the methodology used in the MCMC analysis for Paper I dataset separation tests. All implementations follow standard practices from the cosmological literature.

---

## 1. Type Ia Supernovae (SNe) Analysis

### Data Source
- **Pantheon+ 2022** (Brout et al., arXiv:2202.04077)
- 16 binned distance modulus measurements
- Redshift range: 0.023 - 1.7
- Low-z calibrated with SH0ES Cepheids

### Methodology: Analytic Marginalization over M

SNe Ia measure **relative distances** via the distance modulus:
```
μ = m_B - M
```
where `m_B` is the apparent magnitude and `M` is the absolute magnitude.

Since `M` is degenerate with `H₀`, we use **analytic marginalization** (Conley et al. 2011, Betoule et al. 2014):

```python
# Marginalized chi-squared (flat prior on M)
Delta = μ_obs - μ_pred
inv_var = 1 / σ²
A = Σ(inv_var)
B = Σ(Delta × inv_var)  
C = Σ(Delta² × inv_var)
χ²_marg = C - B²/A
```

This is mathematically equivalent to fitting `M` and profiling over it, but faster.

### Free Parameters
| Parameter | Prior | Description |
|-----------|-------|-------------|
| H₀ | [60, 85] km/s/Mpc | Hubble constant |
| Ωm | [0.15, 0.45] | Matter density |
| w₀ | [-1.5, -0.8] | DE equation of state (today) |
| z_trans | [0.05, 0.5] | Transition redshift |

### References
- Brout et al. 2022, ApJ 938, 110 (Pantheon+ Cosmology)
- Conley et al. 2011, ApJS 192, 1 (SNLS3)
- Betoule et al. 2014, A&A 568, A22 (JLA)

---

## 2. CMB Distance Priors Analysis

### Data Source
- **Planck 2018** (arXiv:1807.06209)  
- Compressed distance priors: R, l_a, ωb

### Observables

**Shift Parameter R** (Efstathiou & Bond 1999):
```
R = √Ωm × D_M(z*) / D_H
```
where `D_M` = comoving distance, `D_H = c/H₀`, `z* = 1089.92`

**Acoustic Scale l_a** (Chen et al. 2017):
```
l_a = π × (1+z*) × D_A(z*) / r_s(z*)
```
where `r_s` = sound horizon at recombination

### Sound Horizon
Using Planck 2018 fitting formula:
```
r_s = 55.154 × exp[-72.3(ωᵥ + 0.0006)²] / (ωm^0.25351 × ωb^0.12807)
```

### Planck 2018 Values
| Observable | Value | σ |
|------------|-------|---|
| R | 1.7502 | 0.0046 |
| l_a | 301.471 | 0.090 |
| ωb | 0.02236 | 0.00015 |
| ρ(R, l_a) | 0.46 | - |

### Note on EU Model
CMB distance priors probe z ~ 1090, far above the EU transition (z_trans ~ 0.2). Therefore:
- CMB cannot distinguish EU from ΛCDM
- Expected result: w₀ unconstrained, Ωm ~ 0.315

### References
- Planck Collaboration 2018, A&A 641, A6 (Cosmological Parameters)
- Chen, Kumar & Ratra 2017, ApJ 835, 86 (Distance Priors)
- Efstathiou & Bond 1999, MNRAS 304, 75 (Shift Parameter)

---

## 3. BAO Analysis

Using BOSS DR12 Ly-α + Galaxy BAO:
- D_H(z)/r_d and D_M(z)/r_d measurements
- r_d = drag epoch sound horizon

### Data Points
| z | Observable | Value | σ |
|---|------------|-------|---|
| 0.122 | D_H/r_d | - | - |
| 0.51 | D_M/r_d | 13.36 | 0.21 |
| 0.61 | D_M/r_d | 14.94 | 0.21 |
| 2.34 | D_H/r_d | 8.86 | 0.29 |

---

## 4. Cosmic Chronometers H(z)

Direct measurements of H(z) from differential age method:
- 28 data points from Moresco et al. (various)
- Cosmology-independent probe

---

## Summary of Approach

| Dataset | Key Feature | EU Sensitivity |
|---------|-------------|----------------|
| SNe | M marginalized analytically | Moderate (probes z < 2) |
| CMB | Distance priors at z* | Low (z* >> z_trans) |
| BAO | D_M/r_d, D_H/r_d | Moderate |
| H(z) | Direct H(z) | High |

All analyses use `emcee` (Foreman-Mackey et al. 2013) with Pool multiprocessing for efficient sampling.

---

*Document prepared for JCAP referee review, January 2026*
