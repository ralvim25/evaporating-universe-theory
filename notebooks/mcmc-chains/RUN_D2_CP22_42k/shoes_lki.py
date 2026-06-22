"""
SH0ES + LKI Likelihood — Evaporating Universe (Run D)
=====================================================
SH0ES H₀ constraint via EU LKI void boost mechanism.

Physics (Linder 2005 + KBC Void):
  H0_LKI = H0_GKI × (1 + 0.102 × Ω_m^0.55)

  where:
    - H0_GKI = global Hubble (from CLASS-EU background)
    - 0.102  = −δ_true/3 = −(−0.306)/3 (KBC void density contrast)
    - f      = Ω_m^0.55 (Linder 2005 growth rate approximation)
    - H0_LKI = local Hubble (what SH0ES actually measures)

Derived quantities exported to chains:
  - H0_LKI:  Local Hubble constant (km/s/Mpc)
  - I_GKI:   Global Kinematic Integral = ∫ ε(z)/(1+z) dz
  - fcdm_z0: CDM survival fraction at z=0 = exp(−λ·I_GKI)

Constraint:
  Gaussian: −0.5 × ((H0_LKI − 73.17) / 0.86)²
  → At Ω_m=0.30, H0=68.9: H0_LKI ≈ 72.5, logp ≈ −0.28 (~0.7σ)
  → ΛCDM at H0=67.4: H0_LKI = 67.4 (no boost), logp ≈ −22.7 (~6.7σ)

Nomenclature (standardized):
  - I_GKI: Global Kinematic Integral = ∫ ε(z)/(1+z) dz
  - H0_GKI: Hubble constant from CLASS-EU background (= H0 param)
  - H0_LKI: Local Hubble constant = H0_GKI × boost

Reference:
  - Alvim (2025, 2026) — EU first-principles EFT
  - Linder (2005) — Growth rate f = Ω_m^γ, γ ≈ 0.55
  - Breuval, Riess et al. (2024) — H₀ = 73.17 ± 0.86 km/s/Mpc
    (4 geometric anchors: MW Cepheids + LMC + NGC 4258 + SMC)
"""

import numpy as np
from cobaya.likelihood import Likelihood


class SH0ES_LKI(Likelihood):
    # ── EU params: FIXED (NB01 UV values, hardcoded as class attributes) ──
    # NOT in Cobaya params: section — passed to CLASS via extra_args.
    # Hardcoded here to compute I_GKI and fcdm_z0 without routing conflicts.
    eu_epsilon_ir: float = 0.04264    # NB01 §6: 1/√(108π) — CW minimum
    eu_z_trans: float = 5.986         # NB01 §5: UV→IR transition
    eu_b: float = 0.52778             # NB01 §3: 19/36 — ghost-graviton vertex

    # ── Structural constant ──
    # λ = (n-1)/n = 2/3 (from eu_derived.py, D-11 consolidation)
    # ❗ Named lambda_LKI because `lambda` is a Python reserved word
    lambda_LKI: float = 2.0 / 3.0

    # ── KBC Void contrast ──
    # δ_true = −0.306 (KBC void real density contrast)
    # boost_coeff = −δ_true/3 = 0.102
    boost_coeff: float = 0.102

    # ── Linder growth rate exponent ──
    # f = Ω_m^γ, γ ≈ 0.55 (Linder 2005)
    gamma_linder: float = 0.55

    # ── SH0ES measurement ──
    # Breuval, Riess et al. (2024)
    H0_shoes: float = 73.17
    H0_shoes_err: float = 0.86

    # ── Derived params exported to chains ──
    params = {
        "H0_LKI": {"derived": True, "latex": r"H_0^{\rm LKI}"},
        "I_GKI": {"derived": True, "latex": r"I_{\rm GKI}"},
        "fcdm_z0": {"derived": True, "latex": r"f_{\rm cdm}(z{=}0)"},
    }

    def get_requirements(self):
        return {"H0": None, "Omega_m": None}

    def logp(self, _derived=None, **params_values):
        # ── H0_GKI: global Hubble from CLASS-EU background ──
        H0_GKI = self.provider.get_param("H0")
        Omega_m = self.provider.get_param("Omega_m")

        # ── I_GKI: Global Kinematic Integral ──
        # ∫₀^z_max ε(z)/(1+z) dz  with logistic coupling ε(z)
        z = np.linspace(0, 50, 2000)
        eps_z = self.eu_epsilon_ir / (
            1.0 + ((1.0 + z) / (1.0 + self.eu_z_trans))**(1.0 / self.eu_b)
        )
        eps_z[z > self.eu_z_trans] = 0.0  # Hard UV cutoff
        _trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
        I_GKI = _trapz(eps_z / (1.0 + z), z)

        # ── fcdm_z0: CDM survival fraction at z=0 ──
        # (consolidated from eu_derived.py — D-11)
        fcdm_z0 = np.exp(-self.lambda_LKI * I_GKI)

        # ── H0_LKI: Local Hubble via KBC Void + Linder growth rate ──
        # H0_LKI = H0_GKI × (1 + 0.102 × Ω_m^0.55)
        # Master Strategy §3.3, Linder (2005)
        f_growth = Omega_m**self.gamma_linder
        boost = 1.0 + self.boost_coeff * f_growth
        H0_LKI = H0_GKI * boost

        # ── Export derived quantities ──
        if _derived is not None:
            _derived["H0_LKI"] = H0_LKI
            _derived["I_GKI"] = I_GKI
            _derived["fcdm_z0"] = fcdm_z0

        # ── Gaussian likelihood: H0_LKI vs SH0ES ──
        return -0.5 * ((H0_LKI - self.H0_shoes) / self.H0_shoes_err)**2
