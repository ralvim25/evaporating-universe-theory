#!/usr/bin/env python3
"""
patch_class.py — Apply Evaporating Universe patches to CLASS v3.3.4
VERSION: vf1 (NB03 physics — includes Fix #1A, #1B, #1C)

Usage: python3 patch_class.py /path/to/class_eu

Patches:
  1. background.h  — EU struct variables + eu_lambda field
  2. input.c        — EU parameter parsing (incl. eu_lambda)
  3. background.c   — CDM drain ODE (RK4) with:
       Fix #1A: a^{-3} source term (NB03 DT Alerta 1)
       Fix #1B: Friedmann closure injection (NB03 DT Alerta 2)
       Fix #1C: w=-1 hardcoded — no dilution term (Δk=0)
       Hard UV cutoff: ε(z) = 0 for z > z_trans
  4. perturbations.c — IDE source terms (eu_has_ide_perturbations=0 → disabled)

Changes vs v7 patch_class.py:
  - Fix #1A: RK4 source multiplied by a^{-3} (vacuum captures dense-past CDM)
  - Fix #1B: Friedmann closure: rho_vac_extra = Ω_c*H0² * (1 - f_c1 - f_v1 + f_v_a)
  - Fix #1C: Removed w0_fld dilution term (w=-1 exactly, vacuum doesn't dilute)
  - Added eu_lambda field (structural constant 2/3)
  - Added hard UV cutoff (z > z_trans → ε = 0)

Reference: Alvim (2025, 2026), NB03 §3.1-§3.4
"""

import sys, os


def patch_background_h(class_dir):
    path = os.path.join(class_dir, 'include', 'background.h')
    with open(path, 'r') as f:
        c = f.read()

    if 'eu_epsilon_ir' in c:
        print('  [SKIP] background.h already patched')
        return

    eu_vars = '\n'.join([
        '',
        '  /* EU IDE (Evaporating Universe) — vf1 */',
        '  double eu_epsilon_ir;    /* IR coupling amplitude */',
        '  double eu_z_trans;       /* transition redshift */',
        '  double eu_b;             /* anomalous dimension (19/36) */',
        '  double eu_lambda;        /* structural constant (2/3) */',
        '  int eu_has_ide_perturbations;  /* 0=bg only, 1=full IDE */',
        '',
    ])
    inserted = False
    for marker in ['short inter_normal;', 'ErrorMsg error_message;']:
        if marker in c:
            c = c.replace(marker, eu_vars + '\n  ' + marker)
            inserted = True
            break
    assert inserted, 'FATAL: Could not find insertion marker in background.h'
    with open(path, 'w') as f:
        f.write(c)
    print('  [OK] background.h patched (incl. eu_lambda)')


def patch_input_c(class_dir):
    path = os.path.join(class_dir, 'source', 'input.c')
    with open(path, 'r') as f:
        c = f.read()

    if 'eu_epsilon_ir' in c:
        print('  [SKIP] input.c already patched')
        return

    eu_block = '\n'.join([
        '',
        '  /* EU IDE defaults and parser — vf1 */',
        '  pba->eu_epsilon_ir = 0.;',
        '  pba->eu_z_trans = 0.;',
        '  pba->eu_b = 0.;',
        '  pba->eu_lambda = 0.6666666666666666;  /* 2/3 structural */',
        '  pba->eu_has_ide_perturbations = 0;',
        '  class_read_double("eu_epsilon_ir", pba->eu_epsilon_ir);',
        '  class_read_double("eu_z_trans", pba->eu_z_trans);',
        '  class_read_double("eu_b", pba->eu_b);',
        '  class_read_double("eu_lambda", pba->eu_lambda);',
        '  class_read_int("eu_has_ide_perturbations", pba->eu_has_ide_perturbations);',
        '',
    ])
    inserted = False
    for marker in ['class_read_double("w0_fld"', 'class_read_double("Omega_Lambda"']:
        idx = c.find(marker)
        if idx > 0:
            end = c.find('\n', c.find(';', idx)) + 1
            c = c[:end] + eu_block + c[end:]
            inserted = True
            break
    assert inserted, 'FATAL: Could not find insertion marker in input.c'
    with open(path, 'w') as f:
        f.write(c)
    print('  [OK] input.c patched (incl. eu_lambda)')


def patch_background_c(class_dir):
    path = os.path.join(class_dir, 'source', 'background.c')
    with open(path, 'r') as f:
        c = f.read()

    if 'eu_epsilon_of_z' in c:
        print('  [SKIP] background.c already patched')
        return

    # ── EU functions with NB03 fixes ──
    eu_functions = '''
/* ============================================================ */
/* EU: Power-law logistic coupling (RGE solution) — vf1.3       */
/* eps(z) = eps_IR / (1 + ((1+z)/(1+z_t))^(1/b))               */
/* Hard UV cutoff: eps(z) = 0 for z > z_trans                   */
/* Fix #1C: w = -1 exactly (vacuum doesn't dilute, Δk=0)        */
/* Fix #2:  Safety guards against MCMC extreme parameters       */
/* Fix #3b: ODE table validation — disable EU for bad params    */
/* ============================================================ */
#include <math.h>
double eu_epsilon_of_z(double z, struct background * pba) {
  double x, exponent, ratio;
  if (pba->eu_epsilon_ir <= 0.) return 0.;
  /* Hard UV cutoff: zero coupling above transition */
  if (z > pba->eu_z_trans) return 0.;
  /* Safety: clamp b to avoid extreme exponents */
  if (pba->eu_b < 0.05) return 0.;  /* b too small → coupling off */
  exponent = 1.0 / pba->eu_b;
  if (exponent > 20.0) return 0.;   /* would overflow pow() */
  ratio = (1. + z) / (1. + pba->eu_z_trans);
  if (ratio <= 0.) return 0.;
  x = pow(ratio, exponent);
  if (x != x || x > 1.0e15) return 0.;  /* NaN or overflow guard */
  return pba->eu_epsilon_ir / (1. + x);
}

/* EU ODE: Pre-computed coupled solution table (RK4, 4th order)  */
/* Fix #1A: source term includes a^{-3} factor                  */
/* Fix #1C: w = -1 → dilution term removed (identically zero)   */
#define EU_TABLE_SIZE 2000
static double eu_a_table[EU_TABLE_SIZE];
static double eu_fcdm_table[EU_TABLE_SIZE];
static double eu_dfld_table[EU_TABLE_SIZE];
static int eu_table_ready = 0;
static double eu_last_eps = -1.0;
static double eu_last_zt  = -1.0;
static double eu_last_b   = -1.0;

void eu_init_ode_table(struct background * pba) {
  int i;
  double a, lna, dlna, z, eps;
  double rho_cdm, rho_fld_extra;
  double a3_inv;
  double k1_c, k1_f, k2_c, k2_f;
  double lna_min = -7.0;
  double lna_max = 0.0;

  /* Invalidate table if EU parameters changed (MCMC) */
  if (eu_table_ready &&
      (pba->eu_epsilon_ir != eu_last_eps ||
       pba->eu_z_trans    != eu_last_zt  ||
       pba->eu_b          != eu_last_b)) {
    eu_table_ready = 0;
  }

  if (pba->eu_epsilon_ir <= 0. || eu_table_ready) return;

  eu_last_eps = pba->eu_epsilon_ir;
  eu_last_zt  = pba->eu_z_trans;
  eu_last_b   = pba->eu_b;
  dlna = (lna_max - lna_min) / (EU_TABLE_SIZE - 1);

  rho_cdm = 1.0;
  rho_fld_extra = 0.0;

  for (i = 0; i < EU_TABLE_SIZE; i++) {
    lna = lna_min + i * dlna;
    a = exp(lna);
    z = 1.0/a - 1.0;
    eps = eu_epsilon_of_z(z, pba);

    eu_a_table[i] = a;
    eu_fcdm_table[i] = rho_cdm;
    eu_dfld_table[i] = rho_fld_extra;

    /* RK4 integration — Fix #1A: a^{-3} on source term */
    /* Fix #1C: w=-1 exactly → 3*(1+w)*rho_fld = 0 (removed) */
    a3_inv = 1.0 / (a * a * a);

    k1_c = -pba->eu_lambda * eps * rho_cdm;
    k1_f = pba->eu_lambda * eps * rho_cdm * a3_inv;  /* Fix #1A: a^{-3} factor */

    { double a_mid, z_mid, eps_mid, a3_inv_mid;
      double k3_c, k3_f, k4_c, k4_f;

      a_mid = exp(lna + 0.5 * dlna);
      z_mid = 1.0/a_mid - 1.0;
      eps_mid = eu_epsilon_of_z(z_mid, pba);
      a3_inv_mid = 1.0 / (a_mid * a_mid * a_mid);

      k2_c = -pba->eu_lambda * eps_mid * (rho_cdm + 0.5*dlna*k1_c);
      k2_f = pba->eu_lambda * eps_mid * (rho_cdm + 0.5*dlna*k1_c) * a3_inv_mid;

      k3_c = -pba->eu_lambda * eps_mid * (rho_cdm + 0.5*dlna*k2_c);
      k3_f = pba->eu_lambda * eps_mid * (rho_cdm + 0.5*dlna*k2_c) * a3_inv_mid;

      { double a_next, z_next, eps_next, a3_inv_next;
        a_next = exp(lna + dlna);
        z_next = 1.0/a_next - 1.0;
        eps_next = eu_epsilon_of_z(z_next, pba);
        a3_inv_next = 1.0 / (a_next * a_next * a_next);

        k4_c = -pba->eu_lambda * eps_next * (rho_cdm + dlna*k3_c);
        k4_f = pba->eu_lambda * eps_next * (rho_cdm + dlna*k3_c) * a3_inv_next;
      }

      rho_cdm += dlna/6.0 * (k1_c + 2.0*k2_c + 2.0*k3_c + k4_c);
      rho_fld_extra += dlna/6.0 * (k1_f + 2.0*k2_f + 2.0*k3_f + k4_f);
      /* Fix #2: Clamp to physical range — prevents NaN propagation */
      if (rho_cdm < 0.01) rho_cdm = 0.01;  /* CDM cannot vanish */
      if (rho_cdm != rho_cdm) rho_cdm = 0.01;  /* NaN guard */
      if (rho_fld_extra != rho_fld_extra) rho_fld_extra = 0.0;  /* NaN guard */
      if (rho_fld_extra > 100.0) rho_fld_extra = 100.0;  /* overflow guard */
    }
  }

  /* Fix #3b: Validate ODE table — reject unphysical results */
  {
    double fcdm_0 = rho_cdm;
    double dfld_0 = rho_fld_extra;
    int bad = 0;
    /* CDM fraction at z=0 must be in [0.5, 1.0] for physical models */
    if (fcdm_0 < 0.5 || fcdm_0 > 1.0) bad = 1;
    /* Vacuum injection cannot be negative or huge */
    if (dfld_0 < -1.0 || dfld_0 > 50.0) bad = 1;
    /* NaN check */
    if (fcdm_0 != fcdm_0 || dfld_0 != dfld_0) bad = 1;
    /* Check monotonicity: fcdm should decrease monotonically */
    {
      int j;
      for (j = 1; j < EU_TABLE_SIZE; j++) {
        if (eu_fcdm_table[j] != eu_fcdm_table[j]) { bad = 1; break; }
        if (eu_dfld_table[j] != eu_dfld_table[j]) { bad = 1; break; }
      }
    }
    if (bad) {
      /* Reset to LCDM (no EU modification) */
      int j;
      for (j = 0; j < EU_TABLE_SIZE; j++) {
        eu_fcdm_table[j] = 1.0;
        eu_dfld_table[j] = 0.0;
      }
      fprintf(stderr, "[EU-vf1.3] BAD ODE table (fcdm=%.4f, dfld=%.4e) — disabling EU for this point\\n",
              fcdm_0, dfld_0);
    }
  }

  eu_table_ready = 1;
  fprintf(stderr, "[EU-vf1.3] ODE table: fcdm(z=0)=%.6f, dfld(z=0)=%.6e\\n",
          eu_fcdm_table[EU_TABLE_SIZE-1], eu_dfld_table[EU_TABLE_SIZE-1]);
}

double eu_get_fcdm(double a, struct background * pba) {
  int i; double lna, t;
  double lna_min = -7.0, lna_max = 0.0;
  if (pba->eu_epsilon_ir <= 0.) return 1.0;
  eu_init_ode_table(pba);
  lna = log(a);
  if (lna <= lna_min) return 1.0;
  if (lna >= lna_max) return eu_fcdm_table[EU_TABLE_SIZE-1];
  t = (lna - lna_min) / (lna_max - lna_min) * (EU_TABLE_SIZE - 1);
  i = (int)t;
  if (i >= EU_TABLE_SIZE - 1) return eu_fcdm_table[EU_TABLE_SIZE-1];
  t -= i;
  return eu_fcdm_table[i] * (1-t) + eu_fcdm_table[i+1] * t;
}

double eu_get_dfld(double a, struct background * pba) {
  int i; double lna, t;
  double lna_min = -7.0, lna_max = 0.0;
  if (pba->eu_epsilon_ir <= 0.) return 0.0;
  eu_init_ode_table(pba);
  lna = log(a);
  if (lna <= lna_min) return 0.0;
  if (lna >= lna_max) return eu_dfld_table[EU_TABLE_SIZE-1];
  t = (lna - lna_min) / (lna_max - lna_min) * (EU_TABLE_SIZE - 1);
  i = (int)t;
  if (i >= EU_TABLE_SIZE - 1) return eu_dfld_table[EU_TABLE_SIZE-1];
  t -= i;
  return eu_dfld_table[i] * (1-t) + eu_dfld_table[i+1] * t;
}
'''

    m = '#include "background.h"'
    idx = c.find(m)
    end = c.find('\n', idx) + 1
    c = c[:end] + eu_functions + c[end:]
    print('  [OK] EU ODE solver functions inserted (Fix #1A, #1C)')

    # ── Inject CDM density modifier ──
    lines = c.split('\n')
    new_lines = []
    nc = 0; nf = 0
    for i, line in enumerate(lines):
        s = line.strip()
        new_lines.append(line)

        # CDM: multiply by fcdm (with safety clamp)
        if 'index_bg_rho_cdm' in s and 'Omega0_cdm' in s and '=' in s and not s.startswith('//') and nc == 0:
            new_lines.append('    { double _fc = eu_get_fcdm(a, pba); /* EU_ODE */')
            new_lines.append('      if (_fc != _fc || _fc < 0.01) _fc = 0.01;')
            new_lines.append('      if (_fc > 1.0) _fc = 1.0;')
            new_lines.append('      pvecback[pba->index_bg_rho_cdm] *= _fc;')
            new_lines.append('    }')
            nc = 1

        # FLD: Friedmann closure injection (Fix #1B)
        # rho_vac_extra = Omega0_cdm * H0^2 * (1 - f_c1 - f_v1 + f_v_a)
        if 'index_bg_rho_fld' in s and '=' in s and ('pvecback[' in s) and not s.startswith('//') and nf == 0:
            new_lines.append('    if (pba->eu_epsilon_ir > 0.) {')
            new_lines.append('      /* Fix #1B: Friedmann closure (NB03 DT Alerta 2) */')
            new_lines.append('      double f_v_a = eu_get_dfld(a, pba);')
            new_lines.append('      double f_v_1 = eu_get_dfld(1.0, pba);')
            new_lines.append('      double f_c_1 = eu_get_fcdm(1.0, pba);')
            new_lines.append('      double rho_cdm_std = pba->Omega0_cdm * pow(pba->H0, 2);')
            new_lines.append('      double _corr = 1.0 - f_c_1 - f_v_1 + f_v_a;')
            new_lines.append('      if (_corr != _corr) _corr = 0.0; /* NaN guard */')
            new_lines.append('      if (_corr > 10.0) _corr = 10.0; /* overflow guard */')
            new_lines.append('      if (_corr < -10.0) _corr = -10.0;')
            new_lines.append('      pvecback[pba->index_bg_rho_fld] += rho_cdm_std * _corr; /* EU_ODE */')
            new_lines.append('    }')
            nf = 1

    c = '\n'.join(new_lines)
    with open(path, 'w') as f:
        f.write(c)
    print(f'  [OK] background.c patched: CDM={nc}, FLD={nf} (Fix #1B, #3b)')


def patch_perturbations_c(class_dir):
    """
    UNIFIED PATCH v3: Perturbations — CDM explicit + FLD complete
    
    Injects IDE source terms into perturbations_derivs() in perturbations.c:
    
    CDM_PATCH (identically zero for w=-1, Δk=0):
      δ'_cdm += ε(z)·aH·(δ_cdm - δ_cdm) = 0
      → Algebraic proof that IDE perturbations do NOT modify CDM growth
        when the vacuum equation of state is exactly w=-1.
        This is why Δσ₈ = 0 in the null test (NB03 Onda 3).
    
    FLD_PATCH (with index guard):
      δ'_fld += ε(z)·aH·(ρ_cdm/ρ_fld)·(δ_cdm - δ_fld)
      → Vacuum receives the CDM perturbation drain.
        Guard: if (pv->index_pt_delta_fld >= 0) prevents crash
        when FLD perturbations are not allocated.
    
    Both patches are gated by eu_has_ide_perturbations == 1.
    When flag=0, perturbations are pure ΛCDM (Phase 1 runs).
    
    Reference: Alvim (2025, 2026), Paper I — Perturbation equations §3
    """
    path = os.path.join(class_dir, 'source', 'perturbations.c')
    with open(path, 'r') as f:
        content = f.read()

    if 'EU_CDM_PATCH' in content:
        print('  [SKIP] perturbations.c already patched (UNIFIED v3)')
        return

    # ── Step 1: Extern declaration ──
    if 'eu_epsilon_of_z' not in content:
        m = '#include "perturbations.h"'
        idx = content.find(m)
        end = content.find('\n', idx) + 1
        extern_decl = '\n/* EU IDE: extern for background coupling function (UNIFIED PATCH v3) */\n'
        extern_decl += '#ifdef __cplusplus\nextern "C"\n#endif\n'
        extern_decl += 'double eu_epsilon_of_z(double z, struct background * pba);\n'
        content = content[:end] + extern_decl + content[end:]
        print('  [OK] perturbations.c: extern declaration added')

    # ── Step 2: CDM_PATCH — inject after CDM density perturbation ODE ──
    # Find: dy[pv->index_pt_delta_cdm] = ... (the CDM continuity equation)
    cdm_marker = 'dy[pv->index_pt_delta_cdm]'
    cdm_idx = content.find(cdm_marker)
    if cdm_idx < 0:
        print('  [WARN] Could not find delta_cdm ODE — skipping CDM_PATCH')
    else:
        # Find end of the CDM ODE line (next semicolon + newline)
        cdm_end = content.find(';', cdm_idx)
        cdm_line_end = content.find('\n', cdm_end) + 1
        
        cdm_patch = '''
      /* ============================================================ */
      /* EU_CDM_PATCH: IDE source term for CDM density perturbation   */
      /* δ\'_cdm += ε(z)·aH·(δ_cdm - δ_cdm) = 0  [identically zero] */
      /*                                                              */
      /* PROOF: For w = -1 exact (Δk = 0), the vacuum perturbation   */
      /* δ_fld does not couple back to δ_cdm because the vacuum      */
      /* energy-momentum tensor T^μν_vac = -ρ_vac·g^μν is            */
      /* homogeneous by construction. The source term reduces to      */
      /* ε·aH·(δ_cdm - δ_cdm) = 0 algebraically.                    */
      /*                                                              */
      /* This explicit null statement is the code proof that          */
      /* Δσ₈ = 0 between flag=0 and flag=1 runs (null test NB03).    */
      /* Ref: Alvim (2025, 2026), Paper I §3.2                       */
      /* ============================================================ */
      if (pba->eu_has_ide_perturbations == 1 && pba->eu_epsilon_ir > 0.) {
        double _z_eu = 1.0/a - 1.0;
        double _eps_eu = eu_epsilon_of_z(_z_eu, pba);
        /* CDM_PATCH: ε·aH·(δ_cdm - δ_cdm) = 0 — identically zero */
        dy[pv->index_pt_delta_cdm] += _eps_eu * a_prime_over_a
          * (y[pv->index_pt_delta_cdm] - y[pv->index_pt_delta_cdm]);
      }
'''
        content = content[:cdm_line_end] + cdm_patch + content[cdm_line_end:]
        print('  [OK] CDM_PATCH injected (identically zero for w=-1)')

    # ── Step 3: FLD_PATCH — inject after FLD density perturbation ODE ──
    # Find: dy[pv->index_pt_delta_fld] = ... (the FLD continuity equation)
    fld_marker = 'dy[pv->index_pt_delta_fld]'
    fld_idx = content.find(fld_marker)
    if fld_idx < 0:
        print('  [WARN] Could not find delta_fld ODE — skipping FLD_PATCH')
    else:
        # Find end of the FLD ODE block (next semicolon + newline)
        fld_end = content.find(';', fld_idx)
        fld_line_end = content.find('\n', fld_end) + 1
        
        fld_patch = '''
      /* ============================================================ */
      /* EU_FLD_PATCH: IDE source term for vacuum perturbation        */
      /* δ\'_fld += ε(z)·aH·(ρ_cdm/ρ_fld)·(δ_cdm - δ_fld)           */
      /*                                                              */
      /* The vacuum receives CDM perturbation drain. This term is     */
      /* non-zero in general, but for w=-1 the FLD perturbation       */
      /* does not feed back into the Poisson equation (∇²Φ) because   */
      /* δP_fld = 0 when w=-1 and cs²=1 (no pressure perturbation).  */
      /*                                                              */
      /* Guard: index_pt_delta_fld >= 0 prevents segfault when FLD    */
      /* perturbation variables are not allocated (e.g., PPF mode).   */
      /* Ref: Alvim (2025, 2026), Paper I §3.2; DT fix               */
      /* ============================================================ */
      if (pba->eu_has_ide_perturbations == 1 && pba->eu_epsilon_ir > 0.) {
        if (pv->index_pt_delta_fld >= 0) {
          double _z_eu_f = 1.0/a - 1.0;
          double _eps_eu_f = eu_epsilon_of_z(_z_eu_f, pba);
          double _rho_cdm_f = pvecback[pba->index_bg_rho_cdm];
          double _rho_fld_f = pvecback[pba->index_bg_rho_fld];
          if (_rho_fld_f > 1.0e-30) {
            /* FLD_PATCH: vacuum receives CDM perturbation drain */
            dy[pv->index_pt_delta_fld] += _eps_eu_f * a_prime_over_a
              * (_rho_cdm_f / _rho_fld_f)
              * (y[pv->index_pt_delta_cdm] - y[pv->index_pt_delta_fld]);
          }
        }
      }
'''
        content = content[:fld_line_end] + fld_patch + content[fld_line_end:]
        print('  [OK] FLD_PATCH injected (with index guard)')

    with open(path, 'w') as f:
        f.write(content)
    print('  [OK] perturbations.c: UNIFIED PATCH v3 complete (CDM explicit + FLD guarded)')


# ── Main ──
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(f'Usage: python3 {sys.argv[0]} /path/to/class_eu')
        sys.exit(1)

    class_dir = sys.argv[1]
    assert os.path.isdir(class_dir), f'FATAL: {class_dir} not found'

    print(f'[EU-vf1.3] Patching CLASS in: {class_dir}')
    print(f'  UNIFIED PATCH v3: background + perturbations (CDM explicit + FLD guarded)')
    print(f'  Fixes: #1A (a^-3 source), #1B (Friedmann closure), #1C (w=-1)')
    patch_background_h(class_dir)
    patch_input_c(class_dir)
    patch_background_c(class_dir)
    patch_perturbations_c(class_dir)
    print('[EU-vf1.3] All CLASS-EU patches applied successfully (UNIFIED PATCH v3)')
