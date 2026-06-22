# NB05 — Nomenclatura das Runs MCMC

**Data:** 2026-05-28  
**Aplicar quando:** Migração AWS → Colab (versão definitiva)  
**Convenção:** Estilo Planck — descritivo por modelo + dados

---

## Mapeamento AWS → Colab

| AWS | Colab (novo) | Prior EU | Estágio | Dados | Propósito | Paper |
|:---:|:-------------|:--------:|:-------:|:-----:|:----------|:-----:|
| A1 | `eu_free_base_s1` | Livre (Δk=3) | Stage 1 (pert OFF) | Base | Covmat generation | App. A |
| A2 | `eu_free_base` | Livre (Δk=3) | Stage 2 (pert ON) | Base | Blind prediction | §4.1 |
| C1 | `eu_fixed_base_s1` | Fixo (Δk=0) | Stage 1 (pert OFF) | Base | Covmat generation | App. A |
| C2 | `eu_fixed_base` | Fixo (Δk=0) | Stage 2 (pert ON) | Base | Resultado principal | §4.2 |
| D1 | `eu_fixed_shoes` | Fixo (Δk=0) | Stage 2 (pert ON) | +SH0ES | H₀ knockout | §4.3 |
| D2 | `eu_fixed_shoes_des` | Fixo (Δk=0) | Stage 2 (pert ON) | +SH0ES+DES | S₈ knockout | §4.3 |

---

## Convenção de nomes

```
eu_{prior}_{dados}[_{stage}]

prior:  free | fixed
dados:  base | shoes | shoes_des
stage:  _s1 (Stage 1, covmat generation) — omitido no Stage 2
```

### Por que `_s1` e não `.1`

Cobaya nomeia chains como `{output_prefix}.1.txt`, `.2.txt`, etc.
Se o nome do run fosse `eu_free_base.1`, as chains seriam `eu_free_base.1.1.txt`
— ambíguo. Com `_s1`, fica `eu_free_base_s1.1.txt` — limpo.

---

## Datasets base

Todas as runs compartilham:

| Likelihood | Sigla |
|:-----------|:-----:|
| Planck PR4 NPIPE CamSpec TTTEEE | CMB |
| Planck PR4 NPIPE lowl TT+EE | CMB |
| Planck PR4 NPIPE lensing φφ | CMB |
| DESI BAO DR2 | BAO |
| Pantheon+ SNe Ia | SNe |

Runs D adicionam:

| Likelihood | Run |
|:-----------|:---:|
| SH0ES + LKI boost (`shoes_lki.py`) | D1, D2 |
| DES-Y3 3×2pt via Cocoa/CosmoLike | D2 |

---

## Referência no paper

No texto do manuscrito, usar nomes curtos:

| Run | Referência no paper |
|:---:|:-------------------:|
| eu_free_base | **EU-blind** |
| eu_fixed_base | **EU** |
| eu_fixed_shoes | **EU+H₀** |
| eu_fixed_shoes_des | **EU+H₀+S₈** |

Stage 1 runs referenciadas no Appendix A como "covmat generation protocol".

---

## Estrutura de pastas no Colab

```
NB05/Colab/
├── eu_free_base_s1/       (ex-A1)
│   ├── chains/
│   ├── covmat/
│   └── yaml/
├── eu_free_base/          (ex-A2)
├── eu_fixed_base_s1/      (ex-C1)
├── eu_fixed_base/         (ex-C2)
├── eu_fixed_shoes/        (ex-D1)
└── eu_fixed_shoes_des/    (ex-D2)
```

> **Nota:** Os nomes A1/A2/C1/C2/D1/D2 permanecem na pasta AWS intocados.
> Esta nomenclatura aplica-se apenas à versão definitiva (Colab + paper).
