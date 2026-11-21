# Pre-Submission Independent Audit & Synchronization Report (Unified Rigorous & Reproducibility Audit Revision)

**Manuscript:** *Observer-bounded quantum genomics: Sub-thermal no-go limits, repair quotients, and cosmic chiral funneling*  
**Author:** Nahuel Aquiles Garcia (`gecorpinfo@gmail.com`)  
**Target Journal:** AAAS *Science Advances* (Research Article)  
**Audit Date:** 2026-11-21 (UTC)  
**Overall Verdict:** **PASS — ALL MATHEMATICAL, DATA-MAPPING, AND STANDALONE REPRODUCIBILITY AUDITS INTEGRATED**

---

## 1. Executive Summary

Following both multi-round biophysics editorial review and the independent code-and-data reproducibility audit (`reproducibility-audit` at `9f36fd3`), an end-to-end synchronization was executed across all computational scripts (`scripts/01`–`09` + `scripts/analysis_paths.py`), portable input shards (`data/human_grch38_mane_promoters_part1.csv.gz` & `part2.csv.gz`), derived datasets (`data/`), vector and 300-dpi raster figures (`assets/fig1`–`fig6`), bibliographic records (`references.bib`), and submission documents (`observer_bounded_dna_evolution.tex`, `observer_bounded_dna_supplementary_materials.tex`, `cover_letter_science_advances.tex`, and `response_to_editorial_review.tex`). Every mathematical scoping refinement, per-library single-cell gene alignment, strict adult-liver tissue filter, descriptive nucleosomal clash reporting, and exploratory cosmic/beamline framing has been verified from scratch.

---

## 2. Clean-Room Computational & Statistical Verification (Modules 1–6)

| Module / Channel Stage | Audited Ground-Truth Quantities Verified on Disk | Verification Status |
| :--- | :--- | :--- |
| **Module 1** (`01_open_quantum_resolvent_and_nogo.py`) | 3-state G:C KMS generator at `310.15 K`: Gibbs stationarity Frobenius error `1.61e-11`; spectral gap `g(L_PP) = 1.0110e+13 s^-1` (`tau_mem = 98.91 fs`); Lumer–Phillips coercivity `g_2(L_PP) = 1.0092e+13 s^-1` (`kappa_2(V) = 1.000093`, unit contraction prefactor `C = 1`); exact trace- and Hermiticity-preserving bare-basis GKSL perturbation on `P` (`Hermiticity err = 0.0`, `||Delta M_0,1,2||_2 = 0.0`, `||Delta M_3||_2 = 4.667e39 s^-3`, rel err `7.77e-16`); `1 uV/m` fixed-bath Floquet Drazin `1.263e-26` (crossover `4.48e5 V/m`), adiabatic KMS fixed bare projector `1.599e-28` (crossover `3.98e6 V/m`), adiabatic KMS instantaneous eigenstate projector `3.484e-30` (crossover `2.70e7 V/m`) vs. asymptotic Kubo scale `7.483e-15` (exact DC `3.63e-15`) and `50 mK` floor `2.533e-3`; 21-bp `8-oxoG` enrichment `1.00x, 10.83x, 15.40x, 17.95x, 16.65x` (`G1`–`G5`). | **VERIFIED** |
| **Module 2** (`02_rate_distortion_and_circadian_kinetics.py`) | Snapshot binary lesion-occupancy entropy `5.09e-3, 1.34e-2, 3.76e-2 bits` (`P_0 = 0.9996, 0.9988, 0.9960` at `1e-6, 3e-6, 1e-5` per G) distinguished from Gaussian Fisher proxy `I_GF(nu)` (`4.68e-3 bits` at `nu=1`; `20.14` vs. `20.16 bits` at `nu=10^5`); prior-averaged van Trees Fisher matrix `E_Theta[F_1(Theta)]` across `N_MC = 2,500` Gauss–Newton MAP trials: at `nu = 10^5`, rhythmic OGG1 `0.0130 > 0.0099` vs. inward-occluded SHL `0.2229 > 0.0737` (`17.1x` MAP MSE penalty); memory-limited onset `D_mem > D_sens` at `K = 1, 1, 2, 5` (`0, 0, 1, 4` sensing-limited modes) for `C_B = 2, 4, 8, 16 bits`. | **VERIFIED** |
| **Module 3** (`03_ecoli_mmr_quotient_and_human_promoters.py`) | `10,467` total `REL606` SNPs across `27` populations (`WT = 471`, `MMR = 5,500`, `mutT = 4,442`, `Ara+1 = 54`); lineage MWU `p = 1.88e-4`, `cos_dist(MMR, mutT) = 0.9818`; portable `GRCh38` MANE Select promoter shards (`4,570` promoters in `part1/2.csv.gz`), Mahalanobis-matched clock promoters (`n=10` vs. `N_ctrl = 1,500`, all BH FDR `q = 0.82`). | **VERIFIED** |
| **Module 4** (`04_alphagenome_circadian_regulome.py`) | Exploratory 1-Mb AlphaGenome scoring conditioned strictly on primary adult human liver (`UBERON:0002107`; zero `EFO:0001187` fallback) across `8` hole-sink vs. `5` isolated `G1` variants (`DNASE 12.74x`, nominal `p = 0.0451`, 6-modality BH `q = 0.186`; `ATAC 7.67x`, `p = 0.0932`, `q = 0.186`; `RNA_SEQ 1.85x`, `p = 0.0653`, `q = 0.186`); `CRY1` `-36 bp` `G5` sink (`chr12:107093585:C>A`) peak 8-mer L2 attribution `1.692`. | **VERIFIED** |
| **Module 5** (`06_upgraded_github_and_astrodynamics_cosmic_chirality.py`) | Exploratory physical modeling: `300`-kbp `REL606` audit (`2,000` windows, `0.0` product QMI); `8,766`-h J2000 solar-neutrino ephemeris (`4.15e-15 yr^-1` diploid DNA [`1.08e19` genomes for `5 sigma`] vs. `6.86e-13 yr^-1` whole `1 ng` cell [`6.53e16` cells / `65 t` water]; `sigma_ion/sigma_nue = 2.71e27`); `240,000`-trajectory 3D helical scattering on horizontal mineral D-DNA (`98.29%` inward at `psi > 26.94 deg`, `-6.883e-19 m`, `1.281x` over isotropic); `500 eV` Dirac SO `|A_SO| = 2.160e-7`. | **VERIFIED** |
| **Module 6** (`08_table6_empirical_wetlab_validation.py`) | **Row 1 (`Chu 2018`):** `n = 20` MA lines/temp (`60` total), `37/28 C = 1.613x` (`P_adj = 6e-8`) vs. KMS `1.599x` (`99.16%`). **Row 2 (`8VWT/8VWV/9EOZ/1KX5`):** `0` extruded clashes vs. descriptive `54.24` outward / `785.32` inward octamer-core clashes (`8.29x` modeled gate, no step-wise `p`-value). **Row 3 (`GSE218316` 24h Cell Ranger):** `GSM6740278/79/82` (`n = 2,465` cells, `3` donors, all `22` genes [`10` clock + `12` BER/NRF2] aligned per library, strictly disjoint pools up to `nu=100`): `85.71%` mean per-gene clock zero-UMI (`30.34%` all-10-zero), `nu=1` LODO `0.486` (clock) / `0.529` (combined) -> `nu=100` Clock LODO `0.347` vs. Combined 22-gene CV `1.000` / LODO `0.960`. **Row 4 (`GSE126550` Stratified):** HepG2 (`n=310` SNVs: non-motif `G3+` `1.03x` `p=0.744`; `SP1/SRE-1 G3+` `9.59x` `p=7.35e-5`, `7.24x` in `LDLR`) & K562 (`n=231` SNVs: motif `G3+` `1.80x` `p=4.58e-5`). **Row 5a:** `14` *E. coli* WGS MA cohorts (`98.62%` `G:C->T:A` in `mutM mutY`, `83.6x` over MMR, `1.41–1.55x` `G3+`–`G5+`). **Row 5b:** Exploratory two-tier beamline (`|A_eff| >= 1e-3` in `4.81 s`/arm vs. `|A_SO| = 2.16e-7` requiring `1.93e14` lesions/arm). | **VERIFIED** |

---

## 3. Editorial, Visual QA & Standalone Archive Checks

- **AAAS *Science Advances* Compliance:**
  - Title: `106` chars (`<= 135`) — **PASS**
  - Short Title: `44` chars (`<= 50`) — **PASS**
  - Teaser: `119` chars (`<= 125`, no `"We"`, ends with period) — **PASS**
  - Abstract: `141` words (`<= 150`, single paragraph, zero citations) — **PASS**
  - Main Body Length: `10,657` prose words / `13,132` math-inclusive tokens before `References and Notes` (`< 15,000` limit) — **PASS**
  - Display Items: `6` figures + `4` tables (`= 10`, at the `<= 10` limit) — **PASS**
  - References: `49` sequential references (`(1)`–`(49)`), `0` missing/uncited keys, `0` `et al.` truncations — **PASS**
  - Solo-Author Voice & Legacy Terms: `0` occurrences of `we/our/ours/us/ourselves`; `0` occurrences of `34 GHz`, `CMB`, `GDS3342`, or `EFO:0001187` — **PASS**
- **PDF & Visual Inspection:**
  - `pdftotext` scan for `??` or `[?]` returned `[]` across all `5` compiled PDFs.
  - Visual inspection of rendered preview PNGs across `observer_bounded_dna_evolution.pdf`, `observer_bounded_dna_supplementary_materials.pdf`, `cover_letter_science_advances.pdf`, and `response_to_editorial_review.pdf` confirmed clean margins, legible figures/tables, and zero text overlap or clipping.
- **Standalone LaTeX Source Bundle (`science_advances_latex_source.zip`):**
  - Unpacked into a clean temporary directory and compiled from scratch (`pdflatex` x2 with pre-built `.bbl` and flat `Fig1.png`–`Fig6.png`): exit code `0`, `0` undefined/multiply-defined warnings.

---

## 4. Final Verified Submission Deliverables (`data/pipeline_execution_audit.json`)

| Deliverable File | Pages | Size (Bytes / MB) | SHA-256 Prefix |
| :--- | :---: | :---: | :--- |
| `observer_bounded_dna_combined_submission.pdf` | 32 | `3,457,002` (`3.46 MB`) | `cdde2b917cd90a24` |
| `observer_bounded_dna_evolution.pdf` | 27 | `2,760,735` (`2.76 MB`) | `2bd6d098169ad75d` |
| `observer_bounded_dna_supplementary_materials.pdf` | 5 | `531,315` (`0.53 MB`) | `90fb7a7f730cc522` |
| `cover_letter_science_advances.pdf` | 2 | `217,586` (`0.22 MB`) | `eb6c617f39a85872` |
| `response_to_editorial_review.pdf` | 4 | `430,643` (`0.43 MB`) | `1c5e242e7977dc7e` |
| `science_advances_latex_source.zip` | — | `1,812,202` (`1.81 MB`) | `9415e128fe310416` |
| `observer_bounded_dna_evolution.tex` | — | `139,278` (`0.14 MB`) | `2cc4fde881801dee` |
| `observer_bounded_dna_supplementary_materials.tex` | — | `28,412` (`0.03 MB`) | `792ccdda186f9d17` |
| `cover_letter_science_advances.tex` | — | `8,073` (`0.01 MB`) | `d00923131e9c5541` |
| `response_to_editorial_review.tex` | — | `19,065` (`0.02 MB`) | `a6af762f872b30bc` |
| `RESPONSE_TO_EDITORIAL_REVIEW.md` | — | `5,876` (`0.01 MB`) | `89e9e9533425fa70` |
| `references.bib` | — | `21,542` (`0.02 MB`) | `dc5d9755b01f8221` |
