# Observer-Bounded Quantum Genomics: Sub-Thermal No-Go Limits, Repair Quotients, and Cosmic Chiral Funneling

**Author:** Nahuel Aquiles Garcia (GECORP & Department of Clinical Biochemistry, Aalborg University Hospital)  
**Manuscript Status:** *Science Advances* Submission Package (Main Manuscript + Supplementary Materials + Reproducible Code, Portable Promoter Shards & Derived Data)

---

## Overview

This repository contains the complete reproducible computational, statistical, structural, genomic, and exploratory astrodynamic/helical-scattering pipeline supporting **"Observer-bounded quantum genomics: Sub-thermal no-go limits, repair quotients, and cosmic chiral funneling."**

Living cells and sequencing assays access fast quantum base-pair dynamics only through a six-stage Markov channel:

$$\Omega \xrightarrow{\mathcal{E}_{\mathrm{chem}}} B \xrightarrow{\Phi_{\mathrm{gate}}} X \xrightarrow{\mathcal{K}_{\mathrm{chrom}}} Z \xrightarrow{\mathcal{W}_{\mathrm{reg}}} R \xrightarrow{\Pi_{\mathcal{O}}} \mathrm{TE}_{\mathcal{O}}$$

Across six computational and empirical modules (`scripts/01` through `09`), this repository provides:
1. **Module 1 (`scripts/01_open_quantum_resolvent_and_nogo.py`):** Three-state G:C KMS-compliant Davies–Lindblad open-quantum resolvent projection ($g(\mathcal{L}_{PP}) = 1.011 \times 10^{13}\ \mathrm{s}^{-1}$, Lumer–Phillips coercivity $g_2 = 1.009 \times 10^{13}\ \mathrm{s}^{-1}$, $\tau_{\mathrm{mem}} = 98.91\ \mathrm{fs}$), exact trace- and Hermiticity-preserving bare-basis GKSL moment annihilation ($\|\Delta M_{0,1,2}\|_2 = 0$, $\|\Delta M_3\|_2 = 4.67 \times 10^{39}\ \mathrm{s}^{-3}$), sub-thermal weak-field bounds across fixed-bath Floquet–Lindblad Drazin ($1.26 \times 10^{-26}$ at $1\ \mu\mathrm{V\ m^{-1}}$) and fixed-projector ($1.60 \times 10^{-28}$) / instantaneous-eigenstate ($3.48 \times 10^{-30}$) adiabatic KMS regimes, and 21-bp duplex guanine hole-sink $\text{8-oxoG}$ enrichment ($17.95\times$ in $\mathrm{G}_4$, $95.66\%$ capture).
2. **Module 2 (`scripts/02_circadian_ogg1_ber_rate_distortion.py`):** Single-cell snapshot Poisson erasure ($99.60\%$–$99.96\%$ zero-lesion promoters; binary occupancy entropy $0.0051$–$0.0376\ \mathrm{bits}$), nonlinear Gauss–Newton Poisson MAP decoding vs. exact prior-averaged van Trees Fisher bounds and Gaussian Fisher information proxy $I_{\mathrm{GF}}(\nu)$ across $\nu \in [1, 10^6]$ cells, $17.1\times$ nucleosomal inward-occlusion MSE penalty, and finite-memory rate-distortion crossovers ($C_B \in \{2,4,8,16\}\ \mathrm{bits}$, memory-limited onset at $K = 1, 1, 2, 5$).
3. **Module 3 (`scripts/03_ecoli_ltee_mae_and_human_promoters.py`):** Empirical uncollapsing of mutational quotients across $27$ independent *E. coli* `REL606` LTEE/MAE populations ($10,467$ SNPs; exact lineage Mann–Whitney $p = 1.88 \times 10^{-4}$; $\cos(\mathrm{MMR}, \mathrm{mutT}) = 0.9818$; flanking SantaLucia stacking Spearman $r = 0.8735$, $p = 9.87 \times 10^{-6}$) and 4-covariate Mahalanobis propensity-score matching of $10$ human circadian clock promoters against $1,500$ `GRCh38` MANE Select CpG-island controls (all outcome BH $q = 0.820$), backed by portable repository shards (`data/human_grch38_mane_promoters_part1.csv.gz` & `part2.csv.gz`).
4. **Module 4 (`scripts/04_alphagenome_circadian_promoters.py`):** Exploratory 1-Mb-context AlphaGenome variant scoring restricted strictly to primary adult human liver (`UBERON:0002107`) across $13$ `GRCh38`-verified $\text{8-oxoG} \to \mathrm{T}$ promoter transversions (`DNASE` $12.74\times$, nominal $p = 0.0451$, 6-modality BH $q = 0.186$; `ATAC` $7.67\times$, $p = 0.0932$, $q = 0.186$; `RNA_SEQ` $1.85\times$, $p = 0.0653$, $q = 0.186$) plus 32-bp *CRY1* single-nucleotide ISM attribution.
5. **Module 5 (`scripts/05_github_quantum_dna_upgrade_and_cosmic_chirality.py`):** Upgraded $300\ \mathrm{kbp}$ (`2,000` windows) *E. coli* `REL606` quantum-DNA codebase audit ($67\times$ trace-normalization GC artifact compression) and **exploratory physical modeling** of $8,766$-hour J2000 Earth–Sun solar-neutrino quadrature ($4.15 \times 10^{-15}\ \mathrm{diploid\text{-}genome^{-1}\ yr^{-1}}$ vs. $6.86 \times 10^{-13}\ \mathrm{whole\text{-}cell^{-1}\ yr^{-1}}$) and $240,000$-trajectory 3D helical B-DNA lepton scattering ($\sigma_{\mathrm{ion}}/\sigma_{\nu e} = 2.71 \times 10^{27}$, $\psi_{\mathrm{crit}} = 26.9^\circ$).
6. **Module 6 (`scripts/08_table6_empirical_wetlab_validation.py` & `scripts/09_generate_fig6.py`):** Stage-specific empirical validation across *E. coli* `REL606` *mutS* thermal MA (`Chu et al. 2018`, $99.16\%$ Arrhenius agreement), cryo-EM nucleosome–hOGG1 structures (`8VWT`, `8VWV`, `9EOZ`, `0` clashes; `1KX5` $54.2$ outward vs. $785.3$ inward mean core clashes, $8.29\times$ descriptive gate), donor-matched 24-hour primary human islet scRNA-seq (`GSE218316`: `GSM6740278`, `GSM6740279`, `GSM6740282`; $n = 2,465$ cells, all $22$ target genes quantified, strictly disjoint cell pools up to $\nu = 100$ giving Combined `CV AUC = 1.000` / `LODO AUC = 0.960` vs. Clock-only `LODO AUC = 0.347`), host-cell-stratified HepG2 ($9.59\times$, $p = 7.35 \times 10^{-5}$) and K562 ($1.80\times$, $p = 4.58 \times 10^{-5}$) saturation MPRA (`GSE126550`), $14$ *E. coli* GO-repair MA VCF cohorts ($\Delta\textit{mutM}\Delta\textit{mutY}$ $98.62\%$ $\mathrm{G{:}C \to T{:}A}$), and an exploratory two-tier polarized-lepton beamline falsification specification.

---

## Reproducing the Full Pipeline

### 1. Python Environment
```bash
pip install -r requirements.txt
```

### 2. Portable Inputs & Raw Public Archives (`OBQG_RAW_DIR`)
- Modules 1 and 2 are pure numerical simulations requiring no external downloads.
- Module 3 automatically reads the bundled `GRCh38` MANE Select promoter shards (`data/human_grch38_mane_promoters_part1.csv.gz` and `part2.csv.gz`) when `hg38.fa` is not staged locally, and reads `LTEE-Ecoli` from `$OBQG_RAW_DIR/LTEE-Ecoli` (`git clone https://github.com/barricklab/LTEE-Ecoli.git`).
- Module 4 requires `ALPHAGENOME_API_KEY` in the environment to re-query live 1-Mb AlphaGenome tracks (`data/module4_alphagenome_circadian_variants.csv` and `data/module4_cry1_reference_bases.csv` are archived in `data/`).
- Module 6 resolves raw PDB (`8VWT.cif`, `8VWV.cif`, `9EOZ.cif`, `1KX5.cif`), GEO (`GSE218316_RAW.tar` or extracted `GSM6740278/79/82` 24-h Cell Ranger tarballs; `GSE126550_RAW.tar`), and `Ecoli_MutationalSignatures` (`git clone https://github.com/garushyants/Ecoli_MutationalSignatures.git`) via `OBQG_RAW_DIR` (see `scripts/analysis_paths.py`).

```bash
export OBQG_RAW_DIR=/path/to/staged_public_archives
python3 -u scripts/01_open_quantum_resolvent_and_nogo.py
python3 -u scripts/02_circadian_ogg1_ber_rate_distortion.py
python3 -u scripts/03_ecoli_ltee_mae_and_human_promoters.py
python3 -u scripts/04_alphagenome_circadian_promoters.py   # requires ALPHAGENOME_API_KEY
python3 -u scripts/05_github_quantum_dna_upgrade_and_cosmic_chirality.py
python3 -u scripts/06_generate_all_figures.py
python3 -u scripts/08_table6_empirical_wetlab_validation.py
python3 -u scripts/09_generate_fig6.py
python3 -u scripts/07_verify_submission_package.py
```

---

## Repository Structure

- `observer_bounded_dna_evolution.tex` / `.pdf` — Main *Science Advances* Research Article (`6` figures, `4` tables, `49` full-author references).
- `observer_bounded_dna_supplementary_materials.tex` / `.pdf` — Supplementary Text S1–S3 and Tables S1–S4.
- `observer_bounded_dna_combined_submission.pdf` — Combined single-PDF reviewer/archive edition (Main + Supplement).
- `cover_letter_science_advances.tex` / `.pdf` — *Science Advances* editorial cover letter.
- `response_to_editorial_review.tex` / `.pdf` & `RESPONSE_TO_EDITORIAL_REVIEW.md` — Point-by-point response to biophysics editorial and reproducibility audits.
- `science_advances_latex_source.zip` — Self-contained Editorial Manager LaTeX source archive.
- `scripts/` — End-to-end computational, statistical, figure-generation, and audit scripts (`01`–`09` + `analysis_paths.py`).
- `data/` — Derived CSV/JSON outputs plus portable `GRCh38` MANE promoter sequence shards (`human_grch38_mane_promoters_part1/2.csv.gz`).
- `assets/` — Publication-ready 300-DPI PNG and vector SVG figures (`fig1`–`fig6`).

## License

MIT License (see `LICENSE`).
