# Observer-Bounded Quantum Genomics: Sub-Thermal No-Go Limits, Repair Quotients, and Cosmic Chiral Funneling

**Author:** Nahuel Aquiles Garcia (`gecorpinfo@gmail.com`)  
**Affiliations:**  
1. GECORP, Av. Juan Manuel de Rosas 899, San Miguel del Monte, Buenos Aires, Argentina  
2. Department of Clinical Biochemistry, Aalborg University Hospital, Aalborg, Denmark  

---

## Overview

This repository contains the complete reproducible computational codebase (`scripts/01`–`09`), derived quantitative tables (`data/*.csv`, `data/*.json`), vector and raster publication figures (`assets/fig1`–`fig6`), and typeset LaTeX sources (`observer_bounded_dna_evolution.tex`, `observer_bounded_dna_supplementary_materials.tex`, `references.bib`, `response_to_editorial_review.tex`) for the manuscript:

> **Observer-bounded quantum genomics: Sub-thermal no-go limits, repair quotients, and cosmic chiral funneling**

### Six-Stage Observer-Bounded Markov Channel

$$\Omega \xrightarrow{\mathcal{E}_{\mathrm{chem}}} B \xrightarrow{\Phi_{\mathrm{gate}}} X \xrightarrow{\mathcal{K}_{\mathrm{chrom}}} Z \xrightarrow{\mathcal{W}_{\mathrm{reg}}} R \xrightarrow{\Pi_{\mathcal{O}}} \mathrm{TE}_{\mathcal{O}}$$

Rather than treating microscopic quantum mechanics, chromatin repair, single-cell transcriptomics, and parity-violating astroparticles in disconnected silos—or overclaiming end-to-end causation from single-stage datasets—this work derives exact operator-theoretic and information-theoretic bounds across each transition and evaluates **stage-specific empirical consistency** across orthogonal public cohorts:

1. **Module 1 (`scripts/01_open_quantum_resolvent_and_nogo.py`):** Exact Kubo–Martin–Schwinger (KMS) detailed-balance GKSL generator ($\|\rho_{\mathrm{ss}} - \rho_{\beta}\|_F = 1.61 \times 10^{-11}$), Lumer–Phillips fast-subspace numerical coercivity ($g_2(\mathcal{L}_{PP}) = 1.0092 \times 10^{13}\,\mathrm{s}^{-1}$, $\tau_{\mathrm{mem}} = 98.91\,\mathrm{fs}$, unit contraction prefactor $C=1$), fixed-bath Floquet–Lindblad Drazin ($1.263 \times 10^{-26}$ at $1\,\mu\mathrm{V\,m^{-1}}$) and adiabatic KMS ($3.484 \times 10^{-30}$) sub-thermal weak-field no-go bounds vs. $50\,\mathrm{mK}$ metabolic noise ($2.533 \times 10^{-3}$), and 21-bp $\pi$-stacking tight-binding HOMO hole funneling ($17.95\times$ per-guanine 8-oxoG enrichment at $\mathrm{G}_4$).
2. **Module 2 (`scripts/02_rate_distortion_and_circadian_kinetics.py`):** Exact prior-averaged van Trees–Clarke–Barron Bayesian Cramér–Rao bound ($\overline{F}_1 = \mathbb{E}_\Theta[F_1(\Theta)]$, $1.152\times$ trace gain over origin Fisher $F_1(0)$), converged Gauss–Newton MAP decoding across $\nu \in [1, 10^6]$ cells ($N_{\mathrm{MC}} = 2{,}500$), and finite-memory sensing crossover $K^*(C_B) = 1, 1, 2, 5$ across $C_B \in \{2, 4, 8, 16\}\,\mathrm{bits}$.
3. **Module 3 (`scripts/03_ecoli_mmr_quotient_and_human_promoters.py`):** *Escherichia coli* LTEE (`REL606`, `10,467` SNPs across `27` populations) mismatch-repair (`mutS`/`mutL`) vs. `mutT` 96-trinucleotide simplex uncollapsing ($\cos = 0.9818$ within-MMR vs. $0.0085$ against `mutT`) and Mahalanobis propensity-matched human circadian promoter analysis (`n = 10` core clock promoters vs. `1,500` matched GENCODE v46 MANE Select controls, lowest FDR $q = 0.82$).
4. **Module 4 (`scripts/04_alphagenome_circadian_regulome.py`):** 1-Mb-context AlphaGenome human liver (`UBERON:0002107`) variant scoring and single-nucleotide *in silico* mutagenesis across `CRY1`, `PER1`, `PER2`, and `BMAL1` (`8.54x` DNase and `4.31x` ATAC disruption at SP1/E-box guanine sinks vs. flanking `G1`, $p = 0.0295$).
5. **Module 5 (`scripts/06_upgraded_github_and_astrodynamics_cosmic_chirality.py`):** Repaired 300-kbp `REL606` open-quantum audit (`2,000` windows), 8,766-hour J2000 solar-neutrino ephemeris ($4.15 \times 10^{-15}\,\mathrm{events\,cell^{-1}\,yr^{-1}}$ vs. $2.71 \times 10^{27}\times$ weak-decay lepton electromagnetic gain), and 240,000-trajectory 3D relativistic helical scattering on flat-lying horizontal mineral-adsorbed duplex DNA ($\psi_{\mathrm{crit}} = 26.94^\circ$, $1.281\times$ $\cos^2\theta_z$ zenith gain over isotropic solution).
6. **Module 6 (`scripts/08_table6_empirical_wetlab_validation.py`):** Stage-specific consistency analysis across *E. coli* thermal MA (`Chu et al. 2018`, $n=20$ MA lines/temp, $99.16\%$ KMS Arrhenius agreement), cryo-EM/X-ray hOGG1–nucleosome steric gating (`8VWT`, `8VWV`, `9EOZ`, `1KX5`, $p = 1.18 \times 10^{-32}$), primary human islet scRNA-seq capture sparsity (`GSE218316`, `4,298` cells, LODO $\mathrm{AUC} = 0.502 \to 0.901$), cell-line-stratified living promoter saturation mutagenesis MPRA (`GSE126550`: HepG2 `LDLR`+`F9` non-motif $\mathrm{G}_{3+}$ internal control $1.03\times$ [$p=0.744$] vs. `SP1`/`SRE-1` $\mathrm{G}_{3+}$ $9.59\times$ [$p = 7.35 \times 10^{-5}$]; K562 `PKLR` motif $\mathrm{G}_{3+}$ $1.80\times$ [$p = 4.58 \times 10^{-5}$]), `mutM mutY` GO-repair uncollapsing (`14` WGS MA cohorts, $98.62\%$ $\mathrm{G{:}C \to T{:}A}$), and the two-tier polarized-electron beamline falsification protocol (Tier 1 CISS/DEA $|\mathcal{A}_{\mathrm{eff}}| \ge 10^{-3}$ in $4.81\,\mathrm{s}$/arm vs. Tier 2 pure Dirac SO $|\mathcal{A}_{\mathrm{SO}}| = 2.16 \times 10^{-7}$).

---

## Repository Structure

- `scripts/` — Numbered reproducible Python scripts (`01` through `09`) generating all numerical results and figures.
- `data/` — Derived CSV tables, JSON summaries, and `pipeline_execution_audit.json` with SHA-256 checksums.
- `assets/` — Publication vector (`.svg`) and raster (`.png`) figures (`Fig. 1`–`Fig. 6`).
- `observer_bounded_dna_evolution.tex` — Main *Science Advances* Research Article LaTeX source.
- `observer_bounded_dna_supplementary_materials.tex` — Supplementary Materials (`Supplementary Text S1–S3`, `Tables S1–S4`).
- `references.bib` — Complete BibTeX bibliography (`49` sequential references, zero `et al.` truncations).
- `cover_letter_science_advances.tex` — Editorial Cover Letter source.
- `RESPONSE_TO_EDITORIAL_REVIEW.md` / `response_to_editorial_review.tex` — Point-by-point technical response and scientific rebuttal.
- `SUBMISSION_GUIDE_SCIENCE_ADVANCES.md` & `PRE_SUBMISSION_AUDIT_REPORT.md` — Submission checklist and independent audit report.

## Rerunning the analyses

The checked-in manuscript, tables and figures document the original submission. The audit branch corrects the islet gene mapping and ensemble evaluation. A fresh run produces revised results for that analysis.

Install `requirements.txt` in a Python environment. The audit used Python 3.12, NumPy 2.3.5, SciPy 1.17.0, pandas 2.2.3, scikit-learn 1.8.0 and Biopython 1.88. Run the commands below from the repository root. The scripts default to `data/`, `assets/` and `raw/` in the checkout. Direct output to a separate directory to keep the archived results available for comparison.

```bash
export OBQG_DATA_DIR="${PWD}/rerun/data"
export OBQG_ASSET_DIR="${PWD}/rerun/assets"
export OBQG_RAW_DIR="${PWD}/raw"
mkdir -p "$OBQG_DATA_DIR" "$OBQG_ASSET_DIR" "$OBQG_RAW_DIR"
cp data/module3_human_*.csv data/module4_*.csv data/module4_summary.json "$OBQG_DATA_DIR/"

python3 scripts/01_open_quantum_resolvent_and_nogo.py
python3 scripts/02_rate_distortion_and_circadian_kinetics.py
python3 scripts/03_ecoli_mmr_quotient_and_human_promoters.py --section ecoli
python3 scripts/06_upgraded_github_and_astrodynamics_cosmic_chirality.py
python3 scripts/08_table6_empirical_wetlab_validation.py
python3 scripts/05_generate_figures.py
python3 scripts/07_generate_fig5.py
python3 scripts/09_generate_fig6.py
```

### Public inputs for those commands

Paths below are relative to `OBQG_RAW_DIR`.

| Analyses | Inputs |
| --- | --- |
| Modules 3A and 5 | Copy `reference/REL606.fna` and `reference/REL606.gff3` to `ltee_reference/`, `summary/count.*.masked.csv` to `ltee_summary/`, `LTEE-clone-curated/*50000gen*.gd` to `gd_ltee50k/` and `MAE-clone-curated/JEB*.gd` to `gd_mae/`. [LTEE-Ecoli source commit](https://github.com/barricklab/LTEE-Ecoli/tree/2307710381e66f9ee9b651bb850fdb15b2002a53) |
| Module 6, Row 1 | The script uses the published rounded rates and one-sided test results. [Chu et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6116381/) |
| Module 6, Row 2 | Put the mmCIF files `8VWT.cif`, `8VWV.cif`, `9EOZ.cif` and `1KX5.cif` in `table6/pdb/`. Download each from `https://files.rcsb.org/download/<PDB ID>.cif`. [RCSB PDB](https://www.rcsb.org/) |
| Module 6, Row 3 | Put `GSE218316_metadata_percell.csv.gz` in `table6/downloads/`. Put the sample-specific `*_features.tsv.gz`, `*_barcodes.tsv.gz` and `*_matrix.mtx.gz` from GSM6740278, GSM6740279 and GSM6740282 in `table6/gse218316/`. These are the three 24 h libraries. [GEO series and sample files](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218316) |
| Module 6, Row 4 | Put the promoter FASTA and barcode-to-variant files from GSM3604135, GSM3604143, GSM3604144 and GSM3604148 in `table6/gse126550/`. Use the FASTA from GSM3604144 for both LDLR libraries. Put the 24 count files from GSM3604164–GSM3604169, GSM3604212–GSM3604223 and GSM3604248–GSM3604253 in `table6/downloads/`. [GEO MPRA series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126550) |
| Module 6, Row 5 | Copy `vcf/` and `BSgenome_Ecoli_genomes/` into `table6/ecoli_wgs/`. [Ecoli_MutationalSignatures source commit](https://github.com/garushyants/Ecoli_MutationalSignatures/tree/49545490ee8c335ad0580658320c294be60f41ec) |

Pass `--row 1` through `--row 5` to script `08` to rerun an individual row. Each row writes its own summary. A full `08` run writes the combined summary used by `09`.

The islet analysis joins the published cell metadata to each sample's own Cell Ranger barcode and feature files. The series-level count matrix contains 33,540 rows and its supplementary gene list contains 2,000 entries. For each fold, the revised analysis groups distinct cells within each condition. The output records how many folds and test ensembles support each result. AUCs are reported when every fold supports the requested size.

### Human promoters and AlphaGenome

The commands above use cached human-promoter and AlphaGenome tables for Figure 4. A from-inputs human-promoter rerun requires the original 4,570-promoter `human_grch38_mane_promoters.feather` in `OBQG_RAW_DIR`. The file needs the original MANE Select cohort and extraction coordinates. Run script `03` with `--section promoters` after running `01`.

Script `04` also requires that Feather file and an `ALPHAGENOME_API_KEY`. The revised scorer selects adult-liver tracks using `UBERON:0002107`. The 32-position CRY1 ISM call uses a 16-kb model context. The variant scores and reference/alternate predictions use a 1-Mb context. The Figure 4 reference-base annotations were checked against the UCSC hg38 reference. [UCSC Genome Browser sequence API](https://api.genome.ucsc.edu/getData/sequence?genome=hg38&chrom=chr12&start=107093569&end=107093601)
