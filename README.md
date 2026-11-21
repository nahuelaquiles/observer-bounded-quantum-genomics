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
2. **Module 2 (`scripts/02_rate_distortion_and_circadian_kinetics.py`):** Single-cell snapshot Poisson erasure ($99.60\%$–$99.96\%$ zero-lesion promoters; binary occupancy entropy $0.0051$–$0.0376\ \mathrm{bits}$), nonlinear Gauss–Newton Poisson MAP decoding vs. exact prior-averaged van Trees Fisher bounds and Gaussian Fisher information proxy $I_{\mathrm{GF}}(\nu)$ across $\nu \in [1, 10^6]$ cells, $17.1\times$ nucleosomal inward-occlusion MSE penalty, and finite-memory rate-distortion crossovers ($C_B \in \{2,4,8,16\}\ \mathrm{bits}$, memory-limited onset at $K = 1, 1, 2, 5$).
3. **Module 3 (`scripts/03_ecoli_mmr_quotient_and_human_promoters.py`):** Empirical uncollapsing of mutational quotients across $27$ independent *E. coli* `REL606` LTEE/MAE populations ($10,467$ SNPs; exact lineage Mann–Whitney $p = 1.88 \times 10^{-4}$; $\cos(\mathrm{MMR}, \mathrm{mutT}) = 0.9818$; flanking SantaLucia stacking Spearman $r = 0.8735$, $p = 9.87 \times 10^{-6}$) and 4-covariate Mahalanobis propensity-score matching of $10$ human circadian clock promoters against $1,500$ `GRCh38` MANE Select CpG-island controls (all outcome BH $q = 0.820$), backed by portable repository shards (`data/human_grch38_mane_promoters_part1.csv.gz` & `part2.csv.gz`).
4. **Module 4 (`scripts/04_alphagenome_circadian_regulome.py`):** Exploratory 1-Mb-context AlphaGenome variant scoring restricted strictly to primary adult human liver (`UBERON:0002107`) across $13$ `GRCh38`-verified $\text{8-oxoG} \to \mathrm{T}$ promoter transversions (`DNASE` $12.74\times$, nominal $p = 0.0451$, 6-modality BH $q = 0.186$; `ATAC` $7.67\times$, $p = 0.0932$, $q = 0.186$; `RNA_SEQ` $1.85\times$, $p = 0.0653$, $q = 0.186$) plus 32-bp *CRY1* single-nucleotide ISM attribution.
5. **Module 5 (`scripts/06_upgraded_github_and_astrodynamics_cosmic_chirality.py`):** Upgraded $300\ \mathrm{kbp}$ (`2,000` windows) *E. coli* `REL606` quantum-DNA codebase audit ($67\times$ trace-normalization GC artifact compression) and **exploratory physical modeling** of $8,766$-hour J2000 Earth–Sun solar-neutrino quadrature ($4.15 \times 10^{-15}\ \mathrm{diploid\text{-}genome^{-1}\ yr^{-1}}$ vs. $6.86 \times 10^{-13}\ \mathrm{whole\text{-}cell^{-1}\ yr^{-1}}$) and $240,000$-trajectory 3D helical B-DNA lepton scattering ($\sigma_{\mathrm{ion}}/\sigma_{\nu e} = 2.71 \times 10^{27}$, $\psi_{\mathrm{crit}} = 26.9^\circ$).
6. **Module 6 (`scripts/08_table6_empirical_wetlab_validation.py` & `scripts/09_generate_fig6.py`):** Stage-specific empirical validation across *E. coli* `REL606` *mutS* thermal MA (`Chu et al. 2018`, comparison across both temperature intervals), cryo-EM nucleosome–hOGG1 structures (`8VWT`, `8VWV`, `9EOZ`, `0` clashes; `1KX5` $54.2$ outward vs. $785.3$ inward mean core clashes, $8.29\times$ descriptive gate), donor-matched 24-hour primary human islet scRNA-seq (`GSE218316`: `GSM6740278`, `GSM6740279`, `GSM6740282`; $n = 2,465$ cells, all $22$ target genes quantified and $19$ used by the combined classifier, strictly disjoint cell pools up to $\nu = 100$ giving Combined `CV AUC = 1.000` / `LODO AUC = 0.960` vs. Clock-only `LODO AUC = 0.347`), host-cell-stratified HepG2 ($9.59\times$, $p = 7.35 \times 10^{-5}$) and K562 ($1.80\times$, $p = 4.58 \times 10^{-5}$) saturation MPRA (`GSE126550`), $14$ *E. coli* GO-repair MA VCF cohorts ($\Delta\textit{mutM}\Delta\textit{mutY}$ $98.62\%$ $\mathrm{G{:}C \to T{:}A}$), and an exploratory two-tier polarized-lepton beamline falsification specification.

---

## Rerunning the analyses

Use Python 3.12 and install the pinned computational dependencies.

```bash
python -m pip install -r requirements.txt
```

Stage public inputs in the layout below. Run these commands from the repository root. The output directories hold the new run and the archived AlphaGenome tables used by Figure 4.

```bash
export OBQG_RAW_DIR="${PWD}/raw"
export OBQG_DATA_DIR="${PWD}/rerun/data"
export OBQG_ASSET_DIR="${PWD}/rerun/assets"
mkdir -p "$OBQG_RAW_DIR" "$OBQG_DATA_DIR" "$OBQG_ASSET_DIR"
cp data/module4_*.csv data/module4_summary.json "$OBQG_DATA_DIR/"

python scripts/01_open_quantum_resolvent_and_nogo.py
python scripts/02_rate_distortion_and_circadian_kinetics.py
python scripts/03_ecoli_mmr_quotient_and_human_promoters.py
python scripts/06_upgraded_github_and_astrodynamics_cosmic_chirality.py
python scripts/08_table6_empirical_wetlab_validation.py
python scripts/05_generate_figures.py
python scripts/07_generate_fig5.py
python scripts/09_generate_fig6.py
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

The islet script quantifies 22 genes and evaluates the clock panel and a combined classifier with 19 features. The Row 3 JSON lists the combined classifier inputs. At each ensemble size it reports the number of eligible folds and test ensembles.

Module 3 uses the two bundled 4,570-promoter GRCh38 shards. Pass `--section promoters` or `--section ecoli` to run either part individually. To use the legacy Feather input, install `pyarrow` and stage `human_grch38_mane_promoters.feather` in `OBQG_RAW_DIR`.

To request live AlphaGenome scores, install the AlphaGenome SDK, configure `ALPHAGENOME_API_KEY` in the environment, and run `python scripts/04_alphagenome_circadian_regulome.py` before generating Figure 4. Variant scores and reference/alternate tracks use a 1-Mb model context. The 32-position CRY1 ISM uses a 16-kb context.

---

## Repository Structure

- `observer_bounded_dna_evolution.tex` — Main *Science Advances* Research Article (`6` figures, `4` tables, `49` full-author references).
- `observer_bounded_dna_supplementary_materials.tex` / `.pdf` — Supplementary Text S1–S3 and Tables S1–S4.
- `cover_letter_science_advances.tex` / `.pdf` — *Science Advances* editorial cover letter.
- `response_to_editorial_review.tex` / `.pdf` & `RESPONSE_TO_EDITORIAL_REVIEW.md` — Point-by-point response to biophysics editorial and reproducibility audits.
- `scripts/` — End-to-end computational, statistical, figure-generation, and audit scripts (`01`–`09` + `analysis_paths.py`).
- `data/` — Derived CSV/JSON outputs plus portable `GRCh38` MANE promoter sequence shards (`human_grch38_mane_promoters_part1/2.csv.gz`).
- `assets/` — Publication-ready 300-DPI PNG and vector SVG figures (`fig1`–`fig6`).

## License

MIT License (see `LICENSE`).
