# AAAS *Science Advances* Complete Submission Guide

**Manuscript Title (106 / 135 characters):**  
`Observer-bounded quantum genomics: Sub-thermal no-go limits, repair quotients, and cosmic chiral funneling`

**Short Title (44 / 50 characters):**  
`Observer-bounded quantum and chiral genomics`

**Teaser (119 / 125 characters):**  
`Physical channel bounds reveal which quantum, regulatory, and cosmic-chiral signals living cells can sense or remember.`

**Author:** Nahuel Aquiles Garcia (`gecorpinfo@gmail.com`) — Solo & Corresponding Author  
**Submission Portal:** AAAS Content Tracking System (CTS) — <https://cts.sciencemag.org/>  
*(Log in with your ORCID-linked AAAS account → Select **Science Advances** → **Submit a new manuscript** → Article Type: **Research Article**)*

---

## 1. Exact Files to Upload to AAAS CTS (`cts.sciencemag.org`)

All files below have been generated, verified, and formatted specifically for **AAAS *Science Advances*** inside `observer_bounded_dna_paper/`:

| Upload Slot in AAAS CTS | Exact File to Upload | Size / Specs | Purpose in Review & Production |
| :--- | :--- | :--- | :--- |
| **1. Cover Letter** | `cover_letter_science_advances.pdf` *(editable source: `cover_letter_science_advances.tex`)* | 2 pages, 0.22 MB PDF | Read first by the Deputy & Associate Editors; states broad advance, stage-specific consistency across five terrestrial cohorts plus two-tier polarized-lepton experimental test, no dual submission, data availability, and AI workflow note *(review the italicized bracketed note on page 2 before uploading)*. |
| **2. Combined Evaluation PDF** *(Primary Evaluation File)* | `observer_bounded_dna_combined_submission.pdf` | 31 pages, 5.88 MB (<25 MB limit) | **Single merged PDF** containing the Main Manuscript (`Figs. 1–6` + `Tables 1–4` in place, 26 pages) followed immediately by the Supplementary Materials (`Supplementary Text S1–S3` + `Tables S1–S4`, 5 pages). Editors and peer reviewers read this file. |
| **3. Main Manuscript PDF** *(if CTS requests split PDFs for LaTeX)* | `observer_bounded_dna_evolution.pdf` | 26 pages, 5.20 MB (`6 Figs + 4 Tables = 10` display items) | Standalone *Science Advances* Research Article (`141`-word Abstract, `14,749` words before `References and Notes` `< 15,000` limit, `49` full-author sequential references, solo-author `"I"` voice). |
| **4. Supplementary Materials PDF** | `observer_bounded_dna_supplementary_materials.pdf` | 5 pages, 0.53 MB (<25 MB limit) | Formatted per *Science Advances* Supplementary Materials layout (`Supplementary Text S1 to S3` and `Tables S1 to S4`, synchronized to the main reference numbers `(1)`–`(49)` without a separate bibliography). |
| **5. Manuscript Source File (LaTeX Bundle)** | `science_advances_latex_source.zip` | 4.95 MB flat `.zip` archive | Self-contained flat archive containing `observer_bounded_dna_evolution.tex`, pre-compiled `observer_bounded_dna_evolution.bbl` (required by AAAS LaTeX system), `references.bib`, `observer_bounded_dna_supplementary_materials.tex`, and `Fig1.png`–`Fig6.png` (300 dpi). |
| **6. Individual High-Resolution Figures** *(mostly requested at Revision; optional Auxiliary at Initial)* | `assets/fig1_...svg` / `.png` through `assets/fig6_...svg` / `.png` | 300 dpi PNG + vector SVG | Vector `.svg` and 300-dpi `.png` originals for `Fig. 1` to `Fig. 6`, plus derived CSV/JSON tables in `data/` and reproducible Python scripts `01`–`09` in `scripts/`. |
| **Companion Pre-Submission Editorial Rebuttal** *(for internal / AI-editor second-round review; do not upload on initial fresh submission unless replying to review)* | `response_to_editorial_review.pdf` *(+ `RESPONSE_TO_EDITORIAL_REVIEW.md`, `response_to_editorial_review.tex`)* | 5 pages, 0.50 MB PDF | Point-by-point technical response documenting all 6 code/dataset fixes (Part I) and the 4 quantitative counter-arguments preserving the unified six-stage Markov channel (Part II). |

> **Important Note on Permanent Code DOI Before Final Acceptance:**  
> *Science Advances* requires custom analysis code (`scripts/01_...py` through `scripts/09_...py` and `data/*.csv`) to be archived in a permanent public repository with a DOI (such as **Zenodo**) by the revision/acceptance stage. All updated scripts, datasets, figures, and LaTeX sources are synced to <https://github.com/nahuelaquiles/observer-bounded-quantum-genomics>; link that GitHub repository to <https://zenodo.org/> to mint a `10.5281/zenodo.XXXXXXX` DOI upon release and insert it into the **Data and materials availability** statement at revision/acceptance.

---

## 2. Copy-Paste Metadata for the AAAS CTS Web Form

When filling out the online submission steps on `https://cts.sciencemag.org/`, copy-paste the verified blocks below directly:

### Step A — Article Type & Disciplinary Classification
- **Article Type:** `Research Article`
- **Primary Section / Grouping:** `Physical and Materials Sciences` (or `Biomedicine and Life Sciences` — both co-handle interdisciplinary biophysics)
- **Suggested Subject Areas / Keywords:**
  1. `Biophysics` / `Biological Physics`
  2. `Genetics and Genomics`
  3. `Chemical Physics` / `Quantum Information`
  4. `Evolutionary Biology`
  5. `Astrobiology`

### Step B — Title, Short Title, Teaser & Abstract
- **Full Title (`106` chars; limit `135`):**
  ```text
  Observer-bounded quantum genomics: Sub-thermal no-go limits, repair quotients, and cosmic chiral funneling
  ```
- **Short Title (`44` chars; limit `50`):**
  ```text
  Observer-bounded quantum and chiral genomics
  ```
- **Teaser (`119` chars; limit `125`, no "We", ends with period):**
  ```text
  Physical channel bounds reveal which quantum, regulatory, and cosmic-chiral signals living cells can sense or remember.
  ```
- **Abstract (`141` words; limit `150` words, single paragraph, zero citations):**
  ```text
  Quantum proton transfer, purine radical-cation migration, and parity-violating astroparticles continually perturb DNA, yet cells access these microscopic events only after femtosecond hydration dephasing, chromatin-gated repair, and sparse molecular recording. Without unified physical channel bounds, sub-thermal electromagnetic claims and cosmic-chiral hypotheses remain disconnected from quantitative genomics. Here, I formulate mutagenesis and lesion readout as a six-stage observer-bounded Markov channel and test stage-specific consistency across bacterial evolution, cryo-electron microscopy nucleosomes, human islet single-cell transcriptomes, and HepG2 and K562 promoter assays. Mori–Zwanzig resolvent and Floquet–Lindblad bounds prove that sub-thermal oscillatory fields (1 microvolt per meter) shift mutation hazards at least 23 orders below 50 mK thermal noise, whereas guanine stacking concentrates 8-oxoguanine 17.95-fold and weak-decay polarized leptons amplify chiral ionization 2.71 × 10^27-fold over solar neutrinos inside a 26.9° helical cone. Stage-resolved cohorts support each terrestrial bottleneck while isolating polarized-lepton mutagenesis for laboratory testing.
  ```

### Step C — Author, ORCID, Funding, CRediT, Competing Interests & Data Availability
- **Corresponding & Sole Author:** Nahuel Aquiles Garcia (`gecorpinfo@gmail.com`) — **Make sure your ORCID iD is linked in your CTS profile** (mandatory for corresponding authors at *Science Advances*).
- **Funding Statement:**
  ```text
  This research was supported by GECORP and institutional research infrastructure; the author received no specific external grant funding for this work.
  ```
- **Author Contributions (CRediT Taxonomy):**
  ```text
  Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Data curation, Writing—original draft, Writing—review & editing, and Visualization: N.A.G.
  ```
- **Competing Interests:**
  ```text
  The author declares that he has no competing interests.
  ```
- **Data and Materials Availability:**
  ```text
  All data needed to evaluate the conclusions in the paper are present in the paper and/or the Supplementary Materials. Public experimental datasets analyzed in this work are available from barricklab/LTEE-Ecoli, garushyants/Ecoli_MutationalSignatures, NCBI Gene Expression Omnibus accessions GSE218316, GDS3342, and GSE126550, RCSB Protein Data Bank accessions 8VWT, 8VWV, 9EOZ, and 1KX5, and UCSC/Ensembl GRCh38 GENCODE v46 MANE Select. Reproducible Python scripts (01 through 09), derived CSV/JSON data tables, and typeset sources are archived at GitHub (https://github.com/nahuelaquiles/observer-bounded-quantum-genomics) and Zenodo [insert minted Zenodo DOI upon upload].
  ```

### Step D — Five Suggested Reviewers (Required by *Science Advances* CTS)
*Science Advances* requires **5 conflict-free suggested reviewers** (no co-authorship, collaboration, or shared institution within the past 3 years) spanning the paper's subfields:
1. **Open Quantum Biology / Proton Tunneling in DNA** *(e.g., Marco Sacchi or Louie Slocombe, University of Surrey — open-quantum Caldeira–Leggett / QM-MM DNA tautomerism)*
2. **Cosmic-Ray & Spin-Polarized Homochirality** *(e.g., Noémie Globus, UC Santa Cruz / Flatiron / UNAM, or S. Furkan Ozturk, Harvard / Caltech — chiral leptons and spin-polarized prebiotic chemistry)*
3. **Nucleosomal 8-Oxoguanine Repair & Genome-Wide Damage Mapping** *(e.g., Bret D. Freudenthal, University of Kansas Medical Center, or Shana J. Sturla, ETH Zürich — nucleosomal OGG1/8-oxoG structure and single-nucleotide oxidation sequencing)*
4. **Bacterial Experimental Evolution & Mutational Signatures** *(e.g., Patricia L. Foster, Indiana University, or Deepa Agashe, NCBS Bangalore — whole-genome MA and repair-knockout mutational spectra in E. coli)*
5. **Quantum Stochastic Processes / Information Bounds** *(e.g., Kavan Modi, Monash University, or Felix A. Pollock — non-Markovian quantum processes and operational memory bounds)*

*(Verify current institutional email addresses and absence of personal conflicts before entering into CTS.)*

---

## 3. Compliance Summary Checked Against *Science Advances* Rules
- **Title:** `106` characters (`<= 135` max) — **PASS**
- **Short Title:** `44` characters (`<= 50` max) — **PASS**
- **Teaser:** `119` characters (`<= 125` max, does not start with "We", ends with period) — **PASS**
- **Abstract:** `141` words (`<= 150` max, 1 paragraph, 0 citations) — **PASS**
- **Body Length:** `14,749` words before `References and Notes` (`< 15,000` max) — **PASS**
- **Main-Text Display Items:** `6 Figures (Figs. 1–6) + 4 Tables (Tables 1–4) = 10` (`<= 10` max) — **PASS**
- **References:** `49` references (`<= 80` max), numbered in strict citation order `(1)`–`(49)`, **zero `et al.` truncations** — **PASS**
- **Section Order:** `INTRODUCTION` → `RESULTS` → `DISCUSSION` → `MATERIALS AND METHODS` (`Experimental design` → `Modules 1–6` → `Statistical analysis`) → `References and Notes` → `Acknowledgments` → `Supplementary Materials` — **PASS**
- **Voice & Legacy Terms:** Strictly solo-author `"I / my"` (`0` occurrences of `we/our/ours/us`), `0` occurrences of `34 GHz` or `CMB` — **PASS**
