# Observer-Bounded Quantum Genomics

Code, derived datasets, and publication figures accompanying the manuscript:

> **Observer-bounded quantum genomics: Sub-thermal no-go limits, repair quotients, and cosmic chiral funneling**  
> Nahuel Aquiles Garcia (`gecorpinfo@gmail.com`)  
> GECORP (Buenos Aires, Argentina) & Department of Clinical Biochemistry, Aalborg University Hospital (Aalborg, Denmark)

---

## Repository Contents

- `scripts/` — All 9 deterministic Python 3.11 analysis and figure-generation scripts (`01` through `09`):
  1. `01_open_quantum_resolvent_and_nogo.py`: KMS G:C Davies–Lindblad diagonalization, Mori–Zwanzig Schur resolvent projection, moment annihilation, 14-decade Floquet–Lindblad weak-field scan, and 21-bp guanine hole migration (Module 1).
  2. `02_rate_distortion_and_circadian_kinetics.py`: 72-hour multi-harmonic circadian OGG1-BER Poisson ensembles, Gauss–Newton MAP solver, and Clarke–Barron bounds (Module 2).
  3. `03_ecoli_mmr_quotient_and_human_promoters.py`: 27-population *E. coli* `REL606` SBS-96 repair-quotient uncollapsing and 4,570-gene human `GRCh38` MANE Select promoter Mahalanobis matching (Module 3).
  4. `04_alphagenome_circadian_regulome.py`: 1-Mb-context AlphaGenome human liver (`UBERON:0002107`) circadian variant scoring and *CRY1* in silico mutagenesis (Module 4; requires `ALPHAGENOME_API_KEY`).
  5. `05_generate_figures.py`: Generates publication Figures 1–4 (`assets/fig1_*` through `assets/fig4_*`).
  6. `06_upgraded_github_and_astrodynamics_cosmic_chirality.py`: 300-kbp `REL606` quantum-DNA audit, 8,766-hour J2000 Earth–Sun solar-neutrino ephemeris, and 240,000-trajectory 3D helical scattering (Module 5).
  7. `07_generate_fig5.py`: Generates publication Figure 5 (`assets/fig5_*`).
  8. `08_table6_empirical_wetlab_validation.py`: Cross-cohort wet-lab empirical validations across *E. coli* thermal/RF MA, cryo-EM/X-ray nucleosomes (`8VWT`, `8VWV`, `9EOZ`, `1KX5`), `GSE218316` primary human islets, `GDS3342` human monocytes, `GSE126550` HepG2 MPRA, and 14 *E. coli* repair-knockout WGS MA cohorts (Module 6).
  9. `09_generate_fig6.py`: Generates publication Figure 6 (`assets/fig6_*`).
- `data/` — All 35 derived CSV/JSON summary tables plus `pipeline_execution_audit.json` containing SHA-256 verification checksums.
- `assets/` — Publication Figures 1–6 in both 300-dpi `.png` and vector `.svg` formats.
- `PRE_SUBMISSION_AUDIT_REPORT.md` — Independent clean-room verification and pre-submission audit summary.

## Public Experimental Datasets Analyzed

- `barricklab/LTEE-Ecoli` (*E. coli* `REL606` 27-population whole-genome evolution)
- `garushyants/Ecoli_MutationalSignatures` (14 *E. coli* repair-knockout WGS MA VCF cohorts)
- NCBI GEO `GSE218316` (primary human pancreatic islet scRNA-seq ± H2O2), `GDS3342` (human monocyte microarray), and `GSE126550` (HepG2 saturation-mutagenesis MPRA)
- RCSB PDB `8VWT`, `8VWV`, `9EOZ`, `1KX5` (hOGG1–nucleosome cryo-EM and X-ray structures)
- UCSC/Ensembl `GRCh38` GENCODE v46 MANE Select promoters

## Environment & Reproducibility Notes

```bash
pip install -r requirements.txt
```

Scripts `01`–`09` were executed under Python 3.11; top-of-file path constants (`OUT_DATA`, `DATA_DIR`, `ASSET_DIR`, `RAW_DIR`) point to the project workspace layout used during clean-room verification and can be adjusted to local paths before re-running. Script `04_alphagenome_circadian_regulome.py` requires an `ALPHAGENOME_API_KEY` environment variable to query the cloud AlphaGenome API.

## License

Code is released under the MIT License (`LICENSE`).
