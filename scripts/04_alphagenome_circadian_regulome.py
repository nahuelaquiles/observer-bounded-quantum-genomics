#!/usr/bin/env python3
"""Module 4: AlphaGenome 1-Mb-Context Regulatory Distortion & Hole-Sink vs
Flanking Guanine Variant Analysis Across Core Human Circadian Clock Promoters
(CRY1, PER1, PER2, BMAL1).

1. Programmatically identifies GRCh38-verified G-run / regulatory motif sites
   (SP1 GGGCGG / E-box CACGTG) vs flanking isolated guanines in [-500, +500] around
   MANE Select TSSs of CRY1, PER1, PER2, and BMAL1.
2. Scores 8-oxoG -> T oxidative transversions (G>T on coding strand, i.e., G>T on
   plus strand for '+' genes and C>A on plus strand for '-' genes) across 6
   AlphaGenome modalities (DNASE, ATAC, CAGE, RNA_SEQ, CHIP_TF, CHIP_HISTONE)
   in human liver (UBERON:0002107).
3. Saves numerical 32-bp ISM attribution matrices (DNASE & CAGE) and 1-bp local
   REF vs ALT regulatory tracks around the CRY1 -36 bp SP1/G5 hole-sink variant
   (chr12:107093585:C>A), plus clickable AlphaGenome Atlas deep links.
"""

import json
import os
import re
import urllib.parse
import numpy as np
import pandas as pd

from analysis_paths import DATA_DIR as OUT_DATA, load_human_mane_promoters

os.makedirs(OUT_DATA, exist_ok=True)
LIVER_ONTOLOGY = "UBERON:0002107"


def select_verified_promoter_variants():
    prom_df = load_human_mane_promoters()
    targets = ["CRY1", "PER1", "PER2", "BMAL1"]
    selected = []

    for gene in targets:
        row = prom_df[prom_df["gene_name"] == gene].iloc[0]
        chrom = row["chrom"]
        tss = int(row["tss_1based"])
        strand = row["strand"]
        s_plus = row["seq_plus_1001bp"]  # index 0 = tss - 500; index 500 = tss
        start_1b = tss - 500

        # Target base on '+' strand that corresponds to a sense-strand Guanine (8-oxoG -> T transversion):
        # If gene is '+' strand: sense G is '+' G, oxidative G->T is G>T on '+' strand.
        # If gene is '-' strand: sense G is '+' C, oxidative G->T on '-' strand is C>A on '+' strand!
        ref_target = "G" if strand == "+" else "C"
        alt_target = "T" if strand == "+" else "A"

        # Find Category A: Regulatory hole-sink motifs (E-box CACGTG or SP1 GGGCGG/CCGCCC or G>=4 tract within 250 bp of TSS)
        motif_hits = []
        for mtype, pat in [
            ("Ebox_CACGTG", r"CACGTG"),
            ("SP1_GCbox", r"GGGCGG" if strand == "+" else r"CCGCCC"),
            ("G4plus_sink", r"G{4,}" if strand == "+" else r"C{4,}"),
        ]:
            for m in re.finditer(pat, s_plus):
                # Find position inside match carrying ref_target closest to center of motif
                idxs = [m.start() + k for k, ch in enumerate(m.group(0)) if ch == ref_target]
                if not idxs:
                    continue
                mid_idx = idxs[len(idxs) // 2]
                pos_1b = start_1b + mid_idx
                rel_tss = (pos_1b - tss) if strand == "+" else (tss - pos_1b)
                if -350 <= rel_tss <= 150:
                    motif_hits.append((abs(rel_tss), mtype, pos_1b, mid_idx, rel_tss, m.group(0)))

        motif_hits.sort(key=lambda x: (0 if "Ebox" in x[1] else (1 if "SP1" in x[1] else 2), x[0]))
        # Pick top 2 distinct positions (>= 10 bp apart) for Category A
        picked_pos = []
        for _, mtype, pos_1b, mid_idx, rel_tss, mseq in motif_hits:
            if all(abs(pos_1b - p) >= 10 for p in picked_pos):
                assert s_plus[mid_idx] == ref_target
                picked_pos.append(pos_1b)
                selected.append({
                    "gene": gene,
                    "chrom": chrom,
                    "pos_1based": pos_1b,
                    "ref": ref_target,
                    "alt": alt_target,
                    "strand": strand,
                    "rel_tss_bp": rel_tss,
                    "category": "Hole_Sink_Motif",
                    "motif_class": mtype,
                    "local_context_11bp": s_plus[max(0, mid_idx - 5) : min(1001, mid_idx + 6)],
                })
                if len(picked_pos) == 2:
                    break

        # Find Category B: Isolated single guanines (flanked by non-ref_target on both sides, outside SP1/E-box)
        iso_hits = []
        for idx in range(250, 750):
            if s_plus[idx] != ref_target:
                continue
            if s_plus[idx - 1] == ref_target or s_plus[idx + 1] == ref_target:
                continue
            win11 = s_plus[idx - 5 : idx + 6]
            if any(mot in win11 for mot in ("CACGTG", "GGGCGG", "CCGCCC")):
                continue
            pos_1b = start_1b + idx
            rel_tss = (pos_1b - tss) if strand == "+" else (tss - pos_1b)
            if -350 <= rel_tss <= 150 and all(abs(pos_1b - p) >= 15 for p in picked_pos):
                iso_hits.append((abs(rel_tss + 80), pos_1b, idx, rel_tss, win11))

        iso_hits.sort(key=lambda x: x[0])
        for _, pos_1b, idx, rel_tss, win11 in iso_hits[:2]:
            if all(abs(pos_1b - p) >= 15 for p in picked_pos):
                assert s_plus[idx] == ref_target
                picked_pos.append(pos_1b)
                selected.append({
                    "gene": gene,
                    "chrom": chrom,
                    "pos_1based": pos_1b,
                    "ref": ref_target,
                    "alt": alt_target,
                    "strand": strand,
                    "rel_tss_bp": rel_tss,
                    "category": "Isolated_Flanking_G",
                    "motif_class": "Isolated_G1",
                    "local_context_11bp": win11,
                })

    return pd.DataFrame(selected)


def main():
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if not api_key:
        raise RuntimeError("ALPHAGENOME_API_KEY not found in environment")

    from alphagenome.data import genome
    from alphagenome.interpretation import ism
    from alphagenome.models.v1 import dna_client, variant_scorers

    dna_model = dna_client.create(
        api_key=api_key,
        address="dns:///gdmscience.googleapis.com:443",
    )

    cand_df = select_verified_promoter_variants()
    print(f"Selected {len(cand_df)} GRCh38-verified circadian promoter 8-oxoG->T variants:")
    print(cand_df[["gene", "chrom", "pos_1based", "ref", "alt", "rel_tss_bp", "category", "motif_class", "local_context_11bp"]].to_string(index=False))

    scorer_names = ["DNASE", "ATAC", "CAGE", "RNA_SEQ", "CHIP_TF", "CHIP_HISTONE"]
    scorers = [variant_scorers.RECOMMENDED_VARIANT_SCORERS[k] for k in scorer_names]

    all_rows = []
    for _, vrow in cand_df.iterrows():
        var = genome.Variant(
            chromosome=str(vrow["chrom"]),
            position=int(vrow["pos_1based"]),
            reference_bases=str(vrow["ref"]),
            alternate_bases=str(vrow["alt"]),
        )
        interval_1mb = var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
        v_id = f"{vrow['chrom']}:{vrow['pos_1based']}:{vrow['ref']}>{vrow['alt']}"
        atlas_url = f"https://deepmind.google.com/science/alphagenome/atlas?q={urllib.parse.quote(v_id)}&m=variant"

        score_list = dna_model.score_variant(
            interval=interval_1mb,
            variant=var,
            variant_scorers=scorers,
        )
        tidy = variant_scorers.tidy_scores(score_list, match_gene_strand=True)

        rec = dict(vrow)
        rec["variant_id"] = v_id
        rec["atlas_url"] = atlas_url

        # Summarize liver (UBERON:0002107) and pan-track peak quantile/raw scores per modality
        for mod in scorer_names:
            sub_mod = tidy[tidy["variant_scorer"].astype(str).str.upper().str.contains(mod, na=False)]
            if mod == "RNA_SEQ" and "gene_name" in sub_mod.columns:
                sub_gene = sub_mod[sub_mod["gene_name"] == vrow["gene"]]
                if not sub_gene.empty:
                    sub_mod = sub_gene
            if sub_mod.empty:
                rec[f"{mod}_liver_raw"] = np.nan
                rec[f"{mod}_liver_quantile"] = np.nan
                rec[f"{mod}_max_abs_quantile"] = np.nan
                continue

            sub_liv = sub_mod[sub_mod["ontology_curie"].eq(LIVER_ONTOLOGY)].dropna(
                subset=["quantile_score", "raw_score"]
            )
            if not sub_liv.empty:
                idx_top = sub_liv["quantile_score"].abs().idxmax()
                rec[f"{mod}_liver_raw"] = float(sub_liv.loc[idx_top, "raw_score"])
                rec[f"{mod}_liver_quantile"] = float(sub_liv.loc[idx_top, "quantile_score"])
            else:
                rec[f"{mod}_liver_raw"] = np.nan
                rec[f"{mod}_liver_quantile"] = np.nan

            quantiles = sub_mod["quantile_score"].dropna()
            rec[f"{mod}_max_abs_quantile"] = (
                float(quantiles.iloc[quantiles.abs().argmax()]) if not quantiles.empty else np.nan
            )

        all_rows.append(rec)
        print(f"  Scored {v_id} ({vrow['gene']} {vrow['category']}): DNASE_liv_q={rec.get('DNASE_liver_quantile')}, CAGE_liv_q={rec.get('CAGE_liver_quantile')}, RNA_liv_q={rec.get('RNA_SEQ_liver_quantile')}")

    res_df = pd.DataFrame(all_rows)
    res_df.to_csv(f"{OUT_DATA}/module4_alphagenome_circadian_variants.csv", index=False)

    # Also save numerical 32-bp ISM matrix & local +/-250 bp REF vs ALT DNASE & CAGE tracks at CRY1 chr12:107093585:C>A
    cry1_var = genome.Variant("chr12", 107093585, "C", "A")
    ism_interval = cry1_var.reference_interval.resize(32)
    seq_interval_16kb = ism_interval.resize(dna_client.SEQUENCE_LENGTH_16KB)

    for mod_name in ["DNASE", "CAGE"]:
        sc = variant_scorers.RECOMMENDED_VARIANT_SCORERS[mod_name]
        ism_res = dna_model.score_ism_variants(
            interval=seq_interval_16kb,
            ism_interval=ism_interval,
            variant_scorers=[sc],
        )
        scores_flat, vars_flat = [], []
        for adata, *_ in ism_res:
            v_obj = adata.uns["variant"]
            col_mask = adata.var["ontology_curie"].eq(LIVER_ONTOLOGY)
            if not col_mask.any():
                raise ValueError(f"AlphaGenome returned no adult-liver {mod_name} tracks for CRY1 ISM")
            val = float(np.nanmean(adata.X[:, col_mask]))
            scores_flat.append(val)
            vars_flat.append(v_obj)
        mat = ism.ism_matrix(scores_flat, variants=vars_flat)
        ism_df = pd.DataFrame(mat, columns=["A", "C", "G", "T"])
        ism_df["pos_1based"] = np.arange(ism_interval.start + 1, ism_interval.end + 1)
        ism_df.to_csv(f"{OUT_DATA}/module4_cry1_ism_{mod_name.lower()}_liver.csv", index=False)

    # Predict raw 1-bp REF vs ALT tracks across +/- 300 bp around CRY1 TSS (chr12:107093549)
    pred_int = cry1_var.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    v_out = dna_model.predict_variant(
        interval=pred_int,
        variant=cry1_var,
        requested_outputs=[dna_client.OutputType.DNASE, dna_client.OutputType.CAGE],
        ontology_terms=["UBERON:0002107"],
    )
    win = genome.Interval("chr12", 107093249, 107093849)
    ref_dnase = v_out.reference.dnase.slice_by_interval(win)
    alt_dnase = v_out.alternate.dnase.slice_by_interval(win)
    ref_cage = v_out.reference.cage.slice_by_interval(win)
    alt_cage = v_out.alternate.cage.slice_by_interval(win)

    track_df = pd.DataFrame({
        "pos_0based": np.arange(win.start, win.end, ref_dnase.resolution),
        "ref_dnase_liver": np.mean(ref_dnase.values, axis=1),
        "alt_dnase_liver": np.mean(alt_dnase.values, axis=1),
    })
    track_df.to_csv(f"{OUT_DATA}/module4_cry1_ref_alt_dnase_tracks.csv", index=False)

    cage_df = pd.DataFrame({
        "pos_0based": np.arange(win.start, win.end, ref_cage.resolution),
        "ref_cage_liver": np.mean(ref_cage.values, axis=1),
        "alt_cage_liver": np.mean(alt_cage.values, axis=1),
    })
    cage_df.to_csv(f"{OUT_DATA}/module4_cry1_ref_alt_cage_tracks.csv", index=False)

    # Statistical comparison: Hole_Sink_Motif vs Isolated_Flanking_G across |liver quantile| & |raw_score|
    from scipy import stats
    sink = res_df[res_df["category"] == "Hole_Sink_Motif"]
    iso = res_df[res_df["category"] == "Isolated_Flanking_G"]
    mod_list = ["DNASE", "ATAC", "CAGE", "RNA_SEQ", "CHIP_TF", "CHIP_HISTONE"]
    comp_summary = {}
    for mod in mod_list:
        x1 = sink[f"{mod}_liver_raw"].dropna().abs().values
        x2 = iso[f"{mod}_liver_raw"].dropna().abs().values
        q1 = sink[f"{mod}_liver_quantile"].dropna().abs().values
        q2 = iso[f"{mod}_liver_quantile"].dropna().abs().values
        if len(x1) > 0 and len(x2) > 0:
            _, p_raw = stats.mannwhitneyu(x1, x2, alternative="two-sided")
            comp_summary[mod] = {
                "liver_ontology_curie": LIVER_ONTOLOGY,
                "n_sink_variants": int(len(x1)),
                "n_isolated_variants": int(len(x2)),
                "sink_mean_abs_raw": float(np.mean(x1)),
                "isolated_mean_abs_raw": float(np.mean(x2)),
                "fold_raw_ratio": float(np.mean(x1) / max(np.mean(x2), 1e-9)),
                "sink_mean_abs_quantile": float(np.mean(q1)),
                "isolated_mean_abs_quantile": float(np.mean(q2)),
                "mannwhitney_p_raw": float(p_raw),
            }

    valid_mods = [m for m in mod_list if m in comp_summary]
    p_arr = np.array([comp_summary[m]["mannwhitney_p_raw"] for m in valid_mods], dtype=float)
    m_tot = len(p_arr)
    ord_idx = np.argsort(p_arr)
    q_sorted = np.minimum.accumulate((p_arr[ord_idx] * m_tot / np.arange(1, m_tot + 1))[::-1])[::-1]
    q_arr = np.empty_like(p_arr)
    q_arr[ord_idx] = np.clip(q_sorted, 0.0, 1.0)
    for m, qv in zip(valid_mods, q_arr):
        comp_summary[m]["bh_fdr_q_6modalities"] = float(qv)

    with open(f"{OUT_DATA}/module4_summary.json", "w", encoding="utf-8") as f:
        json.dump(comp_summary, f, indent=2)

    print("\n=== MODULE 4 ALPHAGENOME SUMMARY (HOLE-SINK MOTIFS vs ISOLATED FLANKING GUANINES, UBERON:0002107) ===")
    print(json.dumps(comp_summary, indent=2))


if __name__ == "__main__":
    main()
