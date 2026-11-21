#!/usr/bin/env python3
"""Module 3: Empirical Observational Quotient Uncollapsing in Escherichia coli
(WT vs MMR-deficient vs mutT-deficient Genomes) and Propensity-Matched Human
GRCh38 Circadian Promoter Analysis (N_ctrl = 1,500).

Resolves exact repair-gene mutations in E. coli REL606 directly against
REL606.gff3 CDS coordinates (mutS, mutL, mutH, uvrD, mutT, mutY, mutM, dnaQ)
so MMR-deficient lineages (mutS/mutL/uvrD) are cleanly separated from
8-oxo-dGTPase (mutT) lineages and repair-proficient WT lineages.
"""

import argparse
import glob
import json
import os
import re
import numpy as np
import pandas as pd
from scipy import stats

from analysis_paths import DATA_DIR as OUT_DATA, RAW_DIR

os.makedirs(OUT_DATA, exist_ok=True)

COMP = str.maketrans("ACGT", "TGCA")

def revcomp(s: str) -> str:
    return s[::-1].translate(COMP)

PYR_SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
BASES = ["A", "C", "G", "T"]
SBS96_KEYS = [f"{b5}[{sub}]{b3}" for sub in PYR_SUBS for b5 in BASES for b3 in BASES]
SBS96_IDX = {k: i for i, k in enumerate(SBS96_KEYS)}

NN_DG37 = {
    "AA": -1.00, "TT": -1.00,
    "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45,
    "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28,
    "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24,
    "GG": -1.84, "CC": -1.84,
}


def load_rel606_genome_and_repair_loci():
    fna_path = f"{RAW_DIR}/ltee_reference/REL606.fna"
    with open(fna_path, encoding="utf-8") as f:
        seq = "".join(ln.strip() for ln in f if not ln.startswith(">")).upper()
    L = len(seq)
    pyr_tri_counts = {f"{b5}{mid}{b3}": 0 for mid in ["C", "T"] for b5 in BASES for b3 in BASES}
    for i in range(1, L - 1):
        tri = seq[i - 1 : i + 2]
        if any(ch not in "ACGT" for ch in tri):
            continue
        if tri[1] in ("C", "T"):
            pyr_tri_counts[tri] += 1
        else:
            pyr_tri_counts[revcomp(tri)] += 1
    opp_96 = np.array([pyr_tri_counts[f"{k[0]}{k[2]}{k[6]}"] for k in SBS96_KEYS], dtype=float)

    # Parse REL606.gff3 for exact coordinates of DNA repair genes (+/- 150 bp promoter margin)
    target_genes = {"mutS", "mutL", "mutH", "uvrD", "mutT", "mutY", "mutM", "fpg", "dnaQ"}
    repair_intervals = {}
    with open(f"{RAW_DIR}/ltee_reference/REL606.gff3", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("#") or "\tgene\t" not in ln:
                continue
            parts = ln.rstrip("\n").split("\t")
            start, end = int(parts[3]), int(parts[4])
            attrs = parts[8]
            m = re.search(r"Name=([^;]+)", attrs)
            if m and m.group(1) in target_genes:
                gname = "mutM" if m.group(1) == "fpg" else m.group(1)
                repair_intervals[gname] = (start - 150, end + 150)
    print("REL606 DNA repair gene intervals:", repair_intervals)
    return seq, opp_96, repair_intervals


def parse_ecoli_gd_cohorts(rel606_seq, opp_96, repair_intervals):
    ltee_meta = pd.read_csv(f"{RAW_DIR}/ltee_summary/count.LTEE.masked.csv")
    mae_meta = pd.read_csv(f"{RAW_DIR}/ltee_summary/count.MAE.masked.csv")
    print("Distinct LTEE MUTATOR_STATUS values:", ltee_meta["MUTATOR_STATUS"].unique().tolist())

    meta_map = {}
    for _, r in mae_meta.iterrows():
        meta_map[r["file"]] = (f"MAE_pop{r['population']}", int(r["time"]), str(r["MUTATOR_STATUS"]))
    for _, r in ltee_meta.iterrows():
        meta_map[r["file"]] = (str(r["population"]), int(r["time"]), str(r["MUTATOR_STATUS"]))

    pop_snps = {}
    pop_repair_hits = {}
    pop_status = {}

    gd_files = sorted(glob.glob(f"{RAW_DIR}/gd_mae/*.gd") + glob.glob(f"{RAW_DIR}/gd_ltee50k/*.gd"))
    for gpath in gd_files:
        fname = os.path.basename(gpath)
        if fname not in meta_map:
            continue
        pop, gen, mstat = meta_map[fname]
        if not pop.startswith("MAE_") and gen != 50000:
            continue
        pop_snps.setdefault(pop, set())
        pop_repair_hits.setdefault(pop, set())
        pop_status[pop] = mstat

        with open(gpath, encoding="utf-8") as f:
            for ln in f:
                if ln.startswith("#") or not ln.strip():
                    continue
                parts = ln.rstrip("\n").split("\t")
                rec_type = parts[0]
                if rec_type in ("SNP", "INS", "DEL", "SUB", "MOB", "AMP") and len(parts) >= 5 and parts[4].isdigit():
                    pos = int(parts[4])
                    span_end = pos
                    if rec_type == "DEL" and len(parts) >= 6 and parts[5].isdigit():
                        span_end = pos + int(parts[5])
                    for gname, (gs, ge) in repair_intervals.items():
                        if not (span_end < gs or pos > ge):
                            pop_repair_hits[pop].add(f"{rec_type}:{gname}")
                if rec_type == "SNP" and len(parts) >= 6 and parts[4].isdigit():
                    pos = int(parts[4])
                    new_base = parts[5].strip().upper()
                    ref_base = rel606_seq[pos - 1]
                    if ref_base in "ACGT" and new_base in "ACGT" and ref_base != new_base:
                        pop_snps[pop].add((pos, ref_base, new_base))

    pop_records = []
    pop_sbs96 = {}
    for pop in sorted(pop_snps.keys()):
        snps = pop_snps[pop]
        mstat = pop_status[pop]
        rhits = sorted(pop_repair_hits[pop])
        vec = np.zeros(96, dtype=int)
        for pos, ref, alt in snps:
            if pos <= 1 or pos >= len(rel606_seq):
                continue
            tri = rel606_seq[pos - 2 : pos + 1]
            if ref in ("C", "T"):
                key = f"{tri[0]}[{ref}>{alt}]{tri[2]}"
            else:
                rc_tri = revcomp(tri)
                key = f"{rc_tri[0]}[{ref.translate(COMP)}>{alt.translate(COMP)}]{rc_tri[2]}"
            if key in SBS96_IDX:
                vec[SBS96_IDX[key]] += 1

        total_snps = int(np.sum(vec))
        ts_count = int(np.sum(vec[16 * 2 : 16 * 3]) + np.sum(vec[16 * 4 : 16 * 5]))
        tv_at_cg = int(np.sum(vec[16 * 5 : 16 * 6]))
        tv_cg_at = int(np.sum(vec[0:16]))

        has_mutT_struct = any(h in ("INS:mutT", "DEL:mutT", "MOB:mutT", "SUB:mutT") for h in rhits) or (mstat == "mutT")
        has_mmr_struct = any(
            any(h.startswith(p) and h.endswith(g) for p in ("INS:", "DEL:", "MOB:", "SUB:") for g in ("mutS", "mutL", "mutH", "uvrD"))
            for h in rhits
        )
        has_mmr_any = any(any(g in h for g in ("mutS", "mutL", "mutH", "uvrD")) for h in rhits)

        if mstat == "non-mutator":
            regime = "WT_repair_proficient"
        elif mstat == "IS-mutator":
            regime = "IS_mutator"
        elif has_mutT_struct:
            regime = "mutT_8oxodGTP_mutator"
        elif has_mmr_struct or (mstat in ("point-mutator", "nan") and has_mmr_any):
            regime = "MMR_deficient"
        else:
            regime = f"mutator_{mstat}"

        pop_sbs96[pop] = vec
        pop_records.append({
            "population": pop,
            "cohort": "MAE_13750gen" if pop.startswith("MAE_") else "LTEE_50000gen",
            "mutator_status": mstat,
            "repair_locus_hits": ";".join(rhits) if rhits else "none",
            "regime": regime,
            "unique_snps": total_snps,
            "transitions_CT_TC": ts_count,
            "transition_fraction": ts_count / max(total_snps, 1),
            "AT_to_CG_transversions": tv_at_cg,
            "AT_to_CG_fraction": tv_at_cg / max(total_snps, 1),
            "GC_to_TA_transversions": tv_cg_at,
        })

    pop_df = pd.DataFrame(pop_records)
    print("\nVerified E. coli Independent Populations & Repair Genotypes:")
    print(pop_df.to_string(index=False))

    # Build SBS-96 Dirichlet-multinomial posteriors per regime AND equal-weight per-lineage normalized profiles
    alpha0 = 0.5
    alpha_lin = 0.05
    pop_opp_norm = {}
    for p, vec in pop_sbs96.items():
        h_p = (vec + alpha_lin) / opp_96
        pop_opp_norm[p] = h_p / np.sum(h_p)

    regime_vecs = {}
    for reg in pop_df["regime"].unique():
        sub_pops = pop_df[pop_df["regime"] == reg]["population"].tolist()
        counts = np.sum([pop_sbs96[p] for p in sub_pops], axis=0)
        p_raw = (counts + alpha0) / np.sum(counts + alpha0)
        haz = (counts + alpha0) / opp_96
        haz_prob = haz / np.sum(haz)
        p_lineage_mean = np.mean([pop_opp_norm[p] for p in sub_pops], axis=0)
        p_lineage_mean /= np.sum(p_lineage_mean)
        regime_vecs[reg] = {
            "n_independent_lineages": len(sub_pops),
            "populations": sub_pops,
            "counts": counts,
            "p_raw": p_raw,
            "p_opp_norm": haz_prob,
            "p_lineage_mean": p_lineage_mean,
        }

    def cos_dist(u, v):
        return float(1.0 - np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))

    def js_div_bits(p, q):
        m = 0.5 * (p + q)
        return float(0.5 * (np.sum(p * np.log2(p / m)) + np.sum(q * np.log2(q / m))))

    def kl_div_bits(p, q):
        return float(np.sum(p * np.log2(p / q)))

    def shannon_bits(p):
        return float(-np.sum(p * np.log2(p)))

    wt_key = "WT_repair_proficient"
    mmr_key = "MMR_deficient"
    mutt_key = "mutT_8oxodGTP_mutator"

    wt_pops = regime_vecs[wt_key]["populations"]
    mmr_pops = regime_vecs[mmr_key]["populations"]
    obs_cos_wt_mmr = cos_dist(regime_vecs[wt_key]["p_opp_norm"], regime_vecs[mmr_key]["p_opp_norm"])
    obs_js_wt_mmr = js_div_bits(regime_vecs[wt_key]["p_opp_norm"], regime_vecs[mmr_key]["p_opp_norm"])
    obs_kl_mmr_wt = kl_div_bits(regime_vecs[mmr_key]["p_opp_norm"], regime_vecs[wt_key]["p_opp_norm"])

    obs_lin_cos_wt_mmr = cos_dist(regime_vecs[wt_key]["p_lineage_mean"], regime_vecs[mmr_key]["p_lineage_mean"])
    obs_lin_js_wt_mmr = js_div_bits(regime_vecs[wt_key]["p_lineage_mean"], regime_vecs[mmr_key]["p_lineage_mean"])

    # Exact Mann-Whitney U test on per-lineage transition fractions (independent populations)
    wt_ts_fracs = pop_df[pop_df["regime"] == wt_key]["transition_fraction"].values
    mmr_ts_fracs = pop_df[pop_df["regime"] == mmr_key]["transition_fraction"].values
    _, p_mwu_ts_exact = stats.mannwhitneyu(mmr_ts_fracs, wt_ts_fracs, alternative="two-sided", method="exact")

    # Equal-weight population-level permutation test across the (n_wt + n_mmr) independent lineages
    rng = np.random.default_rng(20261121)
    pool_pops = wt_pops + mmr_pops
    n_wt = len(wt_pops)
    N_PERM = 10000
    perm_cos = np.zeros(N_PERM)
    perm_js = np.zeros(N_PERM)
    for b in range(N_PERM):
        perm = rng.permutation(pool_pops)
        h1 = np.mean([pop_opp_norm[p] for p in perm[:n_wt]], axis=0); h1 /= np.sum(h1)
        h2 = np.mean([pop_opp_norm[p] for p in perm[n_wt:]], axis=0); h2 /= np.sum(h2)
        perm_cos[b] = cos_dist(h1, h2)
        perm_js[b] = js_div_bits(h1, h2)
    p_emp_cos = float((1 + np.sum(perm_cos >= obs_lin_cos_wt_mmr)) / (1 + N_PERM))
    p_emp_js = float((1 + np.sum(perm_js >= obs_lin_js_wt_mmr)) / (1 + N_PERM))

    # Neutral Synonymous vs Intergenic test (normalized by exact base-substitution opportunities!)
    # Note: In bacterial genomes under different GC% in coding vs intergenic regions, normalize per AT and CG opportunity class!
    def compute_gc_matched_syn_vs_igd(df_sub):
        # Expected synonymous count if intergenic per-spectrum rates applied to synonymous opportunities:
        spectra = ["AT.GC", "AT.CG", "AT.TA", "CG.TA", "CG.AT", "CG.GC"]
        syn_obs = sum(int(df_sub[f"OBSERVED.SYNONYMOUS.{s}"].sum()) for s in spectra)
        igd_obs = sum(int(df_sub[f"OBSERVED.INTERGENIC.{s}"].sum()) for s in spectra)
        syn_expected_from_igd = 0.0
        syn_poss_tot = 0
        igd_poss_tot = 0
        for s in spectra:
            o_igd = float(df_sub[f"OBSERVED.INTERGENIC.{s}"].sum())
            p_igd = float(df_sub[f"POSSIBLE.INTERGENIC.{s}"].sum())
            p_syn = float(df_sub[f"POSSIBLE.SYNONYMOUS.{s}"].sum())
            syn_expected_from_igd += (o_igd / max(p_igd, 1.0)) * p_syn
            syn_poss_tot += int(p_syn)
            igd_poss_tot += int(p_igd)
        raw_ratio = (syn_obs / syn_poss_tot) / (igd_obs / igd_poss_tot)
        spectrum_adj_ratio = syn_obs / max(syn_expected_from_igd, 1e-9)
        # Exact Poisson test of syn_obs vs syn_expected_from_igd
        p_pois = float(2.0 * min(stats.poisson.cdf(syn_obs, syn_expected_from_igd), stats.poisson.sf(syn_obs - 1, syn_expected_from_igd)))
        return {
            "syn_obs": syn_obs,
            "igd_obs": igd_obs,
            "syn_rate_per_opp": syn_obs / syn_poss_tot,
            "igd_rate_per_opp": igd_obs / igd_poss_tot,
            "raw_syn_to_igd_ratio": float(raw_ratio),
            "spectrum_matched_syn_to_igd_ratio": float(spectrum_adj_ratio),
            "poisson_two_sided_p": min(p_pois, 1.0),
        }

    mae_neutral = compute_gc_matched_syn_vs_igd(mae_meta)
    ltee_50k = ltee_meta[ltee_meta["time"] == 50000].drop_duplicates(subset=["population"])
    mmr_neutral = compute_gc_matched_syn_vs_igd(ltee_50k[ltee_50k["population"].isin(mmr_pops)])

    # Curtin-Hammett & flanking base-stacking thermodynamic regression on MMR-deficient transitions
    mmr_counts = regime_vecs[mmr_key]["counts"]
    ts_rows = []
    for sub in ["C>T", "T>C"]:
        for b5 in BASES:
            for b3 in BASES:
                k = f"{b5}[{sub}]{b3}"
                idx = SBS96_IDX[k]
                c_obs = int(mmr_counts[idx])
                opp = float(opp_96[idx])
                mid = sub[0]
                dg_stack = NN_DG37[b5 + mid] + NN_DG37[mid + b3]
                gc_flank = int(b5 in ("G", "C")) + int(b3 in ("G", "C"))
                rate_per_site = (c_obs + 0.5) / opp
                ts_rows.append({
                    "channel": k,
                    "sub_type": sub,
                    "flank_5p": b5,
                    "flank_3p": b3,
                    "gc_flanks": gc_flank,
                    "dG_stack_kcal_mol": dg_stack,
                    "mmr_def_count": c_obs,
                    "rel606_opportunities": int(opp),
                    "posterior_rate_per_opp": rate_per_site,
                    "log10_rate": float(np.log10(rate_per_site)),
                })
    ts_df = pd.DataFrame(ts_rows)
    r_all, p_all = stats.spearmanr(-ts_df["dG_stack_kcal_mol"], ts_df["posterior_rate_per_opp"])
    r_ct, p_ct = stats.spearmanr(-ts_df[ts_df["sub_type"] == "C>T"]["dG_stack_kcal_mol"], ts_df[ts_df["sub_type"] == "C>T"]["posterior_rate_per_opp"])
    r_tc, p_tc = stats.spearmanr(-ts_df[ts_df["sub_type"] == "T>C"]["dG_stack_kcal_mol"], ts_df[ts_df["sub_type"] == "T>C"]["posterior_rate_per_opp"])

    metrics = {
        "wt_independent_lineages": len(wt_pops),
        "mmr_independent_lineages": len(mmr_pops),
        "mmr_populations": mmr_pops,
        "mutT_populations": regime_vecs[mutt_key]["populations"] if mutt_key in regime_vecs else [],
        "wt_total_unique_snps": int(np.sum(regime_vecs[wt_key]["counts"])),
        "mmr_total_unique_snps": int(np.sum(regime_vecs[mmr_key]["counts"])),
        "mutT_total_unique_snps": int(np.sum(regime_vecs[mutt_key]["counts"])) if mutt_key in regime_vecs else 0,
        "wt_transition_fraction_pooled": float(np.sum(regime_vecs[wt_key]["counts"][32:48]) + np.sum(regime_vecs[wt_key]["counts"][64:80])) / max(int(np.sum(regime_vecs[wt_key]["counts"])), 1),
        "mmr_transition_fraction_pooled": float(np.sum(regime_vecs[mmr_key]["counts"][32:48]) + np.sum(regime_vecs[mmr_key]["counts"][64:80])) / max(int(np.sum(regime_vecs[mmr_key]["counts"])), 1),
        "wt_lineage_transition_mean": float(np.mean(wt_ts_fracs)),
        "wt_lineage_transition_sd": float(np.std(wt_ts_fracs, ddof=1)),
        "mmr_lineage_transition_mean": float(np.mean(mmr_ts_fracs)),
        "mmr_lineage_transition_sd": float(np.std(mmr_ts_fracs, ddof=1)),
        "p_mwu_exact_lineage_transition_WT_vs_MMR": float(p_mwu_ts_exact),
        "mutT_AT_to_CG_fraction": float(np.sum(regime_vecs[mutt_key]["counts"][80:96])) / max(int(np.sum(regime_vecs[mutt_key]["counts"])), 1) if mutt_key in regime_vecs else None,
        "cos_dist_WT_vs_MMR": obs_cos_wt_mmr,
        "js_div_bits_WT_vs_MMR": obs_js_wt_mmr,
        "kl_div_bits_MMR_vs_WT": obs_kl_mmr_wt,
        "lineage_equal_cos_dist_WT_vs_MMR": obs_lin_cos_wt_mmr,
        "lineage_equal_js_div_bits_WT_vs_MMR": obs_lin_js_wt_mmr,
        "p_emp_lineage_cos_WT_vs_MMR_10k_perm": p_emp_cos,
        "p_emp_lineage_js_WT_vs_MMR_10k_perm": p_emp_js,
        "cos_dist_MMR_vs_mutT": cos_dist(regime_vecs[mmr_key]["p_opp_norm"], regime_vecs[mutt_key]["p_opp_norm"]) if mutt_key in regime_vecs else None,
        "js_div_bits_MMR_vs_mutT": js_div_bits(regime_vecs[mmr_key]["p_opp_norm"], regime_vecs[mutt_key]["p_opp_norm"]) if mutt_key in regime_vecs else None,
        "cos_dist_WT_vs_mutT": cos_dist(regime_vecs[wt_key]["p_opp_norm"], regime_vecs[mutt_key]["p_opp_norm"]) if mutt_key in regime_vecs else None,
        "shannon_entropy_bits_WT": shannon_bits(regime_vecs[wt_key]["p_opp_norm"]),
        "shannon_entropy_bits_MMR": shannon_bits(regime_vecs[mmr_key]["p_opp_norm"]),
        "shannon_entropy_bits_mutT": shannon_bits(regime_vecs[mutt_key]["p_opp_norm"]) if mutt_key in regime_vecs else None,
        "entropy_compression_WT_to_MMR_bits": shannon_bits(regime_vecs[wt_key]["p_opp_norm"]) - shannon_bits(regime_vecs[mmr_key]["p_opp_norm"]),
        "mae_neutral_syn_vs_intergenic": mae_neutral,
        "mmr_neutral_syn_vs_intergenic": mmr_neutral,
        "stacking_spearman_all32_r": float(r_all),
        "stacking_spearman_all32_p": float(p_all),
        "stacking_spearman_CT_r": float(r_ct),
        "stacking_spearman_CT_p": float(p_ct),
        "stacking_spearman_TC_r": float(r_tc),
        "stacking_spearman_TC_p": float(p_tc),
    }

    return {
        "pop_df": pop_df,
        "regime_vecs": regime_vecs,
        "ts_df": ts_df,
        "opp_96": opp_96,
        "metrics": metrics,
    }


def g4hunter_max_score(seq: str, win: int = 25) -> float:
    n = len(seq)
    scores = np.zeros(n, dtype=float)
    i = 0
    while i < n:
        ch = seq[i]
        if ch in ("G", "C"):
            j = i
            while j < n and seq[j] == ch:
                j += 1
            run_len = j - i
            scores[i:j] = min(run_len, 4) * (1.0 if ch == "G" else -1.0)
            i = j
        else:
            i += 1
    csum = np.concatenate([[0.0], np.cumsum(scores)])
    return float(np.max(np.abs((csum[win:] - csum[:-win]) / win)))


def compute_10_4bp_helical_psd(seq: str) -> tuple[float, float]:
    n = len(seq)
    dinuc_sig = np.zeros(n - 1, dtype=float)
    for i in range(n - 1):
        di = seq[i : i + 2]
        if di in ("AA", "TT", "TA"):
            dinuc_sig[i] = 1.0
        elif di in ("GC", "CG", "GG", "CC"):
            dinuc_sig[i] = -1.0
    dinuc_sig -= np.mean(dinuc_sig)
    freqs = np.fft.rfftfreq(len(dinuc_sig), d=1.0)
    psd_di = np.abs(np.fft.rfft(dinuc_sig)) ** 2
    band = (freqs >= 1.0 / 10.8) & (freqs <= 1.0 / 10.0)
    bg = (freqs >= 1.0 / 25.0) & (freqs <= 1.0 / 5.0)
    snr_dinuc = float(np.mean(psd_di[band]) / max(np.mean(psd_di[bg]), 1e-12))

    g_sig = np.array([1.0 if ch in ("G", "C") else 0.0 for ch in seq], dtype=float)
    g_sig -= np.mean(g_sig)
    freqs_g = np.fft.rfftfreq(len(g_sig), d=1.0)
    psd_g = np.abs(np.fft.rfft(g_sig)) ** 2
    band_g = (freqs_g >= 1.0 / 10.8) & (freqs_g <= 1.0 / 10.0)
    bg_g = (freqs_g >= 1.0 / 25.0) & (freqs_g <= 1.0 / 5.0)
    snr_g = float(np.mean(psd_g[band_g]) / max(np.mean(psd_g[bg_g]), 1e-12))
    return snr_dinuc, snr_g


def analyze_human_promoters_matched_1500():
    prom_df = pd.read_feather(f"{RAW_DIR}/human_grch38_mane_promoters.feather")
    mod1_trap = pd.read_csv(f"{OUT_DATA}/module1_guanine_hole_trapping.csv")
    phi_by_len = {int(r["g_run_length"]): float(r["mean_phi_8oxoG_per_tract_G"]) for _, r in mod1_trap.iterrows()}

    feats = []
    for _, row in prom_df.iterrows():
        s = row["seq_sense_1001bp"]
        L = len(s)
        n_g = s.count("G")
        n_c = s.count("C")
        gc_frac = (n_g + n_c) / L
        cpg = sum(1 for i in range(L - 1) if s[i : i + 2] == "CG")
        cpg_oe = cpg / max((n_c * n_g) / L, 1e-6)

        g_runs = [len(m.group(0)) for m in re.finditer(r"G+", s)] + [len(m.group(0)) for m in re.finditer(r"C+", s)]
        max_grun = max(g_runs) if g_runs else 0
        g3plus_tracts = sum(1 for rl in g_runs if rl >= 3)
        g4h = g4hunter_max_score(s)

        trap_sum = sum(rl * phi_by_len[min(rl, 5)] for rl in g_runs)
        mean_phi_per_G = trap_sum / max(n_g + n_c, 1)
        snr_dinuc_10_4, snr_g_10_4 = compute_10_4bp_helical_psd(s)
        ebox_cnt = len(re.findall(r"CACGTG", s))
        sp1_cnt = len(re.findall(r"GGGCGG|CCGCCC", s))

        feats.append({
            "gene_name": row["gene_name"],
            "chrom": row["chrom"],
            "tss_1based": int(row["tss_1based"]),
            "strand": row["strand"],
            "gc_frac": gc_frac,
            "cpg_oe": cpg_oe,
            "cpg_count": cpg,
            "total_guanines_duplex": n_g + n_c,
            "max_g_run": max_grun,
            "g3plus_tracts": g3plus_tracts,
            "g4hunter_max": g4h,
            "hole_trap_burden": trap_sum,
            "mean_phi_8oxoG_per_G": mean_phi_per_G,
            "psd_10_4bp_dinuc_snr": snr_dinuc_10_4,
            "psd_10_4bp_gc_snr": snr_g_10_4,
            "ebox_CACGTG_count": ebox_cnt,
            "sp1_GGGCGG_count": sp1_cnt,
        })

    fdf = pd.DataFrame(feats)
    clock_set = {"BMAL1", "CLOCK", "CRY1", "CRY2", "DBP", "NPAS2", "NR1D1", "NR1D2", "PER1", "PER2"}
    fdf["is_clock"] = fdf["gene_name"].isin(clock_set)

    clock_df = fdf[fdf["is_clock"]].copy().sort_values("gene_name")
    pool_df = fdf[~fdf["is_clock"]].copy()
    cgi_pool = pool_df[(pool_df["gc_frac"] >= 0.55) & (pool_df["cpg_oe"] >= 0.60)].copy()

    cov_cols = ["gc_frac", "cpg_oe", "max_g_run", "g4hunter_max"]
    mu_clk = clock_df[cov_cols].mean().values
    cov_mat = pool_df[cov_cols].cov().values + 1e-6 * np.eye(len(cov_cols))
    inv_cov = np.linalg.inv(cov_mat)

    diffs = cgi_pool[cov_cols].values - mu_clk
    maha_d2 = np.sum((diffs @ inv_cov) * diffs, axis=1)
    cgi_pool = cgi_pool.assign(mahalanobis_d2=maha_d2).sort_values("mahalanobis_d2")
    matched_1500 = cgi_pool.head(1500).copy()

    test_metrics = [
        "gc_frac", "cpg_oe", "max_g_run", "g4hunter_max",
        "g3plus_tracts", "hole_trap_burden", "mean_phi_8oxoG_per_G",
        "psd_10_4bp_dinuc_snr", "psd_10_4bp_gc_snr", "sp1_GGGCGG_count", "ebox_CACGTG_count",
    ]
    comp_rows = []
    rng = np.random.default_rng(20261122)
    N_MC_NULL = 10000
    for col in test_metrics:
        x_clk = clock_df[col].values
        x_ctl = matched_1500[col].values
        pooled_sd = np.sqrt(0.5 * (np.var(x_clk, ddof=1) + np.var(x_ctl, ddof=1)))
        delta_std = float((np.mean(x_clk) - np.mean(x_ctl)) / max(pooled_sd, 1e-12))
        _, p_mwu = stats.mannwhitneyu(x_clk, x_ctl, alternative="two-sided")
        null_means = np.array([np.mean(rng.choice(x_ctl, size=len(x_clk), replace=False)) for _ in range(N_MC_NULL)])
        obs_diff = abs(np.mean(x_clk) - np.mean(x_ctl))
        p_emp = float((1 + np.sum(np.abs(null_means - np.mean(x_ctl)) >= obs_diff)) / (1 + N_MC_NULL))
        comp_rows.append({
            "feature": col,
            "clock_mean": float(np.mean(x_clk)),
            "clock_sd": float(np.std(x_clk, ddof=1)),
            "matched1500_mean": float(np.mean(x_ctl)),
            "matched1500_sd": float(np.std(x_ctl, ddof=1)),
            "all4560_bg_mean": float(pool_df[col].mean()),
            "delta_std_vs_matched1500": delta_std,
            "mannwhitney_p": float(p_mwu),
            "finite_sample_p_emp": p_emp,
        })

    comp_df = pd.DataFrame(comp_rows)
    out_mask = comp_df["feature"].isin([
        "g3plus_tracts", "hole_trap_burden", "mean_phi_8oxoG_per_G",
        "psd_10_4bp_dinuc_snr", "psd_10_4bp_gc_snr", "sp1_GGGCGG_count", "ebox_CACGTG_count"
    ])
    pvals = comp_df.loc[out_mask, "mannwhitney_p"].values
    m_tests = len(pvals)
    order = np.argsort(pvals)
    qvals = np.minimum.accumulate((pvals[order] * m_tests / np.arange(1, m_tests + 1))[::-1])[::-1]
    q_full = np.empty_like(pvals); q_full[order] = np.clip(qvals, 0.0, 1.0)
    comp_df["bh_fdr_q"] = np.nan
    comp_df.loc[out_mask, "bh_fdr_q"] = q_full

    return clock_df, matched_1500, comp_df


def run_ecoli_analysis():
    rel_seq, opp_96, repair_intervals = load_rel606_genome_and_repair_loci()
    ecoli_res = parse_ecoli_gd_cohorts(rel_seq, opp_96, repair_intervals)

    ecoli_res["pop_df"].to_csv(f"{OUT_DATA}/module3_ecoli_populations.csv", index=False)
    ecoli_res["ts_df"].to_csv(f"{OUT_DATA}/module3_ecoli_mmr_stacking_regression.csv", index=False)

    sbs_rows = []
    for i, k in enumerate(SBS96_KEYS):
        row = {"sbs96_channel": k, "sub_type": k[2:5], "rel606_opportunities": int(opp_96[i])}
        for reg, d in ecoli_res["regime_vecs"].items():
            row[f"{reg}_count"] = int(d["counts"][i])
            row[f"{reg}_opp_norm_prob"] = float(d["p_opp_norm"][i])
        sbs_rows.append(row)
    pd.DataFrame(sbs_rows).to_csv(f"{OUT_DATA}/module3_ecoli_sbs96_spectra.csv", index=False)

    with open(f"{OUT_DATA}/module3_summary.json", "w", encoding="utf-8") as f:
        json.dump(ecoli_res["metrics"], f, indent=2)

    print("\n=== MODULE 3A KEY METRICS ===")
    for k, v in ecoli_res["metrics"].items():
        print(f"  {k}: {v}")


def run_promoter_analysis():
    clock_df, matched_1500, comp_df = analyze_human_promoters_matched_1500()
    clock_df.to_csv(f"{OUT_DATA}/module3_human_circadian_10_promoters.csv", index=False)
    matched_1500.to_csv(f"{OUT_DATA}/module3_human_matched_1500_controls.csv", index=False)
    comp_df.to_csv(f"{OUT_DATA}/module3_human_promoter_comparison.csv", index=False)

    print("\n=== MODULE 3B CIRCADIAN PROMOTERS (n=10) ===")
    print(clock_df[["gene_name", "chrom", "tss_1based", "strand", "gc_frac", "cpg_oe", "max_g_run", "g3plus_tracts", "g4hunter_max", "psd_10_4bp_dinuc_snr", "sp1_GGGCGG_count", "ebox_CACGTG_count"]].to_string(index=False))
    print("\n=== MODULE 3B MATCHED COMPARISON (n=10 vs n=1,500) ===")
    print(comp_df.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description="Run the Module 3 analyses.")
    parser.add_argument("--section", choices=("all", "ecoli", "promoters"), default="all")
    args = parser.parse_args()
    if args.section in ("all", "ecoli"):
        run_ecoli_analysis()
    if args.section in ("all", "promoters"):
        run_promoter_analysis()


if __name__ == "__main__":
    main()
