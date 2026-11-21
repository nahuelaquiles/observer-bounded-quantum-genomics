#!/usr/bin/env python3
"""Generate Figure 6 (PNG 300 dpi + SVG):
Empirical Wet-Lab Resolution of the Table 6 Falsification Matrix (Modules 1-6).
"""

import json
import math
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

DATA_DIR = "/home/sandbox/workspace/observer_bounded_dna_paper/data"
ASSET_DIR = "/home/sandbox/workspace/observer_bounded_dna_paper/assets"
os.makedirs(ASSET_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.6,
    "axes.titlesize": 9.3,
    "axes.titleweight": "bold",
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.6,
    "ytick.labelsize": 7.8,
    "legend.fontsize": 7.0,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.9,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

COL_NAVY = "#1B4F72"
COL_CRIMSON = "#C0392B"
COL_TEAL = "#117A65"
COL_AMBER = "#D35400"
COL_PURPLE = "#6C3483"
COL_SLATE = "#566573"


def make_figure6():
    with open(f"{DATA_DIR}/module6_summary.json") as fh:
        m6 = json.load(fh)
    clash_df = pd.read_csv(f"{DATA_DIR}/module6_nucleosome_rotational_clashes.csv")
    sc_df = pd.read_csv(f"{DATA_DIR}/module6_sc_h2o2_ensemble_scaling.csv")
    mpra_df = pd.read_csv(f"{DATA_DIR}/module6_hepg2_mpra_variants.csv")
    ecoli_df = pd.read_csv(f"{DATA_DIR}/module6_ecoli_mutator_wgs.csv")

    fig = plt.figure(figsize=(15.6, 10.8))
    gs_outer = GridSpec(
        2, 2, figure=fig,
        left=0.052, right=0.958, bottom=0.105, top=0.945,
        hspace=0.42, wspace=0.22,
    )

    # =========================================================================
    # PANEL A (gs_outer[0, 0]): Row 1 (Thermal vs AC) + Row 2 (Nucleosome Gate)
    # =========================================================================
    gs_A = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[0, 0], width_ratios=[0.95, 1.30], wspace=0.33)

    # --- A1: Row 1 E. coli Thermal Arrhenius vs Module 1 KMS Prediction ---
    axA1 = fig.add_subplot(gs_A[0, 0])
    r1 = m6["row1_thermal_vs_ac"]
    t_grid = np.linspace(23.5, 38.5, 120)
    kb = 8.617333262e-5
    mod1_curve = np.exp(-0.42 / (kb * (t_grid + 273.15))) / np.exp(-0.42 / (kb * (37.0 + 273.15)))

    axA1.plot(t_grid, mod1_curve, color=COL_NAVY, lw=2.1, label="Mod. 1 KMS (0.42 eV)")
    emp_t = np.array([25.0, 28.0, 37.0])
    emp_y = np.array(r1["emp_rel_rates_25_28_37C"])
    axA1.plot(
        emp_t, emp_y, "o", color=COL_CRIMSON,
        markersize=6.5, markeredgecolor="white", markeredgewidth=0.7,
        label="REL606 mutS ratios (Chu 2018)", zorder=5,
    )
    axA1.axhline(1.0, color=COL_SLATE, ls=":", lw=0.8, alpha=0.6)
    axA1.text(
        24.0, 0.96,
        f"28→37 °C: {r1['emp_ratio_37C_vs_28C']:.2f}× obs\nvs {r1['mod1_predicted_ratio_37C_vs_28C']:.2f}× pred\n({r1['arrhenius_ratio_agreement_pct']:.1f}% match, P = {r1['chu2018_28_vs_37C_padj']:.0e})",
        fontsize=6.8, color=COL_NAVY, fontweight="bold", va="top",
        bbox=dict(boxstyle="round,pad=0.24", fc="#F4F6F7", ec="#BDC3C7", lw=0.7),
    )
    axA1.text(
        29.2, 0.44,
        "2.4/5 GHz Wi-Fi:\nnull (<10⁻¹⁴ pred)",
        fontsize=6.7, color=COL_TEAL, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.22", fc="#E8F8F5", ec=COL_TEAL, lw=0.6),
    )
    axA1.set_xlim(23.5, 38.5)
    axA1.set_ylim(0.39, 1.16)
    axA1.set_xticks([25, 28, 32, 37])
    axA1.set_xlabel("Growth temperature (°C)")
    axA1.set_ylabel("Relative mutation rate (37 °C = 1)")
    axA1.set_title("A1. Row 1: Thermal KMS Gate", loc="left", pad=6)
    axA1.legend(loc="upper left", bbox_to_anchor=(-0.06, -0.18), ncol=1, frameon=False, fontsize=6.8)

    # --- A2: Row 2 Nucleosomal Minor-Groove Rotational Clash Gate ---
    axA2 = fig.add_subplot(gs_A[0, 1])
    r2 = m6["row2_nucleosome_rotational_gate"]
    axA2.axvspan(0, 60, color="#E8F8F5", alpha=0.55, zorder=0)
    axA2.axvspan(120, 180, color="#FDEDEC", alpha=0.55, zorder=0)
    axA2.scatter(
        clash_df["inward_angle_deg"], clash_df["core_clashes_lt_2p5A"],
        s=14, alpha=0.45, color=COL_SLATE, edgecolors="none", rasterized=True, label=f"1KX5 {r2['n_nucleosome_positions_scanned']} bp steps",
    )
    bins = np.linspace(0, 180, 10)
    clash_df["ang_bin"] = pd.cut(clash_df["inward_angle_deg"], bins=bins, include_lowest=True)
    bstats = clash_df.groupby("ang_bin", observed=False)["core_clashes_lt_2p5A"].agg(["mean", "count", "std"]).reset_index()
    bcenters = 0.5 * (bins[:-1] + bins[1:])
    bsem = (bstats["std"] / np.sqrt(np.maximum(bstats["count"], 1))).fillna(0.0)
    axA2.plot(bcenters, bstats["mean"], color=COL_NAVY, lw=2.2, marker="s", markersize=4.2, label="Octamer core mean")
    axA2.fill_between(bcenters, np.maximum(0, bstats["mean"] - bsem), bstats["mean"] + bsem, color=COL_NAVY, alpha=0.18)

    cryo_ang_dict = r2["cryo_em_inward_angles_deg"]
    cryo_x = [cryo_ang_dict["8VWV"], cryo_ang_dict["8VWT"], cryo_ang_dict["9EOZ"]]
    axA2.scatter(
        cryo_x, [0.0, 0.0, 0.0],
        marker="D", s=44, color=COL_CRIMSON, edgecolor="white", lw=0.7, zorder=6,
        label="Cryo-EM 8VWV/8VWT/9EOZ (0 clashes)",
    )
    axA2.text(4.0, 1015, f"Outward\n({100*r2['outward_core_zero_clash_fraction']:.1f}% 0-clash)", fontsize=6.9, color=COL_TEAL, fontweight="bold")
    axA2.text(114.0, 1015, f"Inward occluded\n(mean {r2['inward_mean_core_clashes_2p5A']:.0f})", fontsize=6.9, color=COL_CRIMSON, fontweight="bold")
    axA2.text(
        62.0, 790,
        f"p = {r2['mwu_p_inward_vs_outward_core_clashes']:.1e}\n{r2['outward_vs_inward_accessibility_ratio']:.1f}× gate",
        fontsize=7.0, color=COL_NAVY, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#BDC3C7", lw=0.7),
    )
    axA2.set_xlim(-3, 183)
    axA2.set_ylim(-45, 1180)
    axA2.set_xticks([0, 60, 120, 180])
    axA2.set_xlabel("Minor-groove inward phase θ_in (deg)")
    axA2.set_ylabel("hOGG1–histone core clashes (<2.5 Å)")
    axA2.set_title("A2. Row 2: Nucleosomal SHL Gate", loc="left", pad=6)
    axA2.legend(loc="upper left", bbox_to_anchor=(-0.02, -0.18), ncol=2, frameon=False, fontsize=6.7)

    # =========================================================================
    # PANEL B (gs_outer[0, 1]): Row 3 Primary Human Islet scRNA-seq Scaling
    # =========================================================================
    axB = fig.add_subplot(gs_outer[0, 1])
    axB_r = axB.twinx()
    axB_r.spines["top"].set_visible(False)

    r3 = m6["row3_single_cell_erasure_and_scaling"]
    nu = sc_df["ensemble_cells_nu"].values
    axB.axhline(0.50, color=COL_SLATE, ls="--", lw=0.9, alpha=0.7)
    l1 = axB.plot(
        nu, sc_df["clock10_roc_auc"], color=COL_NAVY, lw=2.3, marker="o", markersize=5.4,
        label="Circadian clock ROC-AUC (NR1D1/DBP)", zorder=5,
    )
    l2 = axB.plot(
        nu, sc_df["combined_clock_redox_roc_auc"], color=COL_TEAL, lw=2.1, marker="s", markersize=4.8,
        label="Combined clock + BER/redox ROC-AUC", zorder=4,
    )
    l3 = axB_r.plot(
        nu, 100.0 * sc_df["all_10_clock_genes_zero_frac"], color=COL_CRIMSON, lw=1.9, ls="-.", marker="^", markersize=4.8,
        label="Zero-UMI clock erasure (%)", zorder=3,
    )

    axB.set_xscale("log")
    axB.set_xlim(0.8, 1900)
    axB.set_ylim(0.46, 1.04)
    axB_r.set_ylim(-4, 108)
    axB.set_xlabel("Pooled primary human pancreatic islet ensemble size ν (cells, log scale)")
    axB.set_ylabel("Cell-disjoint 5-fold CV H₂O₂ ROC-AUC", color=COL_NAVY)
    axB_r.set_ylabel("All-clock zero-UMI erasure (%)", color=COL_CRIMSON)
    axB_r.tick_params(axis="y", colors=COL_CRIMSON)

    axB.annotate(
        f"ν = 1 single cell:\n{100*r3['mean_all10_clock_gene_zero_erasure']:.2f}% per-gene clock 0-UMI\n(OGG1/MUTYH below filter)\nClock AUC = {r3['nu1_clock10_auc']:.3f}",
        xy=(1.0, r3["nu1_clock10_auc"]), xytext=(1.55, 0.59),
        arrowprops=dict(arrowstyle="->", color=COL_CRIMSON, lw=1.0),
        fontsize=7.0, color=COL_CRIMSON, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.24", fc="#FDEDEC", ec=COL_CRIMSON, lw=0.7),
    )
    axB.annotate(
        f"ν = 1,000 human cells:\nClock AUC = {r3['nu1000_clock10_auc']:.3f} (LODO {r3['nu1000_clock10_lodo_auc']:.3f})\nCombined AUC = {r3['nu1000_combined_auc']:.3f}\nCapacity = {r3['nu1000_combined_bits']:.3f} bits",
        xy=(1000.0, r3["nu1000_clock10_auc"]), xytext=(11.0, 0.855),
        arrowprops=dict(arrowstyle="->", color=COL_NAVY, lw=1.0),
        fontsize=7.0, color=COL_NAVY, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.24", fc="#EBF5FB", ec=COL_NAVY, lw=0.7),
    )
    axB.set_title(f"B. Row 3: Single-Cell Erasure vs. Ensemble Scaling (GSE218316, n = {r3['gse218316_total_cells_analyzed']:,})", loc="left", pad=6)

    lines_B = l1 + l2 + l3
    labels_B = [ln.get_label() for ln in lines_B]
    axB.legend(lines_B, labels_B, loc="upper left", bbox_to_anchor=(0.0, -0.18), ncol=2, frameon=False, fontsize=6.9)

    # =========================================================================
    # PANEL C (gs_outer[1, 0]): Row 4 Living Human HepG2 MPRA (GSE126550)
    # =========================================================================
    gs_C = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[1, 0], width_ratios=[1.05, 1.25], wspace=0.31)

    dedup = mpra_df.copy()
    dedup["promoter"] = dedup["library"].replace({"LDLR.2": "LDLR"})
    prom_ox = (
        dedup[dedup["is_8oxoG_to_T_transversion"]]
        .groupby(["promoter", "pos_1based", "ref", "alt", "duplex_G_run_length", "in_SP1_or_Ebox_motif"], as_index=False)
        .agg({"effect_log2_RNA_DNA": "mean"})
    )
    prom_ox["abs_effect_log2"] = prom_ox["effect_log2_RNA_DNA"].abs()

    # --- C1: Stratified |Delta log2(RNA/DNA)| across G-run length & SP1/SRE sinks ---
    axC1 = fig.add_subplot(gs_C[0, 0])
    r4 = m6["row4_hepg2_mpra_hole_sinks"]

    g_f9 = prom_ox[(prom_ox["promoter"] == "F9") & (prom_ox["duplex_G_run_length"] >= 3)]["abs_effect_log2"]
    g_g1 = prom_ox[(prom_ox["duplex_G_run_length"] == 1) & (~prom_ox["in_SP1_or_Ebox_motif"])]["abs_effect_log2"]
    g_g3 = prom_ox[prom_ox["duplex_G_run_length"] == 3]["abs_effect_log2"]
    g_g4p = prom_ox[prom_ox["duplex_G_run_length"] >= 4]["abs_effect_log2"]
    g_sp1 = prom_ox[(prom_ox["duplex_G_run_length"] >= 3) & (prom_ox["in_SP1_or_Ebox_motif"])]["abs_effect_log2"]

    groups = [g_f9, g_g1, g_g3, g_g4p, g_sp1]
    g_labels = [
        f"F9 ctrl\nG3+\n(n={len(g_f9)})",
        f"Isol.\nG1\n(n={len(g_g1)})",
        f"Duplex\nG3\n(n={len(g_g3)})",
        f"HOMO\nG4+\n(n={len(g_g4p)})",
        f"SP1/SRE\nG3+\n(n={len(g_sp1)})",
    ]
    g_means = [float(g.mean()) for g in groups]
    g_sems = [float(g.std(ddof=1) / math.sqrt(max(len(g), 1))) for g in groups]
    g_cols = [COL_SLATE, "#5DADE2", COL_TEAL, COL_AMBER, COL_CRIMSON]

    x_pos = np.arange(len(groups))
    axC1.bar(x_pos, g_means, yerr=g_sems, capsize=3.2, color=g_cols, edgecolor="#2C3E50", lw=0.7, width=0.60)
    for xi, (mv, sv) in enumerate(zip(g_means, g_sems)):
        axC1.text(xi, mv + sv + 0.04, f"{mv:.2f}", ha="center", va="bottom", fontsize=7.0, fontweight="bold", color="#1C2833")

    axC1.text(
        -0.35, 1.62,
        f"SP1/SRE G3+ vs G1:\n{r4['pooled_SP1_sink_vs_G1_fold_ratio']:.2f}× pooled (p={r4['mwu_p_SP1_G3plus_vs_G1']:.1e})\n{r4['ldlr_SP1_sink_vs_G1_fold_ratio']:.2f}× in LDLR (p={r4['ldlr_mwu_p_SP1_sink_vs_G1']:.1e})",
        fontsize=6.7, color=COL_CRIMSON, fontweight="bold", va="top",
        bbox=dict(boxstyle="round,pad=0.22", fc="#FDEDEC", ec=COL_CRIMSON, lw=0.65),
    )
    axC1.set_xticks(x_pos)
    axC1.set_xticklabels(g_labels, fontsize=6.6)
    axC1.set_ylim(0, 1.72)
    axC1.set_ylabel("HepG2 |Δlog₂(RNA/DNA)| per 8-oxoG→T")
    axC1.set_title("C1. Row 4: HepG2 Hole-Sink Penalty", loc="left", pad=6)

    # --- C2: Base-resolution LDLR promoter 8-oxoG -> T profile ---
    axC2 = fig.add_subplot(gs_C[0, 1])
    ldlr_sub = prom_ox[prom_ox["promoter"] == "LDLR"].sort_values("pos_1based")
    is_sink = (ldlr_sub["duplex_G_run_length"] >= 3) & (ldlr_sub["in_SP1_or_Ebox_motif"])
    is_g3_other = (ldlr_sub["duplex_G_run_length"] >= 3) & (~ldlr_sub["in_SP1_or_Ebox_motif"])
    is_g12 = ldlr_sub["duplex_G_run_length"] <= 2

    axC2.axhline(0.0, color=COL_SLATE, lw=0.8, ls="--")
    axC2.axvspan(181, 191, color="#FCF3CF", alpha=0.65, zorder=0)
    axC2.axvspan(199, 209, color="#FDEDEC", alpha=0.65, zorder=0)

    axC2.scatter(
        ldlr_sub.loc[is_g12, "pos_1based"], ldlr_sub.loc[is_g12, "effect_log2_RNA_DNA"],
        s=15, color="#85929E", alpha=0.7, label="G1–G2 flank", zorder=2,
    )
    axC2.scatter(
        ldlr_sub.loc[is_g3_other, "pos_1based"], ldlr_sub.loc[is_g3_other, "effect_log2_RNA_DNA"],
        s=25, color=COL_TEAL, alpha=0.85, label="Other G3+", zorder=3,
    )
    axC2.scatter(
        ldlr_sub.loc[is_sink, "pos_1based"], ldlr_sub.loc[is_sink, "effect_log2_RNA_DNA"],
        s=36, color=COL_CRIMSON, edgecolor="black", lw=0.4, label="SP1/SRE G4–G5 sink", zorder=5,
    )

    sre_p = r4["ldlr_sre1_peak"]["pos"]
    sre_v = r4["ldlr_sre1_peak"]["log2"]
    sre_loss = 100.0 * (1.0 - 2.0 ** sre_v)
    sp1_p = r4["ldlr_sp1_peak"]["pos"]
    sp1_v = r4["ldlr_sp1_peak"]["log2"]
    sp1_loss = 100.0 * (1.0 - 2.0 ** sp1_v)
    axC2.annotate(
        f"Repeat 2 SRE-1 (G4)\n{sre_v:.2f} log₂ ({sre_loss:.0f}% loss)",
        xy=(sre_p, sre_v), xytext=(55, -2.68),
        arrowprops=dict(arrowstyle="->", color=COL_AMBER, lw=0.9),
        fontsize=6.7, color=COL_AMBER, fontweight="bold",
    )
    axC2.annotate(
        f"Repeat 3 SP1 (G5)\n{sp1_v:.2f} log₂ ({sp1_loss:.0f}% loss)",
        xy=(sp1_p, sp1_v), xytext=(216, -3.38),
        arrowprops=dict(arrowstyle="->", color=COL_CRIMSON, lw=0.9),
        fontsize=6.7, color=COL_CRIMSON, fontweight="bold",
    )
    axC2.text(
        215, 0.56, f"LDLR vs LDLR.2\nr = {r4['ldlr_vs_ldlr2_snv_pearson_r']:.3f}",
        fontsize=6.7, color=COL_NAVY,
        bbox=dict(boxstyle="round,pad=0.2", fc="#EBF5FB", ec="#AEB6BF", lw=0.6),
    )
    axC2.set_xlim(45, 345)
    axC2.set_ylim(-3.85, 1.05)
    axC2.set_xlabel("LDLR promoter position (bp)")
    axC2.set_ylabel("Δlog₂(RNA/DNA) for 8-oxoG→T")
    axC2.set_title("C2. LDLR Base-Resolution Track", loc="left", pad=6)
    axC2.legend(loc="upper left", bbox_to_anchor=(-0.06, -0.18), ncol=3, frameon=False, fontsize=6.6)

    # =========================================================================
    # PANEL D (gs_outer[1, 1]): Row 5a GO Uncollapsing + Row 5b Beamline Spec
    # =========================================================================
    gs_D = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[1, 1], width_ratios=[1.15, 1.05], wspace=0.33)

    # --- D1: Row 5a E. coli WGS MA G:C -> T:A Uncollapsing ---
    axD1 = fig.add_subplot(gs_D[0, 0])
    class_order = [
        ("MMR_deficient", "MMR⁻ (mutS/L/H)"),
        ("proofreading_deficient", "Pol III⁻ (dnaQ)"),
        ("wild_type", "Wild-type (3 str.)"),
        ("pyrimidine_BER_deficient", "Pyr BER⁻ (nth nei)"),
        ("AP_endonuclease_deficient", "AP endo⁻ (xthA nfo)"),
        ("GO_glycosylase_deficient", "GO⁻ (mutY)"),
        ("GO_double_deficient", "GO²⁻ (mutM mutY)"),
    ]
    y_labels = []
    pct_vals = []
    bar_colors = []
    for rclass, lbl in class_order:
        sub = ecoli_df[ecoli_df["repair_class"] == rclass]
        frac = 100.0 * sub["GC_to_TA_count"].sum() / max(sub["total_snvs"].sum(), 1)
        y_labels.append(lbl)
        pct_vals.append(frac)
        bar_colors.append(COL_CRIMSON if "GO" in rclass else (COL_NAVY if rclass == "wild_type" else COL_SLATE))

    y_idx = np.arange(len(y_labels))
    axD1.barh(y_idx, pct_vals, color=bar_colors, edgecolor="#2C3E50", lw=0.65, height=0.62)
    for yi, pv in enumerate(pct_vals):
        axD1.text(pv + 1.8, yi, f"{pv:.1f}%", va="center", ha="left", fontsize=6.9, fontweight="bold")

    r5 = m6["row5_ecoli_go_and_beamline"]
    gr = r5["go_deficient_G_run_normalized"]
    axD1.text(
        37.0, 2.35,
        f"GO⁻ G-run per-bp rate:\nG1: {gr['G1']['per_bp_relative_to_G1']:.2f}× → G3: {gr['G3']['per_bp_relative_to_G1']:.2f}×\nG4: {gr['G4']['per_bp_relative_to_G1']:.2f}× → G5+: {gr['G5+']['per_bp_relative_to_G1']:.2f}×\n({r5['mutMmutY_vs_MMR_fraction_ratio']:.1f}× vs MMR⁻)",
        fontsize=6.6, color=COL_CRIMSON, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.22", fc="#FDEDEC", ec=COL_CRIMSON, lw=0.65),
    )
    axD1.set_yticks(y_idx)
    axD1.set_yticklabels(y_labels, fontsize=7.0)
    axD1.set_xlim(0, 119)
    axD1.set_xlabel("Whole-genome G:C→T:A fraction (%)")
    axD1.set_title("D1. Row 5a: GO Quotient Uncollapsing", loc="left", pad=6)

    # --- D2: Row 5b Sole Unexecuted Polarized Lepton Beamline Budget ---
    axD2 = fig.add_subplot(gs_D[0, 1])
    asym_grid = np.logspace(-4, -2, 200)
    n_3sig = 9.0 / (asym_grid ** 2)
    n_5sig = 25.0 / (asym_grid ** 2)

    bm = r5["row5b_beamline_protocol"]
    a_tgt = float(bm["target_asymmetry_threshold"])
    n_min_bm = float(bm["min_lesion_events_per_helicity_arm_3sigma"])
    t_arm_s = float(bm["estimated_irradiation_seconds_per_arm"])

    axD2.loglog(asym_grid, n_3sig, color=COL_NAVY, lw=2.1, label="3σ helicity-flip threshold")
    axD2.loglog(asym_grid, n_5sig, color=COL_PURPLE, lw=1.8, ls="--", label="5σ discovery threshold")
    axD2.axvspan(1.0e-3, 1.88e-3, color="#FCF3CF", alpha=0.85, label="Mod. 5 B-DNA window")
    axD2.scatter([a_tgt], [n_min_bm], color=COL_CRIMSON, s=46, zorder=6)
    axD2.annotate(
        f"Row 5b target:\n{n_min_bm/1e6:.1f}×10⁶ lesions/arm\n(~{t_arm_s:.1f} s at {bm['standard_GaAs_e_gun_current_nA']:.0f} nA e⁻,\nψ = {bm['recommended_tilt_pair_deg'][0]:.0f}° vs {bm['recommended_tilt_pair_deg'][1]:.0f}°)",
        xy=(a_tgt, n_min_bm), xytext=(1.15e-4, 1.5e5),
        arrowprops=dict(arrowstyle="->", color=COL_CRIMSON, lw=1.0),
        fontsize=6.7, color=COL_CRIMSON, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.22", fc="white", ec=COL_CRIMSON, lw=0.65),
    )
    axD2.set_xlim(8e-5, 1.2e-2)
    axD2.set_ylim(5e4, 4e9)
    axD2.set_xlabel("Chiral yield asymmetry |A_eff|")
    axD2.set_ylabel("Required lesions per arm N_min")
    axD2.set_title("D2. Row 5b: Polarized Beam Spec", loc="left", pad=6)
    axD2.legend(loc="upper left", bbox_to_anchor=(-0.06, -0.18), ncol=2, frameon=False, fontsize=6.5)

    png_path = f"{ASSET_DIR}/fig6_table6_empirical_wetlab_validation.png"
    svg_path = f"{ASSET_DIR}/fig6_table6_empirical_wetlab_validation.svg"
    fig.savefig(png_path, dpi=300)
    fig.savefig(svg_path)
    plt.close(fig)
    print(f"Saved {png_path} ({os.path.getsize(png_path)} bytes) and {svg_path} ({os.path.getsize(svg_path)} bytes)")


if __name__ == "__main__":
    make_figure6()
