#!/usr/bin/env python3
"""Generate publication-grade multi-panel figures (PNG 300 dpi + SVG)
for the Observer-Bounded DNA Access & Rate-Distortion manuscript.
All panels read directly from verified CSV/JSON outputs in data/.
"""

import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import stats

from analysis_paths import ASSET_DIR, DATA_DIR, REPO_ROOT

os.makedirs(ASSET_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.2,
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
COL_LIGHT = "#D5DBDB"


def save_pair(fig, stem):
    png_path = os.path.join(ASSET_DIR, f"{stem}.png")
    svg_path = os.path.join(ASSET_DIR, f"{stem}.svg")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {png_path} ({os.path.getsize(png_path)} bytes) and {svg_path} ({os.path.getsize(svg_path)} bytes)")


def make_figure1():
    res_df = pd.read_csv(f"{DATA_DIR}/module1_resolvent_convergence.csv")
    wf_df = pd.read_csv(f"{DATA_DIR}/module1_weak_field_nogo.csv")
    trap_df = pd.read_csv(f"{DATA_DIR}/module1_guanine_hole_trapping.csv")
    with open(f"{DATA_DIR}/module1_summary.json", encoding="utf-8") as f:
        m1 = json.load(f)

    fig = plt.figure(figsize=(12.2, 8.8), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, hspace=0.14, wspace=0.12)

    # Panel A: Resolvent semigroup error & open vs unitary coherence decay
    axA = fig.add_subplot(gs[0, 0])
    t_ps = res_df["time_s"].values * 1e12
    err = np.maximum(res_df["op_norm_error"].values, 1e-16)
    col_open = "coherence_norm_open" if "coherence_norm_open" in res_df.columns else "open_coherence_l2"
    col_unit = "coherence_norm_unitary" if "coherence_norm_unitary" in res_df.columns else "closed_coherence_l2"
    coh_open = np.maximum(res_df[col_open].values, 1e-16)
    coh_unit = np.maximum(res_df[col_unit].values, 1e-16)
    axA.loglog(t_ps, err, color=COL_NAVY, lw=2.2, label=r"Resolvent error $\|\mathcal{Q}e^{t\mathcal{L}}\mathcal{Q} - e^{t\mathcal{L}_{\mathrm{eff}}}\mathcal{Q}\|_2$")
    axA.loglog(t_ps, coh_open, color=COL_PURPLE, ls="-.", lw=1.5, alpha=0.85, label="Open-system quantum coherence")
    axA.loglog(t_ps, coh_unit, color=COL_AMBER, ls=":", lw=1.3, alpha=0.75, label="Unitary (bath-free) coherence")
    axA.axvline(m1["tau_mem_fs"] * 1e-3, color=COL_CRIMSON, ls="--", lw=1.3, label=rf"Bath memory $\tau_{{\mathrm{{mem}}}} = {m1['tau_mem_fs']:.1f}\ \mathrm{{fs}}$")
    axA.axvline(1.0, color=COL_TEAL, ls=":", lw=1.3, label=r"Picosecond kinetic gate ($1\ \mathrm{ps}$)")
    axA.set_xlabel("Observation interval $t$ (ps)")
    axA.set_ylabel("Operator-norm projection error")
    axA.set_title("A   Mori-Zwanzig Schur Resolvent Convergence", loc="left")
    axA.set_ylim(1e-16, 2.0)
    axA.grid(True, which="major", alpha=0.22, lw=0.6)
    axA.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.92)

    # Panel B: Accessible Moment Annihilation (m = 0..5)
    axB = fig.add_subplot(gs[0, 1])
    moments = m1["moment_annihilation_norms"]
    ms = np.arange(6)
    raw_vals = np.array([float(moments[f"m={m}"]) for m in ms])
    # Normalize m>=3 relative to spectral gap powers g^m so physical dimensionless contrast is plotted cleanly alongside exact zeros!
    g_gap = float(m1["spectral_gap_s_inv"])
    dimless = np.array([raw_vals[m] / (g_gap ** m) if raw_vals[m] > 0 else 1e-22 for m in ms])
    colors_B = [COL_TEAL if raw_vals[m] == 0.0 else COL_CRIMSON for m in ms]
    bars = axB.bar(ms, dimless, color=colors_B, width=0.58, edgecolor="#2C3E50", lw=0.8)
    axB.set_yscale("log")
    axB.set_ylim(1e-24, 1e4)
    axB.axhline(1e-16, color=COL_SLATE, ls="--", lw=1.0, label="IEEE-754 machine-precision floor ($10^{-16}$)")
    for m_idx in [0, 1, 2]:
        axB.text(m_idx, 3e-22, "Exact 0", ha="center", va="bottom", fontsize=7.8, fontweight="bold", color=COL_TEAL)
    axB.annotate(
        r"GKSL pure-dephasing cubic onset: $\mathcal{L}_{QP}\,\delta\mathcal{L}_{PP}\,\mathcal{L}_{PQ} \neq 0$",
        xy=(2.72, dimless[3] * 0.35), xytext=(0.02, 4e-4),
        arrowprops=dict(arrowstyle="->", color=COL_CRIMSON, lw=1.15),
        fontsize=7.7, color=COL_CRIMSON, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor=COL_CRIMSON, lw=0.6, alpha=0.94)
    )
    axB.set_xticks(ms)
    axB.set_xticklabels([f"$m={m}$" for m in ms])
    axB.set_xlabel("Short-time accessible derivative order $m$")
    axB.set_ylabel(r"Scaled moment shift $\|\Delta M_m\|_2 / g(\mathcal{L}_{PP})^m$")
    axB.set_title("B   Accessible Moment Annihilation Theorem", loc="left")
    axB.grid(True, axis="y", alpha=0.22, lw=0.6)
    axB.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.92)

    # Panel C: Universal Sub-Thermal Coherent Wave No-Go Bound
    axC = fig.add_subplot(gs[1, 0])
    E_grid = wf_df["E0_V_per_m"].values
    axC.loglog(E_grid, np.maximum(wf_df["thermo_ceiling_dmu_over_mu"].values, 1e-36), color=COL_NAVY, lw=2.0, label=r"Static Kubo ceiling $2\lambda / (k_B T)$ ($[\Pi,H_0]=0$)")
    axC.loglog(E_grid, np.maximum(wf_df["dc_gibbs_dmu_over_mu"].values, 1e-36), color=COL_TEAL, ls="--", lw=1.7, label=r"Exact DC Gibbs shift ($\Pi_{\varepsilon_1(0)}$)")
    axC.loglog(E_grid, np.maximum(wf_df["drazin_floquet_subgap_dmu_over_mu"].values, 1e-36), color=COL_PURPLE, lw=2.1,
               label=r"Fixed-bath Floquet Drazin ($\Pi_{\varepsilon_1(0)}$, $\omega \ll g$)")
    axC.loglog(E_grid, np.maximum(wf_df["drazin_floquet_thz_dmu_over_mu"].values, 1e-36), color=COL_AMBER, ls=":", lw=1.8,
               label=r"Near-gap THz Floquet ($\Pi_{\varepsilon_1(0)}$, $10^{13}\,\mathrm{rad/s}$)")
    if "adiabatic_kms_fixed_proj_dmu_over_mu" in wf_df.columns:
        axC.loglog(E_grid, np.maximum(wf_df["adiabatic_kms_fixed_proj_dmu_over_mu"].values, 1e-36), color="#1F618D", ls=(0, (3, 1, 1, 1)), lw=1.6,
                   label=r"Adiabatic KMS ($\Pi_{\varepsilon_1(0)}$ fixed eigenproj.)")
    axC.loglog(E_grid, np.maximum(wf_df["adiabatic_kms_ac_dmu_over_mu"].values, 1e-36), color="#2874A6", ls="-.", lw=1.7,
               label=r"Adiabatic KMS ($|\varepsilon_1(t)\rangle$ inst. eigenstate)")
    if "adiabatic_kms_bare_tautomer_dmu_over_mu" in wf_df.columns:
        axC.loglog(E_grid, np.maximum(wf_df["adiabatic_kms_bare_tautomer_dmu_over_mu"].values, 1e-36), color=COL_SLATE, ls=(0, (1, 1.5)), lw=1.4,
                   label=r"Cycle-avg KMS ($|\mathrm{G}^*\text{:}\mathrm{C}^*\rangle$ bare taut.)")
    m50 = float(m1["metabolic_50mK_noise_floor"])
    m50_exp = int(np.floor(np.log10(m50)))
    m50_man = m50 / (10.0 ** m50_exp)
    axC.axhline(m50, color=COL_CRIMSON, ls=(0, (5, 2)), lw=1.6,
                label=rf"50 mK thermal floor (${m50_man:.2f}\times 10^{{{m50_exp}}}$)")
    base_pt = m1["ambient_subgap_1uVm_reference"]
    fl1 = float(base_pt["drazin_floquet_subgap_dmu_over_mu"])
    fl1_exp = int(np.floor(np.log10(fl1)))
    fl1_man = fl1 / (10.0 ** fl1_exp)
    kms1 = float(base_pt["adiabatic_kms_ac_dmu_over_mu"])
    kms1_exp = int(np.floor(np.log10(kms1)))
    kms1_man = kms1 / (10.0 ** kms1_exp)
    axC.scatter([1e-6], [fl1], color=COL_PURPLE, s=54, zorder=6, edgecolor="black", lw=0.8,
                label=rf"$1\ \mu\mathrm{{V/m}}$ Floquet $\Pi_{{\varepsilon_1(0)}}$ (${fl1_man:.2f}\times 10^{{{fl1_exp}}}$)")
    axC.scatter([1e-6], [kms1], color="#2874A6", marker="D", s=44, zorder=6, edgecolor="black", lw=0.7,
                label=rf"$1\ \mu\mathrm{{V/m}}$ KMS $|\varepsilon_1(t)\rangle$ (${kms1_man:.2f}\times 10^{{{kms1_exp}}}$)")
    axC.set_xlabel(r"Non-ionizing electric field amplitude $E_0$ ($\mathrm{V\,m^{-1}}$)")
    axC.set_ylabel(r"Fractional mutation-hazard shift $|\Delta\mu / \mu|$")
    axC.set_title("C   Universal Sub-Thermal Wave No-Go Bound", loc="left")
    axC.set_ylim(1e-35, 1e3)
    axC.grid(True, which="major", alpha=0.22, lw=0.6)
    axC.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.92, fontsize=6.0)

    # Panel D: Duplex Guanine Stack Hole-Sink Enrichment
    axD = fig.add_subplot(gs[1, 1])
    gruns = trap_df["g_run_length"].values
    enrich = trap_df["fold_enrichment_per_G_vs_G1"].values
    homo_mev = -trap_df["homo_sink_eV"].values * 1000.0
    x_pos = np.arange(len(gruns))
    axD.bar(x_pos - 0.20, enrich, width=0.36, color=COL_NAVY, edgecolor="#1B2631", lw=0.7, label="Per-G 8-oxoG fold enrichment vs. $G_1$")
    axD.set_ylabel("Per-guanine 8-oxoG fold enrichment", color=COL_NAVY)
    axD.tick_params(axis="y", labelcolor=COL_NAVY)
    axD.set_ylim(0, 24.5)
    for i, val in enumerate(enrich):
        axD.text(x_pos[i] - 0.20, val + 0.55, f"{val:.1f}x", ha="center", va="bottom", fontsize=7.7, color=COL_NAVY, fontweight="bold")

    axD2 = axD.twinx()
    axD2.spines["top"].set_visible(False)
    axD2.plot(x_pos + 0.22, homo_mev, color=COL_AMBER, marker="o", lw=2.0, ms=6, label=r"HOMO stabilization $-\Delta E_{\mathrm{HOMO}}$ (meV)")
    axD2.set_ylabel(r"HOMO hole-sink depth $-\Delta E_{\mathrm{HOMO}}$ (meV)", color=COL_AMBER)
    axD2.tick_params(axis="y", labelcolor=COL_AMBER)
    axD2.set_ylim(0, 860)

    axD.set_xticks(x_pos)
    axD.set_xticklabels([f"$G_{g}$" if g > 1 else "$G_1$ (isolated)" for g in gruns])
    axD.set_xlabel("Contiguous duplex guanine tract length")
    axD.set_title("D   Duplex Guanine Hole Migration & 8-oxoG Trapping", loc="left")
    h1, l1 = axD.get_legend_handles_labels()
    h2, l2 = axD2.get_legend_handles_labels()
    axD.legend(h1 + h2, l1 + l2, loc="upper left", frameon=True, facecolor="white", framealpha=0.92, fontsize=7.8)

    save_pair(fig, "fig1_resolvent_and_nogo")


def make_figure2():
    sp_df = pd.read_csv(f"{DATA_DIR}/module2_single_cell_sparsity.csv")
    rd_df = pd.read_csv(f"{DATA_DIR}/module2_ensemble_rate_distortion.csv")
    cx_df = pd.read_csv(f"{DATA_DIR}/module2_memory_sensing_crossover.csv")

    fig = plt.figure(figsize=(12.2, 8.8), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, hspace=0.14, wspace=0.12)

    # Panel A: Single-cell erasure vs ensemble detection probability
    axA = fig.add_subplot(gs[0, 0])
    nu_dense = np.logspace(0, 6, 200)
    styles = [(COL_TEAL, "-"), (COL_NAVY, "--"), (COL_CRIMSON, "-.")]
    for (_, r), (col, ls) in zip(sp_df.iterrows(), styles):
        lam1 = float(r["mean_lesions_per_promoter_per_cell"])
        base_r = float(r["baseline_8oxoG_per_G"])
        b_exp = int(round(np.log10(base_r)))
        b_man = base_r / (10.0 ** b_exp)
        rate_str = rf"10^{{{b_exp}}}" if abs(b_man - 1.0) < 1e-6 else rf"{b_man:.0f}\times 10^{{{b_exp}}}"
        lbl = rf"${rate_str}$ / G ($\lambda_1={lam1:.4f}$)"
        p_hit = 1.0 - np.exp(-nu_dense * lam1)
        axA.semilogx(nu_dense, p_hit, color=col, ls=ls, lw=2.1, label=lbl)
    p0_max_pct = 100.0 * float(sp_df["prob_zero_lesions_single_cell"].min())
    axA.axvspan(1, 50, color=COL_LIGHT, alpha=0.45, label="Single-cell / oligocellular erasure regime")
    axA.axhline(0.95, color=COL_SLATE, ls=":", lw=1.1)
    axA.text(1.5, 0.90, rf"$P(N_{{\mathrm{{lesion}}}}=0) \geq {p0_max_pct:.1f}\%$ at $\nu=1$", fontsize=8.3, fontweight="bold", color="#2C3E50")
    axA.set_xlabel(r"Synchronized cell ensemble size $\nu$")
    axA.set_ylabel(r"Active promoter occupancy $P(N_{\mathrm{lesion}} \geq 1)$")
    axA.set_ylim(-0.03, 1.06)
    axA.set_title("A   Single-Cell Poisson Erasure vs. Ensemble Pooling", loc="left")
    axA.grid(True, which="major", alpha=0.22, lw=0.6)
    axA.legend(loc="center right", frameon=True, facecolor="white", framealpha=0.93, fontsize=7.8)

    # Panel B: Normalized waveform distortion and the Monte Carlo van Trees estimate.
    axB = fig.add_subplot(gs[0, 1])
    information_col = (
        "gaussian_fisher_information_proxy_bits"
        if "gaussian_fisher_information_proxy_bits" in rd_df else "mutual_info_bits"
    )
    rd_1e5 = rd_df[rd_df["nu_cells"] == 100000].set_index("schedule")
    mi_rhy = float(rd_1e5.loc["rhythmic_OGG1", information_col])
    mi_flt = float(rd_1e5.loc["flat_equal_OGG1", information_col])
    sched_cfg = [
        ("rhythmic_OGG1", COL_NAVY, "o", "Open SHL (-6), circadian OGG1 (8VWT)"),
        ("flat_equal_OGG1", COL_TEAL, "s", "Open SHL (-6), flat equal-mean OGG1"),
        ("inward_occluded_SHL", COL_CRIMSON, "^", r"Inward-occluded minor groove (1KX5 $\geq 120^\circ$, $18\times$ gate)"),
    ]
    for skey, col, mk, lbl in sched_cfg:
        sub = rd_df[rd_df["schedule"] == skey].sort_values("nu_cells")
        nu = sub["nu_cells"].values
        emp = sub["normalized_distortion_D_over_sigma2"].values
        se = sub["empirical_mse_se"].values / 0.0324
        bcrb = sub["bcrb_D_over_sigma2"].values
        axB.loglog(nu, bcrb, color=col, ls="--", lw=1.5, alpha=0.85)
        axB.errorbar(nu, emp, yerr=se, fmt=mk, color=col, ms=5.5, capsize=2.5, lw=1.2, label=lbl)
    axB.plot([], [], color=COL_SLATE, ls="--", lw=1.5, label="Van Trees bound (Monte Carlo Fisher estimate)")
    axB.set_xlabel(r"Synchronized hepatocyte ensemble size $\nu$")
    axB.set_ylabel(r"Normalized waveform MSE $D / \sigma^2$")
    axB.set_ylim(7e-4, 2.2)
    axB.set_title("B   Circadian Waveform Reconstruction & SHL Occlusion", loc="left")
    axB.grid(True, which="major", alpha=0.22, lw=0.6)
    axB.legend(loc="lower left", frameon=True, facecolor="white", framealpha=0.93, fontsize=7.4)

    # Panel C: Gaussian Fisher approximation to information scaling.
    axC = fig.add_subplot(gs[1, 0])
    c_labels = {
        "rhythmic_OGG1": "Open SHL (-6), circadian OGG1 (8VWT)",
        "flat_equal_OGG1": "Open SHL (-6), flat equal-mean OGG1",
        "inward_occluded_SHL": r"Inward-occluded minor groove (1KX5 $\geq 120^\circ$)",
    }
    for skey, col, mk, _ in sched_cfg:
        sub = rd_df[rd_df["schedule"] == skey].sort_values("nu_cells")
        axC.semilogx(sub["nu_cells"].values, sub[information_col].values, color=col, marker=mk, lw=2.0, ms=5.5, label=c_labels[skey])
    axC.annotate(f"Rhythmic and flat OGG1\n({mi_rhy:.2f} and {mi_flt:.2f} bits at $\\nu=10^5$)",
                 xy=(1e5, mi_rhy), xytext=(8e1, 21.0),
                 arrowprops=dict(arrowstyle="->", color=COL_NAVY, lw=1.1),
                 fontsize=8.0, color=COL_NAVY)
    axC.set_xlabel(r"Synchronized hepatocyte ensemble size $\nu$")
    axC.set_ylabel("Gaussian Fisher approximation (bits)")
    axC.set_ylim(-1.0, 32.5)
    axC.set_title("C   Information Scaling in the Poisson Model", loc="left")
    axC.grid(True, which="major", alpha=0.22, lw=0.6)
    axC.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.93, fontsize=7.8)

    # Panel D: Two-stage Sensing vs Classical-Memory Rate-Distortion Crossover
    axD = fig.add_subplot(gs[1, 1])
    cb_palette = [(2.0, COL_CRIMSON), (4.0, COL_AMBER), (8.0, COL_PURPLE), (16.0, COL_TEAL)]
    sub0 = cx_df[cx_df["C_B_bits"] == 16.0].sort_values("K_modes")
    axD.semilogy(sub0["K_modes"].values, sub0["D_sensing_norm"].values, color="#1C2833", ls="--", lw=2.0, label=r"Front-end sensing floor $D_{\mathrm{sens}}(\nu=10^5, K)$")
    for cb, col in cb_palette:
        sub = cx_df[cx_df["C_B_bits"] == cb].sort_values("K_modes")
        axD.semilogy(sub["K_modes"].values, sub["D_effective_lower_bound_norm"].values, color=col, marker="o", ms=3.8, lw=1.8, label=rf"Register bound $C_B = {int(cb)}$ bits")
    axD.set_xlabel("Decoded harmonic dimension $K$ (Fourier modes)")
    axD.set_ylabel(r"Two-stage lower bound on $D / \sigma^2$")
    axD.set_xlim(0.5, 24.5)
    axD.set_ylim(7e-3, 1.4)
    axD.set_title("D   Sensing-Limited vs. Memory-Limited Crossover", loc="left")
    axD.grid(True, which="major", alpha=0.22, lw=0.6)
    axD.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.93, fontsize=7.7)

    save_pair(fig, "fig2_circadian_rate_distortion")


def make_figure3():
    sbs_df = pd.read_csv(f"{DATA_DIR}/module3_ecoli_sbs96_spectra.csv")
    pop_df = pd.read_csv(f"{DATA_DIR}/module3_ecoli_populations.csv")
    ts_df = pd.read_csv(f"{DATA_DIR}/module3_ecoli_mmr_stacking_regression.csv")
    with open(f"{DATA_DIR}/module3_summary.json", encoding="utf-8") as f:
        m3 = json.load(f)

    fig = plt.figure(figsize=(13.0, 9.6), constrained_layout=True)
    gs = GridSpec(4, 2, figure=fig, height_ratios=[1.0, 1.0, 1.0, 1.35], hspace=0.18, wspace=0.14)

    sub_Order = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    sub_Colors = {
        "C>A": "#03BCE4",
        "C>G": "#1C2833",
        "C>T": "#E42926",
        "T>A": "#999999",
        "T>C": "#9FCE63",
        "T>G": "#ECC6C5",
    }
    x96 = np.arange(96)
    bar_colors = [sub_Colors[st] for st in sbs_df["sub_type"].values]

    spectra_panels = [
        (0, "WT_repair_proficient_opp_norm_prob",
         f"A   Repair-Proficient WT Quotient (n = {m3['wt_independent_lineages']} lineages, {m3['wt_total_unique_snps']} SNPs; H = {m3['shannon_entropy_bits_WT']:.2f} bits)"),
        (1, "MMR_deficient_opp_norm_prob",
         f"B   MMR-Deficient mutS/mutL Uncollapsed Spectrum (n = {m3['mmr_independent_lineages']} lineages, {m3['mmr_total_unique_snps']} SNPs; H = {m3['shannon_entropy_bits_MMR']:.2f} bits)"),
        (2, "mutT_8oxodGTP_mutator_opp_norm_prob",
         f"C   8-oxo-dGTPase mutT-Deficient Spectrum (n = {len(m3['mutT_populations'])} lineages, {m3['mutT_total_unique_snps']} SNPs; H = {m3['shannon_entropy_bits_mutT']:.2f} bits)"),
    ]

    for row_idx, col_name, title in spectra_panels:
        ax = fig.add_subplot(gs[row_idx, :])
        vals = sbs_df[col_name].values
        ax.bar(x96, vals, color=bar_colors, width=0.82, rasterized=True)
        ymax = max(np.max(vals) * 1.24, 0.04)
        ax.set_ylim(0, ymax)
        ax.set_xlim(-0.8, 95.8)
        for b_idx, st in enumerate(sub_Order):
            x0, x1 = b_idx * 16 - 0.4, (b_idx + 1) * 16 - 0.6
            if b_idx > 0:
                ax.axvline(b_idx * 16 - 0.5, color="#BFC9CA", lw=0.7, ls=":")
            if row_idx == 0:
                ax.text(b_idx * 16 + 7.5, ymax * 0.87, st, ha="center", va="center",
                        fontsize=8.6, fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.18", facecolor=sub_Colors[st], alpha=0.22, edgecolor="none"))
        ax.set_ylabel("Opp.-norm. prob.")
        ax.set_title(title, loc="left", fontsize=9.8)
        if row_idx < 2:
            ax.set_xticks([])
        else:
            tick_idx = np.arange(0, 96, 4)
            ax.set_xticks(tick_idx)
            ax.set_xticklabels([sbs_df["sbs96_channel"].values[i][0] + sbs_df["sbs96_channel"].values[i][2] + sbs_df["sbs96_channel"].values[i][6] for i in tick_idx],
                               rotation=90, fontsize=7.2)
            ax.set_xlabel("Trinucleotide opportunity channel (5'-base, pyrimidine reference, 3'-base)")

    # Panel D: Lineage-level Orthogonal Separation across 27 E. coli Populations
    axD = fig.add_subplot(gs[3, 0])
    reg_cfg = [
        ("WT_repair_proficient", COL_NAVY, "o", f"WT non-mutator (n={m3['wt_independent_lineages']})"),
        ("IS_mutator", COL_SLATE, "D", "IS-mutator control (Ara+1, n=1)"),
        ("MMR_deficient", COL_CRIMSON, "s", f"MMR-deficient mutS/mutL (n={m3['mmr_independent_lineages']})"),
        ("mutT_8oxodGTP_mutator", COL_AMBER, "^", f"mutT 8-oxo-dGTPase (n={len(m3['mutT_populations'])})"),
    ]
    for rkey, col, mk, lbl in reg_cfg:
        sub = pop_df[pop_df["regime"] == rkey]
        if len(sub) == 0:
            continue
        axD.scatter(sub["transition_fraction"].values * 100.0,
                    sub["AT_to_CG_fraction"].values * 100.0,
                    color=col, marker=mk, s=58, alpha=0.88, edgecolor="black", lw=0.6, label=lbl)
    axD.annotate(
        f"Exact Mann-Whitney p = {m3['p_mwu_exact_lineage_transition_WT_vs_MMR']:.2e}\n"
        f"cos(MMR, mutT) dist = {m3['cos_dist_MMR_vs_mutT']:.3f}",
        xy=(95.5, 2.0), xytext=(36.0, 38.0),
        arrowprops=dict(arrowstyle="->", color=COL_CRIMSON, lw=1.1),
        fontsize=7.9, color="#1C2833",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#FDEDEC", edgecolor=COL_CRIMSON, lw=0.7)
    )
    axD.set_xlabel("Lineage C>T + T>C transition fraction (%)")
    axD.set_ylabel("Lineage A:T > C:G transversion fraction (%)")
    axD.set_xlim(-4, 104)
    axD.set_ylim(-6, 106)
    axD.set_title("D   Orthogonal Repair-Pathway Separation Across 27 Lineages", loc="left", fontsize=9.8)
    axD.grid(True, alpha=0.22, lw=0.6)
    axD.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.92, fontsize=7.5)

    # Panel E: MMR-deficient C>T mutation rate vs nearest-neighbor stacking stability
    axE = fig.add_subplot(gs[3, 1])
    ct_sub = ts_df[ts_df["sub_type"] == "C>T"].copy()
    x_stack = -ct_sub["dG_stack_kcal_mol"].values
    y_rate = ct_sub["posterior_rate_per_opp"].values * 1e3
    gc_fl = ct_sub["gc_flanks"].values
    for g_val, col, lbl in [(0, COL_TEAL, "0 G/C flanks"), (1, COL_NAVY, "1 G/C flank"), (2, COL_CRIMSON, "2 G/C flanks")]:
        msk = gc_fl == g_val
        axE.scatter(x_stack[msk], y_rate[msk], color=col, s=62, edgecolor="black", lw=0.6, label=lbl, zorder=4)
    slope, intercept = np.polyfit(x_stack, np.log10(y_rate), 1)
    x_line = np.linspace(np.min(x_stack) - 0.12, np.max(x_stack) + 0.12, 100)
    axE.plot(x_line, 10 ** (intercept + slope * x_line), color="#2C3E50", ls="--", lw=1.6,
             label=f"Spearman r = {m3['stacking_spearman_CT_r']:.3f} (p = {m3['stacking_spearman_CT_p']:.2e})")
    axE.set_yscale("log")
    top2 = ct_sub.nlargest(2, "posterior_rate_per_opp")
    offsets_e = [(-0.46, 1.22), (-0.48, 1.28)]
    for (_, r), (dx, fy) in zip(top2.iterrows(), offsets_e):
        axE.annotate(r["channel"], xy=(-r["dG_stack_kcal_mol"], r["posterior_rate_per_opp"] * 1e3),
                     xytext=(-r["dG_stack_kcal_mol"] + dx, r["posterior_rate_per_opp"] * 1e3 * fy),
                     fontsize=7.6, color=COL_CRIMSON, fontweight="bold")
    axE.set_xlabel(r"Flanking base-stacking stability $-\Delta G^\circ_{37,\mathrm{stack}}$ ($\mathrm{kcal\,mol^{-1}}$)")
    axE.set_ylabel(r"MMR-deficient C>T rate ($\times 10^{-3}$ per site)")
    axE.set_ylim(np.min(y_rate) * 0.55, np.max(y_rate) * 2.35)
    axE.set_title("E   Active-Site Stacking Stabilization of Unrepaired C>T", loc="left", fontsize=9.8)
    axE.grid(True, which="major", alpha=0.22, lw=0.6)
    axE.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.92, fontsize=7.6)

    save_pair(fig, "fig3_ecoli_quotient_uncollapsing")


def make_figure4():
    comp_df = pd.read_csv(f"{DATA_DIR}/module3_human_promoter_comparison.csv")
    ctl_df = pd.read_csv(f"{DATA_DIR}/module3_human_matched_1500_controls.csv")
    clk_df = pd.read_csv(f"{DATA_DIR}/module3_human_circadian_10_promoters.csv")
    ism_dnase = pd.read_csv(f"{DATA_DIR}/module4_cry1_ism_dnase_liver.csv")
    trk_dnase = pd.read_csv(f"{DATA_DIR}/module4_cry1_ref_alt_dnase_tracks.csv")
    ag_df = pd.read_csv(f"{DATA_DIR}/module4_alphagenome_circadian_variants.csv")
    with open(f"{DATA_DIR}/module4_summary.json", encoding="utf-8") as f:
        m4 = json.load(f)

    fig = plt.figure(figsize=(12.6, 9.2), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, hspace=0.16, wspace=0.14)

    # Panel A: Propensity-matched forest plot (10 clock promoters vs 1,500 matched CGI controls vs 4,560 10-chr MANE background)
    axA = fig.add_subplot(gs[0, 0])
    label_map = {
        "gc_frac": "GC fraction (matched)",
        "cpg_oe": "CpG O/E ratio (matched)",
        "max_g_run": "Max G-run length (matched)",
        "g4hunter_max": "G4Hunter max (matched)",
        "g3plus_tracts": "Contiguous G>=3 tracts",
        "hole_trap_burden": "8-oxoG hole-trap burden",
        "mean_phi_8oxoG_per_G": "Mean 8-oxoG yield per G",
        "psd_10_4bp_dinuc_snr": "10.4-bp dinucleotide PSD",
        "psd_10_4bp_gc_snr": "10.4-bp G/C rotational PSD",
        "sp1_GGGCGG_count": "SP1 GGGCGG motif count",
        "ebox_CACGTG_count": "Canonical E-box CACGTG",
    }
    sub_c = comp_df.iloc[::-1].reset_index(drop=True)
    y_pos = np.arange(len(sub_c))
    d_matched = sub_c["delta_std_vs_matched1500"].values
    n_clk = float(len(clk_df))
    n_ctl = float(len(ctl_df))
    se_hedges = np.sqrt(1.0 / n_clk + 1.0 / n_ctl + (d_matched ** 2) / (2.0 * (n_clk + n_ctl)))
    d_unadj = (sub_c["clock_mean"].values - sub_c["all4560_bg_mean"].values) / np.maximum(sub_c["clock_sd"].values, 1e-6)
    min_bh_q = float(sub_c["bh_fdr_q"].dropna().min())

    axA.axvline(0.0, color="#1C2833", lw=1.0)
    axA.scatter(d_unadj, y_pos + 0.14, color=COL_AMBER, marker="D", s=38, alpha=0.88, label="Unadjusted vs. 4,560 10-chr MANE promoters")
    axA.errorbar(d_matched, y_pos - 0.14, xerr=1.96 * se_hedges, fmt="o", color=COL_NAVY, ms=5.2, lw=1.3, capsize=2.5,
                 label=f"Mahalanobis-matched vs. {int(n_ctl):,} CGI controls (all BH q ≥ {min_bh_q:.2f})")
    axA.axhline(6.5, color=COL_SLATE, ls=":", lw=0.9)
    axA.set_yticks(y_pos)
    axA.set_yticklabels([label_map[f] for f in sub_c["feature"].values], fontsize=8.0)
    axA.set_xlabel(r"Standardized difference $\Delta_{\mathrm{std}}$ (95% CI)")
    axA.set_xlim(-1.15, 1.85)
    axA.set_ylim(-2.6, 10.7)
    axA.set_title("A   Non-Teleological Inheritance in Human Clock Promoters", loc="left")
    axA.grid(True, axis="x", alpha=0.22, lw=0.6)
    axA.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.95, fontsize=7.1)

    # Panel B: AlphaGenome CRY1 32-bp ISM Attribution Profile in Human Liver
    axB = fig.add_subplot(gs[0, 1])
    rel_x = ism_dnase["pos_1based"].values - 107093585
    acgt_mat = ism_dnase[["A", "C", "G", "T"]].values
    l2_imp = np.linalg.norm(acgt_mat, axis=1)
    # CRY1 is on the minus strand. Complement the GRCh38 plus-strand bases at each plotted position.
    reference = pd.read_csv(REPO_ROOT / "data" / "module4_cry1_reference_bases.csv").set_index("pos_1based")
    reference_bases = reference.loc[ism_dnase["pos_1based"], "ref_plus"].tolist()
    complement = str.maketrans("ACGT", "TGCA")
    sense_bases = [base.translate(complement) for base in reference_bases]
    b_cols = [COL_CRIMSON if sb == "G" else (COL_NAVY if sb == "C" else COL_SLATE) for sb in sense_bases]
    axB.bar(rel_x, l2_imp, color=b_cols, width=0.78, edgecolor="#1C2833", lw=0.45)
    for x_i, imp_i, sb in zip(rel_x, l2_imp, sense_bases):
        if imp_i >= 0.015 or (-9 <= x_i <= 1):
            axB.text(x_i, imp_i + 0.008, sb, ha="center", va="bottom", fontsize=7.2, fontweight="bold",
                     color=COL_CRIMSON if sb == "G" else "#1C2833")
    axB.axvspan(-8.5, 1.5, color="#FDEDEC", alpha=0.50, zorder=0, label="SP1/KLF guanine hole-sink (sense 5'-GGGGGCGGGGC-3')")
    axB.set_xlabel("Genomic offset relative to CRY1 -36 bp C>A (sense G>T) locus (bp)")
    axB.set_ylabel("AlphaGenome liver DNase ISM attribution")
    axB.set_ylim(0, np.max(l2_imp) * 1.25)
    axB.set_title("B   AlphaGenome 32-bp ISM Across the CRY1 Motif", loc="left")
    axB.grid(True, axis="y", alpha=0.22, lw=0.6)
    axB.legend(loc="upper left", frameon=True, facecolor="white", framealpha=0.93, fontsize=7.6)

    # Panel C: Local +/- 250 bp AlphaGenome REF vs ALT Liver DNase Tracks at CRY1 TSS
    axC = fig.add_subplot(gs[1, 0])
    tss_cry1 = 107093549
    pos_1b = trk_dnase["pos_0based"].values + 1
    mask_t = (pos_1b >= tss_cry1 - 250) & (pos_1b <= tss_cry1 + 250)
    sub_t = trk_dnase[mask_t].copy()
    rel_tss = -(pos_1b[mask_t] - tss_cry1)  # negative = upstream of minus-strand CRY1 TSS
    order_t = np.argsort(rel_tss)
    rel_tss = rel_tss[order_t]
    ref_y = sub_t["ref_dnase_liver"].values[order_t]
    alt_y = sub_t["alt_dnase_liver"].values[order_t]
    axC.plot(rel_tss, ref_y, color=COL_NAVY, lw=1.9, label="Reference allele (intact G5 SP1 sink)")
    axC.plot(rel_tss, alt_y, color=COL_CRIMSON, ls="--", lw=1.8, label="8-oxoG > T transversion (chr12:107093585:C>A)")
    axC.fill_between(rel_tss, alt_y, ref_y, where=(ref_y > alt_y), color=COL_CRIMSON, alpha=0.20, label="Predicted DNase accessibility loss")
    axC.axvline(0, color=COL_TEAL, ls=":", lw=1.2, label="CRY1 MANE TSS")
    axC.axvline(-36, color=COL_AMBER, ls="-.", lw=1.1, label="Mutated G5 hole-sink (-36 bp)")
    axC.set_xlabel("Position relative to CRY1 transcription start site (bp, 5' -> 3')")
    axC.set_ylabel("Predicted human liver DNase signal")
    axC.set_title("C   Local Chromatin Accessibility Collapse at CRY1 TSS", loc="left")
    axC.grid(True, alpha=0.22, lw=0.6)
    axC.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.93, fontsize=7.3)

    # Panel D: Multi-modality AlphaGenome Disruption (Hole-Sink Motifs vs Flanking Isolated Guanines)
    axD = fig.add_subplot(gs[1, 1])
    mods = ["DNASE", "ATAC", "CHIP_TF", "CHIP_HISTONE", "CAGE", "RNA_SEQ"]
    mod_names = {"DNASE": "DNase", "ATAC": "ATAC", "CHIP_TF": "ChIP-TF",
                 "CHIP_HISTONE": "ChIP-Hist.", "CAGE": "CAGE", "RNA_SEQ": "RNA-seq"}
    mod_labels = [
        f"{mod_names[m]}\n({m4[m]['fold_raw_ratio']:.1f}x{'*' if m4[m]['mannwhitney_p_raw'] < 0.05 else ''})"
        for m in mods
    ]
    sink_sub = ag_df[ag_df["category"] == "Hole_Sink_Motif"]
    iso_sub = ag_df[ag_df["category"] == "Isolated_Flanking_G"]

    xm = np.arange(len(mods))
    sink_means = [np.mean(np.abs(sink_sub[f"{m}_liver_raw"].values)) for m in mods]
    iso_means = [np.mean(np.abs(iso_sub[f"{m}_liver_raw"].values)) for m in mods]
    sink_se = [stats.sem(np.abs(sink_sub[f"{m}_liver_raw"].values)) for m in mods]
    iso_se = [stats.sem(np.abs(iso_sub[f"{m}_liver_raw"].values)) for m in mods]

    axD.bar(xm - 0.18, sink_means, yerr=sink_se, capsize=2.5, width=0.35, color=COL_CRIMSON, edgecolor="#1C2833", lw=0.7,
            label=f"SP1 / E-box hole-sink guanines (n = {len(sink_sub)})")
    axD.bar(xm + 0.18, iso_means, yerr=iso_se, capsize=2.5, width=0.35, color=COL_SLATE, edgecolor="#1C2833", lw=0.7,
            label=f"Flanking isolated guanines (n = {len(iso_sub)})")
    rng = np.random.default_rng(42)
    for idx_m, m in enumerate(mods):
        ys = np.abs(sink_sub[f"{m}_liver_raw"].values)
        yi = np.abs(iso_sub[f"{m}_liver_raw"].values)
        axD.scatter(np.full_like(ys, xm[idx_m] - 0.18) + rng.uniform(-0.07, 0.07, len(ys)), ys,
                    color="#641E16", s=18, alpha=0.75, zorder=5)
        axD.scatter(np.full_like(yi, xm[idx_m] + 0.18) + rng.uniform(-0.07, 0.07, len(yi)), yi,
                    color="#1C2833", s=18, alpha=0.75, zorder=5)

    axD.set_xticks(xm)
    axD.set_xticklabels(mod_labels, fontsize=8.1)
    axD.set_ylabel("AlphaGenome liver |raw variant score|")
    axD.set_ylim(0, max(np.max(np.abs(ag_df[[f"{m}_liver_raw" for m in mods]].values)) * 1.12, 0.85))
    axD.set_title("D   Cis-Regulatory Vulnerability Across Clock Promoters", loc="left")
    axD.grid(True, axis="y", alpha=0.22, lw=0.6)
    axD.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.93, fontsize=7.8)

    save_pair(fig, "fig4_human_promoters_and_alphagenome")


if __name__ == "__main__":
    make_figure1()
    make_figure2()
    make_figure3()
    make_figure4()
    print("All 4 multi-panel figures generated successfully.")
