#!/usr/bin/env python3
"""Generate Figure 5 (PNG 300 dpi + SVG) with zero 34 GHz / CMB artifacts:
  Panel A: Earth-Sun 2D Quadrature Astrodynamic Time-Encoding (Chiral Solar
           Neutrino Keplerian Flux + 7Be Radial Doppler Quadrature vs.
           Seasonal/27-day Solar-Modulated Spin-Polarized Muon Clock).
  Panel B: Rewritten Quantum-Sequence & 150x150 Tight-Binding HOMO Hole-Sink
           Audit across 2,000 E. coli REL606 Windows (300 kbp).
  Panel C: Explicit 240,000-Trajectory 3D Monte Carlo Helical B-DNA Chiral
           Deflection Convergence (D-DNA m=+2 vs. L-DNA m=-2 & Zenith Sweep).
  Panel D: Chiral Particle Observability & Prebiotic Logistic Homochiralization
           (Direct Solar nu_e vs. Weak-Decay Spin-Polarized mu+/e+ & Young Sun).
"""

import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

DATA_DIR = "/home/sandbox/workspace/observer_bounded_dna_paper/data"
ASSET_DIR = "/home/sandbox/workspace/observer_bounded_dna_paper/assets"
os.makedirs(ASSET_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9.5,
    "axes.titlesize": 10.4,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.3,
    "xtick.labelsize": 8.2,
    "ytick.labelsize": 8.4,
    "legend.fontsize": 7.6,
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


def make_figure5():
    astro_df = pd.read_csv(f"{DATA_DIR}/module5_astrodynamics_doppler_and_flux.csv")
    mc_df = pd.read_csv(f"{DATA_DIR}/module5_monte_carlo_helical_convergence.csv")
    ee_df = pd.read_csv(f"{DATA_DIR}/module5_homochirality_evolution.csv")
    with open(f"{DATA_DIR}/module5_summary.json", "r") as fh:
        m5 = json.load(fh)

    fig = plt.figure(figsize=(13.0, 9.5), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, hspace=0.17, wspace=0.15)

    # -------------------------------------------------------------------------
    # Panel A: Earth-Sun 2D Quadrature Clock (Solar nu_e Flux/Doppler + Muons)
    # -------------------------------------------------------------------------
    axA = fig.add_subplot(gs[0, 0])
    days = astro_df["day_of_year"].values
    nu_flux_pct = astro_df["solar_nu_flux_mod_pct"].values
    muon_mod_pct = astro_df["combined_muon_clock_pct"].values
    be7_doppler_ev = astro_df["be7_nu_orbital_doppler_eV"].values

    nu_amp_pct = 0.5 * float(np.max(nu_flux_pct) - np.min(nu_flux_pct))
    be7_amp_ev = 0.5 * float(np.max(be7_doppler_ev) - np.min(be7_doppler_ev))
    l1, = axA.plot(days, nu_flux_pct, color=COL_NAVY, lw=2.3,
                   label=f"Solar $\\nu_e$ Keplerian flux $(a/r_\\oplus)^2$ ($\\pm {nu_amp_pct:.2f}\\%$, Borexino)")
    l2, = axA.plot(days, muon_mod_pct, color=COL_TEAL, lw=1.45, alpha=0.85,
                   label="Solar-mod. spin-polarized $\\mu^\\pm$ (seasonal + 27-d + diurnal)")
    axA.axhline(0.0, color=COL_SLATE, lw=0.6, alpha=0.5)
    axA.axvline(3.0, color=COL_SLATE, ls=":", lw=1.0, alpha=0.7)
    axA.axvline(185.5, color=COL_SLATE, ls=":", lw=1.0, alpha=0.7)
    axA.text(7.0, 4.45, "Perihelion\n(Jan 3)", fontsize=7.5, color=COL_SLATE, va="top")
    axA.text(189.0, 4.45, "Aphelion\n(Jul 4)", fontsize=7.5, color=COL_SLATE, va="top")

    axA2 = axA.twinx()
    axA2.spines["top"].set_visible(False)
    axA2.spines["right"].set_visible(True)
    axA2.spines["right"].set_color(COL_AMBER)
    l3, = axA2.plot(days, be7_doppler_ev, color=COL_AMBER, ls="--", lw=2.0,
                    label=f"$^7$Be solar $\\nu_e$ radial Doppler $\\Delta E_\\nu$ ($\\pm {be7_amp_ev:.2f}$ eV, $90^\\circ$ quad.)")
    axA2.set_ylabel("$^7$Be ($861.8$ keV) radial Doppler $\\Delta E_\\nu$ (eV)", color=COL_AMBER)
    axA2.tick_params(axis="y", colors=COL_AMBER)
    axA2.set_ylim(-2.15, 2.15)

    axA.set_xlim(0, 365.25)
    axA.set_ylim(-4.8, 4.85)
    axA.set_xlabel("Day of year (8,766-h J2000 Earth-Sun ephemeris)")
    axA.set_ylabel("Astrodynamic particle flux modulation (%)")
    axA.set_title("A   Earth-Sun Solar $\\nu_e$ Quadrature & Polarized $\\mu^\\pm$ Clocks", loc="left")
    axA.grid(True, alpha=0.22, lw=0.6)
    axA.legend(handles=[l1, l3, l2], loc="lower right", frameon=True,
               facecolor="white", framealpha=0.93, fontsize=7.3)

    # -------------------------------------------------------------------------
    # Panel B: Upgraded REL606 2,000-Window Code Audit (1,000 CDS + 1,000 IGR)
    # -------------------------------------------------------------------------
    axB = fig.add_subplot(gs[0, 1])
    gh_audit = m5["github_audit_and_upgrade"]
    w_info = gh_audit["real_rel606_windows"]
    cats = [
        "Orig. Unnorm.\nEntropy $S$\n($\\mathrm{Tr}\\,\\rho \\neq 1$)",
        "Repaired Bloch\nEntropy $S$\n($\\mathrm{Tr}\\,\\rho = 1$)",
        "Repaired\nDinuc. MI\n(bits)",
        "$150\\times 150$ $\\pi$-Stack\n$\\mathrm{G}_{\\geq 3}$ Hole\nCapture",
    ]
    c_means = [
        w_info["orig_unnorm_S_coding_mean"],
        w_info["repaired_bloch_S_coding_mean"],
        w_info["coding_dinuc_MI_mean_bits"],
        w_info["coding_G3plus_hole_frac_mean"],
    ]
    nc_means = [
        w_info["orig_unnorm_S_intergenic_mean"],
        w_info["repaired_bloch_S_intergenic_mean"],
        w_info["intergenic_dinuc_MI_mean_bits"],
        w_info["intergenic_G3plus_hole_frac_mean"],
    ]
    n_c = w_info["n_coding_windows"]
    n_nc = w_info["n_intergenic_windows"]
    c_sem = [
        w_info["orig_unnorm_S_coding_sd"] / np.sqrt(n_c),
        w_info["repaired_bloch_S_coding_sd"] / np.sqrt(n_c),
        w_info["coding_dinuc_MI_sd_bits"] / np.sqrt(n_c),
        w_info["coding_G3plus_hole_frac_sd"] / np.sqrt(n_c),
    ]
    nc_sem = [
        w_info["orig_unnorm_S_intergenic_sd"] / np.sqrt(n_nc),
        w_info["repaired_bloch_S_intergenic_sd"] / np.sqrt(n_nc),
        w_info["intergenic_dinuc_MI_sd_bits"] / np.sqrt(n_nc),
        w_info["intergenic_G3plus_hole_frac_sd"] / np.sqrt(n_nc),
    ]
    p_vals = [
        w_info["mwu_p_orig_unnorm_S"],
        w_info["mwu_p_repaired_bloch_S"],
        w_info["mwu_p_dinuc_MI"],
        w_info["mwu_p_G3plus_hole_frac"],
    ]

    x = np.arange(len(cats))
    w_bar = 0.34
    axB.bar(x - w_bar / 2, c_means, width=w_bar, yerr=c_sem, capsize=3.0,
            color=COL_NAVY, edgecolor="#1B2631", lw=0.8,
            label=f"Coding CDS ($n={n_c}$ windows)")
    axB.bar(x + w_bar / 2, nc_means, width=w_bar, yerr=nc_sem, capsize=3.0,
            color=COL_AMBER, edgecolor="#1B2631", lw=0.8,
            label=f"Intergenic ($n={n_nc}$ windows)")

    for idx in range(len(cats)):
        top_y = max(c_means[idx] + c_sem[idx], nc_means[idx] + nc_sem[idx])
        p_str = f"$p = {p_vals[idx]:.1e}$" if p_vals[idx] < 0.01 else f"$p = {p_vals[idx]:.2f}$"
        prefix = "GC artifact\n" if idx == 0 else ("AT/GC skew\n" if p_vals[idx] < 0.01 else "")
        col_txt = COL_CRIMSON if idx == 0 else "#1B2631"
        axB.text(idx, top_y + 0.025, f"{prefix}{p_str}", ha="center", va="bottom",
                 fontsize=7.3, color=col_txt, fontweight=("bold" if idx == 0 else "normal"))

    qmi_abs = abs(float(gh_audit["product_state_qmi_bits"]))
    neg_val = abs(float(gh_audit["product_state_negativity"]))
    axB.text(2.50, 0.78,
             "Product-state check ($\\rho_1 \\otimes \\rho_2$):\n"
             f"Negativity $\\mathcal{{N}} = {neg_val:.1f}$ (exact),\n"
             f"$|I_Q| = {qmi_abs:.1e}$ bits",
             ha="center", va="center", fontsize=7.4, color="#1B2631",
             bbox=dict(boxstyle="round,pad=0.30", facecolor="#F8F9F9", edgecolor=COL_SLATE, lw=0.7))

    axB.set_xticks(x)
    axB.set_xticklabels(cats)
    axB.set_ylim(0.0, 1.27)
    axB.set_ylabel("Entropy / MI (bits) or hole-sink fraction")
    axB.set_title(f"B   Upgraded Quantum-DNA Audit ({n_c + n_nc} REL606 Windows)", loc="left")
    axB.grid(True, axis="y", alpha=0.22, lw=0.6)
    axB.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.93, fontsize=7.5)

    # -------------------------------------------------------------------------
    # Panel C: 240,000-Trajectory 3D Monte Carlo Helical Chiral Deflection
    # -------------------------------------------------------------------------
    axC = fig.add_subplot(gs[1, 0])
    max_abs_m = float(np.max(np.abs(mc_df["mean_dr_m"].values)))
    exp_scale = int(np.floor(np.log10(max(max_abs_m, 1e-40))))
    scale_dr = 10.0 ** (-exp_scale)
    styles = {
        "D_DNA_m_plus2_isotropic_pos": (COL_NAVY, "o", "-", "D-DNA ($m=+2$, isotropic $\\mu^+/e^+$, $\\mathcal{L}=-0.5$)"),
        "L_DNA_m_minus2_isotropic_pos": (COL_CRIMSON, "s", "-", "L-DNA ($m=-2$, isotropic $\\mu^+/e^+$, $\\mathcal{L}=-0.5$)"),
        "D_DNA_m_plus2_helicity_rev_ctrl": (COL_PURPLE, "^", "--", "D-DNA ($m=+2$, helicity-flipped $\\mathcal{L}=+0.5$ ctrl)"),
        "D_DNA_m_plus2_zenith_oriented": (COL_TEAL, "D", "-.", "D-DNA ($m=+2$, flat terrace $\\cos^2\\theta_z$ shower)"),
    }

    y_peak = 0.0
    for ens_name, (col, mk, ls, lbl) in styles.items():
        sub = mc_df[mc_df["ensemble"] == ens_name].sort_values("n_trajectories")
        n_tr = sub["n_trajectories"].values
        y_val = sub["mean_dr_m"].values * scale_dr
        y_ci = 1.96 * sub["sem_dr_m"].values * scale_dr
        y_peak = max(y_peak, float(np.max(np.abs(y_val) + y_ci)))
        axC.semilogx(n_tr, y_val, color=col, marker=mk, ls=ls, lw=1.9, ms=5.2, label=lbl)
        axC.fill_between(n_tr, y_val - y_ci, y_val + y_ci, color=col, alpha=0.15)

    cc_info = m5["cosmic_ray_chiral_mutagenesis"]
    g4_cap_pct = 100.0 * float(cc_info["g4_hole_sink_capture_fraction"])
    g4_fold_v = float(cc_info["g4_hole_sink_enrichment_fold"])
    n_tot_traj = int(m5["monte_carlo_helical_deflection"]["total_3d_trajectories_integrated"])

    ylim_c = max(1.82 * y_peak, 1.0)
    axC.axhline(0.0, color="#1B2631", lw=0.8, alpha=0.65)
    axC.axhspan(-ylim_c, 0.0, color="#D4EFDF", alpha=0.18)
    axC.axhspan(0.0, ylim_c, color="#FADBD8", alpha=0.18)
    axC.text(540, -0.94 * ylim_c,
             "Inward deflection $\\langle\\Delta r\\rangle < 0$: central nucleobase core\n"
             f"$\\rightarrow$ $\\mathrm{{G}}_4$ HOMO hole sink (${g4_cap_pct:.1f}\\%$ capture, ${g4_fold_v:.1f}\\times$ 8-oxoG)",
             fontsize=7.3, color=COL_TEAL, fontweight="bold", va="bottom")
    axC.text(540, 0.94 * ylim_c,
             "Outward deflection $\\langle\\Delta r\\rangle > 0$: outer sugar-phosphate backbone\n"
             "$\\rightarrow$ lethal single/double-strand scission",
             fontsize=7.3, color=COL_CRIMSON, fontweight="bold", va="top")

    axC.set_xlim(400, 8.0e4)
    axC.set_ylim(-ylim_c, ylim_c)
    axC.set_xlabel(f"3D Monte Carlo trajectories per ensemble $N_{{\\mathrm{{traj}}}}$ (${n_tot_traj:,}$ total)")
    axC.set_ylabel(f"Mean chiral radial deflection $\\langle\\Delta r\\rangle$ ($10^{{{exp_scale}}}$ m)")
    axC.set_title("C   3D Monte Carlo Helical B-DNA Chiral Deflection Convergence", loc="left")
    axC.grid(True, which="major", alpha=0.22, lw=0.6)
    axC.legend(loc="upper right", bbox_to_anchor=(0.99, 0.84), frameon=True, facecolor="white", framealpha=0.95, fontsize=6.9)

    # -------------------------------------------------------------------------
    # Panel D: Solar Neutrino vs. Weak-Decay Polarized Muon/SEP Homochirality
    # -------------------------------------------------------------------------
    axD = fig.add_subplot(gs[1, 1])
    palette_dm = {1.0e-7: COL_SLATE, 1.0e-6: COL_NAVY, 1.0e-5: COL_TEAL}
    for dm_val, col in palette_dm.items():
        sub = ee_df[np.isclose(ee_df["delta_M"], dm_val)]
        exp_int = int(round(np.log10(dm_val)))
        axD.semilogx(sub["generations_tau"], sub["ee_with_chiral_conflict"], color=col, lw=2.2,
                     label=f"$\\delta M = 10^{{{exp_int}}}$ (Young-Sun SEP / conflict)")
        axD.semilogx(sub["generations_tau"], sub["ee_no_conflict"], color=col, ls="--", lw=1.35, alpha=0.8,
                     label=f"$\\delta M = 10^{{{exp_int}}}$ (pure $\\tanh(\\delta M \\tau / 2)$)")

    nu_rate = m5["solar_neutrino_astrodynamics"]["diploid_DNA_nu_e_scatter_rate_per_yr"]
    em_ratio = m5["solar_modulated_muon_astrodynamics"]["electromagnetic_to_weak_cross_section_ratio"]
    axD.axhline(0.99, color=COL_CRIMSON, ls=":", lw=1.2)
    axD.text(1.5e2, 0.94, "Homochiral fixation ($\\mathrm{e.e.} = 0.99$)", fontsize=7.6, color=COL_CRIMSON)
    axD.text(
        1.45e2, 0.50,
        f"Direct solar $\\nu_e$ ($G_F^2$ gated):\n"
        f"  $\\Gamma_{{\\nu,\\mathrm{{DNA}}}} = {nu_rate:.1e}\\ \\mathrm{{cell^{{-1}}\\,yr^{{-1}}}}$\n"
        f"  (needs $\\sim 10^{{17}}$ cells / $100$-t Borexino)\n"
        f"Weak decay $\\pi^+ \\to \\mu^+ + \\nu_\\mu$:\n"
        f"  transfers $\\nu$ chirality to EM $\\mu^+/e^+$\n"
        f"  ($\\sigma_{{\\mathrm{{ion}}}}/\\sigma_{{\\nu e}} \\approx {em_ratio:.1e}\\times$ stronger!)",
        fontsize=7.35, color="#1B2631",
        bbox=dict(boxstyle="round,pad=0.34", facecolor="#F8F9F9", edgecolor=COL_SLATE, lw=0.7),
    )

    axD.set_xlim(1e2, 1e8)
    axD.set_ylim(-0.02, 1.05)
    axD.set_xlabel("Evolutionary radiation-weighted turnovers $\\tau$")
    axD.set_ylabel("Biopolymer enantiomeric excess $\\mathrm{e.e.}(\\tau)$")
    axD.set_title("D   Solar $\\nu_e$ Bottleneck vs. Polarized $\\mu^+/e^+$ Homochiralization", loc="left")
    axD.grid(True, which="major", alpha=0.22, lw=0.6)
    axD.legend(loc="lower right", frameon=True, facecolor="white", framealpha=0.93, fontsize=7.1)

    png_path = f"{ASSET_DIR}/fig5_astrodynamics_and_cosmic_chirality.png"
    svg_path = f"{ASSET_DIR}/fig5_astrodynamics_and_cosmic_chirality.svg"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {png_path} ({os.path.getsize(png_path)} bytes) and {svg_path} ({os.path.getsize(svg_path)} bytes)")


if __name__ == "__main__":
    make_figure5()
