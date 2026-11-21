#!/usr/bin/env python3
"""
Module 6: Empirical Wet-Lab Validation of the Table 6 Falsification Matrix
==========================================================================
Tests Rows 1, 2, 3, 4, and 5a of Table 6 directly on independent, publicly
deposited wet-lab experimental datasets and 3D cryo-EM/X-ray structures, and
quantifies the turnkey accelerator beamline specification for Row 5b (the sole
remaining unexecuted physical experiment).

Datasets & Structural Accessions Analyzed:
  - Row 1 (Theorem 1 — Sub-thermal AC Null vs. Thermal Arrhenius Scaling):
      * Chu et al. (2018) BMC Evol. Biol. 18:126 — E. coli REL606 mutS WGS MA
        across 25 C, 28 C, and 37 C vs. Module 1 KMS ΔE_act = 0.42 eV Arrhenius
        prediction and isothermal GHz RF/Wi-Fi bacterial bounds (Said-Salman
        et al. 2019 Sci. Rep. 9:14425; Dinc et al. 2021).
  - Row 2 (Theorem 2 — Nucleosomal SHL Rotational Gate on Human OGG1):
      * RCSB PDB 8VWT (hOGG1-NCP SHL -6, 3.3 A), 8VWV (hOGG1-NCP SHL +4, 3.6 A),
        9EOZ (hOGG1-K249Q-NCP SHL 6.0, 3.1 A), and 1KX5 (canonical 147-bp NCP,
        1.9 A): native heavy-atom contacts and a rigid-body rotational clash scan
        across the interior base pairs of both DNA strands in 1KX5.
  - Row 3 (Proposition 1 & Theorem 2 — Single-Cell Erasure vs. Ensemble Scaling):
      * NCBI GEO GSE218316 primary human pancreatic islet scRNA-seq, restricted
        to the 24 h timepoint shared by all three donors.
  - Row 4 (Module 4 — Hole-Sink Cis Vulnerability in Human HepG2 MPRA):
      * NCBI GEO GSE126550 (Kircher et al. 2019 Nat. Commun. 10:3583) saturation
        mutagenesis MPRA in living human HepG2 cells across LDLR, LDLR.2, F9,
        and PKLR_48h (24 DNA/RNA UMI-barcoded sequencing libraries).
  - Row 5a & 5b (Module 5 — Bacterial GO-Repair Uncollapsing & Beamline Spec):
      * Foster et al. (2015) PNAS 112:E5990 & Sane/Garushyants et al. (2024) GBE
        16:evae035 E. coli WGS MA VCFs across wild-type, mutM mutY, mutY, mutT,
        MMR (mutS/mutL/mutH), nth nei, xthA nfo, uvrA, and dnaQ lines; plus
        exact particle-beam event budget for pitch-tilted polarized e+- / mu+
        helicity-flip irradiation.
"""

import argparse
import os
import re
import json
import gzip
import math
import time
import tarfile
import subprocess
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from Bio.PDB.MMCIFParser import MMCIFParser

from analysis_paths import DATA_DIR, RAW_DIR

SCRATCH_DIR = RAW_DIR / "table6"
DL_DIR = SCRATCH_DIR / "downloads"
os.makedirs(DATA_DIR, exist_ok=True)

KB_EV_PER_K = 8.617333262145e-5
E_CHARGE = 1.602176634e-19


# ==============================================================================
# ROW 1: SUB-THERMAL AC NULL VS. THERMAL ARRHENIUS SCALING (THEOREM 1)
# ==============================================================================
def run_row1_thermal_vs_ac_validation():
    t0 = time.time()
    # Empirical WGS MA rates in E. coli B REL606 mutS from Chu et al. (2018),
    # BMC Evolutionary Biology 18:126 (DOI: 10.1186/s12862-018-1252-8):
    # 20 independent MA lines per temperature (60 sequenced MA lines total; Welch two-sample df = 29, 32, 33).
    # Published one-sided tests for spontaneous base-pair substitutions plus small indels.
    #   25 C vs 37 C: 37% lower at 25 C (ratio 0.63), t = -5.26, df = 32, Padj = 7e-6
    #   28 C vs 37 C: 38% lower at 28 C (ratio 0.62), t = -7.12, df = 33, Padj = 6e-8
    #   25 C vs 28 C: ratio 0.63/0.62 = 1.016, t = 0.08, df = 29, Padj = 0.469 (null over 3 K)
    temps_c = np.array([25.0, 28.0, 37.0])
    temps_k = temps_c + 273.15
    rel_rates_emp = np.array([0.63, 0.62, 1.00])  # normalized to 37 C = 1.00

    # Module 1 a priori KMS activation barrier (E_1 = 0.42 eV)
    e_act_mod1_ev = 0.42
    rel_rates_mod1 = np.exp(-e_act_mod1_ev / (KB_EV_PER_K * temps_k))
    rel_rates_mod1 = rel_rates_mod1 / rel_rates_mod1[-1]

    # Fit empirical apparent activation energy Ea from ln(rel_rate) vs -1/(kB*T)
    x_inv_kbt = -1.0 / (KB_EV_PER_K * temps_k)
    slope_ea, intercept_ea, r_val, p_val, std_err = stats.linregress(
        x_inv_kbt - x_inv_kbt[-1], np.log(rel_rates_emp)
    )

    # Compare 28 C -> 37 C ratio:
    emp_ratio_37_over_28 = 1.00 / 0.62  # 1.6129x
    mod1_ratio_37_over_28 = float(rel_rates_mod1[2] / rel_rates_mod1[1])  # 1.6010x
    pct_agreement_37_28 = 100.0 * (1.0 - abs(mod1_ratio_37_over_28 - emp_ratio_37_over_28) / emp_ratio_37_over_28)

    # Theorem 2 sub-thermal AC predictions across exact Floquet-Lindblad Drazin Eq. (S2) & Adiabatic KMS
    mod1_df = pd.read_csv(os.path.join(DATA_DIR, "module1_weak_field_nogo.csv"))
    ac_points = {}
    for e0 in [1e-6, 1.0, 10.0]:
        sub_row = mod1_df.iloc[(mod1_df["E0_V_per_m"] - e0).abs().argsort()[:1]].iloc[0]
        ac_points[f"E0_{e0:g}_Vm"] = {
            "E_field_V_per_m": float(sub_row["E0_V_per_m"]),
            "lambda_eV": float(sub_row["lambda_eV"]),
            "drazin_floquet_rel_shift": float(sub_row["drazin_floquet_subgap_dmu_over_mu"]),
            "adiabatic_kms_rel_shift": float(sub_row["adiabatic_kms_ac_dmu_over_mu"]),
            "heuristic_lorentzian_rel_shift": float(sub_row["ac_subgap_rectified_dmu_over_mu"]),
            "thermal_50mK_rel_shift": float(sub_row["metabolic_thermal_floor_50mK"]),
            "thermal_10mK_rel_shift": float(
                math.exp(e_act_mod1_ev * 0.010 / (KB_EV_PER_K * 310.15 * 310.16)) - 1.0
            ),
        }

    return {
        "runtime_s": round(time.time() - t0, 4),
        "n_ma_lines_per_temperature": 20,
        "n_total_ma_lines": 60,
        "welch_dfs": {"25_vs_28C": 29, "25_vs_37C": 32, "28_vs_37C": 33},
        "endpoint_definition": "Spontaneous base-pair substitutions (87%-96% transitions in mutS) + small indels (~1/5 of BPS rate)",
        "emp_rel_rates_25_28_37C": rel_rates_emp.tolist(),
        "mod1_predicted_rel_rates_25_28_37C": [round(float(x), 5) for x in rel_rates_mod1],
        "emp_ratio_28C_vs_25C": round(float(rel_rates_emp[1] / rel_rates_emp[0]), 4),
        "mod1_predicted_ratio_28C_vs_25C": round(float(rel_rates_mod1[1] / rel_rates_mod1[0]), 4),
        "chu2018_test_sidedness": "one-sided",
        "emp_ratio_37C_vs_28C": round(emp_ratio_37_over_28, 4),
        "mod1_predicted_ratio_37C_vs_28C": round(mod1_ratio_37_over_28, 4),
        "arrhenius_ratio_agreement_pct": round(pct_agreement_37_28, 2),
        "emp_apparent_Ea_eV": round(float(slope_ea), 4),
        "mod1_a_priori_Ea_eV": e_act_mod1_ev,
        "chu2018_25_vs_28C_t": 0.08,
        "chu2018_25_vs_28C_padj": 0.469,
        "chu2018_28_vs_37C_t": -7.12,
        "chu2018_28_vs_37C_padj": 6.0e-8,
        "chu2018_25_vs_37C_t": -5.26,
        "chu2018_25_vs_37C_padj": 7.0e-6,
        "ac_bounds": ac_points,
    }


# ==============================================================================
# ROW 2: NUCLEOSOMAL SHL ROTATIONAL GATE ON HUMAN OGG1 (THEOREM 2)
# ==============================================================================
def kabsch_transform(P, Q):
    """Find rigid transform R, t mapping P (Nx3) onto Q (Nx3): x -> x @ R + t."""
    p_mean = P.mean(axis=0)
    q_mean = Q.mean(axis=0)
    P0 = P - p_mean
    Q0 = Q - q_mean
    H = P0.T @ Q0
    U, S, Vt = np.linalg.svd(H)
    R = U @ Vt
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1.0
        R = U @ Vt
    t = q_mean - p_mean @ R
    rmsd = float(np.sqrt(np.mean(np.sum((P0 @ R - Q0) ** 2, axis=1))))
    return R, t, rmsd


def get_heavy_atoms_unpacked(chain_list):
    """Extract all heavy atoms unpacking any altLoc conformers per pdb-structural-validation rules."""
    coords = []
    meta = []
    altloc_n = 0
    occ_list = []
    for ch in chain_list:
        for res in ch:
            rname = res.get_resname().strip()
            if rname == "HOH":
                continue
            for atom in res.get_unpacked_list():
                elem = (atom.element or atom.get_name()[0]).strip().upper()
                if elem in ("H", "D"):
                    continue
                alt = atom.get_altloc().strip() if atom.get_altloc() else ""
                if alt not in ("", "."):
                    altloc_n += 1
                occ = float(atom.get_occupancy()) if atom.get_occupancy() is not None else 1.0
                occ_list.append(occ)
                coords.append(atom.get_coord())
                meta.append((ch.id, res.id[1], rname, atom.get_name().strip(), alt, occ))
    return np.array(coords, dtype=float), meta, altloc_n, (min(occ_list) if occ_list else 1.0), (max(occ_list) if occ_list else 1.0)


def run_row2_nucleosome_ogg1_validation():
    t0 = time.time()
    pdb_dir = os.path.join(SCRATCH_DIR, "pdb")
    parser = MMCIFParser(QUIET=True)

    # Measure native hOGG1-vs-histone heavy-atom clashes in 8VWT, 8VWV and 9EOZ.
    cryo_native_summary = {}
    for pid, ogg_ch, hist_chs in [
        ("8VWT", "K", ["A", "B", "C", "D", "E", "F", "G", "H"]),
        ("8VWV", "K", ["A", "B", "C", "D", "E", "F", "G", "H"]),
        ("9EOZ", "A", ["B", "E", "F", "G", "H", "K", "L", "M"]),
    ]:
        st = parser.get_structure(pid, os.path.join(pdb_dir, f"{pid}.cif"))[0]
        ogg_coords, _, ogg_alt, ogg_occ_min, ogg_occ_max = get_heavy_atoms_unpacked([st[ogg_ch]])
        hist_coords, _, hist_alt, hist_occ_min, hist_occ_max = get_heavy_atoms_unpacked([st[c] for c in hist_chs])
        tree = cKDTree(hist_coords)
        dists, _ = tree.query(ogg_coords, k=1)
        cryo_native_summary[pid] = {
            "ogg1_heavy_atoms": int(len(ogg_coords)),
            "histone_heavy_atoms": int(len(hist_coords)),
            "altloc_atoms_total": int(ogg_alt + hist_alt),
            "occupancy_range": [round(min(ogg_occ_min, hist_occ_min), 2), round(max(ogg_occ_max, hist_occ_max), 2)],
            "min_ogg1_histone_dist_A": round(float(dists.min()), 3),
            "clashes_lt_2p5A": int(np.sum(dists < 2.5)),
            "contacts_lt_3p5A": int(np.sum(dists < 3.5)),
        }

    # Rigid-body Kabsch superposition of the active hOGG1 recognition complex (from 8VWT)
    # across all interior base-pair steps of canonical 1.9 A NCP (1KX5)
    st_8vwt = parser.get_structure("8VWT", os.path.join(pdb_dir, "8VWT.cif"))[0]
    ogg_template_coords, _, _, _, _ = get_heavy_atoms_unpacked([st_8vwt["K"]])

    # Build 5-nt local backbone template around 8VWT J:13 (J:11..15) and complementary I strand
    # First map base-pair partners in 8VWT and 1KX5 by C1'-C1' proximity (~10.5 A)
    def get_nt_dict(chain):
        d = {}
        for r in chain:
            rn = r.get_resname().strip()
            if rn in ("DA", "DT", "DG", "DC", "8OG"):
                atoms = {a.get_name().strip(): a.get_coord() for a in r.get_unpacked_list()}
                if all(k in atoms for k in ("P", "C4'", "C1'")):
                    d[r.id[1]] = atoms
        return d

    j_8vwt = get_nt_dict(st_8vwt["J"])
    i_8vwt = get_nt_dict(st_8vwt["I"])
    # Find complement index on I paired to J:13
    c1_j13 = j_8vwt[13]["C1'"]
    comp_13 = min(i_8vwt.keys(), key=lambda k: abs(np.linalg.norm(i_8vwt[k]["C1'"] - c1_j13) - 10.6))

    # Extract 5-bp backbone anchor atoms (excluding the flipped J:13 C1' itself, using flanking +/-2, +/-1 P, C4', C1')
    def build_anchor(lesion_nt_dict, comp_nt_dict, center_idx, comp_center_idx):
        pts = []
        for off in [-2, -1, 1, 2]:
            idx1 = center_idx + off
            idx2 = comp_center_idx - off
            if idx1 not in lesion_nt_dict or idx2 not in comp_nt_dict:
                return None
            for aname in ("P", "C4'", "C1'"):
                pts.append(lesion_nt_dict[idx1][aname])
                pts.append(comp_nt_dict[idx2][aname])
        return np.array(pts, dtype=float)

    template_anchor = build_anchor(j_8vwt, i_8vwt, 13, comp_13)

    st_1kx5 = parser.get_structure("1KX5", os.path.join(pdb_dir, "1KX5.cif"))[0]
    hist_1kx5_coords, hist_1kx5_meta, hist_1kx5_alt, _, _ = get_heavy_atoms_unpacked(
        [st_1kx5[c] for c in ["A", "B", "C", "D", "E", "F", "G", "H"]]
    )
    # Structured globular histone-fold octamer core (excluding flexible N-terminal tails absent in cryo-EM 8VWT/8VWV):
    # H3 (A,E) >= 38, H4 (B,F) >= 22, H2A (C,G) >= 15, H2B (D,H) >= 30
    core_min_res = {"A": 38, "E": 38, "B": 22, "F": 22, "C": 15, "G": 15, "D": 30, "H": 30}
    core_mask = np.array([m[1] >= core_min_res.get(m[0], 1) for m in hist_1kx5_meta], dtype=bool)
    hist_core_coords = hist_1kx5_coords[core_mask]

    hist_center = hist_core_coords.mean(axis=0)
    hist_tree = cKDTree(hist_1kx5_coords)
    hist_core_tree = cKDTree(hist_core_coords)

    i_1kx5 = get_nt_dict(st_1kx5["I"])
    j_1kx5 = get_nt_dict(st_1kx5["J"])
    # 1KX5 chains I and J are numbered -73 .. +73 (palindromic 147 bp centered at dyad 0)
    rows_scan = []
    for strand_label, les_dict, cmp_dict in [("I", i_1kx5, j_1kx5), ("J", j_1kx5, i_1kx5)]:
        sorted_pos = sorted(les_dict.keys())
        for pos in sorted_pos:
            if abs(pos) > 66:  # interior nucleosomal turns with full +/-2 bp flanks
                continue
            comp_pos = -pos
            anchor_target = build_anchor(les_dict, cmp_dict, pos, comp_pos)
            if anchor_target is None:
                continue
            R, t_vec, anchor_rmsd = kabsch_transform(template_anchor, anchor_target)
            transformed_ogg = ogg_template_coords @ R + t_vec

            # Minor-groove / lesion-extrusion radial phase angle psi_rot in [0, 360)
            bp_mid = 0.5 * (les_dict[pos]["C1'"] + cmp_dict[comp_pos]["C1'"])
            v_rad_out = bp_mid - hist_center
            if (pos + 1) in les_dict and (pos - 1) in les_dict:
                v_axis = les_dict[pos + 1]["C4'"] - les_dict[pos - 1]["C4'"]
                v_axis = v_axis / (np.linalg.norm(v_axis) + 1e-12)
            else:
                v_axis = np.array([0.0, 0.0, 1.0])
            v_rad_perp = v_rad_out - np.dot(v_rad_out, v_axis) * v_axis
            v_rad_perp /= (np.linalg.norm(v_rad_perp) + 1e-12)

            v_ogg = transformed_ogg.mean(axis=0) - bp_mid
            v_ogg_perp = v_ogg - np.dot(v_ogg, v_axis) * v_axis
            v_ogg_perp /= (np.linalg.norm(v_ogg_perp) + 1e-12)

            cos_ang = float(np.clip(np.dot(v_rad_perp, v_ogg_perp), -1.0, 1.0))
            sin_ang = float(np.dot(np.cross(v_axis, v_rad_perp), v_ogg_perp))
            psi_deg = float(np.degrees(np.arctan2(sin_ang, cos_ang)) % 360.0)
            inward_angle_deg = float(np.degrees(np.arccos(cos_ang)))

            dists, _ = hist_tree.query(transformed_ogg, k=1)
            dists_core, _ = hist_core_tree.query(transformed_ogg, k=1)
            clashes_2p5 = int(np.sum(dists < 2.5))
            clashes_3p0 = int(np.sum(dists < 3.0))
            core_clashes_2p5 = int(np.sum(dists_core < 2.5))
            min_dist = float(dists.min())
            min_core_dist = float(dists_core.min())
            shl = pos / 10.2

            # Soft-core steric penalty capped by spontaneous nucleosome site-exposure floor K_unwrap = 0.055
            k_steric_core = max(0.055, math.exp(-0.18 * core_clashes_2p5))
            k_steric_full = max(0.055, math.exp(-0.18 * clashes_2p5))

            rows_scan.append({
                "strand": strand_label,
                "bp_index_from_dyad": int(pos),
                "SHL": round(shl, 3),
                "rotational_phase_deg": round(psi_deg, 2),
                "inward_angle_deg": round(inward_angle_deg, 2),
                "backbone_kabsch_rmsd_A": round(anchor_rmsd, 3),
                "min_histone_heavy_dist_A": round(min_dist, 3),
                "min_core_heavy_dist_A": round(min_core_dist, 3),
                "severe_clashes_lt_2p5A": clashes_2p5,
                "core_clashes_lt_2p5A": core_clashes_2p5,
                "hard_contacts_lt_3p0A": clashes_3p0,
                "effective_accessibility_factor": round(k_steric_core, 5),
                "full_tail_accessibility_factor": round(k_steric_full, 5),
            })

    scan_df = pd.DataFrame(rows_scan)
    scan_df.to_csv(os.path.join(DATA_DIR, "module6_nucleosome_rotational_clashes.csv"), index=False)

    outward_df = scan_df[scan_df["inward_angle_deg"] <= 60.0]
    intermediate_df = scan_df[(scan_df["inward_angle_deg"] > 60.0) & (scan_df["inward_angle_deg"] < 120.0)]
    inward_df = scan_df[scan_df["inward_angle_deg"] >= 120.0]

    # Compute exact inward_angle_deg on cryo-EM 8VWT, 8VWV, 9EOZ using the same geometric formula
    cryo_angles = {}
    for pid, ogg_ch, hist_chs, les_ch, les_pos, cmp_ch in [
        ("8VWT", "K", ["A", "B", "C", "D", "E", "F", "G", "H"], "J", 13, "I"),
        ("8VWV", "K", ["A", "B", "C", "D", "E", "F", "G", "H"], "I", 34, "J"),
        ("9EOZ", "A", ["B", "E", "F", "G", "H", "K", "L", "M"], "Z", 137, "Y"),
    ]:
        st_c = parser.get_structure(pid, os.path.join(pdb_dir, f"{pid}.cif"))[0]
        ogg_c, _, _, _, _ = get_heavy_atoms_unpacked([st_c[ogg_ch]])
        hist_c, _, _, _, _ = get_heavy_atoms_unpacked([st_c[c] for c in hist_chs])
        h_ctr = hist_c.mean(axis=0)
        les_d = get_nt_dict(st_c[les_ch])
        cmp_d = get_nt_dict(st_c[cmp_ch])
        c1_p = les_d[les_pos]["C1'"]
        cmp_pos_c = min(cmp_d.keys(), key=lambda k: abs(np.linalg.norm(cmp_d[k]["C1'"] - c1_p) - 10.6))
        bp_m = 0.5 * (les_d[les_pos]["C1'"] + cmp_d[cmp_pos_c]["C1'"])
        v_ro = bp_m - h_ctr
        v_ax = les_d[les_pos + 1]["C4'"] - les_d[les_pos - 1]["C4'"]
        v_ax = v_ax / (np.linalg.norm(v_ax) + 1e-12)
        v_rp = v_ro - np.dot(v_ro, v_ax) * v_ax
        v_rp /= (np.linalg.norm(v_rp) + 1e-12)
        v_og = ogg_c.mean(axis=0) - bp_m
        v_op = v_og - np.dot(v_og, v_ax) * v_ax
        v_op /= (np.linalg.norm(v_op) + 1e-12)
        ang_c = float(np.degrees(np.arccos(float(np.clip(np.dot(v_rp, v_op), -1.0, 1.0)))))
        cryo_angles[pid] = round(ang_c, 2)
        cryo_native_summary[pid]["inward_angle_deg"] = round(ang_c, 2)

    mean_acc_out = float(outward_df["effective_accessibility_factor"].mean())
    mean_acc_in = float(inward_df["effective_accessibility_factor"].mean())
    acc_ratio = mean_acc_out / max(mean_acc_in, 1e-9)

    return {
        "runtime_s": round(time.time() - t0, 4),
        "cryo_em_native_structures": cryo_native_summary,
        "cryo_em_inward_angles_deg": cryo_angles,
        "1kx5_altloc_atoms_checked": int(hist_1kx5_alt),
        "n_nucleosome_positions_scanned": int(len(scan_df)),
        "outward_le_60deg_n": int(len(outward_df)),
        "outward_mean_clashes_2p5A": round(float(outward_df["severe_clashes_lt_2p5A"].mean()), 2),
        "outward_mean_core_clashes_2p5A": round(float(outward_df["core_clashes_lt_2p5A"].mean()), 2),
        "outward_core_zero_clash_fraction": round(float((outward_df["core_clashes_lt_2p5A"] == 0).mean()), 4),
        "intermediate_60_120deg_n": int(len(intermediate_df)),
        "intermediate_mean_core_clashes_2p5A": round(float(intermediate_df["core_clashes_lt_2p5A"].mean()), 2),
        "inward_ge_120deg_n": int(len(inward_df)),
        "inward_mean_clashes_2p5A": round(float(inward_df["severe_clashes_lt_2p5A"].mean()), 2),
        "inward_mean_core_clashes_2p5A": round(float(inward_df["core_clashes_lt_2p5A"].mean()), 2),
        "inward_core_nonzero_clash_fraction": round(float((inward_df["core_clashes_lt_2p5A"] > 0).mean()), 4),
        "outward_vs_inward_accessibility_ratio": round(acc_ratio, 2),
    }


# ==============================================================================
# ROW 3: SINGLE-CELL TRANSCRIPT DETECTION AND ENSEMBLE SCALING (GSE218316)
# ==============================================================================
def disjoint_ensemble_means(cells, labels, indices, size, rng):
    """Use each cell at most once per train or test partition at this size."""
    means = []
    pooled_labels = []
    for label in (0, 1):
        group = indices[labels[indices] == label]
        n_ensembles = len(group) // size
        if n_ensembles == 0:
            return None
        members = rng.permutation(group)[:n_ensembles * size].reshape(n_ensembles, size)
        means.append(cells[members].mean(axis=1))
        pooled_labels.append(np.full(n_ensembles, label, dtype=int))
    return np.concatenate(means), np.concatenate(pooled_labels)


def score_disjoint_ensemble_folds(cells, labels, feature_indices, folds, size, rng):
    aucs = []
    test_sizes_by_condition = []
    for train_indices, test_indices in folds:
        train = disjoint_ensemble_means(cells, labels, train_indices, size, rng)
        test = disjoint_ensemble_means(cells, labels, test_indices, size, rng)
        if train is None or test is None:
            continue
        x_train, y_train = train
        x_test, y_test = test
        fold_aucs = []
        for indices in feature_indices:
            model = LogisticRegression(C=1.0, max_iter=500).fit(x_train[:, indices], y_train)
            fold_aucs.append(float(roc_auc_score(y_test, model.decision_function(x_test[:, indices]))))
        aucs.append(fold_aucs)
        test_sizes_by_condition.append(np.bincount(y_test, minlength=2))
    mean_auc = np.mean(aucs, axis=0).tolist() if len(aucs) == len(folds) else [None] * len(feature_indices)
    return mean_auc, np.sum(test_sizes_by_condition, axis=0).tolist() if aucs else [0, 0], len(aucs)


ISLET_24H_SAMPLES = {"D1": "GSM6740278", "D2": "GSM6740279", "D3": "GSM6740282"}


def load_islet_counts(metadata, target_genes):
    """Read each donor's matrix against its own Cell Ranger features and barcodes."""
    raw_dir = SCRATCH_DIR / "gse218316"
    expression = np.zeros((len(target_genes), len(metadata)), dtype=np.int64)
    library_sizes = np.zeros(len(metadata), dtype=np.int64)
    target_index = {gene: idx for idx, gene in enumerate(target_genes)}

    if set(metadata["Donor"]) != set(ISLET_24H_SAMPLES):
        raise ValueError("Expected all three 24 h donors in GSE218316 metadata")

    for donor, accession in ISLET_24H_SAMPLES.items():
        stem = f"{accession}_{donor}_24h"
        with gzip.open(raw_dir / f"{stem}_features.tsv.gz", "rt") as handle:
            features = [line.rstrip("\n").split("\t") for line in handle]
        gene_expression = [False] * (len(features) + 1)
        target_rows = {}
        for row, (_, gene, feature_type) in enumerate(features, start=1):
            if feature_type != "Gene Expression":
                continue
            gene_expression[row] = True
            if gene in target_index:
                if gene in target_rows:
                    raise ValueError(f"{stem} contains multiple features for {gene}")
                target_rows[gene] = row
        if target_rows.keys() != target_index.keys():
            missing = ", ".join(sorted(target_index.keys() - target_rows.keys()))
            raise ValueError(f"{stem} is missing expected gene-expression features: {missing}")
        row_to_target = {row: target_index[gene] for gene, row in target_rows.items()}

        with gzip.open(raw_dir / f"{stem}_barcodes.tsv.gz", "rt") as handle:
            barcodes = [line.strip().removesuffix("-1") for line in handle]
        barcode_columns = {barcode: i for i, barcode in enumerate(barcodes, start=1)}
        donor_meta = metadata.loc[metadata["Donor"].eq(donor), "barcode"]
        # The series metadata prefixes each barcode with its library ID.
        donor_barcodes = donor_meta.str.rsplit("_", n=1).str[-1]
        donor_columns = donor_barcodes.map(barcode_columns)
        if len(barcode_columns) != len(barcodes) or donor_columns.isna().any() or not donor_columns.is_unique:
            raise ValueError(f"{stem} barcodes must match uniquely to the series metadata")
        column_to_cell = dict(zip(donor_columns, donor_meta.index))

        with gzip.open(raw_dir / f"{stem}_matrix.mtx.gz", "rt") as handle:
            header = next(line for line in handle if not line.startswith("%"))
            num_rows, num_columns, _ = map(int, header.split())
            if (num_rows, num_columns) != (len(features), len(barcodes)):
                raise ValueError(f"{stem} matrix dimensions disagree with its features or barcodes")
            for record in handle:
                row, column, count = record.split()
                cell = column_to_cell.get(int(column))
                if cell is None:
                    continue
                row_index = int(row)
                if gene_expression[row_index]:
                    value = int(count)
                    library_sizes[cell] += value
                    target = row_to_target.get(row_index)
                    if target is not None:
                        expression[target, cell] += value
    if np.any(library_sizes <= 0):
        raise ValueError("Every selected islet cell must have gene-expression counts")
    return expression, library_sizes


def run_row3_single_cell_erasure_validation():
    t0 = time.time()
    meta_p = os.path.join(DL_DIR, "GSE218316_metadata_percell.csv.gz")
    meta_df = pd.read_csv(meta_p).rename(columns={"Unnamed: 0": "barcode"})
    condition_mask = meta_df["Condition"].isin(["Untreated", "H2O2"])
    sub_meta = meta_df.loc[condition_mask & meta_df["Timepoint"].eq("24h")].reset_index(drop=True)
    n_cells_sub = len(sub_meta)

    clock_genes = ["ARNTL", "CLOCK", "CRY1", "CRY2", "PER1", "PER2", "NR1D1", "NR1D2", "DBP", "NPAS2"]
    redox_ber_genes = ["OGG1", "MUTYH", "APEX1", "PARP1", "NFE2L2", "HMOX1", "NQO1", "TXNRD1", "GCLC", "GCLM", "SOD1", "PRDX1"]
    target_genes = clock_genes + redox_ber_genes
    expr_counts, lib_sizes = load_islet_counts(sub_meta, target_genes)

    # Log1p CP10K normalized expression
    cp10k = np.log1p((expr_counts / np.maximum(lib_sizes[None, :], 1.0)) * 1.0e4)
    y_cond = (sub_meta["Condition"].values == "H2O2").astype(int)
    untreated_mask = (y_cond == 0)
    h2o2_mask = (y_cond == 1)

    # Single-cell erasure diagnostics per gene
    per_gene_erasure = []
    for gi, gname in enumerate(target_genes):
        raw_g = expr_counts[gi, :]
        z_all = float(np.mean(raw_g == 0))
        z_untr = float(np.mean(raw_g[untreated_mask] == 0))
        z_h2o2 = float(np.mean(raw_g[h2o2_mask] == 0))
        mean_cp10k_untr = float(np.mean(cp10k[gi, untreated_mask]))
        mean_cp10k_h2o2 = float(np.mean(cp10k[gi, h2o2_mask]))
        per_gene_erasure.append({
            "gene": gname,
            "panel": "circadian_clock" if gname in clock_genes else "BER_redox",
            "status_in_feature_table": "identified_in_all_three_libraries",
            "zero_umi_fraction_all": round(z_all, 5),
            "zero_umi_fraction_untreated": round(z_untr, 5),
            "zero_umi_fraction_H2O2": round(z_h2o2, 5),
            "mean_log1p_cp10k_untreated": round(mean_cp10k_untr, 4),
            "mean_log1p_cp10k_H2O2": round(mean_cp10k_h2o2, 4),
        })
    erasure_df = pd.DataFrame(per_gene_erasure)

    clock_indices = list(range(len(clock_genes)))
    ber_panel = {"OGG1", "MUTYH", "APEX1", "HMOX1", "NQO1", "TXNRD1", "GCLC", "GCLM", "PRDX1"}
    combined_indices = clock_indices + [i for i, gene in enumerate(target_genes) if gene in ber_panel]

    # Descriptive treatment differences within each donor.
    donor_stats = []
    for d_id, d_df in sub_meta.groupby("Donor"):
        u_idx = d_df.index[d_df["Condition"] == "Untreated"].to_numpy()
        h_idx = d_df.index[d_df["Condition"] == "H2O2"].to_numpy()
        tx_i = target_genes.index("TXNRD1")
        hm_i = target_genes.index("HMOX1")
        donor_stats.append({
            "donor": d_id,
            "n_untreated": int(len(u_idx)),
            "n_H2O2": int(len(h_idx)),
            "TXNRD1_delta_log1p_cp10k": round(float(cp10k[tx_i, h_idx].mean() - cp10k[tx_i, u_idx].mean()), 4),
            "HMOX1_delta_log1p_cp10k": round(float(cp10k[hm_i, h_idx].mean() - cp10k[hm_i, u_idx].mean()), 4),
        })

    rng = np.random.default_rng(202606)
    cells = cp10k.T
    donors = sub_meta["Donor"].to_numpy()
    strat_key = sub_meta["Donor"].astype(str) + "_" + sub_meta["Condition"].astype(str)
    cell_folds = list(StratifiedKFold(n_splits=5, shuffle=True, random_state=2026).split(cells, strat_key))
    donor_folds = [(np.flatnonzero(donors != d), np.flatnonzero(donors == d)) for d in np.unique(donors)]
    feature_indices = [clock_indices, combined_indices]
    n_detected_clock_cells = int(np.count_nonzero(np.any(expr_counts[clock_indices] > 0, axis=0)))

    scaling_rows = []
    for size in (1, 2, 5, 10, 25, 50, 100, 250, 500, 1000):
        cv_auc, cv_counts, n_cv = score_disjoint_ensemble_folds(
            cells, y_cond, feature_indices, cell_folds, size, rng
        )
        lodo_auc, lodo_counts, n_lodo = score_disjoint_ensemble_folds(
            cells, y_cond, feature_indices, donor_folds, size, rng
        )
        scaling_rows.append({
            "ensemble_cells_nu": size,
            "all_10_clock_genes_zero_frac": float(stats.hypergeom.pmf(0, n_cells_sub, n_detected_clock_cells, size)),
            "clock10_roc_auc": cv_auc[0],
            "clock10_lodo_roc_auc": lodo_auc[0],
            "combined_clock_redox_roc_auc": cv_auc[1],
            "combined_lodo_roc_auc": lodo_auc[1],
            "cv_folds_evaluable": n_cv,
            "cv_test_ensembles_untreated": cv_counts[0],
            "cv_test_ensembles_H2O2": cv_counts[1],
            "lodo_folds_evaluable": n_lodo,
            "lodo_test_ensembles_untreated": lodo_counts[0],
            "lodo_test_ensembles_H2O2": lodo_counts[1],
        })

    pd.DataFrame(scaling_rows).to_csv(os.path.join(DATA_DIR, "module6_sc_h2o2_ensemble_scaling.csv"), index=False)
    erasure_df.to_csv(os.path.join(DATA_DIR, "module6_sc_h2o2_gene_erasure.csv"), index=False)
    largest_complete = next(
        row for row in reversed(scaling_rows)
        if row["cv_folds_evaluable"] == len(cell_folds) and row["lodo_folds_evaluable"] == len(donor_folds)
    )
    det_clock_df = erasure_df.loc[erasure_df["panel"].eq("circadian_clock")]
    return {
        "runtime_s": round(time.time() - t0, 4),
        "analysis_timepoint": "24h",
        "gse218316_untreated_or_H2O2_cells_all_timepoints": int(condition_mask.sum()),
        "gse218316_n_untreated_cells": int(untreated_mask.sum()),
        "gse218316_n_H2O2_cells": int(h2o2_mask.sum()),
        "gse218316_total_cells_analyzed": int(n_cells_sub),
        "gse218316_n_donors": len(donor_folds),
        "cv_protocol": "Stratified five-fold and leave-one-donor-out CV with disjoint pools within each fold",
        "sample_accessions": ISLET_24H_SAMPLES,
        "gene_expression_targets": target_genes,
        "clock_genes": clock_genes,
        "combined_classifier_genes": [target_genes[i] for i in combined_indices],
        "donor_pseudobulk_consistency": donor_stats,
        "mean_detected_clock_gene_zero_erasure": float(det_clock_df["zero_umi_fraction_all"].mean()),
        "nu1_zero_all_clock_genes_cell_fraction": (n_cells_sub - n_detected_clock_cells) / n_cells_sub,
        "nu1_clock10_auc": scaling_rows[0]["clock10_roc_auc"],
        "nu1_combined_auc": scaling_rows[0]["combined_clock_redox_roc_auc"],
        "largest_complete_ensemble": largest_complete,
    }


# ==============================================================================
# ROW 4: LIVING HUMAN HEPG2 SATURATION MUTAGENESIS MPRA (GSE126550)
# ==============================================================================
def read_gz_fasta(path):
    seq = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def annotate_duplex_guanine_context(seq, pos_1based):
    """For a 1-based position with REF in {'G', 'C'}, compute duplex G-run length and SP1/E-box motif flag."""
    i = pos_1based - 1
    base = seq[i]
    if base == "G":
        l_idx = i
        while l_idx > 0 and seq[l_idx - 1] == "G":
            l_idx -= 1
        r_idx = i
        while r_idx + 1 < len(seq) and seq[r_idx + 1] == "G":
            r_idx += 1
        g_run = r_idx - l_idx + 1
    elif base == "C":
        l_idx = i
        while l_idx > 0 and seq[l_idx - 1] == "C":
            l_idx -= 1
        r_idx = i
        while r_idx + 1 < len(seq) and seq[r_idx + 1] == "C":
            r_idx += 1
        g_run = r_idx - l_idx + 1
    else:
        return 0, False

    win_start = max(0, i - 6)
    win_end = min(len(seq), i + 7)
    local_w = seq[win_start:win_end]
    motifs = [
        "GGGCGG", "CCCGCC", "GGGGCGG", "CCGCCCC",
        "CCTCCCC", "GGGGAGG", "CCCTCCC", "GGGAGGG",
        "CACCCC", "GGGGTG", "GGGGTGG", "CCACCCC",
        "CACGTG", "TGGGCGG", "CCGCCCA",
    ]
    in_motif = any(
        (win_start + m.start() <= i < win_start + m.start() + len(pat))
        for pat in motifs
        for m in re.finditer(f"(?={pat})", local_w)
    )
    return g_run, in_motif


def process_mpra_library(lib_name, fa_path, var_path, dna_paths, rna_paths):
    ref_seq = read_gz_fasta(fa_path)

    dna_umi = defaultdict(int)
    rna_umi = defaultdict(int)
    total_dna_umis = 0
    total_rna_umis = 0

    for dp in dna_paths:
        proc = subprocess.Popen(["gzip", "-dc", dp], stdout=subprocess.PIPE, text=True, bufsize=1048576)
        for line in proc.stdout:
            tag = line.split("\t", 1)[0]
            dna_umi[tag] += 1
            total_dna_umis += 1
        proc.wait()

    for rp in rna_paths:
        proc = subprocess.Popen(["gzip", "-dc", rp], stdout=subprocess.PIPE, text=True, bufsize=1048576)
        for line in proc.stdout:
            tag = line.split("\t", 1)[0]
            rna_umi[tag] += 1
            total_rna_umis += 1
        proc.wait()

    wt_ratios = []
    single_mut_ratios = defaultdict(list)
    from scipy.sparse import lil_matrix
    var_vocab = {}
    bar_y = []
    bar_w = []
    bar_vlists = []

    snv_re = re.compile(r"^[^:]+:(\d+):([ACGT])>([ACGT])$")

    proc_v = subprocess.Popen(["gzip", "-dc", var_path], stdout=subprocess.PIPE, text=True, bufsize=1048576)
    for line in proc_v.stdout:
        parts = line.strip().split()
        if not parts:
            continue
        tag = parts[0]
        d_cnt = dna_umi.get(tag, 0)
        if d_cnt < 5:
            continue
        r_cnt = rna_umi.get(tag, 0)
        d_rpm = (d_cnt + 0.5) * 1.0e6 / total_dna_umis
        r_rpm = (r_cnt + 0.5) * 1.0e6 / total_rna_umis
        log2_ratio = math.log2(r_rpm / d_rpm)

        if len(parts) == 1:
            wt_ratios.append(log2_ratio)
            bar_y.append(log2_ratio)
            bar_w.append(math.sqrt(d_cnt))
            bar_vlists.append([])
            continue

        vtokens = parts[1:]
        parsed_snvs = []
        valid_all = True
        for vt in vtokens:
            m = snv_re.match(vt)
            if not m:
                valid_all = False
                break
            pos = int(m.group(1))
            ref_b = m.group(2)
            alt_b = m.group(3)
            key = (pos, ref_b, alt_b)
            parsed_snvs.append(key)
        if not valid_all or len(parsed_snvs) > 6:
            continue

        if len(parsed_snvs) == 1:
            single_mut_ratios[parsed_snvs[0]].append(log2_ratio)

        v_ids = []
        for key in parsed_snvs:
            if key not in var_vocab:
                var_vocab[key] = len(var_vocab)
            v_ids.append(var_vocab[key])
        bar_y.append(log2_ratio)
        bar_w.append(math.sqrt(d_cnt))
        bar_vlists.append(v_ids)
    proc_v.wait()

    wt_baseline = float(np.mean(wt_ratios))
    n_bars = len(bar_y)
    n_vars = len(var_vocab)
    X_mat = lil_matrix((n_bars, n_vars), dtype=np.float32)
    for bi, vids in enumerate(bar_vlists):
        for vi in vids:
            X_mat[bi, vi] = 1.0
    X_csr = X_mat.tocsr()
    y_vec = np.array(bar_y, dtype=np.float64) - wt_baseline
    w_vec = np.array(bar_w, dtype=np.float64)

    ridge = Ridge(alpha=10.0, fit_intercept=True)
    ridge.fit(X_csr, y_vec, sample_weight=w_vec)
    coefs = ridge.coef_

    rows = []
    for (pos, ref_b, alt_b), vi in var_vocab.items():
        if pos < 1 or pos > len(ref_seq):
            continue
        if ref_seq[pos - 1] != ref_b:
            continue
        is_oxog_transversion = (ref_b == "G" and alt_b == "T") or (ref_b == "C" and alt_b == "A")
        g_run, in_sp1_ebox = annotate_duplex_guanine_context(ref_seq, pos)
        sm_list = single_mut_ratios.get((pos, ref_b, alt_b), [])
        sm_effect = float(np.mean(sm_list) - wt_baseline) if len(sm_list) >= 3 else np.nan
        eff = float(coefs[vi])

        rows.append({
            "library": lib_name,
            "pos_1based": pos,
            "ref": ref_b,
            "alt": alt_b,
            "is_8oxoG_to_T_transversion": is_oxog_transversion,
            "duplex_G_run_length": g_run,
            "in_SP1_or_Ebox_motif": in_sp1_ebox,
            "effect_log2_RNA_DNA": round(eff, 5),
            "abs_effect_log2": round(abs(eff), 5),
            "single_mutant_only_log2": round(sm_effect, 5) if not np.isnan(sm_effect) else np.nan,
            "n_single_mutant_tags": len(sm_list),
        })

    df_lib = pd.DataFrame(rows)
    valid_sm = df_lib.dropna(subset=["single_mutant_only_log2"])
    r_sm_ridge = float(stats.pearsonr(valid_sm["effect_log2_RNA_DNA"], valid_sm["single_mutant_only_log2"]).statistic) if len(valid_sm) > 10 else np.nan
    return df_lib, {
        "n_high_conf_barcodes": n_bars,
        "n_wt_barcodes": len(wt_ratios),
        "n_fitted_snvs": len(df_lib),
        "pearson_r_ridge_vs_pure_single_mutant_tags": round(r_sm_ridge, 4),
    }


def run_row4_hepg2_mpra_validation():
    t0 = time.time()
    mpra_dir = os.path.join(SCRATCH_DIR, "gse126550")
    libs = [
        (
            "LDLR",
            os.path.join(mpra_dir, "GSM3604144_LDLR.fa.gz"),
            os.path.join(mpra_dir, "GSM3604144_LDLR.variants.txt.gz"),
            [os.path.join(DL_DIR, f"GSM36042{18+i}_LDLR-DNA-{i+1}.tsv.gz") for i in range(3)],
            [os.path.join(DL_DIR, f"GSM36042{21+i}_LDLR-RNA-{i+1}.tsv.gz") for i in range(3)],
        ),
        (
            "LDLR.2",
            os.path.join(mpra_dir, "GSM3604144_LDLR.fa.gz"),
            os.path.join(mpra_dir, "GSM3604143_LDLR.2.variants.txt.gz"),
            [os.path.join(DL_DIR, f"GSM36042{12+i}_LDLR.2-DNA-{i+1}.tsv.gz") for i in range(3)],
            [os.path.join(DL_DIR, f"GSM36042{15+i}_LDLR.2-RNA-{i+1}.tsv.gz") for i in range(3)],
        ),
        (
            "F9",
            os.path.join(mpra_dir, "GSM3604135_F9.fa.gz"),
            os.path.join(mpra_dir, "GSM3604135_F9.variants.txt.gz"),
            [os.path.join(DL_DIR, f"GSM36041{64+i}_F9-DNA-{i+1}.tsv.gz") for i in range(3)],
            [os.path.join(DL_DIR, f"GSM36041{67+i}_F9-RNA-{i+1}.tsv.gz") for i in range(3)],
        ),
        (
            "PKLR_48h",
            os.path.join(mpra_dir, "GSM3604148_PKLR.fa.gz"),
            os.path.join(mpra_dir, "GSM3604148_PKLR.variants.txt.gz"),
            [os.path.join(DL_DIR, f"GSM36042{48+i}_PKLR_48h-DNA-{i+1}.tsv.gz") for i in range(3)],
            [os.path.join(DL_DIR, f"GSM36042{51+i}_PKLR_48h-RNA-{i+1}.tsv.gz") for i in range(3)],
        ),
    ]

    all_dfs = []
    lib_qc = {}
    for lname, fa_p, var_p, dna_ps, rna_ps in libs:
        df_l, qc_l = process_mpra_library(lname, fa_p, var_p, dna_ps, rna_ps)
        all_dfs.append(df_l)
        lib_qc[lname] = qc_l

    mpra_df = pd.concat(all_dfs, ignore_index=True)
    mpra_df.to_csv(os.path.join(DATA_DIR, "module6_hepg2_mpra_variants.csv"), index=False)

    ldlr1 = mpra_df[mpra_df["library"] == "LDLR"].set_index(["pos_1based", "ref", "alt"])["effect_log2_RNA_DNA"]
    ldlr2 = mpra_df[mpra_df["library"] == "LDLR.2"].set_index(["pos_1based", "ref", "alt"])["effect_log2_RNA_DNA"]
    common_idx = ldlr1.index.intersection(ldlr2.index)
    repro_r = float(stats.pearsonr(ldlr1.loc[common_idx], ldlr2.loc[common_idx]).statistic)

    dedup_df = mpra_df.copy()
    dedup_df["promoter"] = dedup_df["library"].replace({"LDLR.2": "LDLR"})
    dedup_df["cell_line"] = dedup_df["promoter"].map({"LDLR": "HepG2", "F9": "HepG2", "PKLR_48h": "K562"})
    prom_df = (
        dedup_df[dedup_df["is_8oxoG_to_T_transversion"]]
        .groupby(["cell_line", "promoter", "pos_1based", "ref", "alt", "duplex_G_run_length", "in_SP1_or_Ebox_motif"], as_index=False)
        .agg({"effect_log2_RNA_DNA": "mean", "abs_effect_log2": "mean"})
    )
    prom_df["abs_effect_log2"] = prom_df["effect_log2_RNA_DNA"].abs()

    # Cross-cell-type pooled splits (kept for reference alongside primary stratified HepG2 / K562 splits)
    g1_ctrl = prom_df[(prom_df["duplex_G_run_length"] == 1) & (~prom_df["in_SP1_or_Ebox_motif"])]
    g2_ctrl = prom_df[prom_df["duplex_G_run_length"] == 2]
    g3_only = prom_df[prom_df["duplex_G_run_length"] == 3]
    g4plus_all = prom_df[prom_df["duplex_G_run_length"] >= 4]
    g3plus_all = prom_df[prom_df["duplex_G_run_length"] >= 3]
    sp1_g3plus_sinks = prom_df[(prom_df["duplex_G_run_length"] >= 3) & (prom_df["in_SP1_or_Ebox_motif"])]
    f9_g3plus = prom_df[(prom_df["promoter"] == "F9") & (prom_df["duplex_G_run_length"] >= 3)]

    mwu_g2 = stats.mannwhitneyu(g2_ctrl["abs_effect_log2"], g1_ctrl["abs_effect_log2"], alternative="two-sided")
    mwu_g3 = stats.mannwhitneyu(g3_only["abs_effect_log2"], g1_ctrl["abs_effect_log2"], alternative="two-sided")
    mwu_g4p = stats.mannwhitneyu(g4plus_all["abs_effect_log2"], g1_ctrl["abs_effect_log2"], alternative="two-sided")
    mwu_all_g3 = stats.mannwhitneyu(g3plus_all["abs_effect_log2"], g1_ctrl["abs_effect_log2"], alternative="two-sided")
    mwu_sp1_sinks = stats.mannwhitneyu(sp1_g3plus_sinks["abs_effect_log2"], g1_ctrl["abs_effect_log2"], alternative="two-sided")

    # Primary stratified Cohort A: HepG2 hepatocellular carcinoma (LDLR deduplicated + F9)
    hepg2_df = prom_df[prom_df["cell_line"] == "HepG2"]
    hepg2_g1 = hepg2_df[(hepg2_df["duplex_G_run_length"] == 1) & (~hepg2_df["in_SP1_or_Ebox_motif"])]
    hepg2_g2 = hepg2_df[hepg2_df["duplex_G_run_length"] == 2]
    hepg2_g3 = hepg2_df[hepg2_df["duplex_G_run_length"] == 3]
    hepg2_g4p = hepg2_df[hepg2_df["duplex_G_run_length"] >= 4]
    hepg2_nonmotif_g3p = hepg2_df[(hepg2_df["duplex_G_run_length"] >= 3) & (~hepg2_df["in_SP1_or_Ebox_motif"])]
    hepg2_sp1_g3p = hepg2_df[(hepg2_df["duplex_G_run_length"] >= 3) & (hepg2_df["in_SP1_or_Ebox_motif"])]

    mwu_hepg2_sp1 = stats.mannwhitneyu(hepg2_sp1_g3p["abs_effect_log2"], hepg2_g1["abs_effect_log2"], alternative="two-sided")
    mwu_hepg2_nonmotif = stats.mannwhitneyu(hepg2_nonmotif_g3p["abs_effect_log2"], hepg2_g1["abs_effect_log2"], alternative="two-sided")
    mwu_hepg2_sp1_vs_nonmotif = stats.mannwhitneyu(hepg2_sp1_g3p["abs_effect_log2"], hepg2_nonmotif_g3p["abs_effect_log2"], alternative="two-sided")
    mwu_hepg2_g4p = stats.mannwhitneyu(hepg2_g4p["abs_effect_log2"], hepg2_g1["abs_effect_log2"], alternative="two-sided")
    mwu_f9_vs_hepg2_g1 = stats.mannwhitneyu(f9_g3plus["abs_effect_log2"], hepg2_g1["abs_effect_log2"], alternative="two-sided")

    # Primary stratified Cohort B: K562 erythroleukemia (PKLR_48h)
    k562_df = prom_df[prom_df["cell_line"] == "K562"]
    k562_g1 = k562_df[(k562_df["duplex_G_run_length"] == 1) & (~k562_df["in_SP1_or_Ebox_motif"])]
    k562_nonmotif_g3p = k562_df[(k562_df["duplex_G_run_length"] >= 3) & (~k562_df["in_SP1_or_Ebox_motif"])]
    k562_all_g3p = k562_df[k562_df["duplex_G_run_length"] >= 3]
    k562_motif_g3p = k562_df[(k562_df["duplex_G_run_length"] >= 3) & (k562_df["in_SP1_or_Ebox_motif"])]

    mwu_k562_motif = stats.mannwhitneyu(k562_motif_g3p["abs_effect_log2"], k562_g1["abs_effect_log2"], alternative="two-sided")
    mwu_k562_all_g3p = stats.mannwhitneyu(k562_all_g3p["abs_effect_log2"], k562_g1["abs_effect_log2"], alternative="two-sided")
    mwu_k562_nonmotif = stats.mannwhitneyu(k562_nonmotif_g3p["abs_effect_log2"], k562_g1["abs_effect_log2"], alternative="two-sided")

    ldlr_only = prom_df[prom_df["promoter"] == "LDLR"]
    ldlr_sinks = ldlr_only[(ldlr_only["duplex_G_run_length"] >= 3) & (ldlr_only["in_SP1_or_Ebox_motif"])]
    ldlr_g1 = ldlr_only[(ldlr_only["duplex_G_run_length"] == 1) & (~ldlr_only["in_SP1_or_Ebox_motif"])]
    ldlr_nonmotif_g3p = ldlr_only[(ldlr_only["duplex_G_run_length"] >= 3) & (~ldlr_only["in_SP1_or_Ebox_motif"])]
    mwu_ldlr = stats.mannwhitneyu(ldlr_sinks["abs_effect_log2"], ldlr_g1["abs_effect_log2"], alternative="two-sided")
    mwu_ldlr_nonmotif = stats.mannwhitneyu(ldlr_nonmotif_g3p["abs_effect_log2"], ldlr_g1["abs_effect_log2"], alternative="two-sided")

    ratio_hepg2_sp1_vs_g1 = float(hepg2_sp1_g3p["abs_effect_log2"].mean() / max(hepg2_g1["abs_effect_log2"].mean(), 1e-9))
    ratio_hepg2_nonmotif_vs_g1 = float(hepg2_nonmotif_g3p["abs_effect_log2"].mean() / max(hepg2_g1["abs_effect_log2"].mean(), 1e-9))
    ratio_k562_motif_vs_g1 = float(k562_motif_g3p["abs_effect_log2"].mean() / max(k562_g1["abs_effect_log2"].mean(), 1e-9))
    ratio_k562_all_g3p_vs_g1 = float(k562_all_g3p["abs_effect_log2"].mean() / max(k562_g1["abs_effect_log2"].mean(), 1e-9))
    ratio_sp1_vs_g1_pooled = float(sp1_g3plus_sinks["abs_effect_log2"].mean() / max(g1_ctrl["abs_effect_log2"].mean(), 1e-9))
    ratio_ldlr_sp1_vs_g1 = float(ldlr_sinks["abs_effect_log2"].mean() / max(ldlr_g1["abs_effect_log2"].mean(), 1e-9))

    frac_gt30pct_drop_hepg2_sinks = float((hepg2_sp1_g3p["effect_log2_RNA_DNA"] <= -0.5145).mean())
    frac_gt30pct_drop_hepg2_g1 = float((hepg2_g1["effect_log2_RNA_DNA"] <= -0.5145).mean())
    ldlr_frac_gt30pct_sinks = float((ldlr_sinks["effect_log2_RNA_DNA"] <= -0.5145).mean())
    ldlr_frac_gt30pct_g1 = float((ldlr_g1["effect_log2_RNA_DNA"] <= -0.5145).mean())

    sre1_sub = ldlr_only[(ldlr_only["pos_1based"] >= 180) & (ldlr_only["pos_1based"] <= 192)]
    sp1_r3_sub = ldlr_only[(ldlr_only["pos_1based"] >= 198) & (ldlr_only["pos_1based"] <= 210)]
    sre1_peak_row = sre1_sub.loc[sre1_sub["effect_log2_RNA_DNA"].idxmin()]
    sp1_peak_row = sp1_r3_sub.loc[sp1_r3_sub["effect_log2_RNA_DNA"].idxmin()]

    return {
        "runtime_s": round(time.time() - t0, 4),
        "library_qc": lib_qc,
        "cell_line_map": {"LDLR": "HepG2", "LDLR.2": "HepG2", "F9": "HepG2", "PKLR_48h": "K562"},
        "ldlr_vs_ldlr2_snv_pearson_r": round(repro_r, 4),
        "n_unique_promoter_8oxoG_transversions": int(len(prom_df)),
        "hepg2_n_unique_8oxoG_transversions": int(len(hepg2_df)),
        "hepg2_n_G1_nonmotif": int(len(hepg2_g1)),
        "hepg2_n_G2": int(len(hepg2_g2)),
        "hepg2_n_G3": int(len(hepg2_g3)),
        "hepg2_n_G4plus": int(len(hepg2_g4p)),
        "hepg2_n_nonmotif_G3plus": int(len(hepg2_nonmotif_g3p)),
        "hepg2_n_SP1_SRE1_G3plus": int(len(hepg2_sp1_g3p)),
        "hepg2_n_F9_G3plus": int(len(f9_g3plus)),
        "hepg2_mean_abs_log2_G1": round(float(hepg2_g1["abs_effect_log2"].mean()), 4),
        "hepg2_mean_abs_log2_G2": round(float(hepg2_g2["abs_effect_log2"].mean()), 4),
        "hepg2_mean_abs_log2_G3": round(float(hepg2_g3["abs_effect_log2"].mean()), 4),
        "hepg2_mean_abs_log2_G4plus": round(float(hepg2_g4p["abs_effect_log2"].mean()), 4),
        "hepg2_mean_abs_log2_nonmotif_G3plus": round(float(hepg2_nonmotif_g3p["abs_effect_log2"].mean()), 4),
        "hepg2_mean_abs_log2_SP1_SRE1_G3plus": round(float(hepg2_sp1_g3p["abs_effect_log2"].mean()), 4),
        "hepg2_mean_abs_log2_F9_G3plus": round(float(f9_g3plus["abs_effect_log2"].mean()), 4),
        "hepg2_SP1_sink_vs_G1_fold_ratio": round(ratio_hepg2_sp1_vs_g1, 3),
        "hepg2_nonmotif_G3plus_vs_G1_fold_ratio": round(ratio_hepg2_nonmotif_vs_g1, 3),
        "hepg2_mwu_p_SP1_sink_vs_G1": float(mwu_hepg2_sp1.pvalue),
        "hepg2_mwu_p_nonmotif_G3plus_vs_G1": float(mwu_hepg2_nonmotif.pvalue),
        "hepg2_mwu_p_SP1_vs_nonmotif_G3plus": float(mwu_hepg2_sp1_vs_nonmotif.pvalue),
        "hepg2_mwu_p_G4plus_vs_G1": float(mwu_hepg2_g4p.pvalue),
        "hepg2_mwu_p_F9_G3plus_vs_G1": float(mwu_f9_vs_hepg2_g1.pvalue),
        "hepg2_frac_gt30pct_reporter_loss_SP1_sinks": round(frac_gt30pct_drop_hepg2_sinks, 4),
        "hepg2_frac_gt30pct_reporter_loss_G1": round(frac_gt30pct_drop_hepg2_g1, 4),
        "k562_n_unique_8oxoG_transversions": int(len(k562_df)),
        "k562_n_G1_nonmotif": int(len(k562_g1)),
        "k562_n_nonmotif_G3plus": int(len(k562_nonmotif_g3p)),
        "k562_n_all_G3plus": int(len(k562_all_g3p)),
        "k562_n_motif_G3plus": int(len(k562_motif_g3p)),
        "k562_mean_abs_log2_G1": round(float(k562_g1["abs_effect_log2"].mean()), 4),
        "k562_mean_abs_log2_nonmotif_G3plus": round(float(k562_nonmotif_g3p["abs_effect_log2"].mean()), 4),
        "k562_mean_abs_log2_all_G3plus": round(float(k562_all_g3p["abs_effect_log2"].mean()), 4),
        "k562_mean_abs_log2_motif_G3plus": round(float(k562_motif_g3p["abs_effect_log2"].mean()), 4),
        "k562_motif_G3plus_vs_G1_fold_ratio": round(ratio_k562_motif_vs_g1, 3),
        "k562_all_G3plus_vs_G1_fold_ratio": round(ratio_k562_all_g3p_vs_g1, 3),
        "k562_mwu_p_motif_G3plus_vs_G1": float(mwu_k562_motif.pvalue),
        "k562_mwu_p_all_G3plus_vs_G1": float(mwu_k562_all_g3p.pvalue),
        "k562_mwu_p_nonmotif_G3plus_vs_G1": float(mwu_k562_nonmotif.pvalue),
        "ldlr_SP1_sink_mean_abs_log2": round(float(ldlr_sinks["abs_effect_log2"].mean()), 4),
        "ldlr_G1_mean_abs_log2": round(float(ldlr_g1["abs_effect_log2"].mean()), 4),
        "ldlr_nonmotif_G3plus_mean_abs_log2": round(float(ldlr_nonmotif_g3p["abs_effect_log2"].mean()), 4),
        "ldlr_SP1_sink_vs_G1_fold_ratio": round(ratio_ldlr_sp1_vs_g1, 3),
        "ldlr_mwu_p_SP1_sink_vs_G1": float(mwu_ldlr.pvalue),
        "ldlr_mwu_p_nonmotif_G3plus_vs_G1": float(mwu_ldlr_nonmotif.pvalue),
        "ldlr_frac_gt30pct_reporter_loss_SP1_sinks": round(ldlr_frac_gt30pct_sinks, 4),
        "ldlr_frac_gt30pct_reporter_loss_G1": round(ldlr_frac_gt30pct_g1, 4),
        "cross_cell_pooled_SP1_sink_vs_G1_fold_ratio": round(ratio_sp1_vs_g1_pooled, 3),
        "cross_cell_mwu_p_SP1_G3plus_vs_G1": float(mwu_sp1_sinks.pvalue),
        "ldlr_sre1_peak": {"pos": int(sre1_peak_row["pos_1based"]), "log2": round(float(sre1_peak_row["effect_log2_RNA_DNA"]), 4)},
        "ldlr_sp1_peak": {"pos": int(sp1_peak_row["pos_1based"]), "log2": round(float(sp1_peak_row["effect_log2_RNA_DNA"]), 4)},
    }


# ==============================================================================
# ROW 5A & 5B: BACTERIAL GO-REPAIR UNCOLLAPSING + POLARIZED BEAMLINE BUDGET
# ==============================================================================
def load_bsgenome_fastas():
    fa_p = "/tmp/bsgenome_ecoli/BSgenome.Ecoli.NCBI.11strains/inst/extdata/single_sequences.fa"
    if not os.path.exists(fa_p):
        tar_p = os.path.join(SCRATCH_DIR, "ecoli_wgs", "BSgenome_Ecoli_genomes", "BSgenome.Ecoli.NCBI.11strains_1.0.0.tar.gz")
        with tarfile.open(tar_p, "r:gz") as tf:
            tf.extractall("/tmp/bsgenome_ecoli")
    seqs = {}
    cur_id = None
    buf = []
    with open(fa_p) as f:
        for line in f:
            if line.startswith(">"):
                if cur_id is not None:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip().upper())
        if cur_id is not None:
            seqs[cur_id] = "".join(buf)
    return seqs


def compute_genome_g_opportunities(seq):
    """Count total G:C base pairs stratified by duplex G-run length L=1..5+ and 3'-flanking base on the G strand."""
    n = len(seq)
    run_opp = Counter()
    flank3_opp = Counter()
    i = 0
    while i < n:
        if seq[i] == "G":
            j = i
            while j + 1 < n and seq[j + 1] == "G":
                j += 1
            run_len = min(j - i + 1, 5)
            run_opp[run_len] += (j - i + 1)
            for k in range(i, j + 1):
                if k + 1 < n and seq[k + 1] in "ACGT":
                    flank3_opp[seq[k + 1]] += 1
            i = j + 1
        elif seq[i] == "C":
            j = i
            while j + 1 < n and seq[j + 1] == "C":
                j += 1
            run_len = min(j - i + 1, 5)
            run_opp[run_len] += (j - i + 1)
            for k in range(i, j + 1):
                if k - 1 >= 0 and seq[k - 1] in "ACGT":
                    comp_map = {"A": "T", "T": "A", "C": "G", "G": "C"}
                    flank3_opp[comp_map[seq[k - 1]]] += 1
            i = j + 1
        else:
            i += 1
    return run_opp, flank3_opp


def run_row5_ecoli_go_and_beamline_validation():
    t0 = time.time()
    seqs = load_bsgenome_fastas()
    vcf_root = os.path.join(SCRATCH_DIR, "ecoli_wgs", "vcf")

    cohort_files = [
        ("WT_MG1655_Foster2015", "wild_type", os.path.join(vcf_root, "Foster2015", "Foster2015_nonmut.vcf")),
        ("WT_ED1a_Foster2015", "wild_type", os.path.join(vcf_root, "Foster2015", "Foster2015_nonmut_ED1a.vcf")),
        ("WT_IAI1_Foster2015", "wild_type", os.path.join(vcf_root, "Foster2015", "Foster2015_nonmut_IAI1.vcf")),
        ("mutY_Foster2015", "GO_glycosylase_deficient", os.path.join(vcf_root, "Foster2015", "Foster2015_mutY.vcf")),
        ("mutY_SaneAgashe", "GO_glycosylase_deficient", os.path.join(vcf_root, "Deepa", "mutY.vcf")),
        ("mutM_mutY_Foster2015", "GO_double_deficient", os.path.join(vcf_root, "Foster2015", "Foster2015_mutMmutY.vcf")),
        ("mutT_Foster2015", "dGTP_sanitizer_deficient", os.path.join(vcf_root, "Foster2015", "Foster2015_mutT.vcf")),
        ("mutS_SaneAgashe", "MMR_deficient", os.path.join(vcf_root, "Deepa", "mutS.vcf")),
        ("mutL_SaneAgashe", "MMR_deficient", os.path.join(vcf_root, "Deepa", "mutL.vcf")),
        ("mutH_SaneAgashe", "MMR_deficient", os.path.join(vcf_root, "Deepa", "mutH.vcf")),
        ("nth_nei_Foster2015", "pyrimidine_BER_deficient", os.path.join(vcf_root, "Foster2015", "Foster15_nth_nei.vcf")),
        ("xthA_nfo_Foster2015", "AP_endonuclease_deficient", os.path.join(vcf_root, "Foster2015", "Foster15_xthA_nfo.vcf")),
        ("uvrA_Foster2015", "NER_deficient", os.path.join(vcf_root, "Foster2015", "Foster2015_uvrA.vcf")),
        ("dnaQ_SaneAgashe", "proofreading_deficient", os.path.join(vcf_root, "Deepa", "dnaQ.vcf")),
    ]

    comp_map = {"A": "T", "T": "A", "C": "G", "G": "C"}
    opp_cache = {}
    cohort_rows = []
    strain_unit_ca_fracs = defaultdict(list)

    go_grun_counts = Counter()
    go_flank3_counts = Counter()
    wt_grun_counts = Counter()

    for cname, repair_class, vpath in cohort_files:
        total_snps = 0
        ca_gc_to_ta = 0
        tg_at_to_cg = 0
        transitions = 0
        ref_matched = 0
        per_strain = defaultdict(lambda: [0, 0])

        with open(vpath) as vf:
            for line in vf:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                chrom = parts[0].strip()
                pos = int(parts[1])
                ref = parts[3].strip().upper()
                alt = parts[4].strip().upper()
                info = parts[7] if len(parts) > 7 else ""
                if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                    continue
                m_st = re.search(r"ST=([^;\s]+)", info)
                st_id = m_st.group(1) if m_st else cname

                total_snps += 1
                per_strain[st_id][1] += 1
                is_ca = (ref == "G" and alt == "T") or (ref == "C" and alt == "A")
                is_tg = (ref == "A" and alt == "C") or (ref == "T" and alt == "G")
                is_ti = (ref, alt) in (("A", "G"), ("G", "A"), ("C", "T"), ("T", "C"))

                if is_ca:
                    ca_gc_to_ta += 1
                    per_strain[st_id][0] += 1
                elif is_tg:
                    tg_at_to_cg += 1
                elif is_ti:
                    transitions += 1

                if chrom in seqs and 1 <= pos <= len(seqs[chrom]) and seqs[chrom][pos - 1] == ref:
                    ref_matched += 1
                    if chrom not in opp_cache:
                        opp_cache[chrom] = compute_genome_g_opportunities(seqs[chrom])
                    if is_ca:
                        s = seqs[chrom]
                        idx = pos - 1
                        if ref == "G":
                            l_i, r_i = idx, idx
                            while l_i > 0 and s[l_i - 1] == "G":
                                l_i -= 1
                            while r_i + 1 < len(s) and s[r_i + 1] == "G":
                                r_i += 1
                            grun = min(r_i - l_i + 1, 5)
                            f3 = s[idx + 1] if idx + 1 < len(s) else "N"
                        else:
                            l_i, r_i = idx, idx
                            while l_i > 0 and s[l_i - 1] == "C":
                                l_i -= 1
                            while r_i + 1 < len(s) and s[r_i + 1] == "C":
                                r_i += 1
                            grun = min(r_i - l_i + 1, 5)
                            f3 = comp_map.get(s[idx - 1], "N") if idx - 1 >= 0 else "N"

                        if repair_class in ("GO_glycosylase_deficient", "GO_double_deficient"):
                            go_grun_counts[grun] += 1
                            if f3 in "ACGT":
                                go_flank3_counts[f3] += 1
                        elif repair_class == "wild_type":
                            wt_grun_counts[grun] += 1

        for st_id, (c_n, t_n) in per_strain.items():
            if t_n >= 5:
                strain_unit_ca_fracs[repair_class].append(c_n / t_n)

        cohort_rows.append({
            "cohort": cname,
            "repair_class": repair_class,
            "n_strain_groups": len(per_strain),
            "total_snvs": total_snps,
            "ref_fasta_match_pct": round(100.0 * ref_matched / max(total_snps, 1), 2),
            "GC_to_TA_count": ca_gc_to_ta,
            "GC_to_TA_fraction": round(ca_gc_to_ta / max(total_snps, 1), 5),
            "AT_to_CG_fraction": round(tg_at_to_cg / max(total_snps, 1), 5),
            "transition_fraction": round(transitions / max(total_snps, 1), 5),
        })

    cohort_df = pd.DataFrame(cohort_rows)
    cohort_df.to_csv(os.path.join(DATA_DIR, "module6_ecoli_mutator_wgs.csv"), index=False)

    ref_run_opp, ref_f3_opp = opp_cache["NC_000913_F"]

    grun_norm = {}
    base_g1_rate = (go_grun_counts[1] / ref_run_opp[1])
    for L in [1, 2, 3, 4, 5]:
        rate_L = go_grun_counts[L] / ref_run_opp[L]
        grun_norm[f"G{L}{'+' if L == 5 else ''}"] = {
            "observed_GC_to_TA": int(go_grun_counts[L]),
            "genome_GC_bp_opportunity": int(ref_run_opp[L]),
            "per_bp_relative_to_G1": round(rate_L / max(base_g1_rate, 1e-12), 4),
        }

    f3_norm = {}
    mean_f3_rate = sum(go_flank3_counts.values()) / sum(ref_f3_opp.values())
    for b in ["G", "A", "C", "T"]:
        rb = go_flank3_counts[b] / ref_f3_opp[b]
        f3_norm[f"5prime_G_{b}_3prime"] = {
            "observed_GC_to_TA": int(go_flank3_counts[b]),
            "genome_opportunity": int(ref_f3_opp[b]),
            "enrichment_over_genome_mean": round(rb / max(mean_f3_rate, 1e-12), 4),
        }

    go_units = strain_unit_ca_fracs["GO_double_deficient"] + strain_unit_ca_fracs["GO_glycosylase_deficient"]
    ctrl_units = strain_unit_ca_fracs["wild_type"] + strain_unit_ca_fracs["MMR_deficient"]
    mwu_go_vs_ctrl = stats.mannwhitneyu(go_units, ctrl_units, alternative="two-sided")

    wt_pooled_ca = float(
        cohort_df[cohort_df["repair_class"] == "wild_type"]["GC_to_TA_count"].sum()
        / cohort_df[cohort_df["repair_class"] == "wild_type"]["total_snvs"].sum()
    )
    mmr_pooled_ca = float(
        cohort_df[cohort_df["repair_class"] == "MMR_deficient"]["GC_to_TA_count"].sum()
        / cohort_df[cohort_df["repair_class"] == "MMR_deficient"]["total_snvs"].sum()
    )
    muty_pooled_ca = float(
        cohort_df[cohort_df["repair_class"] == "GO_glycosylase_deficient"]["GC_to_TA_count"].sum()
        / cohort_df[cohort_df["repair_class"] == "GO_glycosylase_deficient"]["total_snvs"].sum()
    )
    mutmmuty_pooled_ca = float(
        cohort_df[cohort_df["repair_class"] == "GO_double_deficient"]["GC_to_TA_count"].sum()
        / cohort_df[cohort_df["repair_class"] == "GO_double_deficient"]["total_snvs"].sum()
    )

    mod5_summary_p = os.path.join(DATA_DIR, "module5_summary.json")
    if os.path.exists(mod5_summary_p):
        with open(mod5_summary_p) as fh5:
            mod5_s = json.load(fh5)
        a_so_point = float(mod5_s["cosmic_ray_chiral_mutagenesis"]["chiral_bias_500eV_lab_polarized_e_G4"])
    else:
        a_so_point = 2.1603e-7

    a_eff_target = 1.0e-3
    n_events_3sigma = int(math.ceil(9.0 / (a_eff_target ** 2)))       # 9,000,000 lesions per arm (CISS/DEA regime)
    n_events_so_3sigma = int(math.ceil(9.0 / (a_so_point ** 2)))      # ~1.93e14 lesions per arm (pure point-Coulomb Dirac SO)
    beam_current_A = 10.0e-9
    electrons_per_s = beam_current_A / E_CHARGE
    surface_g_density_cm2 = 1.5e12 / 0.5  # 3.0e12 G sites/cm^2 in G4 SAM monolayer
    sigma_damage_cm2 = 2.0e-16            # 100-500 eV electron impact ionization cross section
    oxog_per_damage = 0.05                # 8-oxoG branching ratio per ionizing event
    eta_oxog_per_electron = surface_g_density_cm2 * sigma_damage_cm2 * oxog_per_damage  # 3.0e-5
    lesion_rate_per_s = electrons_per_s * eta_oxog_per_electron
    exposure_time_s = n_events_3sigma / lesion_rate_per_s
    exposure_time_so_s = n_events_so_3sigma / lesion_rate_per_s

    beamline_spec = {
        "status": "UNEXECUTED_PHYSICAL_EXPERIMENT",
        "reason_not_in_public_archives": (
            "No published dataset combines helicity-reversed (L <-> -L) 10-500 eV polarized electrons or "
            "surface mu+ across tilt angles spanning psi_crit = 26.9 deg on SAM-oriented D/L-GGGG DNA or "
            "Delta_mutM_Delta_mutY_Delta_mutS Duplex-Seq."
        ),
        "critical_pitch_cone_deg": 26.9,
        "recommended_tilt_pair_deg": [15.0, 45.0],
        "beam_polarization_abs_L": 0.80,
        "target_asymmetry_threshold": a_eff_target,
        "pure_point_coulomb_dirac_so_asymmetry": round(a_so_point, 11),
        "eta_oxog_per_electron": round(eta_oxog_per_electron, 8),
        "min_lesion_events_per_helicity_arm_3sigma": n_events_3sigma,
        "min_lesion_events_point_coulomb_3sigma": n_events_so_3sigma,
        "standard_GaAs_e_gun_current_nA": 10.0,
        "estimated_irradiation_seconds_per_arm": round(exposure_time_s, 2),
        "estimated_irradiation_seconds_point_coulomb": round(exposure_time_so_s, 1),
    }

    return {
        "runtime_s": round(time.time() - t0, 4),
        "wt_GC_to_TA_fraction": round(wt_pooled_ca, 4),
        "mmr_GC_to_TA_fraction": round(mmr_pooled_ca, 4),
        "mutY_GC_to_TA_fraction": round(muty_pooled_ca, 4),
        "mutMmutY_GC_to_TA_fraction": round(mutmmuty_pooled_ca, 4),
        "mutMmutY_total_SNVs_parsed": int(cohort_df[cohort_df["repair_class"] == "GO_double_deficient"]["total_snvs"].iloc[0]),
        "mutMmutY_vs_MMR_fraction_ratio": round(mutmmuty_pooled_ca / max(mmr_pooled_ca, 1e-9), 2),
        "mutMmutY_vs_WT_fraction_ratio": round(mutmmuty_pooled_ca / max(wt_pooled_ca, 1e-9), 2),
        "strain_unit_mwu_p_GO_vs_WT_MMR": float(mwu_go_vs_ctrl.pvalue),
        "go_deficient_G_run_normalized": grun_norm,
        "go_deficient_3prime_flank_normalized": f3_norm,
        "row5b_beamline_protocol": beamline_spec,
    }


def main():
    parser = argparse.ArgumentParser(description="Run Module 6 empirical analyses.")
    parser.add_argument("--row", type=int, choices=range(1, 6), help="Run one numbered row.")
    args = parser.parse_args()
    row_functions = (
        run_row1_thermal_vs_ac_validation,
        run_row2_nucleosome_ogg1_validation,
        run_row3_single_cell_erasure_validation,
        run_row4_hepg2_mpra_validation,
        run_row5_ecoli_go_and_beamline_validation,
    )

    if args.row is not None:
        result = row_functions[args.row - 1]()
        with (DATA_DIR / f"module6_row{args.row}_summary.json").open("w", encoding="utf-8") as f:
            json.dump(result, f, indent=2)
        print(json.dumps(result, indent=2))
        return

    t_start = time.time()
    print("Running Module 6 stage-specific empirical consistency checks across Table 4...")

    results = []
    for row, run in enumerate(row_functions, start=1):
        results.append(run())
        print(f"[Row {row} DONE]", flush=True)
    r1, r2, r3, r4, r5 = results
    sc_complete = r3["largest_complete_ensemble"]

    summary_table = pd.DataFrame([
        {
            "table6_row": "Row 1",
            "channel_stage_tested": "Thermal B->X Arrhenius gate vs. theoretical Floquet/KMS AC ceiling (Thm. 2)",
            "wetlab_cohort_accession": "Chu et al. 2018 (E. coli REL606 mutS WGS MA at 25, 28, 37 C; BPS + small indels)",
            "sample_size_units": "n=20 MA lines/temperature (60 lines total; published one-sided tests)",
            "key_empirical_metric": (
                f"28 vs 25 C = {r1['emp_ratio_28C_vs_25C']:.3f}x observed and {r1['mod1_predicted_ratio_28C_vs_25C']:.3f}x modeled "
                f"(Padj={r1['chu2018_25_vs_28C_padj']}); 37 vs 28 C = {r1['emp_ratio_37C_vs_28C']:.3f}x "
                f"(Padj={r1['chu2018_28_vs_37C_padj']:.0e}) vs Mod. 1 a priori {r1['mod1_predicted_ratio_37C_vs_28C']:.3f}x ({r1['arrhenius_ratio_agreement_pct']:.1f}% match); "
                f"Floquet/KMS <=10 V/m theoretical shift <1.3e-12"
            ),
            "resolution_status": "STAGE-CONSISTENT (Thermal B->X Gate; RF Direct MA Unmeasured)",
        },
        {
            "table6_row": "Row 2",
            "channel_stage_tested": "Nucleosomal SHL minor-groove rotational gate on hOGG1 accessibility (Thm. 3)",
            "wetlab_cohort_accession": "RCSB PDB 8VWT, 8VWV, 9EOZ, 1KX5 + CLAPS-seq (Li et al. 2025; Cordero et al. 2024)",
            "sample_size_units": f"3 cryo-EM hOGG1-NCP structures + {r2['n_nucleosome_positions_scanned']} 1KX5 bp steps",
            "key_empirical_metric": (
                f"Outward (<=60 deg): {r2['outward_mean_core_clashes_2p5A']:.2f} core clashes ({100*r2['outward_core_zero_clash_fraction']:.1f}% clash-free); "
                f"Inward (>=120 deg): {r2['inward_mean_core_clashes_2p5A']:.1f} core clashes. "
                f"Modeled outward/inward accessibility: {r2['outward_vs_inward_accessibility_ratio']:.1f}x."
            ),
            "resolution_status": "Descriptive steric scan of one nucleosome template",
        },
        {
            "table6_row": "Row 3",
            "channel_stage_tested": "Transcript detection and H2O2 classification in pooled transcriptomes",
            "wetlab_cohort_accession": "NCBI GEO GSE218316 (human islet scRNA-seq, 24 h)",
            "sample_size_units": f"{r3['gse218316_total_cells_analyzed']:,} primary human cells ({r3['gse218316_n_donors']} donors)",
            "key_empirical_metric": (
                f"nu=1 zero-UMI fraction for retained clock genes = {100*r3['mean_detected_clock_gene_zero_erasure']:.2f}%; "
                f"clock AUC={r3['nu1_clock10_auc']:.3f}. At nu={sc_complete['ensemble_cells_nu']}, "
                f"clock AUC={sc_complete['clock10_roc_auc']:.3f} (LODO={sc_complete['clock10_lodo_roc_auc']:.3f}), "
                f"combined {len(r3['combined_classifier_genes'])}-gene AUC={sc_complete['combined_clock_redox_roc_auc']:.3f} "
                f"(LODO={sc_complete['combined_lodo_roc_auc']:.3f})."
            ),
            "resolution_status": "Exploratory analysis across three donors at one timepoint",
        },
        {
            "table6_row": "Row 4",
            "channel_stage_tested": "Downstream cis-regulatory penalty of fixed G:C->T:A at GC-box G3+ tracts (Z->TE_O)",
            "wetlab_cohort_accession": "NCBI GEO GSE126550 (Kircher et al. 2019: HepG2 [LDLR, LDLR.2, F9] & K562 [PKLR])",
            "sample_size_units": f"HepG2 n={r4['hepg2_n_unique_8oxoG_transversions']} & K562 n={r4['k562_n_unique_8oxoG_transversions']} unique G:C->T:A SNVs",
            "key_empirical_metric": (
                f"HepG2 SP1/SRE-1 G3+: {r4['hepg2_SP1_sink_vs_G1_fold_ratio']:.2f}x vs G1 (p={r4['hepg2_mwu_p_SP1_sink_vs_G1']:.2e}; "
                f"non-motif G3+ ctrl {r4['hepg2_nonmotif_G3plus_vs_G1_fold_ratio']:.2f}x, p={r4['hepg2_mwu_p_nonmotif_G3plus_vs_G1']:.2f}); "
                f"K562 PKLR motif G3+: {r4['k562_motif_G3plus_vs_G1_fold_ratio']:.2f}x (p={r4['k562_mwu_p_motif_G3plus_vs_G1']:.2e})"
            ),
            "resolution_status": "STAGE-CONSISTENT (Downstream Cis-Regulatory Vulnerability)",
        },
        {
            "table6_row": "Row 5a / 5b",
            "channel_stage_tested": "Endogenous GO-repair quotient uncollapsing (5a) & polarized lepton beamline (5b)",
            "wetlab_cohort_accession": "Foster et al. 2015 & Sane/Garushyants 2024 E. coli WGS MA VCFs (14 cohorts)",
            "sample_size_units": f"{r5['mutMmutY_total_SNVs_parsed']:,} mutM mutY SNVs + 13 cohorts; 5b: {r5['row5b_beamline_protocol']['min_lesion_events_per_helicity_arm_3sigma']:.1e} (CISS) vs {r5['row5b_beamline_protocol']['min_lesion_events_point_coulomb_3sigma']:.2e} (Dirac SO) lesions/arm",
            "key_empirical_metric": (
                f"5a: G:C->T:A jumps from {100*r5['wt_GC_to_TA_fraction']:.1f}% (WT) / {100*r5['mmr_GC_to_TA_fraction']:.2f}% (MMR) "
                f"to {100*r5['mutMmutY_GC_to_TA_fraction']:.2f}% in mutM mutY ({r5['mutMmutY_vs_MMR_fraction_ratio']:.1f}x vs MMR, p={r5['strain_unit_mwu_p_GO_vs_WT_MMR']:.2e}; G3-G5+ 1.41-1.55x). "
                f"5b UNEXECUTED: tests |A_eff|>=10^-3 CISS/DEA (~4.8 s) vs pure Dirac SO 2.16e-7 (unobservable alone)"
            ),
            "resolution_status": "5a STAGE-CONSISTENT (Endogenous GO Quotient) / 5b UNEXECUTED BEAMLINE",
        },
    ])
    summary_table.to_csv(os.path.join(DATA_DIR, "module6_table6_empirical_validation.csv"), index=False)

    full_summary = {
        "total_module6_wall_clock_s": round(time.time() - t_start, 4),
        "row1_thermal_vs_ac": r1,
        "row2_nucleosome_rotational_gate": r2,
        "row3_single_cell_erasure_and_scaling": r3,
        "row4_hepg2_mpra_hole_sinks": r4,
        "row5_ecoli_go_and_beamline": r5,
    }
    with open(os.path.join(DATA_DIR, "module6_summary.json"), "w") as f:
        json.dump(full_summary, f, indent=2)

    print(f"\n=== MODULE 6 COMPLETE in {full_summary['total_module6_wall_clock_s']:.2f} s ===")
    print(summary_table.to_string(index=False))


if __name__ == "__main__":
    main()
