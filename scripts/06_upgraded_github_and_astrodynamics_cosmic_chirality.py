#!/usr/bin/env python3
"""Module 5: Forensic Audit & Upgrade of Published Quantum-DNA Routines,
Earth-Sun Astrodynamic Time-Encoding via Chiral Solar Neutrinos vs.
Globus-Blandford Weak-Decay Spin-Polarized Leptons, and 240,000-Trajectory
3D Monte Carlo Helical DNA Chiral Mutagenesis.

Eliminates legacy 34 GHz / CMB placeholder assumptions and rigorously evaluates:
  1. Large-scale E. coli REL606 benchmark (2,000 non-overlapping 150-bp genomic
     windows = 1,000 coding CDS + 1,000 intergenic windows = 300 kbp) auditing
     published GitHub routines against trace-preserving Bloch (Tr rho = 1),
     2-qubit dinucleotide mutual information, KMS asymmetric proton transfer,
     and 150x150 tight-binding pi-stack HOMO hole-localization models.
  2. Exact 8,766-hour (1-hour step over 365.25 days) J2000 Earth-Sun Keplerian
     orbital + 24-hour diurnal astrodynamic ephemeris applied to:
       (a) Direct Chiral Solar Neutrinos (pp, monochromatic 7Be 861.8 keV,
           monochromatic pep 1,442.2 keV, 8B): 2D quadrature orbital clock
           combining (a/r)^2 flux modulation (+/-3.34%, Borexino >5-sigma) with
           90-degree phase-shifted radial Doppler shift (+/-1.43 eV on 7Be) +
           diurnal rotational Doppler (+/-1.34 eV), audited against Standard
           Model electroweak (nu-e MSW elastic + CEvNS + Stodolsky WNC) rates.
       (b) Solar-Modulated Spin-Polarized Charged Leptons (Globus & Blandford
           2020, ApJL 895:L11; arXiv:2002.12138v2): weak pion decay
           (pi+ -> mu+ + nu_mu) transducing left-handed neutrino chirality into
           spin-polarized mu+/e+ electromagnetic ionization (~10^27x larger
           cross-section), modulated across diurnal solar time (24 h), orbital
           Compton-Getting Doppler dipole, annual stratospheric temperature wave
           (365.25 d), 27-day Carrington solar rotation, 11-year solar cycle,
           and Young Active Sun solar energetic particles (SEPs).
  3. Explicit 240,000-trajectory 3D Monte Carlo helical B-DNA electrodynamic
     scattering simulation (60,000 trajectories x 4 ensembles x 40 steps =
     9.6e6 3D force evaluations) verifying convergence of inward central-core
     deflection in D-DNA (m = +2) vs outward backbone deflection in L-DNA
     (m = -2), mineral-surface orientation amplification, G4 HOMO hole-sink
     capture (95.66%, 17.95x), prebiotic homochiralization, and repair-deficient
     (Delta mutM Delta mutY Delta mutS) lab beam falsification power.
"""

import glob
import json
import os
import time
import numpy as np
import pandas as pd
import scipy.constants as const
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from scipy.stats import mannwhitneyu

from analysis_paths import DATA_DIR as OUT_DATA, RAW_DIR

os.makedirs(OUT_DATA, exist_ok=True)

HBAR_EV_S = const.hbar / const.e          # 6.582119569e-16 eV*s
HBAR_C_MEV_FM = 197.3269804               # MeV*fm
K_B_EV_K = const.k / const.e              # 8.617333262e-5 eV/K
TEMP_K = 310.15                           # 37.0 C physiological temperature
KBT_EV = K_B_EV_K * TEMP_K                # 0.02672666 eV
C_LIGHT = const.c                         # 299792458.0 m/s
ALPHA_FS = const.alpha                    # 1 / 137.035999084
GF_GEV2 = 1.1663787e-5                    # Fermi weak constant (GeV^-2)
SIN2_THETA_W = 0.23122                    # Weak mixing angle (MS-bar at m_Z)


# =============================================================================
# PART 1: FORENSIC AUDIT & LARGE-SCALE REL606 UPGRADE (2,000 WINDOWS = 300 KBP)
# =============================================================================

SEQ_WAVE_SHIFT_200BP = (
    "TTCGTCGTTTGTTCTTGGAAAGACAGCGCTGGGGATCGTTCGCTGGATACCACCCGCATAACTGGCTCGTC"
    "GCGGTGGGTCAGAGGTCAATGATGAACTTTCAAGTTGACGTGAGAAGCTCTACGGTTGTTGTTCGACTGCT"
    "GTTGCGGCCGTCGTGGCGGGTCACGCGTCATGGGCATTCGTCGTTGGCAGTCCCCACG"
)
SEQ_TUNNELING_200BP = (
    "AGCTGTCGATCGTACGGTACGTTAGGCTAGCTTACGATCGTAGCTGATCGTTAGCTGATCGTAGCATGCTA"
    "GCATCGTAGCTACGATCGTAGCTAGCTGATCGTACGTAGCTAGCTAGCATGCGTAGCTAGCTACGTACGAT"
    "GCGTAGCTACGTACGTAGCATCGTAGCTAGCATCGTAGCTAGCTAGCTAGCTAGCATC"
)


def von_neumann_entropy_bits(rho: np.ndarray) -> float:
    evals = np.linalg.eigvalsh(0.5 * (rho + rho.conj().T))
    evals = evals[evals > 1e-15]
    return float(-np.sum(evals * np.log2(evals)))


def partial_trace_2qubit(rho4: np.ndarray, trace_out: int = 0) -> np.ndarray:
    t = rho4.reshape((2, 2, 2, 2))
    if trace_out == 0:
        return np.trace(t, axis1=0, axis2=2)
    return np.trace(t, axis1=1, axis2=3)


def negativity_peres_horodecki(rho4: np.ndarray) -> float:
    t = rho4.reshape((2, 2, 2, 2))
    rho_pt = np.transpose(t, (0, 3, 2, 1)).reshape((4, 4))
    evals = np.linalg.eigvalsh(0.5 * (rho_pt + rho_pt.conj().T))
    return float(np.sum(np.abs(evals[evals < -1e-15])))


def original_github_density_matrix(seq: str) -> np.ndarray:
    """Exact reproduction of published Von_Neumann_entropy.py (unnormalized C, zero G)."""
    qmap = {
        "A": np.array([1.0, 0.0]),
        "T": np.array([0.0, 1.0]),
        "C": np.array([1.0, 1.0]),   # ||psi||^2 = 2
        "G": np.array([0.0, 0.0]),   # ||psi||^2 = 0
    }
    qs = [qmap[b] for b in seq if b in qmap]
    return sum(np.outer(q, q) for q in qs) / len(qs)


def repaired_bloch_1qubit_density_matrix(seq: str) -> np.ndarray:
    """Trace-preserving 1-qubit Bloch state encoding (Tr rho = 1)."""
    inv_s2 = 1.0 / np.sqrt(2.0)
    qmap = {
        "A": np.array([1.0, 0.0]),
        "T": np.array([0.0, 1.0]),
        "C": np.array([inv_s2, inv_s2]),
        "G": np.array([inv_s2, -inv_s2]),
    }
    qs = [qmap[b] for b in seq if b in qmap]
    return sum(np.outer(q, q) for q in qs) / len(qs)


def repaired_2qubit_ry_ws_dinucleotide_metrics(seq: str) -> dict:
    """Trace-preserving 2-qubit RY x WS chemical basis + adjacent dinucleotide MI."""
    b_idx = {"A": 0, "G": 1, "T": 2, "C": 3}
    idx = np.fromiter((b_idx[b] for b in seq if b in b_idx), dtype=int)
    n = len(idx)
    counts = np.bincount(idx, minlength=4) / n
    rho_1site = np.diag(counts)
    rho_ry = partial_trace_2qubit(rho_1site, trace_out=1)
    rho_ws = partial_trace_2qubit(rho_1site, trace_out=0)
    s_1site = von_neumann_entropy_bits(rho_1site)
    i_ry_ws = von_neumann_entropy_bits(rho_ry) + von_neumann_entropy_bits(rho_ws) - s_1site

    pair_flat = idx[:-1] * 4 + idx[1:]
    dinuc = np.bincount(pair_flat, minlength=16).reshape((4, 4)) / (n - 1)
    p_left = dinuc.sum(axis=1)
    p_right = dinuc.sum(axis=0)
    s_left = -np.sum(p_left[p_left > 0] * np.log2(p_left[p_left > 0]))
    s_right = -np.sum(p_right[p_right > 0] * np.log2(p_right[p_right > 0]))
    s_joint = -np.sum(dinuc[dinuc > 0] * np.log2(dinuc[dinuc > 0]))
    return {
        "shannon_1site_bits": float(s_1site),
        "intra_base_RY_WS_MI_bits": float(i_ry_ws),
        "dinucleotide_MI_bits": float(s_left + s_right - s_joint),
        "gc_fraction": float(counts[1] + counts[3]),
        "g_fraction": float(counts[1]),
    }


def tight_binding_pi_stack_metrics(seq: str) -> dict:
    """Exact 150x150 Hermitian tight-binding pi-stack HOMO hole Hamiltonian
    (Sugiyama-Saito / Voityuk vertical ionization potentials, t_pi = -0.085 eV)."""
    ip_ev = {"G": 7.75, "A": 8.24, "C": 8.87, "T": 9.14}
    e_site = np.array([ip_ev[b] - 7.75 for b in seq], dtype=float)
    n = len(seq)
    h_mat = np.diag(e_site)
    purines = {"A", "G"}
    for i in range(n - 1):
        t_hop = -0.085 if (seq[i] in purines and seq[i + 1] in purines) else -0.045
        h_mat[i, i + 1] = h_mat[i + 1, i] = t_hop
    evals, evecs = np.linalg.eigh(h_mat)
    prob0 = np.abs(evecs[:, 0]) ** 2
    ipr0 = float(np.sum(prob0 ** 2))
    boltz = np.exp(-(evals - evals[0]) / KBT_EV)
    boltz /= np.sum(boltz)
    p_hole_site = (np.abs(evecs) ** 2) @ boltz

    in_g3 = np.zeros(n, dtype=bool)
    i = 0
    while i < n:
        if seq[i] == "G":
            j = i
            while j < n and seq[j] == "G":
                j += 1
            if (j - i) >= 3:
                in_g3[i:j] = True
            i = j
        else:
            i += 1
    return {
        "homo_hole_ground_eV": float(evals[0]),
        "ground_hole_IPR": ipr0,
        "effective_localization_bp": float(1.0 / ipr0),
        "thermal_hole_G3plus_fraction": float(np.sum(p_hole_site[in_g3])),
    }


def audit_and_upgrade_github_scripts() -> tuple[pd.DataFrame, dict, float]:
    t0 = time.perf_counter()
    nc_seq = SEQ_WAVE_SHIFT_200BP[:100]
    c_seq_100 = SEQ_WAVE_SHIFT_200BP[100:]

    rho_nc_orig = original_github_density_matrix(nc_seq)
    rho_c_orig = original_github_density_matrix(c_seq_100)
    tr_nc_orig = float(np.trace(rho_nc_orig))
    tr_c_orig = float(np.trace(rho_c_orig))
    s_nc_orig = von_neumann_entropy_bits(rho_nc_orig)
    s_c_orig = von_neumann_entropy_bits(rho_c_orig)

    rho12_orig = np.kron(rho_nc_orig, rho_c_orig)
    rho12_norm = rho12_orig / np.trace(rho12_orig)
    rho1_norm = partial_trace_2qubit(rho12_norm, trace_out=1)
    rho2_norm = partial_trace_2qubit(rho12_norm, trace_out=0)
    qmi_orig_product = (
        von_neumann_entropy_bits(rho1_norm)
        + von_neumann_entropy_bits(rho2_norm)
        - von_neumann_entropy_bits(rho12_norm)
    )
    neg_orig_product = negativity_peres_horodecki(rho12_norm)

    rho_nc_bloch = repaired_bloch_1qubit_density_matrix(nc_seq)
    rho_c_bloch = repaired_bloch_1qubit_density_matrix(c_seq_100)
    s_nc_bloch = von_neumann_entropy_bits(rho_nc_bloch)
    s_c_bloch = von_neumann_entropy_bits(rho_c_bloch)
    dinuc_nc = repaired_2qubit_ry_ws_dinucleotide_metrics(nc_seq)
    dinuc_c = repaired_2qubit_ry_ws_dinucleotide_metrics(c_seq_100)
    tb_nc = tight_binding_pi_stack_metrics(nc_seq)
    tb_c = tight_binding_pi_stack_metrics(c_seq_100)

    # Large-scale real genome audit across E. coli REL606 (target: 1,000 CDS + 1,000 intergenic 150-bp windows)
    fna_candidates = (
        sorted(glob.glob(f"{RAW_DIR}/**/REL606*.fna", recursive=True))
        + sorted(glob.glob(f"{RAW_DIR}/**/REL606*.fasta", recursive=True))
        + sorted(glob.glob(f"{RAW_DIR}/**/rel606*.fa*", recursive=True))
    )
    gff_candidates = (
        sorted(glob.glob(f"{RAW_DIR}/**/REL606*.gff*", recursive=True))
        + sorted(glob.glob(f"{RAW_DIR}/**/rel606*.gff*", recursive=True))
    )
    gbk_candidates = sorted(glob.glob(f"{RAW_DIR}/**/REL606*.gb*", recursive=True))

    full_genome = ""
    cds_intervals = []
    if fna_candidates and gff_candidates:
        with open(fna_candidates[0], "r") as fh:
            full_genome = "".join(ln.strip().upper() for ln in fh if not ln.startswith(">"))
        with open(gff_candidates[0], "r") as fh:
            for ln in fh:
                if ln.startswith("#"):
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) >= 5 and parts[2] == "CDS":
                    cds_intervals.append((int(parts[3]) - 1, int(parts[4])))
    elif gbk_candidates:
        from Bio import SeqIO
        rec = next(SeqIO.parse(gbk_candidates[0], "genbank"))
        full_genome = str(rec.seq).upper()
        for feat in rec.features:
            if feat.type == "CDS":
                cds_intervals.append((int(feat.location.start), int(feat.location.end)))

    real_window_summary = {}
    if full_genome and cds_intervals:
        cds_intervals.sort()
        # Merge overlapping CDS intervals so intergenic gaps are strictly non-coding!
        merged_cds = []
        for s, e in cds_intervals:
            if not merged_cds or s > merged_cds[-1][1]:
                merged_cds.append([s, e])
            else:
                merged_cds[-1][1] = max(merged_cds[-1][1], e)

        valid_bases = {"A", "C", "G", "T"}
        coding_wins = []
        for s, e in merged_cds:
            pos = s
            while pos + 150 <= e and len(coding_wins) < 1000:
                sub = full_genome[pos:pos + 150]
                if set(sub) <= valid_bases:
                    coding_wins.append(sub)
                pos += 450  # stride across distinct genes to cover the full 4.63 Mb chromosome
            if len(coding_wins) >= 1000:
                break

        igd_wins = []
        for idx_i in range(len(merged_cds) - 1):
            g_start = merged_cds[idx_i][1]
            g_end = merged_cds[idx_i + 1][0]
            pos = g_start
            while pos + 150 <= g_end and len(igd_wins) < 1000:
                sub = full_genome[pos:pos + 150]
                if set(sub) <= valid_bases:
                    igd_wins.append(sub)
                pos += 150  # strictly non-overlapping 150-bp intergenic tiles
            if len(igd_wins) >= 1000:
                break

        n_pair = min(len(coding_wins), len(igd_wins))
        coding_wins = coding_wins[:n_pair]
        igd_wins = igd_wins[:n_pair]

        orig_qmi_zeros = []
        for wc, wnc in zip(coding_wins[:100], igd_wins[:100]):
            r12 = np.kron(original_github_density_matrix(wnc), original_github_density_matrix(wc))
            r12 /= np.trace(r12)
            r1 = partial_trace_2qubit(r12, trace_out=1)
            r2 = partial_trace_2qubit(r12, trace_out=0)
            orig_qmi_zeros.append(
                von_neumann_entropy_bits(r1) + von_neumann_entropy_bits(r2) - von_neumann_entropy_bits(r12)
            )

        orig_s_c = [von_neumann_entropy_bits(original_github_density_matrix(w)) for w in coding_wins]
        orig_s_nc = [von_neumann_entropy_bits(original_github_density_matrix(w)) for w in igd_wins]
        bloch_s_c = [von_neumann_entropy_bits(repaired_bloch_1qubit_density_matrix(w)) for w in coding_wins]
        bloch_s_nc = [von_neumann_entropy_bits(repaired_bloch_1qubit_density_matrix(w)) for w in igd_wins]
        c_dinuc_mi = [repaired_2qubit_ry_ws_dinucleotide_metrics(w)["dinucleotide_MI_bits"] for w in coding_wins]
        nc_dinuc_mi = [repaired_2qubit_ry_ws_dinucleotide_metrics(w)["dinucleotide_MI_bits"] for w in igd_wins]
        c_gc = [repaired_2qubit_ry_ws_dinucleotide_metrics(w)["gc_fraction"] for w in coding_wins]
        nc_gc = [repaired_2qubit_ry_ws_dinucleotide_metrics(w)["gc_fraction"] for w in igd_wins]
        c_g3 = [tight_binding_pi_stack_metrics(w)["thermal_hole_G3plus_fraction"] for w in coding_wins]
        nc_g3 = [tight_binding_pi_stack_metrics(w)["thermal_hole_G3plus_fraction"] for w in igd_wins]

        _, p_orig_s = mannwhitneyu(orig_s_c, orig_s_nc, alternative="two-sided")
        _, p_bloch_s = mannwhitneyu(bloch_s_c, bloch_s_nc, alternative="two-sided")
        _, p_mi = mannwhitneyu(c_dinuc_mi, nc_dinuc_mi, alternative="two-sided")
        _, p_g3 = mannwhitneyu(c_g3, nc_g3, alternative="two-sided")

        real_window_summary = {
            "n_coding_windows": len(coding_wins),
            "n_intergenic_windows": len(igd_wins),
            "total_bp_evaluated": int(150 * (len(coding_wins) + len(igd_wins))),
            "coding_gc_mean": float(np.mean(c_gc)),
            "intergenic_gc_mean": float(np.mean(nc_gc)),
            "orig_product_QMI_max_abs_bits": float(np.max(np.abs(orig_qmi_zeros))),
            "orig_unnorm_S_coding_mean": float(np.mean(orig_s_c)),
            "orig_unnorm_S_coding_sd": float(np.std(orig_s_c, ddof=1)),
            "orig_unnorm_S_intergenic_mean": float(np.mean(orig_s_nc)),
            "orig_unnorm_S_intergenic_sd": float(np.std(orig_s_nc, ddof=1)),
            "mwu_p_orig_unnorm_S": float(p_orig_s),
            "repaired_bloch_S_coding_mean": float(np.mean(bloch_s_c)),
            "repaired_bloch_S_coding_sd": float(np.std(bloch_s_c, ddof=1)),
            "repaired_bloch_S_intergenic_mean": float(np.mean(bloch_s_nc)),
            "repaired_bloch_S_intergenic_sd": float(np.std(bloch_s_nc, ddof=1)),
            "mwu_p_repaired_bloch_S": float(p_bloch_s),
            "coding_dinuc_MI_mean_bits": float(np.mean(c_dinuc_mi)),
            "coding_dinuc_MI_sd_bits": float(np.std(c_dinuc_mi, ddof=1)),
            "intergenic_dinuc_MI_mean_bits": float(np.mean(nc_dinuc_mi)),
            "intergenic_dinuc_MI_sd_bits": float(np.std(nc_dinuc_mi, ddof=1)),
            "mwu_p_dinuc_MI": float(p_mi),
            "coding_G3plus_hole_frac_mean": float(np.mean(c_g3)),
            "coding_G3plus_hole_frac_sd": float(np.std(c_g3, ddof=1)),
            "intergenic_G3plus_hole_frac_mean": float(np.mean(nc_g3)),
            "intergenic_G3plus_hole_frac_sd": float(np.std(nc_g3, ddof=1)),
            "mwu_p_G3plus_hole_frac": float(p_g3),
        }
        pd.DataFrame({
            "window_idx": np.arange(len(coding_wins)),
            "coding_gc": c_gc,
            "intergenic_gc": nc_gc,
            "orig_unnorm_S_coding": orig_s_c,
            "orig_unnorm_S_intergenic": orig_s_nc,
            "repaired_bloch_S_coding": bloch_s_c,
            "repaired_bloch_S_intergenic": bloch_s_nc,
            "coding_dinuc_MI_bits": c_dinuc_mi,
            "intergenic_dinuc_MI_bits": nc_dinuc_mi,
            "coding_G3plus_hole_frac": c_g3,
            "intergenic_G3plus_hole_frac": nc_g3,
        }).to_csv(f"{OUT_DATA}/module5_rel606_coding_vs_intergenic_windows.csv", index=False)

    if not real_window_summary:
        raise FileNotFoundError(f"Provide the REL606 reference FASTA and GFF3 or GenBank under {RAW_DIR}")

    # Tunneling unit-mismatch & Euler wave-shift audit (evaluated purely against generic ambient 1 uV/m dipole)
    lam_si_j = const.e * 1.0e-6 * 1.0e-10                    # 1.6022e-35 J (1 uV/m over 1.0 A)
    lam_ev = 1.0e-16                                         # 1.0e-16 eV
    noise_orig_j = 0.05                                      # 0.05 J in original script
    noise_orig_ev = noise_orig_j / const.e                   # 3.1208e+17 eV
    ratio_noise_to_field = noise_orig_j / lam_si_j           # 3.1208e+33

    p_gc_kms = 1.4601383261440485e-7
    p_at_kms = float(np.exp(-0.56 / KBT_EV))
    tun_nc = SEQ_TUNNELING_200BP[:100]
    tun_c = SEQ_TUNNELING_200BP[100:]
    gc_tun_nc = (tun_nc.count("G") + tun_nc.count("C")) / len(tun_nc)
    gc_tun_c = (tun_c.count("G") + tun_c.count("C")) / len(tun_c)
    p_mut_tun_nc = gc_tun_nc * p_gc_kms + (1.0 - gc_tun_nc) * p_at_kms
    p_mut_tun_c = gc_tun_c * p_gc_kms + (1.0 - gc_tun_c) * p_at_kms

    dt_ws = 1.0e-13  # 100 fs
    euler_courant_factor = 0.08 * dt_ws / HBAR_EV_S          # 12.15 >> 1 (unstable forward Euler)
    # Dimensionless fractional chirp rate |(1/f0)(df/dt)| comparison: Earth diurnal spin = omega_rot^2 * R_E / c
    frac_chirp_true_rot_s_inv = (7.2921159e-5 ** 2) * 6.378137e6 / C_LIGHT  # 1.131e-10 s^-1

    e_hole_200 = np.array([{"G": 0.0, "A": 0.49, "C": 1.12, "T": 1.39}[b] for b in SEQ_WAVE_SHIFT_200BP])
    h_seq = np.diag(e_hole_200)
    for i in range(199):
        h_seq[i, i + 1] = h_seq[i + 1, i] = -0.055
    pos_m = (np.arange(200) - 99.5) * 3.4e-10
    v_dip_ev = 1.0e-6 * pos_m
    u_0 = expm(-1j * h_seq * (1.0e-12 / HBAR_EV_S))
    u_ext = expm(-1j * (h_seq + np.diag(v_dip_ev)) * (1.0e-12 / HBAR_EV_S))
    propagator_diff_1ps = float(np.linalg.norm(u_ext - u_0, ord=2))

    with open(f"{OUT_DATA}/module1_summary.json", "r") as fh:
        m1_ref = json.load(fh)
    ac_ref_val = float(m1_ref["ambient_subgap_1uVm_reference"]["ac_subgap_rectified_dmu_over_mu"])
    tau_mem_ref = float(m1_ref["tau_mem_fs"])

    rows = [
        {
            "script_family": "1_Von_Neumann_Entropy",
            "original_artifact": f"Unnormalized C=[1,1] (norm^2=2) & erased G=[0,0]; Tr(rho_NC)={tr_nc_orig:.3f}, Tr(rho_C)={tr_c_orig:.3f}",
            "original_reported_metric": f"S_orig(NC)={s_nc_orig:.4f}, S_orig(C)={s_c_orig:.4f} (GC-skew trace artifact)",
            "repaired_physical_model": "Trace-preserving Bloch 1-qubit (Tr=1) + 2-qubit RY x WS dinucleotide MI",
            "repaired_metric_value": f"S_Bloch(NC)={s_nc_bloch:.4f}, S_Bloch(C)={s_c_bloch:.4f} bits; Dinuc MI: NC={dinuc_nc['dinucleotide_MI_bits']:.4f}, C={dinuc_c['dinucleotide_MI_bits']:.4f} bits",
        },
        {
            "script_family": "2_Coding_NonCoding_Entanglement",
            "original_artifact": "Outer-product sum builds separable product state rho_12 = rho_1 (x) rho_2",
            "original_reported_metric": f"QMI I(NC:C)={qmi_orig_product:.2e} bits; Negativity N={neg_orig_product:.2e} (exact zero entanglement)",
            "repaired_physical_model": "150x150 tight-binding pi-stack HOMO hole Hamiltonian (Sugiyama-Saito IPs, t_pi=-0.085 eV)",
            "repaired_metric_value": f"Thermal G3+ hole capture: NC={tb_nc['thermal_hole_G3plus_fraction']:.3f} vs C={tb_c['thermal_hole_G3plus_fraction']:.3f}; IPR^-1={tb_nc['effective_localization_bp']:.1f} vs {tb_c['effective_localization_bp']:.1f} bp",
        },
        {
            "script_family": "3_Proton_Tunneling_Stochastic_Resonance",
            "original_artifact": f"SI/dimensionless mismatch: noise=0.05 J ({noise_orig_ev:.2e} eV) vs 1 uV/m dipole ({lam_si_j:.2e} J, {ratio_noise_to_field:.2e}x); symmetric well dE=0",
            "original_reported_metric": "Unphysical ~0.50 right-well saturation driven by 3.12e17 eV noise kicks",
            "repaired_physical_model": "Asymmetric GC (dE=0.42 eV) & AT (dE=0.56 eV) double well with KMS Lindblad bath at 310.15 K",
            "repaired_metric_value": f"p_mut(NC)={p_mut_tun_nc:.3e}, p_mut(C)={p_mut_tun_c:.3e}; sub-thermal AC 1 uV/m |dmu/mu|={ac_ref_val:.2e}",
        },
        {
            "script_family": "4_External_Doppler_Wave_Shift",
            "original_artifact": f"Unstable forward Euler (||H||dt/hbar={euler_courant_factor:.1f}>>1), random non-sequence H, NC-only drive, unphysical chirp",
            "original_reported_metric": "Artificial FFT peak in non-coding block only from asymmetric drive injection",
            "repaired_physical_model": f"Exact matrix-exponential + KMS dephasing (tau_mem={tau_mem_ref:.1f} fs) + Earth-Sun Keplerian/diurnal kinematics",
            "repaired_metric_value": f"||U_ext - U_0||_2 at 1 ps = {propagator_diff_1ps:.3e}; fractional rotational chirp |(1/f)(df/dt)| = {frac_chirp_true_rot_s_inv:.2e} s^-1",
        },
    ]

    elapsed_s = time.perf_counter() - t0
    audit_dict = {
        "trace_orig_nc": tr_nc_orig,
        "trace_orig_c": tr_c_orig,
        "product_state_qmi_bits": qmi_orig_product,
        "product_state_negativity": neg_orig_product,
        "bloch_entropy_nc_bits": s_nc_bloch,
        "bloch_entropy_c_bits": s_c_bloch,
        "dinuc_nc": dinuc_nc,
        "dinuc_c": dinuc_c,
        "tight_binding_nc": tb_nc,
        "tight_binding_c": tb_c,
        "real_rel606_windows": real_window_summary,
        "tunneling_lambda_1uVm_J": lam_si_j,
        "tunneling_lambda_1uVm_eV": lam_ev,
        "tunneling_orig_noise_eV": noise_orig_ev,
        "ratio_orig_noise_to_field": ratio_noise_to_field,
        "euler_instability_courant": euler_courant_factor,
        "fractional_diurnal_chirp_s_inv": frac_chirp_true_rot_s_inv,
        "propagator_norm_diff_1ps_1uVm": propagator_diff_1ps,
    }
    return pd.DataFrame(rows), audit_dict, elapsed_s


# =============================================================================
# PART 2: 8,766-HOUR EARTH-SUN EPHEMERIS, DIRECT CHIRAL SOLAR NEUTRINOS (nu_e),
#         AND SOLAR-MODULATED SPIN-POLARIZED MUON TIME-ENCODING
# =============================================================================

def solve_kepler(mean_anom: np.ndarray, ecc: float, tol: float = 1e-13) -> np.ndarray:
    E = mean_anom.copy()
    for _ in range(14):
        dE = (E - ecc * np.sin(E) - mean_anom) / (1.0 - ecc * np.cos(E))
        E -= dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def solar_neutrino_electroweak_cross_section_cm2(e_nu_mev: float, p_ee: float) -> tuple[float, float]:
    """Exact Standard Model tree-level Weinberg-Salam nu-e elastic scattering cross section
    averaged over MSW solar survival probability P_ee (nu_e CC+NC vs nu_mu/tau NC),
    plus maximum recoil electron kinetic energy T_e_max (MeV)."""
    m_e_mev = 0.51099895
    t_max_mev = (2.0 * (e_nu_mev ** 2)) / (m_e_mev + 2.0 * e_nu_mev)
    z_max = t_max_mev / e_nu_mev

    # sigma_0 = 2 * G_F^2 * m_e / pi in natural units -> convert GeV^-2 * MeV to cm^2
    # (hbar * c) = 1.9732698e-11 MeV*cm; G_F = 1.1663787e-11 MeV^-2
    gf_mev2 = GF_GEV2 * 1.0e-6
    hbar_c_mev_cm = HBAR_C_MEV_FM * 1.0e-13
    prefac_cm2_per_mev = (2.0 * (gf_mev2 ** 2) * m_e_mev / np.pi) * (hbar_c_mev_cm ** 2)
    # Standard numeric check: prefac_cm2_per_mev = 1.723e-44 cm^2 / MeV

    def sigma_flavor(gv: float, ga: float) -> float:
        g_plus = 0.5 * (gv + ga)   # g_L
        g_minus = 0.5 * (gv - ga)  # g_R
        integral = (
            (g_plus ** 2) * z_max
            + (g_minus ** 2) * (1.0 - (1.0 - z_max) ** 3) / 3.0
            - g_plus * g_minus * (m_e_mev / e_nu_mev) * 0.5 * (z_max ** 2)
        )
        return float(prefac_cm2_per_mev * e_nu_mev * integral)

    s_nue = sigma_flavor(0.5 + 2.0 * SIN2_THETA_W, 0.5)
    s_numu = sigma_flavor(-0.5 + 2.0 * SIN2_THETA_W, -0.5)
    return float(p_ee * s_nue + (1.0 - p_ee) * s_numu), float(t_max_mev)


def compute_solar_neutrino_and_muon_astrodynamics() -> tuple[pd.DataFrame, pd.DataFrame, dict, dict, float]:
    """Evaluate 8,766 hourly J2000 Earth-Sun orbital + diurnal rotation kinematics for:
    (1) Direct chiral solar neutrinos (pp, 7Be 861.8 keV, pep 1442.2 keV, 8B) and
    (2) Solar-modulated spin-polarized atmospheric/solar cosmic-ray muons."""
    t0 = time.perf_counter()

    a_au_m = 1.495978707e11
    ecc = 0.0167086
    p_orb_s = 365.256363 * 86400.0
    n_mean = 2.0 * np.pi / p_orb_s
    mu_sun = (n_mean ** 2) * (a_au_m ** 3)
    obliquity = np.deg2rad(23.43929)
    varpi = np.deg2rad(102.937)
    omega_rot = 7.2921159e-5
    r_earth = 6.378137e6
    v_rot_eq = omega_rot * r_earth  # 465.10 m/s

    # Full 8,766-hour hourly ephemeris across 365.25 solar days
    hours = np.arange(8766, dtype=float)
    days_hr = hours / 24.0
    t_sec = (days_hr - 3.0) * 86400.0
    M_anom = n_mean * t_sec
    E_anom = solve_kepler(M_anom, ecc)
    nu_true = 2.0 * np.arctan2(
        np.sqrt(1.0 + ecc) * np.sin(0.5 * E_anom),
        np.sqrt(1.0 - ecc) * np.cos(0.5 * E_anom),
    )
    r_m = a_au_m * (1.0 - ecc * np.cos(E_anom))
    v_rad_sun_m_s = (n_mean * a_au_m * ecc / np.sqrt(1.0 - ecc ** 2)) * np.sin(nu_true)
    v_orb_m_s = np.sqrt(mu_sun * (2.0 / r_m - 1.0 / a_au_m))

    # Solar ecliptic longitude & equatorial declination delta_sun(t)
    lon_sun = varpi + nu_true + np.pi
    sin_dec = np.sin(obliquity) * np.sin(lon_sun)
    cos_dec = np.sqrt(1.0 - sin_dec ** 2)

    # Local hour angle H_solar(t) = 2*pi*(hour_of_day - 12)/24; radial rotational velocity toward Sun:
    # At sunrise (06:00, H = -pi/2), observer rotates TOWARD Sun -> blueshift (+v_rot * cos(lat) * cos(dec))
    # At sunset (18:00, H = +pi/2), observer rotates AWAY from Sun -> redshift (-v_rot * cos(lat) * cos(dec))
    h_angle = 2.0 * np.pi * ((hours % 24.0) - 12.0) / 24.0
    lat_mid = np.deg2rad(41.0)  # Gran Sasso / mid-latitude reference (42.45N / 41N)
    v_los_rot_eq = -v_rot_eq * cos_dec * np.sin(h_angle)  # negative v_los = approaching Sun -> blueshift!
    v_los_rot_mid = -v_rot_eq * np.cos(lat_mid) * cos_dec * np.sin(h_angle)

    # Total line-of-sight recession velocity from Sun: v_rec(t) = v_rad_sun(t) + v_los_rot(t)
    # Relativistic Doppler shift of solar neutrinos: dE / E0 = - v_rec / c + O(v^2/c^2)
    e_be7_ev = 861.8e3
    e_pep_ev = 1442.2e3
    flux_factor = (a_au_m / r_m) ** 2
    flux_mod_pct = 100.0 * (flux_factor - 1.0)

    be7_doppler_orb_ev = -e_be7_ev * (v_rad_sun_m_s / C_LIGHT)
    be7_doppler_tot_eq_ev = -e_be7_ev * ((v_rad_sun_m_s + v_los_rot_eq) / C_LIGHT)
    pep_doppler_orb_ev = -e_pep_ev * (v_rad_sun_m_s / C_LIGHT)

    # Verify 2D quadrature bijection: normalized phase radius R_quad^2 = ((flux_factor-1)/(2e))^2 + (v_rad / v_rad_max)^2 ~ 1!
    v_rad_max = float(np.max(np.abs(v_rad_sun_m_s)))
    quadrature_radius = np.sqrt(((flux_factor - 1.0) / (2.0 * ecc)) ** 2 + (v_rad_sun_m_s / v_rad_max) ** 2)

    # Solar neutrino Standard Model electroweak audit on DNA & whole human cell
    # Standard Solar Model (B16-GS98) fluxes at 1 AU (cm^-2 s^-1) & mean/line energies (MeV):
    nu_branches = [
        ("pp_continuum", 5.98e10, 0.267, 0.55),
        ("7Be_862keV_line", 4.34e9, 0.8618, 0.52),
        ("7Be_384keV_line", 0.50e9, 0.3843, 0.54),
        ("pep_1442keV_line", 1.44e8, 1.4422, 0.50),
        ("CNO_continuum", 4.88e8, 0.740, 0.52),
        ("8B_continuum", 5.16e6, 6.700, 0.36),
    ]
    # Diploid human DNA: 6.4e9 bp, average bp mass ~ 618 Da -> ~ 316 electrons per bp -> N_e_DNA = 2.022e12 electrons
    n_e_diploid_dna = 6.4e9 * 316.0
    # Whole 1-ng mammalian cell (~3.34e14 electrons) & Borexino 100-ton fiducial target (~3.30e31 electrons)
    n_e_whole_cell = 3.34e14
    n_e_borexino_100t = 3.30e31
    sec_per_yr = 365.25 * 86400.0

    branch_rates = []
    tot_rate_dna_s = 0.0
    tot_rate_cell_s = 0.0
    for b_name, flux_cm2_s, e_mev, p_ee in nu_branches:
        sig_cm2, tmax_mev = solar_neutrino_electroweak_cross_section_cm2(e_mev, p_ee)
        r_dna_s = flux_cm2_s * sig_cm2 * n_e_diploid_dna
        r_cell_s = flux_cm2_s * sig_cm2 * n_e_whole_cell
        tot_rate_dna_s += r_dna_s
        tot_rate_cell_s += r_cell_s
        branch_rates.append({
            "branch": b_name,
            "flux_cm2_s": flux_cm2_s,
            "mean_or_line_energy_MeV": e_mev,
            "msw_P_ee": p_ee,
            "sigma_nu_e_cm2": sig_cm2,
            "max_electron_recoil_keV": tmax_mev * 1000.0,
            "rate_per_diploid_DNA_yr": r_dna_s * sec_per_yr,
            "rate_per_cell_yr": r_cell_s * sec_per_yr,
        })

    # Coherent elastic neutrino-nucleus scattering (CEvNS) check on DNA heavy nuclei (C, N, O, P; ~40 heavy nuclei/bp = 2.56e11)
    # Primarily 8B (mean E = 6.7 MeV, mean N_neutron ~ 7): sigma_CEvNS = G_F^2 * E^2 * N^2 / (4*pi)
    gf_mev2 = GF_GEV2 * 1.0e-6
    hbar_c_mev_cm = HBAR_C_MEV_FM * 1.0e-13
    sig_cevns_8b_cm2 = ((gf_mev2 ** 2) * (6.70 ** 2) * (7.0 ** 2) / (4.0 * np.pi)) * (hbar_c_mev_cm ** 2)
    rate_cevns_8b_dna_yr = 5.16e6 * sig_cevns_8b_cm2 * (6.4e9 * 40.0) * sec_per_yr
    # Maximum nuclear recoil of a 861.8 keV neutrino on 12C (M = 11.178 GeV).
    t_nuc_max_be7_12c_ev = (2.0 * (861.8e3 ** 2)) / (11.178e9 + 2.0 * 861.8e3)

    # Stodolsky coherent weak neutral-current spin/chiral energy splitting: dE_WNC = sqrt(2) * G_F * (hbar*c)^3 * n_nu
    n_nu_cm3 = 6.46e10 / (C_LIGHT * 100.0)  # ~ 2.155 cm^-3
    de_stodolsky_ev = np.sqrt(2.0) * (GF_GEV2 * 1.0e-18) * ((HBAR_C_MEV_FM * 1.0e6 * 1.0e-13) ** 3) * n_nu_cm3

    # Minimum cell ensemble nu* needed to detect +/-3.34% Keplerian neutrino flux modulation at SNR = 5 in 1 year:
    # For sinusoidal amplitude A = 2*ecc = 0.03342, Poisson Fisher information over N_tot total events is I_F = N_tot / 2;
    # SNR = A * sqrt(N_tot / 2) >= 5 => N_tot >= 2 * (5 / A)^2 = 44,775 events!
    n_events_5sigma_orbital = 2.0 * ((5.0 / (2.0 * ecc)) ** 2)
    cells_needed_1yr_5sigma = n_events_5sigma_orbital / (tot_rate_cell_s * sec_per_yr)

    # Solar-modulated spin-polarized cosmic-ray / SEP muon astrodynamic clocks across the 8,766-hour year:
    # 1. Orbital Compton-Getting relativistic Doppler dipole: (gamma + 2)*(v_orb / c)*cos(lon - apex)
    gamma_cr = 2.70
    cg_orbital_pct = 100.0 * (gamma_cr + 2.0) * (v_orb_m_s / C_LIGHT) * np.cos(lon_sun)
    # 2. Annual stratospheric temperature pion-decay modulation (Gaisser et al. 2019; Malgin 2015):
    # Sea-level weighted muon seasonal wave ~ +/-1.8% (peaking early July, day ~ 185 in N hemisphere; +/-8% for high-E TeV muons)
    seasonal_strat_pct = 1.80 * np.cos(2.0 * np.pi * (days_hr - 185.0) / 365.25)
    # 3. 27.2753-day synodic Carrington solar rotation Parker-spiral modulation (~ +/-1.5%)
    carrington_27d_pct = 1.50 * np.cos(2.0 * np.pi * days_hr / 27.2753)
    # 4. Diurnal solar-time muon modulation (~ +/-0.15%, 0.30% peak-to-trough dip centered at 11:00 local solar time; Wang et al. 2025)
    diurnal_muon_pct = -0.15 * np.cos(2.0 * np.pi * ((hours % 24.0) - 11.0) / 24.0)

    combined_muon_pct = seasonal_strat_pct + carrington_27d_pct + cg_orbital_pct + diurnal_muon_pct

    # Save daily ephemeris table (366 daily midpoints) for Figure 5A
    daily_idx = np.arange(0, 8766, 24)
    df_astro_daily = pd.DataFrame({
        "day_of_year": days_hr[daily_idx],
        "heliocentric_r_AU": r_m[daily_idx] / a_au_m,
        "v_orbital_km_s": v_orb_m_s[daily_idx] / 1.0e3,
        "v_radial_sun_m_s": v_rad_sun_m_s[daily_idx],
        "solar_nu_flux_mod_pct": flux_mod_pct[daily_idx],
        "be7_nu_orbital_doppler_eV": be7_doppler_orb_ev[daily_idx],
        "pep_nu_orbital_doppler_eV": pep_doppler_orb_ev[daily_idx],
        "seasonal_stratospheric_muon_pct": seasonal_strat_pct[daily_idx],
        "carrington_27d_muon_pct": carrington_27d_pct[daily_idx],
        "compton_getting_muon_pct": cg_orbital_pct[daily_idx],
        "combined_muon_clock_pct": combined_muon_pct[daily_idx],
    })

    # Save 24-hour equinox diurnal clock table
    eq_start = 80 * 24  # ~ March 21 equinox
    sl = slice(eq_start, eq_start + 24)
    hr_24 = np.arange(24, dtype=float)
    # Solar zenith cos(theta_z) at equator on equinox (daytime >= 0):
    cos_zenith_sun = np.maximum(0.0, np.cos(2.0 * np.pi * (hr_24 - 12.0) / 24.0))
    df_diurnal = pd.DataFrame({
        "local_solar_hour": hr_24,
        "be7_diurnal_doppler_equator_eV": -e_be7_ev * (v_los_rot_eq[sl] / C_LIGHT),
        "be7_diurnal_doppler_midlat41N_eV": -e_be7_ev * (v_los_rot_mid[sl] / C_LIGHT),
        "surface_muon_diurnal_pct": diurnal_muon_pct[sl],
        "solar_SEP_zenith_cos2_projection": cos_zenith_sun ** 2,
    })

    nu_summary = {
        "ephemeris_total_hours": int(len(hours)),
        "total_solar_nu_flux_1AU_cm2_s": float(sum(b[1] for b in nu_branches)),
        "keplerian_flux_half_amp_pct": float(0.5 * (np.max(flux_mod_pct) - np.min(flux_mod_pct))),
        "keplerian_flux_peak_to_trough_pct": float(np.max(flux_mod_pct) - np.min(flux_mod_pct)),
        "borexino_measured_2ecc_pct": "3.68 +/- 0.65% (>5 sigma; Appel et al. 2022)",
        "be7_862keV_orbital_doppler_amp_eV": float(np.max(np.abs(be7_doppler_orb_ev))),
        "be7_862keV_diurnal_equator_amp_eV": float(e_be7_ev * v_rot_eq / C_LIGHT),
        "be7_862keV_diurnal_midlat41N_amp_eV": float(e_be7_ev * v_rot_eq * np.cos(lat_mid) / C_LIGHT),
        "be7_total_max_doppler_swing_eV": float(np.max(np.abs(be7_doppler_tot_eq_ev))),
        "pep_1442keV_orbital_doppler_amp_eV": float(np.max(np.abs(pep_doppler_orb_ev))),
        "quadrature_phase_ellipse_min_max_norm": [float(np.min(quadrature_radius)), float(np.max(quadrature_radius))],
        "neutrino_chirality_helicity_7Be": float(-1.0 + 0.5 * ((0.05 / 861.8e3) ** 2)),
        "diploid_DNA_nu_e_scatter_rate_per_yr": float(tot_rate_dna_s * sec_per_yr),
        "whole_cell_nu_e_scatter_rate_per_yr": float(tot_rate_cell_s * sec_per_yr),
        "human_body_37T_cells_nu_e_events_per_yr": float(tot_rate_cell_s * sec_per_yr * 3.7e13),
        "diploid_DNA_CEvNS_8B_rate_per_yr": float(rate_cevns_8b_dna_yr),
        "be7_max_12C_nuclear_recoil_eV": float(t_nuc_max_be7_12c_ev),
        "stodolsky_coherent_WNC_splitting_eV": float(de_stodolsky_ev),
        "whole_cells_required_5sigma_orbital_clock_1yr": float(cells_needed_1yr_5sigma),
        "diploid_DNA_targets_required_5sigma_orbital_clock_1yr": float(
            n_events_5sigma_orbital / (tot_rate_dna_s * sec_per_yr)
        ),
        "branches": branch_rates,
    }

    muon_astro_summary = {
        "v_orb_min_km_s": float(np.min(v_orb_m_s) / 1.0e3),
        "v_orb_max_km_s": float(np.max(v_orb_m_s) / 1.0e3),
        "v_radial_sun_max_m_s": v_rad_max,
        "v_rot_equator_m_s": float(v_rot_eq),
        "compton_getting_orbital_dipole_frac": float(np.max(np.abs(cg_orbital_pct)) / 100.0),
        "compton_getting_diurnal_equator_frac": float((gamma_cr + 2.0) * (v_rot_eq / C_LIGHT)),
        "seasonal_stratospheric_surface_muon_amp_pct": 1.80,
        "seasonal_stratospheric_highE_IceCube_amp_pct": 8.00,
        "carrington_27d_solar_rotation_amp_pct": 1.50,
        "diurnal_solar_time_muon_dip_pct": 0.30,
        "solar_cycle_11yr_muon_amp_pct": 4.50,
        "young_archean_sun_SEP_enhancement_factor_range": [100.0, 10000.0],
        "electromagnetic_to_weak_cross_section_ratio": float(1.0e-17 / branch_rates[1]["sigma_nu_e_cm2"]),
    }

    elapsed_s = time.perf_counter() - t0
    return df_astro_daily, df_diurnal, nu_summary, muon_astro_summary, elapsed_s


# =============================================================================
# PART 3: EXPLICIT 240,000-TRAJECTORY 3D MONTE CARLO HELICAL DNA SCATTERING &
#         GLOBUS-BLANDFORD CHIRAL HOLE FUNNELING
# =============================================================================

def pion_decay_muon_lodacity() -> dict:
    """Exact Globus & Blandford (2020, Appendix A.2, Eqs. A1-A5) pion -> muon + nu_mu
    weak-decay chirality transfer kinematics."""
    m_pi = 139.57039
    m_mu = 105.6583755
    r_pi = (m_mu / m_pi) ** 2
    beta_star = (1.0 - r_pi) / (1.0 + r_pi)

    cos_th_star = np.linspace(-0.9999, 0.9999, 4000)
    z = 0.5 * (1.0 + r_pi) + 0.5 * (1.0 - r_pi) * cos_th_star
    p_long_mu_plus = -((1.0 + r_pi) * z - 2.0 * r_pi) / ((1.0 - r_pi) * z)
    p_spec = 2.7
    weights = z ** (p_spec - 1.0)
    mean_lodacity_mu_plus = float(np.sum(p_long_mu_plus * weights) / np.sum(weights))
    r_mu_ratio = 1.30
    net_charge_excess_mu = (r_mu_ratio - 1.0) / (r_mu_ratio + 1.0)
    return {
        "r_pi": float(r_pi),
        "beta_mu_rest": float(beta_star),
        "spectrum_weighted_lodacity_mu_plus": mean_lodacity_mu_plus,
        "muon_charge_ratio_R": r_mu_ratio,
        "muon_charge_asymmetry_factor": float(net_charge_excess_mu),
        "effective_net_electrostatic_lodacity_mu": float(abs(mean_lodacity_mu_plus) * net_charge_excess_mu),
    }


def run_3d_monte_carlo_helical_scattering(n_per_ensemble: int = 60000) -> tuple[pd.DataFrame, dict, float]:
    """Explicit 3D Monte Carlo trajectory sampling recommended by Globus & Blandford
    (2020, Appendix B.4, p. 16: 'performing a Monte Carlo sampling of isotropic cosmic ray
    trajectories') across N_tot = 4 * 60,000 = 240,000 3D helical B-DNA trajectories
    (2 cylinder sides x 4 rotational quadrature phases x 40 steps = 7.68e7 3D force evaluations).

    Evaluates:
      1. Isotropic spin-polarized e+/mu+ (q=+e, L=-0.5) on right-handed D-DNA (m = +2)
      2. Isotropic spin-polarized e+/mu+ (q=+e, L=-0.5) on mirror L-DNA (m = -2)
      3. Helicity-reversed control e+/mu+ (q=+e, L=+0.5) on right-handed D-DNA (m = +2)
      4. Flat mineral-adsorbed horizontal duplex under cos^2(theta_z) zenith shower on D-DNA (m = +2)
    """
    t0 = time.perf_counter()
    rng = np.random.default_rng(20260419)

    pitch_m = 3.4e-9
    k_helix = 2.0 * np.pi / pitch_m
    r_eval = 0.55e-9
    R1_amp = 0.80

    t_1kev = 1000.0
    m_e_ev = 510998.95
    beta_1kev = np.sqrt(1.0 - 1.0 / ((1.0 + t_1kev / m_e_ev) ** 2))
    v0 = beta_1kev * C_LIGHT
    c_geom = 1.0e-7 / ((ALPHA_FS ** 5) * ((1.0 / beta_1kev) ** 3))
    mu_bohr = const.e * const.hbar / (2.0 * const.m_e)

    n_steps = 40
    s_grid = np.linspace(-1.5 * pitch_m, 1.5 * pitch_m, n_steps)  # (40,)
    ds = s_grid[1] - s_grid[0]
    dt = ds / v0
    # 4-point cylindrical rotational phase quadrature (phi_0 = 0, pi/2, pi, 3*pi/2) per 3D ray
    # eliminates non-chiral O(R1) dipole scatter identically while preserving Monte Carlo sampling
    # over 3D cosmic-ray directions (cos_psi, +/-sin_psi) and entry phase offsets.
    phase_shifts = np.array([0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi])[None, :, None]  # (1, 4, 1)

    def sample_ensemble_dr(m_chiral: int, q_sign: float, lodacity: float, oriented_zenith: bool) -> np.ndarray:
        q_c = q_sign * const.e
        dR1_dr = -abs(m_chiral) * R1_amp / r_eval  # R1(r) ~ r^-|m| outside nucleobase core
        if not oriented_zenith:
            cos_psi = rng.uniform(-1.0, 1.0, size=n_per_ensemble)
        else:
            # Atmospheric cos^2(theta_z) zenith shower incident on flat-lying mineral-adsorbed horizontal duplexes:
            # cos(psi) = sin(theta_z) * cos(phi_h), placing trajectories predominantly outside psi_crit = 26.94 deg.
            cos_tz = rng.uniform(0.0, 1.0, size=n_per_ensemble) ** (1.0 / 3.0)
            sin_tz = np.sqrt(np.maximum(0.0, 1.0 - cos_tz ** 2))
            phi_h = rng.uniform(0.0, 2.0 * np.pi, size=n_per_ensemble)
            cos_psi = sin_tz * np.cos(phi_h)

        sin_mag = np.sqrt(np.maximum(0.0, 1.0 - cos_psi ** 2))
        theta0_base = rng.uniform(0.0, 2.0 * np.pi, size=n_per_ensemble)

        dr_ray_out = np.empty(n_per_ensemble, dtype=float)
        chunk = 10000
        s_3d = s_grid[None, None, :]  # (1, 1, 40)
        for i0 in range(0, n_per_ensemble, chunk):
            i1 = min(i0 + chunk, n_per_ensemble)
            cp = cos_psi[i0:i1, None, None]
            sm = sin_mag[i0:i1, None, None]
            th_b = theta0_base[i0:i1, None, None] + phase_shifts  # (C, 4, 1)

            dr_sides = 0.0
            for az_sgn in (+1.0, -1.0):
                vz = v0 * cp
                vphi = v0 * (az_sgn * sm)
                th_s = th_b + (k_helix * vz - (m_chiral / r_eval) * vphi) * (s_3d / v0)
                cos_th = np.cos(th_s)
                sin_th = np.sin(th_s)

                e_r_1 = -dR1_dr * cos_th
                e_phi_1 = -(m_chiral / r_eval) * R1_amp * sin_th
                e_z_1 = k_helix * R1_amp * sin_th

                dv_r_1 = (q_c / const.m_e) * np.cumsum(e_r_1, axis=2) * dt
                dv_phi_1 = (q_c / const.m_e) * np.cumsum(e_phi_1, axis=2) * dt
                dv_z_1 = (q_c / const.m_e) * np.cumsum(e_z_1, axis=2) * dt

                # Exact Dirac spin-orbit radial force F_r = -d_r U_SO with Thomas factor 1/2:
                # U_SO = +0.5 * sgn(q) * (mu_B * L / c^2) * [v0_hat . (dv_E1 x E_1)]
                # term_z ~ r^(-2|m|-2), term_phi ~ r^(-2|m|-1)
                term_z = (vz / v0) * (dv_r_1 * e_phi_1 - dv_phi_1 * e_r_1)
                term_phi = (vphi / v0) * (dv_z_1 * e_r_1 - dv_r_1 * e_z_1)
                f_chiral_r = 0.5 * q_sign * (mu_bohr * lodacity / (C_LIGHT ** 2)) * (
                    ((2.0 * abs(m_chiral) + 2.0) / r_eval) * term_z
                    + ((2.0 * abs(m_chiral) + 1.0) / r_eval) * term_phi
                )
                vr_pert = np.cumsum(f_chiral_r / const.m_e, axis=2) * dt
                dr_val = np.sum(vr_pert, axis=2) * dt  # (C, 4)
                dr_sides = dr_sides + 0.5 * np.mean(dr_val, axis=1)
            dr_ray_out[i0:i1] = dr_sides
        return dr_ray_out

    ensembles = [
        ("D_DNA_m_plus2_isotropic_pos", +2, +1.0, -0.5, False),
        ("L_DNA_m_minus2_isotropic_pos", -2, +1.0, -0.5, False),
        ("D_DNA_m_plus2_helicity_rev_ctrl", +2, +1.0, +0.5, False),
        ("D_DNA_m_plus2_zenith_oriented", +2, +1.0, -0.5, True),
    ]

    checkpoints = [500, 1500, 5000, 15000, 30000, n_per_ensemble]
    conv_rows = []
    final_stats = {}

    for name, m_val, q_val, l_val, is_oriented in ensembles:
        dr_arr = sample_ensemble_dr(m_val, q_val, l_val, is_oriented)
        for n_ck in checkpoints:
            sub = dr_arr[:n_ck]
            mean_v = float(np.mean(sub))
            sem_v = float(np.std(sub, ddof=1) / np.sqrt(n_ck))
            inward_frac = float(np.mean(sub < 0.0))
            conv_rows.append({
                "ensemble": name,
                "n_trajectories": int(n_ck),
                "mean_dr_m": mean_v,
                "sem_dr_m": sem_v,
                "z_score_vs_zero": float(mean_v / max(sem_v, 1e-50)),
                "inward_trajectory_fraction": inward_frac,
            })
        final_stats[name] = {
            "mean_dr_m": float(np.mean(dr_arr)),
            "sem_dr_m": float(np.std(dr_arr, ddof=1) / np.sqrt(len(dr_arr))),
            "z_score": float(np.mean(dr_arr) / (np.std(dr_arr, ddof=1) / np.sqrt(len(dr_arr)))),
            "inward_fraction": float(np.mean(dr_arr < 0.0)),
        }

    # Exact deterministic quadrature reference check across cos_psi in [-1, 1]
    cos_det = np.linspace(-0.98, 0.98, 64)
    sin_det = np.sqrt(1.0 - cos_det ** 2)

    def exact_quadrature_dr(m_chiral: int) -> float:
        q_c = const.e
        dR1_dr = -abs(m_chiral) * R1_amp / r_eval
        cp = cos_det[:, None, None]
        sm = sin_det[:, None, None]
        s_3d = s_grid[None, None, :]
        tot = 0.0
        for az_sgn in (+1.0, -1.0):
            vz = v0 * cp
            vphi = v0 * (az_sgn * sm)
            th_s = phase_shifts + (k_helix * vz - (m_chiral / r_eval) * vphi) * (s_3d / v0)
            e_r_1 = -dR1_dr * np.cos(th_s)
            e_phi_1 = -(m_chiral / r_eval) * R1_amp * np.sin(th_s)
            e_z_1 = k_helix * R1_amp * np.sin(th_s)
            dv_r_1 = (q_c / const.m_e) * np.cumsum(e_r_1, axis=2) * dt
            dv_phi_1 = (q_c / const.m_e) * np.cumsum(e_phi_1, axis=2) * dt
            dv_z_1 = (q_c / const.m_e) * np.cumsum(e_z_1, axis=2) * dt
            term_z = (vz / v0) * (dv_r_1 * e_phi_1 - dv_phi_1 * e_r_1)
            term_phi = (vphi / v0) * (dv_z_1 * e_r_1 - dv_r_1 * e_z_1)
            f_r = 0.5 * (+1.0) * (mu_bohr * (-0.5) / (C_LIGHT ** 2)) * (
                ((2.0 * abs(m_chiral) + 2.0) / r_eval) * term_z
                + ((2.0 * abs(m_chiral) + 1.0) / r_eval) * term_phi
            )
            vr_p = np.cumsum(f_r / const.m_e, axis=2) * dt
            tot += 0.5 * float(np.mean(np.sum(vr_p, axis=2) * dt))
        return tot

    dr_quad_d = exact_quadrature_dr(+2)
    dr_quad_l = exact_quadrature_dr(-2)

    elapsed_s = time.perf_counter() - t0
    mc_summary = {
        "total_3d_trajectories_integrated": int(4 * n_per_ensemble),
        "total_3d_force_evaluations": int(4 * n_per_ensemble * n_steps * 2 * 4),
        "psi_crit_deg": float(np.degrees(np.arctan(k_helix * r_eval / 2.0))),
        "c_geom_calibrated": float(c_geom),
        "exact_quadrature_D_DNA_m": dr_quad_d,
        "exact_quadrature_L_DNA_m": dr_quad_l,
        "exact_quadrature_antisymmetry_sum_m": float(dr_quad_d + dr_quad_l),
        "monte_carlo_ensembles": final_stats,
        "oriented_to_isotropic_amplification_ratio": float(
            final_stats["D_DNA_m_plus2_zenith_oriented"]["mean_dr_m"]
            / final_stats["D_DNA_m_plus2_isotropic_pos"]["mean_dr_m"]
        ),
    }
    return pd.DataFrame(conv_rows), mc_summary, elapsed_s


def compute_chiral_channels_and_homochiralization(
    nu_summary: dict,
    muon_astro_summary: dict,
    mc_summary: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, dict, float]:
    t0 = time.perf_counter()
    pion_kin = pion_decay_muon_lodacity()

    with open(f"{OUT_DATA}/module1_summary.json", "r") as fh:
        mod1 = json.load(fh)
    g4_row = [r for r in mod1["guanine_hole_trapping"] if r["g_run_length"] == 4][0]
    g4_fold = float(g4_row["fold_enrichment_per_G_vs_G1"])
    g4_capture = float(g4_row["tract_capture_fraction"])
    ac_1uvm_shift = float(mod1["ambient_subgap_1uVm_reference"]["drazin_floquet_subgap_dmu_over_mu"])
    kms_1uvm_shift = float(mod1["ambient_subgap_1uVm_reference"]["adiabatic_kms_ac_dmu_over_mu"])
    met_50mk_floor = float(mod1["metabolic_50mK_noise_floor"])
    tau_mem_fs = float(mod1["tau_mem_fs"])

    m_e_ev = 510998.95
    m_mu_ev = 105.6583755e6
    c_geom = mc_summary["c_geom_calibrated"]
    orient_gain = abs(mc_summary["oriented_to_isotropic_amplification_ratio"])

    def helix_chiral_bias_lepton(t_ev: float, mass_ev: float, l_eff: float) -> float:
        gamma = 1.0 + t_ev / mass_ev
        beta = np.sqrt(1.0 - 1.0 / (gamma ** 2))
        beta_eff = max(beta, 1.5 * ALPHA_FS)
        return float(c_geom * (ALPHA_FS ** 5) * abs(l_eff) * ((1.0 / beta_eff) ** 3))

    l0_sec = 0.20
    bias_1kev_ep = helix_chiral_bias_lepton(1000.0, m_e_ev, l0_sec * np.sqrt(1000.0 / (2.0 * m_e_ev)))
    bias_100ev_ep = helix_chiral_bias_lepton(100.0, m_e_ev, l0_sec * np.sqrt(100.0 / (2.0 * m_e_ev)))
    bias_100ev_ep_oriented = bias_100ev_ep * orient_gain
    bias_2gev_mu = helix_chiral_bias_lepton(2.0e9, m_mu_ev, pion_kin["effective_net_electrostatic_lodacity_mu"])
    bias_lab_500ev = helix_chiral_bias_lepton(500.0, m_e_ev, 0.80) * g4_capture

    # Ionization & 8-oxoG rate budgets per diploid cell (6.4e9 bp ~ 2.6e9 guanines) and per G
    dose_cr_gy_per_yr = 0.34e-3
    dose_cr_gy_per_h = dose_cr_gy_per_yr / (365.25 * 24.0)
    yield_8oxog_per_g_per_gy = 2.0e-7
    k_ox_cr_per_g_h = dose_cr_gy_per_h * yield_8oxog_per_g_per_gy
    k_ox_metabolic_ros_per_g_h = 4.5e-6
    ros_to_cr_ratio = k_ox_metabolic_ros_per_g_h / k_ox_cr_per_g_h

    g_diploid = 2.6e9
    hrs_yr = 365.25 * 24.0
    cr_8oxog_per_cell_yr = k_ox_cr_per_g_h * g_diploid * hrs_yr        # ~ 0.177 8-oxoG/cell/yr
    young_sun_sep_8oxog_per_cell_yr = 1000.0 * cr_8oxog_per_cell_yr   # ~ 177 8-oxoG/cell/yr under active young Sun (10^3x SEP)
    ros_8oxog_per_cell_yr = k_ox_metabolic_ros_per_g_h * g_diploid * hrs_yr

    mod_seasonal_somatic = 0.018 / ros_to_cr_ratio
    mod_solar_11yr_somatic = 0.045 / ros_to_cr_ratio

    channel_rows = [
        {
            "channel": "Ambient_SubThermal_AC_1uVm",
            "quantum_energy_eV": 4.136e-6,
            "regime": "Sub-thermal non-ionizing wave (hv << kBT)",
            "fractional_hazard_or_chiral_shift": ac_1uvm_shift,
            "dna_events_per_diploid_cell_yr": 0.0,
            "dominant_mechanism": f"Exact Floquet-Lindblad Drazin / adiabatic KMS bound (omega << g_PP, tau_mem={tau_mem_fs:.1f} fs)",
            "biological_detectability": "Impossible (~2x10^23x [Floquet] to ~7x10^26x [KMS] below 50 mK floor)",
        },
        {
            "channel": "Solar_Neutrinos_7Be_pp_Doppler",
            "quantum_energy_eV": 8.618e5,
            "regime": "Chiral weak interaction (H_nu = -1, +/-1.43 eV Doppler)",
            "fractional_hazard_or_chiral_shift": float(nu_summary["keplerian_flux_peak_to_trough_pct"] / 100.0),
            "dna_events_per_diploid_cell_yr": float(nu_summary["diploid_DNA_nu_e_scatter_rate_per_yr"]),
            "dominant_mechanism": "Exact 2D Keplerian flux-Doppler clock gated by G_F^2 (sigma ~ 3.6e-45 cm^2)",
            "biological_detectability": (
                f"Idealized one-year count requires {nu_summary['whole_cells_required_5sigma_orbital_clock_1yr']:.2e} whole-cell targets "
                f"or {nu_summary['diploid_DNA_targets_required_5sigma_orbital_clock_1yr']:.2e} diploid DNA targets"
            ),
        },
        {
            "channel": "Intracellular_50mK_Thermal_Noise",
            "quantum_energy_eV": KBT_EV,
            "regime": "Physiological thermal bath (kBT = 26.7 meV)",
            "fractional_hazard_or_chiral_shift": met_50mk_floor,
            "dna_events_per_diploid_cell_yr": float(ros_8oxog_per_cell_yr),
            "dominant_mechanism": "Arrhenius tautomer/repair fluctuation + aerobic metabolic ROS",
            "biological_detectability": "Somatic baseline noise floor",
        },
        {
            "channel": "Modern_SolarModulated_Muons_Somatic",
            "quantum_energy_eV": 2.0e9,
            "regime": "Weak-decay spin-polarized mu+/- (seasonal/27d/11yr)",
            "fractional_hazard_or_chiral_shift": float(mod_solar_11yr_somatic),
            "dna_events_per_diploid_cell_yr": float(cr_8oxog_per_cell_yr),
            "dominant_mechanism": "pi+ -> mu+ + nu_mu EM amplification (10^27x vs nu) diluted 5.8e8x by ROS",
            "biological_detectability": "Masked in modern aerobic cells; active in anoxic deep subsurface",
        },
        {
            "channel": "SpinPolarized_Positron_1keV_Isotropic",
            "quantum_energy_eV": 1000.0,
            "regime": "Supra-thermal chiral ionizing (dE >> I_G = 7.75 eV)",
            "fractional_hazard_or_chiral_shift": float(bias_1kev_ep),
            "dna_events_per_diploid_cell_yr": float(cr_8oxog_per_cell_yr),
            "dominant_mechanism": f"3D helical inward deflection in D-DNA -> central G4 HOMO sink ({g4_fold:.1f}x)",
            "biological_detectability": "Prebiotic homochiralization seed",
        },
        {
            "channel": "YoungSun_SEP_100eV_Mineral_Oriented",
            "quantum_energy_eV": 100.0,
            "regime": "Young-Sun SEP cascade + flat horizontal mineral-adsorbed duplex",
            "fractional_hazard_or_chiral_shift": float(bias_100ev_ep_oriented),
            "dna_events_per_diploid_cell_yr": float(young_sun_sep_8oxog_per_cell_yr),
            "dominant_mechanism": f"{orient_gain:.2f}x flat-terrace cos^2(theta_z) zenith gain + {100.0*g4_capture:.1f}% G4 hole capture + 10^3x Archean SEP flux",
            "biological_detectability": "Prebiotic homochiralization & diurnal solar phase locking",
        },
        {
            "channel": "Lab_Polarized_500eV_Beam_DeltaRepair",
            "quantum_energy_eV": 500.0,
            "regime": "Proposed polarized-lepton experiment",
            "fractional_hazard_or_chiral_shift": float(bias_lab_500ev),
            "dna_events_per_diploid_cell_yr": None,
            "dominant_mechanism": f"Calibrated point-Coulomb spin-orbit estimate {bias_lab_500ev:.2e}",
            "biological_detectability": "Distinguishes CISS/DEA exchange (10^-3) from pure point-Coulomb SO (2.2e-7)",
        },
    ]

    # Prebiotic evolutionary homochiralization (Globus & Blandford Eqs. C34-C35)
    tau_grid = np.logspace(2, 8, 240)
    ee_rows = []
    for dm in [1.0e-7, 1.0e-6, 1.0e-5]:
        ee_no_conflict = np.tanh(0.5 * dm * tau_grid)

        def ode_ee(_t, y, dm_val=dm, kap=12.0 * dm):
            e_val = np.clip(y[0], -0.9999999, 0.9999999)
            return [(0.5 * dm_val + kap * e_val) * (1.0 - e_val ** 2)]

        sol = solve_ivp(
            ode_ee,
            (0.0, float(tau_grid[-1])),
            [0.0],
            t_eval=tau_grid,
            method="RK45",
            rtol=1e-9,
            atol=1e-12,
        )
        for t_val, ee0, eec in zip(tau_grid, ee_no_conflict, sol.y[0]):
            ee_rows.append({
                "delta_M": dm,
                "generations_tau": float(t_val),
                "ee_no_conflict": float(ee0),
                "ee_with_chiral_conflict": float(eec),
            })

    elapsed_s = time.perf_counter() - t0
    cr_summary = {
        "pion_muon_kinematics": pion_kin,
        "g4_hole_sink_enrichment_fold": g4_fold,
        "g4_hole_sink_capture_fraction": g4_capture,
        "ambient_1uVm_floquet_drazin_shift": ac_1uvm_shift,
        "ambient_1uVm_adiabatic_kms_shift": kms_1uvm_shift,
        "chiral_bias_1keV_isotropic_positron": float(bias_1kev_ep),
        "chiral_bias_100eV_isotropic_positron": float(bias_100ev_ep),
        "chiral_bias_100eV_zenith_oriented_positron": float(bias_100ev_ep_oriented),
        "chiral_bias_2GeV_muon_helix": float(bias_2gev_mu),
        "chiral_bias_500eV_lab_polarized_e_G4": float(bias_lab_500ev),
        "sea_level_CR_8oxoG_rate_per_G_h": float(k_ox_cr_per_g_h),
        "metabolic_ROS_8oxoG_rate_per_G_h": float(k_ox_metabolic_ros_per_g_h),
        "metabolic_ROS_to_CR_ratio": float(ros_to_cr_ratio),
        "somatic_8oxoG_shift_seasonal_muon": float(mod_seasonal_somatic),
        "somatic_8oxoG_shift_11yr_solar_cycle": float(mod_solar_11yr_somatic),
        "cr_8oxoG_events_per_diploid_cell_yr": float(cr_8oxog_per_cell_yr),
        "young_sun_1000x_SEP_8oxoG_per_cell_yr": float(young_sun_sep_8oxog_per_cell_yr),
        "tau_99pct_no_conflict_dM_1e6": float(np.arctanh(0.99) * 2.0 / 1.0e-6),
        "min_GtoT_counts_3sigma_at_1e3_asym": float((3.0 / 1.0e-3) ** 2),
        "min_GtoT_counts_3sigma_point_coulomb_500eV": float((3.0 / max(bias_lab_500ev, 1e-30)) ** 2),
    }
    return pd.DataFrame(channel_rows), pd.DataFrame(ee_rows), cr_summary, elapsed_s


def main():
    t_total_start = time.perf_counter()

    df_audit, audit_dict, t_part1 = audit_and_upgrade_github_scripts()
    df_astro, df_diurnal, nu_dict, muon_astro_dict, t_part2 = compute_solar_neutrino_and_muon_astrodynamics()
    df_mc, mc_dict, t_mc = run_3d_monte_carlo_helical_scattering(n_per_ensemble=60000)
    df_chan, df_ee, cr_dict, t_part3 = compute_chiral_channels_and_homochiralization(
        nu_dict, muon_astro_dict, mc_dict
    )

    df_audit.to_csv(f"{OUT_DATA}/module5_github_code_audit_and_upgrade.csv", index=False)
    df_astro.to_csv(f"{OUT_DATA}/module5_astrodynamics_doppler_and_flux.csv", index=False)
    df_diurnal.to_csv(f"{OUT_DATA}/module5_diurnal_24h_solar_clock.csv", index=False)
    df_mc.to_csv(f"{OUT_DATA}/module5_monte_carlo_helical_convergence.csv", index=False)
    df_chan.to_csv(f"{OUT_DATA}/module5_cosmic_ray_chiral_channels.csv", index=False)
    df_ee.to_csv(f"{OUT_DATA}/module5_homochirality_evolution.csv", index=False)

    t_total = time.perf_counter() - t_total_start
    timings = {
        "part1_rel606_2000_windows_s": round(t_part1, 4),
        "part2_8766h_ephemeris_and_neutrinos_s": round(t_part2, 4),
        "part3a_240k_3d_monte_carlo_scattering_s": round(t_mc, 4),
        "part3b_chiral_channels_and_odes_s": round(t_part3, 4),
        "total_module5_wall_clock_s": round(t_total, 4),
    }

    summary = {
        "wall_clock_timings_s": timings,
        "github_audit_and_upgrade": audit_dict,
        "solar_neutrino_astrodynamics": nu_dict,
        "solar_modulated_muon_astrodynamics": muon_astro_dict,
        "monte_carlo_helical_deflection": mc_dict,
        "cosmic_ray_chiral_mutagenesis": cr_dict,
    }
    with open(f"{OUT_DATA}/module5_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print("=== MODULE 5 TIMINGS ===")
    print(json.dumps(timings, indent=2))
    print("\n=== PART 1: REL606 LARGE-SCALE WINDOW AUDIT ===")
    print(json.dumps(audit_dict["real_rel606_windows"], indent=2))
    print("\n=== PART 2A: CHIRAL SOLAR NEUTRINO ASTRODYNAMIC CLOCK & G_F^2 AUDIT ===")
    print(json.dumps({k: v for k, v in nu_dict.items() if k != "branches"}, indent=2))
    print("\n=== PART 2B: SOLAR-MODULATED SPIN-POLARIZED MUON CLOCKS ===")
    print(json.dumps(muon_astro_dict, indent=2))
    print("\n=== PART 3A: 240,000-TRAJECTORY 3D MONTE CARLO HELICAL SCATTERING ===")
    print(json.dumps(mc_dict, indent=2))
    print("\n=== PART 3B: CHIRAL MUTAGENESIS CHANNELS ===")
    print(df_chan[["channel", "quantum_energy_eV", "fractional_hazard_or_chiral_shift", "dna_events_per_diploid_cell_yr"]].to_string(index=False))


if __name__ == "__main__":
    main()
