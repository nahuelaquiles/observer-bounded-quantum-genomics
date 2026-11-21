#!/usr/bin/env python3
"""Module 1: Open-Quantum Liouvillian, Mori-Zwanzig Resolvent Projection,
Accessible Moment Annihilation, Universal Sub-Thermal Coherent Wave No-Go Bound,
and Guanine Radical-Cation Hole-Sink Trapping.

Establishes exact Kubo-Martin-Schwinger (KMS) detailed balance on solvated
Watson-Crick <-> rare tautomer / anionic base pairs at T = 310.15 K, proves
short-time accessible derivative annihilation through quadratic order, derives
a frequency-general sub-thermal electromagnetic perturbation ceiling across
non-ionizing fields (hbar*omega, lambda << k_B*T), and quantifies duplex
oligo-guanine (G_1..G_5) HOMO hole-sink capture.
"""

import json
import os
import time
import numpy as np
import pandas as pd
from scipy import constants as const
from scipy.linalg import expm, inv, eigvals

OUT_DATA = "/home/sandbox/workspace/observer_bounded_dna_paper/data"
os.makedirs(OUT_DATA, exist_ok=True)

HBAR_EV_S = const.hbar / const.e          # 6.582119569e-16 eV*s
KB_EV_K = const.k / const.e               # 8.617333262e-5 eV/K
TEMP_K = 310.15                           # Physiological 37 C
KBT_EV = KB_EV_K * TEMP_K                 # 0.02672666 eV

# Reference proton-transfer barrier & displacement scale (Slocombe et al. 2021/2022)
V0_AT_EV = 0.065
V0_GC_EV = 0.108
LENGTH_A_M = 1.0e-10                      # 1.0 Angstrom proton transfer coordinate
OMEGA_SUBTHERMAL_LOW = 2.0 * np.pi * 1.0e9   # Representative sub-thermal AC (omega << g_PP: Hz-GHz regime)
OMEGA_THZ_NEAR_GAP = 1.0e13                  # Far-IR / THz regime approaching bath spectral gap (1e13 rad/s)
E0_AMBIENT_REF_VM = 1.0e-6                # Reference weak environmental field (1 uV/m)

# Solvated Watson-Crick <-> rare tautomer / anionic open-quantum 3-state model
DELTA_G_TAUT_EV = 0.420    # Slocombe et al. 2021/2022; Kimsey et al. 2018 (~0.40-0.45 eV)
DELTA_G_ION_EV = 0.540     # Kimsey et al. 2018 pKa-shifted anionic wobble state
J_TUNNEL_EV = 0.012        # Coherent proton-transfer tunnel coupling (eV)
J_ION_EV = 0.006           # Secondary coupling (eV)

GAMMA_DEPH = 1.0e13        # Aqueous H-bond pure dephasing rate in H0 eigenbasis (s^-1), tau_deph = 100 fs
GAMMA_RELAX_10 = 2.5e11    # Downhill vibrational relaxation |E_1> -> |E_0> (s^-1)
GAMMA_RELAX_20 = 1.8e11    # Downhill relaxation |E_2> -> |E_0> (s^-1)
GAMMA_RELAX_21 = 1.2e11    # Downhill relaxation |E_2> -> |E_1> (s^-1)


def liouvillian_from_ham_and_jumps(H_ev, jump_ops):
    d = H_ev.shape[0]
    I_d = np.eye(d, dtype=complex)
    H_s = H_ev / HBAR_EV_S
    L = -1j * (np.kron(I_d, H_s) - np.kron(H_s.T, I_d))
    for rate, c in jump_ops:
        cd_c = c.conj().T @ c
        L += rate * (
            np.kron(c.conj(), c)
            - 0.5 * np.kron(I_d, cd_c)
            - 0.5 * np.kron(cd_c.T, I_d)
        )
    return L


def build_basepair_open_system(delta_H_ev=None):
    """Build Davies-type open-quantum generator where bath jumps satisfy exact KMS
    detailed balance with respect to H = H0 + delta_H_ev, with non-trivial Q <-> P
    blocks in the bare chemical configuration basis."""
    H0 = np.array([
        [0.0, J_TUNNEL_EV, J_ION_EV],
        [J_TUNNEL_EV, DELTA_G_TAUT_EV, 0.003],
        [J_ION_EV, 0.003, DELTA_G_ION_EV],
    ], dtype=complex)
    H = H0 if delta_H_ev is None else (H0 + delta_H_ev)

    evals, evecs = np.linalg.eigh(H)
    jumps = []
    for k in range(3):
        vk = evecs[:, k : k + 1]
        jumps.append((GAMMA_DEPH, vk @ vk.conj().T))

    pair_rates = {(0, 1): GAMMA_RELAX_10, (0, 2): GAMMA_RELAX_20, (1, 2): GAMMA_RELAX_21}
    for (lo, hi), r_down in pair_rates.items():
        dE = float(evals[hi] - evals[lo])
        r_up = r_down * np.exp(-dE / KBT_EV)
        v_lo = evecs[:, lo : lo + 1]
        v_hi = evecs[:, hi : hi + 1]
        jumps.append((r_down, v_lo @ v_hi.conj().T))
        jumps.append((r_up, v_hi @ v_lo.conj().T))

    L0 = liouvillian_from_ham_and_jumps(H, jumps)
    L_unitary = liouvillian_from_ham_and_jumps(H, [])
    return H, evals, evecs, jumps, L0, L_unitary


def verify_resolvent_and_moment_annihilation():
    _, evals, evecs, _, L0, L_unitary = build_basepair_open_system()
    d = 3
    q_idx = [0, 4, 8]
    p_idx = [1, 2, 3, 5, 6, 7]

    Q_proj = np.zeros((9, 9), dtype=complex)
    for idx in q_idx:
        Q_proj[idx, idx] = 1.0

    # Verify exact trace preservation & KMS Gibbs stationarity
    trace_bra = np.array([1, 0, 0, 0, 1, 0, 0, 0, 1], dtype=complex)
    col_sums_max = float(np.max(np.abs(trace_bra @ L0)))

    boltz = np.exp(-(evals - evals[0]) / KBT_EV)
    p_eigen_gibbs = boltz / np.sum(boltz)
    rho_gibbs = sum(p_eigen_gibbs[k] * (evecs[:, k:k+1] @ evecs[:, k:k+1].conj().T) for k in range(d))

    rho0 = np.zeros((3, 3), dtype=complex)
    rho0[0, 0] = 1.0
    vec_ss = expm(L0 * 1.0e-9) @ rho0.reshape(-1, order="F")
    rho_ss = vec_ss.reshape((3, 3), order="F")
    gibbs_trace_err = float(np.linalg.norm(rho_ss - rho_gibbs, ord="fro"))
    p_bare_ss = np.real(np.diag(rho_ss))

    # Schur-complement Mori-Zwanzig effective generator on Q
    L_QQ = L0[np.ix_(q_idx, q_idx)]
    L_QP = L0[np.ix_(q_idx, p_idx)]
    L_PQ = L0[np.ix_(p_idx, q_idx)]
    L_PP = L0[np.ix_(p_idx, p_idx)]

    pp_eigs = eigvals(L_PP)
    spectral_gap_s_inv = float(-np.max(np.real(pp_eigs)))
    L_eff_3x3 = L_QQ - L_QP @ inv(L_PP) @ L_PQ

    L_eff_9x9 = np.zeros((9, 9), dtype=complex)
    L_eff_9x9[np.ix_(q_idx, q_idx)] = L_eff_3x3

    times = np.logspace(-15, -6, 73)
    resolvent_rows = []
    psi_init = np.array([np.sqrt(0.85), np.sqrt(0.10), np.sqrt(0.05)], dtype=complex)
    rho_coh0 = np.outer(psi_init, psi_init.conj()).reshape(-1, order="F")

    for t in times:
        exact_prop_Q = Q_proj @ expm(L0 * t) @ Q_proj
        markov_prop_Q = np.zeros((9, 9), dtype=complex)
        markov_prop_Q[np.ix_(q_idx, q_idx)] = expm(L_eff_3x3 * t)
        op_err = float(np.linalg.norm(exact_prop_Q - markov_prop_Q, ord=2))

        rho_t_open = (expm(L0 * t) @ rho_coh0).reshape((3, 3), order="F")
        rho_t_closed = (expm(L_unitary * t) @ rho_coh0).reshape((3, 3), order="F")
        coh_open = float(np.linalg.norm(rho_t_open - np.diag(np.diag(rho_t_open)), ord="fro"))
        coh_closed = float(np.linalg.norm(rho_t_closed - np.diag(np.diag(rho_t_closed)), ord="fro"))

        resolvent_rows.append({
            "time_s": float(t),
            "op_norm_error": op_err,
            "open_coherence_l2": coh_open,
            "closed_coherence_l2": coh_closed,
            "analytic_upper_bound": float(
                min(2.0, (np.linalg.norm(L_QP, 2) * np.linalg.norm(L_PQ, 2)) / (spectral_gap_s_inv ** 2) * (1.0 / max(spectral_gap_s_inv * t, 1.0) + np.exp(-spectral_gap_s_inv * t)))
            ),
        })

    # Accessible Moment Annihilation Theorem check (Theorem 1): Perturbation purely in P-sector
    rng = np.random.default_rng(20260417)
    A_rand = rng.standard_normal((6, 6)) + 1j * rng.standard_normal((6, 6))
    delta_PP = -0.15 * spectral_gap_s_inv * (A_rand @ A_rand.conj().T) / np.linalg.norm(A_rand, 2) ** 2
    delta_L = np.zeros((9, 9), dtype=complex)
    delta_L[np.ix_(p_idx, p_idx)] = delta_PP

    L_pert = L0 + delta_L
    moment_norms = {}
    for m in range(6):
        diff_m = Q_proj @ (np.linalg.matrix_power(L_pert, m) - np.linalg.matrix_power(L0, m)) @ Q_proj
        moment_norms[f"m={m}"] = float(np.linalg.norm(diff_m, ord=2))

    analytical_m3 = float(np.linalg.norm(L_QP @ delta_PP @ L_PQ, ord=2))

    return {
        "spectral_gap_s_inv": spectral_gap_s_inv,
        "tau_mem_fs": 1e15 / spectral_gap_s_inv,
        "col_sums_max": col_sums_max,
        "gibbs_stationary_fro_err": gibbs_trace_err,
        "eigenstate_gibbs_populations": {
            "E0_ground": float(p_eigen_gibbs[0]),
            "E1_tautomer_eigenstate": float(p_eigen_gibbs[1]),
            "E2_anionic_eigenstate": float(p_eigen_gibbs[2]),
        },
        "bare_basis_populations_with_tunnel_dressing": {
            "WC_bare": float(p_bare_ss[0]),
            "tautomer_bare": float(p_bare_ss[1]),
            "anionic_bare": float(p_bare_ss[2]),
        },
        "L_eff_3x3_s_inv": L_eff_3x3.tolist(),
        "moment_annihilation_norms": moment_norms,
        "analytical_m3_norm": analytical_m3,
        "m3_identity_rel_err": float(abs(moment_norms["m=3"] - analytical_m3) / max(analytical_m3, 1e-30)),
        "resolvent_curve": resolvent_rows,
    }


def compute_weak_field_nogo_curves(spectral_gap_s_inv):
    """Universal sub-thermal electromagnetic perturbation bound across E_0 in [10^-9, 10^5] V/m
    for static (DC), universal sub-gap AC waves (omega << g_PP: Hz-GHz), and near-gap 1.59 THz
    (1.0e13 rad/s) waves, including exact Drazin Floquet-Lindblad Eq. (S2) and adiabatic KMS
    second-order evaluation."""
    _, evals0, evecs0, _, L0, _ = build_basepair_open_system()
    q_idx = [0, 4, 8]
    p_idx = [1, 2, 3, 5, 6, 7]
    L_eff_0 = np.real(L0[np.ix_(q_idx, q_idx)] - L0[np.ix_(q_idx, p_idx)] @ inv(L0[np.ix_(p_idx, p_idx)]) @ L0[np.ix_(p_idx, q_idx)])
    k01_base = float(L_eff_0[1, 0])

    boltz0 = np.exp(-(evals0 - evals0[0]) / KBT_EV)
    p_gibbs0 = boltz0 / np.sum(boltz0)
    p1_gibbs_0 = float(p_gibbs0[1])
    rho_beta = sum(p_gibbs0[k] * (evecs0[:, k:k+1] @ evecs0[:, k:k+1].conj().T) for k in range(3))
    vec_beta = rho_beta.reshape(-1, order="F")
    proj_e1 = evecs0[:, 1:2] @ evecs0[:, 1:2].conj().T

    # Exact Drazin pseudoinverse L0^+ on trace-zero subspace
    trace_bra = np.array([[1, 0, 0, 0, 1, 0, 0, 0, 1]], dtype=complex)
    P_tz = np.eye(9, dtype=complex) - vec_beta[:, None] @ trace_bra
    L0_aug = L0 + 1.0e13 * (vec_beta[:, None] @ trace_bra)
    L0_plus = P_tz @ np.linalg.inv(L0_aug) @ P_tz

    # Unit-coupling (1 eV) proton-transfer dipole operator V_unit = diag(0, -1, 0)
    V_unit = np.diag([0.0, -1.0, 0.0]).astype(complex)
    Lv_unit = liouvillian_from_ham_and_jumps(V_unit, [])
    rhs_unit = Lv_unit @ vec_beta

    def exact_drazin_unit_ratio(omega_rad_s):
        A = np.linalg.solve(1j * omega_rad_s * np.eye(9) - L0_aug, rhs_unit)
        B = np.linalg.solve(-1j * omega_rad_s * np.eye(9) - L0_aug, rhs_unit)
        rho2_unit = (-0.25 * (L0_plus @ (Lv_unit @ (A + B)))).reshape((3, 3), order="F")
        return float(abs(np.real(np.trace(proj_e1 @ rho2_unit))) / p1_gibbs_0)

    drazin_unit_subgap = exact_drazin_unit_ratio(OMEGA_SUBTHERMAL_LOW)
    drazin_unit_thz = exact_drazin_unit_ratio(OMEGA_THZ_NEAR_GAP)

    # Exact adiabatic KMS second-order curvature 0.25 * (1/p1) * d^2 p1 / d(eps)^2 at eps=1e-6 eV
    eps_fd = 1.0e-6
    ev_p = np.linalg.eigvalsh(np.array([[0.0, J_TUNNEL_EV, J_ION_EV], [J_TUNNEL_EV, DELTA_G_TAUT_EV - eps_fd, 0.003], [J_ION_EV, 0.003, DELTA_G_ION_EV]]))
    ev_m = np.linalg.eigvalsh(np.array([[0.0, J_TUNNEL_EV, J_ION_EV], [J_TUNNEL_EV, DELTA_G_TAUT_EV + eps_fd, 0.003], [J_ION_EV, 0.003, DELTA_G_ION_EV]]))
    bp = np.exp(-(ev_p - ev_p[0]) / KBT_EV); p1_p = float(bp[1] / np.sum(bp))
    bm = np.exp(-(ev_m - ev_m[0]) / KBT_EV); p1_m = float(bm[1] / np.sum(bm))
    adiabatic_unit_coeff = float(abs(0.25 * (p1_p + p1_m - 2.0 * p1_gibbs_0) / ((eps_fd ** 2) * p1_gibbs_0)))

    N_MUT_BASELINE_COHORT = 1.0e5
    delta_T_K = 0.05
    metabolic_fractional_noise = (DELTA_G_TAUT_EV / KBT_EV) * (delta_T_K / TEMP_K)

    e_fields = np.logspace(-9, 5, 57)
    rows = []
    for E0 in e_fields:
        # Pure proton transfer dipole coupling lambda = e * a * E0 (in eV, with a = 1.0 Angstrom)
        lam_ev = E0 * LENGTH_A_M
        thermo_bound = 2.0 * lam_ev / KBT_EV  # Theorem 2 Eq. (5) static Kubo ceiling
        lam_freq = lam_ev / HBAR_EV_S

        # Coherence-sector Lorentzian envelope 0.5 * (lambda / (hbar * sqrt(g_PP^2 + omega^2)))^2
        ac_subgap_rectified = 0.5 * (lam_freq ** 2) / (spectral_gap_s_inv ** 2 + OMEGA_SUBTHERMAL_LOW ** 2)
        ac_thz_rectified = 0.5 * (lam_freq ** 2) / (spectral_gap_s_inv ** 2 + OMEGA_THZ_NEAR_GAP ** 2)

        # Exact Drazin Floquet-Lindblad Eq. (S2) & Adiabatic KMS second-order AC shifts
        drazin_subgap = (lam_ev ** 2) * drazin_unit_subgap
        drazin_thz = (lam_ev ** 2) * drazin_unit_thz
        adiabatic_ac = (lam_ev ** 2) * adiabatic_unit_coeff

        # Static DC Gibbs equilibrium shift |dp_1 / p_1|
        dH = np.diag([0.0, -min(lam_ev, 0.25), 0.0]).astype(complex)
        _, evals_dc, _, _, _, _ = build_basepair_open_system(dH)
        bdc = np.exp(-(evals_dc - evals_dc[0]) / KBT_EV)
        p1_dc = float(bdc[1] / np.sum(bdc))
        dc_rel_shift = abs(p1_dc - p1_gibbs_0) / p1_gibbs_0

        snr_ac = np.sqrt(N_MUT_BASELINE_COHORT) * ac_subgap_rectified
        snr_dc = np.sqrt(N_MUT_BASELINE_COHORT) * dc_rel_shift

        rows.append({
            "E0_V_per_m": float(E0),
            "lambda_eV": float(lam_ev),
            "thermo_ceiling_dmu_over_mu": float(thermo_bound),
            "dc_gibbs_dmu_over_mu": float(dc_rel_shift),
            "ac_subgap_rectified_dmu_over_mu": float(ac_subgap_rectified),
            "ac_thz_rectified_dmu_over_mu": float(ac_thz_rectified),
            "drazin_floquet_subgap_dmu_over_mu": float(drazin_subgap),
            "drazin_floquet_thz_dmu_over_mu": float(drazin_thz),
            "adiabatic_kms_ac_dmu_over_mu": float(adiabatic_ac),
            "metabolic_thermal_floor_50mK": float(metabolic_fractional_noise),
            "cohort_poisson_snr_ac_subgap": float(snr_ac),
            "cohort_poisson_snr_dc": float(snr_dc),
        })
    return rows, k01_base, metabolic_fractional_noise


def compute_guanine_hole_sink_trapping():
    """21-bp B-DNA segment containing a central G_n tract (n = 1..5) flanked by
    mixed sequence containing isolated background guanines (at 0.0 eV) and pyrimidine/AT steps.
    Quantifies how long-range radical-cation migration concentrates irreversible
    Bronsted N1 deprotonation + hydration onto G_n low-IP hole sinks."""
    tract_profiles = {
        1: [0.00],
        2: [-0.28, -0.06],
        3: [-0.39, -0.44, -0.09],
        4: [-0.43, -0.51, -0.46, -0.10],
        5: [-0.45, -0.55, -0.53, -0.47, -0.10],
    }
    K_HOP = 2.0e10       # Intrastrand hole hopping rate across bp steps (s^-1)
    K_REC_BG = 4.0e8     # Charge recombination to environment (s^-1)
    K_DEP = 1.8e7        # Irreversible N1 Bronsted deprotonation G.+ -> G(-H). (s^-1; Steenken 1992)
    K_HYD = 6.0e4        # Direct G.+ hydration -> 8-oxoG (s^-1; Rokhlenko et al. 2014)
    F_8OXOG = 0.25       # Fraction of trapped radicals maturing to 8-oxoG

    L_SEG = 21
    results = []
    for n_g, site_e in tract_profiles.items():
        energies = np.array([0.0 if (i % 3 == 0) else 0.16 for i in range(L_SEG)], dtype=float)
        is_g = np.array([(i % 3 == 0) for i in range(L_SEG)], dtype=bool)
        start = (L_SEG - n_g) // 2
        energies[start : start + n_g] = site_e
        is_g[start : start + n_g] = True
        if start - 1 >= 0:
            energies[start - 1] = 0.16; is_g[start - 1] = False
        if start + n_g < L_SEG:
            energies[start + n_g] = 0.16; is_g[start + n_g] = False

        H_tb = np.diag(site_e).astype(float)
        for i in range(n_g - 1):
            H_tb[i, i + 1] = H_tb[i + 1, i] = -0.085
        gs_energy_ev = float(np.linalg.eigvalsh(H_tb)[0])

        R = np.zeros((L_SEG, L_SEG), dtype=float)
        for i in range(L_SEG - 1):
            dE = energies[i + 1] - energies[i]
            r_fwd = K_HOP * np.exp(-0.5 * dE / KBT_EV)
            r_rev = K_HOP * np.exp(+0.5 * dE / KBT_EV)
            R[i + 1, i] += r_fwd; R[i, i] -= r_fwd
            R[i, i + 1] += r_rev; R[i + 1, i + 1] -= r_rev

        k_trap = np.where(is_g, K_DEP + K_HYD, 0.0)
        k_rec = K_REC_BG * np.exp(0.25 * np.minimum(energies, 0.0) / KBT_EV)

        M_kin = R - np.diag(k_rec + k_trap)
        p0 = np.where(is_g, 1.0, 0.0)
        p0 /= np.sum(p0)

        tau_occ = -inv(M_kin) @ p0
        site_8oxog = F_8OXOG * k_trap * tau_occ
        mean_per_g_in_tract = float(np.mean(site_8oxog[start : start + n_g]))
        peak_in_tract = float(np.max(site_8oxog[start : start + n_g]))
        fiveprime_share = float(site_8oxog[start] / np.sum(site_8oxog[start : start + n_g]))

        results.append({
            "g_run_length": n_g,
            "homo_sink_eV": gs_energy_ev,
            "mean_phi_8oxoG_per_tract_G": mean_per_g_in_tract,
            "peak_phi_8oxoG_in_tract": peak_in_tract,
            "total_window_8oxoG_yield": float(np.sum(site_8oxog)),
            "tract_capture_fraction": float(np.sum(site_8oxog[start : start + n_g]) / np.sum(site_8oxog)),
            "fiveprime_G_fraction": fiveprime_share,
            "fold_enrichment_per_G_vs_G1": 1.0,
        })

    g1_ref = results[0]["mean_phi_8oxoG_per_tract_G"]
    for r in results:
        r["fold_enrichment_per_G_vs_G1"] = float(r["mean_phi_8oxoG_per_tract_G"] / g1_ref)
    return results


def main():
    t0 = time.perf_counter()
    res_proj = verify_resolvent_and_moment_annihilation()
    nogo_rows, k01_base, met_floor = compute_weak_field_nogo_curves(res_proj["spectral_gap_s_inv"])
    hole_rows = compute_guanine_hole_sink_trapping()
    elapsed_s = time.perf_counter() - t0

    pd.DataFrame(res_proj["resolvent_curve"]).to_csv(f"{OUT_DATA}/module1_resolvent_convergence.csv", index=False)
    pd.DataFrame(nogo_rows).to_csv(f"{OUT_DATA}/module1_weak_field_nogo.csv", index=False)
    pd.DataFrame(hole_rows).to_csv(f"{OUT_DATA}/module1_guanine_hole_trapping.csv", index=False)

    e_ref_row = min(nogo_rows, key=lambda r: abs(np.log10(r["E0_V_per_m"]) - (-6.0)))
    e_cross_ac = min(nogo_rows, key=lambda r: abs(np.log10(max(r["ac_subgap_rectified_dmu_over_mu"], 1e-40)) - np.log10(met_floor)))
    summary = {
        "runtime_seconds": elapsed_s,
        "temperature_K": TEMP_K,
        "kBT_eV": KBT_EV,
        "spectral_gap_s_inv": res_proj["spectral_gap_s_inv"],
        "tau_mem_fs": res_proj["tau_mem_fs"],
        "trace_conservation_err": res_proj["col_sums_max"],
        "gibbs_stationary_fro_err": res_proj["gibbs_stationary_fro_err"],
        "eigenstate_gibbs_populations": res_proj["eigenstate_gibbs_populations"],
        "bare_basis_populations_with_tunnel_dressing": res_proj["bare_basis_populations_with_tunnel_dressing"],
        "k01_effective_s_inv": k01_base,
        "moment_annihilation_norms": res_proj["moment_annihilation_norms"],
        "m3_identity_rel_err": res_proj["m3_identity_rel_err"],
        "resolvent_error_at_1ps": next(r["op_norm_error"] for r in res_proj["resolvent_curve"] if r["time_s"] >= 1e-12),
        "resolvent_error_at_1ns": next(r["op_norm_error"] for r in res_proj["resolvent_curve"] if r["time_s"] >= 1e-9),
        "ambient_subgap_1uVm_reference": e_ref_row,
        "approx_crossover_E0_V_per_m": e_cross_ac["E0_V_per_m"],
        "metabolic_50mK_noise_floor": met_floor,
        "guanine_hole_trapping": hole_rows,
    }
    with open(f"{OUT_DATA}/module1_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"=== MODULE 1 SUMMARY (runtime = {elapsed_s:.3f} s) ===")
    print(f"Spectral gap g(L_PP) = {summary['spectral_gap_s_inv']:.4e} s^-1 (tau_mem = {summary['tau_mem_fs']:.2f} fs)")
    print(f"Gibbs stationary Frobenius error ||rho_ss - rho_Gibbs||_F = {summary['gibbs_stationary_fro_err']:.3e}")
    print(f"Resolvent error at 1 ps = {summary['resolvent_error_at_1ps']:.3e}, at 1 ns = {summary['resolvent_error_at_1ns']:.3e}")
    print("Moment annihilation norms:", summary["moment_annihilation_norms"])
    print(f"m=3 Schur identity rel err = {summary['m3_identity_rel_err']:.3e}")
    print(f"Universal sub-gap AC (omega << g_PP) at 1 uV/m: lambda = {e_ref_row['lambda_eV']:.3e} eV, static ceiling = {e_ref_row['thermo_ceiling_dmu_over_mu']:.3e}, AC rectified = {e_ref_row['ac_subgap_rectified_dmu_over_mu']:.3e}, SNR = {e_ref_row['cohort_poisson_snr_ac_subgap']:.3e}")
    print(pd.DataFrame(hole_rows).to_string(index=False))


if __name__ == "__main__":
    main()
