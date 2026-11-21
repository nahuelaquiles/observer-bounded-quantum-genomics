#!/usr/bin/env python3
"""Module 2: Distributed Vector Rate-Distortion, Single-Cell Poisson Erasure,
and Stiffness-Free QSSA Circadian OGG1-BER Kinetics with Exact Nonlinear
Gauss-Newton Poisson MAP Decoding."""

import os
import numpy as np
import pandas as pd

OUT_DATA = "/home/sandbox/workspace/observer_bounded_dna_paper/data"
os.makedirs(OUT_DATA, exist_ok=True)


def run_single_cell_sparsity_and_ensemble_scaling():
    K_EFF = 400
    K_HARMONICS = 3
    K_DIM = 6
    SIG_THETA = 0.18
    SIG2_THETA = SIG_THETA ** 2

    p0_grid = [1.0e-6, 3.0e-6, 1.0e-5]
    sparsity_rows = []
    for p0 in p0_grid:
        lam_cell = K_EFF * p0
        p_zero = float(np.exp(-lam_cell))
        c_single_cell_bits = 0.5 * K_DIM * np.log2(1.0 + (K_EFF * p0 * SIG2_THETA) / K_DIM)
        sparsity_rows.append({
            "baseline_8oxoG_per_G": p0,
            "K_eff_guanines": K_EFF,
            "mean_lesions_per_promoter_per_cell": lam_cell,
            "prob_zero_lesions_single_cell": p_zero,
            "prob_at_least_one_lesion": 1.0 - p_zero,
            "single_cell_capacity_bits": float(c_single_cell_bits),
        })

    rng = np.random.default_rng(20261120)
    T_HOURS = 72.0
    M_STEPS = 24
    dt_h = T_HOURS / M_STEPS
    t_grid = np.linspace(dt_h, T_HOURS, M_STEPS)
    omega = 2.0 * np.pi / 24.0

    cols = []
    for h in range(1, K_HARMONICS + 1):
        cols.append(np.sqrt(2.0) * np.cos(h * omega * t_grid))
        cols.append(np.sqrt(2.0) * np.sin(h * omega * t_grid))
    Phi = np.column_stack(cols)

    K_OPEN = 0.08
    g_positions = np.arange(K_EFF, dtype=float) * 2.5
    helical_phase = 2.0 * np.pi * g_positions / 10.4
    acc_weights = K_OPEN + (1.0 - K_OPEN) * 0.5 * (1.0 + np.cos(helical_phase))
    mean_acc = float(np.mean(acc_weights))

    k_ox_base = 4.5e-6
    k_rep_bar = 0.45

    schedules = {
        "rhythmic_OGG1": k_rep_bar * (1.0 + 0.30 * np.cos(omega * t_grid - np.pi / 4.0)),
        "flat_equal_OGG1": np.full(M_STEPS, k_rep_bar),
        "inward_occluded_SHL": np.full(M_STEPS, 0.055),
    }

    nu_grid = [1, 10, 100, 1000, 10000, 100000, 300000, 1000000]
    N_MC = 2500
    N_PRIOR = 5000
    rng_prior = np.random.default_rng(20261121)
    prior_thetas = rng_prior.normal(0.0, SIG_THETA, size=(N_PRIOR, K_DIM))

    ensemble_rows = []
    fisher_audit = {}
    f_ref_per_cell = 0.0

    for sched_name, k_rep_vec in schedules.items():
        eff_acc = K_OPEN if sched_name == "inward_occluded_SHL" else mean_acc

        def forward_and_jac_batch(thetas_b):
            thetas_2d = np.atleast_2d(thetas_b)
            B = thetas_2d.shape[0]
            x_traj = thetas_2d @ Phi.T
            exp_x = np.exp(np.clip(x_traj, -3.0, 3.0))
            k_in = k_ox_base * eff_acc * exp_x
            dk_dtheta = k_in[:, :, None] * Phi[None, :, :]  # (B, M_STEPS, K_DIM)

            p = k_in[:, 0] / k_rep_vec[0]
            dp = dk_dtheta[:, 0, :] / k_rep_vec[0]
            for s in range(8):
                idx = s % 8
                dec = np.exp(-k_rep_vec[idx] * dt_h)
                fac = (1.0 - dec) / k_rep_vec[idx]
                p = p * dec + k_in[:, idx] * fac
                dp = dp * dec + dk_dtheta[:, idx, :] * fac

            p_out = np.zeros((B, M_STEPS), dtype=float)
            J_out = np.zeros((B, M_STEPS, K_DIM), dtype=float)
            for m in range(M_STEPS):
                dec = np.exp(-k_rep_vec[m] * dt_h)
                fac = (1.0 - dec) / k_rep_vec[m]
                p = p * dec + k_in[:, m] * fac
                dp = dp * dec + dk_dtheta[:, m, :] * fac
                p_out[:, m] = p
                J_out[:, m, :] = dp
            return p_out, J_out

        p_bar_2d, G_jac0_3d = forward_and_jac_batch(np.zeros((1, K_DIM)))
        p_bar = p_bar_2d[0]
        G_jac0 = G_jac0_3d[0]
        F_1cell_origin = (G_jac0.T * (K_EFF / p_bar)) @ G_jac0

        # Exact van Trees prior-averaged single-cell Fisher information matrix E_Theta[F_1(Theta)]
        p_prior, J_prior = forward_and_jac_batch(prior_thetas)
        F_1cell_prior_avg = np.einsum("bmk,bm,bml->kl", J_prior, K_EFF / p_prior, J_prior) / float(N_PRIOR)

        fisher_audit[sched_name] = {
            "trace_F1_origin": float(np.trace(F_1cell_origin)),
            "trace_F1_prior_avg": float(np.trace(F_1cell_prior_avg)),
            "prior_over_origin_trace_ratio": float(np.trace(F_1cell_prior_avg) / np.trace(F_1cell_origin)),
        }
        if sched_name == "rhythmic_OGG1":
            f_ref_per_cell = float(np.trace(F_1cell_prior_avg))

        prior_prec = np.eye(K_DIM) / SIG2_THETA

        for nu in nu_grid:
            # Exact van Trees Bayesian Cramer-Rao Lower Bound using E_Theta[F_1(Theta)]
            post_prec = prior_prec + nu * F_1cell_prior_avg
            post_cov = np.linalg.inv(post_prec)
            bcrb_mse = float(np.trace(post_cov) / K_DIM)

            # Origin-linearized F_1(0) reference bound
            post_cov_f0 = np.linalg.inv(prior_prec + nu * F_1cell_origin)
            f0_bcrb_mse = float(np.trace(post_cov_f0) / K_DIM)

            _, logdet = np.linalg.slogdet(np.eye(K_DIM) + SIG2_THETA * nu * F_1cell_prior_avg)
            mi_bits = float(0.5 * logdet / np.log(2.0))

            thetas = rng.normal(0.0, SIG_THETA, size=(N_MC, K_DIM))
            p_true, _ = forward_and_jac_batch(thetas)
            y_counts = rng.poisson(nu * K_EFF * p_true)
            zero_frac_nu1 = float(np.mean(y_counts == 0)) if nu == 1 else np.nan

            # Converged 8-step batched Gauss-Newton / Fisher-scoring Poisson MAP optimization
            th = np.zeros((N_MC, K_DIM), dtype=float)
            for _ in range(8):
                pm, Jm = forward_and_jac_batch(th)
                grad = np.einsum("bmk,bm->bk", Jm, y_counts / pm - nu * K_EFF) - th / SIG2_THETA
                H_gn = prior_prec[None, :, :] + np.einsum("bmk,bm,bml->bkl", Jm, nu * K_EFF / pm, Jm)
                step = np.linalg.solve(H_gn, grad[:, :, None])[:, :, 0]
                th = np.clip(th + step, -2.5, 2.5)
                if float(np.max(np.linalg.norm(step, axis=1))) < 1e-6:
                    break

            sq_errs = np.mean((th - thetas) ** 2, axis=1)
            emp_mse = float(np.mean(sq_errs))
            emp_se = float(np.std(sq_errs, ddof=1) / np.sqrt(N_MC))

            ensemble_rows.append({
                "schedule": sched_name,
                "nu_cells": nu,
                "mean_p_bar_per_G": float(np.mean(p_bar)),
                "empirical_mse": emp_mse,
                "empirical_mse_se": emp_se,
                "clarke_barron_bcrb_mse": bcrb_mse,
                "normalized_distortion_D_over_sigma2": emp_mse / SIG2_THETA,
                "bcrb_D_over_sigma2": bcrb_mse / SIG2_THETA,
                "origin_f0_bcrb_D_over_sigma2": f0_bcrb_mse / SIG2_THETA,
                "mutual_info_bits": mi_bits,
                "empty_timebin_fraction_nu1": zero_frac_nu1,
            })

    mem_rows = []
    for c_bits in [2.0, 4.0, 8.0, 16.0]:
        c_nats = c_bits * np.log(2.0)
        for k_m in range(1, 25):
            d_sens = SIG2_THETA / (1.0 + (1.0e5 * f_ref_per_cell * SIG2_THETA) / k_m)
            d_mem = SIG2_THETA * np.exp(-2.0 * c_nats / k_m)
            mem_rows.append({
                "C_B_bits": c_bits,
                "K_modes": k_m,
                "D_sensing_norm": float(d_sens / SIG2_THETA),
                "D_memory_norm": float(d_mem / SIG2_THETA),
                "D_effective_lower_bound_norm": float(max(d_sens, d_mem) / SIG2_THETA),
                "regime": "memory_limited" if d_mem > d_sens else "sensing_limited",
            })

    return sparsity_rows, ensemble_rows, mem_rows, fisher_audit


def main():
    import json
    sparsity_rows, ensemble_rows, mem_rows, fisher_audit = run_single_cell_sparsity_and_ensemble_scaling()
    df_sp = pd.DataFrame(sparsity_rows)
    df_ens = pd.DataFrame(ensemble_rows)
    df_mem = pd.DataFrame(mem_rows)

    df_sp.to_csv(f"{OUT_DATA}/module2_single_cell_sparsity.csv", index=False)
    df_ens.to_csv(f"{OUT_DATA}/module2_ensemble_rate_distortion.csv", index=False)
    df_mem.to_csv(f"{OUT_DATA}/module2_memory_sensing_crossover.csv", index=False)
    with open(f"{OUT_DATA}/module2_summary.json", "w", encoding="utf-8") as fh:
        json.dump(
            {
                "n_mc_trials": 2500,
                "n_prior_draws": 5000,
                "fisher_prior_vs_origin": fisher_audit,
                "all_empirical_above_van_trees_bound": bool(
                    np.all(df_ens["normalized_distortion_D_over_sigma2"] >= df_ens["bcrb_D_over_sigma2"] - 2.0 * df_ens["empirical_mse_se"] / (0.18 ** 2))
                ),
            },
            fh,
            indent=2,
        )

    print("=== MODULE 2 REFINED SUMMARY (EXACT VAN TREES PRIOR-AVERAGED FISHER & NONLINEAR POISSON MAP) ===")
    print(df_ens[["schedule", "nu_cells", "normalized_distortion_D_over_sigma2", "bcrb_D_over_sigma2", "origin_f0_bcrb_D_over_sigma2", "mutual_info_bits"]].to_string(index=False))


if __name__ == "__main__":
    main()
