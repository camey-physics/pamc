#!/usr/bin/env python3

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os

# Parameters
Ls = [4, 6, 8]#, 10]
num_runs = 1
pop_size = 1000
culling_frac = 0.1
beta_max = 0.5

root = "../"
executable = root + "build-release/run_ising"
data_dir = root + "validation/data"
plot_dir = root + "validation/plots"

os.makedirs(data_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

# Known Tc for 3D Ising
Tc = 4.5115
beta_c = 1.0 / Tc

results_binder = {}
results_rho_t = {}

for L in Ls:
    binder_runs = []
    rho_t_runs = []
    step_ref = None
    beta_ref = None

    outfile = os.path.join(data_dir, f"L{L}.txt")

    if os.path.exists(outfile):
        print(f"Loading cached data for L={L}")
        data = np.loadtxt(outfile)
        step_ref = data[:, 0].astype(int)
        beta_ref = data[:, 1]
        binder_mean = data[:, 4]
        rho_t_mean = data[:, 5]
        results_binder[L] = (beta_ref, binder_mean)
        results_rho_t[L] = (step_ref, rho_t_mean)
        continue

    for run_idx in range(num_runs):
        seed = np.random.randint(1e5, 1e9)
        cmd = [
            executable,
            str(L),
            str(pop_size),
            str(culling_frac),
            str(beta_max),
            str(seed)
        ]

        proc = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            universal_newlines=True
        )

        lines = proc.stdout.strip().split("\n")
        data = np.array([
            [float(val) for val in line.split()]
            for line in lines if line.strip()
        ])

        steps = data[:, 0].astype(int)
        beta = data[:, 1]
        binder = data[:, 4]
        rho_t = data[:, 5]

        if beta_ref is None:
            beta_ref = beta
            step_ref = steps
        else:
            assert np.allclose(beta, beta_ref)
            assert np.all(step_ref == steps)

        binder_runs.append(binder)
        rho_t_runs.append(rho_t)

    binder_runs = np.array(binder_runs)
    rho_t_runs = np.array(rho_t_runs)

    binder_mean = np.mean(binder_runs, axis=0)
    rho_t_mean = np.mean(rho_t_runs, axis=0)

    results_binder[L] = (beta_ref, binder_mean)
    results_rho_t[L] = (step_ref, rho_t_mean)

    np.savetxt(outfile, np.column_stack((step_ref, beta_ref, np.mean(binder_runs, axis=0),
                                          np.std(binder_runs, axis=0),
                                          binder_mean, rho_t_mean)),
               header="step beta binder_std binder_std binder_mean rho_t_mean")

# ------------------------------------------
# Plot Binder cumulant crossings
# ------------------------------------------
plt.figure(figsize=(8, 6))

for L in Ls:
    beta, binder_mean = results_binder[L]
    plt.plot(beta, binder_mean, 'o-', label=f"L={L}")

plt.axvline(beta_c, color='k', linestyle='--', label=r"$\beta_c$ (literature)")
plt.xlabel(r"$\beta$")
plt.ylabel(r"Binder Cumulant $U_4$")
plt.title("Binder Cumulant Crossing")
plt.legend()
plt.grid()

binder_plot_file = os.path.join(plot_dir, "binder_crossing.pdf")
plt.savefig(binder_plot_file)
plt.show()

# ------------------------------------------
# Plot rho_t growth verification
# All rho_t comparisons temporarily commented; waiting for full functionality
# ------------------------------------------
# plt.figure(figsize=(8, 6))

# for L in Ls:
#     steps, rho_t_mean = results_rho_t[L]
#     rho_t_theory = 1 + 2 * culling_frac * steps

#     plt.plot(steps, rho_t_mean, 'o-', label=f"L={L}")
#     plt.plot(steps, rho_t_theory, '--', color='gray')

# plt.xlabel("Annealing Step (k)")
# plt.ylabel(r"$\rho_t$")
# plt.title(r"Verification of $\rho_t$ Growth: $\rho_t = 1 + 2 \epsilon k$")
# plt.grid()
# plt.legend()

# rho_t_plot_file = os.path.join(plot_dir, "rho_t_growth.pdf")
# plt.savefig(rho_t_plot_file)
# plt.show()
