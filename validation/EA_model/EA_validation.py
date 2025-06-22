#!/usr/bin/env python3

import subprocess
import os
import numpy as np
import json

# Input configuration
L = 6
pop_size = 50_000
culling_frac = 0.1
beta_max = 5.0
seed = 94723975
num_threads = 4

root = "../../"
instance_dir = f"EA_realizations/L{L}"
bond_path = os.path.join(instance_dir, "bonds.txt")
neighbor_path = os.path.join(instance_dir, "neighbors.txt")
stats_path = os.path.join(instance_dir, f"EA_L{L}_stats.json")

executable = os.path.join(root, "build-release", "run_3D_EA")

# Load benchmark ground truth stats
with open(stats_path, "r") as f:
    benchmark = json.load(f)["ground_truth"]
    gs_energy_ref = benchmark["energy_min"]
    rho_t_ref = benchmark["rho_t"]
    num_fams_ref = benchmark["num_gs_families"]

# Run executable
cmd = [
    executable,
    str(L),
    str(pop_size),
    str(culling_frac),
    str(beta_max),
    str(seed),
    neighbor_path,
    bond_path,
    str(num_threads)
]

print(f"Running: {' '.join(cmd)}")

proc = subprocess.run(
    cmd,
    check=True,
    stdout=subprocess.PIPE,
    universal_newlines=True
)

# Parse final line of output
last_line = proc.stdout.strip().split("\n")[-1]
fields = last_line.split()
E_min = float(fields[3])
rho_t = float(fields[4])
num_gs_families = int(fields[5])

# Compare and report
print("\nValidation Summary:")
print(f"  Reference GS Energy:     {gs_energy_ref:.15f}")
print(f"  Measured GS Energy:      {E_min:.15f}")
print(f"  Match:                   {'YES' if np.isclose(E_min, gs_energy_ref, atol=1e-14) else 'NO'}")
