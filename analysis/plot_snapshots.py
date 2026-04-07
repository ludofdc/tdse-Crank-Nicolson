# plot_snapshots.py
# ------------------------------------------------------------------
# Produces a figure with multiple panels showing the probability
# density |psi(x,t)|^2 at a few chosen instants in time.
# This gives an immediate visual picture of how the wave packet
# travels, hits the barrier, and splits into reflected + transmitted
# components.
#
# HOW TO RUN (from the project root):
#   python analysis/plot_snapshots.py
#
# REQUIREMENTS:
#   pip install matplotlib numpy
# ------------------------------------------------------------------

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# 1. LOCATE SNAPSHOT FILES
# ------------------------------------------------------------------
# All snapshot files live in output/ and are named psi_000000.csv,
# psi_000001.csv, ...  glob.glob finds them all, then we sort them
# so they are in chronological order.
output_dir = "output"
files = sorted(glob.glob(os.path.join(output_dir, "psi_*.csv")))

if len(files) == 0:
    raise FileNotFoundError(
        "No snapshot files found in 'output/'. "
        "Run the C++ solver first with:  make && ./tdse"
    )

# ------------------------------------------------------------------
# 2. CHOOSE WHICH SNAPSHOTS TO PLOT
# ------------------------------------------------------------------
# We pick 5 snapshots spread evenly: t=0, 25%, 50%, 75%, 100%.
# np.linspace gives us 5 equally-spaced indices in [0, N-1].
n_panels = 5
indices = np.linspace(0, len(files) - 1, n_panels, dtype=int)
selected_files = [files[i] for i in indices]

# ------------------------------------------------------------------
# 3. READ THE NORM LOG to get the time axis
# ------------------------------------------------------------------
# norm.csv has two columns: t (time) and norm^2 (should stay ~1).
norm_file = os.path.join(output_dir, "norm.csv")
norm_data = np.loadtxt(norm_file, delimiter=",", skiprows=1)
times_all = norm_data[:, 0]   # all saved times

# Map each selected index to the corresponding simulation time.
# norm.csv includes one row per snapshot file (step 0 included),
# so file k maps directly to times_all[k].
selected_times = []
for idx in indices:
    if idx < len(times_all):
        selected_times.append(times_all[idx])
    else:
        selected_times.append(times_all[-1])

# ------------------------------------------------------------------
# 4. CREATE THE FIGURE
# ------------------------------------------------------------------
# One row of n_panels side-by-side subplots, sharing the y axis so
# the amplitudes are directly comparable across panels.
fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, 4),
                         sharey=True)
fig.suptitle("Wave packet evolution  |ψ(x,t)|²", fontsize=14)

# -- Precompute a single V scale factor shared across all panels --
# Use the global peak over ALL snapshot files (not just the 5 selected)
# so the barrier height is consistent with the animation, which also
# scans all loaded frames.
global_prob_max = max(
    np.loadtxt(f, delimiter=",", skiprows=1)[:, 1].max()
    for f in files
)
global_prob_max = global_prob_max if global_prob_max > 0 else 1.0

# Read V from the first file (it is identical in every snapshot)
V_ref_data = np.loadtxt(selected_files[0], delimiter=",", skiprows=1)
V_max_global = V_ref_data[:, 4].max()

for ax, fpath, t_val in zip(axes, selected_files, selected_times):

    # -- Load CSV columns:  x, |psi|^2, Re(psi), Im(psi), V(x) --
    data = np.loadtxt(fpath, delimiter=",", skiprows=1)
    x    = data[:, 0]   # spatial grid
    prob = data[:, 1]   # |psi|^2  (probability density)
    V    = data[:, 4]   # potential

    # -- Rescale potential for plotting on the same axis --
    # Same formula as the animation: peak of V → 0.6 * global_prob_max,
    # so the barrier height is directly comparable between the two plots.
    V_scale = V / V_max_global * global_prob_max * 0.6 if V_max_global > 0 else V

    # -- Plot probability density --
    ax.plot(x, prob, color="royalblue", lw=1.5, label="|ψ|²")

    # -- Shade the potential barrier --
    ax.fill_between(x, 0, V_scale, color="tomato", alpha=0.35, label="V(x)")

    # -- Labels and title --
    ax.set_title(f"t = {t_val:.3f} a.u.", fontsize=10)
    ax.set_xlabel("x  [bohr]")
    if ax is axes[0]:
        ax.set_ylabel("|ψ(x,t)|²")
    ax.set_xlim(x[0], x[-1])

# Add a single legend outside the last panel
axes[-1].legend(loc="upper right", fontsize=8)

plt.tight_layout()

# ------------------------------------------------------------------
# 5. SAVE AND SHOW
# ------------------------------------------------------------------
out_path = os.path.join(output_dir, "snapshots.pdf")
plt.savefig(out_path)
print(f"Figure saved to {out_path}")
plt.show()
