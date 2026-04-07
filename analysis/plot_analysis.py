# plot_analysis.py
# ------------------------------------------------------------------
# Quantitative analysis of the TDSE simulation.  Four plots on one
# figure:
#
#   (A) Norm conservation:  integral |psi|^2 dx vs time.
#       Should stay = 1 throughout.  Any drift reveals numerical error.
#
#   (B) Transmission T(t) and reflection R(t):
#       T = integral_{barrier}^{+inf} |psi|^2 dx
#       R = integral_{-inf}^{barrier} |psi|^2 dx
#       We also plot T + R to verify T + R ≈ 1  (unitarity check).
#
#   (C) Mean position <x>(t) = integral x |psi|^2 dx
#       Shows the centre of mass moving, slowing at the barrier,
#       and (if partially transmitted) splitting.
#
#   (D) Width sigma(t) = sqrt( <x^2> - <x>^2 )
#       A free Gaussian spreads linearly in time.  The barrier causes
#       additional spreading (interference between R and T components).
#
# HOW TO RUN (from the project root):
#   python analysis/plot_analysis.py
#
# REQUIREMENTS:
#   pip install matplotlib numpy
# ------------------------------------------------------------------

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# 1. LOCATE SNAPSHOT FILES AND NORM LOG
# ------------------------------------------------------------------
output_dir = "output"
files = sorted(glob.glob(os.path.join(output_dir, "psi_*.csv")))

if len(files) == 0:
    raise FileNotFoundError(
        "No snapshot files found in 'output/'. "
        "Run the C++ solver first with:  make && ./tdse"
    )

norm_file = os.path.join(output_dir, "norm.csv")
norm_data = np.loadtxt(norm_file, delimiter=",", skiprows=1)

# norm_data[:,0] = time,  norm_data[:,1] = norm^2 at each snapshot
times  = norm_data[:, 0]
norm2  = norm_data[:, 1]

# ------------------------------------------------------------------
# 2. FIND THE BARRIER POSITION
# ------------------------------------------------------------------
# We read the first snapshot to locate where V > 0.
# The barrier left edge is where the potential first becomes non-zero.
data0 = np.loadtxt(files[0], delimiter=",", skiprows=1)
x_ref = data0[:, 0]
V_ref = data0[:, 4]

# Index of the first grid point inside the barrier
barrier_mask  = V_ref > 0.0
barrier_idx   = np.argmax(barrier_mask)   # first True index
x_bar_left    = x_ref[barrier_idx]        # left edge of barrier

print(f"Barrier detected at x = {x_bar_left:.3f} bohr (grid index {barrier_idx})")

# ------------------------------------------------------------------
# 3. LOOP OVER ALL SNAPSHOTS AND COMPUTE OBSERVABLES
# ------------------------------------------------------------------
# We iterate over every snapshot file and compute T, R, <x>, sigma.
# This loop is the main computational work of this script.

dx = x_ref[1] - x_ref[0]   # grid spacing (uniform)

T_arr    = []   # transmission coefficient at each snapshot
R_arr    = []   # reflection  coefficient at each snapshot
x_mean   = []   # mean position <x>
x_width  = []   # width sigma = sqrt(<x^2> - <x>^2)

for fpath in files:
    data = np.loadtxt(fpath, delimiter=",", skiprows=1)
    x    = data[:, 0]
    prob = data[:, 1]   # |psi|^2

    # ---- Transmission and reflection ----
    # Integrate |psi|^2 to the left (reflection) and right (transmission)
    # of the barrier's left edge.
    R_arr.append(np.sum(prob[:barrier_idx]) * dx)
    T_arr.append(np.sum(prob[barrier_idx:]) * dx)

    # ---- Mean position ----
    # <x> = integral x |psi|^2 dx  ≈  sum_i x_i * prob_i * dx
    xmean = np.sum(x * prob) * dx
    x_mean.append(xmean)

    # ---- Width ----
    # sigma = sqrt( <x^2> - <x>^2 )
    x2mean = np.sum(x**2 * prob) * dx
    x_width.append(np.sqrt(max(x2mean - xmean**2, 0.0)))

T_arr   = np.array(T_arr)
R_arr   = np.array(R_arr)
x_mean  = np.array(x_mean)
x_width = np.array(x_width)

# norm.csv has one row per snapshot (t=0 included), in the same order
# as the psi_*.csv files, so the mapping is direct.
snap_times = times

# ------------------------------------------------------------------
# 4. FINAL VALUES (printed as a summary)
# ------------------------------------------------------------------
T_final = T_arr[-1]
R_final = R_arr[-1]
print(f"\n--- Final state (t = {snap_times[-1]:.4f} a.u.) ---")
print(f"  Transmission T = {T_final:.6f}")
print(f"  Reflection   R = {R_final:.6f}")
print(f"  T + R          = {T_final + R_final:.6f}  (should be ~1)")
print(f"  Norm drift     = {abs(norm2[-1] - norm2[0]):.2e}")

# ------------------------------------------------------------------
# 5. PLOT — 2x2 GRID OF SUBPLOTS
# ------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle("TDSE Crank-Nicolson — quantitative analysis", fontsize=14)

# Unpack axes for readability
ax_norm = axes[0, 0]
ax_TR   = axes[0, 1]
ax_xmn  = axes[1, 0]
ax_wid  = axes[1, 1]

# ---- (A) Norm conservation ----
ax_norm.plot(times, norm2, color="royalblue", lw=1.5)
ax_norm.axhline(1.0, color="gray", lw=0.8, linestyle="--", label="exact = 1")
ax_norm.set_xlabel("t  [a.u.]")
ax_norm.set_ylabel("∫ |ψ|² dx")
ax_norm.set_title("(A) Norm conservation")
ax_norm.legend(fontsize=9)
# Show the absolute drift on the plot
drift = np.abs(norm2 - 1.0).max()
ax_norm.text(0.98, 0.05, f"max |drift| = {drift:.1e}",
             ha="right", va="bottom", transform=ax_norm.transAxes, fontsize=8)

# ---- (B) Transmission and reflection ----
ax_TR.plot(snap_times, T_arr, color="seagreen",  lw=1.5, label="T (transmitted)")
ax_TR.plot(snap_times, R_arr, color="tomato",    lw=1.5, label="R (reflected)")
ax_TR.plot(snap_times, T_arr + R_arr, color="gray", lw=1.0,
           linestyle=":", label="T + R (unitarity)")
ax_TR.axhline(1.0, color="gray", lw=0.6, linestyle="--")
ax_TR.set_xlabel("t  [a.u.]")
ax_TR.set_ylabel("probability")
ax_TR.set_title("(B) Transmission T and reflection R")
ax_TR.legend(fontsize=9)
ax_TR.set_ylim(-0.05, 1.15)

# ---- (C) Mean position <x>(t) ----
ax_xmn.plot(snap_times, x_mean, color="darkorange", lw=1.5)
# Mark the barrier position
ax_xmn.axvline(x=0, color="white", lw=0)   # invisible, just for spacing
ax_xmn.axhline(y=x_bar_left, color="tomato", lw=0.8,
               linestyle="--", label=f"barrier @ x={x_bar_left:.1f}")
ax_xmn.set_xlabel("t  [a.u.]")
ax_xmn.set_ylabel("⟨x⟩  [bohr]")
ax_xmn.set_title("(C) Mean position ⟨x⟩(t)")
ax_xmn.legend(fontsize=9)

# ---- (D) Wave packet width sigma(t) ----
ax_wid.plot(snap_times, x_width, color="mediumpurple", lw=1.5)
ax_wid.set_xlabel("t  [a.u.]")
ax_wid.set_ylabel("σ(t)  [bohr]")
ax_wid.set_title("(D) Packet width σ(t) = √(⟨x²⟩ − ⟨x⟩²)")

plt.tight_layout()

# ------------------------------------------------------------------
# 6. SAVE AND SHOW
# ------------------------------------------------------------------
out_path = os.path.join(output_dir, "analysis.pdf")
plt.savefig(out_path)
print(f"\nFigure saved to {out_path}")
plt.show()
