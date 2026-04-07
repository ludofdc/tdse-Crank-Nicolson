# plot_animation.py
# ------------------------------------------------------------------
# Produces an animated MP4 (or GIF) showing the full time evolution
# of the wave packet.  Two quantities are animated together:
#   - |psi(x,t)|^2   (probability density, filled)
#   - Re(psi(x,t))   (real part, dashed)
# The potential barrier is shown as a static shaded region.
#
# HOW TO RUN (from the project root):
#   python analysis/plot_animation.py
#
# REQUIREMENTS:
#   pip install matplotlib numpy
#   For MP4 output you also need ffmpeg installed on your system.
#   If ffmpeg is missing the script falls back to GIF (Pillow needed):
#   pip install pillow
# ------------------------------------------------------------------

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ------------------------------------------------------------------
# 1. LOAD ALL SNAPSHOT FILES
# ------------------------------------------------------------------
output_dir = "output"
files = sorted(glob.glob(os.path.join(output_dir, "psi_*.csv")))

if len(files) == 0:
    raise FileNotFoundError(
        "No snapshot files found in 'output/'. "
        "Run the C++ solver first with:  make && ./tdse"
    )

# To avoid very large animations we cap the number of frames.
# If there are more snapshots than MAX_FRAMES we subsample them.
# We keep track of the original indices so we can look up the correct
# simulation time for each selected file.
MAX_FRAMES = 120
n_total = len(files)
if n_total > MAX_FRAMES:
    orig_indices = np.linspace(0, n_total - 1, MAX_FRAMES, dtype=int)
    files = [files[i] for i in orig_indices]
else:
    orig_indices = np.arange(n_total)

# ------------------------------------------------------------------
# 2. READ THE TIME AXIS FROM norm.csv
# ------------------------------------------------------------------
norm_file = os.path.join(output_dir, "norm.csv")
norm_data = np.loadtxt(norm_file, delimiter=",", skiprows=1)
# norm_data[:,0] is time — one row per snapshot file (step 0 included),
# so file k maps directly to norm_data[k].
all_times = norm_data[:, 0]
n_all     = len(all_times)

# Build a time array aligned with the selected files.
frame_times = [all_times[min(i, n_all - 1)] for i in orig_indices]

# ------------------------------------------------------------------
# 3. LOAD ALL FRAMES INTO MEMORY
# ------------------------------------------------------------------
# Each frame is a dict with arrays x, prob, re_psi, V.
# Loading everything upfront makes the animation smoother.
frames = []
for fpath in files:
    data = np.loadtxt(fpath, delimiter=",", skiprows=1)
    frames.append({
        "x":      data[:, 0],
        "prob":   data[:, 1],   # |psi|^2
        "re_psi": data[:, 2],   # Re(psi)
        "V":      data[:, 4],   # potential
    })

x = frames[0]["x"]
V = frames[0]["V"]

# ------------------------------------------------------------------
# 4. SET UP THE FIGURE AND STATIC ELEMENTS
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 4))
ax.set_xlabel("x  [bohr]", fontsize=11)
ax.set_ylabel("|ψ(x,t)|²  /  Re(ψ)", fontsize=11)
ax.set_xlim(x[0], x[-1])
ax.set_title("TDSE Crank-Nicolson — wave packet evolution", fontsize=12)

# Determine y limits from the maximum probability over all frames.
prob_max = max(f["prob"].max() for f in frames)
re_max   = max(abs(f["re_psi"]).max() for f in frames)
y_max    = max(prob_max, re_max) * 1.15   # add 15% headroom
ax.set_ylim(-y_max * 0.5, y_max)

# Draw potential barrier (static, plotted once).
# We rescale V so its peak reaches half the maximum probability.
if V.max() > 0:
    V_plot = V / V.max() * prob_max * 0.6
else:
    V_plot = V
ax.fill_between(x, 0, V_plot, color="tomato", alpha=0.30, label="V(x)")

# Create the animated line objects (empty at first).
line_prob, = ax.plot([], [], color="royalblue", lw=2.0,  label="|ψ|²")
line_re,   = ax.plot([], [], color="darkorange", lw=1.2,
                     linestyle="--", label="Re(ψ)")
time_text  = ax.text(0.02, 0.93, "", transform=ax.transAxes, fontsize=10)

ax.legend(loc="upper right", fontsize=9)
ax.axhline(0, color="gray", lw=0.6)   # zero reference line

# ------------------------------------------------------------------
# 5. ANIMATION FUNCTIONS
# ------------------------------------------------------------------

def init_anim():
    """Called once before the animation starts; returns the artists to draw."""
    line_prob.set_data([], [])
    line_re.set_data([], [])
    time_text.set_text("")
    return line_prob, line_re, time_text


def update_frame(frame_idx):
    """Called for each frame; updates the animated artists."""
    f = frames[frame_idx]
    line_prob.set_data(f["x"], f["prob"])
    line_re.set_data(f["x"], f["re_psi"])
    time_text.set_text(f"t = {frame_times[frame_idx]:.4f} a.u.")
    return line_prob, line_re, time_text


# Build the animation object.
# interval = milliseconds between frames in the displayed animation.
anim = animation.FuncAnimation(
    fig,
    update_frame,
    frames=len(frames),
    init_func=init_anim,
    interval=50,       # 50 ms → ~20 fps
    blit=True
)

plt.tight_layout()

# ------------------------------------------------------------------
# 6. SAVE THE ANIMATION
# ------------------------------------------------------------------
# Try MP4 first (requires ffmpeg); fall back to GIF.
mp4_path = os.path.join(output_dir, "evolution.mp4")
gif_path = os.path.join(output_dir, "evolution.gif")

try:
    writer_mp4 = animation.FFMpegWriter(fps=20, bitrate=1200)
    anim.save(mp4_path, writer=writer_mp4)
    print(f"Animation saved to {mp4_path}")
except Exception as e:
    print(f"MP4 failed ({e}), trying GIF...")
    try:
        anim.save(gif_path, writer="pillow", fps=20)
        print(f"Animation saved to {gif_path}")
    except Exception as e2:
        print(f"GIF also failed ({e2}). "
              "Install ffmpeg for MP4 or pillow for GIF.")

plt.show()
