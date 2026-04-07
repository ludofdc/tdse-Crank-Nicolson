# Setup guide

## Requirements

### C++ solver

- C++17-compatible compiler (`g++` or `clang++`)
- Standard library only — no external dependencies

### Python analysis

- Python 3.8+
- `numpy` and `matplotlib`

```bash
pip install numpy matplotlib
```

- For MP4 animation output: [ffmpeg](https://ffmpeg.org) installed on the system
- For GIF animation output (fallback if ffmpeg is absent): `pip install pillow`

---

## Step 1 — Configure the simulation

Edit `input/input.txt` to set the physical parameters:

```
# Domain and grid
L       20.0    # [bohr]      total domain length
dx       0.01   # [bohr]      spatial step
dt       0.001  # [hbar/Eh]   time step
T_tot    2.0    # [hbar/Eh]   total simulation time

# Gaussian wave packet
x0      -5.0    # [bohr]      initial position of packet centre
sigma    0.5    # [bohr]      packet width
k0       5.0    # [bohr^-1]   initial momentum (positive = rightward)

# Potential barrier
x_bar    0.0    # [bohr]      left edge of barrier
w_bar    0.5    # [bohr]      barrier width
V0      20.0    # [Eh]        barrier height
```

The mean kinetic energy of the packet is E = k₀²/2.
- E < V₀ → tunnelling regime (exponentially suppressed transmission)
- E > V₀ → above-barrier regime (partial reflection still occurs due to wave nature)

---

## Step 2 — Build the solver

From the project root:

```bash
make
```

To rebuild from scratch:

```bash
make clean && make
```

---

## Step 3 — Run the solver

```bash
./tdse
```

The program prints a parameter summary and then a line for every saved snapshot:

```
========================================
  TDSE Crank-Nicolson solver
  Atomic units: hbar=1, m=1
========================================
...
step      0  t = 0.0000  norm^2 = 1.00000000e+00
step    200  t = 0.2000  norm^2 = 9.99999999e-01
...
Done. Snapshots written to output/
```

The `norm²` column should remain = 1 throughout.
Any significant drift (> 1e-6) indicates a resolution problem — try reducing `dt` or `dx`.

---

## Step 4 — Run the analysis scripts

All three scripts must be run **from the project root** (not from inside `analysis/`):

```bash
# Static multi-panel snapshot figure
python3 analysis/plot_snapshots.py       # → output/snapshots.png

# Quantitative analysis: norm, T/R, <x>, sigma
python3 analysis/plot_analysis.py        # → output/analysis.png

# Animated full evolution
python3 analysis/plot_animation.py       # → output/evolution.mp4  (or .gif)
```

---

## Output files

| File | Contents |
|---|---|
| `output/psi_XXXXXX.csv` | Columns: `x`, `\|ψ\|²`, `Re(ψ)`, `Im(ψ)`, `V(x)` |
| `output/norm.csv` | Columns: `t`, `norm²` — one row per snapshot |
| `output/snapshots.png` | Static 5-panel evolution figure |
| `output/analysis.png` | Norm conservation, T/R coefficients, ⟨x⟩(t), σ(t) |
| `output/evolution.mp4` | Animated wave packet (requires ffmpeg) |
| `output/evolution.gif` | Animated wave packet fallback (requires pillow) |

---

## Physical parameters: guidance

| Parameter | Effect |
|---|---|
| Increase `k0` | Higher energy → more transmission, less tunnelling |
| Increase `V0` | Higher barrier → less transmission |
| Increase `w_bar` | Wider barrier → exponentially less tunnelling |
| Increase `sigma` | Broader packet → narrower momentum spread (less definite energy) |
| Decrease `dx`, `dt` | Better accuracy, slower run |
| Increase `L` | Larger domain → more room before boundary reflections |

A good starting point for observing tunnelling is `k0 = 5`, `V0 = 20`, `w_bar = 0.5` (default input), where E = k₀²/2 = 12.5 Eh < V₀ = 20 Eh.
