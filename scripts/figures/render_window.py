import argparse
from pathlib import Path
import csv
parser=argparse.ArgumentParser(description="Render Supplementary Figure S2 from locally prepared CGMacros episode CSVs.")
parser.add_argument("input_dir",nargs="?",default="results/application")
parser.add_argument("output_dir",nargs="?",default="results/figures")
args=parser.parse_args()
input_dir=Path(args.input_dir)
output_dir=Path(args.output_dir)
output_dir.mkdir(parents=True,exist_ok=True)
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#1f1f1e", "#55544f", "#8a8984", "#e6e5df"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2})

def read(fn):
    with (input_dir / fn).open() as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}

series = read("schematic_example_series.csv"); native = read("schematic_example_native.csv")
win = read("schematic_example_window.csv"); meta = {k: v[0] for k, v in read("schematic_example_meta.csv").items()}
minute_t, minute_y = series["t_min"].astype(float), series["gl"].astype(float)
native_t, native_y = native["t_min"].astype(float), native["gl"].astype(float)
slots, excess = win["slot_min"].astype(float), win["excess"].astype(float)
flag = win["flag"]; filled = flag == "gapfill"; coll = flag == "collinear"
peak_t, base = float(meta["peak_minute"]), float(meta["baseline_gl"])
win_start, win_end = peak_t, peak_t + 65
peak_y = native_y[np.argmin(np.abs(native_t - peak_t))]

fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], width_ratios=[1.35, 1], hspace=0.42, wspace=0.28)
axA = fig.add_subplot(gs[0, :]); axB = fig.add_subplot(gs[1, 0]); axC = fig.add_subplot(gs[1, 1])

ax = axA
ax.axvspan(0, 120, color=BLUE, alpha=0.07, lw=0)
ax.axvspan(win_start, win_end, color=AQUA, alpha=0.16, lw=0)
ax.plot(minute_t, minute_y, color=MUTED, lw=1.0, label="per-minute file values (linearly interpolated)")
nk = native["kink"] == "TRUE"
ax.plot(native_t[nk], native_y[nk], "o", ms=3.6, color=BLUE, label="grid reading, knot (slope change)")
ax.plot(native_t[~nk], native_y[~nk], "o", ms=4.2, mfc="white", mec=BLUE, mew=1.2, label="grid position, collinear with neighbours")
ax.axhline(base, color=INK2, lw=1, ls=(0, (2, 2)))
ax.axvline(0, color=INK, lw=1.2)
ax.plot([peak_t], [peak_y + 1.5], marker="*", ms=12, color=ORANGE, mec="white", mew=0.8, zorder=5)
ymax = max(minute_y.max(), native_y.max()); ymin = min(minute_y.min(), native_y.min())
top = ymax + 0.30 * (ymax - ymin); bot = ymin - 0.14 * (ymax - ymin)
ax.text(1.5, top - 0.10 * (top - bot), "meal logged\n(onset, t = 0)", color=INK, va="top", fontsize=8)
ax.text(28, bot + 0.02 * (top - bot), "peak-search window: 0–120 min after onset", color=BLUE, fontsize=7.4, ha="center")
ax.text(peak_t + 3, peak_y + 0.05 * (top - bot), "glycemic peak\n(window anchor)", color=ORANGE, fontsize=7.8, va="bottom")
ax.text(win_end + 2, bot + 0.02 * (top - bot), "recovery window: 65 min from peak", color="#127a55", fontsize=7.4, ha="left", va="bottom")
ax.text(200, base + 0.01 * (top - bot), "pre-meal glucose at onset (centering level)", color=INK2, fontsize=7.5, va="bottom", ha="right")
yb = top - 0.045 * (top - bot)
ax.annotate("", xy=(win_end, yb), xytext=(0, yb), arrowprops=dict(arrowstyle="|-|", color=INK2, lw=0.9, mutation_scale=3))
ax.text(win_end / 2, yb + 0.012 * (top - bot), "isolation: no other meal logged from onset to end of window", color=INK2, fontsize=7.3, ha="center", va="bottom")
ax.set_xlim(-32, 202); ax.set_ylim(bot, top)
ax.set_xlabel("Minutes from logged meal onset"); ax.set_ylabel("Glucose (mg/dL)")
ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(loc="upper right", frameon=False, fontsize=7.5, bbox_to_anchor=(1.0, 0.72))
ax.set_title(f"(a) One meal event ({meta['subject']}, {meta['meal_type'].lower()}): anchor the window at the peak, not at meal onset",
             loc="left", fontsize=9, color=INK)

ax = axB
ax.axhline(0, color=INK2, lw=1, ls=(0, (2, 2)))
ax.plot(slots, excess, color=AQUA, lw=1.4, zorder=1)
ax.plot(slots[~coll], excess[~coll], "o", ms=5, color=BLUE, mec="white", mew=0.8, zorder=3, label="knot")
ax.plot(slots[coll], excess[coll], "o", ms=5.5, mfc="white", mec=BLUE, mew=1.5, zorder=3, label="collinear")
ax.text(1, 1.5, r"0 = pre-meal baseline ($\mu_i = 0$)", color=INK2, fontsize=7.5, ha="left", va="bottom")
ax.set_xticks(slots[::2]); ax.set_xlim(-3, 68); ax.set_ylim(min(-6, excess.min() - 5), excess.max() * 1.15)
ax.set_xlabel("Minutes after peak (14 slots, 5 min apart)"); ax.set_ylabel("Excess glucose (mg/dL)")
ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(loc="upper right", frameon=False, fontsize=7.2)
ax.set_title(r"(b) Model input: one unit $\mathbf{Y}_i=(Y_{i1},\ldots,Y_{i,14})$", loc="left", fontsize=9, color=INK)

ax = axC; ax.set_axis_off(); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_title("(c) From meal events to model units", loc="left", fontsize=9, color=INK)
boxes = [
    (0.92, 0.12, "1706 logged meals, 45 participants"),
    (0.575, 0.42, "Complete 2-h search coverage (−63)\nPeak strictly inside the search\n(−86 at onset, −178 at end)\nNo other meal before window end (−204)\nAll slots on one identified grid (−1)"),
    (0.14, 0.26, "1174 recovery windows,\n14 time points each, 44 participants\nhealthy 433 · prediab. 433 · T2D 308"),
]
for y, h, txt in boxes:
    last = y < 0.3
    ax.add_patch(FancyBboxPatch((0.0, y - h / 2), 1.0, h, boxstyle="round,pad=0.01,rounding_size=0.02",
                                fc="#eef6f2" if last else "#f4f3ee", ec=AQUA if last else MUTED, lw=1.0))
    ax.text(0.5, y, txt, ha="center", va="center", fontsize=7.0, color=INK, linespacing=1.4)
for y0, y1 in [(0.85, 0.80), (0.355, 0.285)]:
    ax.add_patch(FancyArrowPatch((0.5, y0), (0.5, y1), arrowstyle="-|>", mutation_scale=10, color=INK2, lw=1))

fig.savefig(str(output_dir / "figure_S2_window_construction.pdf"), bbox_inches="tight")
fig.savefig(str(output_dir / "figure_S2_window_construction.png"), dpi=200, bbox_inches="tight")
print("saved", meta)
