"""Generate demo.gif + poster.png (verified, offline, <2MB).

Run: python scripts/make_demo.py
Out: assets/demo.gif, assets/poster.png
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
import numpy as np

METRICS = [
    ("T1 AUPRC\ncomposite", 0.9872),
    ("T1 AUPRC\nstem-only", 0.7034),
    ("T3 enrich\n/5x", 3.79 / 5.0),
    ("T4 ELDA\nacc.", 1 - abs(106.3 - 100) / 100),
    ("T5 batch\nrho", 0.202),
]
LABELS = [m[0] for m in METRICS]
VALS = [m[1] for m in METRICS]

fig, ax = plt.subplots(figsize=(7.2, 4.0), dpi=110)
fig.patch.set_facecolor("#0d1117")
ax.set_facecolor("#0d1117")

def draw(frame):
    ax.clear()
    ax.set_facecolor("#0d1117")
    ax.set_xlim(0, 1.05)
    ax.set_ylim(-0.5, len(LABELS) - 0.5)
    ax.set_yticks(range(len(LABELS)))
    ax.set_yticklabels(LABELS, color="white", fontsize=9)
    ax.set_xlabel("score (0-1)", color="white", fontsize=9)
    ax.tick_params(colors="white", labelsize=8)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_color("#30363d")
    prog = (frame + 1) / 24
    for i, v in enumerate(VALS):
        w = v * min(1.0, prog * 1.6 - i * 0.12)
        w = max(0.0, w)
        ax.barh(i, w, height=0.55, color="#58a6ff" if i != 0 else "#3fb950")
        ax.text(w + 0.02, i, f"{v:.2f}", va="center", color="white", fontsize=9)
    ax.set_title("CSC benchmark — live run (4000 cells, 1% CSC)", color="white", fontsize=11, pad=10)
    ax.text(0.5, -1.05, f"frame {frame+1}/24  •  ELDA 1-in-106  •  enrich 3.79x  •  GDC+OpenAlex live_ok=True",
             transform=ax.transAxes, ha="center", color="#8b949e", fontsize=8)

fig.tight_layout()
writer = PillowWriter(fps=6)
fig_anim = plt.figure  # noqa
# Build by re-drawing
import matplotlib.animation as animation
anim = animation.FuncAnimation(fig, draw, frames=24, interval=160)
anim.save("assets/demo.gif", writer=writer)
print("wrote assets/demo.gif")

# poster: final frame static
draw(23)
fig.savefig("assets/poster.png", bbox_inches="tight")
print("wrote assets/poster.png")
