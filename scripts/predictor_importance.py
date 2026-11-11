#!/usr/bin/env python3
import os, sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.ticker import MultipleLocator
import sklearn

IN = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = sys.argv[2] if len(sys.argv) > 2 else "."
LWMAX = 1.0
names = {f.name for f in fm.fontManager.ttflist}
plt.rcParams.update({"font.family": "Nimbus Sans" if "Nimbus Sans" in names else "sans-serif",
                     "font.size": 8, "axes.labelsize": 9, "axes.linewidth": 0.6, "pdf.fonttype": 42,
                     "xtick.direction": "out", "ytick.direction": "out",
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
                     "grid.linewidth": 0.4, "patch.linewidth": 0.5, "lines.linewidth": LWMAX,
                     "hatch.linewidth": 0.5})
imp = pd.read_csv(os.path.join(IN, "importance.csv"))
LAB = {"elev": "Elevation", "slope": "Slope", "relief": "Local relief (5 km)", "ksn": "k$_{sn}$",
       "hi": "Hypsometric integral", "ghory": "Ghory thickness", "ghory_grad": "Ghory thickness gradient",
       "seis": "Seismicity density"}
TURBO = plt.get_cmap("turbo")
MOD = {"RF": dict(c=TURBO(0.10), mk="o", hatch=None), "HGB": dict(c=TURBO(0.88), mk="s", hatch="////")}
g = imp.groupby(["model", "feature"]).drop_auc
mean, sd = g.mean().unstack(0), g.std().unstack(0)
order = (mean["RF"] + mean["HGB"]).sort_values().index.tolist()
nf = len(order); folds = sorted(imp.fold.unique())

CM = 1 / 2.54
FW, FH = 17.0, 10.85
fig = plt.figure(figsize=(FW * CM, FH * CM))

ax = fig.add_axes([4.0 / FW, 4.0 / FH, 5.6 / FW, 6.2 / FH])
h = 0.36; ypos = np.arange(nf)
for j, (k, d) in enumerate(MOD.items()):
    yy = ypos + (j - 0.5) * h * 1.08
    ax.barh(yy, mean[k].reindex(order), height=h, color=d["c"], alpha=0.30, edgecolor=d["c"], lw=0.6,
            hatch=d["hatch"], zorder=2, label=f"{k} (mean ± 1 SD)")
    ax.errorbar(mean[k].reindex(order), yy, xerr=sd[k].reindex(order), fmt="none", ecolor="#222222",
                elinewidth=0.6, capsize=1.6, capthick=0.6, zorder=4)
    for i, f in enumerate(order):
        v = imp[(imp.model == k) & (imp.feature == f)].drop_auc.values
        jit = np.linspace(-0.10, 0.10, len(v))
        ax.plot(v, yy[i] + jit, ls="none", marker=d["mk"], ms=2.6, mfc="white", mec=d["c"], mew=0.6, zorder=5)
ax.axvline(0, color="#555555", lw=0.6, zorder=3)
ax.set_yticks(ypos, [LAB[f] for f in order])
ax.set_ylim(-0.6, nf - 0.4)
ax.set_xlabel("Drop in ROC AUC after permutation")
ax.xaxis.set_major_locator(MultipleLocator(0.05)); ax.xaxis.set_minor_locator(MultipleLocator(0.01))
ax.set_xlim(-0.02, 0.21)
ax.grid(axis="x", which="major", color="#dddddd", zorder=0)
ax.tick_params(axis="y", which="minor", left=False)
hd = [mpl.patches.Patch(facecolor=mpl.colors.to_rgba(d["c"], 0.30), edgecolor=d["c"], lw=0.6, hatch=d["hatch"],
                        label=f"{k}: mean ± 1 SD") for k, d in MOD.items()]
hd += [mpl.lines.Line2D([], [], ls="none", marker=d["mk"], ms=3, mfc="white", mec=d["c"], mew=0.6,
                        label=f"{k}: fold value") for k, d in MOD.items()]
ax.legend(handles=hd, loc="lower right", frameon=False, fontsize=7, handlelength=1.3, borderaxespad=0.3)
ax.text(-0.68, 1.03, "(a)", transform=ax.transAxes, fontsize=11, fontweight="bold")

ax2 = fig.add_axes([10.35 / FW, 4.0 / FH, 4.8 / FW, 6.2 / FH])
cols = [(k, f) for k in MOD for f in folds]
M = np.array([[imp[(imp.model == k) & (imp.fold == f) & (imp.feature == ft)].drop_auc.values[0]
               for (k, f) in cols] for ft in order]) * 100
vmax = float(np.ceil(M.max() / 5) * 5)
SPEC = plt.get_cmap("Spectral_r")
NORM = mpl.colors.Normalize(vmin=0, vmax=vmax)
im = ax2.imshow(M, origin="lower", aspect="auto", cmap=SPEC, norm=NORM, interpolation="nearest")
for i in range(nf):
    for j in range(len(cols)):
        v = M[i, j]; c = SPEC(np.clip(NORM(v), 0, 1))
        lum = 0.299 * c[0] + 0.587 * c[1] + 0.114 * c[2]
        ax2.text(j, i, f"{round(v) + 0:.0f}", ha="center", va="center", fontsize=6.3,
                 color="white" if lum < 0.5 else "black")
ax2.axvline(len(folds) - 0.5, color="white", lw=LWMAX)
ax2.set_xticks(range(len(cols)), [str(f + 1) for (_, f) in cols])
ax2.set_yticks(range(nf), [""] * nf)
ax2.tick_params(which="minor", bottom=False, left=False)
ax2.set_xlabel("Cross-validation fold")
for j, k in enumerate(MOD):
    ax2.text((j + 0.5) / 2, 1.015, k, transform=ax2.transAxes, ha="center", va="bottom", fontsize=8)
for s in ax2.spines.values(): s.set_linewidth(0.6)
cax = fig.add_axes([15.35 / FW, 4.0 / FH, 0.25 / FW, 6.2 / FH])
cb = fig.colorbar(im, cax=cax, extend="min")
cb.outline.set_linewidth(0.5); cb.ax.tick_params(width=0.5, labelsize=7)
cb.set_label("Drop in ROC AUC (× 100)", fontsize=8)
ax2.text(-0.10, 1.03, "(b)", transform=ax2.transAxes, fontsize=11, fontweight="bold")

ax0 = fig.add_axes([0, 0, 1, 1]); ax0.axis("off"); ax0.set_xlim(0, FW); ax0.set_ylim(0, FH)
R = fig.canvas.get_renderer(); PX = fig.dpi * CM

def tw(s):
    q = ax0.text(0, 0, s, fontsize=8); w = q.get_window_extent(R).width / PX; q.remove(); return w

top = mean["RF"].sort_values(ascending=False)
def lc(k): s = LAB[k]; return s if s.startswith("Ghory") else s[0].lower() + s[1:]
txt = (f"Depicts: permutation importance of the eight predictors of the blind-thrust susceptibility models, "
       f"measured on the held-out folds of the spatial block cross-validation (5 folds of 0.25° blocks, 3 "
       f"permutation repeats per fold) as the drop in ROC AUC: (a) mean and standard deviation across folds for "
       f"Random Forest (RF) and histogram gradient boosting (HGB), with the individual fold values, predictors "
       f"ordered by their combined importance; (b) drop in AUC (× 100) for every fold and predictor. The three most "
       f"important RF predictors are {lc(top.index[0])}, {lc(top.index[1])} and {lc(top.index[2])}; values near zero or below mean no measurable contribution. Software: Python "
       f"{sys.version.split()[0]}, scikit-learn {sklearn.__version__}, Matplotlib {mpl.__version__}. "
       f"Source: authors.")
x0, wmax = 0.6, FW - 1.2
words = txt.split(); lines = []; cur = ""
for wd in words:
    q = (cur + " " + wd).strip()
    if tw(q) <= wmax: cur = q
    else: lines.append(cur); cur = wd
lines.append(cur); sp0 = tw("a a") - tw("aa"); yy = 2.70
for i, s in enumerate(lines):
    if i < len(lines) - 1:
        ws = s.split(); gap = (wmax - tw(s)) / (len(ws) - 1); x = x0
        for wd in ws:
            ax0.text(x, yy, wd, fontsize=8, va="top", color="#333333"); x += tw(wd) + sp0 + gap
    else:
        ax0.text(x0, yy, s, fontsize=8, va="top", color="#333333")
    yy -= 0.37

for a in fig.axes:
    for ln in a.get_lines():
        ln.set_linewidth(min(ln.get_linewidth(), LWMAX))
    for c in a.collections:
        lw = np.atleast_1d(c.get_linewidth())
        if lw.size and lw.max() > LWMAX: c.set_linewidth(np.minimum(lw, LWMAX))
    for p in a.patches:
        p.set_linewidth(min(p.get_linewidth(), LWMAX))
print("credit lines", len(lines), "bottom", round(yy, 2), "order", order[::-1])
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"predictor_importance.{ext}"), dpi=300)
