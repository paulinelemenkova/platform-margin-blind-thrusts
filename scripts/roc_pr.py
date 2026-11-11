#!/usr/bin/env python3
import json, os, sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.ticker import MultipleLocator
import sklearn
from sklearn.metrics import roc_curve, precision_recall_curve

IN = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = sys.argv[2] if len(sys.argv) > 2 else "."
names = {f.name for f in fm.fontManager.ttflist}
plt.rcParams.update({"font.family": "Nimbus Sans" if "Nimbus Sans" in names else "sans-serif",
                     "font.size": 8, "axes.labelsize": 9, "axes.linewidth": 0.6, "pdf.fonttype": 42,
                     "xtick.direction": "out", "ytick.direction": "out",
                     "xtick.minor.visible": True, "ytick.minor.visible": True})
oof = pd.read_csv(os.path.join(IN, "oof.csv"))
met = json.load(open(os.path.join(IN, "metrics.json")))
y, fold = oof["y"].values, oof["fold"].values
TURBO = plt.get_cmap("turbo")
MOD = {"RF": dict(c=TURBO(0.10), mk="o", name="Random Forest"),
       "HGB": dict(c=TURBO(0.88), mk="s", name="Gradient boosting")}
prev = y.mean()

CM = 1 / 2.54
FW, FH = 17.0, 9.7
fig = plt.figure(figsize=(FW * CM, FH * CM))
SZ = 4.35
AX = [fig.add_axes([(1.15 + i * 5.55) / FW, 4.0 / FH, SZ / FW, SZ / FH]) for i in range(3)]

def style(ax, lab):
    ax.grid(which="major", color="#dddddd", lw=0.4, zorder=0)
    ax.text(-0.25, 1.07, lab, transform=ax.transAxes, fontsize=11, fontweight="bold")

ax = AX[0]
ax.plot([0, 1], [0, 1], ls="--", lw=0.7, color="#888888", zorder=1)
ax.text(0.62, 0.55, "random", rotation=45, rotation_mode="anchor", fontsize=7, color="#777777",
        ha="center", va="bottom", transform=ax.transData)
for k, d in MOD.items():
    p = oof[f"p_{k}"].values; r = met["models"][k]
    for f in np.unique(fold):
        s = fold == f; fp, tp, _ = roc_curve(y[s], p[s])
        ax.plot(fp, tp, color=d["c"], lw=0.5, alpha=0.45, zorder=2)
    fp, tp, _ = roc_curve(y, p)
    ax.plot(fp, tp, color=d["c"], lw=1.6, zorder=3,
            label=f"{k}: AUC {r['auc_mean']:.2f} ± {r['auc_sd']:.2f}")
    i = np.linspace(0, len(fp) - 1, 9).astype(int)[1:-1]
    ax.plot(fp[i], tp[i], ls="none", marker=d["mk"], ms=3.6, mfc="white", mec=d["c"], mew=0.9, zorder=4)
ax.set(xlim=(0, 1), ylim=(0, 1.0), xlabel="False positive rate", ylabel="True positive rate")
ax.set_aspect("equal")
for a in (ax.xaxis, ax.yaxis):
    a.set_major_locator(MultipleLocator(0.2)); a.set_minor_locator(MultipleLocator(0.05))
ax.legend(loc="lower right", frameon=False, fontsize=7, handlelength=1.6, borderaxespad=0.2)
style(ax, "(a)")

ax = AX[1]
ax.axhline(prev, ls="--", lw=0.7, color="#888888", zorder=1)
ax.text(0.03, prev - 0.015, f"no skill ({prev:.2f})", fontsize=7, color="#777777", va="top")
for k, d in MOD.items():
    p = oof[f"p_{k}"].values; r = met["models"][k]
    for f in np.unique(fold):
        s = fold == f; pr, rc, _ = precision_recall_curve(y[s], p[s])
        ax.plot(rc, pr, color=d["c"], lw=0.5, alpha=0.45, zorder=2)
    pr, rc, _ = precision_recall_curve(y, p)
    ax.plot(rc, pr, color=d["c"], lw=1.6, zorder=3,
            label=f"{k}: AP {r['ap_mean']:.2f} ± {r['ap_sd']:.2f}")
    i = np.linspace(0, len(rc) - 1, 9).astype(int)[1:-1]
    ax.plot(rc[i], pr[i], ls="none", marker=d["mk"], ms=3.6, mfc="white", mec=d["c"], mew=0.9, zorder=4)
ax.set(xlim=(0, 1), ylim=(0, 1.0), xlabel="Recall", ylabel="Precision")
ax.set_aspect("equal")
for a in (ax.xaxis, ax.yaxis):
    a.set_major_locator(MultipleLocator(0.2)); a.set_minor_locator(MultipleLocator(0.05))
ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.08), frameon=False, fontsize=7, handlelength=1.6,
          borderaxespad=0.2)
style(ax, "(b)")

ax = AX[2]
best = met["best"]; p = oof[f"p_{best}"].values
edges = np.linspace(0, 1, 26)
br = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
hn, _ = np.histogram(p[y == 0], edges, density=True); hp, _ = np.histogram(p[y == 1], edges, density=True)
ax.stairs(hn, edges, fill=True, color=TURBO(0.10), alpha=0.35, zorder=2)
ax.stairs(hn, edges, color=TURBO(0.10), lw=1.1, zorder=3, label=f"Negative cells (n = {met['n_neg']})")
ax.stairs(hp, edges, fill=True, color=TURBO(0.88), alpha=0.35, zorder=2, hatch="////")
ax.stairs(hp, edges, color=TURBO(0.88), lw=1.1, zorder=3, label=f"Positive cells (n = {met['n_pos']})")
top = 1.12 * max(hn.max(), hp.max())
ax.set(xlim=(0, 1), ylim=(0, top), xlabel=f"Out-of-fold probability ({best})", ylabel="Density")
ax.xaxis.set_major_locator(MultipleLocator(0.2)); ax.xaxis.set_minor_locator(MultipleLocator(0.05))
for b in br[1:-1]:
    ax.axvline(b, color="#999999", lw=0.4, ls=":", zorder=1)
na = met["n_aftershocks"]; fa = met["aftershock_frac"]; cls = int(np.argmax(fa))
ax.annotate(f"{na} aftershocks (2025):\n{100 * fa[cls]:.0f}% in class '{met['class_names'][cls]}'",
            xy=(0.9, 0.06 * top), xytext=(0.62, 0.62 * top), fontsize=7, ha="center",
            arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#333333", mutation_scale=7))
ax.legend(loc="upper center", frameon=False, fontsize=7, handlelength=1.4, borderaxespad=0.2)
style(ax, "(c)")
sb = fig.add_axes([(1.15 + 2 * 5.55) / FW, (4.0 + SZ + 0.08) / FH, SZ / FW, 0.28 / FH])
for j, nm in enumerate(["VL", "L", "M", "H", "VH"]):
    sb.axvspan(br[j], br[j + 1], color=TURBO((j + 0.5) / 5), lw=0)
    sb.text((br[j] + br[j + 1]) / 2, 0.5, nm, ha="center", va="center", fontsize=6.5,
            color="white" if j in (0, 4) else "black")
sb.set(xlim=(0, 1), ylim=(0, 1)); sb.axis("off")

ax0 = fig.add_axes([0, 0, 1, 1]); ax0.axis("off"); ax0.set_xlim(0, FW); ax0.set_ylim(0, FH)
R = fig.canvas.get_renderer(); PX = fig.dpi * CM

def tw(s):
    q = ax0.text(0, 0, s, fontsize=8); w = q.get_window_extent(R).width / PX; q.remove(); return w

rb, hb = met["models"]["RF"], met["models"]["HGB"]
txt = (f"Depicts: validation of the blind-thrust susceptibility models by spatial block cross-validation "
       f"({met['n_blocks']} blocks of 0.25°, 5 stratified group folds, {met['n_labelled']} labelled cells): "
       f"(a) receiver operating characteristic and (b) precision-recall curves of the out-of-fold predictions "
       f"for Random Forest (RF; pooled AUC {rb['auc_pooled']:.2f}, AP {rb['ap_pooled']:.2f}) and histogram "
       f"gradient boosting (HGB; pooled AUC {hb['auc_pooled']:.2f}, AP {hb['ap_pooled']:.2f}); thick lines: all "
       f"folds pooled, thin lines: individual held-out folds, legend values: fold mean ± standard deviation; "
       f"(c) density of out-of-fold {best} probabilities for positive and negative cells, with the five "
       f"susceptibility classes (VL very low to VH very high, strip above); the 2025 aftershock test rests on only {na} ComCat events. Software: "
       f"Python {sys.version.split()[0]}, scikit-learn {sklearn.__version__}, Matplotlib {mpl.__version__}. "
       f"Source: authors.")
x0, wmax = 0.6, FW - 1.2
words = txt.split(); lines = []; cur = ""
for wd in words:
    q = (cur + " " + wd).strip()
    if tw(q) <= wmax: cur = q
    else: lines.append(cur); cur = wd
lines.append(cur); sp0 = tw("a a") - tw("aa"); yy = 2.75
for i, s in enumerate(lines):
    if i < len(lines) - 1:
        ws = s.split(); gap = (wmax - tw(s)) / (len(ws) - 1); x = x0
        for wd in ws:
            ax0.text(x, yy, wd, fontsize=8, va="top", color="#333333"); x += tw(wd) + sp0 + gap
    else:
        ax0.text(x0, yy, s, fontsize=8, va="top", color="#333333")
    yy -= 0.37
print("credit lines", len(lines), "bottom", round(yy, 2))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"roc_pr.{ext}"), dpi=300)
