#!/usr/bin/env python3
import csv, io, math, os, sys, urllib.request, datetime as dt
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
names = {f.name for f in fm.fontManager.ttflist}
plt.rcParams.update({"font.family": "Nimbus Sans" if "Nimbus Sans" in names else "sans-serif",
                     "font.size": 8, "axes.labelsize": 9, "axes.linewidth": 0.6, "pdf.fonttype": 42,
                     "xtick.direction": "out", "ytick.direction": "out",
                     "xtick.minor.visible": True, "ytick.minor.visible": True})
W, S, E, N, DMAX = 65.6, 35.85, 69.4, 37.45, 70
URL = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=1973-01-01&endtime=2026-11-10"
       f"&minlatitude={S}&maxlatitude={N}&minlongitude={W}&maxlongitude={E}&minmagnitude=2.5&maxdepth={DMAX}"
       "&orderby=time-asc")
rows = list(csv.DictReader(io.StringIO(urllib.request.urlopen(URL, timeout=120).read().decode())))
t = np.array([dt.datetime.fromisoformat(r["time"].replace("Z", "")) for r in rows])
m = np.round(np.array([float(r["mag"]) for r in rows]), 1)
mtypes = {}
for r in rows: mtypes[r["magType"]] = mtypes.get(r["magType"], 0) + 1
DM = 0.1
T25 = dt.datetime(2025, 11, 2, 20, 29)

def bmle(mm, mc):
    x = mm[mm >= mc - 1e-9]; n = len(x)
    b = math.log10(math.e) / (x.mean() - (mc - DM / 2))
    sb = 2.3 * b * b * math.sqrt(((x - x.mean()) ** 2).sum() / (n * (n - 1)))
    return b, sb, n

bins = np.round(np.arange(m.min(), m.max() + DM / 2, DM), 1)
inc = np.array([(m == b).sum() for b in bins])
cum = np.array([(m >= b - 1e-9).sum() for b in bins])
mc_maxc = bins[inc.argmax()]; MC = round(mc_maxc + 0.2, 1)
b, sb, nmc = bmle(m, MC)
a = math.log10(nmc) + b * MC
rng = np.random.default_rng(42); x = m[m >= MC - 1e-9]
boot = np.array([bmle(rng.choice(x, len(x)), MC)[0] for _ in range(1000)])
blo, bhi = np.percentile(boot, [2.5, 97.5])

mcs = np.round(np.arange(3.5, 5.41, 0.1), 1)
st = [(c, *bmle(m, c)) for c in mcs if (m >= c - 1e-9).sum() >= 20]
print(f"N={len(m)} types={mtypes} MAXC={mc_maxc} Mc={MC} b={b:.3f}+-{sb:.3f} boot95=[{blo:.3f},{bhi:.3f}] a={a:.2f} N>=Mc={nmc} tmin={t.min()} mmax={m.max()}")

CM = 1 / 2.54
fig = plt.figure(figsize=(17 * CM, 17.4 * CM))
TURBO = plt.get_cmap("turbo")

ax = fig.add_axes([0.075, 0.62, 0.40, 0.34])
ax.bar(bins, inc, width=DM * 0.85, color="#c8c8c8", edgecolor="#6e6e6e", lw=0.4, label="Incremental", zorder=2)
NORM = mpl.colors.Normalize(bins.min(), bins.max())
ax.scatter(bins, cum, s=16, c=bins, cmap=TURBO, norm=NORM, edgecolor="black", lw=0.4, zorder=4,
           label="Cumulative")
xx = np.linspace(MC, bins.max() + 0.2, 50)
ax.plot(xx, 10 ** (a - b * xx), color="black", lw=1.2, zorder=3,
        label=f"log$_{{10}}$N = {a:.2f} − {b:.2f}M")
ax.axvline(MC, color="#c0161c", lw=0.9, ls=(0, (4, 2)), zorder=1)
ax.text(MC - 0.04, 0.70, f"M$_c$ = {MC:.1f}", color="#c0161c", fontsize=8, va="center", ha="right", transform=ax.get_xaxis_transform())
ax.set_yscale("log"); ax.set_xlim(2.4, bins.max() + 0.3); ax.set_ylim(0.7, cum.max() * 60)
ax.set_xlabel("Magnitude"); ax.set_ylabel("Number of events")
ax.xaxis.set_major_locator(MultipleLocator(0.5)); ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.grid(which="major", color="#dddddd", lw=0.4, zorder=0)
hh, ll = ax.get_legend_handles_labels(); oo = [ll.index("Incremental"), ll.index("Cumulative")] + [i for i, q in enumerate(ll) if q.startswith("log")]
ax.legend([hh[i] for i in oo], [ll[i] for i in oo], loc="upper left", frameon=False, fontsize=8, handlelength=1.4)
ax.text(0.97, 0.965, f"b = {b:.2f} ± {sb:.2f} (Shi–Bolt)\n95% bootstrap: {blo:.2f}–{bhi:.2f}\n"
        f"N(M ≥ M$_c$) = {nmc}", transform=ax.transAxes, fontsize=8, va="top", ha="right",
        bbox=dict(fc="white", ec="#999999", lw=0.4, boxstyle="round,pad=0.3"))
ax.text(-0.15, 1.04, "(a)", transform=ax.transAxes, fontsize=11, fontweight="bold")

ax2 = fig.add_axes([0.585, 0.62, 0.315, 0.34])
c_, b_, s_, n_ = map(np.array, zip(*st))
ax2.fill_between(c_, b_ - s_, b_ + s_, color="#d9d9d9", lw=0, zorder=1, label="±1σ (Shi–Bolt)")
ax2.plot(c_, b_, color="black", lw=0.9, zorder=2)
sc = ax2.scatter(c_, b_, c=c_, cmap=TURBO, norm=NORM, s=22,
                 edgecolor="black", lw=0.4, zorder=3)
ax2.axvline(MC, color="#c0161c", lw=0.9, ls=(0, (4, 2)))
ax2.axhspan(blo, bhi, color="#c0161c", alpha=0.10, lw=0, label="95% bootstrap at M$_c$")
ax2.set_xlabel("Cut-off magnitude M$_{min}$"); ax2.set_ylabel("b-value")
ax2.xaxis.set_major_locator(MultipleLocator(0.5)); ax2.xaxis.set_minor_locator(MultipleLocator(0.1))
ax2.yaxis.set_minor_locator(AutoMinorLocator(5))
ax2.grid(which="major", color="#dddddd", lw=0.4, zorder=0)
ax2.set_ylim(min(b_ - s_) - 0.05, max(b_ + s_) + 0.35)
ax2.legend(loc="upper left", frameon=True, framealpha=0.9, edgecolor="none", fontsize=8)
cb = fig.colorbar(sc, cax=fig.add_axes([0.915, 0.62, 0.015, 0.34]))
cb.ax.yaxis.set_major_locator(MultipleLocator(0.5)); cb.ax.yaxis.set_minor_locator(MultipleLocator(0.1))
cb.set_label("Magnitude M (a), cut-off M$_{min}$ (b)", fontsize=8); cb.ax.tick_params(labelsize=8)
ax2.text(-0.15, 1.04, "(b)", transform=ax2.transAxes, fontsize=11, fontweight="bold")

ax3 = fig.add_axes([0.075, 0.22, 0.855, 0.30])
yrs = np.array([x.year + (x.timetuple().tm_yday - 1) / 365.25 for x in t])
ax3.step(yrs, np.arange(1, len(yrs) + 1), where="post", color="#888888", lw=0.8, zorder=2, label="All events (M ≥ 2.5)")
sel = m >= MC - 1e-9
ax3.step(yrs[sel], np.arange(1, sel.sum() + 1), where="post", color="black", lw=1.2, zorder=3,
         label=f"M ≥ M$_c$ ({MC:.1f})")
ax3b = ax3.twinx()
o = np.argsort(m)
sc3 = ax3b.scatter(yrs[o], m[o], c=m[o], cmap=TURBO, vmin=2.5, vmax=6.5, s=4 + 2.2 ** (m[o] * 1.1) / 6,
                   edgecolor="black", lw=0.25, alpha=0.9, zorder=4)
t25 = T25.year + (T25.timetuple().tm_yday - 1) / 365.25
ax3.axvline(t25, color="#c0161c", lw=0.9, ls=(0, (4, 2)), zorder=1)
ax3.text(t25 - 0.4, len(yrs) * 1.17, "2025\nM$_{ww}$ 6.2", color="#c0161c", ha="right", va="top", fontsize=8)
t02 = 2002 + 84 / 365.25
ax3.annotate("2002\nM$_{wc}$ 6.1", xy=(t02, len(yrs) * 0.86), xytext=(t02 - 4.5, len(yrs) * 1.08), fontsize=8, color="#333333", ha="center", va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="#333333"))
t82 = 1982 + 349 / 365.25
ax3.annotate("1982\nM$_s$ 6.6", xy=(t82 + 0.3, len(yrs) * 0.95), xytext=(t82 + 7.5, len(yrs) * 1.10), fontsize=8, color="#333333", ha="center", va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="#333333"))
ax3.set_xlim(1972.5, 2027); ax3.set_ylim(0, len(yrs) * 1.2)
ax3b.set_ylim(2.3, 7.8); ax3b.set_ylabel("Magnitude"); ax3b.yaxis.set_minor_locator(MultipleLocator(0.1))
ax3.set_xlabel("Year"); ax3.set_ylabel("Cumulative number of events")
ax3.xaxis.set_major_locator(MultipleLocator(5)); ax3.xaxis.set_minor_locator(MultipleLocator(1))
ax3.yaxis.set_minor_locator(AutoMinorLocator(5))
ax3.grid(which="major", color="#dddddd", lw=0.4, zorder=0)
ax3.legend(loc="upper left", frameon=False, fontsize=8)
ax3b.yaxis.set_label_position("right"); ax3b.yaxis.tick_right()
ax3.text(-0.07, 1.04, "(c)", transform=ax3.transAxes, fontsize=11, fontweight="bold")

ax0 = fig.add_axes([0, 0, 1, 1]); ax0.axis("off"); ax0.set_xlim(0, 17); ax0.set_ylim(0, 17.4)
R = fig.canvas.get_renderer(); PX = fig.dpi * CM
def tw(s):
    q = ax0.text(0, 0, s, fontsize=8); w = q.get_window_extent(R).width / PX; q.remove(); return w
mt = ", ".join(f"{k} {v}" for k, v in sorted(mtypes.items(), key=lambda z: -z[1]))
txt = (f"Depicts: frequency-magnitude statistics of {len(m)} crustal earthquakes (depth ≤ {DMAX} km, M ≥ 2.5, "
       f"1973-2026; magnitude types: {mt}) in the study area {W}-{E}°E, {S}-{N}°N: (a) incremental and "
       f"cumulative distribution with the Gutenberg-Richter fit above the completeness magnitude M$_c$ (maximum "
       f"curvature + 0.2); (b) maximum-likelihood b-value (Aki-Utsu) against cut-off magnitude; (a) and (b) share the "
       f"magnitude colour scale; (c) cumulative number of events and magnitudes with time. Data: USGS ComCat. Software: Python "
       f"{sys.version.split()[0]}, NumPy {np.__version__}, Matplotlib {mpl.__version__}. Source: authors.")
x0, wmax = 0.6, 16.2
words = txt.split(); lines = []; cur = ""
for wd in words:
    q = (cur + " " + wd).strip()
    if tw(q) <= wmax: cur = q
    else: lines.append(cur); cur = wd
lines.append(cur); sp0 = tw("a a") - tw("aa"); y = 2.75
for i, s in enumerate(lines):
    if i < len(lines) - 1:
        ws = s.split(); gap = (wmax - tw(s)) / (len(ws) - 1); x = x0
        for wd in ws:
            ax0.text(x, y, wd, fontsize=8, va="top", color="#333333"); x += tw(wd) + sp0 + gap
    else:
        ax0.text(x0, y, s, fontsize=8, va="top", color="#333333")
    y -= 0.37
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"magnitude_frequency.{ext}"), dpi=300)
