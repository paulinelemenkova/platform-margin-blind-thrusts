#!/usr/bin/env python3
import math, os, sys
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib import font_manager as fm
from matplotlib.ticker import MultipleLocator

WD = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = sys.argv[2] if len(sys.argv) > 2 else "."
names = {f.name for f in fm.fontManager.ttflist}
plt.rcParams.update({"font.family": "Nimbus Sans" if "Nimbus Sans" in names else "sans-serif",
                     "font.size": 8, "axes.labelsize": 8.5, "axes.linewidth": 0.6, "pdf.fonttype": 42,
                     "xtick.direction": "out", "ytick.direction": "out", "xtick.minor.visible": True,
                     "ytick.minor.visible": True, "hatch.linewidth": 0.4})
WO = [pe.withStroke(linewidth=2.0, foreground="white")]

AZ = 340.0; EQ = (67.3886, 36.5461); EQDEP = 28.0
NP = [(241, 70, 37), (137, 55, 156)]
INS = dict(lon=67.2624, lat=36.6894, depth=7.69, top=3.11, bottom=12.27, strike=270.8, dip=13.24, W=40.0,
           rake=-6.2)
t = np.loadtxt(os.path.join(WD, "topo.txt")); xs, el = t[:, 0], t[:, 1]
k = 9; el_s = np.convolve(np.pad(el, k // 2, mode="edge"), np.ones(k) / k, mode="valid")
AMF_ALL = []

def along(lon, lat):
    kx = 111.32 * math.cos(math.radians(EQ[1])); dx = (lon - EQ[0]) * kx; dy = (lat - EQ[1]) * 111.195
    a = math.radians(AZ); return dx * math.sin(a) + dy * math.cos(a), dx * math.cos(a) - dy * math.sin(a)

def app_dip(strike, dip):
    dd = (strike + 90) % 360; c = math.cos(math.radians(dd - AZ))
    return math.degrees(math.atan(math.tan(math.radians(dip)) * abs(c))), (1 if c > 0 else -1)

segs, cur = [], []
for ln in open(os.path.join(WD, "amf.gmt")):
    if ln.startswith(">"):
        if cur: segs.append(cur)
        cur = []
    elif ln.strip() and not ln.startswith("#"):
        a = ln.split(); cur.append(along(float(a[0]), float(a[1])))
if cur: segs.append(cur)
for s in segs:
    for a_, c_ in s:
        if abs(c_) <= 40 and -70 <= a_ <= 70: AMF_ALL.append(a_)
AMF = float(np.mean(AMF_ALL)) if AMF_ALL else -12.0
AMF_LAB = "Alburz–Marmul Fault (projected)" if AMF_ALL else "range-front thrust (schematic)"
print("AMF projected:", (round(min(AMF_ALL), 1), round(AMF, 1), round(max(AMF_ALL), 1)) if AMF_ALL else None)

X = np.linspace(-70, 70, 1401)
f = (X + 70) / 140
zb = 5.0 + 2.6 * f
z_kug = zb - 1.2
z_gau = z_kug - (0.25 + 0.45 * f)
TH = [("lcr", 1.25), ("ucr", 0.6), ("gho", 0.35), ("eoc", 0.5)]
STR = [dict(c=-18.0, A=0.8, wb=8.0, wf=3.0, tip="ucr_top"),
       dict(c=AMF - 4.5, A=1.7, wb=9.0, wf=3.2, tip=None),
       dict(c=AMF + 17.0, A=0.75, wb=7.0, wf=2.8, tip="ucr_top"),
       dict(c=AMF + 29.0, A=0.4, wb=5.5, wf=2.4, tip="lcr_top")]

def uplift(x):
    u = np.zeros_like(x)
    for s in STR:
        w = np.where(x < s["c"], s["wb"], s["wf"]); u += s["A"] * np.exp(-((x - s["c"]) / w) ** 2)
    return u

U = uplift(X)
H = {"bas_top": zb, "kug_top": z_kug, "gau_top": z_gau - 0.55 * U}
z = z_gau
for key, th in TH:
    z = z - th; H[key + "_top"] = z - U
surf = -np.interp(X, xs, el_s)
UNITS = [("bas", "Palaeozoic–Triassic basement", "#c9c0bb", "xx"),
         ("kug", "Jurassic clastics and carbonates", "#a9d9e8", ".."),
         ("gau", "Gaurdak evaporites (Upper Jurassic)", "#e7c4e8", "++"),
         ("lcr", "Lower Cretaceous (incl. Hauterivian)", "#b8df8f", None),
         ("ucr", "Upper Cretaceous", "#d6eba0", None),
         ("gho", "Ghory limestone (Campanian–Paleocene)", "#f6c58c", None),
         ("eoc", "Eocene–Oligocene", "#fbd3a4", None),
         ("nq", "Neogene–Quaternary molasse", "#fff09a", None)]
TOPS = {"bas": "bas_top", "kug": "kug_top", "gau": "gau_top", "lcr": "lcr_top", "ucr": "ucr_top",
        "gho": "gho_top", "eoc": "eoc_top"}
order = ["bas", "kug", "gau", "lcr", "ucr", "gho", "eoc", "nq"]

def unit_bounds(u, zbot):
    i = order.index(u)
    top = H[TOPS[u]] if u != "nq" else surf
    bot = zbot if u == "bas" else H[TOPS[order[i - 1]]]
    top = np.maximum(top, surf); bot = np.maximum(bot, surf)
    return top, bot

RAMP_DIP = math.radians(30)
ramps = []
z_det = np.interp(AMF, X, z_gau) + 0.15
x0 = AMF - (z_det - np.interp(AMF, X, surf)) / math.tan(RAMP_DIP)
ramps.append(dict(x=[x0, AMF], z=[z_det, np.interp(AMF, X, surf)], blind=False))
for s in STR:
    if s["tip"] is None: continue
    xt = s["c"] + 1.6; zt = float(np.interp(xt, X, H[s["tip"]])) + 0.25
    zd = float(np.interp(xt, X, z_gau)) + 0.15; xr = xt - (zd - zt) / math.tan(RAMP_DIP)
    ramps.append(dict(x=[xr, xt], z=[zd, zt], blind=True))

xc_i, _ = along(INS["lon"], INS["lat"]); dip_i, sg_i = app_dip(INS["strike"], INS["dip"])
hw = 0.5 * INS["W"] * math.cos(math.radians(INS["dip"])) * abs(math.cos(math.radians((INS["strike"] + 90) % 360 - AZ)))
ins_x = [xc_i - sg_i * hw, xc_i + sg_i * hw]; ins_z = [INS["top"], INS["bottom"]]
xm, ym = along(67.110, 36.709)
print(f"InSAR plane x {ins_x[0]:.1f}..{ins_x[1]:.1f} km, apparent dip {dip_i:.1f}; C01 {xm:.1f} km, off {ym:.1f} km")

CM = 1 / 2.54
FW, FH = 17.0, 18.5
fig = plt.figure(figsize=(FW * CM, FH * CM))
L0, WA = 1.35, 15.25
xa = (-70, 70); za = (-3.0, 36.0); sa = WA / (xa[1] - xa[0])
HA = (za[1] - za[0]) * sa
axA = fig.add_axes([L0 / FW, 13.45 / FH, WA / FW, HA / FH])
xb = (-40, 45); zb_lim = (-2.6, 8.4); VE = 3.0; sb = WA / (xb[1] - xb[0])
HB = (zb_lim[1] - zb_lim[0]) * sb * VE
axB = fig.add_axes([L0 / FW, 5.95 / FH, WA / FW, HB / FH])
print("panel heights", round(HA, 2), round(HB, 2))

def draw_units(ax, zbot, hatch=True):
    for u, lab, col, ht in UNITS:
        top, bot = unit_bounds(u, zbot)
        ax.fill_between(X, top, bot, where=bot > top + 1e-6, color=col, lw=0, zorder=1)
        if hatch and ht:
            ax.fill_between(X, top, bot, where=bot > top + 1e-6, facecolor="none", edgecolor="#8a7f78" if u == "bas"
                            else "#7a5a8a" if u == "gau" else "#5d8ea0", hatch=ht, lw=0, zorder=1)
    for key in ("bas_top", "kug_top", "gau_top", "lcr_top", "ucr_top", "gho_top", "eoc_top"):
        zz = np.where(H[key] >= surf, H[key], np.nan)
        ax.plot(X, zz, color="#555555", lw=0.35, zorder=2)
    ax.fill_between(X, -10, surf, color="white", lw=0, zorder=3)
    ax.plot(X, surf, color="black", lw=0.8, zorder=4)

axA.add_patch(mpl.patches.Rectangle((-70, 0), 140, 40, facecolor="#e4e0dc", lw=0, zorder=0))
draw_units(axA, za[1] + 1, hatch=False)
axA.fill_between(X, H["bas_top"], za[1] + 1, color="#e4e0dc", lw=0, zorder=1.5)
axA.fill_between(X, H["bas_top"], za[1] + 1, facecolor="none", edgecolor="#b5aca5", hatch="xx", lw=0, zorder=1.6)
for (s_, d_, r_), ls, nm in zip(NP, ("-", (0, (5, 2.5))), ("NP1", "NP2")):
    ad, sg = app_dip(s_, d_)
    zz = np.array([12.0, 35.0]); xx = sg * (zz - EQDEP) / math.tan(math.radians(ad))
    axA.plot(xx, zz, color="black", lw=1.0, ls=ls, zorder=6)
    xl = sg * (14.5 - EQDEP) / math.tan(math.radians(ad))
    axA.text(xl + (1.2 if sg < 0 else -1.2), 14.5, f"{nm} ({s_}/{d_}/{r_}; apparent dip {ad:.0f}°)", fontsize=6.8,
             ha="left" if sg < 0 else "right", va="center", path_effects=WO, zorder=8)
axA.plot(ins_x, ins_z, color="#c0161c", lw=3.0, alpha=0.55, solid_capstyle="butt", zorder=5)
axA.plot(ins_x, ins_z, color="#c0161c", lw=0.8, zorder=5)
axA.text(ins_x[0] + 0.5, INS["top"] + 4.4, f"InSAR plane ({INS['top']:.0f}–{INS['bottom']:.0f} km)", fontsize=6.8,
         color="#8c0f14", rotation=-math.degrees(math.atan2(ins_z[1] - ins_z[0], ins_x[1] - ins_x[0])),
         rotation_mode="anchor", ha="left", va="top", path_effects=WO, zorder=8)
axA.plot([0], [EQDEP], marker="*", ms=12, mfc="white", mec="black", mew=0.9, zorder=7)
axA.text(2.2, EQDEP + 2.2, "USGS hypocentre 2025 (28 km)", fontsize=6.8, va="top", path_effects=WO, zorder=8)
axA.add_patch(mpl.patches.Rectangle((xb[0], zb_lim[0]), xb[1] - xb[0], zb_lim[1] - zb_lim[0], fill=False,
                                    ec="black", lw=0.6, ls=(0, (3, 2)), zorder=9))
axA.text(xb[0] + 0.8, zb_lim[0] + 0.3, "(b)", fontsize=7, ha="left", va="top", path_effects=WO, zorder=9)
axA.text(-50, 22, "Deformed basement", fontsize=7.5, style="italic", color="#4d4540", ha="center",
         path_effects=WO, zorder=8)
axA.text(-68.5, -3.6, "SSE", fontsize=8, fontweight="bold", ha="left", va="bottom")
axA.text(68.5, -3.6, "NNW", fontsize=8, fontweight="bold", ha="right", va="bottom")
axA.set(xlim=xa, ylim=za[::-1], ylabel="Depth (km)")
axA.set_aspect(1.0)
axA.xaxis.set_major_locator(MultipleLocator(10)); axA.xaxis.set_minor_locator(MultipleLocator(5))
axA.yaxis.set_major_locator(MultipleLocator(10)); axA.yaxis.set_minor_locator(MultipleLocator(5))
axA.set_xlabel("Distance along P–P′ from the 2025 epicentre (km)", labelpad=2)
axA.text(-0.075, 1.06, "(a)", transform=axA.transAxes, fontsize=11, fontweight="bold")
axA.text(0.99, 0.04, "true scale", transform=axA.transAxes, fontsize=7, ha="right", va="bottom", color="#444444",
         path_effects=WO, zorder=9)

draw_units(axB, zb_lim[1] + 1)
for r in ramps:
    axB.plot(r["x"], r["z"], color="black", lw=1.0, zorder=6)
    if r["blind"]:
        axB.plot(r["x"][1], r["z"][1], marker="o", ms=2.6, mfc="white", mec="black", mew=0.7, zorder=7)
        if r is ramps[2]:
            axB.annotate("blind tip", xy=(r["x"][1], r["z"][1]), xytext=(r["x"][1] + 3.5, r["z"][1] + 1.25),
                         fontsize=6.8, style="italic", path_effects=WO, zorder=8,
                         arrowprops=dict(arrowstyle="-", lw=0.5, color="#333333"))

xd = X[(X >= min(r["x"][0] for r in ramps)) & (X <= max(r["x"][0] for r in ramps))]
axB.plot(xd, np.interp(xd, X, z_gau) + 0.15, color="#c0161c", lw=1.0, ls=(0, (4, 2)), zorder=6)
axB.plot(ins_x, ins_z, color="#c0161c", lw=3.0, alpha=0.55, solid_capstyle="butt", zorder=5)
axB.plot(ins_x, ins_z, color="#c0161c", lw=0.8, zorder=5)
xi = ins_x[0] + (5.6 - INS["top"]) * (ins_x[1] - ins_x[0]) / (ins_z[1] - ins_z[0])
axB.annotate("InSAR fault plane", xy=(xi, 5.6), xytext=(xi + 9, 7.0), fontsize=6.8, color="#8c0f14", path_effects=WO, zorder=8,
             arrowprops=dict(arrowstyle="-", lw=0.5, color="#8c0f14"))
xg = -32.0
axB.text(xg, float(np.interp(xg, X, z_gau)) - 0.05, "Gaurdak detachment", fontsize=7, style="italic",
         color="#7a1f6e", ha="center", va="bottom", path_effects=WO, zorder=8)
axB.text(STR[2]["c"] - 1, float(np.interp(STR[2]["c"], X, H["gho_top"])) - 0.25, "fault-propagation fold",
         fontsize=6.8, style="italic", ha="center", va="bottom", path_effects=WO, zorder=8)
sA = float(np.interp(AMF, X, surf))
axB.annotate(AMF_LAB, xy=(AMF, sA + 0.05), xytext=(AMF + 6.5, 0.75), fontsize=7, style="italic",
             ha="center", va="center", path_effects=WO, zorder=8,
             arrowprops=dict(arrowstyle="-|>", lw=0.5, color="#333333", mutation_scale=6))
axB.annotate("epicentre (2025)", xy=(0, float(np.interp(0, X, surf))), xytext=(5, -2.2), fontsize=7,
             ha="center", va="bottom", path_effects=WO, zorder=8,
             arrowprops=dict(arrowstyle="-|>", lw=0.5, color="#333333", mutation_scale=6))
axB.annotate(f"C01 (projected, {abs(ym):.0f} km W)", xy=(xm, float(np.interp(xm, X, surf))),
             xytext=(xm + 4, -1.3), fontsize=7, ha="center", va="bottom", path_effects=WO, zorder=8,
             arrowprops=dict(arrowstyle="-|>", lw=0.5, color="#333333", mutation_scale=6))
axB.text(37, -0.62, "Amu Darya plain", fontsize=7, style="italic", ha="center", va="bottom", path_effects=WO,
         zorder=8)
axB.text(-26, 6.85, "Basement", fontsize=7, style="italic", color="#4d4540", path_effects=WO, zorder=8)
axB.set(xlim=xb, ylim=zb_lim[::-1], ylabel="Depth (km)")
axB.set_aspect(VE)
axB.xaxis.set_major_locator(MultipleLocator(10)); axB.xaxis.set_minor_locator(MultipleLocator(5))
axB.yaxis.set_major_locator(MultipleLocator(2)); axB.yaxis.set_minor_locator(MultipleLocator(1))
axB.set_xlabel("Distance along P–P′ from the 2025 epicentre (km)", labelpad=2)
axB.text(-0.075, 1.03, "(b)", transform=axB.transAxes, fontsize=11, fontweight="bold")
axB.text(1.0, 1.01, "vertical exaggeration 3×", transform=axB.transAxes, fontsize=7, ha="right", va="bottom",
         color="#444444")

hu = [mpl.patches.Patch(facecolor=c, edgecolor="#777777", lw=0.4, hatch=h, label=l) for _, l, c, h in UNITS[::-1]]
hs = [mpl.lines.Line2D([], [], color="black", lw=1.0, label="Thrust ramp; ○ blind tip"),
      mpl.lines.Line2D([], [], color="#c0161c", lw=1.0, ls=(0, (4, 2)), label="Basal detachment"),
      mpl.lines.Line2D([], [], color="#c0161c", lw=3.0, alpha=0.55, label="InSAR fault plane"),
      mpl.lines.Line2D([], [], color="black", lw=1.0, label="NP1 (apparent dip)"),
      mpl.lines.Line2D([], [], color="black", lw=1.0, ls=(0, (5, 2.5)), label="NP2 (apparent dip)"),
      mpl.lines.Line2D([], [], ls="none", marker="*", ms=9, mfc="white", mec="black", mew=0.8,
                       label="2025 hypocentre (USGS)")]
fig.legend(handles=hu, loc="upper left", bbox_to_anchor=(L0 / FW, 5.0 / FH), ncol=2, frameon=False, fontsize=6.8,
           handlelength=1.6, handleheight=0.9, columnspacing=1.0, labelspacing=0.25)
fig.legend(handles=hs, loc="upper left", bbox_to_anchor=(12.3 / FW, 5.0 / FH), ncol=1, frameon=False,
           fontsize=6.8, handlelength=1.8, labelspacing=0.2)

ax0 = fig.add_axes([0, 0, 1, 1]); ax0.axis("off"); ax0.set_xlim(0, FW); ax0.set_ylim(0, FH)
R = fig.canvas.get_renderer(); PX = fig.dpi * CM

def tw(s):
    q = ax0.text(0, 0, s, fontsize=8); w = q.get_window_extent(R).width / PX; q.remove(); return w

txt = ("Depicts: conceptual cross-section of the Northern Platform margin along P–P′ (azimuth 340° through the "
       "2025 epicentre): (a) crust at true scale with the USGS hypocentre, the apparent dips of both nodal "
       f"planes and the best-fitting InSAR fault plane (dip {INS['dip']:.0f}°, mainly left-lateral); (b) the "
       "sedimentary cover at 3× vertical exaggeration, with north-verging fault-propagation folds above "
       "south-dipping ramps rooted in the Gaurdak evaporites, emergent at the range front (Alburz–Marmul Fault, projected onto the section) and blind beneath "
       "the plain. Topography from SRTM15+; horizon depths and thicknesses are schematic, scaled after USGS "
       "stratigraphy and isopachs. Data: USGS ComCat; USGS OFR 2007-1103; COMET-LiCSAR; SRTM15+. Software: GMT, "
       f"Python {sys.version.split()[0]}, Matplotlib {mpl.__version__}. Source: authors.")
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
    fig.savefig(os.path.join(OUT, f"conceptual_model.{ext}"), dpi=300)
