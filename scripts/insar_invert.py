#!/usr/bin/env python3
import sys, os, json, math
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from okada85 import okada85, corners

wd = sys.argv[1]
LON0, LAT0 = 67.3886, 36.5461
KX = 111.32 * math.cos(math.radians(LAT0)); KY = 111.195
BLK = 8
MU = 3.0e10

def km(lon, lat): return (lon - LON0) * KX, (lat - LAT0) * KY

T = {}
for trk in ("asc", "desc"):
    d = np.load(os.path.join(wd, f"{trk}.npz"))
    lon, lat, los, e, n, u, fr = d["lon"], d["lat"], d["los"], d["e"], d["n"], d["u"], d["frame"]
    hg = d["h"] / 1000.0
    fr = np.where(fr >= 0, 0, -1)
    LO, LA = np.meshgrid(lon, lat)
    ny, nx = los.shape; ny2, nx2 = ny // BLK, nx // BLK
    def blk(a, f=np.nanmean):
        a = a[:ny2 * BLK, :nx2 * BLK].reshape(ny2, BLK, nx2, BLK)
        with np.errstate(all="ignore"): return f(f(a, axis=3), axis=1)
    cnt = blk(np.isfinite(los).astype(float), np.mean)
    with np.errstate(all="ignore"):
        bl = blk(los); be, bn, bu = blk(e), blk(n), blk(u); bh = blk(hg)
        bfr = np.round(blk(np.where(np.isfinite(los), fr, np.nan)))
        blo, bla = blk(LO), blk(LA)
    ok = (cnt >= 0.5) & np.isfinite(bl) & np.isfinite(be) & np.isfinite(bfr) & np.isfinite(bh)
    x, y = km(blo[ok], bla[ok])
    T[trk] = dict(x=x, y=y, d=bl[ok], e=be[ok], n=bn[ok], u=bu[ok], h=bh[ok], fr=bfr[ok].astype(int),
                  full=dict(LO=LO, LA=LA, los=los, e=e, n=n, u=u, h=np.nan_to_num(hg), fr=fr))
    print(trk, "points:", ok.sum(), "frames:", np.unique(bfr[ok]))

frames = [(t, f) for t in T for f in np.unique(T[t]["fr"])]
NR = 4

def design(g, trk_data=None, full=False):
    xc, yc, zc, st, dp, L, W = g
    blocks, rows = [], []
    for t in T:
        D = T[t]
        x, y, e, n, u, fr, h = D["x"], D["y"], D["e"], D["n"], D["u"], D["fr"], D["h"]
        cols = []
        for comp in ("ss", "ds"):
            kw = dict(ss=1.0, ds=0.0) if comp == "ss" else dict(ss=0.0, ds=1.0)
            uE, uN, uZ = okada85(x - xc, y - yc, zc, st, dp, L, W, **kw)
            cols.append(uE * e + uN * n + uZ * u)
        R = np.zeros((x.size, NR * len(frames)))
        for k, (tt, f) in enumerate(frames):
            if tt != t: continue
            m = fr == f; R[m, NR * k] = 1; R[m, NR * k + 1] = x[m] / 50; R[m, NR * k + 2] = y[m] / 50
            R[m, NR * k + 3] = h[m]
        blocks.append(np.column_stack(cols + [R])); rows.append(D["d"])
    return blocks, rows

WT = {t: 1 / math.sqrt(T[t]["x"].size) for t in T}

def solve_lin(g):
    blocks, rows = design(g)
    G = np.vstack([b * WT[t] for b, t in zip(blocks, T)])
    dvec = np.concatenate([r * WT[t] for r, t in zip(rows, T)])

    m, *_ = np.linalg.lstsq(G, dvec, rcond=None)
    return m, G, dvec

def resid(g):
    xc, yc, zc, st, dp, L, W = g
    top = zc - math.sin(math.radians(dp)) * W / 2
    m, G, dvec = solve_lin(g)
    r = G @ m - dvec
    pen = 0.0 if top >= 0.2 else (0.2 - top) * 50.0
    return np.append(r, pen)

lb = [-40, -40, 1.0, -720, 10, 3, 3]; ub = [40, 40, 40, 720, 89.5, 70, 40]
starts = []
for st, dp in ((241, 70), (61, 70), (137, 55), (317, 55), (241, 45), (61, 30), (250, 60), (70, 60)):
    for (xc, yc) in ((0, 4),):
        starts.append([xc, yc, 12, st, dp, 25, 12])
best = None
for s0 in starts:
    try:
        r = least_squares(resid, s0, bounds=(lb, ub), x_scale=[5, 5, 5, 20, 10, 5, 5], max_nfev=400,
                          diff_step=1e-3)
    except Exception as ex:
        print("fail", s0, ex); continue
    c = float(np.sum(r.fun ** 2))
    print(f"start st={s0[3]} dp={s0[4]} xy=({s0[0]},{s0[1]}) -> cost {c:.4f} g={np.round(r.x, 2)}")
    if best is None or c < best[0]: best = (c, r)
cost, r = best

def wrms(g):
    m_, G_, d_ = solve_lin(g); rr = G_ @ m_ - d_
    return float(np.sqrt(np.sum(rr ** 2) / len(T))), m_
tests = {}
for name, (fst, fdp) in (("NP1", (241, 70)), ("NP2", (137, 55))):
    def rfix(h):
        return resid([h[0], h[1], h[2], fst, fdp, h[3], h[4]])
    bt = None
    for h0 in ([0, 4, 15, 20, 12], [0, 0, 25, 15, 10]):
        rr_ = least_squares(rfix, h0, bounds=([-40, -40, 1, 3, 3], [40, 40, 40, 70, 40]),
                            x_scale=[5, 5, 5, 5, 5], max_nfev=300, diff_step=1e-3)
        c_ = float(np.sum(rr_.fun ** 2))
        if bt is None or c_ < bt[0]: bt = (c_, rr_)
    h = bt[1].x; gg = [h[0], h[1], h[2], fst, fdp, h[3], h[4]]
    w_, m_ = wrms(gg)
    tests[name] = dict(cost=bt[0], strike=fst, dip=fdp, depth_c=float(h[2]), L=float(h[3]), W=float(h[4]),
                       rake=float(math.degrees(math.atan2(m_[1], m_[0]))), slip_mm=float(math.hypot(m_[0], m_[1])))
    print(name, tests[name])
g = r.x.copy(); g[3] = g[3] % 360
m, G, dvec = solve_lin(g)
ss, ds = m[0], m[1]
slip = math.hypot(ss, ds); rake = math.degrees(math.atan2(ds, ss))

J = r.jac[:-1]; dof = J.shape[0] - J.shape[1] - len(m)
s2 = np.sum(r.fun[:-1] ** 2) / max(dof, 1)
try:
    cov = np.linalg.inv(J.T @ J) * s2; sig = np.sqrt(np.diag(cov))
except np.linalg.LinAlgError:
    sig = np.full(7, np.nan)

covm = np.linalg.inv(G.T @ G) * s2; sm = np.sqrt(np.diag(covm))
xc, yc, zc, st, dp, L, W = g
M0 = MU * L * 1e3 * W * 1e3 * slip * 1e-3; Mw = 2 / 3 * (math.log10(M0) - 9.1)
top = zc - math.sin(math.radians(dp)) * W / 2; bot = zc + math.sin(math.radians(dp)) * W / 2
clon, clat = LON0 + xc / KX, LAT0 + yc / KY

out = {}
for k, t in enumerate(T):
    F = T[t]["full"]; x, y = km(F["LO"], F["LA"])
    uE, uN, uZ = okada85(x - xc, y - yc, zc, st, dp, L, W, ss=ss, ds=ds)
    mod = uE * F["e"] + uN * F["n"] + uZ * F["u"]

    me, mn, mu_ = np.nanmean(F["e"]), np.nanmean(F["n"]), np.nanmean(F["u"])
    modf = uE * np.where(np.isfinite(F["e"]), F["e"], me) + uN * np.where(np.isfinite(F["n"]), F["n"], mn) + \
        uZ * np.where(np.isfinite(F["u"]), F["u"], mu_)
    ramp = np.zeros_like(F["los"])
    for kk, (tt, f) in enumerate(frames):
        if tt != t: continue
        msk = F["fr"] == f
        ramp[msk] = m[2 + NR * kk] + m[3 + NR * kk] * x[msk] / 50 + m[4 + NR * kk] * y[msk] / 50 + \
            m[5 + NR * kk] * F["h"][msk]
    obs = F["los"] - ramp; res = obs - mod
    v = np.isfinite(obs)
    np.savetxt(os.path.join(wd, f"{t}_obs.xyz"), np.column_stack([F["LO"][v], F["LA"][v], obs[v]]), fmt="%.4f %.4f %.2f")
    np.savetxt(os.path.join(wd, f"{t}_res.xyz"), np.column_stack([F["LO"][v], F["LA"][v], res[v]]), fmt="%.4f %.4f %.2f")
    np.savetxt(os.path.join(wd, f"{t}_mod.xyz"), np.column_stack([F["LO"].ravel(), F["LA"].ravel(), modf.ravel()]), fmt="%.4f %.4f %.2f")

    Dt = T[t]; bl, _ = design(g)
    pred = bl[k] @ m
    rr = Dt["d"] - pred
    out[t] = dict(n=int(Dt["x"].size), rms_mm=float(np.sqrt(np.mean(rr ** 2))),
                  var_red=float(1 - np.var(rr) / np.var(Dt["d"] - bl[k][:, 2:] @ m[2:])),
                  obs_min=float(np.nanpercentile(obs, 0.5)), obs_max=float(np.nanpercentile(obs, 99.5)))

C, ztop, zbot = corners(xc, yc, zc, st, dp, L, W)
with open(os.path.join(wd, "fault.txt"), "w") as f:
    for cx, cy in list(C) + [C[0]]: f.write(f"{LON0 + cx / KX:.5f} {LAT0 + cy / KY:.5f}\n")
with open(os.path.join(wd, "fault_top.txt"), "w") as f:
    for cx, cy in C[:2]: f.write(f"{LON0 + cx / KX:.5f} {LAT0 + cy / KY:.5f}\n")

ud = np.array([math.cos(math.radians(st)), -math.sin(math.radians(st))])
dsurf = ztop / math.tan(math.radians(dp))
with open(os.path.join(wd, "fault_surf.txt"), "w") as f:
    for cx, cy in C[:2]:
        sx, sy = np.array([cx, cy]) - ud * dsurf; f.write(f"{LON0 + sx / KX:.5f} {LAT0 + sy / KY:.5f}\n")

res = dict(lon=clon, lat=clat, xc=xc, yc=yc, depth_c=zc, top=ztop, bottom=zbot, strike=st, dip=dp, L=L, W=W,
           ss_mm=ss, ds_mm=ds, slip_mm=slip, rake=rake, M0=M0, Mw=Mw,
           sig=dict(zip(["xc", "yc", "depth_c", "strike", "dip", "L", "W"], [float(v) for v in sig])),
           sig_ss=float(sm[0]), sig_ds=float(sm[1]), tracks=out, cost=cost, tests=tests,
           k_elev_mm_per_km={t: float(m[5 + NR * kk]) for kk, (t, f) in enumerate(frames)})
json.dump(res, open(os.path.join(wd, "best.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
