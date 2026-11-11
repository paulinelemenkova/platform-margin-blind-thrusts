#!/usr/bin/env python3
import sys, os, json, math
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from matplotlib.path import Path
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.inspection import permutation_importance

wd = sys.argv[1]; os.chdir(wd)
W, E, S, N, DX = 65.6, 69.4, 35.85, 37.45, 0.01
lon = np.round(np.arange(W + DX / 2, E, DX), 4); lat = np.round(np.arange(S + DX / 2, N, DX), 4)
LO, LA = np.meshgrid(lon, lat); nyx = LO.shape
LAT0 = 36.65; KX = 111.32 * math.cos(math.radians(LAT0)); KY = 111.195
def xy(lo, la): return np.column_stack([(np.asarray(lo) - 66) * KX, (np.asarray(la) - LAT0) * KY])
XY = xy(LO.ravel(), LA.ravel())

FEAT = [("elev", "Elevation"), ("slope", "Slope"), ("relief", "Local relief (5 km)"), ("ksn", "k_sn"),
        ("hi", "Hypsometric integral"), ("ghory", "Ghory thickness"), ("ghory_grad", "Ghory thickness gradient"),
        ("seis", "Seismicity density")]
df = pd.DataFrame({"lon": LO.ravel(), "lat": LA.ravel()})
for k, _ in FEAT:
    a = np.loadtxt(f"f_{k}.xyz")
    ix = np.round((a[:, 0] - lon[0]) / DX).astype(int); iy = np.round((a[:, 1] - lat[0]) / DX).astype(int)
    g = np.full(nyx, np.nan); ok = (ix >= 0) & (ix < nyx[1]) & (iy >= 0) & (iy < nyx[0])
    g[iy[ok], ix[ok]] = a[ok, 2]; df[k] = g.ravel()

def read_gmt(fn):
    segs, cur = [], []
    for l in open(fn):
        if l.startswith(">"):
            if cur: segs.append(np.array(cur))
            cur = []
        elif l[0] in "#" or not l.strip(): continue
        else:
            t = l.split(); cur.append((float(t[0]), float(t[1])))
    if cur: segs.append(np.array(cur))
    return segs

def densify(segs, step=0.3):
    pts = []
    for s in segs:
        p = xy(s[:, 0], s[:, 1])
        for a, b in zip(p[:-1], p[1:]):
            n = max(1, int(np.hypot(*(b - a)) / step)); t = np.linspace(0, 1, n + 1)[:, None]
            pts.append(a + (b - a) * t)
        if len(p) == 1: pts.append(p)
    return np.vstack(pts)

st_segs = read_gmt("struct.gmt"); qf_segs = read_gmt("qf.gmt")
dS = cKDTree(densify(st_segs)).query(XY)[0]
for s in st_segs:
    if len(s) > 3 and np.allclose(s[0], s[-1]):
        inside = Path(s).contains_points(np.column_stack([LO.ravel(), LA.ravel()])); dS[inside] = 0
dQ = cKDTree(densify(qf_segs)).query(XY)[0]
POS_DS, POS_DQ, NEG_DS, NEG_DQ = 1.0, 5.0, 5.0, 10.0
for a_, b_, c_, d_ in ((2, 10, 10, 10), (1, 5, 5, 10), (1, 3, 5, 10), (0.5, 3, 5, 15)):
    print(f"labels dS<={a_} dQ<={b_} | neg dS>{c_} dQ>{d_}: pos", int(((dS <= a_) & (dQ <= b_)).sum()),
          "neg", int(((dS > c_) & (dQ > d_)).sum()))
from scipy import ndimage
for k in ("ghory", "ghory_grad"):
    g = df[k].values.reshape(nyx).copy(); holes, nh = ndimage.label(~np.isfinite(g))
    sz = ndimage.sum(np.ones(nyx), holes, index=np.arange(1, nh + 1)); small = np.isin(holes, 1 + np.where(sz < 200)[0])
    idx = ndimage.distance_transform_edt(~np.isfinite(g), return_distances=False, return_indices=True)
    g[small] = g[tuple(idx)][small]; df[k] = g.ravel()
    print(k, "filled", int(small.sum()), "cells in", int((sz < 200).sum()), "gaps")
dom = np.isfinite(df["ghory"].values)
y = np.full(len(df), -1); y[(dS <= POS_DS) & (dQ <= POS_DQ)] = 1; y[(dS > NEG_DS) & (dQ > NEG_DQ)] = 0
y[~dom] = -1
df["y"], df["dS"], df["dQ"] = y, dS, dQ
df["block"] = (np.floor((df.lon - W) / 0.25) * 100 + np.floor((df.lat - S) / 0.25)).astype(int)
lab = df[df.y >= 0].reset_index(drop=True)
X, Y, G = lab[[k for k, _ in FEAT]].values, lab.y.values, lab.block.values
print("labelled", len(lab), "positives", int(Y.sum()), "negatives", int((Y == 0).sum()), "blocks", len(np.unique(G)))

models = {
    "HGB": lambda: HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=15,
                                                  min_samples_leaf=40, l2_regularization=1.0,
                                                  class_weight="balanced", random_state=42),
    "RF": lambda: make_pipeline(SimpleImputer(strategy="median"),
                                RandomForestClassifier(n_estimators=200, min_samples_leaf=5, max_features="sqrt",
                                                       class_weight="balanced_subsample", n_jobs=-1, random_state=42)),
}
cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)
oof = pd.DataFrame({"y": Y, "block": G}); res = {}; imps = []
for name, mk in models.items():
    p = np.zeros(len(Y)); fold = np.zeros(len(Y), int); aucs, aps = [], []
    for k, (tr, te) in enumerate(cv.split(X, Y, G)):
        m = mk().fit(X[tr], Y[tr]); p[te] = m.predict_proba(X[te])[:, 1]; fold[te] = k
        aucs.append(roc_auc_score(Y[te], p[te])); aps.append(average_precision_score(Y[te], p[te]))
        if name == "HGB" or name == "RF":
            pi = permutation_importance(m, X[te], Y[te], scoring="roc_auc", n_repeats=3, random_state=k, n_jobs=1)
            imps.append(pd.DataFrame({"model": name, "fold": k, "feature": [f for f, _ in FEAT],
                                      "drop_auc": pi.importances_mean, "sd": pi.importances_std}))
    oof[f"p_{name}"] = p; oof["fold"] = fold
    res[name] = dict(auc_mean=float(np.mean(aucs)), auc_sd=float(np.std(aucs)), ap_mean=float(np.mean(aps)),
                     ap_sd=float(np.std(aps)), auc_pooled=float(roc_auc_score(Y, p)),
                     ap_pooled=float(average_precision_score(Y, p)), auc_folds=aucs, ap_folds=aps)
    print(name, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in res[name].items() if "folds" not in k})
best = max(res, key=lambda k: res[k]["auc_mean"])
final = models[best]().fit(X, Y)
prob = final.predict_proba(df[[k for k, _ in FEAT]].values)[:, 1]
prob[~dom] = np.nan
df["prob"] = prob
np.savetxt("prob.xyz", df.loc[dom, ["lon", "lat", "prob"]].values, fmt="%.4f %.4f %.4f")
oof.to_csv("oof.csv", index=False); pd.concat(imps).to_csv("importance.csv", index=False)

br = [0, 0.2, 0.4, 0.6, 0.8, 1.0001]; names = ["Very low", "Low", "Moderate", "High", "Very high"]
cls = np.where(dom, np.digitize(np.nan_to_num(prob), br[1:-1]), -1); frac = [float(np.mean(cls[dom] == i)) for i in range(5)]

aft = np.loadtxt("eq_after.txt", ndmin=2)
if aft.size: aft = aft[aft[:, 2] < 6.0]
if aft.size:
    ia = np.round((aft[:, 0] - lon[0]) / DX).astype(int); ja = np.round((aft[:, 1] - lat[0]) / DX).astype(int)
    ok = (ia >= 0) & (ia < nyx[1]) & (ja >= 0) & (ja < nyx[0])
    pa = prob.reshape(nyx)[ja[ok], ia[ok]]; pa = pa[np.isfinite(pa)]; ca = np.digitize(pa, br[1:-1])
    aft_frac = [float(np.mean(ca == i)) if ca.size else 0.0 for i in range(5)]; n_aft = int(ca.size)

    d60 = np.hypot(*(XY - xy([67.3886], [36.5461])[0]).T) <= 60
    hv_area60 = float(np.mean(cls[d60 & dom] >= 3))
else:
    aft_frac, n_aft, hv_area60 = [0] * 5, 0, float("nan")

P = prob.reshape(nyx)
with open("struct_cand.gmt", "w") as fc, open("struct_other.gmt", "w") as fo:
    ncand = 0
    for s in st_segs:
        i = np.clip(np.round((s[:, 0] - lon[0]) / DX).astype(int), 0, nyx[1] - 1)
        j = np.clip(np.round((s[:, 1] - lat[0]) / DX).astype(int), 0, nyx[0] - 1)
        inside = (s[:, 0] >= W) & (s[:, 0] <= E) & (s[:, 1] >= S) & (s[:, 1] <= N)
        if not inside.any(): continue
        pv = P[j[inside], i[inside]]; pv = pv[np.isfinite(pv)]
        mp = float(np.mean(pv)) if pv.size else 0.0
        f = fc if mp >= 0.6 else fo; ncand += mp >= 0.6
        f.write(f"> -Z{mp:.3f}\n" + "".join(f"{a:.5f} {b:.5f}\n" for a, b in s))
met = dict(n_labelled=int(len(lab)), n_pos=int(Y.sum()), n_neg=int((Y == 0).sum()), n_blocks=int(len(np.unique(G))),
           models=res, best=best, class_breaks=br[:-1] + [1.0], class_names=names, area_frac=frac,
           n_aftershocks=n_aft, aftershock_frac=aft_frac, hv_area_frac_60km=hv_area60, n_candidates=int(ncand),
           n_structures=len(st_segs), features=[f for f, _ in FEAT],
           domain_frac=float(dom.mean()))
json.dump(met, open("metrics.json", "w"), indent=1)
print(json.dumps({k: v for k, v in met.items() if k != "models"}, indent=1))
