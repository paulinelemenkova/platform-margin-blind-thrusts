#!/usr/bin/env python3
import sys, os, json
import numpy as np, rasterio
from rasterio.warp import reproject, Resampling

W, S, E, N = 66.85, 36.15, 67.95, 36.95
DX = 0.002
LAM = 0.055465763
CCMIN = 0.30
lon = np.arange(W + DX / 2, E, DX); lat = np.arange(N - DX / 2, S, -DX)
nx, ny = lon.size, lat.size
dst_tr = rasterio.transform.from_origin(W, N, DX, DX)

def warp(path, rs=Resampling.average):
    out = np.full((ny, nx), np.nan, np.float32)
    with rasterio.open(path) as src:
        a = src.read(1).astype(np.float32)
        nd = src.nodata
        if nd is not None: a[a == nd] = np.nan
        a[a == 0] = np.nan
        reproject(a, out, src_transform=src.transform, src_crs=src.crs, dst_transform=dst_tr,
                  dst_crs="EPSG:4326", resampling=rs, src_nodata=np.nan, dst_nodata=np.nan)
    return out

wd = sys.argv[1]; tracks = json.loads(sys.argv[2])
info = {}
for tr, frames in tracks.items():
    los_list, e_l, n_l, u_l, h_l = [], [], [], [], []
    for f in frames:
        p = os.path.join(wd, f)
        unw = warp(p + ".unw.tif"); cc = warp(p + ".cc.tif")
        if np.nanmax(cc) > 1.5: cc = cc / 255.0
        los = -unw * LAM / (4 * np.pi) * 1000.0
        los[(cc < CCMIN) | ~np.isfinite(cc)] = np.nan
        e = warp(p + ".E.tif", Resampling.bilinear); n = warp(p + ".N.tif", Resampling.bilinear)
        u = warp(p + ".U.tif", Resampling.bilinear); h = warp(p + ".hgt.tif", Resampling.bilinear)
        los_list.append(los); e_l.append(e); n_l.append(n); u_l.append(u); h_l.append(h)

    mos = los_list[0].copy(); E_, N_, U_, H_ = e_l[0].copy(), n_l[0].copy(), u_l[0].copy(), h_l[0].copy()
    FR = np.where(np.isfinite(mos), 0, -1).astype(np.int16)
    offs = []
    for k in range(1, len(los_list)):
        b = los_list[k]; ov = np.isfinite(mos) & np.isfinite(b)
        off = float(np.nanmedian(mos[ov] - b[ov])) if ov.sum() > 200 else 0.0
        offs.append((int(ov.sum()), off)); b = b + off
        fill = ~np.isfinite(mos) & np.isfinite(b); mos[fill] = b[fill]; FR[fill] = k
        both = np.isfinite(mos) & np.isfinite(b) & ~fill; mos[both] = 0.5 * (mos[both] + b[both])
        for A, B in ((E_, e_l[k]), (N_, n_l[k]), (U_, u_l[k]), (H_, h_l[k])):
            m = ~np.isfinite(A) & np.isfinite(B); A[m] = B[m]
    np.savez_compressed(os.path.join(wd, f"{tr}.npz"), lon=lon, lat=lat, los=mos, e=E_, n=N_, u=U_, h=H_, frame=FR)
    v = mos[np.isfinite(mos)]
    info[tr] = dict(frames=frames, overlap_offsets=offs, valid=int(v.size), p1=float(np.percentile(v, 1)),
                    p99=float(np.percentile(v, 99)), meanE=float(np.nanmean(E_)), meanN=float(np.nanmean(N_)),
                    meanU=float(np.nanmean(U_)))
print(json.dumps(info, indent=1))
