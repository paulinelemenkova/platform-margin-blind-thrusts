# platform-margin-blind-thrusts

Blind ramps beneath an evaporite-detached foreland margin, read from legacy subsurface isopachs,
drainage morphometry, crustal seismicity and coseismic radar interferometry on a 0.01-degree lattice.

## Scope

- Supervised susceptibility of active blind thrusting above a Jurassic evaporite decollement:
  Random Forest vs histogram gradient boosting, 0.25-degree spatial block cross-validation,
  permutation importance.
- Uniform-slip rectangular dislocation (Okada, 1985) inverted from ascending and descending
  Sentinel-1 line-of-sight fields with per-track ramp and phase-elevation terms.
- Frequency-magnitude statistics (maximum curvature, Aki-Utsu, Shi-Bolt, bootstrap) and a
  true-scale hypocentral section across the range front.
- Hypsometric integral and normalised channel steepness (theta_ref = 0.45) as uplift proxies.
- Candidate ramp segments scaled to moment magnitude (Wells and Coppersmith, 1994, reverse RLD).

## Layout

| Path | Role |
|---|---|
| `scripts/grass_*.sh` | drainage basins, HI, k_sn (GRASS GIS) |
| `scripts/susceptibility_ml.py` | labels, predictors, cross-validated classifiers, class map |
| `scripts/seg_inventory.py` | candidate segments, strike, length, M_w |
| `scripts/insar_get.sh`, `insar_prep.py`, `insar_invert.py`, `okada85.py` | LOS retrieval, decimation, dislocation inversion |
| `scripts/*.sh`, `*.py` (others) | maps, sections, validation and importance plots (GMT, Matplotlib) |
| `out/` | cross-validation metrics, out-of-fold probabilities, permutation importance, best-fit source |

## Order

```
grass_susceptibility_features.sh  ->  susceptibility_map.sh (calls susceptibility_ml.py)
                                  ->  seg_inventory.py  ->  roc_pr.py, predictor_importance.py
insar_get.sh -> insar_prep.py -> insar_invert.py  ->  insar_model.sh  ->  segments_3d.sh
magnitude_frequency.py, depth_section.sh  ->  conceptual_model.sh
```

## Environment

GMT >= 6.5, GRASS GIS 8.4, GDAL >= 3.8, Ghostscript; Python >= 3.10 with numpy, scipy, pandas,
scikit-learn, matplotlib, rasterio.

## Inputs

Public: USGS OFR 2007-1103 Quaternary faults, ANSS ComCat, Global CMT, COMET-LiCSAR
interferograms, SRTM 3" and SRTM15+ (GMT remote grids), all fetched at run time.
Legacy subsurface vectors (isopachs, structure outlines, wells): archive identifier in
`GDRIVE_ID`, available on request.

## License

MIT
