#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}; SRC=$(cd "$(dirname "$0")" && pwd)
mkdir -p work_insar && cd work_insar
cp $SRC/okada85.py $SRC/insar_invert.py .

bash $SRC/insar_get.sh > get.log 2>&1
python3 $SRC/insar_prep.py . '{"desc":["151D_05440","151D_05241"],"asc":["144A_05269","144A_05476"]}' > prep.json
python3 insar_invert.py . > invert.log
python3 - > par.sh <<'PY'
import json; b = json.load(open("best.json")); s = b["sig"]; t = b["tracks"]; T = b["tests"]
print(f'STK={b["strike"]:.0f}; DIP={b["dip"]:.0f}; RAKE={b["rake"]:.0f}; LEN={b["L"]:.0f}; WID={b["W"]:.0f}')
print(f'SLIP={b["slip_mm"]:.0f}; ZC={b["depth_c"]:.1f}; ZT={b["top"]:.1f}; ZB={b["bottom"]:.1f}; MW={b["Mw"]:.2f}')
print(f'RA={t["asc"]["rms_mm"]:.1f}; RD={t["desc"]["rms_mm"]:.1f}; VA={100*t["asc"]["var_red"]:.0f}; VD={100*t["desc"]["var_red"]:.0f}')
print(f'C0={b["cost"]:.2f}; C1={T["NP1"]["cost"]:.2f}; C2={T["NP2"]["cost"]:.2f}')
PY
. ./par.sh; cat best.json

curl -s -m 180 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
QSHP=$(find . -name 'faults.shp' | head -1)
g() { s=""; [ -f "${2%.shp}.prj" ] || s="-s_srs EPSG:4326"; ogr2ogr -f OGR_GMT "$1" "$2" $s -t_srs EPSG:4326 "${@:3}"; }
g qf.gmt "$QSHP"; g amf.gmt "$QSHP" -where "Name LIKE '%Alburz%'"

R=66.85/67.95/36.15/36.95; J=L67.4/36.55/36.2/36.9/5.1c
for t in asc desc; do for k in obs mod res; do gmt xyz2grd ${t}_$k.xyz -R$R -I0.002 -r -G${t}_$k.nc; done; done
gmt grdcut @earth_relief_03s -R66.8/68.0/36.1/37.0 -Grel.nc || { sleep 20; gmt grdcut @earth_relief_03s -R66.8/68.0/36.1/37.0 -Grel.nc; }
gmt grdgradient rel.nc -A315 -Nt0.7 -Gint.nc
gmt makecpt -Cgray -T-3000/5000 > gray.cpt
gmt makecpt -Cturbo -T-25/25 -D > los.cpt
GMTV=$(gmt --version); GSV=$(gs --version); WO='=~1p,white'
read MW_ MH_ <<< "$(gmt mapproject -R$R -J$J -W)"
GX=0.25; TW=$(awk -v w=$MW_ -v g=$GX 'BEGIN{print 3*w+2*g}')

panel() {
  gmt grdimage rel.nc -Iint.nc -Cgray.cpt -R$R -J$J
  gmt grdimage $1 -Clos.cpt -Q -nn
  gmt basemap -Bxa30mf10mg15m -Bya20mf5mg20m -B$3 --MAP_GRID_PEN_PRIMARY=0.25p,gray55
  gmt plot qf.gmt -W0.35p,gray20
  gmt plot amf.gmt -W1.4p,white; gmt plot amf.gmt -W0.7p,220/20/60
  gmt plot fault.txt -W1.6p,white; gmt plot fault.txt -W0.8p,black
  gmt plot fault_top.txt -W3.0p,white; gmt plot fault_top.txt -W1.7p,black
  echo "67.3886 36.5461" | gmt plot -Sa0.36c -Gwhite -W0.8p,black
  echo "$2" | gmt text -F+cTL+f8p,Helvetica,black$WO -D0.12c/-0.12c
  if [ -n "$4" ]; then echo "$4" | gmt text -F+cBR+f7.5p,Helvetica,black$WO -D-0.12c/0.12c; fi
}

gmt begin $OUT/insar_model pdf,png E300
gmt set FONT_ANNOT_PRIMARY 7.5p,Helvetica FONT_LABEL 8p,Helvetica MAP_FRAME_TYPE plain \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel MAP_TICK_LENGTH_PRIMARY 3p PS_PAGE_COLOR white

gmt text -R0/$TW/0/2 -JX${TW}c/2c -M -N -F+f8p,Helvetica,gray20+jTL -Y0.2c <<EOF
> 0 2.95 0.35c ${TW}c j
Depicts: Sentinel-1 co-seismic line-of-sight (LOS) displacement of the 2 November 2025 earthquake from COMET-LiCSAR unwrapped interferograms (ascending track 144, frame 144A_05269, 29 October-10 November 2025; descending track 151, frames 151D_05440 and 151D_05241, 30 October-11 November 2025; coherence @~\263@~0.3, 200 m posting): (a, d) observed after removal of a best-fitting offset, planar ramp and linear phase-elevation term per track (descending frames aligned in their overlap), (b, e) Okada (1985) uniform-slip rectangular dislocation (strike ${STK}@., dip ${DIP}@., rake ${RAKE}@., ${LEN} km @~\264@~ ${WID} km, ${SLIP} mm slip, ${ZT}-${ZB} km depth, M@-w@- ${MW}), (c, f) residual. Black rectangle: surface projection of the modelled fault (thick line: upper edge); star: USGS epicentre; red: Alburz-Marmul Fault; grey: Quaternary faults (USGS OFR 2007-1103); shaded relief: SRTM 3 arc-sec. Lambert conformal conic projection (central meridian 67.4@., standard parallels 36.2@.N and 36.9@.N). Software: GMT $GMTV, GDAL, GhostScript $GSV, Python 3 (NumPy, SciPy, rasterio). Source: authors.
EOF
gmt colorbar -Clos.cpt -R$R -J$J -Dx$(awk -v w=$TW 'BEGIN{print w/2}')c/3.95c+w8c/0.3c+h+jBC+e \
   -Bxa5f1+l"LOS displacement (mm), positive towards the satellite" -Y0c
gmt plot -T -Y5.1c
panel desc_obs.nc "(d) Descending 151D, observed" WSen ""
gmt plot -T -X$(awk -v w=$MW_ -v g=$GX 'BEGIN{print w+g}')c
panel desc_mod.nc "(e) Descending 151D, model" wSen ""
gmt plot -T -X$(awk -v w=$MW_ -v g=$GX 'BEGIN{print w+g}')c
panel desc_res.nc "(f) Descending 151D, residual" wSen "RMS ${RD} mm"
gmt plot -T -X-$(awk -v w=$MW_ -v g=$GX 'BEGIN{print 2*(w+g)}')c -Y$(awk -v h=$MH_ 'BEGIN{print h+0.25}')c
panel asc_obs.nc "(a) Ascending 144A, observed" Wsen ""
gmt plot -T -X$(awk -v w=$MW_ -v g=$GX 'BEGIN{print w+g}')c
panel asc_mod.nc "(b) Ascending 144A, model" wsen ""
gmt plot -T -X$(awk -v w=$MW_ -v g=$GX 'BEGIN{print w+g}')c
panel asc_res.nc "(c) Ascending 144A, residual" wsen "RMS ${RA} mm"
gmt end
echo "done $STK/$DIP/$RAKE L=$LEN W=$WID slip=$SLIP Mw=$MW rms $RA/$RD VR $VA/$VD costs $C0 NP1=$C1 NP2=$C2"
