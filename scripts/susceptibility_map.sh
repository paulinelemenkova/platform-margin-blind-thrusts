#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}; SRC=$(cd "$(dirname "$0")" && pwd)
HI=${HI:-$SRC/hi01.tif}; KSN=${KSN:-$SRC/ksn01.tif}
GDRIVE_ID=${GDRIVE_ID:?}
mkdir -p work_su && cd work_su

W=65.6; E=69.4; S=35.85; N=37.45; RA=$W/$E/$S/$N
R=$W/$S/$E/$N+r; J=L67.5/36.65/36/37.3/16c
EQLON=67.3886; EQLAT=36.5461

curl -s -f -m 600 -L -o src.zip "https://drive.usercontent.google.com/download?id=$GDRIVE_ID&export=download&confirm=t"
unzip -q -o src.zip -x '__MACOSX/*'
cpshp() { P=$(find . -name "$1.shp" | head -1); for x in shp shx dbf prj; do [ -f "${P%.shp}.$x" ] && cp "${P%.shp}.$x" "$2.$x"; done; return 0; }
cpshp ghoryisafg gis; cpshp qezelisafg qis; cpshp structafg struct
curl -s -m 180 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
QSHP=$(find . -name 'faults.shp' | head -1)
g() { s=""; [ -f "${2%.shp}.prj" ] || s="-s_srs EPSG:4326"; ogr2ogr -f OGR_GMT "$1" "$2" $s -t_srs EPSG:4326 "${@:3}"; }
g struct.gmt struct.shp; g qf.gmt "$QSHP" -spat 64.6 35.1 70.4 38.2 -spat_srs EPSG:4326
ogrinfo -so -al struct.shp | grep -E "Geometry|Feature Count"
for k in gis qis; do
  s=""; [ -f $k.prj ] || s="-s_srs EPSG:4326"
  ogr2ogr -f CSV $k.csv $k.shp $s -t_srs EPSG:4326 -lco GEOMETRY=AS_WKT
  python3 - $k.csv > $k.xyz <<'PY'
import csv, re, sys
for r in csv.DictReader(open(sys.argv[1])):
    for x, y in re.findall(r'(-?\d+\.?\d*) (-?\d+\.?\d*)', r['WKT']):
        print(x, y, r['THICK'])
PY
done
U=https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv
curl -s -m 180 "$U&starttime=1973-01-01&endtime=2025-11-02T20:28:00&minlatitude=35.4&maxlatitude=37.9&minlongitude=65.1&maxlongitude=69.9&minmagnitude=2.5&maxdepth=70" \
  | awk -F, 'NR>1{print $3,$2,$5}' > eq_pre.txt
curl -s -m 180 "$U&starttime=2025-11-02T20:28:00&endtime=2026-11-10&latitude=$EQLAT&longitude=$EQLON&maxradiuskm=60" \
  | awk -F, 'NR>1 && $5<6.0{print $3,$2,$5,$4}' > eq_after.txt
NPRE=$(wc -l < eq_pre.txt); NAFT=$(wc -l < eq_after.txt); echo "pre=$NPRE after=$NAFT"

for i in 1 2 3; do gmt grdcut @earth_relief_15s -R65.2/69.8/35.5/37.8 -Grel15.nc && break; sleep 15; done
gmt grdgradient rel15.nc -A315 -Nt0.65 -Gint15.nc
gmt grdsample rel15.nc -R$RA -I0.01 -r -Gelev.nc
gmt grdgradient rel15.nc -fg -D -Sgrad15.nc -Gjunk.nc
gmt grdmath grad15.nc ATAND = slope15.nc; gmt grdsample slope15.nc -R$RA -I0.01 -r -Gslope.nc
gmt grdfilter rel15.nc -Fu5 -D2 -Gup.nc; gmt grdfilter rel15.nc -Fl5 -D2 -Glo.nc
gmt grdmath up.nc lo.nc SUB = rlf15.nc; gmt grdsample rlf15.nc -R$RA -I0.01 -r -Grelief.nc
gmt grdconvert $HI -Ghi.nc; gmt grdconvert $KSN -Gksn.nc
for k in gis qis; do
  gmt blockmedian $k.xyz -R62.8/69.8/35.1/37.9 -I0.01 > $k.bm
  gmt surface $k.bm -R62.8/69.8/35.1/37.9 -I0.01 -T0.35 -Ll0 -G${k}_s.nc
  gmt grdmask $k.xyz -R62.8/69.8/35.1/37.9 -I0.01 -NNaN/1/1 -S25k -G${k}_m.nc
  gmt grdmath ${k}_s.nc ${k}_m.nc MUL = ${k}_f.nc
done
gmt grdgradient gis_s.nc -fg -D -Sgg.nc -Gjunk2.nc
gmt grdmath gg.nc 1000 MUL gis_m.nc MUL = ggrad_f.nc
gmt grdsample gis_f.nc -R$RA -I0.01 -r -Gghory.nc
gmt grdsample qis_f.nc -R$RA -I0.01 -r -Gqezel.nc
gmt grdsample ggrad_f.nc -R$RA -I0.01 -r -Gghory_grad.nc
for k in elev slope relief hi ksn ghory ghory_grad; do gmt grd2xyz $k.nc -s > f_$k.xyz; done
python3 - <<'PY'
import numpy as np, math
W, E, S, N, DX = 65.6, 69.4, 35.85, 37.45, 0.01
lon = np.arange(W + DX / 2, E, DX); lat = np.arange(S + DX / 2, N, DX); LO, LA = np.meshgrid(lon, lat)
kx = 111.32 * math.cos(math.radians(36.65)); a = np.loadtxt("eq_pre.txt", ndmin=2)
d = np.zeros(LO.shape)
for x, y, m in a:
    d += np.exp(-(((LO - x) * kx) ** 2 + ((LA - y) * 111.195) ** 2) / (2 * 10.0 ** 2))
d = np.log10(1 + 100 * d / (2 * math.pi * 10.0 ** 2))
np.savetxt("f_seis.xyz", np.column_stack([LO.ravel(), LA.ravel(), d.ravel()]), fmt="%.4f %.4f %.5f")
PY

python3 $SRC/susceptibility_ml.py . | tee ml.log
python3 - > par.sh <<'PY'
import json; m = json.load(open("metrics.json")); b = m["best"]; r = m["models"][b]; o = [k for k in m["models"] if k != b][0]
f = m["area_frac"]; a = m["aftershock_frac"]
print(f'BEST={b}; OTHER={o}; AUC={r["auc_mean"]:.2f}; AUCSD={r["auc_sd"]:.2f}; AP={r["ap_mean"]:.2f}')
print(f'AUCO={m["models"][o]["auc_mean"]:.2f}; NPOS={m["n_pos"]}; NNEG={m["n_neg"]}; NBLK={m["n_blocks"]}; NC={m["n_candidates"]}')
print(f'DOM={100*m["domain_frac"]:.0f}; FHV={100*(f[3]+f[4]):.0f}; NA={m["n_aftershocks"]}; AHV={100*(a[3]+a[4]):.0f}; A60={100*m["hv_area_frac_60km"]:.0f}')
print("F=(" + " ".join(f"{100*v:.0f}" for v in f) + ")")
PY
. ./par.sh; cat par.sh

gmt xyz2grd prob.xyz -R$RA -I0.01 -r -Gprob.nc
gmt grdclip prob.nc -Sb0.0005/0.0005 -Sa0.9995/0.9995 -Gprob.nc
echo "p==0 cells: $(awk '$3<=0' prob.xyz | wc -l)"
gmt grdsample int15.nc -Rprob.nc -Gintp.nc; gmt grdmath intp.nc 0.55 MUL = intp.nc
gmt makecpt -Cgray -T-4000/4000 > gray.cpt
gmt makecpt -Cturbo -T0/1/0.2 -D > p0.cpt
awk -v L="Very low,Low,Moderate,High,Very high" 'BEGIN{split(L,n,",")} /^[0-9]/{k++; print $0" ;"n[k]; next} {print}' p0.cpt > prob.cpt
GMTV=$(gmt --version); GSV=$(gs --version); WO='=~1p,white'
cat > towns.txt <<'EOF'
67.110 36.709 LT C01
68.860 36.728 LM C06
67.686 36.698 LB C02
66.897 36.758 RM C03
65.752 36.667 LM C05
68.018 36.262 LM C04
68.713 35.944 RM C08
67.278 37.224 LM C07
EOF

gmt begin $OUT/susceptibility_map pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 8p,Helvetica MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel MAP_GRID_PEN_PRIMARY 0.25p,gray55 PS_PAGE_COLOR white
gmt text -R0/16/0/2 -JX16c/2c -M -N -F+f8p,Helvetica,gray20+jTL -Y1.45c <<EOF
> 0 2.6 0.35c 16c j
Depicts: probability of active blind thrusting predicted by ${BEST} (spatial block cross-validation, 0.25@. blocks, 5 folds: ROC AUC ${AUC} @~\261@~ ${AUCSD}, average precision ${AP}; alternative model ${OTHER}: AUC ${AUCO}) on a 0.01@. grid from eight predictors (elevation, slope, 5-km relief, k@-sn@-, hypsometric integral, Ghory isopach thickness and its gradient, pre-2025 crustal seismicity density), restricted to the ${DOM}% of the study area covered by Ghory isopachs (grey: no subsurface data). Training labels (class-balanced): ${NPOS} cells within 1 km of USGS subsurface structures and 5 km of Quaternary faults versus ${NNEG} cells over 5 km from structures and 10 km from Quaternary faults. Classes cover ${F[0]}, ${F[1]}, ${F[2]}, ${F[3]} and ${F[4]}% of the modelled area; ${NC} subsurface structures have mean probability @~\263@~0.6; ${NA} of the few ComCat 2025 aftershocks lie in the modelled area (${AHV}% in the high and very high classes). Lambert conformal conic projection (central meridian 67.5@., standard parallels 36@.N and 37.3@.N). Data: USGS petroleum GIS (2006); USGS OFR 2007-1103; USGS ComCat; SRTM. Software: GMT $GMTV, GRASS GIS 8.4, GDAL, GhostScript $GSV, Python 3 (scikit-learn). Source: authors.
EOF
gmt colorbar -Cprob.cpt -R0/1/0/1 -JX1c -Dx0c/3.45c+w7.4c/0.32c+h+jBL -Li0.06c -Y0c \
   --FONT_ANNOT_PRIMARY=7.5p,Helvetica
echo "0 4.25 Susceptibility class (probability bins of 0.2: 0-0.2, ..., 0.8-1)" | gmt text -R0/16/0/6 -JX16c/6c -F+f8p,Helvetica,black+jBL -N
cat > leg.txt <<EOF
N 2
S 0.25c - 0.45c - 1.0p,black 0.6c Candidate structure (p @~\263@~ 0.6)
S 0.25c - 0.45c - 0.5p,gray25 0.6c Other USGS structure
S 0.25c - 0.45c - 0.6p,220/20/60 0.6c Quaternary fault
S 0.25c c 0.14c white 0.4p,black 0.6c 2025 aftershock
S 0.25c a 0.32c white 0.8p,black 0.6c 2025 mainshock
S 0.25c s 0.17c white 0.7p,black 0.6c Town
EOF
gmt legend leg.txt -R0/1/0/1 -JX1c -Dx8.2c/3.2c+w7.8c+jBL --FONT_ANNOT_PRIMARY=7p,Helvetica
gmt basemap -R$R -J$J -B+n -Y5.35c
gmt grdimage rel15.nc -Iint15.nc -Cgray.cpt
gmt grdimage prob.nc -Iintp.nc -Cprob.cpt -Q
gmt coast -Di -N1/0.9p,gray10,- -I1/0.7p,30/110/200 -A20 -Bxa1f0.25g0.5 -Bya0.5f0.25g0.5 -BWSnE
gmt plot struct_other.gmt -W0.5p,gray25
gmt plot qf.gmt -W0.6p,220/20/60
gmt plot struct_cand.gmt -W2.0p,white; gmt plot struct_cand.gmt -W1.0p,black
awk '{print $1,$2}' eq_after.txt | gmt plot -Sc0.14c -Gwhite -W0.4p,black
echo "$EQLON $EQLAT" | gmt plot -Sa0.42c -Gwhite -W0.9p,black
awk '{print $1,$2}' towns.txt | gmt plot -Ss0.18c -Gwhite -W0.7p,black
while read lo la j nm; do case $j in RT|RB|RM) dx=-0.14c;; *) dx=0.14c;; esac
  case $j in RT|LT) dy=0.1c;; RB|LB) dy=-0.1c;; *) dy=0;; esac
  echo "$lo $la $nm" | gmt text -F+f8p,Helvetica,black$WO+j$j -D$dx/$dy; done < towns.txt
echo "68.9 37.03 Amu Darya" | gmt text -F+f7.5p,Helvetica-Oblique,30/90/170$WO+a-8+jCM
gmt basemap -LjBR+c36.6+w25k+f+l"km"+o0.3c/0.35c -F+gwhite@15+p0.4p,gray40+r0.1c
gmt end
echo "done"
