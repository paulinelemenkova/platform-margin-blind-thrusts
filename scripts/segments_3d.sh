#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}; IN=${IN:-$(pwd)}
mkdir -p work_3d && cd work_3d
W=66.6; E=68.35; S=35.95; N=37.1; ZB=-30; ZT=4; DDET=-5
R=$W/$E/$S/$N; R3=$R/$ZB/$ZT; J=M14c; JZ=5.8c; P=-p205/32
EQLON=67.3886; EQLAT=36.5461; EQDEP=28

gmt xyz2grd $IN/prob.xyz -R$R -I0.01 -r -Gprob.nc
gmt grdmath prob.nc 0 MUL $DDET ADD = det.nc
for i in 1 2 3; do gmt grdcut @earth_relief_30s -R$R -Grel.nc && break; sleep 15; done
gmt grdsample rel.nc -I0.025 -Gtopo_c.nc; gmt grdmath topo_c.nc 1000 DIV = topo.nc
gmt grdsample rel.nc -I0.005 -Gtopo_f.nc; gmt grdmath topo_f.nc 1000 DIV = topo_f.nc
gmt makecpt -Cturbo -T0/1/0.2 -D > p0.cpt
awk -v L="Very low,Low,Moderate,High,Very high" 'BEGIN{split(L,n,",")} /^[0-9]/{k++; print $0" ;"n[k]; next} {print}' p0.cpt > prob.cpt

for f in struct_cand struct_other; do gmt select $IN/$f.gmt -R$R -fg > $f.gmt || true
  awk -v z=$DDET '/^>/{print; next} {print $1, $2, z}' $f.gmt > ${f}_z.gmt; done
gmt select $IN/qf.gmt -R$R -fg > qf.gmt || true
gmt sample1d -h0 qf.gmt -fg -I0.5k 2>/dev/null > qf_d.gmt || cp qf.gmt qf_d.gmt
gmt grdtrack qf_d.gmt -Gtopo_f.nc -fg | awk '/^>/{print; next} NF>=3{print $1, $2, $NF+0.05}' > qf_z.gmt
python3 - > insar.txt <<'PY'
import math
lon0, lat0, dc, st, dp, L, W = 67.2624, 36.6894, 7.69, 270.8, 13.24, 51.26, 40.0
kx = 111.32 * math.cos(math.radians(lat0)); ky = 111.195
s = (math.sin(math.radians(st)), math.cos(math.radians(st))); d = (math.sin(math.radians(st + 90)), math.cos(math.radians(st + 90)))
h = 0.5 * W * math.cos(math.radians(dp)); v = 0.5 * W * math.sin(math.radians(dp))
pts = []
for a, b in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
    x = a * 0.5 * L * s[0] + b * h * d[0]; y = a * 0.5 * L * s[1] + b * h * d[1]
    pts.append((lon0 + x / kx, lat0 + y / ky, -(dc + b * v)))
pts.append(pts[0])
for p in pts: print(f"{p[0]:.4f} {p[1]:.4f} {p[2]:.2f}")
PY
ZS=$(echo "$EQLON $EQLAT" | gmt grdtrack -Gtopo_f.nc | awk '{print $3}')
cat > towns.txt <<'EOF'
67.110 36.709 C01
67.686 36.698 C02
66.897 36.758 C03
68.018 36.262 C04
EOF
awk '{print $1, $2}' towns.txt | gmt grdtrack -Gtopo_f.nc > towns_z.txt
paste -d' ' towns_z.txt <(cut -d' ' -f3- towns.txt) > towns_lab.txt
GMTV=$(gmt --version); GSV=$(gs --version); WO='=~1p,white'

gmt begin $OUT/segments_3d pdf,png E300
gmt set FONT_ANNOT_PRIMARY 7.5p,Helvetica FONT_LABEL 8p,Helvetica MAP_FRAME_TYPE plain FORMAT_GEO_MAP ddd:mmF \
        MAP_GRID_PEN_PRIMARY 0.25p,gray60 PS_PAGE_COLOR white MAP_TICK_LENGTH_PRIMARY 3p
gmt text -R0/16/0/2 -JX16c/2c -M -N -F+f8p,Helvetica,gray20+jTL -Y1.25c <<EOF
> 0 2.55 0.35c 16c j
Depicts: perspective view from the SSW (azimuth 205@., elevation 32@.) of the epicentral area of the 2025 M@-w@- 6.2 earthquake: susceptibility of active blind thrusting (Random Forest; classes as in the susceptibility map) draped on a horizontal surface at the approximate depth of the Gaurdak detachment (5 km), with the candidate blind structures (USGS subsurface structures with mean probability @~\263@~0.6) and the other subsurface structures projected onto it; SRTM15+ topography as a semi-transparent mesh with Quaternary faults and towns; the best-fitting uniform-slip InSAR fault plane (3.1-12.3 km, strike 271@., dip 13@. N); the USGS hypocentre (28 km) and the two ComCat aftershocks (fixed depth 10 km). Topography and depth share one vertical scale (about 1.9 times the horizontal). Mercator projection. Data: USGS petroleum GIS; USGS OFR 2007-1103; USGS ComCat; COMET-LiCSAR; SRTM15+. Software: GMT $GMTV, GhostScript $GSV, Python 3 (scikit-learn). Source: authors.
EOF
gmt colorbar -Cprob.cpt -R0/1/0/1 -JX1c -Dx0c/3.35c+w7.2c/0.3c+h+jBL -Li0.06c --FONT_ANNOT_PRIMARY=7p,Helvetica
echo "0 4.05 Susceptibility class at the detachment level (probability bins of 0.2)" | \
  gmt text -R0/16/0/6 -JX16c/6c -F+f7.5p,Helvetica,black+jBL -N
cat > leg.txt <<EOF
N 2
S 0.25c - 0.45c - 1.0p,black 0.6c Candidate structure (p @~\263@~ 0.6)
S 0.25c - 0.45c - 0.4p,gray30 0.6c Other USGS structure
S 0.25c - 0.45c - 0.6p,220/20/60 0.6c Quaternary fault (surface)
S 0.25c r 0.45c/0.25c 220/20/60@30 0.6p,140/15/20 0.6c InSAR fault plane
S 0.25c a 0.32c white 0.8p,black 0.6c 2025 hypocentre (USGS)
S 0.25c c 0.14c white 0.5p,black 0.6c 2025 aftershock (ComCat)
EOF
gmt legend leg.txt -R0/1/0/1 -JX1c -Dx8.0c/3.25c+w8.0c+jBL --FONT_ANNOT_PRIMARY=7p,Helvetica
gmt basemap -R$R3 -J$J -JZ$JZ $P -Bxa0.5f0.25g0.5 -Bya0.5f0.25g0.5 -Bza10f5+l"Elevation / depth (km)" \
    -BWSneZ+b -Y5.6c
gmt grdview det.nc -R$R3 -J$J -JZ$JZ $P -Gprob.nc -Cprob.cpt -Qi250
gmt plot3d struct_other_z.gmt -R$R3 -J$J -JZ$JZ $P -W0.4p,gray30
gmt plot3d struct_cand_z.gmt -R$R3 -J$J -JZ$JZ $P -W2.0p,white
gmt plot3d struct_cand_z.gmt -R$R3 -J$J -JZ$JZ $P -W1.0p,black
gmt grdview topo.nc -R$R3 -J$J -JZ$JZ $P -Qmwhite -Wm0.2p,gray40 -t86
printf "%s %s %s\n%s %s %s\n" $EQLON $EQLAT $ZS $EQLON $EQLAT -$EQDEP | \
  gmt plot3d -R$R3 -J$J -JZ$JZ $P -W0.7p,black,4_2:0p
gmt plot3d insar.txt -R$R3 -J$J -JZ$JZ $P -G220/20/60@30 -W1.0p,140/15/20 -L
echo "$EQLON $EQLAT -$EQDEP" | gmt plot3d -R$R3 -J$J -JZ$JZ $P -Sa0.5c -Gwhite -W0.9p,black
awk '$3<6{print $1, $2, -$4}' $IN/eq_after.txt | gmt plot3d -R$R3 -J$J -JZ$JZ $P -Sc0.16c -Gwhite -W0.5p,black
gmt plot3d qf_z.gmt -R$R3 -J$J -JZ$JZ $P -W0.6p,220/20/60
awk '{print $1, $2, $3+0.05}' towns_z.txt | gmt plot3d -R$R3 -J$J -JZ$JZ $P -Ss0.18c -Gwhite -W0.7p,black
awk '{z=$3+0.6; n=$4; for(i=5;i<=NF;i++) n=n" "$i; print $1, $2, z, n}' towns_lab.txt | \
  gmt text -R$R3 -J$J -JZ$JZ $P -Z -F+f7.5p,Helvetica,black$WO+jCB
printf "%s %s %s InSAR fault plane\n" 67.05 36.80 -9 | \
  gmt text -R$R3 -J$J -JZ$JZ $P -Z -F+f7.5p,Helvetica-Oblique,140/15/20$WO+jCM
printf "%s %s %s USGS hypocentre (28 km)\n" 67.42 $EQLAT -28 | \
  gmt text -R$R3 -J$J -JZ$JZ $P -Z -F+f7.5p,Helvetica,black$WO+jLM
printf "%s %s %s Gaurdak detachment level (5 km)\n" 68.33 35.97 $DDET | \
  gmt text -R$R3 -J$J -JZ$JZ $P -Z -F+f7.5p,Helvetica-Oblique,black$WO+jRB -D0/0.12c
gmt end
echo done
