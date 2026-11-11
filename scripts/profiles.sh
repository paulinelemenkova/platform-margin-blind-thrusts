#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
GDRIVE_ID=${GDRIVE_ID:?}
mkdir -p work_pr && cd work_pr

BW=65.6; BE=69.4; BS=35.85; BN=37.45
R=$BW/$BS/$BE/$BN+r; J=L67.5/36.65/36/37.3/15.2c
LATN=37.15; LATS=35.95; HALF=5
PROFS="A:66.70 B:67.39 C:68.10"
EQLON=67.3886; EQLAT=36.5461

curl -s -f -m 600 -L -o src.zip "https://drive.usercontent.google.com/download?id=$GDRIVE_ID&export=download&confirm=t"
unzip -q -o src.zip -x '__MACOSX/*'
cpshp() { P=$(find . -name "$1.shp" | head -1); for x in shp shx dbf prj; do [ -f "${P%.shp}.$x" ] && cp "${P%.shp}.$x" "$2.$x"; done; return 0; }
cpshp ghoryisafg gis; cpshp structafg struct
curl -s -m 180 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
QSHP=$(find . -name 'faults.shp' | head -1)
g() { s=""; [ -f "${2%.shp}.prj" ] || s="-s_srs EPSG:4326"; ogr2ogr -f OGR_GMT "$1" "$2" $s -t_srs EPSG:4326 "${@:3}"; }
g qf.gmt "$QSHP"; g struct.gmt struct.shp
g qR.gmt "$QSHP" -where "Type='R'"; g qSS.gmt "$QSHP" -where "Type IN ('RL','LL','SS')"; g qN.gmt "$QSHP" -where "Type='N'"
s=""; [ -f gis.prj ] || s="-s_srs EPSG:4326"
ogr2ogr -f CSV gis.csv gis.shp $s -t_srs EPSG:4326 -lco GEOMETRY=AS_WKT
python3 - gis.csv > gis.xyz <<'PY'
import csv, re, sys
for r in csv.DictReader(open(sys.argv[1])):
    for x, y in re.findall(r'(-?\d+\.?\d*) (-?\d+\.?\d*)', r['WKT']): print(x, y, r['THICK'])
PY
gmt blockmedian gis.xyz -R62.8/69.8/35.1/37.9 -I0.01 > gis.bm
gmt surface gis.bm -R62.8/69.8/35.1/37.9 -I0.01 -T0.35 -Ll0 -Ggis_s.nc
gmt grdmask gis.xyz -R62.8/69.8/35.1/37.9 -I0.01 -NNaN/1/1 -S25k -Ggis_m.nc
gmt grdmath gis_s.nc gis_m.nc MUL = gis.nc

gmt grdcut @earth_relief_03s -R65.5/69.5/35.75/37.55 -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.6 -Gint.nc
gmt makecpt -Cetopo1 -G0/8000 -T0/4500 > topo.cpt

for p in $PROFS; do id=${p%%:*}; lon=${p#*:}
  python3 - $lon $LATN $LATS $HALF > sw_$id.pts <<'PY'
import sys, math
lon, latn, lats, half = map(float, sys.argv[1:5])
k = 111.32 * math.cos(math.radians((latn + lats) / 2))
lat = latn
while lat >= lats - 1e-9:
    d = (latn - lat) * 111.195
    for o in range(-int(half), int(half) + 1):
        print(f"{lon + o / k:.6f} {lat:.6f} {d:.3f} {o}")
    lat -= 0.5 / 111.195
PY
  gmt grdtrack sw_$id.pts -Grelief.nc > sw_$id.z
  python3 - sw_$id.z > st_$id.txt <<'PY'
import sys, collections
D = collections.defaultdict(list)
for l in open(sys.argv[1]):
    x, y, d, o, z = l.split(); D[float(d)].append(float(z))
for d in sorted(D):
    v = D[d]; print(f"{d:.2f} {min(v):.1f} {sum(v)/len(v):.1f} {max(v):.1f}")
PY
  awk '{print $1, $2}' st_$id.txt > env_lo_$id.txt
  awk '{print $1, $4}' st_$id.txt | tac > env_hi_$id.txt
  cat env_lo_$id.txt env_hi_$id.txt > env_$id.txt
  awk -v lon=$lon 'NR%11==6{print lon, $2, $3}' sw_$id.pts | gmt grdtrack -Ggis.nc -sa | awk '{print $3, $4}' > gh_$id.txt
  python3 - qf.gmt $lon $LATN $LATS $HALF > fx_$id.txt <<'PY'
import sys, math
f = sys.argv[1]; L, latn, lats, half = map(float, sys.argv[2:6])
dl = half / (111.32 * math.cos(math.radians((latn + lats) / 2)))
feats = []; cur = None
for l in open(f):
    if l.startswith('# @D'):
        p = l[4:].strip().split('|'); cur = [p[1].strip('"'), p[9].strip('"') if len(p) > 9 else '', []]; feats.append(cur)
    elif l.startswith('#') or l.startswith('>') or not l.strip(): continue
    elif cur is not None:
        x, y = map(float, l.split()[:2])
        if abs(x - L) <= dl and lats <= y <= latn: cur[2].append(y)
hits = []
for typ, name, ys in feats:
    if ys:
        t = 'R' if typ == 'R' else ('N' if typ == 'N' else 'S')
        hits.append([(latn - sum(ys) / len(ys)) * 111.195, t, 'AMF' if 'Alburz' in name else '-'])
hits.sort(); out = []
for h in hits:
    if out and out[-1][1] == h[1] and h[0] - out[-1][0] < 4.0:
        if h[2] == 'AMF': out[-1][2] = 'AMF'
        continue
    out.append(h)
for d, t, n in out: print(f"{d:.2f} {t} {n}")
PY
  python3 - $lon $LATN $LATS $HALF > box_$id.txt <<'PY'
import sys, math
lon, latn, lats, half = map(float, sys.argv[1:5]); dl = half / (111.32 * math.cos(math.radians((latn + lats) / 2)))
for x, y in ((lon-dl, latn), (lon+dl, latn), (lon+dl, lats), (lon-dl, lats), (lon-dl, latn)): print(x, y)
PY
done
LEN=$(awk 'END{print $1}' st_A.txt)
ZMAX=$(cat st_*.txt | awk 'BEGIN{m=0}{if($4>m)m=$4}END{print int(m/500+1)*500}')
EQD=$(awk -v a=$LATN -v b=$EQLAT 'BEGIN{print (a-b)*111.195}')
cat fx_*.txt | sort | uniq -c | head -50

cat > cities.txt <<'EOF'
67.110 36.709 CT C01
67.686 36.698 LB C02
68.015 36.265 LM C04
68.860 36.728 LM C06
65.752 36.667 LM C05
EOF
read MW MH <<< "$(gmt mapproject -R$R -J$J -W)"
GMTV=$(gmt --version); GSV=$(gs --version)
WO='=~1p,white'
PH=3.0; PG=1.45

profile() {
  id=$1
  gmt basemap -R0/$LEN/0/$ZMAX -JX${MW}c/${PH}c -Bxa20f5+l"Distance from ${id} (km)" -Bya1000f250+l"Elevation (m)" -BWSn
  gmt plot env_$id.txt -G200/200/200 -W0.25p,gray50 -L
  awk '{print $1,$3}' st_$id.txt | gmt plot -W1.1p,black
  while read d t n; do
    case $t in R) c=220/20/60;; N) c=0/70/255;; *) c=black;; esac
    printf "$d $ZMAX\n" | gmt plot -Si0.25c -G$c -W0.4p,white -N -D0/0.13c
  done < fx_$id.txt
  ad=$(awk '$3=="AMF"{print $1; exit}' fx_$id.txt)
  if [ -n "$ad" ]; then
    printf "$ad 0\n$ad $ZMAX\n" | gmt plot -W0.7p,220/20/60,4_2:0p
    echo "$ad $ZMAX Alburz-Marmul Fault" | gmt text -F+f8p,Helvetica-Oblique,220/20/60$WO+jLB -D0.08c/-0.42c
  fi
  if [ "$id" = "B" ]; then echo "$EQD $ZMAX" | gmt plot -Sa0.42c -Gwhite -W0.9p,black -N -D0/0.13c; fi
  gmt basemap -R0/$LEN/0/1000 -JX${MW}c/${PH}c -Bya200f50+l"Ghory thickness (m)" -BE --MAP_FRAME_PEN=0.6p,200/100/0 \
      --FONT_LABEL=9p,Helvetica,200/100/0 --FONT_ANNOT_PRIMARY=8p,Helvetica,200/100/0 --MAP_TICK_PEN_PRIMARY=0.5p,200/100/0
  gmt plot gh_$id.txt -W1.4p,200/100/0
  echo "$2 ${id}-${id}@~\242@~" | gmt text -F+cTL+f9p,Helvetica,black$WO -D0.2c/-0.15c
  echo "${id}" | gmt text -F+cTL+f9p,Helvetica-Bold,black -D-0.05c/0.42c -N
  echo "${id}@~\242@~" | gmt text -F+cTR+f9p,Helvetica-Bold,black -D0.05c/0.42c -N
}

gmt begin $OUT/profiles pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 8p,Helvetica \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray50 \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel MAP_TICK_LENGTH_PRIMARY 4p PS_PAGE_COLOR white
gmt set MAP_FRAME_TYPE plain
gmt basemap -R0/1/0/1 -JX1c -B+n -Y5.45c
profile C "(d)"
cat > leg.txt <<EOF
N 3
S 0.3c r 0.5c/0.3c 200/200/200 0.25p,gray50 0.75c Swath min-max (10 km)
S 0.3c - 0.55c - 1.1p,black 0.75c Swath mean elevation
S 0.3c - 0.55c - 1.4p,200/100/0 0.75c Ghory Formation thickness
S 0.3c i 0.25c 220/20/60 0.4p,white 0.75c Quaternary reverse fault
S 0.3c i 0.25c black 0.4p,white 0.75c Quaternary strike-slip fault
S 0.3c i 0.25c 0/70/255 0.4p,white 0.75c Quaternary normal fault
S 0.3c a 0.35c white 0.9p,black 0.75c 2025 mainshock
EOF
gmt legend leg.txt -DJBC+w${MW}c+o0/1.1c -F+gwhite+p0.4p,gray40 --FONT_ANNOT_PRIMARY=8.5p,Helvetica
gmt text -R0/$MW/0/2 -JX${MW}c/2c -M -N -F+f8p,Helvetica,gray20+jTL -Y-4.45c <<EOF
> 0 1.75 0.35c ${MW}c j
Depicts: (a) SRTM 3-arc-second shaded relief (ETOPO1 colour ramp) with petroleum-exploration structures, Quaternary faults and three 10-km-wide north-south swaths across the Alburz-Marmul front; (b-d) swath elevation (minimum-maximum envelope and mean, sampled every 0.5 km on 11 parallel lines 1 km apart), Ghory Formation thickness along the swath axis (gridded USGS isopachs) and Quaternary faults within each swath; the 2025 mainshock is projected on B-B@~\242@~; vertical exaggeration about 9. Lambert conformal conic projection (central meridian 67.5@.E, standard parallels 36@.N and 37.3@.N). Data: SRTM 3 arc-sec; USGS petroleum GIS (2006); USGS OFR 2007-1103; GSHHG. Software: GMT $GMTV, GDAL, GhostScript $GSV, Python 3. Source: authors.
EOF
gmt basemap -R0/1/0/1 -JX1c -B+n -Y$(awk -v h=$PH -v g=$PG 'BEGIN{print 4.45+h+g}')c
profile B "(c)"
gmt basemap -R0/1/0/1 -JX1c -B+n -Y$(awk -v h=$PH -v g=$PG 'BEGIN{print h+g}')c
profile A "(b)"
gmt set MAP_FRAME_TYPE fancy
gmt grdimage relief.nc -Iint.nc -Ctopo.cpt -R$R -J$J -Y$(awk -v h=$PH 'BEGIN{print h+1.45}')c
gmt coast -Di -W0.3p,gray25 -N1/0.9p,gray10,- -I1/0.6p,30/110/200 -A20 -Bxa1f10mg30m -Bya30mf10mg30m -BWseN
gmt plot struct.gmt -W0.35p,gray30
for k in "qN:0/70/255" "qSS:black" "qR:220/20/60"; do f=${k%%:*}; c=${k#*:}
  gmt plot $f.gmt -W1.6p,white; gmt plot $f.gmt -W0.8p,$c; done
for p in $PROFS; do id=${p%%:*}; lon=${p#*:}
  gmt plot box_$id.txt -W1.6p,white -Ap; gmt plot box_$id.txt -W0.8p,black -Ap
  printf "$lon $LATN\n$lon $LATS\n" | gmt plot -W1.0p,black,4_2:0p
  echo "$lon $LATN $id" | gmt text -F+f9p,Helvetica-Bold,black$WO+jCB -D0/0.12c
  echo "$lon $LATS ${id}@~\242@~" | gmt text -F+f9p,Helvetica-Bold,black$WO+jCT -D0/-0.12c
done
echo "$EQLON $EQLAT" | gmt plot -Sa0.45c -Gwhite -W0.9p,black
awk '{print $1,$2}' cities.txt | gmt plot -Ss0.2c -Gwhite -W0.8p,black
while read lo la j nm; do
  case $j in CT) dx=0; dy=-0.16c;; LB) dx=0.15c; dy=0.1c;; *) dx=0.15c; dy=0;; esac
  echo "$lo $la $nm" | gmt text -F+f8p,Helvetica,black$WO+j$j -D$dx/$dy
done < cities.txt
echo "(a) Swath profile locations" | gmt text -F+cTL+f9p,Helvetica,black=~1.2p,white -D0.25c/-0.25c
gmt colorbar -Ctopo.cpt -DJRM+w${MH}c/0.32c+o0.3c/0+ef -Bxa1000f100+l"Elevation (m)"
gmt basemap -LjBR+c36.6+w50k+f+l"km"+o0.4c/0.4c -F+gwhite@15+p0.4p,gray40+r0.1c
gmt end
echo "done LEN=$LEN ZMAX=$ZMAX EQD=$EQD"
