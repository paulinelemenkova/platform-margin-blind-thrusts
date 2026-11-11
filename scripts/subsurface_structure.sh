#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
GDRIVE_ID=${GDRIVE_ID:?}
mkdir -p work_ss && cd work_ss

W=63.0; E=69.6; S=35.3; N=37.7
R=$W/$S/$E/$N+r; J=L66.3/36.5/35.8/37.2/15.2c
GR=62.8/69.8/35.1/37.9
BW=65.6; BE=69.4; BS=35.85; BN=37.45
EQLON=67.3886; EQLAT=36.5461

curl -s -f -m 600 -L -o src.zip "https://drive.usercontent.google.com/download?id=$GDRIVE_ID&export=download&confirm=t"
unzip -q -o src.zip -x '__MACOSX/*'
cpshp() { P=$(find . -name "$1.shp" | head -1); for x in shp shx dbf prj; do [ -f "${P%.shp}.$x" ] && cp "${P%.shp}.$x" "$2.$x"; done; return 0; }
cpshp ghoryisafg gis; cpshp qezelisafg qis; cpshp structafg struct; cpshp fldafg fld; cpshp welllogafg well
curl -s -m 180 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
QSHP=$(find . -name 'faults.shp' | head -1)
g() { s=""; [ -f "${2%.shp}.prj" ] || s="-s_srs EPSG:4326"; ogr2ogr -f OGR_GMT "$1" "$2" $s -t_srs EPSG:4326 "${@:3}"; }
g struct.gmt struct.shp; g fld.gmt fld.shp; g well.gmt well.shp
g qR.gmt "$QSHP" -where "Type='R'"; g qSS.gmt "$QSHP" -where "Type IN ('RL','LL','SS')"; g qN.gmt "$QSHP" -where "Type='N'"
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
fc() { ogrinfo -so -al "$@" | awk -F': ' '/Feature Count/{print $2}'; }
NG=$(fc gis.shp); NQ=$(fc qis.shp); NSTR=$(fc struct.shp); NFLD=$(fc fld.shp); NWELL=$(fc well.shp)

gmt grdcut @earth_relief_15s -R$GR -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.6 -Gint15.nc
for k in gis qis; do
  gmt blockmedian $k.xyz -R$GR -I0.01 > $k.bm
  gmt surface $k.bm -R$GR -I0.01 -T0.35 -Ll0 -G${k}_s.nc
  gmt grdmask $k.xyz -R$GR -I0.01 -NNaN/1/1 -S25k -G${k}_m.nc
  gmt grdmath ${k}_s.nc ${k}_m.nc MUL = ${k}.nc
done
gmt grdsample int15.nc -Rgis.nc -Gint.nc
gmt grdmath int.nc 0.75 MUL = int.nc
gmt makecpt -Cgray -T-5000/5000 > gray.cpt
gmt makecpt -CSCM/batlow -T0/900/25 -Z > g.cpt
gmt makecpt -CSCM/batlow -T0/200/5 -Z > q.cpt

cat > cities.txt <<'EOF'
67.110 36.709 RM C01
68.860 36.728 LM C06
65.752 36.667 LM C05
65.122 36.946 LM Andkhoy
65.938 36.216 LM Sar-e Pol
64.746 35.930 LM Maymana
68.708 35.944 LM C08
67.278 37.224 LM C07
EOF
box() { python3 -c "
W,E,S,N=$BW,$BE,$BS,$BN; n=60
for i in range(n+1): print(W+(E-W)*i/n,S)
for i in range(n+1): print(E,S+(N-S)*i/n)
for i in range(n+1): print(E-(E-W)*i/n,N)
for i in range(n+1): print(W,N-(N-S)*i/n)"; }
box > box.txt
read MW MH <<< "$(gmt mapproject -R$R -J$J -W)"
GMTV=$(gmt --version); GDALV=$(gdal-config --version 2>/dev/null || ogr2ogr --version | awk '{print $2}' | tr -d ,); GSV=$(gs --version)
WO='=~0.9p,white'

panel() {
  gmt grdimage relief.nc -Iint15.nc -Cgray.cpt -t35
  gmt grdimage $1.nc -Iint.nc -C$2 -Q
  gmt grdcontour $1.nc -C$3 -A$4+f8p,Helvetica,black$WO -W0.35p,gray20 -Wa0.6p,black -S4 -Gd6c
  gmt coast -Di -W0.3p,gray25 -N1/0.9p,gray10,- -I1/0.6p,30/110/200 -I2/0.4p,30/110/200 -A20 \
            -Bxa1f15mg1 -Bya1f15mg1 -BWSen
  gmt plot struct.gmt -W0.35p,60/60/60
  gmt plot fld.gmt -Gwhite@20 -W0.7p,black
  for k in "qN:0/70/255" "qSS:black" "qR:220/20/60"; do f=${k%%:*}; c=${k#*:}
    gmt plot $f.gmt -W1.6p,white; gmt plot $f.gmt -W0.8p,$c; done
  gmt plot box.txt -W2p,white; gmt plot box.txt -W1p,200/0/0,6_3:0p
  grep -v '^#' well.gmt | grep -v '^>' | gmt plot -Sc0.18c -Gblack -W0.7p,white
  echo "$EQLON $EQLAT" | gmt plot -Sa0.45c -Gred -W0.6p,black
  awk '{print $1,$2}' cities.txt | gmt plot -Ss0.2c -Gwhite -W0.8p,black
  while read lo la j nm; do
    case $j in RM) dx=-0.15c;; *) dx=0.15c;; esac
    echo "$lo $la $nm" | gmt text -F+f8p,Helvetica,black$WO+j$j -D$dx/0
  done < cities.txt
  echo "$5 $6" | gmt text -F+cTL+f9p,Helvetica,black=~1.2p,white -D0.25c/-0.25c
  gmt colorbar -C$2 -DJRM+w${MH}c/0.32c+o0.3c/0+ef -Bxa$7f$(( $7 / 5 ))+l"Thickness (m)" \
            --FONT_ANNOT_PRIMARY=8p,Helvetica --FONT_LABEL=9p,Helvetica
}

gmt begin $OUT/subsurface_structure pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 8p,Helvetica \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray50 \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel MAP_TICK_LENGTH_PRIMARY 4p PS_PAGE_COLOR white
gmt basemap -R$R -J$J -B+n -Y5.6c
panel qis q.cpt 20 40 "(b)" "Hauterivian Qezeltash Formation, thickness" 40
gmt basemap -LjBR+c36.5+w100k+f+l"km"+o0.4c/0.4c -F+gwhite@15+p0.4p,gray40+r0.1c
cat > leg.txt <<EOF
N 3
S 0.3c - 0.55c - 0.35p,60/60/60 0.75c Subsurface structure ($NSTR)
S 0.3c r 0.5c/0.32c white 0.7p,black 0.75c Oil and gas field ($NFLD)
S 0.3c c 0.18c black 0.7p,gray40 0.75c Well with logs ($NWELL)
S 0.3c - 0.55c - 1.2p,220/20/60 0.75c Quaternary reverse fault
S 0.3c - 0.55c - 1.2p,black 0.75c Quaternary strike-slip fault
S 0.3c - 0.55c - 1.2p,0/70/255 0.75c Quaternary normal fault
S 0.3c - 0.55c - 0.6p,black 0.75c Isopach (annotated)
S 0.3c a 0.35c red 0.6p,black 0.75c 2025 mainshock
S 0.3c - 0.55c - 1.0p,200/0/0,6_3:0p 0.75c Study area
EOF
gmt legend leg.txt -DJBC+w${MW}c+o0/1.0c -F+gwhite+p0.4p,gray40 --FONT_ANNOT_PRIMARY=8.5p,Helvetica
gmt text -R0/$MW/0/2 -JX${MW}c/2c -M -N -F+f8p,Helvetica,gray20+jTL -Y-4.4c <<EOF
> 0 1.75 0.35c ${MW}c j
Depicts: thickness of (a) the Upper Campanian-Paleocene Ghory Formation and (b) the Hauterivian Qezeltash Formation, gridded from $NG and $NQ USGS isopachs (GMT surface, tension 0.35, 0.01@. grid, masked beyond 25 km from data) over SRTM15+ shaded relief, with $NSTR petroleum-exploration structures, $NFLD oil and gas fields, $NWELL wells with logs and Quaternary faults; Lambert conformal conic projection (central meridian 66.3@.E, standard parallels 35.8@.N and 37.2@.N). Data: USGS petroleum GIS (2006); USGS OFR 2007-1103; SRTM15+; GSHHG; USGS ComCat. Software: GMT $GMTV, GDAL $GDALV, GhostScript $GSV, Python 3. Source: authors.
EOF
SH=$(awk -v h=$MH 'BEGIN{print 4.4 + h + 1.0}')
gmt basemap -R$R -J$J -B+n -Y${SH}c
panel gis g.cpt 50 100 "(a)" "Ghory Formation (Upper Campanian-Paleocene), thickness" 100
gmt basemap -TdjBR+w1.0c+l,,,N+o0.4c/0.4c
gmt end
echo "done NG=$NG NQ=$NQ"; gmt grdinfo -C gis.nc qis.nc
