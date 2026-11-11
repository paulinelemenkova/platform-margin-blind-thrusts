#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}; HI=${HI:-$(pwd)/hi.tif}; KSN=${KSN:-$(pwd)/ksn.xyz.gz}
GDRIVE_ID=${GDRIVE_ID:?}
mkdir -p work_gi && cd work_gi

BW=65.6; BE=69.4; BS=35.85; BN=37.45
R=$BW/$BS/$BE/$BN+r; J=L67.5/36.65/36/37.3/15.2c
EQLON=67.3886; EQLAT=36.5461

curl -s -f -m 600 -L -o src.zip "https://drive.usercontent.google.com/download?id=$GDRIVE_ID&export=download&confirm=t"
unzip -q -o src.zip -x '__MACOSX/*'
cpshp() { P=$(find . -name "$1.shp" | head -1); for x in shp shx dbf prj; do [ -f "${P%.shp}.$x" ] && cp "${P%.shp}.$x" "$2.$x"; done; return 0; }
cpshp structafg struct; cpshp welllogafg well
curl -s -m 180 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
QSHP=$(find . -name 'faults.shp' | head -1)
g() { s=""; [ -f "${2%.shp}.prj" ] || s="-s_srs EPSG:4326"; ogr2ogr -f OGR_GMT "$1" "$2" $s -t_srs EPSG:4326 "${@:3}"; }
g struct.gmt struct.shp; g well.gmt well.shp; g qf.gmt "$QSHP"

gmt grdcut @earth_relief_03s -R65.5/69.5/35.75/37.55 -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.7 -Gint.nc
gmt makecpt -Cgray -T-7000/4000 > gray.cpt
gmt grdconvert "$HI" -Ghi.nc
gmt grdsample int.nc -Rhi.nc -Ginth.nc
gmt grdmath inth.nc 0.6 MUL = inth.nc
gmt makecpt -Cturbo -T0.1/0.65/0.01 -D > hi.cpt
gzip -dc "$KSN" | gmt blockmean -R$BW/$BE/$BS/$BN -I0.0025 > ksn_bm.txt
sort -k3,3g ksn_bm.txt | awk '$3>0' > ksn.txt
gmt makecpt -Cturbo -T1/3/0.02 -Q -D > ksn.cpt
NKS=$(gzip -dc "$KSN" | wc -l)

cat > cities.txt <<'EOF'
67.110 36.709 CT C01
66.897 36.758 RM C03
67.686 36.698 LB C02
68.015 36.265 LM C04
68.860 36.728 LM C06
65.752 36.667 LM C05
67.278 37.224 LM C07
EOF
read MW MH <<< "$(gmt mapproject -R$R -J$J -W)"
GMTV=$(gmt --version); GSV=$(gs --version)
WO='=~1p,white'

overlay() {
  gmt coast -Di -W0.3p,gray25 -N1/0.9p,gray10,- -I1/0.6p,70/120/190 -A20 -Bxa1f10mg30m -Bya30mf10mg30m -BWSen
  gmt plot struct.gmt -W0.4p,gray20
  gmt plot qf.gmt -W1.5p,white; gmt plot qf.gmt -W0.7p,black
  grep -v '^#' well.gmt | grep -v '^>' | gmt plot -Sc0.17c -Gwhite -W0.7p,black
  echo "$EQLON $EQLAT" | gmt plot -Sa0.5c -Gwhite -W0.9p,black
  awk '{print $1,$2}' cities.txt | gmt plot -Ss0.2c -Gwhite -W0.8p,black
  while read lo la j nm; do
    case $j in RM) dx=-0.15c; dy=0;; CT) dx=0; dy=-0.16c;; LB) dx=0.15c; dy=0.1c;; *) dx=0.15c; dy=0;; esac
    echo "$lo $la $nm" | gmt text -F+f8p,Helvetica,black$WO+j$j -D$dx/$dy
  done < cities.txt
  echo "67.03 36.51 -1 Alburz-Marmul Fault" | gmt text -F+f8p,Helvetica-Oblique,black$WO+a+jCM
}

gmt begin $OUT/geomorphic_indices pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 8p,Helvetica \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray50 \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel MAP_TICK_LENGTH_PRIMARY 4p PS_PAGE_COLOR white
gmt grdimage relief.nc -Iint.nc -Cgray.cpt -R$R -J$J -Y5.4c
gmt plot ksn.txt -Ss0.03c -Cksn.cpt
overlay
echo "(b) Normalised channel steepness k@-sn@- (@~q@~@-ref@- = 0.45)" | gmt text -F+cTL+f9p,Helvetica,black=~1.2p,white -D0.25c/-0.25c
gmt colorbar -Cksn.cpt -DJRM+w${MH}c/0.32c+o0.3c/0+e -Q -Bxa1pf3+l"k@-sn@- (m@+0.9@+)" --FONT_ANNOT_PRIMARY=8p,Helvetica
gmt basemap -LjBR+c36.6+w50k+f+l"km"+o0.4c/0.4c -F+gwhite@15+p0.4p,gray40+r0.1c
cat > leg.txt <<EOF
N 1 0.95 0.8 1.25
S 0.3c - 0.55c - 0.4p,gray20 0.75c Subsurface structure
S 0.3c - 0.55c - 1.0p,black 0.75c Quaternary fault
S 0.3c c 0.17c white 0.7p,black 0.75c Well with logs
S 0.3c a 0.35c white 0.9p,black 0.75c 2025 mainshock
EOF
gmt legend leg.txt -DJBC+w${MW}c+o0/1.0c -F+gwhite+p0.4p,gray40 --FONT_ANNOT_PRIMARY=8.5p,Helvetica
gmt text -R0/$MW/0/2 -JX${MW}c/2c -M -N -F+f8p,Helvetica,gray20+jTL -Y-3.6c <<EOF
> 0 1.75 0.35c ${MW}c j
Depicts: (a) hypsometric integral of drainage basins (GRASS r.watershed, threshold 4000 cells) and (b) normalised channel steepness index k@-sn@- = S A@+0.45@+ on $NKS channel cells (r.stream.extract, threshold 300 cells; slope smoothed over 5 x 5 cells), both on SRTM 3-arc-second shaded relief with petroleum-exploration structures, Quaternary faults, wells and the 2025 mainshock; Lambert conformal conic projection (central meridian 67.5@.E, standard parallels 36@.N and 37.3@.N). Data: SRTM 3 arc-sec; USGS petroleum GIS (2006); USGS OFR 2007-1103; GSHHG. Software: GRASS GIS 8.4.2, GMT $GMTV, GhostScript $GSV. Source: authors.
EOF
SH=$(awk -v h=$MH 'BEGIN{print 3.6 + h + 0.9}')
gmt grdimage relief.nc -Iint.nc -Cgray.cpt -R$R -J$J -Y${SH}c
gmt grdimage hi.nc -Iinth.nc -Chi.cpt -Q -t15
overlay
echo "(a) Hypsometric integral of drainage basins" | gmt text -F+cTL+f9p,Helvetica,black=~1.2p,white -D0.25c/-0.25c
gmt colorbar -Chi.cpt -DJRM+w${MH}c/0.32c+o0.3c/0+e -Bxa0.1f0.01+l"Hypsometric integral" --FONT_ANNOT_PRIMARY=8p,Helvetica
gmt basemap -TdjTR+w1.0c+l,,,N+o0.4c/0.4c
gmt end
echo "done NKS=$NKS"
