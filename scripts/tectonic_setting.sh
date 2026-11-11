#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
GDRIVE_ID=${GDRIVE_ID:?}
mkdir -p work && cd work

RA=59.2/28.6/76.6/39.4+r;   JA=L67.9/34/30/38/14.8c
BW=65.6; BE=69.4; BS=35.85; BN=37.45
RB=$BW/$BS/$BE/$BN+r;       JB=L67.5/36.65/36/37.3/17c
GAP=1.3

curl -s -f -m 600 -L -o src.zip "https://drive.usercontent.google.com/download?id=$GDRIVE_ID&export=download&confirm=t"
unzip -q -o src.zip -x '__MACOSX/*'
cpshp() { P=$(find . -name "$1.shp" | head -1); for x in shp shx dbf prj; do [ -f "${P%.shp}.$x" ] && cp "${P%.shp}.$x" "$2.$x"; done; return 0; }
cpshp fltafg flt; cpshp structafg struct; cpshp fldafg fld
cpshp geoBoundaries-AFG-ADM0_simplified adm0; cpshp geoBoundaries-AFG-ADM1_simplified adm1
curl -s -m 180 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
QSHP=$(find . -name 'faults.shp' | head -1)

g() { S=""; [ -f "${2%.shp}.prj" ] || S="-s_srs EPSG:4326"; ogr2ogr -f OGR_GMT "$1" "$2" $S -t_srs EPSG:4326 "${@:3}"; }
g adm0.gmt adm0.shp; g adm1.gmt adm1.shp
g qR.gmt  "$QSHP" -where "Type='R'"
g qSS.gmt "$QSHP" -where "Type IN ('RL','LL','SS')"
g qN.gmt  "$QSHP" -where "Type='N'"
g am.gmt  "$QSHP" -where "Name='Alburz-Marmul fault'"
g fltP.gmt flt.shp -where "TYPE NOT LIKE '%buried%' AND TYPE NOT LIKE '%inferred%'"
g fltB.gmt flt.shp -where "TYPE LIKE '%buried%' OR TYPE LIKE '%inferred%'"
g struct.gmt struct.shp; g fld.gmt fld.shp
fc() { ogrinfo -so -al "$@" | awk -F': ' '/Feature Count/{print $2}'; }
NSTR=$(fc struct.shp); NFLD=$(fc fld.shp)

curl -s -m 120 "https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=1973-01-01&endtime=2026-11-10&minlatitude=$BS&maxlatitude=$BN&minlongitude=$BW&maxlongitude=$BE&minmagnitude=3.5&orderby=time-asc" > eq.csv
sz='function s(m){return 0.035*2^(m-3)+0.05}'
awk -F, "$sz"' NR>1 && $5<6 && $1<"2025-11-02T20:29"{print $3,$2,s($5)}' eq.csv > eq_old.txt
awk -F, "$sz"' NR>1 && $5<6 && $1>="2025-11-02T20:29"{print $3,$2,s($5)}' eq.csv > eq_new.txt
NOLD=$(wc -l < eq_old.txt); NNEW=$(wc -l < eq_new.txt)
EQLON=67.3886; EQLAT=36.5461
MT=$(curl -s -m 60 "https://earthquake.usgs.gov/fdsnws/event/1/query?eventid=us6000rl31&format=geojson" | python3 -c "
import json,sys; p=json.load(sys.stdin)['properties']['products']['moment-tensor']
m=[x['properties'] for x in p if x['properties'].get('derived-magnitude-type')=='Mww'][0]
print(m['nodal-plane-1-strike'],m['nodal-plane-1-dip'],m['nodal-plane-1-rake'],m['derived-depth'])")
read STK DIP RAK MTD <<< "$MT"
echo "old $NOLD new $NNEW MT $STK $DIP $RAK $MTD"

gmt grdcut @earth_relief_15s -R57.5/78/27.5/41 -Grelief_a.nc
gmt grdgradient relief_a.nc -A315 -Nt0.55 -Gint_a.nc
gmt grdcut @earth_relief_03s -R65.5/69.5/35.75/37.55 -Grelief_b.nc
gmt grdgradient relief_b.nc -A315 -Nt0.6 -Gint_b.nc
gmt makecpt -Cetopo1 -G0/8000 -T0/8000 > topo.cpt

read WA HA <<< "$(gmt mapproject -R$RA -J$JA -W)"; read WB HB <<< "$(gmt mapproject -R$RB -J$JB -W)"
read AMX AMY AMA <<< "$(grep -v '^[#>]' am.gmt | gmt mapproject -R$RB -J$JB | python3 -c "
import sys,math,numpy as n
d=n.loadtxt(sys.stdin); d=d[n.argsort(d[:,0])]; k=n.polyfit(d[:,0],d[:,1],1)
xc=n.median(d[:,0]); print(xc, n.polyval(k,xc)-0.32, math.degrees(math.atan(k[0])))")"
read AMLON AMLAT <<< "$(echo $AMX $AMY | gmt mapproject -R$RB -J$JB -I)"

cat > cities_a.txt <<'EOF'
69.181 34.519 LM C10
62.170 34.330 LM C11
65.695 31.610 LM C12
67.110 36.709 RM C01
70.436 34.442 LM C13
68.860 36.728 LM C06
70.580 37.120 LM C15
68.420 33.550 LM C14
EOF
cat > cities_b.txt <<'EOF'
67.110 36.709 TC C01
66.897 36.758 RM C03
67.686 36.698 LB C02
68.015 36.265 LM C04
68.860 36.728 LM C06
65.752 36.667 LM C05
67.278 37.224 LM C07
68.708 35.944 LM C08
68.712 36.131 LM C09
EOF
cat > blocks.txt <<'EOF'
64.0 35.55 0 Northern Platform
61.75 32.4 0 Farah Block
64.6 31.0 0 Helmand Block
68.25 32.25 30 Katawaz Basin
69.05 34.95 0 Central Block
70.2 35.95 28 Hindu Kush
73.6 38.3 0 Pamir
63.6 38.3 0 Amu Darya Basin
68.75 38.1 0 Eastern Foreland Basin
EOF
cat > countries.txt <<'EOF'
EOF
cat > faults_lab.txt <<'EOF'
64.6 34.55 -3 Hari Rud Fault
67.25 30.7 58 Chaman Fault
EOF
box() { python3 -c "
W,E,S,N=$BW,$BE,$BS,$BN; n=60
for i in range(n+1): print(W+(E-W)*i/n,S)
for i in range(n+1): print(E,S+(N-S)*i/n)
for i in range(n+1): print(E-(E-W)*i/n,N)
for i in range(n+1): print(W,N-(N-S)*i/n)"; }
box > box.txt

GMTV=$(gmt --version); GDALV=$(gdal-config --version 2>/dev/null || ogr2ogr --version | awk '{print $2}' | tr -d ,); GSV=$(gs --version)
YB=5.75

gmt begin $OUT/tectonic_setting pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 10p,Helvetica-Bold \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray50 \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel MAP_TICK_LENGTH_PRIMARY 4p \
        PS_PAGE_COLOR white

gmt grdimage relief_b.nc -Iint_b.nc -Ctopo.cpt -R$RB -J$JB -Y${YB}c
gmt plot struct.gmt -G240/215/160@45 -W0.35p,150/80/30
gmt plot fld.gmt -G0/150/130@25 -W0.45p,0/90/80
gmt coast -Di -W0.3p,gray25 -N1/0.9p,gray10,- -I1/0.6p,30/110/200 -I2/0.4p,30/110/200 -C90/160/230 -A20 \
          -Bxa1f10mg30m -Bya30mf10mg30m -BWSen
gmt plot adm1.gmt -W0.4p,gray20,4_2:0p
gmt plot fltP.gmt -W0.45p,gray20
gmt plot fltB.gmt -W0.45p,gray20,3_2:0p
for k in "qN:0/70/255" "qSS:black" "qR:220/20/60"; do f=${k%%:*}; c=${k#*:}
  gmt plot $f.gmt -W1.9p,white; gmt plot $f.gmt -W1.0p,$c; done
gmt plot eq_old.txt -Sc -G255/255/255@10 -W0.4p,gray15
gmt plot eq_new.txt -Sc -G255/140/0 -W0.4p,black
echo "$EQLON $EQLAT" | gmt plot -Sa0.55c -Gred -W0.6p,black
echo "$EQLON $EQLAT $MTD $STK $DIP $RAK 6.2 66.80 36.28" | gmt meca -Sa1.1c+m -G200/30/30 -Ewhite -W0.6p,black -A+p0.6p,black
echo "66.80 36.08 2025, M@-ww@- 6.2" | gmt text -F+f8p,Helvetica,black=~0.9p,white+jTC
echo "$AMLON $AMLAT $AMA Alburz-Marmul Fault" | gmt text -F+f8p,Helvetica-Oblique,220/20/60=~0.9p,white+a+jCM
echo "66.4 37.27 Amu Darya" | gmt text -F+f8p,Helvetica-Oblique,30/90/180=~0.9p,white+a-12+jCM
awk '{print $1,$2}' cities_b.txt | gmt plot -Ss0.22c -Gwhite -W0.8p,black
while read lo la j nm; do
  case $j in RM|RB) dx=-0.16c;; TC) dx=0;; *) dx=0.16c;; esac; case $j in *B) dy=0.12c;; TC) dy=-0.2c;; *) dy=0;; esac
  echo "$lo $la $nm" | gmt text -F+f8p,Helvetica,black=~0.9p,white+j$j -D$dx/$dy
done < cities_b.txt
echo "(b)" | gmt text -F+cTL+f11p,Helvetica-Bold,black=~1.2p,white -D0.25c/-0.25c
gmt basemap -LjBL+c36.6+w50k+f+l"km"+o0.4c/0.45c -F+gwhite@15+p0.4p,gray40+r0.1c
gmt basemap -TdjTL+w1.0c+l,,,N+o0.4c/1.15c --FONT_TITLE=8p,Helvetica
cat > leg.txt <<EOF
N 3
S 0.3c - 0.5c - 1.2p,220/20/60 0.7c Quaternary reverse fault
S 0.3c - 0.5c - 1.2p,black 0.7c Quaternary strike-slip fault
S 0.3c - 0.5c - 1.2p,0/70/255 0.7c Quaternary normal fault
S 0.3c - 0.5c - 0.5p,gray20 0.7c Mapped fault, proven
S 0.3c - 0.5c - 0.5p,gray20,3_2:0p 0.7c Mapped fault, buried/inferred
S 0.3c r 0.45c/0.3c 240/215/160 0.35p,150/80/30 0.7c Subsurface structure ($NSTR)
S 0.3c r 0.45c/0.3c 0/150/130 0.45p,0/90/80 0.7c Oil and gas field ($NFLD)
S 0.3c c 0.2c white 0.4p,gray15 0.7c Earthquake 1973-2025 ($NOLD)
S 0.3c c 0.2c 255/140/0 0.4p,black 0.7c 2025 sequence ($NNEW)
S 0.3c a 0.35c red 0.6p,black 0.7c 2025 mainshock (USGS)
S 0.3c s 0.22c white 0.8p,black 0.7c City
S 0.3c - 0.5c - 0.9p,gray10,- 0.7c International border
EOF
gmt legend leg.txt -DJBC+w17c+o0/1.0c -F+gwhite+p0.4p,gray40 --FONT_ANNOT_PRIMARY=8.5p,Helvetica
gmt text -R0/17/0/2.1 -JX17c/2.1c -M -N -F+f8p,Helvetica,gray20+jTL -Y-5.55c <<EOF
> 0 2.1 0.35c 17c j
Depicts: (a) tectonic blocks, Quaternary faults and the study area on ETOPO1-coloured SRTM15+ relief, Lambert conformal conic projection (central meridian 67.9@.E, standard parallels 30@.N and 38@.N); (b) study area on SRTM 3-arc-second relief, Lambert conformal conic projection (standard parallels 36@.N and 37.3@.N), with $NSTR petroleum-exploration structures, $NFLD oil and gas fields, mapped and Quaternary faults, ComCat seismicity (M@~\263@~3.5) and the 2025 moment tensor (USGS Mww). Data: USGS petroleum and geology GIS (2006); USGS OFR 2007-1103; USGS ComCat; geoBoundaries; SRTM15+; SRTM 3 arc-sec; GSHHG. Software: GMT $GMTV, GDAL $GDALV, GhostScript $GSV, Python 3. Source: authors.
EOF

SHIFT=$(awk -v h=$HB -v g=$GAP 'BEGIN{print 5.55+h+g}')
gmt grdimage relief_a.nc -Iint_a.nc -Ctopo.cpt -R$RA -J$JA -Y${SHIFT}c
gmt clip adm0.gmt -N
printf "55 25\n80 25\n80 43\n55 43\n" | gmt plot -Gwhite -t42
gmt clip -C
gmt coast -Di -W0.3p,gray25 -N1/0.7p,gray15,- -I1/0.45p,30/110/200 -C90/160/230 -A50 \
          -Bxa2f1g1 -Bya2f1g1 -BWSen
gmt plot adm1.gmt -W0.3p,gray30,4_2:0p
for k in "qN:0/70/255" "qSS:black" "qR:220/20/60"; do f=${k%%:*}; c=${k#*:}
  gmt plot $f.gmt -W1.4p,white; gmt plot $f.gmt -W0.7p,$c; done
gmt plot adm0.gmt -W1.6p,black
gmt plot box.txt -W2.4p,white; gmt plot box.txt -W1.2p,200/0/0
echo "$BE $BN (b)" | gmt text -F+f9p,Helvetica-Bold,200/0/0=~1p,white+jBL -D0.08c/0.06c
echo "$EQLON $EQLAT" | gmt plot -Sa0.4c -Gred -W0.5p,black
gmt text blocks.txt -F+f9p,Helvetica-Oblique,black=~1p,white+a+jCM
gmt text faults_lab.txt -F+f8p,Helvetica-Oblique,220/20/60=~0.9p,white+a+jCM
awk '{print $1,$2}' cities_a.txt | gmt plot -Ss0.2c -Gwhite -W0.8p,black
grep -vw C01 cities_a.txt | awk '{$3="";print}' | gmt text -F+f8p,Helvetica,black=~0.9p,white+jLM -D0.16c/0
grep -w C01 cities_a.txt | awk '{$3="";print}' | gmt text -F+f8p,Helvetica,black=~0.9p,white+jRM -D-0.16c/0
echo "(a)" | gmt text -F+cTL+f11p,Helvetica-Bold,black=~1.2p,white -D0.25c/-0.25c
gmt basemap -LjBR+c34+w200k+f+l"km"+o2.6c/0.45c -F+gwhite@15+p0.4p,gray40+r0.1c
gmt basemap -TdjRM+w1.1c+l,,,N+o0.5c/0 --FONT_TITLE=8p,Helvetica
gmt colorbar -Ctopo.cpt -DJRM+w${HA}c/0.35c+o0.35c/0+ef -Bxa1000f100+l"Elevation (m)" \
          --FONT_ANNOT_PRIMARY=8p,Helvetica --FONT_LABEL=9p,Helvetica

read X1 Y1 <<< "$(echo $BW $BS | gmt mapproject -R$RA -J$JA -Dc)"
read X2 Y2 <<< "$(echo $BE $BS | gmt mapproject -R$RA -J$JA -Dc)"
printf "> \n$X1 $Y1\n0 -$GAP\n> \n$X2 $Y2\n$WB -$GAP\n" | gmt plot -R0/17/0/1 -Jx1c -N -W0.6p,200/0/0,4_2:0p
gmt end
echo done
