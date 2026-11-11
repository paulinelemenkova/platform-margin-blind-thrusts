#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
mkdir -p work_sm && cd work_sm

W=64.5; E=71.5; S=35.0; N=38.0
R=$W/$S/$E/$N+r; J=L68/36.5/35.5/37.5/17c
BW=65.6; BE=69.4; BS=35.85; BN=37.45
DMAX=70

curl -s -m 180 "https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=1973-01-01&endtime=2026-11-10&minlatitude=34.5&maxlatitude=38.5&minlongitude=64&maxlongitude=72&minmagnitude=4&orderby=magnitude-asc" > eq.csv
awk -F, 'NR>1 && $4!="" {printf "%s %s %s %.3f\n", $3, $2, $4, 0.03*1.65^($5-3.5)+0.04}' eq.csv > eq.txt
NEQ=$(wc -l < eq.txt)
B=https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog
curl -s -m 200 $B/jan76_dec20.ndk > cmt.ndk
for y in 2021 2022 2023 2024 2025; do yy=${y:2:2}
  for m in jan feb mar apr may jun jul aug sep oct nov dec; do
    curl -s -f -m 30 $B/NEW_MONTHLY/$y/$m$yy.ndk >> cmt.ndk || true; done; done
curl -s -m 60 $B/NEW_QUICK/qcmt.ndk >> cmt.ndk
python3 - cmt.ndk > cmt_all.txt <<'PY'
import sys, math
L = [l.rstrip("\n") for l in open(sys.argv[1]) if l.strip()]
seen = set()
for i in range(0, len(L) - 4, 5):
    a, b, c, d, e = L[i:i + 5]
    try:
        cid = b.split()[0]
        t = c.split(); lat, lon, dep = float(t[3]), float(t[5]), float(t[7])
        ex = int(d.split()[0]); v = e.split()
        m0 = float(v[10]) * 10 ** ex; s1, d1, r1 = int(v[11]), int(v[12]), int(v[13])
    except (ValueError, IndexError):
        continue
    if cid in seen: continue
    seen.add(cid)
    if 64 <= lon <= 72 and 34.5 <= lat <= 38.5:
        mw = 2 / 3 * (math.log10(m0) - 16.1)
        print(f"{lon} {lat} {dep} {s1} {d1} {r1} {mw:.2f} {a[5:15].replace('/', '-')}")
PY
sort -k7,7n cmt_all.txt | awk -v d=$DMAX '$3<=d' > cmt_crust_d.txt
cut -d' ' -f1-7 cmt_crust_d.txt > cmt_crust.txt
NCMT=$(wc -l < cmt_all.txt); NCR=$(wc -l < cmt_crust.txt)
grep "2025-11-02" cmt_crust_d.txt | cut -d' ' -f1-7 > ev.txt || true
echo "eq=$NEQ cmt=$NCMT crustal=$NCR"; cat ev.txt

curl -s -m 180 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
QSHP=$(find . -name 'faults.shp' | head -1)
ogr2ogr -f OGR_GMT qf.gmt "$QSHP" -t_srs EPSG:4326

gmt grdcut @earth_relief_15s -R63.5/72.5/34/39 -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.7 -Gint.nc
gmt makecpt -Cgray -T-3000/9000 > gray.cpt
gmt makecpt -Cturbo -T0/250/5 > depth.cpt

cat > cities.txt <<'EOF'
67.110 36.709 CB C01
68.860 36.728 LM C06
65.752 36.667 LM C05
69.181 34.519 LM C10
70.580 37.120 LM C15
67.278 37.224 LM C07
68.780 38.560 LM Dushanbe
EOF
box() { python3 -c "
W,E,S,N=$BW,$BE,$BS,$BN; n=60
for i in range(n+1): print(W+(E-W)*i/n,S)
for i in range(n+1): print(E,S+(N-S)*i/n)
for i in range(n+1): print(E-(E-W)*i/n,N)
for i in range(n+1): print(W,N-(N-S)*i/n)"; }
box > box.txt
read MW MH <<< "$(gmt mapproject -R$R -J$J -W)"
GMTV=$(gmt --version); GSV=$(gs --version)
WO='=~1p,white'

gmt begin $OUT/seismicity_mechanisms pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 8p,Helvetica \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray50 \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel MAP_TICK_LENGTH_PRIMARY 4p PS_PAGE_COLOR white
gmt grdimage relief.nc -Iint.nc -Cgray.cpt -R$R -J$J -Y8.4c
gmt coast -Di -W0.35p,gray25 -N1/0.9p,gray10,- -I1/0.6p,70/120/190 -C150/180/220 -A20 \
          -Bxa1f15mg1 -Bya1f15mg1 -BWSen
gmt plot qf.gmt -W0.6p,60/60/60
gmt plot box.txt -W2p,white; gmt plot box.txt -W1.1p,200/0/0,6_3:0p
gmt plot eq.txt -Sc -Cdepth.cpt -W0.25p,gray15
gmt meca cmt_crust.txt -Sa0.33c -Cdepth.cpt -Ewhite -W0.4p,black
[ -s ev.txt ] && awk '{print $1,$2,$3,$4,$5,$6,$7,66.35,36.22}' ev.txt | \
   gmt meca -Sa0.85c+m -Cdepth.cpt -Ewhite -W0.8p,black -A+p0.8p,black
echo "66.35 35.99 2025, M@-w@- 6.3 (GCMT)" | gmt text -F+f8p,Helvetica,black$WO+jTC
awk '{print $1,$2}' cities.txt | gmt plot -Ss0.2c -Gwhite -W0.8p,black
while read lo la j nm; do
  case $j in RM) dx=-0.15c; dy=0;; CB) dx=0; dy=0.16c;; *) dx=0.15c; dy=0;; esac
  echo "$lo $la $nm" | gmt text -F+f8p,Helvetica,black$WO+j$j -D$dx/$dy
done < cities.txt
printf "70.9 35.55 0 Hindu Kush\n71.15 37.85 0 Pamir\n65.3 37.75 0 Amu Darya Basin\n66.6 35.47 0 Northern Platform\n" | \
   gmt text -F+f9p,Helvetica-Oblique,black$WO+a+jCM
echo "$BE $BN Study area" | gmt text -F+f8p,Helvetica,200/0/0$WO+jBR -D0/0.08c
gmt basemap -LjBL+c36.5+w100k+f+l"km"+o0.4c/0.4c -F+gwhite@15+p0.4p,gray40+r0.1c
gmt basemap -TdjTL+w1.0c+l,,,N+o0.4c/0.4c
gmt colorbar -Cdepth.cpt -DJBC+w${MW}c/0.32c+o0/0.85c+h+ef -Bxa50f10+l"Hypocentral / centroid depth (km)"
for m in 4 5 6 7; do echo $m; done | awk '{printf "S 0.35c c %.3f gray80 0.25p,gray15 0.9c M %s\n", 0.03*1.65^($1-3.5)+0.04, $1}' > leg_m.txt
cat > leg.txt <<EOF
N 4
$(cat leg_m.txt)
S 0.35c - 0.6c - 0.6p,60/60/60 0.9c Quaternary fault
S 0.35c - 0.6c - 1.1p,200/0/0,6_3:0p 0.9c Study area
S 0.35c s 0.2c white 0.8p,black 0.9c City
S 0.35c - 0.6c - 0.9p,gray10,- 0.9c International border
G 0.1c
N 1
L 8.5p,Helvetica L Beachballs: Global CMT mechanisms of crustal events (centroid depth @~\243@~ $DMAX km), compressional quadrants coloured by depth
EOF
gmt legend leg.txt -DJBC+w${MW}c+o0/2.45c -F+gwhite+p0.4p,gray40 --FONT_ANNOT_PRIMARY=8.5p,Helvetica
gmt text -R0/$MW/0/2 -JX${MW}c/2c -M -N -F+f8p,Helvetica,gray20+jTL -Y-5.35c <<EOF
> 0 1.0 0.35c ${MW}c j
Depicts: $NEQ earthquakes of M@~\263@~4 (USGS ComCat, 1973-2026) coloured by hypocentral depth with the turbo palette and sized by magnitude, and $NCR Global CMT focal mechanisms of crustal events (centroid depth @~\243@~$DMAX km, 1976-2026) coloured by centroid depth, on SRTM15+ shaded relief with Quaternary faults; the 2025 mechanism is enlarged and offset. Lambert conformal conic projection (central meridian 68@.E, standard parallels 35.5@.N and 37.5@.N). Data: USGS ComCat; Global CMT Project; USGS OFR 2007-1103; SRTM15+; GSHHG. Software: GMT $GMTV, GhostScript $GSV, Python 3. Source: authors.
EOF
gmt end
echo done
