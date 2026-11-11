#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
mkdir -p work_ds && cd work_ds

EQLON=67.3886; EQLAT=36.5461; EQDEP=28
AZ=340; HL=70; HW=40; DMAX=70
NP="241/70/37 137/55/156"
R=66.62/35.80/68.16/37.30; J=L67.4/36.55/36/37.1/6.2c

cat > pq.py <<'PY'
import sys, math
LON0, LAT0, AZ = 67.3886, 36.5461, math.radians(340)
KX = 111.32 * math.cos(math.radians(LAT0)); KY = 111.195
def pq(lon, lat):
    x, y = (lon - LON0) * KX, (lat - LAT0) * KY
    return x * math.sin(AZ) + y * math.cos(AZ), x * math.cos(AZ) - y * math.sin(AZ)
def ll(p, q):
    x = p * math.sin(AZ) + q * math.cos(AZ); y = p * math.cos(AZ) - q * math.sin(AZ)
    return LON0 + x / KX, LAT0 + y / KY
PY
endpt() { python3 -c "import sys; sys.path.insert(0,'.'); from pq import ll; print('%.5f %.5f' % ll($1, $2))"; }
read PLON PLAT <<< "$(endpt -$HL 0)"; read QLON QLAT <<< "$(endpt $HL 0)"
python3 -c "
import sys; sys.path.insert(0,'.'); from pq import ll
pts=[(-$HL,-$HW),($HL,-$HW),($HL,$HW),(-$HL,$HW),(-$HL,-$HW)]
for (p1,q1),(p2,q2) in zip(pts[:-1],pts[1:]):
    for i in range(41):
        t=i/40; print('%.5f %.5f' % ll(p1+(p2-p1)*t, q1+(q2-q1)*t))" > swath_box.txt
python3 -c "
import sys; sys.path.insert(0,'.'); from pq import ll
for d in range(0, 2*$HL+1, 20): print('%.5f %.5f %d' % (*ll(d-$HL, 0), d))" > ticks.txt

curl -s -m 180 "https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=1973-01-01&endtime=2026-11-10&minlatitude=35.5&maxlatitude=37.6&minlongitude=66.2&maxlongitude=68.6&minmagnitude=2.5&maxdepth=$DMAX&orderby=magnitude-asc" > eq.csv
python3 - eq.csv > eq_all.txt <<'PY'
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    d = float(r["depth"]); fx = 1 if d in (10.0, 33.0, 35.0) else 0
    print(r["longitude"], r["latitude"], d, r["mag"], fx)
PY
python3 - eq_all.txt $HL $HW > eq_sec.txt <<'PY'
import sys; sys.path.insert(0, '.'); from pq import pq
HL, HW = float(sys.argv[2]), float(sys.argv[3])
for l in open(sys.argv[1]):
    lon, lat, d, m, fx = l.split(); p, q = pq(float(lon), float(lat))
    if abs(p) <= HL and abs(q) <= HW:
        s = 0.03 * 1.65 ** (float(m) - 2.5) + 0.06
        print(f"{p + HL:.2f} {d} {d} {s:.3f} {fx} {m}")
PY
awk '$5==0{print $1,$2,$3,$4}' eq_sec.txt > eq_free.txt; awk '$5==1{print $1,$2,$3,$4}' eq_sec.txt > eq_fix.txt
awk '{s=0.03*1.65^($4-2.5)+0.06; print $1,$2,$3,s}' eq_all.txt > eq_map.txt
NSEC=$(wc -l < eq_sec.txt); NFIX=$(wc -l < eq_fix.txt); NALL=$(wc -l < eq_all.txt)

B=https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog
curl -s -m 200 $B/jan76_dec20.ndk > cmt.ndk
for y in 2021 2022 2023 2024 2025; do yy=${y:2:2}
  for m in jan feb mar apr may jun jul aug sep oct nov dec; do
    curl -s -f -m 30 $B/NEW_MONTHLY/$y/$m$yy.ndk >> cmt.ndk || true; done; done
curl -s -m 60 $B/NEW_QUICK/qcmt.ndk >> cmt.ndk
python3 - cmt.ndk $HL $HW $DMAX > cmt.txt <<'PY'
import sys, math; sys.path.insert(0, '.'); from pq import pq
HL, HW, DM = map(float, sys.argv[2:5])
L = [l.rstrip("\n") for l in open(sys.argv[1]) if l.strip()]
seen = set()
for i in range(0, len(L) - 4, 5):
    a, b, c, d, e = L[i:i + 5]
    try:
        cid = b.split()[0]; t = c.split(); lat, lon, dep = float(t[3]), float(t[5]), float(t[7])
        ex = int(d.split()[0]); v = e.split(); m0 = float(v[10]) * 10 ** ex
        s1, d1, r1 = int(v[11]), int(v[12]), int(v[13])
    except (ValueError, IndexError):
        continue
    if cid in seen: continue
    seen.add(cid)
    if 66.2 <= lon <= 68.6 and 35.5 <= lat <= 37.6 and dep <= DM:
        p, q = pq(lon, lat); ins = int(abs(p) <= HL and abs(q) <= HW)
        mw = 2 / 3 * (math.log10(m0) - 16.1)
        print(f"{lon} {lat} {dep} {s1} {d1} {r1} {mw:.2f} {ins} {a[5:15]}")
PY
awk '$8==1 && $9!="2025/11/02"{print $1,$2,$3,$4,$5,$6,$7}' cmt.txt > cmt_sec.txt
awk '$9=="2025/11/02"{print $1,$2,$3,$4,$5,$6,$7}' cmt.txt | sort -k7,7nr | head -1 > ev.txt
awk '$9!="2025/11/02"{print $1,$2,$3,$4,$5,$6,$7}' cmt.txt > cmt_map.txt
NCMT=$(wc -l < cmt_sec.txt); echo "eq_all=$NALL sec=$NSEC fixed=$NFIX cmt_sec=$NCMT"; cat ev.txt; cat cmt_sec.txt

curl -s -m 180 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
QSHP=$(find . -name 'faults.shp' | head -1)
g() { s=""; [ -f "${2%.shp}.prj" ] || s="-s_srs EPSG:4326"; ogr2ogr -f OGR_GMT "$1" "$2" $s -t_srs EPSG:4326 "${@:3}"; }
g qf.gmt "$QSHP"; g amf.gmt "$QSHP" -where "Name LIKE '%Alburz%'"
for f in qf amf; do
  awk 'NF>=2 && $1+0==$1 {print $1,$2}' $f.gmt > ${f}_v.txt
  python3 - ${f}_v.txt $HL > ${f}_x.txt <<'PY'
import sys; sys.path.insert(0, '.'); from pq import pq
HL = float(sys.argv[2]); P = []
for l in open(sys.argv[1]):
    lon, lat = map(float, l.split()[:2]); p, q = pq(lon, lat)
    if abs(q) <= 1.5 and abs(p) <= HL: P.append(p + HL)
P.sort(); G = []
for p in P:
    if G and p - G[-1][-1] <= 4: G[-1].append(p)
    else: G.append([p])
for g in G: print(f"{sum(g) / len(g):.2f}")
PY
done
python3 - amf_v.txt $HL $HW > amf_p.txt <<'PY'
import sys; sys.path.insert(0, '.'); from pq import pq
HL, HW = float(sys.argv[2]), float(sys.argv[3]); P = []
for l in open(sys.argv[1]):
    lon, lat = map(float, l.split()[:2]); p, q = pq(lon, lat)
    if abs(q) <= HW and abs(p) <= HL: P.append(p + HL)
if P: print(f"{min(P):.1f} {sum(P) / len(P):.1f} {max(P):.1f}")
PY
read AMIN AMEAN AMAX < amf_p.txt || true
echo "fault crossings:"; cat qf_x.txt | tr '\n' ' '; echo; echo "AMF projected: $AMIN $AMEAN $AMAX"

gmt grdcut @earth_relief_03s -R66.4/68.4/35.6/37.5 -Grelief.nc || { sleep 20; gmt grdcut @earth_relief_03s -R66.4/68.4/35.6/37.5 -Grelief.nc; }
gmt grdgradient relief.nc -A315 -Nt0.65 -Gint.nc
python3 -c "
import sys; sys.path.insert(0,'.'); from pq import ll
p=-$HL
while p<=$HL+1e-9:
    for q in range(-10,11): print('%.6f %.6f %.2f' % (*ll(p,q), p+$HL))
    p+=0.5" > sw.pts
gmt grdtrack sw.pts -Grelief.nc > sw.z
python3 - sw.z > sw_stat.txt <<'PY'
import sys, collections
D = collections.defaultdict(list)
for l in open(sys.argv[1]):
    t = l.split(); D[float(t[2])].append(float(t[3]))
for d in sorted(D): v = D[d]; print(f"{d} {min(v):.0f} {sum(v)/len(v):.0f} {max(v):.0f}")
PY
awk '{print $1,$4}' sw_stat.txt > env.txt; tac sw_stat.txt | awk '{print $1,$2}' >> env.txt
ZMAX=$(awk 'BEGIN{m=0}{if($4>m)m=$4}END{print (int(m/500)+1)*500}' sw_stat.txt)
EQTOPO=$(awk -v x=$HL '$1==x{print $3}' sw_stat.txt)

gmt makecpt -Cgray -T-2500/6000 > gray.cpt
gmt makecpt -Cjet -T0/$DMAX > depth.cpt

python3 - $AZ $HL $EQDEP $NP > np.txt <<'PY'
import sys, math
az, HL, z0 = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
for k, s in enumerate(sys.argv[4:], 1):
    st, dp, rk = map(float, s.split('/'))
    beta = math.radians((st + 90) - az)
    ta = math.tan(math.radians(dp)) * math.cos(beta)
    ad = math.degrees(math.atan(abs(ta)))
    x0 = HL
    xs = x0 - z0 / ta; xe = x0 + (70 - z0) / ta
    print(f"NP{k} {st:.0f} {dp:.0f} {rk:.0f} {ad:.0f} {xs:.2f} {xe:.2f} {'NNW' if ta > 0 else 'SSE'}")
PY
cat np.txt

GMTV=$(gmt --version); GSV=$(gs --version); WO='=~1p,white'
SW=16; L2=$((2*HL)); SH=$(awk -v w=$SW -v l=$L2 -v d=$DMAX 'BEGIN{printf "%.3f", w*d/l}')
TH=1.9
gmt begin $OUT/depth_section pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 8p,Helvetica \
        MAP_FRAME_TYPE plain MAP_GRID_PEN_PRIMARY 0.25p,gray50 FORMAT_GEO_MAP ddd:mmF \
        MAP_ANNOT_OBLIQUE lat_parallel MAP_TICK_LENGTH_PRIMARY 4p PS_PAGE_COLOR white

gmt basemap -R0/$L2/0/$DMAX -JX${SW}c/-${SH}c -Bxa20f5g20+l"Distance along P-P@~\242@~ (km)" \
    -Bya10f5g10+l"Depth (km)" -BWSen -Y4.75c
for z in 10 33 35; do printf "0 $z\n$L2 $z\n" | gmt plot -W0.35p,gray60,2_2:0p; done
while read id st dp rk ad xs xe dir; do
  if [ "$id" = NP1 ]; then pen=1.3p,black; else pen=1.3p,black,6_3:0p; fi
  printf "$xs 0\n$xe $DMAX\n" | gmt plot -W2.6p,white; printf "$xs 0\n$xe $DMAX\n" | gmt plot -W$pen
  if [ "$id" = NP1 ]; then fr=0.21; else fr=0.12; fi
  lx=$(awk -v a=$xs -v b=$xe -v f=$fr 'BEGIN{print a+(b-a)*f}'); lz=$(awk -v f=$fr 'BEGIN{print 70*f}')
  rot=$(awk -v a=$xs -v b=$xe -v w=$SW -v l=$L2 -v h=$SH -v d=$DMAX 'BEGIN{r=atan2(-h, (b-a)*w/l)*57.29578; if(r<-90)r+=180; if(r>90)r-=180; print r}')
  case $dir in NNW) j=CB; dy=0.12c;; *) j=CB; dy=0.12c;; esac
  echo "$lx $lz $rot $id ${st}@./${dp}@./${rk}@." | gmt text -F+f7.5p,Helvetica,black$WO+a+jCB -D0/0.1c -N
done < np.txt
gmt plot eq_fix.txt -Ss -Cdepth.cpt -W0.3p,gray20 -t0
gmt plot eq_free.txt -Sc -Cdepth.cpt -W0.3p,gray20
if [ -s cmt_sec.txt ]; then
  gmt coupe cmt_sec.txt -Aa$PLON/$PLAT/$QLON/$QLAT+d90+w$((2*HW))+z0/$DMAX -Sa0.5c -Cdepth.cpt -Ewhite -W0.4p,black -Q -N
fi
gmt coupe ev.txt -Aa$PLON/$PLAT/$QLON/$QLAT+d90+w400+z0/$DMAX -Sa0.8c -Cdepth.cpt -Ewhite -W0.8p,black -Q -N
echo "$HL $EQDEP" | gmt plot -Sa0.45c -Gwhite -W0.9p,black
printf "$HL $EQDEP\n90 44\n" | gmt plot -W0.5p,black
printf "90.8 44 2025 C02: USGS hypocentre (star, 28 km)\n90.8 47 and GCMT centroid mechanism (27 km)\n" | gmt text -F+f8p,Helvetica,black$WO+jML
echo "(c) Depth section P-P@~\242@~ (no vertical exaggeration)" | gmt text -F+cTL+f9p,Helvetica,black=~1.2p,white -D0.2c/-0.2c
echo "P" | gmt text -F+cTL+f9p,Helvetica-Bold,black -D-0.05c/0.35c -N
echo "P@~\242@~" | gmt text -F+cTR+f9p,Helvetica-Bold,black -D0.05c/0.35c -N
echo "SSE" | gmt text -F+cBL+f8p,Helvetica,gray30 -D0.05c/-0.85c -N
echo "NNW" | gmt text -F+cBR+f8p,Helvetica,gray30 -D-0.05c/-0.85c -N

gmt text -R0/$SW/0/2 -JX${SW}c/2c -M -N -F+f8p,Helvetica,gray20+jTL -Y-4.75c <<EOF
> 0 3.3 0.35c ${SW}c j
Depicts: (a) Global CMT mechanisms and USGS ComCat epicentres (M@~\263@~2.5, depth @~\243@~$DMAX km, 1973-2026) coloured by depth (jet palette) on SRTM shaded relief with Quaternary faults, the ${L2}-km profile P-P@~\242@~ (azimuth ${AZ}@., ticks every 20 km) and its $((2*HW))-km-wide swath; (b) swath topography (20-km-wide minimum-maximum envelope and mean) with Quaternary-fault crossings and the along-profile extent of the Alburz-Marmul Fault traces within the swath; (c) $NSEC hypocentres within the swath ($NFIX at catalogue default depths of 10, 33 or 35 km, squares) and $NCMT GCMT mechanisms projected on the vertical section (back hemisphere), with the apparent dips of both USGS nodal planes of the 2025 mainshock drawn through its hypocentre. Lambert conformal conic projection (central meridian 67.4@., standard parallels 36@.N and 37.1@.N). Data: USGS ComCat; Global CMT Project; USGS OFR 2007-1103; SRTM 3 arc-sec. Software: GMT $GMTV, GDAL, GhostScript $GSV, Python 3. Source: authors.
EOF

gmt basemap -R0/$L2/0/$ZMAX -JX${SW}c/${TH}c -Bxa20f5+l" " -Bya1000f250+l"Elev. (m)" -BWsn -Y$(awk -v h=$SH 'BEGIN{print 4.75+h+0.8}')c
gmt plot env.txt -G205/205/205 -W0.25p,gray50 -L
awk '{print $1,$3}' sw_stat.txt | gmt plot -W1.1p,black
while read x; do echo "$x $ZMAX" | gmt plot -Si0.25c -G60/60/60 -W0.4p,white -N -D0/0.13c; done < qf_x.txt
if [ -n "$AMEAN" ]; then
  printf "$AMIN 2800\n$AMAX 2800\n" | gmt plot -W2.2p,white; printf "$AMIN 2800\n$AMAX 2800\n" | gmt plot -W1.2p,220/20/60
  printf "$AMIN 2800\n$AMAX 2800\n" | gmt plot -Sy0.2c -W1.2p,220/20/60
  echo "$AMAX 2800 Alburz-Marmul Fault (projected)" | gmt text -F+f8p,Helvetica-Oblique,220/20/60$WO+jTR -D0/-0.12c
fi
echo "$HL $ZMAX" | gmt plot -Sa0.4c -Gwhite -W0.9p,black -N -D0/0.15c
echo "(b) Swath topography" | gmt text -F+cTL+f9p,Helvetica,black=~1.2p,white -D0.2c/-0.12c

gmt grdimage relief.nc -Iint.nc -Cgray.cpt -R$R+r -J$J -Y$(awk -v h=$TH 'BEGIN{print h+1.0}')c
gmt coast -Di -W0.3p,gray25 -N1/0.9p,gray10,- -I1/0.6p,30/110/200 -A20 -Bxa30mf10mg30m -Bya30mf10mg30m -BWseN \
    --MAP_FRAME_TYPE=fancy --MAP_FRAME_WIDTH=3p
gmt plot qf.gmt -W0.55p,40/40/40
gmt plot amf.gmt -W2.0p,white; gmt plot amf.gmt -W1.1p,220/20/60
gmt plot swath_box.txt -W1.8p,white; gmt plot swath_box.txt -W0.8p,black,4_2:0p
printf "$PLON $PLAT\n$QLON $QLAT\n" | gmt plot -W2.2p,white; printf "$PLON $PLAT\n$QLON $QLAT\n" | gmt plot -W1.0p,black
awk '{print $1,$2}' ticks.txt | gmt plot -Sc0.11c -Gblack -W0.5p,white
gmt plot eq_map.txt -Sc -Cdepth.cpt -W0.25p,gray15
gmt meca cmt_map.txt -Sa0.38c -Cdepth.cpt -Ewhite -W0.35p,black
awk '{print $1,$2,$3,$4,$5,$6,$7}' ev.txt | gmt meca -Sa0.6c -Cdepth.cpt -Ewhite -W0.6p,black
echo "$EQLON $EQLAT" | gmt plot -Sa0.42c -Gwhite -W0.9p,black
echo "$PLON $PLAT P" | gmt text -F+f9p,Helvetica-Bold,black$WO+jTC -D0/-0.15c
echo "$QLON $QLAT P@~\242@~" | gmt text -F+f9p,Helvetica-Bold,black$WO+jBC -D0/0.15c
printf "67.110 36.709 RM C01\n67.698 36.696 LM C02\n" > cities.txt
awk '{print $1,$2}' cities.txt | gmt plot -Ss0.19c -Gwhite -W0.8p,black
while read lo la j nm; do case $j in RM) dx=-0.15c;; *) dx=0.15c;; esac
  echo "$lo $la $nm" | gmt text -F+f8p,Helvetica,255/215/0+j$j -D$dx/0; done < cities.txt
echo "(a)" | gmt text -F+cTL+f9p,Helvetica,white -D0.15c/-0.15c
gmt basemap -LjBL+c36.5+w20k+f+l"km"+o0.3c/0.3c -F+gwhite@15+p0.4p,gray40+r0.1c
read MW MH <<< "$(gmt mapproject -R$R+r -J$J -W)"

cat > leg.txt <<EOF
N 1
S 0.35c c 0.25c gray70 0.3p,gray20 0.8c Hypocentre / epicentre (size: magnitude 3, 4, 5, 6)
S 0.35c s 0.25c gray70 0.3p,gray20 0.8c Hypocentre at a catalogue default depth
S 0.35c a 0.35c white 0.9p,black 0.8c 2025 mainshock (USGS, M@-ww@- 6.2)
S 0.35c - 0.6c - 1.3p,black 0.8c 2025 nodal plane NP1 (strike/dip/rake), apparent dip
S 0.35c - 0.6c - 1.3p,black,6_3:0p 0.8c 2025 nodal plane NP2, apparent dip
S 0.35c - 0.6c - 1.1p,220/20/60 0.8c Alburz-Marmul Fault (USGS OFR 2007-1103)
S 0.35c - 0.6c - 0.55p,40/40/40 0.8c Other Quaternary faults
S 0.35c - 0.6c - 0.8p,black,4_2:0p 0.8c Swath of P-P@~\242@~ (@~\261@~$HW km)
S 0.35c i 0.25c 60/60/60 0.4p,white 0.8c Quaternary fault crossing the profile
S 0.35c - 0.6c - 1.2p,220/20/60 0.8c Along-profile extent of the Alburz-Marmul Fault in the swath
EOF
LX=$(awk -v w=$MW 'BEGIN{print w+0.45}'); LW=$(awk -v w=$MW -v s=$SW 'BEGIN{print s-w-0.45}')
gmt legend leg.txt -R0/1/0/1 -JX1c -Dx${LX}c/${MH}c+w${LW}c+jTL -F+gwhite+p0.4p,gray40 --FONT_ANNOT_PRIMARY=8p,Helvetica
gmt colorbar -Cdepth.cpt -Dx$(awk -v w=$MW 'BEGIN{print w+0.45}')c/$(awk -v h=$MH 'BEGIN{print h-4.7}')c+w$(awk -v w=$MW -v s=$SW 'BEGIN{print s-w-0.45}')c/0.32c+h+jBL \
   -Bxa10f5+l"Hypocentral / centroid depth (km)" -R0/1/0/1 -JX1c
gmt end
echo "done AMF=$AMIN/$AMEAN/$AMAX ZMAX=$ZMAX MW=$MW MH=$MH"
