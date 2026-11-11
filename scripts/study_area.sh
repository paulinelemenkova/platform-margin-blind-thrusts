#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}
GDRIVE_ID=${GDRIVE_ID:?}
mkdir -p work_sa && cd work_sa

W=64.6; E=69.6; S=35.6; N=37.7
R=$W/$S/$E/$N+r; J=L67.1/36.65/36/37.3/17c
BW=65.6; BE=69.4; BS=35.85; BN=37.45

curl -s -f -m 600 -L -o src.zip "https://drive.usercontent.google.com/download?id=$GDRIVE_ID&export=download&confirm=t"
unzip -q -o src.zip -x '__MACOSX/*'
cpshp() { P=$(find . -name "$1.shp" | head -1); for x in shp shx dbf prj; do [ -f "${P%.shp}.$x" ] && cp "${P%.shp}.$x" "$2.$x"; done; return 0; }
cpshp gregafg greg; cpshp evapfacafg evap; cpshp seismafg seis; cpshp welllogafg well; cpshp fldafg fld
cpshp geoBoundaries-AFG-ADM1_simplified adm1
g() { s=""; [ -f "${2%.shp}.prj" ] || s="-s_srs EPSG:4326"; ogr2ogr -f OGR_GMT "$1" "$2" $s -t_srs EPSG:4326 "${@:3}"; }
CLIP="-clipdst 63.5 35 70.5 38.3"
g gHi.gmt greg.shp $CLIP -where "TYPE='Hi'"; g gLo.gmt greg.shp $CLIP -where "TYPE='Lo'"
g gPz.gmt greg.shp $CLIP -where "TYPE IN ('Pz','Flt')"
g gnames.gmt greg.shp -clipdst $W $S $E $N -where "NAME IS NOT NULL"
g ev10.gmt evap.shp $CLIP -where "FACIES=10"; g ev4.gmt evap.shp $CLIP -where "FACIES=4"
g ev0.gmt evap.shp $CLIP -where "FACIES=0"
g sS.gmt seis.shp -where "SYMBOL='solid'"; g sD.gmt seis.shp -where "SYMBOL<>'solid'"
g well.gmt well.shp; g fld.gmt fld.shp; g adm1.gmt adm1.shp
fc() { ogrinfo -so -al "$@" | awk -F': ' '/Feature Count/{print $2}'; }
NSEIS=$(fc seis.shp); NWELL=$(fc well.shp); NFLD=$(fc fld.shp)

python3 - > glab.txt <<'PY'
import re
OVR={'Murghab Depression':None,'Hindu Kush / Badakhshan':None,
     'N. Karabil-Dawletabad Depressi':(65.22,36.43,'Karabil-Dawletabad Depression'),
     'Obigarm Mega Anticline':(68.85,37.19,None),'Karabekaul Downwarp':(65.38,37.665,None)}
name=None; rings=[]; cur=[]; out={}
def flush():
    global cur
    if cur: rings.append(cur); cur=[]
def area_c(p):
    a=cx=cy=0
    for (x1,y1),(x2,y2) in zip(p,p[1:]+p[:1]):
        k=x1*y2-x2*y1; a+=k; cx+=(x1+x2)*k; cy+=(y1+y2)*k
    return (abs(a/2),cx/(3*a),cy/(3*a)) if a else (0,0,0)
for l in open('gnames.gmt'):
    if l.startswith('# @D'):
        flush()
        if name and rings: out.setdefault(name,[]).extend(rings)
        rings=[]; name=l[4:].strip().split('|')[1].strip('"')
    elif l.startswith('>'): flush()
    elif not l.startswith('#') and l.strip():
        cur.append(tuple(map(float,l.split()[:2])))
flush()
if name and rings: out.setdefault(name,[]).extend(rings)
for n,rs in out.items():
    a,x,y=max(area_c(r) for r in rs)
    if a<=0.02: continue
    o=OVR.get(n,(x,y,n))
    if o is None: continue
    x,y,n=o[0],o[1],(o[2] or n)
    print(f"{x:.3f} {y:.3f} {n}")
PY
cat glab.txt

gmt grdcut @earth_relief_03s -R64.5/69.7/35.5/37.8 -Grelief.nc
gmt grdgradient relief.nc -A315 -Nt0.6 -Gint.nc
gmt makecpt -Csealand -G0/3000 -T0/4500 > topo.cpt

cat > cities.txt <<'EOF'
67.110 36.709 TC C01
66.897 36.758 RM C03
67.686 36.698 LB C02
68.015 36.265 LM C04
68.860 36.728 LM C06
65.752 36.667 LM C05
67.278 37.224 LM C07
68.708 35.944 LM C08
65.938 36.216 LM Sar-e Pol
65.122 36.946 LM Andkhoy
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

gmt begin $OUT/study_area pdf,png E300
gmt set FONT_ANNOT_PRIMARY 8p,Helvetica FONT_LABEL 9p,Helvetica FONT_TITLE 8p,Helvetica \
        MAP_FRAME_TYPE fancy MAP_FRAME_WIDTH 3p MAP_GRID_PEN_PRIMARY 0.25p,gray50 \
        FORMAT_GEO_MAP ddd:mmF MAP_ANNOT_OBLIQUE lat_parallel MAP_TICK_LENGTH_PRIMARY 4p PS_PAGE_COLOR white
gmt grdimage relief.nc -Iint.nc -Ctopo.cpt -R$R -J$J -Y7.2c
gmt plot ev10.gmt -Gp8+r120+f150/80/190+b- -W0.8p,120/40/160 -t35
gmt plot ev4.gmt  -Gp9+r120+f210/110/210+b- -W0.8p,200/90/200,4_2:0p -t35
gmt coast -Di -W0.3p,gray25 -N1/1.0p,gray10,- -I1/0.7p,30/110/200 -I2/0.45p,30/110/200 -C90/160/230 -A20 \
          -Bxa1f10mg30m -Bya30mf10mg30m -BWSEn
gmt plot adm1.gmt -W0.4p,gray25,4_2:0p
gmt plot gPz.gmt -W1.4p,90/90/90
gmt plot gLo.gmt -W1.4p,0/90/170,6_3:0p
gmt plot gHi.gmt -W1.4p,150/60/20
gmt plot fld.gmt -G0/150/130@20 -W0.45p,0/90/80
gmt plot sD.gmt -W0.6p,40/40/40,2_2:0p
gmt plot sS.gmt -W0.6p,40/40/40
gmt plot box.txt -W2.2p,white; gmt plot box.txt -W1.1p,200/0/0,6_3:0p
awk '{print $1,$2}' well.gmt >/dev/null 2>&1 || true
grep -v '^#' well.gmt | grep -v '^>' | gmt plot -Sc0.2c -Gblack -W0.8p,white
echo "67.3886 36.5461" | gmt plot -Sa0.5c -Gred -W0.6p,black
gmt text glab.txt -F+f8.5p,Helvetica-Oblique,60/30/10$WO+jCM
awk '{print $1,$2}' cities.txt | gmt plot -Ss0.22c -Gwhite -W0.8p,black
while read lo la j nm; do
  case $j in RM) dx=-0.16c;; TC) dx=0;; *) dx=0.16c;; esac; case $j in LB) dy=0.12c;; TC) dy=-0.2c;; *) dy=0;; esac
  echo "$lo $la $nm" | gmt text -F+f8p,Helvetica,black$WO+j$j -D$dx/$dy
done < cities.txt
echo "$BE $BN Study area" | gmt text -F+f8p,Helvetica,200/0/0$WO+jBR -D0/0.08c
echo "66.2 37.33 Amu Darya" | gmt text -F+f8p,Helvetica-Oblique,30/90/180$WO+a-12+jCM
gmt basemap -LjBR+c36.6+w50k+f+l"km"+o0.4c/0.45c -F+gwhite@15+p0.4p,gray40+r0.1c
gmt basemap -TdjBR+w1.0c+l,,,N+o0.55c/1.85c
gmt colorbar -Ctopo.cpt -DJBC+w${MW}c/0.3c+o0/0.9c+h+ef -Bxa500f100+l"Elevation (m)"
cat > leg.txt <<EOF
N 3
S 0.3c - 0.55c - 1.4p,150/60/20 0.75c Structural high (uplift)
S 0.3c - 0.55c - 1.4p,0/90/170,6_3:0p 0.75c Structural low (depression)
S 0.3c - 0.55c - 1.4p,90/90/90 0.75c Palaeozoic massif / arch
S 0.3c r 0.5c/0.32c p8+r120+f150/80/190+b- 0.8p,120/40/160 0.75c Salt-anhydrite-limestone (J3)
S 0.3c r 0.5c/0.32c p9+r120+f210/110/210+b- 0.8p,200/90/200 0.75c Anhydrite-dominated (J3)
S 0.3c - 0.55c - 0.6p,40/40/40 0.75c Seismic line ($NSEIS)
S 0.3c c 0.2c black 0.8p,gray40 0.75c Well with logs ($NWELL)
S 0.3c r 0.5c/0.32c 0/150/130 0.45p,0/90/80 0.75c Oil and gas field ($NFLD)
S 0.3c a 0.35c red 0.6p,black 0.75c 2025 mainshock
S 0.3c s 0.22c white 0.8p,black 0.75c City
S 0.3c - 0.55c - 1.0p,gray10,- 0.75c International border
S 0.3c - 0.55c - 1.1p,200/0/0,6_3:0p 0.75c Study area
EOF
gmt legend leg.txt -DJBC+w${MW}c+o0/2.6c -F+gwhite+p0.4p,gray40 --FONT_ANNOT_PRIMARY=8.5p,Helvetica
gmt text -R0/$MW/0/2 -JX${MW}c/2c -M -N -F+f8p,Helvetica,gray20+jTL -Y-6.6c <<EOF
> 0 1.45 0.35c ${MW}c j
Depicts: tectonic regions of the Amu Darya basin margin (structural highs, lows and Palaeozoic massifs), Upper Jurassic (Kimmeridgian-Tithonian) evaporite facies, $NSEIS seismic-reflection line segments, $NWELL wells with scanned logs and $NFLD oil and gas fields on SRTM 3-arc-second relief coloured with the sealand palette; Lambert conformal conic projection (central meridian 67.1@.E, standard parallels 36@.N and 37.3@.N). Data: USGS petroleum GIS (2006); geoBoundaries; SRTM 3 arc-sec; GSHHG; USGS ComCat. Software: GMT $GMTV, GDAL $GDALV, GhostScript $GSV, Python 3. Source: authors.
EOF
gmt end
echo done
