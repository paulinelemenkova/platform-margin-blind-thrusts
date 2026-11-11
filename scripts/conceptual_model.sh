#!/usr/bin/env bash
set -e
OUT=${OUT:-$(pwd)}; SRC=$(cd "$(dirname "$0")" && pwd)
mkdir -p work_cm && cd work_cm
EQLON=67.3886; EQLAT=36.5461; AZ=340; HL=70
gmt project -C$EQLON/$EQLAT -A$AZ -L-$HL/$HL -G0.25 -Q > prof.txt
for i in 1 2 3; do gmt grdcut @earth_relief_15s -R66.6/68.2/35.8/37.3 -Grel.nc && break; sleep 15; done
gmt grdtrack prof.txt -Grel.nc | awk '{print $3, $4/1000}' > topo.txt
curl -s -m 180 -L -o q.zip https://pubs.usgs.gov/of/2007/1103/downloads/1103.zip && unzip -q -o q.zip
QSHP=$(find . -name 'faults.shp' | head -1)
S=""; [ -f "${QSHP%.shp}.prj" ] || S="-s_srs EPSG:4326"
ogr2ogr -f OGR_GMT amf.gmt "$QSHP" $S -t_srs EPSG:4326 -where "Name LIKE '%Alburz%'" -spat 66.6 35.8 68.2 37.3 -spat_srs EPSG:4326
python3 $SRC/conceptual_model.py . "$OUT"
