#!/usr/bin/env bash
set -e
DEM=${1:-dem_sa.tif}; OUTD=${2:-$(pwd)}
cat > /tmp/sf_run.sh <<EOF
set -e
r.in.gdal -o input=$DEM output=dem --q
g.region raster=dem
r.watershed -a elevation=dem threshold=4000 accumulation=acc drainage=dir basin=basin memory=1600 --q
r.stream.extract elevation=dem accumulation=acc threshold=300 stream_raster=str memory=1600 --q
for m in average min max; do r.stats.zonal base=basin cover=dem method=\$m output=b\$m --q; done
r.mapcalc "hi = if(bmax - bmin > 50, (baverage - bmin) / (bmax - bmin), null())" --q
r.slope.aspect elevation=dem slope=slp format=percent --q
r.neighbors input=slp output=slps size=5 method=average --q
r.mapcalc "ksn = if(isnull(str), null(), (slps/100.0) * exp(0.45*log(abs(acc) * (ewres()*111320*cos(y())) * (nsres()*110570))))" --q
g.region n=37.45 s=35.85 w=65.6 e=69.4 res=0:00:36 -a
r.resamp.stats -w input=hi output=hi01 method=average --q
r.resamp.stats -w input=ksn output=ksn01 method=median --q
r.out.gdal -c -f input=hi01 output=$OUTD/hi01.tif type=Float32 --q
r.out.gdal -c -f input=ksn01 output=$OUTD/ksn01.tif type=Float32 --q
r.univar hi01 | grep -E "^n:|mean"; r.univar ksn01 | grep -E "^n:|mean"
EOF
grass -c EPSG:4326 /tmp/sf_loc --exec bash /tmp/sf_run.sh
