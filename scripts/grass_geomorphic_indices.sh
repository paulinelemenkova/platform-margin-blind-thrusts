#!/usr/bin/env bash
set -e
DEM=${1:-dem_sa.tif}; OUTD=${2:-$(pwd)}
cat > /tmp/geomorph_run.sh <<EOF
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
r.out.gdal -c input=hi output=$OUTD/hi.tif type=Float32 createopt=COMPRESS=DEFLATE --q
r.out.xyz input=ksn separator=space | gzip > $OUTD/ksn.xyz.gz
r.univar -e hi percentile=5,50,95; r.univar -e ksn percentile=5,50,95
EOF
grass -c EPSG:4326 /tmp/geomorph_loc --exec bash /tmp/geomorph_run.sh
