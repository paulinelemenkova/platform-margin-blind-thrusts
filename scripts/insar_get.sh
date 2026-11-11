W=66.85; S=36.15; E=67.95; N=36.95
P=https://gws-access.jasmin.ac.uk/public/nceo_geohazards/LiCSAR_products
get() {
  t=$1; f=$2; pr=$3; tn=$((10#${t}))
  for k in unw cc; do
    gdal_translate -q -projwin $W $N $E $S /vsicurl/$P.public/$tn/$f/interferograms/$pr/$pr.geo.$k.tif ${f%_*}.$k.tif || \
    gdal_translate -q -projwin $W $N $E $S /vsicurl/$P/$tn/$f/interferograms/$pr/$pr.geo.$k.tif ${f%_*}.$k.tif || true
  done
  for k in E N U hgt; do
    gdal_translate -q -projwin $W $N $E $S /vsicurl/$P/$tn/$f/metadata/$f.geo.$k.tif ${f%_*}.$k.tif || true
  done
}
get 151 151D_05440_131313 20251030_20251111
get 151 151D_05241_131313 20251030_20251111
get 144 144A_05269_131313 20251029_20251110
get 144 144A_05476_141414 20251029_20251110
ls -la *.tif
