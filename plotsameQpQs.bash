#!/bin/bash
# Plot Qp/Qs from independent Qp and Qs using same stn-evt pairs

#workdir=/P/weisq/attentomo/3D1Dtomo/v25h30/10MPa0.27
paradir=~/GoogleDriveMSU/Work/Lau/Qtomo
srcdir=$paradir/program
gmtdir=~/GoogleDriveMSU/Work/GMT/attentomo
nnzQ=4
nnzQ1D=4
top=2
space="v2550h80top2_QsQps"
#invpara="_0.001_100"
# maskdir="/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/"$space"top"$top"/"sameQps_mask
maskdir=$paradir/3D1Dtomo/QpQs/sameQpQs_mask/

for para in 'sameQpQs' 
do
#       workdir="/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/"$space"top"$top"/"$para
#       workdir="$paradir/3D1Dtomo/QpQs/"$para
    workdir=$paradir"/3D1Dtomo/v2550h80top2_QsQps/"$para

cd $workdir
for resultfl in $(ls Qinv_*.ps)
# for resultfl in $(ls  Qinv_1.5e-02_1e+04.ps)
do
      cd $workdir
      invpara=$(echo $resultfl | awk -F"v" '{print $2}' | awk -F".ps" '{print $1}')
      figdir="figure"$invpara
      echo $resultfl $figdir
if ( ! [ -d $figdir ] );then
		mkdir $figdir
# else
# 		rm $figdir/*
#      continue
fi
cp hitsS hitsSP
### atten4gmt Q-model Velocity-model nnzQ
$srcdir/atten4gmt.py<<EOF
$resultfl Q.grid R $nnzQ1D $top
EOF
$srcdir/hits4gmt.py<<EOF
Q.grid R 10 $top
EOF
mv hits*.dep??? $figdir/
mv mask*.dep??? $figdir/
# cp $maskdir/* $figdir/

mv atten*.dep??? $figdir/
mv Q1D.result $figdir/
cd $figdir
$srcdir/excludeQpQs.py

## Plotting
##for attfile in $(ls attenR.dep060)
##do
##      depth=$(echo $attfile | awk '{print substr($1,11,3)}')
##      $gmtdir/attenmap.gmt5 $depth R
##done
# $gmtdir/attenallmapR.gmt5 R
# $gmtdir/atten_xsec.gmt5 sameR
$gmtdir/attenR_xsec.gmt5
done

done
