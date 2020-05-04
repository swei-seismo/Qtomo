#!/bin/bash

#workdir=/P/weisq/attentomo/3D1Dtomo/v25h30/10MPa0.27
srcdir=~/GoogleDriveMSU/Work/Lau/Qtomo/program
gmtdir=~/GoogleDriveMSU/Work/GMT/attentomo
nnzQ=6
nnzQ1D=4
#invpara="_0.001_100"
maskdir=~/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQps_mask

#for para in 'QpQs' 'QpQs_checker' 'QpQs_synmod' 'QpQs_checkerclean'
for para in 'QpQsjoint'
do
#      workdir="/P/weisq/attentomo/3D1Dtomo/v50h60/"$para
    workdir=~/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/"$para"


cd $workdir
for resultfl in $(ls Qinv_*.sp)
#for resultfl in Qinv_0.01_100000.p Qinv_0.01_10000.p Qinv_0.01_1e+06.p
#for resultfl in Qinv_0.001_0.1.p Qinv_0.001_1.p Qinv_0.0005_0.1.p Qinv_0.0005_1.p Qinv_0.0005_10.p 
#for resultfl in Qinv_1e-03_1e+05.s
do
	cd $workdir
	invpara=$(echo $resultfl | awk -F"v" '{print $2}' | awk -F".sp" '{print $1}')
	figdir="figure"$invpara
	echo $resultfl $figdir
if ( ! [ -d $figdir ] );then
      mkdir $figdir
else
#      rm $figdir/*
      continue
fi

### atten4gmt Q-model Velocity-model nnzQ
#$srcdir/atten4gmt $resultfl Q.grid R $nnzQ1D
#$srcdir/hits4gmt R Q.grid 10
python $srcdir/atten4gmt.py<<EOF
$resultfl Q.grid R $nnzQ1D 2
EOF
#
cp $maskdir/* $figdir

#
mv atten*.dep??? $figdir
mv Q1D.result $figdir
cd $figdir

# Plotting
#for attfile in $(ls attenR.dep100)
#do
#      depth=$(echo $attfile | awk '{print substr($1,11,3)}')
###      $gmtdir/hitsmap.gmt $depth P
#      $gmtdir/attenmap.gmt5 $depth R
#done
#
#$gmtdir/attenallmap.gmt5 S
#$gmtdir/atten_xsec.gmt5 R
##$gmtdir/atten_xsec700.gmt P
##$gmtdir/hits_xsec.gmt P
$gmtdir/attenR_xsec.gmt5

done

done
