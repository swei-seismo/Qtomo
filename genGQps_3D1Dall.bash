#!/bin/bash
## Build G matrix for Qp/Qs from t*(S) and Qp

tstardir=~/GoogleDriveMSU/Work/Lau/Qtomo/alltstar
space="v2550h80top2_QsQps"
phase="PS"
input=~/GoogleDriveMSU/Work/Lau/Qtomo/input

for top in 2
do
for tdir in $tstardir/fcp0.5-20MPa0.27_s_tstar.dat
#for tdir in $tstardir/fcp0.5-20MPa0.27_s_tstar.dat $tstardir/fcs0.5-20MPa0.27_s_tstar.dat $tstardir/fcp0.5-20MPa0.6_s_tstar.dat $tstardir/fcs0.5-20MPa0.6_s_tstar.dat
do
	stress=$(basename $tdir | awk -F"_" '{print $1}')
      spacedir=~/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/"$space"
      if (! [ -d $spacedir ])then
           mkdir $spacedir
      fi
      workdir=$spacedir"/sameQps_top"$top
#   Qpinfl=$spacedir"/indQpQs/Qinv.p"

if (! [ -d $workdir ]);then
      mkdir $workdir
#else
#      rm $workdir/*
fi

#cp $input/"v25h30lau.vel" $workdir/lau.vel
#cp $input/$space"Q.grid" $workdir/Q.grid
#cp $input/Q.prem $workdir/
#cp $input/Q1D100.grid $workdir/Q1D.grid
#cp $input/station.lst $workdir/
##cp $sdat $workdir/s_tstar.dat
#cp $tdir $workdir/s_tstar.dat
#cp $Qpinfl $workdir/Qinv.p

cd $workdir

~/GoogleDriveMSU/Work/Lau/Qtomo/program/src/Gatten3D1D $phase $top
done
done
