#!/bin/bash

space="v2550h80top2_QsQps"
ntop=2
input=~/GoogleDriveMSU/Work/Lau/Qtomo/input
#workdir="/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/"$space"top"$ntop"/sameQpQs"
workdir=~/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/"$space"/sameQpQs
echo $workdir

#python ~/GoogleDriveMSU/Work/Lau/Qtomo/program/sametstar.py

cp $input/"v25h30lau.vel" $workdir/lau.vel
cp $input/$space"Q.grid" $workdir/Q.grid
cp $input/Q.prem $workdir/
cp $input/Q1D100.grid $workdir/Q1D.grid
cp $input/station.lst $workdir/

cd $workdir

~/GoogleDriveMSU/Work/Lau/Qtomo/program/src/Gatten3D1D P $ntop
~/GoogleDriveMSU/Work/Lau/Qtomo/program/src/Gatten3D1D S $ntop
