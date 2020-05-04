#!/bin/bash
# Checker board test for independent Qp and Qs using same stn-evt pairs
# Generate synthetic t*(P) and t*(S), add noise, and then G matrix

space="v125h155"
top=2
input=~/GoogleDriveMSU/Work/Lau/Qtomo/input
spacedir="/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/"$space"top"$top
progdir=~/GoogleDriveMSU/Work/Lau/Qtomo/program
realdatp=$spacedir"/sameQpQs/p_tstar.dat"
realdats=$spacedir"/sameQpQs/s_tstar.dat"

#for ctype in "checker"
#do
#	inmodp=$input"/"$ctype"p_"$space".in"
#   inmods=$input"/"$ctype"s_"$space".in"
#	workdir=$spacedir"/sameQpQs_"$ctype
#
#if (! [ -d $workdir ]);then
#      mkdir $workdir
#else
#      rm -r $workdir/*
#fi
#
#cp $input/"v25h30lau.vel" $workdir/lau.vel
#cp $input/station.lst $workdir/
#cp $input/$space"Q.grid" $workdir/Q.grid
#cp $input/Q.prem $workdir/
#cp $input/Q1D100.grid $workdir/Q1D.grid
#cp $realdatp $workdir/
#cp $realdats $workdir/
#cd $workdir

cd /Users/swei/GoogleDriveMSU/Work/Lau/Qtomo/3D1Dtomo/v2550h80top2_QsQps/sameQpQs_checker2

#cp $inmodp $workdir/checkerp.in
cp checkerp.in checker.in
$progdir/src/syn_tstar P $top
/usr/bin/python $progdir/addnoise.py<<EOF
p_tstar.cleansyn p_tstar.syn P
EOF
$progdir/src/Gatten3D1D PP $top

#cp $inmods $workdir/checkers.in
cp checkers.in checker.in
$progdir/src/syn_tstar S $top
/usr/bin/python $progdir/addnoise.py<<EOF
s_tstar.cleansyn s_tstar.syn S
EOF
$progdir/src/Gatten3D1D SS $top

#done
