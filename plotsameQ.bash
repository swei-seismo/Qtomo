#!/bin/bash
# Plot independent Qp, Qs and Qp/Qs using same stn-evt pairs

#workdir=/P/weisq/attentomo/3D1Dtomo/v25h30/10MPa0.27
paradir=~/GoogleDriveMSU/Work/Lau/Qtomo
srcdir=$paradir/program
gmtdir=~/GoogleDriveMSU/Work/GMT/attentomo
#nnzQ=4
nnzQ1D=4
top=2
space="v2550h80top2_QsQps"
#invpara="_0.001_100"
# maskdir=/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/v125h155top2/sameQps_mask
maskdir=$paradir/3D1Dtomo/v2550h80top2_QsQps/sameQpQs_mask/

for para in 'sameQpQs'
do
#       workdir="/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/"$space"top"$top"/"$para
#       workdir="/Users/sowei/Work/Lau/Qtomo/3D1Dtomo/QpQs/"$para
     workdir=$paradir"/3D1Dtomo/v2550h80top2_QsQps/"$para
#    workdir=$paradir/3D1Dtomo/v2550h80top2_QsQps/"$para"_top"$top"
#		workdir=$paradir/3D1Dtomo_nosite/v125h155top2/sameQpQs_checker/

      cd $workdir
#    for resultfl in $(ls Qinv_*.p)
#      for resultfl in $(ls Qinv_[23]*e-02_1e+*.p)
      for resultfl in Qinv_8.0e-03_1e+04.p
#      for resultfl in Qinv_1.5e-02_1e+04.p
      do
	      cd $workdir
      	invpara=$(echo $resultfl | awk -F"v" '{print $2}' | awk -F".p" '{print $1}')
            resultfl2="Qinv"$invpara".s"
            resultfl3="Qinv"$invpara".ps"
            if ( ! ( [ -f $resultfl2 ] && [ -f $resultfl2 ] ));then
                  echo $invpara "not complete"
                  continue
            fi
	      figdir="figure"$invpara
      	echo $resultfl $figdir
            if ( ! [ -d $figdir ] );then
                  mkdir $figdir
             else
                   rm $figdir/*
 #                   continue
            fi
       cp $maskdir/* $figdir/

            ########################### PLOT Qp ##############################
### atten4gmt Q-model Velocity-model nnzQ
/usr/bin/python $srcdir/atten4gmt.py<<EOF
$resultfl Q.grid P $nnzQ1D $top
EOF
/usr/bin/python $srcdir/hits4gmt.py<<EOF
Q.grid P 10 $top
EOF
      mv hitsP.dep??? $figdir/
      mv maskP.dep??? $figdir/
#      cp $maskdir/* $figdir

		mv atten*.dep??? $figdir/
		mv Q1D.result $figdir/
		cd $figdir
#		rm *000*

		## Plotting
# 		for attfile in $(ls attenP.dep050)
# 		do
# 		      depth=$(echo $attfile | awk '{print substr($1,11,3)}')
# 		##      $gmtdir/hitsmap.gmt $depth P
# 		      $gmtdir/attenmap.gmt5 $depth P
# 		done
		$gmtdir/atten_xsec.gmt5 P

		# Plot model standard deviation
		modresfl=resultfl"_varm"
		if( ! [ -e $modresfl ] );then
			 $gmtdir/attenstd_xsec.gmt5 Pstd
		fi

            ################ PLOTING Qs ###########################
            cd $workdir
            resultfl2="Qinv"$invpara".s"
            echo $resultfl2 $figdir
### atten4gmt Q-model Velocity-model nnzQ
/usr/bin/python $srcdir/atten4gmt.py<<EOF
$resultfl2 Q.grid S $nnzQ1D $top
EOF
/usr/bin/python $srcdir/hits4gmt.py<<EOF
Q.grid S 10 $top
EOF
      mv hitsS.dep??? $figdir/
      mv maskS.dep??? $figdir/
            mv atten*.dep??? $figdir/
		mv Q1D.result $figdir/
		cd $figdir
#		rm *000*
# 		for attfile in $(ls attenS.dep050)
# 		do
# 		      depth=$(echo $attfile | awk '{print substr($1,11,3)}')
# 		      $gmtdir/attenmap.gmt5 $depth S
# 		done
		$gmtdir/atten_xsec.gmt5 S

		# Plot model standard deviation
		modresfl=resultfl"_varm"
		if( ! [ -e $modresfl ] );then
			 $gmtdir/attenstd_xsec.gmt5 SSstd
		fi
##################### PLOT Qp/Qs ######################
#             cd $workdir
#             resultfl3="Qinv"$invpara".ps"
#             cp hitsS hitsSP
# ## atten4gmt Q-model Velocity-model nnzQ
# $srcdir/atten4gmt.py<<EOF
# $resultfl3 Q.grid R $nnzQ1D $top
# EOF
# $srcdir/hits4gmt.py<<EOF
# Q.grid R 10 $top
# EOF
# 		mv atten*.dep??? $figdir/
# 		mv hits*.dep??? $figdir/
# 		mv mask*.dep??? $figdir/
# 		mv Q1D.result $figdir/
# 		cd $figdir
# 		$srcdir/excludeQpQs.py
# # 		$gmtdir/attenallmapR.gmt5 R
# 		$gmtdir/atten_xsec.gmt5 sameR
      done

done
