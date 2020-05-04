ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Main program to invert 3D-1D attenuation tomography.               c
c  Make G-matrix for attenuation using pseudo-bending method          c
c    [Um & Thurber, 1987]. Bending starts with straight lines from    c
c    source to receiver.                                              c
c                                                                     c
c  by S. Wei, June, 2013, based on S. Pozgay's 2-D version            c
c  Modifying to allow different grids for V and Q, October 2014       c
c  Modifying to allow invert 1D Q below V model, December 2014        c
c  Modifying to skip top layer of 3D Q  model, June 2014              c
c  Input: P/S/PP/SS/PSR (P/S, Syn P/S, Qp/Qs based on Qp)             c
c         ntop (0: no top layer skipped, but add 1 layer at 12 km;    c
c               1: first layer at 12 km skipped;                      c
c               2: first layer at 0 km skipped)                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PROGRAM MAIN

      parameter (nstmax = 100)     ! max number of stations
      parameter (nxmax = 45, nymax = 45, nzmax = 30)  ! max nodes
      parameter (maxpts = 2048)       ! max # internal points in raypath
      parameter (maxparam = nxmax*nymax*nzmax+20)      ! max velocity grid
      parameter (maxdata = 21000+maxparam)      ! max measuments + smoothing
      parameter (maxzq = 700) ! max depth of Q model

      real lat0, lon0, beta   ! refrence node
      real lat(nxmax,nymax), lon(nxmax,nymax), temp1, temp2
      real latq(nxmax,nymax), lonq(nxmax,nymax)
      real x(nxmax), y(nymax), z(nzmax), xmax, ymax         ! node locations
      real xq(nxmax), zq(nzmax), yq(nymax), zq1D(nzmax), attin1D(nzmax)
      real stlon(nstmax), stlat(nstmax), scrust(nstmax)     ! station locations
      real xsta(nstmax), ysta(nstmax), zsta(nstmax)
      real ylate, xlone, ze, elat, elon, edep, tstar, tster, grade
      real xe0, ye0, ze0, xs0, ys0, zs0, zmoho, tt, ttP, LttS
      real aveQinv,ttin,maxtt, QpQs, xx, yy, ptstar,pteer
      real*8 vps(2,nzmax,nymax,nxmax), sn(2,nzmax,nymax,nxmax)    ! velocity, slowness
      real*8 attin(nzmax,nymax,nxmax), Qout(2,nzmax), Qtop, toptt     ! 1/Qp input for Qp/Qs
      real grow(maxparam), sumgrow     ! G-matrix row for given raypath
      real gbig(maxdata,maxparam), gbigbak(maxdata,maxparam)
      real datav(maxdata), covdata(maxdata), solnv(maxparam)
      real*8 xarc(0:maxpts+2), yarc(0:maxpts+2), zarc(0:maxpts+2)
C     real*8 xarck(0:maxpts+2), yarck(0:maxpts+2), zarck(0:maxpts+2)
      real wd(maxdata)
      real smdata(maxdata), smwd(maxdata), smgbig(maxdata,maxparam)
      real smlamda, dplamda

      integer nnx, nny, nnz, nnxQ, nnyQ, nnzQ, nnzQ1D, ntop
      integer hitp(maxparam), p, pt, ifskip, lps

      character*3 pors,ntops
      character*4 stname(nstmax), stn, lamda
      character*10 phi(nstmax)

      common /a/ x, y, z, xq, yq, zq
      common /b/ ipolar,nnx,nny,nnz,re,pi,r2d,nnxQ,nnyQ,nnzQ,nnzQ1D
      common /c/ sn
      common /d/ attin, zq1D, attin1D

c     I/O files
      narguments = iargc()
      ntop=1
      if (narguments.lt.1) then
        write(*,*) 'Usage: Gatten3D P/S/PS/SP/PP/SS'
        write(*,*) '(PS: Qp/Qs; SP: S-P diff)'
        stop
      endif
      call getarg(1,pors)
      if (narguments.eq.2) then
        call getarg(2,ntops)
        read(ntops,'(i1)') ntop
      endif
      if (pors.eq.'P' .or. pors.eq.'p') then
        write(*,*) 'Invert Qp with t*(P)'
        lps=1
        icase=1
        open (8, file='p_tstar.dat')
        open (21+lps,file='raypaths.p')     ! output P ray paths
      else if (pors.eq.'S' .or. pors.eq.'s') then
        write(*,*) 'Invert Qs with t*(S)'
        lps=2
        icase=2
        open (8, file='s_tstar.dat')
        open (21+lps,file='raypaths.s')     ! output S ray paths
        open (31+lps,file='raypaths.k')     ! output S ray paths with Qk
      else if (pors.eq.'PP') then
        write(*,*) 'Invert Qp with synthetic t*(P)'
        lps=1
        icase=1
        open (8, file='p_tstar.syn')
        open (21+lps,file='raypaths.p')     ! output P ray paths
      else if (pors.eq.'SS') then
        write(*,*) 'Invert Qs with synthetic t*(S)'
        lps=2
        icase=2
        open (8, file='s_tstar.syn')
        open (21+lps,file='raypaths.s')     ! output S ray paths
        open (31+lps,file='raypaths.k')     ! output S ray paths with Qk
      else if (pors.eq.'PS' .or. pors.eq.'ps' .or.
     &         pors.eq.'Ps' .or. pors.eq.'pS' ) then
        write(*,*) 'Invert Qp/Qs with t*(S) and known Qp'
        lps=2
        icase=3
        open (8, file='s_tstar.dat')
        open (21+lps,file='raypaths.s')     ! output S ray paths
      else if (pors.eq.'PSR' .or. pors.eq.'psr' .or.
     &         pors.eq.'Psr' .or. pors.eq.'pSR' ) then
        write(*,*) 'Invert Qp/Qs with synthetic t*(S) and known Qp'
        lps=2
        icase=3
        open (8, file='s_tstar.syn')
        open (21+lps,file='raypaths.s')     ! output S ray paths
C     else if (pors.eq.'SP' .or. pors.eq.'sp' .or.
C    &         pors.eq.'Sp' .or. pors.eq.'sP' ) then
C       write(*,*) 'Invert Qp/Qs with t*(S)-t*(P)'
C       lps=1
C       icase=4
C       open (8, file='s-p_finaleqa27.dat')
C       open (21+lps,file='raypaths.sp')    ! output S-P ray paths
C     else if (pors.eq.'K' .or. pors.eq.'k' ) then
C       write(*,*) 'Invert Qk with t*(P) and t*(S)'
C       lps=2
C       icase=4
C       open (8, file='s_tstar.dat')
C       open (21+lps,file='raypaths.k')     ! output S ray paths
      else
        write(*,*) 'Usage: Gatten3D P/S/PS/SP/PP/SS'
        write(*,*) '(PS: Qp/Qs; SP: S-P diff)'
        stop
      endif

      open (9, file='station.lst')  ! station lat, lon, z, crust, name
      open (10,file='lau.vel')     ! 3-D velocity model
      open (11,file='Q.grid')     ! 3-D Q model
      open (12,file='Q1D.grid')   ! 1-D Q model
      open (13,file='Q.prem')       ! 1-D Q values of PREM
      open (30,file='event.lst')    ! output used earthquake location
      open (21,file='raypaths.bad')     ! output bad ray paths
      elon=0
      elat=0
      edep=0 

c     Other parameters
      smooth = 8000.
      alpha = 35

      pi = 4.*atan(1.)
      r2d = 180./pi
      re = 6377.        ! radius of Earth ~ 20 deg lat
      ipolar = 1        ! flat Earth = 0, curved Earth = 1
C     iall = 0
      maxtt = 0

ccccc velocity model input ccccc
      read (10,*) nnx, nny, nnz    ! number of nodes in hz, vt directions
      if (nnx.gt.nxmax .or. nny.gt.nymax .or. nnz.gt.nzmax) then
        write (*,*) "ERROR: need larger dimensions of n?max"
        stop
      endif
c     Reference node and beta (counterclockwise angle in rad)
      read (10,*) lat0, lon0, beta
c     Read node lat and lon, then convert to Cartesian coordinates
      do i = 1, nnx
        do j = 1, nny
          read (10,*) lat(i,j), lon(i,j)
          call lonlat2xy(lon0,lat0,beta,lon(i,j),lat(i,j),x(i),y(j))
        enddo
      enddo
c     Find the max x and y
      xmax = x(nnx)
      ymax = y(nny)
      read (10,*) (z(k), k=1,nnz)   ! positive down
      do p = 1,2
        do i = 1, nnx         ! i-th column
          do j = 1, nny       ! j-th row
            read (10,*) (vps(p,k,j,i), k=1,nnz) ! velocity at node
            do k = 1, nnz     ! k-th layer
              sn(p,k,j,i)=1./vps(p,k,j,i) ! slowness
            enddo
          enddo
        enddo
      enddo
      close(10)
      write(*,*) "Read the velocity model successfully"
      open (27,file='cartesian.Vgrid')
      inode = 0
      write(27,*) xmax, ymax
      do i = 1, nnx
        do j = 1, nny
          inode = inode+1
          write(27,*) x(i), y(j), inode
        enddo
      enddo
      close(27)

ccccc Input Q values of PREM for reference outside of Q model
      read(13,*) (Qout(1,k), Qout(2,k), k=1,nnz)
      close(13)
      Qtop=Qout(lps,1)

ccccc 3D Q model input ccccc
      read (11,*) nnxQ, nnyQ, nnzQ    ! number of nodes in hz, vt directions
c     Reference node and beta (counterclockwise angle in rad)
      read (11,*) latq0, lonq0, betaq
c     Read node lat and lon, then convert to velocity Cartesian coordinates
      do i = 1, nnxQ
        do j = 1, nnyQ
          read (11,*) latq(i,j), lonq(i,j)
          call lonlat2xy(lon0,lat0,betaq,lonq(i,j),latq(i,j),
     &             xq(i),yq(j))
        enddo
      enddo
c     Find the max x and y
      xqmax=xq(nnxQ)
      xqmin=xq(1)
      yqmax=yq(nnyQ)
      yqmin=yq(1)
ccc     If skip the top layer of 10 km
      if ((ntop.eq.1).or.(ntop.eq.2)) then
        read (11,*) (zq(k), k=1,nnzQ)   ! positive down
      else if (ntop.eq.0) then
        nnzQ = nnzQ + 1
        zq(1) = 12.0
        read (11,*) (zq(k), k=2,nnzQ)
      else
        write(*,*) "ntop = ",ntop
        stop
      endif
      close(11)
      open (27,file='cartesian.Qgrid')
      inode = 0
      write(27,*) xqmax, yqmax
      do i = 1, nnxQ
        do j = 1, nnyQ
          inode = inode+1
          write(27,*) xq(i), yq(j), inode
        enddo
      enddo
      close(27)
ccccc 1D Q model input ccccc
      read (12,*) nnzQ1D
      do i = 1, nnzQ1D
        read (12,*) zq1D(i)
      enddo
      close(12)
      nparm = nnxQ*nnyQ*nnzQ + nnzQ1D
      write(*,*) "nnxQ, nnyQ, nnzQ, nparm: ", nnxQ, nnyQ, nnzQ, nparm

ccccc Input station info ccccc
      read (9,*) nsta
      do i = 1, nsta
        read (9,*) stname(i),stlon(i),stlat(i),zsta(i),scrust(i),phi(i)
        if (stlon(i) .lt. 0.) stlon(i) = stlon(i) + 360.
        call lonlat2xy(lon0,lat0,beta,stlon(i),stlat(i),xsta(i),ysta(i))
c        if (xsta(i).lt.0 .or. xsta(i).gt.xmax .or. ysta(i).lt.0 
c     &            .or. ysta(i).gt.ymax) then
c          write (*,*) "ERROR: Station lies outside of the grid", stname
c          stop
c        endif
        if (zsta(i).lt.0.) then
          zsta(i) = -zsta(i)          ! positive down
        else
          zsta(i) = 0
        endif
      enddo
      close(9)
      open(28,file='cartesian.sta')
      do i = 1, nsta
        write(28,120) stname(i),xsta(i),ysta(i),zsta(i)
      enddo
      close(28)

cccccc Read 1/Qp to invert Qp/Qs cccccc
      if (icase.eq.3) then
        open(44,file='Qinv.p')
        do i = 1, nnxQ
          do j = 1, nnyQ
            do k = 1, nnzQ
              read(44,*) attin(k,j,i)
              if (attin(k,j,i).eq.0) then
                attin(k,j,i)=0.0001
              endif
            enddo
          enddo
        enddo
        do i = 1, nnzQ1D
          read (44,*) attin1D(i)
c          if (attin1D(i).eq.0) then
            attin1D(i)=1/370
c          endif
        enddo
        close(44)
      endif

cccccc Run main program cccccc
c      do icase = 1, 1

c     find number of data
        nrays = 0
20      read (8,*,end=99)
        nrays = nrays+1
        goto 20
99      continue
        rewind 8
C       if (icase .eq. 1) then
C         nraysp = nrays 
C       else if (icase. eq. 2) then
C         nrayss = nrays 
C       endif

        do i = 1, nrays
          datav(i) = 0.       ! Data-matrix
          covdata(i) = 0.
          do j = 1, nparm
            gbig(i,j) = 0.    ! G-matrix
c            gbigbak(i,j) = 0.
          enddo 
        enddo

cccccc Loop over each raypath cccccccccccccc
        iknt = 0
        iskip = 0
        nbigmis10 = 0
        nbigmis5 = 0
        do iray = 1, nrays
          if (lps.eq.1) then
            read (8,*) stn,ylate,xlone,ze,tstar,tster,xx,aveQinv
          elseif (lps.eq.2) then
            read (8,*) stn,ylate,xlone,ze,tstar,tster,ptstar,pteer,
     &              aveQinv,QpQs
          endif
c         find the event location
          if (xlone .lt. 0) xlone = xlone + 360
          call lonlat2xy(lon0,lat0,beta,xlone,ylate,xe0,ye0)
          if (xe0.lt.0 .or. xe0.gt.xmax .or. ye0.lt.0 
     &            .or. ye0.gt.ymax) then
C            write (*,*) "ERROR: Earthquake lies outside of the grid",
C    &                        ylate, xlone
            iskip = iskip + 1
            go to 513
          endif
          if (ze .le. 35.) then     ! Skip event shallower than 35 km
            iskip = iskip + 1
C            write (*,*) "ERROR: Earthquake is too shallow", ze
            go to 513
          endif
          if (tster .eq. 0) then
            iskip = iskip + 1
C            write (*,*) "ERROR: t* has wrong error", tster
            go to 513
          endif
          if ((lps .eq. 2) .and. (icase .eq. 2)) then
            if ((ptstar .eq. 0) .or. (QpQs .eq. 1.75)) then
              iskip = iskip + 1
C             write (*,*) "ERROR: t*(P) is zero", tster
              go to 513
            else
              ttS0=aveQinv*tstar
              ttP0=aveQinv*ptstar/QpQs
              Lratio=4.*ttP0*ttP0/3./ttS0/ttS0
              if (Lratio .eq. 1) then
                aveQkinv=0.
              else
                aveQkinv=(1/QpQs-Lratio)*aveQinv/(1.-Lratio)
              endif
            endif
          endif
          ze0 = ze
c         find the station
          lsta = 0
          do i = 1, nsta
            if (stn .eq. stname(i)) lsta = i
          enddo
          if (lsta .eq. 0) then
            write (*,*) "Station not found ", stn
            stop
          endif
          if (xsta(lsta).lt.xqmin .or. xsta(lsta).gt.xqmax .or.
     &       ysta(lsta).lt.yqmin .or. ysta(lsta).gt.yqmax) then
c            write (*,*) "ERROR: Station lies outside of the grid ",stn
            iskip = iskip + 1
            go to 513
          endif
          xs0 = xsta(lsta)
          ys0 = ysta(lsta)
          zs0 = zsta(lsta)
c         find Moho depth
          zmoho = zs0 + scrust(lsta)

ccc -- calculate raypath
C      lps = 1 for P & S-P,          =2 for S,
c          write (*,*) lps,xe0,ye0,ze0,xs0,ys0,zs0,zmoho,maxpts
C       if (icase .eq. 4) then   ! Qk needs LttS = Vs^3/Vp^2/len
C         call bendray(3,xe0,ye0,ze0,xs0,ys0,zs0,zmoho,maxpts,
C    &              xarc,yarc,zarc,npts,ttP,LttS)
C         write(*,*) "Finish ray tracing"
C         write(*,*) ttP, LttS
CC         LttS=0.
C       else  ! Others only need xarc,yarc,zarc etc
          call bendray(lps,xe0,ye0,ze0,xs0,ys0,zs0,zmoho,maxpts,
     &              xarc,yarc,zarc,npts,tt,LttS)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Subroutine that does pseudo-bending. Returns tt and     c
c   raypath (xarc, yarc, zarc, npts).  Bending starts with  c
c   straight lines from source to receiver.                 c
c     npts = # interior pts that can bend.                  c
c     kpt = 0 => hypocenter;                                c
c     kpt = npts+1 => moho piercing pt;                     c
c     kpt = npts+2 => station location.                     c
c     INPUT: lps,xe0,ye0,ze0,xs0,ys0,zs0,zmoho,maxpts,iall  c
c     OUTPUT: xarc,yarc,zarc,npts (raypath),tt (travel-time)c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C         write(*,*) "Finish ray tracing"
          if (tt .gt. maxtt) maxtt=tt
          if (aveQinv.eq.0) aveQinv=aveQinv+0.5
          ttin=1000*tstar/aveQinv
c          if (abs(tt-ttin) .gt. 5) then
c            write(*,*) "Warning: ray tracing may fail for ",
c     $                  stn, ylate, xlone, ze
c            write(*,*) "Ray tracing tt - Antelope tt: ", (tt-ttin)
c            write(*,*) "Raytracing: ", tt, " Antelope: ", ttin
c          endif

cccccc    Travel time in the top layer
          if (ntop.eq.1) then
            call layertt(icase,maxpts,xarc,yarc,zarc,npts,toptt,ipt1)
cccccc   Add one point: ipt => piercing point at 25 km   cccccc
          elseif ((ntop.eq.0).or.(ntop.eq.2)) then
            toptt=0
            ipt1=npts+2
          endif
c          write(*,*) npts,toptt
C       endif

          call getgrowQ(icase,maxpts,xarc,yarc,zarc,npts,ipt1,
     &                    grow,ifskip)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Subroutine that builds G-matrix row for given raypath(s) c
c  From the source (0) to 25 km (ipt)                       c
c    id = 1 for P; id = 2 for S                             c
c  INPUT: icase,maxpts,xarc,yarc,zarc,npts                  c
c  OUTPUT: grow,ifskip                                      c
c   Add one point: ipt => piercing point into the Q model   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C         write(*,*) "Finish adding row of G matrix"
          if (ifskip.eq.1) then
c            write (*,*) "SKIP: Earthquake lies outside of the grid",
c     &                        ylate, xlone
            iskip = iskip + 1
            go to 513
          endif
c     Skip INFs and NaNs
          do i = 1, nparm
c            if (.not.(grow(i).le.0.) .and. .not.(grow(i).ge.0.)) then
            if (grow(i) .ne. grow(i)) then
              grow(i) = 0.
              iskip = iskip + 1
c              write(*,*) "Skip NaN"
              go to 513
            endif
            if (abs(grow(i)) .gt. 1.e8) then
              grow(i) = 0.
c              write(*,*) "Skip INF"
              iskip = iskip + 1
              go to 513
            endif
          enddo

c     Output raypath
          iknt = iknt + 1
          if (lps.eq.1) then
            write(21+lps,*) "1234567 ", iknt, aveQinv
C         elseif (icase.eq.4) then
C           write(21+lps,*) "1234567 ", QpQs, aveQkinv
          elseif (lps.eq.2) then
            write(21+lps,*) "1234567 ", QpQs, aveQinv
          endif
          if (icase .eq. 2) then
            write(31+lps,*) "1234567 ", QpQs, aveQkinv
          endif
c          write(21+lps,*) "111 222 ", iknt,stn, ylate, xlone, ze
          do pt = 0, npts+2  ! write out raypath
            write(21+lps,130) xarc(pt), yarc(pt), zarc(pt)
            if (icase .eq. 2) then
              write(31+lps,130) xarc(pt), yarc(pt), zarc(pt)            
            endif
          enddo

C     Sum up each row of G matrix, should be close to travel time
        if ((icase .eq. 1) .or. (icase .eq. 2)) then          ! Qp or Qs
          if (lps .eq. 1) then
            dtcrust=scrust(lsta)/6.1-zmoho/8.0
          else
            dtcrust=scrust(lsta)/3.5-zmoho/4.19
          endif
          sumgrow=toptt+dtcrust
          do i = 1, nparm
            sumgrow=sumgrow+grow(i)
          enddo
          aveQinv3=1000*tstar/sumgrow
c          if (abs(aveQinv-aveQinv3).gt.5) then
          if (abs(sumgrow-tt).gt.10.) then
            write(*,*) "Error: Constructing G row may fail"
            write(*,*) "Difference of tt: ",abs(sumgrow-tt),tt,sumgrow
c            write(*,*) stn,ylate,xlone,ze
            nbigmis10 = nbigmis10 + 1
c            write(*,*) "Average 1000/Q:", aveQinv3, (aveQinv-aveQinv3)
            if (lps.eq.1) then
              write(21,*) "1234567 ", iknt, aveQinv
            elseif (lps.eq.2) then
              write(21,*) "1234567 ", QpQs, aveQinv
            endif
            do pt = 0, npts+2  ! write out raypath
              write(21,130) xarc(pt), yarc(pt), zarc(pt)
            enddo
c            sumgrow1=0
c            do i = 1, nnxQ*nnyQ*nnzQ
c              sumgrow1=sumgrow1+grow(i)
c            enddo
c            sumgrow2=0
c            do i = nnxQ*nnyQ*nnzQ+1, nparm
c              sumgrow2=sumgrow2+grow(i)
c            enddo
c            write(*,*) sumgrow1, sumgrow2
          endif
          if (abs(sumgrow-tt).gt.5) then
            nbigmis5 = nbigmis5 + 1
          endif
        endif

C     Adding row to G-matrix and data-matrix
          do i = 1, nparm
            gbig(iknt,i) = grow(i)
c            gbigbak(iknt,i) = grow(i)
          enddo
c          write(*,*) tster, tstar
C         if (icase .eq. 4) then    ! Qk using t*(P) and t*(S)
Cc     data = t*(P)/ttP-4/3*t*(S)*LttS
C           datav(iknt) = ptstar/ttP - 4./3.*tstar*LttS
C           covdata(iknt) = tster
C           wd(iknt) = 1./tster       ! data weight by tstar errors
C         else
c     data = t* - toptt/Qtop (top layer)
            datav(iknt) = tstar - toptt/Qtop
            covdata(iknt) = tster
c          datavbak(iknt) = tstar
            wd(iknt) = 1./tster       ! data weight by tstar errors
C         endif

c     Output used earthquake location
          if (elon.ne.xlone .or. elat.ne.ylate .or. edep.ne.ze) then
            write(30,130) xlone, ylate, ze
            elon=xlone
            elat=ylate
            edep=ze
          endif

513       continue
c          write(*,*) "Finishi ray: ",iknt,stn, ylate, xlone, ze
        enddo
ccccc End of loop over each raypath (iknt)
        close(21+lps)
        if (icase .eq. 2) close(31+lps)
        close(21)
        close(30)
        ndata = iknt
        write(*,*) "cccccc FINISH BUILDING G-MATRIX cccccc"
        write(*,*) "nparm, ndata, skip", nparm, ndata, iskip
        write(*,*) "Maximum travel time = ",maxtt
        write(*,*) "nbigmis10, nbigmis5: ", nbigmis10, nbigmis5

C      SUM COLUMNS OF GBIG TO GET HIT COUNT
c      hitc is column vector of length nparm
        if (icase .eq. 1) then              ! for P
          open (77,file='hitsP')
        else if (icase.eq.2 .or. icase.eq.3) then              ! for S
          open (77,file='hitsS')
        endif
        if (icase .eq. 1) then              ! for P
          open (78,file='hitsP1D')
        else if (icase.eq.2 .or. icase.eq.3) then              ! for S
          open (78,file='hitsS1D')
        endif
        do j = 1, nparm
          hitp(j) = 0.
          do i = 1, ndata
            if (gbig(i,j) .gt. 0.) then
              hitp(j) = hitp(j) + 1
            endif
          enddo
        enddo
        do i = 1, nnxQ         ! i-th column
          do j = 1, nnyQ       ! j-th row
cccc            write(77,110) (((j-1)*nnx*nnzQ+(i-1)*nnzQ+k), k=1,nnzQ)
            write(77,110) (hitp((i-1)*nnyQ*nnzQ+(j-1)*nnzQ+k), k=1,nnzQ)
          enddo
        enddo
        write(78,110) (hitp(k), k=nnxQ*nnyQ*nnzQ+1, nparm)
        close(77)
        close(78)
        write(*,*) "Finish G matrix"

c     Adding smoothing matrix for 3D Q model
        call smoothg(smooth,maxparam,maxdata,smdata,smwd,smgbig,nkp)
        write(*,*) "Finish smoothing, nkp=",nkp
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   SUBROUTINE TO WRITE CONSTRAINT ROWS TO G-MATRIX  ccc
ccc   2nd ORDER SMOOTHING CONTRAINTS APPEDNING TO G    ccc
ccc    INPUT:  smo,nmod,nd                             ccc
ccc    OUTPUT: DV,WD,GB,nk                             ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        do ikp = 1, nkp
c          do inode = 1, nparm
c            gbig(ikp+iknt,inode) = smgbig(ikp,inode)
c          enddo
c          datav(ikp+iknt) = smdata(ikp)
c          wd(ikp+iknt) = smwd(ikp)
c        enddo
c        iknt = iknt + nkp
c        
cc     Applying the weighting of data
c        do i = 1,iknt          ! G = wd*G; d = wd*d (apply data weights)
c          do j = 1, nparm
c            gbig(i,j) = gbig(i,j)*wd(i)
c          enddo
c          datav(i) = datav(i)*wd(i)
c        enddo

c     Output G-matrix, data matrix, weighting matrix
c     and smoothing matrix for G-matrix, data matrix, weighting matrix
        if (icase .eq. 1) then              ! for P
          open (13,file='GmatP')
          open (14,file='DmatP')
          open (15,file='WDmatP')
          open (16,file='GmatPsm')
          open (17,file='DmatPsm')
          open (18,file='WDmatPsm')
          open (19,file='ErrDmatP')
        else if (icase .eq. 2) then              ! for S
          open (13,file='GmatS')
          open (14,file='DmatS')
          open (15,file='WDmatS')
          open (16,file='GmatSsm')
          open (17,file='DmatSsm')
          open (18,file='WDmatSsm')
          open (19,file='ErrDmatS')
        else if (icase .eq. 3) then              ! for Qp/Qs
          open (13,file='GmatPS')
          open (14,file='DmatPS')
          open (15,file='WDmatPS')
          open (16,file='GmatPSsm')
          open (17,file='DmatPSsm')
          open (18,file='WDmatPSsm')
          open (19,file='ErrDmatPS')
        endif
        do i = 1, iknt
          write(13,100) (gbig(i,j), j=1,nparm)   ! Gmat?
          write(14,*) datav(i)                  ! Dmat?
          write(19,*) covdata(i)                ! COVDmat?
          write(15,*) wd(i)                     ! WDmat?
        enddo
        do ikp = 1, nkp
          write(16,100) (smgbig(ikp,inode), inode=1,nnxQ*nnyQ*nnzQ)   ! Gmat?sm
          write(17,*) smdata(ikp)               ! Dmat?sm
          write(18,*) smwd(ikp)                 ! WDmat?sm
        enddo


        close(13)
        close(14)
        close(15)
        close(16)
        close(17)
        close(18)

100   format (20000(f10.4,1x))
110   format (30(i5,1x))
120   format (a5,3f12.4)
130   format (3(f10.4,1x))

      END
