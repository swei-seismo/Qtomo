ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Forwardly calculate t*(S)/t*(P) based on a given Qp/Qs model       c
c  by S. Wei, August 2015                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      PROGRAM MAIN

      parameter (nstmax = 100)     ! max number of stations
      parameter (nxmax = 45, nymax = 45, nzmax = 30)  ! max nodes
      parameter (maxpts = 2048)       ! max # internal points in raypath
      parameter (maxparam = nxmax*nymax*nzmax+nzmax)      ! max velocity grid
      parameter (maxzq = 700) ! max depth of Q model

      real lat0, lon0, beta   ! refrence node
      real lat(nxmax,nymax), lon(nxmax,nymax), temp1, temp2
      real latq(nxmax,nymax), lonq(nxmax,nymax)
      real x(nxmax), y(nymax), z(nzmax), xmax, ymax         ! node locations
      real xq(nxmax), zq(nzmax), yq(nymax), zq1D(nzmax)
      real stlon(nstmax), stlat(nstmax), scrust(nstmax)     ! station locations
      real xsta(nstmax), ysta(nstmax), zsta(nstmax)
      real ylate, xlone, ze, elat, elon, edep,ststar,ptstar,steer,pteer
      real xe0, ye0, ze0, xs0, ys0, zs0, zmoho, tt, aveQinv, ttin 
      real*8 vps(2,nzmax,nymax,nxmax), sn(2,nzmax,nymax,nxmax)    ! velocity, slowness
      real*8 attin(2,nzmax,nymax,nxmax), Qout(2,nzmax), tstarout     ! 1/Qp input for Qp/Qs
      real*8 Qpin(nzmax,nymax,nxmax), Qpin1D(nzmax) ! 1/Qp input for Qp/Qs
      real*8 Qtop, topttp, toptts
      real growp(maxparam),grows(maxparam)  ! G-matrix row for given raypath
      real syntstarp, syntstars, syntsterp, syntsters
      real*8 xarcp(0:maxpts+2), yarcp(0:maxpts+2), zarcp(0:maxpts+2)
      real*8 xarcs(0:maxpts+2), yarcs(0:maxpts+2), zarcs(0:maxpts+2)
      real pm, QpQs, xx, yy

      integer nnx, nny, nnz, nnxQ, nnyQ, nnzQ, nnzQ1D
      integer p, pt, ifskipp, ifskips

      character*2 pors
      character*4 stname(nstmax), stn, lamda
      character*10 phi(nstmax)

      common /a/ x, y, z, xq, yq, zq
      common /b/ ipolar,nnx,nny,nnz,re,pi,r2d,nnxQ,nnyQ,nnzQ,nnzQ1D
      common /c/ sn
      common /d/ attin, zq1D, Qout, Qpin, Qpin1D

c     I/O files
      write(*,*) 'Calculate t*(S)/t*(P)'
      open (8, file='s_tstar.olddat')
      open (21,file='s_tstar.cleansyn')     ! output synthetic t*(S)/t*(P)

      open (9, file='station.lst')  ! station lat, lon, z, crust, name
      open (10,file='lau.vel')     ! 3-D velocity model
      open (11,file='Q.grid')     ! 3-D velocity model
      open (12,file='Q1D.grid')   ! 1-D Q model
      open (13,file='Q.prem')     ! 1-D Q values of PREM
c      open (30,file='event.lst')    ! output used earthquake location
c      elon=0
c      elat=0
c      edep=0

c     Other parameters
      smooth = 8000.
      alpha = 35

      pi = 4.*atan(1.)
      r2d = 180./pi
      re = 6377.        ! radius of Earth ~ 20 deg lat
      ipolar = 1        ! flat Earth = 0, curved Earth = 1
      iall = 0
      syntsterp = 0.2
      syntsters = 0.3


C     velocity model input
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
c      open (27,file='cartesian.grid')
      inode = 0
      write(27,*) xmax, ymax
      do i = 1, nnx
        do j = 1, nny
          inode = inode+1
c          write(27,*) x(i), y(j), inode
        enddo
      enddo
      close(27)

ccccc Input Q values of PREM for reference outside of Q model
      read(13,*) (Qout(1,k), Qout(2,k), k=1,nnz)
      close(13)
c      Qtop=Qout(lps,1)

ccccc Q grid input ccccc
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
      read (11,*) (zq(k), k=1,nnzQ)   ! positive down
      close(11)
ccccc 1D Q model input ccccc
      read (12,*) nnzQ1D
      do i = 1, nnzQ1D
        read (12,*) zq1D(i)
      enddo
      close(12)
      nparm = nnxQ*nnyQ*nnzQ + nnzQ1D
      write(*,*) "nnxQ, nnyQ, nnzQ, nparm: ", nnxQ, nnyQ, nnzQ, nparm

C     Input station info
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

c     Read synthetic 1/Q model
      open(44,file='checkerp.in')
      do i = 1, nnxQ
        do j = 1, nnyQ
          do k = 1, nnzQ
            read(44,*) attin(1,k,j,i)
          enddo
        enddo
      enddo
      close(44)
      open(44,file='checkers.in')
      do i = 1, nnxQ
        do j = 1, nnyQ
          do k = 1, nnzQ
            read(44,*) attin(2,k,j,i)
          enddo
        enddo
      enddo
      close(44)

c     Run main program
c     find number of data
      nrays = 0
20    read (8,*,end=99)
      nrays = nrays+1
      goto 20
99    continue
      rewind 8

cccccc Loop over each raypath cccccccccccccc
      iknt = 0
      iskip = 0
      do iray = 1, nrays
        read (8,*) stn,ylate,xlone,ze,ststar,steer,ptstar,pteer,
     &              aveQinv,QpQs
        if (ptstar .eq. 0) then
c          write (*,*) "ERROR: t*(P) = 0"
          iskip = iskip + 1
          go to 513
        endif
c         find the event location
        if (xlone .lt. 0) xlone = xlone + 360
        call lonlat2xy(lon0,lat0,beta,xlone,ylate,xe0,ye0)
        if (xe0 .lt. 0 .or. xe0 .gt. xmax .or. ye0 .lt. 0
     &            .or. ye0 .gt. ymax) then
c            write (*,*) "ERROR: Earthquake lies outside of the grid",
c     &                        ylate, xlone
          iskip = iskip + 1
          go to 513
        endif
        if (ze .le. 35.) then     ! Skip event shallower than 35 km
          iskip = iskip + 1
          write (*,*) "ERROR: Earthquake is too shallow", ze
          go to 513
        endif
        if (steer .eq. 0) then
          iskip = iskip + 1
          write (*,*) "ERROR: t* has wrong error", tster
          go to 513
        endif
        ze0 = ze
c       find the station
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
c        write (*,*) lps,xe0,ye0,ze0,xs0,ys0,zs0,zmoho,maxpts
        call bendray(1,xe0,ye0,ze0,xs0,ys0,zs0,zmoho,maxpts,
     &              xarcp,yarcp,zarcp,nptsp,ttp,iall)
        call bendray(2,xe0,ye0,ze0,xs0,ys0,zs0,zmoho,maxpts,
     &              xarcs,yarcs,zarcs,nptss,tts,iall)
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
c          write(*,*) "Finish ray tracing"
        ttin=1000*tstar/aveQinv
c          if (abs(tt-ttin) .gt. 5) then
c            write(*,*) "Warning: ray tracing may fail for ",
c     $                  stn, ylate, xlone, ze
c            write(*,*) "Ray tracing tt - Antelope tt: ", (tt-ttin)
c            write(*,*) "Raytracing: ",tt," Antelope: ",ttin
c          endif
c     Output raypath
c          write(21+lps,*) "1234567 ", iknt, aveQinv
c          write(21+lps,*) "111 222 ", iknt,stn, ylate, xlone, ze
c          do pt = 0, npts+2  ! write out raypath
c            write(21+lps,13) xarc(pt), yarc(pt), zarc(pt)
c          enddo
        nptss=nptsp
        do i = 0, nptss+2
          xarcs(i)=xarcp(i)
          yarcs(i)=yarcp(i)
          zarcs(i)=zarcp(i)
        enddo

cccccc    Travel time in the top layer
        call layertt(1,maxpts,xarcp,yarcp,zarcp,nptsp,topttp,iptp)
        call layertt(2,maxpts,xarcs,yarcs,zarcs,nptss,toptts,ipts)
cccccc   Add one point: ipt => piercing point at 25 km   cccccc
c          write(*,*) npts,toptt

c        write(*,*) stn, ylate, xlone, ze
        call sumtstar(1,maxpts,xarcp,yarcp,zarcp,nptsp,iptp,
     &                  growp,ifskipp)
        call sumtstar(2,maxpts,xarcs,yarcs,zarcs,nptss,ipts,
     &                  grows,ifskips)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Subroutine that calculate l/(v*Q) for each segment       c
c    id = 1 for P; id = 2 for S                             c
c  INPUT: icase,maxpts,xarc,yarc,zarc,npts                  c
c  OUTPUT: grow                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          write(*,*) "Finish adding row of G matrix"
        if ((ifskipp.eq.1).or.(ifskips.eq.1)) then
          iskip = iskip + 1
c          write (*,*) "SKIP: Earthquake lies outside of the grid",
c     &                        ylate, xlone, ze
          go to 513
        endif
c     Skip INFs and NaNs
        do i = 1, nparm
          if (growp(i) .ne. growp(i)) then
            growp(i) = 0.
c              write(*,*) "Skip NaN"
            iskip = iskip + 1
            go to 513
          endif
          if (abs(growp(i)) .gt. 1.e8) then
            growp(i) = 0.
c              write(*,*) "Skip INF"
            iskip = iskip + 1
            go to 513
          endif
        enddo
        do i = 1, nparm
          if (grows(i) .ne. grows(i)) then
            grows(i) = 0.
            iskip = iskip + 1
            go to 513
          endif
          if (abs(grows(i)) .gt. 1.e8) then
            grows(i) = 0.
            iskip = iskip + 1
            go to 513
          endif
        enddo

        iknt = iknt + 1
C     Sum up each row of G matrix, should be close to travel time
        syntstarp = topttp/Qout(1,1)
        do i = 1, nparm
          syntstarp=syntstarp+growp(i)
        enddo
        if (syntstarp .eq. 0.) then
          write (*,*) stn,ylate,xlone,ze
        endif
        syntstars = toptts/Qout(2,1)
        do i = 1, nparm
          syntstars=syntstars+grows(i)
        enddo
        if (syntstars .eq. 0.) then
          write (*,*) stn,ylate,xlone,ze
        endif
C     OUTPUT
        aveQinv3=1000*syntstars/tts
        QpQs = (syntstars*ttp)/(syntstarp*tts)
        write(21,11) stn,ylate,xlone,ze,syntstars,syntsters,syntstarp,
     &              syntsterp,aveQinv3,QpQs
513     continue
c          write(*,*) "Finishi ray: ",iknt,stn, ylate, xlone, ze
      enddo
ccccc End of loop over each raypath (iknt)
      ndata = iknt
      write(*,*) "cccccc FINISH CALCULATING t* cccccc"
      write(*,*) "nparm, iknt, skip", nparm, ndata, iskip

      close(21)

10    format (a5,3f12.4,f12.6,2f5.2,f12.6)
11    format (a5,3f12.4,f12.6,f5.2,f12.6,f5.2,f12.6,f7.2)

      END
