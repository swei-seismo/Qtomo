ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Program to forwardly calculate t* based on a given 1/Q model.      c
c  by S. Wei, May 2014                                                c
c  Input: P/S/PP/SS/PSR (P/S, Syn P/S, Qp/Qs based on Qp)             c
c         ntop (0: no top layer skipped, but add 1 layer at 12 km;    c
c               1: first layer at 12 km skipped;                      c
c               2: first layer at 0 km skipped)                       c
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
      real ylate, xlone, ze, elat, elon, edep, tstar, tster, grade
      real xe0, ye0, ze0, xs0, ys0, zs0, zmoho, tt, aveQinv, ttin 
      real*8 vps(2,nzmax,nymax,nxmax), sn(2,nzmax,nymax,nxmax)    ! velocity, slowness
      real*8 attin(2,nzmax,nymax,nxmax), Qout(2,nzmax), tstarout  ! Syn 1/Q
      real*8 Qpin(nzmax,nymax,nxmax), Qpin1D(nzmax) ! 1/Qp input for Qp/Qs
      real*8 Qtop, toptt
      real grow(maxparam), syntstar, syntsterp, syntsters  ! G-matrix row for given raypath
      real*8 xarc(0:maxpts+1), yarc(0:maxpts+1), zarc(0:maxpts+1)
      real pm, QpQs

      integer nnx, nny, nnz, nnxQ, nnyQ, nnzQ, nnzQ1D, ntop

      character*2 pors, ntops
      character*4 stname(nstmax), stn, lamda
      character*10 phi(nstmax)

      common /a/ x, y, z, xq, yq, zq
      common /b/ ipolar,nnx,nny,nnz,re,pi,r2d,nnxQ,nnyQ,nnzQ,nnzQ1D
      common /c/ sn
      common /d/ attin, zq1D, Qout, Qpin, Qpin1D

c     I/O files
      narguments = iargc()
      if (narguments.lt.1) then
        stop 'Usage: syn_tstar P/S'
      endif
      call getarg(1,pors)
      ntop=0
      if (narguments.eq.2) then
        call getarg(2,ntops)
        read(ntops,'(i1)') ntop
      endif
      if (pors.eq.'P' .or. pors.eq.'p') then
        write(*,*) 'Calculate t*(P)'
        lps=1
        icase=1
        open (8, file='p_tstar.dat')
        open (21,file='p_tstar.cleansyn')     ! output synthetic t*(P)
        open (25,file='raypaths.p')     ! output P ray paths
      else if (pors.eq.'S' .or. pors.eq.'s') then
        write(*,*) 'Calculate t*(S)'
        lps=2
        icase=2
        open (8, file='s_tstar.dat')
        open (21,file='s_tstar.cleansyn')     ! output synthetic t*(S)
        open (25,file='raypaths.s')     ! output S ray paths
      else if (pors.eq.'PS' .or. pors.eq.'ps' .or.
     &         pors.eq.'Ps' .or. pors.eq.'pS' ) then
        write(*,*) 'Calculate t*(S) from Qp and Qp/Qs'
        lps=2
        icase=3
        open (8, file='s_tstar.dat')
        open (21,file='s_tstar.cleansyn')     ! output S ray paths
        open (25,file='raypaths.s')     ! output S ray paths
      else
        stop 'Usage: syn_tstar P/S'
      endif

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
      if (icase.eq.3) then
        Qtop=Qout(1,1)/2.25
      else
        Qtop=Qout(lps,1)
      endif

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

cccccc Read inverted 1/Qp to invert Qp/Qs cccccc
      if (icase.eq.3) then
        open(45,file='Qinv.p')
c        if (ntop.eq.1) then
          do i = 1, nnxQ
            do j = 1, nnyQ
              do k = 1, nnzQ
                read(45,*) Qpin(k,j,i)
              enddo
            enddo
          enddo
c        else if (ntop.eq.0) then
c          do i = 1, nnxQ
c            do j = 1, nnyQ
c              do k = 1, nnzQ+1
c                read(45,*) temp1
c                if (k.gt.1) then
c                  Qpin(k-1,j,i) = temp1
c                endif
c              enddo
c            enddo
c          enddo
c        endif
        do i = 1, nnzQ1D
          read (45,*) Qpin1D(i)
        enddo
        close(45)
      endif

c     Read synthetic 1/Q or Qp/Qs model
      open(44,file='checker.in')
      if ((ntop.eq.1).or.(ntop.eq.2)) then
        do i = 1, nnxQ
          do j = 1, nnyQ
            do k = 1, nnzQ
              attin(1,k,j,i) = 0
              attin(2,k,j,i) = 0
              read(44,*) attin(lps,k,j,i)
c              write(*,*) attin(1,k,j,i)
            enddo
          enddo
        enddo
      else if (ntop.eq.0) then
        do i = 1, nnxQ
          do j = 1, nnyQ
            do k = 1, nnzQ-1
              attin(1,k+1,j,i) = 0
              attin(2,k+1,j,i) = 0
              read(44,*) attin(lps,k+1,j,i)
            enddo
            attin(1,1,j,i) = 1/Qout(1,1)
            attin(2,1,j,i) = 1/Qout(2,1)
          enddo
        enddo
      endif
      close(44)

c     Run main program
c     find number of data
      nrays = 0
20    read (8,*,end=99)
      nrays = nrays+1
      goto 20
99    continue
      rewind 8
      if (icase .eq. 1) then
        nraysp = nrays
      else if (icase. eq. 2) then
        nrayss = nrays
      endif

cccccc Loop over each raypath cccccccccccccc
      iknt = 0
      iskip = 0
      do iray = 1, nrays
        if (icase.eq.1) then
          read (8,*) stn,ylate,xlone,ze,tstar,tster,xx,aveQinv
        elseif ((icase.eq.2).or.(icase.eq.3)) then
          read (8,*) stn,ylate,xlone,ze,tstar,tster,xx,yy,aveQinv,QpQs
        endif
c         find the event location
        if (xlone .lt. 0) xlone = xlone + 360
        call lonlat2xy(lon0,lat0,beta,xlone,ylate,xe0,ye0)
        if (xe0.lt.0 .or. xe0.gt.xmax .or. ye0.lt.0
     &            .or. ye0.gt.ymax) then
c            write (*,*) "ERROR: Earthquake lies outside of the grid",
c     &                        ylate, xlone
          iskip = iskip + 1
          go to 513
        endif
        if (ze .le. 35.) then     ! Skip event shallower than 35 km
          iskip = iskip + 1
c          write (*,*) "ERROR: Earthquake is too shallow", ze
          go to 513
        endif
        if (tster .eq. 0) then
          iskip = iskip + 1
c            write (*,*) "ERROR: t* has wrong error", tster
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
        call bendray(lps,xe0,ye0,ze0,xs0,ys0,zs0,zmoho,maxpts,
     &              xarc,yarc,zarc,npts,tt,iall)
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
        iknt = iknt + 1
c     Output raypath
c          write(21+lps,*) "1234567 ", iknt, aveQinv
c          write(21+lps,*) "111 222 ", iknt,stn, ylate, xlone, ze
c          do pt = 0, npts+2  ! write out raypath
c            write(21+lps,13) xarc(pt), yarc(pt), zarc(pt)
c          enddo

cccccc    Travel time in the top layer
        if (ntop.eq.1) then
          call layertt(icase,maxpts,xarc,yarc,zarc,npts,toptt,ipt1)
cccccc   Add one point: ipt => piercing point at 25 km   cccccc
        elseif ((ntop.eq.0).or.(ntop.eq.2)) then
          toptt=0
          ipt1=npts+2
        endif
c          write(*,*) npts,toptt

c        write(*,*) stn, ylate, xlone, ze
        call sumtstar(icase,maxpts,xarc,yarc,zarc,npts,ipt1,grow,ifskip)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Subroutine that calculate l/(v*Q) for each segment       c
c    id = 1 for P; id = 2 for S                             c
c  INPUT: icase,maxpts,xarc,yarc,zarc,npts                  c
c  OUTPUT: grow                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          write(*,*) "Finish adding row of G matrix"
        if (ifskip .eq. 1) then
          iskip = iskip + 1
c          write (*,*) "SKIP: Earthquake lies outside of the grid",
c     &                        ylate, xlone, ze
          go to 513
        endif
c     Skip INFs and NaNs
        do i = 1, nparm
c          if (.not.(grow(i).le.0.) .and. .not.(grow(i).ge.0.)) then
          if (grow(i) .ne. grow(i)) then
            grow(i) = 0.
c            write(*,*) "Skip NaN"
            go to 513
          endif
          if (abs(grow(i)) .gt. 1.e8) then
            grow(i) = 0.
c            write(*,*) "Skip INF"
            go to 513
          endif
        enddo

C     Sum up each row of G matrix, should be close to travel time
        syntstar = toptt/Qtop
        do i = 1, nparm
          syntstar=syntstar+grow(i)
        enddo
        if (syntstar .eq. 0.) then
          write (*,*) stn,ylate,xlone,ze
        endif
C     OUTPUT
        aveQinv3=1000*syntstar/tt
        if (icase.eq.1) then
          write(21,10) stn,ylate,xlone,ze,syntstar,syntsterp,syntsterp,
     &              aveQinv3
        elseif ((icase.eq.2).or.(icase.eq.3)) then
          write(21,11) stn,ylate,xlone,ze,syntstar,syntsters,syntstar,
     &              syntsterp,aveQinv3,QpQs
        endif
c     Output raypath
        if (lps.eq.1) then
          write(25,*) "1234567 ", iknt, aveQinv3
        elseif (lps.eq.2) then
          write(25,*) "1234567 ", QpQs, aveQinv3
        endif
c          write(21+lps,*) "111 222 ", iknt,stn, ylate, xlone, ze
        do pt = 0, npts+2  ! write out raypath
          write(25,130) xarc(pt), yarc(pt), zarc(pt)
        enddo

513     continue
c          write(*,*) "Finishi ray: ",iknt,stn, ylate, xlone, ze
      enddo
ccccc End of loop over each raypath (iknt)
      ndata = iknt
      write(*,*) "cccccc FINISH CALCULATING t* cccccc"
      write(*,*) "nparm, iknt, skip", nparm, ndata, iskip

      close(21)
      close(25)

10    format (a5,3f12.4,f12.6,2f5.2,f12.6)
11    format (a5,3f12.4,f12.6,f5.2,f12.6,f5.2,f12.6,f7.2)
130   format (3(f10.4,1x))

      END
