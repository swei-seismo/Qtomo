ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Subroutine to calculate the travel time in a layer                c
c   id = 1 for P; id = 2 for S; id = 3 for Qp/Qs from t*(S) and Qp    c
c   Add one point: ipt => piercing point at 25 km                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE layertt(id,maxpts,xarc,yarc,zarc,npts,tt,ipt)
      parameter (nxmax = 45, nymax = 45, nzmax = 30)  ! max nodes

      integer maxpts, pt, ipt, ifskip, casep       ! max # internal points in raypath
      integer nnx, nny, nnz, nnxQ, nnyQ, nnzQ, nnzQ1D, nparm
      real x(nxmax), z(nzmax), y(nymax), xq(nxmax), zq(nzmax), yq(nymax)  ! node locations
      real*8 sn(2,nzmax,nymax,nxmax)
      real*8 xarc(0:maxpts+2), yarc(0:maxpts+2), zarc(0:maxpts+2)
      real*8 savexarc(0:maxpts+2), saveyarc(0:maxpts+2)
      real*8 savezarc(0:maxpts+2)
      real*8 xp1,xp2,xp3,xp4,xp5,yp1,yp2,yp3,yp4,yp5,zp1,zp2,zp3,zp4,zp5
      real*8 wc(8), wd(8), slc, wtc                  ! weighting coefficients
      real*8 xlen, ylen, zlen, sslen, tt, toplen
      real sslength, xlength, ylength, zlength

      common /a/ x, y, z, xq, yq, zq
      common /b/ ipolar,nnx,nny,nnz,re,pi,r2d,nnxQ,nnyQ,nnzQ,nnzQ1D
      common /c/ sn

c     Find the piercing point of the first V layer
      ztop = z(2)
      ipt=0
      do i=0,npts+2
        if (zarc(i).lt.ztop .and. zarc(i-1).gt.ztop) then
          zslope = (zarc(i-1)-z(2))/(zarc(i-1)-zarc(i))
          xtop = xarc(i-1)-zslope*(xarc(i-1)-xarc(i))
          ytop = yarc(i-1)-zslope*(yarc(i-1)-yarc(i))
          ipt = i
          continue
        endif
      enddo
c      write(*,200) ipt,npts,zarc(ipt),zarc(ipt-1),zarc(npts+2)
c200   format (2(i4,1x),3(f10.4,1x))
      do i = 0, npts+2
        savexarc(i) = xarc(i)
        saveyarc(i) = yarc(i)
        savezarc(i) = zarc(i)
      enddo
      npts = npts + 1
      do i = ipt+1, npts+2
        xarc(i) = savexarc(i-1)
        yarc(i) = saveyarc(i-1)
        zarc(i) = savezarc(i-1)
      enddo
      xarc(ipt) = xtop
      yarc(ipt) = ytop
      zarc(ipt) = ztop

      tt = 0.
      toplen = 0.
      x0 = xarc(ipt)
      y0 = yarc(ipt)
      z0 = zarc(ipt)
      x1 = xarc(ipt+1)
      y1 = yarc(ipt+1)
      z1 = zarc(ipt+1)
      pt = ipt    ! End of the segment
      do while (pt .le. npts+2)
        sign = +1.            ! ray moving upwards
        if (z1 .gt. z0) sign = -1.
        if (x1 .eq. x0) then
          xslope = 100000.
          yxslope = 100000.
        else
          xslope = (z1-z0)/(x1-x0)    ! positive means up & left or down & right
          yxslope = (y1-y0)/(x1-x0)
        endif
        if (y1 .eq. y0) then
          yslope = 100000.
          xyslope = 100000.
        else
          yslope = (z1-z0)/(y1-y0)
          xyslope = (x1-x0)/(y1-y0)
        endif
        b = z1 - xslope*x1
        c = z1 - yslope*y1
        do i = 1, nnx   ! i0: node nearest 0 in x direction
          if (x(i).lt.x0 .and. x(i+1).ge.x0 .and. xslope.ge.0.) i0 = i
          if (x(i).le.x0 .and. x(i+1).gt.x0 .and. xslope.lt.0.) i0 = i
        enddo
        if (x0 .eq. x(i0)   .and. sign*xslope .gt. 0.) i0 = i0-1
        if (x0 .eq. x(i0+1) .and. sign*xslope .lt. 0.) i0 = i0+1
        if (xslope*sign .ge. 0.) xi = x(i0)
        if (xslope*sign .lt. 0.) xi = x(i0+1)
        if (x1 .eq. x0) then
          lenz0=10000
          do ii = 1,nnz
            if (((z(ii)-z0)*(z(ii)-z1) .le. 0) .and.
     &            ((abs(z(ii)-z0)) .ge. (abs(z(ii)-z1))) .and.
     &            (abs(z(ii)-z0) .lt. lenz0)) then
              lenz0=abs(z(ii)-z0)
              zi = z(ii) ! between z1 and z0 but close to z0
            endif
          enddo
          leny0=10000
          do ii = 1,nny
            if (((y(ii)-y0)*(y(ii)-y1) .le. 0) .and.
     &            ((abs(y(ii)-y0)) .ge. (abs(y(ii)-y1))) .and.
     &            (abs(y(ii)-y0) .lt. leny0)) then
              leny0=abs(y(ii)-y0)
              yi = y(ii) ! between y1 and y0 but close to y0
            endif
          enddo
        else
          zi = xslope*xi + b
          yi = y0 + yxslope*(xi-x0)
        endif
        do j = 1, nny   ! j0: node nearest 0 in y direction
          if (y(j).lt.y0 .and. y(j+1).ge.y0 .and. yslope.ge.0.) j0 = j
          if (y(j).le.y0 .and. y(j+1).gt.y0 .and. yslope.lt.0.) j0 = j
        enddo
        if (y0 .eq. y(j0)   .and. sign*yslope .gt. 0.) j0 = j0-1
        if (y0 .eq. y(j0+1) .and. sign*yslope .lt. 0.) j0 = j0+1
        if (yslope*sign .ge. 0.) yj = y(j0)
        if (yslope*sign .lt. 0.) yj = y(j0+1)
        if (y1 .eq. y0) then
          lenz0=10000
          do ii = 1,nnz
            if (((z(ii)-z0)*(z(ii)-z1) .le. 0) .and.
     &            ((abs(z(ii)-z0)) .ge. (abs(z(ii)-z1))) .and.
     &            (abs(z(ii)-z0) .lt. lenz0)) then
              lenz0=abs(z(ii)-z0)
              zj = z(ii) ! between z1 and z0 but close to z0
            endif
          enddo
          lenx0=10000
          do ii = 1,nny
            if (((x(ii)-x0)*(x(ii)-x1) .le. 0) .and.
     &            ((abs(x(ii)-x0)) .ge. (abs(x(ii)-x1))) .and.
     &            (abs(x(ii)-x0) .lt. lenx0)) then
              leny0=abs(x(ii)-x0)
              xj = x(ii) ! between x1 and x0 but close to x0
            endif
          enddo
        else
          zj = yslope*yj + c
          xj = x0 + xyslope*(yj-y0)
        endif
        do k = 1, nnz   ! k0: node nearest 0 in z direction
          if (z(k) .lt. z0 .and. z(k+1) .ge. z0) k0 = k
        enddo
        if (z0 .eq. z(k0+1) .and. sign .lt. 0.) k0 = k0+1
        if (sign .gt. 0.) zk = z(k0)
        if (sign .lt. 0.) zk = z(k0+1)
        if (z1 .eq. z0) then
          lenx0=10000
          do ii = 1,nny
            if (((x(ii)-x0)*(x(ii)-x1) .le. 0) .and.
     &            ((abs(x(ii)-x0)) .ge. (abs(x(ii)-x1))) .and.
     &            (abs(x(ii)-x0) .lt. lenx0)) then
              leny0=abs(x(ii)-x0)
              xk = x(ii) ! between x1 and x0 but close to x0
            endif
          enddo
          leny0=10000
          do ii = 1,nny
            if (((y(ii)-y0)*(y(ii)-y1) .le. 0) .and.
     &            ((abs(y(ii)-y0)) .ge. (abs(y(ii)-y1))) .and.
     &            (abs(y(ii)-y0) .lt. leny0)) then
              leny0=abs(y(ii)-y0)
              yk = y(ii) ! between y1 and y0 but close to y0
            endif
          enddo
        else
          xk = (zk - b)/xslope
          yk = y0 + yxslope*(xk-x0)
        endif
        case1length = (x1-x0)**2+(y1-y0)**2+(z1-z0)**2 ! dist between 1 and 0
        xlength = (xi-x0)**2+(yi-y0)**2+(zi-z0)**2  ! dist between i0 and 0
        ylength = (xj-x0)**2+(yj-y0)**2+(zj-z0)**2  ! dist between j0 and 0
        zlength = (xk-x0)**2+(yk-y0)**2+(zk-z0)**2  ! dist between k0 and 0
        icase = 1       ! subsegment ends at arc node
        sslength = case1length
        if (xlength .le. sslength) then   ! subsegment ends at x(i)
          icase = 2
          sslength = xlength
        endif
        if (ylength .le. sslength) then    ! subsegment ends at y(i)
          icase = 3
          sslength = ylength
        endif
        if (zlength .le. sslength) then    ! subsegment ends at z(j)
          icase = 4
        endif
c        write(*,*) icase, sslength, xlength, ylength, zlength
        if (icase .eq. 1) then            ! set subsegment end location
          xss = x1
          yss = y1
          zss = z1
        else if (icase .eq. 2) then
          xss = xi
          yss = yi
          zss = zi
        else if (icase .eq. 3) then
          xss = xj
          yss = yj
          zss = zj
        else if (icase .eq. 4) then
          xss = xk
          yss = yk
          zss = zk
        else
          write (*,*) "icase out of bounds, (", icase, ")"
          write (*,*) 'x1,y1,z1,xi,yi,zi,xj,yj,zj,xk,yk,zk'
          write (*,*) x1,y1,z1,xi,yi,zi,xj,yj,zj,xk,yk,zk
          stop
        endif
        xlen = xss-x0
        ylen = yss-y0
        zlen = zss-z0
        if (ipolar .eq. 1) then     ! x,y distances decr w/ depth
          xlen = xlen*(1.-(zss+z0)/(2.*re))
          ylen = ylen*(1.-(zss+z0)/(2.*re))
        endif
        sslen = sqrt((xlen**2) + (ylen**2) + (zlen**2))
        toplen = toplen + sslen
c        write(*,*) "Length of the subsegment is ",sslen,k0
c        write(*,*) x0,y0,z0,sn(1,k0,j0,i0),pt,npts+2,ipt
ccccc Adding grow elements corresponding to 1D Q model from zq(nnzQ) to z(nnz)
        if (id .eq. 1) then          ! Qp
          tt = tt + sslen*sn(1,k0,j0,i0)
        else if ((id .eq. 2) .or. (id .eq. 3)) then     ! Qs
          tt = tt + sslen*sn(2,k0,j0,i0)
        endif
        x0 = xss
        y0 = yss
        z0 = zss
        if (z0 .eq. zarc(npts+2)) then
          pt = pt+1
        elseif (icase .eq. 1) then      ! Move to next subsegment for case 1
          pt = pt + 1
          x1 = xarc(pt)
          y1 = yarc(pt)
          z1 = zarc(pt)
        endif
c        write(*,*) nnxQ*nnyQ*nnzQ+k0-nnzQ, grow(nnxQ*nnyQ*nnzQ+k0-nnzQ)
      enddo
c      tstarout=0.
c      do i=1,nnz-nnzQ
c        tstarout=tstarout+grow(nnxQ*nnyQ*nnzQ+i)
c      enddo
      if ((toplen+zarc(npts+2)) .lt. 25.) then
        write(*,100) tt, toplen, zarc(npts+1),zarc(npts+2),
     &            (toplen+zarc(npts+2)-25.)
      endif
100   format (5(f10.4,1x))
      END
