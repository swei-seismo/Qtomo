ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Subroutine that builds G-matrix row for given raypath(s)          c
c   id = 1 for P; id = 2 for S; id = 3 for Qp/Qs from t*(S) and Qp    c
c   id = 4 for bulk attenuation t*(P) and t*(S)                       c
c   Add one point: ipt => piercing point into the Q model             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE getgrowQ(id,maxpts,xarc,yarc,zarc,npts,ipt1,
     &                      grow,ifskip)
      parameter (nxmax = 45, nymax = 45, nzmax = 30)  ! max nodes

      integer maxpts, pt, ipt, ifskip, casep      ! max # internal points in raypath
      integer nnx, nny, nnz, nnxQ, nnyQ, nnzQ, nnzQ1D, nparm, lps
      real x(nxmax), z(nzmax), y(nymax)
      real xq(nxmax), zq(nzmax), yq(nymax), zq1D(nzmax), attin1D(nzmax)
      real*8 sn(2,nzmax,nymax,nxmax)
      real*8 attin(nzmax,nymax,nxmax), Qout(2,nzmax)     ! 1/Qp input for Qp/Qs
      real*8 xarc(0:maxpts+2), yarc(0:maxpts+2), zarc(0:maxpts+2)
      real*8 savexarc(0:maxpts+2), saveyarc(0:maxpts+2)
      real*8 savezarc(0:maxpts+2)
      real*8 xp1,xp2,xp3,xp4,xp5,yp1,yp2,yp3,yp4,yp5,zp1,zp2,zp3,zp4,zp5
      real*8 wc(8), wd(8), slc, wtc, snQ(8), snQ2(8)     ! weighting coefficients
      real*8 xlen, ylen, zlen, sslen, Qx, Qy, Qz
      real grow(nxmax*nymax*nzmax+nzmax) ! G-matrix row for given raypaths
      real sslength, xlength, ylength, zlength, alllen

      common /a/ x, y, z, xq, yq, zq
      common /b/ ipolar,nnx,nny,nnz,re,pi,r2d,nnxQ,nnyQ,nnzQ,nnzQ1D
      common /c/ sn
      common /d/ attin, zq1D, attin1D

      ifskip = 0    ! If skip this event because most part of ray is outside
c     check if raypath is in the velocity and Q model space
      do pt = 0, npts+1
        if ((xarc(pt) .lt. x(1)).or.(xarc(pt) .gt. x(nnx)).or.
     &      (zarc(pt) .lt. z(1)).or.(zarc(pt) .gt. z(nnz)).or.
     &      (yarc(pt) .lt. y(1)).or.(yarc(pt) .gt. y(nny))) then
          write (*,*) "ERROR: ray is not in the velocity model space"
          write (*,*) pt, xarc(pt), yarc(pt), zarc(pt)
c          stop
          ifskip = 1
          goto 237
        endif
      enddo

c     Find the first point of the ray within the Q model
      ipt=0
      do i = 0, npts+2
        if ((xarc(i) .le. xq(nnxQ)).and.(xarc(i) .ge. xq(1)).and.
     &      (yarc(i) .le. yq(nnyQ)).and.(yarc(i) .ge. yq(1)).and.
     &      (zarc(i) .le. zq(nnzQ)).and.(zarc(i) .ge. zq(1))) then
          ipt = i
          goto 236
        endif
        ipt = npts+2
      enddo
236   continue
      if (ipt.eq.npts+2) then
        ifskip = 1
        goto 237
      endif
c      write(*,'(6f10.4)') xq(1),xq(nnxQ),yq(1),yq(nnyQ),zq(1),zq(nnzQ)
c      write(*,*) ipt,npts
c      write(*,300) "Earthquake:  ",xarc(0),yarc(0),zarc(0)
c      write(*,300) "First point: ",xarc(ipt),yarc(ipt),zarc(ipt)
c      write(*,300) "Point before:",xarc(ipt-1),yarc(ipt-1),zarc(ipt-1)
      do i = 0, npts+2
        savexarc(i) = xarc(i)
        saveyarc(i) = yarc(i)
        savezarc(i) = zarc(i)
      enddo

ccccc Initiating the row of G matrix
      nparm = nnxQ*nnyQ*nnzQ+nnzQ1D
c      write(*,*) "nparm = ",nparm
      do i = 1, nparm
        grow(i) = 0.
      enddo

ccccc If the whole ray is within the Q model space, jump to buiding G
      if (ipt.eq.0) then
c        ipt = 1
        goto 235
      endif

c      alllen=0.
cccc  -- build G-matrix row (grow) for arc segment in 1D Q model-- cccc
c     find and insert the piercing point
      casep=0
      if (zarc(ipt-1).gt.zq(nnzQ)) then   ! Piercing from bottom
        zp1 = zq(nnzQ)
        xp1 = xarc(ipt)-(zarc(ipt)-zq(nnzQ))/(zarc(ipt)-zarc(ipt-1))
     &                  *(xarc(ipt)-xarc(ipt-1))
        yp1 = yarc(ipt)-(zarc(ipt)-zq(nnzQ))/(zarc(ipt)-zarc(ipt-1))
     &                  *(yarc(ipt)-yarc(ipt-1))
        if ((xp1 .ge. xq(1)).and.(xp1 .le. xq(nnxQ)).and.
     &            (yp1 .ge. yq(1)).and.(yp1 .le. yq(nnyQ))) casep=1
      endif
      if (xarc(ipt-1).lt.xq(1)) then   ! Piercing from left
        xp2 = zq(1)
        zp2 = zarc(ipt)-(xarc(ipt)-xq(1))/(xarc(ipt)-xarc(ipt-1))
     &                  *(zarc(ipt)-zarc(ipt-1))
        yp2 = yarc(ipt)-(xarc(ipt)-xq(1))/(xarc(ipt)-xarc(ipt-1))
     &                  *(yarc(ipt)-yarc(ipt-1))
        if ((zp2 .ge. zq(1)).and.(zp2 .le. zq(nnzQ)).and.
     &            (yp2 .ge. yq(1)).and.(yp2 .le. yq(nnyQ))) casep=2
      endif
      if (xarc(ipt-1).gt.xq(nnxQ)) then   ! Piercing from right
        xp3 = xq(nnxQ)
        zp3 = zarc(ipt)-(xarc(ipt)-xq(nnxQ))/(xarc(ipt)-xarc(ipt-1))
     &                  *(zarc(ipt)-zarc(ipt-1))
        yp3 = yarc(ipt)-(xarc(ipt)-xq(nnxQ))/(xarc(ipt)-xarc(ipt-1))
     &                  *(yarc(ipt)-yarc(ipt-1))
        if ((zp3 .ge. zq(1)).and.(zp3 .le. zq(nnzQ)).and.
     &            (yp3 .ge. yq(1)).and.(yp3 .le. yq(nnyQ))) casep=3
      endif
      if (yarc(ipt-1).lt.yq(1)) then   ! Piercing from front
        yp4 = yq(1)
        zp4 = zarc(ipt)-(yarc(ipt)-yq(1))/(yarc(ipt)-yarc(ipt-1))
     &                  *(zarc(ipt)-zarc(ipt-1))
        xp4 = xarc(ipt)-(yarc(ipt)-yq(1))/(yarc(ipt)-yarc(ipt-1))
     &                  *(xarc(ipt)-xarc(ipt-1))
        if ((zp4 .ge. zq(1)).and.(zp4 .le. zq(nnzQ)).and.
     &            (xp4 .ge. xq(1)).and.(xp4 .le. xq(nnxQ))) casep=4
      endif
      if (yarc(ipt-1).gt.yq(nnyQ)) then   ! Piercing from back
        yp5 = yq(nnyQ)
        zp5 = zarc(ipt)-(yarc(ipt)-yq(nnyQ))/(yarc(ipt)-yarc(ipt-1))
     &                  *(zarc(ipt)-zarc(ipt-1))
        xp5 = xarc(ipt)-(yarc(ipt)-yq(nnyQ))/(yarc(ipt)-yarc(ipt-1))
     &                  *(xarc(ipt)-xarc(ipt-1))
        if ((zp5 .ge. zq(1)).and.(zp5 .le. zq(nnzQ)).and.
     &            (xp5 .ge. xq(1)).and.(xp5 .le. xq(nnxQ))) casep=5
      endif
c      write(*,*) "casep = ",casep
      npts = npts + 1
      ipt1 = ipt1 + 1
      do i = ipt+1, npts+2
        xarc(i) = savexarc(i-1)
        yarc(i) = saveyarc(i-1)
        zarc(i) = savezarc(i-1)
      enddo
      if (casep.eq.0) then
        write(*,*) "Warning: fail to find the piercing point"
        write(*,200) xarc(ipt),yarc(ipt),zarc(ipt)
        write(*,200) xarc(ipt-1),yarc(ipt-1),zarc(ipt-1)
        ifskip = 1
        goto 237
      elseif (casep.eq.1) then
        xarc(ipt) = xp1
        yarc(ipt) = yp1
        zarc(ipt) = zp1
      elseif (casep.eq.2) then
        xarc(ipt) = xp2
        yarc(ipt) = yp2
        zarc(ipt) = zp2
      elseif (casep.eq.3) then
        xarc(ipt) = xp3
        yarc(ipt) = yp3
        zarc(ipt) = zp3
      elseif (casep.eq.4) then
        xarc(ipt) = xp4
        yarc(ipt) = yp4
        zarc(ipt) = zp4
      elseif (casep.eq.5) then
        xarc(ipt) = xp5
        yarc(ipt) = yp5
        zarc(ipt) = zp5
      endif
c      write(*,*) xarc(ipt),yarc(ipt),zarc(ipt),casep,ipt
c      write(*,*) "One ray"

      if (zarc(ipt).lt.zq(nnzQ)) then
c        write(*,*) "Skip outside event piercing too shallow"
c        write(*,*) zarc(ipt-1), zarc(ipt), zarc(ipt+1)
        ifskip = 1
        goto 237
      endif
      x0 = xarc(0)
      y0 = yarc(0)
      z0 = zarc(0)
      x1 = xarc(1)
      y1 = yarc(1)
      z1 = zarc(1)
      pt = 1    ! End of the segment
c      write (*,*) (zarc(i),i=0,ipt)
      do while (pt .le. ipt)
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
c        do k = 1, nnzQ
c          if (zq(k) .lt. z0 .and. zq(k+1) .ge. z0) kq0 = k
c        enddo
        if (z0 .eq. z(k0+1) .and. sign .lt. 0.) k0 = k0+1
c     Determine k1D from 1 to nnzQ1D
c        if (k0 .le. nnzQ) then
c          k1D = nnzQ+1
c        else
c          k1D = k0
c        endif
        k1D = 1
        do k = 1, nnzQ1D
          if (zq1D(k) .lt. z0 .and. zq1D(k+1) .ge. z0) k1D = k
        enddo
        if (z0 .eq. zq1D(k1D+1) .and. sign .lt. 0.) k1D = k1D + 1
c        print *,z0,zq1D(k1D),zq1D(k1D+1),k1D
c        if (z0 .eq. zq(kq0+1) .and. sign .lt. 0.) kq0 = kq0+1
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
        if (ylength .le. sslength) then    ! subsegment ends at y(j)
          icase = 3
          sslength = ylength
        endif
        if (zlength .le. sslength) then    ! subsegment ends at z(k)
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
c        write(*,*) "Length of the subsegment is ",sslen,icase
c        alllen=alllen+sslen
c        write(*,*) x0,y0,z0,sn(id,k0,j0,i0)
ccccc Adding grow elements corresponding to 1D Q model from 1 to nnzQ1D
        if ((id .eq. 1) .or. (id .eq. 2)) then          ! Qp or Qs
          grow(nnxQ*nnyQ*nnzQ+k1D) = grow(nnxQ*nnyQ*nnzQ+k1D)+
     &                                    sslen*sn(id,k0,j0,i0)
        elseif (id .eq. 3) then    ! Qp/Qs
          grow(nnxQ*nnyQ*nnzQ+k1D) = grow(nnxQ*nnyQ*nnzQ+k1D)+
     &                                sslen*sn(2,k0,j0,i0)*attin1D(k1D)
C       elseif (id .eq. 4) then    ! Qk
C         grow(nnxQ*nnyQ*nnzQ+k1D) = grow(nnxQ*nnyQ*nnzQ+k1D)+
C    &                  1-3./4.*(sn(2,k0,j0,i0)**2)/(sn(1,k0,j0,i0)**2)
        endif
        x0 = xss
        y0 = yss
        z0 = zss
        if (z0 .eq. zarc(ipt)) then
          pt = pt+1
        elseif (icase .eq. 1) then      ! Move to next subsegment for case 1
          pt = pt + 1
          x1 = xarc(pt)
          y1 = yarc(pt)
          z1 = zarc(pt)
        endif
c        write(*,*) nnxQ*nnyQ*nnzQ+k0-nnzQ, grow(nnxQ*nnyQ*nnzQ+k0-nnzQ)
      enddo
c      write(*,*) "t* outside = ", tstarout
c      write(*,*) grow(nnxQ*nnyQ*nnzQ+1)
cccc  -- build G-matrix row (grow) for arc segment in 1D Q model-- cccc

cccc  -- build G-matrix row (grow) for arc segment in 3D Q model-- cccc
235   continue
c      write(*,*) 'npts = ', npts
      x0 = xarc(ipt)
      y0 = yarc(ipt)
      z0 = zarc(ipt)
      x1 = xarc(ipt+1)
      y1 = yarc(ipt+1)
      z1 = zarc(ipt+1)
c      write(*,*) x0,y0,z0
c      write(*,*) x1,y1,z1
c      write(*,'(6f10.4)') xq(1),xq(nnxQ),yq(1),yq(nnyQ),zq(1),zq(nnzQ)
c      write(*,*) "Adding G-matrix", ipt1,ipt
cccc  Start segments from bottom
      pt = ipt
      do while (pt .le. ipt1)
c        write(*,*) pt, x1, x0, y1, y0, z1, z0
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
c        write(*,*) "yxslope, xslope, xyslope, yslope",
c     &            yxslope, xslope, xyslope, yslope
         ! sign*xslope => ray moving right-left (+left, -right) 
         ! sign*yslope => ray moving front-back (+front, -back)

        b = z1 - xslope*x1
        c = z1 - yslope*y1

c        do i = 1, nnx   ! i0: node nearest 0 in x direction
c          if (x(i).lt.x0 .and. x(i+1).ge.x0 .and. xslope.ge.0.) i0 = i
c          if (x(i).le.x0 .and. x(i+1).gt.x0 .and. xslope.lt.0.) i0 = i
c        enddo
        iq0=0
        do i = 1, nnxQ
          if (xq(i).lt.x0 .and. xq(i+1).ge.x0 .and. xslope.ge.0.) iq0=i
          if (xq(i).le.x0 .and. xq(i+1).gt.x0 .and. xslope.lt.0.) iq0=i
        enddo
c        if (x0 .eq. x(i0)   .and. sign*xslope .gt. 0.) i0 = i0-1
c        if (x0 .eq. x(i0+1) .and. sign*xslope .lt. 0.) i0 = i0+1
        if (x0 .eq. xq(iq0)   .and. sign*xslope .gt. 0.) iq0 = iq0-1
        if (x0 .eq. xq(iq0+1) .and. sign*xslope .lt. 0.) iq0 = iq0+1
        if (xslope*sign .ge. 0.) xi = xq(iq0)
        if (xslope*sign .lt. 0.) xi = xq(iq0+1)
        if (x1 .eq. x0) then
          lenz0=10000
          do ii = 1,nnzQ
            if (((zq(ii)-z0)*(zq(ii)-z1) .le. 0) .and.
     &            ((abs(zq(ii)-z0)) .ge. (abs(zq(ii)-z1))) .and.
     &            (abs(zq(ii)-z0) .lt. lenz0)) then
              lenz0=abs(zq(ii)-z0)
              zi = zq(ii) ! between z1 and z0 but close to z0
            endif
          enddo
          leny0=10000
          do ii = 1,nnyQ
            if (((yq(ii)-y0)*(yq(ii)-y1) .le. 0) .and.
     &            ((abs(yq(ii)-y0)) .ge. (abs(yq(ii)-y1))) .and.
     &            (abs(yq(ii)-y0) .lt. leny0)) then
              leny0=abs(yq(ii)-y0)
              yi = yq(ii) ! between y1 and y0 but close to y0
            endif
          enddo
        else
          zi = xslope*xi + b
          yi = y0 + yxslope*(xi-x0)
        endif

c        do j = 1, nny   ! j0: node nearest 0 in y direction
c          if (y(j).lt.y0 .and. y(j+1).ge.y0 .and. yslope.ge.0.) j0 = j
c          if (y(j).le.y0 .and. y(j+1).gt.y0 .and. yslope.lt.0.) j0 = j
c        enddo
        jq0=0
        do j = 1, nnyQ
          if (yq(j).lt.y0 .and. yq(j+1).ge.y0 .and. yslope.ge.0.) jq0=j
          if (yq(j).le.y0 .and. yq(j+1).gt.y0 .and. yslope.lt.0.) jq0=j
        enddo
c        if (y0 .eq. y(j0)   .and. sign*yslope .gt. 0.) j0 = j0-1
c        if (y0 .eq. y(j0+1) .and. sign*yslope .lt. 0.) j0 = j0+1
        if (y0 .eq. yq(jq0)   .and. sign*yslope .gt. 0.) jq0 = jq0-1
        if (y0 .eq. yq(jq0+1) .and. sign*yslope .lt. 0.) jq0 = jq0+1
        if (yslope*sign .ge. 0.) yj = yq(jq0)
        if (yslope*sign .lt. 0.) yj = yq(jq0+1)
        if (y1 .eq. y0) then
          lenz0=10000
          do ii = 1,nnzQ
            if (((zq(ii)-z0)*(zq(ii)-z1) .le. 0) .and.
     &            ((abs(zq(ii)-z0)) .ge. (abs(zq(ii)-z1))) .and.
     &            (abs(zq(ii)-z0) .lt. lenz0)) then
              lenz0=abs(zq(ii)-z0)
              zj = zq(ii) ! between z1 and z0 but close to z0
            endif
          enddo
          lenx0=10000
          do ii = 1,nnyQ
            if (((xq(ii)-x0)*(xq(ii)-x1) .le. 0) .and.
     &            ((abs(xq(ii)-x0)) .ge. (abs(xq(ii)-x1))) .and.
     &            (abs(xq(ii)-x0) .lt. lenx0)) then
              leny0=abs(xq(ii)-x0)
              xj = xq(ii) ! between x1 and x0 but close to x0
            endif
          enddo
        else
          zj = yslope*yj + c
          xj = x0 + xyslope*(yj-y0)
        endif

c        do k = 1, nnz   ! k0: node nearest 0 in z direction
c          if (z(k) .lt. z0 .and. z(k+1) .ge. z0) k0 = k
c        enddo
        kq0 = 0
        do k = 1, nnzQ
          if (zq(k) .lt. z0 .and. zq(k+1) .ge. z0) kq0 = k
        enddo
c        if (z0 .eq. z(k0+1) .and. sign .lt. 0.) k0 = k0+1
        if (z0 .eq. zq(kq0+1) .and. sign .lt. 0.) kq0 = kq0+1
        if (kq0 .gt. 0) then
          if (sign .gt. 0.) zk = zq(kq0)
          if (sign .lt. 0.) zk = zq(kq0+1)
        else
          zk = z1
        endif
        if (z1 .eq. z0) then
          lenx0=10000
          do ii = 1,nnyQ
            if (((xq(ii)-x0)*(xq(ii)-x1) .le. 0) .and.
     &            ((abs(xq(ii)-x0)) .ge. (abs(xq(ii)-x1))) .and.
     &            (abs(xq(ii)-x0) .lt. lenx0)) then
              leny0=abs(xq(ii)-x0)
              xk = xq(ii) ! between x1 and x0 but close to x0
            endif
          enddo
          leny0=10000
          do ii = 1,nnyQ
            if (((yq(ii)-y0)*(yq(ii)-y1) .le. 0) .and.
     &            ((abs(yq(ii)-y0)) .ge. (abs(yq(ii)-y1))) .and.
     &            (abs(yq(ii)-y0) .lt. leny0)) then
              leny0=abs(yq(ii)-y0)
              yk = yq(ii) ! between y1 and y0 but close to y0
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
        if (xlength .lt. sslength) then   ! subsegment ends at x(i)
          icase = 2
          sslength = xlength
        endif
        if (ylength .lt. sslength) then    ! subsegment ends at y(j)
          icase = 3
          sslength = ylength
        endif
        if (z0 .le. zq(1)) then
          sslength = 0.
        endif
        if (zlength .lt. sslength) then    ! subsegment ends at z(k)
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
c        write(*,*) x1,y1,z1,xss,yss,zss
c        write(*,*) icase, sslength, xss, yss, zss
c        write(*,*) "start"
c        write(*,*) x0,y0,z0,xlength,ylength
c        write(*,*) xss,yss,zss,case1length,zlength

c     Find slowness at each Q node with trilinear interpolation
cccc          010--------110
cccc          /|         /|
cccc        /  |       /  |
cccc     000--------100   |
cccc      |    |     |    |
cccc      |   011--- | --111
cccc      |   /      |   /
cccc      | /        | /
cccc     001--------101
        iw = 0
        if ((id .eq. 1) .or. (id .eq. 2)) then          ! Qp or Qs
          lps = id
        elseif (id .eq. 3) then    ! Qp/Qs
          lps = 2
        elseif (id .eq. 4) then    ! Qk
          lps = 1
        endif
        do jy = jq0, jq0+1
          do kz = kq0, kq0+1
            do ix = iq0, iq0+1
              iw = iw + 1
              Qx=xq(ix)
              Qy=yq(jy)
              do i = 1, nnx
                if (x(i).lt.Qx .and. x(i+1).ge.Qx) iv = i
              enddo
              do j = 1, nny
                if (y(j).lt.Qy .and. y(j+1).ge.Qy) jv = j
              enddo
              if ((ix.eq.0).or.(jy.eq.0).or.(kz.eq.0)) then
                snQ(iw)=0.
              else
                Qz=zq(kz)
                do k = 1, nnz
                  if (z(k).lt.Qz .and. z(k+1).ge.Qz) kv = k
                enddo
                sn00=sn(lps,kv,jv,iv)+(Qx-x(iv))/(x(iv+1)-x(iv))*
     &                    (sn(lps,kv,jv,iv+1)-sn(lps,kv,jv,iv))
                sn01=sn(lps,kv+1,jv,iv)+(Qx-x(iv))/(x(iv+1)-x(iv))*
     &                    (sn(lps,kv+1,jv,iv+1)-sn(lps,kv+1,jv,iv))
                sn10=sn(lps,kv,jv+1,iv)+(Qx-x(iv))/(x(iv+1)-x(iv))*
     &                    (sn(lps,kv,jv+1,iv+1)-sn(lps,kv,jv+1,iv))
                sn11=sn(lps,kv+1,jv+1,iv)+(Qx-x(iv))/(x(iv+1)-x(iv))*
     &                    (sn(lps,kv+1,jv+1,iv+1)-sn(lps,kv+1,jv+1,iv))
                sn0=sn00+(Qy-y(jv))/(y(jv+1)-y(jv))*(sn10-sn00)
                sn1=sn01+(Qy-y(jv))/(y(jv+1)-y(jv))*(sn11-sn01)
                snQ(iw)=sn0+(Qz-z(kv))/(z(kv+1)-z(kv))*(sn1-sn0)
              endif
            enddo
          enddo
        enddo
        sn000 = snQ(1)
        sn100 = snQ(2)
        sn001 = snQ(3)
        sn101 = snQ(4)
        sn010 = snQ(5)
        sn110 = snQ(6)
        sn011 = snQ(7)
        sn111 = snQ(8)
        
        if (id .eq. 4) then    ! Qk
          do jy = jq0, jq0+1
            do kz = kq0, kq0+1
              do ix = iq0, iq0+1
                iw = iw + 1
                Qx=xq(ix)
                Qy=yq(jy)
                do i = 1, nnx
                  if (x(i).lt.Qx .and. x(i+1).ge.Qx) iv = i
                enddo
                do j = 1, nny
                  if (y(j).lt.Qy .and. y(j+1).ge.Qy) jv = j
                enddo
                if ((ix.eq.0).or.(jy.eq.0).or.(kz.eq.0)) then
                  snQ2(iw)=0.
                else
                  Qz=zq(kz)
                  do k = 1, nnz
                    if (z(k).lt.Qz .and. z(k+1).ge.Qz) kv = k
                  enddo
                  sn00=sn(2,kv,jv,iv)+(Qx-x(iv))/(x(iv+1)-x(iv))*
     &                (sn(2,kv,jv,iv+1)-sn(2,kv,jv,iv))
                  sn01=sn(2,kv+1,jv,iv)+(Qx-x(iv))/(x(iv+1)-x(iv))*
     &                (sn(2,kv+1,jv,iv+1)-sn(2,kv+1,jv,iv))
                  sn10=sn(2,kv,jv+1,iv)+(Qx-x(iv))/(x(iv+1)-x(iv))*
     &                (sn(2,kv,jv+1,iv+1)-sn(2,kv,jv+1,iv))
                  sn11=sn(2,kv+1,jv+1,iv)+(Qx-x(iv))/(x(iv+1)-x(iv))*
     &                (sn(2,kv+1,jv+1,iv+1)-sn(2,kv+1,jv+1,iv))
                  sn0=sn00+(Qy-y(jv))/(y(jv+1)-y(jv))*(sn10-sn00)
                  sn1=sn01+(Qy-y(jv))/(y(jv+1)-y(jv))*(sn11-sn01)
                  snQ2(iw)=sn0+(Qz-z(kv))/(z(kv+1)-z(kv))*(sn1-sn0)
                endif
              enddo
            enddo
          enddo
          sn000 = sn000/snQ2(1)
          sn100 = sn100/snQ2(2)
          sn001 = sn001/snQ2(3)
          sn101 = sn101/snQ2(4)
          sn010 = sn010/snQ2(5)
          sn110 = sn110/snQ2(6)
          sn011 = sn011/snQ2(7)
          sn111 = sn111/snQ2(8)        
        endif
c        write(*,*) sn000,iq0,jq0,kq0

c     Trapezoidal weighting: c - to begin of segment, d - to end of segment
        iw = 0
        slc = 0.
        sld = 0.
        do jy = j0, j0+1
          do kz = k0, k0+1
            do ix = i0, i0+1
              iw = iw + 1
              if ((ix.eq.0).or.(jy.eq.0).or.(kz.eq.0)) then
                wc(iw)=0.
                wd(iw)=0.
              else
                wc(iw) = sqrt(((x(ix)-x0)*(1.-z(kz)/re))**2 +
     &            ((y(jy)-y0)*(1.-z(kz)/re))**2 + (z(kz)-z0)**2)
                wd(iw) =  sqrt(((x(ix)-xss)*(1.-z(kz)/re))**2 +
     &            ((y(jy)-yss)*(1.-z(kz)/re))**2 + (z(kz)-zss)**2)
              endif
              slc = slc + wc(iw)                 ! sum of weighting lengths
              sld = sld + wd(iw)                 ! sum of weighting lengths
            enddo
          enddo
        enddo
        c000 = ((wc(2)+wc(3)+wc(4)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc
        c100 = ((wc(1)+wc(3)+wc(4)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc
        c001 = ((wc(1)+wc(2)+wc(4)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc
        c101 = ((wc(1)+wc(2)+wc(3)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc
        c010 = ((wc(1)+wc(2)+wc(3)+wc(4)+wc(6)+wc(7)+wc(8))/7.)/slc
        c110 = ((wc(1)+wc(2)+wc(3)+wc(4)+wc(5)+wc(7)+wc(8))/7.)/slc
        c011 = ((wc(1)+wc(2)+wc(3)+wc(4)+wc(5)+wc(6)+wc(8))/7.)/slc
        c111 = ((wc(1)+wc(2)+wc(3)+wc(4)+wc(5)+wc(6)+wc(7))/7.)/slc
        d000 = ((wd(2)+wd(3)+wd(4)+wd(5)+wd(6)+wd(7)+wd(8))/7.)/sld
        d100 = ((wd(1)+wd(3)+wd(4)+wd(5)+wd(6)+wd(7)+wd(8))/7.)/sld
        d001 = ((wd(1)+wd(2)+wd(4)+wd(5)+wd(6)+wd(7)+wd(8))/7.)/sld
        d101 = ((wd(1)+wd(2)+wd(3)+wd(5)+wd(6)+wd(7)+wd(8))/7.)/sld
        d010 = ((wd(1)+wd(2)+wd(3)+wd(4)+wd(6)+wd(7)+wd(8))/7.)/sld
        d110 = ((wd(1)+wd(2)+wd(3)+wd(4)+wd(5)+wd(7)+wd(8))/7.)/sld
        d011 = ((wd(1)+wd(2)+wd(3)+wd(4)+wd(5)+wd(6)+wd(8))/7.)/sld
        d111 = ((wd(1)+wd(2)+wd(3)+wd(4)+wd(5)+wd(6)+wd(7))/7.)/sld
        if (abs(1.-(c000+c001+c100+c101+c010+c011+c110+c111)).gt..0001)
     &          write(*,*) "WARNING: cumulative wc = ",
     &          c000+c001+c100+c101+c010+c011+c110+c111

        if (abs(1.-(d000+d001+d100+d101+d010+d011+d110+d111)).gt..0001)
     &          write(*,*) "WARNING: cumulative wd = ",
     &          d000+d001+d100+d101+d010+d011+d110+d111
        xlen = xss-x0
        ylen = yss-y0
        zlen = zss-z0
        if (ipolar .eq. 1) then     ! x,y distances decr w/ depth
          xlen = xlen*(1.-(zss+z0)/(2.*re))
          ylen = ylen*(1.-(zss+z0)/(2.*re))
        endif
        sslen = sqrt((xlen**2) + (ylen**2) + (zlen**2))
        wtc = 0.5
c        write(*,*) "Length of the subsegment is ",sslen,icase
c        write(*,*) "wtc*(d110+c110) = ", wtc*(d110+c110)
c        alllen=alllen+sslen

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Qp: grow(element)=grow(element)+.5*(wgt)*distn*Pslo             C
C      Qs: grow(element)=grow(element)+.5*(wgt)*distn*Sslo             C
C   Qp/Qs: grow(element)=grow(element)+.5*(wgt)*distn*Sslo/Qp          C
C      Qk: grow(element)=grow(element)+.5*(wgt)*(1-3/4*(Pslo/Sslo)**2) C
C           wgt = some "trapezoidal weighting" function                C
C           distn = pathlength in that block                           C
C           ?slo = P or S slowness in that block                       C
C           element = pt-th block, counting in the order of z, x, y    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        if ((id .eq. 1) .or. (id .eq. 2)) then          ! Qp or Qs
          grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0) =
     &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0)+wtc*(d000+c000)*
     &        sslen*sn000
          grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1) =
     &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1)+wtc*(d001+c001)
     &        *sslen*sn001
          grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0) =
     &        grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0) + wtc*(d010+c010)*
     &        sslen*sn010
          grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1) =
     &        grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1) + wtc*(d011+c011)*
     &        sslen*sn011
          grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0) =
     &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0) + wtc*(d100+c100)*
     &        sslen*sn100
          grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0+1) =
     &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0+1) + wtc*(d101+c101)*
     &        sslen*sn101
          grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0) =
     &        grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0) + wtc*(d110+c110)*
     &        sslen*sn110
          grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0+1) =
     &        grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0+1) + wtc*(d111+c111)*
     &        sslen*sn111
        else if (id .eq. 3) then    ! Qp/Qs
          grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0) =
     &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0)+wtc*(d000+c000)*
     &        sslen*sn000*attin(kq0,jq0,iq0)
          grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1) =
     &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1)+wtc*(d001+c001)
     &        *sslen*sn001*attin(kq0+1,jq0,iq0)
          grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0) =
     &        grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0) + wtc*(d010+c010)*
     &        sslen*sn010*attin(kq0,jq0,iq0+1)
          grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1) =
     &        grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1) + wtc*(d011+c011)*
     &        sslen*sn011*attin(kq0+1,jq0,iq0+1)
          grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0) =
     &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0) + wtc*(d100+c100)*
     &        sslen*sn100*attin(kq0,jq0+1,iq0)
          grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0+1) =
     &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0+1) + wtc*(d101+c101)*
     &        sslen*sn101*attin(kq0+1,jq0+1,iq0)
          grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0) =
     &        grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0) + wtc*(d110+c110)*
     &        sslen*sn110*attin(kq0,jq0+1,iq0+1)
          grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0+1) =
     &        grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0+1) + wtc*(d111+c111)*
     &        sslen*sn111*attin(kq0+1,jq0+1,iq0+1)
C       else if (id .eq. 4) then    ! Qp/Qs: sn = Vs/Vp
C         grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0) =
C    &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0)+wtc*(d000+c000)*
C    &        (1-3./4.*sn000*sn000)
C         grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1) =
C    &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1)+wtc*(d001+c001)
C    &        *(1-3./4.*sn001*sn001)
C         grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0) =
C    &        grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0) + wtc*(d010+c010)*
C    &        (1-3./4.*sn010*sn010)
C         grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1) =
C    &        grow(jq0*nnzQ+nnyQ*nnzQ*(iq0-1)+kq0+1) + wtc*(d011+c011)*
C    &        (1-3./4.*sn011*sn011)
C         grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0) =
C    &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0) + wtc*(d100+c100)*
C    &        (1-3./4.*sn100*sn100)
C         grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0+1) =
C    &        grow((jq0-1)*nnzQ+nnyQ*nnzQ*iq0+kq0+1) + wtc*(d101+c101)*
C    &        (1-3./4.*sn101*sn101)
C         grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0) =
C    &        grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0) + wtc*(d110+c110)*
C    &        (1-3./4.*sn110*sn110)
C         grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0+1) =
C    &        grow(jq0*nnzQ+nnyQ*nnzQ*iq0+kq0+1) + wtc*(d111+c111)*
C    &        (1-3./4.*sn111*sn111)
        endif

        x0 = xss
        y0 = yss
        z0 = zss
c        write(*,*) x0, y0, z0
        if (z0 .eq. zarc(ipt1)) then
          pt = pt+1
        elseif (icase .eq. 1) then      ! Move to next subsegment for case 1
          pt = pt + 1
          x1 = xarc(pt)
          y1 = yarc(pt)
          z1 = zarc(pt)
        endif
c        write(*,*) x0,y0,z0,pt,icase
c        write(*,*) "pt, npts+1: ", pt, npts+1
      enddo       ! end of loop over each raypath node (pt -> ipt1)
c      write(*,*) "alllen2 = :",alllen
237   continue
c      write(*,400) ifskip,npts+2,pt-1,ipt,ipt1,
c     &                  zarc(pt-1),zarc(0),zarc(ipt),zarc(ipt1)
c      write(*,*) "Finish adding row of G matrix"
200   format(3(f10.4,1x))
300   format(a15,3(f10.4,1x))
400   format(5(i4,1x),4(f10.4,1x))
      END
