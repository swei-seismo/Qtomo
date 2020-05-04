ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Subroutine that does pseudo-bending in 3D.  Returns tt,   c
c   raypath (xarc, yarc, zarc, npts).                         c
c   Bending starts withstraight lines from source to receiver c
c     npts = # interior pts that can bend.                    c
c     kpt = 0 => hypocenter;                                  c
c     kpt = npts+1 => moho piercing pt;                       c
c     kpt = npts+2 => station location.                       c
c   Add an additional output by S. Wei, July 2018:            c
c     LttS = Vs^3/Vp^2/len
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE bendray(ilps,xe0,ye0,ze0,xs0,ys0,zs0,zmoho,maxpts,
     &     xarc,yarc,zarc,npts,tt,LttS)
      parameter (nxmax = 45, nymax = 45, nzmax = 30)
      integer ilps, maxpts, pt          ! max # internal points in raypath
      real x(nxmax), z(nzmax), y(nymax), xq(nxmax), zq(nzmax), yq(nymax)  ! node locations
      real*8 sn(2,nzmax,nymax,nxmax)    ! velocity, slowness
      real*8 xarc(0:maxpts+2), yarc(0:maxpts+2), zarc(0:maxpts+2)
      real*8 dwx(4), dwy(4), dwz(4), wc(8), wd(8), slc            ! weighting coefficients
      real*8 vl(3)                              ! local velocities for pseudo-bending
      real*8 thetax, thetay, phitmp, hnx, hny, hnz, hlen,alllen
      real*8 a, b, c, a1, b1, c1, alen, rc
      real*8 dvdx, dvdy, dvdz, gradvn, f, f1
      real LttS
      character*2 wtd
      common /a/ x, y, z, xq, yq, zq
      common /b/ ipolar,nnx,nny,nnz,re,pi,r2d,nnxQ,nnyQ,nnzQ,nnzQ1D
      common /c/ sn

C     if (ilps .eq. 3) then    ! Qk needs LttS
C       lps = 1
C       ifqk = 1
C     else
        lps = ilps
        ifqk = 0
C     endif
C     write(*,*) lps, ifqk
      vpcrust = 6.101               ! crustal velocity: P
      vn = 7.78               ! Pn velocity
      if (lps .eq. 2) vpcrust = 3.5 ! crustal velocity: S
      if (lps .eq. 2) vn = 4.28     ! Sn velocity
      zcrust = zmoho - zs0

c starting raypath
c      write(*,*) "start ray tracing"
      if (ipolar .eq. 1) then       ! convert from curved earth to cartesian coord
        thetax = xe0/re             ! angle from x = 0
        thetay = ye0/re             ! angle from y = 0
        xe = (re - ze0)*sin(thetax) 
        ye = (re - ze0)*sin(thetay) 
        ze = (re - (re - ze0)*cos(thetax)*cos(thetay))
        thetax = xs0/re
        thetay = ys0/re
        xs = (re - zmoho)*sin(thetax)
        ys = (re - zmoho)*sin(thetay)
        zs = (re - (re - zmoho)*cos(thetax)*cos(thetay))
      endif
c      write(*,*) "curved earth: ",xs0,ys0,zs0,xs,ys,zs
c      write(*,*) "curved earth: ",xe0,ye0,ze0,xe,ye,ze

      npts = 3                ! Initiating the ray as a line
      xarc(0) = xe0
      yarc(0) = ye0
      zarc(0) = ze0
c      if (iall .ne. 0) then
c        write (21,1) ">", npts
c        write (21,*) xe0, ye0, ze0
c      endif

      do i = 1, npts
        xtemp =  xe+real(i)*(xs-xe)/real(npts+1)
        ytemp =  ye+real(i)*(ys-ye)/real(npts+1)
        ztemp =  ze+real(i)*(zs-ze)/real(npts+1)
        d2 = (xe0-xs0)**2 + (ye0-ys0)**2
c        if (ze0 .lt. 50 .and. d2 .gt. 40000.) then
        if (ze0 .lt. 50) then
          ztemp = ztemp + d2/10000.
        endif
        if ( .true. ) then          ! put back into curved coord
          dtemp = (xtemp*xtemp + ytemp*ytemp)**0.5
          thetax = atan(xtemp/(re-ztemp))
          thetay = atan(ytemp/(re-ztemp))
          thetaz = atan(dtemp/(re-ztemp))
          ztemp = re - (dtemp/sin(thetaz))
          xtemp = re*thetax
          ytemp = re*thetay
        endif
        xarc(i) = xtemp
        yarc(i) = ytemp
        zarc(i) = ztemp
c        if (iall.ne.0) write (21,*) xtemp, ytemp, ztemp
      enddo
      xarc(npts+1) = xs0
      yarc(npts+1) = ys0
      zarc(npts+1) = zmoho
      xarc(npts+2) = xs0
      yarc(npts+2) = ys0
      zarc(npts+2) = zs0
c      if (iall .ne. 0) then
c        write (21,*) xs0, ys0, zmoho
c        write (21,*) xs0, ys0, zs0
c      endif

1     format (a1,1x,i4)
2     format (a1)

      do pt = 0, npts+1
        if (xarc(pt) .le. x(1)) xarc(pt) = x(1) + .005  ! keep arc in model space
        if (xarc(pt) .ge. x(nnx)) xarc(pt) = x(nnx) - .005
c        write(*,*) "pt = ",pt,yarc(pt)
        if (yarc(pt) .le. y(1)) yarc(pt) = y(1) + .005  ! keep arc in model space
        if (yarc(pt) .ge. y(nny)) yarc(pt) = y(nny) - .005
c        write(*,*) "pt = ",pt,yarc(pt)
        if (zarc(pt) .ge. z(nnz)) zarc(pt) = z(nnz) - .005
      enddo

cccc  -- pseudo-bend raypath
      tt0 = 1.e5
      idbl = 0
      epsilon = .001                    ! threshold. Um & thurber: .001
100   ttchng = 10.                      ! tt improvement for subsequent raypaths
      f = 1.01                            ! "enhancement" factor
      do while (abs(ttchng) .gt. epsilon)
cccccccc find xmoho, ymoho
c        write(*,*) "start one loop of psudo-bending"
        if (zmoho .eq. zs0) then
          xmoho = xs0
          ymoho = ys0
          hdist = 0.
          d = 0.
        else
          xm = xarc(npts)
          ym = yarc(npts)
          hdist = sqrt((xm-xs0)**2 + (ym-ys0)**2)
          zl2 = zarc(npts)-zmoho
          d = piercept(hdist,zcrust,zl2,vpcrust,vn)   ! moho piercing pt
          tmpcosphi = abs(xm-xs0)/sqrt((xm-xs0)**2 + (ym-ys0)**2)
          if (tmpcosphi .ge. 1) tmpcosphi = 1
c          write(*,*) 'testtmp = ',tmpcosphi
          phitmp = acos(tmpcosphi)
c          write(*,*) 'phitmp = ',phitmp
          xsign = +1.               ! recheck as sign can flip when phi ~ pi/2 
          if (xm .gt. xs0) xsign = -1.
          ysign = +1.               ! recheck as sign can flip when phi ~ 0 
          if (ym .gt. ys0) ysign = -1.
          xmoho = xm + d*cos(phitmp)*xsign
          ymoho = ym + d*sin(phitmp)*ysign
        endif
c        write(*,*) 'Moho:',xmoho,ymoho,hdist,d
        if (xmoho .gt. x(nnx)) then       ! Only nec for deep, nearby events at Niue
          write (*,*) "WARNING: xmoho adjusted. Was ", xmoho
c         xmoho = x(nnx)
        endif 
        xarc(npts+1) = xmoho
        yarc(npts+1) = ymoho
        tcrust = sqrt((hdist-d)**2 + zcrust**2)/vpcrust
c        write(*,*) 'tcrust',tcrust

ccccccccc find other points
        do ipt = 1, npts              ! pseudo bending loop, need to add do i,j,k - nnx,nny,nnz
          kpt = 1 + (ipt-1)/2         ! work from ends inward
          if (mod(ipt,2) .eq. 0) kpt = npts + 1 - (ipt/2)
          xmid = (xarc(kpt+1) + xarc(kpt-1))/2.
          ymid = (yarc(kpt+1) + yarc(kpt-1))/2.
          zmid = (zarc(kpt+1) + zarc(kpt-1))/2.
          do i = 1, nnx
            if (x(i).lt.xarc(kpt-1).and.x(i+1).ge.xarc(kpt-1)) ik1 = i
            if (x(i).lt.xarc(kpt+1).and.x(i+1).ge.xarc(kpt+1)) ik2 = i
            if (x(i) .lt. xmid .and. x(i+1) .ge. xmid) ik3 = i
          enddo
c          write(*,*) xarc(kpt-1),xarc(kpt+1),xmid,kpt
          do j = 1, nnz
            if (z(j).lt.zarc(kpt-1).and.z(j+1).ge.zarc(kpt-1)) jk1 = j
            if (z(j).lt.zarc(kpt+1).and.z(j+1).ge.zarc(kpt+1)) jk2 = j
            if (z(j) .lt. zmid .and. z(j+1) .ge. zmid) jk3 = j
          enddo
c          write(*,*) zarc(kpt-1),zarc(kpt+1),zmid
          do k = 1, nny
            if (y(k).lt.yarc(kpt-1).and.y(k+1).ge.yarc(kpt-1)) kk1 = k
            if (y(k).lt.yarc(kpt+1).and.y(k+1).ge.yarc(kpt+1)) kk2 = k
            if (y(k) .lt. ymid .and. y(k+1) .ge. ymid) kk3 = k
          enddo
c          write(*,*) yarc(kpt-1),yarc(kpt+1),ymid
c          write(*,*) ik1,ik2,ik3,jk1,jk2,jk3,kk1,kk2,kk3
          do l = 1,3
            if (l .eq. 1) then
              i0 = ik1
              j0 = jk1
              k0 = kk1
              xl = xarc(kpt-1)
              yl = yarc(kpt-1)
              zl = zarc(kpt-1)
            else if (l.eq. 2) then
              i0 = ik2
              j0 = jk2
              k0 = kk2
              xl = xarc(kpt+1)
              yl = yarc(kpt+1)
              zl = zarc(kpt+1)
            else
              i0 = ik3
              j0 = jk3
              k0 = kk3
              xl = xmid
              yl = ymid
              zl = zmid
            endif
                                  ! find velocity at Xmid: trapezoidal weighting
            iw = 0                  ! need to fix these to 8 factors
            slc = 0.                ! only matters for the travel times
            do ky = k0, k0+1
              do jz = j0, j0+1
                do ix = i0, i0+1
                  iw = iw + 1
                  wc(iw) =  sqrt(((x(ix)-xl)*(1.-z(jz)/re))**2 +
     &                ((y(ky)-yl)*(1.-z(jz)/re))**2 + (z(jz)-zl)**2)
                  slc = slc + wc(iw)                 ! sum of weighting lengths
                enddo
              enddo
            enddo
c            write(*,*) "slc = ",slc
            c000=((wc(2)+wc(3)+wc(4)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc ! weighting factors
            c001=((wc(1)+wc(3)+wc(4)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc
            c100=((wc(1)+wc(2)+wc(4)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc
            c101=((wc(1)+wc(2)+wc(3)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc
            c010=((wc(1)+wc(2)+wc(3)+wc(4)+wc(6)+wc(7)+wc(8))/7.)/slc
            c011=((wc(1)+wc(2)+wc(3)+wc(4)+wc(5)+wc(7)+wc(8))/7.)/slc
            c110=((wc(1)+wc(2)+wc(3)+wc(4)+wc(5)+wc(6)+wc(8))/7.)/slc
            c111=((wc(1)+wc(2)+wc(3)+wc(4)+wc(5)+wc(6)+wc(7))/7.)/slc
  
            if (abs(1.-(c000+c001+c100+c101+c010+c011+c110+c111)) 
     &                 .gt. .0001) then
              write(*,*) "WARNING: cumulative wc = ",
     &                     c000+c001+c100+c101+c010+c011+c110+c111
            endif

            stmp = c000*sn(lps,j0,k0,i0) + c001*sn(lps,j0,k0,i0+1) +
     &         c100*sn(lps,j0+1,k0,i0) + c101*sn(lps,j0+1,k0,i0+1) +
     &         c010*sn(lps,j0,k0+1,i0) + c011*sn(lps,j0,k0+1,i0+1) +
     &         c110*sn(lps,j0+1,k0+1,i0) + c111*sn(lps,j0+1,k0+1,i0+1)
            vl(l) = 1./stmp  ! vl(l): average velocity in a small region near l
          enddo           ! enddo l

          vmid = vl(3)
          c = (1./vl(1) + 1./vl(2))/2.    ! c in Um & Thurber Eqa 6
  
          v111 = 1./sn(lps,j0,k0,i0)          ! find velocity gradient at Xmid
          v112 = 1./sn(lps,j0,k0,i0+1)
          v211 = 1./sn(lps,j0+1,k0,i0)
          v212 = 1./sn(lps,j0+1,k0,i0+1)
          v121 = 1./sn(lps,j0,k0+1,i0)
          v122 = 1./sn(lps,j0,k0+1,i0+1)
          v221 = 1./sn(lps,j0+1,k0+1,i0)
          v222 = 1./sn(lps,j0+1,k0+1,i0+1)
  
          dx = x(i0+1) - x(i0) ! need a dy and new indices
          dz = z(j0+1) - z(j0)
          dy = y(k0+1) - y(k0)
  
          distx1 = sqrt((ymid - y(k0))**2 + (zmid - z(j0))**2)
          distx2 = sqrt((ymid - y(k0+1))**2 + (zmid - z(j0))**2)
          distx3 = sqrt((ymid - y(k0))**2 + (zmid - z(j0+1))**2)
          distx4 = sqrt((ymid - y(k0+1))**2 + (zmid - z(j0+1))**2)
          distxsum = distx1+distx2+distx3+distx4
  
          dwx(1) = ((distx2 + distx3 + distx4)/3.)/distxsum
          dwx(2) = ((distx1 + distx3 + distx4)/3.)/distxsum
          dwx(3) = ((distx1 + distx2 + distx4)/3.)/distxsum
          dwx(4) = ((distx1 + distx2 + distx3)/3.)/distxsum
  
          dsx = dwx(1)+dwx(2)+dwx(3)+dwx(4)
  
          disty1 = sqrt((xmid - x(i0))**2 + (zmid - z(j0))**2)
          disty2 = sqrt((xmid - x(i0+1))**2 + (zmid - z(j0))**2)
          disty3 = sqrt((xmid - x(i0))**2 + (zmid - z(j0+1))**2)
          disty4 = sqrt((xmid - x(i0+1))**2 + (zmid - z(j0+1))**2)
          distysum = disty1+disty2+disty3+disty4
  
          dwy(1) = ((disty2 + disty3 + disty4)/3.)/distysum
          dwy(2) = ((disty1 + disty3 + disty4)/3.)/distysum
          dwy(3) = ((disty1 + disty2 + disty4)/3.)/distysum
          dwy(4) = ((disty1 + disty2 + disty3)/3.)/distysum
  
          dsy = dwy(1)+dwy(2)+dwy(3)+dwy(4)
  
          distz1 = sqrt((xmid - x(i0))**2 + (ymid - y(k0))**2)
          distz2 = sqrt((xmid - x(i0+1))**2 + (ymid - y(k0))**2)
          distz3 = sqrt((xmid - x(i0))**2 + (ymid - y(k0+1))**2)
          distz4 = sqrt((xmid - x(i0+1))**2 + (ymid - y(k0+1))**2)
          distzsum = distz1+distz2+distz3+distz4
  
          dwz(1) = ((distz2 + distz3 + distz4)/3.)/distzsum
          dwz(2) = ((distz1 + distz3 + distz4)/3.)/distzsum
          dwz(3) = ((distz1 + distz2 + distz4)/3.)/distzsum
          dwz(4) = ((distz1 + distz2 + distz3)/3.)/distzsum
  
          dsz = dwz(1)+dwz(2)+dwz(3)+dwz(4)
c       write(*,*) 'dsx,dsy,dsz',dsx,dsy,dsz

c      velocity gradient at mid-point
          dvdx = ((v112-v111)*dwx(1) + (v212-v211)*dwx(3) +
     &           (v122-v121)*dwx(2) + (v222-v221)*dwx(4))/dx

          dvdy = ((v121-v111)*dwy(1) + (v122-v112)*dwy(2) +
     &           (v221-v211)*dwy(3) + (v222-v212)*dwy(4))/dy

          dvdz = ((v211-v111)*dwz(1) + (v221-v121)*dwz(3) +
     &           (v212-v112)*dwz(2) + (v222-v122)*dwz(4))/dz
c     write(*,*) dvdx,dvdy,dvdz

          ihead = 1
          if (xe0-xs0.gt.ze0*2.5) ihead = 0
          if (xs0 .ge. x(nnx)) ihead = 0
          if (dvdz.lt.0 .and. ihead.eq.0 .and. z(j0).le.100.) then  ! no head wave
            dvdz = 0.0006
          endif
          a1 = xarc(kpt+1) - xarc(kpt-1) ! length of subsegment in x direc
          b1 = yarc(kpt+1) - yarc(kpt-1) ! "" "" y direc
          c1 = zarc(kpt+1) - zarc(kpt-1) ! in z dirc
          d1 = (a1*dvdx + b1*dvdy + c1*dvdz)/(a1*a1 + b1*b1 + c1*c1)
cccccccccc Um & Thurber Eqa 4
          hnx = dvdx - a1*d1    ! direction to bend (n')
          hny = dvdy - b1*d1
          hnz = dvdz - c1*d1
          hlen = sqrt(hnx*hnx + hny*hny + hnz*hnz)
c                 write (*,*) hnx, hny, hnz, hlen, a1, b1, c1, d1
          hnx = hnx/hlen        ! normalize n'
          hny = hny/hlen
          hnz = hnz/hlen
          alen = sqrt((xarc(kpt+1)-xmid)**2 + (yarc(kpt+1)-ymid)**2 +
     &                (zarc(kpt+1)-zmid)**2)          ! L in Um & Thurber Eqa 6
          a = (c*vmid + 1.)/(4.*c)
          b = alen*alen/(2.*c*vmid)
          gradvn = dvdx*hnx + dvdy*hny + dvdz*hnz                   ! gradVmid dot n
c                write (*,*) dvdx, hnx, dvdy, hny, dvdz, hnz, gradvn
          if (gradvn .eq. 0.) gradvn = 0.00001
cccccccccc Um & Thurber Eqa 6
          rc=-a/gradvn+sqrt(a*a/(gradvn*gradvn)+alen*alen/(2.*c*vmid))
          f1 = f

c     Don't change any points within the crust
          icount = 0
76        ztemp = f1*(zmid + rc*hnz - zarc(kpt)) + zarc(kpt)
          if (ztemp.lt.zmoho+0.4 .and. icount.le.70) then
            f1 = 0.8*f1
            icount = icount + 1
            go to 76
          endif

c     Don't pertubate points out of the model space
77        xarctmp = f1*(xmid + rc*hnx - xarc(kpt)) + xarc(kpt)
          yarctmp = f1*(ymid + rc*hny - yarc(kpt)) + yarc(kpt)
          zarctmp = f1*(zmid + rc*hnz - zarc(kpt)) + zarc(kpt)
          if (zarctmp.gt.z(nnz) .or. (zarctmp.lt.ze0.and.zarctmp.lt.zs0)
     &         .or. xarctmp.gt.x(nnx) .or. xarctmp.lt.x(1)) then
c         write (*,*) "Ray bent out of model space: ", zarctmp
            f1 = 0.8*f1
            go to 77
          endif
          xarc(kpt) = xarctmp
          yarc(kpt) = yarctmp
          zarc(kpt) = zarctmp
          if (zarc(kpt).lt.z(1) .or. zarc(kpt).gt.z(nnz)) then
            write (*,*) "ERROR: raypath bent out of model space ",
     &                        zarc(kpt)
            if (zarc(kpt) .gt. z(nnz)) zarc(kpt) = z(nnz) - 2.
            if (zarc(kpt) .lt. z(1)) zarc(kpt) = zmoho+0.4
c         stop
          endif
        enddo ! ipt

c        if (iall .ne. 0) then
c          write (21,1) ">", npts
          xarc(npts+2) = xs0        ! complete raypath
          yarc(npts+2) = ys0
          zarc(npts+2) = zs0
c          do i = 0, npts+2
c            write (21,*) xarc(i), yarc(i), zarc(i)
c          enddo
c        endif

cccc  -- find travel-time for each arc segment -- cccc
        pt = 1
        tt = 0.
        LttS = 0.
c        alllen = 0.
        x0 = xarc(0)            ! segment ends
        y0 = yarc(0)
        z0 = zarc(0)
        x1 = xarc(1)
        y1 = yarc(1)
        z1 = zarc(1)
c        write(*,*) 'segment endpts:', y0,y1

        jcount=0
        do while (pt .le. npts+1)            ! find subsegments endpts
          sign = +1.            ! ray moving upwards
          if (z1 .gt. z0) sign = -1.
c          write(*,*) 'x1,x0,y1,y0',x1,x0,y1,y0,(x1-x0),(y1-y0),icase
          if (x1 .eq. x0) then
            xslope = 100000.
            yxslope = 100000.
          else
            xslope = (z1-z0)/(x1-x0)    ! positive means up & left(west) or down & right(east)
            yxslope = (y1-y0)/(x1-x0)
          endif
          if (y1 .eq. y0) then
            yslope = 100000.
            xyslope = 100000.
          else
            yslope = (z1-z0)/(y1-y0)      ! positive means up & front(south) or down & back(north)
            xyslope = (x1-x0)/(y1-y0)
          endif
ccccc     sign*xslope => ray moving right-left (+left, -right) 
ccccc     sign*yslope => ray moving front-back (+front, -back)
c          write(*,*) "xslope,yslope:", xslope,yslope
          b = z1 - xslope*x1
          c = z1 - yslope*y1
          do j = 1, nnz
            if (z(j) .lt. z0 .and. z(j+1) .ge. z0) j0 = j
          enddo
          if (z0 .eq. z(j0+1) .and. sign .lt. 0.) j0 = j0+1
          if (sign .gt. 0.) zj = z(j0)
          if (sign .lt. 0.) zj = z(j0+1)
          if (xslope .eq. 0) then
            xj = (zj - b)*100000
          else
            xj = (zj - b)/xslope
          endif
          yj = y0 + yxslope*(xj-x0)
          do i = 1, nnx
            if (x(i).lt.x0 .and. x(i+1).ge.x0 .and. xslope.ge.0) i0 = i
            if (x(i).le.x0 .and. x(i+1).gt.x0 .and. xslope.lt.0) i0 = i
          enddo
          if (x0 .eq. x(i0)   .and. sign*xslope .gt. 0.) i0 = i0-1
          if (x0 .eq. x(i0+1) .and. sign*xslope .lt. 0.) i0 = i0+1
          do k = 1, nny
            if (y(k).lt.y0 .and. y(k+1).ge.y0 .and. yslope.ge.0) k0 = k
            if (y(k).le.y0 .and. y(k+1).gt.y0 .and. yslope.lt.0) k0 = k
          enddo
          if (y0 .eq. y(k0)   .and. sign*yslope .gt. 0.) k0 = k0-1
          if (y0 .eq. y(k0+1) .and. sign*yslope .lt. 0.) k0 = k0+1
c          write(*,*) 'k0:', k0

          if (xslope*sign .ge. 0.) xi = x(i0)
          if (xslope*sign .lt. 0.) xi = x(i0+1)
          zi = xslope*xi + b
          yi = y0 + yxslope*(xi-x0)
          if (yslope*sign .ge. 0.) yk = y(k0)
          if (yslope*sign .lt. 0.) yk = y(k0+1)
          zk = yslope*yk + c
c          write(*,*) 'yk, y0, y1 ',yk, y0, y1
          xk = x0 + xyslope*(yk-y0)

          case1length = (x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2
          xlength = (xi-x0)**2 + (yi-y0)**2 + (zi-z0)**2
          ylength = (xk-x0)**2 + (yk-y0)**2 + (zk-z0)**2
          zlength = (xj-x0)**2 + (yj-y0)**2 + (zj-z0)**2

          icase = 1                   ! arc segment ends at node
          sslength = case1length
          if (xlength .le. sslength) then
            icase = 2                 ! arc segment ends at x(i)
            sslength = xlength
          endif
          if (ylength .le. sslength) then
            icase = 3                 ! arc segment ends at y(k)
            sslength = ylength
          endif
          if (zlength .le. sslength) then
            icase = 4                 ! arc segment ends at z(j)
            sslength = zlength
          endif

c       open(33,file='icase')
c           write(33,*) pt, icase

          if (icase .eq. 1) then            ! fix subsegment end, add icase 4
            xss = x1
            yss = y1
            zss = z1
          else if (icase .eq. 2) then
            xss = xi
            yss = yi
            zss = zi
          else if (icase .eq. 4) then
            xss = xj
            yss = yj
            zss = zj
          else if (icase .eq. 3) then
            xss = xk
            yss = yk
            zss = zk
          else
            write (*,*) "icase out of bounds, (", icase, ")"
            write (*,*) 'x1,y1,z1,xi,yi,zi,xj,yj,zj,xk,yk,zk'
            write (*,*) x1,y1,z1,xi,yi,zi,xj,yj,zj,xk,yk,zk
            goto 666
c            stop
          endif
c          write(*,*) icase, sslength, xss, yss, zss
C        trapezoidal weighting
          iw = 0
          slc = 0.
          sld = 0.
          do ky = k0, k0+1
            do jz = j0, j0+1
              do ix = i0, i0+1
                iw = iw + 1
                wc(iw) =  sqrt(((x(ix)-x0)*(1.-z(jz)/re))**2 +
     &              ((y(ky)-y0)*(1.-z(jz)/re))**2 + (z(jz)-z0)**2)
                wd(iw) =  sqrt(((x(ix)-xss)*(1.-z(jz)/re))**2 +
     &              ((y(ky)-yss)*(1.-z(jz)/re))**2 + (z(jz)-zss)**2)
                slc = slc + wc(iw)                   ! sum of weighting lengths
                sld = sld + wd(iw)                   ! sum of weighting lengths
              enddo
            enddo
          enddo

          c000 = ((wc(2)+wc(3)+wc(4)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc  ! weighting factors
          c001 = ((wc(1)+wc(3)+wc(4)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc
          c100 = ((wc(1)+wc(2)+wc(4)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc
          c101 = ((wc(1)+wc(2)+wc(3)+wc(5)+wc(6)+wc(7)+wc(8))/7.)/slc
          c010 = ((wc(1)+wc(2)+wc(3)+wc(4)+wc(6)+wc(7)+wc(8))/7.)/slc
          c011 = ((wc(1)+wc(2)+wc(3)+wc(4)+wc(5)+wc(7)+wc(8))/7.)/slc
          c110 = ((wc(1)+wc(2)+wc(3)+wc(4)+wc(5)+wc(6)+wc(8))/7.)/slc
          c111 = ((wc(1)+wc(2)+wc(3)+wc(4)+wc(5)+wc(6)+wc(7))/7.)/slc

          d000 = ((wd(2)+wd(3)+wd(4)+wd(5)+wd(6)+wd(7)+wd(8))/7.)/sld  ! weighting factors
          d001 = ((wd(1)+wd(3)+wd(4)+wd(5)+wd(6)+wd(7)+wd(8))/7.)/sld
          d100 = ((wd(1)+wd(2)+wd(4)+wd(5)+wd(6)+wd(7)+wd(8))/7.)/sld
          d101 = ((wd(1)+wd(2)+wd(3)+wd(5)+wd(6)+wd(7)+wd(8))/7.)/sld
          d010 = ((wd(1)+wd(2)+wd(3)+wd(4)+wd(6)+wd(7)+wd(8))/7.)/sld
          d011 = ((wd(1)+wd(2)+wd(3)+wd(4)+wd(5)+wd(7)+wd(8))/7.)/sld
          d110 = ((wd(1)+wd(2)+wd(3)+wd(4)+wd(5)+wd(6)+wd(8))/7.)/sld
          d111 = ((wd(1)+wd(2)+wd(3)+wd(4)+wd(5)+wd(6)+wd(7))/7.)/sld

          if (abs(1.-(c000+c001+c100+c101+c010+c011+c110+c111))
     &            .gt. .0001) then
            write(*,*) "WARNING: cumulative wc = ",
     &                   c000+c001+c100+c101+c010+c011+c110+c111
          endif

          if (abs(1.-(d000+d001+d100+d101+d010+d011+d110+d111))
     &            .gt. .0001) then
            write(*,*) "WARNING: cumulative wd = ",
     &                   d000+d001+d100+d101+d010+d011+d110+d111
          endif

          xlen = xss-x0
          ylen = yss-y0
          zlen = zss-z0
          if (ipolar .eq. 1) then     ! x,y distances decr w/ depth
            xlen = xlen*(1.-(zss+z0)/(2.*re))
            ylen = ylen*(1.-(zss+z0)/(2.*re))
          endif
          sslen = sqrt((xlen**2) + (ylen**2) + (zlen**2))
          sn0  = c000*sn(lps,j0,k0,i0)   + c001*sn(lps,j0,k0,i0+1)  +
     &           c100*sn(lps,j0+1,k0,i0) + c101*sn(lps,j0+1,k0,i0+1) +
     &           c010*sn(lps,j0,k0+1,i0) + c011*sn(lps,j0,k0+1,i0+1) +
     &       c110*sn(lps,j0+1,k0+1,i0) + c111*sn(lps,j0+1,k0+1,i0+1)
          snss  = d000*sn(lps,j0,k0,i0)   + d001*sn(lps,j0,k0,i0+1)  +
     &           d100*sn(lps,j0+1,k0,i0) + d101*sn(lps,j0+1,k0,i0+1) +
     &           d010*sn(lps,j0,k0+1,i0) + d011*sn(lps,j0,k0+1,i0+1) +
     &       d110*sn(lps,j0+1,k0+1,i0) + d111*sn(lps,j0+1,k0+1,i0+1)
          tt = tt + 0.5*(sn0+snss)*sslen
C          write(*,*) tt, sn0, snss, sslen, y0,yss
C          write(*,*) "Length of the subsegment is ",sslen,icase
c          alllen=alllen+sslen
C
Ccccc   Add LttS = Vs^3/Vp^2/len
C         if (ifqk .eq. 1) then
C           sn02  = c000*sn(2,j0,k0,i0)   + c001*sn(2,j0,k0,i0+1) +
C    &           c100*sn(2,j0+1,k0,i0) + c101*sn(2,j0+1,k0,i0+1) +
C    &           c010*sn(2,j0,k0+1,i0) + c011*sn(2,j0,k0+1,i0+1) +
C    &         c110*sn(2,j0+1,k0+1,i0) + c111*sn(2,j0+1,k0+1,i0+1)
C           snss2  = d000*sn(2,j0,k0,i0)   + d001*sn(2,j0,k0,i0+1) +
C    &           d100*sn(2,j0+1,k0,i0) + d101*sn(2,j0+1,k0,i0+1) +
C    &           d010*sn(2,j0,k0+1,i0) + d011*sn(2,j0,k0+1,i0+1) +
C    &         d110*sn(2,j0+1,k0+1,i0) + d111*sn(2,j0+1,k0+1,i0+1)
C         LttS=LttS+((sn0**2)/(sn02**3)+(snss**2)/(snss2**3))/sslen/2.
C          write(*,*) LttS, tt, sn0, snss, sn02, snss2, sslen
C          write(*,*) "Length of the subsegment is ",sslen,icase
C         endif
Ccccc   Add LttS = Vs^3/Vp^2/len

          x0 = xss
          y0 = yss
          z0 = zss
          if (z0 .eq. zarc(npts+1)) pt = pt+1
          if (icase .eq. 1) then      ! only need to change endpt for case 1
            pt = pt + 1
            x1 = xarc(pt)
            y1 = yarc(pt)
            z1 = zarc(pt)
          endif 
          if (jcount .gt. 1000000) pt = pt+1
          jcount=jcount+1
        enddo ! end while for tt summation

        tt = tt + tcrust
C       if (ifqk .eq. 1) then
C         LttS=LttS+3.5**3/(6.101**2)/(sqrt((hdist-d)**2+zcrust**2))
C       endif
        ttchng = tt0 - tt
        if (tt0.lt.tt.and.abs(ttchng).gt.epsilon.and.idbl.eq.0) then
c          write (*,*) "WARNING: raypath not converging", ttchng
          f = f*.98
        endif
        tt0 = tt
        if (abs(ttchng) .gt. epsilon) idbl = 0
c        write(*,*) ttchng,epsilon,idbl
c        write(*,*) "end of one loop of pseudo-bending"
      enddo ! loop for pseudo-bending
c     END OF LOOP FOR PSEUDO-BENDING
c      write(*,*) "end of pseudo-bending with certain points"
c      write(*,*) 'ttchng = ',ttchng,idbl

      if (npts*2 .ge. maxpts) then
        WRITE (*,*) "WARNING: Max pts exceeded. Pseudo-bending stopped"
        idbl = 1
      endif 
      if (idbl .eq. 0) then      ! interpolate by factor of 2
        do i = npts+1, 1, -1
          xarc(i*2) = xarc(i)
          yarc(i*2) = yarc(i)
          zarc(i*2) = zarc(i)
        enddo
        npts = 2*npts + 1
        do i = 1, npts, 2
          xarc(i) = 0.5*(xarc(i-1) + xarc(i+1))
          yarc(i) = 0.5*(yarc(i-1) + yarc(i+1))
          zarc(i) = 0.5*(zarc(i-1) + zarc(i+1))
        enddo
        idbl = 1
c        write (*,*) " doubling # pts to ", npts
        go to 100
      endif

      tvalue=snss
      xarc(npts+2) = xs0            ! complete raypath
      yarc(npts+2) = ys0
      zarc(npts+2) = zs0
c      write(*,*) "alllen1 = :",alllen

666   END
