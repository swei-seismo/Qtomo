cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   SUBROUTINE TO WRITE CONSTRAINT ROWS TO G-MATRIX  ccc
ccc   2nd ORDER SMOOTHING CONTRAINTS APPEDNING TO G    ccc
ccc                                                    ccc
ccc    INPUT:  smo,nmod,nd                             ccc
ccc    OUTPUT: DV,WD,GB,nk                             ccc
ccc                                                    ccc
ccc    S. Wei, June 2013, edited from                  ccc
ccc      J. A. Conder and S.H.Pozgay's 2D version      ccc
ccc    Edited in Nov. 2015 to allow varying spacing    ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE smoothg(smo,nmod,nd,DV,WD,GB,nk)

      parameter (nxmax = 45, nymax = 45, nzmax = 30)  ! max nodes
      integer nmod              ! model dimensions (as in calling prgm)
      integer nd                ! data dimensions (as in calling prgm)
      integer nk                ! number of new rows added here
      integer nnx, nny, nnz, nnxQ, nnyQ, nnzQ, nnzQ1D      ! model nodes
      real x(nxmax), z(nzmax), y(nymax)
      real xq(nxmax), zq(nzmax), yq(nymax), zq1D(nzmax), attin1D(nzmax)
      real smo,smleng,vu,vd,h   ! smoothing weight
      dimension GR(nmod)        ! row of GB
      dimension DV(nd)          ! data vector needed to be added
      dimension WD(nd)          ! weight vector needed to be added
      dimension GB(nd,nmod)     ! smooth matrix needed to be added
      real c0,cu,cd,cw,ce,cs,cn ! smoothing cefficients

      common /a/ x, y, z, xq, yq, zq
      common /b/ ipolar,nnx,nny,nnz,re,pi,r2d,nnxQ,nnyQ,nnzQ,nnzQ1D

      nk = 0
      c0 = 1.           ! smoothing constraints

cccc  INTERIOR NODES  cccc
      smleng = smo
      do i = 2, nnxQ-1
        do j = 2, nnyQ-1
          do k = 2, nnzQ-1
            do node = 1, nnxQ*nnyQ*nnzQ
              GR(node) = 0.
            enddo
            vu = abs(zq(k)-zq(k-1))
            vd = abs(zq(k)-zq(k+1))
            h  = abs(xq(i)-xq(i-1))
            cu = -c0/(h/vu+h/vd+4)*h/vu
            cd = -c0/(h/vu+h/vd+4)*h/vd
            cw = -c0/(h/vu+h/vd+4)
            ce = -c0/(h/vu+h/vd+4)
            cs = -c0/(h/vu+h/vd+4)
            cn = -c0/(h/vu+h/vd+4)
            GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
            GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
            GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
            GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
            GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
            GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
            GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
            nk = nk + 1
            DV(nk) = 0.
            WD(nk) = smleng
            do node = 1,nnxQ*nnyQ*nnzQ
              GB(nk,node) = GR(node)
            enddo
          enddo
        enddo
      enddo
      write(*,*) "Inter nodes"

cccc  SIDE NODES  cccc
      smleng = smo*7/6
      k = 1       ! top side
      do i = 2, nnxQ-1
        do j = 2, nnyQ-1
          do node = 1, nnxQ*nnyQ*nnzQ
            GR(node) = 0.
          enddo
          vd = abs(zq(k)-zq(k+1))
          h  = abs(xq(i)-xq(i-1))
          cd = -c0/(h/vd+4)*h/vd
          cw = -c0/(h/vd+4)
          ce = -c0/(h/vd+4)
          cs = -c0/(h/vd+4)
          cn = -c0/(h/vd+4)
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
          GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
          GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
          GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
          GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
          nk = nk + 1
          DV(nk) = 0.
          WD(nk) = smleng
          do node = 1,nnxQ*nnyQ*nnzQ
            GB(nk,node) = GR(node)
          enddo
        enddo
      enddo

      k = nnzQ    ! bottom side
      do i = 2, nnxQ-1
        do j = 2, nnyQ-1
          do node = 1, nnxQ*nnyQ*nnzQ
            GR(node) = 0.
          enddo
          vu = abs(zq(k)-zq(k-1))
          h  = abs(xq(i)-xq(i-1))
          cu = -c0/(h/vu+4)*h/vu
          cw = -c0/(h/vu+4)
          ce = -c0/(h/vu+4)
          cs = -c0/(h/vu+4)
          cn = -c0/(h/vu+4)
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
          GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
          GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
          GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
          GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
          nk = nk + 1
          DV(nk) = 0.
          WD(nk) = smleng
          do node = 1,nnxQ*nnyQ*nnzQ
            GB(nk,node) = GR(node)
          enddo
        enddo
      enddo

      i = 1       ! west side
      do j = 2, nnyQ-1
        do k = 2, nnzQ-1
          do node = 1, nnxQ*nnyQ*nnzQ
            GR(node) = 0.
          enddo
          vu = abs(zq(k)-zq(k-1))
          vd = abs(zq(k)-zq(k+1))
          h  = abs(xq(i)-xq(i+1))
          cu = -c0/(h/vu+h/vd+3)*h/vu
          cd = -c0/(h/vu+h/vd+3)*h/vd
          ce = -c0/(h/vu+h/vd+3)
          cs = -c0/(h/vu+h/vd+3)
          cn = -c0/(h/vu+h/vd+3)
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
          GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
          GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
          GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
          nk = nk + 1
          DV(nk) = 0.
          WD(nk) = smleng
          do node = 1,nnxQ*nnyQ*nnzQ
            GB(nk,node) = GR(node)
          enddo
        enddo
      enddo

      i = nnxQ     ! east side
      do j = 2, nnyQ-1
        do k = 2, nnzQ-1
          do node = 1, nnxQ*nnyQ*nnzQ
            GR(node) = 0.
          enddo
          vu = abs(zq(k)-zq(k-1))
          vd = abs(zq(k)-zq(k+1))
          h  = abs(xq(i)-xq(i-1))
          cu = -c0/(h/vu+h/vd+3)*h/vu
          cd = -c0/(h/vu+h/vd+3)*h/vd
          cw = -c0/(h/vu+h/vd+3)
          cs = -c0/(h/vu+h/vd+3)
          cn = -c0/(h/vu+h/vd+3)
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
          GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
          GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
          GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
          nk = nk + 1
          DV(nk) = 0.
          WD(nk) = smleng
          do node = 1,nnxQ*nnyQ*nnzQ
            GB(nk,node) = GR(node)
          enddo
        enddo
      enddo

      j = 1       ! south side
      do i = 2, nnxQ-1
        do k = 2, nnzQ-1
          do node = 1, nnxQ*nnyQ*nnzQ
            GR(node) = 0.
          enddo
          vu = abs(zq(k)-zq(k-1))
          vd = abs(zq(k)-zq(k+1))
          h  = abs(xq(i)-xq(i-1))
          cu = -c0/(h/vu+h/vd+3)*h/vu
          cd = -c0/(h/vu+h/vd+3)*h/vd
          cw = -c0/(h/vu+h/vd+3)
          ce = -c0/(h/vu+h/vd+3)
          cn = -c0/(h/vu+h/vd+3)
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
          GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
          GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
          GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
          nk = nk + 1
          DV(nk) = 0.
          WD(nk) = smleng
          do node = 1,nnxQ*nnyQ*nnzQ
            GB(nk,node) = GR(node)
          enddo
        enddo
      enddo

      j = nnyQ     ! north side
      do i = 2, nnxQ-1
        do k = 2, nnzQ-1
          do node = 1, nnxQ*nnyQ*nnzQ
            GR(node) = 0.
          enddo
          vu = abs(zq(k)-zq(k-1))
          vd = abs(zq(k)-zq(k+1))
          h  = abs(xq(i)-xq(i-1))
          cu = -c0/(h/vu+h/vd+3)*h/vu
          cd = -c0/(h/vu+h/vd+3)*h/vd
          cw = -c0/(h/vu+h/vd+3)
          ce = -c0/(h/vu+h/vd+3)
          cs = -c0/(h/vu+h/vd+3)
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
          GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
          GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
          GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
          GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
          nk = nk + 1
          DV(nk) = 0.
          WD(nk) = smleng
          do node = 1,nnxQ*nnyQ*nnzQ
            GB(nk,node) = GR(node)
          enddo
        enddo
      enddo
      write(*,*) "side nodes"

cccc  EDGE NODES  cccc
      smleng = smo*7/5

      j = 1       ! top-south edge
      k = 1
      do i = 2, nnxQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0.
        enddo
        vd = abs(zq(k)-zq(k+1))
        h  = abs(xq(i)-xq(i-1))
        cd = -c0/(h/vd+3)*h/vd
        cw = -c0/(h/vd+3)
        ce = -c0/(h/vd+3)
        cn = -c0/(h/vd+3)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
        GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
        GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
        GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo

      j = nnyQ     ! top-north edge
      k = 1
      do i = 2, nnxQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0.
        enddo
        vd = abs(zq(k)-zq(k+1))
        h  = abs(xq(i)-xq(i-1))
        cd = -c0/(h/vd+3)*h/vd
        cw = -c0/(h/vd+3)
        ce = -c0/(h/vd+3)
        cs = -c0/(h/vd+3)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
        GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
        GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
        GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo

      i = 1       ! top-west edge
      k = 1
      do j = 2, nnyQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0.
        enddo
        vd = abs(zq(k)-zq(k+1))
        h  = abs(xq(i)-xq(i+1))
        cd = -c0/(h/vd+3)*h/vd
        ce = -c0/(h/vd+3)
        cs = -c0/(h/vd+3)
        cn = -c0/(h/vd+3)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
        GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
        GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
        GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo

      i = nnxQ     ! top-east edge
      k = 1
      do j = 2, nnyQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0.
        enddo
        vd = abs(zq(k)-zq(k+1))
        h  = abs(xq(i)-xq(i-1))
        cd = -c0/(h/vd+3)*h/vd
        cw = -c0/(h/vd+3)
        cn = -c0/(h/vd+3)
        cs = -c0/(h/vd+3)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
        GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
        GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
        GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo

      j = 1       ! bottom-south edge
      k = nnzQ
      do i = 2, nnxQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0.
        enddo
        vu = abs(zq(k)-zq(k-1))
        h  = abs(xq(i)-xq(i-1))
        cu = -c0/(h/vu+3)*h/vu
        cw = -c0/(h/vu+3)
        ce = -c0/(h/vu+3)
        cn = -c0/(h/vu+3)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
        GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
        GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
        GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo

      j = nnyQ     ! bottom-north edge
      k = nnzQ
      do i = 2, nnxQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0.
        enddo
        vu = abs(zq(k)-zq(k-1))
        h  = abs(xq(i)-xq(i-1))
        cu = -c0/(h/vu+3)*h/vu
        cw = -c0/(h/vu+3)
        ce = -c0/(h/vu+3)
        cs = -c0/(h/vu+3)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
        GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
        GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
        GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo

      i = 1       ! bottom-west edge
      k = nnzQ
      do j = 2, nnyQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0.
        enddo
        vu = abs(zq(k)-zq(k-1))
        h  = abs(xq(i)-xq(i+1))
        cu = -c0/(h/vu+3)*h/vu
        cs = -c0/(h/vu+3)
        ce = -c0/(h/vu+3)
        cn = -c0/(h/vu+3)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
        GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
        GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
        GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo

      i = nnxQ     ! bottom-east edge
      k = nnzQ
      do j = 2, nnyQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0. 
        enddo
        vu = abs(zq(k)-zq(k-1))
        h  = abs(xq(i)-xq(i-1))
        cu = -c0/(h/vu+3)*h/vu
        cw = -c0/(h/vu+3)
        cs = -c0/(h/vu+3)
        cn = -c0/(h/vu+3)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
        GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
        GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
        GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo

      i = 1       ! south-west edge
      j = 1
      do k = 2, nnzQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0.
        enddo
        vu = abs(zq(k)-zq(k-1))
        vd = abs(zq(k)-zq(k+1))
        h  = abs(xq(i)-xq(i+1))
        cu = -c0/(h/vu+h/vd+2)*h/vu
        cd = -c0/(h/vu+h/vd+2)*h/vd
        ce = -c0/(h/vu+h/vd+2)
        cn = -c0/(h/vu+h/vd+2)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
        GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
        GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo

      i = nnxQ     ! south-east edge
      j = 1
      do k = 2, nnzQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0. 
        enddo
        vu = abs(zq(k)-zq(k-1))
        vd = abs(zq(k)-zq(k+1))
        h  = abs(xq(i)-xq(i-1))
        cu = -c0/(h/vu+h/vd+2)*h/vu
        cd = -c0/(h/vu+h/vd+2)*h/vd
        cn = -c0/(h/vu+h/vd+2)
        cw = -c0/(h/vu+h/vd+2)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
        GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
        GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo

      i = 1       ! north-west edge
      j = nnyQ
      do k = 2, nnzQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0. 
        enddo
        vu = abs(zq(k)-zq(k-1))
        vd = abs(zq(k)-zq(k+1))
        h  = abs(xq(i)-xq(i+1))
        cu = -c0/(h/vu+h/vd+2)*h/vu
        cd = -c0/(h/vu+h/vd+2)*h/vd
        ce = -c0/(h/vu+h/vd+2)
        cs = -c0/(h/vu+h/vd+2)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
        GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
        GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo

      i = nnxQ     ! north-east edge
      j = nnyQ
      do k = 2, nnzQ-1
        do node = 1, nnxQ*nnyQ*nnzQ
          GR(node) = 0. 
        enddo
        vu = abs(zq(k)-zq(k-1))
        vd = abs(zq(k)-zq(k+1))
        h  = abs(xq(i)-xq(i-1))
        cu = -c0/(h/vu+h/vd+2)*h/vu
        cd = -c0/(h/vu+h/vd+2)*h/vd
        cw = -c0/(h/vu+h/vd+2)
        cs = -c0/(h/vu+h/vd+2)
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
        GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
        GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
        GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
        nk = nk + 1
        DV(nk) = 0.
        WD(nk) = smleng
        do node = 1,nnxQ*nnyQ*nnzQ
          GB(nk,node) = GR(node)
        enddo
      enddo
      write(*,*)"edge nodes"

cccc  CORNER NODES  cccc
      smleng = smo*7/4

      i = 1       ! top-south-west corner
      j = 1
      k = 1
      do node = 1, nnxQ*nnyQ*nnzQ
        GR(node) = 0.
      enddo
      vd = abs(zq(k)-zq(k+1))
      h  = abs(xq(i)-xq(i+1))
      cd = -c0/(h/vd+2)*h/vd
      ce = -c0/(h/vd+2)
      cn = -c0/(h/vd+2)
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
      GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
      GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
      nk = nk + 1
      DV(nk) = 0.
      WD(nk) = smleng
      do node = 1,nnxQ*nnyQ*nnzQ
        GB(nk,node) = GR(node)
      enddo

      i = 1       ! bottom-south-west corner
      j = 1
      k = nnzQ
      do node = 1, nnxQ*nnyQ*nnzQ
        GR(node) = 0.
      enddo
      vu = abs(zq(k)-zq(k-1))
      h  = abs(xq(i)-xq(i+1))
      cu = -c0/(h/vu+2)*h/vu
      ce = -c0/(h/vu+2)
      cn = -c0/(h/vu+2)
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
      GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
      GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
      nk = nk + 1
      DV(nk) = 0.
      WD(nk) = smleng
      do node = 1,nnxQ*nnyQ*nnzQ
        GB(nk,node) = GR(node)
      enddo

      i = 1       ! top-north-west corner
      j = nnyQ
      k = 1
      do node = 1, nnxQ*nnyQ*nnzQ
        GR(node) = 0.
      enddo
      vd = abs(zq(k)-zq(k+1))
      h  = abs(xq(i)-xq(i+1))
      cd = -c0/(h/vd+2)*h/vd
      ce = -c0/(h/vd+2)
      cs = -c0/(h/vd+2)
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
      GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
      GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
      nk = nk + 1
      DV(nk) = 0.
      WD(nk) = smleng
      do node = 1,nnxQ*nnyQ*nnzQ
        GB(nk,node) = GR(node)
      enddo

      i = 1       ! bottow-north-west corner
      j = nnyQ
      k = nnzQ
      do node = 1, nnxQ*nnyQ*nnzQ
        GR(node) = 0.
      enddo
      vu = abs(zq(k)-zq(k-1))
      h  = abs(xq(i)-xq(i+1))
      cu = -c0/(h/vu+2)*h/vu
      ce = -c0/(h/vu+2)
      cs = -c0/(h/vu+2)
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
      GR((i-1)*nnzQ*nnyQ+j*nnzQ+k) = ce          ! east node
      GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
      nk = nk + 1
      DV(nk) = 0.
      WD(nk) = smleng
      do node = 1,nnxQ*nnyQ*nnzQ
        GB(nk,node) = GR(node)
      enddo

      i = nnxQ     ! top-south-east corner
      j = 1
      k = 1
      do node = 1, nnxQ*nnyQ*nnzQ
        GR(node) = 0.
      enddo
      vd = abs(zq(k)-zq(k+1))
      h  = abs(xq(i)-xq(i-1))
      cd = -c0/(h/vd+2)*h/vd
      cw = -c0/(h/vd+2)
      cn = -c0/(h/vd+2)
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
      GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
      GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
      nk = nk + 1
      DV(nk) = 0.
      WD(nk) = smleng
      do node = 1,nnxQ*nnyQ*nnzQ
        GB(nk,node) = GR(node)
      enddo

      i = nnxQ     ! bottom-south-east corner
      j = 1
      k = nnzQ
      do node = 1, nnxQ*nnyQ*nnzQ
        GR(node) = 0.
      enddo
      vu = abs(zq(k)-zq(k-1))
      h  = abs(xq(i)-xq(i-1))
      cu = -c0/(h/vu+2)*h/vu
      cw = -c0/(h/vu+2)
      cn = -c0/(h/vu+2)
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
      GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
      GR(i*nnzQ*nnyQ+(j-1)*nnzQ+k) = cn          ! north node
      nk = nk + 1
      DV(nk) = 0.
      WD(nk) = smleng
      do node = 1,nnxQ*nnyQ*nnzQ
        GB(nk,node) = GR(node)
      enddo

      i = nnxQ     ! top-north-east corner
      j = nnyQ
      k = 1
      do node = 1, nnxQ*nnyQ*nnzQ
        GR(node) = 0.
      enddo
      vd = abs(zq(k)-zq(k+1))
      h  = abs(xq(i)-xq(i-1))
      cd = -c0/(h/vd+2)*h/vd
      cw = -c0/(h/vd+2)
      cs = -c0/(h/vd+2)
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k+1) = cd    ! lower node
      GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
      GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
      nk = nk + 1
      DV(nk) = 0.
      WD(nk) = smleng
      do node = 1,nnxQ*nnyQ*nnzQ
        GB(nk,node) = GR(node)
      enddo

      i = nnxQ     ! bottow-north-east corner
      j = nnyQ
      k = nnzQ
      do node = 1, nnxQ*nnyQ*nnzQ
        GR(node) = 0.
      enddo
      vu = abs(zq(k)-zq(k-1))
      h  = abs(xq(i)-xq(i-1))
      cu = -c0/(h/vu+2)*h/vu
      cw = -c0/(h/vu+2)
      cs = -c0/(h/vu+2)
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k) = c0      ! target node
      GR((i-1)*nnzQ*nnyQ+(j-1)*nnzQ+k-1) = cu    ! upper node
      GR((i-1)*nnzQ*nnyQ+(j-2)*nnzQ+k) = cw      ! west node
      GR((i-2)*nnzQ*nnyQ+(j-1)*nnzQ+k) = cs      ! south node
      nk = nk + 1
      DV(nk) = 0.
      WD(nk) = smleng
      do node = 1,nnxQ*nnyQ*nnzQ
        GB(nk,node) = GR(node)
      enddo

      END
