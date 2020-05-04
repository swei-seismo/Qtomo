      PROGRAM atten4gmt
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Program to convert M matrix to gmt-friendly files     c
c     M matrix: nodes are counted in the order of z, y, x   c
c     node=(j-1)*imax*kmax+(i-1)*kmax+k                     c
c     Velocity model: nodes are counted in the order of x, yc
c     by S. Wei, July 2013                                  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      parameter (nxmax = 34, nymax = 14, nzmax = 30)

      integer*4 narguments,iargc,nnx,nny,nnz,dep(nzmax)
      real lat(nxmax,nymax), lon(nxmax,nymax)
      real qinv(nxmax,nymax,nzmax), tomo(nxmax,nymax,nzmax)
      real qinv2(nxmax,nymax,nzmax), QpQs(nxmax,nymax,nzmax)
      character*30 qmodel,vmodel,outfile,qmodel2
      character*70 temp
      character*3 layerid,nnzQ
      character*1 pors

      narguments = iargc()
      if (narguments.lt.4) then
        stop '"atten4gmt Q-model Velocity-model nnzQ P/S/I/O/R"'
      endif
      call getarg(1,qmodel)
      call getarg(2,vmodel)
      call getarg(3,nnzQ)
      call getarg(4,pors)
      if (pors.eq.'p') then
        pors='P'
      elseif (pors.eq.'i') then
        pors='I'
      elseif (pors.eq.'o') then
        pors='O'
      elseif (pors.eq.'r') then
        pors='R'
      endif
      if (pors.eq.'R') then
        call getarg(5,qmodel2)
      endif
      read(nnzQ,'(i2)') nnz

c     Input Velocity model
      open(11,file=vmodel,status='old')
      read(11,*) nnx,nny,temp
      read(11,*) temp
      do i = 1, nnx
        do j = 1, nny
          read(11,*) lat(i,j), lon(i,j)
c          if (lon(i,j).gt.180) then
c            lon(i,j) = lon(i,j) - 360
c          endif
        enddo
      enddo
      read(11,*) (dep(k), k=1,nnz)
      close(11)

c     Input Q model
      open(12,file=qmodel,status='old')
      do i = 1, nnx
        do j = 1, nny
          do k = 1, nnz
            read(12,*) qinv(i,j,k)
            if (qinv(i,j,k).le.0) then
              qinv(i,j,k)=0.00001
            endif
            tomo(i,j,k) = 1000 * qinv(i,j,k)
          enddo
        enddo
      enddo
      close(12)

c     Input Qp model for R (Qp/Qs)
      if (pors.eq.'R') then
        open(22,file=qmodel2,status='old')
        do i = 1, nnx
          do j = 1, nny
            do k = 1, nnz
              read(22,*) qinv2(i,j,k)
              if (qinv2(i,j,k).le.0) then
                qinv2(i,j,k)=0.00001
              endif
              QpQs(i,j,k) = qinv(i,j,k)/qinv2(i,j,k)
            enddo
          enddo
        enddo
        close(22)
      endif

c     Output 1000/Q at each depth
      do k = 1, nnz
        write(layerid,'(i3.3)') dep(k)
        outfile='atten'//pors//'.dep'//layerid
        write(*,*) outfile
        open(13,file=outfile,status='new')
        do i = 1, nnx
          do j = 1, nny
            if (pors.eq.'R') then
              write(13,*) lon(i,j), lat(i,j), QpQs(i,j,k)
            else
              write(13,*) lon(i,j), lat(i,j), tomo(i,j,k)
            endif
          enddo
        enddo
        close(13)
      enddo

      END
