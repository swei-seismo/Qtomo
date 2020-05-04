      PROGRAM hits4gmt
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Program to convert hits file to gmt-friendly files    c
c     hits file: nodes are counted in the order of y, x     c
c     node=(j-1)*imax*kmax+(i-1)*kmax+k                     c
c     Velocity model: nodes are counted in the order of x, yc
c     by S. Wei, July 2013                                  c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      parameter (nxmax = 34, nymax = 14, nzmax = 30)

      integer*4 narguments,iargc,nnx,nny,nnz,dep(nzmax),minhit
      integer*4 hits(nxmax,nymax,nzmax),hits2(nxmax,nymax,nzmax)
      real lat(nxmax,nymax), lon(nxmax,nymax)
      character*30 hitsf,hitsf2,vmodel,outfile,maskfl
      character*70 temp
      character*1 pors
      character*3 layerid,nnzQ,minhits

      narguments = iargc()
      if (narguments.lt.4) then
        stop '"hits4gmt P/S/O/R Velocity-model nnzQ minhits"'
      endif
      call getarg(1,pors)
      call getarg(2,vmodel)
      call getarg(3,nnzQ)
      call getarg(4,minhits)
      if (pors.eq.'R') then
        hitsf='hitsP'
        hitsf2='hitsS'
      else
        hitsf='hits'//pors
      endif
      read(nnzQ,'(i2)') nnz
      read(minhits,'(i3)') minhit

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

c     Input hits file
      open(12,file=hitsf,status='old')
      do i = 1, nnx
        do j = 1, nny
          read(12,*) (hits(i,j,k), k=1,nnz)
          do k=1,nnz
            hits2(i,j,k)=10*hits(i,j,k)
          enddo
        enddo
      enddo
      close(12)
c     Input another hits file for R
      if (pors.eq.'R') then
        open(22,file=hitsf2,status='old')
        do i = 1, nnx
          do j = 1, nny
            read(22,*) (hits2(i,j,k), k=1,nnz)
          enddo
        enddo
        close(22)
      endif

c     Output hits at each depth
      do k = 1, nnz
        write(layerid,'(i3.3)') dep(k)
        outfile='hits'//pors//'.dep'//layerid
        maskfl='mask'//pors//'.dep'//layerid
        write(*,*) outfile
        open(13,file=outfile,status='new')
        open(14,file=maskfl,status='new')
        do i = 1, nnx
          do j = 1, nny
            write(13,*) lon(i,j),lat(i,j),min(hits(i,j,k),hits2(i,j,k))
            if ((hits(i,j,k).ge.minhit).and.(hits2(i,j,k).ge.minhit))
     $                                               then
              write(14,*) lon(i,j), lat(i,j), '100'
            else
              write(14,*) lon(i,j), lat(i,j), 'NaN'
            endif
          enddo
        enddo
        close(13)
        close(14)
      enddo

      END
