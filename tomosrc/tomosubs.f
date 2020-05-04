      SUBROUTINE lonlat2xy(lon0,lat0,beta,lon1,lat1,x1,y1)
c     Convert lon, lat to Cartesian coordinates
c     lon0, lat0: reference location
c     beta: angle of x axis from east (counterclockwise)
      real lon0,lat0,lon1,lat1,x1,y1,beta
      real xx, yy
c      common/b/ ipolar, nnx, nny, nnz, re, pi, r2d, nnzQ
      common /b/ ipolar,nnx,nny,nnz,re,pi,r2d,nnxQ,nnyQ,nnzQ,nnzQ1D
      xx = (lon1-lon0)*re/r2d
      yy = (lat1-lat0)*re/r2d
      x1 = (xx-yy*tan(beta))*cos(beta)
      y1 = x1*tan(beta)+yy/cos(beta)
      if (abs(x1) .lt. 0.01) x1 = 0
      if (abs(y1) .lt. 0.01) y1 = 0
      return
      END

      SUBROUTINE xy2lonlat(lon0,lat0,beta,x2,y2,lon2,lat2)
c     Convert x, y to lon, lat
c     lon0, lat0: reference location
c     beta: angle of x axis from east (counterclockwise)
      real lon0,lat0,lon2,lat2,x2,y2,beta
      real xx, yy
c      common/b/ ipolar, nnx, nny, nnz, re, pi, r2d, nnzQ
      common /b/ ipolar,nnx,nny,nnz,re,pi,r2d,nnxQ,nnyQ,nnzQ,nnzQ1D
      yy = (y2-x2*tan(beta))*cos(beta)
      xx = yy*tan(beta)+x2/cos(beta)
      lon2 = xx*r2d/re+lon0
      lat2 = yy*r2d/re+lat0
      return
      END


      REAL FUNCTION piercept(d,z1,z2,v1,v2)     
ccccc -- find where ray crosses from 1 isotropic layer to next.
ccccc -- d = horizontal distance between source & receiver
ccccc -- z1, z2 = layer thicknesses (receiver at surface of layer 1) 
ccccc -- Uses method of steepest descent.  alpha ~ 20 seems to work well
ccccc -- OUTPUT: distance in horizontal km from source
      real d, z1, z2, v1, v2

      epsilon = 5.e-4   ! dtdr threshold
      alpha = 20.
      u1 = 1./v1        ! slowness of layer 1
      u2 = 1./v2        ! slowness of layer 2
      dr = 0.1          ! needed for dt/dr
      r  = 0.5*d        ! initial guess is halfway between source & receiver
      i = 0

5     t = u2*sqrt(r*r + z2*z2) + u1*sqrt((d-r)**2 + z1*z1)
      dtdr = (u2*sqrt((r+dr)**2 + z2*z2) + 
     &           u1*sqrt((d-(r+dr))**2 + z1*z1) - t )/dr

      r = r - alpha*dtdr
      i = i + 1
      if (mod(i,200) .eq. 0) then
        alpha = 0.5*alpha
        go to 5
      endif
      if (abs(dtdr).gt.epsilon .and. i.lt.5000) go to 5
*     write (*,*) i, " steps"

      piercept = r            ! distance in horizontal km from source
      return
      END
