      subroutine gcpath_e(slat,slon,sbear,dist,ndist,gcloc)
*=====================================================================
* PURPOSE:  To compute locations along a great circle path.
*           assuming elliptical Earth
*=====================================================================
* INPUT ARGUMENTS:
*    SLAT:    Source latitude in degrees, north positive. [f]
*    SLON:    Source longtiude in degrees, east positive. [f]
*    SBEAR:   Source bearing (i.e. local azimuth along path.) [f]
*    DIST:    Array of distances from source in km along path. [f]
*    NDIST:   Length of DIST array. [i]
*=====================================================================
* OUTPUT ARGUMENTS:
*    GCLOC:   Array containing the latitudes and longitudes of GCP
*             line segments.  The first row contains the latitudes,
*             the second row the longitudes. [fa]
*=====================================================================
* MODULE/LEVEL:  MAP/4
*=====================================================================
* GLOBAL INPUT:
*    MAP:     RADIUS, TORAD, TODEG
*=====================================================================
* REFERENCES:
* - Derivation of the spherical geomety equations used in this sub-
*   routine are due to Dave Harris and can be found in the MAP folder.
*=====================================================================

*====== 
* ADDITION TO DEAL WITH GEOCENTRIC COORDINATES
*======

      dimension dist(*), gcloc(2,*)
      real*4 ex, ey, ez, alpha, beta
*
      include 'numerical.h'
*
      e = 1.0/flt
      a = (1.0 - e)**2
      b = 1.0/a

* PROCEDURE:

* - Compute constants that depend only upon the source location.

      theta=drad*slat
*
      beta = atan( a * tan(theta))
      theta = beta
*
      ct=cos(theta)
      st=sin(theta)

      phi=drad*slon
      cp=cos(phi)
      sp=sin(phi)

      bear=drad*sbear
      cb=cos(bear)
      sb=sin(bear)

* - Loop on each element in distance array.

      do 1000 j=1,ndist

        delta=dist(j)/rnk
        cd=cos(delta)
        sd=sin(delta)

        ez=cd*st+sd*cb*ct
        ey=cd*ct*cp+sd*(-cb*st*cp-sb*sp)
        ex=cd*ct*sp+sd*(-cb*st*sp+sb*cp)

c       beta = atan2(ez,sqrt(ex*ex+ey*ey))
c       alpha = atan( b * tan(beta) ) 

        beta = ez/sqrt(ex*ex+ey*ey)
        alpha = atan( b * beta ) 

        gcloc(2,j)=radd*atan2(ex,ey)
 1000   gcloc(1,j)=radd*alpha

 8888 return

*=====================================================================
* MODIFICATION HISTORY:
*    870105:  SUN version.
*    850430:  Combined the latitudes and longitudes of the GCP line
*             segments into one array.
*    840807:  Original version (due to Dave Harris.)
*=====================================================================
* DOCUMENTED/REVIEWED:  850212
*=====================================================================

      end
