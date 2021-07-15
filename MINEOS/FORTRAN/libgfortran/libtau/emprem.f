      subroutine emdlv(r,vp,vs)
c         set up information on earth model (specified by
c         subroutine call emiasp)
c         set dimension of cpr,rd  equal to number of radial
c         discontinuities in model
      save
      character*(*) name
      character*20 modnam
      dimension cpr(13)
      common/emdlc/np,rd(13)
      data np,rd/13,1221.5,3480.0,3630.,5600.,5701.,5771.,5971.,
     1 6151.,6291.,6346.6,6356.,6368.,6371./
      data rn,vn/1.5696123e-4,6.8501006/
      data modnam/'prem'/
c
      call emiask(rn*r,rho,vp,vs)
      vp=vn*vp
      vs=vn*vs
      return
c
      entry emdld(n,cpr,name)
      n=np
      do 1 i=1,np
 1    cpr(i)=rd(i)
      name=modnam
      return
      end
c
      subroutine emiask(x0,ro,vp,vs)
c
c $$$$$ calls no other routine $$$$$
c
c   Emiask returns model parameters for the IASPEI working model 
c   (September 1990.1).  
c   Given non-dimensionalized radius x0, emiasp returns
c   non-dimensionalized density, ro, compressional velocity, vp, and
c   shear velocity, vs.  Non-dimensionalization is according to the
c   scheme of Gilbert in program EOS:  x0 by a (the radius of the
c   Earth), ro by robar (the mean density of the Earth), and velocity
c   by a*sqrt(pi*G*robar) (where G is the universal gravitational
c   constant.
c
c
      save
      dimension r(14),d(13,4),p(13,4),s(13,4)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      data r/0.0,1221.5,3480.0,3630.0,5600.0,5701.0,
     1    5771.0,5971.0,6151.0,6291.0,6346.6,6356.0,6368.0,6371.0/
      data d/13.0885 ,12.5815 , 7.9565 , 7.9565 , 7.9565 , 5.3197 ,
     1 11.2494 , 7.1089 ,  2.691 ,  2.691 , 2.90   , 2.6     , 1.02,
     2        0.     ,-1.2638 ,-6.4761 ,-6.4761 ,-6.4761 ,-1.4836 ,
     2 -8.0298 ,-3.8045 , 0.6924 ,  0.6924, 0.0    , 2*0.,
     3       -8.8381 ,-3.6426 , 5.5283 , 5.5283 , 5.5283 , 8*0.,
     4        0.     ,-5.5281 ,-3.0807 ,-3.0807 ,-3.0807 , 8*0./
      data p/11.2622 ,11.0487 ,15.3891 ,24.9520 ,29.2766 ,19.0957 ,
     1 39.7027 ,20.3926 , 4.1875 , 4.1875 , 6.8    , 5.8     , 5.80,
     2        0.     ,-4.0362 ,-5.3181 ,-40.4673,-23.6027,-9.8672 ,
     2-32.6166,-12.2569, 3.9382 , 3.9382 , 0.     , 2*0.,
     3       -6.3640 , 4.8023 , 5.5242 ,51.4832 , 5.5242 , 8*0.,
     4        0.     ,-13.5732,-2.5514 ,-26.6419,-2.5514 , 8*0./
      data s/ 3.6678 ,  0.0   , 6.9254 ,11.1671 ,22.3459 , 9.9839 ,
     122.3512 , 8.9496 , 2.1519 , 2.1519 , 3.9    , 3.2     , 3.2,
     2        0.0    ,  0.0   , 1.4672 ,-13.7818,-17.2473,-4.9324 ,
     2-18.5856,-4.4597 , 2.3481 , 2.3481 , 0.0    , 2*0.,
     3       -4.4475 ,  0.0   ,-2.0834 ,17.4575 ,-2.0834 , 8*0.,
     4        0.0    ,  0.0   , 0.9783 ,-9.2777 , 0.9783 , 8*0./
      data xn,rn,vn/6371.,.18125793,.14598326/,i/1/
c
c
c     do ii = 1, 13
c       do kk = 1, 4
c         write(6,'(4f10.5)') r(ii), d(ii,kk), p(ii,kk), s(ii,kk)
c       end do
c     end do
c
      x=amax1(x0,0.)
      x1=xn*x
 2    if(x1.ge.r(i)) go to 1
      i=i-1
      go to 2
 1    if(x1.le.r(i+1).or.i.ge.13) go to 3
      i=i+1
      if(i.lt.13) go to 1
 3    ro=rn*(d(i,1)+x*(d(i,2)+x*(d(i,3)+x*d(i,4))))
      vp=vn*(p(i,1)+x*(p(i,2)+x*(p(i,3)+x*p(i,4))))
      vs=vn*(s(i,1)+x*(s(i,2)+x*(s(i,3)+x*s(i,4))))
      return
      end
