      subroutine cubic(jcom,modes)
c
c     cubic interpolation of eigenfunctions to the source depth
c
      include 'parameter5.h'
c
      real*8 a1,a2,c0,c1,c2,c3,x,y,y2,y3
      common/scl/x1,r0,x2,f(4,3,maxmodes)
      data isw/1/
      if (jcom .eq. 3) then
        itype = 2
      else
        itype = 1
      endif
      do nm = 1, modes
        if(isw.ne.1)go to 10
        isw=0
        y=x2-x1
        y2=y**2
        y3=y*y2
        x=r0-x1
        a1=3./y2
        a2=2./y3
   10   continue
        do 20 i=1,itype
           k=2*i
           j=k-1
           c0=f(j,1,nm)
           c1=f(k,1,nm)
           c3=f(j,3,nm)-c0
           c2=a1*c3-(2.*c1+f(k,3,nm))/y
           c3=(c1+f(k,3,nm))/y2-a2*c3
           f(j,2,nm)=c0+x*(c1+x*(c2+x*c3))
           f(k,2,nm)=c1+x*(2.*c2+3.*x*c3)
   20   continue
      end do
      return
      end
