c
      function fcn(p)
      real   l1,l2,l3,mm1,mm2,mm3
      common/b2/l1,l2,l3,mm1,mm2,mm3,factor
      r=p*factor
      g= abs(l3)
      if(g.lt.0.01) go to 1
      ang=-(l1* cos(r)+l2* sin(r))/l3
      ang= atan(ang)
      fcn=ang/factor
      return
    1 ang=-(mm1* cos(r)+mm2* sin(r))/mm3
      ang= atan(ang)
      fcn=ang/factor
      return
      end
