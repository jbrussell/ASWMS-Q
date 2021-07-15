      subroutine load_crust51
c bastardized version of program getCN2point (top half), from Gabi Laske
c reads in Laske and Mooney crustal model crust2.0 (update of
c crust5.1 -- 2 degree spacing), sets up variables for use 
c in avg_crust20.f.
c
c based on load_crust51, rather than getCN2point directl
c
c retains crust51 as the name for this subroutine (internally),
c the common blocks, and the variable names within the common.
c This allows the routine to be updated without changing the main
c program.  The file name records the version number.
c
c outputs:  common block crust20 with variables:
c               amapvp(8,nlo,nla),amaprho(8,nlo,nla),
c     +          amapvs(8,nlo,nla),amapthi(7,nlo,nla),
c     +          amapele(nlo,nla)
c
c JBG 6/02
c        
c layer one and two flipped, after the read statement!
c layer 1: water
c layer 2: ice


      parameter(maxdel = 180)
      real*4 delt(maxdel),gcloc(2,maxdel)

      parameter(ityp=360)
      parameter(nla=90,nlo=180)
      dimension fvel(ityp,8),fvels(ityp,8),frho(ityp,8)
      dimension fthi(ityp,7)
      character*2 ctype(ityp),line*506,dum*1,dum0*5
      character*2 types(nlo),nsta*4,ntyp*2,atype(nlo,nla)
      character*12 names(7)
      data names/'water','ice','soft sed.','hard sed.',
     +         'upper crust','middle crust','lower crust'/

c     common block of crustal parameters to return

      real*4 amapvp(8,nlo,nla),amaprho(8,nlo,nla),
     +          amapvs(8,nlo,nla),amapthi(7,nlo,nla),
     +          amapele(nlo,nla)

      common /crust51/ amapvp,amaprho,amapvs,amapthi,
     +                 amapele

      include 'numerical.h'

      print*,' '
      print*,'Loading crust2.0 parameters'
      print*,' '

      open(2,file=
     +    '/Users/gaherty/Unix/Data/datalib/crust2.0/CNtype2_key.txt')
      open(7,file=
     +    '/Users/gaherty/Unix/Data/datalib/crust2.0/CNtype2.txt')
      open(8,file=
     +    '/Users/gaherty/Unix/Data/datalib/crust2.0/CNelevatio2.txt')

      dx=360./nlo
c... read in key for crust types -- this is all original Laske code
c...............................
      read(2,890)dum
      print*,' ... reading key file ...'
      do 101 i=1,ityp
         read(2,899)ctype(i)
c        print 899,ctype(i)
         read(2,899)line
         read(line,*)(fvel(i,l),l=1,8)
         read(2,899)line
         read(line,*)(fvels(i,l),l=1,8)
         read(2,899)line
         read(line,*)(frho(i,l),l=1,8)
         read(2,899)line
         read(line,*)(fthi(i,l),l=1,7)
c flip layers
         aux=fvel(i,1)
         fvel(i,1)=fvel(i,2)
         fvel(i,2)=aux
         aux=fvels(i,1)
         fvels(i,1)=fvels(i,2)
         fvels(i,2)=aux
         aux=frho(i,1)
         frho(i,1)=frho(i,2)
         frho(i,2)=aux
         aux=fthi(i,1)
         fthi(i,1)=fthi(i,2)
         fthi(i,2)=aux
 101  continue

c... read CNtype file -- also original Laske code
c...............................
      read(7,*)flons
      print*,' ... reading model ...'
      read(8,899)line
      do 40 j=1,nla
         read(8,*)ilat,(amapele(i,j),i=1,nlo)
         read(7,901)ilat,types
c        print*,ilat
         do 10 i=1,nlo
            do 20 l=1,ityp
            if(types(i).eq.ctype(l))then
              atype(i,j)=ctype(l)
              do 30 k=1,8
              amapvp(k,i,j)=fvel(l,k)
              amapvs(k,i,j)=fvels(l,k)
              amaprho(k,i,j)=frho(l,k)
 30           continue
              do 31 k=1,7
 31           amapthi(k,i,j)=fthi(l,k)
              goto 10
            endif
 20         continue
            print*,' crust type code not found: ',types(i)
            print*,' latitude: ',ilat,' long index: ',i
 10      continue
 40   continue


 890  format(////a)
 899  format(a)
 901  format(i4,1x,180(2x,a2,1x))

 99   continue
      close(2)
      close(7)
      close(8)

      return

      end
