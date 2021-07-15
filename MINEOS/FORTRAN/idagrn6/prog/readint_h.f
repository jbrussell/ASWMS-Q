c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine readint_h(hname,int_file)
c-----------------------------------------------------------------------
c  reads in all radial integrals over perturbation and elastic kernels
c  as well as the ellipticity correction factors
c  loads the common blocks
c-----------------------------------------------------------------------
      include 'parameter5.h'
c-----------------------------------------------------------------------
c declaration for radial integrals and ellipticity correction factor
c spheroidal as well as toroidal

      integer*4 jcom_h,ifanis_h,knot_h,lmin,lmax,imode
      integer*4 nnll(0:maxn,0:maxl)
      real*4 tref_h,rn_h
      real*8 int_ker(maxmodes,0:kmax1)
      real*8 int_cru(maxmodes)
      real*8 wa_t(maxmodes),gvelo(maxmodes)
      character*80 model_h

      common /integ_const_h/model_h,tref_h,rn_h,jcom_h,ifanis_h,imode
      common /integ_h/int_ker,int_cru,wa_t,gvelo,nnll
c-----------------------------------------------------------------------
c declaration for M84 coefficients

      integer*4 icrust,ismax
      complex*16 CA(-nord:nord,0:nord,0:kmax)
      complex*16 CC(-nord:nord,0:nord,0:kmax)
      complex*16 crustCC(-nord:nord,0:nord)

      common /coeff_m84/CA,CC,crustCC,icrust,ismax
c-----------------------------------------------------------------------
c declaration for SH8-12 coefficients

      complex*16 manW(-nord:nord,0:nord,0:kmax1)
      complex*16 manL(-nord:nord,0:nord,0:kmax2)
      complex*16 manU(-nord:nord,0:nord,0:kmax3)
      complex*16 cruC(-nord:nord,0:nord)

      common /coeff_sh8/manW,manL,manU,cruC
c-----------------------------------------------------------------------
c declaration for local variables
      
      integer*4 nn,ll
      real*8 w,q,gv,walpha,integ_kernel(0:kmax),integ_bound,integ_crust
      real*8 integ_k1(0:kmax1)
      logical lsh8w
      character*80 int_file,hname

c-----------------------------------------------------------------------
c sets all arrays to zero 
c reads input file 
c-----------------------------------------------------------------------
      model_h(1:80) = hname(1:80)
      do ii=1,maxmodes
         do kk = 0,kmax1
            int_ker(ii,kk)=0.d0
         enddo
         int_cru(ii)=0.d0
         wa_t(ii)=0.d0
      enddo

      print*,' Reading radial integral files for  ',model_h(1:4)

      if (model_h(1:4).eq.'M84A') icrust=0
      if (model_h(1:4).eq.'M84C') icrust=1
      if (model_h(1:4).eq.'SH8W') icrust=1
      if (model_h(1:4).eq.'SH8U') icrust=1
      if (model_h(1:4).eq.'SH12') icrust=1
      
      if (model_h(1:3).eq.'M84') then

        open(1,file=int_file,form='unformatted',status='old')
        read(1) jcom_h,lmin,lmax,ismax,knot_h,ifanis_h,
     &                 tref_h,rn_h
        do kk =0,kmax
           do iss = 0,ismax
             read(1) (CA(itt,iss,kk),CC(itt,iss,kk),itt=-iss,iss)
           enddo
        enddo
        do iss = 0,ismax
           read(1) (crustCC(itt,iss),itt=-iss,iss)
        enddo

        imode = 0
        do while (imode.lt.maxmodes)
          read(1,end=5) nn,ll,w,q,gv,
     &               (integ_kernel(kk),kk=0,kmax),
     &                integ_crust,integ_bound,walpha
          imode = imode + 1
          nnll(nn,ll) = imode
          do kk=0,kmax
            int_ker(imode,kk) = integ_kernel(kk)
          enddo
          int_cru(imode) = integ_crust-integ_bound
          wa_t(imode) = walpha
          gvelo(imode) = gv
        enddo 

      elseif (model_h(1:2).eq.'SH') then

        open(1,file=int_file,form='unformatted',status='old')
        read(1) jcom_h,lmin,lmax,ismax,knot_h,ifanis_h,
     &                 tref_h,rn_h,lsh8w
        if (lsh8w) then
cc...whole mantle model
          if (model_h(1:4).ne.'SH8W'.and.model_h(1:4).ne.'SH12') then
            print*,' WARNING : the integral file is not Whole mantle'
            print*,"           let's stop this non-sense !!"
            close(1)
            stop
          endif
          do kk =0,kmax1
             do iss = 0,ismax
               read(1) (manW(itt,iss,kk),itt=-iss,iss)
             enddo
          enddo

        else
cc...discontinuous model
          if (model_h(1:4).ne.'SH8U') then
            print*,' WARNING : the integral file is for up/low mantle'
            print*,"           let's stop this non-sense !!"
            close(1)
            stop
          endif
          do kk =0,kmax2
             do iss = 0,ismax
               read(1) (manL(itt,iss,kk),itt=-iss,iss)
             enddo
          enddo
          do kk =0,kmax3
             do iss = 0,ismax
               read(1) (manU(itt,iss,kk),itt=-iss,iss)
             enddo
          enddo
        endif

cc...both models
        do iss = 0,ismax
           read(1) (cruC(itt,iss),itt=-iss,iss)
        enddo

        imode = 0
        do while (imode.lt.maxmodes)
          read(1,end=5) nn,ll,w,q,gv,
     &               (integ_k1(kk),kk=0,kmax1),
     &                integ_crust,integ_bound,walpha
cjbg
c          if (nn.eq.0) then 
c           print*,ll,integ_crust,integ_bound,walpha
c          endif
cjbg
          imode = imode + 1
          nnll(nn,ll) = imode
          do kk=0,kmax1
            int_ker(imode,kk) = integ_k1(kk)
          enddo
          int_cru(imode) = integ_crust-integ_bound
          wa_t(imode) = walpha
          gvelo(imode) = gv
        enddo 

      endif
                     
   5  close(1)
      print*,' Order in spherical harmonics : ',ismax

      return
      end
