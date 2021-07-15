*  FUNCTION TO SWAP BYTES IN 2 BYTE INTEGER
      integer*2 function iswap2(int0)

      integer*2 int0,itemp2
      character*1 char1(2),ctemp
      equivalence (itemp2,char1)

      itemp2=int0
      ctemp=char1(1)
      char1(1)=char1(2)
      char1(2)=ctemp
      iswap2=itemp2

      return
      end

*  FUNCTION TO SWAP BYTES IN 4 BYTE INTEGER
      integer*4 function iswap4(int0)

      integer*4 int0,itemp4
      character*1 char1(4),ctemp
      equivalence (itemp4,char1)

      itemp4=int0
      ctemp=char1(1)
      char1(1)=char1(4)
      char1(4)=ctemp
      ctemp=char1(2)
      char1(2)=char1(3)
      char1(3)=ctemp
      iswap4=itemp4

      return
      end

*  FUNCTION TO SWAP BYTES IN 4 BYTE REAL
      real*4 function rswap4(real0)

      real*4 real0,temp4
      character*1 char1(4),ctemp
      equivalence (temp4,char1)

      temp4=real0
      ctemp=char1(1)
      char1(1)=char1(4)
      char1(4)=ctemp
      ctemp=char1(2)
      char1(2)=char1(3)
      char1(3)=ctemp
      rswap4=temp4

      return
      end

*  FUNCTION TO SWAP BYTES IN 8 BYTE REAL
      real*8 function rswap8(real0)

      real*8 real0,temp8
      integer*4 iswap4,integ4(2),itemp4
      equivalence (temp8,integ4)

      temp8=real0
      integ4(1)=iswap4(integ4(1))
      integ4(2)=iswap4(integ4(2))
      itemp4=integ4(1)
      integ4(1)=integ4(2)
      integ4(2)=itemp4
      rswap8=temp8

      return
	end

