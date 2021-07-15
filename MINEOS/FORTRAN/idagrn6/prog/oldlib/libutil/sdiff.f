        double precision function sdiff (dateo,daten)
c
c returns the difference in the input dates
c J. Borsenberger 17/4/89 modified 13/10/89
c
       implicit double precision (a-h,o-z)
       integer*4 dateo(6),daten(6)
c
       jour=daten(2)-dateo(2)
       if (daten(1).gt.dateo(1)) then
         do iannee=dateo(1),daten(1)-1
           if (mod(iannee,4).eq.0 .and. mod(iannee,2000).ne.0) then
             jour=jour+366
           else
             jour=jour+365
           endif
         enddo
       else if (daten(1).lt.dateo(1)) then
         do iannee=daten(1),dateo(1)-1
           if (mod(iannee,4).eq.0 .and. mod(iannee,2000).ne.0) then
             jour=jour-366
           else
             jour=jour-365
           endif
         enddo
       endif
       sdiff= 1.d-3*(daten(6)-dateo(6))+ daten(5)-dateo(5)+ 
     1          6.d1*(daten(4)-dateo(4) + 6.d1* (daten(3)-dateo(3)+
     2          24.d0* jour ))
       return
       end 
