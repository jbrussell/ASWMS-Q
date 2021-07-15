	subroutine sdcoht(time, deltatime)
c version pour deltatime en secondes modified 13/10/89
      integer*4 time(6), deltatime, modul(6)
      data modul / 0, 365, 24, 2*60, 1000 /
	if (deltatime.lt.0.) then
		write (*,*) '** subr addcoht neg argument',deltatime,time
		stop
	endif
	if (mod(time(1),4).eq.0) then
		modul(2)=366
	else
		modul(2)=365
	endif
	time(2)=time(2)-1
      j1 = deltatime
      do j = 5, 2, -1
		j2 = time(j) + j1
		time(j) = mod(j2,modul(j))
		if (j2 .eq. time(j)) then
			go to 80
		else
			j1 = j2 / modul(j)
		end if
      enddo
c supposes less than a year span, to be changed (be carefull for bissextile)
      time(1) = time(1) + 1
80	time(2)=time(2)+1
      return 
      end
