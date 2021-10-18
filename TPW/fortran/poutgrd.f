c** print out the output grid
c** this can be used as a starting model
      character*70 filename

      print *,'enter beglat,endlat,dlat,beglon,endlon,dlon'
      read(*,*) beglat,endlat,dlat,beglon,endlon,dlon

      print *,'enter velocity'
      read(*,*) v0

      print *,'output file'
      read(*,'(a)') filename
      open(8,file=filename)
      write(8,'(6f10.3)') beglat,endlat,dlat,beglon,endlon,dlon
      write(8,'(f10.3)') v0

      nlat=(endlat-beglat)/dlat + 1.01
      nlon=(endlon-beglon)/dlon + 1.01
      do ilat = 1,nlat
      do ilon = 1,nlon
        prdlon=beglon+(ilon-1)*dlon
        prdlat=beglat+(ilat-1)*dlat
        write(8,'(3f10.3)') prdlon,prdlat,v0
      enddo
      enddo
      close(8)
      end
