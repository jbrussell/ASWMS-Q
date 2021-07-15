      character*1 qual,atemp 
      character*2 atemp2
      character*3 atemp3

      it = 1
      write(atemp3,'(i3)') it
      open(4,file='junk.'//atemp3)
      close(4)

      it = 10
      write(atemp3,'(i3)') it
      open(4,file='junk.'//atemp3)
      close(4)

      it = 100
      write(atemp3,'(i3)') it
      open(4,file='junk.'//atemp3)
      close(4)

      stop
      end
