
      logical function lexist (filename)

c  Function to determine whether "filename" exists.


      character*(*) filename
      integer*4 nblen

      if (nblen(filename) .eq. 0) then
         lexist = .false.
      else
         inquire (file = filename, exist = lexist)
         if (.not. lexist) then
            print*, filename (1:nblen (filename)), ' does not exist.'
         else
         end if
      end if

      return

      end
