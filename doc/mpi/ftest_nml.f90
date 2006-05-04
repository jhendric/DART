! very simple fortran program which reads in an external namelist file.
! if successful, will print a message and exit.

program ftest_nml

integer :: iunit, errcode

integer :: array_size, array_data(10)
namelist / ftest / array_size, array_data

   print *, "program start"
  
   array_size = -1
   array_data = -99

   iunit = 11
   open (iunit, name="ftest_nml.nml", iostat=errcode)
   if (errcode /= 0) then
       print *, "cannot open namelist file, error = ", errcode
       stop
   endif

   read (iunit, nml=ftest, iostat=errcode)
   if (errcode /= 0) then
       print *, "cannot read namelist file, error = ", errcode
       stop
   endif

   close (iunit, iostat=errcode)

   print *, "array size should be 10, value is ", array_size
   print *, "array contents should be 1-10, values are: ", array_data

   print *, "program end"

end program ftest_nml

