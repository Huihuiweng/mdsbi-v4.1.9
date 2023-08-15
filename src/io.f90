! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module io

  ! IO contains variables and routines for input/output
  ! 
  ! Modified: 19 July 2010
  
  use constants, only : pr,pin

  implicit none

  ! Global Parameters:
  ! STDOUT = unit of standard out
  ! STDERROR = unit of standard error
  ! COLLECTIVE = collective reads/writes

  integer(pin),save :: stdout=6,stderr=0
  logical,parameter :: collective=.true.

  type :: file_distributed
     character(64) :: name
     integer :: fh,array,pr,MPI_REAL_PR
  end type file_distributed

  interface write_matlab
     ! WRITE_MATLAB overloads the subroutine name for writing
     ! various types of data in matlab format
     module procedure write_matlab_character
     module procedure write_matlab_logical
     module procedure write_matlab_integer
     module procedure write_matlab_integer8
     module procedure write_matlab_real
  end interface

  interface write_file_distributed
     module procedure write_file_distributed_1d,write_file_distributed_2d, &
          write_file_distributed_3d
  end interface

  interface read_file_distributed
     module procedure read_file_distributed_1d,read_file_distributed_2d, &
          read_file_distributed_3d
  end interface

contains


  subroutine flush_io(unit)
    ! FLUSH_IO flushes the file unit (forces data to be written now, rather
    ! than queued up--whether or not this is done automatically depends
    ! on your compiler).  Note that this is not Fortran 95 standard (but is Fortran 2003)
    ! so you may need to comment it if you get a compiler error.  So you don't
    ! have to track this down in many places, I've centralized the call to this
    ! subroutine.
    !
    ! Modified: 29 May 2007

    implicit none

    ! I/O Parameters:
    ! UNIT = file unit to be flushed

    integer(pin),intent(in) :: unit

    ! uncomment one of the following lines for immediate data output,
    ! comment the following lines for maximum speed or if your compiler complains

    !call flush(unit) ! works with most compilers
    flush(unit) ! Fortran 2003 standard

  end subroutine flush_io


  function new_io_unit() result(next)
    ! NEW_IO_UNIT returns the next available file unit
    ! 
    ! Modified: 20 February 2005

    implicit none

    ! I/O Parameters:
    ! NEXT = unit number of new file

    integer(pin) :: next

    ! Internal Parameters:
    ! MIN = minimum unit number
    ! MAX = maximum unit number
    ! LAST = unit number that was last returned

    integer(pin),parameter :: min=10,max=999
    integer(pin),save :: last=10
    logical :: open,top=.false.

    ! start at last opened unit and move up, checking if file unit
    ! is in use; if not, then return that number

    next = min-1
    if (last>0) then
       next = last+1
       inquire(unit=next,opened=open)
       if (.not.open) last = next
       return
    else
       do
          next = next+1
          inquire(unit=next,opened=open)
          if (.not.open) then
             last = next
             exit
          end if
          if (next==max) then
             if (top) call error('No available file unit','new_io_unit')
             top = .true.
             next = min-1
          end if
       end do
    end if

  end function new_io_unit

  
  subroutine message(str,routine)
    ! MESSAGE writes informational message

    use mpi_routines, only : my_rank

    implicit none

    ! STR = string to be written to standard out
    ! ROUTINE =  subroutine name in which message originated

    character(*),intent(in) :: str
    character(*),intent(in),optional :: routine
   
    character(6) :: id

    ! write message and subroutine name (if given)

    write(id,'(i4)') my_rank
    write(id,'(a)') '(' // trim(adjustl(id)) // ') '

    if (present(routine)) write(stdout,'(a)') id // &
         'Message from subroutine: ' // trim(adjustl(routine))
    write(stdout,'(a)') id // trim(adjustl(str))

  end subroutine message


  subroutine warning(str,routine)
    ! WARNING writes error message, but does not terminate program

    use mpi_routines, only : my_rank

    implicit none

    ! STR = string to be written to standard error
    ! ROUTINE =  subroutine name in which error occurred

    character(*),intent(in) :: str
    character(*),intent(in),optional :: routine

    character(6) :: id

    ! write error message and subroutine name (if given)
    
    write(id,'(i4)') my_rank
    write(id,'(a)') '(' // trim(adjustl(id)) // ') '

    if (present(routine)) write(stderr,'(a)') id // &
         'Warning in subroutine: ' // trim(adjustl(routine))
    write(stderr,'(a)') id // trim(adjustl(str))

  end subroutine warning


  subroutine error(str,routine)
    ! ERROR writes error message and terminates program

    use mpi_routines, only : my_rank
    use mpi

    implicit none

    ! STR = string to be written to standard error
    ! ROUTINE =  subroutine name in which error occurred

    character(*),intent(in) :: str
    character(*),intent(in),optional :: routine
   
    ! IERR = MPI error flag

    integer :: ierr
    character(6) :: id

    ! write error message and subroutine name (if given)
    
    write(id,'(i4)') my_rank
    write(id,'(a)') '(' // trim(adjustl(id)) // ') '

    if (present(routine)) write(stderr,'(a)') id // &
         'Error in subroutine: ' // trim(adjustl(routine))
    write(stderr,'(a)') id // trim(adjustl(str))
    write(stderr,'(a)') id // 'Terminating program'

    ! terminate program

    call MPI_Abort(MPI_COMM_WORLD,0,ierr)

  end subroutine error


  subroutine open_array_data(name,unit,ascii)
    ! OPEN_ARRAY_DATA opens file containing array data
    ! 
    ! Modified: 25 October 2007

    implicit none

    ! I/O Parameters:
    ! NAME = file name
    ! UNIT = file unit
    ! ASCII = input file format
    !      T = ASCII text file
    !      F = binary file

    character(*),intent(in) :: name
    integer(pin),intent(out) :: unit
    logical,intent(in) :: ascii

    ! Internal Parameters:
    ! STAT = I/O error flag
    ! FORM = input file format

    integer(pin) :: stat
    character(64) :: form

    unit = new_io_unit()

    if (ascii) then
       form = 'formatted'
    else
       form = 'unformatted'
    end if

    ! open data file, call error routine if needed
    
    open(unit,file=name,form=form,iostat=stat,status='old')
    if (stat/=0) call error("Error opening file '" // trim(adjustl(name)) // "'",'open_array_data')
       
  end subroutine open_array_data


  subroutine close_array_data(unit)
    ! CLOSE_ARRAY_DATA closes file containing array data
    ! 
    ! Modified: 15 March 2005

    implicit none

    ! I/O Parameters:
    integer(pin),intent(in) :: unit

    ! Internal Parameters:
    ! STAT = I/O error flag

    integer(pin) :: stat

    close(unit,iostat=stat)

    ! error trap closing

    if (stat/=0) call error('Error closing array data file','close_array_data')

  end subroutine close_array_data


  subroutine read_array_data(unit,data,ascii,arrname)
    ! READ_ARRAY_DATA reads array data from file
    ! 
    ! Modified: 26 February 2006

    implicit none

    ! I/O Parameters:
    ! UNIT = file unit
    ! DATA = data array
    ! ASCII = input file format
    !      T = ASCII text file
    !      F = binary file
    ! ARRNAME = name of data array

    integer(pin),intent(in) :: unit
    real(pr),dimension(:,:),intent(inout) :: data
    logical,intent(in) :: ascii
    character(*),intent(in),optional :: arrname

    ! Internal Parameters:
    ! STAT = I/O error flag
    ! NAME = file name
    ! STR = error string

    integer(pin) :: stat
    character(64) :: name
    character(128) :: str

    ! read data array

    if (ascii) then
       read(unit,*,iostat=stat) data
    else
       read(unit,iostat=stat) data
    end if

    ! call error routine if needed
    if (stat/=0) then
       inquire(unit,name=name)
       if (present(arrname)) then
          str = "Error reading array '" // trim(adjustl(arrname)) // &
               "' from file '" // trim(adjustl(name)) // "'"
       else
          str = "Error reading array from file '" // trim(adjustl(name)) // "'"
       end if
       call error(str,'read_array_data')
    end if
    
  end subroutine read_array_data


  subroutine openfile(fname,n,unit)
    ! Modified: 17 June 2010

    implicit none

    character(*),intent(in) :: fname
    integer(pin),intent(in) :: n
    integer(pin),intent(out) :: unit

    integer(pin) :: ierr

    unit = new_io_unit()

    ! record length calculated assuming double precision

    open(unit,file=fname,access='direct',recl=8*n,iostat=ierr,status='old')
    if (ierr/=0) call error("Error opening file '" // trim(adjustl(fname)) // "'",'openfile')

  end subroutine openfile


  subroutine read0dfield(f,unit,n)
    ! Modified: 17 June 2010

    implicit none
    
    integer,intent(in) :: unit
    real,intent(out) :: f
    integer,intent(in),optional :: n

    if (present(n)) then
       read(unit,rec=n) f
    else
       read(unit,rec=1) f
    end if

  end subroutine read0dfield


  subroutine read1dfield(f,unit,n)
    ! Modified: 17 June 2010
    
    implicit none
    
    integer,intent(in) :: unit
    real,dimension(:),intent(out) :: f
    integer,intent(in),optional :: n

    if (present(n)) then
       read(unit,rec=n) f
    else
       read(unit,rec=1) f
    end if

  end subroutine read1dfield


  subroutine read2dfield(f,unit,n)
    ! Modified: 17 June 2010
    
    implicit none
    
    integer,intent(in) :: unit
    real,dimension(:,:),intent(out) :: f
    integer,intent(in),optional :: n

    if (present(n)) then
       read(unit,rec=n) f
    else
       read(unit,rec=1) f
    end if

  end subroutine read2dfield


  subroutine read3dfield(f,unit,n)
    ! Modified: 17 June 2010
    
    implicit none
    
    integer,intent(in) :: unit
    real,dimension(:,:,:),intent(out) :: f
    integer,intent(in),optional :: n

    if (present(n)) then
       read(unit,rec=n) f
    else
       read(unit,rec=1) f
    end if

  end subroutine read3dfield


  subroutine write_matlab_character(nunit,name,str,struct)
    ! WRITE_MATLAB_CHARACTER writes character string to matlab file
    ! 
    ! Modified: 20 February 2005

    implicit none

    ! I/O Parameters:
    ! NUNIT = unit number of output file
    ! NAME = name of variable
    ! STR = value of variable
    ! STRUCT = name of structure in which to write variable

    integer(pin),intent(in) :: nunit
    character(*),intent(in) :: name,str
    character(*),intent(in),optional :: struct

    ! Internal Parameters:
    ! STAT = I/O error flag

    integer(pin) :: stat

    ! write character string to file

    if (present(struct)) then
       write(nunit,'(a)',iostat=stat) &
            trim(adjustl(struct)) // '.' // trim(adjustl(name)) // " = '" // trim(adjustl(str)) // "';"
    else
       write(nunit,'(a)',iostat=stat) &
            trim(adjustl(name)) // " = '" // trim(adjustl(str)) // "';"
    end if

    if (stat/=0) call error('Error writing ' // trim(adjustl(name)) // ' to matlab file','write_matlab_character')

  end subroutine write_matlab_character


  subroutine write_matlab_integer(nunit,name,data,struct)
    ! WRITE_MATLAB_INTEGER writes integer to matlab file
    ! 
    ! Modified: 29 August 2006

    implicit none

    ! I/O Parameters:
    ! NUNIT = unit number of output file
    ! NAME = name of variable
    ! DATA = value of variable
    ! STRUCT = name of structure in which to write variable

    integer(pin),intent(in) :: nunit,data
    character(*),intent(in) :: name
    character(*),intent(in),optional :: struct

    ! Internal Parameters:
    ! STAT = I/O error flag

    integer(pin) :: stat
    
    ! write integer to file

    if (present(struct)) then
       write(nunit,'(a,i12,a)',iostat=stat) &
            trim(adjustl(struct)) // '.' // trim(adjustl(name)) // ' = ',data,';'
    else
       write(nunit,'(a,i12,a)',iostat=stat) &
            trim(adjustl(name)) // ' = ',data,';'
    end if
    
    if (stat/=0) call error('Error writing ' // trim(adjustl(name)) // ' to matlab file','write_matlab_character')

  end subroutine write_matlab_integer


  subroutine write_matlab_integer8(nunit,name,data,struct)
    ! WRITE_MATLAB_INTEGER8 writes 8-byte integer to matlab file
    ! 
    ! Modified: 29 August 2006

    use constants, only : int8
    implicit none

    ! I/O Parameters:
    ! NUNIT = unit number of output file
    ! NAME = name of variable
    ! DATA = value of variable
    ! STRUCT = name of structure in which to write variable

    integer(pin),intent(in) :: nunit
    integer(int8),intent(in) :: data
    character(*),intent(in) :: name
    character(*),intent(in),optional :: struct

    ! Internal Parameters:
    ! STAT = I/O error flag

    integer(pin) :: stat
    
    ! write integer to file

    if (present(struct)) then
       write(nunit,'(a,i12,a)',iostat=stat) &
            trim(adjustl(struct)) // '.' // trim(adjustl(name)) // ' = ',data,';'
    else
       write(nunit,'(a,i12,a)',iostat=stat) &
            trim(adjustl(name)) // ' = ',data,';'
    end if
    
    if (stat/=0) call error('Error writing ' // trim(adjustl(name)) // ' to matlab file','write_matlab_character')

  end subroutine write_matlab_integer8


  subroutine write_matlab_real(nunit,name,data,struct)
    ! WRITE_MATLAB_REAL writes real to matlab file
    ! 
    ! Modified: 20 February 2005

    implicit none

    ! I/O Parameters:
    ! NUNIT = unit number of output file
    ! NAME = name of variable
    ! DATA = value of variable
    ! STRUCT = name of structure in which to write variable

    integer(pin),intent(in) :: nunit
    character(*),intent(in) :: name
    real(pr),intent(in) :: data
    character(*),intent(in),optional :: struct

    ! Internal Parameters:
    ! STAT = I/O error flag

    integer(pin) :: stat

    ! write real to file

    if (present(struct)) then
       write(nunit,'(a,e14.6,a)',iostat=stat) &
            trim(adjustl(struct)) // '.' // trim(adjustl(name)) // ' = ',data,';'
    else
       write(nunit,'(a,e14.6,a)',iostat=stat) &
            trim(adjustl(name)) // ' = ',data,';'
    end if

    if (stat/=0) call error('Error writing ' // trim(adjustl(name)) // ' to matlab file','write_matlab_character')

  end subroutine write_matlab_real


  subroutine write_matlab_logical(nunit,name,data,struct)
    ! WRITE_MATLAB_LOGICAL writes logical to matlab file
    ! 
    ! Modified: 16 June 2007

    implicit none

    ! I/O Parameters:
    ! NUNIT = unit number of output file
    ! NAME = name of variable
    ! DATA = value of variable
    ! STRUCT = name of structure in which to write variable

    integer(pin),intent(in) :: nunit
    character(*),intent(in) :: name
    logical,intent(in) :: data
    character(*),intent(in),optional :: struct

    ! Internal Parameters:
    ! STR = string to hold 'T' or 'F'

    character(1) :: str

    ! assign logical value to string

    if (data) then
       str = 'T'
    else
       str = 'F'
    end if

    ! write string to matlab file

    if (present(struct)) then
       call write_matlab_character(nunit,name,str,struct)
    else
       call write_matlab_character(nunit,name,str)
    end if

  end subroutine write_matlab_logical


  subroutine save_data(nunit,name,data,nx,ny)
    ! SAVE_DATA writes a 2D array to a file in ascii format
    ! 
    ! Modified: 3 November 2005

    implicit none
    
    ! I/O Parameters:
    ! NUNIT = unit number of output file
    ! NAME = filename (to which .dat will be appended)
    ! DATA = data array
    ! NX = number of points in x direction
    ! NY = number of points in y direction

    integer,intent(in) :: nunit,nx,ny
    character(*),intent(in) :: name
    real(kind=pr),dimension(nx,ny),intent(in) :: data

    ! Internal Parameters:
    ! FILENAME = filename
    ! FMT = ascii format
    ! I = index in x direction
    ! J = index in y direction

    character(64) :: filename,fmt
    integer :: j

    filename = trim(adjustl(name)) // '.dat'

    open(nunit,file=filename,form='formatted')
    
    write(fmt,'(i12)') nx
    fmt = '(' // trim(adjustl(fmt)) // 'e16.7)'
    
    do j=1,ny
       write(nunit,trim(fmt)) data(:,j)
    end do
    
    close(nunit)
    
  end subroutine save_data


  subroutine open_file_distributed(fid,name,operation,comm,array,precision)

    use mpi_routines, only : MPI_REAL_PR,MPI_REAL_PSAV,pr,psav
    use mpi

    implicit none

    type(file_distributed),intent(out) :: fid
    character(*),intent(in) :: name,operation
    integer,intent(in) :: comm,array,precision

    integer :: ierr
    integer(MPI_OFFSET_KIND),parameter :: zero = 0
    character(128) :: str

    ! file name

    fid%name = name

    ! precision

    if (precision==pr) then
       fid%pr = pr
       fid%MPI_REAL_PR = MPI_REAL_PR
    else
       fid%pr = psav
       fid%MPI_REAL_PR = MPI_REAL_PSAV
    end if

    ! open file
    
    select case(operation)
    case('read')
       call MPI_File_open(comm,fid%name,MPI_MODE_RDONLY, &
            MPI_INFO_NULL,fid%fh,ierr)
    case('write','append')
       call MPI_File_open(comm,fid%name,MPI_MODE_CREATE+MPI_MODE_WRONLY, &
            MPI_INFO_NULL,fid%fh,ierr)
    case default
       call error('Invalid file operation','open_file_distributed') 
    end select
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_open (on file ' // trim(name) // ')',ierr,str)
       call error(str,'open_file_distributed') 
    end if

    ! if writing, delete file (set size to zero)

    if (operation=='write') then
       call MPI_File_set_size(fid%fh,zero,ierr)
       if (ierr/=MPI_SUCCESS) then
          call io_error_message('MPI_File_set_size',ierr,str)
          call error(str,'open_file_distributed') 
       end if
    end if

    ! set view

    call MPI_File_set_view(fid%fh,zero,fid%MPI_REAL_PR,array,'native',MPI_INFO_NULL,ierr)
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_set_view',ierr,str)
       call error(str,'open_file_distributed') 
    end if

  end subroutine open_file_distributed

  
  subroutine write_file_distributed_1d(fid,data)

    use mpi_routines, only : pr,psav
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real(pr),dimension(:),intent(in) :: data

    integer :: ierr
    character(64) :: str

    if (fid%pr==pr) then
       if (collective) then
          call MPI_File_write_all(fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    else
       if (collective) then
          call MPI_File_write_all(fid%fh,real(data,psav),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,real(data,psav),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    end if

    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_write_all',ierr,str)
       call error(str,'write_file_distributed_1d')
    end if

  end subroutine write_file_distributed_1d


  subroutine write_file_distributed_2d(fid,data)

    use mpi_routines, only : pr,psav
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real(pr),dimension(:,:),intent(in) :: data

    integer :: ierr
    character(64) :: str

    if (fid%pr==pr) then
       if (collective) then
          call MPI_File_write_all(fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    else
       if (collective) then
          call MPI_File_write_all(fid%fh,real(data,psav),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,real(data,psav),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    end if

    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_write_all',ierr,str)
       call error(str,'write_file_distributed_2d')
    end if

  end subroutine write_file_distributed_2d


  subroutine write_file_distributed_3d(fid,data)

    use mpi_routines, only : pr,psav
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real(pr),dimension(:,:,:),intent(in) :: data

    integer :: ierr
    character(64) :: str

    if (fid%pr==pr) then
       if (collective) then
          call MPI_File_write_all(fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,data,size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    else
       if (collective) then
          call MPI_File_write_all(fid%fh,real(data,psav),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       else
          call MPI_File_write    (fid%fh,real(data,psav),size(data), &
               fid%MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
       end if
    end if

    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_write_all',ierr,str)
       call error(str,'write_file_distributed_3d')
    end if

  end subroutine write_file_distributed_3d


  subroutine read_file_distributed_1d(fid,data)

    use mpi_routines, only : MPI_REAL_PR
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real(pr),dimension(:),intent(out) :: data

    integer :: ierr
    character(64) :: str

    if (collective) then
       call MPI_File_read_all(fid%fh,data,size(data),MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
    else
       call MPI_File_read    (fid%fh,data,size(data),MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_read_all',ierr,str)
       call error(str,'read_file_distributed_1d')
    end if

  end subroutine read_file_distributed_1d


  subroutine read_file_distributed_2d(fid,data)

    use mpi_routines, only : MPI_REAL_PR
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real(pr),dimension(:,:),intent(out) :: data

    integer :: ierr
    character(64) :: str

    if (collective) then
       call MPI_File_read_all(fid%fh,data,size(data),MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
    else
       call MPI_File_read    (fid%fh,data,size(data),MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_read_all',ierr,str)
       call error(str,'read_file_distributed_2d')
    end if

  end subroutine read_file_distributed_2d


  subroutine read_file_distributed_3d(fid,data)

    use mpi_routines, only : MPI_REAL_PR
    use mpi

    implicit none
    
    type(file_distributed),intent(in) :: fid
    real(pr),dimension(:,:,:),intent(out) :: data

    integer :: ierr
    character(64) :: str

    if (collective) then
       call MPI_File_read_all(fid%fh,data,size(data),MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
    else
       call MPI_File_read    (fid%fh,data,size(data),MPI_REAL_PR,MPI_STATUS_IGNORE,ierr)
    end if
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_read_all',ierr,str)
       call error(str,'read_file_distributed_3d')
    end if

  end subroutine read_file_distributed_3d


  subroutine close_file_distributed(fid)

    use mpi

    implicit none

    type(file_distributed),intent(inout) :: fid

    integer :: ierr
    character(64) :: str

    call MPI_File_close(fid%fh,ierr)
    if (ierr/=MPI_SUCCESS) then
       call io_error_message('MPI_File_close',ierr,str)
       call error(str,'close_file_distributed')
    end if

  end subroutine close_file_distributed


  subroutine io_error_message(routine,ierr,str)

    use mpi

    implicit none

    character(*),intent(in) :: routine
    integer,intent(in) :: ierr
    character(*),intent(out) :: str

    !select case(ierr)
    !case(MPI_ERR_IO)
    !   write(str,'(a)') 'Problem with ' // trim(routine) // ': MPI_ERR_IO'
    !case(MPI_ERR_NO_SPACE)
    !   write(str,'(a)') 'Problem with ' // trim(routine) // ': MPI_ERR_NO_SPACE'
    !case(MPI_ERR_NO_SUCH_FILE)
    !   write(str,'(a)') 'Problem with ' // trim(routine) // ': MPI_ERR_NO_SUCH_FILE'
    !case default
       write(str,'(a,i6)') 'Problem with ' // trim(routine) // ': ierr=',ierr
    !end select

  end subroutine io_error_message


end module io
