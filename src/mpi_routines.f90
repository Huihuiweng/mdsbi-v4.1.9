! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module mpi_routines

  use constants, only : pin,pr,psav

  implicit none

  save

  integer(pin),parameter ::  stdout=6, stderr=0
  integer(pin),parameter :: master=0
  logical :: is_master
  integer(pin) :: nprocs,my_rank,MPI_REAL_PR,MPI_CPLX_PR,MPI_REAL_PSAV
  real(pr) :: starttime,endtime,oldtime,newtime

  interface subarray
     module procedure subarray1d,subarray2d,subarray3d
  end interface
  

contains

  
  subroutine start_mpi

    use constants, only : real4,real8
    use mpi
    
    implicit none

    integer(pin) :: ierr
    
    call MPI_Init(ierr)

    select case(pr)
    case(real4)
       MPI_REAL_PR = MPI_REAL
       MPI_CPLX_PR = MPI_COMPLEX
    case(real8)
       MPI_REAL_PR = MPI_DOUBLE_PRECISION
       MPI_CPLX_PR = MPI_DOUBLE_COMPLEX
    end select
    
    select case(psav)
    case(real4)
       MPI_REAL_PSAV = MPI_REAL
    case(real8)
       MPI_REAL_PSAV = MPI_DOUBLE_PRECISION
    end select
    
    starttime = MPI_WTime()
    
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
    is_master = (my_rank==master)

  end subroutine start_mpi


  subroutine clock_split(split_time,current_time)

    use mpi

    real(pr),intent(out),optional :: split_time,current_time

    newtime = MPI_WTime()
    if (present(split_time)) split_time = newtime-oldtime
    if (present(current_time)) current_time = newtime
    if (is_master.and.(.not.present(split_time)).and.(.not.present(current_time))) &
         write (stdout,'(/,a,g20.10,a,/)') 'Split time ',newtime-oldtime,' s'
    oldtime = newtime

  end subroutine clock_split


  subroutine finish_mpi

    use mpi

    implicit none

    integer(pin) :: ierr

    endtime = MPI_WTime()

    if (is_master) then
       write (stdout,'(a)') ''
       write (stdout,'(a,g20.10,a)') 'Total MPI time ',endtime-starttime,' s'
    end if

    call MPI_Finalize(ierr)

  end subroutine finish_mpi


  function using_mpi()

    use constants, only : zero
    use mpi

    implicit none

    logical :: using_mpi

    if (MPI_WTime()==zero) then
       using_mpi = .false.
    else
       using_mpi = .true.
    end if

  end function using_mpi


  subroutine decompose1d(n,p,i,l,u,c)

    use mpi

    implicit none
    
    ! N = number of tasks
    ! P = number of processors
    ! I = rank of process (0<=i<=p-1)
    ! L = starting index
    ! U = ending index
    ! C = count of tasks for process i

    integer(pin),intent(in) :: n,p,i
    integer(pin),intent(out) :: l,u,c

    integer(pin) :: m,r,ierr

    if (n<p) then
       write(stderr,'(a)') &
            'Error in decompose1d: number of tasks must not be less than number of processes'
       call MPI_Abort(MPI_COMM_WORLD,0,ierr)
    end if

    m = n/p ! integer division ignores any fractional part
    r = mod(n,p) ! remainder

    ! assign p-r processes m tasks each, r processes m+1 tasks each

    ! case 0, r=0

    if (r==0) then
       c = m
       l = 1+i*m
       u = l+c-1
       return
    end if

    ! case 1, p and r even OR p and r odd, symmetric task distribution
    
    if (mod(p,2)==mod(r,2)) then
       if (i<(p-r)/2) then ! low rank (0:(p-r)/2-1), m tasks
          c = m
          l = 1+i*m
       elseif (i<(p+r)/2) then ! intermediate rank ((p-r)/2:(p+r)/2-1), m+1 tasks
          c = m+1
          l = 1+(p-r)/2*m+(i-(p-r)/2)*(m+1)
       else ! high rank ((p+r)/2:p-1), m tasks
          c = m
          l = 1+(p-r)/2*m+r*(m+1)+(i-(p+r)/2)*m
       end if
    end if

    ! case 2, p odd and r even, symmetric task distribution
    
    if (mod(p,2)/=0.and.mod(r,2)==0) then
       if (i<r/2) then ! low rank (0:r/2-1), m+1 tasks
          c = m+1
          l = 1+i*(m+1)
       elseif (i<p-r/2) then ! intermediate rank (r/2:p-r/2-1), m tasks
          c = m
          l = 1+(r/2)*(m+1)+(i-r/2)*m
       else ! high rank (p-r/2:p-1), m+1 tasks
          c = m+1
          l = 1+(r/2)*(m+1)+(p-r)*m+(i-(p-r/2))*(m+1)
       end if
    end if

    ! case 3, p even and r odd, asymmetric task distribution

    if (mod(p,2)==0.and.mod(r,2)/=0) then
       if (i<p-r) then ! low rank (0:p-r-1), m tasks
          c = m
          l = 1+i*m
       else ! high rank (p-r:p-1), m+1 tasks
          c = m+1
          l = 1+(p-r)*m+(i-p+r)*(m+1)
       end if
    end if

    u = l+c-1

  end subroutine decompose1d


  subroutine subarray1d(n,m,p,precision,array)

    use mpi

    implicit none

    integer,intent(in) :: n,m,p,precision
    integer,intent(out) :: array

    integer,parameter :: dim = 1
    integer,dimension(dim) :: sizes,subsizes,starts
    integer :: ierr

    sizes = n
    subsizes = p-m+1
    starts = m-1 ! subtract 1 since MPI starting index is 0

    call MPI_Type_create_subarray(dim,sizes,subsizes,starts, &
         MPI_ORDER_FORTRAN,precision,array,ierr)
    call MPI_Type_commit(array,ierr)

  end subroutine subarray1d


  subroutine subarray2d(nx,ny,mx,px,my,py,precision,array)

    use mpi

    implicit none

    integer,intent(in) :: nx,ny,mx,px,my,py,precision
    integer,intent(out) :: array

    integer,parameter :: dim = 2
    integer,dimension(dim) :: sizes,subsizes,starts
    integer :: ierr

    sizes = (/ nx,ny /)
    subsizes = (/ px-mx+1,py-my+1 /)
    starts = (/ mx,my /)-1 ! subtract 1 since MPI starting index is 0

    call MPI_Type_create_subarray(dim,sizes,subsizes,starts, &
         MPI_ORDER_FORTRAN,precision,array,ierr)
    call MPI_Type_commit(array,ierr)

  end subroutine subarray2d


  subroutine subarray3d(nx,ny,nz,mx,px,my,py,mz,pz,precision,array)

    use mpi

    implicit none

    integer,intent(in) :: nx,ny,nz,mx,px,my,py,mz,pz,precision
    integer,intent(out) :: array

    integer,parameter :: dim = 3
    integer,dimension(dim) :: sizes,subsizes,starts
    integer :: ierr

    sizes = (/ nx,ny,nz /)
    subsizes = (/ px-mx+1,py-my+1,pz-mz+1 /)
    starts = (/ mx,my,mz /)-1 ! subtract 1 since MPI starting index is 0

    call MPI_Type_create_subarray(dim,sizes,subsizes,starts, &
         MPI_ORDER_FORTRAN,precision,array,ierr)
    call MPI_Type_commit(array,ierr)

  end subroutine subarray3d


  subroutine new_communicator(in_new,comm_new,comm)

    use mpi

    implicit none

    logical,intent(in) :: in_new
    integer,intent(out) :: comm_new
    integer,intent(in),optional :: comm

    integer :: np,n,ierr,comm_old,group_old,group_new,nprocs_new,rank
    integer,dimension(:),allocatable :: ranks,ranks_new

    ! old communicator (defaults to MPI_COMM_WORLD)

    if (present(comm)) then
       comm_old = comm
    else
       comm_old = MPI_COMM_WORLD
    end if

    call MPI_Comm_size(comm_old,np,ierr)
    allocate(ranks(np),ranks_new(np))

    ! determine which processes will comprise new communicator

    if (in_new) then
       call MPI_Comm_rank(comm_old,rank,ierr)
    else
       rank = -1
    end if

    call MPI_Allgather(rank,1,MPI_INTEGER,ranks,1,MPI_INTEGER,comm_old,ierr)

    nprocs_new = 0
    do n = 1,np
       if (ranks(n)==-1) cycle
       nprocs_new = nprocs_new+1
       ranks_new(nprocs_new) = ranks(n)
    end do

    ! create new communicator

    call MPI_Comm_group(comm_old,group_old,ierr)
    call MPI_Group_incl(group_old,nprocs_new, &
         ranks_new(1:nprocs_new),group_new,ierr)
    call MPI_Comm_create(comm_old,group_new,comm_new,ierr)
    if (in_new) call MPI_Group_free(group_new,ierr)
    call MPI_Group_free(group_old,ierr)

    deallocate(ranks,ranks_new)

  end subroutine new_communicator
  

end module mpi_routines
