! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
program main

  ! MAIN is a program to solve several problems as defined in user-
  ! configured input files
  ! 
  ! Modified: 6 August 2010

  use constants, only : pin,pr
  use problem, only : problem_type,init_problem,destroy_problem
  use problem_routines, only : solve_problem
  use io, only : new_io_unit,error,message
  use test_utilities, only : test_newton,test_interp,test_update,test_gauss
  use utilities, only : convert_time
  use mpi_routines, only : start_mpi,finish_mpi,is_master,master,nprocs,clock_split

  implicit none

  ! Parameters:
  ! PB = problem variables
  ! INPUT = file unit for problem list file
  ! NAME = problem name
  ! STAT = I/O error flag
  ! IERR = MPI error flag

  type(problem_type) :: pb
  integer(pin) :: input,stat,hr,mn,ninfo
  character(64) :: name,input_file,temp
  character(256) :: str
  real(pr) :: sc,start_time,end_time,total_time

  namelist /problem_list/ name,ninfo

  ! get problem name

  call get_command_argument(1,input_file,status=stat)
  if (stat/=0) input_file = 'default.in'

  ! start MPI

  call start_mpi
  call clock_split(current_time=start_time)

  !---BEGIN INITIALIZATION---

  input = new_io_unit()
     
  if (is_master) then
     
     ! write MPI status
     
     if (nprocs/=1) then
        write(temp,'(a,i3,a,i3,a)') 'Start MPI: number of processes = ',nprocs
        call message(temp)
     end if
     
  end if

  ! open input file
  
  open(input,file=input_file,status='old',iostat=stat)
  
  if (stat/=0) &
       call error("Error opening file '" // trim(input_file) // "'",'main')
  
  !---END INITIALIZATION---

  ! defaults
  
  name = 'default'
  ninfo = 1
  
  ! read in problem parameters
  
  read(input,nml=problem_list,iostat=stat)
  if (stat>0) call error('Error in problem_list','main')
  
  if (is_master) call message('Starting problem: ' // trim(name))

  ! solve problem

  select case(name)
     
  case default
     
     ! create derived type pb
        
     call init_problem(pb,name)
     
     ! initialization complete
     
     call clock_split(current_time=end_time)
     total_time = end_time-start_time
     start_time = end_time
     if (is_master) then
        write(str,'(a,f0.6,a)') 'Initialization complete (',total_time,' s)'
        call message(str)
     end if

     ! solve the problem
     
     call solve_problem(pb,input,start_time,ninfo)
     
     ! collect timing information
     
     call clock_split(current_time=end_time)
     total_time = end_time-start_time
     
     ! display timing information
     
     if (is_master) then
        call convert_time(total_time,hr,mn,sc)
        write(str,'(a,f0.8,a,i0,a,f0.4,a)') &
             'total: ',total_time,' s for ',pb%mdl%n,' time steps (', &
             total_time/dble(pb%mdl%n),' s/step)'
        call message(str)
        write(str,'(a,i0,a,i0,a,f0.8,a)') &
             '= ',hr,' hr ',mn,' min ',sc, ' s'
        call message(str)
        call message('Finished problem: ' // trim(adjustl(name)))
     end if
     
     ! destroy the problem (deallocate memory in pointers, etc.)
     
     call destroy_problem(pb)
     
  case ('test_newton')
     
     ! test Newton's method routines
     
     call test_newton
     
  case ('test_interp')
     
     ! test interpolation routines
     
     call test_interp
     
  case ('test_update')
     
     ! test integration update routines
     
     call test_update
     
  case ('test_gauss')
     
     ! test solving linear system with Gaussian elimination
     
     call test_gauss
     
  case('stop')
     
     stop
     
  end select

  close(input)

  ! finish MPI
  
  call finish_mpi
  
  stop

end program main
