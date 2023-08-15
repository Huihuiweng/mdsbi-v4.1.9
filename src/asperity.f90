! MDSBI version 4.1.9 -- Multi-Dimensional Spectral Boundary Integral Code
module asperity

  ! ASPERITY contains derived types and routines for manipulating a 
  ! singly-linked list of asperities; each asperity is characterized by 
  ! several parameters like its size, type, etc. as well as some number 
  ! of additional variables (the so-called data)
  ! 
  ! Modified: 20 July 2010

  use constants, only : pr,pin

  implicit none

  ! ASPERITY_TYPE is a derived type containing asperity variables
  !
  ! Parameters:
  ! TYPE = asperity shape
  !      box = rectangle
  !      ellipse = ellipse (circle is special case)
  !      gaussian = gaussian
  !      fourier = Fourier mode
  ! INSIDE = add perturbation inside or outside region
  !      T = add inside
  !      F = add outside
  ! X0 = center of asperity in x direction
  ! Y0 = center of asperity in y direction
  ! LENGTH_X = length of asperity in x direction
  ! LENGTH_Y = length of asperity in y direction

  type :: asperity_type
     character(20) :: type
     logical :: inside
     real(pr) :: x0,y0,length_x,length_y
  end type asperity_type

  ! ASPERITY_NODE is a derived type for basic node in singly-linked list 
  ! (single asperity variables and pointer to next node)
  !
  ! Parameters:
  ! ASP = asperity variables
  ! DATA = variables containing amplitude of the perturbed 
  !  fields within the asperity region
  ! NEXT = pointer to next node

  type :: asperity_node
     private
     type(asperity_type) :: asp
     real(pr),dimension(:),allocatable :: data
     type(asperity_node),pointer :: next=>null()
  end type asperity_node

  ! ASPERITY_LIST is a derived type containing singly-linked list of asperities
  !
  ! Parameters:
  ! HEAD = pointer to first node in list

  type :: asperity_list
     private
     type(asperity_node),pointer :: head=>null()
  end type asperity_list


contains


  subroutine create_asperity_node(list,asp,data)
    ! CREATE_ASPERITY_NODE adds node to list
    ! 
    ! Modified: 26 November 2006

    implicit none

    ! I/O Parameters:
    ! LIST = asperity list
    ! ASP = asperity variables
    ! DATA = perturbed fields amplitude data
    
    type(asperity_list),intent(inout) :: list
    type(asperity_type),intent(in) :: asp
    real(pr),dimension(:),intent(in),optional :: data

    ! Internal Parameters:
    ! CURRENT = pointer to current node in list

    type(asperity_node),pointer :: current=>null()

    ! move to end of list, create new node, point current to new 
    ! node, nullify current%next

    if (.not.associated(list%head)) then
       ! allocate memory for first node in list
       allocate(list%head)
       ! nullify pointer to next node
       nullify(list%head%next)
       ! point current to head
       current => list%head
    else
       ! start at head of list
       current => list%head
       ! move down list until empty node is reached
       do while(associated(current%next))
          ! end of list not found, move down one node
          current => current%next
       end do
       ! end of list found, add new node
       ! allocate memory for new node
       allocate(current%next)
       ! point current to new node
       current => current%next
       ! nullify pointer to next node
       nullify(current%next)
    end if

    ! store data in new node (pointed to by current)

    current%asp = asp
    if (present(data)) then
       allocate(current%data(size(data)))
       current%data = data
    end if

  end subroutine create_asperity_node


  subroutine destroy_asperity_node(asp)
    ! DESTROY_ASPERITY_NODE destroys asperity node
    ! 
    ! Modified: 26 November 2006

    implicit none

    ! I/O Parameters:
    ! ASP = asperity variables

    type(asperity_node),intent(inout) :: asp

    if (allocated(asp%data)) deallocate(asp%data)

  end subroutine destroy_asperity_node

  
  subroutine get_asperity_node(list,i_asp,asp,data,idata)
    ! GET_ASPERITY_NODE gets asperity data from node
    ! 
    ! Modified: 26 November 2006

    implicit none

    ! I/O Parameters:
    ! LIST = asperity list
    ! I_ASP = index of node from which data is extracted
    ! ASP = asperity variables
    ! DATA = perturbed fields amplitude data
    ! IDATA = index of data to be returned

    type(asperity_list),intent(in) :: list
    integer(pin),intent(in) :: i_asp
    type(asperity_type),intent(out) :: asp
    real(pr),dimension(:),intent(out),optional :: data    
    integer(pin),intent(in),optional :: idata

    ! Internal Parameters:
    ! CURRENT = pointer to current node in list
    ! I = index of node

    type(asperity_node),pointer :: current=>null()
    integer(pin) :: i

    ! start at beginning of list, move down it until i=i_asp

    current => list%head
    i = 1
    do while(i/=i_asp)
       i = i+1
       current => current%next
    end do

    ! move data from this node to output variables

    asp = current%asp

    if (allocated(current%data)) then
       if (present(idata)) then
          data = current%data(idata)
       else
          data = current%data
       end if
    end if

  end subroutine get_asperity_node


  function present_asperity_list(list)
    ! PRESENT_ASPERITY_LIST checks if list exists
    ! 
    ! Modified: 20 February 2005

    implicit none

    ! I/O Parameters:
    ! LIST = asperity list
    ! PRESENT_ASPERITY_LIST = whether or not list exists
    !      T = list exists
    !      F = list does not exist

    type(asperity_list),intent(in) :: list
    logical :: present_asperity_list

    present_asperity_list = associated(list%head)

  end function present_asperity_list


  function get_num_asperity_nodes(list) result(n)
    ! GET_NUM_ASPERITY_NODES returns number of nodes in list
    ! 
    ! Modified: 26 November 2006

    implicit none

    ! I/O Parameters:
    ! LIST = asperity list
    ! N = number of nodes in list

    type(asperity_list),intent(in) :: list
    integer(pin) :: n

    ! Internal Parameters:
    ! CURRENT = pointer to current node in list

    type(asperity_node),pointer :: current=>null()

    ! start at beginning of list, point current to head
    current => list%head

    ! move down list one node at a time, counting number of nodes
    if (associated(current)) then
       n = 1
       do while(associated(current%next))
          current => current%next
          n = n+1
       end do
    else
       n = 0
    end if

  end function get_num_asperity_nodes


  subroutine destroy_asperity_list(list)
    !  DESTROY_ASPERITY_LIST destroys all nodes in a list
    ! 
    ! Modified: 26 November 2006

    implicit none

    ! I/O Parameters:
    ! LIST = asperity list

    type(asperity_list),intent(inout) :: list

    ! Internal Parameters:
    ! CURRENT = pointer to current node in list
    ! PREVIOUS = pointer to previous node in list

    type(asperity_node),pointer :: current=>null(),previous=>null()

    ! check if list exists, return if it doesn't

    if (.not.present_asperity_list(list)) return

    ! start at beginning, point previous to first node in list, point 
    ! current to next one down

    previous => list%head
    current => list%head%next

    ! move down the list, removing second node and replacing it with 
    ! the third node

    do while(associated(previous%next))
       ! point previous%next to third node
       previous%next => current%next
       ! nullify pointer from second to third node
       nullify(current%next)
       ! destroy objects in second node
       call destroy_asperity_node(current)
       ! deallocate second node
       deallocate(current)
       ! point current to what was the third
       ! and is now the second node
       current => previous%next
    end do

    ! now remove the head node

    call destroy_asperity_node(list%head)

    deallocate(list%head)

  end subroutine destroy_asperity_list


  function ampl(x,y,dx,dy,asp,dim)
    ! AMPL returns the amplitude of the perturbations as a normalized 
    ! value between 0 and 1, given asperity variables and a position
    ! 
    ! Modified: 16 October 2008

    use constants, only : zero,half,one,two,twopi
    use io, only : error
    implicit none

    ! I/O Parameters:
    ! X = position in x direction
    ! Y = position in y direction
    ! DX = grid spacing in x direction
    ! DY = grid spacing in y direction
    ! ASP = asperity variables
    ! AMPL = amplitude of perturbation
    ! DIM = dimension
    
    real(pr),intent(in) :: x,y,dx,dy
    type(asperity_type),intent(in) :: asp
    integer(pin),intent(in) :: dim
    real(pr) :: ampl
    
    ! Internal Parameters:
    ! RAD = radius, used for several asperity types
    ! X1 = left side of grid in x direction
    ! X2 = right side of grid in x direction
    ! Y1 = bottom side of grid in y direction
    ! Y2 = top side of grid in y direction
    ! WEIGHT_X = fraction of grid cell covered by asperity in x direction
    ! WEIGHT_Y = fraction of grid cell covered by asperity in y direction

    real(pr) :: rad,x1,x2,y1,y2,weight_x,weight_y

    ! select amongst types of asperities

    select case(asp%type)

    case default ! wrong asperity type

       call error('Invalid asperity type','ampl')

    case('none') ! no asperity

       ampl = zero
       
       ! asperities with smooth edges

    case('uniform') ! uniform perturbation

       ampl = one

    case('fourier') ! 2D Fourier mode

       ampl = cos(twopi*(x-asp%x0)/asp%length_x)
       if (dim==3) ampl = ampl*cos(twopi*(y-asp%y0)/asp%length_y)

    case('gaussian') ! gaussian

       if (dim==3) then
          rad = (x-asp%x0)**2/(two*asp%length_x**2)+(y-asp%y0)**2/(two*asp%length_y**2)
       else
          rad = (x-asp%x0)**2/(two*asp%length_x**2)
       end if
       ampl = exp(-rad)
       
    ! asperities with discontinuous edges

    case('ellipse') ! ellipse

       if (dim==3) then
          rad = ((x-asp%x0)/asp%length_x)**2+((y-asp%y0)/asp%length_y)**2
       else
          rad = ((x-asp%x0)/asp%length_x)**2
       end if
       if (rad<one) then
          ampl = one
       elseif (rad==one) then
          ampl = half
       else
          ampl = zero
       end if

    case('box') ! box

       ! define grid boundaries relative to asperity center
       x1 = x-half*dx-asp%x0
       x2 = x+half*dx-asp%x0
       y1 = y-half*dy-asp%y0
       y2 = y+half*dy-asp%y0

       ! calculate fraction of cell that is inside box in x direction
       if (x2<=-asp%length_x) then 
          ! outside box to left
          weight_x = zero
       elseif (x1>=asp%length_x) then 
          ! outside box to right
          weight_x = zero
       elseif (x1<-asp%length_x.and.x2>-asp%length_x) then 
          ! partially inside box to left
          weight_x = (x2-(-asp%length_x))/dx
       elseif(x1<asp%length_x.and.x2>asp%length_x) then 
          ! partially inside box to right
          weight_x = (asp%length_x-x1)/dx
       else 
          ! inside box
          weight_x = one
       end if

       if (dim==3) then
          ! calculate fraction of cell that is inside box in y direction
          if (y2<=-asp%length_y) then 
             ! outside box to left
             weight_y = zero
          elseif (y1>=asp%length_y) then 
             ! outside box to right
             weight_y = zero
          elseif (y1<-asp%length_y.and.y2>-asp%length_y) then 
             ! partially inside box to left
             weight_y = (y2-(-asp%length_y))/dy
          elseif(y1<asp%length_y.and.y2>asp%length_y) then 
             ! partially inside box to right
             weight_y = (asp%length_y-y1)/dy
          else 
             ! inside box
             weight_y = one
          end if
       else
          weight_y = one
       end if

       ! multiply fractions to get area of cell inside box

       ampl = weight_x*weight_y

    case('smooth') ! all derivative continuous (C_infinity)

       if (dim==3) then
          rad = ((x-asp%x0)/asp%length_x)**2+((y-asp%y0)/asp%length_y)**2
       else
          rad = ((x-asp%x0)/asp%length_x)**2
       end if

       if (rad>=one) then
          ampl = zero
       else
          ampl = exp(rad/(rad-one))
       end if

    end select
    
    ! check whether or not to apply the perturbation inside or outside the 
    ! asperity region - this can be used for making borders of a given shape

    if (.not.asp%inside) ampl = one-ampl

  end function ampl


  function read_asperity_list(unit,list_name,ndata) result(list)
    ! READ_ASPERITY_LIST reads asperity list data
    !
    ! Modified: 29 July 2005

    use io, only : error
    implicit none

    ! I/O Parameters:
    ! UNIT = file unit to be read
    ! LIST_NAME = name of list
    ! NDATA = number of data items
    ! LIST = asperity list

    integer(pin),intent(in) :: unit,ndata
    character(*),intent(in) :: list_name
    type(asperity_list) :: list

    ! Internal Parameters:
    ! STAT = I/O error flag
    ! IDATA = index for data
    ! ASP = asperity variables
    ! DATA = asperity data
    ! NAME = name of asperity
    ! TEMP = temporary string
    ! ERR = error message
    ! START = string denoting start of list
    ! FINISH = string denoting end of list

    integer(pin) :: stat,idata
    type(asperity_type) :: asp
    real(pr) :: data(ndata)
    character(64) :: name,temp,err,start,finish

    start = '!---begin:' // trim(adjustl(list_name)) // '_asperity_list---'
    finish = '!---end:' // trim(adjustl(list_name)) // '_asperity_list---'
    err = "Error reading asperity list '" // trim(adjustl(list_name)) // "' in .in file"

    ! move to start of asperity list

    rewind(unit)
    do
       read(unit,*,end=100,iostat=stat) temp
       if (stat/=0) goto 100
       if (temp==start) exit
    end do
    
    ! read asperities and add nodes to list
    
    do
       ! read item data
       read(unit,*,end=100,iostat=stat) name
       if (stat/=0) goto 100
       if (name==finish) goto 200
       read(unit,*,iostat=stat) asp%type,asp%inside
       if (stat/=0) goto 100
       read(unit,*,iostat=stat) asp%x0,asp%length_x
       if (stat/=0) goto 100
       read(unit,*,iostat=stat) asp%y0,asp%length_y
       read(unit,*,iostat=stat) (data(idata),idata=1,ndata)
       if (stat/=0) goto 100
       ! store as node in list
       call create_asperity_node(list,asp,data)
    end do
    
100 call error(err,'read_asperity_list')
    
    ! come here when end of file is reached
200 continue
    
  end function read_asperity_list


  subroutine assign_list_data(array,list,idata,x,y,dx,dy,dim)
    ! ASSIGN_LIST_DATA uses the list data to perturb the contents of array
    !
    ! Modified: 27 February 2006

    implicit none

    ! I/O Parameters:
    ! ARRAY = array to be perturbed
    ! LIST = asperity list
    ! IDATA = index of data value containing perturbation amplitude
    ! X = coordinate vector in x direction
    ! Y = coordinate vector in y direction
    ! DX = grid spacing in x direction
    ! DY = grid spacing in y direction
    ! DIM = dimension

    real(pr),dimension(:,:),intent(inout) :: array
    type(asperity_list),intent(in) :: list
    integer(pin),intent(in) :: idata,dim
    real(pr),dimension(:),intent(in) :: x,y
    real(pr),intent(in) :: dx,dy

    ! Internal Parameters:
    ! N_ASP = number of asperities
    ! I_ASP = index of asperity
    ! NX = number of grid points in x direction
    ! NY = number of grid points in y direction
    ! ASP = asperity variables
    ! DATA = asperity data
    ! I = index in x direction
    ! J = index in y direction
    ! AMP = amplitude of the perturbation

    integer(pin) :: i_asp,n_asp,nx,ny,i,j
    type(asperity_type) :: asp
    real(pr) :: data(1:1),amp

    if (.not.present_asperity_list(list)) return

    nx = size(array,1)
    ny = size(array,2)

    n_asp = get_num_asperity_nodes(list)

    ! loop over all asperities

    do i_asp = 1,n_asp

       ! get asp and data
       call get_asperity_node(list,i_asp,asp,data,idata)

       ! loop over all points
       do j = 1,ny
          do i = 1,nx
             ! get amplitude of perturbation
             amp = ampl(x(i),y(j),dx,dy,asp,dim)
             ! perturb array
             array(i,j) = array(i,j)+amp*data(1)
          end do
       end do

    end do

  end subroutine assign_list_data


end module asperity
