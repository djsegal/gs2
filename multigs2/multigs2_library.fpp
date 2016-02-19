!>  A library for running multiscale gs2, which
!! runs separate ion and electron scale simulations 
!! and communicates information between them.
!!
!! This is free software released under the MIT license.
!! Written by:
!!             Michael Hardman (mhardman1@users.sourceforge.net)
!!             Michael Barnes ()
!!             Edmund Highcock (edmundhighcock@users.sourceforge.net)
module multigs2_library
  use gs2_main, only: gs2_program_state_type
  implicit none

  type multigs2_library_type
    logical :: mp_comm_external = .false. 
    integer :: mp_comm
    logical :: run_name_external = .false.
    integer :: sub_comm
    integer :: no_electron_boxes,no_simulations ! the number of simulations will be no_electron_boxes +1 
    integer :: world_id,sub_id !processor id
    character(2000) :: run_name
    type(gs2_program_state_type) :: gs2_state
  end type multigs2_library_type

  contains



! This subroutine initializes the multiscale gs2 
    subroutine initialize_multigs2(multigs2_obj)
! the next command is a preprocessor command which tells the compiler to only include these next lines is MPI is specified on build - taken care of by build script MRH
#ifdef MPI 
      use mpi, only: mpi_comm_world
      use mpi, only: mpi_init
#endif
      use command_line, only: cl_getarg
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      integer :: length
      character(2000) :: temp
      integer ::  ierror
      integer :: ntot_proc,procperbox
      integer,allocatable :: sub_group_label(:) ! this will be a label for each processor - telling it which sub_communicator is belongs to
     ! integer :: world_id,sub_id !first is the world rank or id, second will be the rank or id of the processor within the sub_communicator
      integer :: i ! loop variable
     ! logical, parameter :: reorder=.true.
     ! integer, dimension(2) :: dims
     ! logical, dimension(2) :: period
      
     

      ! set the number of electron boxes - this needs to be input via an input file or cmd line later
  
      multigs2_obj%no_electron_boxes = 4 
      multigs2_obj%no_simulations =  multigs2_obj%no_electron_boxes +1 

! this bit of code sets the mpi communicator to be the world communicator if multiscale gs2 is not being run inside something else already.
      if (.not. multigs2_obj%mp_comm_external) then 
#ifdef MPI
        call mpi_init (ierror)  ! calls the mpi module for the first time MRH
        multigs2_obj%mp_comm = mpi_comm_world
#endif
      end if
      call mpi_comm_size (multigs2_obj%mp_comm, ntot_proc, ierror) !given the communicator mp comm, this returns ntot_proc, the number of processors running
    
     ! ncolumns =  multigs2_obj%no_electron_boxes
      procperbox = ntot_proc/multigs2_obj%no_simulations ! this is the number of processess per electron box

     

! check here that we are not wasting processors 
       
      if(ntot_proc /= multigs2_obj%no_simulations*procperbox) then
        ierror = 1
        write(*,*) 'Number of processes must be divisible by number of ionscale (1)  plus electronscale simulations'
        return
      endif

! end check 

         
      call mpi_comm_rank(multigs2_obj%mp_comm,multigs2_obj%world_id,ierror) !this tells us the world_id of this processor MRH
! now arrange for the processors to be grouped together MRH
      allocate(sub_group_label(ntot_proc))

      do i=0,ntot_proc-1,1 ! sub_id (ranks) range from 0 to N-1 where N is ntot_proc - this is a feature of MPI MRH
        sub_group_label(i) = mod(i,multigs2_obj%no_simulations) ! this allocates an equal number of processors per box MRH
      end do

      multigs2_obj%sub_id = mod(multigs2_obj%world_id,multigs2_obj%no_simulations)
      call mpi_comm_split(multigs2_obj%mp_comm,sub_group_label(multigs2_obj%world_id),&
           multigs2_obj%sub_id, multigs2_obj%sub_comm, ierror)

      deallocate(sub_group_label)
 ! call MPI_COMM_SPLIT(MPI_COMM_WORLD,integer :: second_arg, integer :: third_arg, SUB_COMMUNICATOR, ierr) 
 !Here SUB_COMMUNICATOR is the same for all processes with the same second_arg,
!  and the new my_id variable within this new communicator is given by the third_arg MRH 

! now this processor has been assigned a communicator and sub_id within this communicator MRH

       
! ncolumns is the number of electron simulations per ion simulation - presumably this should be linked to the number of points in the ion flux tube at some point MRH.
  !    nrows = ntot_proc/ncolumns
  !    dims=(/ ncolumns, nrows /)     
  !    if(ntot_proc /= ncolumns*nrows) then
  !       ierror = 1
  !       write(*,*) 'Number of processes must be divisible by number of electron simulations'
  !       return
  !    endif
      ! create 2d cartesian topology for processes
  !    period=(/ .false., .false. /)  !! no circular shift
  !    ndim = 2
  !    call mpi_cart_create(multigs2_obj%mp_comm, ndim, dims, period, reorder, multigs2_obj%sub_comm, ierror)

      if (.not. multigs2_obj%run_name_external) then 
        call cl_getarg(1, temp, length, ierror)
        multigs2_obj%run_name = temp(1:length-3) ! this takes characters 1 to the 4th till last as the run name somehow MRH
      end if 
  
     
    ! <doc> cl_getarg (k,arg,length,ierror)
    !  gets k-th argument string and its length
    !  using intrinsic getarg or POSIX pxfgetarg
    ! </doc>

 

   end subroutine initialize_multigs2



! this subroutine runs a global ion flux tube using gs2 MRH 
! almost identical copy of gs2.f90 here
    subroutine run_ion_box(multigs2_obj)
      use gs2_main, only: gs2_program_state_type
      use gs2_main, only: initialize_wall_clock_timer
      use gs2_main, only: initialize_gs2
      use gs2_main, only: initialize_equations
      use gs2_main, only: initialize_diagnostics
      use gs2_main, only: evolve_equations
      use gs2_main, only: run_eigensolver
      use gs2_main, only: finalize_diagnostics
      use gs2_main, only: finalize_equations
      use gs2_main, only: finalize_gs2
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      multigs2_obj%gs2_state%run_name_external = .true.
      multigs2_obj%gs2_state%run_name = trim(multigs2_obj%run_name)//'_ion'
      multigs2_obj%gs2_state%mp_comm_external = .true.
      multigs2_obj%gs2_state%mp_comm = multigs2_obj%sub_comm
      call initialize_gs2(multigs2_obj%gs2_state)
      call initialize_equations(multigs2_obj%gs2_state)
      call initialize_diagnostics(multigs2_obj%gs2_state)
      multigs2_obj%gs2_state%print_times = .false.
      if (multigs2_obj%gs2_state%do_eigsolve) then 
        call run_eigensolver(multigs2_obj%gs2_state)
      else
        call evolve_equations(multigs2_obj%gs2_state, multigs2_obj%gs2_state%nstep)
      end if
      call finalize_diagnostics(multigs2_obj%gs2_state)
      call finalize_equations(multigs2_obj%gs2_state)
      multigs2_obj%gs2_state%print_times = .true.
      multigs2_obj%gs2_state%print_full_timers = .true.
      call finalize_gs2(multigs2_obj%gs2_state)
    end subroutine run_ion_box



!this subroutine runs an individual electron flux tube within the ion flux tube using gs2 MRH
    subroutine run_electron_box(multigs2_obj)
      ! need to make almost a copy of gs2.f90 here MRH. 
      use gs2_main, only: gs2_program_state_type
      use gs2_main, only: initialize_wall_clock_timer
      use gs2_main, only: initialize_gs2
      use gs2_main, only: initialize_equations
      use gs2_main, only: initialize_diagnostics
      use gs2_main, only: evolve_equations
      use gs2_main, only: run_eigensolver
      use gs2_main, only: finalize_diagnostics
      use gs2_main, only: finalize_equations
      use gs2_main, only: finalize_gs2
      
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      
      
      multigs2_obj%gs2_state%run_name_external = .true.
      multigs2_obj%gs2_state%run_name = trim(multigs2_obj%run_name)//'_electron'
      multigs2_obj%gs2_state%mp_comm_external = .true.
      multigs2_obj%gs2_state%mp_comm = multigs2_obj%sub_comm
      call initialize_gs2(multigs2_obj%gs2_state)
      call initialize_equations(multigs2_obj%gs2_state)
      call initialize_diagnostics(multigs2_obj%gs2_state)
      multigs2_obj%gs2_state%print_times = .false.
      if (multigs2_obj%gs2_state%do_eigsolve) then 
        call run_eigensolver(multigs2_obj%gs2_state)
      else
        call evolve_equations(multigs2_obj%gs2_state, multigs2_obj%gs2_state%nstep)
      end if
      call finalize_diagnostics(multigs2_obj%gs2_state)
      call finalize_equations(multigs2_obj%gs2_state)
      multigs2_obj%gs2_state%print_times = .true.
      multigs2_obj%gs2_state%print_full_timers = .true.
      call finalize_gs2(multigs2_obj%gs2_state)
    end subroutine run_electron_box

    subroutine finalize_multigs2(multigs2_obj)


      type(multigs2_library_type), intent(inout) :: multigs2_obj
      integer :: ierror

      call mpi_finalize(ierror)

    end subroutine finalize_multigs2

end module multigs2_library
