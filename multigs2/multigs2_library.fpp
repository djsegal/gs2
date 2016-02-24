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
    integer :: sub_comm,sub_comm_id !the sub communicator made by mpi split and a label to identify which processor is in which communicator
    integer :: no_electron_boxes,no_simulations ! the number of simulations will be no_electron_boxes +1 
    integer :: world_id,sub_id !processor id, within the world communicator and within in the sub communicator made by mpi_split
    character(2000) :: run_name
    type(gs2_program_state_type) :: gs2_state
  end type multigs2_library_type

  contains



! This subroutine initializes the multiscale gs2 
    subroutine initialize_multigs2(multigs2_obj)
! the next command is a preprocessor command which tells the compiler to only include these next lines is MPI is specified on build - taken care of by build script MRH
#ifdef MPI
! - makes no sense to have multigs2 without installing mpi, so make this unconditional?
      use mpi, only: mpi_comm_world
      use mpi, only: mpi_init
      use mpi, only: mpi_comm_size
      use mpi, only: mpi_comm_rank
      use mpi, only: mpi_comm_split
#endif
      use command_line, only: cl_getarg
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      integer :: length
      character(2000) :: temp
      integer ::  ierror
      integer :: ntot_proc,procperbox
      integer,allocatable :: sub_group_label(:) ! this will be a label for each processor - telling it which sub_communicator is belongs to
     ! integer :: world_id,sub_id !first is the world rank or id, second will be the rank or id of the processor within the sub_communicator
      integer :: i,j ! loop variable
     ! logical, parameter :: reorder=.true.
     ! integer, dimension(2) :: dims
     ! logical, dimension(2) :: period
      
#ifndef MPI

      write(*,*) "multigs2 requires MPI capability to run - exiting"

      stop

#endif      

      ! set the number of electron boxes - this needs to be input via an input file or cmd line later
  
      multigs2_obj%no_electron_boxes =3 
      multigs2_obj%no_simulations =  multigs2_obj%no_electron_boxes +1 

! this bit of code sets the mpi communicator to be the world communicator if multiscale gs2 is not being run inside something else already.
      if (.not. multigs2_obj%mp_comm_external) then 
#ifdef MPI 
        call mpi_init (ierror)  ! calls the mpi module for the first time MRH
        multigs2_obj%mp_comm = mpi_comm_world
#endif
      end if
#ifdef MPI
      call mpi_comm_size (multigs2_obj%mp_comm, ntot_proc, ierror) !given the communicator mp comm, this returns ntot_proc, the number of processors running
#endif    
     ! ncolumns =  multigs2_obj%no_electron_boxes
      procperbox = ntot_proc/multigs2_obj%no_simulations ! this is the number of processess per electron box

     

! check here that we are not wasting processors 
       
      if(ntot_proc /= multigs2_obj%no_simulations*procperbox) then
        ierror = 1
        write(*,*) 'Number of processes must be divisible by number of ionscale (1)  plus electronscale simulations'
        return
      endif

! end check 

#ifdef MPI         
      call mpi_comm_rank(multigs2_obj%mp_comm,multigs2_obj%world_id,ierror) !this tells us the world_id of this processor MRH
#endif

! now arrange for the processors to be grouped together MRH
      allocate(sub_group_label(ntot_proc))

      do i=0,multigs2_obj%no_simulations-1,1
        do j=1,procperbox,1
           sub_group_label(i*procperbox + j) = i 
        end do
      end do
      write(*,*) sub_group_label      
      
!     do i=0,ntot_proc-1,1 ! sub_id (ranks) range from 0 to N-1 where N is ntot_proc - this is a feature of MPI MRH
!        sub_group_label(i) = mod(i,multigs2_obj%no_simulations) ! this allocates an equal number of processors per box MRH
!      end do

      multigs2_obj%sub_id = mod(multigs2_obj%world_id,procperbox)
      multigs2_obj%sub_comm_id = sub_group_label(multigs2_obj%world_id+1) ! plus 1 here because world_id goes from 0 to 1, fortran indices from 1 to ntot
#ifdef MPI
      call mpi_comm_split(multigs2_obj%mp_comm,multigs2_obj%sub_comm_id,&
           multigs2_obj%sub_id, multigs2_obj%sub_comm, ierror)
#endif
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


 
! These subroutines amount to almost identical copy of gs2.f90 here
! there are different routines to intialise the ion scale and electron scale versions of gs2
! however, once the versions are initialised, the code to run and finalize each gs2 version are identical
    subroutine initialize_ion_box(multigs2_obj)
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
   
    end subroutine initialize_ion_box
    
    subroutine initialize_electron_box(multigs2_obj)
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
    
      character(2000)::char_label
! convert label into a string
! here sub_comm_id will go from 1 to N where N is the number of electronscale flux tubes
! this is an example of fortran's internal write statement 
      write(char_label,'(I0)') multigs2_obj%sub_comm_id 
      char_label = trim(adjustl(char_label)) 
 
! now use this to select the correct input file
      multigs2_obj%gs2_state%run_name_external = .true.
      multigs2_obj%gs2_state%run_name = trim(multigs2_obj%run_name)//'_electron'//char_label
      multigs2_obj%gs2_state%mp_comm_external = .true.
      multigs2_obj%gs2_state%mp_comm = multigs2_obj%sub_comm
      call initialize_gs2(multigs2_obj%gs2_state)
      call initialize_equations(multigs2_obj%gs2_state)
      call initialize_diagnostics(multigs2_obj%gs2_state)
      multigs2_obj%gs2_state%print_times = .false.
   
    end subroutine initialize_electron_box
! this subroutine runs an already initialized ocurrance of gs2

    subroutine run_box(multigs2_obj)
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
     

      if (multigs2_obj%gs2_state%do_eigsolve) then 
        call run_eigensolver(multigs2_obj%gs2_state)
      else
        call evolve_equations(multigs2_obj%gs2_state, multigs2_obj%gs2_state%nstep)
      end if
    end subroutine run_box
!this finalizes an occurance of gs2
    subroutine finalize_box(multigs2_obj)
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
      call finalize_diagnostics(multigs2_obj%gs2_state)
      call finalize_equations(multigs2_obj%gs2_state)
      multigs2_obj%gs2_state%print_times = .true.
      multigs2_obj%gs2_state%print_full_timers = .true.
      call finalize_gs2(multigs2_obj%gs2_state)
    end subroutine finalize_box



!this subroutine closes the mpi routine, and will eventually contain the finalize multigs2 routines
    subroutine finalize_multigs2(multigs2_obj)
#ifdef MPI       
      use mpi, only: mpi_finalize
#endif
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      integer :: ierror
#ifdef MPI
      call mpi_finalize(ierror)
#endif
    end subroutine finalize_multigs2

    subroutine communicate_omega(multigs2_obj)

      use diagnostics_omega, only: omega_average    
#ifdef MPI       
       use mpi, only: mpi_complex
       use mpi, only: mpi_status_size 
!      use mpi, only: mpi_send
!      use mpi, only: mpi_recv
! for some reason this statement causes a compile error, although one can use the subroutines without these statements
#endif
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      integer :: ierror
      integer :: mpi_status(mpi_status_size)
      integer :: i ! loop variable
      integer :: tag_omega =1
      complex, dimension(:,:,:), allocatable :: omega_from_all_gs2
      
      allocate(omega_from_all_gs2(size(omega_average,1),size(omega_average,2),&
          multigs2_obj%no_simulations))      
      write(*,*) omega_average

#ifdef MPI
 
      if(sub_comm_id == 0 .and. sub_id == 0) then
        omega_from_all_gs2(:,:,1) = omega_average
        do i = 1,multigs2_obj%no_simulations,1
          call mpi_recv(omega_from_all_gs2(:,:,i), size(omega_from_all_gs2(:,:,i)),&
              mpi_complex,multigs2_obj%world_id,tag_omega,multigs2_obj%mp_comm,mpi_status,ierror)! do a thing
      end if 

!      if(sub_comm_id > 0 .and. sub_id == 0) then 
!        call mpi_send ! do a thing
!      end if 
     
 ! call MPI_RECV(slave_obj%matrix,&
 !          size(slave_obj%matrix),MPI_DOUBLE_PRECISION,0,tag_newtask,sub_COMM,mpi_status,ierr)
!call MPI_SEND(slave_obj%matrix_inv,&
 !              size(slave_obj%matrix_inv),MPI_DOUBLE_PRECISION,0,tag_results,sub_COMM,ierr)
      
#endif

      deallocate(omega_from_all_gs2)
    end subroutine communicate_omega


end module multigs2_library
