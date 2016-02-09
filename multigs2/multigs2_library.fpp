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
    integer :: electron_comm
    character(2000) :: run_name
    type(gs2_program_state_type) :: gs2_state
  end type multigs2_library_type

  contains
    subroutine initialize_multigs2(multigs2_obj)
#ifdef MPI
      use mpi, only: mpi_comm_world
      use mpi, only: mpi_init
#endif
      use command_line, only: cl_getarg
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      integer :: length, ierror, ntot_proc, nrows, ncolumns, ndim
      logical, parameter :: reorder=.true.
      integer, dimension(2) :: dims
      logical, dimension(2) :: period
      
      character(2000) :: temp
      if (.not. multigs2_obj%mp_comm_external) then 
#ifdef MPI
        call mpi_init (ierror)  ! MAB
        multigs2_obj%mp_comm = mpi_comm_world
#endif
      end if
      call mpi_comm_size (multigs2_obj%mp_comm, ntot_proc, ierror)
      ncolumns = 4
      nrows = ntot_proc/ncolumns
      dims=(/ ncolumns, nrows /)     
      if(ntot_proc /= ncolumns*nrows) then
         ierror = 1
         write(*,*) 'Number of processes must be divisible by number of groups'
         return
      endif
      ! create 2d cartesian topology for processes
      period=(/ .false., .false. /)  !! no circular shift
      ndim = 2
      call mpi_cart_create(multigs2_obj%mp_comm, ndim, dims, period, reorder, multigs2_obj%electron_comm, ierror)
      if (.not. multigs2_obj%run_name_external) then 
        call cl_getarg(1, temp, length, ierror)
        multigs2_obj%run_name = temp(1:length-3)
      end if 
    end subroutine initialize_multigs2

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
      multigs2_obj%gs2_state%mp_comm = multigs2_obj%mp_comm
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

    subroutine run_electron_box(multigs2_obj)
      type(multigs2_library_type), intent(inout) :: multigs2_obj
    end subroutine run_electron_box

end module multigs2_library
