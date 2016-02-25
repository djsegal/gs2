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
  use mp, only: comm_type
  implicit none

  type multigs2_library_type
    logical :: multigs2_comm_external = .false. 
!    integer :: multigs2_comm !this is multigs2's communicator
    logical :: run_name_external = .false.
!    integer :: world_id !processor id, within the multigs2 communicator (i.e. multigs2_comm)
    type(comm_type) :: gs2_comm
    type(comm_type) :: multigs2_comm
  ! integer :: sub_comm !the sub communicator made by mpi split and a label to identify which processor is in which communicator
    integer :: gs2_comm_id ! integer 0, 1 - N to identify what sort of simulation is being done by gs2, 0 ionscale 1-N electronscale
    integer :: no_electron_boxes,no_simulations ! the number of simulations will be no_electron_boxes +1 
 !   integer :: world_id,sub_id !processor id, within the world communicator and within in the sub communicator made by mpi_split
    character(2000) :: run_name
    type(gs2_program_state_type) :: gs2_state
  end type multigs2_library_type
! MRH From mp.fpp for convenience (default values shown)
! MRH  type comm_type
! MRH    integer :: id=-1 !The communicator id, used in calls to MPI routines. Contrary to the name, this does not seem to be a unique identifier, rather it is a description of the type of communicator forthe mpi
! MRH    integer :: iproc=-1 !The procs local rank
! MRH    integer :: nproc=-1 !The total number of processors in the communicator
! MRH    logical :: proc0=.false. !Is iproc equal to 0?
! MRH  end type comm_type
  contains



! This subroutine initializes the multiscale gs2 
    subroutine initialize_multigs2(multigs2_obj)
! the next command is a preprocessor command which tells the compiler to only include these next lines is MPI is specified on build - taken care of by build script MRH
!#ifdef MPI
! - makes no sense to have multigs2 without installing mpi, so make this unconditional?
!      use mpi, only: mpi_comm_world
!      use mpi, only: mpi_init
!      use mpi, only: mpi_comm_size
!      use mpi, only: mpi_comm_rank
!      use mpi, only: mpi_comm_split
!#endif One doesn't need preprocessor statements around mp subroutines
      use mp, only : init_mp ! nearly equiv to mpi_initialize
      use mp, only : mp_comm ! holds mpi_comm_world (it seems)
      use mp, only : nproc! give the number of processors
      use mp, only : iproc ! gives the rank
      use mp, only : proc0
  !     use mp, only : split_key_to_commtype ! makes new sub communicators using mpi_split, using full colour and key functionality
      use mp, only : split
      use command_line, only: cl_getarg
      implicit none
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      integer :: length
      character(2000) :: temp
      integer ::  ierror,col,key
     ! integer :: ntot_proc
      integer :: procperbox
      integer,allocatable :: sub_group_label(:) ! this will be a label for each processor - telling it which sub_communicator is belongs to
     ! integer :: world_id,sub_id !first is the world rank or id, second will be the rank or id of the processor within the sub_communicator
      integer :: i,j ! loop variable
     ! logical, parameter :: reorder=.true.
     ! integer, dimension(2) :: dims
     ! logical, dimension(2) :: period
      logical :: multigs2flag = .true. ! this flag merely has to be present - the actual value is not important     
#ifndef MPI

      write(*,*) "multigs2 requires MPI capability to run - exiting"

      stop

#endif      

      ! set the number of electron boxes - this needs to be input via an input file or cmd line later
  
      multigs2_obj%no_electron_boxes =3 
      multigs2_obj%no_simulations =  multigs2_obj%no_electron_boxes +1 

! this bit of code sets the mpi communicator to be the world communicator if multiscale gs2 is not being run inside something else already.
      if (.not. multigs2_obj%multigs2_comm_external) then 
!#ifdef MPI 
        call init_mp(multigs2=multigs2flag)  ! calls the mp module for the first time, with the flag making it with pointers appropriate to multigs2 MRH
        multigs2_obj%multigs2_comm%id = mp_comm ! and sets multigs2_comm to be mpi_comm_world
      else  ! this option is for if multigs2 is being run inside something else and multigs2_comm is already set
        call init_mp(multigs2_obj%multigs2_comm%id,multigs2=multigs2flag)
      end if
  
      multigs2_obj%multigs2_comm%iproc = iproc ! the rank of this proc within the multigs2 communicator
      multigs2_obj%multigs2_comm%nproc = nproc ! the total number of procs running in multigs2
      multigs2_obj%multigs2_comm%proc0 = proc0 ! true if iproc =0

!#ifdef MPI
   !   call nproc_comm(multigs2_obj%multigs2_comm,ntot_proc) !given the communicator mp comm, this returns ntot_proc, the number of processors running
!#endif    
     ! ncolumns =  multigs2_obj%no_electron_boxes
      procperbox = multigs2_obj%multigs2_comm%nproc/multigs2_obj%no_simulations ! this is the number of processess per electron box

     

! check here that we are not wasting processors 
       
      if(multigs2_obj%multigs2_comm%nproc /= multigs2_obj%no_simulations*procperbox) then
        ierror = 1
        write(*,*) 'Number of processes must be divisible by number of ionscale (1)  plus electronscale simulations'
        return
      endif

! end check 

!#ifdef MPI         
    !  call rank_comm(multigs2_obj%multigs2_comm,multigs2_obj%world_id) !this tells us the world_id of this processor MRH
!#endif

! now arrange for the processors to be grouped together MRH
      allocate(sub_group_label(multigs2_obj%multigs2_comm%nproc))

      do i=0,multigs2_obj%no_simulations-1,1
        do j=1,procperbox,1
           sub_group_label(i*procperbox + j) = i 
        end do
      end do
      write(*,*) sub_group_label      
      
!     do i=0,ntot_proc-1,1 ! sub_id (ranks) range from 0 to N-1 where N is ntot_proc - this is a feature of MPI MRH
!        sub_group_label(i) = mod(i,multigs2_obj%no_simulations) ! this allocates an equal number of processors per box MRH
!      end do
!      key = multigs2_obj%world_id
      key = mod(multigs2_obj%multigs2_comm%iproc,procperbox)
      col = sub_group_label(multigs2_obj%multigs2_comm%iproc+1) ! plus 1 here because world_id goes from 0 to 1, fortran indices from 1 to ntot
      multigs2_obj%gs2_comm_id = col
!#ifdef MPI
 !     call mpi_comm_split(multigs2_obj%multigs2_comm,multigs2_obj%sub_comm_id,&
  !         multigs2_obj%sub_id, multigs2_obj%sub_comm, ierror)
!#endif
      call split(col,key,multigs2_obj%gs2_comm) 
     ! call split_key_to_commtype (col,key,multigs2_obj%gs2_comm)
 ! MRH here those processors with the same col are in the same communicator, 
 ! MRH key assigns a rank within the new communicator
 ! MRH gs2_comm keeps track of this for us? so we do not have to store col and key seperately in multigs2_obj?
 ! MRH gs2_comm also handily contains logical which tells us if it is the zero processor. 
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
      implicit none
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      multigs2_obj%gs2_state%run_name_external = .true.
      multigs2_obj%gs2_state%run_name = trim(multigs2_obj%run_name)//'_ion'
      multigs2_obj%gs2_state%mp_comm_external = .true.
      multigs2_obj%gs2_state%mp_comm = multigs2_obj%gs2_comm%id
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
      write(char_label,'(I0)') multigs2_obj%gs2_comm_id 
      char_label = trim(adjustl(char_label)) 
 
! now use this to select the correct input file
      multigs2_obj%gs2_state%run_name_external = .true.
      multigs2_obj%gs2_state%run_name = trim(multigs2_obj%run_name)//'_electron'//char_label
      multigs2_obj%gs2_state%mp_comm_external = .true.
      multigs2_obj%gs2_state%mp_comm = multigs2_obj%gs2_comm%id
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
      implicit none
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      call finalize_diagnostics(multigs2_obj%gs2_state)
      call finalize_equations(multigs2_obj%gs2_state)
      multigs2_obj%gs2_state%print_times = .true.
      multigs2_obj%gs2_state%print_full_timers = .true.
      call finalize_gs2(multigs2_obj%gs2_state)
    end subroutine finalize_box



!this subroutine closes the mpi routine, and will eventually contain the finalize multigs2 routines
    subroutine finalize_multigs2(multigs2_obj)
       
      use mp, only: finish_mp   

      implicit none
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      if (.not. multigs2_obj%multigs2_comm_external) call finish_mp

    end subroutine finalize_multigs2

    subroutine communicate_omega(multigs2_obj,test1,test2)

      use diagnostics_omega, only: omega_average    
      use mp, only: sum_reduce
      use mp, only: barrier
     ! use mp, only:mp_comm,iproc
      use mp, only: scope !confusingly, allprocs is the gs2 proc group label for scope
      use mp, only: multigs2procs
      use mp, only: allprocs
      use unit_tests, only : agrees_with 
      implicit none
      type(multigs2_library_type), intent(inout) :: multigs2_obj
      
      integer :: i ! loop variable
      real :: eps ! the tolerance for comparing numbers
      complex, dimension(:,:,:), allocatable :: omega_from_all_gs2
      logical, intent(inout) :: test1,test2
      !logical :: testr1=.false.,testr2=.true.,testi1=.false.,testi2=.true.,test_all 
      
      allocate(omega_from_all_gs2(size(omega_average,1),size(omega_average,2),&
          multigs2_obj%no_simulations))      
!      write(*,*) omega_average
      omega_from_all_gs2 = 0.0 ! assign array to be all zeros initially
      ! test that assignment     
!      if((multigs2_obj%gs2_comm_id == 0) .and. &
!         ( multigs2_obj%gs2_comm%proc0 .eqv. .true.)) then
!        write(*,*) omega_from_all_gs2
!      end if 
 
      do i=0,multigs2_obj%no_electron_boxes,1
       if( (multigs2_obj%gs2_comm%proc0 .eqv. .true.) .and. &
         (multigs2_obj%gs2_comm_id == i) )  then

            omega_from_all_gs2(:,:,i+1) = omega_average
            write(*,*) omega_from_all_gs2
        end if 
      end do
      call barrier(multigs2_obj%multigs2_comm%id)! stops all processes that reach this point from proceeding until all have reached this point
      call scope(multigs2procs) 
    !  mp_comm = multigs2_obj%multigs2_comm ! assigns the internal mp_comm variable from mp to be the right communicator (here for multigs2 so we can communicate between gs2s)
    !  iproc = multigs2_obj%world_id  ! assigns the the interal iproc variable from mp to be the world_id of that processor (here for multigs2 so we can communicate between gs2s)

      call sum_reduce(omega_from_all_gs2,0) ! the zero here is the world_id of the destination of the result
      call scope(allprocs)

    !  mp_comm = multigs2_obj%gs2_comm%id     ! set mp_comm back now to be the local gs2 communicator for safety
    !  iproc = multigs2_obj%gs2_comm%iproc  ! set iproc back now to be the local gs2 processor rank

      ! test that assignment     
      if(multigs2_obj%multigs2_comm%iproc == 0) then
        write(*,*) omega_from_all_gs2
      end if 
       
      eps = 1.0e-7
      if (precision(eps).lt. 11) eps = eps * 1000.0

      if(multigs2_obj%multigs2_comm%iproc ==0) then
      ! assign the logicals to be sure of their starting value
      ! in here to avoid getting an answer that isn't due to the below test - so only defined on world_id=0
      test1 =.false.
      test2 =.true. 
       do i =2,multigs2_obj%no_simulations,1
         test1 = agrees_with(omega_from_all_gs2(:,:,1),omega_from_all_gs2(:,:,i),eps) .and.  test1 
         test2 = agrees_with(omega_from_all_gs2(:,:,2),omega_from_all_gs2(:,:,i),eps) .and. test2
       end do
       write(*,*) "test1= ",test1,"test2= ",test2       
      end if 
!      if(multigs2_obj%world_id ==0) then 
!        do i=1,multigs2_obj%no_simulations,1
! 
!          omega_from_all_gs2(:,:,i) = omega_from_all_gs2(:,:,i)-omega_average ! fortran supports array operations
!        end do
!        write(*,*) omega_from_all_gs2
!        
!        do i=1,size(omega_average,1),1
!          do j=1, size(omega_average,2),1
!            do k=2,multigs2_obj%no_simulations
!
!              testr1 =( real(omega_from_all_gs2(i,j,k)) < 0.001).and. testr1    
!
!              testr2 =( real(omega_from_all_gs2(i,j,k) - omega_from_all_gs2(i,j,2)) < 0.001) .and. testr2
!           
!              testi1 =( aimag(omega_from_all_gs2(i,j,k)) < 0.001).and. testi1    
!
!              testi2 =( aimag(omega_from_all_gs2(i,j,k) - omega_from_all_gs2(i,j,2)) < 0.001) .and. testr2
!           end do
!          end do 
!        end do
!        test_all = ((.not. testr1) .and. testr2 ).and.  ((.not. testi1) .and. testi2 )
!        write(*,*) "testr1 ",testr1, "testi1= ",testi1, "testr2= ",testr2, "testi2= ",testi2, "test_all= ",test_all
!
!       end if  
 !     if(multigs2_obj%world_id == 1) then
 !       write(*,*) omega_from_all_gs2
 !     end if 
        
      


      deallocate(omega_from_all_gs2)
    end subroutine communicate_omega


end module multigs2_library
