module gs2_gryfx_zonal

  use gs2_main, only: gs2_program_state_type
  use iso_c_binding

  implicit none
  private

  type(gs2_program_state_type) :: state
  public :: state
  public :: init_gs2_gryfx, advance_gs2_gryfx, finish_gs2_gryfx
  public :: allocate_gryfx_zonal_arrays
  public :: deallocate_gryfx_zonal_arrays
  public :: gryfx_parameters_type
 
  public :: test_flag

  logical :: test_flag = .false.

  type, bind(c) :: gryfx_parameters_type
      integer(c_int) :: mpirank

       ! Name of gryfx/gryfx input file
       !character(len=1000) :: input_file
      !Base geometry parameters - not currently set by trinity 
      !See geometry.f90
       integer(c_int) :: equilibrium_type
       !character(len=800) :: eqfile
       integer(c_int) :: irho
       real(c_double) :: rhoc
       integer(c_int) :: bishop
       integer(c_int) :: nperiod
       integer(c_int) :: ntheta

      ! Miller parameters
       real(c_double) :: rgeo_lcfs
       real(c_double) :: rgeo_local
       real(c_double) :: akappa
       real(c_double) :: akappri
       real(c_double) :: tri
       real(c_double) :: tripri
       real(c_double) :: shift
       real(c_double) :: qinp
       real(c_double) :: shat
       real(c_double) :: asym
       real(c_double) :: asympri

       ! Other geometry parameters - Bishop/Greene & Chance
       real(c_double) :: beta_prime_input
       real(c_double) :: s_hat_input

       ! Flow shear
       real(c_double) :: g_exb
       ! Species parameters... I think allowing 20 species should be enough!
       ! Allocating the structs and arrays is tedious and prone to segfaults
       ! and is unnecessary given the tiny memory usage of this data object
       integer(c_int) :: ntspec
       real(c_double):: dens(20)
       real(c_double):: temp(20)
       real(c_double):: fprim(20)
       real(c_double):: tprim(20)
       real(c_double):: nu(20)

       type(c_ptr) :: everything_struct_address

     end type gryfx_parameters_type

contains
  subroutine allocate_gryfx_zonal_arrays
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use nonlinear_terms, only: gryfx_zonal
    
    implicit none

    allocate(gryfx_zonal%NLdens_ky0(naky*ntheta0*2*ntgrid*nspec))
    allocate(gryfx_zonal%NLupar_ky0(naky*ntheta0*2*ntgrid*nspec))
    allocate(gryfx_zonal%NLtpar_ky0(naky*ntheta0*2*ntgrid*nspec))
    allocate(gryfx_zonal%NLtprp_ky0(naky*ntheta0*2*ntgrid*nspec))
    allocate(gryfx_zonal%NLqpar_ky0(naky*ntheta0*2*ntgrid*nspec))
    allocate(gryfx_zonal%NLqprp_ky0(naky*ntheta0*2*ntgrid*nspec))
  end subroutine allocate_gryfx_zonal_arrays    

  subroutine deallocate_gryfx_zonal_arrays
    use nonlinear_terms, only: gryfx_zonal
    
    implicit none

    deallocate(gryfx_zonal%NLdens_ky0)
    deallocate(gryfx_zonal%NLupar_ky0)
    deallocate(gryfx_zonal%NLtpar_ky0)
    deallocate(gryfx_zonal%NLtprp_ky0)
    deallocate(gryfx_zonal%NLqpar_ky0)
    deallocate(gryfx_zonal%NLqprp_ky0)
  end subroutine deallocate_gryfx_zonal_arrays    

  subroutine init_gs2_gryfx_c(strlen, run_name, mp_comm, gryfx_parameters) &
                            bind(c, name='init_gs2')
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: strlen
    character(kind=c_char), intent(in) :: run_name
    integer(c_int), intent(in) :: mp_comm
    type(gryfx_parameters_type), intent(in) :: gryfx_parameters
    call init_gs2_gryfx(strlen, run_name, mp_comm, gryfx_parameters)

  end subroutine init_gs2_gryfx_c

  subroutine init_gs2_gryfx(strlen, file_name, mp_comm, gryfx_parameters )
    use gs2_main, only: initialize_gs2
    use gs2_main, only: initialize_equations
    use gs2_main, only: initialize_diagnostics
    use gs2_main, only: prepare_miller_geometry_overrides
    use gs2_main, only: prepare_profiles_overrides
    use gs2_main, only: prepare_kt_grids_overrides
    use nonlinear_terms, only: gryfx_zonal
    use file_utils, only: run_name
    use mp, only: proc0

    implicit none
    integer, intent(in) :: strlen
    character (len=strlen), intent (in) :: file_name
    integer, intent(in) :: mp_comm
    type(gryfx_parameters_type), intent(in) :: gryfx_parameters

    !gryfx_zonal%on = .true.
    state%mp_comm_external = .true.
    state%mp_comm = mp_comm  

    state%run_name_external = .true.
    state%run_name = file_name

    gryfx_zonal%on = .true.
 

    call initialize_gs2(state)
    if(proc0) write(*,*) 'initialize_gs2 complete. gs2 thinks run_name is ', run_name
    !set overrides of naky, x0, y0, dt, vnewk here. 
    call prepare_kt_grids_overrides(state)
    if(proc0) write(*,*) 'setting overrides'
    state%init%kt_ov%override_gryfx = .true.
    state%init%kt_ov%gryfx = .true.
    !currently this is hard-coded in kt_grids.f90 and run_parameters.f90
    !just search for "GRYFX"
    if(.not. test_flag) then 
    call prepare_miller_geometry_overrides(state)
    call set_miller_geometry_overrides
    call prepare_profiles_overrides(state)
    call set_profiles_overrides
    endif
    call initialize_equations(state)
    if(proc0) write(*,*) 'initialize_equations complete.'
    call initialize_diagnostics(state)
    if(proc0) write(*,*) 'initialize_diagnostics complete.'
    
    call allocate_gryfx_zonal_arrays

contains

      subroutine set_miller_geometry_overrides
        state%init%mgeo_ov%override_rhoc = .true.
        state%init%mgeo_ov%rhoc = gryfx_parameters%rhoc

        state%init%mgeo_ov%override_qinp = .true.
        state%init%mgeo_ov%qinp = gryfx_parameters%qinp ! qval_gs2

        state%init%mgeo_ov%override_shat = .true.
        state%init%mgeo_ov%shat = gryfx_parameters%shat ! shat_gs2

        state%init%mgeo_ov%override_rgeo_lcfs = .true.
        state%init%mgeo_ov%rgeo_lcfs = gryfx_parameters%rgeo_lcfs

        state%init%mgeo_ov%override_rgeo_local = .true.
        state%init%mgeo_ov%rgeo_local = gryfx_parameters%rgeo_local

        state%init%mgeo_ov%override_akappa = .true.
        state%init%mgeo_ov%akappa = gryfx_parameters%akappa

        state%init%mgeo_ov%override_akappri = .true.
        state%init%mgeo_ov%akappri = gryfx_parameters%akappri

        state%init%mgeo_ov%override_tri = .true.
        state%init%mgeo_ov%tri = gryfx_parameters%tri

        state%init%mgeo_ov%override_tripri = .true.
        state%init%mgeo_ov%tripri = gryfx_parameters%tripri

        state%init%mgeo_ov%override_shift = .true.
        state%init%mgeo_ov%shift = gryfx_parameters%shift

        state%init%mgeo_ov%override_betaprim = .true.
        state%init%mgeo_ov%betaprim = gryfx_parameters%beta_prime_input
      end subroutine set_miller_geometry_overrides

      subroutine set_profiles_overrides
        use species, only: nspec
        integer :: is, isg, idens
        !write (*,*) 'SET_PROFILES_OVERRIDES'
        do is = 1,nspec
            isg = is
          state%init%prof_ov%override_dens(isg) = .true.
          state%init%prof_ov%dens(isg) = gryfx_parameters%dens(is) ! dens_ms(1,is)/dens_ms(1,idens)

          state%init%prof_ov%override_temp(isg) = .true.
          state%init%prof_ov%temp(isg) = gryfx_parameters%temp(is) ! temp_ms(1,is)/temp_ms(1,1)

          state%init%prof_ov%override_fprim(isg) = .true.
          state%init%prof_ov%fprim(isg) = gryfx_parameters%fprim(is) ! fprim_gs2(is)

          state%init%prof_ov%override_tprim(isg) = .true.
          state%init%prof_ov%tprim(isg) = gryfx_parameters%tprim(is) ! tprim_gs2(is)

          state%init%prof_ov%override_vnewk(isg) = .true.
          state%init%prof_ov%vnewk(isg) = gryfx_parameters%nu(is) ! nu_gs2(is)
        end do
        !gs2_state%init%prof_ov%override_g_exb = .true.
        !gs2_state%init%prof_ov%g_exb = gexb_gs2
        !gs2_state%init%prof_ov%override_mach = .true.
        !gs2_state%init%prof_ov%mach = mach_gs2

      end subroutine set_profiles_overrides

  end subroutine init_gs2_gryfx

  subroutine advance_gs2_gryfx_c(istep, dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, &
                               qpar_ky0, qprp_ky0, phi_ky0, first_half_step) &
                                    bind(c, name='advance_gs2')
    use iso_c_binding
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    implicit none
    complex (c_float_complex), intent (inout) :: & 
                                dens_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                upar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                tpar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                tprp_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                qpar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                qprp_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (out) :: phi_ky0(naky*ntheta0*2*ntgrid)
    integer(c_int), intent (in) :: first_half_step
    integer(c_int), intent (in) :: istep
    logical :: first_half_step_l

    first_half_step_l = .false.
    if (first_half_step.eq.1) first_half_step_l = .true.

    call advance_gs2_gryfx(istep, dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, &
                               qpar_ky0, qprp_ky0, phi_ky0, first_half_step_l)


  end subroutine advance_gs2_gryfx_c

  subroutine advance_gs2_gryfx(istep, dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, &
                               qpar_ky0, qprp_ky0, phi_ky0, first_half_step)

    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use nonlinear_terms, only: gryfx_zonal
    use gs2_main, only: evolve_equations
    use dist_fn, only: getmoms_gryfx
    use mp, only: proc0, broadcast

    implicit none
    complex*8, dimension (naky*ntheta0*2*ntgrid*nspec), intent (inout) :: &
                    dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, qpar_ky0, qprp_ky0
    complex*8, dimension (naky*ntheta0*2*ntgrid), intent (out) :: phi_ky0
    logical, intent (in) :: first_half_step
    integer, intent (in) :: istep

    !state%istep_end = istep
    state%print_times = .false.

    ! only proc0 knows values of arrays, so broadcast
    call broadcast(dens_ky0)
    call broadcast(upar_ky0)
    call broadcast(tpar_ky0)
    call broadcast(tprp_ky0)
    call broadcast(qpar_ky0)
    call broadcast(qprp_ky0)
    ! all procs know value of first_half_step already

    ! now that all procs know values, copy into gryfx_zonal structure
    gryfx_zonal%first_half_step = first_half_step 
    gryfx_zonal%NLdens_ky0 = dens_ky0
    gryfx_zonal%NLupar_ky0 = upar_ky0
    gryfx_zonal%NLtpar_ky0 = tpar_ky0
    gryfx_zonal%NLtprp_ky0 = tprp_ky0
    gryfx_zonal%NLqpar_ky0 = qpar_ky0
    gryfx_zonal%NLqprp_ky0 = qprp_ky0

    if(first_half_step) then
      state%dont_change_timestep = .true.
      state%skip_diagnostics = .true.
    else 
      state%dont_change_timestep = .false.
      state%skip_diagnostics = .false.
    endif
    ! only advance one timestep
    call evolve_equations(state, 1)
    ! getmoms_gryfx takes moments of g and puts results into arguments.
    ! since these are the same arguments that are passed into 
    ! advance_gs2_gryfx, we don't need to copy them like on the way in
    call getmoms_gryfx(dens_ky0, upar_ky0, tpar_ky0, &
                 tprp_ky0, qpar_ky0, qprp_ky0, phi_ky0)

  end subroutine advance_gs2_gryfx

  subroutine finish_gs2_gryfx_c() bind(c, name='finish_gs2')

    use iso_c_binding
    implicit none
    call finish_gs2_gryfx

  end subroutine finish_gs2_gryfx_c

  subroutine finish_gs2_gryfx
    use gs2_main, only: finalize_diagnostics
    use gs2_main, only: finalize_equations
    use gs2_main, only: finalize_gs2
    use nonlinear_terms, only: gryfx_zonal
    implicit none

    call deallocate_gryfx_zonal_arrays
    call finalize_diagnostics(state)
    call finalize_equations(state)
    state%print_times = .true.
    state%print_full_timers = .true.
    call finalize_gs2(state)
    gryfx_zonal%on = .false.

  end subroutine finish_gs2_gryfx

  subroutine getmoms_gryfx_c (dens_ky0, upar_ky0, tpar_ky0, tprp_ky0, &
                              qpar_ky0, qprp_ky0, phi_ky0) &
                                   bind(c, name='getmoms_gryfx')
    use dist_fn, only: getmoms_gryfx
    use theta_grid, only: ntgrid
    use kt_grids, only: ntheta0, naky
    use species, only: nspec
    use iso_c_binding
    implicit none
    complex (c_float_complex), intent (inout) :: & 
                                dens_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                upar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                tpar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                tprp_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                qpar_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (inout) :: & 
                                qprp_ky0(naky*ntheta0*2*ntgrid*nspec)
    complex (c_float_complex), intent (out) :: phi_ky0(naky*ntheta0*2*ntgrid)
    
    call getmoms_gryfx(dens_ky0, upar_ky0, tpar_ky0, &
            tprp_ky0, qpar_ky0, qprp_ky0, phi_ky0)
  end subroutine getmoms_gryfx_c

  subroutine broadcast_integer_c(a) bind(c, name='broadcast_integer')
    use iso_c_binding
    use mp, only : broadcast
    implicit none
    integer(c_int), intent(in out) :: a

    call broadcast(a) 
  end subroutine broadcast_integer_c
  
end module gs2_gryfx_zonal

