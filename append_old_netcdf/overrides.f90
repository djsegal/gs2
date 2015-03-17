! DO NOT EDIT THIS FILE
! This file is automatically generated by overrides.rb

!> A module which defines the override types. These types
!! are used within the init object (which itself is contained
!! within the gs2_program_state object) to override values 
!! of the specified parameters (i.e. modify their values from
!! what is specified in the input file). The appropriate "prepare_..."
!! function from gs2_main must always be called before setting overrides.
module overrides
!> An object for overriding all or selected 
!! Miller geometry parameters.
  type miller_geometry_overrides_type
    !> DO NOT manually set the value of init.
    !! Nasty things may happen.
    logical :: init = .false.
    logical :: override_rhoc
    logical :: override_qinp
    logical :: override_shat
    logical :: override_rgeo_lcfs
    logical :: override_rgeo_local
    logical :: override_akappa
    logical :: override_akappri
    logical :: override_tri
    logical :: override_tripri
    logical :: override_shift
    logical :: override_betaprim
    real :: rhoc
    real :: qinp
    real :: shat
    real :: rgeo_lcfs
    real :: rgeo_local
    real :: akappa
    real :: akappri
    real :: tri
    real :: tripri
    real :: shift
    real :: betaprim
    
  end type miller_geometry_overrides_type


!> An object for overriding all or selected
!! profile parameters, for example species
!! temps, densities or gradients or the flow gradient or mach
!! number. Note that all species parameters are arrays of 
!! size nspec and you must set the override switches 
!! individually for each species.
  type profiles_overrides_type
    !> DO NOT manually set the value of init.
    !! Nasty things may happen.
    logical :: init = .false.
    logical, dimension(:), pointer :: override_dens
    logical, dimension(:), pointer :: override_temp
    logical, dimension(:), pointer :: override_tprim
    logical, dimension(:), pointer :: override_fprim
    logical, dimension(:), pointer :: override_vnewk
    logical :: override_g_exb
    logical :: override_mach
    real, dimension(:), pointer :: dens
    real, dimension(:), pointer :: temp
    real, dimension(:), pointer :: tprim
    real, dimension(:), pointer :: fprim
    real, dimension(:), pointer :: vnewk
    real :: g_exb
    real :: mach
    
  end type profiles_overrides_type



!> A type for containing overrides to the processor layout
!! and optimisation flags for gs2. 
  type optimisations_overrides_type
    !> DO NOT manually set the value of init.
    !! Nasty things may happen.
    logical :: init = .false.
    logical :: override_nproc
    logical :: override_layout
    integer :: nproc
    character(len=6) :: layout
    integer :: old_comm
  end type optimisations_overrides_type



!> A type for storing overrides of the intial
!! values of the fields and distribution function.
!! This override is different to all the others, 
!! because in order to minimise memory usage and disk writes,
!! this override is used internally during the simulation, and
!! its values can thus change over the course of the simulation.
!! In contrast, no other overrides are modified by running gs2.
!! Also, depending on the value of in_memory, the override
!! values will either be taken from the the arrays within
!! the object, or from the restart files. If you want
!! to externally modify the initial field and dist fn values,
!! you need to use in_memory = .true. If you just want to 
!! use this override to allow you to reinitialise the equations
!! and start from the same values, you can use either memory
!! or restart files. If you want to want to change the number
!! of processors and then reinitialise and then use this override
!! you must use in_memory = .false., because currently the memory
!! is allocated on a processor by processor basis. Changing
!! grid sizes and then using this override is not supported. 
!! This one is too complicated to generate 
!! automatically
type initial_values_overrides_type
  !> True if the object has been initialized.
  logical :: init = .false.
  !> If true, override values are read from the
  !! arrays in this object. If not, they are read
  !! from the restart files. The value of in_memory
  !! should not be changed without reinitializing this  
  !! object (doing so is an excellent way of generating
  !! segmentation faults).
  logical :: in_memory = .true.
  !> Whether to override initial values or not,
  !! i.e., whether or not this override is switched on.
  !! If it is switched on, initial values will be determined
  !! by the values in the arrays or the restart files,
  !! depending on the value of in_memory. If false, 
  !! initial values will be determined by the gs2 input file
  !! (note that of course, this can result in initial values
  !! being taken from the input files).
  logical :: override = .false.
  logical :: force_maxwell_reinit = .true.
  complex, dimension (:,:,:), pointer :: phi
  complex, dimension (:,:,:), pointer :: apar
  complex, dimension (:,:,:), pointer :: bpar
  complex, dimension (:,:,:), pointer :: g
end type initial_values_overrides_type

contains
  subroutine init_miller_geometry_overrides(overrides_obj)
    type(miller_geometry_overrides_type), intent(inout) :: overrides_obj
    if (overrides_obj%init) return 
    overrides_obj%init = .true.
    overrides_obj%override_rhoc = .false.
    overrides_obj%override_qinp = .false.
    overrides_obj%override_shat = .false.
    overrides_obj%override_rgeo_lcfs = .false.
    overrides_obj%override_rgeo_local = .false.
    overrides_obj%override_akappa = .false.
    overrides_obj%override_akappri = .false.
    overrides_obj%override_tri = .false.
    overrides_obj%override_tripri = .false.
    overrides_obj%override_shift = .false.
    overrides_obj%override_betaprim = .false.
  end subroutine init_miller_geometry_overrides


  subroutine finish_miller_geometry_overrides(overrides_obj)
    type(miller_geometry_overrides_type), intent(inout) :: overrides_obj
    if (.not. overrides_obj%init) then
      write (*,*) "ERROR: Called finish_miller_geometry_overrides on an uninitialized object"
      return
    end if
    overrides_obj%init = .false.
    overrides_obj%override_rhoc = .false.
    overrides_obj%override_qinp = .false.
    overrides_obj%override_shat = .false.
    overrides_obj%override_rgeo_lcfs = .false.
    overrides_obj%override_rgeo_local = .false.
    overrides_obj%override_akappa = .false.
    overrides_obj%override_akappri = .false.
    overrides_obj%override_tri = .false.
    overrides_obj%override_tripri = .false.
    overrides_obj%override_shift = .false.
    overrides_obj%override_betaprim = .false.
  end subroutine finish_miller_geometry_overrides


  subroutine init_profiles_overrides(overrides_obj, nspec)
    integer, intent(in) :: nspec
    type(profiles_overrides_type), intent(inout) :: overrides_obj
    if (overrides_obj%init) return 
    overrides_obj%init = .true.
    allocate(overrides_obj%override_dens(nspec), overrides_obj%dens(nspec))
    overrides_obj%override_dens = .false.
    allocate(overrides_obj%override_temp(nspec), overrides_obj%temp(nspec))
    overrides_obj%override_temp = .false.
    allocate(overrides_obj%override_tprim(nspec), overrides_obj%tprim(nspec))
    overrides_obj%override_tprim = .false.
    allocate(overrides_obj%override_fprim(nspec), overrides_obj%fprim(nspec))
    overrides_obj%override_fprim = .false.
    allocate(overrides_obj%override_vnewk(nspec), overrides_obj%vnewk(nspec))
    overrides_obj%override_vnewk = .false.
    overrides_obj%override_g_exb = .false.
    overrides_obj%override_mach = .false.
  end subroutine init_profiles_overrides


  subroutine finish_profiles_overrides(overrides_obj)
    type(profiles_overrides_type), intent(inout) :: overrides_obj
    if (.not. overrides_obj%init) then
      write (*,*) "ERROR: Called finish_profiles_overrides on an uninitialized object"
      return
    end if
    overrides_obj%init = .false.
    overrides_obj%override_dens = .false.
    deallocate(overrides_obj%override_dens, overrides_obj%dens)
    overrides_obj%override_temp = .false.
    deallocate(overrides_obj%override_temp, overrides_obj%temp)
    overrides_obj%override_tprim = .false.
    deallocate(overrides_obj%override_tprim, overrides_obj%tprim)
    overrides_obj%override_fprim = .false.
    deallocate(overrides_obj%override_fprim, overrides_obj%fprim)
    overrides_obj%override_vnewk = .false.
    deallocate(overrides_obj%override_vnewk, overrides_obj%vnewk)
    overrides_obj%override_g_exb = .false.
    overrides_obj%override_mach = .false.
  end subroutine finish_profiles_overrides


  subroutine init_optimisations_overrides(overrides_obj)
    type(optimisations_overrides_type), intent(inout) :: overrides_obj
    if (overrides_obj%init) return 
    overrides_obj%init = .true.
    overrides_obj%override_nproc = .false.
    overrides_obj%override_layout = .false.
  end subroutine init_optimisations_overrides


  subroutine finish_optimisations_overrides(overrides_obj)
    type(optimisations_overrides_type), intent(inout) :: overrides_obj
    if (.not. overrides_obj%init) then
      write (*,*) "ERROR: Called finish_optimisations_overrides on an uninitialized object"
      return
    end if
    overrides_obj%init = .false.
    overrides_obj%override_nproc = .false.
    overrides_obj%override_layout = .false.
  end subroutine finish_optimisations_overrides



!> Warning: You can't change the value of overrides%force_maxwell_reinit 
!! or overrides%in_memory after calling this function
subroutine init_initial_values_overrides(overrides_obj, ntgrid, ntheta0, naky, g_llim, g_ulim, force_maxwell_reinit, in_memory)
  use file_utils, only: error_unit
  implicit none
  type(initial_values_overrides_type), intent(inout) :: overrides_obj
  integer, intent(in) :: ntgrid, ntheta0, naky, g_llim, g_ulim
  logical, intent(in) :: force_maxwell_reinit, in_memory
  integer :: iostat
  if (overrides_obj%init) return
  overrides_obj%init = .true.
  overrides_obj%override = .false.
  !overrides_obj%override_phi = .false.
  !overrides_obj%override_apar = .false.
  !overrides_obj%override_bpar = .false.
  !overrides_obj%override_g = .false.
  overrides_obj%force_maxwell_reinit = force_maxwell_reinit
  overrides_obj%in_memory = in_memory

  write (error_unit(), *) "INFO: changing force_maxwell_reinit or in_memory &
    & after calling initial_values_overrides_type will almost certainly cause &
    & segmentation faults."
  if (overrides_obj%in_memory) then 
    allocate(overrides_obj%g(-ntgrid:ntgrid,2,g_llim:g_ulim), stat=iostat)
    if (overrides_obj%force_maxwell_reinit) then 
      if (iostat.eq.0) allocate(overrides_obj%phi(-ntgrid:ntgrid,ntheta0,naky), stat=iostat)
      if (iostat.eq.0) allocate(overrides_obj%apar(-ntgrid:ntgrid,ntheta0,naky), stat=iostat)
      if (iostat.eq.0) allocate(overrides_obj%bpar(-ntgrid:ntgrid,ntheta0,naky), stat=iostat)
    end if
    if (iostat.ne.0) then
      overrides_obj%in_memory = .false.
      write(error_unit(),*) "WARNING: could not allocate memory for initial_values_overrides. Only restart from file possible (manual setting of initial values not possible)"
      if (associated(overrides_obj%g)) deallocate(overrides_obj%g)
      if (associated(overrides_obj%phi)) deallocate(overrides_obj%phi)
      if (associated(overrides_obj%apar)) deallocate(overrides_obj%apar)
      if (associated(overrides_obj%bpar)) deallocate(overrides_obj%bpar)
    end if
  end if
end subroutine init_initial_values_overrides

subroutine finish_initial_values_overrides(overrides_obj)
  implicit none
  type(initial_values_overrides_type), intent(inout) :: overrides_obj
  if (.not. overrides_obj%init) then
    write (*,*) "WARNING: Called finish_initial_values_overrides on an uninitialized object"
    return 
  end if
  overrides_obj%init = .false.
  overrides_obj%override = .false.
  !overrides%override_phi = .false.
  !overrides%override_apar = .false.
  !overrides%override_bpar = .false.
  !overrides%override_g = .false.
  overrides_obj%force_maxwell_reinit = .true.
  if (overrides_obj%in_memory) then 
    deallocate(overrides_obj%g)
    if (overrides_obj%force_maxwell_reinit) then 
      deallocate(overrides_obj%phi)
      deallocate(overrides_obj%apar)
      deallocate(overrides_obj%bpar)
    end if
  end if
end subroutine finish_initial_values_overrides

end module overrides
