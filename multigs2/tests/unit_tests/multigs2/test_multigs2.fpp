
!> A program that tests the multigs2_library module. It  runs 
!! a linear ion simulation followed by 4 concurrent linear electron simulations.
!!
!! This is free software released under the MIT license
!!   Written by: Edmund Highcock (edmundhighcock@users.sourceforge.net)
!!             Michael Hardman (mhardman1@users.sourceforge.net)

program test_multigs2
  use unit_tests
  use mp, only: init_mp, finish_mp, mp_comm
  use multigs2_library, only: initialize_multigs2
  use multigs2_library, only: multigs2_library_type
  use multigs2_library, only: run_ion_box
  implicit none
  real :: eps
  type(multigs2_library_type) :: multigs2_obj

  eps = 1.0e-7
  if (precision(eps).lt. 11) eps = eps * 1000.0


  call init_mp

  call announce_module_test("multigs2")

  multigs2_obj%mp_comm_external = .true.
  multigs2_obj%mp_comm = mp_comm

  call initialize_multigs2(multigs2_obj)

  call announce_test("run_name")
  call process_test(&
    agrees_with(trim(multigs2_obj%run_name), "test_multigs2"),  &
    "run_name")

  call run_ion_box(multigs2_obj)

  call close_module_test("multigs2")

  call finish_mp

contains



end program test_multigs2
