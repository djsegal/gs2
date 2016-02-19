!This is the main file of multigs2, a program which simulates an electron flux tube at each ion spatial point within the ion flux tube of a standard gs2 run
!
!
! authors: Michael Hardman (mhardman1@users.sourceforge.net) 
! this software is available under the MIT license 

!> multigs2, a program which
!! runs separate ion and electron scale simulations 
!! and communicates information between them.
!!
!! This is free software released under the MIT license.
!! Written by:
!!             Michael Hardman (mhardman1@users.sourceforge.net)
!!             Michael Barnes ()
!!             Edmund Highcock (edmundhighcock@users.sourceforge.net)
program multigs2

  use multigs2_library, only: initialize_multigs2
  use multigs2_library, only: run_ion_box
  use multigs2_library, only: run_electron_box
  use multigs2_library, only: finalize_multigs2
  use multigs2_library, only: multigs2_library_type
  
  implicit none

  type(multigs2_library_type) :: multigs2_obj

  call initialize_multigs2( multigs2_obj) !initialize the program - read in input files, initialize MPI etc
 
  if (multigs2_obj%sub_comm == 0) then ! here if we are sub communicator group 0 then we do an ionscale simulation

    call run_ion_box( multigs2_obj)
  
  else if(multigs2_obj%sub_comm > 0) then ! if we are not then we do an electronscale simulation
 
  call run_electron_box( multigs2_obj)
 
  end if 
 
  call finalize_multigs2( multigs2_obj) ! finalize MPI 


end program multigs2
