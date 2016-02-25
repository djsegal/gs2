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
  use multigs2_library, only: initialize_ion_box
  use multigs2_library, only: initialize_electron_box
  use multigs2_library, only: run_box
  use multigs2_library, only: finalize_box
  use multigs2_library, only: finalize_multigs2
  use multigs2_library, only: multigs2_library_type
  use multigs2_library, only: communicate_omega



  implicit none

  type(multigs2_library_type) :: multigs2_obj
  logical :: test1,test2
  call initialize_multigs2( multigs2_obj) !initialize the program - read in run name, initialize MPI etc


  write(*,*) "for multigs2 this proc has iproc =" ,multigs2_obj%gs2_comm%iproc,&
    "nproc =",multigs2_obj%multigs2_comm%nproc,&    
    "proc0 =",multigs2_obj%multigs2_comm%proc0,&    
    "id =",multigs2_obj%multigs2_comm%id
  write(*,*) "for gs2 this proc has, iproc =",multigs2_obj%gs2_comm%iproc,&    
    "nproc =",multigs2_obj%gs2_comm%nproc,&    
    "proc0 =",multigs2_obj%gs2_comm%proc0,&    
    "id =",multigs2_obj%gs2_comm%id
 
 
  if (multigs2_obj%gs2_comm_id == 0) then ! here if we are sub communicator group 0 then we do an ionscale simulation
    write(*,*) "this proc has gs2_comm_id = ", multigs2_obj%gs2_comm_id,&
      "so runs an ionscale gs2"
    call initialize_ion_box( multigs2_obj)
  
  else if(multigs2_obj%gs2_comm_id > 0) then ! if we are not then we do an electronscale simulation
 
    write(*,*) "this proc has gs2_comm_id = ", multigs2_obj%gs2_comm_id,&
      "so runs an electronscale g2"
    call initialize_electron_box( multigs2_obj)
 
  end if 

  call run_box(multigs2_obj)

  call communicate_omega(multigs2_obj,test1,test2)
  
  call finalize_box(multigs2_obj) 

  call finalize_multigs2( multigs2_obj) ! finalize MPI 


end program multigs2
