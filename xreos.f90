program xreos
!
!-----------------------------------------------------------------------
! This is the MAIN program to control the calculation of relaxed-orbital
! transition probabilities and lifetimes for atomic states which are 
! represented in a determinant basis. Compared with previous versions,
! this program does not read the .csl files for the initial and final
! atomic states anymore. All necessary information are now extracted from
! the corresponding .xpn files. Moreover, this revised version now applies 
! a compact storage of the initial and final state determinants as well as 
! a simplified calculation if the user 'assumes' orthogonality for the
! initial and final state orbital sets. A further acceleration for
! most computations is achieved by calculating all transitions probabilities
! in 'parallel', i.e. each pair of determinants is evaluated only once 
! independent of the number of transitions.
! This module "uses" a few other modules which are specific to the xreos
! program but also modules from the RABS package. 
!
! Calls: reos_collect_input(), reos_get_rwf_final(), reos_get_rwf_initial(),
!        reos_get_xpn_final(), reos_get_xpn_initial(),
!        reos_initialize_rwf_storage(),
!        reos_print_summary(), reos_set_debug(), reos_write_trn_file().
!-----------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_csl
   use rabs_grasp2k
   use rabs_reos
   !
   implicit none
   !
   print *, "REOS: Calculation of relaxed-orbital transitions probabilities"
   print *, " and lifetimes within a determinant basis (Fortran 95 version)"
   print *, " (C) Copyright by S Fritzsche and others, CPC 124 (2000) 340."
   print *, " "
   !
   ! Open the  .sum  file
   call file_open_formatted_stream(24,                                    &
      "Enter a file name for the  reos.sum  file:")
   !
   ! Determine all input information to control the calculation
   call reos_collect_input()
   !
   if (.not.reos_apply_restart) then
      !
      ! Get the eigenvectors of the initial and final states in their
      ! determinant representation
      call file_get_xpn(                                                  &
         "Enter the name of the initial-state CESD .xpn expansion file",  &
	 asf_initial)
      call file_get_xpn(                                                  &
         "Enter the name of the final-state CESD .xpn expansion file",    &
	 asf_final)
      !
      ! Initialize and load the radial wavefunctions for the inital and 
      ! final atomic states
      call reos_initialize_rwf_storage()   
      call file_get_rwf_det(                                              &
         "Enter the name of the initial-state Radial WaveFunction File:", &
	 wave_initial)
      call file_get_rwf_det(                                              &
         "Enter the name of the final-state Radial WaveFunction File:",   &
	 wave_final)
   end if
   !
   ! Append a summary of the inputs to the  .sum  file
   call reos_print_summary()
   !
   ! Carry out the calculation of the transition probabilities and lifetimes
   call reos_calculate_amplitudes()
   !
   ! Output all results to the  .sum  file
   call reos_print_results(6)
   call reos_print_results(24)
   close(24)
   !
   if (reos_write_transition_file) then
      call reos_write_trn_file()
   end if
   !
   print *, " "
   print *, "REOS complete ... ."
   !
end program xreos

