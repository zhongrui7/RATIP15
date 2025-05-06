program xauger
!
!-----------------------------------------------------------------------
! This is the MAIN program to control the calculation of 
!
! Calls: 
!-----------------------------------------------------------------------
   !
   use rabs_auger
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_grasp2k
   use rabs_print
   !
   implicit none
   !
   print *, "AUGER: Calculation of Auger rates, non-radiative lifetimes, "
   print *, " angular distribution parameters, and others (Fortran 95 version)"
   print *, " (C) Copyright by S Fritzsche and others, Kassel (2001)."
   print *, " "
   !
   ! Open the  .sum  file
   call file_open_formatted_stream(24,                                      &
      "Enter a file name for the  auger.sum  file:")
   !
   ! Determine all input information to control the calculation
   call auger_collect_input()
   !
   ! Get the eigenvectors of the initial and final states 
   call file_get_csl_list("Enter the name of the initial-state GRASP2K " // &
      "configuration symmetry list file:",asf_initial%csf_set) 
   call file_get_csl_list("Enter the name of the final-state GRASP2K "   // &
      "configuration symmetry list file:",asf_final%csf_set) 
   !
   ! Get the eigenvectors of the initial and final states 
   call file_get_mix("Enter the name of the initial-state GRASP2K "      // &
      "mixing coefficient file:",asf_initial)
   call file_get_mix("Enter the name of the final-state GRASP2K "        // &
      "mixing coefficient file:",asf_final)
   !
   ! Initialize and load the radial wavefunctions for the inital and 
   ! final atomic states
   !!x call auger_initialize_rwf_storage()   
   call initialize_rwf_storage(asf_initial, wave_initial)   
   call initialize_rwf_storage(asf_final,   wave_final)   
   call file_get_rwf(                                                       &
      "Enter the name of the initial-state Radial WaveFunction File:",      &
      asf_initial,wave_initial,.true.)
   call file_get_rwf(                                                       &
      "Enter the name of the final-state Radial WaveFunction File:",        &
      asf_final,wave_final,.true.)
   !
   ! Append a summary of the inputs to the  .sum  file
   call print_summary("AUGER",24,asf_initial,asf_final,wave_initial, wave_final)
   !
   ! Initialize the computations
   call auger_initialize()
   !
   ! Carry out the calculation of the transition probabilities and lifetimes
   call auger_calculate_amplitudes()
   !
   ! Output all results to the  .sum  file
   call auger_print_results(6)
   call auger_print_results(24)
   close(24)
   !
   if (auger_print_trn_file) then
      call auger_print_amplitudes(27)
      close(27)
   end if
   !
   close(62, status="delete")
   close(63, status="delete")
   close(99, status="delete")
   !
   print *, " "
   print *, "AUGER complete ... ."
   !
end program xauger

