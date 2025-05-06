program xeinstein
!
!-----------------------------------------------------------------------
! This is the MAIN program to control the calculation Einstein coefficients
! but without including relaxation effects into the evaluation of the 
! matrix elements.
!
! Calls: 
!-----------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_grasp2k
   use rabs_einstein
   use rabs_nonorthonormal
   use rabs_print
   !
   implicit none
   !
   print *, "EINSTEIN: Calculation of Einstein A and B coefficients and "
   print *, " oscillator strengths (Fortran 95 version)"
   print *, " (C) Copyright by S Fritzsche, Kassel (2004)."
   print *, " "
   !
   ! Open the  .sum  file
   call file_open_formatted_stream(24,                                      &
      "Enter a file name for the  einstein.sum  file:")
   !
   ! Determine all input information to control the calculation
   call einstein_collect_input()
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
   !!x call einstein_initialize_rwf_storage()   
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
   call print_summary("EINSTEIN",24,asf_initial,asf_final, &
                                    wave_initial, wave_final)
   !
   ! Initialize the computations
   call einstein_initialize()
   !
   ! Carry out the calculation of the transition amplitudes and cross sections
   call einstein_calculate_amplitudes()
   !
   ! Output all results to the  .sum  file
   call einstein_print_results(6)
   call einstein_print_results(24)
   close(24)
   !
   if (einstein_write_transition_file) then
      call einstein_write_trn_file(27)
      close(27)
   end if
   !
   print *, " "
   print *, "EINSTEIN complete ... ."
   !
end program xeinstein

