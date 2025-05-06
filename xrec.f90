program xrec
!
!-----------------------------------------------------------------------
! This is the MAIN program to control the calculation of radiative electron
! capture (REC) cross sections and angular distribution parameters.
!
! Calls: 
!-----------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_grasp2k
   use rabs_nonorthonormal
   use rabs_print
   use rabs_rec
   !
   implicit none
   !
   print *, "REC: Calculation of radiative electron capture cross sections, "
   print *, " angular distribution parameters, and others (Fortran 95 version)"
   print *, " (C) Copyright by S Fritzsche and others, Kassel (2005)."
   print *, " "
   !
   ! Open the  .sum  file
   !!x call rec_open_summary()
   call file_open_formatted_stream(24, &
      "Enter a file name for the  rec.sum  file:")
   !
   ! Determine all input information to control the calculation
   call rec_collect_input()
   !
   ! Get the eigenvectors of the initial and final states 
   call file_get_csl_list("Enter the name of the initial-state GRASP2K "//& 
      "configuration symmetry list file:",asf_initial%csf_set)
   call file_get_csl_list("Enter the name of the final-state GRASP2K "  //& 
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
   !!x call rec_initialize_rwf_storage()   
   call initialize_rwf_storage(asf_initial, wave_initial)   
   call initialize_rwf_storage(asf_final,   wave_final)   
   call file_get_rwf(                                                  &
      "Enter the name of the initial-state Radial WaveFunction File:", &
                     asf_initial,wave_initial,.true.)
   call file_get_rwf(                                                  &
      "Enter the name of the final-state Radial WaveFunction File:",   &
                     asf_final,wave_final,.true.)
   !
   ! Append a summary of the inputs to the  .sum  file
   call print_summary("REC",24,asf_initial,asf_final,wave_initial, wave_final)
   !
   ! Initialize the computations
   call rec_initialize()
   !
   ! Carry out the calculation of the transition amplitudes and cross sections
   call rec_calculate_amplitudes()
   !
   ! Output all results to the  .sum  file
   call rec_print_results(6)
   call rec_print_results(24)
   close(24)
   !
   if (rec_print_chn_file) then
      call rec_print_amplitudes(27)
      close(27)
   end if
   !
   print *, " "
   print *, "REC complete ... ."
   !
end program xrec

