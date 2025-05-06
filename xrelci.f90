program relativistic_CI
!
!-----------------------------------------------------------------------
! This is the MAIN program to control relativistic CI calculations by
! using atomic basis sets. It "uses" a few modules which are specific to
! relativistic CI calculation but also several modules from the RABS
! package. 
!
! Calls: 
!-----------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_hamiltonian
   use rabs_print
   use rabs_relci
   !
   implicit none
   !
   print *, "RELCI: Set-up and diagonalization of a relativistic CI matrix"//&
            " including" 
   print *, " Breit interactions and QED estimates (Fortran 95 version)"
   print *, " (C) Copyright by S Fritzsche and others, Kassel (2001)."
   print *, " "
   !
   ! Open the  .sum  file
   call file_open_formatted_stream(24,                               &
      "Enter a file name for the  relci.sum  file:")
   !
   ! Open, check, load data from, and close, the  .csl  file
   call file_get_csl_list(                                           &
      "Enter the name of the configuration symmetry list file:",     &
      asf_bound%csf_set)
   !
   ! Determine all input information to control the calculation
   call relci_collect_input()
   !
   ! Append a summary of the inputs to the  .sum  file
   call print_summary_c("RELCI",24,asf_bound,wave_bound)
   call relci_print_summary()
   !
   ! Set up the block structure of the Hamiltonian matrix into a derived type
   call hamiltonian_find_blocks(asf_bound%csf_set)
   !
   ! Set up and diagonalize the hamiltonian matrix for all (J,parity) blocks
   ! in turn 
   call relci_diagonalize_matrix(asf_bound,wave_bound)
   !
   ! Output all results to the  .sum  file
   !!x call relci_output_eigenvectors(asf_bound,wave_bound)
   call relci_print_results(6,asf_bound,wave_bound)
   call relci_print_results(24,asf_bound,wave_bound)
   call relci_output_eigenvectors(asf_bound,wave_bound)
   close(24)
   !
   close(62, status="delete")
   close(63, status="delete")
   close(99, status="delete")
   !
   print *, " "
   print *, "RELCI complete ... ."
   !
end program relativistic_CI

