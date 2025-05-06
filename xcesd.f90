program xcesd
!
!-----------------------------------------------------------------------
! This is the MAIN program to control the (complete) expansion of the 
! .csl list files into a determinant basis. This module "uses" several
! other modules which belong to the RABS and/or the RATIP package. 
! Such an expansion applies the cfp-expansion of antisymmetric subshell
! states into uncoupled product states. 
! Two different modes are supported, the expansion of a given set of
! atomic states but also of the whole configuration basis. However, in
! latter case only .csl lists of a small to moderate size is possible 
! since a full rotation of a large subspace of the Hilbert space may 
! require a large amount of memory.
!
! Calls: cesd_collect_input(), get_csl_list(),
!        cesd_open_xpansion(), cesd_print_summary(), cesd_perform_xpansion(),
!        cesd_set_debug(), set_eigenvectors_formatted(),
!        set_eigenvectors_unformatted().
!-----------------------------------------------------------------------
   !
   use rabs_cesd
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_print
   !!x use rabs_utilities
   !
   implicit none
   logical         :: asf_mode
   type(asf_basis) :: asf_set
   !
   print *, "CESD: Complete expansion of jj-coupled symmetry functions"
   print *, " into Slater determinants (Fortran 95 version)"
   print *, " (C) Copyright by S Fritzsche and others, CPC 124 (2000) 353."
   print *, " "
   !
   ! Open the  .sum  file
   call file_open_formatted_stream(24,                             &
      "Enter a file name for the  cesd.sum  file:")
   !
   ! Open, check, load data from, and close, the  .csl  file
   call file_get_csl_list(                                         &
      "Enter the name of the configuration symmetry list file:",   &
      asf_set%csf_set)
   !
   ! Determine all input information to control the calculation
   call cesd_collect_input(asf_mode)
   !
   ! Get eigenvectors of the Atomic State Functions
   if (asf_mode) then
      call file_get_mix("Enter the name of the GRASP2K mixing ",asf_set)
   else
      asf_set%noasf = asf_set%csf_set%nocsf
   end if
   !
   ! Append a summary of the inputs to the  .sum  file
   call print_summary_b("CESD",24,asf_set)
   !
   ! Open the  .xpn  file
   call file_open_formatted_stream(25,                            &
      "Enter a name for the CESD eXPaNsion  .xpn  output File"//  &
      " that is to be created:",cesd_xpn_file)                                   
   !
   ! Perform and output the expansion into Slater determinants
   call cesd_perform_xpansion(asf_mode,asf_set)
   !
   print *, " "
   print *, "CESD complete ... ."
   !
end program xcesd

