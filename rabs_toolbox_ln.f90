module rabs_toolbox_ln
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains procedures which help carry out a number of
! 'utility' tasks within the RATIP environment. 
!-----------------------------------------------------------------------
   !
   use rabs_cesd
   use rabs_constant
   use rabs_file_handling
   use rabs_input_dialog
   use rabs_multipole
   use rabs_print
   use rabs_toolbox_aux
   use rabs_xl
   implicit none
   !
   public  :: toolbox_me_breit0
                 ! Calculates effective Breit0 interaction strengthes for given
		 ! atomic orbitals.
   public  :: toolbox_me_breit0_strength
                 ! Enters quantum numbers and evaluates the interaction
		 ! strengthes.
   public  :: toolbox_me_coulomb
                 ! Calculates effective Coulomb interaction strengthes for given
		 ! atomic orbitals.
   public  :: toolbox_me_coulomb_strength
                 ! Enters quantum numbers and evaluates the interaction
		 ! strengthes.
   public  :: toolbox_me_multipole
                 ! Enters quantum numbers and evaluates the one-electron
		 ! multipole matrix element.
   public  :: toolbox_merge_csl_diffcores
                 ! Merges two .csl lists together.
   public  :: toolbox_merge_csl_lists
                 ! Merges two .csl lists together.
   public  :: toolbox_merge_orbitals_formatted
                 ! Recombines the radial orbitals from two or more (formatted)
                 ! .out files into a new one.
   public  :: toolbox_merge_check_csf_set
                 ! Checks that two CSF lists have proper oder of orbitals and
                 ! creates an unified orbital sequence.
   public  :: toolbox_merge_set_reduced_list
                 ! Packs the information about the CSFs into a more compact 
                 ! format to accelerate the comparison of different CSF.
   public  :: toolbox_merge_two_list
                 ! Compares all (reduced) CSF from list_b with those from 
		 ! list_a and 'eliminates' the doubly defined CSF in list_b.
   public  :: toolbox_merge_write_list
                 ! Writes the 'accepted' CSF from the tow given csl list
                 ! to a new .csl file.
   public  :: toolbox_modify_energies
                 ! Modifies and overwrites the total energies of a given 
		 ! .mix file. 
   public  :: toolbox_nuclear_fermi
                 ! Calculate the Fermi distribution parameters a, c and N 
		 ! as well as the distribution and potential.
   public  :: toolbox_nuclear_radius
                 ! Calculate the nuclear radius from the mass number.
                 !
   type(csf_basis), public :: csf_set_a, csf_set_b, csf_set_c
   !
   !
   ! Define data structures for the 'merge' part
   ! Define number of integer for reduced storage
   integer :: merge_csf_counter, merge_noint, merge_nobit = bit_size(1)
   !
   type, public :: reduced_csf
      logical   :: append
      integer, dimension(:), pointer :: red_occ, red_X
      integer, dimension(:), pointer :: red_shell
   end type reduced_csf
   !
   integer, public :: no_list_a, no_list_b
   type(reduced_csf), dimension(:), pointer :: list_a, list_b
   !
   !
   ! Define an internal structure which stores all information about the
   ! photoionization lines of a given step
   type :: toolbox_p_channel
      integer          :: kappa, totalJ   
      character(len=1) :: parity 
      character(len=2) :: multipole  
      character(len=9) :: gauge 
      real(kind=dp)    :: phase, amplitude_re  
      complex(kind=dp) :: amplitude   
   end type toolbox_p_channel
   !
   type :: toolbox_p_line
      integer          :: asfi, asff, level_i, level_f, totalJ_i, totalJ_f
      integer          :: No_channels
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: p_energy, e_energy, cs_coulomb, cs_babushkin 
      complex(kind=dp) :: beta_b, beta_c, xi_b, xi_c, eta_b, eta_c, &
                          zeta_b, zeta_c, alignment_b, alignment_c
      type(toolbox_p_channel), dimension(:), pointer :: channel
   end type toolbox_p_line
   !
   type(toolbox_p_line), dimension(:), allocatable :: line_step1, line_step2
   !
   type :: toolbox_p_single_step
      integer          :: asfi, asff, totalJ_i, totalJ_f, tr
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: p_energy, e_energy
      real(kind=dp)    :: cs_b, cs_c, beta_b, beta_c, alignment_b, alignment_c 
      ! First-step properties
      real(kind=dp)    :: b2_b, b2_c, b4_b, b4_c
      ! Second-step properties
      real(kind=dp)    :: a0_b, a0_c, a2_b, a2_c, a4_b, a4_c
      real(kind=dp)    :: beta_2_TPDI_b,beta_2_TPDI_c,beta_4_TPDI_b,beta_4_TPDI_c
      real(kind=dp)    :: cs_0if_b, cs_0if_c
   end type toolbox_p_single_step
   !
   type :: toolbox_p_cascade
      ! Cascade properties
      real(kind=dp)    :: beta2_0i_f_b, beta2_0i_f_c, &
                          beta4_0i_f_b, beta4_0i_f_c 
   end type toolbox_p_cascade
   !
   !
   type :: toolbox_p_2step
      integer          :: ni, nn, nf
      real(kind=dp)    :: energy
      type(toolbox_p_single_step), dimension(5,5)   :: step1, step2
      type(toolbox_p_cascade),     dimension(5,5,5) :: cascade12
   end type toolbox_p_2step
   !
   type(toolbox_p_2step) :: photon
   !
   !
   ! Define an internal structure to deal with the coherence transfer and
   ! the calculation of angular correlation functions
   type :: toolbox_coherence_B
       real, dimension(0:6,0:6,0:6) :: value
   end type toolbox_coherence_B
   !
   type :: toolbox_coherence
      integer               :: nb_max = 6, nc_max = 10
      integer               :: ni, nn, nf, asf_i, asf_f
      integer               :: totalJ_i, totalJ_f
      integer, dimension(5) :: asf_n, totalJ_n, tr_in, tr_nf
      character(len=1)      :: parity_i, parity_f 
      !
      character(len=1), dimension(5)                :: parity_n
      real(kind=dp), dimension(5,5)                 :: h
      type(toolbox_coherence_B), dimension(0:5,0:5) :: B, Bbar
      logical, dimension(0:10,0:10,0:10)            :: need_C
      real(kind=dp), dimension(0:10,0:10,0:10)      :: C
      !
      integer                         :: no_angles
      real(kind=dp)                   :: theta_a, phi_a, theta_b, phi_b
      real(kind=dp), dimension(50)    :: theta, phi
      real(kind=dp), dimension(50,50) :: W
   end type toolbox_coherence
   !
   type(toolbox_coherence), save :: coherence
   !
contains
   !
   subroutine toolbox_ln()
   !--------------------------------------------------------------------
   ! Calls: 
   !--------------------------------------------------------------------
      !
      print *, "**************************"
      print *, "*** Not yet implmented ***"
      print *, "**************************"
      !
   end subroutine toolbox_ln
   !
   !
   subroutine toolbox_me_breit0()
   !--------------------------------------------------------------------
   ! Calculates (frequency-independent) effective Breit interaction 
   ! strengthes for a given set of atomic orbitals.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      logical :: yes
      !
      ! Define radial grid parameter and generate the grid
      rnt_grasp2k = 2.0e-6_dp
      h_grasp2k   = 5.0e-2_dp
      n_grasp2k   = 390
      hp_grasp2k  = zero
      !
      ! Default speed of light
      c = c_vacuum
      !
      ! Modify the radial grid parameters
      ! print *, "The default radial grid parameters for this case are:"
      ! print *, " rnt = ",rnt_grasp2k,";"
      ! print *, " h   = ",h_grasp2k,  ";"
      ! print *, " hp  = ",hp_grasp2k, ";"
      ! print *, " n   = ",n_grasp2k,  ";"
      ! print *, " revise these values ?"
      ! yes = get_yes_stream()
      ! if (yes) then
      !    print *, "Enter rnt:"
      !    read *, rnt_grasp2k
      !    print *, "Enter h:"
      !    read *, h_grasp2k
      !    print *, "enter hp:"
      !    read *, hp_grasp2k
      !    print *, "enter n:"
      !    read *, n_grasp2k
      ! end if
      !
      ! accy is an estimate of the accuracy of the numerical procedures
      accy_grasp2k = h_grasp2k**6
      !
      ! Set up the coefficients for the numerical procedures
      call setqic_grasp2k()
      !
      ! Generate the radial grid and all associated arrays
      call radgrd_grasp2k()
      !
      ! Get the occupied shells of the GRASP92 states
      call file_get_csl_list( &
         "Enter the name of the configuration symmetry list file:", &
	 asf_set%csf_set)
      !
      ! Initialize and load the radial wavefunctions 
      !!x call toolbox_initialize_rwf_storage(.true.,.false.,.false.)
      call initialize_rwf_storage(asf_set, wave)
      call file_get_rwf("Enter the name of the Radial WaveFunction File:", &
	 asf_set,wave,.true.)
      !
      ! Enter the quantum numbers and evaluate the Coulomb strength
      call toolbox_me_breit0_strength()
      !
      !
   end subroutine toolbox_me_breit0
   !
   !
   subroutine toolbox_me_breit0_strength()
   !--------------------------------------------------------------------
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer       :: ia, ib, ic, id, na, nb, nc, nd, &
                       kapa, kapb, kapc, kapd, L
      real(kind=dp) :: wx
      !
    1 print *, "Enter the (9!) quantum numbers na, kapa, nb, kapb, "// &
               "nc, kapc, nd, kapd, L:"
      read(*,*,err=5)  na, kapa, nb, kapb, nc, kapc, nd, kapd, L
      !
      ia = 0;   ib = 0;   ic = 0;   id = 0;   
      do  i = 1,wave%number_of_rwf
         if (wave%rwf(i)%orbital%n     == na   .and.  &
	     wave%rwf(i)%orbital%kappa == kapa)  ia = i
         if (wave%rwf(i)%orbital%n     == nb   .and.  &
	     wave%rwf(i)%orbital%kappa == kapb)  ib = i
         if (wave%rwf(i)%orbital%n     == nc   .and.  &
	     wave%rwf(i)%orbital%kappa == kapc)  ic = i
         if (wave%rwf(i)%orbital%n     == nd   .and.  &
	     wave%rwf(i)%orbital%kappa == kapd)  id = i
      end do
      !
      if (ia == 0  .or.  ib == 0  .or.  ic == 0  .or.  id == 0) then
         print *, "Radial functions are not found for all the given orbitals:" 
	 print *, "   ",orbital_name(na,kapa), orbital_name(nb,kapb), &
                        orbital_name(nc,kapc), orbital_name(nd,kapd)
         print *, "Input error; redo ..."
         goto 1
      end if
      !
      wx = XL_Breit0_strength_grasp2k(L,wave%rwf(ia),wave%rwf(ib),   &
                                        wave%rwf(ic),wave%rwf(id))
      !
      print *, " "
      print *, "X^" // trim(angular_momentum_string(L+L,4)) // "_Breit0 ("//  &
               orbital_name(na,kapa) // "," // orbital_name(nb,kapb) //";"//  &
               orbital_name(nc,kapc) // "," // orbital_name(nd,kapd) //       &
	       ") = ",wx
      !       
      print *, "Continue (y/n) ?";   if (get_yes_stream()) goto 1
      return
      !
    5 print *, "Input error; redo ..."
      goto 1
      !     
   end subroutine toolbox_me_breit0_strength
   !
   !
   subroutine toolbox_me_coulomb()
   !--------------------------------------------------------------------
   ! Calculates effective Coulomb interaction strengthes for a given set
   ! of atomic orbitals.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      logical :: yes
      !
      ! Define radial grid parameter and generate the grid
      rnt_grasp2k = 2.0e-6_dp
      h_grasp2k   = 5.0e-2_dp
      n_grasp2k   = 390
      hp_grasp2k  = zero
      !
      ! Default speed of light
      c = c_vacuum
      !
      ! Modify the radial grid parameters
      ! print *, "The default radial grid parameters for this case are:"
      ! print *, " rnt = ",rnt_grasp2k,";"
      ! print *, " h   = ",h_grasp2k,  ";"
      ! print *, " hp  = ",hp_grasp2k, ";"
      ! print *, " n   = ",n_grasp2k,  ";"
      ! print *, " revise these values ?"
      ! yes = get_yes_stream()
      ! if (yes) then
      !    print *, "Enter rnt:"
      !    read *, rnt_grasp2k
      !    print *, "Enter h:"
      !    read *, h_grasp2k
      !    print *, "enter hp:"
      !    read *, hp_grasp2k
      !    print *, "enter n:"
      !    read *, n_grasp2k
      ! end if
      !
      ! accy is an estimate of the accuracy of the numerical procedures
      accy_grasp2k = h_grasp2k**6
      !
      ! Set up the coefficients for the numerical procedures
      call setqic_grasp2k()
      !
      ! Generate the radial grid and all associated arrays
      call radgrd_grasp2k()
      !
      ! Get the occupied shells of the GRASP92 states
      call file_get_csl_list( &
         "Enter the name of the configuration symmetry list file:", &
	 asf_set%csf_set)
      !
      ! Initialize and load the radial wavefunctions 
      !!x call toolbox_initialize_rwf_storage(.true.,.false.,.false.)
      call initialize_rwf_storage(asf_set, wave)
      call file_get_rwf("Enter the name of the Radial WaveFunction File:", &
	 asf_set,wave,.true.)
      !
      ! Enter the quantum numbers and evaluate the Coulomb strength
      call toolbox_me_coulomb_strength()
      !
      !
   end subroutine toolbox_me_coulomb
   !
   !
   subroutine toolbox_me_coulomb_strength()
   !--------------------------------------------------------------------
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer       :: ia, ib, ic, id, na, nb, nc, nd, &
                       kapa, kapb, kapc, kapd, L
      real(kind=dp) :: wx
      !
    1 print *, "Enter the (9!) quantum numbers na, kapa, nb, kapb, "// &
               "nc, kapc, nd, kapd, L:"
      read(*,*,err=5)  na, kapa, nb, kapb, nc, kapc, nd, kapd, L
      !
      ia = 0;   ib = 0;   ic = 0;   id = 0;   
      do  i = 1,wave%number_of_rwf
         if (wave%rwf(i)%orbital%n     == na   .and.  &
	     wave%rwf(i)%orbital%kappa == kapa)  ia = i
         if (wave%rwf(i)%orbital%n     == nb   .and.  &
	     wave%rwf(i)%orbital%kappa == kapb)  ib = i
         if (wave%rwf(i)%orbital%n     == nc   .and.  &
	     wave%rwf(i)%orbital%kappa == kapc)  ic = i
         if (wave%rwf(i)%orbital%n     == nd   .and.  &
	     wave%rwf(i)%orbital%kappa == kapd)  id = i
      end do
      !
      if (ia == 0  .or.  ib == 0  .or.  ic == 0  .or.  id == 0) then
         print *, "Radial functions are not found for all the given orbitals:" 
	 print *, "   ",orbital_name(na,kapa), orbital_name(nb,kapb), &
                        orbital_name(nc,kapc), orbital_name(nd,kapd)
         print *, "Input error; redo ..."
         goto 1
      end if
      !
      wx = XL_Coulomb_strength_grasp2k(L,wave%rwf(ia),wave%rwf(ib),   &
                                         wave%rwf(ic),wave%rwf(id))
      !
      print *, " "
      print *, "X^" // trim(angular_momentum_string(L+L,4)) // "_Coulomb ("// &
               orbital_name(na,kapa) // "," // orbital_name(nb,kapb) //";"//  &
               orbital_name(nc,kapc) // "," // orbital_name(nd,kapd) //       &
	       ") = ",wx
      !       
      print *, "Continue (y/n) ?";   if (get_yes_stream()) goto 1
      return
      !
    5 print *, "Input error; redo ..."
      goto 1
      !     
   end subroutine toolbox_me_coulomb_strength
   !
   !
   subroutine toolbox_me_multipole()
   !--------------------------------------------------------------------
   ! Calculates the one-electron multipole matrix elements of the electron-
   ! photon interaction for a given set of atomic orbitals.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer          :: sel, ia, ib, na, nb, kapa, kapb, L, p
      real(kind=dp)    :: energy, wx1, wx2, wx3, wx4
      character(len=9) :: gauge
      
      !
      ! Define radial grid parameter and generate the grid
      rnt_grasp2k = 2.0e-6_dp
      h_grasp2k   = 5.0e-2_dp
      n_grasp2k   = 390
      hp_grasp2k  = zero
      !
      ! Default speed of light
      c = c_vacuum
      !
      ! accy is an estimate of the accuracy of the numerical procedures
      accy_grasp2k = h_grasp2k**6
      !
      ! Set up the coefficients for the numerical procedures
      call setqic_grasp2k()
      !
      ! Generate the radial grid and all associated arrays
      call radgrd_grasp2k()
      !
      ! Get the occupied shells of the GRASP92 states
      call file_get_csl_list( &
         "Enter the name of the configuration symmetry list file:", &
	 asf_set%csf_set)
      !
      ! Initialize and load the radial wavefunctions 
      call initialize_rwf_storage(asf_set, wave)
      call file_get_rwf("Enter the name of the Radial WaveFunction File:", &
	 asf_set,wave,.true.)
      !
    1 print *, "Select the gauge form of the matrix elements:"
      print *, "  1 - Coulomb gauge."
      print *, "  2 - Babushkin gauge."
      read  *, sel
      select case(sel)
      case(1)
         gauge = "Coulomb  "
      case(2)
         gauge = "Babushkin"
      case default
         print *, "Unable to decode selection; redo ... "
         goto 1
      end select
      !
      print *, "Enter the (6!) quantum numbers na, kapa, nb, kapb, L:"
      print *, "  magnetic (p=0), electric (p=1) "
      read(*,*,err=5)  na, kapa, nb, kapb, L, p
      !
      ia = 0;   ib = 0
      do  i = 1,wave%number_of_rwf
         if (wave%rwf(i)%orbital%n     == na   .and.  &
	     wave%rwf(i)%orbital%kappa == kapa)  ia = i
         if (wave%rwf(i)%orbital%n     == nb   .and.  &
	     wave%rwf(i)%orbital%kappa == kapb)  ib = i
      end do
      !
      if (ia == 0  .or.  ib == 0) then
         print *, "Radial functions are not found for all the given orbitals:" 
	 print *, "   ",orbital_name(na,kapa), orbital_name(nb,kapb)
         print *, "Input error; redo ..."
         goto 1
      end if
      !
      energy = abs(wave%rwf(ia)%energy - wave%rwf(ib)%energy)
      wx1    = multipole_reduced_M_integral(energy,L,p,gauge,&
                                            wave%rwf(ia),wave%rwf(ib))
      wx2    = multipole_reduced_M_andreys(energy,L,p,gauge,&
                                            wave%rwf(ia),wave%rwf(ib))
      wx3    = multipole_reduced_M_capture(energy,L,p,gauge,&
                                            wave%rwf(ia),wave%rwf(ib))
      wx4    = multipole_reduced_M_multiphoton(energy,L,p,gauge,&
                                            wave%rwf(ia),wave%rwf(ib))
      !    
      print *, " "
      print *, "M^" // trim(angular_momentum_string(L+L,4)) // "_integral ("// &
               orbital_name(na,kapa) // "," // orbital_name(nb,kapb) //        &
	       ") = ",wx1
      print *, "M^" // trim(angular_momentum_string(L+L,4)) // "_capture  ("// &
               orbital_name(na,kapa) // "," // orbital_name(nb,kapb) //        &
	       ") = ",wx2
      print *, "M^" // trim(angular_momentum_string(L+L,4)) // "_andreys  ("// &
               orbital_name(na,kapa) // "," // orbital_name(nb,kapb) //        &
	       ") = ",wx3
      print *, "M^" // trim(angular_momentum_string(L+L,4)) // "_multiphtn("// &
               orbital_name(na,kapa) // "," // orbital_name(nb,kapb) //        &
	       ") = ",wx4
      !       
      print *, "Continue (y/n) ?";   if (get_yes_stream()) goto 1
      return
      !
    5 print *, "Input error; redo ..."
      goto 1
      !     
      !
   end subroutine toolbox_me_multipole
   !
   !
   function toolbox_merge_corestring(n,kappa)               result(cadd)
   !--------------------------------------------------------------------
   ! Returns the corestring for a given subshell.
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: n, kappa
      character(len=9)    :: cadd
      !
      if      (n == 1 .and. kappa == -1) then;   cadd = "  1s ( 2)"
      else if (n == 2 .and. kappa == -1) then;   cadd = "  2s ( 2)"
      else if (n == 2 .and. kappa ==  1) then;   cadd = "  2p-( 2)"
      else if (n == 2 .and. kappa == -2) then;   cadd = "  2p ( 4)"
      else if (n == 3 .and. kappa == -1) then;   cadd = "  3s ( 2)"
      else if (n == 3 .and. kappa ==  1) then;   cadd = "  3p-( 2)"
      else if (n == 3 .and. kappa == -2) then;   cadd = "  3p ( 4)"
      else if (n == 3 .and. kappa ==  2) then;   cadd = "  3d-( 4)"
      else if (n == 3 .and. kappa == -3) then;   cadd = "  3d ( 6)"
      else if (n == 4 .and. kappa == -1) then;   cadd = "  4s ( 2)"
      else if (n == 4 .and. kappa ==  1) then;   cadd = "  4p-( 2)"
      else if (n == 4 .and. kappa == -2) then;   cadd = "  4p ( 4)"
      else if (n == 4 .and. kappa ==  2) then;   cadd = "  4d-( 4)"
      else if (n == 4 .and. kappa == -3) then;   cadd = "  4d ( 6)"
      else if (n == 4 .and. kappa ==  3) then;   cadd = "  4f-( 6)"
      else if (n == 4 .and. kappa == -4) then;   cadd = "  4f ( 8)"
      else if (n == 5 .and. kappa == -1) then;   cadd = "  5s ( 2)"
      else if (n == 5 .and. kappa ==  1) then;   cadd = "  5p-( 2)"
      else if (n == 5 .and. kappa == -2) then;   cadd = "  5p ( 4)"
      else if (n == 5 .and. kappa ==  2) then;   cadd = "  5d-( 4)"
      else if (n == 5 .and. kappa == -3) then;   cadd = "  5d ( 6)"
      else if (n == 5 .and. kappa ==  3) then;   cadd = "  5f-( 6)"
      else if (n == 5 .and. kappa == -4) then;   cadd = "  5f ( 8)"
      else if (n == 6 .and. kappa == -1) then;   cadd = "  6s ( 2)"
      else if (n == 6 .and. kappa ==  1) then;   cadd = "  6p-( 2)"
      else if (n == 6 .and. kappa == -2) then;   cadd = "  6p ( 4)"
      else if (n == 6 .and. kappa ==  2) then;   cadd = "  6d-( 4)"
      else if (n == 6 .and. kappa == -3) then;   cadd = "  6d ( 6)"
      else if (n == 6 .and. kappa ==  3) then;   cadd = "  6f-( 6)"
      else if (n == 6 .and. kappa == -4) then;   cadd = "  6f ( 8)"
      else
         stop "toolbox_merge_corestring(): program stop A."
      end if
      !
   end function toolbox_merge_corestring
   !
   !
   subroutine toolbox_merge_csl_diffcores()
   !--------------------------------------------------------------------
   ! Controls the 'merging' of two .csl lists with different cores. 
   ! In contrast to mrgcsl from GRASP92 program, this procedure is much 
   ! faster and does not require that the same list of orbitals is defined.
   ! However, the sequence of orbitals must agree for those orbitals
   ! which are defined in both .csl lists.
   !
   ! Calls: file_get_csl_list(), merge_check_csf_set(), 
   !        toolbox_merge_set_reduced_list(),
   !        toolbox_merge_two_list(), toolbox_merge_write_list().
   !--------------------------------------------------------------------
      !
      integer            :: nocsf_a, nocsf_b, ii, iflag
      character(len=256) :: csl_old_file_a, csl_old_file_b
      character(len=80)  :: recore, cadd
      !
      ! Open, check, load data from, and close, the  first .csl  file
      call file_get_csl_list(                                              &
         "Enter the name of the first configuration symmetry list file:",  &
         csf_set_a,csl_old_file_a)
      !
      ! Open, check, load data from, and close, the  second .csl  file
      call file_get_csl_list(                                              &
         "Enter the name of the second configuration symmetry list file:", &
         csf_set_b,csl_old_file_b)
      !
      ! Check the consistency of the file and create the 'new' sequence of
      ! orbitals
      call toolbox_merge_check_csf_set(.false.)
      csf_set_c%nwcore  = min(csf_set_a%nwcore, csf_set_b%nwcore)
      !
      ! Construct the missing `core' string for those orbitals which are 
      ! not yet defined in the list with larger nwcore
      iflag = 0
      if (csf_set_a%nwcore > csf_set_b%nwcore) then
         iflag = 1
         do  i = csf_set_b%nwcore + 1,csf_set_a%nwcore
	    cadd(1:9) = toolbox_merge_corestring(csf_set_b%subshell(i)%n, &
	                                         csf_set_b%subshell(i)%kappa)
	    ii = 9*(i-csf_set_b%nwcore-1)			      
	    recore(ii+1:ii+9) = cadd(1:9)
	    print *, "i, ii, cadd, recore = ",i, ii, cadd(1:9),recore(1:ii+9)
	 end do
      else if (csf_set_b%nwcore > csf_set_a%nwcore) then
         iflag = 2
         do  i = csf_set_a%nwcore + 1,csf_set_b%nwcore
	    cadd(1:9) = toolbox_merge_corestring(csf_set_a%subshell(i)%n, &
	                                         csf_set_a%subshell(i)%kappa)
	    ii = 9*(i-csf_set_a%nwcore-1)			      
	    recore(ii+1:ii+9) = cadd(1:9)
	    print *, "i, ii, cadd, recore = ",i, ii, cadd(1:9),recore(1:ii+9)
	 end do
      end if 
      !
      recore(ii+10:80) = " "
      print *, "recore = ",recore(1:ii+9)
      !
      merge_noint = csf_set_c%nwshells / merge_nobit + 1
      !
      nocsf_a = csf_set_a%nocsf
      allocate( list_a(1:nocsf_a) )
      call toolbox_merge_set_reduced_list(csf_set_a,list_a)
      call deallocate_csf_basis(csf_set_a)
      !
      nocsf_b = csf_set_b%nocsf
      allocate( list_b(1:nocsf_b) )
      call toolbox_merge_set_reduced_list(csf_set_b,list_b)
      call deallocate_csf_basis(csf_set_b)
      !
      print *, " ... packing complete."
      !
      ! Compare the lists and eliminate all doubly defined CSF from list_b
      call toolbox_merge_two_list(list_a,nocsf_a,list_b,nocsf_b)
      !
      ! Write out the combined list
      merge_csf_counter = 0
      print *, "csf_set_a%nwcore, csf_set_b%nwcore = ", &
                csf_set_a%nwcore, csf_set_b%nwcore
      if (iflag == 1) then
         call toolbox_merge_write_list(csl_old_file_a,list_a,nocsf_a,.false., &
                                       .true.,recore)
         call toolbox_merge_write_list(csl_old_file_b,list_b,nocsf_b,.true.,  &
                                       .false.,recore)
      else if (iflag == 2) then
         call toolbox_merge_write_list(csl_old_file_a,list_a,nocsf_a,.false., &
                                       .false.,recore)
         call toolbox_merge_write_list(csl_old_file_b,list_b,nocsf_b,.true.,  &
                                       .true.,recore)
      else 
         call toolbox_merge_write_list(csl_old_file_a,list_a,nocsf_a,.false., &
                                       .false.,recore)
         call toolbox_merge_write_list(csl_old_file_b,list_b,nocsf_b,.true.,  &
                                       .false.,recore)
      end if			       
      close(22)
      !
      print *, "Write out of the new .csl file complete;"
      print *, " there are ",csf_set_c%nwshells," subshells"
      print *, " in ",merge_csf_counter," relativistic CSF in this list."
      !
   end subroutine toolbox_merge_csl_diffcores
   !
   !
   subroutine toolbox_merge_csl_lists()
   !--------------------------------------------------------------------
   ! Controls the 'merging' of two .csl lists. In contrast to a previous
   ! versions of the GRASP92 program, mrgcsl, this program is much faster 
   ! and also supports a larger variety of .csl files which need not to 
   ! have the same orbitals.
   !
   ! Calls: file_get_csl_list(), merge_check_csf_set(), 
   !        toolbox_merge_set_reduced_list(),
   !        toolbox_merge_two_list(), toolbox_merge_write_list().
   !--------------------------------------------------------------------
      !
      integer            :: nocsf_a, nocsf_b
      character(len=256) :: csl_old_file_a, csl_old_file_b
      !
      ! Open, check, load data from, and close, the  first .csl  file
      call file_get_csl_list(                                              &
         "Enter the name of the first configuration symmetry list file:",  &
         csf_set_a,csl_old_file_a)
      !
      ! Open, check, load data from, and close, the  second .csl  file
      call file_get_csl_list(                                              &
         "Enter the name of the second configuration symmetry list file:", &
         csf_set_b,csl_old_file_b)
      !
      ! Check the consistency of the file and create the 'new' sequence of
      ! orbitals
      call toolbox_merge_check_csf_set(.true.)
      !
      merge_noint = csf_set_c%nwshells / merge_nobit + 1
      !
      nocsf_a = csf_set_a%nocsf
      allocate( list_a(1:nocsf_a) )
      call toolbox_merge_set_reduced_list(csf_set_a,list_a)
      call deallocate_csf_basis(csf_set_a)
      !
      nocsf_b = csf_set_b%nocsf
      allocate( list_b(1:nocsf_b) )
      call toolbox_merge_set_reduced_list(csf_set_b,list_b)
      call deallocate_csf_basis(csf_set_b)
      !
      print *, " ... packing complete."
      !
      ! Compare the lists and eliminate all doubly defined CSF from list_b
      call toolbox_merge_two_list(list_a,nocsf_a,list_b,nocsf_b)
      !
      ! Write out the combined list
      merge_csf_counter = 0
      call toolbox_merge_write_list(csl_old_file_a,list_a,nocsf_a,.false.,  &
                                    .false.," ")
      call toolbox_merge_write_list(csl_old_file_b,list_b,nocsf_b,.true.,   &
                                    .false.," ")
      close(22)
      !
      print *, "Write out of the new .csl file complete;"
      print *, " there are ",csf_set_c%nwshells," subshells"
      print *, " in ",merge_csf_counter," relativistic CSF in this list."
      !
   end subroutine toolbox_merge_csl_lists
   !
   !
   subroutine toolbox_merge_check_csf_set(same_core)
   !--------------------------------------------------------------------
   ! Checks that the two CSF lists in csf_set_a and csf_set_b have the
   ! same number of electrons, the same 'core' as well as the same order
   ! of subshells for those which are defined in both lists.
   ! It prints an ERROR message and terminates the executation if this
   ! is not the case.
   ! This subroutine also merges the two orbital list and generates a
   ! unified orbital list for the derived structure csf_set_c.
   !--------------------------------------------------------------------
      !
      logical, intent(in) :: same_core
      !
      integer                                :: i, j, last, last_in_b, nw
      integer, dimension(csf_set_a%nwshells) :: in_b
      !
      if (csf_set_a%number_of_electrons /= csf_set_b%number_of_electrons) then
         print *, "Different number of electrons in the two .csl list;"
         print *, " "
         stop     "toolbox_merge_check_csf_set(): terminates for this reason ... ."
      else if (same_core  .and.  csf_set_a%nwcore /= csf_set_b%nwcore) then
         print *, "Different 'cores' are defined in the two .csl list;"
         print *, " "
         stop     "toolbox_merge_check_csf_set(): terminates for this reason ... ."
      end if
      !
      in_b(:) = 0
      !
      do  i = 1,csf_set_a%nwshells
         do  j = 1,csf_set_b%nwshells
            if (csf_set_b%subshell(j)%n == csf_set_a%subshell(i)%n  .and.   &
                csf_set_b%subshell(j)%kappa == csf_set_a%subshell(i)%kappa) &
                then
               in_b(i) = j
               exit
            end if
         end do
      end do
      !
      ! Now determine that the numbers in_b(:) are increasing (if not zero)
      last = 0
      do  i = 1,csf_set_a%nwshells
         if (in_b(i) /= 0) then
            if (in_b(i) <= last) then
               print *, "Those (valence) subshells which are defined" //&
                        " in both lists must appear in the same order;"
               print *, " "
               stop     "toolbox_merge_check_csf_set(): terminates for " //&
	                "this reason ... ."
            else
                  last = in_b(i)
            end if
         end if
      end do
      !!x print *, "in_b(:) = ",in_b(:)
      !
      ! Merge the two orbital lists in a correct way
      nw = 0
      do  i = 1,csf_set_a%nwshells
         if (in_b(i) /= 0) nw = nw + 1
      end do
      nw = nw + (csf_set_a%nwshells - nw) + (csf_set_b%nwshells - nw)
      allocate( csf_set_c%subshell(1:nw) )
      csf_set_c%number_of_electrons = csf_set_a%number_of_electrons
      csf_set_c%nwcore              = csf_set_a%nwcore
      csf_set_c%nwshells            = nw
      !
      last_in_b = 0;   nw = 0
      do  i = 1,csf_set_a%nwshells
         if (in_b(i) /= 0) then
          1 if (last_in_b + 1 < in_b(i)) then
               last_in_b = last_in_b + 1
               nw = nw + 1
               csf_set_c%subshell(nw) = nkappa(csf_set_b%subshell(i)%n,   &
                                               csf_set_b%subshell(i)%kappa)
               goto 1
            end if   
            last_in_b = in_b(i)
            nw = nw + 1
            csf_set_c%subshell(nw) = nkappa(csf_set_a%subshell(i)%n,   &
                                            csf_set_a%subshell(i)%kappa)
         else
            nw = nw + 1
            csf_set_c%subshell(nw) = nkappa(csf_set_a%subshell(i)%n,   &
                                            csf_set_a%subshell(i)%kappa)
         end if
      end do
      !
    2 if (last_in_b < csf_set_b%nwshells) then
         last_in_b = last_in_b + 1
         nw = nw + 1
         csf_set_c%subshell(nw) = nkappa(csf_set_b%subshell(last_in_b)%n,   &
                                         csf_set_b%subshell(last_in_b)%kappa)
         goto 2
      end if   
      !
   end subroutine toolbox_merge_check_csf_set
   !
   !
   subroutine toolbox_merge_orbitals_formatted()
   !--------------------------------------------------------------------
   ! Re-collects and combines radial orbitals from two or more (formatted)
   ! .out files into a new one. It assumes that the radial grids are the
   ! same, so all information is copied as strings.
   !--------------------------------------------------------------------
      !
      logical                 :: yes, fail
      logical, dimension(200) :: orbitals_still_needed
      character(len=26)       :: g92rwf
      character(len=256)      :: file_name, record
      character(len=500)      :: pstr
      integer                 :: i, ierr, ios, j, k, No_osn, nc, no_orbitals, &
                                 pqn, kappa, pqnx, kappax, mtpx
      real(kind=dp)           :: wax
      type(csf_basis)         :: csf_list
      type(nkappa), dimension(500) :: orbital
      !
      ! Read in a .csl file to determine the list of orbitals which need
      ! to be collected from the different .out files
      call file_get_csl_list(                                       &
         "Enter the name of the configuration symmetry list file:", &
         csf_list)
      !
      orbitals_still_needed(:) = .true.
      !
      ! Open a new formatted file for the orbtials to be combined in
    2 print *, "Enter a file name for the 'formatted' .out Radial "// &
               "Wavefunction File that is to be generated:"
      read (*,"(a)") file_name
      call file_open(26,file_name,"formatted  ","new",ierr)
      if (ierr /= 0  .or.  len_trim(file_name) == 0) goto 2
      !
      ! Write out a file header
      write (26,"(a)") "G92RWF (formatted file version).  "
      !
      !
    3 continue
      pstr  = trim(" ")
      !
      !!x print *, "csf_list%nwshells = ",csf_list%nwshells
      do  i = 1,csf_list%nwshells
         if (orbitals_still_needed(i)) then
            !!x print *, "n, kappa, name = ", csf_list%subshell(i)%n, &
            !!x          csf_list%subshell(i)%kappa,                  &
            !!x          orbital_name(csf_list%subshell(i)%n,         &
            !!x                       csf_list%subshell(i)%kappa)
            !
            ios = len(trim(pstr))
            if (pstr(ios:ios) == "-")  then
               pstr = trim(pstr)//" "//orbital_name(csf_list%subshell(i)%n,   &
                                                    csf_list%subshell(i)%kappa)
            else
               pstr = trim(pstr)//"  "//orbital_name(csf_list%subshell(i)%n,  &
                                                     csf_list%subshell(i)%kappa)
            end if
         end if  
      end do
      !
      print 21, pstr(1:101)
      if (len(trim(pstr)) > 101)  print 22, pstr(102:202)
      if (len(trim(pstr)) > 202)  print 22, pstr(202:303)
      !
      ! Attempt to open the formatted file and read in all data
    4 print *, "Enter the name of another (formatted) GRASP92 Radial "//&
               "Wavefunction File:"
      read (*,"(a)") file_name
      call file_open(27,file_name,"formatted  ","old",ierr)
      if (ierr /= 0  .or.  len_trim(file_name) == 0) goto 4
      !
      ! Check the header of the file; if not as expected, try again
      read (27,"(a)") g92rwf
      if (g92rwf(1:17) /= "G92RWF (formatted") then
         close (27)
         print *, "g92rwf = ",g92rwf
         print *, "Not a formatted GRASP92 radial wavefunction file;"
         goto 4
      end if
      !
    5 print *, "Enter the list of orbitals to be `used' from this file:"
      read (*,"(a)") record
      print *, "main: record = ",record
      call interprete_orbitals(record,orbital,no_orbitals,500,fail)
      if (fail) goto 5
      !
      do  j = 1,no_orbitals
         pqn   = orbital(j)%n
         kappa = orbital(j)%kappa
         !
         do  i = 1,csf_list%nwshells
            if (orbitals_still_needed(i)            .and.  &
                pqn   == csf_list%subshell(i)%n     .and.  &  
                kappa == csf_list%subshell(i)%kappa) then  
               print *, "Now search for this orbital; n,kappa = ",pqn,kappa
               !
               rewind(27)
               read(27,"(a)") record
               !
             6 read (27,*,end=7) pqnx, kappax, wax, mtpx
               !!x print *, "pqnx, kappax, wax, mtpx = ",pqnx, kappax, wax, mtpx
               if (pqn == pqnx  .and.  kappa == kappax) then
                  backspace (27)
                  do  k = 1,mtpx+2
                     read (27,"(a)") record
                     nc = len(trim(record))
                     write(26,"(a)") record(1:nc)
                  end do
                  orbitals_still_needed(i) = .false.
                  goto 8
               else
                  do  k = 1,mtpx+1
                     read (27,"(a)") record
                  end do
                  !
                  goto 6
               end if
            end if
         end do
         goto 8
         !
       7 print *,"Orbital not found; n,kappa = ",pqn,kappa 
       8 continue
      end do
      !
      close (27)
      !
      do i = 1,csf_list%nwshells
         if (orbitals_still_needed(i)) goto 3
      end do
      !
   21 format(/," Orbitals still be needed for the given .csl list:", &
             /," -------------------------------------------------", &
             /,"  ", a101)
   22 format(  "  ", a101)         
      !
   end subroutine toolbox_merge_orbitals_formatted
   !
   !
   subroutine toolbox_merge_set_reduced_list(csf_set,list)
   !--------------------------------------------------------------------
   ! Packs the information about the CSFs into a much more compact format
   ! to accelerate the comparison of different CSF in the later run.
   !--------------------------------------------------------------------
      !
      type(csf_basis), intent(in)     :: csf_set 
      type(reduced_csf), dimension(:) :: list
      !
      integer :: i, int, j, jj, pos
      integer, dimension(1:csf_set_c%nwshells) :: ndx
      !
      ! First map the orbital indices in csf_set onto those in csf_set_c
      ndx(:) = 0
      i_loop: do  i = 1,csf_set_c%nwshells
         do  j = 1,csf_set%nwshells
            if (csf_set_c%subshell(i)%n == csf_set%subshell(j)%n  .and.   &
                csf_set_c%subshell(i)%kappa == csf_set%subshell(j)%kappa) then
               ndx(i) = j
               cycle i_loop
            end if
         end do
      end do  i_loop
      !
      ! Now pack the CSF list into a reduced format
      do  i = 1,csf_set%nocsf
         !
         list(i)%append = .true.
         allocate(list(i)%red_shell(1:csf_set_c%nwshells),                    &
                  list(i)%red_occ(1:merge_noint), list(i)%red_X(1:merge_noint))
         !
         list(i)%red_occ(:) = 0;   list(i)%red_X(:) = 0
         int = 1
         do  j = 1,csf_set_c%nwshells
            pos = mod(j-1,merge_nobit)
            if (ndx(j) /= 0) then
               jj = ndx(j)
               if (csf_set%csf(i)%occupation(jj) /= 0) then
                  list(i)%red_occ(int) = ibset(list(i)%red_occ(int),pos)
               end if
               if (mod(csf_set%csf(i)%subshellX(jj)*1,4) == 2   .or.  &
                   mod(csf_set%csf(i)%subshellX(jj)*1,4) == 3) then
                  list(i)%red_X(int) = ibset(list(i)%red_X(int),pos)
               end if
            else
               if (mod(csf_set%csf(i)%subshellX(csf_set%nwshells)*1,4)==2 .or.&
                   mod(csf_set%csf(i)%subshellX(csf_set%nwshells)*1,4)==3) then
                  list(i)%red_X(int) = ibset(list(i)%red_X(int),pos)
               end if
            end if
            if (mod(j,merge_nobit) == 0) int = int + 1
         end do
         !
         do  j = 1,csf_set_c%nwshells
            if (ndx(j) /= 0) then
               jj = ndx(j)
               list(i)%red_shell(j) = ((csf_set%csf(i)%occupation(jj) *  10 + &
                                        csf_set%csf(i)%seniority(jj)) * 100 + &
                                        csf_set%csf(i)%subshellJ(jj)) * 100 + &
                                        csf_set%csf(i)%subshellX(jj)
            else
               list(i)%red_shell(j) = csf_set%csf(i)%subshellX(csf_set%nwshells)
            end if
         end do
      end do
      !
   end subroutine toolbox_merge_set_reduced_list
   !
   !
   subroutine toolbox_merge_two_list(list_a,n_a,list_b,n_b)
   !--------------------------------------------------------------------
   ! Compares all (reduced) CSF from list_b with those from list_a and
   ! 'eliminates' the doubly defined CSF in list_b (i.e. set
   ! list_b(i)%append = .false.
   !--------------------------------------------------------------------
      !
      integer, intent(in)                            :: n_a, n_b
      type(reduced_csf), dimension(:), intent(in)    :: list_a
      type(reduced_csf), dimension(:), intent(inout) :: list_b
      !
      integer :: i, int, j, k, count
      !
      count = n_a + n_b
      i_loop:  do  i = 1,n_b
         if (mod(i,1000) == 0) then
            print *, " ... merge complete for ",i,"CSF."
         end if
         j_loop:  do  j = 1,n_a
            do  int = 1,merge_noint
               if (list_a(j)%red_occ(int) /= list_b(i)%red_occ(int)) then
                  cycle  j_loop
               else if (list_a(j)%red_X(int) /= list_b(i)%red_X(int)) then
                  cycle  j_loop
               end if
            end do
            !
            ! Now check the two CSF in detail
            do  k = 1,csf_set_c%nwshells
               if (list_a(j)%red_shell(k) /= list_b(i)%red_shell(k)) then
                  ! print *, "Merge - c"
                  cycle  j_loop
               end if
            end do
            !
            ! This point is reached only if both (reduced) CSF are the same
            list_b(i)%append = .false.
            count            = count - 1
            cycle  i_loop
         end do  j_loop
      end do i_loop
      !
      print *, " ... there remain",count,"CSF in the merged list."
      !
   end subroutine toolbox_merge_two_list
   !
   !
   subroutine toolbox_merge_write_list(csl_old_file,list,nocsf,append_to_file, &
                                       append_corestring,corestring)
   !--------------------------------------------------------------------
   ! Writes the 'accepted' CSF from csl_old_file to a new .csl file.
   ! For append_to_file == .false., it opens a new file and starts with 
   ! the header of .csl file; otherwise, it append to a 'previously opened'
   ! file and starts with the first 'accepted' CSF.
   !--------------------------------------------------------------------
      !
      integer, intent(in)             :: nocsf
      logical, intent(in)             :: append_to_file, append_corestring
      character(len= *),  intent(in)  :: corestring
      character(len=256), intent(in)  :: csl_old_file
      type(reduced_csf), dimension(:) :: list
      !
      integer            :: i, ierr, ilast, na
      character(len=512) :: csl_new_file, record1, record2, record3, &
                            core, peel, ca 
      !
      ! Reopen the 'old' .csl file; no problems should arise for the given
      ! name since this file has been opened successful previously
      call file_open(21,csl_old_file,"formatted  ","old",ierr)
      if (rabs_use_stop   .and.  ierr == 1) then
         stop "merge_write_csl(): program stop A."
      end if
      !
      !
      if (.not.append_to_file) then
       1 print *, "Enter the name of the configuration symmetry list file" //&
                  " that is to be created:"
         read (*,"(a)") csl_new_file
         if (len(trim(csl_new_file)) == 0) goto 1
         !
         call file_open(22,csl_new_file,"formatted  ","new",ierr)
         if (ierr == 1) goto 1
         !
         ! Create the head of the .csl file
         read (21,"(a)") record1
         read (21,"(a)") 
         read (21,"(a)") record2
         read (21,"(a)") 
         read (21,"(a)") record3
         !
         core(:) = " ";   peel(:) = " "
         ilast = 0
         do  i = 1,csf_set_c%nwcore
            core(ilast+1:ilast+1) = " "
            core(ilast+2:ilast+5) = orbital_name(csf_set_c%subshell(i)%n,   &
                                                 csf_set_c%subshell(i)%kappa)
            ilast = ilast + 5
         end do
         ilast = 0
         do  i = csf_set_c%nwcore+1,csf_set_c%nwshells
            peel(ilast+1:ilast+1) = " "
            peel(ilast+2:ilast+5) = orbital_name(csf_set_c%subshell(i)%n,   &
                                                 csf_set_c%subshell(i)%kappa)
            ilast = ilast + 5
         end do
         !
         write(22,"(a)") trim(record1)
         write(22,"(a)") trim(core)
         write(22,"(a)") trim(record2)
         write(22,"(a)") trim(peel)
         write(22,"(a)") trim(record3)
      else
         read (21,"(a)") record1;  read (21,"(a)") record1
         read (21,"(a)") record1;  read (21,"(a)") record1
         read (21,"(a)") record1 
      end if
      !
      ! Now read 'one by one' CSF and write it out if it belongs to the
      ! accepted list
      do  i = 1,nocsf
         read(21,"(a)") record1
         read(21,"(a)") record2
         read(21,"(a)") record3
         if (list(i)%append) then
            merge_csf_counter = merge_csf_counter + 1
	    !
	    if (append_corestring) then
	       print *, "trim(corestring) = ",trim(corestring)
	       ca = "                                                      "
	       record1 = trim(corestring) // record1
	       na = len( trim(corestring) )
	       !!x print *, "na, corestring = ",na, corestring
	       record2 = ca(1:na) // record2
	       record3 = ca(1:na) // record3
	    end if
	    !
            write(22,"(a)") trim(record1)
            write(22,"(a)") trim(record2)
            write(22,"(a)") trim(record3)
         end if
      end do
      !
      close(21)
      !
   end subroutine toolbox_merge_write_list
   !
   !
   subroutine toolbox_modify_energies()
   !--------------------------------------------------------------------
   ! Modifies and overwrites the total energies of a given .mix file.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      logical                  :: yes
      character(len=30)        :: g92mix
      character(len=256)       :: file_name
      integer                  :: i, ierr, ios, j, number_of_levels
      integer, dimension(1000) :: levels
      real(kind=dp)            :: shift, shift_au, refen, refen_au
      !
      ! Determine the units for the printout and further optional input
      call input_energy_unit()
      !
      ! Default speed of light
      c = c_vacuum
      !
      !
      call file_get_csl_list("Enter a file name of the .csl file:", &
                             asf_set%csf_set)
      !
      call file_get_mix("Enter a file name for the (formatted) GRASP92 .mix "//&
                        "coefficient file:",asf_set)
      !
      !
      print *, "Enter energy shifts for all levels of the given .mix file ?"
      yes = get_yes_stream()
      if (yes) then
         ! Here, modify the energies
         do  i = 1,asf_set%noasf
            write(*,1) "Enter the shift (in ",trim(energy_unit),") for level ",&
	               asf_set%asf(i)%level_No,":"
            read(*,*)  shift
	    !
            if (energy_inverse) then
               shift_au = energy_factor / shift
            else
               shift_au = shift / energy_factor 
            end if
	    print *, "i, shift_au = ",i, shift_au
	    asf_set%asf(i)%energy = asf_set%asf(i)%energy + shift_au
	    !
         end do
       1 format(1x,a,a,a,i6,a)
         goto 3
      end if
      !
      !
      print *, "Enter total energies for all levels of the given .mix file ?"
      yes = get_yes_stream()
      if (yes) then
         ! 
         do  i = 1,asf_set%noasf
            write(*,1) "Enter the total energy (in ",trim(energy_unit), &
	               ") for level ",asf_set%asf(i)%level_No,";"
	    write(*,1) "give zero if nothing should be changed."
            read(*,*)  shift
	    !
	    if (abs(shift) < eps10) cycle
	    !
            if (energy_inverse) then
               shift_au = energy_factor / shift
            else
               shift_au = shift / energy_factor 
            end if
	    asf_set%asf(i)%energy = shift_au
	    !
         end do
         goto 3
      end if
      !
      !
      print *, "Enter excitation energies relative to a given (reference) "//&
               "energy for all levels of the given .mix file ?"
      yes = get_yes_stream()
      if (yes) then
         ! 
         write(*,1) "Enter the (total) reference energy (in ",                &
                    trim(energy_unit),") with regard to which all energies "//&
                    "are to be defined:"
       4 read(*,*)  refen
	 !
	 if (abs(refen) < eps10) then
            print *, "The reference energy is a total energy, abs(E)>> 0; "// &
                     "redo ..."
            goto 4
         end if
	 !
         if (energy_inverse) then
            refen_au = energy_factor / refen
         else
            refen_au = refen / energy_factor 
         end if
	 !
         do  i = 1,asf_set%noasf
            write(*,1) "Enter the excitation energy (in ",trim(energy_unit), &
	               ") for level ",asf_set%asf(i)%level_No,";"
	    write(*,1) "give zero if nothing should be changed."
            read(*,*)  shift
	    !
	    if (abs(shift) < eps10) cycle
	    !
            if (energy_inverse) then
               shift_au = energy_factor / shift
            else
               shift_au = shift / energy_factor 
            end if
	    asf_set%asf(i)%energy = refen_au + shift_au
	    !
         end do
         goto 3
      end if
      !
      !
      ! Open a 'formatted' file and write out all data
      !
    3 continue
      !
      number_of_levels = asf_set%noasf
      do  i = 1,number_of_levels
         levels(i) = asf_set%asf(i)%level_No
      end do
      !
      call  file_write_mix(asf_set,levels,number_of_levels)
      !
      !
   end subroutine toolbox_modify_energies
   !
   !
   subroutine toolbox_nuclear_fermi()
   !--------------------------------------------------------------------
   ! Calculate the Fermi distribution parameters c and N as well as the 
   ! Fermi distribution and the associated potential (if requested). These 
   ! computations follow the formulas by Tupitsyn et al., PRA 68 (2003) 
   ! 022511, Eqs. (8-11). For a given Fermi parameter a, in particular, the 
   ! parameter c and normalization N are calculated from analytical formulas.
   ! Using these parameters, the distribution: 
   !                 rho(r; R) = N / ( 1 + exp[(r-c)/a] )  
   ! as well as the associated nuclear potential V_n (r; R) can be printed 
   ! out on request.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer       :: i, k
      logical       :: yes
      real(kind=dp) :: R, a, a_au, c_fermi, c_au, N
      real(kind=dp), dimension(:), allocatable :: rho, Vn, ta
      !
      print *, "***** WARNING *****"
      print *, "***** WARNING *****"
      print *, "***** WARNING *****"
      print *, " "
      print *, "There are difficulties with the given formulas in this reference."
      print *, " "
      print *, " "
      !
      a = 2.3_dp / (two*two * log(three))
      print *, " log(ten) = ",log(ten)
      !
      print *, "Enter the nuclear radius R = sqrt( <r^2> ) [in fm]:"
      read  *, R
      !
      print *, "The Fermi nuclear skin is set to  a [in fm] = ",a
      print *, "  revise this setting ? "
      yes = get_yes_stream()
      if (yes) then
         print *, "Enter the Fermi nuclear skin a [in fm]:"
	 read  *, a
      end if
      !
      c_fermi = 5.0_dp / three * R*R  -  7.0_dp/three *pi*pi * a*a
      c_fermi = sqrt( c_fermi )
      !
      a_au = a * convert_fermi_to_bohr
      c_au = c_fermi * convert_fermi_to_bohr
      !
      N = three / (two*two * pi *c_au*c_au*c_au) &
                / (one + pi*pi*a_au*a_au / (c_au*c_au) )
      !
      print *, " "
      print *, "Fermi a [in fm]        = ",a
      print *, "Fermi c [in fm]        = ",c_fermi
      print *, "Normalization constant = ",N
      !
      print *, " "
      print *, "Calculate the associated nuclear contribution and potential ?"
      yes = get_yes_stream()
      if (yes) then
         !
         ! Determine the parameters controlling the radial grid
         call input_grid_parameters("standard")
         !
         ! Modify the radial grid parameters
         call input_grid_parameters("modify")
         !
         ! Set up the coefficients for the numerical procedures
         call setqic_grasp2k()
         !
         ! Generate the radial grid and all associated arrays
         call radgrd_grasp2k()
         !
	 allocate( rho(1:n_grasp2k) , Vn(1:n_grasp2k), ta(1:n_grasp2k+10) )
	 !
	 do  i = 1,n_grasp2k
	    rho(i) = N / (one + exp( (r_grasp2k(i) - c_au)/a_au ) )
	 end do
	 !
	 c = c_vacuum
	 !
	 do  i = 2,n_grasp2k
	    ta(:) = zero
	    do  k = 1,n_grasp2k
	       ta(k) = r_grasp2k(k)*r_grasp2k(k) * rho(k) &
	               / max( r_grasp2k(k), r_grasp2k(i) )
	    end do
            ta(2:n_grasp2k) = ta(2:n_grasp2k) * rp_grasp2k(2:n_grasp2k)
	    !
	    Vn(i)  = -two*two*pi/c * quad_grasp2k(ta,n_grasp2k)
	 end do
	 !
	 write(*,*) "      i      r(i)     rho(i)      Vn(i) "
	 write(*,*) "----------------------------------------"
	 do  i = 1,n_grasp2k
	    write(*,5) i, r_grasp2k(i),  rho(i), Vn(i) 
	 end do
       5 format(i6,1pe12.5,2x,1pe12.5,2x,1pe12.5)
	 !
      end if
      !
   end subroutine toolbox_nuclear_fermi
   !
   !
   subroutine toolbox_nuclear_radius()
   !--------------------------------------------------------------------
   ! Calculate the nuclear radius from the mass number by using the simple 
   ! formula: R = sqrt( <r^2> ) = 0.836 * A^1/3  +  0.570.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      real(kind=dp) :: A, R
      !
    1 print *, "Enter the nuclear mass A:"
      read *, A
      !
      R = 0.836_dp * (A**(one/three))  +  0.570_dp
      !
      print *, "R [in fm] = sqrt( <r^2> ) = ", R
      print *, "<r^2> [in fm^2]           = ",R*R
      !
   end subroutine toolbox_nuclear_radius
   !
end module rabs_toolbox_ln
