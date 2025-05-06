module rabs_cowf
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains all procedures which are specific to the COWF
! program for generating continuum orbitals for the GRASP-92 package.
! This includes the calculation of the potentials, the integration of
! the continuum orbitals as well as several procedures for file handling 
! and the storage of global data.
!-----------------------------------------------------------------------
   !
   use rabs_anco
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_grasp2k
   use rabs_input_dialog
   implicit none
   !
   private :: cowf_append_csp
                 ! Appends a continuum spinor to the .csp output file.
   private :: cowf_bronstein_260
                 ! Calculates an estimate of the definite integral No. 260
		 ! from Bronstein (1983).
   public  :: cowf_collect_input
                 ! Collects and proceeds all input for the generation of
		 ! continuum orbitals for the GRASP-92 package.
   private :: cowf_fk_r
                 ! Calculates the f^k_r(ab) angular coefficient.
   public  :: cowf_generate_csp
                 ! Generates all continuum spinors as selected by the input.
   private :: cowf_gk_r
                 ! Calculates the g^k_r(ab) angular coefficient.
   public  :: cowf_integrate_csp
                 ! Integrates the radial Dirac-Fock equations for a free 
		 ! electron continuum spinor with energy > 0.
   public  :: cowf_iterate_csp
                 ! Iterates the radial Dirac-Fock equations for a free electron
		 ! continuum spinor with energy > 0.
   private :: cowf_normalize
                 ! Collects some previous code that is not used currently.
   public  :: cowf_normalize_nonrel
                 ! Normalizes a free electron continuum spinor with energy > 0
		 ! and calculates the phase shift by applying a nonrelativistic
		 ! scheme for the large component.
   private :: cowf_normalize_wkb
                 ! Normalizes a free electron continuum spinor with energy > 0
		 ! and calculates the phase shift by applying the WKB method.
   public  :: cowf_open_csp
                 ! Opens a .csp Continuum Spinor file for the COWF program on
                 ! stream 25.
   private :: cowf_schmidt_orthogonalization
                 ! Carries out a Schmidt orthogonalization for a given
                 ! continuum spinor with respect to all bound state orbitals.
   private :: cowf_set_channels
                 ! Set all selected continuum channels to an appropriate
		 ! type(cowf_channel) data structure.
   public  :: cowf_set_debug
                 ! Open if appropriate a .dbg file and sets flags for debug
                 ! print out of the COWF program.
   private :: cowf_set_direct_potential
                 ! Calculates the 'direct' potential Y for a given continuum
		 ! spinor.
   public  :: cowf_set_drs_coefficients
                 ! Calculates the d_rs (alpha) weight coefficients for the
		 ! potentials.
   private :: cowf_set_exchange_potential
                 ! Calculates the 'exchange' potentials X_p and X_q for a 
		 ! given continuum spinor.
   private :: cowf_set_lagrange_multipliers
                 ! Calculates the 'off-diagonal' Lagrange multipliers for a 
		 ! given continuum spinor with respect to the bound-state 
		 ! functions.
   public  :: cowf_set_yk_coefficients
                 ! Calculates the y_k (ab) weight coefficients of the
		 ! direct potential terms.
   public  :: cowf_set_xk_coefficients
                 ! Calculates the x_k (abcd) weight coefficients of the 
		 ! exchange potential terms.
   private :: cowf_start_integration
                 ! Sets up P(1:6) and Q(1:6) to start the integration.
   private :: cowf_total_phase_shift
                 ! Calculates the total phase shift from a pure sinodial
                 ! behavior of the large component P as appropriate for
                 ! the neutral case
   private :: cowf_Vk_rs
                 ! Calculates the V^k_rs(abcd) angular coefficient.
   !
   ! Storage for the initial and final atomic states and wave functions
   type(asf_basis), public       :: asf_bound, asf_cont
   type(grasp2k_orbital), public :: wave_bound
   !
   ! Define an internal structure type(cowf_channel) which stores the
   ! information about the continuum channels to be calculated
   type :: cowf_channel
      integer          :: asf, level_No, kappa_c, totalJ, number_of_energies
      character(len=1) :: parity
      real(kind=dp), dimension(:), pointer :: energy
   end type cowf_channel
   !
   integer :: number_of_channels
   type(cowf_channel), dimension(:), allocatable :: channel
   !
   real(kind=dp), dimension(:), allocatable   :: cowf_gen_occupation
   real(kind=dp), dimension(:,:), allocatable :: cowf_drs
   !
   type :: cowf_yk_coefficient
      integer          :: k, a, b, d
      real(kind=dp)    :: y
   end type cowf_yk_coefficient
   !
   type :: cowf_xk_coefficient
      integer          :: k, a, b, c, d
      real(kind=dp)    :: x
   end type cowf_xk_coefficient
   !
   type :: cowf_eab_lagrange
      integer          :: a, b
      real(kind=dp)    :: epsilon
   end type cowf_eab_lagrange
   !
   integer :: number_of_yk, number_of_xk, number_of_epsilon_ab
   type(cowf_yk_coefficient), dimension(:), allocatable :: cowf_yk
   type(cowf_xk_coefficient), dimension(:), allocatable :: cowf_xk
   type(cowf_eab_lagrange), dimension(200)              :: cowf_eab
   !
   ! Define storage and parameters to control the iteration of the
   ! continuum spinors
   integer, parameter     :: cowf_maximal_iteration         = 48, &
                             cowf_min_points_per_wavelength = 12
   integer, public        :: cowf_extent_mtp                = 0
   !
   type(orbital_function) :: cowf_csp
   real(kind=dp), dimension(:), allocatable :: cowf_nuc_pot, cowf_y_pot, &
                                               cowf_xp_pot, cowf_xq_pot
   !
   ! Define global logical flags for the control of the COWF program; the
   ! default values for these flags may be overwritten interactively during 
   ! input time
   logical, public ::  cowf_average_normalization     = .false., &
                       cowf_norm_nonrel               = .false., &
                       cowf_norm_wkb                  = .false., &
                       cowf_norm_wkb_old              = .false., &
                       cowf_phaseshift_wkb            = .false., &
                       cowf_phaseshift_zero_potential = .true.,  &
                       cowf_phaseshift_coulomb        = .false., &
                       cowf_print_always_csf_scheme   = .false., &
                       cowf_schmidt_always            = .false., &
                       cowf_schmidt_final             = .true.,  &
                       cowf_solve_homogeneous_eqn     = .false., &
                       cowf_specify_r_max             = .false., &
                       cowf_start_homogeneous         = .false., &
                       cowf_use_formatted_mix_file    = .true.,  &
                       cowf_use_formatted_rwf_file    = .true.,  &
		       cowf_use_lagrange_multipliers  = .false.
   !
   real(kind=dp)    :: cowf_r_max = zero
   !
   ! Define logical flags for debugging individual procedures
   logical, public ::  debug_drs_coefficients         = .false., &
                       debug_chi_pq                   = .false., &
                       debug_integrate_csp            = .false., &
                       debug_iterate_csp              = .false., &
                       debug_lagrange_epsilon_ab      = .false., &
		       debug_normalize_csp            = .false., &
		       debug_start_integration        = .false., &
                       debug_yx_potentials            = .false., &
		       debug_yk_coefficients          = .false., &
                       debug_xk_coefficients          = .false. 
   !
   integer :: number_of_selected_energies, number_of_selected_channels
   real(kind=dp), dimension(100)      :: energy_selection
   type(cowf_channel), dimension(100) :: channel_selection
   type(orbital_function)             :: csp_new, csp_sav
   !
contains
   !
   subroutine cowf_aatemplate()
   !--------------------------------------------------------------------
   ! Calculates 
   ! 
   ! Calls: spherical_Bessel_jL().
   !--------------------------------------------------------------------
      !
      !
   end subroutine cowf_aatemplate
   !
   !
   subroutine cowf_append_csp()
   !--------------------------------------------------------------------
   ! Appends the generated continuum spinor in cowf_csp to the .csp
   ! output file on stream 25. 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer :: i
      !
      write(25,*) cowf_csp%orbital%n, cowf_csp%orbital%kappa,    &
                  cowf_csp%energy, cowf_csp%mtp, cowf_csp%phase, &
		  " = n, kappa, energy, mtp, phase"
      write(25,*) cowf_csp%pz, " = pz"
      do  i = 1,cowf_csp%mtp
         write(25,*) r_grasp2k(i),cowf_csp%P(i),cowf_csp%Q(i)
      end do
      !
   end subroutine cowf_append_csp
   !
   !
   function cowf_bronstein_260(a,b,c,xa,xb)            result(value_260)
   !--------------------------------------------------------------------
   ! Returns the value of the integral No. 260 as defined by Bronstein 
   ! and Semendjajew : Taschenbuch der Mathematik (Leipzig 1983) over 
   ! the interval XA ... XB.  
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      real(kind=dp), intent(in) :: a, b, c, xa, xb
      real(kind=dp) :: value_260, wa, wb
      !
      wa = sqrt(a * xa * xa  +  b * xa  +  c)
      wb = sqrt(a * xb * xb  +  b * xb  +  c)
      value_260 = wb  -  wa  -  (b / two) *  cowf_bronstein_241(a,b,c,xa,xb)
      value_260 = value_260  +  c * cowf_bronstein_258(a,b,c,xa,xb)
      !
      contains
      !
      function cowf_bronstein_241(a,b,c,xa,xb)         result(value_241)
      !-----------------------------------------------------------------
      ! Returns the value of the integral No. 241 as defined by Bronstein 
      ! and Semendjajew : Taschenbuch der Mathematik (Leipzig 1983) over 
      ! the interval XA ... XB.
      !-----------------------------------------------------------------
         !
         real(kind=dp), intent(in) :: a, b, c, xa, xb
         real(kind=dp) :: value_241
	 !
	 real(kind=dp) :: wa, wb, wca, wcb, wd
         !
         if (a > zero) then
            wa  = a * xa * xa  +  b * xa  +  c
            wb  = a * xb * xb  +  b * xb  +  c
            wca = dsqrt(a * wa)
            wca = two * wca  +  two * a * xa  +  b
            wca = dlog( wca ) / dsqrt( a )
            wcb = dsqrt(a * wb)
            wcb = two * wcb  +  two * a * xb  +  b
            wcb = log(wcb) / sqrt(a)
            value_241 = wcb - wca
         else if (a < zero) then
            wd  = two * two * a * c  -  b * b
            if (rabs_use_stop   .and.   wd >= zero) then
               stop "cowf_bronstein_241(): program stop A."
            end if
            wd  = sqrt(-wd)
            !
            wca = two * a * xa  +  b
            wca = wca / wd
            wca = asin(wca)
            wca = -wca / dsqrt(-a)
            wcb = two * a * xb  +  b
            wcb = wcb / wd
            wcb = asin(wcb)
            wcb = -wcb / dsqrt(-a)
            value_241 = wcb - wca
         else if (rabs_use_stop) then
            stop "cowf_bronstein_241(): program stop B."
         end if
         !
      end function cowf_bronstein_241
      !
      !
      function cowf_bronstein_258(a,b,c,xa,xb)         result(value_258)
      !-----------------------------------------------------------------
      ! Returns the value of the integral No. 258 as defined by Bronstein 
      ! and Semendjajew : Taschenbuch der Mathematik (Leipzig 1983) over 
      ! the interval XA ... XB. 
      !-----------------------------------------------------------------
         !
         real(kind=dp), intent(in) :: a, b, c, xa, xb
         real(kind=dp) :: value_258
         !
	 real(kind=dp) :: wa, wb, wca, wcb, wd
         !
         if (c > zero) then
            wa  = a * xa * xa  +  b * xa  +  c
            wb  = a * xb * xb  +  b * xb  +  c
            wca = sqrt(c * wa)
            wca = c + c  +  b * xa  -  two * wca
            wca = wca / (two * xa)
            wca = log(wca) / sqrt(c)
            wcb = sqrt(c * wb)
            wcb = c + c + b * xb - two * wcb
            wcb = wcb / (two * xb)
            wcb = log(wcb) / sqrt(c)
            value_258 = wcb - wca
         else if(c < zero) then
            wd  = two * two * a * c  -  b * b
            if (rabs_use_stop   .and.   wd >= zero) then
               stop "cowf_bronstein_258(): program stop A."
            end if
            wd = sqrt(-wd)
            !
            wca = b * xa  +  c  +  c
            wca = wca / ( xa * wd )
            wca = dasin( wca ) / dsqrt( -c )
            wcb = b * xb  +  c  +  c
            wcb = wcb / ( xb * wd )
            wcb = dasin( wcb ) / dsqrt( -c )
            value_258 = wcb - wca
         else if (rabs_use_stop) then
            stop "cowf_bronstein_258(): program stop B."
         end if
         !
      end function cowf_bronstein_258
      !
   end function cowf_bronstein_260
   !
   !
   subroutine cowf_collect_input()
   !--------------------------------------------------------------------
   ! Collect and proceeds all input for the COWF program.
   !
   ! Calls: get_yes_stream().
   !--------------------------------------------------------------------
      !
      integer            :: arrow_position, blank_position, i, ios, ipc, j, l,&
                            kappa_c, level_No, mode, totalJ
      logical            :: yes
      real(kind=dp)      :: energy, en_lower, en_upper, delta_en, &
                            maximal_energy, wavelength, wavenumber
      character(len=1  ) :: parity
      character(len=20 ) :: string
      character(len=256) :: record
      !
      ! Open, check, load data from, and close the  .iso  file
      call file_get_isodat_grasp2k()
      !
      ! Get eigenvectors of the Atomic State Functions
      if (cowf_use_formatted_mix_file) then
         call file_get_eigenvectors_formatted(asf_bound)
      else
         call file_get_eigenvectors_unformatd(asf_bound)
      end if   
      !
      call input_energy_unit()
      !
      ! Print a brief summary of all bound ASF to select proper bound-state
      ! levels
      print *, "Summary of all (ionic) ASF to facilitate the selection of "//&
               " proper bound-state levels"
      print *, " to which additional continuum orbitals are to be coupled."
      print *, " "
      print *, "   i)  Level    J^P   Energy (a.u.)  "// &
               "  Energy ("//trim(energy_unit)//")"
      print *, "----------------------------------------------------"
      do  i = 1,asf_bound%noasf
         if (energy_inverse) then
            energy = energy_factor / asf_bound%asf(i)%energy
         else
            energy = energy_factor * asf_bound%asf(i)%energy
         end if
         print 2, i,asf_bound%asf(i)%level_No,                          &
	            angular_momentum_string(1*asf_bound%asf(i)%totalJ), &
	            asf_bound%asf(i)%parity,asf_bound%asf(i)%energy,    &
		    energy
       2 format(2x,i3,") ",i4,2x,a,a2,3x,1pe12.5,4x,1pe12.5)
      end do
      !
      print *, " "
      print *, "Continuum orbitals are generated in the mean field of an"// &
               " ASF with (given) Level number "
      print *, " which has total angular momentum J_c and parity P_c."
      print *, " A continuum orbitals with (given) Symmetry is coupled"//  &
               " to this bound ASF so, that a new ASF"
      print *, " with total angular momentum J_t and parity P_t"//       &
               " is generated."
      print *, " Enter (now in turn) the Level numbers, the Symmetry"// &
               " of the continuum orbital as well as J_t and P_t."
      print *, " "
      !
      ! Determine the symmetries and energies of the continuum orbitals
      ! to be generated
      number_of_selected_channels = 0
    3 print *, "Enter another Level_No, Symmetry, J_t, P_t, e.g."// &
                  " 3   p-  ==>  3/2 + , "
      print *, " 1  f  ==>  2 - , ...; <cr> if done."
      read (*, "(a)") record
      if (len_trim(record) > 0) then
         arrow_position = scan(record,"=")
         if (arrow_position == 0   .or.   &
	     record(arrow_position+1:arrow_position+2) /= "=>") then
            print *, "Unable to decode the coupling of an additional"// &
	             " continuum orbital; reenter ..."
            goto 3
         else 
            string   = adjustl(record(1:arrow_position-1))
            blank_position = scan(string," ")
            level_No = get_integer_from_string(string(1:blank_position-1))
	    print *, "level_No = ",level_No
            string   = adjustl(string(blank_position+1:))
	    print *, "string(1:2) = ",string(1:2)
	    kappa_c  = get_kappa_from_name(string(1:2))
	    print *, "kappa_c = ",kappa_c
            string   = adjustl(record(arrow_position+3:256))
            blank_position = scan(string," ")
	    totalJ   = get_dinteger_from_string(string(1:blank_position-1))
            parity   = string(blank_position+1:blank_position+1)
	    if (parity /= "+"  .and.   parity /= "-") then
               print *, "Unable to decode the parity of the channel;"// &
	                " reenter ..."
               goto 3
	    end if
	    !
	    ! Check the consistency of the input quantum numbers
	    l = angular_momentum_l(kappa_c)
	    j = angular_momentum_j(kappa_c)
	    do  i = 1,asf_bound%noasf
	       if (level_No == asf_bound%asf(i)%level_No) goto 4
	    end do
	    stop "cowf_collect_input(): program stop A."
	  4 if (asf_bound%asf(i)%parity == "+") then
	       ipc = 1
	    else 
	       ipc = -1
	    end if
	    if ((ipc * ((-1)**l) ==  1  .and.  parity == "-")  .or. &
	        (ipc * ((-1)**l) == -1  .and.  parity == "+")  .or. &
		(asf_bound%asf(i)%totalJ - J > totalJ)         .or. &
		(asf_bound%asf(i)%totalJ + J < totalJ)      )  then
	       print *, "Angular momenta and parities of the continuum "// &
	                " orbital and the selected atomic levels"
               print *, " are not consistent with each other; reenter ..."
	       goto 3
	    end if
	    number_of_selected_channels = number_of_selected_channels + 1
	    channel_selection(number_of_selected_channels)%asf      = i
	    channel_selection(number_of_selected_channels)%level_No = level_No
	    channel_selection(number_of_selected_channels)%kappa_c  = kappa_c
	    channel_selection(number_of_selected_channels)%totalJ   = totalJ
	    channel_selection(number_of_selected_channels)%parity   = parity
	    goto 3
	 end if
      else if (number_of_selected_channels == 0) then
         print *, "At least one continuum channel has to be selected; redo ..."
	 goto 3
      end if
      !
      ! Select the energies of the continuum orbitals
      maximal_energy = zero
      print *, "The energies of the continuum orbitals which are to be "// &
               " generated can be entered either"
      print *, "  by their individual energies or by an appropriate"//     &
               " interval and step size."
      print *, " The default is a list of individual energies; revise this ? " 
      yes = get_yes_stream()
      if (.not.yes) then
       5 print *, "Enter another (positive) energy in ",      &
                  trim(energy_unit)," for generating "//      &
                  " continuum orbitals; enter a negative number if done."
         read  *, energy
         if (energy > zero) then
            if (energy_inverse) then
               energy = energy_factor * energy
            else
               energy = energy /energy_factor
            end if
	    number_of_selected_energies = number_of_selected_energies + 1
	    if (number_of_selected_energies > 100) then
	       print *, "The maximum number of continuum orbital energies"// &
	                " is 100; program continues ..."
	       goto 7
	    end if
	    energy_selection(number_of_selected_energies) = energy
	    if (energy > maximal_energy)   maximal_energy = energy
	    goto 5
	 end if
      else
       6 print *, "Enter an interval of energies and a corresponding"// &
	          " stepsize for generating continuum orbitals; "
         print *, " E_lower  E_upper  delta-E:"
         read  *, en_lower, en_upper, delta_en
	 if (en_lower < zero   .or.   (en_upper - en_lower) < zero   .or. &
	     delta_en <= zero) then
	    print *, "All Energies must be greater than zero and "//&
                     " E_upper > E_lower; reenter ..."
	 else
	    energy = en_lower;   energy_selection(1) = energy
	    number_of_selected_energies = 1
	    do 
	       energy = energy + delta_en
	       if (energy > en_upper) exit
	       number_of_selected_energies = number_of_selected_energies + 1
	       if (number_of_selected_energies > 100) then
	          print *, "The maximum number of continuum orbital energies"//&
	                   " is 100; choose another interval or stepsize."
	          goto 6
	       end if
	       energy_selection(number_of_selected_energies) = energy
	    end do
         endif
      endif
      !
    7 continue
      !
      ! Determine the parameters controlling the radial grid
      call input_grid_parameters("standard")
      !
      ! Default speed of light
      c = c_vacuum
      !
      wavenumber = sqrt( two * maximal_energy + maximal_energy / (c*c) )
      wavelength = two * pi / wavenumber
      hp_grasp2k = wavelength / 100.0_dp
      n_grasp2k  = 50.0_dp / hp_grasp2k + 1000
      !
      ! Now 'overwrite' defaults only if required
      print *, "Modify default set-up and printout of the program ?"
      yes = get_yes_stream()
      if (.not.yes) goto 10
      !
    8 print *, "Select one mode below for normalization and to determine"//&
               " the phase shift d_phi:"
      print *, "  1  - Determine phase shift with respect to zero"       //&
               " potential; a short range potential is assumed for"      //&
               " normalization (default)."
      print *, "  2  - Use WKB method; P(r --> infinity) "               //&
               "~ cos(kr - (l+1)pi/2 + d_phi)."
      print *, "  3  - WKB normalization; phase shift with respect to"   //&
               " pure Coulomb potential."
      read(*,*,iostat=ios) mode
      if (ios /= 0   .or.   mode < 1   .or.   mode > 3) then
         print *, "Unable to interprete phase shift mode; redo ..."
         goto 8
      end if
      !
      cowf_phaseshift_wkb               = .false.
      cowf_phaseshift_zero_potential    = .false.
      cowf_phaseshift_coulomb           = .false.
      select case(mode)
      case(1)
         cowf_phaseshift_zero_potential = .true.
      case(2)
         cowf_phaseshift_wkb            = .true.
      case(3)
         cowf_phaseshift_coulomb        = .true.
      end select
      !
    9 print *, "Select a maximal radius R_max for the integration of"    //&
               " the continuum orbitals ?"
      yes = get_yes_stream()
      if (yes) then
         cowf_specify_r_max = .true.
         print *, "Enter R_max:"
         read *,  cowf_r_max
         if (cowf_r_max < zero   .or.   cowf_r_max  > 1.0e4) then
            print *, "Unable to interprete R_max: This must be in the"   //&
                     " interval  0. < R_max < 10000. (a.u.); redo ..."
            goto 9
         end if
         n_grasp2k  = (cowf_r_max + 5.0_dp) / hp_grasp2k + 1000
      end if
      !
      print *, "Print the CSF scheme each time a new one has been set up ?"
      yes = get_yes_stream()
      if (yes) cowf_print_always_csf_scheme = .true.
      !
      ! Determine the physical effects specifications
      print *, "The physical speed of light in atomic units is",c_vacuum,";"
      print *, " revise this value ?"
      yes = get_yes_stream()
      if (yes) then
         print *, "Enter the revised value:"
         read  *, c
      else
         c = c_vacuum
      endif
      !
      ! Modify the radial grid parameters
      call input_grid_parameters("modify")
      !
      ! Continue if the default set-up has not been modified
   10 continue
      !
      print *, "hp_grasp2k, n_grasp2k = ",hp_grasp2k, n_grasp2k
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
      ! Initialize and load the bound-state radial wavefunctions 
      call file_get_rwf("Enter the name of the bound-state Radial"// &
                        " WaveFunction File:",asf_bound,wave_bound,.false.)
      !
      ! Open a .csp Continuum SPinor output file
      call cowf_open_csp()
      !
   end subroutine cowf_collect_input
   !
   !
   function cowf_fk_r(k,r,ia,ib,csf_set)                   result(value)
   !--------------------------------------------------------------------
   ! Calculates the f^k_r(ab) angular coefficient for the direct potential
   ! in the MC Dirac-Fock eqquations. It follows mainly the paper of 
   ! Dyall et al. CPC 55, p. 425 (1989); Eqs. (24 - 26).          
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in)         :: ia, ib, k, r
      type(csf_basis), intent(in) :: csf_set
      real(kind=dp)               :: value
      !
      integer :: ja, jb
      !
      ja = angular_momentum_j(csf_set%subshell(ia)%kappa)
      jb = angular_momentum_j(csf_set%subshell(ib)%kappa)
      !
      ! Calculate f^k_r(ab) coefficients
      if (.false.   .and.                             &
         (csf_set%csf(r)%occupation(ia) == ja+1  .or. &
          csf_set%csf(r)%occupation(ib) == jb+1))   then
         if (k == 0   .and.   ia == ib) then
	    value = half * csf_set%csf(r)%occupation(ia) *                   &
	                  (csf_set%csf(r)%occupation(ia) - one)
         else if (k == 0   .and.   ia /= ib) then
	    value = csf_set%csf(r)%occupation(ia) *                          &
                    csf_set%csf(r)%occupation(ib)
         else if (k > 0   .and.  ia == ib) then  
            value = - half * (wigner_3j_symbol(ja,k+k,jb,1,0,-1)**2)
         else if (k > 0   .and.  ia /= ib) then 
            value = zero
         else if (rabs_use_stop) then
            stop "cowf_fk_r(): program stop A."
         end if
      else
         value = cowf_Vk_rs(r,r,k,csf_set%subshell(ia),csf_set%subshell(ib), &
                                  csf_set%subshell(ia),csf_set%subshell(ib))
      end if
      !
   end function cowf_fk_r
   !
   !
   subroutine cowf_generate_csp()
   !--------------------------------------------------------------------
   ! Controls the generation of continuum spinors in an atomic mean field 
   ! as given by the wave functions of the GRASP92 atomic structure 
   ! program. It starts from a converged run for one or several atomic 
   ! states in GRASP and calculates the (radial) continuum orbitals 
   ! for given energies, symmetries, and coupling conditions.   
   ! 
   ! Calls: spherical_Bessel_jL().
   !--------------------------------------------------------------------
      !
      integer          :: i, j, totalJ_c
      character(len=1) :: parity_c
      !
      integer       :: nw, nocsf
      real(kind=dp) :: energy
      type(nkappa)  :: subshell
      integer, dimension(1:asf_bound%csf_set%nocsf) :: ndx
      !
      ! Allocate for a "first time"; it is first dellocated before any usage
      allocate( cowf_csp%P(1:10), cowf_csp%Q(1:10) )
      !
      ! Create a list of all "continuum channels" as selected by the input
      ! and cycle through this list 
      ! (J_c, P_c  +  kappa_c  [energy_1,...,energy_n] ---->  J_t, P_t)
      call cowf_set_channels()
      energy = -one
      print *, "number_of_channels = ",number_of_channels
      do  i = 1,number_of_channels
         print *, "i, energy, channel(i)%energy = ", &
	           i, energy, channel(i)%energy 
         !
         ! Set the extended CSF including an continuum orbital
         totalJ_c = asf_bound%asf(channel(i)%asf)%totalJ
         parity_c = asf_bound%asf(channel(i)%asf)%parity
         call set_configuration_scheme(asf_bound%csf_set,asf_cont%csf_set, &
                        -1,channel(i)%kappa_c,totalJ_c,parity_c,           &
                        channel(i)%totalJ,channel(i)%parity,append=.false.,&
		        index=ndx)
         !
         nw = asf_cont%csf_set%nwshells
         if (rabs_use_stop  .and.  nw /= asf_bound%csf_set%nwshells + 1) then
            stop "cowf_generate_csp(): program stop A."
         end if
         ! Calculate the MCP coefficients for the current coupling scheme
         ! as well as the d_rs,  y_k(ab), and x_k(abcd) coefficients
	 nocsf = asf_cont%csf_set%nocsf
         call anco_calculate_csf_matrix(asf_cont%csf_set,1,nocsf,1,nocsf)
         !!x call mcp_generate_coefficients(1,nocsf,1,nocsf,asf_cont%csf_set)
         call cowf_set_drs_coefficients(channel(i)%asf,asf_cont%csf_set,ndx)
	 subshell = nkappa(-1,channel(i)%kappa_c)
         call cowf_set_yk_coefficients(subshell,asf_cont%csf_set)
         call cowf_set_xk_coefficients(subshell,asf_cont%csf_set)
         !
         ! Now iterate the continuum spinors for all energies of this channel
	 do  j = 1,channel(i)%number_of_energies
	    energy = channel(i)%energy(j)
            cowf_solve_homogeneous_eqn = .true.
            cowf_start_homogeneous     = .true.
            call cowf_iterate_csp(energy,nkappa(-1,channel(i)%kappa_c))
            !
            ! Append the orbital to the .csp file
            call cowf_append_csp()
	 end do
      end do
      !
      deallocate( cowf_csp%P, cowf_csp%Q )
      !
   end subroutine cowf_generate_csp
   !
   !
   function cowf_gk_r(k,r,ia,ib,csf_set)                    result(value)
   !--------------------------------------------------------------------
   ! Calculates the g^k_r(ab) angular coefficient for the exchange potential
   ! in the MC Dirac-Fock eqquations. It follows mainly the paper of 
   ! Dyall et al. CPC 55, p. 425 (1989); Eqs. (24 - 26).          
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in)         :: ia, ib, k, r
      type(csf_basis), intent(in) :: csf_set
      real(kind=dp)               :: value
      !
      integer :: ja, jb
      !
      ! Calculate g^k_r(ab) coefficients
      if (csf_set%csf(r)%occupation(ia) ==                                   &
          angular_momentum_j(csf_set%subshell(ia)%kappa) + 1  .or.           &
          csf_set%csf(r)%occupation(ib) ==                                   &
          angular_momentum_j(csf_set%subshell(ib)%kappa) + 1) then
	 !
	 ja = angular_momentum_j(csf_set%subshell(ia)%kappa)
	 jb = angular_momentum_j(csf_set%subshell(ib)%kappa)
	 value = - csf_set%csf(r)%occupation(ia) *                           &
	           csf_set%csf(r)%occupation(ib) *                           &
		   (wigner_3j_symbol(ja,k+k,jb,1,0,-1)**2)
      else if (csf_set%csf(r)%occupation(ia) <                               &
               angular_momentum_j(csf_set%subshell(ia)%kappa) + 1  .and.     &
               csf_set%csf(r)%occupation(ib) <                               &
               angular_momentum_j(csf_set%subshell(ib)%kappa) + 1) then
         value = cowf_Vk_rs(r,r,k,csf_set%subshell(ia),csf_set%subshell(ib), &
                                  csf_set%subshell(ib),csf_set%subshell(ia))
      else if (rabs_use_stop) then
         print *, "k = ",k
         stop "cowf_gk_r(): program stop A."
      end if
      !
   end function cowf_gk_r
   !
   !
   subroutine cowf_integrate_csp(csp)
   !--------------------------------------------------------------------
   ! Integrates the radial DF equations for the continuum spinor csp.
   ! This routine first set up the 'potential' functions A_p, chi_p, and
   ! A_q, chi_q which contain the left-  and right-hand side of the
   ! DF equations for the large (L) and small (S) components, respectively.
   !
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(orbital_function), intent(inout) :: csp
      !
      integer       :: i, ia, ib, j, mtp
      real(kind=dp) :: Q0, dx, accy, wavenumber, wavelength
      real(kind=dp), dimension(1:6)         :: P, Q
      real(kind=dp), dimension(1:n_grasp2k) :: chi_p, chi_q, P_p, Q_p
      !
      integer, parameter               :: steps = 10
      real(kind=dp), dimension(4)      :: rz, yz, chi_pz, chi_qz, Pz, Qz, &
                                          P_pz, Q_pz
      real(kind=dp), dimension(-5:200) :: rx, yx, chi_px, chi_qx, Px, Qx, &
                                          P_px, Q_px
      !
      ! First set up the left-  and right-hand sides of the DF equations
      mtp = csp%mtp
      chi_p(2:mtp) = cowf_xp_pot(2:mtp); chi_q(2:mtp) = cowf_xq_pot(2:mtp)
      !
      do  i = 1,number_of_epsilon_ab
         if (cowf_eab(i)%a == asf_cont%csf_set%nwshells  .and. &
             cowf_eab(i)%b <= wave_bound%number_of_rwf   .and. &
             cowf_eab(i)%a /= cowf_eab(i)%b)  then
            ia  = cowf_eab(i)%a;   ib = cowf_eab(i)%b
            mtp = min(csp%mtp, wave_bound%rwf(ib)%mtp)
	    chi_p(2:mtp) = chi_p(2:mtp) + (r_grasp2k(2:mtp) *   &
	                   cowf_eab(i)%epsilon *                &
	                   wave_bound%rwf(ib)%P(2:mtp)) /       &
			   (c * cowf_gen_occupation(ia) )
	    chi_q(2:mtp) = chi_q(2:mtp) + (r_grasp2k(2:mtp)  *  &
	                   cowf_eab(i)%epsilon *                &
	                   wave_bound%rwf(ib)%Q(2:mtp)) /       &
			   (c * cowf_gen_occupation(ia) )
	 else if (rabs_use_stop) then
	    stop "cowf_integrate_csp(): program stop A."
	 end if
      end do
      !
      if (debug_chi_pq) then
         write(99,*) " "
         write(99,*) "Direct potential and chi_p(r), chi_q(r) :"
         write(99,*) "   i)         R(i)              Y(i)         "// &
	             "  chi_p(i)          chi_q(i)"
         do  i = 1,mtp,100
            write(99,1) i, r_grasp2k(i), cowf_y_pot(i), &
	                chi_p(i),chi_q(i)
         end do
      end if
      !
      ! Generate 6 points to start the integration
      call cowf_start_integration(csp,csp%pz,Q0,P,Q)
      !
      wavenumber = sqrt( two * (-csp%energy) + (-csp%energy) / (c*c) )
      wavelength = two * pi / wavenumber
      !! print *, "csp%orbital%kappa, energy, wavelength = ",               &
      !!           csp%orbital%kappa, csp%energy, wavelength
      !
      ! Calculate derivatives for the first 6 points
      csp%P(1:6) = P(1:6);   csp%Q(1:6) = Q(1:6)
      P_p(1)     = zero;     Q_p(1)     = zero
      P_p(2:6)   = - (csp%orbital%kappa * csp%P(2:6)/r_grasp2k(2:6)      &
                       - (two * c - (csp%energy / c) +                   &
	 	          (cowf_y_pot(2:6) / (c*r_grasp2k(2:6)))         &
	 	         ) * csp%Q(2:6) ) - chi_p(2:6) / r_grasp2k(2:6)
      Q_p(2:6)   = - (- csp%orbital%kappa * csp%Q(2:6)/r_grasp2k(2:6)    &
                       + (- (csp%energy / c) +                           &
                          (cowf_y_pot(2:6)/(c*r_grasp2k(2:6)))           &
	 	         ) * csp%P(2:6) ) + chi_q(2:6) / r_grasp2k(2:6)
      !
      do  i = 6,csp%mtp-1
         ! P_p(i) = - (csp%orbital%kappa * csp%P(i)/r_grasp2k(i)         &
         !                - (two * c - (csp%energy / c) +                &
	 ! 	           (cowf_y_pot(i) / (c*r_grasp2k(i)))            &
	 ! 	          ) * csp%Q(i) ) - chi_p(i) / r_grasp2k(i)
         ! Q_p(i) = - (- csp%orbital%kappa * csp%Q(i)/r_grasp2k(i)       &
         !                + (- (csp%energy / c) +                        &
         !                   (cowf_y_pot(i)/(c*r_grasp2k(i)))            &
	 ! 	          ) * csp%P(i) ) + chi_q(i) / r_grasp2k(i)
         ! csp%P(i+1) = csp%P(i) + P_p(i) * &
         !                         1.0_dp*(r_grasp2k(i+1)-r_grasp2k(i-1))/two
         ! csp%Q(i+1) = csp%Q(i) + Q_p(i) * &
         !                         1.0_dp*(r_grasp2k(i+1)-r_grasp2k(i-1))/two
         !
         ! Prepare integration for one radial step
         dx = (r_grasp2k(i+1)-r_grasp2k(i)) / steps
         do  j = -5,steps
            rx(j) = r_grasp2k(i) + j*dx
         end do
         rz(1:4)     = r_grasp2k(i-2:i+1);   yz(1:4)     = cowf_y_pot(i-2:i+1)
         chi_pz(1:4) = chi_p(i-2:i+1);       chi_qz(1:4) = chi_q(i-2:i+1)
         Pz(1:3)     = csp%P(i-2:i);         Qz(1:3)     = csp%Q(i-2:i)
         P_pz(1:3)   = P_p(i-2:i);           Q_pz(1:3)   = Q_p(i-2:i)
         !
         accy = 1.0e-2_dp
         do  j = -2,0
            !
            call interpolation_lagrange(rz,Pz,   rx(j),Px(j),   3)
            call interpolation_lagrange(rz,Qz,   rx(j),Qx(j),   3)
            call interpolation_lagrange(rz,P_pz, rx(j),P_px(j), 3)
            call interpolation_lagrange(rz,Q_pz, rx(j),Q_px(j), 3)
         end do
         do  j = -2,steps
            !
            call interpolation_lagrange(rz,yz,    rx(j),yx(j),    4)
            call interpolation_lagrange(rz,chi_pz,rx(j),chi_px(j),4)
            call interpolation_lagrange(rz,chi_qz,rx(j),chi_qx(j),4)
         end do
         !
         Px(0) = csp%P(i);   Qx(0) = csp%Q(i) 
         do  j = 0,steps-1
            if (.false.) then
               ! 
               ! Simplest stepwise integration
               P_px(j) = - (csp%orbital%kappa * Px(j)/rx(j)                   &
                         - (two * c - (csp%energy / c) +  (yx(j) / (c*rx(j))) &
	 	            ) * Qx(j) ) - chi_px(j) / rx(j)
               Q_px(j) = - (- csp%orbital%kappa * Qx(j)/rx(j)                 &
                         + (- (csp%energy / c) +  (yx(j)/(c*rx(j)))           &
	 	            ) * Px(j) ) + chi_qx(j) / rx(j)
               Px(j+1) = Px(j) + P_px(j) * (rx(j+1)-rx(j))
               Qx(j+1) = Qx(j) + Q_px(j) * (rx(j+1)-rx(j))
            else if (.true.) then
               !
               ! Adams-Bashfort-Moulton scheme; predictor
               Px(j+1) = Px(j) + dx/12.0_dp * (23.0_dp*P_px(j)                &
                          - 16.0_dp*P_px(j-1) + 5.0_dp*P_px(j-2))
               Qx(j+1) = Qx(j) + dx/12.0_dp * (23.0_dp*Q_px(j)                &
                          - 16.0_dp*Q_px(j-1) + 5.0_dp*Q_px(j-2))
               !
               P_px(j+1) = - (csp%orbital%kappa * Px(j+1)/rx(j+1)             &
                         - (two * c - (csp%energy / c) + (yx(j+1)/(c*rx(j+1)))&
	 	              ) * Qx(j+1) ) - chi_px(j+1) / rx(j+1)
               Q_px(j+1) = - (- csp%orbital%kappa * Qx(j+1)/rx(j+1)           &
                              + (- (csp%energy / c) +  (yx(j+1)/(c*rx(j+1)))  &
	 	              ) * Px(j+1) ) + chi_qx(j+1) / rx(j+1)
               ! 
               ! Corrector
               Px(j+1) = Px(j) + dx/12.0_dp * (5.0_dp*P_px(j+1)               &
                          + 8.0_dp*P_px(j) - P_px(j-1))
               Qx(j+1) = Qx(j) + dx/12.0_dp * (5.0_dp*Q_px(j+1)               &
                          + 8.0_dp*Q_px(j) - Q_px(j-1))
               !
               P_px(j+1) = - (csp%orbital%kappa * Px(j+1)/rx(j+1)             &
                         - (two * c - (csp%energy / c) + (yx(j+1)/(c*rx(j+1)))&
	 	              ) * Qx(j+1) ) - chi_px(j+1) / rx(j+1)
               Q_px(j+1) = - (- csp%orbital%kappa * Qx(j+1)/rx(j+1)           &
                              + (- (csp%energy / c) +  (yx(j+1)/(c*rx(j+1)))  &
	 	              ) * Px(j+1) ) + chi_qx(j+1) / rx(j+1)
            end if   
         end do
         !
         csp%P(i+1) = Px(steps);     csp%Q(i+1) = Qx(steps)
         P_p(i+1)   = P_px(steps);   Q_p(i+1)   = Q_px(steps)
      end do
      !
      ! Print debug information if required
      if (debug_integrate_csp) then
         write(99,*) " "
         write(99,*) "Newly integrated continuum spinor:"
         write(99,*) "   i)        R(i)              P(i)              Q(i)"//&
	             "             P_p(i)            Q_p(i)"
         do  i = 1,10,1
            write(99,1) i, r_grasp2k(i), csp%P(i), csp%Q(i), P_p(i), Q_p(i)
         end do
         do  i = 11,mtp,50
            write(99,1) i, r_grasp2k(i), csp%P(i), csp%Q(i), P_p(i), Q_p(i)
         end do
       1 format(i5,")  ",7(1pe16.9,2x))
      end if
      !
   end subroutine cowf_integrate_csp
   !
   !
   subroutine cowf_iterate_csp(energy,subshell)
   !--------------------------------------------------------------------
   ! Controls the iteration of the radial Dirac-Fock equations for a 
   ! (free electron) continuum spinor with energy > 0. The maximal number
   ! of iterations is given by the parameter cowf_maximal_iteration.
   ! No off-diagonal terms are included in the current version.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      real(kind=dp), intent(in) :: energy
      type(nkappa), intent(in)  :: subshell
      !
      logical, save       :: first_call = .true.
      real(kind=dp), save :: energy_sav = 0.0_dp
      !
      integer       :: i, ii, i1, i2, iteration, kappa, l, mtp, mtp_max, &
                       nr0, pqn_b, pqn_c, imax
      real(kind=dp) :: d, deviation, max_deviation,  max_value, value,   &
                       wavelength, wavenumber, norm_sum, norm_f, phase_sum,  &
                       x1val, x1loc, x2val, x2loc, y, cn
      real(kind=dp), dimension(500)         :: test_norm, test_phase
      real(kind=dp), dimension(1:n_grasp2k) :: tax
      !
      ! Allocate storage for the nuclear potential as well as for the
      ! direct and exchange potentials
      if (.not.allocated(cowf_nuc_pot)) then
         allocate( cowf_nuc_pot(1:n_grasp2k) )
	 cowf_nuc_pot(1) = zero
	 do  i = 2,n_grasp2k
	    cowf_nuc_pot(i) = nuclear_potential(r_grasp2k(i))
	 end do
      end if
      if (.not.allocated(cowf_y_pot))  allocate( cowf_y_pot(1:n_grasp2k) )
      if (.not.allocated(cowf_xp_pot)) allocate( cowf_xp_pot(1:n_grasp2k) )
      if (.not.allocated(cowf_xq_pot)) allocate( cowf_xq_pot(1:n_grasp2k) )
      cowf_y_pot = zero;   cowf_xp_pot = zero;   cowf_xq_pot = zero   
      !
      csp_new%orbital  = subshell
      csp_new%energy   = -energy
      csp_new%pz       = 1.0e-0_dp
      csp_new%phase    = zero
      select case(nuclear_model)
      case("point")
         csp_new%gamma = gamma_r(subshell%kappa)
      case default
         csp_new%gamma = abs(subshell%kappa)
      end select
      !
      ! Determine the radius r_bound of the bound-state density and compare the
      ! radial step size at r_bound with the 'wavelength' of the current csp;
      ! first estimate 'asymtotic wavelength' for given energy
      wavenumber = sqrt( two * energy + energy / (c*c) )
      wavelength = sqrt( two * pi * pi / energy )
      mtp_max = 0
      do  i = 1,wave_bound%number_of_rwf
         mtp_max = max( mtp_max, wave_bound%rwf(i)%mtp )
      end do
      mtp_max = mtp_max + cowf_extent_mtp
      !
      if (first_call) then
         print 1, mtp_max,r_grasp2k(mtp_max), &
                  r_grasp2k(mtp_max)-r_grasp2k(mtp_max-1)
       1 format( " Bound-state electron density vanishes at R(",i5,") = ", &
                 f9.4," a.u. at a step size Delta-R = ",f9.4," a.u. ")
         first_call = .false.
      end if
      !
      if (energy /= energy_sav) then
         energy_sav = energy
         !! print 4, energy,wavelength
       4 format(" The step size at large R has to be compared to the",     &
                " asymtotic wavelength of the continuum electron ",        &
              /," Lambda(",1pe10.4," Hartree) = ",1pe10.4," a.u.")
      end if
      !
      if (wavelength < cowf_min_points_per_wavelength * &
                       (r_grasp2k(mtp_max)-r_grasp2k(mtp_max-1))) then
         print *, "*** WARNING: Radial grid might be inappropriate"// &
                  "for the given wavelength. ***"
      end if
      !
      csp_new%mtp = mtp_max
      allocate( csp_new%P(1:mtp_max), csp_new%Q(1:mtp_max) )
      csp_new%P(1:mtp_max) = zero;    csp_new%Q(1:mtp_max) = zero
      !
      ! Assign a second csp to keep the last iteration
      allocate( csp_sav%P(1:mtp_max), csp_sav%Q(1:mtp_max) )
      csp_sav%P(1:mtp_max) = zero;    csp_sav%Q(1:mtp_max) = zero
      !
      ! The iteration of the continuum spinors cycles over three steps
      !   (1) Calculate the direct and exchange potentials
      !   (2) Calculate Lagrange multipliers (not in the current version)
      !   (3) Integrate the radial DF equations
      !
      print *, " "
      iteration = 0
    2 iteration = iteration + 1
      !
      ! print *, "cowf_iterate_csp - b"
      call cowf_set_direct_potential(csp_new)
      number_of_epsilon_ab = 0
      if (.not.cowf_solve_homogeneous_eqn) then
         if (iteration > 1) call cowf_set_exchange_potential(csp_new)
         if (cowf_use_lagrange_multipliers   .and.   iteration > 1) then
            call cowf_set_lagrange_multipliers(subshell,csp_new)
         end if
      end if
      ! print *, "cowf_iterate_csp - c"
      !
      if (debug_yx_potentials) then
         write(99,*) " "
         write(99,*) "Direct and exchange potentials after ",iteration-1, &
	             " iterations:"
         write(99,*) "   i)         R(i)              Y(i)         "// &
	             "   X_p(i)            X_q(i)"
         do  i = 1,mtp_max,100
            write(99,3) i, r_grasp2k(i), cowf_y_pot(i), &
	                cowf_xp_pot(i),cowf_xq_pot(i)
         end do
       3 format(i5,")  ",4(1pe16.9,2x))
      end if
      !
      ! Now integrate the continuum spinor 
      call cowf_integrate_csp(csp_new)
      ! print *, "cowf_iterate_csp - d"
      !
      ! Carry out a Schmidt orthogonalization for each iteration if required
      if (cowf_schmidt_always) then 
         call cowf_schmidt_orthogonalization(csp_new)
      end if
      ! print *, "cowf_iterate_csp - e"
      !
      ! Compare with previous iteration which is kept in csp_sav
      max_deviation = zero;   max_value = zero
      do  i = 1,csp_new%mtp
         deviation = abs(csp_new%P(i)-csp_sav%P(i)) + &
	             abs(csp_new%Q(i)-csp_sav%Q(i))
         value     = abs(csp_new%P(i))+abs(csp_new%Q(i))
         if (deviation > max_deviation) max_deviation = deviation
         if (value >  max_value) max_value = value
      end do
      !
      print "(a,i5,a,2(1pe14.8,1x))",                                 &
              " Maximal deviation after ",iteration," iterations = ", &
	        max_deviation,max_deviation/max_value
      if (debug_iterate_csp) then
         write(99,*) "Maximal deviation after ",iteration," iterations = ", &
	             max_deviation
      end if
      !
      ! Use a damping factor for higher iterations
      d = 0.05_dp
      if (iteration > 1) then
         csp_sav%P(:) = d * csp_sav%P(:) + (one-d) * csp_new%P(:)
         csp_sav%Q(:) = d * csp_sav%Q(:) + (one-d) * csp_new%Q(:)
      else
         csp_sav%P(:) = csp_new%P(:);   csp_sav%Q(:) = csp_new%Q(:)
      end if
      !
      if (max_deviation/max_value < 10*accy_grasp2k) then
      else if (iteration < cowf_maximal_iteration) then
         goto 2
      else
         print *, "cowf_iterate_csp(): Convergence failure for continuum"// &
	          " spinor |epsilon,",orbital_symmetry(subshell%kappa),     &
		  ">, the max. deviation after ",cowf_maximal_iteration,    &
		  " iterations is ",max_deviation,"."
      end if
      !
      if (debug_normalize_csp) then
         ii = 0
         do  i = csp_new%mtp-2,csp_new%mtp-20,-4
	    ii = ii + 1
            call cowf_normalize_wkb(csp_new,i,.false.,test_norm(ii))
	    test_phase(ii) = csp_new%phase
	    !!x print *, "i, r(i), test_norm(i), test_phase(i) = ",  &
	    !!x           i, r_grasp2k(i), test_norm(ii), test_phase(ii)
	 end do
      end if
      !
      ! Carry out a Schmidt orthogonalization before the normalization
      if (cowf_schmidt_final) then 
         call cowf_schmidt_orthogonalization(csp_new)
      end if
      !
      ! Normalize the continuum spinor and determine the phase shift
      if (cowf_norm_wkb) then
         call cowf_normalize_wkb(csp_new,csp_new%mtp-40,.true.)
         csp_new%phase = mod(csp_new%phase,two*pi) 
      else if (cowf_norm_wkb_old) then
         !!x print *, "+++++ Use previous normalization scheme. +++++"
         do  i = 100,csp_new%mtp
            if (csp_new%P(i-1) - csp_new%P(i-2) > zero  .and. &
      	     	csp_new%P(i)   - csp_new%P(i-1) < zero  .and. &
      	     	i > 0.8* mtp_max) then
      	     	!! r_grasp2k(i) > 6.0_dp) then
      	       nr0 = i
            end if
         end do
         if (nr0 == 0) then
            stop "cowf_iterate_csp - stop B"
         end if
         !!x print *, "nr0, mtp_max, csp_new%mtp = ",nr0, mtp_max, csp_new%mtp
         call cowf_normalize_wkb(csp_new,nr0,.true.)
         csp_new%phase = csp_new%phase +				    &
      	     	     (one+angular_momentum_l(csp_new%orbital%kappa))/two * pi
         csp_new%phase = csp_new%phase - wavenumber*1.0e4_dp + ten*ten*two*pi
         y = cowf_y_pot(nr0) * (energy+c*c) / (c*c*wavenumber)
         csp_new%phase = csp_new%phase - y*log(two*wavenumber*1.0e4_dp)
         csp_new%phase = mod(csp_new%phase,two*pi) 
      else if (cowf_norm_nonrel) then
         !!x print *, "+++++ Use nonrelativistic scheme for normalization. +++++"
	 nr0 = csp_new%mtp-100
         call cowf_normalize_nonrel(csp_new,nr0,cowf_y_pot(nr0))
      else
         stop "cowf_iterate_csp - stop A: No normalization scheme chosen."
      end if
      !
      ! Set result into cowf_csp
      deallocate( cowf_csp%P, cowf_csp%Q )
      allocate( cowf_csp%P(1:mtp_max), cowf_csp%Q(1:mtp_max) )
      cowf_csp%P       = csp_new%P;        cowf_csp%Q     = csp_new%Q
      cowf_csp%orbital = csp_new%orbital;  cowf_csp%gamma = csp_new%gamma
      cowf_csp%energy  = csp_new%energy;   cowf_csp%pz    = csp_new%pz
      cowf_csp%phase   = csp_new%phase;    cowf_csp%mtp   = csp_new%mtp
      !
      !
      ! Deallocate all working arrays
      !
      deallocate( csp_new%P )
      deallocate( csp_new%Q )
      !!!!
      deallocate( csp_sav%P ) 
      deallocate( csp_sav%Q )
      !
      ! Test Lagrange multipliers
      ! -------------------------
      if (debug_lagrange_epsilon_ab) then
      !! if (.false.) then
         write(99,*) " "
         write(99,*) "Evaluation of overlap integrals:"
         do  i = 1,wave_bound%number_of_rwf
            if (wave_bound%rwf(i)%orbital%kappa == cowf_csp%orbital%kappa) then
               kappa = wave_bound%rwf(i)%orbital%kappa
               pqn_b = wave_bound%rwf(i)%orbital%n
               pqn_c = cowf_csp%orbital%n
	       !	   
               ! Determine the maximum tabulation point for the integrand
               mtp = min(wave_bound%rwf(i)%mtp, cowf_csp%mtp)
               !
               ! Tabulate the integrand as required for subroutine quad-grasp2k; 
               ! the value at the first tabulation point is arbitrary
               tax(1) = zero
               do  l = 2,mtp
                  tax(l) = (wave_bound%rwf(i)%P(l) * cowf_csp%P(l)                &
                        + wave_bound%rwf(i)%Q(l) * cowf_csp%Q(l))*rp_grasp2k(l)
               end do
               !
               ! Perform the quadrature
               value = quad_grasp2k(tax,mtp)
               !
	       write(99,*) "  "//orbital_symmetry(kappa)//               &
		           ":  < n_c = ",pqn_c," | n_b = ",pqn_b,"> = ", &
			   value
	       write(*,*)  "  "//orbital_symmetry(kappa)//               &
		           ":  < n_c = ",pqn_c," | n_b = ",pqn_b,"> = ", &
			   value
            end if
	    do  ii = 1,wave_bound%number_of_rwf 
            if (wave_bound%rwf(i)%orbital%kappa == wave_bound%rwf(ii)%orbital%kappa) then
               kappa = wave_bound%rwf(i)%orbital%kappa
               pqn_b = wave_bound%rwf(i)%orbital%n
               pqn_c = wave_bound%rwf(ii)%orbital%n
	       !	   
               ! Determine the maximum tabulation point for the integrand
               mtp = min(wave_bound%rwf(i)%mtp, wave_bound%rwf(ii)%mtp)
               !
               ! Tabulate the integrand as required for subroutine quad-grasp2k; 
               ! the value at the first tabulation point is arbitrary
               tax(1) = zero
               do  l = 2,mtp
                  tax(l) = (wave_bound%rwf(i)%P(l) * wave_bound%rwf(ii)%P(l)                &
                        + wave_bound%rwf(i)%Q(l) * wave_bound%rwf(ii)%Q(l))*rp_grasp2k(l)
               end do
               !
               ! Perform the quadrature
               value = quad_grasp2k(tax,mtp)
               !
	       write(99,*) "  "//orbital_symmetry(kappa)//               &
		           ":  < n_c = ",pqn_c," | n_b = ",pqn_b,"> = ", &
			   value
	       write(*,*)  "  "//orbital_symmetry(kappa)//               &
		           ":  < n_c = ",pqn_c," | n_b = ",pqn_b,"> = ", &
			   value
            end if
	    end do
         end do
      end if
      !
   end subroutine cowf_iterate_csp
   !
   !
   subroutine cowf_normalize()
   !--------------------------------------------------------------------
   ! This procedure collect some of the previous code that has been used 
   ! for normalizing the continuum orbitals. This code is presently not 
   ! used in the program.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
   end subroutine cowf_normalize
   !  
   !  
   subroutine cowf_normalize_nonrel(csp,nrp,za)
   !--------------------------------------------------------------------
   ! Normalizes the continuum orbital csp by using a non-relativistic scheme
   ! for the large component. The amplitude and phase is determined from
   ! the large component at the points nrp and nrp+1. 
   !  
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(orbital_function), intent(inout) :: csp
      integer, intent(in)                   :: nrp
      real(kind=dp), intent(in)             :: za    ! effective charge at nrp
      !
      integer                      :: l, l1, j
      real(kind=dp)                :: ac, as, an, anorm, det, dw1, dw2, e, &
                                      phase, se, x1, x2
      real(kind=dp), dimension(50) :: f1, f2, g1, g2, fc1, fc2, gc1, gc2
      !
      phase = zero
      anorm = zero
      !
      l   = angular_momentum_l(csp%orbital%kappa)
      l1  = l + 1
      e   = - csp%energy * two
      se  = sqrt(e)
      dw1 = csp%P(nrp)
      dw2 = csp%P(nrp+1)
      an  = - za/se
      x1  = r_grasp2k(nrp)   * se
      x2  = r_grasp2k(nrp+1) * se
      call coul_cern_kabachnik(x1,an,l,l,f1,fc1,g1,gc1,1.0e-6_dp,100.0_dp)
      call coul_cern_kabachnik(x2,an,l,l,f2,fc2,g2,gc2,1.0e-6_dp,100.0_dp)
      !
      det = f1(l1)*g2(l1) - f2(l1)*g1(l1)
      if(det == zero) goto 10
      ac = (dw1*g2(l1) - dw2*g1(l1)) / det
      as = (dw2*f1(l1) - dw1*f2(l1)) / det
      anorm = ac*ac + as*as
      phase = pi / two
      !
      if (ac == zero) goto 9
      phase = phase + atan(as/ac) + facouz_kabachnik(e,l,za)
      !
      if (ac < zero) phase = phase + pi   
    9 anorm = one / sqrt(anorm*pi*se)
      !
   10 continue
      csp%P(:)  = anorm * csp%P(:) * sqrt(two)
      csp%Q(:)  = anorm * csp%Q(:) * sqrt(two)
      csp%phase = phase
      !
      !!x print *, "*** phase = ",phase
      !
      contains
         !
         function facouz_kabachnik(e,l,z)                   result(facouz)
         !--------------------------------------------------------------
         ! Calculates the `Coulomb' phase sigma_l.
         !--------------------------------------------------------------
         !
         integer, intent(in)       :: l
	 real(kind=dp), intent(in) :: e, z
	 real(kind=dp)             :: facouz
	 !
	 integer                   :: k, m
	 real(kind=dp)             :: al, be, gam
	 !
         gam = -z/sqrt(e)
         m   = 51-l
         do 19  k=1,m
            al = atan(gam/(52-k))
            if(k == 1) go to 18
            facouz = facouz - al
            go to 19
         18 be = sqrt(gam*gam+(52-k)**two)
            facouz = al*50.5_dp+gam*(log(be)- one) + (-sin(al)/12.0_dp &
	             + sin(three*al)/(360.0_dp*be*be))/be
      19 continue
         !
         end function facouz_kabachnik
         !
      !
   end subroutine cowf_normalize_nonrel
   !
   !
   subroutine cowf_normalize_wkb(csp,nr0,do_normalize,norm_factor)
   !--------------------------------------------------------------------
   ! Normalizes a continuum orbital in an asymtotic coulomb potential. 
   ! This routine uses a WKB-approach which allows for to estimate
   ! the asymtotic behaviour of the continuum orbitals at rather an
   ! arbitrary point r0 in the asymtotic coulomb potential; of course,
   ! r0 should be beyond the extent of the bound-electron charge density. 
   ! The WKB-method applied by this routine is appropriate particularly 
   ! for highly-ionized  atoms with a (remaininn and screened) nuclear 
   ! potential outside the bound-electron charge density.         
   ! This routine follows basically the paper of W. Ong and A. Russek, 
   ! Phys. Rev. A17, 120-124 (1979) which is refered to as I below. 
   ! It should be  noted, however, that the equations from this paper 
   ! needed to be adapted to the GRASP formalism and are not identical 
   ! with the form given in I.
   ! Normalization of continuum spinors is carried out 'per unit energy'
   ! for given (asymptotic) energy as appropriate for the description 
   ! of most ionization  and scattering processes.
   !
   ! Calculates 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(orbital_function), intent(inout) :: csp
      integer, intent(in)                   :: nr0
      logical, intent(in)                   :: do_normalize
      real(kind=dp), optional, intent(out)  :: norm_factor
      !
      integer       :: i
      real(kind=dp) :: a, ae, b, be, cnorm, d, &
                       e, es, phase, ps, te, u, ypot_p,      &
		       ypot_p1, ypot_p2, wa, wb, wc
      real(kind=dp) :: aa, aq, bb, bq, cc, dq, dchi, fq, xa, xb, snphi              
      !
      ! Radius r0 must be chosen far enough from the origin to allow for 
      ! the application of the WKB-approach (I, Eq.(6)).
      b   = (cowf_y_pot(nr0) / r_grasp2k(nr0)) - csp%energy
      b   = b / c
      e   = c + c + b
      !
      ! Calculate the derivative of the potential and the wave function at r0
      ypot_p1 = (cowf_y_pot(nr0)   - cowf_y_pot(nr0-1)) / &
                (r_grasp2k(nr0)    - r_grasp2k(nr0-1))
      ypot_p2 = (cowf_y_pot(nr0+1) - cowf_y_pot(nr0)) / &
                (r_grasp2k(nr0+1)  - r_grasp2k(nr0))
      ypot_p  = (ypot_p1 + ypot_p2) / (two * c)
      !
      if (rabs_use_stop   .and.  abs(ypot_p1 - ypot_p2)  > 8.0e-3_dp) then
	 print *, "cowf_normalize_wkb(): abs(ypot_p1 - ypot_p2) = ", &
	                                 abs(ypot_p1 - ypot_p2)
         stop "cowf_normalize_wkb(): program stop A."
      endif
      !
      be = csp%energy / c
      ae = c + c - be
      ps = -csp%orbital%kappa * csp%P(nr0) / r_grasp2k(nr0) +    &
           (ae + cowf_y_pot(nr0)/r_grasp2k(nr0)) * csp%Q(nr0) - &
	   cowf_xp_pot(nr0)/r_grasp2k(nr0)
      es = ypot_p / r_grasp2k(nr0) - cowf_y_pot(nr0) /           &
           (c*r_grasp2k(nr0)*r_grasp2k(nr0))
      !
      ! Calculate the U-function (I, Eq. 12) at point r0
      wa = b * e  -  csp%orbital%kappa * (csp%orbital%kappa + one) /  &
           (r_grasp2k(nr0)*r_grasp2k(nr0))
      wb = es * (e + b) +                                             &
           two * csp%orbital%kappa * (csp%orbital%kappa + one) /      &
           (r_grasp2k(nr0)*r_grasp2k(nr0)*r_grasp2k(nr0))
      wc = two * two * wa
      !
      u  = half * es * csp%P(nr0) / e - ps - wb * csp%P(nr0) / wc
      u  = u / sqrt(wa)
      !
      ! Calculate the A - amplitude (I, Eq. 13) at the point r0;
      ! it should be independent from r0
      a  = (csp%P(nr0) * csp%P(nr0)  +  u * u) * sqrt(wa) / e
      a  = sqrt(a)
      !!x print *, "*** a: wa, wb, wc, a = ",wa, wb, wc, a
      !
      ! Calculate the normalization constant
      te = - csp%energy
      d  = pi *c*c * sqrt(two * te  +  (te * te) / (c * c))
      d  = (two * c * c  +  te) / d
      d  = dsqrt(d)
      !
      cnorm = sqrt(te / (two * c * c  +  te))
      cnorm = sqrt(cnorm) * d / a
      !
      phase = sqrt(sqrt(a)) * csp%P(nr0) / (a * sqrt(e))
      phase = atan(u/csp%P(nr0)) 
      !
      !
      phase = csp%orbital%kappa * (csp%orbital%kappa + one)
      phase = b * e  -  phase / (r_grasp2k(nr0)*r_grasp2k(nr0))
      phase = sqrt( phase )
      phase = sqrt( phase ) / a
      phase = phase / sqrt( e )
      phase = phase * csp%P(nr0)
      phase = acos( phase )
      snphi = u / a
      if (snphi .lt. zero) then
         phase = two * pi  - phase
      endif
      !
      !
      ! Calculate the phase-integral for chi(r) - chi(r0), I, Eq (9)
      aq = - csp%energy / c
      bq = cowf_y_pot(nr0) / c
      dq = c + c
      fq = csp%orbital%kappa * (csp%orbital%kappa + one)
      !
      aa = aq * dq  +  aq * aq
      bb = two * aq * bq  +  dq * bq
      cc = bq * bq  -  fq
      xa = r_grasp2k(nr0)
      xb = 1.0e4_dp
      dchi  = cowf_bronstein_260(aa,bb,cc,xa,xb)
      phase = phase + dchi
      !
      if (debug_normalize_csp) then
         write(99,*)  " "
         write(99,*)  "Normalization at grid point nr0 = ",nr0,":"
         write(99,*)  "  kappa, energy, te)  = ",csp%orbital%kappa, &
	                                        csp%energy, te
         write(99,*)  "  ps, ypot_p, d, u, a = ",ps, ypot_p, d, u, a
         write(99,*)  "  cnorm, phase        = ",cnorm, phase
      end if
      !
      ! Now normalize the spinor
      ! Phase modified for tests on EIMEX (October 2010)
      csp%phase = phase     - angular_momentum_l(csp%orbital%kappa)*pi/two
      if (do_normalize) then
         csp%P(1:csp%mtp) = csp%P(1:csp%mtp) * cnorm
         csp%Q(1:csp%mtp) = csp%Q(1:csp%mtp) * cnorm
         if (debug_normalize_csp) then
            write(99,*) " "
            write(99,*) "Normalized continuum spinor:"
            write(99,*) "   i)        R(i)              P(i)              Q(i)"
            do  i = 1,10,1
               write(99,1) i, r_grasp2k(i), csp%P(i), csp%Q(i)
            end do
            do  i = 11,csp%mtp,100
               write(99,1) i, r_grasp2k(i), csp%P(i), csp%Q(i)
            end do
          1 format(i5,")  ",7(1pe16.9,2x))
            do  i = 100,csp%mtp
               if (csp%P(i-1) - csp%P(i-2) > zero  .and. &
                   csp%P(i)   - csp%P(i-1) < zero  .and. &
                   r_grasp2k(i) > 20.0_dp	   .and. &
                   r_grasp2k(i) < 30.0_dp)	   then
                  !! print *, "Maxima: i, r_grasp2k(i), csp%P(i) = ", &
                  !!                   i, r_grasp2k(i), csp%P(i)
               end if
            end do
         end if
      end if
      !
      if (present(norm_factor)) norm_factor = cnorm
      !
   end subroutine cowf_normalize_wkb
   !
   !
   subroutine cowf_open_csp()
   !--------------------------------------------------------------------
   ! Opens a .csp  Continuum SPinor output File on stream 25. 
   !--------------------------------------------------------------------
      !
      integer :: ierr
      character(len=3)   :: month
      character(len=8)   :: cdate
      character(len=10)  :: ctime
      character(len=256) :: cowf_csp_file
      !
    1 print *, "Enter a file name for the  cowf.csp  file:"
      read *,  cowf_csp_file
      call file_open(25,cowf_csp_file,"formatted  ","new",ierr)
      !
      if (ierr /= 0) goto 1
      !
      call date_and_time(date=cdate,time=ctime)
      month = get_month(cdate(5:6))
      write(25,*) "COWF continuum spinor output file opened at "//   &
                  ctime(1:2)//":"//ctime(3:4)//":"//ctime(5:6)//     &
		  " on "//month//" "//cdate(7:8)//" "//cdate(1:4)//"."
      write(25,*)     " "
      !
   end subroutine cowf_open_csp
   !
   !
   subroutine cowf_schmidt_orthogonalization(csp)
   !--------------------------------------------------------------------
   ! Orthogonalizes the continuum spinor csp with respect to all bound-
   ! state orbital using Schmidt' procedure. 
   ! 
   ! Calls: .
   !--------------------------------------------------------------------
      !
      type(orbital_function), intent(inout)  :: csp
      !
      integer       :: ib, mtp
      real(kind=dp) :: overlap
      real(kind=dp), dimension(n_grasp2k+10) :: ta
      !
      ta(1)  = zero
      do  ib = 1,wave_bound%number_of_rwf
         if (wave_bound%rwf(ib)%orbital%kappa == csp%orbital%kappa) then
            !
            ! Compute overlap; determine the maximum tabulation point 
            mtp = min(csp%mtp, wave_bound%rwf(ib)%mtp)
            !
            ! Tabulate the integrand as required for subroutine ; the
            ! value at the first tabulation point is arbitrary
            ta(2:mtp) = (  wave_bound%rwf(ib)%P(2:mtp) * csp%P(2:mtp)      &
                         + wave_bound%rwf(ib)%Q(2:mtp) * csp%Q(2:mtp)  )   &
                         * rp_grasp2k(2:mtp)
            !
            ! Perform the quadrature
            overlap = quad_grasp2k(ta,mtp)
            !
            ! Schmidt orthogonalization
            csp%pz       = csp%pz - overlap * wave_bound%rwf(ib)%pz
            csp%P(2:mtp) = csp%P(2:mtp) - overlap * wave_bound%rwf(ib)%P(2:mtp)
            csp%Q(2:mtp) = csp%Q(2:mtp) - overlap * wave_bound%rwf(ib)%Q(2:mtp)
         end if
      end do
      !
   end subroutine cowf_schmidt_orthogonalization
   !
   !
   subroutine cowf_set_channels()
   !--------------------------------------------------------------------
   ! Defines all necessary "channels" as selected at input into an 
   ! appropriate type(cowf_channel) data structure.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer :: i, j
      !
      number_of_channels = number_of_selected_channels
      allocate( channel(1:number_of_channels) )
      !
      do  i = 1,number_of_selected_channels
         channel(i)%asf      = channel_selection(i)%asf
         channel(i)%level_No = channel_selection(i)%level_No
         channel(i)%kappa_c  = channel_selection(i)%kappa_c 
         channel(i)%totalJ   = channel_selection(i)%totalJ
         channel(i)%parity   = channel_selection(i)%parity
	 channel(i)%number_of_energies = number_of_selected_energies
	 allocate( channel(i)%energy(1:number_of_selected_energies) )
         do  j = 1,number_of_selected_energies
            channel(i)%energy(j) = energy_selection(j)
         end do
      end do
      !
   end subroutine cowf_set_channels
   !
   !
   subroutine cowf_set_debug()
   !--------------------------------------------------------------------
   ! Open if appropriate a .dbg file and sets flags for debug print out 
   ! from the COWF program.
   !
   ! Calls: file_open(), print_conversion_factors(), 
   !        print_physical_constants().
   !--------------------------------------------------------------------
      !
      integer :: ierr
      logical :: yes
      character(len=256) :: cowf_dbg_file
      !
      print *, "Generate debug printout ?"
      yes = get_yes_stream()
      if (yes) then
         !
         ! The .dbg file is formatted; open it on unit 99
    1    print *, "Enter a file name for the COWF DeBuG Printout "// &
                  "cowf.dbg  file:"
         read *,  cowf_dbg_file
         call file_open(99,cowf_dbg_file,"formatted  ","new",ierr)
         !
         if (ierr /= 0) goto 1
         !
         ! Set options for general printout 
         print *, "Print out the physical constants used ?"
         yes = get_yes_stream()
         if (yes) then
            call print_physical_constants(99)
            debug_physical_constants = .true.
         end if
         print *, "Print out the conversion factors used ?"
         yes = get_yes_stream()
         if (yes) then
            call print_conversion_factors(99)
            debug_conversion_factors = .true.
         end if
         !
         ! Set debug options for radial modules
         print *, "Printout from radial modules ?"
         yes = get_yes_stream()
         if (yes) then
            print *, "Printout from radgrd_grasp2k() ?"
            yes = get_yes_stream()
            if (yes) debug_radgrd_grasp2k = .true.
            print *, "Printout from load_rwf_file_grasp2k() ?"
            yes = get_yes_stream()
            if (yes) debug_load_rwf_grasp2k = .true.
            print *, "Printout from schmidt_orthogonalize_grasp2k() ?"
            yes = get_yes_stream()
            if (yes) debug_schmidt_grasp2k = .true.
            print *, "Printout from interpolate_rwf_grasp2k() ?"
            yes = get_yes_stream()
            if (yes) debug_interpolate_grasp2k = .true.
         end if
      end if
      !
   end subroutine cowf_set_debug
   !
   !
   subroutine cowf_set_direct_potential(csp)
   !--------------------------------------------------------------------
   ! Tabulates the 'direct' potential function Y(r) for the continuum 
   ! spinor csp (Eq (14) in  I P Grant, B J McKenzie, P H Norrington, 
   ! D F Mayers, and N C Pyper, Computer  Phys  Commun  21 (1980) 211).
   ! The function is tabulated in the array cowf_y_pot(:).
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(orbital_function), intent(in) :: csp
      !
      integer :: i, ib, id, kappa, mtp
      real(kind=dp), dimension(1:n_grasp2k+10) :: yk
      !
      kappa = csp%orbital%kappa
      !
      ! Initialize array yp with the nuclear potential piece
      cowf_y_pot(1:n_grasp2k) =                                               &
                              -r_grasp2k(1:n_grasp2k)*cowf_nuc_pot(1:n_grasp2k)
      !
      do  i = 1,number_of_yk
         if (cowf_yk(i)%a == asf_cont%csf_set%nwshells  .and. &
	     cowf_yk(i)%b <= wave_bound%number_of_rwf   .and. &
	     cowf_yk(i)%d == 0)  then
	    ib = cowf_yk(i)%b
	    call yz_k_grasp2k(cowf_yk(i)%k,wave_bound%rwf(ib),       &
	                                   wave_bound%rwf(ib),yk,mtp)
            cowf_y_pot(1:mtp) = cowf_y_pot(1:mtp) - cowf_yk(i)%y * yk(1:mtp)
	    !
	    !!x print *, "b: i, cowf_yk(i)%y, yk(10) = ",i, cowf_yk(i)%y, yk(10)
            cowf_y_pot(mtp+1:n_grasp2k) = cowf_y_pot(mtp+1:n_grasp2k) - &
	                                  cowf_yk(i)%y * yk(mtp)
            print *, "k,b,y,mtp = ",cowf_yk(i)%k,ib,cowf_yk(i)%y,mtp
         else if (cowf_yk(i)%a == asf_cont%csf_set%nwshells  .and. &
	          cowf_yk(i)%b <= wave_bound%number_of_rwf   .and. &
	          cowf_yk(i)%d <= wave_bound%number_of_rwf)  then
	    ib = cowf_yk(i)%b;   id = cowf_yk(i)%d
	    call yz_k_grasp2k(cowf_yk(i)%k,wave_bound%rwf(ib),       &
	                                   wave_bound%rwf(id),yk,mtp)
            cowf_y_pot(1:mtp) = cowf_y_pot(1:mtp) - cowf_yk(i)%y * yk(1:mtp)
	    !
	    !!x print *, "b: i, cowf_yk(i)%y, yk(10) = ",i, cowf_yk(i)%y, yk(10)
            cowf_y_pot(mtp+1:n_grasp2k) = cowf_y_pot(mtp+1:n_grasp2k) - &
	                                  cowf_yk(i)%y * yk(mtp)
            !!x print *, "k,b,d,y,mtp = ",cowf_yk(i)%k,ib,id,cowf_yk(i)%y,mtp
	 else if (abs(csp%orbital%kappa) > 9) then
	    ! do nothing 
	 else if (rabs_use_stop) then 
	    print *, "cowf_yk(i)%a = ",cowf_yk(i)%a
	    print *, "cowf_yk(i)%b = ",cowf_yk(i)%b
	    print *, "cowf_yk(i)%d = ",cowf_yk(i)%d
	    stop "cowf_set_direct_potential(): program stop A."
	 end if
      end do
      !
   end subroutine cowf_set_direct_potential
   !
   !
   subroutine cowf_set_drs_coefficients(asf_No,csf_set,index)
   !--------------------------------------------------------------------
   ! Calculates the d_rs coefficients and the generalized occupation
   ! numbers q_a for the current configuration scheme as needed for the
   ! computation of individual potential contributions. 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in)               :: asf_No
      type(csf_basis), intent(in)       :: csf_set
      integer, dimension(:), intent(in) :: index
      !
      integer       :: i, r, s
      real(kind=dp) :: sum 
      !
      ! Allocate storage for the matrix of the d_rs coefficients and the
      ! generalized occupation numbers
      if (allocated(cowf_drs))  deallocate( cowf_drs )
      allocate( cowf_drs(1:csf_set%nocsf,1:csf_set%nocsf) )
      if (allocated(cowf_gen_occupation))  deallocate( cowf_gen_occupation )
      allocate( cowf_gen_occupation(1:csf_set%nwshells) )
      !
      ! Calculate the d_rs coefficients in the OL mode
      do  r = 1,csf_set%nocsf
         do  s = 1,csf_set%nocsf
	    cowf_drs(r,s) = asf_bound%asf(asf_No)%eigenvector(index(r)) * &
	                    asf_bound%asf(asf_No)%eigenvector(index(s))
         end do
      end do
      !
      ! Calculate the generalized occupation numbers
      do  i = 1,csf_set%nwshells
         cowf_gen_occupation(i) = zero
         do  r = 1,csf_set%nocsf
	    cowf_gen_occupation(i) = cowf_gen_occupation(i) + &
	                       cowf_drs(r,r) * csf_set%csf(r)%occupation(i)
         end do
      end do
      !
      ! Test the set-up of d_rs coefficients
      sum = zero
      do  r = 1,csf_set%nocsf
         sum = sum + cowf_drs(r,r)
      end do
      !
      if (rabs_use_stop       .and.                    &
         (sum > 1.000001_dp   .or.   sum < 0.999999_dp)) then
         print *, " nocsf, sum = ",csf_set%nocsf,sum
         stop "cowf_set_drs_coefficients(): program stop A."
      end if
      !
      if (debug_drs_coefficients) then
         write(99,*) " "
         write(99,*) "Coefficients d_rs: "
         write(99,*) "    r)    s = 1,...,csf_set%nocsf  "
         write(99,*) " "
         do  r = 1,csf_set%nocsf
	    write(99,2) r,(cowf_drs(r,s),s=1,csf_set%nocsf)
         end do
       2 format(3x,i3,")   ",50(f8.5,3x))
         write(99,*) " "
         write(99,*) "Generalized occupation numbers: "
         write(99,*) "    i)  Orbital   Occupation  "
         do  i = 1,csf_set%nwshells
	    write(99,3) i,orbital_name(csf_set%subshell(i)%n,      &
	                               csf_set%subshell(i)%kappa), &
	                               cowf_gen_occupation(i)
         end do
       3 format(3x,i3,")    ",a4,3x,e12.5)
      end if
      !
   end subroutine cowf_set_drs_coefficients
   !
   !
   subroutine cowf_set_exchange_potential(csp)
   !--------------------------------------------------------------------
   ! Calculates the 'exchange' potentials X_p and X_q for the continuum 
   ! spinor csp. 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(orbital_function), intent(in) :: csp
      !
      integer :: i, ib, ic, id, mtp, mtpc
      real(kind=dp), dimension(1:n_grasp2k+10) :: yk
      !
      ! Initialize arrays X^p and X^q to zero
      cowf_xp_pot(1:n_grasp2k) = zero;   cowf_xq_pot(1:n_grasp2k) = zero
      !
      do  i = 1,number_of_xk
         if (cowf_xk(i)%a == asf_cont%csf_set%nwshells  .and. &
	     cowf_xk(i)%b <= wave_bound%number_of_rwf   .and. &
	     cowf_xk(i)%c == 0 .and. cowf_xk(i)%d == 0  .and. &
	     cowf_xk(i)%a /= cowf_xk(i)%b)  then
	    ib = cowf_xk(i)%b
	    call yz_k_grasp2k(cowf_xk(i)%k,csp,wave_bound%rwf(ib),yk,mtp)
            cowf_xp_pot(1:mtp) = cowf_xp_pot(1:mtp) + cowf_xk(i)%x *  &
	                         yk(1:mtp) * wave_bound%rwf(ib)%Q(1:mtp)
            cowf_xq_pot(1:mtp) = cowf_xq_pot(1:mtp) + cowf_xk(i)%x *  &
	                         yk(1:mtp) * wave_bound%rwf(ib)%P(1:mtp)
         else if (cowf_xk(i)%a == asf_cont%csf_set%nwshells  .and. &
	          cowf_xk(i)%b <= wave_bound%number_of_rwf   .and. &
	          cowf_xk(i)%c /= 0 .and. cowf_xk(i)%d /= 0  .and. &
	          cowf_xk(i)%c <= wave_bound%number_of_rwf   .and. &
	          cowf_xk(i)%d <= wave_bound%number_of_rwf   .and. &
		  cowf_xk(i)%a /= cowf_xk(i)%c)  then
	    ib = cowf_xk(i)%b;   ic = cowf_xk(i)%c;   id = cowf_xk(i)%d
	    call yz_k_grasp2k(cowf_xk(i)%k,wave_bound%rwf(ib),     &
	                                   wave_bound%rwf(id),yk,mtp)
            mtpc = min( mtp, wave_bound%rwf(ic)%mtp )
            cowf_xp_pot(1:mtpc) = cowf_xp_pot(1:mtpc) + cowf_xk(i)%x *   &
	                         yk(1:mtpc) * wave_bound%rwf(ic)%Q(1:mtpc)
            cowf_xq_pot(1:mtpc) = cowf_xq_pot(1:mtpc) + cowf_xk(i)%x *   &
	                         yk(1:mtpc) * wave_bound%rwf(ic)%P(1:mtpc)
         else if (cowf_xk(i)%a == asf_cont%csf_set%nwshells  .and. &
	          cowf_xk(i)%b <= wave_bound%number_of_rwf   .and. &
	          cowf_xk(i)%c /= 0 .and. cowf_xk(i)%d /= 0  .and. &
	          cowf_xk(i)%c <= wave_bound%number_of_rwf   .and. &
	          cowf_xk(i)%d == asf_cont%csf_set%nwshells  .and. &
		  cowf_xk(i)%a /= cowf_xk(i)%c)  then
	    ib = cowf_xk(i)%b;   ic = cowf_xk(i)%c;   id = cowf_xk(i)%d
	    call yz_k_grasp2k(cowf_xk(i)%k,wave_bound%rwf(ib),csp,yk,mtp)
            mtpc = min( mtp, wave_bound%rwf(ic)%mtp )
            cowf_xp_pot(1:mtpc) = cowf_xp_pot(1:mtpc) + cowf_xk(i)%x *    &
	                          yk(1:mtpc) * wave_bound%rwf(ic)%Q(1:mtpc)
            cowf_xq_pot(1:mtpc) = cowf_xq_pot(1:mtpc) + cowf_xk(i)%x *    &
	                          yk(1:mtpc) * wave_bound%rwf(ic)%P(1:mtpc)
	 else if (rabs_use_stop) then
            print *, "k,a,b,c,d,x = ",cowf_xk(i)%k,cowf_xk(i)%a,cowf_xk(i)%b, &
                                      cowf_xk(i)%c,cowf_xk(i)%d,cowf_xk(i)%x
	    stop "cowf_set_exchange_potential(): program stop A."
	 end if
      end do
      cowf_xp_pot(:) =  cowf_xp_pot(:) / c 
      cowf_xq_pot(:) =  cowf_xq_pot(:) / c 
      !
      ! print *, "cowf_xp_pot1:5) = ",cowf_xp_pot(1:5)
      !
   end subroutine cowf_set_exchange_potential
   !
   !
   subroutine cowf_set_lagrange_multipliers(subshell,csp)
   !--------------------------------------------------------------------
   ! Calculates the 'off-diagonal' Lagrange multipliers for the continuum 
   ! spinor csp. It cycles over all bound-state orbitals and, thus,
   ! b /= a is fulfilled automatically. Also, q_bar(a) = 1 is assumed in
   ! this procedure.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(nkappa), intent(in)           :: subshell
      type(orbital_function), intent(in) :: csp
      !
      integer       :: i, ib, mtp
      real(kind=dp) :: epsilon_ab
      real(kind=dp), dimension(n_grasp2k+10) :: ta
      !
      number_of_epsilon_ab = 0
      ta(1)  = zero
      do  ib = 1,wave_bound%number_of_rwf
         if (wave_bound%rwf(ib)%orbital%kappa /= subshell%kappa) cycle
	 mtp = wave_bound%rwf(ib)%mtp
	 ta(2:mtp)  = wave_bound%rwf(ib)%P(2:mtp) * cowf_xq_pot(2:mtp) &
	            - wave_bound%rwf(ib)%Q(2:mtp) * cowf_xp_pot(2:mtp)
	 ta(2:mtp)  = ta(2:mtp) * rp_grasp2k(2:mtp) / r_grasp2k(2:mtp)
	 epsilon_ab = c * quad_grasp2k(ta,mtp)
	 !
	 mtp = min( wave_bound%rwf(ib)%mtp, csp%mtp )
	 ta(2:mtp)  = cowf_y_pot(2:mtp) + r_grasp2k(2:mtp)*cowf_nuc_pot(2:mtp)
	 ta(2:mtp)  = ta(2:mtp) * rp_grasp2k(2:mtp) / r_grasp2k(2:mtp) *   &
	              (wave_bound%rwf(ib)%P(2:mtp)*csp%P(2:mtp) +          &
	 	       wave_bound%rwf(ib)%Q(2:mtp)*csp%Q(2:mtp))
	 epsilon_ab = epsilon_ab + quad_grasp2k(ta,mtp)  &                  
	            - I_ab_grasp2k(csp,wave_bound%rwf(ib))
	 !
	 number_of_epsilon_ab = number_of_epsilon_ab + 1
	 if (rabs_use_stop   .and.                                         &
            (asf_cont%csf_set%subshell(asf_cont%csf_set%nwshells)%n     /= &
	     csp%orbital%n   .or.                                          &
	     asf_cont%csf_set%subshell(asf_cont%csf_set%nwshells)%kappa /= &
	     csp%orbital%kappa)) then
	    stop "cowf_set_lagrange_multipliers(): program stop A."
	 end if
         !!x print *, "a, b, epsilon = ",asf_cont%csf_set%nwshells,ib,epsilon_ab
	 cowf_eab(number_of_epsilon_ab)%a       = asf_cont%csf_set%nwshells
	 cowf_eab(number_of_epsilon_ab)%b       = ib
	 cowf_eab(number_of_epsilon_ab)%epsilon = epsilon_ab
      end do
      !
      if (debug_lagrange_epsilon_ab) then
         write(99,*) " "
         write(99,*) "Lagrange multipliers epsilon_ab: "
         write(99,*) "    i)   ia   ib      epsilon_ab "
         do  i = 1,number_of_epsilon_ab
            write(99,2) i, cowf_eab(i)%a, cowf_eab(i)%b, cowf_eab(i)%epsilon  
         end do
       2 format(3x,i3,")  ",2(i3,2x),3x,e12.5)
      end if
      !
   end subroutine cowf_set_lagrange_multipliers
   !
   !
   subroutine cowf_set_yk_coefficients(subshell,csf_set)
   !--------------------------------------------------------------------
   ! Calculates the y_k(ab) and y_k(abad) coefficients which are required 
   ! for the set-up of the direct part of the MC Dirac-Fock potentials 
   ! of the current configuration scheme.
   ! The calculation of these coefficients follows the paper of
   ! Dyall et al. CPC 55, p. 425 (1989). Only non-zero coefficients
   ! will be stored. 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(nkappa), intent(in)    :: subshell
      type(csf_basis), intent(in) :: csf_set
      !
      integer       :: i, ia, ib, id, k, k_max, number_of_yk_max, r, s
      real(kind=dp) :: sum, wa, ykab, ykabad
      !
      number_of_yk_max = 0;   number_of_yk = 0
    1 number_of_yk_max = number_of_yk_max + 1000
      if (allocated(cowf_yk))  deallocate( cowf_yk ) 
      allocate( cowf_yk(1:number_of_yk_max) ) 
      !
      ia = 0
      do  ib = 1,csf_set%nwshells
         if (subshell%n     == csf_set%subshell(ib)%n .and. &
	     subshell%kappa == csf_set%subshell(ib)%kappa) then
	    ia = ib
	    goto 3
	 end if
      end do
      stop "cowf_set_yk_coefficients(): program stop A."
      !
      ! Go directly to 3 and use the V^k_Coulomb coefficients only
      !-----------------------------------------------------------
      !
      ! Calculate y_k (ab) coefficients for the direct potential
    2 do  ib = 1,csf_set%nwshells
         if (ia == ib) then;   k_max = 0;   else;   
	    k_max = angular_momentum_j(subshell%kappa) - 1
	 end if
         do  k = 0,k_max,2
	    sum = zero
	    do  r = 1,csf_set%nocsf
	       wa = cowf_fk_r(k,r,ia,ib,csf_set)
	       sum = sum + cowf_drs(r,r) * wa
	    end do
	    if (ia == ib) then
	       ykab = two * sum / cowf_gen_occupation(ia)
	    else
	       ykab = sum / cowf_gen_occupation(ia)
	    end if
	    !
	    if (abs(ykab) > eps10) then
	       number_of_yk = number_of_yk + 1
	       if (number_of_yk > number_of_yk_max) then
	          print *, "cowf_set_yk_coefficients():"//&
		           " number_of_yk_max need to be increased."
	          goto 1
	       end if
	       cowf_yk(number_of_yk)%k = k;    cowf_yk(number_of_yk)%y = ykab
	       cowf_yk(number_of_yk)%a = ia;   cowf_yk(number_of_yk)%b = ib
	       cowf_yk(number_of_yk)%d = 0
	    end if
	 end do
      end do	
      !
      ! Calculate y_k (abad) coefficients for the direct potential
    3 do  ib = 1,csf_set%nwshells
         do  id = ib,csf_set%nwshells
            !! do  k = 0,8,2
            do  k = 0,0,2
	       sum = zero
	       do  r = 1,csf_set%nocsf
	             wa = cowf_Vk_rs(r,r,k,subshell,csf_set%subshell(ib),  &
                                           subshell,csf_set%subshell(id))
	             sum = sum + cowf_drs(r,r) * wa
	       end do
	       ykabad = sum / cowf_gen_occupation(ia)
	       !
	       if (abs(ykabad) > eps10) then
	          number_of_yk = number_of_yk + 1
	          if (number_of_yk > number_of_yk_max) then
	             print *, "cowf_set_yk_coefficients():"//&
		              " number_of_yk_max need to be increased."
	             goto 1
	          end if
	          cowf_yk(number_of_yk)%k = k;  cowf_yk(number_of_yk)%y = ykabad
	          cowf_yk(number_of_yk)%a = ia; cowf_yk(number_of_yk)%b = ib
	          cowf_yk(number_of_yk)%d = id
		  !!x print *, "a: number_of_yk, cowf_yk(number_of_yk)%d = ", &
		  !!x              number_of_yk, cowf_yk(number_of_yk)%d
	       end if
	    end do
	 end do
      end do	
      !
      return
      !
      ! Calculate y_k (abad) coefficients for the direct potential
      do  ib = 1,csf_set%nwshells
         do  id = ib,csf_set%nwshells
            do  k = 0,8,2
	       sum = zero
	       do  r = 1,csf_set%nocsf
	          do  s = 1,csf_set%nocsf
	             wa = cowf_Vk_rs(r,s,k,subshell,csf_set%subshell(ib),  &
                                           subshell,csf_set%subshell(id))
	             sum = sum + cowf_drs(r,s) * wa
		  end do
	       end do
	       ykabad = sum / cowf_gen_occupation(ia)
	       !
	       if (abs(ykabad) > eps10) then
	          number_of_yk = number_of_yk + 1
	          if (number_of_yk > number_of_yk_max) then
	             print *, "cowf_set_yk_coefficients():"//&
		              " number_of_yk_max need to be increased."
	             goto 1
	          end if
	          cowf_yk(number_of_yk)%k = k;  cowf_yk(number_of_yk)%y = ykabad
	          cowf_yk(number_of_yk)%a = ia; cowf_yk(number_of_yk)%b = ib
	          cowf_yk(number_of_yk)%d = id
		  print *, "b: number_of_yk, cowf_yk(number_of_yk)%d = ", &
		               number_of_yk, cowf_yk(number_of_yk)%d
	       end if
	    end do
	 end do
      end do	
      !
      if (debug_yk_coefficients) then
         write(99,*) " "
         write(99,*) "Coefficients y_k(ab): "
         write(99,*) "    i)    k   ia   ib   id      y_k(ab/ad) "
         do  i = 1,number_of_yk
            write(99,4) i,cowf_yk(i)%k, cowf_yk(i)%a, cowf_yk(i)%b, &
	                cowf_yk(i)%d, cowf_yk(i)%y
         end do
       4 format(3x,i3,")  ",4(i3,2x),3x,e12.5)
      end if
      !
   end subroutine cowf_set_yk_coefficients
   !
   !
   subroutine cowf_set_xk_coefficients(subshell,csf_set)
   !--------------------------------------------------------------------
   ! Calculates the x_k(ab) and x_k(abcd) coefficients which are required 
   ! for the set-up of the exchange part of the MC Dirac-Fock potentials 
   ! of the current configuration scheme.
   ! The calculation of these coefficients follows the paper of
   ! Dyall et al. CPC 55, p. 425 (1989). Only non-zero coefficients
   ! will be stored. The x_k(ab) coefficients are stored with
   ! c = d = 0. 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(nkappa), intent(in)    :: subshell
      type(csf_basis), intent(in) :: csf_set
      !
      integer       :: i, ia, ib, ic, id, ja, jb, k, kl, ku, &
                       number_of_xk_max, r, s
      real(kind=dp) :: sum, wa, xkab, xkabcd
      !
      number_of_xk_max = 0;   number_of_xk = 0
      !
    1 number_of_xk_max = number_of_xk_max + 1000
      if (allocated(cowf_xk))  deallocate( cowf_xk ) 
      allocate( cowf_xk(1:number_of_xk_max) ) 
      !
      ia = 0
      do  ib = 1,csf_set%nwshells
         if (subshell == csf_set%subshell(ib)) then
	    ia = ib
	    goto 3
	 end if
      end do
      !
      ! Go directly to 3 and use the V^k_Coulomb coefficients only
      !-----------------------------------------------------------
      !
      ! Calculate x_k (ab) coefficients for the exchange potential
      do  ib = 1,csf_set%nwshells
         ja = angular_momentum_j(subshell%kappa)
         jb = angular_momentum_j(csf_set%subshell(ib)%kappa)
         if (subshell%kappa * csf_set%subshell(ib)%kappa > 0) then
	    kl = abs(ja-jb) / 2
         else
	    kl = abs(ja-jb) / 2 + 1
	 end if
         if (mod((ja+jb)/2-kl,2) == 0) then
	    ku = (ja+jb) / 2
	 else
	    ku = (ja+jb) / 2 - 1
	 end if
	 ! 
         do  k = kl,ku,2
	    if (ib == ia) cycle
            sum = zero
            do  r = 1,csf_set%nocsf
               wa  = cowf_gk_r(k,r,ia,ib,csf_set)
               sum = sum + cowf_drs(r,r) * wa
            end do
            xkab = sum / cowf_gen_occupation(ia)
	    !
	    if (abs(xkab) > eps10) then
	       number_of_xk = number_of_xk + 1
	       if (number_of_xk > number_of_xk_max) then
	          print *, "cowf_set_xk_coefficients(): number_of_xk_max"//&
		           " need to be increased."
	          goto 1
	       end if
	       cowf_xk(number_of_xk)%k = k 
	       cowf_xk(number_of_xk)%a = ia;  cowf_xk(number_of_xk)%b = ib
	       cowf_xk(number_of_xk)%c = 0;   cowf_xk(number_of_xk)%d = 0  
	       cowf_xk(number_of_xk)%x = xkab
	    end if
	 end do
      end do
      !
      ! Now add x_k (abcd) coefficients for the exchange potential
    3 do  k = 0,6
         do  ib = 1,csf_set%nwshells
            do  ic = 1,csf_set%nwshells
	       if (ic == ia) cycle
               do  id = ib,csf_set%nwshells
                  sum = zero
                  do  s = 1,csf_set%nocsf
                     do  r = 1,s
                        if (abs(cowf_drs(r,s)) > 1.0e-2_dp) then
                           sum = sum + cowf_drs(r,s) * &
			      cowf_Vk_rs(r,s,k,subshell,csf_set%subshell(ib), &
                                   csf_set%subshell(ic),csf_set%subshell(id))
                        end if
                     end do
                  end do
                  xkabcd = sum / cowf_gen_occupation(ia)
	          !
	          if (abs(xkabcd) > eps10) then
	             number_of_xk = number_of_xk + 1
	             if (number_of_xk > number_of_xk_max) then
	                print *, "cowf_set_xk_coefficients():"//&
		                 " number_of_xk_max need to be increased."
	                goto 1
	             end if
	             cowf_xk(number_of_xk)%k = k 
	             cowf_xk(number_of_xk)%a = ia 
		     cowf_xk(number_of_xk)%b = ib
	             cowf_xk(number_of_xk)%c = ic   
		     cowf_xk(number_of_xk)%d = id  
		     cowf_xk(number_of_xk)%x = xkabcd
	          end if
	          !
               end do
            end do
         end do
      end do
      !
      if (debug_xk_coefficients) then
         write(99,*) " "
         write(99,*) "Coefficients x_k(abcd): "
         write(99,*) "    i)    k   ia   ib   ic   id      x_k(ab/cd)  "
         do  i = 1,number_of_xk
            write(99,4) i,cowf_xk(i)%k, cowf_xk(i)%a, cowf_xk(i)%b, &
	                cowf_xk(i)%c, cowf_xk(i)%d, cowf_xk(i)%x
         end do
       4 format(3x,i3,")  ",5(i3,2x),3x,e12.5)
      end if
      !
   end subroutine cowf_set_xk_coefficients
   !
   !
   subroutine cowf_start_integration(csp,P0,Q0,P,Q)
   !--------------------------------------------------------------------
   ! This subroutine sets up P(1:6) and Q(1:6) to enable the start of
   ! the integration.
   !
   ! Arguments: 
   !
   !   csp  : (Input) continuum spinor.
   !   P0   : (Input) Slope parameter 
   !   P    : (Output) P(1:6) are tabulated by this program 
   !   Q0   : (Output) First term in the series expansion of Q
   !   Q    : (Output) Q(1:6) are tabulated by this program 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(orbital_function), intent(in)         :: csp
      real(kind=dp), intent(in)                  :: P0
      real(kind=dp), intent(out)                 :: Q0
      real(kind=dp), dimension(:), intent(inout) :: P, Q
      !
      integer, parameter :: mxiter = 96
      integer            :: i, ib, j, niter
      real(kind=dp)      :: coefij, difmaw, difmax, enefac, factor,  &
                            pi, pzero, qi, p1, q1, ri, rj, rpi, &
                            rirpi, rsep1, rseq1, sump, sumq, ypirpi
      real(kind=dp), dimension(1:6) :: rdp, rdq, rsep, rseq, spest, sqest
      !
      ! Determine p(1), q(1): these store  r**(-gamma)*(p(1), q(1));
      ! set up  rsep  and  rseq , the inhomogeneous terms
      select case(nuclear_model)
      case("point")
         p(1) = p0
	 if (csp%orbital%kappa < 0) then
            q(1) = -p0 * nuclear_charge / &
	           (c * (csp%gamma - csp%orbital%kappa))
         else
            q(1) = p0 * c * (csp%gamma + csp%orbital%kappa) / nuclear_charge
         end if
      case default
         if (csp%orbital%kappa < 0) then
            p(1) = p0;   q(1) = zero
         else
            p(1) = zero   
            q(1) = p0 * c * (csp%gamma + csp%orbital%kappa) / nuclear_charge
         end if
      end select
      !
      if (cowf_start_homogeneous) then
         rsep = zero;   rseq = zero
      else
         rsep1 = zero;   rseq1 = zero
         do  i = 1,number_of_epsilon_ab
	    ib    = cowf_eab(i)%b
	    pzero = wave_bound%rwf(ib)%pz
            select case(nuclear_model)
            case("point")
               p(1) = pzero
	       if (csp%orbital%kappa < 0) then
                  q(1) = -pzero * nuclear_charge/                           &
		         (c*(csp%gamma -csp%orbital%kappa))
               else  
                  q(1) = pzero *c *                                         &
		         (csp%gamma + csp%orbital%kappa) / nuclear_charge
               end if
               sump =  nuclear_charge /c * q1
               sumq = -nuclear_charge /c * p1
            case default
               if (csp%orbital%kappa < 0) then
                  p(1) = pzero;   q(1) = zero
               else
                  p(1) = zero   
                  q(1) = pzero * c * (csp%gamma + csp%orbital%kappa) /         &
		         nuclear_charge
               end if
	       sump = zero;   sumq = zero
            end select
            rsep1 = rsep1 + cowf_eab(i)%epsilon *                          &
	            (sump + (one-csp%gamma - csp%orbital%kappa) * p1)
            rseq1 = rseq1 + cowf_eab(i)%epsilon *                          &
	            (sump + (one-csp%gamma + csp%orbital%kappa) * p1)
	 end do
	 !
         factor  = rp_grasp2k(1)
         rsep(1) = factor * rsep1
         rseq(1) = factor * rseq1
         do  i = 2,6
            factor = -rp_grasp2k(i) * r_grasp2k(i)**(-csp%gamma)
            rsep(i) = factor * cowf_xp_pot(i)
            rseq(i) = factor * cowf_xq_pot(i)
         end do
      end if
      q0 = q(1)
      !
      ! Set up  rdp  and  rdq
      !!x csq    = c*c
      !!x twocsq = csq + csq
      enefac = two*c*c - csp%energy
      do  i = 1,6
         ri     = r_grasp2k(i);   rpi = rp_grasp2k(i)
         rirpi  = ri * rpi
         ypirpi = cowf_y_pot(i) * rpi
         rdp(i) = -(enefac*rirpi     + ypirpi) / c
         rdq(i) = -(csp%energy*rirpi - ypirpi) / c
      end do
      !
      ! Determine  P(2:6) , Q(2:6);
      ! initilizations for the iterations
      niter  = 0
      p1     = p(1)
      q1     = q(1)
      difmaw = max(abs(p1), abs(q1))
      p(2:6) = p1;   q(2:6) = q1
      !
      ! This is the largest factor by which any result will be multiplied
      factor = r_grasp2k(6)**csp%gamma
      !
      ! Now iterate
    7 niter   = niter + 1
      difmax  = zero
      do j = 2,6
         sump = zero;   sumq = zero
         do  i = 1,6
            coefij = cnc6c_grasp2k(i,j)
            rpi    = rp_grasp2k(i)
            pi     = p(i)
            qi     = q(i)
            sump   = sump + coefij * ((one-csp%gamma-csp%orbital%kappa)* &
	                              rpi*pi-rdp(i)*qi+rsep(i))
            sumq   = sumq + coefij * ((one-csp%gamma+csp%orbital%kappa)* &
	                              rpi*qi-rdq(i)*pi+rseq(i))
         end do
         rj   = r_grasp2k(j)
         sump = sump / rj;   sumq = sumq / rj
         spest(j) = sump;    sqest(j) = sumq
         difmax   = max(difmax, abs(sump-p(j)) )
         difmax   = max(difmax, abs(sumq-q(j)) )
      end do
      if (difmax < difmaw) then
         p(2:6) = spest(2:6);   q(2:6) = sqest(2:6)
         difmaw   = difmax
         difmax   = difmax * factor
         if (difmax > accy_grasp2k) then
            if (niter < mxiter) then
               goto 7
            else
               write (*,8) orbital_name(csp%orbital%n,csp%orbital%kappa),  &
	                   difmax,niter,accy_grasp2k
            end if
         end if
      else
         difmax = difmax * factor
         if (difmax > accy_grasp2k) then
            write (*,8) orbital_name(csp%orbital%n,csp%orbital%kappa),  &
	                difmax,niter,accy_grasp2k
          8 format("cowf_start(): ",a," subshell: accuracy ",1pe7.1,    &
                   " attained after ",i3," iterations; this fails the", &
                   " accuracy criterion ",e7.1,".")
         end if
      end if
      !
      ! All done
      ! this is always true in GRASP92
      p(1) = zero;   q(1) = zero
      do  i = 2,6
         factor = r_grasp2k(i)**csp%gamma
         p(i)   = factor * p(i)
         q(i)   = factor * q(i)
      end do
      !
      if (debug_start_integration) then
         write(99,*) " "
         write(99,*) "First points of P and Q;   P0, Q0 =",P0, Q0
         write(99,*) "   i)        R(i)              P(i)              Q(i)"
         do  i = 1,6
            write(99,9) i, r_grasp2k(i), P(i), Q(i)
          9 format(i5,")  ",3(1pe16.9,2x))
         end do
      end if
      !
   end subroutine cowf_start_integration
   !
   !
   subroutine cowf_total_phase_shift(csp_unp,csp_full,n0,phase)
   !--------------------------------------------------------------------
   ! Calculates 
   ! 
   ! Calls: .
   !--------------------------------------------------------------------
      !
      type(orbital_function), intent(in) :: csp_unp, csp_full
      integer, intent(in)                :: n0
      real(kind=dp), intent(out)         :: phase
      !
      integer       :: i, ii, n_max
      real(kind=dp) :: P_i, P_ip1, P_max, X_max, phase2
      !
      ! Find the 'last positive zero' of the unperturbed solution, i.e.
      ! the point r where P(r+dr) ~ A sin(dr)
      P_ip1 = csp_unp%P(n0)
      do  i = n0-900,n0-1800,-1
         P_i   = csp_unp%P(i)
         if (P_i < zero   .and.   P_ip1 > zero) goto 1
         P_ip1 = P_i
      end do
      !
    1 continue
      !
      if (abs(P_i) > abs(P_ip1)) then
         ii = i + 1
      else
         ii = i
      end if
      !
      ! Determine the constant Amplitude A
      P_max = abs(csp_full%P(n0))
      do  i = ii,ii+900
         if (P_max < abs(csp_full%P(i))) P_max = abs(csp_full%P(i))
      end do
      !
      X_max = zero
      do  i = ii,ii+500
         if (X_max < csp_unp%P(i)) then
            X_max = csp_unp%P(i);   n_max = i
         end if
      end do
      !
      phase  = asin(csp_full%P(ii)/P_max)
      phase2 = asin(csp_full%P(n_max) / P_max)
      print *, "phase, phase2 = ",phase, phase2
      if (phase > zero   .and.   phase2 > zero) then
         phase2 = phase2 + two * (pi*half - phase2) - half*pi
      else if (phase < zero  .and.  phase2 < zero  .and. phase2 > phase) then
         phase  = phase  - two * abs(-pi*half - phase)
         phase2 = phase2 - pi*half
      else if (phase < zero  .and.  phase2 < zero  .and. phase2 < phase) then
         phase  = phase  - two * abs(-pi*half - phase)
         phase2 = phase2 - half*pi
      else if (phase > zero  .and.  phase2 < zero) then
         phase  = phase  + two * (pi*half - phase)
         phase2 = phase2 - two * abs(-pi*half - phase2) + three*half*pi
      else
         phase2 = phase2 - half*pi
      end if
      !   
      if (phase  < zero) phase  = phase  + two*pi
      if (phase2 < zero) phase2 = phase2 + two*pi
      print *, "A: n0, ii, csp_full%P(ii), P_max, phase = ", &
                   n0, ii, csp_full%P(ii), P_max, phase
      print *, "B: n_max, csp_full%P(n_max), X_max, P_max, phase2 = ", &
                   n_max, csp_full%P(n_max), X_max, P_max, phase2
      !
   end subroutine cowf_total_phase_shift
   !
   !
   function cowf_Vk_rs(r,s,k,a,b,c,d)                      result(value)
   !--------------------------------------------------------------------
   ! This routine looks for a Coulomb angular coefficient as specified by 
   ! the indices r,s of the CSF and the four orbital quantum numbers. 
   ! Before a call to this routine should can be made, one has to ensure
   ! that the ANCO component was called for the current configuration
   ! scheme.           
   ! 
   ! Calls: angular_momentum_j(), angular_momentum_l(), CL_reduced_me(),
   !        triangle().
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: k, r, s
      type(nkappa)        :: a, b, c, d
      real(kind=dp)       :: value
      !
      integer       :: i, ja, jb, jc, jd, la, lb, lc, ld, t
      real(kind=dp) :: xc
      !
      value = zero
      !
      ja = angular_momentum_j(a%kappa);   la = angular_momentum_l(a%kappa)
      jb = angular_momentum_j(b%kappa);   lb = angular_momentum_l(b%kappa)
      jc = angular_momentum_j(c%kappa);   lc = angular_momentum_l(c%kappa)
      jd = angular_momentum_j(d%kappa);   ld = angular_momentum_l(d%kappa)
      !
      if (triangle(ja+1,jc+1,k+k+1) * triangle(jb+1,jd+1,k+k+1) == 0  .or. &
          mod(la+lc+k,2) == 1   .or.   mod(lb+ld+k,2) == 1) then
         return
      end if 
      !
      xc = CL_reduced_me(a%kappa,k,c%kappa) *  CL_reduced_me(b%kappa,k,d%kappa)
      if (mod(k,2) == 1) then
         xc = - xc
      end if
      !
      if (abs(xc) .lt. eps10) then
         return
      end if
      !
      do  i = 1,number_of_pair_list
         if (anco_pair_list(i)%r /= r   .or.   anco_pair_list(i)%s /= s) cycle
         do  t = 1,anco_pair_list(i)%no_V_coeff
	    !!x print *, "c: i, t, xc, anco_pair_list(i)%V_coeff(t)%V = ", &
	    !!x              i, t, xc, anco_pair_list(i)%V_coeff(t)%V
            if (anco_pair_list(i)%V_coeff(t)%nu /= k)          cycle
            if       (a == anco_pair_list(i)%V_coeff(t)%a   .and. &
                      b == anco_pair_list(i)%V_coeff(t)%b   .and. &
                      c == anco_pair_list(i)%V_coeff(t)%c   .and. &
                      d == anco_pair_list(i)%V_coeff(t)%d)  then
               value = anco_pair_list(i)%V_coeff(t)%V * xc
               return
            else if  (c == anco_pair_list(i)%V_coeff(t)%a   .and. &
                      d == anco_pair_list(i)%V_coeff(t)%b   .and. &
                      a == anco_pair_list(i)%V_coeff(t)%c   .and. &
                      b == anco_pair_list(i)%V_coeff(t)%d)  then
               value = anco_pair_list(i)%V_coeff(t)%V * xc
               return
            else if  (b == anco_pair_list(i)%V_coeff(t)%a   .and. &
                      a == anco_pair_list(i)%V_coeff(t)%b   .and. &
                      d == anco_pair_list(i)%V_coeff(t)%c   .and. &
                      c == anco_pair_list(i)%V_coeff(t)%d)  then
               value = anco_pair_list(i)%V_coeff(t)%V * xc
               return
            else if  (d == anco_pair_list(i)%V_coeff(t)%a   .and. &
                      c == anco_pair_list(i)%V_coeff(t)%b   .and. &
                      b == anco_pair_list(i)%V_coeff(t)%c   .and. &
                      a == anco_pair_list(i)%V_coeff(t)%d)  then
               value = anco_pair_list(i)%V_coeff(t)%V * xc
               return
            end if
         end do
      end do
      !
   end function cowf_Vk_rs
   !
end module rabs_cowf
