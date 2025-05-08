module rabs_auger
!
!***** July 2014 *****
!-----------------------------------------------------------------------
! This module contains the procedures which are specific to the AUGER
! program. This program supports the calculation of nonradiative rates,
! relative intensities, and angular distribution parameters.
! The continuum orbitals which are required for these computations are
! generated automatically by a call to the COWF component.
! At the present, the program supports nonorthogonality between the
! orbitals only via the form of the radial functions. In the evaluation
! of the many-electron matrix elements, by contrast, orthogonality is
! 'assumed' and  standard Racah algebra techniques used. See the component
! ANCO for calculating the angular integrals for one- and two-particle
! matrix elements.
! In this module below, there are several procedures also for the file
! handling and the intermediate storage to accelerate some of the
! computations.
!-----------------------------------------------------------------------
   !
   use rabs_anco
   use rabs_constant
   use rabs_csl
   use rabs_cowf
   use rabs_file_handling
   use rabs_grasp2k
   use rabs_input_dialog
   use rabs_nonorthonormal
   use rabs_print
   use rabs_xl
   use omp_lib
   implicit none
   !
   private :: auger_angular_parameter
                 ! Returns the alpha_k angular distribution parameter for a
		 ! given transition.
   private :: auger_anisotropy_Ak_kabachnik
                 ! Calculates the generalized anisotropy coefficient due to
                 ! Kabachnik (1994).
   private :: auger_anisotropy_Bbar_kabachnik
                 ! Calculates the generalized coefficient B_bar_k1,k,k0 ()
                 ! due to Kabachnik (1994).
   public  :: auger_calculate_amplitudes
                 ! Calculates for all selected transitions (in turn) the
                 ! required continuum spinors and Auger amplitudes.
   private :: auger_channel_amplitude
                 ! Calculates the Auger amplitude of a channel from the 'pure'
		 ! Auger matrix and the corresponding mixing coefficients.
   public  :: auger_collect_input
                 ! Collects and proceeds all input for the calculation of
		 ! Auger rates, lifetimes, and angular distribution parameters.
   private :: auger_delta_energy
                 ! Returns the energy difference between two levels of the
                 ! same multiplett of a given 'initial' level.
   public  :: auger_initialize
                 ! Set up the selected transitions and initializes some
		 ! arrays as required.
   !!x public  :: auger_initialize_rwf_storage
   !!x               ! Initializes the arrays of type(grasp2k_orbital) for the
   !!x               ! storage of the radial wave functions.
   public  :: auger_print_results
                 ! Writes the transition rates, lifetimes, and others to
                 ! the .sum file.
   private :: auger_print_transitions
                 ! Prints all selected transitions in a neat format before
                 ! the computation starts.
   public  :: auger_print_amplitudes
                 ! Prints the information about all (selected) transitions
                 ! and the amplitudes to a (.trn) transition amplitude file.
   private :: auger_pure_matrix
                 ! Calculates the 'pure' Auger matrix for the given
		 ! configuration scheme.
   private :: auger_rho_after_dipole
   private :: auger_rho_density
   private :: auger_rho_dipole_excitation
                 ! Calculates the density matrix of the coherently-excited
                 ! initial states of an Auger cascade by linearly polarized
                 ! light along the z-axis.
   private :: auger_set_overlaps
                 ! Initializes the array of overlap integrals for the
                 ! calculation of relaxed transition rates.
   private :: auger_set_transitions
                 ! Determines which Auger transitions need to be calculated
                 ! and initializes storage for an appropriate data type for
                 ! them.
   private :: auger_spin_parameter
                 ! Returns the eta_k spin polarization coefficient for a
		 ! given transition.
   private :: auger_transition_properties
                 ! Calculates all selected Auger properties of transition i
		 ! from the amplitudes of the individual channels.
   private :: auger_width
                 ! Returns the width of a given 'initial' level.
   !
   ! Define some global data of the AUGER program; most of these data are
   ! read in during the interactive control and may overwrite existing
   ! default values
   !
   ! Define a 'configuration scheme' to be built up during execution
   type(grasp2k_orbital), public :: wave_continuum
   type(orbital_function)        :: auger_csp
   !
   ! Define an internal structure type(auger_transition) which stores all
   ! necessary information for an radiative transition line
   type :: auger_channel
      integer          :: kappa
      real(kind=dp)    :: phase, amplitude_re
      complex(kind=dp) :: amplitude
   end type auger_channel
   !
   type :: auger_transition
      integer          :: asfi, asff, level_i, level_f, totalJ_i, totalJ_f
      integer          :: number_of_channels
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: energy, probability
      real(kind=dp)    :: alpha_2, alpha_4, eta_2, eta_4
      real(kind=dp), dimension(:), pointer       :: bessel
      type(auger_channel), dimension(:), pointer :: channel
   end type auger_transition
   !
   type(auger_transition), dimension(:), allocatable :: transition
   !
   type :: auger_matrix
      integer :: no_f, no_i
      integer, dimension(:), pointer :: ndx_f, ndx_i
      real(kind=dp), dimension(:,:), pointer :: matrix
   end type auger_matrix
   !
   type(auger_matrix) :: auger
   !
   integer          :: number_of_transitions                 = 0
   integer, public  :: auger_maximal_kappa                   = 20
   !
   integer, private                            :: number_radial_integrals
   real(kind=dp),     dimension(1:500), private :: radial_value
   character(len=25), dimension(1:500), private :: radial_string
   !
   ! Define global logical flags for the control of the AUGER program; the
   ! default values for these flags may be overwritten interactively during
   ! input time
   logical, public :: auger_add_1                     = .false.,  &
                      auger_add_2                     = .false.,  &
                      auger_add_3                     = .false.,  &
                      auger_add_4                     = .false.,  &
                      auger_add_5                     = .false.,  &
                      auger_add_6                     = .false.,  &
                      auger_add_7                     = .false.,  &
                      auger_add_8                     = .false.,  &
                      auger_add_9                     = .false.,  &
                      auger_apply_exp_energies        = .false.,  &
                      auger_calc_angular_parameter    = .true.,   &
                      auger_calc_spin_polarization    = .false.,  &
                      auger_include_breit             = .false.,  &
                      auger_nonorthogonal_eval        = .false.,  &
		      auger_print_csf_scheme          = .false.,  &
		      auger_print_each_transition     = .true.,   &
		      auger_print_main_csf_me         = .false.,  &
		      auger_print_radial_integrals    = .false.,  &
		      auger_print_rates_in_hartree    = .true.,   &
                      auger_print_selected_trans      = .true.,   &
                      auger_print_trn_file            = .false.,  &
                      auger_sort_transition_energy    = .true.,   &
                      auger_use_formatted_mix_file    = .true.,   &
                      auger_use_formatted_rwf_file    = .true.
   !
   ! Energy unit for the output of all energies
   real(kind=dp)    :: auger_rate_factor = zero,                  &
                       auger_print_cut   = 0.001_dp,               &
                       auger_energy_shift   = zero,               &
                       auger_maximal_energy = 30.0_dp,            &
                       !! auger_mean_transition   = 30.0_dp,                    &
                       auger_lowest_transition    = zero
   character(len=7) :: auger_rate_unit
   !
   ! Define storage for the overlap integrals: (kappa,pqn_f,pqn_i)
   real(kind=dp), dimension(:,:,:), allocatable :: auger_overlap
   !
contains
   !
   function auger_angular_parameter(k,trans)               result(alpha)
   !--------------------------------------------------------------------
   ! Returns the alpha_k angular distribution parameter for the transition
   ! trans. It uses the explicit expression in terms of the partial wave
   ! transition amplitudes as given by J Tulkki, H Aksela, and
   ! N M Kabachnik, Phys. Rev. A 48, 1277 (1993), Eqs. A26-27.
   !
   ! Calls: angular_momentum_j(), angular_momentum_l().
   !--------------------------------------------------------------------
      !
      integer, intent(in)    :: k
      type(auger_transition) :: trans
      real(kind=dp)          :: alpha
      !
      integer          :: iphase, ja, jb, la, lb, ma, mb
      real(kind=dp)    :: norm, sum, wa
      !
      ! Calculate the normalization and general phase factor
      norm   = trans%probability / (two * pi)
      iphase = trans%totalJ_i + trans%totalJ_f - 1 + 32
      if (mod(iphase,4) == 2) then
         norm = -norm
      else if (rabs_use_stop       .and.                    &
              (mod(iphase,4) == 1  .or.  mod(iphase,4) == 3)) then
         stop "auger_spin_parameter(): program stop A."
      end if
      !
      sum = zero
      !
      ! Now cycle about the different Auger amplitudes including proper
      ! phase factors and the products of partial-wave expansion coefficients;
      ! Here, b denotes the prime quantum numbers in eq. A27.
      do  ma = 1,trans%number_of_channels
         do  mb = 1,trans%number_of_channels
	    ja = angular_momentum_j(trans%channel(ma)%kappa)
	    jb = angular_momentum_j(trans%channel(mb)%kappa)
	    la = angular_momentum_l(trans%channel(ma)%kappa)
	    lb = angular_momentum_l(trans%channel(mb)%kappa)
	    !
	    wa = (la+la+one)*(lb+lb+one)*(ja+one)*(jb+one)*(trans%totalJ_i+one)
	    wa = sqrt(wa)
	    wa = wa * Clebsch_Gordan(la+la,0,lb+lb,0,k+k,0)
	    wa = wa * wigner_6j_symbol(trans%totalJ_i,trans%totalJ_i,         &
	                               k+k,ja,jb,trans%totalJ_f)              &
	            * wigner_6j_symbol(ja,jb,k+k,lb+lb,la+la,1)               &
	     * trans%channel(ma)%amplitude * conjg(trans%channel(mb)%amplitude)
	    sum = sum + wa
            !!x print *, "sum, norm = ",sum, norm, sum/norm
	 end do
      end do
      !
      alpha = sum / norm
      !
   end function auger_angular_parameter
   !
   !
   function auger_angular_parameter_old(k,trans)               result(alpha)
   !--------------------------------------------------------------------
   ! Returns the alpha_k angular distribution parameter for the transition
   ! trans. It uses the explicit expression in terms of the partial wave
   ! transition amplitudes as given by J Tulkki, H Aksela, and
   ! N M Kabachnik, Phys. Rev. A 48, 1277 (1993), Eqs. A26-27.
   !
   ! Calls: angular_momentum_j(), angular_momentum_l().
   !--------------------------------------------------------------------
      !
      integer, intent(in)    :: k
      type(auger_transition) :: trans
      real(kind=dp)          :: alpha
      !
      integer       :: iphase, ja, jb, la, lb, ma, mb
      real(kind=dp) :: norm, sum, wa
      !
      ! Calculate the normalization and general phase factor
      norm   = trans%probability / (two * pi)
      iphase = trans%totalJ_i + trans%totalJ_f - 1 + 32
      if (mod(iphase,4) == 2) then
         norm = -norm
      else if (rabs_use_stop       .and.                    &
              (mod(iphase,4) == 1  .or.  mod(iphase,4) == 3)) then
         stop "auger_spin_parameter(): program stop A."
      end if
      !
      sum = zero
      !
      ! Now cycle about the different Auger amplitudes including proper
      ! phase factors and the products of partial-wave expansion coefficients;
      ! Here, b denotes the prime quantum numbers in eq. A27.
      do  ma = 1,trans%number_of_channels
         do  mb = 1,trans%number_of_channels
	    ja = angular_momentum_j(trans%channel(ma)%kappa)
	    jb = angular_momentum_j(trans%channel(mb)%kappa)
	    la = angular_momentum_l(trans%channel(ma)%kappa)
	    lb = angular_momentum_l(trans%channel(mb)%kappa)
	    !
	    wa = (la+la+one)*(lb+lb+one)*(ja+one)*(jb+one) &
                *(trans%totalJ_i+one)*(k+k+1)
	    wa = sqrt(wa)
	    iphase = la - lb + 32
            if (mod(iphase,4) == 2) then
                wa = -wa
            end if
	    !
	    wa = wa * wigner_6j_symbol(trans%totalJ_i,trans%totalJ_i,         &
	                               k+k,jb,ja,trans%totalJ_f)              &
	            * wigner_6j_symbol(ja,jb,k+k,lb+lb,la+la,1)               &
		    * cos(trans%channel(ma)%phase - trans%channel(mb)%phase)  &
		    * trans%channel(ma)%amplitude * trans%channel(mb)%amplitude
	    !
	    wa = wa * wigner_3j_symbol(lb+lb,la+la,k+k,0,0,0)
	    !
	    sum = sum + wa
            !!x print *, "sum, norm = ",sum, norm, sum/norm
	 end do
      end do
      !
      alpha = sum / norm
      !
   end function auger_angular_parameter_old
   !
   !
   function auger_anisotropy_Ak_kabachnik(k,J1_level,J1p_level,J2_level) &
                                                                result(Ak)
   !--------------------------------------------------------------------
   ! Calculates the density matrix of the intermediate states from the
   ! statistical tensor of the initial (resonant) state and the amplitudes
   ! of the first-step Auger decay [cf. Ueda et al., JPB (1999), Eq. (6)].
   ! The dummy parameter J1_level, ... refer to the level numbers of the
   ! corresponding transition array
   !
   ! Calls: angular_momentum_j(), angular_momentum_l(),
   !        auger_return_amplitude().
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: k,J1_level,J1p_level,J2_level
      !
      integer          :: i, j, jp, J1, J1p, J2, kappa, kappa_p, l , lp
      complex(kind=dp) :: Ak, wa
      !
      Ak = cmplx(zero,zero)
      !
      ! Determine the total angular momenta of the initial and final states
      J1 = -99;   J1p = -99;   J2 = -99
      do  i = 1,number_of_transitions
         if (transition(i)%level_i == J1_level) then
            J1  = transition(i)%totalJ_i
         end if
         if (transition(i)%level_i == J1p_level) then
            J1p = transition(i)%totalJ_i
         end if
         if (transition(i)%level_f == J2_level) then
            J2  = transition(i)%totalJ_f
         end if
      end do
      !
      if (rabs_use_stop                .and.           &
         (J1 == -99  .or.  J1p == -99  .or.  J2 == -99)) then
         print *, "J1, J1p, J2 = ",J1, J1p, J2
         stop "auger_anisotropy_Ak_kabachnik(): program stop A"
      end if
      !
      ! Sum over kappa and kappa_p
      do  kappa = -12,12
         if (kappa == 0) cycle
         j  = angular_momentum_j(kappa)
         l  = angular_momentum_l(kappa)
         do  kappa_p = -12,12
            if (kappa_p == 0) cycle
            jp  = angular_momentum_j(kappa_p)
            lp  = angular_momentum_l(kappa_p)
            wa  = sqrt( (l+l+one)*(lp+lp+one)*(j+one)*(jp+one)          * &
                               (J1+one)*(J1p+one) )
            wa  = wa * Clebsch_Gordan(l+l,0,lp+lp,0,k+k,0)              * &
                  wigner_6j_symbol(l+l,j,1,jp,lp+lp,k+k)                * &
                  wigner_6j_symbol(J1,j,J2,jp,J1p,k+k)                  * &
                  auger_return_amplitude(J1_level,J2_level,kappa)       * &
                  conjg(auger_return_amplitude(J1p_level,J2_level,kappa_p))
            Ak  = Ak + wa
         end do
      end do
      !
      ! Add phase information
      if (mod(J1+J2-1+32,4) == 0) then
      else if (mod(J1+J2-1+32,4) == 2) then
         Ak = -Ak
      else if (rabs_use_stop) then
         stop "auger_anisotropy_Ak_kabachnik(): program stop B."
      end if
      !
   end function auger_anisotropy_Ak_kabachnik
   !
   !
   function auger_anisotropy_Bbar_kabachnik(k1,k,k0,J0_level,J1_level, &
                                            J1p_level)         result(Bbar)
   !--------------------------------------------------------------------
   ! Calculates the generalized coefficient Bar{B}_k1,k,k0 (J0,J1,J1').
   ! The dummy parameters J0_level, ... refer to the level numbers of the
   ! corresponding transition array
   !
   ! Calls: angular_momentum_j(), angular_momentum_l(),
   !        auger_return_amplitude().
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: k, k0, k1, J0_level, J1_level,J1p_level
      !
      integer          :: i, j, jp, J0, J1, J1p, kappa, kappa_p, l , lp
      complex(kind=dp) :: Bbar, wa
      !
      Bbar = cmplx(zero,zero)
      !
      ! Determine the total angular momenta of the initial and final states
      J0 = -99;   J1 = -99;   J1p = -99
      do  i = 1,number_of_transitions
         if (transition(i)%level_i == J0_level) then
            J0  = transition(i)%totalJ_i
         end if
         if (transition(i)%level_f == J1_level) then
            J1  = transition(i)%totalJ_f
         end if
         if (transition(i)%level_f == J1p_level) then
            J1p = transition(i)%totalJ_f
         end if
      end do
      !
      if (rabs_use_stop                .and.           &
         (J0 == -99  .or.  J1 == -99  .or.  J1p == -99)) then
         print *, "J0, J1, J1p = ",J0, J1, J1p
         stop "auger_anisotropy_Bbar_kabachnik(): program stop A"
      end if
      !
      ! Sum over kappa and kappa_p
      do  kappa = -12,12
         if (kappa == 0) cycle
         j  = angular_momentum_j(kappa)
         l  = angular_momentum_l(kappa)
         do  kappa_p = -12,12
            if (kappa_p == 0) cycle
            jp  = angular_momentum_j(kappa_p)
            lp  = angular_momentum_l(kappa_p)
            wa  = sqrt( (l+l+one)*(lp+lp+one)*(j+one)*(jp+one) )
            wa  = wa * Clebsch_Gordan(l+l,0,lp+lp,0,k1+k1,0)            * &
                  wigner_6j_symbol(l+l,j,1,jp,lp+lp,k1)                 * &
                  wigner_9j_symbol(J1,j,J0,J1p,jp,J0,k+k,k1+k1,k0+k0)   * &
                  auger_return_amplitude(J0_level,J1_level,kappa)       * &
                  conjg(auger_return_amplitude(J0_level,J1p_level,kappa_p))
            if (mod(jp+1,4) == 0) then
               Bbar = Bbar + wa
            else if (mod(jp+1,4) == 2) then
               Bbar = Bbar - wa
            else if (rabs_use_stop) then
               stop "auger_anisotropy_Bbar_kabachnik(): program stop B."
            end if
         end do
      end do
      !
   end function auger_anisotropy_Bbar_kabachnik
   !
   !
   subroutine auger_calculate_amplitudes()
   !--------------------------------------------------------------------
   ! Calculates for all transitions in turn the required continuum
   ! spinors and Auger amplitudes.
   !
   ! Calls: add_csf_to_basis(), anco_calculate_csf_matrix(),
   ! auger_channel_amplitude(), auger_pure_matrix(),
   ! auger_transition_properties() cowf_iterate_csp(),
   ! cowf_set_xk_coefficients(), cowf_set_yk_coefficients(),
   ! print_configuration_scheme()
   ! set_configuration_scheme().
   !--------------------------------------------------------------------
      integer       :: i, j, n, nw, nocsf
      real(kind=dp) :: energy
      type(nkappa)  :: subshell
      integer, dimension(:), allocatable :: ndx

      ! Variables that need to be private to each thread
      integer, dimension(:), allocatable :: local_ndx_f, local_ndx_i
      real(kind=dp), dimension(:,:), allocatable :: local_matrix

      n = asf_final%csf_set%nocsf + asf_initial%csf_set%nocsf

      ! Allocate global arrays outside parallel region
      allocate( auger_csp%P(1:n_grasp2k), auger_csp%Q(1:n_grasp2k) )
      allocate( cowf_csp%P(1:10), cowf_csp%Q(1:10) )
      allocate( ndx(1:n) )

      ! Parallelize the main loop over transitions
      !$OMP PARALLEL PRIVATE(i, j, energy, subshell, local_ndx_f, local_ndx_i, &
      !$OMP                  local_matrix, nw, nocsf) SHARED(ndx)
      !$OMP DO SCHEDULE(DYNAMIC)
      do  i = 1,number_of_transitions
         if (transition(i)%energy < zero) then
            !$OMP CRITICAL
            transition(i)%probability = zero
            transition(i)%alpha_2     = zero
            transition(i)%alpha_4     = zero
            transition(i)%eta_2       = zero
            transition(i)%eta_4       = zero
            !$OMP END CRITICAL
            cycle
         end if

         do  j = 1,transition(i)%number_of_channels
            energy = transition(i)%energy
            call set_configuration_scheme(asf_final%csf_set,asf_cont%csf_set,&
                     -1,transition(i)%channel(j)%kappa,                      &
                     transition(i)%totalJ_f,transition(i)%parity_f,          &
                     transition(i)%totalJ_i,transition(i)%parity_i,          &
                     append=.false.,index=ndx)

            ! Allocate thread-private arrays
            auger%no_f = asf_cont%csf_set%nocsf
            allocate( local_ndx_f(auger%no_f) )
            local_ndx_f(1:auger%no_f) = ndx(1:auger%no_f)

            nw = asf_cont%csf_set%nwshells
            if (rabs_use_stop  .and. nw /= asf_final%csf_set%nwshells + 1) then
               stop "auger_calculate_amplitudes(): program stop A."
            end if

            nocsf = asf_cont%csf_set%nocsf
            call anco_calculate_csf_matrix(asf_cont%csf_set,1,nocsf,1,nocsf)
            call cowf_set_drs_coefficients(transition(i)%asff,             &
                                           asf_cont%csf_set,ndx)
            subshell = nkappa(-1,transition(i)%channel(j)%kappa)
            call cowf_set_yk_coefficients(subshell,asf_cont%csf_set)
            call cowf_set_xk_coefficients(subshell,asf_cont%csf_set)

            ! Set continuum orbital calculation parameters
            cowf_start_homogeneous         = .true.
            cowf_phaseshift_wkb            = .true.
            cowf_phaseshift_zero_potential = .false.
            cowf_phaseshift_coulomb        = .false.
            cowf_norm_wkb                  = .true.
            call cowf_iterate_csp(energy,subshell)

            !$OMP CRITICAL
            auger_csp = cowf_csp
            transition(i)%channel(j)%phase = auger_csp%phase
            !$OMP END CRITICAL

            ! Define extended configuration scheme
            call add_csf_to_basis(asf_initial%csf_set,asf_cont%csf_set,      &
                    transition(i)%totalJ_i,transition(i)%parity_i,index=ndx)
            if (auger_print_csf_scheme) then
               call print_configuration_scheme(6,asf_cont%csf_set)
            end if

            auger%no_i = asf_cont%csf_set%nocsf - auger%no_f
            allocate( local_ndx_i(auger%no_i) )
            local_ndx_i(1:auger%no_i) = ndx(1+auger%no_f:asf_cont%csf_set%nocsf)
            allocate( local_matrix(1:auger%no_f,1:auger%no_i) )

            ! Calculate Auger matrix
            call auger_pure_matrix(asf_cont%csf_set,i)

            ! Move matrix results to shared storage with critical section
            !$OMP CRITICAL
            allocate( auger%ndx_f(auger%no_f) )
            allocate( auger%ndx_i(auger%no_i) )
            allocate( auger%matrix(1:auger%no_f,1:auger%no_i) )
            auger%ndx_f = local_ndx_f
            auger%ndx_i = local_ndx_i
            auger%matrix = local_matrix
            call auger_channel_amplitude(i,j)
            deallocate( auger%ndx_f, auger%ndx_i, auger%matrix )
            !$OMP END CRITICAL

            deallocate( local_ndx_f, local_ndx_i, local_matrix )
            call deallocate_csf_basis(asf_cont%csf_set)
         end do

         ! Calculate transition properties
         call auger_transition_properties(transition(i))
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      deallocate( ndx, auger_csp%P, auger_csp%Q)
   end subroutine auger_calculate_amplitudes
   !
   !
   subroutine auger_channel_amplitude(i,j)
   !--------------------------------------------------------------------
   ! Calculates the Auger amplitude of channel j of transition i
   ! by summing over the 'pure' Auger matrix using the proper weights of
   ! transition i.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: i, j
      !
      integer       :: asfi, asff, l, r, rr, s, ss
      real(kind=dp) :: phase, value
      !
      if (auger_print_main_csf_me) then
         print *, " "
         print *, "Main contribution from initial- and final-state CSF "// &
                  "(abs(c_i*c_f) > 0.01)"
         print *, "----------------------------------------------------"// &
                  "---------------------"
         print *, " "
         print *, "   I-CSF     F-CSF   kappa    c_i     c_f       c_i*c_f"//&
                  "   c_i*c_f*A_if  "
         print *, "-------------------------------------------------------"//&
                  "-----------------"
      end if
      !
      asfi  = transition(i)%asfi;  asff = transition(i)%asff
      value = zero
      do  r = 1,auger%no_f
         rr = auger%ndx_f(r)
         do  s = 1,auger%no_i
            ss = auger%ndx_i(s)
	    value = value + asf_final%asf(asff)%eigenvector(rr) * &
	            auger%matrix(r,s) * asf_initial%asf(asfi)%eigenvector(ss)
            !
            if (auger_print_main_csf_me) then
               if (abs(asf_final%asf(asff)%eigenvector(rr)*                   &
                       asf_initial%asf(asfi)%eigenvector(ss)) > 0.01_dp  .and.&
                   abs(asf_final%asf(asff)%eigenvector(rr) *                  &
	               auger%matrix(r,s) *                                    &
                   asf_initial%asf(asfi)%eigenvector(ss)) > 0.000001_dp)   then
                  !
                  ! Determine first the radial integrals from the occupation
                  ! of the CSF
                  write(*,1) ss,rr,                                         &
                     orbital_symmetry(transition(i)%channel(j)%kappa),      &
                             asf_initial%asf(asfi)%eigenvector(ss),         &
                             asf_final%asf(asff)%eigenvector(rr),           &
                             asf_final%asf(asff)%eigenvector(rr)*           &
                             asf_initial%asf(asfi)%eigenvector(ss),         &
                             asf_final%asf(asff)%eigenvector(rr)*           &
                             asf_initial%asf(asfi)%eigenvector(ss)*         &
                             auger%matrix(r,s)
                1 format(1x,i7,i10,6x,a2,3x,f6.3,2x,f6.3,5x,f8.5,4x,f9.6)
                end if
            end if
            !
	 end do
      end do
      !
      if (auger_print_main_csf_me) then
         print *, "-------------------------------------------------------"//&
                  "-----------------"
      end if
      !
      l     = angular_momentum_l(transition(i)%channel(j)%kappa)
      phase = transition(i)%channel(j)%phase
      !
      transition(i)%channel(j)%amplitude_re = value
      transition(i)%channel(j)%amplitude    = cmplx(zero,one)**l *           &
                              exp( -cmplx(zero,one)*phase) * cmplx(value,zero)
      !
      print *, "i,j,transition(i)%channel(j)%amplitude = ",  &
                i,j,transition(i)%channel(j)%amplitude
      !
   end subroutine auger_channel_amplitude
   !
   !
   subroutine auger_collect_input()
   !--------------------------------------------------------------------
   ! Collect and proceeds all input for the AUGER program.
   !
   ! Calls: get_yes_stream().
   !--------------------------------------------------------------------
      !
      integer            :: level_i, level_f, score_position, kappa, ierr
      logical            :: yes
      real(kind=dp)      :: maximal_energy, wavelength
      character(len=20 ) :: string
      character(len=256) :: record, auger_trn_file
      !
      !
      ! Open, check, load data from, and close the  .iso  file
      call file_get_isodat_grasp2k()
      !
      ! Determine the units for the printout and further optional input
      call input_energy_unit()
      !
      !
    2 print *, "Enter the maximal energy of the Auger transitions "//&
	       " (in " // trim(energy_unit)		           //&
               ") to built-up the radial grid;"
      read (*, *, err=2) auger_maximal_energy
      !
      if (energy_inverse) then
         auger_maximal_energy = energy_factor * auger_maximal_energy
      else
         auger_maximal_energy = auger_maximal_energy / energy_factor
      end if
      !
      if (auger_maximal_energy < one) goto 2
      !
      ! Determine grid parameters
      call input_grid_parameters("standard")
      !
      wavelength = sqrt( two * pi * pi / auger_maximal_energy )
      hp_grasp2k = wavelength / 30.0_dp
      !! n_grasp2k  = 10.0_dp / hp_grasp2k + 3100
      n_grasp2k  = 10.0_dp / hp_grasp2k + 9100
      n_grasp2k  = 10.0_dp / hp_grasp2k + 19100
      !
      !! n_grasp2k  = 10.0_dp / hp_grasp2k + 9100
      !
      ! Default speed of light
      c = c_vacuum
      !
      ! Now 'overwrite' defaults only if required
      print *, "Modify default set-up and printout of the program ?"
      yes = get_yes_stream()
      if (.not.yes) goto 10
      !
      ! Select individual pairs of transitions
      call input_transition_pairs(number_of_transitions)
      !
      print *, "Include exchange interactions into the generation of the"//&
	       " continuum waves ?"
      yes = get_yes_stream()
      cowf_solve_homogeneous_eqn = .not.yes
      !
      auger_nonorthogonal_eval = .false.
      !
      print *, "Include Breit interactions to the Auger matrix ?"
      auger_include_breit = get_yes_stream()
      !
      print *, "Calculate angular distribution parameters ?"
      auger_calc_angular_parameter = get_yes_stream()
      !
      print *, "Calculate spin polarization parameters ?"
      auger_calc_spin_polarization = get_yes_stream()
      !
      print *, "Sort transitions in ascending order of energy ?"
      auger_sort_transition_energy = get_yes_stream()
      !
      print *, "Read in and apply experimental energies for the calculation"//&
               " of Auger rates and other properties ?"
      auger_apply_exp_energies = get_yes_stream()
      !
      print *, "Auger rates are printed in SI units;"// &
               " use Hartree atomic units instead ?"
      auger_print_rates_in_hartree = get_yes_stream()
      !
      print *, "Print all selected transitions and their energies"// &
               " before the computation starts (this is the default) ?"
      auger_print_selected_trans = get_yes_stream()
      !
      print *, "Print the CSF scheme each time a new one has been built ?"
      auger_print_csf_scheme =  get_yes_stream()
      !
      print *, "Print the results for each individual transition immediatly"//&
               " after its computation ?"
      auger_print_each_transition = get_yes_stream()
      !
      print *, "Print the final results to a (.trn) transition amplitude file ?"
      auger_print_trn_file = get_yes_stream()
      if (auger_print_trn_file) then
    3    print *, "Enter a file name for the  auger.trn  file:"
         read *,  auger_trn_file
         call file_open(27,auger_trn_file,"formatted  ","new",ierr)
         !
         if (ierr /= 0) goto 3
      end if
      !
      print *, "Enter an (overall) shift of the Auger energies which applies"//&
               " to all transitions (in"//trim(energy_unit)//"):"
      print *, " Use  0.  or   <cr>  if no shift need to be applied."
      read (*, *, err=5) auger_energy_shift
      !
      if (energy_inverse) then
         auger_energy_shift = energy_factor * auger_energy_shift
      else
         auger_energy_shift = auger_energy_shift / energy_factor
      end if
      !
    5 print *, "Enter a minimal energy (> =0.) of the free electron"        //&
               " (in"//trim(energy_unit)//"):"
      print *, " All other transitions are neglected from the computations;"
      print *, " use  0.  or   <cr>  if all possible transitions are to"    //&
               " be taken into account."
      read (*, *, err=6) auger_lowest_transition
      !
      if (energy_inverse) then
         auger_lowest_transition = energy_factor * auger_lowest_transition
      else
         auger_lowest_transition = auger_lowest_transition / energy_factor
      end if
      !
    6 print *, "Enter a maximal (-)kappa symmetry up to which continuum"//&
               " spinors are taken into account ?"
      print *, " 2 (up to p-waves), 4(f), 6(i), 8(k), ...;"//&
               " 0 or <cr> to include all possible waves."
      read (*, "(i2)", err=8) kappa
      if (kappa == 0) then
      else if (abs(kappa) < auger_maximal_kappa) then
         auger_maximal_kappa = kappa
      end if
      !
      ! Determine the physical effects specifications
    8 print *, "The physical speed of light in atomic units is",c_vacuum,";"
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
      if (auger_print_rates_in_hartree) then
         auger_rate_factor = one
         auger_rate_unit   = "a.u.   "
      else
         auger_rate_factor = convert_au_to_per_sec
         auger_rate_unit   = "1/s    "
      end if
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
   end subroutine auger_collect_input
   !
   !
   function auger_delta_energy(level_1,level_2,asf_set)  result(delta_e)
   !--------------------------------------------------------------------
   ! Returns the energy difference between two given levels of the same
   ! multiplett.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in)         :: level_1, level_2
      type(asf_basis), intent(in) :: asf_set
      real(kind=dp)               :: delta_e, e1, e2
      !
      integer :: i
      !
      do  i = 1,asf_set%noasf
         if (level_1 == asf_set%asf(i)%level_No) then
            e1 = asf_set%asf(i)%energy
            goto 1
         end if
      end do
      !
      print *, "Level not included in the given set; level_1 = ",level_1
      stop "auger_delta_energy(): program stop A."
      !
    1 do  i = 1,asf_set%noasf
         if (level_2 == asf_set%asf(i)%level_No) then
            e2 = asf_set%asf(i)%energy
            goto 2
         end if
      end do
      !
      print *, "Level not included in the given set; level_2 = ",level_2
      stop "auger_delta_energy(): program stop A."
      !
    2 delta_e = (e1 - e2)
      !
   end function auger_delta_energy
   !
   !
   subroutine auger_initialize()
   !--------------------------------------------------------------------
   ! Initializes the computation of Auger rates, lifetimes, angular
   ! distribution parameters and others.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer :: i, mtp, noasf, nocsf, number_of_rwf
      !
      print *, "Initialize the set-up of the overlap integrals and "// &
               "transitions ..."
      call auger_set_transitions()
      call auger_set_overlaps()
      !
      print *, "   ... initialization complete."
      print *, " "
      !
      ! Print the selected transitions before the computation starts
      if (auger_print_selected_trans) then
         call auger_print_transitions(6)
         call auger_print_transitions(24)
      end if
      !
      ! Make a 'copy' of asf_final into asf_bound and wave_final into
      ! wave_bound (kept in rabs_cowf) to use the COWF component
      noasf = asf_final%noasf;   nocsf = asf_final%csf_set%nocsf
      allocate( asf_bound%asf(1:noasf) )
      do  i = 1,noasf
         allocate( asf_bound%asf(i)%eigenvector(1:nocsf) )
      end do
      asf_bound%noasf          = asf_final%noasf
      asf_bound%average_energy = asf_final%average_energy
      asf_bound%asf            = asf_final%asf
      !
      number_of_rwf            = wave_final%number_of_rwf
      wave_bound%number_of_rwf = wave_final%number_of_rwf
      allocate( wave_bound%rwf(1:number_of_rwf) )
      do  i = 1,number_of_rwf
         mtp = wave_final%rwf(i)%mtp
         print *, "i, mpt = ",i, mtp
         allocate( wave_bound%rwf(i)%P(1:mtp), wave_bound%rwf(i)%Q(1:mtp) )
         wave_bound%rwf(i) = wave_final%rwf(i)
         print *, "i,wave_bound%rwf(i)%mtp = ",i,wave_bound%rwf(i)%mtp
      end do
      !
   end subroutine auger_initialize
   !
   !
   subroutine auger_print_amplitudes(stream)
   !--------------------------------------------------------------------
   ! Prints the information about all (selected) transitions and amplitudes
   ! to a (.trn) transition amplitude file on stream.
   !
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      !
      integer            ::  i, j, channels
      !
      write(stream,*) &
         "Transition data and amplitudes for Auger transitions are printed "
      write(stream,*) &
         "in the format (one line per transition): "
      write(stream,*) " "
      write(stream,*) "   ASF_i, ASF_f, Level_i, Level_f, 2*J_i, 2*J_f, ..  ."
      write(stream,*) "        ... P_i, P_f, energy, kappa, phase, amplitude "
      write(stream,*) &
         "================================================================ "
      write(stream,*) " "
      !
      channels = 0
      do  i = 1,number_of_transitions
         channels = channels + transition(i)%number_of_channels
      end do
      write(stream,*) channels, "= Number_of_channels"
      write(stream,*) " "
      !
      do  i = 1,number_of_transitions
         do  j = 1,transition(i)%number_of_channels
            write(stream,1) transition(i)%asfi,transition(i)%asff,         &
                            transition(i)%level_i,transition(i)%level_f,   &
                            transition(i)%totalJ_i,transition(i)%totalJ_f, &
                            transition(i)%parity_i,transition(i)%parity_f, &
                            transition(i)%energy,                          &
                            transition(i)%channel(j)%kappa,                &
                            transition(i)%channel(j)%phase,                &
                            transition(i)%channel(j)%amplitude
         end do
      end do
      !
    1 format(2i5,3x,2i5,3x,2i3,2x,2a3,1x,e14.7,i4,2x,e14.7,2x,e14.7,1x,e14.7)
      !
   end subroutine auger_print_amplitudes
   !
   !
   subroutine auger_print_results(stream)
   !--------------------------------------------------------------------
   ! Writes the Auger rates, lifetimes, angular distribution parameters,
   ! and further properties  in a neat summary to the .sum file.
   !
   ! Calls:  angular_momentum_string().
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      !
      integer       :: i
      real(kind=dp) :: energy, rate, ratio, total, total_rate, width
      !
      write(stream,1)
    1 format(/ &
         /10x,"==================================================================", &
         /10x,"|  Summary of all Auger Rates, Angular Parameters and Lifetimes  |", &
         /10x,"==================================================================", &
             // )
    2 format( / 1x,85("-"),                                                  &
              / 2x,"LevI-LevF   I- J / Parity -F     ",                      &
                   "Energy (",a,")     Rate (a.u.) ",                        &
                   "  Total Rate (a.u.)  ",                                  &
              / 1x,85("-") )
    3 format( / 1x,85("-"),                                                  &
              / 2x,"LevI-LevF   I- J / Parity -F     ",                      &
                   "Energy (",a,")     Rate (1/s)  ",                        &
                   "  Amplitude ratio(%) ",                                  &
              / 1x,85("-") )
    4 format(   2x,i3,' -',i3,3x,a4,a2,4x,a4,a2,5x,1pe12.5,5x,               &
                   (1pe10.3,6x), 1pe10.3)
    7 format( / 1x,85("-"),                                                  &
              / 2x,"LevI-LevF   I- J / Parity -F     ",                      &
                   "Energy (",a,")       Alpha_2      Alpha_4     ",         &
              / 1x,85("-") )
    8 format(   2x,i3,' -',i3,3x,a4,a2,4x,a4,a2,5x,1pe12.5,5x,               &
                   4(1pe10.3,3x) )
    9 format(   1x,85("-") )
   10 format( / 1x,85("-"),                                                  &
              / 2x,"LevI-LevF   I- J / Parity -F     ",                      &
                   "Energy (",a,")        Eta_2        Eta_4      ",         &
              / 1x,85("-") )
   11 format(   1x,115("-"),                                                 &
              / 2x," LeveL",15x,"Lifetime","           Total rate",          &
               32x,"Width",                                                  &
              /23x,"--------           ----------",8x,                       &
                   "-----------------------------------------------------",  &
              /24x,"Seconds              1/s   ",12x,"Hartrees",             &
               12x,"Kaysers",15x,"eV",                                       &
              / 1x,115("-") )
   12 format(3x,i4,6x,5(1pd20.7))
   13 format(   1x,115("-"))
      !
      ! Print the Auger rates of this calculation
      !
      write(stream,*) "Individual rate and Amplitude ratio :"
      write(stream,*) "----------------------------------"
      if (auger_print_rates_in_hartree) then
         write(stream,2) trim(energy_unit)
      else
         write(stream,3) trim(energy_unit)
      end if
      !
      total_rate = zero
      do  i = 1,number_of_transitions
         total_rate = total_rate + transition(i)%probability
      end do
      !
      rate = zero; ratio = zero; total = zero
      do  i = 1,number_of_transitions
         if (energy_inverse) then
            energy = energy_factor / transition(i)%energy
         else
            energy = energy_factor * transition(i)%energy
         end if
         rate = transition(i)%probability * auger_rate_factor
         !
         !        auger_print_only_gt_1percent                              &
        ratio = rate/(total_rate *  auger_rate_factor)
        total = total + ratio
        if (abs(ratio) > auger_print_cut) then
            write(stream,4) transition(i)%level_i, transition(i)%level_f,   &
                  trim(angular_momentum_string(transition(i)%totalJ_i,4)),  &
                            transition(i)%parity_i,                         &
                  trim(angular_momentum_string(transition(i)%totalJ_f,4)),  &
                            transition(i)%parity_f,                         &
                            energy,rate,ratio
         end if
      end do
      write(stream,9)
      !
      write(stream,*) "selected total ratio = ", total
      write(stream,*) "Total rate (a.u.)   = ",total_rate
      if (auger_rate_unit /= "a.u.   ") then
         write(stream,*) "Total rate (",trim(auger_rate_unit),") = ", &
                         total_rate * auger_rate_factor
      end if
      write(stream,*) " "
      !
      !
      ! Print Auger angular parameters
      !
      if (auger_calc_angular_parameter) then
         write(stream,*) " "
         write(stream,*) "Auger angular parameters :"
         write(stream,*) "--------------------------"
         write(stream,7) trim(energy_unit)
         do  i = 1,number_of_transitions
            if (energy_inverse) then
               energy = energy_factor / transition(i)%energy
            else
               energy = energy_factor * transition(i)%energy
            end if
            !
            write(stream,8) transition(i)%level_i, transition(i)%level_f,   &
                  trim(angular_momentum_string(transition(i)%totalJ_i,4)),  &
                            transition(i)%parity_i,                         &
                  trim(angular_momentum_string(transition(i)%totalJ_f,4)),  &
                            transition(i)%parity_f,                         &
                            energy,                                         &
                            transition(i)%alpha_2, transition(i)%alpha_4
         end do
         write(stream,9)
      end if
      !
      !
      ! Print Auger spin-polarization parameters
      !
      if (auger_calc_spin_polarization) then
         write(stream,*) " "
         write(stream,*) "Spin-polarization parameters :"
         write(stream,*) "------------------------------"
         write(stream,10) trim(energy_unit)
         do  i = 1,number_of_transitions
            if (energy_inverse) then
               energy = energy_factor / transition(i)%energy
            else
               energy = energy_factor * transition(i)%energy
            end if
            !
            write(stream,8) transition(i)%level_i, transition(i)%level_f,   &
                  trim(angular_momentum_string(transition(i)%totalJ_i,4)),  &
                            transition(i)%parity_i,                         &
                  trim(angular_momentum_string(transition(i)%totalJ_f,4)),  &
                            transition(i)%parity_f,                         &
                            energy,                                         &
                            transition(i)%eta_2,   transition(i)%eta_4
         end do
         write(stream,9)
      end if
      !
      !
      ! Print Auger widths and lifetimes
      !
      write(stream,*) " "
      write(stream,*)  "Auger lifetimes, total rates and widths :"
      write(stream,*)  "-----------------------------------------"
      write(stream,*) " "
      write(stream,11)
      do  i = 1,asf_initial%noasf
         width = auger_width(asf_initial%asf(i)%level_No)
         if (abs(width) > eps10) then
            write(stream,12) i,one / (width * convert_au_to_per_sec),&
	                     width * convert_au_to_per_sec,          &
                             width, width * convert_au_to_kaysers,   &
                             width * convert_au_to_ev
         end if
      end do
      write(stream,13)
      !
   end subroutine auger_print_results
   !
   !
   subroutine auger_print_transitions(stream)
   !--------------------------------------------------------------------
   ! Prints a neat table of all selected transitions on stream before
   ! the actual computation starts; only the quantum numbers of the atomic
   ! states and the transition energies are displayed.
   !
   ! Calls: angular_momentum_string(), get_kappa_channels().
   !--------------------------------------------------------------------
      !
      integer, intent(in)    :: stream
      integer		     :: i, j, m, nchannel, mchannel
      real(kind=dp)	     :: energy
      integer, dimension(20) :: kappa_channel
      !
      write(stream,1) number_of_transitions,trim(energy_unit)
      do  i = 1,number_of_transitions
	 if (energy_inverse) then
            energy = energy_factor / transition(i)%energy
         else
            energy = energy_factor * transition(i)%energy
         end if
	 call get_kappa_channels(                                          &
	          transition(i)%totalJ_i,transition(i)%parity_i,           &
	          transition(i)%totalJ_f,transition(i)%parity_f,           &
		  kappa_channel,nchannel)
         mchannel = 0
         do  m = 1,nchannel
           if (abs(kappa_channel(m)) <= auger_maximal_kappa) then
              mchannel = mchannel + 1
           end if
         end do
         write(stream,2) transition(i)%level_i,transition(i)%level_f,      &
                  trim(angular_momentum_string(transition(i)%totalJ_i,4)), &
                  transition(i)%parity_i,                                  &
                  trim(angular_momentum_string(transition(i)%totalJ_f,4)), &
                  transition(i)%parity_f,energy,                           &
		  (orbital_symmetry(kappa_channel(j)),j=1,mchannel)
      end do
      write(stream,3)
    1 format(/,"The following ",i5," transitions are selected:",           &
        //,"     I-level-F     I--J^P--F      Transition Energy       ",   &
           "Orbital symmetries (emitted)",                                 &
         /,"                                     (in ",a4,")  ",           &
         /,4x,83("-") )
    2 format(4x,i4," -",i4,3x,a4,a1,3x,a4,a1,5x,1pe14.7,9x,10(a2,1x))
    3 format(4x,83("-") )
      !
   end subroutine auger_print_transitions
   !
   !
   subroutine auger_pure_matrix(csf_set,ii)
   !--------------------------------------------------------------------
   ! Calculates the 'pure' Auger matrix for the given configuration scheme
   ! csf_set. The first no_f CSF belong to the final-state representation
   ! and the following no_f+1,...,no_f+no_i to the initial states.
   ! The procedure takes into account the Coulomb interaction to the
   ! Auger matrix and, if required, also the Breit interaction.
   !
   ! Calls:
   !--------------------------------------------------------------------
      type(csf_basis), intent(in) :: csf_set
      integer, intent(in)         :: ii

      integer           :: i, ia, ib, ic, id, j, no_T_coeff, no_V_coeff, r, s, &
                           ss, t, xnu, Vnu, rrr, sss, asfi, asff
      character(len=30) :: xstring
      real(kind=dp)     :: aweight, xvalue, xweight
      type(nkappa)      :: aa, bb, cc, dd
      type(nkappa)      :: Va, Vb, Vc, Vd

      ! Local thread-private arrays for reduction
      real(kind=dp), dimension(:,:), allocatable :: local_matrix
      integer :: local_number_radial_integrals
      real(kind=dp), dimension(1:500) :: local_radial_value
      character(len=25), dimension(1:500) :: local_radial_string

      if (rabs_use_stop   .and.   csf_set%nocsf /= auger%no_f+auger%no_i) then
         stop "auger_pure_matrix(): program stop A."
      end if

      ! Initialize global variables
      number_radial_integrals = 0
      auger%matrix = zero

      ! Allocate thread-private matrix
      allocate(local_matrix(1:auger%no_f,1:auger%no_i))
      local_matrix = zero
      local_number_radial_integrals = 0
      local_radial_value = zero
      local_radial_string = ""

      if (auger_nonorthogonal_eval) then
         call nonorth_initialize_ci(csf_set,wave_final,wave_initial,auger_csp)
      end if

      ! Parallelize the outer loops over r and s
      !$OMP PARALLEL PRIVATE(r, s, rrr, ss, sss, t, Va, Vb, Vc, Vd, Vnu, &
      !$OMP                  aweight, ia, ib, ic, id, i, j, xweight, xnu, &
      !$OMP                  aa, bb, cc, dd, xvalue, xstring) &
      !$OMP          SHARED(csf_set, auger, wave_final, wave_initial, auger_csp)
      !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC) REDUCTION(+:local_matrix)
      do r = 1,auger%no_f
         do s = auger%no_f+1,auger%no_f+auger%no_i
            rrr = auger%ndx_f(r)
            ss  = s - auger%no_f
            sss = auger%ndx_i(ss)

            if (auger_nonorthogonal_eval) then
               Aoperator%particle = 2
               Aoperator%rank     = 0
               Aoperator%parity   = 1
               call nonorth_calculate_csf_pair(csf_set,Aoperator,r,s,no_V_coeff)
            else
               call anco_calculate_csf_pair(csf_set,r,s,no_T_coeff,no_V_coeff)
            end if

            ! Cycle over all angular coefficients
            do t = 1,no_V_coeff
               if (auger_nonorthogonal_eval) then
                  Va    = nonorth_V_list(t)%a
                  Vb    = nonorth_V_list(t)%b
                  Vc    = nonorth_V_list(t)%c
                  Vd    = nonorth_V_list(t)%d
                  Vnu   = nonorth_V_list(t)%nu
                  aweight = nonorth_V_list(t)%V
               else
                  Va = anco_V_list(t)%a
                  Vb = anco_V_list(t)%b
                  Vc = anco_V_list(t)%c
                  Vd = anco_V_list(t)%d
                  Vnu     = anco_V_list(t)%nu
                  aweight = anco_V_list(t)%V
               end if

               ia = 0; ib = 0; ic = 0; id = 0
               do i = 1,wave_final%number_of_rwf
                  if (wave_final%rwf(i)%orbital == Va) ia = i
                  if (wave_final%rwf(i)%orbital == Vb) ib = i
               end do
               do i = 1,wave_initial%number_of_rwf
                  if (wave_initial%rwf(i)%orbital == Vc) ic = i
                  if (wave_initial%rwf(i)%orbital == Vd) id = i
               end do

               if (Va%n < 0 .and. Vb%n > 0) then
                  local_matrix(r,ss) = local_matrix(r,ss) + aweight * &
                     XL_Coulomb_strength_grasp2k(Vnu,auger_csp, &
                        wave_final%rwf(ib),wave_initial%rwf(ic), &
                        wave_initial%rwf(id),.false.)

                  if (auger_print_radial_integrals) then
                     xweight = aweight
                     xnu     = Vnu
                     aa%n     = auger_csp%orbital%n
                     aa%kappa = auger_csp%orbital%kappa
                     bb%n     = wave_final%rwf(ib)%orbital%n
                     bb%kappa = wave_final%rwf(ib)%orbital%kappa
                     cc%n     = wave_initial%rwf(ic)%orbital%n
                     cc%kappa = wave_initial%rwf(ic)%orbital%kappa
                     dd%n     = wave_initial%rwf(id)%orbital%n
                     dd%kappa = wave_initial%rwf(id)%orbital%kappa
                     xvalue = XL_Coulomb_strength_grasp2k(Vnu,auger_csp, &
                        wave_final%rwf(ib),wave_initial%rwf(ic), &
                        wave_initial%rwf(id),.false.)
                     xstring = Xk_name(xnu,aa%n,aa%kappa,bb%n,bb%kappa, &
                        cc%n,cc%kappa,dd%n,dd%kappa)
                     !$OMP CRITICAL(print_radial)
                     print *, "(a) r, s:",trim(xstring),r,s,xweight,xvalue, &
                        xweight*xvalue,local_matrix(r,ss)
                     !$OMP END CRITICAL
                  end if

                  if (auger_include_breit) then
                     local_matrix(r,ss) = local_matrix(r,ss) + aweight * &
                        XL_Breit0_strength_grasp2k(Vnu,auger_csp, &
                           wave_final%rwf(ib),wave_initial%rwf(ic), &
                           wave_initial%rwf(id))
                  end if
               else if (Va%n > 0 .and. Vb%n < 0) then
                  local_matrix(r,ss) = local_matrix(r,ss) + aweight * &
                     XL_Coulomb_strength_grasp2k(Vnu,wave_final%rwf(ia), &
                        auger_csp,wave_initial%rwf(ic), &
                        wave_initial%rwf(id),.false.)

                  if (auger_print_radial_integrals) then
                     xweight = aweight
                     xnu     = Vnu
                     aa%n     = wave_final%rwf(ia)%orbital%n
                     aa%kappa = wave_final%rwf(ia)%orbital%kappa
                     bb%n     = auger_csp%orbital%n
                     bb%kappa = auger_csp%orbital%kappa
                     cc%n     = wave_initial%rwf(ic)%orbital%n
                     cc%kappa = wave_initial%rwf(ic)%orbital%kappa
                     dd%n     = wave_initial%rwf(id)%orbital%n
                     dd%kappa = wave_initial%rwf(id)%orbital%kappa
                     xvalue = XL_Coulomb_strength_grasp2k(Vnu, &
                        wave_final%rwf(ia),auger_csp, &
                        wave_initial%rwf(ic),wave_initial%rwf(id),.false.)
                     xstring = Xk_name(xnu,aa%n,aa%kappa,bb%n,bb%kappa, &
                        cc%n,cc%kappa,dd%n,dd%kappa)
                     !$OMP CRITICAL(print_radial)
                     print *, "(b) r, s:",trim(xstring),r,s,xweight,xvalue, &
                        xweight*xvalue,local_matrix(r,ss)
                     !$OMP END CRITICAL
                  end if

                  if (auger_include_breit) then
                     local_matrix(r,ss) = local_matrix(r,ss) + aweight * &
                        XL_Breit0_strength_grasp2k(Vnu,wave_final%rwf(ia), &
                           auger_csp,wave_initial%rwf(ic), &
                           wave_initial%rwf(id))
                  end if
               else
                  stop "auger_pure_matrix(): program stop B."
               end if

               ! Handle radial integrals printing
               if (auger_print_radial_integrals) then
                  if (abs(xvalue) >= 1.0e-4) then
                     !$OMP CRITICAL(radial_count)
                     do j = 1,number_radial_integrals
                        if (radial_string(j) == xstring) goto 1
                     end do
                     number_radial_integrals = number_radial_integrals + 1
                     radial_string(number_radial_integrals) = xstring
                     radial_value(number_radial_integrals) = xvalue
                   1 continue
                     !$OMP END CRITICAL
                  end if
               end if
            end do
         end do
      end do
      !$OMP END DO

      ! Reduce local matrices to global matrix
      !$OMP CRITICAL
      auger%matrix = auger%matrix + local_matrix
      !$OMP END CRITICAL
      !$OMP END PARALLEL

      deallocate(local_matrix)

      ! Print radial integrals (outside parallel region)
      if (auger_print_radial_integrals) then
         print *, " "
         print *, "Effective interactions strengths (|X^k()| > 10^-4) for "// &
                  "the current transition; continuum orbital has n=11"
         print *, "-------------------------------------------------------"// &
                  "--------------------------------------------------"
         print *, " "
         write(*,2) (trim(radial_string(j)),radial_value(j), &
                    j=1,number_radial_integrals)
       2 format(4(3x,a25," = ",f6.4,";"))
      end if
   end subroutine auger_pure_matrix
   !
   !
   function auger_return_amplitude(level_i,level_f,kappa)    result(amp)
   !--------------------------------------------------------------------
   ! Return the amplitude of the Auger channel (J_f kappa |V| J_i) from
   ! the (precalculated) Auger channels or zero, if this amplitude is not
   ! defined or calculated.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: level_i,level_f,kappa
      complex(kind=dp)    :: amp
      !
      integer :: i, j
      !
      amp = cmplx(zero,zero)
      !
      do  i = 1,number_of_transitions
         if (transition(i)%level_i == level_i  .and.  &
             transition(i)%level_f == level_f) then
            do  j = 1,transition(i)%number_of_channels
               if (transition(i)%channel(j)%kappa == kappa) then
                  amp = transition(i)%channel(j)%amplitude
               end if
            end do
         end if
      end do
      !
   end function auger_return_amplitude
   !
   !
   subroutine auger_rho_after_dipole(k,nl,J1_levels,matrix)
   !--------------------------------------------------------------------
   ! Calculates the density matrix of the intermediate states from the
   ! statistical tensor of the initial (resonant) state and the amplitudes
   ! of the first-step Auger decay. In the present case, however, the
   ! density matrix of the initial states assumes already a coherent
   ! dipole excitation of two J0 = 1 levels following the excitation with
   ! linearly polarized light along the z-axis.
   !
   ! Due to this assumption, the summation over J0 and J0' runs over two
   ! levels both with J0 = 1; the corresponding density is obtained from
   ! a call to the procedure auger_rho_dipole_excitation().
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in)                           :: k, nl
      integer, dimension(:), intent(in)             :: J1_levels
      complex(kind=dp), dimension(:,:), intent(out) :: matrix
      !
      integer                            :: i, i0, i0p, i1, i1p, kappa,   &
                                            j, J0, J0p, J1, J1p, J0_lev,  &
                                            J0p_lev, l, nl0
      integer, dimension(1:10)           :: J0_levels
      real(kind=dp)                      :: Gamma_J0_J0p, delta_e
      complex(kind=dp)                   :: wa
      complex(kind=dp), dimension(10,10) :: rho_k0_J0_J0p
      !
      J0_levels(:) = 0
      J0_levels(1) = 2;   J0_levels(2) = 3;   nl0 = 2
      J0           = 2;   J0p          = 2
      !
      ! Set-up the statistical tensor of the coherently-excited initial
      ! states following a dipole exciation with linearly-polarized light
      call auger_rho_dipole_excitation(k,nl0,J0_levels,rho_k0_J0_J0p)
      !
      matrix(:,:) = cmplx(zero,zero)
      !
      do  i1 = 1,nl
         do  i1p = 1,nl
            do  i = 1,number_of_transitions
               if (transition(i)%level_f == J1_levels(i1)) then
                  J1  = transition(i)%totalJ_f
               end if
               if (transition(i)%level_f == J1_levels(i1p)) then
                  J1p = transition(i)%totalJ_f
               end if
            end do
            !!x print *, "J1, J1p = ",J1, J1p
            !
            do  kappa = -12,12
               if (kappa == 0) cycle
               j  = angular_momentum_j(kappa)
               l  = angular_momentum_l(kappa)
               ! Sum over all contribution from J0 and J0p of the initial-state
               ! density matrix
               do  i0 = 1,nl0
                  do  i0p = 1,nl0
                  J0_lev = J0_levels(i0);   J0p_lev = J0_levels(i0p)
                  Gamma_J0_J0p = half *                                       &
                                 (auger_width(J0_lev)+auger_width(J0p_lev))
                  delta_e      = auger_delta_energy(J0_lev,J0p_lev,asf_initial)
                  !
                  wa = wigner_6j_symbol(J0,J1,j,J1p,J0p,k+k)                * &
                  auger_return_amplitude(J0_levels(i0),J1_levels(i1),kappa) * &
                  conjg(auger_return_amplitude(J0_levels(i0p),                &
                                               J1_levels(i1p),kappa))       * &
                  sqrt((J0+one)*(J0p+1)) *(Gamma_J0_J0p-cmplx(zero,delta_e))/ &
                  (delta_e*delta_e + Gamma_J0_J0p*Gamma_J0_J0p)
                  !
                  if (abs(wa) > eps10) then
                     wa = rho_k0_J0_J0p(i0,i0p) * wa
                     if (mod(j+J0+J1,4) == 0) then
                        matrix(i1,i1p) = matrix(i1,i1p) + wa
                     else if (mod(j+J0+J1,4) == 2) then
                        matrix(i1,i1p) = matrix(i1,i1p) - wa
                     else if (rabs_use_stop) then
                        stop "auger_rho_after_dipole(): program stop A."
                     end if
                  end if
                  !
                  end do
               end do
               !
            end do
            !
         end do
      end do
      !
   end subroutine auger_rho_after_dipole
   !
   !
   subroutine auger_rho_density(k,J0_level,nl,J1_levels,matrix)
   !--------------------------------------------------------------------
   ! Calculates the density matrix of the intermediate states from the
   ! statistical tensor of the initial (resonant) state and the amplitudes
   ! of the first-step Auger decay [cf. Ueda et al., JPB (1999), Eq. (6)].
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in)                           :: k, J0_level, nl
      integer, dimension(:), intent(in)             :: J1_levels
      complex(kind=dp), dimension(:,:), intent(out) :: matrix
      !
      integer          :: i, i1, i1p, kappa, j, J0, J1, J1p, l
      real(kind=dp)    :: rho_k0_J0
      complex(kind=dp) :: wa
      !
      ! Define the statistical tensor of the initial (resonant) state;
      ! this has to be modified in other cases
      if (k == 0) then
         rho_k0_J0 = 1 / sqrt(three)
      else if (k == 2) then
         rho_k0_J0 = - sqrt(two/three)
      else if (rabs_use_stop) then
         stop "auger_rho_density(): program stop A."
      end if
      !
      matrix(:,:) = cmplx(zero,zero)
      !
      do  i1 = 1,nl
         do  i1p = 1,nl
            do  i = 1,number_of_transitions
               if (transition(i)%level_i == J0_level) then
                  J0  = transition(i)%totalJ_i
               end if
               if (transition(i)%level_f == J1_levels(i1)) then
                  J1  = transition(i)%totalJ_f
               end if
               if (transition(i)%level_f == J1_levels(i1p)) then
                  J1p = transition(i)%totalJ_f
               end if
            end do
            !!x print *, "J0, J1, J1p = ",J0, J1, J1p
            !
            do  kappa = -12,12
               if (kappa == 0) cycle
               j  = angular_momentum_j(kappa)
               l  = angular_momentum_l(kappa)
               wa = wigner_6j_symbol(J0,J1,j,J1p,J0,k+k)                     * &
                    auger_return_amplitude(J0_level,J1_levels(i1),kappa)     * &
                    conjg(auger_return_amplitude(J0_level,J1_levels(i1p),kappa))
               if (abs(wa) > eps10) then
                  wa = rho_k0_J0 * (J0 + one) * wa
                  if (mod(j+J0+J1,4) == 0) then
                     matrix(i1,i1p) = matrix(i1,i1p) + wa
                  else if (mod(j+J0+J1,4) == 2) then
                     matrix(i1,i1p) = matrix(i1,i1p) - wa
                  else if (rabs_use_stop) then
                     stop "auger_rho_density(): program stop B."
                  end if
               end if
            end do
            !
         end do
      end do
      !
   end subroutine auger_rho_density
   !
   !
   subroutine auger_rho_dipole_excitation(k,nl,J0_levels,matrix)
   !--------------------------------------------------------------------
   ! Calculates the density matrix of the coherently-excited initial
   ! states of an Auger cascade by linearly polarized light along the
   ! z-axis. Here, we assume that the excitation occurs from a noble
   ! gas configuration with J_g = 0 and that all J0_levels levels have
   ! total angular momentum J0 = 1. The amplitudes must be set
   ! explicitly inside of this procedure.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in)                           :: k, nl
      integer, dimension(:), intent(in)             :: J0_levels
      complex(kind=dp), dimension(:,:), intent(out) :: matrix
      !
      integer          :: i0, i0p
      complex(kind=dp),dimension(1:10) :: dipole_amplitudes
      complex(kind=dp) :: di, dip
      real(kind=dp)    :: awd2, awd3
      logical, save    :: first = .true.
      !
      awd2 = zero;    awd3 = zero
      dipole_amplitudes(:) = zero
      !
      ! Xenon
      ! Babushkin gauge
      dipole_amplitudes(2) = cmplx(-1.3686616e-3_dp,zero)
      dipole_amplitudes(3) = cmplx( 8.1854684e-5_dp,zero)
      !
      ! Coulomb gauge
      ! dipole_amplitudes(2) = cmplx( 1.6357207e-3_dp,zero)
      ! dipole_amplitudes(3) = cmplx(-9.5924369e-5_dp,zero)
      !
      if (first) then
      !t awd2 = auger_width_dipole(2)
      !t awd3 = auger_width_dipole(3)
      !
      print *,    " "
      print *,    "Excitation amplitudes and widths of the initial states"
      print *,    "------------------------------------------------------"
      print *,    "Level 2: ",dipole_amplitudes(2), awd2
      print *,    "Level 3: ",dipole_amplitudes(3), awd3
      print *,    " "
      !
      write(24,*) " "
      write(24,*) "Excitation amplitudes and  of the initial states"
      write(24,*) "------------------------------------------------------"
      write(24,*) "Level 2: ",dipole_amplitudes(2), awd2
      write(24,*) "Level 3: ",dipole_amplitudes(3), awd3
      write(24,*) " "
      first = .false.
      end if
      !
      matrix(:,:) = zero
      !
      do  i0 = 1,nl
         do  i0p = 1,nl
            di  = dipole_amplitudes(J0_levels(i0))
            dip = dipole_amplitudes(J0_levels(i0p))
            if (di == zero  .or.   dip == zero) then
               print *, "di, dip = ",di, dip
               ! stop "auger_rho_dipole_excitation(): program stop A."
            end if
            !
            if (k == 0) then
               matrix(i0,i0p) = sqrt(one/three) * di * conjg(dip)
            else if (k == 2) then
               matrix(i0,i0p) = - sqrt(two/three) * di * conjg(dip)
            else
               stop "auger_rho_dipole_excitation(): program stop B."
            end if
            !
         end do
      end do
      !
   end subroutine auger_rho_dipole_excitation
   !
   !
   subroutine auger_set_overlaps()
   !--------------------------------------------------------------------
   ! Calculates the non-orthogonal and overlap integrals between all
   ! bound orbitals of the corresponding final- and initial-state arrays.
   !
   ! Calls: rk_integral_grasp2k_ab().
   !--------------------------------------------------------------------
      !
      integer :: i1, i2, kappa, kappa_min, kappa_max, pqn_i, pqn_i_max, &
                 pqn_f, pqn_f_max
      !
      ! Allocate an appropriate array for the overlap integrals
      kappa_min = 100;  kappa_max = -100;  pqn_i_max = -100;  pqn_f_max = -100
      do  i1 = 1,wave_initial%number_of_rwf
         kappa_min = min(kappa_min, wave_initial%rwf(i1)%orbital%kappa)
         kappa_max = max(kappa_max, wave_initial%rwf(i1)%orbital%kappa)
         pqn_i_max = max(pqn_i_max, wave_initial%rwf(i1)%orbital%n)
      end do
      do  i1 = 1,wave_final%number_of_rwf
         kappa_min = min(kappa_min, wave_final%rwf(i1)%orbital%kappa)
         kappa_max = max(kappa_max, wave_final%rwf(i1)%orbital%kappa)
         pqn_f_max = max(pqn_f_max, wave_final%rwf(i1)%orbital%n)
      end do
      !
      allocate( auger_overlap(kappa_min:kappa_max,pqn_f_max,pqn_i_max) )
      !
      ! Calculate the non-orthogonal and overlap integrals
      do  i1 = 1,wave_final%number_of_rwf
         do  i2 = 1,wave_initial%number_of_rwf
            if (wave_final%rwf(i1)%orbital%kappa == &
                wave_initial%rwf(i2)%orbital%kappa) then
               kappa = wave_final%rwf(i1)%orbital%kappa
               pqn_f = wave_final%rwf(i1)%orbital%n
               pqn_i = wave_initial%rwf(i2)%orbital%n
               auger_overlap(kappa,pqn_f,pqn_i) = &
                  rk_integral_grasp2k_ab(wave_final,wave_initial,0,i1,i2)
            end if
         end do
      end do
      !
      ! Print the overlap integrals
      do  i1 = 1,wave_final%number_of_rwf
         do  i2 = 1,wave_initial%number_of_rwf
            if (wave_final%rwf(i1)%orbital%kappa == &
                wave_initial%rwf(i2)%orbital%kappa) then
               kappa = wave_final%rwf(i1)%orbital%kappa
               pqn_f = wave_final%rwf(i1)%orbital%n
               pqn_i = wave_initial%rwf(i2)%orbital%n
            end if
         end do
      end do
      !
   end subroutine auger_set_overlaps
   !
   !
   subroutine auger_set_transitions()
   !--------------------------------------------------------------------
   ! Determines how many and which transitions need to be calculated and
   ! initializes the array transitions of type(auger_transitions).
   ! The default is (for number_of_transitions == 0) that all transitions
   ! with a positive energy are calculated; for number_of_transitions /= 0,
   ! individual transitions were selected during input time and will be
   ! initialized instead; in this case also negative transition energies
   ! are allowed and may be overritten by appropriate experimental energies.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer       :: asfi, asff, i, imin, j, k, m, nchannels, nt
      real(kind=dp) :: energy, energy_exp, energy_min
      integer, dimension(20)       :: kappa_channel
      !x real(kind=dp), dimension(1:n_grasp2k) :: bessel
      !
      ! Determine the total number of transitions
      nt = 0
      if (number_of_transitions == 0) then
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               if (asf_initial%asf(i)%energy - asf_final%asf(j)%energy &
                       + auger_energy_shift > auger_lowest_transition) then
                  nt = nt + 1
               end if
            end do
         end do
      else
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_transitions
                  if ( (select_level_i(k) == asf_initial%asf(i)%level_No .or. &
                        select_level_i(k) == 0)                         .and. &
                       (select_level_f(k) == asf_final%asf(j)%level_No   .or. &
                        select_level_f(k) == 0)                         .and. &
                       (asf_initial%asf(i)%energy - asf_final%asf(j)%energy   &
                         + auger_energy_shift > auger_lowest_transition) ) then
                     nt = nt + 1
                  end if
               end do
            end do
         end do
      end if
      !
      ! Allocate the derived data structure transition
      allocate( transition(1:nt) )
      !
      ! Now initialize all transitions
      nt = 0
      if (number_of_transitions == 0) then
         !$OMP PARALLEL DO REDUCTION(+:nt)
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               if (asf_initial%asf(i)%energy - asf_final%asf(j)%energy  &
                        + auger_energy_shift > auger_lowest_transition) then
                  nt = nt + 1
                  transition(nt)%asfi     = i;   transition(nt)%asff = j
                  transition(nt)%level_i  = asf_initial%asf(i)%level_No
                  transition(nt)%level_f  = asf_final%asf(j)%level_No
                  transition(nt)%energy   = asf_initial%asf(i)%energy - &
                                            asf_final%asf(j)%energy
                  !
                  ! Read in experimental energies if appropriate
                  if (auger_apply_exp_energies) then
                     if (energy_inverse) then
                        energy = energy_factor / transition(nt)%energy
                     else
                        energy = energy_factor * transition(nt)%energy
                     end if
                   1 print *, "Transition",transition(nt)%level_i,"-",                 &
                      transition(nt)%level_f,"has ab-initio energy E_theo = ",&
                      energy,trim(energy_unit),";"
                     print *, " enter E_exp (in ",trim(energy_unit), &
                      ") or 0.0 to use E_theo:"
                     read(*,*,err=1) energy_exp
                     if (energy_exp /= 0) then
                        if (energy_inverse) then
                        transition(nt)%energy = energy_factor/energy_exp
                        else
                        transition(nt)%energy = energy_exp/energy_factor
                        end if
                     end if
                  end if
               end if
            end do
         end do
         !$OMP END PARALLEL DO
      else
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_transitions
                  if ( (select_level_i(k) == asf_initial%asf(i)%level_No .or. &
                        select_level_i(k) == 0)                         .and. &
                       (select_level_f(k) == asf_final%asf(j)%level_No   .or. &
                        select_level_f(k) == 0)                         .and. &
                       (asf_initial%asf(i)%energy - asf_final%asf(j)%energy   &
                        + auger_energy_shift  > auger_lowest_transition) ) then
                     nt = nt + 1
                     transition(nt)%asfi     = i;   transition(nt)%asff = j
                     transition(nt)%level_i  = asf_initial%asf(i)%level_No
                     transition(nt)%level_f  = asf_final%asf(j)%level_No
                     transition(nt)%energy   = asf_initial%asf(i)%energy - &
                                               asf_final%asf(j)%energy
                     !
                     ! Read in experimental energies if appropriate
                     if (auger_apply_exp_energies) then
                        if (energy_inverse) then
                           energy = energy_factor / transition(nt)%energy
                        else
                           energy = energy_factor * transition(nt)%energy
                        end if
                      2 print *, "Transition",transition(nt)%level_i,"-",                 &
                           transition(nt)%level_f,            &
                           "has ab-initio energy E_theo = ",  &
                           energy,trim(energy_unit),";"
                        print *, " enter E_exp (in ",trim(energy_unit), &
                           ") or 0.0 to use E_theo:"
                        read(*,*,err=2) energy_exp
                        if (energy_exp /= 0) then
                           if (energy_inverse) then
                           transition(nt)%energy=energy_factor/energy_exp
                           else
                           transition(nt)%energy=energy_exp/energy_factor
                           end if
                        end if
                     end if
                     !
                  end if
               end do
            end do
         end do
      end if
      !
      number_of_transitions = nt
      !
      ! Order the transitions in ascending order of transition energies
      if (auger_sort_transition_energy) then
         do  i = 1,number_of_transitions-1
            imin       = i
            energy_min = 1.0e20
            do  j = i,number_of_transitions
               if (transition(j)%energy < energy_min) then
                  imin = j;   energy_min = transition(j)%energy
               end if
            end do
            if (imin /= i) then
               asfi     = transition(i)%asfi
               asff     = transition(i)%asff
               energy   = transition(i)%energy
               transition(i)%asfi      = transition(imin)%asfi
               transition(i)%asff      = transition(imin)%asff
               transition(i)%energy    = transition(imin)%energy
               transition(imin)%asfi   = asfi
               transition(imin)%asff   = asff
               transition(imin)%energy = energy
            end if
         end do
      end if
      !
      ! Now assign the full transition specifications and calculate the
      ! necessary Bessel functions
      do  nt = 1,number_of_transitions
         i = transition(nt)%asfi;   j = transition(nt)%asff
         transition(nt)%level_i  = asf_initial%asf(i)%level_No
         transition(nt)%level_f  = asf_final%asf(j)%level_No
         transition(nt)%totalJ_i = asf_initial%asf(i)%totalJ
         transition(nt)%totalJ_f = asf_final%asf(j)%totalJ
         transition(nt)%parity_i = asf_initial%asf(i)%parity
         transition(nt)%parity_f = asf_final%asf(j)%parity
         !
         ! Modify the transition energies if required
         transition(nt)%energy   = transition(nt)%energy + auger_energy_shift
	 call get_kappa_channels(                                   &
	          transition(nt)%totalJ_i,transition(nt)%parity_i,  &
	          transition(nt)%totalJ_f,transition(nt)%parity_f,  &
		  kappa_channel,nchannels)
         if (rabs_use_stop   .and.   nchannels == 0) then
            stop "auger_set_transitions(): program stop A."
         end if
         !
         ! Determine the number of channels
         i = 0
         do  m = 1,nchannels
            if (abs(kappa_channel(m)) <= auger_maximal_kappa) then
               i = i + 1
            end if
         end do
         !
         transition(nt)%number_of_channels = i
         allocate( transition(nt)%channel(1:i) )
         i = 0
         do  m = 1,nchannels
            if (abs(kappa_channel(m)) <= auger_maximal_kappa) then
               i = i + 1
               transition(nt)%channel(i)%kappa      = kappa_channel(m)
               transition(nt)%channel(i)%amplitude  = cmplx(zero,zero)
            end if
         end do
         !
      end do
      !
      print *, "  ",number_of_transitions,                         &
               " transitions have been initialized and will be "// &
               "calculated in this run of the program."
      !
   end subroutine auger_set_transitions
   !
   !
   function auger_spin_parameter(k,trans)                    result(eta)
   !--------------------------------------------------------------------
   ! Returns the eta_k spin angular coefficient for the transition
   ! which is given in terms of a number of individual channels.
   ! It uses the explicit expression in terms of the partial wave
   ! transition amplitudes as given by J Tulkki, H Aksela, and
   ! N M Kabachnik (Phys. Rev. A 48, 1277 (1993), eq. A27).
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in)    :: k
      type(auger_transition) :: trans
      real(kind=dp)          :: eta
      !
      integer          :: iphase, ja, jb, la, lb, ma, mb
      real(kind=dp)    :: norm, sum, wa
      !
      ! Calculate the normalization and general phase factor
      norm   = trans%probability / (two * pi)
      iphase = trans%totalJ_i + trans%totalJ_f - 1 + 32
      if (mod(iphase,4) == 2) then
         norm = -norm
      else if (rabs_use_stop       .and.                    &
              (mod(iphase,4) == 1  .or.  mod(iphase,4) == 3)) then
         stop "auger_spin_parameter(): program stop A."
      end if
      !
      sum = zero
      !
      ! Now cycle about the different Auger amplitudes including proper
      ! phase factors and the products of partial-wave expansion coefficients;
      ! Here, b denotes the prime quantum numbers in eq. A27.
      do  ma = 1,trans%number_of_channels
         do  mb = 1,trans%number_of_channels
	    ja = angular_momentum_j(trans%channel(ma)%kappa)
	    jb = angular_momentum_j(trans%channel(mb)%kappa)
	    la = angular_momentum_l(trans%channel(ma)%kappa)
	    lb = angular_momentum_l(trans%channel(mb)%kappa)
	    !
	    wa = (la+la+one)*(lb+lb+one)*(ja+one)*(jb+one) &
                *(trans%totalJ_i+one)
	    wa = sqrt(wa)
	    iphase = 1 + la + la - jb + 32
            if (mod(iphase,4) == 2) then
               wa = -wa
            else if (rabs_use_stop       .and.                    &
                    (mod(iphase,4) == 1  .or.  mod(iphase,4) == 3)) then
               stop "auger_spin_parameter(): program stop B."
            end if
	    !
	    wa = wa * Clebsch_Gordan(la+la,0,lb+lb,0,k+k,0)                   &
                    * wigner_6j_symbol(trans%totalJ_i,trans%totalJ_i,         &
	                               k+k,ja,jb,trans%totalJ_f)              &
	            * wigner_9j_symbol(la+la,lb+lb,k+k,1,1,2,ja,jb,k+k)       &
		    * aimag(trans%channel(ma)%amplitude *                     &
                            conjg(trans%channel(mb)%amplitude))
	    sum = sum + wa
	 end do
      end do
      !
      eta = sqrt(12.0_dp) * sum / norm
      !
   end function auger_spin_parameter
   !
   !
   subroutine auger_transition_properties(trans)
   !--------------------------------------------------------------------
   ! Calculates all selected Auger properties of transition i from the
   ! amplitudes of the individual channels. This concerns the total
   ! probability as well as angular distribution and spin polarization
   ! parameters.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      type(auger_transition), intent(inout) :: trans
      !
      integer :: j, stream
      real(kind=dp) :: energy, rate, sum
      !
      sum = zero
      do  j = 1,trans%number_of_channels
         sum = sum + trans%channel(j)%amplitude * &
               conjg(trans%channel(j)%amplitude)
      end do
      trans%probability = two * pi * sum
      !
      ! Calculate the angular distribution and spin-polarization coefficients
      if (auger_calc_angular_parameter) then
         trans%alpha_2 = auger_angular_parameter(2,trans)
         trans%alpha_4 = auger_angular_parameter(4,trans)
      end if
      !
      if (auger_calc_spin_polarization) then
         trans%eta_2 = auger_spin_parameter(2,trans)
         trans%eta_4 = auger_spin_parameter(4,trans)
      end if
      !
      ! Print a short summary of the transition if required
      !
      if (auger_print_each_transition) then
         stream = 6
       3 write(stream,*) " "
	 write(stream,*) " "
         write(stream,*) " Results for Auger transition ",               &
	                 trans%level_i," - ",trans%level_f,": ",         &
                         trim(angular_momentum_string(trans%totalJ_i)),  &
		         " ",trans%parity_i,"   ----> ",                 &
	                 trim(angular_momentum_string(trans%totalJ_f)),  &
		         " ",trans%parity_f
	 write(stream,*) " --------------------------------------------"//  &
	                 "------------------"
	 write(stream,*) " "
         if (energy_inverse) then
            energy = energy_factor / trans%energy
         else
            energy = energy_factor * trans%energy
         end if
	 rate = trans%probability * auger_rate_factor
	 write(stream,*) "   Energy     = ",energy,trim(energy_unit)
	 write(stream,*) "   Total rate = ",rate,trim(auger_rate_unit)
         !
         if (auger_calc_angular_parameter) then
	    write(stream,*) "   Alpha_2    = ",trans%alpha_2
	    write(stream,*) "   Alpha_4    = ",trans%alpha_4
	 end if
	 !
         if (auger_calc_spin_polarization) then
	    write(stream,*) "   Eta_2      = ",trans%eta_2
	    write(stream,*) "   Eta_4      = ",trans%eta_4
	 end if
	 write(stream,*) " "
	 write(stream,1) trim(auger_rate_unit)
	 write(stream,*) &
                " ------------------------------------------------------",   &
		"---------------------"
	 do  j = 1,trans%number_of_channels
	    rate = two * pi * trans%channel(j)%amplitude *                   &
	                conjg(trans%channel(j)%amplitude)* auger_rate_factor
	    write(stream,2) orbital_symmetry(trans%channel(j)%kappa),        &
                            trans%channel(j)%amplitude,                      &
			    trans%channel(j)%amplitude_re,rate,              &
		            mod(trans%channel(j)%phase,two*pi)
	 end do
	 write(stream,*) " "
       1 format("  Kappa            Amplitude        Real-Amplitude ",       &
                "  Rate (",a,")     Phase")
       2 format(3x,a4,5x,1pe10.3,2x,1pe10.3,4x,1pe10.3,4x,1pe10.3,4x,1pe10.3)
         !
         ! Re-cycle once on stream 24
         if (stream == 6) then
            stream = 24;   goto 3
         end if
      end if
      !
   end subroutine auger_transition_properties
   !
   !
   function auger_width(level_i)                           result(width)
   !--------------------------------------------------------------------
   ! Return the width of the 'initial' level_i as been derived from the
   ! precalculated transition probabilities into lower-lying levels.
   ! It assumes that all important decay branches are included in the
   ! computations.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: level_i
      real(kind=dp)       :: width
      !
      integer :: i
      !
      width = zero
      !
      do  i = 1,number_of_transitions
         if (transition(i)%level_i == level_i  .and.  &
             transition(i)%energy > zero) then
            width = width + transition(i)%probability
         end if
      end do
      !
   end function auger_width
   !
end module rabs_auger
