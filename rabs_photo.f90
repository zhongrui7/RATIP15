module rabs_photo
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains the procedures which are specific to the PHOTO
! program. This program supports the calculation of photoionization cross
! sections, relative cross sections angular distribution parameters.
! The continuum orbitals which are required for these computations are
! generated automatically by a call to the COWF component.
! At the present, the program supports nonorthogonality between the
! orbitals only via the form of the radial functions. In the evaluation
! of the many-electron matrix elements, by contrast, orthogonality is
! 'assumed' and, hence, standard Racah algebra techniques applied. See
! the component ANCO for calculating the angular integrals for non-scalar
! one-particle matrix elements.
! In this module below, there are several procedures also for the file
! handling and the intermediate storage to accelerate some of the
! computations.
!-----------------------------------------------------------------------
   !
   use rabs_anco
   use rabs_angular
   use rabs_constant
   use rabs_cowf
   use rabs_csl
   use rabs_determinant
   use rabs_dirac_orbital
   use rabs_grasp2k
   use rabs_input_dialog
   use rabs_mcp
   use rabs_nonorthonormal
   use rabs_nucleus
   use rabs_print
   use rabs_rcfp
   use rabs_recoupling
   implicit none
   private
   !
   !
   private :: photo_alignment_parameter
                 ! Returns the alignment parameter for a
		 ! given transition.
   private :: photo_angular_parameter
                 ! Returns the beta_k angular distribution parameter for a
		 ! given transition.
   private :: photo_C_factor
                 ! Calculates the C(kappa,kappa',L) factor as defined by
                 ! Huang (1980).
   public  :: photo_calculate_amplitudes
                 ! Calculates for all selected 'photon energies' (in turn) the
                 ! required continuum spinors and photoionization amplitudes.
   private :: photo_calculate_Bessel
                 ! Calculates the Bessel function over the radial grid
                 ! for a given factor.
   private :: photo_channel_amplitude
                 ! Calculates the photoionization amplitude of a channel from
                 ! the 'pure' photoionization matrix and the corresponding
                 ! mixing coefficients.
   public  :: photo_collect_input
                 ! Collects and proceeds all input for the calculation of
		 ! photoionization cross sections and angular distribution
                 ! parameters.
   public  :: photo_initialize
		 ! Set up the selected transitions and initializes some
		 ! arrays as required.
   private :: photo_kabachnik_B
                 ! Returns the B(k1,ki,k_gamma) parameter due to Kabachnik
		 ! (2007).
   private :: photo_kabachnik_5k_append_table
                 ! Append a table of Kabachnik's B-5k coefficients to the
		 ! results.
   private :: photo_kabachnik_5k_coefficient
                 ! Returns the B-5k coefficient as defined by Nicolai Kabachnik
		 ! (2006) for the collaboration with Fabiana (Italy).
   private :: photo_line_properties
                 ! Calculates all selected photoionization properties of line i
		 ! from the amplitudes of the individual channels.
   public  :: photo_print_amplitudes
                 ! Writes out the information about all (selected)
                 ! photoionzation channels and amplitudes to a (.chn) file
   private :: photo_print_lines
                 ! Prints all selected 'photon' or 'electron' lines in a neat
                 ! format before the computation starts.
   public  :: photo_print_results
                 ! Writes the cross sections and angular distribution
                 ! parameters to the .sum file.
   private :: photo_pure_matrix
                 ! Calculates the 'pure' photoionization matrix for the given
		 ! configuration scheme.
   private :: photo_reduced_M_integral
                 ! Calculates a M integral for the coupling of the radiation
                 ! field for a given multipolarity and gauge form.
   private :: photo_set_lines
                 ! Determines which ionization lines need to be calculated
                 ! and initializes storage for an appropriate data type for
                 ! them.
   private :: photo_set_overlaps
                 ! Initializes the array of overlap integrals for the
                 ! calculation of relaxed ionization cross sections.
   private :: photo_spin_parameter
                 ! Returns the eta_k spin polarization coefficient for a
		 ! given photoionization line.
   !
   ! Define some global data of the PHOTO program; most of these data are
   ! read in during the interactive control and may overwrite existing
   ! default values
   !
   ! Storage for the initial and final atomic states and wave functions
   !
   ! Define a 'configuration scheme' to be built up during execution
   type(grasp2k_orbital), public :: wave_continuum
   type(orbital_function)        :: photo_csp
   !
   ! Define an internal structure type(photo_transition) which stores all
   ! necessary information for an photoionization line
   type :: photo_channel
      integer          :: kappa, totalJ
      character(len=1) :: parity
      character(len=2) :: multipole
      character(len=9) :: gauge
      real(kind=dp)    :: phase, amplitude_re
      complex(kind=dp) :: amplitude
   end type photo_channel
   !
   type :: photo_line
      integer          :: asfi, asff, level_i, level_f, totalJ_i, totalJ_f
      integer          :: No_channels
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: p_energy, e_energy, cs_coulomb, cs_babushkin
      complex(kind=dp) :: beta_b, beta_c, xi_b, xi_c, eta_b, eta_c, &
                          zeta_b, zeta_c, alignment_b, alignment_c
      real(kind=dp), dimension(:), pointer :: bessel0, bessel1, bessel2, &
                                              bessel3, bessel4, bessel5
      type(photo_channel), dimension(:), pointer :: channel
   end type photo_line
   !
   type(photo_line), dimension(:), allocatable :: line
   !
   type :: photo_matrix
      integer :: No_f, No_i
      integer, dimension(:), pointer         :: ndx_f, ndx_i
      real(kind=dp), dimension(:,:), pointer :: matrix
   end type photo_matrix
   !
   type(photo_matrix) :: photo
   !
   integer          :: number_of_lines                       = 0,  &
                       number_of_multipoles                  = 0,  &
                       number_of_photo_energies              = 0,  &
                       number_of_atransitions                = 0
   !
   integer, public  :: photo_maximal_kappa                   = 20
   !
   ! Define global logical flags for the control of the PHOTO program; the
   ! default values for these flags may be overwritten interactively during
   ! input time
   logical, public :: photo_add_1                     = .false.,  &
                      photo_apply_exp_energies        = .false.,  &
                      photo_calc_angular_parameter    = .true.,   &
                      photo_calc_spin_polarization    = .false.,  &
		      photo_kabachnik_5k_coeff        = .true.,   &
		      photo_nonorthogonal_eval        = .false.,  &
		      photo_print_csf_scheme          = .false.,  &
		      photo_print_each_line           = .true.,   &
		      photo_print_cs_in_hartree       = .true.,   &
                      photo_print_selected_lines      = .true.,   &
                      photo_print_chn_file            = .false.,  &
                      photo_sort_line_energy          = .true.,   &
                      photo_use_formatted_mix_file    = .true.,   &
                      photo_use_formatted_rwf_file    = .true.
   !
   ! Energy unit for the output of all energies
   real(kind=dp)    :: photo_cs_factor = zero, &
                       photo_print_cut = 0.001_dp, &
                       photo_energy_shift  = zero, photo_maximal_energy
   character(len=7) :: photo_cs_unit
   !
   ! Define storage for the overlap integrals: (kappa,pqn_f,pqn_i)
   real(kind=dp), dimension(:,:,:), allocatable :: photo_overlap
   !
   ! Define some variables and arrays for processing input data from
   ! photo_collect_input()
   !! integer, dimension(200)         :: select_level_i, select_level_f
   real(kind=dp), dimension(1000)   :: photo_energy_selection
   character(len=2), dimension(20) :: photo_multipole  = "  "
   !
contains
   !
   function photo_alignment_parameter(g,keyword,tline)           result(A_2)
   !--------------------------------------------------------------------
   ! Returns the A_2 alignment parameter for the photoionization
   ! line tline.
   !
   ! Calls: photo_kabachnik_B
   !--------------------------------------------------------------------
      !
      character(len=1), intent(in) :: g
      character(len=9), intent(in) :: keyword
      type(photo_line), intent(in) :: tline
      complex(kind=dp)             :: A_2
      !
      A_2 = cmplx(zero,zero)
      !
      select case (keyword)
      case ("alignment")
         !
	 A_2 = - sqrt(two) * photo_kabachnik_B(g,tline,0,2,2) /  &
	                     photo_kabachnik_B(g,tline,0,0,0)
         !
      case ("beta     ")
         !
	 A_2 = - sqrt(two) * photo_kabachnik_B(g,tline,2,0,2) /  &
	                     photo_kabachnik_B(g,tline,0,0,0)
         !
      case default
         stop "photo_alignment_parameter(): program stop B."
      end select
   end function photo_alignment_parameter
   !
   !
   function photo_angular_parameter(g,keyword,tline)             result(beta)
   !--------------------------------------------------------------------
   ! Returns the beta angular distribution parameter for the photoionization
   ! line tline. It uses the explicit expression in terms of the partial
   ! wave photoionization amplitudes as given by Huang, PRA 22 (1980) 223.
   !
   ! Calls: angular_momentum_j(), angular_momentum_l().
   !--------------------------------------------------------------------
      !
      character(len=1), intent(in) :: g
      character(len=4), intent(in) :: keyword
      type(photo_line), intent(in) :: tline
      complex(kind=dp)             :: beta
      !
      integer                      :: i, ip, j, jp
      complex(kind=dp)             :: N
      !
      beta = cmplx(zero,zero)
      !
      select case (keyword)
      case ("beta")
         !
         select case(g)
         case("B")
            N = cmplx(zero,zero)
            do  i = 1,tline%No_channels
               if (tline%channel(i)%multipole == "E1"   .and.    &
                   tline%channel(i)%gauge     == "Babushkin") then
                  N = N + tline%channel(i)%amplitude *  &
                    conjg(tline%channel(i)%amplitude)
               end if
            end do
            !
            do  i = 1,tline%No_channels
               if (tline%channel(i)%multipole == "E1"   .and.    &
                   tline%channel(i)%gauge     == "Babushkin") then
                  j = angular_momentum_j(tline%channel(i)%kappa)
                  do  ip = 1,tline%No_channels
                     if (tline%channel(ip)%multipole == "E1"   .and.    &
                         tline%channel(ip)%gauge     == "Babushkin") then
                        jp = angular_momentum_j(tline%channel(ip)%kappa)
                        beta = beta + wigner_3j_symbol(j,jp,4,1,-1,0) * &
                               photo_C_factor(2,ip,i,tline)
                     end if
                  end do
               end if
            end do
            !
         case("C")
            N = cmplx(zero,zero)
            do  i = 1,tline%No_channels
               if (tline%channel(i)%multipole == "E1"   .and.    &
                   tline%channel(i)%gauge     == "Coulomb  ") then
                  N = N + tline%channel(i)%amplitude *  &
                    conjg(tline%channel(i)%amplitude)
               end if
            end do
            !
            do  i = 1,tline%No_channels
               if (tline%channel(i)%multipole == "E1"   .and.    &
                   tline%channel(i)%gauge     == "Coulomb  ") then
                  j = angular_momentum_j(tline%channel(i)%kappa)
                  do  ip = 1,tline%No_channels
                     if (tline%channel(ip)%multipole == "E1"   .and.    &
                         tline%channel(ip)%gauge     == "Coulomb  ") then
                        jp = angular_momentum_j(tline%channel(ip)%kappa)
                        beta = beta + wigner_3j_symbol(j,jp,4,1,-1,0) * &
                               photo_C_factor(2,ip,i,tline)
                     end if
                  end do
               end if
            end do
            !
         case default
            stop "photo_angular_parameter(): program stop A."
         end select
         !
         beta = -beta * sqrt(6.0_dp/5.0_dp) / N
         !
      case default
         stop "photo_angular_parameter(): program stop B."
      end select
   end function photo_angular_parameter
   !
   !
   function photo_C_factor(L,ip,i,tline)                     result(C)
   !--------------------------------------------------------------------
   ! Returns the C(L,kappa',kappa) factor for the photoionization
   ! line tline. It uses the explicit expression in terms of the partial
   ! wave photoionization amplitudes as given by Huang, PRA 22 (1980) 223.
   !
   ! Calls: angular_momentum_j(), angular_momentum_l().
   !--------------------------------------------------------------------
      !
      integer, intent(in)          :: L, i, ip
      type(photo_line), intent(in) :: tline
      complex(kind=dp)             :: C
      !
      integer                      :: j, jp
      !
      j  = angular_momentum_j(tline%channel(i)%kappa)
      jp = angular_momentum_j(tline%channel(ip)%kappa)
      !
      C  = tline%channel(ip)%amplitude * conjg(tline%channel(i)%amplitude) * &
           wigner_6j_symbol(tline%channel(i)%totalJ,tline%channel(ip)%totalJ,&
                           L+L,jp,j,tline%totalJ_f) *                        &
           wigner_6j_symbol(tline%channel(i)%totalJ,tline%channel(ip)%totalJ,&
                           L+L,2,2,tline%totalJ_i)  *                        &
           (L+L+one)*sqrt((j+one) *(tline%channel(i)%totalJ+one) *           &
                          (jp+one)*(tline%channel(ip)%totalJ+one))
      !
      if (mod(tline%totalJ_i-tline%totalJ_f+1+32,4) == 0) then
      else if (mod(tline%totalJ_i-tline%totalJ_f+1+32,4) == 2) then
         C = -C
      else
         print *, "tline%totalJ_i, tline%totalJ_f = ", &
                   tline%totalJ_i, tline%totalJ_f
         stop "photo_C_factor(): program stop A."
      end if
      !
   end function photo_C_factor
   !
   !
   subroutine photo_calculate_amplitudes()
   !--------------------------------------------------------------------
   ! Calculates in turn the photoionization amplitudes for all channels and
   ! photoionization lines .
   !
   ! Calls:
   !--------------------------------------------------------------------
      integer       :: i, j, n, nw, nu, nocsf
      real(kind=dp) :: energy
      type(nkappa)  :: subshell
      integer, dimension(:), allocatable :: ndx

      ! Thread-private variables
      integer, dimension(:), allocatable :: local_ndx_f, local_ndx_i
      real(kind=dp), dimension(:,:), allocatable :: local_matrix

      logical           :: first = .true.
      integer           :: ii, ix
      real(kind=dp)     :: ri
      real(kind=dp), dimension(12000), save :: ta_kari

      n = asf_final%csf_set%nocsf + asf_initial%csf_set%nocsf

      ! Allocate global arrays outside parallel region
      allocate( photo_csp%P(1:n_grasp2k), photo_csp%Q(1:n_grasp2k) )
      allocate( cowf_csp%P(1:10), cowf_csp%Q(1:10) )
      allocate( ndx(1:n) )

      ! Parallelize the main loop over lines
      !$OMP PARALLEL PRIVATE(i, j, energy, subshell, local_ndx_f, local_ndx_i, &
      !$OMP                  local_matrix, nw, nocsf, ii, ix, ri) &
      !$OMP          SHARED(ndx, photo_csp, cowf_csp, line, photo, asf_cont, &
      !$OMP                 asf_final, asf_initial, first, ta_kari, n_grasp2k, &
      !$OMP                 r_grasp2k)
      !$OMP DO SCHEDULE(DYNAMIC)
      do i = 1,number_of_lines
         if (line(i)%e_energy < zero) then
            !$OMP CRITICAL
            line(i)%cs_babushkin  = zero
            line(i)%cs_coulomb    = zero
            line(i)%beta_b        = zero
            line(i)%beta_c        = zero
            line(i)%alignment_b   = zero
            line(i)%alignment_c   = zero
            line(i)%xi_b          = zero
            line(i)%xi_c          = zero
            line(i)%eta_b         = zero
            line(i)%eta_c         = zero
            line(i)%zeta_b        = zero
            line(i)%zeta_c        = zero
            !$OMP END CRITICAL
            cycle
         end if

         do j = 1,line(i)%No_channels
            !$OMP CRITICAL(write88)
            write(88,*) "  "
            write(88,2) &
            "Transition: level_i, level_f, channel(kappa,J_tot,parity) = ", &
            line(i)%level_i,line(i)%level_f,line(i)%channel(j)%kappa,       &
            trim(angular_momentum_string(line(i)%channel(j)%totalJ,4)),     &
            line(i)%channel(j)%parity
            !$OMP END CRITICAL

            energy = line(i)%e_energy
            call set_configuration_scheme(asf_final%csf_set,asf_cont%csf_set,&
                     -1,line(i)%channel(j)%kappa,                            &
                     line(i)%totalJ_f,line(i)%parity_f,                      &
                     line(i)%channel(j)%totalJ,line(i)%channel(j)%parity,    &
                     append=.false.,index=ndx)

            ! Allocate thread-private arrays
            photo%No_f = asf_cont%csf_set%nocsf
            allocate( local_ndx_f(photo%No_f) )
            local_ndx_f(1:photo%No_f) = ndx(1:photo%No_f)

            nw = asf_cont%csf_set%nwshells
            if (rabs_use_stop .and. nw /= asf_final%csf_set%nwshells + 1) then
               stop "photo_calculate_amplitudes(): program stop A."
            end if

            nocsf = asf_cont%csf_set%nocsf
            call anco_calculate_csf_matrix(asf_cont%csf_set,1,nocsf,1,nocsf)
            call cowf_set_drs_coefficients(line(i)%asff,asf_cont%csf_set,ndx)
            subshell = nkappa(-1,line(i)%channel(j)%kappa)
            call cowf_set_yk_coefficients(subshell,asf_cont%csf_set)
            call cowf_set_xk_coefficients(subshell,asf_cont%csf_set)

            ! Handle first-time radial grid output
            !$OMP CRITICAL(first_output)
            if (first) then
               first = .false.
               write(47,*) "Radial grid to define the additional potential"
               write(47,*) "----------------------------------------------"
               write(47,*) " "
               do ix = 1,n_grasp2k/3
                  write(47,*) ix, r_grasp2k(ix)
               end do
            end if
            !$OMP END CRITICAL

            ! Set continuum orbital calculation parameters
            cowf_start_homogeneous         = .true.
            cowf_phaseshift_wkb            = .true.
            cowf_phaseshift_zero_potential = .false.
            cowf_phaseshift_coulomb        = .false.
            cowf_norm_nonrel              = .true.
            call cowf_iterate_csp(energy,subshell)

            !$OMP CRITICAL
            photo_csp = cowf_csp
            line(i)%channel(j)%phase = photo_csp%phase
            !$OMP END CRITICAL

            ! Define extended configuration scheme
            call add_csf_to_basis(asf_initial%csf_set,asf_cont%csf_set,      &
                    line(i)%totalJ_i,line(i)%parity_i,index=ndx)
            if (photo_print_csf_scheme) then
               !$OMP CRITICAL(print_csf)
               call print_configuration_scheme(6,asf_cont%csf_set)
               !$OMP END CRITICAL
            end if

            photo%No_i = asf_cont%csf_set%nocsf - photo%No_f
            allocate( local_ndx_i(photo%No_i) )
            local_ndx_i(1:photo%No_i) = ndx(1+photo%No_f:asf_cont%csf_set%nocsf)
            allocate( local_matrix(1:photo%No_f,1:photo%No_i) )

            call photo_pure_matrix(i,j,asf_cont%csf_set)

            ! Move matrix results to shared storage
            !$OMP CRITICAL(photo_matrix)
            allocate( photo%ndx_f(photo%No_f) )
            allocate( photo%ndx_i(photo%No_i) )
            allocate( photo%matrix(1:photo%No_f,1:photo%No_i) )
            photo%ndx_f = local_ndx_f
            photo%ndx_i = local_ndx_i
            photo%matrix = local_matrix
            call photo_channel_amplitude(i,j)
            deallocate( photo%ndx_f, photo%ndx_i, photo%matrix )
            !$OMP END CRITICAL

            deallocate( local_ndx_f, local_ndx_i, local_matrix )
         end do

         ! Calculate line properties
         call photo_line_properties(line(i))
      end do
      !$OMP END DO
      !$OMP END PARALLEL

      deallocate( ndx, photo_csp%P, photo_csp%Q)

    2 format(1x,a,2i5,3x,i4,2x,a4,a)
   end subroutine photo_calculate_amplitudes
   !
   !
   subroutine photo_calculate_Bessel(L,omega_over_c,bessel)
   !--------------------------------------------------------------------
   ! Calculates the Bessel function j_L (w/c*r) over the radial grid
   ! r_grasp2k.
   !
   ! Calls: spherical_Bessel_jL().
   !--------------------------------------------------------------------
      !
      integer, intent(in)       :: L
      real(kind=dp), intent(in) :: omega_over_c
      real(kind=dp), dimension(1:n_grasp2k), intent(out) :: bessel
      !
      integer :: i
      !
      ! Calculate the Bessel function
      bessel(1) = zero
      do  i = 2,n_grasp2k
         bessel(i) = spherical_Bessel_jL(L,omega_over_c*r_grasp2k(i))
      end do
      !
   end subroutine photo_calculate_bessel
   !
   !
   subroutine photo_channel_amplitude(i,j)
   !--------------------------------------------------------------------
   ! Calculates the Photoionization amplitude of channel j of line i
   ! by summing over the 'pure' photoionization matrix using the proper
   ! weights of line i.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: i, j
      !
      integer       :: asfi, asff, l, r, rr, s, ss
      real(kind=dp) :: phase, value
      !
      asfi  = line(i)%asfi;  asff = line(i)%asff
      value = zero
      do  r = 1,photo%no_f
         rr = photo%ndx_f(r)
         do  s = 1,photo%no_i
            ss = photo%ndx_i(s)
	    value = value + asf_final%asf(asff)%eigenvector(rr) * &
	            photo%matrix(r,s) * asf_initial%asf(asfi)%eigenvector(ss)
            !
	 end do
      end do
      !
      l     = angular_momentum_l(line(i)%channel(j)%kappa)
      phase = line(i)%channel(j)%phase
      !
      line(i)%channel(j)%amplitude_re = value
      line(i)%channel(j)%amplitude = cmplx(zero,one)**(-l) *              &
                                exp( cmplx(zero,one)*phase) * cmplx(value,zero)
      !
      !
      print "(a,2i5,2x,2(1pe16.8,1x))", " i,j,line(i)%channel(j)%amplitude = ",&
                                          i,j,line(i)%channel(j)%amplitude
      !
   end subroutine photo_channel_amplitude
   !
   !
   subroutine photo_collect_input()
   !--------------------------------------------------------------------
   ! Collect and proceeds all input for the PHOTO program.
   !
   ! Calls: get_yes_stream().
   !--------------------------------------------------------------------
      !
      integer            :: level_i, level_f, score_position, kappa, ierr
      logical            :: yes
      real(kind=dp)      :: en_lower, en_upper, delta_en, &
                            energy, maximal_energy, wavelength
      character(len=20 ) :: string
      character(len=256) :: record, photo_chn_file
      !
      !
      ! Open, check, load data from, and close the  .iso  file
      call file_get_isodat_grasp2k()
      !
      ! Determine the photoionization multipoles
      call input_transition_multipoles(number_of_multipoles,photo_multipole)
      !
      ! Determine the units for the printout and further optional input
      call input_energy_unit()
      !
      !
   12 print *, "Enter the maximal energy of the photo lines "//&
               " (in " // trim(energy_unit)            //&
               ") to built-up the radial grid;"
      read (*, *, err=12) photo_maximal_energy
      !
      if (energy_inverse) then
         photo_maximal_energy = energy_factor * photo_maximal_energy
      else
         photo_maximal_energy = photo_maximal_energy / energy_factor
      end if
      !
      if (photo_maximal_energy < one) goto 12
      !
      ! Determine grid parameters
      call input_grid_parameters("standard")
      !
      wavelength = sqrt( two * pi * pi / photo_maximal_energy )
      ! Calculations for Kari
      !! hp_grasp2k = wavelength / 50.0_dp
      !! n_grasp2k  = 30.0_dp / hp_grasp2k + 4200
      !! cowf_extent_mtp = 700
      !
      ! Calculations for Nicolai Kabachnik
      hp_grasp2k = wavelength / 177.0_dp
      hp_grasp2k = wavelength / 30.0_dp
      n_grasp2k  = 30.0_dp / hp_grasp2k + 6800
      cowf_extent_mtp = 800
      cowf_extent_mtp = 110
      !
      ! Default speed of light
      c = c_vacuum
      !
      ! Select the photon energies
      maximal_energy = zero
      print *, "The photon energies can be entered either individually or by"
      print *, " an appropriate interval and step size."
      print *, " The default is a list of individual energies; revise this ? "
      yes = get_yes_stream()
      if (.not.yes) then
       4 print *, "Enter another (positive) photon energy in ",      &
		  trim(energy_unit)," enter a negative number if done."
	 read  *, energy
	 if (energy > zero) then
	    if (energy_inverse) then
	       energy = energy_factor * energy
	    else
               energy = energy /energy_factor
            end if
	    number_of_photo_energies = number_of_photo_energies + 1
	    if (number_of_photo_energies > 1000) then
	       print *, "The maximum number of photon energies"// &
	                " is 1000; the program continues ..."
	       goto 6
	    end if
	    photo_energy_selection(number_of_photo_energies) = energy
	    if (energy > maximal_energy)   maximal_energy = energy
	    goto 4
	 end if
      else
       5 print *, "Enter an interval of photon energies and a corresponding"// &
	          " stepsize; "
         print *, " E_lower  E_upper  delta-E:"
         read  *, en_lower, en_upper, delta_en
	 if (en_lower < zero   .or.   (en_upper - en_lower) < zero   .or. &
	     delta_en <= zero) then
	    print *, "All Energies must be greater than zero and "//&
                     " E_upper > E_lower; reenter ..."
	 else
            if (energy_inverse) then
               en_lower = energy_factor * en_lower
               en_upper = energy_factor * en_upper
               delta_en = energy_factor * delta_en
            else
               en_lower = en_lower /energy_factor
               en_upper = en_upper /energy_factor
               delta_en = delta_en /energy_factor
            end if
	    energy = en_lower;   photo_energy_selection(1) = energy
	    number_of_photo_energies = 1
	    do
	       energy = energy + delta_en
	       if (energy > en_upper) exit
	       number_of_photo_energies = number_of_photo_energies + 1
	       if (number_of_photo_energies > 1000) then
	          print *, "The maximum number of photon energies"//&
	                   " is 1000; choose another interval or stepsize."
	          goto 5
	       end if
	       photo_energy_selection(number_of_photo_energies) = energy
	    end do
         endif
      endif
      !
    6 continue
      !
      ! Now 'overwrite' defaults only if required
      print *, "Modify default set-up and printout of the program ?"
      yes = get_yes_stream()
      if (.not.yes) goto 11
      !
      ! Select individual pairs of transitions
      call input_transition_pairs(number_of_atransitions)
      !
      print *, "Include exchange interactions into the generation of the"//&
               " continuum waves ?"
      yes = get_yes_stream()
      cowf_solve_homogeneous_eqn = .not.yes
      !
      photo_nonorthogonal_eval = .false.
      !
      print *, "Calculate angular distribution parameters ?"
      photo_calc_angular_parameter = get_yes_stream()
      !
      print *, "Calculate spin polarization parameters ?"
      photo_calc_spin_polarization = get_yes_stream()
      !
      print *, "Sort all photoionization lines in ascending order of energy ?"
      photo_sort_line_energy = get_yes_stream()
      !
      print *, "Read in and apply experimental energies for the calculation"//&
               " of Auger rates and other properties ?"
      photo_apply_exp_energies = get_yes_stream()
      !
      print *, "Photoionization cross sections are printed in SI units;"// &
               " use Hartree atomic units instead ?"
      photo_print_cs_in_hartree = get_yes_stream()
      !
      print *, "Print all selected photoionization lines"// &
               " before the computation starts (this is the default) ?"
      photo_print_selected_lines = get_yes_stream()
      !
      print *, "Print the CSF scheme each time a new one has been built ?"
      photo_print_csf_scheme =  get_yes_stream()
      !
      print *, "Print the results for each individual line immediatly"//&
               " after its computation ?"
      photo_print_each_line = get_yes_stream()
      !
      print *, "Print the final results to a (.chn) channel amplitude file ?"
      photo_print_chn_file = get_yes_stream()
      if (photo_print_chn_file) then
       8 print *, "Enter a file name for the  photo.chn  file:"
         read *,  photo_chn_file
         call file_open(27,photo_chn_file,"formatted  ","new",ierr)
         !
         if (ierr /= 0) goto 8
      end if
      !
      print *, "Enter an (overall) shift for the atomic transition energies"//&
               " which applies to all transitions (in"//                      &
               trim(energy_unit)//"):"
      print *, " Use  0.  or   <cr>  if no shift need to be applied."
      read (*, *, err=6) photo_energy_shift
      !
      if (energy_inverse) then
         photo_energy_shift = energy_factor * photo_energy_shift
      else
         photo_energy_shift = photo_energy_shift / energy_factor
      end if
      !
    9 print *, "Enter a maximal (-)kappa symmetry up to which continuum"//&
               " spinors are taken into account ?"
      print *, " 2 (up to p-waves), 4(f), 6(i), 8(k), ...;"//&
               " 0 or <cr> to include all possible waves."
      read (*, "(i2)", err=9) kappa
      if (kappa == 0) then
      else if (abs(kappa) < photo_maximal_kappa) then
         photo_maximal_kappa = abs(kappa)
      end if
      !!x write(24,*) "*** photo_maximal_kappa = ",photo_maximal_kappa
      !
      ! Determine the physical effects specifications
   10 print *, "The physical speed of light in atomic units is",c_vacuum,";"
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
   11 continue
      !
      if (photo_print_cs_in_hartree) then
	 photo_cs_factor = one
	 photo_cs_unit   = "a.u.   "
      else
	 photo_cs_factor = bohr_radius_in_cm*bohr_radius_in_cm / 1.0e-18_dp
	 photo_cs_unit   = "Mb     "
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
   end subroutine photo_collect_input
   !
   !
   subroutine photo_initialize()
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
               "photo-ion transitions ..."
      call photo_set_lines()
      call photo_set_overlaps()
      !
      print *, "   ... initialization complete."
      print *, " "
      !
      ! Print the selected transitions before the computation starts
      if (photo_print_selected_lines) then
         call photo_print_lines(6)
         call photo_print_lines(24)
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
   end subroutine photo_initialize
   !
   !
   function photo_kabachnik_B(g,tline,k1,ki,kgamma)                result(B)
   !--------------------------------------------------------------------
   ! Returns the B(k1,ki,k_gamma) parameter due to Kabachnik (2007).
   !
   ! Calls: angular_momentum_j(), angular_momentum_l().
   !--------------------------------------------------------------------
      !
      character(len=1), intent(in) :: g
      type(photo_line), intent(in) :: tline
      integer, intent(in)          :: k1,ki,kgamma
      complex(kind=dp)             :: B
      !
      integer                      :: i, ip, j, jp, totalJ, totalJp, J0, Jf, &
                                      l, lp
      real(kind=dp)                :: phase
      !
      B  = cmplx(zero,zero)
      J0 = tline%totalJ_i
      Jf = tline%totalJ_f
      !
      select case(g)
      case("B")
         do  i = 1,tline%No_channels
            if (tline%channel(i)%multipole == "E1"   .and.    &
        	tline%channel(i)%gauge     == "Babushkin") then
               j      = angular_momentum_j(tline%channel(i)%kappa)
               l      = angular_momentum_l(tline%channel(i)%kappa)
               totalJ = tline%channel(i)%totalJ
               do  ip = 1,tline%No_channels
        	  if (tline%channel(ip)%multipole == "E1"   .and.    &
        	      tline%channel(ip)%gauge	  == "Babushkin") then
        	     jp      = angular_momentum_j(tline%channel(ip)%kappa)
        	     lp      = angular_momentum_l(tline%channel(ip)%kappa)
                     totalJp = tline%channel(ip)%totalJ
		     !
		     if (mod(totalJ+J0+kgamma+kgamma+jp-1+32,4) == 0) then
		        phase = one
		     else if (mod(totalJ+J0+kgamma+kgamma+jp-1+32,4) == 2) then
		        phase = -one
		     else
                        stop "photo_kabachnik_B(): program stop B."
		     end if
		     !
		     B = B + sqrt( (l+l+one) * (lp+lp+one) * (j+one)          &
		                 * (jp+one)  * (totalJ+one) * (totalJp+one) ) &
		       * phase * Clebsch_Gordan(l+l,0,lp+lp,0,k1+k1,0)        &
		       * wigner_6j_symbol(j,l+l,1,lp+lp,jp,k1+k1)	      &
		       * wigner_6j_symbol(2,totalJ,J0,totalJp,2,kgamma+kgamma)&
		       * wigner_9j_symbol(Jf,j,totalJ,Jf,jp,totalJp,          &
		                          ki+ki,k1+k1,kgamma+kgamma)          &
		       * tline%channel(i)%amplitude                           &
		       * conjg(tline%channel(ip)%amplitude)
                  end if
               end do
            end if
         end do
         !
      case("C")
         do  i = 1,tline%No_channels
            if (tline%channel(i)%multipole == "E1"   .and.    &
        	tline%channel(i)%gauge     == "Coulomb  ") then
               j      = angular_momentum_j(tline%channel(i)%kappa)
               l      = angular_momentum_l(tline%channel(i)%kappa)
               totalJ = tline%channel(i)%totalJ
               do  ip = 1,tline%No_channels
        	  if (tline%channel(ip)%multipole == "E1"   .and.    &
        	      tline%channel(ip)%gauge	  == "Coulomb") then
        	     jp      = angular_momentum_j(tline%channel(ip)%kappa)
        	     lp      = angular_momentum_l(tline%channel(ip)%kappa)
                     totalJp = tline%channel(ip)%totalJ
		     !
		     if (mod(totalJ+J0+kgamma+kgamma+jp-1+32,4) == 0) then
		        phase = one
		     else if (mod(totalJ+J0+kgamma+kgamma+jp-1+32,4) == 2) then
		        phase = -one
		     else
                        stop "photo_kabachnik_B(): program stop B."
		     end if
		     !
		     B = B + sqrt( (l+l+one) * (lp+lp+one) * (j+one)          &
		                 * (jp+one)  * (totalJ+one) * (totalJp+one) ) &
		       * phase * Clebsch_Gordan(l+l,0,lp+lp,0,k1+k1,0)        &
		       * wigner_6j_symbol(j,l+l,1,lp+lp,jp,k1+k1)	      &
		       * wigner_6j_symbol(2,totalJ,J0,totalJp,2,kgamma+kgamma)&
		       * wigner_9j_symbol(Jf,j,totalJ,Jf,jp,totalJp,          &
		                          ki+ki,k1+k1,kgamma+kgamma)          &
		       * tline%channel(i)%amplitude                           &
		       * conjg(tline%channel(ip)%amplitude)
                  end if
               end do
            end if
         end do
         !
      case default
         stop "photo_kabachnik_B(): program stop B."
      end select
      !
      B = sqrt( (ki+ki+one) * (k1+k1+one) ) * B
      !
   end function photo_kabachnik_B
   !
   !
   subroutine photo_kabachnik_5k_append_table(stream,tline)
   !--------------------------------------------------------------------
   ! Append a table of Kabachnik's B-5k coefficients to the results
   ! of the PHOTO component on stream.
   !--------------------------------------------------------------------
      !
      integer, intent(in)          :: stream
      type(photo_line), intent(in) :: tline
      !
      integer          :: j0, J1, k, k0, k1, k2, k_gamma
      complex(kind=dp) :: B
      !
      write(stream,*) "The B-5k coefficients are given in the form: "
      write(stream,*) " "
      write(stream,*) "  2*J0  2*J1   k1  k2     k0   k  k_gamma    B (Re, Im) "
      write(stream,*) " "
      write(stream,*) " "
      !
      J0 = tline%totalJ_i
      J1 = tline%totalJ_f
      !
      do  k0 = 0,2
         do  k_gamma = 0,2
            do  k2 = 0,J1
	       !
               do  k = 0,8
	          if (.not.is_triangle(k0+k0, k_gamma+k_gamma, k+k)) cycle
                  do  k1 = 0,8
	             if (.not.is_triangle(k2+k2, k1+k1, k+k)) cycle
		     !
		     B = photo_kabachnik_5k_coefficient(k1,k2,k0,k, &
		                                              k_gamma,tline)
		     !
		     if (abs(B) > eps10*eps10) then
		        write(stream,10) J0, J1, k1, k2, k0, k, k_gamma, B
	             end if
                  end do
               end do
               !
            end do
         end do
      end do
      !
   10 format(2i4,3x,2i4,2x,3i4,3x,1pe12.5,1x,1pe12.5)
      !
   end subroutine photo_kabachnik_5k_append_table
   !
   !
   function photo_kabachnik_5k_coefficient(k1,k2,k0,k, k_gamma,tline)   &
                                                                result(B)
   !--------------------------------------------------------------------
   ! Returns the B-5k coefficient as defined by Nicolai Kabachnik,
   ! April 2006 for the collaboration with Fabiana (Italy).
   !
   ! Calls: angular_momentum_j(), angular_momentum_l().
   !--------------------------------------------------------------------
      !
      integer, intent(in)          :: k1, k2, k0, k, k_gamma
      type(photo_line), intent(in) :: tline
      complex(kind=dp)             :: B
      !
      integer                      :: i, ip, J0, J1, j, jp, l, lp, &
                                      totalJ, totalJp
      real(kind=dp)                :: alpha, factor
      !
      J0 = tline%totalJ_i
      J1 = tline%totalJ_f
      !
      ! Calculate the `factor' between the internal PI amplitudes and
      ! those amplitudes D_ljJ as used by NMK
      alpha  = one / one_over_alpha
      factor = two * pi * alpha * tline%p_energy * sqrt(two * alpha)
      !
      B = zero
      !
      do  i = 1,tline%No_channels
         do  ip = 1,tline%No_channels
            j       = angular_momentum_j(tline%channel(i)%kappa)
            jp      = angular_momentum_j(tline%channel(ip)%kappa)
            l       = angular_momentum_l(tline%channel(i)%kappa)
            lp      = angular_momentum_l(tline%channel(ip)%kappa)
            totalJ  = tline%channel(i)%totalJ
            totalJp = tline%channel(ip)%totalJ
	    !
	    B = B + sqrt( (totalJ+one)*(totalJp+one)*(l+l+one)*(lp+lp+one) &
	                 *(j+one)*(jp+one) )                               &
		    * ( (-1)**((jp+1)/2) )                                 &
		    * Clebsch_Gordan(l+l,0,lp+lp,0,k1+k1,0)                &
		    * wigner_6j_symbol(l+l,j,1,jp,lp+lp,k1+k1)             &
		    * wigner_9j_symbol(J0,2,totalJ,J0,2,totalJp,           &
		                       k0+k0,k_gamma+k_gamma,k+k,.false.)  &
		    * wigner_9j_symbol(J1,j,totalJ,J1,jp,totalJp,          &
		                       k2+k2,k1+k1,k+k,.false.)            &
	            * tline%channel(i)%amplitude                           &
	            * conjg(tline%channel(ip)%amplitude)                   &
		    / (factor * factor)
	    !
         end do
      end do
      !
      if (mod(k1+k2-k_gamma,2) == 0) then
         B = sqrt(k+k+one) * B
      else
         B = - sqrt(k+k+one) * B
      end if
      !
   end function photo_kabachnik_5k_coefficient
   !
   !
   subroutine photo_line_properties(tline)
   !--------------------------------------------------------------------
   ! Calculates all selected photoionization properties of line i from
   ! the amplitudes of the individual photoionization channels. This
   ! concerns the total cross sections as well as angular distribution and
   ! spin polarization parameters.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      type(photo_line), intent(inout) :: tline
      !
      integer       :: j, stream
      real(kind=dp) :: energy, p_energy, e_energy, cs, csc, csb, sumb, sumc
      real(kind=dp) :: wb, wc
      !
      sumc = zero;   sumb = zero
      do  j = 1,tline%No_channels
         if (tline%channel(j)%gauge == "Babushkin"  .or. &
             tline%channel(j)%gauge == "Magnetic ")   then
            sumb = sumb + tline%channel(j)%amplitude *      &
                    conjg(tline%channel(j)%amplitude)
         end if
         if (tline%channel(j)%gauge == "Coulomb  "  .or. &
             tline%channel(j)%gauge == "Magnetic ")   then
            sumc = sumc + tline%channel(j)%amplitude *      &
                    conjg(tline%channel(j)%amplitude)
         end if
      end do
      tline%cs_babushkin = half*c*c/tline%p_energy * sumb &
                           / (three*(tline%totalJ_i + one))
      tline%cs_coulomb   = half*c*c/tline%p_energy * sumc &
                           / (three*(tline%totalJ_i + one))
      !
      ! Calculate the angular distribution and spin-polarization coefficients
      if (photo_calc_angular_parameter) then
         tline%beta_b      = photo_angular_parameter("B","beta",tline)
         tline%beta_c      = photo_angular_parameter("C","beta",tline)
         wb                = photo_alignment_parameter("B","beta     ",tline)
         wc                = photo_alignment_parameter("C","beta     ",tline)
         tline%alignment_b = photo_alignment_parameter("B","alignment",tline)
         tline%alignment_c = photo_alignment_parameter("C","alignment",tline)
      end if
      !
      if (photo_calc_spin_polarization) then
         tline%xi_b   = photo_spin_parameter("B","xi  ",tline)
         tline%xi_c   = photo_spin_parameter("C","xi  ",tline)
         tline%eta_b  = photo_spin_parameter("B","eta ",tline)
         tline%eta_c  = photo_spin_parameter("C","eta ",tline)
         tline%zeta_b = photo_spin_parameter("B","zeta",tline)
         tline%zeta_c = photo_spin_parameter("C","zeta",tline)
      end if
      !
      ! Print a short summary of the line if required
      !
      if (photo_print_each_line) then
         stream = 6
      13 write(stream,*) " "
         write(stream,3) " Results for PI line from the transition ",    &
                         tline%level_i," - ",tline%level_f,": ",         &
                         trim(angular_momentum_string(tline%totalJ_i)),  &
                         " ",tline%parity_i,"   ----> ",                 &
                         trim(angular_momentum_string(tline%totalJ_f)),  &
                         " ",tline%parity_f
         write(stream,*) " --------------------------------------------"//  &
                         "------------------------------------"
         write(stream,*) " "
         if (energy_inverse) then
            p_energy = energy_factor / tline%p_energy
            e_energy = energy_factor / tline%e_energy
         else
            p_energy = energy_factor * tline%p_energy
            e_energy = energy_factor * tline%e_energy
         end if
         csb = tline%cs_babushkin * photo_cs_factor
         csc = tline%cs_coulomb   * photo_cs_factor
         write(stream,4) "   Photon energy       = ",  &
                         p_energy,trim(energy_unit)
         write(stream,4) "   Electron energy     = ",  &
                         e_energy,trim(energy_unit)
         write(stream,4) "   Total cross section = ",csb,trim(photo_cs_unit),&
                         "  (Babushkin gauge) "
         write(stream,4) "                       = ",csc,trim(photo_cs_unit),&
                         "  (Coulomb gauge) "
         !
         if (photo_calc_angular_parameter) then
            write(stream,5) "   Beta                = ",real(tline%beta_b),  &
                            ";   ",real(tline%beta_c),"  (Babushkin; Coulomb) "
            write(stream,5) "   Alignment (final)   = ",real(tline%alignment_b),&
                            ";   ",real(tline%alignment_c),                  &
			    "  (Babushkin; Coulomb) "
         end if
         !
         if (photo_calc_spin_polarization) then
            write(stream,5) "   Delta               = ",                    &
                            real(tline%zeta_b-two*tline%xi_b)/three,";   ", &
                            real(tline%zeta_c-two*tline%xi_c)/three,        &
                            "  (Babushkin; Coulomb) "
            write(stream,5) "   Xi                  = ",real(tline%xi_b),   &
                            ";   ",real(tline%xi_c)
            write(stream,5) "   Eta                 = ",real(tline%eta_b),  &
                            ";   ",real(tline%eta_c)
            write(stream,5) "   Zeta                = ",real(tline%zeta_b), &
                            ";   ",real(tline%zeta_c)
         end if
         write(stream,*) " "
         write(stream,1) trim(photo_cs_unit)
         write(stream,*) &
              " -----------------------------------------------------------"//&
              "---------------------------------------------------"
         do  j = 1,tline%No_channels
            cs = tline%channel(j)%amplitude *                                &
                          conjg(tline%channel(j)%amplitude)* photo_cs_factor
            write(stream,2) orbital_symmetry(tline%channel(j)%kappa),        &
                    trim(angular_momentum_string(tline%channel(j)%totalJ,4)),&
                            tline%channel(j)%parity,                         &
                            tline%channel(j)%multipole,                      &
                            tline%channel(j)%gauge(1:9),                     &
                            tline%channel(j)%amplitude,                      &
			    tline%channel(j)%amplitude_re,cs,                &
                            mod(tline%channel(j)%phase,two*pi)
         end do
         write(stream,*) " "
       1 format("  Kappa    Total J^P  Mp     Gauge            Amplitude ",  &
                "       Real-Amplitude    Rate (",a,")        Phase")
       2 format(3x,a4,5xa4,a1,5x,a2,3x,a9,3x,1pe10.3,2x,1pe10.3,4x,1pe10.3,4x,&
                1pe12.5,4x,1pe10.3)
       3 format(1x,a,i5,a,i5,a,a,a,a,a,a,a,a,a)
       4 format(a,1pe16.8,2x,a,a)
       5 format(a,1pe16.8,2x,a,1pe16.8,2x,a)
         !
         ! Re-cycle once on stream 24
         if (stream == 6) then
            stream = 24;   goto 13
         end if
      end if
      !
   end subroutine photo_line_properties
   !
   !
   subroutine photo_print_amplitudes(stream)
   !--------------------------------------------------------------------
   ! Prints the information about all (selected) transitions and amplitudes
   ! to a (.chn) channel amplitude file on stream.
   !
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      !
      integer             ::  i, j, channels
      real(kind=dp)       :: wa
      !
      write(stream,*) &
         "Transition data and amplitudes for photoionization lines are printed "
      write(stream,*) &
         "in the format (one line per photoionization channel): "
      write(stream,*) " "
      write(stream,*) "   ASF_i, ASF_f, Level_i, Level_f, 2*J_i, 2*J_f, ..."
      write(stream,*) "     ... P_i, P_f, p_energy, e_energy, kappa, ...   "
      write(stream,*) "     ... 2*J_t, P_t, multipole, gauge, phase, amplitude"
      write(stream,*) &
         "================================================================"//&
         "================================================================"
      write(stream,*) " "
      !
      channels = 0
      do  i = 1,number_of_lines
         channels = channels + line(i)%No_channels
      end do
      write(stream,*) channels, "= Number_of_channels"
      write(stream,*) " "
      !
      do  i = 1,number_of_lines
         do  j = 1,line(i)%No_channels
	    wa = abs(line(i)%channel(j)%amplitude)
	    wa = wa*wa
            write(stream,1) line(i)%asfi,line(i)%asff,               &
                            line(i)%level_i,line(i)%level_f,         &
                            line(i)%totalJ_i,line(i)%totalJ_f,       &
                            line(i)%parity_i,line(i)%parity_f,       &
                            line(i)%p_energy,line(i)%e_energy,       &
                            line(i)%channel(j)%kappa,                &
                            line(i)%channel(j)%totalJ,               &
                            line(i)%channel(j)%parity,               &
                            line(i)%channel(j)%multipole,            &
                            line(i)%channel(j)%gauge,                &
                            line(i)%channel(j)%phase,                &
                            line(i)%channel(j)%amplitude,            &
			    wa,line(i)%p_energy*27.21_dp
         end do
      end do
      !
      1 format(2i5,3x,2i5,3x,2i3,2x,2a3,1x,e14.7,1x,e14.7,i4,i5,a3,2x, &
               a2,2x,a9,2x,e14.7,2x,e14.7,1x,e14.7,        4x,e14.7,1x,e14.7)
      !
   end subroutine photo_print_amplitudes
   !
   !
   subroutine photo_print_lines(stream)
   !--------------------------------------------------------------------
   ! Prints a neat table of all selected 'photon' or 'electron' lines
   ! on stream before the actual computation starts; only the quantum
   ! numbers of the atomic states and the transition energies are
   ! displayed.
   !
   ! Calls: angular_momentum_string(), get_kappa_channels().
   !--------------------------------------------------------------------
      !
      integer, intent(in)    :: stream
      integer                :: i, j, m, nchannel, mchannel
      real(kind=dp)          :: e_energy, p_energy
      integer, dimension(20) :: kappa_channel
      !
      write(stream,1) number_of_lines,trim(energy_unit)
      do  i = 1,number_of_lines
         if (energy_inverse) then
            e_energy = energy_factor / line(i)%e_energy
            p_energy = energy_factor / line(i)%p_energy
         else
            e_energy = energy_factor * line(i)%e_energy
            p_energy = energy_factor * line(i)%p_energy
         end if
         write(stream,2) line(i)%level_i,line(i)%level_f,            &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)), &
                  line(i)%parity_i,                                  &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)), &
                  line(i)%parity_f,p_energy,e_energy,                &
                  (orbital_symmetry(line(i)%channel(j)%kappa),       &
                   line(i)%channel(j)%multipole,                     &
                   line(i)%channel(j)%gauge(1:1),                    &
                   trim(angular_momentum_string(line(i)%channel(j)%totalJ,4)),&
                   line(i)%channel(j)%parity,j=1,line(i)%No_channels)
      end do
      write(stream,3)
    1 format(/,"The following ",i5," lines are selected:",                 &
        //,"     I-level-F     I--J^P--F        Photon Energy    Electron",&
           " energy      Photoionization channels",                        &
         /,"                                     (in ",a4,")  ",           &
         /,4x,133("-") )
    2 format(4x,i4," -",i4,3x,a4,a1,3x,a4,a1,5x,1pe14.7,5x,1pe14.7,6x,     &
                   4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x), &
             /,74x,4(a2,"(",a2,";",a1,") J=",a,a1,3x) )
    3 format(4x,133("-") )
      !
   end subroutine photo_print_lines
   !
   !
   subroutine photo_print_results(stream)
   !--------------------------------------------------------------------
   ! Writes the photoionization cross sections and angular distribution
   ! parameters, and further properties in a neat summary to the .sum file.
   !
   ! Calls:  angular_momentum_string().
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      !
      integer       :: i
      real(kind=dp) :: p_energy, e_energy, csc, csb, csb_tot, csc_tot, &
                       ratio_c, ratio_b, total_B, total_C
      !
      write(stream,1)
    1 format(/ &
         /30x,"==========================================================================", &
         /30x,"|  Summary of all Photoionization Cross sections and Angular Parameters  |", &
         /30x,"==========================================================================", &
             // )
    2 format( / 1x,128("-"),                                                            &
              / 2x,"LevI-LevF   I- J / Parity -F      ",                                &
                   "  omega           e-Energy     CS-Babushkin   ",                    &
                   "  CS-Coulomb   Total Babushkin   Total Coulomb    ",                &
              / 2x,"                                  ",                                &
                   "  (",a,")             (",a,")           (a.u.)    ",                &
                   "      (a.u.)          (a.u.)          (a.u.) "   ,                  &
              / 1x,128("-") )
    3 format( / 1x,128("-"),                                                            &
              / 2x,"LevI-LevF   I- J / Parity -F      ",                                &
                   "  omega           e-Energy     CS-Babushkin   ",                    &
                   "  CS-Coulomb   Ratio Babushkin   Ratio Coulomb    ",                &
              / 2x,"                                  ",                                &
                   "  (",a,")             (",a,")            (",a,")     ",             &
                   "       (",a,")            (",a,")            (",a,") "   ,          &
              / 1x,128("-") )
    4 format(   2x,i3,' -',i3,3x,a4,a2,4x,a4,a2,5x,2(1pe12.5,5x),                       &
                   4(1pe10.3,6x),     8x,1pe10.3,2x,1pe10.3 )
    6 format( / 1x,105("-"),                                                             &
              / 2x,"LevI-LevF   I- J / Parity -F      ",                                &
                   "omega (",a,")     e-Energy (",a,")       Beta        ", &
                   " Beta         Beta   ",                                 &
              /68x,"(Babushkin)   (Coulomb)    (Average) "                                         &
              / 1x,105("-") )
    7 format( / 1x,105("-"),                                                             &
              / 2x,"LevI-LevF   I- J / Parity -F      ",                                &
                   "omega (",a,")     e-Energy (",a,")     Alignment    ",  &
                   "Alignment    Alignment ",                               &
              /68x,"(Babushkin)   (Coulomb)    (Average) "                                         &
              / 1x,105("-") )
    8 format(   2x,i3,' -',i3,3x,a4,a2,4x,a4,a2,5x,2(1pe12.5,5x),                       &
                   8(1pe10.3,3x) )
    9 format(   1x,128("-") )
   10 format( / 1x,118("-"),                                                             &
              / 2x,"LevI-LevF   I- J / Parity -F      ",                                &
                   "omega (",a,")     e-Energy (",a,")       Delta          Xi          Eta         Zeta   ",    &
              / 1x,118("-") )
   11 format(   1x,105("-") )
   12 format(   1x,118("-") )
      !
      ! Print the photoionization cross sections of this calculation
      !
      write(stream,*) "Individual photoionization cross sections and ratio:"
      write(stream,*) "-----------------------------------------------------"
      if (photo_print_cs_in_hartree) then
         write(stream,2) trim(energy_unit), trim(energy_unit)
      else
         write(stream,3) trim(energy_unit),  trim(energy_unit),            &
                         trim(photo_cs_unit),trim(photo_cs_unit),          &
                         trim(photo_cs_unit),trim(photo_cs_unit)
      end if
      !
      csc_tot = zero;   csb_tot = zero
      do  i = 1,number_of_lines
         csb_tot = csb_tot + line(i)%cs_babushkin
         csc_tot = csc_tot + line(i)%cs_coulomb
      end do
      !
       ratio_b = zero;  ratio_c = zero;  total_B = zero; total_C=zero;
      do  i = 1,number_of_lines
         if (energy_inverse) then
            p_energy = energy_factor / line(i)%p_energy
            e_energy = energy_factor / line(i)%e_energy
         else
            p_energy = energy_factor * line(i)%p_energy
            e_energy = energy_factor * line(i)%e_energy
         end if
         csb      = line(i)%cs_babushkin * photo_cs_factor
         csc      = line(i)%cs_coulomb   * photo_cs_factor
         !
         ratio_b = csb/(csb_tot * photo_cs_factor)
         ratio_c = csc/(csc_tot * photo_cs_factor)
         total_B = total_B + ratio_b
         total_C = total_C + ratio_c
       !
       !  photo_print_only_gt_1percent                          &
         if(abs(ratio_b) > photo_print_cut) then
           write(stream,4) line(i)%level_i, line(i)%level_f,   &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)),  &
                            line(i)%parity_i,                         &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)),  &
                            line(i)%parity_f,                         &
                            p_energy,e_energy,csb,csc,ratio_b,ratio_c
         end if
      end do
      write(stream,9)
      write(stream,*) "Selected total B-cross section ratio= ", total_B
      write(stream,*) "Total Babushkin-cross section (Mb)= ", csb_tot*photo_cs_factor
      write(stream,*) "Selected total C-cross section ratio= ", total_C
      write(stream,*) "Total Coulomb-cross section (Mb)= ", csc_tot*photo_cs_factor

      !
      ! Print photoionization angular parameters
      !

      if (photo_calc_angular_parameter) then
         write(stream,*) " "
         write(stream,*) "Photoionization angular parameters:"
         write(stream,*) "-----------------------------------"
         write(stream,6) trim(energy_unit), trim(energy_unit)
         do  i = 1,number_of_lines
            if (energy_inverse) then
               p_energy = energy_factor / line(i)%p_energy
               e_energy = energy_factor / line(i)%e_energy
            else
               p_energy = energy_factor * line(i)%p_energy
               e_energy = energy_factor * line(i)%e_energy
            end if
            !
            write(stream,8) line(i)%level_i, line(i)%level_f,               &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)),        &
                            line(i)%parity_i,                               &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)),        &
                            line(i)%parity_f,                               &
                            p_energy,e_energy,                              &
                            real(line(i)%beta_b),real(line(i)%beta_c),      &
                            real(line(i)%beta_b+line(i)%beta_c)/two
         end do
         write(stream,11)
         write(stream,*) " "
         write(stream,*) " "
         write(stream,7) trim(energy_unit), trim(energy_unit)
         do  i = 1,number_of_lines
            if (energy_inverse) then
               p_energy = energy_factor / line(i)%p_energy
               e_energy = energy_factor / line(i)%e_energy
            else
               p_energy = energy_factor * line(i)%p_energy
               e_energy = energy_factor * line(i)%e_energy
            end if
            !
            write(stream,8) line(i)%level_i, line(i)%level_f,               &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)),        &
                            line(i)%parity_i,                               &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)),        &
                      line(i)%parity_f, 			            &
                      p_energy,e_energy,			            &
                      real(line(i)%alignment_b),real(line(i)%alignment_c),	&
                      real(line(i)%alignment_b+line(i)%alignment_c)/two
         end do
         write(stream,11)
         write(stream,*) " "
      end if
      !
      !
      ! Print photoionization spin-polarization parameters
      !
      if (photo_calc_spin_polarization) then
         write(stream,*) " "
         write(stream,*) "Spin-polarization parameters (Babushkin gauge):"
         write(stream,*) "-----------------------------------------------"
         write(stream,10) trim(energy_unit), trim(energy_unit)
         do  i = 1,number_of_lines
            if (energy_inverse) then
               p_energy = energy_factor / line(i)%p_energy
               e_energy = energy_factor / line(i)%e_energy
            else
               p_energy = energy_factor * line(i)%p_energy
               e_energy = energy_factor * line(i)%e_energy
            end if
            !
            write(stream,8) line(i)%level_i, line(i)%level_f,             &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)),      &
                            line(i)%parity_i,                             &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)),      &
                            line(i)%parity_f,                             &
                            p_energy,e_energy,                            &
                            real(line(i)%zeta_b-two*line(i)%xi_b)/three,  &
                            real(line(i)%xi_b),real(line(i)%eta_b),       &
                            real(line(i)%zeta_b)
         end do
         write(stream,12)
         !
         write(stream,*) " "
         write(stream,*) "Spin-polarization parameters (Coulomb gauge):"
         write(stream,*) "---------------------------------------------"
         write(stream,10) trim(energy_unit), trim(energy_unit)
         do  i = 1,number_of_lines
            if (energy_inverse) then
               p_energy = energy_factor / line(i)%p_energy
               e_energy = energy_factor / line(i)%e_energy
            else
               p_energy = energy_factor * line(i)%p_energy
               e_energy = energy_factor * line(i)%e_energy
            end if
            !
            write(stream,8) line(i)%level_i, line(i)%level_f,             &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)),      &
                            line(i)%parity_i,                             &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)),      &
                            line(i)%parity_f,                             &
                            p_energy,e_energy,                            &
                            real(line(i)%zeta_c-two*line(i)%xi_c)/three,  &
                            real(line(i)%xi_c),real(line(i)%eta_c),       &
                            real(line(i)%zeta_c)
         end do
         write(stream,12)
         !
         write(stream,*) " "
         write(stream,*) "Spin-polarization parameters (average):"
         write(stream,*) "---------------------------------------"
         write(stream,10) trim(energy_unit), trim(energy_unit)
         do  i = 1,number_of_lines
            if (energy_inverse) then
               p_energy = energy_factor / line(i)%p_energy
               e_energy = energy_factor / line(i)%e_energy
            else
               p_energy = energy_factor * line(i)%p_energy
               e_energy = energy_factor * line(i)%e_energy
            end if
            !
            write(stream,8) line(i)%level_i, line(i)%level_f,             &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)),      &
                            line(i)%parity_i,                             &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)),      &
                            line(i)%parity_f,                             &
                            p_energy,e_energy,                            &
                            real(line(i)%zeta_b-two*line(i)%xi_b          &
                          +      line(i)%zeta_c-two*line(i)%xi_c)/6.0_dp, &
                            real(line(i)%xi_b + line(i)%xi_c)/two,        &
                            real(line(i)%eta_b + line(i)%eta_c)/two,      &
                            real(line(i)%zeta_b + line(i)%zeta_c)/two
         end do
         write(stream,12)
         !
      end if
      !
      ! Carry out and print some 'non-standard' additional quantities which are
      ! not part of the PHOTO program
      if (photo_kabachnik_5k_coeff  .and.  stream == 6) then
         do  i = 1,number_of_lines
            call photo_kabachnik_5k_append_table(73,line(i))
	 end do
      end if
      !
   end subroutine photo_print_results
   !
   !
   subroutine photo_pure_matrix(tr,ch,csf_set)
   !--------------------------------------------------------------------
   ! Calculates the 'pure' PHOTO matrix for the given configuration scheme
   ! csf_set. The first no_f CSF belong to the final-state representation
   ! and the following no_f+1,...,no_f+no_i to the initial states.
   ! The procedure takes into account ...
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in)	  :: tr, ch
      type(csf_basis), intent(in) :: csf_set
      !
      integer			  :: i, ia, ib, j, ja, jb, no_T_coeff, r, s, &
				     ss, t, nu, par, parity, ii
      type(nkappa)		  :: aa, bb, cc, dd
      real(kind=dp)		  :: aweight, wy
      type(nkappa)                :: Ta, Tb
      !
      select case(line(tr)%channel(ch)%multipole(1:2))
      case("E1");   nu = 1;   par = 1;   parity = -1
      case("M1");   nu = 1;   par = 0;   parity =  1
      case("E2");   nu = 2;   par = 1;   parity =  1
      case("M2");   nu = 2;   par = 0;   parity = -1
      case("E3");   nu = 3;   par = 1;   parity = -1
      case("M3");   nu = 3;   par = 0;   parity =  1
      case("E4");   nu = 4;   par = 1;   parity =  1
      case("M4");   nu = 4;   par = 0;   parity = -1
      case default
         stop "photo_pure_matrix(): program stop A."
      end select
      !
      if (rabs_use_stop   .and.   csf_set%nocsf /= photo%no_f+photo%no_i) then
         stop "photo_pure_matrix(): program stop B."
      end if
      !
      photo%matrix = zero
      !
      if (photo_nonorthogonal_eval) then
         call nonorth_initialize_ci(csf_set,wave_final,wave_initial,photo_csp)
      end if
      !
      !$OMP PARALLEL DO PRIVATE(r, s, ss, t) SCHEDULE(STATIC)
      do  r = 1,photo%No_f
         do  s = photo%No_f+1,photo%No_f+photo%No_i
            ss = s - photo%No_f
	    !
	    if (photo_nonorthogonal_eval) then
	       !
	       ! Use an expansion into Slater determinants to evaluate the
	       ! photoionization matrix
	       Aoperator%particle = 1
	       Aoperator%rank     = nu
	       Aoperator%parity   = parity
	       !
               call nonorth_calculate_csf_pair(csf_set,Aoperator,r,s,no_T_coeff)
	    else
	       !
	       ! Use standard Racah algebra to evaluate the many-electron
	       ! photoionization matrix
               !! call anco_calculate_csf_pair_1p(nu,csf_set,r,s,no_T_coeff)
	       !
	       number_of_mct_coefficients = 0
               call mct_generate_coefficients(r,r,s,s,csf_set,parity,nu)
	       no_T_coeff = number_of_mct_coefficients
	       !!x print *, "*** r, s, number_of_mct_coefficients = ",   &
	       !!x 	         r, s, number_of_mct_coefficients
	    end if
            !
            ! Cycle over all angular coefficients
            do  t = 1,no_T_coeff
	       if (photo_nonorthogonal_eval) then
	          Ta      = nonorth_T_list(t)%a
	          Tb      = nonorth_T_list(t)%b
                  ja      = angular_momentum_j(Ta%kappa)
                  aweight = nonorth_T_list(t)%T * sqrt( ja+one )
		  !
                  ia = -1;   ib = 0
		  if (Ta == photo_csp%orbital) then
		     ia = csf_set%nwshells
		     ia = 0
		  else
                     do  i = 1,wave_final%number_of_rwf
                        if (wave_final%rwf(i)%orbital == Ta) then
                           ia = i
                           exit
                        end if
                     end do
		  end if
		  !
                  do  i = 1,wave_initial%number_of_rwf
                     if (wave_initial%rwf(i)%orbital == Tb) then
                        ib = i
                        exit
                     end if
                  end do
		  !
		  !!x print *, " "
		  !!x print *, "NonO: ia, ib, Ta, Tb, aweight = ", &
		  !!x 		  ia, ib, Ta, Tb, aweight
	       else
	          ia       = mct_list(t)%a
	          ib       = mct_list(t)%b
	          aweight  = mct_list(t)%T
		  if (ia  /= csf_set%nwshells) then
		     Ta    = wave_final%rwf(ia)%orbital
		  else
		     Ta    = photo_csp%orbital
		     ia    = 0
		  end if
		  Tb       = wave_initial%rwf(ib)%orbital
		  !!x print *, "Orth: ia, ib, Ta, Tb, aweight = ", &
		  !!x 		      ia, ib, Ta, Tb, aweight
	       end if
	       !
               if (ia == 0  .and.  ib /= 0) then
                  jb = angular_momentum_j(Tb%kappa)
                  ja = angular_momentum_j(Ta%kappa)
		  !
                  photo%matrix(r,ss) = photo%matrix(r,ss) + aweight *         &
                  photo_reduced_M_integral(tr,line(tr)%channel(ch)%multipole, &
                  		    line(tr)%channel(ch)%gauge, 	      &
                  		    photo_csp,wave_initial%rwf(ib))
		  !!x    rk_integral_grasp2k_cd(photo_csp,wave_initial%rwf(ib),1)
               else
                  jb = angular_momentum_j(Tb%kappa)
                  ja = angular_momentum_j(Ta%kappa)
		  !
		  ! Introduced on 10/7/09 while working on xamplitude
		  print *, "photo_pure_matrix() stop xx"
		  !
                  photo%matrix(r,ss) = photo%matrix(r,ss) + aweight *         &
                  photo_reduced_M_integral(tr,line(tr)%channel(ch)%multipole, &
                  		    line(tr)%channel(ch)%gauge, 	      &
                  		    wave_final%rwf(ia),wave_initial%rwf(ib))
               end if
               !
            end do
         end do
      end do
      !$OMP END PARALLEL DO
      !
      if (photo_nonorthogonal_eval) then
         call nonorth_finalize("ci")
      end if
      !
      !
   end subroutine photo_pure_matrix
   !
   !
   function photo_reduced_M_integral(tr,multipole,gauge,rwf_f,rwf_i) &
                                                          result(red_me)
   !--------------------------------------------------------------------
   ! Calculates the value of an (one-electron) M integral for the coupling
   ! of electrons with the radiation field for given multipole and
   ! gauge. It applies the set of Bessel functions which are associated
   ! with the line tr; the radial wave functions are provided by the
   ! derived data structures rwf_f and rwf_i.
   !
   ! Calls: angular_momentum_j(), orbital_name(), quad_grasp2k(),
   !        wigner_3j_symbol().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: tr
      character(len=2), intent(in)       :: multipole
      character(len=9), intent(in)       :: gauge
      type(orbital_function), intent(in) :: rwf_f, rwf_i
      real(kind=dp)                      :: red_me
      !
      integer       :: ja, jb, kapa, kapb, L, mtp, phase
      real(kind=dp) :: factor, red_Coulomb, red_Gauge, wa, wb
      real(kind=dp), dimension(:), allocatable :: ta
      !
      red_me = zero
      !
      if (.false.) then
         write (99,1) orbital_name(rwf_f%orbital%n,rwf_f%orbital%kappa), &
                      orbital_name(rwf_i%orbital%n,rwf_i%orbital%kappa)
       1 format(/1x,30("+"),"   photo_reduced_M_integral() called for ", &
                "orbitals", a4, " and ", a4, ":",3x,30("+"))
      end if
      !
      kapa = rwf_f%orbital%kappa;        kapb = rwf_i%orbital%kappa
      ja   = angular_momentum_j(kapa);   jb   = angular_momentum_j(kapb)
      select case(multipole)
      case("E1", "M1");   L = 1
      case("E2", "M2");   L = 2
      case("E3", "M3");   L = 3
      case("E4", "M4");   L = 4
      case default
         stop "photo_reduced_M_integral(): program stop A."
      end select
      !
      ! Test parity rules
      select case(multipole)
      case("M1", "M2", "M3", "M4")
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 0)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 2) ) then
         else
            return
         end if
      case("E1", "E2", "E3", "E4")
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 2)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 0) ) then
         else
            return
         end if
      end select
      !
      ! Evaluate factor multiplying Mbar(a,b); use one-particle matrix elements
      ! of Brink-Satchler type
      phase  = (ja + 1)/2 + jb
      factor = ((-one)**phase) *sqrt(jb+one)*wigner_3j_symbol(ja,L+L,jb,1,0,-1)
      if (abs(factor) < eps10) then
         red_me = zero
         return
      end if
      !
      mtp = min( rwf_f%mtp, rwf_i%mtp)
      allocate( ta(1:mtp+10) )
      !
      select case(multipole)
      case("M1", "M2", "M3", "M4")
         ta(1) = zero
         ta(2:mtp) = ( rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +                   &
                       rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) ) * rp_grasp2k(2:mtp)
         select case(multipole)
         case("M1")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel1(2:mtp)
         case("M2")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel2(2:mtp)
         case("M3")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel3(2:mtp)
         case("M4")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel4(2:mtp)
         end select
         red_me = -(L+L+one)*(kapa+kapb)/sqrt(L*(L+one)) * quad_grasp2k(ta,mtp)
         red_me = factor * red_me
      case("E1", "E2", "E3", "E4")
         ta(1) = zero
         !
         wa    = sqrt(L/(L+one))*(kapa-kapb);  wb = sqrt(L/(L+one))*(L+one)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         select case(multipole)
         case("E1")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel2(2:mtp)
         case("E2")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel3(2:mtp)
         case("E3")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel4(2:mtp)
         case("E4")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel5(2:mtp)
         end select
         red_Coulomb = factor * quad_grasp2k(ta,mtp)
         !
         wa    = - sqrt((L+one)/L)*(kapa-kapb);  wb = sqrt((L+one)/L)*L
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         select case(multipole)
         case("E1")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel0(2:mtp)
         case("E2")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel1(2:mtp)
         case("E3")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel2(2:mtp)
         case("E4")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel3(2:mtp)
         end select
         red_Coulomb = red_Coulomb + factor * quad_grasp2k(ta,mtp)
         !
         if (gauge == "Coulomb  ") then
            red_me = red_Coulomb
            if (.false.) &
               write(99,*) "factor, red_Coulomb = ",factor, red_me
            deallocate( ta )
            return
         else if (gauge == "Babushkin") then
         else if (rabs_use_stop) then
            stop "photo_reduced_M_integral(): program stop B."
         end if
         !
         wa    = (kapa-kapb);  wb = (L+1)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         select case(multipole)
         case("E1")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel2(2:mtp)
         case("E2")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel3(2:mtp)
         case("E3")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel4(2:mtp)
         case("E4")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel5(2:mtp)
         end select
         red_Gauge = factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         !
         wa    = (kapa-kapb);  wb = -L
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         select case(multipole)
         case("E1")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel0(2:mtp)
         case("E2")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel1(2:mtp)
         case("E3")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel2(2:mtp)
         case("E4")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel3(2:mtp)
         end select
         red_Gauge = red_Gauge + factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         !
	 ! Sign due to K.G. Dyall et al., CPC (1989), used original and
	 ! up to august 2007
	 !
         wa    = -(L+L+one)
	 !
	 ! New sign due to I.P. Grant, JPB (1974), advised by Sami and used
	 ! since august 2007
         !! wa    = (L+L+one)
	 !
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%P(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%Q(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         select case(multipole)
         case("E1")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel1(2:mtp)
         case("E2")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel2(2:mtp)
         case("E3")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel3(2:mtp)
         case("E4")
            ta(2:mtp) = ta(2:mtp) * line(tr)%bessel4(2:mtp)
         end select
	 !
         red_Gauge = red_Gauge + factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         !
	 ! Gauge test due to Sami; September 2007
	 !
	 !! red_Gauge = red_Gauge + factor * quad_grasp2k(ta,mtp)
         red_me    = red_Coulomb + red_Gauge
         !
         !!x red_me    = red_Coulomb + red_Gauge / ten
         !
      case default
         stop "photo_reduced_M_integral(): program stop C."
      end select
      !
      deallocate( ta )
      !
      if (.true.) then
         write (88,2) orbital_name(rwf_f%orbital%n,rwf_f%orbital%kappa),     &
                      orbital_name(rwf_i%orbital%n,rwf_i%orbital%kappa),red_me
       2 format(1x,3("+")," Reduced_M_integral() for orbitals  ", &
                    a4, " and ", a4, ":",3x,1pe10.3)
      end if
      !
   end function photo_reduced_M_integral
   !
   !
   subroutine photo_set_lines()
   !--------------------------------------------------------------------
   ! Determines how many and which photoionization lines need to be
   ! calculated and initializes the array line of type(photo_lines).
   ! The default is (for number_of_lines == 0) that all lines
   ! with a positive electron energy are calculated; for number_of_lines /= 0,
   ! lines for individual atomic transitions were selected during input time
   ! and will be initialized instead; in this case also negative line energies
   ! are allowed and may be overritten by appropriate experimental energies.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer       :: asfi, asff, i, imin, j, k, l, m, nchannels, nt, nch, &
                       totalJ
      real(kind=dp) :: arg, energy, e_energy, p_energy, energy_exp, energy_min
      logical       :: yes0, yes1, yes2, yes3, yes4, yes5
      integer, dimension(20)             :: kappa_channel
      character(len=1) :: parity
      !
      logical, dimension(:), allocatable    :: ask_for
      integer, dimension(200)               :: wa_kappa, wa_totalJ
      character(len=1), dimension(200)      :: wa_parity
      character(len=2), dimension(200)      :: wa_multipole
      character(len=9), dimension(200)      :: wa_gauge
      real(kind=dp), dimension(1:n_grasp2k) :: bessel
      !
      ! Determine the total number of lines
      nt = 0
      if (number_of_atransitions == 0) then
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_photo_energies
                  if (asf_initial%asf(i)%energy - asf_final%asf(j)%energy &
                                      + photo_energy_selection(k) > zero) then
                     nt = nt + 1
                  end if
               end do
            end do
         end do
      else
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_atransitions
                  if ( (select_level_i(k) == asf_initial%asf(i)%level_No .or. &
                        select_level_i(k) == 0)                         .and. &
                       (select_level_f(k) == asf_final%asf(j)%level_No   .or. &
                        select_level_f(k) == 0) ) then
                     do  l = 1,number_of_photo_energies
                        nt = nt + 1
                     end do
                  end if
               end do
            end do
         end do
      end if
      !
      ! Allocate the derived data structure line
      allocate( line(1:nt) )
      !
      ! Now initialize all lines
      nt = 0
      if (number_of_atransitions == 0) then
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_photo_energies
                  if (asf_initial%asf(i)%energy - asf_final%asf(j)%energy &
                                      + photo_energy_selection(k) > zero) then
                     nt = nt + 1
                     line(nt)%asfi     = i;   line(nt)%asff = j
                     line(nt)%level_i  = asf_initial%asf(i)%level_No
                     line(nt)%level_f  = asf_final%asf(j)%level_No
                     line(nt)%p_energy = photo_energy_selection(k)
                     line(nt)%e_energy = photo_energy_selection(k) +         &
                           asf_initial%asf(i)%energy - asf_final%asf(j)%energy
                  end if
               end do
            end do
         end do
      else
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_atransitions
                  if ( (select_level_i(k) == asf_initial%asf(i)%level_No .or. &
                        select_level_i(k) == 0)                         .and. &
                       (select_level_f(k) == asf_final%asf(j)%level_No   .or. &
                        select_level_f(k) == 0) ) then
                     do  l = 1,number_of_photo_energies
                        nt = nt + 1
                        line(nt)%asfi     = i;   line(nt)%asff = j
                        line(nt)%level_i  = asf_initial%asf(i)%level_No
                        line(nt)%level_f  = asf_final%asf(j)%level_No
                        line(nt)%p_energy = photo_energy_selection(l)
                        line(nt)%e_energy = photo_energy_selection(l) +      &
                           asf_initial%asf(i)%energy - asf_final%asf(j)%energy
                     end do
                  end if
               end do
            end do
         end do
      end if
      !
      number_of_lines = nt
      !
      ! Apply experimental energies if requested
      if (photo_apply_exp_energies) then
         allocate( ask_for(1:nt) )
         ask_for = .false.
         !
         do k = 1,nt
            if (.not.ask_for(k)) then
               i = line(k)%asfi;   j = line(k)%asff
               energy  = asf_initial%asf(i)%energy - asf_final%asf(j)%energy
               if (energy_inverse) then
                  energy = energy_factor / energy
               else
                  energy = energy_factor * energy
               end if
             1 print *, "Atomic transition",line(k)%level_i,"-",                 &
                line(k)%level_f,"has ab-initio energy E_theo = ",&
                energy,trim(energy_unit),";"
               print *, " enter E_exp (in ",trim(energy_unit), &
                ") or 0.0 to use E_theo:"
               read(*,*,err=1) energy_exp
               if (energy_exp /= 0) then
                  if (energy_inverse) then
                     energy = energy_factor/energy_exp
                  else
                     energy = energy_exp/energy_factor
                  end if
               end if
               !
               do l = k,nt
                  if (i == line(l)%asfi   .and.   j == line(l)%asff) then
                     line(l)%e_energy = line(l)%p_energy + energy
                     ask_for(l) = .true.
                  end if
               end do
            end if
         end do
         deallocate( ask_for )
      end if
      !
      !
      ! Order the photoionization lines in ascending order of electron energies
      if (photo_sort_line_energy) then
         do  i = 1,number_of_lines-1
            imin       = i
            energy_min = 1.0e20
            do  j = i,number_of_lines
               if (line(j)%e_energy < energy_min) then
                  imin = j;   energy_min = line(j)%e_energy
               end if
            end do
            if (imin /= i) then
               asfi     = line(i)%asfi
               asff     = line(i)%asff
               e_energy = line(i)%e_energy
               p_energy = line(i)%p_energy
               line(i)%asfi        = line(imin)%asfi
               line(i)%asff        = line(imin)%asff
               line(i)%e_energy    = line(imin)%e_energy
               line(i)%p_energy    = line(imin)%p_energy
               line(imin)%asfi     = asfi
               line(imin)%asff     = asff
               line(imin)%e_energy = e_energy
               line(imin)%p_energy = p_energy
            end if
         end do
      end if
      !
      ! Now assign the full line specifications and calculate the
      ! necessary Bessel functions
      do  nt = 1,number_of_lines
         i = line(nt)%asfi;   j = line(nt)%asff
         line(nt)%level_i     = asf_initial%asf(i)%level_No
         line(nt)%level_f     = asf_final%asf(j)%level_No
         line(nt)%totalJ_i    = asf_initial%asf(i)%totalJ
         line(nt)%totalJ_f    = asf_final%asf(j)%totalJ
         line(nt)%parity_i    = asf_initial%asf(i)%parity
         line(nt)%parity_f    = asf_final%asf(j)%parity
         line(nt)%beta_b      = zero
         line(nt)%beta_c      = zero
         line(nt)%alignment_b = zero
         line(nt)%alignment_c = zero
         line(nt)%xi_b        = zero
         line(nt)%xi_c        = zero
         line(nt)%eta_b       = zero
         line(nt)%eta_c       = zero
         line(nt)%zeta_b      = zero
         line(nt)%zeta_c      = zero
         !
         !
         ! Modify the transition energies if required
         line(nt)%e_energy = line(nt)%e_energy + photo_energy_shift
         !
         nch = 0
         do k = 1,number_of_multipoles
            select case(photo_multipole(k))
            case ("E1")
               if (line(nt)%parity_i == "+") then
                  parity = "-"
               else
                  parity = "+"
               end if
               !
               do  totalJ = line(nt)%totalJ_i-2,line(nt)%totalJ_i+2,2
                  if (totalJ < 0) cycle
                  !! if (totalJ < 0  .or.                                 &
		  !!     .not.is_triangle(line(nt)%totalJ_i,2,totalJ)) cycle
                  call get_kappa_channels(totalJ,parity,               &
                                 line(nt)%totalJ_f,line(nt)%parity_f,  &
                                 kappa_channel,nchannels)
                  if (rabs_use_stop   .and.   nchannels == 0) then
                     stop "photo_set_lines(): program stop A."
                  end if
                  !
                  ! Determine the number of channels
                  do  m = 1,nchannels
		     if (line(nt)%totalJ_i == 0  .and.  totalJ == 0) cycle
                     if (abs(kappa_channel(m)) <= photo_maximal_kappa) then
                        nch = nch + 1
                        wa_kappa(nch)     = kappa_channel(m)
                        wa_totalJ(nch)    = totalJ
                        wa_parity(nch)    = parity
                        wa_multipole(nch) = "E1"
                        wa_gauge(nch)     = "Babushkin"
                        nch = nch + 1
                        wa_kappa(nch)     = kappa_channel(m)
                        wa_totalJ(nch)    = totalJ
                        wa_parity(nch)    = parity
                        wa_multipole(nch) = "E1"
                        wa_gauge(nch)     = "Coulomb  "
                     end if
                  end do
               end do
            case ("M1")
               if (line(nt)%parity_i == "+") then
                  parity = "+"
               else
                  parity = "-"
               end if
               !
               do  totalJ = line(nt)%totalJ_i-2,line(nt)%totalJ_i+2,2
                  if (totalJ < 0  .or.                                 &
		      .not.is_triangle(line(nt)%totalJ_i,2,totalJ)) cycle
                  call get_kappa_channels(totalJ,parity,               &
                                 line(nt)%totalJ_f,line(nt)%parity_f,  &
                                 kappa_channel,nchannels)
                  if (rabs_use_stop   .and.   nchannels == 0) then
                     stop "photo_set_lines(): program stop C."
                  end if
                  !
                  ! Determine the number of channels
                  do  m = 1,nchannels
		     if (line(nt)%totalJ_i == 0  .and.  totalJ == 0) cycle
                     if (abs(kappa_channel(m)) <= photo_maximal_kappa) then
                        nch = nch + 1
                        wa_kappa(nch)     = kappa_channel(m)
                        wa_totalJ(nch)    = totalJ
                        wa_parity(nch)    = parity
                        wa_multipole(nch) = "M1"
                        wa_gauge(nch)     = "Magnetic "
                     end if
                  end do
               end do
            case ("E2")
               if (line(nt)%parity_i == "+") then
                  parity = "+"
               else
                  parity = "-"
               end if
               !
               do  totalJ = line(nt)%totalJ_i-4,line(nt)%totalJ_i+4,2
                  if (totalJ < 0  .or.                                 &
		      .not.is_triangle(line(nt)%totalJ_i,4,totalJ)) cycle
                  call get_kappa_channels(totalJ,parity,               &
                                 line(nt)%totalJ_f,line(nt)%parity_f,  &
                                 kappa_channel,nchannels)
                  if (rabs_use_stop   .and.   nchannels == 0) then
                     stop "photo_set_lines(): program stop D."
                  end if
                  !
                  ! Determine the number of channels
                  do  m = 1,nchannels
		     if (line(nt)%totalJ_i == 0  .and.  totalJ == 0) cycle
                     if (abs(kappa_channel(m)) <= photo_maximal_kappa) then
                        nch = nch + 1
                        wa_kappa(nch)     = kappa_channel(m)
                        wa_totalJ(nch)    = totalJ
                        wa_parity(nch)    = parity
                        wa_multipole(nch) = "E2"
                        wa_gauge(nch)     = "Babushkin"
                        nch = nch + 1
                        wa_kappa(nch)     = kappa_channel(m)
                        wa_totalJ(nch)    = totalJ
                        wa_parity(nch)    = parity
                        wa_multipole(nch) = "E2"
                        wa_gauge(nch)     = "Coulomb  "
                     end if
                  end do
               end do
            case ("M2")
               if (line(nt)%parity_i == "+") then
                  parity = "-"
               else
                  parity = "+"
               end if
               !
               do  totalJ = line(nt)%totalJ_i-4,line(nt)%totalJ_i+4,2
                  if (totalJ < 0  .or.                                 &
		      .not.is_triangle(line(nt)%totalJ_i,4,totalJ)) cycle
                  call get_kappa_channels(totalJ,parity,               &
                                 line(nt)%totalJ_f,line(nt)%parity_f,  &
                                 kappa_channel,nchannels)
                  if (rabs_use_stop   .and.   nchannels == 0) then
                     stop "photo_set_lines(): program stop E."
                  end if
                  !
                  ! Determine the number of channels
                  do  m = 1,nchannels
		     if (line(nt)%totalJ_i == 0  .and.  totalJ == 0) cycle
                     if (abs(kappa_channel(m)) <= photo_maximal_kappa) then
                        nch = nch + 1
                        wa_kappa(nch)     = kappa_channel(m)
                        wa_totalJ(nch)    = totalJ
                        wa_parity(nch)    = parity
                        wa_multipole(nch) = "M2"
                        wa_gauge(nch)     = "Magnetic "
                     end if
                  end do
               end do
            case ("E3")
               if (line(nt)%parity_i == "+") then
                  parity = "-"
               else
                  parity = "+"
               end if
               !
               do  totalJ = line(nt)%totalJ_i-6,line(nt)%totalJ_i+6,2
                  if (totalJ < 0) cycle
                  call get_kappa_channels(totalJ,parity,               &
                                 line(nt)%totalJ_f,line(nt)%parity_f,  &
                                 kappa_channel,nchannels)
                  if (rabs_use_stop   .and.   nchannels == 0) then
                     stop "photo_set_lines(): program stop F."
                  end if
                  !
                  ! Determine the number of channels
                  do  m = 1,nchannels
		     if (line(nt)%totalJ_i == 0  .and.  totalJ == 0) cycle
                     if (abs(kappa_channel(m)) <= photo_maximal_kappa) then
                        nch = nch + 1
                        wa_kappa(nch)     = kappa_channel(m)
                        wa_totalJ(nch)    = totalJ
                        wa_parity(nch)    = parity
                        wa_multipole(nch) = "E3"
                        wa_gauge(nch)     = "Babushkin"
                        nch = nch + 1
                        wa_kappa(nch)     = kappa_channel(m)
                        wa_totalJ(nch)    = totalJ
                        wa_parity(nch)    = parity
                        wa_multipole(nch) = "E3"
                        wa_gauge(nch)     = "Coulomb  "
                     end if
                  end do
               end do
            case default
               stop "photo_set_lines(): program stop G."
            end select
         end do
         !
         line(nt)%No_channels = nch
         allocate( line(nt)%channel(1:nch) )
         !
         line(nt)%channel(1:nch)%kappa      = wa_kappa(1:nch)
         line(nt)%channel(1:nch)%totalJ     = wa_totalJ(1:nch)
         line(nt)%channel(1:nch)%parity     = wa_parity(1:nch)
         line(nt)%channel(1:nch)%multipole  = wa_multipole(1:nch)
         line(nt)%channel(1:nch)%gauge      = wa_gauge(1:nch)
         line(nt)%channel(1:nch)%phase      = zero
         line(nt)%channel(1:nch)%amplitude  = cmplx(zero,zero)
         !
         ! Calculate the Bessel function for the given transition
         ! energy; first allocate the array
         yes0 = .false.;   yes1 = .false.    ! yes == Needs to be
         yes2 = .false.;   yes3 = .false.    !        calculated ?
         yes4 = .false.;   yes5 = .false.
         do  m = 1,nch
            select case(wa_multipole(m))
            case("E1", "M1"); yes0 = .true.; yes1 = .true.; yes2 = .true.
            case("E2", "M2"); yes1 = .true.; yes2 = .true.; yes3 = .true.
            case("E3", "M3"); yes2 = .true.; yes3 = .true.; yes4 = .true.
            case("E4", "M4"); yes3 = .true.; yes4 = .true.; yes5 = .true.
            case default
               stop "photo_set_lines(): program stop B."
            end select
         end do
         !
         arg  = line(nt)%p_energy / c
         if (yes0) then
            allocate( line(nt)%bessel0(1:n_grasp2k) )
            call photo_calculate_Bessel(0,arg,bessel)
            line(nt)%bessel0 = bessel
         end if
         if (yes1) then
            allocate( line(nt)%bessel1(1:n_grasp2k) )
            call photo_calculate_Bessel(1,arg,bessel)
            line(nt)%bessel1 = bessel
         end if
         if (yes2) then
            allocate( line(nt)%bessel2(1:n_grasp2k) )
            call photo_calculate_Bessel(2,arg,bessel)
            line(nt)%bessel2 = bessel
         end if
         if (yes3) then
            allocate( line(nt)%bessel3(1:n_grasp2k) )
            call photo_calculate_Bessel(3,arg,bessel)
            line(nt)%bessel3 = bessel
         end if
         if (yes4) then
            allocate( line(nt)%bessel4(1:n_grasp2k) )
            call photo_calculate_Bessel(4,arg,bessel)
            line(nt)%bessel4 = bessel
         end if
         if (yes5) then
            allocate( line(nt)%bessel5(1:n_grasp2k) )
            call photo_calculate_Bessel(5,arg,bessel)
            line(nt)%bessel5 = bessel
         end if
      end do
      !
      print *, "  ",number_of_lines,                         &
               " lines have been initialized and will be "// &
               "calculated in this run of the program."
      !
   end subroutine photo_set_lines
   !
   !
   subroutine photo_set_overlaps()
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
      allocate( photo_overlap(kappa_min:kappa_max,pqn_f_max,pqn_i_max) )
      !
      ! Calculate the non-orthogonal and overlap integrals
      do  i1 = 1,wave_final%number_of_rwf
         do  i2 = 1,wave_initial%number_of_rwf
            if (wave_final%rwf(i1)%orbital%kappa == &
                wave_initial%rwf(i2)%orbital%kappa) then
               kappa = wave_final%rwf(i1)%orbital%kappa
               pqn_f = wave_final%rwf(i1)%orbital%n
               pqn_i = wave_initial%rwf(i2)%orbital%n
               photo_overlap(kappa,pqn_f,pqn_i) = &
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
   end subroutine photo_set_overlaps
   !
   !
   function photo_spin_parameter(g,keyword,tline)             result(wa)
   !--------------------------------------------------------------------
   ! Returns the spin angular parameters for the photoionization
   ! line tline. It uses the explicit expression in terms of the partial
   ! wave photoionization amplitudes as given by Huang, PRA 22 (1980) 223.
   !
   ! Calls: angular_momentum_j(), angular_momentum_l().
   !--------------------------------------------------------------------
      !
      character(len=1), intent(in) :: g
      character(len=4), intent(in) :: keyword
      type(photo_line)             :: tline
      complex(kind=dp)             :: wa
      !
      integer                      :: i, ip, j, jp, lp
      complex(kind=dp)             :: N
      !
      wa = cmplx(zero,zero)
      !
      N = cmplx(zero,zero)
      !
      select case(g)
      case("B")
         do  i = 1,tline%No_channels
            if (tline%channel(i)%multipole == "E1"   .and.    &
                tline%channel(i)%gauge     == "Babushkin") then
               N = N + tline%channel(i)%amplitude *  &
                 conjg(tline%channel(i)%amplitude)
            end if
         end do
      case("C")
         do  i = 1,tline%No_channels
            if (tline%channel(i)%multipole == "E1"   .and.    &
                tline%channel(i)%gauge     == "Coulomb  ") then
               N = N + tline%channel(i)%amplitude *  &
                 conjg(tline%channel(i)%amplitude)
            end if
         end do
      case default
         stop "photo_spin_parameter(): program stop A."
      end select
      !
      select case (keyword)
      case ("xi  ")
         !
         select case(g)
         case("B")
            do  i = 1,tline%No_channels
               if (tline%channel(i)%multipole == "E1"   .and.    &
                   tline%channel(i)%gauge     == "Babushkin") then
                  j = angular_momentum_j(tline%channel(i)%kappa)
                  do  ip = 1,tline%No_channels
                     if (tline%channel(ip)%multipole == "E1"   .and.    &
                         tline%channel(ip)%gauge     == "Babushkin") then
                        jp = angular_momentum_j(tline%channel(ip)%kappa)
                        lp = angular_momentum_l(tline%channel(ip)%kappa)
                        if (mod(jp+lp+lp-1+32,4) == 0) then
                           wa = wa + wigner_3j_symbol(j,jp,2,1,1,-2) * &
                                photo_C_factor(1,ip,i,tline)
                        else if (mod(jp+lp+lp-1+32,4) == 2) then
                           wa = wa - wigner_3j_symbol(j,jp,2,1,1,-2) * &
                                photo_C_factor(1,ip,i,tline)
                        else
                           stop "photo_spin_parameter(): program stop B."
                        end if
                     end if
                  end do
               end if
            end do
            !
         case("C")
            do  i = 1,tline%No_channels
               if (tline%channel(i)%multipole == "E1"   .and.    &
                   tline%channel(i)%gauge     == "Coulomb  ") then
                  j = angular_momentum_j(tline%channel(i)%kappa)
                  do  ip = 1,tline%No_channels
                     if (tline%channel(ip)%multipole == "E1"   .and.    &
                         tline%channel(ip)%gauge     == "Coulomb  ") then
                        jp = angular_momentum_j(tline%channel(ip)%kappa)
                        lp = angular_momentum_l(tline%channel(ip)%kappa)
                        if (mod(jp+lp+lp-1+32,4) == 0) then
                           wa = wa + wigner_3j_symbol(j,jp,2,1,1,-2) * &
                                photo_C_factor(1,ip,i,tline)
                        else if (mod(jp+lp+lp-1+32,4) == 2) then
                           wa = wa - wigner_3j_symbol(j,jp,2,1,1,-2) * &
                                photo_C_factor(1,ip,i,tline)
                        else
                           stop "photo_spin_parameter(): program stop C."
                        end if
                     end if
                  end do
               end if
            end do
         end select
         !
         wa = wa * sqrt(three) / (two*N)
         !
      case ("eta ")
         !
         select case(g)
         case("B")
            do  i = 1,tline%No_channels
               if (tline%channel(i)%multipole == "E1"   .and.    &
                   tline%channel(i)%gauge     == "Babushkin") then
                  j = angular_momentum_j(tline%channel(i)%kappa)
                  do  ip = 1,tline%No_channels
                     if (tline%channel(ip)%multipole == "E1"   .and.    &
                         tline%channel(ip)%gauge     == "Babushkin") then
                        jp = angular_momentum_j(tline%channel(ip)%kappa)
                        lp = angular_momentum_l(tline%channel(ip)%kappa)
                        if (mod(jp+lp+lp+1+32,4) == 0) then
                           wa = wa + wigner_3j_symbol(j,jp,4,1,1,-2) * &
                                photo_C_factor(2,ip,i,tline)
                        else if (mod(jp+lp+lp+1+32,4) == 2) then
                           wa = wa - wigner_3j_symbol(j,jp,4,1,1,-2) * &
                                photo_C_factor(2,ip,i,tline)
                        else
                           stop "photo_spin_parameter(): program stop D."
                        end if
                     end if
                  end do
               end if
            end do
            !
         case("C")
            do  i = 1,tline%No_channels
               if (tline%channel(i)%multipole == "E1"   .and.    &
                   tline%channel(i)%gauge     == "Coulomb  ") then
                  j = angular_momentum_j(tline%channel(i)%kappa)
                  do  ip = 1,tline%No_channels
                     if (tline%channel(ip)%multipole == "E1"   .and.    &
                         tline%channel(ip)%gauge     == "Coulomb  ") then
                        jp = angular_momentum_j(tline%channel(ip)%kappa)
                        lp = angular_momentum_l(tline%channel(ip)%kappa)
                        if (mod(jp+lp+lp+1+32,4) == 0) then
                           wa = wa + wigner_3j_symbol(j,jp,4,1,1,-2) * &
                                photo_C_factor(2,ip,i,tline)
                        else if (mod(jp+lp+lp+1+32,4) == 2) then
                           wa = wa - wigner_3j_symbol(j,jp,4,1,1,-2) * &
                                photo_C_factor(2,ip,i,tline)
                        else
                           stop "photo_spin_parameter(): program stop E."
                        end if
                     end if
                  end do
               end if
            end do
         end select
         !
         wa = wa * three / (sqrt(5.0_dp)*two*N) * cmplx(zero,one)
         !
      case ("zeta")
         !
         select case(g)
         case("B")
            do  i = 1,tline%No_channels
               if (tline%channel(i)%multipole == "E1"   .and.    &
                   tline%channel(i)%gauge     == "Babushkin") then
                  j = angular_momentum_j(tline%channel(i)%kappa)
                  do  ip = 1,tline%No_channels
                     if (tline%channel(ip)%multipole == "E1"   .and.    &
                         tline%channel(ip)%gauge     == "Babushkin") then
                        jp = angular_momentum_j(tline%channel(ip)%kappa)
                        wa = wa + wigner_3j_symbol(j,jp,2,1,-1,0) * &
                                photo_C_factor(2,ip,i,tline)
                     end if
                  end do
               end if
            end do
            !
         case("C")
            do  i = 1,tline%No_channels
               if (tline%channel(i)%multipole == "E1"   .and.    &
                   tline%channel(i)%gauge     == "Coulomb  ") then
                  j = angular_momentum_j(tline%channel(i)%kappa)
                  do  ip = 1,tline%No_channels
                     if (tline%channel(ip)%multipole == "E1"   .and.    &
                         tline%channel(ip)%gauge     == "Coulomb  ") then
                        jp = angular_momentum_j(tline%channel(ip)%kappa)
                        wa = wa + wigner_3j_symbol(j,jp,2,1,-1,0) * &
                                photo_C_factor(2,ip,i,tline)
                     end if
                  end do
               end if
            end do
         end select
         !
         wa = wa * sqrt(three/two) / N
         !
      case default
         stop "photo_angular_parameter(): program stop F."
      end select
   end function photo_spin_parameter
   !
end module rabs_photo
