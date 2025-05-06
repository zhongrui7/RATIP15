module rabs_rec
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains the procedures which are specific to the REC program. 
! This program supports the calculation of radiative electron capture cross 
! sections and angular distribution parameters. 
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
   use rabs_constant
   use rabs_cowf
   use rabs_csl
   use rabs_determinant
   use rabs_dirac_orbital
   use rabs_file_handling
   use rabs_grasp2k
   use rabs_input_dialog
   use rabs_mcp
   use rabs_multipole
   use rabs_nonorthonormal
   use rabs_print
   use rabs_rcfp
   use rabs_recoupling
   implicit none
   private
   !
   !
   private :: rec_angular_parameter
                 ! Calculates the REC angular parameter beta^RR_nu(J_i,J_f)
                 ! from the many-electron REC amplitudes
   public  :: rec_calculate_amplitudes
                 ! Calculates for all selected 'electron energies' (in turn) the 
                 ! required continuum spinors and REC amplitudes.
   private :: rec_calculate_Bessel
                 ! Calculates the Bessel function over the radial grid
                 ! for a given factor.
   private :: rec_channel_amplitude
                 ! Calculates the REC amplitude of a channel from the 'pure' 
                 ! REC matrix and the corresponding mixing coefficients.
   public  :: rec_collect_input
                 ! Collects and proceeds all input for the calculation of the
		 ! REC capture cross sections and angular distribution 
                 ! parameters.
   public  :: rec_initialize
                 ! Set up the selected transitions and initializes some
		 ! arrays as required.
   private :: rec_line_properties
                 ! Calculates all selected REC properties of line i 
		 ! from the amplitudes of the individual channels.
   public  :: rec_print_amplitudes
                 ! Writes out the information about all (selected)
                 ! REC channels and amplitudes to a (.chn) file
   private :: rec_print_lines
                 ! Prints all selected 'photon' or 'electron' lines in a neat 
                 ! format before the computation starts.
   public ::  rec_print_results
                 ! Prints the REC capture cross sections and angular 
		 ! distribution parameters as well as further properties in a 
		 ! neat summary.
   private :: rec_pure_matrix
                 ! Calculates the 'pure' REC matrix for the given
		 ! configuration scheme.
   private :: rec_reduced_M_integral_old
                 ! Calculates the (one-electron) M integral.
   private :: rec_set_lines
                 ! Determines which ionization lines need to be calculated
                 ! and initializes storage for an appropriate data type for
                 ! them.
   private :: rec_set_overlaps
                 ! Initializes the array of overlap integrals for the
                 ! calculation of relaxed REC cross sections.
   private :: rec_statistical_tensor
                 ! Calculates the REC statistical tensors rho^RR_k(J_i,J_f)
                 ! from the many-electron REC amplitudes
   !
   ! Define some global data of the REC program; most of these data are
   ! read in during the interactive control and may overwrite existing
   ! default values
   !
   ! Storage for the initial and final atomic states and wave functions
   !
   ! Define a 'configuration scheme' to be built up during execution 
   type(grasp2k_orbital), public :: wave_continuum
   type(orbital_function)        :: rec_csp
   !
   ! Define an internal structure type(rec_line) which stores all
   ! necessary information for a REC capture line
   type :: rec_channel
      integer          :: kappa, totalJ, multipole_L, multipole_pi   
      character(len=1) :: parity 
      character(len=9) :: gauge 
      real(kind=dp)    :: phase 
      complex(kind=dp) :: amplitude   
   end type rec_channel
   !
   type :: rec_line
      integer          :: asfi, asff, level_i, level_f, totalJ_i, totalJ_f
      integer          :: No_channels
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: p_energy, e_energy, cs_coulomb, cs_babushkin
      complex(kind=dp),  dimension(40)         :: beta_b, beta_c
      complex(kind=dp),  dimension(0:6)        :: rho_b,  rho_c
      type(rec_channel), dimension(:), pointer :: channel
   end type rec_line
   !
   type(rec_line), dimension(:), allocatable :: line
   !
   type :: rec_matrix
      integer :: No_f, No_i  
      integer, dimension(:), pointer         :: ndx_f, ndx_i  
      real(kind=dp), dimension(:,:), pointer :: matrix   
   end type rec_matrix
   !
   type(rec_matrix) :: rec
   !
   integer          :: number_of_lines                       = 0,  &
                       number_of_electron_energies           = 0,  &
                       number_of_ctransitions                = 0
   !
   integer, public  :: rec_maximal_kappa                     = 0,  &
                       rec_maximal_nu                        = 40, & 
                       rec_maximal_tensor                    = 6,  &
                       rec_anisotropy_J0                     = 0
   !
   ! Define global logical flags for the control of the REC program; the default
   ! values for these flags may be overwritten interactively during input time
   logical, public ::  rec_add_1                     = .false.,  &
                       rec_calc_angular_parameter    = .true.,   &
                       rec_calc_statistical_tensor   = .true.,   &
                       rec_calc_anisotropy_parameter = .true.,   &
		       rec_nonorthogonal_eval        = .false.,  &
		       rec_print_csf_scheme          = .false.,  &
		       rec_print_each_line           = .true.,   &
                       rec_print_only_gt_1percent    = .false.,  &     
		       rec_print_cs_in_hartree       = .true.,   &
                       rec_print_selected_lines      = .true.,   &
                       rec_print_chn_file            = .false.,  &
                       rec_sort_line_energy          = .true.
   !
   ! Energy unit for the output of all energies
   real(kind=dp)    :: rec_cs_factor = zero, &
                       rec_energy_shift  = zero, rec_maximal_energy
   character(len=7) :: rec_cs_unit
   !
   ! Define storage for the overlap integrals: (kappa,pqn_f,pqn_i) 
   real(kind=dp), dimension(:,:,:), allocatable :: rec_overlap
   !
   ! Define some variables and arrays for processing input data from 
   ! rec_collect_input()
   real(kind=dp), dimension(100)   :: rec_energy_selection
   !
contains
   !
   !
   function rec_angular_parameter(G,nu,tline)               result(beta)
   !--------------------------------------------------------------------
   ! Calculates the REC angular parameter beta_nu (alpha_i J_i,alpha_f J_f) 
   ! from the many-electron REC amplitudes.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      character(len=1), intent(in) :: G
      integer, intent(in)          :: nu
      type(rec_line), intent(in)   :: tline
      complex(kind=dp)             :: beta
      !
      integer          :: i, ip, J, Jp, jk, jkp, J_i, J_f,  &
                          L, Lp, lk, lkp, p, pp
      complex(kind=dp) :: D, N, me, mep
      !
      J_i = tline%totalJ_i
      J_f = tline%totalJ_f
      !
      D   = zero;   N = zero
      !  
      do  i = 1,tline%No_channels
         if (tline%channel(i)%gauge(1:1) == G(1:1)   .or. &
             tline%channel(i)%gauge(1:1) == "M") then
            me = tline%channel(i)%amplitude
         else
            cycle
         end if
         !
         D = D + me * conjg(me)
         !
         J   = tline%channel(i)%totalJ
         L   = tline%channel(i)%multipole_L
         p   = tline%channel(i)%multipole_pi
         lk  = angular_momentum_l(tline%channel(i)%kappa)   
         jk  = angular_momentum_j(tline%channel(i)%kappa) 
         !  
         do  ip = 1,tline%No_channels
            if (tline%channel(ip)%gauge(1:1) == G(1:1)   .or. &
                tline%channel(ip)%gauge(1:1) == "M") then
               mep = tline%channel(ip)%amplitude
            else
               cycle
            end if
            !
            Jp  = tline%channel(ip)%totalJ
            Lp  = tline%channel(ip)%multipole_L  
            pp  = tline%channel(ip)%multipole_pi
            lkp = angular_momentum_l(tline%channel(ip)%kappa)   
            jkp = angular_momentum_j(tline%channel(ip)%kappa)
            !
            N  = N + cmplx(zero,one)**(L+p-Lp-pp+16) * (-1)**((J_i-1-J_f+16)/2) &
               * sqrt( (L+L+one)*(Lp+Lp+one)*(lk+lk+one)*(lkp+lkp+one)    &
                      *(jk+one)*(jkp+one)*(J+one)*(Jp+one))               &
               * Clebsch_Gordan(lk+lk,0,lkp+lkp,0,nu+nu,0)                &
               * Clebsch_Gordan(L+L,2,Lp+Lp,-2,nu+nu,0)                   &
               * (one + (-1)**(L+p+Lp+pp-nu+16))                          &
               * wigner_6j_symbol(J,Jp,nu+nu,Lp+Lp,L+L,J_f)               &
               * wigner_6j_symbol(J,Jp,nu+nu,jkp,jk,J_i)                  &
               * wigner_6j_symbol(jk,jkp,nu+nu,lkp+lkp,lk+lk,1)           &
               * me * conjg(mep)
         end do
      end do
      !
      beta = -half * N / D
      !
   end function rec_angular_parameter
   !
   !
   subroutine rec_calculate_amplitudes()
   !--------------------------------------------------------------------
   ! Calculates in turn the REC capture amplitudes for all channels and 
   ! REC lines. 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer       :: i, j, n, nw, nu, nocsf
      real(kind=dp) :: energy
      type(nkappa)  :: subshell
      integer, dimension(:), allocatable :: ndx
      !
      n = asf_final%csf_set%nocsf + asf_initial%csf_set%nocsf 
      !
      ! Allocate for a "first time"; it is first dellocated before any usage
      allocate( rec_csp%P(1:n_grasp2k), rec_csp%Q(1:n_grasp2k) )
      allocate( cowf_csp%P(1:10), cowf_csp%Q(1:10) )
      allocate( ndx(1:n) )
      !
      do  i = 1,number_of_lines
         !
         do  j = 1,line(i)%No_channels
            !
            !
            energy = line(i)%e_energy
            call set_configuration_scheme(asf_initial%csf_set,asf_cont%csf_set,&
                     -1,line(i)%channel(j)%kappa,                              &
                     line(i)%totalJ_i,line(i)%parity_i,                        &
                     line(i)%channel(j)%totalJ,line(i)%channel(j)%parity,      &
                     append=.false.,index=ndx)
            !
            rec%No_i = asf_cont%csf_set%nocsf
            allocate( rec%ndx_i(rec%No_i) )
            rec%ndx_i(1:rec%No_i) = ndx(1:rec%No_i)
            !
            nw = asf_cont%csf_set%nwshells
            if (rabs_use_stop  .and. nw /= asf_initial%csf_set%nwshells +1) then
               stop "rec_calculate_amplitudes(): program stop A."
            end if
            !
            ! Calculate the MCP coefficients for the current coupling scheme
            ! as well as the d_rs,  y_k(ab), and x_k(abcd) coefficients
            nocsf = asf_cont%csf_set%nocsf
            call anco_calculate_csf_matrix(asf_cont%csf_set,1,nocsf,1,nocsf)
            call cowf_set_drs_coefficients(line(i)%asfi,asf_cont%csf_set,ndx)
            subshell = nkappa(-1,line(i)%channel(j)%kappa)
            call cowf_set_yk_coefficients(subshell,asf_cont%csf_set)
            call cowf_set_xk_coefficients(subshell,asf_cont%csf_set)
            !
            ! Now iterate the continuum spinors for this channel
            !! cowf_solve_homogeneous_eqn     = .true.
            cowf_start_homogeneous         = .true.
            cowf_phaseshift_wkb            = .true.
            cowf_phaseshift_zero_potential = .false.
            cowf_phaseshift_coulomb        = .false.
	    !
	    cowf_norm_nonrel               = .false.
            cowf_norm_wkb                  = .false.
            cowf_norm_wkb_old              = .true.
	    cowf_average_normalization     = .false.
            call cowf_iterate_csp(energy,subshell)
            !
            rec_csp = cowf_csp
            line(i)%channel(j)%phase = rec_csp%phase
            !
            ! Define the 'extended' configuration scheme for calculating
            ! the REC matrix and allocate memory
            call add_csf_to_basis(asf_final%csf_set,asf_cont%csf_set,      &
                     line(i)%totalJ_f,line(i)%parity_f,index=ndx)
            if (rec_print_csf_scheme) then
               call print_configuration_scheme(6,asf_cont%csf_set)
            end if
            !
            rec%No_f = asf_cont%csf_set%nocsf - rec%No_i
            allocate( rec%ndx_f(rec%No_f) )
            rec%ndx_f(1:rec%No_f) = ndx(1+rec%No_i:asf_cont%csf_set%nocsf)
            allocate( rec%matrix(1:rec%No_i,1:rec%No_f) )
            !
            ! Calculate the 'pure' REC matrix in the given CSF scheme
            ! (i.e. without including mixing coefficients)
            call rec_pure_matrix(i,j,line(i)%channel(j)%multipole_L, &
                                 asf_cont%csf_set)
            !
            call rec_channel_amplitude(i,j)
            !
            deallocate( rec%ndx_f, rec%ndx_i, rec%matrix  )
            call deallocate_csf_basis( asf_cont%csf_set )          
         end do
         !
         ! Calculates all selected properties for the selected transition
         call rec_line_properties(line(i))
      end do
      deallocate( ndx, rec_csp%P, rec_csp%Q)
      !
   end subroutine rec_calculate_amplitudes
   !
   !
   subroutine rec_calculate_Bessel(L,omega_over_c,bessel)
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
   end subroutine rec_calculate_bessel
   !
   !
   subroutine rec_channel_amplitude(i,j)
   !--------------------------------------------------------------------
   ! Calculates the REC amplitude of channel j of line i by summing over 
   ! the 'pure' REC matrix using the proper weights of line i.
   ! 
   ! Calls: angular_momentum_l().
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: i, j
      !
      integer       :: asfi, asff, l, r, rr, s, ss
      real(kind=dp) :: phase, ph, value
      !
      asfi  = line(i)%asfi;  asff = line(i)%asff
      value = zero
      do  r = 1,rec%no_i
         rr = rec%ndx_i(r)
         do  s = 1,rec%no_f
            ss = rec%ndx_f(s)
	    value = value + asf_initial%asf(asfi)%eigenvector(rr) * &
	            rec%matrix(r,s) * asf_final%asf(asff)%eigenvector(ss)
            !
	 end do
      end do
      !
      l     = angular_momentum_l(line(i)%channel(j)%kappa)
      phase = line(i)%channel(j)%phase
      !
      line(i)%channel(j)%amplitude = cmplx(zero,value) 
      line(i)%channel(j)%amplitude = line(i)%channel(j)%amplitude *           &
                                     (cmplx(zero,one)**(-l)) *                &
                                     exp(-cmplx(zero,one)*phase) 
      !
      print *, "***** i,j,gauge, line(i)%channel(j)%amplitude = ",  &
            i,j,line(i)%channel(j)%gauge,line(i)%channel(j)%amplitude
      !
   end subroutine rec_channel_amplitude
   !
   !
   subroutine rec_collect_input()
   !--------------------------------------------------------------------
   ! Collect and proceeds all input for the REC program.
   !
   ! Calls: get_yes_stream().
   !--------------------------------------------------------------------
      !
      integer            :: i, level_i, level_f, score_position, kappa, ierr
      logical            :: yes
      real(kind=dp)      :: en_lower, en_upper, delta_en, &
                            energy, maximal_energy, wavelength
      character(len=20 ) :: string
      character(len=256) :: record, rec_chn_file
      !
      ! Set the maximal number of multipoles
      number_of_multipoles = 30
      !
      ! Open, check, load data from, and close the  .iso  file
      call file_get_isodat_grasp2k()
      !
      ! Determine the units for the printout and further optional input
      ! Determine the units for the printout and further optional input
      call input_energy_unit()
      !
   12 print *, "Enter the maximal energy of the free electrons"//&
               " (in " // trim(energy_unit)                 //&
               ") to built-up the radial grid;"
      read (*, *, err=12) rec_maximal_energy
      !
      if (energy_inverse) then
         rec_maximal_energy = energy_factor * rec_maximal_energy
      else
         rec_maximal_energy = rec_maximal_energy / energy_factor
      end if
      !
      if (rec_maximal_energy < one) goto 12
      !
      ! Determine grid parameters
      call input_grid_parameters("standard")
      !
      wavelength = sqrt( two * pi * pi / rec_maximal_energy )
      !! hp_grasp2k = wavelength / 2323.0_dp
      hp_grasp2k = wavelength / 1247.0_dp
      n_grasp2k  = 30.0_dp / hp_grasp2k  !!  + 100000
      cowf_extent_mtp = 140	 !! Falk
      cowf_extent_mtp = 140000   !! Andrey
      !! cowf_extent_mtp = 9000     !! Andrey 2008
      !
      ! Default speed of light
      c = c_vacuum
      !
      ! Select the electron energies
      maximal_energy = zero
      print *, "The electron energies can be entered either individually or"//&
	       " by an appropriate interval and step size."
      print *, "At storage rings, the electron energies are obtained from"  //&
	       " the projectile energies by: "
      print *, "       E^(electron) [keV] = E^(projectile) [MeV/u]"	    //&
	       " / 1.8228885."   
      print *, " "     
      print *, "The default is a list of individual electron energies;"     //& 
	       " revise this ? "
      yes = get_yes_stream()
      if (.not.yes) then
       2 print *, "Enter another (positive) electron energy in ",      &
		  trim(energy_unit)," enter a negative number if done."
	 read  *, energy
	 if (energy > zero) then
	    if (energy_inverse) then
	       energy = energy_factor * energy
	    else
               energy = energy /energy_factor
            end if
            number_of_electron_energies = number_of_electron_energies + 1
            if (number_of_electron_energies > 100) then
               print *, "The maximum number of electron energies"// &
                        " is 100; the program continues ..."
               goto 4
            end if
            rec_energy_selection(number_of_electron_energies) = energy
            if (energy > maximal_energy)   maximal_energy = energy
            goto 2
         end if
      else
       3 print *, "Enter an interval of electron energies and a corresponding"// &
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
            energy = en_lower;   rec_energy_selection(1) = energy
            number_of_electron_energies = 1
            do
               energy = energy + delta_en
               if (energy > en_upper) exit
               number_of_electron_energies = number_of_electron_energies + 1
               if (number_of_electron_energies > 100) then
                  print *, "The maximum number of electron energies"//&
                           " is 100; choose another interval or stepsize."
                  goto 3
               end if
               rec_energy_selection(number_of_electron_energies) = energy
            end do
         endif
      endif
      !
    4 continue
      !
      ! Enter the maximal number of partial waves
    5 print *, "Enter the maximal (positive) kappa value to determine the "// &
               " partial waves"
      print *, "in the expansion of the electron wave, "// & 
               "i.e. kappa = -1, 1, -2, ..., |max_kappa|:"
      read  *, rec_maximal_kappa
      rec_maximal_kappa = abs(rec_maximal_kappa)
      if (rec_maximal_kappa < 1  .or.  rec_maximal_kappa > 30) then
         print *, "The maximum allowed kappa value must be in the range "// & 
                  "1 <= kappa_max <= 30; reenter ..."
         goto 5
      end if
      !
      ! Now 'overwrite' defaults only if required
      print *, "Modify default set-up and printout of the program ?"
      yes = get_yes_stream()
      if (.not.yes) goto 11
      !
      ! Select individual pairs of transitions
      call input_transition_pairs(number_of_ctransitions)
      !
      print *, "Include exchange interactions into the generation of the"//&
               " continuum waves ?"
      yes = get_yes_stream()
      cowf_solve_homogeneous_eqn = .not.yes
      !
      rec_nonorthogonal_eval = .false.
      !
      print *, "Calculate angular distribution coefficients ?"
      rec_calc_angular_parameter = get_yes_stream()
      !
      print *, "Calculate statistical tensors for the magnetic population"//&
               " of the final states ?"
      rec_calc_statistical_tensor = get_yes_stream()
      !
      print *, "Calculate the anisotropy parameter for a subsequent decay"//&
               " of the final states ?"
      rec_calc_anisotropy_parameter = get_yes_stream()
      !
      if (rec_calc_anisotropy_parameter) then
         print *, "Enter the total angular momentum J_0 for the lower level"//&
                  " of a subsequent decay:"
         read (*, "(a)") record
         record  = adjustl(record)
         rec_anisotropy_J0 = get_dinteger_from_string(record(1:4))
      end if
      !
      print *, "Sort all REC lines in ascending order of energy ?"
      rec_sort_line_energy = get_yes_stream()
      !
      print *, "REC capture cross sections are printed in SI units;"// &
               " use Hartree atomic units instead ?"
      rec_print_cs_in_hartree = get_yes_stream()
      !
      print *, "Print all selected REC lines"// &
               " before the computation starts (this is the default) ?"
      rec_print_selected_lines = get_yes_stream()
      !
      print *, "Print the CSF scheme each time a new one has been built ?"
      rec_print_csf_scheme =  get_yes_stream()
      !
      print *, "Print the results for each individual line immediatly"//&
               " after its computation ?"
      rec_print_each_line = get_yes_stream()
      !
      print *, "Print the final results to a (.chn) channel amplitude file ?"
      rec_print_chn_file = get_yes_stream()
      if (rec_print_chn_file) then
       8 print *, "Enter a file name for the  rec.chn  file:"
         read *,  rec_chn_file
         call file_open(27,rec_chn_file,"formatted  ","new",ierr)
         !
         if (ierr /= 0) goto 8
      end if
      !
      print *, "Enter an (overall) shift for the atomic transition energies"//&
               " which applies to all transitions (in"//                      &
               trim(energy_unit)//"):"
      print *, " Use  0.  or   <cr>  if no shift need to be applied."
      read (*, *, err=10) rec_energy_shift
      !
      if (energy_inverse) then
         rec_energy_shift = energy_factor * rec_energy_shift
      else
         rec_energy_shift = rec_energy_shift / energy_factor
      end if
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
      if (rec_print_cs_in_hartree) then
	 rec_cs_factor = one
	 rec_cs_unit   = "a.u.   "
      else
	 rec_cs_factor = bohr_radius_in_cm*bohr_radius_in_cm / 1.0e-24_dp
	 rec_cs_unit   = "b	 "
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
   end subroutine rec_collect_input
   !
   !
   subroutine rec_initialize()
   !--------------------------------------------------------------------
   ! Initializes the computation of REC capture cross sections, angular
   ! distributions and others. 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer :: i, mtp, noasf, nocsf, number_of_rwf
      !
      print *, "Initialize the set-up of the overlap integrals and "// &
               "REC lines ..."
      call rec_set_lines()
      call rec_set_overlaps(6)
      !
      print *, "   ... initialization complete."
      print *, " "
      !
      ! Print the selected REC lines before the computation starts
      if (rec_print_selected_lines) then
         call rec_print_lines(6)
         call rec_print_lines(24)
      end if
      !
      ! Make a 'copy' of asf_initial into asf_bound and wave_initial into 
      ! wave_bound (kept in rabs_cowf) to use the COWF component
      noasf = asf_initial%noasf;   nocsf = asf_initial%csf_set%nocsf
      allocate( asf_bound%asf(1:noasf) )
      do  i = 1,noasf
         allocate( asf_bound%asf(i)%eigenvector(1:nocsf) )
      end do
      asf_bound%noasf          = asf_initial%noasf
      asf_bound%average_energy = asf_initial%average_energy
      asf_bound%asf            = asf_initial%asf
      !
      number_of_rwf            = wave_initial%number_of_rwf
      wave_bound%number_of_rwf = wave_initial%number_of_rwf
      allocate( wave_bound%rwf(1:number_of_rwf) )
      do  i = 1,number_of_rwf
         mtp = wave_initial%rwf(i)%mtp
         allocate( wave_bound%rwf(i)%P(1:mtp), wave_bound%rwf(i)%Q(1:mtp) )
         wave_bound%rwf(i) = wave_initial%rwf(i)
      end do
      !
   end subroutine rec_initialize
   !
   !
   subroutine rec_line_properties(tline)
   !--------------------------------------------------------------------
   ! Calculates all selected REC properties of line i from the amplitudes 
   ! of the individual REC channels. This concerns the total cross sections 
   ! as well as angular distribution coefficients.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(rec_line), intent(inout) :: tline
      !
      integer       :: j, nu, stream
      real(kind=dp) :: amp2, energy, p_energy, e_energy, cs, csc, csb, sumb, &
                       sumc, lorentz_beta, lorentz_gamma, wa, weight
      !
      character(len=1), dimension(0:1) :: mult
      !
      !
      character(len=256) :: test_file
      real(kind=dp)      :: wc, wd
      integer            :: ierr
      !
      mult(0) = "M";   mult(1) = "E"
      !
      lorentz_gamma   = one + tline%e_energy / (c*c)
      lorentz_beta    = sqrt( one - (one/(lorentz_gamma*lorentz_gamma)) )
      !
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
      !
      ! Add proper weight coefficients
      weight = two*two*two*pi*pi*pi * tline%p_energy / (137.036_dp**3)                    &
               / ((lorentz_gamma**2)*(lorentz_beta**2)*(tline%totalJ_i + one)) 
      !
      tline%cs_babushkin = weight * sumb    
      tline%cs_coulomb   = weight * sumc  
      !
      ! Calculate the angular distribution coefficients
      if (rec_calc_angular_parameter) then
         do  nu = 1,rec_maximal_nu
            tline%beta_b(nu) = rec_angular_parameter("B",nu,tline)
            tline%beta_c(nu) = rec_angular_parameter("C",nu,tline)
         end do
      end if
      !
      ! Calculate the statistical tensors for the final sublevel population
      if (rec_calc_statistical_tensor) then
         do  nu = 0,rec_maximal_tensor,2
            tline%rho_b(nu) = rec_statistical_tensor("B",nu,tline)
            tline%rho_c(nu) = rec_statistical_tensor("C",nu,tline)
         end do
      end if
      !
      ! Print a short summary of the line if required
      !
      if (rec_print_each_line) then
         stream = 6
       3 write(stream,*) " "
         write(stream,14) " Results for REC line from the transition ",  &
                         tline%level_i," - ",tline%level_f,": " //       &
                         trim(angular_momentum_string(tline%totalJ_i)) //&
                         " " // tline%parity_i // "   ----> " //         &
                         trim(angular_momentum_string(tline%totalJ_f)) //&
                         " " // tline%parity_f               
         write(stream,*) " --------------------------------------------"//  &
                         "----------------------------------"
         write(stream,*) " "
         if (energy_inverse) then
            p_energy = energy_factor / tline%p_energy
            e_energy = energy_factor / tline%e_energy
         else
            p_energy = energy_factor * tline%p_energy
            e_energy = energy_factor * tline%e_energy
         end if
         csb = tline%cs_babushkin * rec_cs_factor
         csc = tline%cs_coulomb   * rec_cs_factor
         write(stream,*) "   Photon energy       = ",  &
                         p_energy,trim(energy_unit)
         write(stream,*) "   Electron energy     = ",  &
                         e_energy,trim(energy_unit)
         write(stream,*) "   Total cross section = ",csb,trim(rec_cs_unit),&
                         "  (Babushkin gauge) "
         write(stream,*) "                       = ",csc,trim(rec_cs_unit),&
                         "  (Coulomb gauge) "
         write(stream,*) " "
         write(stream,*) "   Lorentz_beta        = ",lorentz_beta
         write(stream,*) "   Lorentz_gamma       = ",lorentz_gamma
         !
         if (rec_calc_angular_parameter) then
            write(stream,*) " "
            do  nu = 1,min(9,rec_maximal_nu)
               if (abs(tline%beta_b(nu)) > 1.0e-15  .or.  &
                   abs(tline%beta_c(nu)) > 1.0e-15) then
                  write(stream,4) nu,tline%beta_b(nu)
                  write(stream,5) nu,tline%beta_c(nu)
               end if
            end do
	    !
	    if (rec_maximal_nu > 9) then
               do  nu = 10,rec_maximal_nu
                  if (abs(tline%beta_b(nu)) > 1.0e-15  .or.  &
                      abs(tline%beta_c(nu)) > 1.0e-15) then
                     write(stream,6) nu,tline%beta_b(nu)
                     write(stream,7) nu,tline%beta_c(nu)
                  end if
               end do
	    end if
         end if
         !
         if (rec_calc_statistical_tensor) then
            write(stream,*) " "
            do  nu = 0,rec_maximal_tensor,2
                  write(stream,12) nu,tline%rho_b(nu)
                  write(stream,13) nu,tline%rho_c(nu)
            end do
         end if
         !
         if (rec_calc_anisotropy_parameter) then
            wa = ((-1)**((2+tline%totalJ_f+rec_anisotropy_J0)/2)) *          &
                 sqrt( (three*tline%totalJ_f + three)/two ) *                &
                 wigner_6j_symbol(2,2,4,tline%totalJ_f,tline%totalJ_f,       &
                                  rec_anisotropy_J0)
            !
            write(stream,*) " "
            write(stream,*) "   Anisotropy parameter for the decay to a"//   &
                            " subsequent level with J0 = ",                  &
                            trim(angular_momentum_string(rec_anisotropy_J0,4))
            write(stream,*) " "
            write(stream,*) "   b^an_2 (Babushkin)  = ",real(tline%rho_b(2)/ &
                                                             tline%rho_b(0))*wa
            write(stream,*) "   b^an_2 (Coulomb)    = ",real(tline%rho_c(2)/ &
                                                             tline%rho_c(0))*wa
         end if
         !
         write(stream,*) " "
         write(stream,1) trim(rec_cs_unit)
         write(stream,*) &
              " -----------------------------------------------------------"//&
              "---------------------------------------------"
         do  j = 1,tline%No_channels
            amp2 = tline%channel(j)%amplitude *                              &
                      conjg(tline%channel(j)%amplitude) 
            cs   = tline%channel(j)%amplitude *                              &
                    conjg(tline%channel(j)%amplitude) * weight * rec_cs_factor        
            write(stream,2) orbital_symmetry(tline%channel(j)%kappa),        &
                    trim(angular_momentum_string(tline%channel(j)%totalJ,4)),&
                            tline%channel(j)%parity,                         & 
                            mult(tline%channel(j)%multipole_pi),             &
                            tline%channel(j)%multipole_L,                    &
                            tline%channel(j)%gauge(1:9),                     &
                            tline%channel(j)%amplitude,amp2,cs,              &
                            mod(tline%channel(j)%phase,two*pi)
         end do
         write(stream,*) " "
       1 format("  Kappa   Total J^P  Mp     Gauge             Amplitude ",  &
                "         Amplitude^2       CS (",a,")        Phase")
       2 format(4x,a2,5x,a4,a1,5x,a1,i2,3x,a9,3x,1pe10.3,2x,1pe10.3,3x,      &
                1pe12.5,3x,1pe12.5,3x,1pe10.3)
       4 format(4x,"beta_",i1,"  (Babushkin) = ",1pe15.8,2x,1pe15.8) 
       5 format(4x,"beta_",i1,"  (Coulomb)   = ",1pe15.8,2x,1pe15.8) 
       6 format(4x,"beta_",i2, " (Babushkin) = ",1pe15.8,2x,1pe15.8) 
       7 format(4x,"beta_",i2, " (Coulomb)   = ",1pe15.8,2x,1pe15.8) 
      12 format(4x,"rho_",i1,"   (Babushkin) = ",1pe15.8,2x,1pe15.8) 
      13 format(4x,"rho_",i1,"   (Coulomb)   = ",1pe15.8,2x,1pe15.8) 
      14 format(a,i4,a,i4, a)
         !
         ! Re-cycle once on stream 24
         if (stream == 6) then
            stream = 24;   goto 3
         end if
      end if
      !
   end subroutine rec_line_properties
   !  
   !
   subroutine rec_print_amplitudes(stream)
   !--------------------------------------------------------------------
   ! Prints the information about all (selected) REC lines and amplitudes 
   ! to a (.chn) channel amplitude file on stream.
   !
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      !
      integer             ::  i, j, channels
      !
      write(stream,*) &
         "Transition data and amplitudes for REC capture lines are printed "
      write(stream,*) &
         "in the format (one line per REC channel): "
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
            write(stream,1) line(i)%asfi,line(i)%asff,               &
                            line(i)%level_i,line(i)%level_f,         &
                            line(i)%totalJ_i,line(i)%totalJ_f,       &
                            line(i)%parity_i,line(i)%parity_f,       &
                            line(i)%p_energy,line(i)%e_energy,       &
                            line(i)%channel(j)%kappa,                &
                            line(i)%channel(j)%totalJ,               &
                            line(i)%channel(j)%parity,               &
                            line(i)%channel(j)%multipole_L,          &
                            line(i)%channel(j)%multipole_pi,         &
                            line(i)%channel(j)%gauge,                &
                            line(i)%channel(j)%phase,                &
                            line(i)%channel(j)%amplitude
         end do
      end do
      !
      1 format(2i5,3x,2i5,3x,2i3,2x,2a3,1x,e14.7,1x,e14.7,i4,i5,a3,2x, &
               i2,"^(pi=",i1,")",2x,a9,2x,e14.7,2x,e14.7,1x,e14.7)
      !
   end subroutine rec_print_amplitudes
   !
   !
   subroutine rec_print_lines(stream)
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
      logical                :: test_gt9
      real(kind=dp)          :: e_energy, p_energy
      integer, dimension(20) :: kappa_channel
      !
      character(len=1), dimension(0:1) :: mult
      !
      mult(0) = "M";   mult(1) = "E"
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
         !
         test_gt9 = .false.
         do  j = 1,line(i)%No_channels
            if (line(i)%channel(j)%multipole_L > 9)  test_gt9 = .true.
         end do
         !
         if (.not.test_gt9) then
         write(stream,2) line(i)%level_i,line(i)%level_f,            &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)), &
                  line(i)%parity_i,                                  &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)), &
                  line(i)%parity_f,p_energy,e_energy,                &
                  (orbital_symmetry(line(i)%channel(j)%kappa),       &
                   trim(mult(line(i)%channel(j)%multipole_pi)),      &
                   line(i)%channel(j)%multipole_L,                   &
                   line(i)%channel(j)%gauge(1:1),                    &
                   trim(angular_momentum_string(line(i)%channel(j)%totalJ,4)),&
                   line(i)%channel(j)%parity,j=1,line(i)%No_channels)
         else
         write(stream,3) line(i)%level_i,line(i)%level_f,            &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)), &
                  line(i)%parity_i,                                  &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)), &
                  line(i)%parity_f,p_energy,e_energy,                &
                  (orbital_symmetry(line(i)%channel(j)%kappa),       &
                   trim(mult(line(i)%channel(j)%multipole_pi)),      &
                   line(i)%channel(j)%multipole_L,                   &
                   line(i)%channel(j)%gauge(1:1),                    &
                   trim(angular_momentum_string(line(i)%channel(j)%totalJ,4)),&
                   line(i)%channel(j)%parity,j=1,line(i)%No_channels)
         end if
      end do
      write(stream,4) 
    1 format(/," The following ",i5," lines are selected:",                &
        //,"     I-level-F     I--J^P--F       Photon Energy    Electron", &
           " energy        REC capture channels  ",                        &
         /,"                                     (in ",a4,")  ",           &  
         /,4x,133("-") )
    2 format(4x,i4," -",i4,3x,a4,a1,3x,a4,a1,4x,1pe14.7,4x,1pe14.7,3x,     &
                   4(a2,"(",a1,i1,";",a1,") J=",a,a1,2x),  &
         200(/,69x,4(a2,"(",a1,i1,";",a1,") J=",a,a1,2x)), &
             /,69x,4(a2,"(",a1,i1,";",a1,") J=",a,a1,2x))
    3 format(4x,i4," -",i4,3x,a4,a1,3x,a4,a1,4x,1pe14.7,4x,1pe14.7,3x,     &
                   4(a2,"(",a1,i2,";",a1,") J=",a,a1,2x),  &
         200(/,69x,4(a2,"(",a1,i2,";",a1,") J=",a,a1,2x)), &
             /,69x,4(a2,"(",a1,i2,";",a1,") J=",a,a1,2x))
    4 format(4x,133("-") )
      !
   end subroutine rec_print_lines
   !
   !
   subroutine rec_print_results(stream)
   !--------------------------------------------------------------------
   ! Writes the REC capture cross sections and angular distribution 
   ! parameters as well as further properties in a neat summary to the 
   ! .sum file.
   !
   ! Calls:  angular_momentum_string().
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      !
      integer       :: i
      real(kind=dp) :: p_energy, e_energy, csc, csb, csc_tot, csb_tot, wa, wb
      !
      write(stream,1)
    1 format(/ &
         /18x,"==============================================================", &
         /18x,"|  Summary of all REC Cross sections and Angular Parameters  |", &
         /18x,"==============================================================", &
             // )
    2 format( / 1x,95("-"),                                                     &
              / 2x,"LevI-LevF   I- J / Parity -F     ",                         &
                   "   e-Energy          omega       CS-Babushkin   ",          &
                   "  CS-Coulomb   ",                                           &
              / 2x,"                                  ",                        &
                   "  (",a,")             (",a,")            (a.u.)    ",       &
                   "       (a.u.)         "   ,                                 &
              / 1x,95("-") )
    3 format( / 1x,95("-"),                                                     &
              / 2x,"LevI-LevF   I- J / Parity -F    ",                          &
                   "   e-Energy          omega        CS-Babushkin   ",         &
                   "  CS-Coulomb    ",                                          &
              / 2x,"                                  ",                        &
                   "  (",a,")             (",a,")             (",a,")     ",    &
                   "        (",a,")           "   ,                             &
              / 1x,95("-") )
    4 format(   2x,i3,' -',i3,3x,a4,a2,4x,a4,a2,6x,2(1pe12.5,5x),               &
                   2(1pe10.3,6x) )
    7 format( / 1x,127("-"),                                                    &
              / 2x,"LevI-LevF   I- J / Parity -F    ",                          &
                   "   e-Energy          omega       Gauge  ",                  &
                   "      Beta_1       Beta_2       Beta_3       Beta_4   ",    &
              / 37x,"  (",a,")            (",a,")    ",                         &
              / 1x,127("-") )
    8 format(   2x,i3,' -',i3,3x,a4,a2,4x,a4,a2,4x,2(2x,1pe12.5,3x),a9,2x,      &
                           4(1pe10.3,3x) )
    9 format(   67x,a9,2x, 4(1pe10.3,3x) )
   10 format(   1x,95("-") )
   11 format(   1x,127("-") )
   12 format( / 1x,127("-"),                                                    &
              / 2x,"LevI-LevF   I- J / Parity -F    ",                          &
                   "   e-Energy          omega       Gauge  ",                  &
                   "       rho_0        rho_2        rho_4        rho_6     ",  &
              / 37x,"  (",a,")            (",a,")                           ",  &
              / 1x,127("-") )
   13 format(   2x,i3,' -',i3,3x,a4,a2,4x,a4,a2,4x,2(2x,1pe12.5,3x),a9,2x,      &
                           4(1pe10.3,3x) )
   14 format( / 1x,88("-"),                                                     &
              / 2x,"LevI-LevF   I- J / Parity -F    ",                          &
                   "   e-Energy          omega       Gauge         beta_2",     &
                   "     beta_4",   &
              / 37x,"  (",a,")            (",a,")                           ",  &
              / 1x,88("-") )
   15 format(   2x,i3,' -',i3,3x,a4,a2,4x,a4,a2,4x,2(2x,1pe12.5,3x),a9,2x,      &
                           4(1pe10.3,3x) )
   16 format(   1x,88("-") )
      !
      ! Print the REC cross sections of this calculation
      !
      write(stream,*) "Individual and total REC cross sections :"
      write(stream,*) "-----------------------------------------"
      if (rec_print_cs_in_hartree) then
         write(stream,2) trim(energy_unit), trim(energy_unit)
      else
         write(stream,3) trim(energy_unit), trim(energy_unit),      &
                         trim(rec_cs_unit),trim(rec_cs_unit)
      end if
      !
      csc = zero;   csb = zero
      do  i = 1,number_of_lines
         csb_tot = csb_tot + line(i)%cs_babushkin
         csc_tot = csc_tot + line(i)%cs_coulomb
      end do
      !
      do  i = 1,number_of_lines
         if (energy_inverse) then
            p_energy = energy_factor / line(i)%p_energy
            e_energy = energy_factor / line(i)%e_energy
         else
            p_energy = energy_factor * line(i)%p_energy
            e_energy = energy_factor * line(i)%e_energy
         end if
         csb = line(i)%cs_babushkin * rec_cs_factor
         csc = line(i)%cs_coulomb   * rec_cs_factor
	 !
	 csb = csb / (line(i)%totalJ_f+one)
	 csc = csc / (line(i)%totalJ_f+one)
         !
         if (.not.rec_print_only_gt_1percent  .or.                        &
            (rec_print_only_gt_1percent .and.                             &
               abs(csb/csb_tot) > 0.01_dp)) then
            write(stream,4) line(i)%level_i, line(i)%level_f,   &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)),  &
                            line(i)%parity_i,                         &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)),  &
                            line(i)%parity_f,                         &
                            e_energy,p_energy,csb,csc
         end if
      end do
      write(stream,10)
      write(stream,*) " "
      !
      !
      ! Print REC angular parameters
      !
      if (rec_calc_angular_parameter) then
         write(stream,*) " "
         write(stream,*) "REC angular parameters:"
         write(stream,*) "-----------------------"
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
                            line(i)%parity_f,                               &
                            e_energy,p_energy,                              &
                            "Babushkin",                                    &             
                            real(line(i)%beta_b(1)),real(line(i)%beta_b(2)),&
                            real(line(i)%beta_b(3)),real(line(i)%beta_b(4))
            write(stream,9) "Coulomb  ",                                    &             
                            real(line(i)%beta_c(1)),real(line(i)%beta_c(2)),&
                            real(line(i)%beta_c(3)),real(line(i)%beta_c(4))
                            
         end do
         write(stream,11)
         write(stream,*) " "
      end if
      !
      !
      ! Print REC statistical tensors
      !
      if (rec_calc_statistical_tensor) then
         write(stream,*) " "
         write(stream,*) "REC statistical tensors:"
         write(stream,*) "------------------------"
         write(stream,12) trim(energy_unit), trim(energy_unit)
         do  i = 1,number_of_lines
            if (energy_inverse) then
               p_energy = energy_factor / line(i)%p_energy
               e_energy = energy_factor / line(i)%e_energy
            else
               p_energy = energy_factor * line(i)%p_energy
               e_energy = energy_factor * line(i)%e_energy
            end if
            !
            write(stream,13) line(i)%level_i, line(i)%level_f,              &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)),        &
                            line(i)%parity_i,                               &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)),        &
                            line(i)%parity_f,                               &
                            e_energy,p_energy,                              &
                            "Babushkin",                                    &             
                            real(line(i)%rho_b(0)),real(line(i)%rho_b(2)),  &
                            real(line(i)%rho_b(4)),real(line(i)%rho_b(6))
            write(stream,9) "Coulomb  ",                                    &             
                            real(line(i)%rho_c(0)),real(line(i)%rho_c(2)),  &
                            real(line(i)%rho_c(4)),real(line(i)%rho_c(6))
         end do
         write(stream,11)
         write(stream,*) " "
      end if
      !
      !
      ! Print REC anisotropy parameters
      !
      if (rec_calc_anisotropy_parameter) then
         write(stream,*) " "
         write(stream,*) "REC anisotropy parameters for the subsequent"    //&
                         " decay to a lower level with J0 = "              //&
                         trim(angular_momentum_string(rec_anisotropy_J0,4))//&
                         ":"
         write(stream,*) "--------------------------------------------"//&
                         "---------------------------------------"
         write(stream,14) trim(energy_unit), trim(energy_unit)
         do  i = 1,number_of_lines
            if (energy_inverse) then
               p_energy = energy_factor / line(i)%p_energy
               e_energy = energy_factor / line(i)%e_energy
            else
               p_energy = energy_factor * line(i)%p_energy
               e_energy = energy_factor * line(i)%e_energy
            end if
            !
            wa = ((-1)**((2+line(i)%totalJ_f+rec_anisotropy_J0)/2)) *       &
                 sqrt( (three*line(i)%totalJ_f + three)/two ) *             &
                 wigner_6j_symbol(2,2,4,line(i)%totalJ_f,line(i)%totalJ_f,  &
                 rec_anisotropy_J0)
            !
            wb = ((-1)**((2+line(i)%totalJ_f+rec_anisotropy_J0)/2)) *       &
                 sqrt( (three*line(i)%totalJ_f + three)/two ) *             &
                 wigner_6j_symbol(4,4,8,line(i)%totalJ_f,line(i)%totalJ_f,  &
                 rec_anisotropy_J0)
            !
            write(stream,13) line(i)%level_i, line(i)%level_f,              &
                  trim(angular_momentum_string(line(i)%totalJ_i,4)),        &
                            line(i)%parity_i,                               &
                  trim(angular_momentum_string(line(i)%totalJ_f,4)),        &
                            line(i)%parity_f,                               &
                            e_energy,p_energy,                              &
                            "Babushkin",                                    &
                            wa*real(line(i)%rho_b(2)/line(i)%rho_b(0)),     &
                            wb*real(line(i)%rho_b(4)/line(i)%rho_b(0))
            write(stream,9) "Coulomb  ",                                    &
                            wa*real(line(i)%rho_c(2)/line(i)%rho_c(0)),     &
                            wb*real(line(i)%rho_c(4)/line(i)%rho_c(0))
         end do
         write(stream,16)
         write(stream,*) " "
      end if
      !
      !
   end subroutine rec_print_results
   !
   !
   subroutine rec_pure_matrix(tr,ch,nu,csf_set)
   !--------------------------------------------------------------------
   ! Calculates the 'pure' REC matrix for the given configuration scheme
   ! csf_set. The first no_i CSF belong to the initial-state representation
   ! and the following no_i+1,...,no_i+no_f to the final states.
   ! The procedure takes into account ...
   !
   ! Calls:  
   !--------------------------------------------------------------------
      !
      integer, intent(in)	  :: tr, ch, nu
      type(csf_basis), intent(in) :: csf_set
      !
      integer			  :: i, ia, ib, j, jb, ja, no_T_coeff, &
				     r, s, ss, t, parity     
      type(nkappa)		  :: aa, bb, cc, dd
      real(kind=dp)		  :: aweight, wa, wb, ph
      type(nkappa)                :: Ta, Tb
      !
      if (mod(nu,2) == 1) then
         parity = (-1)**line(tr)%channel(ch)%multipole_pi
      else if (mod(nu,2) == 0) then
         parity = (-1)**(line(tr)%channel(ch)%multipole_pi+1)
      else
         stop "rec_pure_matrix(): program stop A."
      end if
      !
      if (rabs_use_stop   .and.   csf_set%nocsf /= rec%no_i+rec%no_f) then
         stop "rec_pure_matrix(): program stop B."
      end if
      !
      rec%matrix = zero
      !
      if (rec_nonorthogonal_eval) then
         call nonorth_initialize(csf_set,wave_final,wave_initial)
      end if
      !
      do  r = 1,rec%No_i
         do  s = rec%No_i+1,rec%No_i+rec%No_f
            ss = s - rec%No_i
	    if (rec_nonorthogonal_eval) then
	       !
	       ! Use an expansion into Slater determinants to evaluate the
	       ! Einstein matrix
	       Aoperator%particle = 1
	       Aoperator%rank     = nu
	       Aoperator%parity   = parity
	       !
               call nonorth_calculate_csf_pair(csf_set,Aoperator,r,s,no_T_coeff)
	    else
	       !
	       ! Use standard Racah algebra to evaluate the many-electron
	       ! Einstein matrix
               !! call anco_calculate_csf_pair_1p(nu,csf_set,r,s,no_T_coeff)
	       !
	       number_of_mct_coefficients = 0
               call mct_generate_coefficients(r,r,s,s,csf_set,parity,nu)
	       no_T_coeff = number_of_mct_coefficients
	    end if
            !
            ! Cycle over all angular coefficients
            do  t = 1,no_T_coeff
	       if (rec_nonorthogonal_eval) then
	          Ta      = nonorth_T_list(t)%a
	          Tb      = nonorth_T_list(t)%b
                  ja      = angular_momentum_j(Ta%kappa)
                  aweight = nonorth_T_list(t)%T * sqrt( ja+one )
		  !
		  ia = 0
                  ib = 0
                  do  i = 1,wave_final%number_of_rwf
                     if (wave_final%rwf(i)%orbital == Tb) then
                        ib = i
                        exit
                     end if
                  end do
	       else
	          ia       = mct_list(t)%a
	          ib       = mct_list(t)%b
		  Ta       = rec_csp%orbital
		  Tb       = wave_final%rwf(ib)%orbital
                  ja       = angular_momentum_j(Ta%kappa)
	          aweight  = mct_list(t)%T / sqrt( ja+one )
		  !!x print *, "Orth: ia, ib, Ta, Tb, aweight = ", &
		  !!x                 ia, ib, Ta, Tb, aweight
	       end if
               !
               if (Ta%n < 0  .and.  ib /= 0) then
                  jb = angular_momentum_j(Tb%kappa)
                  ja = angular_momentum_j(Ta%kappa)
                  wb = multipole_reduced_M_capture(line(tr)%p_energy,        &
                                        line(tr)%channel(ch)%multipole_L,    &
                                        line(tr)%channel(ch)%multipole_pi,   &
                                        line(tr)%channel(ch)%gauge,          & 
                                        rec_csp,wave_final%rwf(ib))
                  !
                  rec%matrix(r,ss) = rec%matrix(r,ss)                        &
                        + aweight * sqrt(line(tr)%channel(ch)%totalJ+one) * wb
		  !!x print *, "*** (na,kappa_a),(nb,kappa_b), me, aweight, sq = ",&
		  !!x          Ta,Tb,wb,aweight,                                   &
		  !!x 	   sqrt(line(tr)%channel(ch)%totalJ+one)
                  !
                  !! if (line(tr)%channel(ch)%gauge /= "Coulomb  ") then
                  write(97,111)  orbital_name(Ta%n,Ta%kappa),                &
                                 orbital_name(Tb%n,Tb%kappa),                &
                                        line(tr)%channel(ch)%multipole_L,    &
                                        line(tr)%channel(ch)%multipole_pi,   &
                                        line(tr)%p_energy,                   &
                                        Ta%n,Ta%kappa,Tb%n,Tb%kappa,         &
                                        line(tr)%channel(ch)%multipole_L,    &
                                        line(tr)%channel(ch)%multipole_pi,wb       
                  !
              111 format("*** Reduced ME orb_a, orb_b, ",                    &
                         "L^pi, (n,kappa)_a, (n,kappa)_b, L, pi, rme = ",    &
                         a4,3x,a4,3x,i2,"^(",i1,")",2x,1e12.5,2x,            &
                         2i3,2x,2i3,2x,2i3,4x,1e19.12) 
                  !! end if      
                  !
               else
                  stop "rec_pure_matrix(): program stop C."
               end if
               !
               !!x print *, "c: r, s, me = ",r,s,rec%matrix(r,ss)
            end do
	    !
         end do
      end do
      !
   end subroutine rec_pure_matrix
   !
   ! 
   function rec_reduced_M_integral_old(tr,L,p,gauge,rwf_f,rwf_i) &
                                                          result(red_me)
   !--------------------------------------------------------------------
   ! Calculates the value of an (one-electron) M integral for the coupling
   ! of electrons with the radiation field for a given multipole and
   ! gauge. It calculates all of the required Bessel functions; 
   ! the radial wave functions are provided by the derived data structures 
   ! rwf_f and rwf_i.
   !
   ! Calls: angular_momentum_j(), orbital_name(), quad_grasp2k(),
   !        wigner_3j_symbol().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: tr, L, p
      character(len=9), intent(in)       :: gauge
      type(orbital_function), intent(in) :: rwf_f, rwf_i
      real(kind=dp)                      :: red_me
      !
      integer       :: ja, jb, kapa, kapb, mtp, phase
      real(kind=dp) :: arg, factor, red_Coulomb, red_Gauge, wa, wb
      real(kind=dp), dimension(:), allocatable :: ta
      real(kind=dp), dimension(1:n_grasp2k)    :: bessel
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
      !
      ! Test parity rules
      if (p == 0) then
         ! Magnetic case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 0)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 2) ) then
         else
            return
         end if
      else
         ! Electric case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 2)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 0) ) then
         else
            return
         end if
      end if
      !
      !
      ! Evaluate factor multiplying Mbar(a,b); use one-particle matrix elements
      ! of Brink-Satchler type
      phase  = (ja + 1)/2
      factor = ((-one)**phase) *sqrt(ja+one)*Clebsch_Gordan(ja,1,L+L,0,jb,1)
      !
      if (abs(factor) < eps10) then
         red_me = zero
         return
      end if
      !
      mtp = min( rwf_f%mtp, rwf_i%mtp)
      allocate( ta(1:mtp+10) )
      !
      if (p == 0) then
         ! Magnetic case
         ta(1) = zero
         ta(2:mtp) = ( rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +                   &
                       rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) ) * rp_grasp2k(2:mtp)
         !
         arg  = line(tr)%p_energy / c         
         call rec_calculate_Bessel(L,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_me = -(L+L+one)*(kapa+kapb)/sqrt(L*(L+one)) * quad_grasp2k(ta,mtp)
         red_me = factor * red_me
      else if (p == 1) then
         ! Electric case
         ta(1) = zero
         !
         ! Tabulate coulomb integrand
         wa    = sqrt(L/(L+one))*(kapa-kapb);  wb = sqrt(L/(L+one))*(L+one)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = line(tr)%p_energy / c         
         call rec_calculate_Bessel(L+1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Coulomb = factor * quad_grasp2k(ta,mtp)
         !!x print *, "a: red_Coulomb = ",red_Coulomb 
         !
         wa    = - sqrt((L+one)/L)*(kapa-kapb);  wb = sqrt((L+one)/L)*L
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = line(tr)%p_energy / c         
         call rec_calculate_Bessel(L-1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Coulomb = red_Coulomb + factor * quad_grasp2k(ta,mtp)
         !
         if (gauge == "Coulomb  ") then
            red_me = red_Coulomb
            if (.false.) &
               write(99,*) "factor, red_Coulomb = ",factor, red_me
            deallocate( ta )
            !!x print *, "c: red_me = ",red_me
            return
         else if (gauge == "Babushkin") then
         else if (rabs_use_stop) then
            stop "rec_reduced_M_integral(): program stop B."
         end if
         !
         wa    = (kapa-kapb);  wb = (L+1)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = line(tr)%p_energy / c         
         call rec_calculate_Bessel(L+1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Gauge = factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         !
         wa    = (kapa-kapb);  wb = -L
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = line(tr)%p_energy / c         
         call rec_calculate_Bessel(L-1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Gauge = red_Gauge + factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         !
         wa    = -(L+L+one)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%P(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%Q(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = line(tr)%p_energy / c         
         call rec_calculate_Bessel(L,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Gauge = red_Gauge + factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         red_me    = red_Coulomb + red_Gauge
         !
      else
         stop "rec_reduced_M_integral(): program stop C."
      end if
      !
      deallocate( ta ) 
      !
   end function rec_reduced_M_integral_old
   !
   !
   subroutine rec_set_lines()
   !--------------------------------------------------------------------
   ! Determines how many and which of the REC lines need to be calculated 
   ! to determine the cross sections and angular distribution coefficients.
   ! It also initializes the array line of type(rec_lines).
   ! The default is that , for number_of_lines == 0 all lines with a positive 
   ! photon energy are calculated; for number_of_lines /= 0,
   ! lines for individual atomic transitions were selected during input time
   ! and will be initialized instead; in this case also negative line energies
   ! are allowed and may be overritten by an appropriate shift of the total
   ! energies.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer          :: asfi, asff, i, imin, j, k, l, lkappa, m, nchannels, & 
                          ka, kas, kappa, nt, nch, p, totalJ
      real(kind=dp)    :: arg, e_energy, p_energy, energy_exp, energy_min
      character(len=1) :: parity
      !
      logical, dimension(:), allocatable    :: ask_for
      integer, dimension(900)               :: wa_kappa, wa_totalJ,          &
                                               wa_multipole_L, wa_multipole_pi 
      character(len=1), dimension(900)      :: wa_parity
      character(len=9), dimension(900)      :: wa_gauge
      real(kind=dp), dimension(1:n_grasp2k) :: bessel
      !
      ! Determine the total number of lines
      nt = 0
      if (number_of_ctransitions == 0) then
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_electron_energies
                  if (asf_initial%asf(i)%energy - asf_final%asf(j)%energy &
                                        + rec_energy_selection(k) > zero) then
                     nt = nt + 1
                  end if
               end do
            end do
         end do
      else
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_ctransitions
                  if ( (select_level_i(k) == asf_initial%asf(i)%level_No .or. &
                        select_level_i(k) == 0)                         .and. &
                       (select_level_f(k) == asf_final%asf(j)%level_No   .or. &
                        select_level_f(k) == 0) ) then
                     do  l = 1,number_of_electron_energies
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
      !!x print *, "A - nt = ",nt
      !
      ! Now initialize all lines
      nt = 0
      if (number_of_ctransitions == 0) then
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_electron_energies
                  if (asf_initial%asf(i)%energy - asf_final%asf(j)%energy &
                                        + rec_energy_selection(k) > zero) then
                     nt = nt + 1
                     line(nt)%asfi     = i;   line(nt)%asff = j
                     line(nt)%level_i  = asf_initial%asf(i)%level_No
                     line(nt)%level_f  = asf_final%asf(j)%level_No
                     line(nt)%e_energy = rec_energy_selection(k)
                     line(nt)%p_energy = rec_energy_selection(k) +      &
                           asf_initial%asf(i)%energy - asf_final%asf(j)%energy
                  end if
               end do
            end do
         end do
      else
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_ctransitions
                  if ( (select_level_i(k) == asf_initial%asf(i)%level_No .or. &
                        select_level_i(k) == 0)                         .and. &
                       (select_level_f(k) == asf_final%asf(j)%level_No   .or. &
                        select_level_f(k) == 0) ) then
                     do  l = 1,number_of_electron_energies
                        nt = nt + 1
                        line(nt)%asfi     = i;   line(nt)%asff = j
                        line(nt)%level_i  = asf_initial%asf(i)%level_No
                        line(nt)%level_f  = asf_final%asf(j)%level_No
                        line(nt)%e_energy = rec_energy_selection(l)
                        line(nt)%p_energy = rec_energy_selection(l) +    & 
                           asf_initial%asf(i)%energy - asf_final%asf(j)%energy
                     end do
                  end if
               end do
            end do
         end do
      end if
      !
      number_of_lines = nt
      !!x print *, "B - nt = ",nt
      !
      ! Order the REC lines in ascending order of photon energies
      if (rec_sort_line_energy) then
         do  i = 1,number_of_lines-1
            imin       = i
            energy_min = 1.0e20
            do  j = i,number_of_lines
               if (line(j)%p_energy < energy_min) then
                  imin = j;   energy_min = line(j)%p_energy
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
      ! Now assign the full line specifications 
      do  nt = 1,number_of_lines
      !!x print *, "C - nt = ",nt
         i = line(nt)%asfi;   j = line(nt)%asff
         line(nt)%level_i   = asf_initial%asf(i)%level_No
         line(nt)%level_f   = asf_final%asf(j)%level_No
         line(nt)%totalJ_i  = asf_initial%asf(i)%totalJ
         line(nt)%totalJ_f  = asf_final%asf(j)%totalJ
         line(nt)%parity_i  = asf_initial%asf(i)%parity
         line(nt)%parity_f  = asf_final%asf(j)%parity
         line(nt)%beta_b(:) = zero
         line(nt)%beta_c(:) = zero
         !
         !
         ! Modify the transition energies if required
         line(nt)%p_energy = line(nt)%p_energy + rec_energy_shift
         !
         nch = 0
         do  ka = 1,rec_maximal_kappa
            do  kas = -1,1,2
               kappa  = kas * ka
               j      = angular_momentum_j(kappa)
               lkappa = angular_momentum_l(kappa)
               !
               if ((line(nt)%parity_i == "+"  .and. mod(lkappa,2) == 0)  .or. &
                   (line(nt)%parity_i == "-"  .and. mod(lkappa,2) == 1)) then
                  parity = "+"
               else 
                  parity = "-"
               end if
               !
               ! Cycle over all multipoles
               do L = 1,number_of_multipoles
                  do  p = 1,0,-1
                     if ((parity == "+"  .and.  line(nt)%parity_f == "+") .or.&
                         (parity == "-"  .and.  line(nt)%parity_f == "-")) then
                        if (mod(L,2) == 1) then
                           if (p == 1) cycle
                        else
                           if (p == 0) cycle
                        end if
                     else
                        if (mod(L,2) == 1) then
                           if (p == 0) cycle
                        else
                           if (p == 1) cycle
                        end if
                     end if
                     !
                     do  totalJ = line(nt)%totalJ_i-j,line(nt)%totalJ_i+j,2
                        if (totalJ < 0  .or.                                   &
                            .not.is_triangle(totalJ,L+L,line(nt)%totalJ_f) .or.&
                            .not.is_triangle(line(nt)%totalJ_i,j,totalJ))  cycle
                        !
                        nch = nch + 1
                        wa_kappa(nch)        = kappa
                        wa_totalJ(nch)       = totalJ
                        wa_parity(nch)       = parity
                        wa_multipole_L(nch)  = L
                        wa_multipole_pi(nch) = p
                        if (p == 0) then
                           ! Magnetic case
                           wa_gauge(nch)     = "Magnetic "
                        else
                           wa_gauge(nch)        = "Babushkin"
                           nch                  = nch + 1
                           wa_kappa(nch)        = kappa
                           wa_totalJ(nch)       = totalJ
                           wa_parity(nch)       = parity
                           wa_multipole_L(nch)  = L
                           wa_multipole_pi(nch) = p
                           wa_gauge(nch)        = "Coulomb  "
                        end if
                     end do
                  end do
               end do
               !
            end do
         end do
         !
         line(nt)%No_channels = nch
         allocate( line(nt)%channel(1:nch) )
         !
         line(nt)%channel(1:nch)%kappa        = wa_kappa(1:nch)
         line(nt)%channel(1:nch)%totalJ       = wa_totalJ(1:nch)
         line(nt)%channel(1:nch)%parity       = wa_parity(1:nch)
         line(nt)%channel(1:nch)%multipole_L  = wa_multipole_L(1:nch)
         line(nt)%channel(1:nch)%multipole_pi = wa_multipole_pi(1:nch)
         line(nt)%channel(1:nch)%gauge        = wa_gauge(1:nch)
         line(nt)%channel(1:nch)%phase        = zero
         line(nt)%channel(1:nch)%amplitude    = cmplx(zero,zero)
      end do
      !
      print *, "  ",number_of_lines,                             &
               " REC lines have been initialized and will be "// &
               "calculated in this run of the program."
      !
   end subroutine rec_set_lines
   !
   !
   subroutine rec_set_overlaps(stream)
   !--------------------------------------------------------------------
   ! Calculates the non-orthogonal and overlap integrals between all
   ! bound orbitals of the corresponding final- and initial-state arrays.
   !
   ! Calls: rk_integral_grasp2k_ab().
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream
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
      allocate( rec_overlap(kappa_min:kappa_max,pqn_f_max,pqn_i_max) )
      !
      ! Calculate the non-orthogonal and overlap integrals
      do  i1 = 1,wave_final%number_of_rwf
         do  i2 = 1,wave_initial%number_of_rwf
            if (wave_final%rwf(i1)%orbital%kappa == &
                wave_initial%rwf(i2)%orbital%kappa) then
               kappa = wave_final%rwf(i1)%orbital%kappa
               pqn_f = wave_final%rwf(i1)%orbital%n
               pqn_i = wave_initial%rwf(i2)%orbital%n
               rec_overlap(kappa,pqn_f,pqn_i) = &
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
               !
	       !! write(stream,*) "  "//orbital_symmetry(kappa)//               &
	       !!                 ":  < n_f = ",pqn_f," | n_i = ",pqn_i,"> = ", &
	       !!                 rec_overlap(kappa,pqn_f,pqn_i)
            end if
         end do
      end do
      !
   end subroutine rec_set_overlaps
   !
   !
   function rec_statistical_tensor(G,k,tline)               result(rho)
   !--------------------------------------------------------------------
   ! Calculates the REC statistical tensor rho_k (alpha_i J_i, alpha_f J_f) 
   ! from the many-electron REC amplitudes.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      character(len=1), intent(in) :: G
      integer, intent(in)          :: k
      type(rec_line), intent(in)   :: tline
      complex(kind=dp)             :: rho
      !
      integer          :: i, ip, J, Jp, jk, jkp, J_i, J_f,  &
                          L, Lp, lk, lkp, p, pp
      complex(kind=dp) :: N, me, mep
      !
      J_i = tline%totalJ_i
      J_f = tline%totalJ_f
      !
      N = zero
      ! 
      do  i = 1,tline%No_channels
         if (tline%channel(i)%gauge(1:1)  == G   .or. &
             tline%channel(i)%gauge(1:1)  == "M") then
            me = tline%channel(i)%amplitude
         else
            cycle
         end if
         !
         J   = tline%channel(i)%totalJ
         L   = tline%channel(i)%multipole_L
         p   = tline%channel(i)%multipole_pi
         lk  = angular_momentum_l(tline%channel(i)%kappa)   
         jk  = angular_momentum_j(tline%channel(i)%kappa) 
         !  
         do  ip = 1,tline%No_channels
            if (tline%channel(ip)%gauge(1:1) == G   .or. &
                tline%channel(ip)%gauge(1:1) == "M") then
               mep = tline%channel(ip)%amplitude
            else
               cycle
            end if
            !
            Jp  = tline%channel(ip)%totalJ
            Lp  = tline%channel(ip)%multipole_L  
            pp  = tline%channel(ip)%multipole_pi
            lkp = angular_momentum_l(tline%channel(ip)%kappa)   
            jkp = angular_momentum_j(tline%channel(ip)%kappa)
            !
            if (L == Lp   .and.   p == pp) then  
               N  = N  + ((-1)**((J_i+L+L-J_f+J-Jp-1)/2))                    & 
                  * sqrt( (lk+lk+one)*(lkp+lkp+one)*(jk+one)*(jkp+one)       &
                         *(J+one)*(Jp+one) )                                 &
                  * Clebsch_Gordan(lk+lk,0,lkp+lkp,0,k+k,0)                  &
                  * wigner_6j_symbol(jk,jkp,k+k,lkp+lkp,lk+lk,1)             &
                  * wigner_6j_symbol(jk,jkp,k+k,Jp,J,J_i)                    &
                  * wigner_6j_symbol(J,Jp,k+k,J_f,J_f,L+L)                   &
                  * me * conjg(mep)
             end if
         end do
      end do
      !
      rho = N
      !
   end function rec_statistical_tensor
   !
end module rabs_rec
