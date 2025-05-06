module rabs_einstein
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains the procedures which are specific to the EINSTEIN
! program for the computation of transition probabilities and lifetimes.
! Unlike to the REOS component, however, orthogonality of the initial and
! final state orbitals is assumed in the evaluation of the (many-electron)
! matrix elements and, hence, Racah's algebra is utilized to 'divide' the
! computations into angular coefficients (from ANCO) and radial integrals.
! Since the initial and final states are considered to be independent,
! some of the non-orthogonality can be taken into accout in the radial
! integrals.
!
! The EINSTEIN is expected to be much faster than the REOS program and
! has been developed for mass production of transition probability data 
! as well as for studies on the angular distributions of the emitted photon
! for which the determinant approach in REOS may become to cumbersome.
!-----------------------------------------------------------------------
   !
   use rabs_anco
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_grasp2k
   use rabs_input_dialog
   use rabs_mcp
   use rabs_multipole
   use rabs_nonorthonormal
   use rabs_print
   implicit none
   private
   !
   public  :: einstein_calculate_amplitudes
                 ! Calculates the multipole amplitudes for all the selected 
                 ! transitions.
   public  :: einstein_collect_input
                 ! Collects and proceeds all input for the calculation of
		 ! multipole amplitudes, transition probabilities, line
                 ! strengths, and angular distributions.
   public  :: einstein_initialize
                 ! Set up the selected transitions and initializes some
		 ! arrays as required.
   private :: einstein_multipole_amplitude
                 ! Calculates a single multipole amplitude from the 'pure'
                 ! multipole matrix and the corresponding mixing coefficients.
   public  :: einstein_print_results
                 ! Writes the transition probabilities, lifetimes, and 
                 ! angular distribution parameters to the .sum file.
   private :: einstein_print_transitions
                 ! Prints all selected transitions in a neat format before
                 ! the computation starts.
   private :: einstein_pure_matrix
                 ! Calculates the 'pure' multipole matrix for the given
		 ! configuration scheme.
   private :: einstein_set_transitions
                 ! Determines which radiative transitions need to be calculated
                 ! and initializes storage for an appropriate data type for
                 ! them.
   private :: einstein_transition_properties
                 ! Calculates all selected properties of the transition i 
		 ! from the amplitudes of the individual multipole components.
   public  :: einstein_write_trn_file
                 ! Writes out transition energies and amplitudes for further
                 ! processing to a .trn file.
   !
   ! Define some global data of the EINSTEIN program; most of these data are
   ! read in during the interactive control and may overwrite existing
   ! default values
   !
   ! Storage for the initial and final atomic states and wave functions
   type(asf_basis), public       :: asf_bound
   !
   ! Define an internal structure type(einstein_transition) which stores all
   ! necessary information for a multipole transition
   type :: multipole_line
      character(len=2) :: multipole  
      character(len=9) :: gauge 
      real(kind=dp)    :: amplitude   
   end type multipole_line
   !
   type :: einstein_transition
      integer          :: asfi, asff, level_i, level_f, totalJ_i, totalJ_f
      integer          :: number_of_mlines
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: energy, einstein_a
      type(multipole_line), dimension(:), pointer :: mline
   end type einstein_transition
   !
   type(einstein_transition), dimension(:), allocatable :: transition
   !
   type :: einstein_matrix
      integer :: No_f, No_i  
      integer, dimension(:), pointer         :: ndx_f, ndx_i  
      real(kind=dp), dimension(:,:), pointer :: matrix   
   end type einstein_matrix
   !
   type(einstein_matrix) :: einstein
   !
   integer          :: number_of_transitions                 = 0
   !
   ! Define global logical flags for the control of the PHOTO program; the
   ! default values for these flags may be overwritten interactively during 
   ! input time
   logical, public :: einstein_apply_exp_energies        = .false.,  &
                      einstein_calc_angular_parameter    = .false.,  &
		      einstein_nonorthogonal_eval        = .false.,  &
		      einstein_print_csf_scheme          = .false.,  &
		      einstein_print_each_line           = .false.,  &
                      einstein_print_only_gt_1percent    = .false.,  &
                      einstein_print_selected_trans      = .false.,  &     
                      einstein_print_trn_file            = .false.,  &
                      einstein_sort_transition_energy    = .true.,   &
                      einstein_use_formatted_mix_file    = .true.,   &
                      einstein_use_formatted_rwf_file    = .true.,   &
                      einstein_write_transition_file
   !
   ! Energy unit for the output of all energies
   real(kind=dp)    :: einstein_rate_factor      = zero,             &
                       einstein_energy_shift     = zero
   !
   character(len=7) :: einstein_rate_unit
   !
   ! Define some variables and arrays for processing input data from 
   ! einstein_collect_input()
   !
   !
contains
   !
   subroutine einstein_calculate_amplitudes()
   !--------------------------------------------------------------------
   ! Calculates in turn the transition amplitudes for all the multipole
   ! components of the selected radiative transitions.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer       :: i, j, k, n, nw, nu, nocsf, par
      real(kind=dp) :: energy
      type(nkappa)  :: subshell
      integer, dimension(:), allocatable :: ndx
      !
      n = asf_final%csf_set%nocsf + asf_initial%csf_set%nocsf 
      allocate( ndx(1:n) )
      !
      do  i = 1,number_of_transitions
         do  j = 1,transition(i)%number_of_mlines
            energy = transition(i)%energy
            asf_bound%csf_set%nocsf = 0
            !
            asf_bound%csf_set%nwshells = asf_final%csf_set%nwshells
            asf_bound%csf_set%nwcore   = asf_final%csf_set%nwcore
            asf_bound%csf_set%number_of_electrons   =                         &
                                          asf_final%csf_set%number_of_electrons
            allocate( asf_bound%csf_set%subshell(asf_final%csf_set%nwshells) )
            do  k = 1,asf_final%csf_set%nwshells 
               asf_bound%csf_set%subshell(k) = asf_final%csf_set%subshell(k)
            end do
            !
            nw = asf_final%csf_set%nwshells
            allocate( asf_bound%csf_set%csf(1:1) )
	    allocate( asf_bound%csf_set%csf(1)%occupation(1:nw),              &
	              asf_bound%csf_set%csf(1)%seniority(1:nw),               &
	              asf_bound%csf_set%csf(1)%subshellJ(1:nw),               &
		      asf_bound%csf_set%csf(1)%subshellX(1:nw) )
            !
            call add_csf_to_basis(asf_final%csf_set,asf_bound%csf_set,        &
                     transition(i)%totalJ_f,transition(i)%parity_f,index=ndx)
            !
            einstein%No_f = asf_bound%csf_set%nocsf
            allocate( einstein%ndx_f(einstein%No_f) )
            einstein%ndx_f(1:einstein%No_f) = ndx(1:einstein%No_f)
            !
            nw = asf_bound%csf_set%nwshells
            if (rabs_use_stop  .and. nw /= asf_initial%csf_set%nwshells) then
               stop "einstein_calculate_amplitudes(): program stop A."
            end if
            !
            ! Define the 'extended' configuration scheme for calculating
            ! the transition probabilities and allocate the memory
            call add_csf_to_basis(asf_initial%csf_set,asf_bound%csf_set,      &
                     transition(i)%totalJ_i,transition(i)%parity_i,index=ndx)
            if (einstein_print_csf_scheme) then
               call print_configuration_scheme(6,asf_bound%csf_set)
            end if
            !
            einstein%No_i = asf_bound%csf_set%nocsf - einstein%No_f
            allocate( einstein%ndx_i(einstein%No_i) )
            einstein%ndx_i(1:einstein%No_i) =                                 &
                                   ndx(1+einstein%No_f:asf_bound%csf_set%nocsf)
            allocate( einstein%matrix(1:einstein%No_f,1:einstein%No_i) )
            !!x print *, "einstein%No_i, einstein%ndx_i(k) = ",               &
            !!x           einstein%No_i, (einstein%ndx_i(k),k=1,einstein%No_i) 
            !
            ! Calculate the 'pure' multipole matrix in the given CSF scheme
            ! (not including mixing coefficients)
            call einstein_pure_matrix(i,j,asf_bound%csf_set)
            !
            call einstein_multipole_amplitude(i,j)
            !
            deallocate( einstein%ndx_f, einstein%ndx_i, einstein%matrix  )
         end do
         !
         ! Calculates all selected properties for the selected transition
         call einstein_transition_properties(transition(i))
         !
      end do
      deallocate( ndx )
      !
   end subroutine einstein_calculate_amplitudes
   !
   !
   subroutine einstein_collect_input()
   !--------------------------------------------------------------------
   ! Collect and proceeds all input for the EINSTEIN program.
   !
   ! Calls: get_yes_stream().
   !--------------------------------------------------------------------
      !
      integer            :: level_i, level_f, score_position
      logical            :: yes
      character(len=20 ) :: string
      character(len=256) :: record
      !
      !
      ! Open, check, load data from, and close the  .iso  file
      call file_get_isodat_grasp2k()
      !
      ! Determine the transition multipoles
      call input_transition_multipoles(number_of_multipoles,multipoles)    
      !
      ! Determine the units for the printout and further optional input
      call input_energy_unit()
      !
      ! Determine the parameters controlling the radial grid
      call input_grid_parameters("standard")
      !
      hp_grasp2k = zero
      !
      ! Default speed of light
      c = c_vacuum
      !
      ! Now 'overwrite' defaults only if required
      print *, "Modify default set-up and printout of the program ?"
      yes = get_yes_stream()
      if (.not.yes) goto 6
      !
      ! Select individual pairs of transitions 
      call input_transition_pairs(number_of_transitions)
      !
      einstein_nonorthogonal_eval = .false.
      !
      print *, "Sort transitions in ascending order of energy ?"
      yes = get_yes_stream()
      if (yes) einstein_sort_transition_energy = .true.
      !
      print *, "Read in and apply experimental energies for the calculation"//&
               " of transition probabilities ?"
      einstein_apply_exp_energies = get_yes_stream()
      !
      print *, "Einstein A and B coefficients are printed in SI units;"
      print *, " use Hartree atomic units instead ?"
      multipole_print_AB_in_hartree = get_yes_stream()
      !
      print *, "Print all selected transitions and their energies"// &
               " before the computation starts (this is the default) ?"
      einstein_print_selected_trans = get_yes_stream()
      !
      print *, "Print the results for each individual line immediatly"//&
               " after its computation ?"
      einstein_print_each_line = get_yes_stream()
      !
      print *, "Write out the transition energies and amplitudes to an"// &
               " .trn file for further data processing,"
      print *, " for instance, to adopt them to experimental transition"//&
               " energies ?"
      yes = get_yes_stream()
      if (yes) then
         einstein_write_transition_file = .true.
	 !
         call file_open_formatted_stream(27, &
	    "Enter a file name for the  einstein.trn  file:")
         write(27, "(a)") "EINSTEIN Transition energy and amplitude file"
         write(27,*)
      end if
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
    6 continue
      !
      if (multipole_print_AB_in_hartree) then
	 einstein_rate_factor = one
	 einstein_rate_unit   = "a.u.	"
      else
	 einstein_rate_factor = convert_einstein_a_to_si
	 einstein_rate_unit   = "1/s	"
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
    7 continue
      !
      !
   end subroutine einstein_collect_input
   !
   !
   subroutine einstein_initialize()
   !--------------------------------------------------------------------
   ! Initializes the computation of radiative transition probabilities,
   ! Einstein coefficients, lifetimes and others. 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer :: i, mtp, noasf, nocsf, number_of_rwf
      !
      print *, "Initialize the set-up of the transitions probabilities ..."
      call einstein_set_transitions()
      !
      print *, "   ... initialization complete."
      print *, " "
      !
      ! Print the selected transitions before the computation starts
      if (einstein_print_selected_trans) then
         call einstein_print_transitions(6)
         call einstein_print_transitions(24)
      end if
      !
   end subroutine einstein_initialize
   !
   !
   subroutine einstein_multipole_amplitude(i,j)
   !--------------------------------------------------------------------
   ! Calculates the transition amplitude for the transition i and multipole
   ! channel j by the summation over the 'pure' multipole matrix and by
   ! using the proper weights for the transition i.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: i, j
      !
      integer       :: asfi, asff, l, r, rr, s, ss
      real(kind=dp) :: phase, value
      !
      asfi  = transition(i)%asfi;  asff = transition(i)%asff
      value = zero
      do  r = 1,einstein%no_f
         rr = einstein%ndx_f(r)
         do  s = 1,einstein%no_i
            ss = einstein%ndx_i(s)
            value = value + asf_final%asf(asff)%eigenvector(rr) * &
                    einstein%matrix(r,s) * asf_initial%asf(asfi)%eigenvector(ss)
            !
         end do
      end do
      !
      transition(i)%mline(j)%amplitude = value
      !
      !!x print *, "i,j,transition(i)%mline(j)%amplitude = ",  &
      !!x           i,j,transition(i)%mline(j)%amplitude
      !
   end subroutine einstein_multipole_amplitude
   !
   !
   subroutine einstein_print_results(stream)
   !--------------------------------------------------------------------
   ! Writes the transition probabilities, Einstein coefficients, and other
   ! parameters in a neat summary to stream.
   !
   ! Calls:  angular_momentum_string().
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      !
      integer           :: i, j, istate_low, istate_up, nm, rank
      character(len=12) :: multipoles
      real(kind=dp)     :: einstein_A, einstein_B, oscillator, decay_width,    &
                           energy_au, energy, tb, tb_au, tb_cm, tb_ev, tb_sec, &
                           tc, tc_au, tc_cm, tc_ev, tc_sec, f2, f4
      !
      real(kind=dp), dimension(:), allocatable :: total_babushkin,total_coulomb
      logical, dimension(:), allocatable       :: pure_magnetic  
      !
      istate_low = 10000;   istate_up = -10000
      do  i = 1,number_of_transitions
         istate_low = min(istate_low, transition(i)%level_i )
         istate_up  = max(istate_up,  transition(i)%level_i )
      end do
      allocate( total_babushkin(istate_low:istate_up), &
                total_coulomb(istate_low:istate_up),   &
                pure_magnetic(istate_low:istate_up)    )
      !
      do i = istate_low,istate_up
         total_babushkin(i) = zero;   total_coulomb(i) = zero
         pure_magnetic(i) = .true.
      end do
      !
      do i = 1,number_of_transitions
         do j = 1,transition(i)%number_of_mlines
            if (transition(i)%mline(j)%gauge == "Babushkin"   .or. &
                transition(i)%mline(j)%gauge == "Coulomb  ") then
               pure_magnetic(transition(i)%level_i) = .false.
            end if
         end do
      end do
      !
      write(stream,1)
    1 format(/ &
         /32x,"===========================================================", &
         /32x,"|  Summary of all Transition Probabilities and Lifetimes  |", &
         /32x,"==========================================================="  )
      !
      if (energy_unit == "A      "  .and.  &
          multipole_print_AB_in_hartree) then
         write(stream,2)
         write(stream,4) 
      else if (energy_unit == "A      ") then
         write(stream,3)
         write(stream,5) 
      else if (multipole_print_AB_in_hartree) then
         write(stream,2)
         write(stream,6) trim(energy_unit)
      else 
         write(stream,3)
         write(stream,7) trim(energy_unit)
      endif
    2 format(// 1x,122("-"),                                                 &
              / 2x,"LevI-LevF  I- J / Parity -F      Energy   ",             &
                   "Multipol   Gauge         Einstein coefficients",         &
                   "       Oscillator    Decay width  ",                     &
              /71x,"                     " )
    3 format(// 1x,122("-"),                                                 &
              / 2x,"LevI-LevF  I- J / Parity -F      Energy   ",             &
                   "Multipol   Gauge         Einstein coefficients",         &
                   "       Oscillator    Decay width  ",                     &
              /71x,"-1           3 -2 -1 " )   
    4 format(  30x,"   (Angstroms)",                                         &
               23x,"A (a.u.)    gB (a.u.)        strength GF       (eV) ")
    5 format(  30x,"   (Angstroms)",                                         &
               23x,"A (s  )     gB (m s  J  )    strength GF       (eV) ")
    6 format(  30x,"     (",a4,")   ",                                       &
               23x,"A (a.u.)    gB (a.u.)        strength GF       (eV) ")
    7 format(  30x,"     (",a4,")   ",                                       &
               23x,"A (s  )     gB (m s  J  )    strength GF       (eV) ")
      !
      write(stream,8)
    8 format(  /1x,122('-') )
      !
      do i = 1,number_of_transitions
         do j = 1,transition(i)%number_of_mlines
            energy_au = transition(i)%energy
            if (energy_inverse) then
               energy = energy_factor / transition(i)%energy
            else
               energy = energy_factor * transition(i)%energy
            end if
            select case(transition(i)%mline(j)%multipole)
            case("E1", "M1");   rank = 1
            case("E2", "M2");   rank = 2
            case("E3", "M3");   rank = 3
            case("E4", "M4");   rank = 4
            case("E5", "M5");   rank = 5
            case default; stop "einstein_print_results(): program stop A."
            end select
            !
            call multipole_convert_probability(transition(i)%totalJ_i,   &
                                               transition(i)%totalJ_f,   &
                                               rank,energy_au,           &
                             transition(i)%mline(j)%amplitude,           &
                             einstein_A,einstein_B,oscillator,decay_width)
            !
            write(stream,9) transition(i)%level_i,transition(i)%level_f, &
               trim(angular_momentum_string(transition(i)%totalJ_i,4)),  &
               transition(i)%parity_i,                                   &
               trim(angular_momentum_string(transition(i)%totalJ_f,4)),  &
               transition(i)%parity_f,energy,                            &
               transition(i)%mline(j)%multipole,                         &
               transition(i)%mline(j)%gauge,einstein_A,einstein_B,       &
               oscillator,decay_width
            !
            if (transition(i)%mline(j)%gauge == "Babushkin"  .and.       &
                energy_au > zero) then
               total_babushkin(transition(i)%level_i) =                  &
                  total_babushkin(transition(i)%level_i) + einstein_A
            else if (transition(i)%mline(j)%gauge == "Coulomb  "  .and.  &
                energy_au > zero) then
               total_coulomb(transition(i)%level_i) =                    &
                  total_coulomb(transition(i)%level_i)   + einstein_A
            else if (transition(i)%mline(j)%gauge == "Magnetic "  .and.  &
                energy_au > zero) then
               total_babushkin(transition(i)%level_i) =                  &
                  total_babushkin(transition(i)%level_i) + einstein_A
               total_coulomb(transition(i)%level_i) =                    &
                  total_coulomb(transition(i)%level_i)   + einstein_A
            else if (rabs_use_stop) then
               !! stop "einstein_print_results(): program stop B."
            end if
         end do
      end do
    9 format(2x,i3," -",i3,3x,a4,1x,a1,4x,a4,1x,a1,3x,1pd12.5,3x,a2,     &
             4x,a9,3x,1pd12.5,3x,1pd12.5,3x,1pd12.5,3x,1pd12.5)
             !
      write(stream,8)
      !
      !
      ! Print lifetimes and width of levels
      !
      write(stream,*)
      write(stream,10)
   10 format(//"Radiative lifetimes and widths"                          &
              /"------------------------------"                          &
            ///" LeveL",6x,"Gauge",15x,"Lifetime",33x,"Width"            &
              /" -----",6x,"-----",15x,"--------",12x,                   &
               "---------------------------------------------",          &
          /32x,"Seconds",13x,"Hartrees",12x,"Kaysers",16x,"eV"/)
      !
      do  i = istate_low,istate_up
         if( total_babushkin(i) == zero  .and.  total_coulomb(i) == zero) then
            cycle
         end if 
         tb = total_Babushkin(i)
         tc = total_Coulomb(i)
         if (.not.multipole_print_AB_in_hartree) then
            tb_cm  = tb / c_vacuum_in_cm_per_s 
            tc_cm  = tc / c_vacuum_in_cm_per_s 
            tc_au  = tc_cm / convert_au_to_kaysers
            tb_au  = tb_cm / convert_au_to_kaysers
            tc_sec = one / tc
            tb_sec = one / tb
            tc_ev  = tc_au * convert_au_to_ev
            tb_ev  = tb_au * convert_au_to_ev
         else
            tc_au  = tc
            tb_au  = tb
            tc_cm  = tc * convert_au_to_kaysers
            tb_cm  = tb * convert_au_to_kaysers
            tc_sec = one / ( tc_cm * c_vacuum_in_cm_per_s )
            tb_sec = one / ( tb_cm * c_vacuum_in_cm_per_s )
            tc_ev  = tc_au * convert_au_to_ev
            tb_ev  = tb_au * convert_au_to_ev
         end if
         !
         if (pure_magnetic(i)) then
            write(stream,11) i,tc_sec,tc_au,tc_cm,tc_ev
         else
            write(stream,12) i,tb_sec,tb_au,tb_cm,tb_ev
            write(stream,13)   tc_sec,tc_au,tc_cm,tc_ev
         end if
      11 format(1x,i4,6x,"Magnetic:  ",4(1pd20.7)/) 
      12 format(1x,i4,6x,"Babushkin: ",4(1pd20.7))
      13 format(     11x,"Coulomb:   ",4(1pd20.7)/)
         !
      end do
      deallocate( total_babushkin, total_coulomb, pure_magnetic )
      !
      !
   21 format(//"Structure functions f_2() and f_4()"                          &
              /"-----------------------------------" )
   22 format( / 1x,108("-"),                                                  &
              / 2x,"LevI-LevF  I- J / Parity -F      Energy      ",           &
                   "Multipoles      Gauge           f_2            f_4" )
   24 format(  30x,"   (Angstroms)  ")
   26 format(  30x,"     (",a4,")   ")
   28 format(   1x,108('-') )
   29 format(2x,i3," -",i3,3x,a4,1x,a1,4x,a4,1x,a1,3x,1pd12.5,3x,a12,         &
             4x,a9,3x,1pd12.5,3x,1pd12.5)
      !
   end subroutine einstein_print_results
   !
   !
   subroutine einstein_print_transitions(stream)
   !--------------------------------------------------------------------
   ! Prints a neat table of all selected transitions on stream before 
   ! the actual computation starts; only the quantum numbers of the atomic 
   ! states and the transition energies are displayed.
   !
   ! Calls: angular_momentum_string().
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      integer		  :: i, j, nmult
      real(kind=dp)	  :: energy
      character(len=2)    :: mult(20)
      !
      write(stream,1) number_of_transitions,trim(energy_unit)
      do  i = 1,number_of_transitions
	 if (energy_inverse) then
            energy = energy_factor / transition(i)%energy
         else
            energy = energy_factor * transition(i)%energy
         end if
         nmult = 1; mult(1) = transition(i)%mline(1)%multipole
         do j = 2,transition(i)%number_of_mlines
            if (mult(nmult) /= transition(i)%mline(j)%multipole) then
               nmult = nmult + 1
               mult(nmult) = transition(i)%mline(j)%multipole
            end if
         end do
         write(stream,2) transition(i)%level_i,transition(i)%level_f,        &
                   trim(angular_momentum_string(transition(i)%totalJ_i,4)),  &
                         transition(i)%parity_i,                             &
                   trim(angular_momentum_string(transition(i)%totalJ_f,4)),  &
                         transition(i)%parity_f,energy,(mult(j),j=1,nmult)
      end do
      write(stream,3) 
      !
    1 format( "The following ",i5," transitions are selected:",              &
        //,"     I-level-F     I--J^P--F      Transition Energy       ",     &
           "Multipoles ",                                                    &
         /,"                                     (in ",a4,")  ",             &  
         /,4x,66("-") )
    2 format(4x,i4," -",i4,3x,a4,a1,3x,a4,a1,5x,1pe14.7,9x,10(a2,1x))
    3 format(4x,66("-") )
      !
   end subroutine einstein_print_transitions
   !
   !
   subroutine einstein_pure_matrix(tr,ch,csf_set)
   !--------------------------------------------------------------------
   ! Calculates the 'pure' multipole matrix for the given configuration 
   ! scheme csf_set. The first no_f CSF belong to the final-state 
   ! representation and the following no_f+1,...,no_f+no_i to the 
   ! initial states.
   !
   ! Calls:  
   !--------------------------------------------------------------------
      !
      integer, intent(in)         :: tr, ch
      type(csf_basis), intent(in) :: csf_set
      !
      integer                     :: i, ia, ib, j, ja, jb, no_T_coeff, nu, &
                                     par, r, s, ss, t, parity      
      type(nkappa)                :: aa, bb, cc, dd
      real(kind=dp)               :: aweight
      type(nkappa)                :: Ta, Tb
      !
      select case(transition(tr)%mline(ch)%multipole(1:2))
      case("E1");   nu = 1;   par = 1;   parity = -1
      case("M1");   nu = 1;   par = 0;   parity =  1
      case("E2");   nu = 2;   par = 1;   parity =  1
      case("M2");   nu = 2;   par = 0;   parity = -1
      case("E3");   nu = 3;   par = 1;   parity = -1
      case("M3");   nu = 3;   par = 0;   parity =  1
      case("E4");   nu = 4;   par = 1;   parity =  1
      case("M4");   nu = 4;   par = 0;   parity = -1
      case default
         stop "einstein_pure_matrix(): program stop A."
      end select
      !
      if (rabs_use_stop   .and.                            &
          csf_set%nocsf /= einstein%no_f+einstein%no_i) then
         stop "einstein_pure_matrix(): program stop B."
      end if
      !
      einstein%matrix = zero
      !
      if (einstein_nonorthogonal_eval) then
         call nonorth_initialize(csf_set,wave_final,wave_initial)
      end if
      !
      do  r = 1,einstein%No_f
         do  s = einstein%No_f+1,einstein%No_f+einstein%No_i
            ss = s - einstein%No_f
	    if (einstein_nonorthogonal_eval) then
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
	       if (einstein_nonorthogonal_eval) then
	          Ta      = nonorth_T_list(t)%a
	          Tb      = nonorth_T_list(t)%b
                  ja      = angular_momentum_j(Ta%kappa)
                  aweight = nonorth_T_list(t)%T * sqrt( ja+one )
		  !!x  * sqrt(transition(tr)%totalJ_f+one) 
		  !
                  ia = 0;   ib = 0
                  do  i = 1,wave_final%number_of_rwf
                     if (wave_final%rwf(i)%orbital == Ta) then
                        ia = i
                        exit
                     end if
                  end do
                  do  i = 1,wave_initial%number_of_rwf
                     if (wave_initial%rwf(i)%orbital == Tb) then
                        ib = i
                        exit
                     end if
                  end do
		  !
	       else
	          ia       = mct_list(t)%a
	          ib       = mct_list(t)%b
	          aweight  = mct_list(t)%T
		  Ta       = wave_final%rwf(ia)%orbital
		  Tb       = wave_initial%rwf(ib)%orbital
		  !!x print *, "Orth: ia, ib, Ta, Tb, aweight = ", &
		  !!x                 ia, ib, Ta, Tb, aweight
	       end if
               !
               if (ia /= 0  .and.  ib /= 0) then
                  einstein%matrix(r,ss) = einstein%matrix(r,ss) + aweight *   &
                  multipole_reduced_M_integral(abs(transition(tr)%energy),    &
                                       nu,par,transition(tr)%mline(ch)%gauge, &
                                       wave_final%rwf(ia),wave_initial%rwf(ib))
               else
                  stop "einstein_pure_matrix(): program stop C."
               end if
               !
            end do
         end do
      end do
      !
   end subroutine einstein_pure_matrix
   !
   !
   subroutine einstein_set_transitions()
   !--------------------------------------------------------------------
   ! Determines how many and which transitions need to be calculated and 
   ! initializes the array transitions of type(reos_transitions).
   ! The default is (for number_of_transitions == 0) that all transitions
   ! with a positive energy are calculated; for number_of_transitions /= 0,
   ! individual transitions were selected during input time and will be
   ! initialized instead, in this case also negative transition energies
   ! (absorption lines) are allowed.
   !
   ! Calls: multipole_select().
   !--------------------------------------------------------------------
      !
      integer       :: asfi, asff, i, imin, j, k, m, nmult, nt
      real(kind=dp) :: arg, energy, energy_exp, energy_min
      character(len=2), dimension(20)       :: mult
      character(len=9), dimension(20)       :: gauge
      !
      ! Determine the total number of transitions
      nt = 0
      if (number_of_transitions == 0) then
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               if (asf_initial%asf(i)%energy - asf_final%asf(j)%energy &
                                                                 > zero) then
                  call multipole_select(1*asf_initial%asf(i)%totalJ,           &
                                                asf_initial%asf(i)%parity,     &
                        1*asf_final%asf(j)%totalJ,asf_final%asf(j)%parity,     &
                        mult,gauge,nmult)
                  if (nmult == 0) cycle
                  nt = nt + 1
               end if
            end do
         end do
      else
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_transitions
                  if ( (select_level_i(k) == asf_initial%asf(i)%level_No .or.  &
                        select_level_i(k) == 0)                         .and.  &
                       (select_level_f(k) == asf_final%asf(j)%level_No   .or.  &
                        select_level_f(k) == 0) ) then 
                     call multipole_select(1*asf_initial%asf(i)%totalJ,        &
                                                   asf_initial%asf(i)%parity,  &
                           1*asf_final%asf(j)%totalJ,asf_final%asf(j)%parity,  &
                           mult,gauge,nmult)
                     if (nmult == 0) then
                        print *, "Transition ",asf_initial%asf(i)%level_No,"-",&
                                 asf_final%asf(j)%level_No,"was selected at ", &
                                 "input time but does not allow any of ",      &
                                 "the given multipoles; "
                        print *, "this transition is therefore not considered."
                        cycle
                     end if
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
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               if (asf_initial%asf(i)%energy - asf_final%asf(j)%energy         &
                                                                    > zero) then
                  call multipole_select(1*asf_initial%asf(i)%totalJ,           &
                                                asf_initial%asf(i)%parity,     &
                        1*asf_final%asf(j)%totalJ,asf_final%asf(j)%parity,     &
                        mult,gauge,nmult)
                  if (nmult == 0) cycle
                  nt = nt + 1
                  transition(nt)%asfi     = i;   transition(nt)%asff = j
                  transition(nt)%level_i  = asf_initial%asf(i)%level_No
                  transition(nt)%level_f  = asf_final%asf(j)%level_No
                  transition(nt)%energy   = asf_initial%asf(i)%energy -        &
                                            asf_final%asf(j)%energy
                  !
                  ! Read in experimental energies if appropriate
                  if (einstein_apply_exp_energies) then
                     if (energy_inverse) then
                        energy = energy_factor / transition(nt)%energy
                     else
                        energy = energy_factor * transition(nt)%energy
                     end if
                  99 print *, "Transition",transition(nt)%level_i,"-",         &
                      transition(nt)%level_f,"has ab-initio energy E_theo = ", &
                      energy,trim(energy_unit),";"
                     print *, " enter E_exp (in ",trim(energy_unit),  &
                      ") or 0.0 to use E_theo:"
                     read(*,*,err=99) energy_exp
                     if (energy_exp /= 0) then
                        if (energy_inverse) then
                           transition(nt)%energy = energy_factor / energy_exp
                        else
                           transition(nt)%energy = energy_exp / energy_factor
                        end if
                     end if
                  end if
               end if
            end do
         end do
      else
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               do  k = 1,number_of_transitions
                  if ( (select_level_i(k) == asf_initial%asf(i)%level_No .or.  &
                        select_level_i(k) == 0)                         .and.  &
                       (select_level_f(k) == asf_final%asf(j)%level_No   .or.  &
                        select_level_f(k) == 0) ) then 
                     call multipole_select(1*asf_initial%asf(i)%totalJ,        &
                                                   asf_initial%asf(i)%parity,  &
                           1*asf_final%asf(j)%totalJ,asf_final%asf(j)%parity,  &
                           mult,gauge,nmult)
                     if (nmult == 0) cycle
                     nt = nt + 1
                     transition(nt)%asfi     = i;   transition(nt)%asff = j
                     transition(nt)%level_i  = asf_initial%asf(i)%level_No
                     transition(nt)%level_f  = asf_final%asf(j)%level_No
                     transition(nt)%energy   = asf_initial%asf(i)%energy -     &
                                               asf_final%asf(j)%energy
                     !
                     ! Read in experimental energies if appropriate
                     if (einstein_apply_exp_energies) then
                        if (energy_inverse) then
                        energy = energy_factor / transition(nt)%energy
                        else
                        energy = energy_factor * transition(nt)%energy
                        end if
                      2 print *, "Transition",transition(nt)%level_i,"-",      &
                           transition(nt)%level_f,                             &
                           "has ab-initio energy E_theo = ",                   &
                           energy,trim(energy_unit),";"
                        print *," enter E_exp (in ",trim(energy_unit),&
                           ") or 0.0 to use E_theo:"
                        read(*,*,err=2) energy_exp
                        if (energy_exp /= 0) then
                           if (energy_inverse) then
                              transition(nt)%energy = energy_factor/energy_exp
                           else
                              transition(nt)%energy = energy_exp/energy_factor
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
      if (einstein_sort_transition_energy) then
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
         call multipole_select(1*asf_initial%asf(i)%totalJ,          &
                                       asf_initial%asf(i)%parity,    &
               1*asf_final%asf(j)%totalJ,asf_final%asf(j)%parity,    &
               mult,gauge,nmult)
         if (rabs_use_stop   .and.   nmult == 0) then
            stop "einstein_set_transitions(): program stop A."
         end if
         !
         transition(nt)%number_of_mlines = nmult
         allocate( transition(nt)%mline(1:nmult) )
         do  m = 1,nmult
            transition(nt)%mline(m)%multipole  = mult(m)    
            transition(nt)%mline(m)%gauge      = gauge(m)    
            transition(nt)%mline(m)%amplitude  = zero 
         end do
      end do
      !
      print *, "  ",number_of_transitions,                          &
               " transitions have been initialized and will be "//  &
               "calculated in this run of the program."
      !
   end subroutine einstein_set_transitions
   !
   !
   subroutine einstein_transition_properties(tline)
   !--------------------------------------------------------------------
   ! Calculates all selected transition properties of line i from 
   ! the amplitudes of the individual multipole lines. 
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(einstein_transition), intent(inout) :: tline
      !
      integer       :: j, L, stream
      real(kind=dp) :: energy, sumb, sumc, einstein_A, einstein_B,          &
                       oscillator,decay_width
      !
      ! Print a short summary of the line if required
      if (einstein_print_each_line) then
         stream = 6
       3 write(stream,*) " "
         write(stream,4) " Results for the transition ",                    &
                         tline%level_i," - ",tline%level_f,": "//           &
                         trim(angular_momentum_string(tline%totalJ_i)) //   &
                         " " // tline%parity_i // "   ----> " //            &
                         trim(angular_momentum_string(tline%totalJ_f)) //   &
                         " " // tline%parity_f               
         write(stream,*) " --------------------------------------------"//  &
                         "--------------------"
         write(stream,*) " "
         if (energy_inverse) then
            energy = energy_factor / tline%energy
         else
            energy = energy_factor * tline%energy
         end if
         !
         sumc = zero;   sumb = zero
         do  j = 1,tline%number_of_mlines
            select case(tline%mline(j)%multipole(1:2))
            case("E1","M1")
               L = 1
            case("E2","M2")
               L = 2
            case("E3","M3")
               L = 3
            case("E4","M4")
               L = 4
            case default
               stop "einstein_transition_properties(): program stop A."
            end select
            !
            call multipole_convert_probability(tline%totalJ_i,tline%totalJ_f,L,&
                               tline%energy,tline%mline(j)%amplitude,          &
                               einstein_A,einstein_B,oscillator,decay_width)
            !
            if (tline%mline(j)%gauge == "Babushkin"  .or.    &
                tline%mline(j)%gauge == "Magnetic ")   then
               sumb = sumb + einstein_A
            end if
            if (tline%mline(j)%gauge == "Coulomb  "  .or.    &
                tline%mline(j)%gauge == "Magnetic ")   then
               sumc = sumc + einstein_A
            end if
         end do
         !
         !
         write(stream,*) "   Transition energy = ",                    &
                         energy,trim(energy_unit)
         write(stream,*) "   Total rate        = ",sumb,               &
                         trim(einstein_rate_unit),"  (Babushkin gauge) "
         write(stream,*) "                     = ",sumc,               &
                         trim(einstein_rate_unit),"  (Coulomb gauge) "
         !
         write(stream,*) " "
         write(stream,1) trim(einstein_rate_unit)
         write(stream,*) " ------------------------------------------------"
         do  j = 1,tline%number_of_mlines
            select case(tline%mline(j)%multipole(1:2))
            case("E1","M1")
               L = 1
            case("E2","M2")
               L = 2
            case("E3","M3")
               L = 3
            case("E4","M4")
               L = 4
            case default
               stop "einstein_transition_properties(): program stop B."
            end select
            !
            call multipole_convert_probability(tline%totalJ_i,tline%totalJ_f,L,&
                               tline%energy,tline%mline(j)%amplitude,          &
                               einstein_A,einstein_B,oscillator,decay_width)
            write(stream,2) tline%mline(j)%multipole,                          &
                            tline%mline(j)%gauge(1:9),                         &
                            tline%mline(j)%amplitude,einstein_A
         end do
         write(stream,*) " "
       1 format("   Mp   Gauge          Amplitude       A (",a,")")
       2 format(3x,a2,3x,a9,4x,1pe12.5,4x,1pe12.5)
       4 format(a,i4,a,i4, a)
         !
         ! Re-cycle once on stream 24
         if (stream == 6) then
            stream = 24;   goto 3
         end if
      end if
      !
   end subroutine einstein_transition_properties
   !  
   !
   subroutine einstein_write_trn_file(stream)
   !--------------------------------------------------------------------
   ! Writes out a  .trn transition energy and amplitude file for further 
   ! data processing with the EINSTEIN program on stream.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream 
      integer             :: i, j 
      character(len=3)    :: month
      character(len=8)    :: cdate
      character(len=10)   :: ctime
      !
      write(stream, "(i6,a)") number_of_transitions," = Number of transitions"
      write(stream,*)
      do  i = 1,number_of_transitions
         write(stream,1) transition(i)%level_i,  transition(i)%level_f,       &
                         transition(i)%totalJ_i, transition(i)%parity_i,      &
                         transition(i)%totalJ_f, transition(i)%parity_f,      &
                         transition(i)%number_of_mlines, transition(i)%energy,&
                         asf_initial%asf(transition(i)%asfi)%energy,          &
                         (transition(i)%mline(j)%multipole,                   &
                          transition(i)%mline(j)%gauge,                       &
                          transition(i)%mline(j)%amplitude,                   &
                          j=1,transition(i)%number_of_mlines)
      end do
    1 format(i4," -",i4,2x,2(i4,a1,1x),2x,i2,1pe14.7," a.u.   (",             &
             1pe14.7,")",3x,50(a2,1x,a9,1x,1pe14.7,3x))
      write(stream,*)
      !
      ! Get the date and time of day; append this information to the 
      ! einstein.trn transition energy and amplitude file to facilitate 
      ! identification of this file; this information is not to be read in
      !
      call date_and_time(date=cdate,time=ctime)
      month = get_month(cdate(5:6))
      write(stream,*) "EINSTEIN run at "//ctime(1:2)//":"//ctime(3:4)//":"//      &
                  ctime(5:6)//" on "//month//                                 &
                  " "//cdate(7:8)//" "//cdate(1:4)//"."
      !
   end subroutine einstein_write_trn_file
   !
end module rabs_einstein
