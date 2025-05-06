module rabs_toolbox_or
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains procedures which help carry out a number of
! 'utility' tasks within the RATIP environment. 
!-----------------------------------------------------------------------
   !
   use rabs_constant
   !!x use rabs_cesd
   use rabs_input_dialog
   use rabs_file_handling
   use rabs_print
   !!x use rabs_marco
   !!x use rabs_spectrum
   use rabs_toolbox_aux
   !!x use rabs_xl
   implicit none
   !
   public  :: toolbox_oscl_A_tables
                 ! Controls the compilation of transition probability
                 ! tabulations from REOS/EINSTEIN output.
   public  :: toolbox_oscl_branching_ratios
                 ! Controls the calculation of branching ratios from REOS/
                 ! EINSTEIN output.
   !!x private :: toolbox_oscl_gather_transitions
   !!x               ! Collects all transition information from one or several .trn
   !!x               ! REOS/EINSTEIN transition data files. 
   public  :: toolbox_oscl_transform_rate
                 ! Transforms some given transition amplitude or rate into
		 ! an alternative (but equivalent) form.
   public  :: toolbox_overlap
                 ! Controls the calculation of the overlap part.
   private :: toolbox_overlap_set_integrals
                 ! Initializes and calculates the array of overlap integrals
                 ! and prints these overlaps in a neat format.
   public  :: toolbox_radial_properties
                 ! Calculates and writes out the energies and radial expectation
                 ! values of the radial orbitals.
   public  :: toolbox_reduce_mix_file
                 ! Reduces a .mix file based on ASF level numbers or J 
		 ! and P values given by the user.
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
contains
   !
   subroutine toolbox_or()
   !--------------------------------------------------------------------
   ! Calls: 
   !--------------------------------------------------------------------
      !
      print *, "**************************"
      print *, "*** Not yet implmented ***"
      print *, "**************************"
      !
   end subroutine toolbox_or
   !
   !
   subroutine toolbox_oscl_A_tables()
   !--------------------------------------------------------------------
   ! Controls the generation of transition probability tables from REOS 
   ! output. It assumes the input from one or several REOS computations 
   ! which have been carried out independently and whose results have been 
   ! written to  .trn REOS transition data files. 
   !
   ! The procedure first reads in the information from all .trn files
   ! and, then, prints it in neat format. The main purpose is to combine 
   ! the results from independent runs of REOS and to provide a somewhat
   ! larger flexibility in presenting the output.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer            :: ce, i, j, k, no_transitions
      character(len=32)  :: record
      logical            :: yes, ascending_order, first
      real(kind=dp)      :: current, last, energy, energy_au, fa,fb,        &
                            wa, wb, wc, einstein_A, einstein_B, oscillator, &
                            decay_width, tb, tb_au, tb_cm, tb_ev, tb_sec,   &
                            tc, tc_au, tc_cm, tc_ev, tc_sec
      type(transition_from_reos), dimension(:), allocatable ::                &
                            lines_ordered, lines_unordered
      !
      ! Determine the units for the printout and further optional input
      call input_energy_unit()
      !
      !!x ! Determine the units for the printout and further optional input
      !!x 1 print *, "Which units are to be used to print the transition energies ?"
      !!x print *, "    A       : Angstrom;"
      !!x print *, "    eV      : electron volts;"
      !!x print *, "    Hartree : Hartree atomic units;"
      !!x print *, "    Hz      : Hertz;"
      !!x print *, "    Kayser  : [cm**(-1)];"
      !!x read (*, "(a)") record
      !!x if (len_trim(record) > 0) then
      !!x    record = adjustl(record)
      !!x    select case(record(1:7))
      !!x    case("A      ", "eV     ", "Hartree", "Hz     ", "Kayser ")
      !!x       toolbox_energy_unit = record(1:7)      
      !!x    case default
      !!x       print *, "Unable to decode the units '"//record(1:7)// &
      !!x                "' to be used for the transition energies; reenter ..."
      !!x       goto 1
      !!x    end select
      !!x else 
      !!x    toolbox_energy_unit = "eV     "
      !!x end if
      !
      ! Default speed of light
      c = c_vacuum
      !
      !!x ! Determine conversion factor for calculating energies in this unit
      !!x select case(toolbox_energy_unit)
      !!x case("A      ")
      !!x    toolbox_energy_factor  = 1.0e8_dp / convert_au_to_kaysers
      !!x    toolbox_energy_inverse = .true.
      !!x case("eV     ")
      !!x    toolbox_energy_factor  = convert_au_to_ev
      !!x case("Hartree")
      !!x    toolbox_energy_factor  = one
      !!x case("Hz     ")
      !!x    toolbox_energy_factor  = convert_au_to_kaysers * c_vacuum_in_cm_per_s
      !!x case("Kayser ")
      !!x    toolbox_energy_factor  = convert_au_to_kaysers
      !!x end select
      !
      print *, "Transition energies are printed in ascending order; " //&
               "use a decending order instead ?"
      yes = get_yes_stream()
      if (yes) then
         ascending_order = .false.
      else
         ascending_order = .true.
      end if
      !
      print *, "Einstein A and B coefficients are printed in SI units;"
      print *, " use Hartree atomic units instead ?"
      yes = get_yes_stream()
      if (yes) toolbox_print_AB_in_hartree = .true.
      !
      allocate( lines_unordered(1:2000), lines_ordered(1:2000)  )
      call toolbox_gather_transitions(lines_unordered,no_transitions,2000)
      !
      ! Now order all energies and set the 'ordered list'
      ce   = 0  ! current energy
      last = -1.0e99_dp
      !
    2 current =  1.0e99_dp
      !
      do  i = 1,no_transitions
         if (lines_unordered(i)%energy > last  .and.  &
             lines_unordered(i)%energy < current)  then
            current = lines_unordered(i)%energy
         end if
      end do
      !
      do  i = 1,no_transitions
         if (lines_unordered(i)%energy == current) then
            ce = ce + 1
            lines_ordered(ce)%level_i   = lines_unordered(i)%level_i
            lines_ordered(ce)%level_f   = lines_unordered(i)%level_f
            lines_ordered(ce)%totalJ_i  = lines_unordered(i)%totalJ_i
            lines_ordered(ce)%totalJ_f  = lines_unordered(i)%totalJ_f
            lines_ordered(ce)%parity_i  = lines_unordered(i)%parity_i
            lines_ordered(ce)%parity_f  = lines_unordered(i)%parity_f
            lines_ordered(ce)%energy    = current
            lines_ordered(ce)%file_name = lines_unordered(i)%file_name
            lines_ordered(ce)%number_of_mlines       =                      &
                                          lines_unordered(i)%number_of_mlines
            lines_ordered(ce)%initial_energy         =                      &
                                          lines_unordered(i)%initial_energy
            lines_ordered(ce)%total_A_rate_Coulomb   = zero 
            lines_ordered(ce)%total_A_rate_Babushkin = zero 
            !
            allocate( lines_ordered(ce)%mline(                             &
                                     1:lines_ordered(ce)%number_of_mlines) )
            do  j = 1,lines_ordered(ce)%number_of_mlines
               lines_ordered(ce)%mline%amplitude =                         &
                                          lines_unordered(i)%mline%amplitude
               lines_ordered(ce)%mline%multipole =                         &
                                          lines_unordered(i)%mline%multipole
               lines_ordered(ce)%mline%gauge     =                         &
                                          lines_unordered(i)%mline%gauge
               !
               select case(lines_ordered(ce)%mline(j)%multipole)
               case("E1", "M1");   lines_ordered(ce)%mline(j)%rank = 1
               case("E2", "M2");   lines_ordered(ce)%mline(j)%rank = 2
               case("E3", "M3");   lines_ordered(ce)%mline(j)%rank = 3
               case("E4", "M4");   lines_ordered(ce)%mline(j)%rank = 4
               case("E5", "M5");   lines_ordered(ce)%mline(j)%rank = 5
               case default
                  stop "toolbox_control_branching_ratios(): program stop A."
               end select
               !
               call toolbox_convert_probability(lines_ordered(ce)%totalJ_i,  &
                                             lines_ordered(ce)%totalJ_f,     &
                                             lines_ordered(ce)%mline(j)%rank,&
                                             lines_ordered(ce)%energy,       &
                               lines_ordered(ce)%mline(j)%amplitude,         &
                               lines_ordered(ce)%mline(j)%einstein_A,wa,wb,wc)
               lines_ordered(ce)%mline(j)%branching_ratio = zero
            end do
         end if
      end do
      !
      last = current
      if (ce < no_transitions) goto 2
      !
      ! Determine the total A rates of the upper levels; this total rate
      ! is assigned to the first line which belongs to this level while it
      ! is set to -one for all other lines from this level
      do  i = 1,no_transitions
         wa = lines_ordered(i)%initial_energy;  wb = zero;   wc = zero
         do  k = 1,no_transitions
            if (wa == lines_ordered(k)%initial_energy) then
               do  j = 1,lines_ordered(k)%number_of_mlines
                  select case(lines_ordered(k)%mline(j)%gauge)
                  case("Babushkin")
                     wb = wb + lines_ordered(k)%mline(j)%einstein_A
                  case("Coulomb  ")
                     wc = wc + lines_ordered(k)%mline(j)%einstein_A
                  case("Magnetic ")
                     wb = wb + lines_ordered(k)%mline(j)%einstein_A
                     wc = wc + lines_ordered(k)%mline(j)%einstein_A
                  case default
                     stop "toolbox_control_branching_ratios(): program stop B."
                  end select
               end do
            end if
         end do
         !
         first = .true.
         do  k = 1,no_transitions
            if (wa == lines_ordered(k)%initial_energy) then
               if (first) then
                  lines_ordered(k)%total_A_rate_Babushkin = wb
                  lines_ordered(k)%total_A_rate_Coulomb   = wc
                  first = .false.
               else
                  lines_ordered(k)%total_A_rate_Babushkin = -one
                  lines_ordered(k)%total_A_rate_Coulomb   = -one
               end if
            end if
         end do
      end do
      !
      ! Print the results in a neat table similar to the REOS summary file;
      ! Only the transition probabilities and branching ratios are printed.
      write(*,3)
    3 format(/ &
         /40x,"===========================================================",&
         /40x,"|  Summary of all Transition Probabilities and Lifetimes  |",&
         /40x,"===========================================================" )
      !
      if (energy_unit == "A           "  .and.  &
          toolbox_print_AB_in_hartree)       then
         write(*,4)
         write(*,6) 
      else if (energy_unit == "A      ") then
         write(*,5)
         write(*,7) 
      else if (toolbox_print_AB_in_hartree) then
         write(*,4)
         write(*,8) trim(energy_unit)
      else 
         write(*,5)
         write(*,9) trim(energy_unit)
      endif
    4 format(// 1x,140("-"),                                                 &
              / 2x,"LevI-LevF  I- J / Parity -F      Energy   ",             &
                   "Multipol   Gauge         Einstein coefficients",         &
                   "       Oscillator    Decay width    .trn File Name",     &
              /71x,"                     " )
    5 format(// 1x,140("-"),                                                 &
              / 2x,"LevI-LevF  I- J / Parity -F      Energy   ",             &
                   "Multipol   Gauge         Einstein coefficients",         &
                   "       Oscillator    Decay width    .trn File Name",     &
              /71x,"-1           3 -2 -1 " )   
    6 format(  30x,"   (Angstroms)",                                         &
               23x,"A (a.u.)    gB (a.u.)        strength GF       (eV) ")
    7 format(  30x,"   (Angstroms)",                                         &
               23x,"A (s  )     gB (m s  J  )    strength GF       (eV) ")
    8 format(  30x,"     (",a4,")   ",                                       &
               23x,"A (a.u.)    gB (a.u.)        strength GF       (eV) ")
    9 format(  30x,"     (",a4,")   ",                                       &
               23x,"A (s  )     gB (m s  J  )    strength GF       (eV) ")
      !
      write(*,10)
   10 format(  1x,140('-') )
      !
      do i = 1,no_transitions
         do j = 1,lines_ordered(i)%number_of_mlines
            energy_au = lines_ordered(i)%energy
            if (energy_inverse) then
               energy = energy_factor / lines_ordered(i)%energy
            else
               energy = energy_factor * lines_ordered(i)%energy
            end if
            !
            if (toolbox_print_AB_in_hartree) then
               fa = one;   fb = one
            else
               fa = convert_einstein_a_to_si
               fb = convert_einstein_b_to_si
            end if
            call toolbox_convert_probability(lines_ordered(i)%totalJ_i,     &
                                          lines_ordered(i)%totalJ_f,     &
                                          lines_ordered(i)%mline(j)%rank,&
                                          lines_ordered(i)%energy,       &
                            lines_ordered(i)%mline(j)%amplitude,         &
                            einstein_A,einstein_B,oscillator,decay_width)
            !
            write(*,11)     lines_ordered(i)%level_i,lines_ordered(i)%level_f, &
               trim(angular_momentum_string(lines_ordered(i)%totalJ_i,4)),     &
               lines_ordered(i)%parity_i,                                      &
               trim(angular_momentum_string(lines_ordered(i)%totalJ_f,4)),     &
               lines_ordered(i)%parity_f,energy,                               &
               lines_ordered(i)%mline(j)%multipole,                            &
               lines_ordered(i)%mline(j)%gauge,                                &
               einstein_A*fa, einstein_B*fb, oscillator, decay_width,          &
               lines_ordered(i)%file_name
         end do
      end do
      write(*,10)
      !
   11 format(2x,i3," -",i3,3x,a4,1x,a1,4x,a4,1x,a1,3x,1pd12.5,3x,a2,     &
             4x,a9,3x,1pd12.5,3x,1pd12.5,3x,1pd12.5,3x,1pd12.5,4x,a20)
      !
      !
      ! Print lifetimes and width of levels
      write(*,*)
      write(*,12)
   12 format(//"Radiative lifetimes and widths"                          &
              /"------------------------------"                          &
            ///" LeveL",6x,"Gauge",15x,"Lifetime",33x,"Width",28x,"File name" &
              /" -----",6x,"-----",15x,"--------",12x,                   &
               "---------------------------------------------",9x,       &
               "---------",&
          /32x,"Seconds",13x,"Hartrees",12x,"Kaysers",16x,"eV"/)
      !
      do  i = 1,no_transitions
         if (lines_ordered(i)%total_A_rate_Coulomb   == -one   .and.     &
             lines_ordered(i)%total_A_rate_Babushkin == -one) then
            cycle
         end if
         tb = lines_ordered(i)%total_A_rate_Babushkin
         tc = lines_ordered(i)%total_A_rate_Coulomb
         if (.not.toolbox_print_AB_in_hartree) then
            tb =lines_ordered(i)%total_A_rate_Babushkin*convert_einstein_a_to_si
            tc =lines_ordered(i)%total_A_rate_Coulomb  *convert_einstein_a_to_si
            tb_cm  = tb / c_vacuum_in_cm_per_s 
            tc_cm  = tc / c_vacuum_in_cm_per_s 
            tc_au  = tc_cm / convert_au_to_kaysers
            tb_au  = tb_cm / convert_au_to_kaysers
            tc_sec = one / tc
            tb_sec = one / tb
            tc_ev  = tc_au * convert_au_to_ev
            tb_ev  = tb_au * convert_au_to_ev
         else
            tc_au  = lines_ordered(i)%total_A_rate_Coulomb
            tb_au  = lines_ordered(i)%total_A_rate_Babushkin
            tc_cm  = tc * convert_au_to_kaysers
            tb_cm  = tb * convert_au_to_kaysers
            tc_sec = one / ( tc_cm * c_vacuum_in_cm_per_s )
            tb_sec = one / ( tb_cm * c_vacuum_in_cm_per_s )
            tc_ev  = tc_au * convert_au_to_ev
            tb_ev  = tb_au * convert_au_to_ev
         end if
         !
         if (abs((tb-tc)/tb) < eps10) then
            write(*,13) lines_ordered(i)%level_i,tc_sec,tc_au,tc_cm,tc_ev, &
                        lines_ordered(i)%file_name
         else
            write(*,14) lines_ordered(i)%level_i,tb_sec,tb_au,tb_cm,tb_ev, &
                        lines_ordered(i)%file_name
            write(*,15) tc_sec,tc_au,tc_cm,tc_ev,lines_ordered(i)%file_name
         end if
      13 format(1x,i4,6x,"Magnetic:  ",4(1pd20.7),4x,a20/) 
      14 format(1x,i4,6x,"Babushkin: ",4(1pd20.7),4x,a20)
      15 format(     11x,"Coulomb:   ",4(1pd20.7),4x,a20/)
         !
      end do
      !
   end subroutine toolbox_oscl_A_tables
   !
   !
   subroutine toolbox_oscl_branching_ratios()
   !--------------------------------------------------------------------
   ! Controls the calculations of branching ratios from REOS output.
   ! It assumes the input from one or several REOS computations which have 
   ! been carried out independently and whose results have been written 
   ! to  .trn REOS Transition data files. 
   !
   ! The procedure first reads in all information and samples the 
   ! total energies of the initial states in order to 'determine' 
   ! branching ratios.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer            :: ce, i, j, k, no_transitions
      character(len=32)  :: record
      logical            :: yes, ascending_order
      real(kind=dp)      :: current, last, energy, energy_au, factor, wa, wb, wc
      type(transition_from_reos), dimension(:), allocatable ::                &
                            lines_ordered, lines_unordered
      !
      ! Determine the units for the printout and further optional input
      call input_energy_unit()
      !
      !!x ! Determine the units for the printout and further optional input
      !!x 1 print *, "Which units are to be used to print the transition energies ?"
      !!x print *, "    A       : Angstrom;"
      !!x print *, "    eV      : electron volts;"
      !!x print *, "    Hartree : Hartree atomic units;"
      !!x print *, "    Hz      : Hertz;"
      !!x print *, "    Kayser  : [cm**(-1)];"
      !!x read (*, "(a)") record
      !!x if (len_trim(record) > 0) then
      !!x    record = adjustl(record)
      !!x    select case(record(1:7))
      !!x    case("A      ", "eV     ", "Hartree", "Hz     ", "Kayser ")
      !!x       toolbox_energy_unit = record(1:7)      
      !!x    case default
      !!x       print *, "Unable to decode the units '"//record(1:7)// &
      !!x                "' to be used for the transition energies; reenter ..."
      !!x       goto 1
      !!x    end select
      !!x else 
      !!x    toolbox_energy_unit = "eV     "
      !!x end if
      !!x !
      ! Default speed of light
      c = c_vacuum
      !
      !!x ! Determine conversion factor for calculating energies in this unit
      !!x select case(toolbox_energy_unit)
      !!x case("A      ")
      !!x    toolbox_energy_factor  = 1.0e8_dp / convert_au_to_kaysers
      !!x    toolbox_energy_inverse = .true.
      !!x case("eV     ")
      !!x    toolbox_energy_factor  = convert_au_to_ev
      !!x case("Hartree")
      !!x    toolbox_energy_factor  = one
      !!x case("Hz     ")
      !!x    toolbox_energy_factor  = convert_au_to_kaysers * c_vacuum_in_cm_per_s
      !!x case("Kayser ")
      !!x    toolbox_energy_factor  = convert_au_to_kaysers
      !!x end select
      !
      print *, "Transition energies are printed in ascending order; " //&
               "use a decending order instead ?"
      yes = get_yes_stream()
      if (yes) then
         ascending_order = .false.
      else
         ascending_order = .true.
      end if
      !
      allocate( lines_unordered(1:2000), lines_ordered(1:2000) )
      call toolbox_gather_transitions(lines_unordered,no_transitions,2000)
      !
      ! Now order all energies and set the 'ordered list'
      ce   = 0  ! current energy
      last = -1.0e99_dp
      !
    2 current =  1.0e99_dp
      !
      do  i = 1,no_transitions
         if (lines_unordered(i)%energy > last  .and.  &
             lines_unordered(i)%energy < current)  then
            current = lines_unordered(i)%energy
         end if
      end do
      !
      do  i = 1,no_transitions
         if (lines_unordered(i)%energy == current) then
            ce = ce + 1
            lines_ordered(ce)%level_i   = lines_unordered(i)%level_i
            lines_ordered(ce)%level_f   = lines_unordered(i)%level_f
            lines_ordered(ce)%totalJ_i  = lines_unordered(i)%totalJ_i
            lines_ordered(ce)%totalJ_f  = lines_unordered(i)%totalJ_f
            lines_ordered(ce)%parity_i  = lines_unordered(i)%parity_i
            lines_ordered(ce)%parity_f  = lines_unordered(i)%parity_f
            lines_ordered(ce)%energy    = current
            lines_ordered(ce)%file_name = lines_unordered(i)%file_name
            lines_ordered(ce)%number_of_mlines   =                         &
                                          lines_unordered(i)%number_of_mlines
            lines_ordered(ce)%initial_energy     =                         &
                                          lines_unordered(i)%initial_energy
            allocate( lines_ordered(ce)%mline(                             &
                                     1:lines_ordered(ce)%number_of_mlines) )
            do  j = 1,lines_ordered(ce)%number_of_mlines
               lines_ordered(ce)%mline%amplitude =                         &
                                          lines_unordered(i)%mline%amplitude
               lines_ordered(ce)%mline%multipole =                         &
                                          lines_unordered(i)%mline%multipole
               lines_ordered(ce)%mline%gauge     =                         &
                                          lines_unordered(i)%mline%gauge
               !
               select case(lines_ordered(ce)%mline(j)%multipole)
               case("E1", "M1");   lines_ordered(ce)%mline(j)%rank = 1
               case("E2", "M2");   lines_ordered(ce)%mline(j)%rank = 2
               case("E3", "M3");   lines_ordered(ce)%mline(j)%rank = 3
               case("E4", "M4");   lines_ordered(ce)%mline(j)%rank = 4
               case("E5", "M5");   lines_ordered(ce)%mline(j)%rank = 5
               case default
                  stop "toolbox_control_branching_ratios(): program stop A."
               end select
               !
               call toolbox_convert_probability(lines_ordered(ce)%totalJ_i,  &
                                             lines_ordered(ce)%totalJ_f,     &
                                             lines_ordered(ce)%mline(j)%rank,&
                                             lines_ordered(ce)%energy,       &
                               lines_ordered(ce)%mline(j)%amplitude,         &
                               lines_ordered(ce)%mline(j)%einstein_A,wa,wb,wc)
               !!x print *, "i,j,energy,amplitude, A = ",ce,j,      &
               !!x           lines_ordered(ce)%energy,              &
               !!x           lines_ordered(ce)%mline(j)%amplitude,  &
               !!x           lines_ordered(ce)%mline(j)%einstein_A
               lines_ordered(ce)%mline(j)%branching_ratio = zero
            end do
         end if
      end do
      !
      last = current
      if (ce < no_transitions) goto 2
      !
      ! Determine the branching ratios for the individual multipole components
      do  i = 1,no_transitions
         wa = lines_ordered(i)%initial_energy;  wb = zero;   wc = zero
         do  k = 1,no_transitions
            if (wa == lines_ordered(k)%initial_energy) then
               do  j = 1,lines_ordered(k)%number_of_mlines
                  select case(lines_ordered(k)%mline(j)%gauge)
                  case("Babushkin")
                     wb = wb + lines_ordered(k)%mline(j)%einstein_A
                  case("Coulomb  ")
                     wc = wc + lines_ordered(k)%mline(j)%einstein_A
                  case("Magnetic ")
                     wb = wb + lines_ordered(k)%mline(j)%einstein_A
                     wc = wc + lines_ordered(k)%mline(j)%einstein_A
                  case default
                     stop "toolbox_control_branching_ratios(): program stop B."
                  end select
               end do
            end if
         end do
         !
         do  j = 1,lines_ordered(i)%number_of_mlines
            !!x print *, "wb,wc,i,j,einstein_A = ",wb,wc,i,j, &
            !!x          lines_ordered(i)%mline(j)%einstein_A
            select case(lines_ordered(i)%mline(j)%gauge)
            case("Babushkin", "Magnetic ")
               lines_ordered(i)%mline(j)%branching_ratio =                    &
                                      lines_ordered(i)%mline(j)%einstein_A / wb
            case("Coulomb  ")
               lines_ordered(i)%mline(j)%branching_ratio =                    &
                                      lines_ordered(i)%mline(j)%einstein_A / wc
            end select
         end do
      end do
      !
      ! Print the results in a neat table similar to the REOS summary file;
      ! Only the transition probabilities and branching ratios are printed.
      write(*,3)
    3 format(/ &
      /23x,"==================================================================",&
      /23x,"|  Summary of all Transition Probabilities and Branching Ratios  |",&
      /23x,"==================================================================" )
      !
      if (energy_unit == "Hartree") then
         write(*,4)
         write(*,5) 
      else 
         write(*,6)
         write(*,7) trim(energy_unit)
      end if
      ! 
    4 format(// 1x,112("-"),                                                 &
              / 2x,"LevI-LevF  I- J / Parity -F      Energy   ",             &
                   "Multipol   Gauge       Einstein-A",                      &
                   "     Branching      .trn File Name" )
    5 format(  30x,"     (Hartree)", 27x,"(a.u.)         Ratio")
    6 format(// 1x,112("-"),                                                 &
              / 2x,"LevI-LevF  I- J / Parity -F      Energy   ",             &
                   "Multipol   Gauge       Einstein-A",                      &
                   "     Branching      .trn File Name" )
    7 format(  30x,"     (",a4,") ", 27x,"(1/s)          Ratio")
    8 format(  1x,112('-') )
    9 format(2x,i3," -",i3,3x,a4,1x,a1,4x,a4,1x,a1,3x,1pd12.5,3x,a2,         &
             4x,a9,4x,1pd12.5,3x,1pd12.5,4x,a20)
      !
      write(*,8)
      do i = 1,no_transitions
         do j = 1,lines_ordered(i)%number_of_mlines
            energy_au = lines_ordered(i)%energy
            if (energy_inverse) then
               energy = energy_factor / lines_ordered(i)%energy
            else
               energy = energy_factor * lines_ordered(i)%energy
            end if
            if (energy_unit == "Hartree") then
               factor = one
            else
               factor = convert_einstein_a_to_si
            end if
            write(*,9)      lines_ordered(i)%level_i,lines_ordered(i)%level_f, &
               trim(angular_momentum_string(lines_ordered(i)%totalJ_i,4)),     &
               lines_ordered(i)%parity_i,                                      &
               trim(angular_momentum_string(lines_ordered(i)%totalJ_f,4)),     &
               lines_ordered(i)%parity_f,energy,                               &
               lines_ordered(i)%mline(j)%multipole,                            &
               lines_ordered(i)%mline(j)%gauge,                                &
               lines_ordered(i)%mline(j)%einstein_A*factor,                    &
               lines_ordered(i)%mline(j)%branching_ratio,                      &
               lines_ordered(i)%file_name
         end do
      end do
      write(*,8)
      !
   end subroutine toolbox_oscl_branching_ratios
   !
!!x    !
!!x    subroutine toolbox_oscl_gather_transitions(lines,no_transitions,lines_max)
!!x    !--------------------------------------------------------------------
!!x    ! Collects all transition information from one or several .trn 
!!x    ! REOS transition data files.
!!x    !
!!x    ! Calls: file_open(). 
!!x    !--------------------------------------------------------------------
!!x       !
!!x       integer, intent(in)  :: lines_max
!!x       integer, intent(out) :: no_transitions
!!x       type(transition_from_reos), dimension(:), intent(out) :: lines
!!x       !
!!x       integer            :: i, j, k, m, n, no_files, ierr, ios 
!!x 	  character(len=15)  :: file_id
!!x 	  character(len=512) :: record
!!x 	  character(len=256), dimension(10)	   :: toolbox_trn_file
!!x 	  type(multipole_from_reos), dimension(50) :: mlines
!!x 	  !
!!x 	  ! Request a number of .trn files and read the data
!!x 	1 print *, "Enter one (or several)  .trn REOS transition data file(s):"
!!x 	  read(*,"(a)")  record
!!x 	  !
!!x 	  ! Determine file names
!!x 	  no_files = 0
!!x 	2 record = adjustl(record)
!!x 	  !!x print *, "record = ",record
!!x 	  i = scan(record," ")
!!x 	  !!x print *, "i = ",i
!!x 	  if (i /= 1) then
!!x 	     no_files = no_files + 1
!!x 	     toolbox_trn_file(no_files) = record(1:i-1)
!!x 	     record(:) = record(i:)
!!x 	     goto 2
!!x 	  end if
!!x 	  !
!!x 	  ! Try to open these files and to gather all necessary information
!!x 	  no_transitions = 0;	m = 0
!!x 	  do  i = 1,no_files
!!x 	     call file_open(25,toolbox_trn_file(i),"formatted  ","old",ierr)
!!x 	     if (ierr /= 0  .or.  len_trim(toolbox_trn_file(i)) == 0) goto 1
!!x 	     !
!!x 	     ! Check the header of the file; if not as expected for unformatted
!!x 	     ! files, check formatted form
!!x 	     read (25,"(a15)",iostat=ios)  file_id
!!x 	     if (ios /= 0   .or.   file_id /= "REOS Transition") then
!!x 		print *, "ios, file header = ",ios, file_id
!!x 		print *, "Not a REOS Transition Data File;"
!!x 		close (25)
!!x 		goto 1
!!x 	     end if
!!x 	     !
!!x 	     read(25,*)
!!x 	     read(25, "(i6,a)") n
!!x 	     no_transitions = no_transitions + n
!!x 	     read(25,*)
!!x 	     do   k = 1,n
!!x 		m = m + 1
!!x 		if (m > lines_max) then
!!x 		   stop "toolbox_gather_transitions_reos(): program stop A."
!!x 		end if 
!!x 		read(25,4) lines(m)%level_i,  lines(m)%level_f,        &
!!x 			   lines(m)%totalJ_i, lines(m)%parity_i,       &
!!x 			   lines(m)%totalJ_f, lines(m)%parity_f,       &
!!x 			   lines(m)%number_of_mlines, lines(m)%energy, &
!!x 			   lines(m)%initial_energy,		       &
!!x 			   (mlines(j)%multipole, mlines(j)%gauge,      &
!!x 			    mlines(j)%amplitude, j=1,lines(m)%number_of_mlines)
!!x 		allocate( lines(m)%mline(1:lines(m)%number_of_mlines) )
!!x 		do  j = 1,lines(m)%number_of_mlines
!!x 		   lines(m)%mline(j)%multipole = mlines(j)%multipole
!!x 		   lines(m)%mline(j)%gauge     = mlines(j)%gauge
!!x 		   lines(m)%mline(j)%amplitude = mlines(j)%amplitude
!!x 		end do
!!x 		lines(m)%file_name = toolbox_trn_file(i)
!!x 	     end do
!!x 	     !
!!x 	     close (25)
!!x 	  end do
!!x 	4 format(i4,2x,i4,2x,2(i4,a1,1x),2x,i2,1pe14.7,9x,	     &
!!x 		 1pe14.7,4x,50(a2,1x,a9,1x,1pe14.7,3x))
!!x 	  !
!!x 	  if (m /= no_transitions) then
!!x 	     stop "toolbox_gather_transitions_reos(): program stop B."
!!x       end if 
!!x       !
!!x    end subroutine toolbox_oscl_gather_transitions
!!x    !
   !
   subroutine toolbox_oscl_transform_rate()
   !--------------------------------------------------------------------
   ! Transforms some given transition amplitude or rate into an alternative 
   ! (but equivalent) form.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer           :: sel, Ji, Jf, rank
      real(kind=dp)     :: fa, fb, energy, energy_au, amplitude, rate, rate_au,&
                           wa, einstein_A, einstein_B, oscillator, decay_width,&
			   osc, width
      character(len=6)  :: record
      !
      ! Determine the units for the printout and further optional input
      call input_energy_unit()
      !
      ! Default speed of light
      c = c_vacuum
      !
      print *, "Enter the transition energy in "//trim(energy_unit)//":"
      read  *, energy
      !
      energy = abs(energy)
      if (energy_inverse) then
         energy_au = energy_factor / energy
      else
         energy_au = energy / energy_factor
      end if
      print *, "energy (in "//trim(energy_unit)//"), energy_au = ", &
                energy, energy_au
      !
      print *, "Enter the initial-state angular momentum, J_i:"
      read(*,"(a)") record
      Ji = get_dinteger_from_string( trim(record) )
      !
      print *, "Enter the final-state angular momentum, J_f:"
      read(*,"(a)") record
      Jf = get_dinteger_from_string( trim(record) )
      !
      print *, "Enter the (+integer) rank of the multipole transition:"
      read *, rank
      rank = abs(rank)
      !
      print *, "Ji, Jf, rank = ",Ji, Jf
      !
    5 print *, " "
      print *, "Select a particular task:"
      print *, "   1 - Transform a transition amplitude in a.u. "
      print *, "   2 - Transform a transition rate A in a.u. "
      print *, "   3 - Transform a transition rate A in 1/s. "
      print *, "   4 - Transform an oscillator strength. "
      print *, "   5 - Transform a level width. "
      read  *, sel
      !
      select case(sel)
      case(1)
         print *, "Enter the transition amplitude in a.u.:"
	 read  *, amplitude
	 amplitude = abs(amplitude)
      case(2)
         print *, "Enter the transition rate in a.u.:"
	 read  *, rate_au
	 rate_au = abs(rate_au)
      case(3)
         print *, "Enter the Einstein A value in 1/s:"
	 read  *, rate
	 rate    = abs(rate)
         rate_au = rate / convert_einstein_a_to_si
      case(4)
         print *, "Enter oscillator strength:"
	 read  *, osc
	 osc     = abs(osc)
      case(5)
         print *, "Enter level widths in "//trim(energy_unit)//":"
	 read  *, width
	 width   = abs(width)
	 !
         if (energy_inverse) then
            width = energy_factor / width
         else
            width = width / energy_factor
         end if
	 width = width * convert_au_to_ev
         print *, " width_eV = ", width
	 !
      case default
         print *, "Selection cannot be recognized ... redo "
	 goto 5
      end select
      !
      select case(sel)
      case(2,3)
         call toolbox_convert_probability(Ji, Jf, rank, energy_au, one,  &
                                     wa,einstein_B,oscillator,decay_width)
         amplitude = sqrt( rate_au/wa )
      case(4)
         call toolbox_convert_probability(Ji, Jf, rank, energy_au, one,  &
                                     einstein_A,einstein_B,wa,decay_width)
         amplitude = sqrt( osc/wa )
      case(5)
         call toolbox_convert_probability(Ji, Jf, rank, energy_au, one,  &
                                      einstein_A,einstein_B,oscillator,wa)
         amplitude = sqrt( width/wa )
      end select
      !
      call toolbox_convert_probability(Ji, Jf, rank, energy_au, amplitude,  &
                                einstein_A,einstein_B,oscillator,decay_width)
      !
      print *, " "
      print *, "Alternative (but equivalent) information about the "// &
               "transition probability:"
      print *, "---------------------------------------------------"// &
               "-----------------------"
      write(*,7) "  Einstein A (a.u)    = ",einstein_A
      write(*,7) "  Einstein A (1/s)    = ",einstein_A*convert_einstein_a_to_si
      write(*,7) "  Einstein B (a.u)    = ",einstein_B
      write(*,7) "  Einstein B (SI)     = ",einstein_B*convert_einstein_b_to_si
      write(*,7) "  Oscillator strength = ",oscillator
      write(*,7) "  Decay withs (eV)    = ",decay_width
    7 format(a,1pe12.5)
      !
   end subroutine toolbox_oscl_transform_rate
   !
   !
   subroutine toolbox_overlap()
   !--------------------------------------------------------------------
   ! Controls the computation of overlap intagrals using two not quite 
   ! orthogonal sets of one-electron orbitals.
   !
   ! Calls: radgrd_grasp2k(), setqic_grasp2k(), file_get_csl_list(), 
   !        file_get_rwf(), toolbox_overlap_set_integrals().
   !--------------------------------------------------------------------
      !
      ! logical :: yes
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
      ! Get the occupied shells of the initial and final states
      call file_get_csl_list("Enter the name of the initial-state GRASP92 "//&
         "configuration symmetry list file:",asf_initial%csf_set)
      call file_get_csl_list("Enter the name of the final-state GRASP92 "  //&
         "configuration symmetry list file:",asf_final%csf_set)
      !
      ! Initialize and load the radial wavefunctions for the inital and
      ! final atomic states
      !!x call toolbox_initialize_rwf_storage(.false.,.true.,.true.)
      call initialize_rwf_storage(asf_initial, wave_initial)   
      call initialize_rwf_storage(asf_final,   wave_final)   
      call file_get_rwf(                                                  &
         "Enter the name of the initial-state Radial WaveFunction File:", &
	 asf_initial,wave_initial,.true.)
      call file_get_rwf(                                                  &
         "Enter the name of the final-state Radial WaveFunction File:",   &
	 asf_final,wave_final,.true.)
      !
      ! Calculate the overlaps and print results
      call toolbox_overlap_set_integrals()
      !
   end subroutine toolbox_overlap
   !
   !
   subroutine toolbox_overlap_set_integrals()
   !--------------------------------------------------------------------
   ! Calculates the non-orthogonal and overlap integrals between all
   ! bound orbitals of the corresponding final- and initial-state arrays.
   !
   ! Calls: rk_integral_grasp2k_ab().
   !--------------------------------------------------------------------
      !
      integer :: i, i1, i2, kappa, kappa_min, kappa_max, pqn_i, pqn_i_max, &
                 pqn_f, pqn_f_max, no_int
      integer, dimension(100)       :: pqn_ff, pqn_ii
      real(kind=dp), dimension(100) :: value
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
      allocate( overlap(kappa_min:kappa_max,pqn_f_max,pqn_i_max) )
      !
      ! Calculate the non-orthogonal and overlap integrals
      do  i1 = 1,wave_final%number_of_rwf
         do  i2 = 1,wave_initial%number_of_rwf
            if (wave_final%rwf(i1)%orbital%kappa == &
                wave_initial%rwf(i2)%orbital%kappa) then
               kappa = wave_final%rwf(i1)%orbital%kappa
               pqn_f = wave_final%rwf(i1)%orbital%n
               pqn_i = wave_initial%rwf(i2)%orbital%n
               overlap(kappa,pqn_f,pqn_i) = &
                  rk_integral_grasp2k_ab(wave_final,wave_initial,0,i1,i2)
            end if
         end do
      end do
      !
      write(*,*) " "
      write(*,*) "Overlap integrals <n_f kappa | n_i kappa> :"
      write(*,*) "-------------------------------------------"
      !
      ! Print the overlap integrals
      do  kappa = 1,30
         no_int = 0
         do  i1 = 1,wave_final%number_of_rwf
            do  i2 = 1,wave_initial%number_of_rwf
               if (wave_final%rwf(i1)%orbital%kappa   == -kappa   .and. &
                   wave_initial%rwf(i2)%orbital%kappa == -kappa) then
                  no_int = no_int + 1
                  pqn_ff(no_int) = wave_final%rwf(i1)%orbital%n
                  pqn_ii(no_int) = wave_initial%rwf(i2)%orbital%n
                  value(no_int)  = overlap(-kappa,pqn_ff(no_int),pqn_ii(no_int))
               end if
            end do
         end do
         !
         if (no_int > 0) then
            write(*,*) " "
            write(*,1) orbital_symmetry(-kappa),                         &
                      (orbital_name(pqn_ff(i),-kappa),                   &
                       orbital_name(pqn_ii(i),-kappa),value(i)," ",      &
                       i=1,no_int-1),                                    &
                       orbital_name(pqn_ff(no_int),-kappa),              &
                       orbital_name(pqn_ii(no_int),-kappa),value(i)
            if (no_int < 3) write(*,*) "  ------------ "
         end if
         !
         no_int = 0
         do  i1 = 1,wave_final%number_of_rwf
            do  i2 = 1,wave_initial%number_of_rwf
               if (wave_final%rwf(i1)%orbital%kappa   == kappa   .and. &
                   wave_initial%rwf(i2)%orbital%kappa == kappa) then
                  no_int = no_int + 1
                  pqn_ff(no_int) = wave_final%rwf(i1)%orbital%n
                  pqn_ii(no_int) = wave_initial%rwf(i2)%orbital%n
                  value(no_int)  = overlap(kappa,pqn_ff(no_int),pqn_ii(no_int))
               end if
            end do
         end do
         !
         if (no_int > 0) then
            write(*,*) " "
            write(*,1) orbital_symmetry(kappa),                             &
                      (orbital_name(pqn_ff(i),kappa),                       &
                       orbital_name(pqn_ii(i),kappa),value(i)," ",          &
                       i=1,no_int-1),                                       &
                       orbital_name(pqn_ff(no_int),kappa),                  &
                       orbital_name(pqn_ii(no_int),kappa),value(i)
            if (no_int < 3) write(*,*) "  ------------ "
         end if
      end do
      !
    1 format(" + ",a," symmetry: ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /  "   ------------ ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /  "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /  "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /  "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /  "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /  "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /  "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /  "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a))
      !
   end subroutine toolbox_overlap_set_integrals
   !
   !
   subroutine toolbox_radial_properties()
   !--------------------------------------------------------------------
   ! Calculates and writes out the one-electron energies and radial 
   ! expectation values of a set of orbitals which are obtained from a
   ! .mix Mixing Coefficient File. The expectation values of the operators 
   ! r^k are printed for the powers k = -3, -1,  1, and 2 which are 
   ! related to various physical operators describing the structure and 
   ! properties of atoms.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer            :: i, ierr, ios
      logical            :: yes
      character(len=31)  :: record
      character(len=256) :: orb_rwf_file
      real(kind=dp)      :: rm3, rm1, rp1, rp2
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
      print *, "The default radial grid parameters for this case are:"
      print *, " rnt = ",rnt_grasp2k,";"
      print *, " h   = ",h_grasp2k,  ";"
      print *, " hp  = ",hp_grasp2k, ";"
      print *, " n   = ",n_grasp2k,  ";"
      print *, " revise these values ?"
      yes = get_yes_stream()
      if (yes) then
         print *, "Enter rnt:"
         read *, rnt_grasp2k
         print *, "Enter h:"
         read *, h_grasp2k
         print *, "enter hp:"
         read *, hp_grasp2k
         print *, "enter n:"
         read *, n_grasp2k
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
      ! Get the configuration symmetry list
      call file_get_csl_list(                                       &
         "Enter the name of the configuration symmetry list file:", &
         asf_set%csf_set)
      !
      ! Initialize and load the radial wavefunctions; allow for both,
      ! formatted and unformatted orbital files 
      !!x call toolbox_initialize_rwf_storage(.true.,.false.,.false.)
      call initialize_rwf_storage(asf_set, wave)
      !
    1 print *, "Enter the name of the Radial WaveFunction File:"
      read (*,"(a)") orb_rwf_file
      if (len(trim(orb_rwf_file)) == 0) goto 1
      !
      call file_open(21,orb_rwf_file,"unformatted","old",ierr)
      if (ierr == 1) goto 1
      !
      read (21,iostat = ios) record(1:6)
      if (record(1:6) /= 'G92RWF') then
         close (21)
         !
         call file_open(21,orb_rwf_file,"formatted  ","old",ierr)
         if (ierr == 1) goto 1
         read (21,"(a31)",iostat = ios) record
         if (ios /= 0   .or.   &
             record(1:31) /= "G92RWF (formatted file version)") then
            print *, "ios, record(1:31) = ",ios, record(1:31)
            print *, "Not a G92RWF Radial WaveFunction File;"
            close (21)
            goto 1
         end if
         toolbox_use_formatted_rwf_file = .true.
      else
         toolbox_use_formatted_rwf_file = .false.
      endif
      !
      ! Load data from the  .rwf  file
      call load_rwf_file_grasp2k(wave,toolbox_use_formatted_rwf_file,ierr)
      if (ierr /= 0) then
         goto 1
      end if
      !
      ! Close the  .rwf  file
      close (21)
      !
      print *, "                                        -3     " //&
               "               -1                                          2  "
      print *, "Subshell    Energy (a.u.)             <r  >    " //&
               "             <r  >                  <r>                  <r > " 
      print *, " "
      do  i = 1,wave%number_of_rwf
         rm3 = rk_integral_grasp2k(wave,-3,i,i)
         rm1 = rk_integral_grasp2k(wave,-1,i,i)
         rp1 = rk_integral_grasp2k(wave, 1,i,i)
         rp2 = rk_integral_grasp2k(wave, 2,i,i)
         print 2, orbital_name(wave%rwf(i)%orbital%n,                         &
                  wave%rwf(i)%orbital%kappa),wave%rwf(i)%energy,rm3,rm1,rp1,rp2
      end do
    2 format(2x,a4,3x,5(1p,e19.12,3x))
      !
   end subroutine toolbox_radial_properties
   !
  !
  subroutine toolbox_reduce_mix_file(level_numbers)
   !--------------------------------------------------------------------
   ! Reduces a .mix file based either on the ASF level numbers or the
   ! J^P values as given by the user. All ASF from the old .mix file, 
   ! which match with the given specification, are printed in the new 
   ! .mix file. 
   !
   !
   ! Calls: get_yes_stream(), file_open(), load_mix_file_grasp2k, 
   !        util2_interprete_levels(), util2_write_mix_file
   !--------------------------------------------------------------------
      !
      logical, intent(in):: level_numbers
      !
      integer :: ierr, ios, i, j, number_of_levels, nojp, parity_position
      integer, dimension(1:5000)    :: levels
      integer, dimension(1:20)      :: j_chosen    
      logical                       :: yes, fail, is_selected, level_found 
      !
      character(len=20)  :: string
      character(len=120) :: record
      character(len=1), dimension(1:20)        :: parity_chosen
      !
      nojp                  = 0
      number_of_levels      = 0
      levels(:)             = 0
      !
      ! Open, check, load data from, and close, the  .csl  file
      call file_get_csl_list(                                         &
         "Enter the name of the configuration symmetry list file:",   &
         asf_set%csf_set)
      !
      ! Read in the old .mix file 
      call file_get_mix(                                                      &
         "Enter the name of the .mix mixing coefficients file to be reduced:",&
	 asf_set)
      !
      if (level_numbers) then
         !
         ! Read in the ASF serial numbers, which are to be printed in the 
	 ! new .mix file       
       4 print *, "Enter the serial numbers of the ASFs to be selected,"
         print *, "e.g. 1 3 4  7 - 20  48  69 - 85;"
         read (*, "(a)") record
         call toolbox_interprete_levels(record,levels,number_of_levels,fail)
         if (fail) then
            print *, "Unable to interprete the serial level numbers; redo ..."
            goto 4
         else if (number_of_levels == 0) then
            print *, "Unable to interprete the serial level numbers; redo ..."
            goto 4
         end if
         !       
         ! Check that the all givel ASF serial numbers are defined in the 
	 ! old .mix file
         do  i = 1,number_of_levels
            level_found = .false.
            do  j = 1,asf_set%noasf
               if (levels(i) == asf_set%asf(j)%level_No) then
                  level_found = .true.
                  exit
               endif
            end do
	    !              
            if (.not.level_found) then
               print *, "One or more of the chosen levels are not present "//&
	                "in the given .mix file"
               print *, "redo ..." 
               goto 4
            end if
         end do      
         goto 6
      else
         !
         ! Read in the JP information ...
       5 print *, "Enter the J and P values of the ASFs to be selected," 
         print *, "e.g. 3/2 -; 2 +; <cr> if done:"     
         read (*, "(a)") record
         if (len_trim(record) > 0) then
            !
            ! Check the input data
            parity_position = scan(record,"-")    
            if (parity_position == 0) then
               parity_position = scan(record,"+")   
               if (parity_position == 0) then
                   print *, "Unable to decode the states; reenter ..."
                   goto 5
               end if
            end if 
            !
            ! Number of defined JP pairs
            nojp                = nojp + 1 
            parity_chosen(nojp) = record(parity_position:parity_position + 1)
            string         = adjustl(record(1:parity_position-1))
            j_chosen(nojp) = get_dinteger_from_string(string)      
            goto 5
         end if       
         !
         if (nojp == 0) then
            print *, "Enter at least one JP pair; redo ..."
            goto 5
         else
            do i = 1,asf_set%noasf
               do j = 1,nojp
                  if (j_chosen(j)     == asf_set%asf(i)%totalJ .and. &
                      parity_chosen(j)== asf_set%asf(i)%parity)   then   
                     number_of_levels = number_of_levels + 1
                     levels(number_of_levels) = asf_set%asf(i)%level_No
                  end if
               end do
            end do
         end if 
         !
         if (number_of_levels == 0) then
            print *, "No ASFs in the .mix file with given J and P values;"
            print *, "reenter ..."
            goto 5
         end if
      end if
      !
      ! Write out the chosen ASFs in the new .mix file
      !
    6 call file_write_mix(asf_set,levels,number_of_levels)
      !
      do i = 1, asf_set%noasf
         deallocate( asf_set%asf(i)%eigenvector )
      end do
      deallocate( asf_set%asf ) 
      !      
  end subroutine toolbox_reduce_mix_file
  !
end module rabs_toolbox_or
