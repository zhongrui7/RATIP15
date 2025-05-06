module rabs_input_dialog
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module includes `frequent input dialogs' which are needed to
! read in the transition multipoles, units and other properties.
!-----------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_functions_string   
   use rabs_nucleus
   implicit none
   !
   public  :: input_energy_unit
                 ! Reads in and determines the energy unit for input/output. 
   public  :: input_grid_parameters
                 ! Reads in and determines the parameters for the radial grid.
   public  :: input_levels
                 ! Reads in an selected set of levels numbers.
   public  :: input_transition_multipoles
                 ! Reads in the transition multipoles from a given record. 
   public  :: input_transition_pairs
                 ! Reads in and selects pairs of transitions.
   public  :: input_transition_triples
                 ! Reads in and selects triples of transitions.
   !
   integer, dimension(200), public :: select_level_i, select_level_f, &
                                      select_level_m
   !
   !
contains
   !
   !
   subroutine input_energy_unit()
   !--------------------------------------------------------------------
   ! Reads in and determines the energy unit for input/output.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      character(len=256) :: record
      !
      ! Determine the units for the printout and further optional input
    1 print *, "Which units are to be used to enter and to print the"// &
               " energies of the continuum orbitals ?"
      print *, "    A       : Angstrom;"
      print *, "    eV      : electron volts;"
      print *, "    Hartree : Hartree atomic units;"
      print *, "    Hz      : Hertz;"
      print *, "    Kayser  : [cm**(-1)];"
      read (*, "(a)") record
      if (len_trim(record) > 0) then
         record = adjustl(record)
         select case(record(1:7))
         case("A      ", "eV     ", "Hartree", "Hz     ", "Kayser ")
            energy_unit = record(1:7)      
         case default
            print *, "Unable to decode the units '"//record(1:7)// &
                     "' to be used for the transition energies; reenter ..."
            goto 1
         end select
      else 
         energy_unit = "eV     "
      end if
      energy_inverse = .false.
      !
      call save_input(record,.false.)
      !
      ! Determine conversion factor for calculating energies in this unit
      select case(energy_unit)
      case("A      ")
         energy_factor  = 1.0e8_dp / convert_au_to_kaysers
         energy_inverse = .true.
      case("eV     ")
         energy_factor  = convert_au_to_ev
      case("Hartree")
         energy_factor  = one
      case("Hz     ")
         energy_factor  = convert_au_to_per_sec
      case("Kayser ")
         energy_factor  = convert_au_to_kaysers
      end select
      !
   end subroutine input_energy_unit
   !
   !
   subroutine input_grid_parameters(keystring)
   !--------------------------------------------------------------------
   ! Determines of reads in the new parameters for specifying the radial
   ! grid.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      character(len=*), intent(in) :: keystring
      logical                      :: yes
      character(len=256) :: record
      !
      select case(keystring)
      case("standard")
         !
         if (nuclear_model == "point") then
	    rnt_grasp2k = exp(-65.0_dp/16.0_dp) / nuclear_charge
	    h_grasp2k   = half**4
	    n_grasp2k   = 220
         else
	    rnt_grasp2k = 2.0e-6_dp
	    h_grasp2k   = 5.0e-2_dp
	    n_grasp2k   = 390
         endif
	 !
      case("modify")
	 !
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
	    !
	    write(stream_input,"(1pe12.5)") rnt_grasp2k
	    backspace (stream_input)
	    read (stream_input,*) record
            call save_input(trim(record)//" :: rnt_grasp2k",.false.)
	    !
            print *, "Enter h:"
            read *, h_grasp2k
	    !
	    write(stream_input,"(1pe12.5)") h_grasp2k
	    backspace (stream_input)
	    read (stream_input,*) record
            call save_input(trim(record)//" :: h_grasp2k",.false.)
	    !
            print *, "enter hp:"
            read *, hp_grasp2k
	    !
	    write(stream_input,"(1pe12.5)") hp_grasp2k
	    backspace (stream_input)
	    read (stream_input,*) record
            call save_input(trim(record)//" :: hp_grasp2k",.false.)
	    !
            print *, "enter n:"
            read *, n_grasp2k
	    !
	    write(stream_input,"(i6)") n_grasp2k
	    backspace (stream_input)
	    read (stream_input,*) record
            call save_input(trim(record)//"        :: n_grasp2k",.false.)
	    !
         end if
	 !
     case default
         stop "input_grid_parameters(): program stop A."
      end select
      !
   end subroutine input_grid_parameters
   !
   !
   subroutine input_levels(number_of_levels)
   !--------------------------------------------------------------------
   ! Reads in some selected level numbers.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(inout) :: number_of_levels
      !
      integer            :: level_i, blank_position
      logical            :: yes
      character(len=256) :: record
      character(len=20 ) :: string
      !
      number_of_levels  = 0
      !
      print *, "Select individual levels ?"
      yes = get_yes_stream()
      if (yes) then
    4    print *, "Enter a list of levels numbers:"
         read (*, "(a)") record
         call save_input(record,.false.)
	 record         = adjustl(record)
         !
    5    blank_position = scan(record," ")
	 string         = adjustl(record(1:blank_position-1))
	 number_of_levels  = number_of_levels + 1
	 select_level_i(number_of_levels) = get_integer_from_string(string)
	 record         = adjustl(record(blank_position+1:256))
         if (len_trim(record) > 0)  goto 5
      end if
      !
      if (number_of_levels == 0) then
         print *, "No individual level number has been selected; the "// &
        	  " program continues the default path."
      end if
      !
   end subroutine input_levels
   !
   !
   subroutine input_transition_multipoles(number_of_multipoles,multipole)
   !--------------------------------------------------------------------
   ! Reads in the transition multipoles from a given record.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(inout)                        :: number_of_multipoles
      character(len=2), dimension(:), intent(inout) :: multipole
      !
      character(len=256) :: record
      !
    1 print *, "Enter the transition multipoles, e.g.  E1 M2 ... :"
      read (*, "(a)") record
      call save_input(record,.false.)
      !
    2 if (len_trim(record) > 0) then
         record = adjustl(record)
         select case(record(1:2))
         case("E1", "E2", "E3", "E4", "E5", &
              "M1", "M2", "M3", "M4", "M5")
            number_of_multipoles = number_of_multipoles + 1
            if (number_of_multipoles > 20) then
               stop "input_transition_multipoles(): program stop A."
            end if
            multipole(number_of_multipoles) = record(1:2)
            record(1:2) = "  "
            goto 2
         case default
            number_of_multipoles = 0
            print *, "Unable to decode the transition multipole '"// &
                     record(1:2)//"'; reenter ..." 
            goto 1
         end select
      else if (number_of_multipoles == 0) then
         print *, "At least one transition multipole must be specified; "// &
                  "reenter ..."
         goto 1
      end if
      !
   end subroutine input_transition_multipoles
   !
   !
   subroutine input_transition_pairs(number_of_transitions)
   !--------------------------------------------------------------------
   ! Reads in and selects pairs of transitions.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(inout) :: number_of_transitions
      !
      integer            :: level_i, level_f, score_position
      logical            :: yes
      character(len=256) :: record
      character(len=20 ) :: string
      !
      print *, "Select individual transitions ?"
      yes = get_yes_stream()
      if (yes) then
    4    print *, "Enter one pair |i> - |f> of level numbers, e.g."// &
                  " 2 - 3; 2 - 0; 0 - 17 ... "
         print *, " (0 is here equivalent to all); <cr> if done."
         read (*, "(a)") record
         call save_input(record,.false.)
	 !
         if (len_trim(record) > 0) then
            score_position = scan(record,"-")
            if (score_position == 0) then
               print *, "Unable to decode the transition; reenter ..."
               goto 4
            else 
               string  = adjustl(record(1:score_position-1))
               level_i = get_integer_from_string(string)
               string  = adjustl(record(score_position+1:256))
               level_f = get_integer_from_string(string)
               print *, "level_i, level_f = ",level_i, level_f
               if (level_i < 0   .or.   level_f < 0) then
                  print *, "All level numbers must greater or equivalent 0;"//&
                           " reenter ...";   goto 4
               elseif (level_i == 0   .and.   level_f == 0) then
                  print *, "The pair of level numbers 0 - 0 is not allowed;"//&
                           " reenter ...";   goto 4
               endif
               !
               number_of_transitions = number_of_transitions + 1
               if (number_of_transitions > 200) then
                  stop "input_transition_pairs(): program stop A."
               end if
               select_level_i(number_of_transitions) = level_i
               select_level_f(number_of_transitions) = level_f
               goto 4
            end if
         else if (number_of_transitions == 0) then
            print *, "No pair of level numbers has been selected; the "// &
                     " program continues the default path."
         end if
      end if
      !
   end subroutine input_transition_pairs
   !
   !
   subroutine input_transition_triples(number_of_transitions)
   !--------------------------------------------------------------------
   ! Reads in and selects triples of transitions.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(inout) :: number_of_transitions
      !
      integer            :: level_i, level_m, level_f, score_position
      logical            :: yes
      character(len=256) :: record
      character(len=20 ) :: string
      !
      !
      print *, "Select individual transitions |i> --> |m> --> |f> ?"
      yes = get_yes_stream()
      if (yes) then
    6    print *, "Enter one triple |i> - |m> - |f> of level numbers, e.g."// &
                  " 2 - 3 - 1; 2 - 0 - 5; 0 - 17 - 3 ... "
         print *, " (0 is here equivalent to all); <cr> if done."
         read (*, "(a)") record
         if (len_trim(record) > 0) then
            score_position = scan(record,"-")
            if (score_position == 0) then
               print *, "Unable to decode the transition; reenter ..."
               goto 6
            else 
               string  = adjustl(record(1:score_position-1))
               level_i = get_integer_from_string(string)
               record  = adjustl(record(score_position+1:256))
               score_position = scan(record,"-")
               if (score_position == 0) then
                  print *, "Unable to decode the transition; reenter ..."
                  goto 6
               else 
                  string  = adjustl(record(1:score_position-1))
                  level_m = get_integer_from_string(string)
                  string  = adjustl(record(score_position+1:256))
                  level_f = get_integer_from_string(string)
               end if
               !
               print *, "level_i, level_m, level_f = ",level_i,level_m,level_f
               if (level_i < 0   .or.   level_m < 0   .or.   level_f < 0) then
                  print *, "All level numbers must greater or equivalent 0;"//&
                           " reenter ...";   goto 6
               elseif (level_i == 0  .and. level_i == 0 .and. level_f == 0) then
                  print *, "The pair of level numbers 0 - 0 - 0 is not"//&
                           " allowed; reenter ...";   goto 6
               endif
               !
               number_of_transitions = number_of_transitions + 1
               if (number_of_transitions > 200) then
                  stop "dierec_collect_input(): program stop B."
               end if
               select_level_i(number_of_transitions) = level_i
               select_level_m(number_of_transitions) = level_m
               select_level_f(number_of_transitions) = level_f
               goto 6
            end if
         else if (number_of_transitions == 0) then
            print *, "No pair of level numbers has been selected; the "// &
                     " program continues the default path."
         end if
      end if
      !
   end subroutine input_transition_triples
   !
end module rabs_input_dialog
