module rabs_toolbox_df
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains procedures which help carry out a number of
! 'utility' tasks within the RATIP environment. 
!-----------------------------------------------------------------------
   !
   use rabs_anco
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_input_dialog
   use rabs_lsj
   use rabs_lsj_data
   use rabs_toolbox_aux
   implicit none
   !
   public  :: toolbox_display_levels
                 ! Displays the energy levels and level splittings
                 ! (extended version).
   public  :: toolbox_display_levels_lsj
                 ! Displays the energy levels and level splittings
                 ! from one .inp and .mix file together with the corresponding
		 ! (two) leading LSJ-coupled configurations.
   public  :: toolbox_display_levels_simple
                 ! Displays the energy levels and level splittings
                 ! (simple version).
   public  :: toolbox_display_weights
                 ! Displays the major weight contributions in a jj-coupled 
                 ! CSF basis.
   public  :: toolbox_excl_generate_pairlist
                 ! Generates a pair-correlation list from a .csl list with
                 ! respect to a given set of reference CSF.
   public  :: toolbox_exclude_csf
                 ! Excludes CSF from a .csl list due to 'given rules'.
   private :: toolboxe_excl_apply_restrict
                 ! Applies a given set of 'restriction rules' to the
                 ! currently accepted .csl list.
   private :: toolbox_excl_communicate_restr
                 ! Reads in and applies the 'restriction rules' to the 
                 ! currently accepted .csl list.
   private :: toolbox_excl_get_reference_list
                 ! Opens a .csl file, checks it and loads the list of
                 ! CSF as the reference list for generating a pair-correlation
                 ! expansion.
   private :: toolbox_excl_interpret_restr
                 ! Attempts to interpret a given restriction rules and defines
                 ! and corresponding entry in the list of rules.
   private :: toolbox_excl_pair_correlation
                 ! Determines which CSF are not part of a pair-correlation list
                 ! with respect to a given reference list.
   private :: toolbox_exclude_write_csl
                 ! Writes the 'accepted' CSF from csl_old_file to a new 
                 ! .csl file.
   public  :: toolbox_format_mix
                 ! Formats a GRASP92 .mix  file.
   public  :: toolbox_format_out
                 ! Formats a GRASP92 .out  file.
   !
   ! Data structure for the 'exclude' procedures
   type(csf_basis), public :: exclude_csf_set, reference_csf_set, ref_csf_set
   !
   logical, dimension(:), allocatable, public :: exclude_accepted, &
                                                 exclude_prelim
   !
   type, public :: restriction_rule
      integer          :: number_of_terms, limitation
      character(len=2) :: relation
      type(nkappa), dimension(:), pointer :: term
   end type restriction_rule
   !
   integer, public :: number_of_restrictions = 0
   type(restriction_rule), dimension(20), public :: restriction
   !
   !
contains
   !
   subroutine toolbox_df()
   !--------------------------------------------------------------------
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      print *, "**************************"
      print *, "*** Not yet implmented ***"
      print *, "**************************"
      !
   end subroutine toolbox_df
   !
   !
   subroutine toolbox_display_levels()
   !--------------------------------------------------------------------
   ! Prints out the level energies, the level splitting as well as
   ! the excitation energies with respect to a given level. This 'extended 
   ! version' enables the user to select the energy unit, the order of 
   ! printing as well as the level which is taken as the 'zero' reference
   ! level. The information about the levels must be supplied by one or 
   ! several .mix Mixing Coefficient Files from GRASP92 or RELCI.
   !
   ! Calls: get_yes_stream(), toolbox_gather_energies().
   !--------------------------------------------------------------------
      !
      integer            :: ce, i, j, no_files, no_levels, zlevel
      character(len=32)  :: record
      logical            :: yes, ascending_order
      real(kind=dp)      :: current, last, splitting, excitation
      integer, dimension(:,:), allocatable :: history
      !
      ! Determine the units for the printout and further optional input
      call input_energy_unit()
      !
      !
      print *, "Energies are printed in ascending order; " //&
               "use a decending order instead ?"
      yes = get_yes_stream()
      if (yes) then
         ascending_order = .false.
      else
         ascending_order = .true.
      end if
      !
      print *, "Enter the level number of the reference (zero-) level:"
      read *, zlevel
      !
      call toolbox_gather_energies(no_files)
      !
      no_levels = 0
      do  i = 1,no_files
         no_levels = no_levels + asf_mix(i)%noasf
      end do
      allocate( asf_set%asf(1:no_levels),   history(1:no_levels,1:2) )
      !
      if (zlevel <= 0   .or.   zlevel > no_levels) then
         zlevel = 1
         print *, " "
         print *, "The reference level does not match with any of the given "//&
                  "levels; instead, the lowest is taken for reference."
      end if
      !
      ! Now order all energies and print out the results
      ce   = 0  ! current energy
      last = -1.0e99_dp
      !
    4 current =  1.0e99_dp
      !
      do  i = 1,no_files
         do  j = 1,asf_mix(i)%noasf
            if (asf_mix(i)%asf(j)%energy > last  .and.  &
                asf_mix(i)%asf(j)%energy < current)  then
               current = asf_mix(i)%asf(j)%energy
            end if
         end do
      end do
      !
      do  i = 1,no_files
         do  j = 1,asf_mix(i)%noasf
            if (asf_mix(i)%asf(j)%energy == current) then
               ce = ce + 1
               asf_set%asf(ce)%level_No = ce
               asf_set%asf(ce)%totalJ   = asf_mix(i)%asf(j)%totalJ
               asf_set%asf(ce)%parity   = asf_mix(i)%asf(j)%parity
               asf_set%asf(ce)%energy   = asf_mix(i)%asf(j)%energy
               ! Keep information about original files
               history(ce,1)            = i
               history(ce,2)            = asf_mix(i)%asf(j)%level_No
            end if
         end do
      end do
      !
      last = current
      if (ce < no_levels) goto 4
      !
      ! Now print the energy levels and level splittings in a neat format
      !
      if (ascending_order) then
         print *,   " "
         write(*,8)
         write(*,6) energy_unit, energy_unit
       6 format(   16x,"                      Total             Level    ",  &
                       "   Excitation energy  ",                             &
                  /16x,"     [File /          energy          splitting  ",  &
                       "   from the reference  ",                            &
                  /4x,                                                       &
           "Level  J Parity  [Level]        (Hartrees)         (",a ,")  ",  &
                       "       (",a ,")" )
         write(*,8)
         !
         do  i = 1,no_levels
            splitting = zero
            if (i > 1) splitting = (asf_set%asf(i)%energy -                  &
                            asf_set%asf(i-1)%energy)
            excitation = (asf_set%asf(i)%energy-asf_set%asf(zlevel)%energy) 
            if (energy_inverse) then
               splitting  = energy_factor / splitting
               excitation = energy_factor / excitation
            else
               splitting  = energy_factor * splitting
               excitation = energy_factor * excitation
            end if
            !                       
            write(*,7)                                                       &
               asf_set%asf(i)%level_No,                                      &
               trim(angular_momentum_string(1*asf_set%asf(i)%totalJ,4)),     &
               asf_set%asf(i)%parity,history(i,1),history(i,2),              &
               asf_set%asf(i)%energy,splitting,excitation
         7 format(4x,1i3,2x,2a4,4x,"[",i2,",",i3,"]",1x,1p,3e18.8)
         end do
         write(*,8)
       8 format(3x,83("-"))
      else
         print *,   " "
         write(*,8)
         write(*,6) energy_unit, energy_unit
         write(*,8)
         !
         do  i = no_levels,1,-1
            splitting = zero
            if (i < no_levels) splitting = (asf_set%asf(i)%energy -          &
                            asf_set%asf(i+1)%energy)
            excitation = (asf_set%asf(i)%energy-asf_set%asf(zlevel)%energy) 
            if (energy_inverse) then
               splitting  = energy_factor / splitting
               excitation = energy_factor / excitation
            else
               splitting  = energy_factor * splitting
               excitation = energy_factor * excitation
            end if
            !                       
            write(*,7)                                                       &
               asf_set%asf(i)%level_No,                                      &
               trim(angular_momentum_string(1*asf_set%asf(i)%totalJ,4)),     &
               asf_set%asf(i)%parity,history(i,1),history(i,2),              &
               asf_set%asf(i)%energy,splitting,excitation
         end do
         write(*,8)
      end if
      !
   end subroutine toolbox_display_levels
   !
   !
   subroutine toolbox_display_levels_lsj()
   !--------------------------------------------------------------------
   ! Prints out the level energies, level splittings, and the excitation
   ! energies from one set of .inp and .mix files with respect to the lowest,
   ! together with the corresponding (two) leading LSJ-coupled 
   ! configurations.
   
   ! Calls: toolbox_gather_energies().
   !--------------------------------------------------------------------
      !
      integer            :: i, ierr, ios, k, k1, k2, k3, k4, k5
      real(kind=dp)      :: excitation, splitting
      character(len=6)   :: g92mix
      character(len=256) :: toolbox_mix_file
      !
      type(lsj_level_weight), dimension(100,9) :: asf_weight
      type(lsj_string_csf)  , dimension(30000) :: asf_string
      !
      ! Read in a .csl and .mix file for further processing
      ! Get the occupied shells of the initial and final states
      call file_get_csl_list(                                       &
         "Enter the name of the configuration symmetry list file:", &
         asf_set_jj%csf_set)
      !
    1 print *, "Enter the name of corresponding .mix mixing coefficient file:"
      read (*,"(a)") toolbox_mix_file
      if (len(trim(toolbox_mix_file)) == 0) goto 1
      !
      call file_open(25,toolbox_mix_file,"unformatted","old",ierr)
      if (ierr /= 0) goto 1
      !
      ! Check the header of the file; if not as expected for unformatted
      ! files, check formatted form
      read (25,iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         print *, "ios, g92mix = ",ios, g92mix
         close (25)
         goto 2
      else
         call load_mix_file_grasp2k(asf_set,.false.,ierr)
         if (ierr /= 0) then
            print *, "Not a proper .mix mixing coefficient file for the "//&
                     "given .csl list; reenter ..."
            close (25)
            goto 1
         end if
         close (25)
         goto 1
      end if
      !
      ! Try formatted file format; check the header of the file; 
      ! if not as expected for formatted files, try again
    2 call file_open(25,toolbox_mix_file,"formatted  ","old",ierr)
      if (ierr /= 0) goto 2
      read (25,"(a)",iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         print *, "ios, g92mix = ",ios, g92mix
         print *, "Not a GRASP92 Mixing Coefficients File;"
         close (25)
         goto 1
      else
         call load_mix_file_grasp2k(asf_set_jj,.true.,ierr)
         if (ierr /= 0) then
            close (25)
            goto 1
         end if
         close (25)
      end if
      !
      ! Next, perform an LSJ expansion of the wave functions
      call lsj_initialization_LS_jj()
      call lsj_control_from_toolbox(toolbox_mix_file,asf_weight,asf_string)   
      !
      !
      ! Now print the energy levels and level splittings in a neat format
      write(*,*) " "
      write(*,*) " "
      print *,   " "
      write(*,8)
      write(*,6)
    6 format(   16x,"            Total             Level    ",                &
                    "   Excitation energy  ",                                 &
               /16x,"            energy          splitting  ",                &
                    "    from the lowest",                                    &
               /4x,                                                           &
        "Level  J Parity       (Hartrees)         (Kayser)  ",                &
                    "       (Kayser)")
      write(*,8)
      !
      do  i = 1,asf_set_jj%noasf
         splitting = zero
         if (i > 1) splitting = (asf_set_jj%asf(i)%energy -                      &
                            asf_set_jj%asf(i-1)%energy) *convert_au_to_kaysers
         excitation = (asf_set_jj%asf(i)%energy-asf_set_jj%asf(1)%energy)           &
                                * convert_au_to_kaysers
         write(*,7)                                                           &
            asf_set_jj%asf(i)%level_No,                                          &
            trim(angular_momentum_string(1*asf_set_jj%asf(i)%totalJ,4)),         &
            asf_set_jj%asf(i)%parity,                  &
            asf_set_jj%asf(i)%energy,splitting,excitation
      7 format(4x,1i3,2x,2a4,3x,1p,3e18.8)
      end do
      write(*,8)
    8 format(3x,83("-"))
      !
      !
      write(*,*) " "
      write(*,*) " "
      print *,   " "
      write(*,18)
      write(*,6)
      write(*,18)
      !
      if (asf_set_jj%noasf > 30000) then
         stop "toolbox_display_levels_lsj(): program stop A." 
      end if
      !
      do  i = 1,asf_set_jj%noasf
         splitting = zero
         if (i > 1) splitting = (asf_set_jj%asf(i)%energy -                      &
                            asf_set_jj%asf(i-1)%energy) *convert_au_to_kaysers
         excitation = (asf_set_jj%asf(i)%energy-asf_set_jj%asf(1)%energy)           &
                                * convert_au_to_kaysers
         write(*,7)                                                           &
            asf_set_jj%asf(i)%level_No,                                          &
            trim(angular_momentum_string(1*asf_set_jj%asf(i)%totalJ,4)),         &
            asf_set_jj%asf(i)%parity,                  &
            asf_set_jj%asf(i)%energy,splitting,excitation
	 !
	 k1 = 0;   k2 = 0;   k3 = 0;   k4 = 0;   k5 = 0
	 do  k = 1,asf_set_jj%csf_set%nocsf
	    if (asf_weight(i,1)%csf == k) k1 = k
	    if (asf_weight(i,2)%csf == k) k2 = k
	    if (asf_weight(i,3)%csf == k) k3 = k
	    if (asf_weight(i,4)%csf == k) k4 = k
	    if (asf_weight(i,5)%csf == k) k5 = k
	 end do
	 !
         write(*,10) asf_weight(i,1)%weight,asf_string(k1)%s1(12:100),    &
	                                    asf_string(k1)%s2(12:100)
	 if (k2 /= 0)                                                     &				    
         write(*,10) asf_weight(i,2)%weight,asf_string(k2)%s1(12:100),    &
	                                    asf_string(k2)%s2(12:100)
	 if (k3 /= 0)                                                     &				    
         write(*,10) asf_weight(i,3)%weight,asf_string(k3)%s1(12:100),    &
	                                    asf_string(k3)%s2(12:100)
	 if (k4 /= 0)                                                     &				    
         write(*,10) asf_weight(i,4)%weight,asf_string(k4)%s1(12:100),    &
	                                    asf_string(k4)%s2(12:100)
	 if (k5 /= 0)                                                     &				    
         write(*,10) asf_weight(i,5)%weight,asf_string(k5)%s1(12:100),    &
	                                    asf_string(k5)%s2(12:100)
      end do
      write(*,18)
   10 format(16x,f6.4,"  of  ",a,  /,28x,a)
   18 format(3x,103("-"))
      !
      !
   end subroutine toolbox_display_levels_lsj
   !
   !
   subroutine toolbox_display_levels_simple()
   !--------------------------------------------------------------------
   ! Prints out the level energies, level splittings, and the excitation
   ! energies with respect to the lowest. This 'simple version' only
   ! requires the input of one or several .mix files.
   !
   ! Calls: toolbox_gather_energies().
   !--------------------------------------------------------------------
      !
      integer       :: ce, i, j, no_files, no_levels
      real(kind=dp) :: current, last, splitting, excitation
      integer, dimension(:,:), allocatable :: history
      !
      call toolbox_gather_energies(no_files)
      !
      no_levels = 0
      do  i = 1,no_files
         no_levels = no_levels + asf_mix(i)%noasf
      end do
      allocate( asf_set%asf(1:no_levels),   history(1:no_levels,1:2) )
      !
      ! Now order all energies and print out the results
      ce   = 0  ! current energy
      last = -1.0e99_dp
      !
    4 current =  1.0e99_dp
      !
      do  i = 1,no_files
         do  j = 1,asf_mix(i)%noasf
            if (asf_mix(i)%asf(j)%energy > last  .and.  &
                asf_mix(i)%asf(j)%energy < current)  then
               current = asf_mix(i)%asf(j)%energy
            end if
         end do
      end do
      !
      do  i = 1,no_files
         do  j = 1,asf_mix(i)%noasf
            if (asf_mix(i)%asf(j)%energy == current) then
               ce = ce + 1
               asf_set%asf(ce)%level_No = ce
               asf_set%asf(ce)%totalJ   = asf_mix(i)%asf(j)%totalJ
               asf_set%asf(ce)%parity   = asf_mix(i)%asf(j)%parity
               asf_set%asf(ce)%energy   = asf_mix(i)%asf(j)%energy
               ! Keep information about original files
               history(ce,1)            = i
               history(ce,2)            = asf_mix(i)%asf(j)%level_No
            end if
         end do
      end do
      !
      last = current
      if (ce < no_levels) goto 4
      !
      ! Now print the energy levels and level splittings in a neat format
      print *,   " "
      write(*,8)
      write(*,6)
    6 format(   16x,"                      Total             Level    ",    &
                    "   Excitation energy  ",                               &
               /16x,"     [File /          energy          splitting  ",    &
                    "    from the lowest",                                  &
               /4x,                                                         &
        "Level  J Parity  [Level]        (Hartrees)         (Kayser)  ",    &
                    "       (Kayser)")
      write(*,8)
      !
      do  i = 1,no_levels
         splitting = zero
         if (i > 1) splitting = (asf_set%asf(i)%energy -                      &
                                 asf_set%asf(i-1)%energy) *convert_au_to_kaysers
         excitation = (asf_set%asf(i)%energy-asf_set%asf(1)%energy)           &
                                * convert_au_to_kaysers
         write(*,7)                                                           &
            asf_set%asf(i)%level_No,                                          &
            trim(angular_momentum_string(1*asf_set%asf(i)%totalJ,4)),         &
            asf_set%asf(i)%parity,history(i,1),history(i,2),                  &
            asf_set%asf(i)%energy,splitting,excitation
      7 format(4x,1i3,2x,2a4,4x,"[",i2,",",i3,"]",1x,1p,3e18.8)
      end do
      write(*,8)
   8 format(3x,83("-"))
      !
   end subroutine toolbox_display_levels_simple
   !
   !
   subroutine toolbox_display_weights()
   !--------------------------------------------------------------------
   ! Displays the main CSF and their weights for one or several levels
   ! of interest. The .csl and .mix files are required.
   !
   ! Calls: get_yes_stream()
   !--------------------------------------------------------------------
      !
      integer            :: ierr, ios, i, j, lev, number_of_levels, &
                            number_of_leading_csf, nocsf_min
      character(len=6)   :: g92mix
      character(len=120) :: record
      character(len=256) :: toolbox_csl_file, toolbox_mix_file
      logical            :: yes, fail, is_selected
      real(kind=dp)      :: wa, wb  
      !
      integer, dimension(1:100)       :: iw
      integer, dimension(1:1000)      :: levels    
      logical, dimension(1:100000)    :: csf   
      real(kind=dp), dimension(1:100) :: weights
      !
      number_of_leading_csf = 5
      number_of_levels      = 0
      levels(:)             = 0
      !
      print *, "Select individual levels or a particular number of " //&
               "leading CSF ? " 
      yes = get_yes_stream()
      if (yes) then
    1    print *, "Enter the individual level numbers to be selected:"
         print *, " e.g. 1 3 4  7 - 20  48  69 - 85;"
         read (*, "(a)") record
         call toolbox_interprete_levels(record,levels,number_of_levels,fail)
         if (fail) then
            print *, "Unable to interprete the serial level numbers; redo ..."
            goto 1
         end if
         !
         print *, "Enter the (maximal) number of leading CSF to be printed:"
         read *, number_of_leading_csf
      end if
      !
      ! Open, check, load data from, and close, the  .csl  file
    2 print *, "Enter the name of the GRASP92 configuration "//&
               "symmetry list file:"
      read (*,"(a)") toolbox_csl_file
      if (len(trim(toolbox_csl_file)) == 0) goto 2
      !
      call file_open(21,toolbox_csl_file,"formatted  ","old",ierr)
      if (ierr == 1) goto 2
      !
      ! Check the first record of the file; if not as expected, try again
      read (21,"(1a15)",iostat = ios) record
      if (ios /= 0   .or.   record(1:15) /= "Core subshells:") then
         print *, "ios, record(1:15) = ",ios, record(1:15)
         print *, "Not a configuration symmetry list file;"
         close (21)
         goto 2
      end if
      !
      ! Load data from the  .csl  file
      call load_csl_from_file(asf_set%csf_set)
      !
      ! Close the  .csl  file
      close (21)
      !
    3 print *, "Enter the name of .mix mixing coefficient file:"
      read (*,"(a)") toolbox_mix_file
      if (len(trim(toolbox_mix_file)) == 0) goto 3
      !
      call file_open(25,toolbox_mix_file,"unformatted","old",ierr)
      if (ierr /= 0) goto 3
      !
      ! Check the header of the file; if not as expected for unformatted
      ! files, check formatted form
      read (25,iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         close (25)
         goto 4
      else
         call load_mix_file_grasp2k(asf_set,.false.,ierr)
         if (ierr /= 0) then
            close (25)
            goto 3
         end if
         goto 5
      end if
      !
      ! Try formatted file format; check the header of the file; 
      ! if not as expected for formatted files, try again
    4 call file_open(25,toolbox_mix_file,"formatted  ","old",ierr)
      if (ierr /= 0) goto 3
      read (25,"(a)",iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         print *, "ios, g92mix = ",ios, g92mix
         print *, "Not a GRASP92 Mixing Coefficients File;"
         close (25)
         goto 3
      else
         call load_mix_file_grasp2k(asf_set,.true.,ierr)
      end if
      !
    5 csf(:) = .false.
      !
      ! Now print the weights and the explicit coupling of the leading CSF
      print *, " "
      print *, "Weights of major contributors to ASF:"
      print *, " "
      print *, " Level  J Parity      CSF contributions"
      print *, " "
      !
      do  lev = 1,asf_set%noasf
         !
         ! Check that information about this CSF should be printed
         if (number_of_levels == 0) then
            is_selected = .true.
         else
            is_selected = .false.
            do  i = 1,number_of_levels
               if (levels(i) == asf_set%asf(lev)%level_No) then
                  is_selected = .true.
                  exit
               end if
            end do
         end if
         !
         if (is_selected) then
            weights(1:100) = zero;   iw(1:100) = 0
            wb = zero
            do  i = 1,asf_set%csf_set%nocsf
               wa = asf_set%asf(lev)%eigenvector(i) * &
                    asf_set%asf(lev)%eigenvector(i)
               wb = wb + asf_set%asf(lev)%eigenvector(i) * &
                         asf_set%asf(1)%eigenvector(i)
               do  j = 1,99
                  if (wa > weights(j)) then
                     weights(j+1:100) = weights(j:99)
                     weights(j)       = wa
                     iw(j+1:100)      = iw(j:99)
                     iw(j)            = i
                     exit
                  end if
               end do
               ! print *, "iw(:) = ",iw
            end do
            !
            if (lev > 1   .and.   abs(wb) > 0.0001) then
               print *, "level, wb = ",lev,wb
               stop "toolbox_print_weights(): program stop A." 
            end if
            !
            nocsf_min = number_of_leading_csf
            do  j = 1,number_of_leading_csf
               if (weights(j) == zero) then
                  nocsf_min = j - 1
                  exit
               end if
            end do
            !
            print 6, asf_set%asf(lev)%level_No,                              &
                 trim(angular_momentum_string(1*asf_set%asf(lev)%totalJ,4)), &
                      asf_set%asf(lev)%parity,(weights(j),iw(j),j=1,nocsf_min)
            !
            do  j = 1,nocsf_min 
               csf(iw(j)) = .true.
            end do
         end if
      end do
    6 format(1x,i4,3x,2a4,100(3x,f8.5," of",i5))
      !
      ! Print the CSF explicitly
      call file_open(21,toolbox_csl_file,"formatted  ","old",ierr)
      read(21,*);   read(21,*);   read(21,*);   read(21,*);   read(21,*) 
      !
      print *, " "
      print *, "Definition of leading CSF:"
      print *, " "
      !
      do  j = 1,asf_set%csf_set%nocsf
         if (csf(j)) then
            read(21,"(a)") record
            i = len(trim(record))
            write(*,7) j, record 
            read(21,"(a)") record
            i = len(trim(record))
            write(*,8) record 
            read(21,"(a)") record
            i = len(trim(record))
            write(*,8) record 
         else
            read(21,*);   read(21,*);   read(21,*) 
         end if
      end do
      !
      close(21)
    7 format(i6,")"   ,a)
    8 format("       ",a)
      !
   end subroutine toolbox_display_weights
   !
   !
   subroutine toolbox_excl_generate_pairlist()
   !--------------------------------------------------------------------
   ! Utility program to generates a pair-correlation list in order to
   ! exclude all CSF in a given list which are not coupled by the
   ! Hamiltonian to some set of reference configuration.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer            :: i, nocsf, current_number, previous_number
      character(len=256) :: csl_old_file
      !
      ! Open, check, load data from, and close, the  .csl  file
      call file_get_csl_list(                                       &
         "Enter the name of the configuration symmetry list file:", &
         exclude_csf_set,csl_old_file)
      !
      nocsf = exclude_csf_set%nocsf
      allocate( exclude_accepted(1:nocsf), exclude_prelim(1:nocsf) )
      exclude_accepted = .true.;   exclude_prelim = .true.
      !
      ! Generate the pair-correlation list with respect to a given reference
      ! set of CSF.
      call toolbox_excl_get_reference_list(reference_csf_set)
      call toolbox_excl_pair_correlation()
      !!x end if
      !
      !  
      call toolbox_exclude_write_csl(csl_old_file) 
      current_number = 0; previous_number = 0
      do  i = 1,nocsf
         if (exclude_prelim(i))   current_number  = current_number + 1
         if (exclude_accepted(i)) previous_number = previous_number + 1
      end do
      !
   end subroutine toolbox_excl_generate_pairlist
   !
   !
   subroutine toolbox_exclude_csf()
   !--------------------------------------------------------------------
   ! Utility program to 'exclude' CSF from a given GRASP92 .csl file;
   ! it applies a number of user-supplied 'restriction rules' in order 
   ! to exclude unnecessary CSFs in a 'physically' well-defined way.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer            :: i, nocsf, current_number, previous_number
      character(len=256) :: csl_old_file
      !
      ! Open, check, load data from, and close, the  .csl  file
      call file_get_csl_list(                                       &
         "Enter the name of the configuration symmetry list file:", &
         exclude_csf_set,csl_old_file)
      !
      nocsf = exclude_csf_set%nocsf
      allocate( exclude_accepted(1:nocsf), exclude_prelim(1:nocsf) )
      exclude_accepted = .true.;   exclude_prelim = .true.
      !
      call toolbox_excl_communicate_restr(current_number,previous_number)
      !  
      call toolbox_exclude_write_csl(csl_old_file) 
      current_number = 0; previous_number = 0
      do  i = 1,nocsf
         if (exclude_prelim(i))   current_number  = current_number + 1
         if (exclude_accepted(i)) previous_number = previous_number + 1
      end do
      !
   end subroutine toolbox_exclude_csf
   !
   !
   subroutine toolboxe_excl_apply_restrict(relativistic)
   !--------------------------------------------------------------------
   ! Applies a given set of 'restriction rules' to the currently accepted 
   ! .csl list.
   !--------------------------------------------------------------------
      !
      logical, intent(in) :: relativistic
      !
      integer :: i, icsf, iorb, j, kappa, pqn, particle_number
      logical :: all_shells_contribute
      !
      do  icsf = 1,exclude_csf_set%nocsf
         if (exclude_accepted(icsf)) then
            do  i = 1,number_of_restrictions
               particle_number = 0
               all_shells_contribute = .true.
               do  j = 1,restriction(i)%number_of_terms
                  pqn   = restriction(i)%term(j)%n
                  kappa = restriction(i)%term(j)%kappa
                  do  iorb = exclude_csf_set%nwcore+1,exclude_csf_set%nwshells
                     if (pqn   == exclude_csf_set%subshell(iorb)%n  .and. &
                         kappa == exclude_csf_set%subshell(iorb)%kappa) then
                        particle_number = particle_number + &
                           exclude_csf_set%csf(icsf)%occupation(iorb)
                        if (relativistic) then
                           if (exclude_csf_set%csf(icsf)%occupation(iorb)==0) &
                                                                           then
                              all_shells_contribute = .false.
                           end if
                        else
                           if (                                        &
                           exclude_csf_set%csf(icsf)%occupation(iorb)==0) then
                              all_shells_contribute = .false.
                           end if
                        end if
                        exit
                     end if
                  end do
               end do
               !
               ! Test total particle number vs. the given limitation
               select case(restriction(i)%relation)
               case("< ")
                  if (.not.(particle_number <  restriction(i)%limitation)) then
                     exclude_prelim(icsf) = .false.; exit
                  end if
               case("<=")
                  if (.not.(particle_number <= restriction(i)%limitation)) then
                     exclude_prelim(icsf) = .false.; exit
                  end if
               case("<$")
                  if (.not.(particle_number <  restriction(i)%limitation) .and.&
                      all_shells_contribute) then
                     exclude_prelim(icsf) = .false.; exit
                  end if
               case("= ")
                  if (.not.(particle_number == restriction(i)%limitation)) then
                     exclude_prelim(icsf) = .false.; exit
                  end if
               case("> ")
                  if (.not.(particle_number >  restriction(i)%limitation)) then
                     exclude_prelim(icsf) = .false.; exit
                  end if
               case(">=")
                  if (.not.(particle_number >= restriction(i)%limitation)) then
                     exclude_prelim(icsf) = .false.; exit
                  end if
               case default
                  stop "toolboxe_excl_apply_restrict(): program stop A."
               end select
            end do
         end if
      end do 
      !
   end subroutine toolboxe_excl_apply_restrict
   !
   !
   subroutine toolbox_excl_communicate_restr(current_number, &
                                                   previous_number)
   !--------------------------------------------------------------------
   ! Reads in and applies a given set of 'restriction rules' to the 
   ! currently accepted .csl list.
   !--------------------------------------------------------------------
      !
      integer, intent(out) :: current_number, previous_number
      !
      integer            :: i
      logical            :: yes, relativistic, fail
      character(len=100) :: record
      !
      print *, " "
      print *, "Use a nonrelativistic notation for the orbitals ?"
      yes = get_yes_stream()
      if (yes) then;   relativistic = .false.
      else;            relativistic = .true.
      end if
      !
    1 print *, "Enter another restriction rule (null if done):"
      read(*,"(a)") record
      !!x print *, trim(record)//":"
      if (len_trim(record) /= 0) then
         call toolbox_excl_interpret_restr(record,relativistic,fail)
         if (fail) then
            print *, "Unable to interpret string '"//trim(record)//"' as"// &
                     " valid restriction rule; redo ..."
         end if
         goto 1
         if (number_of_restrictions == 20) then
            stop "exclude_communicate_restrictions(): program stop A."
         end if
      end if
      if (number_of_restrictions == 0) return
      !
      ! Process restrictions on the (still) accepted CSF list
      call toolboxe_excl_apply_restrict(relativistic)
      current_number = 0; previous_number = 0
      do  i = 1,exclude_csf_set%nocsf
         if (exclude_prelim(i))   current_number  = current_number + 1
         if (exclude_accepted(i)) previous_number = previous_number + 1
      end do
      print *, "current, previous = ",current_number,previous_number
      !
      print *, "The current .csl list includes ",current_number, &
               " relativistic CSF;"
      print *, " Apply this reduction of the CSFs ?"
      yes = get_yes_stream()
      if (yes) then
         exclude_accepted(1:exclude_csf_set%nocsf) =  &
                            exclude_prelim(1:exclude_csf_set%nocsf)
         do  i = 1,number_of_restrictions
            deallocate(restriction(i)%term)
         end do
         number_of_restrictions = 0
         goto 1
      else
         exclude_prelim(:) = exclude_accepted(:)
         do  i = 1,number_of_restrictions
            deallocate(restriction(i)%term)
         end do
         number_of_restrictions = 0
         print *, "The CSF list has not been modified by the last"//&
                  " restriction rules;"
         print *, " currently, there are ",previous_number," relativistic CSFs;"
         print *, " proceed with further rules ?"
         goto 1
      end if
      !
   end subroutine toolbox_excl_communicate_restr
   !
   !
   subroutine toolbox_excl_get_reference_list(csf_set)
   !--------------------------------------------------------------------
   ! Opens a .csl file, checks it format and number of functions, and 
   ! loads the list of CSF into the internal arrays. This .csl list should
   ! represent the reference list with respect to which the pair-correlation
   ! list will be generated. 
   ! This file is always attached to stream 21.
   !
   ! Calls: load_csl_from_file(), file_open().
   !--------------------------------------------------------------------
      !
      type(csf_basis), intent(inout) :: csf_set
      !
      integer            :: ierr, ios
      character(len=15)  :: record
      character(len=256) :: csl_list_file
      !
    1 print *, "Enter the name of the configuration symmetry list file"//&
               " for the reference configurations:"
      read (*,"(a)") csl_list_file
      if (len(trim(csl_list_file)) == 0) goto 1
      !
      call file_open(21,csl_list_file,"formatted  ","old",ierr)
      if (ierr == 1) goto 1
      !
      ! Check the first record of the file; if not as expected, try again
      read (21,"(1a15)",iostat = ios) record
      if (ios /= 0   .or.   record(1:15) /= "Core subshells:") then
         print *, "ios, record(1:15) = ",ios, record(1:15)
         print *, "Not a configuration symmetry list file;"
         close (21)
         goto 1
      end if
      !
      ! Load data from the  .csl  file
      call load_csl_from_file(csf_set)
      !
      ! Close the  .csl  file
      close (21)
      !
   end subroutine toolbox_excl_get_reference_list
   !
   !
   subroutine toolbox_excl_interpret_restr(record,relativistic,fail)
   !--------------------------------------------------------------------
   ! Attempts to interpret a given restriction rules and defines a
   ! corresponding entry in the list of rules.
   !--------------------------------------------------------------------
      !
      character(len=*), intent(in) :: record
      logical, intent(in)          :: relativistic
      logical, intent(out)         :: fail
      !
      integer :: i, il, it, j, kappa, limitation, pqn, term_number
      logical :: fail1
      character(len=2)   :: relation
      character(len=100) :: reca, recb
      type(nkappa), dimension(100) :: term
      !
      fail = .true.
      !
      ! First find the relation operator and the limitation
      i = scan(record,"<")
      if (i /= 0) then
         if (record(i:i+1) == "<=") then
            relation = "<="; il = i+2; goto 1
         else if (record(i:i+1) == "<$") then
            relation = "<$"; il = i+2; goto 1
         else 
            relation = "< "; il = i+1; goto 1
         end if
      end if
      i = scan(record,">")
      if (i /= 0) then
         if (record(i:i+1) == ">=") then
            relation = ">="; il = i+2; goto 1
         else
            relation = "> "; il = i+1; goto 1
         end if
      end if
      i = scan(record,"=")
      if (i /= 0) then
         relation = "= "; il = i+1; goto 1
      end if
      !!x print *, "relation, limitation = ",relation, limitation
      !
      ! No proper relation operator found; return with fail
      return
      !
    1 reca = adjustl(record(il:))
      limitation = get_integer_from_string(reca(1:2),fail1)
      !!x print *, "relation, limitation = ",relation, limitation
      if (fail1) return
      !
      ! Determine the number of terms
      reca = record(1:i-1); term_number = 0
      do
         it = scan(reca,"+")
         if (it /= 0) then
            recb = adjustl(reca(1:it-1))
            if (recb(2:2) == "*"   .or.   recb(3:3) == "*") then
               kappa = 0
               pqn   = get_integer_from_string("   "//recb(1:1),fail1)
            else if (recb(3:3) == " "   .or.   recb(3:3) == "-") then
               kappa = get_kappa_from_name(recb(2:3),fail1)
               pqn   = get_integer_from_string("   "//recb(1:1),fail1)
            else 
               kappa = get_kappa_from_name(recb(3:4),fail1)
               pqn   = get_integer_from_string("  "//recb(1:2),fail1)
            end if
            !!x print *, "it, recb, fail1 = ",it, recb, fail1
            if (fail1) return
            if (.not.relativistic   .and.   kappa /= 0) then
               term_number = term_number + 1
               term(term_number)%kappa = kappa;   term(term_number)%n = pqn
               if (kappa /= -1) then
                  term_number = term_number + 1
                  term(term_number)%kappa = abs(kappa) - 1   
                  term(term_number)%n     = pqn
               end if
            else if (relativistic   .and.   kappa /= 0) then
               term_number = term_number + 1
               term(term_number)%kappa = kappa;   term(term_number)%n = pqn
            else
               do  j = -pqn,pqn-1
                  if (j == 0) cycle
                  term_number = term_number + 1
                  term(term_number)%kappa = j;   term(term_number)%n = pqn
               end do
            end if
            reca = reca(it+1:)
         else
            recb = adjustl(reca(it+1:))
            if (recb(2:2) == "*"   .or.   recb(3:3) == "*") then
               kappa = 0
               pqn   = get_integer_from_string("   "//recb(1:1),fail1)
            else if (recb(3:3) == " "   .or.   recb(3:3) == "-") then
               kappa = get_kappa_from_name(recb(2:3),fail1)
               pqn   = get_integer_from_string("   "//recb(1:1),fail1)
            else 
               kappa = get_kappa_from_name(recb(3:4),fail1)
               pqn   = get_integer_from_string("  "//recb(1:2),fail1)
            end if
            print *, "it, recb, fail1 = ",it, recb, fail1
            if (fail1) return
            if (.not.relativistic   .and.   kappa /= 0) then
               term_number = term_number + 1
               term(term_number)%kappa = kappa;   term(term_number)%n = pqn
               if (kappa /= -1) then
                  term_number = term_number + 1
                  term(term_number)%kappa = abs(kappa) - 1   
                  term(term_number)%n     = pqn
               end if
            else if (relativistic   .and.   kappa /= 0) then
               term_number = term_number + 1
               term(term_number)%kappa = kappa;   term(term_number)%n = pqn
            else
               do  j = -pqn,pqn-1
                  if (j == 0) cycle
                  term_number = term_number + 1
                  term(term_number)%kappa = j;   term(term_number)%n = pqn
               end do
            end if
            !
            exit
         end if
      end do
      !
      ! Ensure that no orbital appears twice in the list
      do i = 1,term_number-1
         do j = i+1,term_number
            if (term(i)%n     == term(j)%n  .and.  &
                term(i)%kappa == term(j)%kappa) then
               print *, " Orbitals may not appear more than ones in each"//&
                        " restriction rule;"
               fail = .true.;   return
            end if
         end do
      end do
      !
      ! Store this restriction
      number_of_restrictions = number_of_restrictions + 1
      restriction(number_of_restrictions)%number_of_terms = term_number
      allocate( restriction(number_of_restrictions)%term(term_number) )
      do  i = 1,term_number
         restriction(number_of_restrictions)%term(i) = &
                                              nkappa(term(i)%n,term(i)%kappa)
      end do
      restriction(number_of_restrictions)%relation   = relation
      restriction(number_of_restrictions)%limitation = limitation
      fail = .false.  
      !
   end subroutine toolbox_excl_interpret_restr
   !
   !
   subroutine toolbox_excl_pair_correlation()
   !--------------------------------------------------------------------
   ! Determines and excludes those CSF from a csl list which are not
   ! directly 'connected' to the CSF of a given reference list, i.e. which
   ! does not have non-zero Hamiltonian matrix elements. The selection can
   ! be made by means of the Dirac-Coloumb or Dirac-Coulomb-Breit 
   ! Hamiltonian.
   !
   ! Calls: load_csl_from_file(), file_open().
   !--------------------------------------------------------------------
      !
      integer :: icsf, i, j, ja, jb, jc, jd, la, lb, lc, ld, nu, nx, &
                 no_T_coeff, no_V_coeff
      logical :: H_DC, H_DCB, is_coupled, yes
      real(kind=dp) :: xc
      !
    1 print *, "Generate pair-correlation list based on the"//&
               " Dirac-Coulomb Hamiltonian ?"
      yes = get_yes_stream()
      if (yes) then
         H_DC = .true.;    H_DCB = .false.
      else
         print *, "Generate pair-correlation list based on the"//&
                  " Dirac-Coulomb-Breit Hamiltonian ?"
         yes = get_yes_stream()
         if (yes) then  
            H_DCB = .true.;    H_DC = .false.
         else
            print *, "Use either one H_DC or H_DCB; redo ..."
            goto 1
         end if      
      end if
      !
      ! Create a 'new' configuration scheme which includes all reference
      ! CSF and one CSF from the current list;
      ! here, we assume that the number of orbitals in the reference list
      ! is smaller or equivalent to the number of the exclude list and that
      ! the sequence of orbitals does agree for these common orbitals.
      ref_csf_set%number_of_electrons = exclude_csf_set%number_of_electrons
      ref_csf_set%nocsf               = reference_csf_set%nocsf + 1
      ref_csf_set%nwshells            = exclude_csf_set%nwshells
      ref_csf_set%nwcore              = exclude_csf_set%nwcore
      !
      if (reference_csf_set%nwshells > ref_csf_set%nwshells) then
         print *, "The reference list should have less or an equivalent "//&
                  "number of orbitals in the same order "
         print *, "like the list which is to reduced."
         stop 
      else if (reference_csf_set%number_of_electrons /=  &
                     ref_csf_set%number_of_electrons) then
         print *, "Non-consistent number of electrons in the reference " //&
                  "list the list which is to reduced."
         stop 
      end if
      !
      ! Allocate memory for the subshells and the CSF of ref_csf_set
      ! and take over the subshell definition
      allocate( ref_csf_set%subshell(1:ref_csf_set%nwshells),  &
                ref_csf_set%csf(1:ref_csf_set%nocsf) )
      do  i = 1,ref_csf_set%nwshells
         ref_csf_set%subshell(i) = nkappa(exclude_csf_set%subshell(i)%n,   &
                                          exclude_csf_set%subshell(i)%kappa)
      end do
      !
      do  i = 1,reference_csf_set%nwshells
         if (reference_csf_set%subshell(i)%n     /=  &
                   ref_csf_set%subshell(i)%n    .or. &
	     reference_csf_set%subshell(i)%kappa /=  &
                   ref_csf_set%subshell(i)%kappa) then
            print *, "The orbitals of the reference list should have the "//&
                     "same order like the list which is to reduced"
            print *, "for all orbitals of the reference list."
            stop 
         end if
      end do
      !
      nx = reference_csf_set%nwshells
      do  i = 1,ref_csf_set%nocsf-1
         ref_csf_set%csf(i)%totalJ = reference_csf_set%csf(i)%totalJ
         ref_csf_set%csf(i)%parity = reference_csf_set%csf(i)%parity
         allocate( ref_csf_set%csf(i)%occupation(1:ref_csf_set%nwshells), &
                   ref_csf_set%csf(i)%seniority(1:ref_csf_set%nwshells),  &
                   ref_csf_set%csf(i)%subshellJ(1:ref_csf_set%nwshells),  &
                   ref_csf_set%csf(i)%subshellX(1:ref_csf_set%nwshells) )
         ref_csf_set%csf(i)%occupation(:) = 0
         ref_csf_set%csf(i)%seniority(:)  = 0
         ref_csf_set%csf(i)%subshellJ(:)  = 0
         ref_csf_set%csf(i)%subshellX(:)  = 0
         ref_csf_set%csf(i)%occupation(1:nx) =  &
                                   reference_csf_set%csf(i)%occupation(1:nx)
         ref_csf_set%csf(i)%seniority(1:nx)  =  &
                                   reference_csf_set%csf(i)%seniority(1:nx)
         ref_csf_set%csf(i)%subshellJ(1:nx)  =  &
                                   reference_csf_set%csf(i)%subshellJ(1:nx)
         ref_csf_set%csf(i)%subshellX(1:nx)  =  &
                                   reference_csf_set%csf(i)%subshellX(1:nx)
         ref_csf_set%csf(i)%subshellX(nx+1:)  =  &
                                   reference_csf_set%csf(i)%subshellX(nx)
      end do
      !
      nx = reference_csf_set%nocsf
      do  icsf = 1,exclude_csf_set%nocsf
         is_coupled = .false.
         ref_csf_set%csf(nx+1)%totalJ = exclude_csf_set%csf(icsf)%totalJ
         ref_csf_set%csf(nx+1)%parity = exclude_csf_set%csf(icsf)%parity
         ref_csf_set%csf(nx+1)%occupation=> exclude_csf_set%csf(icsf)%occupation
         ref_csf_set%csf(nx+1)%seniority => exclude_csf_set%csf(icsf)%seniority
         ref_csf_set%csf(nx+1)%subshellJ => exclude_csf_set%csf(icsf)%subshellJ
         ref_csf_set%csf(nx+1)%subshellX => exclude_csf_set%csf(icsf)%subshellX
         !!x call print_configuration_scheme(6,ref_csf_set)
         j_loop:  do  j = 1,reference_csf_set%nocsf
            call anco_calculate_csf_pair(ref_csf_set,j,nx+1,  &
                                         no_T_coeff,no_V_coeff)
            if (no_T_coeff > 0  .or.  no_V_coeff > 0) then
               if (H_DCB) then
                  is_coupled = .true.
                  exit j_loop
               else if (H_DC) then
	          do  i = 1,no_V_coeff
	             nu = anco_V_list(i)%nu
                     ja = angular_momentum_j(anco_V_list(i)%a%kappa)
                     jb = angular_momentum_j(anco_V_list(i)%b%kappa)
                     jc = angular_momentum_j(anco_V_list(i)%c%kappa)
                     jd = angular_momentum_j(anco_V_list(i)%d%kappa) 
                     !
                     la = angular_momentum_l(anco_V_list(i)%a%kappa)
                     lb = angular_momentum_l(anco_V_list(i)%b%kappa) 
                     lc = angular_momentum_l(anco_V_list(i)%c%kappa)
                     ld = angular_momentum_l(anco_V_list(i)%d%kappa) 
                     if (triangle(ja+1,jc+1,nu+nu+1) *                       &
	                triangle(jb+1,jd+1,nu+nu+1)== 0  .or.                &
                        mod(la+lc+nu,2) == 1 .or. mod(lb+ld+nu,2) == 1) then
                        cycle
                     end if               
                     xc= CL_reduced_me                                       &
		         (anco_V_list(i)%a%kappa,nu,anco_V_list(i)%c%kappa)* &
                         CL_reduced_me                                       &
		         (anco_V_list(i)%b%kappa,nu,anco_V_list(i)%d%kappa)
                     if (mod(nu,2) == 1) then
                        xc = - xc
                     end if
                     !
                     if (abs(xc) .gt. eps10) then
                        is_coupled = .true.
                        exit j_loop
                     end if
                  end do
               end if
            end if
         end do j_loop
         exclude_prelim(icsf)   = is_coupled
         exclude_accepted(icsf) = is_coupled
      end do
      !
   end subroutine toolbox_excl_pair_correlation
   !
   !
   subroutine toolbox_exclude_write_csl(csl_old_file)
   !--------------------------------------------------------------------
   ! Writes the 'accepted' CSF from csl_old_file to a new .csl file.
   !--------------------------------------------------------------------
      !
      character(len=256), intent(in) :: csl_old_file
      !
      integer :: csf_number, csf_counter, i, ierr
      character(len=256) :: csl_new_file, record1, record2, record3 
      !
      ! Reopen the 'old' .csl file; no problems should arise for the given
      ! name since this file has been opened successful previously
      call file_open(21,csl_old_file,"formatted  ","old",ierr)
      if (rabs_use_stop   .and.  ierr == 1) then
         stop "toolbox_exclude_write_csl(): program stop A."
      end if
      !
      !
    1 print *, "Enter the name of the configuration symmetry list file"//&
               " to be created:"
      read (*,"(a)") csl_new_file
      if (len(trim(csl_new_file)) == 0) goto 1
      !
      call file_open(22,csl_new_file,"formatted  ","new",ierr)
      if (ierr == 1) goto 1
      !
      ! Copy the head of the .csl file
      read (21,"(a)") record1
      write(22,"(a)") trim(record1)
      read (21,"(a)") record1
      write(22,"(a)") trim(record1)
      read (21,"(a)") record1
      write(22,"(a)") trim(record1)
      read (21,"(a)") record1
      write(22,"(a)") trim(record1)
      read (21,"(a)") record1
      write(22,"(a)") trim(record1)
      !
      ! Now read 'one by one' CSF and write it out if it belongs to the
      ! accepted list
      csf_number = 0;   csf_counter = 0
      do  i = 1,exclude_csf_set%nocsf
         read(21,"(a)") record1
         read(21,"(a)") record2
         read(21,"(a)") record3
         if (exclude_accepted(i)) then
            csf_counter = csf_counter + 1
            write(22,"(a)") trim(record1)
            write(22,"(a)") trim(record2)
            write(22,"(a)") trim(record3)
         end if
      end do
      !
      print *, "csf_number = ",csf_number
      !
      print *, "Write out of the new .csl file complete;"
      print *, " there are ",csf_counter," relativistic CSF in this list."
      !
   end subroutine toolbox_exclude_write_csl
   !
   !
   subroutine toolbox_format_mix()
   !--------------------------------------------------------------------
   ! Converts 'unformatted' GRASP92 .mix files into a 'formatted' one.
   !
   ! Calls: get_yes_stream(), file_open().
   !--------------------------------------------------------------------
      !
      logical            :: yes
      character(len=6)   :: g92mix
      character(len=256) :: file_name
      integer            :: i, ierr, ios, j, mx, npr, nak
      real(kind=dp)      :: e, pz
      integer, dimension(:), allocatable       :: iatjpo, iaspar
      real(kind=dp), dimension(:), allocatable :: eval, grid
      !
      ! Open the unformatted file and read in all data
      !
    2 print *, "Enter a file name for the GRASP92 .mix coefficient file:"
      read (*,"(a)") file_name
      call file_open(25,file_name,"unformatted","old",ierr)
      if (ierr /= 0  .or.  len_trim(file_name) == 0) goto 2
      !
      ! Check the header of the file; if not as expected, try again
      read (25,iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         print *, "ios, g92mix = ",ios, g92mix
         print *, "Not a GRASP92 Mixing Coefficients File;"
         close (25)
         goto 2
      end if
      !
      read (25) asf_set%csf_set%number_of_electrons, asf_set%csf_set%nocsf, &
                asf_set%csf_set%nwshells
      !
      ! Allocate memory for the set of ASF and load the data from the .mix file
      print *, "Loading Mixing Coefficients File ..."
      read (25) asf_set%noasf
      allocate( asf_set%asf(1:asf_set%noasf) )
      do  i = 1,asf_set%noasf
         allocate( asf_set%asf(i)%eigenvector(1:asf_set%csf_set%nocsf) )
      end do
      allocate( iatjpo(1:asf_set%noasf), iaspar(1:asf_set%noasf), &
                eval(1:asf_set%noasf) )
      !
      read (25) (asf_set%asf(i)%level_No,i = 1,asf_set%noasf)
      read (25) (iatjpo(i),iaspar(i),i = 1,asf_set%noasf)
      do  i = 1,asf_set%noasf
         asf_set%asf(i)%totalJ = iatjpo(i) - 1
         if (iaspar(i) == 1) then
            asf_set%asf(i)%parity = "+"
         else if (iaspar(i) == -1) then
            asf_set%asf(i)%parity = "-"
         else
            stop "toolbox_format(): program stop A."
         end if
      end do
      !
      read (25) asf_set%average_energy,(eval(i),i = 1,asf_set%noasf)
      do  i = 1,asf_set%noasf
         asf_set%asf(i)%energy = eval(i)
      end do
      !x read (25) ((evec(i+(j-1)*ncf),i = 1,ncf),j = 1,nvec)
      read (25) ((asf_set%asf(i)%eigenvector(j),j=1,asf_set%csf_set%nocsf),&
                 i=1,asf_set%noasf)
      !
      deallocate( iatjpo, iaspar, eval )
      print *, ' ... load complete;'
      close (25)
      !
      ! Open a 'formatted' file and write out all data
      !
    3 print *, "Enter a file name for the 'formatted' .mix Coefficient File:"
      read (*,"(a)") file_name
      call file_open(25,file_name,"formatted  ","new",ierr)
      if (ierr /= 0  .or.  len_trim(file_name) == 0) goto 3
      !
      ! Write out a file header
      write (25,"(a)") g92mix//" (formatted file version)."
      !
      print *, "Write formatted Mixing Coefficients File ..."
      write (25,"(3i6)") asf_set%csf_set%number_of_electrons, &
                         asf_set%csf_set%nocsf, asf_set%csf_set%nwshells
      write (25,"(i6)")  asf_set%noasf
      !
      write (25,"(100i6)") (asf_set%asf(i)%level_No,i = 1,asf_set%noasf)
      write (25,"(100(i5,a1))") (asf_set%asf(i)%totalJ,asf_set%asf(i)%parity, &
                          i = 1,asf_set%noasf)
      !
      write (25,"(101e16.9)") asf_set%average_energy, &
                             (asf_set%asf(i)%energy,i = 1,asf_set%noasf)
      do  j = 1,asf_set%csf_set%nocsf
         write (25,"(100e16.9)") (asf_set%asf(i)%eigenvector(j), &
                                  i=1,asf_set%noasf)
      end do
      !
      print *, ' ... write out complete;'
      close (25)
      !
      do  i = 1,asf_set%noasf
         deallocate( asf_set%asf(i)%eigenvector )
      end do
      deallocate( asf_set%asf )
      close(25);   close(26)
      !
   end subroutine toolbox_format_mix
   !
   !
   subroutine toolbox_format_out()
   !--------------------------------------------------------------------
   ! Converts 'unformatted' GRASP92 .out file into a 'formatted' file.
   !
   ! Calls: get_yes_stream(), file_open().
   !--------------------------------------------------------------------
      !
      logical            :: yes
      character(len=6)   :: g92rwf
      character(len=256) :: file_name
      integer            :: i, ierr, ios, j, mx, npr, nak
      real(kind=dp)      :: e, pz
      integer, dimension(:), allocatable       :: iatjpo, iaspar
      real(kind=dp), dimension(:), allocatable :: eval, grid, p_component, &
                                                  q_component
      !
      ! Attempt to open the unformatted file and read in all data
      !
    4 print *, "Enter the name of the GRASP92 Radial Wavefunction File:"
      read (*,"(a)") file_name
      call file_open(25,file_name,"unformatted","old",ierr)
      if (ierr /= 0  .or.  len_trim(file_name) == 0) goto 4
      !
      ! Check the header of the file; if not as expected, try again
      read (25) g92rwf
      if (g92rwf /= 'G92RWF') then
         close (25)
         print *, "g92rwf = ",g92rwf
         print *, "Not a GRASP92 radial wavefunction file;"
         goto 4
      endif
      !
      ! At this point the  .rwf  file has been opened and checked;
      ! Open a formatted file
    5 print *, "Enter a file name for the 'formatted' .rwf Radial "// &
               "Wavefunction File:"
      read (*,"(a)") file_name
      call file_open(26,file_name,"formatted  ","new",ierr)
      if (ierr /= 0  .or.  len_trim(file_name) == 0) goto 5
      !
      ! Write out a file header
      write (26,"(a)") g92rwf//" (formatted file version)."
      !
      ! Read (unformatted) and write (formatted) the parameters of each orbital
    6 read (25,end =7) npr, nak, e, mx
      write (26,*) npr, nak, e, mx
      allocate( p_component(1:mx), q_component(1:mx), grid(1:mx) )
      read  (25)   pz, (p_component(i),i = 1,mx),(q_component(i),i = 1,mx)
      read  (25)   (grid(i),i = 1,mx)
      write (26,*) pz
      do  i = 1,mx
         write (26,*) grid(i),p_component(i),q_component(i)
      end do
      deallocate( p_component, q_component, grid )
      goto 6
    7 print *, "Conversion of Radial Wavefunction File complete."
      !
      close(25);   close(26)
      !
   end subroutine toolbox_format_out
   !
   !
   subroutine toolbox_gather_energies(no_files)
   !--------------------------------------------------------------------
   ! Collects all level information from one or several .mix files.
   !
   ! Calls: file_open(). 
   !--------------------------------------------------------------------
      !
      integer, intent(out) :: no_files
      !
      character(len=6)   :: g92mix
      character(len=512) :: record
      integer            :: i, ierr, ios, j, k, noasf
      real(kind=dp)      :: wa
      character(len=256), dimension(10)        :: level_mix_file
      integer, dimension(:), allocatable       :: iatjpo, iaspar
      real(kind=dp), dimension(:), allocatable :: eval
      !
    1 print *, "Enter one (or several)  .mix  mixing coefficient file(s):"
      read(*,"(a)")  record
      !
      ! Determine file names
      no_files = 0
    2 record = adjustl(record)
      !!x print *, "record = ",record
      i = scan(record," ")
      !!x print *, "i = ",i
      if (i /= 1) then
         no_files = no_files + 1
         level_mix_file(no_files) = record(1:i-1)
         record(:) = record(i:)
         goto 2
      end if
      !
      ! Try to open these files and to gather all necessary information
      do  i = 1,no_files
         call file_open(25,level_mix_file(i),"unformatted","old",ierr)
         if (ierr /= 0  .or.  len_trim(level_mix_file(i)) == 0) goto 4
         !
         ! Check the header of the file; if not as expected for unformatted
         ! files, check formatted form
         read (25,iostat=ios)  g92mix
         if (ios /= 0   .or.   g92mix /= "G92MIX") then
            close (25)
            goto 3
         end if
         !
         ! Collect the energies from the file
         read (25) asf_mix(i)%csf_set%number_of_electrons,                 &
                   asf_mix(i)%csf_set%nocsf, asf_mix(i)%csf_set%nwshells
         read (25) asf_mix(i)%noasf
         allocate( asf_mix(i)%asf(1:asf_mix(i)%noasf) )
         allocate( iatjpo(1:asf_mix(i)%noasf), iaspar(1:asf_mix(i)%noasf), &
                   eval(1:asf_mix(i)%noasf) )
         !
         !!x print *, "asf_mix(i)%noasf = ",asf_mix(i)%noasf
         read (25) (asf_mix(i)%asf(j)%level_No,j = 1,asf_mix(i)%noasf)
         read (25) (iatjpo(j),iaspar(j),j = 1,asf_mix(i)%noasf)
         do  j = 1,asf_mix(i)%noasf
            asf_mix(i)%asf(j)%totalJ = iatjpo(j) - 1
            if (iaspar(j) == 1) then
               asf_mix(i)%asf(j)%parity = "+"
            else if (iaspar(j) == -1) then
               asf_mix(i)%asf(j)%parity = "-"
            else
               stop "toolbox_control_levels_simple(): program stop A."
            end if
         end do
         !
         read (25) asf_mix(i)%average_energy,(eval(j),j = 1,asf_mix(i)%noasf)
         do  j = 1,asf_mix(i)%noasf
            asf_mix(i)%asf(j)%energy = asf_mix(i)%average_energy + eval(j)
         end do
         deallocate( iatjpo, iaspar, eval )
         !
         ! Take the next .mix mixing file
         cycle
         !
         ! Try formatted file format
       3 call file_open(25,level_mix_file(i),"formatted  ","old",ierr)
         if (ierr /= 0  .or.  len_trim(level_mix_file(i)) == 0) goto 4
         !
         ! Check the header of the file; if not as expected for formatted
         ! files, try again
         read (25,"(a)",iostat=ios)  g92mix
         if (ios /= 0   .or.   g92mix /= "G92MIX") then
            print *, "ios, g92mix = ",ios, g92mix
            print *, "Not a GRASP92 Mixing Coefficients File;"
            close (25)
            goto 5
         end if
         !
         ! Collect the energies from the file
         read (25,"(3i6)") asf_mix(i)%csf_set%number_of_electrons,             &
                           asf_mix(i)%csf_set%nocsf, asf_mix(i)%csf_set%nwshells
         read (25,"(i6)") asf_mix(i)%noasf
         allocate( asf_mix(i)%asf(1:asf_mix(i)%noasf) )
         allocate( iatjpo(1:asf_mix(i)%noasf), iaspar(1:asf_mix(i)%noasf), &
                   eval(1:asf_mix(i)%noasf) )
         !
         !!x print *, "asf_mix(i)%noasf = ",asf_mix(i)%noasf
         noasf = min(60,asf_mix(i)%noasf)
         read (25,"(60i6)")      (asf_mix(i)%asf(j)%level_No,j = 1,noasf)
         read (25,"(60(i5,a1))") (asf_mix(i)%asf(j)%totalJ,               &
                                  asf_mix(i)%asf(j)%parity,j = 1,noasf)
         read (25,"(61e26.19)")  asf_mix(i)%average_energy,(eval(j),j = 1,noasf)
         do  j = 1,noasf
            asf_mix(i)%asf(j)%energy = asf_mix(i)%average_energy + eval(j)
         end do
	 !
         do  j = 1,asf_mix(i)%csf_set%nocsf
            read (25,"(60e16.9)") wa
         end do
         !
         if (asf_mix(i)%noasf <= 60) goto 4
         noasf = min(120,asf_mix(i)%noasf)
         read (25,"(60i6)")      (asf_mix(i)%asf(j)%level_No,j = 61,noasf)
         read (25,"(60(i5,a1))") (asf_mix(i)%asf(j)%totalJ,              &
                                  asf_mix(i)%asf(j)%parity,j = 61,noasf)
         read (25,"(61e26.19)")   (eval(j),j = 61,noasf)
         do  j = 61,noasf
            asf_mix(i)%asf(j)%energy = asf_mix(i)%average_energy + eval(j)
         end do
	 !
         do  j = 1,asf_mix(i)%csf_set%nocsf
            read (25,"(60e16.9)") wa
         end do
         !
         if (asf_mix(i)%noasf <= 120) goto 4
         noasf = min(180,asf_mix(i)%noasf)
         read (25,"(60i6)")      (asf_mix(i)%asf(j)%level_No,j = 121,noasf)
         read (25,"(60(i5,a1))") (asf_mix(i)%asf(j)%totalJ,               &
                                  asf_mix(i)%asf(j)%parity,j = 121,noasf)
         read (25,"(61e26.19)")   (eval(j),j = 121,noasf)
         do  j = 121,noasf
            asf_mix(i)%asf(j)%energy = asf_mix(i)%average_energy + eval(j)
         end do
	 !
         do  j = 1,asf_mix(i)%csf_set%nocsf
            read (25,"(60e16.9)") wa
         end do
         !
         if (asf_mix(i)%noasf <= 180) goto 4
         noasf = min(240,asf_mix(i)%noasf)
         read (25,"(60i6)")      (asf_mix(i)%asf(j)%level_No,j = 181,noasf)
         read (25,"(60(i5,a1))") (asf_mix(i)%asf(j)%totalJ,               &
                                  asf_mix(i)%asf(j)%parity,j = 181,noasf)
         read (25,"(61e26.9)")   (eval(j),j = 181,noasf)
         do  j = 181,noasf
            asf_mix(i)%asf(j)%energy = asf_mix(i)%average_energy + eval(j)
         end do
	 !
         do  j = 1,asf_mix(i)%csf_set%nocsf
            read (25,"(60e16.9)") wa
         end do
         !
         if (asf_mix(i)%noasf <= 240) goto 4
         noasf = min(300,asf_mix(i)%noasf)
         read (25,"(60i6)")      (asf_mix(i)%asf(j)%level_No,j = 241,noasf)
         read (25,"(60(i5,a1))") (asf_mix(i)%asf(j)%totalJ,               &
                                  asf_mix(i)%asf(j)%parity,j = 241,noasf)
         read (25,"(61e26.19)")   (eval(j),j = 241,noasf)
         do  j = 241,noasf
            asf_mix(i)%asf(j)%energy = asf_mix(i)%average_energy + eval(j)
         end do
	 !
         do  j = 1,asf_mix(i)%csf_set%nocsf
            read (25,"(60e16.9)") wa
         end do
         !
         if (asf_mix(i)%noasf <= 300) goto 4
         noasf = min(360,asf_mix(i)%noasf)
         read (25,"(60i6)")      (asf_mix(i)%asf(j)%level_No,j = 301,noasf)
         read (25,"(60(i5,a1))") (asf_mix(i)%asf(j)%totalJ,               &
                                  asf_mix(i)%asf(j)%parity,j = 301,noasf)
         read (25,"(61e26.19)")   (eval(j),j = 301,noasf)
         do  j = 301,noasf
            asf_mix(i)%asf(j)%energy = asf_mix(i)%average_energy + eval(j)
         end do
	 !
         do  j = 1,asf_mix(i)%csf_set%nocsf
            read (25,"(60e16.9)") wa
         end do
         !
         if (asf_set%noasf <= 360) goto 4
         noasf = min(420,asf_mix(i)%noasf)
         read (25,"(60i6)")      (asf_mix(i)%asf(j)%level_No,j = 361,noasf)
         read (25,"(60(i5,a1))") (asf_mix(i)%asf(j)%totalJ,               &
                                  asf_mix(i)%asf(j)%parity,j = 361,noasf)
         read (25,"(61e26.19)")   (eval(j),j = 361,noasf)
         do  j = 361,noasf
            asf_mix(i)%asf(j)%energy = asf_mix(i)%average_energy + eval(j)
         end do
	 !
         do  j = 1,asf_mix(i)%csf_set%nocsf
            read (25,"(60e16.9)") wa
         end do
         !
         if (asf_set%noasf <= 420) goto 4
         noasf = min(480,asf_mix(i)%noasf)
         read (25,"(60i6)")      (asf_mix(i)%asf(j)%level_No,j = 421,noasf)
         read (25,"(60(i5,a1))") (asf_mix(i)%asf(j)%totalJ,               &
                                  asf_mix(i)%asf(j)%parity,j = 421,noasf)
         read (25,"(61e26.19)")   (eval(j),j = 421,noasf)
         do  j = 421,noasf
            asf_mix(i)%asf(j)%energy = asf_mix(i)%average_energy + eval(j)
         end do
	 !
         do  j = 1,asf_mix(i)%csf_set%nocsf
            read (25,"(60e16.9)") wa
         end do
         !
         if (asf_set%noasf <= 480) goto 4
         noasf = min(540,asf_mix(i)%noasf)
         read (25,"(60i6)")      (asf_mix(i)%asf(j)%level_No,j = 481,noasf)
         read (25,"(60(i5,a1))") (asf_mix(i)%asf(j)%totalJ,               &
                                  asf_mix(i)%asf(j)%parity,j = 481,noasf)
         read (25,"(61e26.19)")   (eval(j),j = 481,noasf)
         do  j = 481,noasf
            asf_mix(i)%asf(j)%energy = asf_mix(i)%average_energy + eval(j)
         end do
	 !
         do  j = 1,asf_mix(i)%csf_set%nocsf
            read (25,"(60e16.9)") wa
         end do
         !
         if (asf_set%noasf <= 540) goto 4
         stop "toolbox_control_levels_simple(): program stop B."
         !
       4 deallocate( iatjpo, iaspar, eval )
         !
         cycle
         !
         ! Treat the case than one of the file names is inappropriate
       5 do  j = 1,i-1
            do  k = 1,asf_mix(j)%noasf
               deallocate( asf_mix(j)%asf(k)%eigenvector )
            end do
            deallocate( asf_mix(j)%asf )
         end do
         goto 1
      end do
      !
   end subroutine toolbox_gather_energies
   !
end module rabs_toolbox_df
