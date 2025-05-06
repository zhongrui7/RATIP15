program xtoolbox
!
!-----------------------------------------------------------------------
! This 'toolbox' program supports a number of task which frequently
! occur in atomic structure calculations. I basically comprises a set
! of small components within a common environment to keep the RATIP program
! clear and easy to use. The necesssary input and output files depend on 
! the task to be performed. Each task has its own independent dialog
! and terminates the program after its executation. Thus, to carry out two
! such tasks, the program has to be started twice. This is in line with
! most applications.
! The task which are presently implemented in xutility is seen from the
! menu below; after this menu, the execution continues, as usual, by an
! task-dependent, interactive dialog.
!-----------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_file_handling
   use rabs_grasp2k
   use rabs_toolbox_aux
   use rabs_toolbox_mess
   use rabs_toolbox_ac
   use rabs_toolbox_df
   use rabs_toolbox_gk
   use rabs_toolbox_ln
   use rabs_toolbox_or
   use rabs_toolbox_sz
   use rabs_toolbox_tpdi
   !
   implicit none
   character(len=1) :: sel1
   character(len=2) :: sel
   !
   print *, "XTOOLBOX: Supports a number of frequently used tasks"//&
            " (Fortran 95 version);"
   print *, " (C) Copyright by S Fritzsche, Kassel (2005)."
   print *, " "
   !
 1 print *, "Select from the menu:"
   print *, " "
   print *, "   A - Level energies, notations and weights."
   print *, "   C - Manipulation of GRASP2K .csl lists, .mix and .out files."
   print *, "   D - Format or unformat a GRASP2K file."
   print *, "   E - Properties of atomic orbitals."
   print *, " "
   print *, "   M - Miscellaneous."
   print *, "   N - Nuclear distributions and parameters."
   print *, " "
   print *, "   Q - Quit the program."
   !
   read *, sel1
   select case(sel1)
   case("A", "a");  goto 11
   case("C", "c");  goto 17
   case("D", "d");  goto 20
   case("E", "e");  goto 23
   case("M", "m");  goto 42
   case("N", "n");  goto 44
   case("Q", "q");  goto 61
   case default
      print *, "Unable to recognize the selection; redo ..."
      goto 1
   end select
   !
   !
   ! ===========================================
   ! A  -  Level energies, notations and weights
   ! ===========================================
   !
11 print *, "Select from the menu:"
   print *, " "
   print *, "   1 - Energy levels and level splittings from one or several .mix files."
   print *, "   2 - Energy levels and level splittings, extended form."
   print *, "   3 - Energy levels and level splittings from one  .mix file with LSJ notations."
   print *, "   4 - Display the major jj-coupled CSF and their weights to atomic levels."
   print *, "   5 - Display the major LS-coupled CSF and their weights to atomic levels."
   print *, " "
   print *, "   r - Return to the main menue."
   !
   read *, sel
   select case(sel)
   case("1 ", " 1")
      do i = 1,xmessage_na01;	print *, xmessage_a01(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 11
      call toolbox_display_levels_simple();   stop
   case("2 ", " 2")
      do i = 1,xmessage_na02;	print *, xmessage_a02(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 11
      call toolbox_display_levels();   stop
   case("3 ", " 3")
      do i = 1,xmessage_na03;	print *, xmessage_a03(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 11
      call toolbox_display_levels_lsj();   stop
   case("4 ", " 4")
      do i = 1,xmessage_na04;	print *, xmessage_a04(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 11
      call toolbox_display_weights();   stop
   case("5 ", " 5")
      do i = 1,xmessage_na05;	print *, xmessage_a05(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 11
      call toolbox_control_weights_lsj();   stop
   case("r ", " r", "R ");  goto 1
   case default
      print *, "Unable to recognize the selection; redo ..."
      goto 11
   end select
   !
   !
   ! ========================================
   ! C  -  Manipulation of GRASP2K .csl lists
   ! ========================================
   !
17 print *, "Select from the menu:"
   print *, " "
   print *, "   1 - Exclude a number of CSF from a GRASP2K .csl list."
   print *, "   2 - Split a GRASP2K .csl list into J^P level groups."
   print *, "   3 - Merge two .csl lists from GRASP2K with the same core."
   print *, "   4 - Merge two .csl lists from GRASP2K with different cores."
   print *, "   5 - Condense a .csl list from GRASP2K on a single weight criterium."
   print *, "   6 - Condense a .csl list from GRASP2K as above but for some given levels only."
   print *, "   7 - Reduce a .mix file from GRASP2K using level numbers."
   print *, "   8 - Reduce a .mix file from GRASP2K using symmetry properties (J and P)."
   print *, "   9 - Generate a pair-list from a given .csl list with regard to reference list."
   print *, " " 
   print *, "   r - Return to the main menue."
   !
   read *, sel
   select case(sel)
   case("1 ", " 1")
      do i = 1,xmessage_nc01;	print *, xmessage_c01(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 17
      call toolbox_exclude_csf();   stop
   case("2 ", " 2")
      do i = 1,xmessage_nc02;	print *, xmessage_c02(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 17
      call toolbox_split_csl_list();   stop
   case("3 ", " 3")
      do i = 1,xmessage_nc03;	print *, xmessage_c03(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 17
      call toolbox_merge_csl_lists();   stop
   case("4 ", " 4")
      do i = 1,xmessage_nc04;	print *, xmessage_c04(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 17
      call toolbox_merge_csl_diffcores();   stop
   case("5 ", " 5")
      do i = 1,xmessage_nc05;	print *, xmessage_c05(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 17
      call toolbox_condense_csl_list();   stop
   case("6 ", " 6")
      do i = 1,xmessage_nc06;	print *, xmessage_c06(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 17
      call toolbox_condense_on_levels();	stop
   case("7 ", " 7")
      do i = 1,xmessage_nc07;	print *, xmessage_c07(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 17
      call toolbox_reduce_mix_file(.true.);	stop
   case("8 ", " 8")
      do i = 1,xmessage_nc08;	print *, xmessage_c08(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 17
      call toolbox_reduce_mix_file(.false.);	stop
   case("9 ", " 9")
      do i = 1,xmessage_nc09;	print *, xmessage_c09(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 17
      call toolbox_excl_generate_pairlist();   stop
   case("r ", " r", "R ");  goto 1
   case default
      print *, "Unable to recognize the selection; redo ..."
      goto 17
   end select
   !
   !
   ! =======================================
   ! D  -  Format or unformat a GRASP2K file
   ! =======================================
   !
20 print *, "Select from the menu:"
   print *, " "
   print *, "   1 - Format  .mix  mixing coefficient files from GRASP2K."
   print *, "   2 - Format  .out  radial orbital files from GRASP2K."
   print *, "   3 - Unformat  .mix  mixing coefficient files from GRASP2K."
   print *, "   4 - Unformat  .out  radial orbital files from GRASP2K."
   print *, " "
   print *, "   r - Return to the main menue."
   !
   read *, sel
   select case(sel)
   case("1 ", " 1")
      do i = 1,xmessage_nd01;	print *, xmessage_d01(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 20
      call toolbox_format_mix();	goto 20
   case("2 ", " 2")
      do i = 1,xmessage_nd02;	print *, xmessage_d02(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 20
      call toolbox_format_out();	goto 20
   case("3 ", " 3")
      do i = 1,xmessage_nd03;	print *, xmessage_d03(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 20
      call toolbox_unformat_mix();   goto 20
   case("4 ", " 4")
      do i = 1,xmessage_nd04;	print *, xmessage_d04(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 20
      call toolbox_unformat_out();   goto 20
   case("r ", " r", "R ");  goto 1
   case default
      print *, "Unable to recognize the selection; redo ..."
      goto 20
   end select
   !
   !
   ! ===================================
   ! E  -  Properties of atomic orbitals
   ! ===================================
   !
23 print *, "Select from the menu:"
   print *, " "
   print *, "   1 - Radial properties of atomic orbitals."
   print *, "   2 - Overlaps between two not quite orthogonal orbital sets."
   print *, " "
   print *, "   r - Return to the main menue."
   !
   read *, sel
   select case(sel)
   case("1 ", " 1")
      do i = 1,xmessage_ne01;	print *, xmessage_e01(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 23
      call toolbox_radial_properties();   stop
   case("2 ", " 2")
      do i = 1,xmessage_ne02;	print *, xmessage_e02(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 23
      call toolbox_overlap();   stop
   case("r ", " r", "R ");  goto 1
   case default
      print *, "Unable to recognize the selection; redo ..."
      goto 23
   end select
   !
   !
   ! ==================
   ! M  -  Miscelaneous
   ! ==================
   !
42 print *, "Select from the menu:"
   print *, " "
   print *, "   1 - Effective charge for a given orbital."
   print *, "   2 - Effective radial charge or charge density of a selected ASF."
   print *, " "
   print *, "   4 - GRASP2K grid-parameter calculator."
   print *, " "
   print *, "   r - Return to the main menue."
   !
   read *, sel
   select case(sel)
   case("1 ", " 1")
      do i = 1,xmessage_nm01;	print *, xmessage_m01(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 42
      call toolbox_charge_density_orbital();   stop
   case("2 ", " 2")
      do i = 1,xmessage_nm02;	print *, xmessage_m02(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 42
      call toolbox_charge_density();   stop
   case("4 ", " 4")
      do i = 1,xmessage_nm04;	print *, xmessage_m04(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 42
      call toolbox_grid_calculator();   stop
   case("r ", " r", "R ");  goto 1
   case default
      print *, "Unable to recognize the selection; redo ..."
      goto 42  
   end select
   !
   !
   ! =======================================
   ! N  -  Nuclear parameters and potentials
   ! =======================================
   !
44 print *, "Select from the menu:"
   print *, " "
   print *, "   1 - Calculate the nuclear radius from the mass number."
   print *, "   2 - Calculate the Fermi distribution parameters c and N as well as the distribution and potential."
   print *, " "
   print *, "   r - Return to the main menue."
   !
   read *, sel
   select case(sel)
   case("1 ", " 1")
      do i = 1,xmessage_nn01;	print *, xmessage_n01(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 44
      call toolbox_nuclear_radius();   stop
   case("2 ", " 2")
      do i = 1,xmessage_nn02;	print *, xmessage_n02(i);   end do;
      print *, "Continue (y/n) ?";   if (.not.get_yes_stream()) goto 44
      call toolbox_nuclear_fermi();   stop
   case("r ", " r", "R ");  goto 1
   case default
      print *, "Unable to recognize the selection; redo ..."
   goto 44  
   end select
   !
   !
61 print *, " "
   print *, "XTOOLBOX complete ... ."
   !
end program xtoolbox

