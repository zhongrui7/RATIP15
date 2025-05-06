module rabs_print
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module prints some `frequent output' for the various components.
!-----------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_csl
   use rabs_functions_string 
   use rabs_grasp2k  
   !
   public  :: initialize_rwf_storage
                 ! Initializes the arrays of type(grasp2k_orbital) for the 
		 ! storage of the radial wave functions.
   public  :: print_summary
                 ! Prints the summary information for an initial and final
		 ! set of states.
   public  :: print_summary_b
                 ! Prints the summary information for a single set of bound
                 ! states.
   public  :: print_summary_c
                 ! Prints the summary information for a single set of bound
                 ! orbitals only.
   public  :: print_summary_imf
                 ! Prints the summary information for a given module
		 ! if one need to distinguish between initial, intermediate
		 ! and final states.
   !
   !
contains
   !
   !
   subroutine initialize_rwf_storage(asf_initial,wave_initial)
   !--------------------------------------------------------------------
   ! Initializes the arrays of type(grasp2k_orbital) for the storage of 
   ! the radial wave functions.
   !--------------------------------------------------------------------
      !
      type(asf_basis), intent(in)        :: asf_initial
      type(grasp2k_orbital), intent(out) :: wave_initial
      !
      integer :: i
      !
      ! Initialize storage for initial-state wave functions
      !
      wave_initial%number_of_rwf = asf_initial%csf_set%nwshells
      allocate( wave_initial%rwf(1:wave_initial%number_of_rwf) )
      do  i = 1,asf_initial%csf_set%nwshells
         wave_initial%rwf(i)%orbital = asf_initial%csf_set%subshell(i)
         wave_initial%rwf(i)%mtp     = 0
         wave_initial%rwf(i)%energy  = zero
         wave_initial%rwf(i)%gamma   = zero
         wave_initial%rwf(i)%pz      = zero
         wave_initial%rwf(i)%phase   = zero
      end do
      !
   end subroutine initialize_rwf_storage
   !
   !
   subroutine print_summary(progstring,stream,asf_initial,asf_final, &
                            wave_initial, wave_final)
   !--------------------------------------------------------------------
   ! Appends a summary of the input data from the initial and final state
   ! .mix files to the .sum file.
   !
   ! Calls: angular_momentum_string(), get_month(), orbital_symmetry().
   !--------------------------------------------------------------------
      !
      character(len=*), intent(in)      :: progstring
      integer, intent(in)               :: stream
      type(asf_basis), intent(in)       :: asf_initial, asf_final
      type(grasp2k_orbital), intent(in) :: wave_initial, wave_final
      !
      integer           :: i, totalJ, n_R
      real(kind=dp)     :: Rnuc      
      character(len=3)  :: month
      character(len=8)  :: cdate
      character(len=10) :: ctime
      !
      ! Get the date and time of day; make this information the header of 
      ! the  auger.sum  summary file
      call date_and_time(date=cdate,time=ctime)
      month = get_month(cdate(5:6))
      write(stream,*) progstring//                                             &
                      " run at "//ctime(1:2)//":"//ctime(3:4)//":"//           &
                      ctime(5:6)//" on "//month//                              &
                      " "//cdate(7:8)//" "//cdate(1:4)//"."
      !
      ! Append the input to the .sum file
      call save_input(" ",.true.,stream)
      !
      ! Write out the basic dimensions of the initial and final-state
      ! electron clouds
      write(stream,*)
      write(stream,9)  "There are ",asf_initial%csf_set%number_of_electrons,   &
                       " electrons in the cloud of the initial ion"
      write(stream,9)  "       in ",asf_initial%csf_set%nocsf," CSF based on ",&
                       asf_initial%csf_set%nwshells," subshells ..."
      !
      write(stream,*)
      write(stream,9)  "  ... and ",asf_final%csf_set%number_of_electrons,     &
                       " electrons in the cloud of the final ion"
      write(stream,9)  "       in ",asf_final%csf_set%nocsf," CSF based on ",  &
                       asf_final%csf_set%nwshells," subshells."
      !
      ! Write out the nuclear parameters and the speed of light
      write(stream,*)
      write(stream,8) "The atomic number is Z =",nuclear_charge,	       &
                      " distributed in a "//trim(nuclear_model)//              &
                      "-like nucleus."  
      write(stream,8) "Speed of light         =",c," atomic units."
      !
      select case(nuclear_model)
      case("fermi")
         write(24,*) " "
         write(24,*) " Fermi nucleus:"
         write(24,*) "    c =",fermi_c_parameter,"Bohr radii,"
         write(24,*) "    a =",fermi_a_parameter,"Bohr radii;"
         Rnuc = nuclear_extent()
         do  i = 1,n_grasp2k
            if (r_grasp2k(i) > Rnuc) exit
         end do
         n_R = i
         write(24,*) " there are ",n_R,"tabulation points inside the nucleus."
      end select
      !
      ! Write out the parameters of the radial grid
      write(stream,*)
      if (hp_grasp2k == zero) then
         write(stream,1) rnt_grasp2k,h_grasp2k,n_grasp2k
      else
         write(stream,2) rnt_grasp2k,h_grasp2k,hp_grasp2k,n_grasp2k
      endif
      write(stream,3) r_grasp2k(1),r_grasp2k(2),r_grasp2k(n_grasp2k)
      !
    1 format(" Radial grid: R(I) = RNT*(exp((I-1)*H)-1), I = 1, ..., N;", &
            /" --------------------------------------------------------", &
           //"  RNT  = ",1p,d19.12," Bohr radii;",                        &
            /"  H    = ",   d19.12," Bohr radii;" ,                       &
            /"  N    = ",1i4,";")
    2 format(" Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,",       &
             " I = 1, ..., N;",                                           &
            /" --------------------------------------------------",       &
             "---------------",                                           &
           //"  RNT  = ",1p,d19.12," Bohr radii;",                        &
            /"  H    = ",   d19.12," Bohr radii;",                        &
            /"  HP   = ",   d19.12," Bohr radii;",                        &
            /"  N    = ",1i4,";")
    3 format("  R(1) = ",1p,1d19.12," Bohr radii;",                       &
            /"  R(2) = ",   1d19.12," Bohr radii;",                       &
            /"  R(N) = ",   1d19.12," Bohr radii.")
      !
      ! Write out the orbital properties
      write(stream,*)
      write(stream,*) "Initial-state subshell radial wavefunction summary:"
      write(stream,*) "---------------------------------------------------"
      write(stream,*)
      write(stream,4)
      write(stream,*)
      do  i = 1,wave_initial%number_of_rwf
         write(stream,5) wave_initial%rwf(i)%orbital%n,                   &
                     orbital_symmetry(wave_initial%rwf(i)%orbital%kappa), &
                     wave_initial%rwf(i)%energy,wave_initial%rwf(i)%pz,   &
                     wave_initial%rwf(i)%gamma,                           &
                     wave_initial%rwf(i)%P(2),wave_initial%rwf(i)%Q(2),   &
                     wave_initial%rwf(i)%mtp
      end do
    4 format(" Subshell",11x,"E",20x,"p0",18x,                            &
             "gamma",19x,"P(2)",18x,"Q(2)",11x,"mtp")
    5 format(3x,1i2,1a2,1x,1p,5(3x,1e19.12),3x,1i5)
    8 format(1x,a,f10.5,a)
    9 format(1x,a,i7,a,i7,a)
      !
      write(stream,*)
      write(stream,*) "Final-state subshell radial wavefunction summary:"
      write(stream,*) "-------------------------------------------------"
      write(stream,*)
      write(stream,4)
      write(stream,*)
      do  i = 1,wave_final%number_of_rwf
         write(stream,5) wave_final%rwf(i)%orbital%n,                     &
                     orbital_symmetry(wave_final%rwf(i)%orbital%kappa),   &
                     wave_final%rwf(i)%energy,wave_final%rwf(i)%pz,       &
                     wave_final%rwf(i)%gamma,                             &
                     wave_final%rwf(i)%P(2),wave_final%rwf(i)%Q(2),       &
                     wave_final%rwf(i)%mtp
      end do
      !
      ! Write the list of eigenpair indices
      write(stream,*)
      write(stream,*) 'Initial-state ASF:'
      write(stream,*) '------------------'
      write(stream,*)
      write(stream,6)
      write(stream,*)
      do  i = 1,asf_initial%noasf
         totalJ = asf_initial%asf(i)%totalJ
         write(stream,7) asf_initial%asf(i)%level_No,                     &
                     angular_momentum_string(1*totalJ,4),                 &
                     asf_initial%asf(i)%parity,asf_initial%asf(i)%energy
      end do
    6 format(" Level_No",6x,"J^P",10x,"Energy (a.u.)")
    7 format(2x,i4,8x,a4,a1,6x,1e19.12)
      !
      write(stream,*)
      write(stream,*) 'Final-state ASF:'
      write(stream,*) '----------------'
      write(stream,*)
      write(stream,6)
      write(stream,*)
      do  i = 1,asf_final%noasf
         totalJ = asf_final%asf(i)%totalJ
         write(stream,7) asf_final%asf(i)%level_No,                       &
                     angular_momentum_string(1*totalJ,4),                 &
                     asf_final%asf(i)%parity,asf_final%asf(i)%energy
      end do
      write(stream,*)
      !
   end subroutine print_summary
   !
   !
   subroutine print_summary_b(progstring,stream,asf_initial,wave_initial)
   !--------------------------------------------------------------------
   ! Appends a summary of the input data from the initial and final state
   ! .mix files to the .sum file.
   !
   ! Calls: angular_momentum_string(), get_month(), orbital_symmetry().
   !--------------------------------------------------------------------
      !
      character(len=*), intent(in)                :: progstring
      integer, intent(in)                         :: stream
      type(asf_basis), intent(in)                 :: asf_initial
      type(grasp2k_orbital), intent(in), optional :: wave_initial
      !
      integer           :: i, totalJ, n_R
      real(kind=dp)     :: Rnuc      
      character(len=3)  :: month
      character(len=8)  :: cdate
      character(len=10) :: ctime
      !
      ! Get the date and time of day; make this information the header of 
      ! the  auger.sum  summary file
      call date_and_time(date=cdate,time=ctime)
      month = get_month(cdate(5:6))
      write(stream,*) progstring//                                             &
                      " run at "//ctime(1:2)//":"//ctime(3:4)//":"//           &
                      ctime(5:6)//" on "//month//                              &
                      " "//cdate(7:8)//" "//cdate(1:4)//"."
      !
      ! Append the input to the .sum file
      call save_input(" ",.true.,stream)
      !
      ! Write out the basic dimensions of the initial and final-state
      ! electron clouds
      write(stream,*)
      write(stream,9)  "There are ",asf_initial%csf_set%number_of_electrons,   &
                       " electrons in the cloud of the bound ion"
      write(stream,9)  "       in ",asf_initial%csf_set%nocsf," CSF based on ",&
                       asf_initial%csf_set%nwshells," subshells."
      !
      !
      if (present(wave_initial)) then
      !
      ! Write out the nuclear parameters and the speed of light
      write(stream,*)
      write(stream,8) "The atomic number is Z =",nuclear_charge,	       &
                      " distributed in a "//trim(nuclear_model)//              &
                      "-like nucleus."  
      write(stream,8) "Speed of light         =",c," atomic units."
      !
      select case(nuclear_model)
      case("fermi")
         write(24,*) " "
         write(24,*) " Fermi nucleus:"
         write(24,*) "    c =",fermi_c_parameter,"Bohr radii,"
         write(24,*) "    a =",fermi_a_parameter,"Bohr radii;"
         Rnuc = nuclear_extent()
         do  i = 1,n_grasp2k
            if (r_grasp2k(i) > Rnuc) exit
         end do
         n_R = i
         write(24,*) " there are ",n_R,"tabulation points inside the nucleus."
      end select
      !
      ! Write out the parameters of the radial grid
      write(stream,*)
      if (hp_grasp2k == zero) then
         write(stream,1) rnt_grasp2k,h_grasp2k,n_grasp2k
      else
         write(stream,2) rnt_grasp2k,h_grasp2k,hp_grasp2k,n_grasp2k
      endif
      write(stream,3) r_grasp2k(1),r_grasp2k(2),r_grasp2k(n_grasp2k)
      !
    1 format(" Radial grid: R(I) = RNT*(exp((I-1)*H)-1), I = 1, ..., N;", &
            /" --------------------------------------------------------", &
           //"  RNT  = ",1p,d19.12," Bohr radii;",                        &
            /"  H    = ",   d19.12," Bohr radii;" ,                       &
            /"  N    = ",1i4,";")
    2 format(" Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,",       &
             " I = 1, ..., N;",                                           &
            /" --------------------------------------------------",       &
             "---------------",                                           &
           //"  RNT  = ",1p,d19.12," Bohr radii;",                        &
            /"  H    = ",   d19.12," Bohr radii;",                        &
            /"  HP   = ",   d19.12," Bohr radii;",                        &
            /"  N    = ",1i4,";")
    3 format("  R(1) = ",1p,1d19.12," Bohr radii;",                       &
            /"  R(2) = ",   1d19.12," Bohr radii;",                       &
            /"  R(N) = ",   1d19.12," Bohr radii.")
      !
      ! Write out the orbital properties
      write(stream,*)
      write(stream,*) "Bound-state subshell radial wavefunction summary:"
      write(stream,*) "-------------------------------------------------"
      write(stream,*)
      write(stream,4)
      write(stream,*)
      do  i = 1,wave_initial%number_of_rwf
         write(stream,5) wave_initial%rwf(i)%orbital%n,                   &
                     orbital_symmetry(wave_initial%rwf(i)%orbital%kappa), &
                     wave_initial%rwf(i)%energy,wave_initial%rwf(i)%pz,   &
                     wave_initial%rwf(i)%gamma,                           &
                     wave_initial%rwf(i)%P(2),wave_initial%rwf(i)%Q(2),   &
                     wave_initial%rwf(i)%mtp
      end do
    4 format(" Subshell",11x,"E",20x,"p0",18x,                            &
             "gamma",19x,"P(2)",18x,"Q(2)",11x,"mtp")
    5 format(3x,1i2,1a2,1x,1p,5(3x,1e19.12),3x,1i5)
    8 format(1x,a,f10.5,a)
    9 format(1x,a,i7,a,i7,a)
      !
      ! Write the list of eigenpair indices
      write(stream,*)
      write(stream,*) 'Bound-state ASF:'
      write(stream,*) '----------------'
      write(stream,*)
      write(stream,6)
      write(stream,*)
      do  i = 1,asf_initial%noasf
         totalJ = asf_initial%asf(i)%totalJ
         print *, "abc"
         write(stream,7) asf_initial%asf(i)%level_No,                     &
                     trim(angular_momentum_string(1*totalJ,4)),           &
                     asf_initial%asf(i)%parity,asf_initial%asf(i)%energy
      end do
    6 format(" Level_No",6x,"J^P",10x,"Energy (a.u.)")
    7 format(2x,i4,8x,a4,a1,6x,1e19.12)
      !
      end if
      !
   end subroutine print_summary_b
   !
   !
   subroutine print_summary_c(progstring,stream,asf_initial,wave_initial)
   !--------------------------------------------------------------------
   ! Appends a summary of the input data from the initial and final state
   ! .mix files to the .sum file.
   !
   ! Calls: angular_momentum_string(), get_month(), orbital_symmetry().
   !--------------------------------------------------------------------
      !
      character(len=*), intent(in)                :: progstring
      integer, intent(in)                         :: stream
      type(asf_basis), intent(in)                 :: asf_initial
      type(grasp2k_orbital), intent(in), optional :: wave_initial
      !
      integer           :: i, totalJ, n_R
      real(kind=dp)     :: Rnuc      
      character(len=3)  :: month
      character(len=8)  :: cdate
      character(len=10) :: ctime
      !
      ! Get the date and time of day; make this information the header of 
      ! the  auger.sum  summary file
      call date_and_time(date=cdate,time=ctime)
      month = get_month(cdate(5:6))
      write(stream,*) progstring//                                             &
                      " run at "//ctime(1:2)//":"//ctime(3:4)//":"//           &
                      ctime(5:6)//" on "//month//                              &
                      " "//cdate(7:8)//" "//cdate(1:4)//"."
      !
      ! Append the input to the .sum file
      call save_input(" ",.true.,stream)
      !
      ! Write out the basic dimensions of the bound-state electron cloud
      write(stream,*)
      write(stream,9)  "There are ",asf_initial%csf_set%number_of_electrons,   &
                       " electrons in the cloud of the bound ion"
      write(stream,9)  "       in ",asf_initial%csf_set%nocsf," CSF based on ",&
                       asf_initial%csf_set%nwshells," subshells."
      !
      !
      if (present(wave_initial)) then
      !
      ! Write out the nuclear parameters and the speed of light
      write(stream,*)
      write(stream,8) "The atomic number is Z =",nuclear_charge,	       &
                      " distributed in a "//trim(nuclear_model)//              &
                      "-like nucleus."  
      write(stream,8) "Speed of light         =",c," atomic units."
      !
      select case(nuclear_model)
      case("fermi")
         write(24,*) " "
         write(24,*) " Fermi nucleus:"
         write(24,*) "    c =",fermi_c_parameter,"Bohr radii,"
         write(24,*) "    a =",fermi_a_parameter,"Bohr radii;"
         Rnuc = nuclear_extent()
         do  i = 1,n_grasp2k
            if (r_grasp2k(i) > Rnuc) exit
         end do
         n_R = i
         write(24,*) " there are ",n_R,"tabulation points inside the nucleus."
      end select
      !
      ! Write out the parameters of the radial grid
      write(stream,*)
      if (hp_grasp2k == zero) then
         write(stream,1) rnt_grasp2k,h_grasp2k,n_grasp2k
      else
         write(stream,2) rnt_grasp2k,h_grasp2k,hp_grasp2k,n_grasp2k
      endif
      write(stream,3) r_grasp2k(1),r_grasp2k(2),r_grasp2k(n_grasp2k)
      !
    1 format(" Radial grid: R(I) = RNT*(exp((I-1)*H)-1), I = 1, ..., N;", &
            /" --------------------------------------------------------", &
           //"  RNT  = ",1p,d19.12," Bohr radii;",                        &
            /"  H    = ",   d19.12," Bohr radii;" ,                       &
            /"  N    = ",1i4,";")
    2 format(" Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,",       &
             " I = 1, ..., N;",                                           &
            /" --------------------------------------------------",       &
             "---------------",                                           &
           //"  RNT  = ",1p,d19.12," Bohr radii;",                        &
            /"  H    = ",   d19.12," Bohr radii;",                        &
            /"  HP   = ",   d19.12," Bohr radii;",                        &
            /"  N    = ",1i4,";")
    3 format("  R(1) = ",1p,1d19.12," Bohr radii;",                       &
            /"  R(2) = ",   1d19.12," Bohr radii;",                       &
            /"  R(N) = ",   1d19.12," Bohr radii.")
      !
      ! Write out the orbital properties
      write(stream,*)
      write(stream,*) "Bound-state subshell radial wavefunction summary:"
      write(stream,*) "-------------------------------------------------"
      write(stream,*)
      write(stream,4)
      write(stream,*)
      do  i = 1,wave_initial%number_of_rwf
         write(stream,5) wave_initial%rwf(i)%orbital%n,                   &
                     orbital_symmetry(wave_initial%rwf(i)%orbital%kappa), &
                     wave_initial%rwf(i)%energy,wave_initial%rwf(i)%pz,   &
                     wave_initial%rwf(i)%gamma,                           &
                     wave_initial%rwf(i)%P(2),wave_initial%rwf(i)%Q(2),   &
                     wave_initial%rwf(i)%mtp
      end do
    4 format(" Subshell",11x,"E",20x,"p0",18x,                            &
             "gamma",19x,"P(2)",18x,"Q(2)",11x,"mtp")
    5 format(3x,1i2,1a2,1x,1p,5(3x,1e19.12),3x,1i5)
    8 format(1x,a,f10.5,a)
    9 format(1x,a,i7,a,i7,a)
      !
      end if
      !
   end subroutine print_summary_c
   !
   !
   subroutine print_summary_imf(progstring,stream,                   &
                                asf_initial,asf_intermed,asf_final,  &
                                wave_initial,wave_intermed,wave_final)
   !--------------------------------------------------------------------
   ! Appends a summary of the input data from the initial, intermdiate
   ! and final state .mix files to the .sum files.
   !
   ! Calls: angular_momentum_string(), get_month(), orbital_symmetry().
   !--------------------------------------------------------------------
      !
      character(len=*), intent(in)      :: progstring
      integer, intent(in)               :: stream
      type(asf_basis), intent(in)       :: asf_initial, asf_intermed, asf_final
      type(grasp2k_orbital), intent(in) :: wave_initial, wave_intermed, &
                                           wave_final
      !
      integer           :: i, totalJ
      character(len=3)  :: month
      character(len=8)  :: cdate
      character(len=10) :: ctime
      !
      ! Get the date and time of day; make this information the header of 
      ! the  dierec.sum  summary file
      call date_and_time(date=cdate,time=ctime)
      month = get_month(cdate(5:6))
      write(stream,*) progstring//                                        &
                      " run at "//ctime(1:2)//":"//ctime(3:4)//":"//      &
                      ctime(5:6)//" on "//month//                         &
                      " "//cdate(7:8)//" "//cdate(1:4)//"."
      !
      ! Append the input to the .sum file
      call save_input(" ",.true.,stream)
      !
      ! Write out the basic dimensions of the initial, intermediate and 
      ! final-state electron clouds
      write(stream,*)
      write(stream,9)  "There are ",asf_initial%csf_set%number_of_electrons,   &
                       " electrons in the cloud of the initial ion"
      write(stream,9)  "       in ",asf_initial%csf_set%nocsf," CSF based on ",&
                       asf_initial%csf_set%nwshells," subshells ..."
      !
      write(stream,*)
      write(stream,9)  "  ... and ",asf_final%csf_set%number_of_electrons,     &
                       " electrons in the cloud of the intermediate-state ion"
      write(stream,9)  "       in ",asf_final%csf_set%nocsf," CSF based on ",  &
                       asf_final%csf_set%nwshells," subshells."
      !
      write(stream,*)
      write(stream,9)  "  ... and ",asf_intermed%csf_set%number_of_electrons,  &
                       " electrons in the cloud of the final ion"
      write(stream,9)  "       in ",asf_final%csf_set%nocsf," CSF based on ",  &
                       asf_final%csf_set%nwshells," subshells."
      !
      ! Write out the nuclear parameters and the speed of light
      write(stream,*)
      write(stream,8) "The atomic number is Z =",nuclear_charge,	       &
                      " distributed in a "//trim(nuclear_model)//              &
                      "-like nucleus."  
      write(stream,8) "Speed of light         =",c," atomic units."
      !
      select case(nuclear_model)
      case("fermi")
         write(24,*) " "
         write(24,*) " Fermi nucleus:"
         write(24,*) "    c =",fermi_c_parameter,"Bohr radii,"
         write(24,*) "    a =",fermi_a_parameter,"Bohr radii;"
         Rnuc = nuclear_extent()
         do  i = 1,n_grasp2k
            if (r_grasp2k(i) > Rnuc) exit
         end do
         n_R = i
         write(24,*) " there are ",n_R,"tabulation points inside the nucleus."
      end select
      !
      ! Write out the parameters of the radial grid
      write(stream,*)
      if (hp_grasp2k == zero) then
         write(stream,1) rnt_grasp2k,h_grasp2k,n_grasp2k
      else
         write(stream,2) rnt_grasp2k,h_grasp2k,hp_grasp2k,n_grasp2k
      endif
      write(stream,3) r_grasp2k(1),r_grasp2k(2),r_grasp2k(n_grasp2k)
      !
    1 format(" Radial grid: R(I) = RNT*(exp((I-1)*H)-1), I = 1, ..., N;", &
            /" --------------------------------------------------------", &
           //"  RNT  = ",1p,d19.12," Bohr radii;",                        &
            /"  H    = ",   d19.12," Bohr radii;" ,                       &
            /"  N    = ",1i4,";")
    2 format(" Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,",       &
             " I = 1, ..., N;",                                           &
            /" --------------------------------------------------",       &
             "---------------",                                           &
           //"  RNT  = ",1p,d19.12," Bohr radii;",                        &
            /"  H    = ",   d19.12," Bohr radii;",                        &
            /"  HP   = ",   d19.12," Bohr radii;",                        &
            /"  N    = ",1i4,";")
    3 format("  R(1) = ",1p,1d19.12," Bohr radii;",                       &
            /"  R(2) = ",   1d19.12," Bohr radii;",                       &
            /"  R(N) = ",   1d19.12," Bohr radii.")
      !
      !
      ! Write out the orbital properties
      write(stream,*)
      write(stream,*) "Initial-state subshell radial wavefunction summary:"
      write(stream,*) "---------------------------------------------------"
      write(stream,*)
      write(stream,4)
      write(stream,*)
      do  i = 1,wave_initial%number_of_rwf
         write(stream,5) wave_initial%rwf(i)%orbital%n,                   &
                     orbital_symmetry(wave_initial%rwf(i)%orbital%kappa), &
                     wave_initial%rwf(i)%energy,wave_initial%rwf(i)%pz,   &
                     wave_initial%rwf(i)%gamma,                           &
                     wave_initial%rwf(i)%P(2),wave_initial%rwf(i)%Q(2),   &
                     wave_initial%rwf(i)%mtp
      end do
    4 format(" Subshell",11x,"E",20x,"p0",18x,                            &
             "gamma",19x,"P(2)",18x,"Q(2)",11x,"mtp")
    5 format(3x,1i2,1a2,1x,1p,5(3x,1e19.12),3x,1i5)
    8 format(1x,a,f10.5,a)
    9 format(1x,a,i7,a,i7,a)
      !
      write(stream,*)
      write(stream,*) "Intermediate-state subshell radial wavefunction summary:"
      write(stream,*) "--------------------------------------------------------"
      write(stream,*)
      write(stream,4)
      write(stream,*)
      do  i = 1,wave_intermed%number_of_rwf
         write(stream,5) wave_intermed%rwf(i)%orbital%n,                  &
                     orbital_symmetry(wave_intermed%rwf(i)%orbital%kappa),&
                     wave_intermed%rwf(i)%energy,wave_intermed%rwf(i)%pz, &
                     wave_intermed%rwf(i)%gamma,                          &
                     wave_intermed%rwf(i)%P(2),wave_intermed%rwf(i)%Q(2), &
                     wave_intermed%rwf(i)%mtp
      end do
      !
      write(stream,*)
      write(stream,*) "Final-state subshell radial wavefunction summary:"
      write(stream,*) "-------------------------------------------------"
      write(stream,*)
      write(stream,4)
      write(stream,*)
      do  i = 1,wave_final%number_of_rwf
         write(stream,5) wave_final%rwf(i)%orbital%n,                     &
                     orbital_symmetry(wave_final%rwf(i)%orbital%kappa),   &
                     wave_final%rwf(i)%energy,wave_final%rwf(i)%pz,       &
                     wave_final%rwf(i)%gamma,                             &
                     wave_final%rwf(i)%P(2),wave_final%rwf(i)%Q(2),       &
                     wave_final%rwf(i)%mtp
      end do
      !
      ! Write the list of eigenpair indices
      write(stream,*)
      write(stream,*) 'Initial-state ASF:'
      write(stream,*) '------------------'
      write(stream,*)
      write(stream,6)
      write(stream,*)
      do  i = 1,asf_initial%noasf
         totalJ = asf_initial%asf(i)%totalJ
         write(stream,7) asf_initial%asf(i)%level_No,                     &
                     angular_momentum_string(1*totalJ,4),                 &
                     asf_initial%asf(i)%parity,asf_initial%asf(i)%energy
      end do
    6 format(" Level_No",6x,"J^P",10x,"Energy (a.u.)")
    7 format(2x,i4,8x,a4,a1,6x,1e19.12)
      !
      write(stream,*)
      write(stream,*) 'Intermediate-state ASF:'
      write(stream,*) '-----------------------'
      write(stream,*)
      write(stream,6)
      write(stream,*)
      do  i = 1,asf_intermed%noasf
         totalJ = asf_intermed%asf(i)%totalJ
         write(stream,7) asf_intermed%asf(i)%level_No,                    &
                     angular_momentum_string(1*totalJ,4),                 &
                     asf_intermed%asf(i)%parity,asf_intermed%asf(i)%energy
      end do
      !
      write(stream,*)
      write(stream,*) 'Final-state ASF:'
      write(stream,*) '----------------'
      write(stream,*)
      write(stream,6)
      write(stream,*)
      do  i = 1,asf_final%noasf
         totalJ = asf_final%asf(i)%totalJ
         write(stream,7) asf_final%asf(i)%level_No,                       &
                     angular_momentum_string(totalJ,4),                   &
                     asf_final%asf(i)%parity,asf_final%asf(i)%energy
      end do
      !
   end subroutine print_summary_imf
   !
end module rabs_print
