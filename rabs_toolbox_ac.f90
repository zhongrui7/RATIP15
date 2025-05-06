module rabs_toolbox_ac
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
   use rabs_lsj
   use rabs_lsj_data
   use rabs_print
   use rabs_toolbox_aux
   implicit none
   !
   public  :: toolbox_charge_density
                 ! Controls the calculation of the radial charge densities
		 ! for a single orbital.
   public  :: toolbox_charge_density_orbital
                 ! Controls the calculation of the radial charge densities.
   public  :: toolbox_condense_csl_list
                 ! Condensation of a .csl list on a single weight criterium.
   public  :: toolbox_condense_on_levels
                 ! Condensation of a .csl list for some given levels only.
   public  :: toolbox_control_weights_lsj 
                 ! Controls the transformation of jj-coupled ASF into a
		 ! LS-coupled CSF basis.
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
   subroutine toolbox_ac()
   !--------------------------------------------------------------------
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      print *, "**************************"
      print *, "*** Not yet implmented ***"
      print *, "**************************"
      !
   end subroutine toolbox_ac
   !
   !
   subroutine toolbox_charge_density()
   !--------------------------------------------------------------------
   ! Determines the effective radial charge or charge density of some 
   ! selected level. The effective charge is given by the integral over 
   ! the charge density up to a given radius R; it has been applied to 
   ! estimate the 'anti-screening' in the projectile. The routine also
   ! supports to print the charge density of a particular level
   ! to some text file.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer            :: i, ierr, ios, j,k, level, mtp
      logical            :: yes
      character(len=6)   :: g92mix
      character(len=256) :: record
      character(len=256) :: toolbox_mix_file, toolbox_rwf_file, &
                            cdensity_file
      real(kind=dp)      :: wa, R, Zeff
      real(kind=dp), dimension(:), allocatable :: rho, ta
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
      ! Get the occupied shells of the initial and final states
      call file_get_csl_list(                                       &
         "Enter the name of the configuration symmetry list file:", &
         asf_set%csf_set)
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
         goto 3
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
         call load_mix_file_grasp2k(asf_set,.true.,ierr)
         if (ierr /= 0) then
            close (25)
            goto 1
         end if
         close (25)
      end if
      !
      ! Initialize and load the radial wavefunctions; allow for both,
      ! formatted and unformatted orbital files 
    3 call initialize_rwf_storage(asf_set, wave) 
      !!x call toolbox_initialize_rwf_storage(.true.,.false.,.false.)
      !
    4 print *, "Enter the name of the Radial WaveFunction File:"
      read (*,"(a)") toolbox_rwf_file
      if (len(trim(toolbox_rwf_file)) == 0) goto 4
      !
      call file_open(21,toolbox_rwf_file,"unformatted","old",ierr)
      if (ierr == 1) goto 4
      !
      read (21,iostat = ios) record(1:6)
      if (record(1:6) /= 'G92RWF') then
         close (21)
         !
         call file_open(21,toolbox_rwf_file,"formatted  ","old",ierr)
         if (ierr == 1) goto 1
         read (21,"(a31)",iostat = ios) record
         if (ios /= 0   .or.   &
             record(1:31) /= "G92RWF (formatted file version)") then
            print *, "ios, record(1:31) = ",ios, record(1:31)
            print *, "Not a G92RWF Radial WaveFunction File;"
            close (21)
            goto 4
         end if
         toolbox_use_formatted_rwf_file = .true.
      else
         toolbox_use_formatted_rwf_file = .false.
      endif
      !
      ! Load data from the  .rwf  file
      call load_rwf_file_grasp2k(wave,toolbox_use_formatted_rwf_file,ierr)
      if (ierr /= 0) then
         close (21)
         goto 4
      end if
      !
      ! Close the  .rwf  file
      close (21)
      !
      print *, "Enter the level number of the selected level:"
      read *,  level
      do  i = 1,asf_set%noasf
         if (asf_set%asf(i)%level_No == level) exit
      end do
      if (i == asf_set%noasf + 1) then
         stop "toolbox_control_cdensity(): program stop A." 
      end if
      !
      allocate( rho(1:n_grasp2k+10), ta(1:n_grasp2k+10) )
      rho(:) = zero
      do  j = 1,asf_set%csf_set%nwshells
         wa = zero
         do  k = 1,asf_set%csf_set%nocsf
            wa = wa + asf_set%asf(i)%eigenvector(k) *    &
                      asf_set%asf(i)%eigenvector(k) *    &
                      asf_set%csf_set%csf(k)%occupation(j)
         end do
         if (wave%rwf(j)%orbital%n /= asf_set%csf_set%subshell(j)%n  .or. &
         wave%rwf(j)%orbital%kappa /= asf_set%csf_set%subshell(j)%kappa) then
            stop "toolbox_control_cdensity(): program stop B." 
         end if
         mtp = wave%rwf(j)%mtp
         rho(1:mtp) = rho(1:mtp) + wa *                             &
                      (wave%rwf(j)%P(1:mtp)*wave%rwf(j)%P(1:mtp) +  & 
                       wave%rwf(j)%Q(1:mtp)*wave%rwf(j)%Q(1:mtp))  
      end do
      !
      print *, "Calculate the effective charge at one or several points ?"
      yes = get_yes_stream()
      if (yes) then
       5 print *, "Enter another effective radius R > 0 (0. if done):"
         read (*,*,end=6)  R
         if (R < zero   .or.   R > r_grasp2k(n_grasp2k)) then
            print *, "The effective radius must be 0 < R <= ", &
                     r_grasp2k(n_grasp2k),"; redo ..."
            goto 5
         else if (R == zero)  then
            goto 6
         end if
         !
         ! Calculate and print the effective charge
         do  i = 1,n_grasp2k-1
            if (r_grasp2k(i+1) > R)  then
               mtp = i
               exit
            end if
         end do
         ta(1:mtp) = rho(1:mtp) * rp_grasp2k(1:mtp)
         Zeff  = quad_grasp2k(ta,mtp)
         !
         ! Estimate the last part by linear interpolation
         Zeff = Zeff + rho(mtp)*(R-r_grasp2k(mtp))
         print *, "Effective charge for R = ",R,"is Z_eff(R) = ",Zeff
         goto 5
      end if
      !
    6 continue
      print *, "Write the charge density to some ASCII file ?"
      yes = get_yes_stream()
      if (yes) then
    7    print *, "Enter a file name for the charge density: "
         read (*,"(a)") cdensity_file
         call file_open(26,cdensity_file,"formatted  ","new",ierr)
         if (ierr /= 0  .or.  len_trim(cdensity_file) == 0) goto 7
         !
         ! Write out a file header
         write(26,*) "This file contains an approximations to the " //&
                     "charge density"
         write(26,*) " "
         write(26,*) "Format used:       "
         write(26,*) "                   "
         write(26,*) "       I   Mesh(I)   Density(I)   Zeff(I)           "
         write(26,*) "----------------------------------------------------"
         write(26,*) " "
         write(26,"(a,i3)") "Number of grid points = ",n_grasp2k
         write(26,*) " "
         !
         mtp       = n_grasp2k
         ta(1:mtp) = rho(1:mtp) * rp_grasp2k(1:mtp)
         do  j = 1,mtp
            ta(1:mtp) = rho(1:mtp) * rp_grasp2k(1:mtp)
            Zeff = quad_grasp2k(ta,j)
            !!x print *, "j,Zeff = ",j,Zeff
            write(26,8) j,r_grasp2k(j),rho(j),Zeff
         end do
         close(26)
      end if
    8 format(2x,i4,2x,1pe12.5,4x,2(1pe12.5,2x) )
      !
   end subroutine toolbox_charge_density
   !
   !
   subroutine toolbox_charge_density_orbital()
   !--------------------------------------------------------------------
   ! Determines the effective radial charge or charge density of some 
   ! given orbital from a radial orbital file. The effective charge is 
   ! given by the integral over the charge density up to a given radius R; 
   ! it has been applied to estimate the 'anti-screening' in the projectile. 
   ! The routine also supports to print the charge density of an orbital
   ! to some text file.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer            :: i, ierr, ios, j,k, level, mtp, n, kappa
      logical            :: yes
      character(len=256) :: record
      character(len=256) :: toolbox_rwf_file, cdensity_file
      real(kind=dp)      :: wa, R, Zeff
      real(kind=dp), dimension(:), allocatable :: rho, ta
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
      ! Get the occupied shells of the initial and final states
      call file_get_csl_list(                                       &
         "Enter the name of the configuration symmetry list file:", &
         asf_set%csf_set)
      !
      ! Initialize and load the radial wavefunctions; allow for both,
      ! formatted and unformatted orbital files
      ! Define a maximal number of orbitals 
      !!x asf_set%csf_set%nwshells = 100
      call initialize_rwf_storage(asf_set, wave) 
      !
    4 print *, "Enter the name of the Radial WaveFunction File:"
      read (*,"(a)") toolbox_rwf_file
      if (len(trim(toolbox_rwf_file)) == 0) goto 4
      !
      call file_open(21,toolbox_rwf_file,"unformatted","old",ierr)
      if (ierr == 1) goto 4
      !
      read (21,iostat = ios) record(1:6)
      if (record(1:6) /= 'G92RWF') then
         close (21)
         !
         call file_open(21,toolbox_rwf_file,"formatted  ","old",ierr)
         if (ierr == 1) goto 4
         read (21,"(a31)",iostat = ios) record
         if (ios /= 0   .or.   &
             record(1:31) /= "G92RWF (formatted file version)") then
            print *, "ios, record(1:31) = ",ios, record(1:31)
            print *, "Not a G92RWF Radial WaveFunction File;"
            close (21)
            goto 4
         end if
         toolbox_use_formatted_rwf_file = .true.
      else
         toolbox_use_formatted_rwf_file = .false.
      endif
      !
      ! Load data from the  .rwf  file
      call load_rwf_file_grasp2k(wave,toolbox_use_formatted_rwf_file,ierr)
      if (ierr /= 0) then
         close (21)
         goto 4
      end if
      !
      ! Close the  .rwf  file
      close (21)
      !
      print *, "Enter the principal and kappa quantum number of the orbital:"
      read *, n, kappa
      !
      ! Check that this orbital exists in the given orbital file
      do  i = 1,asf_set%csf_set%nwshells
         if (wave%rwf(i)%orbital == nkappa(n, kappa)) goto 10
      end do
      stop "toolbox_charge_density_orbital(): program stop A." 
      !
   10 allocate( rho(1:n_grasp2k+10), ta(1:n_grasp2k+10) )
      rho(:) = zero
      mtp = wave%rwf(i)%mtp
      rho(1:mtp) = (wave%rwf(i)%P(1:mtp)*wave%rwf(i)%P(1:mtp) +  & 
                    wave%rwf(i)%Q(1:mtp)*wave%rwf(i)%Q(1:mtp))  
      !
      print *, "Calculate the effective charge at one or several points ?"
      yes = get_yes_stream()
      if (yes) then
       5 print *, "Enter another effective radius R > 0 (0. if done):"
         read (*,*,end=6)  R
         if (R < zero   .or.   R > r_grasp2k(n_grasp2k)) then
            print *, "The effective radius must be 0 < R <= ", &
                     r_grasp2k(n_grasp2k),"; redo ..."
            goto 5
         else if (R == zero)  then
            goto 6
         end if
         !
         ! Calculate and print the effective charge
         do  i = 1,n_grasp2k-1
            if (r_grasp2k(i+1) > R)  then
               mtp = i
               exit
            end if
         end do
         ta(1:mtp) = rho(1:mtp) * rp_grasp2k(1:mtp)
         Zeff  = quad_grasp2k(ta,mtp)
         !
         ! Estimate the last part by linear interpolation
         Zeff = Zeff + rho(mtp)*(R-r_grasp2k(mtp))
         print *, "Effective charge for R = ",R,"is Z_eff(R) = ",Zeff
         goto 5
      end if
      !
    6 continue
      print *, "Write the charge density to some ASCII file ?"
      yes = get_yes_stream()
      if (yes) then
    7    print *, "Enter a file name for the charge density: "
         read (*,"(a)") cdensity_file
         call file_open(26,cdensity_file,"formatted  ","new",ierr)
         if (ierr /= 0  .or.  len_trim(cdensity_file) == 0) goto 7
         !
         ! Write out a file header
         write(26,*) "This file contains an approximations to the " //&
                     "charge density"
         write(26,*) " "
         write(26,*) "Format used:       "
         write(26,*) "                   "
         write(26,*) "       I   Mesh(I)   Density(I)   Zeff(I)           "
         write(26,*) "----------------------------------------------------"
         write(26,*) " "
         write(26,"(a,i3)") "Number of grid points = ",n_grasp2k
         write(26,*) " "
         !
         mtp       = n_grasp2k
         ta(1:mtp) = rho(1:mtp) * rp_grasp2k(1:mtp)
         do  j = 1,mtp
            ta(1:mtp) = rho(1:mtp) * rp_grasp2k(1:mtp)
            Zeff = quad_grasp2k(ta,j)
            !!x print *, "j,Zeff = ",j,Zeff
            write(26,8) j,r_grasp2k(j),rho(j),Zeff
         end do
         close(26)
      end if
    8 format(2x,i4,2x,1pe12.5,4x,2(1pe12.5,2x) )
      !
   end subroutine toolbox_charge_density_orbital
   !
   !
   subroutine toolbox_condense_csl_list()
   !--------------------------------------------------------------------
   ! Condenses a .csl list due to a given .mix mixing coefficient file.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer            :: counter, i, ierr, ios, level
      character(len=6)   :: g92mix
      character(len=256) :: record
      character(len=256) :: csl_new_file, toolbox_csl_file, toolbox_mix_file
      real(kind=dp)      :: w, wa
      !
      logical, dimension(:), allocatable :: keep_csf  
      !
      ! Open, check, load data from, and close, the  .csl  file
    1 print *, "Enter the name of the GRASP92 configuration "//&
               "symmetry list file:"
      read (*,"(a)") toolbox_csl_file
      if (len(trim(toolbox_csl_file)) == 0) goto 1
      !
      call file_open(21,toolbox_csl_file,"formatted  ","old",ierr)
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
      ! Load data from the  .csl  file and close the file
      call load_csl_from_file(asf_set%csf_set)
      close (21)
      !
    2 print *, "Enter the name of corresponding .mix mixing coefficient file:"
      read (*,"(a)") toolbox_mix_file
      if (len(trim(toolbox_mix_file)) == 0) goto 2
      !
      call file_open(25,toolbox_mix_file,"unformatted","old",ierr)
      if (ierr /= 0) goto 2
      !
      ! Check the header of the file; if not as expected for unformatted
      ! files, check formatted form
      read (25,iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         close (25)
         goto 3
      else
         call load_mix_file_grasp2k(asf_set,.false.,ierr)
         if (ierr /= 0) then
            print *, "Not a proper .mix mixing coefficient file for the "//&
                     "given .csl list; reenter ..."
            close (25)
            goto 2
         end if
         goto 4
      end if
      !
      ! Try formatted file format; check the header of the file; 
      ! if not as expected for formatted files, try again
    3 call file_open(25,toolbox_mix_file,"formatted  ","old",ierr)
      if (ierr /= 0) goto 2
      read (25,"(a)",iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         print *, "ios, g92mix = ",ios, g92mix
         print *, "Not a GRASP92 Mixing Coefficients File;"
         close (25)
         goto 2
      else
         call load_mix_file_grasp2k(asf_set,.true.,ierr)
         if (ierr /= 0) then
            close (25)
            goto 2
         end if
      end if
      !
    4 allocate( keep_csf(1:asf_set%csf_set%nocsf) )
      keep_csf(:) = .false.
      !
    5 print *, "Enter a cutoff criterion 0 < w <= 1 to remove all CSF"
      print *, " with a weight < w from the .csl list:"
      read *,  w
      !
      if (w <= zero   .or.   w > one) then
         print *, "The cutoff criterion must be in the range 0 < w <= 1; " //&
                  "reenter ..."
         goto 5
      end if
      !
      do  i = 1,asf_set%csf_set%nocsf
         do  level = 1,asf_set%noasf
            wa = asf_set%asf(level)%eigenvector(i) * &
                 asf_set%asf(level)%eigenvector(i)
            if (wa >= w)  keep_csf(i) = .true.
         end do
      end do
      !
      ! Reopen the original .csl list
      call file_open(21,toolbox_csl_file,"formatted  ","old",ierr)
      !
      ! Create the condensed .csl list
    6 print *, "Enter the name of the configuration symmetry list file" //&
               " that is to be created:"
      read (*,"(a)") csl_new_file
      if (len(trim(csl_new_file)) == 0) goto 1
      !
      call file_open(22,csl_new_file,"formatted  ","new",ierr)
      if (ierr == 1) goto 6
      !
      ! Create the head of the .csl file
      read (21,"(a)") record
      write(22,"(a)") trim(record)
      read (21,"(a)") record
      write(22,"(a)") trim(record)
      read (21,"(a)") record
      write(22,"(a)") trim(record)
      read (21,"(a)") record
      write(22,"(a)") trim(record)
      read (21,"(a)") record
      write(22,"(a)") trim(record)
      !
      counter = 0
      do  i = 1,asf_set%csf_set%nocsf
         if (keep_csf(i)) then
            read (21,"(a)") record
            write(22,"(a)") trim(record)
            read (21,"(a)") record
            write(22,"(a)") trim(record)
            read (21,"(a)") record
            write(22,"(a)") trim(record)
            counter = counter + 1
         else
            read (21,*) 
            read (21,*) 
            read (21,*) 
         end if
      end do
      !
      close(21)
      close(22)
      deallocate( keep_csf )
      !
      print *, "... creates a condensed .csl list with a total number of", &
               counter,"CSF."
      print *, " "
      !
   end subroutine toolbox_condense_csl_list
   !
   !
   subroutine toolbox_condense_on_levels()
   !--------------------------------------------------------------------
   ! Condenses a .csl list due to a given .mix mixing coefficient file.
   ! An extended input may allow to select one or several ASF (which are
   ! used for condensation) and to create a new .mix file for the
   ! condensed list. Although the eigenvectors are 're-normalized' in this
   ! list, they are in general not fully orthogonal to each other.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer            :: counter, i, ierr, ios, j,k, ki, level, &
                            number_of_levels
      logical            :: yes, fail
      character(len=6)   :: g92mix
      character(len=256) :: record
      character(len=256) :: csl_new_file, toolbox_csl_file, toolbox_mix_file, &
                            mix_new_file
      real(kind=dp)      :: w, wa
      !
      integer, dimension(1:1000)         :: levels    
      logical, dimension(:), allocatable :: keep_csf  
      !
      ! Open, check, load data from, and close, the  .csl  file
    1 print *, "Enter the name of the GRASP92 configuration "//&
               "symmetry list file:"
      read (*,"(a)") toolbox_csl_file
      if (len(trim(toolbox_csl_file)) == 0) goto 1
      !
      call file_open(21,toolbox_csl_file,"formatted  ","old",ierr)
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
      ! Load data from the  .csl  file and close the file
      call load_csl_from_file(asf_set%csf_set)
      close (21)
      !
    2 print *, "Enter the name of corresponding .mix mixing coefficient file:"
      read (*,"(a)") toolbox_mix_file
      if (len(trim(toolbox_mix_file)) == 0) goto 2
      !
      call file_open(25,toolbox_mix_file,"unformatted","old",ierr)
      if (ierr /= 0) goto 2
      !
      ! Check the header of the file; if not as expected for unformatted
      ! files, check formatted form
      read (25,iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         close (25)
         goto 3
      else
         call load_mix_file_grasp2k(asf_set,.false.,ierr)
         if (ierr /= 0) then
            print *, "Not a proper .mix mixing coefficient file for the "//&
                     "given .csl list; reenter ..."
            close (25)
            goto 2
         end if
         goto 4
      end if
      !
      ! Try formatted file format; check the header of the file; 
      ! if not as expected for formatted files, try again
    3 call file_open(25,toolbox_mix_file,"formatted  ","old",ierr)
      if (ierr /= 0) goto 2
      read (25,"(a)",iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         print *, "ios, g92mix = ",ios, g92mix
         print *, "Not a GRASP92 Mixing Coefficients File;"
         close (25)
         goto 2
      else
         call load_mix_file_grasp2k(asf_set,.true.,ierr)
         if (ierr /= 0) then
            close (25)
            goto 2
         end if
      end if
      !
    4 allocate( keep_csf(1:asf_set%csf_set%nocsf) )
      keep_csf(:) = .false.
      !
    5 print *, "Enter a cutoff criterion 0 < w <= 1 to remove all CSF"
      print *, " with a weight < w from the .csl list:"
      read *,  w
      !
      if (w <= zero   .or.   w > one) then
         print *, "The cutoff criterion must be in the range 0 < w <= 1; " //&
                  "reenter ..."
         goto 5
      end if
      !
    6 print *, "Enter the level numbers of the selected ASF, "
      print *, " e.g. 1 3 4  7 - 20  48  69 - 85 :"
      read (*, "(a)") record
      call toolbox_interprete_levels(record,levels,number_of_levels,fail)
      if (fail) then
         print *, "Unable to interprete the serial level numbers; redo ..."
         goto 6
      end if
      print *, "number_of_levels, levels(:) = ", &
                number_of_levels, levels(1:number_of_levels)
      !
      do  i = 1,asf_set%csf_set%nocsf
         do j = 1,number_of_levels
            do  level = 1,asf_set%noasf
               if (levels(j) == asf_set%asf(level)%level_No) then
                  wa = asf_set%asf(level)%eigenvector(i) * &
                       asf_set%asf(level)%eigenvector(i)
                  if (wa >= w)  keep_csf(i) = .true.
               end if
            end do
         end do
      end do
      !
      ! Reopen the original .csl list
      call file_open(21,toolbox_csl_file,"formatted  ","old",ierr)
      !
      ! Create the condensed .csl list
    7 print *, "Enter the name of the configuration symmetry list file" //&
               " that is to be created:"
      read (*,"(a)") csl_new_file
      if (len(trim(csl_new_file)) == 0) goto 1
      !
      call file_open(22,csl_new_file,"formatted  ","new",ierr)
      if (ierr == 1) goto 7
      !
      ! Create the head of the .csl file
      read (21,"(a)") record
      write(22,"(a)") trim(record)
      read (21,"(a)") record
      write(22,"(a)") trim(record)
      read (21,"(a)") record
      write(22,"(a)") trim(record)
      read (21,"(a)") record
      write(22,"(a)") trim(record)
      read (21,"(a)") record
      write(22,"(a)") trim(record)
      !
      counter = 0
      do  i = 1,asf_set%csf_set%nocsf
         if (keep_csf(i)) then
            read (21,"(a)") record
            write(22,"(a)") trim(record)
            read (21,"(a)") record
            write(22,"(a)") trim(record)
            read (21,"(a)") record
            write(22,"(a)") trim(record)
            counter = counter + 1
         else
            read (21,*) 
            read (21,*) 
            read (21,*) 
         end if
      end do
      !
      close(21)
      close(22)
      !
      print *, "... creates a condensed .csl list with a total number of", &
               counter,"CSF."
      !!x print *, " "
      !
      print *, "Generate a 'corresponding' .mix file with re-normalized " //&
               "eigenvectors ?"
      yes = get_yes_stream()
      if (yes) then
         asf_final%csf_set%number_of_electrons =                           &
                                         asf_set%csf_set%number_of_electrons
         asf_final%csf_set%nocsf    = counter              
         asf_final%csf_set%nwshells = asf_set%csf_set%nwshells              
         asf_final%noasf            = number_of_levels 
         asf_new%average_energy     = zero
         allocate( asf_final%asf(1:number_of_levels) )
         do  j = 1,number_of_levels          
            allocate( asf_final%asf(j)%eigenvector(1:counter) ) 
         end do         
         asf_final%asf(1:number_of_levels)%level_No = levels(1:number_of_levels)
         print *, "number_of_levels = ",number_of_levels 
         ! 
         iloop: do  i = 1,asf_final%noasf 
            do  j = 1,asf_set%noasf
               if (asf_final%asf(i)%level_No == asf_set%asf(j)%level_No) then
                  asf_final%asf(i)%totalJ = asf_set%asf(j)%totalJ
                  asf_final%asf(i)%parity = asf_set%asf(j)%parity
                  asf_final%asf(i)%energy = asf_set%asf(j)%energy
                  ki = 0;  wa = zero
                  do  k = 1,asf_set%csf_set%nocsf
                     if (keep_csf(k)) then
                        ki = ki + 1
                        asf_final%asf(i)%eigenvector(ki) =                   &
                                                 asf_set%asf(j)%eigenvector(k)
                        wa = wa + asf_final%asf(i)%eigenvector(ki)*          &
                                  asf_final%asf(i)%eigenvector(ki)
                     end if
                  end do
                  print *, "ki, wa = ",ki,wa
                  !
                  ! Now re-normalize
                  wa = one / sqrt(wa)
                  do  ki = 1,asf_final%csf_set%nocsf
                     asf_final%asf(i)%eigenvector(ki) = &
                                         wa * asf_final%asf(i)%eigenvector(ki)
                  end do
                  cycle iloop
               end if
            end do
            !
            print *, "asf_final%asf(i)%level_No = ",asf_final%asf(i)%level_No
            stop "toolbox_control_condense(): program stop A." 
         end do  iloop            
                                         
         !
       8 print *, "Enter a file name for a 'formatted' .mix Coefficient File:"
         read (*,"(a)") mix_new_file
         call file_open(25,mix_new_file,"formatted  ","new",ierr)
         if (ierr /= 0  .or.  len_trim(mix_new_file) == 0) goto 8
         !
         ! Write out a file header
         write (25,"(a)") g92mix//" (formatted file version)."
         !
         print *, "Write formatted Mixing Coefficients File ..."
         write (25,"(3i6)") asf_final%csf_set%number_of_electrons, &
                            asf_final%csf_set%nocsf, asf_final%csf_set%nwshells
         write (25,"(i6)")  asf_final%noasf
         !
         write (25,"(100i6)") (asf_final%asf(i)%level_No,i = 1,asf_final%noasf)
         write (25,"(100(i5,a1))") (asf_final%asf(i)%totalJ,                  &
                                 asf_final%asf(i)%parity,i = 1,asf_final%noasf)
         !
         write (25,"(101e16.9)") asf_final%average_energy, &
                                (asf_final%asf(i)%energy,i = 1,asf_final%noasf)
         do  j = 1,asf_final%csf_set%nocsf
            write (25,"(100e16.9)") (asf_final%asf(i)%eigenvector(j), &
                                     i=1,asf_final%noasf)
         end do
         !
         print *, ' ... write out complete;'
         close (25)
      end if
      deallocate( keep_csf )
      !
   end subroutine toolbox_condense_on_levels
   !
   !
   subroutine toolbox_control_weights_lsj()
   !--------------------------------------------------------------------
   ! Controls the transformation of atomic states from a jj- into a 
   ! LS-coupled CSF basis.
   !
   ! Calls: lsj_deallocate_asf_basis_LS, lsj_form_csf_basis_LS,
   ! lsj_interprete_levels, lsj_print_conf_scheme_LS,
   ! lsj_print_single_config_jj, lsj_print_single_config_LS,
   ! lsj_transformation_ASF, load_csl_from_file, load_mix_file_grasp2k, 
   ! print_configuration_scheme, file_open, util_csl_file.
   !--------------------------------------------------------------------
      !
      integer            :: i, j, jj, number_of_levels, ierr, ios
      logical            :: yes, fail, printlevel_gg1
      character(len=4)   :: wLS, XLS
      character(len=6)   :: g92mix
      character(len=256) :: record
      character(len=256) :: toolbox_csl_file, toolbox_mix_file
      !
      character(len=4), dimension(1:100) :: xxLS
      integer, dimension(1:10000)    :: levels
      integer, dimension(0:10000)    :: leading_LS
      integer, dimension(1:100)      :: iw
      integer                        :: level, nocsf_min, lev, string_length
      integer                        :: nocsf_max, sum_nocsf_min
      real(kind=dp), dimension(1:100):: weights, mcoeffs
      real(kind=dp)                  :: wa, wb
      !
      call lsj_initialization_LS_jj()
      !
      ! Open, check, load data from, and close, the  .csl  file
      call file_get_csl_list(                                             &
         "Enter the name of the configuration symmetry list file:",       &
         asf_set_jj%csf_set)
      !
      ! Get eigenvectors of the Atomic State Functions
      call file_get_mix("Enter the name of the mixing coefficient file:", &
         asf_set_jj)
      !
    4 print *, "Maximum number of considered ASF is:",asf_set_jj%noasf
      print *, "Enter the level numbers of the ASF which are to be transformed,"
      print *, " e.g. 1 3 4  7 - 20  48  69 - 85 :"
      read (*, "(a)") record
      call lsj_interprete_levels(record,levels,number_of_levels,fail)
      if (fail) then
         print *, "Unable to interprete the serial level numbers; redo ..."
         goto 4
      end if
      if (asf_set_jj%noasf < number_of_levels) then
         print *, "There are to much ASF:", number_of_levels
         go to 4
      end if
      !
      print *, " Enter the number of the leading CSF to be printed:"
      read (*, "(i2)") nocsf_max
      call lsj_form_csf_basis_LS()
      !
      ! Determine th amount of `additional printout'
      printlevel_gg1 = .true.
      print *, "Include additional printout due to Gediminas ?"
      printlevel_gg1 = get_yes_stream()
      !
      !
      string_length    = index(toolbox_csl_file,'.')
      print *, "string_length = ",string_length
      toolbox_csl_file = toolbox_csl_file(1:string_length)//'LS'
      toolbox_csl_file = "work.LS"
      call file_open(25,toolbox_csl_file,"formatted  ","new",ierr)
      call lsj_print_conf_scheme_LS(25,asf_set_LS%csf_set_LS)
      close(25)
      !
      allocate(asf_set_LS%asf(1:asf_set_jj%noasf))
      do i = 1, asf_set_jj%noasf
         allocate(asf_set_LS%asf(i)%eigenvector(1:asf_set_LS%csf_set_LS%nocsf))
         asf_set_LS%asf(i)%level_No = i
         asf_set_LS%asf(i)%energy   = asf_set_jj%asf(i)%energy
         asf_set_LS%asf(i)%totalJ   = asf_set_jj%asf(i)%totalJ
         asf_set_LS%asf(i)%parity   = asf_set_jj%asf(i)%parity
      end do
      asf_set_LS%noasf = asf_set_jj%noasf
      !
      jj = 0
      do  lev = 1, number_of_levels
         level = levels(lev)
	 !
	 ! Additional printout due to the original LSJ program
	 if (printlevel_gg1) then
            if (lev > 1) then
               print *, " "
               print *, ".  .  .  .  .  .  .  .  .  .  .  .  .  ."
               print *, " "
               print *, " "
               print *, " The new level is under the investigation."
            end if
            print *, " "
            print *, "Weights of major contributors to ASF in jj-coupling:"
            print *, " "
            print *, " Level  J Parity      CSF contributions"
            print *, " "
            !
            weights(1:100) = zero;   iw(1:100) = 0
            wb = zero
            do  i = 1,asf_set_jj%csf_set%nocsf
               wa = asf_set_jj%asf(level)%eigenvector(i) * &
                    asf_set_jj%asf(level)%eigenvector(i)
               wb = wb + asf_set_jj%asf(level)%eigenvector(i) * &
                         asf_set_jj%asf(level)%eigenvector(i)
	       !
               do  j = 1,99
                  if (wa > weights(j)) then
                     weights(j+1:100) = weights(j:99)
                     weights(j)       = wa
                     iw(j+1:100)      = iw(j:99)
                     iw(j)            = i
                     exit
                  end if
               end do
            end do
            !
            if (level > 1   .and.   abs(wb) > 1.0001) then
               print *, "level, wb = ",level,wb
               stop "toolbox_control_weights_lsj(): program stop A." 
            end if
            !
            nocsf_min = 5
            do j = 1,5
               if(abs(weights(j)) < 0.00001) then
                  nocsf_min = j - 1
                  exit
               end if
            end do
	    !
            print 18, asf_set_jj%asf(level)%level_No,                           &
            trim(angular_momentum_string(1*asf_set_jj%asf(level)%totalJ,4)),    &
                 asf_set_jj%asf(level)%parity,(weights(j),iw(j),j=1,nocsf_min)
            print *, " "
            print *, "Definition of leading CSF:"
            print *, " "
            call lsj_print_single_config_jj(-1,asf_set_jj%csf_set,iw(1))
            print *, " "
            print*, " Total sum over   weight    is:  ",wb
         18 format(1x,i4,1x,2a4,2x,100(3x,f8.5," of",i5))
	 end if
         !
         call lsj_transformation_ASF(level)
      end do
      !
      !
      !
      print *, " "
      print *, " "
      print *, "Mixing coefficients of the major LS-coupled CSF to ASF"
      print *, "------------------------------------------------------"
      print *, " "
      print *, "  Level  J/Parity   1 - sum(wa)      CSF contributions"
      print *, " "
      !
      do  lev = 1, number_of_levels
         level = levels(lev)
         weights(1:100) = zero;   mcoeffs(1:100) = zero;   iw(1:100) = 0
         wb = zero
         do  i = 1,asf_set_LS%csf_set_LS%nocsf
            wa =      asf_set_LS%asf(level)%eigenvector(i) * &
                      asf_set_LS%asf(level)%eigenvector(i)
            wb = wb + asf_set_LS%asf(level)%eigenvector(i) * &
                      asf_set_LS%asf(level)%eigenvector(i)
            do  j = 1,99
               if (wa > weights(j)) then
                  weights(j+1:100) = weights(j:99)
                  weights(j)       = wa
                  mcoeffs(j+1:100) = mcoeffs(j:99)
		  mcoeffs(j)       = asf_set_LS%asf(level)%eigenvector(i)
                  iw(j+1:100)      = iw(j:99)
                  iw(j)            = i
                  exit
               end if
            end do
         end do
         !
         if (level > 1   .and.   abs(wb) > 1.0001) then
            print *, "level, wb = ",level,wb
            stop "toolbox_control_weights_lsj(): program stop B." 
         end if
         !
         nocsf_min = nocsf_max
         do j = 1,nocsf_max
            if(abs(weights(j)) < 0.00001) then
               nocsf_min = j - 1
               exit
            end if
            jj = jj + 1
            leading_LS(jj) = iw(j)
         end do
	 !
         sum_nocsf_min = sum_nocsf_min + nocsf_min
         if (nocsf_min <= 0) then
            stop "toolbox_control_weights_lsj(): program stop C." 
         end if
	 !
         print 6,  asf_set_LS%asf(level)%level_No,                   &
                   trim(angular_momentum_string                      &
                   (1*asf_set_LS%csf_set_LS%csf(iw(1))%totalJ,4)),   &
                   asf_set_LS%csf_set_LS%csf(iw(1))%parity,wb,       &
                   (mcoeffs(j),iw(j),J=1,nocsf_min)
         !!x print*, "                     Total sum over  weight  is:",wb
         !!x print *, " "
         asf_set_LS%asf(i)%max_csf_No = iw(1)
       6 format(1x,i4,6x,a4,1x,a1,1x,e14.6,2x,5(3x,f8.5," of",i5),   &
                /34x,5(3x,f8.5," of",i5),                            &
                /34x,5(3x,f8.5," of",i5))
      end do
      !
      !
      !
      print *, " "
      print *, " "
      print *, "Definition of the main LS-coupled CSF"
      print *, "-------------------------------------"
      print *, " "
      do  i = 1,asf_set_jj%csf_set%nocsf
         do  j = 1, sum_nocsf_min
            if (i == leading_LS(j)) then
               call lsj_print_single_config_LS                  &
                       (-1,asf_set_LS%csf_set_LS,leading_LS(j))
               exit
            end if
         end do
      end do
      !
      !
      print *, " "
      print *, " "
      print *, "Weights of the major LS-coupled CSF to ASF"
      print *, "------------------------------------------"
      print *, " "
      print *, "  Level  J/Parity   1 - sum(wa)      CSF contributions"
      print *, " "
      !
      do  lev = 1, number_of_levels
         level = levels(lev)
	 !
         weights(1:100) = zero;   iw(1:100) = 0
         wb = zero
	 !
         do  i = 1,asf_set_LS%csf_set_LS%nocsf
            wa =      asf_set_LS%asf(level)%eigenvector(i) * &
                      asf_set_LS%asf(level)%eigenvector(i)
            wb = wb + asf_set_LS%asf(level)%eigenvector(i) * &
                      asf_set_LS%asf(level)%eigenvector(i)
	    !
            do  j = 1,99
               if (wa > weights(j)) then
                  weights(j+1:100) = weights(j:99)
                  weights(j)       = wa
                  iw(j+1:100)      = iw(j:99)
                  iw(j)            = i
                  exit
               end if
            end do
         end do
         !
         if (level > 1   .and.   abs(wb) > 1.0001) then
            print *, "level, wb = ",level,wb
            stop "toolbox_control_weights_lsj(): program stop B." 
         end if
         !
         nocsf_min = nocsf_max
         do j = 1,nocsf_max
            if(abs(weights(j)) < 0.00001) then
               nocsf_min = j - 1
               exit
            end if
            jj = jj + 1
            leading_LS(jj) = iw(j)
         end do
         sum_nocsf_min = sum_nocsf_min + nocsf_min
         if (nocsf_min <= 0) then
            stop "lsj_control_transformation(): program stop C."
         end if
	 !
	 do  j = 1,nocsf_min
         call lsj_spectroscopic_LS(iw(j),asf_set_LS%csf_set_LS%nwshells,wLS,XLS)
	 xxLS(j) = XLS
	 end do
	 !
	 print 16, asf_set_LS%asf(level)%level_No,		  &
         	trim(angular_momentum_string			  &
         	(1*asf_set_LS%csf_set_LS%csf(iw(1))%totalJ,4)),   &
         	asf_set_LS%csf_set_LS%csf(iw(1))%parity,wb,	  &
         	(weights(j),iw(j),trim(xxLS(j)),J=1,nocsf_min)
         asf_set_LS%asf(i)%max_csf_No = iw(1)
      16 format(1x,i4,6x,a4,1x,a1,1x,e14.6,2x,                        &
                     4(3x,f8.5," of",i5," (",a3,")"),                 &
                /34x,4(3x,f8.5," of",i5," (",a3,")"),                 &
                /34x,4(3x,f8.5," of",i5," (",a3,")"))
      end do
      !
      call lsj_deallocate_asf_basis_LS(asf_set_LS)
      if (debuging /= 0 )then
         close(57)
      end if
      !
   end subroutine toolbox_control_weights_lsj
   !
end module rabs_toolbox_ac
