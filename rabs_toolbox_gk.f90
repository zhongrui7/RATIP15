module rabs_toolbox_gk
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains procedures which help carry out a number of
! 'utility' tasks within the RATIP environment. 
!-----------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_input_dialog
   use rabs_file_handling
   use rabs_nucleus
   use rabs_grasp2k
   use rabs_multipole
   use rabs_toolbox_aux
   implicit none
   !
   public  :: toolbox_gather_energies
                 ! Collects all level information from one or several 
		 ! .mix files.
   public  :: toolbox_grid_calculator
                 ! Determines proper grid parameters if different information
		 ! about the grid size or free-electron energy is given.
   !
   type :: toolbox_isotope
      real(kind=dp) :: Z, A, mass, rms
      real(kind=dp) :: nu
   end type toolbox_isotope
   type(toolbox_isotope), dimension(20) :: isotope
   !
   type :: toolbox_isotope_triples
      integer	    :: i1, i2, i3
      real(kind=dp) :: M, F
   end type toolbox_isotope_triples
   type(toolbox_isotope_triples), dimension(30) :: triples
   !
   type(asf_hfs_basis) :: asf_hfs_i, asf_hfs_f 
   !
contains
   !
   subroutine toolbox_gk()
   !--------------------------------------------------------------------
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      print *, "**************************"
      print *, "*** Not yet implmented ***"
      print *, "**************************"
      !
   end subroutine toolbox_gk
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
         read (25,"(61e26.19)")   (eval(j),j = 181,noasf)
         do  j = 301,noasf
            asf_mix(i)%asf(j)%energy = asf_mix(i)%average_energy + eval(j)
         end do
	 !
         do  j = 1,asf_mix(i)%csf_set%nocsf
            read (25,"(60e16.9)") wa
         end do
         !
         if (asf_set%noasf <= 360) goto 4
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
   !
   subroutine toolbox_grid_calculator()
   !--------------------------------------------------------------------
   ! Determines proper grid parameters if different information about the 
   ! grid size or free-electron energy is given.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer       :: i, sel, dn
      real(kind=dp) :: r_max, e_kin, wavelength
      !
    5 print *, " "
      print *, "Select a particular task:"
      print *, "   1 - Display the grid for given parameters."
      print *, "   2 - Determine  R( N; input: rnt, h, hp, N)."
      print *, "   3 - Determine  N( R_max; input: rnt, h, hp)."
      print *, "   4 - Determine  hp & N     = f(input: rnt, h, R_max, E_kin)."
      print *, "   5 - Determine  hp & R_max = f(input: rnt, h, N, E_kin)."
      print *, "   6 - Determine  E_kin & N  = f(input: rnt, h, hp, R_max)."
      print *, "   7 - Determine  E_kin & R_max = f(input: rnt, h, hp, N)."
      read  *, sel
      !
      ! Determine grid parameters
      call input_grid_parameters("standard")
      !
      select case(sel)
      case(1)
         !
	 ! Display the grid for given parameters 
	 ! -------------------------------------
         !
         call input_grid_parameters("modify")
         !
         ! Generate the radial grid and all associated arrays
         call radgrd_grasp2k()
         !
	 ! Write out the selected grid
	 print *, "Enter the stepsize Delta-n to print the grid: "
	 read  *, dn
	 !
	 write (*,7)
	 do  i = 1,n_grasp2k,dn
            write(*,8) i, r_grasp2k(i),rp_grasp2k(i),rpor_grasp2k(i)
	 end do
       7 format( /"Grid is defined as follows:",                            &
                 /"---------------------------",                            &
                //" --i-- -------- r -------- -------- r' -------",         &
                  " ------- r'/r ------")
       8 format(i6, 1p,6(1x,1d19.12))
      case(2)
         !
	 ! Determine  R( N; input: rnt, h, hp, N)
	 ! --------------------------------------
         !
         call input_grid_parameters("modify")
         !
         ! Generate the radial grid and all associated arrays
         call radgrd_grasp2k()
         !
	 write(*,*) " "
	 write(*,*) " R_max(N) = ",r_grasp2k(n_grasp2k)
      case(3)
         !
	 ! Determine  N( R_max; input: rnt, h, hp)
	 ! ---------------------------------------
         !
         print *, "Enter rnt:"
         read  *, rnt_grasp2k
         print *, "Enter h:"
         read  *, h_grasp2k
         print *, "enter hp:"
         read  *, hp_grasp2k
	 !
         print *, "enter R_max:"
         read  *, r_max
	 !
	 n_grasp2k = 30000
	 !
         ! Generate the radial grid and all associated arrays
         call radgrd_grasp2k()
         !
	 ! determine a proper N
	 do  i = 1,n_grasp2k
	    if (r_grasp2k(i) > r_max) then
	       write(*,*) " "
	       write(*,*) " N(R_max~), R_max~ = ", i,r_grasp2k(i)
	       return
	    end if
	 end do
	 !
         stop "toolbox_grid_calculator(): program stop A."
      case(4)
         !
	 ! Determine  hp & N     = f(input: rnt, h, R_max, E_kin)
	 ! ------------------------------------------------------
         !
         print *, "Enter rnt:"
         read  *, rnt_grasp2k
         print *, "Enter h:"
         read  *, h_grasp2k
	 !
         print *, "enter R_max:"
         read  *, r_max
	 !
         ! Determine the units for the printout and further optional input
         call input_energy_unit()
         print *, "enter E_kin (in "//trim(energy_unit)//"):"
         read  *, e_kin
	 !
         if (energy_inverse) then
            e_kin = energy_factor * e_kin
         else
	    e_kin = e_kin / energy_factor
         end if
	 !!x print *, "e_kin(a.u.) = ",e_kin
	 !
         wavelength = sqrt( two * pi * pi / e_kin )
	 hp_grasp2k = wavelength / 100.0_dp
	 n_grasp2k  = 30000
	 !
         ! Generate the radial grid and all associated arrays
         call radgrd_grasp2k()
         !
	 ! determine a proper N
	 do  i = 1,n_grasp2k
	    if (r_grasp2k(i) > r_max) then
	       write(*,*) " "
	       write(*,9) " N(R_max~), hp, R_max~ = ", i,hp_grasp2k,r_grasp2k(i)
	       return
	    end if
	 end do
       9 format(a, i6, 1p,6(1x,1d19.12))
	 !
         stop "toolbox_grid_calculator(): program stop B."
      case(5)
         !
	 ! Determine  hp & R_max = f(input: rnt, h, N, E_kin)
	 ! --------------------------------------------------
         !
         print *, "Enter rnt:"
         read  *, rnt_grasp2k
         print *, "Enter h:"
         read  *, h_grasp2k
         print *, "enter N:"
         read  *, n_grasp2k
	 !
         ! Determine the units for the printout and further optional input
         call input_energy_unit()
         print *, "enter E_kin (in "//trim(energy_unit)//"):"
         read  *, e_kin
	 !
         if (energy_inverse) then
            e_kin = energy_factor * e_kin
         else
	    e_kin = e_kin / energy_factor
         end if
	 !!x print *, "e_kin(a.u.) = ",e_kin
	 !
         wavelength = sqrt( two * pi * pi / e_kin )
	 hp_grasp2k = wavelength / 100.0_dp
	 !
         ! Generate the radial grid and all associated arrays
         call radgrd_grasp2k()
         !
	 write(*,*) " "
	 write(*,9) " N, R_max(N), hp = ", n_grasp2k, r_grasp2k(n_grasp2k), &
	                                   hp_grasp2k
      case(6)
         !
	 ! Determine  E_kin & N  = f(input: rnt, h, hp, R_max)
	 ! ---------------------------------------------------
         !
         print *, "Enter rnt:"
         read  *, rnt_grasp2k
         print *, "Enter h:"
         read  *, h_grasp2k
         print *, "enter hp:"
         read  *, hp_grasp2k
	 !
         print *, "enter R_max:"
         read  *, r_max
	 !
	 n_grasp2k  = 30000
	 !
         ! Generate the radial grid and all associated arrays
         call radgrd_grasp2k()
	 !
	 ! determine a proper N
	 do  i = 1,n_grasp2k
	    if (r_grasp2k(i) > r_max) then
	       exit
	    end if
	 end do
	 !
	 if (hp_grasp2k > eps10) then
	    wavelength = ten*ten * hp_grasp2k
	 else
	    wavelength = ten * ten * (r_grasp2k(i)-r_grasp2k(i-1))
	 end if
	 e_kin         = two * pi * pi / (wavelength * wavelength)
	 !
	 write(*,*) " "
	 write(*,9) " N(R_max~), R_max~, E_kin (in Hartree) = ",            &
	              i,r_grasp2k(i), e_kin
      case(7)
         !
	 ! Determine  E_kin & R_max = f(input: rnt, h, hp, N)
	 ! ---------------------------------------------------
         !
         call input_grid_parameters("modify")
	 !
         ! Generate the radial grid and all associated arrays
         call radgrd_grasp2k()
	 !
	 if (hp_grasp2k > eps10) then
	    wavelength = ten*ten * hp_grasp2k
	 else
	    wavelength = ten*ten * (r_grasp2k(n_grasp2k)-r_grasp2k(n_grasp2k-1))
	 end if
	 e_kin         = two * pi * pi / (wavelength * wavelength)
	 !
	 write(*,*) " "
	 write(*,9) " N, R_max(N), E_kin (in Hartree) = ",                 &
	              n_grasp2k, r_grasp2k(n_grasp2k), e_kin
      case default
         print *, "Selection cannot be recognized ... redo "
	 goto 5
      end select
      !
      !
   end subroutine toolbox_grid_calculator
   !
end module rabs_toolbox_gk
