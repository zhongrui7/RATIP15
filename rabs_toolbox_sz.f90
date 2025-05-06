module rabs_toolbox_sz
!
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains all procedures which help carry out a number of
! 'utility' tasks within the RATIP environment. A few of these procedures 
! occur in a similar form also in other modules of the RATIP package 
! but have been adapted for the present purpose so that the utility 
! program does by itself not depend on too many modules.
!-----------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_hamiltonian
   use rabs_print
   use rabs_relci
   use rabs_toolbox_aux
   implicit none
   !
   public  :: toolbox_self_energy
                 ! Estimates the self-energy corrections for a set of ASF by
                 ! calling the appropriate procedures from RELCI.
   public  :: toolbox_split_csl_list
                 ! Splits a .csl list into individual level groups owing to
                 ! their J^P symmmetry.
   public  :: toolbox_unformat_mix
                 ! Unformats a GRASP92 .mix  file.
   public  :: toolbox_unformat_out
                 ! UNformats a GRASP92 .out  file.
   !
   type(csf_basis), public :: csf_set_a
   !
   type :: toolbox_amplitude_all
      character(len=4) :: typ
      integer          :: asfi, asff, level_i, level_f, totalJ_i, totalJ_f
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: energy, amplitude, amplitude_mj
   end type toolbox_amplitude_all
   !
   integer                                                :: no_amplitude_all
   type(toolbox_amplitude_all), dimension(:), allocatable :: amplitude_all
   !
contains
   !
   subroutine toolbox_self_energy()
   !--------------------------------------------------------------------
   ! Estimates the self-energy corrections due to Yong-Ki Kim's procedures.
   ! The procedure starts from a complete set of ASF which are read in 
   ! by its .csl, .mix as well as .out files.
   !
   ! Calls: file_get_csl_list().
   !--------------------------------------------------------------------
      !
      integer            :: i, ierr, ios, stream
      logical            :: yes
      character(len=6)   :: g92mix
      character(len=7)   :: relci_energy_unit
      character(len=31)  :: record
      character(len=256) :: orb_rwf_file, toolbox_mix_file
      real(kind=dp)      :: rm3, rm1, rp1, rp2
      !
      real(kind=dp), dimension(:), allocatable :: shift
      !
      ! Open, check, load data from, and close the  .iso  file
      call file_get_isodat_grasp2k()
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
    4 continue
      !
      ! Now proceed with the computations
      !
      allocate( shift(1:asf_set%noasf) )
      shift(:) = zero
      !
      relci_energy_unit = "eV     "
      hamiltonian_no_eigenpairs = asf_set%noasf
      !
      call relci_print_energies(stream,3,asf_set,shift)
      !
      stream = 6
      !!x nuclear_charge = 102
      if (stream == 6) then
         write(stream,*) " "
         write(stream,*) "Entering QED ..."
         write(stream,*) " "
         write(stream,*) "   Orbital     Self-energy estimate    " //&
                         "  charge density ratio       F(alphaZ) "
         write(stream,*) "                    (Hartree)          " //&
                         " this orbital / H--like                "
         write(stream,*) " "
         !
         call relci_selfenergy(asf_set,wave,shift,.true.)
         !
         write(stream,*) " "
         write(stream,*) "... QED complete."
      end if
      !
      write(stream,*) " "
      !
      !!! if (relci_add_qed_to_mix) then
      !!! write(stream,*) "Self-energy corrections estimated --- they are "//&
      !!!                 " added to the total energies "
      !!! write(stream,*) " in the RELCI mixing coefficients file."
      !!! else
      !!! write(stream,*) "Self-energy corrections estimated --- these do "//&
      !!!                 "not influence the data "
      !!! write(stream,*) " in the RELCI mixing coefficients file."
      !!! end if
      !
      !!x relci_energy_unit = "eV     "
      !!x hamiltonian_no_eigenpairs = asf_set%noasf
      !
      call relci_print_energies(stream,3,asf_set,shift)
      !
   end subroutine toolbox_self_energy
   !
   !
   subroutine toolbox_split_csl_list()
   !--------------------------------------------------------------------
   ! 'Splits' a given .csl lists into individual lists (level groups)
   ! owing to their J^P symmetry. For each symmetry block, the program 
   ! prompts for a file name but is omitted if no file name is given.
   !
   ! Calls: file_get_csl_list().
   !--------------------------------------------------------------------
      !
      integer            :: i, ierr, nocsf_a, totalJ, split_csf_counter
      character(len=1)   :: parity
      character(len=256) :: csl_old_file, csl_new_file
      character(len=512) :: record1, record2, record3
      !
      logical, dimension(0:40) :: p_symmetry, m_symmetry
      !
      p_symmetry(:) = .false.;   m_symmetry(:) = .false.
      ! 
      ! Open, check, load data from, and close, the  first .csl  file
      call file_get_csl_list(                                                &
         "Enter the name of the original configuration symmetry list file:", &
         csf_set_a,csl_old_file)
      !
      ! Determine the symmetry blocks in the .csl list
      nocsf_a = csf_set_a%nocsf
      do  i = 1,nocsf_a
         totalJ = csf_set_a%csf(i)%totalJ
         parity = csf_set_a%csf(i)%parity
         if (parity == "+") then
            p_symmetry(totalJ) = .true.
         else
            m_symmetry(totalJ) = .true.
         end if
      end do
      !
      ! Re-open the original .csl list again for copying the CSF
      call file_open(21,csl_old_file,"formatted  ","old",ierr)
      if (ierr == 1) then
         stop "toolbox_split_csl_list(): program stop A."
      end if
      !
      ! Now, in turn, prompt for a file name and 'copy' the selected CSF
      do  totalJ = 0,40
         if (p_symmetry(totalJ)) then
            print *, " "
       1    print *, "Enter the name of the .csl list for the symmetry"      //&
                     " block J^P = "//trim(angular_momentum_string(totalJ,4))//&
                     "+ ; use <cr> to skip this block"
            read (*,"(a)") csl_new_file
            if (len(trim(csl_new_file)) == 0) cycle
            !
            call file_open(22,csl_new_file,"formatted  ","new",ierr)
            if (ierr == 1) goto 1
            !
            split_csf_counter = 0
            !
            rewind 21
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
            do  i = 1,nocsf_a
               if (totalJ == csf_set_a%csf(i)%totalJ   .and.   &
                   parity == csf_set_a%csf(i)%parity) then
                  read (21,"(a)") record1
                  read (21,"(a)") record2
                  read (21,"(a)") record3
                  write(22,"(a)") trim(record1)
                  write(22,"(a)") trim(record2)
                  write(22,"(a)") trim(record3)
                  split_csf_counter = split_csf_counter + 1
               else
                  read (21,"(a)") record1
                  read (21,"(a)") record2
                  read (21,"(a)") record3
               end if
            end do
            !
            print *, "Write out complete; there are ",split_csf_counter, &
                     " relativistic CSF in this list."
            !
            close (22)
         end if
      end do
      !
      do  totalJ = 0,40
         if (m_symmetry(totalJ)) then
            print *, " "
       2    print *, "Enter the name of the .csl list for the symmetry"      //&
                     " block J^P = "//trim(angular_momentum_string(totalJ,4))//&
                     "- ; use <cr> to skip this block"
            read (*,"(a)") csl_new_file
            if (len(trim(csl_new_file)) == 0) cycle
            !
            call file_open(22,csl_new_file,"formatted  ","new",ierr)
            if (ierr == 1) goto 2
            !
            split_csf_counter = 0
            !
            rewind 21
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
            do  i = 1,nocsf_a
               if (totalJ == csf_set_a%csf(i)%totalJ   .and.   &
                   parity == csf_set_a%csf(i)%parity) then
                  read (21,"(a)") record1
                  read (21,"(a)") record2
                  read (21,"(a)") record3
                  write(22,"(a)") trim(record1)
                  write(22,"(a)") trim(record2)
                  write(22,"(a)") trim(record3)
                  split_csf_counter = split_csf_counter + 1
               else
                  read (21,"(a)") record1
                  read (21,"(a)") record2
                  read (21,"(a)") record3
               end if
            end do
            !
            print *, "Write out complete; there are ",split_csf_counter, &
                     " relativistic CSF in this list."
            !
            close (22)
         end if
      end do
      !
      close (21)
      !
   end subroutine toolbox_split_csl_list
   !
   !
   subroutine toolbox_unformat_mix()
   !--------------------------------------------------------------------
   ! Converts 'formatted' GRASP92 .mix files into a 'unformatted' one.
   !
   ! Calls: get_yes_stream().
   !--------------------------------------------------------------------
      !
      logical            :: yes
      character(len=30)  :: g92mix
      character(len=256) :: file_name
      integer            :: i, ierr, ios, j, mx, npr, nak
      real(kind=dp)      :: e, pz
      integer, dimension(:), allocatable       :: iatjpo, iaspar
      real(kind=dp), dimension(:), allocatable :: eval
      !
      ! Open the formatted file and read in all data
    2 print *, "Enter a file name for the (formatted) GRASP92 .mix "//&
               "coefficient file:"
      read (*,"(a)") file_name
      call file_open(25,file_name,"formatted  ","old",ierr)
      if (ierr /= 0  .or.  len_trim(file_name) == 0) goto 2
      !
      ! Check the header of the file; if not as expected, try again
      read (25,"(a)",iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX (formatted file version") then
         print *, "ios, g92mix = ",ios, g92mix
         print *, "Not a formatted GRASP92 Mixing Coefficients File;"
         close (25)
         goto 2
      end if
      !
      ! Allocate memory for the set of ASF and load the data from the 
      ! .mix file
      print *, "Loading Mixing Coefficients File ..."
      read (25,"(3i6)") asf_set%csf_set%number_of_electrons, &
                        asf_set%csf_set%nocsf, asf_set%csf_set%nwshells
      read (25,"(i6)")  asf_set%noasf
      allocate( asf_set%asf(1:asf_set%noasf) )
      do  i = 1,asf_set%noasf
         allocate( asf_set%asf(i)%eigenvector(1:asf_set%csf_set%nocsf) )
      end do
      !
      read (25,"(100i6)") (asf_set%asf(i)%level_No,i = 1,asf_set%noasf)
      read (25,"(100(i5,a1))") (asf_set%asf(i)%totalJ,asf_set%asf(i)%parity, &
                          i = 1,asf_set%noasf)
      !
      read (25,"(101e16.9)") asf_set%average_energy,                   &
                             (asf_set%asf(i)%energy,i = 1,asf_set%noasf)
      do  j = 1,asf_set%csf_set%nocsf
         read (25,"(100e16.9)") (asf_set%asf(i)%eigenvector(j), &
                                 i=1,asf_set%noasf)
      end do
      print *, ' ... load complete;'
      close (25)
      !
      ! Open an 'uformatted' file and write out all data
      !
    3 print *, "Enter a file name for the 'unformatted' .mix "//&
               "Coefficient File:"
      read (*,"(a)") file_name
      call file_open(25,file_name,"unformatted  ","new",ierr)
      if (ierr /= 0  .or.  len_trim(file_name) == 0) goto 3
      !
      ! Write out the file header
      write (25) "G92MIX"
      !
      print *, "Write unformatted Mixing Coefficients File ..."
      write (25) asf_set%csf_set%number_of_electrons, &
                 asf_set%csf_set%nocsf, asf_set%csf_set%nwshells
      write (25) asf_set%noasf         
      !
      allocate( iatjpo(1:asf_set%noasf), iaspar(1:asf_set%noasf), &
                eval(1:asf_set%noasf) ) 
      do  i = 1,asf_set%noasf
         iatjpo(i) = asf_set%asf(i)%totalJ + 1
         eval(i)   = asf_set%asf(i)%energy
         if (asf_set%asf(i)%parity == "+") then
            iaspar(i) = 1
         else if (asf_set%asf(i)%parity == "-") then
            iaspar(i) = -1
         else
            stop "toolbox_unformat(): program stop A."
         end if
      end do
      !
      write (25) (asf_set%asf(i)%level_No,i = 1,asf_set%noasf)
      write (25) (iatjpo(i),iaspar(i),i = 1,asf_set%noasf)
      write (25) asf_set%average_energy,(eval(i),i = 1,asf_set%noasf)
      !x write (25) ((evec(i+(j-1)*ncf),i = 1,ncf),j = 1,nvec)
      write (25) ((asf_set%asf(i)%eigenvector(j),j=1,asf_set%csf_set%nocsf),&
                 i=1,asf_set%noasf)
      !
      deallocate( iatjpo, iaspar, eval )
      print *, ' ... write out complete;'
      close (25)
      !
      do  i = 1,asf_set%noasf
         deallocate( asf_set%asf(i)%eigenvector )
      end do
      deallocate( asf_set%asf )
      close(25);   close(26)
      !
   end subroutine toolbox_unformat_mix
   !
   !
   subroutine toolbox_unformat_out()
   !--------------------------------------------------------------------
   ! Converts 'formatted' GRASP92 .out file into 'unformatted' files.
   !
   ! Calls: get_yes_stream().
   !--------------------------------------------------------------------
      !
      logical            :: yes
      character(len=30)  :: g92rwf
      character(len=256) :: file_name
      integer            :: i, ierr, j, mx, npr, nak
      real(kind=dp)      :: e, pz
      real(kind=dp), dimension(:), allocatable :: eval, grid, p_component, &
                                                  q_component
         !
      ! Attempt to open the formatted file and read in all data
    4 print *, "Enter the name of the (formatted) GRASP92 Radial "//&
               "Wavefunction File:"
      read (*,"(a)") file_name
      call file_open(25,file_name,"formatted  ","old",ierr)
      if (ierr /= 0  .or.  len_trim(file_name) == 0) goto 4
      !
      ! Check the header of the file; if not as expected, try again
      read (25,"(a)") g92rwf
      if (g92rwf /= "G92RWF (formatted file version") then
         close (25)
         print *, "g92rwf = ",g92rwf
         print *, "Not a formatted GRASP92 radial wavefunction file;"
         goto 4
      end if
      !
      ! At this point the  .rwf  file has been opened and checked;
      ! Open a formatted file
    5 print *, "Enter a file name for the 'unformatted' .out Radial "// &
               "Wavefunction File:"
      read (*,"(a)") file_name
      call file_open(26,file_name,"unformatted","new",ierr)
      if (ierr /= 0  .or.  len_trim(file_name) == 0) goto 5
      !
      ! Write out a file header
      write (26) "G92RWF"
      !
      ! Read (formatted) and write (unformatted) the parameters of each 
      ! orbital
    6 read (25,*,end =7) npr, nak, e, mx
      write (26) npr, nak, e, mx
      allocate( p_component(1:mx), q_component(1:mx), grid(1:mx) )
      read (25,*) pz
      do  i = 1,mx
         read (25,*) grid(i),p_component(i),q_component(i)
      end do
      write (26)   pz, (p_component(i),i = 1,mx),(q_component(i),i = 1,mx)
      write (26) (grid(i),i = 1,mx)
      deallocate( p_component, q_component, grid )
      goto 6
    7 print *, "Conversion of Radial Wavefunction File complete."
      !
      close (25);   close(26)
      !
   end subroutine toolbox_unformat_out
   !
end module rabs_toolbox_sz
