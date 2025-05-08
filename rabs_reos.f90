module rabs_reos
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains all procedures which are specific to the REOS
! program. This includes the calculation of the one-particle reduced
! matrix elements but also several procedures for file handling and the
! intermediate storage of several quantities to accelerate the
! computations.
!-----------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   use rabs_csl
   use rabs_determinant
   use rabs_dirac_orbital
   use rabs_file_handling
   use rabs_grasp2k
   use rabs_input_dialog
   use rabs_nucleus
   use omp_lib
   implicit none
   !
   public  :: reos_test_initial
                 ! Opens a .xpn file for the initial states and reads in all
		 ! necessary data from this file.
   private :: reos_amplitude_different_MM
                 ! Calculates the multipole amplitudes for the case that
		 ! Slater determinants with different M-projection occur.
   private :: reos_amplitude_same_MM
                 ! Calculates the multipole amplitudes for the case that
		 ! Slater determinants with the same M-projection occur.
   private :: reos_calculate_Bessel
                 ! Calculates the Bessel function over the radial grid
                 ! for a given factor.
   public  :: reos_calculate_amplitudes
                 ! Initializes the computation of transition probabilities
                 ! and lifetimes and calculates the transition amplitudes for
                 ! all transitions; several computational modes are supported.
   private :: reos_calculate_me
                 ! Calculates the contributions of the one-particle multipole
                 ! operators for all defined transitions for a given pair of
                 ! Slater determinants.
   private :: reos_calculate_me_orth
                 ! Calculates the contributions of the one-particle multipole
                 ! operators for a given pair of determinants by assuming
                 ! orthogonal orbital sets for the initial and final states.
   public  :: reos_collect_input
                 ! Collects and proceeds all input for the expansion of the
                 ! symmetry-adapted functions into determinants.
   private :: reos_convert_probability
                 ! Returns the oscillator strength, the Einstein A and B
                 ! coefficients, and the decay width for a given transition
                 ! amplitude.
   public  :: reos_initialize_rwf_storage
   		 ! Initializes the arrays of type(grasp2k_orbital) for the
   		 ! storage of the radial wave functions.
   private :: reos_load_trn_file
		 ! Opens a .trn  transition amplitude file and reads in data
		 ! for all individual transitions.
   private :: reos_open_restart
                 ! Opens the two .res restart files for restarting the REOS
                 ! from a previously aborted run on stream 31 and 32.
   public  :: reos_print_results
                 ! Writes the transition probabilities and lifetimes to
                 ! the .sum file.
   public  :: reos_print_summary
		 ! Appends a summary of the input data to the .sum file.
   private :: reos_print_transitions
                 ! Prints all selected transitions in a neat format before
                 ! the computation starts.
   private :: reos_readwrite_dump
                 ! Reads or writes all necessary information from the REOS
		 ! program from or to a file except of the transitions.
   private :: reos_readwrite_transitions
                 ! Reads or writes the array transition() of
		 ! type(reos_transition) from or to a file.
   private :: reos_reduced_M_integral
                 ! Calculates a M integral for the coupling of the radiation
                 ! field for a given multipolarity and gauge form.
   public  :: reos_select_multipoles
                 ! Determines a list of multipole line components which
                 ! contribute and are considered for a given transition.
   private :: reos_set_determinant_sym
                 ! Allocates storage and 'initializes' information about the
                 ! symmetry blocks of each determinant to accelerate the
                 ! computations.
   private :: reos_set_overlaps
                 ! Initializes the array of overlap integrals for the
                 ! calculation of relaxed transition probabilities.
   private :: reos_set_radial_integrals
                 ! Initializes 'storage' for all required radial multipole
                 ! integrals if these integrals are calculated only ones
                 ! at the beginning of the computations.
   private :: reos_set_transitions
                 ! Determines which radiative transitions need to be calculated
                 ! and initializes storage for an appropriate data type for
                 ! them.
   public  :: reos_write_trn_file
                 ! Writes out transition energies and amplitudes for further
                 ! processing to a .trn file.
   !
   ! Define some global data of the REOS program; most of these data are
   ! read in during the interactive control and may overwrite existing
   ! default values
   !
   ! Storage for the initial and final atomic states and wave functions
   type(asf_det_basis), public   :: asf_initial, asf_final, asf_test
   type(grasp2k_orbital), public :: wave_initial, wave_final
   !
   ! Define an internal structure type(reos_transition) which stores all
   ! necessary information for an radiative transition line
   type :: multipole_line
      character(len=2) :: multipole
      character(len=9) :: gauge
      real(kind=dp)    :: amplitude
      real(kind=dp), dimension(:,:,:,:), pointer :: radial_int
   end type multipole_line
   !
   type :: reos_transition
      integer          :: asfi, asff, level_i, level_f, totalJ_i, totalJ_f
      integer          :: number_of_mlines
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: energy
      real(kind=dp), dimension(:), pointer :: bessel0, bessel1, bessel2, &
                                              bessel3, bessel4, bessel5
      type(multipole_line), dimension(:), pointer :: mline
   end type reos_transition
   !
   type(reos_transition), dimension(:), allocatable :: transition
   !
   ! Define global counters for the control of the REOS program; the cutoff
   ! defines a criterium to neglect the calculation of matrix elements for
   ! two determinants with a small admixture only
   integer       :: number_of_multipoles      = 0,        &
                    number_of_transitions     = 0,        &
                    number_of_symmetry_blocks = 0,        &
                    reos_determinant_counter  = 100000000,&
                    !reos_determinant_counter = 100000,   &
                    reos_dump_period          = 100000000,&
                    reos_dump_billion         = 0,        &
		    reos_dump_billion_done    = 0,        &
                    reos_dump_counter         = 0,        &
                    reos_dump_counter_done    = 0
   !
   real(kind=dp) :: reos_determinant_cutoff   = 1.0e-8_dp
   !
   ! Define global logical flags for the control of the REOS program; the
   ! default values for these flags may be overwritten interactively during
   ! input time
   logical, public :: reos_use_compact_xpn_file      = .false.,  &
                      reos_apply_exp_energies        = .false.,  &
                      reos_apply_restart             = .false.,  &
                      reos_assume_orthogonality      = .false.,  &
                      reos_enable_restart            = .false.,  &
                      reos_precalculate_integrals    = .true.,   &
                      reos_sort_transition_energy    = .true.,   &
                      reos_print_AB_in_hartree       = .false.,  &
                      reos_print_selected_trans      = .true.,   &
                      reos_use_formatted_rwf_file    = .true.,   &
                      reos_write_transition_file     = .false.
   !
   ! Define storage for the overlap integrals: (kappa,pqn_f,pqn_i)
   real(kind=dp), dimension(:,:,:), allocatable :: reos_overlap
   !
   ! Define storage for symmetries of the determinants;
   ! reos_detsym(determinant_index,symmetry_index) where symmetry_index if
   ! found in sym_block(kappa,mm); the present version therefore only
   ! allows symmetries with |kappa| <= 7
   integer, dimension(-19:19,-37:37) :: sym_block  = 0
   !!x integer, dimension(-7:7,-13:13) :: sym_block  = 0
   integer(kind=i1b), dimension(:,:), allocatable ::                    &
                    reos_detsym_i, reos_detsym_f
   !
   ! Define logical flags for debugging individual procedures
   logical, private :: debug_reduced_M_integral       = .false.
   !
   ! Define some variables and arrays for processing input data from
   ! reos_collect_input()
   character(len=2), dimension(20) :: reos_multipole  = "  "
   !
contains
   !
   subroutine reos_test_initial()
   !--------------------------------------------------------------------
   ! Opens a .xpn file for the initial states and reads in all necessary
   ! data from this file. It also close the file finally.
   !
   ! Calls: file_open(), reos_load_xpn_file().
   !--------------------------------------------------------------------
      !
      integer                 :: i, ia, ib, ierr, ios, j, noint, nobit
      integer, dimension(20)  :: iparity
      integer, dimension(200) :: occupation
      character(len=26)       :: record
      character(len=256)      :: reos_xpn_file
      !
    1 print *, "Enter the name of the corresponding CESD .xpn-old file:"
      read (*,"(a)") reos_xpn_file
      if (len(trim(reos_xpn_file)) == 0) goto 1
      !
      call file_open(21,reos_xpn_file,"formatted  ","old",ierr)
      if (ierr == 1) goto 1
      !
      ! Check the first record of the file; if not as expected, try again
      read (21,"(a26)",iostat = ios) record
      if (ios /= 0   .or.   record(1:26) /= "ASF-based CESD output file") then
         print *, "ios, record(1:26) = ",ios, record(1:26)
         print *, "Not a ASF-based CESD output file;"
         close (21)
         goto 1
      end if
      print *, "a"
      !
      ! Load data from the  .xpn-old  file
      read(21,*);   read(21,*);   read(21,*);    read(21,*)
      read(21,*) asf_test%noasf
      allocate( asf_test%asf(1:asf_test%noasf) )
      print *, "b"
      read(21,*) asf_test%det_set%norbital
      allocate( asf_test%det_set%orbital(1:asf_test%det_set%norbital) )
      print *, "c"
      do  i = 1,asf_test%det_set%norbital
         read(21,"(5x,3i4)")  asf_test%det_set%orbital(i)%n,  &
            asf_test%det_set%orbital(i)%kappa, asf_test%det_set%orbital(i)%mm
      end do
      print *, "d"
      !
      read(21,*);   read(21,*);   read(21,*)
      read(21,*) asf_test%det_set%nod
      print *, "e"
      allocate( asf_test%det_set%determinant(1:asf_test%det_set%nod) )
      nobit = bit_size(noint)
      if (mod(asf_test%det_set%norbital,nobit) == 0) then
         noint = asf_test%det_set%norbital/nobit
      else
         noint = asf_test%det_set%norbital/nobit + 1
      end if
      asf_test%det_set%noint = noint
      do  i = 1,asf_test%det_set%nod
         allocate( asf_test%det_set%determinant(i)%occupation(1:noint) )
         read(21,"(10x,15( 10i2,1x ))")  &
            (occupation(j),j=1,asf_test%det_set%norbital)
         call pack_occupation_in_integer(asf_test%det_set%determinant(i), &
                              noint,occupation,asf_test%det_set%norbital)
      end do
      !
      read(21,*);   read(21,*);   read(21,*);    read(21,*)
      read(21,*) asf_test%average_energy
      read(21,*);   read(21,*);   read(21,*)
      if (asf_test%noasf > 8) asf_test%noasf = 8
      read(21,*) (asf_test%asf(i)%level_No,i=1,asf_test%noasf)
      read(21,*) (asf_test%asf(i)%energy,i=1,asf_test%noasf)
      asf_test%asf(1:asf_test%noasf)%energy = &
         asf_test%asf(1:asf_test%noasf)%energy + asf_test%average_energy
      read(21,*) (asf_test%asf(i)%totalJ,i=1,asf_test%noasf)
      asf_test%asf(1:asf_test%noasf)%totalJ  = &
         asf_test%asf(1:asf_test%noasf)%totalJ - 1
      read(21,*)
      read(21,*) (iparity(i),i=1,asf_test%noasf)
      do  i = 1,asf_test%noasf
         if (iparity(i) == 1) then
            asf_test%asf(i)%parity = "+"
         else
            asf_test%asf(i)%parity = "-"
         end if
         allocate( asf_test%asf(i)%eigenvector(1:asf_test%det_set%nod) )
      end do
      !
      do  i = 1,asf_test%det_set%nod
         read(21,*) (asf_test%asf(j)%eigenvector(i),j=1,asf_test%noasf)
      end do
      !
      ! Close the  .xpn  file
      close (21)
      !
      if (asf_initial%noasf            <  asf_test%noasf             .or. &
          asf_initial%det_set%norbital /= asf_test%det_set%norbital  .or. &
          asf_initial%det_set%nod      /= asf_test%det_set%nod       .or. &
          asf_initial%det_set%noint    /= asf_test%det_set%noint) then
         print *, "asf_initial%noasf, asf_test%noasf = ", &
                   asf_initial%noasf, asf_test%noasf
         print *, "asf_initial%det_set%norbital, asf_test%det_set%norbital = ",&
                   asf_initial%det_set%norbital, asf_test%det_set%norbital
         print *, "asf_initial%det_set%nod, asf_test%det_set%nod = ",&
                   asf_initial%det_set%nod, asf_test%det_set%nod
         print *, "asf_initial%det_set%noint, asf_test%det_set%noint = ",&
                   asf_initial%det_set%noint, asf_test%det_set%noint
         stop "reos_test_initial() - a"
      end if
      !
      do  i = 1,asf_initial%det_set%norbital
      if (asf_initial%det_set%orbital(i)%n /= asf_test%det_set%orbital(i)%n .or. &
          asf_initial%det_set%orbital(i)%kappa /= asf_test%det_set%orbital(i)%kappa .or. &
          asf_initial%det_set%orbital(i)%mm /= asf_test%det_set%orbital(i)%mm) &
          then
         print *, "i, na, nb, kapa, kapb, ma, mb = ",                        &
            asf_initial%det_set%orbital(i)%n, asf_test%det_set%orbital(i)%n, &
            asf_initial%det_set%orbital(i)%kappa,                            &
            asf_test%det_set%orbital(i)%kappa,                               &
            asf_initial%det_set%orbital(i)%mm, asf_test%det_set%orbital(i)%mm
         stop "reos_test_initial() - b"
      end if
      end do
      !
      ialoop: do  ia = 1,asf_initial%det_set%nod
         ibloop: do  ib = 1,asf_test%det_set%nod
            do i = 1,asf_initial%det_set%noint
               if (asf_initial%det_set%determinant(ia)%occupation(i) /= &
                      asf_test%det_set%determinant(ib)%occupation(i))   &
                      cycle ibloop
            end do
            print *, "ia, ib = ", ia, ib
            do  i = 1,asf_test%noasf
               if (abs(asf_initial%asf(i)%eigenvector(ia) - &
                       asf_test%asf(i)%eigenvector(ib)) > ten*ten*eps10) then
                  print *, "   eigenvector: i, weighta, weightb, dw = ", &
                     i, asf_initial%asf(i)%eigenvector(ia),  &
                           asf_test%asf(i)%eigenvector(ib),  &
                        asf_initial%asf(i)%eigenvector(ia) - &
                        asf_test%asf(i)%eigenvector(ib)
               end if
            end do
            cycle ialoop
         end do ibloop
         print *, "*****WARNING******ia, ib = ", ia, ib
      end do ialoop
      !
   end subroutine reos_test_initial
   !
   !
   subroutine reos_amplitude_different_MM()
   !--------------------------------------------------------------------
   ! Calculates the multipole amplitudes for the case that Slater
   ! determinants with different M-projection occurs in the wave function
   ! expansion. This often occurs for the calculation of transition
   ! probabilities between different multipletts.
   !
   ! Calls:
   !--------------------------------------------------------------------
      integer :: diff_symmetry, i, ia, j, ki, kf, nocc_i, nocc_f,                  &
                 MM_i, MM_f, MM_photon, norb_i, norb_f, rank,                  &
                 reos_pair_billion, reos_pair_counter, reos_calc_billion,      &
                 reos_calc_counter, reos_neglect_billion, reos_neglect_counter
      integer, dimension(:), allocatable   :: detsym_i, detsym_f, &
                                              occupation_i, occupation_f
      real(kind=dp) :: geom_factor, wa, total
      type(nkappam), dimension(:), allocatable :: orbital_i, orbital_f

      ! Thread-private variables for parallel region
      integer :: local_reos_calc_billion, local_reos_calc_counter, &
                 local_reos_neglect_billion, local_reos_neglect_counter

      ! Allocate storage for working arrays outside parallel region
      nocc_f = asf_final%det_set%noint * &
               bit_size(asf_final%det_set%determinant(1)%occupation(1))
      nocc_i = asf_initial%det_set%noint * &
               bit_size(asf_initial%det_set%determinant(1)%occupation(1))
      allocate( occupation_f(1:nocc_f), occupation_i(1:nocc_i),         &
                detsym_f(number_of_symmetry_blocks),                    &
                detsym_i(number_of_symmetry_blocks),                    &
                orbital_i(1:asf_initial%det_set%number_of_electrons),   &
                orbital_f(1:asf_final%det_set%number_of_electrons) )

      ! Initialize global counters
      reos_dump_billion    = 0;   reos_dump_counter    = 0
      reos_pair_billion    = 0;   reos_pair_counter    = 1
      reos_calc_billion    = 0;   reos_calc_counter    = 0
      reos_neglect_billion = 0;   reos_neglect_counter = 0

      total = (asf_final%det_set%nod+zero) * asf_initial%det_set%nod * &
              number_of_transitions
      print *, "Less than ",total," matrix elements between final- and"// &
               " initial-state determinants need to be calculated."
      print *, " This includes zero contributions due to different"//     &
               " M projections in the expansion of the wavefunctions."

      ! Parallelize the outer loop over transitions
      !$OMP PARALLEL PRIVATE(i, kf, ki, j, wa, MM_f, MM_i, MM_photon, &
      !$OMP                  diff_symmetry, norb_i, norb_f, detsym_i, detsym_f, &
      !$OMP                  occupation_i, occupation_f, orbital_i, orbital_f, &
      !$OMP                  local_reos_calc_billion, local_reos_calc_counter, &
      !$OMP                  local_reos_neglect_billion, local_reos_neglect_counter, &
      !$OMP                  ia) &
      !$OMP          SHARED(reos_dump_billion, reos_dump_counter, &
      !$OMP                  reos_pair_billion, reos_pair_counter, &
      !$OMP                  reos_calc_billion, reos_calc_counter, &
      !$OMP                  reos_neglect_billion, reos_neglect_counter, &
      !$OMP                  asf_final, asf_initial, transition)
      ! Initialize thread-local counters
      local_reos_calc_billion    = 0
      local_reos_calc_counter    = 0
      local_reos_neglect_billion = 0
      local_reos_neglect_counter = 0

      !$OMP DO SCHEDULE(DYNAMIC)
      do i = 1, number_of_transitions
         ! Determine M projections for this transition
         MM_f = -999
         do kf = 1, asf_final%det_set%nod
            if (asf_final%asf(transition(i)%asff)%eigenvector(kf) /= zero) then
               MM_f = asf_final%det_set%determinant(kf)%totalM
               exit
            end if
         end do
         do kf = 1, asf_final%det_set%nod
            if (asf_final%asf(transition(i)%asff)%eigenvector(kf) /= zero .and. &
                asf_final%det_set%determinant(kf)%totalM /= MM_f) then
               stop "reos_amplitude_different_MM(): program stop A."
            end if
         end do

         MM_i = -999
         do ki = 1, asf_initial%det_set%nod
            if (asf_initial%asf(transition(i)%asfi)%eigenvector(ki) /= zero) then
               MM_i = asf_initial%det_set%determinant(ki)%totalM
               exit
            end if
         end do
         do ki = 1, asf_initial%det_set%nod
            if (asf_initial%asf(transition(i)%asfi)%eigenvector(ki) /= zero .and. &
                asf_initial%det_set%determinant(ki)%totalM /= MM_i) then
               stop "reos_amplitude_different_MM(): program stop B."
            end if
         end do
         MM_photon = MM_f - MM_i

         ! Nested loop over determinants for this transition
         do kf = 1, asf_final%det_set%nod
            do ki = 1, asf_initial%det_set%nod
               ! Update pair counters atomically
               !$OMP ATOMIC
               reos_pair_counter = reos_pair_counter + 1
               if (reos_pair_counter >= 1000000000) then
                  !$OMP CRITICAL(pair_billion)
                  reos_pair_billion = reos_pair_billion + 1
                  reos_pair_counter = 0
                  !$OMP END CRITICAL
               end if

               !$OMP ATOMIC
               reos_dump_counter = reos_dump_counter + 1
               if (reos_dump_counter >= 1000000000) then
                  !$OMP CRITICAL(dump_billion)
                  reos_dump_billion = reos_dump_billion + 1
                  reos_dump_counter = 0
                  !$OMP END CRITICAL
               end if

               ! Periodic status output
               if (mod(reos_pair_counter, reos_determinant_counter) == 0) then
                  !$OMP CRITICAL(print_status)
                  print *, "  ",reos_calc_billion,"billion",reos_calc_counter, &
                           "pairs calculated;",reos_neglect_billion,"billion", &
                           reos_neglect_counter,"neglected or zero contribution;", &
                           reos_pair_billion,"billion",reos_pair_counter,"in total;"
                  !$OMP END CRITICAL
               end if

               ! Skip if restarting and pair already processed
               if (reos_apply_restart) then
                  if (reos_dump_billion < reos_dump_billion_done .or. &
                      (reos_dump_billion == reos_dump_billion_done .and. &
                       reos_dump_counter < reos_dump_counter_done)) cycle
               end if

               ! Dump transition information
               if (reos_enable_restart .and. &
                   mod(reos_dump_counter, reos_dump_period) == 0) then
                  !$OMP CRITICAL(dump_transitions)
                  rewind (32)
                  call reos_readwrite_transitions(32, .false.)
                  rewind (32)
                  print *, "dump transitions after ",reos_dump_billion, &
                           " billion and ",reos_dump_counter," pairs on stream 32."
                  !$OMP END CRITICAL
               end if

               ! Check weight threshold
               wa = asf_final%asf(transition(i)%asff)%eigenvector(kf) * &
                    asf_initial%asf(transition(i)%asfi)%eigenvector(ki)
               if (abs(wa) < reos_determinant_cutoff) then
                  local_reos_neglect_counter = local_reos_neglect_counter + 1
                  if (local_reos_neglect_counter >= 1000000000) then
                     local_reos_neglect_billion = local_reos_neglect_billion + 1
                     local_reos_neglect_counter = 0
                  end if
                  cycle
               else
                  local_reos_calc_counter = local_reos_calc_counter + 1
                  if (local_reos_calc_counter >= 1000000000) then
                     local_reos_calc_billion = local_reos_calc_billion + 1
                     local_reos_calc_counter = 0
                  end if
               end if

               ! Symmetry check
               detsym_i(:) = reos_detsym_f(kf,:) - reos_detsym_i(ki,:)
               diff_symmetry = sum(abs(detsym_i))
               if (diff_symmetry > 2) cycle

               ! Define overlap determinant
               call unpack_occupation_from_integer(asf_final%det_set%determinant(kf), &
                                                   occupation_f, nocc_f)
               call unpack_occupation_from_integer(asf_initial%det_set%determinant(ki), &
                                                   occupation_i, nocc_i)
               norb_i = 0; norb_f = 0
               do ia = 1, asf_final%det_set%norbital
                  if (occupation_f(ia) == 1) then
                     norb_f = norb_f + 1
                     orbital_f(norb_f) = asf_final%det_set%orbital(ia)
                  end if
               end do
               do ia = 1, asf_initial%det_set%norbital
                  if (occupation_i(ia) == 1) then
                     norb_i = norb_i + 1
                     orbital_i(norb_i) = asf_initial%det_set%orbital(ia)
                  end if
               end do
               if (rabs_use_stop .and. &
                   (norb_i /= norb_f .or. norb_i /= asf_initial%det_set%number_of_electrons)) then
                  stop "reos_amplitude_different_MM(): program stop C."
               end if

               ! Calculate matrix elements
               if (reos_assume_orthogonality) then
                  call reos_calculate_me_orth(kf, ki, orbital_f, orbital_i, norb_i, &
                                              MM_photon, i)
               else
                  detsym_f(:) = reos_detsym_f(kf,:)
                  detsym_i(:) = reos_detsym_i(ki,:)
                  call reos_calculate_me(kf, ki, orbital_f, orbital_i, norb_i, &
                                         detsym_f, detsym_i, MM_photon, i)
               end if
            end do
         end do
      end do
      !$OMP END DO

      ! Reduce local counters to global counters
      !$OMP CRITICAL(reduce_counters)
      reos_calc_billion    = reos_calc_billion + local_reos_calc_billion
      reos_calc_counter    = reos_calc_counter + local_reos_calc_counter
      reos_neglect_billion = reos_neglect_billion + local_reos_neglect_billion
      reos_neglect_counter = reos_neglect_counter + local_reos_neglect_counter
      !$OMP END CRITICAL

      ! Calculate reduced parts of amplitudes for this thread's transitions
      do i = 1, number_of_transitions
         MM_f = -999
         do kf = 1, asf_final%det_set%nod
            if (asf_final%asf(transition(i)%asff)%eigenvector(kf) /= zero) then
               MM_f = asf_final%det_set%determinant(kf)%totalM
               exit
            end if
         end do
         MM_i = -999
         do ki = 1, asf_initial%det_set%nod
            if (asf_initial%asf(transition(i)%asfi)%eigenvector(ki) /= zero) then
               MM_i = asf_initial%det_set%determinant(ki)%totalM
               exit
            end if
         end do
         MM_photon = MM_f - MM_i

         do j = 1, transition(i)%number_of_mlines
            select case(transition(i)%mline(j)%multipole)
            case("E1", "M1");   rank = 1
            case("E2", "M2");   rank = 2
            case("E3", "M3");   rank = 3
            case("E4", "M4");   rank = 4
            case("E5", "M5");   rank = 5
            case default; stop "reos_amplitude_different_MM(): program stop D."
            end select
            geom_factor = Wigner_Eckardt_geometry(transition(i)%totalJ_f, MM_f, &
                                                  rank+rank, MM_photon, transition(i)%totalJ_i, MM_i)
            if (abs(geom_factor) < eps10) then
               !$OMP CRITICAL(error_print)
               print *, "i, Jf, MMf, rank, MMp, Ji, MMi = ", &
                        i, transition(i)%totalJ_f, MM_f, &
                        rank, MM_photon, transition(i)%totalJ_i, MM_i
               print *, "geom_factor = ", geom_factor
               stop "reos_amplitude_different_MM(): program stop E."
               !$OMP END CRITICAL
            end if
            !$OMP CRITICAL(update_amplitude)
            transition(i)%mline(j)%amplitude = &
            transition(i)%mline(j)%amplitude / geom_factor
            !$OMP END CRITICAL
         end do
      end do
      !$OMP END PARALLEL
   end subroutine reos_amplitude_different_MM
   !
   !
   subroutine reos_amplitude_same_MM()
   !--------------------------------------------------------------------
   ! Calculates the multipole amplitudes for all transitions for the
   ! case that all Slater determinants have the same M projections
   ! in the initial and final states.
   !
   ! Calls: spherical_Bessel_jL().
   !--------------------------------------------------------------------
      integer :: diff_symmetry, i, j, ki, kf, nocc_i, nocc_f,                  &
                 MM_i, MM_f, MM_photon, norb_i, norb_f, rank,                  &
                 reos_pair_billion, reos_pair_counter, reos_calc_billion,      &
                 reos_calc_counter, reos_neglect_billion, reos_neglect_counter
      integer, dimension(:), allocatable   :: detsym_i, detsym_f, &
                                              occupation_i, occupation_f
      real(kind=dp) :: geom_factor, wa, weight, total
      type(nkappam), dimension(:), allocatable :: orbital_i, orbital_f

      ! Thread-private variables for parallel region
      integer :: local_reos_calc_billion, local_reos_calc_counter, &
                 local_reos_neglect_billion, local_reos_neglect_counter

      ! Allocate storage for working arrays outside parallel region
      nocc_f = asf_final%det_set%noint * &
               bit_size(asf_final%det_set%determinant(1)%occupation(1))
      nocc_i = asf_initial%det_set%noint * &
               bit_size(asf_initial%det_set%determinant(1)%occupation(1))
      allocate( occupation_f(1:nocc_f), occupation_i(1:nocc_i),         &
                detsym_f(number_of_symmetry_blocks),                    &
                detsym_i(number_of_symmetry_blocks),                    &
                orbital_i(1:asf_initial%det_set%number_of_electrons),   &
                orbital_f(1:asf_final%det_set%number_of_electrons) )

      ! Set M quantum numbers for the transition amplitudes
      MM_i      = asf_initial%det_set%determinant(1)%totalM
      MM_f      = asf_final%det_set%determinant(1)%totalM
      MM_photon = MM_f - MM_i

      ! Initialize global counters
      reos_dump_billion    = 0;   reos_dump_counter    = 0
      reos_pair_billion    = 0;   reos_pair_counter    = 1
      reos_calc_billion    = 0;   reos_calc_counter    = 0
      reos_neglect_billion = 0;   reos_neglect_counter = 0

      total = (asf_final%det_set%nod+zero) * asf_initial%det_set%nod
      print *, "About ",total," pairs of final- and initial-state "// &
               "need to be calculated."

      ! Parallelize the outer loop over final-state determinants
      !$OMP PARALLEL PRIVATE(kf, ki, i, j, weight, wa, diff_symmetry, &
      !$OMP                  norb_i, norb_f, detsym_i, detsym_f, &
      !$OMP                  occupation_i, occupation_f, orbital_i, orbital_f, &
      !$OMP                  local_reos_calc_billion, local_reos_calc_counter, &
      !$OMP                  local_reos_neglect_billion, local_reos_neglect_counter) &
      !$OMP          SHARED(reos_dump_billion, reos_dump_counter, &
      !$OMP                  reos_pair_billion, reos_pair_counter, &
      !$OMP                  reos_calc_billion, reos_calc_counter, &
      !$OMP                  reos_neglect_billion, reos_neglect_counter, &
      !$OMP                  asf_final, asf_initial, transition, MM_photon)
      ! Initialize thread-local counters
      local_reos_calc_billion    = 0
      local_reos_calc_counter    = 0
      local_reos_neglect_billion = 0
      local_reos_neglect_counter = 0

      !$OMP DO SCHEDULE(DYNAMIC)
      do kf = 1,asf_final%det_set%nod
         do ki = 1,asf_initial%det_set%nod
            ! Update pair counter (atomic to avoid race condition)
            !$OMP ATOMIC
            reos_pair_counter = reos_pair_counter + 1
            if (reos_pair_counter >= 1000000000) then
               !$OMP CRITICAL
               reos_pair_billion = reos_pair_billion + 1
               reos_pair_counter = 0
               !$OMP END CRITICAL
            end if

            !$OMP ATOMIC
            reos_dump_counter = reos_dump_counter + 1
            if (reos_dump_counter >= 1000000000) then
               !$OMP CRITICAL
               reos_dump_billion = reos_dump_billion + 1
               reos_dump_counter = 0
               !$OMP END CRITICAL
            end if

            ! Periodic status output and dump (only one thread at a time)
            if (mod(reos_pair_counter, reos_determinant_counter) == 0) then
               !$OMP CRITICAL(print_status)
               print *, "  ",reos_calc_billion,"billion",reos_calc_counter, &
                        "pairs calculated;",reos_neglect_billion,"billion", &
                        reos_neglect_counter,"neglected;",                 &
                        reos_pair_billion,"billion",reos_pair_counter,     &
                        "in total;"
               !$OMP END CRITICAL
            end if

            ! Dump transition information if required
            if (reos_enable_restart .and. &
                mod(reos_dump_counter, reos_dump_period) == 0) then
               !$OMP CRITICAL(dump_transitions)
               rewind (32)
               call reos_readwrite_transitions(32, .false.)
               rewind (32)
               print *, "dump transitions after ",reos_dump_billion, &
                        " billion and ",reos_dump_counter," pairs on stream 32."
               !$OMP END CRITICAL
            end if

            ! Skip if restarting and this pair was already processed
            if (reos_apply_restart) then
               if (reos_dump_billion < reos_dump_billion_done .or. &
                   (reos_dump_billion == reos_dump_billion_done .and. &
                    reos_dump_counter < reos_dump_counter_done)) cycle
            end if

            ! Calculate maximum weight across transitions
            weight = zero
            do i = 1, number_of_transitions
               wa = asf_final%asf(transition(i)%asff)%eigenvector(kf) * &
                    asf_initial%asf(transition(i)%asfi)%eigenvector(ki)
               weight = max(weight, abs(wa))
            end do

            if (weight < reos_determinant_cutoff) then
               local_reos_neglect_counter = local_reos_neglect_counter + 1
               if (local_reos_neglect_counter >= 1000000000) then
                  local_reos_neglect_billion = local_reos_neglect_billion + 1
                  local_reos_neglect_counter = 0
               end if
               cycle
            else
               local_reos_calc_counter = local_reos_calc_counter + 1
               if (local_reos_calc_counter >= 1000000000) then
                  local_reos_calc_billion = local_reos_calc_billion + 1
                  local_reos_calc_counter = 0
               end if
            end if

            ! Check symmetry blocks
            detsym_i(:) = reos_detsym_f(kf,:) - reos_detsym_i(ki,:)
            diff_symmetry = sum(abs(detsym_i))
            if (diff_symmetry > 2) cycle

            ! Define the 'overlap determinant'
            call unpack_occupation_from_integer(asf_final%det_set%determinant(kf), &
                                                occupation_f, nocc_f)
            call unpack_occupation_from_integer(asf_initial%det_set%determinant(ki), &
                                                occupation_i, nocc_i)
            norb_i = 0; norb_f = 0
            do i = 1, asf_final%det_set%norbital
               if (occupation_f(i) == 1) then
                  norb_f = norb_f + 1
                  orbital_f(norb_f) = asf_final%det_set%orbital(i)
               end if
            end do
            do i = 1, asf_initial%det_set%norbital
               if (occupation_i(i) == 1) then
                  norb_i = norb_i + 1
                  orbital_i(norb_i) = asf_initial%det_set%orbital(i)
               end if
            end do
            if (rabs_use_stop .and. &
                (norb_i /= norb_f .or. norb_i /= asf_initial%det_set%number_of_electrons)) then
               stop "reos_calculate_amplitudes(): program stop A."
            end if

            ! Calculate matrix elements
            if (reos_assume_orthogonality) then
               call reos_calculate_me_orth(kf, ki, orbital_f, orbital_i, norb_i, &
                                           MM_photon, 0)
            else
               detsym_f(:) = reos_detsym_f(kf,:)
               detsym_i(:) = reos_detsym_i(ki,:)
               call reos_calculate_me(kf, ki, orbital_f, orbital_i, norb_i, &
                                      detsym_f, detsym_i, MM_photon, 0)
            end if
         end do
      end do
      !$OMP END DO

      ! Reduce local counters to global counters
      !$OMP CRITICAL(reduce_counters)
      reos_calc_billion    = reos_calc_billion + local_reos_calc_billion
      reos_calc_counter    = reos_calc_counter + local_reos_calc_counter
      reos_neglect_billion = reos_neglect_billion + local_reos_neglect_billion
      reos_neglect_counter = reos_neglect_counter + local_reos_neglect_counter
      !$OMP END CRITICAL
      !$OMP END PARALLEL

      ! Calculate the 'reduced parts' of the amplitudes (serial section)
      do i = 1, number_of_transitions
         do j = 1, transition(i)%number_of_mlines
            select case(transition(i)%mline(j)%multipole)
            case("E1", "M1");   rank = 1
            case("E2", "M2");   rank = 2
            case("E3", "M3");   rank = 3
            case("E4", "M4");   rank = 4
            case("E5", "M5");   rank = 5
            case default; stop "reos_calculate_amplitudes(): program stop B."
            end select
            geom_factor = Wigner_Eckardt_geometry(transition(i)%totalJ_f, MM_f, &
                              rank+rank, MM_photon, transition(i)%totalJ_i, MM_i)
            if (abs(geom_factor) < eps10) then
               print *, "i, Jf, MMf, rank, MMp, Ji, MMi = ", &
                        i, transition(i)%totalJ_f, MM_f, &
                        rank, MM_photon, transition(i)%totalJ_i, MM_i
               print *, "geom_factor = ", geom_factor
               stop "reos_calculate_amplitudes(): program stop C."
            end if
            !$OMP CRITICAL(update_amplitude)
            transition(i)%mline(j)%amplitude = &
            transition(i)%mline(j)%amplitude / geom_factor
            !$OMP END CRITICAL
         end do
      end do
   end subroutine reos_amplitude_same_MM
   !
   !
   subroutine reos_calculate_Bessel(L,omega_over_c,bessel)
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
         bessel(i) = spherical_Bessel_jL(L, abs(omega_over_c*r_grasp2k(i)))
      end do
      !
   end subroutine reos_calculate_bessel
   !
   !
   subroutine reos_calculate_amplitudes()
   !--------------------------------------------------------------------
   ! Initializes the computation of relaxed-orbital transition probabilities
   ! and lifetimes and calculates the transition amplitudes for all
   ! transitions; several computational modes are supported due to to set
   ! up of various logical flags. The main modes are
   !    (i)   full relaxed-orbital transition probabilities,
   !    (ii)  transition probabilities assuming orthogonal orbitals,
   !    (iii) transition probabilities with the 'same' set of radial
   !          integrals for all transitions, for both modes: relaxed-
   !          orbitals or assumed orthogonal orbitals.
   !
   ! Calls: call reos_calculate_me_orth(), call reos_calculate_me(),
   !        reos_print_transitions(), reos_set_determinant_sym(),
   !        reos_set_overlaps(), reos_set_radial_integrals(),
   !        reos_set_transitions(),unpack_occupation_from_integer(),
   !        Wigner_Eckardt_geometry().
   !--------------------------------------------------------------------
      !
      logical :: same_MM
      integer :: i, MM_i, MM_f
      !
      if (reos_apply_restart) then
         reos_dump_billion_done = reos_dump_billion
         reos_dump_counter_done = reos_dump_counter
      else
         !
         print *, "Initialize the set-up of the overlap integrals and "// &
                  "transitions ..."
         call reos_set_transitions()
         call reos_set_overlaps()
         if (.not.reos_assume_orthogonality) then
            call reos_set_determinant_sym()
         end if
         if (reos_precalculate_integrals) then
            call reos_set_radial_integrals()
         end if
         !
         print *, "   ... initialization complete."
         print *, " "
         !
         ! Dump all information on stream 31 if required
         if (reos_enable_restart) then
            print *, "Dump restart information on stream 31 and close "//&
	             "this file."
            call reos_readwrite_dump(31,.false.)
            close (31)
         end if
      end if
      !
      ! Print the selected transitions before the computation starts
      if (reos_print_selected_trans) then
         call reos_print_transitions(6)
         call reos_print_transitions(24)
      end if
      !
      ! If all determinents of the initial respectively the final states
      ! have the same M-projection, all selected transitions can be calculated
      ! within a single run through all pairs of determinant.
      same_MM = .true.;
      MM_i    = asf_initial%det_set%determinant(1)%totalM
      MM_f    = asf_final%det_set%determinant(1)%totalM
      do  i = 1,asf_initial%det_set%nod
         if (asf_initial%det_set%determinant(i)%totalM /= MM_i) then
            same_MM = .false.
            exit
         end if
      end do
      do  i = 1,asf_final%det_set%nod
         if (asf_final%det_set%determinant(i)%totalM /= MM_f) then
            same_MM = .false.
            exit
         end if
      end do
      !
      if (same_MM) then
         call reos_amplitude_same_MM()
      else
         call reos_amplitude_different_MM()
      end if
      !
   end subroutine reos_calculate_amplitudes
   !
   !
   subroutine reos_calculate_me(detf,deti,orbital_f,orbital_i,norb, &
                                detsym_f,detsym_i,MM_photon,trans)
   !--------------------------------------------------------------------
   ! Calculates the contributions of the one-particle multipole operators
   ! for all defined transitions for a given pair of Slater determinants.
   ! It first calculates the matrix of co-factors D(K/L); then the
   ! individual matrix elements are obtained from
   !
   !         Sum_{k,l}  <k|op|l> * d(k;l) * (-1)**(k+l)
   !
   ! where <k|op|l> is the corresponding one-particle matrix element.
   ! Here, k denotes a final and l an initial state orbital, i.e. in
   ! general <k|op|l>  =/=  <l|op|k>. the d(k/l) are the co-determinants
   ! to the overlap-matrix of the Slater determinants.
   !
   ! Calls: angular_momentum_j(), cofactor_1_of_overlap_matrix(),
   !        Wigner_Eckardt_geometry().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                     :: detf, deti, norb, MM_photon, &
                                                 trans
      integer, dimension(:), intent(in)       :: detsym_f, detsym_i
      type(nkappam), dimension(:), intent(in) :: orbital_f, orbital_i
      !
      integer       :: diff_symmetry, i, j, ki, kf, L, m, rf, ri,   &
                       j_f, kappa_f, mm_f, pqn_f, j_i, kappa_i, mm_i, pqn_i, &
                       i_lower, i_upper
      real(kind=dp) :: gf, me, minor
      logical       :: need_calculation
      logical, dimension(norb,norb)       :: cofactors_1_mask
      real(kind=dp), dimension(norb,norb) :: overlap, cofactors_1
      integer, dimension(number_of_symmetry_blocks) :: detsym_ir, detsym_fr
      !
      ! Determine the logical 'mask' which of the co-factors need to be
      ! calculated and set up the overlap matrix
      cofactors_1_mask = .true.
      need_calculation = .false.
      do  kf = 1,norb
         do  ki = 1,norb
            kappa_f = orbital_f(kf)%kappa;   mm_f  = orbital_f(kf)%mm
            pqn_f   = orbital_f(kf)%n
            kappa_i = orbital_i(ki)%kappa;   mm_i  = orbital_i(ki)%mm
            pqn_i   = orbital_i(ki)%n
            detsym_fr(:) = detsym_f(:)
            detsym_ir(:) = detsym_i(:)
            detsym_fr(sym_block(kappa_f,mm_f)) =                             &
                                        detsym_fr(sym_block(kappa_f,mm_f)) - 1
            detsym_ir(sym_block(kappa_i,mm_i)) =                             &
                                        detsym_ir(sym_block(kappa_i,mm_i)) - 1
            detsym_ir(:) = detsym_fr(:) - detsym_ir(:)
            diff_symmetry = sum( abs(detsym_ir) )
            if (diff_symmetry > 0) then
               cofactors_1_mask(kf,ki) = .false.
            else
               need_calculation = .true.
            end if
            if (kappa_f == kappa_i   .and.   mm_f == mm_i) then
               overlap(kf,ki) = reos_overlap(kappa_f,pqn_f,pqn_i)
            else
               overlap(kf,ki) = zero
            end if
         end do
      end do
      !
      if (.not.need_calculation) then
         ! stop "reos_calculate_me(): program stop A."
      end if
      !
      call cofactor_1_of_overlap_matrix(norb,overlap,cofactors_1, &
                                                     cofactors_1_mask)
      !
      ! Cycle over all transitions and multipole lines;
      ! in an inner loop then cycle over the individual co-factors
      if (trans == 0) then
         i_lower = 1;  i_upper = number_of_transitions
      else
         i_lower = trans;  i_upper = trans
      end if
      !
      do  i = i_lower,i_upper
         !x print *, "line i = ",i
         do  j = 1,transition(i)%number_of_mlines
            select case(transition(i)%mline(j)%multipole)
            case("E1", "M1");   L = 1
            case("E2", "M2");   L = 2
            case("E3", "M3");   L = 3
            case("E4", "M4");   L = 4
            case("E5", "M5");   L = 5
            case default; stop "reos_calculate_me(): program stop B."
            end select
            !
            me = zero
            do  kf = 1,norb
               do  ki = 1,norb
                  minor = cofactors_1(kf,ki)
                  if (abs(minor) < eps20) then
                     cycle
                  else if (mod(kf+ki,2) == 1) then
                     minor = - minor
                  end if
                  kappa_f = orbital_f(kf)%kappa;   mm_f  = orbital_f(kf)%mm
                  j_f     = angular_momentum_j(kappa_f)
                  kappa_i = orbital_i(ki)%kappa;   mm_i  = orbital_i(ki)%mm
                  j_i     = angular_momentum_j(kappa_i)
                  !
                  ! Calculate the 'geometrical factor' to the complete ME
                  gf = Wigner_Eckardt_geometry(j_f,mm_f,L+L,MM_photon,j_i,mm_i)
                  if (abs(gf) > eps10) then
                     !
                     ! Select orbital functions
                     rf = -1;   ri = -1
                     do m = 1,wave_final%number_of_rwf
                        if (orbital_f(kf)%n == wave_final%rwf(m)%orbital%n     &
                              .and. kappa_f == wave_final%rwf(m)%orbital%kappa &
                              ) then
                           rf = m; exit
                        end if
                     end do
                     do m = 1,wave_initial%number_of_rwf
                        if (orbital_i(ki)%n == wave_initial%rwf(m)%orbital%n   &
                            .and. kappa_i == wave_initial%rwf(m)%orbital%kappa &
                            ) then
                           ri = m; exit
                        end if
                     end do
                     if (rabs_use_stop       .and.  &
                        (rf == -1   .or.   ri == -1)) then
                        stop "reos_calculate_me(): program stop C."
                     end if
                     if (reos_precalculate_integrals) then
                        me = me + minor * gf *                                &
                             transition(i)%mline(j)%radial_int(               &
                               orbital_f(kf)%n,kappa_f,orbital_i(ki)%n,kappa_i)
                     else
                        me = me + minor * gf * reos_reduced_M_integral(i,     &
                                          transition(i)%mline(j)%multipole,   &
                                          transition(i)%mline(j)%gauge,       &
                                       wave_final%rwf(rf),wave_initial%rwf(ri))
                     end if
                  end if
               end do
            end do
            !
            ! print *, "i, j, gf, me = ",i, j, gf, me
            transition(i)%mline(j)%amplitude =                           &
               transition(i)%mline(j)%amplitude +                        &
               asf_final%asf(transition(i)%asff)%eigenvector(detf) *     &
               asf_initial%asf(transition(i)%asfi)%eigenvector(deti) * me
         end do
      end do
      !
   end subroutine reos_calculate_me
   !
   !
   subroutine reos_calculate_me_orth(detf,deti,orbital_f,orbital_i,norb,     &
                                                              MM_photon,trans)
   !--------------------------------------------------------------------
   ! Calculates the contributions of the one-particle multipole operators
   ! for all defined transitions for a given pair of Slater determinants.
   ! It assumes that all orbitals are orthogonal to each other.
   ! Therefore, the calculation of cofactors can be avoided; they vanish
   ! identically if the non-active electrons do not occupy the same set
   ! of orbitals. The calculation of matrix elements is then the same
   ! like in the non-orthogonal approach
   !
   !         Sum_{k,l}  <k|op|l> * d(k;l) * (-1)**(k+l)
   !
   ! where <k|op|l> is the corresponding one-particle matrix element.
   ! Here, k denotes a final and l an initial state orbital, i.e. in
   ! general <k|op|l>  =/=  <l|op|k>. the d(k/l) are the co-determinants
   ! to the overlap-matrix of the Slater determinants.
   !
   ! Calls: angular_momentum_j(), reos_return_M_integral(),
   !        Wigner_Eckardt_geometry().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                     :: detf, deti, norb, MM_photon, &
                                                 trans
      type(nkappam), dimension(:), intent(in) :: orbital_f, orbital_i
      !
      integer       :: diff_occ, i, j, ki, kf, L, m, phase, p, rf, ri,   &
                       j_f, kappa_f, mm_f, j_i, kappa_i, mm_i,           &
                       i_lower, i_upper
      real(kind=dp) :: gf, minor, me
      real(kind=dp), dimension(norb,norb) :: cofactors_1
      integer, dimension(norb)            :: occupation_i, occupation_f
      !
      cofactors_1 = zero
      !
      ! The contribution from this pair of determinants vanishs identically
      ! if occupations differ by more than one
      occupation_f = 1;   occupation_i = 0
      do  i = 1,norb
         do  j = 1,norb
            if (orbital_f(i)%n     == orbital_i(j)%n       .and.  &
                orbital_f(i)%kappa == orbital_i(j)%kappa   .and.  &
                orbital_f(i)%mm    == orbital_i(j)%mm  ) then
               occupation_i(i) = 1
               exit
            end if
         end do
      end do
      diff_occ = sum(occupation_f - occupation_i)
      !
      select case(diff_occ)
      case(:-1)
         stop "reos_calculate_me_orth(): program stop A."
      case(2:)
         return   ! Occupation numbers differ by 2 or more
      case(1)
         ! Determine kf, ki, and phase
         occupation_f = 0;   occupation_i = 1
         iloop: do  i = 1,norb
            do  j = 1,norb
               if (orbital_f(i)%n     == orbital_i(j)%n       .and.  &
                   orbital_f(i)%kappa == orbital_i(j)%kappa   .and.  &
                   orbital_f(i)%mm    == orbital_i(j)%mm  ) then
                  occupation_i(j) = 0
                  cycle iloop
               end if
            end do
            occupation_f(i) = 1
         end do iloop
         if (rabs_use_stop            .and.                         &
            (sum(occupation_f) /= 1   .or.   sum(occupation_i) /= 1)) then
            stop "reos_calculate_me_orth(): program stop B."
         end if
         do  i = 1,norb
            if (occupation_f(i) == 1)  kf = i
            if (occupation_i(i) == 1)  ki = i
         end do
         phase = kf - ki
         cofactors_1(kf,ki) = (-1)**phase
         j = 0
         do  i = 1,norb
            if (i == kf) cycle
          1 j = j + 1
            if (j == ki) goto 1
            if (orbital_f(i)%n     /= orbital_i(j)%n        .or.  &
                orbital_f(i)%kappa /= orbital_i(j)%kappa    .or.  &
                orbital_f(i)%mm    /= orbital_i(j)%mm  ) then
               print *, "kf, ki, i, j, nf, ni, kappaf, kappi, mmf, mmi = ",&
                         kf, ki, i, j, orbital_f(i)%n,orbital_i(j)%n,      &
                         orbital_f(i)%kappa,orbital_i(j)%kappa,            &
                         orbital_f(i)%mm,orbital_i(j)%mm
               print *, "(p,orbital_f(p)%n,orbital_f(p)%kappa,"//         &
                                          "orbital_f(p)%mm,p=1,norb) = ", &
                         (p,orbital_f(p)%n,orbital_f(p)%kappa,            &
                                           orbital_f(p)%mm,p=1,norb)
               print *, "(p,orbital_i(p)%n,orbital_i(p)%kappa,"//         &
                                          "orbital_i(p)%mm,p=1,norb) = ", &
                         (p,orbital_i(p)%n,orbital_i(p)%kappa,            &
                                           orbital_i(p)%mm,p=1,norb)
               stop "reos_calculate_me_orth(): program stop C."
            end if
         end do
      case(0)
         ! The final- and initial-state orbitals are assumed to be in the same
         ! order; terminate the program if not
         if (rabs_use_stop) then
         do  i = 1,norb
            if (orbital_f(i)%n     /= orbital_i(i)%n        .or.  &
                orbital_f(i)%kappa /= orbital_i(i)%kappa    .or.  &
                orbital_f(i)%mm    /= orbital_i(i)%mm  ) then
               stop "reos_calculate_me_orth(): program stop D."
            end if
         end do
         end if
         !
         do  i = 1,norb
            cofactors_1(i,i) = one
         end do
      case default
         stop "reos_calculate_me_orth(): program stop E."
      end select
      !
      ! Cycle over all transitions and multipole lines;
      ! in an inner loop then cycle over the individual co-factors
      if (trans == 0) then
         i_lower = 1;  i_upper = number_of_transitions
      else
         i_lower = trans;  i_upper = trans
      end if
      !
      do  i = i_lower,i_upper
         do  j = 1,transition(i)%number_of_mlines
            select case(transition(i)%mline(j)%multipole)
            case("E1", "M1");   L = 1
            case("E2", "M2");   L = 2
            case("E3", "M3");   L = 3
            case("E4", "M4");   L = 4
            case("E5", "M5");   L = 5
            case default; stop "reos_calculate_me_orth(): program stop F."
            end select
            !
            me = zero
            do  kf = 1,norb
               do  ki = 1,norb
                  minor = cofactors_1(kf,ki)
                  if (abs(minor) < eps20) then
                     cycle
                  else if (mod(kf+ki,2) == 1) then
                     minor = - minor
                  end if
                  kappa_f = orbital_f(kf)%kappa;   mm_f  = orbital_f(kf)%mm
                  j_f     = angular_momentum_j(kappa_f)
                  kappa_i = orbital_i(ki)%kappa;   mm_i  = orbital_i(ki)%mm
                  j_i     = angular_momentum_j(kappa_i)
                  !
                  ! Calculate the 'geometrical factor' to the complete ME
                  gf = Wigner_Eckardt_geometry(j_f,mm_f,L+L,MM_photon,j_i,mm_i)
                  if (abs(gf) > eps10) then
                     !
                     ! Select orbital functions
                     rf = -1;   ri = -1
                     do m = 1,wave_final%number_of_rwf
                        if (orbital_f(kf)%n == wave_final%rwf(m)%orbital%n     &
                              .and. kappa_f == wave_final%rwf(m)%orbital%kappa &
                              ) then
                           rf = m; exit
                        end if
                     end do
                     do m = 1,wave_initial%number_of_rwf
                        if (orbital_i(ki)%n == wave_initial%rwf(m)%orbital%n   &
                            .and. kappa_i == wave_initial%rwf(m)%orbital%kappa &
                            ) then
                           ri = m; exit
                        end if
                     end do
                     if (rabs_use_stop      .and.   &
                        (rf == -1   .or.   ri == -1)) then
                        stop "reos_calculate_me_orth(): program stop G."
                     end if
                     if (reos_precalculate_integrals) then
                        me = me + minor * gf *                                &
                             transition(i)%mline(j)%radial_int(               &
                               orbital_f(kf)%n,kappa_f,orbital_i(ki)%n,kappa_i)
                     else
                        me = me + minor * gf * reos_reduced_M_integral(i,     &
                                          transition(i)%mline(j)%multipole,   &
                                          transition(i)%mline(j)%gauge,       &
                                       wave_final%rwf(rf),wave_initial%rwf(ri))
                     end if
                  end if
               end do
            end do
            !
            transition(i)%mline(j)%amplitude =                           &
               transition(i)%mline(j)%amplitude +                        &
               asf_final%asf(transition(i)%asff)%eigenvector(detf)   *   &
               asf_initial%asf(transition(i)%asfi)%eigenvector(deti) * me
         end do
      end do
      !
   end subroutine reos_calculate_me_orth
   !
   !
   subroutine reos_collect_input()
   !--------------------------------------------------------------------
   ! Collect and proceeds all input for the REOS program.
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
      print *, "Restart a previously aborted calculation ?"
      yes = get_yes_stream()
      if (yes) then
         !x reos_apply_restart = .true.
         call reos_open_restart(.false.)
         call reos_readwrite_dump(31,.true.)
         !x print *, "after readwrite_dump(31,.true.)"
         call reos_readwrite_transitions(32,.true.)
         !x print *, "after readwrite_transitions(32,.true.)"
         reos_apply_restart = .true.
	 return
      end if
      !
      ! Open, check, load data from, and close the  .iso  file
      call file_get_isodat_grasp2k()
      !
      ! Determine the transition multipoles
      call input_transition_multipoles(number_of_multipoles,reos_multipole)
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
      ! Adopt previously calculated probabilities and lifetimes to
      ! experimental transition energies
      !
      ! Select individual pairs of transitions
      call input_transition_pairs(number_of_transitions)
      !
      print *, "Assume orthogonality between the orbital sets of the"// &
               " initial and final atomic states ?"
      print *, " This feature accelerates the computations but also "// &
               " neglects all relaxation effects"
      print *, " on the transitions probabilities."
      yes = get_yes_stream()
      if (yes) reos_assume_orthogonality = .true.
      !
      print *, "Calculate all radial integrals only once at the beginning"// &
               " of the computation ?"
      print *, " This feature accelerates the computations but requires"// &
               " considerable more memory, "
      print *, " in particular, if many transitions need to be calculated."
      yes = get_yes_stream()
      if (yes) reos_precalculate_integrals = .true.
      !
      print *, "Enable the restart of the calculation if the run cannot"// &
               " be completed successfully ?"
      print *, " This feature dumps all necessary information to a restart"// &
               " file from which the calculation can be continued. "
      yes = get_yes_stream()
      if (yes) then
         reos_enable_restart = .true.
         !
         ! Open the  .res-dump and .res-tran  files
         call reos_open_restart(.true.)
      end if
      !
      print *, "Sort transitions in ascending order of energy ?"
      yes = get_yes_stream()
      if (yes) reos_sort_transition_energy = .true.
      !
      print *, "Read in and apply experimental energies for the calculation"//&
               " of transition probabilities ?"
      yes = get_yes_stream()
      if (yes) reos_apply_exp_energies = .true.
      !
      print *, "Einstein A and B coefficients are printed in SI units;"
      print *, " use Hartree atomic units instead ?"
      yes = get_yes_stream()
      if (yes) reos_print_AB_in_hartree = .true.
      !
      print *, "Print all selected transitions and their energies"// &
               " before the computation starts (this is the default) ?"
      yes = get_yes_stream()
      if (.not.yes) reos_print_selected_trans = .false.
      !
      print *, "Write out the transition energies and amplitudes to an"// &
               " .trn file for further data processing,"
      print *, " for instance, to adopt them to experimental transition"//&
               " energies ?"
      yes = get_yes_stream()
      if (yes) then
         reos_write_transition_file = .true.
	 !
         call file_open_formatted_stream(26, &
	    "Enter a file name for the  reos.trn  file:")
         write(26, "(a)") "REOS Transition energy and amplitude file"
         write(26,*)
      end if
      !
      ! Fix a cutoff criterium to neglect small admixtures
      print *, "Fix a cutoff criterium other than '10e-8' to neglect small"
      print *, " admixtures to the computation of transition probabilities ?"
      yes = get_yes_stream()
      if (yes) then
         print *, "All pairs of determinants with "// &
                  " 0 < |d_r(alpha)*d_s(beta)| < cutoff "
         print *, " will be neglected; enter cutoff: "
         read *,  reos_determinant_cutoff
    5    print *, "Enter an integer counter >= 1000 to provide a short"// &
                  " printout about the computational process: "
         read *,  reos_determinant_counter
         if (reos_determinant_counter < 1000) then
            print *, "Expect  >= 1000; reenter ..."
            goto 5
         end if
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
   end subroutine reos_collect_input
   !
   !
   subroutine reos_convert_probability(totalJ_i,totalJ_f,L,energy_au, &
		   amplitude,einstein_A,einstein_B,oscillator,decay_width)
   !--------------------------------------------------------------------
   ! Returns the basic oscillator strength, Einstein A and B coefficients,
   ! and the decay width for a given transition amplitude.
   !--------------------------------------------------------------------
      !
      integer, intent(in)        :: totalJ_i, totalJ_f, L
      real(kind=dp), intent(in)  :: energy_au, amplitude
      real(kind=dp), intent(out) :: einstein_A, einstein_B, oscillator, &
                                    decay_width
      !
      real(kind=dp) :: stat_factor, omega_over_c, fosc, faau, fbau
      !
      ! Evaluate statistical factors and constants
      stat_factor  = (totalJ_f + one) / (L + L + one)
      omega_over_c = energy_au / c
      fosc         = c * stat_factor / omega_over_c
      faau         = two * omega_over_c * omega_over_c * fosc / &
                     (c * (totalJ_i + one))
      fbau         = pi * fosc / energy_au
      !
      ! Calculate Einstein A and B coefficients and oscillator strengths
      oscillator  = fosc * amplitude * amplitude
      einstein_A  = faau * amplitude * amplitude
      einstein_B  = fbau * amplitude * amplitude
      decay_width = einstein_A * convert_au_to_ev
      if (.not.reos_print_AB_in_hartree) then
         einstein_A = einstein_A  * convert_einstein_a_to_si
         einstein_B = einstein_B  * convert_einstein_b_to_si
      end if
      !
   end subroutine reos_convert_probability
   !
   !
   subroutine reos_initialize_rwf_storage()
   !--------------------------------------------------------------------
   ! Initializes the arrays of type(grasp2k_orbital) for the storage of
   ! the radial wave functions.
   !--------------------------------------------------------------------
      !
      integer :: i, iorb, kappa, n
      !
      ! Initialize storage for initial-state wave functions
      !
      n = 0;   kappa = 0;   wave_initial%number_of_rwf = 0
      do  iorb = 1,asf_initial%det_set%norbital
         if (asf_initial%det_set%orbital(iorb)%n     /= n     .or.   &
             asf_initial%det_set%orbital(iorb)%kappa /= kappa) then
            wave_initial%number_of_rwf = wave_initial%number_of_rwf + 1
            n     = asf_initial%det_set%orbital(iorb)%n
            kappa = asf_initial%det_set%orbital(iorb)%kappa
         end if
      end do
      !
      allocate( wave_initial%rwf(1:wave_initial%number_of_rwf) )
      n = 0;   kappa = 0;   i = 0
      do  iorb = 1,asf_initial%det_set%norbital
         if (asf_initial%det_set%orbital(iorb)%n     /= n     .or.   &
             asf_initial%det_set%orbital(iorb)%kappa /= kappa) then
            i     = i + 1
            n     = asf_initial%det_set%orbital(iorb)%n
            kappa = asf_initial%det_set%orbital(iorb)%kappa
            wave_initial%rwf(i)%orbital%n     =  &
                                asf_initial%det_set%orbital(iorb)%n
            wave_initial%rwf(i)%orbital%kappa =  &
                                asf_initial%det_set%orbital(iorb)%kappa
            wave_initial%rwf(i)%mtp    = 0
            wave_initial%rwf(i)%energy = zero
            wave_initial%rwf(i)%gamma  = zero
            wave_initial%rwf(i)%pz     = zero
         end if
      end do
      !
      if (rabs_use_stop   .and.  i /= wave_initial%number_of_rwf) then
         stop "reos_initialize_rwf_storage(): program stop A."
      end if
      !
      ! Initialize storage for final-state wave functions
      !
      n = 0;   kappa = 0;   wave_final%number_of_rwf = 0
      do  iorb = 1,asf_final%det_set%norbital
         if (asf_final%det_set%orbital(iorb)%n     /= n     .or.   &
             asf_final%det_set%orbital(iorb)%kappa /= kappa) then
            wave_final%number_of_rwf = wave_final%number_of_rwf + 1
            n     = asf_final%det_set%orbital(iorb)%n
            kappa = asf_final%det_set%orbital(iorb)%kappa
         end if
      end do
      !
      allocate( wave_final%rwf(1:wave_final%number_of_rwf) )
      n = 0;   kappa = 0;   i = 0
      do  iorb = 1,asf_final%det_set%norbital
         if (asf_final%det_set%orbital(iorb)%n     /= n     .or.   &
             asf_final%det_set%orbital(iorb)%kappa /= kappa) then
            i     = i + 1
            n     = asf_final%det_set%orbital(iorb)%n
            kappa = asf_final%det_set%orbital(iorb)%kappa
            wave_final%rwf(i)%orbital%n     = &
                              asf_final%det_set%orbital(iorb)%n
            wave_final%rwf(i)%orbital%kappa = &
                              asf_final%det_set%orbital(iorb)%kappa
            wave_final%rwf(i)%mtp    = 0
            wave_final%rwf(i)%energy = zero
            wave_final%rwf(i)%gamma  = zero
            wave_final%rwf(i)%pz     = zero
         end if
      end do
      !
      if (rabs_use_stop   .and.   i /= wave_final%number_of_rwf) then
         stop "reos_initialize_rwf_storage(): program stop B."
      end if
      !
   end subroutine reos_initialize_rwf_storage
   !
   !
   subroutine reos_load_trn_file()
   !--------------------------------------------------------------------
   ! Opens a .trn transition energy and amplitude file on stream 26 and
   ! reads in the information for all individual transitions.
   !
   ! Calls: file_open().
   !--------------------------------------------------------------------
      !
      integer            :: i, j, ierr
      character(len=15)  :: record
      character(len=256) :: reos_trn_file
      !
    1 print *, "Enter a file name for an (existing)  reos.trn  file:"
      read *,  reos_trn_file
      call file_open(21,reos_trn_file,"formatted  ","old",ierr)
      !
      if (ierr /= 0) goto 1
      !
      ! Check the header of the file
      read(21, "(a)") record(1:15)
      if (record(1:15) /= "REOS Transition") then
         print *, "Not a .trn reos transition energy and amplitude file;"
         goto 1
      end if
      read(21,*)
      !
      read(21, "(i6)") number_of_transitions
      allocate( transition(1:number_of_transitions) )
      !
      read(21,*)
      do  i = 1,number_of_transitions
         read(21,2) transition(i)%asfi,     transition(i)%asff,        &
                    transition(i)%level_i,  transition(i)%level_f,     &
                    transition(i)%totalJ_i, transition(i)%parity_i,    &
                    transition(i)%totalJ_f, transition(i)%parity_f,    &
                    transition(i)%number_of_mlines, transition(i)%energy
         allocate( transition(i)%mline(1:transition(i)%number_of_mlines) )
         read(21,3) (transition(i)%mline(j)%multipole,                 &
                     transition(i)%mline(j)%gauge,                     &
                     transition(i)%mline(j)%amplitude,                 &
                     j=1,transition(i)%number_of_mlines)
      end do
    2 format(2i4,2x,i4,2x,i4,2x,2(i4,a1,1x),2x,i2,1pe14.7)
    3 format(3x,20(a2,1x,a9,1x,1pe14.7,3x))
      !
   end subroutine reos_load_trn_file
   !
   !
   subroutine reos_open_restart(first_time)
   !--------------------------------------------------------------------
   ! Opens two restart file for the REOS program and writes out a
   ! header to these files;  .res-dump  contains a compact representation
   ! of the determinant basis, mixing coefficients, and radial integrals,
   ! and .res-tran contains a list of all transitions and transition
   ! amplitudes from which the calculation starts again.
   !
   ! Calls: file_open().
   !--------------------------------------------------------------------
      !
      logical, intent(in) :: first_time
      !
      integer :: ierr
      character(len=256) :: reos_res_file, reos_res_dump_file, &
                            reos_res_tran_file
      !
      if (first_time) then
       1 print *, "Enter a file name for the  reos.res  file(s):"
         read *,  reos_res_file
         reos_res_dump_file = trim(reos_res_file)//"-dump"
         reos_res_tran_file = trim(reos_res_file)//"-tran"
         call file_open(31,reos_res_dump_file,"formatted  ","new",ierr)
         if (ierr /= 0) goto 1
         call file_open(32,reos_res_tran_file,"formatted  ","new",ierr)
         if (ierr /= 0) goto 1
      else
       2 print *, "Enter a file name for the  reos.res  file(s):"
         read *,  reos_res_file
         reos_res_dump_file = trim(reos_res_file)//"-dump"
         reos_res_tran_file = trim(reos_res_file)//"-tran"
         call file_open(31,reos_res_dump_file,"formatted  ","old",ierr)
         if (ierr /= 0) goto 2
         call file_open(32,reos_res_tran_file,"formatted  ","old",ierr)
         if (ierr /= 0) goto 2
      end if
      !
   end subroutine reos_open_restart
   !
   !
   subroutine reos_print_results(stream)
   !--------------------------------------------------------------------
   ! Writes the transition probabilities and lifetimes in a neat summary
   ! to the .sum file.
   !
   ! Calls: angular_momentum_string(), reos_convert_probability().
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: stream
      !
      integer :: i, j, istate_low, istate_up, rank
      real(kind=dp) :: einstein_A, einstein_B, oscillator, decay_width,    &
                       energy_au, energy, tb, tb_au, tb_cm, tb_ev, tb_sec, &
                       tc, tc_au, tc_cm, tc_ev, tc_sec, tb_inv, tc_inv
      real(kind=dp), dimension(:), allocatable :: total_babushkin, &
                                                  total_coulomb
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
      if (energy_unit == "A      "  .and.  reos_print_AB_in_hartree) then
         write(stream,2)
         write(stream,4)
      else if (energy_unit == "A      ") then
         write(stream,3)
         write(stream,5)
      else if (reos_print_AB_in_hartree) then
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
            case default; stop "reos_print_results(): program stop A."
            end select
            call reos_convert_probability(transition(i)%totalJ_i,        &
                                          transition(i)%totalJ_f,        &
                                          rank,energy_au,                &
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
               ! stop "reos_print_results(): program stop B."
            end if
         end do
      end do
    9 format(2x,i3," -",i3,3x,a4,1x,a1,4x,a4,1x,a1,3x,1pd12.5,3x,a2,     &
             4x,a9,3x,1pd12.5,3x,1pd12.5,3x,1pd12.5,3x,1pd12.5)
             !
      write(stream,8)
      !
      ! Print lifetimes and width of levels
      write(stream,*)
      write(stream,10)
   10 format(//"Radiative lifetimes, total rates, and widths:"                 &
              /"---------------------------------------------"                 &
            ///" LeveL",6x,"Gauge",15x,"Lifetime",11x,"Total rate",32x,"Width" &
              /" -----",6x,"-----",15x,"--------",11x,"----------",8x,         &
               "-----------------------------------------------------",        &
          /32x,"Seconds",14x,"1/sec",14x,"Hartrees",12x,"Kaysers",16x,"eV"/)
      !
      do  i = istate_low,istate_up
         if( total_babushkin(i) == zero  .and.  total_coulomb(i) == zero) then
            cycle
         end if
         tb = total_Babushkin(i)
         tc = total_Coulomb(i)
         if (.not.reos_print_AB_in_hartree) then
            tb_cm  = tb / c_vacuum_in_cm_per_s
            tc_cm  = tc / c_vacuum_in_cm_per_s
            tc_au  = tc_cm / convert_au_to_kaysers
            tb_au  = tb_cm / convert_au_to_kaysers
            tc_sec = one / tc
            tb_sec = one / tb
            tc_ev  = tc_au * convert_au_to_ev
            tb_ev  = tb_au * convert_au_to_ev
	    tb_inv = one / tb_sec
	    tc_inv = one / tc_sec
         else
            tc_au  = tc
            tb_au  = tb
            tc_cm  = tc * convert_au_to_kaysers
            tb_cm  = tb * convert_au_to_kaysers
            tc_sec = one / ( tc_cm * c_vacuum_in_cm_per_s )
            tb_sec = one / ( tb_cm * c_vacuum_in_cm_per_s )
            tc_ev  = tc_au * convert_au_to_ev
            tb_ev  = tb_au * convert_au_to_ev
	    tb_inv = one / tb_sec
	    tc_inv = one / tc_sec
         end if
         !
         if (pure_magnetic(i)) then
            write(stream,11) i,tc_sec,tc_inv,tc_au,tc_cm,tc_ev
         else
            write(stream,12) i,tb_sec,tb_inv,tb_au,tb_cm,tb_ev
            write(stream,13)   tc_sec,tc_inv,tc_au,tc_cm,tc_ev
         end if
      11 format(1x,i4,6x,"Magnetic:  ",5(1pd20.7)/)
      12 format(1x,i4,6x,"Babushkin: ",5(1pd20.7))
      13 format(     11x,"Coulomb:   ",5(1pd20.7)/)
         !
      end do
      deallocate( total_babushkin, total_coulomb, pure_magnetic )
      !
   end subroutine reos_print_results
   !
   !
   subroutine reos_print_summary()
   !--------------------------------------------------------------------
   ! Appends a summary of the input data from the initial and final state
   ! .xpn files to the .sum file.
   !
   ! Calls: angular_momentum_string(), get_month(), orbital_symmetry().
   !--------------------------------------------------------------------
      !
      integer :: i, totalJ
      character(len=3)  :: month
      character(len=8)  :: cdate
      character(len=10) :: ctime
      !
      ! Get the date and time of day; make this information the header of
      ! the  reos.sum  summary file
      call date_and_time(date=cdate,time=ctime)
      month = get_month(cdate(5:6))
      write(24,*) "REOS run at "//ctime(1:2)//":"//ctime(3:4)//":"//	  &
		  ctime(5:6)//" on "//month//				    &
		  " "//cdate(7:8)//" "//cdate(1:4)//"."
      !
      ! Write out the basic dimensions of the initial and final-state
      ! electron clouds
      write(24,*)
      write(24,*) "There are ",asf_initial%det_set%number_of_electrons,     &
		  " electrons in the cloud of the initial ion"
      write(24,*) " in ",asf_initial%det_set%nod," determinants based on ", &
		  asf_initial%det_set%norbital, 			    &
		  " (n,kappa,m)-orbital functions ..."
      !
      write(24,*)
      write(24,*) "... and ",asf_final%det_set%number_of_electrons,	    &
		  " electrons in the cloud of the final ion"
      write(24,*) " in ",asf_final%det_set%nod," determinants based on ",   &
		  asf_final%det_set%norbital," (n,kappa,m)-orbital functions."
      !
      ! Write out the nuclear parameters and the speed of light
      write(24,*)
      write(24,*) "The atomic number is Z =",nuclear_charge,	  &
		  " distributed in a "//trim(nuclear_model)//	  &
		  "-like nucleus."
      write(24,*) "Speed of light =",c," atomic units."
      !
      ! Write out the parameters of the radial grid
      write(24,*)
      if (hp_grasp2k == zero) then
	 write(24,1) rnt_grasp2k,h_grasp2k,n_grasp2k
      else
	 write(24,2) rnt_grasp2k,h_grasp2k,hp_grasp2k,n_grasp2k
      endif
      write(24,3) r_grasp2k(1),r_grasp2k(2),r_grasp2k(n_grasp2k)
      !
    1 format(" Radial grid: R(I) = RNT*(exp((I-1)*H)-1), I = 1, ..., N;", &
             /" --------------------------------------------------------", &
            //"  RNT  = ",1p,d19.12," Bohr radii;",                       &
             /"  H    = ",   d19.12," Bohr radii;" ,                      &
             /"  N    = ",1i4,";")
    2 format(" Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,",       &
             " I = 1, ..., N;",                                           &
             /" --------------------------------------------------",      &
             /"---------------",                                          &
            //"  RNT  = ",1p,d19.12," Bohr radii;",                       &
             /"  H    = ",   d19.12," Bohr radii;",                       &
             /"  HP   = ",   d19.12," Bohr radii;",                       &
             /"  N    = ",1i4,";")
    3 format( "  R(1) = ",1p,1d19.12," Bohr radii;",                      &
             /"  R(2) = ",   1d19.12," Bohr radii;",                      &
             /"  R(N) = ",   1d19.12," Bohr radii.")
      !
      ! Write out the orbital properties
      write(24,*)
      write(24,*) "Initial-state subshell radial wavefunction summary:"
      write(24,*) "---------------------------------------------------"
      write(24,*)
      write(24,4)
      write(24,*)
      do  i = 1,wave_initial%number_of_rwf
	 write(24,5) wave_initial%rwf(i)%orbital%n,			  &
		     orbital_symmetry(wave_initial%rwf(i)%orbital%kappa), &
		     wave_initial%rwf(i)%energy,wave_initial%rwf(i)%pz,   &
		     wave_initial%rwf(i)%gamma, 			  &
		     wave_initial%rwf(i)%P(2),wave_initial%rwf(i)%Q(2),   &
		     wave_initial%rwf(i)%mtp
      end do
    4 format(" Subshell",11x,"E",20x,"p0",18x,                            &
             "gamma",19x,"P(2)",18x,"Q(2)",10x,"mtp")
    5 format(3x,1i2,1a2,1x,1p,5(3x,1e19.12),3x,1i3)
      !
      write(24,*)
      write(24,*) "Final-state subshell radial wavefunction summary:"
      write(24,*) "-------------------------------------------------"
      write(24,*)
      write(24,4)
      write(24,*)
      do  i = 1,wave_final%number_of_rwf
	 write(24,5) wave_final%rwf(i)%orbital%n,			  &
		     orbital_symmetry(wave_final%rwf(i)%orbital%kappa),   &
		     wave_final%rwf(i)%energy,wave_final%rwf(i)%pz,	  &
		     wave_final%rwf(i)%gamma,				  &
		     wave_final%rwf(i)%P(2),wave_final%rwf(i)%Q(2),	  &
		     wave_final%rwf(i)%mtp
      end do
      !
      ! Write the list of eigenpair indices
      write(24,*)
      write(24,*) 'Initial-state ASF:'
      write(24,*) '------------------'
      write(24,*)
      write(24,6)
      write(24,*)
      do  i = 1,asf_initial%noasf
	 totalJ = asf_initial%asf(i)%totalJ
	 write(24,7) asf_initial%asf(i)%level_No,			  &
		     angular_momentum_string(totalJ,4), 		  &
		     asf_initial%asf(i)%parity,asf_initial%asf(i)%energy
      end do
    6 format(" Level_No",6x,"J^P",10x,"Energy (a.u.)")
    7 format(2x,i4,8x,a4,a1,6x,1e19.12)
      !
      write(24,*)
      write(24,*) 'Final-state ASF:'
      write(24,*) '----------------'
      write(24,*)
      write(24,6)
      write(24,*)
      do  i = 1,asf_final%noasf
	 totalJ = asf_final%asf(i)%totalJ
	 write(24,7) asf_final%asf(i)%level_No, 			  &
		     angular_momentum_string(totalJ,4), 		  &
		     asf_final%asf(i)%parity,asf_final%asf(i)%energy
      end do
      !
   end subroutine reos_print_summary
   !
   !
   subroutine reos_print_transitions(stream)
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
         write(stream,2) transition(i)%level_i,transition(i)%level_f,       &
                   trim(angular_momentum_string(transition(i)%totalJ_i,4)), &
                         transition(i)%parity_i,                            &
                   trim(angular_momentum_string(transition(i)%totalJ_f,4)), &
                         transition(i)%parity_f,energy,(mult(j),j=1,nmult)
      end do
      write(stream,3)
    1 format( "The following ",i5," transitions are selected:",              &
        //,"     I-level-F     I--J^P--F      Transition Energy       ",     &
           "Multipoles ",                                                    &
         /,"                                     (in ",a4,")  ",             &
         /,4x,66("-") )
    2 format(4x,i4," -",i4,3x,a4,a1,3x,a4,a1,5x,1pe14.7,9x,10(a2,1x))
    3 format(4x,66("-") )
      !
   end subroutine reos_print_transitions
   !
   !
   subroutine reos_readwrite_dump(stream,read_from)
   !--------------------------------------------------------------------
   ! Reads or writes all necessary information (apart from the transition()
   ! array) for restarting REOS from or to a file on stream.
   !--------------------------------------------------------------------
      !
      integer, intent(in)            :: stream
      logical, intent(in)            :: read_from
      !
      integer :: i1, i2, kappa, pqn_f, pqn_i, bl1, bu1, bl2, bu2, bl3, bu3
      logical :: need_detsym, need_overlap
      !
      call readwrite_asf_det_basis(stream,read_from,asf_initial)
      call readwrite_asf_det_basis(stream,read_from,asf_final)
      !
      call readwrite_grasp2k_grid(stream,read_from)
      call readwrite_grasp2k_orbital(stream,read_from,wave_initial)
      call readwrite_grasp2k_orbital(stream,read_from,wave_final)
      !
      if (read_from) then
         read(stream,*)  c
         read(stream,*)  reos_use_compact_xpn_file   ,  &
                         reos_apply_exp_energies     ,  &
                         reos_apply_restart          ,  &
                         reos_assume_orthogonality   ,  &
                         reos_enable_restart         ,  &
                         reos_precalculate_integrals ,  &
                         reos_sort_transition_energy ,  &
                         reos_print_AB_in_hartree    ,  &
                         reos_print_selected_trans   ,  &
                         reos_use_formatted_rwf_file ,  &
                         reos_write_transition_file  ,  &
                         energy_inverse
         read(stream,*)  number_of_multipoles        ,  &
                         number_of_transitions       ,  &
                         number_of_symmetry_blocks   ,  &
                         reos_determinant_counter    ,  &
                         reos_dump_period            ,  &
                         reos_dump_billion           ,  &
                         reos_dump_counter           ,  &
                         reos_determinant_cutoff
         read(stream,*)  energy_factor               ,  &
                         energy_unit
	 read(stream,*)  need_overlap, bl1, bu1, bl2, bu2, bl3, bu3
	 if (need_overlap) then
	    allocate( reos_overlap(bl1:bu1,bl2:bu2,bl3:bu3) )
            ! Calculate the non-orthogonal and overlap integrals
            do  i1 = 1,wave_final%number_of_rwf
               do  i2 = 1,wave_initial%number_of_rwf
                  if (wave_final%rwf(i1)%orbital%kappa == &
                      wave_initial%rwf(i2)%orbital%kappa) then
                     kappa = wave_final%rwf(i1)%orbital%kappa
                     pqn_f = wave_final%rwf(i1)%orbital%n
                     pqn_i = wave_initial%rwf(i2)%orbital%n
                     reos_overlap(kappa,pqn_f,pqn_i) = &
                        rk_integral_grasp2k_ab(wave_final,wave_initial,0,i1,i2)
                  end if
               end do
            end do
	 end if
	 !
         read(stream,*)  sym_block
	 read(stream,*)  need_detsym
	 if (need_detsym) then
            call reos_set_determinant_sym()
         end if
         read(stream,*)  debug_reduced_M_integral
      else
         write(stream,*) c
         write(stream,*) reos_use_compact_xpn_file   ,  &
                         reos_apply_exp_energies     ,  &
                         reos_apply_restart          ,  &
                         reos_assume_orthogonality   ,  &
                         reos_enable_restart         ,  &
                         reos_precalculate_integrals ,  &
                         reos_sort_transition_energy ,  &
                         reos_print_AB_in_hartree    ,  &
                         reos_print_selected_trans   ,  &
                         reos_use_formatted_rwf_file ,  &
                         reos_write_transition_file  ,  &
                         energy_inverse
         write(stream,*) number_of_multipoles        ,  &
                         number_of_transitions       ,  &
                         number_of_symmetry_blocks   ,  &
                         reos_determinant_counter    ,  &
                         reos_dump_period            ,  &
                         reos_dump_billion           ,  &
                         reos_dump_counter           ,  &
                         reos_determinant_cutoff
         write(stream,*) energy_factor               ,  &
                         energy_unit
         if (allocated(reos_overlap)) then
            need_overlap = .true.
            bl1 = lbound(reos_overlap,1);   bu1 = ubound(reos_overlap,1)
            bl2 = lbound(reos_overlap,2);   bu2 = ubound(reos_overlap,2)
            bl3 = lbound(reos_overlap,3);   bu3 = ubound(reos_overlap,3)
	 else
	    need_overlap = .false.
	    bl1 = 0;   bu1 = 0;   bl2 = 0;   bu2 = 0
            bl3 = 0;   bu3 = 0
	 end if
	 write(stream,*) need_overlap, bl1, bu1, bl2, bu2, bl3, bu3
	 !
         write(stream,*) sym_block
         !
         if (allocated(reos_detsym_i)  .or.  allocated(reos_detsym_f)) then
            need_detsym = .true.
         else
            need_detsym = .false.
         end if
	 write(stream,*) need_detsym
         write(stream,*) debug_reduced_M_integral
      end if
      !
   end subroutine reos_readwrite_dump
   !
   !
   subroutine reos_readwrite_transitions(stream,read_from)
   !--------------------------------------------------------------------
   ! Reads or writes the array transitions()) of type(reos_transition)
   ! from or to a file on stream.
   !--------------------------------------------------------------------
      !
      integer, intent(in)            :: stream
      logical, intent(in)            :: read_from
      !
      integer :: i, i1, i2, j, int_counter, m, kappa_f, kappa_i, pqn_f, pqn_i,&
                 bl1, bu1, bl2, bu2, bl3, bu3, bl4, bu4
      logical :: needb0, needb1, needb2, needb3, needb4, needb5, &
                 need_radial_int
      real(kind=dp)                         :: arg
      real(kind=dp), dimension(1:n_grasp2k) :: bessel
      !
      if (read_from) then
         read(stream,*) reos_dump_billion, reos_dump_counter
         read(stream,*) number_of_transitions
	 allocate( transition(1:number_of_transitions) )
      else
         write(stream,*) reos_dump_billion, reos_dump_counter
         write(stream,*) number_of_transitions
      end if
      !
      int_counter = 0
      do  i = 1,number_of_transitions
         if (read_from) then
	    read(stream,*)  transition(i)%asfi,     transition(i)%asff,     &
	                    transition(i)%level_i,  transition(i)%level_f,  &
			    transition(i)%totalJ_i, transition(i)%totalJ_f, &
			    transition(i)%number_of_mlines,                 &
			    transition(i)%energy
            read(stream,*)  transition(i)%parity_i
	    read(stream,*)  transition(i)%parity_f
 	    allocate( transition(i)%mline(1:transition(i)%number_of_mlines) )
	    read(stream,*)  needb0, needb1, needb2, needb3, needb4, needb5
            arg = transition(i)%energy / c
            if (needb0) then
	       allocate( transition(i)%bessel0(1:n_grasp2k) )
               call reos_calculate_Bessel(0,arg,bessel)
               transition(i)%bessel0 = bessel
	    end if
            if (needb1) then
	       allocate( transition(i)%bessel1(1:n_grasp2k) )
               call reos_calculate_Bessel(1,arg,bessel)
               transition(i)%bessel1 = bessel
	    end if
            if (needb2) then
	       allocate( transition(i)%bessel2(1:n_grasp2k) )
               call reos_calculate_Bessel(2,arg,bessel)
               transition(i)%bessel2 = bessel
	    end if
            if (needb3) then
	       allocate( transition(i)%bessel3(1:n_grasp2k) )
               call reos_calculate_Bessel(3,arg,bessel)
               transition(i)%bessel3 = bessel
	    end if
            if (needb4) then
	       allocate( transition(i)%bessel4(1:n_grasp2k) )
               call reos_calculate_Bessel(4,arg,bessel)
               transition(i)%bessel4 = bessel
	    end if
            if (needb5) then
	       allocate( transition(i)%bessel5(1:n_grasp2k) )
               call reos_calculate_Bessel(5,arg,bessel)
               transition(i)%bessel5 = bessel
	    end if
	    !
	    do  j = 1,transition(i)%number_of_mlines
	       read(stream,*)  transition(i)%mline(j)%multipole
	       read(stream,*)  transition(i)%mline(j)%gauge
	       read(stream,*)  transition(i)%mline(j)%amplitude
	       read(stream,*)  need_radial_int, &
	                       bl1, bu1, bl2, bu2, bl3, bu3, bl4, bu4
	       if (need_radial_int) then
	          allocate( transition(i)%mline(j)%radial_int(bl1:bu1, &
		                              bl2:bu2,bl3:bu3,bl4:bu4) )
                  transition(i)%mline(j)%radial_int = zero
                  !
                  do  i1 = 1,wave_final%number_of_rwf
                     do  i2 = 1,wave_initial%number_of_rwf
                        pqn_f   = wave_final%rwf(i1)%orbital%n
                        kappa_f = wave_final%rwf(i1)%orbital%kappa
                        pqn_i   = wave_initial%rwf(i2)%orbital%n
                        kappa_i = wave_initial%rwf(i2)%orbital%kappa
                        transition(i)%mline(j)%radial_int(pqn_f,kappa_f,      &
                                                          pqn_i,kappa_i) =    &
                        reos_reduced_M_integral(i,                            &
			          transition(i)%mline(j)%multipole,           &
                                  transition(i)%mline(j)%gauge,               &
                                  wave_final%rwf(i1),wave_initial%rwf(i2))
                        if (abs(transition(i)%mline(j)%radial_int(            &
			                                     pqn_f,kappa_f,   &
                          pqn_i,kappa_i)) > zero) int_counter = int_counter + 1
                     end do
                  end do
	      end if
	    end do
	 else
	    write(stream,*) transition(i)%asfi,     transition(i)%asff,     &
	                    transition(i)%level_i,  transition(i)%level_f,  &
			    transition(i)%totalJ_i, transition(i)%totalJ_f, &
			    transition(i)%number_of_mlines,                 &
			    transition(i)%energy
            write(stream,*) transition(i)%parity_i
	    write(stream,*) transition(i)%parity_f
	    needb0 = .false.; needb1 = .false.; needb2 = .false.
	    needb3 = .false.; needb4 = .false.; needb5 = .false.
            do  m = 1,transition(i)%number_of_mlines
            select case(transition(i)%mline(m)%multipole)
            case("E1", "M1"); needb0 = .true.; needb1 = .true.; needb2 = .true.
            case("E2", "M2"); needb1 = .true.; needb2 = .true.; needb3 = .true.
            case("E3", "M3"); needb2 = .true.; needb3 = .true.; needb4 = .true.
            case("E4", "M4"); needb3 = .true.; needb4 = .true.; needb5 = .true.
            case default
               stop "reos_readwrite_transitions(): program stop A."
            end select
            end do
	    write(stream,*) needb0, needb1, needb2, needb3, needb4, needb5
	    !
	    do  j = 1,transition(i)%number_of_mlines
	       write(stream,*) transition(i)%mline(j)%multipole
	       write(stream,*) transition(i)%mline(j)%gauge
	       write(stream,*) transition(i)%mline(j)%amplitude
	       if (reos_precalculate_integrals) then
	          need_radial_int = .true.
		  bl1 = lbound(transition(i)%mline(j)%radial_int,1)
		  bu1 = ubound(transition(i)%mline(j)%radial_int,1)
		  bl2 = lbound(transition(i)%mline(j)%radial_int,2)
		  bu2 = ubound(transition(i)%mline(j)%radial_int,2)
		  bl3 = lbound(transition(i)%mline(j)%radial_int,3)
		  bu3 = ubound(transition(i)%mline(j)%radial_int,3)
		  bl4 = lbound(transition(i)%mline(j)%radial_int,4)
		  bu4 = ubound(transition(i)%mline(j)%radial_int,4)
	       else
	          need_radial_int = .false.
	          bl1 = 0;   bu1 = 0;   bl2 = 0;   bu2 = 0
		  bl3 = 0;   bu3 = 0;   bl4 = 0;   bu4 = 0
	       end if
	       write(stream,*) need_radial_int, &
	                       bl1, bu1, bl2, bu2, bl3, bu3, bl4, bu4
	    end do
	 end if
      end do
      !
      if (read_from) then
         write(*,"(a,i7,a)") &
            "    re-calculation of ",int_counter," radial integrals complete;"
      end if
      !
   end subroutine reos_readwrite_transitions
   !
   !
   function reos_reduced_M_integral(tr,multipole,gauge,rwf_f,rwf_i) &
                                                          result(red_me)
   !--------------------------------------------------------------------
   ! Calculates the value of an (one-electron) M integral for the coupling
   ! of electrons with the radiation field for given multipole and
   ! gauge. It applies the set of Bessel functions which are associated
   ! with the transition tr; the radial wave functions are provided by the
   ! derived data structures rwf_f and rwf_i.
   !
   ! Calls: angular_momentum_j(), orbital_name(), quad_grasp2k(),
   !        wigner_3j_symbol().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: tr
      character(len=2), intent(in)       :: multipole
      character(len=9), intent(in)       :: gauge
      type(orbital_function), intent(in) :: rwf_f, rwf_i
      real(kind=dp)                      :: red_me
      !
      integer       :: ja, jb, kapa, kapb, L, mtp, phase
      real(kind=dp) :: factor, red_Coulomb, red_Gauge, wa, wb
      real(kind=dp), dimension(:), allocatable :: ta
      !
      red_me = zero
      !
      if (debug_reduced_M_integral) then
         write (99,1) orbital_name(rwf_f%orbital%n,rwf_f%orbital%kappa), &
                      orbital_name(rwf_i%orbital%n,rwf_i%orbital%kappa)
       1 format(/1x,30("+"),"   reos_reduced_M_integral() called for ",&
                "orbitals", a4, " and ", a4, ":",3x,30("+"))
      end if
      !
      kapa = rwf_f%orbital%kappa;        kapb = rwf_i%orbital%kappa
      ja   = angular_momentum_j(kapa);   jb   = angular_momentum_j(kapb)
      select case(multipole)
      case("E1", "M1");   L = 1
      case("E2", "M2");   L = 2
      case("E3", "M3");   L = 3
      case("E4", "M4");   L = 4
      case default
         stop "reos_reduced_M_integral(): program stop A."
      end select
      !
      ! Test parity rules
      select case(multipole)
      case("M1", "M2", "M3", "M4")
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 0)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 2) ) then
         else
            return
         end if
      case("E1", "E2", "E3", "E4")
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 2)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 0) ) then
         else
            return
         end if
      end select
      !
      ! Evaluate factor multiplying Mbar(a,b); use one-particle matrix elements
      ! of Brink-Satchler type
      phase  = (ja + 1)/2 + jb
      factor = ((-one)**phase) *sqrt(jb+one)*wigner_3j_symbol(ja,L+L,jb,1,0,-1)
      if (abs(factor) < eps10) then
         red_me = zero
         return
      end if
      !
      mtp = min( rwf_f%mtp, rwf_i%mtp)
      allocate( ta(1:mtp+10) )
      !
      select case(multipole)
      case("M1", "M2", "M3", "M4")
         ta(1) = zero
         ta(2:mtp) = ( rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +                   &
                       rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) ) * rp_grasp2k(2:mtp)
         select case(multipole)
         case("M1")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel1(2:mtp)
         case("M2")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel2(2:mtp)
         case("M3")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel3(2:mtp)
         case("M4")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel4(2:mtp)
         end select
         red_me = -(L+L+one)*(kapa+kapb)/sqrt(L*(L+one)) * quad_grasp2k(ta,mtp)
         red_me = factor * red_me
      case("E1", "E2", "E3", "E4")
         ta(1) = zero
         !
         ! Tabulate coulomb integrand
         !
         wa    = sqrt(L/(L+one))*(kapa-kapb);  wb = sqrt(L/(L+one))*(L+one)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         select case(multipole)
         case("E1")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel2(2:mtp)
         case("E2")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel3(2:mtp)
         case("E3")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel4(2:mtp)
         case("E4")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel5(2:mtp)
         end select
         red_Coulomb = factor * quad_grasp2k(ta,mtp)
         !
         wa    = - sqrt((L+one)/L)*(kapa-kapb);  wb = sqrt((L+one)/L)*L
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !!x print *, "-BA"
         select case(multipole)
         case("E1")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel0(2:mtp)
         case("E2")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel1(2:mtp)
         case("E3")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel2(2:mtp)
         case("E4")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel3(2:mtp)
         end select
         red_Coulomb = red_Coulomb + factor * quad_grasp2k(ta,mtp)
         !
         if (gauge == "Coulomb  ") then
            red_me = red_Coulomb
            if (debug_reduced_M_integral) &
               write(99,*) "factor, red_Coulomb = ",factor, red_me
            deallocate( ta )
            return
         else if (gauge == "Babushkin") then
         else if (rabs_use_stop) then
            stop "reos_reduced_M_integral(): program stop B."
         end if
         !
         wa    = (kapa-kapb);  wb = (L+1)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         select case(multipole)
         case("E1")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel2(2:mtp)
         case("E2")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel3(2:mtp)
         case("E3")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel4(2:mtp)
         case("E4")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel5(2:mtp)
         end select
         red_Gauge = factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         !
         wa    = (kapa-kapb);  wb = -L
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         select case(multipole)
         case("E1")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel0(2:mtp)
         case("E2")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel1(2:mtp)
         case("E3")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel2(2:mtp)
         case("E4")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel3(2:mtp)
         end select
         red_Gauge = red_Gauge + factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         !
         !!x print *, "C"
         wa    = -(L+L+one)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%P(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%Q(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         select case(multipole)
         case("E1")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel1(2:mtp)
         case("E2")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel2(2:mtp)
         case("E3")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel3(2:mtp)
         case("E4")
            ta(2:mtp) = ta(2:mtp) * transition(tr)%bessel4(2:mtp)
         end select
          red_Gauge = red_Gauge + factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         red_me    = red_Coulomb + red_Gauge
         !
      case default
         stop "reos_reduced_M_integral(): program stop C."
      end select
      !
      deallocate( ta )
      !
   end function reos_reduced_M_integral
   !
   !
   subroutine reos_select_multipoles(totalJ_i,parity_i,totalJ_f,parity_f, &
                                     mult,gauge,nmult)
   !--------------------------------------------------------------------
   ! Determines the number and type of the multipole line components
   ! which contribute and need to be considered for a given transition.
   ! It applies the standard selection rules and compares allowed line
   ! components with the list of multipoles selected at input time.
   !
   ! Calls: is_triangle().
   !--------------------------------------------------------------------
      !
      integer, intent(in)           :: totalJ_i, totalJ_f
      character(len=1), intent(in)  :: parity_i, parity_f
      integer, intent(out)          :: nmult
      character(len=2), dimension(20), intent(out) :: mult
      character(len=9), dimension(20), intent(out) :: gauge
      !
      integer :: i
      !
      nmult = 0
      do  i = 1,number_of_multipoles
         select case( reos_multipole(i) )
         case("E1")
            if (parity_i /= parity_f   .and. &
                is_triangle(totalJ_i,totalJ_f,2)) then
               nmult = nmult + 1
               gauge(nmult) = "Babushkin";   mult(nmult) = "E1"
               nmult = nmult + 1
               gauge(nmult) = "Coulomb  ";   mult(nmult) = "E1"
            end if
         case("M1")
            if (parity_i == parity_f   .and. &
                is_triangle(totalJ_i,totalJ_f,2)) then
               nmult = nmult + 1
               gauge(nmult) = "Magnetic ";   mult(nmult) = "M1"
            end if
         case("E2")
            if (parity_i == parity_f   .and. &
                is_triangle(totalJ_i,totalJ_f,4)) then
               nmult = nmult + 1
               gauge(nmult) = "Babushkin";   mult(nmult) = "E2"
               nmult = nmult + 1
               gauge(nmult) = "Coulomb  ";   mult(nmult) = "E2"
            end if
         case("M2")
            if (parity_i /= parity_f   .and. &
                is_triangle(totalJ_i,totalJ_f,4)) then
               nmult = nmult + 1
               gauge(nmult) = "Magnetic ";   mult(nmult) = "M2"
            end if
         case("E3")
            if (parity_i /= parity_f   .and. &
                is_triangle(totalJ_i,totalJ_f,6)) then
               nmult = nmult + 1
               gauge(nmult) = "Babushkin";   mult(nmult) = "E3"
               nmult = nmult + 1
               gauge(nmult) = "Coulomb  ";   mult(nmult) = "E3"
            end if
         case("M3")
            if (parity_i == parity_f   .and. &
                is_triangle(totalJ_i,totalJ_f,6)) then
               nmult = nmult + 1
               gauge(nmult) = "Magnetic ";   mult(nmult) = "M3"
            end if
         case("E4")
            if (parity_i == parity_f   .and. &
                is_triangle(totalJ_i,totalJ_f,8)) then
               nmult = nmult + 1
               gauge(nmult) = "Babushkin";   mult(nmult) = "E4"
               nmult = nmult + 1
               gauge(nmult) = "Coulomb  ";   mult(nmult) = "E4"
            end if
         case("M4")
            if (parity_i /= parity_f   .and. &
                is_triangle(totalJ_i,totalJ_f,8)) then
               nmult = nmult + 1
               gauge(nmult) = "Magnetic ";   mult(nmult) = "M4"
            end if
         case default
            stop "reos_select_multipoles(): program stop A."
         end select
      end do
      !
      if (nmult > 20) then
         stop "reos_select_multipoles(): program stop B."
      end if
      !
   end subroutine reos_select_multipoles
   !
   !
   subroutine reos_set_determinant_sym()
   !--------------------------------------------------------------------
   ! Allocates storage for information about the symmetry blocks of each
   ! determinant of the initial- and final-state determinant basis.
   ! This information later facilitates to 'recognize' easily which pairs
   ! of determinant will have non-zero contributions on transition
   ! probabilities. This information is kept in the (private and
   ! allocatable) arrays reos_detsym_i and reos_detsym_f.
   !
   ! Calls: angular_momentum_j(), unpack_occupation_from_integer().
   !--------------------------------------------------------------------
      !
      integer :: i, j, jmax, kappa, mm, norb
      integer, dimension(:), allocatable :: detsym, occupation
      !
      type(determinant) :: determ           ! For test reasons
      allocate( determ%occupation(3) )
      ! Determine the number of 'symmetry blocks' in the initial- and
      ! final-state determinant basis
      number_of_symmetry_blocks = 0
      do  kappa = -19,19
         if (kappa == 0) cycle
         jmax = angular_momentum_j(kappa)
         mm_loop: do  mm = -jmax,jmax,2
            do  i = 1,asf_initial%det_set%norbital
               if (kappa == asf_initial%det_set%orbital(i)%kappa   .and. &
                   mm    == asf_initial%det_set%orbital(i)%mm) then
                  number_of_symmetry_blocks = number_of_symmetry_blocks + 1
                  sym_block(kappa,mm) = number_of_symmetry_blocks
                  cycle mm_loop
               end if
            end do
            do  i = 1,asf_final%det_set%norbital
               if (kappa == asf_final%det_set%orbital(i)%kappa   .and. &
                   mm    == asf_final%det_set%orbital(i)%mm) then
                  number_of_symmetry_blocks = number_of_symmetry_blocks + 1
                  sym_block(kappa,mm) = number_of_symmetry_blocks
                  cycle mm_loop
               end if
            end do
         end do mm_loop
      end do
      !
      allocate( reos_detsym_i(asf_initial%det_set%nod,number_of_symmetry_blocks) )
      allocate( reos_detsym_f(asf_final%det_set%nod,number_of_symmetry_blocks) )
      !
      ! Cycle trough all determinants and determine the size of each
      ! 'symmetry block' kappa
      !
      norb = asf_initial%det_set%noint * &
             bit_size(asf_initial%det_set%determinant(1)%occupation(1))
      allocate( occupation(1:norb), detsym(number_of_symmetry_blocks) )
      do  i = 1,asf_initial%det_set%nod
         detsym = 0
         call unpack_occupation_from_integer(asf_initial%det_set%determinant(i), &
                                             occupation,norb)
         !x print *, "ADa"
         do  j = 1,asf_initial%det_set%norbital
            if (occupation(j) == 1) then
               kappa = asf_initial%det_set%orbital(j)%kappa
               mm    = asf_initial%det_set%orbital(j)%mm
               detsym(sym_block(kappa,mm)) =  detsym(sym_block(kappa,mm)) + 1
            end if
         end do
         reos_detsym_i(i,:) = detsym(:)
         if (rabs_use_stop                                  .and.  &
             sum(detsym) /= asf_initial%det_set%number_of_electrons) then
            stop "reos_set_determinant_sym(): program stop A."
         end if
      end do
      deallocate( occupation, detsym )
      print *, "   set-up of 'symmetry blocks' for the initial-states complete;"
      !
      norb = asf_final%det_set%noint * &
             bit_size(asf_final%det_set%determinant(1)%occupation(1))
      allocate( occupation(1:norb), detsym(number_of_symmetry_blocks) )
      do  i = 1,asf_final%det_set%nod
         detsym = 0
         call unpack_occupation_from_integer(asf_final%det_set%determinant(i), &
                                             occupation,norb)
         do  j = 1,asf_final%det_set%norbital
            if (occupation(j) == 1) then
               kappa = asf_final%det_set%orbital(j)%kappa
               mm    = asf_final%det_set%orbital(j)%mm
               detsym(sym_block(kappa,mm)) =  detsym(sym_block(kappa,mm)) + 1
            end if
         end do
         reos_detsym_f(i,:) = detsym(:)
         if (sum(detsym) /= asf_final%det_set%number_of_electrons) then
            print *, " "
            print *, "i, sum(detsym), asf_final%det_set%number_of_electrons = ",&
                      i, sum(detsym), asf_final%det_set%number_of_electrons
            print *, "asf_final%det_set%determinant(i)%occupation(:) = ", &
                      asf_final%det_set%determinant(i)%occupation(:)
            print *, "occupation = ",occupation(:)
            call pack_occupation_in_integer(determ,3,occupation,norb)
            print *, "determ%occupation(:) = ",determ%occupation(:)
            !! stop "reos_set_determinant_sym(): program stop B."
         end if
      end do
      deallocate( occupation, detsym )
      print *, "   set-up of 'symmetry blocks' for the final-states complete;"
      !
   end subroutine reos_set_determinant_sym
   !
   !
   subroutine reos_set_overlaps()
   !--------------------------------------------------------------------
   ! Calculates the non-orthogonal and overlap integrals between all
   ! bound orbitals of the corresponding final- and initial-state arrays.
   !
   ! Calls: rk_integral_grasp2k_ab().
   !--------------------------------------------------------------------
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
      allocate( reos_overlap(kappa_min:kappa_max,pqn_f_max,pqn_i_max) )
      !
      ! Calculate the non-orthogonal and overlap integrals
      do  i1 = 1,wave_final%number_of_rwf
         do  i2 = 1,wave_initial%number_of_rwf
            if (wave_final%rwf(i1)%orbital%kappa == &
                wave_initial%rwf(i2)%orbital%kappa) then
               kappa = wave_final%rwf(i1)%orbital%kappa
               pqn_f = wave_final%rwf(i1)%orbital%n
               pqn_i = wave_initial%rwf(i2)%orbital%n
               reos_overlap(kappa,pqn_f,pqn_i) = &
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
            end if
         end do
      end do
      !
   end subroutine reos_set_overlaps
   !
   !
   subroutine reos_set_radial_integrals()
   !--------------------------------------------------------------------
   ! Initializes a proper 'storage' for the necessary radial multipole
   ! integrals.
   !
   ! Calls: reos_calculate_Bessel(), reos_reduced_M_integral().
   !--------------------------------------------------------------------
      !
      integer       :: i, i1, i2, int_counter, j, kappa_f, kappa_i, &
                       kappa_f_min, kappa_f_max, kappa_i_min,       &
                       kappa_i_max, pqn_f, pqn_f_max, pqn_i, pqn_i_max
      !
      ! Allocate an appropriate array for the radial integrals
      kappa_i_min =  100;  kappa_i_max = -100
      kappa_f_min =  100;  kappa_f_max = -100
      pqn_i_max   = -100;  pqn_f_max   = -100
      do  i1 = 1,wave_initial%number_of_rwf
         kappa_i_min = min(kappa_i_min, wave_initial%rwf(i1)%orbital%kappa)
         kappa_i_max = max(kappa_i_max, wave_initial%rwf(i1)%orbital%kappa)
         pqn_i_max   = max(pqn_i_max,   wave_initial%rwf(i1)%orbital%n)
      end do
      do  i1 = 1,wave_final%number_of_rwf
         kappa_f_min = min(kappa_f_min, wave_final%rwf(i1)%orbital%kappa)
         kappa_f_max = max(kappa_f_max, wave_final%rwf(i1)%orbital%kappa)
         pqn_f_max   = max(pqn_f_max,   wave_final%rwf(i1)%orbital%n)
      end do
      !
      int_counter = 0
      do  i = 1,number_of_transitions
         do  j = 1,transition(i)%number_of_mlines
            allocate(transition(i)%mline(j)%radial_int(             &
                               1:pqn_f_max,kappa_f_min:kappa_f_max, &
                               1:pqn_i_max,kappa_i_min:kappa_i_max) )
            transition(i)%mline(j)%radial_int = zero
            !
            do  i1 = 1,wave_final%number_of_rwf
               do  i2 = 1,wave_initial%number_of_rwf
                  pqn_f   = wave_final%rwf(i1)%orbital%n
                  kappa_f = wave_final%rwf(i1)%orbital%kappa
                  pqn_i   = wave_initial%rwf(i2)%orbital%n
                  kappa_i = wave_initial%rwf(i2)%orbital%kappa
                  transition(i)%mline(j)%radial_int(pqn_f,kappa_f,            &
                                                    pqn_i,kappa_i) =          &
                   reos_reduced_M_integral(i,transition(i)%mline(j)%multipole,&
                                  transition(i)%mline(j)%gauge,               &
                                  wave_final%rwf(i1),wave_initial%rwf(i2))
                   if (abs(transition(i)%mline(j)%radial_int(pqn_f,kappa_f,   &
                       pqn_i,kappa_i)) > zero) int_counter = int_counter + 1
                   !
                   ! Deallocate the Bessel function arrays
                   !! if (allocated(transition(i)%bessel0)) then
                   !!    deallocate( transition(i)%bessel0 )
                   !! end if
                   !! if (allocated(transition(i)%bessel1)) then
                   !!    deallocate( transition(i)%bessel1 )
                   !! end if
                   !! if (allocated(transition(i)%bessel2)) then
                   !!    deallocate( transition(i)%bessel2 )
                   !! end if
                   !! if (allocated(transition(i)%bessel3)) then
                   !!    deallocate( transition(i)%bessel3 )
                   !! end if
                   !! if (allocated(transition(i)%bessel4)) then
                   !!    deallocate( transition(i)%bessel4 )
                   !! end if
                   !! if (allocated(transition(i)%bessel5)) then
                   !!    deallocate( transition(i)%bessel5 )
                   !! end if
               end do
            end do
         end do
      end do
      !
      write(*,"(a,i7,a)") &
         "    pre-calculation of ",int_counter," radial integrals complete;"
      !
   end subroutine reos_set_radial_integrals
   !
   !
   subroutine reos_set_transitions()
   !--------------------------------------------------------------------
   ! Determines how many and which transitions need to be calculated and
   ! initializes the array transitions of type(reos_transitions).
   ! The default is (for number_of_transitions == 0) that all transitions
   ! with a positive energy are calculated; for number_of_transitions /= 0,
   ! individual transitions were selected during input time and will be
   ! initialized instead, in this case also negative transition energies
   ! (absorption lines) are allowed.
   !
   ! Calls: reos_calculate_Bessel(), reos_select_multipoles().
   !--------------------------------------------------------------------
      !
      integer       :: asfi, asff, i, imin, j, k, m, nmult, nt
      real(kind=dp) :: arg, energy, energy_exp, energy_min
      logical       :: yes0, yes1, yes2, yes3, yes4, yes5
      character(len=2), dimension(20)       :: mult
      character(len=9), dimension(20)       :: gauge
      real(kind=dp), dimension(1:n_grasp2k) :: bessel
      !
      ! Determine the total number of transitions
      nt = 0
      if (number_of_transitions == 0) then
         do  i = 1,asf_initial%noasf
            do  j = 1,asf_final%noasf
               if (asf_initial%asf(i)%energy - asf_final%asf(j)%energy &
                                                                 > zero) then
                  call reos_select_multipoles(1*asf_initial%asf(i)%totalJ,   &
                                                asf_initial%asf(i)%parity,   &
                        1*asf_final%asf(j)%totalJ,asf_final%asf(j)%parity,   &
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
                  if ( (select_level_i(k) == asf_initial%asf(i)%level_No .or. &
                        select_level_i(k) == 0)                         .and. &
                       (select_level_f(k) == asf_final%asf(j)%level_No   .or. &
                        select_level_f(k) == 0) ) then
                     call reos_select_multipoles(1*asf_initial%asf(i)%totalJ, &
                                                   asf_initial%asf(i)%parity, &
                           1*asf_final%asf(j)%totalJ,asf_final%asf(j)%parity, &
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
               if (asf_initial%asf(i)%energy - asf_final%asf(j)%energy     &
                                                                > zero) then
                  call reos_select_multipoles(1*asf_initial%asf(i)%totalJ, &
                                                asf_initial%asf(i)%parity, &
                        1*asf_final%asf(j)%totalJ,asf_final%asf(j)%parity, &
                        mult,gauge,nmult)
                  if (nmult == 0) cycle
                  nt = nt + 1
                  transition(nt)%asfi     = i;   transition(nt)%asff = j
                  transition(nt)%level_i  = asf_initial%asf(i)%level_No
                  transition(nt)%level_f  = asf_final%asf(j)%level_No
                  transition(nt)%energy   = asf_initial%asf(i)%energy - &
                                            asf_final%asf(j)%energy
                  !
                  ! Read in experimental energies if appropriate
                  if (reos_apply_exp_energies) then
                     if (energy_inverse) then
                        energy = energy_factor / transition(nt)%energy
                     else
                        energy = energy_factor * transition(nt)%energy
                     end if
                   1 print *, "Transition",transition(nt)%level_i,"-",        &
                      transition(nt)%level_f,"has ab-initio energy E_theo = ",&
                      energy,trim(energy_unit),";"
                     print *, " enter E_exp (in ",trim(energy_unit), &
                      ") or 0.0 to use E_theo:"
                     read(*,*,err=1) energy_exp
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
                  if ( (select_level_i(k) == asf_initial%asf(i)%level_No .or. &
                        select_level_i(k) == 0)                         .and. &
                       (select_level_f(k) == asf_final%asf(j)%level_No   .or. &
                        select_level_f(k) == 0) ) then
                     call reos_select_multipoles(1*asf_initial%asf(i)%totalJ, &
                                                   asf_initial%asf(i)%parity, &
                           1*asf_final%asf(j)%totalJ,asf_final%asf(j)%parity, &
                           mult,gauge,nmult)
                     if (nmult == 0) cycle
                     nt = nt + 1
                     transition(nt)%asfi     = i;   transition(nt)%asff = j
                     transition(nt)%level_i  = asf_initial%asf(i)%level_No
                     transition(nt)%level_f  = asf_final%asf(j)%level_No
                     transition(nt)%energy   = asf_initial%asf(i)%energy - &
                                               asf_final%asf(j)%energy
                     !
                     ! Read in experimental energies if appropriate
                     if (reos_apply_exp_energies) then
                        if (energy_inverse) then
                           energy = energy_factor / transition(nt)%energy
                        else
                           energy = energy_factor * transition(nt)%energy
                        end if
                      2 print *, "Transition",transition(nt)%level_i,"-",     &
                           transition(nt)%level_f,            &
                           "has ab-initio energy E_theo = ",  &
                           energy,trim(energy_unit),";"
                        print *, " enter E_exp (in ",trim(energy_unit), &
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
      if (reos_sort_transition_energy) then
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
         call reos_select_multipoles(1*asf_initial%asf(i)%totalJ,    &
                                       asf_initial%asf(i)%parity,    &
               1*asf_final%asf(j)%totalJ,asf_final%asf(j)%parity,    &
               mult,gauge,nmult)
         if (rabs_use_stop   .and.   nmult == 0) then
            stop "reos_set_transitions(): program stop A."
         end if
         !
         transition(nt)%number_of_mlines = nmult
         allocate( transition(nt)%mline(1:nmult) )
         do  m = 1,nmult
            transition(nt)%mline(m)%multipole  = mult(m)
            transition(nt)%mline(m)%gauge      = gauge(m)
            transition(nt)%mline(m)%amplitude  = zero
         end do
         !
         ! Calculate the Bessel function for the given transition
         ! energy; first allocate the array
         yes0 = .false.;   yes1 = .false.    ! yes == Needs to be
         yes2 = .false.;   yes3 = .false.    !        calculated ?
         yes4 = .false.;   yes5 = .false.
         do  m = 1,nmult
         select case(mult(m))
         case("E1", "M1"); yes0 = .true.; yes1 = .true.; yes2 = .true.
         case("E2", "M2"); yes1 = .true.; yes2 = .true.; yes3 = .true.
         case("E3", "M3"); yes2 = .true.; yes3 = .true.; yes4 = .true.
         case("E4", "M4"); yes3 = .true.; yes4 = .true.; yes5 = .true.
         case default
            stop "reos_set_transitions(): program stop B."
         end select
         end do
         arg  = transition(nt)%energy / c
         if (yes0) then
            allocate( transition(nt)%bessel0(1:n_grasp2k) )
            call reos_calculate_Bessel(0,arg,bessel)
            transition(nt)%bessel0 = bessel
         end if
         if (yes1) then
            allocate( transition(nt)%bessel1(1:n_grasp2k) )
            call reos_calculate_Bessel(1,arg,bessel)
            transition(nt)%bessel1 = bessel
         end if
         if (yes2) then
            allocate( transition(nt)%bessel2(1:n_grasp2k) )
            call reos_calculate_Bessel(2,arg,bessel)
            transition(nt)%bessel2 = bessel
         end if
         if (yes3) then
            allocate( transition(nt)%bessel3(1:n_grasp2k) )
            call reos_calculate_Bessel(3,arg,bessel)
            transition(nt)%bessel3 = bessel
         end if
         if (yes4) then
            allocate( transition(nt)%bessel4(1:n_grasp2k) )
            call reos_calculate_Bessel(4,arg,bessel)
            transition(nt)%bessel4 = bessel
         end if
         if (yes5) then
            allocate( transition(nt)%bessel5(1:n_grasp2k) )
            call reos_calculate_Bessel(5,arg,bessel)
            transition(nt)%bessel5 = bessel
         end if
      end do
      !
      print *, "  ",number_of_transitions,                         &
               " transitions have been initialized and will be "// &
               "calculated in this run of the program."
      !
   end subroutine reos_set_transitions
   !
   !
   subroutine reos_write_trn_file()
   !--------------------------------------------------------------------
   ! Writes out a .trn transition energy and amplitude file for further
   ! data processing with the REOS program on stream 26.
   !
   ! Calls: file_open().
   !--------------------------------------------------------------------
      !
      integer           :: i, j
      character(len=3)  :: month
      character(len=8)  :: cdate
      character(len=10) :: ctime
      !
      write(26, "(i6,a)") number_of_transitions," = Number of transitions"
      write(26,*)
      do  i = 1,number_of_transitions
         write(26,1) transition(i)%level_i,  transition(i)%level_f,       &
                     transition(i)%totalJ_i, transition(i)%parity_i,      &
                     transition(i)%totalJ_f, transition(i)%parity_f,      &
                     transition(i)%number_of_mlines, transition(i)%energy,&
                     asf_initial%asf(transition(i)%asfi)%energy,          &
                     (transition(i)%mline(j)%multipole,                   &
                      transition(i)%mline(j)%gauge,                       &
                      transition(i)%mline(j)%amplitude,                   &
                      j=1,transition(i)%number_of_mlines)
      end do
    1 format(i4," -",i4,2x,2(i4,a1,1x),2x,i2,1pe14.7," a.u.   (",         &
             1pe14.7,")",3x,50(a2,1x,a9,1x,1pe14.7,3x))
      write(26,*)
      !
      ! Get the date and time of day; append this information to the reos.trn
      ! transition energy and amplitude file to facilitate identification of
      ! this file; this information is not to be read in
      !
      call date_and_time(date=cdate,time=ctime)
      month = get_month(cdate(5:6))
      write(26,*) "REOS run at "//ctime(1:2)//":"//ctime(3:4)//":"//      &
                  ctime(5:6)//" on "//month//                             &
                  " "//cdate(7:8)//" "//cdate(1:4)//"."
      !
   end subroutine reos_write_trn_file
   !
end module rabs_reos
