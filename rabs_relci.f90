module rabs_relci
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module controls the set-up and the diagonalization of the 
! relativistic Hamiltonian matrix. It comprises several routines which
! support the interactive input of control parameters and the output of
! the final results. This modules also prepares the different 
! computational modes concerning the storage, using either RAM memory or
! disc storage. In principle, is also supports the 'direct' calculation
! of the Hamiltonian matrix but, typically, this mode is not competitive
! in time with any of the other modes, mainly due to a rather slow 
! computation of the angular coefficients.
!-----------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_grasp2k
   use rabs_hamiltonian
   use rabs_input_dialog
   use rabs_nucleus
   implicit none
   !
   public  :: relci_collect_input
                 ! Collect and proceeds all input for the relativistic_CI
                 ! program.
   public  :: relci_diagonalize_block
                 ! Diagonalize a single block of the Hamiltonian matrix.
   public  :: relci_diagonalize_matrix
		 ! Controls the diagonallization of the Hamiltonian matrix 
		 ! for all blocks. 
   private :: relci_interprete_levels
                 ! Attempts to interprete the serial numbers of those energy
                 ! levels which are to be calculated from a given string.
   public  :: relci_output_eigenvectors
                 ! Writes out the calculated eigenvectors in a .mix file.
   public  :: relci_print_eigenvectors
                 ! Prints out selected eigenvectors for debugging purpose.
   public  :: relci_print_energies
                 ! Prints out the energies, splittings, and the transition
                 ! energies relative to the lowest.
   public  :: relci_print_results
                 ! Prints out the results of the RELCI program to stream.
   public  :: relci_print_summary
                 ! Appends a summary of the input data to the .sum file.
   private :: relci_print_weights
                 ! Print  the  weights of the largest five CSF contributors 
                 ! to each ASF.
   private :: relci_qed_F
                 ! Estimates the function  F (Z*\alpha) by using an 
                 ! interpolation of tabulated data.
   private :: relci_qed_F_Klarsfeld
                 ! Estimates the function  F (Z*\alpha) by using an 
                 ! interpolation of tabulated data.from Klarsfeld and Maquet
                 ! (1973); currently not in use.
   private :: relci_qed_F_Mohr
                 ! Computes the function  F (Z*\alpha) for the  1s  2s  2p-  
                 ! and 2p symmetries from tables due to Mohr, ADNDT (1983).
   private :: relci_qed_F_Mohr_Kim
                 ! Computes the function  F (Z*\alpha) for the  3nlj, 4nlj, and
                 ! 5nlj symmetries (lj = s, p-, p, and d-) due to  Mohr
                 ! and Kim, PRA (1992).
   public ::  relci_selfenergy
                 ! Calculates the QED estimations for all selected ASF from 
                 ! a hydrogen-like model which utilizes the charge-density
                 ! ratio (with respect to a pure hydrogenic orbital) near the
                 ! nucleus.
   private :: relci_set_eigenpairs
                 ! Set the selected eigenpairs from the energies obtained for
                 ! the individual J^P blocks of the Hamiltonian matrix.
   !
   ! Define global logical flags for the control of the RELCI program; the
   ! default values for these flags may be overwritten interactively during 
   ! input time
   logical, public ::  relci_print_qed                   = .false.,  &
                       relci_use_formatted_rwf_file      = .false.,  &
                       relci_use_formatted_mix_file      = .true.,   &                      
                       relci_add_qed_to_mix              = .false.
   !
contains
   !
   subroutine relci_collect_input()
   !--------------------------------------------------------------------
   ! Collect and proceeds all input for the set-up and the diagonaization
   ! of the (relativistic) Hamiltonian matrix.
   !
   ! Calls: get_yes_stream().
   !--------------------------------------------------------------------
      !
      integer            :: i
      logical            :: yes, fail
      character(len=256) :: record
      !
      ! Open, check, load data from, and close the  .iso  file
      call file_get_isodat_grasp2k()
      !
      energy_unit    = "eV     "
      energy_factor  = convert_au_to_ev
      energy_inverse = .false.
      !
      ! Determine the parameters controlling the radial grid
      call input_grid_parameters("standard")
      !
      ! Default speed of light
      c = c_vacuum
      !
      ! Now 'overwrite' defaults only if required
      print *, "Modify default set-up and printout of the program ?"
      yes = get_yes_stream()
      if (.not.yes) goto 6
      !
      hamiltonian_XL_coulomb = .true.
      !
      hamiltonian_XL_tbreit = .false.
	 print *, "Include contributions of the frequency-independent Breit"//&
		  " interaction ?"
	 yes = get_yes_stream()
	 if (yes) then
	    hamiltonian_XL_breit0 = .true.
	 else
	    print *, "Include contributions of the Gaunt interaction ?"
	    yes = get_yes_stream()
	    if (yes) hamiltonian_XL_gaunt = .true.
	 end if
      !
      print *, "Include vacuum polarization contributions to H ?"
      hamiltonian_vacuum_pol = get_yes_stream()
      !
      hamiltonian_normal_ms = .false.
      !
      print *, "Include specific mass shift contributions to H ?"
      hamiltonian_specific_ms = get_yes_stream()
      !
      print *, "Estimate contributions from self-energy ?"
      yes = get_yes_stream()
      if (yes) then
         hamiltonian_self_energy = .true.
         !
         print *, "Add self-energy estimates to the total energies in the "//&
                  ".mix output file ?"
         print *, " This can be useful in the computation of transition "  //&
                  "energies."
         relci_add_qed_to_mix = get_yes_stream()
      end if
      !
      ! Determine the storage management for diagonalizing the matrix
      print *, "Store the effective interaction strengths in memory ?"
      print *, " This should always be true for a on-fly calculation"// &
               " of the Hamiltonian matrix."
      hamiltonian_use_storage = get_yes_stream()
      !
    2 print *, "Precalculate and keep Hamiltonian matrix in memory ?"
      yes = get_yes_stream()
      if (yes) then
         hamiltonian_use_memory = .true.
         hamiltonian_use_disc   = .false.
      else
         hamiltonian_use_memory = .false.
         print *, "Precalculate and use disc storage for the Hamiltonian"// &
                  " matrix ?"
         yes = get_yes_stream()
         if (yes) then
            hamiltonian_use_disc = .true.
            !
            ! Open the  .res file
            call file_open_unformatted_stream(27,             &
	       "Enter a file name for the  relci.dsc  file:")
         else
            print *, "Recalculate the Hamiltonian matrix all times"// &
                     " without storage ?"
            yes = get_yes_stream()
            if (yes) then
               hamiltonian_recalculate_alltime = .true.
            end if
         end if
      end if
      !
      if (.not.hamiltonian_use_memory  .and.  .not.hamiltonian_use_disc  .and.& 
          .not.hamiltonian_recalculate_alltime) then
         print *, "The Hamiltonian matrix has to be stored in memory or on "//&
                  "disc or has to be re-calculate all times;"
         print *, " choose on of these modes ..."
         goto 2
      end if
      !
      print *, "Diagonalize the full Hamiltonian matrix ?"
      print *, "(I.e. independent of the parameter hamiltonian_fullmatrix.)"
      hamiltonian_use_full_diagonal = get_yes_stream()
      !
      ! Determine the units for the printout and further optional input
      call input_energy_unit()
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
      ! Get radial orbitals of the Atomic State Functions
      call file_get_rwf(					    &
	 "Enter the name of the GRASP2K Radial WaveFunction File:", &
	 asf_bound,wave_bound,.false.)
      !
      ! Open a .mix Mixing Coefficient Output file
      call file_open_formatted_stream(26, &
	 "Enter a file name for the  relci.mix  file:")
      !
      !
    7 print *, "Enter the serial number(s) of the level(s) to be calculated;"
      print *, " e.g. 1 3 4  7 - 20  48  69 - 85;"
      read (*, "(a)") record
      call relci_interprete_levels(record,hamiltonian_eigenpair,              &
                      hamiltonian_no_eigenpairs,hamiltonian_max_eigenpair,fail)
      if (fail) then
         print *, "Unable to interprete the serial level numbers; redo ..."
         goto 7
      else if (hamiltonian_max_eigenpair > asf_bound%csf_set%nocsf) then
         print *, "Maximal serial level number must be smaller than the"// &
                  " number of CSF; redo ..."
         goto 7
      end if
      !
      call save_input(record,.false.)
      !
      ! Allocate memory and initialize all eigenvectors to zero
      asf_bound%noasf = hamiltonian_no_eigenpairs
      allocate( asf_bound%asf(1:asf_bound%noasf) )
      do  i = 1,asf_bound%noasf
         allocate( asf_bound%asf(i)%eigenvector(1:asf_bound%csf_set%nocsf) )
         asf_bound%asf(i)%eigenvector(1:asf_bound%csf_set%nocsf) = zero
         asf_bound%asf(i)%level_No = hamiltonian_eigenpair(i)
      end do
      !
   end subroutine relci_collect_input
   !
   !
   subroutine relci_diagonalize_block(hblock)
   !--------------------------------------------------------------------
   ! Manages the diagonalization of a single block of the Hamiltonian
   ! matrix.
   !--------------------------------------------------------------------
      !
      type(hamiltonian_block), intent(inout)     :: hblock
      !
      integer       :: ierr, info, iiwsz, iwrsz, m, mblock, maxiter, lim,   &
                       ndx_low, ndx_high, no_eigenvectors, no_estimates,    &
                       nloops, nume, nmv, p, q, s, ss
      real(kind=dp) :: crite, critr, ortho
      logical       :: hiend
      !
      integer, dimension(1:hblock%nocsf)         :: ndx_select
      integer, dimension(:), allocatable         :: iwork, ifail
      real(kind=dp),dimension(1:hblock%nocsf)    :: diagonal_elements
      real(kind=dp), dimension(:), allocatable   :: ap, eigenvalue, work,   &
                                                    full_column
      real(kind=dp), dimension(:,:), allocatable :: z
      external relci_op
      !
      if (hblock%nocsf == 1) then
         !
         ! Trivial case if hblock%nocsf = 1
         ! ================================
         if (abs(H_matrix(1)%me(1)) < eps10) then
            stop "hamiltonian_diagonalize_block(): program stop A."
         end if
         hblock%eigenvalue(1)    = H_matrix(1)%me(1)
         hblock%eigenvector(1,1) = one
         hamiltonian_average     = hamiltonian_average + H_matrix(1)%me(1)
         !
         print *,    " "
         print *,    "Block J^P = ("//                                       &
                     trim(angular_momentum_string(hblock%totalJ*1,4))//","// &
                     hblock%parity//"):  "//                                 &
                     "Trivial case of dimension 1."
         write(24,*) " "
         write(24,*) "Block J^P = ("//                                       &
                     trim(angular_momentum_string(hblock%totalJ*1,4))//","// &
                     hblock%parity//"):  "//                                 &
                     "Trivial case of dimension 1."
         !
      else if (hblock%nocsf <= hamiltonian_fullmatrix   .or. &
               hamiltonian_use_full_diagonal) then
         !
         ! Use a full diagonalization of the corresponding block
         ! =====================================================
         allocate( ap(1:(hblock%nocsf*(hblock%nocsf+1)/2)) )
         allocate( eigenvalue(1:hblock%nocsf), work(1:8*hblock%nocsf) )
         allocate( ifail(1:hblock%nocsf), iwork(1:5*hblock%nocsf) )
         allocate( full_column(1:hblock%nocsf) )
         allocate( z(1:hblock%nocsf,1:hblock%nocsf) )
         !
         ss = 0
         do  s = 1,hblock%nocsf
            full_column(:) = zero
            do  p = 1,H_matrix(s)%non_zero_me
               q              = H_matrix(s)%ndx(p)
               full_column(q) = H_matrix(s)%me(p)
            end do
            ap(ss+1:ss+hblock%nocsf-s+1) = full_column(s:hblock%nocsf)
            ss = ss + hblock%nocsf-s+1
            hamiltonian_average = hamiltonian_average + full_column(s)
         end do
         !
         if (rabs_use_stop .and. ss /= (hblock%nocsf*(hblock%nocsf+1)/2) ) then
            stop "hamiltonian_diagonalize_block(): program stop B."
         end if
         !
         call dspevx("V", "A", "L", hblock%nocsf, ap, zero, zero, 0, 0,  &
                     eps10, m, eigenvalue, z, hblock%nocsf, work, iwork, &
                     ifail, info)
         if (info /= 0) then
            print *, "hblock%nocsf, hblock%totalJ, hblock%parity, info = ", &
                      hblock%nocsf, hblock%totalJ, hblock%parity, info
         end if
         !
         no_eigenvectors = min(hamiltonian_max_eigenpair,hblock%nocsf)
         do  s = 1,no_eigenvectors
            print *, "eigenvalue = ",eigenvalue(s)
            hblock%eigenvalue(s) = eigenvalue(s)
            hblock%eigenvector(s,1:hblock%nocsf) = z(1:hblock%nocsf,s)
         end do
         !
         deallocate( ifail, iwork, ap, eigenvalue, work, z )
         !
         print *,    " "
         print *,    "Block J^P = ("//                                         &
                     trim(angular_momentum_string(hblock%totalJ*1,4))//","//   &
                     hblock%parity//"):  "//                                   &
                     "LAPACK routine DSPEVX selected for eigenvalue problem "//&
                     "of dimension ",hblock%nocsf
         write(24,*) " "
         write(24,*) "Block J^P = ("//                                         &
                     trim(angular_momentum_string(hblock%totalJ*1,4))//","//   &
                     hblock%parity//"):  "//                                   &
                     "LAPACK routine DSPEVX selected for eigenvalue problem "//&
                     "of dimension ",hblock%nocsf
         !
      else
         !
         ! Call Davidson eigensolver
         ! =========================
         ! Allocate storage for workspace; see the header of DVDSON for the 
         ! expression below; the value of LIM can be reduced to 
         ! hamiltonian_max_eigenpair plus a smaller number if storage is 
         ! severely constrained
         lim  = min(hblock%nocsf,hamiltonian_max_eigenpair+20)
         nume = min(hblock%nocsf,hamiltonian_max_eigenpair)
         ndx_low       = 1
         ndx_high      = nume
         ndx_select(:) = 0
         no_estimates  = 0
         mblock        = nume
         crite         = 1.0e-17_dp
         critr         = 1.0e-09_dp
         ortho         = 1.0e-09_dp
         maxiter       = max(1000,lim+50)
         iwrsz	       = 2*hblock%nocsf*lim + lim*lim + &
                         (nume+10)*lim + nume + 1000
         iiwsz	       = 6*lim + nume + 100         
         nloops        = zero 
         ierr          = 0                       
         allocate( work(1:iwrsz), iwork(1:iiwsz) )
         !
         hiend         = .false.
         !
         ! Determine diagonal elements of the Hamiltonian matrix
         do  s = 1,hblock%nocsf
            diagonal_elements(s) = hamiltonian_full_me(s,s)
            hamiltonian_average  = hamiltonian_average + diagonal_elements(s)
         end do
         !
         call dvdson(relci_op,hblock%nocsf,lim,diagonal_elements,      &
                     ndx_low,ndx_high,ndx_select,no_estimates,mblock,  &
                     crite,critr,ortho,maxiter,work,iwrsz,iwork,iiwsz, &
                     hiend,nloops,nmv,ierr)
         !
         if (ierr /= 0) then
            print *, "DVDSON returned from diagonalization of the matrix "//&
                     "with ERROR code ",ierr
            print *, " "
            print *, "RELCI terminates ..."
            stop
         end if
         !
         no_eigenvectors = min(hamiltonian_max_eigenpair,hblock%nocsf)
         do  s = 1,no_eigenvectors
            hblock%eigenvalue(s) = work(nume*hblock%nocsf+s)
            hblock%eigenvector(s,1:hblock%nocsf) = work((s-1)*hblock%nocsf+1: &
                                                        s*hblock%nocsf)
            !! print *, "eigenvalue = ",hblock%eigenvalue(s)
         end do
         deallocate( work, iwork )
         !
         print *,    " "
         print *,    "Block J^P = ("//                                       &
                     trim(angular_momentum_string(hblock%totalJ*1,4))//","// &
                     hblock%parity//"):  "//                                 &
                     "DVDSON routine selected for eigenvalue problem of "  //&
                     "dimension ",hblock%nocsf
         write(24,*) " "
         write(24,*) "Block J^P = ("//                                       &
                     trim(angular_momentum_string(hblock%totalJ*1,4))//","// &
                     hblock%parity//"):  "//                                 &
                     "DVDSON routine selected for eigenvalue problem of "  //&
                     "dimension ",hblock%nocsf
         !
      end if
      !
   end subroutine relci_diagonalize_block
   !
   !
   subroutine relci_diagonalize_matrix(asf_set,wave)
   !--------------------------------------------------------------------
   ! Calculates and diagonalize all (J,parity) blocks of the Hamiltonian 
   ! matrix in turn. For relci_precalculate_matrix = .false., the raws of 
   ! the corresponding block are re-calculated each time they are needed.
   !--------------------------------------------------------------------
      !   
      type(asf_basis), intent(inout)    :: asf_set
      type(grasp2k_orbital), intent(in) :: wave
      !
      integer :: column, i, iblock, mtp, nme
      integer, dimension(:), allocatable       :: ndx
      real(kind=dp), dimension(:), allocatable :: me
      type(hamiltonian_column)                 :: H_column
      !
      ! call cpu_time(time1)
      !
      hamiltonian_average = zero
      !
      ! Initialize the storage management and the vacuum polarization 
      ! potential
      if (hamiltonian_use_storage) then
         call hamiltonian_initialize(wave_bound)
      end if
      !
      if (hamiltonian_vacuum_pol) then
         call ncharg_grasp2k(mtp)
         allocate(vacpol_2_grasp2k(1:n_grasp2k), vacpol_4_grasp2k(1:n_grasp2k))
         call hamiltonian_set_vacpol(mtp)
         zdist_grasp2k(2:n_grasp2k) = vacpol_2_grasp2k(2:n_grasp2k) * &
                                      rp_grasp2k(2:n_grasp2k)
         deallocate( vacpol_2_grasp2k, vacpol_4_grasp2k )
      end if
      !
      ! Set up and diagonalize the Hamiltonian matrix for all (J,parity) blocks
      do  iblock = 1,hamiltonian_noblock
         hamiltonian_iblock = iblock
         if (hamiltonian_use_memory  .or. hamiltonian_use_full_diagonal  .or. &
             H_block(iblock)%nocsf <= hamiltonian_fullmatrix) then
            ! 
            ! Precalculate the full matrix of this (J,parity) block
            allocate( H_matrix(1:H_block(iblock)%nocsf) )
            call hamiltonian_set_matrix(H_block(iblock),asf_set%csf_set,wave)
         else if (hamiltonian_use_disc) then
            ! 
            ! Calculate each column in turn and write to disc
            rewind(27)
            allocate( ndx(1:H_block(iblock)%nocsf), &
                      me(1:H_block(iblock)%nocsf) )
            !
            do  column = 1,H_block(iblock)%nocsf
               call hamiltonian_set_column(H_block(iblock),column, &
                                           asf_set%csf_set,wave,nme,ndx,me)
               ! Pack results into a 'reduced' column
               allocate( H_column%ndx(1:nme), H_column%me(1:nme) )
               write(27) column, nme
               do  i = 1,nme
                  write(27) ndx(i), me(i)
               end do
               deallocate( H_column%ndx, H_column%me )
            end do
            rewind(27)
            deallocate( ndx, me )
         end if
         !
         ! Now diagonalize this block
         call relci_diagonalize_block(H_block(iblock))
         !
         if (hamiltonian_use_memory  .or.  hamiltonian_use_full_diagonal  .or.&
             H_block(iblock)%nocsf <= hamiltonian_fullmatrix) then
            do  i = 1,H_block(iblock)%nocsf
               deallocate( H_matrix(i)%ndx )
               deallocate( H_matrix(i)%me )
            end do
            deallocate( H_matrix )
         end if
      end do
      !
      ! Set the selected eigenpairs in asf_set
      call relci_set_eigenpairs(asf_set)
      !  
      !
   end subroutine relci_diagonalize_matrix
   !
   !
   subroutine relci_interprete_levels(record,eigenpair,no_eigenpairs,  &
                                                      max_eigenpair,fail)
   !--------------------------------------------------------------------
   ! Attempts to interprete the serial level numbers which are to be 
   ! calculated in the present run from a given string of numbers and
   ! 'intervals'. Levels can be given in the format:
   !                1 3 4  7 - 20  48  69 - 85
   ! The level numbers can be given in any order and 'overlapping' 
   ! intervals are also allowed.
   ! The procedure returns with fail = .true. if the level numbers cannot
   ! be interpreted properly (fail = .false. otherwise).
   !
   ! Calls: .
   !--------------------------------------------------------------------
      !
      character(len=*), intent(in)        :: record
      logical, intent(out)                :: fail
      integer, intent(out)                :: no_eigenpairs, max_eigenpair
      integer, dimension(:), pointer      :: eigenpair
      !
      logical, dimension(200)             :: low_to
      character(len=500)                  :: string
      integer                             :: a, i, lower, n
      integer, dimension(200)             :: low
      integer(kind=i1b), dimension(10000) :: eigpair
      !
      fail       = .true.
      !
      string = adjustl(record)
      !
      n = 0;   lower = 0
      low_to(:) = .false.
      !
    1 string = adjustl(string)
      if (string(1:1) == " ") then
         if (n == 0  .or.  low_to(n)) then
            return
         else
            goto 10
         end if
      else if (string(1:1) == "-") then
         if (n == 0)  return
         low_to(n)   = .true.
         string(1:1) = " "
         goto 1
      else if (string(1:1) == "0"   .or.   string(1:1) == "1" .or.  &
               string(1:1) == "2"   .or.   string(1:1) == "3" .or.  &
               string(1:1) == "4"   .or.   string(1:1) == "5" .or.  &
               string(1:1) == "6"   .or.   string(1:1) == "7" .or.  &
               string(1:1) == "8"   .or.   string(1:1) == "9") then
         a = get_integer_from_string(string(1:1))
         lower = 10*lower + a
         if (string(2:2) == " "   .or.   string(2:2) == "-") then
            n = n + 1
            low(n)      = lower
            lower       = 0
         end if
         string(1:1) = " "
         goto 1
      end  if
      !
      ! Determine no_eigenpairs and max_eigenpair
   10 eigpair(:) = 0
      do  i = 1,n
         if (low_to(i)) then
            if (low(i) <= low(i+1)) then;   eigpair(low(i):low(i+1)) = 1
            else;                           eigpair(low(i+1):low(i)) = 1
            end if
            !   
            if (low_to(i+1)) then;   return
            else;                    cycle
            end if
         else
            eigpair(low(i)) = 1
         end if
      end do
      !
      no_eigenpairs = 0;  max_eigenpair = 0
      do  i = 1,10000
         if (eigpair(i) == 1) then
            no_eigenpairs = no_eigenpairs + 1
            max_eigenpair = max( i, max_eigenpair )
         end if
      end do
      !
      allocate( eigenpair(1:no_eigenpairs) )     
      !
      no_eigenpairs = 0
      do  i = 1,10000
         if (eigpair(i) == 1) then
            no_eigenpairs = no_eigenpairs + 1
            eigenpair(no_eigenpairs) = i
         end if
      end do
      !
      fail = .false.
      !
   end subroutine relci_interprete_levels
   !
   !
   subroutine relci_output_eigenvectors(asf_set,wave)
   !--------------------------------------------------------------------
   ! Writes out the calculated eigenvectors in a .mix file and close
   ! this file.
   !--------------------------------------------------------------------
      !
      type(asf_basis), intent(inout)    :: asf_set
      type(grasp2k_orbital), intent(in) :: wave
      !
      integer                                   :: i, j, noasf, noasf_sofar
      integer, dimension(1:asf_set%noasf)       :: iaspar, iatjpo
      real(kind=dp), dimension(1:asf_set%noasf) :: shift, eval
      !   
      ! Add QED estimates to the total energies if required      
      if (hamiltonian_self_energy  .and.   relci_add_qed_to_mix) then
         call relci_selfenergy(asf_set,wave,shift,.false.)
         !
         do  i = 1,asf_set%noasf
            asf_set%asf(i)%energy = asf_set%asf(i)%energy + shift(i)
         end do
      else
         shift(:) = zero
      end if
      !
      if (relci_use_formatted_mix_file) then
         !
         ! Write out the file header
         write (26,"(a)") "G92MIX (formatted file version)."
         !
         write (26,"(3i6)") asf_set%csf_set%number_of_electrons, &
                            asf_set%csf_set%nocsf, asf_set%csf_set%nwshells
         write (26,"(i6)")  asf_set%noasf
         !
         noasf = min(60,asf_set%noasf)
         write (26,"(60i6)") (asf_set%asf(i)%level_No,i = 1,noasf)
         write (26,"(60(i5,a1))") (asf_set%asf(i)%totalJ,asf_set%asf(i)%parity, &
                                   i = 1,noasf)
         !
         write (26,"(61e26.19)") asf_set%average_energy, &
                                (asf_set%asf(i)%energy,i = 1,noasf)
         do  j = 1,asf_set%csf_set%nocsf
            write (26,"(60e16.9)") (asf_set%asf(i)%eigenvector(j), &
                                     i=1,noasf)
         end do
         !
         if (asf_set%noasf <= 60) goto 10
         !
         noasf_sofar = 0
       1 noasf_sofar = noasf_sofar + 60
         !
         noasf = min(noasf_sofar+60,asf_set%noasf)
         write (26,"(60i6)") (asf_set%asf(i)%level_No,i = noasf_sofar+1,noasf)
         write (26,"(60(i5,a1))") (asf_set%asf(i)%totalJ,asf_set%asf(i)%parity, &
                                   i = noasf_sofar+1,noasf)
         !
         write (26,"(60e26.19)") (asf_set%asf(i)%energy,i = noasf_sofar+1,noasf)
         do  j = 1,asf_set%csf_set%nocsf
            write (26,"(60e16.9)") (asf_set%asf(i)%eigenvector(j), &
                                     i=noasf_sofar+1,noasf)
         end do
         !
         if (asf_set%noasf > noasf_sofar+60) goto 1
         !

      else
         ! Write out the file header
         write (26) "G92MIX"
         write (26) asf_set%csf_set%number_of_electrons, &
                    asf_set%csf_set%nocsf, asf_set%csf_set%nwshells
         write (26) asf_set%noasf         
         !
         do  i = 1,asf_set%noasf
            iatjpo(i) = asf_set%asf(i)%totalJ + 1
            eval(i)   = asf_set%asf(i)%energy - asf_set%average_energy
            if (asf_set%asf(i)%parity == "+") then
               iaspar(i) =  1
            else if (asf_set%asf(i)%parity == "-") then
               iaspar(i) = -1
            end if
         end do
         !
         write (26) (asf_set%asf(i)%level_No,i = 1,asf_set%noasf)
         write (26) (iatjpo(i),iaspar(i),i = 1,asf_set%noasf)
         write (26) asf_set%average_energy,              &
                    (eval(i),i = 1,asf_set%noasf)
         !x write (26) ((evec(i+(j-1)*ncf),i = 1,ncf),j = 1,nvec)
         write (26) ((asf_set%asf(i)%eigenvector(j),j=1,asf_set%csf_set%nocsf),&
                    i=1,asf_set%noasf)
         !
      end if
      !
   10 continue
      !
      close(26)
      !
      print *, " "
      print *, "RELCI Mixing Coefficient File generated."  
      !
   end subroutine relci_output_eigenvectors
   !
   !
   subroutine relci_print_eigenvectors()
   !--------------------------------------------------------------------
   ! Prints out selected eigenvectors for debugging purpose.
   ! Start from rci92/evcout.f if necessary.
   !--------------------------------------------------------------------
      !
      !
   end subroutine relci_print_eigenvectors
   !
   !
   subroutine relci_print_energies(stream,mode,asf_set,shift)
   !--------------------------------------------------------------------
   ! Prints the energy levels, their splitting, and transition energies
   ! relative to the lowest in  Hartrees, Kaysers, and  eV, using the
   ! reduced mass corrected value for the Rydberg. This procedure is 
   ! similar to the GRASP2K routine ENGOUT in that different modes are 
   ! supported.
   !
   ! If mode == 0, only eigenenergies are printed.
   !    mode == 1, print all selected eigenenergies and separations.
   !    mode == 2, eigenenergies and energies relative to level 1 are 
   !               printed.
   !    mode == 3, the eigenenergies, separations, and energies relative 
   !               to level 1 are printed.
   !
   !--------------------------------------------------------------------
      !
      integer, intent(in)         :: stream, mode
      type(asf_basis), intent(in) :: asf_set
      real(kind=dp), dimension(:) :: shift
      !
      integer          :: i
      logical          :: en_inverse
      real(kind=dp)    :: energy, energy_un, conversion
      character(len=7) :: en_unit
      !
      if (energy_unit == "eV	 "  .or.  &
          energy_unit == "Hartree")	  then
         conversion  = convert_au_to_kaysers
         en_inverse  = .false.
         en_unit     = "Kaysers"
      else
         conversion  = energy_factor
         en_inverse  = energy_inverse
         en_unit     = energy_unit
      end if
      !
      ! Always print the eigenenergies
      write(stream,*)
      write(stream,*) "Eigenenergies:"
      write(stream,1) en_unit
      do  i = 1,hamiltonian_no_eigenpairs
         energy = asf_set%asf(i)%energy + shift(i)
         if (en_inverse) then
            energy_un = conversion / energy
         else
            energy_un = conversion * energy
         end if
         !
         write(stream,2) asf_set%asf(i)%level_No,                           &
                  trim(angular_momentum_string(1*asf_set%asf(i)%totalJ,4)), &
                     asf_set%asf(i)%parity, energy,                         &
                     energy*convert_au_to_ev, energy_un
      end do
    1 format(/2x,"Level  J Parity",10x,"Hartrees",20x,"eV",20x,a7 /)
    2 format( 2x,1i3,4x,2a4,1p,3d25.15)
      !
      if (hamiltonian_no_eigenpairs > 1) then
         !
         ! Energy separations
         if (mode == 1   .or.   mode == 3) then
            write(stream,*)
            write(stream,*) "Energy of each level relative to immediately lower"//&
                        " level:"
            write(stream,1) en_unit
            do  i = 2,hamiltonian_no_eigenpairs
               energy = asf_set%asf(i)%energy - asf_set%asf(i-1)%energy +   &
                        shift(i) - shift(i-1)
               if (en_inverse) then
                  energy_un = conversion / energy
               else
                  energy_un = conversion * energy
               end if
               !
               write(stream,2) asf_set%asf(i)%level_No,                     &
                  trim(angular_momentum_string(1*asf_set%asf(i)%totalJ,4)), &
                           asf_set%asf(i)%parity, energy,                   &
                           energy*convert_au_to_ev, energy_un
            end do
         end if
         !
         ! Energies relative to level 1
         if (mode == 2   .or.   mode == 3) then
            write(stream,*)
            write(stream,*) "Energy of each level relative to lowest level:"
            write(stream,1) en_unit
            do  i = 2,hamiltonian_no_eigenpairs
               energy = asf_set%asf(i)%energy - asf_set%asf(1)%energy +     &
                        shift(i) - shift(1)
               if (en_inverse) then
                  energy_un = conversion / energy
               else
                  energy_un = conversion * energy
               end if
               !
               write(stream,2) asf_set%asf(i)%level_No,                     &
                  trim(angular_momentum_string(1*asf_set%asf(i)%totalJ,4)), &
                           asf_set%asf(i)%parity, energy,                   &
                           energy*convert_au_to_ev, energy_un
            end do
         end if
      end if
      !
   end subroutine relci_print_energies
   !
   !
   subroutine relci_print_results(stream,asf_set,wave)
   !--------------------------------------------------------------------
   ! Prints out the results of the RELCI program to stream.
   !--------------------------------------------------------------------
      !
      integer, intent(in)               :: stream
      type(asf_basis), intent(in)       :: asf_set
      type(grasp2k_orbital), intent(in) :: wave
      !
      real(kind=dp), dimension(1:asf_set%noasf) :: shift
      !
      ! Print a summary about the calculated integrals
      write(stream,*) " "
      write(stream,*) "Number of Dirac-Coulomb one-electron integrals " //&
                      "computed       = ",hamiltonian_oneint_ener_comput
      if (hamiltonian_vacuum_pol) then
      write(stream,*) "Number of one-electron vacuum-polarization inte" //&
                      "grals computed = ",hamiltonian_oneint_vp_comput
      end if
      if (hamiltonian_normal_ms) then
      write(stream,*) "Number of one-electron normal mass shift integr" //&
                      "als computed   = ",hamiltonian_oneint_nms_comput
      end if
      !
      if (hamiltonian_use_storage) then
      write(stream,*) "Number of (full) one-electron matrix elements s" //&
                      "tored          = ",hamiltonian_one_me_stored
      write(stream,*) "Number of (full) one-electron matrix elements r" //&
                      "e-used         = ",hamiltonian_one_me_reused
      end if
      !
      write(stream,*) " "
      write(stream,*) "Number of (full) two-electron X^k strengths com" //&
                      "puted          = ",hamiltonian_twoint_Xk_comput
      if (hamiltonian_use_storage) then
      write(stream,*) "Number of (full) two-electron X^k strengths sto" //&
                      "red            = ",hamiltonian_twoint_Xk_stored
      write(stream,*) "Number of (full) two-electron X^k strengths re-" //&
                      "used           = ",hamiltonian_twoint_Xk_reused
      end if
      !
      if (hamiltonian_H_me_refered > 0) then
      write(stream,*) " "
      write(stream,*) "Total number of Hamiltonian matrix elmenents re" //&
                      "fered to during the diagonalization procedure = ", &
                      hamiltonian_H_me_refered
      end if
      !
      ! Write out the average energy
      write(stream,*) " "
      write(stream,*) "Average energy = ",hamiltonian_average / &
                                          asf_bound%csf_set%nocsf
      !
      ! Write out eigenvalues and the dominant components of the eigenvectors
      shift(:) = zero
      call relci_print_energies(stream,3,asf_set,shift)
      call relci_print_weights(stream,asf_set)
      !
      ! Estimate diagonal self-energy corrections; the stored eigenvalues and 
      ! eigenvectors are not modified by these estimates
      if (hamiltonian_self_energy) then
         if (stream == 6) then
            write(stream,*) " "
            write(stream,*) "Entering QED ..."
            write(stream,*) " "
            write(stream,*) "   Orbital     Self-energy estimate    " //&
                            "  charge density ratio       F(alphaZ) "
            write(stream,*) "                    (Hartree)          " //&
                            " this orbital / H--like                "
            write(stream,*) " "
            call relci_selfenergy(asf_set,wave,shift,.true.)
            write(stream,*) " "
            write(stream,*) "... QED complete."
         else 
            call relci_selfenergy(asf_set,wave,shift,.false.)
         end if
         !
         write(stream,*) " "
         !
         if (relci_add_qed_to_mix) then
            write(stream,*) "Self-energy corrections estimated --- they are "//&
                            " added to the total energies "
            write(stream,*) " in the RELCI mixing coefficients file."
         else
            write(stream,*) "Self-energy corrections estimated --- these do "//&
                            "not influence the data "
            write(stream,*) " in the RELCI mixing coefficients file."
         end if
         !
         call relci_print_energies(stream,3,asf_set,shift)
      end if
      !
   end subroutine relci_print_results
   !
   !
   subroutine relci_print_summary()
   !--------------------------------------------------------------------
   ! Appends further information to the  relci.sum  file which is open 
   ! on stream 24.
   !--------------------------------------------------------------------
      !
      implicit none
      integer           :: i      
      !
      ! Write out the physical effects specifications
      write(24,*)
      write(24,*) "Speed of light = ",c," atomic units."
      write(24,*)
      write(24,*) "The Hamiltonian matrix includes the following contributions:"
      if (hamiltonian_XL_coulomb) then
         write(24,*) " H (Dirac Coulomb)"
      end if
      if (hamiltonian_XL_breit0) then
         write(24,*) " H (Breit - low_frequency limit; omega --> 0)"
      end if
      if (hamiltonian_XL_tbreit) then
         write(24,*) " H (Transverse Breit)"
      end if
      if (hamiltonian_XL_Gaunt) then
         write(24,*) " H (Gaunt)"
      end if
      if (hamiltonian_vacuum_pol) then
         write(24,*) " H (Vacuum Polarization)"
      end if
      if (hamiltonian_normal_ms) then
         write(24,*) " H (Normal Mass Shift)"
      end if
      if (hamiltonian_specific_ms) then
         write(24,*) " H (Specific Mass Shift)"
      end if
      write(24,*) " the total will be diagonalised."
      !
      write(24,*)
      if (hamiltonian_self_energy) then
         write(24,*) "Diagonal contributions from H (Self Energy) will be"//&
                     " estimated from a screened"
         write(24,*) " hydrogenic approximation."           
      end if
      !
      ! Write the list of eigenpair indices
      write(24,*)
      if (hamiltonian_no_eigenpairs == 1) then
         write(24,*) "Level ",hamiltonian_eigenpair(1)," will be computed."
      else
         write(24,*) "A total of ",hamiltonian_no_eigenpairs," levels"// &
                     " will be computed; their indices are:"
         write(24,*) "  "           
         write(24,*) "  ",(hamiltonian_eigenpair(i),  &
                     i=1,hamiltonian_no_eigenpairs),"."
      end if
      !
   end subroutine relci_print_summary
   !
   !
   subroutine relci_print_weights(stream,asf_set)
   !--------------------------------------------------------------------
   ! Print  the  weights of the largest five CSF contributors to each ASF.
   !--------------------------------------------------------------------
      !
      integer, intent(in)            :: stream
      type(asf_basis), intent(in)    :: asf_set
      !
      integer                        :: i, j, level, nocsf_min
      integer, dimension(1:10)       :: iw
      real(kind=dp)                  :: wa, wb
      real(kind=dp), dimension(1:10) :: weights
      !
      write(stream,*)
      write(stream,*)
      write(stream,*) "Weights of major contributors to ASF:"
      write(stream,*)
      write(stream,*) " Level  J Parity      CSF contributions"
      write(stream,*)
      !
      do  level = 1,asf_set%noasf
         weights(1:10) = zero;   iw(1:10) = 0
         wb = zero
         do  i = 1,asf_set%csf_set%nocsf
            wa = asf_set%asf(level)%eigenvector(i) * &
                 asf_set%asf(level)%eigenvector(i)
            wb = wb + asf_set%asf(level)%eigenvector(i) * &
                      asf_set%asf(1)%eigenvector(i)
            do  j = 1,9
               if (wa > weights(j)) then
                  weights(j+1:10) = weights(j:9)
                  weights(j)      = wa
                  iw(j+1:10)      = iw(j:9)
                  iw(j)           = i
                  exit
               end if
            end do
         end do
         if (level > 1   .and.   abs(wb) > 0.0001) then
            print *, "level, wb = ",level,wb
            stop "relci_print_weights(): program stop A." 
         end if
         !
         nocsf_min = 5
         do  j = 1,5
            if (weights(j) == zero) then
               nocsf_min = j - 1
               exit
            end if
         end do
         !
         write(stream,1) asf_set%asf(level)%level_No,                       &
              trim(angular_momentum_string(1*asf_set%asf(level)%totalJ,4)), &
                   asf_set%asf(level)%parity,(weights(j),iw(j),j=1,nocsf_min) 
      end do
    1 format(1x,i4,3x,2a4,5(3x,f8.5," of",i5))
      !
   end subroutine relci_print_weights
   !
   !
   function relci_qed_F(n,kappa,Z)                             result(F)
   !--------------------------------------------------------------------
   ! Estimates the function  F (Z*\alpha) by using an interpolation of 
   ! tabulated data from Mohr (1983) or Mohr and Kim (1992).
   !
   ! Calls: relci_qed_F_Mohr(), relci_qed_F_Mohr_Kim().
   !--------------------------------------------------------------------
      !
      integer, intent(in)               :: n, kappa
      real(kind=dp), intent(in)         :: Z
      real(kind=dp)                     :: F
      !
      if (n <= 2) then
         F = relci_qed_F_Mohr(n,kappa,z)
      else if (3 <= n   .and.   n <= 7) then
         !! F = relci_qed_F_Mohr_Kim(n,kappa,z)
         F = relci_qed_F_Klarsfeld(n,kappa,z)
      else
         F = zero
      end if
      !
   end function relci_qed_F
   !
   !
   function relci_qed_F_Klarsfeld(n,kappa,Z)                   result(F)
   !--------------------------------------------------------------------
   ! Estimates the function  F (Z*\alpha) by using a series expansion
   ! from S Klarsfeld and A Maquet, Physics Letters  43B (1973) 201,
   ! Eqs (1) and (2) and the table of Bethe logarithms. The 
   ! vacuum-polarization contribution in Eq (2) is omitted. 
   ! This procedure is adapted from RCI92 of GRASP2K, written
   ! by Farid A Parpia, to the Fortran 95 standard. 
   !
   ! This procedure is not used in the current version.
   !--------------------------------------------------------------------
      !
      integer, intent(in)               :: n, kappa
      real(kind=dp), intent(in)         :: Z
      real(kind=dp)                     :: F
      !
      real(kind=dp), dimension(36), parameter :: bethe = &
         (/ 2.9841285_dp,   2.8117699_dp,  -0.0300167_dp,   2.7676636_dp, &
           -0.0381902_dp,  -0.0052321_dp,   2.7498118_dp,  -0.0419549_dp, &
           -0.0067409_dp,  -0.0017337_dp,   2.7408237_dp,  -0.0440347_dp, &
           -0.0076008_dp,  -0.0022022_dp,  -0.0007721_dp,   2.7356642_dp, &
           -0.0453122_dp,  -0.0081472_dp,  -0.0025022_dp,  -0.0009628_dp, &
           -0.0004079_dp,   2.7324291_dp,  -0.0461552_dp,  -0.0085192_dp, &
           -0.0027091_dp,  -0.0010945_dp,  -0.0004997_dp,  -0.0002409_dp, &
            2.7302673_dp,  -0.0467413_dp,  -0.0087850_dp,  -0.0028591_dp, &
           -0.0011904_dp,  -0.0005665_dp,  -0.0002904_dp,  -0.0001539_dp /) 
      !
      real(kind=dp), parameter :: C401 = 11.0_dp/24.0_dp,                 &
                               C402 = 3.0_dp/8.0_dp, ovlfac = 4.0_dp/3.0_dp
      !
      integer       :: l, loc
      real(kind=dp) :: bethel, factor, term
      !
      ! Ensure that the principal quantum number is in range
      if (n < 1   .or.   n > 8) then
         print *, "Principal quantum number,",n,", should be in the range 1-8."
         stop     "relci_qed_F_Klarsfeld(): program stop A."
      end if
      !
      l = angular_momentum_l(kappa)
      if (l > n-1) then
         print *, "Kappa = ",kappa," is out of range for n = ",n,"."
         stop     "relci_qed_F_Klarsfeld(): program stop B."
      end if
      !
      ! Find the appropriate entry in the table
      loc    = (n*n-n)/2+l+1
      bethel = bethe(loc)
      !
      ! Determine the quantity in square brackets in eq.(1) of
      ! Klarsfeld and Maquet
      term = -bethel
      !
      if (kappa > 0) then
         term = term - c402 / (l*(l+l+one))
      else
         term = term + c402 / ((l+one)*(l+l+one))
         if (kappa == -1) then
            factor = log (Z/c)
            factor = - (factor + factor)
            term   = term + factor + c401
         end if
      end if
      !
      F = ovlfac * term
      !
   end function relci_qed_F_Klarsfeld
   !
   !
   function relci_qed_F_Mohr(n,kappa,Z)                        result(F)
   !--------------------------------------------------------------------
   ! Computes the function  F (Z*\alpha) for the  1s  2s  2p-  2p  
   ! symmetries by interpolating in, or extrapolating from, the table 
   ! due to  P J Mohr. See  P J Mohr, At Data Nucl Data Tables 29 
   ! (1983) 453. 
   ! This procedure is adapted from RCI92 of GRASP2K, written
   ! by Farid A Parpia, to the Fortran 95 standard.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in)       :: n, kappa
      real(kind=dp), intent(in) :: Z
      real(kind=dp)             :: F, value
      !
      ! Number of data points
      integer, parameter        :: numval = 12
      real(kind=dp), parameter  :: accy   = 1.0e-3_dp
      !
      ! 1s data
      real(kind=dp), dimension(numval), parameter :: val1s =     &
         (/   10.3168_dp,  4.6540_dp,   3.2460_dp,   2.5519_dp,  &
               2.1351_dp,  1.8644_dp,   1.6838_dp,   1.5675_dp,  &
               1.5032_dp,  1.4880_dp,   1.5317_dp,   1.6614_dp  /)
      ! 2s data
      real(kind=dp), dimension(numval), parameter :: val2s =     &
         (/   10.5468_dp,  4.8930_dp,   3.5063_dp,   2.8391_dp,  &
               2.4550_dp,  2.2244_dp,   2.0948_dp,   2.0435_dp,  &
               2.0650_dp,  2.1690_dp,   2.3870_dp,   2.7980_dp  /)
      ! 2p- data
      real(kind=dp), dimension(numval), parameter :: val2p1 =    &
         (/   -0.1264_dp, -0.1145_dp,  -0.0922_dp,  -0.0641_dp,  &
              -0.0308_dp,  0.0082_dp,   0.0549_dp,   0.1129_dp,  &
               0.1884_dp,  0.2934_dp,   0.4530_dp,   0.7250_dp  /)
      ! 2p data
      real(kind=dp), dimension(numval), parameter :: val2p3 =    &
         (/    0.1235_dp,  0.1303_dp,   0.1436_dp,   0.1604_dp,  &
               0.1794_dp,  0.1999_dp,   0.2215_dp,   0.2440_dp,  &
               0.2671_dp,  0.2906_dp,   0.3141_dp,   0.3367_dp  /)
      ! Z data
      real(kind=dp), dimension(numval), parameter :: arg =       &
         (/    1.0_dp,    10.0_dp,     20.0_dp,     30.0_dp,     &
              40.0_dp,    50.0_dp,     60.0_dp,     70.0_dp,     &
              80.0_dp,    90.0_dp,    100.0_dp,    110.0_dp     /)
      !
      ! Interpolate or issue error message as appropriate
      if (n == 1) then
         select case(kappa)
         case(-1)
            call interpolation_aitken(arg,val1s,numval,z,value,accy)
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr(): program stop A."
         end select
      else if (n == 2) then
         select case(kappa)
         case(-1)
            call interpolation_aitken(arg,val2s, numval,z,value,accy)
         case(1)
            call interpolation_aitken(arg,val2p1,numval,z,value,accy)
         case(-2)
            call interpolation_aitken(arg,val2p3,numval,z,value,accy)
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr(): program stop B."
         end select
      else
         print *, "Principal quantum number, ",n,"should be either 1 or 2."
         stop     "relci_qed_F_Mohr(): program stop B."
      end if
      !
      F = value
      !
   end function relci_qed_F_Mohr
   !
   !
   function relci_qed_F_Mohr_Kim(n,kappa,Z)                    result(F)
   !--------------------------------------------------------------------
   ! Computes the function  F (Z*\alpha) for the  1s  2s  2p-  2p  
   ! symmetries by interpolating in, or extrapolating from, the table 
   ! due to  P J Mohr and Y-K Kim. See P J Mohr and Y-K Kim, 
   ! Phys. Rev. A45 (1992) 2723.
   !
   ! Since no values are given for Z = 1, these values were estimated 
   ! below by (graphical) extrapolation.
   !
   ! Calls: interpolation_aitken().
   !--------------------------------------------------------------------
      !
      integer, intent(in)               :: n, kappa
      real(kind=dp), intent(in)         :: Z
      real(kind=dp)                     :: F, value
      !
      ! Number of data points
      integer, parameter       :: numval = 12
      real(kind=dp), parameter :: accy   = 1.0e-3_dp
      !
      ! 3s data
      real(kind=dp), dimension(numval), parameter :: val3s =     &
         (/   10.5000_dp,  4.9524_dp,   3.5633_dp,   2.8940_dp,  &
               2.5083_dp,  2.2757_dp,   2.1431_dp,   2.0874_dp,  &
               2.1018_dp,  2.1935_dp,   2.3897_dp,   2.7609_dp  /)
      ! 3p- data
      real(kind=dp), dimension(numval), parameter :: val3p1 =    &
         (/   -0.1300_dp, -0.1021_dp,  -0.0760_dp,  -0.0430_dp,  &
               0.0041_dp,  0.0414_dp,   0.0956_dp,   0.1623_dp,  &
               0.2483_dp,  0.3660_dp,   0.5408_dp,   0.8322_dp  /)
      ! 3p data
      real(kind=dp), dimension(numval), parameter :: val3p3 =    &
         (/    0.1250_dp,  0.1421_dp,   0.1572_dp,   0.1761_dp,  &
               0.1977_dp,  0.2214_dp,   0.2470_dp,   0.2745_dp,  &
               0.3038_dp,  0.3350_dp,   0.3679_dp,   0.4020_dp  /)
      ! 3d- data
      real(kind=dp), dimension(numval), parameter :: val3d3 =    &
         (/   -0.0440_dp, -0.0428_dp,  -0.0420_dp,  -0.0410_dp,  &
              -0.0396_dp, -0.0378_dp,  -0.0353_dp,  -0.0321_dp,  &
              -0.0279_dp, -0.0225_dp,  -0.0154_dp,  -0.0062_dp  /)
      !
      ! 4s data
      real(kind=dp), dimension(numval), parameter :: val4s =     &
         (/   10.5000_dp,  4.9749_dp,   3.5834_dp,   2.9110_dp,  &
               2.5215_dp,  2.2842_dp,   2.1455_dp,   2.0814_dp,  &
               2.0840_dp,  2.1582_dp,   2.3262_dp,   2.6484_dp  /)
      ! 4p- data
      real(kind=dp), dimension(numval), parameter :: val4p1 =    &
         (/   -0.1200_dp, -0.0963_dp,  -0.0690_dp,  -0.0344_dp,  &
               0.0064_dp,  0.0538_dp,   0.1098_dp,   0.1780_dp,  &
               0.2649_dp,  0.3819_dp,   0.5525_dp,   0.8311_dp  /)
      ! 4p data
      real(kind=dp), dimension(numval), parameter :: val4p3 =    &
         (/    0.1250_dp,  0.1477_dp,   0.1630_dp,   0.1827_dp,  &
               0.2052_dp,  0.2299_dp,   0.2568_dp,   0.2858_dp,  &
               0.3170_dp,  0.3507_dp,   0.3868_dp,   0.4247_dp  /)
      ! 4d- data
      real(kind=dp), dimension(numval), parameter :: val4d3 =    &
         (/   -0.0410_dp, -0.0403_dp,  -0.0399_dp,  -0.0387_dp,  &
              -0.0371_dp, -0.0348_dp,  -0.0317_dp,  -0.0276_dp,  &
              -0.0222_dp, -0.0149_dp,  -0.0053_dp,   0.0074_dp  /)
      !
      ! 5s data
      real(kind=dp), dimension(numval), parameter :: val5s =     &
         (/   10.5000_dp,  4.9858_dp,   3.5923_dp,   2.9173_dp,  &
               2.5246_dp,  2.2833_dp,   2.1395_dp,   2.0686_dp,  &
               2.0619_dp,  2.1225_dp,   2.2696_dp,   2.5566_dp  /)
      ! 5p- data
      real(kind=dp), dimension(numval), parameter :: val5p1 =    &
         (/   -0.1200_dp, -0.0933_dp,  -0.0652_dp,  -0.0299_dp,  &
               0.0116_dp,  0.0597_dp,   0.1161_dp,   0.1843_dp,  &
               0.2703_dp,  0.3848_dp,   0.5497_dp,   0.8150_dp  /)
      ! 5p data
      real(kind=dp), dimension(numval), parameter :: val5p3 =    &
         (/    0.1300_dp,  0.1502_dp,   0.1662_dp,   0.1861_dp,  &
               0.2089_dp,  0.2341_dp,   0.2614_dp,   0.2910_dp,  &
               0.3229_dp,  0.3574_dp,   0.3946_dp,   0.4338_dp  /)
      ! 5d- data
      real(kind=dp), dimension(numval), parameter :: val5d3 =    &
         (/   -0.0405_dp, -0.0396_dp,  -0.0387_dp,  -0.0374_dp,  &
              -0.0356_dp, -0.0331_dp,  -0.0297_dp,  -0.0252_dp,  &
              -0.0190_dp, -0.0108_dp,   0.0001_dp,   0.0145_dp  /)
      ! Z data
      real(kind=dp), dimension(numval), parameter :: arg =       &
         (/    1.0_dp,    10.0_dp,     20.0_dp,     30.0_dp,     &
              40.0_dp,    50.0_dp,     60.0_dp,     70.0_dp,     &
              80.0_dp,    90.0_dp,    100.0_dp,    110.0_dp     /)
      !
      ! Interpolate or issue error message as appropriate
      if (n == 3) then
         select case(kappa)
         case(-1)
            call interpolation_aitken(arg,val3s,numval,z,value,accy)
         case(1)
            call interpolation_aitken(arg,val3p1,numval,z,value,accy)
         case(-2)
            call interpolation_aitken(arg,val3p3,numval,z,value,accy)
         case(2)
            call interpolation_aitken(arg,val3d3,numval,z,value,accy)
         case(-3)
            call interpolation_aitken(arg,val3d3,numval,z,value,accy)
            value = 0.9_dp * value
            value = zero
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr_Kim(): program stop A."
         end select
      else if (n == 4) then
         select case(kappa)
         case(-1)
            call interpolation_aitken(arg,val4s,numval,z,value,accy)
         case(1)
            call interpolation_aitken(arg,val4p1,numval,z,value,accy)
         case(-2)
            call interpolation_aitken(arg,val4p3,numval,z,value,accy)
         case(2)
            call interpolation_aitken(arg,val4d3,numval,z,value,accy)
         case(-3)
            call interpolation_aitken(arg,val4d3,numval,z,value,accy)
            value = 0.9_dp * value
            value = zero
         case(3)
            value = zero
         case(-4)
            value = zero
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr_Kim(): program stop B."
         end select
      else if (5 <= n  .and.   n <= 7) then
         select case(kappa)
         case(-1)
            call interpolation_aitken(arg,val5s,numval,z,value,accy)
         case(1)
            call interpolation_aitken(arg,val5p1,numval,z,value,accy)
         case(-2)
            call interpolation_aitken(arg,val5p3,numval,z,value,accy)
         case(2)
            call interpolation_aitken(arg,val5d3,numval,z,value,accy)
         case(-3)
            call interpolation_aitken(arg,val4d3,numval,z,value,accy)
            value = 0.9_dp * value
            value = zero
         case(3:6)
            value = zero
         case(-7:-4)
            value = zero
         case default
            print *, "Principal quantum number is ",n," and kappa",kappa
            stop     "relci_qed_F_Mohr_Kim(): program stop C."
         end select
      else
         print *, "Principal quantum number, ",n,"should be in the interval "//&
                  "3 <= n <= 7."
         stop     "relci_qed_F_Mohr(): program stop D."
      end if
      !
      F = value
      !
   end function relci_qed_F_Mohr_Kim
   !
   !
   subroutine relci_selfenergy(asf_set,wave,shift,lprint)
   !--------------------------------------------------------------------
   ! Calculates the QED estimations for all selected ASF from a hydrogen-
   ! like model
   !--------------------------------------------------------------------
      !
      type(asf_basis), intent(in)              :: asf_set
      type(grasp2k_orbital), intent(in)        :: wave
      real(kind=dp), dimension(:), intent(out) :: shift
      logical, intent(in)                      :: lprint
      !
      integer       :: i, j
      real(kind=dp) :: F_alphaZ, se, ratio
      real(kind=dp), dimension(1:asf_set%csf_set%nwshells) :: se_shell
      !
      se_shell(:) = zero
      shift(:)    = zero
      !!x print *, "*** QED contributions in meV ***"
      !!x print *, "*** QED contributions in meV ***"
      !!x print *, "*** QED contributions in meV ***"
      do  i = 1,asf_set%csf_set%nwshells
         do  j = 1,wave%number_of_rwf
            if (asf_set%csf_set%subshell(i)%n == wave%rwf(j)%orbital%n   .and. &
               asf_set%csf_set%subshell(i)%kappa == wave%rwf(j)%orbital%kappa) &
               goto 1
         end do
         stop "relci_selfenergy(): program stop A."
         !
       1 ratio       = selfenergy_ratio_grasp2k(wave%rwf(j))
         F_alphaZ    = relci_qed_F(wave%rwf(j)%orbital%n,           &
                                   wave%rwf(j)%orbital%kappa,nuclear_charge)
         se_shell(i) = ratio * F_alphaZ * nuclear_charge**4 / (pi*c**3) /     &
                                    (asf_set%csf_set%subshell(i)%n**3 *one) 
         if (lprint) then
            print 2,                                                          &
               orbital_name(wave%rwf(j)%orbital%n,wave%rwf(j)%orbital%kappa), &
               se_shell(i), ratio, F_alphaZ
               !!x se_shell(i)*27210., ratio, F_alphaZ
          2 format(5x,a,11x,e11.4,17x,f7.3,11x,e14.7)
         end if
      end do
      !
      do  j = 1,asf_set%csf_set%nocsf
         se = zero
         do  i = 1,asf_set%csf_set%nwshells
            se = se + asf_set%csf_set%csf(j)%occupation(i)*se_shell(i)
         end do
         !
         do  i = 1,asf_set%noasf
            shift(i) = shift(i) + asf_set%asf(i)%eigenvector(j)* &
                                  asf_set%asf(i)%eigenvector(j)* se
         end do
      end do
      !
   end subroutine relci_selfenergy
   !
   !
   subroutine relci_set_eigenpairs(asf_set)
   !--------------------------------------------------------------------
   ! Set the selected eigenpairs from the energies obtained for the
   ! individual (J,parity) blocks of the Hamiltonian matrix.
   !--------------------------------------------------------------------
      !   
      type(asf_basis), intent(inout) :: asf_set
      !
      integer       :: ce, iblock, j, k, no_eigenvectors, s
      real(kind=dp) :: current, last
      real(kind=dp), dimension(1:asf_set%csf_set%nocsf) :: eigenvector
      !
      ce = 0  ! current eigenpair
      !
      last    = -1.0e99_dp
      !
    1 current =  1.0e99_dp
      !
      do  iblock = 1,hamiltonian_noblock 
         no_eigenvectors = min(hamiltonian_no_eigenpairs,H_block(iblock)%nocsf)
         do  s = 1,no_eigenvectors
            if (H_block(iblock)%eigenvalue(s) > last  .and.  &
                H_block(iblock)%eigenvalue(s) < current)  then
               current = H_block(iblock)%eigenvalue(s)
            end if
         end do
      end do
      !
      do  iblock = 1,hamiltonian_noblock 
         no_eigenvectors = min(hamiltonian_no_eigenpairs,H_block(iblock)%nocsf)
         do  s = 1,no_eigenvectors
            if (H_block(iblock)%eigenvalue(s) == current) then
               ce = ce + 1
               do  k = 1,asf_set%noasf
                  if (ce == asf_set%asf(k)%level_No) then
                     asf_set%asf(k)%totalJ = H_block(iblock)%totalJ
                     asf_set%asf(k)%parity = H_block(iblock)%parity
                     asf_set%asf(k)%energy = H_block(iblock)%eigenvalue(s)
                     eigenvector(:)        = zero
                     do  j = 1,H_block(iblock)%nocsf
                        eigenvector(H_block(iblock)%csf_ndx(j)) =   &
                                    H_block(iblock)%eigenvector(s,j)
                     end do
                     asf_set%asf(k)%eigenvector(:) = eigenvector(:)
                  end if 
               end do
            end if
         end do
      end do
      !
      last = current
      if (ce < hamiltonian_max_eigenpair) goto 1
      !
   end subroutine relci_set_eigenpairs
   !
end module rabs_relci
