module rabs_cesd
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains a number of procedures which are specific to the
! CESD program, i.e. to the complete expansion of a list of symmetry
! functions into a determinant basis. It performs this expansion by
! decomposing the individual antisymmetric subshell states into uncoupled
! product functions. Several procedures also deal with the correct file
! handling for the CESD program; a short list of all procedures is given
! below.
!-----------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   use rabs_csl
   use rabs_determinant
   use rabs_dirac_orbital
   use rabs_file_handling
   !
   private :: append_determinants
                 ! 'Appends' the determinant expansion of a given CSF i to the
                 ! previously obtained expansion.
   public  :: cesd_collect_input
                 ! Collects and proceeds all input for the expansion of the
                 ! symmetry-adapted functions into determinants.
   public  :: cesd_expand_csf
                 ! Returns the determinant expansion of a single CSF.
   public  :: cesd_perform_xpansion
                 ! Carries out the complete expansion of the atomic and/or
                 ! configuration state functions into a determinant basis.
   public  :: cesd_perform_xpansion_all
                 ! Carries out the complete expansion of the atomic and/or
                 ! configuration state functions into a determinant basis.
                 ! All possible M substates are taken into account.
   private :: determine_weight_coefficients
                 ! Determines the 'weights' of a given CSF in a (also given)
                 ! determinant basis.
   private :: expand_csf_into_determinants
                 ! Returns the determinant expansion of a single CSF.
   private :: expand_parent_states
                 ! Calculates the expansion coefficients for the decomposition 
                 ! of a coupled antisymmetric parent state.
   private :: expand_subshell_state
                 ! Carries out the expansion of a coupled antisymmetric 
                 ! subshell state for an open shell.
   private :: grant_ipoidw
                 ! Auxiliarity routine which provides 'pointer' to the
                 ! function grant_npoidw.
   private :: grant_npoidw
                 ! Returns subshell state qantum numbers for a proper choice
                 ! of integer parameters; replaces some 'block data' definition
                 ! from GRASP.
   private :: index_combination
                 ! Calculates the index combination (k1,...,k9) for lists of
                 ! subshell determinants for up to nine open shells.
   private :: select_determinant
                 ! Determines whether a list of subshell determinants may
                 ! be combined to a total determinant with proper M quantum
                 ! number.    
   private :: set_orbital_reference
                 ! Generates an orbital reference list of (n, kappa, m) 
                 ! orbitals for a determinant basis.
   private :: write_expansion_to_file
                 ! Writes the (newly calculated) determinant basis to 
                 ! the CESD output file.
   private :: write_expansion_to_file_compact
                 ! Writes the (newly calculated) determinant basis to 
                 ! the CESD output file in the new compact format.
   !
   ! Define some global input/output flags of the CESD program which can be
   ! set to .true. at input time interactively
   logical, private :: cesd_print_orbital_reference  = .false.,  &
                       cesd_print_csf_expansion      = .false.,  &
                       cesd_print_complete_expansion = .false.
   logical, public  :: cesd_use_formatted_mix_file   = .false.,  &
                       cesd_use_compact_xpn_file     = .true.
   character(len=256), public :: cesd_xpn_file
   !
   integer :: i, maximal_no_determinants
   integer, parameter, private :: max_open_shells = 9
   integer, dimension(0:40), private :: total_MM = (/ (i,i=0,40) /)
   !
   type, private :: permutation
      integer :: nperm
      integer, dimension(:,:), pointer :: perm
      real(kind=dp), dimension(:),  pointer :: weight 
   end type permutation
   type(permutation), dimension(:), pointer, private :: permshell
   !
contains
   !
   subroutine append_determinants(asf_mode,asf_set,asf_new,i,MM, &
                                  det_list,det_weight,nod)
   !--------------------------------------------------------------------
   ! This routine 'appends' the expansion of the CSF i (given as a list
   ! of determinants and weights) to the determinants which were 
   ! already included in the expansion of some previous CSF to built up 
   ! the full expansion of some given ASF.
   ! For asf_mode == .false. (CSF-based), the expansion is purely copied 
   ! to store it for later output. 
   ! For asf_mode == .true.  (ASF-based), the configuration mixing of the
   ! ASF in the CSF basis is taken into account and the contribution 
   ! are added for all selected atomic states. 
   !--------------------------------------------------------------------
      !
      implicit none
      logical, intent(in) :: asf_mode
      integer, intent(in) :: i, MM, nod
      type(asf_basis), intent(in)        :: asf_set
      type(asf_det_basis), intent(inout) :: asf_new
      type(determinant), dimension(1:nod), intent(in) :: det_list
      real(kind=dp), dimension(1:nod), intent(in)     :: det_weight
      !
      integer       :: int, j, j_sav, k, kasf, kasf_old, k_new, k_old, &
                       l, l_sav, noint, p
      real(kind=dp) :: sum
      type(determinant), dimension(:), pointer :: determinant_sav
      type(as_function), dimension(:), pointer :: asf_sav
      !
      noint = asf_new%det_set%noint
      !
      ! Add or copy the contribution of each determinant into asf_new;
      ! first check whether the determinant is already defined therein
      sum = zero
      do  k = 1,nod
         if (det_list(k)%parity /= '+' .and. det_list(k)%parity /= '-') then
            print *, "i,k,nod,size(det_list(k)%occupation),weight = ", &
                      i,k,nod,size(det_list(k)%occupation(:)),det_weight(k)
         end if
         sum = sum + det_weight(k) * det_weight(k)
      end do
      !
      if (abs(sum - one) > 1.0e-8_dp ) then
         print 1,     nod, sum
         write (24,1) nod, sum
       1 format(//1x,'(Number of determinants, norm) = ',i4,3x,1pd15.8)
         print *,     "append_determinants(): Incomplete decomposition of "// &
                      "the CSF i = ",i
         write (24,*) "append_determinants(): Incomplete decomposition of "// &
                      "the CSF i = ",i
         stop
      end if
      !
      if (mod(i,50) == 0) then
         print *, I,"CSF complete ...  (with maximal",                     &
                  maximal_no_determinants,"determinants in a single CSF;", &
                  asf_new%det_set%nod,"total No. of determinants)"
      end if
      !
      do  k = 1,nod
         if (abs(det_weight(k)) < eps10) cycle
         do  k_old = asf_new%det_set%nod,1,-1
            do  int = 1,asf_new%det_set%noint
               if (det_list(k)%occupation(int) /= &
                   asf_new%det_set%determinant(k_old)%occupation(int)) then
                  goto 2
               end if
            end do 
            !
            if (asf_mode) then
               do  kasf = 1,asf_new%noasf
                  if (asf_new%asf(kasf)%totalM == MM) then
                  do  j = 1,asf_set%noasf
                  if (asf_new%asf(kasf)%level_No==asf_set%asf(j)%level_No) then
                     kasf_old = j
                     exit
                  end if
                  end do
                  asf_new%asf(kasf)%eigenvector(k_old) =    &
                     asf_new%asf(kasf)%eigenvector(k_old) + &
                     asf_set%asf(kasf_old)%eigenvector(i) * det_weight(k)
                  end if
               end do 
            else 
               asf_new%asf(i)%eigenvector(k_old) =  det_weight(k) 
            end if
            goto 3
       2    continue
         end do
         !
         ! Append a "new" determinant to the total representation;
         ! first check if there is still enough memory available, re-allocate
         ! fresh storage otherwise
         !
         if (asf_new%det_set%nod + 1 > asf_new%det_set%nodmax) then
            ! 
            ! Save and expand the list of determinants, i.e. the array
            ! asf_new%det_set%determinant
            j_sav = asf_new%det_set%nodmax
            print *, "Re-allocate memory for the determinant basis ..."
            print *, "   current number of determinants is ",j_sav
            !
            allocate( determinant_sav(1:asf_new%det_set%nodmax) )
            do  j = 1,j_sav
              allocate( determinant_sav(j)%occupation(1:noint) )
              determinant_sav(j)%totalM = asf_new%det_set%determinant(j)%totalM
              determinant_sav(j)%parity = asf_new%det_set%determinant(j)%parity
              determinant_sav(j)%occupation(1:noint) = &
                         asf_new%det_set%determinant(j)%occupation(1:noint)
              deallocate( asf_new%det_set%determinant(j)%occupation )
            end do
            deallocate( asf_new%det_set%determinant )
            !
            asf_new%det_set%nodmax = asf_new%det_set%nodmax + &
                                     asf_set%csf_set%nocsf
            print *, "   new number of determinants is     ", &
                         asf_new%det_set%nodmax
            allocate( asf_new%det_set%determinant(1:asf_new%det_set%nodmax) )
            do  j = 1,j_sav
              allocate( asf_new%det_set%determinant(j)%occupation(1:noint) )
              asf_new%det_set%determinant(j)%totalM = determinant_sav(j)%totalM
              asf_new%det_set%determinant(j)%parity = determinant_sav(j)%parity
              asf_new%det_set%determinant(j)%occupation(1:noint) = &
                                         determinant_sav(j)%occupation(1:noint)
              deallocate( determinant_sav(j)%occupation )
            end do
            deallocate( determinant_sav )
            !
            ! Save and expand the corresponding set of eigenvectors
            l_sav = asf_new%noasf
            allocate( asf_sav(1:l_sav) )
            do  l = 1,l_sav
               allocate( asf_sav(l)%eigenvector(1:j_sav) )
               asf_sav(l)%level_No = asf_new%asf(l)%level_No
               asf_sav(l)%totalJ   = asf_new%asf(l)%totalJ
               asf_sav(l)%parity   = asf_new%asf(l)%parity
               asf_sav(l)%energy   = asf_new%asf(l)%energy
               asf_sav(l)%eigenvector(1:j_sav) = &
                          asf_new%asf(l)%eigenvector(1:j_sav)
               deallocate( asf_new%asf(l)%eigenvector )
               allocate( asf_new%asf(l)%eigenvector(1:asf_new%det_set%nodmax) )
               asf_new%asf(l)%level_No = asf_sav(l)%level_No
               asf_new%asf(l)%totalJ   = asf_sav(l)%totalJ
               asf_new%asf(l)%parity   = asf_sav(l)%parity
               asf_new%asf(l)%energy   = asf_sav(l)%energy
               asf_new%asf(l)%eigenvector(1:j_sav) = &
                          asf_sav(l)%eigenvector(1:j_sav)
               deallocate( asf_sav(l)%eigenvector )
               do  p = j_sav+1,asf_new%det_set%nodmax
                  asf_new%asf(l)%eigenvector(p) = zero
               end do
            end do
            deallocate( asf_sav )
            !
            print *, "   ... re-allocation complete."
         end if
         !
         ! Append one determinant
         asf_new%det_set%nod = asf_new%det_set%nod + 1
         k_new               = asf_new%det_set%nod
         allocate( asf_new%det_set%determinant(k_new)%occupation(1:noint) )
         asf_new%det_set%determinant(k_new)%occupation(1:noint) =  &
                                det_list(k)%occupation(1:noint)
         asf_new%det_set%determinant(k_new)%totalM = det_list(k)%totalM
         asf_new%det_set%determinant(k_new)%parity = det_list(k)%parity
         !
         if (asf_mode) then
            do  kasf = 1,asf_new%noasf
               if (asf_new%asf(kasf)%totalM == MM) then
               do  j = 1,asf_set%noasf
                  if (asf_new%asf(kasf)%level_No==asf_set%asf(j)%level_No) then
                     kasf_old = j
                     exit
                  end if
               end do
               asf_new%asf(kasf)%eigenvector(k_new) =    &
                  asf_new%asf(kasf)%eigenvector(k_new) + &
                  asf_set%asf(kasf_old)%eigenvector(i) * det_weight(k)
               end if
            end do 
         else 
            !x print *, "pass F; k_new =  ",k_new
            asf_new%asf(i)%eigenvector(k_new) =  det_weight(k) 
            !x print *, "pass G; k_new =  ",k_new
         end if
    3    continue
         !x print *, "pass H"
      end do
      !
   end subroutine append_determinants
   !
   !
   subroutine cesd_collect_input(asf_mode)
   !--------------------------------------------------------------------
   ! Collect and proceeds all input for the closedshell program.
   !
   ! Calls: get_integer_from_string(), get_yes_stream().
   !--------------------------------------------------------------------
      !
      implicit none
      logical, intent(inout) :: asf_mode
      !
      integer           :: jjx, mmx, pos
      logical           :: yes
      character(len=20) :: string
      !
      ! Select the mode of the cesd calculation
    1 print *, "Expand the atomic state functions (ASF) ?"
      yes = get_yes_stream()
      if (yes) then
         asf_mode = .true.
      else
         print *, "Expand the configuration state functions (CSF) ?"
         yes = get_yes_stream()
         if (yes) then
            asf_mode = .false.
         else
            print *, "You must choose either one, an ASF or CSF expansion !"
            goto 1
         end if
      end if
      !
      ! Print non-standard output if required
      print *, "Print non-standard output on the cesd expansion ?"
      yes = get_yes_stream()
      if (yes) then
         print *, "Print the reference list of the one-electron orbitals"
         print *, " each given by its quantum numbers (n kappa m) ?"
         yes = get_yes_stream()
         if (yes) then
            cesd_print_orbital_reference = .true.
         end if
         print *, "Print the expansion of each CSF each time it is" //  &
                  " carried out ?"
         print *, " This prints the occupation and all weights of the"// & 
                  " determinenants."
         yes = get_yes_stream()
         if (yes) then
            cesd_print_csf_expansion = .true.
         end if
         print *, "Print the 'complete' expansion into Slater"// &
                  " determinants after"
         print *, " all symmetry  functions (either ASF or CSF) have"// & 
                  " been projected ? "
         yes = get_yes_stream()
         if (yes) then
            cesd_print_complete_expansion = .true.
         end if
      end if
      !
      ! Specify total angular momentum quantum number M =/= J for the
      ! Slater determinants in the expansion of the symmetry functions 
      !
      print *, "Select non-standard total M values for the "// &
               " expansion into Slater determinants ? "
      print *, " The standard is to use M = +J for each individual"// &
               " determinant."
      yes = get_yes_stream()
      if (yes) then
    2    print *, "Enter a pair of quantum numbers (2*J, 2*M);"
         print *, " half-integers are not allowd here (null if done)."
         read  (*,"(a)") string
         if (len_trim(string) /= 0) then
            string = adjustl(string)
            pos    = scan(string," ") - 1
            jjx    = get_integer_from_string(string(1:pos))
            pos    = pos + 1
            string = adjustl(string(pos:20))
            mmx    = get_integer_from_string(string(1:20))
            if ((mod(jjx,2) == 0   .and.   mod(mmx,2) == 1)   .or. &
                (mod(jjx,2) == 1   .and.   mod(mmx,2) == 0)   .or. &
                jjx < 0           .or.  jjx < abs(mmx))  then
               print *, "2*J, 2*M = ",jjx, mmx
               print *, "Enter a valid pair (2*J,2*M) of total angular"// &
                        " momenta."
               goto 2
            else
               total_MM(jjx) = mmx
               goto 2
            end if
         end if
      end if
      !
   end subroutine cesd_collect_input
   !
   !
   subroutine cesd_expand_csf(csf_set,det_set,i,MM,det_occupation, &
                                                   det_weight,nod)
   !--------------------------------------------------------------------
   ! Expands the configuration state function i from csf_set [which is of
   ! type(csf_basis)] into a series of Slater determinants with the 
   ! total magnetic projection MM. It returns a list of nod 
   ! determinants and weights.
   !
   ! Calls: distribute_m_in_ncells(), pack_occupation_in_integer().
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)         :: i, MM
      integer, intent(out)        :: nod
      type(csf_basis), intent(in) :: csf_set
      type(det_basis), intent(in) :: det_set
      integer, dimension(:,:), pointer         :: det_occupation
      real(kind=dp),	 dimension(:), pointer :: det_weight
      !
      logical :: append
      integer :: d, dmax, j, nelec, nelmax, idet, ndet, noint, np, nq, &
                 no_determinants, no_open_shells, shell
      integer, dimension(1:max_open_shells) :: open_index
      integer, dimension(:), pointer        :: det_occ
      integer, dimension(:,:), pointer      :: perm_work
      real(kind=dp)                         :: unc_weight
      !
      ! Allocate memory to create a list of all possible 'determinants'
      ! for CSF i with given number of open shells;
      ! For all open shells, also find all distribution of the electrons
      ! among the (n, kappa, m) orbitals of these subshells
      !
      allocate( permshell(1:csf_set%nwshells) )
      do  shell = 1,csf_set%nwshells
         nelmax = angular_momentum_j(csf_set%subshell(shell)%kappa) + 1
         nelec  = csf_set%csf(i)%occupation(shell)
         if (nelec == 0) then
            permshell(shell)%nperm = 0
         else if (nelec == nelmax) then
            permshell(shell)%nperm = 1
         else if (0 < nelec  .and.   nelec < nelmax) then
            np = 1;   do  j = nelmax-nelec+1,nelmax;  np = np * j;   end do
            nq = 1;   do  j = 1,nelec;                nq = nq * j;   end do
            permshell(shell)%nperm     = np / nq
            allocate( permshell(shell)%perm(nelmax,permshell(shell)%nperm), &
                      perm_work(nelmax,permshell(shell)%nperm),             &
                      permshell(shell)%weight(permshell(shell)%nperm) )
            do  j = 1,permshell(shell)%nperm  
               permshell(shell)%weight(j) = zero
            end do 
            dmax = permshell(shell)%nperm
            call distribute_m_in_ncells(nelec,nelmax,d,dmax,perm_work)
            permshell(shell)%perm = perm_work
            if (rabs_use_stop   .and.  permshell(shell)%nperm /= d) then
               print *, "d, nperm = ", &
                         d, permshell(shell)%nperm
               stop "expand_csf_into_determinants(): program stop nperm."
            end if
            deallocate( perm_work )
         else if (rabs_use_stop) then
            print *, "i, shell-kappa, nwshells, nelec, nelmax = ", i, &
                     csf_set%subshell(shell)%kappa, csf_set%nwshells, &
                     nelec, nelmax
            stop "expand_csf_into_determinants(): program stop A."
         end if
      end do
      !
      ! Determine the total number of Slater determinants which,
      ! however, might be larger than the number of necessary determinants
      ! because of the total M projection.
      ndet = 1;   no_open_shells = 0
      do  shell = 1,csf_set%nwshells
         if (permshell(shell)%nperm /= 0) then
            ndet = ndet * permshell(shell)%nperm
         end if
         if (permshell(shell)%nperm /= 0  .and.  &
             permshell(shell)%nperm /= 1) then
            no_open_shells = no_open_shells + 1
            open_index(no_open_shells) = permshell(shell)%nperm
         end if
      end do
      !
      do  j = no_open_shells+1,max_open_shells;  open_index(j) = 1;  end do
      !
      ! Construct the complete list of all possible Slater determinants; 
      ! check the magnetic quantum numbers of the individual shells as 
      ! well as MM; allocate memory for this 'intermediate' list of
      ! determinants
      !
      allocate( det_occupation(1:det_set%norbital,1:ndet), &
        	det_weight(1:ndet), det_occ(1:det_set%norbital) )
      call index_combination(open_index,initialization=.true.)
      !
      no_determinants = 0
      do idet = 1,ndet
         call index_combination(open_index,index=idet)
         call select_determinant(open_index,csf_set,det_set,i,MM, &
                                 append,det_occ,unc_weight)
         if (append) then
            no_determinants = no_determinants + 1
            det_occupation(1:det_set%norbital,no_determinants) =  &
               det_occ(1:det_set%norbital)
            det_weight(no_determinants) = unc_weight
         end if
      end do
      nod = no_determinants
      !
      ! Determine the weight coefficients of the different CSF
      !
      call determine_weight_coefficients(csf_set,det_set,i,MM,      &
                                         det_occupation,det_weight, &
                               det_set%norbital,no_determinants,ndet)
      !
      ! Deallocate all working arrays of this routine
      deallocate( det_occ )
      do  shell = 1,csf_set%nwshells
         if (permshell(shell)%nperm > 1) then
            deallocate( permshell(shell)%perm, permshell(shell)%weight)
         end if
      end do
      deallocate( permshell )
      !
   end subroutine cesd_expand_csf
   !
   !
   subroutine cesd_perform_xpansion(asf_mode,asf_set)
   !--------------------------------------------------------------------
   ! This routine controls the CESD expansion. It calculates a complete 
   ! expansion of the jj-coupled CSF or ASF in a determinants basis. 
   ! For both, a Slater determinant as well as a full determinant basis
   ! abstract data types are defines. This definition also includes an 
   ! (ordered) one-particle orbital list where each orbital is given by
   ! its by its quantum numbers (n  kappa  m). 
   !
   ! For invoking this module two modes are allowed. 
   ! For asf_mode == .false., each CSF of the currently defined 
   ! configuration basis (of GRASP2K) is expanded independently into a 
   ! set of determinants and stored on disk.
   ! For asf_mode == .true., the atomic state functions (ASF) are 
   ! transformed into a representation of determinants, i.e. ASF as
   ! obtained from a self-consistent-field calculation. For some given
   ! (complete) active space this is equivalent to an 'orthonormal' 
   ! transformation of the ASF from an representation in terms of 
   ! jj-coupled symmetry functions into (uncoupled) determinants.  
   !
   ! Calls: append_determinants(), expand_csf_into_determinants(),
   !        set_orbital_reference(), write_expansion_to_file().        
   !--------------------------------------------------------------------
      !
      implicit none
      logical, intent(in)         :: asf_mode
      type(asf_basis), intent(in) :: asf_set
      !
      integer             :: i, iorb, j, MM, nobit, nod, noint, pos
      character(len=3)    :: month
      character(len=8)    :: cdate
      character(len=10)   :: ctime
      type(asf_det_basis) :: asf_new
      type(determinant), dimension(:), pointer :: det_list
      real(kind=dp), dimension(:), pointer     :: weight
      !
      ! Allocte initial storage for the 'new' ASF detrminant basis, 
      ! i.e. for asf_new of type(asf_det_basis)
      ! Also, take over the total-J, parity, and energies;
      !
      asf_new%det_set%nod      = 0
      asf_new%det_set%nodmax   = 4 * asf_set%csf_set%nocsf
      asf_new%det_set%norbital = 0
      asf_new%det_set%noint    = 0
      asf_new%det_set%number_of_electrons = asf_set%csf_set%number_of_electrons
      !
      print *, "Allocate memory for ",asf_new%det_set%nodmax," determinants."
      allocate( asf_new%det_set%determinant(1:asf_new%det_set%nodmax) )
      !
      if (asf_mode) then
         asf_new%noasf          = asf_set%noasf
         asf_new%average_energy = asf_set%average_energy
         allocate( asf_new%asf(1:asf_set%noasf) )
         do  i = 1,asf_set%noasf
            asf_new%asf(i)%level_No = asf_set%asf(i)%level_No
            asf_new%asf(i)%totalJ   = asf_set%asf(i)%totalJ
            asf_new%asf(i)%totalM   = total_MM(asf_set%asf(i)%totalJ)
            asf_new%asf(i)%parity   = asf_set%asf(i)%parity
            asf_new%asf(i)%energy   = asf_set%asf(i)%energy
            allocate( asf_new%asf(i)%eigenvector(1:asf_new%det_set%nodmax) )
            asf_new%asf(i)%eigenvector(1:asf_new%det_set%nodmax) = zero
         end do
      else
         asf_new%noasf          = asf_set%csf_set%nocsf
         asf_new%average_energy = zero
         allocate( asf_new%asf(1:asf_set%noasf) )
         do  i = 1,asf_new%noasf
            asf_new%asf(i)%level_No = i
            asf_new%asf(i)%totalJ   = asf_set%csf_set%csf(i)%totalJ
            asf_new%asf(i)%parity   = asf_set%csf_set%csf(i)%parity
            asf_new%asf(i)%energy   = zero
            allocate( asf_new%asf(i)%eigenvector(1:asf_new%det_set%nodmax) )
            asf_new%asf(i)%eigenvector(1:asf_new%det_set%nodmax) = zero
         end do
      end if
      !
      ! Append a header to the cesd expansion output file on stream 25
      if (asf_mode) then
         write(25,"(a)") "ASF-based CESD output file"
      else
         write(25,"(a)") "CSF-based CESD output file"
      end if
      !x print *, "Header appended !"
      !
      ! Define an (ordered) one-electron orbital reference list which 
      ! includes all orbitals of of the 'old' ASF basis (in asf_set). 
      ! Write this list to the CESD output file on stream 25
      !
      call set_orbital_reference(asf_set,asf_new%det_set)
      !
      write (25,301) 
      write (25,302) asf_set%csf_set%nocsf
      write (25,303) asf_new%noasf
      write (25,304) asf_new%det_set%norbital
      if (cesd_use_compact_xpn_file) then
         write (25,305) ( iorb,asf_new%det_set%orbital(iorb)%n,     &
                               asf_new%det_set%orbital(iorb)%kappa, &
                               asf_new%det_set%orbital(iorb)%mm,    &
                               iorb=1,asf_new%det_set%norbital)
      else
         write (25,306) ( iorb,asf_new%det_set%orbital(iorb)%n,     &
                               asf_new%det_set%orbital(iorb)%kappa, &
                               asf_new%det_set%orbital(iorb)%mm,    &
                               iorb=1,asf_new%det_set%norbital)
      end if
      !
      ! Allocate memory for the occupation of determinants in the 
      ! determinant basis; first determine the number of 'integers' 
      ! which is required for storing this occupation and later initialize
      ! all (bits of) occupation to 0.
      !
      noint = bit_size(noint)
      if (mod(asf_new%det_set%norbital,noint) == 0) then
         noint = asf_new%det_set%norbital/noint
      else
         noint = asf_new%det_set%norbital/noint + 1
      end if
      asf_new%det_set%noint = noint
      do  i = 1,asf_new%det_set%nodmax
         allocate( asf_new%det_set%determinant(i)%occupation(1:noint) )
      end do
      !
      nobit = bit_size(asf_new%det_set%determinant(1)%occupation(1))
      do  i = 1,noint;   do  pos = 0,nobit-1
         asf_new%det_set%determinant(1)%occupation(i) = &
            ibclr(asf_new%det_set%determinant(1)%occupation(i),pos)
      end do;   end do
      !
      do  i = 2,asf_new%det_set%nodmax 
         asf_new%det_set%determinant(i)%occupation(1:noint) = &
                     asf_new%det_set%determinant(1)%occupation(1:noint)
      end do
      !
      ! Now expand the CSF from the original ASF basis one by one into
      ! determinants
      !
      do  i = 1,asf_set%csf_set%nocsf
         MM = total_MM(asf_set%csf_set%csf(i)%totalJ)
         call expand_csf_into_determinants(asf_set%csf_set,asf_new%det_set, &
                                           i,MM,det_list,weight,nod)
         call append_determinants(asf_mode,asf_set,asf_new,i,MM, &
                                           det_list,weight,nod)
         do  j = 1,nod
            deallocate( det_list(j)%occupation )
         end do
         deallocate( det_list, weight )
      end do
      !
      ! Write results to the CESD expansion file and close it
      if (cesd_use_compact_xpn_file) then
         call write_expansion_to_file_compact(asf_mode,asf_new)
      else
         call write_expansion_to_file(asf_mode,asf_new)
      end if
      !
      call date_and_time(date=cdate,time=ctime)
      month = get_month(cdate(5:6))
      close (24);   close (25)
      !
      ! Deallocate memory for the new generated basis asf_new
      ! call deallocate_asf_det_basis(asf_new)
      !
      301 format( /,"Orbital reference list:", &
                  /,"-----------------------" )
      302 format(i6,"  = Number of CSF in the 'original' CSF basis")
      303 format(i6,"  = Number of ASF/CSF in 'both' bases, in the original", &
                    " one and in the (new) determinant basis")
      304 format(i6,"  = Number of (n,kappa,m) orbitals in the orbital ", &
                    "reference list")
      305 format(6(i4,') ',i3,1x,i3,1x,i3,4x))
      306 format(i4,') ',i3,1x,i3,1x,i3 )
      !
   end subroutine cesd_perform_xpansion
   !
   !
   subroutine cesd_perform_xpansion_all(asf_set,asf_new)
   !--------------------------------------------------------------------
   ! This routine carries out a 'complete expansion' of the jj-coupled 
   ! ASF in a determinants basis. 
   ! For both, a Slater determinant as well as a full determinant basis
   ! abstract data types are defines. This definition also includes an 
   ! (ordered) one-particle orbital list where each orbital is given by
   ! its by its quantum numbers (n  kappa  m). 
   !
   ! In contrast to the procedure cesd_perform_xpansion(), this procedure
   ! generates the determinant representation for all possible M substates.
   ! Such 'complete' representations are needed, for instance, for
   ! collision calculations where the matrix elements are typically not
   ! independent of total M.
   !
   ! Calls: append_determinants(), expand_csf_into_determinants(),
   !        set_orbital_reference(), write_expansion_to_file().        
   !--------------------------------------------------------------------
      !
      implicit none
      type(asf_basis), intent(in)        :: asf_set
      type(asf_det_basis), intent(inout) :: asf_new
      !
      integer             :: i, j, MM, nobit, nod, noint, pos
      type(determinant), dimension(:), pointer :: det_list
      real(kind=dp), dimension(:), pointer     :: weight
      !
      ! Allocte initial storage for the 'new' ASF detrminant basis, 
      ! i.e. for asf_new of type(asf_det_basis)
      ! Also, take over the total-J, parity, and energies;
      !
      asf_new%det_set%nod      = 0
      asf_new%det_set%nodmax   = 4 * asf_set%csf_set%nocsf
      asf_new%det_set%norbital = 0
      asf_new%det_set%noint    = 0
      asf_new%det_set%number_of_electrons = asf_set%csf_set%number_of_electrons
      !
      print *, "Allocate memory for ",asf_new%det_set%nodmax," determinants."
      allocate( asf_new%det_set%determinant(1:asf_new%det_set%nodmax) )
      !
      asf_new%noasf = 0
      do  i = 1,asf_set%noasf
         do  MM = -asf_set%asf(i)%totalJ,asf_set%asf(i)%totalJ,2
            asf_new%noasf = asf_new%noasf + 1
         end do
      end do
      asf_new%average_energy = asf_set%average_energy
      allocate( asf_new%asf(1:asf_new%noasf) )
      !
      j = 0
      do  i = 1,asf_set%noasf
         do  MM = -asf_set%asf(i)%totalJ,asf_set%asf(i)%totalJ,2
            j = j + 1
            asf_new%asf(j)%level_No = asf_set%asf(i)%level_No
            asf_new%asf(j)%totalJ   = asf_set%asf(i)%totalJ
            asf_new%asf(j)%totalM   = MM
            asf_new%asf(j)%parity   = asf_set%asf(i)%parity
            asf_new%asf(j)%energy   = asf_set%asf(i)%energy
            allocate( asf_new%asf(j)%eigenvector(1:asf_new%det_set%nodmax) )
            asf_new%asf(j)%eigenvector(1:asf_new%det_set%nodmax) = zero
         end do
      end do
      !
      ! Define an (ordered) one-electron orbital reference list which 
      ! includes all orbitals of the 'old' ASF basis (in asf_set). 
      !
      call set_orbital_reference(asf_set,asf_new%det_set)
      !
      ! Allocate memory for the occupation of determinants in the 
      ! determinant basis; first determine the number of 'integers' 
      ! which is required for storing this occupation and later initialize
      ! all (bits of) occupation to 0.
      !
      noint = bit_size(noint)
      if (mod(asf_new%det_set%norbital,noint) == 0) then
         noint = asf_new%det_set%norbital/noint
      else
         noint = asf_new%det_set%norbital/noint + 1
      end if
      asf_new%det_set%noint = noint
      do  i = 1,asf_new%det_set%nodmax
         allocate( asf_new%det_set%determinant(i)%occupation(1:noint) )
      end do
      !
      nobit = bit_size(asf_new%det_set%determinant(1)%occupation(1))
      do  i = 1,noint;   do  pos = 0,nobit-1
         asf_new%det_set%determinant(1)%occupation(i) = &
            ibclr(asf_new%det_set%determinant(1)%occupation(i),pos)
      end do;   end do
      !
      do  i = 2,asf_new%det_set%nodmax 
         asf_new%det_set%determinant(i)%occupation(1:noint) = &
                     asf_new%det_set%determinant(1)%occupation(1:noint)
      end do
      !
      ! Now expand the CSF from the original ASF basis one by one into
      ! determinants
      !
      do  i = 1,asf_set%csf_set%nocsf
         do  MM = -asf_set%csf_set%csf(i)%totalJ, &
                   asf_set%csf_set%csf(i)%totalJ,2
            call expand_csf_into_determinants(asf_set%csf_set,asf_new%det_set,&
                                              i,MM,det_list,weight,nod) 
            call append_determinants(.true.,asf_set,asf_new,i,MM, &
                                              det_list,weight,nod)
            do  j = 1,nod
               deallocate( det_list(j)%occupation )
            end do
            deallocate( det_list, weight )
         end do
      end do
      !
   end subroutine cesd_perform_xpansion_all
   !
   !
   subroutine determine_weight_coefficients(csf_set,det_set,i,MM, &
                    det_occupation,det_weight,norb,no_determinants,ndet)
   !--------------------------------------------------------------------
   ! This procedure determines the 'weights' of the determinants in
   ! det_occupation for CSF i from csf_set. These weights are formally
   ! the Fourier coefficients of the CSF with respect to antisymmetrized 
   ! uncoupled product functions of one-electron orbitals.
   !
   ! Calls: angular_momentum_j(), expand_subshell_states(),
   !        wigner_3j_symbol().
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)         :: i, MM, norb, ndet, no_determinants
      type(csf_basis), intent(in) :: csf_set
      type(det_basis), intent(in) :: det_set
      integer, dimension(1:norb,1:ndet), intent(in)   :: det_occupation
      real(kind=dp), dimension(1:ndet), intent(inout) :: det_weight
      !
      integer :: iu, io, i1, i2, i3, j12r, j1l, j2l, k, &
                 mel, m1l, m2l, m12r, m12rm, mm_test, nueges, j2ges, kapidx
      integer, dimension(-60:60) :: mqz
      real(kind=dp)              :: cg, w, sum,  ws
      !
      do i1 = 1,no_determinants
         w = one
         !
         ! Initialize j12r and m12r which occur as left-hand quantum numbers
         ! in the Clebsch-Gordan coefficients for the coupling of open shells
         j12r  = 0;   m12r  = 0;   mm_test    = 0
         do  k = 1,det_set%norbital
            mm_test = mm_test + det_set%orbital(k)%mm * det_occupation(k,i1)
         end do
         if (rabs_use_stop   .and.   mm_test /= MM) then
            stop "determine_weight_coefficients(): program stop A."
         endif
         !
         ! Couple further open shells
         do  i2 = 1,csf_set%nwshells
            if (csf_set%csf(i)%occupation(i2) /= 0) then
               mqz = 0                              ! Set whole array to 0
               !
               j1l = j12r;   m1l = m12r
               j2l = csf_set%csf(i)%subshellJ(i2)
               m2l = 0
               do  iu = 1,det_set%norbital
                  if (det_set%orbital(iu)%n == csf_set%subshell(i2)%n  .and.&
                      det_set%orbital(iu)%kappa == &
                      csf_set%subshell(i2)%kappa) then
                     goto 2
                  end if
               end do
               print *, "determine_weight_coefficients(): program stop B."
               stop
               !
             2 io  = iu  + angular_momentum_j(csf_set%subshell(i2)%kappa)+1 - 1
               do  i3 = iu,io
                  if ( det_occupation(i3,i1) == 1) then
                     mel = det_set%orbital(i3)%mm
                     m2l = m2l + mel
                     mqz(mel) = 1
                  end if
               end do
               !
               ! Calculate coefficient for the relativistic subshell i2
               j2ges  = csf_set%csf(i)%subshellJ(i2)
               kapidx = (angular_momentum_j(csf_set%subshell(i2)%kappa)+1)/ 2
               nueges = csf_set%csf(i)%seniority(i2)
               call expand_subshell_states(mqz,kapidx,nueges,j2ges,ws)
               !
               ! Calculate coefficient for the coupling of the next open shell
               j12r = csf_set%csf(i)%subshellX(i2);   m12r = 0
               do  i3 = 1,io
                  m12r = m12r + det_set%orbital(i3)%mm * det_occupation(i3,i1)
               end do
               !
               ! Calculate the Clebsch-Gordan coefficient
               m12rm = - m12r
               cg    = wigner_3j_symbol(j1l,j2l,j12r,m1l,m2l,m12rm)
               cg    = cg * sqrt(j12r+one)
               if (mod(j1l-j2l+j12r+j12r-m12r,4) /= 0) then
                  cg = - cg
               endif
               w = w * ws * cg
            end if
         end do
         ! det_weight(i1) = det_weight(i1) * w
         det_weight(i1) = w
      end do
      !
      ! Check the completness of the decomposition, i.e.
      ! <sum over i> det_weight(i) * det_weight(i)  =?=  1
      sum = zero
      do  i1 = 1,no_determinants
         sum = sum + det_weight(i1) * det_weight(i1)
      end do
      !
      if (abs(sum - one) > 1.0e-8_dp ) then
         print *,     " i-csf, cnorm, no_determinants = ",i,sum,no_determinants
         write (24,*) " i-csf, cnorm, no_determinants = ",i,sum,no_determinants
         print *,     "determine_weight_coefficients: "// &
                      "Incomplete decomposition of the given CSF."
         write (24,*) "determine_weight_coefficients: "// &
                      "Incomplete decomposition of the given CSF."
         stop
      end if
      !
   end subroutine determine_weight_coefficients
   !
   !
   subroutine expand_csf_into_determinants(csf_set,det_set,i,MM, &
                                           det_list,weight,nod)
   !--------------------------------------------------------------------
   ! Expands the configuration state function i from csf_set [which is of
   ! type(csf_basis)] into a series of Slater determinants with the 
   ! total magnetic projection MM. It returns a list of determinants
   ! for which the occupation numbers are packed 'bitwise' into standard
   ! integers; for each determinant, also a weight of the projection
   ! is returned.
   !
   ! Calls: distribute_m_in_ncells(), pack_occupation_in_integer().
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)         :: i, MM
      integer, intent(out)        :: nod
      type(csf_basis), intent(in) :: csf_set
      type(det_basis), intent(in) :: det_set
      type(determinant), dimension(:), pointer :: det_list
      real(kind=dp),     dimension(:), pointer :: weight
      !
      logical :: append
      integer :: d, dmax, j, nelec, nelmax, idet, ndet, noint, np, nq, &
                 no_determinants, no_open_shells, shell
      integer, dimension(1:max_open_shells) :: open_index
      integer, dimension(:), pointer        :: det_occ
      integer, dimension(:,:), pointer      :: det_occupation, perm_work
      real(kind=dp)                         :: unc_weight
      real(kind=dp), dimension(:), pointer  :: det_weight
      !
      noint = det_set%noint
      !
      ! Allocate memory to create a list of all possible 'determinants'
      ! for CSF i with given number of open shells;
      ! For all open shells, also find all distribution of the electrons
      ! among the (n, kappa, m) orbitals of these subshells
      !
      allocate( permshell(1:csf_set%nwshells) )
      do  shell = 1,csf_set%nwshells
         nelmax = angular_momentum_j(csf_set%subshell(shell)%kappa) + 1
         nelec  = csf_set%csf(i)%occupation(shell)
         if (nelec == 0) then
            permshell(shell)%nperm = 0
         else if (nelec == nelmax) then
            permshell(shell)%nperm = 1
         else if (0 < nelec  .and.   nelec < nelmax) then
            np = 1;   do  j = nelmax-nelec+1,nelmax;  np = np * j;   end do
            nq = 1;   do  j = 1,nelec;                nq = nq * j;   end do
            permshell(shell)%nperm     = np / nq
            allocate( permshell(shell)%perm(nelmax,permshell(shell)%nperm), &
                      perm_work(nelmax,permshell(shell)%nperm),             &
                      permshell(shell)%weight(permshell(shell)%nperm) )
            do  j = 1,permshell(shell)%nperm  
               permshell(shell)%weight(j) = zero
            end do 
            dmax = permshell(shell)%nperm
            call distribute_m_in_ncells(nelec,nelmax,d,dmax,perm_work)
            permshell(shell)%perm = perm_work
            if (rabs_use_stop   .and.  permshell(shell)%nperm /= d) then
               print *, "d, nperm = ", &
                         d, permshell(shell)%nperm
               stop "expand_csf_into_determinants(): program stop nperm."
            end if
            deallocate( perm_work )
         else if (rabs_use_stop) then
            print *, "i, shell-kappa, nwshells, nelec, nelmax = ", i, &
                     csf_set%subshell(shell)%kappa, csf_set%nwshells, &
                     nelec, nelmax
            stop "expand_csf_into_determinants(): program stop A."
         end if
      end do
      !
      ! Determine the total number of Slater determinants which,
      ! however, might be larger than the number of necessary determinants
      ! because of the total M projection.
      ndet = 1;   no_open_shells = 0
      do  shell = 1,csf_set%nwshells
         if (permshell(shell)%nperm /= 0) then
            ndet = ndet * permshell(shell)%nperm
         end if
         if (permshell(shell)%nperm /= 0  .and.  &
             permshell(shell)%nperm /= 1) then
            no_open_shells = no_open_shells + 1
            open_index(no_open_shells) = permshell(shell)%nperm
         end if
      end do
      !
      do  j = no_open_shells+1,max_open_shells;  open_index(j) = 1;  end do
      !
      ! Construct the complete list of all possible Slater determinants; 
      ! check the magnetic quantum numbers of the individual shells as 
      ! well as MM; allocate memory for this 'intermediate' list of
      ! determinants
      !
      allocate( det_occupation(1:det_set%norbital,1:ndet), &
                det_weight(1:ndet), det_occ(1:det_set%norbital) )
      call index_combination(open_index,initialization=.true.)
      !
      no_determinants = 0
      do idet = 1,ndet
         call index_combination(open_index,index=idet)
         call select_determinant(open_index,csf_set,det_set,i,MM, &
                                 append,det_occ,unc_weight)
         if (append) then
            no_determinants = no_determinants + 1
            det_occupation(1:det_set%norbital,no_determinants) =  &
               det_occ(1:det_set%norbital)
            det_weight(no_determinants) = unc_weight
         end if
      end do
      !
      ! Determine the weight coefficients of the different CSF
      !
      call determine_weight_coefficients(csf_set,det_set,i,MM,      &
                                         det_occupation,det_weight, &
                               det_set%norbital,no_determinants,ndet)
      !
      ! Allocate memory for det_list and 'pack' the determinant occupation
      ! into this list
      if (mod(i,50) == 1) then
         maximal_no_determinants = no_determinants
      else
         maximal_no_determinants = max(maximal_no_determinants,no_determinants)
      end if 
      nod = no_determinants
      allocate( det_list(1:nod), weight(1:nod) )
      weight(1:nod) = det_weight(1:nod)
      do  idet = 1,nod
         allocate( det_list(idet)%occupation(noint) )
         det_occ(1:det_set%norbital) = det_occupation(1:det_set%norbital,idet)
         call pack_occupation_in_integer(det_list(idet),det_set%noint, &
                                         det_occ,det_set%norbital)
         det_list(idet)%totalM = MM                               
         det_list(idet)%parity = csf_set%csf(i)%parity                              
      end do 
      !
      ! Deallocate all working arrays of this routine
      deallocate( det_occupation, det_weight, det_occ )
      do  shell = 1,csf_set%nwshells
         if (permshell(shell)%nperm > 1) then
            deallocate( permshell(shell)%perm, permshell(shell)%weight)
         end if
      end do
      deallocate( permshell )
      !
   end subroutine expand_csf_into_determinants
   !
   !
   subroutine expand_parent_states(nel,kappa,mkappa,nstl,nstr,nqzl,nqzr,wseq)
   !--------------------------------------------------------------------
   ! This routine calculates expansion coefficients for the decomposition 
   ! of the coupled antisymmetric parent states (i.e. for equivalent 
   ! electrons) into uncoupled product wave functions. Here, a single
   ! step of the CSF expansion is carried out; the full decomposition is 
   ! then equivalent to the successive application of an appropriate
   ! annihilation operator a_(n, kappa, m).
   !
   ! In order to transform the coupled states into uncoupled wave functions, 
   ! a summation over a complete set of quantum states in the k-particle
   ! Fock-space must be carried out. In a subsequent step, each parent
   ! state will become a 'new daugther state' for the following step;
   ! therefore, after a call to this routine the intermediate states have 
   ! to be changed from nqzl ---> nqzr.
   !
   ! The parameters are the following:
   !
   !   nel   :  Number of electrons of the right-hand side daugther states. 
   !   kappa :  kappa-value of the electron shell.
   !   mkappa:  Magnetic quantum number of the annihilation operator
   !            a_(n, kappa, m) which corresponds to the decoupled electron.
   !   nstr/nstl: Nnumber of the nel-electron dauther states respectively
   !              the (nel-1) parent states.
   !   nqzl/nqzr: Quantum numbers of the above states, i.e. 
   !              (k = 1..3: nue, Jsub, Msub).
   !   wseq:      Returns the weight factors for these decomposition.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)    :: nel, kappa, mkappa, nstl
      integer, intent(inout) :: nstr
      integer, dimension(90,3), intent(in)      :: nqzl
      integer, dimension(90,3), intent(inout)   :: nqzr
      real(kind=dp), dimension(90), intent(out) :: wseq
      !
      integer       :: jk, nelmax, nuel, jl, ml, nuer, jr, mr, &
                       i1, i2, iwd, iwp, mrm
      real(kind=dp) :: cg, wcfp
      real(kind=dp), dimension(90) :: wsl
      !
      jk     = angular_momentum_j(kappa)
      nelmax = jk  +  1
      !
      if (rabs_use_stop   .and.   nelmax <= 2) then
         stop "expand_parent_states(): program stop A."
      end if
      !
      do  i1 = 1,nstl
         nuel    = nqzl(i1,1);   jl = nqzl(i1,2);   ml = nqzl(i1,3)
         wsl(i1) = zero
         do  i2 = 1,nstr
            if (abs(wseq(i2)) < eps10) cycle
            nuer = nqzr(i2,1);   jr = nqzr(i2,2);   mr = nqzr(i2,3)
            !
            ! Calculate the weights for the decoupling of the electron
            if (ml+mkappa-mr /= 0 ) cycle
            mrm = - mr
            cg  = wigner_3j_symbol(jl,jk,jr,ml,mkappa,mrm)
            cg  = cg * sqrt(jr+one)
            if (mod(jl-jk+jr+jr-mr,4) /= 0) then
               cg = - cg
            end if
            call cfp_coefficient(nelmax,nel,jr,nuer,iwd,jl,nuel,iwp,wcfp)
            wsl(i1) = wsl(i1)  +  (cg * wcfp * wseq(i2) )
         end do
      end do
      !
      ! Rearrange the intermediate states from the left- to the right-hand
      ! side for the decoupling of the next orbital;
      ! the parent states now become the 'new' daugther states.
      nstr = nstl
      do  i1 = 1,nstr
         nqzr(i1,1) = nqzl(i1,1)
         nqzr(i1,2) = nqzl(i1,2)
         nqzr(i1,3) = nqzl(i1,3)
         wseq(i1)   = wsl(i1)
      end do
      !
   end subroutine expand_parent_states
   !
   !
   subroutine expand_subshell_states(mqz,kappa,nueges,j2ges,ws)
   !--------------------------------------------------------------------
   ! Calculates the Fourier coefficients in the projection of 
   ! antisymmetric coupled subshell states onto uncoupled (product-)
   ! functions (subshell determinants). There exist a unique projection
   ! for all subshell states for j <= 7/2 without the need to specify
   ! further quantum numbers. See the long write-up for details.
   !
   ! The parameters are the following:
   !
   !   mqz(i) :  Occupation scheme in a given subshell specified by their
   !             as specified by  magnetic quantum numbers. 
   !   kappa  :  kappa-value of the considered shell. 
   !   nueges :  Senority-quantum number of the electron shell state.
   !   j2ges  :  Total angular momentum in this shell.
   ! 
   !   ws     :  (Fourier-) weights for the uncoupled product functions which 
   !             are contained in the subshell state of the equivalent 
   !             electrons.
   !
   !   grant_ipoidw :  Quantum numbers and locations of all allowed.
   !   grant_npoidw    coupled states of equivalent electrons.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)        :: kappa, nueges, j2ges
      real(kind=dp), intent(out) :: ws
      integer, dimension(-60:60), intent(in) :: mqz
      !
      integer       :: i, ifail, ikap, im, i2, i2max, jim, j2, jt2, jt3, &
                       kk, l1, mkappa, mt3, mt3m, m2frst, m2ges, nel, nelidx, &
                       nelmha, nelmax, nelsto, &
                       nstl, nstr, nanz, nstart
      integer, dimension(90,3)     :: nqzl, nqzr
      real(kind=dp) :: cgs, fakn, wss
      real(kind=dp), dimension(90) :: wseq
      !
      ws = one
      !
      ! Determine the number of electrons and m2ges
      nel     = 0;   nelmax  = 2 * abs(kappa);   nelmha  = abs(kappa)
      ikap    = abs(kappa);   j2 = nelmax  -  1
      m2ges   = 0
      !
      do  i = -j2,j2,2
         if (mqz(i) == 1) then
            nel = nel + 1;   m2ges = m2ges + i
            if (nel == 1) then
               m2frst = i
            end if
         end if
      end do
      !
      nelsto = nel
      !
      if (nel /= 0   .and.   nel <= 8) then
         if (rabs_use_naglib)  then
            fakn = nag_s14aaf(nel+one,ifail)
         else
            fakn = factorial(nel)
         end if
      else if (rabs_use_stop) then
         stop "expand_subshell_states(): program stop A."
      end if
      ! 
      if (nel == 1) goto 35
      !
      if (rabs_use_stop) then
      if (j2 > 7   .and.  nel > 2) then
         print *, "There are maximal two electrons allowed in subshells"// &
                  " with j > 7/2."
         stop "expand_subshell_states(): program stop B."
      else if (j2  > 30) then
         stop "expand_subshell_states(): program stop C."
      endif
      endif
      !
      ! Use the cfp-expansion for the decomposition of the coupled subshell
      ! states
      nstr      = 1
      nqzr(1,1) = nueges;   nqzr(1,2) = j2ges;   nqzr(1,3) = m2ges
      wseq(1)   = one
      !
      do  i = j2,-j2,-2
         if (mqz(i) == 1) then
            mkappa = i
            !
            ! Select a complete set of intermediate states using the quantum
            ! numbers of the left-hand states
            if (nel == 2) then
               wss = zero
               do  l1 = 1,nstr
                  jt2  = j2;           jt3  = nqzr(l1,2)
                  mt3  = nqzr(l1,3);   mt3m = - mt3
                  if (m2frst+mkappa+mt3m  /= 0 ) cycle
                  !x print *, "before A"
                  cgs = wigner_3j_symbol(j2,jt2,jt3,m2frst,mkappa,mt3m)
                  !x print *, "after  A"
                  cgs = cgs * sqrt(jt3+one)
                  if (mod(jt3+jt3-mt3,4) /=  0) then
                     cgs = - cgs
                  end if
                  wss = wss +  cgs * wseq(l1)
               end do
               !
               ! Normalize with sqrt(nel-factorial)
               ws = wss * sqrt(fakn)
               goto 35
            else
               if (nel-1 > nelmha) then
                  nelidx = nelmax - nel + 1
               else
                  nelidx = nel - 1
               end if
               !
               kk = 0
               nstart = grant_ipoidw(ikap,nelidx,1)
               nanz   = grant_ipoidw(ikap,nelidx,2)
               i2max  = nstart + nanz - 1
               do  i2 = nstart,i2max
                  jim = grant_npoidw(i2,2)
                  !x print *, "nstart,i2max,i2,jim = ", nstart,i2max,i2,jim 
                  do  im = -jim,jim,2
                     kk = kk + 1
                     nqzl(kk,1) = grant_npoidw(i2,1)
                     nqzl(kk,2) = jim
                     nqzl(kk,3) = im
                  end do
               end do
               nstl = kk
            end if
            call expand_parent_states(nel,kappa,mkappa,nstl,nstr,nqzl,nqzr,wseq)
            nel = nel - 1
         endif
      end do
      !
      print *, "expand_subshell_states(): program stop D."
      stop
      !
   35 return
      !
   end subroutine expand_subshell_states
   !
   function grant_ipoidw(i,k,l)                            result(index)
   !--------------------------------------------------------------------
   ! ipoidw(jel+1/2,nel=1...4,k) supplies a pointer for grant_npoidw 
   !   k = 1 :  first index m in npoidw(m,..) 
   !   k = 2 :  number of the allowed electron states
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in) :: i, k, l
      integer :: index
      integer, dimension(4,4,2) :: ipoidw
      !
      ! I = 1
      ipoidw(1,1,1) = 1;   ipoidw(1,1,2) = 1
      ipoidw(1,2,1) = 0;   ipoidw(1,2,2) = 0
      ipoidw(1,3,1) = 0;   ipoidw(1,3,2) = 0
      ipoidw(1,4,1) = 0;   ipoidw(1,4,2) = 0
      !
      ! I = 2
      ipoidw(2,1,1) = 2;   ipoidw(2,1,2) = 1
      ipoidw(2,2,1) = 3;   ipoidw(2,2,2) = 2
      ipoidw(2,3,1) = 0;   ipoidw(2,3,2) = 0
      ipoidw(2,4,1) = 0;   ipoidw(2,4,2) = 0
      !
      ! I = 3
      ipoidw(3,1,1) = 5;   ipoidw(3,1,2) = 1
      ipoidw(3,2,1) = 6;   ipoidw(3,2,2) = 3
      ipoidw(3,3,1) = 9;   ipoidw(3,3,2) = 3
      ipoidw(3,4,1) = 0;   ipoidw(3,4,2) = 0
      !
      ! I = 4
      ipoidw(4,1,1) = 12;   ipoidw(4,1,2) = 1
      ipoidw(4,2,1) = 13;   ipoidw(4,2,2) = 4
      ipoidw(4,3,1) = 17;   ipoidw(4,3,2) = 6
      ipoidw(4,4,1) = 23;   ipoidw(4,4,2) = 8
      !
      index = ipoidw(i,k,l)
      !
   end function grant_ipoidw
   !
   !
   function grant_npoidw(i,k)                              result(index)
   !--------------------------------------------------------------------
   ! npoidw(i,k)   k = 1 :  senority-number
   !               k = 2 :  total angular momentum
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in) :: i, k
      integer :: index
      integer, dimension(30,2) :: npoidw
      !
      npoidw( 1,1) = 1;   npoidw( 1,2) = 1
      npoidw( 2,1) = 1;   npoidw( 2,2) = 3
      npoidw( 3,1) = 0;   npoidw( 3,2) = 0
      npoidw( 4,1) = 2;   npoidw( 4,2) = 4
      npoidw( 5,1) = 1;   npoidw( 5,2) = 5
      npoidw( 6,1) = 0;   npoidw( 6,2) = 0
      npoidw( 7,1) = 2;   npoidw( 7,2) = 4
      npoidw( 8,1) = 2;   npoidw( 8,2) = 8
      npoidw( 9,1) = 1;   npoidw( 9,2) = 5
      npoidw(10,1) = 3;   npoidw(10,2) = 3
      !
      npoidw(11,1) = 3;   npoidw(11,2) = 9
      npoidw(12,1) = 1;   npoidw(12,2) = 7
      npoidw(13,1) = 0;   npoidw(13,2) = 0
      npoidw(14,1) = 2;   npoidw(14,2) = 4
      npoidw(15,1) = 2;   npoidw(15,2) = 8
      npoidw(16,1) = 2;   npoidw(16,2) = 12
      npoidw(17,1) = 1;   npoidw(17,2) = 7
      npoidw(18,1) = 3;   npoidw(18,2) = 3
      npoidw(19,1) = 3;   npoidw(19,2) = 5
      npoidw(20,1) = 3;   npoidw(20,2) = 9
      !
      npoidw(21,1) = 3;   npoidw(21,2) = 11
      npoidw(22,1) = 3;   npoidw(22,2) = 15
      npoidw(23,1) = 0;   npoidw(23,2) = 0
      npoidw(24,1) = 2;   npoidw(24,2) = 4
      npoidw(25,1) = 2;   npoidw(25,2) = 8
      npoidw(26,1) = 2;   npoidw(26,2) = 12
      npoidw(27,1) = 4;   npoidw(27,2) = 4
      npoidw(28,1) = 4;   npoidw(28,2) = 8
      npoidw(29,1) = 4;   npoidw(29,2) = 10
      npoidw(30,1) = 4;   npoidw(30,2) = 16
      !
      index = npoidw(i,k)
      !
   end function grant_npoidw
   !
   !
   subroutine index_combination(open_index,index,initialization)
   !--------------------------------------------------------------------
   ! Calculates the index combination (k1,...,k9) for the index-th call.
   ! If  initialization = .true., this routine initializes the values  
   ! kimax = open_index(i); kimax denotes the (maximal) number of 
   ! 'distributions' of the i-th open shell. 
   ! Only one optional paramter, either index or initialization may occur
   ! at the same time.
   !--------------------------------------------------------------------
      !
      implicit none
      logical,  optional, intent(in) :: initialization
      integer,  optional, intent(in) :: index
      integer, dimension(1:max_open_shells), intent(inout) :: open_index
      !
      integer       :: i1, i2, i3, i4, i5, i6, i7, i8, i9
      integer, save :: i1u, i2u, i3u, i4u, i5u, i6u, i7u, i8u, i9u, ktest
      integer, save :: k1max,k2max,k3max,k4max,k5max,k6max,k7max,k8max,k9max
      !
      if (present(initialization)) then
         if (initialization) then
            i1u = 1;   i2u = 1;   i3u = 1;   i4u = 1;   i5u = 1
            i6u = 1;   i7u = 1;   i8u = 1;   i9u = 1
            ktest = 1
            k1max = open_index(1);   k2max = open_index(2)
            k3max = open_index(3);   k4max = open_index(4)
            k5max = open_index(5);   k6max = open_index(6)
            k7max = open_index(7);   k8max = open_index(8)
            k9max = open_index(9)
         else if (rabs_use_stop) then
            stop "index_combination(): program stop A."
         end if
      else if (present(index)) then   
         if (index > 0) then
            ktest = ktest - 1
            do  i1 = i1u,k1max;   open_index(1) = i1;   i1u = i1
            do  i2 = i2u,k2max;   open_index(2) = i2;   i2u = i2
            do  i3 = i3u,k3max;   open_index(3) = i3;   i3u = i3
            do  i4 = i4u,k4max;   open_index(4) = i4;   i4u = i4
            do  i5 = i5u,k5max;   open_index(5) = i5;   i5u = i5
            do  i6 = i6u,k6max;   open_index(6) = i6;   i6u = i6
            do  i7 = i7u,k7max;   open_index(7) = i7;   i7u = i7
            do  i8 = i8u,k8max;   open_index(8) = i8;   i8u = i8
            do  i9 = i9u,k9max;   open_index(9) = i9;   i9u = i9
               ktest = ktest + 1
               if (ktest == index) then
                  return
               end if
            end do;   i9u = 1;   end do;   i8u = 1;   end do;   i7u = 1
            end do;   i6u = 1;   end do;   i5u = 1;   end do;   i4u = 1
            end do;   i3u = 1;   end do;   i2u = 1;   end do
            !
            write (24,*) "index_combination(): program stop B."  
            stop "index_combination(): program stop B."
         else if (rabs_use_stop) then
            stop "index_combination(): program stop C."
         end if
      else if (rabs_use_stop) then
         stop "index_combination(): program stop D."
      endif
      !
   end subroutine index_combination
   !
   !
   subroutine select_determinant(open_index,csf_set,det_set,i,MM, &
                                 append,det_occ,weight)
   !--------------------------------------------------------------------
   ! Selects one determinant from the distribution of the subshell    
   ! determinants for the open shells.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, dimension(1:max_open_shells), intent(in)   :: open_index
      type(csf_basis), intent(in) :: csf_set
      type(det_basis), intent(in) :: det_set
      integer, intent(in)         :: i, MM 
      logical, intent(out)        :: append
      real(kind=dp), intent(out)  :: weight 
      integer, dimension(1:det_set%norbital), intent(out) :: det_occ
      !
      integer :: shell, index, iorb, j, k, MMsub, MM_test, nelmax, &
                 no_open_shell, kapidx, j2ges, nueges
      integer, dimension(-60:60) :: mqz
      real(kind=dp) :: ws
      !
      iorb          = 0
      no_open_shell = 0
      MM_test       = 0
      weight        = one
      append        = .false.
      !
      do  shell = 1,csf_set%nwshells
         nelmax = angular_momentum_j(csf_set%subshell(shell)%kappa) + 1
         do  j = 1,det_set%norbital
            if (det_set%orbital(j)%n     == csf_set%subshell(shell)%n  .and. &
                det_set%orbital(j)%kappa == csf_set%subshell(shell)%kappa) then
               exit
            end if
         end do
         if (j /= iorb+1) then
            print *, "j, iorb+1 = ",j, iorb+1
         end if
         !
         ! j is now the first index of shell in the orbital reference list
         if (permshell(shell)%nperm == 0) then
            do  k = 1,nelmax
               iorb  = iorb + 1
               det_occ(iorb) = 0
            end do
         else if (permshell(shell)%nperm == 1) then
            do  k = 1,nelmax
               iorb  = iorb + 1
               det_occ(iorb) = 1
            end do
         else 
            no_open_shell = no_open_shell + 1
            index = open_index(no_open_shell)
            MMsub = 0
            do  k = 1,nelmax
               iorb  = iorb + 1
               det_occ(iorb) = permshell(shell)%perm(k,index)
               MMsub = MMsub + permshell(shell)%perm(k,index) * &
                               det_set%orbital(iorb)%mm
               mqz(det_set%orbital(iorb)%mm) = permshell(shell)%perm(k,index)
            end do
            if (abs(MMsub) > csf_set%csf(i)%subshellJ(shell)) return
            MM_test = MM_test + MMsub
            if (abs(MM_test) > csf_set%csf(i)%subshellX(shell)) return
            !
            if (permshell(shell)%weight(index) == zero) then
               j2ges  = csf_set%csf(i)%subshellJ(shell)
               ! This should be kappa or |kappa|
               kapidx = (angular_momentum_j(csf_set%subshell(shell)%kappa)+1)/2
               nueges = csf_set%csf(i)%seniority(shell)
               call expand_subshell_states(mqz,kapidx,nueges,j2ges,ws)
               permshell(shell)%weight(index) = ws
            end if
            weight = weight * permshell(shell)%weight(index)
         end if
      end do
      !
      if (rabs_use_stop   .and.  iorb /= det_set%norbital) then
         stop "index_combination(): program stop A." 
      else if (MM_test == MM) then
         append = .true.
      end if
      !
   end subroutine select_determinant
   !
   !
   subroutine set_orbital_reference(asf_set,det_set)
   !--------------------------------------------------------------------
   ! This routine generates an orbital reference list for a determinant
   ! basis in det_set which includes all (n, kappa, m) orbitals from any
   ! CSF in the configuration scheme as given in asf_set. This reference
   ! list is used to define the occupation of determinants by its 
   ! occupation numbers (0 or 1) of the corresponding one-particle 
   ! orbitals.
   !--------------------------------------------------------------------
      !
      implicit none
      type(asf_basis), intent(in)    :: asf_set
      type(det_basis), intent(inout) :: det_set
      !
      integer :: iorb, j2, norb, norb_sub
      !
      ! The present version only allows electrons with j <= 7/2.
      ! Determine first the number of (n, kappa, m) orbitals
      !
      norb = 0
      do  iorb  = 1,asf_set%csf_set%nwshells
         norb_sub = angular_momentum_j(asf_set%csf_set%subshell(iorb)%kappa)+1
         norb     = norb + norb_sub
      end do
      !
      ! Allocate memory for the orbital reference list
      det_set%norbital = norb
      allocate( det_set%orbital(1:norb) )
      !
      ! Create the orbital reference list
      norb = 0
      do  iorb  = 1,asf_set%csf_set%nwshells
         j2 = angular_momentum_j(asf_set%csf_set%subshell(iorb)%kappa)
         do  i = -j2,j2,2
            norb = norb + 1
            det_set%orbital(norb)%n     = asf_set%csf_set%subshell(iorb)%n
            det_set%orbital(norb)%kappa = asf_set%csf_set%subshell(iorb)%kappa
            det_set%orbital(norb)%mm    = i
         end do
      end do
      !
      if(rabs_use_stop   .and.  norb /= det_set%norbital) then
         stop "set_orbital_reference(): program stop A."
      end if
      !
      write(24,300)
      write(24,301) det_set%norbital
      write(24,302) (iorb,det_set%orbital(iorb)%n,det_set%orbital(iorb)%kappa,&
                          det_set%orbital(iorb)%mm,iorb=1,det_set%norbital)
      !
      ! Print the orbital reference list on standard output if required
      if (cesd_print_orbital_reference) then
         print 300
         print 301, det_set%norbital
         print 302, (iorb,det_set%orbital(iorb)%n,det_set%orbital(iorb)%kappa,&
                          det_set%orbital(iorb)%mm,iorb=1,det_set%norbital)
      end if
      !
      300 format(//1x,"++++++++++ set_orbital_reference() ++++++++++")
      301 format( /1x,"The orbital reference list is defined by ",i3, &
                      " orbital functions.",                          &
                  /1x,"they are given in the format  i) n kappa m   ..." /)
      302 format(     3( i4, ") " ,i3,1x,i3,1x,i3, "/2   " ) )
      !
   end subroutine set_orbital_reference
   !
   !
   subroutine write_expansion_to_file(asf_mode,asf_new)
   !--------------------------------------------------------------------
   ! Writes the (newly calculated) determinant basis to the CESD output
   ! file on stream 25. The format of the output is the same for both
   ! modes (asf_mode = .true./.false.).
   ! Before the results are dumped this routine checkes that it is a
   ! "complete" and "orthogonal" representation of the symmetry functions. 
   !--------------------------------------------------------------------
      !
      implicit none
      logical, intent(in)                :: asf_mode
      type(asf_det_basis), intent(inout) :: asf_new
      !
      integer :: i, iblock, idet, isu, iso, j, nblock, norb
      integer, dimension(10)         :: iparity
      integer, dimension(:), pointer :: occupation
      real(kind=dp) :: sum
      !
      ! Test the "completeness" and "orthogonality" of the expansion;
      ! set negligible coefficients to zero
      do  i = 1,asf_new%noasf
         do  idet = 1,asf_new%det_set%nod
            if (abs(asf_new%asf(i)%eigenvector(idet)) < eps10) then
               asf_new%asf(i)%eigenvector(idet) = zero
            end if
         end do
      end do
      !
      do  i = 1,asf_new%noasf
         print *, 'Check orthogonality for ASF/CSF = ',i
         do  j = 1,asf_new%noasf
            !x print *, "j, asf_new%det_set%nod = ",j, asf_new%det_set%nod 
            sum = zero
            do  idet = 1,asf_new%det_set%nod
               sum = sum + asf_new%asf(i)%eigenvector(idet) * &
                           asf_new%asf(j)%eigenvector(idet)
            end do
            if (i == j)   sum = sum - one
            if (abs(sum) > ten * ten * eps10) then
               print 315, asf_mode,i,j,sum
            end if
         end do
      end do
      !
      ! Write the results to the output file on stream 25;
      ! first print the list of determinants in terms of their occupation
      ! numbers which refer to the orbital reference list
      !
      allocate( occupation(1:asf_new%det_set%norbital) )
      norb = asf_new%det_set%norbital
      !
      write (25,300)
      write (25,301) asf_new%det_set%nod
      do  idet = 1,asf_new%det_set%nod
         call unpack_occupation_from_integer( &
            asf_new%det_set%determinant(idet),occupation,norb)
         write (25,312) idet,(occupation(i),i=1,norb)
      enddo
      !
      if (asf_mode) then
         write (25,304)
      else
         write (25,303)
      end if
      write (25,299) asf_new%average_energy
      !
      ! Output the eigenvectors in 'blocks' of maximal 8 eigenvectors
      if (mod(asf_new%noasf,8) == 0) then
         nblock = asf_new%noasf / 8
      else
         nblock = asf_new%noasf / 8 + 1
      end if
      !
      do  iblock = 1,nblock
         isu = (iblock-1)*8 + 1
         iso = min( ((iblock-1)*8 + 8),asf_new%noasf )
         write (25,305) (i,i=isu,iso)
         write (25,*) "             "
         if (asf_mode) then
            write (25,306) (asf_new%asf(i)%level_No,i=isu,iso)
            write (25,307) (asf_new%asf(i)%energy-asf_new%average_energy, &
                            i=isu,iso)
         else 
            write (25,306) (0,   i=isu,iso)
            write (25,307) (zero,i=isu,iso)
         end if
         write (25,306) ((asf_new%asf(i)%totalJ+1),i=isu,iso)
         write (25,306) (total_MM(asf_new%asf(i)%totalJ),i=isu,iso)
         do  i = isu,iso
            if (asf_new%asf(i)%parity == "+") then;      iparity(i-isu+1) = 1
            else if (asf_new%asf(i)%parity == "-") then; iparity(i-isu+1) = -1
            end if
         end do
         write (25,306) (iparity(i-isu+1),i=isu,iso)
         do  idet = 1,asf_new%det_set%nod
            write (25,307) (asf_new%asf(i)%eigenvector(idet),i=isu,iso)
         end do
      enddo
      !
      ! Print a summary of the results in the cesd.sum file in a neat format;
      ! 
      write (24,309) asf_new%det_set%nod,trim(cesd_xpn_file)
      !
      if (cesd_print_complete_expansion) then
         print 308,     asf_new%det_set%norbital
         write (24,308) asf_new%det_set%norbital
         print 309,     asf_new%det_set%nod
         do  idet = 1,asf_new%det_set%nod
            call unpack_occupation_from_integer( &
               asf_new%det_set%determinant(idet),occupation,norb)
            print 302,     idet,(occupation(i),i=1,norb)
            write (24,302) idet,(occupation(i),i=1,norb)
         end do
         !
         if (asf_mode) then
            print 304;   write (24,304)
         else
            print 303;   write (24,303)
         end if
         !
         do  iblock = 1,nblock
            isu = (iblock-1)*8 + 1
            iso = min( ((iblock-1)*8 + 8),asf_new%noasf )
            print 305,     (i,i=isu,iso)
            write (24,305) (i,i=isu,iso)
            if( iblock == nblock   .and.   mod(asf_new%noasf,8) /= 0) then
               print *,     " "
               write (24,*) " "
            endif
            do  idet = 1,asf_new%det_set%nod
               print 307,     (asf_new%asf(i)%eigenvector(idet),i=isu,iso)
               write (24,307) (asf_new%asf(i)%eigenvector(idet),i=isu,iso)
            end do
         end do      
      end if
      !
      ! Print a summary of the results in the cesd.sum file in a neat format
      ! 
      write (24,309) asf_new%det_set%nod,trim(cesd_xpn_file)
      !
      deallocate( occupation )
      !
      299 format( /,2x,1pe15.8," = averaged energy")
      300 format( /,"List of determinants:",  &
                  /,"---------------------" )
      301 format(i6,"  = Number of determinants.")
      302 format(  3x,i4,")  ",3( 10i2,3x )                      &
                  /        10x,3( 10i2,3x )  /  10x,3( 10i2,3x ) &
                  /        10x,3( 10i2,3x )  /  10x,3( 10i2,3x ) &
                  /        10x,3( 10i2,3x )  /  10x,3( 10i2,3x ) &
                  /        10x,3( 10i2,3x ) )
      303 format( /,"CESD expansion of the CSF:",  &
                  /,"--------------------------" )
      304 format( /,"CESD expansion of the ASF:",  &
                  /,"--------------------------" )
      305 format( /,10(6x,i5,6x) /)
      306 format(   10(6x,i5,6x)  )
      307 format(   10(2x,1pe15.8)  )
      308 format( /1x,"++++++++++ write_expansion_to_file() ++++++++++",  &
                 //1x,"Occupation numbers in the determinants refer to ", &
                      "the orbital reference list with ",i5,              &
                      " orbital functions." /)
      309 format( /1x,"The complete expansion contains ",i6," determinants;",&
                  /1x,"it has been written to the file: ",a,"."/)
      312 format(  1x,i6,')  ',33( 10i2,1x ) ) 
      315 format( /1x,"Mode =",l2,": i, j, sum-delta(i,j) = ",i4,1x,i4,3x,e11.4)
      !
   end subroutine write_expansion_to_file
   !
   !
   subroutine write_expansion_to_file_compact(asf_mode,asf_new)
   !--------------------------------------------------------------------
   ! Writes the (newly calculated) determinant basis to the CESD output
   ! file on stream 25. It uses a new more compact format compared with
   ! then earlier versions of the CESD program; however, it remains the
   ! same output for both modes (asf_mode = .true./.false.) 
   ! Before the results are dumped this routine checkes that it is a
   ! "complete" and "orthogonal" representation of the symmetry functions. 
   !--------------------------------------------------------------------
      !
      implicit none
      logical, intent(in)                :: asf_mode
      type(asf_det_basis), intent(inout) :: asf_new
      !
      integer :: i, iblock, idet, isu, iso, j, nblock, norb, nobit
      integer, dimension(:), pointer :: occupation
      real(kind=dp) :: suma
      !
      ! Test the "completeness" and "orthogonality" of the expansion;
      ! set negligible coefficients to zero
      do  i = 1,asf_new%noasf
         do  idet = 1,asf_new%det_set%nod
            if (abs(asf_new%asf(i)%eigenvector(idet)) < eps10) then
               asf_new%asf(i)%eigenvector(idet) = zero
            end if
         end do
      end do
      !
      do  i = 1,asf_new%noasf
         print *, 'Check orthogonality for ASF/CSF = ',i
         do  j = 1,asf_new%noasf
            !x print *, "j, asf_new%det_set%nod = ",j, asf_new%det_set%nod 
            suma = zero
            do  idet = 1,asf_new%det_set%nod
               suma = suma + asf_new%asf(i)%eigenvector(idet) * &
                             asf_new%asf(j)%eigenvector(idet)
            end do
            if (i == j)   suma = suma - one
            if (abs(suma) > ten * ten * eps10) then
               print 315, asf_mode,i,j,suma
            end if
         end do
      end do
      !
      ! Write the results to the output file on stream 25;
      ! first print the list of determinants in terms of their occupation
      ! numbers which refer to the orbital reference list
      !
      allocate( occupation(1:asf_new%det_set%norbital) )
      norb = asf_new%det_set%norbital
      !
      write (25,300)
      write (25,301) asf_new%det_set%nod
      nobit = bit_size(asf_new%det_set%determinant(1)%occupation(1))
      do  idet = 1,asf_new%det_set%nod
         call unpack_occupation_from_integer(asf_new%det_set%determinant(idet),&
                                             occupation,norb)
         if (nobit <= 32) then
            write (25,312) idet, asf_new%det_set%determinant(idet)%totalM, &
                                 asf_new%det_set%determinant(idet)%parity, &
               (asf_new%det_set%determinant(idet)%occupation(i),           &
                                                  i=1,asf_new%det_set%noint)
         else
            write (25,313) idet, asf_new%det_set%determinant(idet)%totalM, &
                                 asf_new%det_set%determinant(idet)%parity, &
               (asf_new%det_set%determinant(idet)%occupation(i),           &
                                                  i=1,asf_new%det_set%noint)
         end if
      end do
      !
      if (asf_mode) then
         write (25,304)
      else
         write (25,303)
      end if
      write (25,299) asf_new%average_energy
      !
      ! Output the eigenvectors in 'blocks' of maximal 13 eigenvectors
      if (mod(asf_new%noasf,13) == 0) then
         nblock = asf_new%noasf / 13
      else
         nblock = asf_new%noasf / 13 + 1
      end if
      !x print *, "nblock = ",nblock
      !
      do  iblock = 1,nblock
         isu = (iblock-1)*13 + 1
         iso = min( ((iblock-1)*13 + 13),asf_new%noasf )
         write (25,305) (i,i=isu,iso)
         if (asf_mode) then
            write (25,306,advance="no") (asf_new%asf(i)%level_No,i=isu,iso)
            write (25,*) "    Atomic level numbers."
            write (25,317,advance="no") (asf_new%asf(i)%energy,  i=isu,iso)
            write (25,*) "        Atomic level energies."
         else 
            write (25,306,advance="no") (0,   i=isu,iso)
            write (25,*) "    No atomic level numbers for a CSF expansion."
            write (25,317,advance="no") (zero,i=isu,iso)
            write (25,*) "        No atomic level energies for a CSF expansion."
         end if
         write (25,306,advance="no") ((asf_new%asf(i)%totalJ),i=isu,iso)
         write (25,*) "    Total 2*J values"
         write (25,306,advance="no") (total_MM(asf_new%asf(i)%totalJ),i=isu,iso)
         write (25,*) "    Total 2*M values"
         write (25,310,advance="no") (asf_new%asf(i)%parity,i=isu,iso)
         write (25,*) "Total parities"
         do  idet = 1,asf_new%det_set%nod
            write (25,307) (asf_new%asf(i)%eigenvector(idet),i=isu,iso)
         end do
      enddo
      !
      if (cesd_print_complete_expansion) then
         print 308,     asf_new%det_set%norbital
         write (24,308) asf_new%det_set%norbital
         print 309,     asf_new%det_set%nod
         do  idet = 1,asf_new%det_set%nod
            call unpack_occupation_from_integer( &
               asf_new%det_set%determinant(idet),occupation,norb)
            print 302,     idet,(occupation(i),i=1,norb)
            write (24,302) idet,(occupation(i),i=1,norb)
         end do
         !
         if (asf_mode) then
            print 304;   write (24,304)
         else
            print 303;   write (24,303)
         end if
         !
         do  iblock = 1,nblock
            isu = (iblock-1)*8 + 1
            iso = min( ((iblock-1)*8 + 8),asf_new%noasf )
            print 305,     (i,i=isu,iso)
            write (24,305) (i,i=isu,iso)
            if( iblock == nblock   .and.   mod(asf_new%noasf,8) /= 0) then
               print *,     " "
               write (24,*) " "
            endif
            do  idet = 1,asf_new%det_set%nod
               print 307,     (asf_new%asf(i)%eigenvector(idet),i=isu,iso)
               write (24,307) (asf_new%asf(i)%eigenvector(idet),i=isu,iso)
            end do
         end do      
      end if
      !
      ! Print a summary of the results in the cesd.sum file in a neat format;
      ! 
      write (24,309) asf_new%det_set%nod,trim(cesd_xpn_file)
      !
      deallocate( occupation )
      !
      299 format( /,2x,e26.19," = averaged energy")
      300 format( /,"List of determinants:",  &
                  /,"---------------------" )
      301 format(i6,"  = Number of determinants.")
      302 format(  3x,i4,")  ",3( 10i2,3x )                      &
                  /        10x,3( 10i2,3x )  /  10x,3( 10i2,3x ) &
                  /        10x,3( 10i2,3x )  /  10x,3( 10i2,3x ) &
                  /        10x,3( 10i2,3x )  /  10x,3( 10i2,3x ) &
                  /        10x,3( 10i2,3x ) )
      303 format( /,"CESD expansion of the CSF:",  &
                  /,"--------------------------" )
      304 format( /,"CESD expansion of the ASF:",  &
                  /,"--------------------------" )
      305 format( /,15(6x,i5,6x) /)
      306 format(   15(6x,i5,6x)  )
      307 format(   15(2x,1pe15.8)  )
      317 format(   15(2x,e26.19)  )
      308 format( /1x,"++++++++++ write_expansion_to_file() ++++++++++",  &
                 //1x,"Occupation numbers in the determinants refer to ", &
                      "the orbital reference list with ",i5,              &
                      " orbital functions." /)
      309 format( /1x,"The complete expansion contains ",i6," determinants;",&
                  /1x,"it has been written to the file: ",a,"."/)
      310 format(   15(6x,4x,a1,6x)  )
      312 format(  1x,i6,")  ",i4,a1,2x,40(i12) ) 
      313 format(  1x,i6,")  ",i4,a1,2x,40(i24) ) 
      315 format( /1x,"Mode =",l2,": i, j, sum-delta(i,j) = ",i4,1x,i4,3x,e11.4)
      !
   end subroutine write_expansion_to_file_compact
   !
end module rabs_cesd
