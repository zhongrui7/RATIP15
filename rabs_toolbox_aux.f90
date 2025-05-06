module rabs_toolbox_aux
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains procedures which help carry out a number of
! 'utility' tasks within the RATIP environment. 
!-----------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_grasp2k
   implicit none
   !
   public ::  toolbox_convert_probability
                 ! Returns the oscillator strength, the Einstein A and B 
                 ! coefficients, and the decay width for a given transition 
                 ! amplitude.
   private :: toolbox_expand_dirac_orbital  
                 ! Carries out the expansion of an (n,kappa,m) Dirac orbital
                 ! in terms of predefined (n,l,m_l,m_s) non-relativistic
                 ! orbitals.
   public ::  toolbox_expand_rel_detbasis
                 ! Carries out the expansion a two- or three-electron 
                 ! relativistic determinant basis into a nonrelativistic 
                 ! (n,l,m_l,m_s) spin-orbital basis.
   private :: toolbox_expand_rel_determinant
                 ! Carries out the expansion a single two- or three-electron 
                 ! relativistic determinant in a nonrelativistic (n,l,m_l,m_s)
                 ! spin-orbital basis.
   private :: toolbox_expand_nonrel_basis
                 ! Determines all possible nonrelativistic determinants,
                 ! based on (n,l,m_l,m_s) spin-orbitals.
   !!x public  :: toolbox_initialize_rwf_storage
   !!x               ! Initializes the arrays of type(grasp2k_orbital) for the 
   !!x               ! storage of the radial wave functions.
   public  :: toolbox_interpolate_orbital
                 ! Interpolates a non-relativistic orbital.
   public  :: toolbox_interprete_levels
                 ! Attempts to interpret the serial level numbers from a 
                 ! string.
   public  :: toolbox_gather_transitions
                 ! Collects all transition information from one or several .trn
                 ! REOS/EINSTEIN transition data files. 
   public  :: toolbox_select_reos_amp
                 ! Returns (if available) a requested transition amplitude
		 ! from a list of REOS/EINSTEIN transition amplitudes.
   !
   ! Storage for the initial and final atomic states and wave functions
   type(asf_basis), public                :: asf_set, asf_initial, asf_final
   type(asf_det_basis), public            :: asf_new
   type(asf_basis), dimension(10), public :: asf_mix
   type(grasp2k_orbital), public :: wave, wave_initial, wave_final
   !
   ! Define storage for the overlap integrals: (kappa,pqn_f,pqn_i) 
   real(kind=dp), dimension(:,:,:), allocatable :: overlap
   !
   ! Define a few useful data types for internal use inside of this module
   type :: multipole_from_hfi
      integer          :: rank
      character(len=2) :: multipole  
      character(len=9) :: gauge 
      real(kind=dp)    :: amplitude, einstein_A  
   end type multipole_from_hfi
   !
   type, public :: transition_from_hfi
      integer          :: totalF_i, totalF_f, totalI, &
                          levelJ_i, levelJ_f, totalJ_i, totalJ_f, &
			  asf_i, asf_f
      integer          :: number_of_mlines
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: energy, initial_energy
      type(multipole_from_hfi), dimension(:), pointer :: mline
   end type transition_from_hfi
   !
   !
   type :: multipole_from_reos
      integer          :: rank
      character(len=2) :: multipole  
      character(len=9) :: gauge 
      real(kind=dp)    :: amplitude, einstein_A, branching_ratio   
   end type multipole_from_reos
   !
   type, public :: transition_from_reos
      integer          :: level_i, level_f, totalJ_i, totalJ_f
      integer          :: number_of_mlines
      character(len=1) :: parity_i, parity_f
      character(len=20):: file_name
      real(kind=dp)    :: energy, initial_energy, total_A_rate_Coulomb, &
                          total_A_rate_Babushkin
      type(multipole_from_reos), dimension(:), pointer :: mline
   end type transition_from_reos
   !
   ! Define global logical flags for the control of the UTILITY program; the
   ! default values for these flags may be overwritten interactively during 
   ! input time
   !!x logical, public :: toolbox_energy_inverse            = .false.,  &
   logical, public :: toolbox_use_formatted_rwf_file    = .false.,  &
                      toolbox_print_AB_in_hartree       = .false.
   !
   ! Energy unit for the output of all energies
   !!x real(kind=dp)    :: toolbox_energy_factor      = zero
   !!x character(len=7) :: toolbox_energy_unit
   !
   ! Define logical flags for debugging individual procedures
   logical, private :: debug_set_integrals              = .false.
   !
contains
   !
   subroutine toolbox_aux()
   !--------------------------------------------------------------------
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      print *, "**************************"
      print *, "*** Not yet implmented ***"
      print *, "**************************"
      !
   end subroutine toolbox_aux
   !
   !
   subroutine toolbox_convert_probability(totalJ_i,totalJ_f,L,energy_au, &
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
      !
   end subroutine toolbox_convert_probability
   !
   !
   subroutine toolbox_expand_dirac_orbital(n,kappa,mm,orblist,norbital,coeff)
   !--------------------------------------------------------------------
   ! Carries out the expansion of an (n,kappa,m) Dirac orbital in a  
   ! predefined basis of (n,l,m_l,m_s) non-relativistic orbitals in orblist.
   !
   !--------------------------------------------------------------------
      !
      integer, intent(in)                            :: n, kappa, mm, norbital
      type(nlmms), dimension(1:norbital), intent(in) :: orblist
      real(kind=dp), dimension(:), intent(out)       :: coeff
      !
      integer       :: counter, i, j, l, phase
      real(kind=dp) :: sum
      !
      counter = 0;   sum = zero;   coeff(1:norbital) = zero
      !
      j = angular_momentum_j(kappa);   l = angular_momentum_l(kappa)
      do  i = 1,norbital
         if (n  == orblist(i)%n         .and.  l == orblist(i)%l  .and. &
             mm == 2*orblist(i)%m-1     .and. -1 == orblist(i)%mms)  then
            counter = counter + 1
            phase   = l + l - 1 + mm + 32
            if (mod(phase,4) == 0) then
               coeff(i) =   wigner_3j_symbol(l+l,1,j,2*orblist(i)%m,-1,-mm) &
                            *sqrt(j+one)
            else if (mod(phase,4) == 2) then
               coeff(i) = - wigner_3j_symbol(l+l,1,j,2*orblist(i)%m,-1,-mm) &
                            *sqrt(j+one)
            else
               stop "toolbox_expand_dirac_orbital(): program stop A."
            end if
            sum = sum + coeff(i)*coeff(i)
         else if (n  == orblist(i)%n      .and.  l == orblist(i)%l  .and. &
                  mm == 2*orblist(i)%m+1  .and.  1 == orblist(i)%mms)  then
            counter = counter + 1
            phase   = l + l - 1 + mm + 32
            if (mod(phase,4) == 0) then
               coeff(i) =   wigner_3j_symbol(l+l,1,j,2*orblist(i)%m,1,-mm) &
                            *sqrt(j+one)
            else if (mod(phase,4) == 2) then
               coeff(i) = - wigner_3j_symbol(l+l,1,j,2*orblist(i)%m,1,-mm) &
                            *sqrt(j+one)
            else
               stop "toolbox_expand_dirac_orbital(): program stop B."
            end if
            sum = sum + coeff(i)*coeff(i)
         end if
      end do
      !
      if (abs(sum - one) > 1.0e-8_dp) then
         print *, "n, kappa, mm, counter, sum = ",n, kappa, mm, counter, sum
         stop "toolbox_expand_dirac_orbital(): program stop C."
      end if
      !
   end subroutine toolbox_expand_dirac_orbital
   !
   !
   subroutine toolbox_expand_rel_detbasis(asf_set,asf_nr)
   !--------------------------------------------------------------------
   ! Carries out the expansion a two- or three-electron relativistic 
   ! determinant basis in a given nonrelativistic (n,l,m_l,m_s) spin-orbital 
   ! basis.
   !
   !--------------------------------------------------------------------
      !
      type(asf_det_basis), intent(in)      :: asf_set
      type(asf_nrdet_basis), intent(inout) :: asf_nr
      !
      integer       :: i, i1, i2, i3, k, n1, kappa1, mm1, n2, kappa2, mm2,  &
                       n3, kappa3, mm3, n4, kappa4, mm4, level, no_electrons
      real(kind=dp) :: sum, weight
      !
      integer, dimension(1:asf_set%det_set%norbital) :: occupation
      real(kind=dp), dimension(:), allocatable       :: eigenvector
      !
      n4 = 0;   kappa4 = 0;   mm4 = 0
      !
      ! The number of electrons is the same
      no_electrons = asf_set%det_set%number_of_electrons
      !
      if (no_electrons == 2) then
         !
         ! Two-electron case
         !
         n3 = 0;   kappa3 = 0;   mm3 = 0
         ialoop: do  i = 1,asf_set%det_set%nod
            call unpack_occupation_from_integer(asf_set%det_set%determinant(i),&
                                            occupation,asf_set%det_set%norbital)
            do  i1 = 1,asf_set%det_set%norbital-1
               do  i2 = i1+1,asf_set%det_set%norbital
                  if (occupation(i1) == 1  .and.  occupation(i2) == 1) then
                     n1     = asf_set%det_set%orbital(i1)%n
                     kappa1 = asf_set%det_set%orbital(i1)%kappa
                     mm1    = asf_set%det_set%orbital(i1)%mm
                     n2     = asf_set%det_set%orbital(i2)%n
                     kappa2 = asf_set%det_set%orbital(i2)%kappa
                     mm2    = asf_set%det_set%orbital(i2)%mm
                     write(*,"(a,3i3,2x,3i3,2x,3i3)")                   &
		              "n1, kappa1, mm1, n2, kappa2, mm2 = ",    &
                               n1, kappa1, mm1, n2, kappa2, mm2
                     !
                     call toolbox_expand_nonrel_basis(no_electrons,     &
		            n1,kappa1,mm1,n2,kappa2,mm2,                &
			    n3,kappa3,mm3,n4,kappa4,mm4,asf_nr%nrdet_set)
                     cycle ialoop
                  end if
               end do
            end do
            !
            ! The two electrons must be found somewhere
            print *, "no_electrons = ",no_electrons
            stop "toolbox_expand_rel_detbasis(): program stop A."
         end do  ialoop
         ! 
         allocate (eigenvector(1:asf_nr%nrdet_set%nodmax) )
         !
         ! Now expand all ASF from the relativistic determinant basis
         do  level = 1,asf_set%noasf
            eigenvector(:) = zero
            ibloop: do  i = 1,asf_set%det_set%nod                 
            call unpack_occupation_from_integer(asf_set%det_set%determinant(i),&
                                            occupation,asf_set%det_set%norbital)
               do  i1 = 1,asf_set%det_set%norbital
                  do  i2 = 1,asf_set%det_set%norbital
                     if (i1 == i2) cycle
                     if (occupation(i1) == 1  .and.  occupation(i2) == 1) then
                        n1     = asf_set%det_set%orbital(i1)%n
                        kappa1 = asf_set%det_set%orbital(i1)%kappa
                        mm1    = asf_set%det_set%orbital(i1)%mm
                        n2     = asf_set%det_set%orbital(i2)%n
                        kappa2 = asf_set%det_set%orbital(i2)%kappa
                        mm2    = asf_set%det_set%orbital(i2)%mm
                        weight = asf_set%asf(level)%eigenvector(i)
                        call toolbox_expand_rel_determinant(                  &
			      no_electrons,                                   &
                              n1,kappa1,mm1,n2,kappa2,mm2,n3,kappa3,mm3,      &
                          level,kappa4,mm4,asf_nr%nrdet_set,eigenvector,weight)
                        cycle ibloop
                     end if
                  end do
               end do
            end do  ibloop 
            !
            ! Test normalization in the new representation and store this
            ! eigenvector
            sum = zero
            do  k = 1,asf_nr%nrdet_set%nod
               sum = sum + eigenvector(k)*eigenvector(k)
            end do
            if (abs(sum - one) > 1.0e-6_dp) then
               print *, " level, sum = ",level, sum
               ! stop "toolbox_expand_rel_detbasis(): program stop B."
            end if
            !
            asf_nr%asf(level)%energy         = asf_set%asf(level)%energy
            asf_nr%asf(level)%parity         = asf_set%asf(level)%parity
            asf_nr%asf(level)%totalJ         = asf_set%asf(level)%totalJ
            asf_nr%asf(level)%totalM         = asf_set%asf(level)%totalM
            asf_nr%asf(level)%eigenvector(:) = eigenvector(:)
         end do
      else if (no_electrons == 3) then
         !
         ! Three-electron case
         !
         icloop: do  i = 1,asf_set%det_set%nod
            call unpack_occupation_from_integer(asf_set%det_set%determinant(i),&
                                            occupation,asf_set%det_set%norbital)
            do  i1 = 1,asf_set%det_set%norbital-2
               do  i2 = i1+1,asf_set%det_set%norbital-1
                  do  i3 = i2+1,asf_set%det_set%norbital
                     if (occupation(i1) == 1  .and.  occupation(i2) == 1     &
                                              .and.  occupation(i3) == 1) then
                        n1     = asf_set%det_set%orbital(i1)%n
                        kappa1 = asf_set%det_set%orbital(i1)%kappa
                        mm1    = asf_set%det_set%orbital(i1)%mm
                        n2     = asf_set%det_set%orbital(i2)%n
                        kappa2 = asf_set%det_set%orbital(i2)%kappa
                        mm2    = asf_set%det_set%orbital(i2)%mm
                        n3     = asf_set%det_set%orbital(i3)%n
                        kappa3 = asf_set%det_set%orbital(i3)%kappa
                        mm3    = asf_set%det_set%orbital(i3)%mm
                        write(*,"(a,3i3,2x,3i3,2x,3i3)")                   &
                            "n1,kappa1,mm1,n2,kappa2,mm2,n3,kappa3,mm3 = ",&
                             n1,kappa1,mm1,n2,kappa2,mm2,n3,kappa3,mm3
                        !
                        call toolbox_expand_nonrel_basis(no_electrons,     &
			       n1,kappa1,mm1,n2,kappa2,mm2,                &
			       n3,kappa3,mm3,n4,kappa4,mm4,asf_nr%nrdet_set)
                        cycle icloop
                     end if
                  end do
               end do
            end do  
            !
            ! The two electrons must be found somewhere
            print *, "no_electrons = ",no_electrons
            stop "toolbox_expand_rel_detbasis(): program stop C."
         end do  icloop
         ! 
         allocate (eigenvector(1:asf_nr%nrdet_set%nodmax) )
         !
         ! Now expand all ASF from the relativistic determinant basis
         do  level = 1,asf_set%noasf
            eigenvector(:) = zero
            idloop: do  i = 1,asf_set%det_set%nod
            call unpack_occupation_from_integer(asf_set%det_set%determinant(i),&
                                            occupation,asf_set%det_set%norbital)
               do  i1 = 1,asf_set%det_set%norbital-2
                  do  i2 = i1+1,asf_set%det_set%norbital-1
                     do  i3 = i2+1,asf_set%det_set%norbital
                        if (occupation(i1) == 1 .and.  occupation(i2) == 1  &
                                                .and.  occupation(i3) == 1) then
                           n1     = asf_set%det_set%orbital(i1)%n
                           kappa1 = asf_set%det_set%orbital(i1)%kappa
                           mm1    = asf_set%det_set%orbital(i1)%mm
                           n2     = asf_set%det_set%orbital(i2)%n
                           kappa2 = asf_set%det_set%orbital(i2)%kappa
                           mm2    = asf_set%det_set%orbital(i2)%mm
                           n3     = asf_set%det_set%orbital(i3)%n
                           kappa3 = asf_set%det_set%orbital(i3)%kappa
                           mm3    = asf_set%det_set%orbital(i3)%mm
                           weight = asf_set%asf(level)%eigenvector(i)
                           call toolbox_expand_rel_determinant(                &
			      no_electrons,                                    &
                              n1,kappa1,mm1,n2,kappa2,mm2,n3,kappa3,mm3,       &
                              n4,kappa4,mm4,asf_nr%nrdet_set,eigenvector,weight)
                           cycle idloop
                        end if
                     end do
                  end do
               end do
            end do  idloop
            !
            ! Test normalization in the new representation and store this
            ! eigenvector
            sum = zero
            do  k = 1,asf_nr%nrdet_set%nod
               sum = sum + eigenvector(k)*eigenvector(k)
            end do
            if (abs(sum - one) > 1.0e-6_dp) then
               print *, " level, sum, J, M, P = ",level, sum, &
	                  asf_set%asf(level)%totalJ,          &
			  asf_set%asf(level)%totalM,asf_set%asf(level)%parity
               !! stop "toolbox_expand_rel_detbasis(): program stop D."
            end if
            !
            asf_nr%asf(level)%energy         = asf_set%asf(level)%energy
            asf_nr%asf(level)%parity         = asf_set%asf(level)%parity
            asf_nr%asf(level)%totalJ         = asf_set%asf(level)%totalJ
            asf_nr%asf(level)%totalM         = asf_set%asf(level)%totalM
            asf_nr%asf(level)%eigenvector(:) = eigenvector(:)
         end do
      else 
         !
         ! More-electron case
         !
         print *, "no_electrons = ",no_electrons
         stop "toolbox_expand_rel_detbasis(): program stop E."
      end if
      !
      print *, "aaax"
      !
   end subroutine toolbox_expand_rel_detbasis
   !
   !
   subroutine toolbox_expand_rel_determinant(no_electrons,                &
                 n1,kappa1,mm1,n2,kappa2,mm2,n3,kappa3,mm3,n4,kappa4,mm4, &
                                              nrdet_set,eigenvector,weight)
   !--------------------------------------------------------------------
   ! Carries out the expansion a two- or three-electron relativistic 
   ! determinant in a given nonrelativistic (n,l,m_l,m_s) spin-orbital basis.
   ! The relativistic determinant has a given weight. The current expansion
   ! is 'added' to the (inout) representation in eigenvector(:) which is
   ! defined in the (n,l,m_l,m_s) spin-orbital basis nrdet_set.
   !
   ! The expansion of four-electron determinants is currently not supported
   ! but already 'defined' by the parameter list.
   !
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: no_electrons, n1, kappa1, mm1, n2, kappa2, mm2, &
                             n3, kappa3, mm3, n4, kappa4, mm4 
      type(nrdet_basis), intent(in)              :: nrdet_set
      real(kind=dp),intent(in)                   :: weight
      real(kind=dp), dimension(:), intent(inout) :: eigenvector
      !
      integer                                        :: i, i1, i2, i3, MM
      !!x real(kind=dp)                                  :: wa
      integer, dimension(1:nrdet_set%norbital)       :: occupation 
      real(kind=dp), dimension(1:nrdet_set%norbital) :: coeff1, coeff2, coeff3
      !
      i = n4;   i = kappa4;   i = mm4   ! To use these variables formally
      !
      if (no_electrons == 2) then
         !
         ! Two-electron case
         !
         MM = mm1 + mm2
         call toolbox_expand_dirac_orbital(n1,kappa1,mm1,nrdet_set%orbital,    &
                                                    nrdet_set%norbital,coeff1)
         call toolbox_expand_dirac_orbital(n2,kappa2,mm2,nrdet_set%orbital,    &
                                                    nrdet_set%norbital,coeff2)
         do  i1 = 1,nrdet_set%norbital-1
            i2loop: do  i2 = 1,nrdet_set%norbital
               if (i1 == i2) cycle i2loop
               if (abs(coeff1(i1)*coeff2(i2)) > 1.0e-8_dp) then
                  do  i = 1,nrdet_set%nod
                     if (MM /= 2*nrdet_set%determinant(i)%totalM + &
                                 nrdet_set%determinant(i)%totalMs) cycle 
                     call unpack_occupation_from_int_nr(               &
                                  nrdet_set%determinant(i),occupation, &
                                  nrdet_set%norbital)
                     if (occupation(i1) == 1  .and.  occupation(i2) == 1) then
                        eigenvector(i) = eigenvector(i) + &
                                         weight * coeff1(i1) * coeff2(i2)
                        cycle i2loop
                     end if
                  end do
                  !
                  print *, "no_electrons = ",no_electrons
                  stop "toolbox_expand_rel_det(): program stop A."
               end if
            end do  i2loop
         end do
      else if (no_electrons == 3) then
         !
         ! Three-electron case
         !
         MM = mm1 + mm2 + mm3
         call toolbox_expand_dirac_orbital(n1,kappa1,mm1,nrdet_set%orbital,    &
                                                    nrdet_set%norbital,coeff1)
         call toolbox_expand_dirac_orbital(n2,kappa2,mm2,nrdet_set%orbital,    &
                                                    nrdet_set%norbital,coeff2)
         call toolbox_expand_dirac_orbital(n3,kappa3,mm3,nrdet_set%orbital,    &
                                                    nrdet_set%norbital,coeff3)
         do  i1 = 1,nrdet_set%norbital
            do  i2 = 1,nrdet_set%norbital
               if (i1 == i2) cycle 
               i3loop: do  i3 = 1,nrdet_set%norbital
                  if (i3 == i2   .or.   i3 == i1) cycle 
                  if (abs(coeff1(i1)*coeff2(i2)*coeff3(i3)) > 1.0e-12_dp) then
                     do  i = 1,nrdet_set%nod
                        if (MM /= 2*nrdet_set%determinant(i)%totalM + &
                                    nrdet_set%determinant(i)%totalMs) cycle 
                        call unpack_occupation_from_int_nr(               &
                                     nrdet_set%determinant(i),occupation, &
                                     nrdet_set%norbital)
                        if (occupation(i1) == 1 .and. occupation(i2) == 1 .and.&
                                                      occupation(i3) == 1)  then
                           eigenvector(i) = eigenvector(i) + weight *          &
                                          coeff1(i1) * coeff2(i2) * coeff3(i3)
                           cycle i3loop
                        end if
                     end do
                     !
                     print *, "no_electrons = ",no_electrons
                     stop "toolbox_expand_rel_det(): program stop B."
                  end if
               end do  i3loop
            end do
         end do
      else 
         !
         ! More-electron case
         !
         print *, "no_electrons = ",no_electrons
         stop "toolbox_expand_rel_determinant(): program stop C."
      end if
      !
   end subroutine toolbox_expand_rel_determinant
   !
   !
   subroutine toolbox_expand_nonrel_basis(no_electrons,n1,kappa1,mm1,     &
                               n2,kappa2,mm2,n3,kappa3,mm3,n4,kappa4,mm4, &
                               nrdet_set)
   !--------------------------------------------------------------------
   ! This routine determines all possible nonrelativistic determinants,
   ! based on (n,l,m_l,m_s) spin-orbitals which contribute to an 
   ! expansion of a relativistic determinant, built from two or more
   ! (n,kappa,m) Dirac orbitals. These determinants are appended to the
   ! list of determinants in (inout) nrdet_set if they are not already 
   ! included in that basis.
   ! The input of four-electron determinants is currently not supported
   ! but already 'defined' by the parameter list.
   !
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: no_electrons, n1, kappa1, mm1, n2, kappa2, mm2, &
                             n3, kappa3, mm3, n4, kappa4, mm4 
      type(nrdet_basis), intent(inout) :: nrdet_set
      !
      integer :: i, i1, i2, i3, j
      type(nrdeterminant)                            :: nrdeterm
      integer, dimension(1:nrdet_set%norbital)       :: occupation
      real(kind=dp), dimension(1:nrdet_set%norbital) :: coeff1, coeff2, coeff3
      !
      i = n4;   i = kappa4;   i = mm4   ! To use these variables formally
      !
      allocate( nrdeterm%occupation(1:4) )
      !
      if (no_electrons == 2) then
         !
         ! Two-electron case
         !
         call toolbox_expand_dirac_orbital(n1,kappa1,mm1,nrdet_set%orbital,    &
                                                    nrdet_set%norbital,coeff1)
         call toolbox_expand_dirac_orbital(n2,kappa2,mm2,nrdet_set%orbital,    &
                                                    nrdet_set%norbital,coeff2)
         do  i1 = 1,nrdet_set%norbital
            i2loop: do  i2 = 1,nrdet_set%norbital
               if (i1 == i2) cycle i2loop
               if (abs(coeff1(i1)*coeff2(i2)) > 1.0e-8_dp) then
                  occupation(:)  = 0
                  occupation(i1) = 1;   occupation(i2) = 1
                  call pack_occupation_in_integer_nr(nrdeterm,nrdet_set%noint, &
                                                occupation,nrdet_set%norbital)
                  !
                  ! Compare the 'current' determinant with what is already keept
                  ialoop: do  i = 1,nrdet_set%nod
                     do j = 1,nrdet_set%noint
                        if (nrdeterm%occupation(j) /= &
                            nrdet_set%determinant(i)%occupation(j)) then
                           cycle ialoop
                        end if
                     end do
                     !
                     ! Determinant already available, cycle i2
                     cycle i2loop
                  end do ialoop
                  !
                  ! Not yet available, append determinant
                  nrdet_set%nod = nrdet_set%nod + 1
                  if (nrdet_set%nod > nrdet_set%nodmax) then
                     stop "toolbox_expand_nonrel_basis(): program stop A."
                  end if
                  allocate( nrdet_set%determinant(nrdet_set%nod)%occupation( &
                                                1:nrdet_set%noint) )
                  nrdet_set%determinant(nrdet_set%nod)%occupation(:) =       &
                                              nrdeterm%occupation(:)
                  nrdet_set%determinant(nrdet_set%nod)%totalM  =             &
                             nrdet_set%orbital(i1)%m + nrdet_set%orbital(i2)%m
                  nrdet_set%determinant(nrdet_set%nod)%totalMs =             &
                         nrdet_set%orbital(i1)%mms + nrdet_set%orbital(i2)%mms
               end if
            end do i2loop
         end do  
      else if (no_electrons == 3) then
         !
         ! Three-electron case
         !
         call toolbox_expand_dirac_orbital(n1,kappa1,mm1,nrdet_set%orbital,    &
                                                    nrdet_set%norbital,coeff1)
         call toolbox_expand_dirac_orbital(n2,kappa2,mm2,nrdet_set%orbital,    &
                                                    nrdet_set%norbital,coeff2)
         call toolbox_expand_dirac_orbital(n3,kappa3,mm3,nrdet_set%orbital,    &
                                                    nrdet_set%norbital,coeff3)
	 !!					    
	 print *, "norbital, nodmax, noint = ",	nrdet_set%norbital, &
	 	   nrdet_set%nodmax,nrdet_set%noint
	 !!	   		    
         do  i1 = 1,nrdet_set%norbital
            do  i2 = i1+1,nrdet_set%norbital
               if (i1 == i2) cycle 
               i3loop: do  i3 = 1,nrdet_set%norbital
                  if (i3 == i2   .or.   i3 == i1) cycle 
                  if (abs(coeff1(i1)*coeff2(i2)*coeff3(i3)) > 1.0e-12_dp) then
                     occupation(:)  = 0
                     occupation(i1) = 1; occupation(i2) = 1; occupation(i3) = 1
                     call pack_occupation_in_integer_nr(nrdeterm,            &
                                nrdet_set%noint,occupation,nrdet_set%norbital)
                     !
                     ! Compare the 'current' determinant with what is 
                     ! already keept
                     ibloop: do  i = 1,nrdet_set%nod
                        do j = 1,nrdet_set%noint
                           if (nrdeterm%occupation(j) /= &
                               nrdet_set%determinant(i)%occupation(j)) then
                              cycle ibloop
                           end if
                        end do
                        !
                        ! Determinant already available, cycle i3
                        cycle i3loop
                     end do ibloop
                     !
                     ! Not yet available, append determinant
                     nrdet_set%nod = nrdet_set%nod + 1
                     if (nrdet_set%nod > nrdet_set%nodmax) then
                        stop "toolbox_expand_nonrel_basis(): program stop B."
                     end if
                     allocate( nrdet_set%determinant(nrdet_set%nod)%occupation(&
                                                   1:nrdet_set%noint) )
                     nrdet_set%determinant(nrdet_set%nod)%occupation(:) =      &
                                                 nrdeterm%occupation(:)
                     nrdet_set%determinant(nrdet_set%nod)%totalM        =      &
                             nrdet_set%orbital(i1)%m + nrdet_set%orbital(i2)%m &
			                             + nrdet_set%orbital(i3)%m
                     nrdet_set%determinant(nrdet_set%nod)%totalMs       =      &
                         nrdet_set%orbital(i1)%mms + nrdet_set%orbital(i2)%mms &
			                           + nrdet_set%orbital(i3)%mms
                  end if
               end do i3loop
            end do 
         end do  
      else 
         !
         ! More-electron case
         !
         print *, "no_electrons = ",no_electrons
         stop "toolbox_expand_nonrel_basis(): program stop C."
      end if
      !
      deallocate( nrdeterm%occupation )
      !
   end subroutine toolbox_expand_nonrel_basis
   !
   !!x !
   !!x subroutine toolbox_initialize_rwf_storage(normal,initial,final)
   !!x !--------------------------------------------------------------------
   !!x ! Initializes the arrays of type(grasp2k_orbital) for the storage of 
   !!x ! the radial wave functions.
   !!x !--------------------------------------------------------------------
   !!x    !
   !!x    logical, intent(in) :: normal, initial, final
   !!x    integer :: i
   !!x    !
   !!x    if (normal) then
   !!x       ! Initialize storage for wave functions
   !!x       wave%number_of_rwf = asf_set%csf_set%nwshells
   !!x       allocate( wave%rwf(1:wave%number_of_rwf) )
   !!x       do  i = 1,asf_set%csf_set%nwshells
   !!x          wave%rwf(i)%orbital%n     = asf_set%csf_set%subshell(i)%n
   !!x          wave%rwf(i)%orbital%kappa = asf_set%csf_set%subshell(i)%kappa
   !!x          wave%rwf(i)%mtp    = 0
   !!x          wave%rwf(i)%energy = zero
   !!x          wave%rwf(i)%gamma  = zero
   !!x          wave%rwf(i)%pz     = zero
   !!x          wave%rwf(i)%phase  = zero
   !!x       end do
   !!x    end if
   !!x    !
   !!x    if (initial) then
   !!x       ! Initialize storage for initial-state wave functions
   !!x       wave_initial%number_of_rwf = asf_initial%csf_set%nwshells
   !!x       allocate( wave_initial%rwf(1:wave_initial%number_of_rwf) )
   !!x       do  i = 1,asf_initial%csf_set%nwshells
   !!x          wave_initial%rwf(i)%orbital%n = asf_initial%csf_set%subshell(i)%n
   !!x          wave_initial%rwf(i)%orbital%kappa = &
   !!x 	                                asf_initial%csf_set%subshell(i)%kappa
   !!x          wave_initial%rwf(i)%mtp    = 0
   !!x          wave_initial%rwf(i)%energy = zero
   !!x          wave_initial%rwf(i)%gamma  = zero
   !!x          wave_initial%rwf(i)%pz     = zero
   !!x          wave_initial%rwf(i)%phase  = zero
   !!x       end do
   !!x    end if
   !!x    !
   !!x    if (final) then
   !!x       ! Initialize storage for final-state wave functions
   !!x       wave_final%number_of_rwf = asf_final%csf_set%nwshells
   !!x       allocate( wave_final%rwf(1:wave_final%number_of_rwf) )
   !!x       do  i = 1,asf_final%csf_set%nwshells
   !!x          wave_final%rwf(i)%orbital%n     = asf_final%csf_set%subshell(i)%n
   !!x          wave_final%rwf(i)%orbital%kappa = &
   !!x                                          asf_final%csf_set%subshell(i)%kappa
   !!x          wave_final%rwf(i)%mtp      = 0
   !!x          wave_final%rwf(i)%energy   = zero
   !!x          wave_final%rwf(i)%gamma    = zero
   !!x          wave_final%rwf(i)%pz       = zero
   !!x       end do
   !!x    end if
   !!x    !
   !!x end subroutine toolbox_initialize_rwf_storage
   !!x !
   !
   subroutine toolbox_interpolate_orbital(r,f,x,fx,nr,nx,ior)
   !--------------------------------------------------------------------
   ! Lagrange-interpolation of the (input) array f(i), defined on the
   ! (input) grid  r(i), i = 1...nr, into the new (input) grid x(k), 
   ! k = 1...nx. This interpolation scheme does not assume any dependence
   ! between the grids r(i) and x(i). For values x < r(ior+1) and 
   ! x > r(n-ior-1), the routine calculates extrapolates the values.
   !
   ! ior2 = 2 * ior : Number of the (surrounding) grid points for the 
   !                  Lagrange interpolation (i.e. order of interpolation).
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in)                      :: nr, nx, ior
      real(kind=dp), dimension(:), intent(in)  :: r, f, x
      real(kind=dp), dimension(:), intent(out) :: fx
      !
      integer                        :: i1, i2, ior1, ior2, iorl, k
      real(kind=dp)                  :: xa, y                 
      real(kind=dp), dimension(1:30) :: arg, val
      !
      ior1 = ior + 1
      ior2 = 2 * ior
      iorl = nr - ior + 1
      !
      k = 1
      do  i1 = ior1,iorl
         do  i2 = 1,ior2
            arg(i2) = r(i1 - ior1 + i2)
            val(i2) = f(i1 - ior1 + i2)
         end do
         do i2 = 1,nx
            if (x(k) > r(i1)) exit 
            xa    = x(k)
            call interpolate(arg,val,xa,y,ior2)
            fx(k) = y
            if (k+1 > nx) return
            k = k + 1
         end do
      end do
      !
      ! Calculate the values for x > r(iorl) using a rough extrapolation
      do  i2 = 1,nx
         xa = x(k)
         call interpolate(arg,val,xa,y,ior2)
         fx(k) = y
         if (k+1 > nx) return
         k = k + 1
      end do
      !
      contains
      !
      subroutine interpolate(arg,val,x,y,n)
      !-----------------------------------------------------------------
      ! Interpolates the value val(x) from the values arg(i),val(i),
      ! i=1,n by using a Lagrange-interpolation formulae.   
      !-----------------------------------------------------------------
      !
      integer                     :: n
      real(kind=dp)               :: x, y
      real(kind=dp), dimension(:) :: arg, val
      !
      integer       :: j, l
      real(kind=dp) :: pl
      !
      y = zero
      do  l = 1,n
         pl = one
         do  j = 1,n
            if (l-j /= 0) then
               pl = (x - arg(j)) * pl / (arg(l) - arg(j))
            end if
         end do
         y = y + pl * val(l)
      end do
      !
      end subroutine interpolate
      !
   end subroutine toolbox_interpolate_orbital
   !
   !
   subroutine toolbox_interprete_levels(record,levels,number_of_levels,fail)
   !--------------------------------------------------------------------
   ! Attempts to interprete the serial level numbers which are given on
   ! record. These level numbers can be given in the format:
   !                1 3 4  7 - 20  48  69 - 85
   ! Any order and 'overlapping' intervals are also supported.
   ! The procedure returns with fail = .true. if the level numbers cannot
   ! be interpreted properly (fail = .false. otherwise).
   !
   ! The level numbers are returned in the vector levels(1:number_of_levels).
   ! The procedure assumes that this vector has a sufficient dimension
   ! to store all level numbers.
   !
   ! Calls: .
   !--------------------------------------------------------------------
      !
      character(len=*), intent(in)         :: record
      logical, intent(out)                 :: fail
      integer, intent(out)                 :: number_of_levels
      integer, dimension(:), intent(inout) :: levels
      !
      logical, dimension(200)              :: low_to
      character(len=500)                   :: string
      integer                              :: a, i, lower, n
      integer, dimension(200)              :: low
      integer(kind=i1b), dimension(10000)  :: run = 0
      !
      fail = .true.
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
   10 levels(:) = 0
      do  i = 1,n
         if (low_to(i)) then
            if (low(i) <= low(i+1)) then;   run(low(i):low(i+1)) = 1
            else;                           run(low(i+1):low(i)) = 1
            end if
            !   
            if (low_to(i+1)) then;   return
            else;                    cycle
            end if
         else
            run(low(i)) = 1
         end if
      end do
      !
      number_of_levels = 0
      do  i = 1,10000
         if (run(i) == 1) then
            number_of_levels         = number_of_levels + 1
            levels(number_of_levels) = i
         end if
      end do
      !
      fail = .false.
      !
   end subroutine toolbox_interprete_levels
   !
   !
   subroutine toolbox_gather_transitions(lines,no_transitions,lines_max)
   !--------------------------------------------------------------------
   ! Collects all transition information from one or several .trn 
   ! REOS transition data files.
   !
   ! Calls: file_open(). 
   !--------------------------------------------------------------------
      !
      integer, intent(in)  :: lines_max
      integer, intent(out) :: no_transitions
      type(transition_from_reos), dimension(:), intent(out) :: lines
      !
      integer            :: i, j, k, m, n, no_files, ierr, ios 
      character(len=15)  :: file_id
      character(len=512) :: record
      character(len=256), dimension(10)        :: toolbox_trn_file
      type(multipole_from_reos), dimension(50) :: mlines
      !
      ! Request a number of .trn files and read the data
    1 print *, "Enter one (or several)  .trn REOS transition data file(s):"
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
         toolbox_trn_file(no_files) = record(1:i-1)
         record(:) = record(i:)
         goto 2
      end if
      !
      ! Try to open these files and to gather all necessary information
      no_transitions = 0;   m = 0
      do  i = 1,no_files
         call file_open(25,toolbox_trn_file(i),"formatted  ","old",ierr)
         if (ierr /= 0  .or.  len_trim(toolbox_trn_file(i)) == 0) goto 1
         !
         ! Check the header of the file; if not as expected for unformatted
         ! files, check formatted form
         read (25,"(a15)",iostat=ios)  file_id
         if (ios /= 0   .or.   file_id /= "REOS Transition") then
            print *, "ios, file header = ",ios, file_id
            print *, "Not a REOS Transition Data File;"
            close (25)
            goto 1
         end if
         !
         read(25,*)
         read(25, "(i6,a)") n
         no_transitions = no_transitions + n
         read(25,*)
         do   k = 1,n
            m = m + 1
            if (m > lines_max) then
               stop "toolbox_gather_transitions_reos(): program stop A."
            end if 
            read(25,4) lines(m)%level_i,  lines(m)%level_f,        &
                       lines(m)%totalJ_i, lines(m)%parity_i,       &
                       lines(m)%totalJ_f, lines(m)%parity_f,       &
                       lines(m)%number_of_mlines, lines(m)%energy, &
                       lines(m)%initial_energy,                    &
                       (mlines(j)%multipole, mlines(j)%gauge,      &
                        mlines(j)%amplitude, j=1,lines(m)%number_of_mlines)
            allocate( lines(m)%mline(1:lines(m)%number_of_mlines) )
            do  j = 1,lines(m)%number_of_mlines
               lines(m)%mline(j)%multipole = mlines(j)%multipole
               lines(m)%mline(j)%gauge     = mlines(j)%gauge
               lines(m)%mline(j)%amplitude = mlines(j)%amplitude
            end do
            lines(m)%file_name = toolbox_trn_file(i)
         end do
         !
         close (25)
      end do
    4 format(i4,2x,i4,2x,2(i4,a1,1x),2x,i2,1pe14.7,9x,           &
             1pe14.7,4x,50(a2,1x,a9,1x,1pe14.7,3x))
      !
      if (m /= no_transitions) then
         stop "toolbox_gather_transitions_reos(): program stop B."
      end if 
      !
   end subroutine toolbox_gather_transitions
   !
   !
   function toolbox_select_reos_amp(levelJ_i,levelJ_f,mult,               &
                                    gauge,lines,no_lines)       result(amp)
   !--------------------------------------------------------------------
   ! Returns (if available) a requested transition amplitude from a list 
   ! of REOS/EINSTEIN transition amplitudes.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in)          :: levelJ_i, levelJ_f, no_lines
      character(len=2), intent(in) :: mult
      character(len=9), intent(in) :: gauge
      type(transition_from_reos), dimension(1:no_lines), intent(in) :: lines
      !
      real(kind=dp)   :: amp
      integer         :: i,m
      !
      amp = zero
      !
      do  i = 1,no_lines
         if (lines(i)%level_i == levelJ_i  .and.                       &
	     lines(i)%level_f == levelJ_f)   then
	    do  m = 1,lines(i)%number_of_mlines
	       if (lines(i)%mline(m)%multipole == mult  .and.          &
	           lines(i)%mline(m)%gauge     == gauge)  then
		  amp = lines(i)%mline(m)%amplitude
		  return
	       end if
	    end do
	 end if
      end do
      !
   end function toolbox_select_reos_amp
   !
   !
   function toolbox_select_reos_energy(levelJ_i,levelJ_f,lines,no_lines)  &
                                                                 result(en)
   !--------------------------------------------------------------------
   ! Returns (if available) a requested transition energy from a list 
   ! of REOS/EINSTEIN transition amplitudes.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in)          :: levelJ_i, levelJ_f, no_lines
      type(transition_from_reos), dimension(1:no_lines), intent(in) :: lines
      !
      real(kind=dp)   :: en
      integer         :: i,m
      !
      en = zero
      !
      do  i = 1,no_lines
         if (lines(i)%level_i == levelJ_i  .and.                       &
	     lines(i)%level_f == levelJ_f)   then
	    en = lines(i)%energy
            return
	 end if
      end do
      !
   end function toolbox_select_reos_energy
   !
end module rabs_toolbox_aux
