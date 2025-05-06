module rabs_nonorthonormal
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains the procedures for calculating 'angular coefficients'
! for two set of (initial- and final-state) CSF which are not quite 
! orthogonal to each other. To include the 'overlap' of the one-electron
! orbitals of the same symmetry, therefore, a decomposition of the 
! (symmetry-adapted) CSF into Slater determinants is made and the expressions
! of Loewdin (1955) are used to evaluate the matrix elements.
!
! The module rabs_nonorthonormal provides procedures to calculate the
! 'generalized angular coefficients' for various zero-, one-, and 
! two-particle operators, in analogy to the ANCO components. Apart from 
! the computation of the coefficients for two full set of (initial- and
! final-state) CSF, it provides subprocedures to evaluate the angular
! coefficients for just a pair of CSF as required frequently in the
! computation of atomic properties.
!-----------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   use rabs_cesd
   use rabs_csl
   use rabs_determinant
   use rabs_dirac_orbital
   use rabs_file_handling
   use rabs_grasp2k
   use rabs_input_dialog
   use rabs_nucleus
   use rabs_print
   implicit none
   !
   public  :: nonorth_calculate_anco_pair
                 ! Calculates the angular coefficients for a given operator
		 ! and determinant expansion of the CSF.
   private :: nonorth_calculate_me_0p
                 ! Calculates the angular coefficients for a scalar
		 ! 0-particle operator and a given pair of determinants.
   private :: nonorth_calculate_me_1p
                 ! Calculates the angular coefficients for a 1-particle
		 ! operator and a given pair of determinants.
   private :: nonorth_calculate_me_2p
                 ! Calculates the angular coefficients for a (scalar) 2-particle
		 ! operator and a given pair of determinants 
   public  :: nonorth_calculate_csf_pair
                 ! Calculates the angular coefficients for a given pair r,s
		 ! of CSF.
   public  :: nonorth_calculate_csf_pair_fi
                 ! Calculates the angular coefficients for a given pair r,s
		 ! of CSF which are taken from asf_final%csf_set%csf(r) and
		 ! asf_initial%csf_set%csf(s), respectively.
   public  :: nonorth_calculate_csf_scheme_fi
                 ! Calculates the angular coefficients for two set of CSF 
		 ! as given by asf_final and asf_initial, respectively.
   public  :: nonorth_collect_input
                 ! Collects and proceeds all input for the calculation of
		 ! 'non-orthogonal' angular coefficients.
   public  :: nonorth_finalize
                 ! Deallocates all arrays which where requested at the time
		 ! of initalization.
   public  :: nonorth_initialize
                 ! Initializes all required arrays if the computations start
		 ! with a given arrays of 'overlaps'.
   public  :: nonorth_initialize_ci
                 ! Initializes all required arrays if the computations start
		 ! with a given arrays of 'overlaps' and with one electron in
		 ! the continuum.
   public  :: nonorth_initialize_fi
                 ! Initializes all required arrays if the computations start
		 ! for two given set of initial- and final-state CSF and 
		 ! radial wave functions
   public  :: nonorth_open_vnu
		 ! Opens a .vnu  Angular Coefficient file for the 
		 ! NON-ORTHOGONAL program on stream 29.
   !
   ! Define some global data of the NON-ORTHONORMAL program; most of these data 
   ! are read in during the interactive control and may overwrite existing
   ! default values
   !
   type, public :: nonorth_determinant
      integer   :: J_r, MM_r, r, nod_r, J_s, MM_s, s, nod_s
      integer   :: number_of_orbitals, number_of_electrons
      type(nkappam), dimension(:), pointer     :: orbital
      real(kind=dp), dimension(:), pointer     :: weight_r, weight_s
      type(determinant), dimension(:), pointer :: determinant_r, determinant_s
   end type nonorth_determinant
   !
   type(nonorth_determinant), private :: pair
   type(det_basis), private           :: det_orbitals
   !
   !
   ! Define data structure for an n-particle operator to keep together
   ! its tensorial properties
   type, public :: nonorth_operator
      integer :: particle, rank, parity
   end type nonorth_operator
   !
   type(nonorth_operator) :: Aoperator
   !
   !
   type, public :: nonorth_matrix_element
      integer                              :: norb, number_of_blocks, MM_r, MM_s
      integer, dimension(:,:), pointer     :: symmetry_index
      integer, dimension(:), pointer       :: block_r, block_s
      type(nkappam), dimension(:), pointer :: orbital_r, orbital_s
      integer                              :: rank, parity, MM_op
      real(kind=dp)                        :: value
   end type nonorth_matrix_element
   !
   !
   ! Define an internal structure to stores the 'pure' zero-, one- and 
   ! two-particle coefficients
   type :: nonorth_Z_coeff
      real(kind=dp)    :: Z
   end type nonorth_Z_coeff
   !
   type :: nonorth_T_coeff
      integer          :: nu    ! Should always be zero for scalar interaction.
      type(nkappa)     :: a, b
      real(kind=dp)    :: T
   end type nonorth_T_coeff
   !
   type :: nonorth_V_coeff
      integer          :: nu
      type(nkappa)     :: a, b, c, d
      real(kind=dp)    :: V
   end type nonorth_V_coeff
   !
   type :: nonorth_csf_pair
      integer          :: r, s
      integer          :: no_T_coeff, no_V_coeff
      type(nonorth_T_coeff), dimension(:), pointer :: T_coeff
      type(nonorth_V_coeff), dimension(:), pointer :: V_coeff
   end type nonorth_csf_pair
   !
   type(nonorth_Z_coeff)                 	   :: nonorth_Z_list
   type(nonorth_T_coeff), dimension(1000) 	   :: nonorth_T_list
   type(nonorth_V_coeff), dimension(1000)	   :: nonorth_V_list
   !
   !
   real(kind=dp), dimension(:,:,:), pointer :: nonorth_overlap
   !
   !
   ! Storage for the initial and final atomic states and wave functions
   type(asf_basis)       :: asf_initial, asf_final
   type(grasp2k_orbital) :: wave_initial, wave_final
   !
   ! Define global logical flags for the control of the NON-ORTHONORMAL program; 
   ! the default values for these flags may be overwritten interactively during 
   ! input time
   logical, public :: nonorth_use_formatted_mix_file      = .true.,   &
                      nonorth_use_formatted_rwf_file      = .true.,   &
		      nonorth_use_formatted_vnu_file      = .true.
   !
   ! Define some variables and arrays for processing input data from 
   ! nonorth_collect_input()
   !!x real(kind=dp), dimension(100)   :: einstein_energy_selection
   !
   !
contains
   !
   !
   subroutine nonorth_calculate_anco_pair(Aop,no_coeff)
   !--------------------------------------------------------------------
   ! Calculates the angular coefficients for a given operator Aop and
   ! pair of determinants which are (already) in terms of Slater determinants.
   ! The determinant expansion of the CSF are provided by the data 
   ! structure pair. The 'generalized angular coefficients' are returned
   ! in the data structures Z_list, T_list(:), and V_list(:) in dependence
   ! of the particle type of the operator.
   !
   ! Calls:  nonorth_calculate_me_0p, nonorth_calculate_me_1p,
   !         nonorth_calculate_me_2p, Wigner_Eckardt_geometry.
   !--------------------------------------------------------------------
      !
      type(nonorth_operator), intent(in)  :: Aop
      integer, intent(out)                :: no_coeff
      !
      integer                             :: i, j, kappa, mm, number_of_blocks,&
                                             MM_r, MM_s, nod_r, nod_s, norb_r, &
					     norb_s, r, s
      real(kind=dp)                       :: weight
      type(nonorth_matrix_element)        :: me
      integer, dimension(-10:10,-23:23)   :: index_r, index_s
      !
      integer :: no_T_me, no_T_coeff, no_V_me, no_V_coeff
      type(nonorth_Z_coeff)                  :: Z_list, Z_me
      type(nonorth_T_coeff), dimension(1000) :: T_list, T_me
      type(nonorth_V_coeff), dimension(1000) :: V_list, V_me
      !
      !
      no_T_coeff = 0;   no_V_coeff = 0;   Z_list%Z = zero
      !
      !
      do  r = 1,pair%nod_r
         do  s = 1,pair%nod_s
	    me%MM_r  = pair%MM_r;   me%MM_s = pair%MM_s
	    me%norb  = pair%number_of_electrons
	    me%value = zero
	    me%rank  = Aop%rank;    me%parity  = Aop%parity
	    me%MM_op = pair%MM_r - pair%MM_s
	    !
	    allocate( me%orbital_r(1:me%norb), me%orbital_s(1:me%norb) )
	    !
	    index_r(:,:) = zero;  index_s(:,:) = zero;  norb_r = 0;  norb_s = 0
	    do  i = 1,pair%number_of_orbitals
	       if (pair%determinant_r(r)%occupation(i) == 1) then
	          index_r(pair%orbital(i)%kappa,pair%orbital(i)%mm) =         &
		          index_r(pair%orbital(i)%kappa,pair%orbital(i)%mm) + 1
		  norb_r = norb_r + 1
		  me%orbital_r(norb_r) = pair%orbital(i)
	       end if
	       if (pair%determinant_s(s)%occupation(i) == 1) then
	          index_s(pair%orbital(i)%kappa,pair%orbital(i)%mm) =         &
		          index_s(pair%orbital(i)%kappa,pair%orbital(i)%mm) + 1
		  norb_s = norb_s + 1
		  me%orbital_s(norb_s) = pair%orbital(i)
	       end if
	    end do
	    !
	    if (norb_r /= me%norb  .or.  norb_s /= me%norb) then
	       print *, "me%norb, norb_r, norb_s = ",me%norb, norb_r, norb_s
               stop "nonorth_calculate_anco_pair(): program stop A."
	    end if
	    !
	    number_of_blocks = 0
	    do  kappa = -10,10
	       do  mm = -23,23,2
	          if (index_r(kappa,mm) > 0  .or.  index_s(kappa,mm) > 0) then
		     number_of_blocks = number_of_blocks  + 1
		  end if
	       end do
	    end do
	    !
	    me%number_of_blocks = number_of_blocks  
	    allocate( me%block_r(1:number_of_blocks),                        &
	              me%block_s(1:number_of_blocks),                        &
		      me%symmetry_index(-10:10,-23:23) )
	    me%symmetry_index(:,:) = 0
	    me%block_r(:) = 0
	    me%block_s(:) = 0
            !
	    number_of_blocks    = 0
	    do  kappa = -10,10
	       do  mm = -23,23,2
	          if (index_r(kappa,mm) > 0  .or.  index_s(kappa,mm) > 0) then
		     number_of_blocks = number_of_blocks  + 1
		     me%symmetry_index(kappa,mm) = number_of_blocks
		     !!x print *, "aa: number_of_blocks, kappa, mm = ", &
		     !!x               number_of_blocks, kappa, mm
		     !
		     me%block_r(number_of_blocks) = index_r(kappa,mm)
		     me%block_s(number_of_blocks) = index_s(kappa,mm)
		  end if
	       end do
	    end do
	    !
	    !
	    ! =============================================================
	    ! Distinguish the computation of the matrix elements due to the
	    ! character of the operator
	    ! =============================================================
	    !
	    if (Aop%particle == 0) then
	       weight = pair%weight_r(r) * pair%weight_s(s) 
	       call nonorth_calculate_me_0p(me,Z_me)
	       !
	       Z_list%Z = Z_list%Z + weight * Z_me%Z
	       !
	    else if (Aop%particle == 1) then
	       !
	       ! One-particle operators
	       ! ----------------------
	       weight = pair%weight_r(r) * pair%weight_s(s)  /               &
	             Wigner_Eckardt_geometry(pair%J_r,pair%MM_r,             &
		                         2*me%rank,me%MM_op,pair%J_s,pair%MM_s)
	       !
	       call nonorth_calculate_me_1p(me,no_T_me,T_me)
	       !
	       do  i = 1,no_T_me
	          do  j = 1,no_T_coeff
		     if (T_me(i)%a  == T_list(j)%a  .and.    &
		         T_me(i)%b  == T_list(j)%b  .and.    & 
		         T_me(i)%nu == T_list(j)%nu)      then
			T_list(j)%T = T_list(j)%T + weight * T_me(i)%T
			goto 1
		     end if
		  end do
		  !
		  no_T_coeff = no_T_coeff + 1
		  T_list(no_T_coeff)%a  = T_me(i)%a
		  T_list(no_T_coeff)%b  = T_me(i)%b
		  T_list(no_T_coeff)%nu = T_me(i)%nu
		  T_list(no_T_coeff)%T  = weight * T_me(i)%T
		  !
		1 continue
	       end do
	       !
	    else if (Aop%particle == 2) then
	       !
	       ! Two-particle operators
	       ! ----------------------
	       if (me%rank == 0  .and.  me%parity == 1) then
	          if (me%MM_r /= me%MM_s) cycle
	          weight = pair%weight_r(r) * pair%weight_s(s)
	       else
                  stop "nonorth_calculate_anco_pair(): program stop B1."
	       end if
	       !
	       call nonorth_calculate_me_2p(me,no_V_me,V_me)
	       !
	       do  i = 1,no_V_me
	          do  j = 1,no_V_coeff
		     if (V_me(i)%a  == V_list(j)%a  .and.    &
		         V_me(i)%b  == V_list(j)%b  .and.    & 
		         V_me(i)%c  == V_list(j)%c  .and.    & 
		         V_me(i)%d  == V_list(j)%d  .and.    & 
		         V_me(i)%nu == V_list(j)%nu)      then
			V_list(j)%V = V_list(j)%V + weight * V_me(i)%V
			goto 2
		     end if
		  end do
		  !
		  no_V_coeff = no_V_coeff + 1
		  V_list(no_V_coeff)%a  = V_me(i)%a
		  V_list(no_V_coeff)%b  = V_me(i)%b
		  V_list(no_V_coeff)%c  = V_me(i)%c
		  V_list(no_V_coeff)%d  = V_me(i)%d
		  V_list(no_V_coeff)%nu = V_me(i)%nu
		  V_list(no_V_coeff)%V  = weight * V_me(i)%V
		  !
		2 continue
	       end do
	       !
	    else 
               stop "nonorth_calculate_anco_pair(): program stop B."
	    end if
	    !
	    deallocate( me%orbital_r, me%orbital_s, me%block_r, me%block_s, &
	                me%symmetry_index )
	 end do
      end do
      !
      ! Copy the collected coefficients
      !
      if (Aop%particle == 0) then
	 nonorth_Z_list%Z = Z_list%Z
      else if (Aop%particle == 1) then
         no_coeff = no_T_coeff
	 if (no_coeff > 1000) then
            stop "nonorth_calculate_anco_pair(): program stop C1."
	 end if
	 nonorth_T_list(1:no_coeff) = T_list(1:no_coeff)
      else if (Aop%particle == 2) then
         no_coeff = no_V_coeff
	 if (no_coeff > 1000) then
            stop "nonorth_calculate_anco_pair(): program stop C2."
	 end if
	 nonorth_V_list(1:no_coeff) = V_list(1:no_coeff)
      else 
         stop "nonorth_calculate_anco_pair(): program stop C."
      end if
      !
   end subroutine nonorth_calculate_anco_pair
   !
   !
   subroutine nonorth_calculate_me_0p(me,Z_me)
   !--------------------------------------------------------------------
   ! Calculates the contributions of the zero-particle scalar operator "1"
   ! for a given pair of determinants, i.e. the 'overlap' of the many-electron
   ! determinants. This overlap is simply given by the determinant of the
   ! one-particle 'overlaps'. 
   !
   ! Calls: calculate_determinant().
   !--------------------------------------------------------------------
      !
      type(nonorth_matrix_element), intent(inout) :: me
      type(nonorth_Z_coeff)                       :: Z_me
      !
      integer       :: diff_symmetry, i, kr, ks, kapr, kaps, mmr, mms, nr, ns
      real(kind=dp) :: deter
      integer, dimension(me%number_of_blocks)     :: detsym_s
      real(kind=dp), dimension(me%norb,me%norb)   :: overlap
      !
      detsym_s(:) = me%block_r(:) - me%block_s(:)
      diff_symmetry = sum( abs(detsym_s) )
      if (diff_symmetry > 0) then
         Z_me%Z = zero
	 return
      end if
      !
      do  kr = 1,me%norb
         do  ks = 1,me%norb
            kapr = me%orbital_r(kr)%kappa;   mmr  = me%orbital_r(kr)%mm
	    nr   = me%orbital_r(kr)%n
            kaps = me%orbital_s(ks)%kappa;   mms  = me%orbital_s(ks)%mm
	    ns   = me%orbital_s(ks)%n
	    !
            if (kapr == kaps   .and.   mmr == mms   .and.       &
	        nr > 0         .and.   ns > 0) then
               overlap(kr,ks) = nonorth_overlap(kapr,nr,ns)
            else
               overlap(kr,ks) = zero
            end if
         end do
      end do
      !
      call calculate_determinant(overlap,deter,me%norb)
      Z_me%Z = deter
      !
   end subroutine nonorth_calculate_me_0p
   !
   !
   subroutine nonorth_calculate_me_1p(me,no_T_me,T_me)
   !--------------------------------------------------------------------
   ! Calculates the contributions of a one-particle multipole operator
   ! for a given pair of determinants. It first calculates the matrix of 
   ! co-factors D(k,l) from which the matrix elements is obtained by 
   !
   !         Sum_{k,l}  <k|op|l> * D(k,l) * (-1)**(k+l) ,
   !
   ! where <k|op|l> is the corresponding one-particle matrix element. 
   ! Here, k denotes a final and l an initial state orbital, i.e. in 
   ! general <k|op|l>  =/=  <l|op|k>. the D(k/l) are the co-determinants 
   ! to the overlap-matrix of the two Slater determinants.
   !
   ! The  no_T_me 'generalized angular coefficients' are returned in 
   ! the list T_me(:).
   !
   ! Calls: angular_momentum_j(), cofactor_1_of_overlap_matrix(),
   !        Wigner_Eckardt_geometry().
   !--------------------------------------------------------------------
      !
      type(nonorth_matrix_element), intent(inout) :: me
      integer, intent(out)                        :: no_T_me
      type(nonorth_T_coeff), dimension(:)         :: T_me
      !
      integer       :: diff_symmetry, i, jr, js, kr, ks, kapr, kaps, lr, ls, &
                       mmr, mms, nr, ns
      real(kind=dp) :: gf, minor           
      logical       :: need_calculation 
      logical, dimension(me%norb,me%norb)       :: cofactors_1_mask
      real(kind=dp), dimension(me%norb,me%norb) :: overlap, cofactors_1
      integer, dimension(me%number_of_blocks)   :: detsym_r, detsym_s
      !
      no_T_me = 0
      !
      ! Determine the logical 'mask' which of the co-factors need to be
      ! calculated and set up the overlap matrix
      cofactors_1_mask = .true.
      need_calculation = .false.
      do  kr = 1,me%norb
         do  ks = 1,me%norb
            kapr = me%orbital_r(kr)%kappa;   mmr  = me%orbital_r(kr)%mm
	    nr   = me%orbital_r(kr)%n
            kaps = me%orbital_s(ks)%kappa;   mms  = me%orbital_s(ks)%mm
	    ns   = me%orbital_s(ks)%n
            detsym_r(:) = me%block_r(:)
            detsym_s(:) = me%block_s(:)
	    detsym_r(me%symmetry_index(kapr,mmr)) =                          &
	                             detsym_r(me%symmetry_index(kapr,mmr)) - 1
	    detsym_s(me%symmetry_index(kaps,mms)) =                          &
	                             detsym_s(me%symmetry_index(kaps,mms)) - 1
	    !! print *, "  ",detsym_r(:)," | ",detsym_s(:)
            detsym_s(:) = detsym_r(:) - detsym_s(:)
            diff_symmetry = sum( abs(detsym_s) )
	    !! print *, "kr, ks, diff_symmetry = ",kr, ks, diff_symmetry
	    !
            if (diff_symmetry > 0) then
               cofactors_1_mask(kr,ks) = .false.
            else
               need_calculation = .true.
            end if
	    !
            if (kapr == kaps   .and.   mmr == mms) then
               overlap(kr,ks) = nonorth_overlap(kapr,nr,ns)
            else
               overlap(kr,ks) = zero
            end if
         end do
      end do
      !
      if (.not.need_calculation) then
         return
         print *, "nonorth_calculate_me_1p(): program return A."
      end if
      !
      call cofactor_1_of_overlap_matrix(me%norb,overlap,cofactors_1, &
        					     cofactors_1_mask)
      !
      do  kr = 1,me%norb
	 do  ks = 1,me%norb
	    minor = cofactors_1(kr,ks)
            if (abs(minor) < eps20) then
               cycle
            else if (mod(kr+ks,2) == 1) then
               minor = - minor
            end if
            kapr = me%orbital_r(kr)%kappa;   mmr  = me%orbital_r(kr)%mm
	    jr   = angular_momentum_j(kapr); lr   = angular_momentum_l(kapr)
	    nr   = me%orbital_r(kr)%n
            kaps = me%orbital_s(ks)%kappa;   mms  = me%orbital_s(ks)%mm
	    js   = angular_momentum_j(kaps); ls   = angular_momentum_l(kaps)
	    ns   = me%orbital_s(ks)%n
            !
            ! Calculate the 'geometrical factor' to the complete ME
            gf = Wigner_Eckardt_geometry(jr,mmr,2*me%rank,me%MM_op,js,mms)
	    gf = gf / sqrt(jr + one)
            if (abs(gf) > eps10   .and.  (-1)**(lr+ls) == me%parity) then
	       !
	       do  i = 1,no_T_me
		  if (T_me(i)%a   == nkappa(nr,kapr)  .and.  &
		      T_me(i)%b   == nkappa(ns,kaps)  .and.  &
		      T_me(i)%nu  == me%rank)         then
		     T_me(i)%T = T_me(i)%T + gf * minor
		     goto 1
		  end if
	       end do
	       !
	       no_T_me = no_T_me + 1
	       T_me(no_T_me)%a  = nkappa(nr,kapr)
	       T_me(no_T_me)%b  = nkappa(ns,kaps)
	       T_me(no_T_me)%nu = me%rank
	       T_me(no_T_me)%T  = gf * minor
	       !
	     1 continue
	    end if
         end do
      end do
      !
   end subroutine nonorth_calculate_me_1p
   !
   !
   subroutine nonorth_calculate_me_2p(me,no_V_me,V_me)
   !--------------------------------------------------------------------
   ! Calculates the contributions of a two-particle scalar operator
   ! for a given pair of determinants. It first calculates the matrix of 
   ! co-factors D(k1,k2;l1,l2) from which the matrix elements is obtained by 
   !
   !         Sum_{k1,k2,l1,l2}  <k1,k2|op|l1,l2> * D(k1,k2;l1,l2) 
   !                            * (-1)**(k1+l1 + ??) ,
   !
   ! where <k1,k2|op|l1,l2> is the corresponding two-particle matrix element. 
   ! Here, k1, k2 denote the final-state and l1, l2 the initial-state orbital, 
   ! i.e. in general <k1,k2|op|l1,l2>  =/=  <l1,l2|op|k1,k2>.
   ! The D(k1,k2;l1,l2) are the co-determinants to the overlap-matrix of the 
   ! two Slater determinants.
   !
   ! The  no_V_me 'generalized angular coefficients' are returned in 
   ! the list V_me(:).
   !
   ! Calls: angular_momentum_j(), cofactor_2_of_overlap_matrix(),
   !        Wigner_Eckardt_geometry().
   !--------------------------------------------------------------------
      !
      type(nonorth_matrix_element), intent(inout) :: me
      integer, intent(out)                        :: no_V_me
      type(nonorth_V_coeff), dimension(:)         :: V_me
      !
      integer       :: diff_symmetry, i, k, kkmax, q,              &
                       jr1, js1, kr1, ks1, kapr1, kaps1, lr1, ls1, &
                       mmr1, mms1, nr1, ns1,                       &
                       jr2, js2, kr2, ks2, kapr2, kaps2, lr2, ls2, &
                       mmr2, mms2, nr2, ns2                        
      real(kind=dp) :: gf, minor, wa           
      logical       :: need_calculation 
      logical, dimension(me%norb,me%norb,me%norb,me%norb) :: cofactors_2_mask
      real(kind=dp), dimension(me%norb,me%norb,me%norb,me%norb) :: cofactors_2
      real(kind=dp), dimension(me%norb,me%norb) :: overlap
      integer, dimension(me%number_of_blocks)   :: detsym_r, detsym_s
      !
      no_V_me = 0
      !
      ! Determine the logical 'mask' which of the co-factors need to be
      ! calculated and set up the overlap matrix
      cofactors_2_mask = .true.
      need_calculation = .false.
      do  kr1 = 1,me%norb
      do  kr2 = 1,me%norb
         if (kr1 == kr2) cycle
         do  ks1 = 1,me%norb
         do  ks2 = 1,me%norb
            if (ks1 == ks2) cycle
            kapr1 = me%orbital_r(kr1)%kappa;   mmr1  = me%orbital_r(kr1)%mm
	    nr1   = me%orbital_r(kr1)%n
            kapr2 = me%orbital_r(kr2)%kappa;   mmr2  = me%orbital_r(kr2)%mm
	    nr2   = me%orbital_r(kr2)%n
            kaps1 = me%orbital_s(ks1)%kappa;   mms1  = me%orbital_s(ks1)%mm
	    ns1   = me%orbital_s(ks1)%n
            kaps2 = me%orbital_s(ks2)%kappa;   mms2  = me%orbital_s(ks2)%mm
	    ns2   = me%orbital_s(ks2)%n
            detsym_r(:) = me%block_r(:)
            detsym_s(:) = me%block_s(:)
	    detsym_r(me%symmetry_index(kapr1,mmr1)) =                        &
	                           detsym_r(me%symmetry_index(kapr1,mmr1)) - 1
	    detsym_r(me%symmetry_index(kapr2,mmr2)) =                        &
	                           detsym_r(me%symmetry_index(kapr2,mmr2)) - 1
	    detsym_s(me%symmetry_index(kaps1,mms1)) =                        &
	                           detsym_s(me%symmetry_index(kaps1,mms1)) - 1
	    detsym_s(me%symmetry_index(kaps2,mms2)) =                        &
	                           detsym_s(me%symmetry_index(kaps2,mms2)) - 1
            detsym_s(:) = detsym_r(:) - detsym_s(:)
            diff_symmetry = sum( abs(detsym_s) )
	    !
            if (diff_symmetry > 0) then
               cofactors_2_mask(kr1,kr2,ks1,ks2) = .false.
            else
               need_calculation = .true.   
            end if
	    !
            if (kapr1 == kaps1   .and.   mmr1 == mms1   .and.       &
	        nr1   > 0        .and.   ns1  > 0) then
               overlap(kr1,ks1) = nonorth_overlap(kapr1,nr1,ns1)
            else
               overlap(kr1,ks1) = zero
            end if
            if (kapr2 == kaps2   .and.   mmr2 == mms2   .and.       &
	        nr2   > 0        .and.   ns2  > 0) then
               overlap(kr2,ks2) = nonorth_overlap(kapr2,nr2,ns2)
            else
               overlap(kr2,ks2) = zero
            end if
         end do
         end do
      end do
      end do
      !
      if (.not.need_calculation) then
         return
         print *, "nonorth_calculate_me_2p(): program return A."
      end if
      !
      call cofactor_2_of_overlap_matrix(me%norb,overlap,cofactors_2, &
                                                     cofactors_2_mask)
      !
      do  kr1 = 1,me%norb
      do  kr2 = 1,me%norb
         if (kr1 == kr2) cycle
	 do  ks1 = 1,me%norb
	 do  ks2 = 1,me%norb
            if (ks1 == ks2) cycle
	    minor = cofactors_2(kr1,kr2,ks1,ks2)
            if (abs(minor) < eps10) cycle
	    !! minor = minor / two
	    if (kr1 > kr2) cycle
	    !
            if (mod(kr1+kr2+ks1+ks2,2) == 1) then
               minor = - minor 
	    end if
	    if (kr1 > kr2) then
               !! minor = - minor 
            end if
	    if (ks1 > ks2) then
               minor = - minor 
            end if
	    !! print *, "kr1,kr2,ks1,ks2, minor = ",kr1,kr2,ks1,ks2, minor
            kapr1 = me%orbital_r(kr1)%kappa;   mmr1  = me%orbital_r(kr1)%mm
	    jr1   = angular_momentum_j(kapr1); lr1   = angular_momentum_l(kapr1)
	    nr1   = me%orbital_r(kr1)%n
            kapr2 = me%orbital_r(kr2)%kappa;   mmr2  = me%orbital_r(kr2)%mm
	    jr2   = angular_momentum_j(kapr2); lr2   = angular_momentum_l(kapr2)
	    nr2   = me%orbital_r(kr2)%n
            kaps1 = me%orbital_s(ks1)%kappa;   mms1  = me%orbital_s(ks1)%mm
	    js1   = angular_momentum_j(kaps1); ls1   = angular_momentum_l(kaps1)
	    ns1   = me%orbital_s(ks1)%n
            kaps2 = me%orbital_s(ks2)%kappa;   mms2  = me%orbital_s(ks2)%mm
	    js2   = angular_momentum_j(kaps2); ls2   = angular_momentum_l(kaps2)
	    ns2   = me%orbital_s(ks2)%n
            !
	    if (me%rank == 0  .and.  me%parity == 1) then
	      !
	      kkmax = min((jr1+js1)/2, (jr2+js2)/2)
	      do  k = 0,kkmax
	         gf = zero
	         do  q = -k,k
		    wa = wigner_3j_symbol(jr1,k+k,js1,-mmr1, q+q,mms1) &
		       * wigner_3j_symbol(jr2,k+k,js2,-mmr2,-q-q,mms2) 
		    if (mod(jr1 - mmr1 + jr2 - mmr2 +k+k- q-q+32,4) == 0) then
		       gf = gf + wa
		    else if                                                    &
		       (mod(jr1 - mmr1 + jr2 - mmr2 +k+k- q-q+32,4) == 2) then
		       gf = gf - wa
		    else
                       stop "nonorth_calculate_me_2p(): program stop A."
		    end if
		 end do
		 !
		 if (abs(gf) > eps10) then
	            do  i = 1,no_V_me
		       if (V_me(i)%a   == nkappa(nr1,kapr1)  .and.  &
		           V_me(i)%b   == nkappa(nr2,kapr2)  .and.  &
		           V_me(i)%c   == nkappa(ns1,kaps1)  .and.  &
		           V_me(i)%d   == nkappa(ns2,kaps2)  .and.  &
		           V_me(i)%nu  == k)         then
		          V_me(i)%V = V_me(i)%V + gf * minor
		          goto 1
		       else if (V_me(i)%a   == nkappa(nr2,kapr2)  .and.  &
		        	V_me(i)%b   == nkappa(nr1,kapr1)  .and.  &
		        	V_me(i)%c   == nkappa(ns2,kaps2)  .and.  &
		        	V_me(i)%d   == nkappa(ns1,kaps1)  .and.  &
		        	V_me(i)%nu  == k)	  then
		          V_me(i)%V = V_me(i)%V + gf * minor
		          goto 1
		       end if
	            end do
	            !
	            no_V_me = no_V_me + 1
	            V_me(no_V_me)%a  = nkappa(nr1,kapr1)
	            V_me(no_V_me)%b  = nkappa(nr2,kapr2)
	            V_me(no_V_me)%c  = nkappa(ns1,kaps1)
	            V_me(no_V_me)%d  = nkappa(ns2,kaps2)
	            V_me(no_V_me)%nu = k
	            V_me(no_V_me)%V  = gf * minor
	            !
	          1 continue
		 end if
	      end do
	    else
               stop "nonorth_calculate_me_2p(): program stop B."
	    end if
         end do
         end do
      end do
      end do
      !
   end subroutine nonorth_calculate_me_2p
   !
   !
   subroutine nonorth_calculate_csf_pair(csf_set,Aop,r,s,no_coeff)
   !--------------------------------------------------------------------
   ! Calculates the 'generalized angular coefficients' for a given pair 
   ! r, s of CSF where  r  refers to the final-state CSF (asf_final) and
   ! s  to the initial-state CSF (asf_initial). It carries out the
   ! decomposition of the CSF in terms of Slater determinants and 'stores'
   ! them into the data structure pair%.
   !
   ! Calls: cesd_expand_csf(), nonorth_calculate_anco_pair().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: r, s
      type(csf_basis), intent(in)        :: csf_set
      type(nonorth_operator), intent(in) :: Aop
      integer, intent(out)               :: no_coeff
      !
      integer       :: J_r, J_s, MM_r, MM_s, nod
      real(kind=dp) :: sum
      !
      integer, dimension(:,:), pointer         :: det_occupation
      real(kind=dp),	 dimension(:), pointer :: weight
      ! 
      J_r = csf_set%csf(r)%totalJ;   MM_r = J_r
      call cesd_expand_csf(csf_set,det_orbitals,r,MM_r, 		   &
        					  det_occupation,weight,nod)
      !
      pair%r = r;   pair%J_r = J_r;   pair%MM_r = MM_r;   pair%nod_r = nod
      allocate( pair%determinant_r(1:nod), pair%weight_r(1:nod) )
      do  i = 1,nod
         allocate( pair%determinant_r(i)%occupation(1:pair%number_of_orbitals) )
	 pair%determinant_r(i)%occupation(:) = det_occupation(:,i)
	 pair%weight_r(i) = weight(i)
      end do
      !
      deallocate( det_occupation, weight )
      !
      J_s = csf_set%csf(s)%totalJ;   MM_s = J_s
      call cesd_expand_csf(csf_set,det_orbitals,s,MM_s, 		   &
        					  det_occupation,weight,nod)
      !
      pair%s = s;   pair%J_s = J_s;   pair%MM_s = MM_s;   pair%nod_s = nod
      allocate( pair%determinant_s(1:nod), pair%weight_s(1:nod) )
      do  i = 1,nod
         allocate( pair%determinant_s(i)%occupation(1:pair%number_of_orbitals) )
	 pair%determinant_s(i)%occupation(:) = det_occupation(:,i)
	 pair%weight_s(i) = weight(i)
      end do
      !
      deallocate( det_occupation, weight )
      !
      ! Now calculate the 'generalized' angular coefficients for the given 
      ! operator and the expansion of the CSF in terms of Slater determinants
      call nonorth_calculate_anco_pair(Aop,no_coeff)
      !
      do  i = 1,nod
         !! deallocate( pair%determinant_r(i)%occupation )
      end do
      deallocate( pair%determinant_r, pair%weight_r )
      !
      do  i = 1,nod
         !! deallocate( pair%determinant_s(i)%occupation )
      end do
      deallocate( pair%determinant_s, pair%weight_s )
      !
   end subroutine nonorth_calculate_csf_pair
   !
   !
   subroutine nonorth_calculate_csf_pair_fi(Aop,r,s,no_coeff)
   !--------------------------------------------------------------------
   ! Calculates the 'generalized angular coefficients' for a given pair 
   ! r, s of CSF where  r  refers to the final-state CSF (asf_final) and
   ! s  to the initial-state CSF (asf_initial). It carries out the
   ! decomposition of the CSF in terms of Slater determinants and 'stores'
   ! them into the data structure pair%.
   !
   ! Calls: cesd_expand_csf(), nonorth_calculate_anco_pair().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: r, s
      type(nonorth_operator), intent(in) :: Aop
      integer, intent(out)               :: no_coeff
      !
      integer       :: J_r, J_s, MM_r, MM_s, nod
      real(kind=dp) :: sum
      !
      integer, dimension(:,:), pointer         :: det_occupation
      real(kind=dp),     dimension(:), pointer :: weight
      ! 
      J_r = asf_final%csf_set%csf(r)%totalJ;   MM_r = J_r
      call cesd_expand_csf(asf_final%csf_set,det_orbitals,r,MM_r, &
                           det_occupation,weight,nod)
      !
      pair%r = r;   pair%J_r = J_r;   pair%MM_r = MM_r;   pair%nod_r = nod
      allocate( pair%determinant_r(1:nod), pair%weight_r(1:nod) )
      do  i = 1,nod
         allocate( pair%determinant_r(i)%occupation(1:pair%number_of_orbitals) )
	 pair%determinant_r(i)%occupation(:) = det_occupation(:,i)
	 pair%weight_r(i) = weight(i)
      end do
      !
      deallocate( det_occupation, weight )
      !
      sum = zero
      do  i = 1,pair%nod_r
         sum = sum + pair%weight_r(i)*pair%weight_r(i)
      end do
      !!x print *, " "
      !!x print *, "r, MM_r, pair%nod_r, sum_r  = ",r,MM_r, pair%nod_r, sum
      !
      J_s = asf_initial%csf_set%csf(s)%totalJ;   MM_s = J_s
      call cesd_expand_csf(asf_initial%csf_set,det_orbitals,s,MM_s, &
                           det_occupation,weight,nod)
      !
      pair%s = s;   pair%J_s = J_s;   pair%MM_s = MM_s;   pair%nod_s = nod
      allocate( pair%determinant_s(1:nod), pair%weight_s(1:nod) )
      do  i = 1,nod
         allocate( pair%determinant_s(i)%occupation(1:pair%number_of_orbitals) )
	 pair%determinant_s(i)%occupation(:) = det_occupation(:,i)
	 pair%weight_s(i) = weight(i)
      end do
      !
      deallocate( det_occupation, weight )
      !
      sum = zero
      do  i = 1,pair%nod_s
         sum = sum + pair%weight_s(i)*pair%weight_s(i)
      end do
      !!x print *, "s, MM_s, pair%nod_s, sum_s  = ",s,MM_s, pair%nod_s, sum
      !!x print *, "========================"
      !
      ! Now calculate the 'generalized' angular coefficients for the given 
      ! operator and the expansion of the CSF in terms of Slater determinants
      call nonorth_calculate_anco_pair(Aop,no_coeff)
      !
      do  i = 1,nod
         deallocate( pair%determinant_r(i)%occupation )
      end do
      deallocate( pair%determinant_r, pair%weight_r )
      !
      do  i = 1,nod
         deallocate( pair%determinant_s(i)%occupation )
      end do
      deallocate( pair%determinant_s, pair%weight_s )
      !
   end subroutine nonorth_calculate_csf_pair_fi
   !
   !
   subroutine nonorth_calculate_csf_scheme_fi()
   !--------------------------------------------------------------------
   ! Calculates the 'generalized angular coefficients' for all pairs 
   ! (r,s) of CSF where  r  refers to the set of final-state CSF (asf_final) 
   ! and  s  to the set of initial-state CSF (asf_initial). All the 
   ! coefficients are printed to stream 29.
   !
   ! Calls: call nonorth_calculate_csf_pair_fi().
   !--------------------------------------------------------------------
      !
      integer          :: no_coeff, i, nu, r, s, no_Tx_coeff,no_Vx_coeff
      character(len=4) :: o_a, o_b, o_c, o_d
      !
      write(29,*) " "
      !
      do  r = 1,asf_final%csf_set%nocsf
         do  s = 1,asf_initial%csf_set%nocsf
            call nonorth_calculate_csf_pair_fi(Aoperator,r,s,no_coeff)
            !
	    ! Distinguish the particle type of the operator for the proper
	    ! printout of the 'generalized angular coefficients'.
            !
	    if (Aoperator%particle == 0) then
               write(29,"(a,i3,a,i3,a,es20.12)")			      &
                  "  pure zero-particle [",r,",",s,"] =",nonorth_Z_list%Z
	    else if (Aoperator%particle == 1) then
               do  i = 1,no_coeff
	          if (abs(nonorth_T_list(i)%T) > eps10) then
		     nu  = nonorth_T_list(i)%nu
                     o_a = orbital_name(nonorth_T_list(i)%a%n,                &
		                        nonorth_T_list(i)%a%kappa)
                     o_b = orbital_name(nonorth_T_list(i)%b%n,                &
		                        nonorth_T_list(i)%b%kappa)
                     write(29,"(a,i3,a,i3,a,i2,5a,es20.12)")                  &
                        "  pure one-particle [",r,",",s,"]  nu=",nu,          &
                        " (",o_a,",",o_b,") =",nonorth_T_list(i)%T
		  end if
	       end do
	    else if (Aoperator%particle == 2) then
	       !
	       write(29,*)
	       write(29,*) "NON-ORTH"
               do  i = 1,no_coeff
	          if (abs(nonorth_V_list(i)%V) > eps10) then
		     nu  = nonorth_V_list(i)%nu
                     o_a = orbital_name(nonorth_V_list(i)%a%n,                &
		                        nonorth_V_list(i)%a%kappa)
                     o_b = orbital_name(nonorth_V_list(i)%b%n,                &
		                        nonorth_V_list(i)%b%kappa)
                     o_c = orbital_name(nonorth_V_list(i)%c%n,                &
		                        nonorth_V_list(i)%c%kappa)
                     o_d = orbital_name(nonorth_V_list(i)%d%n,                &
		                        nonorth_V_list(i)%d%kappa)
                     write(29,"(a,i1,a,i3,a,i3,a,4a,3a,a,es20.12)")           &
                        "  pure two-particle [( ",nu,")]_[",r,",",s,"] (",    &
		        o_a,",",o_b,";",o_c,",",o_d,") =",nonorth_V_list(i)%V
		  end if
	       end do
	       !
	       !
	       write(29,*)
	       write(29,*) "ANCO"
	       !
	    else
               stop "nonorth_calculate_csf_scheme_fi(): program stop A."
	    end if
	 end do
      end do
      !
   end subroutine nonorth_calculate_csf_scheme_fi
   !
   !
   subroutine nonorth_collect_input()
   !--------------------------------------------------------------------
   ! Collect and proceeds all input for the NON-ORTHONORMAL program.
   !
   ! Calls: get_yes_stream().
   !--------------------------------------------------------------------
      !
      integer            :: blank_position, totalJ
      logical            :: yes
      character(len=200) :: string
      character(len=256) :: record
      !
      !
      ! Open, check, load data from, and close the  .iso  file
      call file_get_isodat_grasp2k()
      !
      ! Define the operator for which the angular coefficients need to be
      ! calculated
    1 print *, "Define the properties of the n-particle operator for which"// &
               " for which angular coefficients need to be calculated;"
      print *, " enter the particle number n, the rank, and parity,"//        &
               " for instance:    1 3 -1   or  2 0 1  or  ... "
      read *,  Aoperator%particle, Aoperator%rank, Aoperator%parity
      !
      if (Aoperator%particle < 0  .or.  Aoperator%particle > 2) then
         print *, "The program only supports zero-, one-, and two-particle"// &
	          " operators; redo ..."
         goto 1
      else if (Aoperator%rank < 0) then
         print *, "The rank of the operator must be rank >= 0; redo ..."
         goto 1
      else if (Aoperator%parity /= 1  .and.  Aoperator%parity /= -1) then
         print *, "The parity of the operator must be +1 or -1; redo ..."
         goto 1
      else if (Aoperator%particle == 0  .and.                     &
              (Aoperator%rank /= 0 .or.  Aoperator%parity /= 1)) then
         print *, "For zero-particle operators, only scalar interactions" //  &
	          " (rank=0,parity=1) are currently supported; redo ..."
         goto 1
      else if (Aoperator%particle == 2  .and.                     &
              (Aoperator%rank /= 0 .or.  Aoperator%parity /= 1)) then
         print *, "For two-particle operators, only scalar interactions" //   &
	          " (rank=0,parity=1) are currently supported; redo ..."
         goto 1
      end if
      !      
      !      
      ! Determine the parameters which determine the radial grid
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
      call nonorth_open_vnu(nonorth_use_formatted_vnu_file)
      !
   end subroutine nonorth_collect_input
   !
   !
   subroutine nonorth_finalize(typ)
   !--------------------------------------------------------------------
   ! Deallocates all arrays which where requested at the time of 
   ! initalization.
   ! 
   ! Calls: 
   !--------------------------------------------------------------------
      !
      character(len=2), intent(in) :: typ
      !
      if (typ == "ci") then
         deallocate ( pair%orbital, det_orbitals%orbital )
         deallocate ( nonorth_overlap )
      else
         stop "nonorth_finalize(): program stop A."
      end if
      !
   end subroutine nonorth_finalize
   !
   !
   subroutine nonorth_initialize(csf_set,wave_f,wave_i)
   !--------------------------------------------------------------------
   ! Initializes the computation of the 'generalized angular coefficients'
   ! if it is carried out for a single CSF scheme. It defines the 'reference
   ! list' of orbitals for the data structure pair% and calculates all
   ! the necessary one-electron orbitals. 
   ! 
   ! Calls: angular_momentum_j(), rk_integral_grasp2k_ab().
   !--------------------------------------------------------------------
      !
      type(csf_basis), intent(in)       :: csf_set
      type(grasp2k_orbital), intent(in) :: wave_i, wave_f
      !
      integer :: i, j2, iorb, kapi, kapf, kappa, kappa_max, kappa_min, n_max, &
                 ni, nf, norb, norb_sub, forb
      !
      print *, "Initialize the calculation of angular coefficients for ..."
      !
      ! Count the number of (n-kappa-m) orbitals and allocate memory for the
      ! orbital reference list
      norb = 0
      do  i = 1,csf_set%nwshells
         norb_sub = angular_momentum_j(csf_set%subshell(i)%kappa) + 1
         norb     = norb + norb_sub
      end do
      !
      pair%number_of_orbitals   = norb
      pair%number_of_electrons  = csf_set%number_of_electrons       
      allocate( pair%orbital(1:norb) )
      !
      ! Now define the reference list of the (n-kappa-m) orbitals
      norb = 0
      do  iorb = 1,csf_set%nwshells
         j2 = angular_momentum_j(asf_initial%csf_set%subshell(iorb)%kappa)
         do  i = -j2,j2,2
            norb = norb + 1
            pair%orbital(norb)%n     = asf_initial%csf_set%subshell(iorb)%n
            pair%orbital(norb)%kappa = asf_initial%csf_set%subshell(iorb)%kappa
            pair%orbital(norb)%mm    = i
         end do
      end do
      !
      print 2,  pair%number_of_orbitals
      print 3,  (iorb,pair%orbital(iorb)%n,pair%orbital(iorb)%kappa,	  &
                      pair%orbital(iorb)%mm,iorb=1,pair%number_of_orbitals)
      !
     2 format( /1x,"The orbital reference list is defined by ",i3, &
                      " orbital functions.",                          &
                  /1x,"they are given in the format  i) n kappa m   ..." /)
     3 format(     3( i4, ") " ,i3,1x,i3,1x,i3, "/2   " ) )
      !
      !
      det_orbitals%noint               = 0
      det_orbitals%nod                 = 0
      det_orbitals%nodmax              = 0
      det_orbitals%norbital            = pair%number_of_orbitals
      det_orbitals%number_of_electrons = csf_set%number_of_electrons
      allocate ( det_orbitals%orbital(1:det_orbitals%norbital) )
      !
      do  i = 1,pair%number_of_orbitals
         det_orbitals%orbital(i) = pair%orbital(i)
      end do
      !
      ! Calculate the overlap matrix of the initial- and final-state orbitals
      ! first determine kappa_min, kappa_max, and n_max
      kappa_min = 1000;   kappa_max = -1000; n_max = 0
      do  iorb = 1,pair%number_of_orbitals
         if (pair%orbital(iorb)%n > n_max) then
	    n_max = pair%orbital(iorb)%n 
	 end if
         if (pair%orbital(iorb)%kappa < kappa_min) then
            kappa_min = pair%orbital(iorb)%kappa 
	 end if
         if (pair%orbital(iorb)%kappa > kappa_max) then
            kappa_max = pair%orbital(iorb)%kappa 
	 end if
      end do
      !
      allocate ( nonorth_overlap(kappa_min:kappa_max,1:n_max,1:n_max) )
      nonorth_overlap(:,:,:) = zero
      !
      do  iorb = 1,csf_set%nwshells 
         do  forb = 1,csf_set%nwshells
	    kapi  = csf_set%subshell(iorb)%kappa
	    ni    = csf_set%subshell(iorb)%n
	    kapf  = csf_set%subshell(forb)%kappa
	    nf    = csf_set%subshell(forb)%n
	    !
	    if (kapi /= kapf) cycle
	    !
	    if (wave_i%rwf(iorb)%orbital%n     /= ni	.or.   &
	        wave_i%rwf(iorb)%orbital%kappa /= kapi)     then
               stop "nonorth_initialize(): program stop A."
	    end if
	    !
	    if (wave_f%rwf(forb)%orbital%n     /= nf	.or.	 &
	        wave_f%rwf(forb)%orbital%kappa /= kapf) then
               stop "nonorth_initialize(): program stop B."
	    end if
	    !
	    nonorth_overlap(kapi,nf,ni) =                                    &
	           rk_integral_grasp2k_ab(wave_f,wave_i,0,forb,iorb)
	 end do
      end do
      !
      call nonorth_print_overlap(kappa_min,kappa_max,n_max)
      !
      print *, "   ... initialization complete."
      print *, " "
      !
   end subroutine nonorth_initialize
   !
   !
   subroutine nonorth_initialize_ci(csf_set,wave_f,wave_i,csp)
   !--------------------------------------------------------------------
   ! Initializes the computation of the 'generalized angular coefficients'
   ! if it is carried out for a single CSF scheme but with one orbital
   ! in the continuum. It defines the 'reference
   ! list' of orbitals for the data structure pair% and calculates all
   ! the necessary one-electron orbitals. 
   ! 
   ! Calls: angular_momentum_j(), rk_integral_grasp2k_ab().
   !--------------------------------------------------------------------
      !
      type(csf_basis), intent(in)        :: csf_set
      type(grasp2k_orbital), intent(in)  :: wave_i, wave_f
      type(orbital_function), intent(in) :: csp
      !
      integer :: i, j2, iorb, kapi, kapf, kappa, kappa_max, kappa_min, n_max, &
                 ni, nf, norb, norb_sub, forb
      !
      print *, "Initialize the calculation of angular coefficients for ..."
      !
      ! Count the number of (n-kappa-m) orbitals and allocate memory for the
      ! orbital reference list
      norb = 0
      do  i = 1,csf_set%nwshells
         norb_sub = angular_momentum_j(csf_set%subshell(i)%kappa) + 1
         norb     = norb + norb_sub
      end do
      !
      pair%number_of_orbitals   = norb
      pair%number_of_electrons  = csf_set%number_of_electrons       
      allocate( pair%orbital(1:norb) )
      !
      ! Now define the reference list of the (n-kappa-m) orbitals
      norb = 0
      do  iorb = 1,csf_set%nwshells
         j2 = angular_momentum_j(csf_set%subshell(iorb)%kappa)
         do  i = -j2,j2,2
            norb = norb + 1
            pair%orbital(norb)%n     = csf_set%subshell(iorb)%n
            pair%orbital(norb)%kappa = csf_set%subshell(iorb)%kappa
            pair%orbital(norb)%mm    = i
	    !!x print *, "iorb, norb, n, kappa, mm = ",iorb, norb, &
	    !!x  pair%orbital(norb)%n,pair%orbital(norb)%kappa,pair%orbital(norb)%mm
         end do
      end do
      !
      print 2,  pair%number_of_orbitals
      print 3,  (iorb,pair%orbital(iorb)%n,pair%orbital(iorb)%kappa,	  &
                      pair%orbital(iorb)%mm,iorb=1,pair%number_of_orbitals)
      !
     2 format( /1x,"The orbital reference list is defined by ",i3, &
                      " orbital functions.",                          &
                  /1x,"they are given in the format  i) n kappa m   ..." /)
     3 format(     3( i4, ") " ,i3,1x,i3,1x,i3, "/2   " ) )
      !
      !
      det_orbitals%noint               = 0
      det_orbitals%nod                 = 0
      det_orbitals%nodmax              = 0
      det_orbitals%norbital            = pair%number_of_orbitals
      det_orbitals%number_of_electrons = csf_set%number_of_electrons
      allocate ( det_orbitals%orbital(1:det_orbitals%norbital) )
      !
      do  i = 1,pair%number_of_orbitals
         det_orbitals%orbital(i) = pair%orbital(i)
      end do
      !
      ! Calculate the overlap matrix of the initial- and final-state orbitals
      ! first determine kappa_min, kappa_max, and n_max
      kappa_min = 1000;   kappa_max = -1000; n_max = 0
      do  iorb = 1,pair%number_of_orbitals
         if (pair%orbital(iorb)%n > n_max) then
	    n_max = pair%orbital(iorb)%n 
	 end if
         if (pair%orbital(iorb)%kappa < kappa_min) then
            kappa_min = pair%orbital(iorb)%kappa 
	 end if
         if (pair%orbital(iorb)%kappa > kappa_max) then
            kappa_max = pair%orbital(iorb)%kappa 
	 end if
      end do
      !
      allocate ( nonorth_overlap(kappa_min:kappa_max,-1:n_max,-1:n_max) )
      nonorth_overlap(:,:,:) = zero
      !
      do  forb = 1,csf_set%nwshells - 1
         do  iorb = 1,csf_set%nwshells - 1
	    kapi  = csf_set%subshell(iorb)%kappa
	    ni    = csf_set%subshell(iorb)%n
	    kapf  = csf_set%subshell(forb)%kappa
	    nf    = csf_set%subshell(forb)%n
	    !
	    if (kapi /= kapf) cycle
	    !
	    if (wave_i%rwf(iorb)%orbital%n     /= ni	.or.   &
	        wave_i%rwf(iorb)%orbital%kappa /= kapi)     then
               stop "nonorth_initialize_cont(): program stop A."
	    end if
	    !
	    if (wave_f%rwf(forb)%orbital%n     /= nf	.or.	 &
	        wave_f%rwf(forb)%orbital%kappa /= kapf) then
               stop "nonorth_initialize_cont(): program stop B."
	    end if
	    !
	    nonorth_overlap(kapi,nf,ni) =                                    &
	           rk_integral_grasp2k_ab(wave_f,wave_i,0,forb,iorb)
	 end do
      end do
      !
      kapf  = csf_set%subshell(csf_set%nwshells)%kappa
      nf    = csf_set%subshell(csf_set%nwshells)%n
      !
      do  iorb = 1,csf_set%nwshells - 1
         kapi  = csf_set%subshell(iorb)%kappa
         ni    = csf_set%subshell(iorb)%n
         !
         if (kapi /= kapf) cycle
         !
         if (wave_i%rwf(iorb)%orbital%n     /= ni    .or.   &
             wave_i%rwf(iorb)%orbital%kappa /= kapi)	     then
            stop "nonorth_initialize_cont(): program stop C."
         end if
         !
         if (csp%orbital%n     /= nf	.or.	 &
             csp%orbital%kappa /= kapf) then
            stop "nonorth_initialize_cont(): program stop D."
         end if
         !
         nonorth_overlap(kapi,nf,ni) =  			      &
                 rk_integral_grasp2k_cd(csp,wave_i%rwf(iorb),0) / (one)
         !
         print "(a,3i5,2x,2(1pe16.8))",                               &
	    " *** overlap: kapi, nf, ni, ovl,ovl' = ",kapi, nf, ni,   &
        	  nonorth_overlap(kapi,nf,ni),                        &
		  rk_integral_grasp2k_cd(wave_i%rwf(iorb),csp,0)
      end do
      !
      call nonorth_print_overlap(kappa_min,kappa_max,n_max)
      !
      print *, "   ... initialization complete."
      print *, " "
      !
   end subroutine nonorth_initialize_ci
   !
   !
   subroutine nonorth_initialize_fi()
   !--------------------------------------------------------------------
   ! Initializes the computation of the 'generalized angular coefficients'
   ! for the main program xnonorthogonal, i.e. if they start from two
   ! given set of initial- and final-state CSF. It defines the 'reference
   ! list' of orbitals for the data structure pair% and calculates all
   ! the necessary one-electron orbitals. 
   ! 
   ! Calls: angular_momentum_j(), rk_integral_grasp2k_ab().
   !--------------------------------------------------------------------
      !
      integer :: i, iorb, j2, forb, kapi, kapf, kappa, kappa_max, kappa_min, &
                 ni, nf, n_max, norb, norb_sub
      !
      print *, "Initialize the calculation of angular coefficients for ..."
      !
      ! Define a reference list of (n-kappa-m) orbitals for the Slater
      ! determinants; first check that the orbitals in the same order
      ! for the initial- and final-state CSF
      if (asf_initial%csf_set%number_of_electrons /= &
          asf_final%csf_set%number_of_electrons) then
         print *, " "
         stop "nonorth_initialize_fi(): program stop A."
      else if (asf_initial%csf_set%nwshells /= asf_final%csf_set%nwshells) then
         print *, " "
         stop "nonorth_initialize_fi(): program stop B."
      end if
      !
      do  i = 1,asf_initial%csf_set%nwshells
         if (asf_initial%csf_set%subshell(i) ==                              &
	     asf_final%csf_set%subshell(i)) then
	 else
            print *, " "
            stop "nonorth_initialize_fi(): program stop C."
	 end if
      end do
      !
      ! Count the number of (n-kappa-m) orbitals and allocate memory for the
      ! orbital reference list
      norb = 0
      do  i = 1,asf_initial%csf_set%nwshells
         norb_sub = angular_momentum_j(asf_initial%csf_set%subshell(i)%kappa)+1
         norb     = norb + norb_sub
      end do
      !
      pair%number_of_orbitals   = norb
      pair%number_of_electrons  = asf_initial%csf_set%number_of_electrons       
      allocate( pair%orbital(1:norb) )
      !
      ! Now define the reference list of the (n-kappa-m) orbitals
      norb = 0
      do  iorb = 1,asf_initial%csf_set%nwshells
         j2 = angular_momentum_j(asf_initial%csf_set%subshell(iorb)%kappa)
         do  i = -j2,j2,2
            norb = norb + 1
            pair%orbital(norb)%n     = asf_initial%csf_set%subshell(iorb)%n
            pair%orbital(norb)%kappa = asf_initial%csf_set%subshell(iorb)%kappa
            pair%orbital(norb)%mm    = i
         end do
      end do
      !
      write(24,2) pair%number_of_orbitals
      write(24,3) (iorb,pair%orbital(iorb)%n,pair%orbital(iorb)%kappa,&
                        pair%orbital(iorb)%mm,iorb=1,pair%number_of_orbitals)
      !
      print 2,    pair%number_of_orbitals
      print 3,    (iorb,pair%orbital(iorb)%n,pair%orbital(iorb)%kappa,      &
                        pair%orbital(iorb)%mm,iorb=1,pair%number_of_orbitals)
      !
     2 format( /1x,"The orbital reference list is defined by ",i3, &
                      " orbital functions.",                          &
                  /1x,"they are given in the format  i) n kappa m   ..." /)
     3 format(     3( i4, ") " ,i3,1x,i3,1x,i3, "/2   " ) )
      !
      !
      det_orbitals%noint               = 0
      det_orbitals%nod                 = 0
      det_orbitals%nodmax              = 0
      det_orbitals%norbital            = pair%number_of_orbitals
      det_orbitals%number_of_electrons = asf_initial%csf_set%number_of_electrons
      allocate ( det_orbitals%orbital(1:det_orbitals%norbital) )
      !
      do  i = 1,pair%number_of_orbitals
         det_orbitals%orbital(i) = pair%orbital(i)
      end do
      !
      ! Calculate the overlap matrix of the initial- and final-state orbitals
      ! first determine kappa_min, kappa_max, and n_max
      kappa_min = 1000;   kappa_max = -1000; n_max = 0
      do  iorb = 1,pair%number_of_orbitals
         if (pair%orbital(iorb)%n > n_max) then
	    n_max = pair%orbital(iorb)%n 
	 end if
         if (pair%orbital(iorb)%kappa < kappa_min) then
            kappa_min = pair%orbital(iorb)%kappa 
	 end if
         if (pair%orbital(iorb)%kappa > kappa_max) then
            kappa_max = pair%orbital(iorb)%kappa 
	 end if
      end do
      !
      allocate ( nonorth_overlap(kappa_min:kappa_max,1:n_max,1:n_max) )
      nonorth_overlap(:,:,:) = zero
      !
      do  iorb = 1,asf_initial%csf_set%nwshells 
         do  forb = 1,asf_final%csf_set%nwshells
	    kapi  = asf_initial%csf_set%subshell(iorb)%kappa
	    ni    = asf_initial%csf_set%subshell(iorb)%n
	    kapf  = asf_final%csf_set%subshell(forb)%kappa
	    nf    = asf_final%csf_set%subshell(forb)%n
	    !
	    if (kapi /= kapf) cycle
	    !
	    if (wave_initial%rwf(iorb)%orbital%n     /= ni    .or.   &
	        wave_initial%rwf(iorb)%orbital%kappa /= kapi) then
               stop "nonorth_initialize_fi(): program stop D."
	    end if
	    !
	    if (wave_final%rwf(forb)%orbital%n     /= nf    .or.     &
	        wave_final%rwf(forb)%orbital%kappa /= kapf) then
               stop "nonorth_initialize_fi(): program stop E."
	    end if
	    !
	    nonorth_overlap(kapi,nf,ni) =                                    &
	           rk_integral_grasp2k_ab(wave_final,wave_initial,0,forb,iorb)
	 end do
      end do
      !
      call nonorth_print_overlap(kappa_min,kappa_max,n_max)
      !
      print *, "   ... initialization complete."
      print *, " "
      !
   end subroutine nonorth_initialize_fi
   !
   !
   subroutine nonorth_open_vnu(file_formatted)
   !--------------------------------------------------------------------
   ! Opens a .vnu  Angular Coefficient output File on stream 29. 
   !--------------------------------------------------------------------
      !
      logical, intent(in) :: file_formatted
      !
      integer :: ierr
      character(len=256) :: nonorth_vnu_file
      !
    1 print *, "Enter a file name for the  nonorth.vnu  file:"
      read *,  nonorth_vnu_file
      if (file_formatted) then
	 call file_open(29,nonorth_vnu_file,"formatted  ","new",ierr)
      else
	 call file_open(29,nonorth_vnu_file,"unformatted","new",ierr)
      end if
      !
      if (ierr /= 0) goto 1
      !
      if (file_formatted) then
	 write(29,*) "NON-ORTHOGONAL"
      else
	 write(29) "NON-ORTHOGONAL"
      end if
      !
   end subroutine nonorth_open_vnu
   !
   !
   subroutine nonorth_print_overlap(kappa_min,kappa_max,n_max)
   !--------------------------------------------------------------------
   ! Prints the overlap integrals from the initialization procedure
   ! in a neat format.
   !
   ! Calls: angular_momentum_string(), orbital_symmetry().
   !--------------------------------------------------------------------
      !
      integer, intent(in)           :: kappa_min, kappa_max, n_max
      !
      integer                       :: i, nf, ni, kappa, no_int
      integer, dimension(100)       :: pqn_ff, pqn_ii
      real(kind=dp), dimension(100) :: value
      !
      !
      ! Print the overlap integrals
      do  kappa = kappa_min,kappa_max
         i = 0
	 !
         do  nf = 1,n_max
            do  ni = 1,n_max
	       if (nonorth_overlap(kappa,nf,ni) /= zero) then
	          i = i + 1
		  pqn_ff(i) = nf;   pqn_ii(i) = ni
		  value(i)  = nonorth_overlap(kappa,nf,ni)
	       end if
	    end do
	 end do
	 !
	 if (i > 0) then
	    no_int = i
            write(*,4) orbital_symmetry(kappa),                             &
                      (orbital_name(pqn_ff(i),kappa),                       &
                       orbital_name(pqn_ii(i),kappa),value(i)," ",          &
                       i=1,no_int-1),                                       &
                       orbital_name(pqn_ff(no_int),kappa),                  &
                       orbital_name(pqn_ii(no_int),kappa),value(i)
            if (no_int < 3) write(*,*) "  ------------ "
	 end if
      end do
      !
    4 format(/" + ",a," symmetry: ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /   "   ------------ ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /   "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /   "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /   "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /   "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /   "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /   "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a),  &
             /   "                ",3(2x,"<",a,"|",a,"> = ",1p,e11.4,2x,a))
      !
   end subroutine nonorth_print_overlap
   !
end module rabs_nonorthonormal
















