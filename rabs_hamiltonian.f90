module rabs_hamiltonian
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module comprises all global data and routines for the set-up and
! diagonalization of the Hamiltonian matrix. The 'physical interactions'
! which are incorporated into this matrix are defined by various logical
! flags, see below. -- This module also maintains a data storage for the
! effective (one- and two-particle) interaction strengths; this storage
! has first to be initialized before it can be used in course of the
! calculation. This storage looks up tables whether these effective 
! interaction strengths have been calculated before and 'add' new data
! to the storage as they occur. Thus, the initialization of the storage
! need not to determine which data are indeed needed.
!-----------------------------------------------------------------------
   !
   use rabs_anco
   use rabs_constant
   use rabs_csl
   use rabs_nucleus
   use rabs_XL
   implicit none
   !
   public  :: hamiltonian_find_blocks
                 ! Returns the number of blocks as well as the angular 
                 ! momentum, parity and the total number of CSF in each block.
   public  :: hamiltonian_full_me
                 ! Returns a single (full) matrix elements of the Hamiltonian.
   public  :: hamiltonian_half_column
                 ! Return a (lower-half) column of the Hamiltonian matrix.
   public  :: hamiltonian_initialize
                 ! Initializes the calculation and storage of the effective
                 ! interaction strengths for the current choice of the
                 ! electron-electron interaction and additional contributions,
                 ! ir required.
   private :: hamiltonian_matrix_element
                 ! Calculates a single element of the Hamiltonian matrix from
                 ! a list of angular coefficients.
   public  :: hamiltonian_pack_column
                 ! Packs a column of the Hamiltonian into a 'reduced' vector of
                 ! nonzero mixing coefficients of type(hamiltonian_column).
   public  :: hamiltonian_set_column
                 ! Calculates and returns one column (column) of the 
		 ! Hamiltonian matrix for a given (J,parity) block.
   public  :: hamiltonian_set_matrix
                 ! Set up the Hamiltonian matrix for a given (J,parity) block.
   public  :: hamiltonian_set_vacpol
                 ! Set up the vacuum polarization potential for the given 
                 ! nuclear charge distribution at each grid point.
   private :: hamiltonian_set_vacpol_2
                 ! Set up the second-order vacuum polarization potential.
   private :: hamiltonian_set_vacpol_4
                 ! Set up the fourth-order vacuum polarization potential.
   public  :: hamiltonian_T0_storage
                 ! Returns the one-particle effective interaction strength
                 ! of zero rank.   
   public  :: hamiltonian_unpack_column
                 ! Unpacks a 'reduced' column of the Haamiltonian into its
                 ! original vector structure.
   private :: hamiltonian_vacpol_Kn
                 ! Evaluates the K_N(X) functions of Fullerton and Rinker in
                 ! order to set up the vacuum polarization potential.
   private :: hamiltonian_vacpol_Lk
                 ! Evaluates the L_K(X) functions of Fullerton and Rinker in
                 ! order to set up the vacuum polarization potential.
   public  :: hamiltonian_XL_storage
                 ! Returns the effective interaction strength X^L(abcd)
                 ! for a given set of interacting orbitals  and type of the
                 ! interaction using a storage of these integrals.   
                 !  
   ! Storage for the initial and final atomic states and wave functions
   type(asf_basis), public       :: asf_bound
   type(grasp2k_orbital), public :: wave_bound
   !
   integer, public, parameter :: hamiltonian_nobit      = bit_size(1),    &
                                 hamiltonian_fullmatrix = 8
   !
   integer, public :: hamiltonian_noblock, hamiltonian_iblock,            &
                      hamiltonian_no_eigenpairs, hamiltonian_max_eigenpair 
   !
   integer, dimension(:), pointer, public :: hamiltonian_eigenpair
   !
   ! Define a derived structure to keep information about all (J,parity) blocks
   type, public :: hamiltonian_block
      integer           :: nocsf   ! Number of CSF in this (J,parity) block.
      integer(kind=i1b) :: totalJ
      character(len=1)  :: parity
      integer, dimension(:), pointer         :: csf_ndx
      real(kind=dp), dimension(:), pointer   :: eigenvalue
      real(kind=dp), dimension(:,:), pointer :: eigenvector
   end type hamiltonian_block
   !
   type(hamiltonian_block), dimension(:), pointer :: H_block
   !
   ! Define a derived structure for one (reduced) column of the Hamiltonian
   ! matrix
   type, public :: hamiltonian_column
      integer   :: column     ! The column of the current (J,parity) block.
      integer   :: non_zero_me
      integer, dimension(:), pointer       :: ndx
      real(kind=dp), dimension(:), pointer :: me
   end type hamiltonian_column
   !
   type(hamiltonian_column), dimension(:), pointer :: H_matrix
   !
   integer, parameter, private :: XL_max_number = 32000000
   integer, private            :: XL_number
   type, private :: T0storage
      real(kind=dp), dimension(:), pointer :: me
   end type T0storage
   !
   ! Define a derived structure for storing the effective one- and two-particle
   ! interaction strengths
   type, private :: XLstorage
      !! integer(kind=i2b), dimension(:,:,:,:), pointer   :: ndx_bcdk
      integer, dimension(:,:,:,:), pointer   :: ndx_bcdk
   end type XLstorage
   !
   type(T0storage), dimension(:), allocatable, private :: T0_store
   type(XLstorage), dimension(:), allocatable, private :: XL_index
   real(kind=dp), dimension(:), allocatable ,  private :: XL_store
   !
   ! Define global logical flags for the control of the RELCI program; the
   ! default values for these flags may be overwritten interactively during 
   ! input time
   logical, public :: hamiltonian_normal_ms                 = .false.,  &
                      hamiltonian_specific_ms               = .false.,  &
                      hamiltonian_vacuum_pol                = .false.,  &
                      hamiltonian_self_energy               = .false.,  &
                      hamiltonian_use_memory                = .true.,   &
                      hamiltonian_use_disc                  = .false.,  &
                      hamiltonian_recalculate_alltime       = .false.,  &
		      hamiltonian_nlarger                   = .false.,  &
		      hamiltonian_damped_ci                 = .false.,  &
                      hamiltonian_use_storage               = .true.,   &
                      hamiltonian_use_full_diagonal         = .false.,  &
                      hamiltonian_use_grasp2k               = .true.,   &
                      hamiltonian_XL_coulomb                = .true.,   &
                      hamiltonian_XL_debye                  = .false.,  &
                      hamiltonian_XL_gaunt                  = .false.,  &
                      hamiltonian_XL_breit0                 = .false.,  &
                      hamiltonian_XL_tbreit                 = .false.
   !
   ! Define several counter to enable a summary of the computation
   integer, public :: hamiltonian_oneint_ener_comput        = 0, &
                      hamiltonian_oneint_nms_comput         = 0, &
                      hamiltonian_oneint_vp_comput          = 0, &
                      hamiltonian_one_me_stored             = 0, &
                      hamiltonian_one_me_reused             = 0, &
                      hamiltonian_twoint_Xk_comput          = 0, &
                      hamiltonian_twoint_Xk_reused          = 0, &
                      hamiltonian_twoint_Xk_stored          = 0, &
                      hamiltonian_H_me_calculated           = 0, &
                      hamiltonian_H_me_refered              = 0, &
                      hamiltonian_H_me                      = 0, &
                      hamiltonian_nlarger_zero              = 1000
   !
   real(kind=dp), public :: hamiltonian_cutoff = eps10, hamiltonian_density, &
                            hamiltonian_average, hamiltonian_debye_lambda,   &
			    hamiltonian_ci_damping
   !
   ! Define logical flags for debugging individual procedures
   logical, public :: debug_find_blocks                     = .false.,   &
                      debug_print_Hmatrix                   = .false.,   &
                      debug_set_vacpol_4                    = .false.
   !
contains
   !
   subroutine hamiltonian_find_blocks(csf_set)
   !--------------------------------------------------------------------
   ! Determines and set up the individual blocks of the Hamiltonian matrix.
   ! See the derived type(H_block) for the definition of a single
   ! (J,parity) block of the Hamiltonian matrix.
   !--------------------------------------------------------------------
      !
      type(csf_basis), intent(inout) :: csf_set
      !
      integer                   :: i, ndx, nocsf, no_eigenvectors
      integer, dimension(0:100) :: block_m, block_p, ndx_m, ndx_p
      !
      block_m(:) = 0;   block_p(:) = 0;   ndx_m(:) = 0; ndx_p(:) = 0
      !
      do  i = 1,csf_set%nocsf
         if (csf_set%csf(i)%parity == "+"   .and. &
             csf_set%csf(i)%totalJ <= 100) then
            block_p(csf_set%csf(i)%totalJ) = block_p(csf_set%csf(i)%totalJ) + 1
         else if (csf_set%csf(i)%parity == "-"   .and. &
                  csf_set%csf(i)%totalJ <= 100) then
            block_m(csf_set%csf(i)%totalJ) = block_m(csf_set%csf(i)%totalJ) + 1
         else if (rabs_use_stop) then
            stop "hamiltonian_find_blocks(): program stop A."
         end if   
      end do 
      !
      hamiltonian_noblock = 0;   nocsf = 0;   
      do i = 0,100
         if (block_m(i) > 0) then
            hamiltonian_noblock = hamiltonian_noblock + 1
            ndx_m(i)            = hamiltonian_noblock
            nocsf               = nocsf + block_m(i)
         end if
         if (block_p(i) > 0) then
            hamiltonian_noblock = hamiltonian_noblock + 1
            ndx_p(i)            = hamiltonian_noblock
            nocsf               = nocsf + block_p(i)
         end if
      end do
      !
      if (rabs_use_stop   .and.  nocsf /= csf_set%nocsf) then
         print *, "nocsf, csf_set%nocsf = ",nocsf, csf_set%nocsf
         stop "hamiltonian_find_blocks(): program stop B."
      end if
      !
      ! Allocate memory for all blocks as well as for the eigenvalues and
      ! eigenvectors of each block
      allocate ( H_block(1:hamiltonian_noblock) )
      do i = 0,100
         if (block_m(i) > 0) then
             ndx = ndx_m(i)
             H_block(ndx)%nocsf  = block_m(i)
             H_block(ndx)%totalJ = i
             H_block(ndx)%parity = "-"
             allocate( H_block(ndx)%csf_ndx(1:block_m(i)) )
             no_eigenvectors = min(hamiltonian_no_eigenpairs,H_block(ndx)%nocsf)
             allocate( H_block(ndx)%eigenvalue(1:no_eigenvectors) )
             allocate( H_block(ndx)%eigenvector(1:no_eigenvectors,            &
                                                1:H_block(ndx)%nocsf) )
         end if
         if (block_p(i) > 0) then
             ndx = ndx_p(i)
             H_block(ndx)%nocsf  = block_p(i)
             H_block(ndx)%totalJ = i
             H_block(ndx)%parity = "+"
             allocate( H_block(ndx)%csf_ndx(1:block_p(i)) )
             no_eigenvectors = min(hamiltonian_no_eigenpairs,H_block(ndx)%nocsf)
             allocate( H_block(ndx)%eigenvalue(1:no_eigenvectors) )
             allocate( H_block(ndx)%eigenvector(1:no_eigenvectors,            &
                                                1:H_block(ndx)%nocsf) )
         end if
      end do
      !
      ! Store the CSF indices into this structure
      block_m(:) = 0;   block_p(:) = 0
      do  i = 1,csf_set%nocsf
         if (csf_set%csf(i)%parity == "+") then
            ndx = ndx_p(csf_set%csf(i)%totalJ)
            block_p(ndx) = block_p(ndx) + 1
            H_block(ndx)%csf_ndx(block_p(ndx)) = i
         else if (csf_set%csf(i)%parity == "-") then
            ndx = ndx_m(csf_set%csf(i)%totalJ)
            block_m(ndx) = block_m(ndx) + 1
            H_block(ndx)%csf_ndx(block_m(ndx)) = i
         else if (rabs_use_stop) then
            stop "hamiltonian_find_blocks(): program stop C."
         end if   
      end do 
      !
      if (debug_find_blocks) then
         write(99,*) " "
         write(99,*) "Number of diagonal (J,parity) blocks in H = ", &
                     hamiltonian_noblock
         do  i = 1,hamiltonian_noblock
            write(99,*) " "
            write(99,*) "Block ",i,"   (2*totalJ,parity = ",         &
                        H_block(i)%totalJ, H_block(i)%parity,")"
            write(99,*) "nocsf_block = ",H_block(i)%nocsf
            write(99,"(2x,20i6)") H_block(i)%csf_ndx(1:H_block(i)%nocsf)
         end do
      end if
      !
   end subroutine hamiltonian_find_blocks
   !
   !
   function hamiltonian_full_me(i,k)                          result(me)
   !--------------------------------------------------------------------
   ! Returns the (full) Hamiltonian matrix element H(i,k) for the currently
   ! considered (J,parity) block. 
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: i, k
      real(kind=dp)       :: me
      !
      integer :: ii, kk, m, no_T_coeff, no_V_coeff, rr, ss
      !
      !!x print *, "cccc"
      !!x print *, "cccc"
      me = zero
      !
      if (i >= k) then
         ii = i;   kk = k
      else
         ii = k;   kk = i
      end if
      !
      if (hamiltonian_use_memory) then
         me = zero
         do  m = 1,H_matrix(kk)%non_zero_me
            if (H_matrix(kk)%ndx(m) <  ii) cycle
            if (H_matrix(kk)%ndx(m) == ii) then
               me = H_matrix(kk)%me(m)
               exit
            end if
         end do
         hamiltonian_H_me_refered = hamiltonian_H_me_refered + 1
      else if (hamiltonian_use_disc .or. hamiltonian_recalculate_alltime) then
         ss = H_block(hamiltonian_iblock)%csf_ndx(kk)
         rr = H_block(hamiltonian_iblock)%csf_ndx(ii)
         !
         ! Calculate angular coefficients for CSF pair (rr,ss)
         call anco_calculate_csf_pair(asf_bound%csf_set,rr,ss, &
                                      no_T_coeff,no_V_coeff)
         if (no_T_coeff > 0   .or.   no_V_coeff > 0) then
            ! Calculate the Hamiltonian matrix element
            me = hamiltonian_matrix_element(anco_T_list,no_T_coeff, &
                                            anco_V_list,no_V_coeff,wave_bound)
         end if 
         hamiltonian_H_me_calculated = hamiltonian_H_me_calculated + 1
      else if (rabs_use_stop) then
         stop "hamiltonian_full_me(): program stop B."
      end if                      
      !
      hamiltonian_H_me = hamiltonian_H_me + 1
      !
   end function hamiltonian_full_me
   !
   !
   subroutine hamiltonian_half_column(j,n,ndx,me,non_zero)   
   !--------------------------------------------------------------------
   ! Returns the (lower half) column j of the Hamiltonian matrix, i.e.
   ! H(j,j), ... H(n,j) for the currently considered (J,parity) block.
   ! The upper-half column H(1,j), ... H(j-1,j) is set to zero.
   ! The non-zero indices in this column are returned in ndx(1:non_zero)
   ! along with their actual values.
   !--------------------------------------------------------------------
      !
      integer, intent(in)                      :: j, n
      integer, intent(out)                     :: non_zero
      integer, dimension(n), intent(out)       :: ndx
      real(kind=dp), dimension(n), intent(out) :: me
      !
      integer                  :: column, i, no_T_coeff, no_V_coeff, r, rr, ss
      real(kind=dp)            :: hme
      !
      if (hamiltonian_use_memory) then
         non_zero        = H_matrix(j)%non_zero_me
         ndx(1:non_zero) = H_matrix(j)%ndx(1:non_zero)
         me(1:non_zero)  = H_matrix(j)%me(1:non_zero)
      else if (hamiltonian_use_disc) then
       1 read(27,end=2) column, non_zero
         if (column < j) then
            read(27) 
            read(27) 
            goto 1
         else if (column == j) then
            do  i = 1,non_zero
               read(27) ndx(i), me(i)
            end do
            hamiltonian_H_me_refered = hamiltonian_H_me_refered + n - j + 1
            return
         else if (column > j) then
            goto 2
         end if
         !
       2 rewind 27
         goto 1
      else if (hamiltonian_recalculate_alltime) then
         ss = H_block(hamiltonian_iblock)%csf_ndx(j)
         do  r = j,n
            rr = H_block(hamiltonian_iblock)%csf_ndx(r)
            !
            ! Calculate angular coefficients for CSF pair (rr,ss)
            call anco_calculate_csf_pair(asf_bound%csf_set,rr,ss, &
                                               no_T_coeff,no_V_coeff)
            if (no_T_coeff > 0   .or.   no_V_coeff > 0) then
               ! Calculate the Hamiltonian matrix element
               hme = hamiltonian_matrix_element(anco_T_list,no_T_coeff,     &
                                           anco_V_list,no_V_coeff,wave_bound)
               if (abs(hme) > hamiltonian_cutoff) then
                  non_zero      = non_zero + 1
                  ndx(non_zero) = r
                  me(non_zero)  = hme
                  hamiltonian_H_me_calculated = hamiltonian_H_me_calculated + 1
               end if
            end if 
         end do
      else if (rabs_use_stop) then 
         stop "hamiltonian_half_column(): program stop A."
      end if                      
      !
      hamiltonian_H_me = hamiltonian_H_me + n - j + 1
      !
   end subroutine hamiltonian_half_column
   !
   !
   subroutine hamiltonian_initialize(wave)
   !--------------------------------------------------------------------
   ! Initializes a proper amount of storage for the one-electron interactions
   ! and the effective interaction strengths X^L (abcd). If  
   ! hamiltonian_use_storage = .true., these interactions are stored owing 
   ! to the indices of the orbital functions (in wave) and due to the 
   ! rank L for the effective X^L (abcd). 
   !
   ! One-electron interals: 
   ! ----------------------
   !   Initialization is T0_store(ia)%me(ib) = 1.e20_dp; a value which is
   !   smaller has been calculated before for the corresponding integral
   !   T(ab) = I(ab) + VP(ab) + ...
   !
   ! Two-electron effective strengths X^L (abcd):
   ! --------------------------------------------
   !   Initialization is  XL_index(ia)%ndx_bcd(ib,ic,id) = 0; 
   !   for XL_index(ia)%ndx_bcd(ib,ic,id) < 0, the interaction strength 
   !   is 'zero' (effectively < eps10);
   !   for XL_index(ia)%ndx_bcd(ib,ic,id) > 0, it gives the index i in
   !   XL_store(i) where the interaction strength is stored. 
   !--------------------------------------------------------------------
      !
      type(grasp2k_orbital), intent(in) :: wave
      !
      integer :: i, nrwf, nxl
      !
      if (hamiltonian_use_storage) then
         !!x print *, "aa: wave%number_of_rwf = ",wave%number_of_rwf
         nrwf = wave%number_of_rwf
         !
         ! Allocate the T0 storage
         allocate( T0_store(1:nrwf) )
         do  i = 1,nrwf
            allocate( T0_store(i)%me(1:i) )
            T0_store(i)%me(1:i) = 1.0e20
         end do
         !
         ! Allocate the XL storage
         nxl  = XL_max_number
         nxl  = min( nxl, nrwf*nrwf*nrwf*(nrwf+1)*7/2 )
         allocate( XL_index(1:nrwf),  XL_store(1:nxl) )
         XL_store(:) = zero;   XL_number = 0
         do  i = 1,nrwf
            allocate( XL_index(i)%ndx_bcdk(1:i,1:nrwf,1:nrwf,0:6) )
            XL_index(i)%ndx_bcdk(1:i,1:nrwf,1:nrwf,0:6) = 0
         end do
         !
         print *, "Storage initialization for the one- and two--partice " //&
                  "effective interaction strengths complete."
      end if
      !
   end subroutine hamiltonian_initialize
   !
   !
   function hamiltonian_matrix_element(tlist,nt,vlist,nv,wave) result(me)
   !--------------------------------------------------------------------
   ! Calculates a single element of the Hamiltonian matrix from the 
   ! lists of angular coefficients tlist(1:nt) and vlist(1:nv).
   !--------------------------------------------------------------------
      ! 
      integer, intent(in)                          :: nt, nv
      type(anco_T_coeff), dimension(:), intent(in) :: tlist
      type(anco_V_coeff), dimension(:), intent(in) :: vlist
      type(grasp2k_orbital), intent(in)            :: wave
      real(kind=dp)                                :: me
      !
      integer       :: i
      real(kind=dp) :: aweight
      !
      me = zero
      !
      ! Calculate contributions from T coefficients
      do  i = 1,nt
         aweight = tlist(i)%T * sqrt(angular_momentum_j(tlist(i)%a%kappa)+one)
         me = me + aweight * hamiltonian_T0_storage(tlist(i)%a,tlist(i)%b,wave)
      end do
      !
      ! Calculate contributions from V coefficients
      do  i = 1,nv                  
         if (hamiltonian_nlarger)  then
            if (vlist(i)%a%n > hamiltonian_nlarger_zero  .or. &
                vlist(i)%b%n > hamiltonian_nlarger_zero  .or. &
                vlist(i)%c%n > hamiltonian_nlarger_zero  .or. &
                vlist(i)%d%n > hamiltonian_nlarger_zero)   then
               aweight = zero  
            else
               aweight = vlist(i)%V
            end if
         else
            aweight = vlist(i)%V
         end if
         !
         me = me + aweight * hamiltonian_XL_storage(vlist(i)%nu,             &
                             vlist(i)%a,vlist(i)%b,vlist(i)%c,vlist(i)%d,wave)
      end do  
      !
   end function hamiltonian_matrix_element
   !
   !
   subroutine hamiltonian_pack_column(s,ntot,full_column,bit_is_nonzero, &
                                                                   red_me)
   !--------------------------------------------------------------------
   ! Packs a column of the Hamiltonian matrix into a reduced vector 
   ! structure. Packing is done for all non-zero elements from 
   ! s ... ntot. This procedure is currently not in use here.
   !--------------------------------------------------------------------
      !
      integer, intent(in)                      :: ntot, s
      real(kind=dp), dimension(:), intent(in)  :: full_column
      integer, dimension(:), intent(out)       :: bit_is_nonzero
      real(kind=dp), dimension(:), intent(out) :: red_me
      !
      integer :: count_nonzero, i, int, noint, pos
      !
      ! Clear up the integer storage of the non-zero matrix elements
      noint = size(bit_is_nonzero)
      do  int = 1,noint;   do  pos = 0,hamiltonian_nobit-1
         bit_is_nonzero(int) = ibclr(bit_is_nonzero(int),pos)
      end do;   end do
      !
      int   = 0; count_nonzero = 0
      do  i = 1,ntot-s+1
         if (mod(i-1,hamiltonian_nobit) == 0) then
            int = int + 1
         end if
         pos = mod(i,hamiltonian_nobit) - 1
         if (pos == -1) pos = hamiltonian_nobit -1
         !
         if (abs(full_column(i+s-1)) > eps10) then
            bit_is_nonzero(int) = ibset(bit_is_nonzero(int),pos)
            count_nonzero = count_nonzero + 1
            red_me(count_nonzero) = full_column(i+s-1) 
         end if
      end do
      !
      if (rabs_use_stop   .and.  count_nonzero /= size(red_me)) then
         stop "hamiltonian_pack_column(): program stop A."
      end if
      !
   end subroutine hamiltonian_pack_column
   !
   !
   subroutine hamiltonian_set_column(hblock,column,csf_set,wave,non_zero,ndx,me)
   !--------------------------------------------------------------------
   ! Calculates and returns one column (column) of the Hamiltonian matrix 
   ! for a given (J,parity) block. At output, nme is the number of non-zero 
   ! matrix elements in the vector full_column(column:hblock%nocsf).
   !--------------------------------------------------------------------
      !
      integer, intent(in)                      :: column
      type(hamiltonian_block), intent(in)      :: hblock
      type(csf_basis), intent(in)              :: csf_set
      type(grasp2k_orbital), intent(in)        :: wave
      !
      integer, intent(out)                     :: non_zero
      integer, dimension(:), intent(out)       :: ndx
      real(kind=dp), dimension(:), intent(out) :: me
      !
      real(kind=dp)                            :: H_rr_ss
      !
      integer :: no_T_coeff, no_V_coeff, r, rr, ss 
      !
      ss = hblock%csf_ndx(column)
      non_zero = 0
      do  r = column,hblock%nocsf
         !
         rr = hblock%csf_ndx(r)
         !
         ! Calculate angular coefficients for CSF pair (rr,ss)
         call anco_calculate_csf_pair(csf_set,rr,ss,no_T_coeff,no_V_coeff)
         if (no_T_coeff > 0   .or.   no_V_coeff > 0) then
            ! Calculate the Hamiltonian matrix element
            H_rr_ss = hamiltonian_matrix_element(anco_T_list,no_T_coeff, &
                                                 anco_V_list,no_V_coeff,wave)
            if (abs(H_rr_ss) > hamiltonian_cutoff) then
               non_zero      = non_zero + 1
               ndx(non_zero) = r
               me(non_zero)  = H_rr_ss
            end if
         end if 
      end do
      !
   end subroutine hamiltonian_set_column
   !
   !
   subroutine hamiltonian_set_matrix(hblock,csf_set,wave)
   !--------------------------------------------------------------------
   ! Set up the Hamiltonian matrix for a given (J,parity) block of the
   ! Hamiltonian matrix.
   !--------------------------------------------------------------------
      !
      type(hamiltonian_block), intent(in)      :: hblock
      type(csf_basis), intent(in)              :: csf_set
      type(grasp2k_orbital), intent(in)        :: wave
      !
      integer                                  :: non_zero
      integer, dimension(1:hblock%nocsf)       :: ndx
      real(kind=dp)                            :: H_rr_ss
      real(kind=dp), dimension(1:hblock%nocsf) :: me
      !
      integer :: no_T_coeff, no_V_coeff, r, rr, s, ss 
      !
      do  s = 1,hblock%nocsf
         ss = hblock%csf_ndx(s)
         non_zero = 0
         do  r = s,hblock%nocsf
            !
            rr = hblock%csf_ndx(r)
            !
            ! Calculate angular coefficients for CSF pair (rr,ss)
            call anco_calculate_csf_pair(csf_set,rr,ss,no_T_coeff,no_V_coeff)
            if (no_T_coeff > 0   .or.   no_V_coeff > 0) then
               !
               ! Calculate the Hamiltonian matrix element
               H_rr_ss = hamiltonian_matrix_element(anco_T_list,no_T_coeff, &
                                                    anco_V_list,no_V_coeff,wave)
               if (abs(H_rr_ss) > hamiltonian_cutoff) then
                  non_zero      = non_zero + 1
                  ndx(non_zero) = r
		  !
                  if (hamiltonian_damped_ci   .and.  r /= s)  then
                     me(non_zero)  = H_rr_ss * hamiltonian_ci_damping
		  else
                     me(non_zero)  = H_rr_ss
		  end if
               end if
            end if 
         end do
         !
         H_matrix(s)%non_zero_me = non_zero
         allocate( H_matrix(s)%ndx(1:non_zero), H_matrix(s)%me(1:non_zero) )
         H_matrix(s)%ndx(1:non_zero) = ndx(1:non_zero)
         H_matrix(s)%me(1:non_zero)  = me(1:non_zero)
      end do
      !
   end subroutine hamiltonian_set_matrix
   !
   !
   subroutine hamiltonian_set_vacpol(mtp)
   !--------------------------------------------------------------------
   ! Controls the setting up of the vacuum polarization potentials for 
   ! the given nuclear charge distribution at each grid point using the 
   ! analytic functions defined by  L Wayne Fullerton and G A Rinker Jr 
   ! in Phys Rev A 13 (1976) 1283. The potential is accumulated in 
   ! array vacpol_2_grasp2k(1:n_grasp2k).
   ! This routine is similar to vacpol from GRASP92 [RCI92] which was 
   ! written by F A Parpia; it has been adapted here to the Fortran90/95
   ! standard.
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: mtp
      !
      ! Redefine zdist_grasp2k to be  rho*r*r'
      zdist_grasp2k(1:mtp) = zdist_grasp2k(1:mtp) * r_grasp2k(1:mtp) * &
                             rp_grasp2k(1:mtp)
      !
      ! Second-order vacuum polarisation potential; returned in array 
      ! vacpol_2_grasp2k
      call hamiltonian_set_vacpol_2(mtp)
      !
      ! Fourth-order vacuum polarisation potential; returned in array 
      ! vacpol_4_grasp2k
      call hamiltonian_set_vacpol_4(mtp)
      !
   end subroutine hamiltonian_set_vacpol
   !
   !
   subroutine hamiltonian_set_vacpol_2(mtp)
   !--------------------------------------------------------------------
   ! Sets up the second-order vacuum polarization potential using  
   ! equations (1) and (4) of L Wayne Fullerton and  G A Rinker, Jr,  
   ! Phys Rev A  13 (1976) 1283-1287. The potential is accumulated in  
   ! array  vacpol_2_grasp2k(1:n_grasp2k).
   ! This routine is similar to vac2 from GRASP92 [RCI92] which was 
   ! written by F A Parpia; it has been adapted to the Fortran90/95
   ! standard.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: mtp
      !
      integer             :: i, k
      real(kind=dp)       :: epsi, factor, tbi, twocv, &
                             x, xi, xm, xk, xp
      real(kind=dp), dimension(1:n_grasp2k+10) :: ta
      !
      ! Overall initialization
      epsi  = eps10*eps10*eps10
      twocv = c_vacuum +c_vacuum
      !
      ! Potential for a point nucleus: equation (1)
      ! (this is also the asymptotic form for a finite nucleus)
      factor = -(two * nuclear_charge)/(three * pi * c_vacuum)
      vacpol_2_grasp2k(1) = zero
      !
      i   = 1
    1 i   = i + 1
      x   = twocv * r_grasp2k(i)
      tbi = (factor/r_grasp2k(i)) * hamiltonian_vacpol_Kn(x,1)
      if (abs(tbi) >= epsi) then
         vacpol_2_grasp2k(i) = tbi
         if (i < n_grasp2k) go to 1
      else
         vacpol_2_grasp2k(i:n_grasp2k) = zero
      end if
      !
      ! Potential for a finite nucleus: equation (4)
      if (nuclear_model(1:5) == "fermi") then
         factor = -two / (three * c_vacuum**2)
         !
         ! Set up integrand
         vacpol_2_grasp2k(1) = zero
         !
         k     = 1
       3 k     = k+1
         xk    = twocv * r_grasp2k(k)
         ta(1) = zero
         do  i = 2,mtp
            xi = twocv * r_grasp2k(i)
            xm = abs(xk - xi)
            xp = xk + xi
            ta(i) = (hamiltonian_vacpol_Kn(xm,0) -  &
                     hamiltonian_vacpol_Kn(xp,0)) * zdist_grasp2k(i)
         end do
         !
         x = quad_grasp2k(ta,mtp)
         x = x * factor / r_grasp2k(k)
         !
         ! Get out of loop if the asymptotic value has been attained
         if (abs(x) >= epsi) then
            if (abs((x-vacpol_2_grasp2k(k))/x) > 1.0e-5_dp) then
               vacpol_2_grasp2k(k) = x
               if (k < n_grasp2k) go to 3
            end if
         end if
      end if
      !
   end subroutine hamiltonian_set_vacpol_2
   !
   !
   subroutine hamiltonian_set_vacpol_4(mtp)
   !--------------------------------------------------------------------
   ! Sets up the fourth-order vacuum polarization potential using 
   ! equations (11) and (12) of L Wayne Fullerton and  G A Rinker, Jr,  
   ! Phys  Rev  A 13 (1976) 1283-1287. The potential is accumulated in 
   ! array  vacpol_4_grasp2k(1:n_grasp2k).
   ! This routine is similar to vac4 from GRASP92 [RCI92] which was 
   ! written by F A Parpia; it has been adapted here to the Fortran90/95
   ! standard.
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: mtp
      !
      integer             :: i, ii, ii1, ii2, k, nb2, nrows
      real(kind=dp)       :: epsi, factor, tci, twocv, x, xk, xi, xm, xp
      real(kind=dp), dimension(1:n_grasp2k+10) :: ta
      !
      ! Overall initialization
      epsi  = eps10*eps10*eps10
      twocv = c_vacuum + c_vacuum
      !
      ! Potential for point nucleus: equation (12)
      factor = -nuclear_charge / (pi*c_vacuum)**2
      !
      vacpol_4_grasp2k(1) = zero
      i   = 1
    1 i   = i+1
      x   = twocv * r_grasp2k(i)
      tci = (factor/r_grasp2k(i)) * hamiltonian_vacpol_Lk(x,1)
      if (abs(tci) >= epsi) then
         vacpol_4_grasp2k(i) = tci
         if (i < n_grasp2k) goto 1
      else
         vacpol_4_grasp2k(i:n_grasp2k) = zero
      end if
      !
      ! Potential for finite nucleus: equation (11)
      if (nuclear_model(1:5) == "fermi") then
         factor = -one / (pi*c_vacuum**3)
         vacpol_4_grasp2k(1) = zero
         !
         k     = 1
       3 k     = k+1
         xk    = twocv * r_grasp2k(k)
         ta(1) = zero
         !
         do  i = 2,mtp
            xi = twocv * r_grasp2k(i)
            xm = abs(xk - xi)
            xp = xk + xi
            ta(i) = (hamiltonian_vacpol_Lk(xm,0) -                  &
                     hamiltonian_vacpol_Lk(xp,0) ) * zdist_grasp2k(i)
         end do
         !
         x = quad_grasp2k(ta,mtp)
         x = x * factor / r_grasp2k(k)
         !
         ! Get out of the loop if the asymptotic region has been reached
         if (abs(x) >= epsi) then
            if (abs((vacpol_4_grasp2k(k)-x)/x) > 1.0e-3_dp) then
               vacpol_4_grasp2k(k) = x
               if (k < n_grasp2k) goto 3
            end if
         end if
      end if
      !
      if (debug_set_vacpol_4) then
         write (99,5)
       5 format( /// " ++++++++++ hamiltonian_set_vacpol_4 ++++++++++",  &
          // 2(" -------- r -------- ----- VV2 (r) ----- ----- VV4 (r) -----"))
         nb2 = n_grasp2k / 2
         if (2 * nb2 == n_grasp2k) then
            nrows = nb2
         else
            nrows = nb2 + 1
         end if
         do  ii = 1,nrows
            ii1 = ii;   ii2 = ii1 + nrows
            if (ii2 <= n_grasp2k) then
               write(99,6) r_grasp2k(ii1),vacpol_2_grasp2k(ii1),  &
                                          vacpol_4_grasp2k(ii1),  &
                           r_grasp2k(ii2),vacpol_2_grasp2k(ii2),  &
                                          vacpol_4_grasp2k(ii2)
            else if (ii1 <= n_grasp2k) then
               write(99,6) r_grasp2k(ii1),vacpol_2_grasp2k(ii1),  &
                                          vacpol_4_grasp2k(ii1)
             6 format(1p,6(1x,1d19.12))
            end if
         end do
      end if
      !
      ! Generate total vacuum-polarization potential
      vacpol_2_grasp2k(1:n_grasp2k) = vacpol_4_grasp2k(1:n_grasp2k) + &
                                      vacpol_2_grasp2k(1:n_grasp2k)
      !
   end subroutine hamiltonian_set_vacpol_4
   !
   !
   function hamiltonian_T0_storage(a,b,wave)                  result(T0)
   !--------------------------------------------------------------------
   ! Returns the effective one-particle interaction (strengths) T^0 (ab)
   ! for orbital quantum numbers a and b. 
   ! For hamiltonian_use_storage = .true., this procedure looks up
   ! the storage and, if it needed, calculates T^0 (ab); this is then
   ! stored for later use.
   ! The interactions which are included in T^0 (ab) are derived from 
   ! logical variables like hamiltonian_vacuum_pol, ...
   ! It is assumed that these variables do not change during the use of 
   ! this storage management; no test is made whether the T^0 (ab) 
   ! in the storage have been calculated with the same logical flags.
   ! This routine is useful if T^0 (ab) are needed again and again like
   ! in relativistic CI calculation or MBPT.
   !--------------------------------------------------------------------
      !
      implicit none
      type(nkappa), intent(in)          :: a, b
      type(grasp2k_orbital), intent(in) :: wave
      real(kind=dp)                     :: T0
      !
      integer :: i, ia, iaa, ib, ibb
      !
      ia = 0;   ib = 0 
      do  i = 1,wave%number_of_rwf
         !!x print *, "wave%rwf(i)%orbital%n = ",wave%rwf(i)%orbital%n
         if (wave%rwf(i)%orbital%n == a%n    .and. &
             wave%rwf(i)%orbital%kappa == a%kappa) then
            ia = i
         end if
         if (wave%rwf(i)%orbital%n == b%n    .and. &
             wave%rwf(i)%orbital%kappa == b%kappa) then
            ib = i
         end if
      end do
      !
      if (rabs_use_stop   .and.  &
         (ia == 0  .or.  ib == 0)) then
         stop "hamiltonian_T0_storage(): program stop A."
      else if (ia < ib)  then
         iaa = ib;   ibb = ia
      else
         iaa = ia;   ibb = ib
      end if
      !
      if (hamiltonian_use_storage) then
         if (T0_store(iaa)%me(ibb) < 1.0e19) then
            T0 = T0_store(iaa)%me(ibb)
            hamiltonian_one_me_reused = hamiltonian_one_me_reused + 1
            return
         end if
      end if
      !
      ! If not yet returned, calculate the T^0 (ab)
      T0 = zero
      !
      ! The one-particle energies are always added
      if (hamiltonian_nlarger) then
         if (a%n > hamiltonian_nlarger_zero  .or.  &
             b%n > hamiltonian_nlarger_zero)       then
            if (iaa == ibb) then
               T0 = wave%rwf(iaa)%energy
            else
               T0 = zero
            end if
         else
            T0 = I_ab_grasp2k(wave%rwf(iaa),wave%rwf(ibb))
         end if
      else
         ! Standard route in xrelci
         T0 = I_ab_grasp2k(wave%rwf(iaa),wave%rwf(ibb))
         hamiltonian_oneint_ener_comput = hamiltonian_oneint_ener_comput + 1
      end if
      !
      ! Add vacuum polarization contributions
      if (hamiltonian_vacuum_pol) then
         T0 = T0 + vpintf_grasp2k(wave%rwf(iaa),wave%rwf(ibb)) 
         hamiltonian_oneint_vp_comput = hamiltonian_oneint_vp_comput + 1
      end if
      !
      ! Add normal specific mass shift contributions
      if (hamiltonian_normal_ms) then
         T0 = T0 + (one*electron_mass_in_amu/nuclear_mass) *                              &
                   I_ab_grasp2k(wave%rwf(iaa),wave%rwf(ibb),"kinetic")
         print *, "K_ab = ",I_ab_grasp2k(wave%rwf(iaa),wave%rwf(ibb),"kinetic")
         hamiltonian_oneint_nms_comput = hamiltonian_oneint_nms_comput + 1
      end if
      !
      ! Store the result if appropriate
      if (hamiltonian_use_storage) then
         T0_store(iaa)%me(ibb)     = T0
         hamiltonian_one_me_stored = hamiltonian_one_me_stored + 1
      end if
      !
   end function hamiltonian_T0_storage
   !
   !
   subroutine hamiltonian_unpack_column(s,ntot,bit_is_nonzero,red_me,  &
                                                            full_column)
   !--------------------------------------------------------------------
   ! Unpacks a column of the Hamiltonian matrix from its reduced vector 
   ! structure. Packing is done for all non-zero elements from 
   ! s ... ntot.
   !--------------------------------------------------------------------
      !
      integer, intent(in)                      :: ntot, s
      integer, dimension(:), intent(in)        :: bit_is_nonzero
      real(kind=dp), dimension(:), intent(in)  :: red_me
      real(kind=dp), dimension(:), intent(out) :: full_column
      !
      integer :: count_nonzero, i, int, pos
      !
      ! Clear up the vector storage of the full_column
      full_column(:) = zero
      !
      int   = 0; count_nonzero = 0
      do  i = s,ntot
         if (mod(i-s,hamiltonian_nobit) == 0)  int = int + 1
         pos = mod(i-s,hamiltonian_nobit)
         if (btest(bit_is_nonzero(int),pos)) then
            count_nonzero = count_nonzero + 1
            full_column(i) = red_me(count_nonzero)
         end if
      end do
      !
      if (rabs_use_stop   .and.  count_nonzero /= size(red_me)) then
         print *, "count_nonzero, size(red_me) = ",count_nonzero, size(red_me)
         stop "hamiltonian_unpack_column(): program stop A."
      end if  
      !
   end subroutine hamiltonian_unpack_column
   !
   !
   function hamiltonian_vacpol_Kn(x,n)                        result(Kn)
   !--------------------------------------------------------------------
   ! Evaluates the K_N(X) functions using the analytic functions defined 
   ! in tables 1 and 3 of Fullerton and Rinker, Phys  Rev  A 13 (1976) 
   ! 1283-1287.
   ! This routine is similar to  from GRASP92 [RCI92] which was 
   ! written by F A Parpia; it has been adapted here to the Fortran90/95
   ! standard.
   !--------------------------------------------------------------------
      !
      integer, intent(in)       :: n
      real(kind=dp), intent(in) :: x
      real(kind=dp)             :: Kn
      !
      real(kind=dp), dimension(10,4), parameter :: p = reshape(source =   & 
        (/  8.8357293375e-1_dp, -2.8259817381e-1_dp, -5.8904879578e-1_dp, &
            1.2500133434e-1_dp, -3.2729913852e-2_dp,  8.2888574511e-3_dp, &
           -1.0327765800e-5_dp,  6.3643668900e-5_dp,  0.0_dp,             &
            0.0_dp,                                                       &
           -7.1740181754e-1_dp,  1.1780972274_dp,    -3.7499963087e-1_dp, &
            1.3089675530e-1_dp, -3.8258286439e-2_dp, -2.4297287300e-5_dp, &
           -3.5920148670e-4_dp, -1.7170090700e-5_dp,  0.0_dp,             &
            0.0_dp,                                                       &
            9.9999999987e-1_dp,  1.9770200000e-8_dp, -7.5000050190e-1_dp, &
            7.8540306316e-1_dp, -3.4988601655e-1_dp,  6.4596333000e-5_dp, &
           -9.8189080747e-3_dp,  8.6513145800e-5_dp, -2.3969236620e-4_dp, &
            0.0_dp,                                                       &
            6.0000000002_dp,    -6.4305200000e-8_dp,  2.1049413000e-6_dp, &
           -2.6711271500e-5_dp, -1.3705236152e-1_dp, -6.3476104090e-4_dp, &
           -7.8739801501e-2_dp, -1.9641740173e-3_dp, -3.4752369349e-3_dp, &
           -7.3145316220e-4_dp  /), shape=(/10,4/))
      !
      real(kind=dp), dimension(2,4), parameter :: b = reshape(source =    &
        (/ -3.19999594323e+2_dp,    2.53900995981_dp,                     &
           -6.40514843293e+1_dp,    7.11722714285e-1_dp,                  &
            5.19010136460e+3_dp,    8.28495496200e+1_dp,                  &
            3.18150793824e+2_dp,    4.33898867347e+1_dp  /), shape=(/2,4/))
      !
      real(kind=dp), dimension(3,4), parameter :: c = reshape(source =    &
        (/ -3.19999594333e+2_dp,    2.53901020662_dp,                     &
            0.0_dp,                 6.40514843287e+1_dp,                  &
           -7.11722686403e-1_dp,    8.04220774800e-4_dp,                  &
            2.76805406060e+4_dp,   -3.27039477790e+2_dp,                  &
            0.0_dp,                 8.48402116837e+2_dp,                  &
           -2.56939867765e+1_dp,    3.20844906346e-1_dp  /), shape=(/3,4/))
      !
      real(kind=dp), dimension(5,4), parameter :: d = reshape(source =    &
        (/  5.018065179_dp,         7.1518912620e+1_dp,                   &
            2.116209929e+2_dp,      3.1403274780e+1_dp,                   &
           -1.0_dp,                 2.1723864090e+2_dp,                   &
            1.643364528e+3_dp,      2.1222445120e+3_dp,                   &
           -4.512004044e+1_dp,      1.0_dp,                               &
            8.540770444_dp,         6.0762427660e+1_dp,                   &
            9.714630584e+1_dp,      3.1549735930e+1_dp,                   &
            1.0_dp,                 5.9243015865e-1_dp,                   &
            2.0596312871_dp,        3.7785190424_dp,                      &
            3.5614853214_dp,        1.0_dp               /), shape=(/5,4/))
      !
      real(kind=dp), dimension(5,4), parameter :: e = reshape(source =    &
        (/  2.669207401_dp,         5.172549669e+1_dp,                    &
            2.969809720e+2_dp,      5.364324164e+2_dp,                    &
            1.535335924e+2_dp,      1.155589983e+2_dp,                    &
            1.292191441e+3_dp,      3.831198012e+3_dp,                    &
            2.904410075e+3_dp,      0.0_dp,                               &
            4.543392478_dp,         3.514920169e+1_dp,                    &
            6.019668656e+1_dp,      8.468839579_dp,                       &
            0.0_dp,                 3.1511867816e-1_dp,                   &
            3.473245222e-1_dp,      3.8791936870e-2_dp,                   &
           -1.3059741497e-3_dp,     0.0_dp               /), shape=(/5,4/))
      !
      integer, dimension(4), parameter :: np = (/ 8, 8, 9, 10 /)
      !
      integer       :: i, k, nn
      real(kind=dp) :: bsum, csum, dsum, esum, sum, x2, xm, xn
      !
      if (x == zero) goto 11
      if (rabs_use_stop   .and.                               &
         (n < 0   .or.   n == 2   .or.   n == 4   .or.  n > 5)) then
         print *, "Attempt to calculate FUNK (X,N) for N other than "// &
                  "0, 1, 3 and 5."  
         stop     "hamiltonian_vacpol_Kn(): program stop A."   
      end if
      !
      select case(n-3)
      case(:-1)
         k  = n+1
         xn = one
      case(0)
         k  = n
         xn = one / (x**2)
      case(1:)
         k  = n-1
         xn = one / (x**4)
      case default
         stop "hamiltonian_vacpol_Kn(): program stop B."   
      end select
      if (x > one) goto 9
      !
      ! Calculate function for x < = 1
      nn  = np(k)
      sum = zero
      do  i = 1,nn
         sum = sum + p(i,k) * xn
         xn  = xn * x
      end do
      !
      x2   = x * x
      bsum = b(1,k) + x2*(b(2,k) + x2)
      csum = c(1,k) + x2*(c(2,k) + x2*c(3,k))
      !
      select case(k)
      case(1)
         bsum = bsum * x
      case(2,4)
      case(3)
         bsum = bsum * x2
      case default
         stop "hamiltonian_vacpol_Kn(): program stop C."   
      end select
      sum = sum+bsum*log (x)/csum
      !
      Kn = sum
      return
      !
      ! Calculate function for x > 1
    9 xn   = one
      dsum = zero
      esum = zero
      do  i = 1,5
         dsum = dsum + d(i,k) * xn
         esum = esum + e(i,k) * xn
         xn   = xn / x
      end do
      xm  = -x
      sum = dsum * exp(xm) / (esum*sqrt(x**3))
      !
      Kn  = sum
      return
      !
   11 if (rabs_use_stop   .and.   n /= 0) then
         print *, "Attempt to calculate Kn (0,N) for N > 0."
         stop     "hamiltonian_vacpol_Kn(): program stop D."   
      end if
      !
      Kn  = p(1,1)
      !
   end function hamiltonian_vacpol_Kn
   !
   !
   function hamiltonian_vacpol_Lk(x,k)                        result(Lk)
   !--------------------------------------------------------------------
   ! Evaluates the LK(X) functions using the analytic functions defined  
   ! in table 5  and equations (20) and  (21) of Fullerton and Rinker, 
   ! Phys  Rev  A 13 (1976) 1283-1287. 
   ! This routine is similar to  from GRASP92 [RCI92] which was 
   ! written by F A Parpia; it has been adapted here to the Fortran90/95
   ! standard.
   !--------------------------------------------------------------------
      !
      integer, intent(in)       :: k
      real(kind=dp), intent(in) :: x
      real(kind=dp)             :: Lk
      !
      real(kind=dp), dimension(6,2), parameter :: f = reshape(source =      &
        (/ 2.008188_dp,    -2.397605_dp,    1.046471_dp,    -3.670660e-1_dp,&
           6.374000e-2_dp, -3.705800e-2_dp, 1.646407_dp,    -2.092942_dp,   &
           9.623100e-1_dp, -2.549600e-1_dp, 1.644040e-1_dp,  0.0_dp         &
                                                           /), shape=(/6,2/))
      !
      real(kind=dp), dimension(3,2), parameter :: g = reshape(source =      &
        (/   7.51198e-1_dp,  1.38889e-1_dp,  2.0886e-2_dp,                  &
             1.37691e-1_dp, -4.16667e-1_dp, -9.7486e-2_dp                   &   
                                                           /), shape=(/3,2/))
      !
      real(kind=dp), dimension(2,2), parameter :: h = reshape(source =      &
        (/  -4.44444e-1_dp, -3.472e-3_dp,                                   &
             4.44444e-1_dp,  1.7361e-2_dp                  /), shape=(/2,2/))
      !
      real(kind=dp), parameter :: a = 2.2_dp,   b = -1.72_dp
      !
      integer       :: i, k1
      real(kind=dp) :: sum, sumg, sumh, x2, xm, xn
      !
      if (rabs_use_stop .and.  &
         (k < 0   .or.   k > 1)) then
         print *, "K must be either 0 or 1."
         stop     "hamiltonian_vacpol_Lk(): program stop A."
      end if
      if (x  > two)  goto 3
      if (x == zero) goto 6
      !
      ! Use rational approximation for x < 2
      k1  = k + 1
      sum = zero
      xn  = one
      do  i = 1,6
         sum = sum + xn * f(i,k1)
         xn  = xn * x
      end do
      x2   = x * x
      sumg = g(1,k1) + x2 * (g(2,k1) + x2*g(3,k1))
      sumh = h(1,k1) + x2 * x2 * h(2,k1)
      xn   = log(x)
      sumg = xn*(sumg + xn*sumh)
      if (k == 0) goto 2
      sum  = sum + sumg
      goto 7
    2 sum  = sum + x * sumg
      goto 7
      !
    3 continue
      sum  = a + b / x
      if (k == 0) goto 4
      sum  = sum + (sum + b/x) / x
    4 sum  = sum / x
      xm   = -x
      sum  = sum * exp(xm)
      goto 7
    6 if (rabs_use_stop   .and.  k == 1) then
         print *, "Attempt to calculate function for zero argument "// &
                  "and K value of 1."
         stop     "hamiltonian_vacpol_Lk(): program stop B."
      end if
      sum  = f(1,1)
    7 Lk   = sum
      !
   end function hamiltonian_vacpol_Lk
   !
   !
   function hamiltonian_XL_storage(L,a,b,c,d,wave)            result(XL)
   !--------------------------------------------------------------------
   ! Returns the effective interaction strengths X^L (abcd) for a given 
   ! type of 'interaction' and for given rank L and orbital quantum numbers 
   ! a, b, c, and d. For XL_use_storage = .true., this procedure looks up
   ! the storage and, if it needs to calculate X^L (abcd), it is stored for 
   ! later use.
   ! The mode of the 'interaction' are derived from the logical variables
   ! hamiltonian_XL_coulomb,  hamiltonian_XL_breit0, ... It is assumed 
   ! that these variables do not change during the use of the storage; 
   ! no test is made whether the X^L (abcd) in the storage have been 
   ! calculated with the same logical flags.
   ! This routine is useful if X^L (abcd) are needed again and again like
   ! in relativistic CI calculation of MBPT.
   ! 
   ! For hamiltonian_specific_ms = .true., it also adds the contributions
   ! from the specific mass shift to the total effective strength.
   !
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)               :: L
      type(nkappa), intent(in)          :: a, b, c, d
      type(grasp2k_orbital), intent(in) :: wave
      real(kind=dp)                     :: XL
      !
      integer       :: i, ia, iaa, ib, ibb, ic, icc, id, idd, index
      real(kind=dp) :: wa
      !
      ia = 0;   ib = 0;   ic = 0;   id = 0   
      do  i = 1,wave%number_of_rwf
         if (wave%rwf(i)%orbital%n == a%n    .and. &
             wave%rwf(i)%orbital%kappa == a%kappa) then
            ia = i
         end if
         if (wave%rwf(i)%orbital%n == b%n    .and. &
             wave%rwf(i)%orbital%kappa == b%kappa) then
            ib = i
         end if
         if (wave%rwf(i)%orbital%n == c%n    .and. &
             wave%rwf(i)%orbital%kappa == c%kappa) then
            ic = i
         end if
         if (wave%rwf(i)%orbital%n == d%n    .and. &
             wave%rwf(i)%orbital%kappa == d%kappa) then
            id = i
         end if
      end do
      !
      if (rabs_use_stop   .and.                                &
         (ia == 0  .or.  ib == 0  .or.  ic == 0  .or.  id == 0)) then
         !!! stop "XL_strength_fix(): program stop A."
         print *,  "XL_strength_fix(): program stop A."
      else if (ia < ib)  then
         iaa = ib;   ibb = ia;   icc = id;   idd = ic
      else
         iaa = ia;   ibb = ib;   icc = ic;   idd = id
      end if
      !
      if (hamiltonian_use_storage   .and.   L <= 6) then
         !
         if (XL_index(iaa)%ndx_bcdk(ibb,icc,idd,L) > 0) then
            index = XL_index(iaa)%ndx_bcdk(ibb,icc,idd,L)
            XL = XL_store(index)
            hamiltonian_twoint_Xk_reused = hamiltonian_twoint_Xk_reused + 1
            return
         else if (XL_index(iaa)%ndx_bcdk(ibb,icc,idd,L) < 0) then
            XL = zero
            hamiltonian_twoint_Xk_reused = hamiltonian_twoint_Xk_reused + 1
            return
         end if                      
      end if
      !
      ! If not yet returned, calculate the X^L (abcd) 
      XL = zero
      if (hamiltonian_XL_coulomb) then
      if (wave%rwf(ia)%orbital%n == 77 .or. wave%rwf(ib)%orbital%n == 77 .or.&
          wave%rwf(ic)%orbital%n == 77 .or. wave%rwf(id)%orbital%n == 77) then
         print *, "*** hamiltonian_XL_storage() modified, apr 2010 ***"
      else
         wa = XL_Coulomb_strength_grasp2k(L,wave%rwf(ia),wave%rwf(ib),       &
                                               wave%rwf(ic),wave%rwf(id))
         XL = XL + XL_Coulomb_strength_grasp2k(L,wave%rwf(ia),wave%rwf(ib), &
                                                 wave%rwf(ic),wave%rwf(id))
      end if
      end if
      !
      if (hamiltonian_XL_breit0) then
         wa = XL_Breit0_strength_grasp2k(L,wave%rwf(ia),wave%rwf(ib),       &
                                                wave%rwf(ic),wave%rwf(id))
         XL = XL + XL_Breit0_strength_grasp2k(L,wave%rwf(ia),wave%rwf(ib),  &
                                                wave%rwf(ic),wave%rwf(id))
      end if
      !
      if (hamiltonian_XL_gaunt) then
         XL = XL + XL_Gaunt_strength_grasp2k(L,wave%rwf(ia),wave%rwf(ib),   &
                                               wave%rwf(ic),wave%rwf(id))
      end if
      !
      if (hamiltonian_XL_tbreit) then
         wa = XL_Breit_strength_grasp2k(L,wave%rwf(ia),wave%rwf(ib),       &
                                               wave%rwf(ic),wave%rwf(id))
         XL = XL + XL_Breit_strength_grasp2k(L,wave%rwf(ia),wave%rwf(ib),  &
                                               wave%rwf(ic),wave%rwf(id))
         if (abs(wa) > eps10) then
            print *, "Breit: ia, ib, ic, id, wa, xl = ",ia, ib, ic, id, wa, xl
         end if
      end if
      ! 
      if (hamiltonian_XL_debye) then
         XL = XL + XL_Debye_strength_grasp2k(L,hamiltonian_debye_lambda,  &
                       wave%rwf(ia),wave%rwf(ib),wave%rwf(ic),wave%rwf(id))
      end if
      !
      if (hamiltonian_specific_ms) then     
         XL = XL - (one*electron_mass_in_amu/nuclear_mass) *                                 &
                   XL_sms_strength_grasp2k(L,wave%rwf(ia),wave%rwf(ib),   &
                                             wave%rwf(ic),wave%rwf(id))
      end if
      hamiltonian_twoint_Xk_comput = hamiltonian_twoint_Xk_comput + 1
      !
      !!x print *, "hamiltonian_XL_storage(): XL = ", XL
      !
      ! Store the result if appropriate
      if (hamiltonian_use_storage   .and.   L <= 6) then
         !
         if (abs(XL) < eps10) then
            XL_index(iaa)%ndx_bcdk(ibb,icc,idd,L) = -1
         else
            XL_number = XL_number + 1
            XL_index(iaa)%ndx_bcdk(ibb,icc,idd,L) = XL_number
            XL_store(XL_number) = XL
         end if
         hamiltonian_twoint_Xk_stored = hamiltonian_twoint_Xk_stored + 1
      end if
      !
   end function hamiltonian_XL_storage
   !
end module rabs_hamiltonian
