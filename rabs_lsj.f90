module rabs_lsj
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains the procedures which are used to control the
! jj-LS transformation of a given set of ASF in a jj-coupled basis.
!-----------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   use rabs_csl
   use rabs_file_handling
   use rabs_lsj_data
   implicit none
   !
   public   :: lsj_coefficient_LS_jj
                 ! Returns the value of the LS-jj transformation matrix
                 ! for a given set of quantum numbers.
   private  :: lsj_coefficient_LS_jj_occ_2
                 ! Returns the value of the LS-jj transformation matrix 
                 ! (l^2 LSJ| j_1 j_2 J).
   public   :: lsj_control_from_toolbox
                 ! Controls the transformation of atomic states from a jj-
                 ! into a LS-coupled CSF basis from xtoolbox.
   public   :: lsj_control_transformation
                 ! Controls the transformation of atomic states from a jj-
                 ! into a LS-coupled CSF basis.
   public   :: lsj_deallocate_asf_basis_LS
                 ! Dellocates the storage of asf_set_LS.
   public   :: lsj_form_csf_basis_LS
                 ! This subroutine fills up the variable asf_set_LS%csf_set_LS
                 ! with data generated using the one from asf_set_jj%csf_set
                 ! .....................................................
                 ! This subroutine contains the following internal subroutines:
		 !
                 !  private  :: lsj_form_csf_basis_LS_action
     	         ! The subroutine defines the "action" of subroutine 
   	         ! lsj_form_csf_basis_LS_job_count: whether it counts
   	         ! the number of csfs_LS (asf_set_LS%csf_set_LS%novcsf)
   	         ! or fills the arrays of wave functions in LS coupling
   	         ! with asf_set_LS%csf_set_LS%csf(...) with
   	         ! the corresponding quantum nubers.
		 !
                 !  private  :: lsj_form_csf_basis_LS_add_qn
    	         ! The subroutine adds quantum numbers stored 
   	         ! in temprorary arrays Li, Si, L_i, S_i, w, Q to 
   	         ! the corresponding arrays of asf_set_LS%csf_set_LS%csf().
		 !
                 !  private  :: lsj_form_csf_basis_LS_job_count
   	         ! Recursive subroutine for the calculation of the
   	         ! number of csfs_LS and corresponding quantum numbers.
		 !
                 !  private  :: lsj_form_csf_basis_LS_equiv
     	         ! This subroutine defines the "equivalency" of two csfs_jj
   	         ! in the sence of generation of csfs_LS
   	         ! number of csfs_LS and corresponding quantum numbers.
                 ! ......................................................
   public   :: lsj_get_subshell_term_LS
                 ! This procedure return all allowed subshell terms 
                 ! (l, w, Q, L, S) for given l^N which must be 0, 1, 2 or 3. 
   public   :: lsj_interprete_levels
                 ! Attempts to interpret the serial level numbers from a 
                 ! string.
   public   :: lsj_print_conf_scheme_LS
                 ! Prints all information about the single CSF scheme csf_set
                 ! in LS-coupling
   public   :: lsj_print_single_config_jj
                 ! Prints all information about the single CSF in jj-coupling
   public   :: lsj_print_single_config_LS
                 ! Prints all information about the single CSF in LS-coupling
   public   :: lsj_spectroscopic_LS
                 ! A spectroscopic notation of shell in LS coupling is return.
   public   :: lsj_transformation_ASF
                 ! Expands an atomic state function, which is represented
                 ! in a jj-coupling CSF basis into a basis of LS-coupling CSF.
   public   :: lsj_transformation_LS_jj_gen
                 !  Return the value of the transformation matrix 
                 !  from jj- to LS-coupling scheme in the case of any 
		 !  number of open shells.
   public   :: lsj_transformation_LS_jj_insade
                 !  Return the value of main part of the transformation 
		 !  matrix from jj- to LS-coupling scheme in the 
		 !  case of any number of open shells.
   public   :: lsj_triangle
                 ! Returns .true. if the lengths i1, i2, and i3 may 
		 ! form a triangle and .false. otherwise.
   !
   type, public :: cs_function_LS
      integer(kind=i1b) :: totalJ
      character(len=1)  :: parity
      integer(kind=i1b), dimension(:), pointer :: occupation
      integer(kind=i1b), dimension(:), pointer :: seniority
      integer(kind=i1b), dimension(:), pointer :: w
      integer(kind=i1b), dimension(:), pointer :: shellL
      integer(kind=i1b), dimension(:), pointer :: shellS
      integer(kind=i1b), dimension(:), pointer :: shellLX
      integer(kind=i1b), dimension(:), pointer :: shellSX
   end type cs_function_LS
   !
   type, public :: csf_basis_LS
      integer :: nocsf         ! Number of CSF in the basis.
      integer :: nwshells      ! Number of (nonrelativistic) shells.
      integer :: nwcore        ! Number of (closed) core shells.
      integer :: number_of_electrons
      type(nl), dimension(:), pointer             :: shell
      type(cs_function_LS), dimension(:), pointer :: csf
      type(parent_from_jj), dimension(:), pointer :: parent
   end type csf_basis_LS
   !
   type, public :: as_function_LS
      integer           :: level_No
      integer           :: max_csf_No
      integer(kind=i1b) :: totalL, totalS, totalJ
      character(len=1)  :: parity
      real(kind=dp)     :: energy
      real(kind=dp), dimension(:), pointer :: eigenvector
   end type as_function_LS
   !
   type, public :: asf_basis_LS
      integer         :: noasf           ! Number of considered ASF.
      real(kind=dp)   :: average_energy  ! Averaged energy of this set of ASF.
      type(as_function_LS), dimension(:), pointer :: asf
      type(csf_basis_LS) :: csf_set_LS
   end type asf_basis_LS
   !
   type, public :: parent_from_jj
      integer         :: parent_minus
      integer         :: parent_plius
   end type parent_from_jj
   !
   type, public :: lsj_list
      integer::list_size                    !number items in a list
      integer, dimension(:),pointer:: items !serial numbers of lists items
   end type lsj_list
   !
   ! Define a global jj- and LS-coupled basis
   type(asf_basis), public    :: asf_set_jj
   type(asf_basis_LS), public :: asf_set_LS
   !
   character(len=1), dimension(0:20), private, parameter :: L_string = &
      (/ "S", "P", "D", "F", "G", "H", "I", "K", "L", "M", "N", "O", "Q", &
         "R", "T", "U", "V", "W", "X", "Y", "Z" /)
   integer, public ::  debuging = 0
   !
   type, public :: lsj_level_weight
      integer	    :: csf
      real(kind=dp) :: weight
   end type lsj_level_weight
   !
   type, public :: lsj_string_csf
      integer	    :: csf
      character(len=120) :: s1, s2
   end type lsj_string_csf
   !
contains
   !
   function lsj_coefficient_LS_jj(l_shell,N,w,Q,L,S,J,jm_shell,Nm,Qm,Jm,   &
                                               jp_shell,Qp,Jp) result(wa)
   !--------------------------------------------------------------------
   ! Returns the value of the LS-jj transformation matrix for a given 
   ! set of quantum numbers.
   !
   ! Note that all (generalized) angular momentum quantum numbers except
   ! of l must be given twice the original numbers, i.e. for the quantum 
   ! numbers Q, L, S, J, jm_shell, Qm, Jm, jp_shell, Qp, Jp.
   !
   ! Calls: lsj_coefficient_LS_jj_occ_2,
   !        lsj_coefficient_LS_jj_search.
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: l_shell, N, w, Q, L, S, J,                &
                             jm_shell, Nm, Qm, Jm, jp_shell, Qp, Jp
      !
      real(kind=dp)       :: wa
      integer             :: NN, N1, N2
      !
      wa = zero
      if (jm_shell < jp_shell .or. jm_shell == jp_shell)  then
         if (l_shell > 0 .and. N > 2*l_shell +1)  then
            NN = 4*l_shell + 2 - N;   N1 = jm_shell + 1 - Nm
            N2 = jp_shell + 1 - N + Nm
         else
            NN = N;   N1 = Nm;   N2 = N - Nm;
         end if
         if (NN == 1  .or.  NN == 0)  then
            if (NN == 0 .and. N1 == 0 .and. N2 == 0) then
               wa = one
            else if (N1 == 1 .and. N2 == 0 .and. J == Jm) then
               wa = one
            else if (N1 == 0 .and. N2 == 1 .and. J == Jp) then
               wa = one
            else
               wa = zero
            end if
         else if (NN == 2)  then
            if (J > Jm + Jp  .or.  J < abs(Jm - Jp)) then
               wa = zero
            else
               if (N1 == 2 .and. N2 == 0) then
                  wa = lsj_coefficient_LS_jj_occ_2             &
		               (l_shell,L,S,J,jm_shell,jm_shell)
               else if (N1 == 1 .and. N2 == 1) then
                  wa = lsj_coefficient_LS_jj_occ_2             &
		               (l_shell,L,S,J,jm_shell,jp_shell)
               else if (N1 == 0 .and. N2 == 2) then
                  wa = lsj_coefficient_LS_jj_occ_2             &
		               (l_shell,L,S,J,jp_shell,jp_shell)
               end if;
            end if
         else if (l_shell == 1  .or.  l_shell == 2  .or.  l_shell == 3)  then
            wa = lsj_coefficient_LS_jj_search                              &
	                                (l_shell,NN,w,Q,L,S,J,N1,Qm,Jm,Qp,Jp)
         end if
      else if (l == 0)  then
         if (S == J .and. jm_shell == 1 .and. NN == N1)  then
            wa = one
         else if (S == J .and. jp_shell == 1 .and. NN == N2)  then
            wa = one
         else
            wa = zero
         end if
      else
         stop "lsj_coefficient_LS_jj(): program stop A."
      end if
      !
   end function lsj_coefficient_LS_jj
   !
   !
   function lsj_coefficient_LS_jj_occ_2(l_shell,L,S,J,jm_shell,jp_shell)  &
                                                                 result(wa)
   !--------------------------------------------------------------------
   ! Returns the value of the LS-jj transformation matrix 
   !                                                (l^2 LSJ| j_1 j_2 J).
   !
   ! Note that all (generalized) angular momentum quantum numbers except
   ! of l must be given twice the original numbers, i.e. for the quantum 
   ! numbers L, S, J, jm_shell, jp_shell.
   !
   ! Calls: wigner_9j_triangle, wigner_9j_symbol.
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: l_shell, L, S, J, jm_shell, jp_shell
      real(kind=dp)       :: wa
      !
      integer :: delta_J
      !
      wa = zero
      if (mod(L+S,4) /= 0)  then
         wa = zero
      elseif (mod(L+S+J,2) /= 0)  then
         wa = zero
      elseif (mod(jm_shell+jp_shell+J,2) /= 0)  then
         wa = zero
      else
         if (jm_shell == jp_shell) then
            if (mod(J,4) /= 0)  then
               wa = zero
            else
               wa = (jm_shell+one) * sqrt((L+one) * (S+one))
            end if
         else
            wa = sqrt(2 * (L+one) * (S+one) * (jm_shell+one) * (jp_shell+one))
         end if
         delta_J = wigner_9j_triangle(2*l_shell,2*l_shell,L,1,1,S,jm_shell, &
	                                                          jp_shell,J)
         if (delta_J /= 0)  then
            wa = wa * wigner_9j_symbol(2*l_shell,2*l_shell,L,1,1,S,jm_shell,&
	                                                   jp_shell,J,.true.)

         end if
      end if
      !
   end function lsj_coefficient_LS_jj_occ_2
   !
   !
   subroutine lsj_control_from_toolbox(csl_file,asf_weight,asf_string)
   !--------------------------------------------------------------------
   ! Controls the transformation of atomic states from a jj- into a 
   ! LS-coupled CSF basis from xtoolbox.
   !
   ! Calls: lsj_deallocate_asf_basis_LS, lsj_form_csf_basis_LS,
   ! lsj_interprete_levels, lsj_print_conf_scheme_LS,
   ! lsj_print_single_config_jj, lsj_print_single_config_LS,
   ! lsj_transformation_ASF, load_csl_from_file, load_mix_file_grasp2k, 
   ! print_configuration_scheme, file_open, util_csl_file.
   !--------------------------------------------------------------------
      !
      character(len=256), intent(in) :: csl_file
      type(lsj_level_weight), dimension(100,9), intent(out) :: asf_weight
      type(lsj_string_csf)  , dimension(:), intent(out)     :: asf_string
      !
      !
      integer            :: i, j, jj, number_of_levels, ierr, ios
      logical            :: yes, fail
      character(len=6)   :: g92mix
      character(len=256) :: record, util_csl_file
      !
      integer, dimension(1:10000)    :: levels
      integer, dimension(0:10000)    :: leading_LS
      integer, dimension(1:100)      :: iw
      integer                        :: level, nocsf_min, lev, string_length
      integer                        :: nocsf_max, sum_nocsf_min
      real(kind=dp), dimension(1:100):: weights, mcoeffs
      real(kind=dp)                  :: wa, wb
      !
      number_of_levels = asf_set_jj%noasf
      do lev = 1,number_of_levels
         levels(lev) = asf_set_jj%asf(lev)%level_No
      end do
      !
      print *, " Enter the number of the leading CSF to be printed:"
      read (*, "(i2)") nocsf_max
      call lsj_form_csf_basis_LS()
      !
      string_length = index(csl_file,'.')
      util_csl_file = csl_file(1:string_length)//'LS'
      call file_open(25,util_csl_file,"formatted  ","new",ierr)
      call lsj_print_conf_scheme_LS(25,asf_set_LS%csf_set_LS)
      close(25)
      allocate(asf_set_LS%asf(1:asf_set_jj%noasf))
      do i = 1, asf_set_jj%noasf
         allocate(asf_set_LS%asf(i)%eigenvector(1:asf_set_LS%csf_set_LS%nocsf))
         asf_set_LS%asf(i)%level_No = i
         asf_set_LS%asf(i)%energy = asf_set_jj%asf(i)%energy
         asf_set_LS%asf(i)%totalJ = asf_set_jj%asf(i)%totalJ
         asf_set_LS%asf(i)%parity = asf_set_jj%asf(i)%parity
      end do
      asf_set_LS%noasf = asf_set_jj%noasf
      jj = 0
      do  lev = 1, number_of_levels
         level = levels(lev)
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
            stop "lsj_control_transformation(): program stop A." 
         end if
         !
         nocsf_min = 5
         do j = 1,5
            if(abs(weights(j)) < 0.00001) then
               nocsf_min = j - 1
               exit
            end if
         end do
         print 16, asf_set_jj%asf(level)%level_No,                        &
         trim(angular_momentum_string(1*asf_set_jj%asf(level)%totalJ,4)), &
              asf_set_jj%asf(level)%parity,(weights(j),iw(j),j=1,nocsf_min)
         print *, " "
         print *, "Definition of leading CSF:"
         print *, " "
         call lsj_print_single_config_jj(-1,asf_set_jj%csf_set,iw(1))
         print *, " "
         print*, " Total sum over   weight    is:  ",wb
         !
         call lsj_transformation_ASF(level)
      end do
      !
      !
      !
      if (.true.) then
      !
      !
      print *, " "
      print *, " "
      print *, " "
      print *, " "
      print *, "Weights of major contributors to ASF in LS-coupling:"
      print *, " "
      print *, " Level  J Parity      CSF contributions"
      print *, " "
      !
      do  lev = 1, number_of_levels
         level = levels(lev)
         weights(1:100) = zero;   iw(1:100) = 0
         wb = zero
         do  i = 1,asf_set_LS%csf_set_LS%nocsf
            wa = asf_set_LS%asf(level)%eigenvector(i) * &
                 asf_set_LS%asf(level)%eigenvector(i)
            wb = wb + asf_set_LS%asf(level)%eigenvector(i) * &
                      asf_set_LS%asf(level)%eigenvector(i)
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
            stop "lsj_control_transformation(): program stop B." 
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
	    !
	    asf_weight(lev,j)%csf    = iw(j)
	    asf_weight(lev,j)%weight = weights(j)
         end do
         sum_nocsf_min = sum_nocsf_min + nocsf_min
         if (nocsf_min <= 0) then
            stop "lsj_control_transformation(): program stop C."
         end if
         print 16, asf_set_LS%asf(level)%level_No,                   &
                   trim(angular_momentum_string                      &
                   (1*asf_set_LS%csf_set_LS%csf(iw(1))%totalJ,4)),   &
                   asf_set_LS%csf_set_LS%csf(iw(1))%parity,          &
                   (weights(j),iw(j),J=1,nocsf_min)
         print*, "                     Total sum over  weight  is:",wb
         print *, " "
         asf_set_LS%asf(lev)%max_csf_No = iw(1)
      end do
      !
      end if
      !
      if (.false.) then
      !
      ! Print the mixing coefficients if needed; added by S.F., March 2005
      !
      print *, " "
      print *, " "
      print *, " "
      print *, " "
      print *, "Mixing coefficients of major contributors to ASF in LS-coupling:"
      print *, " "
      print *, " Level  J Parity      CSF contributions"
      print *, " "
      !
      do  lev = 1, number_of_levels
         level = levels(lev)
         weights(1:100) = zero;   mcoeffs(1:100) = zero;   iw(1:100) = 0
         wb = zero
         do  i = 1,asf_set_LS%csf_set_LS%nocsf
            wa = asf_set_LS%asf(level)%eigenvector(i) * &
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
            stop "lsj_control_transformation(): program stop Bx." 
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
            stop "lsj_control_transformation(): program stop Cx."
         end if
         print 16, asf_set_LS%asf(level)%level_No,                   &
                   trim(angular_momentum_string                      &
                   (1*asf_set_LS%csf_set_LS%csf(iw(1))%totalJ,4)),   &
                   asf_set_LS%csf_set_LS%csf(iw(1))%parity,          &
                   (mcoeffs(j),iw(j),J=1,nocsf_min)
         print*, "                     Total sum over  weight  is:",wb
         print *, " "
         asf_set_LS%asf(lev)%max_csf_No = iw(1)
      end do
      !
      end if
      !
      !
      print *, " "
      print *, "Definition of leading CSF:"
      print *, " "
      do  i = 1,asf_set_jj%csf_set%nocsf
         do  j = 1, sum_nocsf_min
            if (i == leading_LS(j)) then
               call lsj_print_single_config_LS                  &
                       (-1,asf_set_LS%csf_set_LS,leading_LS(j))
	       !
	       rewind 59
               call lsj_print_single_config_LS                  &
                       (59,asf_set_LS%csf_set_LS,leading_LS(j))
	       rewind 59
	       asf_string(i)%csf = i
	       read(59,"(a)") asf_string(i)%s1
	       read(59,"(a)") asf_string(i)%s2
	       !       
               exit
            end if
         end do
      end do
   16 format(1x,i4,1x,2a4,2x,100(3x,f8.5," of",i5))
      !
      call lsj_deallocate_asf_basis_LS(asf_set_LS)
      if (debuging /= 0 )then
         close(57)
      end if
   end subroutine lsj_control_from_toolbox
   !
   !
   subroutine lsj_control_transformation()
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
      logical            :: yes, fail
      character(len=6)   :: g92mix
      character(len=256) :: record
      character(len=256) :: util_csl_file, util_mix_file
      !
      integer, dimension(1:10000)    :: levels
      integer, dimension(0:10000)    :: leading_LS
      integer, dimension(1:100)      :: iw
      integer                        :: level, nocsf_min, lev, string_length
      integer                        :: nocsf_max, sum_nocsf_min
      real(kind=dp), dimension(1:100):: weights, mcoeffs
      real(kind=dp)                  :: wa, wb
      !
      !
      character(len=*), parameter :: text = &
         "Transform one or several ASF from a GRASP92 calculation into a " //&
         "LS-coupled CSF basis. The transformation starts from the given " //&
         ".csl and .mix files and is carried out for the n leading CSF in "//&
         "the jj-coupled basis; the new representation in the LS basis is "//&
         "printed similar as in a standard GRASP92 computation.           "
      !
      ! Write a header about the purpose of this branch of the program
      print *, text(  1: 85);    print *, text( 86: 165)
      print *, text(166: 251);   print *, text(252: 315);   print *, " "
      !
      ! Open, check, load data from, and close, the  .csl  file
    1 print *, "Enter the name of the GRASP92 configuration "//&
               "symmetry list file:"
      read (*,"(a)") util_csl_file
      if (len(trim(util_csl_file)) == 0) goto 1
      !
      call file_open(21,util_csl_file,"formatted  ","old",ierr)
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
      call load_csl_from_file(asf_set_jj%csf_set)
      close (21)
      !
    2 print *, "Enter the name of corresponding .mix mixing coefficient file:"
      read (*,"(a)") util_mix_file
      if (len(trim(util_mix_file)) == 0) goto 2
      !
      call file_open(25,util_mix_file,"unformatted","old",ierr)
      if (ierr /= 0) goto 2
      !
      ! Check the header of the file; if not as expected for unformatted
      ! files, check formatted form
      read (25,iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         close (25)
         goto 3
      else
         call load_mix_file_grasp2k(asf_set_jj,.false.,ierr)
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
    3 call file_open(25,util_mix_file,"formatted  ","old",ierr)
      if (ierr /= 0) goto 2
      read (25,"(a)",iostat=ios)  g92mix
      if (ios /= 0   .or.   g92mix /= "G92MIX") then
         print *, "ios, g92mix = ",ios, g92mix
         print *, "Not a GRASP92 Mixing Coefficients File;"
         close (25)
         goto 2
      else
         call load_mix_file_grasp2k(asf_set_jj,.true.,ierr)
         if (ierr /= 0) then
            close (25)
            goto 2
         end if
      end if
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
         print*, "There are to much ASF:", number_of_levels
         go to 4
      end if
      if (debuging /= 0 )then
         string_length = index(util_csl_file,'.')
         util_csl_file = util_csl_file(1:string_length)//'jj'
         call file_open(57,util_csl_file,"formatted  ","new",ierr)
         call print_configuration_scheme(57,asf_set_jj%csf_set)
      end if
      !
      print *, " Enter the number of the leading CSF to be printed:"
      read (*, "(i2)") nocsf_max
      call lsj_form_csf_basis_LS()
      !
      string_length = index(util_csl_file,'.')
      util_csl_file = util_csl_file(1:string_length)//'LS'
      call file_open(25,util_csl_file,"formatted  ","new",ierr)
      call lsj_print_conf_scheme_LS(25,asf_set_LS%csf_set_LS)
      close(25)
      allocate(asf_set_LS%asf(1:asf_set_jj%noasf))
      do i = 1, asf_set_jj%noasf
         allocate(asf_set_LS%asf(i)%eigenvector(1:asf_set_LS%csf_set_LS%nocsf))
         asf_set_LS%asf(i)%level_No = i
         asf_set_LS%asf(i)%energy = asf_set_jj%asf(i)%energy
         asf_set_LS%asf(i)%totalJ = asf_set_jj%asf(i)%totalJ
         asf_set_LS%asf(i)%parity = asf_set_jj%asf(i)%parity
      end do
      asf_set_LS%noasf = asf_set_jj%noasf
      jj = 0
      do  lev = 1, number_of_levels
         level = levels(lev)
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
            stop "lsj_control_transformation(): program stop A." 
         end if
         !
         nocsf_min = 5
         do j = 1,5
            if(abs(weights(j)) < 0.00001) then
               nocsf_min = j - 1
               exit
            end if
         end do
         print 16, asf_set_jj%asf(level)%level_No,                        &
         trim(angular_momentum_string(1*asf_set_jj%asf(level)%totalJ,4)), &
              asf_set_jj%asf(level)%parity,(weights(j),iw(j),j=1,nocsf_min)
         print *, " "
         print *, "Definition of leading CSF:"
         print *, " "
         call lsj_print_single_config_jj(-1,asf_set_jj%csf_set,iw(1))
         print *, " "
         print*, " Total sum over   weight    is:  ",wb
         !
         call lsj_transformation_ASF(level)
      end do
      !
      !
      !
      if (.true.) then
      !
      !
      print *, " "
      print *, " "
      print *, " "
      print *, " "
      print *, "Weights of major contributors to ASF in LS-coupling:"
      print *, " "
      print *, " Level  J Parity      CSF contributions"
      print *, " "
      !
      do  lev = 1, number_of_levels
         level = levels(lev)
         weights(1:100) = zero;   iw(1:100) = 0
         wb = zero
         do  i = 1,asf_set_LS%csf_set_LS%nocsf
            wa = asf_set_LS%asf(level)%eigenvector(i) * &
                 asf_set_LS%asf(level)%eigenvector(i)
            wb = wb + asf_set_LS%asf(level)%eigenvector(i) * &
                      asf_set_LS%asf(level)%eigenvector(i)
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
            stop "lsj_control_transformation(): program stop B." 
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
         print 16, asf_set_LS%asf(level)%level_No,                   &
                   trim(angular_momentum_string                      &
                   (1*asf_set_LS%csf_set_LS%csf(iw(1))%totalJ,4)),   &
                   asf_set_LS%csf_set_LS%csf(iw(1))%parity,          &
                   (weights(j),iw(j),J=1,nocsf_min)
         print*, "                     Total sum over  weight  is:",wb
         print *, " "
         asf_set_LS%asf(lev)%max_csf_No = iw(1)
      end do
      !
      end if
      !
      if (.false.) then
      !
      ! Print the mixing coefficients if needed; added by S.F., March 2005
      !
      print *, " "
      print *, " "
      print *, " "
      print *, " "
      print *, "Mixing coefficients of major contributors to ASF in LS-coupling:"
      print *, " "
      print *, " Level  J Parity      CSF contributions"
      print *, " "
      !
      do  lev = 1, number_of_levels
         level = levels(lev)
         weights(1:100) = zero;   mcoeffs(1:100) = zero;   iw(1:100) = 0
         wb = zero
         do  i = 1,asf_set_LS%csf_set_LS%nocsf
            wa = asf_set_LS%asf(level)%eigenvector(i) * &
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
            stop "lsj_control_transformation(): program stop Bx." 
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
            stop "lsj_control_transformation(): program stop Cx."
         end if
         print 16, asf_set_LS%asf(level)%level_No,                   &
                   trim(angular_momentum_string                      &
                   (1*asf_set_LS%csf_set_LS%csf(iw(1))%totalJ,4)),   &
                   asf_set_LS%csf_set_LS%csf(iw(1))%parity,          &
                   (mcoeffs(j),iw(j),J=1,nocsf_min)
         print*, "                     Total sum over  weight  is:",wb
         print *, " "
         asf_set_LS%asf(lev)%max_csf_No = iw(1)
      end do
      !
      end if
      !
      !
      print *, " "
      print *, "Definition of leading CSF:"
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
   16 format(1x,i4,1x,2a4,2x,100(3x,f8.5," of",i5))
      !
      call lsj_deallocate_asf_basis_LS(asf_set_LS)
      if (debuging /= 0 )then
         close(57)
      end if
   end subroutine lsj_control_transformation
   !
   !
   subroutine lsj_deallocate_asf_basis_LS(asf_set_LS)
   !--------------------------------------------------------------------
   ! Dellocates the storage of asf_set_LS.
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      type(asf_basis_LS), intent(inout)  :: asf_set_LS
      !
      integer :: i
      !
      deallocate(asf_set_LS%csf_set_LS%shell)
      deallocate(asf_set_LS%csf_set_LS%parent)
      do i = 1, asf_set_LS%csf_set_LS%nocsf, 1
         deallocate(asf_set_LS%csf_set_LS%csf(i)%occupation)
         deallocate(asf_set_LS%csf_set_LS%csf(i)%seniority)
         deallocate(asf_set_LS%csf_set_LS%csf(i)%shellL)
         deallocate(asf_set_LS%csf_set_LS%csf(i)%shellS)
         deallocate(asf_set_LS%csf_set_LS%csf(i)%shellLX)
         deallocate(asf_set_LS%csf_set_LS%csf(i)%shellSX)
      end do
      deallocate(asf_set_LS%csf_set_LS%csf)
      do i = 1, asf_set_jj%noasf
         deallocate(asf_set_LS%asf(i)%eigenvector)
      end do
      deallocate(asf_set_LS%asf)
      !
   end subroutine lsj_deallocate_asf_basis_LS
   !
   !
   subroutine lsj_form_csf_basis_LS()
   !--------------------------------------------------------------------
   ! This subroutine fills up the variable asf_set_LS%csf_set_LS
   ! with data generated using the one from asf_set_jj%csf_set.
   ! 
   ! This subroutine contains the following internal subroutines:
   !   * subroutine           lsj_form_csf_basis_LS_action 
   !   * subroutine           lsj_form_csf_basis_LS_add_qn 
   !   * recursive subroutine lsj_form_csf_basis_LS_job_count
   !   * function             lsj_form_csf_basis_LS_equiv
   !
   ! Calls: angular_momentum_l, lsj_form_csf_basis_LS_equiv.
   !--------------------------------------------------------------------
      !
      implicit none
      !
      type(nl), dimension(:), pointer :: shell_temp
      type(lsj_list)                 :: nonequiv_csfs_jj
      !
      integer :: isubc, isubc2, icsf_jj, icsf_jj2, icsf_jj_real, icsf_LS, N
      integer :: action_type 
      logical :: new_one, found_parent_minus, found_parent_plius
      !
      integer, dimension(:), pointer           :: all_occupation
      integer(kind=i1b), dimension(:), pointer :: Li, Si, L_i, S_i, w, Q
      integer(kind=i1b)                        :: J
      !
      asf_set_LS%csf_set_LS%number_of_electrons =                          &
                                      asf_set_jj%csf_set%number_of_electrons
      !
      !1. define  nl, parent
      allocate(shell_temp(asf_set_jj%csf_set%nwshells))
      !
      asf_set_LS%csf_set_LS%nwshells = 0
      do isubc = 1, asf_set_jj%csf_set%nwshells, 1
         new_one = .true.
         do isubc2 = 1, isubc-1, 1
            if(asf_set_jj%csf_set%subshell(isubc)%n == 			       &
               asf_set_jj%csf_set%subshell(isubc2)%n  .and.                    &
               angular_momentum_l(asf_set_jj%csf_set%subshell(isubc)%kappa) == &
               angular_momentum_l(asf_set_jj%csf_set%subshell(isubc2)%kappa))  & 
			   new_one = .false.
	 end do
	 if(new_one) then
            asf_set_LS%csf_set_LS%nwshells = asf_set_LS%csf_set_LS%nwshells + 1
            shell_temp(asf_set_LS%csf_set_LS%nwshells)%n =                    &
                                           asf_set_jj%csf_set%subshell(isubc)%n
            shell_temp(asf_set_LS%csf_set_LS%nwshells)%l =                    &
                   angular_momentum_l(asf_set_jj%csf_set%subshell(isubc)%kappa)
         end if
      end do
      !
      allocate(asf_set_LS%csf_set_LS%shell(asf_set_LS%csf_set_LS%nwshells))
      allocate(asf_set_LS%csf_set_LS%parent(asf_set_LS%csf_set_LS%nwshells))
      !
      do isubc=1, asf_set_LS%csf_set_LS%nwshells, 1
         asf_set_LS%csf_set_LS%shell(isubc) = shell_temp(isubc)
      end do
      !
      deallocate(shell_temp)
      !
      ! begin find parent
      do isubc = 1, asf_set_LS%csf_set_LS%nwshells, 1
         found_parent_minus = .false.
	 found_parent_plius = .false.
	 do isubc2 = 1, asf_set_jj%csf_set%nwshells , 1
            if(asf_set_LS%csf_set_LS%shell(isubc)%l ==                       & 
               angular_momentum_l(asf_set_jj%csf_set%subshell(isubc2)%kappa) &
               .and.                                                         &
               asf_set_LS%csf_set_LS%shell(isubc)%n ==                       &
                                   asf_set_jj%csf_set%subshell(isubc2)%n) then 
               if(asf_set_jj%csf_set%subshell(isubc2)%kappa > 0) then
		  found_parent_minus = .true.
		  asf_set_LS%csf_set_LS%parent(isubc)%parent_minus = isubc2
               else
		  found_parent_plius = .true.
		  asf_set_LS%csf_set_LS%parent(isubc)%parent_plius = isubc2
               end if
            end if
	 end do
         if(.not.found_parent_plius)  & 
            asf_set_LS%csf_set_LS%parent(isubc)%parent_plius = 0
         if(.not.found_parent_minus)  &
            asf_set_LS%csf_set_LS%parent(isubc)%parent_minus = 0
      end do
      !
      ! end find parent  --------------------
      !
      !2. define the number of "core" shells
      !   (the LS shell is supposed to be "core" if: 
      !        1. l=0 and corresponding jj subshell is "core"
      !        2. l<>0 l+ and l - "core" subshells            )
      asf_set_LS%csf_set_LS%nwcore=0
      do isubc=1, asf_set_LS%csf_set_LS%nwshells, 1
         if(asf_set_LS%csf_set_LS%parent(isubc)%parent_minus .le. &
            asf_set_jj%csf_set%nwcore 		                  &
            .and.                                                 &
            asf_set_LS%csf_set_LS%parent(isubc)%parent_plius .le. &
            asf_set_jj%csf_set%nwcore) 			          &
            asf_set_LS%csf_set_LS%nwcore =                        &
            asf_set_LS%csf_set_LS%nwcore + 1
      end do
      !
      !3. form the list of "nonequivalent" csfs_jj 
      !   (i.e. csfs_jj different in J,parity, or
      !   some l's occupation numbers Ni = N_(i+) + N_(i-))
      allocate(nonequiv_csfs_jj%items(asf_set_jj%csf_set%nocsf))
      nonequiv_csfs_jj%list_size = 0
      do icsf_jj = 1, asf_set_jj%csf_set%nocsf, 1
         new_one = .true.
         do icsf_jj2 = 1 , icsf_jj - 1 !asf_set_jj%csf_set%nocsf
            if(lsj_form_csf_basis_LS_equiv(icsf_jj, icsf_jj2)) then
               new_one = .false.
               exit
            end if
	 end do
	 if(new_one) then
	    nonequiv_csfs_jj%list_size = nonequiv_csfs_jj%list_size + 1
            nonequiv_csfs_jj%items(nonequiv_csfs_jj%list_size) = icsf_jj
         end if
      end do
      !
      !4. for each nonequivalent csf_jj find all the csfs_LS 
      !   
      !	  To avoid the dependency on the number of subshells 
      !   the recursive subroutine is used 
      allocate(Li(asf_set_LS%csf_set_LS%nwshells))
      allocate(L_i(asf_set_LS%csf_set_LS%nwshells))
      allocate(Si(asf_set_LS%csf_set_LS%nwshells))
      allocate(S_i(asf_set_LS%csf_set_LS%nwshells))
      allocate(Q(asf_set_LS%csf_set_LS%nwshells))
      allocate(w(asf_set_LS%csf_set_LS%nwshells))
      allocate(all_occupation(asf_set_LS%csf_set_LS%nwshells))
      !
      !4.1 - find the number of csfs_LS
      asf_set_LS%csf_set_LS%nocsf = 0
      !
      do icsf_jj=1, nonequiv_csfs_jj%list_size, 1
         !
	 ! set the variable for conviniency ...
	 icsf_jj_real = nonequiv_csfs_jj%items(icsf_jj)
         J = asf_set_jj%csf_set%csf(icsf_jj_real)%totalJ
	 !
	 !define the occupation numbers
	 !
	 do isubc = 1, asf_set_LS%csf_set_LS%nwshells, 1
	    all_occupation(isubc) = 0
            do isubc2 = 1, asf_set_jj%csf_set%nwshells, 1
               if(asf_set_LS%csf_set_LS%shell(isubc)%l ==                      & 
                  angular_momentum_l(asf_set_jj%csf_set%subshell(isubc2)%kappa)&
                 .and.	                                                       &
                 asf_set_LS%csf_set_LS%shell(isubc)%n ==                       &
                 asf_set_jj%csf_set%subshell(isubc2)%n)                        &
                 all_occupation(isubc) = all_occupation(isubc) +               &
                 asf_set_jj%csf_set%csf(icsf_jj_real)%occupation(isubc2)
            end do  !isubc2
            !
         end do  !isubc
         N = 0;    action_type = 1
         call lsj_form_csf_basis_LS_job_count(1, N)
         asf_set_LS%csf_set_LS%nocsf=asf_set_LS%csf_set_LS%nocsf + N
      end do
      !
      ! 4.2 - fill up the csf_LS arrays with the corresponding quantum numbers
      ! 4.2.1  allocate the arrays
      allocate(asf_set_LS%csf_set_LS%csf(asf_set_LS%csf_set_LS%nocsf))
      do N = 1, asf_set_LS%csf_set_LS%nocsf, 1
         allocate(asf_set_LS%csf_set_LS%csf(N)%occupation(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%seniority(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%w(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%shellL(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%shellS(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%shellLX(asf_set_LS%csf_set_LS%nwshells))
         allocate(asf_set_LS%csf_set_LS%csf(N)%shellSX(asf_set_LS%csf_set_LS%nwshells))
      end do
      !
      ! 4.2.2 - fill them with quantum numbers
      icsf_LS = 0
      !
      do icsf_jj=1, nonequiv_csfs_jj%list_size, 1
         !
         ! set the variable for conviniency ...
         icsf_jj_real = nonequiv_csfs_jj%items(icsf_jj)
         J = asf_set_jj%csf_set%csf(icsf_jj_real)%totalJ
         do isubc = 1, asf_set_LS%csf_set_LS%nwshells, 1
            all_occupation(isubc)=0
            do isubc2 = 1, asf_set_jj%csf_set%nwshells, 1
               if(asf_set_LS%csf_set_LS%shell(isubc)%l ==                      & 
                  angular_momentum_l(asf_set_jj%csf_set%subshell(isubc2)%kappa)&
                  .and.						               &
                  asf_set_LS%csf_set_LS%shell(isubc)%n ==                      &
                  asf_set_jj%csf_set%subshell(isubc2)%n)                       &
                  all_occupation(isubc) = all_occupation(isubc) +              &
	          asf_set_jj%csf_set%csf(icsf_jj_real)%occupation(isubc2)
            end do  !isubc2
         end do  !isubc
         action_type = 2
         call lsj_form_csf_basis_LS_job_count(1, N)
      end do
      !
      deallocate(nonequiv_csfs_jj%items)
      !
      deallocate(Q);   deallocate(w);    deallocate(all_occupation)
      deallocate(Li);  deallocate(L_i);  deallocate(Si);
      deallocate(S_i)
      !
   contains
      !
      subroutine lsj_form_csf_basis_LS_action(tip, irez)
      !--------------------------------------------------------------------
      ! The subroutine defines the "action" of subroutine 
      ! lsj_form_csf_basis_LS_job_count: whether it counts the number 
      ! of csfs_LS (asf_set_LS%csf_set_LS%novcsf) or fills the
      ! arrays of wave functions in LS coupling with
      ! asf_set_LS%csf_set_LS%csf(...) with the corresponding
      ! quantum nubers.
      !
      ! Calls: lsj_form_csf_basis_LS_add_qn.
      !--------------------------------------------------------------------
         !
         integer tip, irez
         !
         if (tip .eq. 1) then
            irez = irez + 1
         else
            call lsj_form_csf_basis_LS_add_qn()
         end if
         !
      end subroutine lsj_form_csf_basis_LS_action
      !
      !
      subroutine lsj_form_csf_basis_LS_add_qn()
      !--------------------------------------------------------------------
      ! The subroutine adds quantum numbers stored in temprorary arrays
      ! Li, Si, L_i, S_i, w, Q to the corresponding arrays
      ! of asf_set_LS%csf_set_LS%csf() (i.e. to the arrays of
      ! corresponding quantum numbers of the wave function in LS
      ! coupling).
      !
      ! Calls: 
      !--------------------------------------------------------------------
         !
         integer :: isubcx
         !
         icsf_LS = icsf_LS + 1
         ! 
         if(icsf_LS  .gt.  asf_set_LS%csf_set_LS%nocsf) then
            stop 'lsj_form_csf_basis_LS_add_qn(): program stop A.'
         end if
         !
         asf_set_LS%csf_set_LS%csf(icsf_LS)%totalJ = J               
         !
         asf_set_LS%csf_set_LS%csf(icsf_LS)%parity = & 
            asf_set_jj%csf_set%csf(icsf_jj_real)%parity
         !
         do isubcx = 1, asf_set_LS%csf_set_LS%nwshells, 1
            asf_set_LS%csf_set_LS%csf(icsf_LS)%occupation(isubcx) =        &
                                               all_occupation(isubcx)
            asf_set_LS%csf_set_LS%csf(icsf_LS)%shellL(isubcx)  = Li(isubcx)
            asf_set_LS%csf_set_LS%csf(icsf_LS)%shellS(isubcx)  = Si(isubcx)
            asf_set_LS%csf_set_LS%csf(icsf_LS)%shellLX(isubcx) = L_i(isubcx)  
            asf_set_LS%csf_set_LS%csf(icsf_LS)%shellSX(isubcx) = S_i(isubcx)
            asf_set_LS%csf_set_LS%csf(icsf_LS)%w(isubcx)       = w(isubcx)
            asf_set_LS%csf_set_LS%csf(icsf_LS)%seniority(isubcx) =         &
                     2*asf_set_LS%csf_set_LS%shell(isubcx)%l + 1 - Q(isubcx)
         end do
         !
      end subroutine lsj_form_csf_basis_LS_add_qn 
      !
      !
      recursive subroutine lsj_form_csf_basis_LS_job_count(isubc, rez)
      !--------------------------------------------------------------------
      ! Recursive subroutine for the calculation of the
      ! number of csfs_LS and corresponding quantum numbers.
      !
      ! Calls: lsj_form_csf_basis_LS_action,lsj_form_csf_basis_LS_job_count,
      !        lsj_get_subshell_term_LS, lsj_triangle.
      !--------------------------------------------------------------------
         integer :: isubc, rez
         integer :: iterm, nr_terms, suma
         integer                                :: l_shell, N
         type(subshell_term_LS), dimension(120) :: LS_terms
         integer                                :: number
         integer(kind=i1b) :: tempLmax,tempLmin,tempSmax,tempSmin,tempS,tempL
         !
         if(isubc.gt.(asf_set_LS%csf_set_LS%nwshells) .or. isubc.lt.1) then
            print *, 'isubc = ', isubc
            stop "lsj_form_csf_basis_LS_job_count(): program stop A."
         end if
         !
         if(isubc.le.asf_set_LS%csf_set_LS%nwshells) then
            if(all_occupation(isubc).eq.0) then
               if(isubc.gt.1) then
                  Li(isubc)  = 0;          L_i(isubc) = L_i(isubc-1)
                  Si(isubc)  = 0;          S_i(isubc) = S_i(isubc-1)
                  if(isubc .lt. asf_set_LS%csf_set_LS%nwshells) then
                     call lsj_form_csf_basis_LS_job_count(isubc + 1, rez)
                  else
                     if(lsj_triangle(S_i(isubc), L_i(isubc), J))  &
	   	          call lsj_form_csf_basis_LS_action(action_type, rez) 
			  !rez=rez+1
                  end if 
               else
                  Li(isubc)  = 0;          L_i(isubc) = 0
                  Si(isubc)  = 0;          S_i(isubc) = 0
               end if
            else
               N = all_occupation(isubc)
               l_shell = asf_set_LS%csf_set_LS%shell(isubc)%l
               call lsj_get_subshell_term_LS(l_shell,N,LS_terms,number)
               do iterm=1, number, 1
                  Li(isubc)  = LS_terms(iterm)%LL
                  Si(isubc)  = LS_terms(iterm)%S
                  w(isubc)   = LS_terms(iterm)%w
                  Q(isubc)   = LS_terms(iterm)%Q
                  if(isubc.eq.1) then
                     L_i(isubc) = LS_terms(iterm)%LL
                     S_i(isubc) = LS_terms(iterm)%S
                     if(asf_set_LS%csf_set_LS%nwshells.gt.1) then
                        call lsj_form_csf_basis_LS_job_count(isubc + 1, rez)
                     else 
                        if(lsj_triangle(S_i(isubc), L_i(isubc), J)) &
	   	            call lsj_form_csf_basis_LS_action(action_type, rez) !rez=rez+1
                     end if
                  else
                     tempLmax=L_i(isubc-1)+Li(isubc)
                     tempLmin=abs(L_i(isubc-1)-Li(isubc))
                     tempSmax=S_i(isubc-1)+Si(isubc)
                     tempSmin=abs(S_i(isubc-1)-Si(isubc))
                     do tempL = tempLmin, tempLmax, 2
                        L_i(isubc) = tempL
                        do tempS = tempSmin, tempSmax, 2
                           S_i(isubc) = tempS
                           if(isubc.lt.asf_set_LS%csf_set_LS%nwshells) then
                              call lsj_form_csf_basis_LS_job_count(isubc+1,rez)
                           else
                              if(lsj_triangle(S_i(isubc), L_i(isubc), J)) &
                                  call lsj_form_csf_basis_LS_action       &
				               (action_type, rez) ! rez=rez+1
                           end if
                        end do   ! tempS
                     end do      ! tempL
                  end if
               end do            ! iterm
            end if
         end if
         ! 
      end subroutine lsj_form_csf_basis_LS_job_count
      !
      !
      function lsj_form_csf_basis_LS_equiv(ncsf1, ncsf2)  result(rez)
      !--------------------------------------------------------------------
      ! This subroutine defines the "equivalency" of two csfs_jj
      ! in the sence of generation of csfs_LS
      ! number of csfs_LS and corresponding quantum numbers.
      !
      ! Calls: angular_momentum_l.
      !--------------------------------------------------------------------
         !
         integer :: ncsf1, ncsf2, isubc_LS, isubc_jj, NLS1, NLS2
         logical rez
         !
         rez= .true.
         !
         if(asf_set_jj%csf_set%csf(ncsf1)%totalJ .ne.                          &
                                       asf_set_jj%csf_set%csf(ncsf2)%totalJ) then
            rez = .false.
         else if (asf_set_jj%csf_set%csf(ncsf1)%parity .ne.                    &
		                       asf_set_jj%csf_set%csf(ncsf2)%parity) then
            rez = .false.
         else
            do isubc_LS = asf_set_LS%csf_set_LS%nwcore + 1,                    &
	                                 asf_set_LS%csf_set_LS%nwshells, 1
               NLS1 = 0;   NLS2	= 0
               do isubc_jj = 1, asf_set_jj%csf_set%nwshells, 1
                  if(asf_set_LS%csf_set_LS%shell(isubc_LS)%l ==                &
                   angular_momentum_l(asf_set_jj%csf_set%subshell(isubc_jj)%kappa)&
                   .and. asf_set_LS%csf_set_LS%shell(isubc_LS)%n ==            &
			              asf_set_jj%csf_set%subshell(isubc_jj)%n) then
   		      NLS1 = NLS1+asf_set_jj%csf_set%csf(ncsf1)%occupation(isubc_jj)
                      NLS2 = NLS2+asf_set_jj%csf_set%csf(ncsf2)%occupation(isubc_jj)
	           end if
               end do
               if(NLS1.ne.NLS2) rez=.false.
	    end do
         end if
         !
      end function lsj_form_csf_basis_LS_equiv
      !
   end subroutine lsj_form_csf_basis_LS
   !
   !
   subroutine lsj_get_subshell_term_LS(l_shell,N,LS,number)
   !--------------------------------------------------------------------
   !  This procedure return all allowed subshell terms 
   !                                                  (l, w, Q, L, S)
   !  for given l^N which must be 0, 1, 2 or 3. 
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      integer, intent(in)                                 :: l_shell, N
      type(subshell_term_LS), dimension(120), intent(out) :: LS
      integer, intent(out)                                :: number
      !
      integer  :: M_Q, i, j
      integer, save :: k = 0
      !
      M_Q = N - 2* l_shell - 1;  j = 0
      select case (l_shell)
      case (0)
         do i = 1,2
            if (mod(M_Q + term_LS_s(i)%Q,2) == 0) then
	       if (abs(M_Q) <= term_LS_s(i)%Q) then
	          j = j + 1
	          LS(j)%l_shell  = term_LS_s(i)%l_shell
                  LS(j)%w	 = term_LS_s(i)%w 
	          LS(j)%Q	 = term_LS_s(i)%Q
                  LS(j)%LL	 = term_LS_s(i)%LL
                  LS(j)%S	 = term_LS_s(i)%S
	       end if
	    end if
         end do
      case (1)
         do i = 1,6
            if (mod(M_Q + term_LS_p(i)%Q,2) == 0) then
	       if (abs(M_Q) <= term_LS_p(i)%Q) then
	          j = j + 1
	          LS(j)%l_shell  = term_LS_p(i)%l_shell
                  LS(j)%w	 = term_LS_p(i)%w 
	          LS(j)%Q	 = term_LS_p(i)%Q
                  LS(j)%LL	 = term_LS_p(i)%LL
                  LS(j)%S	 = term_LS_p(i)%S
	       end if
	    end if
         end do
      case (2)
         do i = 1,32
            if (mod(M_Q + term_LS_d(i)%Q,2) == 0) then
	       if (abs(M_Q) <= term_LS_d(i)%Q) then
	          j = j + 1
	          LS(j)%l_shell  = term_LS_d(i)%l_shell
                  LS(j)%w	 = term_LS_d(i)%w 
	          LS(j)%Q	 = term_LS_d(i)%Q
                  LS(j)%LL	 = term_LS_d(i)%LL
                  LS(j)%S	 = term_LS_d(i)%S
	       end if
	    end if
         end do
      case (3)
         do i = 1,238
            if (mod(M_Q + term_LS_f(i)%Q,2) == 0) then
	       if (abs(M_Q) <= term_LS_f(i)%Q) then
	          j = j + 1
	          LS(j)%l_shell  = term_LS_f(i)%l_shell
                  LS(j)%w	 = term_LS_f(i)%w 
	          LS(j)%Q	 = term_LS_f(i)%Q
                  LS(j)%LL	 = term_LS_f(i)%LL
                  LS(j)%S	 = term_LS_f(i)%S
	       end if
	    end if
         end do
      case default
         k = k + 1
         if (k < 6) then
            !!print *, "l_shell = ",l_shell
            print *, "lsj_get_subshell_term_LS(): program stop A."
	 end if
	 j = 0
         !! stop  "lsj_get_subshell_term_LS(): program stop A."
      end select
      !
      number = j
      !
   end subroutine lsj_get_subshell_term_LS
   !
   !
   subroutine lsj_interprete_levels(record,levels,number_of_levels,fail)
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
   ! Calls: get_integer_from_string.
   !--------------------------------------------------------------------
      !
      character(len=*), intent(in)               :: record
      logical, intent(out)                       :: fail
      integer, intent(out)                       :: number_of_levels
      integer, dimension(1:10000), intent(inout) :: levels
      !
      logical, dimension(200)              :: low_to
      character(len=500)                   :: string
      integer                              :: a, i, lower, n
      integer, dimension(200)              :: low
      integer(kind=i1b), dimension(10000)  :: run
      !
      fail = .true.;  levels(1:10000) = 0;  number_of_levels = 0
      run(1:10000) = 0
      !
      string = adjustl(record);
      !
      n = 0;   lower = 0;  low_to(:) = .false.
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
   end subroutine lsj_interprete_levels
   !
   !
   subroutine lsj_print_conf_scheme_LS(stream,csf_set_LS)
   !--------------------------------------------------------------------
   ! Print all information about the CSF scheme csf_set in a nead
   ! format on stream.
   !
   ! Calls: lsj_print_single_config_LS.
   !--------------------------------------------------------------------
      !
      integer, intent(in)            :: stream
      type(csf_basis_LS), intent(in) :: csf_set_LS
      !
      integer :: i, j
      !
      write(stream,*) " "
      write(stream,*) "The current configuration scheme with",         &
                      csf_set_LS%nocsf,"CSF in LS coupling" 
      write(stream,*)                                                  &
           "                                     is defined as follows:"
      write(stream,*) " "
      write(stream,*) "Number of (nonrelativistic) subshells:",        &
                                          csf_set_LS%nwshells
      write(stream,*) "Core shells:"
      write(stream,1)(orbital_name(csf_set_LS%shell(j)%n,              &
                      -1-csf_set_LS%shell(j)%l),j=1,csf_set_LS%nwcore)
      write(stream,*) "Peel shells:"
      write(stream,1)(orbital_name(csf_set_LS%shell(j)%n,              &
                      -1-csf_set_LS%shell(j)%l),                       &
		          j=csf_set_LS%nwcore+1,csf_set_LS%nwshells)
      write(stream,*) "CSF(s):"
      do i = 1,csf_set_LS%nocsf
         call lsj_print_single_config_LS(stream,csf_set_LS,i)
      end do
    1 format(1x,100(a4,3x))
      !
   end subroutine lsj_print_conf_scheme_LS
   !
   !
   subroutine lsj_print_single_config_jj(stream,csf_set,csf_number)
   !--------------------------------------------------------------------
   ! Print all information about the single CSF scheme csf_set in a nead
   ! format on stream.
   !
   ! Calls: angular_momentum_string, orbital_name.
   !--------------------------------------------------------------------
      !
      integer, intent(in)         :: stream, csf_number
      type(csf_basis), intent(in) :: csf_set
      !
      integer                                  :: j, counter
      integer, dimension(:), pointer           :: occupation
      character(len=4), dimension(:), pointer  :: string_J, string_X
      !
      counter = 0
      allocate(occupation(asf_set_jj%csf_set%nwshells))
      allocate(string_J(asf_set_jj%csf_set%nwshells))
      allocate(string_X(asf_set_jj%csf_set%nwshells))
      do j = csf_set%nwcore+1, csf_set%nwshells
         if (csf_set%csf(csf_number)%occupation(j) > 0) then
            counter = counter +1
            occupation(counter) = j
            if (csf_set%csf(csf_number)%subshellJ(j) == 0) then
               string_J(counter) = "    "
            else
               string_J(counter) = angular_momentum_string     &
	       (1*csf_set%csf(csf_number)%subshellJ(j),4)
            end if
            if (csf_set%csf(csf_number)%subshellX(j) == 0) then
               string_X(counter) = "     "
            else
               string_X(counter) = angular_momentum_string     &
	       (1*csf_set%csf(csf_number)%subshellX(j),4)
            end if
         end if
      end do
      if (stream == -1) then
         write(*,1)csf_number,                                           &
	    (orbital_name(csf_set%subshell(occupation(j))%n,             &
            csf_set%subshell(occupation(j))%kappa),                      &
            csf_set%csf(csf_number)%occupation(occupation(j)),j=1,counter)
         write(*,2)(string_J(j),j=1,counter)
         write(*,3)(string_X(j),j=2,counter-1),angular_momentum_string   &
         (1*csf_set%csf(csf_number)%totalJ,4)
      end if
      deallocate(occupation)
      deallocate(string_J)
      deallocate(string_X)
    1 format(4x,i9,')',100(a4,'(',i2,')'2x))
    2 format(18x,a4,100(6X,a4))
    3 format(31x,a4,100(6X,a5))
      !
   end subroutine lsj_print_single_config_jj
   !
   !
   subroutine lsj_print_single_config_LS(stream,csf_set_LS,csf_number)
   !--------------------------------------------------------------------
   ! Print all information about the single CSF scheme csf_set in a nead
   ! format on stream.
   !
   ! Calls: angular_momentum_string, lsj_spectroscopic_LS, orbital_name.
   !--------------------------------------------------------------------
      !
      integer, intent(in)            :: stream, csf_number
      type(csf_basis_LS), intent(in) :: csf_set_LS
      !
      integer                                 :: j, counter
      integer, dimension(:), pointer          :: occupation
      character(len=4)                        :: LS, XLS
      character(len=4), dimension(:),pointer  :: string_LS, string_XLS
      !
      counter = 0
      allocate(occupation(csf_set_LS%nwshells))
      allocate(string_LS(csf_set_LS%nwshells))
      allocate(string_XLS(csf_set_LS%nwshells))
      do j = csf_set_LS%nwcore+1, csf_set_LS%nwshells     
         if (csf_set_LS%csf(csf_number)%occupation(j) > 0) then
            counter = counter +1;        occupation(counter) = j
            call lsj_spectroscopic_LS(csf_number,j,LS,XLS)
            string_LS(counter)  = LS;    string_XLS(counter) = XLS
         end if
      end do
      if (stream == -1) then
         write(*,1)csf_number,                                            &
            (orbital_name(csf_set_LS%shell(occupation(j))%n,              &
            -1-csf_set_LS%shell(occupation(j))%l),                        &
            csf_set_LS%csf(csf_number)%occupation(occupation(j)),j=1,counter)
         write(*,2)(string_LS(j),j=1,counter),(string_XLS(j),j=2,counter),&
             angular_momentum_string(1*csf_set_LS%csf(csf_number)%totalJ,4)
      else
         write(stream,1)csf_number,                                       &
	    (orbital_name(csf_set_LS%shell(occupation(j))%n,              &
            -1-csf_set_LS%shell(occupation(j))%l),                        &
            csf_set_LS%csf(csf_number)%occupation(occupation(j)),j=1,counter)
         write(stream,2)(string_LS(j),j=1,counter),(string_XLS(j),j=2,counter),&
             angular_momentum_string(1*csf_set_LS%csf(csf_number)%totalJ,4)
      end if
      deallocate(occupation)
      deallocate(string_LS)
      deallocate(string_XLS)
    1 format(i10,')',4x,100(a3,'(',i2,')'2x))
    2 format(19x,a4,100(5X,a4))
      !
   end subroutine lsj_print_single_config_LS
   !
   !
   subroutine lsj_spectroscopic_LS(csf_number,shell_number,string_LS,  &
                                                             string_XLS)
   !--------------------------------------------------------------------
   ! A spectroscopic notation of shell in LS coupling is return.
   !
   ! Calls: angular_momentum_string.
   !--------------------------------------------------------------------
      !
      integer, intent(in)            :: csf_number,shell_number
      character(len=4), intent(out)  :: string_LS, string_XLS
      !
      character(len=1)               :: string_S, string_v
      !
      string_S = angular_momentum_string                                      &
           (2*(1+asf_set_LS%csf_set_LS%csf(csf_number)%shellS(shell_number)),4)
      if (asf_set_LS%csf_set_LS%shell(shell_number)%l < 3)    then
         string_v = angular_momentum_string                                   &
            (2*asf_set_LS%csf_set_LS%csf(csf_number)%seniority(shell_number),1)
      else
!      else if                                                                &
!       (asf_set_LS%csf_set_LS%csf(csf_number)%occupation(shell_number)<4) then
!         string_v = angular_momentum_string                                   &
!            (2*asf_set_LS%csf_set_LS%csf(csf_number)%seniority(shell_number),1)
!      else
         string_v = angular_momentum_string                                   &
            (2*asf_set_LS%csf_set_LS%csf(csf_number)%w(shell_number),1)
      end if
      string_LS = string_S //                     &
        L_string(asf_set_LS%csf_set_LS%csf(csf_number)%shellL(shell_number)/2)&
        //string_v
      !
      string_S = angular_momentum_string                                      &
         (2*(1+asf_set_LS%csf_set_LS%csf(csf_number)%shellSX(shell_number)),1)
      string_XLS = string_S //                                                &
         L_string(asf_set_LS%csf_set_LS%csf(csf_number)%shellLX(shell_number)/2)
      !
   end subroutine lsj_spectroscopic_LS
   !
   !
   subroutine lsj_transformation_ASF(level_jj)
   !--------------------------------------------------------------------
   ! Expands an atomic state function, which is represented 
   ! in a jj-coupling CSF basis into a basis of LS-coupling CSF.
   !
   ! Calls: lsj_transformation_LS_jj_gen.
   !--------------------------------------------------------------------
      !
      integer, intent(in) :: level_jj
      !
      real(kind=dp)       :: wa, wa_transformation
      !
      integer  :: jj_number, LS_number, i
      integer  :: shell_number
      !
!GG debug begin (4 lines)
      do LS_number = 1, asf_set_LS%csf_set_LS%nocsf                 !real run
!            do jj_number = 1, asf_set_jj%csf_set%nocsf             !for debug
         if (asf_set_LS%csf_set_LS%csf(LS_number)%parity ==  &      !real run
             asf_set_jj%asf(level_jj)%parity) then                  !real run
! debug end
            wa = zero
!GG debug begin (4 lines)
!            do LS_number = 1, asf_set_LS%csf_set_LS%nocsf          !for debug
!              if (asf_set_LS%csf_set_LS%csf(LS_number)%parity == & !for debug
!	          asf_set_jj%asf(level_jj)%parity) then             !for debug
            do jj_number = 1, asf_set_jj%csf_set%nocsf              !real run
! debug end
               if ((asf_set_LS%csf_set_LS%csf(LS_number)%totalJ ==      &
!GG debug begin (5 lines)
!	             asf_set_jj%csf_set%csf(jj_number)%totalJ)          &!for debug
	             asf_set_jj%csf_set%csf(jj_number)%totalJ)  .and.   &!real run
		     abs(asf_set_jj%asf(level_jj)%eigenvector(jj_number))>&!real run
                                                        eps20)  then!real run
!		                                              )  then!for debug
! debug end
                  wa_transformation =                                    &
		        lsj_transformation_LS_jj_gen(jj_number,LS_number)
                  if (debuging /= 0) then
                     if(abs(wa_transformation) > eps20) then
                      write(57,*)"wa_transformation", "jj_number", "LS_number"
                      write(57,*) wa_transformation, jj_number, LS_number
                      write(57,*)"asf_set_jj%asf(level_jj)%eigenvector(jj_number)"
                      write(57,*) asf_set_jj%asf(level_jj)%eigenvector(jj_number)
                      write(57,*)" "
                      write(57,*)" "
                     end if
                  end if
                  wa = wa +                                       &
!GG debug begin (4 lines)
                  asf_set_jj%asf(level_jj)%eigenvector(jj_number) & !real run
                  * wa_transformation                               !real run
!	          wa_transformation*wa_transformation               !for debug
! debug end
               end if
!GG debug begin (1 lines)
!              end if           !for debug
! debug end
            end do
            asf_set_LS%asf(level_jj)%eigenvector(LS_number) = wa
            asf_set_LS%asf(level_jj)%level_No = asf_set_jj%asf(level_jj)%level_No
!GG debug begin (9 lines)
!            if(abs(wa) > eps10) then                               !for debug
!               if(abs(abs(wa) - one) > eps10) then                 !for debug
!                    print*, " Wa*wa=",wa, "jj_number=",jj_number   !for debug
!                    print*, " Wa*wa=",wa, "LS_number=",LS_number   !for debug
!               end if                                              !for debug
!            end if                                                 !for debug
!             if(abs(wa_transformation) > eps10)  &        !for debug
!             print*, "wa_tra=",wa_transformation,"LS_number=",LS_number !for debug
         end if                                                    !real run
!GG debug end
      end do
      !
   end subroutine lsj_transformation_ASF
   !
   !
   function lsj_transformation_LS_jj_gen(jj_number,LS_number)  result(wa)
   !--------------------------------------------------------------------
   !  This procedure return the value of the transformation matrix 
   !  from jj- to LS-coupling scheme in the case of any number of
   !  open shells.
   !
   !  Calls: lsj_coefficient_LS_jj, LS_jj_transformation_insade.
   !--------------------------------------------------------------------
      !
      integer, intent(in)   :: jj_number,LS_number
      !
      real(kind=dp)         :: wa
      !
      integer   :: total_number, shell_number, number_minus, number_plius, &
                   jj_minus, N_minus, Q_minus, J_1_minus,                  &
	 	   jj_plius, N_plius, Q_plius, J_1_plius, J_1,             &
                   l_shell, N_LS, W_1, Q_1, L_1, S_1
      !
      shell_number = asf_set_LS%csf_set_LS%nwcore+1
      number_minus = asf_set_LS%csf_set_LS%parent(shell_number)%parent_minus
      number_plius = asf_set_LS%csf_set_LS%parent(shell_number)%parent_plius
      if (number_minus+1 /= number_plius .and.  &
          number_minus*number_plius  /= 0) then
        stop "lsj_transformation_LS_jj_gen(): program stop A."
      end if
      !
      if (number_minus == 0) then
         jj_minus = iabs(asf_set_jj%csf_set%subshell(number_plius)%kappa)*2 - 3
         N_minus  = 0;    J_1_minus = 0;   Q_minus = (jj_minus + 1)/2
      else
         jj_minus = iabs(asf_set_jj%csf_set%subshell(number_minus)%kappa)*2 - 1
         N_minus  = asf_set_jj%csf_set%csf(jj_number)%occupation(number_minus)
         Q_minus  = (jj_minus +1)/2 -                                 &
                 asf_set_jj%csf_set%csf(jj_number)%seniority(number_minus)
         J_1_minus = asf_set_jj%csf_set%csf(jj_number)%subshellJ(number_minus)
         J_1 = asf_set_jj%csf_set%csf(jj_number)%subshellX(number_minus)
      end if
      !
      if (number_plius == 0) then
         jj_plius = iabs(asf_set_jj%csf_set%subshell(number_minus)%kappa)*2 + 1
         N_plius  = 0;    J_1_plius = 0;   Q_plius = (jj_plius + 1)/2
      else
         jj_plius = iabs(asf_set_jj%csf_set%subshell(number_plius)%kappa)*2 - 1
         N_plius  =  asf_set_jj%csf_set%csf(jj_number)%occupation(number_plius)
         Q_plius  = (jj_plius + 1)/2 -                                 &
                 asf_set_jj%csf_set%csf(jj_number)%seniority(number_plius)
         J_1_plius = asf_set_jj%csf_set%csf(jj_number)%subshellJ(number_plius)
         J_1 = asf_set_jj%csf_set%csf(jj_number)%subshellX(number_plius)
      end if
      !
      l_shell = asf_set_LS%csf_set_LS%shell(shell_number)%l
      N_LS    = asf_set_LS%csf_set_LS%csf(LS_number)%occupation(shell_number)
      W_1     = asf_set_LS%csf_set_LS%csf(LS_number)%w(shell_number)
      Q_1     = 2*l_shell+1-                                               & 
                asf_set_LS%csf_set_LS%csf(LS_number)%seniority(shell_number)
      !
      L_1     = asf_set_LS%csf_set_LS%csf(LS_number)%shellL(shell_number)
      S_1     = asf_set_LS%csf_set_LS%csf(LS_number)%shellS(shell_number)
      !
      wa = lsj_coefficient_LS_jj(l_shell,N_LS,W_1,Q_1,L_1,S_1,J_1,      &
        	                  jj_minus,N_minus,Q_minus,J_1_minus,    &
        	                  jj_plius,Q_plius,J_1_plius)
      if (abs(wa) > eps20) then
         total_number =                                                  &
	     asf_set_LS%csf_set_LS%nwshells - asf_set_LS%csf_set_LS%nwcore
         if (total_number == 1) then
            wa = wa
         else if (total_number >= 2) then
            do shell_number = asf_set_LS%csf_set_LS%nwcore + 2,          &
	                      asf_set_LS%csf_set_LS%nwshells
               if (abs(wa) > eps20) then
                  wa = wa * lsj_transformation_LS_jj_insade             &
                         (shell_number,jj_number,LS_number)
               end if
            end do
         else
            wa = zero
         end if
      else
         wa = zero
      end if
      !
   end function lsj_transformation_LS_jj_gen
   !
   !
   function lsj_transformation_LS_jj_insade(shell_number,jj_number,LS_number)&
                                                                     result(wa)
   !--------------------------------------------------------------------
   ! Return the value of main part of the transformation matrix from 
   ! jj- to LS-coupling scheme in the case of any number of open shells.
   !
   ! Calls: lsj_coefficient_LS_jj, wigner_6j_symbol wigner_6j_triangle,
   ! wigner_9j_symbol, wigner_9j_triangle.
   !--------------------------------------------------------------------
      !
      integer, intent(in)    :: shell_number,jj_number,LS_number
      !
      real(kind=dp) :: wa, wa_sum
      integer       :: delta_J, number_minus, number_plius, number_plius_1,&
                       number_minus_1,                                     &
                       jj_minus, N_minus, Q_minus,                         &
		       jj_plius, N_plius, Q_plius,                         &
		       J_i_min, J_i_max, J_i_minus,                        &
		       J_i_plius, J_i, J_1_i, Jp_1_i, J_1_i1,              &
                       l_shell, N_LS, W_i, Q_i, L_i, S_i,                  &
		       L_1_i, S_1_i, L_1_i1, S_1_i1
      !
      wa     = zero;      wa_sum = zero
      !
      number_minus = asf_set_LS%csf_set_LS%parent(shell_number)%parent_minus
      number_plius = asf_set_LS%csf_set_LS%parent(shell_number)%parent_plius
      if (number_minus+1 /= number_plius .and.  &
          number_minus*number_plius  /= 0) then
         stop "lsj_transformation_LS_jj_insade(): program stop A."
      end if
      !
      number_plius_1=asf_set_LS%csf_set_LS%parent(shell_number-1)%parent_plius
      if (number_plius_1 == 0) then
         number_plius_1 =                                            &
	    asf_set_LS%csf_set_LS%parent(shell_number-1)%parent_minus
      end if
      !
      if (number_minus == 0) then
         jj_minus  = iabs(asf_set_jj%csf_set%subshell(number_plius)%kappa)*2 - 1
         N_minus   = 0;    J_i_minus = 0;    Q_minus = (jj_minus + 1)/2
         Jp_1_i    = asf_set_jj%csf_set%csf(jj_number)%subshellX(number_plius_1)
      else
         jj_minus  = iabs(asf_set_jj%csf_set%subshell(number_minus)%kappa)*2 - 1
         N_minus   = asf_set_jj%csf_set%csf(jj_number)%occupation(number_minus)
         Q_minus   = (jj_minus + 1)/2 -                                 &
                 asf_set_jj%csf_set%csf(jj_number)%seniority(number_minus)
         J_i_minus = asf_set_jj%csf_set%csf(jj_number)%subshellJ(number_minus)
         Jp_1_i  = asf_set_jj%csf_set%csf(jj_number)%subshellX(number_minus)
      end if
      !
      if (number_plius == 0) then
         jj_plius  = iabs(asf_set_jj%csf_set%subshell(number_minus)%kappa)*2 - 1
         Q_plius   = (jj_plius + 1)/2
         N_plius   = 0;         J_i_plius = 0;         J_1_i   = Jp_1_i
      else
         jj_plius  = iabs(asf_set_jj%csf_set%subshell(number_plius)%kappa)*2 - 1
         N_plius   =  asf_set_jj%csf_set%csf(jj_number)%occupation(number_plius)
         Q_plius   = (jj_plius + 1)/2 -                                  &
                 asf_set_jj%csf_set%csf(jj_number)%seniority(number_plius)
         J_i_plius = asf_set_jj%csf_set%csf(jj_number)%subshellJ(number_plius)
         J_1_i   = asf_set_jj%csf_set%csf(jj_number)%subshellX(number_plius)
      end if
      !
      J_i_min = iabs(J_i_minus - J_i_plius);   J_i_max = J_i_minus + J_i_plius
      !
      l_shell = asf_set_LS%csf_set_LS%shell(shell_number)%l
      N_LS    = asf_set_LS%csf_set_LS%csf(LS_number)%occupation(shell_number)
      W_i     = asf_set_LS%csf_set_LS%csf(LS_number)%w(shell_number)
      Q_i     = 2 * l_shell + 1 -                                          &
                asf_set_LS%csf_set_LS%csf(LS_number)%seniority(shell_number)
      if (N_LS == N_minus + N_plius) then
         L_i = asf_set_LS%csf_set_LS%csf(LS_number)%shellL(shell_number)
         S_i = asf_set_LS%csf_set_LS%csf(LS_number)%shellS(shell_number)
         !
         L_1_i = asf_set_LS%csf_set_LS%csf(LS_number)%shellLX(shell_number)
         S_1_i = asf_set_LS%csf_set_LS%csf(LS_number)%shellSX(shell_number)
         if (shell_number == 2) then
            number_minus_1  =                                           &
	       asf_set_LS%csf_set_LS%parent(shell_number-1)%parent_minus
            if (number_minus_1 == 0) then
               number_minus_1 =                                         &
	       asf_set_LS%csf_set_LS%parent(shell_number-1)%parent_plius
            end if
            J_1_i1 =asf_set_jj%csf_set%csf(jj_number)%subshellJ(number_minus_1)
            L_1_i1 =asf_set_LS%csf_set_LS%csf(LS_number)%shellL(shell_number-1)
            S_1_i1 =asf_set_LS%csf_set_LS%csf(LS_number)%shellS(shell_number-1)
         else
            J_1_i1 =asf_set_jj%csf_set%csf(jj_number)%subshellX(number_plius_1)
            L_1_i1 =asf_set_LS%csf_set_LS%csf(LS_number)%shellLX(shell_number-1)
            S_1_i1 =asf_set_LS%csf_set_LS%csf(LS_number)%shellSX(shell_number-1)
         end if
         !
         do J_i = J_i_min, J_i_max, 2
            delta_J = wigner_6j_triangle(J_i_minus,J_i_plius,J_i,      &
		                         J_1_i,    J_1_i1,   Jp_1_i)
            if (delta_J /= 0) then
               delta_J = wigner_9j_triangle(L_1_i1,S_1_i1,J_1_i1,      &
                                            L_i,   S_i,   J_i,         &
                                            L_1_i, S_1_i, J_1_i)
               if (delta_J /= 0) then
                  wa_sum = wa_sum + (J_i + one) *                         &
                      wigner_6j_symbol(J_i_minus,J_i_plius,J_i,           &
      		                       J_1_i,    J_1_i1,   Jp_1_i,.true.)*&
                      wigner_9j_symbol(L_1_i1,S_1_i1,J_1_i1,              &
                                       L_i,   S_i,   J_i,                 &
                                       L_1_i, S_1_i, J_1_i,.true.)       *&
                      lsj_coefficient_LS_jj(l_shell,N_LS,W_i,Q_i,L_i,S_i,J_i, &
		                       jj_minus,N_minus,Q_minus,J_i_minus,&
                                       jj_plius,Q_plius,J_i_plius)
               end if
            end if
         end do
         wa = wa_sum * sqrt((Jp_1_i+one)*(J_1_i1+one)*(L_1_i+one)*(S_1_i+one))
         if (mod(J_i_minus+J_i_plius+J_1_i1+J_1_i,4) /= 0) wa = - wa
      else
         wa = zero
      end if
      !
   end function lsj_transformation_LS_jj_insade
   !
   !
   function lsj_triangle(i1,i2,i3)                            result(yes)
   !--------------------------------------------------------------------
   ! Returns .true. if the lengths i1, i2, and i3 may form a triangle 
   ! and .false. otherwise.
   !
   ! Calls:
   !--------------------------------------------------------------------
      !
      integer(kind=i1b), intent(in) :: i1, i2, i3
      !
      logical                       :: yes
      !
      if (i1 <= i2+i3  .and.  i2 <= i3+i1  .and.  i3 <= i1+i2) then
         yes = .true.
      else
         yes = .false.
      end if
      !
   end function lsj_triangle
   !
end module rabs_lsj
