module rabs_XL
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module implements a variety of interaction strength which are used
! in different programs. The main differences of the given subprocedures refer
! to different interactions as well as orbital sets used in the computations.
!-----------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   use rabs_grasp2k
   use rabs_nucleus
   implicit none
   !
   public  :: W_integral  
                 ! Calculates W_mu^nu basic integral for given type mu, 
		 ! rank nu, and quantum numbers a, b, c, and d.
   public  :: XL_Breit0_coefficients  
                 ! Determines and calculates the coefficients of the basic
                 ! integral arrays for the effective Breit strength 
                 ! X^L_Breit(abcd;omega=0).
   private :: XL_Breit0_strength
                 ! Calculates the ordinary Breit effective interaction strength 
                 ! X^L_Breit(abcd).   
   public  :: XL_Breit0_strength_grasp2k
                 ! Calculates the ordinary Breit effective interaction strength 
                 ! X^L_Breit(abcd;omega=0) using GRASP92 wave functions.   
   public  :: XL_Breit_coefficients  
                 ! Determines and calculates the coefficients of the basic
                 ! integral arrays for the effective (full transverse) Breit 
                 ! strength X^L_Breit(abcd).
   public  :: XL_Breit_strength_grasp2k
                 ! Calculates the (full transverse) Breit effective interaction
                 ! strength  X^L_Breit(abcd) using GRASP92 wave functions.   
   public  :: XL_Coulomb_coefficients  
                 ! Determines and calculates the coefficients of the basic
                 ! integral arrays for the effective Coulomb strength 
                 ! X^L_Coulomb(abcd).
   private :: XL_Coulomb_strength
                 ! Calculates the Coulomb effective interaction strength 
                 ! X^L_Coulomb(abcd).   
   public  :: XL_Coulomb_strength_grasp2k
                 ! Calculates the Coulomb effective interaction strength 
                 ! X^L_Coulomb(abcd) using GRASP92 wave functions.   
   public  :: XL_Debye_strength_grasp2k
                 ! Calculates the Debye effective interaction strength 
                 ! X^L_Debye(abcd) for the description of plasma interactions
                 ! using GRASP92 wave functions.   
   private :: XL_Gaunt_coefficients  
                 ! Determines and calculates the coefficients of the basic
                 ! integral arrays for the effective Gaunt strength 
                 ! X^L_Gaunt(abcd).
   private :: XL_Gaunt_strength
                 ! Calculates the Gaunt effective interaction strength 
                 ! X^L_Gaunt(abcd).   
   public  :: XL_Gaunt_strength_grasp2k
                 ! Calculates the Gaunt effective interaction strength 
                 ! X^L_Gaunt(abcd) using GRASP92 wave functions.   
   public  :: XL_sms_strength_grasp2k
                 ! Calculates the specific mass shift effective interaction 
                 ! strength X^0_sms(abcd) using GRASP92 wave functions.   
   private :: XL_strength
                 ! Returns the effective interaction strength X^L(abcd)
                 ! for a given set of interacting orbitals  and type of the
                 ! interaction.   
   public  :: XL_strength_grasp2k
                 ! Returns the effective interaction strength X^L(abcd)
                 ! for a given set of interacting orbitals  and type of the
                 ! interaction.   
   !
   type, public :: XLcoefficient
      integer       :: mu, nu
      type(nkappa)  :: a, c, b, d
      real(kind=dp) :: coeff
   end type XLcoefficient
   !
   !
   ! Define global logical flags for the control of the RELCI program; the
   ! default values for these flags may be overwritten interactively during 
   ! input time
   !
contains
   !
   function W_integral(mu,nu,a,c,b,d,silent)               result(Wmunu)
   !--------------------------------------------------------------------
   ! Calculates for given type mu, rank nu, and quantum numbers a, b, c,
   ! and d the value of the corresponding W_mu^nu basic integral.
   ! The orbital quantum numbers are a = (pqna, kapa), b = (pqnb, kapb),
   ! etc, all of type(nkappa).
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)      :: mu, nu
      type(nkappa), intent(in) :: a, b, c, d
      logical, optional        :: silent
      real(kind=dp)            :: Wmunu
      !
      Wmunu = one
      print *, "mu,nu,a,c,b,d,silent = ",mu,nu,a,c,b,d,silent
      print *, " "
      !
      print *, "rabs_use_basis_set = ",rabs_use_basis_set
      stop "W_integral() in rabs_xl: program stop A."
      !
   end function W_integral
   !
   !
   subroutine XL_Breit0_coefficients(L,a,b,c,d,xlcoeff,nx,nxmax)
   !--------------------------------------------------------------------
   ! Determines and calculates the number, weights, and identifiers of 
   ! the basic W_mu^nu (...) integral arrays which contribute to effective
   ! Breit interaction strengths X^L_Breit (abcd).
   ! This routine follows the paper of I P Grant and B J McKenzie, 
   ! J. Phys. B13, 2671-81 (1980), eq.(5).
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)      :: L, nxmax
      integer, intent(out)     :: nx
      type(nkappa), intent(in) :: a, b, c, d
      type(XLcoefficient), dimension(1:nxmax), intent(out) :: xlcoeff
      !
      integer       :: ja, jb, jc, jd, la, lb, lc, ld, nu
      real(kind=dp) :: wa, wb, xc, xcc
      !
      ! Ensure a big enough dimension of the xlcoeff array
      if (nxmax < 40) then
         print *, "nx, nxmax = ",nx, nxmax
         stop     "XL_Breit0_coefficients(): program stop A."
      end if
      !
      ! Initializes the number of XL_Breit coefficients
      nx = 0
      !
      ja = angular_momentum_j(a%kappa);  jb = angular_momentum_j(b%kappa) 
      jc = angular_momentum_j(c%kappa);  jd = angular_momentum_j(d%kappa) 
      !
      la = angular_momentum_l(a%kappa);  lb = angular_momentum_l(b%kappa) 
      lc = angular_momentum_l(c%kappa);  ld = angular_momentum_l(d%kappa) 
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0) then
         return
      end if 
      !
      xc = CL_reduced_me(a%kappa,L,c%kappa) * CL_reduced_me(b%kappa,L,d%kappa)
      if (mod(L,2) == 1) then
         xc = - xc
      end if
      !
      if (abs(xc) .lt. eps10) then
         return
      end if
      !
      ! Consider the individual contributions from sum_nu and sum_mu.
      ! First, take T^(nu,L)_mu = R^(nu,L)_mu
      ! 
      nu = L - 1
      !
      if (mod(la+lc+nu,2) == 1   .and.   &
          mod(lb+ld+nu,2) == 1   .and.   L /= 0) then
         wa = (L+one) / ( L*(L+L-one)*(L+L+one) )
         !
         ! mu = 1
         xcc = xc * wa * (c%kappa-a%kappa+L) * (d%kappa-b%kappa+L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 2
         xcc = xc * wa * (c%kappa-a%kappa-L) * (d%kappa-b%kappa-L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 3
         xcc = xc * wa * (c%kappa-a%kappa+L) * (d%kappa-b%kappa-L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 4
         xcc = xc * wa * (c%kappa-a%kappa-L) * (d%kappa-b%kappa+L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = xcc
         end if
      end if
      !
      nu = L
      !
      if (mod(la+lc+nu,2) == 1   .and.   &
          mod(lb+ld+nu,2) == 1   .and.   L /= 0) then
         wa = - real(a%kappa+c%kappa,kind=dp) * (b%kappa+d%kappa ) / &
                (L*(L+one))
         !
         ! mu = 1
         xcc = xc * wa 
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 2
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 3
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = xcc
         endif
         !
         ! mu = 4
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = xcc
         end if
      end if
      !
      nu = L + 1
      !
      if (mod(la+lc+nu,2) == 1   .and.   &
          mod(lb+ld+nu,2) == 1   .and.    L /= 0) then
         wa = real(L,kind=dp) / ( (L+one)*(L+L+one)*(L+L+three) )
         !
         ! mu = 1
         xcc = xc * wa * (c%kappa-a%kappa-L-one) * (d%kappa-b%kappa-L-one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 2
         xcc = xc * wa * (c%kappa-a%kappa+L+one) * (d%kappa-b%kappa+L+one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 3
         xcc = xc * wa * (c%kappa-a%kappa-L-one) * (d%kappa-b%kappa+L+one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 4
         xcc = xc * wa * (c%kappa-a%kappa+L+one) * (d%kappa-b%kappa-L-one) 
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = xcc
         end if
      end if
      !
      ! Add contributions of the S^k_mu integrals
      !
      ! mu = 1
      if (mod(la+lc+L-1,2) == 1   .and.   mod(lb+ld+L+1,2) == 1) then
         wb  = one / ( (L+L+one)*(L+L+one) ) 
         xcc = xc * wb * (c%kappa-a%kappa+L) * (d%kappa-b%kappa-L-one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 2
         xcc = xc * wb * (d%kappa-b%kappa+L) * (c%kappa-a%kappa-L-one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 3
         xcc = xc * wb * (d%kappa-b%kappa+L+one) * (c%kappa-a%kappa-L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 4
         xcc = xc * wb * (d%kappa-b%kappa-L) * (c%kappa-a%kappa+L+one) 
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 5
         xcc = xc * wb * (d%kappa-b%kappa+L+one) * (c%kappa-a%kappa+L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 6
         xcc = xc * wb * (d%kappa-b%kappa-L) * (c%kappa-a%kappa-L-one) 
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 7
         xcc = xc * wb * (d%kappa-b%kappa-L-one) * (c%kappa-a%kappa-L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 8
         xcc = xc * wb * (d%kappa-b%kappa+L) * (c%kappa-a%kappa+L+one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
      end if
      !
      if (nx > nxmax) then
         stop "XL_Breit0_coefficients(): program stop B."
      end if
      !
   end subroutine XL_Breit0_coefficients
   !
   !
   function XL_Breit0_strength(L,a,b,c,d)                result(XL_Breit)
   !--------------------------------------------------------------------
   ! Calculates the effective Breit interaction strengths 
   ! X^L_Breit (abcd) for given rank L and orbital quantum numbers
   ! a, b, c, and d.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)      :: L
      type(nkappa), intent(in) :: a, b, c, d
      real(kind=dp)            :: XL_Breit
      !
      integer :: i, ja, jb, jc, jd, la, lb, lc, ld, nx, nxmax
      type(XLcoefficient), dimension(:), pointer :: xlcoeff
      !
      XL_Breit = zero
      !
      ja = angular_momentum_j(a%kappa);  jb = angular_momentum_j(b%kappa) 
      jc = angular_momentum_j(c%kappa);  jd = angular_momentum_j(d%kappa) 
      !
      la = angular_momentum_l(a%kappa);  lb = angular_momentum_l(b%kappa) 
      lc = angular_momentum_l(c%kappa);  ld = angular_momentum_l(d%kappa) 
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0) then
         return
      end if 
      !
      ! Allocate an array for calculating coefficients and calculate them
      nxmax = 42;   nx = 0
      allocate( xlcoeff(1:nxmax) )
      call XL_Breit0_coefficients(L,a,b,c,d,xlcoeff,nx,nxmax)
      !
      !
      if (rabs_use_basis_set) then
         do  i = 1,nx
            XL_Breit = XL_Breit + xlcoeff(i)%coeff *           &
                       W_integral(xlcoeff(i)%mu,xlcoeff(i)%nu, &
                                  xlcoeff(i)%a,xlcoeff(i)%c,   &
                                  xlcoeff(i)%b,xlcoeff(i)%d,.true.)
            !!                      xlcoeff(i)%b,xlcoeff(i)%d,silent=.true.)
         end do
      end if
      !
      deallocate( xlcoeff )
      !
   end function XL_Breit0_strength
   !
   !
   function XL_Breit0_strength_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d)       &
                                                        result(XL_Breit)
   !--------------------------------------------------------------------
   ! Calculates the effective Breit interaction strengths 
   ! X^L_Gaunt (abcd) for given rank L and orbital quantum numbers
   ! a, b, c, and d using GRASP92 wave functions.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)                        :: L
      type(orbital_function), intent(in), target :: rwf_a, rwf_b, rwf_c, rwf_d
      real(kind=dp)                              :: XL_Breit
      !
      integer, parameter :: nxmax = 100
      integer            :: i, ja, jb, jc, jd, nx
      real(kind=dp)      :: wa
      type(orbital_function), pointer         :: rwf_aa, rwf_bb, rwf_cc, rwf_dd
      type(XLcoefficient), dimension(1:nxmax) :: xlcoeff
      !
      XL_Breit = zero
      !
      ja = angular_momentum_j(rwf_a%orbital%kappa)
      jb = angular_momentum_j(rwf_b%orbital%kappa) 
      jc = angular_momentum_j(rwf_c%orbital%kappa)
      jd = angular_momentum_j(rwf_d%orbital%kappa) 
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0   .or.  &
          L == 0) then
         return
      end if 
      !
      ! Calculate coefficients
      nx = 0
      call XL_Breit0_coefficients(L,rwf_a%orbital,rwf_b%orbital,rwf_c%orbital, &
                                    rwf_d%orbital,xlcoeff,nx,nxmax)
      !
      do  i = 1,nx
         if (xlcoeff(i)%a%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%a%kappa == rwf_a%orbital%kappa) then
            rwf_aa => rwf_a
         else if (xlcoeff(i)%a%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%a%kappa == rwf_b%orbital%kappa) then
            rwf_aa => rwf_b
         else if (xlcoeff(i)%a%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%a%kappa == rwf_c%orbital%kappa) then
            rwf_aa => rwf_c
         else if (xlcoeff(i)%a%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%a%kappa == rwf_d%orbital%kappa) then
            rwf_aa => rwf_d
         else
            stop "XL_Breit_strength_grasp2k(): program stop A."
         end if
         !
         if (xlcoeff(i)%b%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%b%kappa == rwf_a%orbital%kappa) then
            rwf_bb => rwf_a
         else if (xlcoeff(i)%b%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%b%kappa == rwf_b%orbital%kappa) then
            rwf_bb => rwf_b
         else if (xlcoeff(i)%b%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%b%kappa == rwf_c%orbital%kappa) then
            rwf_bb => rwf_c
         else if (xlcoeff(i)%b%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%b%kappa == rwf_d%orbital%kappa) then
            rwf_bb => rwf_d
         else
            stop "XL_Breit_strength_grasp2k(): program stop B."
         end if
         !
         if (xlcoeff(i)%c%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%c%kappa == rwf_a%orbital%kappa) then
            rwf_cc => rwf_a
         else if (xlcoeff(i)%c%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%c%kappa == rwf_b%orbital%kappa) then
            rwf_cc => rwf_b
         else if (xlcoeff(i)%c%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%c%kappa == rwf_c%orbital%kappa) then
            rwf_cc => rwf_c
         else if (xlcoeff(i)%c%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%c%kappa == rwf_d%orbital%kappa) then
            rwf_cc => rwf_d
         else
            stop "XL_Breit_strength_grasp2k(): program stop C."
         end if
         !
         if (xlcoeff(i)%d%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%d%kappa == rwf_a%orbital%kappa) then
            rwf_dd => rwf_a
         else if (xlcoeff(i)%d%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%d%kappa == rwf_b%orbital%kappa) then
            rwf_dd => rwf_b
         else if (xlcoeff(i)%d%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%d%kappa == rwf_c%orbital%kappa) then
            rwf_dd => rwf_c
         else if (xlcoeff(i)%d%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%d%kappa == rwf_d%orbital%kappa) then
            rwf_dd => rwf_d
         else
            stop "XL_Breit_strength_grasp2k(): program stop D."
         end if
         !
         wa = W_integral_grasp2k(xlcoeff(i)%mu,xlcoeff(i)%nu,             &
                                               rwf_aa,rwf_cc,rwf_bb,rwf_dd)
         !
         XL_Breit = XL_Breit + xlcoeff(i)%coeff * wa
      end do 
      !
   end function XL_Breit0_strength_grasp2k
   !
   !
   subroutine XL_Breit_coefficients(L,a,b,c,d,xlcoeff,nx,nxmax)
   !--------------------------------------------------------------------
   ! Determines and calculates the number, weights, and identifiers of 
   ! the basic W_mu^nu (...) integral arrays which contribute to effective
   ! (full transverse) Breit interaction strengths X^L_Breit (abcd).
   ! This routine follows the paper of I P Grant and B J McKenzie, 
   ! J. Phys. B13, 2671-81 (1980), eq.(5).
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)      :: L, nxmax
      integer, intent(out)     :: nx
      type(nkappa), intent(in) :: a, b, c, d
      type(XLcoefficient), dimension(1:nxmax), intent(out) :: xlcoeff
      !
      integer       :: ja, jb, jc, jd, la, lb, lc, ld, nu
      real(kind=dp) :: wa, wb, xc, xcc
      !
      ! Ensure a big enough dimension of the xlcoeff array
      if (nxmax < 40) then
         print *, "nx, nxmax = ",nx, nxmax
         stop     "XL_Breit_coefficients(): program stop A."
      end if
      !
      ! Initializes the number of XL_Breit coefficients
      nx = 0
      !
      ja = angular_momentum_j(a%kappa);  jb = angular_momentum_j(b%kappa) 
      jc = angular_momentum_j(c%kappa);  jd = angular_momentum_j(d%kappa) 
      !
      la = angular_momentum_l(a%kappa);  lb = angular_momentum_l(b%kappa) 
      lc = angular_momentum_l(c%kappa);  ld = angular_momentum_l(d%kappa) 
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0) then
         return
      end if 
      !
      xc = CL_reduced_me(a%kappa,L,c%kappa) * CL_reduced_me(b%kappa,L,d%kappa)
      if (mod(L,2) == 1) then
         xc = - xc
      end if
      !
      if (abs(xc) .lt. eps10) then
         return
      end if
      !
      ! Consider the individual contributions from sum_nu and sum_mu.
      ! First, take T^(nu,L)_mu = R^(nu,L)_mu
      ! 
      nu = L - 1
      !
      if (mod(la+lc+nu,2) == 1   .and.   &
          mod(lb+ld+nu,2) == 1   .and.   L /= 0) then
         wa = (L+one) / ( L*(L+L-one)*(L+L+one) )
         !
         ! mu = 1
         xcc = xc * wa * (c%kappa-a%kappa+L) * (d%kappa-b%kappa+L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 2
         xcc = xc * wa * (c%kappa-a%kappa-L) * (d%kappa-b%kappa-L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 3
         xcc = xc * wa * (c%kappa-a%kappa+L) * (d%kappa-b%kappa-L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 4
         xcc = xc * wa * (c%kappa-a%kappa-L) * (d%kappa-b%kappa+L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
         end if
      end if
      !
      nu = L
      !
      if (mod(la+lc+nu,2) == 1   .and.   &
          mod(lb+ld+nu,2) == 1   .and.   L /= 0) then
         wa = - real(a%kappa+c%kappa,kind=dp) * (b%kappa+d%kappa ) / &
                (L*(L+one))
         !
         ! mu = 1
         xcc = xc * wa 
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 2
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 3
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
         endif
         !
         ! mu = 4
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
         end if
      end if
      !
      nu = L + 1
      !
      if (mod(la+lc+nu,2) == 1   .and.   &
          mod(lb+ld+nu,2) == 1   .and.    L /= 0) then
         wa = real(L,kind=dp) / ( (L+one)*(L+L+one)*(L+L+three) )
         !
         ! mu = 1
         xcc = xc * wa * (c%kappa-a%kappa-L-one) * (d%kappa-b%kappa-L-one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 2
         xcc = xc * wa * (c%kappa-a%kappa+L+one) * (d%kappa-b%kappa+L+one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 3
         xcc = xc * wa * (c%kappa-a%kappa-L-one) * (d%kappa-b%kappa+L+one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 4
         xcc = xc * wa * (c%kappa-a%kappa+L+one) * (d%kappa-b%kappa-L-one) 
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 6;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
         end if
      end if
      !
      goto 1
      !---------------------------
      ! Add contributions of the S^k_mu integrals
      !
      ! mu = 1
      if (mod(la+lc+L-1,2) == 1   .and.   mod(lb+ld+L+1,2) == 1) then
         wb  = one / ( (L+L+one)*(L+L+one) ) 
         xcc = xc * wb * (c%kappa-a%kappa+L) * (d%kappa-b%kappa-L-one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 2
         xcc = xc * wb * (d%kappa-b%kappa+L) * (c%kappa-a%kappa-L-one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 3
         xcc = xc * wb * (d%kappa-b%kappa+L+one) * (c%kappa-a%kappa-L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 4
         xcc = xc * wb * (d%kappa-b%kappa-L) * (c%kappa-a%kappa+L+one) 
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 5
         xcc = xc * wb * (d%kappa-b%kappa+L+one) * (c%kappa-a%kappa+L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 6
         xcc = xc * wb * (d%kappa-b%kappa-L) * (c%kappa-a%kappa-L-one) 
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 7
         xcc = xc * wb * (d%kappa-b%kappa-L-one) * (c%kappa-a%kappa-L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
         !
         ! mu = 8
         xcc = xc * wb * (d%kappa-b%kappa+L) * (c%kappa-a%kappa+L+one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L+1  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = (L+L+one) / two * xcc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L-1  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = -(L+L+one) / two * xcc
         end if
      end if
      !
      goto 1
      !---------------------------
      ! Add contributions of the S^k_mu integrals
      !
      ! mu = 1
      if (mod(la+lc+L-1,2) == 1   .and.   mod(lb+ld+L+1,2) == 1) then
         wb  = one / ( (L+L+one)*(L+L+one) ) 
         xcc = xc * wb * (c%kappa-a%kappa+L) * (d%kappa-b%kappa-L-one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 2
         xcc = xc * wb * (d%kappa-b%kappa+L) * (c%kappa-a%kappa-L-one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 3
         xcc = xc * wb * (d%kappa-b%kappa+L+one) * (c%kappa-a%kappa-L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L 
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 4
         xcc = xc * wb * (d%kappa-b%kappa-L) * (c%kappa-a%kappa+L+one) 
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 5
         xcc = xc * wb * (d%kappa-b%kappa+L+one) * (c%kappa-a%kappa+L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 6
         xcc = xc * wb * (d%kappa-b%kappa-L) * (c%kappa-a%kappa-L-one) 
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 7
         xcc = xc * wb * (d%kappa-b%kappa-L-one) * (c%kappa-a%kappa-L)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L 
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = xcc
         end if
         !
         ! mu = 8
         xcc = xc * wb * (d%kappa-b%kappa+L) * (c%kappa-a%kappa+L+one)
         if (abs(xcc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 7;   xlcoeff(nx)%nu = L  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = xcc
         end if
      end if
      !
    1 if (nx > nxmax) then
         stop "XL_Breit_coefficients(): program stop B."
      end if
      !
   end subroutine XL_Breit_coefficients
   !
   !
   function XL_Breit_strength_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d)       &
                                                        result(XL_Breit)
   !--------------------------------------------------------------------
   ! Calculates the effective Breit interaction strengths 
   ! X^L_Gaunt (abcd) for given rank L and orbital quantum numbers
   ! a, b, c, and d using GRASP92 wave functions.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)                        :: L
      type(orbital_function), intent(in), target :: rwf_a, rwf_b, rwf_c, rwf_d
      real(kind=dp)                              :: XL_Breit
      !
      integer, parameter :: nxmax = 100
      integer            :: i, ja, jb, jc, jd, nx
      real(kind=dp)      :: wa
      type(orbital_function), pointer         :: rwf_aa, rwf_bb, rwf_cc, rwf_dd
      type(XLcoefficient), dimension(1:nxmax) :: xlcoeff
      !
      XL_Breit = zero
      !
      ja = angular_momentum_j(rwf_a%orbital%kappa)
      jb = angular_momentum_j(rwf_b%orbital%kappa) 
      jc = angular_momentum_j(rwf_c%orbital%kappa)
      jd = angular_momentum_j(rwf_d%orbital%kappa) 
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0   .or.  &
          L == 0) then
         return
      end if 
      !
      ! Calculate coefficients
      nx = 0
      call XL_Breit_coefficients(L,rwf_a%orbital,rwf_b%orbital,rwf_c%orbital, &
                                   rwf_d%orbital,xlcoeff,nx,nxmax)
      !
      do  i = 1,nx
         if (xlcoeff(i)%a%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%a%kappa == rwf_a%orbital%kappa) then
            rwf_aa => rwf_a
         else if (xlcoeff(i)%a%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%a%kappa == rwf_b%orbital%kappa) then
            rwf_aa => rwf_b
         else if (xlcoeff(i)%a%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%a%kappa == rwf_c%orbital%kappa) then
            rwf_aa => rwf_c
         else if (xlcoeff(i)%a%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%a%kappa == rwf_d%orbital%kappa) then
            rwf_aa => rwf_d
         else
            stop "XL_Breit_strength_grasp2k(): program stop A."
         end if
         !
         if (xlcoeff(i)%b%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%b%kappa == rwf_a%orbital%kappa) then
            rwf_bb => rwf_a
         else if (xlcoeff(i)%b%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%b%kappa == rwf_b%orbital%kappa) then
            rwf_bb => rwf_b
         else if (xlcoeff(i)%b%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%b%kappa == rwf_c%orbital%kappa) then
            rwf_bb => rwf_c
         else if (xlcoeff(i)%b%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%b%kappa == rwf_d%orbital%kappa) then
            rwf_bb => rwf_d
         else
            stop "XL_Breit_strength_grasp2k(): program stop B."
         end if
         !
         if (xlcoeff(i)%c%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%c%kappa == rwf_a%orbital%kappa) then
            rwf_cc => rwf_a
         else if (xlcoeff(i)%c%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%c%kappa == rwf_b%orbital%kappa) then
            rwf_cc => rwf_b
         else if (xlcoeff(i)%c%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%c%kappa == rwf_c%orbital%kappa) then
            rwf_cc => rwf_c
         else if (xlcoeff(i)%c%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%c%kappa == rwf_d%orbital%kappa) then
            rwf_cc => rwf_d
         else
            stop "XL_Breit_strength_grasp2k(): program stop C."
         end if
         !
         if (xlcoeff(i)%d%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%d%kappa == rwf_a%orbital%kappa) then
            rwf_dd => rwf_a
         else if (xlcoeff(i)%d%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%d%kappa == rwf_b%orbital%kappa) then
            rwf_dd => rwf_b
         else if (xlcoeff(i)%d%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%d%kappa == rwf_c%orbital%kappa) then
            rwf_dd => rwf_c
         else if (xlcoeff(i)%d%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%d%kappa == rwf_d%orbital%kappa) then
            rwf_dd => rwf_d
         else
            stop "XL_Breit_strength_grasp2k(): program stop D."
         end if
         !
         print *, " "
         print *, rwf_a%orbital%kappa,rwf_c%orbital%kappa, &
                  rwf_b%orbital%kappa,rwf_d%orbital%kappa
         print *, rwf_aa%orbital%kappa,rwf_cc%orbital%kappa, &
                  rwf_bb%orbital%kappa,rwf_dd%orbital%kappa
         wa = W_integral_grasp2k(xlcoeff(i)%mu,xlcoeff(i)%nu,          &
                                            rwf_aa,rwf_cc,rwf_bb,rwf_dd)
         !
         XL_Breit = XL_Breit + xlcoeff(i)%coeff * wa
      end do 
      !
   end function XL_Breit_strength_grasp2k
   !
   !
   subroutine XL_Coulomb_coefficients(L,a,b,c,d,xlcoeff,nx,nxmax,      &
                                                       Slater_integrals)
   !--------------------------------------------------------------------
   ! Determines and calculates the number, weights, and identifiers of 
   ! the basic W_mu^nu (...) integral arrays which contribute to effective
   ! Coulomb interaction strengths X^L_Coulomb (abcd).
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)      :: L, nxmax
      integer, intent(inout)   :: nx
      type(nkappa), intent(in) :: a, b, c, d
      type(XLcoefficient), dimension(1:nxmax), intent(out) :: xlcoeff
      logical, intent(in), optional                        :: Slater_integrals
      !
      integer       :: ja, jb, jc, jd, la, lb, lc, ld, mu, nxc
      real(kind=dp) :: xc
      !
      ! Initializes the number of XL_Coulomb coefficients
      nx = 0
      !
      ja = angular_momentum_j(a%kappa);  jb = angular_momentum_j(b%kappa) 
      jc = angular_momentum_j(c%kappa);  jd = angular_momentum_j(d%kappa) 
      !
      la = angular_momentum_l(a%kappa);  lb = angular_momentum_l(b%kappa) 
      lc = angular_momentum_l(c%kappa);  ld = angular_momentum_l(d%kappa) 
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0  .or. &
          mod(la+lc+L,2) == 1   .or.   mod(lb+ld+L,2) == 1) then
         return
      end if 
      !
      xc = CL_reduced_me(a%kappa,L,c%kappa) * CL_reduced_me(b%kappa,L,d%kappa)
      if (mod(L,2) == 1) then
         xc = - xc
      end if
      !
      if (abs(xc) .lt. eps10) then
         return
      end if
      !
      ! Collect the contributions from the basic integrals 
      ! W_1^L(..), ..., W_4^L(..)
      !
      if (present(Slater_integrals)   .and.   Slater_integrals) then
         nx = 1
      else
         nx = 8
      end if
      !
      if (nx > nxmax) then 
         stop "XL_Coulomb_coeffients(): program stop A."
      end if
      !
      if (present(Slater_integrals)   .and.   Slater_integrals) then
         xlcoeff(1)%mu    = 0
         xlcoeff(1)%nu    = L
         xlcoeff(1)%a     = a
         xlcoeff(1)%c     = c
         xlcoeff(1)%b     = b
         xlcoeff(1)%d     = d
         xlcoeff(1)%coeff = xc
      else
         nxc = 0
         do  mu = 1,4
            nxc = nxc + 1
            xlcoeff(nxc)%mu    = mu
            xlcoeff(nxc)%nu    = L
            xlcoeff(nxc)%a     = a
            xlcoeff(nxc)%c     = c
            xlcoeff(nxc)%b     = b
            xlcoeff(nxc)%d     = d
            xlcoeff(nxc)%coeff = xc
         end do
         !
         do  mu = 1,4
            nxc = nxc + 1
            xlcoeff(nxc)%mu    = mu
            xlcoeff(nxc)%nu    = L
            xlcoeff(nxc)%a     = b
            xlcoeff(nxc)%c     = d
            xlcoeff(nxc)%b     = a
            xlcoeff(nxc)%d     = c
            xlcoeff(nxc)%coeff = xc
            !! print *, "XL_Coulomb_coeff: nxc,mu,coeff = ",nxc,xlcoeff(nxc)%mu,  &
            !!                             xlcoeff(nxc)%coeff
         end do
      end if
      !! print *, "XL_Coulomb_coefficients: nxc = ",nxc
      !
   end subroutine XL_Coulomb_coefficients
   !
   !
   function XL_Coulomb_strength(L,a,b,c,d)            result(XL_Coulomb)
   !--------------------------------------------------------------------
   ! Calculates the effective Coulomb interaction strengths 
   ! X^L_Coulomb (abcd) for given rank L and orbital quantum numbers
   ! a, b, c, and d.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)      :: L
      type(nkappa), intent(in) :: a, b, c, d
      real(kind=dp)            :: XL_Coulomb
      !
      integer :: i, ja, jb, jc, jd, la, lb, lc, ld, nx, nxmax
      type(XLcoefficient), dimension(:), pointer :: xlcoeff
      !
      XL_Coulomb = zero
      !
      ja = angular_momentum_j(a%kappa);  jb = angular_momentum_j(b%kappa) 
      jc = angular_momentum_j(c%kappa);  jd = angular_momentum_j(d%kappa) 
      !
      la = angular_momentum_l(a%kappa);  lb = angular_momentum_l(b%kappa) 
      lc = angular_momentum_l(c%kappa);  ld = angular_momentum_l(d%kappa) 
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0  .or. &
          mod(la+lc+L,2) == 1   .or.   mod(lb+ld+L,2) == 1) then
         return
      end if 
      !
      ! Allocate an array for calculating coefficients and calculate them
      nxmax = 10;   nx = 0
      allocate( xlcoeff(1:nxmax) )
      call XL_Coulomb_coefficients(L,a,b,c,d,xlcoeff,nx,nxmax,.false.)
      !!x print *, "XL_Coulomb_strength: nx = ",nx
      !
      !! do  i = 1,nx
      !!    print *, "XL_Coulomb_strength: i,mu,nu,a,c,b,d,coeff = ", &
      !!             i,xlcoeff(i)%mu,xlcoeff(i)%nu,xlcoeff(i)%a,xlcoeff(i)%c,   &
      !!                                           xlcoeff(i)%b,xlcoeff(i)%d,   &
      !!                                           xlcoeff(i)%coeff
      !! end do
      !
      if (rabs_use_basis_set) then
      do  i = 1,nx
         !! print *, "XL_Coulomb_strength: i,mu,nu,a,c,b,d,coeff = ", &
         !!          i,xlcoeff(i)%mu,xlcoeff(i)%nu,xlcoeff(i)%a,xlcoeff(i)%c,   &
         !!                                        xlcoeff(i)%b,xlcoeff(i)%d,   &
         !!                                        xlcoeff(i)%coeff
         XL_Coulomb = XL_Coulomb + xlcoeff(i)%coeff *         &
                      W_integral(xlcoeff(i)%mu,xlcoeff(i)%nu, &
                                 xlcoeff(i)%a,xlcoeff(i)%c,   &
                                 xlcoeff(i)%b,xlcoeff(i)%d,.true.)
         !                        xlcoeff(i)%b,xlcoeff(i)%d,silent=.true.)
         !! print *, "XL_Coulomb_strength: XL_Coulomb = ",XL_Coulomb
      end do
      end if
      !
      deallocate( xlcoeff )
      !
   end function XL_Coulomb_strength
   !
   !
   function XL_Coulomb_strength_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d,is_bound) &
                                                      result(XL_Coulomb)
   !--------------------------------------------------------------------
   ! Calculates the effective Coulomb interaction strengths 
   ! X^L_Coulomb (abcd) for given rank L and orbital functions rwf_a, 
   ! rwf_b, rwf_c, and rwf_d using wave functions from GRASP92.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)                :: L
      type(orbital_function), intent(in) :: rwf_a, rwf_b, rwf_c, rwf_d
      logical, optional                  :: is_bound
      real(kind=dp)                      :: XL_Coulomb
      !
      integer       :: ja, jb, jc, jd, la, lb, lc, ld
      real(kind=dp) :: xc
      !
      XL_Coulomb = zero
      !
      ja = angular_momentum_j(rwf_a%orbital%kappa)
      jb = angular_momentum_j(rwf_b%orbital%kappa) 
      jc = angular_momentum_j(rwf_c%orbital%kappa)
      jd = angular_momentum_j(rwf_d%orbital%kappa) 
      !
      la = angular_momentum_l(rwf_a%orbital%kappa)
      lb = angular_momentum_l(rwf_b%orbital%kappa)
      lc = angular_momentum_l(rwf_c%orbital%kappa)
      ld = angular_momentum_l(rwf_d%orbital%kappa)
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0  .or. &
          mod(la+lc+L,2) == 1   .or.   mod(lb+ld+L,2) == 1) then
         return
      end if 
      !
      xc = CL_reduced_me(rwf_a%orbital%kappa,L,rwf_c%orbital%kappa) *      &
           CL_reduced_me(rwf_b%orbital%kappa,L,rwf_d%orbital%kappa)
      if (mod(L,2) == 1) then
         xc = - xc
      end if
      !
      if (abs(xc) .lt. eps10) then
         return
      end if
      !
      !!! XL_Coulomb = xc * slater_integral_explicit(L,rwf_a,rwf_b,rwf_c,rwf_d)
      !!! return
      !
      if (present(is_bound)) then
         XL_Coulomb = xc * slater_integral_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d,&
                                                                     is_bound)
      else
         XL_Coulomb = xc * slater_integral_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d,&
                                                                       .true.)
      end if
      !                                                                 
   end function XL_Coulomb_strength_grasp2k
   !
   !
   function XL_Debye_strength_grasp2k(L,lambda,rwf_a,rwf_b,rwf_c,rwf_d, &
                                             is_bound)   result(XL_Debye)
   !--------------------------------------------------------------------
   ! Calculates the effective Debye interaction strengths 
   ! X^L_Debye (abcd) for given rank L and orbital functions rwf_a, 
   ! rwf_b, rwf_c, and rwf_d using wave functions from GRASP92.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)                :: L
      real(kind=dp), intent(in)          :: lambda
      type(orbital_function), intent(in) :: rwf_a, rwf_b, rwf_c, rwf_d
      logical, optional                  :: is_bound
      real(kind=dp)                      :: XL_Debye, wa, wb, wc
      !
      integer       :: ja, jb, jc, jd, la, lb, lc, ld
      real(kind=dp) :: xc
      !
      XL_Debye = zero
      !
      ja = angular_momentum_j(rwf_a%orbital%kappa)
      jb = angular_momentum_j(rwf_b%orbital%kappa) 
      jc = angular_momentum_j(rwf_c%orbital%kappa)
      jd = angular_momentum_j(rwf_d%orbital%kappa) 
      !
      la = angular_momentum_l(rwf_a%orbital%kappa)
      lb = angular_momentum_l(rwf_b%orbital%kappa)
      lc = angular_momentum_l(rwf_c%orbital%kappa)
      ld = angular_momentum_l(rwf_d%orbital%kappa)
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0  .or. &
          mod(la+lc+L,2) == 1   .or.   mod(lb+ld+L,2) == 1) then
         return
      end if 
      !
      xc = CL_reduced_me(rwf_a%orbital%kappa,L,rwf_c%orbital%kappa) *      &
           CL_reduced_me(rwf_b%orbital%kappa,L,rwf_d%orbital%kappa)
      if (mod(L,2) == 1) then
         xc = - xc
      end if
      !
      if (abs(xc) .lt. eps10) then
         return
      end if
      !
      wa = slater_integral_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d,.true.)
      wb = slater_integral_explicit(L,rwf_a,rwf_b,rwf_c,rwf_d)
      wc = debye_integral_explicit(L,lambda,rwf_a,rwf_b,rwf_c,rwf_d)
      !
      print *, "wa, wb, wc = ",wa, wb, wc
      !
      XL_Debye = xc * debye_integral_explicit(L,lambda,rwf_a,rwf_b,rwf_c,rwf_d)
      !                                                                 
   end function XL_Debye_strength_grasp2k
   !
   !
   subroutine XL_Gaunt_coefficients(L,a,b,c,d,xlcoeff,nx,nxmax)
   !--------------------------------------------------------------------
   ! Determines and calculates the number, weights, and identifiers of 
   ! the basic W_mu^nu (...) integral arrays which contribute to effective
   ! Gaunt interaction strengths X^L_Gaunt (abcd).
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)      :: L, nxmax
      integer, intent(out)     :: nx
      type(nkappa), intent(in) :: a, b, c, d
      type(XLcoefficient), dimension(1:nxmax), intent(out) :: xlcoeff
      !
      integer       :: ja, jb, jc, jd, la, lb, lc, ld, nu
      real(kind=dp) :: wa, xc
      !
      ! Initializes the number of XL_Gaunt coefficients
      nx = 0
      !
      ja = angular_momentum_j(a%kappa);  jb = angular_momentum_j(b%kappa) 
      jc = angular_momentum_j(c%kappa);  jd = angular_momentum_j(d%kappa) 
      !
      la = angular_momentum_l(a%kappa);  lb = angular_momentum_l(b%kappa) 
      lc = angular_momentum_l(c%kappa);  ld = angular_momentum_l(d%kappa) 
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0) then
         return
      end if 
      !
      xc = CL_reduced_me(a%kappa,L,c%kappa) * CL_reduced_me(b%kappa,L,d%kappa)
      if (mod(L,2) == 1) then
         !! xc = - xc
      end if
      !
      if (abs(xc) .lt. eps10) then
         return
      end if
      !
      ! Consider nu = L-1, nu = L, and nu = L+1 in turn
      ! 
      do  nu = L-1,L+1
         if (mod(la+lc+nu+4,2) /= 1   .or.   &
             mod(lb+ld+nu+4,2) /= 1   .or.   nu < 0) then
            cycle
         end if
         !
         wa = - angular_E_beta_nu( 1,nu,a%kappa,c%kappa,L) * &
                angular_E_beta_nu( 1,nu,b%kappa,d%kappa,L) * (-1)**nu
         if (abs(wa*xc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = wa*xc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = wa*xc
         end if
         !
         wa = angular_E_beta_nu( 1,nu,a%kappa,c%kappa,L) * &
              angular_E_beta_nu(-1,nu,b%kappa,d%kappa,L) * (-1)**nu
         if (abs(wa*xc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = c;   xlcoeff(nx)%c  = a
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = wa*xc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = c;   xlcoeff(nx)%d  = a
            xlcoeff(nx)%coeff = wa*xc
         end if
         !
         wa = angular_E_beta_nu(-1,nu,a%kappa,c%kappa,L) * &
              angular_E_beta_nu( 1,nu,b%kappa,d%kappa,L) * (-1)**nu
         if (abs(wa*xc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = d;   xlcoeff(nx)%d  = b
            xlcoeff(nx)%coeff = wa*xc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = d;   xlcoeff(nx)%c  = b
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = wa*xc
         end if
         !
         wa = - angular_E_beta_nu(-1,nu,a%kappa,c%kappa,L) * &
                angular_E_beta_nu(-1,nu,b%kappa,d%kappa,L) * (-1)**nu
         if (abs(wa*xc) > eps10) then
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = a;   xlcoeff(nx)%c  = c
            xlcoeff(nx)%b  = b;   xlcoeff(nx)%d  = d
            xlcoeff(nx)%coeff = wa*xc
            !
            nx = nx + 1
            xlcoeff(nx)%mu = 5;   xlcoeff(nx)%nu = nu  
            xlcoeff(nx)%a  = b;   xlcoeff(nx)%c  = d
            xlcoeff(nx)%b  = a;   xlcoeff(nx)%d  = c
            xlcoeff(nx)%coeff = wa*xc
         end if
      end do
      !
      if (nx > nxmax) then 
         stop "XL_Gaunt_coeffients(): program stop A."
      end if
      !
   end subroutine XL_Gaunt_coefficients
   !
   !
   function XL_Gaunt_strength(L,a,b,c,d)                result(XL_Gaunt)
   !--------------------------------------------------------------------
   ! Calculates the effective Gaunt interaction strengths 
   ! X^L_Gaunt (abcd) for given rank L and orbital quantum numbers
   ! a, b, c, and d.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)      :: L
      type(nkappa), intent(in) :: a, b, c, d
      real(kind=dp)            :: XL_Gaunt
      !
      integer :: i, ja, jb, jc, jd, nx, nxmax
      type(XLcoefficient), dimension(:), pointer :: xlcoeff
      !
      XL_Gaunt = zero
      !
      ja = angular_momentum_j(a%kappa);  jb = angular_momentum_j(b%kappa) 
      jc = angular_momentum_j(c%kappa);  jd = angular_momentum_j(d%kappa) 
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0   .or.  &
          L == 0) then
         return
      end if 
      !
      ! Allocate an array for calculating coefficients and calculate them
      nxmax = 26;   nx = 0
      allocate( xlcoeff(1:nxmax) )
      call XL_Gaunt_coefficients(L,a,b,c,d,xlcoeff,nx,nxmax)
      !
      if (rabs_use_basis_set) then
      do  i = 1,nx
         XL_Gaunt = XL_Gaunt + xlcoeff(i)%coeff *         &
                    W_integral(xlcoeff(i)%mu,xlcoeff(i)%nu, &
                               xlcoeff(i)%a,xlcoeff(i)%c,   &
                               xlcoeff(i)%b,xlcoeff(i)%d,.true.)
         !                      xlcoeff(i)%b,xlcoeff(i)%d,silent=.true.)
      end do
      end if
      !
      deallocate( xlcoeff )
      !
   end function XL_Gaunt_strength
   !
   !
   function XL_Gaunt_strength_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d) &
                                                        result(XL_Gaunt)
   !--------------------------------------------------------------------
   ! Calculates the effective Gaunt interaction strengths 
   ! X^L_Gaunt (abcd) for given rank L and orbital quantum numbers
   ! a, b, c, and d using GRASP92 wave functions.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)                :: L
      type(orbital_function), intent(in), target :: rwf_a, rwf_b, rwf_c, rwf_d
      real(kind=dp)                      :: XL_Gaunt
      !
      integer       :: i, ja, jb, jc, jd, nx, nxmax
      real(kind=dp) :: wa
      type(orbital_function), pointer    :: rwf_aa, rwf_bb, rwf_cc, rwf_dd
      type(XLcoefficient), dimension(:), pointer :: xlcoeff
      !
      XL_Gaunt = zero
      !
      ja = angular_momentum_j(rwf_a%orbital%kappa)
      jb = angular_momentum_j(rwf_b%orbital%kappa) 
      jc = angular_momentum_j(rwf_c%orbital%kappa)
      jd = angular_momentum_j(rwf_d%orbital%kappa) 
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0   .or.  &
          L == 0) then
         return
      end if 
      !
      ! Allocate an array for calculating coefficients and calculate them
      nxmax = 26;   nx = 0
      allocate( xlcoeff(1:nxmax) )
      call XL_Gaunt_coefficients(L,rwf_a%orbital,rwf_b%orbital,rwf_c%orbital, &
                                                 rwf_d%orbital,xlcoeff,nx,nxmax)
      !
      do  i = 1,nx
         if (xlcoeff(i)%a%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%a%kappa == rwf_a%orbital%kappa) then
            rwf_aa => rwf_a
         else if (xlcoeff(i)%a%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%a%kappa == rwf_b%orbital%kappa) then
            rwf_aa => rwf_b
         else if (xlcoeff(i)%a%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%a%kappa == rwf_c%orbital%kappa) then
            rwf_aa => rwf_c
         else if (xlcoeff(i)%a%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%a%kappa == rwf_d%orbital%kappa) then
            rwf_aa => rwf_d
         else
            stop "XL_Gaunt_strength_grasp2k(): program stop A."
         end if
         !
         if (xlcoeff(i)%b%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%b%kappa == rwf_a%orbital%kappa) then
            rwf_bb => rwf_a
         else if (xlcoeff(i)%b%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%b%kappa == rwf_b%orbital%kappa) then
            rwf_bb => rwf_b
         else if (xlcoeff(i)%b%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%b%kappa == rwf_c%orbital%kappa) then
            rwf_bb => rwf_c
         else if (xlcoeff(i)%b%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%b%kappa == rwf_d%orbital%kappa) then
            rwf_bb => rwf_d
         else
            stop "XL_Gaunt_strength_grasp2k(): program stop B."
         end if
         !
         if (xlcoeff(i)%c%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%c%kappa == rwf_a%orbital%kappa) then
            rwf_cc => rwf_a
         else if (xlcoeff(i)%c%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%c%kappa == rwf_b%orbital%kappa) then
            rwf_cc => rwf_b
         else if (xlcoeff(i)%c%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%c%kappa == rwf_c%orbital%kappa) then
            rwf_cc => rwf_c
         else if (xlcoeff(i)%c%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%c%kappa == rwf_d%orbital%kappa) then
            rwf_cc => rwf_d
         else
            stop "XL_Gaunt_strength_grasp2k(): program stop C."
         end if
         !
         if (xlcoeff(i)%d%n == rwf_a%orbital%n  .and.      &
             xlcoeff(i)%d%kappa == rwf_a%orbital%kappa) then
            rwf_dd => rwf_a
         else if (xlcoeff(i)%d%n == rwf_b%orbital%n  .and. &
             xlcoeff(i)%d%kappa == rwf_b%orbital%kappa) then
            rwf_dd => rwf_b
         else if (xlcoeff(i)%d%n == rwf_c%orbital%n  .and. &
             xlcoeff(i)%d%kappa == rwf_c%orbital%kappa) then
            rwf_dd => rwf_c
         else if (xlcoeff(i)%d%n == rwf_d%orbital%n  .and. &
             xlcoeff(i)%d%kappa == rwf_d%orbital%kappa) then
            rwf_dd => rwf_d
         else
            stop "XL_Gaunt_strength_grasp2k(): program stop A."
         end if
         !
         wa = W_integral_grasp2k(xlcoeff(i)%mu,xlcoeff(i)%nu,             &
                                               rwf_aa,rwf_cc,rwf_bb,rwf_dd)
         !
         XL_Gaunt = XL_Gaunt + xlcoeff(i)%coeff * wa
      end do 
      !
      deallocate( xlcoeff )
      !
   end function XL_Gaunt_strength_grasp2k
   !
   !
   function XL_sms_strength_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d) result(X0_sms)
   !--------------------------------------------------------------------
   ! Calculates the effective specific-mass shift interaction strengths 
   ! X^0_sms (abcd) for given orbital functions rwf_a, rwf_b, rwf_c, and 
   ! rwf_d using wave functions from GRASP92.
   !--------------------------------------------------------------------
      !
      implicit none
      integer, intent(in)                :: L
      type(orbital_function), intent(in) :: rwf_a, rwf_b, rwf_c, rwf_d
      real(kind=dp)                      :: X0_sms
      !
      integer       :: ja, jb, jc, jd, la, lb, lc, ld
      real(kind=dp) :: xc
      !
      X0_sms = zero
      !
      if (L /= 1) return
      !
      ja = angular_momentum_j(rwf_a%orbital%kappa)
      jb = angular_momentum_j(rwf_b%orbital%kappa) 
      jc = angular_momentum_j(rwf_c%orbital%kappa)
      jd = angular_momentum_j(rwf_d%orbital%kappa) 
      !
      la = angular_momentum_l(rwf_a%orbital%kappa)
      lb = angular_momentum_l(rwf_b%orbital%kappa)
      lc = angular_momentum_l(rwf_c%orbital%kappa)
      ld = angular_momentum_l(rwf_d%orbital%kappa)
      !
      if (triangle(ja+1,jc+1,L+L+1) * triangle(jb+1,jd+1,L+L+1) == 0  .or. &
          mod(la+lc+L,2) == 1   .or.   mod(lb+ld+L,2) == 1) then
         return
      end if 
      !
      xc = CL_reduced_me(rwf_a%orbital%kappa,L,rwf_c%orbital%kappa) *      &
           CL_reduced_me(rwf_b%orbital%kappa,L,rwf_d%orbital%kappa)
      if (mod(L,2) == 1) then
         xc = - xc
      end if
      !
      if (abs(xc) .lt. eps10) then
         return
      end if
      !
      X0_sms = xc * V_ab_grasp2k(rwf_a,rwf_c) * V_ab_grasp2k(rwf_b,rwf_d)
      !                                                                 
   end function XL_sms_strength_grasp2k
   !
   !
   function XL_strength(mode,L,a,b,c,d)                       result(XL)
   !--------------------------------------------------------------------
   ! Calculates the effective interaction strengths X^L (abcd) for 
   ! a given type of 'interaction' and for given rank L and orbital 
   ! quantum numbers a, b, c, and d.
   ! Allowed modes are mode = { "Coulomb", "Breit0" , "Gaunt", 
   ! "transverseBreit", "Coulomb+Breit0", "Coulomb+Gaunt", 
   ! "Coulomb+transverseBreit". }
   !
   !--------------------------------------------------------------------
      !
      implicit none
      character(len=*), intent(in) :: mode
      integer, intent(in)          :: L
      type(nkappa), intent(in)     :: a, b, c, d
      real(kind=dp)                :: XL
      !
      if (mode == "Coulomb") then
         XL = XL_Coulomb_strength(L,a,b,c,d) 
      else if (mode == "Breit0") then  
         XL = XL_Breit0_strength(L,a,b,c,d) 
      else if (mode == "Gaunt") then  
         XL = XL_Gaunt_strength(L,a,b,c,d) 
      else if (mode == "CB-complete") then  
         XL = XL_Coulomb_strength(L,a,b,c,d) +  XL_Breit0_strength(L,a,b,c,d)
      else
         stop "XL_strength: program stop A."
      end if
      !
   end function XL_strength
   !
   !
   function XL_strength_grasp2k(mode,L,rwf_a,rwf_b,rwf_c,rwf_d) result(XL)
   !--------------------------------------------------------------------
   ! Calculates the effective interaction strengths X^L (abcd) for 
   ! a given type of 'interaction' and for given rank L and orbital 
   ! quantum numbers a, b, c, and d.
   ! Allowed modes are mode = { "Coulomb", "Breit0" , "Gaunt", 
   ! "transverseBreit", "Coulomb+Breit0", "Coulomb+Gaunt", 
   ! "Coulomb+transverseBreit". }
   !
   !--------------------------------------------------------------------
      !
      implicit none
      character(len=*), intent(in)       :: mode
      integer, intent(in)                :: L
      type(orbital_function), intent(in) :: rwf_a, rwf_b, rwf_c, rwf_d
      real(kind=dp)                      :: XL
      !
      if (mode == "Coulomb") then
         XL = XL_Coulomb_strength_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d) 
      else if (mode == "Breit0") then  
         XL = XL_Breit0_strength_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d) 
      else if (mode == "Gaunt") then  
         XL = XL_Gaunt_strength_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d) 
      else if (mode == "CB-complete") then  
         XL = XL_Coulomb_strength_grasp2k(L,rwf_a,rwf_b,rwf_c,rwf_d) +  &
	      XL_Breit0_strength_grasp2k( L,rwf_a,rwf_b,rwf_c,rwf_d)
      else
         stop "XL_strength_grasp2k: program stop A."
      end if
      !
   end function XL_strength_grasp2k
   !
end module rabs_XL
