module rabs_multipole
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains the procedures which are specific to the computation
! of multipole transitions and reduced matrix elements. It is used for
! calculating transition probabilites, photoionization and 
! photo-recombination amplitudes and elsewhere.
!-----------------------------------------------------------------------
   !
   use rabs_angular
   use rabs_constant
   use rabs_grasp2k
   use rabs_nucleus
   implicit none
   private
   !
   public :: multipole_calculate_Bessel
                 ! Calculates the Bessel function over the radial grid
                 ! for a given factor.
   public :: multipole_convert_probability
                 ! Returns the oscillator strength, the Einstein A and B 
                 ! coefficients, and the decay width for a given transition 
                 ! amplitude.
   public :: multipole_grant_Ipm
                 ! Calculates the I^(+,-) integral as defined by Grant.
   public :: multipole_grant_J
                 ! Calculates the J_L integral as defined by Grant.
   public :: multipole_reduced_F1_integral
                 ! Calculates a M-F1 integral as it arises due to the 
		 ! Coulomb excitation of an electron.
   public :: multipole_reduced_F23_integral
                 ! Calculates a M-F23 integral as it arises due to the 
		 ! Coulomb excitation of an electron.
   public :: multipole_reduced_M_andreys
                 ! Calculates a M integral for the coupling of the radiation
                 ! field for a given multipolarity and gauge form for the
                 ! capture of an electron due to Andrey's formula.
   public :: multipole_reduced_M_capture
                 ! Calculates a M integral for the coupling of the radiation
                 ! field for a given multipolarity and gauge form for the
                 ! capture of an electron.
   public :: multipole_reduced_M_integral
                 ! Calculates a M integral for the coupling of the radiation
                 ! field for a given multipolarity and gauge form due to
                 ! Grant (1989).
   public :: multipole_reduced_M_multiphoton
                 ! Calculates a M integral for the coupling of the radiation
                 ! field for a given multipolarity and gauge form due to
                 ! Andrey (2009).
   public :: multipole_reduced_M_operator
                 ! Calculates a M operator for the coupling of the radiation
                 ! field for a given multipolarity and gauge form for the
                 ! generation of modified orbitals.
   public :: multipole_mphoton_select
                 ! Determines the number and type of the multipole transition
                 ! components which contribute to a given multi-photon
		 ! transition.
   public :: multipole_msymmetry_select
                 ! Determines the total angular momenta J_nu and parity_nu of 
		 ! allowed intermediate states for a given pair of multipoles.
   public :: multipole_select
                 ! Determines the number and type of the multipole transition
                 ! components which contribute to a given transition.
   !
   !
   integer, public                         :: number_of_multipoles = 0
   character(len=2), dimension(20), public :: multipoles           = "  "
   !
   ! Define global logical flags for the control of the PHOTO program; the
   ! default values for these flags may be overwritten interactively during 
   ! input time
   logical, public :: multipole_print_AB_in_hartree       = .false.
   !
contains
   !
   subroutine multipole_calculate_Bessel(L,omega_over_c,bessel)
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
         bessel(i) = spherical_Bessel_jL(L,omega_over_c*r_grasp2k(i))
      end do
      !
   end subroutine multipole_calculate_bessel
   !
   !
   subroutine multipole_convert_probability(totalJ_i,totalJ_f,L,energy_au, &
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
      if (.not.multipole_print_AB_in_hartree) then
         einstein_A = einstein_A  * convert_einstein_a_to_si
         einstein_B = einstein_B  * convert_einstein_b_to_si
      end if
      !
   end subroutine multipole_convert_probability
   !
   !
   function multipole_grant_Ipm(pm,energy,L,rwf_f,rwf_i)   result(Ipm)
   !--------------------------------------------------------------------
   ! Returns the I^(+,-) integral as defined by Grant.
   !--------------------------------------------------------------------
      !
      character(len=1), intent(in)       :: pm
      integer, intent(in)                :: L
      real(kind=dp), intent(in)          :: energy
      type(orbital_function), intent(in) :: rwf_f, rwf_i
      real(kind=dp)                      :: Ipm, arg
      !
      integer                                  :: mtp
      real(kind=dp), dimension(:), allocatable :: ta
      real(kind=dp), dimension(1:n_grasp2k)    :: bessel
      !
      mtp = min( rwf_f%mtp, rwf_i%mtp)
      allocate( ta(1:mtp+10) )
      !
      if (pm == "+") then
         ta(1)     = zero
         ta(2:mtp) = (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +	  &
                      rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) ) * rp_grasp2k(2:mtp)   
         !		    
         arg  = energy / c         
         call multipole_calculate_bessel(L,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         Ipm = quad_grasp2k(ta,mtp)
      else if (pm == "-") then
         ta(1)     = zero
         ta(2:mtp) = (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -	  &
                      rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) ) * rp_grasp2k(2:mtp)   
         !		    
         arg  = energy / c         
         call multipole_calculate_bessel(L,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         Ipm = quad_grasp2k(ta,mtp)
      else
         Ipm = 0._dp
      end if
      !
      deallocate( ta ) 
      !
   end function multipole_grant_Ipm
   !
   !
   function multipole_grant_J(energy,L,rwf_f,rwf_i)         result(JL)
   !--------------------------------------------------------------------
   ! Returns the J_L integral as defined by Grant.
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: L
      real(kind=dp), intent(in)          :: energy
      type(orbital_function), intent(in) :: rwf_f, rwf_i
      real(kind=dp)                      :: JL, arg
      !
      integer                                  :: mtp
      real(kind=dp), dimension(:), allocatable :: ta
      real(kind=dp), dimension(1:n_grasp2k)    :: bessel
      !
      mtp = min( rwf_f%mtp, rwf_i%mtp)
      allocate( ta(1:mtp+10) )
      !
      ta(1)	= zero
      ta(2:mtp) = (rwf_f%P(2:mtp)*rwf_i%P(2:mtp) +     &
        	   rwf_f%Q(2:mtp)*rwf_i%Q(2:mtp) ) * rp_grasp2k(2:mtp)   
      ! 		 
      arg  = energy / c 	
      call multipole_calculate_bessel(L,arg,bessel)
      ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
      !
      JL = quad_grasp2k(ta,mtp)
      !
      deallocate( ta ) 
      !
   end function multipole_grant_J
   !
   ! 
   function multipole_reduced_F1_integral(q,L,p,rwf_f,rwf_i)   result(red_me)
   !--------------------------------------------------------------------
   ! Calculates the value of an (one-electron) M_F1 integral as it arises
   ! in the Coulomb excitation of electrons in the field of the target 
   ! nucleus. It calculates all of the required Bessel functions; 
   ! the radial wave functions are provided by the derived data structures 
   ! rwf_f and rwf_i.
   !
   ! Calls: angular_momentum_j(), orbital_name(), quad_grasp2k(),
   !        wigner_3j_symbol().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: L, p
      real(kind=dp), intent(in)          :: q
      type(orbital_function), intent(in) :: rwf_f, rwf_i
      real(kind=dp)                      :: red_me
      !
      integer       :: ja, jb, kapa, kapb, mtp
      !
      real(kind=dp), dimension(:), allocatable :: ta
      real(kind=dp), dimension(1:n_grasp2k)    :: bessel
      !
      kapa = rwf_f%orbital%kappa;        kapb = rwf_i%orbital%kappa
      ja   = angular_momentum_j(kapa);   jb   = angular_momentum_j(kapb)
      !
      mtp = min( rwf_f%mtp, rwf_i%mtp)
      allocate( ta(1:mtp+10) )
      !
      ta(1)     = zero
      ta(2:mtp) = (rwf_f%P(2:mtp)*rwf_i%P(2:mtp) +	  &
                   rwf_f%Q(2:mtp)*rwf_i%Q(2:mtp) ) * rp_grasp2k(2:mtp)   
      !	
      call multipole_calculate_bessel(L,q,bessel)
      ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
      !
      red_me = sqrt( (L+L+1) / (two*two*pi) ) * CL_reduced_me(kapa,L,kapb)  &
               * quad_grasp2k(ta,mtp)    
      red_me = red_me / sqrt(ja + one)
      !    
      deallocate( ta ) 
      !
   end function multipole_reduced_F1_integral
   !
   ! 
   function multipole_reduced_F23_integral(q,L,t,p,rwf_f,rwf_i) result(red_me)
   !--------------------------------------------------------------------
   ! Calculates the value of an (one-electron) M_F23 integral as it arises
   ! in the Coulomb excitation of electrons in the field of the target 
   ! nucleus. It calculates all of the required Bessel functions; 
   ! the radial wave functions are provided by the derived data structures 
   ! rwf_f and rwf_i.
   !
   ! Calls: angular_momentum_j(), orbital_name(), quad_grasp2k(),
   !        wigner_3j_symbol().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: L, t, p
      real(kind=dp), intent(in)          :: q
      type(orbital_function), intent(in) :: rwf_f, rwf_i
      real(kind=dp)                      :: red_me
      !
      integer       :: ja, jb, kapa, kapb, mtp
      !
      real(kind=dp), dimension(:), allocatable :: ta
      real(kind=dp), dimension(1:n_grasp2k)    :: bessel
      !
      red_me = zero
      !
      kapa = rwf_f%orbital%kappa;        kapb = rwf_i%orbital%kappa
      ja   = angular_momentum_j(kapa);   jb   = angular_momentum_j(kapb)
      !
      mtp = min( rwf_f%mtp, rwf_i%mtp)
      allocate( ta(1:mtp+10) )
      !
      ta(1)     = zero
      ta(2:mtp) = (rwf_f%Q(2:mtp)*rwf_i%P(2:mtp)) * rp_grasp2k(2:mtp)   
      !		    
      call multipole_calculate_bessel(L,q,bessel)
      ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
      !
      red_me = red_me + sigma_TtL_reduced_me(-kapa,L,t,kapb)           &
                        * quad_grasp2k(ta,mtp)        
      !
      ta(1)     = zero
      ta(2:mtp) = (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp)) * rp_grasp2k(2:mtp)   
      !		    
      call multipole_calculate_bessel(L,q,bessel)
      ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
      !
      red_me = red_me - sigma_TtL_reduced_me(kapa,L,t,-kapb)           &
                        * quad_grasp2k(ta,mtp)        
      red_me = red_me / sqrt(ja + one)
      !
      deallocate( ta ) 
      !
   end function multipole_reduced_F23_integral
   !
   ! 
   function multipole_reduced_M_andreys(energy,L,p,gauge,rwf_f,rwf_i) &
                                                          result(red_me)
   !--------------------------------------------------------------------
   ! Calculates the value of an (one-electron) M integral for the coupling
   ! of electrons with the radiation field for a given multipole and
   ! gauge, following Andrey's one-particle formula. It calculates all of 
   ! the required Bessel functions; the radial wave functions are provided 
   ! by the derived data structures rwf_f and rwf_i.
   !
   ! Calls: angular_momentum_j(), orbital_name(), quad_grasp2k(),
   !        wigner_3j_symbol().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: L, p
      real(kind=dp), intent(in)          :: energy
      character(len=9), intent(in)       :: gauge
      type(orbital_function), intent(in) :: rwf_f, rwf_i
      real(kind=dp)                      :: red_me
      !
      integer       :: ja, jb, kapa, kapb, mtp, phase
      real(kind=dp) :: arg, factor, red_Coulomb, red_Gauge, wa, wb
      !
      kapa = rwf_f%orbital%kappa;        kapb = rwf_i%orbital%kappa
      ja   = angular_momentum_j(kapa);   jb   = angular_momentum_j(kapb)
      !
      red_me = zero
      !
      ! Test parity rules
      if (p == 0) then
         ! Magnetic case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 0)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 2) ) then
         else
            return
         end if
      else
         ! Electric case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 2)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 0) ) then
         else
            return
         end if
      end if
      !
      factor = ((-1)**L) * sqrt(ja+one) * Clebsch_Gordan(ja,1,L+L,0,jb,1) &
                / (- sqrt(two*two*pi))
      !
      if (p == 0) then
         ! Magnetic case
	 red_me = factor * (kapa+kapb) * sqrt((L+L+one)/(L*(L+one))) &
	          * multipole_grant_Ipm("+",energy,L,rwf_f,rwf_i)
      else if (p == 1) then
         ! Electric case
	 red_me = sqrt( L/((L+one)*(L+L+one)) )                                &
	        * ( (kapa-kapb) *                                              &
		     multipole_grant_Ipm("+",energy,L+1,rwf_f,rwf_i)           &
		  + (L+one) * multipole_grant_Ipm("-",energy,L+1,rwf_f,rwf_i) )&
              -   sqrt( (L+one)/(L*(L+L+one)) )                                &
	        * ( (kapa-kapb) *                                              &
		     multipole_grant_Ipm("+",energy,L-1,rwf_f,rwf_i)           &
		  -  L * multipole_grant_Ipm("-",energy,L-1,rwf_f,rwf_i) )
	 !
	 red_me = factor * red_me   
      else
         stop "multipole_reduced_M_andreys(): program stop A."
      end if
      !
   end function multipole_reduced_M_andreys
   !
   ! 
   function multipole_reduced_M_capture(energy,L,p,gauge,rwf_f,rwf_i) &
                                                          result(red_me)
   !--------------------------------------------------------------------
   ! Calculates the value of an (one-electron) M integral for the coupling
   ! of electrons with the radiation field for a given multipole and
   ! gauge. It calculates all of the required Bessel functions; 
   ! the radial wave functions are provided by the derived data structures 
   ! rwf_f and rwf_i.
   !
   ! This procedure follows our conventions with Andrey for calculating
   ! capture cross sections.
   !
   ! Calls: angular_momentum_j(), orbital_name(), quad_grasp2k(),
   !        wigner_3j_symbol().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: L, p
      real(kind=dp), intent(in)          :: energy
      character(len=9), intent(in)       :: gauge
      type(orbital_function), intent(in) :: rwf_f, rwf_i
      real(kind=dp)                      :: red_me
      !
      integer       :: ja, jb, kapa, kapb, mtp, phase
      real(kind=dp) :: arg, factor, red_Coulomb, red_Gauge, wa, wb
      real(kind=dp), dimension(:), allocatable :: ta
      real(kind=dp), dimension(1:n_grasp2k)    :: bessel
      !
      kapa = rwf_f%orbital%kappa;        kapb = rwf_i%orbital%kappa
      ja   = angular_momentum_j(kapa);   jb   = angular_momentum_j(kapb)
      !
      red_me = zero
      !
      ! Test parity rules
      if (p == 0) then
         ! Magnetic case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 0)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 2) ) then
         else
            return
         end if
      else
         ! Electric case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 2)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 0) ) then
         else
            return
         end if
      end if
      !
      !
      ! Evaluate factor multiplying Mbar(a,b); use one-particle matrix elements
      ! of Brink-Satchler type
      !
      factor = ((-one)**L) *sqrt(ja+one)*Clebsch_Gordan(ja,1,L+L,0,jb,1) &
               / (- sqrt(two*two*pi))
      !
      if (abs(factor) < eps10*eps10) then
         red_me = zero
         return
      end if
      !
      mtp = min( rwf_f%mtp, rwf_i%mtp)
      allocate( ta(1:mtp+10) )
      !
      if (p == 0) then
         ! Magnetic case
         ta(1) = zero
         ta(2:mtp) = ( rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +                   &
                       rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) ) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_Bessel(L,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         !!x red_me = -(L+L+one)*(kapa+kapb)/sqrt(L*(L+one)) * quad_grasp2k(ta,mtp)
         red_me = (kapa+kapb)*sqrt((L+L+one)/(L*(L+one)))*quad_grasp2k(ta,mtp)
         red_me = factor * red_me
      else if (p == 1) then
         ! Electric case
         ta(1) = zero
         !
         ! Tabulate coulomb integrand
         wa    = sqrt(L/((L+one)*(L+L+one)))*(kapa-kapb)
         wb    = sqrt(L/((L+one)*(L+L+one)))*(L+one)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_bessel(L+1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Coulomb = factor * quad_grasp2k(ta,mtp)
         !
         wa    = - sqrt((L+one)/(L*(L+L+one)))*(kapa-kapb)
         wb    =   sqrt((L+one)/(L*(L+L+one)))*L
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_bessel(L-1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Coulomb = red_Coulomb + factor * quad_grasp2k(ta,mtp)
         !
         if (gauge == "Coulomb  ") then
            red_me = red_Coulomb
	    !!x print *, "*** gauge,red_me = ",gauge,red_me 
            if (.false.) &
               write(99,*) "factor, red_Coulomb = ",factor, red_me
            deallocate( ta )
            !!x print *, "c: red_me = ",red_me
            return
         else if (gauge == "Babushkin") then
         else if (rabs_use_stop) then
            stop "multipole_reduced_M_capture(): program stop A."
         end if
         !
         wa    = (kapa-kapb);  wb = (L+1)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_bessel(L+1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Gauge = factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         !
         wa    = (kapa-kapb);  wb = -L
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_bessel(L-1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Gauge = red_Gauge + factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         !
         wa    = -(L+L+one)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%P(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%Q(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_bessel(L,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Gauge = red_Gauge + factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         red_me    = red_Coulomb + red_Gauge * sqrt(one/(L+L+one))
         !
      else
         stop "multipole_reduced_M_capture(): program stop B."
      end if
      !
      deallocate( ta ) 
      !
   end function multipole_reduced_M_capture
   !
   ! 
   function multipole_reduced_M_integral(energy,L,p,gauge,rwf_f,rwf_i) &
                                                          result(red_me)
   !--------------------------------------------------------------------
   ! Calculates the value of an (one-electron) M integral for the coupling
   ! of electrons with the radiation field for given multipole and
   ! gauge. It applies the set of Bessel functions which are associated
   ! with the transition tr; the radial wave functions are provided by the
   ! derived data structures rwf_f and rwf_i.
   !
   ! This procedure follows Grant's conventions which is used also in REOS
   ! and EINSTEIN.
   !
   ! Calls: angular_momentum_j(), orbital_name(), quad_grasp2k(),
   !        wigner_3j_symbol().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: L, p
      real(kind=dp), intent(in)          :: energy
      character(len=9), intent(in)       :: gauge
      type(orbital_function), intent(in) :: rwf_f, rwf_i
      real(kind=dp)                      :: red_me
      !
      integer       :: ja, jb, kapa, kapb, mtp, phase
      real(kind=dp) :: arg, factor, red_Coulomb, red_Gauge, wa, wb
      real(kind=dp), dimension(:), allocatable :: ta
      real(kind=dp), dimension(1:n_grasp2k)    :: bessel
      !
      red_me = zero
      !
      kapa = rwf_f%orbital%kappa;        kapb = rwf_i%orbital%kappa
      ja   = angular_momentum_j(kapa);   jb   = angular_momentum_j(kapb)
      !
      ! Test parity rules
      if (p == 0) then
         ! Magnetic case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 0)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 2) ) then
         else
            return
         end if
      else
         ! Electric case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 2)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 0) ) then
         else
            return
         end if
      end if
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
      if (p == 0) then
         ! Magnetic case
         ta(1) = zero
         ta(2:mtp) = ( rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +                   &
                       rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) ) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_Bessel(L,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_me = -(L+L+one)*(kapa+kapb)/sqrt(L*(L+one)) * quad_grasp2k(ta,mtp)
         red_me = factor * red_me
      else if (p == 1) then
         ! Electric case
         ta(1) = zero
         !
         ! Tabulate coulomb integrand
         wa    = sqrt(L/(L+one))*(kapa-kapb);  wb = sqrt(L/(L+one))*(L+one)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_bessel(L+1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Coulomb = factor * quad_grasp2k(ta,mtp)
         !
         wa    = - sqrt((L+one)/L)*(kapa-kapb);  wb = sqrt((L+one)/L)*L
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_bessel(L-1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Coulomb = red_Coulomb + factor * quad_grasp2k(ta,mtp)
         !
         if (gauge == "Coulomb  ") then
            red_me = red_Coulomb
            return
         else if (gauge == "Babushkin") then
         else if (rabs_use_stop) then
            stop "multipole_reduced_M_integral(): program stop A."
         end if
         !
         wa    = (kapa-kapb);  wb = (L+1)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_bessel(L+1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Gauge = factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         !
         wa    = (kapa-kapb);  wb = -L
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )        &
                   + wb * (rwf_f%P(2:mtp)*rwf_i%Q(2:mtp) -        &
                           rwf_f%Q(2:mtp)*rwf_i%P(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_bessel(L-1,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Gauge = red_Gauge + factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         !
         wa    = -(L+L+one)
         ta(2:mtp) = wa * (rwf_f%P(2:mtp)*rwf_i%P(2:mtp) +        &
                           rwf_f%Q(2:mtp)*rwf_i%Q(2:mtp) )
         ta(2:mtp) = ta(2:mtp) * rp_grasp2k(2:mtp)
         !
         arg  = energy / c         
         call multipole_calculate_bessel(L,arg,bessel)
         ta(2:mtp) = ta(2:mtp) * bessel(2:mtp)
         !
         red_Gauge = red_Gauge + factor * sqrt((L+one)/L) * quad_grasp2k(ta,mtp)
         red_me    = red_Coulomb + red_Gauge
         !
      else
         stop "multipole_reduced_M_integral(): program stop B."
      end if
      !
      deallocate( ta ) 
      !
   end function multipole_reduced_M_integral
   !
   ! 
   function multipole_reduced_M_multiphoton(energy,L,p,gauge,rwf_f,rwf_i) &
                                                             result(red_me)
   !--------------------------------------------------------------------
   ! Calculates the value of an (one-electron) M integral for the coupling
   ! of electrons with the radiation field for a given multipole and
   ! gauge, following Andrey's one-particle formula for the multiphoton
   ! decay. It calculates all of the required Bessel functions; the radial 
   ! wave functions are provided by the derived data structures rwf_f 
   ! and rwf_i.
   !
   ! Calls: angular_momentum_j(), orbital_name(), quad_grasp2k(),
   !        wigner_3j_symbol().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: L, p
      real(kind=dp), intent(in)          :: energy
      character(len=9), intent(in)       :: gauge
      type(orbital_function), intent(in) :: rwf_f, rwf_i
      complex(kind=dp)                   :: red_me
      !
      integer       :: ja, jb, kapa, kapb, mtp, phase, na, nb
      real(kind=dp) :: arg, factor, G, red_Coulomb, red_Gauge, wb, wbp, wbm
      complex(kind=dp)                   :: wa, wg
      !
      kapa = rwf_f%orbital%kappa;        kapb = rwf_i%orbital%kappa
      ja   = angular_momentum_j(kapa);   jb   = angular_momentum_j(kapb)
      !
      red_me = zero
      !
      ! Test parity rules
      if (p == 0) then
         ! Magnetic case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 0)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 2) ) then
         else
            return
         end if
      else
         ! Electric case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 2)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 0) ) then
         else
            return
         end if
      end if
      !
      factor = ((-1)**L) * sqrt(ja+one) * Clebsch_Gordan(ja,1,L+L,0,jb,1) 
      !
      if (p == 0) then
         ! Magnetic case
	 red_me = cmplx(zero, -one) / sqrt(two*two*pi)                         &
	          * (kapa+kapb) * sqrt((L+L+one)/(L*(L+one)))                  &
	          * multipole_grant_Ipm("+",energy,L,rwf_f,rwf_i)              
		  !!x * CL_reduced_me(kapa,L,-kapb)
	 wa     = red_me
	 red_me = conjg(red_me * factor)
      else if (p == 1) then
         ! Electric case
	 red_me = cmplx(zero, -one) / sqrt(two*two*pi)  		    &
	     * ( sqrt( L/((L+one)*(L+L+one)) )  			    &
	     * ( (kapa-kapb) *  					    &
	     	  multipole_grant_Ipm("+",energy,L+1,rwf_f,rwf_i)	    &
	       + (L+one) * multipole_grant_Ipm("-",energy,L+1,rwf_f,rwf_i) )&
           -   sqrt( (L+one)/(L*(L+L+one)) )				    &
	     * ( (kapa-kapb) *  					    &
	     	  multipole_grant_Ipm("+",energy,L-1,rwf_f,rwf_i)	    &
	       -  L * multipole_grant_Ipm("-",energy,L-1,rwf_f,rwf_i) ) )  
	 wa	= red_me
	 wbp	= multipole_grant_Ipm("+",energy,L+1,rwf_f,rwf_i)
	 wbm	= multipole_grant_Ipm("+",energy,L-1,rwf_f,rwf_i)
	 !
	 if (gauge == "Babushkin") then
	    G = - sqrt( (L+one) / L )
	    !
	    wg = cmplx(zero, one) / sqrt(two*two*pi)  		            &
	       * ( sqrt( one/(L+L+one) ) )        			    &
	       * ( (kapa-kapb) *  					    &
	     	   (multipole_grant_Ipm("+",energy,L+1,rwf_f,rwf_i)	    &
	          + multipole_grant_Ipm("+",energy,L-1,rwf_f,rwf_i))        &
	       +  (L+one) * multipole_grant_Ipm("-",energy,L+1,rwf_f,rwf_i) &
	       -   L * multipole_grant_Ipm("-",energy,L-1,rwf_f,rwf_i)      &
	       +  (L+L+1) * multipole_grant_J(energy,L,rwf_f,rwf_i) )  
	    red_me = red_me + G*wg
	 end if
	 !
	 red_me = conjg(red_me * factor)
      else
         stop "multipole_reduced_M_multiphoton(): program stop A."
      end if
      !
      if (mod(ja-jb+p+p + 32,4) == 2) then
         red_me = - red_me
      else if (mod(ja-jb+p+p + 32,4) == 0) then
      else
         stop "multipole_reduced_M_multiphoton(): program stop A."
      end if
      !
   end function multipole_reduced_M_multiphoton
   !
   ! 
   subroutine multipole_reduced_M_operator(energy,L,p,gauge,rwf_b, &
                                           kapa,Pa,Qa,mtp) 
   !--------------------------------------------------------------------
   ! Calculates the (one-electron) M operator
   ! 
   !                  <kappa_a || M^L || n_b, kappa_b >
   ! 
   ! for the coupling of electrons with the radiation field for a given 
   ! multipole and gauge. It calculates all of the required Bessel functions; 
   ! the radial wave functions are provided by the derived data structures 
   ! rwf_b while the two components of the operator are returned in the 
   ! arrays Pa(:) and Qa(:).
   !
   ! This procedure follows our conventions with Andrey for calculating
   ! capture cross sections.
   !
   ! Calls: angular_momentum_j(), orbital_name(), quad_grasp2k(),
   !        wigner_3j_symbol().
   !--------------------------------------------------------------------
      !
      integer, intent(in)                :: L, p, kapa
      real(kind=dp), intent(in)          :: energy
      character(len=9), intent(in)       :: gauge
      type(orbital_function), intent(in) :: rwf_b
      integer, intent(out)               :: mtp
      real(kind=dp), dimension(1:n_grasp2k), intent(out) :: Pa, Qa
      !
      integer       :: ja, jb, kapb, phase
      real(kind=dp) :: arg, factor, red_Coulomb, red_Gauge, wa, wb
      !!x real(kind=dp), dimension(:), allocatable :: ta
      real(kind=dp), dimension(1:n_grasp2k)    :: bessel
      !
      kapb = rwf_b%orbital%kappa
      ja   = angular_momentum_j(kapa)
      jb   = angular_momentum_j(kapb)
      !
      Pa(:) = zero;   Qa(:) = zero
      mtp   = rwf_b%mtp
      !
      ! Test parity rules
      if (p == 0) then
         ! Magnetic case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 0)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 2) ) then
         else
            return
         end if
      else
         ! Electric case
         if ( (kapa*kapb > 0 .and. mod(ja+jb+L+L,4) == 2)   .or. &
              (kapa*kapb < 0 .and. mod(ja+jb+L+L,4) == 0) ) then
         else
            return
         end if
      end if
      !
      !
      ! Evaluate factor multiplying Mbar(a,b); use one-particle matrix elements
      ! of Brink-Satchler type
      !
      factor = ((-one)**L) *sqrt(ja+one)*Clebsch_Gordan(ja,1,L+L,0,jb,1) &
               / (- sqrt(two*two*pi))
      !
      if (abs(factor) < eps10*eps10) then
         return
      end if
      !
      if (p == 0) then
         ! Magnetic case
	 wa   = factor * (kapa+kapb)*sqrt((L+L+one)/(L*(L+one)))
         arg  = energy / c         
         call multipole_calculate_Bessel(L,arg,bessel)
	 Pa(2:mtp) = wa * rwf_b%Q(2:mtp) * rp_grasp2k(2:mtp) * bessel(2:mtp)
	 Qa(2:mtp) = wa * rwf_b%P(2:mtp) * rp_grasp2k(2:mtp) * bessel(2:mtp)
      else if (p == 1) then
         ! Electric case
         wa   = factor * sqrt(L/((L+one)*(L+L+one)))*(kapa-kapb)
         wb   = factor * sqrt(L/((L+one)*(L+L+one)))*(L+one)
         arg  = energy / c         
         call multipole_calculate_bessel(L+1,arg,bessel)
	 Pa(2:mtp) = (wa+wb) * rwf_b%Q(2:mtp) * rp_grasp2k(2:mtp) &
	                     * bessel(2:mtp) 
	 Qa(2:mtp) = (wa-wb) * rwf_b%P(2:mtp) * rp_grasp2k(2:mtp) &
	                     * bessel(2:mtp)  
         !
         wa    = - factor * sqrt((L+one)/(L*(L+L+one)))*(kapa-kapb)
         wb    =   factor * sqrt((L+one)/(L*(L+L+one)))*L
         call multipole_calculate_bessel(L-1,arg,bessel)
         Pa(2:mtp) = Pa(2:mtp) + (wa+wb) * rwf_b%Q(2:mtp) * rp_grasp2k(2:mtp)&
	                                 * bessel(2:mtp)
         Qa(2:mtp) = Qa(2:mtp) + (wa-wb) * rwf_b%P(2:mtp) * rp_grasp2k(2:mtp)&
	                                 * bessel(2:mtp)
         !
         if (gauge == "Coulomb  ") then
            return
         else if (gauge == "Babushkin") then
         else 
            stop "multipole_reduced_M_operator(): program stop A."
         end if
         !
         wa    = factor * sqrt((L+one)/L) * (kapa-kapb) * sqrt(one/(L+L+one))
         wb    = factor * sqrt((L+one)/L) * (L+1) * sqrt(one/(L+L+one))
         call multipole_calculate_bessel(L+1,arg,bessel)
         Pa(2:mtp) = Pa(2:mtp) + (wa+wb) * rwf_b%Q(2:mtp) * rp_grasp2k(2:mtp)&
	                                 * bessel(2:mtp)
         Qa(2:mtp) = Qa(2:mtp) + (wa-wb) * rwf_b%P(2:mtp) * rp_grasp2k(2:mtp)&
	                                 * bessel(2:mtp)
         !
         wa    = factor * sqrt((L+one)/L) * (kapa-kapb) * sqrt(one/(L+L+one))
         wb    = - factor * sqrt((L+one)/L) * L * sqrt(one/(L+L+one))
         call multipole_calculate_bessel(L-1,arg,bessel)
         Pa(2:mtp) = Pa(2:mtp) + (wa+wb) * rwf_b%Q(2:mtp) * rp_grasp2k(2:mtp)&
	                                 * bessel(2:mtp)
         Qa(2:mtp) = Qa(2:mtp) + (wa-wb) * rwf_b%P(2:mtp) * rp_grasp2k(2:mtp)&
	                                 * bessel(2:mtp)
         !
         wa    = - factor * sqrt((L+one)/L) * (L+L+one) * sqrt(one/(L+L+one))
         call multipole_calculate_bessel(L,arg,bessel)
         Pa(2:mtp) = Pa(2:mtp) + wa * rwf_b%Q(2:mtp) * rp_grasp2k(2:mtp)&
	                            * bessel(2:mtp)
         Qa(2:mtp) = Qa(2:mtp) + wa * rwf_b%P(2:mtp) * rp_grasp2k(2:mtp)&
	                            * bessel(2:mtp)
      else
         stop "multipole_reduced_M_operator(): program stop B."
      end if
      !
   end subroutine multipole_reduced_M_operator
   !
   !
   subroutine multipole_mphoton_select(totalJ_i,parity_i,totalJ_f,parity_f, &
                                calc_coulomb,calc_babushkin,mult,gauge,nmult)
   !--------------------------------------------------------------------
   ! Determines the number and type of the multi-photon multipole components
   ! which contribute and need to be considered for a given transition.
   ! It applies the standard selection rules and compares allowed line
   ! components with the list of multipoles selected at input time.
   !
   ! Calls: is_triangle().
   !--------------------------------------------------------------------
      !
      integer, intent(in)           :: totalJ_i, totalJ_f 
      character(len=1), intent(in)  :: parity_i, parity_f
      logical, intent(in)           :: calc_coulomb,calc_babushkin
      integer, intent(out)          :: nmult
      character(len=2), dimension(120,3), intent(out) :: mult
      character(len=9), dimension(120), intent(out) :: gauge
      !
      integer :: i, j
      !
      nmult = 0
      do  i = 1,number_of_multipoles
         do  j = 1,number_of_multipoles
	    if (multipoles(i) == "E1"  .and.  multipoles(j) == "E1") then
        	if (parity_i == parity_f   .and.           &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "E1"  .and.  multipoles(j) == "M1") then
               if (parity_i /= parity_f   .and.            &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "M1"  .and.  multipoles(j) == "E1") then
               if (parity_i /= parity_f   .and.            &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "E1"  .and.  multipoles(j) == "E2") then
               if (parity_i /= parity_f   .and.            &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "E2"  .and.  multipoles(j) == "E1") then
               if (parity_i /= parity_f   .and.            &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "E1"  .and.  multipoles(j) == "M2") then
               if (parity_i == parity_f   .and. &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "M2"  .and.  multipoles(j) == "E1") then
               if (parity_i == parity_f   .and. &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "M1"  .and.  multipoles(j) == "M1") then
               if (parity_i == parity_f   .and.            &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)) ) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Magnetic "  
               end if
	       !
	    else if (multipoles(i) == "M1"  .and.  multipoles(j) == "E2") then
               if (parity_i == parity_f   .and. &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "E2"  .and.  multipoles(j) == "M1") then
               if (parity_i == parity_f   .and. &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "M1"  .and.  multipoles(j) == "M2") then
               if (parity_i /= parity_f   .and. &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)) ) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Magnetic "  
               end if
	       !
	    else if (multipoles(i) == "M2"  .and.  multipoles(j) == "M1") then
               if (parity_i /= parity_f   .and. &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)) ) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Magnetic "  
               end if
	       !
	    else if (multipoles(i) == "E2"  .and.  multipoles(j) == "E2") then
               if (parity_i == parity_f   .and. &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)  .or. &
		    is_triangle(totalJ_i,totalJ_f,8)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "E2"  .and.  multipoles(j) == "M2") then
               if (parity_i /= parity_f   .and. &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)  .or. &
		    is_triangle(totalJ_i,totalJ_f,8)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "M2"  .and.  multipoles(j) == "E2") then
               if (parity_i /= parity_f   .and. &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)  .or. &
		    is_triangle(totalJ_i,totalJ_f,8)) ) then
		  if (calc_babushkin) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Babushkin"  
		  end if
		  !
		  if (calc_coulomb) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Coulomb  "
		  end if
               end if
	       !
	    else if (multipoles(i) == "M2"  .and.  multipoles(j) == "M2") then
               if (parity_i == parity_f   .and. &
                   (is_triangle(totalJ_i,totalJ_f,0)  .or. &
		    is_triangle(totalJ_i,totalJ_f,2)  .or. &
		    is_triangle(totalJ_i,totalJ_f,4)  .or. &
		    is_triangle(totalJ_i,totalJ_f,6)  .or. &
		    is_triangle(totalJ_i,totalJ_f,8)) ) then
                  nmult = nmult + 1
		  mult(nmult,1) = multipoles(i); mult(nmult,2) = multipoles(j)
                  gauge(nmult)  = "Magnetic "  
               end if
	       !
	    else
	       print *, "The pair ", multipoles(i), multipoles(j),          &
	                "  is not supported by the present program version. " 
            end if
	    !
	 end do
      end do
      !
      if (nmult > 120) then
         stop "multipole_mphoton_select(): program stop A."
      end if
      !
   end subroutine multipole_mphoton_select
   !
   !
   subroutine multipole_msymmetry_select(totalJ_i,parity_i,totalJ_f,parity_f, &
                                        mult1,mult2,symmetry_j,symmetry_p,njnu)
   !--------------------------------------------------------------------
   ! Determines the total angular momenta J_nu and parity_nu of those 
   ! intermediate states that can connect the initial and final states
   ! via the two multipoles as given in mult1 and mult2.
   !
   ! Calls: is_triangle().
   !--------------------------------------------------------------------
      !
      integer, intent(in)           :: totalJ_i, totalJ_f 
      character(len=1), intent(in)  :: parity_i, parity_f
      character(len=2), intent(in)  :: mult1, mult2
      integer, intent(out)          :: njnu
      character(len=1), dimension(:), intent(out) :: symmetry_p
      integer, dimension(:), intent(out)          :: symmetry_j
      !
      integer          :: j, jlow, jnu
      logical          :: add, a,b 
      character(len=1) :: jparity
      !
      njnu = 0
      !
      if (mod(totalJ_i,2) == 1)  then
         jlow = 1
      else
         jlow = 0
      end if
      if (mod(totalJ_f,2) /= jlow) then
         stop "multipole_msymmetry_select(): program stop A."
      end if
      !
      do  j = jlow,40,2
         add = .true.
	 jnu = j;   jparity = "+"
	 !
	 select case(mult1)
	 case("E1")
            if (parity_i == jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,2)) add = .false.
	 case("M1")
            if (parity_i /= jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,2)) add = .false.
	 case("E2")
            if (parity_i /= jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,4)) add = .false.
	 case("M2")
            if (parity_i == jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,4)) add = .false.
	 case("E3")
            if (parity_i == jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,6)) add = .false.
	 case("M3")
            if (parity_i /= jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,6)) add = .false.
	 case default
            stop "multipole_msymmetry_select(): program stop B."
	 end select
	 !
	 select case(mult2)
	 case("E1")
            if (parity_f == jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,2)) add = .false.
	 case("M1")
            if (parity_f /= jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,2)) add = .false.
	 case("E2")
            if (parity_f /= jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,4)) add = .false.
	 case("M2")
            if (parity_f == jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,4)) add = .false.
	 case("E3")
            if (parity_f == jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,6)) add = .false.
	 case("M3")
            if (parity_f /= jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,6)) add = .false.
	 case default
            stop "multipole_msymmetry_select(): program stop C."
	 end select
	 !
	 if (add) then
	    njnu             = njnu + 1
	    symmetry_j(njnu) = jnu
	    symmetry_p(njnu) = "+"
	 end if
	 !
	 !
	 add = .true.
	 jnu = j;   jparity = "-"
	 !
	 select case(mult1)
	 case("E1")
            if (parity_i == jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,2)) add = .false.
	 case("M1")
            if (parity_i /= jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,2)) add = .false.
	 case("E2")
            if (parity_i /= jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,4)) add = .false.
	 case("M2")
            if (parity_i == jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,4)) add = .false.
	 case("E3")
            if (parity_i == jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,6)) add = .false.
	 case("M3")
            if (parity_i /= jparity   .or.  &
                .not.is_triangle(totalJ_i,jnu,6)) add = .false.
	 case default
            stop "multipole_msymmetry_select(): program stop D."
	 end select
	 !
	 select case(mult2)
	 case("E1")
            if (parity_f == jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,2)) add = .false.
	 case("M1")
            if (parity_f /= jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,2)) add = .false.
	 case("E2")
            if (parity_f /= jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,4)) add = .false.
	 case("M2")
            if (parity_f == jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,4)) add = .false.
	 case("E3")
            if (parity_f == jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,6)) add = .false.
	 case("M3")
            if (parity_f /= jparity   .or.  &
                .not.is_triangle(totalJ_f,jnu,6)) add = .false.
	 case default
            stop "multipole_msymmetry_select(): program stop E."
	 end select
	 !
	 if (add) then
	    njnu             = njnu + 1
	    symmetry_j(njnu) = jnu
	    symmetry_p(njnu) = "-"
	 end if
	 !
      end do
      !
   end subroutine multipole_msymmetry_select
   !
   !
   subroutine multipole_select(totalJ_i,parity_i,totalJ_f,parity_f,    &
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
         select case( multipoles(i) )
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
            stop "multipole_select(): program stop A."
         end select
      end do
      !
      if (nmult > 20) then
         stop "multipole_select(): program stop B."
      end if
      !
   end subroutine multipole_select
   !
end module rabs_multipole
