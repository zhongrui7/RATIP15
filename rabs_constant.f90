module rabs_constant
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module defines all processor-dependent constants as well as
! a large variety of physical constants and conversion factors.
! All processor dependent procedures should be placed in this module.
! It is used by almost all other modules. 
!
! This module further includes a variety of useful type definitions for
! relativistic atomic structure calculations.
!-----------------------------------------------------------------------
   !
   private :: is_equal_nkappa
                 ! "Overloads" the equivalence (==) of two variables of 
                 ! type(nkappa).
   private :: is_equal_nkappam
                 ! "Overloads" the equivalence (==) of two variables of 
                 ! type(nkappam).
   private :: is_equal_nl
                 ! "Overloads" the equivalence (==) of two variables of 
                 ! type(nl).
   private :: is_equal_nlmms
                 ! "Overloads" the equivalence (==) of two variables of 
                 ! type(nlmms).
   public  :: print_physical_constants             
                 ! Prints out a short summary about all physical constants.
   public  :: print_conversion_factors
                 ! Prints out a short summary about all physical conversion 
                 ! factors.                  
   public  :: get_yes_stream        
                 ! Supports input.
   public  :: save_input       
                 ! Helps collect and print all the input data of some
		 ! computation.
   !
   integer, parameter, public :: i6b = selected_int_kind(12)
   integer, parameter, public :: i4b = selected_int_kind(9)
   integer, parameter, public :: i2b = selected_int_kind(4)
   integer, parameter, public :: i1b = selected_int_kind(2)
   !
   integer, parameter, public :: sp  = kind(1.0)
   integer, parameter, public :: dp  = selected_real_kind(2*precision(1.0_sp))
   integer, parameter, public :: qp  = selected_real_kind(2*precision(1.0_sp))
   !! integer, parameter, public :: qp  = selected_real_kind(2*precision(1.0_dp))
   !
   type, public :: matrix_dp
      real(kind=dp), dimension(:,:), pointer :: matrix
   end type matrix_dp
   type, public :: vector_dp
      real(kind=dp), dimension(:),   pointer :: vector
   end type vector_dp
   !
   ! Define derived structures for (relativistic) orbital indices
   type, public :: nkappa
      integer :: n, kappa
   end type nkappa
   !
   type, public :: nkappam
      integer :: n, kappa, mm
   end type nkappam
   !
   type, public :: nl
      integer :: n, l
   end type nl
   !
   type, public :: nlmms
      integer :: n, l, m, mms
   end type nlmms
   !
   interface operator(==)
      module procedure is_equal_nkappa
      module procedure is_equal_nkappam
      module procedure is_equal_nl
      module procedure is_equal_nlmms
   end interface
   !
   ! Define simple constants
   real(kind=dp), parameter, public :: one = 1.0_dp, two = 2.0_dp,  &
       three = 3.0_dp, four = 4.0_dp, five = 5.0_dp, ten = 10.0_dp, &
       zero  = 0.0_dp, half = one/two, third = one/three,           &
       eps10 = 1.0e-10_dp, eps20 = 1.0e-20_dp
   !
   ! Define a logical array to control the calculation and printout
   logical, dimension(1:100), public :: options
   !
   ! Define a logical variable and a corresponding array to control 
   ! a possible debugging of the program
   integer, public                   :: rabs_dbg = 99
   logical, public                   :: rabs_debug = .false.
   logical, dimension(1:100), public :: debugging
   !
   ! Define a logical flag for including a (large) number of additional
   ! tests and 'stops' if data seem not to be appropriate; this helps in
   ! the development and debugging of the code
   logical, parameter, public :: rabs_use_stop = .true.
   !
   ! Define a logical flag for using the NAG library
   logical, parameter, public :: rabs_use_naglib = .false.
   !
   ! Define a logical flag for using finite-difference features (GRASP92) 
   logical, parameter, public :: rabs_use_finite_difference = .true.
   !
   ! Define a logical flag for using basis sets 
   logical, parameter, public :: rabs_use_basis_set = .false.
   !
   ! Define some global variables to deal with (energy) units
   logical, public          :: energy_inverse = .false.
   real(kind=dp), public    :: energy_factor  = zero
   character(len=7), public :: energy_unit
   !
   ! Define some global variables as used in GRASP92
   integer, public       :: n_grasp2k    = 390, kmp1_grasp2k = 60 
   real(kind=dp), public :: rnt_grasp2k  = 2.0e-6_dp, &
                            h_grasp2k    = 5.0e-2_dp, &
                            hp_grasp2k   = zero
   !
   ! Define an input stream to convert some data into the string format
   integer, parameter, public :: stream_input = 57
   !
   ! Define some useful physical constants
   !
   real(kind=dp), parameter, public ::                 &
      bohr_radius_in_cm      =    0.529177249e-8_dp,   &
      ! one_over_alpha       =  137.0373_dp,           & !Old value from Grasp-2
      one_over_alpha         =  137.0359895_dp,        &
      c_vacuum_in_cm_per_s   =    2.99792458e10_dp,    &
      electron_charge_in_esu =    4.80320680e-10_dp,   &
      electron_mass_in_g     =    9.1093897e-28_dp,    &
      electron_mass_in_amu   =    5.48579903e-4_dp,    &
      proton_mass_in_amu     =    1.007276470_dp,      &
      hbar_in_ergs           =    1.05457266e-27_dp,   &
      rydberg_in_ev          =   13.6056981_dp,        &
      rydberg_in_kaysers     =  109.73731534e3_dp
   !   
   ! The values of the physical constants used here are taken from
   ! E R Cohen and B N Taylor, The 1986 Adjstment of the Fundamental
   ! Physical Constants, Report of the CODATA Task Group on Funda-
   ! mental COnstants, CODATA Bulletin 63, Pergamon, Elmsford, NY (1986)
   !
   ! In the context of new results, Taylor warns that ... since the
   ! output values of a least-squares adjustment are related in a
   ! complex way and a change in the measured value of one constant
   ! usually leads to corresponding changes in the adjusted values
   ! of others, one must be cautious in carrying out calculations using both 
   ! the [above] values and the results of more recent measurements.
   !
   ! convert_au_to_kaysers and convert_au_to_ev assume an infinitely-heavy
   ! nucleus; convert_einstein_a_to_si and convert_einstein_b_to_si convert
   ! to SI units
   !
   real(kind=dp), parameter, public ::                                      &
      convert_au_to_kaysers    = two * rydberg_in_kaysers,                  &
      convert_au_to_ev         = two * rydberg_in_ev,                       &
      convert_einstein_a_to_si = (electron_mass_in_g/hbar_in_ergs) *        &
         ((electron_charge_in_esu*electron_charge_in_esu/hbar_in_ergs)**2), &
      convert_einstein_b_to_si = (ten * bohr_radius_in_cm**3/hbar_in_ergs)  &
         * convert_einstein_a_to_si,                                        &
      convert_fermi_to_bohr    = 1.0e-13_dp / bohr_radius_in_cm,            &
      convert_au_to_per_sec    = 4.134138e16
   !
   real(kind=dp), parameter, public :: pi = two * 1.570796326794896619231_dp
   real(kind=dp), parameter, public :: c_vacuum = one_over_alpha
   real(kind=dp), save, public      :: c
   !
   logical, public :: debug_physical_constants = .false.,  &
                      debug_conversion_factors = .false.
   ! 
contains
   !
   !
   function is_equal_nkappa(a,b)                             result(yes)
   !--------------------------------------------------------------------
   ! This function "defines" the equivalence of two variables of 
   ! type(nkappa).
   !--------------------------------------------------------------------
      !
      implicit none
      type(nkappa), intent(in) :: a, b
      logical                  :: yes
      !
      if (a%n == b%n   .and.   a%kappa == b%kappa) then
         yes = .true.
      else
         yes = .false.
      end if
      !
   end function is_equal_nkappa
   !
   !
   function is_equal_nkappam(a,b)                            result(yes)
   !--------------------------------------------------------------------
   ! This function "defines" the equivalence of two variables of 
   ! type(nkappam).
   !--------------------------------------------------------------------
      !
      implicit none
      type(nkappam), intent(in) :: a, b
      logical                   :: yes
      !
      if (a%n == b%n   .and.   a%kappa == b%kappa  .and.  a%mm == b%mm) then
         yes = .true.
      else
         yes = .false.
      end if
      !
   end function is_equal_nkappam
   !
   !
   function is_equal_nl(a,b)                                 result(yes)
   !--------------------------------------------------------------------
   ! This function "defines" the equivalence of two variables of 
   ! type(nl).
   !--------------------------------------------------------------------
      !
      implicit none
      type(nl), intent(in) :: a, b
      logical              :: yes
      !
      if (a%n == b%n   .and.   a%l == b%l) then
         yes = .true.
      else
         yes = .false.
      end if
      !
   end function is_equal_nl
   !
   !
   function is_equal_nlmms(a,b)                              result(yes)
   !--------------------------------------------------------------------
   ! This function "defines" the equivalence of two variables of 
   ! type(nlmms).
   !--------------------------------------------------------------------
      !
      implicit none
      type(nlmms), intent(in) :: a, b
      logical                 :: yes
      !
      if (a%n == b%n   .and.   a%l == b%l  .and.  a%m == b%m  .and.  & 
          a%mms == b%mms) then
         yes = .true.
      else
         yes = .false.
      end if
      !
   end function is_equal_nlmms
   !
   !
   subroutine print_physical_constants(unit)
   !--------------------------------------------------------------------
   ! Print a short summary about all physical constants as currently
   ! defined in the rabs package.
   !--------------------------------------------------------------------
      !
      integer :: unit
      !
      write(unit,*) "Physical constants are currently defined as: "
      write(unit,*) "-------------------------------------------- "
      write(unit,*) " "
      write(unit,"(a,es15.9)")  "   Bohr radius (cm) :                      ",&
                                    bohr_radius_in_cm
      write(unit,"(a,es15.9)")  "   Inverse of the fine-structure constant: ",&
                                    one_over_alpha
      write(unit,"(a,es14.8)")  "   Speed of light (cm/s) :                 ",&
                                    c_vacuum_in_cm_per_s                           
      write(unit,"(a,es14.8)")  "   Electron charge (esu) :                 ",&
                                    electron_charge_in_esu                           
      write(unit,"(a,es13.7)")  "   Electron mass (g) :                     ",&
                                    electron_mass_in_g                           
      write(unit,"(a,es15.9)")  "   Electron mass (u) :                     ",&
                                    electron_mass_in_amu                           
      write(unit,"(a,es15.9)")  "   Proton mass (u) :                       ",&
                                    proton_mass_in_amu                           
      write(unit,"(a,es14.8)")  "   Rationalized Planck constant (erg s) :  ",&
                                    hbar_in_ergs                           
      write(unit,"(a,es14.8)")  "   Rydberg (eV) :                          ",&
                                    rydberg_in_ev                           
      write(unit,"(a,es16.10)") "   Rydberg (Kaysers) :                     ",&
                                    rydberg_in_kaysers 
      write(unit,*) " "
      !
   end subroutine print_physical_constants 
   !                        
   !
   subroutine print_conversion_factors(unit)
   !--------------------------------------------------------------------
   ! Print a short summary about all physical conversion factors as 
   ! currently defined in the rabs package.
   !--------------------------------------------------------------------
      !
      integer :: unit
      !
      write(unit,*) "Conversion factors are currently defined as: "
      write(unit,*) "-------------------------------------------- "
      write(unit,*) " "
      write(unit,"(a,es15.9)")  "   a.u. to Kaysers :      ", &
                                    convert_au_to_kaysers
      write(unit,"(a,es15.9)")  "   a.u. to eV :           ",convert_au_to_ev
      write(unit,"(a,es15.9)")  "   Einstein-A to SI :     ", &
                                    convert_einstein_a_to_si
      write(unit,"(a,es15.9)")  "   Einstein-B to SI :     ", &
                                    convert_einstein_b_to_si
      write(unit,"(a,es15.9)")  "   Fermi to a.u. (Bohr) : ", &
                                    convert_fermi_to_bohr
      write(unit,*) " "
      !
   end subroutine print_conversion_factors
   !                        
   !
   function get_yes_stream()   result(yes)
   !--------------------------------------------------------------------
   ! This function reads a response on the default input unit; this 
   ! response must be either 'y' or 'n'. It returns  .true.  if 'y' is   
   ! entered and  .false.  if 'n' is entered.  
   !--------------------------------------------------------------------
      !
      logical          :: yes
      integer          :: ios
      character(len=1) :: response
      !
    2 read (unit=*,fmt="(a)",iostat=ios) response
      if (ios /= 0) then
         print *, "Expecting <y><cr> or <n><cr> ..."
         read (unit= *,fmt="(a)",iostat=ios) response
         if (ios /= 0) then
            print *, "Expecting <y><cr> or <n><cr> ..."
            read (unit=*,fmt="(a)",iostat=ios) response
         endif
      end if
      !
      if (response == "y") then
         yes = .true.
      else if (response == "n") then
         yes = .false.
      else
         print *, "Expecting <y><cr> or <n><cr> ..."
         goto 2
      end if
      !
      call save_input(response,.false.)
      !
   end function get_yes_stream
   !
   !
   subroutine save_input(string,prnt,stream)
   !--------------------------------------------------------------------
   ! Helps collect all the input data of some computation in order to
   ! print it later to .sum files.
   !--------------------------------------------------------------------
      !
      character(len=*), intent(in)  :: string
      logical, intent(in)           :: prnt
      integer, optional, intent(in) :: stream
      !    
      integer, save                            :: No_lines = 0
      character(len=100), dimension(200), save :: sinput = " "
      !
      if (prnt) then
         write(stream,*) " "
         write(stream,*) "The following input data were used for the "//&
	                 "computations:"
         write(stream,*) "-------------------------------------------"//&
	                 "-------------"
	 do  i = 1,No_lines
	    write(stream,*) sinput(i)(1:100)
	 end do
         write(stream,*) "--- end of input ---"
      else
         No_lines = No_lines + 1
	 sinput(No_lines)(:) = trim(string)
      end if
      !
   end subroutine save_input
   !                        
end module rabs_constant

