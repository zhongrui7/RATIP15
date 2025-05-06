module rabs_toolbox_tpdi
!
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains all procedures which help carry out a number of
! 'utility' tasks within the RATIP environment. A few of these procedures 
! occur in a similar form also in other modules of the RATIP package 
! but have been adapted for the present purpose so that the utility 
! program does by itself not depend on too many modules.
!-----------------------------------------------------------------------
   !
   use rabs_constant
   use rabs_file_handling
   use rabs_toolbox_aux
   use rabs_xl
   implicit none
   !
   type(csf_basis), public :: csf_set_a, csf_set_b, csf_set_c
   !
   !
   ! Define an internal structure which stores all information about the
   ! photoionization lines of a given step
   type :: toolbox_p_channel
      integer          :: kappa, totalJ   
      character(len=1) :: parity 
      character(len=2) :: multipole  
      character(len=9) :: gauge 
      real(kind=dp)    :: phase, amplitude_re  
      complex(kind=dp) :: amplitude   
   end type toolbox_p_channel
   !
   type :: toolbox_p_line
      integer          :: asfi, asff, level_i, level_f, totalJ_i, totalJ_f
      integer          :: No_channels
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: p_energy, e_energy, cs_coulomb, cs_babushkin 
      complex(kind=dp) :: beta_b, beta_c, xi_b, xi_c, eta_b, eta_c, &
                          zeta_b, zeta_c, alignment_b, alignment_c
      type(toolbox_p_channel), dimension(:), pointer :: channel
   end type toolbox_p_line
   !
   type(toolbox_p_line), dimension(:), allocatable :: line_step1, line_step2
   !
   type :: toolbox_p_single_step
      integer          :: asfi, asff, totalJ_i, totalJ_f, tr
      character(len=1) :: parity_i, parity_f
      real(kind=dp)    :: p_energy, e_energy
      real(kind=dp)    :: cs_b, cs_c, beta_b, beta_c, alignment_b, alignment_c 
      ! First-step properties
      real(kind=dp)    :: b2_b, b2_c, b4_b, b4_c
      ! Second-step properties
      real(kind=dp)    :: a0_b, a0_c, a2_b, a2_c, a4_b, a4_c
      real(kind=dp)    :: beta_2_TPDI_b,beta_2_TPDI_c,beta_4_TPDI_b,beta_4_TPDI_c
      real(kind=dp)    :: cs_0if_b, cs_0if_c
   end type toolbox_p_single_step
   !
   type :: toolbox_p_cascade
      ! Cascade properties
      real(kind=dp)    :: beta2_0i_f_b, beta2_0i_f_c, &
                          beta4_0i_f_b, beta4_0i_f_c 
   end type toolbox_p_cascade
   !
   !
   type :: toolbox_p_2step
      integer          :: ni, nn, nf
      real(kind=dp)    :: energy
      type(toolbox_p_single_step), dimension(5,5)   :: step1, step2
      type(toolbox_p_cascade),     dimension(5,5,5) :: cascade12
   end type toolbox_p_2step
   !
   type(toolbox_p_2step) :: photon
   !
   !
   ! Define an internal structure to deal with the coherence transfer and
   ! the calculation of angular correlation functions
   type :: toolbox_coherence_B
       real, dimension(0:6,0:6,0:6) :: value
   end type toolbox_coherence_B
   !
   type :: toolbox_coherence
      integer               :: nb_max = 6, nc_max = 10
      integer               :: ni, nn, nf, asf_i, asf_f
      integer               :: totalJ_i, totalJ_f
      integer, dimension(5) :: asf_n, totalJ_n, tr_in, tr_nf
      character(len=1)      :: parity_i, parity_f 
      !
      character(len=1), dimension(5)                :: parity_n
      real(kind=dp), dimension(5,5)                 :: h
      type(toolbox_coherence_B), dimension(0:5,0:5) :: B, Bbar
      logical, dimension(0:10,0:10,0:10)            :: need_C
      real(kind=dp), dimension(0:10,0:10,0:10)      :: C
      !
      integer                         :: no_angles
      real(kind=dp)                   :: theta_a, phi_a, theta_b, phi_b
      real(kind=dp), dimension(50)    :: theta, phi
      real(kind=dp), dimension(50,50) :: W
   end type toolbox_coherence
   !
   type(toolbox_coherence), save :: coherence
   !
contains
   !
end module rabs_toolbox_tpdi
