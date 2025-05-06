module rabs_toolbox_mess
!
!***** July 2010 *****
!-----------------------------------------------------------------------
! This module contains all the messages which explain the individual tasks
! and which are printed to screen during the executation.
!-----------------------------------------------------------------------
   !
   use rabs_constant
   implicit none
   !
   !
   !
   integer, public :: xmessage_na01 = 2
   character(len=100), dimension(2) :: xmessage_a01 = (/                                                  &
      "Returns the energy levels and level splitting for all given ASF.                             ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_na02 = 4
   character(len=100), dimension(4) :: xmessage_a02 = (/                                                  &
      "Returns the energy levels and excitation energies of some ASF from one or several .mix files.",    &
      "For the printout of the excitation energies, the energy units, an ascending/decending order, ",    &
      "and the level, relative to which the excitation energies are taken, can be selected below.   ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_na03 = 3
   character(len=100), dimension(3) :: xmessage_a03 = (/                                                  &
      "Returns the energy levels and level splitting for all ASF from a single .inp and mix file    ",    &
      "together with the two leading LSJ configurations.                                            ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_na04 = 4
   character(len=100), dimension(4) :: xmessage_a04 = (/                                                  &
      "Displays the main CSF and their weights to atomic levels. In the printout,                   ",    &
      "the level numbers and the number of leading CSF can be chosen explicitly. As default,        ",    &
      "the five leading CSF are displayed for all levels of the given .mix mixing coefficient file. ",    &     
      "                                                                                             " /)
   !
   integer, public :: xmessage_na05 = 5
   character(len=100), dimension(5) :: xmessage_a05 = (/                                                  &
      "Transform one or several ASF in a jj-coupled CSF basis from a GRASP2K calculation into a     ",    &
      "LS-coupled CSF basis in order to determine the (main) LS symmetry terms. The transformation  ",    &
      "starts from the given .csl and .mix files and calculates the mixing coefficients and weigths ",    &
      "in LS-coupling, similarly as in standard GRASP2K computations.                               ",    &
      "                                                                                             " /)
   !
   !
   !
   integer, public :: xmessage_nc01 = 15
   character(len=100), dimension(15) :: xmessage_c01 = (/                                                 &
      "Reduce a given .csl list from GRASP2K due to a set of 'restriction rules' which can be       ",    &
      "supplied interatively. These 'rules' can be defined in several ways; their general form      ",    &
      "is:                                                                                          ",    &   
      "          subshell_1  +  subshell_2  +  ...  +  subshell_n     rel_op    N                   ",    &  
      "where  subshell_i  is 1s, 2p, ...,  and  rel_op  denotes a relational operator,              ",    &
      "i.e. one of {<, <=, <$, =, >, >=} while N (> 0) is the number of allowed particles in        ",    &
      "these subshells. A  *  denotes as usual all subshells of a given layer; each restriction     ",    &
      "rule must be given on a separate line and may contain only peel shell orbitals.              ",    &
      "The relational operator  <$  means that all shells of this restriction must contribute       ",    &
      "to it with an occupation > 0.                                                                ",    & 
      "                                                                                             ",    & 
      "A few examples are:                                                                          ",    &
      "            4s + 4p- + 4d <= 3      or      5* < 2   or   3p + 3d > 12                       ",    &
      "            4s + 4p <$ 2   ... in order to omit 4s(1) 4p(1) but not 4s(2) or 4p(2).          ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nc02 = 4
   character(len=100), dimension(4) :: xmessage_c02 = (/                                                  &
      "'Split' a given  .csl  files from GRASP2K into separate lists owing to their J^P symmetry.   ",    &
      "For each of these symmetries, the program prompts for a file name. It is assumed, moreover,  ",    &
      "that the peel orbitals are 'occupied' in all the symmetry blocks.                            ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nc03 = 8
   character(len=100), dimension(8) :: xmessage_c03 = (/                                                  &
      "'Merge' two independent .csl  files from GRASP2K together. From the 'combined' list,         ",    &
      "all doubly defined CSF are removed so that no linear dependence may arise.                   ",    &
      "When compared with the GRASP2K component 'mrgcsl', this subtask proceeds much faster and     ",    &
      "supports a larger variety of different  .csl  files which need, for instance, not to have    ",    &
      "the same set of subshells in both lists. However, there is still a restrictions at present   ",    &
      "in that both .csl files must have the same 'core' and the same order of those orbitals       ",    &
      "which occur in both lists.                                                                   ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nc04 = 7
   character(len=100), dimension(7) :: xmessage_c04 = (/                                                  &
      "'Merge' two independent .csl  files from GRASP2K with different cores together. From the     ",    &
      "'combined' list, all doubly defined CSF are removed so that no linear dependence may arise.  ",    &
      "When compared with the GRASP2K component 'mrgcsl', this subtask proceeds much faster and     ",    &
      "supports a larger variety of different  .csl  files which need, for instance, not to have    ",    &
      "the same set of subshells in both lists. Althouhg the two .csl files need not to have the    ",    &
      "same core, the same order of orbitals is still required.                                     ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nc05 = 4
   character(len=100), dimension(4) :: xmessage_c05 = (/                                                  &
      "Condense a GRASP2K  .csl  list due to a given cutoff 0 < w <= 1. All CSF with a              ",    &
      "weight = [c_r (alpha)]^2 < w (for any of the atomic levels alpha in the .mix file) will      ",    &
      "be removed from the configuration list.                                                      ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nc06 = 7
   character(len=100), dimension(7) :: xmessage_c06 = (/                                                  &
      "Condense a GRASP2K  .csl  list due to a given cutoff 0 < w <= 1. All CSF with a              ",    &
      "weight = [c_r (alpha)]^2 < w for any of the specified atomic levels in the  .mix  file will  ",    &
      "be removed from the configuration list. The user can also create a 'new'  .mix  mixing       ",    &
      "coefficient file which contains the re-normalized eigenvectors of the specified levels.      ",    &
      "Note that this is of course only an approximation to the exact eigenvectors in the given     ",    &
      "basis.                                                                                       ",    & 
      "                                                                                             " /)
   !
   integer, public :: xmessage_nc07 = 3
   character(len=100), dimension(3) :: xmessage_c07 = (/                                                  &
      "Reduces a .mix file based on ASF serial numbers. The program reads in one .mix file and      ",    &	
      "prints those ASF as specified by the user; a formatted .mix file is generated.               ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nc08 = 4
   character(len=100), dimension(4) :: xmessage_c08 = (/                                                  &
      "Reduces a .mix file based on symmetry (J and P). The program reads in one .mix file and      ",    &	
      "prints those ASF for one or several J^P pairs as specified by the user;                      ",    &
      "a formatted .mix file is generated.                                                          ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nc09 = 6
   character(len=100), dimension(6) :: xmessage_c09 = (/                                                  &
      "Generate a pair-correlation list in order to exclude all CSF in a .csl list which are        ",    &
      "not coupled directly by the Hamiltonian to a given set of reference configurations.          ",    &
      "Such a  .csl  list can be generated either on the basis of the Dirac-Coulomb or              ",    &
      "Dirac-Coulomb-Breit Hamiltonian. The subshells in the reference list must have the           ",    &
      "same order as the list which is to be reduced.                                               ",    &
      "                                                                                             " /)
   !
   !
   !
   integer, public :: xmessage_nd01 = 2
   character(len=100), dimension(2) :: xmessage_d01 = (/                                                  &
      "Format an 'unformatted' GRASP2K  .mix  Mixing Coefficient File.                              ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nd02 = 2
   character(len=100), dimension(2) :: xmessage_d02 = (/                                                  &
      "Format an 'unformatted' GRASP2K  .out Radial Wavefunction File.                              ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nd03 = 2
   character(len=100), dimension(2) :: xmessage_d03 = (/                                                  &
      "Unformat a 'formatted' GRASP2K  .mix  Mixing Coefficient File.                               ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nd04 = 2
   character(len=100), dimension(2) :: xmessage_d04 = (/                                                  &
      "Unformat a 'formatted' GRASP2K  .out Radial Wavefunction File.                               ",    &
      "                                                                                             " /)
   !
   !
   !
   integer, public :: xmessage_ne01 = 5
   character(len=100), dimension(5) :: xmessage_e01 = (/                                                  &
      "Calculates and writes out the one-electron energies and radial expectation values of a set   ",    &
      "of orbitals as given by a  .mix  GRASP2K mixing coefficient file. The expectation values     ",    &
      "are printed for the operators  r^k  with powers k = -3, -1,  1, and 2 which are related to   ",    &
      "various physical operators in atomic structure and collision theory.                         ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_ne02 = 3
   character(len=100), dimension(3) :: xmessage_e02 = (/                                                  &
      "Calculates the 'overlaps' of orbitals with the same symmetry for two not quite orthogonal    ",    &
      "sets of orbitals. However, the orbitals within each set are supposed to be orthogonal.       ",    &
      "                                                                                             " /)
   !
   !
   !
   integer, public :: xmessage_nm01 = 4
   character(len=100), dimension(4) :: xmessage_m01 = (/                                                  &
      "Determines the effective radial charge or the charge density of a given radial orbital.      ",    &
      "The effective charge is given by the integral over the charge density up to a given          ",    &
      "radius R. The charge density of a particular level can also be printed to some text file.    ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nm02 = 5
   character(len=100), dimension(5) :: xmessage_m02 = (/                                                  &
      "Determines the effective radial charge or the charge density of a selected level.            ",    &
      "The effective charge is given by the integral over the charge density up to a given          ",    &
      "radius R. The charge density of a particular level can also be printed to some text file.    ",    &
      "To calculate the charge density, a full representation of the wave functions is required.    ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nm04 = 5
   character(len=100), dimension(5) :: xmessage_m04 = (/                                                  &
      "Supports the determination and interpretation of proper grid parameters for various          ",    &
      "information about the grid, the free-electron energy, or others.                             ",    &
      "If the kinetic energy, E_kin, of a free-electron is given, hp is determined so that one      ",    &
      "obtains a grid with 100 point per wavelength, and vice versa for determining E_kin.          ",    &
      "                                                                                             " /)
   !
   integer, public :: xmessage_nm05 = 5
   character(len=100), dimension(5) :: xmessage_m05 = (/                                                  &
      "Estimate the self-energy shift for a given set of ASF; for that, the program requires the    ",    &
      ".csl and .mix files as well as the .out radial orbital file. The self-energy estimate is     ",    &
      "based on a procedure as suggested by Yong-Ki Kim and described in the long write-up for      ",    &
      "RELCI.                                                                                       ",    &
      "                                                                                             " /)
   !
   !
   !
   integer, public :: xmessage_nn01 = 2
   character(len=100), dimension(2) :: xmessage_n01 = (/                                                  &
      "Calculate the nuclear radius from the mass number by using the simple formula:               ",    &
      "     R = sqrt( <r^2> ) = 0.836 * A^1/3  +  0.570.                                            " /)
   !
   integer, public :: xmessage_nn02 = 7
   character(len=100), dimension(7) :: xmessage_n02 = (/                                                  &
      "Calculate the Fermi distribution parameters c and N as well as the Fermi distribution and the",    &
      "associated potential (if appropriate). These computations follow the formulas by             ",    &
      "Tupitsyn et al., PRA 68 (2003) 022511, Eqs. (8-11). For a given Fermi parameter a,           ",    &
      "in particular, the parameter c and normalization N are calculated from analytical formulas.  ",    &
      "Using these parameters, the distribution:                                                    ",    &
      "     rho(r; R) = N / ( 1 + exp[(r-c)/a] )                                                    ",    &
      "as well as the associated nuclear potential V_n (r; R) can be printed out on request.        " /)
   !  
   !
   !
contains
   !
   subroutine toolbox_mess()
   !--------------------------------------------------------------------
   !
   ! Calls: 
   !--------------------------------------------------------------------
      !
      print *, "**************************"
      print *, "*** Not yet implmented ***"
      print *, "**************************"
      !
   end subroutine toolbox_mess
   !
end module rabs_toolbox_mess
