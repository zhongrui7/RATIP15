!
!***** July 2010 *****
!
subroutine relci_op(n,m,current,new)
!--------------------------------------------------------------------
! Up-dates the Hamiltonian matrix in course of the Davidson diagonalization.
!--------------------------------------------------------------------
   ! 
   use rabs_constant 
   use rabs_hamiltonian
   ! 
   implicit none
   integer, intent(in)                        :: n, m
   real(kind=dp), dimension(n,m), intent(in)  :: current
   real(kind=dp), dimension(n,m), intent(out) :: new
   !
   integer, save               :: count = 0
   !
   integer                     :: i, j, k, p, non_zero
   integer, dimension(n)       :: ndx
   real(kind=dp), dimension(n) :: me
   !
   count = count + 1
   !
   if (mod(count,10) == 0) then
      print *, count, "-th call of relci_op"
   end if
   !
   new(:,:) = zero
   !
   do  j = 1,n
      call hamiltonian_half_column(j,n,ndx,me,non_zero)
      do  k = 1,m
         do  p = 2,non_zero
            i        = ndx(p)
            new(i,k) = new(i,k) + me(p) * current(j,k)
            new(j,k) = new(j,k) + me(p) * current(i,k)
         end do
         if (rabs_use_stop   .and.   ndx(1) /= j) then
            stop "relci_op() - program stop A"
         end if
         new(j,k) = new(j,k) + me(1) * current(j,k)
      end do
   end do
   !
end subroutine relci_op
