*=======================================================================
	subroutine addabs(op,n,lim,hiend,kpass,ncnv,basis,ab,s)
*=======================================================================
*	Called by: DVDRVR, SETUP
*
*	Calculates the new column in the D matrix and the new column
*	in the S matrix. The new D column is D(new)=AB(new). S has a 
*	new row and column, but being symmetric only the new column is 
*	stored. S(i,kpass+1)=B(i)^T D(kpass+1) for all i.
*
*	subroutines called:
*	OP, DDOT, DSCAL
*	
*-----------------------------------------------------------------------
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION BASIS(N*LIM),AB(N*LIM)
	DIMENSION S(LIM*(LIM+1)/2)
	LOGICAL HIEND
	EXTERNAL OP
*-----------------------------------------------------------------------
*   on entry
*   -------
*   N		The order of the matrix A
*   kpass	The current dimension of the expanding sub-basis
*   NCNV	Number of new basis vectors.
*   Basis	the basis vectors, including the new NCNV ones.
*   on exit
*   -------
*   AB		The new matrix D=AB. (with new NCNV columns)
*   S           The small matrix with NCNV new columns at the last part
*-----------------------------------------------------------------------
*
* The user specified matrix-vector routine is called with the new 
* basis vector B(*,kpass+1) and the result is assigned to AB(idstart) 
*
        IDSTART=KPASS*N+1
*        print *, "IDSTART= ",IDSTART
*        print *, "BASIS = ",basis
*        print *, "ab = ",ab
	CALL OP(N,NCNV,BASIS(IDSTART),AB(IDSTART))
*
* If highest pairs are sought, use the negative of the matrix
*
	IF (HIEND) CALL DSCAL(N*NCNV,-1.D0,AB(IDSTART),1)
*
* The new S is calculated by adding the new last columns
* S(new)=B^T D(new).
*
	ISSTART=KPASS*(KPASS+1)/2
	DO 20 IV=1,NCNV
	   IBSTART=1
	   DO 10 IBV=1,KPASS+IV
	       SS=DDOT(N,BASIS(IBSTART),1,AB(IDSTART),1)
	       S(ISSTART + IBV)=SS
	       IBSTART=IBSTART+N
 10        CONTINUE
	   ISSTART=ISSTART+KPASS+IV
	   IDSTART=IDSTART+N
 20	CONTINUE

	RETURN
	END
