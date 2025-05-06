!***********************************************************************

! Code converted using TO_F90 by Alan Miller
! Date: 2020-02-25  Time: 21:28:06

!                                                                      *

SUBROUTINE bubble (jpol,fail)
!                                                                      *
!   Reduces a circuit of order 2 , giving  delta  function and phase   *
!   factors.                                                           *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: jpol
LOGICAL, INTENT(IN OUT)                  :: fail
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME,namsub
INTEGER :: arr,tab1
LOGICAL :: sumvar

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/nam/namsub  &
    /graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc  &
    /argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp

DATA NAME/'BUBBLE'/

namsub = NAME
k2 = 2
k23 = 3
i1 = 1
i2 = 1
it1 = npoint(1)
it2 = npoint(2)

IF (it2 == ilast) THEN
  IF (it1 /= ifirst) THEN
    it2 = it1
    it1 = ilast
  END IF
  i1 = -1
  k23 = 2
  i2 = 2
END IF

CALL phase (it1,jdiag,m4trd)
k = ABS ((3*arr(it2,1)+2*arr(it2,2)+arr(it2,3))/2)+1
IF (k /= 4) CALL phase2 (jdiag(it2,k))
IF (nbnode == 2) RETURN
il1 = il(it2)+i1
it = ih(il1)
arr(it,k23) = arr(it1,k23)
l = jdiag(it1,k23)
l1 = jdiag(it,k23)
jdiag(it,k23) = l

IF (jpol /= 1) THEN
  CALL delta (l,l1,fail)
  IF (fail) RETURN
ELSE
  mp = mp-1
  kw(2,jwc) = l
  j6(j6c-1) = l
  j6(j6c) = l
  IF (k == 2) j8(j8c) = l
END IF

tab1(l,i2) = it

IF (it1 /= ilast) THEN
  IF (it2 == ilast) THEN
    tab1(l,1) = ih(2)
    il1 = 2
    k2 = 1
  END IF
  
  DO  i = il1,nbnode
    it = ih(i)
    il(it) = i-k2
    ih(i-k2) = it
  END DO
  
END IF

j9(j9c+1) = l
j9c = j9c+2
j9(j9c) = l

RETURN
END SUBROUTINE bubble
!***********************************************************************
!                                                                      *

SUBROUTINE change (l,k)
!                                                                      *
!   Exchanges the free ends in either first or last triad of JDIAG.    *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: l
INTEGER, INTENT(IN OUT)                  :: k
IMPLICIT doubleprecision (a-h, o-z)
INTEGER :: arr,tab1

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc

CALL phase (l,jdiag,m4trd)
jp = jdiag(l,k)
jdiag(l,k) = jdiag(l,1)
jdiag(l,1) = jp
jar = arr(l,k)
arr(l,k) = arr(l,1)
arr(l,1) = jar

RETURN
END SUBROUTINE change
!***********************************************************************
!                                                                      *

SUBROUTINE chklp1 (fail)
!                                                                      *
!   This routine checks if  there are active triads with two identi-   *
!   cal  arguments.  This is a loop of  order 1 (a "lollypop").  The   *
!   other argument  must then be zero;  i.e. j1(j) = 1 in 2j+1 nota-   *
!   tion. Suppression of the loop introduces factors and phases. Two   *
!   Two triads  become inactive.  All this is  performed by invoking   *
!   ZERO with first argument 1.                                        *
!                                                                      *
!   Written by  Marcel Klapisch for correcting  an error detected by   *
!   Charlotte F Fischer.    This version includes Marcel's fix of 12   *
!   June 1992.                                                         *
!                                         Last revision: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


LOGICAL, INTENT(OUT)                     :: fail
IMPLICIT doubleprecision (a-h, o-z)
INTEGER :: arrow
LOGICAL :: free,tabs
CHARACTER (LEN=6) :: NAME,namsub

INTEGER, PARAMETER :: mangm = 60
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad

COMMON/couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3), free(mangm)  &
    /debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6 /nam/namsub  &
    /tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),  &
    lcol(mangm,2),tabs(m2trd),nbtr

DATA NAME/'CHKLP1'/

namsub = NAME
nbtr1 = 2*(n-1)
DO  l = 1,nbtr1
  IF (.NOT. tabs(l)) THEN
    jdif = 0
    IF     (j23(l,1) == j23(l,2)) THEN
      jdif = j23(l,3)
    ELSE IF (j23(l,1) == j23(l,3)) THEN
      jdif = j23(l,2)
    ELSE IF (j23(l,2) == j23(l,3)) THEN
      jdif = j23(l,1)
    END IF
    IF (jdif /= 0) THEN
      
!   Putting the link to 0. ZERO changes NBTR
      
      fail = .false.
      IF (j1(jdif) /= 1.AND. .NOT. free(jdif)) THEN
        fail = .true.
        IF (ibug3 == 1) WRITE (99,300) jdif,j1(jdif)
        RETURN
      ELSE
        CALL zero (1,jdif,fail)
        IF (fail) RETURN
      END IF
    END IF
  END IF
END DO
IF (jdif /= 0) CALL printj (NAME,4)

RETURN

300 FORMAT (1X,'JDIF = ',1I2,'; should be 0; J1(JDIF) = ',1I2,  &
    '; RECUP -> 0.')

END SUBROUTINE chklp1
!***********************************************************************
!                                                                      *

SUBROUTINE chvar (jp,nbc,kbc,jt,jinv,nsum)
!                                                                      *
!   Change  the  order  of  summation variable to be able to perform   *
!   separately the summations in GENSUM.                               *
!                                                                      *
!                                           Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: jp(nbc)
INTEGER, INTENT(IN)                      :: nbc
INTEGER, INTENT(IN OUT)                  :: kbc
LOGICAL, INTENT(IN OUT)                  :: jt(nsum)
INTEGER, INTENT(IN)                      :: jinv(nsum)
INTEGER, INTENT(IN OUT)                  :: nsum
IMPLICIT doubleprecision (a-h, o-z)




kb = kbc+1

IF (kb <= nbc) THEN
  DO  i = kb,nbc
    jk = jp(i)
    IF (jt(jk)) THEN
      kbc = kbc+1
      jp(i) = jp(kbc)
      jp(kbc) = jinv(jk)
    END IF
  END DO
END IF

RETURN
END SUBROUTINE chvar
!***********************************************************************
!                                                                      *

SUBROUTINE cut1l (fail)
!                                                                      *
!   Cut  on  one  line, that  was  left as a free end in JDIAG. Puts   *
!   corresponding delta in J23.                                        *
!                                                                      *
!                                           Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************


LOGICAL, INTENT(IN OUT)                  :: fail
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME
INTEGER :: arr,tab1
LOGICAL :: sumvar,free

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc  &
    /argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)

DATA NAME/'CUT1L '/

it = itfree(1)
j0 = jdiag(it,1)
CALL delta (j0,m,fail)
IF (fail) GO TO 2
CALL delta (jdiag(it,3),jdiag(it,2),fail)
IF (fail) GO TO 2
jdiag(it+1,3) = jdiag(it,3)

IF (arr(it,2) == arr(it,3)) THEN
  arr(it+1,3) = 1
  arr(it-1,2) = -1
ELSE IF (arr(it,2) < arr(it,3)) THEN
  arr(it+1,3) = -1
  arr(it-1,2) = 1
END IF

j9c = j9c+1
j9(j9c) = jdiag(it,3)
j = 2
CALL zero (j,j0,fail)
IF (fail) GO TO 2
il1 = il(it+1)

DO  i = il1,nbnode
  it = ih(i)
  ilp = i-1
  il(it) = ilp
  ih(ilp) = it
END DO

nbnode = nbnode-1

2 CALL printj (NAME,12)
RETURN

END SUBROUTINE cut1l
!***********************************************************************
!                                                                      *

SUBROUTINE cut2l (fail)
!                                                                      *
!   Cut on two lines that were left as free ends in JDIAG. Puts cor-   *
!   responding delta in J23.                                           *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


LOGICAL, INTENT(IN OUT)                  :: fail
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME
INTEGER :: arr,tab1,arrow
LOGICAL :: tabs,sumvar

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc  &
    /tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),  &
    lcol(mangm,2),tabs(m2trd),nbtr  &
    /argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp

DATA NAME/'CUT2L '/


it1 = itfree(1)
it2 = itfree(2)
jt1 = jdiag(it1,1)
jt2 = jdiag(it2,1)
CALL delta (jt1,jt2,fail)
IF (fail) GO TO 1
IF (arr(it1,1) == arr(it2,1)) CALL phase2 (jt1)
arr(it2,1) = -arr(it1,1)
jdiag(it2,1) = jt1
tab1(jt1,2) = it2
j9(j9c+1) = jt1
j9c = j9c+2
j9(j9c) = jt1
CALL otherj (0,jt1,l1,lc1,k1)
CALL otherj (0,jt2,l2,lc2,k2)
j23(l2,lc2) = jt1
line(jt1,k1) = l2
lcol(jt1,k1) = lc2
arrow(l2,lc2) = -arrow(l1,lc1)

1 CALL printj (NAME,12)

RETURN
END SUBROUTINE cut2l
!***********************************************************************
!                                                                      *

SUBROUTINE cutnl (fail)
!                                                                      *
!   This subroutine  examines the case where there are more than two   *
!   free ends, but they are contiguous, so that the graph can be cut   *
!   without destroying the flat structure.                             *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


LOGICAL, INTENT(OUT)                     :: fail
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME
INTEGER :: arrow,arr,tab1
LOGICAL :: tabs,sumvar

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),  &
    lcol(mangm,2),tabs(m2trd),nbtr  &
    /graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc  &
    /KEEP/jkp(2,3),jarr(2,3),it2,it3,it5  &
    /argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp

DATA NAME/'CUTNL '/

ntf = itfree(nfree)-itfree(1)
IF (ntf > nfree) GO TO 8
it2 = itfree(1)
it3 = itfree(nfree)
it1 = it2-1
it4 = it3+1

IF (ntf /= nfree) THEN
  
  jt = jdiag(it2,3)
  CALL delta (jt,jdiag(it3,2),fail)
  
  IF (fail) GO TO 9
  
  IF (arr(it2,3) == arr(it3,2)) THEN
    CALL phase2 (jt)
    arr(it2,3) = -arr(it2,3)
    arr(it1,2) = -arr(it1,2)
  END IF
  
  jdiag(it3,2) = jt
  jdiag(it4,3) = jt
  j9(j9c+1) = jt
  j9c = j9c+2
  j9(j9c) = jt
  nbtr = nbtr+nfree
  it5 = 0
  
ELSE
  
  nfr = 0
  
  DO  it5 = it2,it3
    nfr = nfr+1
    IF (itfree(nfr) > it5) GO TO 4
  END DO
  
  4    jkp(1,1) = jdiag(it5,1)
  jarr(1,1) = -arr(it5,1)
  jkp(1,2) = jdiag(it2,3)
  jarr(1,2) = -arr(it2,3)
  jkp(1,3) = jdiag(it3,2)
  jarr(1,3) = -arr(it3,2)
  
  DO  j = 1,3
    jkp(2,j) = jdiag(it5,j)
    jarr(2,j) = arr(it5,j)
  END DO
  
  jdiag(it5,2) = jdiag(it3,2)
  arr(it5,2) = arr(it3,2)
  jdiag(it5,3) = jdiag(it2,3)
  arr(it5,3) = arr(it2,3)
  ilp = il(it2)
  il(it5) = ilp
  ih(ilp) = it5
  nbtr = nbtr+nfree+2
  CALL phase (it5,jdiag,m4trd)
  k = ABS ((3*arr(it5,1)+2*arr(it5,2)+arr(it5,3))/2+1)
  IF (k /= 4) CALL phase2 (jdiag(it5,k))
  
END IF

il1 = il(it4)

DO  i = il1,nbnode
  it = ih(i)
  ilp = i-nfree
  il(it) = ilp
  ih(ilp) = it
END DO

nbnode = nbnode-nfree
nfin = 0
GO TO 8

9 fail = .true.
8 CALL printj (NAME,8)

RETURN

END SUBROUTINE cutnl
!***********************************************************************
!                                                                      *

SUBROUTINE delta (ja,jb,fail)
!                                                                      *
!   Test for delta(JA,JB). If they are summation variables, the sec-   *
!   ond  is  changed  into  the first everywhere. if they are fixed,   *
!   their value is checked, and fail put to .TRUE. if they differ.     *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN)                      :: ja
INTEGER, INTENT(IN)                      :: jb
LOGICAL, INTENT(OUT)                     :: fail
IMPLICIT doubleprecision (a-h, o-z)
LOGICAL :: cut,sumvar,free

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /cutdig/cut /debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6  &
    /dim/j6cc,j7cc,j8cc,j9cc,jwcc,jdelc  &
    /couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)

IF (ibug3 == 1) WRITE (99,300) ja,sumvar(ja),jb,sumvar(jb)
IF (sumvar(ja) .AND. sumvar(jb)) GO TO 2
IF (free(ja) .OR. free(jb)) THEN
  jdel = jdel+1
  ldel(jdel,1) = ja
  ldel(jdel,2) = jb
  sumvar(ja) = .false.
  sumvar(jb) = .false.
  RETURN
END IF

IF (j1(ja) /= j1(jb)) fail = .true.
cut = .true.
RETURN

2 IF (j6c /= j6cc) THEN
  j61 = j6cc+1
  
  DO  i = j61,j6c
    IF (j6(i) == jb) j6(i) = ja
  END DO
  
END IF

IF (j7c /= j7cc) THEN
  j71 = j7cc+1
  
  DO  i = j71,j7c
    IF (j7(i) == jb) j7(i) = ja
  END DO
END IF

IF (j8c /= j8cc) THEN
  j81 = j8cc+1
  
  DO  i = j81,j8c
    IF (j8(i) == jb) j8(i) = ja
  END DO
END IF

IF (j9c /= j9cc) THEN
  j91 = j9cc+1
  
  DO  i = j91,j9c
    IF (j9(i) == jb) j9(i) = ja
  END DO
END IF

IF (jwc /= jwcc) THEN
  jw1 = jwcc+1
  
  DO  i = jw1,jwc
    DO  j = 1,6
      IF (kw(j,i) == jb) kw(j,i) = ja
    END DO
  END DO
END IF

IF (jdel /= jdelc) THEN
  jdel1 = jdelc+1
  
  DO  i = jdel1,jdel
    DO  j = 1,2
      IF (ldel(i,j) == jb) ldel(i,j) = ja
    END DO
  END DO
  
  sumvar(jb) = .false.
END IF

RETURN

300 FORMAT (/'From DELTA: JA = ',i2,l2,5X,'JB = ',i2,l2)

END SUBROUTINE delta
!***********************************************************************
!                                                                      *

SUBROUTINE diagrm (jump)
!                                                                      *
!   This subroutine builds up a flat diagram from the triads J23 and   *
!   places them in JDIAG . Arrows  are in ARR (INTEGER). The diagram   *
!   is built so as to maximize the number of triads involved, within   *
!   a one-step-forward-check process. If the diagram does not inclu-   *
!   de all the NBTR triads, it will have 'free ends'. JDIAG has dim-   *
!   ension double that of  J23 , because the path may proceed either   *
!   way.                                                               *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: jump
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME
INTEGER :: arr,tab1,arrow
LOGICAL :: tabs,free,sumvar

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),  &
    lcol(mangm,2),tabs(m2trd),nbtr  &
    /graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc  &
    /build/ial(m4trd),if1,if2,node  &
    /couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)

DATA NAME/'DIAGRM'/

DATA nb/0/

!   Initialization

IF (jump > 2) GO TO 17
IF (jump < 2) nb = 0
1 nb = nb+1
IF (tabs(nb)) GO TO 1
node = nbtr
ilast = nbtr

DO  j = 1,3
  jdiag(node,j) = j23(nb,j)
  arr(node,j) = arrow(nb,j)
END DO

tabs(nb) = .true.

DO  i = 1,mp
  ial(i) = 0
END DO

if1 = jdiag(node,1)
if2 = jdiag(node,3)
ial(if1) = 1
ial(if2) = 1
17 ntime = 0
i1 = 1
k1 = 1
k2 = 2
k3 = 3
3 jb = jdiag(node,k2)
CALL otherj (0,jb,l,lc,kp)
CALL neibor (lc,l1,l2)

!   Check consistency of triads

IF (tabs(l)) THEN
  WRITE (*,300)
  STOP
END IF

CALL way (l,l1,l2,ich,nd)
node = node+i1
tabs(l) = .true.
jdiag(node,k3) = j23(l,lc)
arr(node,k3) = arrow(l,lc)
ict = ich*i1

IF (ich <= 0) THEN
  lp = l1
  l1 = l2
  l2 = lp
END IF

IF (ict <= 0) CALL phase (l,j23,m2trd)
jdiag(node,k1) = j23(l,l1)
arr(node,k1) = arrow(l,l1)
jdiag(node,k2) = j23(l,l2)
arr(node,k2) = arrow(l,l2)
j = j23(l,l1)
ial(j) = ial(j)+1
j = j23(l,l2)
ial(j) = ial(j)+1
IF (nd < 1) GO TO 3
ntime = ntime+1
ilast = MAX (node,ilast)
ifirst = MIN (node,nbtr)
nbp = ial(if1)+ial(if2)
IF ((nbp > 3) .OR. (ntime > 1)) THEN
  nbnode = ilast-ifirst+1
  nbtr = nbtr-nbnode
  
!   Definition of free ends and other quantities.
  
  CALL intab
  CALL printj (NAME,12)
  GO TO 50
END IF

IF (nbp > 2) THEN
  IF (ial(if1) <= ial(if2)) THEN
    jt = jdiag(nbtr,1)
    jar = arr(nbtr,1)
    jdiag(nbtr,1) = jdiag(nbtr,3)
    arr(nbtr,1) = arr(nbtr,3)
    jdiag(nbtr,3) = jt
    arr(nbtr,3) = jar
    CALL phase (nbtr,jdiag,m4trd)
  END IF
END IF

node = nbtr
i1 = -1
k2 = 3
k3 = 2
GO TO 3

50 RETURN

300 FORMAT ('DIAGRM: Flat graph impossible to build.')

END SUBROUTINE diagrm
!***********************************************************************
!                                                                      *

SUBROUTINE dracah (i,j,k,l,m,n,rac)
!                                                                      *
!   SUBROUTINE  to calculate Racah coefficients. The arguments I, J,   *
!   K, L, M, N should be twice their actual value. Works for integer   *
!   and  half-integer  values of  angular momenta. The routine makes   *
!   use of the GAM  array, thus  SUBROUTINE FACTT must be called be-   *
!   fore this routine is used.                                         *
!                                                                      *
!   Written by N S Scott                    Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN)                      :: i
INTEGER, INTENT(IN)                      :: j
INTEGER, INTENT(IN)                      :: k
INTEGER, INTENT(IN)                      :: l
INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN)                      :: n
doubleprecision, INTENT(OUT)             :: rac
IMPLICIT doubleprecision (a-h, o-z)

INTEGER, PARAMETER :: mfact = 500

COMMON/facts/gam(mfact)

j1 = i+j+m
j2 = k+l+m
j3 = i+k+n
j4 = j+l+n
IF (((2*MAX (i,j,m)-j1) > 0) .OR. (MOD (j1,2) /= 0)) GO TO 2
IF (((2*MAX (k,l,m)-j2) > 0) .OR. (MOD (j2,2) /= 0)) GO TO 2
IF (((2*MAX (i,k,n)-j3) > 0) .OR. (MOD (j3,2) /= 0)) GO TO 2
IF (((2*MAX (j,l,n)-j4) > 0) .OR. (MOD (j4,2) /= 0)) GO TO 2
GO TO 1
2  rac = 0.0D 00
RETURN

1  CONTINUE
j1 = j1/2
j2 = j2/2
j3 = j3/2
j4 = j4/2
j5 = (i+j+k+l)/2
j6 = (i+l+m+n)/2
j7 = (j+k+m+n)/2
numin = MAX (j1,j2,j3,j4)+1
numax = MIN (j5,j6,j7)+1
rac = 1.0D 00
icount = 0

IF (numin == numax) GO TO 4
numin = numin+1

DO  kk = numin,numax
  ki = numax-icount
  rac = 1.0D 00 -(rac*DBLE(ki*(j5-ki+2)*(j6-ki+2)*(j7-ki+2))/  &
      DBLE((ki-1-j1)*(ki-1-j2)*(ki-1-j3)*(ki-1-j4)))
  icount = icount+1
END DO

numin = numin-1
4  rac = rac*((-1.0D 00)**(j5+numin+1)) *EXP( (gam(numin+1)-gam(numin-j1)  &
    -gam(numin  -j2)-gam(numin  -j3)-gam(numin  -j4)-gam(j5+2-numin)  &
    -gam(j6+2-numin)-gam(j7+2-numin))+((gam(j1+1-i)+gam(j1+1-j)  &
    +gam(j1+1-m)-gam(j1+2)+gam(j2+1-k)+gam(j2+1-l)+gam(j2+1-m)  &
    -gam(j2+2)+gam(j3+1-i)+gam(j3+1-k)+gam(j3+1-n)-gam(j3+2)  &
    +gam(j4+1-j)+gam(j4+1-l)+gam(j4+1-n)-gam(j4+2)) *0.5D 00  ))

RETURN
END SUBROUTINE dracah
!***********************************************************************
!                                                                      *

SUBROUTINE gensum (j6c,j7c,j8c,j9c,jwc,j6,j7,j8,j9,jw,jdel,  &
    ldel,sumvar,mp,j6p,j7p,j8p,j9p,jword,nlsum,  &
    nbj,nb6j,k6cp,k7cp,k8cp,k9cp,jsum4,jsum5, jsum6,inv6j,recup)
!                                                                      *
!   Carries  out the summation over coefficients defined by the arr-   *
!   ays J6, J7, J8, LDEL and  JW  to give RECUP. The entry is either   *
!   made from NJGRAF or directly  assuming that the arrays J6,...,JW   *
!   have already been determined  by  a previous entry to NJGRAf and   *
!   that the summation is required for another set of j values defi-   *
!   ned by the array J1. RECUP is the recoupling coefficient.          *
!                                                                      *
!   Call(s) to: [NJGRAF]: DRACAH, RDIAG.                               *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN)                      :: j6c
INTEGER, INTENT(IN)                      :: j7c
INTEGER, INTENT(IN)                      :: j8c
INTEGER, INTENT(IN)                      :: j9c
INTEGER, INTENT(IN)                      :: jwc
INTEGER, INTENT(IN)                      :: j6(m3mngm)
INTEGER, INTENT(IN)                      :: j7(m3mngm)
INTEGER, INTENT(IN)                      :: j8(m3mngm)
INTEGER, INTENT(IN)                      :: j9(mangmp)
INTEGER, INTENT(IN)                      :: jw(6,m6j)
INTEGER, INTENT(IN)                      :: jdel
INTEGER, INTENT(IN)                      :: ldel(m6j,2)
LOGICAL, INTENT(IN OUT)                  :: sumvar(mangm)
INTEGER, INTENT(IN)                      :: mp
INTEGER, INTENT(IN)                      :: j6p(mangmp)
INTEGER, INTENT(IN)                      :: j7p(mangmp)
INTEGER, INTENT(IN)                      :: j8p(mangmp)
INTEGER, INTENT(IN)                      :: j9p(mangmp)
INTEGER, INTENT(IN)                      :: jword(6,m6j)
INTEGER, INTENT(IN OUT)                  :: nlsum
INTEGER, INTENT(IN)                      :: nbj(msum)
INTEGER, INTENT(IN)                      :: nb6j(msum)
INTEGER, INTENT(IN)                      :: k6cp(msum)
INTEGER, INTENT(IN)                      :: k7cp(msum)
INTEGER, INTENT(IN)                      :: k8cp(msum)
INTEGER, INTENT(IN)                      :: k9cp(msum)
INTEGER, INTENT(IN)                      :: jsum4(mtriad,m6j)
INTEGER, INTENT(IN)                      :: jsum5(mtriad,m6j)
INTEGER, INTENT(IN)                      :: jsum6(mtriad)
INTEGER, INTENT(IN OUT)                  :: inv6j(m6j)
doubleprecision, INTENT(OUT)             :: recup
IMPLICIT doubleprecision (a-h, o-z)
LOGICAL :: ldiag,noel,free

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10
INTEGER, PARAMETER :: mfact = 500

doubleprecision, PARAMETER :: epsil = 1.0D-10

DIMENSION mat(mtriad,mtriad),noel(mtriad),  &
    maxlp(mtriad),jsum2(mtriad),jsum3(mtriad),  &
    jsum(2,m6j),jwtest(m6j),wstor(m6j),ipair(2,2),ldiag(mtriad)
DIMENSION xj1(mangm),ist(6)
DIMENSION j12(4,mtriad,mtriad)



COMMON/debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6 /facts/gam(mfact)  &
    /couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)

DATA mxcsvr/4/

!   evaluates all terms in J6, J7, J8, J9, LDEL, JW which do not
!   involve a summation. The result is stored in RECUP and IASTOR

IF (ibug3 == 1) THEN
  
  DO  i = 1,m
    xj1(i) = 0.5D 00*DBLE (j1(i)-1)
  END DO
  
  WRITE (99,400) (xj1(i),i = 1,m)
  WRITE (99,306) nlsum
  WRITE (99,401)
END IF

mm = m+1
j1(mm) = 1

!   Test delta functions

j1(mm) = 1
IF (jdel <= 0) GO TO 180

DO  i = 1,jdel
  i1 = ldel(i,1)
  i2 = ldel(i,2)
  IF ((i1 > mm) .OR. (i2 > mm))THEN
    IF (i1 > mm) j1(i1) = j1(i2)
    IF (i2 > mm) j1(i2) = j1(i1)
  ELSE
    IF (j1(i1) /= j1(i2)) THEN
      recup = 0.0D 00
      RETURN
    END IF
  END IF
END DO

180 recup = 1.0D 00
IF (jwc /= 0) THEN
  
!   Multiply RECUP by all Racah coefficients which do not involve a
!   summation
  
  IF (ibug3 == 1) WRITE (99,309)
  
  DO  i = 1,jwc
    IF (inv6j(i) > 0) CYCLE
    DO  j = 1,6
      i1 = jw(j,i)
      ist(j) = j1(i1) - 1
    END DO
    
    CALL dracah (ist(1),ist(2),ist(3),ist(4),ist(5),ist(6),x1)
    IF (ibug3 == 1) WRITE (99,305) (xj1(jw(k,i)),k = 1,6),x1
    recup = recup*x1
    
  END DO
  
END IF

sqr = 1.0D 00

IF (j6c /= 0) THEN
  DO  i = 1,j6c
    i1 = j6(i)
    sqr = sqr*j1(i1)
  END DO
END IF

spr = 1.0D 00

IF (j9c /= 0) THEN
  DO  i = 1,j9c
    i1 = j9(i)
    spr = spr*j1(i1)
  END DO
END IF

recup = recup*SQRT (sqr/spr)
IF (ABS(recup) < epsil) GO TO 145
iastor = 0

IF (j7c /= 0) THEN
  DO  i = 1,j7c
    i1 = j7(i)
    iastor = iastor + j1(i1) -1
  END DO
END IF

IF (j8c /= 0) THEN
  DO  i = 1,j8c
    i1 = j8(i)
    iastor = iastor +2*(j1(i1)-1)
  END DO
END IF

IF (nlsum <= 0) THEN
  iastor = iastor/2
  
!   No summation involved. End of computation
  
  stor1 = 1.0D 00
  stor = 1.0D 00
  IF (MOD (iastor,2) == 1) recup = -recup
  IF (ibug3 == 1) WRITE (99,303) recup
  RETURN
  
END IF

!   Evaluation of the part involving summations.

nfs = 0
jwr = 0
j6f = 0
j7f = 0
j8f = 0
j9f = 0
nps = 0
25 nps = nps+1
IF (ibug3 == 1) WRITE (99,302) nps

!   Loop on the disconnected summations

ias = 0
nsum = nbj(nps)-nfs
jwrd = nb6j(nps)-jwr
j6cp = k6cp(nps)
j7cp = k7cp(nps)
j8cp = k8cp(nps)
j9cp = k9cp(nps)

!   The range of values of each summation variable is defined by
!   establishing a matrix of the links between variables.
!   MAT(I,J) contains:
!       I = J  Number of possible values of I due to triangular
!              relations with non-variables, i.e. constants.
!       I > J  Number of links between I and J through constants
!       I < J  Value of the constant, if the above is 1. If not,
!              these values are srored in J12(L,I,J) where there
!              is room for MXCSVR such values (L .LE. 4)

DO  i = 1,nsum
  DO  j = 1,nsum
    mat(i,j) = 0
  END DO
END DO

DO  i1 = 1,nsum
  i1t = i1+nfs
  i2 = jsum6(i1t)
  DO  i3 = 1,i2
    i = jsum5(i1t,i3)
    j = jsum4(i1t,i3)
    SELECT CASE ( J )
      CASE (    1)
        GO TO 54
      CASE (    2)
        GO TO 55
      CASE (    3)
        GO TO 56
      CASE (    4)
        GO TO 57
      CASE (    5)
        GO TO 58
      CASE (    6)
        GO TO 59
    END SELECT
    
!   The rows of the IPAIR arrays give limits of summation imposed
    
    54       ipair(1,1) = jword(2,i)
    ipair(1,2) = jword(5,i)
    ipair(2,1) = jword(3,i)
    ipair(2,2) = jword(6,i)
    GO TO 60
    
    55       ipair(1,1) = jword(1,i)
    ipair(1,2) = jword(5,i)
    ipair(2,1) = jword(4,i)
    ipair(2,2) = jword(6,i)
    GO TO 60
    
    56       ipair(1,1) = jword(1,i)
    ipair(1,2) = jword(6,i)
    ipair(2,1) = jword(4,i)
    ipair(2,2) = jword(5,i)
    GO TO 60
    
    57       ipair(1,1) = jword(2,i)
    ipair(1,2) = jword(6,i)
    ipair(2,1) = jword(3,i)
    ipair(2,2) = jword(5,i)
    GO TO 60
    
    58       ipair(1,1) = jword(1,i)
    ipair(1,2) = jword(2,i)
    ipair(2,1) = jword(3,i)
    ipair(2,2) = jword(4,i)
    GO TO 60
    
    59       ipair(1,1) = jword(1,i)
    ipair(1,2) = jword(3,i)
    ipair(2,1) = jword(2,i)
    ipair(2,2) = jword(4,i)
    
    loop63:  60       DO  i4 = 1,2
      km = 0
      DO  i5 = 1,2
        IF (ipair(i4,i5) > mp) km = km+1
      END DO
      
      jj1 = ipair(i4,1)
      jj2 = ipair(i4,2)
      IF (km == 1) GO TO 67
      IF (km > 1) CYCLE loop63
      
!   One variable linked to two constants. Fix the diagonal MAT(I,I)
      
      jt1 = j1(jj1)-1
      jt2 = j1(jj2)-1
      jmin = ABS (jt1-jt2)
      jmax = jt1+jt2
      
      IF (mat(i1,i1) > 1) THEN
        
!   If there are several couples of constants, take the more
!   stringent combination
        
        jmin = MAX (jmin,jsum(1,i1))
        jmax = MIN (jmax,jsum(2,i1))
        IF (jmax >= jmin) THEN
          jsum(1,i1) = jmin
          jsum(2,i1) = jmax
          mat(i1,i1) = (jmax-jmin)/2+1
          CYCLE loop63
        ELSE
          recup = 0.0D 00
          GO TO 110
        END IF
      ELSE IF (mat(i1,i1) < 1) THEN
        
!   First time
        
        mat(i1,i1) = (jmax-jmin)/2+1
        jsum(1,i1) = jmin
        jsum(2,i1) = jmax
      END IF
      
      CYCLE loop63
      
!   One variable linked to one constant and one variable  non diagonal
!   element
      
      67          jt1 = MIN (jj1,jj2)
      jt2 = MAX (jj1,jj2)-mp
      IF (jt2 > i1) CYCLE loop63
      jt4 = j1(jt1)-1
      k = mat(i1,jt2)
      IF (k == 0) GO TO 107
      
      DO  ll = 1,k
        IF (jt4 == j12(ll,jt2,i1)) CYCLE loop63
      END DO
      
      107          k = k+1
      IF (k > mxcsvr) CYCLE loop63
      mat(i1,jt2) = k
      j12(k,jt2,i1) = jt4
      
    END DO loop63
  END DO
END DO

!   Reduce the diagonal elements by taking into account the non
!   diagonal elements, and keep the latter only if needed

150 ichan = 0

DO  i = 1,nsum
  noel(i) = .true.
  i1 = i-1
  IF (i1 == 0) GO TO 170
  DO   j = 1,i1
    IF ((mat(i,j) == 0) .OR. (mat(j,j) == 0)) CYCLE
    ik1 = i
    ik2 = j
    CALL rdiag (i,j,ik1,ik2,ichan,mat,jsum,j12)
    noel(i) = .false.
  END DO
  170    IF (i == nsum) CYCLE
  i2 = i+1
  
  DO  j = i2,nsum
    IF ((mat(j,i) == 0) .OR. (mat(j,j) == 0)) CYCLE
    ik1 = j
    ik2 = i
    CALL rdiag (i,j,ik1,ik2,ichan,mat,jsum,j12)
  END DO
END DO

IF (ichan /= 0) GO TO 150
GO TO 220

!   Carry out the summations.

220 DO  i = 1,nsum
  jsum3(i) = 1
  ldiag(i) = .false.
  IF (mat(i,i) == 1) ldiag(i) = .true.
END DO

DO  i = 1,jwrd
  jwtest(i) = 1
END DO

stor = 0.0D 00
stor1 = 1.0D 00
nolp = 0
ip = 1
240 nolp = nolp+1

!   Find the range of JSUM2(NOLP)
!   NOLP is the index  of the summation variable

jmin = jsum(1,nolp)
jmax = jsum(2,nolp)
IF (noel(nolp)) GO TO 241
no1 = nolp-1

DO  nj = 1,no1
  IF (mat(nolp,nj) == 1) THEN
    jj1 = mat(nj,nolp)
    jj2 = jsum2(nj)
    jmin = MAX (jmin,IABS(jj2-jj1))
    jmax = MIN (jmax,jj1+jj2)
  ELSE IF (mat(nolp,nj) > 1) THEN
    k = mat(nolp,nj)
    jj2 = jsum2(nj)
    
    DO  i = 1,k
      jj1 = j12(i,nj,nolp)
      jmin = MAX (jmin,IABS(jj2-jj1))
      jmax = MIN (jmax,jj1+jj2)
    END DO
    
  END IF
  
END DO

241 jsum2(nolp) = jmin
maxlp(nolp) = jmax
IF (ldiag(nolp)) jsum3(nolp) = 0
IF (nolp < nsum) GO TO 240

DO  jj = jmin,jmax,2
  jsum2(nsum) = jj
  
!   Determine which RACAH coefficients need re-evaluating and
!   set JWTEST appropriately
  
  DO  j = ip,nsum
    IF (jsum3(j) <= 0) CYCLE
    i2 = jsum6(j)
    
    DO  i1 = 1,i2
      i3 = jsum5(j,i1)
      jwtest(i3) = 1
    END DO
  END DO
  
  DO  j = 1,jwrd
    IF (jwtest(j) == 0) CYCLE
    jwj = j+jwr
    
    DO  i = 1,6
      IF (jword(i,jwj) <= mp) THEN
        i1 = jword(i,jwj)
        ist(i) = j1(i1) - 1
      ELSE
        i1 = jword(i,jwj)-mp-nfs
        ist(i) = jsum2(i1)
      END IF
    END DO
    
    CALL dracah (ist(1),ist(2),ist(3),ist(4),ist(5),ist(6),x1)
    wstor(j) = x1
    IF (ibug3 == 1) THEN
      DO  i = 1,6
        xj1(i) = 0.5D 00*DBLE (ist(i))
      END DO
      
      WRITE (99,305) (xj1(i), i = 1,6),x1
    END IF
  END DO
  
!   Form product of Racah coefficients, (2J+1) factors and (-1)
!   factors in STOR1
  
  DO  i = 1,jwrd
    stor1 = stor1*wstor(i)
  END DO
  
!   IASTOR contains the power of (-1) which is common to all terms
  
  ix2 = 0
  ij6cp = 1
  IF (j6cp /= j6f) THEN
    jb = j6f+1
    
    DO  i = jb,j6cp
      i1 = j6p(i)-nfs
      ij6cp = ij6cp*(jsum2(i1)+1)
    END DO
  END IF
  
  IF (j9cp /= j9f) THEN
    jb = j9f+1
    
    DO  i = jb,j9cp
      i1 = j9p(i)-nfs
      ij6cp = ij6cp/(jsum2(i1)+1)
    END DO
  END IF
  
  stor1 = stor1*SQRT (DBLE (ij6cp))
  
  IF (j7cp /= j7f) THEN
    jb = j7f+1
    
    DO  i = jb,j7cp
      i1 = j7p(i)-nfs
      ix2 = ix2 + jsum2(i1)
    END DO
  END IF
  
  IF (j8cp /= j8f) THEN
    jb = j8f+1
    
    DO  i = jb,j8cp
      i1 = j8p(i)-nfs
      ix2 = ix2 + 2*(jsum2(i1))
    END DO
  END IF
  
  IF (MOD(ix2,2) == 1) THEN
    ias = -1
    ix2 = ix2+1
  END IF
  
  ix2 = ix2/2
  
!   Add term into STOR and reset STOR1 to 1 ready for next term
  
  IF (MOD(ix2,2) == 1) stor1 = -stor1
  stor = stor + stor1
  stor1 = 1.0D 00
  nsum1 = nsum-1
  IF (nsum1 == 0) CYCLE
  
  DO  ik = 1,nsum1
    jsum3(ik) = 0
  END DO
  
  DO  ik = 1,jwrd
    jwtest(ik) = 0
  END DO
  
END DO

250 nolp = nolp-1

IF (nolp /= 0) THEN
  IF (ldiag(nolp)) GO TO 250
  jsum3(nolp) = 1
  jsum2(nolp) = jsum2(nolp)+2
  IF (jsum2(nolp) > maxlp(nolp)) GO TO 250
  ip = nolp
  
!   Proceed to next variable
  
  GO TO 240
  
END IF

recup = recup*stor
IF (ibug3 == 1) WRITE (99,307) nps,stor,recup
IF (ABS(recup) < epsil) GO TO 145
jwr = jwrd+jwr
nfs = nsum+nfs
j6f = j6cp
j7f = j7cp
j8f = j8cp
j9f = j9cp
iastor = iastor+ias

!   Proceed to next sum

IF (nps < nlsum) GO TO 25
iastor = iastor/2
IF (MOD (iastor,2) /= 0) recup = -recup
IF (ibug3 == 1) WRITE (99,304) recup
110 RETURN

!   No summations. Check that there are no inconsistencies. Then
!   multiply by (-1) factor and exit

145 recup = 0.0D 00
RETURN

302 FORMAT (' Sum Nr.',i3)
303 FORMAT (' No summation. Recoupling coefficient = ',g15.8)
304 FORMAT (' Recoupling coefficient = ',g15.8)
305 FORMAT (6F5.1,10X,g15.8)
306 FORMAT (' Number of independent sums:',i3)
307 FORMAT (' Sum Nr.',i2,' Sum value = ',g15.8,' RECUP = ',g15.8)
309 FORMAT (' Not involving summation variable')
400 FORMAT (//' Printout from SUBROUTINE GENSUM'  &
    //' Values of angular momenta in *REAL* FORMAT' /(14F5.1))
401 FORMAT (/' Racah W functions (6J)'  &
    /' Arguments in *REAL* FORMAT',18X,'value')

END SUBROUTINE gensum
!***********************************************************************
!                                                                      *

SUBROUTINE intab
!                                                                      *
!   This SUBROUTINE called at the end of DIAGRM, fixes the arrays IH   *
!   and IL - so to speak hardware and logical addresses of triads in   *
!   JDIAG . Also  determines the number of free ends NFREE and their   *
!   location ITFREE.                                                   *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************

IMPLICIT doubleprecision (a-h, o-z)
INTEGER :: arr,tab1
LOGICAL :: free

PARAMETER ( mangm = 60,m3mngm = 3*mangm,  &
    mtriad = 12,m2trd = 2*mtriad,m4trd = 4*mtriad, m6j = 20,msum = 10)

COMMON/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc  &
    /build/ial(m4trd),if1,if2,node  &
    /couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)

DO  i = 1,m
  ial(i) = 1
END DO

DO  i = ifirst,ilast
  j = jdiag(i,1)
  k = ial(j)
  tab1(j,k) = i
  ial(j) = k+1
END DO

ifr = ifirst-1

DO  i = ifirst,ilast
  it = i-ifr
  il(i) = it
  ih(it) = i
END DO

j = jdiag(ifirst,3)
k = ial(j)
IF (k > 1) tab1(j,2) = tab1(j,1)
tab1(j,1) = ifirst
ial(j) = 3
j = jdiag(ilast,2)
tab1(j,2) = ilast
ial(j) = 3
nfree = 0

DO  i = ifirst,ilast
  j = jdiag(i,1)
  IF (ial(j) /= 3) THEN
    nfree = nfree+1
    itt = ilast+nfree
    tab1(j,2) = itt
    il(itt) = nfree*1000
    itfree(nfree) = i
  END IF
END DO

RETURN
END SUBROUTINE intab
!***********************************************************************
!                                                                      *

SUBROUTINE lolpop (fail)
!                                                                      *
!   Reduces a loop with one line and one node in the flat graph.       *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


LOGICAL, INTENT(IN OUT)                  :: fail
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME,namsub
INTEGER :: arr,tab1
LOGICAL :: sumvar,free

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

DIMENSION kp(3),ks(3)

COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3), free(mangm)  &
    /graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc /nam/namsub

DATA NAME/'LOLPOP'/
DATA kp/2,3,1/
DATA ks/0,1,-1/

namsub = NAME
i1 = npoint(1)
k3 = 2
IF (i1 == ilast) k3 = 3
l = jdiag(i1,k3)
CALL delta (l,mp,fail)
IF (fail) RETURN
k = kp(k3)
IF (arr(i1,k) < 0) CALL phase2 (jdiag(i1,k))
k1 = ks(k3)
il1 = il(i1)+k1
i2 = ih(il1)
l1 = jdiag(i2,1)
CALL delta (l1,jdiag(i2,k3),fail)
IF (fail) RETURN
IF (arr(i2,k3) == k1) CALL phase2 (l1)
il2 = il(i2)+k1
i3 = ih(il2)
k2 = k3+k1
jdiag(i3,k2) = l1
arr(i3,k2) = arr(i2,1)
j9c = j9c+1
j9(j9c) = l1
j6c = j6c+1
j6(j6c) = jdiag(i1,1)
IF (k3 == 3) RETURN

DO  i = 3,nbnode
  it = ih(i)
  ilp = i-2
  il(it) = ilp
  ih(ilp) = it
END DO

RETURN
END SUBROUTINE lolpop
!***********************************************************************
!                                                                      *

SUBROUTINE neibor (lc,l1,l2)
!                                                                      *
!   Gives the positions of the other two arguments in the triad.       *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: lc
INTEGER, INTENT(OUT)                     :: l1
INTEGER, INTENT(OUT)                     :: l2
IMPLICIT doubleprecision (a-h, o-z)

IF (lc < 2) THEN
  l1 = 2
  l2 = 3
ELSE IF (lc == 2) THEN
  l1 = 3
  l2 = 1
ELSE
  l1 = 1
  l2 = 2
END IF

RETURN
END SUBROUTINE neibor
!***********************************************************************
!                                                                      *

SUBROUTINE njgraf (recup,igen,fail)
!                                                                      *
!   Program to calculate a general recoupling coefficient. This ver-   *
!   sion is slightly modified  (Colin T Johnson, Oxford University).   *
!   The changes are as follows:                                        *
!                                                                      *
!      1. PARAMETER IGEN has been included in the argument list:       *
!            IGEN =  0  normal call to NJGRAF                          *
!            IGEN = -1  GENSUM is not called                           *
!      2. The contents of COMMON blocks /ARGU/ and /SUMARG/ are used   *
!         by  GENSUM  to calculate the recoupling coefficient. These   *
!         COMMON  blocks  have  been removed from GENSUM. Their con-   *
!         tents are passed to  GENSUM  through the argument list in-   *
!         stead, so that NJGRAF can be called to set up formulae for   *
!         both the direct and exchange cases in COR and BREIT.         *
!      3. Extra  dimension  tests  have  been  included  in routines   *
!         NJGRAF, PRINTJ, SPRATE, VAR and ZERO. These are  discussed   *
!         below.                                                       *
!      4. An  extra routine  RDIAG  has been introduced to remove an   *
!         extended  DO loop from GENSUM, to conform with the FORTRAN   *
!         77 standard.                                                 *
!                                                                      *
!                                                                      *
!   Description of some COMMON blocks. A full discussion is given in   *
!   the NJGRAF program description (Bar-Shalom and Klapisch op cit).   *
!                                                                      *
!      COMMON block COUPLE                                             *
!                                                                      *
!         M                The total number of angular momentum val-   *
!                          ues in the initial and final states         *
!         N                The number of basic angular momentum val-   *
!                          ues that are coupled                        *
!         J1(I),           The angular momentum values stored as 2J+1  *
!            I = 1,M                                                   *
!         J2(I,J),         The position in the J1 array of the init-   *
!            I = 1,(N-1),  ial state triads                            *
!            J = 1,3                                                   *
!         J3(I,J),         The position in the J1 array of the final   *
!            I = 1,(N-1),  state triads                                *
!            J = 1,3                                                   *
!         FREE(I),         If FREE(I) = .TRUE., no reference is made   *
!            I = 1,M       to the value of J1(I) when establishing a   *
!                          formula in  NJGRAF .  GENSUM  may then be   *
!                          called  for  repeated  occurences of this   *
!                          formula  with  differing values of J1(I).   *
!                          If J1(I) does  not  vary between calls to   *
!                          GENSUM then FREE(I) should be set .FALSE.   *
!                          so that zero branches  can be removed be-   *
!                          fore the formula is established.            *
!                                                                      *
!      COMMON block DEBUG                                              *
!                                                                      *
!         IBUG1            Not used                                    *
!         IBUG2            Not used                                    *
!         IBUG3            Debug prints in NJGRAF and GENSUM if 1      *
!         IBUG4            Not used                                    *
!         IBUG5            Not used                                    *
!         IBUG6            Not used                                    *
!                                                                      *
!      COMMON block ARGU                                               *
!                                                                      *
!         J6C              The number of elements in the K6 array      *
!         J7C              The number of elements in the K7 array      *
!         J8C              The number of elements in the K8 array      *
!         J9C              The number of elements in the K9 array      *
!         JWC              The number of columns in the KW array       *
!         J6(I),           Each entry corresponds to a  factor  SQRT   *
!            I = 1,J6C     (2J+1) in RECUP. The value  of  J6  GIVES   *
!                          position in  J1  array where  J  value is   *
!                          found                                       *
!         J7(I),           Each entry corresponds to factor  (-1)**J   *
!            I = 1,J7C     in RECUP                                    *
!         J8(I),           Each entry corresponds to a factor (-1)**   *
!            I = 1,J8C     (2J) in RECUP                               *
!         J9(I),           Each entry corresponds to a factor (2J+1)   *
!            I = 1,J9C     **(-0.5) in RECUP                           *
!         KW(I,J),         Each column corresponds to a Racah coeff-   *
!            I = 1,6,      icient in RECUP                             *
!            J = 1,JWC                                                 *
!         JDEL             The number of delta functions               *
!         LDEL(I,J),       The arguments of the delta functions        *
!              J = 1,2                                                 *
!         SUMVAR(I)        .TRUE. for ang. mom. I (a summation vari-   *
!                          able                                        *
!         MP               The index of the last variable              *
!                                                                      *
!   The arrays  J6, J7, J8, J9 and  KW, Are evaluated by NJGRAF. The   *
!   summation over the variables in  J6, J7, J8, J9 and  KW, and the   *
!   evaluation of RECUP is carried out in GENSUM. GENSUM  can be re-   *
!   entered directly to evaluate different  recoupling  coefficients   *
!   with the same structure  by just  altering the numbers in the J1   *
!   array.                                                             *
!                                                                      *
!   This is the main program. It handles all the analysis of the re-   *
!   coupling  coefficient without referring explicitly to the values   *
!   of angular  momenta  which  are in J1(J),except for zero in case   *
!   FREE = .FALSE. . Like NJSYM it  prepares arrays of arguments for   *
!   phase factors, (2*J+1) factors and  6j-coefficients to be compu-   *
!   ted in GENSUM, which can also be called separately when only the   *
!   numerical values of angular momenta change. These variable angu-   *
!   lar momenta should be declared  FREE(J)  = .TRUE. , so  that the   *
!   formula prepared for GENSUM should be correct when J1 is not ze-   *
!   ro. FAIL will be TRUE when the  recoupling  coefficient  is zero   *
!   because of unsatisfied delta or other similar causes.              *
!                                                                      *
!   This version holds the array dimensions in parameter statements.   *
!   The dimensions are labelled:                                       *
!                                                                      *
!      MANGM  : Dimension of the J1 and FREE arrays in /COUPLE/, and   *
!               the  first  dimension of the LINE and LCOL arrays in   *
!               /TREE/. Also  the  dimension  of the SUMVAR array in   *
!               /ARGU/, AND OF THE INVER array in routine SPRATE. It   *
!               is tested for  M  on entry to  NJGRAF, and for MP in   *
!               routine SPRATE.                                        *
!      MTRIAD : Dimension of the  J2 and  J3 arrays in /COUPLE/. The   *
!               dimensions of these  arrays  are checked on entry to   *
!               NJGRAF in addition  MTRIAD sets the dimension of the   *
!               JSUM6 array and the first dimension of the JSUM4 and   *
!               JSUM5  arrays in /SUMARG/. Also gives the dimensions   *
!               of some  temporary working arrays in SPRATE and GEN-   *
!               SUM. In these  cases  mtriad sets the maximum number   *
!               of summation variables  in any particular sum, which   *
!               is tested in SPRATE.                                   *
!      M2TRD  : (=2*MTRIAD) Dimension of the J23 ,  ARROW  and  TABS   *
!               arrays in /TREE/. Also  the  dimension of the npoint   *
!               array in /GRAPH/.                                      *
!      M4TRD  : (=4*MTRIAD) Dimension of the  JDIAG,  ARR, IL and IH   *
!               arrays in /GRAPH/, and of the IAL array in /BUILD/.    *
!      M3MNGM : Dimension of the J6 array in /ARGU/, tested in SPRATE  *
!               Dimension of the J7 array in /ARGU/, tested in SPRATE  *
!               Dimension of the J8 array in /ARGU/, tested in SPRATE  *
!      MANGMP : Dimension of the J9 array in /ARGU/, tested in SPRATE  *
!               MANGMP also sets the dimension of the J6P,  J7P, J8P   *
!               and J9P arrays in /SUMARG/, And of the JNS  array in   *
!               routine VAR. The dimension of the JNS array is  tes-   *
!               ted in VAR.                                            *
!      M6J    : Dimension of the JW(or KW) and LDEL arrays in /ARGU/,  *
!               and of the JWORD and INV6J arrays in /SUMARG/.  Also   *
!               the second dimension of the  JSUM4 and  JSUM5 arrays   *
!               in /SUMARG/. In addition it  gives the dimensions of   *
!               a  number of  temporary  working  arrays in routines   *
!               SPRATE and GENSUM. M6J is tested in SPRATE.            *
!      MFACT  : The dimension of the factorial array GAM in /FACTS /.  *
!      MSUM   : Dimension of the NBJ, NB6J, K6CP, K7CP, K8CP and K9CP  *
!               arrays in /SUMARG/. MSUM is the  maximum  number  of   *
!               sums allowed, and is tested in routine SPRATE.         *
!      MTAB   : The dimension of the JTAB array in  routine  PRINTJ.   *
!               MTAB is tested in PRINTJ.                              *
!      MZERO  : Dimension of the JZERO array in /ZER/. MZERO is tes-   *
!               ted in routine ZERO.                                   *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


doubleprecision, INTENT(OUT)             :: recup
INTEGER, INTENT(IN OUT)                  :: igen
LOGICAL, INTENT(OUT)                     :: fail
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME,namsub
LOGICAL :: find,tabs,cut,free,sumvar
INTEGER :: arrow,arr,tab1

INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mfact = 500
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /cons/zro,half,tenth,one,two,three,ten  &
    /couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm) /cutdig/cut  &
    /debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6 /facts /gam(mfact)  &
    /graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc /nam/namsub  &
    /sumarg/j6p(mangmp),j7p(mangmp),j8p(mangmp),j9p(mangmp),  &
    jword(6,m6j),nlsum,nbj(msum),nb6j(msum),k6cp(msum),  &
    k7cp(msum),k8cp(msum),k9cp(msum),jsum6(mtriad),  &
    jsum4(mtriad,m6j),jsum5(mtriad,m6j),inv6j(m6j)  &
    /tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),  &
    lcol(mangm,2),tabs(m2trd),nbtr

DATA NAME/'NJGRAF'/

!   Debug printout

IF (ibug3 == 1) THEN
  WRITE (99,300)
  WRITE (99,301) m,n-1
  WRITE (99,302)
  WRITE (99,303) (j1(i),i = 1,m)
  WRITE (99,304) (free(i),i = 1,m)
  WRITE (99,305)
  WRITE (99,306) ((j2(i,j),j = 1,3), (j3(i,j),j = 1,3),i = 1,n-1)
END IF

!   Test the dimension of the J1 array

IF (m+1 > mangm) THEN
  WRITE (*,307)
  WRITE (*,308) m+1,mangm
  STOP
END IF

!   Test the dimensions of the J2 and J3 arrays

IF (n-1 > mtriad) THEN
  WRITE (*,307)
  WRITE (*,309) n-1,mtriad
  STOP
END IF

!   Initializations

DO  i = n,mtriad
  DO  j = 1,3
    j2(i,j) = 0
    j3(i,j) = 0
  END DO
END DO

fail = .false.
j6c = 0
j7c = 0
j8c = 0
j9c = 0
jwc = 0
jdel = 0
CALL setdm
nfin = 0
cut = .false.

!   Building up of the unstructured graph

CALL settab (fail)

!   Exit with RECUP set to zero if any angular momentum is
!   impossible

m = m+1
IF (fail) GO TO 7

m = m-1

!   Locate and eliminate any zero angular momentum; simplify the
!   graph

jf = 0
jf1 = 0
CALL zero (jf1,jf,fail)
IF (fail) GO TO 7

mp = m
IF (nbtr == 0) GO TO 6
jump = 1

1 CALL chklp1 (fail)
IF (fail) GO TO 7

!   Build a flat diagram out of the unstructured graph; several flat
!   diagrams may constitute the original graph, in which case there
!   are possible cuts; the flat diagrams will have free ends if cut

CALL diagrm (jump)
nfin = MAX (0,nfree-2)

IF (nfin /= 0) THEN
  jump = 3
  
!   Handling of free ends if a cut was found
  
  CALL cutnl (fail)
  IF (fail) GO TO 7
ELSE
  jump = 2
  IF (nfree == 1) THEN
    CALL cut1l (fail)
    IF (fail) GO TO 7
  ELSE IF (nfree > 1) THEN
    CALL cut2l (fail)
    IF (fail) GO TO 7
  END IF
END IF

nbtr = nbtr+nfin
IF (nbtr /= 0) cut = .true.

!   Analysis of the flat diagram.
!   Closed circuits of increasing order NC are searched, analysed,
!   and taken out of the flat diagram, thus reducing the number of
!   nodes, NBNODE.

nc = 0
10 nc = nc+1
CALL search (find)
IF (.NOT. find) GO TO 10
ncp = nc-2
jpol = 0
IF ((m == mp) .AND. (nc > 3)) CALL setdm
IF (ipartl > 2) CALL polygn (jpol)
SELECT CASE ( nc )
  CASE (    1)
    GO TO 11
  CASE (    2)
    GO TO 12
  CASE (    3)
    GO TO 13
  CASE (    4)
    GO TO 14
END SELECT
11 CALL lolpop (fail)
IF (fail) GO TO 7
GO TO 15
12 CALL bubble (jpol,fail)
IF (fail) GO TO 7
GO TO 15
13 CALL triang (fail)
IF (fail) GO TO 7
GO TO 15
14 CALL square
15 nbnode = nbnode-2
IF (nbnode == 0) GO TO 9
ifirst = ih(1)
ilast = ih(nbnode)

!   PRINTJ is an all purpose printing SUBROUTINE called from many
!   places

CALL printj (namsub,8)
IF (nbnode == nfin) GO TO 9
nc = ncp

!   Proceed to other circuits of order NC-1

GO TO 10
9 IF (nbtr == 0) GO TO 6
IF (jump == 3) CALL ordtri

!   At this stage, the flat diagram has been reduced to nodes
!   involving free ends. Proceed to build other flat diagrams
!   if necessary.

GO TO 1

!   All parts of the original graph have been reduced.

7 recup = zro
m = m-1
RETURN
6 CALL printj (NAME,0)

!   Preparation of the results, and separation in several sums
!   if cuts have been detected, also in the flat diagram itself

CALL sprate (m)
m = m-1

!   GENSUM computes the numerical value of the recoupling
!   coefficient

IF (igen /= -1) CALL gensum (j6c,j7c,j8c,j9c,jwc,j6,j7,j8,j9,kw,  &
    jdel,ldel,sumvar,mp,j6p,j7p,j8p, j9p,jword,nlsum,nbj,nb6j,k6cp,k7cp,  &
    k8cp,k9cp,jsum4,jsum5,jsum6,inv6j, recup)

RETURN

300 FORMAT (//' ++++++++++ NJGRAF ++++++++++'/)
301 FORMAT (' Total number of angular momenta (M) = ',1I3  &
    //' Number of triads in each of the ''LEFT-HAND'' AND',  &
    ' ''RIGHT-HAND'' STATES (N-1) = ',1I3)
302 FORMAT (/' (2J+1)-value for each angular momentum:')
303 FORMAT (1X,42I3)
304 FORMAT (1X,42L3)
305 FORMAT (/' ''LEFT-HAND'' TRIADS',10X,'''RIGHT-HAND'' TRIADS')
306 FORMAT (1X,3I3,19X,3I3)
307 FORMAT (/' ***** Error in NJGRAF *****'/)
308 FORMAT (' M+1 = ',1I3,', exceeds PARAMETER MANGM = ',1I3)
309 FORMAT (' N-1 = ',1I3,', exceeds PARAMETER MTRIAD = ',1I3)

END SUBROUTINE njgraf
!***********************************************************************
!                                                                      *

SUBROUTINE ordtri
!                                                                      *
!   This subroutine orders the triads which were left with free ends   *
!   as consequence of cutting,so that the new graph will start there.  *
!                                                                      *
!                                           Last update: 16 Aug 1992   *
!                                                                      *
!***********************************************************************

IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME
INTEGER :: arrow,arr,tab1
LOGICAL :: sumvar,tabs

PARAMETER ( mangm = 60,m3mngm = 3*mangm,mangmp = 2*(mangm/3),  &
    mtriad = 12,m2trd = 2*mtriad,m4trd = 4*mtriad, m6j = 20,msum = 10)

COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),  &
    lcol(mangm,2),tabs(m2trd),nbtr  &
    /graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc  &
    /build/ial(m4trd),if1,if2,node /KEEP/jkp(2,3),jarr(2,3),it2,it3,it5

DATA NAME/'ORDTRI'/

DO  i = 1,mp
  ial(i) = 0
END DO

IF (nfin /= 0) THEN
  nbt1 = nbtr-1
  nbt = nbt1+nfin
  nbtt = nbt+1
  nb = 0
  GO TO 31
END IF

nf = nbtr-itfree(1)

IF (it5 == 0) THEN
  nbt1 = nbtr-1
  n0 = 0
  nft = nfree
  isw = 2
  GO TO 100
END IF

nft = it5-it2
nm = nft+nbtr+1
nbt1 = nbtr

DO  j = 1,3
  jdiag(nbtr,j) = jkp(1,j)
  arr(nbtr,j) = jarr(1,j)
END DO

jt = jdiag(nm,1)
n0 = 0
isw = 1
GO TO 100

22 n0 = nft

DO  j = 1,3
  jdiag(nm,j) = jkp(2,j)
  arr(nm,j) = jarr(2,j)
END DO

nbt1 = nbt1+1
nft = it3-it5
isw = 3
GO TO 100

24 nbt1 = k-nft

23 node = nbt1+nft
CALL change (node,2)
GO TO 40

31 DO  i = 1,nbnode
  i1 = ih(i)
  IF (il(i1) > ilast) CYCLE
  i2 = nbt1+i
  IF (i1 > nbtt) GO TO 33
  IF (i1 == i2) GO TO 32
  IF (il(i2) <= nbnode) CYCLE
  
  33    DO  j = 1,3
    jdiag(i2,j) = jdiag(i1,j)
    arr(i2,j) = arr(i1,j)
  END DO
  
  il(i1) = ilast+i
  32    nb = nb+1
  il(i2) = 0
  
END DO

IF (nb /= nfin) GO TO 31
node = nbt
40 if1 = jdiag(nbtr,1)
if2 = jdiag(nbtr,3)

DO  i = nbtr,node
  DO  k = 1,3
    j = jdiag(i,k)
    ial(j) = ial(j)+1
  END DO
END DO

ilast = node
CALL printj (NAME,8)

RETURN

100 IF (nf <= 0) THEN
  nfr = n0
  i1 = 1
ELSE
  nfr = nft+1
  i1 = -1
END IF

DO  i = 1,nft
  ik = nfr+i1*i
  it = itfree(ik)
  k = nbt1+ik
  
  DO  j = 1,3
    jdiag(k,j) = jdiag(it,j)
    arr(k,j) = arr(it,j)
  END DO
  
END DO

SELECT CASE ( isw )
  CASE (    1)
    GO TO 22
  CASE (    2)
    GO TO 23
  CASE (    3)
    GO TO 24
END SELECT

END SUBROUTINE ordtri
!***********************************************************************
!                                                                      *

SUBROUTINE otherj (lin,j,lo,lco,k)
!                                                                      *
!   Gives the other triad where a given J occurs and its position.     *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: lin
INTEGER, INTENT(IN OUT)                  :: j
INTEGER, INTENT(OUT)                     :: lo
INTEGER, INTENT(OUT)                     :: lco
INTEGER, INTENT(OUT)                     :: k
IMPLICIT doubleprecision (a-h, o-z)
INTEGER :: arrow
LOGICAL :: tabs

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),  &
    lcol(mangm,2),tabs(m2trd),nbtr

lo = line(j,1)
IF ((lo == lin) .OR. (tabs(lo))) THEN
  k = 1
  lo = line(j,2)
  lco = lcol(j,2)
ELSE
  k = 2
  lco = lcol(j,1)
END IF

RETURN
END SUBROUTINE otherj
!***********************************************************************
!                                                                      *

SUBROUTINE phase (l,jm,ndim)
!                                                                      *
!   Phase factor arising from non-cyclic permutation of arguments in   *
!   triad L. JM may be either J23 or JDIAG.                            *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: l
INTEGER, INTENT(IN)                      :: jm(ndim,3)
INTEGER, INTENT(IN OUT)                  :: ndim
IMPLICIT doubleprecision (a-h, o-z)
LOGICAL :: sumvar

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10



COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp

j7(j7c+1) = jm(l,1)
j7(j7c+2) = jm(l,2)
j7c = j7c+3
j7(j7c) = jm(l,3)

RETURN
END SUBROUTINE phase
!***********************************************************************
!                                                                      *

SUBROUTINE phase2 (j)
!                                                                      *
!   Adds a phase factor (-1)**2J .                                     *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN)                      :: j
IMPLICIT doubleprecision (a-h, o-z)
LOGICAL :: sumvar

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp

j8c = j8c+1
j8(j8c) = j

RETURN
END SUBROUTINE phase2
!***********************************************************************
!                                                                      *

SUBROUTINE polygn (jpol)
!                                                                      *
!   This routine reduces a circuit of arbitrary order NC. It exchan-   *
!   ges nodes on the flat diagram until the distance on the axis be-   *
!   tween nodes equals one. Each exchange introduces a summation va-   *
!   riable  and  a 6j-symbol. The circuit has a maximum of NPART = 2   *
!   disconnected parts on the axis.                                    *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(OUT)                     :: jpol
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME
INTEGER :: arr,tab1
LOGICAL :: sumvar

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc  &
    /argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp

DATA NAME/'POLYGN'/

nc1 = nc+1
nc2 = nc
nbc = ipartl-2

10 DO  i = 1,nbc
  it2 = npoint(nc1-i)
  it1 = npoint(nc2-i)
  jb = jdiag(it1,1)
  jc = jdiag(it2,1)
  jdiag(it1,1) = jc
  jdiag(it2,1) = jb
  jar = arr(it1,1)
  arr(it1,1) = arr(it2,1)
  arr(it2,1) = jar
  je = jdiag(it1,2)
  mp = mp+1
  sumvar(mp) = .true.
  jdiag(it1,2) = mp
  jdiag(it2,3) = mp
  
  IF (tab1(jb,1) == it1) THEN
    tab1(jb,1) = it2
  ELSE
    tab1(jb,2) = it2
  END IF
  
  IF (tab1(jc,1) == it2) THEN
    tab1(jc,1) = it1
  ELSE
    tab1(jc,2) = it1
  END IF
  
  IF (arr(it1,2) <= 0) THEN
    CALL phase2 (je)
    arr(it1,2) = 1
    arr(it2,3) = -1
  END IF
  
  jwc = jwc+1
  kw(1,jwc) = jb
  kw(2,jwc) = mp
  kw(3,jwc) = je
  kw(4,jwc) = jc
  kw(5,jwc) = jdiag(it2,2)
  kw(6,jwc) = jdiag(it1,3)
  j6(j6c+1) = mp
  j6c = j6c+2
  j6(j6c) = mp
END DO

nc = nc-nbc

IF (nc > 4) THEN
  nbc = iparts-2
  nc1 = iparts+1
  nc2 = iparts
  GO TO 10
END IF

IF (npart /= 1) THEN
  npoint(3) = npoint(nc1)
  npoint(4) = npoint(nc1+1)
END IF

IF (nc == 2) jpol = 1
CALL printj (NAME,10)

RETURN
END SUBROUTINE polygn
!***********************************************************************
!                                                                      *

SUBROUTINE printj (names,jp)
!                                                                      *
!   This  SUBROUTINE  prints  intermediate  results in standard form   *
!   from wherever it is called.                                        *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


CHARACTER (LEN=6), INTENT(IN OUT)        :: names
INTEGER, INTENT(IN)                      :: jp
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=8) :: iblank,ifree,ifr
CHARACTER (LEN=6) :: nsettb
CHARACTER (LEN=4) :: i6,i7,i8,i9,ij1
CHARACTER (LEN=1) :: im,ip,is(3)
INTEGER :: arr,tab1,arrow
LOGICAL :: tabs,sumvar,free

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10
INTEGER, PARAMETER :: mtab = 30
INTEGER, PARAMETER :: mzero = 20

DIMENSION ix(7),jtab(mtab,3)

COMMON/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),  &
    lcol(mangm,2),tabs(m2trd),nbtr  &
    /argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /const/i6c,i7c,i8c,i9c,idel,iwc /zer/nzero,jzero(mzero)  &
    /debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6  &
    /couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)  &
    /graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc

EQUIVALENCE(i6c,ix(1))

DATA iblank,ifree,ip,im/'        ','FREE END','+','-'/
DATA nsettb/'SETTAB'/

DATA i6,i7,i8,i9,ij1/'I6=','I7=','I8=','I9=','J1='/

IF (ibug3 /= 1) RETURN
WRITE (99,1050) names

!   Initialise variables

i6c = 1
i7c = 1
i8c = 1
i9c = 1
idel = 1
iwc = 1

jump = jp
IF (jump == 0) THEN
  
  DO  i = 1,7
    ix(i) = 1
  END DO
  
  WRITE (99,1020) ij1,(j1(i),i = 1,m)
END IF

IF (jump < 8) GO TO 20
WRITE (99,1000) nbnode,nbtr,nfin,ifirst,ilast,nfree
jump = jump-8
WRITE (99,1001)
k = 0

DO  i = 1,nbnode
  it = ih(i)
  ifr = iblank
  jt = jdiag(it,1)
  
  IF ((tab1(jt,2) /= it) .OR. (jt == jdiag(ifirst,3))) THEN
    k = k+1
    IF (k > mtab) THEN
      WRITE (*,100) k,mtab
      STOP
    END IF
    jtab(k,1) = jt
    jtab(k,2) = tab1(jt,1)
    jtab(k,3) = tab1(jt,2)
  END IF
  
  IF (tab1(jt,2) > ilast) ifr = ifree
  
  DO  j = 1,3
    is(j) = ip
    IF (arr(it,j) < 1) is(j) = im
  END DO
  
  WRITE (99,1002) (is(j),j = 1,3)
  WRITE (99,1003) il(it),it,ifr,(jdiag(it,j),j = 1,3)
  
END DO

WRITE (99,1004)
ntime = 0
jt = jdiag(ifirst,3)
IF (jt /= jdiag(ilast,2)) THEN
  IF (tab1(jt,2) < 1000) GO TO 5
END IF
4 k = k+1
IF (k > mtab) THEN
  WRITE (*,101) k,mtab
  STOP
END IF
jtab(k,1) = jt
jtab(k,2) = tab1(jt,1)
jtab(k,3) = tab1(jt,2)
5 ntime = ntime+1

IF (ntime /= 2) THEN
  jt = jdiag(ilast,2)
  IF (tab1(jt,2) == 1000) GO TO 4
END IF

WRITE (99,1005) ((jtab(i,j),j = 1,3),i = 1,k)
WRITE (99,1006) (i,sumvar(i),i = 1,mp)
20 IF (jump < 4) GO TO 30
jump = jump-4
nbtr1 = 2*n-2
WRITE (99,1010) nbtr1
k = 0

DO  i = 1,nbtr1
  IF (tabs(i)) CYCLE
  k = k+1
  
  DO  j = 1,3
    is(j) = ip
    IF (arrow(i,j) < 1) is(j) = im
  END DO
  
  WRITE (99,1012) (is(j),j = 1,3)
  WRITE (99,1013) k,i,(j23(i,j),j = 1,3)
  
END DO

WRITE (99,1014)
mm = m
IF (names /= nsettb) mm = m-1
WRITE (99,1015) (i,(line(i,j),lcol(i,j),j = 1,2),i = 1,mm)

30 IF (jump >= 2) THEN
  jump = jump-2
  WRITE (99,1030) nc,npart,ipartl,iparts,icross, (npoint(i),i = 1,nc)
END IF

IF (jump >= 1) WRITE (99,1040) nzero,(i,jzero(i),i = 1,nzero)
IF (j6c >= i6c) WRITE (99,1020) i6,(j6(i),i = i6c,j6c)
IF (j7c >= i7c) WRITE (99,1020) i7,(j7(i),i = i7c,j7c)
IF (j8c >= i8c) WRITE (99,1020) i8,(j8(i),i = i8c,j8c)
IF (j9c >= i9c) WRITE (99,1020) i9,(j9(i),i = i9c,j9c)
IF (jdel >= idel) WRITE (99,1021) ((ldel(i,j),j = 1,2),i = idel,jdel)
IF (jwc >= iwc) WRITE (99,1022) ((kw(j,i),j = 1,6),i = iwc,jwc)
i6c = j6c+1
i7c = j7c+1
i8c = j8c+1
i9c = j9c+1
idel = jdel+1
iwc = jwc+1
RETURN

100 FORMAT (' Dimension error in PRINTJ. K = ',i5,' MTAB = ',i5)
101 FORMAT (' Dimension error IN PRINTJ. K = ',i5,' MTAB = ',i5)
1000 FORMAT (/10X,'NBNODE = ',i3,10X,'NBTR = ',i3,10X,'NFIN = ',i3,  &
    /10X,'IFIRST = ',i3,10X,'ILAST = ',i3,9X,'NFREE = ',i3)
1001 FORMAT (//7X,'IL',3X,'IH',14X,'JDIAG'//)
1002 FORMAT (28X,3(a1,2X))
1003 FORMAT (7X,i2,3X,i2,2X,a8,2X,3I3/)
1004 FORMAT (/5X,'TAB1'/)
1005 FORMAT (4(i3,1H),2X,i3,i5,5X))
1006 FORMAT (/2X,'SUMVAR = ',15(i3,l1))
1010 FORMAT (//10X,'J23',10X,'NBTR1 = ',i3//)
1012 FORMAT (18X,3(a1,2X))
1013 FORMAT (i9,i5,2X,3I3/)
1014 FORMAT (/3X,'J  L1 K1  L2 K2')
1015 FORMAT (4(i4,1H),i3,i3,i4,i3))
1020 FORMAT (/3X,a4,3X,3(20I3/))
1021 FORMAT (/3X,'DELTA = ',7(i5,i3))
1022 FORMAT (/3X,'KW(ARG. OF 6J)',6I3)
1030 FORMAT (//2X,'NC = ',i2,4X,'NPART = ',i2,4X,'IPARTL = ',i2,4X,  &
    'IPARTS = ',i2,4X,'ICROSS = ',i2,4X,/2X,'NPOINT = ',20I3)
1040 FORMAT (//2X,'NZERO = ',i2,5X,12(i4,1H),i3))
1050 FORMAT (///3X,'Printout after calling SUBROUTINE ',a7)

END SUBROUTINE printj
!***********************************************************************
!                                                                      *

SUBROUTINE rdiag (i,j,ik1,ik2,ichan,mat,jsum,j12)
!                                                                      *
!   Called by  GENSUM to establish the range of values of the summa-   *
!   tion variables.  This routine replaces an extended range do loop   *
!   in GENSUM, to conform with the FORTRAN 77 standard.                *
!                                                                      *
!                                           Last update: 02 Sep 1987   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: i
INTEGER, INTENT(IN OUT)                  :: j
INTEGER, INTENT(IN OUT)                  :: ik1
INTEGER, INTENT(IN OUT)                  :: ik2
INTEGER, INTENT(OUT)                     :: ichan
INTEGER, INTENT(IN OUT)                  :: mat(mtriad,mtriad)
INTEGER, INTENT(IN OUT)                  :: jsum(2,m6j)
INTEGER, INTENT(IN OUT)                  :: j12(4,mtriad)
IMPLICIT doubleprecision (a-h, o-z)

INTEGER, INTENT(IN OUT), PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m6j = 20

DIMENSION  jmnp(5),jmxp(5)


jmin1 = 0
jmax1 = 1000
k = mat(ik1,ik2)

DO  l1 = 1,k
  
  l3 = mat(j,j)
  jj1 = jsum(1,j)
  jnd = j12(l1,ik2,ik1)
  jmin = 1000
  jmax = 0
  jmnp(l1) = 0
  jmxp(l1) = 1000
  
  DO  l2 = 1,l3
    
    jmn = IABS(jnd-jj1)
    jmx = jnd+jj1
    jmin = MIN (jmn,jmin)
    jmax = MAX (jmx,jmax)
    jmnp(l1) = MAX (jmn,jmnp(l1))
    jmxp(l1) = MIN (jmx,jmxp(l1))
    jj1 = jj1+2
    
  END DO
  
  jmin1 = MAX (jmin1,jmin)
  jmax1 = MIN (jmax1,jmax)
  
END DO

IF (mat(i,i) == 0) THEN
  jsum(1,i) = jmin1
  jsum(2,i) = jmax1
  mat(i,i) = (jmax1-jmin1)/2+1
  ichan = ichan+1
  GO TO 3
END IF

IF (jsum(1,i) < jmin1) THEN
  jsum(1,i) = jmin1
  ichan = ichan+1
END IF

IF (jsum(2,i) > jmax1) THEN
  jsum(2,i) = jmax1
  ichan = ichan+1
END IF

3 k1 = 0

DO  l1 = 1,k
  IF ((jmnp(l1) <= jsum(1,i)) .AND. (jmxp(l1) >= jsum(2,i))) CYCLE
  k1 = k1+1
  j12(k1,ik2,ik1) = j12(l1,ik2,ik1)
END DO

IF (k1 /= k) THEN
  mat(ik1,ik2) = k1
  ichan = ichan+1
END IF

mat(ik2,ik1) = j12(1,ik2,ik1)

RETURN
END SUBROUTINE rdiag
!***********************************************************************
!                                                                      *

SUBROUTINE search (find)
!                                                                      *
!   This  routine locates circuits or loops of order  NC. NPOINT(NC)   *
!   are the  indices  of the points (triads) pertaining to the first   *
!   such loop found.  NPART  is the number of separate parts (groups   *
!   of contiguous points) on  the  axis of the flat graph. IPARTS is   *
!   the number of points in the smallest  part. IPARTL is the number   *
!   of points in  the  largest  part.  The  SUBROUTINE finds all the   *
!   possible loops  of  order  3  and 4. For NC .GE. 5, it looks for   *
!   only those who are partitionned in NPART .LE. 2. which can even-   *
!   tually reduce to a loop of  order  4  without breaking the basic   *
!   structure of the flat graph. ICROSS = -1, if lines cross.          *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


LOGICAL, INTENT(OUT)                     :: find
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME
INTEGER :: arr,tab1


doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc

DATA NAME/'SEARCH'/

!   Initialization

find = .false.
ncm1 = nc-1
ncm = nc-2
icross = 0

!   First treat two cases that do not involve do loops:

!   1. One isolated point, either the first or the last

npart = 1
ipartl = nc-1
iparts = 1

!   A. First

i1 = ifirst
k3 = 3
k2 = 2
200 ja = jdiag(i1,1)
jc = jdiag(i1,k3)

IF (ja == jc) THEN
  IF (nc > 1) THEN
    WRITE (*,300) i1,k3,ja,jc,nc
    STOP
  END IF
  npoint(1) = i1
  GO TO 900
END IF

i2 = tab1(ja,k2)
i3 = tab1(jc,k2)

IF (ABS(il(i3)-il(i2))-ncm < 0) THEN
  WRITE (*,301) i2,i3,ja,jc,k2,nc
  STOP
END IF

IF (ABS(il(i3)-il(i2))-ncm > 0) THEN
  
!   B. Last
  
  IF (i1 /= ifirst) GO TO 250
  i1 = ilast
  k3 = 2
  k2 = 1
  GO TO 200
END IF

ic = 1
npoint(ic) = i1
i20 = MIN (i2,i3)
i21 = il(i20)
i31 = i21+ncm1

DO  ii = i21,i31
  ic = ic+1
  npoint(ic) = ih(ii)
END DO

IF (nc <= 2) THEN
  IF (jdiag(ifirst,1) /= jdiag(ilast,1)) CALL phase (i1,jdiag,m4trd)
  GO TO 900
END IF

IF (i1 /= ilast) THEN
  it = i2
  jt = jdiag(ilast,2)
  k4 = 2
  i4 = ilast
ELSE
  it = i3
  jt = jdiag(ifirst,3)
  k4 = 3
  i4 = ifirst
END IF

IF (it == i20) CALL phase (i1,jdiag,m4trd)
IF ((jt == ja) .OR. (jt == jc)) CALL change (i4,k4)
GO TO 900

!   2. Two isolated points,first and last

250 IF (nc == 1) RETURN
IF (nc <= 3) GO TO 100
ipartl = nc-2
iparts = 1
i1 = ifirst
i2 = ilast
ja = jdiag(i1,1)
jb = jdiag(i1,3)

IF (tab1(ja,2) /= i2) THEN
  ja = jdiag(i1,3)
  jb = jdiag(i1,1)
  IF (tab1(ja,2) /= i2) GO TO 100
END IF

IF (ja == jdiag(i2,1)) THEN
  jc = jdiag(i2,2)
ELSE
  jc = jdiag(ilast,1)
END IF

i3 = tab1(jb,2)
i4 = tab1(jc,1)
idist = il(i4)-il(i3)

IF (ABS(idist)-(ncm-1) < 0) THEN
  WRITE (*,302) i3,i4,jb,jc,idist,nc
  STOP
END IF
IF (ABS(idist)-(ncm-1) == 0) THEN
  npoint(1) = ilast
  npoint(2) = ifirst
  icross = SIGN (1,idist)
  ic = 2
  i20 = MIN (i3,i4)
  i21 = il(i20)
  i31 = i21+ncm
  
  DO  ii = i21,i31
    ic = ic+1
    npoint(ic) = ih(ii)
  END DO
  
  IF (ja == jdiag(ifirst,1)) CALL change (ifirst,3)
  IF (ja == jdiag(ilast,1)) CALL change (ilast,2)
  GO TO 900
END IF

!   First general case: all points in one group

100 npart = 1
iparts = 0
ipartl = nc
k3 = 1

DO  in = 1,nbnode
  i = ih(in)
  108    ja = jdiag(i,k3)
  IF (i /= tab1(ja,2))THEN
    i2 = tab1(ja,2)
    
    IF (il(i2)-in-ncm1 < 0) THEN
      WRITE (*,303) in,i,i2,il(i2),ja,nc
      STOP
    END IF
    IF (il(i2)-in-ncm1 == 0) THEN
      i21 = il(i2)
      ic = 0
      
      DO  ii = in,i21
        ic = ic+1
        npoint(ic) = ih(ii)
      END DO
      
      IF (ja == jdiag(ifirst,3)) CALL change (ifirst,3)
      IF (ja == jdiag(ilast,2)) CALL change (ilast,2)
      GO TO 900
    END IF
  END IF
  
  IF (in == 1) THEN
    IF (k3 /= 3) THEN
      k3 = 3
      GO TO 108
    ELSE
      k3 = 1
    END IF
  END IF
  
END DO

!   Search did not find loop NC .LE. 3

IF (nc <= 3) RETURN

!   General case of loop partitionned in 2 groups. DO loop
!   on IPARTS

npart = 2
nc2 = nc/2
k3 = 1
k2 = 1

DO  ips = 2,nc2
  jps = ips-1
  nbn = nbnode-jps
  
  DO  i1 = 1,nbn
    i = ih(i1)
    i2 = ih(i1+jps)
    402       ja = jdiag(i,k3)
    jd = jdiag(i2,k2)
    
    IF (i == tab1(ja,1)) THEN
      ii2 = tab1(jd,2)
      ii1 = tab1(ja,2)
    ELSE
      ii1 = tab1(ja,1)
      ii2 = tab1(jd,1)
    END IF
    
    idist = il(ii1)-il(ii2)
    
    IF (ABS (idist)-(ncm-jps) < 0) THEN
      WRITE (*,304) jps,i1,i,i2,ja,jd,ii1,ii2,idist,nc
      STOP
    END IF
    IF (ABS (idist)-(ncm-jps) > 0) GO TO 420
    icross = SIGN (1,idist)
    ic = 0
    i21 = il(i2)
    
    DO  ii = i1,i21
      ic = ic+1
      npoint(ic) = ih(ii)
    END DO
    
    i20 = MIN (ii1,ii2)
    i30 = MAX (ii1,ii2)
    i21 = il(i20)
    i31 = il(i30)
    
    DO  ii = i21,i31
      ic = ic+1
      npoint(ic) = ih(ii)
    END DO
    
    iparts = ips
    ipartl = nc-ips
    IF ((jdiag(ifirst,3) == ja) .OR.  &
        (jdiag(ifirst,3) == jd)) CALL change (ifirst,3)
    IF ((jdiag(ilast,2) == ja) .OR.  &
        (jdiag(ilast,2) == jd)) CALL change (ilast,2)
    GO TO 900
    
    420       IF (i1 == 1) THEN
      IF (k3 == 3) THEN
        k3 = 1
        CYCLE
      ELSE
        k3 = 3
        GO TO 402
      END IF
    END IF
    
    IF (i2 == ilast) THEN
      IF (k2 /= 2) THEN
        k2 = 2
        GO TO 402
      END IF
    END IF
    
  END DO
END DO

!   SEARCH did not find circuit of order NC

RETURN

!   Loop found

900 find = .true.
CALL printj (NAME,10)

RETURN

!   Error printout

300 FORMAT (' Error in SEARCH. I1,K3,JA,JC,NC = ',5I5)
301 FORMAT (' Error in SEARCH. I2,I3,JA,JC,K2,NC = ',6I5)
302 FORMAT (' Error in SEARCH. I3,I4,JB,JC,IDIST,NC = ',6I5)
303 FORMAT (' Error in SEARCH. IN,I,I2,IL(I2),JA,NC = ',6I5)
304 FORMAT (' Error in SEARCH. JPS,I1,I,I2,JA,JD,II1,II2,IDIST,NC = ' ,10I5)

END SUBROUTINE search
!***********************************************************************
!                                                                      *

SUBROUTINE setdm
!                                                                      *
!   Sets dimensions of arrays.                                         *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************

IMPLICIT doubleprecision (a-h, o-z)
LOGICAL :: sumvar

PARAMETER ( mangm = 60,m3mngm = 3*mangm,mangmp = 2*(mangm/3),  &
    mtriad = 12,m2trd = 2*mtriad,m4trd = 4*mtriad, m6j = 20,msum = 10)

COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /dim/j6cc,j7cc,j8cc,j9cc,jwcc,jdelc

jwcc = jwc
jdelc = jdel
j6cc = j6c
j7cc = j7c
j8cc = j8c
j9cc = j9c

RETURN
END SUBROUTINE setdm
!***********************************************************************
!                                                                      *

SUBROUTINE settab (fail)
!                                                                      *
!   Builds up the unstructured graph. Sets the array J23, containing   *
!   the  two lists of original triads J2 and J3, and the correspond-   *
!   ing arrows  on the  angular  momenta lines. Also establishes the   *
!   numerical and phase factors  connecting  recoupling  coefficient   *
!   and graphs, according to Yutsis, Levinson, and Vanagas. For this   *
!   purpose determines the total J.                                    *
!                                                                      *
!                                           Last update: 16 Ocy 1992   *
!                                                                      *
!***********************************************************************


LOGICAL, INTENT(IN OUT)                  :: fail
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME
INTEGER :: arrow
LOGICAL :: tabs,free,sumvar

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),  &
    lcol(mangm,2),tabs(m2trd),nbtr  &
    /couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)  &
    /build/ial(m4trd),if1,if2,node  &
    /argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp

DATA NAME/'SETTAB'/

DO  i = 1,m2trd
  DO  j = 1,3
    j23(i,j) = 0
  END DO
END DO

ipr = n-1
nbtr = ipr+ipr

DO  i = 1,ipr
  DO  j = 1,2
    j23(i,j) = j2(i,j)
    arrow(i,j) = 1
  END DO
  tabs(i) = .false.
  j23(i,3) = j2(i,3)
  arrow(i,3) = -1
END DO

ipr1 = ipr+1

DO  i = ipr1,nbtr
  ii = i-ipr
  DO  j = 1,2
    j23(i,j) = j3(ii,j)
    arrow(i,j) = -1
  END DO
  tabs(i) = .false.
  j23(i,3) = j3(ii,3)
  arrow(i,3) = 1
END DO

DO  j = 1,nbtr
  j8(j) = j23(j,1)
END DO

j8c = nbtr+ipr
nb1 = nbtr+1

DO  j = nb1,j8c
  i = j-ipr
  j8(j) = j23(i,3)
END DO

j6c = nbtr

DO  j = 1,j6c
  j6(j) = j23(j,3)
END DO

DO  i = 1,m
  sumvar(i) = .false.
  ial(i) = 1
END DO

DO  i = 1,nbtr
  DO  j = 1,3
    ji = j23(i,j)
    k = ial(ji)
    line(ji,k) = i
    lcol(ji,k) = j
    ial(ji) = k+1
  END DO
END DO

it = 0

DO  i = 1,nbtr
  
  jt = j23(i,3)
  
  IF (ial(jt) == 3) THEN
    
    CALL otherj (i,jt,l,lc,k)
    IF (lc == 3) GO TO 19
    
  ELSE
    
    IF (it == 1) THEN
      CALL delta (jt1,jt,fail)
      IF (fail) GO TO 20
      k = line(jt,1)
      kc = lcol(jt,1)
      line(jt1,2) = k
      lcol(jt1,2) = kc
      line(jt,2) = line(jt1,1)
      lcol(jt,2) = lcol(jt1,1)
      j23(k,kc) = jt1
      ial(jt) = 1
      GO TO 19
    END IF
    
    jt1 = jt
    it = 1
    
  END IF
  
END DO

19 j9(j9c+1) = jt
j9c = j9c+2
j9(j9c) = jt

20 CALL printj (NAME,4)

RETURN
END SUBROUTINE settab
!***********************************************************************
!                                                                      *

SUBROUTINE sprate (m)
!                                                                      *
!   This  subroutine  prepares  the  information to be transfered to   *
!   GENSUM for numerical evaluation.                                   *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN)                      :: m
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME
CHARACTER (LEN=5) :: nme
LOGICAL :: sum6j,t6j,jt,js,sumvar,cut

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

DIMENSION jtem4(mtriad,m6j),jtem5(mtriad,m6j),jtem6(mtriad),  &
    nsum6j(m6j),j6sum(m6j)
DIMENSION sum6j(m6j),t6j(m6j),jt(mtriad),js(mtriad),  &
    inver(mangm),jnsum(mtriad),jinv(mtriad),n6jn(m6j),in6j(m6j), jsumt(m6j,6)

COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),jw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /cutdig/cut /dim/j6cc,j7cc,j8cc,j9cc,jwcc,jdelc  &
    /sumarg/j6p(mangmp),j7p(mangmp),j8p(mangmp),j9p(mangmp),  &
    jword(6,m6j),nlsum,nbj(msum),nb6j(msum),k6cp(msum),  &
    k7cp(msum),k8cp(msum),k9cp(msum),jsum6(mtriad),  &
    jsum4(mtriad,m6j),jsum5(mtriad,m6j),inv6j(m6j)

!   Test that array dimensions have not been exceeded.

IF (mp > mangm) THEN
  nmx = mangm
  npx = mp
  NAME = 'MANGM '
  nme  = 'MP   '
ELSE IF (jwc > m6j) THEN
  nmx = m6j
  npx = jwc
  NAME = 'M6J   '
  nme  = 'JWC  '
ELSE IF (j6c > m3mngm) THEN
  nmx = m3mngm
  npx = j6c
  NAME = 'M3MNGM'
  nme  = 'J6C  '
ELSE IF (j7c > m3mngm) THEN
  nmx = m3mngm
  npx = j7c
  NAME = 'M3MNGM'
  nme  = 'J7C  '
ELSE IF (j8c > m3mngm) THEN
  nmx = m3mngm
  npx = j8c
  NAME = 'M3MNGM'
  nme  = 'J8C  '
ELSE
  IF (j9c <= mangmp) GO TO 54
  nmx = mangmp
  npx = j9c
  NAME = 'MANGMP'
  nme  = 'J9C  '
END IF

60 WRITE (*,300) NAME,nme,npx,nmx
STOP

!   Determination of effective summation variables and their
!   relationships with 6j coefficients.

54 DO  i = 1,jwc
  inv6j(i) = 0
  sum6j(i) = .false.
END DO

nsum = 0
nlsum = 0
IF (mp == m) RETURN
m1 = m+1

DO  i = m1,mp
  IF (sumvar(i)) THEN
    nsum = nsum+1
    jsum6(nsum) = 0
    inver(i) = nsum
  END IF
END DO

IF (nsum == 0) RETURN

IF (nsum > mtriad) THEN
  nmx = mtriad
  npx = nsum
  NAME = 'MTRIAD'
  nme  = 'NSUM '
  GO TO 60
END IF

kt = 0

DO  i = 1,jwc
  DO  j = 1,6
    ik = jw(j,i)
    IF (.NOT. sumvar(ik)) CYCLE
    
    IF (.NOT. sum6j(i)) THEN
      sum6j(i) = .true.
      kt = kt+1
      j6sum(kt) = 0
      nsum6j(kt) = i
      inv6j(i) = kt
    END IF
    
    isk = inver(ik)
    i2 = jsum6(isk)+1
    jsum6(isk) = i2
    jsum4(isk,i2) = j
    jsum5(isk,i2) = kt
    i3 = j6sum(kt)+1
    j6sum(kt) = i3
    jsumt(kt,i3) = isk
  END DO
END DO

CALL var (j6,j6p,j6c,j6cp,j6cc,sumvar,mp,m,inver)
CALL var (j7,j7p,j7c,j7cp,j7cc,sumvar,mp,m,inver)
CALL var (j8,j8p,j8c,j8cp,j8cc,sumvar,mp,m,inver)
CALL var (j9,j9p,j9c,j9cp,j9cc,sumvar,mp,m,inver)

IF (.NOT. cut) THEN
  nlsum = 1
  nbj(1) = nsum
  nb6j(1) = kt
  k6cp(1) = j6cp
  k7cp(1) = j7cp
  k8cp(1) = j8cp
  k9cp(1) = j9cp
  
  DO  i = 1,kt
    i1 = nsum6j(i)
    DO  j = 1,6
      jword(j,i) = jw(j,i1)
    END DO
  END DO
  
  DO  i = 1,nsum
    isu = jsum6(i)
    DO  j = 1,isu
      i1 = jsum5(i,j)
      j1 = jsum4(i,j)
      jword(j1,i1) = mp+i
    END DO
  END DO
  
  RETURN
END IF

!   Separation of variables and sums in case a cut was detected.

k6c = 0
k7c = 0
k8c = 0
k9c = 0
nj = 0
n6j = 0

DO  i = 1,kt
  t6j(i) = .false.
END DO

DO  i = 1,nsum
  jt(i) = .false.
  js(i) = .false.
END DO

j = 1

10 nj = nj+1
jnsum(nj) = j
jinv(j) = nj
jt(j) = .true.
18 js(j) = .true.
js6 = jsum6(j)

DO  i = 1,js6
  i6j = jsum5(j,i)
  
  IF (.NOT. t6j(i6j)) THEN
    t6j(i6j) = .true.
    n6j = n6j+1
    n6jn(n6j) = nsum6j(i6j)
    in6j(i6j) = n6j
  END IF
  
  j6j = j6sum(i6j)
  
  DO  k = 1,j6j
    jk = jsumt(i6j,k)
    IF (.NOT. jt(jk)) THEN
      nj = nj+1
      jnsum(nj) = jk
      jinv(jk) = nj
      jt(jk) = .true.
    END IF
  END DO
  
END DO

DO  jj = 1,nsum
  j = jj
  IF ((.NOT. js(jj)) .AND. jt(jj)) GO TO 18
END DO

nlsum = nlsum+1

IF (nlsum > msum) THEN
  nmx = msum
  npx = nlsum
  NAME = 'MSUM  '
  nme = 'NLSUM'
  GO TO 60
END IF
nbj(nlsum) = nj
nb6j(nlsum) = n6j

IF (j6cp /= 0) CALL chvar (j6p,j6cp,k6c,jt,jinv,nsum)
k6cp(nlsum) = k6c
IF (j7cp /= 0) CALL chvar (j7p,j7cp,k7c,jt,jinv,nsum)
k7cp(nlsum) = k7c
IF (j8cp /= 0) CALL chvar (j8p,j8cp,k8c,jt,jinv,nsum)
k8cp(nlsum) = k8c
IF (j9cp /= 0) CALL chvar (j9p,j9cp,k9c,jt,jinv,nsum)
k9cp(nlsum) = k9c

IF (nj /= nsum) THEN
  DO  jj = 1,nsum
    j = jj
    IF (.NOT. jt(jj)) GO TO 10
  END DO
END IF

DO  i = 1,kt
  i1 = n6jn(i)
  DO  j = 1,6
    jword(j,i) = jw(j,i1)
  END DO
END DO

DO  i = 1,nsum
  ik = jnsum(i)
  i2 = jsum6(ik)
  jtem6(i) = i2
  DO  j = 1,i2
    jtem4(i,j) = jsum4(ik,j)
    k = jsum5(ik,j)
    jtem5(i,j) = in6j(k)
  END DO
END DO

DO  i = 1,nsum
  i2 = jtem6(i)
  jsum6(i) = i2
  DO  j = 1,i2
    i1 = jtem5(i,j)
    j1 = jtem4(i,j)
    jsum4(i,j) = j1
    jsum5(i,j) = i1
    jword(j1,i1) = i+mp
  END DO
END DO

RETURN

300 FORMAT (' Dimension error for ',a6  &
    /2X,a5,' = ',i5,' is out of allowed range',i4)

END SUBROUTINE sprate
!***********************************************************************
!                                                                      *

SUBROUTINE square
!                                                                      *
!   Reduces  a  circuit  of  order 4 in the two cases which are left   *
!   over by POLYGN, namely two disconnected groups of two points and   *
!   one group of two points plus  the  two  ends of the axis. In the   *
!   latter, the end of the axis is transferred  to the beginning. In   *
!   this  process,  one summation variable and two  6j  symbols  are   *
!   introduced.                                                        *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************

IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME,namsub
INTEGER :: arr,tab1
LOGICAL :: sumvar

PARAMETER ( mangm = 60,m3mngm = 3*mangm,mangmp = 2*(mangm/3),  &
    mtriad = 12,m2trd = 2*mtriad,m4trd = 4*mtriad, m6j = 20,msum = 10)

COMMON/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc  &
    /argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /nam/namsub

DATA NAME/'SQUARE'/

namsub = NAME
mp = mp+1
sumvar(mp) = .true.
it1 = npoint(1)
it2 = npoint(2)

IF (icross == 1) THEN
  it3 = npoint(3)
  it4 = npoint(4)
  k23 = 3
  k32 = 2
ELSE
  it3 = npoint(4)
  it4 = npoint(3)
  k23 = 2
  k32 = 3
END IF

l4 = jdiag(it2,1)

IF (arr(it2,1) <= 0) THEN
  CALL phase2 (l4)
  arr(it2,1) = 1
  arr(it3,1) = -1
END IF

l2 = jdiag(it1,1)
IF (arr(it1,1) > 0) CALL phase2 (l2)
jwc = jwc+1
kw(1,jwc) = l4
kw(2,jwc) = l2
kw(3,jwc) = jdiag(it2,2)
jj1 = jdiag(it1,3)
kw(4,jwc) = jj1
kw(5,jwc) = mp
kw(6,jwc) = jdiag(it1,2)
IF (arr(it1,2) < 0) CALL phase2 (jdiag(it1,2))
jwc = jwc+1
kw(1,jwc) = l4
kw(2,jwc) = l2
jj3 = jdiag(it3,k23)
jj2 = jdiag(it4,k32)
kw(3,jwc) = jj3
kw(4,jwc) = jj2
kw(5,jwc) = mp
kw(6,jwc) = jdiag(it3,k32)
IF (arr(it3,k32) < 0) CALL phase2 (jdiag(it3,k32))
j6(j6c+1) = mp
j6c = j6c+2
j6(j6c) = mp

IF (npart == 1) THEN
  itmin = it2
  itmax = it3
ELSE
  itmin = MIN (it2,it3)
  itmax = MAX (it2,it3)
END IF
itmn = MIN (it1,it4)
itmx = MAX (it1,it4)

tab1(mp,1) = itmin
tab1(mp,2) = itmax
jdiag(it2,1) = mp
jdiag(it3,1) = mp
jdiag(it2,3) = jj1
arr(it2,3) = arr(it1,3)
jdiag(it3,k32) = jj2
arr(it3,k32) = arr(it4,k32)

IF (icross == 1) THEN
  j7(j7c+1) = l2
  j7(j7c+2) = l4
  CALL phase2 (l4)
  j7c = j7c+3
  j7(j7c) = mp
ELSE
  CALL phase2 (jj2)
END IF

itll = il(itmn)
ithl = il(itmx)

DO  i = itll+1,ithl-1
  it = ih(i)
  ilp = i-1
  il(it) = ilp
  ih(ilp) = it
END DO
IF (ithl /= nbnode) THEN
  DO  i = ithl+1,nbnode
    it = ih(i)
    ilp = i-2
    il(it) = ilp
    ih(ilp) = it
  END DO
END IF

IF (npart /= 2) THEN
  tab1(jj1,1) = ih(1)
  tab1(jj1,2) = ih(nbnode-2)
END IF

RETURN
END SUBROUTINE square
!***********************************************************************
!                                                                      *

SUBROUTINE trdel (jj1,jj2,jj3,nbn,fail)
!                                                                      *
!   Test for triangular delta. If not satisfied FAIL = .TRUE. .        *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: jj1
INTEGER, INTENT(IN OUT)                  :: jj2
INTEGER, INTENT(IN OUT)                  :: jj3
INTEGER, INTENT(IN OUT)                  :: nbn
LOGICAL, INTENT(OUT)                     :: fail
IMPLICIT doubleprecision (a-h, o-z)
LOGICAL :: sumvar,cut,free

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /cutdig/cut /couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)

IF (sumvar(jj1) .OR. sumvar(jj2) .OR. sumvar(jj3)) RETURN
IF (nbn > 4) cut = .true.
IF ((.NOT. free(jj1)) .AND. (.NOT. free(jj2)) .AND. (.NOT.free(jj3))) THEN
  i1 = j1(jj1)
  i2 = j1(jj2)
  i3 = j1(jj3)
  IF ((i1 < (ABS (i2-i3)+1)) .OR. (i1 > (i2+i3-1))) fail = .true.
END IF

RETURN
END SUBROUTINE trdel
!***********************************************************************
!                                                                      *

SUBROUTINE triang (fail)
!                                                                      *
!   Reduces  a triangle having one apex at either end of the axis of   *
!   the flat  diagram.  This introduces one 6j symbol and some phase   *
!   factors.                                                           *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


LOGICAL, INTENT(IN OUT)                  :: fail
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME,namsub
INTEGER :: arr,tab1
LOGICAL :: sumvar

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/graph/jdiag(m4trd,3),arr(m4trd,3),tab1(mangm,2),il(m4trd),  &
    ih(m4trd),npoint(m2trd),nbnode,ifirst,ilast,iparts,  &
    ipartl,npart,icross,nfree,itfree(m6j),nfin,nc  &
    /argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp  &
    /nam/namsub

DATA NAME/'TRIANG'/

namsub = NAME
it1 = npoint(1)
it2 = npoint(2)
it3 = npoint(3)
jwc = jwc+1
kw(1,jwc) = jdiag(it3,2)
kw(2,jwc) = jdiag(it2,3)
kw(3,jwc) = jdiag(it3,1)
IF (arr(it3,1) > 0) CALL phase2 (kw(3,jwc))
kw(4,jwc) = jdiag(it2,1)
IF (arr(it2,1) < 0) CALL phase2 (kw(4,jwc))
k23 = 3
IF (it1 == ifirst) k23 = 2
kw(5,jwc) = jdiag(it1,k23)
kw(6,jwc) = jdiag(it3,3)
CALL trdel (kw(1,jwc),kw(2,jwc),kw(5,jwc),nbnode,fail)
IF (fail) GO TO 15
IF (arr(it3,3) > 0) CALL phase2 (kw(6,jwc))
jt1 = kw(5,jwc)
jdiag(it3,1) = jt1
jdiag(it3,3) = kw(2,jwc)
arr(it3,1) = arr(it1,k23)
arr(it3,3) = arr(it2,3)

IF (it1 /= ifirst) THEN
  tab1(jt1,1) = it3
  tab1(jt1,2) = ih(nbnode-1)
  k12 = 1
ELSE
  tab1(jt1,1) = ih(2)
  tab1(jt1,2) = it3
  k12 = 2
END IF

il3 = il(it3)

IF (it1 /= ilast) THEN
  il2 = il(it2)-1
  
  DO  i = 2,il2
    it = ih(i)
    ilp = i-1
    il(it) = ilp
    ih(ilp) = it
  END DO
END IF

DO  i = il3,nbnode
  it = ih(i)
  ilp = i-k12
  il(it) = ilp
  ih(ilp) = it
END DO

15 RETURN
END SUBROUTINE triang
!***********************************************************************
!                                                                      *

SUBROUTINE var (jn,jns,jnc,jnsc,jbc,sumvar,mp,m,inver)
!                                                                      *
!   Test  for  variable  character and put in JNS if yes, and JN now   *
!   contains 0.                                                        *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: jn(jnc)
INTEGER, INTENT(OUT)                     :: jns(mangmp)
INTEGER, INTENT(IN)                      :: jnc
INTEGER, INTENT(OUT)                     :: jnsc
INTEGER, INTENT(IN)                      :: jbc
LOGICAL, INTENT(IN OUT)                  :: sumvar(mp)
INTEGER, INTENT(IN OUT)                  :: mp
INTEGER, INTENT(IN)                      :: m
INTEGER, INTENT(IN)                      :: inver(mp)
IMPLICIT doubleprecision (a-h, o-z)


doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10



jnsc = 0
IF (jbc /= jnc) THEN
  jbbc = jbc+1
  
  DO  i = jbbc,jnc
    i1 = jn(i)
    IF (sumvar(i1)) THEN
      jnsc = jnsc+1
      IF (jnsc > mangmp) THEN
        WRITE (*,300) jnsc,mangmp
        STOP
      END IF
      j = inver(i1)
      jns(jnsc) = j
      jn(i) = m
    END IF
  END DO
END IF

RETURN

300 FORMAT (' Dimension error in VAR. JNSC = ',i5,' MANGMP = ',i5)

END SUBROUTINE var
!***********************************************************************
!                                                                      *

SUBROUTINE way (l,ka,kb,ich,nb)
!                                                                      *
!   Tests  one  step  forward  if  the way is free. First and second   *
!   arguments are interchanged or not according to ICH = -1, or +1.    *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: l
INTEGER, INTENT(IN OUT)                  :: ka
INTEGER, INTENT(IN OUT)                  :: kb
INTEGER, INTENT(OUT)                     :: ich
INTEGER, INTENT(OUT)                     :: nb
IMPLICIT doubleprecision (a-h, o-z)
INTEGER :: arrow
LOGICAL :: tabs

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10

COMMON/tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),  &
    lcol(mangm,2),tabs(m2trd),nbtr /build/ial(m4trd),if1,if2,node

k1 = j23(l,ka)
k2 = j23(l,kb)
nb = ial(k1)+ial(k2)-1
IF (nb < 0) THEN
  GO TO     3
ELSE IF (nb == 0) THEN
  GO TO     2
ELSE
  GO TO     8
END IF
2 nb1 = ial(k1)-ial(k2)
IF (nb1 < 0) THEN
  GO TO     9
ELSE
  GO TO     8
END IF
3 CALL otherj (l,k1,l1,lc1,la)
CALL otherj (l,k2,l2,lc2,lb)
CALL neibor (lc1,i1,i2)
CALL neibor (lc2,i3,i4)
ji1 = j23(l1,i1)
ji2 = j23(l1,i2)
ji3 = j23(l2,i3)
ji4 = j23(l2,i4)
ia = ial(ji1)+ial(ji2)
ib = ial(ji3)+ial(ji4)
nbp = ib+ia+1
nbm = ib-ia
SELECT CASE ( nbp )
  CASE (    1)
    GO TO 8
  CASE (    2)
    GO TO 4
  CASE (    3)
    GO TO 5
  CASE (    4)
    GO TO 4
  CASE (    5)
    GO TO 6
END SELECT
4 IF (nbm < 0) THEN
  GO TO     9
ELSE
  GO TO     8
END IF
5 IF (nbm < 0) THEN
  GO TO     9
ELSE IF (nbm == 0) THEN
  GO TO     6
ELSE
  GO TO     8
END IF
6 IF ((ji3 == if1) .OR. (ji3 == if2) .OR.  &
    (ji4 == if1) .OR. (ji4 == if2)) GO TO 9
8 ich = 1
GO TO 10
9 ich = -1
10 RETURN

END SUBROUTINE way
!***********************************************************************
!                                                                      *

SUBROUTINE zero (j,jz,fail)
!                                                                      *
!   Suppresses  one  line  and  two  nodes of the unstructured graph   *
!   introduces  zeros in the triads  J23. As a consequence the other   *
!   two arguments of the triad are put equal. If there was already a   *
!   zero in the triad which is changed, it is a special case.          *
!                                                                      *
!                                           Last update: 16 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(OUT)                     :: j
INTEGER, INTENT(OUT)                     :: jz
LOGICAL, INTENT(IN OUT)                  :: fail
IMPLICIT doubleprecision (a-h, o-z)
CHARACTER (LEN=6) :: NAME
INTEGER :: arrow
LOGICAL :: tabs,free,sumvar,cut,nocut

doubleprecision, PARAMETER ::mangm = 60
INTEGER, PARAMETER :: m3mngm = 3*mangm
INTEGER, PARAMETER :: mangmp = 2*(mangm/3)
INTEGER, PARAMETER :: mtriad = 12
INTEGER, PARAMETER :: m2trd = 2*mtriad
INTEGER, PARAMETER :: m4trd = 4*mtriad
INTEGER, PARAMETER :: m6j = 20
INTEGER, PARAMETER :: msum = 10
INTEGER, PARAMETER :: mzero = 20

COMMON/zer/nzero,jzero(mzero) /cutdig/cut  &
    /couple/m,n,j1(mangm),j2(mtriad,3),j3(mtriad,3),free(mangm)  &
    /KEEP/jkp(2,3),jarr(2,3),it2,it3,it5 /build/ial(m4trd),if1,if2,node  &
    /tree/j23(m2trd,3),arrow(m2trd,3),line(mangm,2),  &
    lcol(mangm,2),tabs(m2trd),nbtr  &
    /argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm), mp

DATA NAME/'ZERO  '/

nocut = .false.
nzero = 0

IF (j >= 1) THEN
  CALL otherj (0,jz,lin,lc,k1)
  i = nzero
  GO TO 8
END IF

DO  i = 1,m
  IF ((j1(i) /= 1) .OR. free(i) .OR. (ial(i) <= 1)) CYCLE
  nzero = nzero+1
  IF (nzero > mzero) THEN
    WRITE (*,300) nzero,mzero
    STOP
  END IF
  jzero(nzero) = i
END DO

nocut = .true.
m = m+1
j1(m) = 1
sumvar(m) = .false.
free(m) = .false.
IF (nzero == 0) GO TO 7
CALL printj (NAME,1)
i = 0
1 i = i+1
jz = jzero(i)
j = 0
13 j = j+1
lin = line(jz,j)
IF (tabs(lin)) GO TO 2
lc = lcol(jz,j)
8 CALL neibor (lc,l1,l2)
jj1 = j23(lin,l1)
jj2 = j23(lin,l2)

IF (jj1 == jj2) THEN
  j6c = j6c+1
  j6(j6c) = jj1
  lo1=lin
  lo2=lin
  lco1=l1
  lco2=l2
  GO TO 10
END IF

CALL delta (jj1,jj2,fail)
IF (fail) GO TO 7

IF ((j1(jj1) /= 1) .AND. (j1(jj2) /= 1)) GO TO 15
IF (j1(jj1) < j1(jj2)) GO TO 15
IF (j1(jj1) > j1(jj2)) GO TO 19

IF (nzero /= 0) THEN
  DO  jjx = i,nzero
    jjz = jzero(jjx)
    IF (jj1 == jjz) GO TO 15
    IF (jj2 == jjz) GO TO 19
  END DO
END IF

GO TO 15

19 jjz = jj2
jj2 = jj1
jj1 = jjz

15 CALL otherj (lin,jj1,lo1,lco1,k1)
CALL otherj (lin,jj2,lo2,lco2,k2)
j9c = j9c+1
j9(j9c) = jj1
j23(lo2,lco2) = jj1
line(jj1,k1) = lo2
lcol(jj1,k1) = lco2

10 IF     (arrow(lin,l1) < arrow(lin,l2)) THEN
  CALL phase2 (jj1)
ELSE IF (arrow(lin,l1) == arrow(lin,l2)) THEN
  arrow(lo1,lco1) = 1
  arrow(lo2,lco2) = -1
END IF

tabs(lin) = .true.
nbtr = nbtr-1
IF (nbtr == 0) GO TO 7
IF (lo1 == lo2) THEN
  l = 6-lco1-lco2
  jt = j23(lo1,l)
  IF ((j1(jt) == 1) .AND. (.NOT.free(jt))) GO TO 2
  CALL delta (jt,m,fail)
  IF (fail) GO TO 7
  nzero = nzero+1
  jzero(nzero) = jt
END IF
2 IF (j == 1) GO TO 13

IF (nbtr /= 0) THEN
  IF (i < nzero) GO TO 1
END IF

7 CALL printj (NAME,4)
IF (nocut) cut = .false.

RETURN

300 FORMAT (' Dimension error in ZERO. NZERO = ',i5,' MZERO = ',i5)

END SUBROUTINE zero
