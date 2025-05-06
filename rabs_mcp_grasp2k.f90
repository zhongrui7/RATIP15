!***********************************************************************

! Code converted using TO_F90 by Alan Miller
! Date: 2020-02-25  Time: 21:27:54

! This file contains several procedures from the GRASP-92 package
! to enable the calculation of angular coefficients. Since parts of this
! code date back to the late sixties, this code is (almost) impossible
! to modify; this has hampared the implementation of a more efficient
! computation of angular coefficients in the past.

! No attempt has been made to 'improve' the readability of this code;
! a few modifications need to be done as are clearly indicated beloww.
! This file just list all required procedures (as a list of external
! Fortran-77 procedures without any explicit interface) in alphabetic
! order. -- In a long term, we intent to replace this file and to
! implement a more efficient scheme using the features of Fortran-90.

!***********************************************************************
!                                                                      *

SUBROUTINE breid (ja,jb,ja1,ipca,jb1)
!                                                                      *
!   Computes closed shell contributions - aaaa and exchange only.      *
!                                                                      *
!   Call(s) to: [LIB92]: CLRX, CXK, ITRIG, TALK, SNRC.                 *
!                                                                      *
!                                           LAST UPDATE: 09 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: ja
INTEGER, INTENT(IN OUT)                  :: jb
INTEGER, INTENT(IN OUT)                  :: ja1
INTEGER, INTENT(IN OUT)                  :: ipca
INTEGER, INTENT(IN OUT)                  :: jb1
IMPLICIT doubleprecision (a-h, o-z)

DIMENSION cone(7,20),js(4),kaps(4),ks(4),s(12)

COMMON/bcore/icore(149) /cons/zero,half,tenth,one,two,three,ten  &
    /debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6 /m1/nq1(149),nq2(149)  &
    /m3/jlist(149),klist(149),npeel,ncore /orb4/np(149),nak(149)

INTEGER, PARAMETER :: numax = 20

!   1.0  Initialization

IF (ipca == 2) THEN
  ia1 = klist(ja1)
ELSE
  ia1 = jlist(ja1)
END IF
ib1 = klist(jb1)

isg = 1
IF (ja == jb) THEN
  IF ((icore(ia1) /= 0) .AND. (icore(ib1) /= 0)) THEN
    IF (ja > 1) RETURN
    isg = -1
  END IF
END IF

js(1) = ia1
js(2) = ib1
js(3) = ia1
js(4) = ib1
nqs1 = nq1(ia1)
nqs2 = nq2(ib1)
DO  i = 1,4
  kaps(i) = 2*nak(js(i))
  ks(i) = IABS (kaps(i))
END DO
const = nqs1*nqs2
IF (ibug2 /= 0) WRITE (99,300) ia1,ib1

!   2.0  Set range of tensor indices

CALL snrc (js,kaps,ks,nd1,nd2,ne1,ne2,ibrd,ibre)
IF (ibug2 /= 0) WRITE (99,301) nd1,nd2,ne1,ne2,ibrd,ibre
IF (ia1 /= ib1) GO TO 3

!   3.0 Calculate aaaa interaction

DO  n = 1,nd2
  nu = nd1+2*(n-1)
  k = nu
  IF (MOD (k,2) /= 1) RETURN
  kap1 = kaps(1)/2
  gam = clrx (kap1,nu,kap1)
  dksks = ks(1)*ks(1)
  dnunu1 = nu*(nu+1)
  coef = const*two*dksks*gam*gam/dnunu1
  IF (ibug2 /= 0) WRITE (99,302) nu,gam,coef
  itype = isg*4
  CALL talk (ja,jb,nu,ia1,ia1,ia1,ia1,itype,coef)
END DO
RETURN

!   Calculate exchange interactions

3 CONTINUE
IF (ibre < 0) RETURN
IF (ne2 > numax) THEN
  WRITE (*,304)
  STOP
END IF

DO  n = 1,ne2
  DO  mu = 1,7
    cone(mu,n) = zero
  END DO
END DO

proc = -const/DBLE (ks(1)*ks(2))

!   Negative sign arises from Pauli phase factor

DO  n = 1,ne2
  nu = ne1+2*(n-1)
  k = nu
  ip = (ks(1)-ks(2))/2+k
  ipp = ip+1
  IF (nu == 0) GO TO 8
  kk = k+k+1
  IF (itrig(ks(1),ks(2),kk) == 0) GO TO 6
  prod = proc
  IF (MOD (ip,2) /= 0) prod = -prod
  CALL cxk (s,js,kaps,nu,k,ibre,2)
  IF (ibug2 /= 0) WRITE (99,303) prod,(s(mu),mu = 1,3)
  DO  mu = 1,3
    cone (mu,n) = cone(mu,n)+prod*s(mu)
  END DO
  
  6    k = nu-1
  kk = k+k+1
  IF (itrig(ks(1),ks(2),kk) == 0) GO TO 8
  prod = proc
  IF (MOD (ipp,2) /= 0) prod = -prod
  CALL cxk (s,js,kaps,nu,k,ibre,2)
  IF (ibug2 /= 0) WRITE (99,303) prod,(s(mu),mu = 1,3)
  DO  mu = 1,3
    cone(mu,n) = cone(mu,n)+prod*s(mu)
  END DO
  
  8    IF (n == ne2) GO TO 11
  k = nu+1
  kk = k+k+1
  prod = proc
  IF (MOD (ipp,2) /= 0) prod = -prod
  CALL cxk (s,js,kaps,nu,k,ibre,2)
  IF (ibug2 /= 0) WRITE (99,303) prod,(s(mu),mu = 1,7)
  DO  mu = 1,7
    cone (mu,n) = cone(mu,n)+prod*s(mu)
  END DO
END DO

!   4.0  Output results

11 CONTINUE

DO  n = 1,ne2
  nu = ne1+2*(n-1)
  itype = isg*5
  CALL talk (ja,jb,nu,ib1,ia1,ib1,ia1,itype,cone(1,n))
  CALL talk (ja,jb,nu,ia1,ib1,ib1,ia1,itype,cone(2,n))
  CALL talk (ja,jb,nu,ia1,ib1,ia1,ib1,itype,cone(3,n))
  IF (n == ne2) CYCLE
  nup1 = nu+1
  itype = isg*6
  CALL talk (ja,jb,nup1,ia1,ib1,ia1,ib1,itype,cone(4,n))
  CALL talk (ja,jb,nup1,ib1,ia1,ib1,ia1,itype,cone(5,n))
  CALL talk (ja,jb,nup1,ia1,ib1,ib1,ia1,itype,cone(6,n))
  CALL talk (ja,jb,nup1,ib1,ia1,ia1,ib1,itype,cone(7,n))
END DO
RETURN

300 FORMAT ('BREID: orbitals ',2I3)
301 FORMAT (2X,'ND1 ND2 NE1 NE2 IBRD IBRE ',6I5)
302 FORMAT (2X,'aaaa contribution: NU,GAM,COEF',i5,2(3X,1PD15.8))
303 FORMAT (2X,'PROD = ',1PD15.8 /' S',7D15.8)
304 FORMAT ('BREID: Dimension error for NUMAX.')

END SUBROUTINE breid
!***********************************************************************
!                                                                      *

SUBROUTINE breit (ja,jb,ja1,jb1,ja2,jb2)
!                                                                      *
!   Computes  the  coefficients  appearing in EQS. 5, 8, 9 AND 10 OF   *
!   I P Grant and B J McKenzie,  J Phys B 13 (1980) 2671--2681.  The   *
!   coefficients for each choice of  orbitals JA1, JB1, JA2, and JB2   *
!   depend on two further parameters  NU and K;  there are IMU inte-   *
!   grals for each such choice, where:                                 *
!                                                                      *
!                  IMU = 4          TYPE = 1                           *
!                        8                 2                           *
!                        1                 3                           *
!                        1                 4                           *
!                        3                 5                           *
!                        4                 6                           *
!                                                                      *
!   See the paper cited above for details.                             *
!                                                                      *
!   Call(s) to: [LIB92]: GENSUM, ITRIG, KNJ, LTAB, MODJ23, MUMDAD,     *
!                        NJGRAF, OCON, SETJ                            *
!               [RCI92]: CXK, SNRC, TALK.                              *
!                                                                      *
!                                           Last update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: ja
INTEGER, INTENT(IN OUT)                  :: jb
INTEGER, INTENT(IN)                      :: ja1
INTEGER, INTENT(IN)                      :: jb1
INTEGER, INTENT(IN)                      :: ja2
INTEGER, INTENT(IN)                      :: jb2
IMPLICIT doubleprecision (a-h, o-z)
INTEGER :: pntriq

doubleprecision, PARAMETER ::mangm=60
INTEGER, PARAMETER :: m3mngm=3*mangm
INTEGER, PARAMETER :: mangmp=2*(mangm/3)
INTEGER, PARAMETER :: mtriad=12
INTEGER, PARAMETER :: m2trd=2*mtriad
INTEGER, PARAMETER :: m4trd=4*mtriad
INTEGER, PARAMETER :: m6j=20
INTEGER, PARAMETER :: msum=10

LOGICAL :: free,dsumvr,esumvr,faild,faile

DIMENSION cond(12,20),cone(12,20),s(12)
DIMENSION is(4),kaps(4),ks(4),nqs(4),ils(4),lls(4),it1(4),irows(4)
DIMENSION js(4)

DIMENSION jd6(m3mngm),jd7(m3mngm),jd8(m3mngm),  &
    jd9(mangmp),kdw(6,m6j),lddel(m6j,2),dsumvr(mangm)
DIMENSION jd6p(mangmp),jd7p(mangmp),jd8p(mangmp),jd9p(mangmp),  &
    jdword(6,m6j), ndbj(msum),ndb6j(msum),kd6cp(msum),kd7cp(msum),kd8cp(msum),  &
    kd9cp(msum),jdsum6(mtriad),jdsum4(mtriad,m6j),jdsum5(mtriad,m6j), invd6j(m6j)
DIMENSION je6(m3mngm),je7(m3mngm),je8(m3mngm),  &
    je9(mangmp),kew(6,m6j),ledel(m6j,2),esumvr(mangm)
DIMENSION je6p(mangmp),je7p(mangmp),je8p(mangmp),je9p(mangmp),  &
    jeword(6,m6j), nebj(msum),neb6j(msum),ke6cp(msum),ke7cp(msum),ke8cp(msum),  &
    ke9cp(msum),jesum6(mtriad),jesum4(mtriad,m6j),jesum5(mtriad,m6j), inve6j(m6j)

COMMON/couple/mja,nja,j1(mangm),j2(mtriad,3),j3(mtriad,3), free(mangm)  &
    /debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6  &
    /l1/jbq1(3,149),jbq2(3,149),jtq1(3),jtq2(3) /l2/j2s(mtriad,3),j3s(mtriad,3)  &
    /m1/nq1(149),nq2(149) /m2/jjq1(3,149),jjq2(3,149)  &
    /m3/jlist(149),klist(149),npeel,ncore /orb2/ncf,nw,pntriq  &
    /orb4/np(149),nak(149) /terms/nrows,itab(31),jtab(32),ntab(327)

doubleprecision, PARAMETER :: eps = 1.0D-10
INTEGER, PARAMETER :: numax = 20

!   1.0  Initialize pointers and flags and set any
!        tables required.

!        In this segment, the array IS points to the
!        full list of orbitals, the array JS to the
!        array JLIST of peel orbital pointerS.

!   1.1  Initialization

js(1) = ja1
js(2) = jb1
js(3) = ja2
js(4) = jb2
DO  i = 1,4
  is(i) = jlist(js(i))
  kaps(i) = 2*nak(is(i))
  ks(i) = ABS (kaps(i))
END DO
ia1 = is(1)
ib1 = is(2)
ia2 = is(3)
ib2 = is(4)
nqs(1) = nq1(ia1)
nqs(2) = nq1(ib1)
nqs(3) = nq2(ia2)
nqs(4) = nq2(ib2)

kj23 = 0
isnj = 0

faild = .false.
faile = .false.
nbrj = 3*npeel + 7
DO  i = 1,(nbrj-1)
  free(i) = .false.
END DO
free(nbrj) = .true.

!   2.0  Set quantum numbers of spectator shells.

DO  j = 1,nw
  DO  k = 1,3
    jbq1(k,j) = 0
    jbq2(k,j) = 0
  END DO
END DO

DO  jj = 1,npeel
  j = jlist(jj)
  IF ((j /= ia1) .AND. (j /= ib1)) THEN
    DO  k = 1,3
      jbq1(k,j) = jjq1(k,j)
    END DO
  END IF
  IF ((j /= ia2) .AND. (j /= ib2)) THEN
    DO  k = 1,3
      jbq2(k,j) = jjq2(k,j)
    END DO
  END IF
  
!   2.1  Examine spectator shells for orthogonality
  
  IF ((j /= ia1) .AND. (j /= ib1) .AND. (j /= ia2) .AND. (j /= ib2)) THEN
    DO  k = 1,3
      IF (jbq1(k,j) /= jbq2(k,j) ) GO TO  98
    END DO
  END IF
END DO

!   3.0  Start main calculation
!        Begin with common factors

const = ocon (ia1,ib1,ia2,ib2)
IF (ibug2 /= 0) WRITE (99,307) const

!   3.1  Set range of tensor index NU

IF (ibug2 /= 0) WRITE (99,302) ia1,ib1,ia2,ib2
CALL snrc (is,kaps,ks,nd1,nd2,ne1,ne2,ibrd,ibre)
IF (ibug2 /= 0) WRITE (99,303) nd1,nd2,ne1,ne2,ibrd,ibre
IF ((ibrd < 0) .AND. (ibre < 0)) RETURN
IF ((nd2 > numax) .OR. (ne2 > numax)) THEN
  kk = MAX (nd2,ne2)
  WRITE (*,301) kk
  STOP
END IF
IF (ibrd >= 0) THEN
  DO  n = 1,nd2
    DO  mu = 1,12
      cond(mu,n) = 0.0D 00
    END DO
  END DO
END IF
IF (ibre >= 0) THEN
  DO  n = 1,ne2
    DO  mu = 1,12
      cone(mu,n) = 0.0D 00
    END DO
  END DO
END IF

!   3.2  Set parameters of summation over parent
!        (barred) terms in Eq. 2 (loc cit). The array
!        IROWS is formed to point to the list of
!        allowed parents of active shells in the
!        array NTAB

CALL ltab (is,nqs,ks,irows)

DO  i = 1,4
  ii = irows(i)
  lls(i) = itab(ii)
  ils(i) = jtab(ii)
END DO

!   4.0  Sum over all parent terms permitted by
!        angular momentum and seniority selection rules

lls1 = lls(1)
IF (lls1 /= 1) free(ja1) = .true.
lls2 = lls(2)
lls3 = lls(3)
lls4 = lls(4)

ls2 = ils(2)
DO  lb1 = 1,lls2
  ls2 = ls2+3
  it1(2) = ntab(ls2)
  it12 = it1(2)
  it2 = ks(2)
  it3 = jjq1(3,ib1)
  IF (itrig (it12,it2,it3) == 0) CYCLE
  IF (ABS (ntab(ls2-2)-jjq1(1,ib1) ) /= 1) CYCLE
  
  ls1 = ils(1)
  DO  la1 = 1,lls1
    ls1 = ls1+3
    it1(1) = ntab(ls1)
    it11 = it1(1)
    it2 = ks(1)
    IF (ia1 == ib1) THEN
      
!   Treat IA1 .EQ. IB1 as a special case
      
      it3 = it1(2)
      IF (itrig (it11,it2,it3) == 0) CYCLE
      IF (ABS (ntab(ls1-2)-ntab(ls2-2)) /= 1) CYCLE
      IF (lls2 /= 1) free(nbrj-8) = .true.
    ELSE
      it3 = jjq1(3,ia1)
      IF (itrig (it11,it2,it3) == 0) CYCLE
      IF (ABS (ntab(ls1-2)-jjq1(1,ia1)) /= 1) CYCLE
      IF (lls2 /= 1) free(jb1) = .true.
    END IF
    
    ls4 = ils(4)
    DO  lb2 = 1,lls4
      ls4 = ls4+3
      it1(4) = ntab(ls4)
      it14 = it1(4)
      it2 = ks(4)
      it3 = jjq2(3,ib2)
      IF (itrig(it14,it2,it3) == 0) CYCLE
      IF (ABS (ntab(ls4-2)-jjq2(1,ib2)) /= 1) CYCLE
      
      ls3 = ils(3)
      loop26:  DO  la2 = 1,lls3
        ls3 = ls3+3
        it1(3) = ntab(ls3)
        it13 = it1(3)
        it2 = ks(3)
        IF (ia2 == ib2) THEN
          
!   TREAT IA2 .EQ. IB2 as a special case
          
          it3 = it1(4)
          IF (lls4 /= 1) free(nbrj-6) = .true.
          IF (itrig (it13,it2,it3) == 0) CYCLE loop26
          IF (ABS (ntab(ls3-2)-ntab(ls4-2)) /= 1) CYCLE loop26
        ELSE
          it3 = jjq2(3,ia2)
          IF (itrig (it13,it2,it3) == 0) CYCLE loop26
          IF (ABS (ntab(ls3-2)-jjq2(1,ia2)) /= 1) CYCLE loop26
        END IF
        
!   At this point the current parent has been completely defined,
!   and its quantum numbers can now be set.  The JTQ arrays must
!   be set if IA1 .EQ. IB1 or IA2 .EQ. IB2. The matrix element should be
!   diagonal in barred quantum numbers.
        
        DO  k = 1,3
          jbq1(k,ia1) = ntab(ls1+k-3)
          jbq2(k,ia2) = ntab(ls3+k-3)
          jtq1(k) = 0
          IF (ib1 == ia1) THEN
            jtq1(k) = ntab(ls2+k-3)
          ELSE
            jbq1(k,ib1) = ntab(ls2+k-3)
          END IF
          jtq2(k) = 0
          IF (ib2 == ia2) THEN
            jtq2(k) = ntab(ls4+k-3)
          ELSE
            jbq2(k,ib2) = ntab(ls4+k-3)
          END IF
          DO  kk = 1,4
            IF (jbq1(k,is(kk)) /= jbq2(k,is(kk))) CYCLE loop26
          END DO
        END DO
        
!   4.1 Evaluate product of 4 CFPs
        
        CALL mumdad (is,kaps,prod)
        IF (ABS (prod) < eps) CYCLE loop26
        
!    4.2  Set arrays for defining the recoupling
!         coefficient
        
        CALL setj (is,js,ks,npeel,kj23)
        
        IF (isnj == 0) THEN
          
!******************** N J G R A F   V E R S I O N **********************
          
!     Set up the arrays and variables for the direct case.
          
          IF (ibrd >= 0) THEN
            CALL njgraf (recup,-1,faild)
            isnj = 1
            IF (.NOT. faild) THEN
              CALL knj(jd6c,jd7c,jd8c,jd9c,jdwc,jd6,jd7,jd8,jd9,  &
                  kdw,jddel,lddel,dsumvr,mdp,  &
                  jd6p,jd7p,jd8p,jd9p,jdword,ndlsum,ndbj,  &
                  ndb6j,kd6cp,kd7cp,kd8cp,kd9cp, jdsum4,jdsum5,jdsum6,invd6j)
            END IF
          END IF
          
!   Set up the arrays and variables for the exchange case.
          
          IF (ibre >= 0) THEN
            CALL modj23
            CALL njgraf (recup,-1,faile)
            isnj = 2
            IF (.NOT. faile) THEN
              CALL knj(je6c,je7c,je8c,je9c,jewc,je6,je7,je8,je9,  &
                  kew,jedel,ledel,esumvr,mep,  &
                  je6p,je7p,je8p,je9p,jeword,nelsum,nebj,  &
                  neb6j,ke6cp,ke7cp,ke8cp,ke9cp, jesum4,jesum5,jesum6,inve6j)
            END IF
          END IF
          
          IF (faild .AND. faile) GO TO 30
          
        END IF
        
!   4.3.1 Summation for direct terms
        
        IF ((ibrd >= 0). AND. (.NOT. faild)) THEN
          
          imud = 4
          IF (ibrd > 1) imud = 1
          ncode = 0
          DO  n = 1,nd2
            nu = nd1+2*(n-1)
            nud = nu+nu+1
            
            IF (nu /= 0) THEN
              
              IF ((itrig (ks(1),ks(3),nud) /= 0) .AND.  &
                    (itrig (ks(2),ks(4),nud) /= 0)) THEN
                
                k = nu
                j1(mja) = nud
                CALL gensum (jd6c,jd7c,jd8c,jd9c,jdwc,jd6,jd7,jd8,  &
                    jd9,kdw,jddel,lddel,dsumvr,mdp,  &
                    jd6p,jd7p,jd8p,jd9p,jdword,ndlsum,ndbj,ndb6j,  &
                    kd6cp,kd7cp,kd8cp,kd9cp,jdsum4,jdsum5,jdsum6, invd6j,x)
                IF (ibug2 /= 0) WRITE (99,304) nu,k,x
                
                IF (ABS (x) >= eps) THEN
                  x = x*prod
                  CALL cxk(s,is,kaps,nu,k,ibrd,1)
                  IF (ibug2 /= 0) WRITE (99,305) (s(iii),iii = 1,imud)
                  DO  mu = 1,imud
                    cond(mu,n) = cond(mu,n)+x*s(mu)
                  END DO
                END IF
                
              END IF
              
!   K = NU-1
              
              IF (ibrd > 1) CYCLE
              
              k = nu-1
              
              IF (ncode == n) THEN
                x=xcode
              ELSE
                itkmo = nud-2
                IF (itrig (ks(1),ks(3),itkmo) == 0) GO TO 18
                IF (itrig (ks(2),ks(4),itkmo) == 0) GO TO 18
                j1(mja) = itkmo
                CALL gensum(jd6c,jd7c,jd8c,jd9c,jdwc,jd6,jd7,jd8,  &
                    jd9,kdw,jddel,lddel,dsumvr,mdp,  &
                    jd6p,jd7p,jd8p,jd9p,jdword,ndlsum,ndbj,ndb6j,  &
                    kd6cp,kd7cp,kd8cp,kd9cp,jdsum4,jdsum5,jdsum6, invd6j,x)
              END IF
              
              IF (ibug2 /= 0) WRITE (99,304) nu,k,x
              
              IF (ABS (x) >= eps) THEN
                x = x*prod
                CALL cxk (s,is,kaps,nu,k,ibrd,1)
                IF (ibug2 /= 0) WRITE (99,305) (s(iii),iii = 1,4)
                DO  mu = 1,4
                  cond (mu,n) = cond(mu,n)+x*s(mu)
                END DO
              END IF
              
            END IF
            
!   K = NU+1
            
            18           IF ((ibrd > 1) .OR. (n == nd2)) CYCLE
            
            ncode = n+1
            xcode = 0.0D 00
            itkmo = nud+2
            
            IF ((itrig(ks(1),ks(3),itkmo) /= 0) .AND.  &
                  (itrig(ks(2),ks(4),itkmo) /= 0)) THEN
              k = nu+1
              j1(mja) = itkmo
              CALL gensum (jd6c,jd7c,jd8c,jd9c,jdwc,jd6,jd7,jd8,  &
                  jd9,kdw,jddel,lddel,dsumvr,mdp,  &
                  jd6p,jd7p,jd8p,jd9p,jdword,ndlsum,ndbj,ndb6j,  &
                  kd6cp,kd7cp,kd8cp,kd9cp,jdsum4,jdsum5,jdsum6, invd6j,x)
              xcode = x
              IF (ibug2 /= 0) WRITE (99,304) nu,k,x
              
              IF (ABS (x) >= eps) THEN
                x = x*prod
                CALL cxk (s,is,kaps,nu,k,ibrd,1)
                IF (ibug2 /= 0) WRITE (99,305) (s(iii),iii = 1,12)
                DO  mu = 1,12
                  cond(mu,n) = cond(mu,n)+x*s(mu)
                END DO
              END IF
              
            END IF
            
          END DO
          
        END IF
        
!   4.3.2 Summation for exchange terms
        
        IF ((ibre >= 0) .AND. (.NOT. faile)) THEN
          
          ncode = 0
          
          DO  n = 1,ne2
            imue = 4
            IF (ibre == 2) imue = 1
            IF (ibre == 4) imue = 3
            nu = ne1+2*(n-1)
            nud = nu+nu+1
            
            IF (nu /= 0) THEN
              
              IF ((itrig(ks(1),ks(4),nud) /= 0) .AND.  &
                    (itrig(ks(2),ks(3),nud) /= 0)) THEN
                k = nu
                j1(mja) = nud
                CALL gensum(je6c,je7c,je8c,je9c,jewc,je6,je7,je8,  &
                    je9,kew,jedel,ledel,esumvr,mep,  &
                    je6p,je7p,je8p,je9p,jeword,nelsum,nebj,neb6j,  &
                    ke6cp,ke7cp,ke8cp,ke9cp,jesum4,jesum5,jesum6, inve6j,x)
                IF (ibug2 /= 0) WRITE (99,306) nu,k,x
                
                IF (ABS (x) >= eps) THEN
                  x = x*prod
                  CALL cxk (s,is,kaps,nu,k,ibre,2)
                  IF (ibug2 /= 0) WRITE (99,305) (s(iii),iii = 1,imue)
                  DO  mu = 1,imue
                    cone(mu,n) = cone(mu,n)+x*s(mu)
                  END DO
                END IF
                
              END IF
              
!   K = NU-1
              
              IF (ibre == 2) CYCLE
              
              imue = 4
              IF (ibre == 4) imue = 3
              k = nu-1
              
              IF (ncode == n) THEN
                x=xcode
              ELSE
                itkmo = nud-2
                IF (itrig (ks(1),ks(4),itkmo) == 0) GO TO 23
                IF (itrig (ks(2),ks(3),itkmo) == 0) GO TO 23
                j1(mja) = itkmo
                CALL gensum (je6c,je7c,je8c,je9c,jewc,je6,je7,je8,  &
                    je9,kew,jedel,ledel,esumvr,mep,  &
                    je6p,je7p,je8p,je9p,jeword,nelsum,nebj,neb6j,  &
                    ke6cp,ke7cp,ke8cp,ke9cp,jesum4,jesum5,jesum6, inve6j,x)
              END IF
              
              IF (ibug2 /= 0) WRITE (99,306) nu,k,x
              
              IF (ABS (x) >= eps) THEN
                x = x*prod
                CALL cxk (s,is,kaps,nu,k,ibre,2)
                IF (ibug2 /= 0) WRITE (99,305) (s(iii),iii = 1,imue)
                DO  mu = 1,imue
                  cone(mu,n) = cone(mu,n)+x*s(mu)
                END DO
              END IF
              
            END IF
            
!   K = NU+1
            
            23             IF ((ibre == 2) .OR. (n == ne2)) CYCLE
            
            ncode = n+1
            xcode = 0.0D 00
            imue = 12
            IF (ibre == 4) imue = 7
            itkmo = nud+2
            
            IF ((itrig (ks(1),ks(4),itkmo) /= 0) .AND.  &
                  (itrig (ks(2),ks(3),itkmo) /= 0)) THEN
              k = nu+1
              j1(mja) = itkmo
              CALL gensum (je6c,je7c,je8c,je9c,jewc,je6,je7,je8,  &
                  je9,kew,jedel,ledel,esumvr,mep,  &
                  je6p,je7p,je8p,je9p,jeword,nelsum,nebj,neb6j,  &
                  ke6cp,ke7cp,ke8cp,ke9cp,jesum4,jesum5,jesum6, inve6j,x)
              xcode = x
              IF (ibug2 /= 0) WRITE (99,306) nu,k,x
              
              IF (ABS (x) >= eps) THEN
                x = x*prod
                CALL cxk (s,is,kaps,nu,k,ibre,2)
                IF (ibug2 /= 0) WRITE (99,305) (s(iii),iii = 1,imue)
                DO  mu = 1,imue
                  cone (mu,n) = cone(mu,n)+x*s(mu)
                END DO
              END IF
              
            END IF
            
          END DO
          
        END IF
        
      END DO loop26
    END DO
  END DO
END DO

!   4.4 Insert outside factors

30 IF (ibrd >= 0) THEN
  prodd = ks(1)*ks(4)
  prodd = const/SQRT (prodd)
  IF ((ia1 == ib1) .AND. (ia2 == ib2)) prodd = 0.5D 00*prodd
  DO  n = 1,nd2
    DO  mu = 1,12
      cond(mu,n) = cond(mu,n)*prodd
    END DO
  END DO
END IF

IF (ibre >= 0) THEN
  prode = ks(1)*ks(3)
  prode = -const/SQRT(prode)
  DO  n = 1,ne2
    DO  mu = 1,12
      cone(mu,n) = cone(mu,n)*prode
    END DO
  END DO
END IF

!   5.0 Output results

IF (ibrd >= 0) THEN
  DO  n = 1,nd2
    nu = nd1+2*(n-1)
    itype = 1
    IF (ibrd == 2) itype = 3
    IF (ibrd == 3) itype = 4
    CALL talk (ja,jb,nu,ia1,ia2,ib1,ib2,itype,cond(1,n))
    IF (ibrd > 1) CYCLE
    CALL talk (ja,jb,nu,ia2,ia1,ib2,ib1,itype,cond(2,n))
    CALL talk (ja,jb,nu,ia1,ia2,ib2,ib1,itype,cond(3,n))
    CALL talk (ja,jb,nu,ia2,ia1,ib1,ib2,itype,cond(4,n))
    IF (n == nd2) CYCLE
    nup1 = nu+1
    itype = 2
    CALL talk (ja,jb,nup1,ia1,ia2,ib1,ib2,itype,cond(5,n))
    CALL talk (ja,jb,nup1,ib1,ib2,ia1,ia2,itype,cond(6,n))
    CALL talk (ja,jb,nup1,ia2,ia1,ib2,ib1,itype,cond(7,n))
    CALL talk (ja,jb,nup1,ib2,ib1,ia2,ia1,itype,cond(8,n))
    CALL talk (ja,jb,nup1,ia1,ia2,ib2,ib1,itype,cond(9,n))
    CALL talk (ja,jb,nup1,ib2,ib1,ia1,ia2,itype,cond(10,n))
    CALL talk (ja,jb,nup1,ia2,ia1,ib1,ib2,itype,cond(11,n))
    CALL talk (ja,jb,nup1,ib1,ib2,ia2,ia1,itype,cond(12,n))
  END DO
END IF

IF (ibre < 0) RETURN

DO  n = 1,ne2
  nu = ne1+2*(n-1)
  IF (ibre /= 4) THEN
    itype = 1
    IF (ibre == 2) itype = 3
    CALL talk (ja,jb,nu,ia1,ib2,ib1,ia2,itype,cone(1,n))
    IF (ibre == 2) CYCLE
    CALL talk (ja,jb,nu,ib2,ia1,ia2,ib1,itype,cone(2,n))
    CALL talk (ja,jb,nu,ia1,ib2,ia2,ib1,itype,cone(3,n))
    CALL talk (ja,jb,nu,ib2,ia1,ib1,ia2,itype,cone(4,n))
    IF (n == ne2) CYCLE
    nup1 = nu+1
    itype = 2
    CALL talk (ja,jb,nup1,ia1,ib2,ib1,ia2,itype,cone(5,n))
    CALL talk (ja,jb,nup1,ib1,ia2,ia1,ib2,itype,cone(6,n))
    CALL talk (ja,jb,nup1,ib2,ia1,ia2,ib1,itype,cone(7,n))
    CALL talk (ja,jb,nup1,ia2,ib1,ib2,ia1,itype,cone(8,n))
    CALL talk (ja,jb,nup1,ia1,ib2,ia2,ib1,itype,cone(9,n))
    CALL talk (ja,jb,nup1,ia2,ib1,ia1,ib2,itype,cone(10,n))
    CALL talk (ja,jb,nup1,ib2,ia1,ib1,ia2,itype,cone(11,n))
    CALL talk (ja,jb,nup1,ib1,ia2,ib2,ia1,itype,cone(12,n))
  ELSE
    itype = 5
    CALL talk (ja,jb,nu,ib1,ia1,ib1,ia1,itype,cone(1,n))
    CALL talk (ja,jb,nu,ia1,ib1,ib1,ia1,itype,cone(2,n))
    CALL talk (ja,jb,nu,ia1,ib1,ia1,ib1,itype,cone(3,n))
    IF (n == ne2) CYCLE
    nup1 = nu+1
    itype = 6
    CALL talk (ja,jb,nup1,ia1,ib1,ia1,ib1,itype,cone(4,n))
    CALL talk (ja,jb,nup1,ib1,ia1,ib1,ia1,itype,cone(5,n))
    CALL talk (ja,jb,nup1,ia1,ib1,ib1,ia1,itype,cone(6,n))
    CALL talk (ja,jb,nup1,ib1,ia1,ia1,ib1,itype,cone(7,n))
  END IF
END DO

RETURN

!   6.0 Fault diagnostic prints

98 IF (ibug2 /= 0) WRITE (99,300)
RETURN

300 FORMAT ('BREIT: Spectator quantum numbers not diagonal for',  &
    ' non-interacting shells')
301 FORMAT ('BREIT: Increase second dimension of arrays',  &
    ' COND(MU,N) and CONE(MU,N) to the new value of NUMAX,'  &
    /' (at least ',1I3,').')
302 FORMAT ('BREIT: Subshells ',4I5)
303 FORMAT ('  ND1 ND2 NE1 NE2 IBRD IBRE',6I5)
304 FORMAT ('  Direct NU K recoupling coef ',2I5,1P,d20.9)
305 FORMAT (' S',1P,8D15.7)
306 FORMAT ('  Exchange NU K recoupling coef ',2I5,1P,d20.9)
307 FORMAT ('  Statistical factor ',1P,d20.9)

END SUBROUTINE breit
!***********************************************************************
!                                                                      *

FUNCTION clrx (kappaa,k,kappab)
!                                                                      *
!   The value of CLRX is the 3-j symbol:                               *
!                                                                      *
!                    ( JA        K        JB  )                        *
!                    ( 1/2       0       -1/2 )                        *
!                                                                      *
!   The  K'S are kappa angular quantum numbers. The formula is taken   *
!   from D M Brink and G R Satchler, <Angular Momentum>, second edi-   *
!   tion (Oxford: Clarendon press, 1968), p 138.   The logarithms of   *
!   the first  MFACT  factorials must be available in  COMMON/FACTS/   *
!   for this program to function correctly. Note that  N!  is stored   *
!   in FACT(N+1)                                                       *
!                                                                      *
!   No subroutines called.                                             *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: kappaa
INTEGER, INTENT(IN)                      :: k
INTEGER, INTENT(IN OUT)                  :: kappab
IMPLICIT doubleprecision (a-h, o-z)

INTEGER, PARAMETER :: mfact = 500

COMMON/facts/gam(mfact)

!   Determine the absolute values of the kappas

ka = ABS (kappaa)
kb = ABS (kappab)

!   Perform the triangularity check

IF ((ABS (ka-kb) <= k) .AND. (ka+kb-1 >= k)) THEN
  
!   Triangularity satisfied; compute the 3j coefficient
  
!   Begin with the logarithm of the square of the leading term
  
  exptrm = -LOG (DBLE (ka*kb))
  
!   Compute the logarithm of the square root of the leading term
!   and the factorial part that doesn't depend on the parity of
!   KA+KB+K (the delta factor)
  
  kapkb = ka+kb
  kabkp = kapkb+k
  kamkb = ka-kb
  kbmka = kb-ka
  exptrm = 0.5D 00 *(exptrm+gam(kapkb-k  )+gam(kamkb+k+1)  &
      +gam(kbmka+k+1)-gam(kabkp  +1) )
  
!   The remainder depends on the parity of KA+KB+K
  
  IF (MOD (kabkp,2) == 0) THEN
    
!   Computation for even parity case
    
!   Include the phase factor: a minus sign if necessary
    
    IF (MOD (3*kabkp/2,2) == 0) THEN
      clrx =  1.0D 00
    ELSE
      clrx = -1.0D 00
    END IF
    
!   Include the contribution from the factorials
    
    exptrm = exptrm+gam((kabkp  +2)/2)-gam((kapkb-k  )/2)  &
        -gam((kamkb+k+2)/2)-gam((kbmka+k+2)/2)
    
  ELSE
    
!   Computation for odd parity case
    
!   Include the phase factor: a minus sign if necessary
    
    IF (MOD ((3*kabkp-1)/2,2) == 0) THEN
      clrx =  1.0D 00
    ELSE
      clrx = -1.0D 00
    END IF
    
!   Include the contribution from the factorials
    
    exptrm = exptrm+gam((kabkp  +1)/2)-gam((kapkb-k+1)/2)  &
        -gam((kamkb+k+1)/2)-gam((kbmka+k+1)/2)
    
  END IF
  
!   Final assembly
  
  clrx = clrx*EXP (exptrm)
  
ELSE
  
!   Triangularity violated; set the coefficient to zero
  
  clrx = 0.0D 00
  
END IF

RETURN

END FUNCTION clrx
!***********************************************************************
!                                                                      *

SUBROUTINE cor (ja,jb,ja1,jb1,ja2,jb2)
!                                                                      *
!   Computes  the  MCP  coefficients.  Equation numbers are those of   *
!   Computer Phys Commun 5 (1973) 263                                  *
!                                                                      *
!   Call(s) to: [LIB92]: CRE, ITRIG, LTAB, MODJ23, MUMDAD, OCON,       *
!                        SETJ, SKRC, SPEAK, KNJ.                       *
!               [NJGRAF]: NJGRAF, GENSUM.                              *
!                                                                      *
!                                           Last update: 02 Nov 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: ja
INTEGER, INTENT(IN OUT)                  :: jb
INTEGER, INTENT(IN)                      :: ja1
INTEGER, INTENT(IN)                      :: jb1
INTEGER, INTENT(IN)                      :: ja2
INTEGER, INTENT(IN)                      :: jb2
IMPLICIT doubleprecision (a-h, o-z)
INTEGER :: pntriq

doubleprecision, PARAMETER ::mangm=60
INTEGER, PARAMETER :: m3mngm=3*mangm
INTEGER, PARAMETER :: mangmp=2*(mangm/3)
INTEGER, PARAMETER :: mtriad=12
INTEGER, PARAMETER :: m6j=20
INTEGER, PARAMETER :: msum=10

doubleprecision, PARAMETER :: eps=1.0D-10

LOGICAL :: free,dsumvr,esumvr,faild,faile

INTEGER, PARAMETER :: idim = 11
DIMENSION cond(idim),cone(idim)

DIMENSION kaps(4),ks(4),nqs(4),ils(4),lls(4),irows(4)
DIMENSION is(4),js(4)

DIMENSION jd6(m3mngm),jd7(m3mngm),jd8(m3mngm),  &
    jd9(mangmp),kdw(6,m6j),lddel(m6j,2),dsumvr(mangm)
DIMENSION jd6p(mangmp),jd7p(mangmp),jd8p(mangmp),jd9p(mangmp),  &
    jdword(6,m6j), ndbj(msum),ndb6j(msum),kd6cp(msum),kd7cp(msum),kd8cp(msum),  &
    kd9cp(msum),jdsum6(mtriad),jdsum4(mtriad,m6j),jdsum5(mtriad,m6j), invd6j(m6j)
DIMENSION je6(m3mngm),je7(m3mngm),je8(m3mngm),  &
    je9(mangmp),kew(6,m6j),ledel(m6j,2),esumvr(mangm)
DIMENSION je6p(mangmp),je7p(mangmp),je8p(mangmp),je9p(mangmp),  &
    jeword(6,m6j), nebj(msum),neb6j(msum),ke6cp(msum),ke7cp(msum),ke8cp(msum),  &
    ke9cp(msum),jesum6(mtriad),jesum4(mtriad,m6j),jesum5(mtriad,m6j), inve6j(m6j)

COMMON/cons/zero,half,tenth,one,two,three,ten  &
    /couple/mja,nja,j1(mangm),j2(mtriad,3),j3(mtriad,3), free(mangm)  &
    /debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6  &
    /l1/jbq1(3,149),jbq2(3,149),jtq1(3),jtq2(3) /m1/nq1(149),nq2(149)  &
    /m2/jjq1(3,149),jjq2(3,149) /m3/jlist(149),klist(149),npeel,ncore  &
    /orb2/ncf,nw,pntriq /orb4/np(149),nak(149)  &
    /terms/nrows,itab(31),jtab(32),ntab(327)

!   1.0  Initialize pointers and flags and set any
!        tables required.

!        In this segment, the array IS points to the
!        full list of orbitals, the array JS to the
!        array JLIST of peel orbital pointers.

!   1.1  Initialization

js(1) = ja1
js(2) = jb1
js(3) = ja2
js(4) = jb2
DO  i = 1,4
  is(i) = jlist(js(i))
  kaps(i) = 2*nak(is(i))
  ks(i) = ABS (kaps(i))
END DO
ia1 = is(1)
ib1 = is(2)
ia2 = is(3)
ib2 = is(4)

kj23 = 0
isnj = 0
faild = .false.
faile = .false.

!   Initialize arrays

DO  j = 1,nw
  DO  k = 1,3
    jbq1(k,j) = 0
    jbq2(k,j) = 0
  END DO
END DO

nbrj = 3*npeel+7
DO  i = 1,nbrj
  free(i) = .false.
END DO

!   2.0 Set tables of quantum numbers of spectator
!       shells.

DO  jj = 1,npeel
  j = jlist(jj)
  IF ((j /= ia1) .AND. (j /= ib1)) THEN
    DO  k = 1,3
      jbq1(k,j) = jjq1(k,j)
    END DO
  END IF
  IF ((j /= ia2) .AND. (j /= ib2)) THEN
    DO  k = 1,3
      jbq2(k,j) = jjq2(k,j)
    END DO
  END IF
  
!   2.1 Examine quantum numbers of spectator
!       shells for orthogonality and
!       exit if found.
  
  IF ((j /= ia1) .AND. (j /= ib1) .AND. (j /= ia2) .AND. (j /= ib2)) THEN
    DO  k = 1,3
      IF (jbq1(k,j) /= jbq2(k,j)) THEN
        IF (ibug2 /= 0) WRITE (99,300)
        GO TO 41
      END IF
    END DO
  END IF
END DO

!   3.1  Set range of the parameter k for Coulomb
!        integrals.
!        Terminate run if buffer store dimension
!        IDIM is too small.

CALL skrc (is,kaps,ks,kd1,kd2,ke1,ke2)
IF ((kd2 == 0) .AND. (ke2 == 0)) GO TO 41
IF ((kd2. GT. idim) .OR. (ke2 > idim)) THEN
  kk = MAX (ke2,kd2)
  WRITE (*,301) kk
  STOP
END IF

IF (kd2 /= 0) THEN
  DO  k = 1,kd2
    cond(k) = zero
  END DO
END IF
IF (ke2 /= 0) THEN
  DO  k = 1,ke2
    cone(k) = zero
  END DO
END IF

nqs(1) = nq1(ia1)
nqs(2) = nq1(ib1)
nqs(3) = nq2(ia2)
nqs(4) = nq2(ib2)

!   3.3  Set parameters of summation over parent
!        (barred) terms in Eq. (5). The array IROWS
!        is formed to point at the list of allowed
!        parents of active shells in the array
!        NTAB.

CALL ltab (is,nqs,ks,irows)

DO  i = 1,4
  ii = irows(i)
  lls(i) = itab(ii)
  ils(i) = jtab(ii)
END DO

!   4.0  Sum contributions over all parent terms
!        permitted by angular momentum and seniority
!        selection rules.

lls1 = lls(1)
IF (lls1 /= 1) free (ja1) = .true.
lls2 = lls(2)
lls3 = lls(3)
lls4 = lls(4)

ls2 = ils(2)
DO  lb1 = 1,lls2
  ls2 = ls2+3
  it12 = ntab(ls2)
  it2 = ks(2)
  it3 = jjq1(3,ib1)
  IF (itrig (it12,it2,it3) == 0) CYCLE
  IF (ABS (ntab(ls2-2)-jjq1(1,ib1)) /= 1) CYCLE
  
  ls1 = ils(1)
  DO  la1 = 1,lls1
    ls1 = ls1+3
    it11 = ntab(ls1)
    it2 = ks(1)
    
    IF (ia1 == ib1) THEN
      
!   Treat IA1 .EQ. IB1 as special case.
      
      it3 = it12
      IF (itrig (it11,it2,it3) == 0) CYCLE
      IF (ABS (ntab(ls1-2)-ntab(ls2-2)) /= 1) CYCLE
      IF (lls2 /= 1) free(nbrj-8) = .true.
    ELSE
      it3 = jjq1(3,ia1)
      IF (itrig (it11,it2,it3) == 0) CYCLE
      IF (ABS (ntab(ls1-2)-jjq1(1,ia1)) /= 1) CYCLE
      IF (lls2 /= 1) free(jb1) = .true.
    END IF
    
    ls4 = ils(4)
    DO  lb2 = 1,lls4
      ls4 = ls4+3
      it14 = ntab(ls4)
      it2 = ks(4)
      it3 = jjq2(3,ib2)
      IF (itrig (it14,it2,it3) == 0) CYCLE
      IF (ABS (ntab(ls4-2)-jjq2(1,ib2)) /= 1) CYCLE
      
      ls3 = ils(3)
      loop35:  DO  la2 = 1,lls3
        ls3 = ls3+3
        it13 = ntab(ls3)
        it2 = ks(3)
        
        IF (ia2 == ib2) THEN
          
!   Treat IA2 .EQ. IB2 as special case.
          
          it3 = it14
          IF (lls4 /= 1) free(nbrj-6) = .true.
          IF (itrig (it13,it2,it3) == 0) CYCLE loop35
          IF (ABS (ntab(ls3-2)-ntab(ls4-2)) /= 1) CYCLE loop35
        ELSE
          it3 = jjq2(3,ia2)
          IF (itrig (it13,it2,it3) == 0) CYCLE loop35
          IF (ABS (ntab(ls3-2)-jjq2(1,ia2)) /= 1) CYCLE loop35
        END IF
        
!   At this point the current parent has been completely defined,
!   and its quantum numbers can now be set. The JTQ arrays must
!   be set if IA1 .EQ. IB1 or IA2 .EQ. IB2. The matrix element
!   should be diagonal in barred quantum numbers.
        
        DO  k = 1,3
          jbq1(k,ia1) = ntab(ls1+k-3)
          jbq2(k,ia2) = ntab(ls3+k-3)
          jtq1(k) = 0
          IF (ib1 == ia1) THEN
            jtq1(k) = ntab(ls2+k-3)
          ELSE
            jbq1(k,ib1) = ntab(ls2+k-3)
          END IF
          jtq2(k) = 0
          IF (ib2 == ia2) THEN
            jtq2(k) = ntab(ls4+k-3)
          ELSE
            jbq2(k,ib2) = ntab(ls4+k-3)
          END IF
          DO  kk = 1,4
            IF (jbq1(k,is(kk)) /= jbq2(k,is(kk))) CYCLE loop35
          END DO
        END DO
        
!   4.1  Evaluate product of 4 CFPs, Eq. (5).
        
        CALL mumdad (is,kaps,prod)
        IF (ABS (prod) < eps) CYCLE loop35
        
!   4.2  Set arrays for defining the recoupling
!        coefficient.
        
        CALL setj (is,js,ks,npeel,kj23)
        
        IF (isnj == 0) THEN
          
!******************** N J G R A F   V e r s i o n **********************
          
!   Set up the arrays and variables for the direct case.
!   J1(NBRJ) ( = J1(MJA) ) is set to (2*KD1+1) so that NJGRAF is
!   called correctly.
          
          IF (kd2 /= 0) THEN
            IF (kd2 > 1) free(nbrj) = .true.
            j1(nbrj) = kd1+kd1+1
            CALL njgraf (recup,-1,faild)
            isnj = 1
            IF (.NOT. faild) THEN
              CALL knj (jd6c,jd7c,jd8c,jd9c,jdwc,jd6,jd7,jd8,jd9,  &
                  kdw,jddel,lddel,dsumvr,mdp,  &
                  jd6p,jd7p,jd8p,jd9p,jdword,ndlsum,ndbj,  &
                  ndb6j,kd6cp,kd7cp,kd8cp,kd9cp, jdsum4,jdsum5,jdsum6,invd6j)
            END IF
          END IF
          
          
!   Set up the arrays and variables for the exchange case.
!   J1(NBRJ) ( = J1(MJA) ) is set to (2*KE1+1) so that NJGRAF is
!   called correctly.
          
          IF (ke2 /= 0) THEN
            CALL modj23
            free(nbrj) = .false.
            IF (ke2 > 1) free(nbrj) = .true.
            j1(nbrj) = ke1+ke1+1
            CALL njgraf (recup,-1,faile)
            isnj = 2
            IF (.NOT. faile) THEN
              CALL knj (je6c,je7c,je8c,je9c,jewc,je6,je7,je8,je9,  &
                  kew,jedel,ledel,esumvr,mep,  &
                  je6p,je7p,je8p,je9p,jeword,nelsum,nebj,  &
                  neb6j,ke6cp,ke7cp,ke8cp,ke9cp, jesum4,jesum5,jesum6,inve6j)
            END IF
          END IF
        END IF
        
!   4.3  Calculate AD, Eq. (6),
!        without the phase factor.
        
        IF ((kd2 /= 0) .AND. (.NOT. faild)) THEN
          kk = kd1-2
          DO  k = 1,kd2
            kk = kk+2
            j1(mja) = kk+kk+1
            CALL gensum (jd6c,jd7c,jd8c,jd9c,jdwc,jd6,jd7,jd8,jd9,  &
                kdw,jddel,lddel,dsumvr,mdp,  &
                jd6p,jd7p,jd8p,jd9p,jdword,ndlsum,ndbj,ndb6j,  &
                kd6cp,kd7cp,kd8cp,kd9cp,jdsum4,jdsum5,jdsum6,invd6j, x)
            IF (ABS (x) >= eps) cond(k)=cond(k)+x*prod
          END DO
        END IF
        
!   4.4  Calculate AE, Eq. (6),
!        without the phase factor.
        
        IF ((ke2 /= 0) .AND. (.NOT. faile)) THEN
          kk = ke1-2
          DO  k = 1,ke2
            kk = kk+2
            j1(mja) = kk+kk+1
            CALL gensum (je6c,je7c,je8c,je9c,jewc,je6,je7,je8,je9,  &
                kew,jedel,ledel,esumvr,mep,  &
                je6p,je7p,je8p,je9p,jeword,nelsum,nebj,neb6j,  &
                ke6cp,ke7cp,ke8cp,ke9cp,jesum4,jesum5,jesum6,inve6j, y)
            IF (ABS(y) >= eps) cone(k) = cone(k)+y*prod
          END DO
        END IF
        
        IF (faild .AND. faile) GO TO 500
        
      END DO loop35
    END DO
  END DO
END DO

!   4.5  Insert factors independent of barred
!        quantum numbers.
!        Output results

!        Begin with common statistical factors, Eq. (5).

500 const = ocon (ia1,ib1,ia2,ib2)

kap1 = nak(ia1)
kap2 = nak(ib1)
kap3 = nak(ia2)
kap4 = nak(ib2)

!   4.6  Compute products of reduced matrix
!        elements, Eq. (7).
!        CRED for direct terms
!        CREE for exchange terms

IF (kd2 /= 0) THEN
  prodd = const/SQRT (DBLE (ks(1)*ks(4)))
  IF (MOD (kd1,2) /= 0) prodd = -prodd
  IF ((ia1 == ib1) .AND. (ia2 == ib2)) prodd = prodd*half
  kk = kd1-2
  DO  k = 1,kd2
    kk = kk+2
    cred = cre (kap1,kk,kap3)*cre (kap2,kk,kap4)
    x = prodd*cond(k)*cred
    IF (ABS (x) >= eps) CALL speak (ja,jb,ia1,ib1,ia2,ib2,kk,x)
  END DO
END IF

IF (ke2 /= 0) THEN
  prode = const/SQRT(DBLE (ks(1)*ks(3)))
  IF (MOD (ke1,2) /= 0) prode = -prode
  prode = -prode
  kk = ke1-2
  DO  k = 1,ke2
    kk = kk+2
    cree = cre (kap1,kk,kap4)*cre (kap2,kk,kap3)
    x = prode*cone(k)*cree
    IF (ABS (x) >= eps) CALL speak (ja,jb,ia1,ib1,ib2,ia2,kk,x)
  END DO
END IF

41 RETURN

300 FORMAT('COR: Spectator quantum numbers not diagonal for',  &
    ' non-interacting shells')
301 FORMAT('COR: Dimension error: reset PARAMETER IDIM to at least ',  &
    1I2,' and recompile.')

END SUBROUTINE cor
!***********************************************************************
!                                                                      *

SUBROUTINE cord (ja,jb,ja1,ipca,jb1)
!                                                                      *
!   Computes the MCP coefficients for contributions involving closed   *
!   shells.  The  standard formulae are given in I P Grant, Advances   *
!   in Physics  19 (1970) 747, Eq. (8.33).  In this segment JA1, JB1   *
!   point to the JLIST array, IA1, IB1 to the full list of orbitals.   *
!                                                                      *
!   Call(s) to: CLRX, SPEAK.                                           *
!                                                                      *
!                                           Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: ja
INTEGER, INTENT(IN OUT)                  :: jb
INTEGER, INTENT(IN OUT)                  :: ja1
INTEGER, INTENT(IN OUT)                  :: ipca
INTEGER, INTENT(IN OUT)                  :: jb1
IMPLICIT doubleprecision (a-h, o-z)

doubleprecision, PARAMETER :: eps = 1.0D-10

COMMON/orb4/np(149),nak(149) /m1/nq1(149),nq2(149)  &
    /m3/jlist(149),klist(149),npeel,ncore

!   Set quantum numbers required.

IF (ipca == 2) THEN
  ia1 = klist(ja1)
ELSE
  ia1 = jlist(ja1)
END IF
ib1 = klist(jb1)

!   Force IA1 to be greater than IB1

IF (ia1 > ib1) THEN
  ns = ia1
  ia1 = ib1
  ib1 = ns
END IF

kap1 = nak(ia1)
j1 = IABS (kap1)
nqs1 = nq1(ia1)

IF (ia1 == ib1) THEN
  
!   Case when IA1 .EQ. IB1
  
  x = DBLE (nqs1*(nqs1-1)/2)
  CALL speak (ja,jb,ia1,ib1,ia1,ib1,0,x)
  numax = j1+j1-2
  IF (numax  <=  0) RETURN
  const = DBLE (nqs1*nqs1/2)
  DO  nu = 2,numax,2
    gam = clrx (kap1,nu,kap1)
    x = -const*gam*gam
    IF (ABS (x) >= eps) CALL speak (ja,jb,ia1,ib1,ia1,ib1,nu,x)
  END DO
  
!   Case when IA1 .NE. IB1
  
ELSE
  
  kap2 = nak(ib1)
  j2 = ABS (kap2)
  nqs2 = nq1(ib1)
  const = DBLE (nqs1*nqs2)
  CALL speak (ja,jb,ia1,ib1,ia1,ib1,0,const)
  numin = ABS (j1-j2)
  numax = j1+j2-1
  IF (kap1*kap2 < 0) numin = numin+1
  DO  nu = numin,numax,2
    gam = clrx (kap1,nu,kap2)
    x = -const*gam*gam
    IF (ABS (x) >= eps) CALL speak (ja,jb,ia1,ib1,ib1,ia1,nu,x)
  END DO
  
END IF

RETURN
END SUBROUTINE cord
!***********************************************************************
!                                                                      *

FUNCTION cre (kap1,k,kap2)
!                                                                      *
!   Computes the relativistic reduced matrix element                   *
!                                                                      *
!                         (j1 || C(K) || j2),                          *
!                                                                      *
!   Eq. (5.15) of I P Grant, Advances in Physics 19 (1970) 762. KAP1,  *
!   KAP2 are the kappa values corresponding to j1, j2.  The triangle   *
!   conditions are tested by the routine CLRX.                         *
!                                                                      *
!   Call(s) to: [LIB92] CLRX.                                          *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************



INTEGER, INTENT(IN OUT)                  :: kap1
INTEGER, INTENT(IN)                      :: k
INTEGER, INTENT(IN OUT)                  :: kap2
IMPLICIT doubleprecision (a-h, o-z)

k1 = ABS (kap1)
dk1k2 = DBLE (4*k1*IABS (kap2))
cre = SQRT (dk1k2)*clrx (kap1,k,kap2)
IF (MOD (k1,2) == 1) cre  = -cre

RETURN
END FUNCTION cre
!***********************************************************************
!                                                                      *

SUBROUTINE cxk (s,is,kaps,nu,k,ibr,iex)
!                                                                      *
!   Computes  the  coefficients of radial integrals in the expansion   *
!   of the effective interaction strength: X(K,IA1,IB1,IA2,IB2).       *
!                                                                      *
!   Input variables:                                                   *
!                                                                      *
!      IS  : Orbital labels                                            *
!      KAPS: Values of 2*kappa                                         *
!      NU  : Order of radial integral                                  *
!      K   : Index of tensor operator                                  *
!      IEX : 1 for direct, 2 for exchange terms                        *
!      IBR : Classifies type of radial integral.There are 4 distinct   *
!            cases:                                                    *
!            IBR = 1 A. All states distinct                            *
!                    B. ((IA .EQ. IB) .AND. (IC .NE. ID)), or          *
!                       ((IA .NE. IB) .AND. (IC .EQ. ID))              *
!                    These give 12 distinct  radial  integrals, with   *
!                    values of K and NU limited only by angular mom-   *
!                    entum and parity                                  *
!            IBR = 2 ((IA .EQ. IC) .AND. (IB .NE. ID)) or              *
!                    ((IA .NE. IC) .AND. (IB .EQ. ID))                 *
!                    This case gives one non-zero integral when K =    *
!                    NU is ODD                                         *
!            IBR = 3 ((IA .EQ. IC) .AND. (IB .EQ. ID)) AND             *
!                    (IA .NE. IB)                                      *
!                    Integrals of magnetic F-type when K = NU is odd   *
!            IBR = 4 ((IA .EQ. ID) .AND. (IB .EQ. IC)) gives 3  mag-   *
!                    netic G-type integrals and  four  H-TYPE  inte-   *
!                    grals                                             *
!                                                                      *
!   Output:                                                            *
!                                                                      *
!      S   : Coefficients S(MU) MU = 1,12                              *
!                                                                      *
!                                                                      *
!   Call(s) to: [LIB92] CRE.                                           *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************


doubleprecision, INTENT(OUT)             :: s(12)
INTEGER, INTENT(IN)                      :: is(4)
INTEGER, INTENT(IN)                      :: kaps(4)
INTEGER, INTENT(IN OUT)                  :: nu
INTEGER, INTENT(IN)                      :: k
INTEGER, INTENT(IN OUT)                  :: ibr
INTEGER, INTENT(IN OUT)                  :: iex
IMPLICIT doubleprecision (a-h, o-z)



!   1.0  Initialization

DO  mu = 1,12
  s(mu) = 0.0D 00
END DO

ia = is(1)
ib = is(2)
ic = is(3)
id = is(4)
ka = kaps(1)/2
kb = kaps(2)/2
kc = kaps(3)/2
kd = kaps(4)/2
IF (iex /= 2) GO TO 2
kk = kd
ik = id
kd = kc
id = ic
kc = kk
ic = ik
2 SELECT CASE ( ibr )
  CASE (    1)
    GO TO 3
  CASE (    2)
    GO TO 8
  CASE (    3)
    GO TO 11
  CASE (    4)
    GO TO 12
END SELECT
GO TO 17

!   2.0  IBR = 1 --- The general case

3 CONTINUE
IF (nu-k < 0) THEN
  GO TO     7
ELSE IF (nu-k == 0) THEN
  GO TO     4
ELSE
  GO TO     6
END IF

!   2.1  NU = K .GT. 0

4 CONTINUE
s(1) = -(ka+kc)*(kd+kb)
IF (k == 0) GO TO 16
d = k*(k+1)
h = cre (ka,k,kc)*cre(kb,k,kd)
IF (MOD (k,2) /= 0) h = -h
s(1) = s(1)*h/d
DO  mu = 2,4
  s(mu) = s(1)
END DO
RETURN

!   2.2  NU = K+1

6 CONTINUE
dk1 = kc-ka
dk2 = kd-kb
fk = k
gk = k+1
g1 = dk1-gk
g2 = dk1+gk
g3 = dk2-gk
g4 = dk2+gk
kk = k+k+1
h = cre (ka,k,kc)*cre(kb,k,kd)
IF (MOD (k,2) /= 0) h = -h
a = h*fk/gk/DBLE (kk*(kk+2))
s(1) = a*g1*g3
s(2) = a*g2*g4
s(3) = a*g1*g4
s(4) = a*g2*g3
RETURN

!   2.2  NU = K-1

7 CONTINUE
dk1 = kc-ka
dk2 = kd-kb
fk = k
gk = k+1
f1 = dk1-fk
f2 = dk1+fk
f3 = dk2-fk
f4 = dk2+fk
g1 = dk1-gk
g2 = dk1+gk
g3 = dk2-gk
g4 = dk2+gk
kk = k+k+1
h = cre (ka,k,kc)*cre(kb,k,kd)
IF (MOD (k,2) /= 0) h = -h
a = h*gk/fk/DBLE (kk*(kk-2))
s(1) = a*f2*f4
s(2) = a*f1*f3
s(3) = a*f2*f3
s(4) = a*f1*f4
b = h/DBLE (kk*kk)
s(5) = b*f2*g3
s(6) = b*f4*g1
s(7) = b*f1*g4
s(8) = b*f3*g2
s(9) = b*f2*g4
s(10) = b*f3*g1
s(11) = b*f1*g3
s(12) = b*f4*g2
RETURN

!   3.0  IBR = 2  Degenerate case: only one non-zero R-integral

8 CONTINUE
IF ((ia == ic) .AND. (ib /= id)) GO TO 10
IF ((ia /= ic) .AND. (ib == id)) GO TO 9
GO TO 17

9 ik = ib
ib = ia
ia = ik
ik = id
id = ic
ic = ik

kk = kb
kb = ka
ka = kk
kk = kd
kd = kc
kc = kk

10 IF (MOD (k,2) /= 1) RETURN
dk = k*(k+1)
h = cre (ka,k,kc)*cre(kb,k,kd)/dk
s(1) = h*DBLE (4*ka*(kb+kd))
RETURN

!   4.0  IBR = 3. Direct magnetic F-integrals

11 CONTINUE
IF ((ia /= ic) .OR. (ib /= id)) GO TO 17
IF (MOD (k,2) /= 1) RETURN
dk = k*(k+1)
h = cre(ka,k,ka)*cre(kb,k,kb)/dk
s(1) = h*DBLE (16*ka*kb)
RETURN

!   5.0   IBR = 4. Exchange magnetic G- and H-integrals

12 CONTINUE
IF ((ia /= id) .OR. (ib /= ic) )GO TO 17
IF (nu-k < 0) THEN
  GO TO    15
ELSE IF (nu-k == 0) THEN
  GO TO    13
ELSE
  GO TO    14
END IF

!   5.1  NU = K

13 CONTINUE
s(1) = DBLE (ka+kb)*cre(ka,k,kb)
ip = ABS (ka)-ABS (kb)+k+1
s(1) = s(1)*s(1)/DBLE(k*(k+1))
IF (MOD (ip,2) /= 0) s(1) = -s(1)
s(3) = s(1)
s(2) = s(1)+s(1)
RETURN

!   5.2  NU = K+1

14 CONTINUE
dk = kb-ka
gk = k+1
fk = k
g1 = dk+gk
g2 = dk-gk
kk = k+k+1
h = cre (ka,k,kb)**2
IF (ka*kb < 0) h = -h
a = h*fk/gk/DBLE(kk*(kk+2))
s(1) = -a*g1*g1
s(2) = -2.0D 00*a*g1*g2
s(3) = -a*g2*g2
RETURN

!   5.3  NU = K-1

15 CONTINUE
dk = kb-ka
fk = k
gk = k+1
f1 = dk+fk
f2 = dk-fk
g1 = dk+gk
g2 = dk-gk
kk = k+k+1
h = cre (ka,k,kb)**2
IF (ka*kb < 0) h = -h
a = h*gk/fk/DBLE (kk*(kk-2))
s(1) = -a*f2*f2
s(2) = -2.0D 00*a*f1*f2
s(3) = -a*f1*f1
b = h/DBLE (kk*kk)
b = b+b
s(4) = -b*f1*g2
!     S(5) = S(4)
s(5) = -b*f2*g1
s(6) = -b*f1*g1
s(7) = -b*f2*g2
RETURN

!   6.0  Special cases and errors

!   Illegal zero value of K in Type 1

16 WRITE (*,300) is(1),is(2),is(3),is(4),nu,ibr,iex
STOP

!   Illegal combination of states in Type 3 or 4

17 WRITE (*,301) ibr,is(1),is(2),is(3),is(4),nu,k,iex
STOP

300 FORMAT ('CXK: Illegal value K = 0 -' /1X,4I3,2X,i3,2X,2I2)
301 FORMAT ('CXK: Type ',i2,'-' /1X,i2,3X,4I3,2X,2I3,2X,i2)

END SUBROUTINE cxk
!***********************************************************************
!                                                                      *

SUBROUTINE fixj (ja1,ja2,ka,is,ks,ns,kj23)
!                                                                      *
!   Sets up the arrays J1, J2, J3 required by the recoupling package   *
!   NJSYM.                                                             *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN)                      :: ja1
INTEGER, INTENT(IN)                      :: ja2
INTEGER, INTENT(IN)                      :: ka
INTEGER, INTENT(IN)                      :: is(2)
INTEGER, INTENT(IN)                      :: ks(2)
INTEGER, INTENT(IN)                      :: ns
INTEGER, INTENT(IN OUT)                  :: kj23
INTEGER, PARAMETER :: mangm=60
INTEGER, PARAMETER :: mtriad=12

LOGICAL :: free



COMMON/m0/jjc1(149),jjc2(149) /m2/jjq1(3,149),jjq2(3,149)  &
    /m3/jlist(149),klist(149),npeel,ncore  &
    /couple/mja,nja,j1(mangm),j2(mtriad,3),j3(mtriad,3), free(mangm)  &
    /l1/jbq1(3,149),jbq2(3,149),jtq1(3),jtq2(3)

!   Set up the J2 and J3 arrays

nm1 = ns-1
IF (kj23 == 1) GO TO 5
ns1 = ns+1
n2 = ns+ns
n3 = n2+ns

j2(1,1) = n3+2
j2(1,2) = n3+3
j2(1,3) = n3+1
j2(2,1) = ja1
j2(2,2) = n3+1
j2(2,3) = n3-1

j3(1,1) = ja2
j3(1,2) = n3+2
j3(1,3) = n3

IF (ns == 1) GO TO 3

DO  jw = 1,nm1
  jj = jw+2
  j2(jj,1) = ns+jw-1
  j2(jj,2) = jw+1
  j2(jj,3) = ns+jw
  
  jk = jw+1
  j3(jk,1) = n2+jw-2
  j3(jk,2) = jw+1
  j3(jk,3) = n2+jw-1
END DO

j2(3,1) = 1
IF (ja1 == 1) j2(3,1) = n3-1

j3(2,1) = 1
IF (ja2 == 1) j3(2,1) = n3

j2(ns1,3) = n2-1

j3(ns1,1) = n3-2
j3(ns1,2) = n3+3
j3(ns1,3) = n2-1

IF (ja1 == 1) GO TO 2
jaf1 = ja1+1
j2(jaf1,2) = n3-1

2 IF (ja2 == 1) GO TO 4
j3(ja2,2) = n3

IF (ns > 1) GO TO 4
3 j3(2,1) = n3
j3(2,2) = n3+3
j3(2,3) = n3-1

4 CONTINUE

!   Set the J1 array

5 CONTINUE
ii = 0

DO  jw = 1,ns
  ij = jlist(jw)
  ii = ii+1
  j1(ii) = jbq2(3,ij)
END DO

IF (ns == 1) GO TO 9

DO  jw = 1,nm1
  ii = ii+1
  j1(ii) = jjc1(jw)
END DO

DO  jw = 1,nm1
  ii = ii+1
  j1(ii) = jjc2(jw)
END DO

9 CONTINUE
ii = ii+1
ij = is(1)
j1(ii) = jjq1(3,ij)
j1(ii+2) = ks(1)
ii = ii+1
ij = is(2)
j1(ii) = jjq2(3,ij)
j1(ii+2) = ks(2)

ii = ii+3
j1(ii) = ka+ka+1
mja = ii
nja = ns+2

RETURN

END SUBROUTINE fixj
!***********************************************************************
!                                                                      *

FUNCTION irow1 (nelc,ksi)
!                                                                      *
!   Locate the row position of configuration j(**n) in table NTAB.     *
!                                                                      *
!                                           Last update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************

IF ((nelc <= 0) .OR. (nelc > ksi)) THEN
  WRITE (*,300) nelc,ksi
  STOP
END IF

kq1 = nelc-1
kq2 = ksi-kq1
kql = MIN (kq1,kq2)+1
IF (kql == 1) THEN
  irow1 = 1
ELSE
  irow1 = (ksi*(ksi-2))/8+kql
END IF

RETURN

300 FORMAT ('IROW1: ',i3,' electrons in shell with 2j+1 = ',i3)

END FUNCTION irow1
!***********************************************************************
!                                                                      *

FUNCTION itrig (i1,i2,i3)
!                                                                      *
!   The  triangular delta. Input: Values of 2*J+1; Output: 1, IF J'S   *
!   form a triangle; 0, otherwise.                                     *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************

i4 = i2-i3
IF ((i1 >= (ABS (i4)+1)) .AND. ((i1 <= (i2+i3-1)))) THEN
  itrig = 1
ELSE
  itrig = 0
END IF

RETURN
END FUNCTION itrig
!***********************************************************************
!                                                                      *

SUBROUTINE knj (jd6c,jd7c,jd8c,jd9c,jdwc,jd6,jd7,jd8,jd9,kdw,  &
    jddel,lddel,dsumvr,mdp, jd6p,jd7p,jd8p,jd9p,jdword,ndlsum,ndbj,ndb6j,  &
    kd6cp,kd7cp,kd8cp,kd9cp,jdsum4,jdsum5,jdsum6,invd6j)
!                                                                      *
!   This routine stores data for future calls to GENSUM.               *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(OUT)                     :: jd6c
INTEGER, INTENT(OUT)                     :: jd7c
INTEGER, INTENT(OUT)                     :: jd8c
INTEGER, INTENT(OUT)                     :: jd9c
INTEGER, INTENT(OUT)                     :: jdwc
INTEGER, INTENT(OUT)                     :: jd6(m3mngm)
INTEGER, INTENT(OUT)                     :: jd7(m3mngm)
INTEGER, INTENT(OUT)                     :: jd8(m3mngm)
INTEGER, INTENT(OUT)                     :: jd9(mangmp)
INTEGER, INTENT(OUT)                     :: kdw(6,m6j)
INTEGER, INTENT(OUT)                     :: jddel
INTEGER, INTENT(OUT)                     :: lddel(m6j,2)
LOGICAL, INTENT(OUT)                     :: dsumvr(mangm)
INTEGER, INTENT(OUT)                     :: mdp
INTEGER, INTENT(OUT)                     :: jd6p(mangmp)
INTEGER, INTENT(OUT)                     :: jd7p(mangmp)
INTEGER, INTENT(OUT)                     :: jd8p(mangmp)
INTEGER, INTENT(OUT)                     :: jd9p(mangmp)
INTEGER, INTENT(OUT)                     :: jdword(6,m6j)
INTEGER, INTENT(OUT)                     :: ndlsum
INTEGER, INTENT(OUT)                     :: ndbj(msum)
INTEGER, INTENT(OUT)                     :: ndb6j(msum)
INTEGER, INTENT(OUT)                     :: kd6cp(msum)
INTEGER, INTENT(OUT)                     :: kd7cp(msum)
INTEGER, INTENT(OUT)                     :: kd8cp(msum)
INTEGER, INTENT(OUT)                     :: kd9cp(msum)
INTEGER, INTENT(OUT)                     :: jdsum4(mtriad,m6j)
INTEGER, INTENT(OUT)                     :: jdsum5(mtriad,m6j)
INTEGER, INTENT(OUT)                     :: jdsum6(mtriad)
INTEGER, INTENT(OUT)                     :: invd6j(m6j)
IMPLICIT doubleprecision (a-h, o-z)

doubleprecision, PARAMETER ::mangm=60
INTEGER, PARAMETER :: m3mngm=3*mangm
INTEGER, PARAMETER :: mangmp=2*(mangm/3)
INTEGER, PARAMETER :: mtriad=12
INTEGER, PARAMETER :: m6j=20
INTEGER, PARAMETER :: msum=10

LOGICAL :: sumvar




COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp  &
    /sumarg/j6p(mangmp),j7p(mangmp),j8p(mangmp),j9p(mangmp),  &
    jword(6,m6j),nlsum,nbj(msum),nb6j(msum),k6cp(msum),  &
    k7cp(msum),k8cp(msum),k9cp(msum),jsum6(mtriad),  &
    jsum4(mtriad,m6j),jsum5(mtriad,m6j),inv6j(m6j)

jd6c = j6c
jd7c = j7c
jd8c = j8c
jd9c = j9c
jdwc = jwc
jddel = jdel
mdp = mp
ndlsum = nlsum
IF (j6c /= 0) THEN
  DO  i = 1,j6c
    jd6(i) = j6(i)
  END DO
END IF
IF (j7c /= 0) THEN
  DO  i = 1,j7c
    jd7(i) = j7(i)
  END DO
END IF
IF (j8c /= 0) THEN
  DO  i = 1,j8c
    jd8(i) = j8(i)
  END DO
END IF
IF (j9c /= 0) THEN
  DO  i = 1,j9c
    jd9(i) = j9(i)
  END DO
END IF
IF (jwc /= 0) THEN
  DO  i = 1,6
    DO  j = 1,jwc
      kdw(i,j) = kw(i,j)
    END DO
  END DO
  DO  i = 1,jwc
    invd6j(i) = inv6j(i)
  END DO
END IF
IF (jdel /= 0) THEN
  DO  i = 1,2
    DO  j = 1,jdel
      lddel(j,i) = ldel(j,i)
    END DO
  END DO
END IF
IF (mp /= 0) THEN
  DO  i = 1,mp
    dsumvr(i) = sumvar(i)
  END DO
END IF
IF (nlsum /= 0) THEN
  DO  i = 1,nlsum
    ndbj(i) = nbj(i)
    ndb6j(i) = nb6j(i)
    kd6cp(i) = k6cp(i)
    kd7cp(i) = k7cp(i)
    kd8cp(i) = k8cp(i)
    kd9cp(i) = k9cp(i)
  END DO
END IF
DO  i = 1,mangmp
  jd6p(i) = j6p(i)
  jd7p(i) = j7p(i)
  jd8p(i) = j8p(i)
  jd9p(i) = j9p(i)
END DO
DO  i = 1,mtriad
  jdsum6(i) = jsum6(i)
  DO  j = 1,m6j
    jdsum4(i,j) = jsum4(i,j)
    jdsum5(i,j) = jsum5(i,j)
  END DO
END DO
DO  i = 1,6
  DO  j = 1,m6j
    jdword(i,j) = jword(i,j)
  END DO
END DO

RETURN
END SUBROUTINE knj
!***********************************************************************
!                                                                      *

SUBROUTINE ltab (is,nqs,ks,irows)
!                                                                      *
!   locates rows of possible parents of active shell states for acc-   *
!   essing  NTAB. It is assumed that empty shells have been elimina-   *
!   ted from consideration by SUBROUTINE RKCO.                         *
!                                                                      *
!                                           Last update: 15 Oct 1992   *
!                                                                      *
!***********************************************************************



INTEGER, INTENT(IN)                      :: is(4)
INTEGER, INTENT(OUT)                     :: nqs(4)
INTEGER, INTENT(IN)                      :: ks(4)
INTEGER, INTENT(OUT)                     :: irows(4)
DIMENSION  kq(4)

!old  COMMON/TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
!new
nrows = 31
!new

IF (is(1) == is(2)) nqs(1) = nqs(2)-1
IF (is(3) == is(4)) nqs(3) = nqs(4)-1

DO  i = 1,4
  
!   Check that input data are consistent
  
  IF ((nqs(i) <= 0) .OR. (nqs(i) > ks(i))) THEN
    WRITE (*,300) nqs(i),is(i),ks(i)
    STOP
  END IF
  
  kq1 = nqs(i)-1
  kq2 = ks(i)-kq1
  kq(i) = MIN (kq1,kq2)+1
  IF (kq(i) /= 1) THEN
    irows(i) = (ks(i)*(ks(i)-2))/8+kq(i)
  ELSE
    irows(i) = 1
  END IF
  
  IF (irows(i) > nrows) THEN
    WRITE (*,301)
    STOP
  END IF
  
END DO

RETURN

300 FORMAT ('LTAB: ',1I3,' Electrons in shell ',1I3, ' with 2j+1 = ',1I3)
301 FORMAT ('LTAB: Extend COMMON block TERMS')

END SUBROUTINE ltab
!***********************************************************************
!                                                                      *

SUBROUTINE modj23
!                                                                      *
!   Restores  COMMON  block  /COUPLE/ from saved values for exchange   *
!   case.                                                              *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************

PARAMETER (mangm = 60,mtriad = 12)

LOGICAL :: free

COMMON/l2/j2s(mtriad,3),j3s(mtriad,3)  &
    /couple/mja,nja,j1(mangm),j2(mtriad,3),j3(mtriad,3), free(mangm)

ns2 = nja-1
DO  j = 1,3
  DO  i = 1,ns2
    j2(i,j) = j2s(i,j)
    j3(i,j) = j3s(i,j)
  END DO
END DO

i = j3(1,3)
j3(1,3) = j2(1,1)
j2(1,1) = i

RETURN
END SUBROUTINE modj23
!***********************************************************************
!                                                                      *

SUBROUTINE mumdad (is,kaps,x)
!                                                                      *
!   Evaluate the product of 4 CFPs.                                    *
!                                                                      *
!   Call(s) to: [LIB92]: CFP.                                          *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN)                      :: is(2,2)
INTEGER, INTENT(IN)                      :: kaps(2,2)
doubleprecision, INTENT(OUT)             :: x
IMPLICIT doubleprecision (a-h, o-z)

doubleprecision, PARAMETER :: eps = 1.0D-10



COMMON/debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6 /m1/nq1(149),nq2(149)  &
    /m2/jjq1(3,149),jjq2(3,149) /l1/jbq1(3,149),jbq2(3,149),jtq1(3),jtq2(3)

x = 1.0D 00

!   First index

lock = kaps(1,1)
IF (ABS (lock) == 2) GO TO 4
ii = is(1,1)
nel = nq1(ii)
ivp = jbq1(1,ii)
iwp = jbq1(2,ii)
ijp = jbq1(3,ii)-1

!   IA1 .NE. IB1 and IA2 .NE. IB2; use JJQ array.

IF (is(1,1) == is(2,1)) GO TO 1
ivd = jjq1(1,ii)
iwd = jjq1(2,ii)
ijd = jjq1(3,ii)-1
GO TO 2

!   IA1 .EQ. IB1 or IA2 .EQ. IB2; JTQ array needed.

1 nel = nel-1
ivd = jtq1(1)
iwd = jtq1(2)
ijd = jtq1(3)-1
2 CALL cfp (lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,c)
IF (ibug2 /= 0) WRITE (99,300) lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,c
IF (ABS (c) < eps) GO TO 17
x = x*c

4 lock = kaps(2,1)
IF (IABS (lock) == 2) GO TO 8
ii = is(2,1)
nel = nq1(ii)
ivd = jjq1(1,ii)
iwd = jjq1(2,ii)
ijd = jjq1(3,ii)-1
IF (is(1,1) == is(2,1)) GO TO 5
ivp = jbq1(1,ii)
iwp = jbq1(2,ii)
ijp = jbq1(3,ii)-1
GO TO 6
5 ivp = jtq1(1)
iwp = jtq1(2)
ijp = jtq1(3)-1
6 CALL cfp (lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,c)
IF (ibug2 /= 0) WRITE (99,300) lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,c
IF (ABS (c) < eps) GO TO 17
x = x*c
8 CONTINUE

!   Second index

lock = kaps(1,2)
IF (ABS (lock) == 2) GO TO 12
ii = is(1,2)
nel = nq2(ii)
ivp = jbq2(1,ii)
iwp = jbq2(2,ii)
ijp = jbq2(3,ii)-1

!   IA1 .NE. IB1 and IA2 .NE. IB2; use JJQ array.

IF (is(1,2) == is(2,2)) GO TO 9
ivd = jjq2(1,ii)
iwd = jjq2(2,ii)
ijd = jjq2(3,ii)-1
GO TO 10

!   IA1 .EQ. IB1 or IA2 .EQ. IB2; JTQ array needed.

9 nel = nel-1
ivd = jtq2(1)
iwd = jtq2(2)
ijd = jtq2(3)-1
10 CALL cfp (lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,c)
IF (ibug2 /= 0) WRITE (99,300) lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,c
IF (ABS (c) < eps) GO TO 17
x = x*c

12 lock = kaps(2,2)
IF (ABS (lock) == 2) GO TO 16
ii = is(2,2)
nel = nq2(ii)
ivd = jjq2(1,ii)
iwd = jjq2(2,ii)
ijd = jjq2(3,ii)-1
IF (is(1,2) == is(2,2)) GO TO 13
ivp = jbq2(1,ii)
iwp = jbq2(2,ii)
ijp = jbq2(3,ii)-1
GO TO 14
13 ivp = jtq2(1)
iwp = jtq2(2)
ijp = jtq2(3)-1
14 CALL cfp (lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,c)
IF (ibug2 /= 0) WRITE (99,300) lock,nel,ijd,ivd,iwd,ijp,ivp,iwp,c
IF (ABS (c) < eps) GO TO 17
x = x*c
16 CONTINUE
RETURN

17 x = 0.0D 00
RETURN

300 FORMAT ('MUMDAD: CFP ',i3,i4,i7,2I4,i7,2I4,1P,1D19.12)

END SUBROUTINE mumdad
!***********************************************************************
!                                                                      *

FUNCTION ocon (ia1,ib1,ia2,ib2)
!                                                                      *
!   Evaluates the  multiplicative statistical  factor. It is assumed   *
!   that states are ordered so that IA1 .LE. IB1, IA2 .LE. IB2.        *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: ia1
INTEGER, INTENT(IN OUT)                  :: ib1
INTEGER, INTENT(IN OUT)                  :: ia2
INTEGER, INTENT(IN OUT)                  :: ib2
IMPLICIT doubleprecision (a-h, o-z)

COMMON/m1/nq1(149),nq2(149)

wa = DBLE (nq1(ia1)*nq1(ib1))
IF (ia1 == ib1) wa = wa-DBLE (nq1(ia1))
wb = DBLE (nq2(ia2)*nq2(ib2))
IF (ia2 == ib2) wb = wb-DBLE (nq2(ib2))
wc = wa*wb
ocon = SQRT (wc)

!   Set phase factor (-1)**(DELTA P)

lrd1 = MIN (ia2,ib2)+1
lrd2 = MAX (ia2,ib2)
IF (lrd1 > lrd2) THEN
  idr = 0
ELSE
  idr = 1
  DO  k = lrd1,lrd2
    idr = idr+nq2(k)
  END DO
END IF

lld1 = MIN (ia1,ib1)+1
lld2 = MAX (ia1,ib1)
IF (lld1 > lld2) THEN
  idl = 0
ELSE
  idl = 1
  DO  k = lld1,lld2
    idl = idl+nq1(k)
  END DO
END IF

iphas = idr-idl
IF (MOD (iphas,2) /= 0) ocon = -ocon

RETURN
END FUNCTION ocon
!***********************************************************************
!                                                                      *

SUBROUTINE rkco (ja,jb,cor,cord,incor)
!                                                                      *
!   Configurations JA, JB. Analyse the tables of quantum numbers set   *
!   in the COMMON  blocks M0 , M1, M2, M3  to determine all possible   *
!   sets of interacting  orbitals which give a non-vanishing Coulomb   *
!   matrix element,  and  initiates the calculation of coefficients.   *
!   The following conventions are in force: (1) labels 1, 2 refer to   *
!   left, right sides of matrix element respectively;   (2) pointers   *
!   JA1, JB1, JA2, JB2 point to the JLIST array of active  orbitals;   *
!   IA1, IB1, IA2, IB2 point to the complete list of orbitals.         *
!                                                                      *
!   Call(s) to: [LIB92]: COR, CORD, ISPAR, ITJPO, SETQNA, VIJOUT.      *
!                                                                      *
!                                           Last update: 02 Nov 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: ja
INTEGER, INTENT(IN OUT)                  :: jb
REAL, INTENT(IN OUT)                     :: cor
REAL, INTENT(IN OUT)                     :: cord
INTEGER, INTENT(IN OUT)                  :: incor
INTEGER :: pntriq

COMMON/debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6  &
    /dumx/jlis(149),jc1s(149),jc2s(149) /m0/jjc1(149),jjc2(149)  &
    /m1/nq1(149),nq2(149) /m2/jjq1(3,149),jjq2(3,149)  &
    /m3/jlist(149),klist(149),npeel,ncore /orb2/ncf,nw,pntriq  &
    /orb4/np(149),nak(149)

!   The Hamiltonian is an even scalar operator

IF (itjpo (ja) /= itjpo (jb)) RETURN
IF (ispar (ja) /= ispar (jb)) RETURN

CALL setqna (ja,jb)
IF (ibug2 == 1) CALL vijout (ja,jb)

!   1.0 Analyse peel shell interactions

!   1.1 Analyse electron distribution in peel. (The full procedure is
!       needed only if the number of peel orbitals NPEEL .GE. 2)

IF (nw < 1) THEN
  PRINT *, 'RKCO: No subshells.'
  STOP
END IF
IF (npeel == 0) GO TO 48
IF (npeel == 1) GO TO 43

!   Find differences in occupations, NDQ, for each peel orbital in
!   turn and use to set up labels of active orbitals maintaining the
!   convention JA1 .LE. JB1, JA2 .LE. JB2.

idq = 0
ja1 = 0
jb1 = 0
ja2 = 0
jb2 = 0
DO  jw = 1,npeel
  j = jlist(jw)
  ndq = nq1(j) - nq2(j)
  IF (IABS (ndq) > 2) RETURN
  IF (ndq < 0) GO TO 5
  IF (ndq-1 < 0) THEN
    GO TO    10
  ELSE IF (ndq-1 == 0) THEN
    GO TO     1
  ELSE
    GO TO     4
  END IF
  1    IF (ja1 > 0) GO TO 2
  ja1 = jw
  GO TO 3
  2    jb1 = jw
  3    idq = idq+1
  CYCLE
  4    ja1 = jw
  idq = idq+2
  CYCLE
  5    IF (ndq+1 < 0) THEN
    GO TO     9
  ELSE IF (ndq+1 == 0) THEN
    GO TO     6
  ELSE
    GO TO    10
  END IF
  6    IF (ja2 > 0) GO TO 7
  ja2 = jw
  GO TO 8
  7    jb2 = jw
  8    idq = idq+1
  CYCLE
  9    ja2 = jw
  idq = idq+2
END DO

!   1.2 Calculate coefficients for all possible sets of active shells.

!   There are 4 cases, depending on the value of IDQ, the sum of the
!   absolute differences NDQ:

!   1.2.1 IDQ .GT. 4: matrix element null

IF (idq > 4) RETURN
IF (idq == 4) GO TO 12
IF (idq /= 2) GO TO 11
klast = 1
GO TO 16
11 IF (idq /= 0) GO TO 54
IF (ja == jb) GO TO 43
klast = npeel
GO TO 16

!   1.2.2 IDQ .EQ. 4: matrix element uniquely defined

12 IF (jb1 /= 0) GO TO 13
jb1 = ja1
13 IF (jb2 /= 0) GO TO 14
jb2 = ja2
14 IF (ibug2 /= 0) WRITE (99,301) ja1,jb1,ja2,jb2
CALL cor (ja,jb,ja1,jb1,ja2,jb2)
RETURN

!   1.2.3 IDQ .EQ. 2: One orbital fixed each side include all
!                     possible spectators.

!   Also IDQ .EQ. 0 for a matrix element off-diagonal in coupling
!   only. Must sum over all pairs of orbitals excluding core-core
!   terms

16 DO  kwa = 1,klast
  IF (idq == 2) GO TO 17
  ja1 = kwa
  ja2 = kwa
  17    jt1 = ja1
  jt2 = ja2
  it1 = jlist(ja1)
  it2 = jlist(ja2)
  DO  kw = kwa,npeel
    k1 = jlist(kw)
    IF (nq1(k1)*nq2(k1) == 0) CYCLE
    jb1 = kw
    jb2 = kw
    ja1 = jt1
    ja2 = jt2
    
!   Interchange JA1 and JB1 and/or JA2 and JB2 if necessary
    
    IF (ja1-jb1 < 0) THEN
      GO TO    20
    ELSE IF (ja1-jb1 == 0) THEN
      GO TO    19
    END IF
    18       jt3 = jb1
    jb1 = ja1
    ja1 = jt3
    GO TO 20
    19       ib1 = jlist(jb1)
    IF (nq1(ib1) <= 1) CYCLE
    20       IF (ja2-jb2 < 0) THEN
      GO TO    23
    ELSE IF (ja2-jb2 == 0) THEN
      GO TO    22
    END IF
    21       jt3 = jb2
    jb2 = ja2
    ja2 = jt3
    GO TO 23
    22       ib2 = jlist(jb2)
    IF (nq2(ib2) <= 1) CYCLE
    23       IF (ibug2 /= 0) WRITE (99,301) ja1,jb1,ja2,jb2
    CALL cor (ja,jb,ja1,jb1,ja2,jb2)
  END DO
  IF ((idq == 0) .AND. (ncore == 0)) CYCLE
  IF ((ncore == 0) .OR. (nak(it1) /= nak(it2))) RETURN
  
!   This section calculates the terms arising from active electrons
!   which are in closed shells
  
  npeelm = npeel-1
  DO  i = 1,npeel
    jlis(i) = jlist(i)
  END DO
  DO  i = 1,npeelm
    jc1s(i) = jjc1(i)
    jc2s(i) = jjc2(i)
  END DO
  DO  kw = 1,ncore
    ijw = klist(kw)
    DO  i = 1,npeel
      ij = jlist(i)
      IF (ijw < ij) GO TO 29
    END DO
    i = npeel+1
    GO TO 31
    29       im = npeel-i+1
    DO  ii = 1,im
      jlist(npeel+2-ii) = jlist(npeel+1-ii)
      IF (npeel == ii) EXIT
      jjc1(npeel+1-ii) = jjc1(npeel-ii)
      jjc2(npeel+1-ii) = jjc2(npeel-ii)
    END DO
    31       CONTINUE
    IF (i < 3) GO TO 32
    jjc1(i-1) = jjc1(i-2)
    jjc2(i-1) = jjc2(i-2)
    GO TO 33
    32       i1 = jlist(1)
    jjc1(1) = jjq1(3,i1)
    jjc2(1) = jjq2(3,i1)
    33       jlist(i) = ijw
    ja1 = jt1
    IF (jt1 >= i) ja1 = ja1+1
    jb1 = i
    ja2 = jt2
    IF (jt2 >= i) ja2 = ja2+1
    jb2 = i
    IF (ja1-jb1 > 0) THEN
      GO TO    34
    ELSE
      GO TO    35
    END IF
    34       jt3 = jb1
    jb1 = ja1
    ja1 = jt3
    35       CONTINUE
    IF (ja2-jb2 > 0) THEN
      GO TO    36
    ELSE
      GO TO    37
    END IF
    36       jt3 = jb2
    jb2 = ja2
    ja2 = jt3
    37       CONTINUE
    npeel = npeel+1
    IF (ibug2 /= 0) THEN
      npeelm = npeel-1
      WRITE (99,302) ja1,jb1,ja2,jb2,kw,klist(kw)
      WRITE (99,303) (jlist(i),i = 1,npeel)
      WRITE (99,304) (jjc1(i),i = 1,npeelm)
      WRITE (99,305) (jjc2(i),i = 1,npeelm)
    END IF
    CALL cor (ja,jb,ja1,jb1,ja2,jb2)
    npeel = npeel-1
    npeelm = npeel-1
    DO  i = 1,npeel
      jlist(i) = jlis(i)
    END DO
    DO  i = 1,npeelm
      jjc1(i)  = jc1s(i)
      jjc2(i)  = jc2s(i)
    END DO
  END DO
END DO
RETURN

!   1.2.4 IDQ .EQ. 0 - diagonal case. Include all pairs with
!         JA1 = JA2, JB1 = JB2.

43 DO  kw1 = 1,npeel
  k1 = jlist(kw1)
  jb1 = kw1
  jb2 = kw1
  DO  kw2 = 1,kw1
    ja1 = kw2
    IF (ja1 /= jb1) GO TO 44
    IF (nq1(k1) <= 1) CYCLE
    44       ja2 = ja1
    IF (ibug2 /= 0) WRITE (99,301) ja1,jb1,ja2,jb2
    CALL cor (ja,jb,ja1,jb1,ja2,jb2)
  END DO
END DO
48 IF (incor < 1) RETURN
IF (ncore == 0) RETURN

!   2.0 The diagonal case. deal with contributions from core orbitals
!       if INCOR .EQ. 1.

DO  kw1 = 1,ncore
  jb1 = kw1
  jb2 = kw1
  
!   2.1 Calculate contribution from core/core terms
  
  ipca = 2
  DO  kw2 = 1,kw1
    ja1 = kw2
    ja2 = kw2
    IF (ibug2 /= 0) WRITE (99,301) ja1,jb1,ja1,jb1
    CALL cord (ja,jb,ja1,ipca,jb1)
  END DO
  
!   2.2 Calculate contribution from peel/core terms
  
  IF (npeel == 0) CYCLE
  ipca = 1
  DO  kw2 = 1,npeel
    ja1 = kw2
    ja2 = kw2
    IF (ibug2 /= 0) WRITE (99,301) ja1,jb1,ja1,jb1
    CALL cord (ja,jb,ja1,ipca,jb1)
  END DO
END DO
RETURN

!   3.0 Diagnostic print - NW .LT. 1

54 WRITE (*,300)
STOP

300 FORMAT ('RKCO: Error.')
301 FORMAT ('From RKCO:' /10X,' JA1 = ',i3,4X,' JB1 = ',i3,4X,' JA2 = ',i3,  &
    ' JB2 = ',i3)
302 FORMAT ('From RKCO:' /10X,' JA1 = ',i3,4X,' JB1 = ',i3,4X,' JA2 = ',i3,  &
    ' JB2 = ',i3,' K2  = ',i3,   ' KW  = ',i3)
303 FORMAT (1X,'JLIST : ',25I4)
304 FORMAT (1X,'JJC1  : ',25I4)
305 FORMAT (1X,'JJC2  : ',25I4)

END SUBROUTINE rkco
!***********************************************************************
!                                                                      *

SUBROUTINE setj (is,js,ks,ns,kj23)
!                                                                      *
!   Sets the tables required by  the recoupling  coefficient package   *
!   NJGRAF. This routine loads  the COMMON block /COUPLE/ with para-   *
!   meters for the first call  of NJGRAF involving direct integrals.   *
!   Subsequent exchange calls  of NJGRAF must be preceeded by a call   *
!   of MODJ23 to restore these arrays to their correct initial state.  *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN)                      :: is(2,2)
INTEGER, INTENT(IN)                      :: js(2,2)
INTEGER, INTENT(IN)                      :: ks(2,2)
INTEGER, INTENT(IN)                      :: ns
INTEGER, INTENT(OUT)                     :: kj23
INTEGER, PARAMETER :: mangm = 60
INTEGER, PARAMETER :: mtriad = 12

LOGICAL :: free



COMMON/l1/jbq1(3,149),jbq2(3,149),jtq1(3),jtq2(3)  &
    /l2/j2s(mtriad,3),j3s(mtriad,3) /m0/jjc1(149),jjc2(149)  &
    /m2/jjq1(3,149),jjq2(3,149) /m3/jlist(149),klist(149),npeel,ncore  &
    /couple/mja,nja,j1(mangm),j2(mtriad,3),j3(mtriad,3), free(mangm)

!   1.0  Set J1 array

ii = 0
DO  ij = 1,ns
  i = jlist(ij)
  ii = ii+1
  j1(ii) = jbq1(3,i)
END DO
IF (ns == 1) GO TO 4
ns1 = ns-1
DO  i = 1,ns1
  ii = ii+1
  j1(ii) = jjc1(i)
END DO
DO  i = 1,ns1
  ii = ii+1
  j1(ii) = jjc2(i)
END DO
4 CONTINUE
DO  i = 1,2
  ii = ii+1
  ij = is(i,1)
  j1(ii) = jjq1(3,ij)
  IF ((i == 1) .AND. (is(1,1) == is(2,1))) j1(ii) = jtq1(3)
  j1(ii+4) = ks(i,1)
END DO
DO  i = 1,2
  ii = ii+1
  ij = is(i,2)
  j1(ii) = jjq2(3,ij)
  IF ((i == 1) .AND. (is(1,2) == is(2,2))) j1(ii) = jtq2(3)
  j1(ii+4) = ks(i,2)
END DO

!   2.0  Set J2, J3 arrays if not already available

ns2 = MAX (4,ns+2)
IF (kj23 > 0) GO TO 14

DO  i = 4,ns2
  j2(i,1) = ns+i-4
  j2(i,2) = i-2
  j2(i,3) = ns+i-3
  j3(i,1) = j2(i,1)+ns-1
  j3(i,2) = i-2
  j3(i,3) = j2(i,3)+ns-1
END DO
j2(4,1) = 1
j3(4,1) = 1

!   At this stage, the entries in rows corresponding to active
!   shells are set incorrectly.

!   3.0  Set rows 1 through 3

ns3 = 3*ns
j2(1,1) = ns3+5
j2(1,2) = ns3+7
j2(1,3) = ns3+3
j2(2,1) = js(1,1)
j2(2,2) = ns3+3
j2(2,3) = ns3-1
j2(3,1) = js(2,1)
j2(3,2) = ns3+4
j2(3,3) = ns3

j3(1,1) = ns3+7
j3(1,2) = ns3+4
j3(1,3) = ns3+6
j3(2,1) = js(1,2)
j3(2,2) = ns3+5
j3(2,3) = ns3+1
j3(3,1) = js(2,2)
j3(3,2) = ns3+6
j3(3,3) = ns3+2

!   4.0  Set remaining resultants

ij1 = js(1,1)
ij2 = js(2,1)
IF (ij2 > 1) j2(ij2+2,2) = j2(3,3)
IF (ij2 == 1) j2(4,1) = j2(3,3)
IF (ij1 /= ij2) GO TO 8
j2(3,1) = j2(2,3)
GO TO 9

8 IF (ij1 > 1) j2(ij1+2,2) = j2(2,3)
IF (ij1 == 1) j2(4,1) = j2(2,3)

9 ij1 = js(1,2)
ij2 = js(2,2)
IF (ij2 > 1) j3(ij2+2,2) = j3(3,3)
IF (ij2 == 1) j3(4,1) = j3(3,3)
IF (ij1 /= ij2) GO TO 10
j3(3,1) = j3(2,3)
GO TO 11

10 IF (ij1 > 1) j3(ij1+2,2) = j3(2,3)
IF (ij1 == 1) j3(4,1) = j3(2,3)

!   All arrays now set. Put up flag KJ23.

11 kj23 = 1
mja = ns3+7
nja = ns+3

!   5.0  Save J2, J3 and return

DO  j = 1,3
  DO  i = 1,ns2
    j2s(i,j) = j2(i,j)
    j3s(i,j) = j3(i,j)
  END DO
END DO
RETURN

!   6.0  Reset J2, J3 from buffers if KJ23 has been set

14 DO  j = 1,3
  DO  i = 1,ns2
    j2(i,j) = j2s(i,j)
    j3(i,j) = j3s(i,j)
  END DO
END DO
RETURN

END SUBROUTINE setj
!***********************************************************************
!                                                                      *

SUBROUTINE setqna (ja,jb)
!                                                                      *
!   This generates the  arrays  defining  the quantum numbers of the   *
!   states involved in the  matrix  element  linking  configurations   *
!   labelled by JA, JB.                                                *
!                                                                      *
!   Call(s) to: [LIB92]: ICHOP, IQ, JCUP, JQS.                         *
!                                                                      *
!                                           Last update: 30 Oct 1987   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: ja
INTEGER, INTENT(IN OUT)                  :: jb
IMPLICIT doubleprecision (a-h, o-z)
INTEGER :: pntriq

COMMON/m0/jjc1(149),jjc2(149) /m1/nq1(149),nq2(149)  &
    /m2/jjq1(3,149),jjq2(3,149) /m3/jlist(149),klist(149),npeel,ncore  &
    /orb2/ncf,nw,pntriq /orb4/np(149),nak(149)

!   List parameters defining all shells in both configurations, whether
!   participating or not

DO  j = 1,nw
  nq1(j) = iq (j,ja)
  nq2(j) = iq (j,jb)
  DO  k = 1,3
    jjq1(k,j) = jqs (k,j,ja)
    jjq2(k,j) = jqs (k,j,jb)
  END DO
END DO

!   Define coupling schemes: set JLIST array to define those shells
!   which are open in either configuration, and KLIST array to locate
!   the rest. Exclude shells which are empty in both configurations

npeel = 0
ncore = 0
DO  j = 1,nw
  IF ((ichop (j,ja) /= -1) .OR. (ichop (j,jb) /= -1)) THEN
    IF ((ichop (j,ja) == 1) .AND. (ichop (j,jb) == 1)) THEN
      ncore = ncore+1
      klist(ncore) = j
    ELSE
      npeel = npeel+1
      jlist(npeel) = j
    END IF
  END IF
END DO

!   Return if not more than one shell is open

IF (npeel <= 1) RETURN

!   Set arrays of coupling angular momenta interpolating closed
!   shells where necessary. Left hand side first ...

jcnt = 1
jcntop = 0
jw1 = jlist(1)
jw2 = jlist(2)
IF (ichop (jw1,ja) /= 0) THEN
  jjc1(1) = jqs (3,jw2,ja)
  IF (ichop (jw2,ja) == 0) jcntop = 1
ELSE
  jcntop = 1
  IF (ichop (jw2,ja) == 0) THEN
    jjc1(1) = jcup (jcnt,ja)
    jcnt = jcnt+1
  ELSE
    jjc1(1) = jqs (3,jw1,ja)
  END IF
END IF

DO  j = 3,npeel
  jw = jlist(j)
  IF (ichop (jw,ja) /= 0) THEN
    jjc1(j-1) = jjc1(j-2)
  ELSE
    IF (jcntop /= 0) THEN
      jjc1(j-1) = jcup (jcnt,ja)
      jcnt = jcnt+1
    ELSE
      jjc1(j-1) = jqs (3,jw,ja)
    END IF
    jcntop = jcntop+1
  END IF
END DO

!   ... and repeat for right hand side

jcnt = 1
jcntop = 0
jw1 = jlist(1)
jw2 = jlist(2)
IF (ichop (jw1,jb) /= 0) THEN
  jjc2(1) = jqs (3,jw2,jb)
  IF (ichop (jw2,jb) == 0) jcntop = 1
ELSE
  jcntop = 1
  IF (ichop (jw2,jb) == 0) THEN
    jjc2(1) = jcup (jcnt,jb)
    jcnt = jcnt+1
  ELSE
    jjc2(1) = jqs (3,jw1,jb)
  END IF
END IF

DO  j = 3,npeel
  jw = jlist(j)
  IF (ichop (jw,jb) /= 0) THEN
    jjc2(j-1) = jjc2(j-2)
  ELSE
    IF (jcntop /= 0) THEN
      jjc2(j-1) = jcup (jcnt,jb)
      jcnt = jcnt+1
    ELSE
      jjc2(j-1) = jqs (3,jw,jb)
    END IF
    jcntop = jcntop+1
  END IF
END DO

RETURN
END SUBROUTINE setqna
!***********************************************************************
!                                                                      *

SUBROUTINE skrc (is,kaps,ks,kd1,kd2,ke1,ke2)
!                                                                      *
!   Determines the range of the tensor rank k for Coulomb integral.    *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************



INTEGER, INTENT(IN)                      :: is(4)
INTEGER, INTENT(IN)                      :: kaps(4)
INTEGER, INTENT(IN)                      :: ks(4)
INTEGER, INTENT(OUT)                     :: kd1
INTEGER, INTENT(OUT)                     :: kd2
INTEGER, INTENT(OUT)                     :: ke1
INTEGER, INTENT(OUT)                     :: ke2


kd2 = 0
ke2 = 0

!   Direct terms --- KD1 = minimum k , KD2 = number of terms

isd1 = 1
IF (kaps(1)*kaps(3) < 0) isd1 = -1
isd2 = 1
IF (kaps(2)*kaps(4) < 0) isd2 = -1
kd1a = ABS (ks(1)-ks(3))
IF (isd1 < 0) kd1a = kd1a+2
kd1b = ABS (ks(2)-ks(4))
IF (isd2 < 0) kd1b = kd1b+2
IF (MOD ((kd1a-kd1b)/2,2) /= 0) GO TO 1
kd2a = ks(1)+ks(3)-2
IF (isd1 > 0) kd2a = kd2a-2
kd2b = ks(2)+ks(4)-2
IF (isd2 > 0) kd2b = kd2b-2
kd1 = MAX (kd1a,kd1b)/2
kd2 = MIN (kd2a,kd2b)/2
kd2 = (kd2-kd1)/2+1

!   Exchange terms --- KE1 = minimum k , KE2 = number of terms

1 CONTINUE
IF ((is(1) == is(2)) .OR. (is(3) == is(4))) RETURN
ise1 = 1
IF (kaps(1)*kaps(4) < 0) ise1 = -1
ise2 = 1
IF (kaps(2)*kaps(3) < 0) ise2 = -1
ke1a = ABS (ks(1)-ks(4))
IF (ise1 < 0) ke1a = ke1a+2
ke1b = ABS (ks(2)-ks(3))
IF (ise2 < 0) ke1b = ke1b+2
IF (MOD ((ke1a-ke1b)/2,2) /= 0) RETURN
ke2a = ks(1)+ks(4)-2
IF (ise1 > 0) ke2a = ke2a-2
ke2b = ks(2)+ks(3)-2
IF (ise2 > 0) ke2b = ke2b-2
ke1 = MAX (ke1a,ke1b)/2
ke2 = MIN (ke2a,ke2b)/2
ke2 = (ke2-ke1)/2+1

RETURN
END SUBROUTINE skrc
!***********************************************************************
!                                                                      *

SUBROUTINE snrc (is,kaps,ks,nd1,nd2,ne1,ne2,ibrd,ibre)
!                                                                      *
!   Determines the range of tensor rank NU for direct/exchange terms,  *
!   and classifies the types of radial integral.                       *
!                                                                      *
!   Input variables:                                                   *
!                                                                      *
!      IS      : Orbital labels                                        *
!      KAPS    : Values of 2*kappa                                     *
!      KS      : Values of 2*J+1                                       *
!                                                                      *
!   Outputs:                                                           *
!                                                                      *
!      ND1/NE1 : Lowest NU value for direct/exchange types             *
!      ND2/NE2 : Corresponding number of contributing NU values: NU    *
!                = ND1, ND1+2 ,..., ND1+2*(ND2-1) etc                  *
!      IBRD    : Classify types of  radial  integrals  contributing;   *
!      IBRE      negative value implies null contribution              *
!                                                                      *
!                                           Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************



INTEGER, INTENT(IN)                      :: is(4)
INTEGER, INTENT(IN)                      :: kaps(4)
INTEGER, INTENT(IN)                      :: ks(4)
INTEGER, INTENT(OUT)                     :: nd1
INTEGER, INTENT(OUT)                     :: nd2
INTEGER, INTENT(OUT)                     :: ne1
INTEGER, INTENT(OUT)                     :: ne2
INTEGER, INTENT(OUT)                     :: ibrd
INTEGER, INTENT(OUT)                     :: ibre


nd2 = 0
ne2 = 0

!   2.0  Form limits for direct terms

iac = 1
IF ((kaps(1)*kaps(3)) < 0) iac = -1
iad = 1
IF ((kaps(2)*kaps(4)) < 0) iad = -1
nd1 = ABS (ks(1)-ks(3))/2-1
IF (iac == -1) nd1 = nd1+1
IF (nd1 == -1) nd1 = 1
nd1a = ABS (ks(2)-ks(4))/2-1
IF (iad == -1) nd1a = nd1a+1
IF (nd1a == -1) nd1a = 1
IF (MOD (nd1-nd1a,2) == 0) GO TO 1
ibrd = -1
GO TO 2
1 nd2 = ABS (ks(1)+ks(3))/2
IF (iac == -1) nd2 = nd2+1
nd2a = ABS (ks(2)+ks(4))/2
IF (iad == -1) nd2a = nd2a+1
nd1 = MAX (nd1,nd1a)
nd2 = MIN (nd2,nd2a)
nd2 = (nd2-nd1)/2+1

!   2.1  Identify type of radial integrals

ibrd = 1
IF (is(1) == is(3) .AND. is(2) /= is(4)) ibrd = 2
IF (is(1) /= is(3) .AND. is(2) == is(4)) ibrd = 2
IF (is(1) == is(3) .AND. is(2) == is(4)) ibrd = 3

!   3.0  Form limits for exchange terms

2 IF ((is(1) /= is(2)) .AND. (is(3) /= is(4))) GO TO 3
ibre = -1
RETURN
3 CONTINUE
iac = 1
IF ((kaps(1)*kaps(4)) < 0) iac = -1
iad = 1
IF ((kaps(2)*kaps(3)) < 0) iad = -1
ne1 = IABS (ks(1)-ks(4))/2-1
IF (iac == -1) ne1 = ne1+1
IF (ne1 == -1) ne1 = 1
ne1a = ABS (ks(2)-ks(3))/2-1
IF (iad == -1) ne1a = ne1a+1
IF (ne1a == -1) ne1a = 1
IF (MOD (ne1-ne1a,2) == 0) GO TO 4
ibre = -1
RETURN

4 ne2 = ABS (ks(1)+ks(4))/2
IF (iac == -1) ne2 = ne2+1
ne2a = ABS (ks(2)+ks(3))/2
IF (iad == -1) ne2a = ne2a+1
ne1 = MAX (ne1,ne1a)
ne2 = MIN (ne2,ne2a)
ne2 = (ne2-ne1)/2+1

!   3.1  Identify type of radial integrals

ibre = 1
IF ((is(1) == is(4)) .AND. (is(2) /= is(3))) ibre = 2
IF ((is(1) /= is(4)) .AND. (is(2) == is(3))) ibre = 2
IF ((is(1) == is(3)) .AND. (is(2) == is(4))) ibre = 4
RETURN

END SUBROUTINE snrc
!***********************************************************************
!                                                                      *
BLOCK DATA term
!                                                                      *
!   Taken from GRASP2.                                                 *
!                                                                      *
!   Table of subshell quantum numbers.    Symmetry  of the table for   *
!   particle/hole configurations is used to compress it.               *
!                                                                      *
!                                          Last updated: 30 Sep 1992   *
!                                                                      *
!***********************************************************************

COMMON/terms/nrows,itab(31),jtab(32),ntab(327)

DATA nrows/31/

!   A row is defined by a subshell angular momentum and an occupation
!   number

!   Each entry ITAB gives the number of terms in a row

!   Each entry JTAB gives the starting location -1 of the first triad
!   in a row

!   Each triad in NTAB is (v,w,2J+1); here v is the seniority,
!   w resolves any degeneracy in the seniority scheme, and J is the
!   subshell total angular momentum

!   Empty subshell or full subshell

DATA (itab(i),i =   1,  1)/ 1/
DATA (jtab(i),i =   1,  1)/ 0/
DATA (ntab(i),i =   1,  3)/ 0,0, 1/

!   s, p-   (j = 1/2)

DATA (itab(i),i =   2,  2)/ 1/
DATA (jtab(i),i =   2,  2)/ 3/
DATA (ntab(i),i =   4,  6)/ 1,0, 2/

!   p, d-   (j = 3/2)

DATA (itab(i),i =   3,  4)/ 1,2/
DATA (jtab(i),i =   3,  4)/ 6,  9/
DATA (ntab(i),i =   7, 15)/ 1,0, 4,  &
    0,0, 1, 2,0, 5/

!  d, f-   (j = 5/2)

DATA (itab(i),i =   5,  7)/ 1,3,3/
DATA (jtab(i),i =   5,  7)/ 15, 18, 27/
DATA (ntab(i),i =  16, 36)/ 1,0, 6,  &
    0,0, 1, 2,0, 5, 2,0, 9,  &
    1,0, 6, 3,0, 4, 3,0,10/

!   f, g-   (j = 7/2)

DATA (itab(i),i =   8, 11)/ 1,4,6,8/
DATA (jtab(i),i =   8, 11)/ 36, 39, 51, 69/
DATA (ntab(i),i =  37, 93)/ 1,0, 8,  &
    0,0, 1, 2,0, 5, 2,0, 9, 2,0,13,  &
    1,0, 8, 3,0, 4, 3,0, 6, 3,0,10, 3,0,12, 3,0,16,  &
    0,0, 1, 2,0, 5, 2,0, 9, 2,0,13,  &
    4,0, 5, 4,0, 9, 4,0,11, 4,0,17/

!   g, h-   (j = 9/2)

DATA (itab(i),i =  12, 16)/ 1,5,10,18,20/
DATA (jtab(i),i =  12, 16)/ 93, 96, 111,141,195/
DATA (ntab(i),i =  94,255)/ 1,0,10,  &
    0,0, 1, 2,0, 5, 2,0, 9, 2,0,13, 2,0,17,  &
    1,0,10, 3,0, 4, 3,0, 6, 3,0, 8, 3,0,10, 3,0,12, 3,0,14, 3,0,16, 3,0,18,  &
    3,0,22, 0,0, 1,  &
    2,0,5, 2,0,9, 2,0,13, 2,0,17,  &
    4,0, 1, 4,0, 5, 4,0, 7, 4,0, 9, 4,1, 9, 4,0,11, 4,0,13, 4,1,13,  &
    4,0,15, 4,0,17, 4,0,19, 4,0,21, 4,0,25, 1,0,10,  &
    3,0, 4, 3,0, 6, 3,0, 8, 3,0,10, 3,0,12, 3,0,14, 3,0,16, 3,0,18, 3,0,22,  &
    5,0, 2, 5,0, 6, 5,0, 8, 5,0,10, 5,0,12, 5,0,14, 5,0,16, 5,0,18,  &
    5,0,20, 5,0,26/

!   h, i-   (j = 11/2)

!   First two rows only

DATA (itab(i),i =  17, 18)/ 1,6/
DATA (jtab(i),i =  17, 19)/ 255,258,277/
DATA (ntab(i),i = 256,276)/ 1,0,12,  &
    0,0, 1, 2,0, 5, 2,0, 9, 2,0,13, 2,0,17, 2,0,21/

!   i, k-   (j = 13/2)

!   First two rows only

DATA (itab(i),i =  23, 24)/ 1,7/
DATA (jtab(i),i =  23, 25)/ 276,279,301/
DATA (ntab(i),i = 277,300)/ 1,0,14,  &
    0,0, 1, 2,0, 5, 2,0, 9, 2,0,13, 2,0,17, 2,0,21, 2,0,25/

!   k, l-   (j = 15/2)

!   First two rows only

DATA (itab(i),i =  30, 31)/ 1,8/
DATA (jtab(i),i =  30, 32)/ 300,303,328/
DATA (ntab(i),i = 301,327)/ 1,0,16,  &
    0,0, 1, 2,0, 5, 2,0, 9, 2,0,13, 2,0,17, 2,0,21, 2,0,25, 2,0,29/

END
!***********************************************************************
!                                                                      *

SUBROUTINE tnsrjj (ka,iopar,ja,jb,ia1,ia2,vshell)
!                                                                      *
!   The  main  program for evaluating the reduced matrix elements of   *
!   a one particle operator for configurations in jj-coupling.         *
!                                                                      *
!   Call(s) to: [LIB92]: CFP, FIXJ, GENSUM, ICHOP, IROW1, ISPAR,       *
!                        ITJPO, ITRIG, SETQNA, VIJOUT.                 *
!               [NJGRAF]: NJGRAF.                                      *
!                                                                      *
!                                           Last update: 02 Nov 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN)                      :: ka
INTEGER, INTENT(IN)                      :: iopar
INTEGER, INTENT(IN OUT)                  :: ja
INTEGER, INTENT(IN OUT)                  :: jb
INTEGER, INTENT(OUT)                     :: ia1
INTEGER, INTENT(OUT)                     :: ia2
doubleprecision, INTENT(OUT)             :: vshell(149)
IMPLICIT doubleprecision (a-h,o-z)
INTEGER :: pntriq

doubleprecision, PARAMETER ::mangm=60
INTEGER, PARAMETER :: m3mngm=3*mangm
INTEGER, PARAMETER :: mangmp=2*(mangm/3)
INTEGER, PARAMETER :: mtriad=12
INTEGER, PARAMETER :: m6j=20
INTEGER, PARAMETER :: msum=10

doubleprecision, PARAMETER :: eps=1.0D-10

LOGICAL :: free,sumvar,fail


DIMENSION is(2),ks(2)

COMMON/argu/j6c,j7c,j8c,j9c,jwc,j6(m3mngm),j7(m3mngm),j8(m3mngm),  &
    j9(mangmp),kw(6,m6j),jdel,ldel(m6j,2),sumvar(mangm),mp  &
    /cons/zero,half,tenth,one,two,three,ten  &
    /couple/mja,nja,j1(mangm),j2(mtriad,3),j3(mtriad,3), free(mangm)  &
    /dumx/jlis(149),jc1s(149),jc2s(149)  &
    /debug/ibug1,ibug2,ibug3,ibug4,ibug5,ibug6  &
    /l1/jbq1(3,149),jbq2(3,149),jtq1(3),jtq2(3) /m0/jjc1(149),jjc2(149)  &
    /m1/nq1(149),nq2(149) /m2/jjq1(3,149),jjq2(3,149)  &
    /m3/jlist(149),klist(149),npeel,ncore /orb2/ncf,nw,pntriq  &
    /orb4/np(149),nak(149)
COMMON/sumarg/j6p(mangmp),j7p(mangmp),j8p(mangmp),j9p(mangmp),  &
    jword(6,m6j),nlsum,nbj(msum),nb6j(msum),k6cp(msum),  &
    k7cp(msum),k8cp(msum),k9cp(msum),jsum6(mtriad),  &
    jsum4(mtriad,m6j),jsum5(mtriad,m6j),inv6j(m6j)  &
    /terms/nrows,itab(31),jtab(32),ntab(327)

ia1 = 0
kk = ka+ka+1
IF (itrig (itjpo (ja),itjpo (jb),kk) == 0) RETURN
IF ((iopar /= 0) .AND. (ispar (ja)*ispar (jb)*iopar /= 1)) RETURN

CALL setqna (ja,jb)
IF (ibug4 /= 0) CALL vijout (ja,jb)

DO  ij = 1,nw
  vshell(ij) = zero
END DO

!   Analyse peel shell interactions

idq = 0
ja1 = 0
ja2 = 0

IF (npeel /= 0) THEN
  
  DO  jw = 1,npeel
    ij = jlist(jw)
    ndq = nq1(ij)-nq2(ij)
    IF (ABS (ndq) > 1) GO TO 39
    IF (ndq > 0) THEN
      ja1 = jw
      idq = idq+1
    ELSE IF (ndq < 0) THEN
      ja2 = jw
      idq = idq+1
    END IF
  END DO
  
  IF (idq > 2) GO TO 39
  
!   Evaluate the array VSHELL
  
!   Then there are two possibilities IDQ = 0 or IDQ = 2
!   if IDQ = 0, then loop over all shells by index ISH
!   if IDQ = 2, then one orbital fixed on each side
  
  ns = npeel
END IF

IF (idq == 2) GO TO 19

!   Loop over shells when IDQ = 0

ish = 0
IF (npeel == 0) GO TO 9
DO  i = 1,npeel
  jlis(i) = jlist(i)
END DO
IF (npeel == 1) GO TO 9
npeelm = npeel-1
DO  i = 1,npeelm
  jc1s(i) = jjc1(i)
  jc2s(i) = jjc2(i)
END DO

!   If ISH .GT. NW, then loop is over and return

9 ish = ish+1
IF (ish > nw) RETURN
IF (ichop (ish,ja) == -1) GO TO 9
IF (ibug6 /= 0) WRITE (99,308) ish
IF (ichop (ish,ja) == 0) GO TO 16

!   Case one --- the ISH-th shell is in the core or in the peel and
!   closed for both sides

i = 1
IF (npeel == 0) GO TO 15
DO  i = 1,npeel
  ij = jlist(i)
  IF (ish < ij) GO TO 11
END DO
i = npeel+1
GO TO 13
11 im = npeel-i+1
DO  ii = 1,im
  jlist(npeel+2-ii) = jlist(npeel+1-ii)
  IF (npeel == ii) EXIT
  jjc1(npeel+1-ii) = jjc1(npeel-ii)
  jjc2(npeel+1-ii) = jjc2(npeel-ii)
END DO
13 CONTINUE
IF (i < 3) GO TO 14
jjc1(i-1) = jjc1(i-2)
jjc2(i-1) = jjc2(i-2)
GO TO 15
14 i1 = jlist(1)
jjc1(1) = jjq1(3,i1)
jjc2(1) = jjq2(3,i1)
15 jlist(i) = ish
ja1 = i
ja2 = i
ns = npeel+1
GO TO 19

!   Case two --- the ISH-th shell is in the peel and open for either
!   side

16 ns = npeel
DO  jw = 1,npeel
  nx = ish-jlist(jw)
  IF (nx == 0) EXIT
END DO
18 ja1 = jw
ja2 = jw

!   Main computation

!     JA1, JA2 are the indices of interacting shells in JLIST
!     IA1, IA2 are the indices of interacting shells in NW

19 ia1 = jlist(ja1)
ia2 = jlist(ja2)
ks1 = 2*ABS (nak(ia1))
ks2 = 2*ABS (nak(ia2))

!   Check triangular condition for the active shells

IF (itrig (ks1,ks2,kk) == 1) GO TO 99
IF (idq == 2) RETURN
GO TO 100

!   Set tables of quantum numbers of non-interacting spectator shells

99 CONTINUE

DO  jw = 1,ns
  ij = jlist(jw)
  IF (ij == ia1) GO TO 23
  DO  k = 1,3
    jbq1(k,ij) = jjq1(k,ij)
  END DO
  
  23   IF (ij == ia2) CYCLE
  DO  k = 1,3
    jbq2(k,ij) = jjq2(k,ij)
  END DO
  IF ((ij == ia1) .OR. (ij == ia2)) CYCLE
  DO  k = 1,3
    IF (jbq1(k,ij) /= jbq2(k,ij)) GO TO 40
  END DO
END DO

!   Loop over parent states

is(1) = ia1
is(2) = ia2
ks(1) = ks1
ks(2) = ks2
val = zero
kj23 = 0
ix = 0
fail = .false.

nelcts = nq2(ia2)
l2 = irow1(nelcts,ks2)
lls2 = itab(l2)
ls2 = jtab(l2)

DO  lb = 1,lls2
  ls2 = ls2+3
  it1 = ntab(ls2)
  it2 = ks2
  it3 = jjq2(3,ia2)
  IF (itrig (it1,it2,it3) == 0) CYCLE
  IF (ABS (ntab(ls2-2)-jjq2(1,ia2)) /= 1) CYCLE
  DO  k = 1,3
    jbq2(k,ia2) = ntab(ls2+k-3)
  END DO
  
  nelcts = nq1(ia1)
  l1 = irow1(nelcts,ks1)
  lls1 = itab(l1)
  ls1 = jtab(l1)
  
  loop33:  DO  la = 1,lls1
    ls1 = ls1+3
    it1 = ntab(ls1)
    it2 = ks1
    it3 = jjq1(3,ia1)
    IF (itrig (it1,it2,it3) == 0) CYCLE loop33
    IF (ABS (ntab(ls1-2)-jjq1(1,ia1)) /= 1) CYCLE loop33
    DO  k = 1,3
      jbq1(k,ia1) = ntab(ls1+k-3)
    END DO
    
    DO  k = 1,3
      IF (jbq1(k,ia1) /= jbq2(k,ia1)) CYCLE loop33
      IF (jbq1(k,ia2) /= jbq2(k,ia2)) CYCLE loop33
    END DO
    
!   Parent shells now defined
    
    CALL fixj (ja1,ja2,ka,is,ks,ns,kj23)
    kj23 = 1
    
    IF (ibug6 /= 0) THEN
      mn1 = mja
      ns1 = nja-1
      WRITE (99,302)
      WRITE (99,303) (j1(j),j = 1,mn1)
      WRITE (99,304)
      DO  jw = 1,ns1
        WRITE (99,305) (j2(jw,k),k = 1,3),(j3(jw,k),k = 1,3)
      END DO
    END IF
    
!   Evaluate recoupling coefficient
    
    IF (ix == 0) THEN
      DO  i = 1,mja
        free(i) = .false.
      END DO
      IF (lls2 /= 1) free(ja2) = .true.
      CALL njgraf (recups,-1,fail)
      ix = 1
      IF (fail) GO TO 501
    END IF
    CALL gensum (j6c,j7c,j8c,j9c,jwc,j6,j7,j8,j9,kw,jdel, ldel,sumvar,mp,  &
        j6p,j7p,j8p,j9p,jword,nlsum,nbj,nb6j,  &
        k6cp,k7cp,k8cp,k9cp,jsum4,jsum5,jsum6,inv6j, recups)
    IF (ibug6 /= 0) WRITE (99,307) recups
    IF (ABS(recups) < eps) CYCLE loop33
    
!   Evaluates 2 CFPs
    
    IF (ks1 == 2) GO TO 31
    ii = ia1
    nel = nq1(ii)
    ivp = jbq1(1,ii)
    iwp = jbq1(2,ii)
    ijp = jbq1(3,ii)-1
    ivd = jjq1(1,ii)
    iwd = jjq1(2,ii)
    ijd = jjq1(3,ii)-1
    CALL cfp (ks1,nel,ijd,ivd,iwd,ijp,ivp,iwp,c)
    IF (ibug6 /= 0) WRITE (99,306) ks1,nel,ijd,ivd,iwd,ijp,ivp, iwp,c
    IF (ABS(c) < eps) CYCLE loop33
    recups = recups*c
    
    31        IF (ks2 == 2) GO TO 32
    ii = ia2
    nel = nq2(ii)
    ivd = jjq2(1,ii)
    iwd = jjq2(2,ii)
    ijd = jjq2(3,ii)-1
    ivp = jbq2(1,ii)
    iwp = jbq2(2,ii)
    ijp = jbq2(3,ii)-1
    CALL cfp (ks2,nel,ijd,ivd,iwd,ijp,ivp,iwp,c)
    IF (ibug6 /= 0) WRITE (99,306) ks2,nel,ijd,ivd,iwd, ijp,ivp,iwp,c
    IF (ABS(c) < eps) CYCLE loop33
    recups = recups*c
    
    32     CONTINUE
    val = val+recups
  END DO loop33
END DO

!   End of loop over parent states

501 IF (idq == 2) GO TO 37

!   IDQ = 0 CASE

vshell(ish) = val*DBLE (nq1(ia1))

!   Loop over all shells when IDQ = 0

100   CONTINUE
IF (npeel == 0) GO TO 9
DO  i = 1,npeel
  jlist(i) = jlis(i)
END DO
IF (npeel == 1) GO TO 9
npeelm = npeel-1
DO  i = 1,npeelm
  jjc1(i)  = jc1s(i)
  jjc2(i)  = jc2s(i)
END DO
GO TO 9

!   IDQ = 2 Case

!       Permutation factor for IDQ = 2

37 CONTINUE
val = val*SQRT (DBLE (nq1(ia1)*nq2(ia2)))
lld1 = MIN (ia1,ia2)+1
lld2 = MAX (ia1,ia2)
idl = 1
IF (ia1 < ia2) idl = 0
DO  k = lld1,lld2
  idl = idl+nq1(k)
END DO
IF (MOD(idl,2) /= 0) val = -val
vshell(1) = val
RETURN

39 IF (ibug6 /= 0) WRITE (99,300)
RETURN
40 IF (ibug6 /= 0) WRITE (99,301)
RETURN

300 FORMAT (' One side has more than one interacting electron')
301 FORMAT (' Spectator quantum numbers not diagonal for non-interact'  &
    ,'ing shells')
302 FORMAT (/' J1')
303 FORMAT (24I5)
304 FORMAT (' J2                   J3')
305 FORMAT (3I5,i10,2I5)
306 FORMAT(' CFP  ',i3,i4,i7,2I4,i7,2I4,1P,d20.9)
307 FORMAT(/' Recoupling coefficient = ',1P,d19.12)
308 FORMAT(//' ISH = ',i3)

END SUBROUTINE tnsrjj
!***********************************************************************
!                                                                      *

SUBROUTINE vijout (ja,jb)
!                                                                      *
!   Prints  out tables of configurational quantum numbers defined by   *
!   SETQNA for current matrix element.                                 *
!                                                                      *
!                                           Last update: 14 Oct 1992   *
!                                                                      *
!***********************************************************************


INTEGER, INTENT(IN OUT)                  :: ja
INTEGER, INTENT(IN OUT)                  :: jb
CHARACTER (LEN=2) :: ic,ij

DIMENSION ij(2),jc(2),ic(2)

COMMON/m0/jjc1(149),jjc2(149) /m1/nq1(149),nq2(149)  &
    /m2/jjq1(3,149),jjq2(3,149) /m3/jlist(149),klist(149),npeel,ncore

DATA ij/'/2','  '/

IF (npeel > 0) THEN
  
!   Identify CSFs
  
  WRITE (99,300) ja,jb
  
!   Print active shell quantum numbers from JLIST table
  
  WRITE (99,301)
  
  DO  j = 1,npeel
    
    jw = jlist(j)
    
    jc(1) = jjq1(3,jw)-1
    IF (MOD (jc(1),2) == 1) THEN
      ic(1) = ij(1)
    ELSE
      jc(1) = jc(1)/2
      ic(1) = ij(2)
    END IF
    
    jc(2) = jjq2(3,jw)-1
    
    IF (MOD (jc(2),2) == 1) THEN
      ic(2) = ij(1)
    ELSE
      jc(2) = jc(2)/2
      ic(2) = ij(2)
    END IF
    
    WRITE (99,302) jw,nq1(jw),jjq1(1,jw),jjq1(2,jw),jc(1),ic(1)  &
        ,nq2(jw),jjq2(1,jw),jjq2(2,jw),jc(2),ic(2)
  END DO
  
!   Print coupling angular momenta if NPEEL .GE. 2
  
  IF (npeel > 2) THEN
    
    WRITE (99,303)
    DO  j = 2,npeel
      
      jc(1) = jjc1(j-1)-1
      
      IF (MOD (jc(1),2) == 1) THEN
        ic(1) = ij(1)
      ELSE
        jc(1) = jc(1)/2
        ic(1) = ij(2)
      END IF
      
      jc(2) = jjc2(j-1)-1
      IF (MOD (jc(2),2) == 1) THEN
        ic(2) = ij(1)
      ELSE
        jc(2) = jc(2)/2
        ic(2) = ij(2)
      END IF
      
      WRITE (99,304) (jc(i),ic(i),i = 1,2)
    END DO
  END IF
  
END IF

WRITE (99,305) ncore

RETURN

300 FORMAT (/'From VIJOUT: CSF ',1I2,35X,'CSF ',1I2)
301 FORMAT (3X,'subshell',4X,'q',4X,'v',2X,'w',2X,'J',  &
    19X,'q',4X,'v',2X,'w',2X,'J')
302 FORMAT (7X,i3,i6,i5,2I3,a2,15X,i3,i5,2I3,a2)
303 FORMAT (' coupling schemes:')
304 FORMAT (14X,i2,a2,27X,i2,a2)
305 FORMAT (' there are ',i3,' inactive closed shells.'/)

END SUBROUTINE vijout
