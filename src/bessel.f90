      SUBROUTINE RIBESL(X,ALPHA,NB,IZE,B,NCALC)
!-------------------------------------------------------------------
!
!  Source: https://www.netlib.org/specfun/ribesl
!
!  This routine calculates Bessel functions I SUB(N+ALPHA) (X)
!  for non-negative argument X, and non-negative order N+ALPHA,
!  with or without exponential scaling.
!
!
! Explanation of variables in the calling sequence
!
! X     - Working precision non-negative real argument for which
!         I's or exponentially scaled I's (I*EXP(-X))
!         are to be calculated.  If I's are to be calculated,
!         X must be less than EXPARG (see below).
! ALPHA - Working precision fractional part of order for which
!         I's or exponentially scaled I's (I*EXP(-X)) are
!         to be calculated.  0 .LE. ALPHA .LT. 1.0.
! NB    - Integer number of functions to be calculated, NB .GT. 0.
!         The first function calculated is of order ALPHA, and the
!         last is of order (NB - 1 + ALPHA).
! IZE   - Integer type.  IZE = 1 if unscaled I's are to calculated,
!         and 2 if exponentially scaled I's are to be calculated.
! B     - Working precision output vector of length NB.  If the routine
!         terminates normally (NCALC=NB), the vector B contains the
!         functions I(ALPHA,X) through I(NB-1+ALPHA,X), or the
!         corresponding exponentially scaled functions.
! NCALC - Integer output variable indicating possible errors.
!         Before using the vector B, the user should check that
!         NCALC=NB, i.e., all orders have been calculated to
!         the desired accuracy.  See error returns below.
!
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   beta   = Radix for the floating-point system
!   minexp = Smallest representable power of beta
!   maxexp = Smallest power of beta that overflows
!   it     = Number of bits in the mantissa of a working precision
!            variable
!   NSIG   = Decimal significance desired.  Should be set to
!            INT(LOG10(2)*it+1).  Setting NSIG lower will result
!            in decreased accuracy while setting NSIG higher will
!            increase CPU time without increasing accuracy.  The
!            truncation error is limited to a relative error of
!            T=.5*10**(-NSIG).
!   ENTEN  = 10.0 ** K, where K is the largest integer such that
!            ENTEN is machine-representable in working precision
!   ENSIG  = 10.0 ** NSIG
!   RTNSIG = 10.0 ** (-K) for the smallest integer K such that
!            K .GE. NSIG/4
!   ENMTEN = Smallest ABS(X) such that X/4 does not underflow
!   XLARGE = Upper limit on the magnitude of X when IZE=2.  Bear
!            in mind that if ABS(X)=N, then at least N iterations
!            of the backward recursion will be executed.  The value
!            of 10.0 ** 4 is used on every machine.
!   EXPARG = Largest working precision argument that the library
!            EXP routine can handle and upper limit on the
!            magnitude of X when IZE=1; approximately
!            LOG(beta**maxexp)
!
!
!     Approximate values for some important machines are:
!
!                        beta       minexp      maxexp       it
!
!  CRAY-1        (S.P.)    2        -8193        8191        48
!  Cyber 180/855
!    under NOS   (S.P.)    2         -975        1070        48
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)    2         -126         128        24
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)    2        -1022        1024        53
!  IBM 3033      (D.P.)   16          -65          63        14
!  VAX           (S.P.)    2         -128         127        24
!  VAX D-Format  (D.P.)    2         -128         127        56
!  VAX G-Format  (D.P.)    2        -1024        1023        53
!
!
!                        NSIG       ENTEN       ENSIG      RTNSIG
!
! CRAY-1        (S.P.)    15       1.0E+2465   1.0E+15     1.0E-4
! Cyber 180/855
!   under NOS   (S.P.)    15       1.0E+322    1.0E+15     1.0E-4
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)    16       1.0D+308    1.0D+16     1.0D-4
! IBM 3033      (D.P.)     5       1.0D+75     1.0D+5      1.0D-2
! VAX           (S.P.)     8       1.0E+38     1.0E+8      1.0E-2
! VAX D-Format  (D.P.)    17       1.0D+38     1.0D+17     1.0D-5
! VAX G-Format  (D.P.)    16       1.0D+307    1.0D+16     1.0D-4
!
!
!                         ENMTEN      XLARGE   EXPARG
!
! CRAY-1        (S.P.)   1.84E-2466   1.0E+4    5677
! Cyber 180/855
!   under NOS   (S.P.)   1.25E-293    1.0E+4     741
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)   4.70E-38     1.0E+4      88
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)   8.90D-308    1.0D+4     709
! IBM 3033      (D.P.)   2.16D-78     1.0D+4     174
! VAX           (S.P.)   1.17E-38     1.0E+4      88
! VAX D-Format  (D.P.)   1.17D-38     1.0D+4      88
! VAX G-Format  (D.P.)   2.22D-308    1.0D+4     709
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  In case of an error,  NCALC .NE. NB, and not all I's are
!  calculated to the desired accuracy.
!
!  NCALC .LT. 0:  An argument is out of range. For example,
!     NB .LE. 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE. EXPARG.
!     In this case, the B-vector is not calculated, and NCALC is
!     set to MIN0(NB,0)-1 so that NCALC .NE. NB.
!
!  NB .GT. NCALC .GT. 0: Not all requested function values could
!     be calculated accurately.  This usually occurs because NB is
!     much larger than ABS(X).  In this case, B(N) is calculated
!     to the desired accuracy for N .LE. NCALC, but precision
!     is lost for NCALC .LT. N .LE. NB.  If B(N) does not vanish
!     for N .GT. NCALC (because it is too small to be represented),
!     and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
!     significant figures of B(N) can be trusted.
!
!
! Intrinsic functions required are:
!
!     DBLE, EXP, DGAMMA, GAMMA, INT, MAX, MIN, REAL, SQRT
!
!
! Acknowledgement
!
!  This program is based on a program written by David J.
!  Sookne (2) that computes values of the Bessel functions J or
!  I of real argument and integer order.  Modifications include
!  the restriction of the computation to the I Bessel function
!  of non-negative real argument, the extension of the computation
!  to arbitrary positive order, the inclusion of optional
!  exponential scaling, and the elimination of most underflow.
!  An earlier version was published in (3).
!
! References: "A Note on Backward Recurrence Algorithms," Olver,
!              F. W. J., and Sookne, D. J., Math. Comp. 26, 1972,
!              pp 941-947.
!
!             "Bessel Functions of Real Argument and Integer Order,"
!              Sookne, D. J., NBS Jour. of Res. B. 77B, 1973, pp
!              125-132.
!
!             "ALGORITHM 597, Sequence of Modified Bessel Functions
!              of the First Kind," Cody, W. J., Trans. Math. Soft.,
!              1983, pp. 242-245.
!
!  Latest modification: May 30, 1989
!
!  Modified by: W. J. Cody and L. Stoltz
!               Applied Mathematics Division
!               Argonne National Laboratory
!               Argonne, IL  60439
!
!-------------------------------------------------------------------
      INTEGER IZE,K,L,MAGX,N,NB,NBMX,NCALC,NEND,NSIG,NSTART
!    REAL              GAMMA,
      DOUBLE PRECISION DGAMMA,ALPHA,B,CONST,CONV,EM,EMPAL
      DOUBLE PRECISION EMP2AL,EN,ENMTEN,ENSIG,POLD,PSAVE,PSAVEL
      DOUBLE PRECISION ENTEN,EXPARG,FUNC,HALF,HALFX,ONE,P,PLAST,X,XLARGE
      DOUBLE PRECISION RTNSIG,SUM,TEMPA,TEMPB,TEMPC,TEST,TOVER,TWO,ZERO
      DIMENSION B(NB)
!-------------------------------------------------------------------
!  Mathematical constants
!-------------------------------------------------------------------
!    DATA ONE,TWO,ZERO,HALF,CONST/1.0E0,2.0E0,0.0E0,0.5E0,1.585E0/
      DATA ONE,TWO,ZERO,HALF,CONST/1.0D0,2.0D0,0.0D0,0.5D0,1.585D0/
!-------------------------------------------------------------------
!  Machine-dependent parameters
!-------------------------------------------------------------------
!    DATA NSIG,XLARGE,EXPARG /8,1.0E4,88.0E0/
!    DATA ENTEN,ENSIG,RTNSIG/1.0E38,1.0E8,1.0E-2/
!    DATA ENMTEN/4.7E-38/
      DATA NSIG,XLARGE,EXPARG /16,1.0D4,709.0D0/
      DATA ENTEN,ENSIG,RTNSIG/1.0D308,1.0D16,1.0D-4/
      DATA ENMTEN/8.9D-308/
!-------------------------------------------------------------------
!  Statement functions for conversion
!-------------------------------------------------------------------
!    CONV(N) = REAL(N)
!    FUNC(X) = GAMMA(X)
      CONV(N) = DBLE(N)
      FUNC(X) = DGAMMA(X)
!-------------------------------------------------------------------
! Check for X, NB, OR IZE out of range.
!-------------------------------------------------------------------
      IF ((NB.GT.0) .AND. (X .GE. ZERO) .AND. &
         (ALPHA .GE. ZERO) .AND. (ALPHA .LT. ONE) .AND. &
         (((IZE .EQ. 1) .AND. (X .LE. EXPARG)) .OR. &
          ((IZE .EQ. 2) .AND. (X .LE. XLARGE)))) THEN
!-------------------------------------------------------------------
! Use 2-term ascending series for small X
!-------------------------------------------------------------------
            NCALC = NB
            MAGX = INT(X)
            IF (X .GE. RTNSIG) THEN
!-------------------------------------------------------------------
! Initialize the forward sweep, the P-sequence of Olver
!-------------------------------------------------------------------
                  NBMX = NB-MAGX
                  N = MAGX+1
                  EN = CONV(N+N) + (ALPHA+ALPHA)
                  PLAST = ONE
                  P = EN / X
!-------------------------------------------------------------------
! Calculate general significance test
!-------------------------------------------------------------------
                  TEST = ENSIG + ENSIG
                  IF (2*MAGX .GT. 5*NSIG) THEN
                        TEST = SQRT(TEST*P)
                     ELSE
                        TEST = TEST / CONST**MAGX
                  END IF
                  IF (NBMX .GE. 3) THEN
!-------------------------------------------------------------------
! Calculate P-sequence until N = NB-1.  Check for possible overflow.
!-------------------------------------------------------------------
                     TOVER = ENTEN / ENSIG
                     NSTART = MAGX+2
                     NEND = NB - 1
                     DO 100 K = NSTART, NEND
                        N = K
                        EN = EN + TWO
                        POLD = PLAST
                        PLAST = P
                        P = EN * PLAST/X + POLD
                        IF (P .GT. TOVER) THEN
!-------------------------------------------------------------------
! To avoid overflow, divide P-sequence by TOVER.  Calculate
! P-sequence until ABS(P) .GT. 1.
!-------------------------------------------------------------------
                           TOVER = ENTEN
                           P = P / TOVER
                           PLAST = PLAST / TOVER
                           PSAVE = P
                           PSAVEL = PLAST
                           NSTART = N + 1
   60                      N = N + 1
                              EN = EN + TWO
                              POLD = PLAST
                              PLAST = P
                              P = EN * PLAST/X + POLD
                           IF (P .LE. ONE) GO TO 60
                           TEMPB = EN / X
!-------------------------------------------------------------------
! Calculate backward test, and find NCALC, the highest N
! such that the test is passed.
!-------------------------------------------------------------------
                           TEST = POLD*PLAST / ENSIG
                           TEST = TEST*(HALF-HALF/(TEMPB*TEMPB))
                           P = PLAST * TOVER
                           N = N - 1
                           EN = EN - TWO
                           NEND = MIN0(NB,N)
                           DO 80 L = NSTART, NEND
                              NCALC = L
                              POLD = PSAVEL
                              PSAVEL = PSAVE
                              PSAVE = EN * PSAVEL/X + POLD
                              IF (PSAVE*PSAVEL .GT. TEST) GO TO 90
   80                      CONTINUE
                           NCALC = NEND + 1
   90                      NCALC = NCALC - 1
                           GO TO 120
                        END IF
  100                CONTINUE
                     N = NEND
                     EN = CONV(N+N) + (ALPHA+ALPHA)
!-------------------------------------------------------------------
! Calculate special significance test for NBMX .GT. 2.
!-------------------------------------------------------------------
                     TEST = MAX(TEST,SQRT(PLAST*ENSIG)*SQRT(P+P))
                  END IF
!-------------------------------------------------------------------
! Calculate P-sequence until significance test passed.
!-------------------------------------------------------------------
  110             N = N + 1
                     EN = EN + TWO
                     POLD = PLAST
                     PLAST = P
                     P = EN * PLAST/X + POLD
                  IF (P .LT. TEST) GO TO 110
!-------------------------------------------------------------------
! Initialize the backward recursion and the normalization sum.
!-------------------------------------------------------------------
  120             N = N + 1
                  EN = EN + TWO
                  TEMPB = ZERO
                  TEMPA = ONE / P
                  EM = CONV(N) - ONE
                  EMPAL = EM + ALPHA
                  EMP2AL = (EM - ONE) + (ALPHA + ALPHA)
                  SUM = TEMPA * EMPAL * EMP2AL / EM
                  NEND = N - NB
                  IF (NEND .LT. 0) THEN
!-------------------------------------------------------------------
! N .LT. NB, so store B(N) and set higher orders to zero.
!-------------------------------------------------------------------
                        B(N) = TEMPA
                        NEND = -NEND
                        DO 130 L = 1, NEND
  130                      B(N+L) = ZERO
                     ELSE
                        IF (NEND .GT. 0) THEN
!-------------------------------------------------------------------
! Recur backward via difference equation, calculating (but
! not storing) B(N), until N = NB.
!-------------------------------------------------------------------
                           DO 140 L = 1, NEND
                              N = N - 1
                              EN = EN - TWO
                              TEMPC = TEMPB
                              TEMPB = TEMPA
                              TEMPA = (EN*TEMPB) / X + TEMPC
                              EM = EM - ONE
                              EMP2AL = EMP2AL - ONE
                              IF (N .EQ. 1) GO TO 150
                              IF (N .EQ. 2) EMP2AL = ONE
                              EMPAL = EMPAL - ONE
                              SUM = (SUM + TEMPA*EMPAL) * EMP2AL / EM
  140                      CONTINUE
                        END IF
!-------------------------------------------------------------------
! Store B(NB)
!-------------------------------------------------------------------
  150                   B(N) = TEMPA
                        IF (NB .LE. 1) THEN
                           SUM = (SUM + SUM) + TEMPA
                           GO TO 230
                        END IF
!-------------------------------------------------------------------
! Calculate and Store B(NB-1)
!-------------------------------------------------------------------
                        N = N - 1
                        EN = EN - TWO
                        B(N)  = (EN*TEMPA) / X + TEMPB
                        IF (N .EQ. 1) GO TO 220
                        EM = EM - ONE
                        EMP2AL = EMP2AL - ONE
                        IF (N .EQ. 2) EMP2AL = ONE
                        EMPAL = EMPAL - ONE
                        SUM = (SUM + B(N)*EMPAL) * EMP2AL / EM
                  END IF
                  NEND = N - 2
                  IF (NEND .GT. 0) THEN
!-------------------------------------------------------------------
! Calculate via difference equation and store B(N), until N = 2.
!-------------------------------------------------------------------
                     DO 200 L = 1, NEND
                        N = N - 1
                        EN = EN - TWO
                        B(N) = (EN*B(N+1)) / X +B(N+2)
                        EM = EM - ONE
                        EMP2AL = EMP2AL - ONE
                        IF (N .EQ. 2) EMP2AL = ONE
                        EMPAL = EMPAL - ONE
                        SUM = (SUM + B(N)*EMPAL) * EMP2AL / EM
  200                CONTINUE
                  END IF
!-------------------------------------------------------------------
! Calculate B(1)
!-------------------------------------------------------------------
                  B(1) = TWO*EMPAL*B(2) / X + B(3)
  220             SUM = (SUM + SUM) + B(1)
!-------------------------------------------------------------------
! Normalize.  Divide all B(N) by sum.
!-------------------------------------------------------------------
  230             IF (ALPHA .NE. ZERO) THEN
                     SUM = SUM * FUNC(ONE+ALPHA) * (X*HALF)**(-ALPHA)
                  END IF
                  IF (IZE .EQ. 1) SUM = SUM * EXP(-X)
                  TEMPA = ENMTEN
                  IF (SUM .GT. ONE) TEMPA = TEMPA * SUM
                  DO 260 N = 1, NB
                     IF (B(N) .LT. TEMPA) B(N) = ZERO
                     B(N) = B(N) / SUM
  260             CONTINUE
                  RETURN
!-------------------------------------------------------------------
! Two-term ascending series for small X.
!-------------------------------------------------------------------
               ELSE
                  TEMPA = ONE
                  EMPAL = ONE + ALPHA
                  HALFX = ZERO
                  IF (X .GT. ENMTEN) HALFX = HALF * X
                  IF (ALPHA .NE. ZERO) TEMPA = HALFX**ALPHA /FUNC(EMPAL)
                  IF (IZE .EQ. 2) TEMPA = TEMPA * EXP(-X)
                  TEMPB = ZERO
                  IF ((X+ONE) .GT. ONE) TEMPB = HALFX * HALFX
                  B(1) = TEMPA + TEMPA*TEMPB / EMPAL
                  IF ((X .NE. ZERO) .AND. (B(1) .EQ. ZERO)) NCALC = 0
                  IF (NB .GT. 1) THEN
                     IF (X .EQ. ZERO) THEN
                           DO 310 N = 2, NB
                              B(N) = ZERO
  310                      CONTINUE
                        ELSE
!-------------------------------------------------------------------
! Calculate higher-order functions.
!-------------------------------------------------------------------
                           TEMPC = HALFX
                           TOVER = (ENMTEN + ENMTEN) / X
                           IF (TEMPB .NE. ZERO) TOVER = ENMTEN / TEMPB
                           DO 340 N = 2, NB
                              TEMPA = TEMPA / EMPAL
                              EMPAL = EMPAL + ONE
                              TEMPA = TEMPA * TEMPC
                              IF (TEMPA .LE. TOVER*EMPAL) TEMPA = ZERO
                              B(N) = TEMPA + TEMPA*TEMPB / EMPAL
                              IF ((B(N) .EQ. ZERO) .AND. (NCALC .GT. N))&
                                  NCALC = N-1
  340                      CONTINUE
                     END IF
                  END IF
            END IF
         ELSE
            NCALC = MIN0(NB,0)-1
      END IF
      RETURN
!---------- Last line of RIBESL ----------
      END
