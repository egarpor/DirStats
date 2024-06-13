      SUBROUTINE RIBESL(X,ALPHA,NB,IZE,B,NCALC)
!-------------------------------------------------------------------
!
!  Adapted from https://www.netlib.org/specfun/ribesl
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
! B     - Working precision output vector of length NB. If the routine
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
!     DBLE, EXP, DGAMMA2, GAMMA, INT, MAX, MIN, REAL, SQRT
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
      DOUBLE PRECISION DGAMMA2, &
       ALPHA,B,CONST,CONV,EM,EMPAL,EMP2AL,EN,ENMTEN,ENSIG, &
       ENTEN,EXPARG,FUNC,HALF,HALFX,ONE,P,PLAST,POLD,PSAVE,PSAVEL, &
       RTNSIG,SUM,TEMPA,TEMPB,TEMPC,TEST,TOVER,TWO,X,XLARGE,ZERO
      DIMENSION B(NB)
!-------------------------------------------------------------------
!  Mathematical constants
!-------------------------------------------------------------------
      DATA ONE,TWO,ZERO,HALF,CONST/1.0D0,2.0D0,0.0D0,0.5D0,1.585D0/
!-------------------------------------------------------------------
!  Machine-dependent parameters
!-------------------------------------------------------------------
      DATA NSIG,XLARGE,EXPARG /16,1.0D4,709.0D0/
      DATA ENTEN,ENSIG,RTNSIG/1.0D308,1.0D16,1.0D-4/
      DATA ENMTEN/8.9D-308/
!-------------------------------------------------------------------
!  Statement functions for conversion
!-------------------------------------------------------------------
      CONV(N) = DBLE(N)
      FUNC(X) = DGAMMA2(X)
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
! Calculate P-sequence until N = NB-1. Check for possible overflow.
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
! To avoid overflow, divide P-sequence by TOVER. Calculate
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
                           B(N+L) = ZERO
  130                   CONTINUE
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
! Normalize. Divide all B(N) by sum.
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
                           IF (TEMPB .NE. ZERO) &
                              TOVER = ENMTEN / TEMPB
                           DO 340 N = 2, NB
                              TEMPA = TEMPA / EMPAL
                              EMPAL = EMPAL + ONE
                              TEMPA = TEMPA * TEMPC
                              IF (TEMPA .LE. TOVER*EMPAL) &
                                 TEMPA = ZERO
                              B(N) = TEMPA + TEMPA*TEMPB / EMPAL
                              IF ((B(N) .EQ. ZERO) &
                                  .AND. (NCALC .GT. N)) NCALC = N-1
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

      DOUBLE PRECISION FUNCTION DGAMMA2(X)
!----------------------------------------------------------------------
!
! This routine calculates the GAMMA function for a real argument X.
!   Computation is based on an algorithm outlined in reference 1.
!   The program uses rational functions that approximate the GAMMA
!   function to at least 20 significant decimal digits.  Coefficients
!   for the approximation over the interval (1,2) are unpublished.
!   Those for the approximation for X .GE. 12 are from reference 2.
!   The accuracy achieved depends on the arithmetic system, the
!   compiler, the intrinsic functions, and proper selection of the
!   machine-dependent constants.
!
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows
! XBIG   - the largest argument for which GAMMA(X) is representable
!          in the machine, i.e., the solution to the equation
!                  GAMMA(XBIG) = beta**maxexp
! XINF   - the largest machine representable floating-point number;
!          approximately beta**maxexp
! EPS    - the smallest positive floating-point number such that
!          1.0+EPS .GT. 1.0
! XMININ - the smallest positive floating-point number such that
!          1/XMININ is machine representable
!
!     Approximate values for some important machines are:
!
!                            beta       maxexp        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! Cyber 180/855
!   under NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-Format   (D.P.)        2          127        34.844
! VAX G-Format   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! Cyber 180/855
!   under NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, etc.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, etc.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-Format   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-Format   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The program returns the value XINF for singularities or
!     when overflow would occur.  The computation is believed
!     to be free of underflow and overflow.
!
!
!  Intrinsic functions required are:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! References: "An Overview of Software Development for Special
!              Functions", W. J. Cody, Lecture Notes in Mathematics,
!              506, Numerical Analysis Dundee, 1975, G. A. Watson
!              (ed.), Springer Verlag, Berlin, 1976.
!
!              Computer Approximations, Hart, Et. Al., Wiley and
!              sons, New York, 1968.
!
!  Latest modification: October 12, 1989
!
!  Authors: W. J. Cody and L. Stoltz
!           Applied Mathematics Division
!           Argonne National Laboratory
!           Argonne, IL 60439
!
!----------------------------------------------------------------------
      INTEGER I,N
      LOGICAL PARITY
      DOUBLE PRECISION &
          C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
!----------------------------------------------------------------------
!  Mathematical constants
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/, &
           SQRTPI/0.9189385332046727417803297D0/ &
           PI/3.1415926535897932384626434D0/
!----------------------------------------------------------------------
!  Machine dependent parameters
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,XINF/1.79D308/
!----------------------------------------------------------------------
!  Numerator and denominator coefficients for rational minimax
!     approximation over (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811D+0, &
             2.47656508055759199108314D+1, &
             -3.79804256470945635097577D+2, &
             6.29331155312818442661052D+2, &
             8.66966202790413211295064D+2, &
             -3.14512729688483675254357D+4, &
             -3.61444134186911729807069D+4, &
             6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1, &
             3.15350626979604161529144D+2, &
             -1.01515636749021914166146D+3, &
             -3.10777167157231109440444D+3, &
              2.25381184209801510330112D+4, &
              4.75584627752788110767815D+3, &
              -1.34659959864969306392456D+5, &
              -1.15132259675553483497211D+5/
!----------------------------------------------------------------------
!  Coefficients for minimax approximation over (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728D-03,8.4171387781295D-04, &
           -5.952379913043012D-04,7.93650793500350248D-04, &
           -2.777777777777681622553D-03, &
           8.333333333333333331554247D-02,5.7083835261D-03/
!----------------------------------------------------------------------
!  Statement functions for conversion between integer and float
!----------------------------------------------------------------------
      CONV(I) = DBLE(I)
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .LE. ZERO) THEN
!----------------------------------------------------------------------
!  Argument is negative
!----------------------------------------------------------------------
            Y = -X
            Y1 = AINT(Y)
            RES = Y - Y1
            IF (RES .NE. ZERO) THEN
                  IF (Y1 .NE. AINT(Y1 * HALF) * TWO) PARITY = .TRUE.
                  FACT = -PI / SIN(PI * RES)
                  Y = Y + ONE
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
!----------------------------------------------------------------------
!  Argument is positive
!----------------------------------------------------------------------
      IF (Y .LT. EPS) THEN
!----------------------------------------------------------------------
!  Argument .LT. EPS
!----------------------------------------------------------------------
            IF (Y .GE. XMININ) THEN
                  RES = ONE / Y
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
         ELSE IF (Y .LT. TWELVE) THEN
            Y1 = Y
            IF (Y .LT. ONE) THEN
!----------------------------------------------------------------------
!  0.0 .LT. argument .LT. 1.0
!----------------------------------------------------------------------
                  Z = Y
                  Y = Y + ONE
               ELSE
!----------------------------------------------------------------------
!  1.0 .LT. argument .LT. 12.0, reduce argument if necessary
!----------------------------------------------------------------------
                  N = INT(Y) - 1
                  Y = Y - CONV(N)
                  Z = Y - ONE
            END IF
!----------------------------------------------------------------------
!  Evaluate approximation for 1.0 .LT. argument .LT. 2.0
!----------------------------------------------------------------------
            XNUM = ZERO
            XDEN = ONE
            DO 260 I = 1, 8
               XNUM = (XNUM + P(I)) * Z
               XDEN = XDEN * Z + Q(I)
  260       CONTINUE
            RES = XNUM / XDEN + ONE
            IF (Y1 .LT. Y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  0.0 .LT. argument .LT. 1.0
!----------------------------------------------------------------------
                  RES = RES / Y1
               ELSE IF (Y1 .GT. Y) THEN
!----------------------------------------------------------------------
!  Adjust result for case  2.0 .LT. argument .LT. 12.0
!----------------------------------------------------------------------
                  DO 290 I = 1, N
                     RES = RES * Y
                     Y = Y + ONE
  290             CONTINUE
            END IF
         ELSE
!----------------------------------------------------------------------
!  Evaluate for argument .GE. 12.0,
!----------------------------------------------------------------------
            IF (Y .LE. XBIG) THEN
                  YSQ = Y * Y
                  SUM = C(7)
                  DO 350 I = 1, 6
                     SUM = SUM / YSQ + C(I)
  350             CONTINUE
                  SUM = SUM/Y - Y + SQRTPI
                  SUM = SUM + (Y - HALF) * LOG(Y)
                  RES = EXP(SUM)
               ELSE
                  RES = XINF
                  GO TO 900
            END IF
      END IF
!----------------------------------------------------------------------
!  Final adjustments and return
!----------------------------------------------------------------------
      IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
  900 DGAMMA2 = RES
      RETURN
! ---------- Last line of DGAMMA2 ----------
      END
