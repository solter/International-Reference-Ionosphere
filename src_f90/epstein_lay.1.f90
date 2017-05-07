! *********************************************************************
! ************************ EPSTEIN FUNCTIONS **************************
! *********************************************************************
! REF:  H. G. BOOKER, J. ATMOS. TERR. PHYS. 39, 619-623, 1977
!       K. RAWER, ADV. SPACE RES. 4, #1, 11-15, 1984
! *********************************************************************

! RAWER  LAYER
REAL FUNCTION  RLAY ( X, XM, SC, HX )
    Y1  = EPTR ( X , SC, HX )
    Y1M = EPTR ( XM, SC, HX )
    Y2M = EPST ( XM, SC, HX )
    RLAY = Y1 - Y1M - (X-XM)*Y2M/SC
    RETURN
END FUNCTION RLAY


! dLAY/dX
REAL FUNCTION D1LAY ( X, XM, SC, HX )
    D1LAY = ( EPST(X,SC,HX) - EPST(XM,SC,HX) ) /  SC
    RETURN
END FUNCTION D1LAY


! d2LAY/dX2
REAL FUNCTION D2LAY ( X, XM, SC, HX )
    D2LAY = EPLA(X,SC,HX) / SC**2
    RETURN
END FUNCTION D2LAY


! TRANSITION
REAL FUNCTION EPTR ( X, SC, HX )
    COMMON/ARGEXP/ARGMAX

    D1 = ( X - HX ) / SC
    IF (ABS(D1) >= ARGMAX) THEN
        IF (D1 > 0.0) THEN
            EPTR = D1
        ELSE
            EPTR = 0.0
        END IF
    ELSE
        EPTR = LOG ( 1. + EXP( D1 ))
    END IF
    RETURN
END FUNCTION EPTR


! STEP
REAL FUNCTION EPST ( X, SC, HX )
    COMMON/ARGEXP/ARGMAX

    D1 = ( X - HX ) / SC
    IF (ABS(D1) >= ARGMAX) THEN
        IF (D1 > 0.0) THEN
            EPST = 1.
        ELSE
            EPST = 0.
        ENDIF
    ELSE
        EPST = 1. / ( 1. + EXP( -D1 ))
    END IF

    RETURN
END FUNCION EPST


! STEP FROM Y1 TO Y2      
REAL FUNCTION EPSTEP ( Y2, Y1, SC, HX, X)
    EPSTEP = Y1 + ( Y2 - Y1 ) * EPST ( X, SC, HX)
    RETURN
END FUNCTION EPSTEP 


! PEAK 
REAL FUNCTION EPLA ( X, SC, HX )
    COMMON/ARGEXP/ARGMAX

    D1 = ( X - HX ) / SC
    IF (ABS(D1) >= ARGMAX) THEN
        EPLA = 0
    ELSE
        D0 = EXP ( D1 )
        D2 = 1. + D0
        EPLA = D0 / ( D2 * D2 )
    END IF
    RETURN
END FUNCTION EPLA 


!----------------------------------------------------------------------
! NORMALIZED ELECTRON DENSITY (N/NMF2) FOR THE MIDDLE IONOSPHERE FROM 
! HME TO HMF2 USING LAY-FUNCTIONS.
!----------------------------------------------------------------------
FUNCTION XE2TO5(H,HMF2,NL,HX,SC,AMP)

    DIMENSION       HX(NL),SC(NL),AMP(NL)

    SUM = 1.0
    DO I = 1,NL
        YLAY = AMP(I) * RLAY( H, HMF2, SC(I), HX(I) )
        zlay = 10.**ylay
        sum = sum*zlay
    end do

    XE2TO5 = sum

    RETURN
END FUNCTION XE2TO5


!----------------------------------------------------------------------
! ELECTRON DENSITY WITH NEW MIDDLE IONOSPHERE
!----------------------------------------------------------------------
REAL FUNCTION XEN(H,HMF2,XNMF2,HME,NL,HX,SC,AMP)
    DIMENSION       HX(NL),SC(NL),AMP(NL)
    !
    IF(H < HMF2) GOTO 100
    XEN = XE1(H)
    RETURN
    100     IF(H < HME) GOTO 200
    XEN = XNMF2 * XE2TO5(H,HMF2,NL,HX,SC,AMP)
    RETURN
    200     XEN = XE6(H)
    RETURN
END FUNCTION XEN


! --------------------------------------------------------------------- 
!   CALCULATES RATIO H0.5/HMF2 FOR HALF-DENSITY POINT (NE(H0.5)=0.5*
!   NMF2) T. GULYAEVA, ADVANCES IN SPACE RESEARCH 7, #6, 39-48, 1987.
!
!       INPUT:  IDAY    DAY OF YEAR
!               XHI     SOLAR ZENITH ANGLE [DEGREE]
!       
!       OUTPUT: GRO     RATIO OF HALF DENSITY HEIGHT TO F PEAK HEIGHT
!               SX      SMOOTHLY VARYING SEASON PARAMTER (SX=1 FOR 
!                       DAY=1; SX=3 FOR DAY=180; SX=2 FOR EQUINOX)
! ---------------------------------------------------------------------
SUBROUTINE ROGUL(IDAY,XHI,SX,GRO)

    common  /const1/humr,dumr

    SX = 2. - COS ( IDAY * dumr )
    XS = ( XHI - 20.*SX )/15.
    GRO = 0.8 - 0.2/( 1. + EXP(XS) )
    ! same as gro=0.6+0.2/(1+exp(-xs))
    RETURN
END SUBROUTINE ROGUL


! --------------------------------------------------------------------
! SOLVES QUADRATIC SYSTEM OF LINEAR EQUATIONS:
!
!       INPUT:  N       NUMBER OF EQUATIONS (= NUMBER OF UNKNOWNS)
!               A(N,N)  MATRIX (LEFT SIDE OF SYSTEM OF EQUATIONS)
!               B(N)    VECTOR (RIGHT SIDE OF SYSTEM)
!
!       OUTPUT: AUS     =.TRUE.   NO SOLUTION FOUND
!                       =.FALSE.  SOLUTION IS IN  A(N,J) FOR J=1,N
! --------------------------------------------------------------------
SUBROUTINE LNGLSN ( N, A, B, AUS)

    DIMENSION :: A(5,5), B(5), AZV(10)
    LOGICAL   :: AUS, break=.false.

    NN = N - 1
    AUS = .FALSE.
    DO K=1,N-1
        IMAX = K
        L    = K
        IZG  = 0
        AMAX = ABS( A(K,K) )
        HSP  = 0.0

        DO WHILE (HSP <= AMAX)
            L = L + 1
            IF (L > N) exit

            HSP = ABS( A(L,K) )
            IF (HSP < 1.E-8) IZG = IZG + 1
        END DO

        IF (ABS(AMAX) < 1.E-10) THEN
            AUS = .TRUE.
            RETURN
        END IF

        IF (IMAX /= K) THEN
            DO L=K,N
                AZV(L+1)  = A(IMAX,L)
                A(IMAX,L) = A(K,L)
                A(K,L)    = AZV(L+1)
            end do
            AZV(1)  = B(IMAX)
            B(IMAX) = B(K)
            B(K)    = AZV(1)
        END IF

        IF (IZG == (N-K)) exit

        AMAX = 1. / A(K,K)
        AZV(1) = B(K) * AMAX

        DO M=K+1,N
            AZV(M+1) = A(K,M) * AMAX
        end do

        DO L=K+1,N
            AMAX = A(L,K)
            IF (ABS(AMAX) < 1.E-8) THEN
                break = .true.
                exit
            END IF

            A(L,K) = 0.0
            B(L) = B(L) - AZV(1) * AMAX
            DO M=K+1,N
                A(L,M) = A(L,M) - AMAX * AZV(M+1)
            end do
        end do
        if (break) exit
    end do

    DO K=N,1,-1
        AMAX = 0.0
        IF (K < N) THEN
            DO L=K+1,N
                AMAX = AMAX + A(K,L) * A(N,L)
            end do
        ENDIF
        IF (ABS(A(K,K)) < 1.E-6) THEN
            A(N,K) = 0.0
        ELSE
            A(N,K) = ( B(K) - AMAX ) / A(K,K)
        ENDIF
    end do

    RETURN
END SUBROUTINE LNGLSN 


! --------------------------------------------------------------------
!   DETERMINES LAY-FUNCTIONS AMPLITUDES FOR A NUMBER OF CONSTRAINTS:
!
!       INPUT:  N       NUMBER OF AMPLITUDES ( LAY-FUNCTIONS)
!               M       NUMBER OF CONSTRAINTS
!               M0      NUMBER OF POINT CONSTRAINTS
!               M1      NUMBER OF FIRST DERIVATIVE CONSTRAINTS
!               HM      F PEAK ALTITUDE  [KM]
!               SC(N)   SCALE PARAMETERS FOR LAY-FUNCTIONS  [KM]
!               HX(N)   HEIGHT PARAMETERS FOR LAY-FUNCTIONS  [KM]
!               W(M)    WEIGHT OF CONSTRAINTS
!               X(M)    ALTITUDES FOR CONSTRAINTS  [KM]
!               Y(M)    LOG(DENSITY/NMF2) FOR CONSTRAINTS
!
!       OUTPUT: VAR(M)  AMPLITUDES
!               SING    =.TRUE.   NO SOLUTION
! ---------------------------------------------------------------------
SUBROUTINE LSKNM ( N, M, M0, M1, HM, SC, HX, W,X,Y,VAR,SING)
    
    LOGICAL         SING
    DIMENSION       VAR(N), HX(N), SC(N), W(M), X(M), Y(M),&
                    BLI(5), ALI(5,5), XLI(5,10)
    
    M01 = M0 + M1
    SCM = 0

    BLI = 0.
    ALI = 0.

    DO I=1,N
        DO K=1,M0
            XLI(I,K) = RLAY( X(K), HM, SC(I), HX(I) )
        end do
        DO K=M0+1,M01
            XLI(I,K) = D1LAY( X(K), HM, SC(I), HX(I) )
        end do
        DO K=M01+1,M
            XLI(I,K) = D2LAY( X(K), HM, SC(I), HX(I) )
        end do
    end do

    BLI(1:N) = matmul( XLI(1:N, 1:M), W(1:M)*Y(1:M) )
    ALI = matmul( XLI(:,:M)*SPREAd(w(:M), 1, 2), transpose(XLI(:,:M)) )

    CALL LNGLSN( N, ALI, BLI, SING )

    IF (.NOT.SING) THEN
        VAR(1:N) = ALI(N, 1:N)
    ENDIF

    RETURN
END SUBROUTINE LSKNM 


!-------------------------------------------------------------------
! CALCULATES AMPLITUDES FOR LAY FUNCTIONS
! D. BILITZA, DECEMBER 1988
!
! INPUT:        NIGHT   LOGICAL VARIABLE FOR DAY/NIGHT DISTINCTION
!               F1REG   LOGICAL VARIABLE FOR F1 OCCURRENCE
!               XNMF2   F2 PEAK ELECTRON DENSITY [M-3]
!               XNMF1   F1 PEAK ELECTRON DENSITY [M-3]
!               XNME    E  PEAK ELECTRON DENSITY [M-3]
!               VNE     ELECTRON DENSITY AT VALLEY BASE [M-3]
!               HMF2    F2 PEAK ALTITUDE [KM]
!               HMF1    F1 PEAK ALTITUDE [KM]
!               HME     E  PEAK ALTITUDE [KM]
!               HV1     ALTITUDE OF VALLEY TOP [KM]
!               HV2     ALTITUDE OF VALLEY BASE [KM]
!               HHALF   ALTITUDE OF HALF-F2-PEAK-DENSITY [KM]
!
! OUTPUT:       HXL(4)  HEIGHT PARAMETERS FOR LAY FUNCTIONS [KM] 
!               SCL(4)  SCALE PARAMETERS FOR LAY FUNCTIONS [KM]
!               AMP(4)  AMPLITUDES FOR LAY FUNCTIONS
!               IQUAL   =0 ok, =1 ok using second choice for HXL(1)
!                       =2 NO SOLUTION
!---------------------------------------------------------------  
SUBROUTINE INILAY(NIGHT,F1REG,XNMF2,XNMF1,XNME,VNE,HMF2,HMF1, &
                  HME,HV1,HV2,HHALF,HXL,SCL,AMP,IQUAL)

    DIMENSION       XX(8),YY(8),WW(8),AMP(4),HXL(4),SCL(4)
    LOGICAL         SSIN,NIGHT,F1REG
    
    ! constants --------------------------------------------------------
    NUMLAY = 4
    NC1 = 2
    ALG102 = LOG10(2.)
    
    ! constraints: xx == height     yy == log(Ne/NmF2)    ww == weights
    ! -----------------------------------------------------------------
    ALOGF  = LOG10(XNMF2)
    ALOGEF = LOG10(XNME) - ALOGF
    XHALF = XNMF2/2.
    XX(1) = HHALF
    XX(2) = HV1
    XX(3) = HV2
    XX(4) = HME
    XX(5) = HME - ( HV2 - HME )
    YY(1) = -ALG102
    YY(2) = ALOGEF
    YY(3) = LOG10(VNE) - ALOGF
    YY(4) = ALOGEF
    YY(5) = YY(3)
    YY(7) = 0.0
    WW(2) = 1.
    WW(3) = 2.
    WW(4) = 5.
    
    ! geometric paramters for LAY -------------------------------------
    ! difference to earlier version:  HXL(3) = HV2 + SCL(3)
    SCL0   = 0.7*( 0.216*( HMF2 - HHALF ) + 56.8 )
    SCL(1) = 0.8 * SCL0
    SCL(2) = 10.
    SCL(3) = 9.
    SCL(4) = 6.
    HXL(3) = HV2
    HFFF   = HHALF
    XFFF   = XHALF
    
    IF(NIGHT) THEN
        ! NIGHT CONDITION
        ! different HXL,SCL values were tested including: 
        !       SCL(1) = HMF2 * 0.15 - 27.1     HXL(2) = 200.   
        !       HXL(2) = HMF1 + SCL(2)          HXL(3) = 140.
        !       SCL(3) = 5.                     HXL(4) = HME + SCL(4)
        !       HXL(4) = 105.                   
        NUMCON = 7
        HXL(1) = HHALF
        HXL1T  = 0.4 * HMF2 + 30.
        HXL(2) = ( HMF2 + HV1 ) / 2.
        HXL(4) = HME
        XX(6) = HV2
        XX(7) = HME
        YY(6) = 0.0
        WW(1) = 1.
        WW(3) = 3.
        WW(5) = 0.5
        WW(6) = 50.
        WW(7) = 500.
        HFFF=HHALF
        XFFF=XHALF
    ELSE
        ! DAY CONDITION
        ! earlier tested:       HXL(2) = HMF1 + SCL(2)
        NUMCON = 8
        HXL(1) = 0.9 * HMF2
        HXL1T  = HHALF
        HXL(2) = HMF1
        HXL(4) = HME - SCL(4)
        XX(6) = HMF1
        XX(7) = HV2
        XX(8) = HME
        YY(8) = 0.0
        WW(5) = 1.
        WW(7) = 50.
        WW(8) = 500.

        IF(F1REG) THEN
            ! with F-region
            YY(6) = LOG10(XNMF1) - ALOGF
            WW(6) = 3.

            IF( (XNMF1-XHALF)*(HMF1-HHALF) < 0.0) THEN
                WW(1)=0.5
            ELSE
                ZET = YY(1) - YY(6)
                WW(1) = EPST( ZET, 0.1, 0.15)
            ENDIF

            IF(HHALF > HMF1) THEN
                HFFF = HMF1
                XFFF = XNMF1
            ELSE
                HFFF = HHALF
                XFFF = XHALF
            ENDIF
        ELSE
            ! without F-region
            HXL(2) = ( HMF2 + HHALF )/2.
            YY(6)  = 0.
            WW(6)  = 0.
            WW(1)  = 1.
        END IF
    END IF
    
    ! are valley-top and bottomside point compatible ?
    IF((HV1-HFFF)*(XNME-XFFF) < 0.0) WW(2) = 0.5
    IF(HV1 <= HV2+5.0) WW(2)=0.5
    
    ! DETERMINE AMPLITUDES
    NC0=NUMCON-NC1
    IQUAL=0

    do
        CALL LSKNM(NUMLAY,NUMCON,NC0,NC1,HMF2,SCL,&
                   HXL,WW,XX,YY, AMP,SSIN)

        IF(IQUAL > 0) exit

        IF((ABS(AMP(1)) > 10.0).OR.(SSIN)) THEN
            IQUAL=1
            HXL(1) = HXL1T
        ELSE
            exit
        END IF
    enddo

    IF(SSIN) IQUAL=2

    RETURN
END SUBROUTINE INILAY
