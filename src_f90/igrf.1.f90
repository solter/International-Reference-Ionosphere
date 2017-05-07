! igrf.for, version number can be found at the end of this comment.
!-----------------------------------------------------------------------
!
! Subroutines to compute IGRF parameters for IRI and all functions and
! subroutines required for this computation, including:
!     IGRF_SUB, IGRF_DIP, FINDB0, SHELLG, STOER, FELDG, FELDCOF, GETSHC,
!     INTERSHC, EXTRASHC, GEODIP, fmodip
!
! CGM coordinates : GEOCGM01, OVL_ANG, CGMGLA, CGMGLO, DFR1DR,
!   AZM_ANG, MLTUT, MFC, FTPRNT, GEOLOW, CORGEO, GEOCOR, SHAG, RIGHT,
!   IGRF, RECALC, SPHCAR, BSPCAR, GEOMAG, MAGSM, SMGSM
!
! MLT: CLCMLT, DPMTRX
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Required i/o units:
!  KONSOL= 6 Program messages (used when jf(12)=.true. -> konsol)
!  KONSOL=11 Program messages (used when jf(12)=.false. -> MESSAGES.TXT)
!
!     COMMON/iounit/konsol,mess is used to pass the value of KONSOL from
!     IRISUB to IRIFUN and IGRF. If mess=false then messages are turned off.
!
!  UNIT=14 IGRF/GETSHC: IGRF coeff. (DGRF%%%%.DAT or IGRF%%%%.DAT, %%%%=year)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!-----------------------------------------------------------------------
! INPUT:
!    xlat      geodatic latitude in degrees
!    xlong     geodatic longitude in degrees
!    year      decimal year (year+(month-0.5)/12.0-0.5 or
!                  year+day-of-year/365 or ../366 if leap year)
!    height    height in km
! OUTPUT:
!    xl        L value
!    icode      =1  L is correct; =2  L is not correct;
!               =3  an approximation is used
!    dipl      dip latitude in degrees
!    babs      magnetic field strength in Gauss
!-----------------------------------------------------------------------
subroutine igrf_sub(xlat, xlong, year, height, xl, icode, dipl, babs)

    REAL              LATI,LONGI
    COMMON /CONST/UMR,PI

    lati  = xlat
    longi = xlong
    CALL FELDG(LATI, LONGI, HEIGHT, BNORTH, BEAST, BDOWN, BABS)
    CALL SHELLG(LATI, LONGI, HEIGHT, XL, ICODE, BAB1)
    DIPL = ATAN( BDOWN/2.0/sqrt(BNORTH**2 + BEAST**2) )/umr
    RETURN
END SUBROUTINE igrf_sub


!-----------------------------------------------------------------------
! INPUT:
!    xlat      geodatic latitude in degrees
!    xlong     geodatic longitude in degrees
!    year      decimal year (year+month/12.0-0.5 or
!                  year+day-of-year/365 or ../366 if leap year)
!    height    height in km
! OUTPUT:
!    dec       magnetic declination in degrees
!    dip       magnetic inclination (dip) in degrees
!    dipl      dip latitude in degrees
!    ymodip    modified dip latitude = asin{dip/sqrt[dip^2+cos(LATI)]}
!-----------------------------------------------------------------------
subroutine igrf_dip(xlat, xlong, year, height, dec, dip, dipl, ymodip)

    COMMON /CONST/UMR,PI

    xlati  = xlat
    xlongi = xlong
    h = height

    CALL FELDG(XLATI, XLONGI, H, BNORTH, BEAST, BDOWN, BABS)
    DECARG = BEAST/SQRT(BEAST**2 + BNORTH**2)
    IF(ABS(DECARG) > 1.) DECARG = SIGN(1.,DECARG)

    DEC = ASIN(DECARG)
    BDBA = BDOWN/BABS
    IF(ABS(BDBA) > 1.) BDBA = SIGN(1.,BDBA)

    DIP = ASIN(BDBA)
    dipdiv = DIP/SQRT( DIP**2 + cos(XLATI*UMR) )
    IF(ABS(dipdiv) > 1.) dipdiv = SIGN(1.,dipdiv)

    SMODIP = ASIN(dipdiv)
    DIPL = ATAN( BDOWN/2.0/sqrt( BNORTH**2 + BEAST**2) )/umr
    YMODIP = SMODIP/UMR
    DEC = DEC/UMR
    DIP = DIP/UMR
    RETURN
END SUBROUTINE igrf_dip


!*********************************************************************
!  SUBROUTINES SHELLG, SHELLC,  STOER, FELDG, FELDCOF, GETSHC,       *
!       INTERSHC, EXTRASHC                                           *
!*********************************************************************
!*********************************************************************


!-----------------------------------------------------------------------
!      SUBROUTINE SHELLG(GLAT,GLON,ALT,FL,ICODE,B0)
!-----------------------------------------------------------------------
! CALCULATES L-VALUE FOR SPECIFIED GEODAETIC COORDINATES, ALTITUDE
! AND GEMAGNETIC FIELD MODEL.
! REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTER, INTERNAL NOTE
!      NO. 67, 1970.
!      G. KLUGE, COMPUTER PHYSICS COMMUNICATIONS 3, 31-35, 1972
!-----------------------------------------------------------------------
!  INPUT:     GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
!             GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
!             ALT   ALTITUDE IN KM ABOVE SEA LEVEL
!
!          DIMO     DIPOL MOMENT IN GAUSS (NORMALIZED TO EARTH RADIUS)
!
!          COMMON
!             X(3)    NOT USED
!             H(144)  FIELD MODEL COEFFICIENTS ADJUSTED FOR SHELLG
!-----------------------------------------------------------------------
!  OUTPUT: FL           L-VALUE
!          ICODE        =1 NORMAL COMPLETION
!                       =2 UNPHYSICAL CONJUGATE POINT (FL MEANINGLESS)
!                       =3 SHELL PARAMETER GREATER THAN LIMIT UP TO
!                          WHICH ACCURATE CALCULATION IS REQUIRED;
!                          APPROXIMATION IS USED.
!          B0           MAGNETIC FIELD STRENGTH IN GAUSS
!-----------------------------------------------------------------------
SUBROUTINE SHELLG(GLAT, GLON, ALT, FL, ICODE, B0)

    COMMON/IGRF2/     X(3), H(196)
    COMMON/CONST/     UMR, PI
    COMMON/IGRF1/     ERA, AQUAD, BQUAD, DIMO

    DIMENSION :: V(3), U(3,3), P(8,100), SP(3)
    
    ! Convert Geodetic to Cartesian
    RLAT = GLAT*UMR
    CT = SIN(RLAT)
    ST = COS(RLAT)
    D = SQRT( AQUAD - ( AQUAD - BQUAD )*CT**2 )
    X(1) = ( ALT + AQUAD/D )*ST/ERA
    X(3) = ( ALT + BQUAD/D )*CT/ERA
    RLON = GLON*UMR
    X(2) = X(1)*SIN(RLON)
    X(1) = X(1)*COS(RLON)

    ! Call into the ShellC subroutine
    call SHELLC(X, DIMO, FL, ICODE, B0)

END SUBROUTINE SHELLG


!-----------------------------------------------------------------------
!      SUBROUTINE SHELLC(V, DIMO, FL, ICODE, B0)
!-----------------------------------------------------------------------
! CALCULATES L-VALUE FOR SPECIFIED CARTESIAN COORDINATES, ALTITUDE
! AND GEMAGNETIC FIELD MODEL.
! REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTER, INTERNAL NOTE
!      NO. 67, 1970.
!      G. KLUGE, COMPUTER PHYSICS COMMUNICATIONS 3, 31-35, 1972
!-----------------------------------------------------------------------
!  INPUT:       V(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
!                     X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
!                     Y-AXIS POINTING TO EQUATOR AT 90 LONG.
!                     Z-AXIS POINTING TO NORTH POLE
!               DIMO  Dipole Moment (prevents a common block need.
!
!          COMMON
!             X(3)    NOT USED
!             H(144)  FIELD MODEL COEFFICIENTS ADJUSTED FOR SHELLG
!-----------------------------------------------------------------------
!  OUTPUT: FL           L-VALUE
!          ICODE        =1 NORMAL COMPLETION
!                       =2 UNPHYSICAL CONJUGATE POINT (FL MEANINGLESS)
!                       =3 SHELL PARAMETER GREATER THAN LIMIT UP TO
!                          WHICH ACCURATE CALCULATION IS REQUIRED;
!                          APPROXIMATION IS USED.
!          B0           MAGNETIC FIELD STRENGTH IN GAUSS
!-----------------------------------------------------------------------
SUBROUTINE SHELLC(V, DIMO, FL, ICODE, B0)

    logical           finished
    DIMENSION         V(3), U(3, 3), P(8, 100), SP(3), ZZ(7)
    COMMON/IGRF2/     X(3), H(196)
    COMMON/FIDB0/     SP    

    !-- RMIN, RMAX ARE BOUNDARIES FOR IDENTIFICATION OF ICODE=2 AND 3
    !-- STEP IS STEP SIZE FOR FIELD LINE TRACING
    !-- STEQ IS STEP SIZE FOR INTEGRATION
    
    DATA RMIN,RMAX    /0.05,1.01/
    DATA STEP,STEQ    /0.20,0.03/
    DATA U/+0.3511737,-0.9148385,-0.1993679,&
           +0.9335804,+0.3583680,+0.0000000,&
           +0.0714471,-0.1861260,+0.9799247/

    BEQU = 1.E10

    !*****CONVERT TO DIPOL-ORIENTED CO-ORDINATES
    X = V
    RQ = 1./SUM(X**2)
    R3H = SQRT( RQ*SQRT(RQ) )

    P(1:3,2) = MATMUL(X, U) * (/ R3H, R3H, RQ /)

    !*****FIRST THREE POINTS OF FIELD LINE
    STEP = -SIGN(STEP, P(3,2))
    CALL STOER(P(1:7,2), BQ2, R2)

    B0 = SQRT(BQ2)
    P(1,3) = P(1,2) + 0.5*STEP*P(4,2)
    P(2,3) = P(2,2) + 0.5*STEP*P(5,2)
    P(3,3) = P(3,2) + 0.5*STEP
    CALL STOER(P(1:7,3), BQ3, R3)

    P(1,1) = P(1,2) - STEP*( 2.*P(4,2) - P(4,3) )
    P(2,1) = P(2,2) - STEP*( 2.*P(5,2) - P(5,3) )
    P(3,1) = P(3,2) - STEP
    CALL STOER(P(1:7,1), BQ1, R1)

    P(1,3) = P(1,2) + STEP*( 20.*P(4,3) - 3.*P(4,2) + P(4,1) )/18.
    P(2,3) = P(2,2) + STEP*( 20.*P(5,3) - 3.*P(5,2) + P(5,1) )/18.
    P(3,3) = P(3,2) + STEP
    CALL STOER(P(1:7,3), BQ3, R3)

    !*****INVERT SENSE IF REQUIRED
    IF(BQ3 > BQ1) THEN
        STEP = -STEP
        R3 = R1
        BQ3 = BQ1
        ZZ = P(1:7,1)
        P(1:7,1) = P(1:7,3)
        P(1:7,3) = ZZ
    end if

    !*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
    IF(BQ1 < BEQU) THEN
        BEQU = BQ1
        IEQU = 1
    ENDIF
    IF(BQ2 < BEQU) THEN
        BEQU = BQ2
        IEQU = 2
    ENDIF
    IF(BQ3 < BEQU) THEN
        BEQU = BQ3
        IEQU = 3
    ENDIF

    !*****INITIALIZATION OF INTEGRATION LOOPS
    STEP12 = STEP/12.
    STEP2 = 2*STEP
    STEQ = SIGN(STEQ, STEP)
    FI = 0.
    ICODE = 1
    ORADIK = 0.
    OTERM = 0.
    STP = R2*STEQ
    Z = P(3,2) + STP
    STP = STP/0.75
    P(8,1) = STEP2*( P(1,1)*P(4,1) + P(2,1)*P(5,1) )
    P(8,2) = STEP2*( P(1,2)*P(4,2) + P(2,2)*P(5,2) )

    !*****MAIN LOOP (FIELD LINE TRACING)
    finished = .FALSE.
    DO N=3,99
        !*****CORRECTOR (FIELD LINE TRACING)
        P(1,N) = P(1,N-1) + STEP12*( 5.*P(4,N) + 8.*P(4,N-1) - P(4,N-2) )
        P(2,N) = P(2,N-1) + STEP12*( 5.*P(5,N) + 8.*P(5,N-1) - P(5,N-2))

        !*****PREPARE EXPANSION COEFFICIENTS FOR INTERPOLATION
        !*****OF SLOWLY VARYING QUANTITIES
        P(8,N) = STEP2*( P(1,N)*P(4,N) + P(2,N)*P(5,N) )
        C0 =   P(1,N-1)**2 + P(2,N-1)**2
        C1 =   P(8,N-1)
        C2 = ( P(8,N) - P(8,N-2) )*0.25
        C3 = ( P(8,N) + P(8,N-2) )/6.0 - C1/3.
        D0 = P(6,N-1)
        D1 = ( P(6,N) - P(6,N-2) )*0.5
        D2 = ( P(6,N) + P(6,N-2) )*0.5 - D0
        E0 = P(7,N-1)
        E1 = ( P(7,N) - P(7,N-2) )*0.5
        E2 = ( P(7,N) + P(7,N-2) )*0.5 - E0

        !*****INNER LOOP (FOR QUADRATURE)
        do
            T = ( Z - P(3,N-1) )/STEP
            ! Terminate the loop
            IF(T > 1.) EXIT

            HLI = C3*T + C2
            HLI = T*HLI + C1
            HLI = T*HLI + C0
            HLI = .5*HLI
            ZQ = Z**2
            R = HLI + SQRT( HLI**2 + ZQ )

            IF(R <= RMIN) then
                ! APPROXIMATION FOR HIGH VALUES OF L.
                ICODE=3
                T = -P(3,N-1)/STEP
                FL = 1./(ABS(((C3*T+C2)*T+C1)*T+C0)+1E-15)
                FL = C3*T + C2
                FL = T*FL + C1
                FL = FL*T + C0
                FL = ABS(FL) + 1E-15
                FL = 1./FL
                RETURN
            endif

            RQ = R**2
            FF = SQRT( 1. + 3.*ZQ/RQ )
            RADIK = B0 - ( (D2*T + D1)*T + D0 )*R*RQ*FF
            IF(R > RMAX) then
                ICODE=2
                RADIK = RADIK - 12.*(R - RMAX)**2
            endif
            IF(2*RADIK <= ORADIK) then
                finished = .TRUE.
                EXIT
            endif

            TERM = SQRT(RADIK)*FF*( (E2*T+E1)*T + E0 )/(RQ + ZQ)
            FI = FI+STP*( OTERM + TERM )
            ORADIK = RADIK
            OTERM = TERM
            STP = R*STEQ
            Z = Z + STP
        enddo
        if (finished) EXIT

        !*****PREDICTOR (FIELD LINE TRACING)
        P(1,N+1) = P(1,N) + STEP12*( 23.*P(4,N) - 16.*P(4,N-1) + 5.*P(4,N-2) )
        P(2,N+1) = P(2,N) + STEP12*( 23.*P(5,N) - 16.*P(5,N-1) + 5.*P(5,N-2) )
        P(3,N+1) = P(3,N) + STEP
        CALL STOER(P(1:7,N+1), BQ3, R3)

        !*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
        IF(BQ3 < BEQU) THEN
            IEQU = N + 1
            BEQU = BQ3
        ENDIF
    enddo

    IF(IEQU < 2) IEQU = 2
    SP(1:3) = P(1:3, IEQU-1)
    IF(ORADIK >= 1E-15) FI = FI + STP/0.75*OTERM*ORADIK/( ORADIK - RADIK )
    
    !-- The minimal allowable value of FI was changed from 1E-15 to 1E-12,
    !-- because 1E-38 is the minimal allowable arg. for ALOG in our envir.
    !-- D. Bilitza, Nov 87.
    FI = 0.5*ABS(FI)/SQRT(B0) + 1E-12
    
    !*****COMPUTE L FROM B AND I.  SAME AS CARMEL IN INVAR.
    !-- Correct dipole moment is used here. D. Bilitza, Nov 87.
    DIMOB0 = DIMO/B0
    arg1 = LOG(FI)
    arg2 = LOG(DIMOB0)
    XX = 3*arg1 - arg2

    if (XX <= -22.) then
        GG = 3.33338E-1
        GG = GG*XX + 3.0062102E-1
    else if (XX <= -3.0) then
        GG = -8.1537735E-14
        GG = GG*XX + 8.3232531E-13
        GG = GG*XX + 1.0066362E-9
        GG = GG*XX + 8.1048663E-8
        GG = GG*XX + 3.2916354E-6
        GG = GG*XX + 8.2711096E-5
        GG = GG*XX + 1.3714667E-3
        GG = GG*XX + 1.5017245E-2
        GG = GG*XX + 4.3432642E-1
        GG = GG*XX + 6.2337691E-1
    else if( XX <= 3.0) then
        GG = 2.6047023E-10
        GG = GG*XX + 2.3028767E-9
        GG = GG*XX - 2.1997983E-8
        GG = GG*XX - 5.3977642E-7
        GG = GG*XX - 3.3408822E-6
        GG = GG*XX + 3.8379917E-5
        GG = GG*XX + 1.1784234E-3
        GG = GG*XX + 1.4492441E-2
        GG = GG*XX + 4.3352788E-1
        GG = GG*XX + 6.228644E-1
    else if( XX <= 11.7 ) then
        GG = 6.3271665E-10
        GG = GG*XX - 3.958306E-8
        GG = GG*XX + 9.9766148E-7
        GG = GG*XX - 1.2531932E-5
        GG = GG*XX + 7.9451313E-5
        GG = GG*XX - 3.2077032E-4
        GG = GG*XX + 2.1680398E-3
        GG = GG*XX + 1.2817956E-2
        GG = GG*XX + 4.3510529E-1
        GG = GG*XX + 6.222355E-1
    else if ( XX < 23.0 ) then
        GG = 2.8212095E-8
        GG = GG*XX - 3.8049276E-6
        GG = GG*XX + 2.170224E-4
        GG = GG*XX - 6.7310339E-3
        GG = GG*XX + 1.2038224E-1
        GG = GG*XX - 1.8461796E-1
        GG = GG*XX + 2.0007187E0
    else
        GG = XX - 3.0460681E0
    endif

    FL = EXP( LOG( (1. + EXP(GG))*DIMOB0 )/3.0 )
    RETURN
END SUBROUTINE SHELLC


!*******************************************************************
!* SUBROUTINE USED FOR FIELD LINE TRACING IN SHELLG                *
!                                                                  *
!*******************************************************************
SUBROUTINE STOER(P,BQ,R)
    DIMENSION         P(7), U(3,3), DX(3), DM(3), XM(3)
    COMMON/IGRF2/     XI(3), H(196)

    !*****XM,YM,ZM  ARE GEOMAGNETIC CARTESIAN INVERSE CO-ORDINATES
    FLI= P(1)**2 + P(2)**2 + 1E-15
    R  = 0.5*( FLI + SQRT( FLI**2 + 4*ZM**2 ) )
    RQ = R**2
    WR = SQRT(R)

    XM(1) = P(1)*WR
    XM(2) = P(2)*WR
    XM(3) = P(3)

    !*****TRANSFORM TO GEOGRAPHIC CO-ORDINATE SYSTEM
    DATA U/+0.3511737, -0.9148385, -0.1993679,&
           +0.9335804, +0.3583680, +0.0000000,&
           +0.0714471, -0.1861260, +0.9799247/
    XI = matmul(U, XM)

    !*****COMPUTE DERIVATIVES
    !      CALL FELDI(XI,H)
    CALL FELDI
    Q = H(1)/RQ
    DX(1) = 2*H(3) + Q*XI(1)
    DX(2) = 2*H(4) + Q*XI(2)
    DX(3) = 2*H(2) + Q*XI(3)

    !*****TRANSFORM BACK TO GEOMAGNETIC CO-ORDINATE SYSTEM
    DM = matmul(DX, U)
    DR  = dot_product(XM, DM)/R

    !*****FORM SLOWLY VARYING EXPRESSIONS
    P(4) = ( WR*DM(1) - 0.5*P(1)*DR )/(R*DM(3))
    P(5) = ( WR*DM(2) - 0.5*P(2)*DR )/(R*DM(3))
    DSQ  = RQ*dot_product(DM, DM)
    BQ   = DSQ*RQ**2
    P(6) = SQRT( DSQ/( RQ + 3.*XM(3)**2 ) )
    P(7) = P(6)*( RQ + XM(3)**2 )/(RQ*DM(3))
    RETURN
END SUBROUTINE STOER


!-----------------------------------------------------------------------
! CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL
! REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61,
!      1970.
!-----------------------------------------------------------------------
!  INPUT:  
!         GLAT  GEODETIC LATITUDE IN DEGREES (NORTH)
!         GLON  GEODETIC LONGITUDE IN DEGREES (EAST)
!         ALT   ALTITUDE IN KM ABOVE SEA LEVEL
!
!
!          COMMON /IGRF1/ and /CONST/
!               UMR     = ATAN(1.0)*4./180.   <DEGREE>*UMR=<RADIANT>
!               ERA     EARTH RADIUS FOR NORMALIZATION OF CARTESIAN
!                       COORDINATES (6371.2 KM)
!               AQUAD, BQUAD   SQUARE OF MAJOR AND MINOR HALF AXIS OF
!                       EARTH ELLIPSOID AS RECOMMENDED BY INTERNAT.
!                       ASTRONOMICAL UNION (6378.160, 6356.775 KM).
!-----------------------------------------------------------------------
!  OUTPUT: 
!          BNORTH, BEAST, BDOWN   COMPONENTS OF THE FIELD WITH RESPECT
!                 TO THE LOCAL GEODETIC COORDINATE SYSTEM, WITH AXIS
!                 POINTING IN THE TANGENTIAL PLANE TO THE NORTH, EAST
!                 AND DOWNWARD.
!
!           BABS   MAGNETIC FIELD STRENGTH IN GAUSS
!-----------------------------------------------------------------------
SUBROUTINE FELDG(GLAT, GLON, ALT, BNORTH, BEAST, BDOWN, BABS)
    DIMENSION         X(3),B(3)

    COMMON/IGRF1/ERA,AQUAD,BQUAD,DIMO    /CONST/UMR,PI

    RLAT = GLAT*UMR
    CT = SIN(RLAT)
    ST = COS(RLAT)
    D  = SQRT( AQUAD - (AQUAD-BQUAD)*CT**2 )

    RLON = GLON*UMR
    CP = COS(RLON)
    SP = SIN(RLON)

    RHO = ( ALT + AQUAD/D )*ST/ERA

    X(1) = RHO*CP
    X(2) = RHO*SP
    X(3) = ( ALT + BQUAD/D )*CT/ERA

    call FELDC(X, B)

    BABS  = SQRT(dot_product(B,B))
    BEAST =  B(2)*CP - B(1)*SP
    BRHO  =  B(2)*SP + B(1)*CP
    BNORTH=  B(3)*ST - BRHO*CT
    BDOWN = -B(3)*CT - BRHO*ST

    return
    
END SUBROUTINE FELDG


!-----------------------------------------------------------------------
! CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL
! REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61,
!      1970.
!-----------------------------------------------------------------------
!  INPUT:  
!         X(3)  CARTESIAN COORDINATES IN EARTH RADII (6371.2 KM)
!                   X-AXIS POINTING TO EQUATOR AT 0 LONGITUDE
!                   Y-AXIS POINTING TO EQUATOR AT 90 LONG.
!                   Z-AXIS POINTING TO NORTH POLE
!
!          COMMON /IGRF2/
!-----------------------------------------------------------------------
!  OUTPUT: 
!          B Components of the field with respect to the the 
!            cartesian coordinate system.
!-----------------------------------------------------------------------
SUBROUTINE FELDC(X, B)
    DIMENSION         X(3),B(3)
    CHARACTER*13      NAME

    COMMON/IGRF2/XI(3),H(196)

    RQ = 1./dot_product(X, X)
    XI = X*RQ

    call FELDI

    S = .5*H(1) + 2.*( H(2)*XI(3) + H(3)*XI(1) + H(4)*XI(2) )
    T = 2*RQ*SQRT(RQ)
    B(1) = T*( H(3) - S*X(1) )
    B(2) = T*( H(4) - S*X(2) )
    B(3) = T*( H(2) - S*X(3) )

    return
END SUBROUTINE FELDC


!-----------------------------------------------------------------------
! CALCULATES EARTH MAGNETIC FIELD FROM SPHERICAL HARMONICS MODEL
! REF: G. KLUGE, EUROPEAN SPACE OPERATIONS CENTRE, INTERNAL NOTE 61,
!      1970.
!-----------------------------------------------------------------------
!  INPUT:  
!
!          COMMON /MODEL/ AND /IGRF2/
!               G(M)    NORMALIZED FIELD COEFFICIENTS (SEE FELDCOF)
!                       M=NMAX*(NMAX+2)
!-----------------------------------------------------------------------
!  Results stored in the common block /IGRF2/ H array.
!-----------------------------------------------------------------------
SUBROUTINE FELDI
    CHARACTER*13      NAME
    COMMON/IGRF2/XI(3),H(196)
    COMMON/MODEL/NMAX,TIME,G(196),NAME

    IHMAX = NMAX**2 + 1
    LAST  = IHMAX + 2*NMAX
    IMAX  = 2*NMAX - 1
    H(IHMAX:LAST) = G(IHMAX:LAST)

    DO K=1,3,2
        I  = IMAX
        IH = IHMAX
        DO
            IL = IH-I
            F = 2./FLOAT(I-K+2)
            X = XI(1)*F
            Y = XI(2)*F
            Z = XI(3)*2.*F
            I = I-2

            IF (I > 1) then
                DO M=3,I,2
                    H(IL+M+1) = G(IL+M+1) &
                                + X*( H(IH+M+3) - H(IH+M-1) ) &
                                - Y*( H(IH+M+2) + H(IH+M-2) ) &
                                + Z*H(IH+M+1)
                    H(IL+M)   = G(IL+M) &
                                + X*( H(IH+M+2) - H(IH+M-2) ) &
                                + Y*( H(IH+M+3) + H(IH+M-1) ) &
                                + Z*H(IH+M)
                END DO
            END IF

            IF (I == 1) then
                H(IL+2) = G(IL+2) &
                          + X*H(IH+4) &
                          - Y*( H(IH+3) + H(IH) ) &
                          + Z*H(IH+2)
                H(IL+1) = G(IL+1) &
                          + X*( H(IH+3) - H(IH) ) &
                          + Y*H(IH+4) &
                          + Z*H(IH+1)
            END IF

            H(IL) = G(IL) &
                    + 2.*X*H(IH+1) &
                    + 2.*Y*H(IH+2) &
                    +    Z*H(IH)

            IH = IL
            IF(I < K) exit
        ENDDO
    END DO
    RETURN
END SUBROUTINE FELDI


!-----------------------------------------------------------------------
!  DETERMINES COEFFICIENTS AND DIPOL MOMENT FROM IGRF MODELS
!
!       INPUT:  YEAR    DECIMAL YEAR FOR WHICH GEOMAGNETIC FIELD IS TO
!                       BE CALCULATED
!                        COMMON/IGRF1/ERAD,AQUAD,BQUAD,DIMO /CONST/UMR,PI
!       OUTPUT:         COMMON/MODEL/NMAX,TIME,GH1,FIL1
!                        COMMON/DIPOL/GHI1,GHI2,GHI3
!
! THE GEOMAGNETIC DIPOL MOMENT (DIMO) IN GAUSS (NORMALIZED TO EARTH'S
! RADIUS) AT THE TIME (YEAR) IS COMPUTED BUT NOT USED.
!
! 05/31/2000 updated to IGRF-2000 version (###)
! 03/24/2000 updated to IGRF-2005 version (###)
! 07/22/2009 NMAX=13 for DGRF00 and IGRF05; H/G-arrays(195) ! 02/26/2010 update to IGRF-11 (2010) (###)
! 10/05/2011 added COMMON/DIPOL/ for MLT computation in DPMTRX (IRIFUN)
! 02/10/2015 update to IGRF-12 (2015) (###)
!-----------------------------------------------------------------------
SUBROUTINE FELDCOF(YEAR)
    CHARACTER*13    FILMOD, FIL1, FIL2
    ! ### FILMOD, DTEMOD array-size is number of IGRF maps
    DIMENSION       GH1(196),GH2(196),GHA(196),FILMOD(16)
    DIMENSION        DTEMOD(16)
    DOUBLE PRECISION X,F0,F

    COMMON/MODEL/   NMAX,TIME,GH1,FIL1
    COMMON/IGRF1/   ERAD,AQUAD,BQUAD,DIMO /CONST/UMR,PI
    COMMON/DIPOL/    GHI1,GHI2,GHI3

    ! ### updated coefficient file names and corresponding years
    DATA  FILMOD   / 'dgrf1945.dat','dgrf1950.dat','dgrf1955.dat',&
        'dgrf1960.dat','dgrf1965.dat','dgrf1970.dat','dgrf1975.dat',&
        'dgrf1980.dat','dgrf1985.dat','dgrf1990.dat','dgrf1995.dat',&
        'dgrf2000.dat','dgrf2005.dat','dgrf2010.dat','igrf2015.dat',&
        'igrf2015s.dat'/
    DATA  DTEMOD / 1945., 1950., 1955., 1960., 1965.,&
        1970., 1975., 1980., 1985., 1990., 1995., 2000.,2005.,&
        2010., 2015., 2020./
    
    ! ### numye is number of IGRF coefficient files minus 1
    NUMYE=15
    
    !  IS=0 FOR SCHMIDT NORMALIZATION   IS=1 GAUSS NORMALIZATION
    !  IU  IS INPUT UNIT NUMBER FOR IGRF COEFFICIENT SETS
    IU = 14
    IS = 0

    !-- DETERMINE IGRF-YEARS FOR INPUT-YEAR
    TIME = YEAR
    IYEA = INT(YEAR/5.)*5
    L = (IYEA - 1945)/5 + 1
    IF(L < 1) L = 1
    IF(L > NUMYE) L = NUMYE
    DTE1 = DTEMOD(L)
    FIL1 = FILMOD(L)
    DTE2 = DTEMOD(L+1)
    FIL2 = FILMOD(L+1)

    !-- GET IGRF COEFFICIENTS FOR THE BOUNDARY YEARS
    CALL GETSHC (IU, FIL1, NMAX1, ERAD, GH1, IER)
    IF (IER  /=  0) STOP !TODO - do something more intelligent...
    CALL GETSHC (IU, FIL2, NMAX2, ERAD, GH2, IER)
    IF (IER  /=  0) STOP !TODO - do something more intelligent...

    !-- DETERMINE IGRF COEFFICIENTS FOR YEAR
    IF (L <= NUMYE-1) THEN
        CALL INTERSHC (YEAR, DTE1, NMAX1, GH1, DTE2, NMAX2, GH2, NMAX, GHA)
    ELSE
        CALL EXTRASHC (YEAR, DTE1, NMAX1, GH1, NMAX2, GH2, NMAX, GHA)
    ENDIF

    !-- DETERMINE MAGNETIC DIPOL MOMENT AND COEFFIECIENTS G
    DIMO = sqrt(dot_product(GHA(1:3), GHA(1:3))) * 1.D-5

    GHI1 = GHA(1)
    GHI2 = GHA(2)
    GHI3 = GHA(3)
    GH1(1) = 0.0
    I = 2
    F0 = 1.D-5
    IF(IS == 0) F0 = -F0

    DO N=1,NMAX
        X = N

        F0 = F0*X**2/( 4.D0*X - 2.D0 )
        IF(IS == 0) F0 = F0*( 2.D0*X - 1.D0 )/X
        F = F0*0.5D0
        IF(IS == 0) F = F*SQRT(2.)

        GH1(I) = GHA(I-1)*F0

        I = I + 1
        DO M=1,N
            F = F*( X + M )/( X - M + 1.D0 )
            IF(IS == 0) F = F*SQRT( (X-M+1.D0)/(X+M) )

            GH1(I)   = GHA(I-1)*F
            GH1(I+1) = GHA(I)*F

            I = I + 2
        END DO
    END DO
    RETURN
END SUBROUTINE FELDCOF


! ===============================================================
!       Reads spherical harmonic coefficients from the specified
!       file into an array.
!       Input:
!           IU    - Logical unit number
!           FSPEC - File specification
!       Output:
!           NMAX  - Maximum degree and order of model
!           ERAD  - Earth's radius associated with the spherical
!                   harmonic coefficients, in the same units as
!                   elevation
!           GH    - Schmidt quasi-normal internal spherical
!                   harmonic coefficients
!           IER   - Error number: =  0, no error
!                                 = -2, records out of order
!                                 = FORTRAN run-time error number
! ===============================================================
SUBROUTINE GETSHC (IU, FSPEC, NMAX, ERAD, GH, IER)
    CHARACTER  FSPEC*(*), FOUT*80
    DIMENSION       GH(196)
    LOGICAL        mess

    COMMON/iounit/konsol,mess

    GH = 0.0

    ! ---------------------------------------------------------------
    !       Open coefficient file. Read past first header record.
    !       Read degree and order of model and Earth's radius.
    ! ---------------------------------------------------------------
    WRITE(FOUT,'(A13)') FSPEC

    OPEN (IU, FILE=FOUT, STATUS='OLD', IOSTAT=IER, ERR=999)
    READ (IU, *, IOSTAT=IER, ERR=999)
    READ (IU, *, IOSTAT=IER, ERR=999) NMAX, ERAD, XMYEAR
    nm = nmax*(nmax+2)
    READ (IU, *, IOSTAT=IER, ERR=999) (GH(i),i=1,nm)
    CLOSE (IU)

    RETURN

    ! File IO Error
999 if (mess) write(konsol,"('Error while reading ',A13)") FOUT
    RETURN
END SUBROUTINE GETSHC


! ===============================================================
!
!       Version 1.01
!
!       Interpolates linearly, in time, between two spherical
!       harmonic models.
!
!       Input:
!           DATE  - Date of resulting model (in decimal year)
!           DTE1  - Date of earlier model
!           NMAX1 - Maximum degree and order of earlier model
!           GH1   - Schmidt quasi-normal internal spherical
!                   harmonic coefficients of earlier model
!           DTE2  - Date of later model
!           NMAX2 - Maximum degree and order of later model
!           GH2   - Schmidt quasi-normal internal spherical
!                   harmonic coefficients of later model
!
!       Output:
!           GH    - Coefficients of resulting model
!           NMAX  - Maximum degree and order of resulting model
!
!       A. Zunde
!       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225
!
! ===============================================================
SUBROUTINE INTERSHC (DATE, DTE1, NMAX1, GH1, DTE2, NMAX2, GH2, NMAX, GH)
    DIMENSION       GH1(*), GH2(*), GH(*)
    ! ---------------------------------------------------------------
    !       The coefficients (GH) of the resulting model, at date
    !       DATE, are computed by linearly interpolating between the
    !       coefficients of the earlier model (GH1), at date DTE1,
    !       and those of the later model (GH2), at date DTE2. If one
    !       model is smaller than the other, the interpolation is
    !       performed with the missing coefficients assumed to be 0.
    ! ---------------------------------------------------------------
    FACTOR = (DATE - DTE1) / (DTE2 - DTE1)

    IF (NMAX1  ==  NMAX2) THEN
        K = NMAX1 * (NMAX1 + 2)
        NMAX = NMAX1
    ELSE IF (NMAX1  >  NMAX2) THEN
        K = NMAX2 * (NMAX2 + 2)
        L = NMAX1 * (NMAX1 + 2)
        GH(K+1:L) = GH1(K+1:L) + FACTOR * (-GH1(K+1:L))
        NMAX = NMAX1
    ELSE
        K = NMAX1 * (NMAX1 + 2)
        L = NMAX2 * (NMAX2 + 2)
        GH(K+1:L) = FACTOR * GH2(K+1:L)
        NMAX = NMAX2
    ENDIF

    GH(1:K) = GH1(1:K) + FACTOR * (GH2(1:K) - GH1(1:K))
    RETURN
END SUBROUTINE INTERSHC


! ===============================================================
!
!       Version 1.01
!
!       Extrapolates linearly a spherical harmonic model with a
!       rate-of-change model.
!
!       Input:
!           DATE  - Date of resulting model (in decimal year)
!           DTE1  - Date of base model
!           NMAX1 - Maximum degree and order of base model
!           GH1   - Schmidt quasi-normal internal spherical
!                   harmonic coefficients of base model
!           NMAX2 - Maximum degree and order of rate-of-change
!                   model
!           GH2   - Schmidt quasi-normal internal spherical
!                   harmonic coefficients of rate-of-change model
!
!       Output:
!           GH    - Coefficients of resulting model
!           NMAX  - Maximum degree and order of resulting model
!
!       A. Zunde
!       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225
!
! ===============================================================
SUBROUTINE EXTRASHC (DATE, DTE1, NMAX1, GH1, NMAX2, GH2, NMAX, GH)
    DIMENSION       GH1(*), GH2(*), GH(*)
    ! ---------------------------------------------------------------
    !       The coefficients (GH) of the resulting model, at date
    !       DATE, are computed by linearly extrapolating the coef-
    !       ficients of the base model (GH1), at date DTE1, using
    !       those of the rate-of-change model (GH2), at date DTE2. If
    !       one model is smaller than the other, the extrapolation is
    !       performed with the missing coefficients assumed to be 0.
    ! ---------------------------------------------------------------
    FACTOR = (DATE - DTE1)

    IF (NMAX1  ==  NMAX2) THEN
        K = NMAX1 * (NMAX1 + 2)
        NMAX = NMAX1
    ELSE IF (NMAX1  >  NMAX2) THEN
        K = NMAX2 * (NMAX2 + 2)
        L = NMAX1 * (NMAX1 + 2)
        GH(K+1:L) = GH1(K+1:L)
        NMAX = NMAX1
    ELSE
        K = NMAX1 * (NMAX1 + 2)
        L = NMAX2 * (NMAX2 + 2)
        GH(K+1:L) = FACTOR * GH2(K+1:L)
        NMAX = NMAX2
    ENDIF

    GH(1:K) = GH1(1:K) + FACTOR * GH2(1:K)
    RETURN
END SUBROUTINE EXTRASHC


! ===============================================================
!  Calculates dipole geomagnetic coordinates from geocentric coordinates
!  or vice versa.
!                     J=0           J=1
!        INPUT:     J,SLA,SLO     J,DLA,DLO
!        OUTPUT:     DLA,DLO       SLA,SLO
!  Last revision: November 2005 (Vladimir Papitashvili)
!  The code is modifed from GEOCOR written by V.Popov and V.Papitashvili
!  in mid-1980s.
! ===============================================================
SUBROUTINE GEODIP(IYR,SLA,SLO,DLA,DLO,J)
    COMMON /CONST/UMR,PI
    !  Earth's radius (km) RE = 6371.2
    !  The radius of the sphere to compute the coordinates (in Re)
    !        RH = (RE + HI)/RE
    R = 1.
    if(j <= 0) then
        COL = (90.- SLA)*UMR
        RLO = SLO*UMR
        CALL SPHCAR(R, COL, RLO, X, Y, Z, 1)
        CALL GEOMAG(X, Y, Z, XM, YM, ZM, 1, IYR)
        CALL SPHCAR(RM, TH, PF, XM, YM, ZM, -1)
        SZM = ZM
        DLO = PF/UMR
        DCO = TH/UMR
        DLA = 90.- DCO
    ELSE
        COL = (90.- DLA)*UMR
        RLO = DLO*UMR
        CALL SPHCAR(R, COL, RLO, XM, YM, ZM, 1)
        CALL GEOMAG(X, Y, Z, XM, YM, ZM, -1, IYR)
        CALL SPHCAR(RM, TH, PF, X, Y, Z, -1)
        SZM = ZM
        SLO = PF/UMR
        SCO = TH/UMR
        SLA = 90.- SCO
    END IF
    RETURN
END SUBROUTINE GEODIP


FUNCTION fmodip(xlat)
    common/findRLAT/xlong,year

    call igrf_dip(xlat,xlong,year,300.,dec,dip,dipl,ymodip)
    fmodip=ymodip
    return
END FUNCTION fmodip


!  =====================================================================
!  Input parameters:
!     ICOR = +1    geo to cgm
!            -1    cgm to geo
!     IYEAR= year
!     HI   = altitude in km
!  Input/Output parameters:
!     DAT(1,i)=slar geocentric latitude (input/output if icor=+1/-1)
!     DAT(2,i)=slor geocentric longitude (input/output if icor=+1/-1)
!     DAT(3,i)=clar CGM latitude (input/output if icor=-1/+1)
!     DAT(4,i)=clor CGM longitude (input/output if icor=-1/+1)
!  Output parameters:
!     DAT(5,i)=rbm apex of the magnetic field line in Re (Re=6371.2 km)
!            (this parameter approximately equals the McIlwain L-value)
!     DAT(6,i)=btr IGRF Magnetic field H (nT)
!     DAT(7,i)=brr IGRF Magnetic field D (deg)
!     DAT(8,i)=ovl oval_angle as the azimuth to "magnetic north":
!                + east in Northern Hemisphere
!                + west in Southern Hemisphere
!     DAT(9,i)=azm meridian_angle as the azimuth to the CGM pole:
!                + east in Northern Hemisphere
!                + west in Southern Hemisphere
!     DAT(10,i)=utm magnetic local time (MLT) midnight in UT hours
!              i=1    for the start point
!              i=2    for the conjugate point of the start point (slac, sloc)
!             i=3    for the footprint at 1-Re of the start point (slaf,slof)
!             i=4    for the conjugate footprint at 1-Re of the start point
!     PLA(1)    geocentric latitude of the CGM pole in the Northern hemisphere
!     PLO(1)    geocentric longitude of the CGM pole in the Northern hemisphere
!     PLA(2)    geocentric latitude of the CGM pole in the Southern hemisphere
!     PLO(2)    geocentric longitude of the CGM pole in the Southern hemisphere
!     PLA(3)    geoce lati CGM North pole at the Earth's surface 1-Re or zero alt.
!     PLO(3)    geoce long CGM North pole at the Earth's surface 1-Re or zero alt.
!     PLA(4)    geoce lati CGM South pole at the Earth's surface 1-Re or zero alt.
!     PLO(4)    geoce long CGM South pole at the Earth's surface 1-Re or zero alt.
!
SUBROUTINE GEOCGM01(ICOR,IYEAR,HI,DAT,PLA,PLO)
    DIMENSION DAT(11,4),PLA(4),PLO(4)
    CHARACTER STR*12
    ! dla  = dipole latitude
    ! dlo  = dipole longitude

    COMMON /IYR/ IYR
    COMMON /NM/ NM

    !  Year (for example, as for Epoch 1995.0 - no fraction of the year)
    IYR = iyear
    !  Earth's radius (km)
    RE = 6371.2
    !  NM is the number of harmonics
    NM = 10
    !  The radius of the sphere to compute the coordinates (in Re)
    RH = (RE + HI)/RE

    !  Correction of latitudes and longitudes if they are entered beyond of
    !  the limits (this actually does not affect coordinate calculations
    !  but the oval/meridian angles and MLT midnight cannot be computed)
    IF (DAT(1,1) >  90.) DAT(1,1) =  180. - DAT(1,1)
    IF (DAT(1,1) < -90.) DAT(1,1) = -180. - DAT(1,1)
    IF (DAT(3,1) >  90.) DAT(3,1) =  180. - DAT(3,1)
    IF (DAT(3,1) < -90.) DAT(3,1) = -180. - DAT(3,1)
    IF (DAT(2,1) >  360.) DAT(2,1) = DAT(2,1) - 360.
    IF (DAT(2,1) < -360.) DAT(2,1) = DAT(2,1) + 360.
    IF (DAT(4,1) >  360.) DAT(4,1) = DAT(4,1) - 360.
    IF (DAT(4,1) < -360.) DAT(4,1) = DAT(4,1) + 360.

    !  Computation of CGM coordinates from geocentric ones at high- and
    !  middle latitudes
    IF (ICOR ==  1) THEN
        SLAR = DAT(1,1)
        SLOR = DAT(2,1)
        IF (ABS(SLAR) == 90.) SLOR = 360.
        CALL GEOCOR(SLAR,SLOR,RH,DLA,DLO,CLAR,CLOR,PMR)
        DAT(3,1) = CLAR
        DAT(4,1) = CLOR
    ELSE
        !  Computation of geocentric coordinates from CGM ones at high- and
        !  middle latitudes
        CLAR = DAT(3,1)
        CLOR = DAT(4,1)
        IF (ABS(CLAR) == 90.) CLOR = 360.
        CALL CORGEO(SLAR,SLOR,RH,DLA,DLO,CLAR,CLOR,PMR)
        DAT(1,1) = SLAR
        DAT(2,1) = SLOR
    ENDIF

    !  PMI is L-shell parameter for the magnetic field line; limit to 16 Re
    IF(PMR >= 16.) PMR = 999.99
    DAT(5,1) = PMR

    !  Check if CGM_Lat has been calculated, then go for the conjugate point
    IF(CLAR > 999.) THEN
        !  CGM_Lat has NOT been calculated, call GEOLOW for computation of the
        !  CGM coordinates at low latitudes using the CBM approach (see the
        !  reference in GEOLOW)
        CALL GEOLOW(SLAR,SLOR,RH,CLAR,CLOR,RBM,SLAC,SLOC)
        DAT(3,1) = CLAR
        DAT(4,1) = CLOR
        IF(RBM >= 16.) RBM = 999.99
        DAT(5,1) = RBM
        !  Conjugate point coordinates at low latitudes
        WRITE(STR,'(2F6.2)') SLAC,SLOC
        READ (STR,'(2F6.2)') SLAC,SLOC
        DAT(1,2) = SLAC
        DAT(2,2) = SLOC
        CALL GEOCOR(SLAC,SLOC,RH,DAA,DOO,CLAC,CLOC,RBM)
        IF(CLAC > 999.)&
            CALL GEOLOW(SLAC,SLOC,RH,CLAC,CLOC,RBM,SLAL,SLOL)
        DAT(3,2) = CLAC
        DAT(4,2) = CLOC
        DAT(5,2) = RBM
    ELSE
        !  Computation of the magnetically conjugated point at high- and
        !  middle latitudes
        CLAC = -CLAR
        CLOC =  CLOR
        DAT(3,2) = CLAC
        DAT(4,2) = CLOC
        CALL CORGEO(SLAC,SLOC,RH,DAA,DOO,CLAC,CLOC,PMC)
        DAT(1,2) = SLAC
        DAT(2,2) = SLOC
        IF(PMC >= 16.) PMC = 999.99
        DAT(5,2) = PMC
    ENDIF

    !  Same RBM for footprints as for the starting and conjugate points
    DAT(5,3) = DAT(5,1)
    DAT(5,4) = DAT(5,2)

    !  Calculation of the magnetic field line footprint at the
    !  Earth's surface for the starting point
    IF(RH > 1..and.CLAR < 999..and.CLAR < 999.) THEN
        CALL FTPRNT(RH,SLAR,SLOR,CLAR,CLOR,ACLAR,ACLOR,SLARF,SLORF,1.)
        DAT(1,3) = SLARF
        DAT(2,3) = SLORF
        DAT(3,3) = ACLAR
        DAT(4,3) = ACLOR
        !  and for the conjugate point
        CALL FTPRNT(RH,SLAC,SLOC,CLAC,CLOC,ACLAC,ACLOC,SLACF,SLOCF,1.)
        DAT(1,4) = SLACF
        DAT(2,4) = SLOCF
        DAT(3,4) = ACLAC
        DAT(4,4) = ACLOC
    ELSE
        do i = 1,4
            do j = 3,4
                DAT(i,j) = 999.99
            enddo
        enddo
    ENDIF

    !  Computation of geocentric coordinates of the North or South CGM
    !  poles for a given year at the altitude RH and Earth's surface (1-Re)
    CALL CORGEO(PLAN,PLON,RH,DAA,DOO, 90.,360.,PMP)
    PLAN1 = PLAN
    PLON1 = PLON
    CALL CORGEO(PLAS,PLOS,RH,DAA,DOO,-90.,360.,PMP)
    PLAS1 = PLAS
    PLOS1 = PLOS
    IF(RH > 1.) THEN
        CALL CORGEO(PLAN1,PLON1,1.,DAA,DOO, 90.,360.,PMP)
        CALL CORGEO(PLAS1,PLOS1,1.,DAA,DOO,-90.,360.,PMM)
    ENDIF
    IF(CLAR < 0.) THEN
        PLA(1) = PLAS
        PLO(1) = PLOS
    ELSE
        PLA(1) = PLAN
        PLO(1) = PLON
    ENDIF
    IF(ACLAR < 0.) THEN
        PLA(3) = PLAS1
        PLO(3) = PLOS1
    ELSE
        PLA(3) = PLAN1
        PLO(3) = PLON1
    ENDIF
    IF(CLAC < 0.) THEN
        PLA(2) = PLAS
        PLO(2) = PLOS
    ELSE
        PLA(2) = PLAN
        PLO(2) = PLON
    ENDIF
    IF(ACLAC < 0.) THEN
        PLA(4) = PLAS1
        PLO(4) = PLOS1
    ELSE
        PLA(4) = PLAN1
        PLO(4) = PLON1
    ENDIF

    DAT( 6,:) = 99999.
    DAT( 7,:) = 999.99
    DAT( 8,:) = 99999.
    DAT( 9,:) = 999.99
    DAT(10,:) = 999.99
    DAT(11,:) =  99.99

    icount = 2
    if(RH > 1.) icount = 4
    RJ = RH
    do j = 1,icount
        if(j > 2) RJ = 1.
        PLAJ = PLA(j)
        PLOJ = PLO(j)
        SLAJ = DAT(1,j)
        SLOJ = DAT(2,j)
        CLAJ = DAT(3,j)
        CLOJ = DAT(4,j)
        !  Computation of the IGRF components
        CALL MFC(SLAJ,SLOJ,RJ,BTR,BFR,BRR)
        DAT(6,j) = BTR
        DAT(7,j) = BFR
        DAT(8,j) = BRR
        !  Computation of the oval_angle (OVL) between the tangents to
        !  geographic and CGM latitudes at a given point (the code is slightly
        !  modified from the source provided by Therese Morreto in 1994). Note
        !  that rotation of OVL on 90 deg anticlockwise provides the azimuth
        !  to the local "magnetic" north (south) measured from the local
        !  geographic meridian. The OVL_ANG can be calculated only at middle
        !  and high latitudes where CGM --> GEO is permitted.
        OVL = OVL_ANG(SLAJ,SLOJ,CLAJ,CLOJ,RJ)
        DAT(9,j) = OVL
        !  Computation of the meridian_angle (AZM) between the geographic
        !  meridian and direction (azimuth along the great-circle arc) to
        !  the North (South) CGM pole
        AZM = AZM_ANG(SLAJ,SLOJ,CLAJ,PLAJ,PLOJ)
        DAT(10,j) = AZM
        !  Computation of the MLT midnight (in UT)
        CALL MLTUT(SLAJ,SLOJ,CLAJ,PLAJ,PLOJ,UT)
        DAT(11,j) = UT
        !  End of loop j = 1,icount
    enddo
    RETURN
END SUBROUTINE GEOCGM01


!  *********************************************************************
!  This function returns an estimate at the given location of the angle
!  (oval_angle) between the directions (tangents) along the constant
!  CGM and geographic latitudes by utilizing the function DFRIDR from
!  Numerical Recipes for FORTRAN.
!  This angle can be taken as the azimuth to the local "magnetic" north
!  (south) if the eastward (westward) tangent to the local CGM latitude
!  points south (north) from the local geographic latitude.
!  Written by Therese Moretto in August 1994 (revised by V. Papitashvili
!  in January 1999).
!  *********************************************************************
real function OVL_ANG(sla,slo,cla,clo,rr)
    real cgmgla,cgmglo,dfridr
    logical cr360,cr0

    external cgmgla,cgmglo,dfridr
    common/cgmgeo/clat,cr360,cr0,rh
    common/const/umr, pi

    !  Ignore points which nearly coincide with the geographic or CGM poles
    !  within 0.01 degree in latitudes; this also takes care if SLA or CLA
    !  are dummy values (e.g., 999.99)
    if ( abs(sla) >= 89.99     &
         .or. abs(cla) >= 89.99 &
         .or. abs(sla) < 30.) then
        OVL_ANG = 999.99
        return
    endif

    !  Initialize values for the cgmglo and cgmgla functions
    rh = rr
    clat = cla
    cr360 = .false.
    cr0 = .false.

    !  Judge if SLO may be crossing the 360-0 limit. If geocentric
    !  longitude of the location is larger than 270 deg, then cr360 is
    !  set "true"; if it is less than 90 deg, then cr0 is set "true".
    if(slo >= 270.) cr360 = .true.
    if(slo <=  90.)   cr0 = .true.

    !  An initial stepsize (in degrees)
    step = 10.
    !  Note that in the near-pole region the functions CGMGLA and CGMGLO
    !  could be called from DFRIDR with the CGM latitudes exceeded 90 or
    !  -90 degrees (e.g., 98 or -98) when STEP is added or subtracted to a
    !  given CGM latitude (CLA). This does not produce discontinuities in
    !  the functions because GEOCOR calculates GEOLAT smoothly for the
    !  points lying behind the pole (e.g., as for 82 or - 82 deg. in the
    !  above-mentioned example). However, it could be discontinuity in
    !  GEOLON if |GEOLAT| = 90 deg. - see CGMGLO for details.
    hom = dfridr(cgmgla,clo,step,err1)
    denom = dfridr(cgmglo,clo,step,err2)
    denom = denom*cos(sla*pi/180.)
    OVL_ANG = -atan2(hom,denom)
    OVL_ANG = OVL_ANG*57.2957751
    return
end function ovl_ang


!  *********************************************************************
!  This function returns the geocentric latitude as a function of CGM
!  longitude with the CGM latitude held in common block CGMGEO.
!  Essentially this function just calls the subroutine CORGEO.
!  *********************************************************************
real function cgmgla(clon)
    logical cr360,cr0

    common/cgmgeo/cclat,cr360,cr0,rh

    rr = rh
    if(clon > 360.) clon = clon - 360.
    if(clon < 0.) clon = clon + 360.
    call CORGEO(geolat,geolon,rr,dla,dlo,cclat,clon,pmi)
    cgmgla = geolat
    return
end function cgmgla


! *********************************************************************
!  Same as the function CGMGLA but this returns the geocentric
!  longitude. If cr360 is true, geolon+360 deg is returned when geolon
!  is less than 90 deg. If cr0 is true, geolon-360 deg is returned
!  when geolon is larger than 270 degrees.
! *********************************************************************
real function cgmglo(clon)
    logical cr360,cr0

    common/cgmgeo/cclat,cr360,cr0,rh

    rr = rh
    if(clon > 360.) clon = clon - 360.
    if(clon < 0.) clon = clon + 360.

    do
        call CORGEO(geolat,geolon,rr,dla,dlo,cclat,clon,pmi)

        !  Geographic longitude geolon could be any number (e.g., discontinued)
        !  when geolat is the geographic pole
        if (abs(geolat) < 89.99) exit

        clon = clon - 0.01
    end do

    if(cr360.and.(geolon <= 90.)) then
        cgmglo = geolon + 360.
    else
        if (cr0.and.(geolon >= 270.)) then
            cgmglo = geolon - 360.
        else
            cgmglo = geolon
        endif
    endif
    return
end function cgmglo


! **********************************************************************
!  Numerical Recipes Fortran 77 Version 2.07
!  Copyright (c) 1986-1995 by Numerical Recipes Software
! **********************************************************************
FUNCTION DFRIDR(func,x,h,err)
    INTEGER NTAB
    REAL dfridr,err,h,x,func,CON,CON2,BIG,SAFE
    LOGICAL mess
    PARAMETER (CON=1.4,CON2=CON*CON,BIG=1.E30,NTAB=10,SAFE=2.)
    INTEGER i,j
    REAL errt,fac,hh,a(NTAB,NTAB)

    EXTERNAL func
    COMMON/iounit/konsol,mess

    if(h == 0.) then
        if (mess) write(konsol,100)
        return
    endif

    hh = h
    a(1,1) = ( func(x+hh) - func(x-hh) )/(2.0*hh)
    err = BIG

    DO i=2,NTAB
        hh = hh/CON
        a(1,i) = ( func(x+hh) - func(x-hh) )/(2.0*hh)
        fac = CON2
        do j=2,i
            a(j,i) = ( a(j-1,i)*fac - a(j-1,i-1) )/( fac - 1. )
            fac = CON2*fac
            errt = max( abs( a(j,i) - a(j-1,i  ) ), &
                        abs( a(j,i) - a(j-1,i-1) ) )
            if (errt <= err) then
                err = errt
                dfridr = a(j,i)
            endif
        END DO
        if( abs( a(i,i) - a(i-1,i-1) ) >= SAFE*err) return
    END DO

    return
100 FORMAT('h must be nonzero in dfridr')
END FUNCTION DFRIDR


!  *********************************************************************
!  Computation of an angle between the north geographic meridian and
!  direction to the North (South) CGM pole: positive azimuth is
!  measured East (West) from geographic meridian, i.e., the angle is
!  measured between the great-circle arc directions to the geographic
!  and CGM poles. In this case the geomagnetic field components in
!  XYZ (NEV) system can be converted into the CGM system in both
!  hemispheres as:
!                           XM = X cos(alf) + Y sin(alf)
!                           YM =-X sin(alf) + Y cos(alf)
!  Written by V. O. Papitashvili in mid-1980s; revised in February 1999
!  Ignore points which nearly coincide with the geographic or CGM poles
!  within 0.01 degree in latitudes; this also takes care if SLA or CLA
!  are dummy values (e.g., 999.99)
!  *********************************************************************
real function AZM_ANG(sla,slo,cla,pla,plo)
    common/const/ umr, pi

    if(abs(sla) >= 89.99.or.abs(cla) >= 89.99) then
        AZM_ANG = 999.99
        return
    endif

    sp = 1.
    ss = 1.
    if(sign(sp,pla) /= sign(ss,cla)) then
        write(7,2) pla,cla
    endif

    RAD = pi/180.
    am = (90. - abs(pla))*rad
    if(sign(sp,pla) == sign(ss,sla)) then
        cm = (90. - abs(sla))*rad
    else
        cm = (90. + abs(sla))*rad
    endif

    if(sla >= 0.) then
        bet = (plo - slo)*rad
    else
        bet = (slo - plo)*rad
    endif

    sb = sin(bet)
    st = sin(cm)/tan(am) - cos(cm)*cos(bet)
    alfa = atan2(sb,st)
    AZM_ANG = alfa/rad
    RETURN
2   Format(/'WARNING - The CGM pole PLA = ',f6.2,' and station CLAT = ',&
            f6.2,' are not in the same hemisphere: AZM_ANG is incorrect!')
END FUNCTION AZM_ANG


!  *********************************************************************
!  Calculates the MLT midnight in UT hours
!  Definition of the MLT midnight (MLTMN) here is different from the
!  approach described elsewhere. This definition does not take into
!  account the geomagnetic meridian of the subsolar point which causes
!  seasonal variations of the MLTMN in UT time. The latter approach is
!  perfectly applicable to the dipole or eccentric dipole magnetic
!  coordinates but it fails with the CGM coordinates because there are
!  forbidden areas near the geomagnetic equator where CGM coordinates
!  cannot be calculated by definition [e.g., Gustafsson et al., JATP,
!  54, 1609, 1992].
!  In this code the MLT midnight is defined as location of a given point
!  on (or above) the Earth's surface strictly behind the North (South)
!  CGM pole in such the Sun, the pole, and the point are lined up.
!  This approach was originally proposed and coded by Boris Belov
!  sometime in the beginning of 1980s; here it is slightly edited by
!  Vladimir Papitashvili in February 1999.
!  Ignore points which nearly coincide with the geographic or CGM poles
!  within 0.01 degree in latitudes; this also takes care if SLA or CLA
!  are dummy values (e.g., 999.99)
!  *********************************************************************
SUBROUTINE MLTUT(SLA,SLO,CLA,PLA,PLO,UT)
    common/const/umr,pi

    if(abs(sla) >= 89.99.or.abs(cla) >= 89.99) then
        UT = 99.99
        return
    endif

    TPI = 2.*pi 
    RAD = pi/180.
    sp = 1.
    ss = 1.
    if(sign(sp,pla) /= sign(ss,cla)) then
        write(7,2) pla,cla
    endif

    !  Solve the spherical triangle
    QQ = PLO*RAD
    CFF = 90. - abs(PLA)
    CFF = CFF*RAD
    IF(CFF < 0.0000001) CFF=0.0000001
    if(sign(sp,pla) == sign(ss,sla)) then
        CFT = 90. - abs(SLA)
    else
        CFT = 90. + abs(SLA)
    endif
    CFT = CFT*RAD
    IF(CFT < 0.0000001) CFT=0.0000001
    QT = SLO*RAD
    A = SIN(CFF)/SIN(CFT)
    Y = A*SIN(QQ) - SIN(QT)
    X = COS(QT) - A*COS(QQ)
    UT = ATAN2(Y,X)
    IF(UT < 0.) UT = UT + TPI
    QQU = QQ + UT
    QTU = QT + UT
    BP = SIN(CFF)*COS(QQU)
    BT = SIN(CFT)*COS(QTU)
    UT = UT/RAD
    UT = UT/15.

    IF(BP >= BT) then
        IF(UT < 12.) UT = UT + 12.
        IF(UT > 12.) UT = UT - 12.
    END IF

    RETURN

2   format(/'WARNING - The CGM pole PLA = ',f6.2,' and station CLAT = ',&
            f6.2,' are not in the same hemisphere: MLTMN is incorrect!')
END SUBROUTINE MLTUT


!  *********************************************************************
!  Computation of the IGRF magnetic field components
!  Extracted as a subroutine from the earlier version of GEO-CGM.FOR
!  V. Papitashvili, February 1999
!  *********************************************************************
SUBROUTINE MFC(SLA,SLO,R,H,D,Z)
    COMMON /NM/NM
    COMMON /IYR/IYR
    COMMON /CONST/ UMR, PI

    !  This takes care if SLA or CLA are dummy values (e.g., 999.99)
    if(sla >= 999.) then
        X = 99999.
        Y = 99999.
        Z = 99999.
        H = 99999.
        D = 999.99
        I = 999.99
        F = 99999.
        return
    endif

    !  Computation of all geomagnetic field components
    RLA = (90.-SLA)*pi/180.
    RLO = SLO*pi/180.
    CALL IGRF(IYR,NM,R,RLA,RLO,BR,BT,BF)
    X = -BT
    Y =  BF
    Z = -BR
    H = SQRT(X**2 + Y**2)
    D = 180./pi*ATAN2(Y,X)
    I = 180./pi*ATAN2(Z,H)
    F = SQRT(H**2 + Z**2)
    RETURN
END SUBROUTINE MFC


!  *********************************************************************
!  Calculation of the magnetic field line footprint at the Earth's
!  (or any higher) surface.
!  Extracted as a subroutine from the earlier version of GEO-CGM.FOR by
!  V. Papitashvili in February 1999 but then the subroutine was revised
!  to obtain the Altitude Adjusted CGM coordinates. The AACGM approach
!  is proposed by Kile Baker of the JHU/APL, see their World Wide Web
!  site http://sd-www.jhuapl.edu/RADAR/AACGM/ for details.
!  If RF = 1-Re (i.e., at the Earth's surface), then the footprint
!  location is defined as the Altitude Adjusted (AA) CGM coordinates
!  for a given point (ACLA, ACLO).
!
!  If RF = 1.xx Re (i.e., at any altitude above or below the starting
!  point), then the conjunction between these two points can be found
!  along the field line.
!  *********************************************************************
SUBROUTINE FTPRNT(RH,SLA,SLO,CLA,CLO,ACLA,ACLO,SLAF,SLOF,RF)
    logical noloop


    COMMON /NM/NM
    COMMON /IYR/IYR
    COMMON /CONST/ UMR, PI

    !  This takes care if SLA or CLA are dummy values (e.g., 999.99)
    if(sla > 999..or.cla > 999.or.RF == RH) then
        ACLA = 999.99
        ACLO = 999.99
        SLAF = 999.99
        SLOF = 999.99
        return
    endif

    !  Defining the Altitude Adjusted CGM coordinates for a given point
    COL = (90. - CLA)*pi/180.
    SN2 = (SIN(COL))**2
    DECARG=SQRT( SN2*RF/RH )
    IF(ABS(DECARG) > 1.) DECARG=SIGN(1.,DECARG)
    ACOL = ASIN(DECARG)
    ACLA = 90. - ACOL*180./pi
    IF(CLA < 0.) ACLA = -ACLA
    ACLO = CLO
    CALL CORGEO(SLAF,SLOF,RF,DLAF,DLOF,ACLA,ACLO,PMIF)
    IF(SLAF < 999.) RETURN

    !  Tracing the magnetic field line down to the Earth's surface at low
    !  latitudes if CORGEO failed to calculate geocentric coordinates SLAF
    !  and SLOF
    IF(SN2 < 0.0000001) SN2 = 0.0000001
    RL = RH/SN2
    FRAC = 0.03/( 1. + 3./( RL - 0.6 ) )
    !  Checking direction of the magnetic field-line, so the step along
    !  the field-line will go down, to the Earth surface
    IF(CLA >= 0.) FRAC = -FRAC
    DS = RH*FRAC

    do
        !  Start from an initial point
        R = RH
        RSLA = (90. - SLA)*pi/180.
        RSLO = SLO*pi/180.
        CALL SPHCAR(R,RSLA,RSLO,XF,YF,ZF,1)
        RF1 = R
        XF1 = XF
        YF1 = YF
        ZF1 = ZF
        noloop = .false.
        do
            CALL SHAG(XF,YF,ZF,DS)
            RR = SQRT( XF**2 + YF**2 + ZF**2 )
            IF (RR > RH) THEN
                DS = -DS
                XF = XF1
                YF = YF1
                ZF = ZF1
                exit ! inner loop
            ELSE IF (RR > RF) THEN
                RF1 = RR
                XF1 = XF
                YF1 = YF
                ZF1 = ZF
            ELSE
                noloop = .true.
                exit ! inner loop
            ENDIF
        enddo
        if (noloop) exit ! outer loop
    end do

    DR1 = ABS(RF1 - RF)
    DR0 = ABS( RF - RR)
    DR10 = DR1 + DR0
    IF(DR10 /= 0.) THEN
        DS = DS*(DR1/DR10)
        CALL SHAG(XF1,YF1,ZF1,DS)
    ENDIF
    CALL SPHCAR(RR,SLAF,SLOF,XF1,YF1,ZF1,-1)
    SLAF = 90. - SLAF*180./pi
    SLOF = SLOF*180./pi
    RETURN
END SUBROUTINE FTPRNT


!  *********************************************************************
!  Calculates CGM coordinates from geocentric ones at low latitudes
!  where the DGRF/IGRF magnetic field lines may never cross the dipole
!  equatorial plane and, therefore, the definition of CGM coordinates
!  becomes invalid.
!
!  The code is written by Natalia and Vladimir Papitashvili as a part
!  of the earlier versions of GEO-CGM.FOR; extracted as a subroutine by
!  V. Papitashvili in February 1999.
!
!  Apr 11, 2001  GEOLOW is modified to account for interpolation of
!                CGM meridians near equator across the 360/0 boundary
!
!  See the paper by  Gustafsson, G., N. E. Papitashvili, and V. O.
!  Papitashvili, A revised corrected geomagnetic coordinate system for
!  Epochs 1985 and 1990 [J. Atmos. Terr. Phys., 54, 1609-1631, 1992]
!  for detailed description of the B-min approach utilized here.
!  *********************************************************************
SUBROUTINE GEOLOW(SLAR,SLOR,RH,CLAR,CLOR,RBM,SLAC,SLOC)

    DIMENSION BC(2),ARLAT(181),ARLON(181)
    REAL*8 BM,B2,B3
    logical :: NOBM, NOBM_IHEM
    
    COMMON /NM/NM
    COMMON /IYR/IYR
    COMMON /CONST/UMR,PI

    !  This takes care if SLA is a dummy value (e.g., 999.99)
    if(slar > 999.) then
        CLAR = 999.99
        CLOR = 999.99
        SLAC = 999.99
        SLOC = 999.99
        RBM = 999.99
        return
    endif

    !  HH is an error (nT) to determine B-min along the magnetic field line
    DHH = 0.5
    
    !  Filling the work arrays of CGM latitudes and longitudes with 999.99
    !  Note that at certain geocentric longitudes in the very near-equator
    !  region no "geomagnetic equator" can be defined at all.
    ARLAT(61:121) = 999.99
    ARLON(61:121) = 999.99
    SLO  = SLOR

    !  Finding the geomagnetic equator as a projection of the B-min point
    !  found for the field lines started from the last latitude in each
    !  hemisphere where the CGM coordinates were obtained from geocentric
    !  ones (GEO --> CGM). First the CGM coordinates are calculated in the
    !  Northern  and then in the Southern hemispheres 

    !! Northern hemisphere
    !  Program works from 30 deg. latitude down to the geographic equator
    !  in the Northern Hemisphere
    DO JC = 61,91
        SLA = 90. - (JC-1)
        CALL GEOCOR(SLA,SLO,RH,DAA,DOO,CLA,CLO,PMM)
        IF(CLA > 999.) exit
        ARLAT(JC) = CLA
        ARLON(JC) = CLO
    ENDDO

    !! Southern hemisphere
    !  Program works from -30 deg. latitude down to the geographic equator
    !  in the Southern Hemisphere
    DO JC = 121,92,-1
        SLA = 90. - (JC-1)
        CALL GEOCOR(SLA,SLO,RH,DAA,DOO,CLA,CLO,PMM)
        IF(CLA > 999.) exit
        ARLAT(JC) = CLA
        ARLON(JC) = CLO
    ENDDO

    !  Finding last geographic latitudes along SLO where CGM coordinates
    !  can be calculated
    n999=0
    ndir=0
    do jc = 61,121
        if(arlat(jc) > 999.) then
            if(ndir == 0) then
                jcn = jc - 1
                rnlat = arlat(jcn)
                rnlon = arlon(jcn)
                ndir = 1
                n999 = 1
            endif
        endif
        if(arlat(jc) < 999.) then
            if(ndir == 1) then
                jcs = jc
                rslat = arlat(jc)
                rslon = arlon(jc)
                ndir = 0
                exit
            endif
        endif
    enddo

    !  If there is no points with 999.99 found along the SLO meridian,
    !  then the IHEM loop will start from 3; otherwise it starts from 1
    if(n999 == 0) then
        ih = 3
    else
        ih = 1
        !  Interpolation of the appropriate CGM longitudes between last
        !  geocentric latitudes along SLO where CGM coordinates were defined
        ! (modified by Freddy Christiansen of DMI to account for interpolation
        !  across the 360/0 boundary - April 11, 2001)
        rdel = jcs - jcn
        if(rdel == 0.) then
            delon = 0.
        else
            if(rslon > 270. .and. rnlon < 90.) then
                delon = (rslon - (rnlon + 360.))/rdel
            else
                if(rslon < 90..and.rnlon > 270.) then
                    delon = (rslon - (rnlon - 360.))/rdel
                else
                    delon = (rslon - rnlon)/rdel
                endif
            endif
        endif
        do jc = jcn+1,jcs-1
            arlon(jc) = rnlon + delon*(jc-jcn)
            if (arlon(jc) < 0.) arlon(jc) = arlon(jc) + 360.
        enddo
    endif

    !  Finding the CGM equator at SLO on the sphere with radius RH
    NOBM = .FALSE.
    do ihem = ih,3
        RM = RH
        if(ihem == 1) then
            !  Defining the real equator point from the Northern Hemisphere
            CLA = rnlat
            SLA = 90. - (jcn - 1.)
            SLAN = SLA
        else if(ihem == 2) then
            !  Defining the real equator point from the Southern Hemisphere
            CLA = rslat
            SLA = 90. - (jcs - 1)
            SLAS = SLA
        else if(ihem == 3) then
            !  Defining the apex of the current magnetic field line
            CLA = 0.
            SLA = SLAR
        endif

        !  Here CLA is used only to calculate FRAC
        COL = (90. - CLA)*PI/180.
        SLM = (90. - SLA)*PI/180.
        SLL = SLO*PI/180.
        CALL IGRF(IYR,NM,RM,SLM,SLL,BR,BT,BF)
        SZ = -BR
        CALL SPHCAR(RM,SLM,SLL,XGEO,YGEO,ZGEO,1)
        BM = SQRT(BR**2 + BT**2 + BF**2)
        XBM = XGEO
        YBM = YGEO
        ZBM = ZGEO
        RL = 1./(SIN(COL))**2
        FRAC = 0.03/(1. + 3./(RL - 0.6))
        IF(SZ <= 0.) FRAC = -FRAC
        DSD = RL*FRAC
        DS = DSD

        NOBM_IHEM = .FALSE.
        do
            !  Keep two consequently computed points to define B-min
            DO I = 1,2
                DD = DS
                CALL SHAG(XGEO,YGEO,ZGEO,DD)
                IF(I == 1) then
                    XBM1 = XGEO
                    YBM1 = YGEO
                    ZBM1 = ZGEO
                    RBM1 = SQRT(XBM1**2 + YBM1**2 + ZBM1**2)
                endif
                CALL SPHCAR(RM,SLM,SLL,XGEO,YGEO,ZGEO,-1)
                CALL IGRF(IYR,NM,RM,SLM,SLL,BR,BT,BF)
                !  Go and compute the conjugate point if no B-min was found at this
                !  magnetic field line (could happen at very near geomagnetic equator)
                if(RM < RH) then
                    NOBM = .TRUE.
                    NOBM_IHEM = .TRUE.
                    exit
                endif
                BC(I) = SQRT(BR**2 + BT**2 + BF**2)
            end do
            if (NOBM_IHEM) exit

            B2 = BC(1)
            B3 = BC(2)
            IF(B2 >= min(BM, B3)) then 
                IF(B2 >= max(BM, B3) .and. B2 > min(BM, B3)) then
                    XGEO = XBM1
                    YGEO = YBM1
                    ZGEO = ZBM1
                    XBM  = XBM1
                    YBM  = YBM1
                    ZBM  = ZBM1
                    CYCLE
                endif
            else
                BB3 = ABS(B3 - B2)
                BB2 = ABS(BM - B2)
            endif
            IF(BB2 >= DHH .OR. BB3 >= DHH) then
                BM = BM
                XGEO = XBM
                YGEO = YBM
                ZGEO = ZBM
                DS = DS/2.
            else
                exit
            endif
        enddo
        if (NOBM_IHEM) CYCLE
        
        CALL SPHCAR(RBM1,RLA,RLO,XBM1,YBM1,ZBM1,-1)
        RLA = 90. - RLA*180./PI
        RLO = RLO*180./PI
        if(ihem == 1) then
            rlan = rla
        else if(ihem == 2) then
            rlas = rla
        else if(ihem == 3) then
            !  Computation of the magnetically conjugate point at low latitudes
            RBM = RBM1
            RM = RBM
            DS = DSD
            do
                CALL SHAG(XBM1,YBM1,ZBM1,DS)
                RR = SQRT(XBM1**2 + YBM1**2 + ZBM1**2)
                IF (RR > RH) THEN
                    R1 = RR
                    X1 = XBM1
                    Y1 = YBM1
                    Z1 = ZBM1
                ELSE
                    DR1 = ABS(RH - R1)
                    DR0 = ABS(RH - RR)
                    DR10 = DR1 + DR0
                    IF(DR10 /= 0.) THEN
                        DS = DS*(DR1/DR10)
                        RM = R1
                        CALL SHAG(X1,Y1,Z1,DS)
                    ENDIF
                    CALL SPHCAR(RR,SLAC,SLOC,X1,Y1,Z1,-1)
                    SLAC = 90. - SLAC*180./PI
                    SLOC = SLOC*180./PI
                    exit
                ENDIF
            enddo
        endif
    enddo !  End of loop IHEM

    if (n999 /= 0) then
        IF (NOBM) THEN
            !  Interpolation of CGM latitudes if there is no B-min at this
            !  magnetic field line
            rdel = jcs - jcn
            if(rdel == 0.) then
                delat = 0.
            else
                delat = (rslat - rnlat)/rdel
            endif
            jdel = 0
            do jc=jcn+1,jcs-1
                jdel = jdel + 1
                arlat(jc) = rnlat + delat*jdel
            enddo
            RBM = 999.99
            SLAC = 999.99
            SLOC = 999.99
        ELSE
            !  Geocentric latitude of the CGM equator
            rla = (rlan + rlas)/2.
            
            !  Interpolation of the CGM latitudes in the Northern hemisphere
            rdel = SLAN - rla
            if(rdel == 0.) then
                delat = 0.
            else
                delat = rnlat/rdel
            endif
            jdn = abs(rdel)
            jdel = 0
            do jc = jcn+1,jcn+jdn
                jdel = jdel + 1
                arlat(jc) = rnlat - delat*jdel
            enddo

            !  Interpolation of the CGM latitudes in the Southern hemisphere
            rdel = SLAS - rla
            if(rdel == 0.) then
                delat = 0.
            else
                delat = rslat/rdel
            endif
            jds = abs(rdel)
            jdel = 0
            do jc = jcs-1,jcs-jds,-1
                jdel = jdel + 1
                arlat(jc) = rslat + delat*jdel
            enddo
        ENDIF
    endif

    !  Defining by interpolation the exact values of the CGM latitude
    !  and longitude between two adjacent values
    L1 = 90. - SLAR + 1.
    IF(SLAR < 0.) THEN
        L2 = L1-1
    ELSE
        L2 = L1+1
    ENDIF

    DSLA   = ABS(SLAR - INT(SLAR))
    DELCLA = ARLAT(L2) - ARLAT(L1)
    DELCLO = ARLON(L2) - ARLON(L1)
    CLAR   = ARLAT(L1) + DELCLA*DSLA
    CLOR   = ARLON(L1) + DELCLO*DSLA

    RETURN
END SUBROUTINE GEOLOW


!  *********************************************************************
!  Calculates geocentric coordinates from corrected geomagnetic ones.
!  The code is written by Vladimir Popov and Vladimir Papitashvili
!  in mid-1980s; revised by V. Papitashvili in February 1999
!  *********************************************************************
SUBROUTINE CORGEO(SLA,SLO,RH,DLA,DLO,CLA,CLO,PMI)

    COMMON /NM/NM
    COMMON /IYR/IYR
    COMMON /CONST/UMR,PI

    !  This takes care if CLA is a dummy value (e.g., 999.99)
    jc = 0
    if(abs(cla) < 0.1) then
        write(7,2)
        jc = 1
    endif

    if(cla > 999..or.jc == 1) then
        SLA = 999.99
        SLO = 999.99
        DLA = 999.99
        DLO = 999.99
        PMI = 999.99
        return
    endif

    NG = NM
    COL = 90. - CLA
    R = 10.
    R1 = R
    R0 = R
    COL = COL*PI/180.
    RLO = CLO*PI/180.
    SN = SIN(COL)
    SN2 = SN**2

    !  The CGM latitude should be at least 0.01 deg. away of the CGM pole
    SN2 = MAX(SN2, (.01*PI/180.)**2)

    RFI = RH/SN2
    PMI = RFI
    IF(PMI > 99.999) PMI = 999.99
    AA10 = R/RFI

    !  RFI = R if COL = 90 deg.
    IF(RFI > R) then
        SAA = AA10/(1.-AA10)
        SAQ = SQRT(SAA)
        SCLA = ATAN(SAQ)
        IF(CLA < 0) SCLA = PI - SCLA
    else
        SCLA = PI/2
        R0 = RFI
    endif
    CALL SPHCAR(R0,SCLA,RLO,XM,YM,ZM,1)
    CALL GEOMAG(X,Y,Z,XM,YM,ZM,-1,IYR)
    RL = R0
    FRAC = -0.03/(1. + 3./(RL - 0.6))
    IF(CLA < 0.) FRAC = -FRAC
    R = R0
    do
        DS = R*FRAC
        NM = (1. + 9./R) + 0.5
        CALL SHAG(X,Y,Z,DS)
        R = SQRT(X**2 + Y**2 + Z**2)
        IF(R <= RH) exit
        R1 = R
        X1 = X
        Y1 = Y
        Z1 = Z
    enddo

    !  Define intersection with the start surface
    DR1 = ABS(RH - R1)
    DR0 = ABS(RH - R)
    DR10 = DR1 + DR0
    IF(DR10 /= 0.) THEN
        DS = DS*(DR1/DR10)
        CALL SHAG(X1,Y1,Z1,DS)
    ENDIF
    CALL SPHCAR(R,GTET,GXLA,X1,Y1,Z1,-1)
    GTH = GTET*180./PI
    SLO = GXLA*180./PI
    SLA = 90. - GTH
    CALL GEOMAG(X1,Y1,Z1,XM,YM,ZM,1,IYR)
    CALL SPHCAR(RM,TH,PF,XM,YM,ZM,-1)
    DLO = PF*180./PI
    DLA = 90. - TH*180./PI
    NM = NG
    !  Because CORGEO cannot check if the CGM --> GEO transformation is
    !  performed correctly in the equatorial area (that is, where the IGRF
    !  field line may never cross the dipole equatorial plane). Therefore,
    !  the backward check is required for geocentric latitudes lower than
    !  30 degrees (see the paper referenced in GEOLOW)
    IF(ABS(SLA) < 30..OR.ABS(CLA) < 30.) THEN
        CALL GEOCOR(SLA,SLO,RH,DLS,DLS,CLAS,CLOS,PMS)
        IF(CLAS > 999.) CALL GEOLOW(SLA,SLO,RH,CLAS,CLOS,RBM,SLAC,SLOC)
        IF(ABS(ABS(CLA)-ABS(CLAS)) >= 1.) THEN
            write(7,22) CLA
            SLA = 999.99
            SLO = 999.99
            PMI = 999.99
        ENDIF
    ENDIF
    RETURN

! FORMAT STATEMENTS
2    format(/'WARNING - No calculations within +/-0.1 degree near CGM equator')
22   format(/'WARNING - Selected CGM_Lat.=',f6.2,' is too close to geomagnetic'&
            /'          equator where CGM coordinates are not defined')
END SUBROUTINE CORGEO


!  *********************************************************************
!  Calculates corrected geomagnetic coordinates from geocentric ones
!  The code is written by Vladimir Popov and Vladimir Papitashvili
!  in mid-1980s; revised by V. Papitashvili in February 1999
!  *********************************************************************
SUBROUTINE GEOCOR(SLA,SLO,RH,DLA,DLO,CLA,CLO,PMI)

    COMMON /NM/NM
    COMMON /IYR/IYR
    COMMON /CONST/UMR,PI

    !  This takes care if SLA is a dummy value (e.g., 999.99)
    if(sla > 999.) then
        CLA = 999.99
        CLO = 999.99
        DLA = 999.99
        DLO = 999.99
        PMI = 999.99
        return
    endif
    NG = NM
    COL = 90. - SLA
    R = RH
    R1 = R
    COL = COL*PI/180.
    RLO = SLO*PI/180.
    CALL SPHCAR(R,COL,RLO,X,Y,Z,1)
    CALL GEOMAG(X,Y,Z,XM,YM,ZM,1,IYR)
    CALL SPHCAR(RM,TH,PF,XM,YM,ZM,-1)
    SZM = ZM
    DLO = PF*180./PI
    DCO = TH*180./PI
    DLA = 90. - DCO
    RL = R/(SIN(TH))**2
    FRAC = 0.03/(1. + 3./(RL - 0.6))
    IF(SZM < 0.) FRAC = -FRAC
    !  Error to determine the dipole equtorial plane: aprox. 0.5 arc min
    HHH = PI/180. * 30./3600.

    !  Trace the IGRF magnetic field line to the dipole equatorial plane
    RZM = SZM
    do while(SZM*RZM > 0.)
        DS = R*FRAC
        do
            NM = (1. + 9./R) + 0.5
            R1 = R
            X1 = X
            Y1 = Y
            Z1 = Z
            CALL SHAG(X,Y,Z,DS)
            CALL GEOMAG(X,Y,Z,XM,YM,ZM,1,IYR)
            CALL SPHCAR(R,C,S,XM,YM,ZM,-1)

            !  As tracing goes above (RH+10_Re), use the dipole field line
            IF(R > 10.+RH) then
                RZM = -SZM ! Make sure the outer loop terminates
                exit
            endif

            !  If the field line returns to the start surface without crossing the
            !  dipole equatorial plane, no CGM coordinates can be calculated
            IF(R <= RH) then
                CLA = 999.99
                CLO = 999.99
                PMI = 999.99
                NM = NG
                RETURN
            endif

            DCL = C - PI/2
            IF(ABS(DCL) <= HHH) then
                RZM = -SZM ! Make sure the outer loop terminates
                exit 
            endif

            RZM = ZM
            ! If SZM and RZM share a sign, skip assignment and recalculate.
            IF(SZM*RZM > 0.) exit
            R = R1
            X = X1
            Y = Y1
            Z = Z1
            DS = DS/2.
        enddo
    enddo

    CALL GEOMAG(X,Y,Z,XM,YM,ZM,1,IYR)
    CALL SPHCAR(R,GTET,GXLA,XM,YM,ZM,-1)
    ST = ABS(SIN(GTET))
    RRH = ABS(RH/(R - RH*ST**2))
    CLA = PI/2 - ATAN(ST*SQRT(RRH))
    CLA = CLA*180./PI
    CLO = GXLA*180./PI
    IF(SZM < 0.) CLA = -CLA
    SSLA = 90. - CLA
    SSLA = SSLA*PI/180.
    SN = SIN(SSLA)
    PMI = RH/SN**2
    NM = NG
    RETURN
END


!  *********************************************************************
!  Similar to SUBR STEP from GEOPACK-1996 but SHAG takes into account
!  only internal sources
!  The code is re-written from Tsyganenko's subroutine STEP by
!  Natalia and Vladimir Papitashvili in mid-1980s
!  *********************************************************************
SUBROUTINE SHAG(X,Y,Z,DS)

    COMMON/A5/DS3

    REAL :: NX, NY, NZ

    DS3 = -DS/3.
    CALL RIGHT(X, Y, Z, R11, R12, R13)
    NX = X + R11
    NY = Y + R12
    NZ = Z + r13
    CALL RIGHT(NX, NY, NZ, R21, R22, R23)
    NX = X + .5*(R11 + R21)
    NY = Y + .5*(R12 + R22)
    NZ = Z + .5*(R13 + R23)
    CALL RIGHT(NX, NY, NZ, R31, R32, R33)
    NX = X + .375*(R11 + 3.*R31)
    NY = Y + .375*(R12 + 3.*R32)
    NZ = Z + .375*(R13 + 3.*R33)
    CALL RIGHT(NX, NY, NZ, R41, R42, R43)
    NX = X + 1.5*(R11 - 3.*R31 + 4.*R41)
    NY = Y + 1.5*(R12 - 3.*R32 + 4.*R42)
    NZ = Z + 1.5*(R13 - 3.*R33 + 4.*R43)
    CALL RIGHT(NX, NY, NZ, R51, R52, R53)

    X = X + .5*(R11 + 4.*R41 + R51)
    Y = Y + .5*(R12 + 4.*R42 + R52)
    Z = Z + .5*(R13 + 4.*R43 + R53)
    RETURN
END


!  *********************************************************************
!  Similar to SUBR RHAND from GEOPACK-1996 but RIGHT takes into account
!  only internal sources
!  The code is re-written from Tsyganenko's subroutine RHAND
!  by Natalia and Vladimir Papitashvili in mid-1980s
!  *********************************************************************
SUBROUTINE RIGHT(X,Y,Z,R1,R2,R3)

    COMMON /A5/DS3
    COMMON /NM/NM
    COMMON /IYR/IYR

    CALL SPHCAR(R,T,F,X,Y,Z,-1)
    CALL IGRF(IYR,NM,R,T,F,BR,BT,BF)
    CALL BSPCAR(T,F,BR,BT,BF,BX,BY,BZ)

    B  = DS3/SQRT(BX**2 + BY**2 + BZ**2)
    R1 = BX*B
    R2 = BY*B
    R3 = BZ*B
    RETURN
END


!  *********************************************************************
!     CALCULATES COMPONENTS OF THE MAIN (INTERNAL) GEOMAGNETIC FIELD IN SPHERICAL
!     GEOGRAPHICAL COORDINATE SYSTEM, USING IAGA INTERNATIONAL GEOMAGNETIC REFERENCE
!     MODEL COEFFICIENTS (e.g., http://www.ngdc.noaa.gov/IAGA/wg8/igrf2000.html)
!
!     UPDATING THE COEFFICIENTS TO A GIVEN EPOCH IS MADE AUTOMATICALLY UPON THE FIRST
!     CALL AND AFTER EVERY CHANGE OF THE PARAMETER IY.
!
!-----INPUT PARAMETERS:
!
!     IY  -  YEAR NUMBER (FOUR-DIGIT; 1965 &LE IY &LE 2005)
!     NM  -  HIGHEST ORDER OF SPHERICAL HARMONICS IN THE SCALAR POTENTIAL (NM &LE 10)
!     R,T,F -  SPHERICAL COORDINATES (RADIUS R IN UNITS RE=6371.2 KM, GEOGRAPHIC
!                COLATITUDE  T  AND LONGITUDE  F  IN RADIANS)
!
!-----OUTPUT PARAMETERS:
!
!     BR,BT,BF - SPHERICAL COMPONENTS OF THE MAIN GEOMAGNETIC FIELD IN NANOTESLA
!
!     LAST MODIFICATION:  JANUARY 5, 2001, BY: N. A. TSYGANENKO
!     THE CODE WAS MODIFIED TO ACCEPT DATES THROUGH 2005.
!     IT HAS ALSO BEEN SLIGHTLY SIMPLIFIED BY TAKING OUT SOME REDUNDANT STATEMENTS,
!     AND A "SAVE" STATEMENT WAS ADDED, TO AVOID POTENTIAL PROBLEMS WITH SOME
!     FORTRAN COMPILERS.
!
!     MODIFIED TO DGRF TO ACCEPT YEARS FROM 1900 THROUGH 2005
!     BY SCOTT BOARDSEN, NASA GSFC, OCTOBER 2004
!     MODIFIED TO IGRF-10 WITH YEARS THROUGH 2010
!     BY V. PAPITASHVILI, NOVEMBER 2005
!     MODIFIED TO IGRF-11 WITH YEARS THROUGH 2015
!     BY V. PAPITASHVILI, January 2011
!
!  *********************************************************************
SUBROUTINE IGRF(IY,NM,R,T,F,BR,BT,BF)

    SAVE FIRST_RUN,IYR,G,H,REC

    DIMENSION A(11),B(11),G(66),H(66),REC(66)
    logical    mess
    logical :: FIRST_RUN = .TRUE.
    integer :: IYR = 0

    common /iounit/konsol,mess

    ! data from all the years
    REAL, DIMENSION(66) :: G1900 = [&
              0., -31543.,  -2298.,   -677.,   2905.,    924.,   1022.,&
          -1469.,   1256.,    572.,    876.,    628.,    660.,   -361.,&
            134.,   -184.,    328.,    264.,      5.,    -86.,    -16.,&
             63.,     61.,    -11.,   -217.,    -58.,     59.,    -90.,&
             70.,    -55.,      0.,     34.,    -41.,    -21.,     18.,&
              6.,     11.,      8.,     -4.,     -9.,      1.,      2.,&
             -9.,      5.,      8.,      8.,     10.,      1.,    -11.,&
             12.,      1.,     -2.,      2.,     -1.,     -1.,     -3.,&
             -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,&
              2.,      2.,      0.]
    REAL, DIMENSION(66) :: H1900 = [&
              0.,      0.,   5922.,      0.,  -1061.,   1121.,      0.,&
           -330.,      3.,    523.,      0.,    195.,    -69.,   -210.,&
            -75.,      0.,   -210.,     53.,    -33.,   -124.,      3.,&
              0.,     -9.,     83.,      2.,    -35.,     36.,    -69.,&
              0.,    -45.,    -13.,    -10.,     -1.,     28.,    -12.,&
            -22.,      0.,      8.,    -14.,      7.,    -13.,      5.,&
             16.,     -5.,    -18.,      0.,    -20.,     14.,      5.,&
             -3.,     -2.,      8.,     10.,     -2.,      2.,      0.,&
              2.,      1.,      2.,      6.,     -4.,      0.,     -2.,&
              4.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1905 = [&
              0., -31464.,  -2298.,   -728.,   2928.,   1041.,   1037.,&
          -1494.,   1239.,    635.,    880.,    643.,    653.,   -380.,&
            146.,   -192.,    328.,    259.,     -1.,    -93.,    -26.,&
             62.,     60.,    -11.,   -221.,    -57.,     57.,    -92.,&
             70.,    -54.,      0.,     33.,    -41.,    -20.,     18.,&
              6.,     11.,      8.,     -4.,     -9.,      1.,      2.,&
             -8.,      5.,      8.,      8.,     10.,      1.,    -11.,&
             12.,      1.,     -2.,      2.,      0.,     -1.,     -3.,&
             -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,&
              2.,      2.,      0.]
    REAL, DIMENSION(66) :: H1905 = [&
              0.,      0.,   5909.,      0.,  -1086.,   1065.,      0.,&
           -357.,     34.,    480.,      0.,    203.,    -77.,   -201.,&
            -65.,      0.,   -193.,     56.,    -32.,   -125.,     11.,&
              0.,     -7.,     86.,      4.,    -32.,     32.,    -67.,&
              0.,    -46.,    -14.,    -11.,      0.,     28.,    -12.,&
            -22.,      0.,      8.,    -15.,      7.,    -13.,      5.,&
             16.,     -5.,    -18.,      0.,    -20.,     14.,      5.,&
             -3.,     -2.,      8.,     10.,     -2.,      2.,      0.,&
              2.,      1.,      2.,      6.,     -4.,      0.,     -2.,&
              4.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1910 = [&
              0., -31354.,  -2297.,   -769.,   2948.,   1176.,   1058.,&
          -1524.,   1223.,    705.,    884.,    660.,    644.,   -400.,&
            160.,   -201.,    327.,    253.,     -9.,   -102.,    -38.,&
             62.,     58.,    -11.,   -224.,    -54.,     54.,    -95.,&
             71.,    -54.,      1.,     32.,    -40.,    -19.,     18.,&
              6.,     11.,      8.,     -4.,     -9.,      1.,      2.,&
             -8.,      5.,      8.,      8.,     10.,      1.,    -11.,&
             12.,      1.,     -2.,      2.,      0.,     -1.,     -3.,&
             -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,&
              2.,      2.,      0.]
    REAL, DIMENSION(66) :: H1910 = [&
              0.,      0.,   5898.,      0.,  -1128.,   1000.,      0.,&
           -389.,     62.,    425.,      0.,    211.,    -90.,   -189.,&
            -55.,      0.,   -172.,     57.,    -33.,   -126.,     21.,&
              0.,     -5.,     89.,      5.,    -29.,     28.,    -65.,&
              0.,    -47.,    -14.,    -12.,      1.,     28.,    -13.,&
            -22.,      0.,      8.,    -15.,      6.,    -13.,      5.,&
             16.,     -5.,    -18.,      0.,    -20.,     14.,      5.,&
             -3.,     -2.,      8.,     10.,     -2.,      2.,      0.,&
              2.,      1.,      2.,      6.,     -4.,      0.,     -2.,&
              4.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1915 = [&
              0., -31212.,  -2306.,   -802.,   2956.,   1309.,   1084.,&
          -1559.,   1212.,    778.,    887.,    678.,    631.,   -416.,&
            178.,   -211.,    327.,    245.,    -16.,   -111.,    -51.,&
             61.,     57.,    -10.,   -228.,    -51.,     49.,    -98.,&
             72.,    -54.,      2.,     31.,    -38.,    -18.,     19.,&
              6.,     11.,      8.,     -4.,     -9.,      2.,      3.,&
             -8.,      6.,      8.,      8.,     10.,      1.,    -11.,&
             12.,      1.,     -2.,      2.,      0.,     -1.,     -3.,&
             -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,&
              1.,      2.,      0.]
    REAL, DIMENSION(66) :: H1915 = [&
              0.,      0.,   5875.,      0.,  -1191.,    917.,      0.,&
           -421.,     84.,    360.,      0.,    218.,   -109.,   -173.,&
            -51.,      0.,   -148.,     58.,    -34.,   -126.,     32.,&
              0.,     -2.,     93.,      8.,    -26.,     23.,    -62.,&
              0.,    -48.,    -14.,    -12.,      2.,     28.,    -15.,&
            -22.,      0.,      8.,    -15.,      6.,    -13.,      5.,&
             16.,     -5.,    -18.,      0.,    -20.,     14.,      5.,&
             -3.,     -2.,      8.,     10.,     -2.,      2.,      0.,&
              2.,      1.,      2.,      6.,     -4.,      0.,     -2.,&
              4.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1920 = [&
              0., -31060.,  -2317.,   -839.,   2959.,   1407.,   1111.,&
          -1600.,   1205.,    839.,    889.,    695.,    616.,   -424.,&
            199.,   -221.,    326.,    236.,    -23.,   -119.,    -62.,&
             61.,     55.,    -10.,   -233.,    -46.,     44.,   -101.,&
             73.,    -54.,      2.,     29.,    -37.,    -16.,     19.,&
              6.,     11.,      7.,     -3.,     -9.,      2.,      4.,&
             -7.,      6.,      8.,      8.,     10.,      1.,    -11.,&
             12.,      1.,     -2.,      2.,      0.,     -1.,     -3.,&
             -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,&
              1.,      3.,      0.]
    REAL, DIMENSION(66) :: H1920 = [&
              0.,      0.,   5845.,      0.,  -1259.,    823.,      0.,&
           -445.,    103.,    293.,      0.,    220.,   -134.,   -153.,&
            -57.,      0.,   -122.,     58.,    -38.,   -125.,     43.,&
              0.,      0.,     96.,     11.,    -22.,     18.,    -57.,&
              0.,    -49.,    -14.,    -13.,      4.,     28.,    -16.,&
            -22.,      0.,      8.,    -15.,      6.,    -14.,      5.,&
             17.,     -5.,    -19.,      0.,    -20.,     14.,      5.,&
             -3.,     -2.,      9.,     10.,     -2.,      2.,      0.,&
              2.,      1.,      2.,      6.,     -4.,      0.,     -2.,&
              4.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1925 = [&
              0., -30926.,  -2318.,   -893.,   2969.,   1471.,   1140.,&
          -1645.,   1202.,    881.,    891.,    711.,    601.,   -426.,&
            217.,   -230.,    326.,    226.,    -28.,   -125.,    -69.,&
             61.,     54.,     -9.,   -238.,    -40.,     39.,   -103.,&
             73.,    -54.,      3.,     27.,    -35.,    -14.,     19.,&
              6.,     11.,      7.,     -3.,     -9.,      2.,      4.,&
             -7.,      7.,      8.,      8.,     10.,      1.,    -11.,&
             12.,      1.,     -2.,      2.,      0.,     -1.,     -3.,&
             -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,&
              1.,      3.,      0.]
    REAL, DIMENSION(66) :: H1925 = [&
              0.,      0.,   5817.,      0.,  -1334.,    728.,      0.,&
           -462.,    119.,    229.,      0.,    216.,   -163.,   -130.,&
            -70.,      0.,    -96.,     58.,    -44.,   -122.,     51.,&
              0.,      3.,     99.,     14.,    -18.,     13.,    -52.,&
              0.,    -50.,    -14.,    -14.,      5.,     29.,    -17.,&
            -21.,      0.,      8.,    -15.,      6.,    -14.,      5.,&
             17.,     -5.,    -19.,      0.,    -20.,     14.,      5.,&
             -3.,     -2.,      9.,     10.,     -2.,      2.,      0.,&
              2.,      1.,      2.,      6.,     -4.,      0.,     -2.,&
              4.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1930 = [&
              0., -30805.,  -2316.,   -951.,   2980.,   1517.,   1172.,&
          -1692.,   1205.,    907.,    896.,    727.,    584.,   -422.,&
            234.,   -237.,    327.,    218.,    -32.,   -131.,    -74.,&
             60.,     53.,     -9.,   -242.,    -32.,     32.,   -104.,&
             74.,    -54.,      4.,     25.,    -34.,    -12.,     18.,&
              6.,     11.,      7.,     -3.,     -9.,      2.,      5.,&
             -6.,      8.,      8.,      8.,     10.,      1.,    -12.,&
             12.,      1.,     -2.,      3.,      0.,     -2.,     -3.,&
             -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,&
              1.,      3.,      0.]
    REAL, DIMENSION(66) :: H1930 = [&
              0.,      0.,   5808.,      0.,  -1424.,    644.,      0.,&
           -480.,    133.,    166.,      0.,    205.,   -195.,   -109.,&
            -90.,      0.,    -72.,     60.,    -53.,   -118.,     58.,&
              0.,      4.,    102.,     19.,    -16.,      8.,    -46.,&
              0.,    -51.,    -15.,    -14.,      6.,     29.,    -18.,&
            -20.,      0.,      8.,    -15.,      5.,    -14.,      5.,&
             18.,     -5.,    -19.,      0.,    -20.,     14.,      5.,&
             -3.,     -2.,      9.,     10.,     -2.,      2.,      0.,&
              2.,      1.,      2.,      6.,     -4.,      0.,     -2.,&
              4.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1935 = [&
              0., -30715.,  -2306.,  -1018.,   2984.,   1550.,   1206.,&
          -1740.,   1215.,    918.,    903.,    744.,    565.,   -415.,&
            249.,   -241.,    329.,    211.,    -33.,   -136.,    -76.,&
             59.,     53.,     -8.,   -246.,    -25.,     25.,   -106.,&
             74.,    -53.,      4.,     23.,    -33.,    -11.,     18.,&
              6.,     11.,      7.,     -3.,     -9.,      1.,      6.,&
             -6.,      8.,      7.,      8.,     10.,      1.,    -12.,&
             11.,      1.,     -2.,      3.,      0.,     -2.,     -3.,&
             -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,&
              2.,      3.,      0.]
    REAL, DIMENSION(66) :: H1935 = [&
              0.,      0.,   5812.,      0.,  -1520.,    586.,      0.,&
           -494.,    146.,    101.,      0.,    188.,   -226.,    -90.,&
           -114.,      0.,    -51.,     64.,    -64.,   -115.,     64.,&
              0.,      4.,    104.,     25.,    -15.,      4.,    -40.,&
              0.,    -52.,    -17.,    -14.,      7.,     29.,    -19.,&
            -19.,      0.,      8.,    -15.,      5.,    -15.,      5.,&
             18.,     -5.,    -19.,      0.,    -20.,     15.,      5.,&
             -3.,     -3.,      9.,     11.,     -2.,      2.,      0.,&
              2.,      1.,      2.,      6.,     -4.,      0.,     -1.,&
              4.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1940 = [&
              0., -30654.,  -2292.,  -1106.,   2981.,   1566.,   1240.,&
          -1790.,   1232.,    916.,    914.,    762.,    550.,   -405.,&
            265.,   -241.,    334.,    208.,    -33.,   -141.,    -76.,&
             57.,     54.,     -7.,   -249.,    -18.,     18.,   -107.,&
             74.,    -53.,      4.,     20.,    -31.,     -9.,     17.,&
              5.,     11.,      7.,     -3.,    -10.,      1.,      6.,&
             -5.,      9.,      7.,      8.,     10.,      1.,    -12.,&
             11.,      1.,     -2.,      3.,      1.,     -2.,     -3.,&
             -4.,      2.,     -5.,     -2.,      6.,      4.,      0.,&
              2.,      3.,      0.]
    REAL, DIMENSION(66) :: H1940 = [&
              0.,      0.,   5821.,      0.,  -1614.,    528.,      0.,&
           -499.,    163.,     43.,      0.,    169.,   -252.,    -72.,&
           -141.,      0.,    -33.,     71.,    -75.,   -113.,     69.,&
              0.,      4.,    105.,     33.,    -15.,      0.,    -33.,&
              0.,    -52.,    -18.,    -14.,      7.,     29.,    -20.,&
            -19.,      0.,      8.,    -14.,      5.,    -15.,      5.,&
             19.,     -5.,    -19.,      0.,    -21.,     15.,      5.,&
             -3.,     -3.,      9.,     11.,     -2.,      2.,      0.,&
              2.,      1.,      2.,      6.,     -4.,      0.,     -1.,&
              4.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1945 = [&
              0., -30594.,  -2285.,  -1244.,   2990.,   1578.,   1282.,&
          -1834.,   1255.,    913.,    944.,    776.,    544.,   -421.,&
            304.,   -253.,    346.,    194.,    -20.,   -142.,    -82.,&
             59.,     57.,      6.,   -246.,    -25.,     21.,   -104.,&
             70.,    -40.,      0.,      0.,    -29.,    -10.,     15.,&
             29.,     13.,      7.,     -8.,     -5.,      9.,      7.,&
            -10.,      7.,      2.,      5.,    -21.,      1.,    -11.,&
              3.,     16.,     -3.,     -4.,     -3.,     -4.,     -3.,&
             11.,      1.,      2.,     -5.,     -1.,      8.,     -1.,&
             -3.,      5.,     -2.]
    REAL, DIMENSION(66) :: H1945 = [&
              0.,      0.,   5810.,      0.,  -1702.,    477.,      0.,&
           -499.,    186.,    -11.,      0.,    144.,   -276.,    -55.,&
           -178.,      0.,    -12.,     95.,    -67.,   -119.,     82.,&
              0.,      6.,    100.,     16.,     -9.,    -16.,    -39.,&
              0.,    -45.,    -18.,      2.,      6.,     28.,    -17.,&
            -22.,      0.,     12.,    -21.,    -12.,     -7.,      2.,&
             18.,      3.,    -11.,      0.,    -27.,     17.,     29.,&
             -9.,      4.,      9.,      6.,      1.,      8.,      0.,&
              5.,      1.,    -20.,     -1.,     -6.,      6.,     -4.,&
             -2.,      0.,     -2.]
    REAL, DIMENSION(66) :: G1950 = [&
              0., -30554.,  -2250.,  -1341.,   2998.,   1576.,   1297.,&
          -1889.,   1274.,    896.,    954.,    792.,    528.,   -408.,&
            303.,   -240.,    349.,    211.,    -20.,   -147.,    -76.,&
             54.,     57.,      4.,   -247.,    -16.,     12.,   -105.,&
             65.,    -55.,      2.,      1.,    -40.,     -7.,      5.,&
             19.,     22.,     15.,     -4.,     -1.,     11.,     15.,&
            -13.,      5.,     -1.,      3.,     -7.,     -1.,    -25.,&
             10.,      5.,     -5.,     -2.,      3.,      8.,     -8.,&
              4.,     -1.,     13.,     -4.,      4.,     12.,      3.,&
              2.,     10.,      3.]
    REAL, DIMENSION(66) :: H1950 = [&
              0.,      0.,   5815.,      0.,  -1810.,    381.,      0.,&
           -476.,    206.,    -46.,      0.,    136.,   -278.,    -37.,&
           -210.,      0.,      3.,    103.,    -87.,   -122.,     80.,&
              0.,     -1.,     99.,     33.,    -12.,    -12.,    -30.,&
              0.,    -35.,    -17.,      0.,     10.,     36.,    -18.,&
            -16.,      0.,      5.,    -22.,      0.,    -21.,     -8.,&
             17.,     -4.,    -17.,      0.,    -24.,     19.,     12.,&
              2.,      2.,      8.,      8.,    -11.,     -7.,      0.,&
             13.,     -2.,    -10.,      2.,     -3.,      6.,     -3.,&
              6.,     11.,      8.]
    REAL, DIMENSION(66) :: G1955 = [&
              0., -30500.,  -2215.,  -1440.,   3003.,   1581.,   1302.,&
          -1944.,   1288.,    882.,    958.,    796.,    510.,   -397.,&
            290.,   -229.,    360.,    230.,    -23.,   -152.,    -69.,&
             47.,     57.,      3.,   -247.,     -8.,      7.,   -107.,&
             65.,    -56.,      2.,     10.,    -32.,    -11.,      9.,&
             18.,     11.,      9.,     -6.,    -14.,      6.,     10.,&
             -7.,      6.,      9.,      4.,      9.,     -4.,     -5.,&
              2.,      4.,      1.,      2.,      2.,      5.,     -3.,&
             -5.,     -1.,      2.,     -3.,      7.,      4.,     -2.,&
              6.,     -2.,      0.]
    REAL, DIMENSION(66) :: H1955 = [&
              0.,      0.,   5820.,      0.,  -1898.,    291.,      0.,&
           -462.,    216.,    -83.,      0.,    133.,   -274.,    -23.,&
           -230.,      0.,     15.,    110.,    -98.,   -121.,     78.,&
              0.,     -9.,     96.,     48.,    -16.,    -12.,    -24.,&
              0.,    -50.,    -24.,     -4.,      8.,     28.,    -20.,&
            -18.,      0.,     10.,    -15.,      5.,    -23.,      3.,&
             23.,     -4.,    -13.,      0.,    -11.,     12.,      7.,&
              6.,     -2.,     10.,      7.,     -6.,      5.,      0.,&
             -4.,      0.,     -8.,     -2.,     -4.,      1.,     -3.,&
              7.,     -1.,     -3.]
    REAL, DIMENSION(66) :: G1960 = [&
              0., -30421.,  -2169.,  -1555.,   3002.,   1590.,   1302.,&
          -1992.,   1289.,    878.,    957.,    800.,    504.,   -394.,&
            269.,   -222.,    362.,    242.,    -26.,   -156.,    -63.,&
             46.,     58.,      1.,   -237.,     -1.,     -2.,   -113.,&
             67.,    -56.,      5.,     15.,    -32.,     -7.,     17.,&
              8.,     15.,      6.,     -4.,    -11.,      2.,     10.,&
             -5.,     10.,      8.,      4.,      6.,      0.,     -9.,&
              1.,      4.,     -1.,     -2.,      3.,     -1.,      1.,&
             -3.,      4.,      0.,     -1.,      4.,      6.,      1.,&
             -1.,      2.,      0.]
    REAL, DIMENSION(66) :: H1960 = [&
              0.,      0.,   5791.,      0.,  -1967.,    206.,      0.,&
           -414.,    224.,   -130.,      0.,    135.,   -278.,      3.,&
           -255.,      0.,     16.,    125.,   -117.,   -114.,     81.,&
              0.,    -10.,     99.,     60.,    -20.,    -11.,    -17.,&
              0.,    -55.,    -28.,     -6.,      7.,     23.,    -18.,&
            -17.,      0.,     11.,    -14.,      7.,    -18.,      4.,&
             23.,      1.,    -20.,      0.,    -18.,     12.,      2.,&
              0.,     -3.,      9.,      8.,      0.,      5.,      0.,&
              4.,      1.,      0.,      2.,     -5.,      1.,     -1.,&
              6.,      0.,     -7.]
    REAL, DIMENSION(66) :: G1965 = [&
              0., -30334.,  -2119.,  -1662.,   2997.,   1594.,   1297.,&
          -2038.,   1292.,    856.,    957.,    804.,    479.,   -390.,&
            252.,   -219.,    358.,    254.,    -31.,   -157.,    -62.,&
             45.,     61.,      8.,   -228.,      4.,      1.,   -111.,&
             75.,    -57.,      4.,     13.,    -26.,     -6.,     13.,&
              1.,     13.,      5.,     -4.,    -14.,      0.,      8.,&
             -1.,     11.,      4.,      8.,     10.,      2.,    -13.,&
             10.,     -1.,     -1.,      5.,      1.,     -2.,     -2.,&
             -3.,      2.,     -5.,     -2.,      4.,      4.,      0.,&
              2.,      2.,      0.]
    REAL, DIMENSION(66) :: H1965 = [&
              0.,      0.,   5776.,      0.,  -2016.,    114.,      0.,&
           -404.,    240.,   -165.,      0.,    148.,   -269.,     13.,&
           -269.,      0.,     19.,    128.,   -126.,    -97.,     81.,&
              0.,    -11.,    100.,     68.,    -32.,     -8.,     -7.,&
              0.,    -61.,    -27.,     -2.,      6.,     26.,    -23.,&
            -12.,      0.,      7.,    -12.,      9.,    -16.,      4.,&
             24.,     -3.,    -17.,      0.,    -22.,     15.,      7.,&
             -4.,     -5.,     10.,     10.,     -4.,      1.,      0.,&
              2.,      1.,      2.,      6.,     -4.,      0.,     -2.,&
              3.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1970 = [&
              0., -30220.,  -2068.,  -1781.,   3000.,   1611.,   1287.,&
          -2091.,   1278.,    838.,    952.,    800.,    461.,   -395.,&
            234.,   -216.,    359.,    262.,    -42.,   -160.,    -56.,&
             43.,     64.,     15.,   -212.,      2.,      3.,   -112.,&
             72.,    -57.,      1.,     14.,    -22.,     -2.,     13.,&
             -2.,     14.,      6.,     -2.,    -13.,     -3.,      5.,&
              0.,     11.,      3.,      8.,     10.,      2.,    -12.,&
             10.,     -1.,      0.,      3.,      1.,     -1.,     -3.,&
             -3.,      2.,     -5.,     -1.,      6.,      4.,      1.,&
              0.,      3.,     -1.]
    REAL, DIMENSION(66) :: H1970 = [&
              0.,      0.,   5737.,      0.,  -2047.,     25.,      0.,&
           -366.,    251.,   -196.,      0.,    167.,   -266.,     26.,&
           -279.,      0.,     26.,    139.,   -139.,    -91.,     83.,&
              0.,    -12.,    100.,     72.,    -37.,     -6.,      1.,&
              0.,    -70.,    -27.,     -4.,      8.,     23.,    -23.,&
            -11.,      0.,      7.,    -15.,      6.,    -17.,      6.,&
             21.,     -6.,    -16.,      0.,    -21.,     16.,      6.,&
             -4.,     -5.,     10.,     11.,     -2.,      1.,      0.,&
              1.,      1.,      3.,      4.,     -4.,      0.,     -1.,&
              3.,      1.,     -4.]
    REAL, DIMENSION(66) :: G1975 = [&
              0., -30100.,  -2013.,  -1902.,   3010.,   1632.,   1276.,&
          -2144.,   1260.,    830.,    946.,    791.,    438.,   -405.,&
            216.,   -218.,    356.,    264.,    -59.,   -159.,    -49.,&
             45.,     66.,     28.,   -198.,      1.,      6.,   -111.,&
             71.,    -56.,      1.,     16.,    -14.,      0.,     12.,&
             -5.,     14.,      6.,     -1.,    -12.,     -8.,      4.,&
              0.,     10.,      1.,      7.,     10.,      2.,    -12.,&
             10.,     -1.,     -1.,      4.,      1.,     -2.,     -3.,&
             -3.,      2.,     -5.,     -2.,      5.,      4.,      1.,&
              0.,      3.,     -1.]
    REAL, DIMENSION(66) :: H1975 = [&
              0.,      0.,   5675.,      0.,  -2067.,    -68.,      0.,&
           -333.,    262.,   -223.,      0.,    191.,   -265.,     39.,&
           -288.,      0.,     31.,    148.,   -152.,    -83.,     88.,&
              0.,    -13.,     99.,     75.,    -41.,     -4.,     11.,&
              0.,    -77.,    -26.,     -5.,     10.,     22.,    -23.,&
            -12.,      0.,      6.,    -16.,      4.,    -19.,      6.,&
             18.,    -10.,    -17.,      0.,    -21.,     16.,      7.,&
             -4.,     -5.,     10.,     11.,     -3.,      1.,      0.,&
              1.,      1.,      3.,      4.,     -4.,     -1.,     -1.,&
              3.,      1.,     -5.]
    REAL, DIMENSION(66) :: G1980 = [&
              0., -29992.,  -1956.,  -1997.,   3027.,   1663.,   1281.,&
          -2180.,   1251.,    833.,    938.,    782.,    398.,   -419.,&
            199.,   -218.,    357.,    261.,    -74.,   -162.,    -48.,&
             48.,     66.,     42.,   -192.,      4.,     14.,   -108.,&
             72.,    -59.,      2.,     21.,    -12.,      1.,     11.,&
             -2.,     18.,      6.,      0.,    -11.,     -7.,      4.,&
              3.,      6.,     -1.,      5.,     10.,      1.,    -12.,&
              9.,     -3.,     -1.,      7.,      2.,     -5.,     -4.,&
             -4.,      2.,     -5.,     -2.,      5.,      3.,      1.,&
              2.,      3.,      0.]
    REAL, DIMENSION(66) :: H1980 = [&
              0.,      0.,   5604.,      0.,  -2129.,   -200.,      0.,&
           -336.,    271.,   -252.,      0.,    212.,   -257.,     53.,&
           -297.,      0.,     46.,    150.,   -151.,    -78.,     92.,&
              0.,    -15.,     93.,     71.,    -43.,     -2.,     17.,&
              0.,    -82.,    -27.,     -5.,     16.,     18.,    -23.,&
            -10.,      0.,      7.,    -18.,      4.,    -22.,      9.,&
             16.,    -13.,    -15.,      0.,    -21.,     16.,      9.,&
             -5.,     -6.,      9.,     10.,     -6.,      2.,      0.,&
              1.,      0.,      3.,      6.,     -4.,      0.,     -1.,&
              4.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1985 = [&
              0., -29873.,  -1905.,  -2072.,   3044.,   1687.,   1296.,&
          -2208.,   1247.,    829.,    936.,    780.,    361.,   -424.,&
            170.,   -214.,    355.,    253.,    -93.,   -164.,    -46.,&
             53.,     65.,     51.,   -185.,      4.,     16.,   -102.,&
             74.,    -62.,      3.,     24.,     -6.,      4.,     10.,&
              0.,     21.,      6.,      0.,    -11.,     -9.,      4.,&
              4.,      4.,     -4.,      5.,     10.,      1.,    -12.,&
              9.,     -3.,     -1.,      7.,      1.,     -5.,     -4.,&
             -4.,      3.,     -5.,     -2.,      5.,      3.,      1.,&
              2.,      3.,      0.]
    REAL, DIMENSION(66) :: H1985 = [&
              0.,      0.,   5500.,      0.,  -2197.,   -306.,      0.,&
           -310.,    284.,   -297.,      0.,    232.,   -249.,     69.,&
           -297.,      0.,     47.,    150.,   -154.,    -75.,     95.,&
              0.,    -16.,     88.,     69.,    -48.,     -1.,     21.,&
              0.,    -83.,    -27.,     -2.,     20.,     17.,    -23.,&
             -7.,      0.,      8.,    -19.,      5.,    -23.,     11.,&
             14.,    -15.,    -11.,      0.,    -21.,     15.,      9.,&
             -6.,     -6.,      9.,      9.,     -7.,      2.,      0.,&
              1.,      0.,      3.,      6.,     -4.,      0.,     -1.,&
              4.,      0.,     -6.]
    REAL, DIMENSION(66) :: G1990 = [&
              0., -29775.,  -1848.,  -2131.,   3059.,   1686.,   1314.,&
          -2239.,   1248.,    802.,    939.,    780.,    325.,   -423.,&
            141.,   -214.,    353.,    245.,   -109.,   -165.,    -36.,&
             61.,     65.,     59.,   -178.,      3.,     18.,    -96.,&
             77.,    -64.,      2.,     26.,     -1.,      5.,      9.,&
              0.,     23.,      5.,     -1.,    -10.,    -12.,      3.,&
              4.,      2.,     -6.,      4.,      9.,      1.,    -12.,&
              9.,     -4.,     -2.,      7.,      1.,     -6.,     -3.,&
             -4.,      2.,     -5.,     -2.,      4.,      3.,      1.,&
              3.,      3.,      0.]
    REAL, DIMENSION(66) :: H1990 = [&
              0.,      0.,   5406.,      0.,  -2279.,   -373.,      0.,&
           -284.,    293.,   -352.,      0.,    247.,   -240.,     84.,&
           -299.,      0.,     46.,    154.,   -153.,    -69.,     97.,&
              0.,    -16.,     82.,     69.,    -52.,      1.,     24.,&
              0.,    -80.,    -26.,      0.,     21.,     17.,    -23.,&
             -4.,      0.,     10.,    -19.,      6.,    -22.,     12.,&
             12.,    -16.,    -10.,      0.,    -20.,     15.,     11.,&
             -7.,     -7.,      9.,      8.,     -7.,      2.,      0.,&
              2.,      1.,      3.,      6.,     -4.,      0.,     -2.,&
              3.,     -1.,     -6.]
    REAL, DIMENSION(66) :: G1995 = [&
              0., -29692.,  -1784.,  -2200.,   3070.,   1681.,   1335.,&
          -2267.,   1249.,    759.,    940.,    780.,    290.,   -418.,&
            122.,   -214.,    352.,    235.,   -118.,   -166.,    -17.,&
             68.,     67.,     68.,   -170.,     -1.,     19.,    -93.,&
             77.,    -72.,      1.,     28.,      5.,      4.,      8.,&
             -2.,     25.,      6.,     -6.,     -9.,    -14.,      9.,&
              6.,     -5.,     -7.,      4.,      9.,      3.,    -10.,&
              8.,     -8.,     -1.,     10.,     -2.,     -8.,     -3.,&
             -6.,      2.,     -4.,     -1.,      4.,      2.,      2.,&
              5.,      1.,      0.]
    REAL, DIMENSION(66) :: H1995 = [&
              0.,      0.,   5306.,      0.,  -2366.,   -413.,      0.,&
           -262.,    302.,   -427.,      0.,    262.,   -236.,     97.,&
           -306.,      0.,     46.,    165.,   -143.,    -55.,    107.,&
              0.,    -17.,     72.,     67.,    -58.,      1.,     36.,&
              0.,    -69.,    -25.,      4.,     24.,     17.,    -24.,&
             -6.,      0.,     11.,    -21.,      8.,    -23.,     15.,&
             11.,    -16.,     -4.,      0.,    -20.,     15.,     12.,&
             -6.,     -8.,      8.,      5.,     -8.,      3.,      0.,&
              1.,      0.,      4.,      5.,     -5.,     -1.,     -2.,&
              1.,     -2.,     -7.]
    REAL, DIMENSION(66) :: G2000 = [&
             0.0,-29619.4, -1728.2, -2267.7,  3068.4,  1670.9,  1339.6,&
         -2288.0,  1252.1,   714.5,   932.3,   786.8,   250.0,  -403.0,&
           111.3,  -218.8,   351.4,   222.3,  -130.4,  -168.6,   -12.9,&
            72.3,    68.2,    74.2,  -160.9,    -5.9,    16.9,   -90.4,&
            79.0,   -74.0,     0.0,    33.3,     9.1,     6.9,     7.3,&
            -1.2,    24.4,     6.6,    -9.2,    -7.9,   -16.6,     9.1,&
             7.0,    -7.9,    -7.0,     5.0,     9.4,     3.0,    -8.4,&
             6.3,    -8.9,    -1.5,     9.3,    -4.3,    -8.2,    -2.6,&
            -6.0,     1.7,    -3.1,    -0.5,     3.7,     1.0,     2.0,&
             4.2,     0.3,    -1.1]
    REAL, DIMENSION(66) :: H2000 = [&
             0.0,     0.0,  5186.1,     0.0, -2481.6,  -458.0,     0.0,&
          -227.6,   293.4,  -491.1,     0.0,   272.6,  -231.9,   119.8,&
          -303.8,     0.0,    43.8,   171.9,  -133.1,   -39.3,   106.3,&
             0.0,   -17.4,    63.7,    65.1,   -61.2,     0.7,    43.8,&
             0.0,   -64.6,   -24.2,     6.2,    24.0,    14.8,   -25.4,&
            -5.8,     0.0,    11.9,   -21.5,     8.5,   -21.5,    15.5,&
             8.9,   -14.9,    -2.1,     0.0,   -19.7,    13.4,    12.5,&
            -6.2,    -8.4,     8.4,     3.8,    -8.2,     4.8,     0.0,&
             1.7,     0.0,     4.0,     4.9,    -5.9,    -1.2,    -2.9,&
             0.0,    -2.2,    -7.4]
    REAL, DIMENSION(66) :: G2005 = [&
            0.00,-29554.63,-1669.05,-2337.24, 3047.69, 1657.76, 1336.30,&
        -2305.83, 1246.39,  672.51,  920.55,  797.96,  210.65, -379.86,&
          100.00, -227.00,  354.41,  208.95, -136.54, -168.05,  -13.55,&
           73.60,   69.56,   76.74, -151.34,  -14.58,   14.58,  -86.36,&
           79.88,  -74.46,   -1.65,   38.73,   12.30,    9.37,    5.42,&
            1.94,   24.80,    7.62,  -11.73,   -6.88,  -18.11,   10.17,&
            9.36,  -11.25,   -4.87,    5.58,    9.76,    3.58,   -6.94,&
            5.01,  -10.76,   -1.25,    8.76,   -6.66,   -9.22,   -2.17,&
           -6.12,    1.42,   -2.35,   -0.15,    3.06,    0.29,    2.06,&
            3.77,   -0.21,   -2.09]
    REAL, DIMENSION(66) :: H2005 = [&
            0.00,    0.00, 5077.99,    0.00,-2594.50, -515.43,    0.00,&
         -198.86,  269.72, -524.72,    0.00,  282.07, -225.23,  145.15,&
         -305.36,    0.00,   42.72,  180.25, -123.45,  -19.57,  103.85,&
            0.00,  -20.33,   54.75,   63.63,  -63.53,    0.24,   50.94,&
            0.00,  -61.14,  -22.57,    6.82,   25.35,   10.93,  -26.32,&
           -4.64,    0.00,   11.20,  -20.88,    9.83,  -19.71,   16.22,&
            7.61,  -12.76,   -0.06,    0.00,  -20.11,   12.69,   12.67,&
           -6.72,   -8.16,    8.10,    2.92,   -7.73,    6.01,    0.00,&
            2.19,    0.10,    4.46,    4.76,   -6.58,   -1.01,   -3.47,&
           -0.86,   -2.31,   -7.93]
    REAL, DIMENSION(66) :: G2010 = [&
             0.0,-29496.5, -1585.9, -2396.6,  3026.0,  1668.6,  1339.7,&
         -2326.3,  1231.7,   634.2,   912.6,   809.0,   166.6,  -357.1,&
            89.7,  -231.1,   357.2,   200.3,  -141.2,  -163.1,    -7.7,&
            72.8,    68.6,    76.0,  -141.4,   -22.9,    13.1,   -77.9,&
            80.4,   -75.0,    -4.7,    45.3,    14.0,    10.4,     1.6,&
             4.9,    24.3,     8.2,   -14.5,    -5.7,   -19.3,    11.6,&
            10.9,   -14.1,    -3.7,     5.4,     9.4,     3.4,    -5.3,&
             3.1,   -12.4,    -0.8,     8.4,    -8.4,   -10.1,    -2.0,&
            -6.3,     0.9,    -1.1,    -0.2,     2.5,    -0.3,     2.2,&
             3.1,    -1.0,    -2.8]
    REAL, DIMENSION(66) :: H2010 = [&
             0.0,     0.0,  4945.1,     0.0, -2707.7,  -575.4,     0.0,&
          -160.5,   251.7,  -536.8,     0.0,   286.4,  -211.2,   164.4,&
          -309.2,     0.0,    44.7,   188.9,  -118.1,     0.1,   100.9,&
             0.0,   -20.8,    44.2,    61.5,   -66.3,     3.1,    54.9,&
             0.0,   -57.8,   -21.2,     6.6,    24.9,     7.0,   -27.7,&
            -3.4,     0.0,    10.9,   -20.0,    11.9,   -17.4,    16.7,&
             7.1,   -10.8,     1.7,     0.0,   -20.5,    11.6,    12.8,&
            -7.2,    -7.4,     8.0,     2.2,    -6.1,     7.0,     0.0,&
             2.8,    -0.1,     4.7,     4.4,    -7.2,    -1.0,    -4.0,&
            -2.0,    -2.0,    -8.3]
    REAL, DIMENSION(45) :: DG = [&
             0.0,    11.4,    16.7,   -11.3,    -3.9,     2.7,     1.3,&
            -3.9,    -2.9,    -8.1,    -1.4,     2.0,    -8.9,     4.4,&
            -2.3,    -0.5,     0.5,    -1.5,    -0.7,     1.3,     1.4,&
            -0.3,    -0.3,    -0.3,     1.9,    -1.6,    -0.2,     1.8,&
             0.2,    -0.1,    -0.6,     1.4,     0.3,     0.1,    -0.8,&
             0.4,    -0.1,     0.1,    -0.5,     0.3,    -0.3,     0.3,&
             0.2,    -0.5,     0.2]
    REAL, DIMENSION(45) :: DH = [&
             0.0,     0.0,   -28.8,     0.0,   -23.0,   -12.9,     0.0,&
             8.6,    -2.9,    -2.1,     0.0,     0.4,     3.2,     3.6,&
            -0.8,     0.0,     0.5,     1.5,     0.9,     3.7,    -0.6,&
             0.0,    -0.1,    -2.1,    -0.4,    -0.5,     0.8,     0.5,&
             0.0,     0.6,     0.3,    -0.2,    -0.1,    -0.8,    -0.3,&
             0.2,     0.0,     0.0,     0.2,     0.5,     0.4,     0.1,&
            -0.1,     0.4,     0.4]
    
    IF(FIRST_RUN) then ! Initialize REC
        FIRST_RUN = .FALSE.

        DO N=1,11
            N2 = 2*N - 1
            N2 = N2*(N2 - 2)
            DO M=1,N
                MN = N*(N - 1)/2 + M
                REC(MN) = FLOAT( (N-M)*(N+M-2) )/FLOAT(N2)
            ENDDO
        ENDDO
    ENDIF

    ! If the year has changed, interpolate to the new year.
    IF (IY /= IYR) then
        IYR = IY
        IF (IYR < 1900) IYR=1900
        IF (IYR > 2015) IYR=2015
        IF (IY /= IYR.AND.mess) WRITE (konsol,999)IY,IYR
     
        F2 = MOD(IYR,5)/5.
        F1 = 1. - F2
        IF (IYR  <  1905) THEN      !INTERPOLATE BETWEEN 1900 - 1905
            G = G1900*F1 + G1905*F2
            H = H1900*F1 + H1905*F2
        ELSE IF (IYR  <  1910) THEN      !INTERPOLATE BETWEEN 1905 - 1910
            G = G1905*F1 + G1910*F2
            H = H1905*F1 + H1910*F2
        ELSE IF (IYR  <  1915) THEN      !INTERPOLATE BETWEEN 1910 - 1915
            G = G1910*F1 + G1915*F2
            H = H1910*F1 + H1915*F2
        ELSE IF (IYR  <  1920) THEN      !INTERPOLATE BETWEEN 1915 - 1920
            G = G1915*F1 + G1920*F2
            H = H1915*F1 + H1920*F2
        ELSE IF (IYR  <  1925) THEN      !INTERPOLATE BETWEEN 1920 - 1925
            G = G1920*F1 + G1925*F2
            H = H1920*F1 + H1925*F2
        ELSE IF (IYR  <  1930) THEN      !INTERPOLATE BETWEEN 1925 - 1930
            G = G1925*F1 + G1930*F2
            H = H1925*F1 + H1930*F2
        ELSE IF (IYR  <  1935) THEN      !INTERPOLATE BETWEEN 1930 - 1935
            G = G1930*F1 + G1935*F2
            H = H1930*F1 + H1935*F2
        ELSE IF (IYR  <  1940) THEN      !INTERPOLATE BETWEEN 1935 - 1940
            G = G1935*F1 + G1940*F2
            H = H1935*F1 + H1940*F2
        ELSE IF (IYR  <  1945) THEN      !INTERPOLATE BETWEEN 1940 - 1945
            G = G1940*F1 + G1945*F2
            H = H1940*F1 + H1945*F2
        ELSE IF (IYR  <  1950) THEN      !INTERPOLATE BETWEEN 1945 - 1950
            G = G1945*F1 + G1950*F2
            H = H1945*F1 + H1950*F2
        ELSE IF (IYR  <  1955) THEN      !INTERPOLATE BETWEEN 1950 - 1955
            G = G1950*F1 + G1955*F2
            H = H1950*F1 + H1955*F2
        ELSE IF (IYR  <  1960) THEN      !INTERPOLATE BETWEEN 1955 - 1960
            G = G1955*F1 + G1960*F2
            H = H1955*F1 + H1960*F2
        ELSE IF (IYR  <  1965) THEN      !INTERPOLATE BETWEEN 1960 - 1965
            G = G1960*F1 + G1965*F2
            H = H1960*F1 + H1965*F2
        ELSE IF (IYR  <  1970) THEN      !INTERPOLATE BETWEEN 1965 - 1970
            G = G1965*F1 + G1970*F2
            H = H1965*F1 + H1970*F2
        ELSE IF (IYR  <  1975) THEN      !INTERPOLATE BETWEEN 1970 - 1975
            G = G1970*F1 + G1975*F2
            H = H1970*F1 + H1975*F2
        ELSE IF (IYR  <  1980) THEN      !INTERPOLATE BETWEEN 1975 - 1980
            G = G1975*F1 + G1980*F2
            H = H1975*F1 + H1980*F2
        ELSE IF (IYR  <  1985) THEN      !INTERPOLATE BETWEEN 1980 - 1985
            G = G1980*F1 + G1985*F2
            H = H1980*F1 + H1985*F2
        ELSE IF (IYR  <  1990) THEN      !INTERPOLATE BETWEEN 1985 - 1990
            G = G1985*F1 + G1990*F2
            H = H1985*F1 + H1990*F2
        ELSE IF (IYR  <  1995) THEN      !INTERPOLATE BETWEEN 1990 - 1995
            G = G1990*F1 + G1995*F2
            H = H1990*F1 + H1995*F2
        ELSE IF (IYR  <  2000) THEN      !INTERPOLATE BETWEEN 1995 - 2000
            G = G1995*F1 + G2000*F2
            H = H1995*F1 + H2000*F2
        ELSE IF (IYR  <  2005) THEN      !INTERPOLATE BETWEEN 2000 - 2005
            G = G2000*F1 + G2005*F2
            H = H2000*F1 + H2005*F2
        ELSE IF (IYR  <  2010) THEN      !INTERPOLATE BETWEEN 2005 - 2010
            G = G2005*F1 + G2010*F2
            H = H2005*F1 + H2010*F2
        ELSE ! EXTRAPOLATE BEYOND 2010:
            DT=FLOAT(IYR)-2010.
            G(:45) = G2010(:45) + DG*DT
            H(:45) = H2010(:45) + DH*DT
            G(46:) = G2010(46:)
            H(46:) = H2010(46:)
        ENDIF

        !   COEFFICIENTS FOR A GIVEN YEAR HAVE BEEN CALCULATED; NOW MULTIPLY
        !   THEM BY SCHMIDT NORMALIZATION FACTORS:
        S = 1.
        DO N=2,11
            MN = N*(N-1)/2+1
            S = S*FLOAT( 2*N - 3 )/FLOAT( N - 1 )
            G(MN) = G(MN)*S
            H(MN) = H(MN)*S
            P = S
            DO M=2,N
                AA = 1.
                IF (M == 2) AA = 2.

                P = P*SQRT( AA*FLOAT(N-M+1)/FLOAT(N+M-2) )
                MNN = MN + M - 1
                G(MNN) = G(MNN)*P
                H(MNN) = H(MNN)*P
            ENDDO
        ENDDO
    ENDIF

    ! NOW CALCULATE THE FIELD COMPONENTS
    ! (IN CASE OF MULTIPLE INVOCATIONS WITH THE SAME VALUES OF IY AND NM,
    ! CALCULATIONS START RIGHT HERE):
    PP = 1./R
    P = PP
    K = NM+1
    DO N=1,K
        P = P*PP
        A(N) = P
        B(N) = P*N
    END DO

    P  = 1.
    D  = 0.
    BBR= 0.
    BBT= 0.
    BBF= 0.
    U  = T
    CF = COS(F)
    SF = SIN(F)
    C  = COS(U)
    S  = SIN(U)

    DO M=1,K
        IF(M /= 1) then
            MM = M - 1
            W  = X
            X  = W*CF + Y*SF
            Y  = Y*CF - W*SF
        else
            X = 0.
            Y = 1.
        endif

        Q  = P
        Z  = D
        BI = 0.
        P2 = 0.
        D2 = 0.
        DO N = M,K
            AN = A(N)
            MN = N*(N - 1)/2 + M
            E  = G(MN)
            HH = H(MN)
            W  = E*Y + HH*X
            BBR = BBR + B(N)*W*Q
            BBT = BBT - AN*W*Z
            IF(M /=  1) then
                QQ = Q
                IF(S < 1.E-5) QQ = Z
                BI = BI + AN*(E*X - HH*Y)*QQ
            endif
            XK = REC(MN)
            DP = C*Z - S*Q - XK*D2
            PM = C*Q - XK*P2
            D2 = Z
            P2 = Q
            Z = DP
            Q = PM
        ENDDO

        D = S*D + C*P
        P = S*P
        IF(M /=  1) then
            BI = BI*MM
            BBF = BBF + BI
        endif
    END DO

    BR = BBR
    BT = BBT
    IF(S >= 1.E-5) then
        BF = BBF/S
        RETURN
    endif

    IF(C < 0.) BBF = -BBF
    BF = BBF
    RETURN

! FORMAT statemetns
999 FORMAT(/'   IGRF: GIVEN YEAR',I5,' IS OUT OF INTERVAL 1900-2015'/,&
            '   *** CALCULATIONS WILL BE DONE FOR YEAR =',I5,' ***'/)
END SUBROUTINE IGRF


!  *********************************************************************
!  If only IYR is given then CALL RECALC(IYR,0,25,0,0)
!  THIS IS A MODIFIED VERSION OF THE SUBROUTINE RECOMP WRITTEN BY
!  N. A. TSYGANENKO. SINCE I WANT TO USE IT IN PLACE OF SUBROUTINE
!  RECALC, I HAVE RENAMED THIS ROUTINE RECALC AND ELIMINATED THE
!  ORIGINAL RECALC FROM THIS VERSION OF THE <GEOPACK.FOR> PACKAGE.
!  THIS WAY ALL ORIGINAL CALLS TO RECALC WILL CONTINUE TO WORK WITHOUT
!  HAVING TO CHANGE THEM TO CALLS TO RECOMP.
!
!  AN ALTERNATIVE VERSION OF THE SUBROUTINE RECALC FROM THE GEOPACK
!  PACKAGE BASED ON A DIFFERENT APPROACH TO DERIVATION OF ROTATION
!  MATRIX ELEMENTS
!
!  THIS SUBROUTINE WORKS BY 20% FASTER THAN RECALC AND IS EASIER TO
!  UNDERSTAND
!  #####################################################
!  #  WRITTEN BY  N.A. TSYGANENKO ON DECEMBER 1, 1991  #
!  #####################################################
!  Modified by Mauricio Peredo, Hughes STX at NASA/GSFC Code 695,
!  September 1992
!
!  Modified to accept years up to year 2000 and updated IGRF coeficients
!     from 1945 (updated by V. Papitashvili, February 1995)
!
!  Modified to accept years up to 2005 (V. Papitashvili, January 2001)
!
!  Modified to accept years from 1900 through 2010 using the DGRF &
!     IGRF-10 coeficients (updated by V. Papitashvili, November 2005)
!
!  Modified to accept years up to 2015 (V. Papitashvili, January 2011)
!
!  Modified to accept years up to 2020 (D. Bilitza, October 2015)
!
!   OTHER SUBROUTINES CALLED BY THIS ONE: SUN
!
!     IYR = YEAR NUMBER (FOUR DIGITS)
!     IDAY = DAY OF YEAR (DAY 1 = JAN 1)
!     IHOUR = HOUR OF DAY (00 TO 23)
!     MIN = MINUTE OF HOUR (00 TO 59)
!     ISEC = SECONDS OF DAY(00 TO 59)
!  *********************************************************************
SUBROUTINE RECALC(IYR,IDAY,IHOUR,MIN,ISEC)

    IMPLICIT NONE

    REAL ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,CPS,&
        SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,&
        A33,DS3,F2,F1,G10,G11,H11,DT,SQ,SQQ,SQR,S1,S2,&
        S3,CGST,SGST,DIP1,DIP2,DIP3,Y1,Y2,Y3,Y,Z1,Z2,Z3,DJ,&
        T,OBLIQ,DZ1,DZ2,DZ3,DY1,DY2,DY3,EXMAGX,EXMAGY,EXMAGZ,&
        EYMAGX,EYMAGY,GST,SLONG,SRASN,SDEC,BA(8),DECARG, PI, UMR

    INTEGER IYR,IDAY,IHOUR,MIN,ISEC,K,IY,IDE,IYE,konsol

    logical mess

    COMMON/C1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,CPS,&
        SHI,CHI,HI,PSI,XMUT,A11,A21,A31,A12,A22,A32,A13,A23,A33,DS3,&
        K,IY,BA
    common/iounit/konsol,mess
    COMMON /CONST/UMR,PI

    DATA IYE,IDE/2*0/

    !  IYE AND IDE ARE THE CURRENT VALUES OF YEAR AND DAY NUMBER
    IF (IYR /= IYE .OR. IDAY /= IDE) then
        IY  = IYR
        IDE = IDAY
        IF(IY < 1900) IY=1900
        IF(IY > 2020) IY=2020
        !  WE ARE RESTRICTED BY THE INTERVAL 1900-2020, FOR WHICH THE DGRF & IGRF-11
        !  COEFFICIENTS ARE KNOWN; IF IYR IS OUTSIDE THIS INTERVAL, THE
        !  SUBROUTINE GIVES A WARNING (BUT DOES NOT REPEAT IT AT THE NEXT CALLS)
        IF(IY /= IYR.AND.mess) write(konsol,10) IYR,IY

        IYE = IY
        !  LINEAR INTERPOLATION OF THE GEODIPOLE MOMENT COMPONENTS BETWEEN THE
        !  VALUES FOR THE NEAREST EPOCHS:
        F2 = MOD(FLOAT(IY) + FLOAT(IDAY)/365., 5.)/5.
        F1 = 1.D0 - F2
        IF (IY < 1905) THEN                             !1900 - 1905
            G10 = 31543.*F1 + 31464.*F2
            G11 = -2298.*F1 - 2298.*F2
            H11 =  5922.*F1 + 5909.*F2
        ELSEIF (IY < 1910) THEN                         !1905 - 1910
            G10 = 31464.*F1 + 31354.*F2
            G11 = -2298.*F1 - 2297.*F2
            H11 =  5909.*F1 + 5898.*F2
        ELSEIF (IY < 1915) THEN                         !1910 - 1915
            G10 = 31354.*F1 + 31212.*F2
            G11 = -2297.*F1 - 2306.*F2
            H11 =  5898.*F1 + 5875.*F2
        ELSEIF (IY < 1920) THEN                         !1915 - 1920
            G10 = 31212.*F1 + 31060.*F2
            G11 = -2306.*F1 - 2317.*F2
            H11 = 5875.*F1 + 5845.*F2
        ELSEIF (IY < 1925) THEN                         !1920 - 1925
            G10 = 31060.*F1 + 30926.*F2
            G11 = -2317.*F1 - 2318.*F2
            H11 =  5845.*F1 + 5817.*F2
        ELSEIF (IY < 1930) THEN                         !1925 - 1930
            G10 = 30926.*F1 + 30805.*F2
            G11 = -2318.*F1 - 2316.*F2
            H11 =  5817.*F1 + 5808.*F2
        ELSEIF (IY < 1935) THEN                        !1930 - 1935
            G10 = 30805.*F1 + 30715.*F2
            G11 = -2316.*F1 - 2306.*F2
            H11 =  5808.*F1 + 5812.*F2
        ELSEIF (IY < 1940) THEN                        !1935 - 1940
            G10 = 30715.*F1 + 30654.*F2
            G11 = -2306.*F1 - 2292.*F2
            H11 =  5812.*F1 + 5821.*F2
        ELSEIF (IY < 1945) THEN                        !1940 - 1945
            G10 = 30654.*F1 + 30594.*F2
            G11 = -2292.*F1 - 2285.*F2
            H11 =  5821.*F1 + 5810.*F2
        ELSEIF (IY < 1950) THEN                        !1945 - 1950
            G10 = 30594.*F1 + 30554.*F2
            G11 = -2285.*F1 - 2250.*F2
            H11 =  5810.*F1 + 5815.*F2
        ELSEIF (IY < 1955) THEN                        !1950 - 1955
            G10 = 30554.*F1 + 30500.*F2
            G11 = -2250.*F1 - 2215.*F2
            H11 =  5815.*F1 + 5820.*F2
        ELSEIF (IY < 1960) THEN                        !1955 - 1960
            G10 = 30500.*F1 + 30421.*F2
            G11 = -2215.*F1 - 2169.*F2
            H11 =  5820.*F1 + 5791.*F2
        ELSEIF (IY < 1965) THEN                        !1960 - 1965
            G10 = 30421.*F1 + 30334.*F2
            G11 = -2169.*F1 - 2119.*F2
            H11 =  5791.*F1 + 5776.*F2
        ELSEIF (IY < 1970) THEN                        !1965 - 1970
            G10 = 30334.*F1 + 30220.*F2
            G11 = -2119.*F1 - 2068.*F2
            H11 =  5776.*F1 + 5737.*F2
        ELSEIF (IY < 1975) THEN                        !1970 - 1975
            G10 = 30220.*F1 + 30100.*F2
            G11 = -2068.*F1 - 2013.*F2
            H11 =  5737.*F1 + 5675.*F2
        ELSEIF (IY < 1980) THEN                        !1975 - 1980
            G10 = 30100.*F1 + 29992.*F2
            G11 = -2013.*F1 - 1956.*F2
            H11 =  5675.*F1 + 5604.*F2
        ELSEIF (IY < 1985) THEN                        !1980 - 1985
            G10 = 29992.*F1 + 29873.*F2
            G11 = -1956.*F1 - 1905.*F2
            H11 =  5604.*F1 + 5500.*F2
        ELSEIF (IY < 1990) THEN                        !1985 - 1990
            G10 = 29873.*F1 + 29775.*F2
            G11 = -1905.*F1 - 1848.*F2
            H11 =  5500.*F1 + 5406.*F2
        ELSEIF (IY < 1995) THEN                        !1990 - 1995
            G10 = 29775.*F1 + 29692.*F2
            G11 = -1848.*F1 - 1784.*F2
            H11 =  5406.*F1 + 5306.*F2
        ELSEIF (IY < 2000) THEN                        !1995 - 2000
            G10 = 29692.*F1 + 29619.4*F2
            G11 = -1784.*F1 - 1728.2*F2
            H11 =  5306.*F1 + 5186.1*F2
        ELSEIF (IY < 2005) THEN                        !2000 - 2005
            G10 = 29619.4*F1 + 29554.63*F2
            G11 = -1728.2*F1 - 1669.05*F2
            H11 =  5186.1*F1 + 5077.99*F2
        ELSEIF (IY < 2010) THEN                        !2005 - 2010
            G10 = 29554.63*F1 + 29496.57*F2
            G11 = -1669.05*F1 - 1586.42*F2
            H11 =  5077.99*F1 + 4944.26*F2
        ELSEIF (IY < 2015) THEN                        !2010 - 2015
            G10 = 29496.57*F1 + 29442.0*F2
            G11 = -1586.42*F1 - 1501.0*F2
            H11 =  4944.26*F1 + 4797.1*F2
        ELSE                                            !2015 - 2020
            DT = FLOAT(IY) + FLOAT(IDAY)/365. - 2015.
            G10 = 29442.0 - 10.3*DT
            G11 = -1501.0 + 18.1*DT
            H11 =  4797.1 - 26.6*DT
        ENDIF

        !  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EzMAG IN GEO COORD
        !  SYSTEM:
        !  SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
        !         ST0 * CL0                ST0 * SL0                CT0
        SQ   = G11**2 + H11**2
        SQQ  = SQRT(SQ)
        SQR  = SQRT(G10**2 + SQ)
        SL0  = - H11/SQQ
        CL0  = - G11/SQQ
        ST0  = SQQ/SQR
        CT0  = G10/SQR
        STCL = ST0*CL0
        STSL = ST0*SL0
        CTSL = CT0*SL0
        CTCL = CT0*CL0
    endif

    !  THE CALCULATIONS ARE TERMINATED IF ONLY GEO-MAG TRANSFORMATION
    !  IS TO BE DONE  (IHOUR>24 IS THE AGREED CONDITION FOR THIS CASE):
    IF (IHOUR > 24) RETURN

    CALL SUN(IY,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
    !  S1,S2, AND S3 ARE THE COMPONENTS OF THE UNIT VECTOR EXGSM=EXGSE
    !  IN THE SYSTEM GEI POINTING FROM THE EARTH'S CENTER TO THE SUN:
    S1   = COS(SRASN)*COS(SDEC)
    S2   = SIN(SRASN)*COS(SDEC)
    S3   = SIN(SDEC)
    CGST = COS(GST)
    SGST = SIN(GST)

    !  DIP1, DIP2, AND DIP3 ARE THE COMPONENTS OF THE UNIT VECTOR
    !  EZSM = EZMAG IN THE SYSTEM GEI:
    DIP1 = STCL*CGST - STSL*SGST
    DIP2 = STCL*SGST + STSL*CGST
    DIP3 = CT0

    !  NOW CALCULATE THE COMPONENTS OF THE UNIT VECTOR EYGSM IN THE SYSTEM
    !  GEI BY TAKING THE VECTOR PRODUCT D x S AND NORMALIZING IT TO UNIT
    !  LENGTH:
    Y1 = DIP2*S3 - DIP3*S2
    Y2 = DIP3*S1 - DIP1*S3
    Y3 = DIP1*S2 - DIP2*S1
    Y  = SQRT(Y1**2 + Y2**2 + Y3**2)
    Y1 = Y1/Y
    Y2 = Y2/Y
    Y3 = Y3/Y

    !  THEN IN THE GEI SYSTEM THE UNIT VECTOR Z = EZGSM = EXGSM x EYGSM = S x Y
    !  HAS THE COMPONENTS:
    Z1 = S2*Y3 - S3*Y2
    Z2 = S3*Y1 - S1*Y3
    Z3 = S1*Y2 - S2*Y1

    !  THE VECTOR EZGSE (HERE DZ) IN GEI HAS THE COMPONENTS (0, - SIN(DELTA),
    !  COS(DELTA))  =  (0., - 0.397823,0.917462); HERE DELTA  =  23.44214 DEG FOR
    !  THE EPOCH 1978 (SEE THE BOOK BY GUREVICH OR OTHER ASTRONOMICAL
    !  HANDBOOKS). HERE THE MOST ACCURATE TIME - DEPENDENT FORMULA IS USED:
    DJ = FLOAT(365*(IY - 1900) + (IY - 1901)/4  + IDAY) - 0.5 + FLOAT(ISEC)/86400.
    T  = DJ/36525.
    OBLIQ = (23.45229 - 0.0130125*T)*PI/180
    DZ1 = 0.
    DZ2 =  - SIN(OBLIQ)
    DZ3 = COS(OBLIQ)
    !  THEN THE UNIT VECTOR EYGSE IN GEI SYSTEM IS THE VECTOR PRODUCT DZ x S
    DY1 = DZ2*S3 - DZ3*S2
    DY2 = DZ3*S1 - DZ1*S3
    DY3 = DZ1*S2 - DZ2*S1

    !  THE ELEMENTS OF THE MATRIX GSE TO GSM ARE THE SCALAR PRODUCTS:
    !  CHI = EM22 = (EYGSM,EYGSE), SHI = EM23 = (EYGSM,EZGSE),
    !  EM32 = (EZGSM,EYGSE) =  - EM23, AND EM33 = (EZGSM,EZGSE) = EM22
    CHI = Y1*DY1 + Y2*DY2 + Y3*DY3
    SHI = Y1*DZ1 + Y2*DZ2 + Y3*DZ3
    DECARG = SHI
    IF(ABS(DECARG) > 1.) DECARG = SIGN(1.,DECARG)
    HI = ASIN(DECARG)
    !  TILT ANGLE: PSI = ARCSIN(DIP,EXGSM)
    SPS = DIP1*S1 + DIP2*S2 + DIP3*S3
    CPS = SQRT(1. - SPS**2)
    DECARG = SPS
    IF(ABS(DECARG) > 1.) DECARG = SIGN(1.,DECARG)
    PSI = ASIN(DECARG)

    !  THE ELEMENTS OF THE MATRIX MAG TO SM ARE THE SCALAR PRODUCTS:
    !  CFI = GM22 = (EYSM,EYMAG), SFI = GM23 = (EYSM,EXMAG); THEY CAN BE DERIVED
    !  AS FOLLOWS:
    !  IN GEO THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS
    !  (CT0*CL0,CT0*SL0, - ST0) AND ( - SL0,CL0,0), RESPECTIVELY. HENCE, IN
    !  GEI SYSTEM THE COMPONENTS ARE:
    !  EXMAG:    CT0*CL0*COS(GST) - CT0*SL0*SIN(GST)
    !            CT0*CL0*SIN(GST) + CT0*SL0*COS(GST)
    !             - ST0
    !  EYMAG:     - SL0*COS(GST) - CL0*SIN(GST)
    !             - SL0*SIN(GST) + CL0*COS(GST)
    !             0
    !  THE COMPONENTS OF EYSM IN GEI WERE FOUND ABOVE AS Y1, Y2, AND Y3;
    !  NOW WE ONLY HAVE TO COMBINE THE QUANTITIES INTO SCALAR PRODUCTS:
    EXMAGX = CT0*(CL0*CGST - SL0*SGST)
    EXMAGY = CT0*(CL0*SGST + SL0*CGST)
    EXMAGZ =-ST0
    EYMAGX =-(SL0*CGST + CL0*SGST)
    EYMAGY =-(SL0*SGST - CL0*CGST)
    CFI    = Y1*EYMAGX + Y2*EYMAGY
    SFI    = Y1*EXMAGX + Y2*EXMAGY + Y3*EXMAGZ
    XMUT   = (ATAN2(SFI,CFI) + PI)*12./PI

    !  THE ELEMENTS OF THE MATRIX GEO TO GSM ARE THE SCALAR PRODUCTS:
    !  A11 = (EXGEO,EXGSM), A12 = (EYGEO,EXGSM), A13 = (EZGEO,EXGSM),
    !  A21 = (EXGEO,EYGSM), A22 = (EYGEO,EYGSM), A23 = (EZGEO,EYGSM),
    !  A31 = (EXGEO,EZGSM), A32 = (EYGEO,EZGSM), A33 = (EZGEO,EZGSM),
    !  ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:
    !  EXGEO = (CGST,SGST,0), EYGEO = ( - SGST,CGST,0), EZGEO = (0,0,1)
    !  EXGSM = (S1,S2,S3),  EYGSM = (Y1,Y2,Y3),   EZGSM = (Z1,Z2,Z3)
    !  AND  THEREFORE:
    A11 = S1*CGST + S2*SGST
    A12 =-S1*SGST + S2*CGST
    A13 = S3
    A21 = Y1*CGST + Y2*SGST
    A22 =-Y1*SGST + Y2*CGST
    A23 = Y3
    A31 = Z1*CGST + Z2*SGST
    A32 =-Z1*SGST + Z2*CGST
    A33 = Z3
    RETURN

! FORMAT statemetns
10  FORMAT(/' RECALC: GIVEN YEAR',I5,' IS OUT OF INTERVAL 1900-2010'/,&
            '   *** CALCULATIONS WILL BE DONE FOR YEAR =',I5,' ***'/)
END SUBROUTINE RECALC


!  *********************************************************************
!   CONVERTS SPHERICAL COORDS INTO CARTESIAN ONES AND VICA VERSA
!    (THETA AND PHI IN RADIANS).
!                  J>0            J<0
!-----INPUT:   J,R,THETA,PHI     J,X,Y,Z
!----OUTPUT:      X,Y,Z        R,THETA,PHI
!  AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
!      STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
!      (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
!  *********************************************************************
SUBROUTINE SPHCAR(R,THETA,PHI,X,Y,Z,J)

    IMPLICIT NONE

    REAL R,THETA,PHI,X,Y,Z,SQ, PI, UMR
    INTEGER J

    COMMON /CONST/UMR, PI

    IF(J > 0) then
        SQ = R*SIN(THETA)
        X  = SQ*COS(PHI)
        Y  = SQ*SIN(PHI)
        Z  = R*COS(THETA)
        RETURN
    endif

    SQ = X**2 + Y**2
    R  = SQRT( SQ + Z**2 )
    IF (SQ /= 0.) then
        SQ  = SQRT(SQ)
        PHI = ATAN2(Y,X)
        THETA=ATAN2(SQ,Z)
        IF (PHI < 0.) PHI = PHI + 2.*PI
    else
        PHI = 0.
        IF (Z < 0.) THEN
            THETA = PI
        ELSE
            THETA=0.
        ENDIF
    ENDIF
    RETURN
END SUBROUTINE SPHCAR


!  *********************************************************************
!   CALCULATES CARTESIAN FIELD COMPONENTS FROM SPHERICAL ONES
!-----INPUT:   TETA,PHI - SPHERICAL ANGLES OF THE POINT IN RADIANS
!              BR,BTET,BPHI -  SPHERICAL COMPONENTS OF THE FIELD
!-----OUTPUT:  BX,BY,BZ - CARTESIAN COMPONENTS OF THE FIELD
!  AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
!      STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
!      (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
!  *********************************************************************
SUBROUTINE BSPCAR(THETA,PHI,BR,BTET,BPHI,BX,BY,BZ)
    IMPLICIT NONE

    REAL THETA,PHI,BR,BTET,BPHI,BX,BY,BZ,S,C,SF,CF,BE

    S  = SIN(THETA)
    C  = COS(THETA)
    SF = SIN(PHI)
    CF = COS(PHI)

    BE = BR*S  + BTET*C
    BX = BE*CF - BPHI*SF
    BY = BE*SF + BPHI*CF
    BZ = BR*C  - BTET*S
    RETURN
END SUBROUTINE BSPCAR


!  *********************************************************************
! CONVERTS GEOCENTRIC (GEO) TO DIPOLE (MAG) COORDINATES OR VICA VERSA.
! IYR IS YEAR NUMBER (FOUR DIGITS).
!                           J>0                J<0
!-----INPUT:  J,XGEO,YGEO,ZGEO,IYR   J,XMAG,YMAG,ZMAG,IYR
!-----OUTPUT:    XMAG,YMAG,ZMAG        XGEO,YGEO,ZGEO
!  AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
!      STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
!      (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
!  *********************************************************************
SUBROUTINE GEOMAG(XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,J,IYR)

    IMPLICIT NONE

    REAL XGEO,YGEO,ZGEO,XMAG,YMAG,ZMAG,ST0,CT0,SL0,CL0,CTCL,&
        STCL,CTSL,STSL,AB(19),BB(8)
    INTEGER J,IYR,K,IY,II

    COMMON/C1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,AB,K,IY,BB
    DATA II/1/

    IF(IYR /= II) then
        II=IYR
        CALL RECALC(II,0,25,0,0)
    ENDIF

    IF(J > 0) then
        XMAG = XGEO*CTCL + YGEO*CTSL - ZGEO*ST0
        YMAG = YGEO*CL0  - XGEO*SL0
        ZMAG = XGEO*STCL + YGEO*STSL + ZGEO*CT0
    else
        XGEO = XMAG*CTCL - YMAG*SL0 + ZMAG*STCL
        YGEO = XMAG*CTSL + YMAG*CL0 + ZMAG*STSL
        ZGEO = ZMAG*CT0  - XMAG*ST0
    endif
    RETURN
END SUBROUTINE GEOMAG


!  *********************************************************************
! CONVERTS DIPOLE (MAG) TO SOLAR MAGNETIC (SM) COORDINATES OR VICA VERSA
!                    J>0              J<0
!-----INPUT: J,XMAG,YMAG,ZMAG     J,XSM,YSM,ZSM
!----OUTPUT:    XSM,YSM,ZSM       XMAG,YMAG,ZMAG
!  ATTENTION: SUBROUTINE RECALC MUST BE CALLED BEFORE MAGSM IN TWO CASES
!     /A/  BEFORE THE FIRST USE OF MAGSM
!     /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC ARE
!          DIFFERENT FROM THOSE IN THE PRECEDING CALL OF  MAGSM
!  AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
!      STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
!      (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
!  *********************************************************************
SUBROUTINE MAGSM(XMAG,YMAG,ZMAG,XSM,YSM,ZSM,J)

    IMPLICIT NONE

    REAL XMAG,YMAG,ZMAG,XSM,YSM,ZSM,SFI,CFI,A(8),B(7),&
        AB(10),BA(8)
    INTEGER J,K,IY

    COMMON/C1/ A,SFI,CFI,B,AB,K,IY,BA

    IF (J > 0) then
        XSM = XMAG*CFI - YMAG*SFI
        YSM = XMAG*SFI + YMAG*CFI
        ZSM = ZMAG
    else
        XMAG = XSM*CFI + YSM*SFI
        YMAG = YSM*CFI - XSM*SFI
        ZMAG = ZSM
    endif
    RETURN
END SUBROUTINE MAGSM


!  *********************************************************************
! CONVERTS SOLAR MAGNETIC (SM) TO SOLAR MAGNETOSPHERIC (GSM) COORDINATES
!   OR VICA VERSA.
!                  J>0                 J<0
!-----INPUT: J,XSM,YSM,ZSM        J,XGSM,YGSM,ZGSM
!----OUTPUT:  XGSM,YGSM,ZGSM       XSM,YSM,ZSM
!  ATTENTION: SUBROUTINE RECALC MUST BE CALLED BEFORE SMGSM IN TWO CASES
!     /A/  BEFORE THE FIRST USE OF SMGSM
!     /B/  IF THE CURRENT VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC ARE
!          DIFFERENT FROM THOSE IN THE PRECEDING CALL OF SMGSM
!  AUTHOR: NIKOLAI A. TSYGANENKO, INSTITUTE OF PHYSICS, ST.-PETERSBURG
!      STATE UNIVERSITY, STARY PETERGOF 198904, ST.-PETERSBURG, RUSSIA
!      (now the NASA Goddard Space Fligth Center, Greenbelt, Maryland)
!  *********************************************************************
SUBROUTINE SMGSM(XSM,YSM,ZSM,XGSM,YGSM,ZGSM,J)

    IMPLICIT NONE

    REAL XSM,YSM,ZSM,XGSM,YGSM,ZGSM,SPS,CPS,A(10),B(15),AB(8)
    INTEGER J,K,IY
    COMMON/C1/ A,SPS,CPS,B,K,IY,AB

    IF (J > 0) then
        XGSM = XSM*CPS + ZSM*SPS
        YGSM = YSM
        ZGSM = ZSM*CPS - XSM*SPS
    else
        XSM = XGSM*CPS - ZGSM*SPS
        YSM = YGSM
        ZSM = XGSM*SPS + ZGSM*CPS
    endif
    RETURN
END SUBROUTINE SMGSM


!--------------------------------------------------------------------
!      calculates magnetic local time
!      Inputs:
!             IYYYY..Year as YYYY, e.g. 1998
!             DDD..day of year (1.1. = 0)
!             UTHR..universal time in decimal hours
!             GLAT,GLON..latitude north and longitude east in degrees
!      Output:
!             MLT..magnetic local time in decimal hours
!      Required subroutines: DPMTRX
!--------------------------------------------------------------------
SUBROUTINE CLCMLT(IYYYY,DDD,UTHR,GLAT,GLON,MLT)

    INTEGER IYYYY,DDD
    REAL UTHR,GLAT,GLON,MLT
    REAL DTOR,PI,XG(3)
    REAL XXM(3),YYM(3),ZZM(3)
    INTEGER IHOUR,MIN,ISEC
    REAL GST,SLONG,SRASN,SDEC
    REAL BE,CAL,SA(3),S,C,SG(3),SM(3)
    REAL LAM,LAMS,DELLAM

    COMMON /CONST/DTOR,PI

    XG(1) = COS(GLAT*DTOR)*COS(GLON*DTOR)
    XG(2) = COS(GLAT*DTOR)*SIN(GLON*DTOR)
    XG(2) = SIN(GLAT*DTOR)

    CALL DPMTRX(IYYYY,DDD,XXM,YYM,ZZM)
    !       transform
    XM = dot_product(XXM, XG)
    YM = dot_product(YYM, XG)
    ZM = dot_product(ZZM, XG)

    IHOUR = INT(UTHR)
    MIN   = INT( (UTHR - IHOUR)*60 )
    ISEC  = INT( (UTHR - IHOUR - MIN/60.0)*3600 )
    CALL SUN(IYYYY,DDD+1,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
    BE  = GST
    CAL = COS(SRASN)
    SA(3) = SIN(SDEC)
    SA(1) = COS(SDEC)
    SA(2) = SA(1)*SIN(SRASN)
    SA(1) = SA(1)*CAL
    S = SIN(BE)
    C = COS(BE)
    SG(1) = C*SA(1) + S*SA(2)
    SG(2) = C*SA(2) - S*SA(1)
    SG(3) = SA(3)

    !       transform
    SM(1) = dot_product(XXM, SG)
    SM(2) = dot_product(YYM, SG)
    SM(3) = dot_product(ZZM, SG)

    LAM  = ATAN2(YM, XM)
    LAMS = ATAN2(SM(2), SM(1))
    DELLAM = LAM - LAMS
    IF (DELLAM < 0.) DELLAM = DELLAM + 2*PI
    MLT = MOD(DELLAM/PI*12. + 12., 24.)
    RETURN
END SUBROUTINE CLCMLT


!--------------------------------------------------------------------------
!      calculates othonormal matrix (columns XM,YM,ZM) for transformation
!      from geographic to magnetic coordinates
!      Inputs:
!             IYYYY..year
!               DDD..day of year (1.1 = 0)
!      Outputs:
!               XM,YM,ZM..colums of the matrix
!      Notes:
!      MX(N),MY(N),MZ(N)..coordinates of the B vector in geographic system
!                for years stored in YR(N)
!      N..number of elements of arrays MX,MY,MZ and YR
!--------------------------------------------------------------------------
SUBROUTINE DPMTRX(IYYYY,DDD,XM,YM,ZM)

    INTEGER IYYYY,DDD
    REAL XM(3),YM(3),ZM(3)
    REAL YR(10),MX(10),MY(10),MZ(10)
    REAL INTERP,YEAR
    REAL M,MXI,MYI,MZI,ZM12
    INTEGER N

    COMMON /DIPOL/ GHI1,GHI2,GHI3
    DATA N/10/

    ! IGRF coefficients (dipole) calculated in FELDCOF in IGRF.FOR
    MXI = -GHI2
    MYI = -GHI3
    MZI = -GHI1

    ! normalization of the vector of the dipole exis of the magnetic field
    M   = SQRT(MXI**2 + MYI**2 + MZI**2)
    MYZ = SQRT(MYI**2 + MZI**2)

    ZM(1) = MXI/M
    ZM(2) = MYI/M
    ZM(3) = MZI/M

    ZM12  = SQRT(ZM(1)**2 + ZM(2)**2)

    YM(1) =-ZM(2)/ZM12
    YM(2) = ZM(1)/ZM12
    YM(3) = 0.

    XM(1) = YM(2)*ZM(3) - YM(3)*ZM(2)
    XM(2) = YM(3)*ZM(1) - YM(1)*ZM(3)
    XM(3) = YM(1)*ZM(2) - YM(2)*ZM(1)

    RETURN
END SUBROUTINE DPMTRX
