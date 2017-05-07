!****************************************************************************************
! subroutines for IDC model
!
! includes: main subroutine CHEMION and the following subroutines and functions
!           KEMPPRN.FOR: CN2D, CNO, CN4S, CN2PLS, CNOP, CO2P, COP4S, COP2D, COP2P,
!                        CNPLS, CN2A, CN2P, CNOPV
!           RATES.FOR:   RATS
!           PESIMP.FOR:  SECIPRD, FLXCAL, FACFLX, SIGEXS, TXSION, OXRAT, T_XS_N2,
!                        T_XS_OX, OXSIGS
!           RSPRIM.FOR:  PRIMPR, SCOLUM, PARAMS, PROBS, PROBN2, YLDISS, PROBO2,
!                        SCHUMN, FACEUV, FACSR,
!
! turn on printout of intermediate quantities with JPRINT=1 also in PARAMS, PROBS,
! PROBN2, YLDISS, and PROBO2 with ISW=1.
!
! Richards, P. G., D. Bilitza, and D. Voglozin (2010), Ion density calculator (IDC):
!     A new efficient model of ionospheric ion densities, Radio Sci., 45, RS5007,
!   doi:10.1029/2009RS004332.
!
!****************************************************************************************

!-------------------------------------------------------------------------------------
!
! Subroutine ChemIon
!
! Inputs:
!     JPRINT             -- Turn file output on or off&
!     ALT                --  Altitude(km)&
!     F107, F107A        --  Solar activity indices&
!     TE, TI, TN         --  Electron, ion and neutral temperatures&
!     OXN, O2N, N2N, HEN --  O, O2, N2, and He densities (cm-3)&
!     USER_NO            --  User specified NO density (cm-3)&
!     N4S                --  N4S should be 0.5*MSIS N density (cm-3)&
!     NE                 --  electron density (cm-3)&
!     USER_NPLUS         --  User specified N+ density (cm-3)&
!     SZAD               --  solar zenith angle(D)&
! Outputs:
!     OXPLUS, O2PLUS, NOPLUS --  O+, O2+, NO+ densities (cm-3)&
!     N2PLUS, NPLUS          --  N2+ and N+ densities (cm-3)&
!     NNO, N2D               --  NO and N(2D) density (cm-3)&
!     INEWT                  --  Newton procedure fails if INEWT=0
!
!    This routine was written by Phil Richards April 2010.
!    It takes the specified input electron density and returns O+, O2+, NO+,
!    N2+, N+, NO, and N(2D) densities.These densities generally agree well
!    with the FLIP model densities.
!    All the densities except O+ are calculated from chemical equilibrium.
!    O+ is calculated using a Newton iterative procedure so that the total
!    ion density matches the input electron density.
!    N+ and NO densities can either be user specified or calculated by the model.
!    N+ generally agrees well with AE-C data and the FLIP model during the day,
!    but is inaccurate at night due to lack of diffusion.
!    NO will be very good except below about 130 km where it will be
!    underestimated due to neglect of diffusion. There is an artificial
!    floor on the NO density to prevent it from getting too low below 130 km.
!    If the Newton procedure fails to converge, all ions including O+ are
!    calculated from chemical equilibrium and then normalized to reproduce
!    the input electron density. This generally works well.
!    The Newton procedure usually works if the total calculated molecular ion
!    densities do not approach the input electron density. Difficulties are most
!    likely to happen below ~150 km and especially at night. A Newton solution is
!    usually found when the O+ density is not too much smaller than the modeled
!    molecular ion density.
!    The EUVAC model is used for solar EUV irradiances
!-------------------------------------------------------------------------------------
subroutine ChemIon( jPrint, alt, F107, F107A, TE, TI, TN,&
                    OXN, O2N, N2N, HEN, USER_NO, N4S, NE,&
                    USER_NPLUS, SZAD, OXPLUS, O2PLUS,    &
                    NOPLUS, N2PLUS, NPLUS, NNO, N2D, INEWT) 

    implicit none

    INTEGER JPRINT        !.. Turns on printing of production and loss
    INTEGER INEWT         !.. Signifies when the Newton procedure fails
    INTEGER I,J,K,ITERS   !.. loop control variables
    INTEGER IRATS         !.. Switch for different rates
    INTEGER ITS,JITER     !.. Variables for Newton procedure
    REAL TE,TN,TI         !.. Electron and ion temperatures
    !.. Geophysical parameters
    REAL F107,F107A,ALT,SZAD
    !.. Measured H+, He+, O+, N2+, NO+, O2+, N+, RPA ion density
    REAL HEPLUS,OXPLUS,N2PLUS,NOPLUS,O2PLUS,NPLUS,USER_NPLUS
    !.. O2,O,N2,NO,He,N4S, user specified NO
    REAL O2N,OXN,N2N,NNO,HEN,N4S,USER_NO
    !.. Ne, N(2P),N(2D),O+(2P),O+(2D) densities
    REAL NE,N2P,N2D,OP2D,OP2P
    !.. Total (photon & photoel) production rates O+(4S),O+(2P),O+(2D),O2+
    REAL TPROD1,PDISOP,TPROD2,TPROD3,TPROD5
    !.. Total Production rates from all sources for NO+, O2+,
    REAL TPNOP,O2PPROD
    !.. Production rates hv(e*)+N2->N+, hv+N->N+, Lyman-a -> NO+
    REAL DISNP,PHOTN,PLYNOP
    REAL PSEC                     !.. generic PE production
    REAL RTS(99)                  !.. Reaction rates array
    REAL SECPN2PLUS,EUVN2PLUS     !.. N2+ total production
    REAL H,DEX,FEX(2)             !.. used in Newton solver
    REAL SUMIONS                  !.. Sum of the major ions
    REAL PNO,LNO,PDNOSR           !.. Production and loss of NO
    REAL N2A                      !.. N2(A) density
    REAL VCON                     !.. FLIP N2(v) factor. Not used here
    REAL DISN2D,UVDISN,PN2D,LN2D  !.. Production and loss of N(2D)
    REAL ALTCHEM                  !.. altitude for chemistry calculation
    REAL N2APRD                   !.. PE production rate of N2(A)
    REAL PN4S,LN4S,DISN4S
    REAL OXPLUSAVE
    LOGICAL :: NEWT_SUCCESS
    !.. various ionization and excitation rates by EUV and PE
    REAL EUVION,PEXCIT,PEPION,OTHPR1,OTHPR2
    COMMON/EUVPRD/EUVION(3,12),PEXCIT(3,12),PEPION(3,12),OTHPR1(6)&
        ,OTHPR2(6)
    !.. initialize parameters
    DATA VCON/1.0/K/0/
    DATA PNO,LNO,PDNOSR,PLYNOP,N2A/5*0.0/
    DATA DISN2D,UVDISN/0.0,0.0/
    DATA HEPLUS/0.0/

    ALTCHEM = 150  !.. Initial altitude for O+ for imposing chemistry
    JITER   = 0    !.. Counts the number of Newton iterations
    N2P     = 0.0  !.. N(2P) density
    CALL RATS(0, TE, TI, TN, RTS)  !.. Get the reaction rates
    !.. PRIMPR calculates solar EUV production rates.
    !      print*,'ALT,OXN,N2N,O2N,HEN,SZAD*0.01745,TN,F107,F107A,N4S'
    !      print*,ALT,OXN,N2N,O2N,HEN,SZAD*0.01745,TN,F107,F107A,N4S
    UVDISN = OTHPR1(1)
    CALL PRIMPR(1, ALT, OXN, N2N, O2N, HEN, SZAD*0.01745, TN,&
                F107, F107A, N4S)
    UVDISN = OTHPR1(1)
    !.. Calculate secondary Production from photoelectrons
    CALL SECIPRD(ALT, SZAD, F107, F107A, TE, TN, OXN, O2N,&
                 N2N, NE, N2APRD)
    DISNP = EUVION(3, 4) + EUVION(3, 5) + EUVION(3, 6)&
            + 0.1*( PEPION(3, 1) + PEPION(3, 2) + PEPION(3, 3) )&  !.. Rydberg diss&
            + PEPION(3, 4) + PEPION(3, 5) + PEPION(3, 6)
    DISN2D = 2.0*PEPION(3, 1) + OTHPR2(3)
    DISN4S = 2.0*PEPION(3, 1) + OTHPR2(3)

    NEWT_SUCCESS = .FALSE.
    DO WHILE ( .not. NEWT_SUCCESS )
        !.. initialize variables to avoid using left over values
        HEPLUS = 0.0
        OXPLUS = 0.0
        N2PLUS = 0.0
        NOPLUS = 0.0
        O2PLUS = 0.0
        NPLUS  = 0.0
        N2P    = 0.0
        N2D    = 0.0
        OP2D   = 0.0
        OP2P   = 0.0
        N2A    = 0.0
        ! Assume that Newton is going to succeed.
        ! This will be reset if Newton's method fails.
        NEWT_SUCCESS = .TRUE.

        !.. Iterate through chemistry twice in case O+ is chemical equilibrium
        DO ITERS=1,3
            !.. K counts number of iterations for printing headers in routines
            K = K + 1

            !.. Set initial value for electron density
            !.. O+(2P) Calculate and print densities, production, loss
            PSEC   = PEPION(1, 3)           !.. Photoelectron production
            TPROD3 = EUVION(1, 3) + PSEC    !.. Add EUV and photoelectrons
            CALL COP2P(JPRINT, 7, K, ALT, RTS, OXN, O2N, N2N, NE, &
                       OP2P, TPROD3, PSEC, HEPLUS, N4S, NNO, TE)

            !.. O+(2D) Calculate and print densities, production, loss
            PSEC   = PEPION(1, 2)    !.. Photoelectron production
            TPROD2 = EUVION(1, 2)    !.. EUV
            CALL COP2D(JPRINT, 8, K, ALT, RTS, OXN, O2N, N2N, NE,&
                       OP2D, TPROD2, OP2P, HEPLUS, N4S, NNO, PSEC)

            !.. O+(4S) Calculate and print densities, production, loss.
            TPROD1 = EUVION(1, 1)
            PDISOP = EUVION(2, 4) + EUVION(2, 5) + PEPION(2, 4) + PEPION(2, 5)
            CALL COP4S(JPRINT, 4, K, ALT, RTS, OXN, O2N, N2N, NE, OXPLUS, TPROD1, &
                       OP2D, OP2P, PEPION(1, 1), PDISOP, N2PLUS, N2D, NNO, 1.0, HEPLUS)
            CALL CN2A(JPRINT, 27, K, ALT, RTS, OXN, O2N, N2N, NE, N2A, N2APRD, 0.0, &
                      0.0, 0.0)

            !..CALL CN2P(JPRINT,0,K,ALT,RTS,OXN,O2N,N2N,NE,PN2P,LN2P
            !..> ,N2P,P3X7,UVDISN,O2PLUS,NNO,N2PLUS)
            !.. N(2D) Calculate and print densities, production, loss.
            CALL CN2D(JPRINT, 16, K, ALT, RTS, OXN, O2N, N2N, NOPLUS, NE, PN2D, LN2D,&
                      N2PLUS, DISN2D, UVDISN, NPLUS, N2P, N2D, OXPLUS, NNO, N2A)
            N2D = PN2D / LN2D

            !.. N2+ Calculate and print densities, production, loss.
            CALL CN2PLS(JPRINT, 9, K, ALT, RTS, OXN, O2N, N2N, NE, N2PLUS, EUVION(3, 1),&
                        EUVION(3, 2), EUVION(3, 3), PEPION(3, 1), PEPION(3, 2), PEPION(3, 3),&
                        OP2D, OP2P, HEPLUS, NPLUS, NNO, N4S)

            !.. N+ Calculate and print densities, production, loss.
            !.. Note that N(2D) is turned off in N+ solution
            PHOTN = OTHPR2(3)  !.. N+ prod
            CALL CNPLS(JPRINT, 10, K, ALT, RTS, OXN, O2N, N2N, NE, DISNP, NPLUS, &
                       OXPLUS, N2D, OP2P, HEPLUS, PHOTN, O2PLUS, N4S, OP2D, N2PLUS, NNO)
            IF ( USER_NPLUS > 0 ) NPLUS = USER_NPLUS  !.. User specified N+
            CALL CN4S(JPRINT, 28, K, ALT, RTS, OXN, O2N, N2N, NE, PN4S, LN4S, N4S, &
                      DISN4S, N2D, N2P, OXPLUS, N2PLUS, UVDISN, NOPLUS, NPLUS, NNO, &
                      O2PLUS, PDNOSR, 1.0)

            !..N4S=PN4S/LN4S
            !.. NO Calculate and print densities, production, loss.
            CALL CNO(JPRINT, 15, K, ALT, RTS, OXN, O2N, N2N, NE, PNO, LNO,&
                     N2D, N4S, N2P, NNO, O2PLUS, OXPLUS, OTHPR2(2), OTHPR2(1), N2A, NPLUS)
            NNO = PNO / LNO     !.. NO chemical equilibrium density
            !.. Set a floor on NO density, which is needed below ~150 km at night
            IF ( NNO < 1.0E8*EXP( (100-ALT)/20 ) )&
                                 NNO = 1.0E8 * EXP( (100-ALT)/20 )
            IF ( USER_NO > 1.0 ) NNO = USER_NO  !.. substitute user specified value
            IF ( NNO > 1.5E8 )   NNO = 1.5E8      !.. Don't let NO get too big

            !.. NO+ Calculate and print densities, production, loss.
            CALL CNOP(JPRINT, 11, K, ALT, RTS, OXN, O2N, N2N, NE, TPNOP, NOPLUS, OXPLUS,&
                      N2PLUS, O2PLUS, N4S, NNO, NPLUS, N2P, PLYNOP, VCON, N2D, OP2D)
            !.. O2+ Calculate and print densities,  production,  loss.
            !.. EUV + PE production
            TPROD5  =  EUVION(2, 1)  +  EUVION(2, 2)  +  EUVION(2, 3)  +  PEPION(2, 1) &
                      +  PEPION(2, 2)  +  PEPION(2, 3)
            CALL CO2P(JPRINT, 12, K, ALT, RTS, OXN, O2N, N2N, NE, O2PPROD,&
                      O2PLUS, TPROD5, OXPLUS, OP2D, N2PLUS, NPLUS, N4S, NNO, OP2P)
        END DO

        !.. This section for chemical equilibrium densities for all species
        !.. including O+. It is used when the Newton procedure fails to get O+
        !.. Don't bother if molecular ions greater than Ne/2
        !.. It increments ALTCHEM to force this action. The ion densities are
        !.. normalized to the input NE. N+ is excluded in case it is user specified
        INEWT   = 1
        SUMIONS = OXPLUS + NOPLUS + O2PLUS + NPLUS
        IF ( ALT < ALTCHEM .OR. NOPLUS + O2PLUS > 0.85*NE ) THEN
            OXPLUS = OXPLUS * NE/SUMIONS
            NOPLUS = NOPLUS * NE/SUMIONS
            O2PLUS = O2PLUS * NE/SUMIONS
            NPLUS  = NPLUS * NE/SUMIONS
            INEWT  = 0
            RETURN
        ENDIF

        !.. Newton solution for O+ density given the electron density (NE)
        !.. Go through twice to set up the derivative (F(x+h)-F(x))/h
        !.. First O+ guess for Newton. This is important for high altitudes
        !.. because Newton may converge on first try.
        OXPLUSAVE = OXPLUS
        IF ( NE - NOPLUS - O2PLUS > 100) OXPLUS = NE - NOPLUS - O2PLUS

        SUMIONS = 0 ! Re-initialize SUMIONS
        DO WHILE ( ABS( 1. - SUMIONS/NE ) > 0.05 )
            DO ITS = 1, 2
                IF ( ITS == 1 ) H = OXPLUS * 0.0001         !.. increment for dn/dt
                IF ( ITS == 2 ) OXPLUS = OXPLUS + H       !.. increment N

                !.. N+ Calculate and print densities, production, loss.
                CALL CNPLS(JPRINT, 10, K, ALT, RTS, OXN, O2N, N2N, NE, DISNP, NPLUS, &
                           OXPLUS, N2D, OP2P, HEPLUS, PHOTN, O2PLUS, N4S, OP2D, N2PLUS, NNO)
                IF ( USER_NPLUS > 0 ) NPLUS = USER_NPLUS  !.. User specified N+

                !.. N2+ Calculate and print densities, production, loss.
                CALL CN2PLS(JPRINT, 9, K, ALT, RTS, OXN, O2N, N2N, NE, N2PLUS, EUVION(3, 1),&
                            EUVION(3, 2), EUVION(3, 3), PEPION(3, 1), PEPION(3, 2), PEPION(3, 3),&
                            OP2D, OP2P, HEPLUS, NPLUS, NNO, N4S)

                !.. NO+ Calculate and print densities,  production,  loss.
                CALL CNOP(JPRINT, 11, K, ALT, RTS, OXN, O2N, N2N, NE, TPNOP, NOPLUS, OXPLUS,&
                          N2PLUS, O2PLUS, N4S, NNO, NPLUS, N2P, PLYNOP, VCON, N2D, OP2D)

                !.. O2+ Calculate and print densities, production, loss.
                !.. EUV + PE production
                TPROD5 = EUVION(2,1) + EUVION(2,2) + EUVION(2,3) + PEPION(2,1)&
                         + PEPION(2,2) + PEPION(2,3)
                CALL CO2P(JPRINT, 12, K, ALT, RTS, OXN, O2N, N2N, NE, O2PPROD,&
                          O2PLUS, TPROD5, OXPLUS, OP2D, N2PLUS, NPLUS, N4S, NNO, OP2P)

                !.. N(2D) Calculate and print densities, production, loss.
                CALL CN2D(JPRINT, 16, K, ALT, RTS, OXN, O2N, N2N, NOPLUS, NE, PN2D, LN2D,&
                          N2PLUS, DISN2D, UVDISN, NPLUS, N2P, N2D, OXPLUS, NNO, N2A)
                N2D = PN2D/LN2D

                !.. calculation of F(x) from the total ion concentration
                FEX(ITS) = OXPLUS + NOPLUS + O2PLUS + N2PLUS + NPLUS - NE
            ENDDO

            !.. Test for convergence and add increment to O+ if not
            JITER  = JITER + 1          !.. for stopping the iterations
            DEX    = ( FEX(2) - FEX(1) )/H
            OXPLUS = OXPLUS - H - FEX(1)/DEX

            !.. If Newton fails, go back and calculate O+ chem equilibrium
            IF ( JITER > 6 .OR. OXPLUS < 0.0 .OR. OXPLUS > 1.0E7 ) THEN
                ALTCHEM = ALT + 1  !.. forces chemical equilibrium densities
                NEWT_SUCCESS = .FALSE.
                exit
            END IF

            !.. Update SUMIONS, used in convergence test
            SUMIONS = OXPLUS + NOPLUS + O2PLUS + N2PLUS + NPLUS
        END DO

        IF (NEWT_SUCCESS) THEN
            !.. Normalize ion densities to the input electron density
            OXPLUS = OXPLUS * NE/SUMIONS
            NOPLUS = NOPLUS * NE/SUMIONS
            O2PLUS = O2PLUS * NE/SUMIONS
            NPLUS  = NPLUS  * NE/SUMIONS

            !.. If O+ a very minor ion, Newton solution may not be good, force
            !.. chemical equilibrium solution instead
            IF ( OXPLUS/SUMIONS < 0.1 ) THEN
                NEWT_SUCCESS = .FALSE.
                ALTCHEM = ALT + 1  !.. forces chemical equilibrium densities
            ENDIF
        ENDIF
    END DO
    RETURN
end subroutine ChemIon


!.................... KEMPRN.FOR ......................................
!.. This file contains the chemistry routines for ions and neutrals
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine CN2D(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NOP, NE, P1, L1, &
                N2PLS, DISN2D, UVDISN, NPLUS, N2P, N2D, OPLS, NNO, N2A)
    !.......n(2d)
    IMPLICIT REAL(A-H,L,N-Z)
    DIMENSION RTS(99), LR(22), PR(22)

    PR(1) = NOP     * NE * RTS(50)
    PR(2) = N2PLS   * NE * RTS(32) * RTS(11)
    PR(3) = N2PLS   * ON * RTS(10)
    PR(4) = DISN2D
    PR(5) = RTS(63) * UVDISN
    PR(6) = RTS(65) * NPLUS * O2N
    PR(7) = N2P     * RTS(57)
    PR(8) = RTS(27) * N2A * ON
    LR(1) = ON      * RTS(15)
    LR(2) = O2N     * RTS(16)
    LR(3) = NE      * RTS(8)
    LR(4) = OPLS    * RTS(29)
    LR(5) = RTS(61)
    LR(6) = RTS(41) * NNO
    P1    = SUM( PR(1:8) )
    L1    = SUM( LR(1:6) )
    !....... EF is used to convert production rates to volume emission rates
    EF    = 1.0
    !....... This line in for volume emission rates
    !...      EF=RTS(61)*0.76/L1
    IF ( JPT == 1 .AND. JPR > 0 .AND. INT( EF + 0.1 ) /= 1) WRITE(I, 189)
    IF ( JPT == 1 .AND. JPR > 0 .AND. INT( EF + 0.1 ) == 1) WRITE(I, 191)
    IF ( JPT == 1 .AND. JPR > 0) WRITE(I,96)
    IF ( JPR > 0 ) WRITE(I,7) Z, P1/L1, (PR(K)*EF, K=1, 8), (LR(K)*N2D, K=1, 6)
    RETURN

! Format Statements
7    FORMAT(F6.1,1P,22E8.1)
96   FORMAT(2X,'ALT   [N2D]   NO++e   N2++e   N2++O    e+N2   hv+N2',&
            3X,'N++O2   N(2P)   N2A+O    +O     +O2      +e     +O+',&
            5X,'RAD     +NO')
189  FORMAT(/2X,'N(2D)',25X,'EMISSION',28X,':',20X,'Loss rate')
191  FORMAT(/2X,'N(2D)',25X,'Production',36X,':',20X,'Loss rate')
end subroutine CN2D


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine CNO(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NE, P1, L1, &
               N2D, N4S, N2P, NNO, O2P, OPLS, PDNOSR, PLYNOP, N2A, NPLUS)
    !........no
    IMPLICIT REAL(A-H,L,N-Z)
    DIMENSION RTS(99),LR(22),PR(22)

    PR(1) = RTS(16) * O2N   * N2D
    PR(2) = RTS(7)  * O2N   * N4S
    PR(3) = RTS(38) * N2P   * O2N
    PR(4) = RTS(27) * N2A   * ON
    PR(5) = RTS(22) * NPLUS * O2N          !.. Fox
    LR(1) = RTS(9)  * N4S
    LR(2) = RTS(23) * O2P
    LR(3) = RTS(24) * OPLS
    LR(4) = RTS(41) * N2D
    LR(5) = PDNOSR
    LR(6) = PLYNOP
    P1    = SUM( PR(1:5) )
    L1    = SUM( LR(1:6) )

    IF ( JPT == 1 .AND. JPR > 0 ) WRITE(I,192)
    IF ( JPR > 0 ) WRITE(I,7) Z, NNO, P1/L1, (PR(K), K=1, 5), (LR(K)*NNO, K=1, 6)
    RETURN

! Format Statements
7    FORMAT(F6.1,1P,22E9.2)
192  FORMAT(/2X,'NO',17X,'PRODUCTION',20X,':',10X,'LOSS RATES'/&
    ,4X,'ALT',3X,'[NO]',5X,'[NO]c',3X,'O2+N2D',&
    3X,'O2+N4S   N2P+O2   N2A+O    N++O2    N4S+NO   O2P+NO   O++NO'&
    ,3X,'N2D+NO   hv<1910   Lyman-a')
end subroutine CNO


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine CN4S(JPR,I,JPT,Z,RTS,ON,O2N,N2N,NE,P1,L1,N4S,DISN4S,&
                N2D,N2P,OPLS,N2PLS,UVDISN,NOP,NPLUS,NNO,O2P,PDNOSR,VCON)
    !........N(4S)
    IMPLICIT REAL(A-H,L,N-Z)
    DIMENSION RTS(99),LR(22),PR(22)

    PR(1)  = DISN4S
    PR(2)  = RTS(15) * ON      * N2D
    PR(3)  = RTS(8)  * NE      * N2D
    PR(4)  = VCON    * RTS(3)  * OPLS  * N2N
    PR(5)  = RTS(53) * RTS(11) * N2PLS * NE
    PR(6)  = RTS(62) * UVDISN
    PR(7)  = NOP     * NE      * RTS(49)
    PR(8)  = N2D     * RTS(61)
    PR(9)  = N2P     * RTS(58)
    PR(10) = RTS(25) * NPLUS   * O2N
    PR(11) = PDNOSR  * NNO
    PR(12) = NPLUS   * NNO     * RTS(81)      !..Fox
    LR(1)  = RTS(7)  * O2N
    LR(2)  = RTS(9)  * NNO
    LR(3)  = RTS(21) * O2P
    LR(4)  = RTS(79) * N2PLS     !..Fox
    P1 = SUM( PR(1:12) )
    L1 = SUM( LR(1:4) )

    IF ( JPT == 1 .AND. JPR > 0 ) WRITE(I,193)
    PR(10) = PR(10) + PR(12)  !... for printing fit
    LR(3)  = LR(3)  + LR(4)     !... for printing fit
    IF(JPR > 0) WRITE(I,7) Z, N4S, (PR(K), K=1, 11), (LR(K)*N4S, K=1, 3)
    RETURN

! Format Statements
7    FORMAT(F6.1,1P,22E8.1)
193  FORMAT(/2X,'N(4S)',38X,'PRODUCTION',46X,':',7X,'LOSS RATES'/,&
            3X,'ALT',2X,'[N4S]',2X,'hv->N+',&
            3X,'O+N2D',2X,'e+N2D',3X,'O++N2',3X,'N2++e',4X,'hv->2N',&
            2X,'NO++e',2X,'N(2D)',4X,'N(2P)   N+&X    hv+NO    +O2  ',&
            2X,' +NO  ',2X,'+O2+ , N2+')
end subroutine CN4S


!::::::::::::::::::::::::::::::: CN2PLS :::::::::::::::::::::::::::::::
!..... Simplified chemistry of N2+.  PUN2P* = production of N2+ by euv
!..... in the (X,A,B states). PEN2P* same for p.e.s (X,A,B states)
subroutine CN2PLS(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NE, N2PLS, PUN2PX, &
                  PUN2PA, PUN2PB, PEN2PX, PEN2PA, PEN2PB, OP2D, OP2P, HEPLUS, NPLUS, &
                  NNO, N4S)
    IMPLICIT REAL(A-H,L,N-Z)
    REAL PUN2PX,PUN2PA,PUN2PB,PEN2PX,PEN2PA,PEN2PB
    DIMENSION RTS(99),LR(22),PR(22)

    PR(1)  = PUN2PX
    PR(2)  = PUN2PA
    PR(3)  = PUN2PB
    PR(4)  = PEN2PX
    PR(5)  = PEN2PA
    PR(6)  = PEN2PB
    PR(7)  = RTS(19) * OP2D   * N2N
    PR(8)  = RTS(20) * OP2P   * N2N
    PR(9)  = RTS(44) * HEPLUS * N2N
    PR(10) = RTS(82) * NPLUS  * NNO   !..Fox
    LR(1)  = RTS(10) * ON
    LR(2)  = RTS(11) * NE
    LR(3)  = RTS(17) * O2N
    LR(4)  = RTS(99) * ON
    LR(5)  = RTS(79) * N4S     !..Fox
    LR(6)  = RTS(80) * NNO     !..Fox
    N2PLS  = SUM( PR(1:10) ) / SUM( LR(1:6) )

    IF( JPT == 1 .AND. JPR > 0 ) WRITE(I,95)
    PR(9) = PR(9) + PR(10)     !.. for printing fit
    LR(5) = LR(5) + LR(6)      !.. for printing fit
    IF(JPR > 0) WRITE(I,7) Z, N2PLS, (PR(K), K=1, 9), (LR(K)*N2PLS, K=1, 5)
    RETURN

! Format Statements
7    FORMAT(F6.1,1P,22E8.1)
95   FORMAT(/2X,'N2+',29X,'PRODUCTION',45X,':',12X,'LOSS RATES'/,&
            3X,'ALT  [N2+]  EUV-X   EUV-A    EUV-B   PE-X',&
            5X,'PE-A    PE-B  O+2D+N2  O+2P+N2  He++N2  O+N2+',&
            2X,'e+N2+  O2+N2+  N2++O  Other')
end subroutine CN2PLS


!:::::::::::::::::::::::::::::: CNOP ::::::::::::::::::::::::::::::::::
subroutine CNOP(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NE, P1, NOP, OPLS, &
                N2PLS, O2P, N4S, NNO, NPLUS, N2P, PLYNOP, VCON, N2D, OP2D)
    !........no+
    IMPLICIT REAL(A-H,L,N-Z)
    DIMENSION RTS(99),LR(22),PR(22)

    PR(1)  = VCON  * RTS(3) * N2N * OPLS
    PR(2)  = N2PLS * ON     * RTS(10)
    PR(3)  = O2P   * N4S    * RTS(21)
    PR(4)  = O2P   * NNO    * RTS(23)
    !.. N +   +  O2 -> O2 +   +  N(2D,4S) or NO +   +  O(1S)
    PR(5)  = (RTS(30) + RTS(66) + RTS(59)) * NPLUS * O2N
    PR(6)  = RTS(37) * N2P  * ON
    PR(7)  = RTS(24) * OPLS * NNO
    PR(8)  = PLYNOP  * NNO
    PR(9)  = O2P     * N2D  * RTS(77)    !..Fox
    PR(10) = N2PLS   * NNO  * RTS(80)    !..Fox
    PR(11) = NPLUS   * NNO  * RTS(81)    !..Fox
    PR(12) = RTS(83) * NNO  * OP2D       !..Fox
    PR(13) = OP2D    * RTS(90) * N2N     !.. -> NO +   +  N, Li et al. [1997]
    LR(1)  = NE      * RTS(5)
    P1     = SUM( PR(1:13) )
    NOP    = P1/LR(1)

    IF ( JPT == 1 .AND. JPR > 0 ) WRITE(I,96)
    !PR(9)=PR(9)+PR(10)+PR(11)+PR(12)+PR(13)
    IF ( JPR > 0 ) WRITE(I,7) Z,NOP,(PR(K),K=1,13),LR(1)*NOP
    RETURN

! Format Statements
7    FORMAT(F6.1,1P,22E9.2)
96   FORMAT(/2X,'NO+',31X,'PRODUCTION',48X,':',2X,'LOSS RATES'/,&
            3X,'ALT',3X,'[NO+]',4X,'O++N2',3X,'N2++O',3X,'O2++N4S',&
            3X,'O2++NO',3X,'N++O2',4X,'N2P+O',3X,'O++NO   hv+NO',&
            5X,'O2++N2D   N2++NO   N++NO   OP2D+NO   OP2D+N2  NO++e')
end subroutine CNOP


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine CO2P(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NE, P1, &
                O2P, TPROD5, OPLS, OP2D, N2PLS, NPLUS, N4S, NNO, OP2P)
    !........o2+
    IMPLICIT REAL(A-H,L,N-Z)
    DIMENSION RTS(99),LR(22),PR(22)

    !......... TPROD5=euv @ p.e. production
    PR(1) = TPROD5
    PR(2) = RTS(4)  * O2N   * OPLS
    PR(3) = RTS(43) * OP2D  * O2N
    PR(4) = RTS(17) * O2N   * N2PLS
    PR(5) = RTS(25) * NPLUS * O2N
    PR(6) = RTS(86) * OP2P  * O2N           !.. Fox
    PR(7) = RTS(65) * NPLUS * O2N
    LR(1) = RTS(6)  * NE
    LR(2) = RTS(21) * N4S
    LR(3) = RTS(23) * NNO
    P1    = SUM( PR(1:7) )
    O2P   = P1/SUM( LR(1:3) )

    IF ( JPT == 1 .AND. JPR > 0 ) WRITE(I,97)
    IF(JPR > 0) WRITE(I,7) Z,O2P,(PR(K),K=1,7),(LR(K)*O2P,K=1,3)
    RETURN

! Format Statements
7    FORMAT(F6.1,1P,22E9.2)
97   FORMAT(/2X,'O2+',22X,'PRODUCTION',24X,':',12X,'LOSS RATES'/,&
            3X,'ALT',3X,'[O2+]',3X,'hv+O2',3X,'O++O2',3X,'O+(2D)+O2',&
            4X,'N2++O2   N++O2   O+(2P)+O2  O2++e   O2++N   O2++NO')
end subroutine CO2P


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine COP4S(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NE, OPLS, TPROD1, OP2D, &
                 OP2P, PEPION, PDISOP, N2PLS, N2D, NNO, VCON, HEPLUS)
    !...........o+(4s)
    IMPLICIT REAL(A-H,L,N-Z)
    DIMENSION RTS(99),LR(22),PR(22)

    !.........pr(1)= euv production of o+(4s)
    PR(1)  = TPROD1
    PR(2)  = OP2D    * NE * RTS(12)
    PR(3)  = OP2P    * ON * RTS(26)
    PR(4)  = PEPION
    PR(5)  = PDISOP
    PR(6)  = RTS(99) * N2PLS * ON
    PR(7)  = OP2P    * NE * RTS(14)
    PR(8)  = OP2P    * 0.047
    PR(9)  = RTS(28) * ON * OP2D
    PR(10) = RTS(85) * OP2P * O2N              !.. Fox
    PR(11) = HEPLUS  * O2N * (RTS(91) + RTS(93))  !..Fox
    PR(12) = RTS(95) * NNO * HEPLUS            !..Fox
    !      PR(13) = RTS(22) * NPLUS * O2N             !..Fox
    PR(13) = 0.0
    LR(1)  = N2N * VCON * RTS(3)
    LR(2)  = O2N * RTS(4)
    LR(3)  = NNO * RTS(24)
    LR(4)  = N2D * RTS(29)    !.... small loss?? ..Fox
    !..LR(4) = (LR(1) + LR(2) + LR(3))  !.. total loss for printing
    PR(10) = SUM(PR(10:13))
    PRTOT  = SUM(PR(1:10))
    LRTOT  = SUM(LR(1:4))
    OPLS   = PRTOT/LRTOT

    IF ( JPT == 1 .AND. JPR > 0 ) WRITE(I,98)
    IF ( JPR > 0 ) WRITE(I,7) Z,OPLS,(PR(K),K=1,10),(LR(K)*OPLS,K=1,4)
    RETURN

! Format Statements
7    FORMAT(F6.1,1P,22E8.1)
98   FORMAT(/2X,'O+',41X,'PRODUCTION',39X,':',10X,'LOSS RATES'/,&
            ' ALT    [O+]   hv+O  O+(2D)+e O+(2P)+O   e*+O  O2-diss  ',&
            'N2++O  O+(2P)+e O+(2P) O+O+(2D)   Other  +N2     +O2    ',&
            '+NO   +N2D')
end subroutine COP4S


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine COP2D(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NE, OP2D, TPROD2, OP2P, &
                 HEPLUS, N4S, NNO, PSEC)
    !.......o+(2d)
    IMPLICIT REAL (A-H,L,N-Z)
    DIMENSION RTS(99),LR(22),PR(22)

    PR(1) = TPROD2                ! EUV  prod
    PR(2) = OP2P    * NE * RTS(13)
    PR(3) = OP2P    * 0.171
    PR(4) = HEPLUS  * O2N * RTS(76)     !..Fox
    PR(5) = PSEC
    LR(1) = RTS(19) * N2N
    LR(2) = 7.7E-5                 !.. radiation at 3726 and 3729 A
    LR(3) = NE      * RTS(12)
    LR(4) = ON      * RTS(28)
    LR(5) = RTS(43) * O2N
    LR(6) = RTS(83) * NNO      !..Fox
    LR(7) = RTS(84) * N4S      !..Fox
    LR(8) = RTS(90) * N2N      !.. -> NO +   +  N, Li et al. [1997]
    OP2D  = SUM( PR(1:5) )/SUM( LR(1:8) )

    IF ( JPT == 1 .AND. JPR > 0 ) WRITE(I,99)
    IF ( JPR > 0 ) WRITE(I,7) Z,OP2D,PR(1),PR(5),PR(2),PR(3),PR(4),(LR(K)*OP2D,K=1,8)
    RETURN

! Format Statements
7    FORMAT(F6.1,1P,22E9.2)
99   FORMAT(/2X,'O+(2D)',13X,'PRODUCTION',27X,':',18X,'LOSS RATES'/,&
            3X,'ALT',3X,'[O+2D]',3X,'hv+O',4X,'e*+O',4X,'O+2P+e',3X,&
            'O+2P>hv',2X,'He++O2     +N2    E3726_29    +e       +O',&
            '      +O2      +NO     +N  +N2>NO+')
end subroutine COP2D

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine COP2P(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NE, OP2P, TPROD3, PSEC, &
                 HEPLUS, N4S, NNO, TE)
    !.......o+(2p)
    IMPLICIT REAL(A-H,L,N-Z)
    DIMENSION RTS(99),LR(22),PR(22)

    PR(1) =  0.0
    IF ( TPROD3 >= PSEC ) PR(1) = TPROD3 - PSEC
    PR(2) =  PSEC
    PR(3) =  HEPLUS  * O2N * RTS(92)        !..Fox
    LR(1) =  RTS(26) * ON
    LR(2) =  RTS(20) * N2N
    LR(3) =  RTS(13) * NE
    LR(4) =  0.218
    LR(5) =  RTS(14) * NE
    LR(6) = (RTS(85) + RTS(86)) * O2N    !..Fox
    LR(7) =  RTS(87) * N4S              !..Fox
    LR(8) =  RTS(88) * NNO              !..Fox
    OP2P  = (TPROD3  + PR(3))/SUM( LR(1:8) )

    IF(JPT == 1.AND.JPR > 0) WRITE(I,100)
    IF(JPR > 0) WRITE(I,7) Z,OP2P,PR(1),PR(2),PR(3),(LR(K)*OP2P,K=1,8),ON,N2N,NE,TE,OP2P*0.218*0.781
    RETURN

7    FORMAT(F6.1,1P,22E9.2)
100  FORMAT(/2X,' O+(2P)',6X,'PRODUCTION',10X,':',12X,'LOSS RATES'/,&
            3X,'ALT   [O+2P]    hv+O     e*+O  He++O2      +O',7X,'+N2',&
            6x,'+e       RAD      +e      +O2      +N4S     +NO',&
            6x,'OX       N2        e      Te       E7320')
end subroutine COP2P


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine CNPLS(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NE, DISNP, &
                 NPLUS, OPLS, N2D, OP2P, HEPLUS, PHOTN, O2P, N4S, OP2D, N2PLS, NNO)
    !........n+
    IMPLICIT REAL(A-H,L,N-Z)
    DIMENSION RTS(99),LR(22),PR(22)

    PR(1) = DISNP
    PR(2) = RTS(29) * OPLS * N2D
    PR(3) = 0
    PR(4) = RTS(45) * HEPLUS * N2N
    PR(5) = PHOTN
    PR(6) = O2P     * N2D * RTS(78)    !..Fox
    PR(7) = N2PLS   * N4S * RTS(79)    !..Fox
    PR(8) = OP2D    * N4S * RTS(84)    !..Fox
    PR(9) = RTS(94) * NNO * HEPLUS     !..Fox
    LR(1) = RTS(30) * O2N              !..Fox
    LR(2) = RTS(25) * O2N              !..Fox
    LR(3) = RTS(22) * O2N              !..Fox
    LR(4) = RTS(65) * O2N              !..Fox
    LR(5) = RTS(66) * O2N              !..Fox
    LR(6) = RTS(31) * ON               !..Fox
    CNPLUS = 0.0
    IF ( SUM( LR(1:3) ) > 0.0) CNPLUS = SUM( PR(1:9) ) / SUM( LR(1:6) )
    NPLUS = CNPLUS

    IF(JPT == 1.AND.JPR > 0) WRITE(I,101)
    IF(JPR > 0) WRITE(I,7) Z,NPLUS,CNPLUS,(PR(K),K=1,9),(LR(K)*NPLUS,K=1,6)
    RETURN

7    FORMAT(F6.1,1P,22E9.2)
101  FORMAT(/2X,'N+',20X,'PRODUCTION',71X,':',8X,'LOSS RATES'/,&
            4X,'ALT   [N+]   [N+]c     hv+N2   O++N2D  O+2P+N2',3X,&
            'He++N2',3X,' hv+N   O2++N2D  N2++N4S O+(2D)+N4S  He++NO',&
            3X,'N++O2    N++O2    N++O2    N++O2    N++O2    N++O')
end subroutine CNPLS


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine CN2A(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NE, &
                N2A, P3X1, P3X2, P3X3, P3X4)
    !........n2(a3sigma) and total LBH
    IMPLICIT REAL(A-H,L,N-Z)
    REAL P3X1,P3X2,P3X3,P3X4
    DIMENSION RTS(99),LR(22),PR(22)

    !....... pr(1,2,3)= electron impact excitation of n2(a,b,c) states
    PR(1) = P3X1
    PR(2) = P3X2
    PR(3) = P3X3
    LR(1) = RTS(36) * ON
    LR(2) = RTS(27) * ON
    LR(3) = 0.57
    N2A = SUM( PR(1:3) ) / SUM( LR(1:3) )

    IF(JPT == 1.AND.JPR > 0) WRITE(I,102)
    IF(JPR > 0) WRITE(I,7) Z,N2A,(PR(K),K=1,3),(LR(K)*N2A,K=1,3),P3X4
    RETURN

7    FORMAT(F6.1,1P,22E9.2)
102  FORMAT(/2X,'N2(A)',12X,'PRODUCTION',13X,':',5X,'LOSS RATES',&
            3X,':  Total LBH'/,&
            3X,'ALT',3X,'N2(A)',3X,'e*->N2A',3X,'e*->N2B',3X,'e*->N2C',2X,&
            'N2A>O1S',2X,'N2A>NO',2X,'RAD',5X,'LBH')
end subroutine CN2A


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine CN2P(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NE, P1, L1, &
                N2P, P3X7, UVDISN, O2P, NNO, N2PLUS)
    !....... N(2P). the rates are from Zipf et al JGR, 1980 p687
    !.... 21-AUG-1992. Added N2+ recombination source
    IMPLICIT REAL(A-H,L,N-Z)
    DIMENSION RTS(99),LR(22),PR(22)

    PR(1) = P3X7
    PR(2) = RTS(64) * UVDISN
    PR(3) = RTS(73) * RTS(11) * N2PLUS * NE
    LR(1) = RTS(37) * ON
    LR(2) = RTS(38) * O2N
    LR(3) = RTS(39) * O2P
    LR(4) = RTS(40) * NNO
    LR(5) = RTS(57)
    LR(6) = RTS(58)
    LR(7) = (RTS(96) + RTS(97)) * NE    !..Fox
    P1 = SUM( PR(1:3) )
    L1 = SUM( LR(1:7) )

    IF(JPT == 1.AND.JPR > 0) WRITE(I,103)
    IF(JPR > 0) WRITE(I,7) Z,N2P,(PR(K),K=1,3),(LR(K)*N2P,K=1,7)
    RETURN

7    FORMAT(F6.1,1P,22E9.2)
103  FORMAT(/2X,'N(2P)',9X,'PRODUCTION',17X,':',20X,'LOSS RATES'/,&
            3X,'ALT',3X,'[N2P]',3X,'e+N2',5X,'hv+N2',3X,'e+N2+',6X,'+O  ',&
            3X,'+O2      +O2+      +NO       +2D     +4S      +e')
END


!:::::::::::::::::::::::::::::: CNOPV ::::::::::::::::::::::::::::::::::
!...... This routine calculates the vibrational distribution of no+
!...... Uses AFRL report by Winick et al. AFGL-TR-87-0334, Environmental
!...... Research papers, NO. 991, "An infrared spectral radiance code
!...... for the auroral thermosphere (AARC)
!...... Written by P. Richards in February 2004
subroutine CNOPV(JPR, I, JPT, Z, RTS, ON, O2N, N2N, NE, P1, NOP, OPLS, &
                 N2PLS, O2P, N4S, NNO, NPLUS, N2P, PLYNOP, VCON, N2D, OP2D)

    IMPLICIT NONE

    INTEGER JPR, I, JPT, INV, IJ, IV, K
    PARAMETER (INV=20)              !.. NO+(v) array dimensions
    REAL Z, RTS(99)                 !.. Altitude,  rate coefficients&
    REAL ON, O2N, N2N, NE           !.. O,  N2,  O2,  electron densities&
    REAL P1                         !.. total source output for finding [e]&
    REAL NOP, OPLS, N2PLS, O2P      !.. NO+,  )+, N2+, O2+ densities&
    REAL N4S, NNO, NPLUS, N2P       !.. N(4S),  NO,  N+,  N(2P) densities&
    REAL PLYNOP, VCON               !.. Lyman-a source,  N2(v) rate factor&
    REAL N2D, OP2D                  !.. N(2D),  O+(2D) densities&
    REAL NOPV(INV), NOPTOT          !.. NO+(v) densities and total NO+ density&
    REAL LR(22), PR(22)             !.. storage for NO+ sources and sinks&
    REAL EINSCO1(INV), EINSCO2(INV) !.. Einstein coeffs for delv=1, 2&
    REAL LRV(INV)                   !.. NO+(v) + e rate factors&
    REAL PNOPV, LNOPV               !.. Sources and sinks of NO+(v)&
    REAL PCASC, LRAD                !.. Temp total cascade source,  sink&
    REAL K_N2_Q, P_N2_Q, L_N2_Q     !.. N2 queching rate :- coeff,  source,  sink
    !.. Fractions of each source going to each vib. level. Assume
    !.. N+ + O2 fractions for each source. NEED UPDATE
    REAL PRV1(INV), PRV2(INV), PRV3(INV), PRV4(INV), PRV5(INV), PRV6(INV), &
         PRV7(INV), PRV8(INV), PRV9(INV), PRV10(INV), PRV11(INV), PRV12(INV)
    DATA PRV1/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
    DATA PRV2/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
    DATA PRV3/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
    DATA PRV4/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
    DATA PRV5/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
    DATA PRV6/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
    DATA PRV7/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
    DATA PRV8/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
    DATA PRV9/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
    DATA PRV10/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
    DATA PRV12/0,.05,.07,.09,.11,.13,.14,.17,.07,.01,.02,.06,.08,7*0/
    DATA K_N2_Q/7.0E-12/  !.. Quenching rate coeff. by N2
    DATA EINSCO1/0.0,10.9,20.2,28.4,35.5,41.5,46.6,50.9,54.3,57.0,&
        58.9,60.2,60.8,60.7,59.9,5*60.0/  !.. Einstein coeff delv=1
    DATA EINSCO2/0.0,0.0,.697,1.93,3.61,5.74,8.24,11.1,14.2,17.7,&
        21.3,25.1,29.0,33.2,37.4,5*40.0/ !.. Einstein coeff delv=1
    !.. rate factors for NO+(v)+e -> N + O. Sheehan and St-Maurice 2004
    DATA LRV/1.0,19*0.3333/

    !... Evaluate total production and loss rates
    PR(1)  = VCON  * RTS(3) * N2N * OPLS    !.. N2  +  O + 
    PR(2)  = N2PLS * ON  * RTS(10)          !.. N2 +   +  O
    PR(3)  = O2P   * N4S * RTS(21)          !.. O2 +   +  N(4S)
    PR(4)  = O2P   * NNO * RTS(23)          !.. O2 +   +  NO
    !.. N+   +  O2 -> O2 +   +  N(2D,4S) or NO+   +  O(1S)
    PR(5)  = (RTS(30) + RTS(66) + RTS(59)) * NPLUS * O2N
    PR(6)  = RTS(37) * N2P * ON             !.. N2 +   +  O
    PR(7)  = RTS(24) * OPLS * NNO           !.. O +   +  NO
    PR(8)  = PLYNOP  * NNO                  !.. Lyman-a  +  NO
    PR(9)  = O2P     * N2D * RTS(77)        !.. Fox: O2 +   +  N(2D)
    PR(10) = N2PLS   * NNO * RTS(80)        !.. Fox: N2 +   +  NO
    PR(11) = NPLUS   * NNO * RTS(81)        !.. Fox: N +   +  NO
    PR(12) = RTS(83) * NNO * OP2D           !.. Fox: O + (2D)  +  NO
    PR(13) = OP2D    * RTS(90) * N2N        !.. -> NO +   +  N, Li et al. [1997]
    LR(1)  = NE      * RTS(5)               !.. NO +   +  e

    !..Total source term used in main program to calculate [e]
    P1          = SUM( PR(1:13) )
    NOP         = P1/LR(1)         !.. NO+ density
    NOPV(1:INV) = 0.0
    NOPTOT      = 0.0

    !.. loop down evaluating the vibrational populations. Must start
    !.. less than INV-1 because need cascade from v+2
    DO IV=INV-4,1,-1

        !... chemical production for v=IV = total source * fraction
        PNOPV = PR(1)*PRV1(IV)   + PR(2)*PRV2(IV)   &
               +PR(3)*PRV3(IV)   + PR(4)*PRV4(IV)   &
               +PR(5)*PRV5(IV)   + PR(6)*PRV6(IV)   &
               +PR(7)*PRV7(IV)   + PR(8)*PRV8(IV)   &
               +PR(9)*PRV9(IV)   + PR(10)*PRV10(IV) &
               +PR(11)*PRV11(IV) + PR(12)*PRV12(IV)

        !.. cascade production from v+1, v+2
        PCASC  = NOPV(IV+1)*EINSCO1(IV+1) + NOPV(IV+2)*EINSCO2(IV+2)
        LRAD   = EINSCO1(IV) + EINSCO2(IV)  !.. total radiative loss
        L_N2_Q = K_N2_Q * N2N               !.. sink of quanta by N2 quenching
        IF (IV == 1) L_N2_Q = 0.0
        P_N2_Q = K_N2_Q * N2N * NOPV(IV+1)  !.. source of quanta by N2 quenching
        LNOPV  = LR(1) * LRV(IV)            !.. recombination rate for level IV

        !.. evaluate NO+(iV) population
        NOPV(IV) = ( PNOPV + PCASC + P_N2_Q )/( LNOPV + LRAD + L_N2_Q)
        NOPTOT   = NOPTOT + NOPV(IV)  !.. total NO+ concentration

        !... diagnostic print. Set alt range to invoke
        IF ( JPR > 0 .AND. Z >= 0 .AND. Z < 10 )&
            WRITE(6,'(F10.1,I7,1P,22E10.2)') &
                Z,IV,PNOPV,P1,PCASC,P_N2_Q,LNOPV,LRAD,L_N2_Q,NOPV(IV),NOP,NOPTOT
    ENDDO

    PR(9) = SUM( PR(9:12) )                 !.. for printing only
    IF ( JPT == 1 .AND. JPR > 0 ) WRITE(I,96)
    IF ( JPR > 0 ) WRITE(I,'(F6.1,1P,22E9.2)')  Z,(NOPV(K),K=1,4),(PR(K),K=1,9),LR(1)*NOP
    RETURN

96   FORMAT(/5X,'NO+',34X,'PRODUCTION',69X,':LOSS RATE'/,&
            3X,'ALT NO+(v=0) NO+(v=1) NO+(v=2)  NO+(v=3) O++N2    N2++O',&
            4X,'O2++N4S   O2++NO   N++O2    N2P+O   O++NO   hv+NO',&
            4X,'Other_P   NO++e')
end subroutine CNOPV


!.................................... RATES.FOR .................
!.... This is the reaction rate subroutine for the FLIP model. It takes
!.... temperatures as input and puts the rates in array RTS. It includes
!.... reaction rates, efficiencies for various products, and Einstein
!.... coefficients. For a complete set of references see Fox and Sung JGR
!.... 2001, page 21,305. Rates different from Fox and Sung indicated by PGR
subroutine RATS(J,TE,TI,TN,RTS)

    IMPLICIT REAL(A-H,L,N-Z)
    REAL TE,TI,TN,RTS,T13,TOT_NP_O2_RATE
    DIMENSION RTS(99)

    !.. zero out array
    RTS(1:99) = 0.0

    !.. O + H+ -> O+ + H      Fox and Sung [2001]
    ZED    = 1 + 0.6*exp(-228.0/TN) + 0.2*exp(-326.0/TN)

    RTS(1) = 8 * 6.4E-10/(9 * ZED) * &
             ( exp(-232.1/TI) + 0.6*exp(-228.0/TI) + 0.2*exp(-326.0/TI) )
    !.. O+ + H -> O + H+    Anicich et al. [1993]
    RTS(2) = 6.4E-10
    !.. O+ + N2 --> NO+ + N,   Hierl et al.[1997]
    !.. The Hierl et al. [1997] lab rate is contaminated by N2(v)
    !.. for T > 1300K. Therefore, the Hierl et al. rate is not really
    !.. appropriate  in the ionosphere. The IDC model uses the Hierl et
    !.. al. rate because it does not solve for N2(v). The FLIP model
    !.. solves for N2(v) and uses the St. Maurice and Torr rate (JGR,1978,p969)
    IF(TI <= 1000) RTS(3)=1.2E-12*(300/TI)**0.45     !.. Hierl et al.[1997]
    IF(TI > 1000) RTS(3)=7.0E-13*(TI/1000)**2.12    !.. Hierl et al.[1997]
    !.. O+ + O2 -> O2+ + O,   Lindinger et al. [1974]
    !.. Hierl et al. lists different rates. Hierl et al. [1997] not
    !.. used above 1600 because rates are contaminated by O2(v) for
    !.. T > 1000K. We don't know the vibrational state in the
    !.. thermosphere. This fit was done by PGR May 2009. It is similar
    !.. to Fox and Sung but does not increase sharply above 1000K.
    IF(TI <= 1600) RTS(4) = 1.6E-11*(300/TI)**0.52
    IF(TI > 1600)  RTS(4) = 6.7E-12*(TI/1600)**0.6
    !.. NO+ + e -> N + O    Walls and Dunn [1974)
    !.. Vejby-Christensen et al [1998] gives 4.0E-7*(300/TE)**0.5
    !.. Torr and Torr [1979] gives 4.3E-7*(300/TE)**0.83(-0.16,+.08)
    !.. Sheehan and St. Maurice gives 3.5E-7*(300/TE)**0.65
    RTS(5) = 4.0E-7*(300/TE)**0.85

    !.. O2+ + e -> O + O   Mehr and Biondi (1969)
    IF(TE <= 1200) RTS(6) = 1.953E-7*(300/TE)**0.70
    IF(TE > 1200)  RTS(6) = 7.389E-8*(1200/TE)**0.56
    !..  O2 + N(4S)-> NO + O           Baulch et al.[1994]
    RTS(7) = 1.5E-14*TN*EXP(-3270.0/TN)
    !..   N(2D) + e -> N(4S) + e     Berrington and Burke [1981]
    RTS(8) = 3.86E-10*(TE/300.)**0.81
    !..   NO + N(4S) -> N2 + O      Lee et al. [1978]
    RTS(9) = 3.4E-11
    !.. N2+ + O -> NO+ + N   Scott et al.[1999]
    IF(TI <= 1500) RTS(10)= 1.33E-10*(300/TI)**0.44
    IF(TI > 1500)  RTS(10)= 6.55E-11*(1500/TI)**(-0.2)

    !.. N2+ + e -> N + N  Mehr and Biondi (1969)
    RTS(11) = 2.2E-7*(300/TE)**0.39    !.. Zipf (1980)
    !.. O+(2D) + e -> O+(4S) + e   McLaughlin and Bell (1998)
    !.. Henry [1969] gives 7.8E-8*(300/TE)**0.5
    RTS(12) = 6.03E-8*SQRT(300/TE)
    !.. O+(2P) + e ->  O+(2D) + e   McLaughlin and Bell (1998)
    !.. RTS(13)+RTS(14) agrees with Walker et al (1975) and
    !.. Chang et al (1993)
    RTS(13) = 1.84E-7*SQRT(300/TE)
    !.. O+(2P) + e -> O+(4S) + e  McLaughlin and Bell (1998)
    RTS(14) = 3.03E-8*SQRT(300/TE)
    !.. N(2D) + O ->  N(4S) + O  Fell et al.[1990]. Lin and Kaufman[1971]
    RTS(15) = 6.9E-13

    !.. N(2D) + O2 -> NO + O  Herron[1999]. Shihira et al.[1994]
    RTS(16) = 9.7E-12*EXP(-185.0/TN)
    !.. N2+ + O2 -> O2+ + N2   Scott et al.[1999]
    IF(TI < 1000) THEN
        RTS(17) = 5.1E-11*(300/TI)**1.16
    ELSE IF(TI <= 2000) THEN
        RTS(17) = 1.26E-11*(TI/1000)**0.67
    ELSE
        RTS(17) = 2.39E-11
    ENDIF
    !.. thermal electron excitation of O(1D); Rees et al 1967 pss, p1097 .....
    RTS(18) = 1.1E-10*SQRT(TE)*EXP(-2.27E+4/TE)*&
              ( 0.406 + 0.357E-4*TE - (0.333 + 0.183E-4*TE) * EXP(-1.37E4/TE)&
                - (0.456 + 0.174E-4*TE)*EXP(-2.97E4/TE) )
    !.. N2 + O+(2D) -> N2+ + O
    !..RTS(19) = 8.0E-10                   !.. Johnson and Biondi
    RTS(19) = 1.50E-10*(300/Ti)**(-0.55)   !.. Li et al by PGR
    !.. N2 + O+(2P) -> N2+ + 0    Fox
    !.. RTS(20) = 6.2E-10*EXP(-340/TI)   !.. Li et al from Fox wrong
    RTS(20) = 2.0E-10*(300/Ti)**(-0.55)    !.. Li et al by PGR

    !.. O2+ + N(4S) -> NO+ + 0   Scott et al.[1999]
    RTS(21) = 1.0E-10
    !.. N+ + O2 -> O+ + NO
    !.. Torr and Torr gives 6.0E-10 for total N+ + O2 reaction rate
    !.. Dotan et al [1997] from Fox and Sung gives
    !IF(TI < =  1000) TOT_NP_O2_RATE=2.02E-10*(300/TI)**(-0.45)
    !IF(TI > 1000) TOT_NP_O2_RATE = 3.49E-10
    !.. does not seem to be correct. Probably vibrationally excited O2
    !.. Branching ratios for N+ + O2 from O'Keefe et al J. Chem. Phys. 1986
    !.. NO+ +O(3P)  =  .09, NO+ + O(1D) = .36, O2+ + N(4S) = 0.35,
    !.. O2+ + N(2D)  =  0.15, O+(4S) + NO = .05
    TOT_NP_O2_RATE = 6.0E-10                !.. Total N+ + O2 rate
    RTS(22) = 0.05*TOT_NP_O2_RATE
    !.. O2+ + NO -> NO+ + O2 Midey and Viggiano [1999]
    RTS(23) = 4.5E-10 * 1.0000
    !.. O+ + NO -> O + NO+   Dotan and Viggiano [1999]
    IF(TI <= 300) RTS(24) = 7.0E-13*(300/TI)**0.66
    IF(TI > 300)  RTS(24) = 7.0E-13*(TI/300)**0.87
    !.. N+ + O2 -> O2+ + N(4S)
    RTS(25) = 0.35*TOT_NP_O2_RATE

    !.. O+(2P) + O -> O+(4S) + O
    !..RTS(26)=5.2E-10  !.. Fox appears to be wrong
    !.. (Chang et al., JGR 1993) c.f. 5.2E-11  (Rusch)
    RTS(26) = 4.0E-10
    !.. N2(A3sig) + O -> NO + N(2D)
    RTS(27) = 2.0E-11       !..see Campbell et al. 2006
    RTS(27) = 0.000000      !.. Torr and Torr value
    !.. O+(2D) + O ->  O+(4S) + O  Torr and Torr [1980]
    RTS(28) = 1.0E-11
    !.. O+ + N(2D) -> O + N+  Constantinides et al.[1979].Bates[1989]
    RTS(29) = 1.3E-10
    !.. O2 + N+ -> O(3P) + NO+
    !.. Branching ratio from O'Keefe et al J. Chem. Phys. 1968
    RTS(30) = 0.09*TOT_NP_O2_RATE

    !.. O + N+ -> O+ + N   Constantinides et al.[1979].Bates[1989]
    RTS(31) = 2.2E-12
    !.. Efficiency for   N2+ + e -> N(2D) + N(2D)
    RTS(32) = 1.46
    !.. N2 + O(1D) -> O + NO
    RTS(33) = 1.8E-11*EXP(107.0/TN)
    !.. O2 + O(1D) -> O + O2
    RTS(34) = 3.2E-11*EXP(67/TN)
    !.. O2 + N(4S) -> O(1S) + NO. Kopp et al. 1977, JGR, p4715
    RTS(35) = 2.5E-11

    !.. N2(A3sig) + O -> O(1S) + N2
    RTS(36) = 2.5E-11*EXP(TN/298)**0.55    !..  see Campbell et al. 2006
    RTS(36) = 2.0E-11                      ! .. Torr et al.
    !.. N(2P) + O -> products (N(2D,4S) and NO+) and O(3P,1D)
    !.. from Piper et al 1993, J. Chem. Phys. vol 98 page 8560.
    RTS(37) = 1.7E-11
    !.. N(2P) + O2 -> NO + O
    RTS(38) = 3.9E-12*EXP(-60/TN)
    !.. N(2P) quenching rates(O2+,NO) from Zipf et al jgr 1980 p687
    RTS(39) = 2.2E-11
    RTS(40) = 1.8E-10
    
    !.. N(2D) + NO -> N2 + O
    RTS(41) = 6.7E-11
    !.. efficiency N2+ + O -> N2 + O+(4S)
    IF(TI <= 1500) RTS(42) = 7.0E-12*(300/TI)**0.21
    IF(TI > 1500)  RTS(42) = 4.83E-12*(1500/TI)**(-0.41)
    RTS(42) = RTS(42)/RTS(10)    !.. converts to efficiency
    !.. O+(2D) + O2 -> O2+ + O   Fox
    RTS(43) = 7.0E-10
    !.. He+ + N2 -> He + N2+
    RTS(44) = 5.2E-10
    !.. He+ + N2 -> He + N+
    RTS(45) = 7.8E-10

    !.. O(1S)+ e -> O(1D) + e
    RTS(46) = 8.5E-9
    !.. O(1S)+ e -> O(3P) + e
    RTS(47) = 1.56E-10*(TE/300)**0.94
    !.. O(1S) + O2 -> O2 + O
    RTS(48) = 4.4E-12*EXP(-815.0/TN)
    !.. NO+ + e -> N(4S) + O
    RTS(49) = 0.15 * RTS(5)
    !.. NO+ + e -> N(2D) + O
    RTS(50) = 0.85 * RTS(5)

    !.. O2+ + e -> O(1D) + O
    RTS(51) = 1.11 * RTS(6)
    !.. O2+ + e -> O(1S) + O
    RTS(52) = 0.05 * RTS(6)
    !.. Efficiency for   N2+ + e -> N(4S) + N(2D)
    RTS(53) = 0.46
    !.. O(1D) -> O + 6300 + 6364
    RTS(54) = 0.00934
    !.. O(1S) -> O(1D) + 5577
    RTS(55) =  1.06

    !.. O(1S) -> O(3P) + hv (2972) RTS(56)= 4.5E-2 !.. old value
    RTS(56) = 0.10 * RTS(55)  !.. From Slanger, Spring AGU 2005
    !.. N(2P) -> N(2D) + hv
    RTS(57) = 7.9E-2
    !.. N(2P) -> N(4S) + hv
    RTS(58) = 5.0E-3
    !.. N+ + O2 -> NO+ + O(1S) Langford et al., PSS, 33,1225,1985
    RTS(59) = 1.0E-3*TOT_NP_O2_RATE
    !.. Efficiency for   N2(A3sig) + O -> O(1S) + N2
    RTS(60) = 0.37

    !.. N(2D) -> N(4S) + hv
    RTS(61) = 1.07E-5
    !.. hv(>600A) + N2 -> N(4S) + N   branching ratio
    RTS(62) = 0.5
    !.. hv(>600A) + N2 -> N(2D) + N   branching ratio
    RTS(63) = 0.4
    !.. hv(>600A) + N2 -> N(2P) + N   branching ratio
    RTS(64) = 0.1
    !.. N+ + O2 -> O2+ + N(2D)
    !.. Branching ratio from O'Keefe et al J. Chem. Phys. 1968
    RTS(65) = 0.15*TOT_NP_O2_RATE

    !.. N+ + O2 -> NO+ + O(1D)
    !.. Branching ratio from O'Keefe et al J. Chem. Phys. 1968
    RTS(66) = 0.36*TOT_NP_O2_RATE
    !.. hv(Scum-Runge) + O2 -> O(1S) + O   branching ratio
    RTS(67) = 0.001
    !.. Effic of O2(A3,DEL) + O -> O(1S)
    RTS(68) = 0.1
    !.. O(1D) + O -> O + O   Abreu et al. PSS, p1143, 1986
    RTS(69) = 6.47E-12*(TN/300)**0.14
    !.. hv + N2 -> N+(5S) -> 2143 A emission yield from the 2s sigma g state
    !.. of N2. This was taken as 0.6 the value of Cleary and Barth JGR 1987,
    !.. p13,635 because they did not double EUV below 250 A.
    RTS(70) = 0.06

    !.. hv + N2 -> N+(1D) -> 6584 A emission (guess)
    RTS(71) = 0.3
    !.. hv + N2 -> N+(1S) -> 5755 A emission (guess)
    RTS(72) = 0.03
    !.. efficiency of production of N(2P) from e + N2+ reaction
    RTS(73) = 0.08
    !.. Efficiency for production of O(1D) from N(2D) + O2 reaction
    !.. See Link and Swaminathan, PSS, 1992, page 699
    RTS(74) = 0.1    !??? check
    !.. He+ + O2 -> He + O2+
    RTS(75) = 9.2E-12

    !.. He+ + O2 -> He + O+(2D) + O(3P)
    RTS(76) = 2.37E-10
    !.. O2+ + N(2D) -> NO+ + O
    RTS(77) = 1.8E-10
    !.. O2+ + N(2D) -> N+ + O2
    RTS(78) = 8.65E-10
    !.. N2+ + N(4S) -> N+ + N2
    RTS(79) = 1.0E-11
    !.. N2+ + NO -> NO+ + N2
    RTS(80) = 3.6E-10

    !.. N+ + NO -> N(4S) + NO+
    RTS(81) = 4.72E-10*(300/TI)**0.24
    !.. N+ + NO -> N2+ + O
    RTS(82) = 8.33E-11*(300/TI)**0.24
    !.. O+(2D) + NO -> NO+ + O
    RTS(83) = 1.2E-9
    !.. O+(2D) + N -> N+ + O
    RTS(84) = 1.5E-10
    !.. O+(2P) + O2 -> O+ + O2  Fox
    RTS(85) = 1.3E-10

    !.. O+(2P) + O2 -> O2+ + O
    RTS(86) = 1.3E-10
    !.. O+(2P) + N -> O+ + N(2D)
    RTS(87) = 1.0E-11
    !.. O+(2P) + NO -> NO+ + O
    RTS(88) = 1.2E-9
    !.. H+ + O2 -> O2+ + H
    RTS(89) = 3.8E-9
    !.. O+(2D) + N2 -> NO+ + N  !.. Li et al. (1997).
    !.. From the ratio of the cross sections.
    !.. The branching ratio to O+(4S) + N2 not given by Li et al.
    RTS(90) = 2.5E-11

    !.. He+ + O2 -> He + O+(4S) + O
    RTS(91) = 2.39E-11
    !.. He+ + O2 -> He + O+(2P) + O
    RTS(92) = 6.04E-10
    !.. He+ + O2 -> He + O+(4S) + O(1D)
    RTS(93) = 4.6E-11
    !.. He+ + NO -> He + N+ + O
    RTS(94) = 1.35E-9
    !.. He+ + NO -> He + O+ + N
    RTS(95) = 1.0E-10

    !.. N(2P) + e -> N(4S) + e
    RTS(96) = 2.04E-10*(TE/300)**0.85
    !.. N(2P) + e -> N(2D) + e
    RTS(97) = 9.5E-9
    !.. O(1D) + e -> O(3P) + e
    RTS(98) = 2.87E-10*(TE/300)**0.91
    !.. N2+ + O -> O+ + N2  !.. McFarland et al.(1974)
    !.. From Fox and Sung 2001
    RTS(99) = 0.07*1.0E-10*(300/Ti)**0.23
    RETURN
end subroutine RATS
    
   
!.........................<PESIMP.FOR>......................3 APR 92
!...... This test driver demonstrates how to call the model and
!...... how to calculate electron heating rate and 3371 excitation rate.
subroutine SECIPRD(ALT, SZADEG, F107, F107A, TE, TN, OXN, O2N, N2N, XNE, N2APRD)

    INTEGER I, K, IK          !-- loop control variables
    INTEGER IDIM              !.. Array dimensions
    INTEGER IMAX              !.. Maximum PE energy
    PARAMETER (IDIM=501)
    REAL PEFLUX(IDIM)         !-- PE flux
    REAL SIGIT(3)             !-- Total electron impact cross sections
    REAL SIGEX(22)            !-- Electron impact cross sections for OX
    REAL ALT                  !-- ALT = altitude (km)  { 120 -> 500 }
    REAL SZADEG               !-- solar zenith angle {0 -> 90 degrees}
    REAL F107,  F107A         !-- F107 = Solar 10.7 cm flux
    REAL TE, TN               !-- electron,  neutral temperatures (K)
    REAL XN(3), OXN, O2N, N2N !-- XN,  O,  O2,  N2,  neutral densities  (cm-3)
    REAL XNE                  !-- electron density  (cm-3)
    REAL XN2D                 !-- N(2D) density for N(2D) + e -> 2.5 eV
    REAL XOP2D                !-- O+(2D) density for O+(2D) + e -> 3.3 eV
    REAL SPRD(3, 6)
    REAL SIGOX, SIGN2, SIGEE  !.. Total exciation cross sections for O,  N2,  O2
    REAL N2APRD               !.. Production of N2A
    REAL DE(IDIM), EV(IDIM)
    !.. various ionization and excitation rates by EUV and PE
    REAL EUVION, PEXCIT, PEPION, OTHPR1, OTHPR2
    COMMON/EUVPRD/EUVION(3,12),PEXCIT(3,12),PEPION(3,12),OTHPR1(6)&
        ,OTHPR2(6)
    DATA SPRD/.4,.56,.44, .4,.28,.44, .2,.06,.10, 0.,.05,.00, 0.,.05&
        ,.00, 0.0,0.0,0.02/
    DATA IMAX/0/              !.. Initialize IMAX Reset in FLXCAL

    !.. Transfer neutral densities to the density array
    XN(1) = OXN
    XN(2) = O2N
    XN(3) = N2N
    N2APRD = 0.0

    PEPION(1:3,1:6) = 1.0E-15
    PEXCIT(1:3,1:6) = 1.0E-15

    !.. Cannot calculate PE if no densities
    IF ( OXN < 1.0E5 .OR. N2N < 1.0E5 ) RETURN
    IF ( SZADEG > 105 ) RETURN

    !********************************************************************
    !.. Go and get the photoelectron fluxes
    XN2D  = 0     !.. N(2D) density for calculating N(2D) + e -> 2.5 eV
    XOP2D = 0     !.. O+(2D) density for calculating O+(2D) + e -> 3.3 eV
    CALL FLXCAL(IDIM, ALT, SZADEG, TE, TN, XN, XNE, XN2D,&
                XOP2D, PEFLUX, AFAC, IMAX, DE, EV)
    !***************************************************************

    !........ sample calculation of ion production rates.
    DO I=1,IMAX

        E = EV(I)
        CALL TXSION(E, SIGIT)                          !.. total ion XS
        CALL SIGEXS(E, TE, XNE, SIGOX, SIGN2, SIGEE)   !.. Total excitation XS
        CALL OXSIGS(E, SIGEX, SIGOX)                   !.. OX cross sections
        IF(E < 250) N2APRD = N2APRD + 0.22*PEFLUX(I)*SIGN2*XN(3)*DE(I) !.. N2(A) prod
        PEXCIT(1,1) = PEXCIT(1,1) + PEFLUX(I)*SIGEX(1)*XN(1)*DE(I)  !.. O(1D) prod
        PEXCIT(1,2) = PEXCIT(1,2) + PEFLUX(I)*SIGEX(2)*XN(1)*DE(I)  !.. O(1S) prod

        !.. Evaluate ionization branching ratios for O+
        CALL OXRAT(E,SPRD(1,1),SPRD(1,2),SPRD(1,3))
        !.. Calculate ion production rates
        DO K=1,6
            PEPION(1:3,K) = PEPION(1:3,K) &
                              + PEFLUX(I)*SIGIT(1:3)*XN(1:3)*SPRD(1:3,K)*DE(I)
        END DO
        EP = E + 17
        PEFLX = PEFLUX(I)/12.57
    END DO
    RETURN
end subroutine SECIPRD


!:::::::::::::::::::::::::: PHOTOELECTRON MODEL  ::::::::::::::::::::::::
!....... This subroutine evaluates the photoelectron flux using the concept
!.......  production frequencies developed by Phil Richards at Utah
!....... State University March 1984. It supercedes the model described in
!....... JGR, p2155, 1983. Contact EAST::CSPARA::RICHARDS on SPAN network
!------- Some minor updates in April 1992 indicated by C----
!....... I would appreciate any feedback on bugs or clarity and if it
!....... contributes substantially to a paper, I would appreciate the
!....... appropriate acknowledgement.
!......       **************** WARNING ****************
!...... This program is constructed to produce reasonable agreement with
!...... the Atmosphere Explorer-E PES fluxes of John Doering (Lee et al.
!...... PSS 1980, page 947). It will NOT give good fluxes if the EUV
!...... attenuation is greater than about a factor of 7 (AFAC < 0.14).
!...... The model accurately reproduces the measured fluxes very closely
!...... for the case in the test driver at 148 km SZA=53 when AFAC=0.19.
!...... You should compare the output against the Lee et al. 1980 fluxes
!...... periodically as a check. It is doubtful below 140km during the
!...... day and below 200km near sunset. Between 200km & 350km, it should
!...... be good for solar zenith angles < 90 degrees. Above 350 km there
!...... is considerable uncertainty due to neglect of transport but most
!...... models have similar uncertainties at high altitudes due to the
!...... uncertainty in the conjugate photoelectron flux, and the pitch
!...... angle distribution.
!
!------ ALT = altitude (km)  { 120 -> 500 }
!------ SZADEG = solar zenith angle  {0 -> 90 degrees ? }
!------ TE, TN = electron, neutral temperatures (K)
!------ XN, XNE = O, O2, N2, and electron densities  (cm-3)
!------ XN2D, XOP2D = N(2D) and O+(2D) densities for electron quenching
!------ (cm-3). You may put these to ZERO if not available.
!------ PEFLUX = photoelectron flux to be returned (eV cm2 sec)-1
!------ AFAC = the solar EUV attenuation warning flag
subroutine FLXCAL(IDIM, ALT, SZADEG, TE, TN, XN, XNE, XN2D,&
                  XOP2D, PEFLUX, AFAC, IMAX, DE, EV)

    INTEGER RDIM,IMAX
    REAL EMAX             !.. maximum photoelectron energy.
    PARAMETER (RDIM=84)
    REAL RJOX(RDIM), RJN2(RDIM), XN(3), COLUM(3), PEFLUX(IDIM),&
         DE(RDIM), DELTE(RDIM), EV(RDIM), EN(RDIM), UVFAC(59), EUV
    COMMON/SOL/UVFAC,EUV
    !.. photoelectron production frequencies by 1.0E9. Renormalized below
    !-- O production frequencies
    DATA RJOX/10*19,15,18,14,10,13,9,13,9,7,11,6,26,6,31,6,5,22,4,&
        4,5,4.04,2.56,1.9,2.28,2.12,0.96,0.24,0.14,0.14,0.1,0.1,0.1&
        ,0.1,0.1, 10*.05,10*.04,10*.01,10*.02/
    !-- N2 production frequencies
    DATA RJN2/6*40,43,35,35,28,29,21,25,19,19,13,19,16,12,11,7,18,8,&
        46,27,5,5,5,5,5,5.34,2.92,1.84,2.22,1.62,0.62,0.05,0.05,0.05,&
        0.05,0.05,0.05,0.05,0.044,10*.02,10*.01,10*.01,10*.01/
    !-- PE energy grid
    DATA EN/0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,&
        13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,&
        25.5,26.5,27.5,28.5,29.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,&
        67.5,72.5,77.5,82.5,87.5,92.5,97.5,105,115,125,135,145,155,&
        165,175,185,195,205,215,225,235,245,255,265,275,285,295,305,&
        315,325,335,345,355,365,375,385,395,405,415,425,435,445,455,&
        465,475,485,495/
    !-- PE energy steps
    DATA DELTE/30*1.0,14*5.0,40*10/
    DATA EMAX/286.0/          !..  Maximum PE energy

    SZA = SZADEG/57.29578   !.. convert solar zenith angle to radians
    IF(IMAX < 10) THEN
        !.. transfer energy grid to pass back
        DE(1:RDIM) = DELTE(1:RDIM)
        EV(1:RDIM) = EN(1:RDIM)
        DO IE=1,RDIM
            IF(EN(IE) < EMAX) IMAX=IE
        ENDDO
    ENDIF

    !.. 2.5eV production from electron quenching of N2D
    PN2D  = XN2D*XNE*6.0E-10*SQRT(TE/300.0)
    !.. 3.3eV production from electron quenching of O+(2D)
    POP2D = XOP2D*XNE*6.6E-8*SQRT(300./TE)
    CASEL = 0.0
    !.. evaluate the neutral column density  &&&&&&&&
    CALL SCOLUM(I, SZA, ALT*1.0E5, TN, XN, COLUM)

    !...... begin flux calculation loop............................
    DO IE=1,IMAX
        I = IMAX + 1 - IE
        IF(I < 1) exit

        PEFLUX(I) = 0.0

        !... evaluate energy of photon responsible for electron at energy EE
        EE = EV(I)
        EP = EE + 17
        IF(EE < 22) EP = 45
        IF(EE >= 22 .AND. EE < 28) EP = 41
        IF(EE >= 28 .AND. EE < 38) EP = 49

        !.. evaluate total photoionization cross sections for photon energy EP
        XSOXT = T_XS_OX(EP)         !.. New OX cross section
        XSO2T = 2.2*T_XS_OX(EP)     !.. O2 XS is 2.2* O XS
        XSN2T = T_XS_N2(EP)         !.. New N2 cross section

        !... evaluate EUV attenuation factor AFAC
        TAU  = COLUM(1)*XSOXT + COLUM(2)*XSO2T + COLUM(3)*XSN2T
        AFAC = EXP(-TAU)

        !..... low energy cascade production from O(1D) and N2* impact
        CASOX = 0.0
        IF(EE < 10) CASOX = PEFLUX(I+2)*SIGOX*XN(1)
        CASN2 = 0.0
        IF(EE < 6) CASN2 = PEFLUX(I+1)*SIGN2*XN(3)

        !..... cascade production from thermal electron degradation
        CASEL = 0.0
        IF(I < IMAX) CASEL = PEFLUX(I+1)*SIGEE*XNE

        !... Production from electron quenching of metastables
        EPN2D = 0.0
        IF(NINT(EE) == 3) EPN2D  = PN2D
        EPOP2D = 0.0
        IF(NINT(EE) == 4) EPOP2D = POP2D

        !.... evaluate cross sections (must be after cascade production)
        CALL SIGEXS(EE,TE,XNE,SIGOX,SIGN2,SIGEE)
        !..... adjust EUV production rate for different period of solar cycle
        CALL FACFLX(EE,UVFAC,FFAC)

        !..... Production of pe's at energy EE, taking into account
        !..... attenuation and EUV variation, and renormalize frequencies
        PRODOX = RJOX(I) * XN(1) * AFAC * FFAC * 1.0E-9
        PRODN2 = RJN2(I) * XN(3) * AFAC * FFAC * 1.0E-9
        !..... Sum all the production rates
        PROD = PRODOX + PRODN2 + CASEL + CASOX + CASN2 + EPN2D + EPOP2D
        !.... total loss through collisions
        RLOSS = SIGOX*XN(1) + SIGN2*XN(3) + SIGEE*XNE
        !....... evaluate photoelectron flux
        PEFLUX(I) = PROD/RLOSS
    end do
    RETURN
end subroutine FLXCAL


subroutine FACFLX(EE,UVFAC,FFAC)
    !....... solar UVFAC factors. Correspond to the first 9 wavelengths
    !....... TORR et al.[1979] GRL page 771 table 3. UVFAC(9) is for 304A
    REAL UVFAC(59)

    if (EE <= 30) then
        FFAC=(7*UVFAC(9)+UVFAC(8)+0.2*UVFAC(6))/8.2
    else if(EE <= 38) then
        FFAC = ( 2*UVFAC(7) + .5*UVFAC(5) )/2.5
    else if(EE <= 45) then
        FFAC = UVFAC(4)
    else if(EE <= 66) then
        FFAC = UVFAC(3)
    else if(EE <= 108) then
        FFAC = UVFAC(2)
    else
        FFAC = UVFAC(1)
    end if
    RETURN
end subroutine FACFLX


!..... Program for evaluating the total inelastic cross sections
subroutine SIGEXS(E,TE,XNE,SIGOX,SIGN2,SIGEE)
    !........ loss to thermal electrons ....
    ET    = 8.618E-5*TE
    SIGEE = ( 3.37E-12/E**0.94/XNE**0.03 )&
            *( ( E - ET )/( E - 0.53*ET ) )**2.36
    
    !...... cross section for o(1d)
    SIGO1D = 0.0
    IF(E > 1.96) SIGO1D = 4E-16*( 1 - 1.96/E )**2 / E

    !...... total excitation cross section for O excluding O(1D)
    IF(E < 25) then
        SIGO = ( 0.4*E - 5 )*1.4E-17
    else
        SIGO = 7.0E-17
    end if
    IF(SIGO < 0.0) SIGO=0.0
    
    !...... total excitation cross section for N2......
    if(E < 4.0) then
        SIGN2 = 5.0E-9*( 1 - 1.4/E )**9 * (1.4/E)**16
    else if (E <= 11.5) then
        SIGN2 = ( 15.5*E - 104.8 )*1.7E-18
    else
        SIGN2 = 1.4E-16
    end if
    IF(SIGN2 < 0.0) SIGN2=0.0
    !
    !........ total ionization cross sections from Keiffer and Dunn ....
    SIGION = 0.0
    AL = ALOG10(E)
    if(E < 50) then
        SIGION = 1.0E-16*( 0.068*E - 1.06 )
    else if(1.2 <= AL .AND. AL < 2.7) then
        SIGION = -3.6E-16*( AL - 1.2 )*( AL - 3 )
    else if(2.7 < AL) then
        SIGION = 1.2E-14*EXP(-AL*1.6)
    end if
    IF(SIGION <= 0.0) SIGION=0.0
    
    SIGOX = SIGO1D + SIGO + 0.5*SIGION
    SIGN2 = SIGN2 + SIGION
    RETURN
end subroutine SIGEXS


!::::::::::::::::::::::: TXSION ::::::::::::::::::::::::::::::::::
!..... total ionization cross sections for O, O2, and N2
!..... ionization cross sections keiffer and dunn ........
!..... The N2+ and O2+ cross sections were modified in April 99 to agree
!..... with the Schram et al. cross sections at high energies
subroutine TXSION(E,SIGIT)
    DIMENSION SIGIT(3)

    !... SIGTMP is used for N2+ and O2+ at the high energies
    SIGTMP = 1.0E-13*EXP( -2.303*ALOG10(E) )

    !... N2+ cross section
    SIGIT(3) = 0.0
    IF(E > 15.0)          SIGIT(3) = 1.42E-14 * (1 - 9.0/E)**7.1 * E**(-0.7)
    IF(SIGTMP < SIGIT(3)) SIGIT(3) = SIGTMP
    !... This correction to convert units to cm**2. Keiffer and Dunn page 10
    SIGIT(3) = SIGIT(3) * 0.87972

    !... O2+ cross section
    SIGIT(2) = 0.0
    IF(E > 12.0)          SIGIT(2) = 1.08E-14 * (1 - 7.0/E)**8.6 * E**(-0.65)
    IF(SIGTMP < SIGIT(2)) SIGIT(2) = SIGTMP
    !... This correction to convert units to cm**2. Keiffer and Dunn page 10
    SIGIT(2) = SIGIT(2) * 0.87972

    !... O+ cross section from Brook et al. J. Phys. B. Vol 11 p 3115, 1978
    SIGIT(1) = 0.0
    IF(E > 12.0) SIGIT(1) = 7.33E-15 * (1 - 2.0/E)**34.3 * E**(-0.7)
    RETURN
end subroutine TXSION


!:::::::::::::::::::::::::::::::::::: OXRAT ::::::::::::::::::::::::::::::::
!....... This subroutine returns the electron impact branching ratios
!....... for atomic oxygen from Burnett and Rountree Phys. Rev. A. 20
!....... 1979 page 1468
SUBROUTINE OXRAT(E,R4S,R2D,R2P)

    R4S = 1.0
    R2D = 0.0
    R2P = 0.0
    EV  = E
    IF(E > 100.0) EV  = 100.0
    IF(EV > 17)   R4S = -1.6E-3*EV + 0.56
    IF(EV > 17)   R2D = 1.067E-3*EV + 0.2933
    R2P = 1 - R4S - R2D

    IF(EV < 22) THEN
        R2P  = 0.0
        RTOT = R4S + R2D
        R4S  = R4S/RTOT
        R2D  = R2D/RTOT
    ENDIF
    RETURN
END SUBROUTINE OXRAT


!::::::::::::::::::::: T_XS_N2 :::::::::::::::::::::::::::
!.... This function calculates the N2 total photoionization
!.... cross section. P. Richards 2003-10-04
REAL FUNCTION T_XS_N2(EP)
    IMPLICIT NONE

    REAL EP   !... photon energy
    REAL ESAVE
    DATA ESAVE/0.0/

    IF(EP >= 600.0) THEN
        !.. Wavelength < 20 A, Auger ionization
        T_XS_N2 = 0.5E-18
    ELSE IF(EP >= 400.0) THEN
        !.. Wavelength < 31 A, Auger ionization
        T_XS_N2 = 1.0E-18
    ELSE IF(EP >= 392.0) THEN
        !.. Wavelength 31.62 to 23.70 A
        T_XS_N2 = EXP( 7.9864*ALOG(EP) - 91.6604 )
    ELSE IF(EP >= 55.09) THEN
        !.. Wavelength 225 to 125 A
        T_XS_N2 = EXP( -2.3711*ALOG(EP) - 29.8142 )
    ELSE
        !.. Wavelength > 225 A
        T_XS_N2 = EXP( -1.1077*ALOG(EP) - 34.8787 )
    ENDIF
    ESAVE = EP
    RETURN
END FUNCTION T_XS_N2


!::::::::::::::::::::: T_XS_OX :::::::::::::::::::::::::::
!.... This function calculates the OX total photoionization
!.... cross section. P. Richards 2003-10-04
!.... Samson and Pareek Phys. Rev. A, 31, 1470, 1985
REAL FUNCTION T_XS_OX(EP)
    IMPLICIT NONE
    REAL EP   !... photon energy
    REAL ESAVE
    DATA ESAVE/0.0/

    IF(EP >= 500.0) THEN
        !.. Wavelength shorter than 25 A, Auger ionization
        T_XS_OX = 0.5E-18
    ELSE IF(EP >= 165.26) THEN
        !.. Wavelength shorter than 75 A
        T_XS_OX = EXP( -2.5209*ALOG(EP) - 28.8855 )
    ELSE IF(EP >= 55.09) THEN
        !.. Wavelength between 78 and 256.26 A
        T_XS_OX = EXP( -1.7871*ALOG(EP) - 32.6335 )
    ELSE
        !.. Wavelength longer than 256.26 A
        T_XS_OX = EXP( -1.3077*ALOG(EP) - 34.5556 )
    ENDIF
    ESAVE = EP
    RETURN
END FUNCTION T_XS_OX


!:::::::::::::::::::::: OXSIGS :::::::::::::::::::::::::::::::::::::
!....... Inelastic cross sections for electron impact on atomic oxygen
!....... E=electron energy, SIGEX(22)=array of partial cross sections,
!....... SIGEXT=total excitation cross section, and S
subroutine OXSIGS(E,SIGEX,SIGEXT)
    DIMENSION SIGEX(22),SO1D(7)
    !..- CROSS SECTION FOR O(1D) - New Doering cross section from JGR
    !..- p19531, 1992. Increases production by a factor of 1.13
    DATA SO1D/0.0,0.0,15.0,30.0,44.0,54.0,38.0/
    
    SIGEX = 0.0

    !..... Old cross section of Henry
    IF(E > 1.96) SIGEX(1) = 4E-16 * ( 1 - 1.96/E )**2 / E

    !........ O(1S) cross section: may be double Shyn et al. JGR 1986, 13751
    IF(E > 4.17) SIGEX(2) = 6.54E-17 * ( 1 - SQRT( 4.17/E ) ) / E

    !....... 1304, 1027 A, Zipf and Erdman JGR 1985, 11088 include cascade.
    !....... Direct excitation is half for  1304 (Vaughan and Doering,
    !........ JGR 1986, 13755 and 1987 in press)
    IF(E >= 10) SIGEX(3) = 67.6E-17*( E - 10 )/E**2

    !....... 989 cross section from Doering 1987 (1/2 of Zipf)
    IF(E >= 14) SIGEX(4) = 7.0E-17*( 1 - 14/E )/SQRT(E)

    SIGEX(5) = 0.38*SIGEX(4)

    !....... O(5S) 1356 A Stone And Zipf Corrected By Zipf And Erdman 1985
    !..- reparameterized 1 May 92 using FITXS.FOR (PGR)
    IF(E > 10.0) SIGEX(6) = 4.867E-12 * ( 1.0 - 9.0/E )**2.67 / E**4.0

    SIGEXT = SUM(SIGEX(1:6))
    RETURN
end subroutine OXSIGS
    
   
!.................... RSPRIM.FOR ..................................
!.... This routine evaluates the ionization rates for photon impact
!.... It is based on a FLIP model routine that was modified in August
!.... 2009 for the chemical equilibrium model by P. richards.
SUBROUTINE PRIMPR(IJ, Z, ZOX, ZN2, ZO2, HE, SZA, TN, F107, F107A, N4S)
    IMPLICIT NONE

    INTEGER IVERT, I, IJ, IK, IPROBS, IS, K, L, LMAX, NNI, K1
    REAL EUVION, F107, F107A, F107SV, FNFAC, FREQLY, FREQSR, O2LYXS, &
         O2SRXS, TAUN, UVFAC, ZLAM, EUV, FLUXN, LAMAX, PEPION, PEXCIT, &
         SIGABS, SIGION, TPOT, ZFLUX
    REAL Z, ZOX, ZN2, ZO2, HE, SZA, TN, CHI, ZZ, TNJ, TAU, FLUX, HEPLS, &
         FBSBN, DISN, TAUGAM, FLUXG, ALTG, XNSIGF, DSPECT, GL, N4S
    REAL COLUM(3), OTHPR1, OTHPR2, OTHPR3(6)
    REAL COLUMN(3), XN(3), PROB(3,6,37), XSNPLS(37), FNITE(37), CLNITE(3)
    REAL TPROB(3,6,37)

    !-- common to hold the EUV and photoelectron production rates
    COMMON/EUVPRD/EUVION(3,12), PEXCIT(3,12), PEPION(3,12),&
                  OTHPR1(6), OTHPR2(6)
    COMMON/SIGS/ZFLUX(37), SIGABS(3,37), ZLAM(37), SIGION(3,37),&
                TPOT(3,10), NNI(3), LAMAX
    COMMON/SOL/UVFAC(59), EUV

    SAVE PROB, F107SV, TPROB    !.. Save values that are only calc once

    DATA LMAX/0/, F107SV/0.0/, IPROBS/0/
    !.. Fluxes for nighttime ion production in the 37 wavelength bins of
    !.. Torr et al GRL 1979. The fluxes are set to reproduce the production
    !.. rates in Strobel et al. PSS, p1027, 1980. Note that most bins are
    !.. set to zero and that the Strobel production rates are scaled by
    !.. FNFAC to stabilize the O+ solution below 200 km. Note also that
    !.. the wavelengths in FNITE go from largest (#3=HI) to smallest.
    DATA FNITE/9E5,0.0,9E5,2*0.0,9E6,13*0.0,3E5,8*0.0,3E5,8*0.0/
    DATA FNFAC/1.0/

    !.. UVFAC(58) is left over from FLIP routines for compatibility
    UVFAC(58) = -1.0
    IF( ABS( (F107 - F107SV)/F107 ) > 0.005) THEN
        !.. update UV flux factors
        CALL FACEUV(UVFAC, F107, F107A)
        CALL FACSR(UVFAC, F107, F107A)
        !.. call params to get solar flux data and cross sections
        CALL PARAMS(0, LMAX)
        F107SV = F107
    ENDIF

    !..  find probability for formation of each state  ........
    IF(IPROBS == 0) THEN
        CALL PROBS(0, PROB, ZLAM, LMAX, NNI)
        TPROB = PROB
        IPROBS=1
    ENDIF

    !... initialization of production rates. 1.0E-15 stabilizes
    !... e density evaluation at low altitudes in CMINOR
    EUVION = 1.0E-15
    DISN   = 0.0
    OTHPR1 = 1.0E-15
    OTHPR2 = 1.0E-15
        
    !........ Nighttime He+ production is calculated and stored. Attenuated to
    !........ avoid excess production at low altitudes
    OTHPR1(2) = 8E-11* EXP( -1.0E-11 * ZN2 ) * HE
    COLUM(I) = 1.0E+25

    TNJ   = TN
    XN(1) = ZOX
    XN(2) = ZO2
    XN(3) = ZN2
    ZZ    = Z * 1.0E+5
    CHI   = SZA
    
    !*****  obtain reaction rates from subr rats to get their densities
    !....... determine if sun is below the horizon ...
    !---- Must now do calculation for night production - Feb 93
    ALTG = ( 6371.0 + Z )*SIN(3.1416 - CHI) - 6371.0
    IF(Z > 1500) RETURN

    !...... get column densities for scattered light at night
    CALL SCOLUM(IJ, 0.0E0, ZZ, TNJ, XN, CLNITE)
    
    !...... evaluate the neutral column density  
    CALL SCOLUM(IJ, CHI, ZZ, TNJ, XN, COLUMN)

    !........ Store the column densities for the 2-Stream program
    COLUM = COLUMN
    
    !........ O2 dissociation by Schumann-Runge UV.
    !........ OTHPR1(3)= dissociation rate. OTHPR1(5)= Energy
    CALL SCHUMN(IJ, Z, ZO2, COLUMN, OTHPR1(3), OTHPR1(5))
    
    !---- Calculate hv + NO ion. freq. from Lyman-a (Brasseur & Solomon)
    !---- OTHPR2(2) is photodissociation of NO in the SR bands.
    !---- A small night production from scattered light is included. FREQLY
    !---- varies with solar activity using Richards et al. 1994 page 8981
    !---- LY_a=2.5E11 (Lean), sigi(NO)=2.0E-18 (Brasseur & Solomon page 329)
    DATA O2LYXS,O2SRXS,FREQSR /1.0E-20,1.0E-21,5.0E-6/
    FREQLY = 5.0E-7*( 1 + 4.0E-3*( 0.5*( F107 + F107A ) - 80.0 ) )
    OTHPR2(1) = FREQLY*( EXP( -O2LYXS*COLUMN(2) ) &
                         + 0.001*EXP( -O2LYXS*CLNITE(2) ) )
    OTHPR2(2) = FREQSR*( EXP( -O2SRXS*COLUMN(2) ) &
                         + 0.001*EXP( -O2SRXS*CLNITE(2) ) )
    
    HEPLS=0.0
    DO L=1,LMAX ! wavelength loop
        !..  TAU, TAUN = optical depth for day, night
        TAUN = SUM(SIGABS(:,L)*CLNITE)
        TAU  = SUM(SIGABS(:,L)*COLUMN)

        !.. evaluate nighttime flux and daytime flux
        FLUXN= FNFAC * (F107/75.0) * FNITE(L) * EXP(-TAUN)
        FLUX = ZFLUX(L)*EXP(-TAU) + FLUXN

        !.. he+ production. He+ X-S  = 0.25 N2  X-S. HEPRDN = nite He+
        IF(ZLAM(L) < 500.) HEPLS = HEPLS + HE*0.25*SIGION(3,L)*FLUX

        !-- hv + N -> N+ + e. ion. freq. Richards et al. JGR 1994 page 8989
        DATA XSNPLS/6*0.0,.211,10.294,11.171,10.961,11.244,11.323,12.098&
            ,13.265,12.423,11.951,11.212,11.798,11.758,11.778,11.772,11.503&
            ,11.016,10.578,9.556,8.15,8.302,7.298,6.413,6.399,5.192,5.725&
            ,4.787,3.778,2.3,.878,.286/

        OTHPR2(3) = OTHPR2(3) + XSNPLS(L)*1.0E-18*FLUX*N4S

        IF(ZLAM(L) >= 600.0) THEN
            !...... calculation of total euv absorption-ionization.....
            FBSBN = FLUX*( SIGABS(3,L) - SIGION(3,L) )*XN(3)
            !.. Save energy absorbed in the photodissociative process
            OTHPR1(4) = OTHPR1(4) + 1.24E+4*FBSBN/ZLAM(L)
            !.. production on atomic nitrogen by dissociation
            DISN = DISN + FBSBN

            !.. take into account the large n2 absorption of lyman gamma(972.54)
            IF(NINT(ZLAM(L)) == 975) THEN
                TAUGAM = 370E-18*COLUMN(3)
                IF(TAUGAM > 70.0) TAUGAM = 70.0
                FLUXG  = UVFAC(34) * 0.82E+9 * EXP(-TAUGAM)
                DISN   = DISN + FLUXG*370E-18*XN(3)
            ENDIF
        ENDIF

        !***** species loop begins here *****
        DO I=1,3
            XNSIGF = XN(I) * SIGION(I,L) * FLUX
            K1 = NNI(I)

            DO K=1,K1
                !.. # ions formed by w-l l by ionization of k state of species i
                DSPECT = XNSIGF * PROB(I,K,L)
                !.. store ion production rates .....
                EUVION(I,K)  = EUVION(I,K) + DSPECT
                !.. calculation of ion heating rate......
                EUVION(1,10) = EUVION(1,10) + DSPECT*TPOT(I,K)
            end do
        end do
    end do ! wavelength loop ends here
    
    !.........Store UV disoc of N2 2 atoms produced for every dissociation
    OTHPR1(1) = 2.0 * DISN

    !........ Transfer He+ production to storage
    OTHPR1(2) = OTHPR1(2) + HEPLS

    RETURN
end subroutine PRIMPR


!:::::::::::::::::::::::::::: SCOLUM ::::::::::::::::::::::::::::::::::
!.... this routine evaluates the neutral column density for O, O2, and N2
!.... see Smith & Smith JGR 1972 p 3592
!.... chi=solar zenith angle, RE & GE radius and grav constant for earth
!.... Modified by P. Richards January 2010 to eliminate need for calling
!.... the MSIS model at grazing incidence
subroutine SCOLUM(J,CHI,Z,TN,XN,COLUMN)
    IMPLICIT NONE

    ! Variable declarations
    INTEGER I, J
    REAL ZG, CHI, Z, TNJ, ALTG, GE, GR, RE, RP
    REAL SH, XP, Y, ERFY2, CHAPFN, RG, HG, XG, EM, F, G, A, B, C, D, GL
    REAL TN, XI, TINF, GTN
    REAL XN(3), COLUMN(3), SN(3), M(3), DG(9), T(2), GN(3)

    ! Variable Initialization
    DATA A, B, C, D, F, G&
        /1.0606963,0.55643831,1.0619896,1.7245609,0.56498823,0.06651874/
    DATA SN/0.0,0.0,0.0/
    DATA EM, M(1), M(2), M(3), RE, GE&
        / 1.662E-24 ,   16. ,  32. ,  28. ,6.357E8, 980/
    DATA T, ALTG, ERFY2/0.0,0.0,0.0D0,0.0D0/
    DATA DG/9*0.0/

    SN     = 0.0
    COLUMN = 1.E+25

    TNJ = TN     !.. Avoids changing Tn at grazing incidence
    IF(CHI >= 1.5708) then      !.. is sza>90.0 degrees
        !..Grazing incidence parameters
        ALTG = (6371.0E5 + Z) * SIN(3.1416 - CHI) - 6371.0E5
        IF(ALTG >= 85*1.0E5) THEN
            ZG = ALTG*1.E-5

            !.. Bates-Walker temperature
            XI   = (ZG - 120.0) * (6357.0 + 120.0) / (6357.0 + ZG)
            TINF = MAX(TN, 500.0)   !.. Crude TINF
            GTN  = TINF - ( TINF - 300.0 )*EXP(-0.025*XI)
            if(GTN < 180.0) GTN = 180.0

            !.. Neutral densities are extrapolated from altitude to grazing
            !.. altitude. Weighted average Tn and GTn is used
            GR = GE * ( RE/(RE + Z) )**2   !.. gravity
            GN = XN*EXP( 3*EM*M*GR*(Z - ALTG)/( 1.38E-16*(TN + 2*GTN) ) )

            !.. Make sure that grazing incidence density is not lower than MSIS
            !.. This is for using non MSIS densities like CTIPe
            TNJ = GTN
            SN = MAX(GN,XN) ! Note, this is per element of each array
        ELSE
            RETURN
        ENDIF
    end if

    !.. sn(1)=o , sn(2)=o2 , sn(3)=n2 , tnj=tn,  gr=gravity, rp=distance
    !.. to pt p, sh=scale height, rg=distance to pt g, hg=scale height at g
    GR = GE * ( RE/(RE+Z) )**2
    RP = RE + Z

    !.. Calculate column densities for each species
    DO I=1,3
        SH = (1.38E-16*TNJ)/(EM*M(I)*GR)
        XP = RP/SH
        Y  = SQRT(0.5*XP)*ABS(COS(CHI))

        IF(Y > 100.0) WRITE(6,100) I,Z/1.0E5,CHI*57.3,TNJ,EM,M(I),GR,RP
        IF(Y > 8) ERFY2 = F/(G + Y)
        IF(Y < 8) ERFY2 = (A + B*Y)/(C + D*Y + Y**2)
        IF(CHI <= 1.5708) then
            CHAPFN    = SQRT( 0.5*3.1416*XP )*ERFY2
            COLUMN(I) = XN(I)*SH*CHAPFN
            exit
        end if
        RG = RP*SIN(CHI)
        HG = 1.38E-16*TNJ/( EM*M(I)*980. ) * (1. + ALTG/6371.E5)**2
        XG = RG/HG
        COLUMN(I) = SQRT(0.5*3.1416*XG) * HG * ( 2.0*SN(I) - XN(I)*ERFY2 )
    end do

    RETURN

! format statements
100 FORMAT('WARNING, Y IN COLUMN(I) > 100',I4,1P,9E10.2)
end subroutine SCOLUM


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!........ this program determines the cross sections, solar fluxes, and
!........ given in m. torr et al g.r.l 1979 p771, table 2 and 3. but
!........ the longer wavelengths come first
subroutine PARAMS(ISW,LMAX)
    IMPLICIT NONE

    INTEGER  I, IN, IS, ISW, J, L, LAMAX, LMAX, NNI
    REAL EUV, FFAC,  SIGABS, SIGION, TPOT, UVFAC, ZFLUX, ZLAM
    REAL X1(111), X2(111), X3(18), ZLX(37), ZFX(37)
    COMMON/SIGS/ZFLUX(37), SIGABS(3, 37), ZLAM(37), SIGION(3, 37), &
                TPOT(3, 10), NNI(3), LAMAX
    COMMON/SOL/UVFAC(59), EUV
    
    !....... ionization potentials for o,o2 and n2 see kirby et al note the
    !....... o2 2(pi)u and 2(sigma-)u , and n2 f2(sigma)u pots are guesses
    !....... the sixth n2 potential is for dissociation
    DATA X3/13.6,16.9,18.6,28.5,40.,0.0,12.1,16.5,18.2,20.,&
        25.,0.,15.6,16.7,18.8,25.,29.,37./
    !........ wavelength data. average is taken for bands
    DATA ZLX/1025.,1031.91,1025.72,975.,977.02,925.,875.,825.,775.,&
        789.36,770.41,765.15,725.,703.36,675.,625.,629.73,609.76,575.,&
        584.33,554.31,525.,475.,465.22,425.,375.,368.07,325.,303.78,&
        303.31,275.,284.15,256.3,225.,175.,125.,75./
    !........ fluxes from table 3. these are for 74113. just replace this data
    !........ for other years in the table. note!!!! flux doubled for lambda<250
    !........ shortest wavelenghts have been tripled
    DATA ZFX/2.4665,2.1,3.5,1.4746,4.4,3.,3.537,1.625,.758,.702,&
        .26,.17,.141,.36,.23,.342,1.59,.53,.357,1.27,.72,.452,.285,&
        .29,.383,.314,.65,.965,6.9,.8,1.679,.21,.46,3.1,4.8,.45,1.2/
    !........ absorption cross sections -- o first ,o2, then n2
    DATA X1/5*0.0,1.315,4.554,3.498,5.091,3.749,3.89,4,10.736,11.46&
        ,17.245,13.365,13.4,13.4,13.024,13.09,12.59,12.059,12.127,11.93&
        ,11.496,9.687,9.84,8.693,7.7,7.68,6.461,7.08,6.05,5.202,3.732&
        ,1.839,.73,  1.346,1.0,1.63,21.108,18.73,12.817,8.562,16.631&
        ,22.145,26.668,18.91,20.8,28.535,27.44,21.919,26.017,32.06&
        ,28.07,26.61,22.79,26.04,24.606,23.101,21.91,20.31,18.118&
        ,18.32,17.438,16.81,16.8,14.387,15.79,13.37,10.9,7.509,3.806&
        ,1.316,  3*0.0,50.988,2.24,9.68,20.249,16.992,33.578,16.487&
        ,14.18,120.49,24.662,26.54,31.755,23.339,23.37,22.79,22.787&
        ,22.4,24.13,24.501,23.471,23.16,21.675,16.395,16.91,13.857&
        ,11.7,11.67,10.493,10.9,10.21,8.392,4.958,2.261,0.72/
    !....... ionization cross sections
    DATA X2/5*0.0,1.315,4.554,3.498,5.091,3.749,3.89,4,10.736,11.46&
        ,17.245,13.365,13.4,13.4,13.024,13.09,12.59,12.059,12.127,11.93&
        ,11.496,9.687,9.84,8.693,7.7,7.68,6.461,7.08,6.05,5.202,3.732&
        ,1.839,.73,  .259,0.0,1.05,13.94,15.54,9.374,5.494,6.413,10.597&
        ,10.191,8.47,11.72,23.805,23.75,21.306,24.937,31.1,26.39&
        ,26.61,22.79,26.04,24.606,23.101,21.91,20.31,18.118&
        ,18.32,17.438,16.81,16.8,14.387,15.79,13.37,10.9,7.509,3.806&
        ,1.316,  8*0.0,14.274,8.86,8.5,65.8,15.06,25.48,29.235&
        ,23.339,23.37,22.79,22.787&
        ,22.4,24.13,24.501,23.471,23.16,21.675,16.395,16.91,13.857&
        ,11.7,11.67,10.493,10.9,10.21,8.392,4.958,2.261,0.72/
    
    NNI(1) = 5
    NNI(2) = 5
    NNI(3) = 6
    LMAX   = 37
    IF(ISW /= 0) WRITE(17,95)
    
    DO L=1,LMAX
        ZLAM(L) = ZLX(L)
        FFAC    = UVFAC( LMAX + 1 - L )
        IF(ZFX(L) < 100) ZFLUX(L) = ZFX(L)*1.0E+9*FFAC

        !..- setting up ionization potentials
        IF(L <= 6)THEN
            TPOT(1,L) = X3(L)
            TPOT(2,L) = X3(6+L)
            TPOT(3,L) = X3(12+L)
        ENDIF

        !..- setting up cross sections
        DO IS=1,3
            IN = LMAX*(IS-1) + L
            SIGABS(IS,L) = X1(IN)*1.0E-18
            SIGION(IS,L) = X2(IN)*1.0E-18
            IF( SIGABS(IS,L) < SIGION(IS,L) ) SIGABS(IS,L) = SIGION(IS,L)
        end do
            
        IF(ISW /= 0) then
            WRITE(17,90) L,ZLAM(L),ZFLUX(L),(SIGABS(I,L),I=1,3),&
                         (SIGION(I,L),I=1,3),FFAC
        end if
    end do
    
    IF(ISW == 0) RETURN
    WRITE(17,94)
    WRITE(17,91) ((TPOT(I,J),J=1,6),I=1,3)
    
    RETURN

! format statements
90   FORMAT(1X,I4,F9.2,1P,22E9.2)
91   FORMAT(22F5.1)
94   FORMAT(/5X,' Ionization potentials for O, O2, N2',/&
            /2X,'4S   2D   2P   4P   2P*  -   X2   a+A  b4   B2   dis  -',&
            '  X2   A2   B2   C2   F2   2s')
95   FORMAT(/5X,'EUV fluxes, Photoabsorption, and Photoionization ',&
            'Cross sections',&
            /4X,'I',5X,'lam',5X,'flux',4X,'sigaOX',3X,'sigaO2',&
            3X,'sigaN2',3X,'sigiOX',3X,'sigiO2',3X,'sigiN2',3X,'UVfac')
end subroutine PARAMS


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!.... program for finding branching ratios (probabilities for various ion
!.... and molecular states) of o,o2,n2
!.... ---refs--- m. torr et al grl 1979 page 771, kirby et al atomic data
!.... and nuclear tables 1979 23,page 63
subroutine PROBS(ISW,PROB,ZLAM,LMAX,NNI)
    IMPLICIT NONE

    INTEGER I,IS,ISW,J,L,LL,LMAX,NNI(3)
    REAL YO(37,5),PROB(3,6,37),ZLAM(37),YSUM

    !...... coefficients of o ionization cross sections from torr et al
    !..... table 2
    DATA YO/.19,.486,.952,1.311,1.539,1.77,1.628,1.92,1.925,2.259&
        ,2.559,2.523,3.073,3.34,3.394,3.421,3.65,3.92,3.62,3.61,3.88,4.25&
        ,5.128,4.89,6.739,4.0,3.89,3.749,5.091,3.498,4.554,1.315,5*0.0&
        ,.206,.529,1.171,1.762,2.138,2.62,2.325,2.842,2.849,3.446,3.936&
        ,3.883,4.896,5.37,5.459,5.427,5.67,6.02,5.91,6.17,6.29,6.159&
        ,11.453,6.57,3.997,12*0.0, .134,.345,.768,1.144,1.363,1.63&
        ,1.488,1.92,1.925,2.173,2.558,2.422,2.986,3.22,3.274,3.211,3.27&
        ,3.15,3.494,3.62,3.23,2.956,0.664,14*0.0,  .062,.163,.348,.508&
        ,.598,.71,.637,.691,.693,.815,.787,.859,.541,24*0.0, .049,.13&
        ,.278,.366,.412,.35,.383,.307,.308,28*0.0/

    !....... production of o states from torr et al table 2 (yo array)
    !....... need to reverse order of yo to correspond with lambda
    DO L=1,LMAX
        LL = LMAX + 1 - L
        YSUM = SUM(YO(LL,:))
        DO I=1,5
            PROB(1,I,L)=0.0
            IF(YSUM /= 0.0) PROB(1,I,L) = YO(LL,I)/YSUM
        end do
    end do
            
    !....... call separate subroutines for o2 and n2 probabilities
    DO L=1,LMAX
        CALL PROBO2(1,L,ZLAM(L),PROB,NNI(2))
        CALL PROBN2(1,L,ZLAM(L),PROB,NNI(3))
    end do
        
    IF(ISW == 0) RETURN

    WRITE(17,95)
    DO L=1,LMAX
        WRITE(17,90) ZLAM(L),((PROB(IS,J,L),J=1,6),IS=1,3)
    end do

    RETURN

! format statements
90   FORMAT(F8.2,22F5.2)
95   FORMAT(/5X,' Photoionization branching ratios for O, O2, N2'&
    ,/3X,'Lam    4S   2D   2P   4P   2P*   -   X2   a+A  b4   B2 '&
    ,'  dis   -  X2   A2   B2   C2   F2   2s')
end subroutine PROBS


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!...... the n2 probabilities are taken from kirby et al tables b and c
!...... the yield of n+ is determined first then the remaining portion
!...... of the cross section is distributed amongst the n2+ ion states
!...... (x,a,b). the dissociation yield is divided between the 3 higher
!...... energy states according to wight et al. j.phys. b, 1976
!...... the 2 other states of kirby et al are not included
SUBROUTINE PROBN2(ISW,L,ZLAM,PROB,JPTS)
    IMPLICIT NONE

    INTEGER I,IPTS,ISW,J,JPTS,L
    REAL A(6),B(6),X(14),Y(14,6),PROB(3,6,37),PSUM,YIELD,YLAM,ZLAM

    DATA IPTS/14/
    DATA X/50.,210.,240.,280.,300.,332.,428.,500.,600.,660.,660.01,&
        720.,747.,796./
    DATA Y/5*.32,.3,.46,.404,.308,.308,.308,.42,1.,1.,5*.55,&
        .52,.46,.506,2*.589,.692,.58,2*.0,5*.13,.12,.08,&
        .09,.103,.103,4*0.0,3*0.0,.05,.1,.15,.83,1.,6*.0,3*.0,&
        .3,.4,.79,.17,7*.0,3*1.,.65,.5,.06,8*.0/
    
    !...... if zlam is too big set equal to x(max)
    YLAM = ZLAM

    !.. Prevent divide by zero
    IF(ZLAM > X(14)) YLAM = X(14) - 1
    IF(ZLAM < X(1))  YLAM = X(1)  + 1
    YIELD=0.0

    !...... determine yield of n+, and store in prob array
    CALL YLDISS(1, YLAM, YIELD)
    !
    DO I=1,IPTS
        ! kjh 6/22/92   NOTE:  I realize the following statement is strange
        !   looking, but its purpose is to prevent the CRAY compiler from
        !   vectorizing this loop.  (Which it does incorrectly).
        if(i == 25)write(6,*)' '
        IF(X(I) < YLAM .AND. YLAM <= X(I+1))  exit
    end do

    PSUM = 0.0
    !...... fit straight line between points
    A(1:JPTS) = ( Y(I+1,1:JPTS) - Y(I,1:JPTS) )/( X(I+1) - X(I) )
    B(1:JPTS) = Y(I,1:JPTS) - A(1:JPTS)*X(I)

    !...... determine probabilities of n2+ states
    DO J=1,JPTS
        IF(J <= 3) PROB(3,J,L) = ( A(J)*YLAM + B(J) )*(1-YIELD)
        IF(J > 3)  PROB(3,J,L) = ( A(J)*YLAM + B(J) )*YIELD
        PSUM = PSUM + PROB(3,J,L)
    END DO
        
    IF(PSUM == 0.0) RETURN

    !....... normalise probabilities
    PROB(3,1:JPTS,L) = PROB(3,1:JPTS,L) / PSUM
    RETURN
end subroutine PROBN2

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!..... determination of dissociative yield of n+, refer to kirby et al
!..... page 66 and table b
subroutine YLDISS(ISW,ZLAM,YIELD)
    IMPLICIT NONE
    INTEGER ISW,I,IPTS
    REAL X(9),Y(9),ZLAM,YIELD
    DATA IPTS/9/
    DATA X/50.,210.,240.,302.,387.,477.,496.,509.,2000./
    DATA Y/.36,.36,.346,.202,.033,.041,.024,0.0,0.0/
    
   
    IF(387 < ZLAM .AND. ZLAM < 477) then
        !...... parabolic interpolation see formula page 66 kirby et al
        YIELD = .0329 + 8.13E-6*( ZLAM - 442 )**2
    else
        !....... linear interpolation
        DO I=1,IPTS
            ! kjh 6/22/92   NOTE:  I realize the following statement is strange
            !   looking, but its purpose is to prevent the CRAY compiler from
            !   vectorizing this loop.  (Which it does incorrectly).
            if(i == 25)write(6,*)' '
            IF( X(I) <= ZLAM .AND. ZLAM < X(I+1) ) exit
        end do

        YIELD = ( ZLAM - X(I) )/( X(I+1) - X(I) )*( Y(I+1) - Y(I) ) + Y(I)
    end if
    RETURN
end subroutine YLDISS


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!....... o2 branching ratios are taken from kirby et al table d
!....... columns 4 & 9 are combined. columns 5,6,7&8 are combined
subroutine PROBO2(ISW,L,ZLAM,PROB,JPTS)
    IMPLICIT NONE

    INTEGER ISW,I,IPTS,J,JPTS,L
    REAL A(5),B(5),X(20),Y(20,5),PROB(3,6,37),BRSUM,YLAM,ZLAM
    DATA IPTS,X/20,304.,323.,454.,461.,504.,537.,556.,573.,584.,598.&
        ,610.,637.,645.,662.,684.,704.,720.,737.,774.,1026./
    DATA Y/.365,.374,.432,.435,.384,.345,.356,.365,.306,.23,.235,.245,&
        .34,.27,.482,.675,.565,.565,1.,1.,.205,.21,.243,.245,.27,.29,&
        .23,.27,.33,.295,.385,.35,.305,.385,.518,.325,.435,.435,2*.0,&
        .125,.124,.12,.12,.126,.13,.225,.216,.21,.375,.305,.37,.33,.345,&
        6*.0,.055,.167,.11,.105,.194,.234,.189,.149,.155,.103,.075,.036&
        ,.025,7*0.,.25,.125,.095,.95,.026,15*0./
    
    !...... if zlam is too big set equal to x(max)
    YLAM = ZLAM

    !...... if zlam is outside range of data values set equal to max or min
    IF(ZLAM > X(20)) then
        YLAM = X(20)
    else if(ZLAM <= X(1)) then
        YLAM = X(1)+1.E-3
    end if
    
    DO I=1,IPTS
        ! kjh 6/22/92   NOTE:  I realize the following statement is strange
        !   looking, but its purpose is to prevent the CRAY compiler from
        !   vectorizing this loop.  (Which it does incorrectly).
        if(i == 25)write(6,*)' '
        IF(X(I) < YLAM .AND.YLAM <= X(I+1))  exit
    end do

    A(1:JPTS) = ( Y(I+1,1:JPTS) - Y(I,1:JPTS) )/( X(I+1) - X(I) )
    B(1:JPTS) = Y(I,1:JPTS) - A(1:JPTS)*X(I)
    BRSUM = SUM(A(1:JPTS)*YLAM + B(1:JPTS))
       
    PROB(2,1:JPTS,L) = ( A(1:JPTS)*YLAM + B(1:JPTS) )/BRSUM
            
    RETURN
end subroutine PROBO2


!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!......... production of o(1d) by schumann-runge bands
!......... The fluxes are from Torr et al. GRL 1980 p6063. Scaling is
!......... done using UVFAC which may be set according to F10.7 cm flux
!......... may be done in FACEUV
subroutine SCHUMN(J,Z,ZO2,COLUMN,SCHUPR,SCHUHT)
    IMPLICIT NONE

    INTEGER J,JTI,LMAX,LSR
    REAL Z,ZO2,SCHUPR,SCHUHT,HSRX,FLD,EUV,SRXSCT,UVFAC
    COMMON/SOL/UVFAC(59),EUV
    REAL COLUMN(3),SRFLUX(8),SRXS(8),SRLAM(8)

    DATA SRFLUX/2.4,1.4,.63,.44,.33,.17,.12,.053/
    DATA SRXS/.5,1.5,3.4,6,10,13,15,12/
    DATA SRLAM/1725,1675,1625,1575,1525,1475,1425,1375/
    
    !........ lmax=# of lambdas in sub. primpr: schuht=heating: schupr=o(1d) prod
    LMAX=37

    DO LSR=1,8
        !......... photoabsorption cross section
        SRXSCT = 1.0E-18 * SRXS(LSR)
        HSRX   = SRXSCT * COLUMN(2)
        IF(HSRX > 70) HSRX = 70

        !........ attentuated solar flux
        FLD = UVFAC(LMAX + LSR) * 1.E+11 * SRFLUX(LSR) * EXP(-HSRX)

        !............ neutral heating SCHUHT and photodissociation rate SCHUPR
        SCHUHT = SCHUHT + 1.24E+4*FLD*SRXSCT*ZO2/SRLAM(LSR)
        SCHUPR = SCHUPR + FLD*SRXSCT
    end do
        
    SCHUPR = ZO2*SCHUPR
    JTI = 1

    RETURN

end subroutine SCHUMN


!:::::::::::::::::::::::::::::::: FACEUV :::::::::::::::::::::::
!----- This routine uses the EUV scaling from Richards et al.[1994]
!----- The EUVAC flux model is based on the F74113 solar reference
!----- spectrum and Hinteregger's scaling factors. This subroutine
!----- just provides the scaling factors as a function of the proxy
!----- (F107+F107A)/2
subroutine FACEUV(UVFAC,F107,F107A)
    IMPLICIT NONE

    INTEGER I
    REAL UVFAC(59),HFG200(37),A,B,C,F107,F107A,F107AV

    DATA HFG200/2.202,1.855,2.605,3.334,1.333,17.522,4.176,4.0&
        ,1.4,3.694,1.791,5.385,1.889,1.899,3.427,2.051,1.392,1.619&
        ,1.439,2.941,1.399,2.416,1.512,1.365,1.570,1.462,2.537,1.393&
        ,1.572,1.578,1.681,1.598,1.473,1.530,1.622,1.634,1.525/

    !--  Test to see if need to scale - see DATRD2 subroutine
    IF(NINT(UVFAC(58)) == -1.OR.NINT(UVFAC(58)) == -3) THEN
        !........... EUV scaling
        F107AV = (F107 + F107A)*0.5
        DO I=1,37
            A = (HFG200(I) - 1)/120.0
            B = 1 - A*80.0
            UVFAC(I) = A*F107AV + B
            IF(UVFAC(I) < 0.8) UVFAC(I)=0.8
        end do
    ENDIF
    RETURN
end subroutine FACEUV


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!........ The Schumann-Runge factors are scaled according to F10.7
!........ from Torr et al. GRL 1980 p6063
subroutine FACSR(UVFAC,F107,F107A)
    IMPLICIT NONE

    INTEGER I,LSR
    REAL UVFAC(59),SRFLUX(8),SRA(8),SRB(8),F107,F107A

    !............. Schumann-Runge scaling
    DATA SRFLUX/2.4,1.4,.63,.44,.33,.17,.12,.053/
    !...... first two SRA and SRB values out of order in Marsha's paper
    DATA SRA/25.5,20.7,13.2,11.6,11.3,7.86,7.68,4.56/
    DATA SRB/222.,129.,53.4,36.0,25.0,11.3,6.35,2.05/
    
    DO I=38,50
        LSR = I - 37
        UVFAC(I) = 1.0
        IF(LSR <= 8)then
            UVFAC(I) = SRA(LSR)*1.0E7*F107 + SRB(LSR)*1.0E9
            UVFAC(I) = 1.0E11*UVFAC(I)/SRFLUX(LSR)
        end if
    end do
    RETURN
end subroutine FACSR
