! irisub.for, version number can be found at the end of this comment.
!-----------------------------------------------------------------------        
! Includes subroutines IRI_SUB and IRI_WEB to compute IRI parameters 
! for specified location, date, time, and altitude range and subroutine 
! IRI_WEB to computes IRI parameters for specified location, date, time 
! and variable range; variable can be altitude, latitude, longitude, 
! year, month, day of month, day of year, or hour (UT or LT). 
! IRI_WEB requires IRI_SUB. Both subroutines require linking with the 
! following library files IRIFUN.FOR, IRITEC.FOR, IRIDREG.FOR, 
! CIRA.FOR, IGRF.FOR
!-----------------------------------------------------------------------        
! Programs using subroutine IRI_SUB need to include (see IRITEST.FOR):
!
!        call read_ig_rz
!       call readapf107
!
! Programs using subroutineIRI_WEB need to include (see IRITEST.FOR):
!
!       do i=1,100
!          oar(i,1)=-1.0
!          enddo
!
!        
!-----------------------------------------------------------------------        
! Required i/o units:  
!  KONSOL= 6 IRISUB: Program messages (used when jf(12)=.true. -> konsol)
!  IUCCIR=10 IRISUB: CCIR and URSI coefficients (CCIR%%.ASC, %%=month+10)
!  KONSOL=11 IRISUB: Program messages (used when jf(12)=.false. -> MESSAGES.TXT)
!    KONSOL=6/11 is also used in IRIFUN and IGRF. COMMON/iounit/konsol,mess 
!    is used to pass the value of KONSOL. If mess=false messages are turned off.
!  UNIT=12 IRIFUN/TCON:  Solar/ionospheric indices IG12, R12 (IG_RZ.DAT) 
!  UNIT=13 IRIFUN/APF..: Magnetic indices and F10.7 (APF107.DAT 
!  UNIT=14 IGRF/GETSHC:  IGRF coeff. (DGRF%%%%.DAT or IGRF%%%%.DAT, %%%%=year)
!-----------------------------------------------------------------------        

!*****************************************************************
!********* INTERNATIONAL REFERENCE IONOSPHERE (IRI). *************
!*****************************************************************

!-----------------------------------------------------------------
!
! INPUT:  JF(1:50)      true/false switches for several options
!         JMAG          =0 geographic   = 1 geomagnetic coordinates
!         ALATI,ALONG   LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
!         IYYYY         Year as YYYY, e.g. 1985
!         MMDD (-DDD)   DATE (OR DAY OF YEAR AS A NEGATIVE NUMBER)
!         DHOUR         LOCAL TIME (OR UNIVERSAL TIME + 25) IN DECIMAL 
!                          HOURS
!         HEIBEG,       HEIGHT RANGE IN KM; maximal 100 heights, i.e.
!          HEIEND,HEISTP        int((heiend-heibeg)/heistp)+1.le.100
!
!    JF switches to turn off/on (.true./.false.) several options
!
!    i       .true.                  .false.          standard version
!    -----------------------------------------------------------------
!    1    Ne computed            Ne not computed                     t
!    2    Te, Ti computed        Te, Ti not computed                 t
!    3    Ne & Ni computed       Ni not computed                     t
!    4    B0,B1 - Bil-2000       B0,B1 - other models jf(31)     false
!    5    foF2 - CCIR            foF2 - URSI                     false
!    6    Ni - DS-1995 & DY-1985 Ni - RBV-2010 & TTS-2005        false
!    7    Ne - Tops: f10.7<188   f10.7 unlimited                     t            
!    8    foF2 from model        foF2 or NmF2 - user input           t
!    9    hmF2 from model        hmF2 or M3000F2 - user input        t
!   10    Te - Standard          Te - Using Te/Ne correlation        t
!   11    Ne - Standard Profile  Ne - Lay-function formalism         t
!   12    Messages to unit 6     to messages.txt on unit 11          t
!   13    foF1 from model        foF1 or NmF1 - user input           t
!   14    hmF1 from model        hmF1 - user input (only Lay version)t
!   15    foE  from model        foE or NmE - user input             t
!   16    hmE  from model        hmE - user input                    t
!   17    Rz12 from file         Rz12 - user input                   t
!   18    IGRF dip, magbr, modip old FIELDG using POGO68/10 for 1973 t
!   19    F1 probability model   critical solar zenith angle (old)   t
!   20    standard F1            standard F1 plus L condition        t
!   21    ion drift computed     ion drift not computed          false
!   22    ion densities in %     ion densities in m-3                t
!   23    Te_tops (Bil-1985)     Te_topside (TBT-2012)           false
!   24    D-region: IRI-1990     FT-2001 and DRS-1995                t
!   25    F107D from APF107.DAT  F107D user input (oarr(41))         t
!   26    foF2 storm model       no storm updating                   t
!   27    IG12 from file         IG12 - user                         t
!   28    spread-F probability      not computed                    false
!   29    IRI01-topside          new options as def. by JF(30)   false
!   30    IRI01-topside corr.    NeQuick topside model            false 
! (29,30) = (t,t) IRIold, (f,t) IRIcor, (f,f) NeQuick
!   31    B0,B1 ABT-2009         B0 Gulyaeva-1987 h0.5               t   
! (4,31) = (t,t) Bil-00, (f,t) ABT-09, (f,f) Gul-87, (t,f) not used
!   32    F10.7_81 from file     F10.7_81 - user input (oarr(46))    t
!   33    Auroral boundary model on/off  true/false                 false
!   34    Messages on            Messages off                        t
!   35    foE storm model        no foE storm updating           false
!   36    hmF2 w/out foF2_storm  with foF2-storm                     t
!   37    topside w/out foF2-storm  with foF2-storm                  t
!   38    turn WRITEs off in IRIFLIP   turn WRITEs on                t
!   39    hmF2 (M3000F2)         new models                      false
!   40    hmF2 AMTB-model        Shubin-COSMIC model                 t
!   41    Use COV=F10.7_12       COV=f(IG12) (IRI before Oct 2015)   t
!   42    Te with PF10.7 dep.     w/o PF10.7 dependance               t
!      ....
!   50    
!   ------------------------------------------------------------------
!
!  Depending on the jf() settings additional INPUT parameters may 
!  be required:
!
!       Setting              INPUT parameter
!    -----------------------------------------------------------------
!    jf(8)  =.false.     OARR(1)=user input for foF2/MHz or NmF2/m-3
!    jf(9)  =.false.     OARR(2)=user input for hmF2/km or M(3000)F2
!    jf(10 )=.false.     OARR(15),OARR(16)=user input for Ne(300km),
!       Ne(400km)/m-3. Use OARR()=-1 if one of these values is not 
!       available. If jf(23)=.false. then Ne(300km), Ne(550km)/m-3.
!    jf(13) =.false.     OARR(3)=user input for foF1/MHz or NmF1/m-3 
!    jf(14) =.false.     OARR(4)=user input for hmF1/km
!    jf(15) =.false.     OARR(5)=user input for foE/MHz or NmE/m-3 
!    jf(16) =.false.     OARR(6)=user input for hmE/km
!    jf(17) =.flase.     OARR(33)=user input for Rz12
!    jf(21) =.true.      OARR(41)=user input for daily F10.7 index
!    jf(23) =.false.     OARR(41)=user input for daily F10.7 index
!    jf(24) =.false.     OARR(41)=user input for daily F10.7 index
!          optional for jf(21:24); default is F10.7D=COV
!    jf(25) =.false.     OARR(41)=user input for daily F10.7 index
!          if oarr(41).le.0 then 12-month running mean is 
!          taken from internal file]
!    jf(27) =.false.     OARR(39)=user input for IG12
!    jf(28) =.true.      OARR(41)=user input for daily F10.7 index
!
!
!  OUTPUT:  OUTF(1:20,1:1000)
!               OUTF(1,*)  ELECTRON DENSITY/M-3
!               OUTF(2,*)  NEUTRAL TEMPERATURE/K
!               OUTF(3,*)  ION TEMPERATURE/K
!               OUTF(4,*)  ELECTRON TEMPERATURE/K
!               OUTF(5,*)  O+ ION DENSITY/% or /M-3 if jf(22)=f 
!               OUTF(6,*)  H+ ION DENSITY/% or /M-3 if jf(22)=f
!               OUTF(7,*)  HE+ ION DENSITY/% or /M-3 if jf(22)=f
!               OUTF(8,*)  O2+ ION DENSITY/% or /M-3 if jf(22)=f
!               OUTF(9,*)  NO+ ION DENSITY/% or /M-3 if jf(22)=f
!                 AND, IF JF(6)=.FALSE.:
!               OUTF(10,*)  CLUSTER IONS DEN/% or /M-3 if jf(22)=f
!               OUTF(11,*)  N+ ION DENSITY/% or /M-3 if jf(22)=f
!               OUTF(12,*)  
!               OUTF(13,*)  
!  if(jf(24)    OUTF(14,1:11) standard IRI-Ne for 60,65,..,110km 
!     =.false.)        12:22) Friedrich (FIRI) model at these heights 
!                      23:33) standard Danilov (SW=0, WA=0) 
!                      34:44) for minor Stratospheric Warming (SW=0.5) 
!                      45:55) for major Stratospheric Warming (SW=1) 
!                      56:66) weak Winter Anomaly (WA=0.5) conditions
!                      67:77) strong Winter Anomaly (WA=1) conditions
!               OUTF(15-20,*)  free
!
!            OARR(1:100)   ADDITIONAL OUTPUT PARAMETERS         
!
!      #OARR(1) = NMF2/M-3           #OARR(2) = HMF2/KM
!      #OARR(3) = NMF1/M-3           #OARR(4) = HMF1/KM
!      #OARR(5) = NME/M-3            #OARR(6) = HME/KM
!       OARR(7) = NMD/M-3             OARR(8) = HMD/KM
!       OARR(9) = HHALF/KM            OARR(10) = B0/KM
!       OARR(11) =VALLEY-BASE/M-3     OARR(12) = VALLEY-TOP/KM
!       OARR(13) = TE-PEAK/K          OARR(14) = TE-PEAK HEIGHT/KM
!      #OARR(15) = TE-MOD(300KM)     #OARR(16) = TE-MOD(400KM)/K
!       OARR(17) = TE-MOD(600KM)      OARR(18) = TE-MOD(1400KM)/K
!       OARR(19) = TE-MOD(3000KM)     OARR(20) = TE(120KM)=TN=TI/K
!       OARR(21) = TI-MOD(430KM)      OARR(22) = X/KM, WHERE TE=TI
!       OARR(23) = SOL ZENITH ANG/DEG OARR(24) = SUN DECLINATION/DEG
!       OARR(25) = DIP/deg            OARR(26) = DIP LATITUDE/deg
!       OARR(27) = MODIFIED DIP LAT.  OARR(28) = Geographic latitude
!       OARR(29) = sunrise/dec. hours OARR(30) = sunset/dec. hours
!       OARR(31) = ISEASON (1=spring) OARR(32) = Geographic longitude
!      #OARR(33) = Rz12               OARR(34) = Covington Index
!       OARR(35) = B1                 OARR(36) = M(3000)F2
!      $OARR(37) = TEC/m-2           $OARR(38) = TEC_top/TEC*100.
!      #OARR(39) = gind (IG12)        OARR(40) = F1 probability 
!      #OARR(41) = F10.7 daily        OARR(42) = c1 (F1 shape)
!       OARR(43) = daynr              OARR(44) = equatorial vertical 
!       OARR(45) = foF2_storm/foF2_quiet         ion drift in m/s
!      #OARR(46) = F10.7_81           OARR(47) = foE_storm/foE_quiet 
!       OARR(48) = spread-F probability          
!       OARR(49) = Geomag. latitude   OARR(50) = Geomag. longitude  
!       OARR(51) = ap at current time OARR(52) = daily ap
!       OARR(53) = invdip/degree      OARR(54) = MLT-Te
!       OARR(55) = CGM-latitude       OARR(56) = CGM-longitude
!       OARR(57) = CGM-MLT            OARR(58) = CGM lat eq. aurl bodry
!       OARR(59) = CGM-lati(MLT=0)    OARR(60) = CGM-lati for MLT=1
!       OARR(61) = CGM-lati(MLT=2)    OARR(62) = CGM-lati for MLT=3
!       OARR(63) = CGM-lati(MLT=4)    OARR(64) = CGM-lati for MLT=5
!       OARR(65) = CGM-lati(MLT=6)    OARR(66) = CGM-lati for MLT=7
!       OARR(67) = CGM-lati(MLT=8)    OARR(68) = CGM-lati for MLT=9
!       OARR(69) = CGM-lati(MLT=10)   OARR(70) = CGM-lati for MLT=11
!       OARR(71) = CGM-lati(MLT=12)   OARR(72) = CGM-lati for MLT=13
!       OARR(73) = CGM-lati(MLT=14)   OARR(74) = CGM-lati for MLT=15
!       OARR(75) = CGM-lati(MLT=16)   OARR(76) = CGM-lati for MLT=17
!       OARR(77) = CGM-lati(MLT=18)   OARR(78) = CGM-lati for MLT=19
!       OARR(79) = CGM-lati(MLT=20)   OARR(80) = CGM-lati for MLT=21
!       OARR(81) = CGM-lati(MLT=22)   OARR(82) = CGM-lati for MLT=23
!       OARR(83) = Kp at current time OARR(84) = magnetic declination 
!       OARR(85) = L-value            OARR(86) = dipole moment 
!                # INPUT as well as OUTPUT parameter
!                $ special for IRIWeb (only place-holders)
!-----------------------------------------------------------------------        
!*****************************************************************
!*** THE ALTITUDE LIMITS ARE:  LOWER (DAY/NIGHT)  UPPER        ***
!***     ELECTRON DENSITY         60/80 KM       1500 KM       ***
!***     TEMPERATURES               60 KM        2500/3000 KM  ***
!***     ION DENSITIES             100 KM        1500 KM       ***
!*****************************************************************
!*****************************************************************
!*********            INTERNALLY                    **************
!*********       ALL ANGLES ARE IN DEGREE           **************
!*********       ALL DENSITIES ARE IN M-3           **************
!*********       ALL ALTITUDES ARE IN KM            **************
!*********     ALL TEMPERATURES ARE IN KELVIN       **************
!*********     ALL TIMES ARE IN DECIMAL HOURS       **************
!*****************************************************************
!*****************************************************************
!*****************************************************************
SUBROUTINE IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR, &
                   HEIBEG,HEIEND,HEISTP,OUTF,OARR)

    logical, intent(in) :: JF(50)
    integer   :: DAYNR,DDO,DO2,SEASON,SEADAY
    real      :: LATI,LONGI,MO2,MO,MODIP,NMF2,MAGBR,INVDIP,IAPO,&
                 NMF1,NME,NMD,MM,MLAT,MLONG,NMF2S,NMES,INVDPC
    character :: FILNAM*12

    dimension :: ARIG(3),RZAR(3),F(3),E(4),XDELS(4),DNDS(4),&
        FF0(988),XM0(441),F2(13,76,2),FM3(9,49,2),ddens(5,11),&
        elg(7),FF0N(988),XM0N(441),F2N(13,76,2),FM3N(9,49,2),&
        INDAP(13),AMP(4),HXL(4),SCL(4),XSM(4),MM(5),DTI(4),AHH(7),&
        STTE(6),DTE(5),ATE(7),TEA(6),XNAR(2),param(2),OARR(100),&
        OUTF(20,1000),DDO(4),DO2(2),DION(7),&
        osfbr(25),D_MSIS(9),T_MSIS(2),IAPO(7),SWMI(25),ab_mlat(48),&
        DAT(11,4), PLA(4), PLO(4)

    logical   :: EXT,SCHALT,TECON(2),sam_mon,sam_yea,sam_ut,sam_date,&
        F1REG,FOF2IN,HMF2IN,URSIF2,LAYVER,DY,DREG,rzino,FOF1IN,&
        HMF1IN,FOEIN,HMEIN,RZIN,sam_doy,F1_OCPRO,F1_L_COND,NODEN,&
        NOTEM,NOION,TENEOP,OLD79,URSIFO,igin,igino,mess,&
        dnight,enight,fnight,TOPO,TOPC,fstorm_on,estorm_on,&
        f2_in, use_ursi, spwx_in, rightTime, readSpwx, readOldSpwx

    COMMON /CONST/UMR,PI  /const1/humr,dumr   /ARGEXP/ARGMAX &
        /IGRF1/ERA,AQUAD,BQUAD,DIMO&
        /BLOCK1/HMF2,NMF2,HMF1,F1REG  /BLOCK2/B0,B1,C1  &
        /BLOCK3/HZ,T,HST              /BLOCK4/HME,NME,HEF &
        /BLOCK5/ENIGHT,E              /BLOCK6/HMD,NMD,HDX&
        /BLOCK7/D1,XKK,FP30,FP3U,FP1,FP2&
        /BLOCK8/HS,TNHS,XSM,MM,DTI,MXSM       &
        /BLOTE/AHH,ATE1,STTE,DTE&
        /BLO10/BETA,ETA,DELTA,ZETA    /findRLAT/FLON,RYEAR   &
        /BLO11/B2TOP,TC3,itopn,alg10,hcor1       &
        /iounit/konsol,mess     /CSW/SW(25),ISW,SWC(25)&
        /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU

    EXTERNAL          XE1,XE2,XE3_1,XE4_1,XE5,XE6,FMODIP

    DATA icalls/0/

    save
            
    mess=jf(34)
    
    ! set switches for NRLMSIS00  
    ISW=0
    do KI=1,25
        SWMI(KI)=1.
    end do

    nummax=1000
    do KI=1,20
        do kk=1,nummax
            OUTF(KI,kk)=-1.
        end do
    end do
    
    ! oarr(1:6,15,16,33,39:41) is used for inputs
    do kind=7,14,1
        oarr(kind)=-1.
    end do
    do kind=17,32,1
        oarr(kind)=-1.
    end do
    do kind=34,38,1
        oarr(kind)=-1.
    end do
    oarr(40)=-1.
    do kind=42,100,1
        if(kind /= 46) oarr(kind)=-1.
    end do

    ! PROGRAM CONSTANTS AND INITIALIZATION
    if(icalls < 1) then
        ARGMAX  = 88.0
        pi      = ATAN(1.0)*4.
        UMR     = pi/180.
        humr    = pi/12.
        dumr    = pi/182.5
        ALOG2   = ALOG(2.)
        ALG10   = ALOG(10.)
        ALG100  = ALOG(100.)
        montho  = -1
        nmono   = -1
        iyearo  = -1
        idaynro = -1
        rzino   = .true.
        igino   = .true.
        ut0     = -1
        ursifo  = .true.
        ! Initialize parameters for COMMON/IGRF1/
        !   ERA        EARTH RADIUS (WGS-84: 6371.137 KM) 
        !   EREQU   MAJOR HALF AXIS FOR EARTH ELLIPSOID (6378.160 KM)
        !   ERPOL   MINOR HALF AXIS FOR EARTH ELLIPSOID (6356.775 KM)
        !   AQUAD   SQUARE OF MAJOR HALF AXIS FOR EARTH ELLIPSOID
        !   BQUAD   SQUARE OF MINOR HALF AXIS FOR EARTH ELLIPSOID
        !   EEXC    Eccentricity of Earth's orbit
        !   DIMO    Earth's dipole moment in Gauss 
        ! ERA, EREQU and ERPOL as recommended by the INTERNATIONAL 
        ! ASTRONOMICAL UNION .
        ERA     = 6371.2
        EREQU   = 6378.16
        ERPOL   = 6356.775
        AQUAD   = EREQU**2
        BQUAD   = ERPOL**2
        EEXC    = 0.01675
        dimo    = 0.311653
    end if

    numhei = 1 + int( abs(heiend-heibeg) / abs(heistp) )
    if(numhei > nummax) numhei=nummax

    ! NEW-GUL------------------------------
    Y05     = .6931473
    QF      = 1.
    h05top  = 0.
    ! NEW-GUL------------------------------

    ! Code inserted to aleviate block data problem for PC version.
    ! Thus avoiding DATA statement with parameters from COMMON block.
    XDELS(1)= 5.
    XDELS(2)= 5.
    XDELS(3)= 5.
    XDELS(4)= 10.
    DNDS(1) = .016
    DNDS(2) = .01
    DNDS(3) = .016
    DNDS(4) = .016
    DDO(1)  = 9
    DDO(2)  = 5
    DDO(3)  = 5
    DDO(4)  = 25
    DO2(1)  = 5
    DO2(2)  = 5
    XNAR(1) = 0.0
    XNAR(2) = 0.0
    DTE(1)  = 5.
    DTE(2)  = 5.
    DTE(3)  = 10.
    DTE(4)  = 20.
    DTE(5)  = 20.
    DTI(1)  = 10.
    DTI(2)  = 10.
    DTI(3)  = 20.
    DTI(4)  = 20.

    ! FIRST SPECIFY YOUR COMPUTERS CHANNEL NUMBERS ....................
    ! AGNR=OUTPUT (OUTPUT IS DISPLAYED OR STORED IN FILE OUTPUT.IRI)...
    ! IUCCIR=UNIT NUMBER FOR CCIR COEFFICIENTS ........................
    IUCCIR  = 10

    !-web- special for web version
    !-web- messages should be turned off with mess=jf(34)=.false. 
    KONSOL  = 6
    if( .not. jf(12) .and. mess) then
        konsol=11
        open(11,file='messages.txt')
    end if

    ! selection of density, temperature and ion composition options ......
    NODEN   = (.not. jf(1))
    NOTEM   = (.not. jf(2))
    NOION   = (.not. jf(3))
    if ( .not. NOION ) NODEN = .false.
    DY      = (.not. jf(6))
    LAYVER  = (.not. jf(11))
    OLD79   = (.not. jf(7))
    F1_OCPRO= jf(19)
    F1_L_COND=(.not. jf(20))
    DREG    = jf(24)
    TOPO    = jf(29)
    TOPC    = jf(30)

    ! rz12, IG12, F10.7D, PF10.7 input option ............................
    RZIN    = (.not. jf(17))
    if (RZIN) then
      ARZIN = OARR(33)
    else
      oarr(33) = -1.
    end if
  
    IGIN    = (.not. jf(27))
    if(IGIN) then
      AIGIN=OARR(39)
    else
      oarr(39)=-1.
    end if

    if (.not. jf(25)) then
        f107din=OARR(41)
    else
        oarr(41)=-1.
    endif

    if (.not. jf(32)) then
        f10781in=OARR(46)
    else
        oarr(46)=-1.
    end if

    ! Topside density ....................................................
    if (TOPO) then
        if (TOPC) then 
            itopn=0
        else
            itopn=3
        endif
    else 
        if (TOPC) then
            itopn=1
        else
            itopn=2
        endif
    endif

    ! F2 peak density ....................................................

    FOF2IN  = (.not. jf(8))
    if (FOF2IN) then
        OARR1 = OARR(1)
        AFOF2 = OARR1
        ANMF2 = OARR1
        if(OARR1 <  100.) ANMF2 = 1.24E10 * AFOF2 * AFOF2
        if(OARR1 >= 100.) AFOF2 = SQRT( ANMF2 / 1.24E10 )
    else
        oarr(1)=-1.
    end if
    URSIF2  = (.not. jf(5))

    ! F2 peak altitude ..................................................
    HMF2IN  = (.not. jf(9))
    if (HMF2IN) then
        AHMF2 = OARR(2)
    else
        oarr(2) = -1.
    endif

    ! F1 peak density ...................................................
    FOF1IN  = (.not. jf(13))
    if (FOF1IN) then
        OARR3 = OARR(3)
        AFOF1 = OARR3
        ANMF1 = OARR3
        if (OARR3 <  100.) ANMF1 = 1.24E10 * AFOF1 * AFOF1
        if (OARR3 >= 100.) AFOF1 = SQRT( ANMF1 / 1.24E10 )
    else
        oarr(3)=-1.
    end if

    ! F1 peak altitude ..................................................
    HMF1IN  = (.not. jf(14))
    if (HMF1IN) then
        AHMF1=OARR(4)
        if( .not. layver .and. mess) write(konsol,1939)
    else
        oarr(4)=-1.
    end if

    ! E peak density ....................................................
    FOEIN   = (.not. jf(15))
    if(FOEIN) then
        OARR5 = OARR(5)
        AFOE  = OARR5
        ANME  = OARR5
        if (OARR5 <  100.) ANME=1.24E10*AFOE*AFOE
        if (OARR5 >= 100.) AFOE=SQRT(ANME/1.24E10)
    else
        oarr(5)=-1.
    end if 

    ! E peak altitude ..................................................
    HMEIN   = (.not. jf(16))
    if(HMEIN) then
        AHME = OARR(6)
    else
        oarr(6) = -1.
    end if

    ! TE-NE MODEL OPTION ..............................................
    TENEOP  = (.not. jf(10))
    if (TENEOP) then
        do JXNAR=1,2
            XNAR(JXNAR)=OARR(JXNAR+14)
            TECON(JXNAR)=.FALSE. 
            if(XNAR(JXNAR) > 0.) TECON(JXNAR)=.TRUE. 
        end do
    else
        oarr(15)=-1.
        oarr(16)=-1.
    end if

    ! lists the selected options before starting the table
    if (icalls < 1 .and. mess) then
                            write (konsol, 2911) 
        if (.not. NODEN) then
            if (LAYVER)     write(konsol,9012) 
            if (OLD79)      write(konsol,9014) 
            if (itopn == 0) write(konsol,9207)
            if (itopn == 1) write(konsol,9204)
            if (itopn == 2) write(konsol,9205)
            if (itopn == 3) write(konsol,9206)
            if(FOF2IN) then
                            write(konsol,9015) 
            else
                if(URSIF2) then
                            write(konsol,9016) 
                else
                            write(konsol,9017) 
                end if
                if(HMF2IN)  write(konsol,9018) 
                if(fof1in)  write(konsol,9019) 
                if(HMF1IN.and.LAYVER)&
                            write(konsol,9021) 
                if(jf(4)) then
                            write(konsol,9214)
                else 
                    if (jf(31)) then
                            write(konsol,9216)
                    else
                            write(konsol,9215)
                    endif
                endif
                if(foein)   write(konsol,9022) 
                if(HMEIN)   write(konsol,9023) 
                if(F1_OCPRO) write(konsol,9024) 
                if(F1_L_COND) write(konsol,9025) 
                if(DREG) then 
                            write(konsol,9026) 
                else
                            write(konsol,9027) 
                endif
                if(jf(26)) then
                    if(fof2in) then 
                            write(konsol,9028) 
                            jf(26)=.false.
                    else
                            write(konsol,9029) 
                    endif
                endif
            end if
        end if

        if ( (.not. NOION) .and. (DY) )&
                            write(konsol,9031) 
        if ( (.not. NOION) .and. (.not. DY) )&
                            write(konsol,9039) 

        if (.not. NOTEM) then 
            if (TENEOP)     write(konsol,9032) 
            if (jf(23)) then 
                            write(konsol,9033) 
            else
                            write(konsol,9034) 
            endif

            if (jf(33)) then
                            write(konsol,4031)
            else
                            write(konsol,4032)
            endif

            if (jf(35)) then
                            write(konsol,9128)
            else
                            write(konsol,9129)
            endif
        end if
    end if

    ! CALCULATION OF DAY OF YEAR OR MONTH/DAY AND DECIMAL YEAR 
    ! NRDAYM is the number of days in the current month 
    ! IDAYY is the number of days in the current year
    !
    !  leap year rule: years evenly divisible by 4 are leap years, except
    !  years also evenly divisible by 100 are not leap years, except years 
    !  also evenly divisible by 400 are leap years. The year 2000 is a 100 
    !  and 400 year exception and therefore it is a normal leap year. 
    !  The next 100 year exception will be in the year 2100!

    iyear = iyyyy
    if(iyear < 100) iyear = iyear + 1900
    if(iyear < 30 ) iyear = iyear + 2000
    idayy = 365
    if ( MOD(iyear, 4) == 0 ) then ! leap year
        if ( MOD(iyear, 100) /= 0 .or. MOD(iyear, 400) == 0) idayy = 366
    end if

    if (MMDD < 0) then
        DAYNR=-MMDD
        call MODA(1,iyear,MONTH,IDAY,DAYNR,nrdaym)
    else
        MONTH=MMDD/100
        IDAY=MMDD-MONTH*100
        call MODA(0,iyear,MONTH,IDAY,DAYNR,nrdaym)
    end if

    ryear = iyear + (daynr-1.0)/idayy
    iyd = iyear *1000 + daynr
    amx = pi*(daynr-3.)/182.6
    radj = 1.-eexc*(cos(amx)+eexc*(cos(2*amx)-1.)/2.)

    ! calculate center height for CGM computation
    height_center=(HEIBEG+HEIEND)/2.

    ! CALCULATION OF GEODETIC/GEOMAGNETIC COORDINATES (LATI, LONGI AND 
    ! MLAT, MLONG), MAGNETIC INCLINATION (DIP), DIP LATITUDE (MAGBR) 
    ! AND MODIFIED DIP (MODIP), ALL IN DEGREES
    if(along < 0.) along = along + 360. ! -180/180 to 0-360
        
    if (JMAG > 0) then
        MLAT  = ALATI
        MLONG = ALONG
    else
        LATI  = ALATI
        LONGI = ALONG
    end if

    CALL GEODIP(IYEAR,LATI,LONGI,MLAT,MLONG,JMAG)
    CALL FELDCOF(RYEAR)

    if (jf(18)) then
        call igrf_dip(lati,longi,ryear,300.0,dec,dip,magbr,modip)
    else
        CALL FIELDG(LATI,LONGI,300.0,XMA,YMA,ZMA,BET,DIP,DEC,MODIP)
        MAGBR = ATAN( 0.5 * TAN(DIP*UMR) ) / UMR
    endif

    ! calculate L-value, dip lati, and B_abs needed for invdip computation
    ! calculating invdip at 600 km
    invdip = -100.0
    if( (jf(2) .and. .not. jf(23)) .or. (jf(3) .and. .not. jf(6)) ) then
        call igrf_sub(lati,longi,ryear,600.0,fl,icode,dipl,babs)
        if(fl > 10.) fl = 10.
        invdip = INVDPC(FL,DIMO,BABS,DIPL)
    endif

    ABSLAT=ABS(LATI)
    ABSMLT=ABS(MLAT)
    ABSMDP=ABS(MODIP)
    ABSMBR=ABS(MAGBR)

    ! CALCULATION OF UT/LT and XMLT  ...............
    if(DHOUR <= 24.0) then
        HOUR = DHOUR                       ! dhour =< 24 is LT
        hourut = hour - longi/15.
        if(hourut < 0) hourut = hourut+24.
    else
        hourut = DHOUR - 25.               ! dhour>24 is UT+25
        hour = hourut + longi/15.          ! hour becomes LT
        if(hour > 24.) hour = hour - 24.
    end if
    CALL CLCMLT(IYEAR,DAYNR,HOURUT,LATI,LONGI,XMLT)

    ! SEASON assumes equal length seasons (92 days) with spring 
    ! (SEASON=1) starting at day-of-year=45; for lati < 0 adjustment 
    ! for southern hemisphere is made. Some models require the
    ! seasonal month (ISEAMON) or the seasonal day-of year (SEADAY)
    ! ZMONTH is decimal month (Jan 1 = 1.0 and Dec 31 = 12.97)
    ! SDAY is the day number reduced to a 360 day year (TOPH05)
    ! NRDAYM is the number of days in the current month 
    ! IDAYY is the number of days in the current year
    SEASON = INT( (DAYNR + 45.0) / 92.0 )
    if(SEASON < 1) SEASON = 4
    NSEASN = SEASON                ! Northern hemisphere season
    zmonth = month + (iday-1)*1./nrdaym
    ! NEW-GUL------------------------------
    sday = daynr/idayy*360.             
    ! NEW-GUL------------------------------
    seaday = daynr
    iseamon = month
    if (LATI >= 0.0) then
        SEASON = SEASON - 2            
        if(SEASON < 1) SEASON = SEASON + 4
        iseamon = month + 6
        if(iseamon > 12) iseamon = iseamon-12
        seaday = daynr + idayy/2.
        if(seaday > idayy) seaday = seaday-idayy
        ! NEW-GUL------------------------------
        sday=sday+180.                        
        if (sday > 360.) sday = sday - 360.    
        ! NEW-GUL------------------------------
    end if

    ! 12-month running mean sunspot number (rssn) and Ionospheric Global 
    ! index (gind), daily F10.7 cm solar radio flux (f107d) and monthly 
    ! F10.7 (cov) index   
    sam_mon  = (month == montho)
    sam_yea  = (iyear == iyearo)
    sam_doy  = (daynr == idaynro)
    sam_date = (sam_yea .and. sam_doy)
    sam_ut   = (hourut == ut0)
        
    if (sam_date .and..not. rzino .and..not. rzin &
                 .and..not. igin  .and..not. igino) then
     
        call tcon(iyear,month,iday,daynr,rzar,arig,ttt,nmonth)
        if(nmonth < 0) then ! exit routine
            icalls = icalls + 1
            return
        end if

        if (RZIN) then
            rrr = arzin
            rzar(1) = rrr
            rzar(2) = rrr
            rzar(3) = rrr
            !zi=-12.349154+(1.4683266-2.67690893e-03*rrr)*rrr
            !if(zi.gt.174.0) zi=174.0
            !arig(1) = zi
            !arig(2) = zi
            !arig(3) = zi
        end if
        if (IGIN) then
            zi = aigin
            arig(1) = zi
            arig(2) = zi
            arig(3) = zi
        end if
        rssn = rzar(3)
        gind = arig(3)
        cov  = 63.75 + rssn * (0.728 + rssn*0.00089)
        !rlimit=gind
        !COVSAT=63.75+rlimit*(0.728+rlimit*0.00089)        
        f107d   = cov
        f107pd  = cov
        f10781  = cov
        f107365 = cov
        if (jf(32) .or. jf(25)) then
            call APF_ONLY(iyear,month,iday,F107_daily,F107PD,F107_81,&
                          F107_365,IAP_daily,isdate)
            if (F107_daily > -11.1) then
                f107d   = f107_daily
                f107y   = f107PD
                f10781  = f107_81
                f107365 = f107_365
            end if
        else if (.not. jf(25)) then
            f107d  = f107din         ! F10.7D user input (adjusted expected)
            f107y  = f107din         ! F10.7PD user input
        else
            f10781 = f10781in        ! F10.7_81 user input
        end if
        
        pf107 = (f107d + f10781)/2.

        ! correcting F10.7 adjusted flux from APF107.DAT to flux observed at
        ! Earth that is expected by NRLMSIS00 (GTD7) and CHEMION          

        f_adj   = radj**2
        f107yo  = f107y/f_adj
        f10781o = f10781/f_adj
               
        if(jf(41)) cov = f107365                
        COVSAT  = cov
        if(covsat > 188.) covsat=188
    end if

    ! CALCULATION OF SOLAR ZENITH ANGLE (XHI/DEG), SUN DECLINATION ANGLE 
    ! (SUNDEC),SOLAR ZENITH ANGLE AT NOON (XHINON) AND TIME OF LOCAL 
    ! SUNRISE/SUNSET (SAX, SUX; dec. hours) AT 70 KM (D-REGION), 110 KM
    ! (E-REGION), 200 KM (F1-REGION), AND 500 KM (TE, TI).
    CALL SOCO(daynr,HOUR,LATI,LONGI,80. ,SUNDEC,XHI1  ,SAX80 ,SUX80)
    CALL SOCO(daynr,HOUR,LATI,LONGI,110.,SUD1  ,XHI2  ,SAX110,SUX110)
    CALL SOCO(daynr,HOUR,LATI,LONGI,200.,SUD1  ,XHI   ,SAX200,SUX200)
    CALL SOCO(daynr,HOUR,LATI,LONGI,300.,SUD1  ,XHI3  ,SAX300,SUX300)
    CALL SOCO(daynr,12.0,LATI,LONGI,110.,SUNDE1,XHINON,SAX1  ,SUX1)

    ! Night for D-layer
    dnight = .false.
    if(sax80 < -25.0) dnight = .true.
    if(sax80 <= sux80) then
        if ( hour < sax80 .or.  sux80 < hour ) dnight = .true.
    else
        if ( sux80 < hour .and. hour < sax80 ) dnight = .true.
    end if

    ! Night for E-layer
    enight = .false.
    if(sax110 < -25) enight = .true.
    if(sax110 <= sux110) then
        if ( hour < sax110 .or.  sux110 < hour ) enight = .true.
    else
        if ( sux110 < hour .and. hour < sax110 ) enight = .true.
    end if

    ! Night for F-layer
    fnight = .false.
    if(sax200 < -25) fnight = .true.
    if(sax200 <= sux200) then
        if ( hour < sax200 .or.  sux200 < hour ) fnight = .true.
    else
        if ( sux200 < hour .and. hour < sax200 ) fnight = .true.
    end if

    ! CALCULATION OF ELECTRON DENSITY PARAMETERS................
    ! lower height boundary (HNEA), upper boundary (HNEE)
    HNEA = 65.
    if(dnight) HNEA = 80.
    HNEE = 2000.
    if(.not. NODEN) then
 
        DELA = 4.32
        if(ABSMDP >= 18.) DELA = 1.0 + EXP(-(ABSMDP-30.0)/10.0)
        DELL = 1 + EXP(-(ABSLAT-20.)/10.)
 
        ! E peak critical frequency (foE), density (NmE), and height (hmE)
        if(FOEIN) then
            FOE = AFOE
            NME = ANME
        else
            FOE = FOEEDI(COV,XHI,XHINON,ABSLAT)
            NME = 1.24E10 * FOE**2
        end if
        if(HMEIN) then
            HME = AHME
        else
            HME = 110.0
        end if
 
        ! F2 peak critical frequency foF2, density NmF2, and height hmF2
        f2_in = FOF2IN .and. HMF2IN .and. itopn /= 2 
        use_ursi  = URSIF2 .eqv. URSIFO
        spwx_in = rzin .or. rzino .or. igin .or. igino 
        rightTime = nmonth == nmono .and. sam_yea
        readSpwx = .not. use_ursi .or. (use_ursi.and.spwx_in) .or. (use_ursi.and..not.spwx_in.and..not.sam_mon)
        readOldSpwx = use_ursi .and. .not.spwx_in .and. sam_mon .and..not.rightTime

        if (.not. f2_in) then
            if (readSpwx) then
                ! READ CCIR AND URSI COEFFICIENT SET FOR CHOSEN MONTH 
                URSIFO=URSIF2
                WRITE(FILNAM,104) MONTH+10
                OPEN(IUCCIR,FILE=FILNAM,STATUS='OLD',ERR=8448,FORM='FORMATTED')
                READ(IUCCIR,4689) F2,FM3
                CLOSE(IUCCIR)
         
                ! then URSI if chosen ....................................
                if(URSIF2) then
                    WRITE(FILNAM,1144) MONTH+10
                    OPEN(IUCCIR,FILE=FILNAM,STATUS='OLD',ERR=8448,FORM='FORMATTED')
                    READ(IUCCIR,4689) F2
                    CLOSE(IUCCIR)
                endif
            end if 
            if (readSpwx .or. readOldSpwx) then
     
                ! READ CCIR AND URSI COEFFICIENT SET FOR NMONTH, i.e. previous 
                ! month if day is less than 15 and following month otherwise 
         
                ! first CCIR ..............................................
                WRITE(FILNAM,104) NMONTH+10
                OPEN(IUCCIR,FILE=FILNAM,STATUS='OLD',ERR=8448,FORM='FORMATTED')
                READ(IUCCIR,4689) F2N,FM3N
                CLOSE(IUCCIR)
         
                ! then URSI if chosen .....................................
                if(URSIF2) then
                    WRITE(FILNAM,1144) NMONTH+10
                    OPEN(IUCCIR,FILE=FILNAM,STATUS='OLD',ERR=8448,FORM='FORMATTED')
                    READ(IUCCIR,4689) F2N
                    CLOSE(IUCCIR)
                end if
         
                ! LINEAR INTERPOLATION IN SOLAR ACTIVITY. IG12 used for foF2
                RR2  = ARIG(1)/100.
                RR2N = ARIG(2)/100.
                RR1  = 1. - RR2
                RR1N = 1. - RR2N
                do I=1,76
                    do J=1,13
                        K=J+13*(I-1)
                        FF0N(K)=F2N(J,I,1)*RR1N+F2N(J,I,2)*RR2N
                        FF0(K)=F2(J,I,1)*RR1+F2(J,I,2)*RR2
                    end do
                end do
         
                RR2  = RZAR(1)/100.
                RR2N = RZAR(2)/100.
                RR1  = 1. - RR2
                RR1N = 1. - RR2N
                DO I=1,49
                    DO J=1,9
                        K      = J + 9 * (I - 1)
                        XM0N(K)= FM3N(J,I,1)* RR1N+ FM3N(J,I,2)*RR2N
                        XM0(K) = FM3(J,I,1) * RR1 + FM3(J,I,2) *RR2
                    end do
                end do

            end if 

            zfof2  =  FOUT(MODIP,LATI,LONGI,HOURUT,FF0)
            fof2n  =  FOUT(MODIP,LATI,LONGI,HOURUT,FF0N)
            zm3000 = XMOUT(MODIP,LATI,LONGI,HOURUT,XM0)
            xm300n = XMOUT(MODIP,LATI,LONGI,HOURUT,XM0N)
            midm   = 15
            if(month == 2) midm=14
            if (iday < midm) then
                yfof2  = fof2n + ttt * (zfof2-fof2n)
                xm3000 = xm300n+ ttt * (zm3000-xm300n)
            else
                yfof2  = zfof2 + ttt * (fof2n-zfof2)
                xm3000 = zm3000+ ttt * (xm300n-zm3000)
            end if
        end if 

        if(FOF2IN) then
            FOF2 = AFOF2
            NMF2 = ANMF2
        else
            FOF2 = YFOF2
            NMF2 = 1.24E10 * FOF2**2
        end if
 
        ! stormtime updating for foF2 (foF2s, NmF2s) and foE (foEs,
        ! NmEs) and auroral boundary computation.
        foF2s     = foF2
        foEs      = foE
        NmF2s     = NmF2
        NmEs      = NmE
        stormcorr = -1.
        estormcor = -1.
        fstorm_on = jf(26) .and. jf(8)
        estorm_on = jf(35) .and. jf(15)
        if(fstorm_on .or. jf(33) .or. estorm_on) &
           call apf(isdate,hourut,indap)
           
        ! stormtime updating for foF2 (foF2s, NmF2s) 
        if ( fstorm_on .and. (indap(1) > -1) ) then    
            icoord = 1
            kut    = int(hourut)
            call STORM(indap,lati,longi,icoord,cglat,kut,daynr,stormcorr)
            fof2s  = fof2 * stormcorr
            NMF2S  = 1.24E10 * FOF2S**2
        end if
 
        ! stormtime updating for foE (foEs, NmEs)
        if ( estorm_on .and. (indap(1) > -1) ) then    
            estormcor = STORME_AP(DAYNR,MLAT,INDAP(13)*1.0)
            if(estormcor > -2.0) foes = foe*estormcor
            NMES = 1.24E10 * FOES**2
        end if
 
        ! calculation of equatorward auroral boundary
        if (jf(33)) then
            if( indap(1) > -1 ) then
                xkp = ckp( indap(13) )
            else     
                xkp = 3.0
            end if   
            ! Corrected magnetic latitude CGM of equatorward boundary, 
            ! ab_mlat(48), for MLT=0.0,0.5,1.0 ... 23.5 h and kp=xkp
            call auroral_boundary(xkp,-1.0,cgmlat,ab_mlat)
            DAT(1,1) = lati
            DAT(2,1) = longi
            call GEOCGM01(1,IYEAR,height_center,DAT,PLA,PLO)
            cgm_lat  = DAT(3,3)
            cgm_lon  = DAT(4,3)
            cgm_mlt00_ut = DAT(11,3)
            cgm_mlt  = hourut - cgm_mlt00_ut
            if(cgm_mlt < 0.) cgm_mlt = 24. + hourut - cgm_mlt00_ut
 
            ! CGM latitude of boundary (cgmlat) for present MLT value
            ! 2012.02 12/17/12 Add magnetic declination as oarr(84) output
            zmlt=xmlt
            cgmlat=100.0
            if( zmlt >= 0.0 .and. zmlt <= 24.0) &
                call auroral_boundary(xkp,zmlt,cgmlat,ab_mlat)
        end if
 
        if( ( .not.JF(4) .and. JF(31) ) .or. ( .not.JF(39) .and. JF(40) ) ) then
            FLON = LONGI + 15. * hourut
            if(FLON > 360.) FLON = FLON-360.
            X11  = -90
            X22  = 90
            FX11 = fmodip(x11)
            FX22 = fmodip(x22)
            CALL REGFA1(X11,X22,FX11,FX22,0.001,MODIP,FMODIP,SCHALT,XRLAT)
            if(SCHALT) then
                XRLAT = LATI
                if(mess) WRITE(KONSOL,656) LATI
            end if
            RLAT = XRLAT
        end if
        
        if(HMF2IN) then
            if(AHMF2 < 50.0) then
                  XM3000 = AHMF2
                  HMF2   = HMF2ED(MAGBR,RSSN,FOF2/FOE,XM3000)
              else 
                  HMF2=AHMF2
              end if
        else if(JF(39)) then
            ratf = fof2/foe
            ! set jf(36)=false to use foF2_storm in hmF2 formula  
            if( .not.jf(36) ) ratf = fof2s/foe
            HMF2 = HMF2ED(MAGBR,RSSN,RATF,XM3000)
        else if(JF(40)) then
            ! AMTB digisonde model        
            CALL SHAMDHMF2(RLAT,FLON,ZMONTH,RSSN,HMF2)
        else
            ! SHUBIN-COSMIC model
            CALL SDMF2(hourut,month,F10781,modip,longi,HMF2)
        end if
 
        nmono   = nmonth
        MONTHO  = MONTH
        iyearo  = iyear
        idaynro = daynr
        rzino   = rzin
        igino   = igin
        ut0     = hourut
 
        ! topside profile parameters .............................
        COS2 = COS(MLAT*UMR)
        COS2 = COS2**2
        FLU  = (COVSAT - 40.0)/30.0
 
        ! option to use unlimited F10.7M for the topside
        ! previously: if(OLD79) ETA1=-0.0070305*COS2
        if(OLD79) FLU = (COV - 40.0)/30.0
        FO1  = FOF2S
        if(JF(37)) FO1 = FOF2
        EX   = EXP(-MLAT/15.)
        EX1  = EX + 1
        EPIN = 4. * EX/(EX1**2)
        ETA1 = -0.02 * EPIN
        ETA  = 0.058798 + ETA1 - &
               FLU * (0.014065  - 0.0069724 * COS2) + &
               FO1* (0.0024287 + 0.0042810 * COS2  - 0.0001528 * FO1)
        ZETA = 0.078922 - 0.0046702 * COS2 - &
               FLU * (0.019132  - 0.0076545 * COS2) + &
               FO1* (0.0032513 + 0.0060290 * COS2  - 0.00020872 * FO1)
        BETA = -128.03 + 20.253 * COS2 - &
               FLU * (8.0755  + 0.65896 * COS2) + &
               FO1* (0.44041 + 0.71458 * COS2 - 0.042966 * FO1)
        Z    = EXP(94.5/BETA)
        Z1   = Z + 1
        Z2   = Z/(BETA * Z1**2)
        DELTA= ( ETA/Z1 - ZETA/2.0 )/( ETA*Z2 + ZETA/400.0 )
 
        ! Correction term for topside (Bilitza) depends on modip, hour,
        ! sax300, sux300, and hmF2
        if (itopn == 1) then
            zmp1   = exp(modip / 10.)
            zmp11  = 1. + zmp1
            zmp111 = zmp1 / (zmp11 * zmp11)
            zmp2   = exp(modip / 19.)
            zmp22  = 1. + zmp2
            zmp222 = zmp2 / (zmp22 * zmp22) 
            r2n    = -0.84 - 1.6 * zmp111  
            r2d    = -0.84 - 0.64 * zmp111 
            x1n    = 230. - 700. * zmp222  
            x1d    = 550. - 1900. * zmp222
            r2     = HPOL(HOUR,r2d,r2n,SAX300,SUX300,1.,1.)
            x1     = HPOL(HOUR,x1d,x1n,SAX300,SUX300,1.,1.)
            hcor1  = hmF2 + x1
            x12    = 1500. - x1
            tc3    = r2 / x12
        end if
 
        ! NEW-GUL--------------------------------
        !
        ! Correction term for topside (Gulyaeva)
        if (itopn == 3) then
            hei05  = 0.
            CALL TOPH05(COV,MLAT,HOUR,HMF2,HEI05,SDAY)
            h05top = hei05
            xnetop = XE_1(H05TOP)
        end if
        ! NEW-GUL--------------------------------
 
        ! NeQuick topside parameters (use CCIR-M3000F2 even if user-hmF2)
        if (itopn == 2) then
            fo2    = foF2s
            if(jf(37)) fo2 = foF2
            dNdHmx = -3.467 + 1.714*log(fo2) + 2.02*log(XM3000)
            dNdHmx = exp(dNdHmx) * 0.01
            B2bot  = 0.04774 * fo2**2 / dNdHmx
            b2k    = 3.22 - 0.0538*fo2 - 0.00664*hmF2 + 0.113*hmF2/B2bot + 0.00257*rssn
            ee     = exp(2.0*(b2k-1.0))
            b2k    = (b2k*ee + 1.0)/(ee + 1.0)
            B2TOP  = b2k * B2bot
        endif
 
        ! Bottomside thickness parameter B0 and shape parameters B1
        if(jf(4)) then
            B0 = B0_98(HOUR,SAX200,SUX200,NSEASN,RSSN,LONGI,MODIP)
            B1 = HPOL(HOUR,1.9,2.6,SAX200,SUX200,1.,1.)
        else if (jf(31)) then
            CALL SHAMDB0D (RLAT,FLON,ZMONTH,RSSN,B0)
            CALL SHAB1D (LATI,FLON,ZMONTH,RSSN,B1)
        else
            CALL ROGUL(SEADAY,XHI,SEAX,GRAT)
            IF (FNIGHT) GRAT = 0.91d0 - HMF2/4000.d0
            B1     = HPOL(HOUR,1.9,2.6,SAX200,SUX200,1.,1.)
            BCOEF  = B1*(B1*(0.0046d0*B1 - 0.0548d0) + 0.2546d0) + 0.3606d0
            B0CNEW = HMF2*(1.d0 - GRAT)
            B0     = B0CNEW/BCOEF
        endif
 
        ! F1 layer height hmF1, critical frequency foF1, peak density NmF1
        if(FOF1IN) then
            FOF1 = AFOF1
            NMF1 = ANMF1
        ELSE
            FOF1 = FOF1ED(ABSMBR,RSSN,XHI)
            NMF1 = 1.24E10 * FOF1**2
        ENDIF
 
        ! F1 layer thickness parameter c1
        c1 = f1_c1(modip,hour,sux200,sax200)
 
        ! F1 occurrence probability with Scotto et al. 1997 or Ducharme et al. 
        ! if jf(19)=f1_ocpro=.true. or .false.
        ! If .not.jf(20)=f1_l_cond=.true. then Scotto model with L-condition
        if (f1_ocpro) then
            call f1_prob(xhi,mlat,rssn,f1pbw,f1pbl)
            f1pb  = f1pbw
            if(f1_l_cond) f1pb = f1pbl
            f1reg = .false.
            if( fof1in .or. (f1pb >= 0.5) ) f1reg = .true.
        else            
            f1pb  = 0.0 
            if( fof1in .or. (.not.fnight .and. (fof1 > 0.0)) ) f1pb = 1. 
            f1reg = .false.
            if(f1pb > 0.0) f1reg = .true.
        end if
 
        ! E-valley: DEPTH=(NmE-N_deepest)/NmE*100, WIDTH=HEF-HmE, 
        ! distance of deepest value point above E-peak(HDEEP), 
        ! derivative at valley top divided by NmE (DLNDH), 
        ! and height of valley top (HEF)
        XDEL   = XDELS(SEASON)/DELA
        DNDHBR = DNDS(SEASON)/DELA
        HDEEP  = HPOL(HOUR, 10.5/DELA, 28.         , SAX110,SUX110,1.,1.)
        WIDTH  = HPOL(HOUR, 17.8/DELA, 45.+22./DELA, SAX110,SUX110,1.,1.)
        DEPTH  = HPOL(HOUR, XDEL     , 81.         , SAX110,SUX110,1.,1.)
        DLNDH  = HPOL(HOUR, DNDHBR   , .06         , SAX110,SUX110,1.,1.)
        if(depth >= 1.0) then
            if(enight) depth = -depth
            call TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
            if(EXT .and. mess) WRITE(KONSOL,650)
            if(EXT) WIDTH = .0
        else
            WIDTH = .0
        end if
        HEF    = HME + WIDTH
        hefold = hef
        VNER   = (1. - ABS(DEPTH) / 100.) * NMES
 
        ! Parameters below E  .............................
        hmex = hme - 9.
        NMD  = XMDED(XHI,RSSN,4.0E8)
        HMD  = HPOL(HOUR, 81.0          , 88.0,SAX80,SUX80,1.,1.)
        F(1) = HPOL(HOUR, 0.02+0.03/DELA, 0.05,SAX80,SUX80,1.,1.)
        F(2) = HPOL(HOUR, 4.6           ,  4.5,SAX80,SUX80,1.,1.)
        F(3) = HPOL(HOUR, -11.5         , -4.0,SAX80,SUX80,1.,1.)
        FP1  = F(1)
        FP2  = -FP1**2/2.0
        FP30 = ( -F(2)*FP2 - FP1 + 1.0/F(2) )/(F(2)**2)
        FP3U = ( -F(3)*FP2 - FP1 - 1.0/F(3) )/(F(3)**2)
        HDX  = HMD + F(2)
 
        ! indermediate region between D and E region; parameters xkk
        ! and d1 are found such that the function reaches hdx/xdx/dxdh
        X    = HDX - HMD
        XDX  = NMD*EXP( X*( FP1 + X*( FP2 + X*FP30 ) ) )
        DXDX = XDX*( FP1 + X*(2.0*FP2 + X*3.0*FP30) )
        X    = HME - HDX
        XKK  = -DXDX*X/( XDX*ALOG(XDX/NMES) )
 
        ! if exponent xkk is larger than xkkmax, then xkk will be set to 
        ! xkkmax and d1 will be determined such that the point hdx/xdx is 
        ! reached; derivative is no longer continuous.
        xkkmax = 5.
        if (xkk > xkkmax) then
                xkk = xkkmax
                d1  = -alog(xdx/nmes) / (x**xkk)
        else
                D1  = DXDX/( XDX * XKK * X**(XKK - 1.0) )
        endif
 
        ! compute Danilov et al. (1995) D-region model values
        if(.not. dreg) then
            vKp  = 1.
            f5sw = 0.
            f6wa = 0.
            call DRegion(xhi,month,f107d,vKp,f5SW,f6WA,elg)
            do ii=1,11
              ddens(1,ii) = -1.
              if(ii < 8) ddens(1,ii) = 10**( elg(ii) + 6 )
            end do
            f5sw = 0.5
            f6wa = 0.
            call DRegion(xhi,month,f107d,vKp,f5SW,f6WA,elg)
            do ii=1,11
              ddens(2,ii) = -1.
              if(ii < 8) ddens(2,ii) = 10**( elg(ii) + 6 )
            end do
            f5sw = 1.
            f6wa = 0.
            call DRegion(xhi,month,f107d,vKp,f5SW,f6WA,elg)
            do ii=1,11
              ddens(3,ii) = -1.
              if(ii < 8) ddens(3,ii) = 10**( elg(ii) + 6 )
            end do
            f5sw = 0.
            f6wa = 0.5
            call DRegion(xhi,month,f107d,vKp,f5SW,f6WA,elg)
            do ii=1,11
              ddens(4,ii) = -1.
              if(ii < 8) ddens(4,ii) = 10**( elg(ii) + 6 )
            end do          
            f5sw = 0.
            f6wa = 1.
            call DRegion(xhi,month,f107d,vKp,f5SW,f6WA,elg)
            do ii=1,11
              ddens(5,ii) = -1.
              if(ii < 8) ddens(5,ii) = 10**( elg(ii) + 6 )
            end do
        end if
 
        ! SEARCH FOR HMF1 ..................................................
        if(.not. LAYVER) then
            hmf1 = 0
            if(F1REG) then
                bnmf1 = 0.9 * nmf1
                f1GTe = nmes < bnmf1 
                firstIter = .true.
                do
                    nSchalt = .true. ! this WILL be reset if it's needed
                    if (f1GTe .or. .not. firstIter) then
                        do
                            xe2h = xe2(hef)
                            xeLEf1  = xe2h <= bnmf1
                            if (xeLEf1) exit
                            hef  = hef - 1 
                            if (hef <= hme) exit
                        end do
                        if (xeLEf1) then
                            CALL REGFA1(HEF,HMF2,XE2H,NMF2S,0.001,NMF1,XE2,SCHALT,HMF1) 
                            nSchalt = .not. SCHALT
                        else
                            hef    = hme  
                            width  = 0.0  
                            hefold = hef   
                        end if
                    end if
                    if ( (.not.f1GTe .and. firstIter) .or. &
                        .not. (xeLEf1 .and. nSchalt) ) then
                        if(mess) WRITE(KONSOL,11) 
                        HMF1=0.
                        F1REG=.FALSE.
                    end if
                    if (hef == hefold) exit
                    width  = hef - hme
                    if(ENIGHT) DEPTH = -DEPTH
                    CALL TAL(HDEEP,DEPTH,WIDTH,DLNDH,EXT,E)
                    if (.not. EXT) exit
                    if(mess) WRITE(KONSOL,650)
                    WIDTH  = .0
                    hef    = hme
                    hefold = hef
                    first_iter = .false.
                end do
            end if
     
            ! SEARCH FOR HST [NE3(HST)=NMEs] ......................................
            if(F1REG) then
                hf1 = hmf1
                xf1 = nmf1
            else
                hf1 = (hmf2 + hef)/2.
                xf1 = xe2(hf1)
            endif
     
            hf2 = hef
            xf2 = xe3_1(hf2)
            if(xf2 <= nmes) &
                CALL REGFA1(hf1,HF2,XF1,XF2,0.001,NMES,XE3_1,SCHALT,HST)

            if (xf2 <= nmes .and. .not. schalt) then
                HZ = (HST + HF1)/2.0
                D  = HZ - HST
                T  = D**2/( HZ - HEF - D)
            else
                ! assume linear interpolation between HZ and HEF ..................
                if(mess) WRITE(KONSOL,100)
                HZ    = (HEF + HF1)/2.
                xnehz = xe3_1(hz)
                if(mess) WRITE(KONSOL,901) HZ,HEF
                T     = (XNEHZ - NMES)/(HZ - HEF)
                HST   = -333.
            end if
        else
            ! LAY-functions for middle ionosphere
            if(HMF1IN) then
                HMF1M = AHMF1
            else
                HMF1M = 165. + 0.6428 * XHI
            end if
            HHALF = GRAT * HMF2
            HV1R  = HME + WIDTH
            HV2R  = HME + HDEEP
            HHMF2 = HMF2
            CALL INILAY(FNIGHT,F1REG,NMF2S,NMF1,NMES,VNER,HHMF2,HMF1M,HME,&
                        HV1R,HV2R,HHALF,HXL,SCL,AMP,IIQU)
            if((IIQU == 1).and.mess) WRITE(KONSOL,7733)
            if((IIQU == 2).and.mess) WRITE(KONSOL,7722)
        end if
    end if

    !---------- CALCULATION OF NEUTRAL TEMPERATURE PARAMETER-------
    HTA   = 60.0
    HEQUI = 120.0
    if(.not. NOTEM) then
        SEC   = hourut * 3600.
        CALL APFMSIS(ISDATE,HOURUT,IAPO)
        if(iapo(2) < 0.0) then
            SWMI(9) = 0.
            IAPO(1) = 0.
        else
            SWMI(9) = -1.0
        end if           
        CALL TSELEC(SWMI)
        CALL GTD7(IYD,SEC,HEQUI,LATI,LONGI,HOUR,F10781o,F107Yo,IAPO,&
                  0,D_MSIS,T_MSIS)
        TN120 = T_MSIS(2)
        if(HOUR /= 0.0) then
            secni= (24. - longi/15)*3600.
            CALL GTD7(IYD,SECNI,HEQUI,LATI,LONGI,0.0,F10781o,F107Yo,&
                      IAPO,0,D_MSIS,T_MSIS)
        end if
        TN1NI= T_MSIS(2)
    
        !--------- CALCULATION OF ELECTRON TEMPERATURE PARAMETER--------
        ! Te(120km) = Tn(120km)
        AHH(1) = 120.
        ATE(1) = TN120
    
        ! Te-MAXIMUM based on JICAMARCA and ARECIBO data 
    
        HMAXD  = 60. * EXP( -(MLAT/22.41)**2 ) + 210.
        HMAXN  = 150.
        AHH(2) = HPOL(HOUR,HMAXD,HMAXN,SAX200,SUX200,1.,1.)
        TMAXD  = 800.* EXP( -(MLAT/33.  )**2 ) + 1500.
        CALL GTD7(IYD,SECNI,HMAXN,LATI,LONGI,0.0,F10781o,F107Yo,&
                  IAPO,0,D_MSIS,T_MSIS)
        TMAXN  = T_MSIS(2)
        ATE(2) = HPOL(HOUR,TMAXD,TMAXN,SAX200,SUX200,1.,1.)
    
        ! Te(300km), Te(400km) from AE-C, Te(1400km), Te(3000km) from 
        ! ISIS, Brace and Theis
        DIPLAT=MAGBR
        CALL TEBA(DIPLAT,HOUR,NSEASN,TEA)
    
        icd=0              
        if(jf(23)) then
            ! Te at fixed heights taken from Brace and Theis
            AHH(3) = 300.
            AHH(4) = 400.
            AHH(5) = 600.
            AHH(6) = 1400.
            AHH(7) = 3000.
            hte    = 3000
            ATE(3) = TEA(1)
            ATE(4) = TEA(2)
            ATE(6) = TEA(3)
            ATE(7) = TEA(4)
    
            ! Te(600km) from AEROS, Spenner and Plugge (1979)
            ETT    = EXP(-MLAT/11.35)
            TET    = 2900. - 5600.*ETT/( (ETT+1)**2. )
            TEN    = 839. + 1161./( 1. + EXP( -(ABSMLT-45.)/5. ) )
            ATE(5) = HPOL(HOUR,TET,TEN,SAX300,SUX300,1.5,1.5)
        else
            ! New model with solar activity effects included (Truhlik et al., 2011)
            ! Te at fixed heights 350, 550, 850, 1400, and 2000 km
            AHH(3) = 350.
            AHH(4) = 550.
            AHH(5) = 850.
            AHH(6) = 1400.
            AHH(7) = 2000.
            hte    = 2500
            isa    = 0
            if(jf(42)) isa=1
            do ijk=3,7
                call elteik(isa,invdip,xmlt,ahh(ijk),daynr,&
                            pf107,teh2,sdte)
                ate(ijk)=teh2
            end do
        end if
    
        ! Option to use Te = f(Ne) relation at ahh(3), ahh(4)
        if(TENEOP) then
            DO I=1,2
                if(TECON(I)) ATE(I+2) = TEDE(AHH(I+2), XNAR(I), -COV)
            end do
        end if
    
        ! Te corrected and Te > Tn enforced
        CALL GTD7(IYD,SEC,AHH(2),LATI,LONGI,HOUR,F10781o,F107Yo,&
                  IAPO,0,D_MSIS,T_MSIS)
        TNAHH2 = T_MSIS(2)
        if(ATE(2) < TNAHH2) ATE(2) = TNAHH2
        STTE1  = ( ATE(2) - ATE(1) )/( AHH(2) - AHH(1) )
        do I=2,6
            CALL GTD7(IYD,SEC,AHH(I+1),LATI,LONGI,HOUR,F10781o,&
                      F107Yo,IAPO,0,D_MSIS,T_MSIS)
            TNAHHI = T_MSIS(2)
            if(ATE(I+1) < TNAHHI) ATE(I+1) = TNAHHI
            STTE2  = ( ATE(I+1) - ATE(I) )/( AHH(I+1) - AHH(I) )
            ATE(I) = ATE(I) - (STTE2 - STTE1)*DTE(I-1)*ALOG2
            STTE1  = STTE2
        end do
    
        ! Te gradients STTE are computed for each segment
        do I=1,6
            STTE(I) = ( ATE(I+1) - ATE(I) )/( AHH(I+1) - AHH(I) )
        end do
    
        ATE1 = ATE(1)
    
        !------------ CALCULATION OF ION TEMPERATURE PARAMETERS--------
        ! Ti(430km) during daytime from AEROS data
        XSM1   = 430.0
        XSM(1) = XSM1
        Z1     = EXP(-0.09*MLAT)
        Z2     = Z1 + 1.
        TID1   = 1240.0 - 1400.0 * Z1 / ( Z2 * Z2 )
        MM(2)  = HPOL(HOUR,3.0,0.0,SAX300,SUX300,1.,1.)
    
        ! Ti(430km) duirng nighttime from AEROS data
        Z1   = ABSMLT
        Z2   = Z1*( 0.47 + Z1*0.024 )*UMR
        Z3   = COS(Z2)
        TIN1 = 1200.0 - 300.0 * SIGN(1.0, Z3) * SQRT(ABS(Z3))
    
        ! Ti(430km) for specified time using HPOL
        TI1  = TIN1  
        if(TID1 > TIN1) TI1 = HPOL(HOUR,TID1,TIN1,SAX300,SUX300,1.,1.)
    
        ! Tn < Ti < Te enforced
        TEN1 = ELTE(XSM1)
        CALL GTD7(IYD,SECNI,XSM1,LATI,LONGI,0.0,F10781o,F107Yo,&
                  IAPO,0,D_MSIS,T_MSIS)
        TNN1 = T_MSIS(2)
        if(TEN1 < TNN1) TEN1 = TNN1
        if(TI1  > TEN1) TI1  = TEN1
        if(TI1  < TNN1) TI1  = TNN1
    
        ! Tangent on Tn profile determines HS
        HS   = 200.
        CALL GTD7(IYD,SEC,HS,LATI,LONGI,HOUR,F10781o,F107Yo,&
                  IAPO,0,D_MSIS,T_MSIS)
        TNHS = T_MSIS(2)
        MM(1)= (TI1 - TNHS)/(XSM1 - HS)
        MXSM = 2
    
        ! XTETI is altitude where Te=Ti
        XTTS  = 500.
        X     = 500.
        do
            X   = X + XTTS
            ! if X gets too large, exit loop
            if(X >=  AHH(7)) exit
            TEX = ELTE(X)
            TIX = TI(X)
            if(TIX >= TEX) then
                X    = X - XTTS
                XTTS = XTTS/10.
                ! if increment is too small, exit loop
                if(XTTS <= 0.1) exit
            end if
        end do
    
        if (X < AHH(7)) then
            XTETI = X + XTTS*5.
            ! Ti=Te above XTETI 
            MXSM   = 3
            MM(3)  = STTE(6)
            XSM(2) = XTETI
            if(XTETI <= AHH(6)) then
                MXSM   = 4
                MM(3)  = STTE(5)
                MM(4)  = STTE(6)
                XSM(3) = AHH(6)
                if(XTETI <= AHH(5)) then
                    MXSM   = 5
                    DTI(1) = 5.
                    DTI(2) = 5.
                    MM(3)  = STTE(4)
                    MM(4)  = STTE(5)
                    MM(5)  = STTE(6)
                    XSM(3) = AHH(5)
                    XSM(4) = AHH(6)
                end if
            end if
        end if
    end if 

    ! CALCULATION OF ION DENSITY PARAMETER..................
    if(.not. NOION) then
        HNIA = 75.
        if(DY) HNIA = 80.
        HNIE = 2000.
    end if

    ! CALCULATION FOR THE REQUIRED HEIGHT RANGE.......................
    ! In the absence of an F1 layer hmf1=hz since hmf1 is used in XE
    xhmf1  = hmf1
    if(hmf1 <= 0.0) HMF1 = HZ
    height = heibeg
    kk     = 1

    ! PARAMETER COMPUTATION LOOP
    do
        if (.not. NODEN .or. &
            ( HNEA <= Height .and. Height <= HNEE )) then
                if(LAYVER) then
                    ELEDE = -9.
                    if(IIQU < 2) &
                        ELEDE = XEN(HEIGHT,HMF2,NMF2S,HME,4,HXL,SCL,AMP)
                    outf(1,kk) = elede
                else
                    ELEDE = XE_1(HEIGHT)
                    ! electron density in m-3 in outf(1,*)
                    OUTF(1,kk)=ELEDE
                end if
            end if
        end if

        ! plasma temperatures
        if (.not.NoTem .or. &
            (HTA <= Height .and. Height <= HTE)) then
            CALL GTD7(IYD,SEC,HEIGHT,LATI,LONGI,HOUR,&
                      F10781o,F107Yo,IAPO,0,D_MSIS,T_MSIS)
            TNH = T_MSIS(2)
            TIH = TNH
            if(HEIGHT > HS) then
                TIH = TI(HEIGHT)
                if(TIH < TNH) TIH = TNH
            end if
            TEH = TNH
            if(HEIGHT > HEQUI) then 
                TEH = ELTE(HEIGHT)
                if(TEH < TIH) TEH = TIH
            endif

            OUTF(2,kk) = TNH
            OUTF(3,kk) = TIH
            OUTF(4,kk) = TEH
        end if

        ! ion composition
        if (.not.NOION &
            (HNIA <= HEIGHT .and. HEIGHT <= HNIE )) then

            ROX    =-1.
            RHX    =-1.
            RNX    =-1.
            RHEX   =-1.
            RNOX   =-1.
            RO2X   =-1.
            RCLUST =-1.

            if(DY) then
                if (height > 300.) then
                    call CALION(invdip,xmlt,height,daynr,&
                                f107d,xic_O,xic_H,xic_He,xic_N)
                    rox  = xic_O * 100.
                    rhx  = xic_H * 100.
                    rnx  = xic_N * 100.
                    rhex = xic_He* 100.
                    rnox = 0.
                    ro2x = 0.
                else
                    ! Richards-Bilitza-Voglozin-2010 IDC model
                    CALL GTD7(IYD,SEC,height,lati,longi,&
                              HOUR,f10781o,f107yo,IAPO,&
                              48,D_MSIS,T_MSIS)
                    XN4S   = 0.5 * D_MSIS(8)
                    EDENS  = ELEDE/1.e6
                    jprint = 1
                    if(jf(38)) jprint = 0
                    CALL CHEMION(jprint,height,F107Yo,F10781o,&
                                 TEH,TIH,TNH,D_MSIS(2),D_MSIS(4),&
                                 D_MSIS(3),D_MSIS(1),-1.0,XN4S,&
                                 EDENS,-1.0,xhi,ro,ro2,rno,rn2,&
                                 rn,Den_NO,Den_N2D,INEWT)
                    sumion = edens/100.
                    rox    = ro/sumion
                    rhx    = 0.
                    rhex   = 0.
                    rnx    = rn/sumion
                    rnox   = rno/sumion
                    ro2x   = ro2/sumion
                end if
            else
                ! Danilov-Smirnova-1995 model and Danilov-Yaichnikov-1985 model (upper)
                call iondani(iday,iseamon,height,xhi,&
                             lati,f107365,dion)
                ROX    = DION(1)
                RHX    = DION(2)
                RNX    = DION(3)
                RHEX   = DION(4)
                RNOX   = DION(5)
                RO2X   = DION(6)
                RCLUST = DION(7)
            end if

            ! ion densities are given in percent of total electron density;
            if(jf(22)) then 
                xnorm = 1
            else
                xnorm = elede/100.
            end if
            OUTF(5,kk)  = ROX*xnorm
            OUTF(6,kk)  = RHX*xnorm
            OUTF(7,kk)  = RHEX*xnorm
            OUTF(8,kk)  = RO2X*xnorm
            OUTF(9,kk)  = RNOX*xnorm
            OUTF(10,kk) = RCLUST*xnorm
            OUTF(11,kk) = RNX*xnorm
        end if

        ! D region special: Friedrich&Torkar model in outf(13,*)
        if(.not.dreg .and. height <= 140.) then
            outf(1,kk) = -1.
            call F00(HEIGHT,LATI,DAYNR,XHI,F107D,EDENS,IERROR)
            if(ierror == 0 .or. ierror == 2) outf(1,kk) = edens
        end if

        height = height + heistp
        kk     = kk + 1
        if(kk <= numhei) exit
    end do
    ! END OF PARAMETER COMPUTATION LOOP 

    ! D region special: densities for 11 heights (60,65,70,..,110km)
    if (.not. dreg) then
        do ii=1,11
            Htemp           = 55 + ii*5  
            outf(14, ii)    = -1.     
            if(Htemp >= 65.) outf(14,ii) = XE6(Htemp)     
            outf(14, 11+ii) = -1.
            call F00(Htemp,LATI,DAYNR,XHI,F107D,EDENS,IERROR)
            if(ierror == 0 .or. ierror == 2) outf(14,11+ii)=edens
            outf(14, 22+ii) = ddens(1,ii)      
            outf(14, 33+ii) = ddens(2,ii)      
            outf(14, 44+ii) = ddens(3,ii)      
            outf(14, 55+ii) = ddens(4,ii)      
            outf(14, 66+ii) = ddens(5,ii)      
        end do
    end if

    ! equatorial vertical ion drift
    drift=-1.
    if (jf(21) .and. abs(magbr) < 25.0) then
        param(1) = daynr
        param(2) = f107d
        call vdrift(hour,longi,param,drift)
    end if

    ! spread-F occurrence probability
    spreadf = -1.
    if ( jf(28) .and. &
         ( hour <= 7.25 .or. 17.75 <= hour ) .and. &
         ( abs(lati) <= 25.0 )) then
        spfhour = hour
        daynr1  = daynr
        if(hour < 12.0) then
            spfhour = hour  + 24.0
            daynr1  = daynr - 1
            if(daynr1 < 1) daynr1 = idayy
        end if
        call spreadf_brazil(daynr,idayy,f107d,lati,osfbr)
        ispf = int( (spfhour - 17.75)/0.5 ) + 1
        if(ispf > 0 .and. ispf < 26) spreadf = osfbr(ispf)
    end if

    ! ADDITIONAL PARAMETER FIELD OARR
    if(.not.NODEN) then
        OARR(1)  = NMF2S
        OARR(2)  = HMF2
        if(f1reg) OARR(3) = NMF1
        if(f1reg) OARR(4) = XHMF1
        OARR(5)  = NMES
        OARR(6)  = HME
        OARR(7)  = NMD
        OARR(8)  = HMD
        OARR(9)  = HHALF
        OARR(10) = B0
        OARR(11) = VNER
        OARR(12) = HEF
    end if
    if(.not. NOTEM) then
        OARR(13) = ATE(2)
        OARR(14) = AHH(2)
        OARR(15) = ATE(3)
        OARR(16) = ATE(4)
        OARR(17) = ATE(5)
        OARR(18) = ATE(6)
        OARR(19) = ATE(7)
        OARR(20) = ATE(1)
        OARR(21) = TI1
        OARR(22) = XTETI
    end if
    OARR(23) = XHI
    OARR(24) = SUNDEC
    OARR(25) = DIP
    OARR(26) = MAGBR
    OARR(27) = MODIP
    OARR(28) = LATI        
    OARR(29) = SAX200
    OARR(30) = SUX200
    OARR(31) = SEASON
    OARR(32) = LONGI        
    OARR(33) = rssn
    OARR(34) = COV
    OARR(35) = B1
    OARR(36) = xm3000
    ! OARR(37) used for TEC and 38 for TEC-top
    OARR(39) = gind
    OARR(40) = f1pb
    OARR(41) = f107d
    OARR(42) = c1
    OARR(43) = daynr
    OARR(44) = drift
    OARR(45) = stormcorr
    OARR(46) = f10781
    OARR(47) = estormcor
    OARR(48) = spreadf
    OARR(49) = MLAT
    OARR(50) = MLONG
    OARR(51) = indap(13)*1.0   ! ap for current UT
    OARR(52) = IAP_daily*1.0   ! daily ap
    OARR(53) = invdip
    OARR(54) = XMLT
    OARR(55) = cgm_lat
    OARR(56) = cgm_lon
    OARR(57) = cgm_mlt
    OARR(58) = cgmlat   ! CGM latitude of equatorward boundary
    ! include only every second auroral boundary point (MLT=0,1,2..23)
    jjj=58
    do iii=1,47,2
        jjj       = jjj + 1 
        oarr(jjj) = ab_mlat(iii)
    end do
    OARR(83) = xkp
    OARR(84) = dec
    OARR(85) = fl
    OARR(86) = dimo

    icalls=icalls+1

    return

    ! File IO Error
8448    write(konsol,8449) FILNAM
        icalls = icalls + 1
        return
    ! End File IO Error

END SUBROUTINE IRI_SUB

!-----------------------------------------------------------------------        
! input:   jmag,alati,along,iyyyy,mmdd,dhour  see IRI_SUB
!          height  height in km
!          h_tec_max  =0 no TEC otherwise upper boundary for integral
!          iut     =1 for UT       =0 for LT
!          ivar    =1      altitude
!                  =2,3    latitude,longitude
!                  =4,5,6  year,month,day
!                  =7      day of year
!                  =8      hour (UT or LT)
!          vbeg,vend,vstp  variable range (begin,end,step)
! output:  a       similar to outf in IRI_SUB
!          b       similar to oarr in IRI_SUB
!
!          numstp  number of steps; maximal 1000
!-----------------------------------------------------------------------        
subroutine iri_web(jmag,jf,alati,along,iyyyy,mmdd,iut,dhour,&
                   height,h_tec_max,ivar,vbeg,vend,vstp,a,b)

    dimension   outf(20,1000),oar(100),oarr(100),a(20,1000)
    dimension   xvar(8),b(100,1000)
    logical     jf(50)

    nummax = 1000
    numstp = int((vend - vbeg)/vstp) + 1
    if(numstp > nummax) numstp = nummax

    do i=1,100
        oar(i)=b(i,1) 
    end do

    if(ivar == 1) then
        do i=1,100
            oarr(i)=oar(i) 
        end do

        xhour = dhour + iut*25.
        call IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,XHOUR,&
                     VBEG,VEND,VSTP,a,OARR)

        if (h_tec_max > 50.) then 
            call iri_tec(50.,h_tec_max,2,tec,tect,tecb)
            oarr(37) = tec
            oarr(38) = tect
        end if

        do i=1,100
            b(i,1)=oarr(i)
        end do
    else

        if(height <= 0.0) height = 100
        xvar(2) = alati
        xvar(3) = along
        xvar(4) = iyyyy
        xvar(5) = mmdd/100
        xvar(6) = mmdd - xvar(5)*100
        xvar(7) = abs(mmdd*1.)
        xvar(8) = dhour
        xvar(ivar) = vbeg
        alati = xvar(2)
        along = xvar(3)
        iyyyy = int(xvar(4))

        if(ivar == 7) then
            mmdd = -int(vbeg)
        else
            mmdd = int( xvar(5)*100 + xvar(6) )
        endif

        dhour = xvar(8) + iut*25.

        do i=1,numstp
            do iii=1,100
                oarr(iii)=b(iii,i) 
            end do
            call IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,&
                         height,height,1.,OUTF,OARR)
            if (h_tec_max < 50.) then
                call iri_tec(50.,h_tec_max,2,tec,tect,tecb)
                oarr(37) = tec
                oarr(38) = tect
            end if
            do ii=1,20
                a(ii,i)  = outf(ii,1)
            end do
            do ii=1,100
                b(ii,i) = oarr(ii)
            end do
            xvar(ivar) = xvar(ivar) + vstp
            alati = xvar(2)
            along = xvar(3)
            iyyyy = int(xvar(4))
            if(ivar == 7) then
                mmdd = -xvar(7)
            else
                mmdd = int(xvar(5)*100 + xvar(6))
            endif
            dhour = xvar(8) + iut*25.
        end do
    end if
    return
END SUBROUTINE IRI_WEB

! Format statements
11   format(1X,'*NE* HMF1 IS NOT EVALUATED BY THE FUNCTION XE2'/&
            1X,'CORR.: NO F1 REGION, B1=3, C1=0.0')
100  format(1X,'*NE* HST IS NOT EVALUATED BY THE FUNCTION XE3')
104  format('ccir',I2,'.asc')
650  format(1X,'*NE* E-REGION VALLEY CAN NOT BE MODELLED')
656  format(1X,'*NE* ABT-B0 computed with RLAT=LATI=',F6.2)
901  format(6X,'CORR.: LIN. APP. BETWEEN HZ=',F5.1,' AND HEF=',F5.1)
1144 format('ursi',I2,'.asc')
1939 format(' *Ne* User input of hmF1 is only possible for the LAY-version')
2911 format('*** IRI parameters are being calculated ***')
4689 format(1X,4E15.8)
7722 format('*NE* LAY amplitudes could not be found.')
7733 format('*NE* LAY amplitudes found with 2nd choice of HXL(1).')
8449 format(1X////,' The file ',A30,'is not in your directory.')
9012 format('Ne, E-F: The LAY-Version is prelimenary.', &
            ' Erroneous profile features can occur.')
9014 format('Ne: No upper limit for F10.7 in', &
            ' topside formula.')
9204 format('Ne: Corrected Topside Formula')
9205 format('Ne: NeQuick for Topside')
9206 format('Ne: Gul-h0.5 for Topside')
9207 format('Ne: IRI-2001 for Topside')
9214 format('Ne: B0,B1 Bil-2000')
9215 format('Ne: B0 Gul-1987')
9216 format('Ne: B0,B1-ABT-2009')
9015 format('Ne, foF2/NmF2: provided by user.')
9016 format('Ne, foF2: URSI model is used.')
9017 format('Ne, foF2: CCIR model is used.')
9018 format('Ne, hmF2/M3000F2: provided by user.')
9019 format('Ne, foF1/NmF1: provided by user.')
9021 format('Ne, hmF1: provided by user.')
9022 format('Ne, foE/NmE: provided by user.')
9023 format('Ne, hmE: provided by user.')
9024 format('Ne, foF1: probability function used.')
9025 format('Ne, foF1: L condition cases included.')
9026 format('Ne, D: IRI1990')
9027 format('Ne, D: FT2001; IRI-90, FT-01, DRS-95)')
9028 format('Ne, foF2: Storm model turned off if foF2 or', &
            ' NmF2 user input')
9029 format('Ne, foF2: storm model included')
9128 format('Ne, foE: storm model on')
9129 format('Ne, foE: storm model off')
9039 format('Ion Com.: DS-95 & DY-85')
9031 format('Ion Com.: RBV-10 & TTS-03')
9032 format('Te: Temperature-density correlation is used.')
9033 format('Te: Aeros/AE/ISIS model')
9034 format('Te: TBT-2012 model')
4031 format('Auroral boundary model on')
4032 format('Auroral boundary model off')
