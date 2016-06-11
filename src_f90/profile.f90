!************************************************************                   
!***************** PROFILE PARAMETERS ***********************                   
!************************************************************                 


SUBROUTINE TOPH05(COVI,AMLAT,TIME,HMAX,HT05,SG)
    !---------------------------------------------------------------------------------
    ! Gulyaeva T.L. (2003) Variations in the half-width of the topside ionosphere 
    !    according to the observations by space ionosondes ISIS 1,ISIS 2, and IK19.
    !    International J. of Geomagnetism and Aeronomy, 4(3), 201-207.
    ! Gulyaeva T.L., Titheridge J.E. (2006) Advanced specification of electron density 
    !    and temperature in the IRI ionosphere-plasmasphere model. 
    !    Adv. Space Res. 38(11), 2587-2595, doi:10.1016/j.asr.2005.08.045.
    !
    !  Implementation of empirical RAT=(h05top-hmF2)/hmF2 derived from ISIS and IK19
    !  topside electron density profiles to obtain half peak density topside height
    !  h05top  from the Chebishev polinomial coefficients given for 
    !  (1) 4 levels of solar activity: Rz= 0,  50, 100, 150 replaced by
    !      solar radio flux          covi=60, 106, 152, 198
    !  (2) 10 selected grids of geomagnetic latitude (N=S):0,10,20,30,40,50,60,70,80,90
    !  (3) 5 selected grids of local time: 0, 6, 12, 18, 24.
    !  (4) 4 seasonal grids: 1 equinox(SG=90deg), 2 summer (SG=180), 
    !                        3 equinox (SG=270), 4 winter(SG=360)
    !   SG=season grids=90,180,270,360
    !---------------------------------------------------------------------------------
    DIMENSION CVLEV(4)    
    COMMON   /BLOCK1/HMF2,XNMF2,XHMF1,F1REG         &
        /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU
    DATA CVLEV/60.,106.,152.,198./
    LOGICAL F1REG
    ABMLAT=ABS(AMLAT)
    IR=IFIX((covi-60.)/46.)+1    
    M1=IFIX(ABMLAT/10.)+1
    L1=IFIX(TIME/6.)+1
    M2=M1+1
    IF(M1 == 10) M2=10
    L2=L1+1
    IF(L1 == 5) L2=5
    !
    ! INTERPOLATE RAT FOR GIVEN RZI
    ! Call Chebishev approximation to interpolate for given ABMLAT, HRLT
    !
    call CHEBISH(CVLEV(IR),TIME,ABMLAT,XX,SG)
    IF (IR == 4) THEN
        RAT05=XX
        GOTO 10
    ENDIF
    call CHEBISH(CVLEV(IR+1),TIME,ABMLAT,YY,SG)
    RAT05=XX+(YY-XX)*(COVI-CVLEV(IR))/46.
    10    HT05=HMAX*(1.+RAT05)
    RETURN
END
!
!      
SUBROUTINE CHEBISH(COVS,HOURLT,ABMLAT,RATCH,SG)
    !---------------------------------------------------------------------------------
    ! CHEBISHEV POLINOMIALS FOR ABMLAT(10),HOURLT(5)
    ! CR((C0...C5),(LT=0,6,...24),(SG=season grids=90,180,270,360)
    !                            (COV=60,106,152,198)
    !---------------------------------------------------------------------------------
    !      REAL UK(0:10),CR(0:5,5,3,4),YI(5),YY(5,3)
    REAL BR(6,5,3,4),YI(5),YY(5,3)
    REAL PL1(5),PL2(5),PL3(5),CL(0:3)
    !
    COMMON /CONST/rad,pi  
    DATA PL1/-2.,-1.,0.,1.,2./
    DATA PL2/2.,-1.,-2.,-1.,2./
    DATA PL3/-1.,2.,0.,-2.,1./
    DATA BR/
    ! Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=60 (mlat/10=0,1,...,9)
    ! Equinox    B0MLAT:&
        -1.5859,  3.5789, -3.7884, 2.7094,-1.2962,.6759&
        , -5.3449, 12.8379,-12.0165, 5.9746,-1.9084,.7669&
        ,-12.8000, 35.3084,-38.0043,19.6004,-4.4974,.6975&
        ,  5.8282,-13.3538,  9.1674,-0.9593,-0.8909,.6062&
        , -1.5859,  3.5789, -3.7884, 2.7094,-1.2962,.6759
    ! Summer    B0MLAT    &
        , -7.1103, 21.0389,-24.5539,14.1607,-3.8537,.7266&
        ,  5.5333,-10.6242,  4.8751, 1.0587,-1.0821,.7527&
        ,-15.4487, 42.9269,-45.0314,21.4718,-4.2116,.6026&
        , -6.6436, 16.4533,-15.5142, 6.8287,-1.2871,.4976&
        , -7.1103, 21.0389,-24.5539,14.1607,-3.8537,.7266
    ! Winter   B0MLAT                                                                       &
        , 14.9103,-35.2337, 27.3078,-6.5362,-0.6265,.7509&
        ,  2.3846, -5.8840,  3.7023, 0.8525,-1.2663,.7086&
        , -9.8846, 26.6649,-27.0173,12.6959,-2.6536,.6295&
        ,  1.7692, -2.3578, -0.7945, 2.2477,-0.9691,.5719&
        , 14.9103,-35.2337, 27.3078,-6.5362,-0.6265,.7509
    ! Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=106 (mlat=0,10,...,90)
    ! Equinox    B1MLAT&
        , -4.1218, 10.6136,-11.4922, 6.0470,-1.3620,.5563&
        ,  0.9077,  2.9562, -8.9880, 6.8680,-1.9621,.7737&
        ,-16.2744, 42.8047,-43.7009,20.7965,-4.0697,.6619&
        ,-17.3038, 44.3336,-40.9249,15.9042,-2.1554,.4796&
        , -4.1218, 10.6136,-11.4922, 6.0470,-1.3620,.5563
    ! Summer    B1MLAT  &
        , -4.9692, 16.5753,-21.3543,12.7061,-3.1758,.6446&
        ,  1.9000, -2.8167, -0.9962, 3.0687,-1.3454,.6859&
        ,  7.6769,-14.8343,  6.7030, 1.5578,-1.0626,.4291&
        ,  5.4833,-10.6322,  4.7571, 1.2178,-0.8223,.4615&
        , -4.9692, 16.5753,-21.3543,12.7061,-3.1758,.6446
    ! Winter    B1MLAT  &
        , -4.7282, 13.4491,-15.6931, 8.8388,-1.9732,.5874&
        ,  5.6756,-14.8458, 11.8927,-2.2632,-0.6122,.6948&
        ,-14.2872, 40.0829,-41.2716,18.1696,-2.7203,.4916&
        ,-13.6128, 33.4657,-29.7231,11.0972,-1.2884,.5034&
        , -4.7282, 13.4491,-15.6931, 8.8388,-1.9732,.5874
    ! Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=152 (mlat=0,10,...,90)
    ! Equinox    B2MLAT&
        , -3.3282, 10.4296,-12.4722, 6.7623,-1.5172,.4931&
        , -8.9744, 20.1311,-17.4552, 7.6518,-1.7371,.6702&
        , 12.0462,-27.8932, 20.6241,-4.5781, 0.0814,.3501&
        ,-17.0551, 42.3258,-37.1874,13.3608,-1.4804,.4216&
        , -3.3282, 10.4296,-12.4722, 6.7623,-1.5172,.4931
    ! Summer    B2MLAT  &
        ,  7.3077,-17.1579, 11.6872,-0.7405,-1.0298,.5754&
        , 19.2641,-45.1886, 34.3297,-8.1879,-0.1875,.6562&
        ,  6.0987,-11.0903,  4.3569, 1.4001,-0.7309,.3885&
        ,  5.9295,-13.9205, 10.2347,-2.2818, 0.0853,.3915&
        ,  7.3077,-17.1579, 11.6872,-0.7405,-1.0298,.5754
    ! Winter    B2MLAT  &
        , -1.6821,  8.6010,-13.6570, 8.6307,-1.9846,.5635&
        ,  5.4679,-12.3750,  7.5620, 0.5394,-1.4415,.6659&
        , -8.0821, 21.9288,-21.8597, 9.3455,-1.4644,.3599&
        , -8.3000, 19.3076,-16.3295, 6.1619,-0.9144,.3846&
        , -1.6821,  8.6010,-13.6570, 8.6307,-1.9846,.5635
    ! Polinomial Coefficients B1,B2,B3,B4,B5,B6 for COV=198 (mlat=0,10,...,90)
    ! Equinox    B3MLAT&
        ,-16.4051, 28.2880,-16.0982, 4.6328,-1.0405,.5486&
        , 13.0846,-34.8291, 30.0074,-8.6402, 0.1529,.6165&
        , 19.7474,-42.7116, 28.9430,-6.0487, 0.1492,.3748&
        , 16.2795,-36.6982, 26.5094,-6.3492, 0.2926,.3946&
        ,-16.4051, 28.2880,-16.0982, 4.6328,-1.0405,.5486
    ! Summer    B3MLAT&
        ,  4.6410,-13.7931, 11.6548,-1.9248,-0.7246,.5264&
        , -2.4090,  3.1805, -2.8423, 2.8861,-0.9937,.5234&
        ,  6.3410,-13.9643,  8.2461,-0.0186,-0.7009,.3582&
        ,  9.0987,-20.8618, 14.7262,-2.8798,-0.0512,.3662&
        ,  4.6410,-13.7931, 11.6548,-1.9248,-0.7246,.5264
    ! Winter    B3MLAT&
        , -4.6526, 12.1878,-14.4047, 8.5226,-2.0493,.5903&
        ,  3.9821, -6.9477,  0.8382, 3.4898,-1.5694,.6283&
        , -7.0474, 17.3974,-17.3465, 8.3671,-1.5708,.3759&
        ,  4.2782, -9.9880,  5.9834, 0.0975,-0.4900,.3842&
        , -4.6526, 12.1878,-14.4047, 8.5226,-2.0493,.5903/
    !    DATA UL/-2.,-1.,0.,1.,2./
    do k=0,3
        cl(k)=0.
    enddo
    !
    IR=IFIX((covs-60.)/46.)+1    
    ! Given geomagnetic latitude parameter:
    xi=abmlat/100.
    DO LS=1,3
        DO LL=1,5
            B1=BR(6,LL,LS,IR)      
            B2=BR(5,LL,LS,IR)        
            B3=BR(4,LL,LS,IR)       
            B4=BR(3,LL,LS,IR)        
            B5=BR(2,LL,LS,IR)        
            B6=BR(1,LL,LS,IR)        
            HLT=(LL-1)*6.0
            YY(LL,LS)=B1+xi*(B2+xi*(B3+xi*(B4+xi*(B5+xi*B6))))
        ENDDO
    ENDDO            ! end of season/day cycle
    ! Apply seasonal interpolation
    do i=1,5
        p0=(2.*YY(i,1)+YY(i,2)+YY(i,3))/4.
        p1=(YY(i,3)-YY(i,2))/2.
        p2=(YY(i,2)+YY(i,3)-2.*YY(i,1))/4.
        YI(i)=p0+p1*cos(sg*rad)+p2*cos(2.*sg*rad)
    enddo
    DO K=1,5
        CL(0)=CL(0)+YI(K)
        CL(1)=CL(1)+YI(K)*PL1(K)
        CL(2)=CL(2)+YI(K)*PL2(K)
        CL(3)=CL(3)+YI(K)*PL3(K)
    ENDDO
    CL(0)=CL(0)/5.
    CL(1)=CL(1)/10.
    CL(2)=CL(2)/14.
    CL(3)=CL(3)/12.
    ULL=(HOURLT-12.)/6.
    ZA=CL(0)-2.*CL(2)
    RATCH=ZA+ULL*(CL(1)-3.4*CL(3)+ULL*(CL(2)+ULL*CL(3)))
    RETURN
END    
!
!
SUBROUTINE  SHAMDB0D (RLAT,FLON,T,RZ,B)
    !-------------------------------------------------------------------
    !    COMPUTES THE HOURLY VALUES OF B0 FROM A SET OF SH COEFFICIENTS
    !    IN A POINT OF A GIVEN GEOCENTRIC LATITUDE AND LONGITUDE
    !    OF THE EARTH'S SURFACE FOR A GIVEN MONTH AND A GIVEN SUSPOT NUMER
    !
    ! INPUT:    RLAT    The geogrphic latitude on the meridian given by 
    !                    the local time (FLON), where the modified dip
    !                   latitude is the same as of the orginal site.
    !            FLON    =LONGITUDE+15.*UT(hours)
    !            T        Month as a REAL number (1.0 to 12.0)
    !            RZ        12-month running mean
    ! OUTOUT    B        =B0
    !
    !  Blanch E., D. Arrazola, D. Altadill, D. Buresova, M. Mosert, 
    !     Adv. Space Res. 39, 701-710, 2007.
    !  Altadill, D., D. Arrazola, E. Blanch, D. Buresova, 
    !     Adv. Space Res. 42, 610-616, 2008.
    !  Altadill, D., J.M. Torta, and E. Blanch, 
    !     Adv. Space Res. 43,1825-1834, 2009.
    !-------------------------------------------------------------------
    PARAMETER   (IBO=0,JBO=1,KDIM=6,LDIM=4,L=-1)
    DIMENSION   FN(0:KDIM,0:KDIM), CONST(0:KDIM,0:KDIM)
    DIMENSION   GNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        HNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
    DIMENSION   GANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        HANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        GBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        HBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
    DIMENSION   BINT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        BEXT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
    CHARACTER*1 IE       
    COMMON/ATB/BINT,BEXT,RE,TZERO,IFIT,IB,KINT,LINT,KEXT,&
        LEXT,KMAX,FN
    DATA THETA,RE,TZERO,IFIT,ICEN,IREF,IB,KINT,LINT,KEXT,LEXT&
        /180.,6371.2,1.0,  -1,  0,   0,   2, 6,   4,   0,   -1/
    DATA ((CONST(N,M), M=0,N), N=0,KDIM)&
        /4*1.,1.73205,0.866025,1.,2.44949,1.93649,0.7905691,1.,&
        3.16228,3.35410,2.09165,0.739510,1.,3.87298,5.12348,&
        4.18330,2.21853,0.701561,1.,4.58258,2*7.24569,4.96078,&
        2.32681,0.671693/
    DATA IE /"I"/
    DATA (((GANM(N,M,J),GBNM(N,M,J),HANM(N,M,J),HBNM(N,M,J),&
        J=0,LDIM), M=0,N), N=0,KDIM)&
        /176.746, 0.233,   0.000, 0.000,&
        -109.413, 0.072,   0.000, 0.000, -66.000,-0.004,   0.000, 0.000,&
        36.874,-0.018,   0.000, 0.000, -19.515, 0.040,   0.000, 0.000,&
        94.998, 0.724,   0.000, 0.000,&
        -116.900,-0.971,   0.000, 0.000, -93.849,-0.590,   0.000, 0.000,&
        80.579, 0.425,   0.000, 0.000, -19.205,-0.220,   0.000, 0.000,&
        -94.824, 0.115,   6.055, 0.265,&
        84.720,-0.161,  -7.101,-0.374,  35.200,-0.138,   1.043,-0.350,&
        -23.960, 0.109,  -2.478, 0.133,  25.550,-0.049,  -3.143, 0.003,&
        -29.418,-0.823,   0.000, 0.000,&
        40.070, 0.800,   0.000, 0.000,  24.019, 0.478,   0.000, 0.000,&
        -13.175,-0.265,   0.000, 0.000,   8.799, 0.090,   0.000, 0.000,&
        -31.099,-1.117,  -1.906, 0.498,&
        43.807, 1.406,  -3.216,-0.520,  37.957, 0.902,  -0.717,-0.570,&
        -40.232,-0.583,  12.171, 0.244,  11.595, 0.241,  -4.890, 0.054,&
        -87.665, 1.635, 117.581,-1.765,&
        123.444,-2.119,-146.917, 2.131,  81.137,-1.173, -99.063, 1.548,&
        -42.646, 0.681,  61.263,-0.811,  17.550,-0.408, -24.374, 0.260,&
        54.538,-0.170,   0.000, 0.000,&
        -71.552, 0.361,   0.000, 0.000, -50.565,-0.077,   0.000, 0.000,&
        36.653,-0.071,   0.000, 0.000, -10.816, 0.236,   0.000, 0.000,&
        -31.138, 1.156,  37.307,-1.407,&
        40.390,-1.390, -34.573, 1.730,  41.597,-0.835, -41.318, 1.550,&
        -19.779, 0.404,  15.954,-0.696,  -1.706,-0.220,   5.084, 0.040,&
        -57.671, 0.045,  42.683,-0.800,&
        71.491, 0.048, -49.441, 0.980,  47.893,-0.037, -36.191, 0.562,&
        -26.638,-0.029,  20.346,-0.384,   9.998, 0.067,  -6.787, 0.213,&
        90.187,-1.198, -66.032,-0.056,&
        -119.148, 1.428,  81.202, 0.022, -63.375, 0.754,  53.070, 0.149,&
        39.249,-0.436, -30.898,-0.052, -27.293, 0.301,  12.838,-0.067,&
        -110.261, 1.509,   0.000, 0.000,&
        164.956,-1.761,   0.000, 0.000, 103.699,-1.005,   0.000, 0.000,&
        -55.127, 0.569,   0.000, 0.000,  25.376,-0.315,   0.000, 0.000,&
        -104.655, 1.341, 109.057,-1.367,&
        139.129,-1.730,-127.325, 1.532,  88.526,-1.068,-106.461, 1.397,&
        -38.306, 0.508,  56.240,-0.798,  17.239,-0.267,  -7.766, 0.058,&
        -6.494,-1.253,   5.714, 0.132,&
        3.547, 1.545,  -5.372,-0.106,  -4.343, 1.103,  -3.393,-0.017,&
        2.454,-0.626,  -3.297,-0.025,   5.871, 0.160,   2.040,-0.036,&
        50.814,-0.230, -25.094, 0.817,&
        -65.502, 0.304,  32.267,-1.075, -44.176, 0.019,  14.606,-0.605,&
        27.869,-0.009,  -5.147, 0.387, -11.041, 0.131,   5.922,-0.225,&
        77.825,-0.728, 128.501,-0.810,&
        -87.685, 0.838,-164.016, 1.103, -74.431, 0.807, -95.539, 0.498,&
        40.631,-0.454,  49.950,-0.292,  -4.229, 0.000, -29.666, 0.272,&
        152.380,-1.232,   0.000, 0.000,&
        -192.098, 1.514,   0.000, 0.000,-132.417, 1.370,   0.000, 0.000,&
        82.894,-0.709,   0.000, 0.000, -28.162, 0.050,   0.000, 0.000,&
        -12.633, 1.192,  47.246,-1.193,&
        -5.488,-1.387, -67.206, 1.486,  -9.917,-0.914, -34.438, 0.552,&
        13.185, 0.477,  21.225,-0.387,   0.586,-0.208, -15.426, 0.419,&
        -4.478,-0.118,  17.908, 0.175,&
        -0.417, 0.067, -27.047,-0.241,   7.636, 0.028, -10.075,-0.109,&
        -10.582, 0.005,  14.496, 0.086,   0.421, 0.001, -12.200,-0.041,&
        16.086, 0.321,  47.044,-0.126,&
        -24.823,-0.280, -62.615, 0.210, -12.030,-0.136, -44.003,-0.023,&
        4.929, 0.137,  28.340,-0.009,  -4.688,-0.057,  -9.315, 0.103,&
        28.023,-0.031, -21.535, 0.115,&
        -31.946, 0.011,  24.143,-0.180, -21.019,-0.057,  24.108,-0.116,&
        13.969, 0.004, -13.823, 0.042,  -6.860, 0.031,   0.546,-0.035,&
        20*0.000,&
        -31.994, 0.409,   0.000, 0.000,&
        18.217,-0.458,   0.000, 0.000,  39.280,-0.754,   0.000, 0.000,&
        -20.453, 0.324,   0.000, 0.000,  -8.111, 0.139,   0.000, 0.000,&
        28.765,-0.477, -28.368, 0.516,&
        -50.604, 0.751,  25.725,-0.471, -23.444, 0.283,  29.966,-0.558,&
        2.759,-0.146, -10.824, 0.341,  -7.419, 0.206,  -3.711, 0.056,&
        42.429,-0.415,   1.993, 0.117,&
        -53.162, 0.555,   7.229,-0.246, -19.307, 0.039,  -8.028, 0.028,&
        9.849,-0.035,   6.834, 0.033, -17.010, 0.272,   4.668,-0.129,&
        4.546,-0.359, -57.796, 0.359,&
        0.738, 0.343,  73.027,-0.423,  -7.421, 0.420,  56.067,-0.327,&
        5.093,-0.279, -37.581, 0.226,   3.636,-0.041,  10.910,-0.059,&
        88.440,-0.393, -69.598, 0.643,&
        -109.481, 0.532,  82.266,-0.765, -59.229, 0.182,  55.279,-0.580,&
        28.514,-0.057, -30.282, 0.326, -22.924, 0.164,  11.602,-0.073,&
        40*0.000/                                
    KMAX = MAX(KINT,KEXT)
    IF (KMAX  >  KDIM)  GO TO 9999
    KT = MAX(LINT,LEXT)
    IF (KT  >  LDIM)  GO TO 9999
    DO N=0,KMAX
        DO M=0,N
            DO J=0,KT
                GNM(N,M,J)=GANM(N,M,J)+GBNM(N,M,J)*rz
                HNM(N,M,J)=HANM(N,M,J)+HBNM(N,M,J)*rz
            ENDDO
            
            IF (IE  ==  'I')  THEN
                IF (N  >  KINT)  GO TO 500
                LJ = LINT
            ELSE
                IF (N  >  KEXT)  GO TO 500
                LJ = LEXT
            END IF
            FN(N,M) = FLOAT(N)
            IF (M  >  0)  GO TO 300
            DO J=1-IBO-JBO,KT
                IF (IE  ==  'I')  THEN
                    BINT(N,M,J)   = GNM(N,M,J)
                ELSE
                    BEXT(N,M,J)   = GNM(N,M,J)
                END IF
            end do
            GO TO 500
            300 continue
            DO J=1-IBO-JBO,LJ
                IF (IE  ==  'I')  THEN
                    BINT(N,M,J)   = GNM(N,M,J)
                    BINT(M-1,N,J) = HNM(N,M,J)
                ELSE
                    BEXT(N,M,J)   = GNM(N,M,J)
                    BEXT(M-1,N,J) = HNM(N,M,J)
                END IF
            end do
            !
        end do
    end do
    500 CONTINUE
    !     **********************************************************
    !     SYNTHESIZES THE VALUE OF B0 FROM THE MODEL    
    !     **********************************************************
    CALL SCHNEVPD(RZ,RLAT,FLON,dum,T,L,dum,dum,B)
    RETURN
    9999  STOP
END
!
!
SUBROUTINE  SHAB1D (FLAT,FLON,T,RZ,B)
    !-------------------------------------------------------------------
    !    COMPUTES THE HOURLY VALUES OF B1 FROM A SET OF SH COEFFICIENTS
    !    IN A POINT OF A GIVEN GEOCENTRIC LATITUDE AND LONGITUDE
    !    OF THE EARTH'S SURFACE FOR A GIVEN MONTH AND A GIVEN SUSPOT NUMER
    !
    !   PARAMETERS ARE THE SAME AS IN SHAMDB0D, EXCEPT:
    !        FLAT    Geographic latitude
    !        B        =B1
    !
    !    ***** PARAMS   COEFFS IN DATA SENTENCES *****
    !-------------------------------------------------------------------
    !
    PARAMETER   (IBO=0,JBO=1,KDIM=6,LDIM=4,L=-1)
    DIMENSION   FN(0:KDIM,0:KDIM), CONST(0:KDIM,0:KDIM)
    DIMENSION   GNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        HNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
    DIMENSION   GANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        HANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        GBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        HBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
    DIMENSION   BINT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        BEXT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
    CHARACTER*1 IE       
    COMMON/ATB/BINT,BEXT,RE,TZERO,IFIT,IB,KINT,LINT,KEXT,&
        LEXT,KMAX,FN
    DATA ALT /300./
    DATA THETA,RE,TZERO,IFIT,ICEN,IREF,IB,KINT,LINT,KEXT,LEXT&
        /180.,6371.2,1.0,  -1,  0,   0,   2, 6,   4,   0,   -1/
    DATA ((CONST(N,M), M=0,N), N=0,KDIM)&
        /4*1.,1.73205,0.866025,1.,2.44949,1.93649,0.7905691,1.,&
        3.16228,3.35410,2.09165,0.739510,1.,3.87298,5.12348,&
        4.18330,2.21853,0.701561,1.,4.58258,2*7.24569,4.96078,&
        2.32681,0.671693/
    DATA IE /"I"/
    DATA (((GANM(N,M,J),GBNM(N,M,J),HANM(N,M,J),HBNM(N,M,J),&
        J=0,LDIM), M=0,N), N=0,KDIM)&
        /1.156, 0.039,  0.000, 0.000,&
        1.725,-0.053,  0.000, 0.000,  1.097,-0.032,  0.000, 0.000,&
        -0.579, 0.019,  0.000, 0.000,  0.265,-0.010,  0.000, 0.000,&
        -2.895, 0.023,  0.000, 0.000,&
        3.269,-0.025,  0.000, 0.000,  2.278,-0.015,  0.000, 0.000,&
        -1.789, 0.009,  0.000, 0.000,  0.653,-0.005,  0.000, 0.000,&
        -3.240, 0.052, -2.645, 0.030,&
        4.404,-0.062,  3.283,-0.038,  2.827,-0.038,  2.181,-0.019,&
        -1.496, 0.020, -1.143, 0.011,  0.688,-0.011,  0.512,-0.008,&
        -0.023,-0.025,  0.000, 0.000,&
        -0.370, 0.031,  0.000, 0.000, -0.385, 0.017,  0.000, 0.000,&
        0.150,-0.009,  0.000, 0.000,  0.019, 0.007,  0.000, 0.000,&
        1.704, 0.006, -1.766, 0.022,&
        -2.115,-0.007,  2.309,-0.028, -1.610,-0.003,  1.314,-0.015,&
        0.806, 0.002, -0.873, 0.010, -0.249,-0.002,  0.445,-0.005,&
        0.256,-0.009, -2.959, 0.023,&
        -0.662, 0.013,  3.643,-0.030, -0.208, 0.001,  2.208,-0.015,&
        0.144,-0.001, -1.326, 0.010, -0.205, 0.005,  0.678,-0.007,&
        1.730,-0.030,  0.000, 0.000,&
        -2.072, 0.035,  0.000, 0.000, -1.027, 0.016,  0.000, 0.000,&
        0.946,-0.008,  0.000, 0.000, -0.644, 0.009,  0.000, 0.000,&
        3.925,-0.060,  1.613,-0.015,&
        -4.790, 0.070, -2.021, 0.019, -3.293, 0.048, -1.273, 0.006,&
        1.829,-0.026,  0.797,-0.005, -0.767, 0.011, -0.395, 0.006,&
        -1.988, 0.027, -0.761, 0.003,&
        2.258,-0.031,  0.978,-0.004,  1.772,-0.026,  0.459,-0.003,&
        -0.813, 0.014, -0.257, 0.004,  0.173,-0.003,  0.223,-0.001,&
        20*0.000,&
        -1.356, 0.011,  0.000, 0.000,&
        0.079,-0.003,  0.000, 0.000,  0.415,-0.012,  0.000, 0.000,&
        -0.249, 0.004,  0.000, 0.000, -0.087, 0.004,  0.000, 0.000,&
        -1.155,-0.012,  0.261,-0.026,&
        1.619, 0.011, -0.217, 0.032,  0.989, 0.009, -0.361, 0.019,&
        -0.745,-0.001,  0.009,-0.009,  0.347, 0.000,  0.178, 0.005,&
        4.672,-0.032,  0.562, 0.017,&
        -5.868, 0.040, -0.850,-0.020, -3.798, 0.026, -0.145,-0.020,&
        2.079,-0.014,  0.160, 0.008, -0.943, 0.007, -0.417, 0.001,&
        40*0.000,&
        2.477, 0.000,  0.000, 0.000,&
        -1.815,-0.011,  0.000, 0.000, -1.571, 0.002,  0.000, 0.000,&
        0.551, 0.004,  0.000, 0.000, -0.044,-0.008,  0.000, 0.000,&
        2.160,-0.010, -0.305, 0.012,&
        -2.618, 0.010,  0.529,-0.016, -1.782, 0.006,  0.634,-0.007,&
        0.976,-0.003, -0.449, 0.006, -0.331, 0.001, -0.004,-0.004,&
        -0.394, 0.002,  0.851,-0.020,&
        0.359,-0.003, -1.051, 0.024,  0.357, 0.002, -0.239, 0.012,&
        0.005,-0.001,  0.210,-0.008, -0.028,-0.003, -0.322, 0.005,&
        60*0.000,&
        -6.760, 0.064,  0.000, 0.000,&
        7.700,-0.073,  0.000, 0.000,  5.394,-0.054,  0.000, 0.000,&
        -2.788, 0.026,  0.000, 0.000,  0.923,-0.007,  0.000, 0.000,&
        -2.328, 0.024, -0.463, 0.020,&
        2.923,-0.027,  0.490,-0.025,  1.768,-0.019,  0.711,-0.017,&
        -1.068, 0.009, -0.363, 0.010,  0.596,-0.004, -0.073,-0.004,&
        -1.911, 0.016, -4.519, 0.041,&
        2.644,-0.024,  5.569,-0.050,  1.287,-0.009,  3.707,-0.031,&
        -0.894, 0.007, -2.121, 0.019,  0.669,-0.007,  0.933,-0.010,&
        80*0.000/ 
    KMAX = MAX(KINT,KEXT)
    IF (KMAX  >  KDIM)  GO TO 9999
    KT = MAX(LINT,LEXT)
    IF (KT  >  LDIM)  GO TO 9999
    DO N=0,KMAX
        DO M=0,N
            DO J=0,KT
                GNM(N,M,J)=GANM(N,M,J)+GBNM(N,M,J)*rz
                HNM(N,M,J)=HANM(N,M,J)+HBNM(N,M,J)*rz
            ENDDO
            
            IF (IE  ==  'I')  THEN
                IF (N  >  KINT)  GO TO 500
                LJ = LINT
            ELSE
                IF (N  >  KEXT)  GO TO 500
                LJ = LEXT
            END IF
            FN(N,M) = FLOAT(N)
            IF (M  >  0)  GO TO 300
            255 FORMAT (1X,A1,2I3,F9.4,E15.6,F10.3,4F20.3:/(22X,5F20.3))
            DO J=1-IBO-JBO,KT
                IF (IE  ==  'I')  THEN
                    BINT(N,M,J)   = GNM(N,M,J)
                ELSE
                    BEXT(N,M,J)   = GNM(N,M,J)
                END IF
            end do
            GO TO 500
            300 continue
            260 FORMAT (1X,A1,2I3,F9.4,E15.6,10F10.3:/(32X,10F10.3))
            DO J=1-IBO-JBO,LJ
                IF (IE  ==  'I')  THEN
                    BINT(N,M,J)   = GNM(N,M,J)
                    BINT(M-1,N,J) = HNM(N,M,J)
                ELSE
                    BEXT(N,M,J)   = GNM(N,M,J)
                    BEXT(M-1,N,J) = HNM(N,M,J)
                END IF
            end do
        end do
    end do
    500 CONTINUE
    !     **********************************************************
    !     SYNTHESIZES THE VALUE OF B1 FROM THE MODEL    
    !     **********************************************************
    CALL SCHNEVPD(RZ,FLAT,FLON,dum,T,L,dum,dum,B)
    !
    RETURN
    9999  STOP
END
!
!
SUBROUTINE SCHNEVPD (RZ,FLAT,FLON,R,T,L,BN,BE,BV)
    !-------------------------------------------------------------------
    !     WHEN L IS POSITIVE:
    !     COMPUTES SPHERICAL CAP HARMONIC (GEOCENTRIC) FIELD COMPONENTS
    !     HORIZONTAL NORTH BN,HORIZONTAL EAST BE,AND VERTICAL DOWNWARD BV.
    !     WHEN L IS NEGATIVE:
    !     COMPUTES GENERAL FUNCTION BV, ITS HORIZONTAL NORTH DERIVATIVE BN,
    !     AND ITS HORIZONTAL EAST DERIVATIVE BE, ON SPHERICAL CAP SURFACE.
    !     NOTE THAT THESE ARE METRICAL DERIVATIVES, AND BE IS THE
    !     LONGITUDINAL DERIVATIVE DIVIDED BY SIN(COLATITUDE).
    !     FLAT,FLON,R ARE GEOCENTRIC SPHERICAL CAP LATITUDE,LONGITUDE,RADIAL
    !     DISTANCE; T IS TIME.
    !     L =  0  ON FIRST CALL:  RETURNS SPHERICAL CAP POLE POSITION FLATO,FLONO
    !             AND HALF-ANGLE THETA AS BN,BE, AND BV AFTER INITIALIZATION.
    !             ON SUBSEQUENT CALLS:  ACTS AS L=1.
    !          1  COMPUTES POTENTIAL FIELD COMPONENTS FROM INTERNAL COEFFICIENTS.
    !          2  COMPUTES POTENTIAL FIELD COMPONENTS FROM EXTERNAL COEFFICIENTS.
    !          3  COMPUTES FIELD FROM BOTH INTERNAL AND EXTERNAL COEFFICIENTS.
    !         -1  COMPUTES GENERAL FUNCTION BV AND DERIVATIVES BN WITH RESPECT TO
    !             LATITUDE AND BE WITH RESPECT TO LONGITUDE DIVIDED BY COS(LAT)
    !             (R IS DUMMY VARIABLE IN THIS CASE).
    !     NOTE:   SUBROUTINE IS INITIALIZED DURING FIRST CALL REGARDLESS OF L.
    !     SUBPROGRAM USED:  LEGFUN
    !    ***** PARAMS   COEFFS TRANSFERRED FROM MAIN PROGRAM *****
    !    ADAPTED FROM SUBROUTINE SCHNEV OF G.V. HAINES (COMPUTERS   GEOSCIENCES, 
    !      14, 413-447, 1988)
    !-------------------------------------------------------------------
    PARAMETER   (IBO=0,JBO=1,KDIM=6,LDIM=4)                                     
    DIMENSION   FN(0:KDIM,0:KDIM), CONST(0:KDIM,0:KDIM)
    DIMENSION   CML(KDIM), SML(KDIM)
    DIMENSION   DELT(0:LDIM)
    DIMENSION   BINT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        BEXT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
    COMMON/ATB/BINT,BEXT,RE,TZERO,IFIT,IB,KINT,LINT,KEXT,&
        LEXT,KMAX,FN
    CHARACTER*1 IE,RESP
    !
    DATA ((CONST(N,M), M=0,N), N=0,KDIM)&
        /4*1.,1.73205,0.866025,1.,2.44949,1.93649,0.7905691,1.,&
        3.16228,3.35410,2.09165,0.739510,1.,3.87298,5.12348,&
        4.18330,2.21853,0.701561,1.,4.58258,2*7.24569,4.96078,&
        2.32681,0.671693/
    dfarg=(atan(1.0)*4.)/180.
    !     IBF   =  0   TO USE ORDINARY POLYNOMIALS AS BASIS FUNCTIONS
    !              1          LEGENDRE POLYNOMIALS
    !              2          FOURIER SERIES
    !              3          COSINE SERIES
    !              4          SINE SERIES
    !     NOTE:    TZERO AND THINT MAY DEPEND ON IBF.
    IBF   =  2                                                       
    T1=1.
    T2=12.
    CALL TBFIT (T1,T2,IBF,THINT,TZERO)                                
    !      IF (L  /=  0)  GO TO 100
    !      BN = FLATO
    !      BE = FLONO
    !      BV = THETA
    !      RETURN
    100 IF (L  >=  0)  THEN
    IF (IFIT  <  0)  GO TO 999
    AOR = RE/R
    AR = AOR**2
    IF (L  >  1)  GO TO 107
ELSE
    IF (IFIT  >=  0)  GO TO 999
    AR = -1.
END IF
KT = LINT
GO TO 109
107 IF (KEXT  >  0)  AOR3 = AOR*AR
IF (L  >  2)  GO TO 108
KT = LEXT
GO TO 109
108 KT = MAX (LINT,LEXT)
109 DELT(0) = 1.
IF (KT  <=  0)  GO TO 103
DEL = (T - TZERO)/THINT
DO I=1,KT
    IF (I  ==  1)  THEN
        IF (IBF  <=  1) THEN
            DELT(I) = DEL
        ELSE IF (IBF  ==  2)  THEN
            ST = SIN(DEL)
            DELT(I) = ST
        ELSE IF (IBF  ==  3)  THEN
            DELT(I) = COS(DEL)
        ELSE
            DELT(I) = SIN(DEL)
        ENDIF
        GO TO 102
    ENDIF
    IF (IBF  ==  0)  THEN
        DELT(I) = DELT(I-1)*DEL
    ELSE IF (IBF  ==  1)  THEN
        RECIP = 1./FLOAT(I)
        DELT(I) = (2.-RECIP)*DELT(I-1)*DEL - (1.-RECIP)*DELT(I-2)
    ELSE IF (IBF  ==  2)  THEN
        IF ((I/2)*2  ==  I)  THEN
            IF (I  ==  2)  THEN
                CT = COS(DEL)
                DELT(I) = CT
            ELSE
                DELT(I) = DELT(I-2)*CT - DELT(I-3)*ST
            ENDIF
        ELSE
            DELT(I) = DELT(I-2)*CT + DELT(I-1)*ST
        ENDIF
    ELSE IF (IBF  ==  3)  THEN
        DELT(I) = COS(I*DEL)
    ELSE IF (IBF  ==  4)  THEN
        DELT(I) = SIN(I*DEL)
    ELSE
        GO TO 999
    ENDIF
end do
102 CONTINUE
incept = 0                                                              
if ((ibf == 2 .or. ibf == 3) .and. incept  ==  1)  then
    !     change to intercept form of fourier series.
    do i=2,lint,4-ibf
        delt(i) = 1. - delt(i)
    enddo
endif
103 X = 0.
Y = 0.
Z = 0.
IF (L  ==  2)  GO TO 106
IF (KINT  <  0)  GO TO 106
GTI = 0.
DO I=1-IBO-JBO,LINT
    GTI = GTI + BINT(0,0,I)*DELT(I)
end do
Z = -AR*GTI
N =  0
106 COLAT = 90. - FLAT
DO N=1,KMAX
    IF (N  >  1)  GO TO 115
    CL = COS(FLON*dfarg)
    SL = SIN(FLON*dfarg)
    CML(1) = CL
    SML(1) = SL
    GO TO 120
    115 SML(N) = SL*CML(N-1) + CL*SML(N-1)
    CML(N) = CL*CML(N-1) - SL*SML(N-1)
    120 CONTINUE
    DO M=0,N
        IF (IB  ==  2)  GO TO 121
        NMM = N - M
        IF ((NMM/2)*2  /=  NMM)  GO TO 150
        121 FFN = FN(N,M)
        CALL LEGFUN (M,FFN,CONST(N,M),COLAT,P,DP,PMS,0)
        IF (L  >=  0)  THEN
            AR = AOR**(FFN+2.)
        ELSE
            AR = 1.
            FFN = -2.
            DP = -DP
            PMS = -PMS
        END IF
        IF (M  /=  0)  GO TO 130
        BT1 = 0.
        BT3 = 0.
        BT  = 0.
        IF (L  ==  2)  GO TO 123
        IF (N  >  KINT)  GO TO 123
        GTI = 0.
        DO I=1-IBO-JBO,LINT
            GTI  = GTI  + BINT(N,M,I)*DELT(I)
        end do
        BT1  = AR*GTI
        BT3  = BT1
        123 IF (L  <=  1)  GO TO 125
        IF (N  >  KEXT)  GO TO 125
        GTE = 0.
        DO I=1-IBO-JBO,LEXT
            GTE = GTE + BEXT(N,M,I)*DELT(I)
        end do
        BT  = AOR3/AR*GTE
        BT1 = BT1 + BT
        125 X = X + BT1*DP
        Z = Z - (FFN*(BT3-BT)+BT3)*P
        GO TO 150
        130 BT1 = 0.
        BT2 = 0.
        BT3 = 0.
        BT  = 0.
        IF (L  ==  2)  GO TO 133
        IF (N  >  KINT)  GO TO 133
        GTI = 0.
        HTI = 0.
        DO I=1-IBO-JBO,LINT
            GTI = GTI + BINT(N,M,I)*DELT(I)
            HTI = HTI + BINT(M-1,N,I)*DELT(I)
        end do
        BT1 = AR*(GTI*CML(M) + HTI*SML(M))
        BT2 = AR*(GTI*SML(M) - HTI*CML(M))
        BT3 = BT1
        133 IF (L  <=  1)  GO TO 135
        IF (N  >  KEXT)  GO TO 135
        GTE = 0.
        HTE = 0.
        DO I=1-IBO-JBO,LEXT
            GTE = GTE + BEXT(N,M,I)*DELT(I)
            HTE = HTE + BEXT(M-1,N,I)*DELT(I)
        end do
        RA = AOR3/AR
        BT = RA*(GTE*CML(M) + HTE*SML(M))
        BT1 = BT1 + BT
        BT2 = BT2 + RA*(GTE*SML(M) - HTE*CML(M))
        135 X = X + BT1*DP
        Y = Y + BT2*PMS
        Z = Z - (FFN*(BT3-BT)+BT3)*P
    end do
end do
150 CONTINUE
BN = X
BE = Y
BV = Z
RETURN
999 STOP
END
!
!
SUBROUTINE TBFIT (T1,T2,IBF,THINT,TZERO)
    !-------------------------------------------------------------------
    !    COURTESY OF G.V. HAINES
    !
    !     T2    =  BEGINNING OF TIME INTERVAL.
    !     T1    =  END OF TIME INTERVAL.
    !     IBF   =  0   TO USE ORDINARY POLYNOMIALS AS TEMPORAL BASIS FUNCTIONS
    !              1          LEGENDRE POLYNOMIALS
    !              2          FOURIER SERIES
    !              3          COSINE SERIES
    !              4          SINE SERIES
    !              5          COSINE + SINE SERIES
    !     TZERO =  TIME-TRANSLATION PARAMETER:
    !              FOR IBF <= 1, CHOOSE TZERO = CENTER OF TIME INTERVAL
    !              FOR IBF >= 2, CHOOSE TZERO = BEGINNING OF TIME INTERVAL.
    !     THINT =  TIME-SCALING PARAMETER. THINT = HALF OF TIME INTERVAL T2-T1
    !              FOR IBF <= 1; PI*HALF OF TIME INTERVAL FOR IBF == 2;
    !              AND PI*TIME INTERVAL FOR IBF >= 3.
    !     NOTE:    CHOOSING TZERO AND THINT IN THIS WAY SCALES TIME
    !              TO (-1,1) FOR IBF <= 1;  TO (0,2PI) FOR IBF == 2;
    !              AND TO (0,PI) FOR IBF >= 3.
    !-------------------------------------------------------------------
    IBFI  = IBF
    IF (IBFI  <=  1)  THEN
        TZERO = (T2+T1)/2.D0
    ELSE
        TZERO =  T1
    ENDIF
    THINT = T2 - T1
    IF (IBFI  <=  2)  THINT = THINT/2.D0
    !      IF (IBFI  ==  4)  THEN
    !          DELT(0) = 0.D0
    !         ELSE
    !          DELT(0) = 1.D0
    !         ENDIF
    RETURN
END
!
!
SUBROUTINE LEGFUN (M,FN,CONST,COLAT,P,DP,PMS,IPRT)
    !-------------------------------------------------------------------
    !     SERIES FORM FOR ASSOCIATED LEGENDRE FUNCTION P, ITS DERIVATIVE DP,
    !     AND THE FUNCTION PMS=P*M/SIN(COLAT), IN POWERS OF (1-COS(COLAT))/2.
    !     INTEGRAL ORDER M, REAL DEGREE FN, NORMALIZING CONSTANT CONST.
    !     COLATITUDE COLAT IN DEGREES.
    !     IPRT = 0     NO PRINT-OUT
    !            1     PRINT PARAMETERS AND P SERIES
    !            2     PRINT PARAMETERS AND DP SERIES
    !            3     PRINT PARAMETERS AND BOTH P AND DP SERIES
    !           -1     PRINT PARAMETERS ONLY
    !     INPUT M,FN,CONST,COLAT,IPRT.   OUTPUT P,DP,PMS
    !    ADAPTED FROM G.V. HAINES (COMPUTERS   GEOSCIENCES, 14, 413-447, 1988).
    !-------------------------------------------------------------------
    REAL*8          FNN,AL,A,B,PNM,DPNM
    DIMENSION      AM(60), BM(60)
    LOGICAL        mess
    COMMON/iounit/konsol,mess
    DATA   JMAX/60/
    dfarg=(atan(1.0)*4.)/180.
    FNN = FN*(FN+1.)
    IF (COLAT  <  60.)  THEN
        X = SIN(dfarg*COLAT/2.)**2
        C = 1. - 2.*X
    ELSE
        C = COS(COLAT*dfarg)
        X = (1. - C)/2.
    END IF
    S = SIN(COLAT*dfarg)
    IF (M  > 1)  GO TO 20
    IF (M  <  0)  STOP
    AL = CONST
    GO TO 50
    20 AL = CONST*S**(M-1)
    50 PNM = AL
    DPNM = 0.
    J = 0
    100 J = J + 1
    JPM = J + M
    B = AL*((JPM-1)-FNN/JPM)
    DPNM = DPNM + B
    A = (B*X)/J
    PNM = PNM + A
    AL = A
    !     STORE P OR DP SERIES FOR PRINTOUT.
    IF (IPRT  <=  0)  GO TO 150
    IF (IPRT  ==  2)  GO TO 145
    AM(J) = A
    IF (IPRT  ==  1)  GO TO 150
    145 BM(J) = B
    !     CHECK FOR TERMINATION OF SERIES.
    150 ABSA = ABS(A)
    ABSB = ABS(B)
    IF (ABSB  >=  1.E-7)  GO TO 160
    IF (ABSA  <  1.E-7)  GO TO 110
    160 IF (ABSB  >=  1.E+13)  GO TO 105
    IF (ABSA  >=  1.E+13)  GO TO 105
    !     CHANGE CHECK LIMITS ACCORDING TO ACCURACY DESIRED AND ACCORDING
    !     TO WORD SIZE OF COMPUTER.
    !     FOR 32-BIT WORD, DOUBLE PRECISION, E-8 AND E+8 GIVE 7 DIGITS ACCURACY.
    !     FOR 60-BIT WORD, DOUBLE PRECISION, E-15 AND E+15 GIVE 14 DIGITS ACCURACY.
    !     FOR 60-BIT WORD, SINGLE PRECISION, E-8 AND E+7 GIVE 7 DIGITS ACCURACY.
    !     (DOUBLE OR SINGLE PRECISION REFER TO FNN,AL,A,B,PNM,DPNM)
    IF (J  <  JMAX)  GO TO 100
    !     CONVERGENCE SLOW OR JMAX TOO SMALL
    105 CONTINUE
    !     NUMERICAL ERROR UNACCEPTABLY LARGE DUE TO ADDING OF
    !     LARGE AND SMALL NUMBERS.
    if (mess) WRITE(konsol,108) M,FN,CONST,J,A,B
    108 FORMAT (//12H ** ERROR **/1X,I5,F10.5,E15.7,I5,2D15.7)
    STOP
    !     SERIES TRUNCATED SUCCESSFULLY.
    110 PS = PNM
    DPS = DPNM
    IF (M  /=  0)  GO TO 115
    PMS = 0.
    P = PS
    DP = DPS*S/2.
    GO TO 120
    115 PMS = PS*M
    P = PS*S
    DP = DPS*S*S/2. + C*PMS
    120 CONTINUE
    !     PRINT TERMS OF SERIES
    IF (IPRT  ==  0)  RETURN
    if (mess) WRITE(konsol,125) M,FN,CONST,COLAT,P,DP,PMS,J
    125 FORMAT (/1X,I5,F10.5,E20.12,F10.2,3F25.14,I5)
    IF (IPRT  <  0)  RETURN
    IF (IPRT  ==  2)  GO TO 135
    if (mess) write(konsol,130) (AM(I),I=1,J)
    130 FORMAT (1X,16E8.1)
    IF (IPRT  ==  1)  RETURN
    135 CONTINUE
    !  135 PRINT *, (BM(I),I=1,J)
    RETURN
END
!      
!      
REAL FUNCTION B0_98 ( HOUR, SAX, SUX, NSEASN, R, ZLO, ZMODIP)
    !-----------------------------------------------------------------
    ! Interpolation procedure for bottomside thickness parameter B0.
    ! Array B0F(ILT,ISEASON,IR,ILATI) distinguishes between day and
    ! night (ILT=1,2), four seasons (ISEASON is northern season with
    ! ISEASON=1 northern spring), low and high solar activity Rz12=10,
    ! 100 (IR=1,2), and modified dip latitudes of 0, 18 and 45
    ! degress (ILATI=1,2,3). In the DATA statement the first value
    ! corresponds to B0F(1,1,1,1), the second to B0F(2,1,1,1), the
    ! third to B0F(1,2,1,1) and so on.
    !
    ! input:
    !       hour    LT in decimal hours
    !       SAX     time of sunrise in decimal hours
    !       SUX     time of sunset in decimal hours
    !       nseasn  season in northern hemisphere (1=spring)
    !       R       12-month running mean of sunspot number
    !       ZLO     longitude
    !       ZMODIP  modified dip latitude
    !
    ! JUNE 1989 --------------------------------------- Dieter Bilitza
    !
    ! Updates (B0_new -> B0_98):
    !
    ! 01/98 corrected to include a smooth transition at the modip equator
    !       and no discontinuity at the equatorial change in season.
    ! 09/98 new B0 values incl values at the magnetic equator
    ! 10/98 longitude as input to determine if magnetic equator in northern 
    !         or southern hemisphere
    !-------------------------------------------------------------------
    REAL      NITVAL
    DIMENSION B0F(2,4,2,3),bfr(2,2,3),bfd(2,3),zx(5),g(6),dd(5)
    DATA      B0F/201,68,210,61,192,68,199,67,240,80,245,83,&
        233,71,230,65,108,65,142,81,110,68,77,75,&
        124,98,164,100,120,94,96,112,78,81,94,84,&
        81,81,65,70,102,87,127,91,109,88,81,78/
    DATA      zx/45.,72.,90.,108.,135./,dd/5*3.0/
    num_lat=3
    ! jseasn is southern hemisphere season
    jseasn=nseasn+2
    if(jseasn > 4) jseasn=jseasn-4
    zz = zmodip + 90.
    zz0 = 0.
    ! Interpolation in Rz12: linear from 10 to 100
    DO ISL=1,num_lat
        DO ISD=1,2
            bfr(isd,1,isl) = b0f(isd,nseasn,1,isl) +&
                (b0f(isd,nseasn,2,isl) - b0f(isd,nseasn,1,isl))/90.*(R-10.)
            bfr(isd,2,isl) = b0f(isd,jseasn,1,isl) +&
                (b0f(isd,jseasn,2,isl) - b0f(isd,jseasn,1,isl))/90.*(R-10.)
        end do
        ! Interpolation day/night with transitions at SAX (sunrise)
        ! and SUX (sunset) for northern/southern hemisphere iss=1/2
        do iss=1,2
            DAYVAL = BFR(1,ISS,ISL)
            NITVAL = BFR(2,ISS,ISL)
            BFD(iss,ISL) = HPOL(HOUR,DAYVAL,NITVAL,SAX,SUX,1.,1.)
        end do
    end do
    ! Interpolation with epstein-transitions in modified dip latitude.
    ! Transitions at +/-18 and +/-45 degrees; constant above +/-45.
    !
    ! g(1:5) are the latitudinal slopes of B0;
    !       g(1) is for the region from -90 to -45 degrees
    !       g(2) is for the region from -45 to -18 degrees
    !       g(3) is for the region from -18 to   0 degrees
    !       g(4) is for the region from   0 to  18 degrees
    !       g(5) is for the region from  18 to  45 degrees
    !       g(6) is for the region from  45 to  90 degrees
    !
    ! B0 =  bfd(2,3) at modip = -45,
    !       bfd(2,2) at modip = -18,
    !       bfd(2,1) or bfd(1,1) at modip = 0,
    !       bfd(1,2) at modip = 20,
    !       bfd(1,3) at modip = 45.
    ! If the Longitude is between 200 and 320 degrees then the modip 
    ! equator is in the southern hemisphere and bfd(2,1) is used at the 
    ! equator, otherwise bfd(1,1) is used.
    !
    zx1=bfd(2,3)
    zx2=bfd(2,2)
    zx3=bfd(1,1)
    if(zlo > 200.0.and.zlo < 320) zx3=bfd(2,1)
    zx4=bfd(1,2)
    zx5=bfd(1,3)
    g(1) = 0.
    g(2) = ( zx2 - zx1 ) / 27.
    g(3) = ( zx3 - zx2 ) / 18.
    g(4) = ( zx4 - zx3 ) / 18.
    g(5) = ( zx5 - zx4 ) / 27.
    g(6) = 0.
    !        bb0 = bfd(2,3)
    !      SUM = bb0
    sum=zx1
    DO I=1,5
        aa = eptr(zz ,dd(i),zx(i))
        bb = eptr(zz0,dd(i),zx(i))
        DSUM = (G(I+1) - G(I)) * (AA-BB) * dd(i)
        SUM = SUM + DSUM
    end do
    B0_98 = SUM
    RETURN
END
!
!
SUBROUTINE TAL(SHABR,SDELTA,SHBR,SDTDH0,AUS6,SPT)                         
    !-----------------------------------------------------------
    ! CALCULATES THE COEFFICIENTS SPT FOR THE POLYNOMIAL
    ! Y(X)=1+SPT(1)*X**2+SPT(2)*X**3+SPT(3)*X**4+SPT(4)*X**5               
    ! TO FIT THE VALLEY IN Y, REPRESENTED BY:                
    ! Y(X=0)=1, THE X VALUE OF THE DEEPEST VALLEY POINT (SHABR),                    
    ! THE PRECENTAGE DEPTH (SDELTA), THE WIDTH (SHBR) AND THE                       
    ! DERIVATIVE DY/DX AT THE UPPER VALLEY BOUNDRY (SDTDH0).                        
    ! IF THERE IS AN UNWANTED ADDITIONAL EXTREMUM IN THE VALLEY                     
    ! REGION, THEN AUS6=.TRUE., ELSE AUS6=.FALSE..     
    ! FOR -SDELTA THE COEFF. ARE CALCULATED FOR THE FUNCTION                        
    ! Y(X)=EXP(SPT(1)*X**2+...+SPT(4)*X**5).           
    !-----------------------------------------------------------
    DIMENSION SPT(4)                             
    LOGICAL AUS6    
    AUS6=.FALSE.
    if(SHBR <= 0.0) then
        AUS6=.TRUE.
        RETURN
    ENDIF
    Z1=-SDELTA/(100.0*SHABR*SHABR)               
    IF(SDELTA > 0.) GOTO 500                    
    SDELTA=-SDELTA  
    Z1=ALOG(1.-SDELTA/100.)/(SHABR*SHABR)        
    500   Z3=SDTDH0/(2.*SHBR)                          
    Z4=SHABR-SHBR   
    SPT(4)=2.0*(Z1*(SHBR-2.0*SHABR)*SHBR+Z3*Z4*SHABR)/                        &
        (SHABR*SHBR*Z4*Z4*Z4)                        
    SPT(3)=Z1*(2.0*SHBR-3.0*SHABR)/(SHABR*Z4*Z4)-&
        (2.*SHABR+SHBR)*SPT(4)          
    SPT(2)=-2.0*Z1/SHABR-2.0*SHABR*SPT(3)-3.0*SHABR*SHABR*SPT(4)              
    SPT(1)=Z1-SHABR*(SPT(2)+SHABR*(SPT(3)+SHABR*SPT(4)))                      
    B=4.*SPT(3)/(5.*SPT(4))+SHABR                
    C=-2.*SPT(1)/(5*SPT(4)*SHABR)                
    Z2=B*B/4.-C     
    IF(Z2 < 0.0) GOTO 300                       
    Z3=SQRT(Z2)     
    Z1=B/2.         
    Z2=-Z1+Z3       
    IF(Z2 > 0.0.AND.Z2 < SHBR) AUS6=.TRUE.     
    IF (ABS(Z3) > 1.E-15) GOTO 400              
    Z2=C/Z2         
    IF(Z2 > 0.0.AND.Z2 < SHBR) AUS6=.TRUE.     
    RETURN          
    400   Z2=-Z1-Z3       
    IF(Z2 > 0.0.AND.Z2 < SHBR) AUS6=.TRUE.     
    300   RETURN          
END             
!
!
SUBROUTINE VALGUL(XHI,HVB,VWU,VWA,VDP)
    ! --------------------------------------------------------------------- 
    !   CALCULATES E-F VALLEY PARAMETERS; T.L. GULYAEVA, ADVANCES IN
    !   SPACE RESEARCH 7, #6, 39-48, 1987.
    !
    !       INPUT:  XHI     SOLAR ZENITH ANGLE [DEGREE]
    !       
    !       OUTPUT: VDP     VALLEY DEPTH  (NVB/NME)
    !               VWU     VALLEY WIDTH  [KM]
    !               VWA     VALLEY WIDTH  (SMALLER, CORRECTED BY RAWER)
    !               HVB     HEIGHT OF VALLEY BASE [KM]
    ! -----------------------------------------------------------------------
    !
    COMMON  /CONST/UMR,PI
    !
    CS = 0.1 + COS(UMR*XHI)
    ABC = ABS(CS)
    VDP = 0.45 * CS / (0.1 + ABC ) + 0.55
    ARL = ( 0.1 + ABC + CS ) / ( 0.1 + ABC - CS)
    ZZZ = ALOG( ARL )
    VWU = 45. - 10. * ZZZ
    VWA = 45. -  5. * ZZZ
    HVB = 1000. / ( 7.024 + 0.224 * CS + 0.966 * ABC )
    RETURN
END
!
!
Subroutine DRegion(z,it,f,vKp,f5SW,f6WA,elg)
    !-----------------------------------------------------------------------
    ! Reference: Danilov, Rodevich, and Smirnova, Adv. Space Res.  
    !     15, #2, 165, 1995.
    !
    ! Input:     z    - solar zenith angle in degrees
    !            it   - season (month)
    !            f    - F10.7 solar radio flux (daily)
    !            vKp  - Kp magnetic index (3-hour)
    !            f5SW - indicator for Stratospheric Warming (SW) conditions
    !                   =0 no SW, =0.5 minor SW, =1 major SW
    !            f6WA - indicator for Winter Anomaly (WA) conditions
    !                   =0 no WA, =0.5 weak WA, =1 strong WA
    ! Criteria for SW and WA indicators:
    !      SW minor:  Temperature increase at the 30 hPa level by 10 deg.
    !      SA major:  The same but by 20 degrees.
    !         Temperature data for each year are published  
    !         in Beilage zur Berliner Wetterkarte (K. Labitzke et al.).
    !      WA weak:   An increase of the absorption in the 2-2.8 MHz  
    !                 range at short A3 paths by 15 dB
    !      WA strong: The same by 30 dB.
    ! 
    !       Only for month 12 to 2 (winter).
    !
    ! Output:      elg(7)  alog10 of electron density [cm-3] at h=60,65,
    !                  70,75,80,85, and 90km
    !-----------------------------------------------------------------------
    !            
    !or   dimension h(7),A0(7),A1(7),A2(7),A3(7),A4(7),A5(7),A6(7),elg(7)
    dimension A0(7),A1(7),A2(7),A3(7),A4(7),A5(7),A6(7),elg(7)
    data A0/1.0,1.2,1.4,1.5,1.6,1.7,3.0/
    data A1/0.6,0.8,1.1,1.2,1.3,1.4,1.0/
    data A2/0.,0.,0.08,0.12,0.05,0.2,0./
    data A3/0.,0.,0.,0.,0.,0.,1./
    data A4/0.,0.,-0.30,0.10,0.20,0.30,0.15/
    data A5/0.,-0.10,-0.20,-0.25,-0.30,-.30,0./
    data A6/0.,0.1,0.3,0.6,1.,1.,0.7/
    pi=3.14159265
    if(z <= 45) then
        f1z=1.
    else
        if(z < 90) then
            f1z=1.1892*(cos(z*pi/180))**0.5
        else
            f1z=0.
        endif
    endif
    f4S=1.
    if((it >= 5).and.(it <= 9))then
        f4S=0.
        f5SW=0
        f6WA=0
    endif
    if((it == 3).or.(it == 4).or.(it == 10).or.(it == 11))then
        f4S=0.5
        f5SW=0
        f6WA=0
    endif
    f2Kp=vKp
    if(vKp > 2) f2Kp=2.
    f3F=(f-60.)/300.*f1z
    do i=1,7
        elg(i)=A0(i)+A1(i)*f1z+A2(i)*f2Kp+A3(i)*f3F+A4(i)*f4S&
            +A5(i)*f5SW+A6(i)*f6WA
    end do
end
