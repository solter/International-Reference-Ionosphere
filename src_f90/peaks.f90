!*************************************************************                  
!************* PEAK VALUES ELECTRON DENSITY ******************                  
!*************************************************************                  
!
!
real function FOUT(XMODIP,XLATI,XLONGI,UT,FF0)
    !--------------------------------------------------------------
    ! CALCULATES CRITICAL FREQUENCY FOF2/MHZ USING SUBROUTINE GAMMA1.      
    ! XMODIP = MODIFIED DIP LATITUDE, XLATI = GEOG. LATITUDE, XLONGI=
    ! LONGITUDE (ALL IN DEG.), MONTH = MONTH, UT =  UNIVERSAL TIME 
    ! (DEC. HOURS), FF0 = ARRAY WITH RZ12-ADJUSTED CCIR/URSI COEFF.
    ! D.BILITZA,JULY 85.
    !--------------------------------------------------------------
    DIMENSION FF0(988)
    INTEGER QF(9)
    DATA QF/11,11,8,4,1,0,0,0,0/
    FOUT=GAMMA1(XMODIP,XLATI,XLONGI,UT,6,QF,9,76,13,988,FF0)
    RETURN
END
!
!
real function XMOUT(XMODIP,XLATI,XLONGI,UT,XM0)
    !--------------------------------------------------------------
    ! CALCULATES PROPAGATION FACTOR M3000 USING THE SUBROUTINE GAMMA1.
    ! XMODIP = MODIFIED DIP LATITUDE, XLATI = GEOG. LATITUDE, XLONGI=
    ! LONGITUDE (ALL IN DEG.), MONTH = MONTH, UT =  UNIVERSAL TIME 
    ! (DEC. HOURS), XM0 = ARRAY WITH RZ12-ADJUSTED CCIR/URSI COEFF.
    ! D.BILITZA,JULY 85.
    !--------------------------------------------------------------
    DIMENSION XM0(441)
    INTEGER QM(7)
    DATA QM/6,7,5,2,1,0,0/
    XMOUT=GAMMA1(XMODIP,XLATI,XLONGI,UT,4,QM,7,49,9,441,XM0)
    RETURN
END
!
!
REAL FUNCTION HMF2ED(XMAGBR,R,X,XM3)         
    !--------------------------------------------------------------
    ! CALCULATES THE PEAK HEIGHT HMF2/KM FOR THE MAGNETIC                           
    ! LATITUDE XMAGBR/DEGREE AND THE SMOOTHED ZUERICH SUNSPOT                         
    ! NUMBER R USING CCIR-M3000 XM3 AND THE RATIO X=FOF2/FOE.
    ! FOLLOWING CCIR RECOMMENDATION X IS LIMITED TO VALUE
    ! GREATER OR EQUAL TO 1.7 .                       
    ! [REF. D.BILITZA ET AL., TELECOMM.J., 46, 549-553, 1979]                       
    ! D.BILITZA,1980.     
    !--------------------------------------------------------------
    F1=0.00232*R+0.222                         
    F2=1.2-0.0116*EXP(0.0239*R)            
    F3=0.096*(R-25.0)/150.0                      
    F4=1.0-R/150.0*EXP(-XMAGBR*XMAGBR/1600.0)
    if(x < 1.7) x=1.7
    DELM=F1*F4/(X-F2)+F3                
    HMF2ED=1490.0/(XM3+DELM)-176.0 
    RETURN          
END             
!
!
REAL FUNCTION XM3000HM(XMAGBR,R,X,HMF2)         
    !--------------------------------------------------------------
    ! CALCULATES THE PROPAGATION FACTOR M3000 FOR THE MAGNETIC LATITUDE
    ! XMAGBR/DEG. AND THE SMOOTHED ZUERICH SUNSPOT NUMBER R USING THE                        
    ! PEAK HEIGHT HMF2/KM AND THE RATIO X=FOF2/FOE. Reverse of HMF2ED.                      
    ! [REF. D.BILITZA ET AL., TELECOMM.J., 46, 549-553, 1979]                       
    ! D.BILITZA,1980. ----- no longer used    
    !--------------------------------------------------------------
    F1=0.00232*R+0.222                         
    F2=1.2-0.0116*EXP(0.0239*R)            
    F3=0.096*(R-25.0)/150.0                      
    F4=1.0-R/150.0*EXP(-XMAGBR*XMAGBR/1600.0)
    if(x < 1.7) x=1.7
    DELM=F1*F4/(X-F2)+F3                
    XM3000HM=1490.0/(HMF2+176.0)-DELM
    RETURN          
END             
!
!
SUBROUTINE  SHAMDHMF2 (RLAT,FLON,T,RZ,HMF2)
    !-------------------------------------------------------------------
    !    COMPUTES THE HOURLY VALUES OF hmF2 FROM A SET OF SH COEFFICIENTS
    !    IN A POINT OF A GIVEN GEOCENTRIC LATITUDE AND LONGITUDE
    !    OF THE EARTH'S SURFACE FOR A GIVEN MONTH AND A GIVEN SUSPOT NUMER.
    !   PARAMETERS AND COEFFICIENTS ARE GIVEN IN DATA STATEMENTS.
    !
    ! INPUT:    RLAT    The geogrphic latitude on the meridian given by 
    !                    the local time (FLON), where the modified dip
    !                   latitude is the same as of the orginal site.
    !            FLON    =LONGITUDE+15.*UT(hours)
    !            T        Month as a REAL number (1.0 to 12.0)
    !            RZ        12-month running mean
    ! OUTPUT    HMF2    F2 peak altitude in km
    !
    !  Altadill, D., S. Magdaleno, J.M. Torta, and E. Blanch 
    !     Adv. Space Res. 52, 1756-1769, 2013.
    !-------------------------------------------------------------------
    PARAMETER   (IBO=0,JBO=1,KDIM=8,LDIM=4,L=-1)
    DIMENSION   FN(0:KDIM,0:KDIM), CONST(0:KDIM,0:KDIM)
    DIMENSION   GNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        HNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
    DIMENSION     GANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        HANM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        GBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        HBNM(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
    DIMENSION   BINT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        BEXT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
    CHARACTER*1 IE       
    COMMON/AMTB/BINT,BEXT,RE,TZERO,IFIT,IB,KINT,LINT,KEXT,&
        LEXT,KMAX,FN
    !     ,CONST
    DATA THETA,RE,TZERO,IFIT,ICEN,IREF,IB,KINT,LINT,KEXT,LEXT&
        /180.,6371.2,1.0,  -1,  0,   0,   2, 8,   4,   0,   -1/
    DATA ((CONST(N,M), M=0,N), N=0,KDIM)&
        /4*1.,1.73205,0.866025,1.,2.44949,1.93649,0.7905691,1.,&
        3.16228,3.35410,2.09165,0.739510,1.,3.87298,5.12348,&
        4.18330,2.21853,0.701561,1.,4.58258,2*7.24569,4.96078,&
        2.32681,0.671693,1.,5.29150,9.72111,11.4564,9.49918,&
        5.69951,2.42182,0.647260,1.,6.,12.5499,16.9926,16.4531,&
        11.8645,6.40755,2.50683,0.626707/
    DATA IE /"I"/
    DATA (((GANM(N,M,J),GBNM(N,M,J),HANM(N,M,J),HBNM(N,M,J),&
        J=0,LDIM), M=0,N), N=0,KDIM)&
        / 340.552,-0.376,   0.000, 0.000,&
        -70.823, 1.150,   0.000, 0.000, -49.250, 0.591,   0.000, 0.000,&
        24.037,-0.401,   0.000, 0.000,  -8.780, 0.282,   0.000, 0.000,&
        -38.548, 0.177,   0.000, 0.000,&
        52.811,-0.282,   0.000, 0.000,  35.165,-0.265,   0.000, 0.000,&
        -21.227, 0.273,   0.000, 0.000,  10.183,-0.064,   0.000, 0.000,&
        67.946, 1.288,  67.397, 0.850,&
        -47.987,-1.643, -93.888,-1.257, -40.298,-1.075, -56.623,-0.768,&
        19.471, 0.532,  37.264, 0.400,  -0.830,-0.267, -17.603,-0.243,&
        -139.733,-1.204,   0.000, 0.000,&
        140.567, 1.145,   0.000, 0.000,  93.356, 0.980,   0.000, 0.000,&
        -48.353,-0.427,   0.000, 0.000,  21.370,-0.016,   0.000, 0.000,&
        -64.759, 0.228, -83.194, 0.121,&
        83.857,-0.228, 107.345, 0.013,  63.040,-0.268,  80.004,-0.254,&
        -48.692, 0.109, -52.767, 0.079,  15.922, 0.004,  14.904, 0.172,&
        -134.112, 3.405, 152.622,-2.016,&
        192.632,-4.277,-188.201, 2.383, 114.890,-2.568,-121.624, 1.712,&
        -69.827, 1.414,  76.826,-0.902,  35.669,-0.723, -33.598, 0.290,&
        -241.557,-0.964,   0.000, 0.000,&
        303.719, 1.393,   0.000, 0.000, 198.143, 0.731,   0.000, 0.000,&
        -136.985,-0.557,   0.000, 0.000,  56.365, 0.420,   0.000, 0.000,&
        -173.032, 0.275, -59.684,-0.368,&
        251.737,-0.285,  79.946, 0.667, 184.674,-0.313,  30.894, 0.828,&
        -100.746, 0.190, -14.771,-0.389,  30.175, 0.002,  21.618,-0.069,&
        6.871, 0.189,  70.508,-0.893,&
        -12.804, 0.018, -91.499, 1.128,  -5.759,-0.284, -60.616, 0.712,&
        21.336, 0.036,  37.937,-0.537,  -9.778, 0.175, -15.741, 0.211,&
        28.298, 1.499,  40.261,-2.376,&
        -54.316,-1.889, -45.615, 2.789, -34.102,-1.211, -36.189, 2.006,&
        20.337, 0.617,  15.090,-1.031, -10.873,-0.285,  -4.738, 0.324,&
        -90.022, 1.611,   0.000, 0.000,&
        139.241,-1.741,   0.000, 0.000,  69.750,-0.769,   0.000, 0.000,&
        -28.658, 0.280,   0.000, 0.000,  22.881,-0.282,   0.000, 0.000,&
        -44.702,-0.560,  87.636,-0.530,&
        53.150, 0.510,-119.895, 0.539,  31.695, 0.735, -87.077, 0.547,&
        6.348,-0.394,  54.557,-0.382,   0.875,-0.114, -15.031,-0.032,&
        -159.602, 0.277,-144.746, 0.246,&
        188.366,-0.118, 186.080,-0.124, 123.670,-0.172, 107.109,-0.067,&
        -58.392, 0.093, -64.485,-0.014,  30.929,-0.054,  37.817, 0.010,&
        19.949,-0.483,  29.250, 0.916,&
        -30.468, 0.725, -36.526,-1.320, -13.301, 0.259, -31.395,-0.638,&
        12.609,-0.244,  27.485, 0.374, -10.153, 0.253,  -5.509,-0.324,&
        42.493,-0.379, 138.667,-1.337,&
        -50.716, 0.303,-181.247, 1.744, -47.625, 0.434,-110.876, 1.038,&
        27.315,-0.192,  57.106,-0.545,   0.778,-0.106, -30.096, 0.284,&
        125.336, 0.325,   0.000, 0.000,&
        -164.302,-0.458,   0.000, 0.000,-119.173, 0.120,   0.000, 0.000,&
        78.273,-0.125,   0.000, 0.000, -26.471,-0.234,   0.000, 0.000,&
        -164.084, 1.528,  68.950,-2.145,&
        187.815,-1.754, -97.028, 2.521, 123.385,-1.348, -54.733, 1.333,&
        -58.316, 0.589,  33.701,-0.805,  27.475,-0.118, -22.520, 0.539,&
        30.638,-1.477, -71.140, 1.204,&
        -43.070, 1.448,  86.162,-1.443, -31.966, 1.202,  68.472,-1.014,&
        23.982,-0.552, -29.584, 0.633,  -3.900, 0.025,   4.363,-0.167,&
        81.222, 0.150,  66.889, 0.006,&
        -98.117,-0.095, -82.546, 0.149, -65.414,-0.057, -59.117,-0.113,&
        33.457, 0.113,  34.682,-0.021, -13.056,-0.019, -12.982, 0.178,&
        39.267,-0.337, -33.557, 0.109,&
        -42.755, 0.319,  34.263,-0.255, -26.894, 0.297,  30.837,-0.086,&
        6.342,-0.172, -17.113, 0.082,  -4.386,-0.001,   2.650,-0.097,&
        20*0.000,&
        -56.482,-0.951,   0.000, 0.000,&
        49.870, 1.103,   0.000, 0.000,  47.355, 0.689,   0.000, 0.000,&
        -21.830,-0.374,   0.000, 0.000,   0.708, 0.176,   0.000, 0.000,&
        8.535, 0.422, -56.552, 0.790,&
        -27.737,-0.208,  76.156,-0.834,  -9.129,-0.660,  44.140,-0.477,&
        -5.598, 0.234, -24.694, 0.257,  -2.797, 0.226,   9.278,-0.054,&
        -97.867, 1.578,  44.857,-0.495,&
        116.177,-1.680, -49.651, 0.408,  93.725,-1.524, -42.861, 0.447,&
        -55.150, 0.740,  28.612,-0.210,   7.615, 0.039,  -7.158,-0.040,&
        84.841,-1.243, -74.812, 0.103,&
        -100.335, 1.300,  87.349,-0.037, -78.486, 1.059,  70.399,-0.062,&
        50.717,-0.613, -36.938, 0.044, -12.384, 0.093,   7.046, 0.005,&
        50.573,-0.185,  47.967,-0.154,&
        -59.881, 0.390, -55.625, 0.206, -35.187, 0.074, -39.209, 0.028,&
        12.490,-0.039,  18.874,-0.072, -11.613, 0.191,  -8.454, 0.127,&
        40*0.000,&
        -243.265, 1.567,   0.000, 0.000,&
        306.183,-1.910,   0.000, 0.000, 183.262,-0.922,   0.000, 0.000,&
        -97.132, 0.553,   0.000, 0.000,  53.004,-0.403,   0.000, 0.000,&
        -212.916, 1.108,   8.322,-0.089,&
        257.690,-1.278, -15.853, 0.247, 179.077,-0.962, -15.457, 0.243,&
        -99.826, 0.628,   1.487,-0.098,  34.120,-0.178,   1.112,-0.006,&
        -61.739, 2.205,  62.259,-1.018,&
        69.646,-2.288, -77.640, 1.293,  46.838,-1.850, -65.274, 0.858,&
        -20.796, 0.919,  33.662,-0.542,   7.719,-0.165,  -5.840, 0.227,&
        57.782, 0.176,  30.904,-0.524,&
        -64.106,-0.272, -32.527, 0.468, -42.125,-0.123, -22.160, 0.471,&
        21.006, 0.034,  12.768,-0.214,  -9.326,-0.077,  -3.047,-0.047,&
        -15.644,-0.224, -96.643, 1.478,&
        30.526, 0.117, 116.896,-1.707,   8.932, 0.254,  72.474,-1.045,&
        -8.442,-0.098, -38.557, 0.559,  12.040,-0.100,  20.158,-0.290,&
        60*0.000,&
        -93.235, 1.670,   0.000, 0.000,&
        115.863,-1.904,   0.000, 0.000,  42.475,-0.878,   0.000, 0.000,&
        -16.805, 0.317,   0.000, 0.000,  26.766,-0.353,   0.000, 0.000,&
        -10.938,-0.552, 229.921,-3.050,&
        8.234, 0.386,-289.115, 3.548,   3.405, 0.746,-184.290, 2.115,&
        12.877,-0.410, 116.519,-1.229,  -3.016,-0.147, -49.942, 0.564,&
        -40.203,-0.273,   4.801, 0.004,&
        46.012, 0.180, -11.896, 0.129,  38.085, 0.142,   2.665,-0.054,&
        -6.701,-0.111,  -8.326, 0.065,   2.219,-0.065,  -3.252, 0.071,&
        -5.474, 0.833,  61.060,-0.178,&
        0.739,-0.777, -79.078, 0.293,  17.611,-0.764, -55.985, 0.102,&
        -7.107, 0.359,  40.779,-0.184,  -8.184, 0.030, -14.352, 0.127,&
        -51.003, 0.485,  36.183, 0.144,&
        65.647,-0.752, -42.617,-0.228,  33.590,-0.298, -27.554,-0.093,&
        -15.194, 0.193,  20.247, 0.033,  14.304,-0.227,  -6.789,-0.088,&
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
    !     SYNTHESIZES THE VALUE OF HMF2 FROM THE MODEL    
    !     **********************************************************
    CALL SCHNEVPDH (RZ,RLAT,FLON,dum,T,L,dum,dum,HMF2)
    RETURN
    9999  STOP
END
!
!
SUBROUTINE SCHNEVPDH (RZ,FLAT,FLON,R,T,L,BN,BE,BV)
    !------------------------------------------------------------------------
    ! WHEN L IS POSITIVE:
    ! COMPUTES SPHERICAL CAP HARMONIC (GEOCENTRIC) FIELD COMPONENTS
    ! HORIZONTAL NORTH BN,HORIZONTAL EAST BE,AND VERTICAL DOWNWARD BV.
    ! WHEN L IS NEGATIVE:
    ! COMPUTES GENERAL FUNCTION BV, ITS HORIZONTAL NORTH DERIVATIVE BN,
    ! AND ITS HORIZONTAL EAST DERIVATIVE BE, ON SPHERICAL CAP SURFACE.
    ! NOTE THAT THESE ARE METRICAL DERIVATIVES, AND BE IS THE
    ! LONGITUDINAL DERIVATIVE DIVIDED BY SIN(COLATITUDE).
    !      INPUT:
    ! FLAT,FLON,R ARE GEOCENTRIC SPHERICAL CAP LATITUDE,LONGITUDE,RADIAL
    ! DISTANCE; T IS TIME.
    ! L =  0  ON FIRST CALL:  RETURNS SPHERICAL CAP POLE POSITION FLATO,FLONO
    !         AND HALF-ANGLE THETA AS BN,BE, AND BV AFTER INITIALIZATION.
    !         ON SUBSEQUENT CALLS:  ACTS AS L=1.
    !      1  COMPUTES POTENTIAL FIELD COMPONENTS FROM INTERNAL COEFFICIENTS.
    !      2  COMPUTES POTENTIAL FIELD COMPONENTS FROM EXTERNAL COEFFICIENTS.
    !      3  COMPUTES FIELD FROM BOTH INTERNAL AND EXTERNAL COEFFICIENTS.
    !     -1  COMPUTES GENERAL FUNCTION BV AND DERIVATIVES BN WITH RESPECT TO
    !         LATITUDE AND BE WITH RESPECT TO LONGITUDE DIVIDED BY COS(LAT)
    !         (R IS DUMMY VARIABLE IN THIS CASE).
    ! NOTE:   SUBROUTINE IS INITIALIZED DURING FIRST CALL REGARDLESS OF L.
    !
    ! SUBROUTINE USED:  LEGFUN
    !
    ! PARAMS   COEFFS TRANSFERRED FROM MAIN PROGRAM IN COMMON/AMTB/
    !
    ! ADAPTED FROM SUBROUTINE SCHNEV OF G.V. HAINES (COMPUTERS   GEOSCIENCES, 
    ! 14, 413-447, 1988)
    !------------------------------------------------------------------------
    PARAMETER   (IBO=0,JBO=1,KDIM=8,LDIM=4)                           
    DIMENSION   FN(0:KDIM,0:KDIM), CONST(0:KDIM,0:KDIM)
    DIMENSION   CML(KDIM), SML(KDIM)
    DIMENSION   DELT(0:LDIM)
    DIMENSION   BINT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM),&
        BEXT(0:KDIM,0:KDIM,1-IBO-JBO:LDIM)
    COMMON/AMTB/BINT,BEXT,RE,TZERO,IFIT,IB,KINT,LINT,KEXT,&
        LEXT,KMAX,FN
    !     ,CONST
    CHARACTER*1 IE,RESP
    DATA ((CONST(N,M), M=0,N), N=0,KDIM)&
        /4*1.,1.73205,0.866025,1.,2.44949,1.93649,0.7905691,1.,&
        3.16228,3.35410,2.09165,0.739510,1.,3.87298,5.12348,&
        4.18330,2.21853,0.701561,1.,4.58258,2*7.24569,4.96078,&
        2.32681,0.671693,1.,5.29150,9.72111,11.4564,9.49918,&
        5.69951,2.42182,0.647260,1.,6.,12.5499,16.9926,16.4531,&
        11.8645,6.40755,2.50683,0.626707/
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
incept = 0                                                            set up
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
    CL = COS(FLON*UMR)
    SL = SIN(FLON*UMR)
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
subroutine SDMF2(UT,monthut,F107A,xmodip,long,hmF2)
    !--------------------------------------------------------------------------
    !    Global median model of the F2-layer peak height
    !
    !  Requires the following subroutines and functions:
    !     SDMF2, hmF2_med_SD, read_data_SD, fun_hmF2_SD,
    !     fun_Gk, Legendre, fun_hmF2UT, Koeff_UT, fun_Akp_UT,
    !     fun_Fk_UT, fun_Gk_UT  
    !
    ! Author of the code:
    !         Valentin Shubin
    !         Pushkov Institute of Terrestrial Magnetism,
    !         Ionosphere and Radio wave propagation (IZMIRAN)
    !         Moscow, Troitsk, 142190, Russia
    !         e-mail: shubin@izmiran.ru
    !         
    !     [Ref. V.N. Shubin. Global median model of the F2-layer
    !     peak height based on ionospheric radio-occultation and
    !     ground-based Digisonde observations. Advances in Space Research (2015)
    !     http://dx.doi.org/10.1016/j.asr.2015.05.029]
    !
    !  Input:
    !    UT      - universal time (real)
    !    monthut - month (integer)
    !    xmodip  - modified dip latitude in degrees (real)
    !    long    - geodatic longitude    in degrees (real)
    !    F107A   - F10.7 index averaged over the 3 Sun rotations 
    !              in units of 10^-22 W/(m^2 Hz) (real)
    !
    !  Output:
    !    hmF2  - F2-layer peak height in km (real)
    !--------------------------------------------------------------------------
    implicit none
    !     .. scalar arguments ..
    integer monthut
    real F107A
    real UT, hmF2
    real xmodip, long
    !     .. local scalars ..
    integer i
    double precision T
    !     .. local arrays ..
    double precision xUT(0:23)
    !    .. array in common ..
    double precision hmF2_UT(0:23)
    common/hmF2UT/hmF2_UT
    !     .. function references .
    real hmF2_med_SD, fun_hmF2UT
    !
    hmF2_UT = 0.0
    do i=0,23
        hmF2_UT(i) = hmF2_med_SD(i,monthut,F107A,xmodip,long)
        xUT(i) = dble(i)
    end do
    ! 
    T = dble(UT)
    hmF2 = fun_hmF2UT(T) 
    !
    return
end
!
!
real function hmF2_med_SD(iUT,monthut,F107A,xmodip,long)
    !---------------------------------------------------------------------
    !    Input: 
    !      iUT     - universal time (real)
    !      monthut - month (integer)
    !      F107A   - F10.7 index averaged over the 3 Sun rotations 
    !                in units of 10^-22 W/(m^2 Hz) (real)
    !      xmodip  - modified dip latitude in degrees (real)
    !      long    - geodatic longitude    in degrees (real)
    !
    !    function to interpolate hmF2 between the two levels of solar activity
    !    used the following auxiliary subroutines and functions:
    !    read_data_SD, fun_hmF2_SD
    !---------------------------------------------------------------------
    implicit none
    !    ..   scalar arguments ..
    integer monthut, iUT
    real F107A
    real xmodip, long
    !    real xmodip, long, umr, pi
    !    ..   local scalars ..
    real cov, cov1, cov2
    real a, b, hmF2_1, hmF2_2
    double precision teta
    !    ..   local arrays ..
    double precision coeff_month(0:148,0:47)
    double precision Kf(0:148)
    real ft1(12), ft2(12)
    !
    !    Arrays ft1 (12) and ft2 (12) are the median values of F10.7A,
    !    which were used as the margins for
    !    the low and high solar activity levels      
    !
    !               Jan   Feb   Mar   Apr   May  Jun
    data ft1/ 73.6, 72.3, 71.8, 70.9, 73.6,73.0,&
        71.1, 69.4, 69.1, 70.9, 72.3,74.6/
    data ft2/144.2,142.9,167.2,125.3,124.4,127.9,&
        142.0,165.9,132.6,142.0,145.6,143.0/
    !               Jul   Aug   Sep   Oct   Nov  Dec
    !    .. local in common ..
    double precision umr
    common/constt/umr
    !    common/const/umr,pi
    !     .. function references ..
    real fun_hmF2_SD
    !     .. subroutine references ..
    !     read_data_SD
    !
    umr=atan(1.0)*4./180
    teta = 90.0-xmodip
    !
    call read_data_SD(monthut,coeff_month)
    Kf = coeff_month(0:148,iUT)
    hmF2_1 = fun_hmF2_SD(teta,long,Kf)
    Kf = coeff_month(0:148,iUT+24)
    hmF2_2 = fun_hmF2_SD(teta,long,Kf)
    !
    cov = F107A
    cov1 = ft1(monthut) 
    cov2 = ft2(monthut)
    ! 
    a = (hmF2_2 - hmF2_1)/log(cov2/cov1)
    b =  hmF2_2 - a*log(cov2)
    hmF2_med_SD = a*log(cov) + b
    !
    return
end
!
!
subroutine read_data_SD(month,coeff_month)
    !------------------------------------------------------------------
    !    subroutine to read arrays mcsat11.datÖ mcsat22.dat
    !    with coefficients of hmF2 spatial decomposition
    !    for for 12 month, 24 UT hour and two solar activity levels
    !------------------------------------------------------------------
    implicit none
    !     .. scalar arguments ..
    integer month
    !     .. array arguments ..
    double precision coeff_month(0:148,0:47)
    !     .. local scalars ..
    integer coeff_month_read(1:12)
    character(256) filedata
    integer i, j
    !     .. local arrays ..
    double precision coeff_month_all(0:148,0:47,1:12)
    save coeff_month_all
    data coeff_month_read /12*0/
    !
    if (coeff_month_read(month)  ==  0) then
        write(filedata, 10) month+10
        open(10, File=filedata, status='old')
        do j=0,47
            read(10,20) (coeff_month_all(i,j,month),i=0,148)
        end do
        close(10)
        coeff_month_read(month) = 1
    end if
    !
    coeff_month = coeff_month_all(0:148,0:47,month)    
    !
    !    Previous operator can be replaced by
    !
    !      do i=0,148
    !       do j=0,47 
    !          coeff_month(i,j) = coeff_month_all(i,j,month)
    !        end do
    !    end do       
    !
    return
    
    ! 10   format('data\mcsat',i2,'.dat')
    10   format('mcsat',i2,'.dat')
    20   format(6(d12.5))
end
!
!      
real function fun_hmF2_SD(teta,long,Kf)
    !---------------------------------------------------------------------
    !    Input: 
    !      teta - (90-modip) in degrees
    !      long - geodatic longitude in degrees
    !      Kf   - coefficients of hmF2 spatial decomposition 
    !
    !    function to calculate spherical harmonics decomposition
    !    for the spatial dependencies of hmF2
    !    used the following auxiliary subroutines and functions:
    !    fun_Gk, Legendre
    !---------------------------------------------------------------------
    implicit none
    !    .. scalar arguments ..
    double precision teta
    real long
    !     .. array arguments ..
    double precision Kf(0:148)
    !    .. local scalars ..
    integer k
    double precision hmF2
    !     .. local arrays ..
    double precision Gk(0:148)
    !     .. subroutine references ..
    !     fun_Gk
    !
    call fun_Gk(teta,long,Gk)
    hmF2 = 0.d0
    do k=0,148
        hmF2 = hmF2 + Kf(k)*Gk(k) 
    end do
    fun_hmF2_SD = hmF2
    !
    return
end
!
!
subroutine fun_Gk(teta,long,Gk)
    !---------------------------------------------------------------------
    implicit none
    !    .. scalar arguments ..
    double precision teta
    real long
    !      real long, umr, pi
    !     .. array arguments ..
    double precision Gk(0:148)
    !    .. local scalars ..
    integer mm, nn, m, n, k
    !    .. local arrays ..
    double precision Pl_mn(0:8,0:12)
    !    .. local in common ..
    double precision umr
    common/constt/umr
    !    common/const/umr,pi
    !     .. subroutine references ..
    !     Legendre
    !
    Pl_mn = 0.d0
    mm = 8
    nn = 12
    call Legendre(mm,nn,Pl_mn,teta)
    Gk = 0.d0
    k = 0
    do m=0,mm
        if (m==0) then
            do n=0,nn
                Gk(k) = Pl_mn(m,n)
                k = k + 1
            end do
        else 
            do n=m,nn
                Gk(k)   = Pl_mn(m,n)*cos(m*long*umr)
                Gk(k+1) = Pl_mn(m,n)*sin(m*long*umr)
                k = k + 2
            end do
        end if
    end do
    !
    return
end
!
!
subroutine Legendre(mm,nn,p,teta)
    !---------------------------------------------------------------------
    !     Input:
    !       mm       - harmonics for longitude
    !       nn       - harmonics for the modified dip latitude (modip) 
    !       teta     - (90-modip) in degrees
    !    Output:
    !       P(mm,nn) - associated Legendre function    
    !
    !    subroutine to calculate associated Legendre function P(mm,nn)
    !    with Schmidt normalization
    !---------------------------------------------------------------------
    implicit none
    !     .. scalar arguments ..
    integer mm, nn
    !      real umr, pi
    double precision teta
    !     .. array arguments ..
    double precision p(0:mm,0:nn)
    !     .. local scalars ..
    integer j,l,m,n
    double precision z, s
    !    .. local in common ..
    double precision umr
    common/constt/umr
    !      common/const/umr,pi
    !
    p = 0.0
    z=cos(umr*teta)
    p(0,0)=1.
    p(0,1)=z
    if (mm /= 0) p(1,1)=sin(umr*teta)
    !
    do j=2,mm
        p(j,j)=(2*j-1)*p(j-1,j-1)*p(1,1)
    end do
    !
    do m=0,mm
        do n=1,nn
            if (m > n) then
                p(m,n) = 0.0
                cycle
            end if
            if ((n+1) > nn) exit
            if (n+1 == m) cycle
            if (m > (n-1)) then
                p(m,n+1)= (2*n+1)*z*p(m,n)/(n+1-m)
            else
                p(m,n+1)=((2*n+1)*z*p(m,n)-(n+m)*p(m,n-1))/(n+1-m)
            end if
        end do
    end do
    !
    do n=1,nn
        do m=1,mm
            if (m > n) then
                p(m,n) = 0.0
                exit
            end if
            s=1
            do l=n-m+1,n+m
                s=s*l
            end do
            p(m,n)=p(m,n)*sqrt(2./s)
        end do
    end do
    !
    return
end
!
!
real function fun_hmF2UT(T)
    !---------------------------------------------------------------------
    !    Input: T - universal time        
    ! 
    !    function to calculate Fourier  decomposition
    !    for the temporal variations of hmF2
    !    used the following auxiliary subroutines:
    !    Koeff_UT, fun_Akp_UT, fun_Fk_UT, fun_Gk_UT
    !---------------------------------------------------------------------
    implicit none
    !    .. scalar arguments ..
    double precision T
    !    .. local scalars ..
    integer k, mm, mk
    !      real dtr, dumr
    double precision hmF2
    !   .. local arrays ..    
    double precision Gk_UT(0:6), Kf_UT(0:6)
    !    .. local in common ..    
    double precision dtr
    common/radUT/dtr
    !      common/const1/dtr,dumr
    !   .. subroutine references ..
    !     Koeff_UT, fun_Gk_UT
    !
    dtr=atan(1.0)*4.0/12.0
    mm = 3
    mk = 2*mm
    !
    call Koeff_UT(mm,mk,Kf_UT)
    call fun_Gk_UT(mm,mk,t,Gk_UT)
    hmF2 = 0.d0
    do k=0,mk
        hmF2 = hmF2 + Kf_UT(k)*Gk_UT(k) 
    end do
    fun_hmF2UT = hmF2
    !
    return
end
!
!
subroutine Koeff_UT(mm,mk,Kf_UT)
    !---------------------------------------------------------------------
    implicit none
    !    ..  scalar arguments ..
    integer mm, mk
    !   .. array arguments ..
    double precision Kf_UT(0:mk)
    !    .. local scalars ..
    integer k, m
    double precision sum_D
    !   .. local arrays ..
    double precision Akp_UT(0:mk,0:mk)
    double precision Dk_UT(0:mk)
    !   .. subroutine references ..
    !   fun_Akp_UT 
    !
    call fun_Akp_UT(mm,mk,Akp_UT,Dk_UT)
    Kf_UT = 0.d0
    do k=mk,0,-1
        sum_D = 0.d0
        do m=k+1,mk
            sum_D = sum_D + Akp_UT(m,k)*Kf_UT(m)
        end do
        Kf_UT(k) = sum_D + Dk_UT(k)
    end do   
    return
end
!
!
subroutine fun_Akp_UT(mm,mk,Akp_UT,Dk_UT)
    !---------------------------------------------------------------------
    implicit none
    !    .. scalar arguments ..
    integer mm, mk
    !   .. array arguments ..
    double precision Akp_UT(0:mk,0:mk), Dk_UT(0:mk)
    !    .. local scalars ..
    integer i, k, p
    double precision t
    double precision sum_An, sum_Dn
    double precision sum_Ad, sum_Dd
    !   .. local arrays ..
    double precision Gk_UT(0:mk), Fk_UT(0:mk)
    !    .. array in common ..    
    double precision hmF2_UT(0:23)
    common/hmF2UT/hmF2_UT
    !   .. subroutine references ..
    !    fun_Gk_UT, fun_Fk_UT
    !
    Gk_UT = 0.d0
    Gk_UT(0) = 1.0
    Fk_UT = 0.d0
    Fk_UT(0) = 1.d0
    Akp_UT = 0.d0
    Dk_UT = 0.d0
    do p=0,mk
        sum_Dn=0.d0
        sum_Dd=0.d0
        do k=p+1,mk
            sum_An=0.d0
            sum_Ad=0.d0
            do i=0,23
                t = dble(i)
                call fun_Gk_UT(mm,mk,t,Gk_UT)
                call fun_Fk_UT(mk,Gk_UT,Akp_UT,Fk_UT)
                sum_An = sum_An + Gk_UT(k)*Fk_UT(p)
                sum_Ad = sum_Ad + Fk_UT(p)*Fk_UT(p)
                if (p == (k-1)) then
                    sum_Dn = sum_Dn + hmF2_UT(i)*Fk_UT(p)
                    sum_Dd = sum_Dd + Fk_UT(p)*Fk_UT(p)
                end if
            end do
            Akp_UT(k,p) = - sum_An/sum_Ad
        end do
        if (p < mk) then
            Dk_UT(p) = sum_Dn/sum_Dd
        end if
    end do
    !
    p=mk
    sum_Dn=0.d0
    sum_Dd=0.d0
    do i=0,23
        t = dble(i)
        call fun_Gk_UT(mm,mk,t,Gk_UT)
        call fun_Fk_UT(mk,Gk_UT,Akp_UT,Fk_UT)
        sum_Dn = sum_Dn + hmF2_UT(i)*Fk_UT(p)
        sum_Dd = sum_Dd + Fk_UT(p)*Fk_UT(p)
    end do
    Dk_UT(p) = sum_Dn/sum_Dd
    !
    return
end
!
!
subroutine fun_Fk_UT(mk,Gk_UT,Akp_UT,Fk_UT)
    !---------------------------------------------------------------------
    implicit none
    !    .. scalar arguments ..
    integer mk
    !   .. array arguments ..
    double precision Gk_UT(0:mk)
    double precision Akp_UT(0:mk,0:mk)
    double precision Fk_UT(0:mk)
    !    .. local scalars ..
    integer k, p
    double precision sum_G
    !
    Fk_UT = 0.d0
    do k=0,mk
        sum_G = 0.d0
        do p=0,k
            if (k == p) cycle
            sum_G = sum_G + Akp_UT(k,p)*Fk_UT(p)
        end do
        Fk_UT(k) = sum_G + Gk_UT(k)
    end do
    !
    return
end
!
!
subroutine fun_Gk_UT(mm,mk,t,Gk_UT)
    !---------------------------------------------------------------------
    implicit none
    !    .. scalar arguments ..
    integer mm, mk
    !      real dtr, dumr
    double precision t
    !   .. array arguments ..
    double precision Gk_UT(0:mk)
    !    .. local scalars ..
    integer m, k
    !    .. local in common ..
    double precision dtr
    common/radUT/dtr
    !      common/const1/dtr, dumr
    !
    Gk_UT = 0.d0
    k = 0
    do m=0,mm
        if (m == 0) then
            Gk_UT(k) = 1
            k = k + 1
        else 
            Gk_UT(k)   = cos(m*t*dtr)
            Gk_UT(k+1) = sin(m*t*dtr)
            k = k + 2
        end if
    end do
    !
    return
end
!
!
REAL FUNCTION FOF1ED(YLATI,R,CHI)
    !--------------------------------------------------------------
    ! CALCULATES THE F1 PEAK PLASMA FREQUENCY (FOF1/MHZ)
    ! INPUT:   
    !        YLATI    ABSOLUT VALUE OF DIP-LATITUDE IN DEGREE
    !       R        12-MONTH RUNNING MEAN OF SUNSPOT NUMBER 
    !       CHI        SOLAR ZENITH ANGLE IN DEGREE
    ! REFERENCE: 
    !       E.D.DUCHARME ET AL., RADIO SCIENCE 6, 369-378, 1971
    !                                      AND 8, 837-839, 1973
    !       HOWEVER WITH MAGNETIC DIP LATITUDE INSTEAD OF GEOMAGNETIC
    !       DIPOLE LATITUDE, EYFRIG, 1979                    
    !--------------------------------------------- D. BILITZA, 1988.   
    COMMON/CONST/UMR,PI
    fof1ed=0.0
    if (chi > 90.0) return
    DLA =  YLATI
    F0 = 4.35 + DLA * ( 0.0058 - 1.2E-4 * DLA ) 
    F100 = 5.348 + DLA * ( 0.011 - 2.3E-4 * DLA )
    FS = F0 + ( F100 - F0 ) * R / 100.0
    XMUE = 0.093 + DLA * ( 0.0046 - 5.4E-5 * DLA ) + 3.0E-4 * R
    FOF1 = FS * COS( CHI * UMR ) ** XMUE
    CHI0 = 49.84733 + 0.349504 * DLA
    CHI100 = 38.96113 + 0.509932 * DLA
    CHIM = ( CHI0 + ( CHI100 - CHI0 ) * R / 100. )
    IF(CHI > CHIM) FOF1=-FOF1 
    FOF1ED = FOF1     
    RETURN
END             
!
!
real function f1_c1(xmodip,hour,suxnon,saxnon)
    ! F1 layer shape parameter C1 after Reinisch and Huang, Advances in
    ! Space Research, Volume 25, Number 1, 81-88, 2000.
    common    /const/umr,pi
    pi = umr * 180.
    
    ABSMDP=ABS(XMODIP)
    DELA=4.32
    IF(ABSMDP >= 18.) DELA=1.0+EXP(-(ABSMDP-30.0)/10.0)
    C1OLD = 0.09 + 0.11/DELA
    if(suxnon == saxnon) then
        c1 = 2.5 * c1old
    else
        c1 = 2.5*c1old*cos((HOUR-12.)/(suxnon-saxnon)*pi)
    endif
    if(c1 < 0.0) c1=0.0
    f1_c1=c1
    return
end
!
!
subroutine f1_prob (sza,glat,rz12,f1prob,f1probl)
    !--------------------------------------------------------------------------
    ! Occurrence probability of F1 layer after Scotto et al., Advances in
    ! Space Research, Volume 20, Number 9, 1773-1775, 1997.
    !
    ! Input:     sza        solar zenith angle in degrees 
    !             glat    geomagnetic latitude in degrees
    !            rz12    12-month running mean of sunspot number
    ! Output:     f1prob    F1 occurrence probability without L-condition cases 
    !             f1probl    F1 occurrence probability with L-condition cases
    !--------------------------------------------------------------------------
    !
    common /const/umr,pi
    xarg = 0.5 + 0.5 * cos(sza*umr)
    a = 2.98 + 0.0854 * rz12
    b = 0.0107 - 0.0022 * rz12
    c = -0.000256 + 0.0000147 * rz12
    gamma = a + ( b + c * glat) * glat
    f1pr = xarg ** gamma
    if(f1pr < 1.e-3) f1pr=0.0
    f1prob=f1pr
    f1prl = xarg ** 2.36
    if(f1prl < 1.e-3) f1prl=0.0
    f1probl=f1prl
    return
end
!
!
REAL FUNCTION FOEEDI(COV,XHI,XHIM,XLATI)
    !-------------------------------------------------------
    ! CALCULATES FOE/MHZ BY THE EDINBURGH-METHOD.      
    ! INPUT: 
    !     COV        MONTHLY MEAN 10.7CM SOLAR RADIO FLUX measured at 
    !           ground level  
    !   XHI        SOLAR ZENITH ANGLE IN DEGREE 
    !   XHIM    SOLAR ZENITH ANGLE AT NOON IN DEGREE
    !   XLATI     ABSOLUTE VALUE OF GEOGRAPHIC LATITUDE IN DEGREE, 
    ! REFERENCE: 
    !       KOURIS-MUGGELETON, CCIR DOC. 6/3/07, 1973
    !       TROST, J. GEOPHYS. RES. 84, 2736, 1979 (was used
    !               to improve the nighttime varition)
    !       RAWER AND BILITZA, Adv. Space Res. 10(8), 5-14, 1990
    ! D.BILITZA--------------------------------- AUGUST 1986.    
    COMMON/CONST/UMR,PI
    ! variation with solar activity (factor A) ...............
    A=1.0+0.0094*(COV-66.0)                      
    ! variation with noon solar zenith angle (B) and with latitude (C)
    SL=COS(XLATI*UMR)
    IF(XLATI < 32.0) THEN
        SM=-1.93+1.92*SL                             
        C=23.0+116.0*SL                              
    ELSE
        SM=0.11-0.49*SL                              
        C=92.0+35.0*SL  
    ENDIF
    if(XHIM >= 90.) XHIM=89.999
    B = COS(XHIM*UMR) ** SM
    ! variation with solar zenith angle (D) ..........................        
    IF(XLATI > 12.0) THEN
        SP=1.2
    ELSE
        SP=1.31         
    ENDIF
    ! adjusted solar zenith angle during nighttime (XHIC) .............
    XHIC=XHI-3.*ALOG(1.+EXP((XHI-89.98)/3.))   
    D=COS(XHIC*UMR)**SP       
    ! determine foE**4 ................................................
    R4FOE=A*B*C*D     
    ! minimum allowable foE (foe_min=sqrt[SMIN])...............................
    SMIN=0.121+0.0015*(COV-60.)
    SMIN=SMIN*SMIN
    IF(R4FOE < SMIN) R4FOE=SMIN                     
    FOEEDI=R4FOE**0.25                           
    RETURN          
END   
!
!
REAL FUNCTION XMDED(XHI,R,YW)                
    ! D. BILITZA, 1978, CALCULATES ELECTRON DENSITY OF D MAXIMUM.                   
    ! XHI/DEG. IS SOLAR ZENITH ANGLE, R SMOOTHED ZURICH SUNSPOT NUMBER              
    ! AND YW/M-3 THE ASSUMED CONSTANT NIGHT VALUE.     
    ! [REF.: D.BILITZA, WORLD DATA CENTER A REPORT UAG-82,7,BOULDER,1981]
    ! corrected 4/25/97 - D. Bilitza
    !
    COMMON/CONST/UMR,PI
    !
    if(xhi >= 90) goto 100
    Y = 6.05E8 + 0.088E8 * R
    yy = cos ( xhi * umr )
    yyy = -0.1 / ( yy**2.7 ) 
    if (yyy < -40.) then 
        ymd=0.0
    else
        ymd = y * exp(yyy)
    endif
    if (ymd < yw) ymd = yw
    xmded=ymd
    RETURN          
    100     XMDED=YW        
    RETURN          
END
!
!
REAL FUNCTION GAMMA1(SMODIP,SLAT,SLONG,HOUR,&
        IHARM,NQ,K1,M,MM,M3,SFE)      
    !---------------------------------------------------------------
    ! CALCULATES GAMMA1=FOF2 OR M3000 USING CCIR NUMERICAL MAP                      
    ! COEFFICIENTS SFE(M3) FOR MODIFIED DIP LATITUDE (SMODIP/DEG)
    ! GEOGRAPHIC LATITUDE (SLAT/DEG) AND LONGITUDE (SLONG/DEG)  
    ! AND UNIVERSIAL TIME (HOUR/DECIMAL HOURS). IHARM IS THE MAXIMUM
    ! NUMBER OF HARMONICS USED FOR DESCRIBING DIURNAL VARIATION.
    ! NQ(K1) IS AN INTEGER ARRAY GIVING THE HIGHEST DEGREES IN 
    ! LATITUDE FOR EACH LONGITUDE HARMONIC WHERE K1 GIVES THE NUMBER 
    ! OF LONGITUDE HARMONICS. M IS THE NUMBER OF COEFFICIENTS FOR 
    ! DESCRIBING VARIATIONS WITH SMODIP, SLAT, AND SLONG. MM IS THE
    ! NUMBER OF COEFFICIENTS FOR THE FOURIER TIME SERIES DESCRIBING
    ! VARIATIONS WITH UT.
    ! M=1+NQ(1)+2*[NQ(2)+1]+2*[NQ(3)+1]+... , MM=2*IHARM+1, M3=M*MM  
    ! SHEIKH,4.3.77.      
    !---------------------------------------------------------------
    REAL*8 C(12),S(12),COEF(100),SUM             
    DIMENSION NQ(K1),XSINX(13),SFE(M3)           
    COMMON/CONST/UMR,PI
    HOU=(15.0*HOUR-180.0)*UMR                    
    S(1)=SIN(HOU)   
    C(1)=COS(HOU)   
    DO I=2,IHARM                             
        C(I)=C(1)*C(I-1)-S(1)*S(I-1)                 
        S(I)=C(1)*S(I-1)+S(1)*C(I-1)                 
    end do
    DO I=1,M    
        MI=(I-1)*MM     
        COEF(I)=SFE(MI+1)                            
        DO J=1,IHARM                             
            COEF(I)=COEF(I)+SFE(MI+2*J)*S(J)+SFE(MI+2*J+1)*C(J)                       
        end do
    end do
    SUM=COEF(1)     
    SS=SIN(SMODIP*UMR)                           
    S3=SS           
    XSINX(1)=1.0    
    INDEX=NQ(1)     
    DO J=1,INDEX                             
        SUM=SUM+COEF(1+J)*SS                         
        XSINX(J+1)=SS   
        SS=SS*S3        
    end do
    XSINX(NQ(1)+2)=SS                            
    NP=NQ(1)+1      
    SS=COS(SLAT*UMR)                             
    S3=SS           
    DO J=2,K1   
        S0=SLONG*(J-1.)*UMR                          
        S1=COS(S0)      
        S2=SIN(S0)      
        INDEX=NQ(J)+1   
        DO L=1,INDEX                             
            NP=NP+1         
            SUM=SUM+COEF(NP)*XSINX(L)*SS*S1              
            NP=NP+1         
            SUM=SUM+COEF(NP)*XSINX(L)*SS*S2              
        end do
        SS=SS*S3        
    end do
    
    GAMMA1=SUM      
    RETURN          
END 
!

