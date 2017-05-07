!-----------------------------------------------------------------------
! Functions and subroutines for the International Reference Ionosphere 
! (IRI) model's electron density model. These functions and subroutines 
! are called by the main IRI subroutine IRI_SUB in IRISUB.FOR.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Required i/o units:  
!  KONSOL= 6 Program messages (used when jf(12)=.true. -> konsol)
!  KONSOL=11 Program messages (used when jf(12)=.false. -> MESSAGES.TXT)
!
!     COMMON/iounit/konsol,mess is used to pass the value of KONSOL from 
!     IRISUB to IRIFUN and IGRF. If mess=false then messages are turned off.
!     
!  UNIT=12 TCON: Solar/ionospheric indices IG12, R12 (IG_RZ.DAT) 
!  UNIT=13 APF,APFMSIS,APF_ONLY: Magnetic indices and F10.7 (APF107.DAT) 
!
! I/o Units used in other programs:
!  IUCCIR=10 in IRISUB for CCIR and URSI coefficients (CCIR%%.ASC, %%=month+10)
!  UNIT=14 in IGRF/GETSHC for IGRF coeff. (DGRF%%%%.DAT, %%%%=year)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                  
        
!*************************************************************   
!*************** ELECTRON DENSITY ****************************   
!*************************************************************   

!----------------------------------------------------------------
! DETERMINING ELECTRON DENSITY(M-3) IN THE TOPSIDE IONOSPHERE   
! (H=HMF2....2000 KM) BY HARMONIZED BENT-MODEL ADMITTING 
! VARIABILITY OF THE GLOBAL PARAMETERS BETA,ETA,DELTA,ZETA WITH        
! GEOM. LATITUDE, SMOOTHED SOLAR FLUX AND CRITICAL FREQUENCY.     
! BETA,ETA,DELTA,ZETA are computed in IRISUB program and 
! communicated via COMMON /BLO10. This is the IRI-2001 approach
! [REF.:K.RAWER,S.RAMAKRISHNAN,1978] 
! New options include:
! (1) IRI-corrected: TC3,alg10,hcor1 in COMMON /BLO11. 
!   TC3     correction term divided by (1500-(hcor1-hmF2))
!   alg10   = alog(10.)
!    hcor1    lower height boundary for correction
! (2) NeQuick:  B2TOP  in COMMON /BLO11.
!    B2TOP   is the topside scale height that depends on foF2 and 
!           hmF2. 
! Switch for choosing the desired option is itopn in COMMON /BLO11
!   itopn   =0 IRI-2001, =1 IRI-2001-corrected, =2 NeQuick
!           =3 Gulyaeva-0.5 is not yet implemented. 
!----------------------------------------------------------------
FUNCTION XE1(H)    

    COMMON  /BLOCK1/HMF2,XNMF2,HMF1,F1REG&
        /BLO10/BETA,ETA,DELTA,ZETA&
        /BLO11/B2TOP,TC3,itopn,alg10,hcor1&
        /QTOP/Y05,H05TOP,QF,XNETOP,xm3000,hhalf,tau&
        /ARGEXP/ARGMAX

    logical     f1reg              

    if(itopn == 2) then
        XE1 = TOPQ(H,XNMF2,HMF2,B2TOP)
        RETURN
    end if
    
    DXDH  = ( 1000. - HMF2 )/700.
    x0    = 300. - delta
    xmx0  = (H-HMF2)/DXDH
    x     = xmx0 + x0
    eptr1 = eptr(x, beta, 394.5) - eptr(x0, beta, 394.5)
    eptr2 = eptr(x, 100., 300.0) - eptr(x0, 100., 300.0) 
    y     = BETA * ETA * eptr1 + ZETA * (100. * eptr2 - xmx0)
    y     = y*dxdh

    if ( abs(y) > argmax ) y = sign(argmax, y)

    if ( itopn == 3 ) then
        if ( QF == 1. .and. ABS(H-H05TOP) < 1. ) QF = Y05 / Y
        XE1 = XNMF2 * EXP( -Y * QF)                             
        return          
    end if

    TCOR = 0.
    if ( itopn == 1 .and. h > hcor1 ) then
        xred = h - hcor1
        rco  = tc3 * xred
        tcor = rco * alg10
    end if

    XE1 = XNMF2 * EXP( -Y + TCOR )                             
    return          
end function XE1

!----------------------------------------------------------------
!  NeQuick formula
!----------------------------------------------------------------
real function TOPQ(h, No, hmax, Ho)
    REAL No
    PARAMETER (g=0.125, rfac=100.0)

    dh = h - hmax
    g1 = g * dh
    z  = rfac*Ho + g1
    z  = 1 + rfac * g1 / z
    z  = z*Ho
    z  = dh / z

    if (z > 40) then
        TOPQ = 0.0
        return
    end if

    ee = exp(z)
    if (ee > 1.0e7) then
        ep = 4.0/ee
    else
        ep = 4.0*ee/( 1.0 + ee )**2
    end if

    TOPQ = No*ep
    return
end function TOPQ


! FOR A PEAK AT X0 THE FUNCTION ZERO HAS TO BE EQUAL TO 0.
real function ZERO(DELTA)
    COMMON  /BLO10/         BETA,ETA,DEL,ZETA&
            /ARGEXP/        ARGMAX

    arg1 = delta/100.

    if ( abs(arg1) < argmax ) then
        z1 = 1./( 1. + exp(arg1) )
    else if (arg1 < 0) then
        z1 = 1.
    else
        z1 = 0.
    end if

    arg1 = ( delta + 94.5 )/beta
    if ( abs(arg1) < argmax ) then
        z2 = 1./( 1. + exp(arg1) )
    else if (arg1 < 0) then
        z2 = 1.
    else
        z2 = 0.
    end if

    ZERO = zeta*( 1. - z1 ) - eta*z2
    return
end function ZERO


! LOGARITHMIC DERIVATIVE OF FUNCTION XE1 (KM-1).   
function DXE1N(H)                            
    COMMON  /BLOCK1/HMF2,XNMF2,HMF1,F1REG&
            /BLO10/BETA,ETA,DELTA,ZETA                    
    logical f1reg

    x0 = 300. - delta
    X  = ( H - HMF2 )/( 1000.0 - HMF2 )*700.0 + x0
    epst2 = epst(x, 100.0, 300.0)
    epst1 = epst(x, beta , 394.5)
    DXE1N = -ETA * epst1 + ZETA * (1. - epst2)             
    RETURN          
end function DXE1N


! ELECTRON DENSITY FOR THE BOTTOMSIDE F-REGION (HMF1...HMF2).                   
real function XE2(H)                         
    COMMON  /BLOCK1/HMF2,XNMF2,HMF1,F1REG&
            /BLOCK2/B0,B1,C1    /ARGEXP/ARGMAX
    logical    f1reg
    X = ( HMF2 - H )/B0
    if (x <= 0.0) x = 0.0
    z = x**b1
    if (z > argmax) z = argmax
    XE2 = XNMF2 * exp( -z )/cosh( X )                 
    return          
end function XE2


! ELECTRON DENSITY FOR THE F1-LAYER (HZ.....HMF1)
! USING THE NEW DEFINED F1-LAYER FUNCTION (Reinisch and Huang, Advances 
! in Space Research, Volume 25, Number 1, 81-88, 2000)
real function XE3_1(H)
    COMMON    /BLOCK1/    HMF2,XNMF2,HMF1,F1REG&
              /BLOCK2/    B0,B1,D1F1
    logical    f1reg
    
    h1bar = h
    if ( f1reg ) h1bar = HmF1*( 1.0 - ( ( HmF1 - H )/HmF1 )**( 1.0 + D1F1 ) )
    XE3_1 = XE2( H1BAR )
    return
end function XE3_1


! ELECTRON DENSITY FOR THE INTERMEDIATE REGION (HEF...HZ)
! USING THE NEW DEFINED FUNCTION
real function XE4_1(H)
    COMMON    /BLOCK3/    HZ,T,HST&
              /BLOCK4/    HME,XNME,HEF
    
    if (hst < 0.0) then
        xe4_1 = xnme + t*( h - hef )
        return
    end if

    if ( HST == HEF ) then
        H1BAR = H
    else
        H1BAR = HZ + 0.5*T - SIGN(1.0, T) * SQRT( T*( 0.25*T + HZ - H ) )
    end if
    XE4_1 = XE3_1(H1BAR)
    return
end function XE4_1


! ELECTRON DENSITY FOR THE E AND VALLEY REGION (HME..HEF).   
real function XE5(H)                         
    LOGICAL NIGHT   
    COMMON  /BLOCK4/        HME,XNME,HEF&
            /BLOCK5/        NIGHT,E(4)                    

    T3 = H - HME        
    T1 = T3**2 * ( E(1) + T3*( E(2) + T3*( E(3) + T3*E(4) ) ) )

    if (NIGHT) then
        XE5 = XNME * exp(T1)                              
    else
        XE5 = XNME * ( 1 + T1 )  
    end if

    return          
end function XE5


! ELECTRON DENSITY FOR THE D REGION (HA...HME).    
real function XE6(H)                         
    COMMON  /BLOCK4/        HME,XNME,HEF&
            /BLOCK6/        HMD,XNMD,HDX&
            /BLOCK7/        D1,XKK,FP30,FP3U,FP1,FP2    

    IF (H > HDX) then
        z   = HME - H         
        XE6 = XNME * exp( -D1 * Z**XKK )
    else
        z   = H - HMD         
        FP3 = FP3U        
        IF ( z > 0.0 ) FP3 = FP30                        
        XE6 = XNMD * exp( Z*( FP1 + Z*( FP2 + Z*FP3 ) ) )
    end if
    return
end function XE6


! ELECTRON DENSITY BEETWEEN HA(KM) AND 1000 KM     
! SUMMARIZING PROCEDURES  NE1....6;                
real function XE_1(H)                          
    COMMON  /BLOCK1/HMF2,XNMF2,XHMF1,F1REG &
            /BLOCK3/HZ,T,HST&
            /BLOCK4/HME,XNME,HEF
    logical     f1reg

    if (f1reg) then
        hmf1 = xhmf1
    else
        hmf1 = hmf2
    end if

    if (H >= HMF2) then
        XE_1 = XE1(H)     
    else if(H >= HMF1) then
        XE_1 = XE2(H)       
    else if (H >= HZ) then
        XE_1 = XE3_1(H)       
    else if (H >= HEF) then
        XE_1 = XE4_1(H)       
    else if(H >= HME) then
        XE_1 = XE5(H)       
    else
        XE_1 = XE6(H)       
    end if
    return          
end function XE_1
