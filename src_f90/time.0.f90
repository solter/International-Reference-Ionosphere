!******************************************************************
!********** ZENITH ANGLE, DAY OF YEAR, TIME ***********************
!******************************************************************

subroutine soco (ld,t,flat,Elon,height,&
        DECLIN, ZENITH, SUNRSE, SUNSET)
    !--------------------------------------------------------------------
    !       s/r to calculate the solar declination, zenith angle, and
    !       sunrise   sunset times  - based on Newbern Smith's algorithm
    !       [leo mcnamara, 1-sep-86, last modified 16-jun-87]
    !       {dieter bilitza, 30-oct-89, modified for IRI application}
    !
    ! in:   ld      local day of year
    !       t       local hour (decimal)
    !       flat    northern latitude in degrees
    !       elon    east longitude in degrees
    !        height    height in km
    !
    ! out:  declin      declination of the sun in degrees
    !       zenith      zenith angle of the sun in degrees
    !       sunrse      local time of sunrise in hours 
    !       sunset      local time of sunset in hours 
    !-------------------------------------------------------------------
    !
    common/const/dtr,pi     /const1/humr,dumr
    ! amplitudes of Fourier coefficients  --  1955 epoch.................
    data    p1,p2,p3,p4,p6 /&
        0.017203534,0.034407068,0.051610602,0.068814136,0.103221204 /
    !
    ! s/r is formulated in terms of WEST longitude.......................
    wlon = 360. - Elon
    !
    ! time of equinox for 1980...........................................
    td = ld + (t + Wlon/15.) / 24.
    te = td + 0.9369
    !
    ! declination of the sun..............................................
    dcl = 23.256 * sin(p1*(te-82.242)) + 0.381 * sin(p2*(te-44.855))&
        + 0.167 * sin(p3*(te-23.355)) - 0.013 * sin(p4*(te+11.97))&
        + 0.011 * sin(p6*(te-10.41)) + 0.339137
    DECLIN = dcl
    dc = dcl * dtr
    !
    ! the equation of time................................................
    tf = te - 0.5
    eqt = -7.38*sin(p1*(tf-4.)) - 9.87*sin(p2*(tf+9.))&
        + 0.27*sin(p3*(tf-53.)) - 0.2*cos(p4*(tf-17.))
    et = eqt * dtr / 4.
    !
    fa = flat * dtr
    phi = humr * ( t - 12.) + et
    !
    a = sin(fa) * sin(dc)
    b = cos(fa) * cos(dc)
    cosx = a + b * cos(phi)
    if(abs(cosx) > 1.) cosx=sign(1.,cosx)
    zenith = acos(cosx) / dtr
    !
    ! calculate sunrise and sunset times --  at the ground...........
    ! see Explanatory Supplement to the Ephemeris (1961) pg 401......
    ! sunrise at height h metres is at...............................
    h=height*1000.
    chih = 90.83 + 0.0347 * sqrt(h)
    ! this includes corrections for horizontal refraction and........
    ! semi-diameter of the solar disk................................
    ch = cos(chih * dtr)
    cosphi = (ch -a ) / b
    ! if abs(secphi) > 1., sun does not rise/set.....................
    ! allow for sun never setting - high latitude summer.............
    secphi = 999999.
    if(cosphi /= 0.) secphi = 1./cosphi
    sunset = 99.
    sunrse = 99.
    if(secphi > -1.0.and.secphi <= 0.) return
    ! allow for sun never rising - high latitude winter..............
    sunset = -99.
    sunrse = -99.
    if(secphi > 0.0.and.secphi < 1.) return
    !
    cosx = cosphi
    if(abs(cosx) > 1.) cosx=sign(1.,cosx)
    phi = acos(cosx)
    et = et / humr
    phi = phi / humr
    sunrse = 12. - phi - et
    sunset = 12. + phi - et
    if(sunrse < 0.) sunrse = sunrse + 24.
    if(sunset >= 24.) sunset = sunset - 24.
    !
    return
end
!
!
FUNCTION HPOL(HOUR,TW,XNW,SA,SU,DSA,DSU)            
    !-------------------------------------------------------
    ! PROCEDURE FOR SMOOTH TIME-INTERPOLATION USING EPSTEIN  
    ! STEP FUNCTION AT SUNRISE (SA) AND SUNSET (SU). THE 
    ! STEP-WIDTH FOR SUNRISE IS DSA AND FOR SUNSET DSU.
    ! TW,NW ARE THE DAY AND NIGHT VALUE OF THE PARAMETER TO 
    ! BE INTERPOLATED. SA AND SU ARE TIME OF SUNRIES AND 
    ! SUNSET IN DECIMAL HOURS.
    ! BILITZA----------------------------------------- 1979.
    IF(ABS(SU) > 25.) THEN
        IF(SU > 0.0) THEN
            HPOL=TW
        ELSE
            HPOL=XNW
        ENDIF
        RETURN
    ENDIF
    HPOL=XNW+(TW-XNW)*EPST(HOUR,DSA,SA)+&
        (XNW-TW)*EPST(HOUR,DSU,SU) 
    RETURN          
END       
!      
!
SUBROUTINE MODA(IN,IYEAR,MONTH,IDAY,IDOY,NRDAYMO)
    !-------------------------------------------------------------------
    ! CALCULATES DAY OF YEAR (IDOY, ddd) FROM YEAR (IYEAR, yy or yyyy), 
    ! MONTH (MONTH, mm) AND DAY OF MONTH (IDAY, dd) IF IN=0, OR MONTH 
    ! AND DAY FROM YEAR AND DAY OF YEAR IF IN=1. NRDAYMO is an output 
    ! parameter providing the number of days in the specific month.
    !-------------------------------------------------------------------
    DIMENSION       MM(12)
    DATA            MM/31,28,31,30,31,30,31,31,30,31,30,31/
    IMO=0
    MOBE=0
    !
    !  leap year rule: years evenly divisible by 4 are leap years, except
    !  years also evenly divisible by 100 are not leap years, except years 
    !  also evenly divisible by 400 are leap years. The year 2000 therefore 
    !  is a leap year. The 100 and 400 year exception rule
    !     if((iyear/4*4 == iyear).and.(iyear/100*100 /= iyear)) mm(2)=29
    !  will become important again in the year 2100 which is not a leap 
    !  year.
    !
    mm(2)=28
    if(iyear/4*4 == iyear) mm(2)=29
    IF(IN > 0) GOTO 5
    mosum=0
    if(month > 1) then
        do i=1,month-1 
            mosum=mosum+mm(i)
        end do
    endif
    idoy=mosum+iday
    nrdaymo=mm(month)
    RETURN
    5       IMO=IMO+1
    IF(IMO > 12) GOTO 55
    MOOLD=MOBE
    nrdaymo=mm(imo)
    MOBE=MOBE+nrdaymo
    IF(MOBE < IDOY) GOTO 5
    55              MONTH=IMO
    IDAY=IDOY-MOOLD
    RETURN
END             
!
!
subroutine ut_lt(mode,ut,slt,glong,iyyy,ddd)
    ! -----------------------------------------------------------------
    ! Converts Universal Time UT (decimal hours) into Solar Local Time
    ! SLT (decimal hours) for given date (iyyy is year, e.g. 1995; ddd
    ! is day of year, e.g. 1 for Jan 1) and geodatic longitude in degrees.
    ! For mode=0 UT->LT and for mode=1 LT->UT
    ! Please NOTE that iyyy and ddd are input as well as output parameters
    ! since the determined LT may be for a day before or after the UT day.
    ! ------------------------------------------------- bilitza nov 95
    integer         ddd,dddend
    xlong=glong
    if(glong > 180) xlong=glong-360
    if(mode /= 0) goto 1
    !
    ! UT ---> LT
    !
    SLT=UT+xlong/15.
    if((SLT >= 0.).and.(SLT <= 24.)) goto 2
    if(SLT > 24.) goto 3
    SLT=SLT+24.
    ddd=ddd-1
    if(ddd < 1.) then
        iyyy=iyyy-1
        ddd=365
        !
        ! leap year if evenly divisible by 4 and not by 100, except if evenly
        ! divisible by 400. Thus 2000 will be a leap year.
        !
        if(iyyy/4*4 == iyyy) ddd=366
    endif
    goto 2
    3               SLT=SLT-24.
    ddd=ddd+1
    dddend=365
    if(iyyy/4*4 == iyyy) dddend=366
    if(ddd > dddend) then
        iyyy=iyyy+1
        ddd=1
    endif
    goto 2
    !
    ! LT ---> UT
    !
    1       UT=SLT-xlong/15.
    if((UT >= 0.).and.(UT <= 24.)) goto 2
    if(UT > 24.) goto 5
    UT=UT+24.
    ddd=ddd-1
    if(ddd < 1.) then
        iyyy=iyyy-1
        ddd=365
        if(iyyy/4*4 == iyyy) ddd=366
    endif
    goto 2
    5               UT=UT-24.
    ddd=ddd+1
    dddend=365
    if(iyyy/4*4 == iyyy) dddend=366
    if(ddd > dddend) then
        iyyy=iyyy+1
        ddd=1
    endif
    2       return
end
!
!
SUBROUTINE SUN (IYEAR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
    !-----------------------------------------------------------------------------
    !  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
    !  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
    !
    !-------  INPUT PARAMETERS:
    !  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
    !    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
    !
    !-------  OUTPUT PARAMETERS:
    !  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
    !  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
    !  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
    !  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
    !
    !  LAST MODIFICATION:  MARCH 31, 2003 (ONLY SOME NOTATION CHANGES)
    !
    !     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
    !-----------------------------------------------------------------------------
    !
    DOUBLE PRECISION DJ,FDAY
    COMMON /CONST/UMR,PI
    !
    IF(IYEAR < 1901.OR.IYEAR > 2099) RETURN
    FDAY=DFLOAT(IHOUR*3600+MIN*60+ISEC)/86400.D0
    DJ=365*(IYEAR-1900)+(IYEAR-1901)/4+IDAY-0.5D0+FDAY
    T=DJ/36525.
    VL=DMOD(279.696678+0.9856473354*DJ,360.D0)
    GST=DMOD(279.690983+.9856473354*DJ+360.*FDAY+180.,360.D0)*UMR
    G=DMOD(358.475845+0.985600267*DJ,360.D0)*UMR
    SLONG=(VL+(1.91946-0.004789*T)*SIN(G)+0.020094*SIN(2.*G))*UMR
    IF(SLONG > 6.2831853) SLONG=SLONG-6.2831853
    IF (SLONG < 0.) SLONG=SLONG+6.2831853
    OBLIQ=(23.45229-0.0130125*T)*UMR
    SOB=SIN(OBLIQ)
    SLP=SLONG-9.924E-5
    !
    !   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION  DUE TO
    !   THE ORBITAL MOTION OF THE EARTH
    !
    SIN1=SOB*SIN(SLP)
    COS1=SQRT(1.-SIN1**2)
    SC=SIN1/COS1
    SDEC=ATAN(SC)
    SRASN=3.141592654-ATAN2(COS(OBLIQ)/SOB*SC,-COS(SLP)/COS1)
    RETURN
END
