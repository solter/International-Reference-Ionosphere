!-----------------------------------------------------------------------
! Functions and subroutines for the International Reference Ionosphere 
! (IRI) model. These functions and subroutines are called by the main
! IRI subroutine IRI_SUB in IRISUB.FOR.
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
                  
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! IRI functions and subroutines:
! MAG. FIELD: FIELDG, CONVER(Geom. Corrected Latitude)
! INDICES:  TCON,APF,APFMSIS,APF_ONLY
! STORM:       CONVER, STORM, STORME_AP
! Spread-F:    spreadf_brazil,bspl2f,bspl2l,bspl2s,bspl4t
! Auroral:     auroral_boundary, ckp
! Misc:     REGFA1
!-----------------------------------------------------------------------
                     

!                     


!                     
!                     
!
!
!                     
!************************************************************                   
!*************** EARTH MAGNETIC FIELD ***********************                   
!**************************************************************                 
!
!
SUBROUTINE FIELDG(DLAT,DLONG,ALT,X,Y,Z,F,DIP,DEC,SMODIP)                  
    ! THIS IS A SPECIAL VERSION OF THE POGO 68/10 MAGNETIC FIELD                    
    ! LEGENDRE MODEL. TRANSFORMATION COEFF. G(144) VALID FOR 1973.                  
    ! INPUT: DLAT, DLONG=GEOGRAPHIC COORDINATES/DEG.(-90/90,0/360),                 
    !        ALT=ALTITUDE/KM.                          
    ! OUTPUT: F TOTAL FIELD (GAUSS), Z DOWNWARD VERTICAL COMPONENT                  
    !        X,Y COMPONENTS IN THE EQUATORIAL PLANE (X TO ZERO LONGITUDE).          
    !        DIP INCLINATION ANGLE(DEGREE). SMODIP RAWER'S MODFIED DIP.             
    ! SHEIK,1977.         
    DIMENSION H(144),XI(3),G(144),FEL1(72),FEL2(72)
    COMMON/CONST/UMR,PI                           
    DATA FEL1/0.0, 0.1506723,0.0101742, -0.0286519, 0.0092606,                &
        -0.0130846, 0.0089594, -0.0136808,-0.0001508, -0.0093977,                &
        0.0130650, 0.0020520, -0.0121956, -0.0023451, -0.0208555,                &
        0.0068416,-0.0142659, -0.0093322, -0.0021364, -0.0078910,                &
        0.0045586,  0.0128904, -0.0002951, -0.0237245,0.0289493,                 &
        0.0074605, -0.0105741, -0.0005116, -0.0105732, -0.0058542,               &
        0.0033268, 0.0078164,0.0211234, 0.0099309, 0.0362792,                     &
        -0.0201070,-0.0046350,-0.0058722,0.0011147,-0.0013949,                    &
        -0.0108838,  0.0322263, -0.0147390,  0.0031247, 0.0111986,               &
        -0.0109394,0.0058112,  0.2739046, -0.0155682, -0.0253272,                &
        0.0163782, 0.0205730,  0.0022081, 0.0112749,-0.0098427,                 &
        0.0072705, 0.0195189, -0.0081132, -0.0071889, -0.0579970,                &
        -0.0856642, 0.1884260,-0.7391512, 0.1210288, -0.0241888,                 &
        -0.0052464, -0.0096312, -0.0044834, 0.0201764,  0.0258343,               &
        0.0083033,  0.0077187/                       
    DATA FEL2/0.0586055,0.0102236,-0.0396107,    &
        -0.0167860, -0.2019911, -0.5810815,0.0379916,  3.7508268,                &
        1.8133030, -0.0564250, -0.0557352, 0.1335347, -0.0142641,                &
        -0.1024618,0.0970994, -0.0751830,-0.1274948, 0.0402073,                  &
        0.0386290, 0.1883088,  0.1838960, -0.7848989,0.7591817,                 &
        -0.9302389,-0.8560960, 0.6633250, -4.6363869, -13.2599277,               &
        0.1002136,  0.0855714,-0.0991981, -0.0765378,-0.0455264,                 &
        0.1169326, -0.2604067, 0.1800076, -0.2223685, -0.6347679,               &
        0.5334222, -0.3459502,-0.1573697,  0.8589464, 1.7815990,                  &
        -6.3347645, -3.1513653, -9.9927750,13.3327637, -35.4897308,               &
        37.3466339, -0.5257398,  0.0571474, -0.5421217,  0.2404770,               &
        -0.1747774,-0.3433644, 0.4829708,0.3935944, 0.4885033,                   &
        0.8488121, -0.7640999, -1.8884945, 3.2930784,-7.3497229,                &
        0.1672821,-0.2306652, 10.5782146, 12.6031065, 8.6579742,                 &
        215.5209961, -27.1419220,22.3405762,1108.6394043/                        
    K=0             
    DO I=1,72    
        K=K+1           
        G(K)=FEL1(I)    
        G(72+K)=FEL2(I)                              
    end do
    RLAT=DLAT*UMR   
    CT=SIN(RLAT)    
    ST=COS(RLAT)    
    NMAX=11         
    D=SQRT(40680925.0-272336.0*CT*CT)            
    RLONG=DLONG*UMR                              
    CP=COS(RLONG)   
    SP=SIN(RLONG)   
    ZZZ=(ALT+40408589.0/D)*CT/6371.2             
    RHO=(ALT+40680925.0/D)*ST/6371.2             
    XXX=RHO*CP      
    YYY=RHO*SP      
    RQ=1.0/(XXX*XXX+YYY*YYY+ZZZ*ZZZ)             
    XI(1)=XXX*RQ    
    XI(2)=YYY*RQ    
    XI(3)=ZZZ*RQ    
    IHMAX=NMAX*NMAX+1                            
    LAST=IHMAX+NMAX+NMAX                         
    IMAX=NMAX+NMAX-1                             
    DO I=IHMAX,LAST                          
        H(I)=G(I)       
    end do
    DO K=1,3,2  
        I=IMAX          
        IH=IHMAX        
        300   IL=IH-I         
        F1=2./(I-K+2.)  
        X1=XI(1)*F1     
        Y1=XI(2)*F1     
        Z1=XI(3)*(F1+F1)                             
        I=I-2           
        IF((I-1) < 0) GOTO 400                      
        IF((I-1) == 0) GOTO 500                      
        DO M=3,I,2  
            H(IL+M+1)=G(IL+M+1)+Z1*H(IH+M+1)+X1*(H(IH+M+3)-H(IH+M-1))-                &
                Y1*(H(IH+M+2)+H(IH+M-2))                     
            H(IL+M)=G(IL+M)+Z1*H(IH+M)+X1*(H(IH+M+2)-H(IH+M-2))+                      &
                Y1*(H(IH+M+3)+H(IH+M-1))                     
        end do
        500   H(IL+2)=G(IL+2)+Z1*H(IH+2)+X1*H(IH+4)-Y1*(H(IH+3)+H(IH))                  
        H(IL+1)=G(IL+1)+Z1*H(IH+1)+Y1*H(IH+4)+X1*(H(IH+3)-H(IH))                  
        400   H(IL)=G(IL)+Z1*H(IH)+2.0*(X1*H(IH+1)+Y1*H(IH+2))                          
        700   IH=IL           
        IF(I >= K) GOTO 300                          
    end do
    S=0.5*H(1)+2.0*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))                         
    XT=(RQ+RQ)*SQRT(RQ)                          
    X=XT*(H(3)-S*XXX)                            
    Y=XT*(H(4)-S*YYY)                            
    Z=XT*(H(2)-S*ZZZ)                            
    F=SQRT(X*X+Y*Y+Z*Z)                          
    BRH0=Y*SP+X*CP  
    Y=Y*CP-X*SP     
    X=Z*ST-BRH0*CT  
    Z=-Z*CT-BRH0*ST 
    zdivf=z/f
    IF(ABS(zdivf) > 1.) zdivf=SIGN(1.,zdivf)
    DIP=ASIN(zdivf)
    ydivs=y/sqrt(x*x+y*y)  
    IF(ABS(ydivs) > 1.) ydivs=SIGN(1.,ydivs)
    DEC=ASIN(ydivs)
    dipdiv=DIP/SQRT(DIP*DIP+ST)
    IF(ABS(dipdiv) > 1.) dipdiv=SIGN(1.,dipdiv)
    SMODIP=ASIN(dipdiv)
    DIP=DIP/UMR     
    DEC=DEC/UMR     
    SMODIP=SMODIP/UMR                            
    RETURN          
END             
!
!
!************************************************************                   
!*********** INTERPOLATION AND REST ***************************                 
!**************************************************************                 
!
!
SUBROUTINE REGFA1(X11,X22,FX11,FX22,EPS,FW,F,SCHALT,X) 
    ! REGULA-FALSI-PROCEDURE TO FIND X WITH F(X)-FW=0. X1,X2 ARE THE                
    ! STARTING VALUES. THE COMUTATION ENDS WHEN THE X-INTERVAL                      
    ! HAS BECOME LESS THEN EPS . IF SIGN(F(X1)-FW)= SIGN(F(X2)-FW)                  
    ! THEN SCHALT=.TRUE.  
    LOGICAL L1,LINKS,K,SCHALT                    
    SCHALT=.FALSE.
    EP=EPS  
    X1=X11          
    X2=X22          
    F1=FX11-FW     
    F2=FX22-FW     
    K=.FALSE.       
    NG=2       
    LFD=0  
    IF(F1*F2 <= 0.0) GOTO 200
    X=0.0           
    SCHALT=.TRUE.   
    RETURN
    200   X=(X1*F2-X2*F1)/(F2-F1) 
    GOTO 400        
    300     L1=LINKS        
    DX=(X2-X1)/NG
    IF(.NOT.LINKS) DX=DX*(NG-1)
    X=X1+DX
    400   FX=F(X)-FW
    LFD=LFD+1
    IF(LFD > 20) THEN
        EP=EP*10.
        LFD=0
    ENDIF 
    LINKS=(F1*FX > 0.0)
    K=.NOT.K        
    IF(LINKS) THEN
        X1=X            
        F1=FX           
    ELSE
        X2=X 
        F2=FX 
    ENDIF   
    IF(ABS(X2-X1) <= EP) GOTO 800               
    IF(K) GOTO 300  
    IF((LINKS.AND.(.NOT.L1)).OR.(.NOT.LINKS.AND.L1)) NG=2*NG                  
    GOTO 200        
    800   RETURN          
END             
!
!
  
!
!
subroutine read_ig_rz 
    !----------------------------------------------------------------
    ! Reads the Rz12 and IG12 indices file IG_RZ.DAT from I/O UNIT=12 
    ! and stores the indices in COMMON:
    !        common/igrz/aig,arziyst,iyed   with aig(806),arz(806),
    !                                            start year (iyst)
    !                                           end year (iyed)
    ! 
    ! The indices file IG_RZ.DAT is structured as follows (values are 
    ! separated by comma): 
    !   day, month, year of the last update of this file,
    !   a blank line
    !   start month, start year, end month, end year,
    !   a blank line
    !   the IG index for December of start year minus 1 (this value is 
    !        needed for interpolating from 1st to 15th of first year)
    !   the 12 IG indices (13-months running mean) for start year, 
    !   the 12 IG indices for the second year 
    !       .. and so on until the last year,
    !   the 12 IG indices for the last year 
    !   the IG index for January of end year plus 1 (needed for interpolation)
    !   a blank line
    !   the Rz index for December of start year minus 1 
    !   the 12 Rz indices (13-months running mean) for the start year,
    !   the 12 Rz indices for the second year 
    !       .. and so on until the last year.
    !   the 12 Rz indices for the last year 
    !   the Rz index for January of end year plus 1 
    ! 
    ! A negative Rz index means that the given index is the 13-months-
    ! running mean of the solar radio flux (F10.7). The close correlation 
    ! between (Rz)12 and (F10.7)12 is used to compute the (Rz)12 indices.
    !
    ! An IG index of -111 indicates that no IG values are available for the
    ! time period. In this case a correlation function between (IG)12 and 
    ! (Rz)12 is used to obtain (IG)12.
    !
    ! The computation of the 13-month-running mean for month M requires the
    ! indices for the six months preceeding M and the six months following 
    ! M (month: M-6, ..., M+6). To calculate the current running mean one 
    ! therefore requires predictions of the indix for the next six months. 
    ! Starting from six months before the UPDATE DATE (listed at the top of 
    ! the file) and onward the indices are therefore based on indices 
    ! predictions.
    !----------------------------------------------------------------
    integer    iyst,iyend,iymst,iupd,iupm,iupy,imst,imend
    real        aig(806),arz(806)
    
    common /igrz/aig,arz,iymst,iymend
    open(unit=12,file='ig_rz.dat',FORM='FORMATTED',status='old')
    !-web- special for web version
    !            open(unit=12,file=
    !     *         '/var/www/omniweb/cgi/vitmo/IRI/ig_rz.dat',
    !     *         FORM='FORMATTED',status='old')
    ! Read the update date, the start date and the end date (mm,yyyy), and
    ! get number of data points to read.
    read(12,*) iupd,iupm,iupy
    read(12,*) imst,iyst,imend,iyend
    iymst=iyst*100+imst
    iymend=iyend*100+imend
    ! inum_vals= 12-imst+1+(iyend-iyst-1)*12 +imend + 2
    ! 1st year \ full years       \last y\ before   after
    inum_vals= 3-imst+(iyend-iyst)*12 +imend
    ! read all the IG12 (ionoindx) and Rz12 (indrz) values
    read(12,*) (aig(i),i=1,inum_vals)
    read(12,*) (arz(i),i=1,inum_vals)
    !            do 1 jj=1,inum_vals
    !                rrr=arz(jj)
    !                ggg=aig(jj)
    !                if(rrr < 0.0) then
    !                    covr=abs(rrr)
    !                    rrr=33.52*sqrt(covr+85.12)-408.99
    !                    if(rrr < 0.0) rrr=0.0
    !                    endif
    !                if(ggg <= -90.) then
    !                    zi=-12.349154+(1.4683266-2.67690893e-03*rrr)*rrr
    !                    if(zi > 274.0) zi=274.0
    !                    ggg=zi
    !                    endif
    !                arz(jj)=rrr
    !                aig(jj)=ggg    
    !1               continue
    close(unit=12)
    return
end
!
!
subroutine tcon(yr,mm,day,idn,rz,ig,rsn,nmonth)
    !----------------------------------------------------------------
    ! input:        yr,mm,day       year(yyyy),month(mm),day(dd)
    !               idn             day of year(ddd)
    ! output:       rz(3)           12-month-smoothed solar sunspot number
    !               ig(3)           12-month-smoothed IG index
    !               rsn             interpolation parameter
    !               nmonth          previous or following month depending
    !                               on day
    !
    ! Uses read_ig_rz and common/igrz/ to get indices 
    ! 
    ! rz(1)   ig(1) contain the indices for the month mm and rz(2) & ig(2)
    ! for the previous month (if day less than 15) or for the following
    ! month (if day greater than 15). These indices are for the mid of the 
    ! month. The indices for the given day are obtained by linear 
    ! interpolation and are stored in rz(3) and ig(3).
    !----------------------------------------------------------------
    integer    yr,mm,day,iyst,iyend,iymst
    integer    imst,iymend
    real        ionoindx(806),indrz(806)
    real        ig(3),rz(3)
    logical    mess
    
    common     /iounit/konsol,mess  
    common    /igrz/ionoindx,indrz,iymst,iymend
    iytmp=yr*100+mm
    if (iytmp < iymst.or.iytmp > iymend) then
        if(mess) write(konsol,8000) iytmp,iymst,iymend
        8000           format(1x,I10,'** OUT OF RANGE **'/,5x,&
            'The file IG_RZ.DAT which contains the indices Rz12',&
            ' and IG12'/5x,'currently only covers the time period',&
            ' (yymm) : ',I6,'-',I6)
        nmonth=-1
        return
    endif
    iyst=iymst/100
    imst=iymst-iyst*100
    !       num=12-imst+1+(yr-iyst-1)*12+mm+1
    num=2-imst+(yr-iyst)*12+mm
    rz(1)=indrz(num)
    ig(1)=ionoindx(num)
    midm=15
    if(mm == 2) midm=14
    call MODA(0,yr,mm,midm,idd1,nrdaym)
    if(day < midm) goto 1926
    ! day is at or after mid of month
    imm2=mm+1
    if(imm2 > 12) then
        imm2=1
        iyy2=yr+1
        idd2=380            ! =365+15 mid-January
        !               if((yr/4*4 == yr).and.(yr/100*100 /= yr)) idd2=381
        if(yr/4*4 == yr) idd2=381
    else
        iyy2=yr
        midm=15
        if(imm2 == 2) midm=14
        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
    endif
    rz(2)=indrz(num+1)
    ig(2)=ionoindx(num+1)
    rsn=(idn-idd1)*1./(idd2-idd1)                
    rz(3)=rz(1)+(rz(2)-rz(1))*rsn
    ig(3)=ig(1)+(ig(2)-ig(1))*rsn
    goto 1927
    1926            imm2=mm-1
    if(imm2 < 1) then
        imm2=12
        idd2=-16
        iyy2=yr-1
    else
        iyy2=yr
        midm=15
        if(imm2 == 2) midm=14
        call MODA(0,iyy2,imm2,midm,IDD2,nrdaym)
    endif
    rz(2)=indrz(num-1)
    ig(2)=ionoindx(num-1)
    rsn=(idn-idd2)*1./(idd1-idd2)
    rz(3)=rz(2)+(rz(1)-rz(2))*rsn
    ig(3)=ig(2)+(ig(1)-ig(2))*rsn
    1927    nmonth=imm2
    return
end
!
!
subroutine readapf107
    !-------------------------------------------------------------------------
    ! Reads APF107.DAT file (on UNIT=13) and stores contents in COMMON block:
    !     COMMON/AAP,AF107,N/ with  AAP(23000,9) and AF107(23000,3)
    !        AAP(*,1)    3-hour Ap indices for the UT interval )0-3)
    !        AAP(*,2)    3-hour Ap indices for the UT interval )3-6)
    !          ....                       ....
    !        AAP(*,8)    3-hour Ap indices for the UT interval )21-6)
    !        AAP(*,9)    daily Ap
    !        AF107(*,1)    F10.7 radio flux for the day
    !        AF107(*,2)    81-day average of F10.7 radio flux 
    !        AF107(*,3)    365-day average of F10.7
    !       N           total number of records
    !
    ! APF107.DAT is structured as follows:
    !         JY(I3),JMN(I3),JD(I3)    year, month, day 
    !        IIAP(8)    (8I3)            3-hour Ap indices for the UT intervals 
    !                                (0-3(,(3-6(,(6-9(, .., (18-21(,(21-24(
    !        IAPD (I3)                daily Ap
    !        IR (I3)                    sunspot number for the day (empty)
    !        F107 (F5.1)                F10.7 radio flux for the day
    !        F107_81 (F5.1)            81-day average of F10.7 radio flux 
    !       F107_365 (F5.1)         365-day average of F10.7 centered on 
    !                               the date of interest. At start and end  
    !                                of index file it takes all available  
    !                               indices, e.g. for the first date the 
    !                               average is only over 40 F10.7 values  
    !                               and over 41 values on the 2nd date.  
    !
    ! If date is outside the range of the Ap indices file then IAP(1)=-5  
    !-------------------------------------------------------------------------
    !
    INTEGER        aap(23000,9),iiap(8)
    DIMENSION     af107(23000,3)
    COMMON        /apfa/aap,af107,n
    Open(13,FILE='apf107.dat',FORM='FORMATTED',STATUS='OLD')
    !-web-sepcial vfor web version
    !      OPEN(13,FILE='/var/www/omniweb/cgi/vitmo/IRI/apf107.dat',
    !     *    FORM='FORMATTED',STATUS='OLD')
    i=1
    1       READ(13,10,END=21) JY,JMN,JD,iiap,iapda,IR,F107D,F107_81,&
        F107_365
    10      FORMAT(3I3,9I3,I3,3F5.1)
    !         adate(i)=jy*10000+jmn*100+jd 
    do j=1,8 
        aap(i,j)=iiap(j)
    enddo
    aap(i,9)=iapda
    !        irza(i)=ir
    if(F107_81 < -4.) F107_81=F107D
    if(F107_365 < -4.) F107_365=F107D
    af107(i,1)=f107d 
    af107(i,2)=f107_81 
    af107(i,3)=f107_365
    i=i+1 
    goto 1
    21        n=i-1
    
    CLOSE(13)
    return
end
!
!
SUBROUTINE APF(ISDATE,HOUR,IAP)
    !-----------------------------------------------------------------------
    ! Finds 3-hourly Ap indices for IRI-STORM model
    !    INPUTS:     ISDATE        Array-index from APF_ONLY
    !                HOUR        UT in decimal hours
    !    OUTPUT:    IAP(1:13)    3-hourly Ap index
    !                            IAP(13) Ap index for current UT
    !                            IAP(1) AP index for UT-39 hours.
    !
    ! Gets indices from COMMON/APFA/
    !
    ! If date is outside the range of the Ap indices file than IAP(1)=-5  
    !-----------------------------------------------------------------------
    INTEGER        aap(23000,9),iiap(8),iap(13)
    DIMENSION     af107(23000,3)
    LOGICAL     mess
    COMMON         /iounit/konsol,mess    /apfa/aap,af107,nf107
    
    do i=1,8
        iap(i)=-1
    enddo
    IS = ISDATE
    ihour=int(hour/3.)+1
    if(ihour > 8) ihour=8
    if(is*8+ihour < 13) goto 21   ! less then 13 indices available    
    j1=13-ihour
    do i=1,ihour
        iapi=aap(is,i)
        if(iapi < -2) goto 21
        iap(j1+i)=iapi
    enddo
    if(ihour > 4) then
        do i=1,j1
            iapi=aap(is-1,8-j1+i)
            if(iapi < -2) goto 21
            iap(i)=iapi
        enddo
    else           
        j2=5-ihour
        do i=1,8
            iapi=aap(is-1,i)
            if(iapi < -2) goto 21
            iap(j2+i)=iapi
        enddo
        do i=1,j2
            iapi=aap(is-2,8-j2+i)
            if(iapi < -2) goto 21
            iap(i)=iapi
        enddo
    endif         
    goto 20
    
    21      if(mess) write(konsol,100)
    100     format(1X,'One of the ap indeces is negative.',&
        ' STORM model is turned off.')
    IAP(1)=-5
    
    20    RETURN
END
!
!
SUBROUTINE APFMSIS(ISDATE,HOUR,IAPO)
    !-----------------------------------------------------------------------
    ! Finds 3-hourly Ap indices for NRLMSIS00 model 
    !    INPUTS:     ISDATE        Array-index from APF_ONLY
    !                HOUR        UT in decimal hours
    !    OUTPUT:       IAPO(1:7)    3-hourly Ap index
    !  
    ! IAPO(1) DAILY AP
    ! IAPO(2) 3-HR AP INDEX FOR CURRENT TIME                   
    ! IAPO(3) 3-HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME    
    ! IAPO(4) 3-HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME    
    ! IAPO(5) 3-HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME  
    ! IAPO(6) AVERAGE OF EIGHT 3-HR AP INDICIES FROM 12 TO 33 HRS PRIOR
    !         TO CURRENT TIME
    ! IAPO(7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
    !         TO CURRENT TIME
    !
    ! The 3-hour UT intervals during the day are: (0-3),)3-6),)6-9),)9-12),
    ! )12-15),)15-18),)18-21),)21-24(.
    ! 
    ! If date is outside the range of the Ap indices file then IAPO(2)=-5  
    !-----------------------------------------------------------------------
    !
    REAL         IAPO
    INTEGER        aap(23000,9),iiap(8)
    DIMENSION     af107(23000,3),iap(20),lm(12),iapo(7)
    LOGICAL      mess
    COMMON         /iounit/konsol,mess    /apfa/aap,af107,nf107
    IS=ISDATE
    ihour=int(hour/3.)+1
    if(ihour > 8) ihour=8
    iapo(1)=aap(is,9)        
    ! There must be at least 20 indices available
    if((is-1)*8+ihour < 20) goto 21          
    ! assemble Ap values as needed by MSIS
    j1=ihour+1
    do i=1,ihour
        iap(j1-i)=aap(is,i)
    enddo
    j1=ihour+9
    do i=1,8
        iap(j1-i)=aap(is-1,i)
    enddo
    j1=ihour+17
    j2=8-(20-ihour-8)+1
    if(j2 < 1) j2=1
    do i=j2,8
        iap(j1-i)=aap(is-2,i)
    enddo
    if(ihour < 4) then
        j1=ihour+25
        j2=8-(20-ihour-16)+1
        do i=j2,8
            iap(j1-i)=aap(is-3,i)
        enddo
    endif
    do i=1,4 
        iapo(i+1)=iap(i)*1.0
    end do
    sum1=0.
    sum2=0.
    do i=1,8
        sum1=sum1+iap(4+i)
        sum2=sum2+iap(12+i)
    end do
    !        iapo(6)=int(sum1/8.+.5)
    !        iapo(7)=int(sum2/8.+.5)
    iapo(6)=sum1/8.
    iapo(7)=sum2/8.
    goto 20
    
    21      if(mess) write(konsol,100)
    100     format(1X,'APFMSIS: No Ap dependence because date is not',&
        ' covered by APF107.DAT indices file')
    IAPO(2)=-5.0
    20        continue
    RETURN
END
!
!
SUBROUTINE APF_ONLY(IYYYY,IMN,ID,F107D,F107PD,F107_81,F107_365,&
        IAPDA,ISDATE)
    !-----------------------------------------------------------------------
    ! Finds daily F10.7, daily Ap, and 81-day and 365-day F10.7 index: 
    !
    !    INPUTS:   IYYYY (yyyy)    year 
    !              IMN (mm)        month 
    !               ID (dd)        day 
    !    OUTPUT:   F107D        F10.7 index for the day (adjusted 
    !                                to 1AU)
    !              F107PD          F10.7 index for one day prior (used in MSIS)
    !              F107_81        F10.7 average over 3 solar rotations
    !                               (81 days, centered on the current day) 
    !              F107_365     F10.7 12-month running mean
    !              IAPDA        Daily Ap
    !              ISDATE        Array-index for the specified date (for
    !                                use in APF subroutine.
    ! 
    ! Using COMMON/apfa/ for indices
    !
    ! Is used for vdrift and foeedi.
    !
    ! If date is outside the range of indices file than F107D=F107_81=-11.1  
    !-----------------------------------------------------------------------
    INTEGER        aap(23000,9),iiap(8),lm(12)
    DIMENSION     af107(23000,3)
    LOGICAL     mess
    common         /iounit/konsol,mess /apfa/aap,af107,nf107
    DATA LM/31,28,31,30,31,30,31,31,30,31,30,31/
    IYBEG=1958
    if(iyyyy < IYBEG) goto 21   ! APF107.DAT starts at Jan 1, 1958
    is=0
    do i=IYBEG,iyyyy-1
        nyd=365
        if(i/4*4 == i) nyd=366    ! leap year
        IS=IS+nyd
    enddo
    lm(2)=28
    if(iyyyy/4*4 == iyyyy) lm(2)=29      ! leap year
    do i=1,IMN-1
        IS=IS+LM(i)
    ENDDO
    
    IS=IS+ID
    ISDATE = IS
    if(IS > nf107) goto 21
    F107D = AF107(IS,1)
    F107PD = F107D
    if(IS > 1) F107PD = AF107(IS-1,1)
    F107_81=AF107(IS,2)
    if(F107_81 < -4.) F107_81=F107D
    F107_365=AF107(IS,3)
    if(F107_365 < -4.) F107_365=F107D
    IAPDA=AAP(is,9)
    goto 20
    
    21      if(mess) write(konsol,100)
    100     format(1X,'APF_ONLY: Date is outside range of F10.7D indices',&
        ' file (F10.7D = F10.7_81 = F10.7RM12).')
    F107D = -11.1
    F107_81 = -11.1
    F107_365 = -11.1
    IAPDA = -11
    
    20    RETURN
END
!      
!
!----------------------STORM MODEL --------------------------------
!
SUBROUTINE CONVER(rga,rgo,rgma)
    !     This subroutine converts a geographic latitude and longitude
    !     location to a corrected geomagnetic latitude.
    !
    !     INPUT: 
    !       geographic latitude   -90. to +90.
    !       geographic longitude  0. to 360. positive east from Greenwich.
    !
    !     OUTPUT:
    !       corrected geomagnetic latitude    -90. to +90.
    DIMENSION CORMAG(20,91)      
    DATA ((CORMAG(i,j),i=1,20),j=1,31)/&
        163.68,163.68,163.68,163.68,163.68,163.68,&
        163.68,163.68,163.68,163.68,163.68,163.68,163.68,163.68,&
        163.68,163.68,163.68,163.68,163.68,163.68,162.60,163.12,&
        163.64,164.18,164.54,164.90,165.16,165.66,166.00,165.86,&
        165.20,164.38,163.66,162.94,162.42,162.00,161.70,161.70,&
        161.80,162.14,161.20,162.18,163.26,164.44,165.62,166.60,&
        167.42,167.80,167.38,166.82,166.00,164.66,163.26,162.16,&
        161.18,160.40,159.94,159.80,159.98,160.44,159.80,161.14,&
        162.70,164.50,166.26,167.90,169.18,169.72,169.36,168.24,&
        166.70,164.80,162.90,161.18,159.74,158.60,157.94,157.80,&
        157.98,158.72,158.40,160.10,162.02,164.28,166.64,169.00,&
        170.80,171.72,171.06,169.46,167.10,164.64,162.18,160.02,&
        158.20,156.80,156.04,155.80,156.16,157.02,157.00,158.96,&
        161.24,163.86,166.72,169.80,172.42,173.72,172.82,170.34,&
        167.30,164.22,161.34,158.74,156.60,155.00,154.08,153.90,&
        154.36,155.36,155.50,157.72,160.36,163.32,166.60,170.20,&
        173.70,175.64,174.18,170.80,167.10,163.56,160.24,157.36,&
        154.96,153.10,152.08,151.92,152.46,153.76,154.10,156.52,&
        159.36,162.52,166.24,170.30,174.62,177.48,175.04,170.82,&
        166.60,162.70,159.02,155.88,153.22,151.20,150.08,149.92,&
        150.64,152.20,152.80,155.32,158.28,161.70,165.58,170.00,&
        174.84,178.46,175.18,170.38,165.80,161.64,157.80,154.38,&
        151.52,149.30,148.18,148.02,148.92,150.60,151.40,154.08,&
        157.18,160.68,164.78,169.40,174.34,177.44,174.28,169.44,&
        164.70,160.34,156.30,152.78,149.72,147.40,146.18,146.04,&
        147.12,149.04,150.10,152.88,156.00,159.58,163.78,168.50,&
        173.28,175.60,172.86,168.14,163.40,158.98,154.88,151.10,&
        147.98,145.50,144.18,144.14,145.40,147.48,148.80,151.68,&
        154.88,158.48,162.68,167.40,171.76,173.60,171.12,166.68,&
        162.00,157.48,153.28,149.50,146.18,143.50,142.18,142.24,&
        143.68,145.98,147.50,150.54,153.68,157.28,161.42,166.10,&
        170.10,171.48,169.22,164.98,160.40,155.88,151.68,147.80,&
        144.34,141.60,140.18,140.26,141.98,144.62,146.30,149.34,&
        152.48,155.98,160.08,164.60,168.34,169.38,167.20,163.18,&
        158.60,154.18,149.98,146.02,142.54,139.70,138.18,138.46,&
        140.26,143.16,145.10,148.14,151.18,154.60,158.68,163.10,&
        166.48,167.28,165.18,161.32,156.90,152.48,148.28,144.32,&
        140.74,137.80,136.22,136.48,138.64,141.76,143.90,146.98,&
        149.98,153.30,157.24,161.40,164.52,165.16,162.86,159.42,&
        155.00,150.68,146.48,142.52,138.94,135.90,134.22,134.68,&
        137.02,140.40,142.70,145.84,148.76,151.92,155.74,159.70,&
        162.52,162.96,160.98,157.42,153.10,148.84,144.68,140.82,&
        137.20,134.00,132.32,132.80,135.42,139.10,141.60,144.74,&
        147.46,150.52,154.20,158.00,160.46,160.76,158.86,155.36,&
        151.20,146.94,142.88,139.02,135.40,132.10,130.32,131.00,&
        133.80,137.74,140.50,143.58,146.24,149.12,152.60,156.20,&
        158.40,158.66,156.76,153.36,149.30,145.04,141.08,137.30,&
        133.60,130.30,128.42,129.12,132.28,136.44,139.30,142.48,&
        144.94,147.64,150.48,154.30,156.34,156.36,154.56,151.26,&
        147.30,143.14,139.20,135.50,131.90,128.40,126.52,127.32,&
        130.76,135.18,138.20,141.28,143.72,146.24,149.26,152.40,&
        154.24,154.16,152.36,149.16,145.30,141.24,137.30,133.70,&
        130.10,126.60,124.62,125.54,129.16,133.92,137.10,140.18,&
        142.42,144.66,147.62,150.50,152.18,151.96,150.16,147.10,&
        143.30,139.24,135.50,131.90,128.36,124.80,122.72,123.74,&
        127.64,132.62,135.90,139.02,141.12,143.18,145.92,148.60,&
        149.98,149.76,148.04,145.00,141.20,137.30,133.60,130.10,&
        126.60,123.00,120.86,121.96,126.12,131.36,134.80,137.88,&
        139.80,141.68,144.08,146.60,147.88,147.56,145.84,142.90,&
        139.20,135.30,131.70,128.28,124.86,121.30,118.96,120.18,&
        124.70,130.16,133.60,136.72,138.48,140.10,142.38,144.60,&
        145.72,145.34,143.64,140.80,137.10,133.30,129.72,126.48,&
        123.10,119.50,117.16,118.48,123.18,128.86,132.40,135.42,&
        137.08,138.50,140.54,142.60,143.52,143.06,141.44,138.70,&
        135.10,131.30,127.82,124.58,121.40,117.70,115.26,116.70,&
        121.66,127.60,131.20,134.22,135.66,136.82,138.70,140.60,&
        141.36,140.86,139.24,136.50,133.00,129.30,125.92,122.78,&
        119.60,116.00,113.40,114.92,120.16,126.30,130.00,132.92,&
        134.24,135.14,136.80,138.60,139.16,138.64,137.12,134.40,&
        130.90,127.20,123.92,120.96,117.90,114.20,111.56,113.12,&
        118.64,124.90,128.70,131.56,132.74,133.44,134.90,136.50,&
        137.00,136.36,134.82,132.30,128.70,125.16,121.94,119.06,&
        116.10,112.50,109.70,111.42,117.14,123.60,127.30,130.16,&
        131.22,131.66,133.00,134.50,134.80,134.14,132.62,130.14,&
        126.60,123.06,119.94,117.16,114.30,110.70,107.80,109.64,&
        115.62,122.24,125.90,128.76,129.62,129.96,131.06,132.40,&
        132.60,131.86,130.42,128.00,124.50,120.96,117.96,115.26,&
        112.54,108.90,105.94,107.86,114.02,120.84/
    DATA ((CORMAG(i,j),i=1,20),j=32,61)/&
        124.05,126.79,&
        127.55,127.83,128.90,130.21,130.41,129.71,128.33,125.96,&
        122.49,118.96,115.97,113.26,110.52,106.89,104.01,106.00,&
        112.21,119.06,122.19,124.82,125.48,125.69,126.73,128.03,&
        128.22,127.55,126.23,123.92,120.47,116.97,113.97,111.26,&
        108.50,104.89,102.08,104.14,110.41,117.29,120.34,122.85,&
        123.40,123.56,124.57,125.84,126.03,125.40,124.14,121.88,&
        118.46,114.97,111.98,109.26,106.48,102.88,100.15,102.28,&
        108.60,115.51,118.49,120.88,121.33,121.42,122.40,123.65,&
        123.84,123.24,122.04,119.83,116.45,112.97,109.98,107.26,&
        104.46,100.87,098.22,100.42,106.79,113.74,116.63,118.91,&
        119.26,119.29,120.24,121.47,121.65,121.09,119.95,117.79,&
        114.43,110.98,107.99,105.26,102.44,098.87,096.29,098.56,&
        104.98,111.96,114.78,116.94,117.19,117.15,118.07,119.28,&
        119.46,118.93,117.86,115.75,112.42,108.98,106.00,103.26,&
        100.42,096.86,094.36,096.70,103.18,110.19,112.93,114.97,&
        115.12,115.02,115.91,117.09,117.27,116.78,115.76,113.71,&
        110.41,106.98,104.00,101.26,098.40,094.85,092.43,094.84,&
        101.37,108.41,111.07,113.00,113.04,112.88,113.74,114.91,&
        115.08,114.62,113.67,111.67,108.39,104.99,102.01,099.26,&
        096.38,092.85,090.51,092.97,099.56,106.64,109.22,111.03,&
        110.97,110.75,111.58,112.72,112.89,112.47,111.57,109.63,&
        106.38,102.99,100.01,097.26,094.36,090.84,088.58,091.11,&
        097.75,104.86,107.37,109.06,108.90,108.61,109.41,110.53,&
        110.70,110.31,109.48,107.59,104.37,100.99,098.02,095.26,&
        092.34,088.83,086.65,089.25,095.95,103.09,105.51,107.09,&
        106.83,106.48,107.25,108.35,108.51,108.16,107.39,105.55,&
        102.35,099.00,096.03,093.26,090.32,086.83,084.72,087.39,&
        094.14,101.31,103.66,105.12,104.76,104.34,105.08,106.16,&
        106.32,106.00,105.29,103.50,100.34,097.00,094.03,091.26,&
        088.30,084.82,082.79,085.53,092.33,099.54,101.81,103.15,&
        102.68,102.21,102.92,103.97,104.13,103.85,103.20,101.46,&
        098.33,095.00,092.04,089.26,086.28,082.81,080.86,083.67,&
        090.52,097.76,099.95,101.18,100.61,100.07,100.75,101.79,&
        101.94,101.69,101.10,099.42,096.31,093.01,090.04,087.26,&
        084.26,080.81,078.93,081.81,088.72,095.99,098.10,099.21,&
        098.54,097.94,098.59,099.60,099.75,099.54,099.01,097.38,&
        094.30,091.01,088.05,085.26,082.24,078.80,077.00,079.95,&
        086.91,094.21,096.25,097.24,096.47,095.81,096.43,097.41,&
        097.56,097.39,096.92,095.34,092.29,089.01,086.06,083.26,&
        080.22,076.79,075.07,078.09,085.10,092.43,094.39,095.27,&
        094.40,093.67,094.26,095.23,095.37,095.23,094.82,093.30,&
        090.27,087.02,084.06,081.26,078.20,074.79,073.14,076.23,&
        083.30,090.66,092.54,093.30,092.32,091.54,092.10,093.04,&
        093.18,093.08,092.73,091.26,088.26,085.02,082.07,079.26,&
        076.18,072.78,071.21,074.37,081.49,088.88,090.69,091.33,&
        090.25,089.40,089.93,090.85,090.99,090.92,090.63,089.21,&
        086.25,083.02,080.07,077.26,074.16,070.77,069.28,072.51,&
        079.68,087.11,088.83,089.36,088.18,087.27,087.77,088.67,&
        088.80,088.77,088.54,087.17,084.23,081.03,078.08,075.26,&
        072.14,068.77,067.35,070.65,077.87,085.33,086.98,087.39,&
        086.11,085.13,085.60,086.48,086.61,086.61,086.45,085.13,&
        082.22,079.03,076.09,073.26,070.12,066.76,065.42,068.79,&
        076.07,083.56,085.13,085.42,084.04,083.00,083.44,084.29,&
        084.42,084.46,084.35,083.09,080.21,077.03,074.09,071.26,&
        068.10,064.75,063.49,066.93,074.26,081.78,083.27,083.45,&
        081.96,080.86,081.27,082.11,082.23,082.30,082.26,081.05,&
        078.19,075.04,072.10,069.26,066.08,062.75,061.57,065.06,&
        072.45,080.01,081.42,081.48,079.89,078.73,079.11,079.92,&
        080.04,080.15,080.16,079.01,076.18,073.04,070.10,067.26,&
        064.06,060.74,059.64,063.20,070.64,078.23,079.57,079.51,&
        077.82,076.59,076.94,077.73,077.85,077.99,078.07,076.97,&
        074.17,071.04,068.11,065.26,062.04,058.73,057.71,061.34,&
        068.84,076.46,077.71,077.54,075.75,074.46,074.78,075.55,&
        075.66,075.84,075.98,074.93,072.15,069.05,066.12,063.26,&
        060.02,056.73,055.78,059.48,067.03,074.68,075.86,075.57,&
        073.68,072.32,072.61,073.36,073.47,073.68,073.88,072.88,&
        070.14,067.05,064.12,061.26,058.00,054.72,053.85,057.62,&
        065.22,072.91,074.01,073.60,071.60,070.19,070.45,071.17,&
        071.28,071.53,071.79,070.84,068.13,065.05,062.13,059.26,&
        055.98,052.71,051.92,055.76,063.41,071.13,072.15,071.63,&
        069.53,068.05,068.28,068.99,069.09,069.37,069.69,068.80,&
        066.11,063.06,060.13,057.26,053.96,050.71,049.99,053.90,&
        061.61,069.36,070.30,069.66,067.46,065.92,066.12,066.80,&
        066.90,067.22,067.60,066.76,064.10,061.06,058.14,055.26,&
        051.94,048.70,048.06,052.04,059.80,067.58/
    DATA ((CORMAG(i,j),i=1,20),j=62,91)/&
        067.70,067.06,&
        065.08,063.72,063.98,064.60,064.80,065.12,065.60,064.86,&
        062.40,059.26,056.24,053.18,049.84,046.60,046.12,050.12,&
        057.52,064.80,064.90,064.42,062.70,061.62,061.78,062.40,&
        062.60,063.04,063.58,063.00,060.60,057.46,054.42,051.18,&
        047.70,044.60,044.22,048.02,055.06,061.92,062.10,061.72,&
        060.32,059.50,059.68,060.20,060.46,060.94,061.58,061.00,&
        058.70,055.66,052.52,049.18,045.60,042.50,042.22,046.00,&
        052.60,058.98,059.20,059.18,058.12,057.32,057.48,058.00,&
        058.30,058.84,059.48,059.04,056.90,053.86,050.62,047.10,&
        043.50,040.50,040.28,043.98,050.22,056.18,056.40,056.64,&
        055.84,055.20,055.38,055.80,056.16,056.84,057.48,057.04,&
        055.10,052.06,048.70,045.10,041.40,038.40,038.28,041.88,&
        047.94,053.44,053.70,054.14,053.56,053.10,053.24,053.70,&
        054.06,054.74,055.38,055.14,053.20,050.26,046.80,043.10,&
        039.34,036.40,036.38,039.96,045.56,050.84,051.10,051.70,&
        051.36,051.00,051.14,051.50,051.96,052.64,053.38,053.08,&
        051.30,048.36,044.90,041.02,037.24,034.40,034.38,037.86,&
        043.28,048.20,048.50,049.26,049.18,048.90,049.04,049.40,&
        049.86,050.64,051.28,051.08,049.40,046.46,042.98,039.02,&
        035.14,032.40,032.48,035.72,041.00,045.70,046.00,046.96,&
        046.98,046.80,046.94,047.30,047.76,048.54,049.28,049.08,&
        047.40,044.56,041.08,037.02,033.14,030.40,030.58,033.84,&
        038.72,043.20,043.50,044.62,044.80,044.80,044.94,045.20,&
        045.76,046.54,047.18,046.98,045.50,042.66,039.08,035.02,&
        031.14,028.40,028.58,031.82,036.52,040.80,041.20,042.32,&
        042.54,042.70,042.84,043.20,043.66,044.44,045.08,044.98,&
        043.50,040.76,037.08,033.04,029.04,026.40,026.68,029.82,&
        034.34,038.40,038.80,040.12,040.60,040.70,040.84,041.10,&
        041.62,042.34,042.98,042.88,041.50,038.76,035.18,031.04,&
        027.14,024.50,024.78,027.70,032.14,036.06,036.50,037.88,&
        038.50,038.68,038.84,039.10,039.56,040.34,040.88,040.82,&
        039.40,036.76,033.18,029.12,025.14,022.50,022.88,025.90,&
        029.96,033.86,034.30,035.68,036.42,036.68,036.84,037.10,&
        037.56,038.24,038.88,038.72,037.40,034.76,031.18,027.12,&
        023.14,020.60,020.98,023.90,027.88,031.66,032.10,033.58,&
        034.32,034.68,034.84,035.10,035.56,036.24,036.78,036.62,&
        035.30,032.72,029.18,025.14,021.24,018.70,019.08,021.90,&
        025.88,029.42,029.90,031.48,032.32,032.68,032.84,033.10,&
        033.56,034.22,034.68,034.42,033.20,030.72,027.28,023.22,&
        019.34,016.80,017.24,020.00,023.78,027.32,027.70,029.38,&
        030.24,030.68,030.94,031.20,031.66,032.22,032.58,032.32,&
        031.10,028.62,025.28,021.32,017.48,015.00,015.38,018.18,&
        021.80,025.22,025.70,027.28,028.24,028.78,029.04,029.30,&
        029.66,030.22,030.50,030.22,029.00,026.62,023.30,019.42,&
        015.64,013.10,013.54,016.28,019.80,023.12,023.60,025.24,&
        026.24,026.78,027.14,027.40,027.76,028.22,028.40,028.12,&
        026.80,024.52,021.30,017.52,013.78,011.30,011.74,014.48,&
        017.90,021.12,021.60,023.24,024.34,024.88,025.24,025.50,&
        025.86,026.22,026.40,025.98,024.70,022.48,019.40,015.72,&
        012.04,009.50,009.94,012.58,016.02,019.12,019.60,021.24,&
        022.34,022.98,023.34,023.70,024.00,024.30,024.40,023.88,&
        022.60,020.48,017.52,014.00,010.34,007.80,008.18,010.88,&
        014.22,017.18,017.60,019.34,020.44,021.16,021.54,021.90,&
        022.16,022.40,022.32,021.78,020.60,018.48,015.62,012.20,&
        008.68,006.00,006.44,009.18,012.42,015.28,015.80,017.44,&
        018.54,019.26,019.74,020.10,020.30,020.50,020.32,019.72,&
        018.50,016.54,013.84,010.68,007.14,004.40,004.74,007.58,&
        010.74,013.48,014.00,015.54,016.74,017.46,017.94,018.30,&
        018.50,018.58,018.32,017.72,016.50,014.64,012.24,009.18,&
        005.84,002.90,003.30,006.16,009.14,011.84,012.30,013.78,&
        014.94,015.66,016.24,016.50,016.70,016.70,016.42,005.78,&
        014.60,012.90,010.66,007.86,004.88,001.60,001.72,004.96,&
        007.84,010.24,010.70,012.14,013.24,013.96,014.44,014.80,&
        014.90,014.88,014.52,013.92,012.80,011.30,009.28,006.94,&
        004.32,001.80,001.94,004.34,006.78,008.94,009.40,010.58,&
        011.64,012.36,012.74,013.10,013.20,013.08,012.72,012.12,&
        011.10,009.86,008.30,006.50,004.60,003.10,003.16,004.50,&
        006.20,007.90,008.40,009.42,010.14,010.76,011.14,011.40,&
        011.40,011.38,011.02,010.46,009.70,008.72,007.64,006.46,&
        005.42,004.60,004.70,005.34,006.24,007.36,007.90,008.46,&
        008.92,009.28,009.54,009.70,009.70,009.68,009.42,009.06,&
        008.60,008.08,007.56,007.02,006.56,006.30,006.30,006.52,&
        006.96,007.38,008.15,008.15,008.15,008.15,008.15,008.15,&
        008.15,008.15,008.15,008.15,008.15,008.15,008.15,008.15,&
        008.15,008.15,008.15,008.15,008.15,008.15/
    !     Data Input      
    rlan = rga
    rlo = rgo      
    
    !     From "normal" geographic latitude 
    !     to angle from South Pole.       
    rla = rlan + 90
    IF (rlo  ==  360) THEN
        rlo = 0
    END IF
    !     PROXIMITY
    !     coefficients of the latitudinal points        
    LA1 = (INT(rla/2)+1)
    LA2 = LA1 + 1
    if(la2 > 91) la2=91
    !     coefficients of the longitudinal points        
    LO1 = (INT(rlo/18)+1)
    !orr      LO2 = LO1 + 1
    LO2 = MOD(LO1,20) + 1 
    !     Four points of Geomagnetic Coordinates
    gm1 = CORMAG(LO1,LA1)
    gm2 = CORMAG(LO1,LA2) 
    gm3 = CORMAG(LO2,LA1)
    gm4 = CORMAG(LO2,LA2)
    !     latitudinal points        
    !      X1 = ABS(rla - (INT(rla)))                        
    !      X2 = 2. - X1
    x = (rla/2.0 - (INT(rla/2.0)))
    !     longitudinal points        
    !      Y1 = ABS(rlo - (INT(rlo)))                        
    !      Y2 = 18. - Y1
    y =(rlo/18.0 - (INT(rlo/18.0))) 
    !     X AND Y VALUES
    !      x = X1 / (X1 + X2)
    !      y = Y1 / (Y1 + Y2)
    !     INTERPOLATION
    gmla = gm1 * (1 - x) * (1 - y) + gm2 * (1 - y) * (x) + gm3 * (y)&
        * (1 - x) + gm4 * (x) * (y)
    !     OUTPUT OF THE PROGRAM
    !     From corrected geomagnetic latitude from North Pole
    !     to "normal"  geomagnetic latitude.       
    rgma = 90. - gmla
END
!
!
SUBROUTINE STORM(ap,rga,rgo,coor,rgma,ut,doy,cf)
    !----------------------------------------------------------------------
    !      Fortran code to obtain the foF2 storm-time correction factor at 
    !      a given location and time, using the current and the 12 previous
    !      ap values as input.
    !
    !      ap ---> (13 elements integer array). Array with the preceeding
    !              13 value of the 3-hourly ap index. The 13th value
    !              in the array will contain the ap at the UT of interest,
    !              the 12th value will contain the 1st three hourly interval
    !              preceeding the time of interest, and so on to the 1st
    !              ap value at the earliest time.
    !     coor --> (integer). If coor = 2, rga should contain the 
    !                         geomagnetic latitude.
    !                         If coor = 1, rga should contain the 
    !                         geographic latitude.
    !     rga ---> (real, -90 to 90) geographic or geomagnetic latitude.
    !     rgo ---> (real, 0 to 360, positive east from Greenwich.)
    !                           geographic longitude, only used if coor=1.
    !     ut  ---> (integer, hours 00 to 23) Universal Time of interest.
    !     doy ---> (integer, 1 to 366)Day of the year.
    !     cf  ---> (real) The output; the storm-time correction factor used
    !              to scale foF2, foF2 * cf.
    !     rgma --> corrected magnetic latitude calculated from rga and rgo
    !
    !     This model and computer code was developed by E. Araujo-Pradere,
    !     T. Fuller-Rowell and M. Condrescu, SEC, NOAA, Boulder, USA
    !     Ref: 
    !     T. Fuller-Rowell, E. Araujo-Pradere, and M. Condrescu, An 
    !       Empirical Ionospheric Storm-Time Ionospheric Correction Model,
    !       Adv. Space Res. 8, 8, 15-24, 2000.
    !----------------------------------------------------------------------
    !     DIMENSIONS AND COEFFICIENTS VALUES
    DIMENSION c4(20)
    DATA c4/0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,&
        0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,&
        0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00/
    DIMENSION c3(20)
    DATA c3/0.00E+00,0.00E+00,0.00E+00,0.00E+00,0.00E+00,-9.44E-12,&
        0.00E+00,3.04E-12,0.00E+00,9.32E-12,-1.07E-11,0.00E+00,0.00E+00,&
        0.00E+00,1.09E-11,0.00E+00,0.00E+00,0.00E+00,0.00E+00,-1.01E-11/
    DIMENSION c2(20)
    DATA c2/1.16E-08,0.00E+00,0.00E+00,-1.46E-08,0.00E+00,9.86E-08,&
        2.25E-08,-1.67E-08,-1.62E-08,-9.42E-08,1.17E-07,4.32E-08,3.97E-08,&
        3.13E-08,-8.04E-08,3.91E-08,2.58E-08,3.45E-08,4.76E-08,1.13E-07/
    DIMENSION c1(20)
    DATA c1/-9.17E-05,-1.37E-05,0.00E+00,7.14E-05,0.00E+00,-3.21E-04,&
        -1.66E-04,-4.10E-05,1.36E-04,2.29E-04,-3.89E-04,-3.08E-04,&
        -2.81E-04,-1.90E-04,4.76E-05,-2.80E-04,-2.07E-04,-2.91E-04,&
        -3.30E-04,-4.04E-04/
    DIMENSION c0(20)
    DATA c0/1.0136E+00,1.0478E+00,1.00E+00,1.0258E+00,1.00E+00,&
        1.077E+00,1.0543E+00,1.0103E+00,9.9927E-01,9.6876E-01,1.0971E+00,&
        1.0971E+00,1.0777E+00,1.1134E+00,1.0237E+00,1.0703E+00,1.0248E+00,&
        1.0945E+00,1.1622E+00,1.1393E+00/
    DIMENSION fap(36)
    DATA fap/0.,0.,0.037037037,0.074074074,0.111111111,0.148148148,&
        0.185185185,0.222222222,0.259259259,0.296296296,0.333333333,&
        0.37037037,0.407407407,0.444444444,0.481481481,0.518518519,&
        0.555555556,0.592592593,0.62962963,0.666666667,0.703703704,&
        0.740740741,0.777777778,0.814814815,0.851851852,0.888888889,&
        0.925925926,0.962962963,1.,0.66666667,0.33333334,0.,0.333333,&
        0.666666,1.,0.7/
    integer code(8,6)
    data code/3,4,5,4,3,2,1,2,3,2,1,2,3,4,5,4,8,7,6,7,8,9,10,9,&
        13,12,11,12,13,14,15,14,18,17,16,17,18,19,20,19,18,17,16,17,&
        18,19,20,19/
    INTEGER ape(39)
    INTEGER ap(13)
    INTEGER ut,doy,dayno,coor,s1,s2,l1,l2
    REAL rgma, rap, rga, rgo, rs, rl
    !      CALLING THE PROGRAM TO CONVERT TO GEOMAGNETIC COORDINATES
    IF (coor  ==  1) THEN
        CALL CONVER (rga,rgo,rgma)
    ELSE IF (coor  ==  2) THEN
        rgma = rga
    ELSE
        WRITE (6,*)' '
        WRITE (6,*)' '
        WRITE (6,*)'   Wrong Coordinates Selection -------- >>', coor
        WRITE (6,*)' '
        GOTO 100
    ENDIF
    ! FROM 3-HOURLY TO HOURLY ap (New, interpolates between the three hourly ap values)
    ape(1)=ap(1)
    ape(2)=ap(1)
    ape(38)=ap(13)
    ape(39)=ap(13)
    DO k = 1,13
        i = (k * 3) - 1
        ape(i) = ap(k)
    END DO
    DO k = 1,12
        i = k * 3
        ape(i) = (ap(k)*2 + ap(k+1))/3.0
    END DO
    DO k = 2,13
        i = (k * 3) - 2
        ape(i) = (ap(k-1) + ap(k)*2)/3.0
    END DO
    !     FROM 3-HOURLY TO HOURLY ap (old version without interpolation)
    !      i = 1
    !      DO 10 k = 1,13
    !         DO j = 1,3
    !            ape(i) = ap(k)
    !            i = i + 1
    !            END DO
    !10    CONTINUE
    !     TO OBTAIN THE INTEGRAL OF ap.
    !     INTEGRAL OF ap
    if(ut == 24) ut=0
    IF (ut  ==  0 .OR. ut  ==  3 .OR. ut  ==  6 .OR. ut  ==  9 .OR.&
        ut  ==  12 .OR. ut  ==  15 .OR. ut  ==  18 .OR. ut  ==  21) THEN
        k = 1
    ELSE IF (ut  ==  1 .OR. ut  ==  4 .OR. ut  ==  7 .OR. ut  ==  10&
            .OR.ut  ==  13 .OR. ut  ==  16 .OR. ut  ==  19 .OR. ut  ==  22)&
            THEN
        k = 2
    ELSE IF (ut  ==  2 .OR. ut  ==  5 .OR. ut  ==  8 .OR. ut  ==  11&
            .OR. ut  ==  14 .OR. ut  ==  17 .OR. ut  ==  20 .OR. ut  ==  23)&
            THEN
        k = 3
    ELSE
        WRITE (6,*)' '
        WRITE (6,*)' '
        WRITE (6,*)'  Wrong Universal Time value -------- >>', ut
        WRITE (6,*)' '
        GOTO 100
    END IF
    rap = 0
    DO j = 1,36
        rap = rap + fap(j) * ape(k+j)
    END DO
    if(rap <= 200.)then
        cf=1.0
        goto 100
    end if
    if(doy > 366.or.doy < 1)then
        WRITE (6,*)' '
        WRITE (6,*)' '
        WRITE (6,*)' '
        WRITE (6,*)'      Wrong Day of Year value --- >>', doy
        WRITE (6,*)' '
        GOTO 100
    end if
    if(rgma > 90.0.or.rgma < -90.0)then
        WRITE (6,*)' '
        WRITE (6,*)' '
        WRITE (6,*)' '
        WRITE (6,*)'   Wrong GEOMAGNETIC LATITUDE value --- >>', rgma
        WRITE (6,*)' '
        GOTO 100
    end if
    !      write(6,*)rgma
    dayno=doy
    if(rgma < 0.0)then
        dayno=doy+172
        if(dayno > 365)dayno=dayno-365
    end if
    if (dayno >= 82) rs=(dayno-82.)/45.6+1.
    if (dayno < 82) rs=(dayno+283.)/45.6+1.
    s1=rs
    facs=rs-s1
    s2=s1+1
    if(s2 == 9) s2=1
    !      write(6,*)s1,s2,rs
    rgma = abs(rgma)
    rl=(rgma+10.)/20.+1
    if(rl == 6.0)rl=5.9
    l1=rl
    facl=rl-l1
    l2=l1+1
    !      write(6,*)l1,l2,rl
    !     FACTORS CALCULATIONS
    if(rap < 300.)then
        rapf=300.
        n1=code(s1,l1)
        cf1=c4(n1)*(rapf**4)+c3(n1) * (rapf**3) + c2(n1) * (rapf**2) +&
            c1(n1) * rapf + c0(n1)
        n2=code(s1,l2)
        cf2=c4(n2)*(rapf**4)+c3(n2) * (rapf**3) + c2(n2) * (rapf**2) +&
            c1(n2) * rapf + c0(n2)
        n3=code(s2,l1)
        cf3=c4(n3)*(rapf**4)+c3(n3) * (rapf**3) + c2(n3) * (rapf**2) +&
            c1(n3) * rapf + c0(n3)
        n4=code(s2,l2)
        cf4=c4(n4)*(rapf**4)+c3(n4) * (rapf**3) + c2(n4) * (rapf**2) +&
            c1(n4) * rapf + c0(n4)
        !     INTERPOLATION
        cf300=cf1*(1 - facs) * (1 - facl) + cf2 * (1 - facs) * (facl) +&
            cf3 * (facs) * (1 - facl) + cf4 * (facs) * (facl)
        cf = (cf300-1.0)*rap/100.-2.*cf300+3.
        goto 100
    end if
    n1=code(s1,l1)
    !      write(6,*)n1
    cf1 = c4(n1) * (rap**4) + c3(n1) * (rap**3) + c2(n1) * (rap**2) +&
        c1(n1) * rap + c0(n1)
    n2=code(s1,l2)
    cf2 = c4(n2) * (rap**4) + c3(n2) * (rap**3) + c2(n2) * (rap**2) +&
        c1(n2) * rap + c0(n2)
    n3=code(s2,l1)
    cf3 = c4(n3) * (rap**4) + c3(n3) * (rap**3) + c2(n3) * (rap**2) +&
        c1(n3) * rap + c0(n3)
    n4=code(s2,l2)
    cf4 = c4(n4) * (rap**4) + c3(n4) * (rap**3) + c2(n4) * (rap**2) +&
        c1(n4) * rap + c0(n4)
    !     INTERPOLATION
    cf = cf1 * (1 - facs) * (1 - facl) + cf2 * (1 - facs) * (facl) +&
        cf3 * (facs) * (1 - facl) + cf4 * (facs) * (facl)
    100   CONTINUE
    RETURN
END
!
!
FUNCTION STORME_AP(JDOY,XMLAT,AP)
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !
    ! EMPIRICAL STORM-E MODEL: COMPUTES A STORM-TO-QUIET RATIO (SQR) FACTOR
    ! TO ADJUST THE QUIESCENT E-REGION PEAK ELECTRON DENSITY TO ACCOUNT FOR
    ! ENHANCEMENTS DUE TO GEOMAGNETIC ACTIVITY. THE SQR FACTORS WERE
    ! COMPUTED FROM NO+ 4.3 UM VOLUME EMISSION RATES DERIVED FROM 
    ! TIMED/SABER LIMB RADIANCE MEASUREMENTS. THE SABER-DERIVED SQR FACTORS
    ! WERE FIT TO POWE-LAW IN THE ap INDEX.   
    !
    ! INPUT PARAMETERS:
    !
    !  JDOY      --- DAY OF YEAR (1-365) 
    !  XMLAT     --- MAGNETIC LATITUDE (DEGREES)
    !  AP        --- ap INDEX
    !
    ! OUTPUT PARAMETER
    ! 
    !  STORME_AP --- STORM-TO-QUIET RATIO (SQR) TO ADJUST QUIESCENT E-REGION
    !                PEAK ELECTRON DENSITY TO ACCOUNT FOR GEOMAGNETIC
    !                ENHANCEMENTS. SQR COMPUTED FROM A POWER-LAW FIT
    !                IN AP-INDEX: SQR=C1*AP**C2+C3
    !
    ! REFERENCES:
    ! 
    !  (1) Mertens et al. [submitted to JASR, 2011]
    !  (2) Fernandez et al. [JASR, Vol. 46, 2010]
    !  (3) Mertens et al. [Proc. of SPIE, Vol. 7475, 2009]
    !  (4) Mertens et al. [Proc. of SPIE, Vol. 6745, 2007]
    !  (5) Mertens et al. [JASR, Vol. 39, 2007] 
    ! 
    ! SOFTWARE WRITTEN BY Christopher J. Mertens
    !                     NASA Langley Research Center
    !                     Atmospheric Sciences Competency
    !                     21 Langley Blvd., Mail Stop 401B
    !                     Hampton, VA 23681-2199
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    PARAMETER(NMLG=37,NDBD=5)
    DIMENSION C3(NMLG,NDBD)
    DIMENSION XMLG(NMLG),IDBD(NDBD),C1(NMLG,NDBD),C2(NMLG,NDBD)
    LOGICAL mess
    COMMON /iounit/konsol,mess
    DATA XMLG/-90.0,-85.0,-80.0,-75.0,-70.0,-65.0,-60.0,-55.0,-50.0,&
        -45.0,-40.0,-35.0,-30.0,-25.0,-20.0,-15.0,-10.0,-5.0,&
        0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,&
        55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0/
    DATA IDBD/79,171,264,354,366/
    !
    ! JDOY        | 0-79 ||80-171||172-264||265-354||355-365|
    !
    DATA ((C1(IML,IDB),IDB=1,NDBD),IML=1,NMLG)&
        / 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, & !-90.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, & !-85.0
        0.00000, 0.00508, 0.17360, 0.00000, 0.00000, &!-80.0
        0.00000, 0.31576, 0.31498, 0.00000, 0.00000, &!-75.0
        0.00000, 0.39217, 0.40121, 0.00000, 0.00000, &!-70.0
        0.00000, 0.32634, 0.30179, 0.00000, 0.00000, &!-65.0
        0.11573, 0.06211, 0.16230, 0.20233, 0.11573, &!-60.0
        0.00526, 0.00013, 0.00204, 0.04965, 0.00526, &!-55.0
        0.00011, 0.00013, 0.00018, 0.00040, 0.00011, &!-50.0
        0.00001, 0.00002, 0.00040, 0.00001, 0.00001, &!-45.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-40.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-35.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-30.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-25.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-20.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-15.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-10.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! -5.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!  0.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!  5.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 10.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 15.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 20.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 25.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 30.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 35.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 40.0
        0.00000,-0.02738, 0.00004, 0.00001, 0.00000, &! 45.0
        0.00001,-0.00022, 0.00001, 0.00004, 0.00001, &! 50.0
        0.00126, 0.00062, 0.00011, 0.00345, 0.00126, &! 55.0
        0.14923, 0.05483, 0.07113, 0.18282, 0.14923, &! 60.0
        0.37361, 0.00000, 0.00000, 0.44592, 0.37361, &! 65.0
        0.27792, 0.00000, 0.00000, 0.00804, 0.27792, &! 70.0
        0.06445, 0.00000, 0.00000, 0.10315, 0.06445, &! 75.0
        0.00149, 0.00000, 0.00000, 0.00073, 0.00149, &! 80.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 85.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000/ ! 90.0  
    !
    ! JDOY        | 0-79 ||80-171||172-264||265-354||355-365|
    !
    DATA ((C2(IML,IDB),IDB=1,NDBD),IML=1,NMLG)&
        / 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-90.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-85.0
        0.00000, 1.00000, 0.32415, 0.00000, 0.00000, &!-80.0
        0.00000, 0.36538, 0.35455, 0.00000, 0.00000, &!-75.0
        0.00000, 0.41287, 0.38062, 0.00000, 0.00000, &!-70.0
        0.00000, 0.52224, 0.52810, 0.00000, 0.00000, &!-65.0
        0.73025, 0.90723, 0.68107, 0.64815, 0.73025, &!-60.0
        1.29410, 2.06038, 1.47332, 0.84843, 1.29410, &!-55.0
        1.79442, 1.77511, 1.59906, 1.59141, 1.79442, &!-50.0
        1.84434, 1.70607, 1.03056, 1.92168, 1.84434, &!-45.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-40.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-35.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-30.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-25.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-20.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-15.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!-10.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! -5.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!  0.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &!  5.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 10.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 15.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 20.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 25.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 30.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 35.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 40.0
        2.04797, 0.27034, 1.00000, 1.90792, 2.04797, &! 45.0
        2.24520, 1.00000, 1.88721, 2.01535, 2.24520, &! 50.0
        1.56493, 1.52468, 1.93389, 1.38532, 1.56493, &! 55.0
        0.71442, 0.87492, 0.78890, 0.66828, 0.71442, &! 60.0
        0.53546, 0.00000, 0.00000, 0.42597, 0.53546, &! 65.0
        0.48647, 0.00000, 0.00000, 1.00000, 0.48647, &! 70.0
        0.67340, 0.00000, 0.00000, 0.36809, 0.67340, &! 75.0
        1.44025, 0.00000, 0.00000, 1.13529, 1.44025, &! 80.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000, &! 85.0
        0.00000, 0.00000, 0.00000, 0.00000, 0.00000/ ! 90.0 
    !
    ! JDOY        | 0-79 ||80-171||172-264||265-354||355-365|
    !
    DATA ((C3(IML,IDB),IDB=1,NDBD),IML=1,NMLG)&
        / 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &!-90.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &!-85.0
        1.00000, 1.03132, 0.76703, 1.00000, 1.00000, &!-80.0
        1.00000, 0.57588, 0.56324, 1.00000, 1.00000, &!-75.0
        1.00000, 0.41370, 0.38549, 1.00000, 1.00000, &!-70.0
        1.00000, 0.51704, 0.50217, 1.00000, 1.00000, &!-65.0
        0.55236, 0.80162, 0.60824, 0.46999, 0.55236, &!-60.0
        0.90923, 0.99688, 0.96752, 0.67312, 0.90923, &!-55.0
        0.99338, 0.98486, 0.99503, 0.87473, 0.99338, &!-50.0
        1.00031, 1.00369, 1.00225, 0.91242, 1.00031, &!-45.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &!-40.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &!-35.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &!-30.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &!-25.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &!-20.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &!-15.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &!-10.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &! -5.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &!  0.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &!  5.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &! 10.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &! 15.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &! 20.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &! 25.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &! 30.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &! 35.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &! 40.0
        1.04797, 1.02319, 0.98581, 1.01457, 1.04797, &! 45.0
        1.03332, 0.97016, 0.97807, 0.99044, 1.03332, &! 50.0
        1.00633, 0.94822, 0.96340, 0.95363, 1.00633, &! 55.0
        0.67902, 0.71540, 0.70230, 0.60821, 0.67902, &! 60.0
        0.35017, 1.00000, 1.00000, 0.51033, 0.35017, &! 65.0
        0.63358, 1.00000, 1.00000, 1.37782, 0.63358, &! 70.0
        0.85724, 1.00000, 1.00000, 0.91942, 0.85724, &! 75.0
        0.92703, 1.00000, 1.00000, 1.00502, 0.92703, &! 80.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, &! 85.0
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000/ ! 90.0
    !
    ! ... Find Season-Averaged Coefficient Index 
    !
    IDXS=0
    IF(JDOY <= IDBD(1)) IDXS=1
    DO IS=2,NDBD
        IF((JDOY > IDBD(IS-1)).AND.(JDOY <= IDBD(IS))) IDXS=IS
    ENDDO
    IF(IDXS == 0) THEN
        if(mess) WRITE(konsol,*) 'ERROR IN STORME_AP: ',&
            'PROBLEM FINDING SEASON-AVERAGED COEFFICIENT',&
            'DAY OF YEAR = ',JDOY
        STORME_AP=-5.0
        GOTO 222 
    ENDIF
    !
    ! ... Find Magnetic Latitude Coefficient Index
    ! 
    IDXL=0
    DELG=ABS(XMLG(1)-XMLG(2))
    DELD=DELG/2.0
    YMP=XMLG(1)+DELD
    YMM=XMLG(NMLG)-DELD
    IF((XMLAT >= XMLG(1)).AND.(XMLAT <= YMP)) IDXL=1
    IF((XMLAT > YMM).AND.(XMLAT <= XMLG(NMLG))) IDXL=NMLG
    DO IL=2,NMLG-1
        YMP=XMLG(IL)+DELD
        YMM=XMLG(IL)-DELD
        IF((XMLAT > YMM).AND.(XMLAT <= YMP)) IDXL=IL
    ENDDO
    IF(IDXL == 0) THEN
        if(mess) WRITE(konsol,*) 'ERROR IN STORME_AP: ', &
            'PROBLEM FINDING MAGNETIC LATITUDE COEFFICIENT',&
            'MAGNETIC LATITUDE(DEGREES) = ',XMLAT
        STORME_AP=-5.0
        GOTO 222
    ENDIF
    !
    ! ... COMPUTE E-REGION ELECTRON DENSITY GEOMAGNETIC STORM ENHANCEMET
    ! ... FACTOR (i.e., THE STORM-TO-QUIET RATIO (SQR)) 
    !
    SQR=C1(IDXL,IDXS)*AP**(C2(IDXL,IDXS))+C3(IDXL,IDXS)
    IF(SQR < 1.0) SQR=1.0
    STORME_AP=SQR
    
    222   RETURN
END    
!
!
!***************************************************************************
!
subroutine spreadf_brazil(idoy,idiy,f107,geolat,osfbr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
    !
    !       SUBROUTINE CALCULATES PERCENTAGE OF SPREAD F OCCURRENCE OVER 
    !       BRAZILIAN SECTOR AS DESCRIBED IN:
    !       ABDU ET AL., Advances in Space Research, 31(3), 
    !       703-716, 2003
    !
    !    INPUT:
    !         IDOY: DAY OF YEAR (1 TO 365/366)
    !         IDIY: DAYS IN YEAR (365 OR 366)
    !         F107: F10.7 cm SOLAR FLUX (DAILY VALUE)
    !         GEOLAT: BRAZILIAN GEOGRAPHIC LATITUDE BETWEEN -4 AND -22.5
    !
    !    OUTPUT:         
    !         OSFBR(25): PERCENTAGE OF SPREAD F OCCURRENCE FOR 25 TIME 
    !                    STEPS FROM LT=18 TO LT=7 ON THE NEXT DAY IN
    !                    STEPS OF 0.5 HOURS.
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    dimension param(3),osfbr(25),coef_sfa(684),coef_sfb(684),&
        sosf(2,32,3,12)
    common/mflux/kf,n     
    data coef_sfa/&
        0.07,0.13,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.05,0.04,0.03&
        ,0.06,0.07,0.02,0.03,0.03,0.07,0.06,0.07,0.21,0.28,0.34,0.16&
        ,0.12,0.00,0.02,0.02,0.04,0.05,0.02,0.11,0.19,0.31,0.31,0.11&
        ,0.14,0.16,0.03,0.00,0.00,0.02,0.00,0.00,0.05,0.55,0.61,0.28&
        ,0.26,0.10,0.15,0.23,0.07,0.06,0.03,0.03,0.41,0.88,0.89,0.65&
        ,0.19,0.18,0.17,0.10,0.14,0.15,0.03,0.14,0.46,0.72,0.71,0.53&
        ,0.57,0.38,0.30,0.14,0.00,0.04,0.03,0.02,0.21,0.84,0.87,0.72&
        ,0.79,0.60,0.65,0.70,0.29,0.19,0.19,0.32,0.73,0.96,0.99,0.84&
        ,0.75,0.78,0.79,0.70,0.63,0.24,0.28,0.53,0.75,0.77,0.75,0.85&
        ,0.78,0.51,0.59,0.24,0.00,0.07,0.05,0.06,0.33,0.92,0.96,0.89&
        ,0.90,0.84,0.86,0.81,0.33,0.27,0.23,0.47,0.90,1.00,1.00,0.96&
        ,0.96,0.89,0.92,0.84,0.80,0.27,0.35,0.61,0.81,0.93,0.86,0.97&
        ,0.84,0.65,0.75,0.25,0.00,0.04,0.08,0.06,0.53,0.93,0.96,0.94&
        ,0.95,0.84,0.91,0.71,0.18,0.17,0.21,0.42,0.92,0.99,0.97,0.92&
        ,0.92,0.93,0.92,0.67,0.58,0.21,0.38,0.55,0.83,0.90,0.89,0.97&
        ,0.84,0.71,0.91,0.21,0.02,0.07,0.03,0.03,0.60,0.95,0.96,0.92&
        ,0.97,0.91,0.92,0.67,0.11,0.08,0.09,0.23,0.90,0.99,0.99,0.96&
        ,0.96,0.93,0.98,0.63,0.25,0.08,0.12,0.41,0.79,0.95,0.98,0.99&
        ,0.86,0.80,0.94,0.22,0.02,0.04,0.03,0.03,0.63,0.95,0.96,0.94&
        ,0.98,0.90,0.91,0.59,0.10,0.04,0.07,0.15,0.83,0.97,0.97,0.90&
        ,0.92,0.93,0.95,0.57,0.12,0.03,0.05,0.23,0.74,0.94,0.94,0.99&
        ,0.84,0.84,0.90,0.24,0.02,0.07,0.07,0.03,0.60,0.95,0.96,0.97&
        ,0.93,0.82,0.83,0.51,0.08,0.07,0.09,0.09,0.71,0.95,0.92,0.87&
        ,0.91,0.91,0.89,0.50,0.14,0.03,0.06,0.14,0.61,0.84,0.89,0.94&
        ,0.77,0.82,0.84,0.34,0.10,0.11,0.12,0.06,0.43,0.87,0.94,0.97&
        ,0.91,0.77,0.68,0.42,0.06,0.08,0.10,0.04,0.51,0.78,0.71,0.77&
        ,0.85,0.88,0.77,0.35,0.16,0.05,0.08,0.15,0.53,0.70,0.60,0.89&
        ,0.85,0.71,0.72,0.26,0.16,0.17,0.08,0.15,0.38,0.73,0.91,0.91&
        ,0.89,0.68,0.53,0.26,0.06,0.12,0.08,0.09,0.32,0.63,0.67,0.77&
        ,0.81,0.79,0.59,0.21,0.14,0.03,0.06,0.09,0.23,0.51,0.34,0.79&
        ,0.88,0.66,0.59,0.16,0.18,0.15,0.16,0.16,0.33,0.67,0.75,0.88&
        ,0.80,0.64,0.52,0.16,0.04,0.09,0.04,0.09,0.24,0.47,0.53,0.50&
        ,0.73,0.69,0.48,0.11,0.14,0.03,0.03,0.03,0.20,0.37,0.28,0.54&
        ,0.81,0.64,0.49,0.18,0.12,0.17,0.16,0.19,0.31,0.57,0.70,0.83&
        ,0.76,0.57,0.52,0.13,0.04,0.06,0.05,0.08,0.21,0.49,0.47,0.39&
        ,0.69,0.66,0.43,0.11,0.10,0.02,0.00,0.03,0.16,0.39,0.24,0.35&
        ,0.77,0.45,0.39,0.10,0.10,0.13,0.15,0.18,0.29,0.57,0.70,0.69&
        ,0.71,0.49,0.54,0.20,0.05,0.06,0.05,0.06,0.27,0.42,0.36,0.42&
        ,0.61,0.59,0.50,0.08,0.06,0.02,0.03,0.02,0.16,0.40,0.17,0.31&
        ,0.68,0.30,0.28,0.13,0.10,0.16,0.14,0.08,0.19,0.50,0.63,0.62&
        ,0.63,0.45,0.51,0.13,0.06,0.07,0.04,0.06,0.27,0.42,0.28,0.35&
        ,0.68,0.53,0.57,0.15,0.05,0.00,0.00,0.05,0.31,0.33,0.18,0.22&
        ,0.59,0.32,0.21,0.06,0.10,0.16,0.12,0.10,0.19,0.41,0.55,0.54&
        ,0.69,0.43,0.43,0.15,0.06,0.05,0.05,0.08,0.29,0.39,0.23,0.29&
        ,0.57,0.51,0.56,0.13,0.06,0.00,0.00,0.05,0.34,0.27,0.19,0.24&
        ,0.49,0.16,0.13,0.09,0.04,0.11,0.11,0.05,0.17,0.32,0.49,0.49&
        ,0.60,0.42,0.38,0.11,0.06,0.04,0.07,0.07,0.25,0.36,0.21,0.25&
        ,0.65,0.48,0.53,0.17,0.05,0.00,0.00,0.11,0.29,0.14,0.20,0.22&
        ,0.44,0.16,0.18,0.07,0.04,0.04,0.07,0.03,0.12,0.23,0.39,0.43&
        ,0.57,0.40,0.35,0.14,0.06,0.03,0.04,0.07,0.18,0.27,0.14,0.15&
        ,0.45,0.50,0.50,0.19,0.06,0.00,0.02,0.05,0.26,0.19,0.15,0.18&
        ,0.23,0.09,0.12,0.06,0.04,0.02,0.02,0.02,0.10,0.03,0.14,0.26&
        ,0.39,0.34,0.22,0.07,0.03,0.00,0.04,0.01,0.15,0.01,0.04,0.14&
        ,0.41,0.39,0.35,0.13,0.02,0.00,0.00,0.06,0.17,0.07,0.06,0.14&
        ,0.07,0.02,0.03,0.00,0.00,0.00,0.00,0.00,0.00,0.01,0.03,0.08&
        ,0.19,0.14,0.14,0.00,0.03,0.01,0.02,0.00,0.09,0.00,0.01,0.00&
        ,0.18,0.09,0.16,0.08,0.01,0.00,0.02,0.02,0.15,0.00,0.03,0.04/
    
    data coef_sfb/&
        0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00&
        ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.03,0.00,0.00,0.00,0.00&
        ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.02,0.01,0.00,0.00,0.00&
        ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00&
        ,0.00,0.01,0.00,0.00,0.00,0.00,0.00,0.01,0.01,0.00,0.00,0.00&
        ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.05,0.03,0.00,0.02,0.00&
        ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00&
        ,0.00,0.04,0.00,0.01,0.00,0.00,0.00,0.01,0.01,0.05,0.00,0.00&
        ,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.04,0.00,0.03,0.03,0.00&
        ,0.01,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.01,0.00&
        ,0.01,0.04,0.04,0.03,0.00,0.01,0.00,0.01,0.00,0.27,0.14,0.06&
        ,0.05,0.04,0.02,0.00,0.00,0.00,0.00,0.04,0.09,0.48,0.43,0.27&
        ,0.05,0.04,0.01,0.00,0.00,0.00,0.00,0.00,0.00,0.13,0.16,0.06&
        ,0.26,0.12,0.29,0.04,0.01,0.02,0.00,0.01,0.08,0.65,0.56,0.45&
        ,0.43,0.42,0.42,0.09,0.00,0.02,0.00,0.00,0.34,0.67,0.73,0.72&
        ,0.10,0.05,0.04,0.00,0.01,0.00,0.00,0.00,0.00,0.18,0.39,0.15&
        ,0.61,0.37,0.51,0.06,0.01,0.02,0.01,0.01,0.18,0.72,0.63,0.80&
        ,0.77,0.66,0.70,0.19,0.00,0.02,0.02,0.02,0.41,0.68,0.88,0.85&
        ,0.24,0.11,0.08,0.00,0.01,0.00,0.00,0.00,0.00,0.28,0.51,0.29&
        ,0.75,0.48,0.57,0.11,0.00,0.02,0.01,0.01,0.19,0.77,0.77,0.88&
        ,0.89,0.81,0.74,0.21,0.02,0.02,0.02,0.02,0.42,0.71,0.93,0.95&
        ,0.49,0.30,0.19,0.00,0.00,0.00,0.00,0.01,0.06,0.38,0.64,0.48&
        ,0.86,0.60,0.62,0.12,0.00,0.02,0.01,0.00,0.18,0.81,0.84,0.94&
        ,0.88,0.79,0.70,0.26,0.03,0.02,0.02,0.02,0.36,0.61,0.98,0.93&
        ,0.60,0.46,0.31,0.03,0.00,0.01,0.00,0.00,0.09,0.50,0.71,0.58&
        ,0.90,0.65,0.66,0.10,0.00,0.02,0.01,0.02,0.15,0.69,0.80,0.86&
        ,0.84,0.75,0.64,0.09,0.03,0.00,0.00,0.04,0.26,0.54,0.78,0.92&
        ,0.62,0.59,0.44,0.01,0.00,0.01,0.00,0.00,0.13,0.52,0.77,0.63&
        ,0.84,0.67,0.63,0.11,0.00,0.00,0.03,0.03,0.18,0.65,0.75,0.84&
        ,0.81,0.63,0.47,0.06,0.02,0.00,0.00,0.05,0.14,0.49,0.76,0.91&
        ,0.58,0.63,0.47,0.09,0.00,0.07,0.01,0.04,0.15,0.48,0.68,0.61&
        ,0.79,0.63,0.55,0.12,0.01,0.01,0.02,0.05,0.13,0.57,0.51,0.63&
        ,0.72,0.54,0.43,0.11,0.02,0.00,0.00,0.09,0.16,0.39,0.59,0.72&
        ,0.46,0.55,0.39,0.07,0.01,0.03,0.03,0.06,0.15,0.37,0.51,0.50&
        ,0.61,0.43,0.38,0.11,0.01,0.03,0.02,0.03,0.10,0.38,0.38,0.60&
        ,0.58,0.42,0.38,0.15,0.02,0.00,0.00,0.11,0.13,0.24,0.41,0.51&
        ,0.36,0.36,0.21,0.04,0.04,0.03,0.06,0.05,0.06,0.26,0.39,0.43&
        ,0.43,0.31,0.24,0.09,0.02,0.00,0.02,0.02,0.06,0.24,0.24,0.40&
        ,0.53,0.19,0.28,0.13,0.02,0.02,0.02,0.09,0.13,0.17,0.24,0.40&
        ,0.32,0.27,0.17,0.03,0.04,0.02,0.04,0.03,0.06,0.13,0.34,0.36&
        ,0.42,0.31,0.20,0.09,0.03,0.00,0.02,0.01,0.07,0.19,0.24,0.32&
        ,0.44,0.10,0.23,0.13,0.03,0.02,0.00,0.09,0.12,0.17,0.21,0.33&
        ,0.32,0.23,0.16,0.00,0.02,0.04,0.03,0.03,0.06,0.15,0.29,0.34&
        ,0.36,0.26,0.28,0.07,0.01,0.00,0.01,0.02,0.04,0.19,0.17,0.27&
        ,0.34,0.14,0.26,0.09,0.03,0.02,0.00,0.06,0.13,0.09,0.16,0.22&
        ,0.29,0.21,0.15,0.00,0.02,0.02,0.02,0.03,0.11,0.16,0.26,0.28&
        ,0.29,0.22,0.27,0.05,0.01,0.00,0.01,0.01,0.02,0.14,0.09,0.19&
        ,0.25,0.19,0.25,0.07,0.02,0.02,0.00,0.00,0.09,0.07,0.12,0.15&
        ,0.23,0.20,0.16,0.00,0.03,0.04,0.00,0.00,0.08,0.09,0.21,0.18&
        ,0.22,0.21,0.19,0.02,0.02,0.00,0.01,0.03,0.04,0.08,0.06,0.14&
        ,0.20,0.12,0.23,0.02,0.00,0.02,0.00,0.00,0.05,0.05,0.09,0.11&
        ,0.14,0.16,0.13,0.00,0.03,0.04,0.00,0.00,0.05,0.05,0.04,0.09&
        ,0.09,0.13,0.16,0.03,0.01,0.00,0.01,0.03,0.01,0.03,0.04,0.10&
        ,0.14,0.09,0.17,0.02,0.02,0.00,0.00,0.02,0.04,0.04,0.03,0.07&
        ,0.00,0.11,0.09,0.00,0.02,0.00,0.00,0.00,0.01,0.00,0.02,0.02&
        ,0.02,0.06,0.11,0.00,0.00,0.00,0.00,0.01,0.00,0.00,0.01,0.02&
        ,0.06,0.09,0.13,0.00,0.02,0.00,0.03,0.02,0.03,0.01,0.02,0.01/
    
    param(1)=idoy
    param(2)=f107
    param(3)=geolat
    n=idiy-365
    
    if(param(1) <= 31.)kf=1
    if(param(1) > 31..and.param(1) <= (59+n))kf=2
    if(param(1) > (59+n).and.param(1) <= (90+n))kf=3
    if(param(1) > (90+n).and.param(1) <= (120+n))kf=4
    if(param(1) > (120+n).and.param(1) <= (151+n))kf=5
    if(param(1) > (151+n).and.param(1) <= (181+n))kf=6
    if(param(1) > (181+n).and.param(1) <= (212+n))kf=7
    if(param(1) > (212+n).and.param(1) <= (243+n))kf=8
    if(param(1) > (243+n).and.param(1) <= (273+n))kf=9
    if(param(1) > (273+n).and.param(1) <= (304+n))kf=10
    if(param(1) > (304+n).and.param(1) <= (334+n))kf=11
    if(param(1) > (334+n).and.param(1) <= (365+n))kf=12
    
    do i=1,32
        do j=1,3
            do k=1,12
                sosf(1,i,j,k)=0.
                sosf(2,i,j,k)=0.           
            enddo
        enddo
    enddo
             
    kc=0
    do i=5,23
        do j=1,3
            do k=1,12
                kc=kc+1
                sosf(1,i,j,k)=coef_sfa(kc)
                sosf(2,i,j,k)=coef_sfb(kc) 
            enddo
        enddo
    enddo
    
    kk=0    
    do it=1600,3200,50      
        slt=it/100.
        osft=0.
        do i=1,23
            il=i+3
            if(il > 23)il=il-23
            do j=1,12
                jl=j+2
                if(jl > 12)jl=jl-12
                do m=1,3
                    ml=m+1
                    if(ml > 3)ml=ml-3
                    do  l=1,2
                        bspl4=bspl4t(i,slt)*bspl2s(j,param(1))*&
                            bspl2l(l,param(3))*bspl2f(m,param(2))
                        osft=osft+bspl4*sosf(l,il,ml,jl)           
                    enddo
                enddo
            enddo
        enddo
        if(slt > 17.98.and.slt < 30.01)then
            kk=kk+1
            osfbr(kk)=osft 
        endif
    enddo
   
    
    do iii=1,25 
        if(osfbr(iii) > 1.) osfbr(iii)=1.
        if(osfbr(iii) < 0.) osfbr(iii)=0.
    enddo
    return
end

function bspl4t(i,t1)
    dimension tt(0:78),b(30,30)
    data tt/16.00,16.50,17.00,17.50,18.00,18.50,19.00,19.50,20.00,&
        20.50,21.00,22.00,23.00,24.00,25.00,26.00,27.00,27.50,28.00,&
        28.50,29.00,29.50,30.00,30.50,31.00,32.00,40.00,40.50,41.00,&
        41.50,42.00,42.50,43.00,43.50,44.00,44.50,45.00,46.00,47.00,&
        48.00,49.00,50.00,51.00,51.50,52.00,52.50,53.00,53.50,54.00,&
        54.50,55.00,56.00,64.00,64.50,65.00,65.50,66.00,66.50,67.00,&
        67.50,68.00,68.50,69.00,70.00,71.00,72.00,73.00,74.00,75.00,&
        75.50,76.00,76.50,77.00,77.50,78.00,78.50,79.00,80.00,88.00/
    
    t=t1
    if(i >= 0.and.t < tt(i)) then
        t=t+24.
    endif
    do j=i,i+4-1
        if(t >= tt(j).and.t < tt(j+1)) then
            b(j,1)=1.
        else
            b(j,1)=0.
        endif
    enddo
    do j=2,4
        do k=i,i+4-j
            b(k,j)=(t-tt(k))/(tt(k+j-1)-tt(k))*b(k,j-1)
            b(k,j)=b(k,j)+(tt(k+j)-t)/(tt(k+j)-tt(k+1))*b(k+1,j-1)
        enddo
    enddo
    
    bspl4t=b(i,4)
    
    return
    
end

function bspl2s(i,t1)
    dimension ts(0:36),b(30,30)
    
    data ts/ 15,46,74,105,135,166,196,227,258,288,319,349,&
        380,411,439,470,500,531,561,592,623,653,684,714,&
        745,776,804,835,865,896,926,957,988,1018,1049,&
        1079,1110/
    
    t=t1
    if(i >= 0.and.t < ts(i)) then
        t=t+365.
    endif
    do j=i,i+2-1
        if(t >= ts(j).and.t < ts(j+1)) then
            b(j,1)=1.
        else
            b(j,1)=0.
        endif
    enddo
    
    do j=2,2
        do k=i,i+4-j
            b(k,j)=(t-ts(k))/(ts(k+j-1)-ts(k))*b(k,j-1)
            b(k,j)=b(k,j)+(ts(k+j)-t)/(ts(k+j)-ts(k+1))*b(k+1,j-1)
        enddo
    enddo
    
    bspl2s=b(i,2)
    return
end

function bspl2l(i,t1)
    dimension ts(0:6),b(30,30)
    
    data ts/ 94.,112.5,454.,472.5,814.,832.5,1174./
    
    t=t1
    if(i >= 0.and.t < ts(i)) then
        t=t+360.
    endif
    do j=i,i+2-1
        if(t >= ts(j).and.t < ts(j+1)) then
            b(j,1)=1.
        else
            b(j,1)=0.
        endif
    enddo
    
    do j=2,2
        do k=i,i+2-j
            b(k,j)=(t-ts(k))/(ts(k+j-1)-ts(k))*b(k,j-1)
            b(k,j)=b(k,j)+(ts(k+j)-t)/(ts(k+j)-ts(k+1))*b(k+1,j-1)
        enddo
    enddo
    
    bspl2l=b(i,2)
    
    return
    
end

function bspl2f(i,t1)
    dimension ts(0:9),b(30,30),&
        ifnodes1(12),ifnodes2(12),ifnodes3(12)
    common/mflux/kf,n
    
    data ifnodes1 / 78, 77, 75, 79, 80, 77, 78, 80, 76, 81, 78, 78/
    data ifnodes2 /144,140,139,142,139,146,142,139,150,151,150,157/
    data ifnodes3 /214,211,201,208,213,220,203,209,213,215,236,221/ 
    
    ts(0)=ifnodes1(kf)
    ts(1)=ifnodes2(kf)
    ts(2)=ifnodes3(kf)
    ts(3)=ts(1)+367
    ts(4)=ts(2)+367
    ts(5)=ts(3)+367
    ts(6)=ts(4)+367
    ts(7)=ts(5)+367
    ts(8)=ts(6)+367
    ts(9)=ts(7)+367
    
    t=t1
    if(i >= 0.and.t < ts(i)) then
        t=t+367.
    endif
    do j=i,i+2-1
        if(t >= ts(j).and.t < ts(j+1)) then
            b(j,1)=1.
        else
            b(j,1)=0.
        endif
    enddo
    
    do j=2,2
        do k=i,i+2-j
            b(k,j)=(t-ts(k))/(ts(k+j-1)-ts(k))*b(k,j-1)
            b(k,j)=b(k,j)+(ts(k+j)-t)/(ts(k+j)-ts(k+1))*b(k+1,j-1)
        enddo
    enddo
    
    bspl2f=b(i,2)
    return
end
!
!
function ckp(ap)
    !-----------------------------------------------------------------------
    ! Converts ap index (ap is integer variable varying from 0 to 400) into 
    ! kp index (xkp is real variable varying from 0 to 9). Using standard
    ! tables for deriving the 3-hourly ap index from the 3-hourly Kp index
    ! (e.g., http://www.ngdc.noaa.gov/stp/GEOMAG/kp_ap.shtml) 
    !-----------------------------------------------------------------------
    integer        ap,ap_array
    real        kp_array,ap_log_array
    dimension     ap_array(28),kp_array(28),alap(28)
    data ap_array /0,2,3,4,5,6,7,9,12,15,18,22,27,32,39,48,56,67,&
        80,94,111,132,154,179,207,236,300,400/
    
    do i=2,28
        kp_array(i)=(i-1)/3.
    end do
    
    if(ap == 0) then
        ckp=0.0
        return
    endif
    if(ap == 1) then
        ckp=kp_array(2)/2.
        return
    endif
    if(ap < 8.and.ap > 1) then
        ckp=kp_array(ap)
        return
    endif
    xl_ap=log(ap*1.0)
    
    i=8
    1257    alap(i)=log(ap_array(i)*1.0)
    if(xl_ap > alap(i)) then
        i=i+1
        if(i <= 28) goto 1257
    endif
    slope=(kp_array(i)-kp_array(i-1))/(alap(i)-alap(i-1))
    
    ckp = kp_array(i) + slope * (xl_ap - alap(i))
    
    return
end
!                   
!        
subroutine auroral_boundary(xkp,xmlt,cgmlat,ab_mlat)
    !-----------------------------------------------------------------------
    ! Computes equatorward auroral boundary values for givern kp value.
    ! kp given in units of 0.1 (xkp) for the range from 0.0 to 9.0. Model 
    ! values are only used for kp=0,1,2,3,4,5,6,7,8,9 and a linear inter-
    ! polation is applied for intermediate kp values.
    ! 
    ! The auroral oval boundary is given as an array for corrected magnetic 
    ! latitude CGM (ab_mlat). The 48 values correspond to the MLT values 
    ! of 0.0,0.5,1.0,1.5,2.0 .. 23.5. If the input xmlt is greater than
    ! -1, then the program determines the CGM latitude, cgmlat, that 
    ! corresponds to the given MLT value (xmlt).
    !
    ! Y. Zhang and L.J. Paxton, An empirical Kp-dependent global auroral 
    ! model based on TIMED/GUVI FUV data, Journal of Atmospheric and 
    ! Solar-Terrestrial Physics 70, 1231–1242, 2008.
    ! 
    !----------------------------------------------------------------------- 
    dimension zp_mlat(48,10),ab_mlat(48),ab_mlt(48)
    data zp_mlat/&
        66.1,65.9,65.9,66.0,66.2,66.5,66.7,67.0,67.4,67.8,68.2,68.7,&
        69.4,70.0,70.5,70.8,71.0,71.3,71.9,72.7,73.8,75.4,77.0,77.8,&
        77.3,76.7,76.4,76.2,75.9,75.4,74.7,74.1,73.5,73.0,72.5,71.9,&
        71.6,71.1,70.6,70.0,69.2,68.5,67.9,67.5,67.2,66.9,66.7,66.4,&
        63.0,62.9,63.0,63.0,63.1,63.2,63.3,63.6,64.0,64.5,65.0,65.6,&
        66.2,66.7,67.1,67.6,68.1,68.6,69.4,70.3,71.5,72.8,74.1,74.9,&
        74.8,74.7,74.7,74.5,74.0,73.1,72.2,71.2,70.4,69.7,68.9,68.2,&
        67.5,66.9,66.3,65.8,65.4,65.0,64.5,64.1,63.7,63.3,63.1,62.9,&
        61.1,61.3,61.5,61.7,61.9,62.1,62.2,62.4,62.7,63.0,63.4,64.0,&
        64.4,65.0,65.5,65.9,66.4,66.9,67.5,68.3,69.3,70.5,71.6,72.3,&
        72.7,72.8,72.8,72.4,71.6,70.6,69.7,68.9,68.2,67.5,66.7,65.9,&
        65.0,64.3,63.7,63.4,63.1,62.8,62.4,62.0,61.5,61.2,61.0,61.0,&
        59.6,60.0,60.4,60.7,60.7,60.7,60.5,60.4,60.6,61.1,61.8,62.5,&
        63.0,63.4,63.7,64.1,64.6,65.4,66.2,67.0,67.7,68.5,69.3,70.1,&
        70.6,70.9,70.9,70.4,69.3,68.0,66.9,66.0,65.2,64.5,63.7,63.0,&
        62.3,61.7,61.2,60.9,60.6,60.4,60.2,60.1,59.8,59.6,59.4,59.4,&
        58.5,58.8,59.2,59.4,59.4,59.2,58.9,58.9,59.2,59.7,60.5,61.2,&
        61.7,62.0,62.3,62.6,63.3,64.1,65.1,65.9,66.4,66.9,67.5,68.4,&
        69.0,69.2,68.9,68.1,66.7,65.3,64.2,63.4,62.7,62.0,61.1,60.4,&
        59.7,59.2,58.9,58.6,58.4,58.4,58.4,58.5,58.5,58.4,58.3,58.3,&
        57.6,57.8,57.9,57.9,57.7,57.5,57.4,57.5,57.9,58.6,59.3,59.8,&
        60.3,60.5,60.8,61.2,62.0,63.0,64.1,64.9,65.4,65.7,66.2,66.8,&
        67.2,66.8,66.0,64.8,63.4,62.1,61.3,60.7,60.1,59.3,58.2,57.1,&
        56.6,56.3,56.1,56.0,56.0,56.2,56.6,56.9,57.1,57.3,57.4,57.5,&
        54.3,54.9,55.4,55.6,55.6,55.3,55.1,55.0,55.3,55.8,56.5,57.2,&
        57.8,58.3,58.7,59.4,60.3,61.2,62.0,62.8,63.4,63.9,64.1,64.1,&
        63.8,63.3,62.5,61.5,60.2,58.9,57.7,56.7,56.0,55.5,55.1,54.8,&
        54.6,54.3,53.9,53.4,53.2,53.1,53.3,53.4,53.5,53.5,53.6,53.8,&
        52.9,53.6,54.2,54.6,54.6,54.3,54.0,53.7,53.8,54.1,54.8,55.7,&
        56.4,56.9,57.4,58.0,58.8,59.6,60.2,60.8,61.7,62.4,62.6,62.2,&
        61.5,60.7,59.8,58.9,57.9,56.8,55.5,54.5,53.8,53.4,53.3,53.3,&
        53.5,53.2,52.7,52.1,51.7,51.7,51.7,51.9,51.9,52.0,52.1,52.4,&
        51.8,52.5,53.2,53.6,53.6,53.3,52.9,52.5,52.4,52.6,53.3,54.2,&
        55.0,55.6,56.1,56.7,57.3,57.9,58.1,58.6,59.7,60.7,60.9,60.3,&
        59.2,58.1,57.1,56.3,55.5,54.6,53.5,52.5,51.8,51.4,51.4,51.7,&
        52.0,51.9,51.4,50.8,50.4,50.3,50.4,50.5,50.6,50.7,50.8,51.2,&
        50.9,51.7,52.4,52.9,52.9,52.5,52.0,51.5,51.3,51.4,52.1,53.0,&
        53.9,54.5,55.0,55.5,56.0,56.4,56.4,56.8,58.0,59.3,59.5,58.7,&
        57.4,56.1,54.9,54.1,53.5,52.8,51.8,50.8,50.1,49.7,49.8,50.4,&
        50.9,50.9,50.3,49.7,49.3,49.2,49.3,49.4,49.5,49.6,49.8,50.2/
    
    if(xkp > 9.0) xkp=9.0
    kp1=int(xkp)+1
    xkp1=int(xkp)*1.0
    kp2=kp1+1
    if(kp2 > 10) kp2=10
    do i=1,48 
        ab_mlat(i)=zp_mlat(i,kp1)+(xkp-xkp1)*&
            (zp_mlat(i,kp2)-zp_mlat(i,kp1))
    enddo
    
    cgmlat=-99.99
    if(xmlt < 0.0) return
    
    do i=1,48 
        ab_mlt(i)=(i-1)*.5
    enddo
    i1=int(xmlt/0.5)+1
    if(i1 >= 48) i1=1
    i2=i1+1
    
    s1=(zp_mlat(i2,kp1)-zp_mlat(i1,kp1))/(ab_mlt(i2)-ab_mlt(i1))
    zmlkp1=zp_mlat(i1,kp1)+(xmlt-ab_mlt(i1))*s1
    s2=(zp_mlat(i2,kp2)-zp_mlat(i1,kp2))/(ab_mlt(i2)-ab_mlt(i1))
    zmlkp2=zp_mlat(i1,kp2)+(xmlt-ab_mlt(i1))*s2
    
    cgmlat=zmlkp1+(xkp-xkp1)*(zmlkp2-zmlkp1)
    return
end
!
!
