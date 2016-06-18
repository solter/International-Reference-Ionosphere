! iritec.for, version number can be found at the end of this comment.
!-----------------------------------------------------------------------        
!
! contains IRIT13, IONCORR, IRI_TEC subroutines to computed the 
! total ionospheric electron content (TEC)
!
!-----------------------------------------------------------------------        
!-----------------------------------------------------------------------        
! Program for numerical integration of IRI-94 profiles from h=100km
! to h=alth. 
!       
!  INPUT:  ALATI,ALONG  LATITUDE NORTH AND LONGITUDE EAST IN DEGREES
!          jmag         =0 geographic   =1 geomagnetic coordinates
!          jf(1:50)     =.true./.false. flags; explained in IRISUB.FOR
!          iy,md        date as yyyy and mmdd (or -ddd)
!          hour         decimal hours LT (or UT+25)
!          hbeg,hend    upper and lower integration limits in km
! 
!  OUTPUT: TEC          Total Electron Content in m-2
!          tecb,tect    percentage of bottomside and topside content
!-----------------------------------------------------------------------        
subroutine IRIT13(ALATI,ALONG,jmag,jf,iy,md,hour,hbeg,hend,tec,tecb,tect)

    dimension       outf(20,1000),oarr(100)
    logical         jf(50)

    !  select various options and choices for IRI-94
    tec = -111.
    tect= -111.
    tecb= -111.

    ! initialize IRI parameter in COMMON blocks
    abeg=hbeg
    aend=hend
    astp=hend-hbeg
    call IRI_SUB(JF,JMAG,ALATI,ALONG,IY,MD,HOUR,abeg,aend,astp,OUTF,OARR)

    !  calculate total electron content (TEC) in m-2 using the
    !  stepsize selection 2 (highest accuracy)
    call iri_tec (hbeg,hend,2,tec,tect,tecb)

    return
end subroutine IRIT13


!-----------------------------------------------------------------------        
! computes ionospheric correction IONCORR (in m) for given vertical
! ionospheric electron content TEC (in m-2) and frequency f (in Hz)
!-----------------------------------------------------------------------        
real function ioncorr(tec,f)
    ioncorr = 40.3 * tec / f**2
    return
end


!-----------------------------------------------------------------------        
! subroutine to compute the total ionospheric content
! INPUT:      
!   hstart  altitude (in km) where integration should start
!   hend    altitude (in km) where integration should end
!   istep   =0 [fast, but higher uncertainty <5%]
!           =1 [standard, recommended]
!           =2 [stepsize of 1 km; best TEC, longest CPU time]
! OUTPUT:
!   tectot  total ionospheric content in tec-units (10^16 m^-2)
!   tectop  topside content (in %)
!   tecbot  bottomside content (in %)
!
! The different stepsizes for the numerical integration are 
! defined as follows (h1=100km, h2=hmF2-10km, h3=hmF2+10km, 
! h4=hmF2+150km, h5=hmF2+250km):
!       istep   h1-h2   h2-h3   h3-h4   h4-h5   h5-hend
!       0       2.0km   1.0km   2.5km   exponential approximation
!       1       2.0km   1.0km   2.5km   10.0km  30.0km
!       2       1.0km   0.5km   1.0km   1.0km   1.0km   
!
!-----------------------------------------------------------------------        
subroutine iri_tec (hstart,hend,istep,tectot,tectop,tecbot)

    logical         expo
    dimension       step(5),hr(6)
    logical         f1reg
    common  /block1/hmf2,xnmf2,hmf1,f1reg &
            /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,HHALF,TAU
    ! NEW-GUL------------------------------
    !     &         /QTOP/Y05,H05TOP,QF,XNETOP,XM3000,hht,TAU

    !test   
    save

    expo    = .false.
    numstep = 5
    xnorm   = xnmf2/1000.

    ! NEW-2003: Half-density: XNETOP at htop in [hmf2,1000 km)
    xxx  = xnmf2/2.
    ht1  = hmf2
    xne1 = xnmf2
    ht2  = ht1
    xne2 = xne1
    hht  = 0.0

    ! NEW-2003: Half-density: XNETOP at htop in [hmf2,1000 km)
    hr(1) = 100.
    hr(2) = hmf2-10.
    hr(3) = hmf2+10.
    hr(4) = hmf2+150.
    hr(5) = hmf2+250.
    hr(6) = hend

    do i=2,6 
        if (hr(i) > hend) hr(i) = hend
    end do

    if (istep == 0) then 
        step(1) = 2.0
        step(2) = 1.0
        step(3) = 2.5
        step(4) = 5.
        if (hend > hr(5)) expo = .true.
    end if

    if (istep == 1) then
        step(1) = 2.0
        step(2) = 1.0
        step(3) = 2.5
        step(4) = 10.0
        step(5) = 30.0
    end if

    if (istep == 2) then
        step(1) = 1.0
        step(2) = 0.5
        step(3) = 1.0
        step(4) = 1.0
        step(5) = 1.0
    end if

    sumtop = 0.0
    sumbot = 0.0

    ! find the starting point for the integration
    i  = 0
    ia = 1

    do
        i  = i + 1
        h  = hr(i)
        if (hstart > h) then
            hr(i) = hstart
            ia    = i
        else
            exit
        end if
    end do

    ! start the numerical integration
    i    = ia
    h    = hr(i)
    hu   = hr(i+1)
    delx = step(i)
    do
        h    = h + delx
        hh   = h

        if (h >= hu) then
            delx = hu - h + delx
            hx   = hu - delx/2.
            YNE  = XE_1(hx)
            if( (hx > hmf2).and.(yne > xnmf2) ) yne=xnmf2
            yyy  = yne * delx / xnorm
            i    = i + 1
            if(i < 6) then
                h    = hr(i)
                hu   = hr(i+1)
                delx = step(i)
            end if
        else
            hx = h - delx/2.
            YNE = XE_1(hx)
            if( (hx > hmf2).and.(yne > xnmf2) ) yne = xnmf2
            yyy = yne * delx / xnorm
        end if

        if (hx <= hmf2) then
            sumbot = sumbot + yyy
        else
            sumtop = sumtop + yyy
            ! NEW-GUL: remember xne2 at ht2 :
            ht2  = hx
            xne2 = yne
            ! NEW-GUL------------------------------
        end if

        ! NEW-GUL: interpolate for htop
        if ( (hx > hmf2).and.(hht <= 0.0) ) then
            if ((xxx <= xne1).and.(xxx > xne2)) then
                hht  = ht1 + (ht2 - ht1)/(xne2 - xne1)*(xxx - xne1)
            else
                ht1  = ht2
                xne1 = xne2
            endif
        end if
        ! NEW-GUL------------------------------

        if ( expo .and. hh >= hr(4) ) cycle

        if (hh >= hend .or. i >= 6) then
            zzz    = sumtop + sumbot
            tectop = sumtop / zzz * 100.
            tecbot = sumbot / zzz * 100.
            tectot = zzz * xnmf2    
            return
        end if
    end do

    num_step = 3
    hei_top  = hr(4)
    hei_end  = hend
    top_end  = hei_end - hei_top
    del_hei  = top_end / num_step
    xntop    = xe_1(hei_end)/xnmf2

    if(xntop > 0.9999) then
        ss_t = top_end  
    else
        hei_2 = hei_top
        hei_3 = hei_2 + del_hei
        hei_4 = hei_3 + del_hei
        hei_5 = hei_end

        hss = top_end / 4.
        xkk = exp ( - top_end / hss ) - 1.
        x_2 = hei_2
        x_3 = hei_top - hss * alog( xkk*(hei_3 - hei_top)/top_end + 1. ) 
        x_4 = hei_top - hss * alog( xkk*(hei_4 - hei_top)/top_end + 1. )
        x_5 = hei_end

        ed_2 = min(xe_1(x_2)/xnmf2, 1.)
        ed_3 = min(xe_1(x_3)/xnmf2, 1.)
        ed_4 = min(xe_1(x_4)/xnmf2, 1.)
        ed_5 = xntop

        if(ed_3 == ed_2) then
            ss_2 = ed_3 * (x_3 - x_2)
        else
            ss_2 = ( ed_3 - ed_2 ) * ( x_3 - x_2 ) / alog( ed_3 / ed_2 )
        end if

        if(ed_4 == ed_3) then
            ss_3 = ed_4 * (x_4 - x_3)
        else
            ss_3 = ( ed_4 - ed_3 ) * ( x_4 - x_3 ) / alog( ed_4 / ed_3 )
        end if

        if(ed_5 == ed_4) then
            ss_4 = ed_5 * (x_5 - x_4)
        else
            ss_4 = ( ed_5 - ed_4 ) * ( x_5 - x_4 ) / alog( ed_5 / ed_4 )
        end if

        ss_t = ss_2 + ss_3 + ss_4 

    end if

    sumtop = sumtop + ss_t * 1000.
    zzz    = sumtop + sumbot
    tectop = sumtop / zzz * 100.
    tecbot = sumbot / zzz * 100.
    tectot = zzz * xnmf2

    return
end subroutine iri_tec
