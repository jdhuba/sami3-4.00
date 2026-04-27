!******************************************
!******************************************

!   heating.F90
!
!   Ion and electron temperature equations.
!
!   References:
!     Schunk & Nagy, Ionospheres, Cambridge Press
!     Millward et al., in STEP Handbook (red book), p. 269, 1996
!     Rees, Physics and Chemistry of the Upper Atmosphere, Cambridge, 1989

!******************************************
!******************************************

! -----------------------------------------------------------------------
! convfac  = amu / (bolt * 3)  -- frictional heating conversion factor.
! Derived purely from parameter_mod constants; defined once here and
! used by itemp (htemp / hetemp / otemp wrappers) and etemp.
! twothirds = 2/3 exactly -- replaces the inconsistent mix of 0.6667,
! .66667, and 0.66667 scattered throughout the original.
! -----------------------------------------------------------------------

!******************************************
!******************************************

!   itemp  --  generic ion temperature subroutine
!
!   Replaces the three identical subroutines htemp / hetemp / otemp.
!   The only difference between those routines was the ion species index
!   passed to tisolv and used for the Coulomb logarithm, thermal
!   conductivity, and exchange terms.  That index is now an argument.
!
!   Callers:
!     htemp  -> call itemp(..., ion_neut=pth,  ion_idx=pthp,  ...)
!     hetemp -> call itemp(..., ion_neut=pthe, ion_idx=pthep, ...)
!     otemp  -> call itemp(..., ion_neut=pto,  ion_idx=ptop,  ...)
!
!   ion_neut : neutral species index matching this ion (for Coulomb log)
!   ion_idx  : ion species index (for all exchange and solver terms)

!******************************************
!******************************************

    subroutine itemp(tti, tiold, tvn, nuin, nfl, nll, ion_neut, ion_idx)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use atomic_mod
    use exb_mod
    use grid_mod

    implicit none

    ! Arguments
    real,    intent(inout) :: tti(nz)
    real,    intent(in)    :: tiold(nz)
    real,    intent(in)    :: tvn(nz,nl)
    real,    intent(in)    :: nuin(nz,nion,nneut)
    integer, intent(in)    :: nfl, nll
    integer, intent(in)    :: ion_neut   ! neutral index for Coulomb log (pth/pthe/pto)
    integer, intent(in)    :: ion_idx    ! ion index for exchange terms (pthp/pthep/ptop)

    ! Local variables
    integer :: i, ni, nn
    real    :: kapi(nz), s1i(nz), s2i(nz), s3i(nz), s4i(nz), s5i(nz)
    real    :: s6i(nz), s7i(nz), divvexb(nz)
    real    :: coulomb_log    ! Coulomb logarithm (called 'lambda' in original)
    real    :: schunkfac, redmass, tfac, xs6i, xs7i
    real    :: nzh_r, vexbeq
    integer :: nzh

    ! Derived constants (same value as original convfac = amu/bolt/3.)
    real, parameter :: convfac   = amu / bolt / 3.0
    real, parameter :: twothirds = 2.0 / 3.0

    ! ------------------------------------------------------------------
    ! Initialise source term arrays
    ! ------------------------------------------------------------------
    s1i(:)  = 0.0
    s2i(:)  = 0.0
    s3i(:)  = 0.0
    s4i(:)  = 0.0
    s5i(:)  = 0.0
    s6i(:)  = 0.0
    s7i(:)  = 0.0
    kapi(:) = 0.0

    ! ------------------------------------------------------------------
    ! Compute thermal conductivity and heating/cooling source terms
    ! ------------------------------------------------------------------
    do i = 1, nz

        ! Coulomb logarithm (Schunk & Nagy)
        coulomb_log = 23.0 - &
            0.5 * log( deni(i,nfl,nll,ion_neut) / &
                      (ti(i,nfl,nll,ion_neut) / evtok)**3 )

        ! Thermal conductivity (Schunk factor accounts for other ions)
        schunkfac = 0.0
        do ni = nion1, nion2
            if (ni /= ion_neut) then
                schunkfac = schunkfac + &
                    deni(i,nfl,nll,ni) / deni(i,nfl,nll,ion_neut) * &
                    sqrt( ami(ni) / (ami(ni) + ami(ion_neut))**5 ) * &
                    ( 3.0*ami(ion_neut)**2 + &
                      1.6*ami(ion_neut)*ami(ni) + &
                      1.3*ami(ni)**2 )
            endif
        enddo

        kapi(i) = 15.0 * 3.1e4 * sqrt( ti(i,nfl,nll,ion_neut)**5 ) / &
                  sqrt( ami(ion_neut) ) / &
                  ( 1.0 + 1.75*schunkfac ) / coulomb_log
        kapi(i) = twothirds * kapi(i) * evtok

        ! Neutral heating/cooling (s2i) and frictional heating (s3i)
        do nn = 1, nneut
            redmass = ami(ion_idx) * amn(nn) / ( ami(ion_idx) + amn(nn) )**2
            s2i(i)  = s2i(i) + 2.0 * nuin(i,ion_idx,nn) * redmass
            s3i(i)  = s3i(i) + convfac * amn(nn) &
                      * abs( vsi(i,nfl,nll,ion_idx) - tvn(i,nll) )**2 &
                      * 2.0 * nuin(i,ion_idx,nn) * redmass
        enddo

        s1i(i) = s2i(i) * tn(i,nfl,nll)

        ! Electron exchange (s4i, s5i)
        s4i(i) = 7.7e-6 * ne(i,nfl,nll) / ami(ion_idx) &
                 / te(i,nfl,nll) / sqrt( te(i,nfl,nll) ) &
                 * twothirds * evtok
        s5i(i) = s4i(i) * te(i,nfl,nll)

        ! Other-ion exchange (s6i, s7i)
        do ni = nion1, nion2
            if ( ni /= ion_idx ) then
                tfac  = ti(i,nfl,nll,ion_idx) / ami(ion_idx) &
                      + ti(i,nfl,nll,ni)      / ami(ni)
                xs6i  = 3.3e-4 * deni(i,nfl,nll,ni) / ami(ion_idx) / ami(ni) &
                        / tfac / sqrt(tfac) * twothirds * evtok
                xs7i  = xs6i * ti(i,nfl,nll,ni)
                s6i(i) = s6i(i) + xs6i
                s7i(i) = s7i(i) + xs7i
            endif
        enddo

    enddo

    ! ------------------------------------------------------------------
    ! MS: Neglected term -- divergence of ExB drift
    ! Divergence of the ExB drift; requires equatorial drift
    ! ------------------------------------------------------------------
    nzh    = nz / 2
    vexbeq = vexbp(nzh, nfl, nll)
    do i = 1, nz
        divvexb(i) = 6.0 * vexbeq / &
                     ( ps(i,nfl,nll) * re * 1.0e5 ) * &
                     cos( blats(i,nfl,nll) * po180 )**2 * &
                     ( 1.0 + sin( blats(i,nfl,nll) * po180 )**2 ) / &
                     ( 1.0 + 3.0 * sin( blats(i,nfl,nll) * po180 )**2 )**2
        s2i(i) = s2i(i) - (1.0/3.0) * divvexb(i)
    enddo

    call tisolv(tti, tiold, kapi, s1i, s2i, s3i, s4i, s5i, s6i, s7i, &
                ion_idx, nfl, nll)

    end subroutine itemp


!******************************************
!******************************************

!   htemp  --  H+ ion temperature
!   Thin wrapper around itemp; ion_neut=pth (H neutral), ion_idx=pthp (H+)

!******************************************
!******************************************

    subroutine htemp(tti, tiold, tvn, nuin, nfl, nll)

    use parameter_mod

    implicit none

    real,    intent(inout) :: tti(nz)
    real,    intent(in)    :: tiold(nz)
    real,    intent(in)    :: tvn(nz,nl)
    real,    intent(in)    :: nuin(nz,nion,nneut)
    integer, intent(in)    :: nfl, nll

    call itemp(tti, tiold, tvn, nuin, nfl, nll, pth, pthp)

    end subroutine htemp


!******************************************
!******************************************

!   hetemp  --  He+ ion temperature
!   Thin wrapper around itemp; ion_neut=pthe (He neutral), ion_idx=pthep (He+)

!******************************************
!******************************************

    subroutine hetemp(tti, tiold, tvn, nuin, nfl, nll)

    use parameter_mod

    implicit none

    real,    intent(inout) :: tti(nz)
    real,    intent(in)    :: tiold(nz)
    real,    intent(in)    :: tvn(nz,nl)
    real,    intent(in)    :: nuin(nz,nion,nneut)
    integer, intent(in)    :: nfl, nll

    call itemp(tti, tiold, tvn, nuin, nfl, nll, pthe, pthep)

    end subroutine hetemp


!******************************************
!******************************************

!   otemp  --  O+ ion temperature
!   Thin wrapper around itemp; ion_neut=pto (O neutral), ion_idx=ptop (O+)

!******************************************
!******************************************

    subroutine otemp(tti, tiold, tvn, nuin, nfl, nll)

    use parameter_mod

    implicit none

    real,    intent(inout) :: tti(nz)
    real,    intent(in)    :: tiold(nz)
    real,    intent(in)    :: tvn(nz,nl)
    real,    intent(in)    :: nuin(nz,nion,nneut)
    integer, intent(in)    :: nfl, nll

    call itemp(tti, tiold, tvn, nuin, nfl, nll, pto, ptop)

    end subroutine otemp


!******************************************
!******************************************

!   etemp  --  electron temperature

!******************************************
!******************************************

    subroutine etemp(tte, te_old, phprodr, nfl, nll)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use atomic_mod
    use exb_mod
    use grid_mod

    implicit none

    ! Arguments
    real,    intent(inout) :: tte(nz)
    real,    intent(in)    :: te_old(nz)
    real,    intent(in)    :: phprodr(nz,nion)
    integer, intent(in)    :: nfl, nll

    ! Local variables
    integer :: i, ni, nn, nzh
    integer :: izs, izn
    integer :: iz300s(nf,nl), iz300n(nf,nl)

    real    :: kape(nz), s1e(nz), s2e(nz), s3e(nz), s4e(nz), s5e(nz)
    real    :: qphe(nz), phprod(nz), qen(nz,nneut)
    real    :: ratio(nz), divvexb(nz)

    real    :: fac1, fac2, fac3, akpefac
    real    :: xs3e, xs4e
    real    :: xarg, x, earg, epsi
    real    :: facts, factn
    real    :: ne300s, ne300n
    real    :: o2300, n2300, o300, phprod300
    real    :: xarg300, x300, earg300, epsi300
    real    :: q0s, q0n
    real    :: xbms, xbmn
    real    :: dels300s, dels300n
    real    :: xn, xss, xints, xintn, xqs, xqn
    real    :: vexbeq

    real, parameter :: twothirds = 2.0 / 3.0

    ! ------------------------------------------------------------------
    ! Initialise source term arrays
    ! ------------------------------------------------------------------
    s1e(:)  = 0.0
    s2e(:)  = 0.0
    s3e(:)  = 0.0
    s4e(:)  = 0.0
    s5e(:)  = 0.0
    kape(:) = 0.0
    qen(:,:) = 0.0

    ! ------------------------------------------------------------------
    ! Thermal conductivity and neutral/ion exchange source terms
    ! ------------------------------------------------------------------
    do i = 1, nz

        ! Electron thermal conductivity (modified by neutral collisions)
        fac1 = denn(i,nfl,nll,pto)  * 1.1e-16 &
               * ( 1.0 + 5.7e-4 * te(i,nfl,nll) )
        fac2 = denn(i,nfl,nll,ptn2) * 2.82e-17 &
               * ( 1.0 - 1.2e-4 * te(i,nfl,nll) ) * sqrt( te(i,nfl,nll) )
        fac3 = denn(i,nfl,nll,pto2) * 2.2e-16 &
               * ( 1.0 + 3.6e-2 * sqrt( te(i,nfl,nll) ) )
        akpefac = fac1 + fac2 + fac3

        kape(i) = 7.7e5 * sqrt( te(i,nfl,nll)**5 ) * twothirds * evtok &
                  / ( 1.0 + 3.22e4 * ( te(i,nfl,nll)**2 / &
                  ne(i,nfl,nll) * akpefac ) )

        ! Neutral (Tn-Te) cooling terms

        ! N2 -- vibrational state from red book (p. 269) Millward et al.
        qen(i,ptn2) = twothirds * evtok * denn(i,nfl,nll,ptn2) * &
            ( 1.2e-19 * ( 1.0 - 1.2e-4 * te(i,nfl,nll) ) &
              * te(i,nfl,nll) &
              + 2.0e-14 / sqrt( te(i,nfl,nll) ) &
              + 6.5e-22 * ( tn(i,nfl,nll) - 310.0 )**2 * &
                exp( 0.0023 * ( te(i,nfl,nll) - tn(i,nfl,nll) ) ) )

        ! O2
        qen(i,pto2) = twothirds * evtok * denn(i,nfl,nll,pto2) * &
            ( 7.9e-19 * ( 1.0 + 3.6e-2 * sqrt( te(i,nfl,nll) ) ) &
              * sqrt( te(i,nfl,nll) ) &
              + 7.0e-14 / sqrt( te(i,nfl,nll) ) )

        ! O
        qen(i,pto) = twothirds * 7.2e-18 * evtok * &
                     denn(i,nfl,nll,pto) * sqrt( te(i,nfl,nll) )

        ! H
        qen(i,pth) = twothirds * 6.3e-16 * evtok * denn(i,nfl,nll,pth) * &
                     ( 1.0 - 1.35e-4 * te(i,nfl,nll) ) * sqrt( te(i,nfl,nll) )

        do nn = 1, nneut
            s2e(i) = s2e(i) + qen(i,nn)
        enddo

        s1e(i) = s2e(i) * tn(i,nfl,nll)

        ! Ion (Ti-Te) exchange terms
        do ni = nion1, nion2
            xs3e   = 7.7e-6 * deni(i,nfl,nll,ni) / ami(ni) &
                     / te(i,nfl,nll) / sqrt( te(i,nfl,nll) ) &
                     * twothirds * evtok
            xs4e   = xs3e * ti(i,nfl,nll,ni)
            s3e(i) = s3e(i) + xs3e
            s4e(i) = s4e(i) + xs4e
        enddo

    enddo

    ! ------------------------------------------------------------------
    ! Photoelectron heating
    ! Red book (Millward et al. p. 269)
    ! ------------------------------------------------------------------

    ! Total ion photoproduction (= photoelectron production rate)
    phprod(:) = 0.0
    do i = 1, nz
        do ni = nion1, nion2
            phprod(i) = phprod(i) + phprodr(i,ni) * denn(i,nfl,nll,ni)
        enddo
    enddo

    ! Ionosphere/plasmasphere ratio used to locate ~300 km boundary
    do i = 1, nz
        ratio(i) = ne(i,nfl,nll) / &
                   ( 0.1 * denn(i,nfl,nll,pto) + &
                     denn(i,nfl,nll,pto2) + &
                     denn(i,nfl,nll,ptn2) )
    enddo

    ! Find south 300 km index (ascending from bottom)
    ! Loop guards prevent out-of-bounds: index tested before ratio access
    iz300s(nfl,nll) = 1
    i = 1
    do while ( i < nz )
        if ( ratio(i) > 3.0e-3 ) exit
        iz300s(nfl,nll) = i
        i = i + 1
    enddo

    ! Find north 300 km index (descending from top)
    iz300n(nfl,nll) = nz
    i = nz
    do while ( i > 1 )
        if ( ratio(i) > 3.0e-3 ) exit
        iz300n(nfl,nll) = i
        i = i - 1
    enddo

    if ( iz300s(nfl,nll) > iz300n(nfl,nll) ) then

        ! Entire field line is above 300 km: apply heating everywhere
        do i = 1, nz
            xarg    = ne(i,nfl,nll) / &
                      ( denn(i,nfl,nll,pto2) + &
                        denn(i,nfl,nll,ptn2) + &
                        0.1 * denn(i,nfl,nll,pto) )
            x       = log( xarg )
            earg    = 12.75 + 6.941*x + 1.166*x**2 + 0.08034*x**3 + 0.001996*x**4
            epsi    = exp( -earg )
            qphe(i) = epsi * phprod(i)
        enddo

    else

        ! Field line straddles 300 km boundary: three-region treatment

        ! South low-altitude region
        do i = 1, iz300s(nfl,nll)
            xarg    = ne(i,nfl,nll) / &
                      ( denn(i,nfl,nll,pto2) + &
                        denn(i,nfl,nll,ptn2) + &
                        0.1 * denn(i,nfl,nll,pto) )
            x       = log( xarg )
            earg    = 12.75 + 6.941*x + 1.166*x**2 + 0.08034*x**3 + 0.001996*x**4
            epsi    = exp( -earg )
            qphe(i) = epsi * phprod(i)
        enddo

        ! Smooth at south 300 km boundary
        izs   = iz300s(nfl,nll)
        facts = ( 3.0e-3 - ratio(izs) ) / ( ratio(izs+1) - ratio(izs) )

        ne300s    = ne(izs,nfl,nll)   + ( ne(izs+1,nfl,nll)   - ne(izs,nfl,nll)   ) * facts
        o2300     = denn(izs,nfl,nll,pto2) + ( denn(izs+1,nfl,nll,pto2) - denn(izs,nfl,nll,pto2) ) * facts
        n2300     = denn(izs,nfl,nll,ptn2) + ( denn(izs+1,nfl,nll,ptn2) - denn(izs,nfl,nll,ptn2) ) * facts
        o300      = denn(izs,nfl,nll,pto)  + ( denn(izs+1,nfl,nll,pto)  - denn(izs,nfl,nll,pto)  ) * facts
        phprod300 = phprod(izs) + ( phprod(izs+1) - phprod(izs) ) * facts

        xarg300   = ne300s / ( o2300 + n2300 + 0.1*o300 )
        x300      = log( xarg300 )
        earg300   = 12.75 + 6.941*x300 + 1.166*x300**2 + 0.08034*x300**3 + 0.001996*x300**4
        epsi300   = exp( -earg300 )
        q0s       = epsi300 * phprod300 / ne300s

        ! North low-altitude region
        do i = iz300n(nfl,nll), nz
            xarg    = ne(i,nfl,nll) / &
                      ( denn(i,nfl,nll,pto2) + &
                        denn(i,nfl,nll,ptn2) + &
                        0.1 * denn(i,nfl,nll,pto) )
            x       = log( xarg )
            earg    = 12.75 + 6.941*x + 1.166*x**2 + 0.08034*x**3 + 0.001996*x**4
            epsi    = exp( -earg )
            qphe(i) = epsi * phprod(i)
        enddo

        ! Smooth at north 300 km boundary
        izn   = iz300n(nfl,nll)
        factn = ( 3.0e-3 - ratio(izn) ) / ( ratio(izn-1) - ratio(izn) )

        ne300n    = ne(izn,nfl,nll)   + ( ne(izn-1,nfl,nll)   - ne(izn,nfl,nll)   ) * factn
        o2300     = denn(izn,nfl,nll,pto2) + ( denn(izn-1,nfl,nll,pto2) - denn(izn,nfl,nll,pto2) ) * factn
        n2300     = denn(izn,nfl,nll,ptn2) + ( denn(izn-1,nfl,nll,ptn2) - denn(izn,nfl,nll,ptn2) ) * factn
        o300      = denn(izn,nfl,nll,pto)  + ( denn(izn-1,nfl,nll,pto)  - denn(izn,nfl,nll,pto)  ) * factn
        phprod300 = phprod(izn) + ( phprod(izn-1) - phprod(izn) ) * factn

        xarg300   = ne300n / ( o2300 + n2300 + 0.1*o300 )
        x300      = log( xarg300 )
        earg300   = 12.75 + 6.941*x300 + 1.166*x300**2 + 0.08034*x300**3 + 0.001996*x300**4
        epsi300   = exp( -earg300 )
        q0n       = epsi300 * phprod300 / ne300n

        xbms    = bms(izs,nfl,nll) + ( bms(izs+1,nfl,nll) - bms(izs,nfl,nll) ) * facts
        xbmn    = bms(izn,nfl,nll) + ( bms(izn-1,nfl,nll) - bms(izn,nfl,nll) ) * factn

        dels300s = dels(iz300s(nfl,nll),   nfl,nll) * facts
        dels300n = dels(iz300n(nfl,nll)-1, nfl,nll) * factn

        ! MS: Cleaner xn calculation; set bottom integration bound to 300 km.
        xn = 0.5 * ( ne(iz300n(nfl,nll)-1,nfl,nll) + ne300n ) * &
             ( dels(iz300n(nfl,nll)-1,nfl,nll) - dels300n )
        do i = iz300n(nfl,nll)-2, iz300s(nfl,nll)+1, -1
            xn = xn + 0.5 * ( ne(i,nfl,nll) + ne(i+1,nfl,nll) ) * dels(i,nfl,nll)
        enddo

        if ( q0s < 0.0 .OR. q0n < 0.0 ) then
            print *, ' q0s = ', q0s, ' q0n = ', q0n, ' nfl = ', nfl
        endif

        ! Put in dels (arc length along field line); compute middle region heating
        xss = 0.0
        do i = iz300s(nfl,nll)+1, iz300n(nfl,nll)-1
            if ( i == iz300s(nfl,nll)+1 ) then
                xss = xss + 0.5 * ( ne300s + ne(i,nfl,nll) ) * &
                      ( dels(iz300s(nfl,nll),nfl,nll) - dels300s )
            else
                xss = xss + 0.5 * ( ne(i,nfl,nll) + ne(i-1,nfl,nll) ) * dels(i-1,nfl,nll)
                xn  = xn  - 0.5 * ( ne(i,nfl,nll) + ne(i-1,nfl,nll) ) * dels(i-1,nfl,nll)
            endif

            xints   = cqe * xss
            xintn   = cqe * xn
            xqs     = ne(i,nfl,nll) * q0s * bms(i,nfl,nll) / xbms * exp(-xints)
            xqn     = ne(i,nfl,nll) * q0n * bms(i,nfl,nll) / xbmn * exp(-xintn)
            qphe(i) = xqs + xqn
        enddo

    endif

    ! Photoelectron heating source term
    do i = 1, nz
        s5e(i) = twothirds * evtok * qphe(i) / ne(i,nfl,nll) * s5e_fac
    enddo

    ! ------------------------------------------------------------------
    ! MS: Neglected term -- divergence of ExB drift
    ! Divergence of the ExB drift; requires equatorial drift
    ! ------------------------------------------------------------------
    nzh    = nz / 2
    vexbeq = vexbp(nzh, nfl, nll)
    do i = 1, nz
        divvexb(i) = 6.0 * vexbeq / &
                     ( ps(i,nfl,nll) * re * 1.0e5 ) * &
                     cos( blats(i,nfl,nll) * po180 )**2 * &
                     ( 1.0 + sin( blats(i,nfl,nll) * po180 )**2 ) / &
                     ( 1.0 + 3.0 * sin( blats(i,nfl,nll) * po180 )**2 )**2
        s2e(i) = s2e(i) - (1.0/3.0) * divvexb(i)
    enddo

    call tesolv(tte, te_old, kape, s1e, s2e, s3e, s4e, s5e, nfl, nll)

    end subroutine etemp
