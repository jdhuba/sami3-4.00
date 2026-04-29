!******************************************
!******************************************

!            photprod

!******************************************
!******************************************

! photoproduction rates

    subroutine photprod ( phprodr, nfl, nll )

    use parameter_mod
    use variable_mod
    use namelist_mod
    use photo_production_mod
    use atomic_mod
    use grid_mod
    use chapmanv

    implicit none

    ! Arguments
    real,    intent(out) :: phprodr(nz,nion)
    integer, intent(in)  :: nfl, nll

    ! Local variables
    integer :: iz, i, j, l
    integer :: itheta0, itheta
    real    :: xmass(3)
    integer :: idx(3)
    real(8) :: ch1, xscale, xr   ! fix: real*8 -> real(8); atm_chapman removed (unused)
    real    :: hscale, coschi, rp, rp2
    real    :: exa, flx, pei_rate, ang, del, fluxntt

    ! fix: hcof is a constant derived entirely from parameter_mod constants.
    ! Computed once here as a local parameter rather than re-evaluated per call.
    real, parameter :: hcof = 1.e-5 * bolt / ( gzero * amu * re ** 2 )

    do iz = 1, nz
        coschi = cx(iz,nfl,nll)
        do j = nion1, nion2
            phprodr(iz,j) = 0.0
        enddo

        ! Only consider O, N2, O2 for absorption
        idx(1) = pto
        idx(2) = ptn2
        idx(3) = pto2

        rp  = alts(iz,nfl,nll) + re
        rp2 = rp * rp

        if ( coschi >= coschicrit(iz,nfl,nll) ) then  ! sun is up

            ! Daytime deposition

            do i = 1, 3
                hscale   = hcof * tn(iz,nfl,nll) * rp2 / amn(idx(i))
                xscale   = rp / hscale
                xr       = acos(coschi)
                ch1      = chapmanApprox(xscale,xr)
                xmass(i) = denn(iz,nfl,nll,idx(i)) * hscale * ch1 * 1.e5
            enddo

                ! EUVAC
                do l = 1, linesuv_euvac
                    exa = xmass(1) * sigabsdt_euvac(l,1) &
                        + xmass(2) * sigabsdt_euvac(l,2) &
                        + xmass(3) * sigabsdt_euvac(l,3)
                    if (exa > 22.) exa = 22.
                    flx = flux_euvac(l) * exp(-exa)
                    do j = nion1, nion2
                        phprodr(iz,j) = phprodr(iz,j) + sigidt_euvac(l,j) * flx
                    enddo
                enddo

            ! Photoelectron ionization (Buonsanto et al. 1992 via Siskind)
            ! fix: three sequential ifs replaced with if/else if/else to make
            !      branches mutually exclusive (boundary at 200 km was ambiguous)

            if ( alts(iz,nfl,nll) > 200.0 ) then
                pei_rate = 0.2
            else if ( alts(iz,nfl,nll) >= 140.0 ) then
                pei_rate = exp( -(alts(iz,nfl,nll) - 140.0) / 37.0 )
            else
                pei_rate = exp( -(alts(iz,nfl,nll) - 140.0) / 14.0 )
            end if

            phprodr(iz,ptop)  = phprodr(iz,ptop)  * (1.0 + pei_rate)
            phprodr(iz,ptn2p) = phprodr(iz,ptn2p) * (1.0 + pei_rate)
            phprodr(iz,pto2p) = phprodr(iz,pto2p) * 1.2

            ! Add nighttime ionization (daytime branch also includes scattered flux)

            ang     = acos(coschi)
            itheta0 = int(ang / po180) - 90
            ! fix: amax1(float(...),1.) -> max(real(...),1.)
            itheta  = int( max( real(itheta0), 1.0 ) )
            del     = ang/po180 - int(ang/po180)

            do l = 1, linesnt
                do j = nion1, nion2
                    if (itheta0 < 1) then
                        fluxntt = fluxnt(iz,nfl,nll,itheta,l)
                    else
                        fluxntt = fluxnt(iz,nfl,nll,itheta,  l) * (1.0 - del) &
                                + fluxnt(iz,nfl,nll,itheta+1,l) * del
                    endif
                    phprodr(iz,j) = phprodr(iz,j) + sigint(l,j) * fluxntt
                enddo
            enddo

        else    ! sun is down

            ! Nighttime deposition

            ang     = acos(coschi)
            itheta0 = int(ang / po180) - 90
            ! fix: amax1(float(...),1.) -> max(real(...),1.)
            itheta  = int( max( real(itheta0), 1.0 ) )
            del     = ang/po180 - int(ang/po180)

            do l = 1, linesnt
                do j = nion1, nion2
                    if (itheta0 < 1) then
                        fluxntt = fluxnt(iz,nfl,nll,itheta,l)
                    else
                        fluxntt = fluxnt(iz,nfl,nll,itheta,  l) * (1.0 - del) &
                                + fluxnt(iz,nfl,nll,itheta+1,l) * del
                    endif
                    phprodr(iz,j) = phprodr(iz,j) + sigint(l,j) * fluxntt
                enddo
            enddo

        endif
    enddo

    end subroutine photprod




!******************************************
!******************************************

!   sf_nighttime_flux  (new helper)
!
!   Shared theta-interpolation kernel for sf1026, sf584, sf304, sf1216.
!   Fixes issue #6: ~300 lines of duplicated code across four subroutines
!   consolidated here.
!
!   The caller fills f(1:nz, nfl, nll, known_theta_indices) with the
!   altitude-dependent flux values for the 4 known zenith angles before
!   calling this routine.  sf_nighttime_flux then:
!     1. Floors all known-theta values to >= 1.0
!     2. Interpolates (log10) to all 91 theta values (90 deg to 180 deg)
!
!   The interpolation for the last theta bin (ji==4) differs between
!   sf1026/sf1216 (extrapolate to jip1=5) and sf584/sf304 (interpolate
!   to a fixed endpoint of 1.0 at theta=180).  The logical argument
!   'extrap_last' selects the behaviour:
!     extrap_last = .true.  -> sf1026 / sf1216 style (jip1 = ji+1 always)
!     extrap_last = .false. -> sf584  / sf304  style (clamp at 180 deg)

!******************************************
!******************************************

    subroutine sf_nighttime_flux ( f, line, nfl, nll, extrap_last )

    use parameter_mod
    use variable_mod
    use photo_production_mod
    use grid_mod

    implicit none

    ! Arguments
    real,    intent(inout) :: f(nz,nf,nl,91)
    integer, intent(in)    :: line, nfl, nll
    logical, intent(in)    :: extrap_last

    ! Local variables
    integer :: i, k, j, k90, ji, ki, jip1, kip1

    ! fix: float(...) -> real(...)
    real    :: delk, flog

    ! 1. Floor all 4 known-theta values to >= 1.0
    !    fix: amax1(1., ...) -> max(1.0, ...)
    do k = 1, 4
        do i = 1, nz
            f(i,nfl,nll, int(thetant(line,k))+1-90) = &
                max(1.0, f(i,nfl,nll, int(thetant(line,k))+1-90))
        enddo
    enddo

    ! 2. Interpolate (log10) to all 91 theta values (90 to 180 deg)
    do k = 1, 91
        k90 = 90 + k - 1
        ji  = 1
        ki  = int(thetant(line,1))
        do j = 1, 4
            if ( k90 > int(thetant(line,j)) ) then
                ji = j
                ki = int(thetant(line,ji))
            endif
        enddo

        if ( extrap_last .or. ji /= 4 ) then
            jip1 = ji + 1
            kip1 = int(thetant(line,jip1))
            delk = real( int(thetant(line,jip1)) - int(thetant(line,ji)) )
            do i = 1, nz
                flog = log10(f(i,nfl,nll,ki+1-90)) &
                     + (k90 - ki) / delk &
                     * ( log10(f(i,nfl,nll,kip1+1-90)) &
                       - log10(f(i,nfl,nll,ki  +1-90)) )
                ! fix: 10 ** flog (integer base) -> 10.0 ** flog
                f(i,nfl,nll,k) = 10.0 ** flog
            enddo
        else
            ! ji == 4 and not extrap_last: interpolate toward f=1 at theta=180
            delk = real( 180 - int(thetant(line,ji)) )
            do i = 1, nz
                flog = log10(f(i,nfl,nll,ki+1-90)) &
                     + (k90 - ki) / delk &
                     * ( log10(1.0) - log10(f(i,nfl,nll,ki+1-90)) )
                ! fix: 10 ** flog (integer base) -> 10.0 ** flog
                f(i,nfl,nll,k) = 10.0 ** flog
            enddo
        endif
    enddo

    end subroutine sf_nighttime_flux


!******************************************
!******************************************

!            f1026

!******************************************
!******************************************

! subroutine to calculate the nighttime flux of
! lyman beta (1026) (note: line = 1)

    subroutine sf1026 ( f, line, nfl, nll )

    use parameter_mod
    use variable_mod
    use photo_production_mod
    use grid_mod

    implicit none

    ! Arguments
    real,    intent(inout) :: f(nz,nf,nl,91)
    integer, intent(in)    :: line, nfl, nll

    ! Local variables
    integer :: i, k, imax

    imax = 1

    ! Determine f for the 4 known values of theta
    do i = 1, nz
        if ( alts(i,nfl,nll) < zaltnt(line,1) ) then
            do k = 1, 4
                f(i,nfl,nll, int(thetant(line,k))+1-90) = 1.0
            enddo
        elseif ( zaltnt(line,1) <= alts(i,nfl,nll) .and. &
                 alts(i,nfl,nll) <= zaltnt(line,2) ) then
            f(i,nfl,nll, int(thetant(line,1))+1-90) = &
                1.4e8 * tanh( (alts(i,nfl,nll) - 90.) / 50. )
            f(i,nfl,nll, int(thetant(line,2))+1-90) = &
                3.8e7 * tanh( (alts(i,nfl,nll) - 90.) / 50. )
            f(i,nfl,nll, int(thetant(line,3))+1-90) = &
                1.4e7 * tanh( (alts(i,nfl,nll) - 93.) / 55. )
            f(i,nfl,nll, int(thetant(line,4))+1-90) = &
                9.2e6 * tanh( (alts(i,nfl,nll) - 94.) / 55. )
            imax = i
        else
            do k = 1, 4
                f(i,   nfl,nll, int(thetant(line,k))+1-90) = &
                f(imax,nfl,nll, int(thetant(line,k))+1-90)
            enddo
        endif
    enddo

    ! Interpolate to all values of theta (90 - 180)
    ! sf1026 uses extrap_last=.true. (jip1 = ji+1 always, no endpoint clamping)
    call sf_nighttime_flux(f, line, nfl, nll, .true.)

    end subroutine sf1026


!******************************************
!******************************************

!            f584

!******************************************
!******************************************

! subroutine to calculate the nighttime flux of
! he i (584) (note: line = 2)

    subroutine sf584 ( f, line, nfl, nll )

    use parameter_mod
    use variable_mod
    use photo_production_mod
    use grid_mod

    implicit none

    ! Arguments
    real,    intent(inout) :: f(nz,nf,nl,91)
    integer, intent(in)    :: line, nfl, nll

    ! Local variables
    integer :: i, k, imax

    imax = 1

    ! Determine f for the 4 known values of theta
    do i = 1, nz
        if ( alts(i,nfl,nll) < zaltnt(line,1) ) then
            do k = 1, 4
                f(i,nfl,nll, int(thetant(line,k))+1-90) = 1.0
            enddo
        elseif ( zaltnt(line,1) <= alts(i,nfl,nll) .and. &
                 alts(i,nfl,nll) <= zaltnt(line,2) ) then
            f(i,nfl,nll, int(thetant(line,1))+1-90) = &
                1.85e5 * ( alts(i,nfl,nll) - 170. ) ** 1.20
            f(i,nfl,nll, int(thetant(line,2))+1-90) = &
                2.60e4 * ( alts(i,nfl,nll) - 170. ) ** 1.25
            f(i,nfl,nll, int(thetant(line,3))+1-90) = &
                2.60e3 * ( alts(i,nfl,nll) - 170. ) ** 1.20
            f(i,nfl,nll, int(thetant(line,4))+1-90) = &
                2.60e2 * ( alts(i,nfl,nll) - 170. ) ** 1.20
            imax = i
        else
            do k = 1, 4
                f(i,   nfl,nll, int(thetant(line,k))+1-90) = &
                f(imax,nfl,nll, int(thetant(line,k))+1-90)
            enddo
        endif
    enddo

    ! Interpolate to all values of theta (90 - 180)
    ! sf584 uses extrap_last=.false. (clamp at f=1 when theta=180)
    call sf_nighttime_flux(f, line, nfl, nll, .false.)

    end subroutine sf584


!******************************************
!******************************************

!            f304

!******************************************
!******************************************

! subroutine to calculate the nighttime flux of
! he ii (304) (note: line = 3)

    subroutine sf304 ( f, line, nfl, nll )

    use parameter_mod
    use variable_mod
    use photo_production_mod
    use grid_mod

    implicit none

    ! Arguments
    real,    intent(inout) :: f(nz,nf,nl,91)
    integer, intent(in)    :: line, nfl, nll

    ! Local variables
    integer :: i, k, imax

    imax = 1

    ! Determine f for the 4 known values of theta
    do i = 1, nz
        if ( alts(i,nfl,nll) < zaltnt(line,1) ) then
            do k = 1, 4
                f(i,nfl,nll, int(thetant(line,k))+1-90) = 1.0
            enddo
        elseif ( zaltnt(line,1) <= alts(i,nfl,nll) .and. &
                 alts(i,nfl,nll) <= zaltnt(line,2) ) then
            f(i,nfl,nll, int(thetant(line,1))+1-90) = &
                3.8e6 * tanh( (alts(i,nfl,nll) - 138.) / 80. )
            f(i,nfl,nll, int(thetant(line,2))+1-90) = &
                3.0e6 * tanh( (alts(i,nfl,nll) - 138.) / 80. )
            f(i,nfl,nll, int(thetant(line,3))+1-90) = &
                2.5e6 * tanh( (alts(i,nfl,nll) - 138.) / 80. )
            f(i,nfl,nll, int(thetant(line,4))+1-90) = &
                2.5e6 * tanh( (alts(i,nfl,nll) - 138.) / 80. )
            imax = i
        else
            do k = 1, 4
                f(i,   nfl,nll, int(thetant(line,k))+1-90) = &
                f(imax,nfl,nll, int(thetant(line,k))+1-90)
            enddo
        endif
    enddo

    ! Interpolate to all values of theta (90 - 180)
    ! sf304 uses extrap_last=.false. (clamp at f=1 when theta=180)
    call sf_nighttime_flux(f, line, nfl, nll, .false.)

    end subroutine sf304


!******************************************
!******************************************

!            f1216

!******************************************
!******************************************

!     subroutine to calculate the nighttime flux of
!     lyman alpha (1216) (note: line = 4)

    subroutine sf1216 ( f, line, nfl, nll )

    use parameter_mod
    use variable_mod
    use photo_production_mod
    use grid_mod

    implicit none

    ! Arguments
    real,    intent(inout) :: f(nz,nf,nl,91)
    integer, intent(in)    :: line, nfl, nll

    ! Local variables
    integer :: i, k, imax

    imax = 1

    ! Determine f for the 4 known values of theta
    do i = 1, nz
        if ( alts(i,nfl,nll) < zaltnt(line,1) ) then
            do k = 1, 4
                f(i,nfl,nll, int(thetant(line,k))+1-90) = 1.0
            enddo
        elseif ( zaltnt(line,1) <= alts(i,nfl,nll) .and. &
                 alts(i,nfl,nll) <= zaltnt(line,2) ) then
            f(i,nfl,nll, int(thetant(line,1))+1-90) = &
                1.2e10 * tanh( (alts(i,nfl,nll) - 80.) / 50. ) + 3.e9
            f(i,nfl,nll, int(thetant(line,2))+1-90) = &
                4.0e9  * tanh( (alts(i,nfl,nll) - 80.) / 50. ) + 1.e9
            f(i,nfl,nll, int(thetant(line,3))+1-90) = &
                2.0e9  * tanh( (alts(i,nfl,nll) - 65.) / 50. ) + 1.e8
            f(i,nfl,nll, int(thetant(line,4))+1-90) = &
                1.5e9  * tanh( (alts(i,nfl,nll) - 75.) / 50. ) + 1.e8
            imax = i
        else
            do k = 1, 4
                f(i,   nfl,nll, int(thetant(line,k))+1-90) = &
                f(imax,nfl,nll, int(thetant(line,k))+1-90)
            enddo
        endif
    enddo

    ! Interpolate to all values of theta (90 - 180)
    ! sf1216 uses extrap_last=.true. (jip1 = ji+1 always, no endpoint clamping)
    call sf_nighttime_flux(f, line, nfl, nll, .true.)

    end subroutine sf1216


!******************************************
!******************************************

!            zenith

!******************************************
!******************************************

    subroutine zenith (hrut, nfl, nll, iday)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use photo_production_mod
    use grid_mod

    implicit none

    ! Arguments
    real,    intent(in) :: hrut
    integer, intent(in) :: nfl, nll, iday

    ! Local variables
    integer :: i
    real    :: hrl, sdec, cossdec, sinsdec, clat, slat

! bdec: magnetic declination angle
! sdec: solar declination angle
! cx:   cos of the solar zenith angle

    do i = 1, nz
        hrl = mod(hrut + glons(i,nfl,nll) / 15., 24.)
        sdec    = rtod * asin( sin(2.*pie*(nday(iday)-dayve)/sidyr) &
                             * sin(solinc/rtod) )
        cossdec = cos(po180 * sdec)
        sinsdec = sin(po180 * sdec)
        clat    = cos(po180 * glats(i,nfl,nll))
        slat    = sin(po180 * glats(i,nfl,nll))
        cx(i,nfl,nll) = slat * sinsdec &
                      - clat * cossdec * cos(15.0*po180*hrl)
        ! MS: Guard against round-off error minutely exceeding |cx|=1
        !     before taking acos in photprod.
        if (abs(abs(cx(i,nfl,nll)) - 1.) < 1.e-6) &
            cx(i,nfl,nll) = sign(1., cx(i,nfl,nll))
    enddo

    end subroutine zenith
