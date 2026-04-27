! hardy.f90

!      precipitation ionization rates

!      Rees, Physics and chemistry of the upper atmosphere
!            Cambridge Press, 1989
!            pp. 39 - 41

!      Rees, Auroral ionization and excitation by
!            incident energetic electrons,
!            Planet. Space Sci. 11, 1209, 1963

!   particle flux and average electron energy from Hardy model

    subroutine hardy_sami3(hrut, nfl, nll, iday)

    use hardy_mod
    use namelist_mod
    use grid_mod
    use parameter_mod

    implicit none

    ! Arguments
    real,    intent(in) :: hrut
    integer, intent(in) :: nfl, nll, iday

    ! Local variables
    integer :: i, j, k, ii, ngrids, is, iit, nglam, ngridn, in, nzh
    real    :: dele0, aven, value, part_flux
    real    :: e0ps, rps, r0s
    real    :: e0pn, rpn, r0n
    real    :: thelat, themlt

    ! File units for one-time reads
    integer :: lu1, lu2

    ! ------------------------------------------------------------------
    ! Read atmospheric profile tables on first call only.
    ! iread_hardy lives in hardy_mod with initialisation = 1.
    ! ------------------------------------------------------------------
    if ( iread_hardy == 1 ) then
        open(newunit=lu1, file='lambda_zor.inp')
        open(newunit=lu2, file='atm_params_int.inp')
        do ii = 1, no_lamzor
            read(lu1,*) zor(ii), lambda(ii)
        enddo
        do ii = 1, no_atm_ps
            read(lu2,*) height(ii), z(ii), rho(ii), nNn(ii), nM(ii)
        enddo
        close(lu1)
        close(lu2)
        iread_hardy = 0
    endif

    ! ------------------------------------------------------------------
    ! Initialise all ionization rates to 0
    ! ------------------------------------------------------------------
    preciprs(:,:,:) = 0.0
    preciprn(:,:,:) = 0.0
    iin(:,:) = nz
    iis(:,:) = 1

    ! ------------------------------------------------------------------
    ! Some parameters
    ! ------------------------------------------------------------------
    dele0 = 35.0e-3
    nzh   = nz / 2
    j     = nfl
    k     = nll

    ! ------------------------------------------------------------------
    ! South pole
    ! ------------------------------------------------------------------
    thelat = blats(1,j,k)
    themlt = mod(hrut + blons(nzh,nfl,nll) / 15.0, 24.0)

    call eavekp(kp(iday), thelat, themlt, aven)
    call elekp(1, kp(iday), thelat, themlt, value)
    part_flux = 10.0 ** value          ! explicit real base (was: 10 ** value)

    e0ps = aven
    rps  = 4.57e-6 * e0ps ** 1.75

    ! Search for the atmospheric grid level corresponding to rps.
    ! Guard against rps exceeding the top of the z table.
    ii = 1
    do while ( ii < no_atm_ps .and. z(ii) < rps )
        ii = ii + 1
    enddo
    ngrids = ii
    r0s    = rps / rho(ngrids)

    is = 1
    do while ( is < nzh .and. alts(is,j,k) <= 300.0 )
        iis(j,k) = is
        if ( alts(is,j,k) < height(ngrids) ) then
            preciprs(is,j,k) = 0.0
        else
            iit = ngrids
            do while ( iit > 1 .and. alts(is,j,k) > height(iit) )
                iit = iit - 1
            enddo
            if ( iit < ngrids .and. rps /= 0.0 ) then
                nglam = int( ( z(iit) / rps ) * no_lamzor ) + 1
                ! Clamp nglam to valid lambda index range
                nglam = max(1, min(nglam, no_lamzor))
                preciprs(is,j,k) = ( e0ps / r0s ) / dele0 * &
                                   lambda(nglam) * part_flux * &
                                   nM(iit) / nM(ngrids)
            endif
        endif
        is = is + 1
    enddo

    ! ------------------------------------------------------------------
    ! North pole
    ! ------------------------------------------------------------------
    thelat = blats(nz,j,k)
    themlt = mod(hrut + blons(nzh,nfl,nll) / 15.0, 24.0)

    call eavekp(kp(iday), thelat, themlt, aven)
    call elekp(1, kp(iday), thelat, themlt, value)
    part_flux = 10.0 ** value          ! explicit real base (was: 10 ** value)

    e0pn = aven
    rpn  = 4.57e-6 * e0pn ** 1.75

    ! Search for the atmospheric grid level corresponding to rpn.
    ! Guard against rpn exceeding the top of the z table.
    ii = 1
    do while ( ii < no_atm_ps .and. z(ii) < rpn )
        ii = ii + 1
    enddo
    ngridn = ii
    r0n    = rpn / rho(ngridn)

    in = nz
    do while ( in > nzh .and. alts(in,j,k) <= 300.0 )
        iin(j,k) = in
        if ( alts(in,j,k) < height(ngridn) ) then
            preciprn(in,j,k) = 0.0
        else
            iit = ngridn
            do while ( iit > 1 .and. alts(in,j,k) > height(iit) )
                iit = iit - 1
            enddo
            if ( iit < ngridn .and. rpn /= 0.0 ) then
                nglam = int( ( z(iit) / rpn ) * no_lamzor ) + 1
                ! Clamp nglam to valid lambda index range
                nglam = max(1, min(nglam, no_lamzor))
                preciprn(in,j,k) = ( e0pn / r0n ) / dele0 * &
                                   lambda(nglam) * part_flux * &
                                   nM(iit) / nM(ngridn)
            endif
        endif
        in = in - 1
    enddo

    end subroutine hardy_sami3


    subroutine eavekp(ikp, thelat, themlt, aven)

!******************************************************************************
!
!   INPUTS:
!     ikp          integral Kp value (0 through 6)
!     thelat       corrected geomagnetic latitude
!     themlt       magnetic local time
!
!   OUTPUTS:
!     aven         the electron average energy for this Kp map (keV)
!
!   FILES:
!     eavekp.inp   ASCII file containing coefficients for the
!                  Chebyshev-Fourier expansion
!
!   NOTES:
!     Routine is valid above corrected geomagnetic latitude of 50 deg.
!
!   REVISIONS:
!     16 SEP 98  original
!
!******************************************************************************

    use hardy_mod

    implicit none

    ! Arguments
    integer, intent(in)  :: ikp
    real,    intent(in)  :: thelat, themlt
    real,    intent(out) :: aven

    ! Local variables
    integer :: i, j, k, imod, jmod, ll1, mm1
    integer :: lu_eave
    real    :: ttlat, twothirds, xx, sig
    real    :: av, todeg, degr, fmm
    real    :: y(ll_eave)

    ! ------------------------------------------------------------------
    ! Read coefficient file on first call only.
    ! ifirst_eave lives in hardy_mod with initialisation = 1.
    ! ------------------------------------------------------------------
    if ( ifirst_eave == 1 ) then
        ifirst_eave = 0
        open(newunit=lu_eave, file='eavekp.inp', form='formatted')
        do i = 1, 7
            read(lu_eave, *, err=98) jmod, ll1, mm1
            if ( ll1 /= ll_eave ) then
                print *, 'EAVEKP: not proper L value'
                stop
            endif
            if ( mm1 /= mm_eave ) then
                print *, 'EAVEKP: not proper M value'
                stop
            endif
            read(lu_eave, *, err=98) (eave_a(jmod,j), j=1,ll_eave),    &
            ((eave_b(jmod,j,k), j=1,ll_eave), k=1,mm_eave), &
            ((eave_c(jmod,j,k), j=1,ll_eave), k=1,mm_eave)
        enddo
        close(lu_eave)
    endif

    ! ------------------------------------------------------------------
    ! Validate Kp
    ! ------------------------------------------------------------------
    if ( ikp < 0 ) then
        print *, 'EAVEKP: kp must be >= 0'
        stop
    endif
    if ( ikp > 6 ) then
        print *, 'EAVEKP: kp must be <= 6'
        stop
    endif

    ! Model number is ikp+1
    imod = ikp + 1

    ! ------------------------------------------------------------------
    ! Get x value from latitude
    ! ------------------------------------------------------------------
    ttlat = abs(thelat)
    ttlat = max(ttlat, 50.0)
    ttlat = min(ttlat, 89.0)

    twothirds = 2.0 / 3.0
    xx = 2.0 * (ttlat - 50.0) / 40.0 - 1.0
    if ( xx < 0.0 ) then
        sig = -1.0
    else
        sig =  1.0
    endif
    xx = sig * abs(xx) ** twothirds

    ! ------------------------------------------------------------------
    ! Evaluate Chebyshev polynomials
    ! ------------------------------------------------------------------
    call getchevy3(xx, ll_eave, y)

    ! ------------------------------------------------------------------
    ! Accumulate Chebyshev-Fourier sum
    ! ------------------------------------------------------------------
    av    = 0.0
    todeg = 45.0 / atan(1.0)
    degr  = themlt * 360.0 / 24.0 / todeg

    do i = 1, ll_eave
        av = av + y(i) * eave_a(imod,i)
        do j = 1, mm_eave
            fmm = real(j)
            av  = av + y(i) * eave_b(imod,i,j) * sin(fmm*degr)
            av  = av + y(i) * eave_c(imod,i,j) * cos(fmm*degr)
        enddo
    enddo

    aven = exp(av)

    return

    98 print *, 'EAVEKP: file eavekp.inp is corrupted'
    stop

    end subroutine eavekp


    subroutine getchevy3(x, iord, chevy)

!   Evaluates the direct expansion of the Chebyshev polynomials.

    implicit none

    ! Arguments
    real,    intent(in)  :: x
    integer, intent(in)  :: iord
    real,    intent(out) :: chevy(iord)   ! sized to actual order, not fixed 20

    ! Local variables
    integer :: j

    if ( iord > 20 ) then
        print *, 'GETCHEVY3: not more than 20 terms'
        stop
    endif

    chevy(1) = 1.0
    if ( iord >= 2 ) chevy(2) = x
    do j = 3, iord
        chevy(j) = 2.0 * x * chevy(j-1) - chevy(j-2)
    enddo

    end subroutine getchevy3


    subroutine elekp(iopt, kp, thelat, themlt, value)

!   INPUTS:      iopt    1=electron number flux
!                        2=electron energy flux
!
!                kp      the Kp map number (0 through 6)
!
!                thelat  corrected geomagnetic latitude (deg)
!                themlt  magnetic local time (0->24)
!
!   OUTPUT:      value   log10 of the flux (1/cm^2/sec)
!
!   NOTES:       The poleward values of the fluxes are set to a minimum
!                of the values in prain(), derived from data.
!                Equatorward values are set at unity.
!
!   REVISIONS:   11 SEP 98  original

    use hardy_mod

    implicit none

    ! Arguments
    integer, intent(in)  :: iopt, kp
    real,    intent(in)  :: thelat, themlt
    real,    intent(out) :: value

    ! Local variables
    integer :: i, k, ikp, jopt, jkp, jcoef, jnum
    integer :: lu_elekp
    real    :: pi, arg, aarg, cosarg, sinarg
    real    :: r0, h0, h1, s0, r1, s2, s1
    real    :: b1, b2, h, eh
    character(80) :: line

    ! ------------------------------------------------------------------
    ! Read coefficient file on first call only.
    ! ifirst_elekp lives in hardy_mod with initialisation = 0.
    ! ------------------------------------------------------------------
    if ( ifirst_elekp == 0 ) then
        ifirst_elekp = 1
        open(newunit=lu_elekp, file='elekp.inp', form='formatted')
        do jopt = 1, 2
            do jkp = 1, 7
                read(lu_elekp, '(a80)', err=98) line
                do jcoef = 1, 17
                    read(lu_elekp, *, err=98) jnum,              &
                        elekp_cr0(jopt,jkp,jcoef),               &
                        elekp_ch0(jopt,jkp,jcoef),               &
                        elekp_cr1(jopt,jkp,jcoef),               &
                        elekp_ch1(jopt,jkp,jcoef),               &
                        elekp_cs0(jopt,jkp,jcoef),               &
                        elekp_cs2(jopt,jkp,jcoef)
                enddo
            enddo
        enddo
        close(lu_elekp)
    endif

    ! ------------------------------------------------------------------
    ! Validate Kp before use
    ! ------------------------------------------------------------------
    if ( kp < 0 ) then
        print *, 'ELEKP: kp must be >= 0'
        stop
    endif
    if ( kp > 6 ) then
        print *, 'ELEKP: kp must be <= 6'
        stop
    endif

    ikp = kp + 1

    ! ------------------------------------------------------------------
    ! Accumulate Fourier sum over local time
    ! ------------------------------------------------------------------
    pi  = 4.0 * atan(1.0)
    r0  = elekp_cr0(iopt,ikp,1) / 2.0
    h0  = elekp_ch0(iopt,ikp,1) / 2.0
    h1  = elekp_ch1(iopt,ikp,1) / 2.0
    s0  = elekp_cs0(iopt,ikp,1) / 2.0
    r1  = elekp_cr1(iopt,ikp,1) / 2.0
    s2  = elekp_cs2(iopt,ikp,1) / 2.0

    arg = pi * themlt / 12.0
    k   = 2
    do i = 1, 8
        aarg   = arg * i
        cosarg = cos(aarg)
        sinarg = sin(aarg)
        r0 = r0 + cosarg * elekp_cr0(iopt,ikp,k)
        h0 = h0 + cosarg * elekp_ch0(iopt,ikp,k)
        h1 = h1 + cosarg * elekp_ch1(iopt,ikp,k)
        s0 = s0 + cosarg * elekp_cs0(iopt,ikp,k)
        r1 = r1 + cosarg * elekp_cr1(iopt,ikp,k)
        s2 = s2 + cosarg * elekp_cs2(iopt,ikp,k)
        r0 = r0 + sinarg * elekp_cr0(iopt,ikp,k+8)
        h0 = h0 + sinarg * elekp_ch0(iopt,ikp,k+8)
        h1 = h1 + sinarg * elekp_ch1(iopt,ikp,k+8)
        s0 = s0 + sinarg * elekp_cs0(iopt,ikp,k+8)
        r1 = r1 + sinarg * elekp_cr1(iopt,ikp,k+8)
        s2 = s2 + sinarg * elekp_cs2(iopt,ikp,k+8)
        k  = k + 1
    enddo

    b1 = log((1.0 + exp(h1-h0)) / 2.0)
    b2 = log(2.0 / (1.0 + exp(h0-h1)))
    s1 = (r1 - r0 - s0*(h1-h0) + s0*b1 - s2*b2) / (b1 - b2)

    ! ------------------------------------------------------------------
    ! Evaluate latitude profile
    ! ------------------------------------------------------------------
    h = abs(thelat)
    if ( h < 50.0 ) h = 50.0

    eh = r0 + s0*(h - h0)
    eh = eh + (s1 - s0) * log((1.0 + exp(h-h0)) / 2.0)
    eh = eh + (s2 - s1) * log((1.0 + exp(h-h1)) / &
                               (1.0 + exp(h0-h1)))

    if ( eh < 0.0 ) eh = 0.0

    if ( h > h1 .and. eh < prain(iopt,ikp) ) then
        eh = prain(iopt,ikp)
    endif

    value = eh

    return

    98 print *, 'ELEKP: file elekp.inp is corrupted'
    stop

    end subroutine elekp
