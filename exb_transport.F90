!******************************************
!******************************************

!             EXB

!******************************************
!******************************************

    subroutine exb(hrut, phi)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use time_mod
    use exb_mod
    use misc_mod
    use grid_mod

    implicit none

    ! Arguments
    real,    intent(in) :: hrut
    real,    intent(in) :: phi(nnx,nny)

    ! Local variables
    integer :: i, j, k, ni
    real    :: arg0, fac
    real    :: alt_crit_high, dela_high, dela

! define the e x b drift

    call vexb_phi(phi)

! here we add in exb drift from the potential

! reduce e x b velocities below alt_crit
! with 20 km decay (dela)

! and kill above some high altitude

    alt_crit_high = pcrit * re
    dela_high     = 2.0 * re
    dela          = 20.0

!     vexbp

    do k = 1, nl
        do j = 1, nfp1
            do i = 1, nz
                vexbp(i,j,k) = vexbp_phi(i,j,k)
                if ( baltp(i,j,k) < alt_crit ) then
                    arg0 = ( alt_crit - baltp(i,j,k) ) / dela
                    fac  = exp(-arg0*arg0)
                    vexbp(i,j,k) = vexbp(i,j,k) * fac
                endif
                if ( baltp(i,j,k) > alt_crit_high ) then
                    arg0 = abs( baltp(i,j,k) - alt_crit_high ) / dela_high
                    fac  = exp(-arg0*arg0)
                    vexbp(i,j,k) = vexbp(i,j,k) * fac
                endif
            enddo
        enddo
    enddo

!     vexbs

    do k = 1, nl
        do j = 1, nf
            do i = 1, nzp1
                vexbs(i,j,k) = vexbs_phi(i,j,k)
                if ( baltp(i,j,k) < alt_crit ) then
                    arg0 = ( alt_crit - baltp(i,j,k) ) / dela
                    fac  = exp(-arg0*arg0)
                    vexbs(i,j,k) = vexbs(i,j,k) * fac
                endif
                if ( baltp(i,j,k) > alt_crit_high ) then
                    arg0 = abs( baltp(i,j,k) - alt_crit_high ) / dela_high
                    fac  = exp(-arg0*arg0)
                    vexbs(i,j,k) = vexbs(i,j,k) * fac
                endif
            enddo
        enddo
    enddo

!     vexbh

    do k = 1, nlp1
        do j = 1, nf
            do i = 1, nz
                vexbh(i,j,k) = vexbh_phi(i,j,k)
                if ( baltp(i,j,k) < alt_crit ) then
                    arg0 = ( alt_crit - baltp(i,j,k) ) / dela
                    fac  = exp(-arg0*arg0)
                    vexbh(i,j,k) = vexbh(i,j,k) * fac
                endif
                if ( baltp(i,j,k) > alt_crit_high ) then
                    arg0 = abs( baltp(i,j,k) - alt_crit_high ) / dela_high
                    fac  = exp(-arg0*arg0)
                    vexbh(i,j,k) = vexbh(i,j,k) * fac
                endif
            enddo
        enddo
    enddo

! limit e x b velocities

    do k = 1, nl
        do j = 1, nfp1
            do i = 1, nz
                if (vexbp(i,j,k) > 0.0) &
                vexbp(i,j,k) = min(vexbp(i,j,k),  vexb_max)
                if (vexbp(i,j,k) < 0.0) &
                vexbp(i,j,k) = max(vexbp(i,j,k), -vexb_max)
            enddo
        enddo
    enddo

    do k = 1, nl
        do j = 1, nf
            do i = 1, nzp1
                if (vexbs(i,j,k) > 0.0) &
                vexbs(i,j,k) = min(vexbs(i,j,k),  vexb_max)
                if (vexbs(i,j,k) < 0.0) &
                vexbs(i,j,k) = max(vexbs(i,j,k), -vexb_max)
            enddo
        enddo
    enddo

    do k = 1, nlp1
        do j = 1, nf
            do i = 1, nz
                if (vexbh(i,j,k) > 0.0) &
                vexbh(i,j,k) = min(vexbh(i,j,k),  vexb_max)
                if (vexbh(i,j,k) < 0.0) &
                vexbh(i,j,k) = max(vexbh(i,j,k), -vexb_max)
            enddo
        enddo
    enddo

! output e x b velocities
! NOTE: u1p/u2s/u3h are cell-centred (nz,nf,nl); only interior
!       values are copied. Face values at nfp1/nzp1/nlp1 are not
!       stored here by design.

    do k = 1, nl
        do j = 1, nf
            do i = 1, nz
                u1p(i,j,k) = vexbp(i,j,k)
                u2s(i,j,k) = vexbs(i,j,k)
                u3h(i,j,k) = vexbh(i,j,k)
            enddo
        enddo
    enddo

! calculate conserved particle number: denic
! and 'conserved' temperature: tic, tec

    do ni = nion1, nion2
        do k = 1, nl
            do j = 1, nf
                do i = 1, nz
                    denic(i,j,k,ni) = deni(i,j,k,ni) * vol(i,j,k)
                    tic(i,j,k,ni)   = ti(i,j,k,ni)   * vol(i,j,k)
                enddo
            enddo
        enddo
    enddo

    do k = 1, nl
        do j = 1, nf
            do i = 1, nz
                tec(i,j,k) = te(i,j,k) * vol(i,j,k)
            enddo
        enddo
    enddo

! calculate flux in p-direction at interface
! NOTE: neutral flux condition at outer boundary (JH 11/29/07)
! altered to consider NS pole densities

    do ni = nion1, nion2
        do k = 1, nl
            do j = 2, nf
                do i = 1, nz
                    if ( vexbp(i,j,k) >= 0 ) then
                        fluxnp(i,j,k,ni) = deni(i,j-1,k,ni) * vexbp(i,j,k)
                        fluxtp(i,j,k,ni) = ti(i,j-1,k,ni)   * vexbp(i,j,k)
                    else
                        fluxnp(i,j,k,ni) = deni(i,j,k,ni) * vexbp(i,j,k)
                        fluxtp(i,j,k,ni) = ti(i,j,k,ni)   * vexbp(i,j,k)
                    endif
                enddo
            enddo
        enddo
    enddo

    do k = 1, nl
        do j = 2, nf
            do i = 1, nz
                if ( vexbp(i,j,k) >= 0 ) then
                    fluxtep(i,j,k) = te(i,j-1,k) * vexbp(i,j,k)
                else
                    fluxtep(i,j,k) = te(i,j,k)   * vexbp(i,j,k)
                endif
            enddo
        enddo
    enddo

!      flux at nfp1 (near magnetic north/south poles)

    do ni = nion1, nion2
        do k = 1, nl
            j = nfp1
            do i = 1, nz
                if ( vexbp(i,j,k) >= 0 ) then
                    fluxnp(i,j,k,ni) = deni(i,j-1,k,ni) * vexbp(i,j,k)
                    fluxtp(i,j,k,ni) = ti(i,j-1,k,ni)   * vexbp(i,j,k)
                else
                    fluxnp(i,j,k,ni) = deni_mnp(i,ni) * vexbp(i,j,k)
                    fluxtp(i,j,k,ni) = ti_mnp(i,ni)   * vexbp(i,j,k)
                endif
            enddo
        enddo
    enddo

    do k = 1, nl
        j = nfp1
        do i = 1, nz
            if ( vexbp(i,j,k) >= 0 ) then
                fluxtep(i,j,k) = te(i,j-1,k) * vexbp(i,j,k)
            else
                fluxtep(i,j,k) = te_mnp(i)   * vexbp(i,j,k)
            endif
        enddo
    enddo

! calculate flux in s-direction at interface

    do ni = nion1, nion2
        do k = 1, nl
            do j = 1, nf
                do i = 2, nz
                    if ( vexbs(i,j,k) >= 0 ) then
                        fluxns(i,j,k,ni) = deni(i-1,j,k,ni) * vexbs(i,j,k)
                        fluxts(i,j,k,ni) = ti(i-1,j,k,ni)   * vexbs(i,j,k)
                    else
                        fluxns(i,j,k,ni) = deni(i,j,k,ni) * vexbs(i,j,k)
                        fluxts(i,j,k,ni) = ti(i,j,k,ni)   * vexbs(i,j,k)
                    endif
                enddo
            enddo
        enddo
    enddo

    do k = 1, nl
        do j = 1, nf
            do i = 2, nz
                if ( vexbs(i,j,k) >= 0 ) then
                    fluxtes(i,j,k) = te(i-1,j,k) * vexbs(i,j,k)
                else
                    fluxtes(i,j,k) = te(i,j,k)   * vexbs(i,j,k)
                endif
            enddo
        enddo
    enddo

! calculate flux in h-direction at interface (k > 1)

    do ni = nion1, nion2
        do k = 2, nl
            do j = 1, nf
                do i = 1, nz
                    if ( vexbh(i,j,k) >= 0 ) then
                        fluxnh(i,j,k,ni) = deni(i,j,k-1,ni) * vexbh(i,j,k)
                        fluxth(i,j,k,ni) = ti(i,j,k-1,ni)   * vexbh(i,j,k)
                    else
                        fluxnh(i,j,k,ni) = deni(i,j,k,ni) * vexbh(i,j,k)
                        fluxth(i,j,k,ni) = ti(i,j,k,ni)   * vexbh(i,j,k)
                    endif
                enddo
            enddo
        enddo
    enddo

    do k = 2, nl
        do j = 1, nf
            do i = 1, nz
                if ( vexbh(i,j,k) >= 0 ) then
                    fluxteh(i,j,k) = te(i,j,k-1) * vexbh(i,j,k)
                else
                    fluxteh(i,j,k) = te(i,j,k)   * vexbh(i,j,k)
                endif
            enddo
        enddo
    enddo

! update total particle number and density
! and temperatures
! NOTE: the temperature update is an approximation
!       (probably better than no update but, strictly
!       speaking, not exactly correct)

    do ni = nion1, nion2
        do k = 2, nl-1
            do j = 2, nf
                do i = 2, nz-1
                    denic(i,j,k,ni) = denic(i,j,k,ni) &
                    + dt * ( areap(i,j,k)   * fluxnp(i,j,k,ni) - &
                    areap(i,j+1,k) * fluxnp(i,j+1,k,ni) ) &
                    + dt * ( areas(i,j,k)   * fluxns(i,j,k,ni) - &
                    areas(i+1,j,k) * fluxns(i+1,j,k,ni) ) &
                    + dt * ( areah(i,j,k)   * fluxnh(i,j,k,ni) - &
                    areah(i,j,k+1) * fluxnh(i,j,k+1,ni) )
                    deni(i,j,k,ni)  = denic(i,j,k,ni) / vol(i,j,k)

                ! brazen fix
                    deni(i,j,k,ni)  = max(deni(i,j,k,ni), denmin)

                    tic(i,j,k,ni) = tic(i,j,k,ni) &
                    + dt * ( areap(i,j,k)   * fluxtp(i,j,k,ni) - &
                    areap(i,j+1,k) * fluxtp(i,j+1,k,ni) ) &
                    + dt * ( areas(i,j,k)   * fluxts(i,j,k,ni) - &
                    areas(i+1,j,k) * fluxts(i+1,j,k,ni) ) &
                    + dt * ( areah(i,j,k)   * fluxth(i,j,k,ni) - &
                    areah(i,j,k+1) * fluxth(i,j,k+1,ni) )
                    ti(i,j,k,ni)  = tic(i,j,k,ni) / vol(i,j,k)
                ! brazen fix
                    ti(i,j,k,ni)  = max(ti(i,j,k,ni), 200.0)
                    if (isnan(ti(i,j,k,ni))) then
                        print *, 'Ti fixed', i, j, k, ni
                        ti(i,j,k,ni) = 200.0
                    endif
                enddo
            enddo
        enddo
    enddo

    do k = 2, nl-1
        do j = 2, nf
            do i = 2, nzm1
                tec(i,j,k) = tec(i,j,k) &
                + dt * ( areap(i,j,k)   * fluxtep(i,j,k) - &
                areap(i,j+1,k) * fluxtep(i,j+1,k) ) &
                + dt * ( areas(i,j,k)   * fluxtes(i,j,k) - &
                areas(i+1,j,k) * fluxtes(i+1,j,k) ) &
                + dt * ( areah(i,j,k)   * fluxteh(i,j,k) - &
                areah(i,j,k+1) * fluxteh(i,j,k+1) )
                te(i,j,k)  = tec(i,j,k) / vol(i,j,k)
            ! brazen fix
                te(i,j,k)  = max(te(i,j,k), 200.0)
                if (te(i,j,k) < 0.0) print *, 'i,j,k', &
                i, j, k, te(i,j,k)
            enddo
        enddo
    enddo

! fill cells at j = 1 with j = 2

    do ni = nion1, nion2
        do k = 2, nl-1
            do i = 2, nzm1
                deni(i,1,k,ni) = deni(i,2,k,ni)
                ti(i,1,k,ni)   = ti(i,2,k,ni)
            enddo
        enddo
    enddo

    do k = 2, nl-1
        do i = 2, nzm1
            te(i,1,k) = te(i,2,k)
        enddo
    enddo

! fill cells at i = 1 and nz with i = 2 and nzm1

    do ni = nion1, nion2
        do k = 2, nl-1
            do j = 1, nf
                deni(1,j,k,ni)  = deni(2,j,k,ni)
                deni(nz,j,k,ni) = deni(nzm1,j,k,ni)
                ti(1,j,k,ni)    = ti(2,j,k,ni)
                ti(nz,j,k,ni)   = ti(nzm1,j,k,ni)
            enddo
        enddo
    enddo

    do k = 2, nl-1
        do j = 1, nf
            te(1,j,k)  = te(2,j,k)
            te(nz,j,k) = te(nzm1,j,k)
        enddo
    enddo

    end subroutine exb



!     ********************************************

!     VEXB_PHI

!     ********************************************

    subroutine vexb_phi(phi)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use exb_mod
    use grid_mod

    implicit none

    ! Arguments
    real, intent(in) :: phi(nnx,nny)

    ! Local variables
    integer :: i, j, k
    real    :: sinangx, sinangy, sinangz, sinang
    real    :: vps, vhs, slope
    real    :: phihp(nfp1,nlp1)
    real    :: phish(nfp1,nl)
    real    :: phisp(nf,nlp1)

!      p face

    call phi_nfp1_nlp1(bradss, blonss, phihp, phi)

    do k = 1, nl
        do j = 1, nfp1
            do i = 1, nz
                sinangx =  bdirpy(i,j,k) * ehpz(i,j,k) - &
                bdirpz(i,j,k) * ehpy(i,j,k)
                sinangy = -bdirpx(i,j,k) * ehpz(i,j,k) + &
                bdirpz(i,j,k) * ehpx(i,j,k)
                sinangz =  bdirpx(i,j,k) * ehpy(i,j,k) - &
                bdirpy(i,j,k) * ehpx(i,j,k)
                sinang  =  sqrt( sinangx * sinangx + &
                sinangy * sinangy + &
                sinangz * sinangz  )
                ehp(i,j,k) = -1.0 * &
                ( phihp(j,k+1) - phihp(j,k) ) / delhp(i,j,k) / &
                sinang
            enddo
        enddo
    enddo

    do k = 1, nl
        do j = 1, nfp1
            do i = 1, nz
                vexbp_phi(i,j,k) = ehp(i,j,k) * sol / bmag / bmpf(i,j,k) &
                * tvexb0
            enddo
        enddo
    enddo

!      h face

    do k = 1, nlp1
        do j = 1, nf
            do i = 1, nz
                sinangx =  bdirhy(i,j,k) * ephz(i,j,k) - &
                bdirhz(i,j,k) * ephy(i,j,k)
                sinangy = -bdirhx(i,j,k) * ephz(i,j,k) + &
                bdirhz(i,j,k) * ephx(i,j,k)
                sinangz =  bdirhx(i,j,k) * ephy(i,j,k) - &
                bdirhy(i,j,k) * ephx(i,j,k)
                sinang  =  sqrt( sinangx * sinangx + &
                sinangy * sinangy + &
                sinangz * sinangz  )
                eph(i,j,k) = -1.0 * &
                ( phihp(j+1,k) - phihp(j,k) ) / &
                delph(i,j,k) / sinang
            enddo
        enddo
    enddo

    do k = 1, nlp1
        do j = 1, nf-1
            do i = 1, nz
                vexbh_phi(i,j,k) = -eph(i,j,k) * sol / bmag / bmhf(i,j,k) &
                * tvexb0
            enddo
        enddo
    enddo

!   fix at j = nf (extrapolate except at nz/2+1 - then interpolate)
!   fix re doug drob

    j = nf

! north pole

    do k = 1, nlp1
        do i = nz/2, nz
            slope = ( blatp(i,j,k) - blatp(i,j-1,k) ) / &
                    ( 90.0          - blatp(i,j-1,k) )
            vexbh_phi(i,j,k) = vexbh_phi(i,j-1,k) * (1.0 - slope) &
            * tvexb0
        enddo
    enddo

! south pole

    do k = 1, nlp1
        do i = 1, nz/2-1
            slope = ( blatp(i,j,k)  - blatp(i,j-1,k) ) / &
                    ( -90.0          - blatp(i,j-1,k) )
            vexbh_phi(i,j,k) = vexbh_phi(i,j-1,k) * (1.0 - slope) &
            * tvexb0
        enddo
    enddo

!      s face

    call phi_nfp1_nl (bradsh, blonsh, phish, phi)
    call phi_nf_nlp1 (bradsp, blonsp, phisp, phi)

    do k = 1, nl
        do j = 1, nf
            do i = 1, nzp1
                eps(i,j,k) = -1.0 * &
                ( phish(j+1,k) - phish(j,k) ) / delps(i,j,k)
            enddo
        enddo
    enddo

    do k = 1, nl
        do j = 1, nf
            do i = 1, nzp1
                ehs(i,j,k) = -1.0 * &
                ( phisp(j,k+1) - phisp(j,k) ) / delhs(i,j,k)
            enddo
        enddo
    enddo

    do k = 1, nl
        do j = 1, nf
            do i = 1, nz
                vps = vexbp_phi(i,j,k) * &
                ( vpsnx(i,j,k) * xnorms(i,j,k) + &
                vpsny(i,j,k) * ynorms(i,j,k) + &
                vpsnz(i,j,k) * znorms(i,j,k)   )

                vhs = vexbh_phi(i,j,k) * &
                ( vhsnx(i,j,k) * xnorms(i,j,k) + &
                vhsny(i,j,k) * ynorms(i,j,k) + &
                vhsnz(i,j,k) * znorms(i,j,k)   )

                vexbs_phi(i,j,k) = (vps + vhs) * tvexb0
            enddo
        enddo
    enddo

    do k = 1, nl
        do j = 1, nf
            vexbs_phi(nzp1,j,k) = vexbs_phi(nz,j,k)
        enddo
    enddo

    end subroutine vexb_phi



!     ********************************************

!     PHI_NFP1_NLP1

!     ********************************************


! takes a variable in one grid (xu,zu):(nf,nl)
! and interpolates to another grid (x0,z0):(nx,ny)
! the nf,nf grid need not be orthogonal

    subroutine phi_nfp1_nlp1(brad, blon, phiout, phiin)

    use parameter_mod
    use variable_mod
    use message_passing_mod
    use grid_mod

    implicit none

    ! Arguments
    real, intent(in)  :: brad(nfp1,nlp1), blon(nfp1,nlp1)
    real, intent(out) :: phiout(nfp1,nlp1)
    real, intent(in)  :: phiin(nnx,nny)

    ! Local variables
    integer :: i, j, k, kk, jj
    real    :: slope

    phiout(:,:) = 0.0

    if ( taskid /= 1 .AND. taskid /= numworkers ) then
        do k = 1, nlp1
            do j = nf, 2, -1
                kk = (taskid-1)*(nl-2) + (k-1)
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == 1 ) then
        do k = 1, nlp1
            do j = nf, 2, -1
                kk = k - 1
                if ( k == 1 ) kk = nnx - 1
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == numworkers ) then
        do k = 1, nlp1
            do j = nf, 2, -1
                kk = (taskid - 1)*(nl-2) + (k-1)
                if ( k == nlp1   ) kk = 2
                if ( k == nlp1-1 ) kk = 1
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    do k = 1, nlp1
        phiout(1,k)    = phiout(2,k)
        slope          = ( blatp(nz-1,nf+1,k) - blatp(nz-1,nf-1,k) ) / &
        ( blatp(nz-1,nf,k) - blatp(nz-1,nf-1,k) )
        phiout(nf+1,k) = phiout(nf-1,k) + &
        ( phiout(nf,k) - phiout(nf-1,k) ) * slope
    enddo

    end subroutine phi_nfp1_nlp1

!     ********************************************

!     PHI_NFP1_NL

!     ********************************************


! takes a variable in one grid (xu,zu):(nf,nl)
! and interpolates to another grid (x0,z0):(nx,ny)
! the nf,nf grid need not be orthogonal

    subroutine phi_nfp1_nl(brad, blon, phiout, phiin)

    use parameter_mod
    use variable_mod
    use message_passing_mod
    use grid_mod

    implicit none

    ! Arguments
    real, intent(in)  :: brad(nfp1,nl), blon(nfp1,nl)
    real, intent(out) :: phiout(nfp1,nl)
    real, intent(in)  :: phiin(nnx,nny)

    ! Local variables
    integer :: i, j, k, kk, jj
    real    :: slope

    phiout(:,:) = 0.0

    if ( taskid /= 1 .AND. taskid /= numworkers ) then
        do k = 1, nl
            do j = nf, 2, -1
                kk = (taskid-1)*(nl-2) + (k-1)
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == 1 ) then
        do k = 1, nl
            do j = nf, 2, -1
                kk = k - 1
                if ( k == 1 ) kk = nnx - 1
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == numworkers ) then
        do k = 1, nl
            do j = nf, 2, -1
                kk = (taskid - 1)*(nl-2) + (k-1)
                if ( k == nl ) kk = 1
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    do k = 1, nl
        phiout(1,k)    = phiout(2,k)
        slope          = ( blatp(nz-1,nf+1,k) - blatp(nz-1,nf-1,k) ) / &
        ( blatp(nz-1,nf,k) - blatp(nz-1,nf-1,k) )
        phiout(nf+1,k) = phiout(nf-1,k) + &
        ( phiout(nf,k) - phiout(nf-1,k) ) * slope
    enddo

    end subroutine phi_nfp1_nl



!     ********************************************

!     PHI_NF_NLP1

!     ********************************************


    subroutine phi_nf_nlp1(brad, blon, phiout, phiin)

    use parameter_mod
    use variable_mod
    use message_passing_mod
    use grid_mod

    implicit none

    ! Arguments
    real, intent(in)  :: brad(nf,nlp1), blon(nf,nlp1)
    real, intent(out) :: phiout(nf,nlp1)
    real, intent(in)  :: phiin(nnx,nny)

    ! Local variables
    integer :: i, j, k, kk, jj

    phiout(:,:) = 0.0

    if ( taskid /= 1 .AND. taskid /= numworkers ) then
        do k = 1, nlp1
            do j = nf-1, 2, -1
                kk = (taskid-1)*(nl-2) + (k-1)
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == 1 ) then
        do k = 1, nlp1
            do j = nf-1, 2, -1
                kk = k - 1
                if ( k == 1 ) kk = nnx - 1
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    if ( taskid == numworkers ) then
        do k = 1, nlp1
            do j = nf-1, 2, -1
                kk = (taskid - 1)*(nl-2) + (k-1)
                if ( k == nlp1   ) kk = 2
                if ( k == nlp1-1 ) kk = 1
                jj = j - 1
                phiout(j,k) = phiin(kk,jj)
            enddo
        enddo
    endif

    do k = 1, nlp1
        phiout(nf,k) = phiout(nf-1,k) * &
        ( blatp(nz-1,nf,k)   - 90.0 ) / &
        ( blatp(nz-1,nf-1,k) - 90.0 )
        phiout(1,k)  = phiout(2,k)
    enddo

    end subroutine phi_nf_nlp1
