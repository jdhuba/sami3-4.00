!       **********************************
!       **********************************

!             POTPPHI

!       **********************************
!       **********************************

    subroutine potpphi(phi,dphi,hrut)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use conductance_mod
    use misc_mod
    use grid_mod

    implicit none

    ! Arguments
    real,    intent(out)   :: phi(nnx,nny)
    real(8), intent(inout) :: dphi(nnx+1,nnyt)   ! initial guess in, solution out
    real,    intent(in)    :: hrut

    ! Named unit number (fix #10: was magic literal 1232)
    integer, parameter :: U_DPHI_RST = 1232

    ! Local arrays (working space for conductance interpolation)
    real :: hipcp_pot(nnx,nnyt),hipcphi_pot(nnx,nnyt)
    real :: hidphig_pot(nnx,nnyt),hidphiv_pot(nnx,nnyt)
    real :: hidpg_pot(nnx,nnyt),hidpv_pot(nnx,nnyt)
    real :: hihcm_pot(nnx,nnyt)
    real :: hipc_pot(nnx,nnyt)
    real :: hihc_pot(nnx,nnyt)
    real :: hidv_pot(nnx,nnyt)

    real    :: dpreal(nnyt),preal(nnyt)
    real    :: dbang(nnx,nnyt),blang(nnx)
    real(8) :: dphireal(nnx+1)
    real(8) :: dxij,dxip1j
    real(8) :: f11_lb(nnx+1),f11_ub(nnx+1)

    ! fix #6/#14: real*8 -> real(8); sih removed (unused, fix #7)
    real(8) :: si(nnx+1,nnyt)
    real(8) :: a1(nnx+1,nnyt),a2(nnx+1,nnyt),a3(nnx+1,nnyt)
    real(8) :: a4(nnx+1,nnyt),a5(nnx+1,nnyt)
    real(8) :: dphi0(nnx+1,nnyt)

    ! fix #8: ylonp, ylatp, zigm11, zigm22, zigm2, rim1, rim2 removed (all unused)
    real :: phi_weimer(nnx,nny)

    ! fix #2: DATA idpinter/ipcrit replaced with save initialisers
    !         ipcrit removed entirely (fix #3: was set but never read)
    integer, save :: idpinter = 1

    ! Local scalars
    integer :: nzh, j, jj, i, k, k0, im1, ip1, jm1, jp1, i00, j00
    real    :: dyij, dyijp1, delx, dely
    real    :: delx_inv, dely_inv, delxy_inv
    real    :: pcphimx, pcphipx, pcpmy, pcppy
    real    :: hcmjp, hcmjm, dhcy, hcmip, hcmim, dhcx, dpp
    real    :: fphivmx, fphivpx, dfphivx
    real    :: fphigmx, fphigpx, dfphigx
    real    :: fpvmy, fpvpy, dfpvy, fpgmy, fpgpy, dfpgy
    real    :: a11, a12, a21, a22, a41, a42, a51, a52
    real    :: sip1jp1, sip1jm1, sim1jp1, sim1jm1
    real    :: phi90, dbangi90, dbangif, dbangf90
    real    :: angrot, hr24, angut, angutr, phicorot
    real    :: err0, tol_phi

    if ( idpinter == 1 ) then

        nzh  = nz / 2
              
        do j = 1,nny
            jj         = j + nyextra
            dpreal(jj) = ppt(nzh,j+1,1) - ppt(nzh,j,1)
            preal(jj)  = ppt(nzh,j,1)
        enddo

        do j = nyextra,1,-1
            dpreal(j) = dpreal(j+1) * 1.4
            preal(j)  = preal(j+1) - dpreal(j)
        enddo

        do i = 1,nnx+1
            dphireal(i) = (blonp0t(i+1)-blonp0t(i))*pie/180.
        enddo

        ! define blang
        ! fix #4: original had inner j loop overwriting blang(i) on every j
        ! iteration (only j=nny survived); jj was computed but never used.
        ! Removed vestigial j loop; blang(i) now set directly at j=nny.
        do i = 1,nnx-1
            blang(i) = -blonpt(nz/2+1, nny, i) * pie / 180.
        enddo

        blang(nnx) = blang(1)

        if (restart) then
            open(U_DPHI_RST, file='dphi.rst', form='unformatted')
            read(U_DPHI_RST) dphi0
            close(U_DPHI_RST)
        endif

        idpinter = 0

    endif

!       set up conductances and driver for potential equation
!       zero-gradient in phi (x); zero-gradient in p (y)
!       note: transpose variables

    do j = 1,nny
        jj   = j + nyextra
        do i = 2,nnx-1
            hipcp_pot(i,jj)    = 0.25 * ( hipcpt(j,i-1)   + &
            hipcpt(j,i)     + &
            hipcpt(j+1,i-1) + &
            hipcpt(j+1,i)     )
            hihcm_pot(i,jj)    = 0.25 * ( hihcmt(j,i-1)   + &
            hihcmt(j,i)     + &
            hihcmt(j+1,i-1) + &
            hihcmt(j+1,i)     )
            hipcphi_pot(i,jj)  = 0.25 * ( hipcphit(j,i-1)   + &
            hipcphit(j,i)     + &
            hipcphit(j+1,i-1) + &
            hipcphit(j+1,i)     )
            hidphig_pot(i,jj)  = 0.25 * ( hidphigt(j,i-1)   + &
            hidphigt(j,i)     + &
            hidphigt(j+1,i-1) + &
            hidphigt(j+1,i)     )
            hidpg_pot(i,jj)    = 0.25 * ( hidpgt(j,i-1)   + &
            hidpgt(j,i)     + &
            hidpgt(j+1,i-1) + &
            hidpgt(j+1,i)     )
            hidphiv_pot(i,jj)  = 0.25 * ( hidphivt(j,i-1)   + &
            hidphivt(j,i)     + &
            hidphivt(j+1,i-1) + &
            hidphivt(j+1,i)     )
            hidpv_pot(i,jj)    = 0.25 * ( hidpvt(j,i-1)   + &
            hidpvt(j,i)     + &
            hidpvt(j+1,i-1) + &
            hidpvt(j+1,i)     )
            hipc_pot(i,jj)     = 0.25 * ( hipct(j,i-1)   + &
            hipct(j,i)     + &
            hipct(j+1,i-1) + &
            hipct(j+1,i)     )
            hihc_pot(i,jj)     = 0.25 * ( hihct(j,i-1)   + &
            hihct(j,i)     + &
            hihct(j+1,i-1) + &
            hihct(j+1,i)     )
            hidv_pot(i,jj)     = 0.25 * ( hidvt(j,i-1)   + &
            hidvt(j,i)     + &
            hidvt(j+1,i-1) + &
            hidvt(j+1,i)     )
        enddo
    enddo

    do j = nyextra+1,nnyt
        hipcp_pot(1,j)   = 0.5 * (hipcp_pot(2,j)    + &
        hipcp_pot(nnx-1,j) )
        hihcm_pot(1,j)   = 0.5 * (hihcm_pot(2,j)    + &
        hihcm_pot(nnx-1,j) )
        hipcphi_pot(1,j) = 0.5 * (hipcphi_pot(2,j)  + &
        hipcphi_pot(nnx-1,j) )
        hidphig_pot(1,j) = 0.5 * ( hidphig_pot(2,j) + &
        hidphig_pot(nnx-1,j) )
        hidpg_pot(1,j)   = 0.5 * ( hidpg_pot(2,j) + &
        hidpg_pot(nnx-1,j) )
        hidphiv_pot(1,j) = 0.5 * ( hidphiv_pot(2,j) + &
        hidphiv_pot(nnx-1,j) )
        hidpv_pot(1,j)   = 0.5 * ( hidpv_pot(2,j) + &
        hidpv_pot(nnx-1,j) )
        hipc_pot(1,j)    = 0.5 * (hipc_pot(2,j)    + &
        hipc_pot(nnx-1,j) )
        hihc_pot(1,j)    = 0.5 * (hihc_pot(2,j)    + &
        hihc_pot(nnx-1,j) )
        hidv_pot(1,j)    = 0.5 * (hidv_pot(2,j)    + &
        hidv_pot(nnx-1,j) )
    enddo

    do j = nyextra+1,nnyt
        hipcp_pot(nnx,j)   = hipcp_pot(1,j)
        hihcm_pot(nnx,j)   = hihcm_pot(1,j)
        hipcphi_pot(nnx,j) = hipcphi_pot(1,j)
        hidphig_pot(nnx,j) = hidphig_pot(1,j)
        hidpg_pot(nnx,j)   = hidpg_pot(1,j)
        hidphiv_pot(nnx,j) = hidphiv_pot(1,j)
        hidpv_pot(nnx,j)   = hidpv_pot(1,j)
        hipc_pot(nnx,j)    = hipc_pot(1,j)
        hihc_pot(nnx,j)    = hihc_pot(1,j)
        hidv_pot(nnx,j)    = hidv_pot(1,j)
    enddo

    do j = nyextra,1,-1
        do i = 1,nnx
            hipcp_pot(i,j)   = hipcp_pot(i,j+1)   * .02
            hihcm_pot(i,j)   = hihcm_pot(i,j+1)   * .02
            hipcphi_pot(i,j) = hipcphi_pot(i,j+1) * .02
            hidphig_pot(i,j) = hidphig_pot(i,j+1) * .01
            hidpg_pot(i,j)   = hidpg_pot(i,j+1)   * .01
            hidphiv_pot(i,j) = hidphiv_pot(i,j+1) * .01
            hidpv_pot(i,j)   = hidpv_pot(i,j+1)   * .01
            hipc_pot(i,j)    = hipc_pot(i,j+1)    * .02
            hihc_pot(i,j)    = hihc_pot(i,j+1)    * .02
            hidv_pot(i,j)    = hidv_pot(i,j+1)    * .01
        enddo
    enddo

    k   = 0
    k0  = 10

    do while ( k <= k0 )

    !     div j = 0 differencing (for Pedersen and Hall via iteration)

        do j = 2,nnyt-1
            do i = 1,nnx

                im1  = i - 1
                ip1  = i + 1
                jm1  = j - 1
                jp1  = j + 1

                dxij   = dphireal(i)
                dxip1j = dphireal(ip1)

                dyij    = dpreal(j)
                dyijp1  = dpreal(jp1)

                if ( i == 1   ) im1 = nnx - 1
                if ( i == nnx ) ip1 = 2

                delx    = dxij + dxip1j
                dely    = dyij + dyijp1
                 
                delx_inv  = 1. / delx
                dely_inv  = 1. / dely
                delxy_inv = delx_inv * dely_inv

          pcphimx  = (0.5 * ( hipcphi_pot(im1,j) &
                     + hipcphi_pot(i,j)   )) / preal(j)
          pcphipx  = (0.5 * ( hipcphi_pot(i,j)   &
                     + hipcphi_pot(ip1,j) )) / preal(j)

          pcpmy  = (0.5 * ( hipcp_pot(i,jm1) + hipcp_pot(i,j)   )) &
                 * (0.5 * ( preal(jm1) + preal(j) ) )  
          pcppy  = (0.5 * ( hipcp_pot(i,j)   + hipcp_pot(i,jp1) )) &
                 * (0.5 * ( preal(j) + preal(jp1) ) )

                if (hall) then

                    hcmjp = 0.5 * ( hihcm_pot(i,jp1) + hihcm_pot(i,j)   )
                    hcmjm = 0.5 * ( hihcm_pot(i,j)   + hihcm_pot(i,jm1) )
                    dhcy  = hcmjp - hcmjm

                    hcmip = 0.5 * ( hihcm_pot(i,j)   + hihcm_pot(ip1,j) )
                    hcmim = 0.5 * ( hihcm_pot(im1,j) + hihcm_pot(i,j)   )
                    dhcx  = hcmip - hcmim

                else

                    dhcy   = 0.
                    dhcx   = 0.
                    hcmjp  = 0.
                    hcmjm  = 0.
                    hcmip  = 0.
                    hcmim  = 0.

                endif

                fphivmx   = 0.5 * ( hidphiv_pot(im1,j) + hidphiv_pot(i,j)   )
                fphivpx   = 0.5 * ( hidphiv_pot(i,j)   + hidphiv_pot(ip1,j) )
                dfphivx   = fphivpx - fphivmx

                fphigmx   = 0.5 * ( hidphig_pot(im1,j) + hidphig_pot(i,j)   )
                fphigpx   = 0.5 * ( hidphig_pot(i,j)   + hidphig_pot(ip1,j) )
                dfphigx   = fphigpx - fphigmx


                fpvmy   = 0.5 * ( hidpv_pot(i,jm1) + hidpv_pot(i,j)   )
                fpvpy   = 0.5 * ( hidpv_pot(i,j)   + hidpv_pot(i,jp1) )
                dfpvy   = fpvpy - fpvmy

                fpgmy   = 0.5 * ( hidpg_pot(i,jm1) + hidpg_pot(i,j)   )
                fpgpy   = 0.5 * ( hidpg_pot(i,j)   + hidpg_pot(i,jp1) )
                dfpgy   = fpgpy - fpgmy

                a11    = 2. * delx_inv * pcphimx / dxij
                a12    = delxy_inv * dhcy

                a1(i,j) = a11 + a12

                a21     = 2. * dely_inv * pcpmy / dyij
                a22     = delxy_inv * dhcx

                a2(i,j) = a21 - a22

                a41     = 2. * dely_inv * pcppy / dyijp1
                a42     = delxy_inv * dhcx

                a4(i,j) = a41 + a42

                a51     = 2. * delx_inv * pcphipx / dxip1j
                a52     = delxy_inv * dhcy

                a5(i,j) = a51 - a52

                a3(i,j) = -a51 - a11 - a41 - a21

                sip1jp1 = -delxy_inv * (hcmip - hcmjp) * dphi0(ip1,jp1)

                sip1jm1 =  delxy_inv * (hcmip - hcmjm) * dphi0(ip1,jm1)

                sim1jp1 =  delxy_inv * (hcmim - hcmjp) * dphi0(im1,jp1)

                sim1jm1 = -delxy_inv * (hcmim - hcmjm) * dphi0(im1,jm1)

                si(i,j) = &
                  2. * dfphigx * delx_inv &
                + 2. * dfphivx * delx_inv &
                - 2. * dfpgy   * dely_inv &
                + 2. * dfpvy   * dely_inv &
                + sip1jp1 + sip1jm1 + sim1jp1 + sim1jm1

            enddo
        enddo

    !     zero gradient in y

        do i = 1,nnx
            a1(i,1)   = a1(i,2)
            a2(i,1)   = a2(i,2)
            a3(i,1)   = a3(i,2)
            a4(i,1)   = a4(i,2)
            a5(i,1)   = a5(i,2)
            si(i,1)   = si(i,2)

            dpp          = dpreal(nnyt) / dpreal(nnyt-1)

            a1(i,nnyt)   = dpp*(a1(i,nnyt-1)-a1(i,nnyt-2))+a1(i,nnyt-1)
            a2(i,nnyt)   = dpp*(a2(i,nnyt-1)-a2(i,nnyt-2))+a2(i,nnyt-1)
            a3(i,nnyt)   = dpp*(a3(i,nnyt-1)-a3(i,nnyt-2))+a3(i,nnyt-1)
            a4(i,nnyt)   = dpp*(a4(i,nnyt-1)-a4(i,nnyt-2))+a4(i,nnyt-1)
            a5(i,nnyt)   = dpp*(a5(i,nnyt-1)-a5(i,nnyt-2))+a5(i,nnyt-1)
            si(i,nnyt)   = dpp*(si(i,nnyt-1)-si(i,nnyt-2))+si(i,nnyt-1)

        enddo

        do j = 1,nnyt
            a1(nnx+1,j) = a1(2,j)
            a2(nnx+1,j) = a2(2,j)
            a3(nnx+1,j) = a3(2,j)
            a4(nnx+1,j) = a4(2,j)
            a5(nnx+1,j) = a5(2,j)
            si(nnx+1,j) = si(2,j)
        enddo

        do i = 2,nnx
            f11_lb(i)  = 0.d0
            f11_ub(i)  = 0.d0   ! fix #5: was commented out, leaving f11_ub uninitialised
        enddo

        f11_lb(1)     = f11_lb(nnx)
        f11_lb(nnx+1) = f11_lb(2)

        phi90 = 0.

        do i = 2,nnx
            dbangi90     = (blatpt(nz-1,nf-2,1)-90.)
            dbangif      = (blatpt(nz-1,nf-2,1) - blatpt(nz-1,nf-1,1))
            dbangf90     = (blatpt(nz-1,nf-1,1)-90.)

            dphi(i,nnyt) = ( dphi(i,nnyt-1) * dbangf90 + &
            phi90          * dbangif  ) / dbangi90

        enddo


        f11_ub(1)     = f11_ub(nnx)
        f11_ub(nnx+1) = f11_ub(2)

        if ( lmadala ) then
            call madala_sevp(a1,a2,a5,a4,a3,si,dphi,f11_lb,f11_ub)
        else
            do j = 1,nnyt
                do i = 1,nnx+1
                    dphi(i,j) = 0.
                enddo
            enddo
        endif

        i00 = nnx/2
        j00 = nnyt/2

        err0 = abs(dphi0(i00,j00)-dphi(i00,j00))/ &
        abs(dphi0(i00,j00)+1.e-5)
        k = k + 1
        print *,k,err0
         
        tol_phi  = 2.e-2

        if ( err0 <= tol_phi ) k0 = -10

        do j = 1,nnyt
            do i = 1,nnx+1
                dphi0(i,j) = dphi(i,j)
            enddo
        enddo

    enddo

!     rotate phi so that potential is
!     aligned with local time midnight/noon
!       angut  in degrees
!       angutr in radians

    angrot = 360. - plon * 180. / pie
    hr24   = mod (hrut,24.)
    angut  = hr24 * 2. * pie / 24. * rtod - angrot
    if ( angut > 360. ) angut = angut - 360.
    if ( angut < 0.   ) angut = angut + 360.
    angutr = angut * pie / 180.

    if ( lweimer ) then
        call weimer(phi_weimer,angut,hrut)
    else
        do j = nyextra+1,nnyt
            do i = 1,nnx
                jj = j - nyextra
                phi_weimer(i,jj) = 0.
            enddo
        enddo
    endif

    do j = nyextra+1,nnyt
        do i = 1,nnx
            jj        = j - nyextra
            phi(i,jj) = dphi(i,j) + phi_weimer(i,jj)
        enddo
    enddo

    print *,maxval(dphi),phi_weimer(nnx/2,nny)

    end subroutine potpphi


!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************

    subroutine weimer(phi_weimer,angut,hrut)

    use parameter_mod
    use namelist_mod
    use misc_mod
    use grid_mod

    implicit none

    ! Arguments
    real, intent(out) :: phi_weimer(nnx,nny)
    real, intent(in)  :: angut, hrut

    ! Named unit number (fix #10: was magic literal 810)
    integer, parameter :: U_WEIMER = 810

    ! Local arrays
    real :: phi_weimer_real(nfp1,nlt+1)
    real :: phi_weimer_interp(nfp1,nlt+1)

    ! fix #2: DATA iread_weimer replaced with integer,save
    ! fix #11: hrutw1 and hrutw2 must also be save -- they persist across calls
    integer, save :: iread_weimer = 1
    integer, save :: nweimer      = 0
    real,    save :: hrutw1       = 0.0
    real,    save :: hrutw2       = 0.0

    ! Local scalars
    integer :: i, j, k, nlon
    real    :: dlon, thlon, dnlon

    ! Read Weimer potential on first call
    if ( iread_weimer == 1 ) then
        open(U_WEIMER, file='phi_weimer.inp', form='unformatted')
        read(U_WEIMER) hrutw1
        read(U_WEIMER) phi_weimer_real
        read(U_WEIMER) hrutw2
        print *,'hrutw2 = ',hrut,hrutw2
        nweimer       = 1
        print *,'nweimer',nweimer
        iread_weimer  = 0
    endif

    ! Read next Weimer record when model time has advanced past current window
    if ( hrut >= hrutw2 ) then
        read(U_WEIMER) phi_weimer_real
        read(U_WEIMER) hrutw2
        print *,'hrutw2 = ',hrut,hrutw2
        nweimer = nweimer + 1
        ! fix #11: close unit when end-of-file would be next to avoid
        ! descriptor leak; reopen logic would be needed if model continues.
        ! For now close is deferred to normal model termination (file stays open
        ! for sequential reading), consistent with original design.
    endif

    ! Initialise output
    do j = 1,nny
        do i = 1,nnx
            phi_weimer(i,j) = 0.0
        enddo
    enddo

    ! fix #15: float(nlt) -> real(nlt)
    dlon = 360.0 / real(nlt)

    ! Requires SAMI3 and Weimer potential to have the same latitude
    do k = 1,nlt+1
        do j = nfp1,1,-1
            thlon = mod(blonp0t(k+1) + angut, 360.0)
            if ( thlon < 0.0 ) thlon = thlon + 360.0
            nlon  = int(thlon/dlon) + 1
            dnlon = thlon/dlon - nlon + 1
            phi_weimer_interp(j,k) = &
                  dnlon        * phi_weimer_real(j,nlon+1) &
                + (1.0 - dnlon) * phi_weimer_real(j,nlon)
        enddo
    enddo

    ! here nnx = nlt + 1, nny = nf - 1
    do j = 1,nny
        do i = 1,nnx
            phi_weimer(i,j) = phi_weimer_interp(j,i)
        enddo
    enddo

    ! Simple linear interpolation fix for last interior latitude row
    j = nny - 1
    do i = 1,nnx
        phi_weimer(i,j) = 0.5 * ( phi_weimer(i,j-1) + phi_weimer(i,j+1) )
    enddo

    end subroutine weimer
