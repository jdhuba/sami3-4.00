!****************************************
!******************************************

!           GRID3_MPI-4.00

!******************************************
!******************************************

    subroutine grid3_mpi

    use parameter_mod
    use message_passing_mod
    use grid_mod
    use exb_mod

    implicit none  ! fix #1

    include 'mpif.h'
    integer :: status(MPI_STATUS_SIZE)
    integer :: ierr  ! MPI error code (not in message_passing_mod)

!     Begin MPI stuff

    call da_grid 

!     calculate geometric parameters

!     vol is cell volume
     
    call volume 

!     areas is cell face area in i-direction (s)
!     areap is cell face area in j-direction (p)
!     areah is cell face area in k-direction (phi)

    call area 

!     xdels is line distance in i-direction (s)
!     xdelp is line distance in j-direction (p)
!     xdelh is line distance in k-direction (phi)

    call line 

!     normal: calculates normal to cell face in s-direction

    call  normal 

!     normals are calculated on s-grid
!     vpnormal: calculates e x b direction in p-direction
!     vhnormal: calculates e x b direction in h-direction (phi)

    call vpsnormal
    call vhsnormal

!     normals are calculated on p-grid
!     vpnormal: calculates e x b direction in p-direction
!     vhnormal: calculates e x b direction in h-direction (phi)

    call vppnormal
    call vhpnormal

!     unit direction in geographic latitude (theta)
!                       geographic longitude (phi)

    call gstheta
    call gsphi
    call gsr

    call  bdirs
    call  bdirh
    call  bdirp

    call facesp
    call facesh
    call facess

!    print *,'sending vpsnx '

    call mpi_send(vpsnx, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(vpsny, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(vpsnz, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(vhsnx, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(vhsny, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(vhsnz, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(bdirsx, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(bdirsy, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(bdirsz, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(gsthetax, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(gsthetay, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(gsthetaz, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(gsphix, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(gsphiy, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(gsphiz, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(gsrx, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(gsry, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(gsrz, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(vol, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)

!        print *,'send xyznorms'

    call mpi_send(xnorms, nzp1*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(ynorms, nzp1*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(znorms, nzp1*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)

!        print *,'send xyznormp'

    call mpi_send(xnormp, nz*nfp1*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(ynormp, nz*nfp1*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(znormp, nz*nfp1*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)

!        print *,'send xyznormh'

    call mpi_send(xnormh, nz*nf*nlp1, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(ynormh, nz*nf*nlp1, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(znormh, nz*nf*nlp1, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)

!        print *,'send vppnx/y/z'

    call mpi_send(vppnx, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(vppny, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(vppnz, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)

!        print *,'send vhpnx/y/z'

    call mpi_send(vhpnx, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(vhpny, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(vhpnz, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)

!        print *,'send bdirpx/y/z'

    call mpi_send(bdirpx, nz*nfp1*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(bdirpy, nz*nfp1*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(bdirpz, nz*nfp1*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)

!        print *,'send bdirhx/y/z'

    call mpi_send(bdirhx, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(bdirhy, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(bdirhz, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)

!        print *,'send ehpx/y/z'

    call mpi_send(ehpx, nz*nfp1*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(ehpy, nz*nfp1*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(ehpz, nz*nfp1*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)

!        print *,'send ephx/y/z'

    call mpi_send(ephx, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(ephy, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)
    call mpi_send(ephz, nz*nf*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)

!        print *,'send delhp'

    call mpi_send(delhp, nz*nfp1*nl, MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)

!        print *,'send delph'

    call mpi_send(delph, nz*nf*(nl+1), MPI_REAL, 0, 0, &
                  MPI_COMM_WORLD, ierr)  ! fix #5: was nz*nf*nl+1 (wrong precedence); delph(nz,nf,nlp1)


    end subroutine grid3_mpi


!     ************************
!     ************************

!           da_grid

!     ************************
!     ************************

    subroutine da_grid

!    use apex,only: apex_q2g, apex_sub, apex_mall
    use parameter_mod
    use parameter_apex_mod
    use namelist_mod
    use message_passing_mod
    use misc_mod
    use grid_mod

    implicit none  ! fix #1

    include 'mpif.h'
    integer :: status(MPI_STATUS_SIZE)

!     s grid parameters: gridding for parallel dynamics

    real :: s(nz,nf,nl)

!     p grid parameters: gridding for perpendicular dynamics

    real :: qptmp(nzp1),pval(nfp1),qpm(nfp1)
    real :: qpextra(nze+nseg+1)
    real :: fs(nz0+1),fn(nz0+1),ft(nz0+1),qpnew(nz0+1),f00(nz0+1)
    real :: qp0(nz0+1),fb0

! DEfine APEXrelated parameters (GL)
      real :: clm2
      real :: b(3),bhat(3),d1(3),d2(3),d3(3),e1(3),e2(3),e3(3),f1(2),f2(2)

    ! fix #1: previously implicit locals in da_grid
    ! ierr: MPI error code (not in message_passing_mod)
    ! blon0: local magnetic longitude for current field line
    ! qp_solve: external real function (solves q coordinate)
    ! h0: NOT declared here -- comes from parameter_apex_mod (h0=57.)
    ! mlat,mlon: real scalar outputs from apxg2all (keep original names)
    integer :: i, j, k, kk, ii, iextra, ie, iee, nz0h, iminn, imins, ier
    integer :: ierr
    real    :: blon0
    real    :: qp_solve               ! external function
    real    :: rgmin, ang0, angq, ang, bb, aa, del_blat, pvalue
    real    :: blata, blatb, blatc, rba, rbb, rbc, rga, rgb, rgc
    real    :: glat, glon, delrg, delrc, qminn, qmins, delqp, delqe
    real    :: farg, darg, f10, fa, qminm, delq
    real    :: r, br_norm, theta, brsx, grstmp, glatstmp, glonstmp
    real    :: fac1, fac2, x2, y2, zz2
    real    :: glatout, glonout, mlat, mlon, f, d
    real    :: BNRT, BEAST, BDOWN, BABS

!     defines blon from longitude_grif.f

    call blonp0a

!     define grid via p
!     grid is uniform in each region

!     will not work for IGRF-like field
!     or offset dipole

    rgmin    = altmin + re

!   new grid (arbitrary power laydown)

!   gamp      = 3. (added to namelist)

    bb        = (blat_max - blat_min)/ &
                (blat_max**gamp - blat_min**gamp)
    aa        =  blat_max - bb * blat_max**gamp

    del_blat  = (blat_max - blat_min)/float(nf)

    do j = 1,nfp1
        ang0    = blat_min + (j-1)*del_blat
        angq    = aa + bb * ang0**gamp
        ang     = angq * pie / 180.
        pval(j) = 1./cos(ang)/cos(ang)
!        print *,'j,pval',j,angq,pval(j)
    enddo

!   p grid

!   loop over longitude (phi direction)

    do k = 1,nlp1

    !       specify blon0 on worker (magnetic longitude)

        kk    = (taskid - 1) * (nl - 2) + k
        blon0 = blonp0t(kk)
                
    !        print *,k,taskid,kk,blon0

        do j = 1,nfp1
            pvalue   = pval(j)

        !         north minimum: find q value at northern most point

            blata   = 0.
            blatb   = 89.9
            blatc   = 0.5 * ( blata + blatb )
            rba     = pvalue * re * cos ( blata * po180 ) ** 2
            rbb     = pvalue * re * cos ( blatb * po180 ) ** 2
            rbc     = pvalue * re * cos ( blatc * po180 ) ** 2
            call  btog ( rba,blon0,blata,rga,glat,glon )
            call  btog ( rbb,blon0,blatb,rgb,glat,glon )
            call  btog ( rbc,blon0,blatc,rgc,glat,glon )
            delrg   = .01
            delrc   = abs ( rgc - rgmin )
            qminn   = ( re / rbc ) ** 2 * sin ( blatc * po180 )

            iminn   = 0

            do while ( delrc > delrg .AND. iminn < 20)
                iminn   = iminn + 1
                if ( rgc < rgmin ) blatb = blatc
                if ( rgc > rgmin ) blata = blatc
                blatc   = 0.5 * ( blata + blatb )
                rbc     = pvalue * re * cos ( blatc * po180 ) ** 2
                call  btog ( rbc,blon0,blatc,rgc,glat,glon )
                delrc   = abs ( rgc - rgmin )
                qminn   = ( re / rbc ) ** 2 * sin ( blatc * po180 )
            enddo
            
!            print *,'qminn ',taskid,j,k,pvalue,qminn

        !         south minimum: find q value at southern most point

            blata   = 0.
            blatb   = -89.9
            blatc   = 0.5 * ( blata + blatb )
            rba     = pvalue * re * cos ( blata * po180 ) ** 2
            rbb     = pvalue * re * cos ( blatb * po180 ) ** 2
            rbc     = pvalue * re * cos ( blatc * po180 ) ** 2
            call  btog ( rba,blon0,blata,rga,glat,glon )
            call  btog ( rbb,blon0,blatb,rgb,glat,glon )
            call  btog ( rbc,blon0,blatc,rgc,glat,glon )
            delrg   = .01
            delrc   = abs ( rgc - rgmin )
            qmins   = ( re / rbc ) ** 2 * sin ( blatc * po180 )

            imins   = 0

            do while ( delrc > delrg .AND. imins < 20 )
                imins   = imins + 1
                if ( rgc < rgmin ) blatb = blatc
                if ( rgc > rgmin ) blata = blatc
                blatc   = 0.5 * ( blata + blatb )
                rbc     = pvalue * re * cos ( blatc * po180 ) ** 2
                call  btog ( rbc,blon0,blatc,rgc,glat,glon )
                delrc   = abs ( rgc - rgmin )
                qmins   = ( re / rbc ) ** 2 * sin ( blatc * po180 )
            enddo

!            print *,'qmins ',taskid,j,k,pvalue,qmins

        !         first lay out qp uniformly on reduced grid (e.g., nz0)

            nz0h = nz0/2

            qminm = 0.
            delqp = ( qminn - qmins ) / float(nz0)
            do i = 1,nz0+1
                qp0(i) = qmins + float(i-1) * delqp
            enddo
              
!            print *,'delqp ',delqp

        !         new exponential grid laydown

            delqp = rmin * delqp ! exponential decay

            do i=1,nz0+1
                farg  = float(i-nz0h)
                darg  = float(nz0) / 10.
                f00(i) = gams + (.99-gams) * exp(-(farg/darg)**2)
            enddo

            f10  = 1.0

            do i = 1,nz0h
                fb0    = (f10-f00(i)) / (exp((qminm-qmins)/delqp) - 1.)
                fa    = f00(i) - fb0
                ft(i) = fa + fb0 * exp((qp0(i)-qmins)/delqp)
                fs(i) = ft(i)
                if ( fs(i) > f10 ) then
                    fb0   = (f10-f00(i)) / (exp(-(qminm-qminn)/delqp) - 1.)
                    fa   = f00(i) - fb0
                    ft(i)= fa + fb0* exp(-(qp0(i)-qminn)/delqp)
                endif
            enddo

            do i = nz0h+1,nz0
                fb0   = (f10-f00(i)) / (exp(-(qminm-qminn)/delqp) - 1.)
                fa   = f00(i) - fb0
                ft(i)= fa + fb0* exp(-(qp0(i)-qminn)/delqp)
                fn(i) = ft(i)
                if (fn(i) > f10 ) then
                    fb0   = (f10-f00(i)) / (exp((qminm-qmins)/delqp) - 1.)
                    fa   = f00(i) - fb0
                    ft(i)= fa + fb0* exp((qp0(i)-qmins)/delqp)
                endif
            enddo

            do i = 2,nz0h
                delq     = qp0(i) - qmins
                qpnew(i) = qmins + ft(i) * delq
                if (fs(i) > f10 ) then
                    delq     = qp0(i) - qminn
                    qpnew(i) = qminn + ft(i) * delq
                endif
            enddo

            do i = nz0,nz0h+1,-1
                delq     = qp0(i) - qminn
                qpnew(i) = qminn + ft(i) * delq
                if (fn(i) > f10 ) then
                    delq     = qp0(i) - qmins
                    qpnew(i) = qmins + ft(i) * delq
                endif
            enddo

            qpnew(1)     = qmins
            qpnew(nz0+1) = qminn

            iextra = 0

            do ie = 1,nseg
                do iee = 1,nextra+1
                    iextra = iextra + 1
                    delqe   = qpnew(nz0h+(1-nseg/2)+ie) - &
                    qpnew(nz0h+(1-nseg/2)+(ie-1))
                    qpextra(iextra) = qpnew(nz0h+(1-nseg/2)+(ie-1)) &
                    + delqe/float(nextra+1) * (iee-1)
                enddo
            enddo

            qpextra(nze+nseg+1) = qpnew(nz0h+(1+nseg/2))

            do i = 1,nz0h+(1-nseg/2)-1
                qptmp(i) = qpnew(i)
            enddo

            do i = nz0h+(1-nseg/2),nz0h+(1-nseg/2)+(nze+nseg+1)-1
                qptmp(i) = qpextra(i+1-(nz0h+(1-nseg/2)))
            enddo

            do i = nz0h+(1-nseg/2)+(nze+nseg+1),nz+1
                ii    = i + 1 - (nz0h+(1-nseg/2)+(nze+nseg+1)) &
                + nz0h + (1+nseg/2)
                qptmp(i) = qpnew(ii)
            enddo

            do i = 1,nz+1
                qp(i,j,k) = qptmp(i)
            enddo

        enddo

    !     define p grid

        do j = 1,nfp1
            do i = 1,nzp1

                r = re*qp_solve(qp(i,j,k),pval(j))
                pp(i,j,k)     = pval(j)
                brp(i,j,k)    = r
                br_norm       = r / re
                blatp(i,j,k)  = asin ( qp(i,j,k) * br_norm ** 2 ) * rtod
                blonp(i,j,k)  = blon0
                baltp(i,j,k)  = brp(i,j,k) - re   ! altitude above earth

            !     p grid in cartesian coordinates

                xp(i,j,k)    = brp(i,j,k) * cos( blatp(i,j,k)*po180 ) * &
                sin( blonp(i,j,k)*po180 )
                yp(i,j,k)    = brp(i,j,k) * sin( blatp(i,j,k)*po180 )
                zp(i,j,k)    = brp(i,j,k) * cos( blatp(i,j,k)*po180 ) * &
                cos( blonp(i,j,k) *po180)

            enddo ! end i (along the field)
        enddo ! end j (altitude)
    enddo ! end k (longitude)

! Now we know altp blatp blonp and we need to send them to master

    print *,'sending baltp ',taskid

    call mpi_send(baltp, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(blatp, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(blonp, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)

! Now we know xp yp zp and we need to send them to master

    call mpi_send(xp, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(yp, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(zp, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
     print *,'sending pp',taskid
    call mpi_send(pp, nzp1*nfp1*nlp1, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)

!     now do s grid

!     define qs and ps (average of qp and pp)

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz

                xs(i,j,k)    = 0.125 * &
                ( xp(i,j,k)     + xp(i,j+1,k)   + &
                xp(i+1,j,k)   + xp(i+1,j+1,k) + &
                xp(i,j,k+1)   + xp(i,j+1,k+1) + &
                xp(i+1,j,k+1) + xp(i+1,j+1,k+1) )

                ys(i,j,k)    = 0.125 * &
                ( yp(i,j,k)     + yp(i,j+1,k)   + &
                yp(i+1,j,k)   + yp(i+1,j+1,k) + &
                yp(i,j,k+1)   + yp(i,j+1,k+1) + &
                yp(i+1,j,k+1) + yp(i+1,j+1,k+1) )

                zs(i,j,k)    = 0.125 * &
                ( zp(i,j,k)     + zp(i,j+1,k)   + &
                zp(i+1,j,k)   + zp(i+1,j+1,k) + &
                zp(i,j,k+1)   + zp(i,j+1,k+1) + &
                zp(i+1,j,k+1) + zp(i+1,j+1,k+1) )

                brs(i,j,k)   = sqrt(xs(i,j,k)**2+ys(i,j,k)**2+zs(i,j,k)**2)
                blats(i,j,k) = asin(ys(i,j,k)/brs(i,j,k))/po180
                blons(i,j,k) = atan2(xs(i,j,k),zs(i,j,k))/po180
                if(blons(i,j,k) < 0) blons(i,j,k)=360.+blons(i,j,k)

            !           define qs and ps

                qs(i,j,k)     = ( re / brs(i,j,k) ) ** 2 * &
                sin ( blats(i,j,k) * po180 )
                ps(i,j,k)     = ( brs(i,j,k) / re ) / &
                ( cos ( blats(i,j,k) * po180 ) ** 2 )

            !            if (taskid.eq.1 .and. k.eq.1 .and. j.eq.nf/2)
            !     .   print *,i,qs(i,j,k),ps(i,j,k),pp(i,j,k)

            !            qtmp   = qs(i,j,k)
            !            ptmp   = ps(i,j,k)
            !            r = re*qp_solve(qtmp,ptmp)
            !            brs(i,j,k)    = r
            !            br_norm       = r / re
            !            blats(i,j,k)  = asin ( qs(i,j,k) * br_norm ** 2 ) * rtod
                 
!                call btog ( brs(i,j,k),blons(i,j,k),blats(i,j,k), &
!                grstmp,glatstmp,glonstmp )

! apex btog (GL)
! change re doug

            brsx = brs(i,j,k) - re
             call apxq2g(blats(i,j,k),blons(i,j,k),brsx,1.e-3,&
                glatstmp,glonstmp,ier)
            if (glonstmp < 0.) glonstmp = glonstmp + 360.
            grstmp = brs(i,j,k)

                grs(i,j,k)    = grstmp
                glats(i,j,k)  = glatstmp
                glons(i,j,k)  = glonstmp
                glons0(i,j,k) = glonstmp
                 
            !     theta is in radians
                 
                theta = acos ( qs(i,j,k) * (brs(i,j,k)/re) ** 2 )

            !     get the dimensionless magnetic field: bm = b/b0

                br_norm       = brs(i,j,k) / re


! Recalculate B-field from IGRF (GL) re Doug Drob

  call apxg2all(glats(i,j,k),glons(i,j,k),brsx,h0,1, &
                glatout,glonout,mlat,mlon,f1,f2,f,d1,d2,d3,d,e1,e2,e3)

  call FELDG(1,glats(i,j,k),glons(i,j,k),brsx,BNRT,BEAST,BDOWN,BABS)


!  Convert from nT to Gauss, then normalized by bmag
            bms(i,j,k) = babs/bmag
! Calculate d1 and d2 from Richmond 1995 reference paper
! and normalize d1 and d2
            fac1 = sqrt(d1(1)**2+d1(2)**2+d1(3)**2)
            fac2 = sqrt(d2(1)**2+d2(2)**2+d2(3)**2)
            dvec(i,j,k,1,1) = d1(1) 
            dvec(i,j,k,1,2) = d1(2) 
            dvec(i,j,k,1,3) = d1(3) 
            dvec(i,j,k,2,1) = d2(1) 
            dvec(i,j,k,2,2) = d2(2) 
            dvec(i,j,k,2,3) = d2(3) 
! new re doug
            dvec(i,j,k,3,1) = beast/babs
            dvec(i,j,k,3,2) = bnrt/babs
            dvec(i,j,k,3,3) = -bdown/babs

            call vector(brs(i,j,k),blons(i,j,k),blats(i,j,k), &
                        xrg(i,j,k),xthg(i,j,k),xphig(i,j,k)  )

            !     we define the dimensional grid points given by s = q * re

                s(i,j,k)    = qs(i,j,k)  * re
                alts(i,j,k) = grs(i,j,k) - re   ! altitude above earth

            !     grid in cartesian coordinates

                xs(i,j,k)    = brs(i,j,k) * cos( blats(i,j,k)*po180 ) * &
                sin( blons(i,j,k)*po180 )
                ys(i,j,k)    = brs(i,j,k) * sin( blats(i,j,k)*po180 )
                zs(i,j,k)    = brs(i,j,k) * cos( blats(i,j,k)*po180 ) * &
                cos( blons(i,j,k)*po180 )
            enddo
        enddo
    enddo

!         print *,'sending data',taskid

! sending dvec

    print *,'sending dvec ',taskid

    call mpi_send(dvec(:,:,:,1,1), nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(dvec(:,:,:,1,2), nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(dvec(:,:,:,1,3), nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
!    print *,'done dvec1 ',taskid

    call mpi_send(dvec(:,:,:,2,1), nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(dvec(:,:,:,2,2), nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(dvec(:,:,:,2,3), nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)

    call mpi_send(dvec(:,:,:,3,1), nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(dvec(:,:,:,3,2), nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(dvec(:,:,:,3,3), nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)



! Now we know alts glats glons and we need to send them to master

    call mpi_send(alts, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(glats, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(glons, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)

    call mpi_send(brs, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(blats, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(blons, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)

    call mpi_send(xs, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(ys, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(zs, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)

    call mpi_send(xrg, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(xthg, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)
    call mpi_send(xphig, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)

    call mpi_send(bms, nz*nf*nl, MPI_REAL, 0, 0, &
    MPI_COMM_WORLD, ierr)


!     following distances are in cm
         
    do k = 1,nl
        do j = 1,nf
            do i = 2,nz-1
                ds(i,j,k)   = ( s(i,j,k)   - s(i-1,j,k) ) * 1.e5
                d2s(i,j,k)  = ( s(i+1,j,k) - s(i-1,j,k) ) * 1.e5
                d22s(i,j,k) = .5 * d2s(i,j,k)
            enddo
            ds(1,j,k)    = ds(2,j,k)
            ds(nz,j,k)   = ds(nz-1,j,k)
            d2s(1,j,k)   = d2s(2,j,k)
            d2s(nz,j,k)  = d2s(nz-1,j,k)
            d22s(1,j,k)  = d22s(2,j,k)
            d22s(nz,j,k) = d22s(nz-1,j,k)
        enddo
    enddo

!     calculate dels: grid length along field line using cartesian geometry
!     and convert to cm from km

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz-1
                x2  = ( xs(i+1,j,k) - xs(i,j,k) ) ** 2
                y2  = ( ys(i+1,j,k) - ys(i,j,k) ) ** 2
                zz2 = ( zs(i+1,j,k) - zs(i,j,k) ) ** 2
                dels(i,j,k) = sqrt ( x2 + y2 + zz2 ) * 1.e5
            enddo
            dels(nz,j,k) = dels(nz-1,j,k)
        enddo
    enddo

    return
    end subroutine da_grid


!******************************************
!******************************************

!            gtob

!******************************************
!******************************************

    subroutine gtob(brad,blond,blatd,grad,glatd,glond)

    use parameter_mod
    use grid_mod
    use namelist_mod
    use misc_mod, only: plat, plon  ! plat/plon: magnetic pole coords

    implicit none  ! fix #1

    real, intent(in)  :: grad, glatd, glond
    real, intent(out) :: brad, blond, blatd
    real :: glat, glon, xg, yg, zg, xm, ym, zm

    glat = glatd * po180
    glon = glond * po180

!     rotate in longitude

    xg = grad * cos ( glat ) * cos ( glon - plon )
    yg = grad * cos ( glat ) * sin ( glon - plon )
    zg = grad * sin ( glat )

!     rotate in latitude

    xm = xg * sin ( plat ) - zg * cos ( plat )
    ym = yg
    zm = xg * cos ( plat ) + zg * sin ( plat )

!     magnetic lat and long converted to degrees

    brad  = sqrt ( xm ** 2 + ym ** 2 + zm ** 2 )
    blatd = asin ( zm / brad ) * rtod
    blond = ( atan2 ( ym/brad,xm/brad ) ) * rtod
    if (blond < -0.01) blond = blond+360.

    return
    end subroutine gtob

!******************************************
!******************************************

!            btog

!******************************************
!******************************************


    subroutine btog ( brad,blond,blatd,grad,glatd,glond )

!     conversion from centered dipole to geographic coordinates
!     plat,plon =  coords of north magnetic pole (in param2d.9.14e.inc)
!     the geomagnetic longitude is measured east from the
!     geomagnetic prime meridian at 291 degrees east geographic.
!     brad: radius in the offset dipole frame
!     grad: radius in the geocentric frame
!     (g joyce june 1998)

    use parameter_mod
    use grid_mod
    use namelist_mod

    implicit none  ! fix #1

    real, intent(in)  :: brad, blond, blatd
    real, intent(out) :: grad, glatd, glond
    real :: plat0, plon0, blonr, blatr, xmm, ymm, zmm, xg, yg, zg

!     NO OFFSET (just tilt)

    plat0 = 1.57
    plon0 = 6.28

!     convert magnetic lat and long to radians

    blonr = blond * po180
    blatr = blatd * po180

!     position of point in geomagnetic coords
!     first get xm ym zm in the eccentric dipole system

    xmm = brad *cos ( blatr ) * cos ( blonr )
    ymm = brad *cos ( blatr ) * sin ( blonr )
    zmm = brad *sin ( blatr )

!     r is invariant under the rotations of the tilted dipole

    grad = sqrt ( xmm ** 2 + ymm ** 2 + zmm ** 2 )

!     rotate coords in north-south direction

    xg =  xmm * sin ( plat0 ) + zmm * cos ( plat0 )
    yg =  ymm
    zg = -xmm * cos ( plat0 ) + zmm * sin ( plat0 )

!     geographic latitude and longitude converted to degrees

    glatd = asin ( zg / grad ) * rtod

    if ( glatd > 90. )  glatd =  90.
    if ( glatd < -90. ) glatd = -90.

    glond = ( plon0 + atan2 ( yg/grad,xg/grad ) ) * rtod

    if (glond >= 360.) glond = glond -  360.

!     print *,'plat0/plon0 ',plat0,plon0

    return
    end subroutine btog

!******************************************
!******************************************

!            btog_xyz

!******************************************
!******************************************


    subroutine btog_xyz ( brad,blond,blatd,x,y,z )

!     conversion from centered dipole to geographic coordinates
!     plat,plon =  coords of north magnetic pole (in param2d.9.14e.inc)
!     the geomagnetic longitude is measured east from the
!     geomagnetic prime meridian at 291 degrees east geographic.
!     brad: radius in the offset dipole frame
!     grad: radius in the geocentric frame
!     (g joyce june 1998)

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real, intent(in)  :: brad, blond, blatd
    real, intent(out) :: x, y, z
    real :: blonr, blatr

    blonr = blond * po180
    blatr = blatd * po180

!     position of point in geomagnetic coords
!     first get xm ym zm in the eccentric dipole system

    z = brad *cos ( blatr ) * cos ( blonr )
    x = brad *cos ( blatr ) * sin ( blonr )
    y = brad *sin ( blatr )

    return
    end subroutine btog_xyz


!******************************************
!******************************************

!             volume

!******************************************
!******************************************

    subroutine volume

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: voli(nz,nf,nl),volj(nz,nf,nl),volk(nz,nf,nl)
    integer :: i, j, k
    real :: xmid, ymid, zmid
    real :: ax1, ay1, az1, bx1, by1, bz1, cx1, cy1, cz1
    real :: dx1, dy1, dz1, ax2, ay2, az2
    real :: v1, v2, v3, v4

!       volume from sidei

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz

                xmid  = 0.5 * ( xp(i,j,k) + xp(i+1,j+1,k+1) )
                ymid  = 0.5 * ( yp(i,j,k) + yp(i+1,j+1,k+1) )
                zmid  = 0.5 * ( zp(i,j,k) + zp(i+1,j+1,k+1) )
                              
                ax1 = xp(i,j,k) - xmid
                ay1 = yp(i,j,k) - ymid
                az1 = zp(i,j,k) - zmid

                bx1 = xp(i,j+1,k) - xmid
                by1 = yp(i,j+1,k) - ymid
                bz1 = zp(i,j+1,k) - zmid

                cx1 = xp(i,j,k+1) - xmid
                cy1 = yp(i,j,k+1) - ymid
                cz1 = zp(i,j,k+1) - zmid

                dx1 =    by1 * cz1 - bz1 * cy1
                dy1 = -( bx1 * cz1 - bz1 * cx1 )
                dz1 =    bx1 * cy1 - by1 * cx1

                v1 = (1.0/6.0) * &
                abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))

                ax2 = xp(i,j+1,k+1) - xmid
                ay2 = yp(i,j+1,k+1) - ymid
                az2 = zp(i,j+1,k+1) - zmid

                v2 = (1.0/6.0) * &
                abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

                ax1 = xp(i+1,j,k) - xmid
                ay1 = yp(i+1,j,k) - ymid
                az1 = zp(i+1,j,k) - zmid

                bx1 = xp(i+1,j+1,k) - xmid
                by1 = yp(i+1,j+1,k) - ymid
                bz1 = zp(i+1,j+1,k) - zmid

                cx1 = xp(i+1,j,k+1) - xmid
                cy1 = yp(i+1,j,k+1) - ymid
                cz1 = zp(i+1,j,k+1) - zmid

                dx1 =    by1 * cz1 - bz1 * cy1
                dy1 = -( bx1 * cz1 - bz1 * cx1 )
                dz1 =    bx1 * cy1 - by1 * cx1

                v3 = (1.0/6.0) * &
                abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))
                       
                ax2 = xp(i+1,j+1,k+1) - xmid
                ay2 = yp(i+1,j+1,k+1) - ymid
                az2 = zp(i+1,j+1,k+1) - zmid

                v4 = (1.0/6.0) * &
                abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

                voli(i,j,k) = v1 + v2 + v3 + v4
                               
            enddo
        enddo
    enddo

!       volume from sidej

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz

                xmid = 0.5 * ( xp(i,j,k) + xp(i+1,j+1,k+1) )
                ymid = 0.5 * ( yp(i,j,k) + yp(i+1,j+1,k+1) )
                zmid = 0.5 * ( zp(i,j,k) + zp(i+1,j+1,k+1) )

                ax1 = xp(i,j,k) - xmid
                ay1 = yp(i,j,k) - ymid
                az1 = zp(i,j,k) - zmid

                bx1 = xp(i+1,j,k) - xmid
                by1 = yp(i+1,j,k) - ymid
                bz1 = zp(i+1,j,k) - zmid

                cx1 = xp(i,j,k+1) - xmid
                cy1 = yp(i,j,k+1) - ymid
                cz1 = zp(i,j,k+1) - zmid

                dx1 =    by1 * cz1 - bz1 * cy1
                dy1 = -( bx1 * cz1 - bz1 * cx1 )
                dz1 =    bx1 * cy1 - by1 * cx1

                v1 = (1.0/6.0) * &
                abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))

                ax2 = xp(i+1,j,k+1) - xmid
                ay2 = yp(i+1,j,k+1) - ymid
                az2 = zp(i+1,j,k+1) - zmid

                v2 = (1.0/6.0) * &
                abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

                ax1 = xp(i,j+1,k) - xmid
                ay1 = yp(i,j+1,k) - ymid
                az1 = zp(i,j+1,k) - zmid

                bx1 = xp(i+1,j+1,k) - xmid
                by1 = yp(i+1,j+1,k) - ymid
                bz1 = zp(i+1,j+1,k) - zmid

                cx1 = xp(i,j+1,k+1) - xmid
                cy1 = yp(i,j+1,k+1) - ymid
                cz1 = zp(i,j+1,k+1) - zmid

                dx1 =    by1 * cz1 - bz1 * cy1
                dy1 = -( bx1 * cz1 - bz1 * cx1 )
                dz1 =    bx1 * cy1 - by1 * cx1

                v3 = (1.0/6.0) * &
                abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))
                       
                ax2 = xp(i+1,j+1,k+1) - xmid
                ay2 = yp(i+1,j+1,k+1) - ymid
                az2 = zp(i+1,j+1,k+1) - zmid

                v4 = (1.0/6.0) * &
                abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

                volj(i,j,k) = v1 + v2 + v3 + v4
                 
            enddo
        enddo
    enddo

!       volume from sidek

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz

                xmid = 0.5 * ( xp(i,j,k) + xp(i+1,j+1,k+1) )
                ymid = 0.5 * ( yp(i,j,k) + yp(i+1,j+1,k+1) )
                zmid = 0.5 * ( zp(i,j,k) + zp(i+1,j+1,k+1) )

                ax1 = xp(i,j,k) - xmid
                ay1 = yp(i,j,k) - ymid
                az1 = zp(i,j,k) - zmid

                bx1 = xp(i+1,j,k) - xmid
                by1 = yp(i+1,j,k) - ymid
                bz1 = zp(i+1,j,k) - zmid

                cx1 = xp(i,j+1,k) - xmid
                cy1 = yp(i,j+1,k) - ymid
                cz1 = zp(i,j+1,k) - zmid

                dx1 =    by1 * cz1 - bz1 * cy1
                dy1 = -( bx1 * cz1 - bz1 * cx1 )
                dz1 =    bx1 * cy1 - by1 * cx1

                v1 = (1.0/6.0) * &
                abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))

                ax2 = xp(i+1,j+1,k) - xmid
                ay2 = yp(i+1,j+1,k) - ymid
                az2 = zp(i+1,j+1,k) - zmid

                v2 = (1.0/6.0) * &
                abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

                ax1 = xp(i,j,k+1) - xmid
                ay1 = yp(i,j,k+1) - ymid
                az1 = zp(i,j,k+1) - zmid

                bx1 = xp(i+1,j,k+1) - xmid
                by1 = yp(i+1,j,k+1) - ymid
                bz1 = zp(i+1,j,k+1) - zmid

                cx1 = xp(i,j+1,k+1) - xmid
                cy1 = yp(i,j+1,k+1) - ymid
                cz1 = zp(i,j+1,k+1) - zmid

                dx1 =    by1 * cz1 - bz1 * cy1
                dy1 = -( bx1 * cz1 - bz1 * cx1 )
                dz1 =    bx1 * cy1 - by1 * cx1

                v3 = (1.0/6.0) * &
                abs(( ax1 * dx1 + ay1 * dy1 + az1 * dz1 ))
                       
                ax2 = xp(i+1,j+1,k+1) - xmid
                ay2 = yp(i+1,j+1,k+1) - ymid
                az2 = zp(i+1,j+1,k+1) - zmid

                v4 = (1.0/6.0) * &
                abs(( ax2 * dx1 + ay2 * dy1 + az2 * dz1 ))

                volk(i,j,k) = v1 + v2 + v3 + v4

            enddo
        enddo
    enddo

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                vol(i,j,k) = 1.e15 * &
                ( voli(i,j,k) + volj(i,j,k) + volk(i,j,k) )
                 
            enddo
        enddo
    enddo

    return
    end subroutine volume

!******************************************
!******************************************

!            area

!******************************************
!******************************************


    subroutine area

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    integer :: i, j, k
    real :: ax1,ax2,ay1,ay2,az1,az2,bx1,bx2,by1,by2,bz1,bz2
    real :: cx1,cx2,cy1,cy2,cz1,cz2,a1,a2

!       sidei (s-direction)

    do k = 1,nl
        do j = 1,nf
            do i = 1,nzp1

                ax1 = xp(i,j+1,k) - xp(i,j,k)
                ay1 = yp(i,j+1,k) - yp(i,j,k)
                az1 = zp(i,j+1,k) - zp(i,j,k)

                bx1 = xp(i,j,k+1) - xp(i,j,k)
                by1 = yp(i,j,k+1) - yp(i,j,k)
                bz1 = zp(i,j,k+1) - zp(i,j,k)

                cx1 =    ay1 * bz1 - az1 * by1
                cy1 = -( ax1 * bz1 - az1 * bx1 )
                cz1 =    ax1 * by1 - ay1 * bx1

                a1 = 0.5 * sqrt ( cx1*cx1 + cy1*cy1 + cz1*cz1 )

                ax2 = xp(i,j+1,k) - xp(i,j+1,k+1)
                ay2 = yp(i,j+1,k) - yp(i,j+1,k+1)
                az2 = zp(i,j+1,k) - zp(i,j+1,k+1)

                bx2 = xp(i,j,k+1) - xp(i,j+1,k+1)
                by2 = yp(i,j,k+1) - yp(i,j+1,k+1)
                bz2 = zp(i,j,k+1) - zp(i,j+1,k+1)

                cx2 =    ay2 * bz2 - az2 * by2
                cy2 = -( ax2 * bz2 - az2 * bx2 )
                cz2 =    ax2 * by2 - ay2 * bx2

                a2 = 0.5 * sqrt ( cx2*cx2 + cy2*cy2 + cz2*cz2 )

                areas(i,j,k) = ( a1 + a2 ) * 1.e10

            enddo
        enddo
    enddo

!       sidej (p-direction)

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nz

                ax1 = xp(i+1,j,k) - xp(i,j,k)
                ay1 = yp(i+1,j,k) - yp(i,j,k)
                az1 = zp(i+1,j,k) - zp(i,j,k)

                bx1 = xp(i,j,k+1) - xp(i,j,k)
                by1 = yp(i,j,k+1) - yp(i,j,k)
                bz1 = zp(i,j,k+1) - zp(i,j,k)

                cx1 =    ay1 * bz1 - az1 * by1
                cy1 = -( ax1 * bz1 - az1 * bx1 )
                cz1 =    ax1 * by1 - ay1 * bx1

                a1 = 0.5 * sqrt ( cx1*cx1 + cy1*cy1 + cz1*cz1 )
                       
                ax2 = xp(i+1,j,k) - xp(i+1,j,k+1)
                ay2 = yp(i+1,j,k) - yp(i+1,j,k+1)
                az2 = zp(i+1,j,k) - zp(i+1,j,k+1)

                bx2 = xp(i,j,k+1) - xp(i+1,j,k+1)
                by2 = yp(i,j,k+1) - yp(i+1,j,k+1)
                bz2 = zp(i,j,k+1) - zp(i+1,j,k+1)

                cx2 =    ay2 * bz2 - az2 * by2
                cy2 = -( ax2 * bz2 - az2 * bx2 )
                cz2 =    ax2 * by2 - ay2 * bx2

                a2 = 0.5 * sqrt ( cx2*cx2 + cy2*cy2 + cz2*cz2 )
                       
                areap(i,j,k) = ( a1 + a2 ) * 1.e10

            enddo
        enddo
    enddo

!       sidek (phi-direction)

    do k = 1,nlp1
        do j = 1,nf
            do i = 1,nz

                ax1 = xp(i+1,j,k) - xp(i,j,k)
                ay1 = yp(i+1,j,k) - yp(i,j,k)
                az1 = zp(i+1,j,k) - zp(i,j,k)

                bx1 = xp(i,j+1,k) - xp(i,j,k)
                by1 = yp(i,j+1,k) - yp(i,j,k)
                bz1 = zp(i,j+1,k) - zp(i,j,k)

                cx1 =    ay1 * bz1 - az1 * by1
                cy1 = -( ax1 * bz1 - az1 * bx1 )
                cz1 =    ax1 * by1 - ay1 * bx1

                a1 = 0.5 * sqrt ( cx1*cx1 + cy1*cy1 + cz1*cz1 )

                ax2 = xp(i+1,j,k) - xp(i+1,j+1,k)
                ay2 = yp(i+1,j,k) - yp(i+1,j+1,k)
                az2 = zp(i+1,j,k) - zp(i+1,j+1,k)

                bx2 = xp(i,j+1,k) - xp(i+1,j+1,k)
                by2 = yp(i,j+1,k) - yp(i+1,j+1,k)
                bz2 = zp(i,j+1,k) - zp(i+1,j+1,k)

                cx2 =    ay2 * bz2 - az2 * by2
                cy2 = -( ax2 * bz2 - az2 * bx2 )
                cz2 =    ax2 * by2 - ay2 * bx2

                a2 = 0.5 * sqrt ( cx2*cx2 + cy2*cy2 + cz2*cz2 )

                areah(i,j,k) = ( a1 + a2 ) * 1.e10

            enddo
        enddo
    enddo

    return
    end subroutine area

!******************************************
!******************************************

!            line

!******************************************
!******************************************


    subroutine line

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    integer :: i, j, k
    real :: ax1, ay1, az1

!       xdels (s-direction)

    do k = 1,nlp1
        do j = 1,nfp1
            do i = 1,nz

                ax1 = xp(i+1,j,k) - xp(i,j,k)
                ay1 = yp(i+1,j,k) - yp(i,j,k)
                az1 = zp(i+1,j,k) - zp(i,j,k)

                xdels(i,j,k) = sqrt ( ax1*ax1 + ay1*ay1 + az1*az1 ) &
                * 1.e5

            enddo
        enddo
    enddo

!       xdelp (p-direction)

    do k = 1,nlp1
        do j = 1,nf
            do i = 1,nzp1

                ax1 = xp(i,j+1,k) - xp(i,j,k)
                ay1 = yp(i,j+1,k) - yp(i,j,k)
                az1 = zp(i,j+1,k) - zp(i,j,k)

                xdelp(i,j,k) = sqrt ( ax1*ax1 + ay1*ay1 + az1*az1 ) &
                * 1.e5

            enddo
        enddo
    enddo

!       xdelh (phi-direction)

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nzp1

                ax1 = xp(i,j,k+1) - xp(i,j,k)
                ay1 = yp(i,j,k+1) - yp(i,j,k)
                az1 = zp(i,j,k+1) - zp(i,j,k)

                xdelh(i,j,k) = sqrt ( ax1*ax1 + ay1*ay1 + az1*az1 ) &
                * 1.e5

            enddo
        enddo
    enddo

    return
    end subroutine line

!******************************************
!******************************************

!            normal

!******************************************
!******************************************


    subroutine normal

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    integer :: i, j, k
    real :: ax1,ay1,az1,bx1,by1,bz1,cx1,cy1,cz1,ca1

!       norms (normal to cell face in s-direction)

    do k = 1,nl
        do j = 1,nf
            do i = 1,nzp1
                            
!!$                ax1 = xp(i,j,k+1) - xp(i,j+1,k)
!!$                ay1 = yp(i,j,k+1) - yp(i,j+1,k)
!!$                az1 = zp(i,j,k+1) - zp(i,j+1,k)
!!$
!!$                bx1 = xp(i,j+1,k+1) - xp(i,j,k)
!!$                by1 = yp(i,j+1,k+1) - yp(i,j,k)
!!$                bz1 = zp(i,j+1,k+1) - zp(i,j,k)

! original 

                ax1 = xp(i,j+1,k+1) - xp(i,j,k)
                ay1 = yp(i,j+1,k+1) - yp(i,j,k)
                az1 = zp(i,j+1,k+1) - zp(i,j,k)

                bx1 = xp(i,j,k+1) - xp(i,j+1,k)
                by1 = yp(i,j,k+1) - yp(i,j+1,k)
                bz1 = zp(i,j,k+1) - zp(i,j+1,k)

                cx1 =   ay1 * bz1 - az1 * by1
                cy1 = -(ax1 * bz1 - az1 * bx1)
                cz1 =   ax1 * by1 - ay1 * bx1

                ca1 = sqrt( cx1*cx1 + cy1*cy1 + cz1*cz1 )

                xnorms(i,j,k) = cx1/ca1
                ynorms(i,j,k) = cy1/ca1
                znorms(i,j,k) = cz1/ca1

            enddo
        enddo
    enddo

!       normp (normal to cell face in p-direction)

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nz
                            
!!$                ax1 = xp(i+1,j,k+1) - xp(i,j,k)
!!$                ay1 = yp(i+1,j,k+1) - yp(i,j,k)
!!$                az1 = zp(i+1,j,k+1) - zp(i,j,k)
!!$
!!$                bx1 = xp(i,j,k+1) - xp(i+1,j,k)
!!$                by1 = yp(i,j,k+1) - yp(i+1,j,k)
!!$                bz1 = zp(i,j,k+1) - zp(i+1,j,k)

! original

                ax1 = xp(i,j,k+1) - xp(i+1,j,k)
                ay1 = yp(i,j,k+1) - yp(i+1,j,k)
                az1 = zp(i,j,k+1) - zp(i+1,j,k)

                bx1 = xp(i+1,j,k+1) - xp(i,j,k)
                by1 = yp(i+1,j,k+1) - yp(i,j,k)
                bz1 = zp(i+1,j,k+1) - zp(i,j,k)

                cx1 =   ay1 * bz1 - az1 * by1
                cy1 = -(ax1 * bz1 - az1 * bx1)
                cz1 =   ax1 * by1 - ay1 * bx1

                ca1 = sqrt( cx1*cx1 + cy1*cy1 + cz1*cz1 )

                xnormp(i,j,k) = cx1/ca1
                ynormp(i,j,k) = cy1/ca1
                znormp(i,j,k) = cz1/ca1

            enddo
        enddo
    enddo

!       normh (normal to cell face in phi-direction: longitude)

    do  k = 1,nlp1
        do j = 1,nf
            do i = 1,nz
                              
!!$                ax1 = xp(i+1,j+1,k) - xp(i,j,k)
!!$                ay1 = yp(i+1,j+1,k) - yp(i,j,k)
!!$                az1 = zp(i+1,j+1,k) - zp(i,j,k)
!!$
!!$                bx1 = xp(i,j+1,k) - xp(i+1,j,k)
!!$                by1 = yp(i,j+1,k) - yp(i+1,j,k)
!!$                bz1 = zp(i,j+1,k) - zp(i+1,j,k)

! original

                ax1 = xp(i+1,j,k) - xp(i,j+1,k)
                ay1 = yp(i+1,j,k) - yp(i,j+1,k)
                az1 = zp(i+1,j,k) - zp(i,j+1,k)

                bx1 = xp(i+1,j+1,k) - xp(i,j,k)
                by1 = yp(i+1,j+1,k) - yp(i,j,k)
                bz1 = zp(i+1,j+1,k) - zp(i,j,k)

                cx1 =   ay1 * bz1 - az1 * by1
                cy1 = -(ax1 * bz1 - az1 * bx1)
                cz1 =   ax1 * by1 - ay1 * bx1

                ca1 = sqrt( cx1*cx1 + cy1*cy1 + cz1*cz1 )

                xnormh(i,j,k) = cx1/ca1
                ynormh(i,j,k) = cy1/ca1
                znormh(i,j,k) = cz1/ca1

            enddo
        enddo
    enddo

    return
    end subroutine normal



!       ********************************************************
!       ********************************************************

    function qp_solve(q,p)

!       ********************************************************
!       ********************************************************

    real :: qp_solve,q,p
    real :: term1,term2

!       MS: Functionally unnecessary, but kind of neat.  To get r from
!       (q,p) one needs to find the root of a fourth-order polynomial.
!       The code did this with a standard root-finding method,
!       but, of course, there is a closed-form solution.  Since the
!       polynomial has no second or third order terms, the answer
!       isn't completely ugly, and the code below gives the exact
!       result.  Note the special case for q = 0, i.e. points on
!       the magnetic equator.
      
! MS: The formula is actually my third attempt.  The first had huge
! cancellation problems when q was small, the second had smaller
! problems when term0 was large. This should be well-behaved
! everywhere.  Massive algebra was involved along the way.

    term0 = 256./27.*q**2*p**4
    term1 = ( 1. + sqrt(1. + term0) )**(2./3.)
    term2 = term0**(1./3.)
    term3 = 0.5 * ( (term1**2 + term1*term2 + &
    term2**2)/term1 )**(3./2.)
    qp_solve = p * (4.*term3) / (1.+term3) / &
    ( 1.+sqrt(2.*term3-1.) )


    return
    end function qp_solve


!******************************************
!******************************************

!            vpsnormal

!******************************************
!******************************************


    subroutine vpsnormal

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: xss(nz,nf,nl),yss(nz,nf,nl),zss(nz,nf,nl)
    integer :: i, j, k
    real :: delp, pval2, r, br_norm, blatp2, blonp2, blons2, mag
    real :: ax1, ay1, az1, a1, blat, blon, qtmp, pvalue
    real :: qp_solve  ! external function

    delp  = .01

!     calculate grid at qs, blons, and ps + delp

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                qtmp   = qs(i,j,k)
                pvalue = ps(i,j,k) * ( 1. + delp )
                r = re*qp_solve(qtmp,pvalue)

                br_norm   = r / re
                blat      = asin ( qs(i,j,k) * br_norm ** 2 ) * rtod
                blon      = blons(i,j,k)

            !           grid in cartesian coordinates

                xss(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                yss(i,j,k) = r * sin( blat * po180 )
                zss(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )
            enddo
        enddo
    enddo

!     direction of e x b velocity in p-direction

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                ax1 = xss(i,j,k) - xs(i,j,k)
                ay1 = yss(i,j,k) - ys(i,j,k)
                az1 = zss(i,j,k) - zs(i,j,k)
                a1  = sqrt ( ax1*ax1 + ay1*ay1 + az1*az1 )
                vpsnx(i,j,k) = ax1 / a1
                vpsny(i,j,k) = ay1 / a1
                vpsnz(i,j,k) = az1 / a1
            enddo
        enddo
    enddo

    return
    end subroutine vpsnormal


!******************************************
!******************************************

!            vhsnormal

!******************************************
!******************************************


    subroutine vhsnormal

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: xss(nz,nf,nl),yss(nz,nf,nl),zss(nz,nf,nl)
    integer :: i, j, k
    real :: delblon, blon2, r, br_norm, blatp2, blonp2, mag
    real :: ax1, ay1, az1, a1, blat, blon

    delblon = .1

!     calculate grid at qs, ps, and blons + delblon

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
            !            qtmp   = qs(i,j,k)
            !            pvalue = ps(i,j,k)
            !            r = re*qp_solve(qtmp,pvalue)
                           
            !            br_norm  = r / re
            !            blat     = asin ( qs(i,j,k) * br_norm ** 2 ) * rtod
            !            blon     = blons(i,j,k) + delblon
                 
                r         = brs(i,j,k)
                blat      = blats(i,j,k)
                blon      = blons(i,j,k) + delblon

            !           grid in cartesian coordinates

                xss(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                yss(i,j,k) = r * sin( blat * po180 )
                zss(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )

            enddo
        enddo
    enddo

!     direction of e x b velocity in h-direction

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                ax1 = xss(i,j,k) - xs(i,j,k)
                ay1 = yss(i,j,k) - ys(i,j,k)
                az1 = zss(i,j,k) - zs(i,j,k)
                a1  = sqrt ( ax1*ax1 + ay1*ay1 + az1*az1 )
                vhsnx(i,j,k) = ax1 / a1
                vhsny(i,j,k) = ay1 / a1
                vhsnz(i,j,k) = az1 / a1
            enddo
        enddo
    enddo

    return
    end subroutine vhsnormal


!******************************************
!******************************************

!            vppnormal

!******************************************
!******************************************


    subroutine vppnormal

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: xpp(nzp1,nfp1,nlp1),ypp(nzp1,nfp1,nlp1),zpp(nzp1,nfp1,nlp1)
    integer :: i, j, k
    real :: delp, pval2, r, br_norm, blatp2, blonp2, mag
    real :: ax1, ay1, az1, a1, blat, blon, qtmp, pvalue
    real :: qp_solve  ! external function

!     calculate grid at qp, blonp, and ps + delp

    delp  = .01 

    do k = 1,nlp1
        do j = 1,nfp1
            do i = 1,nzp1
                qtmp   = qp(i,j,k)
                pvalue = pp(i,j,k) * ( 1. + delp )
                r = re*qp_solve(qtmp,pvalue)

                br_norm   = r / re
                blat      = asin ( qp(i,j,k) * br_norm ** 2 ) * rtod
                blon      = blonp(i,j,k)

            !           grid in cartesian coordinates

                xpp(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                ypp(i,j,k) = r * sin( blat * po180 )
                zpp(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )
            enddo
        enddo
    enddo

!     direction of e x b velocity in p-direction

    do k = 1,nlp1
        do j = 1,nfp1
            do i = 1,nzp1
                ax1 = xpp(i,j,k) - xp(i,j,k)
                ay1 = ypp(i,j,k) - yp(i,j,k)
                az1 = zpp(i,j,k) - zp(i,j,k)
                a1  = sqrt ( ax1*ax1 + ay1*ay1 + az1*az1 )
                vppnx(i,j,k) = ax1 / a1
                vppny(i,j,k) = ay1 / a1
                vppnz(i,j,k) = az1 / a1
            enddo
        enddo
    enddo

    return
    end subroutine vppnormal


!******************************************
!******************************************

!            vhpnormal

!******************************************
!******************************************


    subroutine vhpnormal

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: xpp(nzp1,nfp1,nlp1),ypp(nzp1,nfp1,nlp1),zpp(nzp1,nfp1,nlp1)
    integer :: i, j, k
    real :: delblon, blon2, r, br_norm, blatp2, blonp2, mag
    real :: ax1, ay1, az1, a1, blat, blon

    delblon = .1

!     calculate grid at qs, ps, and blons + delblon

    do k = 1,nlp1
        do j = 1,nfp1
            do i = 1,nzp1
                 
                r         = brp(i,j,k)
                blat      = blatp(i,j,k)
                blon      = blonp(i,j,k) + delblon

            !           grid in cartesian coordinates

                xpp(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                ypp(i,j,k) = r * sin( blat * po180 )
                zpp(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )

            enddo
        enddo
    enddo

!     direction of e x b velocity in h-direction

    do k = 1,nlp1
        do j = 1,nfp1
            do i = 1,nzp1
                ax1 = xpp(i,j,k) - xp(i,j,k)
                ay1 = ypp(i,j,k) - yp(i,j,k)
                az1 = zpp(i,j,k) - zp(i,j,k)
                a1  = sqrt ( ax1*ax1 + ay1*ay1 + az1*az1 )
                vhpnx(i,j,k) = ax1 / a1
                vhpny(i,j,k) = ay1 / a1
                vhpnz(i,j,k) = az1 / a1
            enddo
        enddo
    enddo

    return
    end subroutine vhpnormal


!************************************************

!            vector

!************************************************

! MS 9/15/05
! This subroutine finds the coefficients needed to calculate the
! parallel components of vectors.  The math is shamelessly copied from
! the CTIP paper in the red book.  Inputs are the magnetic coordinates
! of a point, outputs are the coefficients at that point.
     
    subroutine vector(brad,blond,blatd,varg,vathg,vaphig)

    use parameter_mod
    use grid_mod
    use namelist_mod
    use misc_mod

    implicit none  ! fix #1

    real, intent(in)  :: brad, blond, blatd
    real, intent(out) :: varg, vathg, vaphig
    real :: blonr, grad, glat, glon
    real :: thetaprime,cosphiprime,sinphiprime
    real :: sinomega,cosomega,coplat,theta
    real :: arm,athm,ax,ay,az,arp,athp,aphip
     
! Convert blatd and plat to colatitude, blon to radians

    theta  = pie/2. - blatd*po180
    coplat = pie/2. - plat
    blonr  = blond*po180
     
    arm = 2.* cos (theta) /  sqrt ( ( 2. * cos(theta) ) ** 2 + &
    sin(theta)   ** 2 )
    athm =     sin (theta) / sqrt ( ( 2. * cos(theta) ) ** 2 + &
    sin(theta)   ** 2 )
     
    ax = arm * sin(theta)*cos(blonr) + athm*cos(theta)*cos(blonr)
    ay = arm * sin(theta)*sin(blonr) + athm*cos(theta)*sin(blonr)
    az = arm * cos(theta) - athm*sin(theta)
     
! Find geographic coordinates of point

    call btog(brad,blond,blatd,grad,glat,glon)
    glat = pie/2. - glat*po180
    glon = glon * po180
     
! See p.243 of red book

    thetaprime = acos(cos(coplat)*cos(glat) + &
    sin(coplat)*sin(glat)*cos(glon-plon))

! Account for possibility of no tilt.

    if (coplat == 0.) then
        cosphiprime = 1.
        sinphiprime = 0.
    else
        cosphiprime = ( cos(coplat)*cos(thetaprime) - cos(glat) ) / &
        (sin(coplat) * sin(thetaprime))
        sinphiprime = sin(glat) * sin(glon-plon) / sin(thetaprime)
    endif
     
    arp = ax * sin(thetaprime)*cosphiprime + &
    ay*sin(thetaprime)*sinphiprime + &
    az*cos(thetaprime)
    athp = ax * cos(thetaprime)*cosphiprime + &
    ay*cos(thetaprime)*sinphiprime - &
    az*sin(thetaprime)
    aphip = -ax*sinphiprime + ay*cosphiprime
     
    sinomega = -sin(coplat)*sin(glon - plon) / sin(thetaprime)

! MS: This formula is wrong in the red book.  Corrected here.

    cosomega = (cos(coplat) - cos(glat)*cos(thetaprime)) / &
    (sin(thetaprime) * sin(glat))
     
    varg  = arp
    vathg = athp*cosomega + aphip*sinomega

! Because of a typographical error this formula does not appear in
! the red book.

    vaphig = -athp*sinomega + aphip*cosomega
     
    return
    end subroutine vector


!******************************************
!******************************************

!           bdirh (on center of h face)

!******************************************
!******************************************

    subroutine bdirh

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: xm(nz,nf,nlp1),ym(nz,nf,nlp1),zm(nz,nf,nlp1)
    real :: xpp(nz,nf,nlp1),ypp(nz,nf,nlp1),zpp(nz,nf,nlp1)
    integer :: i, j, k
    real :: qfac0, r, blatp2, blonp2, mag
    real :: arg, blat, blon, br_norm, d0, dx, dy, dz, ptmp, qfac, qtmp, qtmpm, qtmpp
    real :: qp_solve  ! external function

    qfac0 = 0.001

    do k = 1,nlp1
        do j = 1,nf
            do i = 1,nz

                qtmp   = .25 * ( qp(i,j,k)   + qp(i+1,j,k)  + &
                qp(i,j+1,k) + qp(i+1,j+1,k)  )
                ptmp   = .25 * ( pp(i,j,k)   + pp(i+1,j,k)  + &
                pp(i,j+1,k) + pp(i+1,j+1,k)  )
                blon   = blonp(i,j,k)

                qfac   = max(abs(.01 * qtmp),qfac0)

                qtmpm  = qtmp - qfac
                r = re*qp_solve(qtmpm,ptmp)
                br_norm   = r / re
                arg       = qtmpm * br_norm ** 2
                blat      = asin ( arg ) * rtod
                xm(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                ym(i,j,k) = r * sin( blat * po180 )
                zm(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )

                qtmpp   = qtmp + qfac
                r = re*qp_solve(qtmpp,ptmp)
                br_norm   = r / re
                arg       = qtmpp * br_norm ** 2
                blat      = asin ( arg ) * rtod
                xpp(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                ypp(i,j,k) = r * sin( blat * po180 )
                zpp(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )

            enddo
        enddo
    enddo

!     calculate bdirh

    do k = 1,nlp1
        do j = 1,nf
            do i = 1,nz
                dx = xpp(i,j,k) - xm(i,j,k)
                dy = ypp(i,j,k) - ym(i,j,k)
                dz = zpp(i,j,k) - zm(i,j,k)
                d0 = sqrt( dx*dx + dy*dy + dz*dz )
                bdirhx(i,j,k) = dx / d0
                bdirhy(i,j,k) = dy / d0
                bdirhz(i,j,k) = dz / d0
            enddo
        enddo
    enddo

    return
    end subroutine bdirh

!******************************************
!******************************************

!           bdirp (on center of p face)

!******************************************
!******************************************

    subroutine bdirp

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: xm(nz,nfp1,nl),ym(nz,nfp1,nl),zm(nz,nfp1,nl)
    real :: xpp(nz,nfp1,nl),ypp(nz,nfp1,nl),zpp(nz,nfp1,nl)
    integer :: i, j, k
    real :: qfac0, r, blatp2, blonp2, mag
    real :: arg, blat, blon, br_norm, d0, dx, dy, dz, ptmp, qfac, qtmp, qtmpm, qtmpp
    real :: qp_solve  ! external function

    qfac0 = 0.001

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nz

                qtmp   = .25 * ( qp(i,j,k)   + qp(i+1,j,k)  + &
                qp(i,j,k+1) + qp(i+1,j,k+1)  )
                ptmp   = .25 * ( pp(i,j,k)   + pp(i+1,j,k)  + &
                pp(i,j,k+1) + pp(i+1,j,k+1)  )
                blon   = .25 * ( blonp(i,j,k)   + blonp(i+1,j,k)  + &
                blonp(i,j,k+1) + blonp(i+1,j,k+1)  )

                qfac   = max(abs(.01 * qtmp),qfac0)

                if ( j /= nfp1 ) then
                    blon  = blons(i,j,k)
                else
                    blon  = blons(i,nf,k)
                endif

                qtmpm  = qtmp - qfac
                r = re*qp_solve(qtmpm,ptmp)
                br_norm   = r / re
                arg       = qtmpm * br_norm ** 2
                blat      = asin ( arg ) * rtod
                xm(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                ym(i,j,k) = r * sin( blat * po180 )
                zm(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )

                qtmpp   = qtmp + qfac
                r = re*qp_solve(qtmpp,ptmp)
                br_norm   = r / re
                arg       = qtmpp * br_norm ** 2
                blat      = asin ( arg ) * rtod
                xpp(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                ypp(i,j,k) = r * sin( blat * po180 )
                zpp(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )

            enddo
        enddo
    enddo

!     calculate bdirp

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nz
                dx = xpp(i,j,k) - xm(i,j,k)
                dy = ypp(i,j,k) - ym(i,j,k)
                dz = zpp(i,j,k) - zm(i,j,k)
                d0 = sqrt( dx*dx + dy*dy + dz*dz )
                bdirpx(i,j,k) = dx / d0
                bdirpy(i,j,k) = dy / d0
                bdirpz(i,j,k) = dz / d0
            enddo
        enddo
    enddo

    return
    end subroutine bdirp

!******************************************
!******************************************

!           bdirs (on center of s face)

!******************************************
!******************************************

    subroutine bdirs

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: xm(nz,nf,nl),ym(nz,nf,nl),zm(nz,nf,nl)
    real :: xpp(nz,nf,nl),ypp(nz,nf,nl),zpp(nz,nf,nl)
    integer :: i, j, k
    real :: qfac0, r, blatp2, blonp2, mag
    real :: arg, blat, blon, br_norm, d0, dx, dy, dz, ptmp, qfac, qtmp, qtmpm, qtmpp
    real :: qp_solve  ! external function

    qfac0 = 0.001

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz

            blon = blons(i,j,k)
            qtmp = qs(i,j,k)
            ptmp = ps(i,j,k)
            qfac = max(abs(.01 * qtmp),qfac0)

                qtmpm  = qtmp - qfac
                r = re*qp_solve(qtmpm,ptmp)
                br_norm   = r / re
                arg       = qtmpm * br_norm ** 2
                blat      = asin ( arg ) * rtod
                xm(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                ym(i,j,k) = r * sin( blat * po180 )
                zm(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )

                qtmpp   = qtmp + qfac
                r = re*qp_solve(qtmpp,ptmp)
                br_norm   = r / re
                arg       = qtmpp * br_norm ** 2
                blat      = asin ( arg ) * rtod
                xpp(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                ypp(i,j,k) = r * sin( blat * po180 )
                zpp(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )

            enddo
        enddo
    enddo

!     calculate bdirs

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                dx = xpp(i,j,k) - xm(i,j,k)
                dy = ypp(i,j,k) - ym(i,j,k)
                dz = zpp(i,j,k) - zm(i,j,k)
                d0 = sqrt( dx*dx + dy*dy + dz*dz )
                bdirsx(i,j,k) = dx / d0
                bdirsy(i,j,k) = dy / d0
                bdirsz(i,j,k) = dz / d0
            enddo
        enddo
    enddo


    return
    end subroutine bdirs

!******************************************
!******************************************

!            facesp

!******************************************
!******************************************


    subroutine facesp

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: xpp(nzp1,nf,nlp1),ypp(nzp1,nf,nlp1),zpp(nzp1,nf,nlp1)
    real :: bmtmp(nzp1,nf,nlp1)
    integer :: i, j, k, nzh
    real :: r, br_norm, blatp2, blonp2, brs2, blats2, blons2, mag
    real :: dpphi, dptheta, blat, blon, delx, dely, delz, ptmp, qtmp
    real :: qp_solve  ! external function

!     calculate grid on p face centered on s
!     to obtain blaths, delhs, bmhf, and bmsf

    nzh = nz / 2

    do k = 1,nlp1
        do j = 1,nf
            do i = 1,nzp1
                qtmp   = .5 * ( qp(i,j,k) + qp(i,j+1,k) )
                ptmp   = .5 * ( pp(i,j,k) + pp(i,j+1,k) )
                r = re*qp_solve(qtmp,ptmp)
                br_norm   = r / re
                blat      = asin ( qtmp * br_norm ** 2 ) * rtod
                blon = blonp(i,j,k)
!                theta = acos ( qtmp * br_norm ** 2 )
!                bmtmp(i,j,k) = sqrt ( 1. + 3.*cos(theta)**2 ) / br_norm ** 3
            if (i .le. nz .and. k.le.nl) bmtmp(i,j,k) = bms(i,j,k)
            if (i .eq. nzp1) bmtmp(i,j,k) = bmtmp(nz,j,k)
            if (k .eq. nlp1) bmtmp(i,j,k) = bmtmp(i,j,nl)
                xpp(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                ypp(i,j,k) = r * sin( blat * po180 )
                zpp(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )
                if ( i == nzh) then
                    blatsp(j,k) = blat
                    blonsp(j,k) = blon
                    bradsp(j,k) = r
                endif
            enddo
        enddo
    enddo

!     calculate bmhf

    do k = 1,nlp1
        do j = 1,nf
            do i = 1,nz
                bmhf(i,j,k) = .5 * ( bmtmp(i,j,k) + bmtmp(i+1,j,k) )
            enddo
        enddo
    enddo

!     calculate bmsf

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                bmsf(i,j,k) = .5 * (bmtmp(i,j,k) + bmtmp(i,j,k+1) )
            enddo
        enddo
    enddo

    do k = 1,nl
        do j = 1,nf
            i = nz + 1
            bmsf(i,j,k) = bmsf(nz,j,k)
        enddo
    enddo
          
!     calculate delhs

    do k = 1,nl
        do j = 1,nf
            do i = 1,nzp1
                delx  = xpp(i,j,k+1) - xpp(i,j,k)
                dely  = ypp(i,j,k+1) - ypp(i,j,k)
                delz  = zpp(i,j,k+1) - zpp(i,j,k)
                dpphi = sqrt( delx*delx + dely*dely + delz*delz )
                delhs(i,j,k) = dpphi * 1.e5 ! convert to cm
            enddo
        enddo
    enddo

    return
    end subroutine facesp

!******************************************
!******************************************

!            facesh

!******************************************
!******************************************


    subroutine facesh

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: xpp(nzp1,nfp1,nlp1),ypp(nzp1,nfp1,nlp1),zpp(nzp1,nfp1,nlp1)
    integer :: i, j, k, nzh
    real :: r, br_norm, blatp2, blonp2, brs2, blats2, blons2, mag
    real :: dpphi, dptheta, blat, blon, delx, dely, delz, ptmp, qtmp
    real :: qp_solve  ! external function

!     calculate grid on h face centered on s
!     to obtain blatsh and delps

    nzh = nz / 2

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nzp1
                qtmp   = .5 * ( qp(i,j,k) + qp(i,j,k+1) )
                ptmp   = .5 * ( pp(i,j,k) + pp(i,j,k+1) )
                r = re*qp_solve(qtmp,ptmp)
                br_norm   = r / re
                blat      = asin ( qtmp * br_norm ** 2 ) * rtod
                blon      = blonp(i,j,k)
                xpp(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                ypp(i,j,k) = r * sin( blat * po180 )
                zpp(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )
                if ( i == nzh ) then
                    blatsh(j,k) = blat
                    blonsh(j,k) = blon
                    bradsh(j,k) = r
                endif
            enddo
        enddo
    enddo

!     calculate delps

    do k = 1,nl
        do j = 1,nf
            do i = 1,nzp1
                delx  = xpp(i,j+1,k) - xpp(i,j,k)
                dely  = ypp(i,j+1,k) - ypp(i,j,k)
                delz  = zpp(i,j+1,k) - zpp(i,j,k)
                dpphi = sqrt( delx*delx + dely*dely + delz*delz )
                delps(i,j,k) = dpphi * 1.e5
            enddo
        enddo
    enddo


    return
    end subroutine facesh
      
!******************************************
!******************************************

!            facess

!******************************************
!******************************************


    subroutine facess

    use parameter_mod
    use grid_mod
    use exb_mod

    implicit none  ! fix #1

    real :: xpp(nz,nfp1,nlp1),ypp(nz,nfp1,nlp1),zpp(nz,nfp1,nlp1)
    real :: bmtmp(nz,nfp1,nlp1)
    integer :: i, j, k, nzh
    real :: r, br_norm, blatp2, blonp2, brs2, blats2, blons2, mag
    real :: dpphi, dptheta, blat, blon, delx, dely, delz, ptmp, qtmp
    real :: qp_solve  ! external function

!     calculate grid on s face centered on s
!     to obtain blatss, delph, delhp, bmpf

    nzh = nz / 2

    do k = 1,nlp1
        do j = 1,nfp1
            do i = 1,nz
                qtmp   = 0.5 * ( qp(i,j,k) + qp(i+1,j,k) )
                ptmp   = 0.5 * ( pp(i,j,k) + pp(i+1,j,k) )
                r = re*qp_solve(qtmp,ptmp)
                br_norm   = r / re
                blat      = asin ( qtmp * br_norm ** 2 ) * rtod
                blon = blonp(i,j,k)
!                theta = acos ( qtmp * br_norm ** 2 )
!                bmtmp(i,j,k) = sqrt ( 1. + 3.*cos(theta)**2 ) / br_norm ** 3
            if (j .le. nf .and. k.le.nl) bmtmp(i,j,k) = bms(i,j,k)
            if (j .eq. nfp1) bmtmp(i,j,k) = bmtmp(i,nf,k)
            if (k .eq. nlp1) bmtmp(i,j,k) = bmtmp(i,j,nl)
                xpp(i,j,k) = r * cos( blat * po180 ) * sin( blon * po180 )
                ypp(i,j,k) = r * sin( blat * po180 )
                zpp(i,j,k) = r * cos( blat * po180 ) * cos( blon * po180 )
                if ( i == nzh) then
                    blatss(j,k) = blat
                    blonss(j,k) = blon
                    bradss(j,k) = r
                endif
            enddo
        enddo
    enddo

!     calculate bmpf

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nz
                bmpf(i,j,k) = .5 * (bmtmp(i,j,k) + bmtmp(i,j,k+1) )
            enddo
        enddo
    enddo

!     calculate delph

    do k = 1,nlp1
        do j = 1,nf
            do i = 1,nz
                delx  = xpp(i,j+1,k) - xpp(i,j,k)
                dely  = ypp(i,j+1,k) - ypp(i,j,k)
                delz  = zpp(i,j+1,k) - zpp(i,j,k)
                dpphi = sqrt( delx*delx + dely*dely + delz*delz )
                delph(i,j,k) = dpphi * 1.e5 ! convert to cm
                ephx(i,j,k)  = delx/dpphi
                ephy(i,j,k)  = dely/dpphi
                ephz(i,j,k)  = delz/dpphi
            enddo
        enddo
    enddo

!     calculate delhp

    do k = 1,nl
        do j = 1,nfp1
            do i = 1,nz
                delx  = xp(i,j,k+1) - xp(i,j,k)
                dely  = yp(i,j,k+1) - yp(i,j,k)
                delz  = zp(i,j,k+1) - zp(i,j,k)
                dpphi = sqrt( delx*delx + dely*dely + delz*delz )
                delhp(i,j,k) = dpphi * 1.e5 ! convert to cm
                ehpx(i,j,k)  = delx/dpphi
                ehpy(i,j,k)  = dely/dpphi
                ehpz(i,j,k)  = delz/dpphi
            enddo
        enddo
    enddo

    return
    end subroutine facess


!******************************************
!******************************************

!            gstheta

!******************************************
!******************************************


    subroutine gstheta

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: xm(nz,nf,nl),ym(nz,nf,nl),zm(nz,nf,nl)
    real :: xpp(nz,nf,nl),ypp(nz,nf,nl),zpp(nz,nf,nl)
    integer :: i, j, k
    real :: delglat, glat2, glon2, r2, mag
    real :: d0, dx, dy, dz, glattmp, glon, grad, x, y, z

    delglat = 0.1

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                glattmp = glats(i,j,k) - delglat
                glon    = glons(i,j,k) + 72.
                grad    = grs(i,j,k)
!                call gtob(brad,blonr,blatr,grad,glattmp,glon)
!                call btog_xyz(brad,blonr,blatr,x,y,z)
                z = grad *cos ( glattmp*po180 ) * cos ( glon*po180 )
                x = grad *cos ( glattmp*po180 ) * sin ( glon*po180 )
                y = grad *sin ( glattmp*po180 )
                xm(i,j,k) = x
                ym(i,j,k) = y
                zm(i,j,k) = z

                glattmp = glats(i,j,k) + delglat
                glon    = glons(i,j,k) + 72.
                grad    = grs(i,j,k)
!                call gtob(brad,blonr,blatr,grad,glattmp,glon)
!                call btog_xyz(brad,blonr,blatr,x,y,z)
                z = grad *cos ( glattmp*po180 ) * cos ( glon*po180 )
                x = grad *cos ( glattmp*po180 ) * sin ( glon*po180 )
                y = grad *sin ( glattmp*po180 )
                xpp(i,j,k) = x
                ypp(i,j,k) = y
                zpp(i,j,k) = z
            enddo
        enddo
    enddo

!     calculate gstheta

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                dx = xpp(i,j,k) - xm(i,j,k)
                dy = ypp(i,j,k) - ym(i,j,k)
                dz = zpp(i,j,k) - zm(i,j,k)
                d0 = sqrt( dx*dx + dy*dy + dz*dz )
                gsthetax(i,j,k) = dx / d0
                gsthetay(i,j,k) = dy / d0
                gsthetaz(i,j,k) = dz / d0
            enddo
        enddo
    enddo

    return
    end subroutine gstheta

!******************************************
!******************************************

!            gsphi

!******************************************
!******************************************


    subroutine gsphi

    use parameter_mod
    use grid_mod

    implicit none  ! fix #1

    real :: xm(nz,nf,nl),ym(nz,nf,nl),zm(nz,nf,nl)
    real :: xpp(nz,nf,nl),ypp(nz,nf,nl),zpp(nz,nf,nl)
    integer :: i, j, k
    real :: delglon, glat2, glon2, r2, mag
    real :: d0, dx, dy, dz, glat, glontmp, grad, x, y, z

    delglon = 0.1

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                glontmp = glons(i,j,k) - delglon + 72.
                glat    = glats(i,j,k)
                grad    = grs(i,j,k)
!                call gtob(brad,blonr,blatr,grad,glat,glontmp)
!                call btog_xyz(brad,blonr,blatr,x,y,z)
                z = grad *cos ( glat*po180 ) * cos ( glontmp*po180 )
                x = grad *cos ( glat*po180 ) * sin ( glontmp*po180 )
                y = grad *sin ( glat*po180 )
                xm(i,j,k) = x
                ym(i,j,k) = y
                zm(i,j,k) = z

                glontmp = glons(i,j,k) + delglon + 72.
                glat    = glats(i,j,k)
                grad    = grs(i,j,k)
!                call gtob(brad,blonr,blatr,grad,glat,glontmp)
!                call btog_xyz(brad,blonr,blatr,x,y,z)
                z = grad *cos ( glat*po180 ) * cos ( glontmp*po180 )
                x = grad *cos ( glat*po180 ) * sin ( glontmp*po180 )
                y = grad *sin ( glat*po180 )
                xpp(i,j,k) = x
                ypp(i,j,k) = y
                zpp(i,j,k) = z
            enddo
        enddo
    enddo

!     calculate gsphi

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                dx = xpp(i,j,k) - xm(i,j,k)
                dy = ypp(i,j,k) - ym(i,j,k)
                dz = zpp(i,j,k) - zm(i,j,k)
                d0 = sqrt( dx*dx + dy*dy + dz*dz )
                gsphix(i,j,k) = dx / d0
                gsphiy(i,j,k) = dy / d0
                gsphiz(i,j,k) = dz / d0
            enddo
        enddo
    enddo

    return
    end subroutine gsphi


!******************************************
!******************************************
     
!            gsr
     
!******************************************
!******************************************
     
     
    subroutine gsr

    use parameter_mod
    use grid_mod
    use message_passing_mod

    implicit none  ! fix #1

    real :: xm(nz,nf,nl),ym(nz,nf,nl),zm(nz,nf,nl)
    real :: xpp(nz,nf,nl),ypp(nz,nf,nl),zpp(nz,nf,nl)
    integer :: i, j, k
    real :: delr, glat2, glon2, r2, mag
    real :: d0, dx, dy, dz, glat, glon, gradtmp, x, y, z

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                delr    = .01 * grs(i,j,k)
                glat    = glats(i,j,k)
                glon    = glons(i,j,k) + 72.
                gradtmp = grs(i,j,k) - delr
!                call gtob(brad,blonr,blatr,gradtmp,glat,glon)
!                call btog_xyz(brad,blonr,blatr,x,y,z)
                z = gradtmp * cos ( glat*po180 ) * cos ( glon*po180 )
                x = gradtmp * cos ( glat*po180 ) * sin ( glon*po180 )
                y = gradtmp * sin ( glat*po180 )
                xm(i,j,k) = x
                ym(i,j,k) = y
                zm(i,j,k) = z

                glat    = glats(i,j,k)
                glon    = glons(i,j,k) + 72.
                gradtmp = grs(i,j,k) + delr
!                call gtob(brad,blonr,blatr,gradtmp,glat,glon)
!                call btog_xyz(brad,blonr,blatr,x,y,z)
                z = gradtmp * cos ( glat*po180 ) * cos ( glon*po180 )
                x = gradtmp * cos ( glat*po180 ) * sin ( glon*po180 )
                y = gradtmp * sin ( glat*po180 )
                xpp(i,j,k) = x
                ypp(i,j,k) = y
                zpp(i,j,k) = z
            enddo
        enddo
    enddo
     
!     calculate gsr
     
    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                dx = xpp(i,j,k) - xm(i,j,k)
                dy = ypp(i,j,k) - ym(i,j,k)
                dz = zpp(i,j,k) - zm(i,j,k)
                d0 = sqrt( dx*dx + dy*dy + dz*dz )
                gsrx(i,j,k) = dx / d0
                gsry(i,j,k) = dy / d0
                gsrz(i,j,k) = dz / d0
            enddo
        enddo
    enddo
     
!!$    if (taskid == 1) then
!!$        open(188,file='gsr.dat',form='unformatted')
!!$        write(188) gsrx,gsry,gsrz
!!$        close(188)
!!$    endif
     
    return
    end subroutine gsr


!******************************************
!******************************************

!            blonp0a

!******************************************
!******************************************

    subroutine blonp0a

    use parameter_mod
    use grid_mod
    use misc_mod
    use message_passing_mod

    implicit none  ! fix #1

    real :: theta(nltp1)
    integer :: i,n,nn
    real :: theta_min, theta_max, dtheta

    theta_min = 0.
    theta_max = 360.
    dtheta    = theta_max/real(nlt)  ! fix #9: float->real
                        
    do i = 1,nlt
      theta(i) = float(i-1) * dtheta
     enddo
     
    theta(nlt+1) = 360.
     
    do n = 1,nltp1
        nn = n + 1
        blonp0t(nn) = theta(n)
    enddo
     
    blonp0t(1)     = blonp0t(nlt+1) - blonp0t(nlt+2)
    blonp0t(nlt+3) = blonp0t(nlt+2) + blonp0t(3)

!    if ( taskid == 0) then
!        do i = 1,nlt+3
!            print *,i,blonp0t(i)
!        enddo
!    endif

    return
    end subroutine blonp0a
               

! *********************

!     smoothp

! *********************

    subroutine smoothp(finout)

    use parameter_mod

    implicit none  ! fix #1

    real, intent(inout) :: finout(nfp1)
    real :: tempz(nfp1)
    integer :: i, im1, ip1


! This is the binomial filter (in x space) as described in
! Birdsall appendix C. no compensation

     
! do smoothp in the p direction
     
    do i = 1,nfp1
        tempz(i) = finout(i)
    enddo

    do i = 3,nf-2
        ip1 = i + 1
        im1 = i - 1
        tempz(i) = .25*(finout(im1) +2.*finout(i) &
        +finout(ip1))
    enddo
    do i = 2,nf
        finout(i) = tempz(i)
    enddo
      
    return
    end subroutine smoothp


! *********************

!     smooths

! *********************

    subroutine smooths(finout)

    use parameter_mod

    implicit none  ! fix #1

    real, intent(inout) :: finout(nzp1)
    real :: tempz(nzp1)
    integer :: i, im1, ip1


! This is the binomial filter (in x space) as described in
! Birdsall appendix C. no compensation

     
! do smoothz in the z direction
     
    do i = 2,nz
        ip1 = i +1
        im1 = i -1
        tempz(i) = .25*(finout(im1) +2.*finout(i) &
        +finout(ip1))
    enddo
    do i = 2,nz
        finout(i) = tempz(i)
    enddo
      
    return
    end subroutine smooths


