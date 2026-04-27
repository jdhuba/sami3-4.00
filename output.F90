! ---------------------------------------------------------------
! output_units_mod
!   Named integer parameters for every output unit number used by
!   open_u and output.  Fixes issue #7: previously all unit numbers
!   were magic literals scattered across two subroutines with no
!   central registry.
! ---------------------------------------------------------------
module output_units_mod
    implicit none

    ! Formatted
    integer, parameter :: U_TIME    = 70

    ! Electron / ion temperatures
    integer, parameter :: U_TE      = 75

    ! Height-integrated conductances (2-D, nf x nlt)
    integer, parameter :: U_HIPC    = 93,  U_HIHC    = 94
    integer, parameter :: U_SIGMA0  = 951, U_SIGMAP  = 95
    integer, parameter :: U_SIGMAH  = 96,  U_SIGMAPIC = 97, U_SIGMAHIC = 98

    ! Ion densities (per species)
    integer, parameter :: U_DENI(7) = [711,712,713,714,715,716,717]
    integer, parameter :: U_DENE    = 1718

    ! Ion temperatures (per species)
    integer, parameter :: U_TI(7)   = [811,812,813,814,815,816,817]

    ! Ion field-aligned velocities (per species)
    integer, parameter :: U_VSI(7)  = [911,912,913,914,915,916,917]

    ! Ion loss rates (per species)
    integer, parameter :: U_LOSS(7) = [1811,1812,1813,1814,1815,1816,1817]

    ! Ion collision frequencies (per species)
    integer, parameter :: U_NUIN(7) = [1911,1912,1913,1914,1915,1916,1917]

    ! Neutral densities (per species)
    integer, parameter :: U_DENN(7) = [1711,1712,1713,1714,1715,1716,1717]

    ! RHS eigenvector
    integer, parameter :: U_RHSEGV  = 569

    ! Neutral winds / currents
    integer, parameter :: U_VNQ     = 196, U_VNP     = 197, U_VNPHI   = 198
    integer, parameter :: U_JP      = 201, U_JPHI    = 202

    ! Height-integrated conductance components (2-D)
    integer, parameter :: U_HIPCPU  = 491, U_HIPCPHI = 492, U_HIHCM   = 493
    integer, parameter :: U_HIDPV   = 494, U_HIDPHIV = 495
    integer, parameter :: U_HIDPG   = 496, U_HIDPHIG = 497
    integer, parameter :: U_PHI     = 498

    ! Diagnostic / velocity fields
    integer, parameter :: U_U1P     = 384, U_U2S     = 385, U_U3H     = 386
    integer, parameter :: U_U1      = 84,  U_U2      = 85,  U_U3      = 86
    integer, parameter :: U_U4      = 87,  U_U5      = 88,  U_GP      = 89

end module output_units_mod


!******************************************
!******************************************

!             open_u

!******************************************
!******************************************

    subroutine open_u

! Open all output files (unformatted, except time.dat).
! Unit numbers are drawn from output_units_mod to stay consistent
! with the write statements in subroutine output.

    use output_units_mod
    implicit none

    open ( unit=U_TIME,    file='time.dat',      form='formatted'   )
    open ( unit=U_TE,      file='teu.dat',        form='unformatted' )
    open ( unit=U_HIPC,    file='hipcu.dat',      form='unformatted' )
    open ( unit=U_HIHC,    file='hihcu.dat',      form='unformatted' )
    open ( unit=U_SIGMA0,  file='sigma0u.dat',    form='unformatted' )
    open ( unit=U_SIGMAP,  file='sigmapu.dat',    form='unformatted' )
    open ( unit=U_SIGMAH,  file='sigmahu.dat',    form='unformatted' )
    open ( unit=U_SIGMAPIC,file='sigmapicu.dat',  form='unformatted' )
    open ( unit=U_SIGMAHIC,file='sigmahicu.dat',  form='unformatted' )

    open ( unit=U_DENI(1), file='deni1u.dat',     form='unformatted' )
    open ( unit=U_DENI(2), file='deni2u.dat',     form='unformatted' )
    open ( unit=U_DENI(3), file='deni3u.dat',     form='unformatted' )
    open ( unit=U_DENI(4), file='deni4u.dat',     form='unformatted' )
    open ( unit=U_DENI(5), file='deni5u.dat',     form='unformatted' )
    open ( unit=U_DENI(6), file='deni6u.dat',     form='unformatted' )
    open ( unit=U_DENI(7), file='deni7u.dat',     form='unformatted' )
    open ( unit=U_DENE,    file='deneu.dat',       form='unformatted' )

    open ( unit=U_TI(1),   file='ti1u.dat',       form='unformatted' )
    open ( unit=U_TI(2),   file='ti2u.dat',       form='unformatted' )
    open ( unit=U_TI(3),   file='ti3u.dat',       form='unformatted' )
    open ( unit=U_TI(4),   file='ti4u.dat',       form='unformatted' )
    open ( unit=U_TI(5),   file='ti5u.dat',       form='unformatted' )
    open ( unit=U_TI(6),   file='ti6u.dat',       form='unformatted' )
    open ( unit=U_TI(7),   file='ti7u.dat',       form='unformatted' )

    open ( unit=U_VSI(1),  file='vsi1u.dat',      form='unformatted' )
    open ( unit=U_VSI(2),  file='vsi2u.dat',      form='unformatted' )
    open ( unit=U_VSI(3),  file='vsi3u.dat',      form='unformatted' )
    open ( unit=U_VSI(4),  file='vsi4u.dat',      form='unformatted' )
    open ( unit=U_VSI(5),  file='vsi5u.dat',      form='unformatted' )
    open ( unit=U_VSI(6),  file='vsi6u.dat',      form='unformatted' )
    open ( unit=U_VSI(7),  file='vsi7u.dat',      form='unformatted' )

    open ( unit=U_LOSS(1), file='loss1u.dat',      form='unformatted' )
    open ( unit=U_LOSS(2), file='loss2u.dat',      form='unformatted' )
    open ( unit=U_LOSS(3), file='loss3u.dat',      form='unformatted' )
    open ( unit=U_LOSS(4), file='loss4u.dat',      form='unformatted' )
    open ( unit=U_LOSS(5), file='loss5u.dat',      form='unformatted' )
    open ( unit=U_LOSS(6), file='loss6u.dat',      form='unformatted' )
    open ( unit=U_LOSS(7), file='loss7u.dat',      form='unformatted' )

    open ( unit=U_NUIN(1), file='nuin1u.dat',      form='unformatted' )
    open ( unit=U_NUIN(2), file='nuin2u.dat',      form='unformatted' )
    open ( unit=U_NUIN(3), file='nuin3u.dat',      form='unformatted' )
    open ( unit=U_NUIN(4), file='nuin4u.dat',      form='unformatted' )
    open ( unit=U_NUIN(5), file='nuin5u.dat',      form='unformatted' )
    open ( unit=U_NUIN(6), file='nuin6u.dat',      form='unformatted' )
    open ( unit=U_NUIN(7), file='nuin7u.dat',      form='unformatted' )

    open ( unit=U_DENN(1), file='denn1u.dat',      form='unformatted' )
    open ( unit=U_DENN(2), file='denn2u.dat',      form='unformatted' )
    open ( unit=U_DENN(3), file='denn3u.dat',      form='unformatted' )
    open ( unit=U_DENN(4), file='denn4u.dat',      form='unformatted' )
    open ( unit=U_DENN(5), file='denn5u.dat',      form='unformatted' )
    open ( unit=U_DENN(6), file='denn6u.dat',      form='unformatted' )
    open ( unit=U_DENN(7), file='denn7u.dat',      form='unformatted' )

    open ( unit=U_RHSEGV,  file='rhsegv.dat',      form='unformatted' )

    open ( unit=U_VNQ,     file='vnqu.dat',         form='unformatted' )
    open ( unit=U_VNP,     file='vnpu.dat',         form='unformatted' )
    open ( unit=U_VNPHI,   file='vnphiu.dat',       form='unformatted' )
    open ( unit=U_JP,      file='jpu.dat',          form='unformatted' )
    open ( unit=U_JPHI,    file='jphiu.dat',        form='unformatted' )

    open ( unit=U_HIPCPU,  file='hipcpu.dat',       form='unformatted' )
    open ( unit=U_HIPCPHI, file='hipcphiu.dat',     form='unformatted' )
    open ( unit=U_HIHCM,   file='hihcmu.dat',       form='unformatted' )
    open ( unit=U_HIDPV,   file='hidpvu.dat',        form='unformatted' )
    open ( unit=U_HIDPHIV, file='hidphivu.dat',      form='unformatted' )
    open ( unit=U_HIDPG,   file='hidpgu.dat',        form='unformatted' )
    open ( unit=U_HIDPHIG, file='hidphigu.dat',      form='unformatted' )
    open ( unit=U_PHI,     file='phiu.dat',          form='unformatted' )

! diagnostic files (unformatted)

    open ( unit=U_U1P,     file='u1pu.dat',          form='unformatted' )
    open ( unit=U_U2S,     file='u2su.dat',          form='unformatted' )
    open ( unit=U_U3H,     file='u3hu.dat',          form='unformatted' )
    open ( unit=U_U1,      file='u1u.dat',           form='unformatted' )
    open ( unit=U_U2,      file='u2u.dat',           form='unformatted' )
    open ( unit=U_U3,      file='u3u.dat',           form='unformatted' )
    open ( unit=U_U4,      file='u4u.dat',           form='unformatted' )
    open ( unit=U_U5,      file='u5u.dat',           form='unformatted' )
    open ( unit=U_GP,      file='gpu.dat',           form='unformatted' )

    end subroutine open_u


!******************************************
!******************************************

!             write_ion_slices
!
!   Helper: split the 4th dimension of a (nz,nf,nlt,nsp) array into
!   nsp separate 3D slabs and write each to the corresponding unit.
!   Fixes issue #6: the allocate/copy/write/deallocate pattern was
!   duplicated six times in the original output subroutine.

!******************************************
!******************************************

    subroutine write_ion_slices(arr4d, nsp, units, nz_in, nf_in, nlt_in)

    implicit none

    integer, intent(in) :: nsp, nz_in, nf_in, nlt_in
    real,    intent(in) :: arr4d(nz_in, nf_in, nlt_in, nsp)
    integer, intent(in) :: units(nsp)

    real, allocatable :: slab(:,:,:)
    integer :: s, i, j, k

    allocate( slab(nz_in, nf_in, nlt_in) )

    do s = 1, nsp
        do k = 1, nlt_in
            do j = 1, nf_in
                do i = 1, nz_in
                    slab(i,j,k) = arr4d(i,j,k,s)
                enddo
            enddo
        enddo
        write(units(s)) slab
    enddo

    deallocate(slab)

    end subroutine write_ion_slices


!******************************************
!******************************************

!             output

!******************************************
!******************************************

    subroutine output ( hr,ntm,istep,phi,denit,dennt,vsit,nuintt,&
    losstt,sumvsit, tit,ut,vt,vpit,tet,tnt,u1t, &
    u2t,u3t,u4t,u5t,vnqt,vnpt,vnphit,jpt,jphit, &
    u1pt,u2st,u3ht,sigmapict,sigmahict, &
    sigmapt,sigmaht,sigma0t,gpt)

    use parameter_mod
    use variable_mod
    use namelist_mod
    use time_mod
    use conductance_mod
    use output_units_mod

    implicit none

    ! Arguments
    real,    intent(in) :: hr
    integer, intent(in) :: ntm, istep
    real,    intent(in) :: phi(nnx,nny)
    real,    intent(in) :: denit(nz,nf,nlt,nion)
    real,    intent(in) :: dennt(nz,nf,nlt,nneut)  ! fix: was (nion), correct dim is nneut
    real,    intent(in) :: vsit(nz,nf,nlt,nion)
    real,    intent(in) :: nuintt(nz,nf,nlt,nion),losstt(nz,nf,nlt,nion)
    real,    intent(in) :: sumvsit(nz,nf,nlt,nion)
    real,    intent(in) :: tet(nz,nf,nlt),tit(nz,nf,nlt,nion),tnt(nz,nf,nlt)
    real,    intent(in) :: ut(nz,nf,nlt),vt(nz,nf,nlt),vpit(nz,nf,nlt)
    real,    intent(in) :: u1t(nz,nf,nlt),u2t(nz,nf,nlt)
    real,    intent(in) :: u3t(nz,nf,nlt),u4t(nz,nf,nlt),u5t(nz,nf,nlt)
    real,    intent(in) :: vnqt(nz,nf,nlt),vnpt(nz,nf,nlt),vnphit(nz,nf,nlt)
    real,    intent(in) :: jpt(nz,nf,nlt),jphit(nz,nf,nlt)
    real,    intent(in) :: u1pt(nz,nf,nlt),u2st(nz,nf,nlt),u3ht(nz,nf,nlt)
    real,    intent(in) :: sigmapict(nz,nf,nlt),sigmahict(nz,nf,nlt)
    real,    intent(in) :: sigmapt(nz,nf,nlt),sigmaht(nz,nf,nlt)
    real,    intent(in) :: sigma0t(nz,nf,nlt),gpt(nz,nf,nlt)

    ! Local variables
    integer :: i, j, k, ni, s
    real    :: hr24, totsec, thr, tmin, tsec
    integer :: nthr, ntmin, ntsec

    ! Electron density accumulation array and restart slab
    real, allocatable :: denet(:,:,:)
    real, allocatable :: rst_slab(:,:,:)

    ! ------------------------------------------------------------------
    ! Time-of-day breakdown
    ! ------------------------------------------------------------------
    hr24   = mod(hr, 24.0)
    totsec = hr24 * 3600.0
    thr    = totsec / 3600.0
    nthr   = int(thr)
    tmin   = (thr - nthr) * 60.0
    ntmin  = int(mod(tmin, 60.0))
    tsec   = (tmin - ntmin) * 60.0
    ntsec  = int(tsec)

    print *,'istep = ',istep,' ntm = ',ntm
    print *,' hr = ',hr,' dt = ',dt

    write (U_TIME,100) ntm,nthr,ntmin,ntsec,hr
    flush(U_TIME)

    ! ------------------------------------------------------------------
    ! Ion densities -- split 4D array into per-species 3D files
    ! ------------------------------------------------------------------
    call write_ion_slices(denit, nion, U_DENI, nz, nf, nlt)

    ! Restart file for ion densities (opened/closed here per call)
    open (144, file='denit_cg.rst', form='unformatted')
    allocate( rst_slab(nz,nf,nlt) )
    do s = 1, nion
        do k = 1,nlt; do j = 1,nf; do i = 1,nz
            rst_slab(i,j,k) = denit(i,j,k,s)
        enddo; enddo; enddo
        write(144) rst_slab
    enddo
    deallocate(rst_slab)
    close(144)

    ! Electron density = sum over active ion species
    allocate( denet(nz,nf,nlt) )
    do k = 1, nlt
        do j = 1, nf
            do i = 1, nz
                denet(i,j,k) = 0.0          ! fix: was integer 0
                do ni = nion1, nion2
                    denet(i,j,k) = denet(i,j,k) + denit(i,j,k,ni)
                enddo
            enddo
        enddo
    enddo

    print *,'hr,dene ',hr,denet(nz/2,35,nlt/2)
    write(U_DENE) denet
    deallocate(denet)

    ! ------------------------------------------------------------------
    ! Neutral densities
    ! ------------------------------------------------------------------
    call write_ion_slices(dennt, nneut, U_DENN, nz, nf, nlt)

    ! ------------------------------------------------------------------
    ! Ion temperatures
    ! ------------------------------------------------------------------
    call write_ion_slices(tit, nion, U_TI, nz, nf, nlt)

    write(U_TE) tet

    print *,'hr,te  ',tet(nz/2,50,nlt/2)

    ! ------------------------------------------------------------------
    ! Ion field-aligned velocities
    ! ------------------------------------------------------------------
    call write_ion_slices(vsit, nion, U_VSI, nz, nf, nlt)

    ! ------------------------------------------------------------------
    ! Ion collision frequencies
    ! ------------------------------------------------------------------
    call write_ion_slices(nuintt, nion, U_NUIN, nz, nf, nlt)

    ! ------------------------------------------------------------------
    ! Ion loss rates
    ! ------------------------------------------------------------------
    call write_ion_slices(losstt, nion, U_LOSS, nz, nf, nlt)

    ! ------------------------------------------------------------------
    ! Velocity / E x B drift fields
    ! ------------------------------------------------------------------
    write(U_U1P)   u1pt
    write(U_U2S)   u2st
    write(U_U3H)   u3ht
    write(U_U1)    u1t
    write(U_U2)    u2t
    write(U_U3)    u3t
    write(U_U4)    u4t
    write(U_U5)    u5t
    write(U_GP)    gpt

    ! ------------------------------------------------------------------
    ! Conductances
    ! ------------------------------------------------------------------
    write(U_SIGMA0)   sigma0t
    write(U_SIGMAP)   sigmapt
    write(U_SIGMAH)   sigmaht
    write(U_SIGMAPIC) sigmapict
    write(U_SIGMAHIC) sigmahict

    ! ------------------------------------------------------------------
    ! Neutral winds / current densities
    ! ------------------------------------------------------------------
    write(U_VNQ)   vnqt
    write(U_VNP)   vnpt
    write(U_VNPHI) vnphit
    write(U_JP)    jpt
    write(U_JPHI)  jphit

    ! ------------------------------------------------------------------
    ! Height-integrated conductance components and potential
    ! (2-D arrays from conductance_mod: nf x nlt)
    ! ------------------------------------------------------------------
    write(U_HIPCPU)  hipcpt
    write(U_HIPCPHI) hipcphit
    write(U_HIHCM)   hihcmt
    write(U_HIDPV)   hidpvt
    write(U_HIDPHIV) hidphivt
    write(U_HIDPG)   hidpgt
    write(U_HIDPHIG) hidphigt
    write(U_PHI)     phi

    100 format(1x,4i6,1p1e14.4)

    end subroutine output
