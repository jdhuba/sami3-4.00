module exb_mod

    use parameter_mod

    implicit none

    ! E x B drift velocities (face-centred)
    real :: vexbs(nzp1,nf,nl), vexbp(nz,nfp1,nl), vexbh(nz,nf,nlp1)

    ! E x B drift velocities derived from electric potential
    real :: vexbs_phi(nzp1,nf,nl), vexbp_phi(nz,nfp1,nl), &
            vexbh_phi(nz,nf,nlp1)

    ! Cell-centred output drift velocities
    ! NOTE: these hold interior cell values only; face values at
    !       nfp1/nzp1/nlp1 are intentionally not copied from the
    !       face-centred vexb* arrays.
    real :: u1p(nz,nf,nl), u2s(nz,nf,nl), u3h(nz,nf,nl)

    ! Electric field components (s-face and p/h-face)
    real :: eps(nzp1,nf,nl),  eph(nz,nf,nlp1)
    real :: ehs(nzp1,nf,nl),  ehp(nz,nfp1,nl)

    ! Unit-vector components for cross products
    real :: ehpx(nz,nfp1,nl), ehpy(nz,nfp1,nl), ehpz(nz,nfp1,nl)
    real :: ephx(nz,nf,nlp1), ephy(nz,nf,nlp1), ephz(nz,nf,nlp1)

    ! Current densities
    real :: jp(nz,nf,nl), jphi(nz,nf,nl)

    ! ---------------------------------------------------------------
    ! Large work arrays formerly on the stack in subroutine exb.
    ! Moved here to prevent stack overflow for large grids.
    ! ---------------------------------------------------------------

    ! Conserved particle number and temperature * volume
    real :: denic(nz,nf,nl,nion)
    real :: tic(nz,nf,nl,nion)
    real :: tec(nz,nf,nl)

    ! Upwind fluxes — p direction (nfp1 interfaces)
    real :: fluxnp(nz,nfp1,nl,nion), fluxtp(nz,nfp1,nl,nion)
    real :: fluxtep(nz,nfp1,nl)

    ! Upwind fluxes — s direction (nz interfaces)
    real :: fluxns(nz,nf,nl,nion),  fluxts(nz,nf,nl,nion)
    real :: fluxtes(nz,nf,nl)

    ! Upwind fluxes — h direction (nl interfaces)
    real :: fluxnh(nz,nf,nl,nion),  fluxth(nz,nf,nl,nion)
    real :: fluxteh(nz,nf,nl)

end module exb_mod
