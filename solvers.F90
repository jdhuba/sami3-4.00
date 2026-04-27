

!******************************************
!******************************************

!            densolv2

!******************************************
!******************************************

    subroutine densolv2( ni,tdeni,prod,loss,oldion,nfl,nll )

    use parameter_mod
    use variable_mod
    use namelist_mod
    use time_mod
    use grid_mod

    implicit none

    ! Arguments
    integer, intent(in)    :: ni, nfl, nll
    real,    intent(inout) :: tdeni(nz)
    real,    intent(in)    :: prod(nz), loss(nz), oldion(nz)

    ! Local variables
    real :: a(nz), b(nz), c(nz), d(nz)
    integer :: j
    real    :: ujm1, uj, ujp1, ur, ul, a0, b0, c0

! initialize

    do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
    enddo


    do j = 2,nz-1

        ujm1  = vsi(j-1,nfl,nll,ni)/bms(j-1,nfl,nll)
        uj    = vsi(j,nfl,nll,ni)  /bms(j,nfl,nll)
        ujp1  = vsi(j+1,nfl,nll,ni)/bms(j+1,nfl,nll)
        ur = .5*( uj +ujp1)
        ul = .5*( uj +ujm1)

        ! fix #4: four independent ifs replaced with if/else if chain
        ! to guarantee mutually exclusive upwind flux selection
        if (ur >= 0. .AND. ul >= 0.) then
            a0 = -ul
            b0 =  ur
            c0 =  0.
        else if (ur <= 0. .AND. ul <= 0.) then
            a0 = 0.
            b0 = -ul
            c0 = ur
        else if (ur >= 0. .AND. ul <= 0.) then
            a0 = 0.
            b0 = ur - ul
            c0 = 0.
        else
            a0 = -ul
            b0 = 0.
            c0 = ur
        end if

        a(j) =  a0 * bms(j,nfl,nll) ** 2 / d22s(j,nfl,nll)
        b(j) = 1. / dt + loss(j) + &
        b0 * bms(j,nfl,nll) ** 2 / d22s(j,nfl,nll)
        c(j) = c0 * bms(j,nfl,nll) ** 2 / d22s(j,nfl,nll)
        d(j) = oldion(j) / dt + prod(j)

    enddo

! we will assume that they are determined by the production and loss
! at both ends of the field line

!     lower bc

    a(1) = 0.
    b(1) = 1.
    c(1) = 0.
    d(1) = &
    sqrt ( tdeni(1) * prod(1) / loss(1) ) + denmin

! upper bc

    a(nz) = 0.
    b(nz) = 1.
    c(nz) = 0.
    d(nz) = &
    sqrt ( tdeni(nz) * prod(nz) / loss(nz) ) + denmin

    call rtryds ( a,b,c,d,tdeni,nz )

    end subroutine densolv2



!******************************************
!******************************************

!            vsisolv

!******************************************
!******************************************

    subroutine vsisolv ( vii,vid,viold,snuj,nfl,nll,cs )

    use parameter_mod
    use variable_mod
    use namelist_mod
    use time_mod
    use grid_mod

    implicit none

    ! Arguments
    integer, intent(in)    :: nfl, nll
    real,    intent(inout) :: vii(nz)
    real,    intent(in)    :: vid(nz), viold(nz), snuj(nz), cs(nz)

    ! Local variables
    ! fix #5: bare DIMENSION replaced with typed declaration
    real    :: a(nz), b(nz), c(nz), d(nz)
    integer :: j
    real    :: ujm1, uj, ujp1, ur, ul, a0, b0, c0
    real    :: delcs, alpha_drag_p, alpha_drag_n, anu_drag

! initialize

    do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
    enddo

    do j = 2,nz-1

        ujm1 = vii(j-1)
        uj   = vii(j)
        ujp1 = vii(j+1)
        ur = .25*( uj +ujp1)
        ul = .25*( uj +ujm1)

        ! fix #4: four independent ifs replaced with if/else if chain
        if (ur >= 0. .AND. ul >= 0.) then
            a0 = -ul
            b0 =  ur
            c0 =  0.
        else if (ur <= 0. .AND. ul <= 0.) then
            a0 = 0.
            b0 = -ul
            c0 = ur
        else if (ur >= 0. .AND. ul <= 0.) then
            a0 = 0.
            b0 = ur - ul
            c0 = 0.
        else
            a0 = -ul
            b0 = 0.
            c0 = ur
        end if

    ! anomalous drag to prevent ions from going supersonic

        delcs        = 0.1 * cs(j)
        alpha_drag_p = ( vii(j) - 0.9 * cs(j) )  / delcs
        alpha_drag_n = ( vii(j) + 0.9 * cs(j) )  / delcs
        anu_drag    = 0.5 * anu_drag0 * ( 1. + tanh(alpha_drag_p)) + &
                      0.5 * anu_drag0 * ( 1. - tanh(alpha_drag_n))

                 
        a(j) = a0 / d22s(j,nfl,nll) * bms(j,nfl,nll)
        ! fix #2: 1/dt was integer division (always 0); corrected to 1./dt
        b(j) = 1./dt + snuj(j) + b0 / d22s(j,nfl,nll) * &
        bms(j,nfl,nll) + anu_drag
        c(j) = c0 / d22s(j,nfl,nll) * bms(j,nfl,nll)
        d(j) = viold(j)/dt + vid(j)

    enddo

! we will assume that the bc's are zero
! at both ends of the field line

! lower bc

    a(1) = 0.
    b(1) = 1.
    c(1) = 0.
    d(1) = 0.

! upper bc
     
    a(nz) = 0.
    b(nz) = 1.
    c(nz) = 0.
    d(nz) = 0.

    call rtryds(a,b,c,d,vii,nz)

    end subroutine vsisolv

!******************************************
!******************************************

!            tisolv

!******************************************
!******************************************

    subroutine tisolv(tti,tio,kap,s1,s2,s3,s4,s5,s6,s7,npt,nfl,nll)

    use parameter_mod
    use variable_mod
    use time_mod
    use grid_mod

    implicit none

    ! Arguments
    integer, intent(in)    :: npt, nfl, nll
    real,    intent(inout) :: tti(nz)
    real,    intent(in)    :: tio(nz), kap(nz)
    real,    intent(in)    :: s1(nz), s2(nz), s3(nz), s4(nz)
    real,    intent(in)    :: s5(nz), s6(nz), s7(nz)

    ! Local variables
    real    :: a(nz), b(nz), c(nz), d(nz)
    integer :: j
    real    :: ujm1, uj, ujp1, ur, ul, a0, b0, c0

! initialize

    do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
    enddo


    do j = 2,nz-1
        ujm1 = bms(j-1,nfl,nll)*vsi(j-1,nfl,nll,npt)
        uj   = bms(j,nfl,nll)  *vsi(j,nfl,nll,npt)
        ujp1 = bms(j+1,nfl,nll)*vsi(j+1,nfl,nll,npt)
        ur = .5*( uj +ujp1)
        ul = .5*( uj +ujm1)

        ! fix #4: four independent ifs replaced with if/else if chain
        if (ur >= 0. .AND. ul >= 0.) then
            a0 = -ul
            b0 =  ur
            c0 =  0.
        else if (ur <= 0. .AND. ul <= 0.) then
            a0 = 0.
            b0 = -ul
            c0 = ur
        else if (ur >= 0. .AND. ul <= 0.) then
            a0 = 0.
            b0 = ur - ul
            c0 = 0.
        else
            a0 = -ul
            b0 = 0.
            c0 = ur
        end if

        a(j) =     a0 / d22s(j,nfl,nll) &
        - ( bms(j,nfl,nll)**2 / deni(j,nfl,nll,npt) ) / &
        d22s(j,nfl,nll) &
        *.5 * ( kap(j) + kap(j-1) ) / ds(j,nfl,nll)
        ! fix #6: .333333 (truncated) replaced with (1.0/3.0) for full precision
        b(j) = 1. / dt + b0 / d22s(j,nfl,nll) &
        -(1.0/3.0) * ( bms(j,nfl,nll) &
        * (vsi(j+1,nfl,nll,npt) - &
        vsi(j-1,nfl,nll,npt) ) &
        + 5. * vsi(j,nfl,nll,npt) &
        * (bms(j+1,nfl,nll) - &
        bms(j-1,nfl,nll) ) ) &
        / d2s(j,nfl,nll) &
        +  ( bms(j,nfl,nll)**2 / deni(j,nfl,nll,npt) ) / &
        d22s(j,nfl,nll) &
        *(.5* (kap(j+1) + kap(j) ) / ds(j+1,nfl,nll) &
        +.5 * (kap(j) + kap(j-1) ) / ds(j,nfl,nll)) &
        + s2(j) + s4(j) + s6(j) !+ anu_drag
        c(j) =     c0 / d22s(j,nfl,nll) &
        - ( bms(j,nfl,nll)**2 / deni(j,nfl,nll,npt) ) / &
        d22s(j,nfl,nll) &
        *.5 * (kap(j+1) + kap(j) ) / ds(j+1,nfl,nll)
        d(j) = tio(j)/dt + s1(j) + s3(j) + s5(j) + s7(j)

    enddo

! we will assume that the bc's are the neutral temperature
! at both ends of the field line

! lower bc

    a(1) = 0.
    b(1) = 1.
    c(1) = 0.
    d(1) = tn(1,nfl,nll)

! upper bc
     
    a(nz) = 0.
    b(nz) = 1.
    c(nz) = 0.
    d(nz) = tn(nz,nfl,nll)

    call rtryds ( a,b,c,d,tti,nz )

    end subroutine tisolv

!******************************************
!******************************************

!            tesolv

!******************************************
!******************************************

    subroutine tesolv(tte,te_old,kap,s1,s2,s3,s4,s5,nfl,nll)

    use parameter_mod
    use variable_mod
    use time_mod
    use grid_mod

    implicit none

    ! Arguments
    integer, intent(in)    :: nfl, nll
    real,    intent(inout) :: tte(nz)
    real,    intent(in)    :: te_old(nz), kap(nz)
    real,    intent(in)    :: s1(nz), s2(nz), s3(nz), s4(nz), s5(nz)

    ! Local variables
    ! fix #5: bare DIMENSION replaced with typed declarations
    real    :: a(nz), b(nz), c(nz), d(nz)
    integer :: j

! initialize

    do j = 1,nz
        a(j) = 0.
        b(j) = 0.
        c(j) = 0.
        d(j) = 0.
    enddo

! note: ne used here is from variable_mod

    do j = 2,nz-1

        a(j) = - bms(j,nfl,nll)**2 / ne(j,nfl,nll) / d22s(j,nfl,nll) &
        *.5 * ( kap(j) + kap(j-1) ) / ds(j,nfl,nll)
        b(j) = 1. / dt + bms(j,nfl,nll)**2 / ne(j,nfl,nll) / &
        d22s(j,nfl,nll) &
        *(  .5 * (kap(j+1) + kap(j)   ) /ds(j+1,nfl,nll) &
        +.5 * (kap(j)   + kap(j-1) ) /ds(j,nfl,nll)   ) &
        + s2(j) + s3(j)
        c(j) = - bms(j,nfl,nll)**2 / ne(j,nfl,nll) /d22s(j,nfl,nll) &
        *.5 * ( kap(j+1) + kap(j) )/ ds(j+1,nfl,nll)
        d(j) = te_old(j)/dt + s1(j) + s4(j) + s5(j)

    enddo

! we will assume that the bc's are the neutral temperature
! at both ends of the field line

! lower bc

    a(1) = 0.
    b(1) = 1.
    c(1) = 0.
    d(1) = tn(1,nfl,nll)

! upper bc
     
    a(nz) = 0.
    b(nz) = 1.
    c(nz) = 0.
    d(nz) = tn(nz,nfl,nll)

    call rtryds(a,b,c,d,tte,nz)

    end subroutine tesolv

!******************************************
!******************************************

!            rtryds

!******************************************
!******************************************

    subroutine rtryds(a,b,c,d,x,n)

    implicit none

    ! Arguments
    integer, intent(in)    :: n
    real,    intent(inout) :: a(n), b(n), c(n), d(n)
    real,    intent(out)   :: x(n)

    ! Local variables
    ! fix #5: bare DIMENSION replaced with typed declarations
    ! fix (issue #7): alfa/beta sized (n) not (nz) -- correct for any n
    real    :: alfa(n), beta(n)
    integer :: nm1, k, i
    real    :: ast, dst, rb, z

! arrays a,b,c, and d may be used for storage of alfa, beta and x
! in the actual call of this routine, but remember, whatever you
! use will be lost by the definition of alfa and beta here.
! fix #11: "stoage" -> "storage"; "of of" -> "of"
! form,  a(k)*x(k-1) + b(k)*x(k) + c(k)*x(k+1) = d(k)

! i have modified the input sequence to the routine, but have left it
! otherwise intact.  we may want to eventually change this (gj)

    nm1=n-1

! apply the boundary condition at x(1)
! alfa(1) and beta(1) determined from b(1),c(1),d(1),a(1)=0.

    dst     = d(1)
    rb      = 1. / b(1)
    alfa(1) = -c(1) * rb
    beta(1) =   dst * rb

! calculate the alfas and betas of k on forward sweep

    do k=2,nm1
        ast     =  a(k)
        z       =  1. / ( b(k) + ast * alfa(k-1) )
        dst     =  d(k)
        alfa(k) = -c(k) * z
        beta(k) =  ( dst - ast * beta(k-1) ) * z
    enddo

! apply the boundary condition at x(n)
! x(n) determined from a(n),b(n),d(n),c(n)=0.

    x(n) = ( d(n) - a(n) *beta(nm1) ) / ( b(n) + a(n) * alfa(nm1) )

! calculate x of k from the alfas and betas on backward sweep

    do i=2,n
        k    = n + 1 - i
        x(k) = x(k+1) * alfa(k) + beta(k)
    enddo

    end subroutine rtryds
