module misc_mod

  use parameter_mod

!     diagnostic variables

    implicit none

    real :: u1(nz,nf,nl),u2(nz,nf,nl),u3(nz,nf,nl), &
            u4(nz,nf,nl),u5(nz,nf,nl)

    real :: ppt(nzp1,nfp1,nlt)
    real :: blonp0t(nlt+3)

    real :: deni_mnp(nz,nion),ti_mnp(nz,nion),te_mnp(nz)

    real :: plat,plon

!    character(128) :: dir = './data_dir/'

end module misc_mod
