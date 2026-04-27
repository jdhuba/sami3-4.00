! WEIMER_GRID.F

    subroutine weimer_grid

    use parameter_mod
    use grid_mod          ! fix #2/#5: blatpt, blonpt come from here;
                          ! local ~62 MB stack copies removed

    implicit none         ! fix #1

    ! fix #3: named unit parameters replace magic literals 76, 77, 78
    integer, parameter :: U_BLATP      = 76
    integer, parameter :: U_BLONP      = 77
    integer, parameter :: U_WEIMER_OUT = 78

    ! fix #2/#5: local blatpt(nzp1,nfp1,nlt) and blonpt(nzp1,nfp1,nlt)
    ! declarations removed -- these are now read from grid_mod directly.
    ! The open/read/close for units 76 and 77 are also removed.

    real :: weimer_lat(nf+1), weimer_lon(nlt+1)

    ! fix #1: declare true implicit locals
    integer :: j, k

    open ( unit=U_WEIMER_OUT, file='weimer_grid.dat', form='formatted' )

    do j = 1,nf+1
        weimer_lat(j) = blatpt(nz-1,j,1)
    enddo

    do k = 1,nlt
        weimer_lon(k) = blonpt(1,1,k)
    enddo

    weimer_lon(nlt+1) = weimer_lon(1) + 360.

    write(U_WEIMER_OUT,*) weimer_lat, weimer_lon

    close(U_WEIMER_OUT)

    end subroutine weimer_grid

