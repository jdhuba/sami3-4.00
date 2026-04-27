!     namelist data

module namelist_mod

  use parameter_mod

    implicit none

    logical :: hall,restart
    logical :: lmadala,lweimer

    integer :: imonth_1,iday_1
    integer :: imonth_2,iday_2

    integer :: nion1,nion2
    integer :: maxstep,mmass,kp(8)
    integer :: nday(8),iyear

    real :: snn(nneut),fbar(8),f10p7(8),ap(8)
    real :: hrmax, dthr, hrpr, dt0, &
            rmin, altmin, &                            
            hrinit, tvn0, tvexb0, ver, veh, vw,&
            gams, gamp, alt_crit, cqe, alt_crit_avg
    real :: vexb_max, &
            decay_time, pcrit, anu_drag0, &
            blat_max, stn,denmin,blat_min,xhardy,tphi,&
            tmax,euv_fac,s5e_fac

end module namelist_mod
