

!!! ===========================================================================
!!! NRLMSIS 2.0:
!!! Neutral atmosphere empirical model from the surface to lower exosphere
!!! John Emmert (john.emmert@nrl.navy.mil)
!!! Doug Drob (douglas.drob@nrl.navy.mil)
!!! ===========================================================================
!!!
!!! GTD8D: Legacy wrapper with input and output arguments used in NRLMSISE-00
!
!     PREREQUISITES:
!       Must first run MSISINIT to load parameters and set switches. The
!       MSISCALC subroutine (called by this wrapper) checks for initialization
!       and does a default initialization if necessary. This self-initialization
!       will be removed in future versions.
!
!     CALLING SEQUENCE:
!       CALL GTD8D(IYD, SEC, ALT, GLAT, GLONG, STL, F107A, F107, AP, MASS, D, T)
!
!     INPUT VARIABLES:
!       IYD    Year and day as YYDDD (day of year from 1 to 365 (or 366))
!                (Year is ignored in current model)
!       SEC    Universal time (seconds)
!       ALT    Geodetic altitude (km)
!       GLAT   Geodetic latitude (deg)
!       GLONG  Geodetic longitude (deg)
!       STL    Local solar time (Ignored in 2.0; calculated from SEC and GLONG)
!       F107A  81 day average, centered on input time, of F10.7 solar activity
!                index
!       F107   Daily F10.7 for previous day
!       AP     Geomagnetic activity index array:
!                (1) Daily Ap
!                (2) 3 hr ap index for current time
!                (3) 3 hr ap index for 3 hrs before current time
!                (4) 3 hr ap index for 6 hrs before current time
!                (5) 3 hr ap index for 9 hrs before current time
!                (6) Average of eight 3 hr ap indices from 12 to 33 hrs
!                    prior to current time
!                (7) Average of eight 3 hr ap indices from 36 to 57 hrs
!                    prior to current time
!              AP(2:7) are only used when switch_legacy(9) = -1.0 in MSISINIT
!       MASS   Mass number (Ignored in 2.0)
!
!     NOTES ON INPUT VARIABLES:
!       - If lzalt_type = .false. in the MSISINIT call, then the ALT input
!         argument is treated as geopotential height.
!       - The STL input argument is ignored in NRLMSIS 2.0. Instead, local time
!         is computed from universal time and longitude.
!       - F107 and F107A values are the 10.7 cm radio flux at the Sun-Earth
!         distance, not the radio flux at 1 AU.
!       - The MASS input argument is ignored in NRLMSIS 2.0; species to be
!         calculated are set in MSISINIT.
!
!     OUTPUT VARIABLES:
!       D(1)  He number density (cm-3)
!       D(2)  O number density (cm-3)
!       D(3)  N2 number density (cm-3)
!       D(4)  O2 number density (cm-3)
!       D(5)  Ar number density (cm-3)
!       D(6)  Total mass density (g/cm3)
!       D(7)  H number density (cm-3)
!       D(8)  N number density (cm-3)
!       D(9)  Anomalous oxygen number density (cm-3)
!       T(1)  Exospheric temperature (K)
!       T(2)  Temperature at altitude (K)
!
!     NOTES ON OUTPUT VARIABLES:
!       - Missing density values are returned as 9.999e-38
!       - Species included in mass density calculation are set in MSISINIT
!
!!! =========================================================================

!==================================================================================================
! GTD8D: Legacy wrapper
!==================================================================================================

module MSIS

    contains


    subroutine msis2init()

      use msis_constants, only     : rp,dmissing
      use msis_init, only          : msisinit
      implicit none
  
    !  real(4),intent(in)      :: sw(25)
      call msisinit() !switch_legacy=sw)
  
      return

    end subroutine msis2init


    subroutine gtd8d(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)

      use msis_constants, only     : rp,dmissing
      use msis_calc, only          : msiscalc
!      use tune,only                : f107scale  ! re doug email 10/26/20
      
      implicit none

      ! MSIS Legacy subroutine arguments
      integer, intent(in)         :: iyd
      real(4), intent(in)         :: sec
      real(4), intent(in)         :: alt
      real(4), intent(in)         :: glat
      real(4), intent(in)         :: glong
      real(4), intent(in)         :: stl
      real(4), intent(in)         :: f107a
      real(4), intent(in)         :: f107
      real(4), intent(in)         :: ap(7)
      integer, intent(in)         :: mass
      real(4), intent(out)        :: d(9), t(2)

      ! MSIS 1.97 subroutine arguments
      real(kind=rp)               :: xday
      real(kind=rp)               :: xutsec
      real(kind=rp)               :: xalt
      real(kind=rp)               :: xlat
      real(kind=rp)               :: xlon
      real(kind=rp)               :: xsfluxavg, xsflux
      real(kind=rp)               :: xap(1:7)
      real(kind=rp)               :: xtn
      real(kind=rp)               :: xdn(1:10)
      real(kind=rp)               :: xtex
      real(kind=rp)               :: xdnA(1:10)
      real(kind=rp)               :: xdnB(1:10)

      
      real(kind=rp)               :: slope
      
      ! Convert the legacy input arguments to the new interface values and precision
      xday = mod(iyd,1000) + sec / 86400.0_rp
      xutsec = sec
      xalt = alt
      xlat = glat
      xlon = glong
      xsfluxavg = f107a!*f107scale%now
      xsflux = f107!*f107scale%now
      xap = ap

      ! Call the new subroutine
      call msiscalc(xday,xutsec,xalt,xlat,xlon,xsfluxavg,xsflux,xap,xtn,xdn,tex=xtex)

      ! Convert the output arguments to the legacy format (mks to cgs, re-order species)
      t(1) = sngl(xtex)    ! Expospheric temperature
      t(2) = sngl(xtn)     ! Temperature at altitude
      
      
       ! Patch to fix Hydrogen below 76.5 km for SAMI4 D-region
      if (xalt .lt. 76.5 .and. xdn(6) .eq. dmissing) then
        call msiscalc(xday,xutsec,76.5_rp,xlat,xlon,xsfluxavg,xsflux,xap,xtn,xdnA,tex=xtex)
        call msiscalc(xday,xutsec,77.5_rp,xlat,xlon,xsfluxavg,xsflux,xap,xtn,xdnB,tex=xtex)
        slope = log(xdnB(6)/xdnA(6))
        xdn(6) = exp( log(xdnB(6)) + (xalt - 76.5)*slope )
      endif
      
      if (xdn(6) .eq. dmissing) then
          stop "MSIS Hydrogen Error"
      endif
      
      where (xdn .ne. dmissing) xdn = xdn*1d-6
      if (xdn(1) .ne. dmissing) xdn(1) = xdn(1)*1e3_rp
       
      d(1) = sngl(xdn(5))  ! [He]
      d(2) = sngl(xdn(4))  ! [O]
      d(3) = sngl(xdn(2))  ! [N2]
      d(4) = sngl(xdn(3))  ! [O2]
      d(5) = sngl(xdn(7))  ! [Ar]
      d(6) = sngl(xdn(1))  ! Mass density
      d(7) = sngl(xdn(6))  ! [H]
      d(8) = sngl(xdn(8))  ! [N]
      d(9) = sngl(xdn(9))  ! [Anomalous O]

      return

    end subroutine gtd8d


    !!! ===============================================================
    !!! Calculate the geometric altitude, temperature, and individual
    !!! species number density for a given pressure with MSIS 2.0
    !!! ==============================================================

    subroutine ghp8(ccyyddd,utsec,z0,glat,glon,f107a,f107,ap,press,alt,dn,tn)

          use msis_constants, only: kB,Na,g0,rp
          use msis_calc, only          : msiscalc
  
          implicit none
          
          integer,intent(in)          :: ccyyddd
          real(kind=rp),intent(in)    :: utsec
          real(kind=rp),intent(in)    :: z0       !! first guess
          real(kind=rp),intent(in)    :: glat,glon
          real(kind=rp),intent(in)    :: f107a,f107
          real(kind=rp),intent(in)    :: ap(7)
          real(kind=rp),intent(in)    :: press   !!! pressure in hPa

          real(kind=rp),intent(out)   :: alt
          real(kind=rp),intent(out)   :: dn(10)   !!! # density not different order then MSIS00
          real(kind=rp),intent(out)   :: tn


          integer                     :: iday
          real(kind=rp)               :: day,plog,delta
          real(kind=rp)               :: zkm,pzkm
          real(kind=rp)               :: xn,gz,xmbar,scl

          integer                     :: n

          real(kind=rp),parameter     :: tol = 0.000043_rp
          integer,parameter           :: maxit = 30
          real(8)                     :: alt2gph
          real(8)                     :: xlat,alt0,alt1

          day = mod(ccyyddd,1000) + utsec / 86400.0_rp
          plog = log10(press*100.0_rp)
          zkm = z0
          delta = 1.0_rp
          
          ! Calculation performed in MKS except for the altitude
          ! Very simple Newtons method
          
          n = 0
          do while ( abs(delta) .ge. tol .and. n .le. maxit)

            n = n + 1

            call msiscalc(day,utsec,zkm,glat,glon,f107a,f107,ap,tn,dn)

            xn = sum(dn(2:8))
            pzkm = kB * xn * tn
            delta = plog - log10(pzkm)

            xmbar = dn(1) / xn / 1.66E-24_rp

            xlat = dble(glat)
            alt0 = dble(zkm)
            alt1 = alt0 + 1.0d0
            gz = (alt2gph(xlat,alt1) - alt2gph(xlat,alt0)) * g0

            scl = Na * kB * tn / (xmbar * gz)

            ! difference

            zkm = zkm - scl * delta / 1000.0_rp

          end do

      alt = zkm

      return
      end subroutine 
      
end module MSIS
