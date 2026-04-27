
module parameter_apex_mod

  use parameter_mod

      implicit none

! Geographic grid parameters:

      integer, parameter :: nlat = 90,nlon = 180,nlev = 48,nalt = 85
      integer, parameter :: nlonp1=nlon+1, nlonp2=nlon+2, nlatp1=nlat+1
      real, parameter :: glat1 = -87.5,dlat  = 2.,glon1 = 0.,dlon  = 2.
      real, parameter :: spval = 1.e36, ispval = 999

! Magnetic grid:
      integer, parameter ::  nmlat = 97,  nmlon = 80
      integer, parameter ::  nmlatp1=nmlat+1,nmlath=(nmlat+1)/2
      integer, parameter ::  nmlonp1=nmlon+1,nmlonp2=nmlon+2

! Magnetospheric grid:
      integer, parameter :: nmagphrlat = 31, nmagphrlon = 40
      real, parameter :: magphrlat1=71.97,magphrlat2=10.14,magphrlon1=-180.

! Arrays for APEX

      real,dimension(nz,nf,nl,3,3) :: dvec
      real,parameter :: h0 = 57. ! reference height        

end module parameter_apex_mod
