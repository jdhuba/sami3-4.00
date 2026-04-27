module hardy_mod

    use parameter_mod

    implicit none

    ! ---------------------------------------------------------------
    ! Atmospheric profile table dimensions
    ! ---------------------------------------------------------------
    integer, parameter :: no_lamzor = 250
    integer, parameter :: no_atm_ps = 250

    ! ---------------------------------------------------------------
    ! Chebyshev expansion dimensions (from EAVEKP)
    ! ---------------------------------------------------------------
    integer, parameter :: ll_eave = 10
    integer, parameter :: mm_eave = 6

    ! ---------------------------------------------------------------
    ! Atmospheric profile lookup tables (read once from file)
    ! ---------------------------------------------------------------
    real :: lambda(no_lamzor)
    real :: zor(no_lamzor)
    real :: height(no_atm_ps)
    real :: z(no_atm_ps), nNn(no_atm_ps), nM(no_atm_ps), rho(no_atm_ps)

    ! ---------------------------------------------------------------
    ! Precipitation ionization rates and temperature perturbations
    ! ---------------------------------------------------------------
    real    :: preciprn(nz,nf,nl), preciprs(nz,nf,nl)
    real    :: tpn(nz,nf,nl),      tps(nz,nf,nl)
    integer :: iin(nf,nl),         iis(nf,nl)

    ! ---------------------------------------------------------------
    ! One-time file-read gates
    ! Replaces unreliable DATA/local pattern in each subroutine.
    ! Initialised to 1 (true); set to 0 after first read.
    ! ---------------------------------------------------------------
    integer :: iread_hardy  = 1   ! gates lambda_zor / atm_params read
    integer :: ifirst_eave  = 1   ! gates eavekp.inp read
    integer :: ifirst_elekp = 0   ! gates elekp.inp read (0 = unread)

    ! ---------------------------------------------------------------
    ! EAVEKP coefficient arrays
    ! Replaces COMMON /SAVEKPE/ A, B, C
    ! ---------------------------------------------------------------
    real :: eave_a(7, ll_eave)
    real :: eave_b(7, ll_eave, mm_eave)
    real :: eave_c(7, ll_eave, mm_eave)

    ! ---------------------------------------------------------------
    ! ELEKP coefficient arrays
    ! Replaces COMMON /SAVKPE/ CR0, CH0, CH1, CS0, CR1, CS2, PRAIN
    ! ---------------------------------------------------------------
    real :: elekp_cr0(2,7,17), elekp_ch0(2,7,17), elekp_ch1(2,7,17)
    real :: elekp_cs0(2,7,17), elekp_cr1(2,7,17), elekp_cs2(2,7,17)

    real :: prain(2,7) = reshape( [ &
        7.07292, 7.01250, 7.02709, 7.18542, 7.00417, 6.90417, 6.98750, &
        6.48542, 6.68333, 6.62083, 6.82708, 6.54583, 6.49167, 6.46458  &
        ], shape=[2,7] )

end module hardy_mod
