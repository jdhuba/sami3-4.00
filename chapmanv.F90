! Chapman integral solver
! Written by R. H. Varney for ionospheric models in the Center for Geospace Storms
!
! This code provides two different functions for evaluating the Chapman function
!
! chapmanInt(X,chi) is developed to be highly accurate (~1e-14 error) for all values of input parameters
! This function uses numerical integration, and is relatively slow.
!
! chapmanApprox(X,chi) uses various approximations to rapidly accelerate the computation. For X>36, this
! option is equivalent to the large-X approximation  in Huestis (2001) section 3, but I have written a
! new implementation that accelerates the computation and only requires a single call to erfcx.
! For X< 36 this falls back to using numerical integration, but it can do so with a small number of points
! since the integrals converge more rapidly for smaller X.
!
!
! The chapman grazing incidence integral is related to the atmospheric column depth
! assuming a spherically symmetric and exponential atmosphere
!
! For theory read
!
! David L. Huestis,
! Accurate evaluation of the Chapman function for atmospheric attenuation,
! Journal of Quantitative Spectroscopy and Radiative Transfer,
! Volume 69, Issue 6,
! 2001,
! Pages 709-721,
! ISSN 0022-4073,
! https://doi.org/10.1016/S0022-4073(00)00107-2
!
!
! Ch(X,chi) = 1/(H * n(z0) ) N(chi,z0)
! where
! H is the atmospheric scale height
! chi is the solar zenith angle
! and
! N(chi,z0) is the column total integral of the atmospheric density from altitude z0 upwards towards the sun
!
! For chi <= pi/2, the Chapman Integral is (Huestis equation 1):
!
! Ch(X,chi) = 1/H * integral_z0 ^infinity exp(-(z-z0)/H) * sqrt(1-((R_E+z0)/(R_E+z))*sin^2(chi)) dz
! where X = (R_E + z) / H
!
! For chi > pi/2 (see Banks and Kockarts equations 17.19 and 17.20)
!
! Ch(X,chi) = 2*Chi(Y,pi/2)*exp(X-Y) - Ch(X,pi-chi)
!
! Where Y = X*sin(pi-chi) is the normalized distance of the atmospheric tangent point
!
! Approach for accurate chapmanInt routine:
!
! The numerical approach for this routine begins with a transformed version of the integral
!
! Huestis equation 5 (only valid for chi <= pi/2)
!
! Ch(X,chi) = 1 + X * sin(chi) * integral_0 ^chi exp(X(1-sin(chi)/sin(lambda)))/(1+cos(lambda)) dlambda
!
! For chi > pi/2, the above identity is used instead.
!
! This integrand is skewed towards the right side of the interval.
! Near lambda = 0 this integrand approaches exp(-infinity)
! Near lambda = chi this integrand approaches exp(0)/(1+cos(chi))
! Therefore we want more points in the numerical integration clustered on the right side of the interval
!
! The interval from [0,chi] is divided into two segments from [0,2/3*chi] and [2/3*chi,chi]
! Then the second segment is divided into two segments from [2/3*chi,2/3*(2/3*chi)] and [2/3*(2/3*chi),chi]
! and so forth for nsegs levels of refinement.
!
! Each segment is then evaluated with Gauss-Legendre quadrature.
!
! The parameters nsegs = 4 and nGaussPoints = 16 were chosen after trial-and-error
! These parameters consistently achieve absolute errors < 5e-14 over the entire range
! X in [20,2000]
! chi in [0,pi]
! compared to using scipy.integrate.quad with a very small tolerance
!
!
! Approach for chapmanApprox:
!
! This approach is based on the large-X approximation in Huestis (2001) section 3.
!
! The chapman integral is written as (Huestis Eq. 6)
! 
! Ch(X,chi) = 2 X integral_t0 ^ infinity exp[-X(t^2 - t0^2)] f(chi) dt
!
! With
! t0^2 = 1 - sin(chi)
! a = 1 + sin(chi)
! w = (t^2 - t0^2)/a
! f(chi) = 1/sqrt(1+sin(chi)) * (1+a*w)/sqrt(1+w) 
!
! The function f(chi) is taylor expanded around w=0 out to 4 terms
!
! f(chi) = 1/sqrt(1+sin(chi)) * (1 + (a-1/2)*w + (1/8)*(-4a+3)*w^2 + (1/16)*(6a-5)*w^3)
!
! This series only converges for |w|<1 since the function has a non-analytic branch point at w=-1
! Nonetheless, for large X the assumption is that the rapidly decaying exp(-X(t^2-t0^2)) will squash
! the divergent part of the series for large w. This approach cannot work for small X, no matter
! how many extra terms are included in the series, since the series is divergent.
!
! With this series substitution, the integral can be expanded into a series of upper incomplete Gamma
! functions with half-integer arguments.
!
! Gamma(s,x) = integral_x ^ infinity t^(s-1) exp(-t) dt
!
! For s = 1/2
!
! Gamma(1/2,x) = sqrt(pi) * erfc(sqrt(x))
!
! For larger s, integration by parts yields a recurrence relation that can be applied repeatedly
! to walk down to s=1/2.
!
! Gamma(n+1/2,x) = x^(n-1+1/2)*exp(-x) + (n-1+1/2)*Gamma(n-1+1/2,x)
!
! I have expanded the entire series and collected terms to generate the simplest possible expression.
! The final algorithm only requires one single call to erfc_scaled, which is the special function for
!
! erfcx(x) = exp(x^2)*erfc(x)
!
! At X=36, the Huestis approximation has an absolute error around ~5e-7
! It gets a lot better as X gets larger, and a lot worse as X gets smaller
!
! For X<36, I fall back to using 16-point Gauss-Legendre quadrature
! At X=36, this quadrature has absolute error around ~1e-8,
! and it gets better as X gets smaller.
! Therefore the more elaborate segement integration is not needed.
!
!
! Special case of chi = 90 degrees:
!
! Ch(X,chi=pi/2) = X exp(X) * K1(X)
! Where K1 is the first-order modifed Bessel function of the second kind.
!
! For moderately large X, this Bessel function has an asymptotic form
! (Abramowitz and Stegun 9.7.2)
!
! K1(X) = sqrt(pi/(2*X)) exp(-X) [1 + (4-1)/(8X) + (4-1)(4-9)/(8X)^2 + (4-1)(4-9)(4-25)/(8X)^3 + ...]
!
! This means the chapman function at 90 degrees has a simple series representation
!
! Chi(X,chi=pi/2) = sqrt(pi/2) * sqrt(X) * [1 + (4-1)/(8X) + (4-1)(4-9)/(8X)^2 + (4-1)(4-9)(4-25)/(8X)^3 + ...]
!
! The function chap90approx(X) evaluates this series
!
! The code utilizes this function to evalute the identity for chi>pi/2 more rapidly.
!
!
! TL;DR
! Call chapmanApprox(X,chi) if you want speed
! Call chapmanInt(X,chi) if you want accuracy
! Remember that chi is in radians

module chapmanv
  !use params
  implicit none
  integer, parameter :: rrp = selected_real_kind(15, 307)
  
contains
  
  !Chapman integral top level function for fast approximations
  !X = (R_E + z) / H (dimensionless number)
  !chi0 = solar zenith angle in radians
  real(rrp) function chapmanApprox(x,chi0)
    implicit none
    real(rrp), intent(in) :: x
    real(rrp), intent(in) :: chi0
    real(rrp), parameter :: piover2 = 1.5707963267948966_rrp
    real(rrp), parameter :: pi = 3.141592653589793_rrp

    real(rrp) y

    real(rrp) integral, sinchi !only needed if the approximation needs to fall back to using numerical integration for small X

    if (chi0 .eq. 0.0_rrp) then
       ! Identity by definition for chi=0.0
       chapmanApprox = 1.0_rrp
    elseif (chi0 .eq. piover2) then
       ! Use the fast approximation for 90 degrees
       chapmanApprox = chap90approx(x)
    elseif (chi0 .lt. piover2) then
       if (x.lt.36.0_rrp) then
          ! fall back to 16-point quadrature for the small X case
          sinchi = sin(chi0)
          integral = gaussQ16(0.0_rrp,chi0,x,sinchi)
          chapmanApprox = 1.0_rrp + x*sinchi*integral
       else
          ! For large X, use the fast approximation
          chapmanApprox = chap1approx(x,chi0)
       endif
    else
       ! for solar zenith angles above 90 degrees
       ! apply the identity
       ! Ch(X,chi) = 2*Chi(Y,pi/2)*exp(X-Y) - Ch(X,pi-chi)
       ! Where Y = X*sin(pi-chi) is the normalized radius of the atmospheric tangent point

       !Use the fast approximation for 90 degrees for the first term
       y = x*sin(pi-chi0)
       chapmanApprox = 2*chap90approx(y)*exp(x-y)

       !second term
       if (x.lt.36.0_rrp) then
          ! fall back to 16-point quadrature for the small X case
          sinchi = sin(pi-chi0)
          integral = gaussQ16(0.0_rrp,pi-chi0,x,sinchi)
          chapmanApprox = chapmanApprox - 1.0_rrp - x*sinchi*integral
       else
          ! For large X, use the fast approximation
          chapmanApprox = chapmanApprox - chap1approx(x,pi-chi0)
       endif
    endif
    return

  end function chapmanApprox

  ! Fast implementation of the approximation from Huestis (2001) section 3
  ! This implementation is designed to only call erfcx once, sin once, sqrt 3 times,
  ! and otherwise only use elementary operations
  real(rrp) function chap1approx(x,chi)
    implicit none
    real(rrp), intent(in) :: x
    real(rrp), intent(in) :: chi

    real(rrp), parameter :: rrpi = 1.7724538509055159_rrp !sqrt(pi)                                                                                
    real(rrp) :: sinchi,c0,c1,c2,c3,t0sq,Rsq,R,F,q0,q1,q2,q3,xi

    sinchi = sin(chi)

    ! These coefficients come from the Taylor expansion of f(chi)
    ! Have the same meaning as C0 - C3 in Huestis (2001)
    C0 = 1/sqrt(1+sinchi)
    C1 = 0.5*(1+2*sinchi)*C0**3
    C2 = -0.125*(1+4*sinchi)*C0**5
    C3 = 0.0625*(1+6*sinchi)*C0**7

    ! t0^2 in Huestis (2001)
    t0sq = 1 - sinchi

    ! Convenient parameters to evaluate only once
    Rsq = X*t0sq
    R = sqrt(Rsq)
    F = rrpi*erfc_scaled(R)

    ! Coefficients of a polynomial expansion in 1/X
    ! These do not appear in Huestis' paper,
    ! but they can be derived after considerable algebra
    ! and repeated application of Gamma-function recurrences.
    q0 = (C0 + t0sq*(-C1 + t0sq*(C2-C3*t0sq))) * F

    q1 = (C1 + t0sq*(-2*C2 + 3*C3*t0sq)) * &
        (R + 0.5*F)

    q2 = (C2 - 3*C3*t0sq) * &
        ((Rsq + 1.5)*R + 0.75*F)

    q3 = (C3) * &
        (((Rsq + 2.5)*Rsq + 3.75)*R + 1.875*F)

    ! Final result in terms of a polynomial expansion in 1/X
    xi = 1.0_rrp/X
    chap1approx = sqrt(X)*(q0 + xi*(q1 + xi*(q2 + xi*q3)))

    return
  end function chap1approx

  ! This is a fast approximation of
  ! Ch(X,chi=pi/2) = X exp(X) K1(X)
  ! Using the asymptotic series expansion of the Bessel function
  real(rrp) function chap90approx(x)
    implicit none
    real(rrp), intent(in) :: x

    real(rrp) xi

    real(rrp), parameter :: rrpiover2 = 1.2533141373155001_rrp !sqrt(pi/2)                                                                        
    xi=1.0_rrp/x

    chap90approx = rrpiover2*sqrt(x)*(1+xi*(0.375+xi*(-0.1171875+xi*0.1025390625)))

    return
  end function chap90approx
  
  !Chapman integral top level function for highly accurate numerical integration
  !X = (R_E + z) / H (dimensionless number)
  !chi = solar zenith angle in radians
  real(rrp) function chapmanInt(x,chi0)
    implicit none

    real(rrp), intent(in) :: x
    real(rrp), intent(in) :: chi0
    real(rrp), parameter :: piover2 = 1.5707963267948966_rrp
    real(rrp), parameter :: pi = 3.141592653589793_rrp

    real(rrp) y !only needed for chi > piover2
    
    if (chi0 .lt. piover2) then
       chapmanInt = chap1(x,chi0)
    else
       y = x*sin(pi-chi0)
       chapmanInt = 2*chap1(y,piover2)*exp(x-y)
       chapmanInt = chapmanInt - chap1(x,pi-chi0)
    endif
    return
  end function chapmanInt

  !Evaluate the Ch(X,chi) integral for chi< pi/2 using Huestis equation 5
  real(rrp) function chap1(x,chi)
    implicit none

    real(rrp), intent(in) :: x
    real(rrp), intent(in) :: chi
    real(rrp) integral
    
    integral = segIntegrator(x,chi)
    
    chap1 = 1.0 + x*sin(chi)*integral
    return
  end function chap1

  !Procedure to break the integral into multiple segments using the 2/3 1/3 division
  !This particular trick happens to work well for this integrand since it is skewed to the right
  real(rrp) function segIntegrator(x,chi)
    implicit none

    real(rrp), intent(in) :: x
    real(rrp), intent(in) :: chi
    
    integer, parameter :: nsegs=4
    integer n, m, k
    real(rrp) f1,f2,sinchi

    sinchi = sin(chi)
    segIntegrator = 0
    do n=0,nsegs-1
       m=3**n
       k=3**(n+1)
       
       f1=real(m-1,rrp)/real(m,rrp)
       f2=real(k-1,rrp)/real(k,rrp)
       segIntegrator = segIntegrator + gaussQ16(f1*chi,f2*chi,x,sinchi)
    enddo
    !last bit
    m=3**nsegs
    f1=real(m-1,rrp)/real(m,rrp)
    segIntegrator = segIntegrator + gaussQ16(f1*chi,chi,x,sinchi)
    return
  end function segIntegrator

  !integrand of the integral in Huestis equation 5
  real(rrp) function integrand(lam,x,sinchi)
    implicit none

    real(rrp), intent(in) :: lam
    real(rrp), intent(in) :: x
    real(rrp), intent(in) :: sinchi ! pass sinchi instead of chi to save the work of repeated calls to sin

    if (lam.eq.0.0) then
       integrand = 0.0_rrp
    else
       integrand = exp(x*(1-sinchi/sin(lam)))/(1+cos(lam))
    endif
    !should give integrand = 0 for lam=0
    
    return
  end function integrand

  !16 point Gauss-Legendre quadrature
  real(rrp) function gaussQ16(a,b,x,sinchi)
    implicit none

    real(rrp), intent(in) :: a,b
    real(rrp), intent(in) :: x
    real(rrp), intent(in) :: sinchi ! pass sinchi instead of chi to save the work of repated calls to sin
    
    real(rrp) span, mid
    
    integer, parameter :: n=16
    
    real(rrp), dimension(n) :: wi,xi
    
    integer i
    
    !weights
    data wi/0.1894506104550685, &
         0.1894506104550685, & 
         0.1826034150449236, & 
         0.1826034150449236, & 
         0.1691565193950025, & 
         0.1691565193950025, & 
         0.1495959888165767, & 
         0.1495959888165767, & 
         0.1246289712555339, & 
         0.1246289712555339, & 
         0.0951585116824928, & 
         0.0951585116824928, & 
         0.0622535239386479, & 
         0.0622535239386479, & 
         0.0271524594117541, & 
         0.0271524594117541/
    
    !abcissae
    data xi/-0.0950125098376374, &
         0.0950125098376374, &
         -0.2816035507792589, &
         0.2816035507792589, &
         -0.4580167776572274, &
         0.4580167776572274, &
         -0.6178762444026438, &
         0.6178762444026438, &
         -0.7554044083550030, &
         0.7554044083550030, &
         -0.8656312023878318, &
         0.8656312023878318, &
         -0.9445750230732326, &
         0.9445750230732326, &
         -0.9894009349916499, &
         0.9894009349916499/
    
    
    span = 0.5*(b-a)
    mid = 0.5*(a+b)
    
    gaussQ16 = 0.0_rrp
    
    do i=1,n
       gaussQ16 = gaussQ16 + wi(i)*integrand(span*xi(i)+mid,x,sinchi)
    enddo
    
    gaussQ16 = span*gaussQ16
    
    return
  end function gaussQ16
end module chapmanv
