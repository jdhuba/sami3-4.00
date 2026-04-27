
!******************************************
!******************************************

!             courant

!******************************************
!******************************************

    subroutine courant

    use parameter_mod
    use variable_mod
    use namelist_mod
    use time_mod
    use exb_mod
    use grid_mod

! parallel motion

    dtnew = 1.e6
    do l = nion1,nion2
        do k = 1,nl
            do j = 1,nf
                do i = 1,nz
                    dt1 = dels(i,j,k) / amax1(1.,abs(vsi(i,j,k,l)))
                    if ( dt1 <= dtnew ) then
                        dtnew = dt1
                        i0par    = i
                        j0par    = j
                        k0par    = k
                        l0par    = l
                    endif
                enddo
            enddo
        enddo
    enddo

! perpendicular motion

    do k = 1,nl
        do j = 1,nf
            do i = 1,nz
                dts = xdels(i,j,k) / amax1(1.,abs(vexbs(i,j,k)))
                dtp = xdelp(i,j,k) / amax1(1.,abs(vexbp(i,j,k)))
                dth = xdelh(i,j,k) / amax1(1.,abs(vexbh(i,j,k)))
                dt1 = amin1 ( dts,dtp,dth )
                if ( dt1 <= dtnew ) then
                    dtnew = dt1
                    i0per    = i
                    j0per    = j
                    k0per    = k
                endif
            enddo
        enddo
    enddo

    if ( dtnew <= .01 ) then
        print *,' Time step too small',dtnew,i0par,j0par,k0par,l0par
        print *,' Time step too small',dtnew,i0par,j0per,k0per,l0per
        print *,' vparallel',vsi(i0par,j0par,k0par,l0par)
        stop
    elseif ( dtnew >= 5e4 ) then
        print *,' Time step too big: dtnew',dtnew
        stop
    endif
    dtnew = 1.5 * dtnew
    if ( dtnew <= dt      ) dt = dtnew
    if ( dtnew > dt*1.2  ) dt = dt * 1.2

    return
    end subroutine courant


! *********************

!     smoothz

! *********************

    subroutine smoothz(finout,ncomp)
     
    use parameter_mod
     
    dimension finout(nz), tempz(nz)
     

! This is the binomial filter (in x space) as described in
! Birdsall appendix C.
! We have the choice of a compensating filter or not.
! if ncomp=0, no compensation, else compensation

     
! do smoothz in the z direction
     
    do i = 1,nz
        ip1 = i +1
        if(i == nz) ip1 = 1
        im1 = i -1
        if(i == 1) im1 = nz
        tempz(i) = .25*(finout(im1) +2.*finout(i) &
        +finout(ip1))
    enddo
    do i = 1,nz
        finout(i) = tempz(i)
    enddo
     
    if ( ncomp /= 0 ) then
         
    ! put in compensator
    ! the alogrithm below is equivalent to
    ! fftmp(i)=(1./16.)*(-ff0(i-2)+4.*ff0(i-1)+10.*ff0(i)+4.*ff0(i+1)-ff0(i+2))
         
    ! do compensation in the z direction
         
        const = sqrt(1.4571072)
        do i = 1,nz
            ip1 = i +1
            if(i == nz) ip1 = 1
            finout(i) = const*(finout(i) -.171573*finout(ip1))
        enddo
        do i = nz,1,-1
            im1 = i -1
            if(i == 1) im1 = nz
            finout(i) = const*(finout(i) -.171573*finout(im1))
        enddo
         
    endif
     
    return
    end subroutine smoothz



!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************

    subroutine smoothx(f)
     
    use parameter_mod

    parameter ( nnxp2 = nnx + 2, nnyp2 = nny + 2 )
    parameter ( nnxp1 = nnx + 1, nnyp1 = nny + 1 )
             
    real :: f(nnx,nny),f0(nnx+2,nny+2)
     
    u12 = 1.
     
    do j = 1,nny
        do i = 1,nnx
            f0(i+1,j+1) = f(i,j)
        enddo
    enddo
     
!       zero-gradient in x
     
!        do j = 2,nnyp1
!          f0(1,j)     = f0(2,j)
!          f0(nnx+2,j)  = f0(nnx+1,j)
!        enddo

!       periodic in x
     
    do j = 2,nnyp1
        f0(1,j)     = f0(nnx,j)
        f0(nnx+2,j)  = f0(3,j)
    enddo
     
!       zero gradient in y
     
    do i = 1,nnxp2
        f0(i,1)    = f0(i,2)
        f0(i,nnyp2) = f0(i,nnyp1)
    enddo
     
!       sweep in x (1/2)
        
    do j = 2,nnyp1
        do i = 1,nnxp1
            f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i+1,j) )
        enddo
    enddo
     
    do j = 2,nnyp1
        do i = nnxp2,2,-1
            f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i-1,j) )
        enddo
    enddo
     
!       now get f
     
    do j = 1,nny
        do i = 1,nnx
            f(i,j)  = f0(i+1,j+1)
        enddo
    enddo
     
    return
    end subroutine smoothx

!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************
!*******************************************************


    subroutine smoothy(f)
     
    use parameter_mod

    parameter ( nnxp2 = nnx + 2, nnyp2 = nny + 2 )
    parameter ( nnxp1 = nnx + 1, nnyp1 = nny + 1 )
             
    real :: f(nnx,nny),f0(nnx+2,nny+2)
     
    u12 = 1.
     
    do j = 1,nny
        do i = 1,nnx
            f0(i+1,j+1) = f(i,j)
        enddo
    enddo
            
!       zero-gradient in x
     
!        do j = 2,nnyp1
!          f0(1,j)     = f0(2,j)
!          f0(nnx+2,j)  = f0(nnx+1,j)
!        enddo
     
!       periodic in x
     
    do j = 2,nnyp1
        f0(1,j)     = f0(nnx,j)
        f0(nnx+2,j)  = f0(3,j)
    enddo

!       zero gradient in y
     
    do i = 1,nnxp2
        f0(i,1)    = f0(i,2)
        f0(i,nnyp2) = f0(i,nnyp1)
    enddo
     
!       sweep in y (1/2)
     
    do j = 1,nnyp1
        do i = 1,nnxp2
            f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i,j+1) )
        enddo
    enddo
     
    do j = nnyp2,2,-1
        do i = 1,nnxp2
            f0 (i,j) = 0.5 * ( f0(i,j) + u12*f0(i,j-1) )
        enddo
    enddo
     
!       now get f
     
    do j = 1,nny
        do i = 1,nnx
            f(i,j)  = f0(i+1,j+1)
        enddo
    enddo
     
    return
    end subroutine smoothy
     



!******************************************
!******************************************
!          splinenr
!    (from numerical recipes)
!******************************************
!******************************************

    subroutine splinenr(x,y,n,yp1,ypn,y2)
    parameter (nmax=200)
    dimension x(n),y(n),y2(n),u(nmax)
    if (yp1 > .99e30) then
        y2(1)=0.
        u(1)=0.
    else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    11 END DO
    if (ypn > .99e30) then
        qn=0.
        un=0.
    else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
    12 END DO
    return
    end subroutine splinenr

!******************************************
!******************************************
!          splintnr
!    (from numerical recipes)
!******************************************
!******************************************

    subroutine splintnr(xa,ya,y2a,n,x,y)
    dimension xa(n),ya(n),y2a(n)
    klo=1
    khi=n
    1 if (khi-klo > 1) then
        k=(khi+klo)/2
        if(xa(k) > x)then
            khi=k
        else
            klo=k
        endif
        goto 1
    endif
    h=xa(khi)-xa(klo)
    if (h == 0.) print *, 'bad xa input. from splintnr'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+ &
    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
    return
    end subroutine splintnr



!*******************************************************

    subroutine smoothz_1D(f)
     
    use parameter_mod

    real :: f(nz),f0(nz+2)
     
    do i = 1,nz
      f0(i+1) = f(i)
    enddo
     
!   zero-gradient in z
     
    f0(1)     = f0(2)
    f0(nz+2)  = f0(nz+1)

!   sweep in z (forward)
        
    do i = 1,nz+1
      f0 (i) = 0.5 * ( f0(i) + f0(i+1) )
    enddo

!   sweep in z (backward)
     
    do i = nz+2,2,-1
      f0 (i) = 0.5 * ( f0(i) + f0(i-1) )
    enddo
     
!       now get f
     
    do i = 1,nz
      f(i)  = f0(i+1)
    enddo
     
    return
    end subroutine smoothz_1D


!    pure real(8) function kcc1(nx,ny,x,y,f,xp,yp)
    real(8) function kcc1(nx,ny,x,y,f,xp,yp)      

    implicit none

! Input variables

    integer,intent(in)     :: nx
    integer,intent(in)     :: ny
    real(8),intent(in)     :: x(nx)
    real(8),intent(in)     :: y(ny)
    real(8),intent(in)     :: f(nx,ny)

    real(8),intent(in)     :: xp
    real(8),intent(in)     :: yp

! Internal variables

    real(8)     :: x1,x2,x3,x4
    real(8)     :: y1,y2,y3,y4
    integer     :: i1,i2,i3,i4
    integer     :: j1,j2,j3,j4
    integer     :: i,j

! Apply the Key's Cubic Convolution interpolation algorithm

! Calculate the interpolation kernel for the x direction

    i1 = int(xp) + 1
    i2 = i1 - 1
    i3 = i1 + 1
    i4 = i1 + 2
    x1 = abs(xp - x(i1))
    x2 = x1 + 1.0d0
    x3 = 1.0d0 - x1
    x4 = 2.0d0 - x1
    x1 =  1.5d0*x1*x1*x1 - 2.5d0*x1*x1 + 1.0d0
    x2 = -0.5d0*x2*x2*x2 + 2.5d0*x2*x2 - 4.0d0*x2 + 2.0d0
    x3 =  1.5d0*x3*x3*x3 - 2.5d0*x3*x3 + 1.0d0
    x4 = -0.5d0*x4*x4*x4 + 2.5d0*x4*x4 - 4.0d0*x4 + 2.0d0

    ! Calculate the interpolation kernel for the y direction

    j1 = int(yp) + 1
    j2 = j1 - 1
    j3 = j1 + 1
    j4 = j1 + 2
    y1 = abs(yp - y(j1))
    y2 = y1 + 1.0d0
    y3 = 1.0d0 - y1
    y4 = 2.0d0 - y1
    y1 =  1.5d0*y1*y1*y1 - 2.5d0*y1*y1 + 1.0d0
    y2 = -0.5d0*y2*y2*y2 + 2.5d0*y2*y2 - 4.0d0*y2 + 2.0d0
    y3 =  1.5d0*y3*y3*y3 - 2.5d0*y3*y3 + 1.0d0
    y4 = -0.5d0*y4*y4*y4 + 2.5d0*y4*y4 - 4.0d0*y4 + 2.0d0

        ! Sum the contributions from the eight adjacent data points

    kcc1 = y1*(f(i1,j1)*x1 + f(i2,j1)*x2 + f(i3,j1)*x3 + f(i4,j1)*x4) + &
           y2*(f(i1,j2)*x1 + f(i2,j2)*x2 + f(i3,j2)*x3 + f(i4,j2)*x4) + &
           y3*(f(i1,j3)*x1 + f(i2,j3)*x2 + f(i3,j3)*x3 + f(i4,j3)*x4) + &
           y4*(f(i1,j4)*x1 + f(i2,j4)*x2 + f(i3,j4)*x3 + f(i4,j4)*x4)

    return

end function kcc1
