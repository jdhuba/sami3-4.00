
!-----------------------------------------------------------------------
module read_data
  implicit none
!
! Data read from W05scEpot.inp or W05scBpot.inp:
!
  integer,parameter :: csize=28, d1_pot=15, d2_pot=18
  integer :: ab(csize), ls(csize), ms(csize)
  real :: alschfits(d2_pot,csize), schfits(d1_pot,csize), ex_pot(2)
  integer :: maxl_pot,maxm_pot
!
! Data read from SCHAtable.inp:
!
  integer,parameter :: d1_scha=19, d2_scha=7, d3_scha=68
  real :: allnkm(d1_scha,d2_scha,d3_scha)
  integer :: maxk_scha,maxm_scha
  real :: th0s(d3_scha)
!
! Data read from W05scBndy.inp:
!
  integer,parameter :: na=6, nb=7
  real :: bndya(na),bndyb(nb),ex_bndy(2)
!
  contains 
!-----------------------------------------------------------------------
subroutine read_potential(infile)
!
! Read ascii data file W05scEpot.inp or W05scBpot.inp, written by 
!   pro write_potential (write_data.pro)
!
  implicit none
!
! Args:
  character(len=*),intent(in) :: infile 
!
! Local:
!
  character(len=16) :: fname
  integer :: i
  integer :: csize_rd,d1_rd,d2_rd
  ! fix: lu as named parameter instead of initialised local variable
  integer,parameter :: lu = 20
!
  open(lu,file=infile,status='old')
  read(lu,"(a)") fname
  read(lu,"(28i3)") ab
  read(lu,"(3i3)") csize_rd,d1_rd,d2_rd
  if (csize_rd /= csize) then 
    write(6,"('>>> read_potential: file ',a,': incompatable csize: ',&
      &'csize_rd=',i4,' csize=',i4)") fname,csize_rd,csize
    ERROR STOP 1  ! fix: stop 'string' -> ERROR STOP
  endif
  if (d1_rd /= d1_pot) then 
    write(6,"('>>> read_potential: file ',a,': incompatable d1: ',&
      &'d1_rd=',i4,' d1_pot=',i4)") fname,d1_rd,d1_pot
    ERROR STOP 1
  endif
  if (d2_rd /= d2_pot) then 
    write(6,"('>>> read_potential: file ',a,': incompatable d2: ',&
      &'d2_rd=',i4,' d2_pot=',i4)") fname,d2_rd,d2_pot
    ERROR STOP 1
  endif
  do i=1,csize
    read(lu,"(6e20.9)") alschfits(:,i)
  enddo
  read(lu,"(2f10.3)") ex_pot
  read(lu,"(28i3)") ls
  read(lu,"(2i3)") maxl_pot,maxm_pot
  read(lu,"(28i3)") ms
  do i=1,csize 
    read(lu,"(6e20.9)") schfits(:,i)
  enddo
  close(lu)
end subroutine read_potential
!-----------------------------------------------------------------------
subroutine read_schatable(infile)
!
! Read ascii data file SCHAtable.inp, written by pro write_scha
!   (write_data.pro)
!
  implicit none
!
! Args:
  character(len=*),intent(in) :: infile 
!
! Local:
!
  character(len=16) :: fname
  integer :: i,j
  integer,parameter :: lu = 20  ! fix: lu as named parameter
!
  open(lu,file=infile,status='old')
  read(lu,"(a)") fname
  read(lu,"(2i3)") maxk_scha,maxm_scha
  do i=1,d3_scha
    do j=1,d2_scha
      read(lu,"(6e20.9)") allnkm(:,j,i)
    enddo
  enddo
  read(lu,"(8f10.4)") th0s
  close(lu)
end subroutine read_schatable
!-----------------------------------------------------------------------
subroutine read_bndy(infile)
!
! Read ascii data file W05scBndy.inp, written by pro write_bndy
!   (write_data.pro)
!
  implicit none
!
! Args:
  character(len=*),intent(in) :: infile 
!
! Local:
!
  character(len=16) :: fname
  integer :: rd_na,rd_nb
  integer,parameter :: lu = 20  ! fix: lu as named parameter
!
  open(lu,file=infile,status='old')
  read(lu,"(a)") fname
  read(lu,"(2i3)") rd_na,rd_nb
  if (rd_na /= na) then 
    ! fix: error message corrected from 'read_potential' to 'read_bndy'
    write(6,"('>>> read_bndy: file ',a,': incompatable na: ',&
      &'rd_na=',i4,' na=',i4)") fname,rd_na,na
    ERROR STOP 1  ! fix: stop 'string' -> ERROR STOP
  endif
  if (rd_nb /= nb) then 
    write(6,"('>>> read_bndy: file ',a,': incompatable nb: ',&
      &'rd_nb=',i4,' nb=',i4)") fname,rd_nb,nb
    ERROR STOP 1
  endif
  read(lu,"(8e20.9)") bndya
  read(lu,"(8e20.9)") bndyb
  read(lu,"(8e20.9)") ex_bndy
  close(lu)
end subroutine read_bndy
!-----------------------------------------------------------------------
end module read_data


!-----------------------------------------------------------------------
module w05sc
  use read_data,only: csize
  implicit none
!
  real,parameter :: pi = 4.*atan(1.)  ! fix: module-level constant; eliminates local recomputation
  real :: rad2deg,deg2rad           ! set by setmodel
  real :: bndyfitr                  ! calculated by setboundary
  real :: esphc(csize),bsphc(csize) ! calculated by setmodel
  real :: tmat(3,3),ttmat(3,3)      ! from setboundary
  integer,parameter :: mxtablesize=200
  real :: plmtable(mxtablesize,csize), &
          colattable(mxtablesize)
  real :: nlms(csize)
  contains
!-----------------------------------------------------------------------
  subroutine setmodel(by,bz,tilt,swvel,swden,file_path,model)
    use read_data,only: read_potential,ex_pot,d1_pot,csize,schfits,&
      read_schatable
    implicit none
!
! Args:
!
! file_path: directory in which to find data file (must have "/" at end)
! model: must be either 'epot' or 'bpot' for electric or magnetic potential
!
    real,intent(in) :: by,bz,tilt,swvel,swden
    character(len=*),intent(in) :: file_path,model
!
! Local:
    integer :: i,j
    ! fix: pi removed -- now the module-level parameter
    real :: bt,angle,stilt,stilt2,sw,swp,swe,c0,&
      rang,cosa,sina,cos2a,sin2a
    real :: cfits(d1_pot,csize),a(d1_pot)
!
    if (trim(model) /= 'epot'.and.trim(model) /= 'bpot') then
      write(6,"('>>> model=',a)") trim(model)
      write(6,"('>>> setmodel: model must be either epot or bpot')")
      ERROR STOP 1  ! fix: stop 'string' -> ERROR STOP
    endif
!
! Read data:
    if (trim(model) == 'epot') then
      call read_potential(file_path//'W05scEpot.inp')
    else
      call read_potential(file_path//'W05scBpot.inp')
    endif
    call read_schatable(file_path//'SCHAtable.inp')
!
    ! fix: pi is now the module parameter; rad2deg/deg2rad set from it
    rad2deg = 180./pi
    deg2rad = pi/180.
!
    bt = sqrt(by**2+bz**2)
    angle = atan2(by,bz)*rad2deg
    call setboundary(angle,bt,tilt,swvel,swden,file_path)
!
    stilt = sin(tilt*deg2rad)
    stilt2 = stilt**2
    sw = bt*swvel/1000.
    swe = (1.-exp(-sw*ex_pot(2)))*sw**ex_pot(1)
    c0 = 1.
    swp = swvel**2 * swden*1.6726e-6
    rang = angle*deg2rad
    cosa = cos(rang)
    sina = sin(rang)
    cos2a = cos(2.*rang)
    sin2a = sin(2.*rang)
    if (bt < 1.) then ! remove angle dependency for IMF under 1 nT
      cosa = -1.+bt*(cosa+1.)
      cos2a = 1.+bt*(cos2a-1.)
      sina = bt*sina
      sin2a = bt*sin2a
    endif
    cfits = schfits ! schfits(d1_pot,csize) is in module read_data
    a = (/c0      , swe       , stilt      , stilt2     , swp, &
          swe*cosa, stilt*cosa, stilt2*cosa, swp*cosa,         &
          swe*sina, stilt*sina, stilt2*sina, swp*sina,         &
          swe*cos2a,swe*sin2a/)
    if (trim(model) == 'epot') then
      esphc(:) = 0.
      do j=1,csize
        do i=1,d1_pot  ! fix: int(d1_pot) -> d1_pot (already integer parameter)
          esphc(j) = esphc(j)+cfits(i,j)*a(i)
        enddo
      enddo
    else
      bsphc(:) = 0.
      do j=1,csize
        do i=1,d1_pot  ! fix: int(d1_pot) -> d1_pot
          bsphc(j) = bsphc(j)+cfits(i,j)*a(i)
        enddo
      enddo
    endif
  end subroutine setmodel
!-----------------------------------------------------------------------
  subroutine setboundary(angle,bt,tilt,swvel,swden,file_path)
    use read_data,only: read_bndy,ex_bndy,na,bndya ! na=6
    implicit none
!
! Args:
!
! file_path: directory in which to find data file (must have "/" at end)
    character(len=*),intent(in) :: file_path
    real,intent(in) :: angle,bt,tilt,swvel,swden
!
! Local:
    integer :: i
    real :: swp,xc,theta,ct,st,tilt2,cosa,btx,x(na),c(na)

!   write(6,"('Enter setboundary: angle=',f8.3,' bt=',f8.3)") angle,bt
!
! Read data:
    call read_bndy(file_path//'W05scBndy.inp')
!
! Calculate the transformation matrix to the coordinate system
! of the offset pole.
!
    xc = 4.2
    theta = xc*(deg2rad)
    ct = cos(theta)
    st = sin(theta)
!
    tmat(1,:) = (/ ct, 0., st/) 
    tmat(2,:) = (/ 0., 1., 0./) 
    tmat(3,:) = (/-st, 0., ct/)
!
    ttmat(1,:) = (/ct, 0.,-st/)
    ttmat(2,:) = (/ 0.,1., 0./)
    ttmat(3,:) = (/st, 0., ct/)
!
    swp = swden*swvel**2*1.6726e-6 ! pressure
    tilt2 = tilt**2
    cosa = cos(angle*deg2rad)
    btx = 1.-exp(-bt*ex_bndy(1))
    if (bt > 1.) then
      btx = btx*bt**ex_bndy(2)
    else
      cosa = 1.+bt*(cosa-1.) ! remove angle dependency for IMF under 1 nT
    endif
    x = (/1., cosa, btx, btx*cosa, swvel, swp/)
    c = bndya
    bndyfitr = 0.
    do i=1,na
      bndyfitr = bndyfitr+x(i)*c(i)
    enddo

!   write(6,"('setboundary: cosa=',f8.3,' btx=',f8.3)") cosa,btx
!   write(6,"('setboundary: x=',/,(6e12.4))") x
!   write(6,"('setboundary: c=',/,(6e12.4))") c
!   write(6,"('setboundary: bndyfitr=',e12.4)") bndyfitr

  end subroutine setboundary
!-----------------------------------------------------------------------
  subroutine epotval(lat,mlt,fill,epot)
    use read_data,only: maxm_pot,ms,ab
    implicit none
!
! Args:
    real,intent(in) :: lat,mlt,fill
    real,intent(out) :: epot
!
! Local:
    integer :: inside,j,m,mm,skip
    real :: z,phir,plm,colat,nlm
    real :: phim(2),cospm(2),sinpm(2)
!
! checkinputs returns inside=1 if lat is inside model boundary,
! inside=0 otherwise. Phir and colat are also returned by checkinputs.
!

    call checkinputs(lat,mlt,inside,phir,colat)
    if (inside == 0) then
      epot = fill
      return
    endif

    !!if (colat .le. 0.332906575352760) then
    !!  print *, lat, mlt
    !!endif


!
! IDL code: 
! phim=phir # replicate(1,maxm) * ((indgen(maxm)+1) ## replicate(1,n_elements(phir)))
!   where the '#' operator multiplies columns of first array by rows of second array,
!   and the '##' operator multiplies rows of first array by columns of second array.
! Here, maxm == maxm_pot == 2 (from read_data module), and phir is a scalar. The 
!   above IDL statement then becomes: phim = ([phir] # [1,1]) * ([1,2] ## [phir]) where
!   phim will be dimensioned [1,2]
!
    phim(1) = phir
    phim(2) = phir*2.
    cospm(:) = cos(phim(:))
    sinpm(:) = sin(phim(:))

    z = 0.
    skip = 0  ! fix: initialise before loop (was used uninitialised)
    do j=1,csize
      if (skip == 1) then
        skip = 0
        cycle
      endif

      m = ms(j)
      if (ab(j)==1) then
        plm = scplm(j,colat,nlm) ! scplm function is in this module
        skip = 0
        if (m == 0) then
          z = z+plm*esphc(j)
        else
          z = z+plm*(esphc(j)*cospm(m)+esphc(j+1)*sinpm(m))
          skip = 1
        endif
      endif ! ab(j)
    enddo
    epot = z 
  end subroutine epotval
!-----------------------------------------------------------------------
  subroutine mpfac(lat,mlt,fill,fac)
    use read_data,only: ls,ms,ab
    implicit none
!
! Args:
    real,intent(in) :: lat,mlt,fill
    real,intent(out) :: fac
!
! Local:
    integer :: j,m,inside,skip
    real :: phim(2),cospm(2),sinpm(2),cfactor
    ! fix: pi removed -- now the module-level parameter
    real :: re,z,phir,plm,colat,nlm
!
    re = 6371.2 + 110. ! km radius (allow default ht=110)
!
! checkinputs returns inside=1 if lat is inside model boundary,
! inside=0 otherwise. Phir and colat are also returned by checkinputs.
!
    call checkinputs(lat,mlt,inside,phir,colat)
    if (inside == 0) then
      fac = fill
      return
    endif
!
    phim(1) = phir
    phim(2) = phir*2.
    cospm(:) = cos(phim(:))
    sinpm(:) = sin(phim(:))
!
    z = 0.
    skip = 0  ! fix: initialise before loop (was used uninitialised)
    jloop: do j=1,csize
      if (skip == 1) then
        skip = 0
        cycle
      endif
      if (ls(j) >= 11) exit jloop
      m = ms(j)
      if (ab(j) == 1) then
        plm = scplm(j,colat,nlm) ! colat and nlm are returned (both reals)
        plm = plm*(nlm*(nlm+1.))
!
! bsphc was calculated in setmodel (when setmodel called with 'bpot')
        if (m==0) then
          z = z-plm*bsphc(j)
        else
          z = z-(plm*(bsphc(j)*cospm(m)+bsphc(j+1)*sinpm(m)))
          skip = 1
        endif
      endif
    enddo jloop ! j=1,csize
    cfactor = -1.e5/(4.*pi*re**2) ! convert to uA/m2; pi from module
    fac = z*cfactor
  end subroutine mpfac
!-----------------------------------------------------------------------
  real function scplm(index,colat,nlm)
    use read_data,only: ls,ms,ab
    implicit none
!
! Args:
    integer,intent(in) :: index
    real,intent(in) :: colat
    real,intent(out) :: nlm
!
! Local:
    integer :: i,j,l,m,skip
    real :: th0,out(1),colata(1)
    real :: cth(mxtablesize)
    real,save :: prevth0=1.e36
    integer,save :: tablesize=0  ! fix: initialised (was uninitialised save variable)

    scplm = 0.
    th0 = bndyfitr
    if (prevth0 /= th0) then
      tablesize = 3*nint(th0)
      if (tablesize > mxtablesize) then 
        write(6,"('>>> scplm: tablesize > mxtablesize: tablesize=',i5,&
          &' mxtablesize=',i5,' th0=',e12.4)") tablesize,mxtablesize,th0
        ERROR STOP 1  ! fix: stop 'string' -> ERROR STOP
      endif
!
      do i=1,tablesize
        ! fix: float() -> real()
        colattable(i) = real(i-1)*(th0/real(tablesize-1))
        cth(i) = cos(colattable(i)*deg2rad)
      enddo

      prevth0 = th0
      nlms = 0. ! whole array init
      skip = 0  ! fix: initialise before loop (was used uninitialised)
      do j=1,csize
        if (skip == 1) then
          skip = 0
          cycle
        endif
        l = ls(j)
        m = ms(j)

        nlms(j) = nkmlookup(l,m,th0) ! nkmlookup in this module

        call pm_n(m,nlms(j),cth,plmtable(1:tablesize,j),tablesize)

        skip = 0
        if (m /= 0 .and. ab(j) > 0) then
          plmtable(1,j+1) = plmtable(1,j)
          nlms(j+1) = nlms(j)
          skip = 1
        endif

      enddo ! j=1,csize

    endif ! prevth0
    nlm = nlms(index)
    colata(1) = colat

    call interpol_quad_idx(plmtable(1:tablesize,index),colattable(1:tablesize),&
      colata,out,index)
    scplm = out(1)

  end function scplm
!-----------------------------------------------------------------------
  subroutine pm_n(m,r,cth,plmtable,tablesize)
    implicit none
!
! Args:
    integer,intent(in) :: m,tablesize
    real,intent(in) :: r
    real,intent(in) :: cth(tablesize)
    real,intent(out) :: plmtable(tablesize)
!
! Local:
    integer :: i,k,ii
    real :: rm,rk,div,ans,xn
    real,dimension(tablesize) :: a,x,tmp,table
!
    if (m == 0) then 
      a = 1. ! whole array op
    else
      do i=1,tablesize
        a(i) = sqrt(1.-cth(i)**2)**m
      enddo
    endif
    xn = r*(r+1.)
    x(:) = (1.-cth(:))/2.

!   write(6,"('pm_n: a=',/,(6(1pe12.4))") a
!   write(6,"('pm_n: xn=',1pe12.4") xn
!   write(6,"('pm_n: x=',/,(6(1pe12.4))") x

    table = a ! whole array init
!
    k = 1
!   write(6,"(/)")
    pmn_loop: do         ! repeat-until loop in idl code
      do i=1,tablesize
        ! fix: float() -> real()
        rm = real(m)
        rk = real(k)
        a(i) = a(i)*(x(i)*((rk+rm-1.)*(rk+rm)-xn)/(rk*(rk+rm)))
        table(i) = table(i)+a(i) ! "result" in idl code
      enddo

!     write(6,"('pm_n: k=',i3,' a=',/,(6(1pe12.4)))") k,a
!     write(6,"('pm_n: k=',i3,' table=',/,(6(1pe12.4)))") k,table

      k = k+1
      do i=1,tablesize
        div = abs(table(i))
        if (div <= 1.e-6) div = 1.e-6
        tmp(i) = abs(a(i)) / div
      enddo

!     write(6,"('pm_n: k=',i3,' abs(a)=',/,(6(1pe12.4)))") k,abs(a)
!     write(6,"('pm_n: k=',i3,' abs(table)=',/,(6(1pe12.4)))") k,abs(table)
!     write(6,"('pm_n: k=',i3,' tmp=',/,(6(1pe12.4)))") k,tmp
!     write(6,"('pm_n: k=',i3,' max(tmp)=',1pe12.4)") k,maxval(tmp)

!     write(6,"('pm_n: k=',i5,' min,max table=',2(1pe12.4),' min,max tmp=',2(1pe12.4))")&
!       k,minval(table),maxval(table),minval(tmp),maxval(tmp)

      if (maxval(tmp) < 1.e-6) exit pmn_loop
    enddo pmn_loop
    ans = km_n(m,r)

!   write(6,"('pm_n: ans=',1pe12.4,' table=',/,(6(1pe12.4)))") ans,table(1:tablesize)

    plmtable(:) = table(:)*ans

!   write(6,"('pm_n returning: tablesize=',i4,' plmtable=',/,6(1pe12.4))") &
!     tablesize,plmtable

  end subroutine pm_n
!-----------------------------------------------------------------------
  real function km_n(m,rn)
    implicit none
!
! Args:
    integer,intent(in) :: m
    real,intent(in) :: rn
!
! Local:
    integer :: i,n
    real :: rm
!
    if (m == 0) then 
      km_n = 1.
      return
    endif
    
    ! fix: float() -> real()
    rm = real(m)
    km_n = sqrt(2.*exp(lngamma(rn+rm+1.)-lngamma(rn-rm+1.))) / (2.**m*factorial(m))
!   write(6,"('km_n: m=',i3,' rn=',f8.4,' km_n=',e12.4)") m,rn,km_n

  end function km_n
!-----------------------------------------------------------------------
  real function nkmlookup(k,m,th0)
    use read_data,only: d1_scha,d2_scha,d3_scha ! d1=19, d2=7, d3=68
    use read_data,only: allnkm,th0s             ! allnkm(d1,d2,d3),th0s(d3)
    use read_data,only: maxk_scha,maxm_scha
    implicit none
!
! Args:
    integer,intent(in) :: k,m
    real,intent(in) :: th0
!
! Local:
    integer :: kk,mm
    real :: th0a(1),out(1)

    if (th0 == 90.) then
      ! fix: float(k) -> real(k)
      nkmlookup = real(k)
      return
    endif
    th0a(1) = th0
    kk = k+1
    mm = m+1
    ! fix: out-of-bounds branches now return after the corrected call,
    !      preventing fall-through to the unclamped interpol_quad below
    if (kk > maxk_scha) then
      write(6,"('>>> nkmlookup: kk > maxk: kk=',i4,' maxk=',i4)") kk,maxk_scha
      call interpol_quad(allnkm(maxk_scha,mm,:),th0s,th0a,out)
      nkmlookup = out(1)
      return
    endif
    if (mm > maxm_scha) then
      write(6,"('>>> nkmlookup: mm > maxm: kk=',i4,' maxm=',i4)") kk,maxm_scha
      call interpol_quad(allnkm(kk,maxm_scha,:),th0s,th0a,out)
      nkmlookup = out(1)
      return
    endif
    if (th0 < th0s(1)) then
      write(6,"('>>> nkmlookup: th0 < th0s(1): th0=',e12.4,' th0s(1)=',e12.4)") &
        th0,th0s(1)
    endif

    call interpol_quad(allnkm(kk,mm,:),th0s,th0a,out)
    nkmlookup = out(1)

  end function nkmlookup
!-----------------------------------------------------------------------
  subroutine checkinputs(lat,mlt,inside,phir,colat)
    implicit none
!
! Args:
    real,intent(in) :: lat,mlt
    integer,intent(out) :: inside
    real,intent(out) :: phir,colat
!
! Local:
    real :: lon,tlat,tlon,radii

    lon = mlt*15.
    call dorotation(lat,lon,tlat,tlon)
    radii = 90.-tlat
    inside = 0
    if (radii <= bndyfitr) inside = 1 ! bndyfitr from setboundary
    phir = tlon*deg2rad
    colat = radii

  end subroutine checkinputs
!-----------------------------------------------------------------------
  subroutine dorotation(latin,lonin,latout,lonout)
    implicit none
!
! Args:
    real,intent(in) :: latin,lonin
    real,intent(out) :: latout,lonout
!
! Local:
    real :: latr,lonr,stc,ctc,sf,cf,a,b,pos(3)
    integer :: i

    latr = latin*deg2rad
    lonr = lonin*deg2rad
    stc = sin(latr)
    ctc = cos(latr)
    sf = sin(lonr)
    cf = cos(lonr)
    a = ctc*cf
    b = ctc*sf
!
! IDL code: Pos= TM ## [[A],[B],[STC]]
! The ## operator multiplies rows of first array by columns of second array.
! Currently, TM(3,3) = Tmat (or TTmat if "reversed" was set)
! If called w/ single lat,lon, then a,b,stc are dimensioned (1), and
!   Pos is then (1,3)
!
    do i=1,3
      pos(i) = tmat(1,i)*a + tmat(2,i)*b + tmat(3,i)*stc
    enddo
  
    latout = asin(pos(3))*rad2deg
    lonout = atan2(pos(2),pos(1))*rad2deg

!   write(6,"('dorotation: latin,lonin=',2f9.4,' latout,lonout=',2f9.4)") &
!     latin,lonin,latout,lonout
  
  end subroutine dorotation
!-----------------------------------------------------------------------
  subroutine interpol_quad(v,x,u,p)
!
! f90 translation of IDL function interpol(v,x,u,/quadratic)
!
    implicit none
!
! Args:
    real,intent(in) :: v(:),x(:),u(:)
    real,intent(out) :: p(:)
!
! Local:
    integer :: nv,nx,nu,i,ix
    real :: x0,x1,x2
!
    nv = size(v)
    nx = size(x)
    nu = size(u)
    if (nx /= nv) then
      write(6,"('>>> interpol_quad: nx /= nv: nx=',i4,' nv=',i4)") nx,nv
      p(:) = 0.
      return
    endif

    do i=1,nu
      ix = value_locate(x,u(i))
      if (ix <= 1.or.ix >= nx) then ! bug fix by btf 12/23/09
          p(i) = 0.
          cycle                       ! bug fix by btf 12/23/09
      endif
      x1 = x(ix)
      x0 = x(ix-1)
      x2 = x(ix+1)
      p(i) = v(ix-1) * (u(i)-x1) * (u(i)-x2) / ((x0-x1) * (x0-x2)) + &
             v(ix)   * (u(i)-x0) * (u(i)-x2) / ((x1-x0) * (x1-x2)) + &
             v(ix+1) * (u(i)-x0) * (u(i)-x1) / ((x2-x0) * (x2-x1))
      
    enddo
!   write(6,"('interpol_quad: nu=',i4,' p=',/,(1pe12.4)") nu,p

  end subroutine interpol_quad
!-----------------------------------------------------------------------
  integer function value_locate(vec,val)
!
! f90 translation of IDL function value_locate
! Return index i into vec for which vec(i) <= val <= vec(i+1)
! fix: original comment had inequality backwards (was vec(i) <= val >= vec(i+1))
! Input vec must be monotonically increasing
!
    implicit none
!
! Args:
    real,intent(in) :: vec(:),val
!
! Local:
    integer :: n,i
!
    value_locate = 0
    n = size(vec)
    if (val < vec(1)) return
    if (val > vec(n)) then
      value_locate = n
      return
    endif
    do i=1,n-1
      if (val >= vec(i) .and. val <= vec(i+1)) then
        value_locate = i
        return
      endif
    enddo

  end function value_locate
!-----------------------------------------------------------------------
real function lngamma(xx)
!
! This is an f90 translation from C code copied from 
! www.fizyka.umk.pl/nrbook/c6-1.pdf (numerical recipes gammln)
!
  implicit none
  real,intent(in) :: xx
  real :: x,y,tmp,ser
  real :: cof(6) = (/76.18009172947146, -86.50532032941677, 24.01409824083091,&
                     -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5/)
  integer :: j
!
  y = xx
  x = xx
  tmp = x+5.5
  tmp = tmp-(x+0.5)*log(tmp)
  ser = 1.000000000190015
  do j=1,6  ! fix: was 1,5 -- 6th Lanczos coefficient was never applied
    y = y+1
    ser = ser+cof(j)/y
  enddo
  lngamma = -tmp+log(2.5066282746310005*ser/x)
end function lngamma
!-----------------------------------------------------------------------
real function factorial(n)
  implicit none
  integer,intent(in) :: n
  integer :: m
  if (n <= 0) then
    write(6,"('>>> factorial: n must be positive: n=',i4)") n
    factorial = 0.
    return
  endif
  if (n == 1) then
    factorial = 1.
    return
  endif
  ! fix: float() -> real()
  factorial = real(n)
  do m = n-1,1,-1
    factorial = factorial * real(m)
  enddo
end function factorial
!-----------------------------------------------------------------------
  subroutine interpol_quad_idx(v,x,u,p,idx)
!
! f90 translation of IDL function interpol(v,x,u,/quadratic)
!
    implicit none
!
! Args:
    real,intent(in) :: v(:),x(:),u(:)
    real,intent(out) :: p(:)
    integer, intent(in) :: idx
!
! Local:
    integer :: nv,nx,nu,i,ix
    real :: x0,x1,x2
!
    nv = size(v)
    nx = size(x)
    nu = size(u)
    if (nx /= nv) then
      write(6,"('>>> interpol_quad: nx /= nv: nx=',i4,' nv=',i4)") nx,nv
      p(:) = 0.
      return
    endif


    do i=1,nu
      ix = value_locate(x,u(i))
      if (ix < 1.or.ix >= nx) then ! bug fix by btf 12/23/09
          p(i) = 0.
          cycle                       ! bug fix by btf 12/23/09
      endif
      !! at boundary, use simple linear interpolation
      if (ix == 1) then ! bug fix by tsaiwei 6/12/15
       x1 = x(ix)
       x2 = x(ix+1)
       p(i) =  (v(ix+1) - v(ix))/(x2-x1)*(u(i)-x1) + v(ix)
       cycle
      endif
      x1 = x(ix)
      x0 = x(ix-1)
      x2 = x(ix+1)
      p(i) = v(ix-1) * (u(i)-x1) * (u(i)-x2) / ((x0-x1) * (x0-x2)) + &
             v(ix)   * (u(i)-x0) * (u(i)-x2) / ((x1-x0) * (x1-x2)) + &
             v(ix+1) * (u(i)-x0) * (u(i)-x1) / ((x2-x0) * (x2-x1))
      
    enddo
    !!if (idx .eq. 1 .and. u(1) .lt. 0.33) print *, nu, p(1)
!   write(6,"('interpol_quad: nu=',i4,' p=',/,(1pe12.4)") nu,p

  end subroutine interpol_quad_idx
!-----------------------------------------------------------------------
end module w05sc

 subroutine test_w05sc

! Reads OMNI solar wind data and SAMI3 grid, evaluates the Weimer 2005
! electric potential and field-aligned current model on the SAMI3 grid,
! and writes results to unformatted output files for potential.F90.

  use w05sc,only: setmodel,epotval,mpfac
  use parameter_mod
  implicit none

  ! fix: named unit number parameters replace magic literals
  integer,parameter :: U_GRID   = 121
  integer,parameter :: U_DATA   = 122
  integer,parameter :: U_PHI    = 1101
  integer,parameter :: U_FAC    = 1102
  integer,parameter :: U_POLE   = 1103

  integer,parameter :: nlat = nf+1, nlon = nlt+1

  real,parameter :: fill = 1.e-6  ! points outside boundary return this value

  character(len=*),parameter :: file_path = "./"

  ! Input solar wind parameters
  real :: by,bz,tilt,swvel,swden

  ! Grid
  real :: mlat(nlat),mlon(nlon),mlat_pole,mlt_pole
  real :: mlt(nlon)

  ! Model output on SAMI3 grid
  real :: phi_weimer(nlat,nlon)
  real :: fac_weimer(nlat,nlon)
  real :: epot(nlon,nlat),fac(nlon,nlat)
  real :: epot_pole

  ! Time loop
  integer :: i,j,ntime,year,doy,itime
  real    :: hrutw

  ! Open output files
  open(unit=U_PHI,  file='phi_weimer.inp', form='unformatted')
  open(unit=U_FAC,  file='fac_weimer.inp', form='unformatted')
  open(unit=U_POLE, file='phi_pole.inp',   form='unformatted')

  ! Read SAMI3 grid coordinates
  open(U_GRID, file='weimer_grid.dat', form='formatted')
  read(U_GRID,*) mlat, mlon
  close(U_GRID)

  ! Read OMNI solar wind input file
  open(U_DATA, file='weimer_data.inp', form='formatted')
  read(U_DATA,*) ntime, year, doy

  do itime = 1,ntime

    read(U_DATA,*) hrutw, by, bz, swvel, swden, tilt

    ! Evaluate electric potential
    call setmodel(by,bz,tilt,swvel,swden,file_path,'epot')

    do i = 1,nlon
      mlt(i) = mlon(i)/15.
    enddo

    do j = 1,nlat
      do i = 1,nlon
        call epotval(mlat(j),mlt(i),fill,epot(i,j))
      enddo
    enddo

    ! Potential at magnetic pole
    mlat_pole = 90.
    mlt_pole  =  0.
    call epotval(mlat_pole,mlt_pole,fill,epot_pole)

    ! Evaluate field-aligned current from magnetic potential model
    call setmodel(by,bz,tilt,swvel,swden,file_path,'bpot')

    do j = 1,nlat
      do i = 1,nlon
        call mpfac(mlat(j),mlt(i),fill,fac(i,j))
      enddo
    enddo

    ! Transpose to SAMI3 array layout (nlat,nlon)
    do j = 1,nlon
      do i = 1,nlat
        phi_weimer(i,j) = epot(j,i)
        fac_weimer(i,j) = fac(j,i)
      enddo
    enddo

    ! Write output; convert kV to statvolt for SAMI3 (1 kV = 1e3/300 statvolt)
    write(U_PHI) hrutw
    write(U_PHI) phi_weimer * 1.e3 / 300.
    write(U_POLE) epot_pole * 1.e3 / 300.
    write(U_FAC) hrutw
    write(U_FAC) fac_weimer

  enddo

  close(U_DATA)
  close(U_PHI)
  close(U_FAC)
  close(U_POLE)

end subroutine test_w05sc


