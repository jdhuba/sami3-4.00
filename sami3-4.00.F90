!     *******************************************
!     *******************************************
     
!                  SAMI3_MPI-4.00

!     *******************************************
!     *******************************************

    program main

    use parameter_mod
    use parameter_apex_mod
    use variable_mod
    use namelist_mod
    use message_passing_mod
    use chemistry_mod
    use time_mod
    use atomic_mod
    use conductance_mod
    use exb_mod
    use misc_mod
    use grid_mod
    use hardy_mod
    use photo_production_mod

    implicit none  ! fix #1

!     Some local variables

    real :: denitmp(nz,nf,nl), titmp(nz,nf,nl), vsitmp(nz,nf,nl)
    real :: nuintmp(nz,nf,nl),losstmp(nz,nf,nl)
    real :: deni_mnptmp(nz,nion),ti_mnptmp(nz,nion),te_mnptmp(nz)
    real :: denntmp(nz,nf,nl)
    real :: phi(nnx,nny)
    logical :: tflag,ttflag

    character(len=80) :: charday,str,dir
    character(len=4)  :: cyear,cday
    character(len=3)  :: ccday

! allocatable total matrices

! Total

    real, dimension(:,:,:,:), allocatable :: denit,dennt
    real, dimension(:,:,:,:), allocatable :: vsit, sumvsit, nuintt, losstt
    real, dimension(:,:,:,:), allocatable :: tit
    real, dimension(:,:,:), allocatable   :: ut, vt, vpit
    real, dimension(:,:,:), allocatable   :: tet, tnt

!     height integrated pedersen/hall conductivities

    real, dimension(:,:,:), allocatable :: u1t, u2t, u3t, u4t, u5t
    real, dimension(:,:,:), allocatable :: vnqt, vnpt, vnphit
    real, dimension(:,:,:), allocatable :: jpt, jphit
    real, dimension(:,:,:), allocatable :: u1pt, u2st, u3ht
    real, dimension(:,:,:), allocatable :: sigmapict,sigmahict
    real, dimension(:,:,:), allocatable :: sigmapt,sigmaht,sigma0t
    real, dimension(:,:,:), allocatable :: gpt
!    real, dimension(:,:,:), allocatable :: gpt,xnuin0t,xnuin1t,xnuin2t

    real, dimension(:), allocatable :: dtc

! Output matrices for restart

    real :: deniout(nz,nf,nl,nion,numwork), &
            tiout(nz,nf,nl,nion,numwork), &
            vsiout(nz,nf,nl,nion,numwork)
    real :: teout(nz,nf,nl,numwork)
    real(8) :: dphi(nnx+1,nnyt)   ! fix #10: real*8 -> real(8)

    logical :: lprnt

!     Begin MPI stuff

    include 'mpif.h'
    integer :: status(MPI_STATUS_SIZE)
    integer :: left,right
    integer :: nerrorcode

    ! fix #1: previously implicit locals in program main
    integer :: ierr, numtasks
    integer :: iwrk, jwrk, nntmp, ktmp, jtmp, itmp
    integer :: k, kk, j, nn, ni, nll, nfl
    integer :: icnt10, ifintot, ifintot1, ifintot2
    integer :: ntm, ntmmax, istep, nfile, iday, iday0
    integer :: i,j0
    real    :: xxx, dtmp
    real    :: hrut, timemax, time, tneut, thardy, tpot
    real    :: code_start, code_finish, dt_old
    real    :: deni_mnp0, ti_mnp0, te_mnp0
    real    :: tprnt, timem, timep
    real    :: f74, ai, p


! ************************ initializations ***********************************
! Find out how many tasks are in this partition and what my task id is.  Then
! define the number of worker tasks and the array partition size as chunksize.
! Note:  For this example, the MP_PROCS environment variable should be set
! to an odd number...to insure even distribution of the array to numtasks-1
! worker tasks.
! *****************************************************************************

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)
    write(*,*)'taskid =',taskid

    numworkers = numtasks-1

! Check to see if the number of processors selected agrees with
! the number of divisions in params

    if(taskid == 0) then
        if(numwork /= numworkers) then
            print *, ' numworkers is ',numworkers
            print *, ' numwork (in parameter_mod_mpi) is',numwork
            print *, ' in order for the code to work correctly '
            print *, ' these two numbers must be the same '
            print *, ' Either set np = numwork +1 and rerun or '
            print *, ' change numwork and recompile '
            call mpi_abort(MPI_COMM_WORLD, nerrorcode, ierr)
            call mpi_finalize(ierr)
        endif
    endif

! Determine what is left (down) and right (up)
! Here we assume that taskid=0 is the Master and does nothing but
! deal with handling the data

    print *, ' taskid ',taskid

    if(taskid == numtasks -1) then
        right = 1
    else
        right = taskid +1
    endif
    if(taskid == 1) then
        left = numtasks -1
    else
        left = taskid -1
    endif

! open input files
! only need these on the master

    if(taskid == 0) then
        open ( unit=10, file='sami3-4.00.namelist'  )
        open ( unit=20, file='deni-init.inp'        )
        open ( unit=30, file='ichem.inp'            )
        open ( unit=50, file='phabsdt_euvac.inp'    )  ! euvac
        open ( unit=60, file='phiondt_euvac.inp'    )  ! euvac
        open ( unit=61, file='phionnt.inp'          )
        open ( unit=65, file='euvflux_euvac.inp'    )  ! euvac
        open ( unit=66, file='thetant.inp'          )
        open ( unit=67, file='zaltnt.inp'           )
    endif

    call initial

    call cpu_time(code_start)

! We are out of initial now.
! We will overwrite the values of
! deni, vsi, ti, te if this is a restart (restart = true)

    if(restart) then
        if(taskid == 0) then
            print *,'doing restart'
            open ( unit=210, file='time.rst', form='formatted' )
            open ( unit=211, file='deni.rst', form='unformatted' )
            open ( unit=212, file='vsi.rst', form='unformatted' )
            open ( unit=213, file='ti.rst', form='unformatted' )
            open ( unit=214, file='te.rst', form='unformatted' )

! note: this is good for 24 hr restarts only

!            read(210,*) hrinit  
            read(211) deniout
!            read(212) vsiout
            read(213) tiout
            read(214) teout

            close (210)
            close (211)
            close (212)
            close (213)
            close (214)

            do iwrk = 1,numworkers
                do nntmp = 1,nion
                    do ktmp = 1,nl
                        do jtmp = 1,nf
                            do itmp = 1,nz
                                denitmp(itmp,jtmp,ktmp) &
                                =  deniout(itmp,jtmp,ktmp,nntmp,iwrk)
                                titmp(itmp,jtmp,ktmp) &
                                =  tiout(itmp,jtmp,ktmp,nntmp,iwrk)
                                vsitmp(itmp,jtmp,ktmp) &
                                =  vsiout(itmp,jtmp,ktmp,nntmp,iwrk)
                            enddo
                        enddo
                    enddo

                    call mpi_send(denitmp, nz*nf*nl, MPI_REAL, iwrk, 9, &
                    MPI_COMM_WORLD, ierr)
                    call mpi_send(titmp, nz*nf*nl, MPI_REAL, iwrk, 9, &
                    MPI_COMM_WORLD, ierr)
                    call mpi_send(vsitmp, nz*nf*nl, MPI_REAL, iwrk, 9, &
                    MPI_COMM_WORLD, ierr)
                enddo

                do ktmp = 1,nl
                    do jtmp = 1,nf
                        do itmp = 1,nz
                            te(itmp,jtmp,ktmp) &
                            =  teout(itmp,jtmp,ktmp,iwrk)
                        enddo
                    enddo
                enddo
                call mpi_send(te, nz*nf*nl, MPI_REAL, iwrk, 9, &
                MPI_COMM_WORLD, ierr)

            enddo
        endif

    ! Now let's get those files

        if(taskid > 0 .AND. taskid <= numworkers) then

            do nntmp = 1,nion
                call mpi_recv(denitmp, nz*nf*nl, MPI_REAL, 0, 9, &
                MPI_COMM_WORLD, status, ierr)
                call mpi_recv(titmp, nz*nf*nl, MPI_REAL, 0, 9, &
                MPI_COMM_WORLD, status, ierr)
                call mpi_recv(vsitmp, nz*nf*nl, MPI_REAL, 0, 9, &
                MPI_COMM_WORLD, status, ierr)
                do ktmp = 1,nl
                    do jtmp = 1,nf
                        do itmp = 1,nz
                            deni(itmp,jtmp,ktmp,nntmp) &
                            =  denitmp(itmp,jtmp,ktmp)
                            ti(itmp,jtmp,ktmp,nntmp) &
                            =  titmp(itmp,jtmp,ktmp)
                            vsi(itmp,jtmp,ktmp,nntmp) &
                            =  vsitmp(itmp,jtmp,ktmp)
                        enddo
                    enddo
                enddo
            enddo
            call mpi_recv(te, nz*nf*nl, MPI_REAL, 0, 9, &
            MPI_COMM_WORLD, status, ierr)

        endif

    ! tell the workers the starting time
    ! this call has to be seen by the master and workers

        call mpi_bcast(hrinit,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

    endif

! open output files

    if(taskid == 0) call open_u

    if(taskid == 0) then
        close (10)
        close (20)
        close (30)
        close (50)
        close (54)
        close (60)
        close (64)
        close (61)
        close (65)
        close (66)
        close (67)
        close (68)
    endif

!******************* master task *******************************************
          
    if(taskid == 0) then

    ! allocate the total matrices only on master

        allocate (denit(nz,nf,nlt,nion),dennt(nz,nf,nlt,nneut))
        allocate  (vsit(nz,nf,nlt,nion),sumvsit(nz,nf,nlt,nion))
        allocate  (nuintt(nz,nf,nlt,nion),losstt(nz,nf,nlt,nion))
        allocate  (tit(nz,nf,nlt,nion))
        allocate (ut(nz,nf,nlt),vt(nz,nf,nlt),vpit(nz,nf,nlt))
        allocate (tet(nz,nf,nlt),tnt(nz,nf,nlt))
        allocate (u1t(nz,nf,nlt),u2t(nz,nf,nlt),u3t(nz,nf,nlt), &
                  u4t(nz,nf,nlt),u5t(nz,nf,nlt))
        allocate (vnqt(nz,nf,nlt),vnpt(nz,nf,nlt),vnphit(nz,nf,nlt))
        allocate (jpt(nz,nf,nlt),jphit(nz,nf,nlt))
        allocate (u1pt(nz,nf,nlt),u2st(nz,nf,nlt),u3ht(nz,nf,nlt))
        allocate (sigmapict(nz,nf,nlt),sigmahict(nz,nf,nlt))
        allocate (sigmapt(nz,nf,nlt),sigmaht(nz,nf,nlt))
        allocate (sigma0t(nz,nf,nlt))
        allocate (gpt(nz,nf,nlt))
!        allocate (xnuin0t(nz,nf,nlt),xnuin1t(nz,nf,nlt),xnuin2t(nz,nf,nlt))

        allocate (dtc(numworkers))

        hrut    = hrinit
        timemax = hrmax * sphr
        istep   = 0
        tneut   = 0.
        tpot    = 0.
        thardy  = 0.
        time    = 0.
        ntm     = 0
!        ntmmax  = min(maxstep, int( (hrmax-hrpr) / dthr ))
        ntmmax  = min(maxstep, int( (hrmax-hrpr) / dthr )) + 1

        print *,'max,ntmmax ',max(dt0/3600.,dthr),ntmmax

        ifintot  = numworkers
        ifintot1 = numworkers
        ifintot2 = numworkers

        tflag  = .TRUE. 
        icnt10 =  0

        do while ( tflag )

            do iwrk = 1,numworkers
                call mpi_iprobe(iwrk,10,MPI_COMM_WORLD,flagit10, &
                status,ierr)
                if (flagit10) then
                    icnt10 = icnt10 + 1
                    call mpi_recv(xxx,1,MPI_REAL,iwrk,10, &
                    MPI_COMM_WORLD,status,ierr)
                endif
                if (icnt10 == numworkers) tflag= .FALSE. 
            enddo

        ! Now wait to receive back the results from each worker task

            do  iwrk = 1, numworkers
                source = iwrk
                dest = source

                call mpi_iprobe(source, 2, &
                MPI_COMM_WORLD, flagit, status, ierr)

                if(flagit .AND. ifintot2 > 0) then

                    call mpi_recv(hipcp, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hihcm, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hipcphi, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hidphig, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hidpg, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hidphiv, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hidpv, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hipc, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hihc, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hidv, nf*nl, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)


                    call mpi_recv(hrut, 1, MPI_REAL, iwrk, 2, &
                    MPI_COMM_WORLD, status, ierr)

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            hipcpt(j,kk)   = hipcp(j,k)
                            hihcmt(j,kk)   = hihcm(j,k)
                            hipcphit(j,kk) = hipcphi(j,k)
                            hidphigt(j,kk) = hidphig(j,k)
                            hidpgt(j,kk)   = hidpg(j,k)
                            hidphivt(j,kk) = hidphiv(j,k)
                            hidpvt(j,kk)   = hidpv(j,k)
                            hipct(j,kk)    = hipc(j,k)
                            hihct(j,kk)    = hihc(j,k)
                            hidvt(j,kk)    = hidv(j,k)
                        enddo
                    enddo

                    ifintot2 = ifintot2 - 1

                endif

                if ( ifintot2 == 0 ) then
                    ifintot2 = numworkers
                    do jwrk = 1,numworkers
                        call mpi_send(tpot,1,MPI_REAL,jwrk,3, &
                        MPI_COMM_WORLD,ierr)
                    enddo
                    if ( tpot >= tphi .OR. tpot == 0) then
                      print *,'calling potpphi',tpot,dt
                      call potpphi(phi,dphi,hrut)
                      do jwrk = 1,numworkers
                          call mpi_send(phi,nnx*nny,MPI_REAL,jwrk,3, &
                          MPI_COMM_WORLD,ierr)
                      enddo
                      tpot = 0.
                    endif
                    tpot  = tpot + dt
                endif


                call mpi_iprobe(source, 0, &
                MPI_COMM_WORLD, flagit, status, ierr)

                if(flagit .AND. ifintot > 0) then
                                      
                ! This is just for outputting the data
                ! only sent as often as data dumps are requested

                    call mpi_recv(time, 1, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(hrut, 1, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(istep, 1, MPI_INTEGER, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    do nntmp = 1,nion
                        call mpi_recv(denitmp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                      MPI_COMM_WORLD, status, ierr)
                        call mpi_recv(denntmp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                      MPI_COMM_WORLD, status, ierr)
                        call mpi_recv(titmp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                      MPI_COMM_WORLD, status, ierr)
                        call mpi_recv(vsitmp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                      MPI_COMM_WORLD, status, ierr)
                        call mpi_recv(nuintmp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                      MPI_COMM_WORLD, status, ierr)
                        call mpi_recv(losstmp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                      MPI_COMM_WORLD, status, ierr)
                        do itmp = 1,nz
                            do jtmp = 1,nf
                                do ktmp = 1,nl
                                    deni(itmp,jtmp,ktmp,nntmp) &
                                      =  denitmp(itmp,jtmp,ktmp)
                                    denn(itmp,jtmp,ktmp,nntmp) &
                                      =  denntmp(itmp,jtmp,ktmp)
                                    ti(itmp,jtmp,ktmp,nntmp) &
                                      =  titmp(itmp,jtmp,ktmp)
                                    vsi(itmp,jtmp,ktmp,nntmp) &
                                      =  vsitmp(itmp,jtmp,ktmp)
                                    nuin_tot(itmp,jtmp,ktmp,nntmp) &
                                      =  nuintmp(itmp,jtmp,ktmp)
                                    loss_tot(itmp,jtmp,ktmp,nntmp) &
                                      =  losstmp(itmp,jtmp,ktmp)
                                enddo
                            enddo
                        enddo
                    enddo
                    call mpi_recv(te, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u1p, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u2s, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u3h, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u1, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u2, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u3, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u4, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(u5, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(gp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
!!$                    call mpi_recv(xnuin0, nz*nf*nl, MPI_REAL, iwrk, 0, &
!!$                                  MPI_COMM_WORLD, status, ierr)
!!$                    call mpi_recv(xnuin1, nz*nf*nl, MPI_REAL, iwrk, 0, &
!!$                                  MPI_COMM_WORLD, status, ierr)
!!$                    call mpi_recv(xnuin2, nz*nf*nl, MPI_REAL, iwrk, 0, &
!!$                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(sigma0, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(sigmap, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(sigmah, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(sigmapic, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(sigmahic, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(vnq, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(vnp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(vnphi, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(jp, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(jphi, nz*nf*nl, MPI_REAL, iwrk, 0, &
                                  MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do nn = 1,nion
                        do k = 2,nl-1
                            kk = (iwrk-1)*(nl -2) +k -1
                            if(kk == 0) kk = nlt
                            if(kk == nltp1) kk = 1
                            do j = 1,nf
                                do i = 1,nz
                                    denit(i,j,kk,nn)  = deni(i,j,k,nn)
                                    dennt(i,j,kk,nn)  = denn(i,j,k,nn)
                                    tit(i,j,kk,nn)    = ti(i,j,k,nn)
                                    vsit(i,j,kk,nn)   = vsi(i,j,k,nn)
                                    nuintt(i,j,kk,nn) = nuin_tot(i,j,k,nn)
                                    losstt(i,j,kk,nn) = loss_tot(i,j,k,nn)
                                enddo
                            enddo
                        enddo

                    ! Put the submatrices into the total matrix for restart

                        do k = 1,nl
                            do j = 1,nf
                                do i = 1,nz
                                    deniout(i,j,k,nn,iwrk) = deni(i,j,k,nn)
                                    tiout(i,j,k,nn,iwrk)   = ti(i,j,k,nn)
                                    vsiout(i,j,k,nn,iwrk)  = vsi(i,j,k,nn)
                                enddo
                            enddo
                        enddo
                    enddo  ! for nion loop

                    do k = 1,nl
                        do j = 1,nf
                            do i = 1,nz
                                teout(i,j,k,iwrk) = te(i,j,k)
                            enddo
                        enddo
                    enddo

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                tet(i,j,kk)        = te(i,j,k)
                                u1pt(i,j,kk)       = u1p(i,j,k)
                                u2st(i,j,kk)       = u2s(i,j,k)
                                u3ht(i,j,kk)       = u3h(i,j,k)
                                u1t(i,j,kk)        = u1(i,j,k)
                                u2t(i,j,kk)        = u2(i,j,k)
                                u3t(i,j,kk)        = u3(i,j,k)
                                u4t(i,j,kk)        = u4(i,j,k)
                                u5t(i,j,kk)        = u5(i,j,k)
                                gpt(i,j,kk)        = gp(i,j,k)
!!$                                xnuin0t(i,j,kk)    = xnuin0(i,j,k)
!!$                                xnuin1t(i,j,kk)    = xnuin1(i,j,k)
!!$                                xnuin2t(i,j,kk)    = xnuin2(i,j,k)
                                sigma0t(i,j,kk)    = sigma0(i,j,k)
                                sigmapt(i,j,kk)    = sigmap(i,j,k)
                                sigmaht(i,j,kk)    = sigmah(i,j,k)
                                sigmapict(i,j,kk)  = sigmapic(i,j,k)
                                sigmahict(i,j,kk)  = sigmahic(i,j,k)
                                vnqt(i,j,kk)       = vnq(i,j,k)
                                vnpt(i,j,kk)       = vnp(i,j,k)
                                vnphit(i,j,kk)     = vnphi(i,j,k)
                                jpt(i,j,kk)        = jp(i,j,k)
                                jphit(i,j,kk)      = jphi(i,j,k)
                            enddo
                        enddo
                    enddo
                                     
                    ifintot = ifintot -1

                endif

                call mpi_iprobe(source, 1, &
                MPI_COMM_WORLD, flagit1, status, ierr)
                   
                if(flagit1 .AND. ifintot1 > 0) then
                    call mpi_recv(dtmp, 1, MPI_REAL, iwrk, 1, &
                                  MPI_COMM_WORLD, status, ierr)
                    dtc(iwrk) = dtmp
                    call mpi_recv(time, 1, MPI_REAL, iwrk, 1, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(istep, 1, MPI_INTEGER, iwrk, 1, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(deni_mnptmp,nz*nion,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(ti_mnptmp,nz*nion,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(te_mnptmp,nz,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, status, ierr)

                    if ( ifintot1 == numworkers ) then
                        do ni = nion1,nion2
                            do i = 1,nz
                                deni_mnp(i,ni) = 0.
                                ti_mnp(i,ni)   = 0.
                            enddo
                        enddo

                        do i = 1,nz
                            te_mnp(i)     = 0.
                        enddo
                    endif

                    do ni = nion1,nion2
                        do i = 1,nz
                            deni_mnp(i,ni) = deni_mnp(i,ni) + &
                            deni_mnptmp(i,ni)/numworkers
                            ti_mnp(i,ni)   = ti_mnp(i,ni) + &
                            ti_mnptmp(i,ni)/numworkers
                        enddo
                    enddo

                    do i = 1,nz
                        te_mnp(i)     = te_mnp(i) + te_mnptmp(i)/numworkers
                    enddo

                    ifintot1 = ifintot1 - 1

                endif

            enddo               ! end worker loop


        ! if we are here, we should have gathered up all the data

            if(ifintot == 0) then
                ifintot = numworkers
                ntm = ntm + 1
                call output ( hrut,ntm,istep,phi,denit,dennt,vsit,nuintt, &
                              losstt,sumvsit,tit,ut,vt,vpit,tet,tnt,u1t, &
                              u2t,u3t,u4t,u5t,vnqt,vnpt,vnphit,jpt,jphit, &
                              u1pt,u2st,u3ht,sigmapict,sigmahict, &
                              sigmapt,sigmaht,sigma0t,gpt)
!                             ,xnuin0t,xnuin1t,xnuin2t)
                                         
            ! write the restart files and close those files

                open ( unit=210,  file='time.rst', form='formatted' )
                open ( unit=211,  file='deni.rst', form='unformatted' )
                open ( unit=212,  file='vsi.rst', form='unformatted' )
                open ( unit=213,  file='ti.rst', form='unformatted' )
                open ( unit=214,  file='te.rst', form='unformatted' )
                open ( unit= 2322,file='dphi.rst',form='unformatted')

                write(210,*) hrut
                write(211)   deniout
!                write(212)   vsiout
                write(213)   tiout
                write(214)   teout
                write(2322)  dphi

                close (210)
                close (211)
                close (212)
                close (213)
                close (214)
                close (2322)

            endif

        ! Need to fix up dt calculation

            if(ifintot1 == 0) then
                ifintot1 = numworkers
                dt = min(minval(dtc), dt0)  ! fix #9: amin1 -> min
!                print *,'dt = ',dt
                do  iwrk = 1,numworkers
                    call mpi_send(dt, 1, MPI_REAL, iwrk, 1, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(deni_mnp,nz*nion,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(ti_mnp,nz*nion,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(te_mnp,nz,MPI_REAL,iwrk,1, &
                                  MPI_COMM_WORLD, ierr)
                enddo
            endif

        enddo       ! end while (tflag)
              
        print *, 'MASTER: All Done!'

    ! close files

        close (10)
        close (20)
        close (40)
        close (70)
        close (71)
        close (72)
        close (73)
        close (74)
        close (75)
        close (78)
        close (79)
        close (80)
        close (90)
        close (91)
        close (92)
        close (93)
        close (94)
        close (95)
        close (96)
        close (97)
        close (98)
        close (81)
        close (82)
        close (83)
        close (84)
        close (85)
        close (86)
        close (87)
        close (88)
        close (711)
        close (712)
        close (713)
        close (714)
        close (715)
        close (1712)
        close (1713)
        close (1714)
        close (1715)
        close (569)
        close (716)
        close (717)
        close (1718)
        close (811)
        close (812)
        close (813)
        close (814)
        close (815)
        close (816)
        close (817)
        close (911)
        close (912)
        close (913)
        close (914)
        close (915)
        close (916)
        close (917)
        close (384)
        close (385)
        close (386)
        close (196)
        close (197)
        close (198)
        close (201)
        close (202)
        close (491)
        close (492)
        close (493)
        close (494)
        close (495)
        close (496)
        close (497)
        close (498)

    endif

!******************* end master task ***************************************

!******************** worker task *******************************************

    if(taskid > 0) then

    ! field line loop: actual run

        hrut    = hrinit
        timemax = hrmax * sphr
        istep   = 0
        tneut   = 0.
        thardy  = 0.
        time    = 0.
        ttflag  = .TRUE. 
        lprnt   = .FALSE.
        ntm     = 0
        nfile   = 2
        iday0   = 0
        dt_old  = dt0

        ntmmax  = min(maxstep, int( (hrmax-hrpr) / dthr)) + 1
        print *,'iwrk ',taskid,ntmmax
                 
    ! initialize neutrals
    ! neutral density, temperature, and neutral wind
    ! already done in initialization

        if (restart) then
            do nll = 1,nl
                call neutambt (hrinit,nll,1)
            enddo
        endif

        do while ( istep <= maxstep .AND. ttflag )

    ! parallel transport

        iday  = max(min(ceiling(time / 86400.), 8), 1)  ! fix #9: max0/min0 removed in F95

       ! EUVAC
         
        if ( iday > iday0 ) then   

          print *,'in if parallel - iday,iday0',iday,iday0

          p  = 0.5 * ( f10p7(iday) + fbar(iday) )
      
          do i = 1,linesuv_euvac
            f74   = fluxdat_euvac(i,1) * euv_fac
            ai    = fluxdat_euvac(i,2)
            flux_euvac(i) = max(0.,f74 * ( 1. + ai * ( p - 80.) ) * 1.e9)
          enddo

          print *,'max flux_euvac', maxval(flux_euvac)

          do i = 1,7
            aap(i) = ap(iday)
          enddo

          iday0 = iday
  
        endif

    ! Below is  nll = 2,nl-1 because of guard cells

            do nll = 1,nl
                do nfl = 1,nf
                    call zenith (hrut,nfl,nll,iday)
                    call transprt (nfl,nll)
                enddo
            enddo

    ! Do data exchanges between guard cells

    ! buffer and send to the LEFT

            do k = 1,nion
                do j = 1,nf
                    do i = 1,nz
                        tl1s(i,j,k) = deni(i,j,2,k)
                    enddo
                enddo
            enddo
            do k = nion+1,nion+nion
                do j = 1,nf
                    do i = 1,nz
                        tl1s(i,j,k) = ti(i,j,2,k-nion)
                    enddo
                enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
                do i = 1,nz
                    tl1s(i,j,k) = te(i,j,2)
                enddo
            enddo

            call mpi_sendrecv(tl1s, (nion+nion+1)*nz*nf, MPI_REAL, &
            left, 0, tl1r, (nion+nion+1)*nz*nf, MPI_REAL, &
            right, 0, MPI_COMM_WORLD, status, ierr)

        ! Now everybody receives

            do k = 1,nion
                do j = 1,nf
                    do i = 1,nz
                        deni(i,j,nl,k) = tl1r(i,j,k)
                    enddo
                enddo
            enddo
            do k = nion+1,nion+nion
                do j = 1,nf
                    do i = 1,nz
                        ti(i,j,nl,k-nion) = tl1r(i,j,k)
                    enddo
                enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
                do i = 1,nz
                    te(i,j,nl) = tl1r(i,j,k)
                enddo
            enddo

        ! Buffer and send to the RIGHT

            do k = 1,nion
                do j = 1,nf
                    do i = 1,nz
                        tr1s(i,j,k) = deni(i,j,nl-1,k)
                    enddo
                enddo
            enddo
            do k = nion+1,nion+nion
                do j = 1,nf
                    do i = 1,nz
                        tr1s(i,j,k) = ti(i,j,nl-1,k-nion)
                    enddo
                enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
                do i = 1,nz
                    tr1s(i,j,k) = te(i,j,nl-1)
                enddo
            enddo

            call mpi_sendrecv(tr1s, (nion+nion+1)*nz*nf, MPI_REAL, &
            right, 0, tr1r, (nion +nion +1)*nz*nf, MPI_REAL, &
            left, 0, MPI_COMM_WORLD, status, ierr)

            do k = 1,nion
                do j = 1,nf
                    do i = 1,nz
                        deni(i,j,1,k) = tr1r(i,j,k)
                    enddo
                enddo
            enddo
            do k = nion+1,nion+nion
                do j = 1,nf
                    do i = 1,nz
                        ti(i,j,1,k-nion) = tr1r(i,j,k)
                    enddo
                enddo
            enddo
            k = nion + nion + 1
            do j = 1,nf
                do i = 1,nz
                    te(i,j,1) = tr1r(i,j,k)
                enddo
            enddo

        ! We are now finished exchanging guard cell data

        ! Sending hipcp and hidphig to master to calculate the
        ! potential

            call mpi_send(hipcp, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hihcm, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hipcphi, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hidphig, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hidpg, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hidphiv, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hidpv, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hipc, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hihc, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(hidv, nf*nl, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)


        ! send master current time: hrut

            call mpi_send(hrut, 1, MPI_REAL, 0, 2, &
                          MPI_COMM_WORLD, ierr)

!       now get the potential from master

            call mpi_recv(tpot, 1, MPI_REAL, 0, 3, &
                          MPI_COMM_WORLD, status, ierr)
 
           if ( tpot >= tphi .OR. tpot == 0 ) &
           call mpi_recv(phi, nnx*nny, MPI_REAL, 0, 3, &
                          MPI_COMM_WORLD, status, ierr)

!       perpendicular transport

            call exb  (hrut,phi)
            call neut (hrut)
            call courant

!       average magnetic pole grid values (deni,Ti,Te)

            j0           = nf

            do ni = nion1,nion2
                do i = 1,nz
                    deni_mnp0 = 0.
                    ti_mnp0   = 0.
                    do k = 2,nl-1
                        if ( alts (i,j0,k) < alt_crit_avg) then
                            deni_mnp0      = deni_mnp0 + deni(i,j0,k,ni)
                            ti_mnp0        = ti_mnp0 + ti(i,j0,k,ni)
                        endif
                    enddo
                    deni_mnp(i,ni) = deni_mnp0 / real(nl-2)  ! fix #9
                    ti_mnp(i,ni)   = ti_mnp0 / real(nl-2)  ! fix #9
                enddo
            enddo

            do i = 1,nz
                te_mnp0 = 0.
                do k = 2,nl-1
                    if ( alts (i,j0,k) < alt_crit_avg) then
                        te_mnp0     = te_mnp0 + te(i,j0,k)
                    endif
                enddo
                te_mnp(i)   = te_mnp0 / real(nl-2)  ! fix #9
            enddo


        ! send local dt

            call mpi_send(dt, 1, MPI_REAL, 0, 1, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(time, 1, MPI_REAL, 0, 1, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(istep, 1, MPI_INTEGER, 0, 1, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(deni_mnp,nz*nion,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(ti_mnp,nz*nion,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, ierr)
            call mpi_send(te_mnp,nz,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, ierr)

        ! get global dt

            call mpi_recv(dt, 1, MPI_REAL, 0, 1, &
                          MPI_COMM_WORLD, status, ierr)
            call mpi_recv(deni_mnp,nz*nion,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, status, ierr)
            call mpi_recv(ti_mnp,nz*nion,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, status, ierr)
            call mpi_recv(te_mnp,nz,MPI_REAL,0,1, &
                          MPI_COMM_WORLD, status, ierr)

        ! update neutrals

            if( tneut >= 0.25 ) then
                do nll = 1,nl
                    call neutambt (hrut,nll,iday)
                enddo
                tneut = 0.
                if ( hrut < hrpr+hrinit ) then
                    print *,'No output yet: hr = ',hrut,dt
                endif
            endif

        ! update hardy

            if( thardy >= 0.25 ) then
                do nll = 1,nl
                  do nfl = 1,nf
                    call hardy_sami3 (hrut,nfl,nll,iday)
                    do i = 1,nz
                      tpn(i,nfl,nll) = preciprn(i,nfl,nll)
                      tps(i,nfl,nll) = preciprs(i,nfl,nll)
                    enddo
                  enddo
                enddo
                thardy = 0.
            endif

            if ( (lprnt .AND. hrut >= hrpr+hrinit) .OR. istep == 0 ) then

            ! We no longer call output from here, but send data to the MASTER
            ! The four things we want to send are  deni, ti, vsi, te

                ntm    = ntm + 1

                call mpi_send(time, 1, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(hrut, 1, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(istep, 1, MPI_INTEGER, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                do nntmp = 1,nion
                    do itmp = 1,nz
                        do jtmp = 1,nf
                            do ktmp = 1,nl
                                denitmp(itmp,jtmp,ktmp) &
                                  = deni(itmp,jtmp,ktmp,nntmp)
                                denntmp(itmp,jtmp,ktmp) &
                                  = denn(itmp,jtmp,ktmp,nntmp)
                                titmp(itmp,jtmp,ktmp) &
                                  = ti(itmp,jtmp,ktmp,nntmp)
                                vsitmp(itmp,jtmp,ktmp) &
                                  = vsi(itmp,jtmp,ktmp,nntmp)
                                nuintmp(itmp,jtmp,ktmp) &
                                  = nuin_tot(itmp,jtmp,ktmp,nntmp)
                                losstmp(itmp,jtmp,ktmp) &
                                  = loss_tot(itmp,jtmp,ktmp,nntmp)
                            enddo
                        enddo
                    enddo
                    call mpi_send(denitmp, nz*nf*nl, MPI_REAL, 0, 0, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(denntmp, nz*nf*nl, MPI_REAL, 0, 0, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(titmp, nz*nf*nl, MPI_REAL, 0, 0, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(vsitmp, nz*nf*nl, MPI_REAL, 0, 0, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(nuintmp, nz*nf*nl, MPI_REAL, 0, 0, &
                                  MPI_COMM_WORLD, ierr)
                    call mpi_send(losstmp, nz*nf*nl, MPI_REAL, 0, 0, &
                                  MPI_COMM_WORLD, ierr)
                enddo

                call mpi_send(te, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(u1p, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u2s, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u3h, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u1, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u2, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u3, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u4, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(u5, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
                call mpi_send(gp, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD,  ierr)
!!$                call mpi_send(xnuin0, nz*nf*nl, MPI_REAL, 0, 0, &
!!$                              MPI_COMM_WORLD,  ierr)
!!$                call mpi_send(xnuin1, nz*nf*nl, MPI_REAL, 0, 0, &
!!$                              MPI_COMM_WORLD,  ierr)
!!$                call mpi_send(xnuin2, nz*nf*nl, MPI_REAL, 0, 0, &
!!$                              MPI_COMM_WORLD,  ierr)
                call mpi_send(sigma0, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(sigmap, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(sigmah, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(sigmapic, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(sigmahic, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(vnq, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(vnp, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(vnphi, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(jp, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                call mpi_send(jphi, nz*nf*nl, MPI_REAL, 0, 0, &
                              MPI_COMM_WORLD, ierr)
                 
                dt      = dt_old
                lprnt   = .FALSE. 
                if ( ntm >= ntmmax ) ttflag = .FALSE. 

            endif

        ! time/step advancement

            istep  = istep + 1
            time   = time  + dt
            hrut   = time / sphr + hrinit
            tneut  = tneut + dt / sphr
            thardy = thardy + dt / sphr
            tprnt  = (hrpr + ntm * dthr) * 3600.
            timem  = time  - dt
            timep  = time  + dt
!            dts    = abs(dthr*3600. - tprnt)
            dt_old = dt
            if ( timem <= tprnt .AND. timep >= tprnt .AND. .not. lprnt ) then
!                dt     = dts
!                if ( dts == 0 ) dt = dt_old
                lprnt  = .TRUE. 
!                print *,'to print',ntm,timem,tprnt,timep
            endif

        enddo    ! end time loop

    xxx = 1.
    call mpi_send(xxx, 1, MPI_REAL, 0, 10, &
    MPI_COMM_WORLD, ierr)

    endif ! end worker task


    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    call mpi_finalize(ierr)
    print *,'done finalizing,taskid',taskid

    call cpu_time(code_finish)
    print *,'timing',code_start,code_finish,code_finish-code_start

    stop
    END PROGRAM

!******************************************
!******************************************

!            initial

!******************************************
!******************************************

    subroutine initial

!    use apex_module,only: apxparm ! (GL)
    use parameter_mod
    use variable_mod
    use namelist_mod
    use message_passing_mod
    use chemistry_mod
    use time_mod
    use photo_production_mod
    use atomic_mod
    use exb_mod
    use misc_mod
    use grid_mod
    use hardy_mod

    implicit none  ! fix #1

    include 'mpif.h'

    real, dimension(:,:,:), allocatable :: altstmp,glatstmp,glonstmp,altstmp2
    real, dimension(:,:,:), allocatable :: altstmp1,glatstmp1,glonstmp1
    real, dimension(:,:,:), allocatable :: altst,glatst,glonst
    real, dimension(:,:,:), allocatable :: dvec11,dvec12,dvec13
    real, dimension(:,:,:), allocatable :: dvec21,dvec22,dvec23
    real, dimension(:,:,:), allocatable :: dvec31,dvec32,dvec33
    real, dimension(:,:,:), allocatable :: baltst,blatst,blonst
    real, dimension(:,:,:), allocatable :: xst,yst,zst
    real, dimension(:,:,:), allocatable :: altptmp,blatptmp,blonptmp
    real, dimension(:,:,:), allocatable :: baltpt
    real, dimension(:,:,:), allocatable :: vpsnxt,vpsnyt,vpsnzt
    real, dimension(:,:,:), allocatable :: vhsnxt,vhsnyt,vhsnzt
    real, dimension(:,:,:), allocatable :: xpt,ypt,zpt
    real, dimension(:,:,:), allocatable :: bdirsxt,bdirsyt,bdirszt
    real, dimension(:,:,:), allocatable :: gsthetaxt
    real, dimension(:,:,:), allocatable :: gsthetayt
    real, dimension(:,:,:), allocatable :: gsthetazt
    real, dimension(:,:,:), allocatable :: gsphixt
    real, dimension(:,:,:), allocatable :: gsphiyt
    real, dimension(:,:,:), allocatable :: gsphizt
    real, dimension(:,:,:), allocatable :: gsrxt
    real, dimension(:,:,:), allocatable :: gsryt
    real, dimension(:,:,:), allocatable :: gsrzt
    real, dimension(:,:,:), allocatable :: xrgt
    real, dimension(:,:,:), allocatable :: volt
    real, dimension(:,:,:), allocatable :: xthgt
    real, dimension(:,:,:), allocatable :: xphigt
    real, dimension(:,:,:), allocatable :: bmst

    real, dimension(:,:,:), allocatable :: xnormst,ynormst,znormst
    real, dimension(:,:,:), allocatable :: xnormpt,ynormpt,znormpt
    real, dimension(:,:,:), allocatable :: xnormht,ynormht,znormht
    real, dimension(:,:,:), allocatable :: xstmp,ystmp,zstmp
    real, dimension(:,:,:), allocatable :: xptmp,yptmp,zptmp
    real, dimension(:,:,:), allocatable :: xhtmp,yhtmp,zhtmp
    real, dimension(:,:,:), allocatable :: vppnxt,vppnyt,vppnzt
    real, dimension(:,:,:), allocatable :: vhpnxt,vhpnyt,vhpnzt
    real, dimension(:,:,:), allocatable :: bdirpxt,bdirpyt,bdirpzt
    real, dimension(:,:,:), allocatable :: bdirhxt,bdirhyt,bdirhzt
    real, dimension(:,:,:), allocatable :: ehpxt,ehpyt,ehpzt
    real, dimension(:,:,:), allocatable :: ephxt,ephyt,ephzt
    real, dimension(:,:,:), allocatable :: delhpt
    real, dimension(:,:,:), allocatable :: delpht

    real    :: f1026(nz,nf,nl,91),f584(nz,nf,nl,91), &
               f304 (nz,nf,nl,91),f1216(nz,nf,nl,91)
    integer :: status(MPI_STATUS_SIZE),ierr
    real    :: zi(29),denii(29,7)
    real    :: phionrnt(linesnt,4)
    real    :: phionrdt_euvac(linesuv_euvac,5)

! Some local variables

    ! fix #1: previously implicit locals in initial
    integer :: i, j, k, n, l, nn, ni, nll, nfl
    integer :: j0, jj, iwrk, kk, ifintot
    real    :: slope, amuf, amimn, plon_ang
    real    :: deni_mnp0, ti_mnp0, te_mnp0


    namelist / go / maxstep,hrmax,dthr,hrpr,dt0, &
      rmin, altmin, fbar,f10p7,ap, iyear,nday,&
      imonth_1,iday_1,imonth_2,iday_2,&
      mmass, &
      nion1,nion2,hrinit,tvn0,tvexb0,ver,veh,vw, &
      gams, gamp, alt_crit_avg, blat_max, blat_min,&
      snn,stn,denmin,alt_crit,cqe, &
      hall,restart, &
      vexb_max, &
      lmadala,lweimer,decay_time,pcrit, &
      anu_drag0,kp,xhardy,tphi,&
      tmax,euv_fac,s5e_fac

    if ( taskid == 0 ) then
        allocate (altstmp(nz,nf,nl), glatstmp(nz,nf,nl),glonstmp(nz,nf,nl))
        allocate (altstmp1(nz,nf,nl), glatstmp1(nz,nf,nl),glonstmp1(nz,nf,nl))
        allocate (altstmp2(nz,nf,nl+1))
        allocate (altst(nz,nf,nlt),  glatst(nz,nf,nlt), glonst(nz,nf,nlt))
        allocate (dvec11(nz,nf,nlt), dvec12(nz,nf,nlt), dvec13(nz,nf,nlt))
        allocate (dvec21(nz,nf,nlt), dvec22(nz,nf,nlt), dvec23(nz,nf,nlt))
        allocate (dvec31(nz,nf,nlt), dvec32(nz,nf,nlt), dvec33(nz,nf,nlt))
        allocate (baltst(nz,nf,nlt), blatst(nz,nf,nlt), blonst(nz,nf,nlt))
        allocate (xst(nz,nf,nlt),    yst(nz,nf,nlt), zst(nz,nf,nlt))
        allocate (altptmp(nzp1,nfp1,nlp1), blatptmp(nzp1,nfp1,nlp1), blonptmp(nzp1,nfp1,nlp1))
        allocate (vpsnxt(nz,nf,nlt), vpsnyt(nz,nf,nlt), vpsnzt(nz,nf,nlt))
        allocate (vhsnxt(nz,nf,nlt), vhsnyt(nz,nf,nlt), vhsnzt(nz,nf,nlt))
        allocate (vppnxt(nzp1,nfp1,nlt), vppnyt(nzp1,nfp1,nlt), vppnzt(nzp1,nfp1,nlt))
        allocate (vhpnxt(nzp1,nfp1,nlt), vhpnyt(nzp1,nfp1,nlt), vhpnzt(nzp1,nfp1,nlt))
        allocate (xpt(nzp1,nfp1,nlt), ypt(nzp1,nfp1,nlt), zpt(nzp1,nfp1,nlt))
        allocate (bdirsxt(nz,nf,nlt), bdirsyt(nz,nf,nlt), bdirszt(nz,nf,nlt))
        allocate (bdirpxt(nz,nfp1,nlt), bdirpyt(nz,nfp1,nlt), bdirpzt(nz,nfp1,nlt))
        allocate (bdirhxt(nz,nf,nlt), bdirhyt(nz,nf,nlt), bdirhzt(nz,nf,nlt))
        allocate (ehpxt(nz,nfp1,nlt), ehpyt(nz,nfp1,nlt), ehpzt(nz,nfp1,nlt))
        allocate (ephxt(nz,nf,nlt), ephyt(nz,nf,nlt), ephzt(nz,nf,nlt))
        allocate (delhpt(nz,nfp1,nlt))
        allocate (delpht(nz,nf,nlt))
        allocate (gsthetaxt(nz,nf,nlt), gsthetayt(nz,nf,nlt), gsthetazt(nz,nf,nlt))
        allocate (gsphixt(nz,nf,nlt), gsphiyt(nz,nf,nlt), gsphizt(nz,nf,nlt))
        allocate (gsrxt(nz,nf,nlt), gsryt(nz,nf,nlt), gsrzt(nz,nf,nlt))
        allocate (volt(nz,nf,nlt))
        allocate (xrgt(nz,nf,nlt), xthgt(nz,nf,nlt), xphigt(nz,nf,nlt))
        allocate (bmst(nz,nf,nlt))
        allocate (baltpt(nzp1,nfp1,nlt))

        allocate (xnormst(nzp1,nf,nlt),ynormst(nzp1,nf,nlt),znormst(nzp1,nf,nlt))
        allocate (xnormpt(nz,nfp1,nlt),ynormpt(nz,nfp1,nlt),znormpt(nz,nfp1,nlt))
        allocate (xnormht(nz,nf,nlt+1),ynormht(nz,nf,nlt+1),znormht(nz,nf,nlt+1))
        allocate (xstmp(nzp1,nf,nl),ystmp(nzp1,nf,nl),zstmp(nzp1,nf,nl))
        allocate (xptmp(nz,nfp1,nl),yptmp(nz,nfp1,nl),zptmp(nz,nfp1,nl))
        allocate (xhtmp(nz,nf,nlp1),yhtmp(nz,nf,nlp1),zhtmp(nz,nf,nlp1))

    endif

! read in parameters and initial ion density data

    if(taskid == 0) then
        read(10,go)
    endif

! send the namelist data to all the other processors

!    print *,'start bcast 0'

    call mpi_bcast(maxstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hrmax,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dthr,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hrpr,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(dt0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(rmin,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(altmin,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(fbar,8,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(f10p7,8,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ap,8,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(iyear,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nday,8,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(imonth_1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(iday_1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(imonth_2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(iday_2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(mmass,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nion1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(nion2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(tvn0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(tvexb0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ver,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(veh,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vw,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(gams,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(gamp,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(alt_crit_avg,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(blat_max,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(blat_min,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(snn,nneut,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(stn,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(denmin,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(alt_crit,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(cqe,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hall,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(vexb_max,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lmadala,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(lweimer,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(decay_time,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(pcrit,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(hrinit,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(anu_drag0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(kp,8,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(xhardy,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(tphi,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(tmax,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(euv_fac,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(s5e_fac,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

!    print *,'end bcast 0'

    dt = dt0

    ami(pthp)  = 1.
    ami(pthep) = 4.
    ami(ptnp)  = 14.
    ami(ptop)  = 16.
    ami(ptn2p) = 28.
    ami(ptnop) = 30.
    ami(pto2p) = 32.

    amn(pth)  = 1.
    amn(pthe) = 4.
    amn(ptn)  = 14.
    amn(pto)  = 16.
    amn(ptn2) = 28.
    amn(ptno) = 30.
    amn(pto2) = 32.

    alpha0(pth)  = 0.67
    alpha0(pthe) = 0.21
    alpha0(ptn)  = 1.10
    alpha0(pto)  = 0.79
    alpha0(ptn2) = 1.76
    alpha0(ptno) = 1.74
    alpha0(pto2) = 1.59

    do i = 1,7
        aap(i) = ap(1)
    enddo

    do nn = 1,nneut
      do ni = 1,nion
        amuf    = ami(ni) * amn(nn) / ( ami(ni) + amn(nn) )
        amimn   = amn(nn) / ( ami(ni) + amn(nn) )
        nufacin(ni,nn) = 2.69e-9 / sqrt(amuf) * amimn * sqrt(alpha0(nn))
      enddo
    enddo

! read in initial density data

    if(taskid == 0) then
        do i = 1,29
            read(20,102) zi(i),(denii(i,j),j=1,7)
            102 format(1x,f7.1,1p7e8.1)
        enddo
    endif

    call mpi_bcast(zi,29, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(denii,29*7, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)


!    print *,'end bcast 1'

! read in chemistry data
! in format statement 104 need to 'hardwire' nneut (= 7)

    if(taskid == 0) then
        do k = 1,nchem
            read(30,103) (ichem(k,j),j=1,3)
            103 format(3i3)
        enddo
    endif
    call mpi_bcast(ichem,nchem*3, &
                   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)  ! fix #2: was MPI_REAL; ichem is integer
!    print *,'end bcast 2'

! generate the mesh data by everybody but the Master

! Calculate APEX magnetic field using apex_suns  4/5/2013 (GL)

!      write(6,"('initial call apxparm: iyear = ',f7.0)") float(iyear)
!      call apxparm(float(iyear))

! new re doug

    call loadapxsh(2010.)
    call cofrm(2010.)

    if (taskid == 0) call blonp0a
    if(taskid > 0)  call grid3_mpi

! Now wait to receive back the results from each worker task

    if(taskid == 0) then
                 
        ifintot = numworkers

        do while( ifintot > 0)

            do  iwrk = 1, numworkers
                source = iwrk
                dest = source
                            
                call mpi_iprobe(source, 0, &
                MPI_COMM_WORLD, flagit, status, ierr)
                               
                if(flagit .AND. ifintot > 0) then

                !  The three things we want to receive are  altpt blatpt blonpt

    print *,'rcving baltp ',taskid

                    call mpi_recv(altptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blatptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blonptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +  k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nzp1
                                baltpt(i,j,kk) = altptmp(i,j,k)
                                blatpt(i,j,kk) = blatptmp(i,j,k)
                                blonpt(i,j,kk) = blonptmp(i,j,k)
                            enddo
                        enddo
                    enddo

                  print *,'print recvd balt/blat/blonpt',taskid

                !  The three things we want to receive are  xpt ypt zpt

                    call mpi_recv(altptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blatptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blonptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nzp1
                                xpt(i,j,kk) = altptmp(i,j,k)
                                ypt(i,j,kk) = blatptmp(i,j,k)
                                zpt(i,j,kk) = blonptmp(i,j,k)
                            enddo
                        enddo
                    enddo
!                  print *,'print recvd x/y/zpt',taskid

                 print *,'recv pp '

                ! We want to receive pp

                    call mpi_recv(altptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nzp1
                                ppt(i,j,kk)   = altptmp(i,j,k)
                            enddo
                        enddo
                    enddo

              print *,'start recv dvec',taskid

                ! The three things we want to receive are  dvec11/12/13

                    call mpi_recv(altstmp1, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp1, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp1, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                dvec11(i,j,kk) = altstmp1(i,j,k)
                                dvec12(i,j,kk) = glatstmp1(i,j,k)
                                dvec13(i,j,kk) = glonstmp1(i,j,k)
                            enddo
                        enddo
                    enddo
!              print *,' recv dvec1',taskid


                ! The three things we want to receive are  dvec21/22/23

                    call mpi_recv(altstmp1, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp1, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp1, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                dvec21(i,j,kk) = altstmp1(i,j,k)
                                dvec22(i,j,kk) = glatstmp1(i,j,k)
                                dvec23(i,j,kk) = glonstmp1(i,j,k)
                            enddo
                        enddo
                    enddo

!              print *,' recv dvec2',taskid

                ! The three things we want to receive are  dvec31/32/33

                    call mpi_recv(altstmp1, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp1, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp1, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                dvec31(i,j,kk) = altstmp1(i,j,k)
                                dvec32(i,j,kk) = glatstmp1(i,j,k)
                                dvec33(i,j,kk) = glonstmp1(i,j,k)
                            enddo
                        enddo
                    enddo

!              print *,' recv dvec3',taskid


                ! The three things we want to receive are  altst glatst glonst

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                altst(i,j,kk)  = altstmp(i,j,k)
                                glatst(i,j,kk) = glatstmp(i,j,k)
                                glonst(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                ! The three things we want to receive are  baltst blatst blonst

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                baltst(i,j,kk) = altstmp(i,j,k)
                                blatst(i,j,kk) = glatstmp(i,j,k)
                                blonst(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                ! The three things we want to receive are  xst yst zst

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                xst(i,j,kk) = altstmp(i,j,k)
                                yst(i,j,kk) = glatstmp(i,j,k)
                                zst(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                ! The three things we want to receive are  xrg,xthg,xphig

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                ! Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                xrgt(i,j,kk)   = altstmp(i,j,k)
                                xthgt(i,j,kk)  = glatstmp(i,j,k)
                                xphigt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) +k -1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                bmst(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are vpsnx vpsny vpsnz

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                vpsnxt(i,j,kk) = altstmp(i,j,k)
                                vpsnyt(i,j,kk) = glatstmp(i,j,k)
                                vpsnzt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are vhsnx vhsny vhsnz

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                vhsnxt(i,j,kk) = altstmp(i,j,k)
                                vhsnyt(i,j,kk) = glatstmp(i,j,k)
                                vhsnzt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are bdirsx bdirsy bdirsz

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                bdirsxt(i,j,kk) = altstmp(i,j,k)
                                bdirsyt(i,j,kk) = glatstmp(i,j,k)
                                bdirszt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are gsthetax/y/z

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                gsthetaxt(i,j,kk) = altstmp(i,j,k)
                                gsthetayt(i,j,kk) = glatstmp(i,j,k)
                                gsthetazt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are gsphix/y/z

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                gsphixt(i,j,kk) = altstmp(i,j,k)
                                gsphiyt(i,j,kk) = glatstmp(i,j,k)
                                gsphizt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo


                !  The three things we want to receive are gsrx/y/z

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                gsrxt(i,j,kk) = altstmp(i,j,k)
                                gsryt(i,j,kk) = glatstmp(i,j,k)
                                gsrzt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  receive cell volume

                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                volt(i,j,kk) = altstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are x/y/znorms
        print *,'recv xyznorms',iwrk
                    call mpi_recv(xstmp, nzp1*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(ystmp, nzp1*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(zstmp, nzp1*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nzp1
                                xnormst(i,j,kk) = xstmp(i,j,k)
                                ynormst(i,j,kk) = ystmp(i,j,k)
                                znormst(i,j,kk) = zstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are x/y/znormp
        print *,'recv xyznormp'
                    call mpi_recv(xptmp, nz*nfp1*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(yptmp, nz*nfp1*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(zptmp, nz*nfp1*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nz
                                xnormpt(i,j,kk) = xptmp(i,j,k)
                                ynormpt(i,j,kk) = yptmp(i,j,k)
                                znormpt(i,j,kk) = zptmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are x/y/znormh
!        print *,'recv xyznormh'
                    call mpi_recv(xhtmp, nz*nf*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(yhtmp, nz*nf*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(zhtmp, nz*nf*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                xnormht(i,j,kk) = xhtmp(i,j,k)
                                ynormht(i,j,kk) = yhtmp(i,j,k)
                                znormht(i,j,kk) = zhtmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are vppnx/y/z
        print *,'recv vppnx/y/z'
                    call mpi_recv(altptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blatptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blonptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nzp1
                                vppnxt(i,j,kk) = altptmp(i,j,k)
                                vppnyt(i,j,kk) = blatptmp(i,j,k)
                                vppnzt(i,j,kk) = blonptmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are vhpnx/y/z
        print *,'recv vhpnx/y/z'
                    call mpi_recv(altptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blatptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(blonptmp, nzp1*nfp1*nlp1, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nzp1
                                vhpnxt(i,j,kk) = altptmp(i,j,k)
                                vhpnyt(i,j,kk) = blatptmp(i,j,k)
                                vhpnzt(i,j,kk) = blonptmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are bdirp/y/z
        print *,'recv bdirpx/y/z'
                    call mpi_recv(xptmp, nz*nfp1*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(yptmp, nz*nfp1*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(zptmp, nz*nfp1*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nz
                                bdirpxt(i,j,kk) = xptmp(i,j,k)
                                bdirpyt(i,j,kk) = yptmp(i,j,k)
                                bdirpzt(i,j,kk) = zptmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are bdirh/y/z
        print *,'recv bdirhx/y/z'
                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                bdirhxt(i,j,kk) = altstmp(i,j,k)
                                bdirhyt(i,j,kk) = glatstmp(i,j,k)
                                bdirhzt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are ehp/y/z
        print *,'recv ehpx/y/z'
                    call mpi_recv(xptmp, nz*nfp1*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(yptmp, nz*nfp1*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(zptmp, nz*nfp1*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nz
                                ehpxt(i,j,kk) = xptmp(i,j,k)
                                ehpyt(i,j,kk) = yptmp(i,j,k)
                                ehpzt(i,j,kk) = zptmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are eph/y/z
        print *,'recv ephx/y/z'
                    call mpi_recv(altstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glatstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)
                    call mpi_recv(glonstmp, nz*nf*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                ephxt(i,j,kk) = altstmp(i,j,k)
                                ephyt(i,j,kk) = glatstmp(i,j,k)
                                ephzt(i,j,kk) = glonstmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are delhp'
        print *,'recv delhp'
                    call mpi_recv(xptmp, nz*nfp1*nl, MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nfp1
                            do i = 1,nz
                                delhpt(i,j,kk) = xptmp(i,j,k)
                            enddo
                        enddo
                    enddo

                !  The three things we want to receive are delph
        print *,'recv delph'
                    call mpi_recv(altstmp2, nz*nf*(nl+1), MPI_REAL, &
                                  iwrk, 0, MPI_COMM_WORLD, status, ierr)  ! fix #5: was nz*nf*nl+1

                !  Put the submatrices into the correct matrix

                    do k = 2,nl-1
                        kk = (iwrk-1)*(nl -2) + k - 1
                        if(kk == 0) kk = nlt
                        if(kk == nltp1) kk = 1
                        do j = 1,nf
                            do i = 1,nz
                                delpht(i,j,kk) = altstmp2(i,j,k)
                            enddo
                        enddo
                    enddo


                    ifintot = ifintot -1
                                      
                endif

            enddo                  ! end worker loop

        enddo

!     calculate plon (for angut in potential)

                 plat     = 1.57  ! not used
                 plon_ang = glonst(nz/2,1,1) - blonst(nz/2,1,1)
                 plon     = pie * plon_ang / 180.

                 print *,'plon_ang, plon = ',plon_ang,plon

    ! if we are here, we should have gathered up all the data

    ! output grid data

            open ( unit=69, file='dvec11.dat',form='unformatted' )
            open ( unit=76, file='dvec12.dat',form='unformatted' )
            open ( unit=77, file='dvec13.dat',form='unformatted' )
            write(69) dvec11
            write(76) dvec12
            write(77) dvec13
            close(69)
            close(76)
            close(77)

            open ( unit=69, file='dvec21.dat',form='unformatted' )
            open ( unit=76, file='dvec22.dat',form='unformatted' )
            open ( unit=77, file='dvec23.dat',form='unformatted' )
            write(69) dvec21
            write(76) dvec22
            write(77) dvec23
            close(69)
            close(76)
            close(77)

            open ( unit=69, file='dvec31.dat',form='unformatted' )
            open ( unit=76, file='dvec32.dat',form='unformatted' )
            open ( unit=77, file='dvec33.dat',form='unformatted' )
            write(69) dvec31
            write(76) dvec32
            write(77) dvec33
            close(69)
            close(76)
            close(77)

            open ( unit=69, file='zaltu.dat',form='unformatted' )
            open ( unit=76, file='glatu.dat',form='unformatted' )
            open ( unit=77, file='glonu.dat',form='unformatted' )
            write(69) altst
            write(76) glatst
            write(77) glonst
            close(69)
            close(76)
            close(77)

            open (144,file='glons_cg.rst',form='unformatted')
            write(144) glonst
            close(144)


            open ( unit=69, file='baltu.dat',form='unformatted' )
            open ( unit=76, file='blatu.dat',form='unformatted' )
            open ( unit=77, file='blonu.dat',form='unformatted' )
            write(69) baltst
            write(76) blatst
            write(77) blonst
            close(69)
            close(76)
            close(77)

            open (144,file='blons_cg.rst',form='unformatted')
            write(144) blonst
            close(144)

            open ( unit=69, file='xsu.dat',form='unformatted' )
            open ( unit=76, file='ysu.dat',form='unformatted' )
            open ( unit=77, file='zsu.dat',form='unformatted' )
            write(69) xst
            write(76) yst
            write(77) zst
            close(69)
            close(76)
            close(77)

            open ( unit=169, file='baltpu.dat'   ,form='unformatted' )
            open ( unit=176, file='blatpu.dat'   ,form='unformatted' )
            open ( unit=177, file='blonpu.dat'   ,form='unformatted' )
            write(169) baltpt
            write(176) blatpt
            write(177) blonpt
            close(169)
            close(176)
            close(177)

            open (144,file='blonp_cg.rst',form='unformatted')
            write(144) blonpt
            close(144)

            open ( unit=169, file='xpu.dat'   ,form='unformatted' )
            open ( unit=176, file='ypu.dat'   ,form='unformatted' )
            open ( unit=177, file='zpu.dat'   ,form='unformatted' )
            write(169) xpt
            write(176) ypt
            write(177) zpt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='xrgu.dat'   ,form='unformatted' )
            open ( unit=176, file='xthgu.dat'  ,form='unformatted' )
            open ( unit=177, file='xphigu.dat' ,form='unformatted' )
            write(169) xrgt
            write(176) xthgt
            write(177) xphigt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='bmstu.dat'   ,form='unformatted' )
            write(169) bmst
            close(169)

            open ( unit=169, file='vpsnxu.dat'   ,form='unformatted' )
            open ( unit=176, file='vpsnyu.dat'   ,form='unformatted' )
            open ( unit=177, file='vpsnzu.dat'   ,form='unformatted' )
            write(169) vpsnxt
            write(176) vpsnyt
            write(177) vpsnzt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='vhsnxu.dat'   ,form='unformatted' )
            open ( unit=176, file='vhsnyu.dat'   ,form='unformatted' )
            open ( unit=177, file='vhsnzu.dat'   ,form='unformatted' )
            write(169) vhsnxt
            write(176) vhsnyt
            write(177) vhsnzt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='bdirsxu.dat'   ,form='unformatted' )
            open ( unit=176, file='bdirsyu.dat'   ,form='unformatted' )
            open ( unit=177, file='bdirszu.dat'   ,form='unformatted' )
            write(169) bdirsxt
            write(176) bdirsyt
            write(177) bdirszt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='bdirpxu.dat'   ,form='unformatted' )
            open ( unit=176, file='bdirpyu.dat'   ,form='unformatted' )
            open ( unit=177, file='bdirpzu.dat'   ,form='unformatted' )
            write(169) bdirpxt
            write(176) bdirpyt
            write(177) bdirpzt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='bdirhxu.dat'   ,form='unformatted' )
            open ( unit=176, file='bdirhyu.dat'   ,form='unformatted' )
            open ( unit=177, file='bdirhzu.dat'   ,form='unformatted' )
            write(169) bdirhxt
            write(176) bdirhyt
            write(177) bdirhzt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='ehpxu.dat'   ,form='unformatted' )
            open ( unit=176, file='ehpyu.dat'   ,form='unformatted' )
            open ( unit=177, file='ehpzu.dat'   ,form='unformatted' )
            write(169) ehpxt
            write(176) ehpyt
            write(177) ehpzt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='ephxu.dat'   ,form='unformatted' )
            open ( unit=176, file='ephyu.dat'   ,form='unformatted' )
            open ( unit=177, file='ephzu.dat'   ,form='unformatted' )
            write(169) ephxt
            write(176) ephyt
            write(177) ephzt
            close(169)
            close(176)
            close(177)

         print *,'ehp ',maxval(ehpxt),minval(ehpxt)
         print *,'ehp ',maxval(ehpyt),minval(ehpyt)
         print *,'ehp ',maxval(ehpzt),minval(ehpzt)

            open ( unit=169, file='delhpu.dat'   ,form='unformatted' )
            write(169) delhpt
            close(169)

            open ( unit=169, file='delphu.dat'   ,form='unformatted' )
            write(169) delpht
            close(169)

            open ( unit=169, file='gsthetaxu.dat' ,form='unformatted' )
            open ( unit=176, file='gsthetayu.dat' ,form='unformatted' )
            open ( unit=177, file='gsthetazu.dat' ,form='unformatted' )
            write(169) gsthetaxt
            write(176) gsthetayt
            write(177) gsthetazt
            close(169)
            close(176)
            close(177)

            open ( unit=169, file='gsphixu.dat'   ,form='unformatted' )
            open ( unit=176, file='gsphiyu.dat'   ,form='unformatted' )
            open ( unit=177, file='gsphizu.dat'   ,form='unformatted' )
            write(169) gsphixt
            write(176) gsphiyt
            write(177) gsphizt
            close(169)
            close(176)
            close(177)


            open ( unit=169, file='gsrxu.dat'   ,form='unformatted' )
            open ( unit=176, file='gsryu.dat'   ,form='unformatted' )
            open ( unit=177, file='gsrzu.dat'   ,form='unformatted' )
            write(169) gsrxt
            write(176) gsryt
            write(177) gsrzt
            close(169)
            close(176)
            close(177)


    endif

    100 format (1x,1p10e16.6)

! The rest of the initialization is done also for everthing but the master

    if(taskid > 0) then

    ! MS: chicrit is the zenith angle below which the Sun is visible.
    ! For points on the surface this is just pi/2, but at higher
    ! altitudes it is bigger.

        do k = 1,nl
            do j = 1,nf
                do i = 1,nz
                    coschicrit(i,j,k) = cos(pie - &
                    asin( 1./ (1. + alts(i,j,k)/re) ))
                enddo
            enddo
        enddo


    ! put deni on mesh via linear interpolation
    ! and put on lower limit

    ! initialize all ions

        j0 = 1
        do n = 1,nion
            do k = 1,nl
                do j = 1,nf
                    do i = 1,nz
                        jj = 1
                        do while (  alts(i,j,k) >= zi(jj) .AND. jj <= 28 )
                            j0 = jj
                            jj = jj + 1
                        enddo
                        if ( n == 1 ) nn = pthp
                        if ( n == 2 ) nn = pthep
                        if ( n == 3 ) nn = ptnp
                        if ( n == 4 ) nn = ptop
                        if ( n == 5 ) nn = ptn2p
                        if ( n == 6 ) nn = ptnop
                        if ( n == 7 ) nn = pto2p
                        slope   = ( denii(j0+1,n) - denii(j0,n) ) &
                                / ( zi   (j0+1)   - zi   (j0) )
                        deni(i,j,k,nn) = denii(j0,n) + &
                                       ( alts(i,j,k) - zi(j0) ) * slope
                        deni(i,j,k,nn) = max ( 3. * deni(i,j,k,nn) , denmin )  ! fix #9
                        if ( alts(i,j,k) > zi(29) ) then
                            if ( n == 1 )  then
                                nn = pthp
                                deni(i,j,k,nn) = &
                                   max(denii(29,n)*zi(29)/alts(i,j,k),denmin)  ! fix #9
                            else
                                deni(i,j,k,nn) = denmin
                            endif
                        endif
                    enddo
                enddo
            enddo
        enddo

    !     initialize helium density = 10% hydrogen density

        do k = 1,nl
            do j = 1,nf
                do i = 1,nz
                    deni(i,j,k,pthep) = 0.1 * deni(i,j,k,pthp)
                enddo
            enddo
        enddo

    ! print *,'done initializing deni',taskid

    ! initialize neutrals
    ! neutral density, temperature, and neutral wind

        if ( .NOT. restart ) then
            do nll = 1,nl
                call neutambt (hrinit,nll,1)
            enddo
        endif

    ! electron and ion temperature initialization

        do k = nion1,nion2
            do n = 1,nl
                do j = 1,nf
                    do i = 1,nz
                        ti(i,j,n,k)    = tni(i,j,n)
                    enddo
                enddo
            enddo
        enddo

        do n = 1,nl
            do j = 1,nf
                do i = 1,nz
                    te(i,j,n)      = tni(i,j,n)
                enddo
            enddo
        enddo

    !       average magnetic pole grid values (deni,Ti,Te)

        j0  = nf

        do ni = nion1,nion2
            do i = 1,nz
                deni_mnp0 = 0.
                ti_mnp0   = 0.
                do k = 2,nl-1
                    deni_mnp0      = deni_mnp0 + deni(i,j0,k,ni)
                    ti_mnp0        = ti_mnp0 + ti(i,j0,k,ni)
                enddo
                deni_mnp(i,ni) = deni_mnp0 / real(nl-2)  ! fix #9
                ti_mnp(i,ni)   = ti_mnp0 / real(nl-2)  ! fix #9
            enddo
        enddo

        do i = 1,nz
            te_mnp0 = 0.
            do k = 2,nl-1
                te_mnp0     = te_mnp0 + te(i,j0,k)
            enddo
            te_mnp(i)   = te_mnp0 / real(nl-2)  ! fix #9
        enddo



    ! initialize ion velocity to zero

        do ni = nion1,nion2
            do k = 1,nl
                do j = 1,nf
                    do i = 1,nz
                        vsi(i,j,k,ni)     = 0.
                        sumvsi(i,j,k,ni)  = 0.
                    enddo
                enddo
            enddo
        enddo

    endif

! endif for taskid > 0 initialization

! read in photoabsorption rates

!   EUVAC

    if(taskid == 0) then

        do i = 1,linesuv_euvac
          read (50,105) (sigabsdt_euvac(i,j), j=1,3)
        enddo
        do j = 1,3
          do i = 1,linesuv_euvac
            sigabsdt_euvac(i,j) = tm18 * sigabsdt_euvac(i,j)
          enddo
        enddo

    endif

         
105 format (3f7.2)

      call mpi_bcast(sigabsdt_euvac,linesuv_euvac*3, &
                     MPI_REAL,0,MPI_COMM_WORLD,ierr)
    
!    print *,'end bcast 3'

! read in daytime photoionization line data
! (only n, o, he, n_2, o_2)

    if(taskid == 0) then

! EUVAC

          do j = 1,nneut
            do i = 1,linesuv_euvac
              sigidt_euvac(i,j)  = 0.
            enddo
          enddo

          do i = 1,linesuv_euvac
              read(60,106) (phionrdt_euvac(i,j), j=1,5)
                sigidt_euvac(i,ptn ) = phionrdt_euvac(i,1)
                sigidt_euvac(i,pto ) = phionrdt_euvac(i,2)
                  ! can increase He+ photoproduction rate
                  ! bailey and sellek used 2.5
                  ! JK used 1.5
                sigidt_euvac(i,pthe) = phionrdt_euvac(i,3)
                sigidt_euvac(i,ptn2) = phionrdt_euvac(i,4)
                sigidt_euvac(i,pto2) = phionrdt_euvac(i,5)
          enddo
                 
          do j = 1,nneut
            do i = 1,linesuv_euvac
                sigidt_euvac(i,j) = tm18 * sigidt_euvac(i,j)
            enddo
          enddo

106     format(5f7.2)
                 
    ! read in nighttime photoionization line data
    ! (only o, n_2, n0, o_2)

   ! initialize nighttime photoionization rates to zero

        do j = 1,nneut
          do i = 1,linesnt
              sigint(i,j)  = 0.
          enddo
        enddo

        do i = 1,linesnt
          read(61,106) (phionrnt(i,j), j=1,4)
            sigint(i,pto ) = phionrnt(i,1)
            sigint(i,ptn2) = phionrnt(i,2)
            sigint(i,ptno) = phionrnt(i,3)
            sigint(i,pto2) = phionrnt(i,4)
        enddo
                 
        do j = 1,nneut
            do i = 1,linesnt
                sigint(i,j) = tm18 * sigint(i,j)
            enddo
        enddo

    endif

    call mpi_bcast(sigidt_euvac,linesuv_euvac*nneut, &
                     MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(sigint,linesnt*nneut, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)

!         print *,'end bcast 4'

! read in f74113, ai data and set euv flux
! (from richards et al., jgr 99, 8981, 1994)

    if( taskid == 0 ) then
        do i = 1,linesuv_euvac
            read (65,107) (fluxdat_euvac(i,j),j=1,2)
        enddo
107     format (f6.3,1pe11.4)
    endif

    call mpi_bcast(fluxdat_euvac,2*linesuv_euvac, &
                       MPI_REAL,0,MPI_COMM_WORLD,ierr)

! read in angles for nighttime deposition fluxes

    if(taskid == 0) then
        do i = 1,linesnt
            read(66,108) (thetant(i,j), j=1,4)
        enddo
    endif
    108 format (4f7.1)
    call mpi_bcast(thetant,linesnt*4, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)


! read in min/max altitude for nighttime deposition fluxes
!   zaltnt(i,1): zmin(i)
!   zaltnt(i,2): zmax(i)

    if(taskid == 0) then
        do i = 1,linesnt
            read(67,108) (zaltnt(i,j), j=1,2)
        enddo
    endif

    call mpi_bcast(zaltnt,linesnt*2, &
                   MPI_REAL,0,MPI_COMM_WORLD,ierr)

! Do this for everything but the master

    if (taskid > 0) then

    ! call nighttime euv flux subroutines
    ! (lyman beta 1026, he i 584, he ii 304, lyman alpha 1216)

        do k = 1,nl
            do j = 1,nf
                call sf1026 ( f1026,1,j,k )
                call sf584  ( f584 ,2,j,k )
                call sf304  ( f304 ,3,j,k )
                call sf1216 ( f1216,4,j,k )
                do l = 1,91
                    do i = 1,nz
                        fluxnt(i,j,k,l,1) = f1026(i,j,k,l)
                        fluxnt(i,j,k,l,2) = f584 (i,j,k,l)
                        fluxnt(i,j,k,l,3) = f304 (i,j,k,l)
                        fluxnt(i,j,k,l,4) = f1216(i,j,k,l)
                    enddo
                enddo
            enddo
        enddo

    ! initialize e x b drift to 0

        do k = 1,nl
            do j = 1,nf
                do i = 1,nzp1
                    vexbs(i,j,k) = 0.
                enddo
            enddo
        enddo

        do k = 1,nl
            do j = 1,nfp1
                do i = 1,nz
                    vexbp(i,j,k) = 0.
                enddo
            enddo
        enddo

        do k = 1,nlp1
            do j = 1,nf
                do i = 1,nz
                    vexbh(i,j,k) = 0.
                enddo
            enddo
        enddo

    ! intialize diagnostic variables to 0

        do k = 1,nl
            do j = 1,nf
                do i = 1,nz
                    u1p(i,j,k) = 0.
                    u2s(i,j,k) = 0.
                    u3h(i,j,k) = 0.
                    u1(i,j,k) = 0.
                    u2(i,j,k) = 0.
                    u3(i,j,k) = 0.
                    u4(i,j,k) = 0.
                    u5(i,j,k) = 0.
                    gp(i,j,k) = 0.
!!$                    xnuin0(i,j,k) = 0.
!!$                    xnuin1(i,j,k) = 0.
!!$                    xnuin2(i,j,k) = 0.
                enddo
            enddo
        enddo

!        do k = 1,nion
!            do n = 1,nl
!                do j = 1,nf
!                    do i = 1,nz
!                        t1(i,j,n,k) = 0.
!                        t2(i,j,n,k) = 0.
!                        t3(i,j,n,k) = 0.
!                    enddo
!                enddo
!            enddo
!        enddo
    endif

!!$    ! -----------------------------------------------
!!$    ! Read in the chapman function table lookup
!!$    ! -----------------------------------------------
!!$    
!!$    open(unit=77, file='./chapmanfnc.inp' , access='stream' )
!!$    read(77) cfx
!!$    read(77) cfchi
!!$    read(77) chapfunc
!!$    close(77)
!!$

    if ( taskid == 0 ) then
        deallocate (altstmp,glatstmp,glonstmp)
        deallocate (altstmp1,glatstmp1,glonstmp1)
        deallocate (altst,glatst,glonst)
        deallocate (baltst,blatst,blonst)
        deallocate (xst,yst,zst)
        deallocate (altptmp,blatptmp,blonptmp)
        deallocate (baltpt)
        deallocate (vpsnxt,vpsnyt,vpsnzt)
        deallocate (vhsnxt,vhsnyt,vhsnzt)
        deallocate (xpt,ypt,zpt)
        deallocate (bdirsxt,bdirsyt,bdirszt)
        deallocate (gsthetaxt,gsthetayt,gsthetazt)
        deallocate (gsphixt,gsphiyt,gsphizt)
        deallocate (gsrxt,gsryt,gsrzt)
        deallocate (xrgt,xthgt,xphigt)
        deallocate (bmst)
    endif

! set up weimer by master: weimer grid / data

    if ( taskid == 0 .AND. lweimer) then
      call weimer_omni(imonth_1,iday_1,&
           imonth_2,iday_2,iyear)  ! generates 'weimer_data.inp'
      call weimer_grid  ! generates 'weimer_grid.dat'
      call test_w05sc   ! generates 'phi_weimer.inp'
    endif

    print *,' finished initialization taskid = ',taskid

    return
    end subroutine initial


!******************************************
!******************************************

!            transprt

!******************************************
!******************************************

    subroutine transprt (nfl,nll)

    use parameter_mod
    use parameter_apex_mod
    use message_passing_mod
    use variable_mod
    use namelist_mod
    use chemistry_mod
    use time_mod
    use atomic_mod
    use conductance_mod
    use grid_mod
    use hardy_mod
    use misc_mod

    implicit none  ! fix #1

    ! Arguments
    integer, intent(in) :: nfl, nll

    real :: prod(nz,nion),loss(nz,nion),lossr, &
            phprodr(nz,nion),chrate(nz,nchem), &
            chloss(nz,nion),chprod(nz,nion),relossr(nz,nion)
    real :: deni_old(nz,nion),te_old(nz),ti_old(nz,nion),vsi_old(nz,nion)
    real :: tvn(nz,nl)
    real :: nuin(nz,nion,nneut),nuij(nz,nion,nion),sumnuj(nz,nion)
    real :: vsin(nz,nion),vsidn(nz,nion),denin(nz,nion),cs(nz,nion)
    real :: ten(nz),tin(nz,nion)

    real :: vnq_apex, vnp_apex, vnphi_apex ! (GL)
    real :: denntotn, denntots
    real :: pnout(nz,nion), psout(nz,nion)

    ! fix #1: previously implicit locals in transprt
    integer :: i, j, n, ni
    real    :: gs_apex, gp_apex, gphi_apex
    real    :: fzero, clat, slat, cfac, tn0

! calculation of production and loss
!   phprodr: photo production rates
!   chrate:  chemical rates (ichem)
!   chloss:  chemical loss term
!   chprod:  chemical production term
!   relossr: recombination loss rates
!   hardy: auroral ionization rate

! initialize tvn and gs and cfs (centrifugal force)

    do i = 1,nz
        tvn(i,nll) = 0.
        gs(i,nll)  = 0.
        cfs(i,nll)  = 0.
    enddo

    do i = 1,nz
        ne(i,nfl,nll)   = 1.
        te_old(i)       = te(i,nfl,nll)
        do j = nion1,nion2
            deni_old(i,j) = deni(i,nfl,nll,j)
            ne(i,nfl,nll) = ne(i,nfl,nll) + deni(i,nfl,nll,j)
            ti_old(i,j)   = ti(i,nfl,nll,j)
            vsi_old(i,j)  = vsi(i,nfl,nll,j)
        enddo
    enddo

    call photprod ( phprodr,nfl,nll               ) ! calculates phprodr
    call chemrate ( chrate,nfl,nll                ) ! calculates chrate
    call chempl   ( chrate,chloss,chprod,nfl,nll  ) ! calcualtes chloss,chprod
    call recorate ( relossr,nfl,nll               ) ! calculates relossr

!   ionization rates for N2, O2, O
!   from Rees (1989) p. 44

!   N2,O2,0

!   initiialize precipitaing quantities to 0

    denntotn    = 0. 
    denntots    = 0. 

    do i = 1,nz
      do j = nion1,nion2
        pnout(i,j)  = 0.
        psout(i,j)  = 0.
      enddo
    enddo

    do i = 1,nz
!      u5(i,nfl,nll)  = tpn(i,nfl,nll)
      denntotn = 0.92 * denn(i,nfl,nll,ptn2) + &
                 1.00 * denn(i,nfl,nll,pto2) + &
                 0.56 * denn(i,nfl,nll,pto)
      denntots = 0.92 * denn(i,nfl,nll,ptn2) + &
                 1.00 * denn(i,nfl,nll,pto2) + &
                 0.56 * denn(i,nfl,nll,pto)
      pnout(i,ptn2) = tpn(i,nfl,nll) * &
                      0.92 * denn(i,nfl,nll,ptn2) / denntotn
      psout(i,ptn2) = tps(i,nfl,nll) * &
                      0.92 * denn(i,nfl,nll,ptn2) / denntots
      pnout(i,pto2) = tpn(i,nfl,nll) * &
                      1.00 * denn(i,nfl,nll,pto2) / denntotn
      psout(i,pto2) = tps(i,nfl,nll) * &
                      1.00 * denn(i,nfl,nll,pto2) / denntots
      pnout(i,pto)  = tpn(i,nfl,nll) * &
                      0.56 * denn(i,nfl,nll,pto)  / denntotn
      psout(i,pto)  = tps(i,nfl,nll) * &
                      0.56 * denn(i,nfl,nll,pto)  / denntots

    enddo

    do i = 1,nz
        do j = nion1,nion2
          prod  (i,j) =  phprodr(i,j) * denn(i,nfl,nll,j) &
                         + chprod(i,j) +  (psout(i,j) + pnout(i,j)) * xhardy
          lossr       =  relossr(i,j) * deni(i,nfl,nll,j) * &
                         ne(i,nfl,nll) &
                         + chloss(i,j)
          loss (i,j)  =  lossr / deni(i,nfl,nll,j)
!          loss_tot(i,nfl,nll,j) = loss(i,j)
          loss_tot(i,nfl,nll,j) = relossr(i,j) * ne(i,nfl,nll)
        enddo


!  add loss for NO+ below 90 km
!  NEED TO IMPROVE THIS

    if ( alts(i,nfl,nll) .lt. 90. ) then
      loss(i,ptnop) = loss(i,pto2p)
      loss(i,pthp)  = loss(i,pto2p)
    endif

    !     loss term for hydrogen and helium

        if ( alts(i,nfl,nll) > pcrit*re ) then
            loss(i,pthp)  = loss(i,pthp)  + 1./decay_time
            loss(i,pthep) = loss(i,pthep) + 1./decay_time
        ! add loss term for O+ (JK)
        !          loss(i,ptop) = loss(i,ptop) + 1./decay_time
        endif


    ! K     centrifugal force (see notes 2012/01/04)

        fzero = 3.369
        clat = cos(pie*glats(i,nfl,nll)/180.0)
        slat = sin(pie*glats(i,nfl,nll)/180.0)
        cfs(i,nll)   =  -fzero * &
        (clat*xrg(i,nfl,nll) + slat*xthg(i,nfl,nll)) &
        * (re + alts(i,nfl,nll)) * clat / re

!       gp(i,nfl,nll)   = gzero &
!       * ( re / (re + alts(i,nfl,nll)) ) ** 2 &
!       * ( gsrx(i,nfl,nll)*vpsnx(i,nfl,nll) + &
!       gsry(i,nfl,nll)*vpsny(i,nfl,nll) + &
!       gsrz(i,nfl,nll)*vpsnz(i,nfl,nll)  )
!
!       vnq(i,nfl,nll) = v(i,nfl,nll) * &
!       ( gsthetax(i,nfl,nll) * bdirsx(i,nfl,nll) + &
!       gsthetay(i,nfl,nll) * bdirsy(i,nfl,nll) + &
!       gsthetaz(i,nfl,nll) * bdirsz(i,nfl,nll)   ) + &
!       u(i,nfl,nll) * &
!       ( gsphix(i,nfl,nll) * bdirsx(i,nfl,nll) + &
!       gsphiy(i,nfl,nll) * bdirsy(i,nfl,nll) + &
!       gsphiz(i,nfl,nll) * bdirsz(i,nfl,nll)   )   + &
!       w(i,nfl,nll) * &
!       ( gsrx(i,nfl,nll) * bdirsx(i,nfl,nll) + &
!       gsry(i,nfl,nll) * bdirsy(i,nfl,nll) + &
!       gsrz(i,nfl,nll) * bdirsz(i,nfl,nll)   )
!
!       vnp(i,nfl,nll) = v(i,nfl,nll) * &
!       ( gsthetax(i,nfl,nll) * vpsnx(i,nfl,nll) + &
!       gsthetay(i,nfl,nll) * vpsny(i,nfl,nll) + &
!       gsthetaz(i,nfl,nll) * vpsnz(i,nfl,nll)   ) + &
!       u(i,nfl,nll) * &
!       ( gsphix(i,nfl,nll) * vpsnx(i,nfl,nll) + &
!       gsphiy(i,nfl,nll) * vpsny(i,nfl,nll) + &
!       gsphiz(i,nfl,nll) * vpsnz(i,nfl,nll)   )   + &
!       w(i,nfl,nll) * &
!       ( gsrx(i,nfl,nll) * vpsnx(i,nfl,nll) + &
!       gsry(i,nfl,nll) * vpsny(i,nfl,nll) + &
!       gsrz(i,nfl,nll) * vpsnz(i,nfl,nll)   )
!
!       vnphi(i,nfl,nll) = v(i,nfl,nll) * &
!       ( gsthetax(i,nfl,nll) * vhsnx(i,nfl,nll) + &
!       gsthetay(i,nfl,nll) * vhsny(i,nfl,nll) + &
!       gsthetaz(i,nfl,nll) * vhsnz(i,nfl,nll)   ) + &
!       u(i,nfl,nll) * &
!       ( gsphix(i,nfl,nll) * vhsnx(i,nfl,nll) + &
!       gsphiy(i,nfl,nll) * vhsny(i,nfl,nll) + &
!       gsphiz(i,nfl,nll) * vhsnz(i,nfl,nll)   )   + &
!       w(i,nfl,nll) * &
!       ( gsrx(i,nfl,nll) * vhsnx(i,nfl,nll) + &
!       gsry(i,nfl,nll) * vhsny(i,nfl,nll) + &
!       gsrz(i,nfl,nll) * vhsnz(i,nfl,nll)   )

!c$$$! Project gravity g (downward) using APEX (GL)
!c$$$

! GL

        gs_apex   = -dvec(i,nfl,nll,3,3) * &
                     gzero * ( re / (re + alts(i,nfl,nll)) ) ** 2
        gp_apex   = -dvec(i,nfl,nll,2,3) * &
                     gzero * ( re / (re + alts(i,nfl,nll)) ) ** 2
        gphi_apex =  dvec(i,nfl,nll,1,3) * &
                     gzero * ( re / (re + alts(i,nfl,nll)) ) ** 2

        gs(i,nll)     = gs_apex
        gp(i,nfl,nll) = gp_apex
        gphi(i,nll)   = gphi_apex

! GL

! Projection using APEX parameters

       vnq_apex   = dvec(i,nfl,nll,3,1)*u(i,nfl,nll) + &
                    dvec(i,nfl,nll,3,2)*v(i,nfl,nll) + &
                    dvec(i,nfl,nll,3,3)*w(i,nfl,nll)
       vnp_apex   = dvec(i,nfl,nll,2,1)*u(i,nfl,nll) + &
                    dvec(i,nfl,nll,2,2)*v(i,nfl,nll) + &
                    dvec(i,nfl,nll,2,3)*w(i,nfl,nll)
       vnphi_apex = dvec(i,nfl,nll,1,1)*u(i,nfl,nll) + &
                    dvec(i,nfl,nll,1,2)*v(i,nfl,nll) + &
                    dvec(i,nfl,nll,1,3)*w(i,nfl,nll)

       vnq(i,nfl,nll) = vnq_apex

! Reverse sign for vnp since it is inward/downward in the APEX coordinates

        vnp(i,nfl,nll)   = -vnp_apex
        vnphi(i,nfl,nll) = vnphi_apex

        tvn(i,nll)    = vnq(i,nfl,nll)

    enddo

    call update ( tvn,nuin,sumnuj,nuij,nfl,nll )

! define new arrays for velocity and density

    do ni = nion1,nion2
        do i = 1,nz
            vsin (i,ni) = vsi(i,nfl,nll,ni)
            vsidn(i,ni) = vsid(i,nfl,nll,ni)
            denin(i,ni) = deni(i,nfl,nll,ni)
        enddo
    enddo

! define sound velocity used in vsisolv

    do ni = nion1,nion2
        do i = 1,nz
            cfac     = 1.6667 * 8.6174e-5 * te(i,nfl,nll) / ami(ni)
            cs(i,ni) = 9.79e5 * sqrt(cfac)
        enddo
    enddo

! update variables

    do ni = nion1,nion2

        call vsisolv ( vsin(1,ni),vsidn(1,ni),vsi_old(1,ni), &
                       sumnuj(1,ni),nfl,nll,cs(1,ni) )

    ! compensating filter

        call smoothz ( vsin(1,ni), 1 )

    ! put stuff back into velocity array

        do i = 1,nz
            vsi(i,nfl,nll,ni)  = vsin(i,ni)
            vsid(i,nfl,nll,ni) = vsidn(i,ni)
        enddo

        call densolv2 ( ni,denin(1,ni), &
                        prod(1,ni),loss(1,ni),deni_old(1,ni),nfl,nll )

    ! put stuff back into density array

        do i = 1,nz
            deni(i,nfl,nll,ni) = denin(i,ni)
        enddo

    ! put floor on density

        do i = 1,nz
            deni(i,nfl,nll,ni) = max ( deni(i,nfl,nll,ni), denmin )  ! fix #9
        ! below commented out (JK)
            if ( alts(i,nfl,nll) > pcrit*re .AND. ni == pthp ) &
              deni(i,nfl,nll,ni) = max ( deni(i,nfl,nll,ni), .1 )   ! fix #9
            if ( alts(i,nfl,nll) > pcrit*re .AND. ni == pthep ) &
              deni(i,nfl,nll,ni) = max ( deni(i,nfl,nll,ni), .01 )  ! fix #9
        enddo

    enddo


! define new arrays for temperature

    do ni = nion1,nion2
        do i = 1,nz
            tin(i,ni)  = ti(i,nfl,nll,ni)
        enddo
    enddo

    do i = 1,nz
        ten(i)  = te(i,nfl,nll)
    enddo

! temperatures (with floors and warnings)

    tn0 = 200. ! floor on temperature

    call etemp  (ten,te_old,phprodr,nfl,nll)
    do i = 1,nz
        ! fix #4: second floor must use te(i,nfl,nll) not ten(i), otherwise
        !         the tn floor from line 1 is overwritten by the tn0 floor from ten
        ! fix #9: amax1/amin1 -> max/min
        te(i,nfl,nll)  = max(tn(i,nfl,nll), ten(i))
        te(i,nfl,nll)  = max(tn0, te(i,nfl,nll))
        te(i,nfl,nll)  = min(te(i,nfl,nll), tmax)
        if ( te(i,nfl,nll) < 0 ) then
            print *,' T(e) negative: i,nfl,nll taskid',i,nfl,nll,taskid
            stop
        endif
    enddo

    call htemp  (tin(1,pthp) ,ti_old(1,pthp) ,tvn,nuin,nfl,nll)
    do i = 1,nz
        ! fix #4 + #9
        ti(i,nfl,nll,pthp)  = max(tn(i,nfl,nll), tin(i,pthp))
        ti(i,nfl,nll,pthp)  = max(tn0, ti(i,nfl,nll,pthp))
        ti(i,nfl,nll,pthp)  = min(ti(i,nfl,nll,pthp), tmax)
        if ( ti(i,nfl,nll,pthp) < 0 ) then
            print *,' T(H) negative: i,nfl,nll',i,nfl,nll
            stop
        endif
    enddo

    call hetemp (tin(1,pthep),ti_old(1,pthep),tvn,nuin,nfl,nll)
    do i = 1,nz
        ! fix #4 + #9
        ti(i,nfl,nll,pthep)  = max(tn(i,nfl,nll), tin(i,pthep))
        ti(i,nfl,nll,pthep)  = max(tn0, ti(i,nfl,nll,pthep))
        ti(i,nfl,nll,pthep)  = min(ti(i,nfl,nll,pthep), tmax)
        if ( ti(i,nfl,nll,pthep) < 0 ) then
            print *,' T(He) negative: i,nfl,nll',i,nfl,nll
            stop
        endif
    enddo

    call otemp  (tin(1,ptop) ,ti_old(1,ptop) ,tvn,nuin,nfl,nll)
    do i = 1,nz
        ! fix #4 + #9
        ti(i,nfl,nll,ptop)  = max(tn(i,nfl,nll), tin(i,ptop))
        ti(i,nfl,nll,ptop)  = max(tn0, ti(i,nfl,nll,ptop))
        ti(i,nfl,nll,ptop)  = min(ti(i,nfl,nll,ptop), tmax)
        if ( ti(i,nfl,nll,ptop) < 0 ) then
            print *,' T(O) negative: i,nfl,nll',i,nfl,nll
            stop
        endif
    enddo

    do i = 1,nz
        ti(i,nfl,nll,ptnp )    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,ptn2p)    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,ptnop)    = ti(i,nfl,nll,ptop)
        ti(i,nfl,nll,pto2p)    = ti(i,nfl,nll,ptop)
    enddo

    return
    end subroutine transprt




!******************************************
!******************************************

!            update

!******************************************
!******************************************

    subroutine update ( tvn,nuin,sumnuj,nuij,nfl,nll )

    use parameter_mod
    use message_passing_mod
    use variable_mod
    use namelist_mod
    use atomic_mod
    use conductance_mod
    use grid_mod
    use exb_mod

    implicit none  ! fix #1

    ! Arguments
    integer, intent(in) :: nfl, nll
    real, intent(in)    :: tvn(nz,nl)
    real, intent(out)   :: nuin(nz,nion,nneut), nuij(nz,nion,nion)
    real, intent(out)   :: sumnuj(nz,nion)

    real :: nuint(nz,nion)
    real :: nufacij
    real :: k0, mi
    real :: nuen(nz,nf,nl), nuei(nz,nf,nl)

    ! fix #1: previously implicit locals in update
    integer :: i, iz, nn, ni, nj, nni
    real    :: teff, fac, tfactor
    real    :: alame1, alame2, alame, alam, amufac
    real    :: pip, pim, denid, pep, pem, dened
    real    :: term1, term2, term3
    real    :: Te_ev, dene, oce, sige, cole, denome, sigpe, sighe
    real    :: sigpi, sighi, sighic, sigpic, sigpic0
    real    :: oci, sigi, coli, denomi
    real    :: ang, bang, del, b, fdpv, fdpg, fdphiv, fdphig

! ion-neutral collision frequency

! initialize everything to 0

    do nn = 1,nneut
        do ni = nion1,nion2
            do iz = 1,nz
                nuin (iz,ni,nn) = 0.
                nuint(iz,ni)    = 0.
            enddo
        enddo
    enddo

! collision frequencies/factors

! hydrogen (H)

    ni = pthp
    do nn = 1,nneut
        do i = 1,nz
            if ( nn == pto ) then
                teff    = ti(i,nfl,nll,ni)
                fac     = ( 1.00 - .047 * log10(teff) ) ** 2
                tfactor = sqrt(teff) * fac
                nuin(i,ni,nn)  = 6.61e-11 * denn(i,nfl,nll,nn) * tfactor
            else
                nuin(i,ni,nn) = nufacin(ni,nn) * denn(i,nfl,nll,nn)
            endif
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! helium (He)

    ni = pthep
    do nn = 1,nneut
        do i = 1,nz
            nuin(i,ni,nn) = nufacin(ni,nn) * denn(i,nfl,nll,nn)
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! nitrogen (N)

    ni = ptnp
    do nn = 1,nneut
        do i = 1,nz
            nuin(i,ni,nn) = nufacin(ni,nn) * denn(i,nfl,nll,nn)
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! oxygen (O)

    ni = ptop
    do nn = 1,nneut
        do i = 1,nz
            if ( nn == pto ) then
                teff    = 0.5 * ( ti(i,nfl,nll,ni) + tn(i,nfl,nll) )
                fac     = ( 1.04 - .067 * log10(teff) ) ** 2
                tfactor = sqrt(teff) * fac
                nuin(i,ni,nn)  = 4.45e-11 * denn(i,nfl,nll,nn) * tfactor
            else
                nuin(i,ni,nn) = nufacin(ni,nn) * denn(i,nfl,nll,nn)
            endif
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! nitrogen 2(N2)

    ni = ptn2p
    do nn = 1,nneut
        do i = 1,nz
            if ( nn == ptn2 ) then
                teff    = 0.5 * ( ti(i,nfl,nll,ni) + tn(i,nfl,nll) )
                fac     = ( 1.00 - .069 * log10(teff) ) ** 2
                tfactor = sqrt(teff) * fac
                nuin(i,ni,nn) = 5.14e-11 * denn(i,nfl,nll,nn) * tfactor
            else
                nuin(i,ni,nn) = nufacin(ni,nn) * denn(i,nfl,nll,nn)
            endif
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! nitrous oxide (N0)

    ni = ptnop
    do nn = 1,nneut
        do i = 1,nz
            nuin(i,ni,nn) = nufacin(ni,nn) * denn(i,nfl,nll,nn)
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! oxygen 2(O2)

    ni = pto2p
    do nn = 1,nneut
        do i = 1,nz
            if ( nn == pto2 ) then
                teff    = 0.5 * ( ti(i,nfl,nll,ni) + tn(i,nfl,nll) )
                fac     = ( 1.00 - .073 * log10(teff) ) ** 2
                tfactor = sqrt(teff) * fac
                nuin(i,ni,nn) = 2.59e-11 * denn(i,nfl,nll,nn) * tfactor
            else
                nuin(i,ni,nn) = nufacin(ni,nn) * denn(i,nfl,nll,nn)
            endif
            nuint(i,ni) = nuint(i,ni) + nuin(i,ni,nn)
        enddo
    enddo

! ion-ion collision frequency

      do i = 1,nz
        do ni = nion1,nion2
            do nj = nion1,nion2
                if(ni /= nj) then
                    alame1  = ( ami(ni) + ami(nj) ) * evtok / &
                              ( ami(ni)*ti(i,nfl,nll,nj) + &
                                ami(nj)*ti(i,nfl,nll,ni) )
                    alame2  = deni(i,nfl,nll,ni) * evtok / ti(i,nfl,nll,ni) + &
                              deni(i,nfl,nll,nj) * evtok / ti(i,nfl,nll,nj)
                    if ( alame2 < 0 ) then
                        print *,'ni,i,nj,nfl,nll,tii,tij,alame1,alame2,nii,nij', &
                        ni,i,nj,nfl,nll,ti(i,nfl,nll,ni),ti(i,nfl,nll,nj), &
                        alame1,alame2, &
                        deni(i,nfl,nll,ni),deni(i,nfl,nll,nj)
                        stop
                    endif
                    alame   = alame1 * sqrt(alame2)
                    alam    = 23. - log(alame)
                    amufac  = (ami(nj)/ami(ni))/(ami(ni) +ami(nj))
                    nufacij = 9.2e-2*alam*sqrt(amufac)
                    nuij(i,ni,nj) =  nufacij * deni(i,nfl,nll,nj) &
                                   / sqrt( ti(i,nfl,nll,ni)**3 )
                else
                    nuij(i,ni,nj) = 0.
                endif
            enddo
        enddo
    enddo

!     add this for restart (sarah mcdonald)
!     it's used in term 1 below

    do i = 1,nz
        do nni = nion1,nion2
            sumvsi(i,nfl,nll,nni) = 0.
            do nj = nion1,nion2
                sumvsi(i,nfl,nll,nni) =   sumvsi(i,nfl,nll,nni) + &
                                          nuij(i,nni,nj)*vsi(i,nfl,nll,nj)
            enddo
        enddo
    enddo

! sumnuj: sum of ion-ion coll freq and nuin

    do ni = nion1,nion2
        do i = 1,nz
            sumnuj(i,ni) = 0.
            do nj = nion1,nion2
                sumnuj(i,ni) = sumnuj(i,ni) + nuij(i,ni,nj)
            enddo
            sumnuj(i,ni) = sumnuj(i,ni) + nuint(i,ni)
            nuin_tot(i,nfl,nll,ni) = nuint(i,ni)
        enddo
    enddo

! update ne

    do i = 1,nz
        ne(i,nfl,nll) = 1.
        do ni = nion1,nion2
            ne(i,nfl,nll) = ne(i,nfl,nll) + deni(i,nfl,nll,ni)
        enddo
    enddo

! get a new value for vsid

    do i = 2,nz-1
        do ni = nion1,nion2
            mi    = amu * ami(ni)
            k0    = bolt / mi
            term1 = nuint(i,ni) * tvn(i,nll) + &
                    sumvsi(i,nfl,nll,ni) + gs(i,nll) !+ cfs(i,nll) ! need to modify for apex
            pip   = 0.5 * (   deni(i+1,nfl,nll,ni) * ti(i+1,nfl,nll,ni) &
                  + deni(i,nfl,nll,ni)   * ti(i,nfl,nll,ni)   )
            pim   = 0.5 * (   deni(i,nfl,nll,ni)   * ti(i,nfl,nll,ni) &
                  + deni(i-1,nfl,nll,ni) * ti(i-1,nfl,nll,ni) )
            denid = &
                (        deni(i-1,nfl,nll,ni) &
                  + 4. * deni(i,nfl,nll,ni) &
                  +      deni(i+1,nfl,nll,ni)  ) / 6.
            term2 =  - bms(i,nfl,nll) * k0 /  denid &
                  * ( pip - pim ) / d22s(i,nfl,nll)
            pep   = 0.5 * (   ne(i+1,nfl,nll) * te(i+1,nfl,nll) &
                  + ne(i,nfl,nll)   * te(i,nfl,nll)   )
            pem   = 0.5 * (   ne(i,nfl,nll)   * te(i,nfl,nll) &
                  + ne(i-1,nfl,nll) * te(i-1,nfl,nll) )
            dened = &
                ( ne(i-1,nfl,nll) + 4. * ne(i,nfl,nll) + ne(i+1,nfl,nll) ) / 6.
            term3 =  - bms(i,nfl,nll) * k0 /  dened &
                  * ( pep - pem ) / d22s(i,nfl,nll)

            vsid(i,nfl,nll,ni)  =  term1 + term2 + term3
        enddo
    enddo

! fix up end points for vsid

    do ni = nion1,nion2
        vsid (1,nfl,nll,ni)    = vsid (2,nfl,nll,ni)
        vsid (nz,nfl,nll,ni)   = vsid (nz-1,nfl,nll,ni)
    enddo

! calculate the electron-neutral collision frequency
! nuen = 5.4e-10*n_n*T_e^1/2 (kelley, the earth's ionosphere, p. 462)
! nuei = 2.91e-6*n_e*lambda/T_e^3/2 (nrl formulary p. 28) where Te in ev
!        take lambda = 15

    do i = 1,nz
        Te_ev           = te(i,nfl,nll) / evtok
        nuei(i,nfl,nll) = 2.91e-6*ne(i,nfl,nll)*15./ Te_ev / sqrt(Te_ev)
        nuen(i,nfl,nll) = 0
        do nn = 1,nneut
            nuen(i,nfl,nll) = nuen(i,nfl,nll) + 5.4e-10 * &
            denn(i,nfl,nll,nn) * sqrt(te(i,nfl,nll))
        enddo
    enddo


! calculate pedersen and hall conductivities

    do i = 1,nz
        dene    = ne(i,nfl,nll)
        oce     = 1.76e7 * bmag * bms(i,nfl,nll)
        sige    = dene * charge * sol / ( bmag * bms(i,nfl,nll) )
        cole    = nuen(i,nfl,nll) / oce
        denome  = 1. + cole * cole
        sigpe   = sige * cole / denome
        sighe   = sige * cole * cole / denome
        sigpi   = 0.
        sighi   = 0.
        sighic  = 0.
        sigpic  = 0.
        sigpic0 = 0.
        do ni = nion1,nion2
            oci     = 9580. * bmag * bms(i,nfl,nll) / ami(ni)
            sigi    = deni(i,nfl,nll,ni) * charge * sol / &
                      ( bmag * bms(i,nfl,nll) )
            coli    = nuint(i,ni) / oci
            denomi  = 1. + coli * coli
            sigpi   = sigpi   + sigi * coli / denomi
            sigpic  = sigpic  + sigi * coli / denomi / oci
            sighi   = sighi   + sigi * coli * coli / denomi
            sighic  = sighic  + sigi / denomi / oci
        enddo
        sigmap(i,nfl,nll)    = sigpi + sigpe
        sigmah(i,nfl,nll)    = sighi - sighe
        sigmapic(i,nfl,nll)  = sigpic
        sigmahic(i,nfl,nll)  = sighic
        sigma0(i,nfl,nll)    = 2.e12
        if (alts(i,nfl,nll) <= 1000.) &
        sigma0(i,nfl,nll)    = dene * charge * charge / &
                               emass / (nuen(i,nfl,nll)+nuei(i,nfl,nll))
        if (alts(i,nfl,nll) >= 1.e4) then
            sigmap(i,nfl,nll)   = 0.00
            sigmah(i,nfl,nll)   = 0.00
            sigmapic(i,nfl,nll) = 0.00
            sigmahic(i,nfl,nll) = 0.00
        endif
    enddo


!    print *,'Ln_inv',i,maxval(Ln_inv)

    if ( .NOT. hall ) then
        do i=1,nz
            sigmah(i,nfl,nll)   = 0.
            sigmahic(i,nfl,nll) = 0.
        enddo
    endif

    hipcp(nfl,nll)    = 0.
    hipcphi(nfl,nll)  = 0.
    hihcm(nfl,nll)    = 0.
    hidpv(nfl,nll)    = 0.
    hidphiv(nfl,nll)  = 0.
    hidpg(nfl,nll)    = 0.
    hidphig(nfl,nll)  = 0.

    hipc(nfl,nll)     = 0.
    hihc(nfl,nll)     = 0.
    hidv(nfl,nll)     = 0.

    do i = 1,nz
        ang     = .5 * pie - blats(i,nfl,nll) * pie / 180.
        bang    = blats(nz,nfl,nll) * pie / 180.
        del     = sqrt ( 1. + 3. * cos(ang) * cos(ang) )
        b       = bmag * bms(i,nfl,nll)

        hipcp(nfl,nll) = hipcp(nfl,nll) + &
        sigmap(i,nfl,nll) * del / bms(i,nfl,nll) * &
        d22s(i,nfl,nll)

        hipcphi(nfl,nll) = hipcphi(nfl,nll) + &
        sigmap(i,nfl,nll) / del / bms(i,nfl,nll) * &
        d22s(i,nfl,nll)

        hihcm(nfl,nll) = hihcm(nfl,nll) + &
        sigmah(i,nfl,nll) / bms(i,nfl,nll) * &
        d22s(i,nfl,nll)

        fdpv = (bmag/sol) * ( sigmap(i,nfl,nll) * vnphi(i,nfl,nll) + &
        sigmah(i,nfl,nll) * vnp(i,nfl,nll)     )

        hidpv(nfl,nll) = hidpv(nfl,nll) + &
        brs(i,nfl,nll) * 1.e5  * sin(ang) * &
        fdpv * d22s(i,nfl,nll)

        fdpg = (bmag/sol) * sigmapic(i,nfl,nll) * gp(i,nfl,nll)

        hidpg(nfl,nll) = hidpg(nfl,nll) + &
        brs(i,nfl,nll) * 1.e5  * sin(ang) * &
        fdpg * d22s(i,nfl,nll)

        fdphiv = (bmag/sol) * ( -sigmap(i,nfl,nll) * vnp(i,nfl,nll) + &
        sigmah(i,nfl,nll) * vnphi(i,nfl,nll)  )

        hidphiv(nfl,nll) = hidphiv(nfl,nll) + &
        re * 1.e5 * ( sin(ang) ** 3 ) / del * &
        fdphiv * d22s(i,nfl,nll)

        fdphig = (bmag/sol) * sigmahic(i,nfl,nll) * gp(i,nfl,nll)

        hidphig(nfl,nll) = hidphig(nfl,nll) + &
        re * 1.e5 * ( sin(ang) ** 3 ) / del * &
        fdphig * d22s(i,nfl,nll)


    !       integrated quantities for current

        hipc(nfl,nll) = hipc(nfl,nll) + &
        sigmap(i,nfl,nll) * del / re / sin(ang) ** 3 * &
        d22s(i,nfl,nll) / 1.e5
        hihc(nfl,nll) = hihc(nfl,nll) + &
        sigmah(i,nfl,nll) / brs(i,nfl,nll) / sin(ang) * &
        d22s(i,nfl,nll) / 1.e5
        hidv(nfl,nll) = hidv(nfl,nll) + &
        fdpv * d22s(i,nfl,nll)

    enddo

    return
    end subroutine update



