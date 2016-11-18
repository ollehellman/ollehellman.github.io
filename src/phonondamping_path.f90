
!> calculate the intensity along a path
subroutine get_intensity_along_path(bs,uc,fc,fct,fcf,qp,dr,loto,opts,mw,lsmpi)
    !> the bandstructure
    type(lo_phonon_bandstructure), intent(inout) :: bs
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(in) :: fc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-point mesh
    type(lo_monkhorst_pack_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(in) :: dr
    !> electrostatic corrections
    type(lo_loto), intent(in) :: loto
    !> all settings
    type(lo_opts), intent(in) :: opts
    !> mpi communicator
    type(lo_mpiinfo), intent(in) :: mw
    !> mpi helper
    type(lo_lsmpi), intent(in) :: lsmpi

    type(lo_phonon_selfenergy) :: se
    real(flyt), dimension(:,:,:), allocatable :: rebuf,imbuf,dumbuf
    real(flyt), dimension(:,:), allocatable :: intbuf,lwbuf,shbuf
    real(flyt), dimension(:), allocatable :: dum
    real(flyt), dimension(2) :: lsintx,lsinty
    real(flyt) :: t0,f0
    integer :: q1,nb,lqp,i,j,path,ii,jj,k

    t0=mpi_wtime()
    ! Make space for linewidths
    do q1=1,bs%nptot
        allocate(bs%p(q1)%linewidth(bs%nb))
        allocate(bs%p(q1)%shift(bs%nb))
        allocate(bs%p(q1)%threephononphasespace(bs%nb))        
    enddo
    ! Space for intensity
    allocate(bs%intensity(bs%nptot,opts%nf))
    allocate(bs%selfenergy_real(bs%nptot,opts%nf,dr%nb))
    allocate(bs%selfenergy_imag(bs%nptot,opts%nf,dr%nb))

    bs%intensity=0.0_flyt
    bs%selfenergy_real=0.0_flyt
    bs%selfenergy_imag=0.0_flyt

    ! Dump some general info
    if ( mw%talk ) then
        write(*,*) '        isotope:',opts%isotopescattering
        write(*,*) '    threephonon:',opts%thirdorder
        write(*,*) '     fourphonon:',opts%fourthorder
        write(*,*) '           loto:',opts%loto
        write(*,*) 'integrationtype:',opts%integrationtype
    endif

    ! Turn off openmp
    call omp_set_num_threads(1)


    ! Calculate self energy
    if ( mw%talk ) call lo_progressbar_init()
    allocate(rebuf(bs%nptot,opts%nf,dr%nb))
    allocate(imbuf(bs%nptot,opts%nf,dr%nb))
    allocate(lwbuf(dr%nb,bs%nptot))
    allocate(shbuf(dr%nb,bs%nptot))
    rebuf=0.0_flyt
    imbuf=0.0_flyt
    lwbuf=0.0_flyt
    shbuf=0.0_flyt
    do q1=1,lsmpi%nq
        ! global q-point index
        lqp=lsmpi%ind(q1)
        ! get the actual self-energy
        call se%generate(bs%q(lqp),bs%p(lqp),uc,fc,fct,fcf,qp,dr,loto,opts)
        ! Add it to the intensity
        do j=1,bs%nb
            do i=2,se%nf
                imbuf(lqp,i,j)=se%im_3ph(i,j)+se%im_iso(i,j)
                rebuf(lqp,i,j)=se%re_3ph(i,j)+se%re_4ph(i,j)
            enddo
           ! get the linewidth exactly at the harmonic frequency
           lwbuf(j,lqp)=lo_linear_interpolation(se%faxis,se%im_3ph(:,j)+se%im_iso(:,j),bs%p(lqp)%omega(j))*2.0_flyt
           ! get the anharmonic shift at the harmonic frequency
           shbuf(j,lqp)=lo_linear_interpolation(se%faxis,se%re_3ph(:,j)+se%re_4ph(:,j),bs%p(lqp)%omega(j))
           ! Also store the phase-space volume at this q-point
!          ! bs%p(lqp)%threephononphasespace(j)=lo_linear_interpolation(se%faxis,se%twophonondos(:,j),bs%p(lqp)%omega(j))
        enddo
        if ( mw%talk ) call lo_progressbar(' ... lineshape on path',q1,lsmpi%nq,mpi_wtime()-t0)
    enddo

    ! Add these up!
    call mpi_allreduce(rebuf,bs%selfenergy_real,bs%nptot*opts%nf*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    call mpi_allreduce(imbuf,bs%selfenergy_imag,bs%nptot*opts%nf*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    allocate(intbuf(dr%nb,bs%nptot))
    intbuf=0.0_flyt
    call mpi_allreduce(lwbuf,intbuf,bs%nptot*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    do i=1,bs%nptot
    do j=1,bs%nb
        bs%p(i)%linewidth(j)=intbuf(j,i)
    enddo
    enddo
    intbuf=0.0_flyt
    call mpi_allreduce(shbuf,intbuf,bs%nptot*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    do i=1,bs%nptot
    do j=1,bs%nb
        bs%p(i)%shift(j)=intbuf(j,i)
    enddo
    enddo
    deallocate(intbuf)
    deallocate(lwbuf)
    deallocate(shbuf)
    ! Dump the shifts and widths
    if ( mw%r .eq. 0 ) then
        call bs%write_dispersive_property(opts%enhet,'shift','outfile.dispersion_shifts',.false.)
        call bs%write_dispersive_property(opts%enhet,'linewidth','outfile.dispersion_linewidths',.false.)
!        call bs%write_dispersive_property(opts%enhet,'threephononphasespace','outfile.dispersion_threephononphasespace',.false.)
    endif

    ! dump some files, they are done now:

    rebuf=0.0_flyt
    imbuf=0.0_flyt
    t0=mpi_wtime()
    ! Figure out some neat interpolation of self-energy for really small q
    if ( mw%talk ) call lo_progressbar_init()
    do q1=1,lsmpi%nq
        lqp=lsmpi%ind(q1)
        ! what path am I on?
        path=bs%q(lqp)%path
        ! does it contain gamma?
        if ( norm2(bs%segment(path)%r1) .gt. lo_tol .and. norm2(bs%segment(path)%r2) .gt. lo_tol ) cycle
        ! seems it does, have to fix this, maybe.
        ! Good small number to use
        f0=(se%intensityaxis(2)-se%intensityaxis(1))*0.25_flyt ! smallest selfenergy
        ! Is it in the beginning or the end?
        if ( norm2(bs%segment(path)%r1) .lt. lo_tol ) then
            ! Index of gamma
            ii=(path-1)*bs%npts+1
            ! Fix the acoustic branches
            do j=1,3
                ! Find index of point that is ok
                do i=ii,ii+bs%npts-1
                    if ( bs%p(i)%omega(j) .gt. dr%omega_min*0.5_flyt ) then
                        jj=i
                        exit
                    endif
                enddo
                ! now I know that things are zero at ii, and ok at jj
                bs%selfenergy_imag(ii,:,j)=f0       ! set imaginary at gamma
                bs%selfenergy_real(ii,:,j)=0.0_flyt ! set real at gamma
                lsintx(1)=bs%q_axis(ii)-lo_sqtol
                lsintx(2)=bs%q_axis(jj)+lo_sqtol
                ! Interpolate the missing self-energies at this point
                do k=1,se%nf
                    ! y-axis for interpolation, imaginary part
                    lsinty(1)=bs%selfenergy_imag(ii,k,j)
                    lsinty(2)=bs%selfenergy_imag(jj,k,j)
                    imbuf(lqp,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(lqp))
                    ! y-axis for interpolation, real part
                    lsinty(1)=bs%selfenergy_real(ii,k,j)
                    lsinty(2)=bs%selfenergy_real(jj,k,j)
                    rebuf(lqp,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(lqp))
                enddo
            enddo
        else
            ! Same thing again, but this time gamma is at the end.
            ii=path*bs%npts
            ! loop over the three lowest branches
            do j=1,3
                jj=0
                do i=ii,(path-1)*bs%npts+1,-1
                    if ( bs%p(i)%omega(j) .gt. dr%omega_min*0.5_flyt ) then
                        jj=i
                        exit
                    endif
                enddo
                ! Interpolate this, somehow
                bs%selfenergy_imag(ii,:,j)=f0       ! set imaginary at gamma
                bs%selfenergy_real(ii,:,j)=0.0_flyt ! set real at gamma
                ! x-axis for interpolation
                lsintx(2)=bs%q_axis(ii)-lo_sqtol
                lsintx(1)=bs%q_axis(jj)+lo_sqtol
                ! interpolate to missing points
                do k=1,se%nf
                    ! y-axis for interpolation
                    lsinty(2)=bs%selfenergy_imag(ii,k,j)
                    lsinty(1)=bs%selfenergy_imag(jj,k,j)
                    imbuf(lqp,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(lqp))
                    ! y-axis for interpolation
                    lsinty(2)=bs%selfenergy_real(ii,k,j)
                    lsinty(1)=bs%selfenergy_real(jj,k,j)
                    rebuf(lqp,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(lqp))
                enddo
            enddo
        endif
        !
        if ( mw%talk ) call lo_progressbar(' ... fixing tiny q',q1,lsmpi%nq,mpi_wtime()-t0)
        !
    enddo

    ! Add this together, and add it to the self energy
    allocate(dumbuf(bs%nptot,opts%nf,dr%nb))
    dumbuf=0.0_flyt
    call mpi_allreduce(rebuf,dumbuf,bs%nptot*opts%nf*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    bs%selfenergy_real=bs%selfenergy_real+dumbuf
    dumbuf=0.0_flyt
    call mpi_allreduce(imbuf,dumbuf,bs%nptot*opts%nf*bs%nb,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    bs%selfenergy_imag=bs%selfenergy_imag+dumbuf
    deallocate(dumbuf)

    ! Now all the self-energies are nice, time to get the intensities
    allocate(intbuf(bs%nptot,opts%nf))
    allocate(dum(opts%nf))
    intbuf=0.0_flyt
    t0=mpi_wtime()
    ! Figure out some neat interpolation of self-energy for really small q
    if ( mw%talk ) call lo_progressbar_init()
    do q1=1,lsmpi%nq
        lqp=lsmpi%ind(q1)
        do j=1,bs%nb
            ! Get the lineshape
            dum=0.0_flyt
            if ( bs%p(lqp)%omega(j) .gt. lo_freqtol ) then
                call getintensity(se%faxis,bs%selfenergy_imag(lqp,:,j),bs%selfenergy_real(lqp,:,j),&
                bs%p(lqp)%omega(j),se%intensityaxis,dum)
            else
                ! acoustic branch at Gamma. Add a gaussian at 0 to no make it disappear.
                do i=1,se%nf
                    dum(i)=lo_gauss(se%intensityaxis(i),0.0_flyt,se%intensityaxis(2)-se%intensityaxis(1))
                enddo
            endif
            ! Add it to the intensity
            intbuf(lqp,:)=intbuf(lqp,:)+dum
        enddo
        !
        if ( mw%talk ) call lo_progressbar(' ... intensities',q1,lsmpi%nq,mpi_wtime()-t0)
    enddo
    ! add them up
    call mpi_allreduce(intbuf,bs%intensity,bs%nptot*opts%nf,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    deallocate(intbuf)
    deallocate(dum)

    ! Dump to file
    if ( mw%r .eq. 0 ) then
        write(*,*) 'Writing intensity to file'
        lo_allocate(bs%faxis(se%nf))
        bs%faxis=se%intensityaxis
        call bs%write_intensity(opts%enhet,logscale=.true.)
    endif
    ! And it is done!
end subroutine
