
!> Calculate the spectral function along a path in the BZ
subroutine spectralfunction_along_path(bs,uc,fc,fct,fcf,qp,dr,opts,mw)
    !> the bandstructure
    type(lo_phonon_bandstructure), intent(inout) :: bs
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(in) :: dr
    !> all settings
    type(lo_opts), intent(in) :: opts
    !> mpi communicator
    type(lo_mpi_helper), intent(inout) :: mw

    real(flyt), parameter :: timereport=30.0_flyt
    real(flyt), dimension(:,:,:), allocatable :: imbuf,rebuf
    real(flyt), dimension(:), allocatable :: seax,inax
    real(flyt) :: timer,minsmear

    ! Make some space and things like that
    init: block
        real(flyt) :: f0
        integer :: i,j,k,l
        
        timer=walltime()
        ! Figure out how much space is needed to store the buffers, and maybe print a warning
        ! in case it's some really large amount
        f0=4*dr%nb*opts%nf*bs%nptot*4.0_flyt/1024/1024
        if ( f0 .gt. 20.0_flyt .and. mw%talk ) then
            write(*,*) '... Will use at least ',tochar(int(f0)),'MB memory per rank, probably a lot more.'
        endif
        ! Make space for linewidth, shifts and so on
        do i=1,bs%nptot
            allocate(bs%p(i)%linewidth(bs%nb))
            allocate(bs%p(i)%shift3(bs%nb))
            allocate(bs%p(i)%shift4(bs%nb))
            bs%p(i)%linewidth=0.0_flyt
            bs%p(i)%shift3=0.0_flyt
            bs%p(i)%shift4=0.0_flyt
        enddo
        ! Space for intensity
        allocate(bs%intensity(bs%nptot,opts%nf))
        allocate(bs%selfenergy_real(bs%nptot,opts%nf,dr%nb))
        allocate(bs%selfenergy_imag(bs%nptot,opts%nf,dr%nb))
        allocate(rebuf(bs%nptot,opts%nf,dr%nb))
        allocate(imbuf(bs%nptot,opts%nf,dr%nb))
        allocate(bs%faxis(opts%nf))
        bs%faxis=0.0_flyt
        bs%selfenergy_real=0.0_flyt
        bs%selfenergy_imag=0.0_flyt
        bs%intensity=0.0_flyt
        lo_allocate(seax(opts%nf))
        lo_allocate(inax(opts%nf))
        seax=0.0_flyt
        inax=0.0_flyt
    end block init

    ! Get the self-energies
    selfenergy: block
        type(lo_phonon_selfenergy) :: se
        real(flyt), dimension(:), allocatable :: x,y,z
        real(flyt) :: f0,f1,t0
        integer :: i,j,k,l,ctr,path,q,nq_per_path,band

        t0=walltime()

        ! Count q-points per path, and make space for dummy arrays
        nq_per_path=0
        do q=1,bs%npts,opts%stride
            nq_per_path=nq_per_path+1
        enddo
        lo_allocate(x(nq_per_path))
        lo_allocate(y(nq_per_path))
        lo_allocate(z(nq_per_path))

        rebuf=0.0_flyt
        imbuf=0.0_flyt
        if ( mw%talk ) call lo_progressbar_init()
        do path=1,bs%npath
        do q=1,bs%npts,opts%stride
            i=(path-1)*bs%npts+q
            ! get the self-energies
            call se%generate(bs%q(i),bs%p(i),uc,fc,fct,fcf,qp,dr,opts,mw)
        
            minsmear=(se%faxis(2)-se%faxis(1))*opts%minsmear
            ! store them
            rebuf(i,:,:)=se%re_3ph+se%re_4ph
            imbuf(i,:,:)=se%im_3ph+se%im_iso
            ! add minimum smearing
            do j=1,se%nb
            do k=1,opts%nf
                imbuf(i,k,j)=max(imbuf(i,k,j),minsmear)
            enddo
            enddo

            ! Interpolate the shifts
            do j=1,se%nb
                bs%p(i)%shift3(j)=lo_linear_interpolation( se%faxis,se%re_3ph(:,j),bs%p(i)%omega(j) )
                bs%p(i)%shift4(j)=lo_linear_interpolation( se%faxis,se%re_4ph(:,j),bs%p(i)%omega(j) )
            enddo

            if ( mw%talk ) then
                if ( walltime()-t0 .gt. timereport ) then
                    call lo_looptimer('... spectralfunction along path',timer,walltime(),i,bs%nptot)
                    t0=walltime()
                endif
            endif
        enddo
        enddo

        if ( opts%stride .gt. 1 ) then
            t0=walltime()
            bs%selfenergy_real=0.0_flyt
            bs%selfenergy_imag=0.0_flyt
            ! Interpolate the self-energy to all q
            ctr=0
            if ( mw%talk ) call lo_progressbar_init()
            do path=1,bs%npath
                ! fetch the x-values for this path
                l=0
                do q=1,bs%npts,opts%stride
                    i=(path-1)*bs%npts+q
                    l=l+1
                    x(l)=bs%q_axis(i)
                enddo
                
                do band=1,dr%nb
                    do j=1,opts%nf
                        ! fetch self-energies
                        l=0
                        y=0.0_flyt
                        z=0.0_flyt
                        do q=1,bs%npts,opts%stride
                            i=(path-1)*bs%npts+q
                            l=l+1
                            y(l)=rebuf(i,j,band)
                            z(l)=imbuf(i,j,band)
                        enddo
                        ! interpolate self-energies
                        do q=1,bs%npts
                            i=(path-1)*bs%npts+q
                            f0=lo_linear_interpolation(x,y,bs%q_axis(i))
                            f1=lo_linear_interpolation(x,z,bs%q_axis(i))
                            bs%selfenergy_real(i,j,band)=f0
                            bs%selfenergy_imag(i,j,band)=max(f1,0.0_flyt)
                        enddo
                    enddo
                    ctr=ctr+1
                    if ( mw%talk ) call lo_progressbar(' ... interpolating selfenergy',ctr,bs%npath*dr%nb,walltime()-t0)
                enddo
            enddo
            rebuf=bs%selfenergy_real
            imbuf=bs%selfenergy_imag
        endif
        ! and the different intensity axes
        seax=se%faxis
        inax=se%intensityaxis
    end block selfenergy

    ! Figure out some neat interpolation of self-energy for really small q
    smallq: block
        real(flyt), dimension(:), allocatable :: dumre,dumim
        real(flyt), dimension(3) :: qv1,qv2
        real(flyt), dimension(2) :: lsintx,lsinty
        real(flyt) :: f0,t0
        integer :: q,path,ii,jj,i,j,k

        ! Temporarily store self-energies
        lo_allocate(dumre(opts%nf))
        lo_allocate(dumim(opts%nf))
        ! Set smallest imaginary selfenergy
        dumre=0.0_flyt
        dumim=0.0_flyt
        ! Good small number to use
        f0=(seax(2)-seax(1))*0.5_flyt ! smallest selfenergy
        t0=walltime()
        !
        if ( mw%talk ) call lo_progressbar_init()
        do q=1,bs%nptot
            ! what path am I on?
            path=bs%q(q)%path
            ! The start and end-points
            qv1=bs%segment(path)%r1-uc%bz%gshift( bs%segment(path)%r1 + lo_degenvector )
            qv2=bs%segment(path)%r2-uc%bz%gshift( bs%segment(path)%r2 + lo_degenvector )
            ! does it contain gamma?
            if ( norm2(qv1) .gt. lo_tol .and. norm2(qv2) .gt. lo_tol ) cycle
            ! Is it in the beginning or the end?
            if ( norm2(qv1) .lt. lo_tol ) then
                ! Index of gamma
                ii=(path-1)*bs%npts+1
                ! Fix the acoustic branches
                do j=1,dr%nb
                    ! skip if omega too large
                    if ( bs%p(q)%omega(j) .gt. dr%omega_min*0.5_flyt ) cycle
                    ! Find index of point that is ok
                    jj=ii+bs%npts-1
                    do i=ii,ii+bs%npts-1
                        if ( bs%p(i)%omega(j) .gt. dr%omega_min*0.5_flyt ) then
                            jj=i
                            exit
                        endif
                    enddo
                    ! Fetch real and imaginary at this q
                    dumre=bs%selfenergy_real(jj,:,j)
                    dumim=bs%selfenergy_imag(jj,:,j)
                    ! now I know that things are zero at ii, and ok at jj
                    lsintx(1)=bs%q_axis(ii)-lo_sqtol
                    lsintx(2)=bs%q_axis(jj)+lo_sqtol
                    ! Interpolate the missing self-energies at this point
                    do k=1,opts%nf
                        ! y-axis for interpolation, imaginary part
                        lsinty(1)=f0
                        lsinty(2)=dumim(k) 
                        imbuf(q,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(q))
                        ! y-axis for interpolation, real part
                        lsinty(1)=0.0_flyt
                        lsinty(2)=dumre(k) 
                        rebuf(q,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(q))
                    enddo
                enddo
            else
                ! Same thing again, but this time gamma is at the end.
                ii=path*bs%npts
                ! loop over the three lowest branches
                do j=1,dr%nb
                    ! skip if omega too large
                    if ( bs%p(q)%omega(j) .gt. dr%omega_min*0.5_flyt ) cycle
                    jj=(path-1)*bs%npts+1
                    do i=ii,(path-1)*bs%npts+1,-1
                        if ( bs%p(i)%omega(j) .gt. dr%omega_min*0.5_flyt ) then
                            jj=i
                            exit
                        endif
                    enddo
                    ! Fetch real and imaginary at this q
                    dumre=bs%selfenergy_real(jj,:,j)
                    dumim=bs%selfenergy_imag(jj,:,j)
                    ! x-axis for interpolation
                    lsintx(2)=bs%q_axis(ii)-lo_sqtol
                    lsintx(1)=bs%q_axis(jj)+lo_sqtol
                    ! interpolate to missing points
                    do k=1,opts%nf
                        ! y-axis for interpolation
                        lsinty(2)=f0
                        lsinty(1)=dumim(k) 
                        imbuf(q,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(q))
                        ! y-axis for interpolation
                        lsinty(2)=0.0_flyt 
                        lsinty(1)=dumre(k) 
                        rebuf(q,k,j)=lo_linear_interpolation(lsintx,lsinty,bs%q_axis(q))
                    enddo
                enddo
            endif
            if ( mw%talk ) call lo_progressbar(' ... fixing tiny q',q,bs%nptot,walltime()-t0)
        enddo
        ! Store the self-energies
        bs%selfenergy_real=rebuf
        bs%selfenergy_imag=imbuf
    end block smallq

    ! Smear the self-energies in q-and energy direction, just a little 
    smearse: block
        integer :: i,j,band,i1,i2,j1,j2,ii,jj,iii,jjj,ctr
        real(flyt), dimension(5,5) :: kernel
        real(flyt), dimension(2) :: v0
        !
        do i=1,5
        do j=1,5
            v0=[i-3,j-3]*1.0_flyt
            kernel(j,i)=lo_gauss(norm2(v0),0.0_flyt,2.6_flyt)
        enddo
        enddo
        kernel=kernel/sum(kernel)
        !
        rebuf=0.0_flyt
        imbuf=0.0_flyt
        ctr=0
        if ( mw%talk ) call lo_progressbar_init()
        do band=1,dr%nb
            do i=1,bs%nptot
                do j=1,opts%nf
                    i1=max(1,i-2)
                    i2=min(i+2,bs%nptot)
                    j1=max(1,j-2)
                    j2=min(opts%nf,j+2)
                    do ii=i1,i2
                    do jj=j1,j2
                        iii=ii-i+3
                        jjj=jj-j+3
                        rebuf(ii,jj,band)=rebuf(ii,jj,band)+kernel(iii,jjj)*bs%selfenergy_real(i,j,band)
                        imbuf(ii,jj,band)=imbuf(ii,jj,band)+kernel(iii,jjj)*bs%selfenergy_imag(i,j,band)
                    enddo
                    enddo
                enddo
                ctr=ctr+1
            enddo
        enddo
        bs%selfenergy_real=rebuf
        bs%selfenergy_imag=imbuf
    end block smearse

    ! Get the intensities
    intensities: block
        real(flyt), dimension(:), allocatable :: dum
        real(flyt) :: f0,f1
        integer :: i,j,k
        ! put something at gamma to make the intensities not weird.
        f0=(seax(2)-seax(1))*0.25_flyt
        f1=(seax(2)-seax(1))*opts%minsmear

        lo_allocate(dum(opts%nf))
        do i=1,bs%nptot
        do j=1,dr%nb
            if ( bs%p(i)%omega(j) .gt. lo_freqtol ) then
                call getintensity(seax,imbuf(i,:,j),rebuf(i,:,j),bs%p(i)%omega(j),inax,dum)
            else
                do k=1,opts%nf
                    dum(k)=lo_lorentz(inax(k),0.0_flyt,f0)
                enddo
            endif
            bs%intensity(i,:)=bs%intensity(i,:)+dum/lo_trapezoid_integration(inax,dum)
        enddo
        enddo
        bs%faxis=inax

        ! Also store the linewidth at the harmonic frequencies, as well as
        ! the shifts
        do i=1,bs%nptot
        do j=1,dr%nb
            if ( bs%p(i)%omega(j) .gt. lo_freqtol ) then
                f0=lo_linear_interpolation(seax,imbuf(i,:,j),bs%p(i)%omega(j))
                bs%p(i)%linewidth(j)=f0
            else
                bs%p(i)%linewidth(j)=0.0_flyt
            endif 
        enddo
        enddo
    end block intensities

end subroutine

