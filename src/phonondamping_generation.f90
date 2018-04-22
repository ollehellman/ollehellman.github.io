
!> Generate the self-energy at a single q-point
subroutine generate(se,qpoint,ompoint,uc,fc,fct,fcf,qp,dr,opts,mw)
    !> self energy 
    class(lo_phonon_selfenergy), intent(out) :: se
    !> qpoint of interest
    class(lo_qpoint), intent(in) :: qpoint
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: ompoint
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
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
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw

    type(lo_listofscatteringrates) :: sr

    ! First thing to do is to initialize all arrays, and set all options
    setopts: block
        integer :: i
        se%nf=opts%nf                               ! Number of points on the energy axis
        se%nb=uc%na*3                               ! Number of bands
        se%sigma=opts%sigma*dr%default_smearing()   ! Baseline fix smearing parameter
        se%adaptiveparameter=opts%sigma             ! Scaling for adaptive gaussian
        se%integrationtype=opts%integrationtype     ! How to integrate
        se%blochlcorrections=.false.                ! dunno
        se%verbosity=opts%verbosity                 ! How much to talk
        se%isotope=opts%isotopescattering           ! isotope scattering
        se%thirdorder=opts%thirdorder               ! threephonon term
        se%fourthorder=opts%fourthorder             ! fourphonon term
        se%dos=opts%phonondos                       ! dunno   
        se%slightsmear=opts%slightsmearing          ! smear things a little little bit.
        se%diagonal=opts%diagonal                   ! calculate only the diagonal part of the self-energy
        ! Make space for all the necessary arrays
        if ( allocated(se%faxis) )          deallocate(se%faxis)
        if ( allocated(se%intensityaxis) )  deallocate(se%intensityaxis)
        if ( allocated(se%im_3ph) )         deallocate(se%im_3ph)
        if ( allocated(se%im_iso) )         deallocate(se%im_iso)
        if ( allocated(se%re_3ph) )         deallocate(se%re_3ph)
        if ( allocated(se%re_4ph) )         deallocate(se%re_4ph)
        if ( allocated(se%intensity) )      deallocate(se%intensity)
!        if ( allocated(se%im_nd) )          deallocate(se%im_nd)
!        if ( allocated(se%re_nd) )          deallocate(se%re_nd)
        allocate(se%faxis(se%nf))
        allocate(se%intensityaxis(se%nf))
        allocate(se%im_3ph(se%nf,se%nb))
        allocate(se%im_iso(se%nf,se%nb))
        allocate(se%re_3ph(se%nf,se%nb))
        allocate(se%re_4ph(se%nf,se%nb))
        allocate(se%intensity(se%nf,se%nb))
!        allocate(se%im_nd(se%nf,se%nb,se%nb))
!        allocate(se%re_nd(se%nf,se%nb,se%nb))
        ! And set the to appropriate values
        se%im_3ph=0.0_flyt
        se%im_iso=0.0_flyt
        se%re_3ph=0.0_flyt
        se%re_4ph=0.0_flyt
!        se%re_nd=0.0_flyt
!        se%im_nd=0.0_flyt
        se%intensity=0.0_flyt
        call lo_linspace(0.0_flyt,2.1_flyt*dr%omega_max,se%faxis)
        call lo_linspace(0.0_flyt,opts%maxf*dr%omega_max,se%intensityaxis)
        ! Store the harmonic properties, and the q-point
        se%q%w=qpoint%w
        se%q%v=qpoint%v
        call lo_get_small_group_of_qpoint(se%q,uc)
        se%p=ompoint
        ! Now things should be clean and nice. Maybe say what we are about to do?
        if ( se%verbosity .gt. 1 ) then
            write(*,*) '         liso:',se%isotope
            write(*,*) '  lthirdorder:',se%thirdorder
            write(*,*) ' lfourthorder:',se%fourthorder
            write(*,*) '          dos:',se%dos
            write(*,*) '  slightsmear:',se%slightsmear
            write(*,*) '        sigma:',se%sigma*lo_frequency_Hartree_to_THz
            write(*,*) '  temperature: ',tochar(opts%temperature),' K'
            write(*,*) '  frequencies:'
            do i=1,dr%nb
                write(*,*) '    mode ',tochar(i,-3),', omega: ',tochar(se%p%omega(i)*lo_frequency_Hartree_to_THz,6),' THz'
            enddo
        endif
        ! Divide everything over MPI. Right now, I only divide over qpoint I think.
        if ( qp%mpi%initialized .eqv. .false. ) then
            write(*,*) 'q-grid not distributed over MPI'
            stop
        endif
    end block setopts

    ! Now calculate all the matrix elements
    call sr%generate(se%q,se%p,qp,dr,uc,fc,fct,se%verbosity,se%isotope,se%thirdorder,mw)

    ! Start to actually calculate stuff, isotope things first
    if ( se%isotope ) then
    select case(se%integrationtype)
        case(1)
            call isotope_imaginary_selfenergy_gaussian(qp,dr,opts%temperature,se,sr,mw)
        case(2)
            call isotope_imaginary_selfenergy_gaussian(qp,dr,opts%temperature,se,sr,mw)
        case(3)
            se%blochlcorrections=.false.
            call isotope_imaginary_selfenergy_tetrahedron(qp,dr,opts%temperature,se,sr,mw)
        case(4)
            se%blochlcorrections=.true.
            call isotope_imaginary_selfenergy_tetrahedron(qp,dr,opts%temperature,se,sr,mw)
    end select
    endif

    ! Then three-phonon things
    if ( se%thirdorder ) then
    select case(se%integrationtype)
       case(1)
          call threephonon_imaginary_selfenergy_gaussian(se,sr,qp,dr,opts%temperature,mw)
       case(2)
          call threephonon_imaginary_selfenergy_gaussian(se,sr,qp,dr,opts%temperature,mw)
       case(3)
          se%blochlcorrections=.false.
          call threephonon_imaginary_selfenergy_tetrahedron(qp,sr,dr,opts%temperature,se,mw)
       case(4)
          se%blochlcorrections=.true.
          call threephonon_imaginary_selfenergy_tetrahedron(qp,sr,dr,opts%temperature,se,mw)
    end select        
    endif

    ! Maybe fourthorder things
    if ( se%fourthorder ) then
    fourthorder: block
        real(flyt), dimension(:), allocatable :: delta
        integer :: j
        allocate(delta(dr%nb))
        call fourphonon_selfenergy(se%q,se%p,qp,opts%temperature,dr,uc,fc,fcf,delta,mw,se%verbosity)
        do j=1,se%nb
            se%re_4ph(:,j)=delta(j)
            se%re_4ph(1,j)=0.0_flyt
        enddo
        deallocate(delta)
    end block fourthorder
    endif

    !> finalize to ensure that it's reasonable.
    sanity: block
        real(flyt) :: f0,f1
        integer :: i,j
        ! Make sure the selfenergy is zero where it's supposed to be
        do j=1,se%nb
        do i=1,se%nf
            if ( se%im_3ph(i,j) .lt. 0.0_flyt ) se%im_3ph(i,j)=0.0_flyt 
            if ( se%im_iso(i,j) .lt. 0.0_flyt ) se%im_iso(i,j)=0.0_flyt 
        enddo
        enddo
        ! zero at zero
        se%im_3ph(1,:)=0.0_flyt
        se%im_iso(1,:)=0.0_flyt

        ! Do a slight smearing of the imaginary selfenergy, it can become 
        ! really spiky sometimes. 
        if ( se%slightsmear ) then
            do j=1,se%nb
                if ( se%thirdorder ) then
                    f0=lo_trapezoid_integration(se%faxis,se%im_3ph(:,j))
                    if( f0*se%nf .lt. lo_freqtol ) cycle
                    call slightsmearing(se%im_3ph(:,j),1)
                    se%im_3ph(:,j)=se%im_3ph(:,j)*f0/lo_trapezoid_integration(se%faxis,se%im_3ph(:,j))
                endif
                !
                if ( se%isotope ) then
                    f0=lo_trapezoid_integration(se%faxis,se%im_iso(:,j))
                    if( f0*se%nf .lt. lo_freqtol ) cycle
                    call slightsmearing(se%im_iso(:,j),1)
                    se%im_iso(:,j)=se%im_iso(:,j)*f0/lo_trapezoid_integration(se%faxis,se%im_iso(:,j))
                endif
                se%im_3ph(1,j)=0.0_flyt
                se%im_iso(1,j)=0.0_flyt
            enddo
        endif
    end block sanity

    ! Kramer-Kronig transform the imaginary to get the real part 
    call kramer_kronig_transform_to_get_real_part(se)

    ! Add the isotope part do the non-diagonal self-energy
!    do i=1,se%nf
!        do j=1,se%nb
!            se%im_nd(i,j,j)=se%im_nd(i,j,j)+se%im_iso(i,j)
!        enddo
!    enddo
!    ! and maybe the fourthorder
!    if ( se%fourthorder ) then
!        if ( se%diagonal ) then
!            do i=1,se%nf
!                do j=1,se%nb
!                    se%re_nd(i,j,j)= se%re_nd(i,j,j)+se%re_4ph(i,j)
!                enddo
!            enddo
!        else
!            do i=1,se%nf
!                se%re_nd(i,:,:)=se%re_nd(i,:,:)+nddelta
!            enddo
!        endif
!    endif

end subroutine

!> transform the imaginary part to the real part
subroutine kramer_kronig_transform_to_get_real_part(se)
    !> selfenergy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !
    integer :: i,j,j1,j2
    real(flyt), dimension(:), allocatable :: x,y
    real(flyt) :: sig
    
    allocate(x(se%nf),y(se%nf))
    sig=(se%faxis(2)-se%faxis(1))*lo_sqtol 
    x=se%faxis

    do i=1,se%nf
        do j=1,se%nb
            y=se%im_3ph(:,j) 
            y=y*x
            !y=real(y/(x(i)**2-x**2+lo_imag*sig) )
            y=real(y/(x(i)**2-x**2+lo_imag*sig) )
            !y(i)=0.0_flyt
            se%re_3ph(i,j)=lo_trapezoid_integration(x,y)*2.0_flyt/lo_pi
        enddo
    enddo

!    se%re_nd=0.0_flyt
!    if ( se%diagonal ) then
!        do j=1,se%nb
!            se%re_nd(:,j,j)=se%re_3ph(:,j)
!        enddo
!    else
!        do i=1,se%nf
!            do j1=1,se%nb
!            do j2=1,se%nb
!                y=se%im_nd(:,j1,j2) 
!                y=y*x
!                y=real(y/(x(i)**2-x**2+lo_imag*sig) )
!                !y(i)=0.0_flyt
!                se%re_nd(i,j1,j2)=lo_trapezoid_integration(x,y)*2.0_flyt/lo_pi
!            enddo
!            enddo
!        enddo
!    endif
        
    deallocate(x,y)
end subroutine
