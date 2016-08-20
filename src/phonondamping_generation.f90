
!> Generate the self-energy at a single q-point
subroutine generate(se,qpoint,ompoint,uc,fc,fct,fcf,qp,dr,loto,opts)
    !> self energy 
    class(lo_phonon_selfenergy), intent(out) :: se
    !> qpoint of interest
    class(lo_qpoint), intent(in) :: qpoint
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: ompoint
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
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
    
    type(lo_listofscatteringrates) :: sr
    type(lo_listofscatteringrates_isotopes) :: sri
    real(flyt), dimension(:), allocatable :: delta
    real(flyt) :: f0 
    integer :: i,j

    ! This routine got a little long, so I chopped it up into blocks. First, make sure we
    ! set all the options, and clean all arrays
    se%nf=opts%nf                               ! Number of points on the energy axis
    se%nb=uc%na*3                               ! Number of bands
    se%sigma=opts%sigma*dr%default_smearing()   ! Baseline fix smearing parameter
    se%adaptiveparameter=opts%sigma             ! Scaling for adaptive gaussian
    se%integrationtype=opts%integrationtype     ! How to integrate
    se%blochlcorrections=.false.                ! dunno
    se%fixedsmearing=.false.                    ! dunno
    se%verbosity=opts%verbosity                 ! How much to talk
    se%loto=opts%loto                           ! Electrostatic corrections
    se%isotope=opts%isotopescattering           ! isotope scattering
    se%thirdorder=opts%thirdorder               ! threephonon term
    se%fourthorder=opts%fourthorder             ! fourphonon term
    se%dos=opts%phonondos                       ! dunno   
    se%slightsmear=opts%slightsmearing          ! smear things a little little bit.
    ! Make space for all the necessary arrays
    if ( allocated(se%faxis) )          deallocate(se%faxis)
    if ( allocated(se%intensityaxis) )  deallocate(se%intensityaxis)
    if ( allocated(se%im_3ph) )         deallocate(se%im_3ph)
    if ( allocated(se%im_iso) )         deallocate(se%im_iso)
    if ( allocated(se%re_3ph) )         deallocate(se%re_3ph)
    if ( allocated(se%re_4ph) )         deallocate(se%re_4ph)
    if ( allocated(se%intensity) )      deallocate(se%intensity)
    if ( allocated(se%twophonondos) )   deallocate(se%twophonondos)
    allocate(se%faxis(se%nf))
    allocate(se%intensityaxis(se%nf))
    allocate(se%im_3ph(se%nf,se%nb))
    allocate(se%im_iso(se%nf,se%nb))
    allocate(se%re_3ph(se%nf,se%nb))
    allocate(se%re_4ph(se%nf,se%nb))
    allocate(se%twophonondos(se%nf,se%nb))
    allocate(se%intensity(se%nf,se%nb))
    ! And set the to appropriate values
    se%im_3ph=0.0_flyt
    se%im_iso=0.0_flyt
    se%re_3ph=0.0_flyt
    se%re_4ph=0.0_flyt
    se%intensity=0.0_flyt
    se%twophonondos=0.0_flyt        
    call lo_linspace(0.0_flyt,2.1_flyt*dr%omega_max,se%faxis)
    call lo_linspace(0.0_flyt,opts%maxf*dr%omega_max,se%intensityaxis)
    ! Store the harmonic properties, and the q-point
    se%q=qpoint
    se%p=ompoint

    ! Now things should be clean and nice. Maybe say what we are about to do?
    if ( se%verbosity .gt. 1 ) then
        write(*,*) '         liso:',se%isotope
        write(*,*) '        lloto:',se%loto
        write(*,*) '  lthirdorder:',se%thirdorder
        write(*,*) ' lfourthorder:',se%fourthorder
        write(*,*) '          dos:',se%dos
        write(*,*) '  slightsmear:',se%slightsmear
        write(*,*) '        qgrid:',qp%griddensity
        write(*,*) '  temperature:',opts%temperature
        write(*,*) '  frequencies:'
        do i=1,dr%nb
        write(*,*) '    mode '//trim(int2char(i))//', omega: ',trim(flyt2char(se%p%omega(i)*opts%unitfactor,8)),' ',opts%unitname
        enddo
    endif

    ! Start to actually calculate stuff, isotope things first
    if ( se%isotope ) then
        sri%verbosity=se%verbosity
        call find_isotope_scatteringrates(se%p,qp,dr,sri,uc)
        select case(se%integrationtype)
            case(1)
                call isotope_imaginary_selfenergy_gaussian(qp,dr,opts%temperature,se,sri)
            case(2)
                call isotope_imaginary_selfenergy_gaussian(qp,dr,opts%temperature,se,sri)
            case(3)
                se%blochlcorrections=.false.
                call isotope_imaginary_selfenergy_tetrahedron(qp,dr,opts%temperature,se,sri)
            case(4)
                se%blochlcorrections=.true.
                call isotope_imaginary_selfenergy_tetrahedron(qp,dr,opts%temperature,se,sri)
        end select
    endif

    ! Then three-phonon things
    if ( se%thirdorder ) then
        sr%verbosity=se%verbosity
        call find_threephonon_scatteringrates(se%q,se%p,qp,dr,sr,uc,fc,fct,loto)
        select case(se%integrationtype)
            case(1)
               call threephonon_imaginary_selfenergy_gaussian(se,sr,qp,dr,opts%temperature)
            case(2)
               call threephonon_imaginary_selfenergy_gaussian(se,sr,qp,dr,opts%temperature)
            case(3)
               se%blochlcorrections=.false.
               call threephonon_imaginary_selfenergy_tetrahedron(qp,sr,dr,opts%temperature,se)
            case(4)
               se%blochlcorrections=.true.
               call threephonon_imaginary_selfenergy_tetrahedron(qp,sr,dr,opts%temperature,se)
        end select
    endif

    ! Maybe fourthorder things    
    if ( se%fourthorder ) then
        allocate(delta(dr%nb))
        call fourphonon_selfenergy(se%q,se%p,qp,opts%temperature,dr,uc,fc,fcf,delta,loto)
        do j=1,se%nb
            se%re_4ph(:,j)=delta(j)
            se%re_4ph(1,j)=0.0_flyt
        enddo
        deallocate(delta)
    endif
    
    ! Make sure the selfenergy is zero where it's supposed to be    
    do i=1,se%nf
        ! no contributions at too high frequency
        if ( se%faxis(i) .gt. dr%omega_max*2 ) se%im_3ph(i,:)=0.0_flyt
        if ( se%faxis(i) .gt. dr%omega_max*2 ) se%im_iso(i,:)=0.0_flyt
        ! no negative values
        do j=1,se%nb
            if ( se%im_3ph(i,j) .lt. 0.0_flyt ) se%im_3ph(i,j)=0.0_flyt
            if ( se%im_iso(i,j) .lt. 0.0_flyt ) se%im_iso(i,j)=0.0_flyt
        enddo
    enddo
    ! zero at zero
    se%im_3ph(1,:)=0.0_flyt
    se%im_iso(1,:)=0.0_flyt
    !
    ! Do a slight smearing of the imaginary selfenergy, it can become really spiky
    ! sometimes. 
    !
    if ( se%slightsmear ) then
        do j=1,se%nb
            if ( se%thirdorder ) then
                f0=lo_trapezoid_integration(se%faxis,se%im_3ph(:,j))
                if( f0*1E-24 .lt. lo_sqtol ) cycle
                call slightsmearing(se%im_3ph(:,j),1)
                se%im_3ph(:,j)=se%im_3ph(:,j)*f0/lo_trapezoid_integration(se%faxis,se%im_3ph(:,j))
            endif
            !
            if ( se%isotope ) then
                f0=lo_trapezoid_integration(se%faxis,se%im_iso(:,j))
                if( f0*1E-24 .lt. lo_sqtol ) cycle
                call slightsmearing(se%im_iso(:,j),1)
                se%im_iso(:,j)=se%im_iso(:,j)*f0/lo_trapezoid_integration(se%faxis,se%im_iso(:,j))
            endif
            se%im_3ph(1,j)=0.0_flyt
            se%im_iso(1,j)=0.0_flyt
        enddo
    endif

    ! Kramer-Kronig transform the imaginary to get the real part 
    call kramer_kronig_transform_to_get_real_part(se)

    ! Some cleanup
    if ( allocated(sr%gi2) )       deallocate(sr%gi2)
    if ( allocated(sr%omega3) )    deallocate(sr%omega3)
    if ( allocated(sr%vel3) )      deallocate(sr%vel3)
    if ( allocated(sr%psisquare) ) deallocate(sr%psisquare)
    if ( allocated(sri%gi2) )      deallocate(sri%gi2)
    if ( allocated(sri%psisquare)) deallocate(sri%psisquare)

end subroutine

!> transform the imaginary part to the real part
subroutine kramer_kronig_transform_to_get_real_part(se)
    !> selfenergy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !
    integer :: i,j
    real(flyt), dimension(:), allocatable :: x,y
    real(flyt) :: sig
    !
    !
    allocate(x(se%nf),y(se%nf))
    sig=se%faxis(2)-se%faxis(1)
    do i=1,se%nf
        do j=1,se%nb
            x=se%faxis
            y=se%im_3ph(:,j) !+se%im_iso(:,j)
            y=y*x
            !y=y/(x(i)**2-x**2)
            y=real(y/(x(i)**2-x**2+lo_imag*sig) )
            y(i)=0.0_flyt
            se%re_3ph(i,j)=lo_trapezoid_integration(x,y)*2.0_flyt/lo_pi
        enddo
    enddo
    !
    deallocate(x,y)
end subroutine
