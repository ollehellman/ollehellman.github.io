
!> count three-phonons scattering events with gaussian smearing
subroutine threephonon_gaussian_fft_oneqp(qp,dr,scq,gi1,scsigma,thres,adaptiveparameter,progressbar)
    !type(lo_fft_mesh), intent(in) :: qp
    class(lo_qpoint_mesh), intent(in) :: qp
    type(lo_phonon_dispersions), intent(in) :: dr
    type(lo_3phqp2), intent(inout) :: scq
    integer, intent(in) :: gi1
    real(flyt), intent(in) :: scsigma,thres
    real(flyt), intent(in), optional :: adaptiveparameter
    logical, intent(in), optional :: progressbar
    !
    integer :: i
    integer :: b1,b2,b3
    integer :: gi2,gi3
    integer, dimension(3) :: dims
    integer, dimension(:,:,:), allocatable :: bc1,bc2
    real(flyt), dimension(:,:), allocatable :: vel2,vel3
    real(flyt), dimension(:), allocatable :: omr1,omr2,omr3
    real(flyt) :: deltafunction,om1,om2,om3,sigma,omthres
    logical :: adaptive,progress
    
    ! Adaptive or fixed gaussian    
    if ( present(adaptiveparameter) ) then
        adaptive=.true.
    else
        adaptive=.false.
    endif
    ! Show a progress bar?
    if ( present(progressbar) ) then
        progress=progressbar
    else
        progress=.false.
    endif
    
    ! Threshold for small frequencies
    omthres=dr%omega_min*0.2_flyt    

    ! Some temporary space    
    lo_allocate(bc1(dr%nb,dr%nb,dr%nb))
    lo_allocate(bc2(dr%nb,dr%nb,dr%nb))
    lo_allocate(omr1(dr%nb))
    lo_allocate(omr2(dr%nb))
    lo_allocate(omr3(dr%nb))
    lo_allocate(vel2(3,dr%nb))
    lo_allocate(vel3(3,dr%nb))
    ! grid dimensions
    select type(qp)
    class is(lo_fft_mesh)
        dims=qp%griddensity
    class default
        write(*,*) 'Really need an fft grid for this'
        stop
    end select

    ! Set the q-index and omega for q, and reset counters
    scq%gi1=gi1
    bc1=0
    bc2=0    
    omr1=dr%aq(gi1)%omega

    ! Do the actual counting
    if ( progress ) call lo_progressbar_init() 
    do i=1,qp%nq_tot        
        ! Get q'', This is q1+q2+q3=G        
        gi2=i
        gi3=fft_third_grid_index(gi1,gi2,dims)
        !
        omr2=dr%aq(gi2)%omega
        omr3=dr%aq(gi3)%omega
        vel2=dr%aq(gi2)%vel
        vel3=dr%aq(gi3)%vel
        
        ! Count events        
        do b1=1,dr%nb
        do b2=1,dr%nb
        do b3=1,dr%nb            
            om1=omr1(b1)
            om2=omr2(b2)
            om3=omr3(b3)
            ! plus-events, first get sigma
            if ( adaptive ) then
                sigma=qp%smearingparameter(vel2(:,b2)-vel3(:,b3),scsigma,adaptiveparameter)
            else
                sigma=scsigma
            endif            
            if ( abs(om1-om3+om2) .lt. thres*sigma ) then
            if ( om1 .gt. omthres .and. om2 .gt. omthres .and. om3 .gt. omthres ) then
                bc1(b1,b2,b3)=bc1(b1,b2,b3)+1
            endif
            endif
            ! minus-events
            if ( adaptive ) then
                sigma=qp%smearingparameter(vel2(:,b2)+vel3(:,b3),scsigma,adaptiveparameter)
            else
                sigma=scsigma
            endif            
            if ( abs(om1-om2-om3) .lt. thres*sigma ) then
            if ( om1 .gt. omthres .and. om2 .gt. omthres .and. om3 .gt. omthres ) then
                bc2(b1,b2,b3)=bc2(b1,b2,b3)+1
            endif
            endif
        enddo
        enddo
        enddo
        if ( progress ) call lo_progressbar(' ... counting scattering events',i,qp%nq_tot*2)
    enddo
    
    ! Allocate the storage    
    do b1=1,dr%nb
    do b2=1,dr%nb
    do b3=1,dr%nb
        scq%plus(b1,b2,b3)%n=bc1(b1,b2,b3)
        scq%minus(b1,b2,b3)%n=bc2(b1,b2,b3)
        if ( scq%plus(b1,b2,b3)%n .gt. 0 ) then
            lo_allocate(scq%plus(b1,b2,b3)%e( bc1(b1,b2,b3) ))
        endif
        if ( scq%minus(b1,b2,b3)%n .gt. 0 ) then
            lo_allocate(scq%minus(b1,b2,b3)%e( bc2(b1,b2,b3) ))
        endif
    enddo
    enddo
    enddo
    
    ! Count again and store things
    bc1=0
    bc2=0
    do i=1,qp%nq_tot        
        ! This is q1+q2+q3=G        
        gi2=i
        gi3=fft_third_grid_index(gi1,gi2,dims)        
        omr2=dr%aq(gi2)%omega
        omr3=dr%aq(gi3)%omega
        vel2=dr%aq(gi2)%vel
        vel3=dr%aq(gi3)%vel
        do b1=1,dr%nb
        do b2=1,dr%nb
        do b3=1,dr%nb
            om1=omr1(b1)
            om2=omr2(b2)
            om3=omr3(b3)
            ! plus-events, first get sigma
            if ( adaptive ) then
                sigma=qp%smearingparameter(vel2(:,b2)-vel3(:,b3),scsigma,adaptiveparameter)
            else
                sigma=scsigma
            endif            
            if ( abs(om1-om3+om2) .lt. thres*sigma ) then
            if ( om1 .gt. omthres .and. om2 .gt. omthres .and. om3 .gt. omthres ) then
                deltafunction=lo_gauss(om1,om3-om2,sigma)
                bc1(b1,b2,b3)=bc1(b1,b2,b3)+1
                scq%plus(b1,b2,b3)%e( bc1(b1,b2,b3) )%gi2=gi2
                scq%plus(b1,b2,b3)%e( bc1(b1,b2,b3) )%gi3=gi3
                scq%plus(b1,b2,b3)%e( bc1(b1,b2,b3) )%deltafunction=deltafunction/qp%nq_tot
            endif
            endif
            ! minus-events
            if ( adaptive ) then
                sigma=qp%smearingparameter(vel2(:,b2)+vel3(:,b3),scsigma,adaptiveparameter)
            else
                sigma=scsigma
            endif
            if ( abs(om1-om2-om3) .lt. thres*sigma ) then
            if ( om1 .gt. omthres .and. om2 .gt. omthres .and. om3 .gt. omthres ) then
                bc2(b1,b2,b3)=bc2(b1,b2,b3)+1
                deltafunction=lo_gauss(om1,om2+om3,sigma)
                scq%minus(b1,b2,b3)%e( bc2(b1,b2,b3) )%gi2=gi2
                scq%minus(b1,b2,b3)%e( bc2(b1,b2,b3) )%gi3=gi3
                scq%minus(b1,b2,b3)%e( bc2(b1,b2,b3) )%deltafunction=deltafunction/qp%nq_tot
            endif
            endif
        enddo
        enddo
        enddo
        if ( progress ) call lo_progressbar(' ... counting scattering events',qp%nq_tot+i,qp%nq_tot*2)
    enddo
    ! Cleanup
    lo_deallocate(bc1)
    lo_deallocate(bc2)
    lo_deallocate(omr1)
    lo_deallocate(omr2)
    lo_deallocate(omr3)
    lo_deallocate(vel2)
    lo_deallocate(vel3)
end subroutine

!> Gaussian integration weights for isotope scattering, from one q-point
subroutine iso_gaussian_fft_oneqp(qp,dr,scq,gi1,scsigma,thres,adaptiveparameter)
    !> q-point mesh
    !type(lo_fft_mesh), intent(in) :: qp
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> qpoint in question
    type(lo_iso2), intent(inout) :: scq
    !> gridindex in question
    integer, intent(in) :: gi1
    !> baseline smearing
    real(flyt), intent(in) :: scsigma
    !> threshold to cut of gaussian
    real(flyt), intent(in) :: thres
    !> scaling factor for adaptive gaussian
    real(flyt), intent(in), optional :: adaptiveparameter
    !
    integer :: i,ii,b1,b2
    integer :: gi2
    integer, dimension(:,:), allocatable :: bandcounter
    real(flyt) :: deltafunction,om1,om2,omthres,sigma
    logical :: adaptive
    
    ! Fix or adaptive gaussian?    
    if ( present(adaptiveparameter) ) then
        adaptive=.true.
    else
        adaptive=.false.
    endif
    ! Values for q
    omthres=dr%omega_min*0.2_flyt
    scq%gi1=gi1

    ! count first
    lo_allocate(bandcounter(dr%nb,dr%nb))
    bandcounter=0
    do i=1,qp%nq_tot
        gi2=i
        do b1=1,dr%nb
        do b2=1,dr%nb
            om1=dr%aq(gi1)%omega(b1)
            om2=dr%aq(gi2)%omega(b2)
            if ( om1 .gt. omthres .and. om2 .gt. omthres ) then
                if ( adaptive ) then
                    sigma=qp%smearingparameter(dr%aq(gi2)%vel(:,b2),scsigma,adaptiveparameter)
                else
                    sigma=scsigma
                endif
                if ( abs(om1-om2) .lt. thres*sigma ) then
                    bandcounter(b1,b2)=bandcounter(b1,b2)+1
                endif
            endif
        enddo
        enddo
    enddo
    
    ! Allocate storage
    do b1=1,dr%nb
    do b2=1,dr%nb
        scq%band(b1,b2)%n=bandcounter(b1,b2)
        if ( scq%band(b1,b2)%n .gt. 0 ) then
           lo_allocate(scq%band(b1,b2)%e( scq%band(b1,b2)%n ))
        endif
    enddo
    enddo
    
    ! Count again and store
    bandcounter=0
    do i=1,qp%nq_tot
        gi2=i
        ! It should not bounce to itself, perhaps
        do b1=1,dr%nb
        do b2=1,dr%nb
            om1=dr%aq(gi1)%omega(b1)
            om2=dr%aq(gi2)%omega(b2)
            if ( om1 .gt. omthres .and. om2 .gt. omthres ) then
                if ( adaptive ) then
                    sigma=qp%smearingparameter(dr%aq(gi2)%vel(:,b2),scsigma,adaptiveparameter)
                else
                    sigma=scsigma
                endif
                if ( abs(om1-om2) .lt. thres*sigma ) then
                    bandcounter(b1,b2)=bandcounter(b1,b2)+1
                    ii=bandcounter(b1,b2)
                    deltafunction=lo_gauss(om1,om2,sigma)
                    scq%band(b1,b2)%e(ii)%deltafunction=deltafunction/qp%nq_tot
                    scq%band(b1,b2)%e(ii)%gi2=gi2
                endif
            endif
        enddo
        enddo
    enddo
    lo_deallocate(bandcounter)
end subroutine

