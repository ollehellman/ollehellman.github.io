
!> Tetrahedron integration weights for isotope scattering, from one q-point
subroutine iso_tetrahedron_fft_oneqp(qp,dr,scq,gi1,scsigma,blochlcorrections)
    !> qpoint mesh
    !type(lo_fft_mesh), intent(in) :: qp
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersoins
    type(lo_phonon_dispersions), intent(in) :: dr
    !> current q-point
    type(lo_iso2), intent(inout) :: scq
    !> current grid-index
    integer, intent(in) :: gi1
    !> baseline smearing parameter
    real(flyt), intent(in) :: scsigma
    !> blochl corrections?
    logical, intent(in) :: blochlcorrections
    !
    real(flyt), parameter :: thres_weight=lo_tiny
    integer, dimension(4) :: tetqpoints
    integer :: gi2
    integer :: i,j,ii,b1,b2,i_gamma
    real(flyt), dimension(4) :: cval,wts1,wts2
    real(flyt) :: om1,omthres,minc,maxc
    real(flyt), dimension(:,:,:), allocatable :: qpw
    real(flyt), dimension(:,:), allocatable :: omr2
    real(flyt), dimension(:), allocatable :: omr1
    
    lo_allocate(omr1(dr%nb))
    lo_allocate(omr2(4,dr%nb))
    lo_allocate(qpw(dr%nb,dr%nb,qp%nq_tot))
    omthres=dr%omega_min*0.2_flyt


    ! Calculate all integration weights for the relevant tetrahedra
    scq%gi1=gi1
    qpw=0.0_flyt
    omr1=dr%aq(gi1)%omega
    tetloop: do i=1,qp%ntet
        ! Grab frequencies
        do j=1,4
            gi2=qp%tet(i)%gridind(j)
            omr2(j,:)=dr%aq(gi2)%omega
            tetqpoints(j)=gi2
        enddo
        
        ! Add integration weights
        do b1=1,dr%nb
        do b2=1,dr%nb
            minc=lo_huge
            maxc=-lo_huge
            do j=1,4
                cval(j)=omr2(j,b2)
                minc=min(minc,cval(j))
                maxc=max(maxc,cval(j))
            enddo
            ! add a small delta for safety
            minc=minc-3*scsigma/100.0_flyt
            maxc=maxc+3*scsigma/100.0_flyt
            ! check if it's relevant
            om1=omr1(b1)
            if ( om1 .gt. minc .and. om1 .lt. maxc .and. om1 .gt. omthres ) then
                ! Take care if it's completely degenerate
                if ( lo_stddev(cval)/om1 .lt. 1E-3_flyt ) then
                    call lo_integration_weights_for_one_tetrahedron(&
                    qp%tet(i),cval,om1-scsigma/200.0_flyt,wts1,scsigma,blochlcorrections)
                    call lo_integration_weights_for_one_tetrahedron(&
                    qp%tet(i),cval,om1+scsigma/200.0_flyt,wts2,scsigma,blochlcorrections)
                else
                    call lo_integration_weights_for_one_tetrahedron(&
                    qp%tet(i),cval,om1-scsigma/200.0_flyt,wts1,scsigma,blochlcorrections)
                    call lo_integration_weights_for_one_tetrahedron(&
                    qp%tet(i),cval,om1+scsigma/200.0_flyt,wts2,scsigma,blochlcorrections)
                endif
                ! sum up integration weights
                qpw(b1,b2,tetqpoints)=qpw(b1,b2,tetqpoints)+&
                (wts1+wts2)*0.5_flyt
            endif
        enddo
        enddo
    enddo tetloop
   
    ! Set the integration weights for the acoustic branches to zero
    i_gamma=-1
    do i=1,qp%nq_tot
        if ( lo_sqnorm(qp%ap(i)%w) .lt. lo_sqtol ) then
            i_gamma=i
            exit
        endif
    enddo
    if ( i_gamma .lt. 0 ) call lo_stop_gracefully(['FFT mesh does not contain gamma'],lo_exitcode_symmetry,__FILE__,__LINE__)

    do b2=1,dr%nb
    do b1=1,dr%nb
        if ( dr%aq(i_gamma)%omega(b1) .lt. lo_freqtol .or. dr%aq(i_gamma)%omega(b2) .lt. lo_freqtol ) then
            qpw(b1,b2,i_gamma)=0.0_flyt
        endif
    enddo
    enddo
 
    ! Store weights per q-point, and only the relevant q-points
    do b1=1,dr%nb
    do b2=1,dr%nb
        ! count q-points
        ii=0
        do i=1,qp%nq_tot
            if ( abs(qpw(b1,b2,i)) .gt. thres_weight ) then
                ii=ii+1
            endif
        enddo
        ! make some space
        scq%band(b1,b2)%n=ii
        if ( scq%band(b1,b2)%n .gt. 0 ) then
           lo_allocate(scq%band(b1,b2)%e( ii ))
        endif
        ! store weights and indices
        ii=0
        do i=1,qp%nq_tot
            if ( abs(qpw(b1,b2,i)) .gt. thres_weight ) then
                ii=ii+1
                scq%band(b1,b2)%e(ii)%gi2=i
                scq%band(b1,b2)%e(ii)%deltafunction=qpw(b1,b2,i)
            endif
        enddo
    enddo
    enddo
    
    lo_deallocate(omr1)
    lo_deallocate(omr2)
    lo_deallocate(qpw)
end subroutine

!> Tetrahedron integration weights for threephonon scattering, from one q-point
subroutine threephonon_tetrahedron_fft_oneqp(qp,dr,scq,gi1,scsigma,blochlcorrections)
    !> qpoint mesh
    !type(lo_fft_mesh), intent(in) :: qp
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> q-point in question
    type(lo_3phqp2), intent(inout) :: scq
    !> grid-index in question
    integer, intent(in) :: gi1
    !> baseline smearing
    real(flyt), intent(in) :: scsigma
    !> blochl corrections?
    logical, intent(in) :: blochlcorrections

    integer :: i,j,ii,jj,b1,b2,b3
    integer :: gi2,gi3,i_gamma
    integer, dimension(3) :: dims
    integer, dimension(4) :: tetqpoints
    real(flyt), dimension(:,:,:,:), allocatable :: qpwp,qpwm
    real(flyt), dimension(:,:), allocatable :: omr2,omr3
    real(flyt), dimension(:), allocatable :: omr1
    real(flyt), dimension(4) :: cval_p,cval_m,wts1,wts2
    real(flyt) :: om1,omthres,minc_m,minc_p,maxc_m,maxc_p,wthres,sigmatol,deltaeps

    ! Some temporary space
    lo_allocate(omr1(dr%nb))
    lo_allocate(omr2(4,dr%nb))
    lo_allocate(omr3(4,dr%nb))
    lo_allocate(qpwp(dr%nb,dr%nb,dr%nb,qp%nq_tot))
    lo_allocate(qpwm(dr%nb,dr%nb,dr%nb,qp%nq_tot))
    ! smallest frequency allowed, should allow everything except acoustic modes at Gamma
    omthres=dr%omega_min*0.5_flyt
    ! smallest weight to care about
    wthres=1E-30_flyt !lo_tiny ! 1E-40_flyt
    ! Dimensions of q-point grid
    select type(qp)
    class is(lo_fft_mesh)
        dims=qp%griddensity
    class default
        write(*,*) 'Really need an fft grid for this'
        stop
    end select
    ! The reference q-point, q1 in q1+q2+q3=G
    scq%gi1=gi1
    ! Some tolerances
    sigmatol=scsigma/30.0_flyt
    deltaeps=scsigma/1E2_flyt
    qpwp=0.0_flyt
    qpwm=0.0_flyt
    omr1=dr%aq(gi1)%omega
    tetloop: do i=1,qp%ntet
        ! Grab frequencies        
        do j=1,4
            gi2=qp%tet(i)%gridind(j)
            gi3=fft_third_grid_index(gi1,gi2,dims)
            omr2(j,:)=dr%aq(gi2)%omega
            omr3(j,:)=dr%aq(gi3)%omega
            tetqpoints(j)=gi2
        enddo
        ! Add integration weights        
        do b1=1,dr%nb
        do b2=1,dr%nb
        do b3=1,dr%nb
            minc_p=lo_huge
            maxc_p=-lo_huge
            minc_m=lo_huge
            maxc_m=-lo_huge
            do j=1,4
                cval_p(j)=omr3(j,b3)-omr2(j,b2)   ! plus events
                cval_m(j)=omr3(j,b3)+omr2(j,b2)   ! minus events
                minc_p=min(minc_p,cval_p(j))
                minc_m=min(minc_m,cval_m(j))
                maxc_p=max(maxc_p,cval_p(j))
                maxc_m=max(maxc_m,cval_m(j))
            enddo
            ! add a small delta for safety, to avid pathological cases
            minc_p=minc_p-sigmatol
            minc_m=minc_m-sigmatol
            maxc_p=maxc_p+sigmatol
            maxc_m=maxc_m+sigmatol
            ! check if it's relevant
            om1=omr1(b1)
            if ( om1 .gt. minc_p .and. om1 .lt. maxc_p ) then !.and. om1 .gt. omthres ) then
                if ( blochlcorrections ) then
                    wts1=lo_LV_tetrahedron_weights(cval_p,om1,lo_freqtol,scsigma)
                    qpwp(b1,b2,b3,tetqpoints)=qpwp(b1,b2,b3,tetqpoints)+wts1*qp%tet(i)%weight
                else
                    wts1=0.0_flyt
                    wts2=0.0_flyt
                    ! this plus event is relevant, store weights
                    call lo_integration_weights_for_one_tetrahedron(qp%tet(i),cval_p,om1-deltaeps,wts1,scsigma,blochlcorrections)
                    call lo_integration_weights_for_one_tetrahedron(qp%tet(i),cval_p,om1+deltaeps,wts2,scsigma,blochlcorrections)
                    qpwp(b1,b2,b3,tetqpoints)=qpwp(b1,b2,b3,tetqpoints)+(wts1+wts2)*0.5_flyt
                endif
            endif
            if ( om1 .gt. minc_m .and. om1 .lt. maxc_m ) then !.and. om1 .gt. omthres ) then
                if ( blochlcorrections ) then
                    wts1=lo_LV_tetrahedron_weights(cval_m,om1,lo_freqtol,scsigma)
                    qpwm(b1,b2,b3,tetqpoints)=qpwm(b1,b2,b3,tetqpoints)+wts1*qp%tet(i)%weight
                else
                ! this minus event is relevant, store weights
                wts1=0.0_flyt
                wts2=0.0_flyt
                call lo_integration_weights_for_one_tetrahedron(qp%tet(i),cval_m,om1-deltaeps,wts1,scsigma,blochlcorrections)
                call lo_integration_weights_for_one_tetrahedron(qp%tet(i),cval_m,om1+deltaeps,wts2,scsigma,blochlcorrections)
                qpwm(b1,b2,b3,tetqpoints)=qpwm(b1,b2,b3,tetqpoints)+(wts1+wts2)*0.5_flyt
                endif
            endif
        enddo
        enddo
        enddo        
    enddo tetloop

    i_gamma=-1
    do i=1,qp%nq_tot
        if ( lo_sqnorm(qp%ap(i)%w) .lt. lo_sqtol ) then
            i_gamma=i
            exit
        endif
    enddo
    if ( i_gamma .lt. 0 ) call lo_stop_gracefully(['FFT mesh does not contain gamma'],lo_exitcode_symmetry,__FILE__,__LINE__)
    do b3=1,dr%nb
    do b2=1,dr%nb
    do b1=1,dr%nb
        if ( dr%aq(i_gamma)%omega(b1) .lt. lo_freqtol .or. &
             dr%aq(i_gamma)%omega(b2) .lt. lo_freqtol .or. &
             dr%aq(i_gamma)%omega(b3) .lt. lo_freqtol ) then
            qpwp(b1,b2,b3,i_gamma)=0.0_flyt
            qpwm(b1,b2,b3,i_gamma)=0.0_flyt
        endif
    enddo
    enddo
    enddo
    qpwp(1:3,1:3,1:3,i_gamma)=0.0_flyt
    qpwm(1:3,1:3,1:3,i_gamma)=0.0_flyt
    
    ! Store weights per q-point    
    do b1=1,dr%nb
    do b2=1,dr%nb
    do b3=1,dr%nb
        ! count q-points
        ii=0
        jj=0
        do i=1,qp%nq_tot
            if ( abs(qpwp(b1,b2,b3,i)) .gt. wthres ) then
                ii=ii+1
            endif
            if ( abs(qpwm(b1,b2,b3,i)) .gt. wthres ) then
                jj=jj+1
            endif
        enddo
        ! make some space
        scq%plus(b1,b2,b3)%n=ii
        scq%minus(b1,b2,b3)%n=jj
        if ( scq%plus(b1,b2,b3)%n .gt. 0 ) then
           lo_allocate(scq%plus(b1,b2,b3)%e( ii ))
        endif
        if ( scq%minus(b1,b2,b3)%n .gt. 0 ) then
           lo_allocate(scq%minus(b1,b2,b3)%e( jj ))
        endif
        ! store weights and indices
        ii=0
        jj=0
        do i=1,qp%nq_tot
            if ( abs(qpwp(b1,b2,b3,i)) .gt. wthres ) then
                ii=ii+1
                gi2=i
                scq%plus(b1,b2,b3)%e(ii)%gi2=gi2
                scq%plus(b1,b2,b3)%e(ii)%gi3=fft_third_grid_index(gi1,gi2,dims)
                scq%plus(b1,b2,b3)%e(ii)%deltafunction=qpwp(b1,b2,b3,i)
            endif
            if ( abs(qpwm(b1,b2,b3,i)) .gt. wthres ) then
                jj=jj+1
                gi2=i
                scq%minus(b1,b2,b3)%e(jj)%gi2=gi2
                scq%minus(b1,b2,b3)%e(jj)%gi3=fft_third_grid_index(gi1,gi2,dims)
                scq%minus(b1,b2,b3)%e(jj)%deltafunction=qpwm(b1,b2,b3,i)
            endif
        enddo
    enddo
    enddo
    enddo

    ! Cleanup    
    lo_deallocate(qpwp)
    lo_deallocate(qpwm)
    lo_deallocate(omr1)
    lo_deallocate(omr2)
    lo_deallocate(omr3)    
end subroutine

