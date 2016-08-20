
!> The fourth order self-energy
subroutine fourphonon_selfenergy(qpoint,ompoint,qp,temperature,dr,uc,fc,fcf,delta,loto)
    !> qpoint for q
    type(lo_qpoint), intent(in) :: qpoint
    !> harmonic properties for q
    type(lo_phonon_dispersions_qpoint), intent(in) :: ompoint
    !> grid for q',q'',q'''
    type(lo_monkhorst_pack_mesh), intent(in) :: qp
    !> harmonic properties for q',q'',q''
    type(lo_phonon_dispersions), intent(in) :: dr
    !> cyrstal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(in) :: fc
    !> third order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> electrostatic corrections
    type(lo_loto), intent(in) :: loto
    !> temperature
    real(flyt), intent(in) :: temperature
    !> real four-phonon self-energy
    real(flyt), dimension(:), intent(out) :: delta
    !
    integer :: q1,b1,b2,b3,b4,gi3
    complex(flyt) :: c0
    real(flyt) :: t0,omegathres,enhet
    real(flyt), dimension(3) :: qv1,qv2,qv3,qv4
    real(flyt), dimension(4) :: omega
    complex(flyt), dimension(:,:), allocatable :: egv,dum
    !
    ! For the extra q-points
    !
    type(lo_qpoint) :: dumq2,dumq3,dumq4
    complex(flyt), dimension(:,:,:), allocatable :: Dq
    complex(flyt), dimension(:,:), allocatable :: D
    complex(flyt), dimension(:,:), allocatable :: dumegv2,dumegv3,dumegv4
    real(flyt), dimension(:,:), allocatable :: dumvel2,dumvel3,dumvel4
    real(flyt), dimension(:), allocatable :: dumomega2,dumomega3,dumomega4
    !
    ! Dummy array to hold things 
    !
    t0=walltime()
    allocate(dum(dr%nb,qp%nq_tot))
    dum=0.0_flyt
    !$OMP PARALLEL DEFAULT(private) SHARED(dr,fc,fcf,qp,dum,uc,loto,temperature,qpoint,ompoint)
    omegathres=dr%omega_min*0.2_flyt
    !$OMP CRITICAL
    lo_allocate(dumegv2(dr%nb,dr%nb))
    lo_allocate(dumegv3(dr%nb,dr%nb))
    lo_allocate(dumegv4(dr%nb,dr%nb))
    lo_allocate(dumvel2(3,dr%nb))
    lo_allocate(dumvel3(3,dr%nb))
    lo_allocate(dumvel4(3,dr%nb))
    lo_allocate(dumomega2(dr%nb))    
    lo_allocate(dumomega3(dr%nb))    
    lo_allocate(dumomega4(dr%nb))
    lo_allocate(D(dr%nb,dr%nb))
    lo_allocate(Dq(3,dr%nb,dr%nb))
    ! Get q2
    dumq2%v=-qpoint%w
    dumq2%w=dumq2%v-uc%bz%gshift(dumq2%v)
    call lo_get_small_group_of_qpoint(dumq2,uc)
    call lo_get_dynamical_matrix(fc,uc,dumq2,loto,D,Dq)
    call lo_get_omega_and_velocities(D,Dq,uc,dumomega2,dumegv2,dumvel2,qpoint=dumq2)
    lo_deallocate(dumq2%invariant_operations)
    ! Some other stuff
    lo_allocate(egv(dr%nb,4))
    !$OMP END CRITICAL
    !
    !$OMP DO
    do q1=1,qp%nq_tot        
        !
        ! Get the other q-points
        !
        gi3=q1
        dumq3%v=qp%ap(gi3)%w
        dumq3%w=dumq3%v-uc%bz%gshift(dumq3%v)
        call lo_get_small_group_of_qpoint(dumq3,uc)
        call lo_get_dynamical_matrix(fc,uc,dumq3,loto,D,Dq)
        call lo_get_omega_and_velocities(D,Dq,uc,dumomega3,dumegv3,dumvel3,qpoint=dumq3)
        lo_deallocate(dumq3%invariant_operations)
        dumq4%v=-dumq3%w
        dumq4%w=dumq4%v-uc%bz%gshift(dumq4%v)
        call lo_get_small_group_of_qpoint(dumq4,uc)
        call lo_get_dynamical_matrix(fc,uc,dumq4,loto,D,Dq)
        call lo_get_omega_and_velocities(D,Dq,uc,dumomega4,dumegv4,dumvel4,qpoint=dumq4)
        lo_deallocate(dumq4%invariant_operations)
        !
        qv1=qpoint%w*2*lo_pi
        qv2=-qv1
        qv3=dumq3%w*2*lo_pi 
        qv4=-qv3
        !
        do b1=1,dr%nb
            b2=b1
            do b3=1,dr%nb
                b4=b3
                !
                ! And freqencies
                !
                omega(1)=ompoint%omega(b1)
                omega(2)=dumomega2(b2)
                omega(3)=dumomega3(b3)
                omega(4)=dumomega4(b4)
                !
                ! As well as eigenvectors
                !
                egv(:,1)=ompoint%egv(:,b1)
                egv(:,2)=dumegv2(:,b2)
                egv(:,3)=dumegv3(:,b3)
                egv(:,4)=dumegv4(:,b4)
                !
                ! To finally get the scattering rates
                !
                if ( minval(omega) .gt. omegathres ) then
                    !
                    c0=fcf%scatteringamplitude(uc,omega,egv,qv2,qv3,qv4)
                    c0=c0*(1.0_flyt+lo_planck(temperature,omega(3))+lo_planck(temperature,omega(4)) )
                    !
                else
                    c0=0.0_flyt
                endif
                !
                dum(b1,q1)=dum(b1,q1)+real(c0)
                !
            enddo
        enddo
        !
    enddo
    !$OMP END DO    
    lo_deallocate(dumegv2)
    lo_deallocate(dumegv3)
    lo_deallocate(dumegv4)
    lo_deallocate(dumvel2)
    lo_deallocate(dumvel3)
    lo_deallocate(dumvel4)
    lo_deallocate(dumomega2)    
    lo_deallocate(dumomega3)    
    lo_deallocate(dumomega4)
    lo_deallocate(D)
    lo_deallocate(Dq)
    lo_deallocate(egv)
    !$OMP END PARALLEL
    !
    do b1=1,dr%nb
        delta(b1)=sum(real(dum(b1,:)))
    enddo
    !
    enhet=lo_hbar_J/16.0_flyt/qp%nq_tot
    delta=delta*enhet
    !
    lo_deallocate(dum)
end subroutine

