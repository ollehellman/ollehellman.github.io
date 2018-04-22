!> calculate things as a density of states
subroutine get_intensity_as_dos(pd,qpd,drd,uc,fc,fct,fcf,qp,dr,opts,mw)
    !> phonon density of states
    type(lo_phonon_dos), intent(out) :: pd
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qpd
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(inout) :: drd
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
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw

    ! stuff for the DOS mesh
    real(flyt) :: tt0

    ! Start the timer
    init: block
        tt0=walltime()
        ! make space in the dos
        pd%na=uc%na
        pd%nb=uc%na*3
        pd%dosmin=0.0_flyt
        pd%dosmax=drd%omega_max*opts%maxf
        pd%ndos=opts%nf
        pd%enhet=opts%enhet
        pd%verbosity=opts%verbosity
        pd%integrationtype=-1   ! no choice here
        pd%a=opts%sigma
        pd%dossmear=drd%default_smearing()
        lo_allocate(pd%omega(pd%ndos))
        lo_allocate(pd%dos(pd%ndos))
        lo_allocate(pd%pdos_site(pd%ndos,pd%na))
        lo_allocate(pd%pdos_mode(pd%ndos,pd%nb))
        call lo_linspace(0.0_flyt,pd%dosmax,pd%omega)
        pd%dos=0.0_flyt
        pd%pdos_site=0.0_flyt
        pd%pdos_mode=0.0_flyt
    end block init

    ! We can start by calculating the lineshape at all the points we need.
    lshp: block
        type(lo_phonon_selfenergy) :: se
        complex(flyt), dimension(3) :: cv0
        real(flyt), dimension(:,:), allocatable :: lsbuf
        real(flyt), dimension(uc%na) :: siteproj
        real(flyt) :: t0,sigma
        integer :: i,j,k

        t0=walltime()
        ! Make some space
        lo_allocate(lsbuf(opts%nf,dr%nb))
        lsbuf=0.0_flyt

        if ( mw%talk ) call lo_progressbar_init()
        do i=1,qpd%nq_irr
            ! Get the self-energy
            call se%generate(qpd%ip(i),drd%iq(i),uc,fc,fct,fcf,qp,dr,opts,mw)
            ! Get the intensity
            lsbuf=0.0_flyt
            do j=1,dr%nb
                ! get the smearing parameter
                sigma=qpd%smearingparameter(drd%iq(i)%vel(:,j),pd%dossmear,pd%a)
                ! the site-projections
                do k=1,uc%na
                    cv0=drd%iq(i)%egv( (k-1)*3+1:k*3, j )
                    siteproj(k)=abs(dot_product(cv0,conjg(cv0)))
                enddo
                ! intensity, smeared by this
                call getintensity(pd%omega,se%im_3ph(:,j)+se%im_iso(:,j),se%re_3ph(:,j)+se%re_4ph(:,j),&
                    drd%iq(i)%omega(j),se%faxis,lsbuf(:,j),sigma)
                ! Add this in the right place
                pd%pdos_mode(:,j)=pd%pdos_mode(:,j)+lsbuf(:,j)*qpd%ip(i)%weight
                do k=1,uc%na
                    pd%pdos_site(:,k)=pd%pdos_site(:,k)+lsbuf(:,j)*siteproj(k)*qpd%ip(i)%weight
                enddo
            enddo
            if ( mw%talk ) call lo_progressbar(' ... lineshapes across q-mesh',i,qpd%nq_irr,walltime()-t0)
        enddo
    end block lshp

    ! Clean up the dos so that it makes sense
    cleandos: block
        real(flyt), dimension(:,:), allocatable :: sitebuf
        real(flyt) :: f0,f1
        integer :: i,j,k

        ! I might have gotten a contribution at 0 frequency due to the smearing. That has to
        ! be removed. I subtract the line that goes from the dos at 0 to the max frequency.
        do j=1,pd%nb
            f1=pd%pdos_mode(1,j)
            do i=1,pd%ndos
                f0=(pd%ndos-i)*1.0_flyt/( (pd%ndos-1)*1.0_flyt )
                f0=f0**2
                pd%pdos_mode(i,j)=pd%pdos_mode(i,j)-f0*f1
                if ( pd%pdos_mode(i,j) .lt. 0.0_flyt ) pd%pdos_mode(i,j)=0.0_flyt
            enddo
            pd%pdos_mode(1,j)=0.0_flyt
        enddo
        do j=1,pd%na
            f1=pd%pdos_site(1,j)
            do i=1,pd%ndos
                f0=(pd%ndos-i)*1.0_flyt/( (pd%ndos-1)*1.0_flyt )
                f0=f0**2
                pd%pdos_site(i,j)=pd%pdos_site(i,j)-f0*f1
                if ( pd%pdos_site(i,j) .lt. 0.0_flyt ) pd%pdos_site(i,j)=0.0_flyt
            enddo
            pd%pdos_site(1,j)=0.0_flyt
        enddo

        ! Sum up contributions and normalize things
        pd%dos=0.0_flyt
        do j=1,pd%nb
            f0=1.0_flyt/lo_trapezoid_integration(pd%omega,pd%pdos_mode(:,j))
            pd%dos=pd%dos+pd%pdos_mode(:,j)
            pd%pdos_mode(:,j)=pd%pdos_mode(:,j)*f0
        enddo
        f0=lo_trapezoid_integration(pd%omega,pd%dos)
        pd%dos=pd%dos*dr%nb/f0

        ! Adjust the mode projected so that they sum up to the total
        do i=1,pd%ndos
            if ( pd%dos(i) .gt. lo_tol/lo_twopi/1E12_flyt ) then
                f0=sum(pd%pdos_mode(i,:))
                pd%pdos_mode(i,:)=pd%pdos_mode(i,:)*pd%dos(i)/f0
            endif
        enddo

        ! Fix the degeneracy of the site-projected
        lo_allocate(sitebuf(pd%ndos,uc%na))
        sitebuf=0.0_flyt
        do i=1,pd%na
            ! enfore site degeneracy
            do j=1,uc%sym%degeneracy(i)
                k=uc%sym%degenerate_atom(j,i)
                sitebuf(:,i)=sitebuf(:,i)+pd%pdos_site(:,k)/(1.0_flyt*uc%sym%degeneracy(j))
            enddo
            f0=lo_trapezoid_integration(pd%omega,sitebuf(:,i))
            sitebuf(:,i)=sitebuf(:,i)*3.0_flyt/f0
        enddo
        pd%pdos_site=sitebuf
        lo_deallocate(sitebuf)
        ! normalize the projections so that they add up to the total
        do i=1,pd%ndos
            if ( pd%dos(i) .gt. lo_tol/lo_twopi/1E12_flyt ) then
                f0=sum(pd%pdos_site(i,:))
                pd%pdos_site(i,:)=pd%pdos_site(i,:)*pd%dos(i)/f0
            endif
        enddo
    end block cleandos

end subroutine
