#include "precompilerdefinitions"
module magneticdisorder
    !
    use constants
    use helpers
    use type_crystalstructure
    use type_symmetryoperation
    use type_symmetrytable
    !
    implicit none
    private
    public :: lo_magdisorder

    type lo_magdisorder_shell
        ! How many pairs in the supercell map to this shell?
        integer :: npair
        ! How should it be weighted?
        real(flyt) :: weight
        ! Radius of this coordination shell
        real(flyt) :: rad
        ! indices to the first atom in the pair
        integer, dimension(:), allocatable :: i1
        ! and the second atom
        integer, dimension(:), allocatable :: i2
    end type

    type lo_magdisorder
        !> collinear or noncollinear
        logical :: coll
        !> How many bins of different levels of disorder do we want?
        integer :: nbin
        !> How many configurations per bin
        integer :: nconf
        !> number of magnetic coordination shells
        integer :: nshell
        !> info about the coordination shells
        type(lo_magdisorder_shell), dimension(:), allocatable :: sh
        !> history of configurations
        integer, dimension(:,:,:), allocatable :: collhistory
        real(flyt), dimension(:,:,:,:), allocatable :: noncollhistory
        !> initial configuration
        integer, dimension(:), allocatable :: initial_collinear_configuration
        real(flyt), dimension(:,:), allocatable :: initial_noncollinear_configuration
        !> which sites are switchable?
        integer, dimension(:), allocatable :: sites
        !> initial correlation function
        real(flyt) :: maxcf        
        contains
            !> create the structure
            procedure :: generate    
            !> get the correlation function
            procedure :: correlation_function
            !> generate magnetic sqs
            procedure :: optimize
            !> dump to file
            procedure :: dump_configurations
    end type
    
contains

function random_unit_vector() result(v)
    real(flyt), dimension(3) :: v
    !
    integer :: i
    do i=1,3
        v(i)=lo_random_gaussian_number(0.0_flyt,1.0_flyt)
    enddo
    v=v/norm2(v)
end function

subroutine zerosum(x,rel)
    real(flyt), dimension(:,:), intent(inout) :: x
    logical, dimension(:), intent(in) :: rel
    !
    real(flyt), dimension(3) :: v
    real(flyt) :: f0
    integer :: i,j,n

    v=0.0_flyt
    n=size(x,2)
    j=0
    do i=1,n
        if ( rel(i) ) then
            v=v+x(:,i)
            j=j+1
        endif
    enddo
    v=v/(j*1.0_flyt)
    do i=1,n
        if ( rel(i) ) then
            x(:,i)=x(:,i)-v
            x(:,i)=x(:,i)/norm2(x(:,i))
        endif
    enddo
end subroutine

subroutine dump_configurations(mag)
    !> shells and stuff
    class(lo_magdisorder), intent(in) :: mag
    !
    integer :: i,j,k,l,m,u,ii,jj
    integer :: ncf
    character(len=10000) :: dum,fn
    integer, dimension(mag%nconf) :: ind
    real(flyt) :: f0

!    do l=1,mag%nbin
!    write(*,*) ' '
!    do i=1,mag%nconf
!        do j=1,mag%nconf
!            ind(j)=sum(abs(mag%history(:,j,l)-mag%history(:,i,l)))
!        enddo
!        write(*,'(40(1X,I2))') ind
!    enddo
!    enddo

!    ! I arranged stuff in strange order. Fix that.
!    ind(1)=1
!    do i=1,mag%nbin
!        ind(i+1)=mag%nbin-i+1
!    enddo
!    !

    
    ncf=mag%nbin+1
    do i=1,ncf
        fn='outfile.magmom_'//trim(int2char(i))
        u=open_file('out',trim(fn))
            do k=1,mag%nconf
                dum="MAGMOM = "
                if ( mag%coll ) then
                    do l=1,size(mag%collhistory,1)
                        if ( i .eq. 1 ) then
                            ii=mag%initial_collinear_configuration(l)
                        else
                            ii=mag%collhistory(l,k,i-1)
                        endif
                        ii=ii*3
                        dum=trim(dum)//" "//trim(int2char(ii))
                    enddo
                else
                    do l=1,size(mag%noncollhistory,2)
                    do m=1,3
                        if ( i .eq. 1 ) then
                            f0=mag%initial_noncollinear_configuration(m,l)
                        else
                            f0=mag%noncollhistory(m,l,k,i-1)
                        endif
                        f0=f0*3
                        dum=trim(dum)//" "//trim(flyt2char(f0))
                    enddo
                    enddo
                endif
                write(u,'(1X,A)') trim(dum)
            enddo
        close(u)
    enddo

end subroutine

!> calculate the correlation function
subroutine correlation_function(mag,collconf,noncollconf,cf)
    !> list of shells and stuff
    class(lo_magdisorder), intent(in) :: mag
    !> current collinear magnetic configuration
    integer, dimension(:), intent(in), optional :: collconf 
    !> current noncollinear magnetic configuration
    real(flyt), dimension(:,:), intent(in), optional :: noncollconf
    !> the correlation function per shell
    real(flyt), dimension(:), intent(out) :: cf
    !
    integer :: i,j,k,l,i1,i2
    real(flyt) :: f0
    logical :: coll
    !
    if ( mag%coll ) then
        cf=0.0_flyt
        do i=1,mag%nshell
            l=0
            do j=1,mag%sh(i)%npair
                i1=mag%sh(i)%i1(j)
                i2=mag%sh(i)%i2(j)
                l=l+collconf(i1)*collconf(i2)
            enddo
            cf(i)=(l*1.0_flyt)/(mag%sh(i)%npair*1.0_flyt)
        enddo
    else
        cf=0.0_flyt
        do i=1,mag%nshell
            f0=0.0_flyt
            do j=1,mag%sh(i)%npair
                i1=mag%sh(i)%i1(j)
                i2=mag%sh(i)%i2(j)
                f0=f0+dot_product(noncollconf(:,i1),noncollconf(:,i2))
            enddo
            cf(i)=(f0*1.0_flyt)/(mag%sh(i)%npair*1.0_flyt)
        enddo
    endif
end subroutine

!> generate optimized configuration
subroutine optimize(mag,ss,nconf,nbin)
    !> shells and stuff
    class(lo_magdisorder), intent(inout) :: mag
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: ss
    !> how many configurations do I want in each bin?
    integer, intent(in) :: nconf
    !> how many bins?
    integer, intent(in) :: nbin

    real(flyt), dimension(mag%nshell) :: cf
    real(flyt), dimension(:,:), allocatable :: cftargets
    integer :: bin,i,j,k

    call lo_seed_random_numbers()

    ! Space for the history
    mag%nbin=nbin
    mag%nconf=nconf
    if ( mag%coll ) then
        lo_allocate(mag%collhistory(ss%na,nconf,nbin))
        mag%collhistory=0
    else
        lo_allocate(mag%noncollhistory(3,ss%na,nconf,nbin))
        mag%noncollhistory=0.0_flyt
    endif

    ! Set the target correlation functions
    lo_allocate(cftargets(mag%nshell,mag%nbin))
    call mag%correlation_function(mag%initial_collinear_configuration,mag%initial_noncollinear_configuration,cf)
    cftargets=0.0_flyt
    do i=1,size(cftargets,1)
        cftargets(:,i)=cf*abs(( 1.0_flyt-(i)/(1.0_flyt*mag%nbin) )*mag%maxcf)
    enddo

    ! I have a series of correlation function targets, one minimization for each:
    do bin=1,mag%nbin

    findonetarget: block
        real(flyt), dimension(:), allocatable :: convcheck
        real(flyt) :: trg,cf0,cf1,f0,f1
        real(flyt) :: temperature,tempinc,tempdec,breaktol 
        integer, dimension(ss%na) :: conf0,conf1
        real(flyt), dimension(3,ss%na) :: ncconf0,ncconf1
        integer, dimension(:,:), allocatable :: collhist
        real(flyt), dimension(:,:,:), allocatable :: noncollhist
        integer :: nouter,ninner,na,histcounter
        integer :: outiter,initer,nflip
        
        ! Initial configuration
        if ( mag%coll ) then
            conf0=mag%initial_collinear_configuration
            conf1=0
            call mag%correlation_function(collconf=conf0,cf=cf)
        else
            ncconf0=mag%initial_noncollinear_configuration
            ncconf1=0
            call mag%correlation_function(noncollconf=ncconf0,cf=cf)
        endif
        cf0=sum(abs(cf-cftargets(:,bin)))/mag%nshell !maxcf
        cf1=0.0_flyt
        na=size(mag%sites,1)
        ! Some counters for the minimization
        nouter=100
        ninner=mag%nconf*ss%na
        temperature=0.3_flyt
        tempinc=1.5_flyt
        tempdec=0.5_flyt
        lo_allocate(convcheck(nouter))
        convcheck=0.0_flyt
        breaktol=1E-2_flyt

        write(*,*) 'Simulated annealing '//trim(int2char(bin))//', target: '!//trim(flyt2char(trg))

        ! Start minimizing
        outerloop1: do outiter=1,nouter
            nflip=0
            do initer=1,ninner
                ! Flip a spin, get new correlation function
                if ( mag%coll ) then
                    conf1=conf0
                    do
                        i=lo_random_int(na)
                        j=lo_random_int(na)
                        if ( conf1(mag%sites(i))*conf1(mag%sites(j)) .lt. 0 ) exit
                    enddo
                    conf1(mag%sites(i))=-1*conf1(mag%sites(i))
                    conf1(mag%sites(j))=-1*conf1(mag%sites(j))
                    call mag%correlation_function(collconf=conf1,cf=cf)
                else
                    ncconf1=ncconf0
                    i=lo_random_int(na)
                    ncconf1(:,mag%sites(i))=random_unit_vector()
                    call zerosum(ncconf1,ss%mag%atom_has_moment)
                    call mag%correlation_function(noncollconf=ncconf1,cf=cf)
                endif
                cf1=sum(abs(cf-cftargets(:,bin)))/mag%nshell
                ! MC compare thingy
                f0=exp(-(cf1-cf0)/temperature)
                call random_number(f1)
                ! keep?
                if ( f0 .gt. f1 ) then
                    if ( mag%coll ) then
                        conf0=conf1
                    else
                        ncconf0=ncconf1
                    endif
                    cf0=cf1
                    nflip=nflip+1
                endif 
            enddo
            ! check how many flips there were, and maybe adjust the temperature
            f0=(1.0_flyt*nflip)/(1.0_flyt*ninner)        
            if ( f0 .lt. 0.02_flyt ) then
                temperature=temperature*tempinc
            elseif ( f0 .gt. 0.10_flyt ) then
                temperature=temperature*tempdec
            endif
            ! check for convergence
            convcheck(outiter)=cf0
            if ( outiter .ge. 5 ) then
                f1=lo_mean(convcheck(outiter-4:outiter))
            else
                f1=1.0_flyt
            endif            
            if ( f1 .lt. breaktol ) exit outerloop1
            !
            write(*,'(1X,I8,1X,F10.7,3(1X,F10.7))') outiter,cf0,temperature,real(f0),f1
            !
        enddo outerloop1

        convcheck=0.0_flyt
        histcounter=0
        if ( mag%coll ) then
            lo_allocate(collhist(ss%na,mag%nconf))
            collhist=0
        else
            lo_allocate(noncollhist(3,ss%na,mag%nconf))
            noncollhist=0.0_flyt
        endif
        ! Now gather some statistics for the history
        outerloop2: do outiter=1,nouter
            !
            nflip=0
            do initer=1,ninner
                ! Flip a spin, get new correlation function
                if ( mag%coll ) then
                    conf1=conf0
                    do
                        i=lo_random_int(na)
                        j=lo_random_int(na)
                        if ( conf1(mag%sites(i))*conf1(mag%sites(j)) .lt. 0 ) exit
                    enddo
                    conf1(mag%sites(i))=-1*conf1(mag%sites(i))
                    conf1(mag%sites(j))=-1*conf1(mag%sites(j))
                    call mag%correlation_function(collconf=conf1,cf=cf)
                else
                    ncconf1=ncconf0
                    i=lo_random_int(na)
                    ncconf1(:,mag%sites(i))=random_unit_vector()
                    call zerosum(ncconf1,ss%mag%atom_has_moment)
                    call mag%correlation_function(noncollconf=ncconf1,cf=cf)
                endif
                cf1=sum(abs(cf-cftargets(:,bin)))/mag%nshell
                ! MC compare thingy
                f0=exp(-(cf1-cf0)/temperature)
                call random_number(f1)
                ! keep?
                if ( f0 .gt. f1 ) then
                    if ( mag%coll ) then
                        conf0=conf1
                        ! Is this a new configuration?
                        j=0
                        do i=1,histcounter
                            if ( sum(abs(conf1-collhist(:,i))) .le. 2 ) then
                                j=j+1
                                exit
                            endif
                        enddo
                        ! If it is new, keep it!
                        if ( j .eq. 0 ) then
                            histcounter=histcounter+1
                            collhist(:,histcounter)=conf1
                            ! We might be done now!
                            if ( histcounter .eq. mag%nconf ) exit outerloop2
                        endif
                    else
                        ncconf0=ncconf1
                        ! Is this a new configuration?
                        j=0
                        do i=1,histcounter
                            if ( sum(abs(ncconf1-noncollhist(:,:,i))) .le. 0.2_flyt ) then
                                j=j+1
                                exit
                            endif
                        enddo
                        ! If it is new, keep it!
                        if ( j .eq. 0 ) then
                            histcounter=histcounter+1
                            noncollhist(:,:,histcounter)=ncconf1
                            ! We might be done now!
                            if ( histcounter .eq. mag%nconf ) exit outerloop2
                        endif
                    endif
                    !
                    cf0=cf1
                    nflip=nflip+1
                    !
                endif 
            enddo
            ! check how many flips there were, and maybe adjust the temperature
            f0=(1.0_flyt*nflip)/(1.0_flyt*ninner)        
            if ( f0 .lt. 0.02_flyt ) then
                temperature=temperature*tempinc
            elseif ( f0 .gt. 0.10_flyt ) then
                temperature=temperature*tempdec
            endif
        enddo outerloop2
        ! And store the history
        if ( mag%coll ) then
            mag%collhistory(:,:,bin)=collhist
        else
            mag%noncollhistory(:,:,:,bin)=noncollhist
        endif
    end block findonetarget
    enddo
    !
end subroutine

!> set up all the coordination shells and stuff
subroutine generate(mag,sy,uc,ss)
    !> to keep track of all the coordination shells
    class(lo_magdisorder), intent(out) :: mag
    !> a symmetry table, useful for a bunch of stuff
    type(lo_symcell), intent(in) :: sy
    !> unit cell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss

    ! Count unique atoms with magnetic moments
    findrelevantshells: block
        integer, dimension(:), allocatable :: dumi1,dumi2
        integer :: i,j,k,l
        integer, parameter ::magicfactor=1000000

        ! Count number of relevant coordination shells
        l=0
        do i=1,sy%nun
            do j=1,sy%un(i)%npair
                if ( uc%mag%atom_has_moment( sy%un(i)%pair(j)%i1 ) .and. uc%mag%atom_has_moment( sy%un(i)%pair(j)%i2 ) ) then
                    if ( lo_sqnorm(sy%un(i)%pair(j)%v2) .lt. lo_sqtol ) cycle
                    l=l+1
                endif
            enddo
        enddo
        lo_allocate(dumi1(l))
        l=0
        do i=1,sy%nun
            do j=1,sy%un(i)%npair
                if ( uc%mag%atom_has_moment( sy%un(i)%pair(j)%i1 ) .and. uc%mag%atom_has_moment( sy%un(i)%pair(j)%i2 ) ) then
                    if ( lo_sqnorm(sy%un(i)%pair(j)%v2) .lt. lo_sqtol ) cycle
                    l=l+1
                    ! magic factor 100000 is to differentiate between shells from different atoms
                    dumi1(l)=sy%un(i)%pair(j)%kvalisort_of_pair+i*magicfactor
                endif
            enddo
        enddo
        call lo_return_unique(dumi1,dumi2)
        mag%nshell=size(dumi2,1)
        lo_allocate(mag%sh(mag%nshell))
        do i=1,mag%nshell
            mag%sh(i)%npair=0
            mag%sh(i)%weight=0.0_flyt
        enddo

        ! Now count, from the supercell, where each pair ends up, if it is at all 
        do i=1,sy%nss
        do j=1,sy%ss(i)%npair
            if ( ss%mag%atom_has_moment( sy%ss(i)%pair(j)%i1 ) .and. ss%mag%atom_has_moment( sy%ss(i)%pair(j)%i2 ) ) then
                ! skip self-terms, irrelevant
                if ( lo_sqnorm(sy%ss(i)%pair(j)%v2) .lt. lo_sqtol ) cycle
                ! ensure i1<i2, to avoid doublecounting
                if ( sy%ss(i)%pair(j)%i1 .gt. sy%ss(i)%pair(j)%i2 ) cycle
                ! figure out where it should be added. Same magic factor as above.
                k=sy%ss(i)%pair(j)%kvalisort_of_pair+sy%ss(i)%kvalisort*magicfactor
                do l=1,mag%nshell
                    if ( k .eq. dumi2(l) ) then
                        mag%sh(l)%npair=mag%sh(l)%npair+1
                        exit
                    endif
                enddo
            endif
        enddo
        enddo

        ! make some space
        do i=1,mag%nshell
            lo_allocate(mag%sh(i)%i1( mag%sh(i)%npair ))
            lo_allocate(mag%sh(i)%i2( mag%sh(i)%npair ))
            mag%sh(i)%i1=0
            mag%sh(i)%i2=0
            mag%sh(i)%npair=0
        enddo

        ! Now count, from the supercell, where each pair ends up, if it is at all 
        do i=1,sy%nss
        do j=1,sy%ss(i)%npair
            if ( ss%mag%atom_has_moment( sy%ss(i)%pair(j)%i1 ) .and. ss%mag%atom_has_moment( sy%ss(i)%pair(j)%i2 ) ) then
                ! skip self-terms, irrelevant
                if ( lo_sqnorm(sy%ss(i)%pair(j)%v2) .lt. lo_sqtol ) cycle
                ! ensure i1<i2, to avoid doublecounting
                if ( sy%ss(i)%pair(j)%i1 .gt. sy%ss(i)%pair(j)%i2 ) cycle
                ! figure out where it should be added
                k=sy%ss(i)%pair(j)%kvalisort_of_pair+sy%ss(i)%kvalisort*magicfactor
                do l=1,mag%nshell
                    if ( k .eq. dumi2(l) ) then
                        mag%sh(l)%npair=mag%sh(l)%npair+1
                        mag%sh(l)%i1( mag%sh(l)%npair )=sy%ss(i)%pair(j)%i1
                        mag%sh(l)%i2( mag%sh(l)%npair )=sy%ss(i)%pair(j)%i2
                        exit
                    endif
                enddo
            endif
        enddo
        enddo
    end block findrelevantshells
    write(*,*) '... identified magnetic coordination shells'

    ! Now set up the starting configuration
    setupconf: block
        real(flyt), dimension(:), allocatable :: cf
        integer :: i,j

        ! collinear or noncollinear?
        if ( uc%info%collmag ) then
            mag%coll=.true.
        else
            mag%coll=.false.
        endif

        ! Get a list of switchable sites
        j=0
        do i=1,ss%na
            if ( ss%mag%atom_has_moment(i) ) then
                j=j+1
            endif
        enddo
        lo_allocate(mag%sites(j))
        j=0
        do i=1,ss%na
            if ( ss%mag%atom_has_moment(i) ) then
                j=j+1
                mag%sites(j)=i
            endif
        enddo

        ! Set up the initial configuration
        if ( mag%coll ) then
            lo_allocate(mag%initial_collinear_configuration(ss%na))
            mag%initial_collinear_configuration=0        
            do i=1,ss%na
                if ( ss%mag%atom_has_moment(i) ) then
                    mag%initial_collinear_configuration(i)=int(anint(ss%mag%collinear_moment(i)/abs(ss%mag%collinear_moment(i))))
                endif
            enddo
        else
            lo_allocate(mag%initial_noncollinear_configuration(3,ss%na))
            mag%initial_noncollinear_configuration=0        
            do i=1,ss%na
                if ( ss%mag%atom_has_moment(i) ) then
                    mag%initial_noncollinear_configuration(:,i)=ss%mag%noncollinear_moment(:,i)/norm2(ss%mag%noncollinear_moment(:,i))
                endif
            enddo
        endif

        ! Get the reference correlation function
        lo_allocate(cf(mag%nshell))
        if ( mag%coll ) then
            call mag%correlation_function(collconf=mag%initial_collinear_configuration,cf=cf)
        else
            call mag%correlation_function(noncollconf=mag%initial_noncollinear_configuration,cf=cf)
        endif
        mag%maxcf=sum(abs(cf))
        
    end block setupconf

end subroutine

end module
