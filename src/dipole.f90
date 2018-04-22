#include "precompilerdefinitions"
module dipole
use konstanter, only: flyt,lo_sqtol,lo_tol
use helpers !, only: lo_frobnorm,walltime,lo_trueNtimes,lo_progressbar_init,lo_progressbar
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forcemap, only: lo_forcemap
use type_mdsim, only: lo_mdsim
use quadratures_stencils
implicit none
private
public :: subtract_longrange_interactions

public :: histfit

contains

subroutine histfit(sim)
    !> force-displacement data
    type(lo_mdsim), intent(inout) :: sim
   
    integer, parameter :: nhist=100
    real(flyt), dimension(:), allocatable :: moments
    real(flyt), parameter :: sigmafactor=0.05_flyt
    real(flyt), dimension(nhist) :: hx,hy,cdfy,cdfx

    ! First build the histogram. Not sure how to do it smoothly, but I'll figure it out.
    hist: block
        real(flyt) :: xmin,xmax,sigma,foursigma,stddev,invf
        real(flyt) :: invh
        integer :: i,j,k,l,ii,jj,u

        ! Grab all the moments
        lo_allocate(moments(sim%nt*sim%na))
        moments=0.0_flyt
        l=0
        do i=1,sim%nt
        do j=1,sim%na
            l=l+1
            moments(l)=norm2(sim%m(:,j,i))
        enddo
        enddo

        ! Build the cumulative distribution function
        xmin=minval(moments)*0.97_flyt
        xmax=maxval(moments)*1.03_flyt
        call lo_linspace(xmin,xmax,cdfx)
        do i=1,nhist
            cdfy(i)=sum(moments,moments<cdfx(i))
        enddo
        cdfy=cdfy/cdfy(nhist)

        ! Get the pdf from this.
        hx=cdfx
        hy=0.0_flyt
        invh=1.0_flyt/(hx(2)-hx(1))
        ! get the edges first
        hy(1)=(cdfy(2)-cdfy(1))*invh
        hy(nhist)=( cdfy(nhist)-cdfy(nhist-1) )*invh
        do i=2,nhist-1
            hy(i)= (cdfy(i+1)-cdfy(i-1))*0.5_flyt*invh
        enddo
    end block hist

    ! Now try to fit it to something
    fithist: block
        real(flyt) :: nrm0,nrm1
        real(flyt) :: f0,beta1,beta2,alpha
        real(flyt) :: mu0,sig0,e0,e1,e2
        real(flyt) :: step,a,b,c,d,par_tk,par_theta
        real(flyt), dimension(3) :: y
        real(flyt), dimension(2) :: par0,par1,grad0,grad1,dir,searchdir,bfgs_y,bfgs_s
        real(flyt), dimension(2,2) :: hessian,bfgs_B,bfgs_C,invhessian
        integer :: i,j,k,l,u
        integer :: iter,outer,ctr
        logical :: damped_update

        ! Starting guess for the mean
        mu0=lo_trapezoid_integration(hx,hx*hy)
        ! Starting guess for sigma. Hmm.
        e1=sqrt(lo_trapezoid_integration(hx,hy*(hx-mu0)**2))
        e1=e1*2*Sqrt(2*log(2.0_flyt))
        sig0=-4*Log(2.0_flyt)/( e1**4-4*e1**2*mu0**2)

        ! Get the starting point        
        par0=[mu0,sig0]
        !call cdferr(cdfx,cdfy,par0,grad0,e0)
        e0=lsqerr(hx,hy,par0)
        call lo_identitymatrix(hessian)
        hessian=hessian*1E6_flyt
        par1=par0
        grad1=grad0

        ! First find decent sigma to start from
        step=0.001_flyt
        e1=lo_huge

        write(*,*) 0,par0,e0
        do outer=1,1 !200 !1 !20
            grad0=-lsqgrad(hx,hy,par0)
            step=1E-5_flyt
            ctr=0
            do iter=1,1000 !00
                grad0=-lsqgrad(hx,hy,par0)
                par0=par0+grad0*step
                e1=lsqerr(hx,hy,par0)
                if ( e1 .lt. e0 ) then
                    e0=e1
                    step=step*2_flyt
                    ctr=0
                else
                    par0=par0-grad0*step
                    step=step*0.25_flyt
                    ctr=ctr+1
                    if ( ctr .gt. 5 ) exit
                endif
                write(*,*) iter,par0,e1,step
            enddo
            write(*,*) outer,par0,e0,norm2(grad0)
        enddo

!        do iter=1,500 !-1 !200
!            !call cdferr(cdfx,cdfy,par1,grad1,e1
!            call lsqerr(hx,hy,par1,grad1,e1)
!            if ( iter .gt. 1 ) then
!                bfgs_y = grad1-grad0 
!                bfgs_s = par1-par0
!                bfgs_B = hessian
!                ! Slight damping in the Hessian update ...
!                beta1=dot_product(bfgs_y,bfgs_s)
!                beta2=dot_product(bfgs_s,matmul(bfgs_B,bfgs_s))
!                par_tk=1.2_flyt*beta2
!                if ( dot_product(bfgs_s,bfgs_y) .lt. par_tk ) then
!                    par_theta=0.8_flyt*beta2/(beta2+dot_product(bfgs_s,bfgs_y))
!                    damped_update=.true.
!                else
!                    par_theta=1.0_flyt
!                    damped_update=.false.
!                endif
!
!                bfgs_y=bfgs_y*par_theta+(1.0_flyt-par_theta)*matmul(bfgs_B,bfgs_s)
!                beta1=dot_product(bfgs_y,bfgs_s)
!                beta2=dot_product(bfgs_s,matmul(bfgs_B,bfgs_s))
!                bfgs_s=matmul(bfgs_B,bfgs_s)
!                do i=1,2
!                do j=1,2
!                    bfgs_C(j,i)=bfgs_B(j,i)+bfgs_y(i)*bfgs_y(j)/beta1-bfgs_s(i)*bfgs_s(j)/beta2
!                enddo
!                enddo
!                hessian=bfgs_C*0.2_flyt+hessian*0.8_flyt
!
!                ! Check eigenvalues of hessian
!                a=hessian(1,1)
!                b=hessian(2,1)
!                c=hessian(1,2)
!                d=hessian(2,2)
!                beta1=(a + d - Sqrt(a**2 + 4*b*c - 2*a*d + d**2))*0.5_flyt
!                beta2=(a + d + Sqrt(a**2 + 4*b*c - 2*a*d + d**2))*0.5_flyt
!                a=min(beta1,beta2)
!                if ( a .lt. lo_tol ) then
!                    hessian(1,1)=hessian(1,1)+abs(a)+1.0_flyt
!                    hessian(2,2)=hessian(2,2)+abs(a)+1.0_flyt
!                endif
!!write(*,*) 'bb',beta1,beta2
!            endif
!            e0=e1
!            grad0=grad1
!            par0=par1
!            ! Invert the hessian
!            a=hessian(1,1)
!            b=hessian(2,1)
!            c=hessian(1,2)
!            d=hessian(2,2)
!            invhessian(:,1)=[d,-b]
!            invhessian(:,2)=[-c,-a]
!            invhessian=invhessian/(a*d-b*c)
!            ! New gradient
!
!            if ( damped_update ) then
!                searchdir=grad1/norm2(grad1)
!            else
!                searchdir=-matmul(invhessian,grad1)
!                searchdir=searchdir/norm2(searchdir)
!            endif
!            ! How far?
!!write(*,*) hessian
!
!            beta1=dot_product(grad1,searchdir)
!            beta2=dot_product(searchdir,matmul(hessian,searchdir))
!
!            if ( beta2 .gt. lo_sqtol ) then
!                alpha=abs(beta1/beta2)
!if ( abs(alpha) .lt. 1E-8_flyt ) exit
!            else
!                alpha=abs(beta1)
!            endif
!            par1=par1+alpha*searchdir*1.0_flyt
!            !
!            write(*,*) iter,e0,alpha,norm2(grad1),damped_update
!            !
!        enddo
        
u=open_file('out','dumhist')
    e0=normfactor(par0)
    do i=1,nhist
        write(u,*) hx(i),hy(i),distr(hx(i),par0(1),par0(2))*e0
    enddo
close(u)

u=open_file('out','dumcdf')
    nrm0=normfactor(par0)
    do i=1,nhist
        y=cdffn(par0,nrm0,cdfx(i))

        write(u,*) cdfx(i),cdfy(i),y(1)
    enddo
close(u)

!par0=[1.0_flyt,0.25_flyt]

    do i=3,50
!        e0=1.0_flyt/normfactor(par0,20)
!        e1=cdffn(par0,i,e0,2.5_flyt)
!        write(*,*) i,e0,e1,e1-cdffn(par0,40,e0,2.5_flyt)
    enddo
    end block fithist

stop

contains

function lsqgrad(hx,hy,par) result(grad)
    real(flyt), dimension(:), intent(in) :: hx,hy
    real(flyt), dimension(2), intent(in) :: par
    real(flyt), dimension(2) :: grad

    integer, parameter :: n=6,m=2*n+1
    real(flyt), dimension(m,4) :: sc
    real(flyt), dimension(m) :: fv
    integer :: i

    call lo_centraldifference(n,par(1),1E-5_flyt,sc)
    do i=1,m
        fv(i)=lsqerr(hx,hy,[sc(i,1),par(2)])
    enddo
    grad(1)=sum(sc(:,2)*fv)

    call lo_centraldifference(n,par(2),1E-5_flyt,sc)
    do i=1,m
        fv(i)=lsqerr(hx,hy,[par(1),sc(i,1)])
    enddo
    grad(2)=sum(sc(:,2)*fv)
end function

function lsqerr(hx,hy,par) result(errorsq)
    real(flyt), dimension(:), intent(in) :: hx,hy
    real(flyt), dimension(2), intent(in) :: par
    real(flyt) :: errorsq

    real(flyt) :: mu,sig,gmu,gsig,nrm
    real(flyt) :: f0,f1,f2
    integer :: i

    nrm=normfactor(par)

    ! Squared error and gradient
    mu=par(1)
    sig=par(2)
    errorsq=0.0_flyt
    do i=1,nhist
        f0=distr(hx(i),mu,sig)*nrm
        errorsq=errorsq+( hy(i) - f0 )**2
    enddo        
end function 

function distr(x,mu,sig) result(y)
    real(flyt) :: x,mu,sig
    real(flyt) :: y

    real(flyt) :: arg
    arg=-sig*(x**2-mu**2)**2
    if ( arg .lt. -40.0_flyt ) then
        y=0.0_flyt
    else
        y=exp(arg)
    endif
end function

function energy(x,mu,sig) result(e)
    real(flyt) :: x,mu,sig
    real(flyt) :: e
    e=sig*( x**2-mu**2)**2
end function

!> normalization factor for my distribution
function normfactor(par) result(nrm)
    real(flyt), dimension(2), intent(in) :: par
    real(flyt) :: nrm

    integer, parameter :: nquad=30
    real(flyt), dimension(2,nquad) :: gq
    real(flyt) :: x0,x1,x2,x3,z1,z2,z3,w1,w2,w3
    integer :: i

    ! Get the quadrature weights
    call lo_gaussianquadrature(nquad,0.0_flyt,1.0_flyt,gq)

    x0=0.0_flyt
    z1=sqrt(log(2.0_flyt)/par(2))
    if ( par(1)**2 .gt. z1 ) then
        x1=sqrt(par(1)**2-z1)
    else
        x1=par(1)*0.5_flyt
    endif
    x2=sqrt(par(1)**2+z1)
    x3=sqrt(par(1)**2+sqrt(log(1E14_flyt)/par(2)))

    nrm=0.0_flyt
    do i=1,nquad
        z1=gq(1,i)*(x1-x0)+x0
        z2=gq(1,i)*(x2-x1)+x1
        z3=gq(1,i)*(x3-x2)+x2
        w1=gq(2,i)*(x1-x0)
        w2=gq(2,i)*(x2-x1)
        w3=gq(2,i)*(x3-x2)
        nrm=nrm+distr( z1,par(1),par(2) )*w1
        nrm=nrm+distr( z2,par(1),par(2) )*w2
        nrm=nrm+distr( z3,par(1),par(2) )*w3
    enddo
    nrm=1.0_flyt/nrm
end function

!> cumulative distribution function thing
function cdffn(par,normfactor,x) result(y)
    !> parameters for the distribution
    real(flyt), dimension(2), intent(in) :: par
    !> normalization factor
    real(flyt), intent(in) :: x
    !> value, mu-derivative,sigma-derivative
    real(flyt), dimension(3) :: y 

    integer, parameter :: nquad=50
    real(flyt), dimension(2,nquad) :: gq
    real(flyt) :: normfactor,f0
    integer :: i

    ! Get the quadrature weights
    call lo_gaussianquadrature(nquad,0.0_flyt,x,gq)
    y=0.0_flyt
    do i=1,nquad
        f0=distr(gq(1,i),par(1),par(2))*gq(2,i)
        y(1)=y(1)+f0
        y(2)=y(2)-f0*4*par(1)*(par(1)**2-gq(1,i)**2)*par(2)
        y(3)=y(3)-f0*(par(1)**2-gq(1,i)**2)**2
    enddo
    y=y*normfactor
end function
 
end subroutine

!> remove longrange interactions
subroutine subtract_longrange_interactions(sim,map,uc,ss)
    !> force-displacement data
    type(lo_mdsim), intent(inout) :: sim
    !> forcemap
    type(lo_forcemap), intent(in) :: map
    !> unitcell
    type(lo_crystalstructure), intent(inout) :: uc
    !> supercell
    type(lo_crystalstructure), intent(inout) :: ss

    type(lo_forceconstant_secondorder) :: fc
    real(flyt), dimension(:,:,:,:), allocatable :: forceconstant
    real(flyt), dimension(:,:), allocatable :: f
    real(flyt), dimension(3,3) :: m0,m1
    real(flyt), dimension(3) :: v0
    real(flyt) :: t0,f0
    integer :: a1,a2,t

    write(*,*) ''
    write(*,*) 'SUBTRACTING LONGRANGE ELECTROSTATIC INTERACTIONS'

    ! Get some forceconstant thingy
    call map%get_secondorder_forceconstant(uc,fc)
    lo_allocate(forceconstant(3,3,ss%na,ss%na))
    call fc%supercell_longrange_dynamical_matrix_at_gamma(ss,forceconstant,1E-15_flyt)

    ! small sanity-check, just in case
    f0=0.0_flyt
    do a1=1,ss%na
    do a2=a1,ss%na
        m0=forceconstant(:,:,a1,a2)
        m1=transpose(forceconstant(:,:,a2,a1))
        f0=f0+lo_frobnorm(m0-m1)
    enddo
    enddo
    f0=f0/(ss%na**2)
    if ( f0 .gt. lo_tol ) then
        write(*,*) 'ERROR, non-Hermitian electrostatic forceconstant: ',f0
        write(*,*) 'If this keeps happening I should not hard-code the Ewald tolerance.'
        stop
    endif

    t0=walltime()
    call lo_progressbar_init()
    ! Subtract forces
    lo_allocate(f(3,ss%na))
    do t=1,sim%nt
        ! calculate the long-range forces
        f=0.0_flyt
        do a1=1,ss%na
        do a2=1,ss%na
            v0=matmul(forceconstant(:,:,a1,a2),sim%u(:,a2,t))
            f(:,a1)=f(:,a1)-v0
        enddo
        enddo
        ! sanity check that they add up to zero
        if ( abs(sum(f)) .gt. lo_sqtol ) then
            write(*,*) ''
            write(*,*) 'ERROR: dipole-dipole forces do not add up to zero!'
            stop
        endif
        ! subtract them
        sim%f(:,:,t)=sim%f(:,:,t)-f
        if( lo_trueNtimes(t,50,sim%nt) ) call lo_progressbar(' ... subtracting dipole-dipole forces',t,sim%nt)
    enddo
    call lo_progressbar(' ... subtracting dipole-dipole forces',sim%nt,sim%nt,walltime()-t0)
end subroutine


end module
