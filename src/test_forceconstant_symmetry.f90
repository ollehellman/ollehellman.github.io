#include "precompilerdefinitions"
module test_forceconstant_symmetry
use konstanter, only: flyt,lo_sqtol,lo_status
use helpers, only: tochar,lo_sqnorm
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_firstorder, only: lo_forceconstant_firstorder
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_forcemap, only: lo_forcemap

implicit none
private
public :: test_forceconstant_constraints

contains

subroutine test_forceconstant_constraints(map,fc1,fc2,fc3,fc4,uc)
    !> all the forceconstant metadata
    type(lo_forcemap), intent(in) :: map
    !> firstorder forceconstant
    type(lo_forceconstant_firstorder), intent(in) :: fc1
    !> secondorder forceconstant
    type(lo_forceconstant_secondorder), intent(in) :: fc2
    !> thirdorder forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fc3
    !> fourthorder forceconstant
    type(lo_forceconstant_fourthorder), intent(in) :: fc4
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc

    ! Start testing a lot of sum rules on the forceconstants, start with the trivial
    ! rotational constraints on the first order forceconstant
    if ( map%firstorder ) then
    rot1: block
        real(flyt), dimension(3,3) :: m0,m1
        real(flyt), dimension(3) :: r
        integer :: a1,i,j

        ! Check the constraints
        m0=0.0_flyt
        m1=0.0_flyt
        do a1=1,uc%na
            r=uc%rcart(:,a1) 
            do i=1,3
            do j=1,3
                m0(i,j)=m0(i,j)+fc1%atom(a1)%m(i)*r(j)
                m1(i,j)=m1(i,j)+fc1%atom(a1)%m(j)*r(i)
            enddo
            enddo
        enddo
        write(*,*) '     firstorder rotational constraint:',abs(sum(m0-m1))
    end block rot1
    endif

    ! Then we have the first+second order
    if ( map%firstorder .and. map%secondorder ) then
    rot2: block
        real(flyt), dimension(3,3,3) :: m0,m1
        real(flyt), dimension(3) :: r
        integer :: a1,i,al,be,gm

        do a1=1,uc%na
            m0=0.0_flyt
            m1=0.0_flyt
            do i=1,fc2%atom(a1)%n
                r=fc2%atom(a1)%pair(i)%r
                do al=1,3
                do be=1,3
                do gm=1,3
                    m0(al,be,gm)=m0(al,be,gm)+fc2%atom(a1)%pair(i)%m(al,be)*r(gm)+fc1%atom(a1)%m(be)*kron(al,gm)
                    m1(al,be,gm)=m1(al,be,gm)+fc2%atom(a1)%pair(i)%m(al,gm)*r(be)+fc1%atom(a1)%m(gm)*kron(al,be)
                enddo
                enddo
                enddo
            enddo
            write(*,*) '   second+first rotational constraint:',sum(m0-m1),'(atom ',tochar(a1),')'
        enddo
    end block rot2
    endif

    ! Pure secondorder rotational
    if ( (map%firstorder .eqv. .false.) .and. map%secondorder ) then
    rot3: block
        real(flyt), dimension(3,3,3) :: m0,m1
        real(flyt), dimension(3) :: r
        integer :: a1,i,al,be,gm

        do a1=1,uc%na
            m0=0.0_flyt
            m1=0.0_flyt
            do i=1,fc2%atom(a1)%n
                r=fc2%atom(a1)%pair(i)%r
                do al=1,3
                do be=1,3
                do gm=1,3
                    m0(al,be,gm)=m0(al,be,gm)+fc2%atom(a1)%pair(i)%m(al,be)*r(gm)
                    m1(al,be,gm)=m1(al,be,gm)+fc2%atom(a1)%pair(i)%m(al,gm)*r(be)
                enddo
                enddo
                enddo
            enddo
            write(*,*) '    secondorder rotational constraint:',abs(sum(m0-m1))
        enddo
    end block rot3
    endif

    if ( map%secondorder .and. map%thirdorder ) then
    rot4: block
        real(flyt), dimension(3,3,3,3) :: m0,m1
        real(flyt), dimension(3) :: r
        real(flyt) :: f0
        integer :: a1,i,j,k,al,be,gm,lm,ii,ti
        integer, dimension(:), allocatable :: di

        do a1=1,uc%na
            lo_allocate( di(fc3%atom(a1)%n))
            di=0
            f0=0.0_flyt
            do i=1,fc2%atom(a1)%n
                r=fc2%atom(a1)%pair(i)%r
                ! get a list if triplets with the same j-vector, to sum over
                ii=0                
                do j=1,fc3%atom(a1)%n
                    if ( lo_sqnorm(r-fc3%atom(a1)%triplet(j)%rv3 ) .lt. lo_sqtol ) then
                        ii=ii+1
                        di(ii)=j
                    endif
                enddo
                ! Now loop and sum over these triplets
                m0=0.0_flyt
                m1=0.0_flyt
                do k=1,ii
                    ti=di(k)
                    do al=1,3
                    do be=1,3
                    do gm=1,3
                    do lm=1,3
                        m0(al,be,gm,lm)=m0(al,be,gm,lm)+&
                            fc2%atom(a1)%pair(i)%m(gm,be)*kron(al,lm)+&
                            fc2%atom(a1)%pair(i)%m(gm,be)*kron(al,lm)+&
                            fc3%atom(a1)%triplet(ti)%m(al,be,gm)*r(lm)
                        m1(al,be,gm,lm)=m1(al,be,gm,lm)+&
                            fc2%atom(a1)%pair(i)%m(lm,be)*kron(al,gm)+&
                            fc2%atom(a1)%pair(i)%m(lm,be)*kron(al,gm)+&
                            fc3%atom(a1)%triplet(ti)%m(al,be,lm)*r(gm)
                    enddo
                    enddo
                    enddo
                    enddo
                enddo
                f0=f0+abs(sum(m0-m1))
            enddo
            lo_deallocate(di)
            write(*,*) '     thirdorder rotational constraint:',a1,f0
        enddo
    end block rot4
    endif

    contains
    function kron(i,j) result(k)
        integer :: i,j,k
        if ( i .eq. j ) then
            k=1
        else
            k=0
        endif
    end function
end subroutine

end module
