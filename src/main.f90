#include "precompilerdefinitions"
!! main file for canonical configuration
program canonical_configuration
!!{!src/canonical_configuration/manual.md!}
use konstanter, only: flyt,lo_tol,lo_kb_hartree,lo_bohr_to_A,lo_twopi,lo_frequency_Hartree_to_THz,lo_eV_to_Hartree
use helpers, only: lo_mpi_helper,lo_random_int,tochar,lo_seed_random_numbers,walltime
use options, only: lo_opts
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_jij_secondorder, only: lo_jij_secondorder

implicit none
type(lo_opts) :: opts
type(lo_crystalstructure) :: ss,uc
type(lo_forceconstant_secondorder) :: fc,fcss
type(lo_jij_secondorder) :: jij 
type(lo_mpi_helper) :: mw

! Get the necessary things first
init: block
    ! Get CLI options
    call opts%parse()
    call mw%init()
    ! Seed the random numbers
    call lo_seed_random_numbers()

    ! Read structures
    write(*,*) '... reading infiles'
    call ss%readfromfile('infile.ssposcar',verbosity=opts%verbosity)
    call uc%readfromfile('infile.ucposcar',verbosity=opts%verbosity)
    ! Match the supercell to the unitcell, always a good idea
    call uc%classify('spacegroup',timereversal=.true.)
    call ss%classify('supercell',uc)

    ! get a forceconstant somehow
    if ( opts%debye_temperature .gt. 0.0_flyt ) then
        call fc%fake_forceconstant(uc,ss,debye_temperature=opts%debye_temperature,verbosity=opts%verbosity)
        write(*,*) '... constructed fake forceconstant corresponding to Td = ',tochar(opts%debye_temperature),'K'
        call fc%writetofile(uc,'outfile.fakeforceconstant')
        write(*,*) '... wrote it to "outfile.fakeforceconstant", check that the frequency range is reasonable'
    elseif ( opts%maximum_frequency .gt. 0.0_flyt ) then
        call fc%fake_forceconstant(uc,ss,maximum_frequency=opts%maximum_frequency,verbosity=opts%verbosity)
        write(*,*) '... constructed fake forceconstant corresponding to max(omega) = ',tochar(opts%maximum_frequency*lo_frequency_Hartree_to_THz),' THz'
        call fc%writetofile(uc,'outfile.fakeforceconstant')
        write(*,*) '... wrote it to "outfile.fakeforceconstant", check that dispersions look reasonable'
    else
        ! read the forceconstant from file
        call fc%readfromfile(uc,'infile.forceconstant')
    endif

    ! Create fake magnetic exchange interactions
    if ( abs(opts%exchange_J) .gt. lo_tol ) then
        ! Convert from the meV in the input to Hartree
        opts%exchange_J=opts%exchange_J*lo_eV_to_Hartree/1000
        opts%exchange_J=opts%exchange_J/(opts%mean_moment**2)
        call jij%fake_jij(uc,opts%exchange_J,opts%mean_moment)
        call jij%writetofile(uc,'outfile.fakejij')
    endif

    ! Remap force constant to supercell
    if ( uc%info%alloy ) then
        write(*,*) '... alloy detected'
        write(*,*) 'Not done'
        stop
        ! the unitcell is an alloy, now I have to think a little
        !call sqs%readfromfile('infile.sqs')
        !call alloy%readfromfile('infile.sqs_alloy_supercell')
        ! the remapping is a bit different
        !alloy%r=sqs%r
        !call fc%remap(uc,alloy,fcss)
    else
        ! in normal case just remap it
        call fc%remap(uc,ss,fcss)
    endif
    write(*,*) '... remapped fc'
end block init

write(*,*) '         ek(K)            ep(K)            <ek/ep>           T(K)            <T>(K)       <msd>(A)'
dumpconf: block
    type(lo_crystalstructure) :: p
    real(flyt), dimension(:,:,:,:), allocatable :: polar_fc
    real(flyt) :: ep,ek,temp,avgtemp,avgmsd,msd,ratio,rek,rep,f0
    integer :: i,j,a1,a2
    character(len=1000) :: fname

    ! Clean copy to work with
    p=ss
    ! Get a copy of the polar forceconstant to evaluate potential energy
    if ( fc%polar ) then
        lo_allocate(polar_fc(3,3,p%na,p%na))
        polar_fc=0.0_flyt
        call fc%supercell_longrange_dynamical_matrix_at_gamma(ss,polar_fc,1E-10_flyt)
    endif
    
    ratio=0.0_flyt
    avgtemp=0.0_flyt
    avgmsd=0.0_flyt
    rek=0.0_flyt
    rep=0.0_flyt
    do i=1,opts%nconf
        ! reset the structure
        p%r=ss%r
        p%rcart=ss%rcart
        p%v=0.0_flyt
        p%u=0.0_flyt
        p%f=0.0_flyt
        ! initialize
        call fcss%initialize_cell(p,uc,fc,opts%temperature,opts%zpm,.false.,opts%threshold)

        ! dump to file
        select case(opts%output_format)
        case(1) ! vasp output
            fname='contcar_conf'//tochar(i,4)
            call p%writetofile(trim(fname),opts%output_format,write_velocities=.true.)
        case(2) ! abinit output
            fname='abinput_conf'//tochar(i,4)
            call p%writetofile(trim(fname),opts%output_format,write_velocities=.true.)
        case(3) ! LAMMPS output
            fname='lammps_conf'//tochar(i,4)
            call p%writetofile(trim(fname),opts%output_format,write_velocities=.true.)
        case(4) ! AIMS output
            fname='aims_conf'//tochar(i,4)
            call p%writetofile(trim(fname),opts%output_format,write_velocities=.true.)
        end select

        ! just measure some stuff, for no good reason
        ek=p%kinetic_energy()/(p%na)
        ep=fcss%potential_energy(p%u)/(p%na)
        if ( fc%polar ) then
            f0=0.0_flyt
            do a1=1,p%na
            do a2=1,p%na
                f0=f0+dot_product(matmul(p%u(:,a1),polar_fc(:,:,a1,a2)),p%u(:,a2))*0.5_flyt
            enddo
            enddo
            ep=ep+f0/p%na
        endif

        rek=rek+ek
        rep=rep+ep
        ratio=rek/rep
        temp=(ek+ep)/(3*lo_kb_hartree)
        msd=0
        do j=1,p%na
            msd=msd+norm2(p%u(:,j))*lo_bohr_to_A/p%na
        enddo
        avgmsd=avgmsd+msd
        avgtemp=avgtemp+temp
        write(*,"(1X,5(2X,F15.5),2X,F12.8)") ek/lo_kb_hartree/1.5_flyt,ep/lo_kb_hartree/1.5_flyt,ratio,temp,avgtemp/i,avgmsd/i
    enddo
end block dumpconf

call mw%destroy()


end program
