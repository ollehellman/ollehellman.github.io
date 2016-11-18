#include "precompilerdefinitions"
!! main file for canonical configuration
program canonical_configuration
!!{!src/canonical_configuration/manual.md!}
use constants
use helpers
use options
use type_qpointmesh
use type_phonon_dispersions
use type_crystalstructure
use type_forceconstant_secondorder
use type_qpointmesh
use type_lotosplitting

implicit none
type(lo_opts) :: opts
type(lo_crystalstructure) :: ss,uc,sqs,alloy,p
type(lo_forceconstant_secondorder) :: fc,fcss
class(lo_qpoint_mesh), allocatable :: qp
type(lo_loto) :: loto
type(lo_phonon_dispersions) :: dr

integer :: i
integer, dimension(3) :: qdim
real(flyt) :: ep,ek,temp,avgtemp
real(flyt) :: f0,f1,alpha,step
character(len=1000) :: fname

! Get CLI options
call opts%parse()
! Seed the random numbers
call lo_seed_random_numbers()

! Read structures
write(*,*) '... reading infiles'
call ss%readfromfile('infile.ssposcar',verbosity=opts%verbosity)
call uc%readfromfile('infile.ucposcar',verbosity=opts%verbosity)
! Match the supercell to the unitcell, always a good idea
call ss%classify('supercell',uc)

! get a forceconstant somehow
if ( opts%debye_temperature .gt. 0.0_flyt ) then
    call fc%fake_forceconstant(uc,ss,opts%debye_temperature,verbosity=opts%verbosity)
    write(*,*) '... constructed fake forceconstant corresponding to Td=',trim(flyt2char(opts%debye_temperature))
    ! dump the stupid forceconstant
    call fc%writetofile(uc,'outfile.fakeforceconstant')
    write(*,*) '... wrote it to "outfile.fakeforceconstant", check that the frequency range is reasonable'
    ! now we have a fake forceconstant and can keep going!
else
    ! read the forceconstant from file
    call fc%readfromfile(uc,'infile.forceconstant')
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

! Clean copy to work with
p=ss
avgtemp=0.0_flyt
do i=1,opts%nconf
    ! reset the structure
    p%r=ss%r
    p%v=0.0_flyt
    p%u=0.0_flyt
    ! initialize
    call fcss%initialise_cell(p,opts%temperature,opts%zpm,.false.,opts%threshold)
    ! dump to file
    select case(opts%output_format)
    case(1) ! vasp output
        fname='contcar_conf'//trim(int2char_padded(i,4))
        call p%writetofile(trim(fname),opts%output_format,write_velocities=.true.)
    case(2) ! abinit output
        fname='abinput_conf'//trim(int2char_padded(i,4))
        call p%writetofile(trim(fname),opts%output_format,write_velocities=.true.)
    end select
    ! just measure some stuff, for no good reason
    ek=p%kinetic_energy()/(p%na-1)
    ep=fcss%potential_energy(p%u)/(p%na-1)
    temp=(ek+ep)/(3*lo_kb_ev)
    avgtemp=avgtemp+temp
    write(*,"(1X,'conf:',I4,' ek:',F10.5,' ep:',F10.5,' temp:',F10.2,' mean temp',F10.2)") i,ek,ep,temp,avgtemp/i
enddo

contains

end program
