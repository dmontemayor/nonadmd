!-----------------------------------------------------------------------------
!> \brief
!! Main loop of Program!
!<----------------------------------------------------------------------------
!> \details
!! Here the user can build a Hamiltonian of the form \f$H=hs+hb+hc\f$
!! and perform spectroscopic measurements. The user can also choose an
!! appropriate propagation method to compute time resolved observables. 
!<----------------------------------------------------------------------------
!>\authors
!! Daniel Montemayor
!!
!!\date
!! Jan 2012
!<----------------------------------------------------------------------------
program main
  use type_kinds
  use atomicunits
  use math
  use rand
  use string
  use MPIframework
  use filemanager
  use ErrorLog
  use molreader
  use wallclock

  !hamiltonian  
  use quantum_class
  use classical_class
  use coupling_class
  
  !spectrometer
  use spectrometer_class

  !propagators
  use verlet_class
  use islandmapserial_class
  use LDM_class  

  implicit none
# include "config.h"

  integer(long)::ierr

  type(quantum),target::qs
  type(classical),target::cs
  type(coupling),target::cp

  type(spectrum)::spec

  type(verlet)::prop
  !type(islandmapserial)::prop
  !type(PLDM)::prop
  !type(ILDM)::prop
!------------------------------------------

# ifdef MPI
  call MPI_INIT(ierr)
  call MPI_COMM_rank(MPI_COMM_World,myid,ierr)
  call MPI_COMM_size(MPI_COMM_World,nproc,ierr)
# endif

  !Developers can run system tests here
  !call tests

  call openLog(level=1,file=EXEDIR//'/runtime.log'&
       ,ProgramName=PACKAGE_STRING,BugReport=PACKAGE_BUGREPORT)
  !--------------------   Body of experiment   ----------------------
  if(myid.EQ.master)then
     call Prompt('Creating Quantum Subsystem.')
     call NEW(qs,type='MorseOscillator')        !< Create new qs
     call display(qs)                           !< Display newly created qs
!stop
     call save(qs,file='example.qs')            !< Save qs for later use

     call Prompt('Creating Classial Subsystem.')
     call NEW(cs,type='harmonicbath')           !< Create new cs
     call display(cs)                           !< Display newly created cs
     call save(cs,file='example.cs')            !< Save cs for later use
     
     call Prompt('Creating Coupling Term.')
     call NEW(cp,qs=qs,cs=cs,type='DMBLcoupling')!< Create new cp
     call display(cp)                            !< Display newly created cp
     call save(cp,file='example.cp')             !< Save cp for later use

     call Prompt('Hamiltonian Complete!')

     call Prompt('Making pre-production observatons.')
     !call Observe(spec,qs,cs,cp,file='example.pre.spec')

     call Prompt('Setting Up Propagator.')
     call NEW(prop,qs,cs,cp)              !< Create new propagator
     !call display(prop)                  !< Display propagator
     call save(prop,file='example.prop')  !< Save propagator for later

  end if
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif



  call Prompt('Running Propagator.')
  call RUN(prop,file='example.run') !< Conduct the propagation scheme


  if(myid.EQ.master)then
     call Prompt('Making post-production observatons.')
     !call Observe(spec,qs,cs,cp,file='example.post.spec')
  end if
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif

  !-----------------------------------------------------------------
  call closeLog

# ifdef MPI
  call MPI_Finalize(ierr)
# endif

end program main
