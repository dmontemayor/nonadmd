!-----------------------------------------------------------------------------
!> \brief
!! Parallel PLDM dynamics of a 7 state FMO model
!! with site specific spectral densities.
!<----------------------------------------------------------------------------
!> \details
!! An FMO model is parameterized for a 7 state quantum subsystem. Each state
!! is bi-linearly coupled to an independent dissipative bath. Each bath
!! consists of 60 harmonic oscillators with spectral density defined in 
!! the file 'MK22.46.jw'.
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
  use PLDM_class  

  implicit none
# include "config.h"

  integer(long)::ierr

  type(quantum),target::qs
  type(classical),target::cs
  type(coupling),target::cp

  type(spectrum)::spec

  type(PLDM)::prop
  !------------------------------------------

# ifdef MPI
  call MPI_INIT(ierr)
  call MPI_COMM_rank(MPI_COMM_World,myid,ierr)
  call MPI_COMM_size(MPI_COMM_World,nproc,ierr)
# endif

  !Developers can run system tests here
  !call tests

  call openLog(level=1)!,file=EXEDIR//'/runtime.log'&
       !,ProgramName=PACKAGE_STRING,BugReport=PACKAGE_BUGREPORT)
  !--------------------   Body of experiment   ----------------------

  if(myid.EQ.master)then
     call Prompt('Creating Quantum Subsystem.')
     call NEW(qs,type='hs')                 !< Create new qs
     !call display(qs)                      !< Display newly created qs
     call save(qs,file='fmo.qs')            !< Save qs for later use
     call kill(qs)
  end if
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif
  !create quantum subsystem from file
  call NEW(qs,file='fmo.qs')!< Load qs from save file for each processor

  if(myid.EQ.master)then
     call Prompt('Creating Classial Subsystem.')
     call NEW(cs,type='harmonicbath')       !< Create new cs
     !call display(cs)                      !< Display newly created cs
     call save(cs,file='fmo.cs')            !< Save cs for later use
     call kill(cs)
  end if
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif
  call NEW(cs,file='fmo.cs')!< Load cs from save file for each processor

  if(myid.EQ.master)then
     call Prompt('Creating Coupling Term.')
     call NEW(cp,qs=qs,cs=cs,type='DiagonalCaldeiraLeggett')!< Create new cp
     !call display(cp)                           !< Display newly created cp
     call save(cp,file='fmo.cp')                 !< Save cp for later use
     call kill(cp)
  end if
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif
  call NEW(cp,qs=qs,cs=cs,file='fmo.cp')!< Load cp from save file for each processor
  call Prompt('Hamiltonian Complete!')

  if(myid.EQ.master)then
     call Prompt('Setting Up Propagator.')
     call NEW(prop,qs,cs,cp)                      !< Create new propagator
     !call display(prop)                          !< Display propagator
     call save(prop,file='fmo.prop')              !< Save propagator for later
     call kill(prop)
  end if
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif
  call NEW(prop,qs,cs,cp,file='fmo.prop')!< Load propagator from save file for each processor

  call Prompt('Running Propagator.')
  call RUN(prop,file=trim(int2str(myid))//'.den') !<each processor will run the propagation scheme

  !-----------------------------------------------------------------
  call closeLog

# ifdef MPI
  call MPI_Finalize(ierr)
# endif

end program main
