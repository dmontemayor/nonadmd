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
  use outputdisplay

  !hamiltonian  
  use quantum_class
  use classical_class
  use coupling_class
  
  !spectrometer
  use spectrometer_class

  !propagators
  use PLDM_class  

  implicit none
# include "config.h"

  type(quantum),target::qs
  type(classical),target::cs
  type(coupling),target::cp
  type(spectrum)::spec
  type(PLDM)::prop

  integer(long)::ierr

  !parameters used for fitting
  real(double)::param(2)
  real(double),external::f

!------------------------------------------

!!$# ifdef MPI
!!$  call MPI_INIT(ierr)
!!$  call MPI_COMM_rank(MPI_COMM_World,myid,ierr)
!!$  call MPI_COMM_size(MPI_COMM_World,nproc,ierr)
!!$# endif

  !Developers can run system tests here
  !call tests

  call openLog(level=1_short,file='./runtime.log')
!!$  call openLog(level=1,file=EXEDIR//'/runtime.log'&
!!$       ,ProgramName=PACKAGE_STRING,BugReport=PACKAGE_BUGREPORT)
  call display('display.out')
  !--------------------   Body of experiment   ----------------------

  if(myid.EQ.master)then
     call Prompt('Creating Quantum Subsystem.')
     !call NEW(qs,type='hs')                     !< Create new qs
     !call NEW(qs,file='example.qs')             !< Create qs from save file
     !call display(qs)                           !< Display newly created qs
     !call save(qs,file='example.qs')            !< Save qs for later use

     call Prompt('Creating Classial Subsystem.')
     !call NEW(cs,type='hb')                     !< Create new cs
     !call NEW(cs,file='example.cs')             !< Create cs from save file
     !call display(cs)                           !< Display newly created cs
     !call save(cs,file='example.cs')            !< Save cs for later use
     
     call Prompt('Creating Coupling Term.')
     !call NEW(cp,qs=qs,cs=cs,type='hc')         !< Create new cp
     !call NEW(cp,qs=qs,cs=cs,file='example.cp') !< Create cp from save file
     !call display(cp)                           !< Display newly created cp
     !call save(cp,file='example.cp')            !< Save cp for later use

     call Prompt('Hamiltonian Complete!')

     call Prompt('Making pre-production observatons.')
     !call Observe(spec,qs,cs,cp,file='example.pre.spec')

     call Prompt('Setting Up Propagator.')
     !call NEW(prop,qs,cs,cp)                     !< Create new propagator
     !call NEW(prop,qs,cs,cp,file='example.prop') !< Load propagator from file
     !call display(prop)                          !< Display propagator
     !call save(prop,file='example.prop')         !< Save propagator for later

  end if
!!$# ifdef MPI
!!$  call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$# endif

  call Prompt('Running Propagator.')
  !call RUN(prop,file='example.run') !< Conduct the propagation scheme
  if(myid.EQ.master)then
     call Prompt('Fitting data.')
     param=1.0                !initial guess of parameters: here all params=1.0
     !call praxisfit(param,2)
     write(*,*)param          !print parameters after fitting
  end if
!!$# ifdef MPI
!!$  call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$# endif

  if(myid.EQ.master)then
     call Prompt('Making post-production observatons.')
  end if
!!$# ifdef MPI
!!$  call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$# endif

  !-----------------------------------------------------------------
  call closeLog

!!$# ifdef MPI
!!$  call MPI_Finalize(ierr)
!!$# endif

end program main
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!function to fit
FUNCTION F(param,n)
  use type_kinds
  use filemanager
  use errorlog
  
  INTEGER(long),intent(in)::n
  REAL(double),intent(inout):: param(n)
  INTEGER(long):: i,ierr,unit
  REAL(double):: F,xobs,yobs,ycal
  character(len=path)::file='./example.dat'
  !character(len=path)::file=EXAMPLESDIR//'example.dat'
  F=0.0

!!$  if(check(trim(file)).NE.0)call stop('cannot find file '//trim(file))
!!$  unit=newunit()
!!$  open(unit,file=trim(file))
!!$  ierr=0
!!$  DO while(ierr.eq.0)
!!$     read(unit,*,iostat=ierr)xobs,yobs
!!$     if(ierr.EQ.0)then
!!$        ycal=param(1)*xobs*exp(-param(2)*xobs) !example fit to ohmic with debey cuttoff
!!$        F=F+(ycal-yobs)**2         !minimize error calulated as sum of the diference squared
!!$     end if
!!$  END DO
!!$  close(unit)  
  

END FUNCTION F
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
