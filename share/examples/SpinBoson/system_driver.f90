!-----------------------------------------------------------------------------
!> \brief
!! The Spin-Boson model.
!<----------------------------------------------------------------------------
!> \details
!! Examples from: Bonella and Coker, JCP 122, 194102 (2005)
!! figures 2,3,and 4. Examples are generated from input file input files 
!! provided with this driver. 
!<----------------------------------------------------------------------------
!>\authors
!! Daniel Montemayor
!!
!!\date
!! Jan 2012, rev. Jul 2014
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
  !type(SH)::prop

  integer(long)::ierr

  !parameters used for fitting
  real(double)::param(2),t,pop
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
  !setup quantum subsystem
  if(myid.EQ.master)then
     call Prompt('Creating Quantum Subsystem')
     !call quantum subsystem (qs) creator
     call NEW(qs,type='hs')
     call display(qs)
     !save quantum subsystem to file
     call save(qs,file='hs.qs')
     call kill(qs)
  end if
!!$# ifdef MPI
!!$  call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$# endif
  !create quantum subsystem from file
  call NEW(qs,file='hs.qs')

  !setup classical subsystem
  if(myid.EQ.master)then
     call Prompt('Creating Classical Subsystem')
     !call classical subsystem (cs) creator
     call NEW(cs,type='harmonicbath')
     call display(cs)
     !save classical subsystem to file
     call save(cs,file='harmonicbath.cs')
     call kill(cs)
  end if
!!$# ifdef MPI
!!$  call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$# endif
  !create classical subsystem from file
  call NEW(cs,file='harmonicbath.cs')

  !setup coupling scheme
  if(myid.EQ.master)then
     call Prompt('Creating Coupling Term')
     !create coupling subsystem (cp) and link with qs and cs
     call NEW(cp,qs=qs,cs=cs,type='bilinear')
     call display(cp)
     !save coupling subsystem to file
     call save(cp,file='bilinear.cp')
     call Kill(cp)
  end if
!!$# ifdef MPI
!!$  call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$# endif
  !create coupling subsystem from file and link with qs and cs
  call NEW(cp,qs=qs,cs=cs,file='bilinear.cp')
     
  !make pre-production observatons
  if(myid.EQ.master)then
     !call spectrometer observation 
  end if
!!$# ifdef MPI
!!$  call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$# endif
  
  !setup experiment
  if(myid.EQ.master)then
     call Prompt('Creating Propagator')
     !create propagator and link with qs, cs, and cp
     call NEW(prop,qs,cs,cp)

     !call display(prop)

     !save propagator to file
     call save(prop,file='PLDM.prop')
  end if
!!$# ifdef MPI
!!$  call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$# endif
  !create propagator from file and link with qs, cs, and cp
  call NEW(prop,qs,cs,cp,file='PLDM.prop')

  !run propagator
  call Run(prop,file='SpinBoson.den')

  !make post-production observatons
  if(myid.EQ.master)then
     open(200,file='praxis.in')
     open(100,file='SpinBoson.den',iostat=ierr)
     do while (ierr.EQ.0)
        read(100,*,iostat=ierr)t,pop
        !write population difference
        if(ierr.EQ.0)write(200,*)abs(t*qs%hs%diabat(1,2)),2*pop-1
     end do
     close(100)
     close(200)

     !praxis fit
     param(1)=.5_double !Rate initial guess
     param(2)=0._double !Eq value initial guess
     call praxisfit(param,2)
     write(*,*)'Population difference decay rate= ',param(1)
     write(*,*)'Population difference eq value= ',param(2)
  end if
!!$# ifdef MPI
!!$  call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$# endif

  call closeLog
  !-----------------------------------------------------------------

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
  character(len=path)::file='./praxis.in'
  F=0.0

  if(check(trim(file)).NE.0)call stop('cannot find file '//trim(file))
  unit=newunit()
  open(unit,file=trim(file))
  ierr=0
  DO while(ierr.eq.0)
     read(unit,*,iostat=ierr)xobs,yobs
     if(ierr.EQ.0)then
        ycal=(1._double-param(2))*exp(-param(1)*xobs)+param(2) !example fit to exponetial decay param(1)=rate, param(2)=eq value
        F=F+(ycal-yobs)**2         !minimize error calulated as sum of the diference squared
     end if
  END DO
  close(unit)  
  

END FUNCTION F
