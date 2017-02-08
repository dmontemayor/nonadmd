!-----------------------------------------------------------------------------
!> \brief
!! PLDM dynamics of Spin-Boson model.
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
!! Jan 2012
!<----------------------------------------------------------------------------
program main
  use type_kinds
  use MPIframework
  use wallclock
  use string
  use ErrorLog
  use rand
  use atomicunits

  use quantum_class
  use classical_class
  use coupling_class
  
  use spectrometer_class
  use LDM_class
  use islandmapserial_class  

  implicit none
# include "config.h"

  integer(short)::ierr
  real(double)::tol,dummy

  type(quantum),target::qs
  type(classical),target::cs
  type(coupling),target::cp

  type(spectrum)::spec

 ! type(verlet)::prop
  type(PLDM)::prop
 ! type(ILDM)::prop
 !type(islandmapserial)::prop
!------------------------------------------

# ifdef MPI
  call MPI_INIT(ierr)
  call MPI_COMM_rank(MPI_COMM_World,myid,ierr)
  call MPI_COMM_size(MPI_COMM_World,nproc,ierr)
# endif

  !Developers can run system tests here
  !call tests

  if(myid.EQ.master)then
     write(*,*)"**********************************************"
     write(*,*)"*         B E G I N   P R O G R A M          *"
     write(*,*)"* "//PACKAGE_STRING
     write(*,*)"*                                            *"
     write(*,*)"* Report Bugs to:"
     write(*,*)"* "//PACKAGE_BUGREPORT
     write(*,*)"* Execution Time:                            *"
     write(*,"(A3,8(I4,1X),A4)")" * ",timearray(),       "   *"
     write(*,*)"**********************************************"
  end if
  dummy=ran0()
  write(*,*)'core:'//trim(int2str(myid))//', seed:'//trim(int2str(seed()))
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif

  call openLog(level=1,file=EXEDIR//'/runtime.log')
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
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif
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
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif
  !create classical subsystem from file
  call NEW(cs,file='harmonicbath.cs')

  !setup coupling scheme
  if(myid.EQ.master)then
     call Prompt('Creating Coupling Term')
     !create coupling subsystem (cp) and link with qs and cs
     call NEW(cp,qs=qs,cs=cs,type='DiagonalCaldeiraLeggett')
     call display(cp)
     !save coupling subsystem to file
     call save(cp,file='DiagonalCaldeiraLeggett.cp')
     call Kill(cp)
  end if
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif
  !create coupling subsystem from file and link with qs and cs
  call NEW(cp,qs=qs,cs=cs,file='DiagonalCaldeiraLeggett.cp')
     
  !make pre-production observatons
  if(myid.EQ.master)then
     !call spectrometer observation 
  end if
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif
  
  !setup experiment
  if(myid.EQ.master)then
     call Prompt('Creating Propagator')
     !create propagator and link with qs, cs, and cp
     call NEW(prop,qs,cs,cp)

     !call display(prop)

     !save propagator to file
     call save(prop,file='PLDM.prop')
  end if
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif
  !create propagator from file and link with qs, cs, and cp
  call NEW(prop,qs,cs,cp,file='PLDM.prop')

  !run propagator
  call Run(prop,'den.out')!,tol)
  !call Run(prop,file='den.out',normalize=.true.)

  !make post-production observatons
  if(myid.EQ.master)then
     !call spectrometer observation
  end if
# ifdef MPI
  call MPI_barrier(MPI_COMM_WORLD,ierr)
# endif

  call closeLog
  !-----------------------------------------------------------------
  if(myid.EQ.master)then
     write(*,*)"**********************************************"
     write(*,*)"* S U C C E S S F U L   E N D   P R O G R A M*"
     write(*,*)"*                                            *"
     write(*,*)"* End Time:                                  *"
     write(*,"(A3,8(I4,1X),A4)")" * ",timearray(),       "   *"
     write(*,*)"**********************************************"
  end if


# ifdef MPI
  call MPI_Finalize(ierr)
# endif

end program main
