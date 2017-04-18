!-----------------------------------------------------------------------------
!> \brief
!! Tests Hamiltonian. 
!<----------------------------------------------------------------------------
!> \details
!! Used to ensure that the Hamiltonian  quantum, classical, and coupling
!! subsystems are compatible and have proper NEW, DISPLAY, UPDATE, RESAMPLE,
!! SAVE, and KILL and CHECK functionality.
!<----------------------------------------------------------------------------
!>\authors
!! Daniel Montemayor
!!
!!\date
!! Jan 2012
!<----------------------------------------------------------------------------
!>\todo
!! Add CHECK functionality tests.
!<-----------------------------------------------------------------------------
! !>\bug
! !<-----------------------------------------------------------------------------
!Change Log:
!
! === v1.0 Jan 2012 ===
!-----------------------------------------------------------------------------
subroutine testHamiltonian
  use type_kinds
  !use MPIframework
  use wallclock
  use string
  use RunTimeLog
  use rand
  use atomicunits
  
  use quantum_class
  use classical_class
  use coupling_class

  !use spectrometer_class
  !use LDM_class  

  implicit none
# include "config.h"
  
  integer(short)::unit,ierr
  logical::usedunit,exists
  
  type(quantum),target::qs
  type(classical),target::cs
  type(coupling),target::cp

  !type(ABspectrum)::ABspec
  !type(SpectralDensity)::SpecDen

  !type(LDM)::landmap
  !type(PLDM)::Huomap
  !type(ILDM)::islandmap

  character(len=title)::qstype
  character(len=title)::cstype
  character(len=title)::cptype

  !real(double)::Emin,Emax
  !------------------------------------------
  
  if(myid.EQ.master)then
     call openLog(level=3,file='testHamiltonian.log')

     call Prompt("**********************************************")
     call Prompt("*        B E G I N   Q S   T E S T S         *")
     call Prompt("**********************************************")
     
     Call Prompt('Please enter qs type')
     read(*,*)qstype
     qstype=adjustl(qstype)

     call new(qs,type=trim(qstype))
     call save(qs,file=EXEDIR//'/testHamiltonian.'//trim(qstype))
     call new(qs,file=EXEDIR//'/testHamiltonian.'//trim(qstype))
     call update(qs)
     call resample(qs)
     call display(qs)
     call kill(qs)

     call Prompt("**********************************************")
     call Prompt("*        P A S S   Q S   T E S T S    !      *")
     call Prompt("**********************************************")
     !---------------------------------------------------------------
     call Prompt("**********************************************")
     call Prompt("*        B E G I N   C S   T E S T S         *")
     call Prompt("**********************************************")

     call Prompt('Please enter cs type')
     read(*,*)cstype
     cstype=adjustl(cstype)

     call new(cs,type=trim(cstype))
     call save(cs,file=EXEDIR//'/testHamiltonian.'//trim(cstype))
     call new(cs,file=EXEDIR//'/testHamiltonian.'//trim(cstype))
     call update(cs)
     call resample(cs)
     call display(cs)
     call kill(cs)

     call Prompt("**********************************************")
     call Prompt("*        P A S S   C S   T E S T S    !      *")
     call Prompt("**********************************************")
     !---------------------------------------------------------------
     call Prompt("**********************************************")
     call Prompt("*        B E G I N   C P   T E S T S         *")
     call Prompt("**********************************************")

     !initiate quantum and classical subsystems for coupling subsystem tests
     call new(qs,file=EXEDIR//'/testHamiltonian.'//trim(qstype))
     call new(cs,file=EXEDIR//'/testHamiltonian.'//trim(cstype))

     call Prompt('Please enter cp type')
     read(*,*)cptype
     cptype=adjustl(cptype)

     call new(cp,qs,cs,type=trim(cptype))
     call save(cp,file=EXEDIR//'/testHamiltonian.'//trim(cptype))
     call new(cp,qs,cs,file=EXEDIR//'/testHamiltonian.'//trim(cptype))
     call update(cp)
     call resample(cp)
     call display(cp)
     call kill(cp)

     call Prompt("**********************************************")
     call Prompt("*        P A S S   C P   T E S T S    !      *")
     call Prompt("**********************************************")
     !========================================================================
     call Prompt("**********************************************")
     call Prompt("*     P A S S E D   A L L   T E S T S  !     *")
     call Prompt("**********************************************")
     call closeLog
  endif
  
end subroutine testHamiltonian
