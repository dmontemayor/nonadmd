!-----------------------------------------------------------------------------
!> \brief
!! Loads Hamiltonian and performs spectrometric observations and runs ILDM and PLDM dynamics. 
!<----------------------------------------------------------------------------
!> \details
!! Used to compare ABspectra, PLDM and ILDM dynamics with literature.
!<----------------------------------------------------------------------------
!>\authors
!! Daniel Montemayor
!!
!!\date
!! Jan 2012
!<----------------------------------------------------------------------------
!>\todo
!! Change this file to function as tests driver.
!<-----------------------------------------------------------------------------
! !>\bug
! !<-----------------------------------------------------------------------------
!Change Log:
!
! === v1.0 Jan 2012 ===
!-----------------------------------------------------------------------------
subroutine tests
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

  implicit none
# include "config.h"
  
  integer(short)::unit,ierr
  logical::usedunit,exists
  
  type(quantum),target::qs
  type(classical),target::cs
  type(coupling),target::cp

  !type(ABspectrum)::ABspec
  !type(SpectralDensity)::SpecDen

  type(PLDM)::Huomap
  type(ILDM)::islandmap

  character(len=title)::qstype
  character(len=title)::cstype
  character(len=title)::cptype

  real(double)::Emin,Emax
  !------------------------------------------
  
  if(myid.EQ.master)then
     call openLog(level=3,file='qstest.log')

     call Prompt("**********************************************")
     call Prompt("*        B E G I N   Q S   T E S T S         *")
     call Prompt("**********************************************")
     
     call Prompt('Please enter qs type')
     read(*,*)qstype
     qstype=adjustl(qstype)

     call new(qs,type=trim(qstype))
     call save(qs,file=EXEDIR//'/test.'//trim(qstype))
     call new(qs,file=EXEDIR//'/test.'//trim(qstype))
     call update(qs)
     call resample(qs)
     call display(qs)
     call kill(qs)

     call closeLog
     !---------------------------------------------------------------
     call openLog(level=3,file='cstest.log')
     call Prompt("**********************************************")
     call Prompt("*        B E G I N   C S   T E S T S         *")
     call Prompt("**********************************************")

     call Prompt('Please enter cs type')
     read(*,*)cstype
     cstype=adjustl(cstype)

     call new(cs,type=trim(cstype))
     call save(cs,file=EXEDIR//'/test.'//trim(cstype))
     call new(cs,file=EXEDIR//'/test.'//trim(cstype))
     call update(cs)
     call resample(cs)
     call display(cs)
     call kill(cs)

     call closeLog
     !---------------------------------------------------------------
     call openLog(level=3,file='cptest.log')
     call Prompt("**********************************************")
     call Prompt("*        B E G I N   C P   T E S T S         *")
     call Prompt("**********************************************")


     !initiate quantum and classical subsystems for coupling subsystem tests
     call new(qs,file=EXEDIR//'/test.'//trim(qstype))
     call new(cs,file=EXEDIR//'/test.'//trim(cstype))

     call Prompt('Please enter cp type')
     read(*,*)cptype
     cptype=adjustl(cptype)

     call new(cp,qs,cs,type=trim(cptype))
     call save(cp,file=EXEDIR//'/test.'//trim(cptype))
     call new(cp,qs,cs,file=EXEDIR//'/test.'//trim(cptype))
     call update(cp)
     call resample(cp)
     call display(cp)
     call kill(cp)

     call closeLog
     !---------------------------------------------------------------
     call openLog(level=1,file='obsvtest.log')
     call Prompt("**********************************************")
     call Prompt("*        B E G I N   O B V   T E S T S       *")
     call Prompt("**********************************************")

     !initiate quantum, classical and coupling subsystems for observaton tests 
     call new(qs,file=EXEDIR//'/test.'//trim(qstype))
     call new(cs,file=EXEDIR//'/test.'//trim(cstype))
     call new(cp,qs,cs,file=EXEDIR//'/test.'//trim(cptype))


     call Prompt('Please enter Emin ABspectrum')
     read(*,*)Emin
     call Prompt('Please enter Emax ABspectrum')
     read(*,*)Emax
    ! call observe(ABspec,qs=qs,cs=cs,cp=cp&
    !      ,Emin=Emin,Emax=Emax&
    !      ,N=1000,samples=1000,file=EXEDIR//'/test.ABspec')


!!$     call Prompt('Please enter Emin Spectral Density')
!!$     read(*,*)Emin
!!$     call Prompt('Please enter Emax Spectral Density')
!!$     read(*,*)Emax
!!$     call observe(SpecDen,qs=qs,cs=cs,cp=cp&
!!$          ,Emin=Emin,Emax=Emax&
!!$          ,N=1000,samples=1000,file=EXEDIR//'/test.SpecDen')

     call closeLog
     !---------------------------------------------------------------
     call openLog(level=1,file='proptest.log')
     call Prompt("**********************************************")
     call Prompt("* B E G I N   P R O P A G A T O R  T E S T S *")
     call Prompt("**********************************************")

     !initiate quantum, classical and coupling subsystems for observaton tests 
     call new(qs,file=EXEDIR//'/test.'//trim(qstype))
     call new(cs,file=EXEDIR//'/test.'//trim(cstype))
     call new(cp,qs,cs,file=EXEDIR//'/test.'//trim(cptype))

     !call new(landmap,qs,cs,cp)
     !call save(landmap,file=EXEDIR//'/test.LDM')
     !call new(landmap,qs,cs,cp,file=EXEDIR//'/test.LDM')
     !call run(landmap,file=EXEDIR//'/test.LDM.den')
     !call display(landmap)
     !call kill(landmap)

     call new(Huomap,qs,cs,cp)
     call save(Huomap,file=EXEDIR//'/test.PLDM')
     call new(Huomap,qs,cs,cp,file=EXEDIR//'/test.PLDM')
     call run(Huomap,file=EXEDIR//'/test.PLDM.den')
     call display(Huomap)
     call kill(Huomap)

     call new(islandmap,qs,cs,cp)
     call save(islandmap,file=EXEDIR//'/test.ILDM')
     call new(islandmap,qs,cs,cp,file=EXEDIR//'/test.ILDM')
     call run(islandmap,file=EXEDIR//'/test.ILDM.den')
     call display(islandmap)
     call kill(islandmap)

     call closeLog
     !------------------------------------------------------------------------
     call openLog(level=1,file='systest.log')
     call Prompt("**********************************************")
     call Prompt("*     B E G I N   S Y S T E M   T E S T S    *")
     call Prompt("**********************************************")

     call closeLog     
     !========================================================================
     call Prompt("**********************************************")
     call Prompt("*     P A S S E D   A L L   T E S T S  !     *")
     call Prompt("**********************************************")
  endif
  
end subroutine tests
