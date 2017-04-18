!-----------------------------------------------------------------------------
!> \brief
!! Use Praxis to parameterize spectral density equation from data.
!<----------------------------------------------------------------------------
!> \details
!! Provided with this driver is data set containing spectral density data
!! caluclated from a molecular dynamics simulation. In this example we will
!! calculate the parameters of an analytic form of the spectral density 
!! that best fit the data set provided. The analytic form of the  density can be
!! found in EQ2.14 in Garg et al, J. Chem. Phys. 83 (9), 1 November 1985.  
!<----------------------------------------------------------------------------
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

  ! The form to fit is F(x)=eta*x*Omega**4/((Omega**2-x**2)**2+4*x**2*gamma**2)
  ! The parameters to the function are eta, Omega, and gamma and should be
  !   contained in the param array.
  !
  !1)Edit the declaration of the param array above such that it has 3 elements instead of 2.

  !2) Praxis needs and initial guess for each parameter. Enter some reasonable numbers here.
  param(1)=1.0   ! this is eta     (kondo friction parameter)
  param(2)=100   ! this is Omega   (undamped frequency)
  param(3)=100   ! this is gamma   (related to friction)
  
  !3) Call the fitting function with these initial parameters and specifying the number of parameters
  call praxisfit(param,3) !<- what is wrong here?

  !print the parameters in units of wavenumbers to the screen after fitting
  write(*,*)'eta =',param(1)
  write(*,*)'Omega =',param(2)
  write(*,*)'gamma =',param(3)


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

  ! The form to fit is F(x)=eta*x*Omega**4/((Omega**2-x**2)**2+4*x**2*gamma**2)
  !param(1) is eta
  !param(2) is Omega
  !param(3) is gamma

  if(check(trim(file)).NE.0)call stop('cannot find file '//trim(file))
  unit=newunit()
  open(unit,file=trim(file))
  ierr=0
  DO while(ierr.eq.0)
     read(unit,*,iostat=ierr)xobs,yobs
     if(ierr.EQ.0)then

          !change this equation for ycal to the form we want
        ycal=param(1)*xobs*exp(-param(2)*xobs) !example here is f(x)=eta*x*exp(-x/Omega)


        F=F+(ycal-yobs)**2  !minimizes error calulated as sum of the diference squared


     end if
  END DO
  close(unit)  
  

END FUNCTION F
