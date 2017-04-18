!============================================================
!> \brief
!! Definition of integer, floating point, and string kinds.
!! \details
!! Cross-Platform portability is the goal such that type
!! precsision is maintained from build to build of the software.
!! \authors
!! Author: Daniel Montemayor
!<===========================================================
module type_kinds
  implicit none

  Private

  !Integer Kinds
  Public :: byte,short,long

  !Floating point Kinds
  Public :: single,double

  !String Lengths
  Public :: title,path,line,comment 

  !> string length for a name/title
  integer, parameter:: title=50 
  !> string length for a directory path
  integer, parameter:: path=100
  !> string length for a line in a file
  integer, parameter:: line=100
  !> string length for a message/comment
  integer, parameter:: comment=500

  !> byte integer size 
  integer, parameter:: byte=selected_int_kind(1)
  !> short integer size 
  integer, parameter:: short=selected_int_kind(4)
  !> long integer size 
  integer, parameter:: long=selected_int_kind(8)

  !> single precision float size 
  integer, parameter:: single=selected_real_kind(6)
  !> double precision float size 
  integer, parameter:: double=selected_real_kind(15)

end module type_kinds

!===========================================================
!> \brief
!! Definition of mathematical constants and operators.
!! \authors
!! Author: Daniel Montemayor
!! \todo
!! * Add cross and dot product operators
!<===========================================================
module math
  use type_kinds
  implicit none

  real(double),parameter::pi=3.1415926535898_double
  real(double),parameter::twopi=pi*2.0_double
  real(double),parameter::fourpi=pi*4.0_double
  real(double),parameter::halfpi=pi/2.0_double
  !real(double),parameter::sqrt2=sqrt(2.0_double)
  !real(double),parameter::invsqrt2=1.0_double/sqrt2
  complex(double),parameter::eye=(0.0_double,1.0_double)
  public::iden,trace,cnorm,fnorm

contains
  !> \brief identity matrix
  !! \param[in] N size of matrix
  !!<
  function iden(N)
    integer(long),intent(in)::N
    real(double)::iden(N,N)
    integer(long)::i

    iden=0._double
    do i=1,N
       iden(i,i)=1._double
    end do
  end function iden

  !> \brief trace of a complex matrix
  !! \param[in] A matrix
  !!<
  function trace(A)
    complex(double),intent(in)::A(:,:)
    real(double)::trace
    integer(long)::i,N
    N=size(A,1)

    trace=0._double
    do i=1,N
       trace=trace+sqrt(real(A(i,i)*conjg(A(i,i))))
    end do
  end function trace

  !> \brief Entrywise p-norm of a complex square matrix with p=1
  !! \param[in] A matrix
  !!<
  function Cnorm(A)
    complex(double),intent(in)::A(:,:)
    real(double)::Cnorm
    integer(long)::i,j,N
    N=size(A,1)

    cnorm=0._double
    do i=1,N
       do j=1,N
          cnorm=cnorm+sqrt(real(A(i,j)*conjg(A(i,j))))
       end do
    end do
  end function Cnorm

  !> \brief Entrywise p-norm of a complex matrix
  !! \param[in] A matrix
  !! \param[in] p optional p-norm integer default: p=2 Frobinius normalization 
  !!<
  function Fnorm(A,p)
    complex(double),intent(in)::A(:,:)
    integer(long),intent(in),optional::p
    real(double)::Fnorm,s,sh,si
    integer(long)::i,j,N,M
    N=size(A,1)
    M=size(A,2)

    fnorm=0._double
    s=2._double
    if(present(p))then
       s=real(p,double)
       sh=s/2._double
       si=1._double/s
       do i=1,N
          do j=1,M
             fnorm=fnorm+real(A(i,j)*conjg(A(i,j)))**sh
          end do
       end do
       fnorm=fnorm**si
    else
       do i=1,N
          do j=1,M
             fnorm=fnorm+real(A(i,j)*conjg(A(i,j)))
          end do
       end do
       fnorm=sqrt(fnorm)
    end if
  end function Fnorm

end module math
!===========================================================
!> \brief
!! Definition of atomic units.
!! \details
!! Atomic units are the default used throughout the software.
!! \author
!! Daniel Montemayor
!<===========================================================
module atomicunits
  use type_kinds
  use math
  implicit none

  !Fundamental Units
  real(double),parameter::a0=1.0_double   !<(Fundamental length) Bohr radius
  real(double),parameter::me=1.0_double   !<(Fundamental mass) electron mass
  real(double),parameter::hbar=1.0_double !<(Fundamental action) reduced Planks const
  real(double),parameter::e=1.0_double    !<(Fundamental charge) electron charge
  real(double),parameter::Kel=1.0_double  !<(Fundamental temperature) Kelvin
  real(double),parameter::mol=1.0_double  !<(Fundamental amount) Avogadro's number
  real(double),parameter::rad=1.0_double  !<(Fundamental angle) radian

  !Derived Units
  real(double),parameter::Eh=(hbar/a0)**2/me                      !<Hartree (unit value) 
  real(double),parameter::kC=Eh*a0/(e**2)                         !<Coulomb force (unit value)
  real(double),parameter::mp=1836.152663_double*me                !<proton rest mass
  real(double),parameter::mn=1838.685239_double*me                !<neutron rest mass
  real(double),parameter::kb=3.166815203_double*1E-6*Eh/Kel       !<boltzman constant
  real(double),parameter::eV=3.674932587_double*1E-2*Eh           !<Electon volt
  real(double),parameter::c0=137.0359996287515_double*hbar/(a0*me)!<speed of light
  real(double),parameter::eps0=1.0_double/(fourpi*kC)             !<vacuum permittivity
  real(double),parameter::Debye=.3934302014076827_double*e*a0     !<Dipole moment
  real(double),parameter::Deg=twopi/360_double*rad                !<Degrees

  !Converstion Factors
  real(double),parameter::angstrom=1.889726133921252_double*a0    !<angstrom
  real(double),parameter::fs=41.34137337414122_double*hbar/Eh     !<femtosecond
  real(double),parameter::ps=fs*1_double*1E3                      !<picosecond
  real(double),parameter::invcm=4.5663352672_double*1E-6*Eh/hbar  !<wavenumber(1/cm)

end module atomicunits
!===========================================================
!===========================================================
!> \brief
!! Timming Module
!! \details
!! Supplies the software with a standard method of computing
!! wall clock time in units of mils, sec, hours or even days.
!! \authors
!! Author: Daniel Montemayor
!<===========================================================
module wallclock
  use type_kinds
  implicit none
  
  private
  public:: date,time,zone,timearray
  
     character(len=8)::d
     character(len=10)::t
     character(len=5)::z
     integer(long) :: a(8)
  
contains
  !----------------------------
  function date()!<returns date in 8 char string format YYYYMMDD
    character(len=8)::date
    call date_and_time(d,t,z,a)
    date=d
  end function date
  !----------------------------
  function time()!<returns time in 10 char string format hhmmss.sss
    character(len=10)::time
    call date_and_time(d,t,z,a)
    time=t
  end function time
  !----------------------------
  function zone()!<returns zone relative to Coordinated Univeral Time (UTC) in 5 char string format (+-)hhmm
    character(len=5)::zone
    call date_and_time(d,t,z,a)
    zone=z
  end function zone
  !----------------------------
  function timearray()!<returns 8 element integer array representing time in format Yr,Mo,Day,UTCmin_dif,hr,min,sec,msec
    integer(long) :: timearray(8)
    call date_and_time(d,t,z,a)
    timearray=a
  end function timearray

end module wallclock
!===========================================================
!> \brief
!! Enables MPI
!! \details
!! Checks if MPI is available from preprocessor and sets up appropriate variables
!<
module MPIframework
  use type_kinds
  implicit none
  
# ifdef MPI
  include "mpif.h"
  integer(long)::status(MPI_STATUS_SIZE)
# endif
  
  integer(long), parameter :: master=0
  integer(long)::myid=0
  integer(long)::nproc=1
  
end module MPIFRAMEWORK
!===========================================================
module filemanager
!> \brief
!! IO utility
!! \details 
!! Controls Fortran unit assignments for read and write statements and checks
!! if files exist.
!< --------------------------------------------------------------
  use type_kinds
  use MPIFRAMEWORK
  implicit none
  
  private
  public:: newunit,check

  interface check
     module procedure checkfile
  end interface

contains
  !----------------------------
  integer(long) function newunit()
    logical::usedunit  
    newunit=1000+myid*100
    usedunit=.true.
    do while (usedunit)
       newunit=newunit+1
       inquire(newunit,opened=usedunit)
    end do
    return
  end function newunit
  !-----------------------------
  integer(short) function checkfile(file)
    character*(*),intent(in)::file
    logical::there
    
    checkfile=0
    inquire(file=file,exist=there)
    if(.not.there)then
       checkfile=1
    end if
    return
  end function checkfile

end module filemanager
!===========================================================
module rand
!> \brief
!! Random number generators
!! \todo
!! * correlated random numbers
!<

  use type_kinds
  use wallclock
  use MPIframework
  implicit none

  private
  public :: seed,setseed,ran0,gran

  integer(long) :: idum=0_long

contains
  subroutine setseed(seedin)
    IMPLICIT NONE
    integer(long) :: seedin
    idum=seedin
    RETURN
  END subroutine setseed
  !----------------------------
  function seed()
    IMPLICIT NONE
    integer(long) :: seed
    seed=idum
    RETURN
  END function seed
  !----------------------------
  FUNCTION RAN0(seedin)
    IMPLICIT NONE
    REAL(double):: RAN0
    integer(long), optional:: seedin
    INTEGER(long),PARAMETER:: IA=16807_long
    INTEGER(long),PARAMETER:: IM=2147483647_long
    INTEGER(long),PARAMETER:: IQ=127773_long
    INTEGER(long),PARAMETER:: IR=2836_long
    INTEGER(long),PARAMETER:: MASK=20410876_long
    REAL(double),PARAMETER:: AM=1.0_double/IM
    INTEGER(long):: k

    integer(long)::array(8)

    if(present(seedin))idum=seedin
    if(idum.EQ.0)then
       do idum=0, myid*100000!E5
          array=timearray()
       end do
       idum=myid*1E4
       idum=idum+array(8)
       idum=idum+array(7)*1E3
       idum=idum+array(6)*6E4
       idum=idum+array(5)*36E5
       !write(*,*)"core ",myid,"generating new random seed ",idum)
    end if

    idum=IEOR(idum,MASK)
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    IF (idum.LT.0) idum=idum+IM
    RAN0=AM*idum
    idum=IEOR(idum,MASK)
    RETURN
  END FUNCTION RAN0
  !----------------------------
  FUNCTION GRAN(seed)
    !         FUNCTION TO GENERATE GAUSSIAN RANDOM NUMBERS
    !         USING BOX-MUELLER/WEINER METHOD.
    IMPLICIT NONE
    integer(long), optional:: seed
    REAL(double) :: GRAN,w1,w2
    real(double) , parameter :: PI=3.1415926535898_double

    if(present(seed))idum=seed

    w1=RAN0()
    w2=RAN0()
    gran=SQRT(-2.*LOG(w1))*COS(2.*PI*w2)

  END FUNCTION GRAN

end module rand
!===========================================================---------
!> \brief
!! Get strings from numbers and other string operations
module string
  use type_kinds
  use MPIFRAMEWORK
  use filemanager
  private
  
  Public::quote
  Public::int2str,float2str,int2hex
  Public::str2int,str2float,hex2int

  interface int2str
     module procedure long2str
     module procedure short2str
  end interface
  interface float2str
     module procedure double2str
     module procedure single2str
  end interface

  interface str2int
     module procedure str2long
  end interface
  interface str2float
     module procedure str2double
  end interface
  interface hex2int
     module procedure hex2long
  end interface
  interface int2hex
     module procedure long2hex
  end interface

contains
  !----------------------------
  function quote(str)
    character*(*),intent(in)::str
    character(len=len(str)+2)::quote

    quote="'"//str//"'"
  end function quote
  !----------------------------
  function long2str(I)

    integer(long),intent(in)::I
    character(len=title)::long2str
    integer(long)::ierr

    if((I.NE.I).or.(abs(I).GT.huge(I)))then
       long2str='int2str_ERROR'
    else
       write(long2str,'(I10)',iostat=ierr)I
       if(ierr.EQ.1)long2str='int2str_ERROR'
    end if

    long2str=adjustl(long2str)

  end function long2str
  !----------------------------
  function short2str(I)

    integer(short),intent(in)::I
    integer(long)::longI
    character(len=title)::short2str

    longI=I
    short2str=trim(int2str(longI))
    
  end function short2str
  !----------------------------
  function str2long(string)

    character*(*),intent(in)::string
    integer(long)::str2long,order,i
    integer(long)::ierr,tryhex

    !make sure all characters are appropriate
    order=len(string)
    tryhex=0
    do i=1,order
       select case(string(i:i))
       case(' ')
       case('0')
       case('1')
       case('2')
       case('3')
       case('4')
       case('5')
       case('6')
       case('7')
       case('8')
       case('9')
       case default
          !write(*,*)'str2int error: '//string&
          !     //' contains non integer characters!'
          tryhex=1
       end select
    end do
    if(tryhex.EQ.1)then
       !attempt hexadecimal
       str2long=hex2long(string)
       !write(*,*)'in hex: '//string//' =',str2long
    else
       read(string,*,iostat=ierr)str2long
       if(ierr.EQ.1)then
          write(*,*)'str2int error: '//string//' is not an integer! Program will stop!'
          stop
       end if
    end if
  end function str2long
  !----------------------------
  function double2str(X)

    real(double),intent(in)::X
    character(len=title)::double2str
    integer(long)::ierr

    if((X.NE.X).or.(abs(X).GT.huge(X)))then
       double2str='float2str_ERROR'
    else
       write(double2str,'(E14.6E3)',iostat=ierr)X
       if(ierr.EQ.1)double2str='float2str_ERROR'
    end if

    double2str=adjustl(double2str)
  end function double2str
  !----------------------------
  function str2double(string)

    character*(*),intent(in)::string
    real(double)::str2double
    integer(long)::ierr

    read(string,*,iostat=ierr)str2double
    if(ierr.EQ.1)then
       write(*,*)'str2float error: '//string//' is not a float! Program will stop!'
       stop
    end if

  end function str2double
  !----------------------------
  function single2str(X)
    
    real(single),intent(in)::X
    real(double)::doubleX
    character(len=title)::single2str
    
    doubleX=X
    single2str=trim(double2str(doubleX))
    
  end function single2str
  !----------------------------
  function hex2long(string)

    character*(*),intent(in)::string
    integer(long)::hex2long,order,i,j

    !get length of hex
    order=len(string)

    !make sure all characters are appropriate
    do i=1,order
       select case(string(i:i))
       case('0')
       case('1')
       case('2')
       case('3')
       case('4')
       case('5')
       case('6')
       case('7')
       case('8')
       case('9')
       case('a')
       case('b')
       case('c')
       case('d')
       case('e')
       case('f')
       case('A')
       case('B')
       case('C')
       case('D')
       case('E')
       case('F')
       case default
          write(*,*)'hex2int error: '//string//' is not a hexadecimal number! Program will stop.'
          stop
       end select
    end do

    hex2long=0
    do i=1, order
       select case(string(i:i))
       case('0')
          j=0
       case('1')
          j=1
       case('2')
          j=2
       case('3')
          j=3
       case('4')
          j=4
       case('5')
          j=5
       case('6')
          j=6
       case('7')
          j=7
       case('8')
          j=8
       case('9')
          j=9
       case('a')
          j=10
       case('b')
          j=11
       case('c')
          j=12
       case('d')
          j=13
       case('e')
          j=14
       case('f')
          j=15
       case('A')
          j=10
       case('B')
          j=11
       case('C')
          j=12
       case('D')
          j=13
       case('E')
          j=14
       case('F')
          j=15
       end select
       hex2long=hex2long+j*16**(order-i)
    end do
  end function hex2long
  !----------------------------
  function long2hex(I)

    integer(long),intent(in)::I
    character(len=title)::long2hex
    character::char
    integer(long)::order,j,remainder,digit

    !get length of hex
    order=0
    do while(I.GE.16**order)
       order=order+1
    end do

    remainder=I
    long2hex=''
    do j=1, order
       digit=0
       do while(remainder.GE.(digit+1)*16**(order-j))
          digit=digit+1
       end do
       remainder=remainder-digit*16**(order-j)
       if(digit.LT.10)then
          char=trim(int2str(digit))
       else
          select case(digit)
          case(10)
             char='a'
          case(11)
             char='b'
          case(12)
             char='c'
          case(13)
             char='d'
          case(14)
             char='e'
          case(15)
             char='f'
          case default
             write(*,*)'error hex2int: cannot convert integer to hexadecimal! Program will stop!'
             stop
          end select
       end if 
          long2hex=trim(long2hex)//char
    end do
  end function long2hex
end module string

!===========================================================
!>\brief
!! Error Handleing
!! \details
!! Enables Error, Warning, Note, and Prompt messages by setting
!! a log level. Also allows user to save ErrorLog messages to a file.
!<------------------------------------------------------------------- 
module Errorlog
  !-------------------    ChangeLog    ---------------------
  ! - NOTES and PROMPTS only returned by master processor
  !Author: Daniel Montemayor
  !---------------------------------------------------------
  use type_kinds
  use wallclock
  use MPIframework
  use filemanager
  use string
  implicit none
  private

  public::note,warn,stop,openLog,closeLog,prompt,displayLog

  integer(short)::LOGLevel=0
  integer(long)::LOGunit=6
  character(len=path)::LOGfile='stout'
  logical::logclosed=.True.
  logical::logopen=.False.
  integer(long)::LogCount=0

  interface note
     module procedure log_note
  end interface
  
  interface warn
     module procedure log_warn
  end interface

  interface stop
     module procedure log_stop
  end interface

  interface prompt
     module procedure log_prompt
  end interface

contains
  !----------------------------
  subroutine displaylog
    
    write(6,*)'Log Level= '//trim(int2str(LOGlevel))
    write(6,*)'Log File= '//trim(LOGfile)
    if(logopen)write(6,*)'Log is currently open.'
    if(logclosed)write(6,*)'Log is currently closed.'
    write(6,*)'Log Count= '//trim(int2str(LOGCount))
    
  end subroutine displaylog
  !----------------------------
  subroutine openlog(file,level,ProgramName,BugReport)
    character*(*),intent(in),optional::file
    integer(short),intent(in),optional::level
    character*(*),intent(in),optional::ProgramName
    character*(*),intent(in),optional::BugReport
    logical::fileopened
    if(myid.EQ.master)then    
       
       if(logclosed)then

          LOGcount=0

          LOGLevel=0
          if(present(level))LOGlevel=level
       
          if(present(file))then
             LOGfile=file
             LOGunit=newunit()
             inquire(file=trim(LOGfile),opened=fileopened)
             if(.not.fileopened)open(LOGunit,file=trim(LOGfile))
          else
             LOGfile='the runtime messages above!'
             LOGunit=6
             open(LOGunit)
          end if

              
          write(LOGunit,*)"**********************************************"
          write(LOGunit,*)"*        B E G I N   E R R O R   L O G       *"
          write(LOGunit,*)"*                                            *"
          write(LOGunit,*)"* Log Level = "//trim(int2str(LogLevel))&
               //"                              *"
          write(LOGunit,*)"* Create Time:                               *"
          write(LOGunit,"(A3,8(I4,1X),A4)")" * ",timearray(),       "   *"
          write(LOGunit,*)"*                                            *"
          if(present(ProgramName))&
               write(LOGunit,*)"* Program Title: "//ProgramName
          if(present(BugReport))&
               write(LOGunit,*)"* Report Bugs to: "//BugReport
          write(LOGunit,*)"**********************************************"
          
          logclosed=.false.
       end if
       logopen=.not.logclosed

    end if

  end subroutine openlog
  !----------------------------
  subroutine closelog
    if(myid.EQ.master)then
       if(logopen)then
          
          write(LOGunit,*)"**********************************************"
          write(LOGunit,*)"*          E N D   E R R O R   L O G         *"
          write(LOGunit,*)"*                                            *"
          write(LOGunit,*)"* Close Time:                                *"
          write(LOGunit,"(A3,8(I4,1X),A4)")" * ",timearray(),       "   *"
          write(LOGunit,*)"*                                            *"
          write(LOGunit,*)"**********************************************"
          close(LOGunit)
          
          logclosed=.true.
       end if
       logopen=.not.logclosed
    end if
  end subroutine closelog
  !----------------------------
  subroutine log_note(msg)
    character*(*),intent(in)::msg
    logical::record,halt
    
    call openlog
    
    if(myid.EQ.master)then
       
       record=.false.
       halt=.false.
       
       select case(LOGlevel)
       case(0)
       case(1)
       case(2)
          record=.true.
       case(3)
          record=.true.
       case(4)
          record=.true.
          halt=.true.
       end select
       
       if(record)then
          Logcount=LogCount+1
          write(LOGunit,*)trim(int2str(LogCount))//' NOTE: '//msg
       end if
       
       if(halt)call stop('NOTES will STOP program: '//msg)

    end if

  end subroutine log_note
  !----------------------------
  subroutine log_warn(msg,recover)
    character*(*),intent(in)::msg
    character*(*),intent(in),optional::recover
    
    logical::record,halt
    
    call openlog
    
    record=.false.
    halt=.false.
    
    select case(LOGlevel)
    case(0)
    case(1)
       record=.true.
    case(2)
       record=.true.
    case(3)
       record=.true.
       halt=.true.
    case(4)
       record=.true.
       halt=.true.
    end select
    
    if(myid.NE.master)record=.false.

    if(record)then
       Logcount=LogCount+1
       write(LOGunit,*)trim(int2str(LogCount))// ' WARNING: '//msg       
       if(present(recover))write(LOGunit,*)trim(int2str(LogCount))&
            // ' WARNING: Program will recover by '//recover
    end if
    
    if(halt)call stop('WARNINGS will STOP program: '//msg)
    
  end subroutine log_warn
  !----------------------------
  subroutine log_stop(msg,suggest)
    character*(*),intent(in)::msg
    character*(*),intent(in),optional::suggest
    call openlog

    call Prompt('What did you do?! Check '//trim(LOGfile) )
    
    Logcount=LogCount+1
    
    if(myid.EQ.master)then
       write(LOGunit,*)trim(int2str(LogCount))//' ERROR: '//msg
       if(present(suggest))write(LOGunit,*)trim(int2str(LogCount))&
            //' ERROR: Try '//suggest
       write(LOGunit,*)'Program will stop!'
    else
       write(*,*)trim(int2str(LogCount))//' ERROR: '//msg
       if(present(suggest))write(*,*)trim(int2str(LogCount))&
            //' ERROR: Try '//suggest
       write(*,*)'Program will stop!'
    end if
    call closeLog
    stop

  end subroutine log_stop
  !----------------------------
  subroutine log_prompt(msg)
    character*(*),intent(in)::msg
    call openlog
    if(myid.EQ.master)then
       Logcount=LogCount+1
       write(6,*)trim(int2str(LogCount))//' PROMPT: '//msg
       if(LOGunit.NE.6)write(LOGunit,*)trim(int2str(LogCount))//' PROMPT: '//msg
    end if
  end subroutine log_prompt
  
end module Errorlog
!===========================================================
!===========================================================
!>\brief
!! Enables Display to file 
!! \details
!! Offers the user the ability to control where DISPLAY command will output to.
!<------------------------------------------------------------------- 
module OutputDisplay
  use type_kinds
  use wallclock
  use MPIframework
  use filemanager
  use string
  implicit none
  private

  public::Display,Dunit

  integer(long)::Dunit=6

  interface Display
     module procedure Display_init
  end interface

contains
  !-------------------------------------
  !>\brief redirects display output to newfile
  !-------------------------------------
  subroutine display_init(file)
    character*(*),intent(in)::file

    if(myid.EQ.master)then
       Dunit=newunit()
       open(Dunit,file=file)
    end if
  end subroutine display_init

end module OutputDisplay


!-----------------------------------------------------
!> \brief
!! Procedures to Fit XY data to functional forms
!!\details
!!Returns select function parameters
!<------------------------------------------------------
module fitXYdata
  implicit none
  private
  public::fit2exp
  
  real*8,allocatable::XYdata(:,:)
  
  interface fit2exp
     module procedure fit2exp
  end interface

contains
  !------------------------------------------------------------------
  !> \brief
  !! Fit to exponential function
  !!\details
  !!fits set to exp function (A-B)*exp(-R*t)+B
  !<------------------------------------------------------
  subroutine fit2exp(set,R,A,B)

    use rand
    implicit none

    real*8, intent(in)::set(:,:)
    real*8, intent(inout)::R,A
    real*8,intent(inout),optional::B


    real*8,parameter::tol=1E-10
    real*8 :: dtn,dt,dt2,domain
    real*8 :: V,V0,Vp,Vm
    real*8,dimension(3):: X,Xdot,dV,dX,X0

    integer :: istep,ipt,npt,idof,ndof,minpt(1),maxpt(1)

    logical :: newiteration

    !write(*,*)'begin fit2exp'
    !write(*,*)'size(set)',size(set(:,1))
    !write(*,*)'domain',minval(set(:,1)),maxval(set(:,1))
    !write(*,*)'R',R
    !write(*,*)'A',A
    !if(present(B))then
    !   write(*,*)'B',B
    !else
    !   write(*,*)'B is not present'
    !end if

    ndof=2
    if(present(B))ndof=3

    !setup
    npt=size(set(:,1))
    domain=maxval(set(:,1))-minval(set(:,1))
    dtn=domain/real(npt)
    dt=dtn
    dt2=dt*dt

    !initial guess
    maxpt=maxloc(set(:,2))
    minpt=minloc(set(maxpt(1):,2))

    X(3)=0.0
    if(present(B)) X(3)=set(minpt(1),1)
    X(1)=-log((set(minpt(1),2)+X(3))/(set(maxpt(1),2)+X(3)))&
         /(set(minpt(1),1)-set(maxpt(1),1))
    X(2)=X(3)+(set(maxpt(1),2)+X(3))/exp(-X(1)*set(maxpt(1),1))

    dX=epsilon(X)

    !calibrate initial variables
    V0=sum((set(:,2)-((X(2)-X(3))*exp(-X(1)*set(:,1))+X(3)))**2)
    V=V0

    !write(*,*)
    !write(*,*)'ndof=',ndof
    !write(*,*)'V0=',V0
    !write(*,*)'X=',X(:ndof)
    !write(*,*)'dX=',dX(:ndof)
    !write(*,*)'max point index=',maxpt(1)
    !write(*,*)'min point index=',minpt(1)
    !write(*,*)'maxval=',set(maxpt(1),2)
    !write(*,*)'minval=',set(minpt(1),2)
    !stop

    X0=X
    do idof=1,ndof
       newiteration=.true.
       do while(newiteration)
          X(idof)=X0(idof)+dX(idof)
          Vp=sum((set(:,2)-((X(2)-X(3))*exp(-X(1)*set(:,1))+X(3)))**2)

          X(idof)=X0(idof)-dX(idof)
          Vm=sum((set(:,2)-((X(2)-X(3))*exp(-X(1)*set(:,1))+X(3)))**2)

          if(abs(Vp-Vm).LE.epsilon(V))then
             dX(idof)=dX(idof)*1.05
             !write(*,*)'gradient too small',idof,Vp-Vm,dX(idof)
          else
             newiteration=.false.
             !if(Vp-Vm.LT.0.)dX(idof)=-dX(idof)
             !write(*,*)'gradient is good',idof,Vp-Vm,dX(idof)
          end if
       end do

       !stop

       dV(idof)=0.
       if(abs(dX(idof)).GT.epsilon(dX)) dV(idof)=.5*(Vp-Vm)/dX(idof)
       X(idof)=X0(idof)

       !Calibrate integration time step
       do while(abs(dV(idof))*dt2.GT.abs(X(idof))*1E-3)
          dt=dt/2.0
          dt2=dt*dt
          write(*,*)'dt',dt,dV(idof),X(idof)*1E-3
          if(dt.LE.epsilon(dt))stop
       end do
    end do

    !write(*,*)
    !write(*,*)'initial guess'
    !write(*,*)'X',X(:ndof)
    !write(*,*)'dX',dX(:ndof)
    !write(*,*)'dV',dV(:ndof)
    !stop


    !start from rest
    Xdot=0.

    istep=0
    do while ((abs(V-V0)/V0.GT.tol).or.(istep.EQ.0))
       newiteration=.true.
       do while(V0.GT.V.or.newiteration)
          istep=istep+1
          V0=V

          X=X+Xdot*dt-.5*dV*dt2
          Xdot=Xdot-.5*dV*dt

          !compute  dV
          V=sum((set(:,2)-((X(2)-X(3))*exp(-X(1)*set(:,1))+X(3)))**2)
          do idof=1,ndof
             X0(idof)=X(idof)
             X(idof)=X0(idof)+dX(idof)
             Vp=sum((set(:,2)-((X(2)-X(3))*exp(-X(1)*set(:,1))+X(3)))**2)

             X(idof)=X0(idof)-dX(idof)
             Vm=sum((set(:,2)-((X(2)-X(3))*exp(-X(1)*set(:,1))+X(3)))**2)

             dV(idof)=0.
             if(dX(idof).NE.0.0) dV(idof)=.5*(Vp-Vm)/dX(idof)
             X(idof)=X0(idof)
          end do

          Xdot=Xdot-.5*dV*dt
          newiteration=.false.
          !write(*,*)'---------------------'
          !write(*,*)'X=',X(:ndof)
          !write(*,*)'Xdot=',Xdot(:ndof)
          !write(*,*)'dV=',dV(:ndof)

       end do
       Xdot=-Xdot
       dt=dt*.95
       dt2=dt*dt
       !write(*,*)
       !write(*,*)'new iteration'
       !write(*,*)'istep',istep
       !write(*,*)'V=',V,V0,abs(V-V0)/V0
       !write(*,*)'X=',X(:ndof)
       !write(*,*)'dt=',dt
       !stop
    end do

    !write(*,*)'final result'
    !write(*,*)X

    R=X(1)
    A=X(2)
    if(present(B))  B=X(3)

  end subroutine fit2exp
!----------------------------------------------------------------------

end module fitXYdata
!---------------------------------------------------------------
!> \brief
!! Text based plotting.
!! \details
!! A set of procedures accessed with the DISPLAY call that return
!! text based graphs of real or complex matrix data and 2D vector
!! fields.
!<--------------------------------------------------------------- 
module textgraphs
  use type_kinds
  use string
  use math
  use errorlog
  use outputdisplay
  implicit none
  
  private
  public:: display
  
  interface display
     module procedure displayCMPLXMatrix
     module procedure displayREALMatrix
     module procedure displayVectField
     module procedure displayValues
     module procedure displayXY
  end interface
contains
  !----------------------------
  !> \brief
  !! Plot Complex Matrix
  !! \details
  !! Displays a 2D complex matrix elements as 2 char
  !! text in NxM table, where N and M are the
  !! dimensions of the input matrix.
  !<--------------------------------------------------
  subroutine displayCMPLXMatrix(Z,tol,msg)
    complex(double),intent(in),target::Z(:,:)
    real(double),intent(in),optional::tol
    character*(*),intent(in),optional::msg
    integer(long)::N,M,i,j
    character(len=title)::FMT

    real(double)::mag

    write(Dunit,*)'____________________________________________________'
    write(Dunit,*)'---------------   CMPLX Matrix  --------------------'
    if(present(msg))write(Dunit,*)msg
    write(Dunit,*)'----------------------------------------------------'

    N=size(Z,1)
    M=size(Z,2)

    mag=sqrt(sum(Z**2))/real(N*M)
    if(present(tol))mag=tol
    mag=abs(mag) !make sure magnitute is positive

    FMT='('//trim(int2str(M))//'(A2,1X))'

    write(Dunit,FMT)((textgraphsChar(real(Z(i,j)),tol=mag,alt=.false.)&
         //textgraphsChar(aimag(Z(i,j)),tol=mag,alt=.false.),j=1,M),i=1,N)
    
    write(Dunit,*)'----------------------------------------------------'
    write(Dunit,*)
  end subroutine displayCMPLXMatrix
  !----------------------------
  !>\brief
  !! Plot Real Matrix
  !! \details
  !! Displays a 2D real matrix elements as 1 char 
  !! text in NxM table, where N and M are the
  !! dimensions of the input matrix.
  !<---------------------------------------------
  subroutine displayREALMatrix(A,tol,msg)
    real(double),intent(in),target::A(:,:)
    real(double),intent(in),optional::tol
    character*(*),intent(in),optional::msg
    integer(long)::N,M,i,j
    character(len=title)::FMT
    real(double)::mag

    write(Dunit,*)'____________________________________________________'
    write(Dunit,*)'---------------   REAL Matrix   --------------------'
    if(present(msg))write(Dunit,*)msg
    write(Dunit,*)'----------------------------------------------------'

    N=size(A,1)
    M=size(A,2)
    mag=sqrt(sum(A**2))/real(N*M)
    if(present(tol))mag=tol
    mag=abs(mag) !make sure magnitute is positive

    FMT='('//trim(int2str(M))//'(A1,2X))'

    write(Dunit,FMT)((textgraphsChar(A(i,j),tol=mag),j=1,M),i=1,N)

    write(Dunit,*)'----------------------------------------------------'
    write(Dunit,*)
  end subroutine displayREALMatrix
  !-------------
  !>\brief
  !! Scaled Value character
  !! \details
  !! Returns a character based double precision input.
  !!\n
  !!-----character key------\n
  !! %+ : significant real positive value\n
  !! %- : significant real negative value\n
  !! %p : significant imag positive value\n
  !! %n : significant imag negative value\n
  !! %. : insignificant value\n
  !!------------------------\n
  !<----------------------------------------
  character function textgraphsChar(X,tol,alt)
    real(double),intent(in)::X
    real(double),intent(in),optional::tol
    logical,intent(in),optional::alt

    logical::alttxt
    real(double)::mag

    alttxt=.false.
    if(present(alt))alttxt=alt

    mag=.10_double
    if(present(tol))mag=tol
    mag=abs(mag) !make sure magnitute is positive

    textgraphsChar='.'
    if(X.GT.mag)then
       textgraphsChar='+'
       if(alttxt)textgraphsChar='p'
    end if
    if(X.LT.-mag)then
       textgraphsChar='-'
       if(alttxt)textgraphsChar='n'
    end if
    
  end function textgraphsChar
  !-------------
  !> \brief
  !! Vector character
  !! \details
  !! Returns a character representing the orientation of
  !! an input 2D vector.
  !! \n
  !!-----character key------\n
  !! %. : insignificant magnitute\n
  !! %- : pi15/8 < angle < pi/8   achar(45)\n
  !! %/ : pi/8   < angle < pi3/8  achar(47)\n
  !! %| : pi3/8  < angle < pi5/8  achar(124)\n
  !! %\ : pi5/8  < angle < pi7/8  achar(92)\n
  !! %- : pi7/8  < angle < pi9/8\n
  !! %/ : pi9/8  < angle < pi11/8\n
  !! %| : pi11/8 < angle < pi13/8\n
  !! %\ : pi13/8 < angle < pi15/8\n
  !!------------------------\n
  !<------------------------------------------
  character function textgraphsVect(X,tol)
    real(double),intent(in)::X(2)
    real(double),intent(in),optional::tol

    real(double)::mag,ang,mtol,ethpi

    ethpi=pi/8._double

    mtol=.10_double
    if(present(tol))mtol=tol
    mtol=abs(mtol) !make sure magnitute tol is positive

    mag=sqrt(sum(X**2))
    ang=0
    ang=atan2(X(2),X(1))
    if(ang.LT.0)ang=ang+pi !make sure angle is positive

    textgraphsVect='.'
    if(mag.GT.mtol)then
       if(ang.GT.15*ethpi.or.ang.LE.ethpi)textgraphsVect=achar(45)
       if(ang.GT.ethpi.and.ang.LE.3*ethpi)textgraphsVect=achar(47)
       if(ang.GT.3*ethpi.and.ang.LE.5*ethpi)textgraphsVect=achar(124)
       if(ang.GT.5*ethpi.and.ang.LE.7*ethpi)textgraphsVect=achar(92)
       if(ang.GT.7*ethpi.and.ang.LE.9*ethpi)textgraphsVect=achar(45)
       if(ang.GT.9*ethpi.and.ang.LE.11*ethpi)textgraphsVect=achar(47)
       if(ang.GT.11*ethpi.and.ang.LE.13*ethpi)textgraphsVect=achar(124)
       if(ang.GT.13*ethpi.and.ang.LE.15*ethpi)textgraphsVect=achar(92)
       
    end if

  end function textgraphsVect
  !----------------------------
  !> \brief
  !! Plot Vector Field
  !! \details
  !! Displays a 2D vector field elements as 
  !!nptx x npty char text.
  !<------------------------------
  subroutine displayVectField(Q,Mu,tol,nptx,npty,labelvect,msg)
    real(double),intent(in),target::Q(:,:),Mu(:,:)
    real(double),intent(in),optional::tol
    integer(long),intent(in),optional::nptx,npty
    logical,intent(in),optional::labelvect
    character*(*),intent(in),optional::msg

    character,allocatable::record(:)
    character::vectchar
    character(len=title)::vectlabel
    character(len=title)::FMT
    integer(long)::N,M,i,j,k,nx,ny,ni,nj,nk,pad
    real(double)::mtol,xmin,xmax,ymin,ymax,dx,dy,x,y
    logical::pass,label

    pad=5

    nx=line/2
    if(present(nptx))nx=nptx
    ny=line/4
    if(present(npty))ny=npty

    N=size(Q,1)
    M=size(Mu,2)

    pass=.true.
    if(size(Q,1).NE.size(Mu,1))pass=.false.
    if(size(Q,2).NE.size(Mu,2))pass=.false.
    if(N.LT.1)pass=.false.
    if(M.LT.2)pass=.false.
    if(nx.LT.1)pass=.false.
    if(ny.LT.1)pass=.false.
    if(.not.pass)call Warn('DisplayVectField: cannot plot vector 2d vector field!','displaying nothing.')


    label=.true.
    !if(present(labelvect))label=labelvect

    if(allocated(record))deallocate(record)
    allocate(record(0:nx+pad))
    FMT='(A'//trim(int2str(nx+pad+1))//')'

    if(pass)then
       
       write(Dunit,*)'____________________________________________________'
       write(Dunit,*)'--------------  2D  Vector Field  ------------------'
       if(present(msg))write(Dunit,*)msg
       write(Dunit,*)'----------------------------------------------------'

       !compute lattice spacing
       xmin=minval(Q(:,1))
       xmax=maxval(Q(:,1))
       dx=(xmax-xmin)/real(nx)
       ymin=minval(Q(:,2))
       ymax=maxval(Q(:,2))
       dy=(ymax-ymin)/real(ny)

       !compute significant magnitute
       mtol=0._double
       do i=1,N
          mtol=mtol+sqrt(sum(Mu(i,1:2)**2))
       end do
       mtol=mtol/real(10*N)
       if(present(tol))mtol=tol
       mtol=abs(mtol) !make sure magnitute is positive
       
       do j=0,ny
          y=ymin+j*dy

          record=' '
          do i=0,nx
             x=xmin+i*dx

             vectchar=' '
             vectlabel=' '
             do k=1,N
                if((Q(k,1).GT.x-.5*dx)&
                     .and.(Q(k,1).LE.x+.5*dx)&
                     .and.(Q(k,2).GT.y-.5*dy)&
                     .and.(Q(k,2).LE.y+.5*dy))then
                   if(vectchar.EQ.' ')then
                      vectchar=textgraphsVect(Mu(k,1:2),tol=mtol)
                      nj=len(trim(int2str(k)))
                      vectlabel(1:nj)=''//trim(int2str(k))
                      record(i:i)=vectchar
                      if(label)then
                         do ni=1,nj
                            record(i+ni:i+ni)=vectlabel(ni:ni)
                         end do
                      end if
                   end if
                end if
             end do
          end do
          write(Dunit,*)record
       end do
       write(Dunit,*)'----------------------------------------------------'
    end if
    
    write(Dunit,*)
    if(allocated(record))deallocate(record)
  end subroutine displayVectField
  !----------------------------
  !>\brief
  !!Display Array Vaules
  !!\details
  !!Plots a function contained in Q with text
  !!each element of Q gets a character on a the record
  !!with npt represents the range of values
  !< ---------------------
  subroutine displayValues(Q,npt,msg)
    real(double),intent(in),target::Q(:)
    integer(long),intent(in),optional::npt
    character*(*),intent(in),optional::msg

    character,allocatable::record(:)
    character(len=title)::FMT

    integer(long)::N,ny,i,j
    real(double)::ymin,ymax,dy,y
    logical::pass

    ny=line/4
    if(present(npt))ny=npt

    N=size(Q)
    FMT='('//trim(int2str(N))//'(A1))'

    pass=.true.
    if(ny.LT.1)pass=.false.
    if(.not.pass)call warn('DisplayValues: cannot plot function!','diplaying nothing.')

    if(allocated(record))deallocate(record)
    allocate(record(N))
    
    if(pass)then
       
       write(Dunit,*)'____________________________________________________'
       write(Dunit,*)'--------------    Array Vaules    ------------------'
       if(present(msg))write(Dunit,*)msg
       write(Dunit,*)'----------------------------------------------------'

       !compute lattice spacing
       ymin=minval(Q)
       ymax=maxval(Q)
       dy=(ymax-ymin)/real(ny)
       write(Dunit,*)'ymax= '//trim(float2str(ymax))       
       do j=ny,0,-1
          y=ymin+j*dy
          record=' '
          if(j.EQ.ny.or.j.EQ.0)record='-'
          do i=1,N
             !record(i:i)=' '
             if(i.EQ.1.or.i.EQ.N)record(i:i)='|'
             !if(Q(i).GT.y-.5*dy)record(i:i)=':'
             !if(Q(i).LE.y+.5*dy)record(i:i)='*'
             if((Q(i).GT.y-.5*dy).and.(Q(i).LE.y+.5*dy))record(i:i)='*'
          end do
          write(Dunit,FMT)record
       end do
       write(Dunit,*)'ymin= '//trim(float2str(ymin))       
       write(Dunit,*)'----------------------------------------------------'
    end if

    write(Dunit,*)
    if(allocated(record))deallocate(record)
  end subroutine displayValues
  !----------------------------
  !> \brief
  !! Scatter plot
  !! details/
  !! displays xy data as nptx x npty char text.
  !< -------------------------------------------------
  subroutine displayXY(Xset,Yset,nptx,npty,mask,msg)
    real(double),intent(in),target::Xset(:),Yset(:)
    integer(long),intent(in),optional::nptx,npty
    logical,intent(in),optional::mask(:)
    character*(*),intent(in),optional::msg

    character(len=title)::FMT
    character,allocatable::record(:)
    integer(long)::N,M,i,j,k,nx,ny
    real(double)::xmin,xmax,ymin,ymax,dx,dy,x,y
    logical::pass

    nx=line/2
    if(present(nptx))nx=nptx
    ny=line/4
    if(present(npty))ny=npty
    FMT='('//trim(int2str(nx+1))//'(A1))'

    N=size(Xset)
    M=N
    if(present(mask))M=size(mask)

    pass=.true.
    if(size(Xset).NE.size(Yset))pass=.false.
    if(M.NE.N)pass=.false.
    if(nx.LT.1)pass=.false.
    if(ny.LT.1)pass=.false.
    if(.not.pass)call warn('DisplayXY: cannot plot XY scatter plot!','displaying nothing.')

    if(allocated(record))deallocate(record)
    allocate(record(0:nx))
    
    if(pass)then
       
       write(Dunit,*)'____________________________________________________'
       write(Dunit,*)'--------------      XY plot       ------------------'
       if(present(msg))write(Dunit,*)msg
       write(Dunit,*)'----------------------------------------------------'

       !compute lattice spacing
       xmin=minval(Xset)
       xmax=maxval(Xset)
       dx=(xmax-xmin)/real(nx)
       ymin=minval(Yset)
       ymax=maxval(Yset)
       dy=(ymax-ymin)/real(ny)

       write(Dunit,*)'ymax='//trim(float2str(ymax))
       do j=ny,0,-1
          y=ymin+j*dy
          record=' '
          if(j.EQ.ny.or.j.EQ.0)record='-'
          record(0)='|'
          record(nx)='|'
          do i=0,nx
             x=xmin+i*dx

             do k=1,N
                pass=.true.
                if(present(mask))pass=mask(k)
                if(pass&
                     .and.(Xset(k).GT.x-.5*dx)&
                     .and.(Xset(k).LE.x+.5*dx)&
                     .and.(Yset(k).GT.y-.5*dy)&
                     .and.(Yset(k).LE.y+.5*dy))then
                   record(i:i)='*'
                end if
             end do
          end do
          write(Dunit,FMT)record
       end do
       write(Dunit,*)'ymin='//trim(float2str(ymin))
       write(Dunit,*)'xmin='//trim(float2str(xmin))//', xmax='//trim(float2str(xmax))
       write(Dunit,*)'----------------------------------------------------'
    end if

    write(Dunit,*)
    if(allocated(record))deallocate(record)
  end subroutine displayXY
  !----------------------------  
end module textgraphs
