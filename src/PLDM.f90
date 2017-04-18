!============================================================================
!>\brief
!! Partial Lineraized Path Integral Propagators
!!\details
!!\authors
!! Daniel Montemayor
!!\date
!! June 2012
!!\todo
!! * display output to file option
!<
!============================================================================
!--------------------------------changelog-----------------------------------
! + return of the intermediate file
! + added norming option to pldm_run
! + added popoly option to pldm_run
! + istate now properly sampled
! - remove harmonic oscillator hard code
! - remove diagonalcaldeiraleggett hard code
!============================================================================
module PLDM_class
  use type_kinds
  use math
  use atomicunits
  use wallclock
  use ErrorLog
  use filemanager
  use MPIframework
  use rand
  use string
  use quantum_class
  use classical_class
  use coupling_class
  implicit none
  Private

  Public::run,PLDM
  Public::new,kill,save,display,check

  real(double),allocatable::x0(:),p0(:),xt0(:),pt0(:)

  type PLDM
     logical::initialized=.false.
     type(quantum),pointer::qs
     type(classical),pointer::cs
     type(coupling),pointer::cp

     integer(long)::ntraj
     real(double)::runtime
     real(double)::dtn
     real(double)::dte
     real(double)::dtout
     real(double)::dtsave
     integer(long)::nstep,nbstep,nlstep
     !real(double),dimension(:),pointer::q,p,qt,pt
     !complex(double),dimension(:,:,:,:),pointer::buf
     !logical::bufok=.true.
  end type PLDM
  interface new
     module procedure PLDM_new
  end interface

  interface kill
     module procedure PLDM_kill
  end interface

  interface save
     module procedure PLDM_save
  end interface

  interface display
     module procedure PLDM_display
  end interface

  interface check
     module procedure PLDM_check
  end interface

  interface run
     module procedure PLDM_run
  end interface


contains
  !------------------------------------------------------------------------------
  subroutine mapverlet(this,nstate,nlit,dt,hel)!,x0,p0,xt0,pt0)
    type(PLDM),intent(inout)::this
    integer(long),intent(in)::nstate,nlit
    real(double),intent(in)::dt
    real(double),intent(in)::hel(nstate,nstate)

    !real(double),intent(inout)::x0(nstate),p0(nstate),xt0(nstate),pt0(nstate)
    !real(double)::x2(nstate),p2(nstate),xt2(nstate),pt2(nstate)

    
    real(double)::xrhs(nstate),prhs(nstate),xtrhs(nstate),ptrhs(nstate),force(nstate),forcet(nstate)
    real(double)::xsumdum,psumdum
    
    integer i,n,m
    
    do i=1,nlit
       
       !     Generate current derivatives
       
       !     Forward
       
       do n=1,nstate
          xrhs(n)=hel(n,n)*p0(n)
          prhs(n)=-hel(n,n)*x0(n)
       end do
       
       do n=1,nstate
          xsumdum=0.
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                xsumdum=xsumdum+hel(n,m)*p0(m)
                psumdum=psumdum+hel(n,m)*x0(m)
             end if
          end do
          xrhs(n)=xrhs(n)+xsumdum
          prhs(n)=prhs(n)-psumdum
       end do
       
       !     Backward
       
       do n=1,nstate
          xtrhs(n)=hel(n,n)*pt0(n)
          ptrhs(n)=-hel(n,n)*xt0(n)
       end do
       
       do n=1,nstate
          xsumdum=0.
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                xsumdum=xsumdum+hel(n,m)*pt0(m)
                psumdum=psumdum+hel(n,m)*xt0(m)
             end if
          end do
          xtrhs(n)=xtrhs(n)+xsumdum
          ptrhs(n)=ptrhs(n)-psumdum
       end do
       
       !     Generate current second derivatives
       
       !     Forward
       
       do n=1,nstate
          force(n)=hel(n,n)*prhs(n)
       end do
       
       
       do n=1,nstate
          xsumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                xsumdum=xsumdum+hel(n,m)*prhs(m)
             end if
          end do
          force(n)=force(n)+xsumdum
       end do
       
       !     backward
       
       do n=1,nstate
          forcet(n)=hel(n,n)*ptrhs(n)
       end do
       
       
       do n=1,nstate
          xsumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                xsumdum=xsumdum+hel(n,m)*ptrhs(m)
             end if
          end do
          forcet(n)=forcet(n)+xsumdum
       end do
       
       !     Advance Ver step
       
       do n=1,nstate
          x0(n)=x0(n)+xrhs(n)*dt+0.5*force(n)*dt*dt
          p0(n)=p0(n)+0.5*prhs(n)*dt
          
          xt0(n)=xt0(n)+xtrhs(n)*dt+0.5*forcet(n)*dt*dt
          pt0(n)=pt0(n)+0.5*ptrhs(n)*dt
       end do
       
       !     Compute new first derivatives
       
       !     Forward
       
       do n=1,nstate
          prhs(n)=-hel(n,n)*x0(n)
       end do
       
       do n=1,nstate
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                psumdum=psumdum+hel(n,m)*x0(m)
             end if
          end do
          prhs(n)=prhs(n)-psumdum
       end do
       
       !     Backward
       
       do n=1,nstate
          ptrhs(n)=-hel(n,n)*xt0(n)
       end do
       
       do n=1,nstate
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                psumdum=psumdum+hel(n,m)*xt0(m)
             end if
          end do
          ptrhs(n)=ptrhs(n)-psumdum
       end do
       
       !    Advance let step
       
       do n=1,nstate
          p0(n)=p0(n)+0.5*prhs(n)*dt
          pt0(n)=pt0(n)+0.5*ptrhs(n)*dt
       end do
       
    end do
    
    return
  end subroutine mapverlet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine PLDM_new(this,qs,cs,cp,file)
    implicit none
    type(PLDM),intent(inout)::this
    type(quantum),intent(in),target::qs
    type(classical),intent(in),target::cs
    type(coupling),intent(in),target::cp
    character*(*),intent(in),optional::file
    character(len=title)::filetype
    character(len=path)::filename
    integer(long)::unit

    integer(long)::ubound,lbound
    real(double)::value

    !check inputs
    call Note('Begin PLDM_new.')
    if(present(file))call Note('Input file= '//file)
    if(check(qs).EQ.1)call stop('PLDM_new, failed qs object check.')
    if(check(cs).EQ.1)call stop('PLDM_new, failed cs object check.')
    if(check(cp).EQ.1)call stop('PLDM_new, failed cp object check.')

    if(associated(this%qs))nullify(this%qs)
    this%qs=>qs
    if(associated(this%cs))nullify(this%cs)
    this%cs=>cs
    if(associated(this%cp))nullify(this%cp)
    this%cp=>cp    

    if(check(this%qs).EQ.1)call stop('PLDM_new, failed this%qs object check.')
    if(check(this%cs).EQ.1)call stop('PLDM_new, failed this%cs object check.')
    if(check(this%cp).EQ.1)call stop('PLDM_new, failed this%cp object check.')

    !load input file
    if(present(file))then
       unit=newunit()
       open(unit,file=file)
       read(unit,*)filetype
       filetype=adjustl(filetype)
       if(trim(filetype).NE.'PLDM')call stop('PLDM_init error: not a valid input file.')
       read(unit,*)this%ntraj
       read(unit,*)this%runtime
       read(unit,*)this%dtn
       read(unit,*)this%dte
       read(unit,*)this%dtout
       read(unit,*)this%dtsave
       close(unit)
    else
       this%ntraj=0
       do while(this%ntraj.LE.0)
          write(*,*)'Enter number of trajectories.'
          read(*,*)this%ntraj
          if(this%ntraj.LE.0)call warn('Number of trajectories must be positive.'&
               ,'requesting user to try again.')
       end do
       this%runtime=0._double
       do while(this%runtime.LE.0._double)
          write(*,*)'Enter total runtime.'
          read(*,*)this%runtime
          if(this%runtime.LE.0._double)call warn('total runtime must be positive.'&
               ,'requesting user to try again.')
       end do
       this%dtn=0._double
       do while(this%dtn.LE.0._double.OR.this%dtn.GT.this%runtime)
          write(*,*)'Enter nuclear time step.'
          read(*,*)this%dtn
          if(this%dtn.LE.0._double)call warn('nuclear time step must positive.'&
               ,'requesting user to try again.')
          if(this%dtn.GT.this%runtime)then
             call warn('nuclear time step cannot be larger than total runtime.'&
                  ,'requesting user to try again.')
          end if
       end do
       this%dte=0._double
       do while(this%dte.LE.0._double.or.this%dte.GT.this%dtn)
          write(*,*)'Enter electronic time step.'
          read(*,*)this%dte
          if(this%dte.LE.0._double)call warn('Electronic time step must positive.'&
               ,'requesting user to try again')
          if(this%dte.GT.this%dtn)call warn('Electronic time step cannot be larger than nuclear time step.'&
               ,'requesting user to try again.')
       end do
       this%dtout=0._double
       do while(this%dtout.LE.0._double.or.this%dtout.LT.this%dtn)
          write(*,*)'Enter output time step.'
          read(*,*)this%dtout
          if(this%dtout.LE.0._double)call warn('Output time step must be positive.'&
               ,'requesting user to try again.')
          if(this%dtout.LE.this%dtn)call warn('Output time step cannot be smaller than nuclear time step.'&
               ,'requesting user to try again.')
       end do
       this%dtsave=0._double
       do while(this%dtsave.LE.0._double)
          write(*,*)'Enter time interval in wallclock seconds between Restart files.'
          read(*,*)this%dtsave
          if(this%dtsave.LE.0._double)call warn('Save time step must positive.'&
               ,'requesting user to try again!')
       end do
    end if

    !snap output time step to grid
    value=this%runtime/this%dtout
    ubound=ceiling(value)
    lbound=floor(value)
    if(value.GE.(ubound+lbound)/2._double)then
       this%nstep=ubound
    else
       this%nstep=lbound
    end if
    this%dtout=this%runtime/real(this%nstep)

    !snap nuclear time step to grid
    value=this%runtime/this%dtn
    ubound=ceiling(value)
    lbound=floor(value)
    if(value.GE.(ubound+lbound)/2._double)then
       this%nbstep=ubound
    else
       this%nbstep=lbound
    end if
    this%dtn=this%runtime/real(this%nbstep)
    
    !snap electronic electronic time step to grid
    value=this%dtn/this%dte
    ubound=ceiling(value)
    lbound=floor(value)
    if(value.GE.(ubound+lbound)/2._double)then
       this%nlstep=ubound
    else
       this%nlstep=lbound
    end if
    this%dte=this%dtn/real(this%nlstep)

!!$    !initialize mapping variables
!!$    if(associated(this%q))nullify(this%q)
!!$    if(associated(this%p))nullify(this%p)
!!$    if(associated(this%qt))nullify(this%qt)
!!$    if(associated(this%pt))nullify(this%pt)
!!$
!!$    allocate(this%q(this%qs%hs%nstate),this%p(this%qs%hs%nstate))
!!$    allocate(this%qt(this%qs%hs%nstate),this%pt(this%qs%hs%nstate))

    this%initialized=.true.
    if(check(this).EQ.1)call stop('PLDM_new, failed final object check.')
  end subroutine PLDM_new
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine PLDM_kill(this)
    implicit none
    type(PLDM),intent(inout)::this
    call Note('Begin PLDM_kill.')
    this%initialized=.false.
    if(associated(this%qs))nullify(this%qs)
    if(associated(this%cs))nullify(this%cs)
    if(associated(this%cp))nullify(this%cp)
!!$    if(associated(this%q))nullify(this%q)
!!$    if(associated(this%p))nullify(this%p)
!!$    if(associated(this%qt))nullify(this%qt)
!!$    if(associated(this%pt))nullify(this%pt)
  end subroutine PLDM_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine PLDM_display(this)
    implicit none
    type(PLDM),intent(in)::this
    call Note('Begin PLDM_display.')
    if(check(this).EQ.1)then
       call warn('display_PLDM: object failed check.','displaying nothing.')
    else
       call display(this%qs)
       call display(this%cs)
       call display(this%cp)
       write(*,*)'ntraj=',this%ntraj
       write(*,*)'runtime(fs)=',this%runtime/fs
       write(*,*)'nuclear time step(fs)=',this%dtn/fs
       write(*,*)'electronic time step(fs)=',this%dte/fs
    end if
  end subroutine PLDM_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine PLDM_save(this,file)
    implicit none
    type(PLDM),intent(in)::this
    character*(*),intent(in)::file

    integer(long)::unit    

    call Note('Begin PLDM_save.')
    call Note('Input file= '//file)
    if(check(this).EQ.1)call stop('PLDM_save, failed object check.')

    unit=newunit()
    open(unit,file=file)
    write(unit,*)'PLDM'
    write(unit,*)this%ntraj
    write(unit,*)this%runtime
    write(unit,*)this%dtn
    write(unit,*)this%dte
    write(unit,*)this%dtout
    write(unit,*)this%dtsave
    close(unit)  
    
  end subroutine PLDM_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short) function PLDM_check(this)
    type(PLDM),intent(in)::this
    integer(long)::nstate
    PLDM_check=0

    call Note('Checking PLDM.')

    nstate=this%qs%hs%nstate

    !logical::initialized=.false.
    if(.not.this%initialized)then
       call Warn('PLDM not initialized.')
       PLDM_check=1
       return
    end if

    !type(quantum),pointer::qs
    if(check(this%qs).NE.0)then
       call Warn('PLDM%qs failed check.')
       PLDM_check=1
       return
    end if
    !type(classical),pointer::cs
    if(check(this%cs).NE.0)then
       call Warn('PLDM%cs failed check.')
       PLDM_check=1
       return
    end if
    !type(coupling),pointer::cp
    if(check(this%cp).NE.0)then
       call Warn('PLDM%cp failed check.')
       PLDM_check=1
       return
    end if

    !integer(long)::ntraj
    if(this%ntraj.NE.this%ntraj)then
       call Warn('PLDM%ntraj is NAN')
       PLDM_check=1
       return
    end if
    if(this%ntraj.LE.0)then
       call Warn('PLDM%ntraj is not a positive integer.')
       PLDM_check=1
       return
    end if
    if(this%ntraj.GT.huge(this%ntraj))then
       call Warn('PLDM%ntraj is huge.')
       PLDM_check=1
       return
    end if

    !real(double)::runtime
    if(this%runtime.NE.this%runtime)then
       call Warn('PLDM%runtime is NAN')
       PLDM_check=1
       return
    end if
    if(this%runtime.LE.0._double)then
       call Warn('PLDM%runtime is not a positive value.')
       PLDM_check=1
       return
    end if
    if(this%runtime.LT.epsilon(this%runtime))then
       call Warn('PLDM%runtime is tiny.')
       PLDM_check=1
       return
    end if
    if(this%runtime.GT.huge(this%runtime))then
       call Warn('PLDM%runtime is huge.')
       PLDM_check=1
       return
    end if

    !real(double)::dtn
    if(this%dtn.NE.this%dtn)then
       call Warn('PLDM%dtn is NAN')
       PLDM_check=1
       return
    end if
    if(this%dtn.LE.0._double)then
       call Warn('PLDM%dtn is not a positive value.')
       PLDM_check=1
       return
    end if
    if(this%dtn.LT.epsilon(this%dtn))then
       call Warn('PLDM%dtn is tiny.')
       PLDM_check=1
       return
    end if
    if(this%dtn.GT.huge(this%dtn))then
       call Warn('PLDM%dtn is huge.')
       PLDM_check=1
       return
    end if
    if(this%dtn.GT.this%runtime)then
       call Warn('PLDM%dtn is larger than PLDM%runtime.')
       PLDM_check=1
       return
    end if

    !real(double)::dte
    if(this%dte.NE.this%dte)then
       call Warn('PLDM%dte is NAN')
       PLDM_check=1
       return
    end if
    if(this%dte.LE.0._double)then
       call Warn('PLDM%dte is not a positive value.')
       PLDM_check=1
       return
    end if
    if(this%dte.LT.epsilon(this%dte))then
       call Warn('PLDM%dte is tiny.')
       PLDM_check=1
       return
    end if
    if(this%dte.GT.huge(this%dte))then
       call Warn('PLDM%dte is huge.')
       PLDM_check=1
       return
    end if
    if(this%dte.GT.this%dtn)then
       call Warn('PLDM%dte is larger than PLDM%dtn.')
       PLDM_check=1
       return
    end if

    !real(double)::dtout
    if(this%dtout.NE.this%dtout)then
       call Warn('PLDM%dtout is NAN')
       PLDM_check=1
       return
    end if
    if(this%dtout.LE.0._double)then
       call Warn('PLDM%dtout is not a positive value.')
       PLDM_check=1
       return
    end if
    if(this%dtout.LT.epsilon(this%dtout))then
       call Warn('PLDM%dtout is tiny.')
       PLDM_check=1
       return
    end if
    if(this%dtout.GT.huge(this%dtout))then
       call Warn('PLDM%dtout is huge.')
       PLDM_check=1
       return
    end if
    if(this%dtout.LT.this%dtn)then
       call Warn('PLDM%dtout is less than PLDM%dtn.')
       PLDM_check=1
       return
    end if
    !if((this%dtout/this%dtn)-floor(this%dtout/this%dtn).NE.0_double)then
    !   call Warn('PLDM%dtout is not an integer number of PLDM%dtn.')
    !   PLDM_check=1
    !   return
    !end if

    !real(double)::dtsave
    if(this%dtsave.NE.this%dtsave)then
       call Warn('PLDM%dtsave is NAN')
       PLDM_check=1
       return
    end if
    if(this%dtsave.LE.0._double)then
       call Warn('PLDM%dtsave is not a positive value.')
       PLDM_check=1
       return
    end if
    if(this%dtsave.LT.epsilon(this%dtsave))then
       call Warn('PLDM%dtsave is tiny.')
       PLDM_check=1
       return
    end if
    !if(this%dtsave.GT.86400._double)then
    !   call Warn('PLDM%dtsave is larger than 1 day.')
    !   PLDM_check=1
    !   return
    !end if

    !integer(long)::nbstep
    if(this%nbstep.LE.0)then
       call Warn('PLDM%nbstep is not a positive integer.')
       PLDM_check=1
       return
    end if
    if(this%nbstep.GT.huge(this%nbstep))then
       call Warn('PLDM%nbstep is huge.')
       PLDM_check=1
       return
    end if
    if(this%nbstep.NE.this%runtime/this%dtn)then
       call Warn('PLDM%nbstep is not snapped to grid.')
       PLDM_check=1
       return
    end if

    !integer(long)::nlstep
    if(this%nlstep.LE.0)then
       call Warn('PLDM%nlstep is not a positive integer.')
       PLDM_check=1
       return
    end if
    if(this%nlstep.GT.huge(this%nlstep))then
       call Warn('PLDM%nlstep is huge.')
       PLDM_check=1
       return
    end if
    if(this%nlstep.NE.this%dtn/this%dte)then
       call Warn('PLDM%nlstep is not snapped to grid.')
       PLDM_check=1
       return
    end if
    
    !integer(long)::nstep
    if(this%nstep.LE.0)then
       call Warn('PLDM%nstep is not a positive integer.')
       PLDM_check=1
       return
    end if
    if(this%nstep.GT.huge(this%nstep))then
       call Warn('PLDM%nstep is huge.')
       PLDM_check=1
       return
    end if
    if(this%runtime/this%dtout.NE.this%nstep)then
       call Warn('PLDM%nstep not snapped to grid.')
       PLDM_check=1
       return
    end if
    
  end function PLDM_check
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine PLDM_run(this,file,poponly,normstyle,forcesnap)
    implicit none
    type(PLDM),intent(inout)::this
    character*(*),intent(in),optional::file
    logical,intent(in),optional::poponly
    character*(*),intent(in),optional::normstyle
    logical,intent(in),optional::forcesnap

    !input parameters
    real(double)::runtime,tol,beta
    integer(long)::nbstep,nlit,ntraj,proctraj
    integer(long)::ndof,nstate,istate,istatet

    !Hamiltonian components
    real(double),allocatable::hel(:,:)
    real(double),allocatable::dhel(:,:,:)
    real(double)::x0i,p0i,x0it,p0it
    real(double),allocatable::fnon(:)
    complex(double),allocatable::initden(:,:)
    complex(double)::phase
    real(double)::counterterm

    !derived parameters
    integer(long)::nstep
    real(double)::dt,dtn,dtn2

    !dynamic arrays
    complex(double),allocatable:: redmat(:,:,:),redmat_sum(:,:,:),redmat_temp(:,:,:)
!    real(double),allocatable::x0(:),p0(:),xt0(:),pt0(:)

    !book keeping
    integer(long)::itime,iatt,itraj,unit,trajsum,norming,itraj_sum,jtraj
    integer(long)::i,j,k,n,nt,fs,fst 
    real(double)::norm,wnorm,w,cume
    logical::ierr
    integer(long)::timestart(8),timedif(8)
    real(double)::elapsedtime=0.
    character(len=path)::FMT

    integer(long)::redmatsize,bdof,Esize
    real(double)::time  !Current time                       

    real(double)::value
    integer(long)::ubound,lbound

    if(present(forcesnap))then
       if(forcesnap)then
          !snap output time step to grid
          value=this%runtime/this%dtout
          ubound=ceiling(value)
          lbound=floor(value)
          if(value.GE.(ubound+lbound)/2._double)then
             this%nstep=ubound
          else
             this%nstep=lbound
          end if
          this%dtout=this%runtime/real(this%nstep)
          
          !snap nuclear time step to grid
          value=this%runtime/this%dtn
          ubound=ceiling(value)
          lbound=floor(value)
          if(value.GE.(ubound+lbound)/2._double)then
             this%nbstep=ubound
          else
             this%nbstep=lbound
          end if
          this%dtn=this%runtime/real(this%nbstep)
          
          !snap electronic electronic time step to grid
          value=this%dtn/this%dte
          ubound=ceiling(value)
          lbound=floor(value)
          if(value.GE.(ubound+lbound)/2._double)then
             this%nlstep=ubound
          else
             this%nlstep=lbound
          end if
          this%dte=this%dtn/real(this%nlstep)
       end if
    end if
    !check after snapping to grid
    if(check(this).EQ.1)call stop('PLDM_run: input object failed check!')

    norming=1
    if(present(normstyle))then
       select case(normstyle)
       case('none')
          norming=0
       case('bytraj')
          norming=1
       case('bytrace')
          norming=2
       case('byinitden')
          norming=3
       case default
          call warn('PLDM_run:unknown normstyle','setting normstyle=none.')
       end select
    end if

    !input params
    nstate=this%qs%hs%nstate
    ndof=this%cs%hb%ndof
    bdof=ndof/nstate
    runtime=this%runtime!/ps
    nbstep=this%nbstep
    nlit=this%nlstep
    ntraj=this%ntraj
    beta=1._double/(kb*this%cs%hb%temperature)!/ps
    FMT="("//trim(int2str(nstate**2*2+1))//"(ES18.10E2,1X))"

    !hard coded params
    tol=0.

    !derived params
    nstep=nbstep*nlit
    dt=this%dte
    dtn=this%dtn
    dtn2=0.5*dtn*dtn

    !write(*,*)'runtime(ps)=',runtime/ps
    !write(*,*)'dtn(ps)=',dtn/ps
    !write(*,*)'dte(ps)=',dt/ps
    !write(*,*)'dtout(ps)=',dtout/ps


    if (nproc<ntraj)then
       proctraj=int(ntraj/nproc)
    else
       nproc=ntraj
       proctraj=1
    end if
    ntraj=proctraj*nproc

    !Allocate the dynamic arrays
    allocate(redmat(nstate,nstate,nbstep))
    allocate(hel(nstate,nstate),dhel(ndof,nstate,nstate))
    allocate(x0(nstate),p0(nstate),xt0(nstate),pt0(nstate))
    allocate(fnon(ndof))

    allocate(initden(nstate,nstate))
    allocate(redmat_sum(nstate,nstate,nbstep))
    allocate(redmat_temp(nstate,nstate,nbstep))

    !Set up potential parameters and constants
    initden=this%qs%hs%den
    wnorm=cnorm(initden)
!!$    wnorm=0.
!!$    do i=1,nstate
!!$       do j=1,nstate
!!$          wnorm=wnorm+sqrt(real(initden(i,j)&
!!$               *conjg(initden(i,j))))
!!$       end do
!!$    end do

    !initialize density matrix
    redmat=cmplx(0.,0.)

    !initialize trajectory
    do itraj=1,proctraj

       !sample input reduced density matrix
       w=ran0()*wnorm
       cume=0.
       k=0
       do while(cume.lt.w)
          i=k/nstate+1
          j=mod(k,nstate)+1
          
          cume=cume+sqrt(real(initden(i,j)&
               *conjg(initden(i,j))))

          istate=i
          istatet=j
          k=k+1
       end do
       if(istate.LT.1.or.istate.GT.nstate)call stop('PLDM_run: cannot properly sample initial density.')
       if(istatet.LT.1.or.istatet.GT.nstate)call stop('PLDM_run: cannot properly sample initial density.')
       call Note('trajectory: '//trim(int2str(itraj))//' initial state: '//trim(int2str(istate))// ','//trim(int2str(istate)))


       phase=initden(istate,istatet)/sqrt(initden(istate,istatet)*conjg(initden(istate,istatet)))

       !resample mapping variables
       do i=1,nstate
          x0(i)=gran()
          p0(i)=gran()
          xt0(i)=gran()
          pt0(i)=-gran()
       end do

       !focus initial mapping variables
       x0i=x0(istate)
       p0i=p0(istate)
       x0it=xt0(istatet)
       p0it=pt0(istatet) 


       !Initialize nuclear positions and momenta
       call resample(this%cs)
       call resample(this%cp)


       !!get Hamiltonian and force matrix
       hel=this%qs%hs%diabat+(this%cs%hb%V*iden(nstate))+real(this%cp%hc%V)
       dhel=real(this%cp%hc%dV)

       !Calculate Forces
       fnon=this%cs%hb%F
       do k=1,ndof
          do i=1,nstate
             fnon(k)=fnon(k)-0.25*dhel(k,i,i)*(x0(i)**2+p0(i)**2+xt0(i)**2+pt0(i)**2)
             do j=1,nstate
                if(j.ne.i) then
                   fnon(k)=fnon(k)-0.25*dhel(k,i,j)*(x0(i)*x0(j)+p0(i)*p0(j)+xt0(i)*xt0(j)+pt0(i)*pt0(j))
                end if
             end do
          end do
       end do

       !Main time stepping loop
       do itime=1,nbstep

          !Advance nuclear positions with current forces
          !first half Ver(let) step 
          this%cs%hb%Q=this%cs%hb%Q+this%cs%hb%rmass*this%cs%hb%P*dtn+this%cs%hb%rmass*fnon*dtn2
          this%cs%hb%P=this%cs%hb%P+0.5*fnon*dtn

          !Update Hamiltonian and force matrix
          call update(this%cs)
          call update(this%cp)
          hel=this%qs%hs%diabat+(this%cs%hb%V*iden(nstate))+real(this%cp%hc%V)
          dhel=real(this%cp%hc%dV)

          !Advance mapping variables 
          call mapverlet(this,nstate,nlit,dt,hel)!,x0,p0,xt0,pt0)

          !Calculate Forces
          fnon(:)=this%cs%hb%F
          do k=1,ndof
             do i=1,nstate
                fnon(k)=fnon(k)-0.25*dhel(k,i,i)*(x0(i)**2+p0(i)**2+xt0(i)**2+pt0(i)**2)
                do j=1,nstate
                   if(j.ne.i) then
                      fnon(k)=fnon(k)-0.25*dhel(k,i,j)*(x0(i)*x0(j)+p0(i)*p0(j)+xt0(i)*xt0(j)+pt0(i)*pt0(j))
                   end if
                end do
             end do
          end do

          !second half of (ver)Let
          this%cs%hb%P=this%cs%hb%P+0.5*fnon*dtn

          !Construct redmat
          do n=1,nstate
             do nt=1,nstate
                redmat(n,nt,itime)=redmat(n,nt,itime)&
                     +0.25*(x0(n)+eye*p0(n))&
                     *(xt0(nt)-eye*pt0(nt))*(x0i-eye*p0i)&
                     *(x0it+eye*p0it)&
                     *phase
             end do
          end do


       end do

       !restart files
       timedif=timearray()-timestart
       if(timedif(8).LT.0)then
          timedif(7)=timedif(7)-1
          timedif(8)=timedif(8)+1E3
       end if
       if(timedif(7).LT.0)then
          timedif(6)=timedif(6)-1
          timedif(7)=timedif(7)+60
       end if
       if(timedif(6).LT.0)then
          timedif(5)=timedif(5)-1
          timedif(6)=timedif(6)+60
       end if
       if(timedif(5).LT.0)then
          timedif(4)=timedif(4)-1
          timedif(5)=timedif(5)+24
       end if
       elapsedtime=timedif(8)*1E-3
       elapsedtime=elapsedtime+timedif(7)
       elapsedtime=elapsedtime+timedif(6)*60
       elapsedtime=elapsedtime+timedif(5)*1440
       
       if (elapsedtime.GE.this%dtsave)then
          
          call Prompt('Saving intermediate density '//normstyle//' approx ntraj= '&
               //trim(int2str(itraj*nproc)))
          
          !comunicate and put a barrier here
          redmatsize=2*size(redmat)
          if(myid.ne.master)then
!!$# ifdef MPI
!!$             call MPI_send(redmat,redmatsize,MPI_COMPLEX,master,100,&
!!$                  MPI_COMM_WORLD,ierr) !send redmat to master
!!$             call MPI_send(itraj,1,MPI_INTEGER,master,101,&
!!$                  MPI_COMM_WORLD,ierr) !send itraj to master
!!$# endif
          else
             redmat_sum=redmat
             itraj_sum=itraj
!!$# ifdef MPI
!!$             if(nproc.GT.1)then
!!$                do i=1,nproc-1
!!$                   call MPI_Recv(redmat_temp,redmatsize,MPI_COMPLEX,i,100,&
!!$                        MPI_COMM_WORLD,status,ierr) !some say not safe
!!$                   redmat_sum=redmat_sum+redmat_temp
!!$                   call MPI_Recv(jtraj,1,MPI_INTEGER,i,101,&
!!$                        MPI_COMM_WORLD,status,ierr) !some say not safe
!!$                   itraj_sum=itraj_sum+jtraj
!!$                   !write(*,*)i,jtraj,itraj_sum
!!$                end do
!!$             end if
!!$# endif
             
             if(present(file))then
                if(norming.EQ.1)norm=real(itraj_sum,double)
                if(norming.EQ.2)norm=trace(initden)*real(itraj_sum,double)
                if(norming.EQ.3)norm=cnorm(initden)*real(itraj_sum,double)/2._double
                if(norming.EQ.0)norm=1._double
                
                !call prompt(normstyle//' - '//trim(int2str(norming))//' - '//trim(float2str(norm)))
                
                unit=newunit()
                open(unit,file=file)
                
                ierr=.false.
                if(present(poponly))ierr=poponly
                if(ierr)then
                   do itime=1,nbstep
                      write(unit,FMT) itime*dtn,(redmat_sum(i,i,itime)/norm,i=1,nstate)
                   end do
                else
                   do itime=1,nbstep
                      write(unit,FMT) itime*dtn,((redmat_sum(i,j,itime)/norm,j=1,nstate),i=1,nstate)
                   end do
                end if
                
                close(unit)
             end if
             
             call Note('Overwritting quantum subsystem density with evolved Matrix.')
             this%qs%hs%den=redmat_sum(:,:,nbstep)
             
          end if
          
!!$# ifdef MPI
!!$          call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$# endif
          
          timestart=timearray()
       end if
       
    end do

!!$# ifdef MPI
!!$    call MPI_barrier(MPI_COMM_WORLD,ierr)
!!$# endif
    
    redmatsize=2*size(redmat)

    if(myid.ne.master)then
!!$# ifdef MPI
!!$       call MPI_send(redmat,redmatsize,MPI_COMPLEX,master,100,&
!!$            MPI_COMM_WORLD,ierr) !send redmat to master
!!$             call MPI_send(itraj,1,MPI_INTEGER,master,101,&
!!$                  MPI_COMM_WORLD,ierr) !send itraj to master
!!$# endif
    else
       redmat_sum=redmat
       itraj_sum=itraj
!!$# ifdef MPI
!!$       if(nproc.GT.1)then
!!$          do i=1,nproc-1
!!$             call MPI_Recv(redmat_temp,redmatsize,MPI_COMPLEX,i,100,&
!!$                  MPI_COMM_WORLD,status,ierr) !some say not safe
!!$             redmat_sum=redmat_sum+redmat_temp
!!$             call MPI_Recv(jtraj,1,MPI_INTEGER,i,101,&
!!$                  MPI_COMM_WORLD,status,ierr) !some say not safe
!!$             itraj_sum=itraj_sum+jtraj
!!$             !write(*,*)i,jtraj,itraj_sum
!!$          end do
!!$       end if
!!$# endif
       
       
       if(present(file))then
          if(norming.EQ.1)norm=real(itraj_sum,double)
          if(norming.EQ.2)norm=trace(initden)*real(itraj_sum,double)
          if(norming.EQ.3)norm=cnorm(initden)*real(itraj_sum,double)/2._double
          if(norming.EQ.0)norm=1._double

          !call prompt(normstyle//' - '//trim(int2str(norming))//' - '//trim(float2str(norm)))
             
          unit=newunit()
          open(unit,file=file)

          ierr=.false.
          if(present(poponly))ierr=poponly
          if(ierr)then
             do itime=1,nbstep
                write(unit,FMT) itime*dtn,(redmat_sum(i,i,itime)/norm,i=1,nstate)
             end do
          else
             do itime=1,nbstep
                write(unit,FMT) itime*dtn,((redmat_sum(i,j,itime)/norm,j=1,nstate),i=1,nstate)
             end do
          end if

          close(unit)
       end if

       call Note('Overwritting quantum subsystem density with evolved Matrix.')
       this%qs%hs%den=redmat_sum(:,:,nbstep)

    end if
    
    !Deallocate the dynamical arrays  
    deallocate(redmat)
    deallocate(hel,dhel)
    deallocate(x0,p0,xt0,pt0)
    deallocate(fnon)

    deallocate(initden)
    deallocate(redmat_sum)
    deallocate(redmat_temp)

  end subroutine PLDM_run
  
end module PLDM_class
