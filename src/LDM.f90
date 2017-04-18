!============================================================================
!>\brief
!! Lineraized Path Integral Propagators
!!\details
!! Supplies iterative and partial linearized density matrix propagators
!! ,ILDM and PLDM respectively.
!!\authors
!! Daniel Montemayor
!!\dates
!! August 2011
!!\bugs
!! No population relaxation is observed.
!!\todo
!! Reverse changes to qs%diabat made during RUN at the end of subroutine
!<============================================================================
!-----------------------------------------------------------------------------
!Dr Daniel Montemayor, Aug 2011, ACAM, UCD 
!-----------------------------------------------------------------------------
!----------------------C H A N G E  L O G------------------------------
! - dyn file only recorded by master processor
! + output hs%den is normalized by number of trajectories
! - output buf file is no longer divided by trajectories
! - allowed zero coupling tol
!-------------------v1.0.1 Sept 2011
!-----------------------------------------------------------------------------
module LDM_class
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

  Public::run,PLDM,ILDM!,LDM
  Public::new,kill,save,display,check

  type mappingH
     logical::initialized=.false.
     integer(long)::nstate
     real(double),dimension(:),pointer::q,p,qt,pt
  end type mappingH

  type redmat
     logical::initialized=.false.
     logical::normalize=.true.
     integer(long)::nstate
     integer(long)::ntraj
     real(double)::runtime
     real(double)::dtn
     real(double)::dte
     real(double)::dtout
     real(double)::dtsave
     real(double)::tol

     integer(long)::nbstep
     integer(long)::nlstep
     integer(long)::nstep

     !holds only as much of the reduced density matrix history
     !as needed by method
     complex(double),dimension(:,:,:,:),pointer::buf
     logical::bufok=.true.
  end type redmat

!!$  type LDM
!!$     logical::initialized=.false.
!!$     type(mappingH)::H
!!$     type(redmat)::den
!!$     type(quantum),pointer::qs
!!$     type(classical),pointer::cs
!!$     type(coupling),pointer::cp
!!$  end type LDM

  type PLDM
     logical::initialized=.false.
     type(mappingH)::H
     type(redmat)::den
     type(quantum),pointer::qs
     type(classical),pointer::cs
     type(coupling),pointer::cp
  end type PLDM

  type ILDM
     logical::initialized=.false.
     type(mappingH)::H
     type(redmat)::den
     type(quantum),pointer::qs
     type(classical),pointer::cs
     type(coupling),pointer::cp
     real(double)::shorttime
     integer(long)::nattempt
  end type ILDM

  interface new
     module procedure mappingH_new
     module procedure redmat_new
     !module procedure LDM_new
     module procedure PLDM_new
     module procedure ILDM_new
  end interface

  interface kill
     module procedure mappingH_kill
     module procedure redmat_kill
     !module procedure LDM_kill
     module procedure PLDM_kill
     module procedure ILDM_kill
  end interface

  interface save
     module procedure mappingH_save
     module procedure redmat_save
     !module procedure LDM_save
     module procedure PLDM_save
     module procedure ILDM_save
  end interface

  interface display
     module procedure mappingH_display
     module procedure redmat_display
     !module procedure LDM_display
     module procedure PLDM_display
     module procedure ILDM_display
  end interface

  interface check
     module procedure mappingH_check
     module procedure redmat_check
     !module procedure LDM_check
     module procedure PLDM_check
     module procedure ILDM_check
  end interface

  interface run
     !module procedure LDM_run
     module procedure PLDM_run
     module procedure ILDM_run
  end interface


contains
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine mappingH_new(this,nstate)
    implicit none
    type(mappingH),intent(inout)::this
    integer(long),intent(in),optional::nstate
    call Note('Begin new_mappingH.')
    
    this%nstate=0
    if(present(nstate))this%nstate=nstate
    do while(this%nstate.LT.1)
       write(*,*)'Enter number of quantum states.'
       read(*,*)this%nstate
       if(this%nstate.LT.1)then
          write(*,*)'Number of quantum sates must be a positive integer. try again!'
       end if
    end do
    
    if(associated(this%q))nullify(this%q)
    allocate(this%q(this%nstate))
    if(associated(this%p))nullify(this%p)
    allocate(this%p(this%nstate))
    if(associated(this%qt))nullify(this%qt)
    allocate(this%qt(this%nstate))
    if(associated(this%pt))nullify(this%pt)
    allocate(this%pt(this%nstate))
    
    this%q=0._double
    this%p=0._double
    this%qt=0._double
    this%pt=0._double

    this%initialized=.true.
  end subroutine mappingH_new
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine mappingH_kill(this)
    implicit none
    type(mappingH),intent(inout)::this
    call Note('Begin mappingH_kill.')
    this%initialized=.false.
    if(associated(this%q))nullify(this%q)
    if(associated(this%p))nullify(this%p)
    if(associated(this%qt))nullify(this%qt)
    if(associated(this%pt))nullify(this%pt)
  end subroutine mappingH_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine mappingH_display(this)
    implicit none
    type(mappingH),intent(in)::this
    call Note('Begin mappingH_display.')

    if(check(this).EQ.1)then
       call warn('display_mappingH: object failed check.','displaying nothing.')
    else
       write(*,*)'nstate=',this%nstate
       write(*,*)'q=',this%q
       write(*,*)'p=',this%p
       write(*,*)'qt=',this%qt
       write(*,*)'pt=',this%pt
    end if
  end subroutine mappingH_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine mappingH_save(this,file)
    implicit none
    type(mappingH),intent(in)::this
    character*(*),intent(in)::file

    integer(long)::unit    

    call Note('Begin mappingH_save.')
    call Note('input file='//file)
    if(check(this).EQ.1)call stop('mappingH_save, failed object check.')

    unit=newunit()
    open(unit,file=file)
    write(unit,*)'mappingH'
    write(unit,*)this%nstate
    close(unit)

  end subroutine mappingH_save
!------------------------------------------------------
  integer(short) function mappingH_check(this)
    type(mappingH),intent(in)::this
    mappingH_check=0

    call Note('Checking mappingH.')

    !logical::initialized=.false.
    if(.not.this%initialized)then
       call Warn('mappingH not initialized.')
       mappingH_check=1
       return
    end if

    !integer(long)::nstate
    if(this%nstate.NE.this%nstate)then
       call Warn('mappingH%nstate is NAN')
       mappingH_check=1
       return
    end if
    if(this%nstate.LE.0)then
       call Warn('mappingH%nstate is not a positive integer.')
       mappingH_check=1
       return
    end if
    if(this%nstate.GT.huge(this%nstate))then
       call Warn('mappingH%nstate is huge.')
       mappingH_check=1
       return
    end if

    !real(double),dimension(:),pointer::q
    if(.not.associated(this%q))then
       call Warn('mappingH%q is not associated.')
       mappingH_check=1
       return
    end if
    if(size(this%q).NE.this%nstate)then
       call Warn('mappingH%q has wrong size.')
       mappingH_check=1
       return
    end if
    if(any(this%q.NE.this%q))then
       call Warn('mappingH%q has NAN values.')
       mappingH_check=1
       return
    end if
    if(any(abs(this%q).GT.huge(this%q)))then
       call Warn('mappingH%q has huge values.')
       mappingH_check=1
       return
    end if

    !real(double),dimension(:),pointer::p
    if(.not.associated(this%p))then
       call Warn('mappingH%p is not associated.')
       mappingH_check=1
       return
    end if
    if(size(this%p).NE.this%nstate)then
       call Warn('mappingH%p has wrong size.')
       mappingH_check=1
       return
    end if
    if(any(this%p.NE.this%p))then
       call Warn('mappingH%p has NAN values.')
       mappingH_check=1
       return
    end if
    if(any(abs(this%p).GT.huge(this%p)))then
       call Warn('mappingH%p has huge values.')
       mappingH_check=1
       return
    end if

    !real(double),dimension(:),pointer::qt
    if(.not.associated(this%qt))then
       call Warn('mappingH%qt is not associated.')
       mappingH_check=1
       return
    end if
    if(size(this%qt).NE.this%nstate)then
       call Warn('mappingH%qt has wrong size.')
       mappingH_check=1
       return
    end if
    if(any(this%qt.NE.this%qt))then
       call Warn('mappingH%qt has NAN values.')
       mappingH_check=1
       return
    end if
    if(any(abs(this%qt).GT.huge(this%qt)))then
       call Warn('mappingH%qt has huge values.')
       mappingH_check=1
       return
    end if

    !real(double),dimension(:),pointer::pt
    if(.not.associated(this%pt))then
       call Warn('mappingH%pt is not associated.')
       mappingH_check=1
       return
    end if
    if(size(this%pt).NE.this%nstate)then
       call Warn('mappingH%pt has wrong size.')
       mappingH_check=1
       return
    end if
    if(any(this%pt.NE.this%pt))then
       call Warn('mappingH%pt has NAN values.')
       mappingH_check=1
       return
    end if
    if(any(abs(this%pt).GT.huge(this%pt)))then
       call Warn('mappingH%pt has huge values.')
       mappingH_check=1
       return
    end if

  end function mappingH_check
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine mapverlet(this,dt,nlstep,qs,cs,cp)
    implicit none

    type(mappingH),intent(inout)::this
    type(quantum),intent(in)::qs
    type(classical),intent(in)::cs
    type(coupling),intent(in)::cp
    real(double),intent(in)::dt
    integer(long),intent(in)::nlstep

    real(double),allocatable,dimension(:)::xrhs,prhs,xtrhs,ptrhs,force,forcet
    real(double),allocatable,dimension(:,:)::hel
    real(double):: xsumdum,psumdum
    integer(long):: i,n,m,nstate

    call Note('Begin mapverlet.')

    if(.not.this%initialized)call stop('mapverlet: mapping hamiltonian is not initialized.')

    nstate=this%nstate

    if(allocated(hel))deallocate(hel)
    if(allocated(xrhs))deallocate(xrhs)
    if(allocated(prhs))deallocate(prhs)
    if(allocated(xtrhs))deallocate(xtrhs)
    if(allocated(ptrhs))deallocate(ptrhs)
    if(allocated(force))deallocate(force)
    if(allocated(forcet))deallocate(forcet)
    allocate(hel(nstate,nstate))
    allocate(xrhs(nstate))
    allocate(prhs(nstate))
    allocate(xtrhs(nstate))
    allocate(ptrhs(nstate))
    allocate(force(nstate))
    allocate(forcet(nstate))

    hel=0.0_double
    do n=1,nstate
       do m=1,nstate
          hel(n,m)=hel(n,m)+cp%hc%V(n,m)+qs%hs%diabat(n,m)
       end do
       hel(n,n)=hel(n,n)+cs%hb%V
    end do

    do i=1,nlstep

       !     Generate current derivatives

       !     Forward

       do n=1,nstate
          xrhs(n)=hel(n,n)*this%p(n)
          prhs(n)=-hel(n,n)*this%q(n)
       end do

       do n=1,nstate
          xsumdum=0.
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                xsumdum=xsumdum+hel(n,m)*this%p(m)
                psumdum=psumdum+hel(n,m)*this%q(m)
             end if
          end do
          xrhs(n)=xrhs(n)+xsumdum
          prhs(n)=prhs(n)-psumdum
       end do

       !     Backward

       do n=1,nstate
          xtrhs(n)=hel(n,n)*this%pt(n)
          ptrhs(n)=-hel(n,n)*this%qt(n)
       end do

       do n=1,nstate
          xsumdum=0.
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                xsumdum=xsumdum+hel(n,m)*this%pt(m)
                psumdum=psumdum+hel(n,m)*this%qt(m)
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
          this%q(n)=this%q(n)+xrhs(n)*dt+0.5*force(n)*dt*dt
          this%p(n)=this%p(n)+0.5*prhs(n)*dt

          this%qt(n)=this%qt(n)+xtrhs(n)*dt+0.5*forcet(n)*dt*dt
          this%pt(n)=this%pt(n)+0.5*ptrhs(n)*dt
       end do

       !     Compute new first derivatives

       !     Forward

       do n=1,nstate
          prhs(n)=-hel(n,n)*this%q(n)
       end do

       do n=1,nstate
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                psumdum=psumdum+hel(n,m)*this%q(m)
             end if
          end do
          prhs(n)=prhs(n)-psumdum
       end do

       !     Backward

       do n=1,nstate
          ptrhs(n)=-hel(n,n)*this%qt(n)
       end do

       do n=1,nstate
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                psumdum=psumdum+hel(n,m)*this%qt(m)
             end if
          end do
          ptrhs(n)=ptrhs(n)-psumdum
       end do

       !    Advance let step

       do n=1,nstate
          this%p(n)=this%p(n)+0.5*prhs(n)*dt
          this%pt(n)=this%pt(n)+0.5*ptrhs(n)*dt
       end do

    end do

    if(allocated(hel))deallocate(hel)
    if(allocated(xrhs))deallocate(xrhs)
    if(allocated(prhs))deallocate(prhs)
    if(allocated(xtrhs))deallocate(xtrhs)
    if(allocated(ptrhs))deallocate(ptrhs)
    if(allocated(force))deallocate(force)
    if(allocated(forcet))deallocate(forcet)

    return
  end subroutine mapverlet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine redmat_new(this,nstate,file)
    implicit none
    type(redmat),intent(inout)::this
    character*(*),intent(in),optional::file
    integer(long),intent(in),optional::nstate
    character(len=title)::filetype
    character(len=path)::filename
    integer(long)::n

    integer(long)::unit

    call Note('Begin redmat_new.')
    if(present(file))&
         call Note('input file='//file)

    if(present(file))then
       unit=newunit()
       open(unit,file=file)
       read(unit,*)filetype
       filetype=adjustl(filetype)
       if(trim(filetype).NE.'redmat')&
            call stop('redmat_init error: not a valid input file.')
       read(unit,*)this%nstate
       read(unit,*)this%ntraj
       read(unit,*)this%runtime
       read(unit,*)this%dtn
       read(unit,*)this%dte
       read(unit,*)this%dtout
       read(unit,*)this%dtsave
       read(unit,*)this%tol
    else
       this%nstate=0
       if(present(nstate))then
          this%nstate=nstate
          call Note('input nstate='//trim(int2str(nstate)))
       end if
       do while(this%nstate.LT.1)
          write(*,*)'Enter number of quantum states of trajectories.'
          read(*,*)this%nstate
          if(this%nstate.LT.1)then
             write(*,*)'Number of quantum states  must be a positive integer. Try again!'
          end if
       end do
       this%ntraj=0
       do while(this%ntraj.LT.1)
          write(*,*)'Enter total number of trajectories.'
          read(*,*)this%ntraj
          if(this%ntraj.LT.1)then
             write(*,*)'Number of trajectories must be a positiveinteger . Try again!'
          end if
       end do
       this%runtime=0._double
       do while(this%runtime.LE.0._double)
          write(*,*)'Enter total runtime.'
          read(*,*)this%runtime
          if(this%runtime.LE.0._double)then
             write(*,*)'Tuntime must be positive. try again!'
          end if
       end do
       this%dtn=0._double
       do while(this%dtn.LE.0._double)
          write(*,*)'Enter nuclear time step.'
          read(*,*)this%dtn
          if(this%dtn.LE.0._double)then
             write(*,*)'Nuclear time step must positive. Try again!'
          end if
       end do
       this%dte=0._double
       do while(this%dte.LE.0._double)
          write(*,*)'Enter electronic time step.'
          read(*,*)this%dte
          if(this%dte.LE.0._double)then
             write(*,*)'Electronic time step must positive. Try again!'
          end if
       end do
       this%dtout=0._double
       do while(this%dtout.LE.0._double)
          write(*,*)'Enter output time step.'
          read(*,*)this%dtout
          if(this%dtout.LE.this%dtn)then
             write(*,*)'Output time step cannot be smaller than nuclear time step. Try again!'
          end if
       end do

       this%dtsave=0._double
       do while(this%dtsave.LE.0._double)
          write(*,*)'Enter time interval in wallclock seconds between Restart files.'
          read(*,*)this%dtsave
          if(this%dtsave.LE.0._double)then
             write(*,*)'Save time step must positive. Try again!'
          end if
       end do
       this%tol=-1._double
       do while(this%tol.LT.0._double)
          write(*,*)'Coupling threashold(Off-diagonal elements with absolute energy less than will not be propagated).'
          read(*,*)this%tol
          if(this%tol.LT.0._double)then
             write(*,*)'Threashold can not be negative. Try again.'
          end if
       end do
    end if

    this%nbstep=ceiling(this%runtime/this%dtn)
    this%dtn=this%runtime/real(this%nbstep)
    this%nstep=ceiling(this%dtout/this%dtn)
    do while(mod(this%nbstep,this%nstep).NE.0)
       this%nstep=this%nstep-1
    end do
    this%dtout=this%nstep*this%dtn 
    this%nstep=int(this%nbstep/this%nstep)
    if(associated(this%buf))nullify(this%buf)
    allocate(this%buf(this%nstate,this%nstate,0:this%nstep,0:1))
    this%buf=(0._double,0._double)
    this%bufok=.true.

    !electronic time step
    this%nlstep=ceiling(this%dtn/this%dte)
    this%dte=this%dtn/real(this%nlstep)

    if(present(file))close(unit)
    this%initialized=.true.
    if(check(this).EQ.1)call stop('redmat_new: failed check')
  end subroutine redmat_new
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine redmat_kill(this)
    implicit none
    type(redmat),intent(inout)::this
    call Note('Begin redmat_kill.')
    this%initialized=.false.
    if(associated(this%buf))nullify(this%buf)
  end subroutine redmat_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine redmat_display(this)
    implicit none
    type(redmat),intent(in)::this
    call Note('Begin redmat_display.')
    if(check(this).EQ.1)then
       call warn('display_redmat: object failed check.','displaying nothing.')
    else
       write(*,*)'nstate=',this%nstate
       write(*,*)'ntraj=',this%ntraj
       write(*,*)'runtime(fs)=',this%runtime/fs
       write(*,*)'nuclear time step(fs)=',this%dtn/fs
       write(*,*)'electronic time step(fs)=',this%dte/fs
       write(*,*)'Coupling threashold(1/cm)=',this%tol/invcm
    end if
  end subroutine redmat_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine redmat_save(this,file)
    implicit none
    type(redmat),intent(in)::this
    character*(*),intent(in)::file

    integer(long)::unit    

    call Note('Begin redmat_save.')
    call Note('Input file= '//file)
    if(check(this).EQ.1)call stop('redmat_save, failed object check.')
    unit=newunit()
    open(unit,file=file)
    write(unit,*)'redmat'
    write(unit,*)this%nstate
    write(unit,*)this%ntraj
    write(unit,*)this%runtime
    write(unit,*)this%dtn
    write(unit,*)this%dte
    write(unit,*)this%dtout
    write(unit,*)this%dtsave
    write(unit,*)this%tol
    close(unit)  

  end subroutine redmat_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short) function redmat_check(this)
    type(redmat),intent(in)::this
    redmat_check=0

    call Note('Checking redmat.')

    !logical::initialized=.false.
    if(.not.this%initialized)then
       call Warn('redmat not initialized.')
       redmat_check=1
       return
    end if

    !integer(long)::nstate
    if(this%nstate.NE.this%nstate)then
       call Warn('redmat%nstate is NAN')
       redmat_check=1
       return
    end if
    if(this%nstate.LE.0)then
       call Warn('redmat%nstate is not a positive integer.')
       redmat_check=1
       return
    end if
    if(this%nstate.GT.huge(this%nstate))then
       call Warn('redmat%nstate is huge.')
       redmat_check=1
       return
    end if

    !integer(long)::ntraj
    if(this%ntraj.NE.this%ntraj)then
       call Warn('redmat%ntraj is NAN')
       redmat_check=1
       return
    end if
    if(this%ntraj.LE.0)then
       call Warn('redmat%ntraj is not a positive integer.')
       redmat_check=1
       return
    end if
    if(this%ntraj.GT.huge(this%ntraj))then
       call Warn('redmat%ntraj is huge.')
       redmat_check=1
       return
    end if

    !real(double)::runtime
    if(this%runtime.NE.this%runtime)then
       call Warn('redmat%runtime is NAN')
       redmat_check=1
       return
    end if
    if(this%runtime.LE.0._double)then
       call Warn('redmat%runtime is not a positive value.')
       redmat_check=1
       return
    end if
    if(this%runtime.LT.epsilon(this%runtime))then
       call Warn('redmat%runtime is tiny.')
       redmat_check=1
       return
    end if
    if(this%runtime.GT.huge(this%runtime))then
       call Warn('redmat%runtime is huge.')
       redmat_check=1
       return
    end if

    !real(double)::dtn
    if(this%dtn.NE.this%dtn)then
       call Warn('redmat%dtn is NAN')
       redmat_check=1
       return
    end if
    if(this%dtn.LE.0._double)then
       call Warn('redmat%dtn is not a positive value.')
       redmat_check=1
       return
    end if
    if(this%dtn.LT.epsilon(this%dtn))then
       call Warn('redmat%dtn is tiny.')
       redmat_check=1
       return
    end if
    if(this%dtn.GT.huge(this%dtn))then
       call Warn('redmat%dtn is huge.')
       redmat_check=1
       return
    end if
    if(this%dtn.GT.this%runtime)then
       call Warn('redmat%dtn is larger than redmat%runtime.')
       redmat_check=1
       return
    end if

    !real(double)::dte
    if(this%dte.NE.this%dte)then
       call Warn('redmat%dte is NAN')
       redmat_check=1
       return
    end if
    if(this%dte.LE.0._double)then
       call Warn('redmat%dte is not a positive value.')
       redmat_check=1
       return
    end if
    if(this%dte.LT.epsilon(this%dte))then
       call Warn('redmat%dte is tiny.')
       redmat_check=1
       return
    end if
    if(this%dte.GT.huge(this%dte))then
       call Warn('redmat%dte is huge.')
       redmat_check=1
       return
    end if
    if(this%dte.GT.this%dtn)then
       call Warn('redmat%dte is larger than redmat%dtn.')
       redmat_check=1
       return
    end if

    !real(double)::dtout
    if(this%dtout.NE.this%dtout)then
       call Warn('redmat%dtout is NAN')
       redmat_check=1
       return
    end if
    if(this%dtout.LE.0._double)then
       call Warn('redmat%dtout is not a positive value.')
       redmat_check=1
       return
    end if
    if(this%dtout.LT.epsilon(this%dtout))then
       call Warn('redmat%dtout is tiny.')
       redmat_check=1
       return
    end if
    if(this%dtout.GT.huge(this%dtout))then
       call Warn('redmat%dtout is huge.')
       redmat_check=1
       return
    end if
    if(this%dtout.LT.this%dtn)then
       call Warn('redmat%dtout is less than redmat%dtn.')
       redmat_check=1
       return
    end if
    if((this%dtout/this%dtn)-floor(this%dtout/this%dtn).NE.0_double)then
       call Warn('redmat%dtout is not an integer number of redmat%dtn.')
       redmat_check=1
       return
    end if

    !real(double)::dtsave
    if(this%dtsave.NE.this%dtsave)then
       call Warn('redmat%dtsave is NAN')
       redmat_check=1
       return
    end if
    if(this%dtsave.LE.0._double)then
       call Warn('redmat%dtsave is not a positive value.')
       redmat_check=1
       return
    end if
    if(this%dtsave.LT.epsilon(this%dtsave))then
       call Warn('redmat%dtsave is tiny.')
       redmat_check=1
       return
    end if
    if(this%dtsave.GT.86400._double)then
       call Warn('redmat%dtsave is larger than 1 day.')
       redmat_check=1
       return
    end if

    !real(double)::tol
    if(this%tol.NE.this%tol)then
       call Warn('redmat%tol is NAN')
       redmat_check=1
       return
    end if
    if(abs(this%tol).GT.huge(this%tol))then
       call Warn('redmat%tol is huge.')
       redmat_check=1
       return
    end if

    !integer(long)::nbstep
    if(this%nbstep.LE.0)then
       call Warn('redmat%nbstep is not a positive integer.')
       redmat_check=1
       return
    end if
    if(this%nbstep.GT.huge(this%nbstep))then
       call Warn('redmat%nbstep is huge.')
       redmat_check=1
       return
    end if

    !integer(long)::nlstep
    if(this%nlstep.LE.0)then
       call Warn('redmat%nlstep is not a positive integer.')
       redmat_check=1
       return
    end if
    if(this%nlstep.GT.huge(this%nlstep))then
       call Warn('redmat%nlstep is huge.')
       redmat_check=1
       return
    end if

    !integer(long)::nstep
    if(this%nstep.LE.0)then
       call Warn('redmat%nstep is not a positive integer.')
       redmat_check=1
       return
    end if
    if(this%nstep.GT.huge(this%nstep))then
       call Warn('redmat%nstep is huge.')
       redmat_check=1
       return
    end if

    !complex(double),dimension(:,:,:),pointer::buf
    if(.not.associated(this%buf))then
       call Warn('redmat%buf is not associated.')
       redmat_check=1
       return
    end if
    if(size(this%buf,1).NE.this%nstate)then
       call Warn('redmat%buf does not have proper size.')
       redmat_check=1
       return
    end if
    if(size(this%buf,2).NE.this%nstate)then
       call Warn('redmat%buf does not have proper size.')
       redmat_check=1
       return
    end if
    if(size(this%buf,3).NE.this%nstep+1)then
       call Warn('redmat%buf does not have proper size.')
       redmat_check=1
       return
    end if
    if(size(this%buf,4).NE.2)then
       call Warn('redmat%buf does not have proper size.')
       redmat_check=1
       return
    end if
    if(any(this%buf.NE.this%buf))then
       call Warn('redmat%buf has NAN values.')
       redmat_check=1
       return
    end if
    if(any(abs(this%buf).GT.huge(this%tol)))then
       call Warn('redmat%buf has huge values.')
       redmat_check=1
       return
    end if

  end function redmat_check
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
 

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

    call new(this%H,this%qs%hs%nstate)

    if(present(file))then
       unit=newunit()
       open(unit,file=file)
       read(unit,*)filetype
       filetype=adjustl(filetype)
       if(trim(filetype).NE.'LDM')&
            call stop('PLDM_init error: not a valid input file.')
       read(unit,*)filename
       filename=adjustl(filename)
       call new(this%den,nstate=this%qs%hs%nstate,file=trim(filename))
    else
       call new(this%den,nstate=this%qs%hs%nstate)
    end if

    if(present(file))close(unit)
    this%initialized=.true.
    if(check(this).EQ.1)call stop('PLDM_new, failed final object check.')
  end subroutine PLDM_new
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine PLDM_kill(this)
    implicit none
    type(PLDM),intent(inout)::this
    call Note('Begin PLDM_kill.')
    this%initialized=.false.
    call kill(this%H)
    call kill(this%den)
    if(associated(this%qs))nullify(this%qs)
    if(associated(this%cs))nullify(this%cs)
    if(associated(this%cp))nullify(this%cp)
  end subroutine PLDM_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine PLDM_display(this)
    implicit none
    type(PLDM),intent(in)::this
    call Note('Begin PLDM_display.')
    if(check(this).EQ.1)then
       call warn('display_PLDM: object failed check.','displaying nothing.')
    else
       call display(this%H)
       call display(this%den)
       call display(this%qs)
       call display(this%cs)
       call display(this%cp)
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
    write(unit,*)'LDM'
    write(unit,*)quote(file//'.den')
    call save(this%den,file=file//'.den')
    write(unit,*)this%den%dtn
    close(unit)  

  end subroutine PLDM_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short) function PLDM_check(this)
    type(PLDM),intent(in)::this
    PLDM_check=0

    call Note('Checking PLDM.')

    !logical::initialized=.false.
    if(.not.this%initialized)then
       call Warn('PLDM not initialized.')
       PLDM_check=1
       return
    end if

    if(check(this%H).NE.0)then
       call Warn('PLDM%H failed check.')
       PLDM_check=1
       return
    end if

    if(check(this%den).NE.0)then
       call Warn('PLDM%den failed check.')
       PLDM_check=1
       return
    end if

    !type(quantum),pointer::qs
    if(.not.associated(this%qs))then
       call Warn('PLDM%qs is not associated.')
       PLDM_check=1
       return
    end if
    if(check(this%qs).NE.0)then
       call Warn('PLDM%qs failed check.')
       PLDM_check=1
       return
    end if
    if(this%H%nstate.NE.this%qs%hs%nstate)then
       call Warn('PLDM%H%nstate.NE.PLDM%qs%hs%nstate')
       PLDM_check=1
       return
    end if
    if(this%den%nstate.NE.this%qs%hs%nstate)then
       call Warn('PLDM%den%nstate.NE.PLDM%qs%hs%nstate')
       PLDM_check=1
       return
    end if

    !type(classical),pointer::cs
    if(.not.associated(this%cs))then
       call Warn('PLDM%cs is not associated.')
       PLDM_check=1
       return
    end if
    if(check(this%cs).NE.0)then
       call Warn('PLDM%cs failed check.')
       PLDM_check=1
       return
    end if

    !type(coupling),pointer::cp
    if(.not.associated(this%cp))then
       call Warn('PLDM%cp is not associated.')
       PLDM_check=1
       return
    end if
    if(check(this%cp).NE.0)then
       call Warn('PLDM%cp failed check.')
       PLDM_check=1
       return
    end if

  end function PLDM_check
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine PLDM_run(this,file,normalize)
    implicit none
    type(PLDM),intent(inout)::this
    character*(*),intent(in),optional::file
    logical,intent(in),optional::normalize

    integer(long)::timestart(8),timedif(8)
    real(double)::elapsedtime

    integer(long)::unit

    real(double):: x0i,p0i,x0it,p0it

    complex*16,allocatable::bigA(:,:)
    real(double),allocatable::fnon(:)

    real(double)::QE,CE
    real(double)::dtn,dtn2,dtdump,norm
    integer(long)::nstate,ndof,ntraj
    integer(long)::itraj,istep,ndump
    integer(long)::istate,istatet
    integer(long)::i,j

    logical::recdyn
    character(len=path)::FMT
 

    integer,parameter::sp=100
    real*8::stat(sp),S0,Smin,Smax,dS,S,param


   if(check(this).EQ.1)call stop('PLDM_run, failed object check.')
    call Note('Begin PLDM_run.')
    if(present(file))call Note('Input file= '//file)
    FMT="("//trim(int2str(this%qs%hs%nstate**2*2+1))//"(ES18.10E2,1X))"
    timestart=timearray()


    dtn=this%den%dtn
    dtn2=0.5_double*dtn*dtn
    nstate=this%qs%hs%nstate
    ndof=this%cs%hb%ndof
    ntraj=this%den%ntraj
    dtdump=this%den%dtout
    this%den%normalize=.true.
    if(present(normalize))this%den%normalize=normalize

    if(allocated(bigA))deallocate(bigA)
    allocate(bigA(nstate,nstate))
    if(allocated(fnon))deallocate(fnon)
    allocate(fnon(ndof))
    !write(*,*)'PLDM'
    !write(*,*)this%den%dtout,this%den%dtn,this%den%dte
    !write(*,*)this%den%nstep,this%den%nbstep,this%den%nlstep
    !write(*,*)

    istate=1
    istatet=1
    call warn('run_PLDM: initial state is hard coded: 1, 1')
  
    !initialize reduced density matrix
    this%den%buf=cmplx(0.,0.)
    this%den%bufok=.true.

    if(myid.EQ.master)then

       unit=newunit()
       open(unit,file='fort.2000')

       Smin=0.
       Smax=400.
       dS=(Smax-Smin)/real(sp,8)
       S0=smin
       stat=0.
       do itraj=1,2000

          call resample(this%cs)
          call resample(this%cp)

          !initialize mapping variables
          do i=1,nstate
             this%H%q(i)=gran()
             this%H%p(i)=gran()
             this%H%qt(i)=gran()
             this%H%pt(i)=gran()
          end do
          x0i=this%H%q(istate)
          p0i=this%H%p(istate)
          x0it=this%H%qt(istatet)
          p0it=this%H%pt(istatet)

          !Calculate nonadiabatic Forces
          fnon=this%cs%hb%F
          do i=1,nstate
             do j=1,nstate
                fnon=fnon-0.25*this%cp%hc%dV(:,i,j)&
                     *(this%H%q(i)*this%H%q(j)+this%H%p(i)*this%H%p(j)&
                     +this%H%qt(i)*this%H%qt(j)+this%H%pt(i)*this%H%pt(j))
             end do
          end do


!!$        !Calculate Forces
!!$        !---Diagonal contribution
!!$        fnon(:)=-omega(:)**2*qs(:)
!!$        
!!$        do k=1,ndof
!!$           do i=1,nstate
!!$              fnon(k)=fnon(k)-0.25*dhel(k,i,i)*(x0(i)**2+p0(i)**2+xt0(i)**2+pt0(i)**2)
!!$              do j=1,nstate
!!$                 if(j.ne.i) then
!!$                    fnon(k)=fnon(k)-0.25*dhel(k,i,j)*(x0(i)*x0(j)+p0(i)*p0(j)+xt0(i)*xt0(j)+pt0(i)*pt0(j))
!!$                 end if
!!$              end do
!!$           end do
!!$        end do
!!$



          param=(maxval(this%cs%hb%mode))*ps
          !param=(this%cs%hb%F(1))

          if(itraj.EQ.1)smin=param
          if(itraj.EQ.1)smax=param
       
          if(param.LE.smin)smin=param
          if(param.GE.smax)smax=param
          
          !accumulate stats
          do i=1,sp
             S=((i-1)*dS+S0)
             do j=1,this%cs%hb%ndof
                param=this%cs%hb%mode(j)*ps
                if(param.GT.S.and.param.LE.S+dS)&
                     stat(i)=stat(i)+1.
             end do
             !if(param.GT.S.and.param.LE.S+dS)stat(i)=stat(i)+1.
          end do
       end do
       write(*,*)'smin=',smin
       write(*,*)'smax=',smax
       do i=1,sp
          S=(i-1)*dS+S0 
          write(unit,*)S,stat(i)/maxval(stat)
       end do
    end if
stop

    !open ECON dyn file
    recdyn=.false.
    if(present(file).and.(myid.EQ.master))recdyn=.true.
    if(recdyn)then
       unit=newunit()
       open(unit,file=file//'.dyn')
    end if

    !loop over trajectories
    itraj=0
    do while(itraj.LT.ntraj)

       this%den%bufok=.true.

       !Resample system subsystem
       call resample(this%cs)
       call resample(this%cp)

       !initialize mapping variables
       do i=1,nstate
          this%H%q(i)=gran()
          this%H%p(i)=gran()
          this%H%qt(i)=gran()
          this%H%pt(i)=gran()
       end do
       x0i=this%H%q(istate)
       p0i=this%H%p(istate)
       x0it=this%H%qt(istatet)
       p0it=this%H%pt(istatet)

       !Calculate nonadiabatic Forces
       fnon=this%cs%hb%F
       do i=1,nstate
          do j=1,nstate
             fnon=fnon-0.25*this%cp%hc%dV(:,i,j)&
                  *(this%H%q(i)*this%H%q(j)+this%H%p(i)*this%H%p(j)&
                  +this%H%qt(i)*this%H%qt(j)+this%H%pt(i)*this%H%pt(j))
          end do
       end do
       if(any(fnon.NE.fnon))then
          call warn('PLDM_run: NAN nondaiabatic force.'&
               ,'thowing away trajectory.')
          this%den%bufok=.false.
          goto 31
       end if
       if(any(abs(fnon).GE.huge(fnon)))then
          call warn('PLDM_run: huge nondaiabatic force.'&
               ,'thowing away trajectory.')
          this%den%bufok=.false.
          goto 31
       end if

       !Construct bigA
       do i=1,nstate
          do j=1,nstate
             bigA(i,j)=0.25*(this%H%q(i)+eye*this%H%p(i))&
                  *(this%H%qt(j)-eye*this%H%pt(j))&
                  *(x0i-eye*p0i)*(x0it+eye*p0it)
          end do
       end do

       !Main time stepping loop
       ndump=0
       this%den%buf(:,:,ndump,1)=bigA+this%den%buf(:,:,ndump,1)
       do istep=1,this%den%nbstep 

          !write(*,*)itraj,istep,this%den%nbstep

          !Advance nuclear positions with current forces
          !first half Ver(let) step 
          this%cs%hb%Q=this%cs%hb%Q&
               +this%cs%hb%rmass*this%cs%hb%P*dtn&
               +this%cs%hb%rmass*fnon*dtn2
          this%cs%hb%P=this%cs%hb%P&
               +0.5*fnon*dtn
          
          !Compute New force 
          call update(this%cs)
          call update(this%cp)
          
          !Advance mapping variables 
          call mapverlet(this%H,this%den%dte,this%den%nlstep&
               ,this%qs,this%cs,this%cp)

          !Calculate nonadiabatic Forces
          fnon=this%cs%hb%F
          do i=1,nstate
             do j=1,nstate
                fnon=fnon-0.25*this%cp%hc%dV(:,i,j)&
                     *(this%H%q(i)*this%H%q(j)+this%H%p(i)*this%H%p(j)&
                     +this%H%qt(i)*this%H%qt(j)+this%H%pt(i)*this%H%pt(j))
             end do
          end do
          if(any(fnon.NE.fnon))then
             call warn('PLDM_run: NAN nondaiabatic force.'&
                  ,'thowing away trajectory.')
             this%den%bufok=.false.
             goto 31
          end if
          if(any(abs(fnon).GE.huge(fnon)))then
             call warn('PLDM_run: huge nondaiabatic force.'&
                  ,'thowing away trajectory.')
             this%den%bufok=.false.
             goto 31
          end if
          
          !second half of (ver)Let
          this%cs%hb%P=this%cs%hb%P+0.5*fnon*dtn

          !Construct bigA
          do i=1,nstate
             do j=1,nstate
                bigA(i,j)=0.25*(this%H%q(i)+eye*this%H%p(i))&
                     *(this%H%qt(j)-eye*this%H%pt(j))&
                     *(x0i-eye*p0i)*(x0it+eye*p0it)
             end do
          end do

          !record output
          if(istep*dtn.GE.(ndump+1)*dtdump)then
             ndump=ndump+1
             this%den%buf(:,:,ndump,1)=bigA+this%den%buf(:,:,ndump,1)

             !optional dyn record
             if(itraj.LE.1.and.recdyn)then
                call update(this%cs)
                call update(this%cp)
                QE=0._double
                CE=0._double
                do i=1,nstate
                   QE=QE+this%qs%hs%diabat(i,i)*bigA(i,i)
                   CE=CE+this%cp%hc%V(i,i)*bigA(i,i)
                end do
                write(unit,"(20(ES18.10E2,1X))")&
                     istep*dtn&
                     ,this%qs%hs%diabat&
                     ,this%cs%hb%V&
                     ,this%cs%hb%T&
                     ,this%cp%hc%V&
                     ,this%cs%hb%F(1)&
                     ,this%cp%hc%dV(1,:,:)&
                     ,sum(this%H%q**2+this%H%p**2)&
                     ,sum(this%H%qt**2+this%H%pt**2)&
                     ,this%cs%hb%Q(1),this%cs%hb%P(1)
             end if

          end if

       end do!End of Nuclear Stepping loop
       if(ndump.NE.this%den%nstep)&
            call stop('PLDM_run: mismatch in output count'&
            ,'make sure dtout is an integer number of nuclear timesteps.')

31     continue       
       !save buf if traj passed 
       if(this%den%bufok)then
          this%den%buf(:,:,:,0)=this%den%buf(:,:,:,0)+this%den%buf(:,:,:,1)
          itraj=itraj+1
       end if
       write(1234,*)real(this%den%buf(1,1,100,0))/real(itraj)       

       if(recdyn.and.itraj.GT.1)then
          close(unit) !close dyn file
       else
          write(unit,*)!add blank line
       end if

       if(present(file))then
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
          
          if (elapsedtime.GE.this%den%dtsave)then
             call Prompt('Saving intermediate density ntraj= '&
                  //trim(int2str(itraj)))

             norm=1._double
             if(this%den%normalize)norm=real(itraj)

             unit=newunit()
             open(unit,file=file)
             do istep=0,this%den%nstep
                write(unit,trim(FMT))&
                     istep*dtdump,((real(this%den%buf(i,j,istep,1))/norm&
                     ,aimag(this%den%buf(i,j,istep,1))/norm&
                     ,j=1,nstate),i=1,nstate)
             end do
             close(unit)
             timestart=timearray()
             this%qs%hs%den=this%den%buf(:,:,this%den%nstep,1)/real(itraj)
          end if
       end if
    end do!end ntraj loop
    
    if(present(file))then
       norm=1._double
       if(this%den%normalize)norm=real(itraj)

       unit=newunit()
       open(unit,file=file)
       do istep=0,this%den%nstep
          write(unit,trim(FMT))istep*dtdump,&
               ((real(this%den%buf(i,j,istep,0))/norm&
               &,aimag(this%den%buf(i,j,istep,0))/norm&
               ,j=1,nstate),i=1,nstate)
       end do
       close(unit)
    end if
    
    this%qs%hs%den=this%den%buf(:,:,this%den%nstep,0)/real(ntraj)
    
    if(allocated(bigA))deallocate(bigA)
    if(allocated(fnon))deallocate(fnon)

  end subroutine PLDM_run
 !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine ILDM_new(this,qs,cs,cp,file)
    implicit none
    type(ILDM),intent(inout)::this
    type(quantum),intent(in),target::qs
    type(classical),intent(in),target::cs
    type(coupling),intent(in),target::cp
    character*(*),intent(in),optional::file
    character(len=title)::filetype
    character(len=path)::filename
    integer(long)::unit

    call Note('Begin ILDM_new.')
    if(present(file))call Note('Input file= '//file)
    if(check(qs).EQ.1)call stop('ILDM_new, failed qs check')
    if(check(cs).EQ.1)call stop('ILDM_new, failed cs check')
    if(check(cp).EQ.1)call stop('ILDM_new, failed cp check')

    if(associated(this%qs))nullify(this%qs)
    this%qs=>qs
    if(associated(this%cs))nullify(this%cs)
    this%cs=>cs
    if(associated(this%cp))nullify(this%cp)
    this%cp=>cp

    call new(this%H,this%qs%hs%nstate)

    if(present(file))then
       unit=newunit()
       open(unit,file=file)
       read(unit,*)filetype
       filetype=adjustl(filetype)
       if(trim(filetype).NE.'LDM')&
            call stop('ILDM_init error: not a valid input file.')
       read(unit,*)filename
       filename=adjustl(filename)
       call new(this%den,nstate=this%qs%hs%nstate,file=trim(filename))
    else
       call new(this%den,nstate=this%qs%hs%nstate)
    end if

    if(present(file))then
       read(unit,*)this%shorttime
    else
       this%shorttime=this%den%runtime*2.0_double
       do while(this%shorttime.GT.this%den%runtime)
          write(*,*)'Enter total short time approx.'
          read(*,*)this%shorttime
          if(this%shorttime.GT.this%den%runtime)then
             write(*,*)'short time must be < or = to runtime. try again!'
          end if
       end do
    end if

    this%nattempt=ceiling(this%den%runtime/this%shorttime)
    this%shorttime=this%den%runtime/real(this%nattempt)
    this%den%nbstep=ceiling(this%shorttime/this%den%dtn)
    this%den%dtn=this%shorttime/real(this%den%nbstep)
    this%den%nstep=ceiling(this%den%dtout/this%den%dtn)
    do while(mod(this%den%nbstep,this%den%nstep).NE.0)
       this%den%nstep=this%den%nstep-1
    end do
    this%den%dtout=this%den%nstep*this%den%dtn
    this%den%nstep=int(this%den%nbstep/this%den%nstep)*this%nattempt
    if(associated(this%den%buf))nullify(this%den%buf)
    allocate(this%den%buf(this%den%nstate,this%den%nstate,0:this%den%nstep,0:1))
    this%den%buf=(0._double,0._double)
    this%den%bufok=.true.

    !electronic time step
    this%den%nlstep=ceiling(this%den%dtn/this%den%dte)
    this%den%dte=this%den%dtn/real(this%den%nlstep)

    if(present(file))close(unit)
    this%initialized=.true.
    if(check(this).NE.0)call stop('ILDM_new: failed final object check')

  end subroutine ILDM_new
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine ILDM_kill(this)
    implicit none
    type(ILDM),intent(inout)::this
    call Note('Begin ILDM_kill.')
    this%initialized=.false.
    call kill(this%H)
    call kill(this%den)
    if(associated(this%qs))nullify(this%qs)
    if(associated(this%cs))nullify(this%cs)
    if(associated(this%cp))nullify(this%cp)
  end subroutine ILDM_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine ILDM_display(this)
    implicit none
    type(ILDM),intent(in)::this
    call Note('Begin ILDM_display.')
    if(check(this).NE.0)then
       call warn('display_ILDM: object failed check.','displaying nothing.')
    else
       call display(this%H)
       call display(this%den)
       write(*,*)'shorttime approx=',this%shorttime
       write(*,*)'number of iterations=',this%nattempt
       call display(this%qs)
       call display(this%cs)
       call display(this%cp)
    end if
  end subroutine ILDM_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine ILDM_save(this,file)
    implicit none
    type(ILDM),intent(in)::this
    character*(*),intent(in)::file

    integer(long)::unit    

    call Note('Begin ILDM_save.')
    call Note('Input file= '//file)
    if(check(this).EQ.1)call stop('ILDM_save, failed object check')

    unit=newunit()
    open(unit,file=file)
    write(unit,*)'LDM'
    write(unit,*)quote(file//'.den')
    call save(this%den,file=file//'.den')
    write(unit,*)this%shorttime
    close(unit)  

  end subroutine ILDM_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short) function ILDM_check(this)
    type(ILDM),intent(in)::this
    ILDM_check=0

    call Note('Checking ILDM.')

    !logical::initialized=.false.
    if(.not.this%initialized)then
       call Warn('ILDM not initialized.')
       ILDM_check=1
       return
    end if

    if(check(this%H).NE.0)then
       call Warn('ILDM%H failed check.')
       ILDM_check=1
       return
    end if

    if(check(this%den).NE.0)then
       call Warn('ILDM%den failed check.')
       ILDM_check=1
       return
    end if

    !type(quantum),pointer::qs
    if(.not.associated(this%qs))then
       call Warn('ILDM%qs is not associated.')
       ILDM_check=1
       return
    end if
    if(check(this%qs).NE.0)then
       call Warn('ILDM%qs failed check.')
       ILDM_check=1
       return
    end if
    if(this%H%nstate.NE.this%qs%hs%nstate)then
       call Warn('ILDM%H%nstate.NE.ILDM%qs%hs%nstate')
       ILDM_check=1
       return
    end if
    if(this%den%nstate.NE.this%qs%hs%nstate)then
       call Warn('ILDM%den%nstate.NE.ILDM%qs%hs%nstate')
       ILDM_check=1
       return
    end if

    !type(classical),pointer::cs
    if(.not.associated(this%cs))then
       call Warn('ILDM%cs is not associated.')
       ILDM_check=1
       return
    end if
    if(check(this%cs).NE.0)then
       call Warn('ILDM%cs failed check.')
       ILDM_check=1
       return
    end if

    !type(coupling),pointer::cp
    if(.not.associated(this%cp))then
       call Warn('ILDM%cp is not associated.')
       ILDM_check=1
       return
    end if
    if(check(this%cp).NE.0)then
       call Warn('ILDM%cp failed check.')
       ILDM_check=1
       return
    end if

    !real(double)::shortime
    if(this%shorttime.NE.this%shorttime)then
       call Warn('ILDM%shorttime is NAN')
       ILDM_check=1
       return
    end if
    if(this%shorttime.LE.0._double)then
       call Warn('ILDM%shorttime is not a positive value.')
       ILDM_check=1
       return
    end if
    if(this%shorttime.LT.epsilon(this%shorttime))then
       call Warn('ILDM%shorttime is tiny.')
       ILDM_check=1
       return
    end if
    if(this%shorttime.GT.huge(this%shorttime))then
       call Warn('ILDM%shorttime is huge.')
       ILDM_check=1
       return
    end if
    if(this%shorttime.GT.this%den%runtime)then
       call Warn('ILDM%shortime is larger than runtime.')
       ILDM_check=1
       return
    end if
    if(this%shorttime.LT.this%den%dtn)then
       call Warn('ILDM%shortime is less than nuclear time step.')
       ILDM_check=1
       return
    end if

  end function ILDM_check
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine ILDM_run(this,file,tol,normalize)
    implicit none
    type(ILDM),intent(inout)::this
    character*(*),intent(in),optional::file
    real(double),intent(in),optional::tol
    logical,intent(in),optional::normalize

    integer(long)::timestart(8),timedif(8)
    real(double)::elapsedtime

    integer(long)::unit

    real(double),parameter:: x0i=1.0_double,p0i=1.0_double
    real(double),parameter::x0it=1.0_double,p0it=-1.0_double

    real(double),allocatable::fnon(:)
    real(double),allocatable::qsin(:),psin(:)
    real(double),allocatable::qsnew(:),psnew(:)
    real(double),allocatable::qsjudge(:,:,:),psjudge(:,:,:)
    complex(double),allocatable::bigA(:,:),bigAjudge(:,:,:,:,:)

    real(double)::dtn,dtn2,tol0,dtdump
    real(double)::wmc,wnorm,cume,rsq
    real(double)::w,wt,wn,wtn,norm
    complex(double)::aweight
    integer(long)::nstate,ndof,ntraj
    integer(long)::itraj,istep,ndump
    integer(long)::iatt,nstepA
    integer(long)::istate_read,istatet_read
    integer(long)::istate,istatet
    integer(long)::fstate,fstatet
    integer(long)::i,j

    logical::recdyn
    character(len=path)::FMT
    if(check(this).EQ.1)call stop('ILDM_run, failed object check')
    call Note('Begin ILDM_run.')
    if(present(file))call Note('Input file= '//file)
    if(present(tol))call Note('Input tol= '//trim(float2str(tol)))
    FMT="("//trim(int2str(this%qs%hs%nstate**2*2+1))//"(ES18.10E2,1X))"
    timestart=timearray()

    tol0=1E-5
    if(present(tol))tol0=tol
    dtn=this%den%dtn
    dtn2=0.5_double*dtn*dtn
    nstate=this%qs%hs%nstate
    ndof=this%cs%hb%ndof
    ntraj=this%den%ntraj
    dtdump=this%den%dtout
    nstepA=int(this%shorttime/dtdump)
    this%den%normalize=.true.
    if(present(normalize))this%den%normalize=normalize

    if(allocated(bigA))deallocate(bigA)
    allocate(bigA(nstate,nstate))
    if(allocated(bigAjudge))deallocate(bigAjudge)
    allocate(bigAjudge(nstate,nstate,nstate,nstate,nstepA))
    if(allocated(qsjudge))deallocate(qsjudge)
    allocate(qsjudge(ndof,nstate,nstate))
    if(allocated(psjudge))deallocate(psjudge)
    allocate(psjudge(ndof,nstate,nstate))
    if(allocated(qsnew))deallocate(qsnew)
    allocate(qsnew(ndof))
    if(allocated(psnew))deallocate(psnew)
    allocate(psnew(ndof))
    if(allocated(qsin))deallocate(qsin)
    allocate(qsin(ndof))
    if(allocated(psin))deallocate(psin)
    allocate(psin(ndof))
    if(allocated(fnon))deallocate(fnon)
    allocate(fnon(ndof))

    istate_read=1
    istatet_read=1
    call warn('ILDM_run: initial state is hard coded: 1, 1')

    !initialize quantum subsystem to consider only strongly coupled states
    do fstate=1,nstate
       do fstatet=1,nstate
          if(abs(this%qs%hs%diabat(fstate,fstatet)).LT.this%den%tol&
               .and.(fstate.NE.fstatet)) this%qs%hs%diabat(fstate,fstatet)=0.0_double
       end do
    end do
    call update(this%qs)

    !initialize reduced density matrix
    this%den%buf=cmplx(0.,0.)
    bigAjudge=cmplx(0.,0.)
    this%den%bufok=.true.

    !open ECON dyn file
    recdyn=.false.
    if(present(file).and.myid.EQ.master)recdyn=.true.
    if(recdyn)then
       unit=newunit()
       open(unit,file=file//'.dyn')
    end if

    !loop over trajectories
    itraj=0
    do while(itraj.LT.ntraj)

       !sample initial state
       istatet = istatet_read
       istate = istate_read

       wmc=1._double          !>Initialize trajectory MC weight
       aweight=(1._double,0.) !>Initialize trajectory Phase weight

       !init buffer
       this%den%buf(istate,istatet,0,1)=wmc*aweight
       this%den%bufok=.true.

       !Resample system
       call resample(this%cs)
       call resample(this%cp)
       qsin=this%cs%hb%Q
       psin=this%cs%hb%P

       !Hop attempt loop
       do iatt=1,this%nattempt

          !Loop over possible final states
          do fstate=1,nstate
             if((istate.EQ.fstate).or.&
                  abs(this%qs%hs%diabat(istate,fstate))&
                  .GE.this%den%tol)then
                do fstatet=1,nstate
                   if((istatet.EQ.fstatet).or.&
                        abs(this%qs%hs%diabat(istatet,fstatet))&
                        .GE.this%den%tol)then
                      
                      !resample system at each hop
                      if(iatt.eq.1) then
                         this%cs%hb%Q=qsin
                         this%cs%hb%P=psin
                      else
                         this%cs%hb%Q=qsnew
                         this%cs%hb%P=psnew
                      end if
                      
                      !Compute New force 
                      call update(this%cs)
                      call update(this%cp)

                      !resample mapping variables
                      this%H%q=0.
                      this%H%p=0.
                      this%H%qt=0.
                      this%H%pt=0.
                      this%H%q(istate)=x0i
                      this%H%p(istate)=p0i
                      this%H%qt(istatet)=x0it
                      this%H%pt(istatet)=p0it

                      !Calculate nonadiabatic Forces
                      wn=(this%H%q(fstate)**2+this%H%p(fstate)**2)
                      wtn=(this%H%qt(fstatet)**2+this%H%pt(fstatet)**2)
                      fnon=0.0_double
                      do i=1,nstate
                         w=0.0_double
                         wt=0.0_double
                         if(wn.GT.tol0)w=(this%H%q(fstate)*this%H%q(i)&
                              +this%H%p(fstate)*this%H%p(i))/wn
                         if(i.EQ.fstate)w=1.0
                         if(wtn.GT.tol0)wt=(this%H%qt(fstatet)*this%H%qt(i)&
                              +this%H%pt(fstatet)*this%H%pt(i))/wtn
                         if(i.EQ.fstatet)wt=1.0
                         fnon=fnon-this%cp%hc%dV(:,i,fstate)*w&
                              -this%cp%hc%dV(:,i,fstatet)*wt
                         if(fnon(i).NE.fnon(i))then
                            call warn('ILDM_run: NAN nondaiabatic force.')
                            this%den%bufok=.false.
                            goto 32
                         end if
                         if(abs(fnon(i)).GE.huge(fnon(i)))then
                            call warn('ILDM_run: Huge nondaiabatic force.')
                            this%den%bufok=.false.
                            goto 32
                         end if
                      end do
                      fnon=0.5_double*fnon+this%cs%hb%F

                      !Main time stepping loop
                      ndump=0
                      do istep=1,this%den%nbstep 

                         !Advance nuclear positions with current forces
                         !first half Ver(let) step 
                         this%cs%hb%Q=this%cs%hb%Q&
                              +this%cs%hb%rmass*this%cs%hb%P*dtn&
                              +this%cs%hb%rmass*fnon*dtn2
                         this%cs%hb%P=this%cs%hb%P&
                              +0.5*fnon*dtn
                         
                         !Compute New force 
                         call update(this%cs)
                         call update(this%cp)

                         !Advance mapping variables 
                         call mapverlet(this%H,this%den%dte,this%den%nlstep&
                              ,this%qs,this%cs,this%cp)

                         !Calculate nonadiabatic Forces
                         wn=(this%H%q(fstate)**2+this%H%p(fstate)**2)
                         wtn=(this%H%qt(fstatet)**2+this%H%pt(fstatet)**2)
                         fnon=0.0_double
                         do i=1,nstate
                            w=0.0_double
                            wt=0.0_double
                            if(wn.GT.tol0)w=(this%H%q(fstate)*this%H%q(i)&
                                 +this%H%p(fstate)*this%H%p(i))/wn
                            if(i.EQ.fstate)w=1.0
                            if(wtn.GT.tol0)wt=(this%H%qt(fstatet)*this%H%qt(i)&
                                 +this%H%pt(fstatet)*this%H%pt(i))/wtn
                            if(i.EQ.fstatet)wt=1.0
                            fnon=fnon-this%cp%hc%dV(:,i,fstate)*w&
                                 -this%cp%hc%dV(:,i,fstatet)*wt
                            if(fnon(i).NE.fnon(i))then
                               call warn('ILDM_run: NAN nondaiabatic force.')
                               this%den%bufok=.false.
                               goto 32
                            end if
                            if(abs(fnon(i)).GE.huge(fnon(i)))then
                               call warn('ILDM_run: Huge nondaiabatic force.')
                               this%den%bufok=.false.
                               goto 32
                            end if
                         end do
                         fnon=0.5_double*fnon+this%cs%hb%F

                         !second half of (ver)Let
                         this%cs%hb%P=this%cs%hb%P+0.5*fnon*dtn

                         !record output
                         if(istep*dtn.GE.(ndump+1)*dtdump)then
                            ndump=ndump+1

                            !Construct bigA
                            do i=1,nstate
                               do j=1,nstate
                                  bigA(i,j)=0.25*(this%H%q(i)+eye*this%H%p(i))&
                                       *(this%H%qt(j)-eye*this%H%pt(j))&
                                       *(x0i-eye*p0i)*(x0it+eye*p0it)
                               end do
                            end do
                            bigAjudge(:,:,fstate,fstatet,ndump)=bigA

                            !optional dyn record
                            if(fstate.EQ.1.and.fstatet.EQ.1&
                                 .and.itraj.LE.1.and.recdyn)then
                               call update(this%cs)
                               call update(this%cp)
                               write(unit,"(10(ES18.10E2,1X))")&
                                    istep*dtn&
                                    ,this%qs%hs%diabat(fstate,fstatet)&
                                    ,this%cs%hb%V,this%cs%hb%T&
                                    ,this%cp%hc%V(fstate,fstatet)&
                                    ,sum(this%H%q**2+this%H%p**2)&
                                    ,sum(this%H%qt**2+this%H%pt**2)&
                                    ,this%cs%hb%Q(1),this%cs%hb%P(1)
                            end if!recdyn loop
                         end if !record loop
                         
                      end do!End of Nuclear Stepping loop
                      if(ndump.NE.nstepA)&
                           call stop('ILDM_run: mismatch in output count'&
                           ,'make sure dtout is an integer number'&
                           //' of nuclear timesteps.')
                      qsjudge(:,fstate,fstatet)=this%cs%hb%Q
                      psjudge(:,fstate,fstatet)=this%cs%hb%P

                   end if!end istatet-fstatet coupling condition
                end do!end fstatet loop
             end if!end istate-fstate coupling condition
          end do!end fstate loop

          !close dyn file
          if(recdyn.and.itraj.GT.1)close(unit)

          !Monte Carlo branching
          wnorm=0.
          do i=1,nstate
             do j=1,nstate
                wnorm=wnorm+sqrt(bigAjudge(i,j,i,j,nstepA)&
                     *conjg(bigAjudge(i,j,i,j,nstepA)))
             end do
          end do
          w=ran0()*wnorm
          cume=0.
          do i=1,nstate
             do j=1,nstate
                rsq=sqrt(bigAjudge(i,j,i,j,nstepA)&
                     *conjg(bigAjudge(i,j,i,j,nstepA)))
                cume=cume+rsq
                if(w.le.cume) then
                   istate=i
                   istatet=j
                   goto 10
                end if
             end do
          end do
10        continue

          !save chosen initial state phase-space coords
          qsnew=qsjudge(:,istate,istatet)
          psnew=psjudge(:,istate,istatet)

          !reconstruct redmat up to current hopping time for recording
          do istep=1,nstepA-1
             do i=1,nstate
                do j=1,nstate
                   !this%den%buf(i,j,(iatt-1)*nstepA+istep,1)&
                   !     =this%den%buf(i,j,(iatt-1)*nstepA+istep,1)&
                   !     +bigAjudge(i,j,istate,istatet,istep)*wmc*aweight
                   this%den%buf(i,j,(iatt-1)*nstepA+istep,1)&
                        =bigAjudge(i,j,istate,istatet,istep)*wmc*aweight
                end do
             end do
          end do
          
          !new trajectory MC weight
          wmc=wmc*wnorm

          !new trajectory phase weight
          bigA=bigAjudge(:,:,istate,istatet,nstepA)
          aweight = aweight*bigA(istate,istatet)&
               /sqrt(bigA(istate,istatet)&
               *conjg(bigA(istate,istatet)))

          !accumulate redmat at next hop
          this%den%buf(istate,istatet,iatt*nstepA,1)&
               =wmc*aweight
               !=this%den%buf(istate,istatet,iatt*nstepA,1)&
               !+wmc*aweight
          
       end do!end nattempt loop

32        continue
          !save buf if traj passed 
          if(this%den%bufok)then
             this%den%buf(:,:,:,0)=this%den%buf(:,:,:,0)+this%den%buf(:,:,:,1)
             itraj=itraj+1
          end if

       if(present(file))then
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
          
          if (elapsedtime.GE.this%den%dtsave)then
             call Prompt('Saving intermediate density ntraj= '&
                  //trim(int2str(itraj)))
             norm=1._double
             if(this%den%normalize)norm=real(itraj)

             unit=newunit()
             open(unit,file=file)
             do istep=0,this%nattempt*nstepA!this%den%nstep
                write(unit,trim(FMT))&
                     istep*dtdump,((real(this%den%buf(i,j,istep,0))/norm&
                     ,aimag(this%den%buf(i,j,istep,0))/norm&
                     ,j=1,nstate),i=1,nstate)
             end do
             close(unit)
             timestart=timearray()
             this%qs%hs%den=this%den%buf(:,:,this%den%nstep,0)/real(itraj)

          end if
       end if
    end do!end ntraj loop


    !write(*,*)'***',nstepA,this%den%nstep,this%nattempt*nstepA
    !write(*,*)'***',this%den%runtime,this%den%dtout,dtdump
    !write(*,*)'***',ndump
                      
!stop


    if(present(file))then
       norm=1._double
       if(this%den%normalize)norm=real(itraj)
       
       unit=newunit()
       open(unit,file=file)
       do istep=0,this%nattempt*nstepA!this%den%nstep
          write(unit,trim(FMT))istep*dtdump,&
               ((real(this%den%buf(i,j,istep,0))/norm&
               &,aimag(this%den%buf(i,j,istep,0))/norm&
               ,j=1,nstate),i=1,nstate)
       end do
       close(unit)
    end if
    this%qs%hs%den=this%den%buf(:,:,this%den%nstep,0)/real(itraj)

    if(allocated(qsin))deallocate(qsin)
    if(allocated(psin))deallocate(psin)
    if(allocated(fnon))deallocate(fnon)

  end subroutine ILDM_run
!---------------------------------------------------
!!$ subroutine LDM_new(this,qs,cs,cp,file)
!!$    implicit none
!!$    type(LDM),intent(inout)::this
!!$    type(quantum),intent(in),target::qs
!!$    type(classical),intent(in),target::cs
!!$    type(coupling),intent(in),target::cp
!!$    character*(*),intent(in),optional::file
!!$    character(len=title)::filetype
!!$    character(len=path)::filename
!!$    integer(long)::unit
!!$
!!$    call Note('Begin LDM_new.')
!!$    if(present(file))call Note('Input file= '//file)
!!$    if(check(qs).EQ.1)
!!$    if(check(cs).EQ.1)
!!$    if(check(cp).EQ.1)
!!$
!!$    if(associated(this%qs))nullify(this%qs)
!!$    this%qs=>qs
!!$    if(associated(this%cs))nullify(this%cs)
!!$    this%cs=>cs
!!$    if(associated(this%cp))nullify(this%cp)
!!$    this%cp=>cp
!!$
!!$    call new(this%H,this%qs%hs%nstate)
!!$
!!$    if(present(file))then
!!$       unit=newunit()
!!$       open(unit,file=file)
!!$       read(unit,*)filetype
!!$       filetype=adjustl(filetype)
!!$       if(trim(filetype).NE.'LDM')&
!!$            call stop('LDM_init error: not a valid input file.')
!!$       read(unit,*)filename
!!$       filename=adjustl(filename)
!!$       call new(this%den,nstate=this%qs%hs%nstate,file=trim(filename))
!!$    else
!!$       call new(this%den,nstate=this%qs%hs%nstate)
!!$    end if
!!$
!!$    if(present(file))close(unit)
!!$    this%initialized=.true.
!!$  end subroutine LDM_new
!!$  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!!$  subroutine LDM_kill(this)
!!$    implicit none
!!$    type(LDM),intent(inout)::this
!!$    call Note('Begin LDM_kill.')
!!$    this%initialized=.false.
!!$    call kill(this%H)
!!$    call kill(this%den)
!!$    if(associated(this%qs))nullify(this%qs)
!!$    if(associated(this%cs))nullify(this%cs)
!!$    if(associated(this%cp))nullify(this%cp)
!!$  end subroutine LDM_kill
!!$  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!!$  subroutine LDM_display(this)
!!$    implicit none
!!$    type(LDM),intent(in)::this
!!$    call Note('Begin LDM_display.')
!!$    if(check(this).NE.0)then
!!$       call warn('display_LDM: object failed check.','displaying nothing.')
!!$    else
!!$       call display(this%H)
!!$       call display(this%den)
!!$       call display(this%qs)
!!$       call display(this%cs)
!!$       call display(this%cp)
!!$    end if
!!$  end subroutine LDM_display
!!$  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!!$  subroutine LDM_save(this,file)
!!$    implicit none
!!$    type(LDM),intent(in)::this
!!$    character*(*),intent(in)::file
!!$
!!$    integer(long)::unit    
!!$
!!$    call Note('Begin LDM_save.')
!!$    call Note('Input file= '//file)
!!$    if(check(this).EQ.1)
!!$
!!$    unit=newunit()
!!$    open(unit,file=file)
!!$    write(unit,*)'LDM'
!!$    write(unit,*)quote(file//'.den')
!!$    call save(this%den,file=file//'.den')
!!$    write(unit,*)0._double
!!$    close(unit)  
!!$
!!$  end subroutine LDM_save
!!$  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!!$  integer(short) function LDM_check(this)
!!$    type(LDM),intent(in)::this
!!$    LDM_check=0
!!$
!!$    call Note('Checking LDM.')
!!$
!!$    !logical::initialized=.false.
!!$    if(.not.this%initialized)then
!!$       call Warn('LDM not initialized.')
!!$       LDM_check=1
!!$       return
!!$    end if
!!$
!!$    if(check(this%H).NE.0)then
!!$       call Warn('LDM%H failed check.')
!!$       LDM_check=1
!!$       return
!!$    end if
!!$
!!$    if(check(this%den).NE.0)then
!!$       call Warn('LDM%den failed check.')
!!$       LDM_check=1
!!$       return
!!$    end if
!!$
!!$    !type(quantum),pointer::qs
!!$    if(.not.associated(this%qs))then
!!$       call Warn('LDM%qs is not associated.')
!!$       LDM_check=1
!!$       return
!!$    end if
!!$    if(check(this%qs).NE.0)then
!!$       call Warn('LDM%qs failed check.')
!!$       LDM_check=1
!!$       return
!!$    end if
!!$    if(this%H%nstate.NE.this%qs%hs%nstate)then
!!$       call Warn('LDM%H%nstate.NE.LDM%qs%hs%nstate')
!!$       LDM_check=1
!!$       return
!!$    end if
!!$    if(this%den%nstate.NE.this%qs%hs%nstate)then
!!$       call Warn('LDM%den%nstate.NE.LDM%qs%hs%nstate')
!!$       LDM_check=1
!!$       return
!!$    end if
!!$
!!$    !type(classical),pointer::cs
!!$    if(.not.associated(this%cs))then
!!$       call Warn('LDM%cs is not associated.')
!!$       LDM_check=1
!!$       return
!!$    end if
!!$    if(check(this%cs).NE.0)then
!!$       call Warn('LDM%cs failed check.')
!!$       LDM_check=1
!!$       return
!!$    end if
!!$    !if(.not.this%cs%hb%wignerbath)then
!!$    !   call Warn('LDM%cs%hb%wignerbath not true.')
!!$    !   LDM_check=1
!!$    !   return
!!$    !end if
!!$
!!$    !type(coupling),pointer::cp
!!$    if(.not.associated(this%cp))then
!!$       call Warn('LDM%cp is not associated.')
!!$       LDM_check=1
!!$       return
!!$    end if
!!$    if(check(this%cp).NE.0)then
!!$       call Warn('LDM%cp failed check.')
!!$       LDM_check=1
!!$       return
!!$    end if
!!$
!!$  end function LDM_check
!!$  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!!$  subroutine LDM_run(this,file,tol,normalize)
!!$    implicit none
!!$    type(LDM),intent(inout)::this
!!$    character*(*),intent(in),optional::file
!!$    real(double),intent(in),optional::tol
!!$    logical,intent(in),optional::normalize
!!$
!!$    integer(long)::timestart(8),timedif(8)
!!$    real(double)::elapsedtime
!!$
!!$    integer(long)::unit
!!$
!!$    real(double),parameter:: x0i=1.0_double,p0i=1.0_double
!!$    real(double),parameter::x0it=1.0_double,p0it=-1.0_double
!!$
!!$    real(double),allocatable::qsin(:),psin(:)
!!$    real(double),allocatable::fnon(:)
!!$
!!$    real(double)::dtn,dtn2,tol0,dtdump
!!$    real(double)::w,wn,wt,wtn,norm
!!$    integer(long)::nstate,ndof,ntraj
!!$    integer(long)::itraj,istep,ndump
!!$    integer(long)::istate,istatet
!!$    integer(long)::fstate,fstatet
!!$    integer(long)::i,j
!!$
!!$    logical::recdyn
!!$    character(len=path)::FMT
!!$
!!$    if(check(this).EQ.1)
!!$    call Note('Begin LDM_run.')
!!$    Call PROMPT('LDM_run is incomplete. USE ILDM type with short time approx = total runtime for equivalent dynamics.')
!!$    call warn('LDM_run is incomplete.','use ILDM type with short time approx = total runtime for equivalent dynamics.')
!!$    if(present(file))call Note('Input file= '//file)
!!$    if(present(tol))call Note('Input tol= '//trim(float2str(tol)))
!!$    timestart=timearray()
!!$    FMT="("//trim(int2str(this%qs%hs%nstate**2*2+1))//"(ES18.10E2,1X))"
!!$
!!$    tol0=1E-5
!!$    if(present(tol))tol0=tol
!!$    dtn=this%den%dtn
!!$    dtn2=0.5_double*dtn*dtn
!!$    nstate=this%qs%hs%nstate
!!$    ndof=this%cs%hb%ndof
!!$    ntraj=this%den%ntraj
!!$    dtdump=this%den%dtout
!!$    this%den%normalize=.true.
!!$    if(present(normalize))this%den%normalize=normalize
!!$
!!$    if(allocated(qsin))deallocate(qsin)
!!$    allocate(qsin(ndof))
!!$    if(allocated(psin))deallocate(psin)
!!$    allocate(psin(ndof))
!!$    if(allocated(fnon))deallocate(fnon)
!!$    allocate(fnon(ndof))
!!$
!!$    istate=1
!!$    istatet=1
!!$    call warn('LDM_run: initial state is hard coded: 1, 1')
!!$
!!$    !initialize quantum subsystem
!!$    do fstate=1,nstate
!!$       do fstatet=1,nstate
!!$          if(abs(this%qs%hs%diabat(fstate,fstatet)).LT.this%den%tol&
!!$               .and.(fstate.ne.fstatet))&
!!$               this%qs%hs%diabat(fstate,fstatet)=0.0_double
!!$       end do
!!$    end do
!!$    call update(this%qs)
!!$
!!$    !initialize reduced density matrix
!!$    this%den%buf=cmplx(0.,0.)
!!$    this%den%bufok=.true.
!!$
!!$    !open ECON dyn file
!!$    recdyn=.false.
!!$    if(present(file).and.(myid.EQ.master))recdyn=.true.
!!$    if(recdyn)then
!!$       unit=newunit()
!!$       open(unit,file=file//'.dyn')
!!$    end if
!!$
!!$    !loop over trajectories
!!$    itraj=0
!!$    do while (itraj.LT.ntraj)
!!$
!!$       !current buf = last good buf
!!$       this%den%buf(:,:,:,1)=this%den%buf(:,:,:,0)
!!$       this%den%bufok=.true.
!!$
!!$       !Resample classical subsystem
!!$       call resample(this%cs)
!!$       call resample(this%cp)
!!$
!!$       !initialize mapping variables before time loop
!!$       qsin=this%cs%hb%Q
!!$       psin=this%cs%hb%P
!!$
!!$       !Loop over possible final states
!!$       do fstate=1,nstate
!!$          if((istate.EQ.fstate).or.&
!!$               abs(this%qs%hs%diabat(istate,fstate))&
!!$               .GE.this%den%tol)then
!!$             do fstatet=1,nstate
!!$                if((istatet.EQ.fstatet).or.&
!!$                     abs(this%qs%hs%diabat(istatet,fstatet))&
!!$                     .GE.this%den%tol)then
!!$                   
!!$                   !consider only strongly coupled states
!!$                   this%H%q=0.
!!$                   this%H%p=0.
!!$                   this%H%qt=0.
!!$                   this%H%pt=0.
!!$                   
!!$                   this%H%q(istate)=x0i
!!$                   this%H%p(istate)=p0i
!!$                   this%H%qt(istatet)=x0it
!!$                   this%H%pt(istatet)=p0it
!!$                   
!!$                   this%cs%hb%Q=qsin
!!$                   this%cs%hb%P=psin
!!$                   
!!$                   !Compute New force 
!!$                   call update(this%cs)
!!$                   call update(this%cp)
!!$                   
!!$                   !Calculate nonadiabatic Forces
!!$                   wn=(this%H%q(fstate)**2+this%H%p(fstate)**2)
!!$                   wtn=(this%H%qt(fstatet)**2+this%H%pt(fstatet)**2)
!!$                   fnon=0.0_double
!!$                   do i=1,nstate
!!$                      w=0.0_double
!!$                      wt=0.0_double
!!$                      if(wn.GT.tol0)w=(this%H%q(fstate)*this%H%q(i)&
!!$                           +this%H%p(fstate)*this%H%p(i))/wn
!!$                      if(i.EQ.fstate)w=1.0
!!$                      if(wtn.GT.tol0)wt=(this%H%qt(fstatet)*this%H%qt(i)&
!!$                           +this%H%pt(fstatet)*this%H%pt(i))/wtn
!!$                      if(i.EQ.fstatet)wt=1.0
!!$                      fnon=fnon-this%cp%hc%dV(:,i,fstate)*w&
!!$                           -this%cp%hc%dV(:,i,fstatet)*wt
!!$                         if(fnon(i).NE.fnon(i))then
!!$                            call warn('LDM_run: NAN nondaiabatic force.'&
!!$                                 ,'thowing away trajectory.')
!!$                            this%den%bufok=.false.
!!$                            goto 30
!!$                         end if
!!$                         if(abs(fnon(i)).GE.huge(fnon(i)))then
!!$                            call warn('LDM_run: Huge nondaiabatic force.'&
!!$                                 ,'thowing away trajectory.')
!!$                            this%den%bufok=.false.
!!$                            goto 30
!!$                         end if
!!$                   end do
!!$                   fnon=0.5_double*fnon+this%cs%hb%F
!!$                   
!!$                   !Main time stepping loop
!!$                   ndump=0
!!$                   do istep=1,this%den%nbstep 
!!$                      !Advance nuclear positions with current forces
!!$                      !first half Ver(let) step 
!!$                      this%cs%hb%Q=this%cs%hb%Q&
!!$                           +this%cs%hb%rmass*this%cs%hb%P*dtn&
!!$                           +this%cs%hb%rmass*fnon*dtn2
!!$                      this%cs%hb%P=this%cs%hb%P&
!!$                           +0.5*fnon*dtn
!!$                      
!!$                      !Compute New force 
!!$                      call update(this%cs)
!!$                      call update(this%cp)
!!$                      
!!$                      !Advance mapping variables 
!!$                      call mapverlet(this%H,this%den%dte,this%den%nlstep&
!!$                           ,this%qs,this%cs,this%cp)
!!$                      
!!$                      !Calculate nonadiabatic Forces
!!$                      wn=(this%H%q(fstate)**2+this%H%p(fstate)**2)
!!$                      wtn=(this%H%qt(fstatet)**2+this%H%pt(fstatet)**2)
!!$                      fnon=0.0_double
!!$                      do i=1,nstate
!!$                         w=0.0_double
!!$                         wt=0.0_double
!!$                         if(wn.GT.tol0)w=(this%H%q(fstate)*this%H%q(i)&
!!$                              +this%H%p(fstate)*this%H%p(i))/wn
!!$                         if(i.EQ.fstate)w=1.0
!!$                         if(wtn.GT.tol0)wt=(this%H%qt(fstatet)*this%H%qt(i)&
!!$                              +this%H%pt(fstatet)*this%H%pt(i))/wtn
!!$                         if(i.EQ.fstatet)wt=1.0
!!$                         fnon=fnon-this%cp%hc%dV(:,i,fstate)*w&
!!$                              -this%cp%hc%dV(:,i,fstatet)*wt
!!$                         if(fnon(i).NE.fnon(i))then
!!$                            call warn('LDM_run: NAN nondaiabatic force.'&
!!$                                 ,'thowing away trajectory.')
!!$                            this%den%bufok=.false.
!!$                            goto 30
!!$                         end if
!!$                         if(abs(fnon(i)).GE.huge(fnon(i)))then
!!$                            call warn('LDM_run: Huge nondaiabatic force.'&
!!$                                 ,'thowing away trajectory.')
!!$                            this%den%bufok=.false.
!!$                            goto 30
!!$                         end if
!!$                      end do
!!$                      fnon=0.5_double*fnon+this%cs%hb%F
!!$                      
!!$                      !second half of (ver)Let
!!$                      this%cs%hb%P=this%cs%hb%P+0.5*fnon*dtn
!!$                      
!!$                      !record output
!!$                      if(istep*dtn.GE.(ndump+1)*dtdump)then
!!$                         ndump=ndump+1
!!$
!!$                         !Construct buf
!!$                         do i=1,nstate
!!$                            do j=1,nstate
!!$                               this%den%buf(i,j,ndump,1)=&
!!$                                    this%den%buf(i,j,ndump,1)&
!!$                                    +0.25*(this%H%q(i)+eye*this%H%p(i))&
!!$                                    *(this%H%qt(j)-eye*this%H%pt(j))&
!!$                                    *(x0i-eye*p0i)*(x0it+eye*p0it)      
!!$                            end do
!!$                         end do
!!$                              
!!$                                                 
!!$                         if(fstate.EQ.1.and.fstatet.EQ.1.and.itraj.LE.1)then
!!$                            if(recdyn)then
!!$                               call update(this%cs)
!!$                               call update(this%cp)
!!$                               !call prompt('dyn itraj'//trim(int2str(itraj)))
!!$                               write(unit,"(20(ES18.10E2,1X))")&
!!$                                    istep*dtn,this%qs%hs%diabat(fstate,fstatet)&
!!$                                    ,this%cs%hb%V,this%cs%hb%T&
!!$                                    ,this%cp%hc%V(fstate,fstatet)&
!!$                                    ,sum(this%H%q**2+this%H%p**2)&
!!$                                    ,sum(this%H%qt**2+this%H%pt**2)&
!!$                                    ,this%cs%hb%Q(1),this%cs%hb%P(1)&
!!$                                    ,this%H%q(fstate),this%H%p(fstate)&
!!$                                    ,this%H%qt(fstatet),this%H%pt(fstatet)&
!!$                                    ,x0i,p0i,x0it,p0it
!!$                            end if
!!$                         end if
!!$                      end if
!!$                   end do!End of Nuclear Stepping loop
!!$                   if(ndump.NE.this%den%nstep)&
!!$                        call stop('LDM_run: mismatch in output count'&
!!$                        ,'make sure dtout is an integer number'&
!!$                        //' of nuclear timesteps.')
!!$                end if!end istate-fstate coupling condition
!!$             end do!end fstatet loop
!!$          end if!end of istatet-fstatet coupling condition
!!$       end do!end fstate loop
!!$
!!$       !save buf if traj passed 
!!$30     if(this%den%bufok)then
!!$          this%den%buf(:,:,:,0)=this%den%buf(:,:,:,1)
!!$          itraj=itraj+1
!!$       end if
!!$
!!$       if(recdyn.and.itraj.GT.1)close(unit) !close dyn file
!!$
!!$       if(present(file))then
!!$          timedif=timearray()-timestart
!!$          if(timedif(8).LT.0)then
!!$             timedif(7)=timedif(7)-1
!!$             timedif(8)=timedif(8)+1E3
!!$          end if
!!$          if(timedif(7).LT.0)then
!!$             timedif(6)=timedif(6)-1
!!$             timedif(7)=timedif(7)+60
!!$          end if
!!$          if(timedif(6).LT.0)then
!!$             timedif(5)=timedif(5)-1
!!$             timedif(6)=timedif(6)+60
!!$          end if
!!$          if(timedif(5).LT.0)then
!!$             timedif(4)=timedif(4)-1
!!$             timedif(5)=timedif(5)+24
!!$          end if
!!$          elapsedtime=timedif(8)*1E-3
!!$          elapsedtime=elapsedtime+timedif(7)
!!$          elapsedtime=elapsedtime+timedif(6)*60
!!$          elapsedtime=elapsedtime+timedif(5)*1440
!!$
!!$          if (elapsedtime.GE.this%den%dtsave)then
!!$             call Prompt('Saving intermediate density ntraj= '&
!!$                  //trim(int2str(itraj)))
!!$             norm=1._double
!!$             if(this%den%normalize)norm=real(itraj)
!!$
!!$             unit=newunit()
!!$             open(unit,file=file)
!!$             do istep=1,this%den%nstep
!!$                write(unit,trim(FMT))&
!!$                     istep*dtdump,((real(this%den%buf(i,j,istep,0))/norm&
!!$                     ,aimag(this%den%buf(i,j,istep,0))/norm&
!!$                     ,j=1,nstate),i=1,nstate)
!!$             end do
!!$             close(unit)
!!$             timestart=timearray()
!!$             this%qs%hs%den=this%den%buf(:,:,this%den%nstep,0)/real(itraj)
!!$          end if
!!$       end if
!!$    end do!end ntraj loop
!!$
!!$    if(present(file))then
!!$       norm=1._double
!!$       if(this%den%normalize)norm=real(itraj)
!!$
!!$       unit=newunit()
!!$       open(unit,file=file)
!!$       do istep=1,this%den%nstep
!!$          write(unit,trim(FMT))istep*dtdump,&
!!$               ((real(this%den%buf(i,j,istep,0))/norm&
!!$               &,aimag(this%den%buf(i,j,istep,0))/norm&
!!$               ,j=1,nstate),i=1,nstate)
!!$       end do
!!$       close(unit)
!!$    end if
!!$    this%qs%hs%den=this%den%buf(:,:,this%den%nstep,0)/real(ntraj)
!!$
!!$    if(allocated(qsin))deallocate(qsin)
!!$    if(allocated(psin))deallocate(psin)
!!$    if(allocated(fnon))deallocate(fnon)
!!$
!!$  end subroutine LDM_run
!!$  !------------------------------------------------------------------------


end module LDM_class
