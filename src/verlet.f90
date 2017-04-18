module verlet_class
  use atomicunits
  use filemanager
  use Errorlog
  use quantum_class
  use classical_class
  use coupling_class
  implicit none
  Private

  Public::verlet
  Public::new,save,run

  type verlet
     logical::initialized=.false.
     real*8::runtime
     real*8::dtn
     integer::nstep
     type(quantum),pointer::qs
     type(classical),pointer::cs
     type(coupling),pointer::cp
  end type verlet

  interface new
     module procedure new_verlet
  end interface

  interface save
     module procedure save_verlet
  end interface

  interface run
     module procedure run_verlet
  end interface

contains
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine new_verlet(this,qs,cs,cp,file)
    implicit none
    type(verlet),intent(inout)::this
    type(quantum),intent(in),target::qs
    type(classical),intent(in),target::cs
    type(coupling),intent(in),target::cp
    character*(*),intent(in),optional::file

    integer::ierr,unit,nstate
    logical::usedefault,usedunit
    character(len=50)::filetype      

    if(present(file))then
       inquire(file=file,exist=usedefault)
       if(.not.usedefault)then
          write(*,*)'new_verlet error cannot find inputfile '//file
          stop
       end if
    end if

    unit=1000
    usedunit=.true.
    do while (usedunit)
       unit=unit+1
       inquire(unit,opened=usedunit)
    end do

    if(present(file))then
       open(unit,file=file)
       read(unit,'(A50)')filetype
       filetype=adjustl(filetype)
       if(trim(filetype).NE.'verlet')then
          write(*,*)'verlet_init error: not a valid input file.'
          stop
       end if
       read(unit,*)this%runtime
       read(unit,*)this%dtn
    else
       this%runtime=0.0
       do while(this%runtime.LE.0.0)
          write(*,*)'Enter total runtime.'
          read(*,*)this%runtime
          if(this%runtime.LE.0.0)then
             write(*,*)'runtime must be positive. try again!'
          end if
       end do
       write(*,*)'Enter nuclear timestep.'
       read(*,*)this%dtn
       if(this%dtn.GT.this%runtime)then
          call warn('new_verlet: Nuclear time step greater than runtime','seting runtime to equal nuclear time step.')
          this%runtime=this%dtn
       end if
    end if
    this%nstep=int(this%runtime/this%dtn)

    if(.not.qs%initialized)then
       write(*,*)'verlet_init error: quantum object not initialized.'
       stop
    end if

    if(associated(this%qs))nullify(this%qs)
    this%qs=>qs

    if(.not.cs%initialized)then
       write(*,*)'verlet_init error: classical object not initialized.'
       stop
    end if
    if(associated(this%cs))nullify(this%cs)
    this%cs=>cs

    if(.not.cp%initialized)then
       write(*,*)'verlet_init error: coupling object not initialized.'
       stop
    end if
    if(associated(this%cp))nullify(this%cp)
    this%cp=>cp

    if(present(file))close(unit)

    this%initialized=.true.

  end subroutine new_verlet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine save_verlet(this,file)
    implicit none
    type(verlet),intent(in)::this
    character*(*),intent(in)::file

    integer::unit
    logical::usedunit    

    if(.not.this%initialized)then
       write(*,*)'verlet_save error: dynamics not initialized.'
       stop
    end if

    unit=newunit()
    open(unit,file=file)
    write(unit,'(A50)')'verlet'
    write(unit,*)this%runtime
    write(unit,*)this%dtn
    close(unit)  

  end subroutine save_verlet
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine run_verlet(this,file)
    implicit none
    type(verlet),intent(inout)::this
    character*(*),intent(in),optional::file

    integer::unit
    logical::recording

    real*8:: dtn,dtn2,fnon(this%cs%hb%ndof)
    integer::ndof,nstep
    integer::idof,istep
    
    dtn=this%dtn
    dtn2=0.5*dtn*dtn
    ndof=this%cs%hb%ndof
    nstep=this%nstep

    recording=.false.
    if(present(file))then
       unit=newunit()
       open(unit,file=file)
       recording=.true.
    end if

    !initiate system
    istep=0
    call get_GroundStateAdiabaticForce(this,fnon)
    if(recording)then
       write(unit,"(1000(ES18.10E2,1X))")istep*dtn,this%cs%hb%V,this%cs%hb%T
    end if

    !Main time stepping loop
    do istep=1,nstep
       
       !Advance nuclear positions with current forces
       !First half Verlet step 
       this%cs%hb%Q=this%cs%hb%Q&
            +this%cs%hb%rmass*this%cs%hb%P*dtn&
            +this%cs%hb%rmass*fnon*dtn2
       this%cs%hb%P=this%cs%hb%P+0.5*fnon*dtn

       !update forces
       call get_GroundStateAdiabaticForce(this,fnon)
       !second half of (ver)Let
       this%cs%hb%P=this%cs%hb%P+0.5*fnon*dtn

       if(recording)then
          call update(this%cs)
          write(unit,"(1000(ES18.10E2,1X))")istep*dtn,this%cs%hb%V,this%cs%hb%T
       end if

    end do!End of Nuclear Stepping loop

    if(recording)close(unit)
    
  end subroutine run_verlet

  !------------------------------------------------------------------------
  subroutine get_GroundStateAdiabaticForce(this,force)
    implicit none
    type(verlet),intent(inout)::this
    real(double),intent(inout)::force(:)
    complex(double)::H(this%qs%hs%nstate,this%qs%hs%nstate)
    real*8:: Eng(this%qs%hs%nstate)
    integer(long)::nstate,ndof
    integer(long)::istate,jstate,idof

    nstate=this%qs%hs%nstate
    ndof=this%cs%hb%ndof

    call update(this%qs)
    call update(this%cs)
    call update(this%cp)

    !comupte eigen vectors
    H=this%qs%hs%diabat+this%cp%hc%V
    call diagonalize(nstate,H,Eng)

    !compute coupling term dV/dQ in adiabatic representation
    force=0._double
    do idof=1,ndof
       do istate=1,nstate
          do jstate=1,nstate
             force(idof)=force(idof)-real(H(istate,1)*conjg(H(jstate,1)))&
                  *this%cp%hc%dV(idof,istate,jstate)
          end do
       end do
    end do

    !add bath and coupling term Force
    !force=this%cs%hb%F
    force=force+this%cs%hb%F

    !checks
    if(any(force.NE.force))then
       write(*,*)'Warning NAN nondaiabatic force.'
       write(*,*)'Program will stop.'
       stop
    end if
    if(any(force.GT.Huge(force)))then
       write(*,*)'Warning Huge nondaiabatic force.'
       write(*,*)'Program will stop.'
       stop
    end if

  end subroutine get_GroundStateAdiabaticForce

end module verlet_class
