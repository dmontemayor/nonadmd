!============================================================================
!> \brief
!! Classical Subsystem Primitive
!! \details
!! This class defines the classical subsystem kernal inherited by all derived
!! classical subsystems.
!! \date
!! 12 Nov 2012
!! \authors
!!  Daniel Montemayor
!<============================================================================
module hb_class
  use type_kinds
  use ErrorLog
  use filemanager
  use string
  use atomicunits
  use textgraphs
  use outputdisplay
  implicit none
  private

  public::hb
  public::new,kill,update,resample,display,save,check

  type hb
     logical::initialized
     integer(long)::npart
     real(double),dimension(:),pointer::mass,rmass
     real(double),dimension(:),pointer::charge
     integer(byte)::ndim
     integer(long)::ndof
     real(double),dimension(:),pointer::Qmin,Qmax
     logical,dimension(:),pointer::PBC
     real(double),dimension(:),pointer::Q,P,F
     real(double)::V,T
     real(double)::temperature
     real(double),dimension(:),pointer::mode,W
     real(double),dimension(:,:),pointer::EigenVec
  end type hb

  interface new
     module procedure hb_init
  end interface

  interface kill
     module procedure hb_kill
  end interface

  interface update
     module procedure hb_update
  end interface

  interface resample
     module procedure hb_resample
  end interface

  interface display
     module procedure hb_display
  end interface

  interface save
     module procedure hb_save
  end interface

  interface check
     module procedure hb_check
  end interface

contains
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hb_init(this,file)
    type(hb),intent(inout)::this
    character*(*),intent(in),optional::file
    
    character(len=title)::filetype

    integer(long)::npart,ipart,ndim,idim,ndof,idof

    integer(short)::ierr
    integer(long)::unit
    logical::usedefault

    call Note('Begin hb_init.')
    if(present(file))call Note('input file= '//file)

    if(present(file))then

       if(check(file=file).NE.0)call Stop('hb_init: cannot find input file.')

       unit=newunit()
       open(unit,file=file)
       read(unit,*)filetype
       filetype=adjustl(filetype)
       if(trim(filetype).NE.'hb')call stop('hb_init: input file not valid.')
    end if
    
    if(present(file))then
       read(unit,*)npart
    else
       npart=0
       do while (npart.LE.0)
          write(*,*)'Please enter the number of bath particles.'
          read(*,*)npart
          if(npart.LT.0)write(*,*)'Number of bath particles must be positive integer. Try again.'
       end do
    end if
    this%npart=npart

    if(associated(this%mass))nullify(this%mass)
    allocate(this%mass(npart),stat=ierr)
    if(ierr.NE.0)call stop('hb_init: failed mass memory allocation.')
    this%mass=1.0_double
    if(present(file))then
       read(unit,*)(this%mass(ipart),ipart=1,npart)
    else
       write(*,*)'Use default particle mass settings [all unit masses]? (T/F)'
       read(*,*)usedefault
       if(.not.usedefault)then
          do ipart=1,npart
             write(*,*)'Enter particle',ipart,' mass.'
             read(*,*)this%mass(ipart)
          end do
       end if
    end if


    if(associated(this%charge))nullify(this%charge)
    allocate(this%charge(npart),stat=ierr)
    if(ierr.NE.0)call stop('hb_init: failed particle charge memory allocation.')
    this%charge=0._double
    if(present(file))then
       read(unit,*)(this%charge(ipart),ipart=1,npart)
    else
       write(*,*)'Use default particle charge settings [no charged particles]? (T/F)'
       read(*,*)usedefault
       if(.not.usedefault)then
          do ipart=1,npart
             write(*,*)'Enter particle',ipart,' charge.'
             read(*,*)this%charge(ipart)
          end do
       end if
    end if


    if(present(file))then
       read(unit,*)ndim
    else
       ndim=0
       do while (ndim.LE.0)
          write(*,*)'Please enter the number of phase-space dimensions.'
          read(*,*)ndim
          if(ndim.LE.0)write(*,*)'Number of pase-space dimensions must be positve integer. Try again.'
       end do
    end if
    this%ndim=ndim
    ndof=ndim*npart
    this%ndof=ndof

    if(associated(this%Qmin))nullify(this%Qmin)
    if(associated(this%Qmax))nullify(this%Qmax)
    if(associated(this%PBC))nullify(this%PBC)
    allocate(this%Qmin(ndim),this%Qmax(ndim),this%PBC(ndim),stat=ierr)
    if(ierr.NE.0)call stop('hb_init: failed phase-space domain memory allocation.')
    this%Qmin=-1.0_double
    this%Qmax=1.0_double
    this%PBC=.false.

    if(present(file))then
       read(unit,*)(this%Qmin(idim),idim=1,ndim)
       read(unit,*)(this%Qmax(idim),idim=1,ndim)
       read(unit,*)(this%PBC(idim),idim=1,ndim)
    else
       write(*,*)'Use default phase-space domain and boundary condition settings [Qmin=-1.0, Qmax=1.0, no PBC]? (T/F)'
       read(*,*)usedefault
       if(.not.usedefault)then
          do idim=1,ndim
             write(*,*)'Impose dimension',idim&
                  ,' periodic boundary conditions? (T/F)'
             read(*,*)this%PBC(idim)
             if(this%PBC(idim))then
                write(*,*)'Enter dimension',idim,' minimum value.'
                read(*,*)this%Qmin(idim)
                write(*,*)'Enter dimension',idim,' maximum value.'
                read(*,*)this%Qmax(idim)
             end if
          end do
       end if
    end if
        
    if(associated(this%Q))nullify(this%Q)
    if(associated(this%P))nullify(this%P)
    if(associated(this%F))nullify(this%F)
    allocate(this%Q(ndof),this%P(ndof),this%F(ndof),stat=ierr)
    if(ierr.NE.0)call stop('hb_init: failed phase-space memory allocation.')
    this%Q=0.0_double
    this%P=0.0_double
    this%F=0.0_double
    this%V=0.0_double
    this%T=0.0_double

    if(present(file))then
       read(unit,*)(this%Q(idof),idof=1,ndof)
       read(unit,*)(this%P(idof),idof=1,ndof)
    else
       write(*,*)'Use default phase-space settings [obsolete just enter T]? (T/F)'
       read(*,*)usedefault
       if(.not.usedefault)then
          do idof=1,ndof
             write(*,*)'Enter values for dof',idof
             write(*,*)'Q='
             read(*,*)this%Q(idof)
             write(*,*)'P='
             read(*,*)this%P(idof)
          end do
       end if
    end if

    if(associated(this%rmass))nullify(this%rmass)
    allocate(this%rmass(ndof),stat=ierr)
    if(ierr.NE.0)call stop('hb_init: failed rmass memory allocation.')
    do ipart=1,npart
       do idim=1,ndim
          idof=(ipart-1)*ndim+idim
          this%rmass(idof)=1.0_double/this%mass(ipart)
       end do
    end do

    this%temperature=298.15
    if(present(file))then
       read(unit,*)this%temperature
    else
       this%temperature=-1.0_double
       do while(this%temperature.LT.0)
          write(*,*)'Enter bath temperature in Kelvin.'
          read(*,*)this%temperature
          if(this%temperature.LT.0)write(*,*)'Temperature cannot be negative. Try again!'
       end do
    end if

    !normal modes
    if(associated(this%mode))nullify(this%mode)
    allocate(this%mode(ndof),stat=ierr)
    if(ierr.NE.0)call stop('hb_init: failed normal mode memory allocation.')
    this%mode=0.0_double

    !normal mode coords
    if(associated(this%W))nullify(this%W)
    allocate(this%W(ndof),stat=ierr)
    if(ierr.NE.0)call stop('hb_init: failed normal mode coordinate memory allocation.')
    this%W=0.0_double

    !Eigen vectors that transform Q into W coordinates
    if(associated(this%EigenVec))nullify(this%EigenVec)
    allocate(this%EigenVec(ndof,ndof),stat=ierr)
    if(ierr.NE.0)call stop('hb_init: failed Eigen vector memory allocation.')
    this%EigenVec=0.0_double
    do idof=1,ndof
       this%EigenVec(idof,idof)=1.0_double
    end do

    if(present(file))close(unit)
    this%initialized=.true.

    if(check(this).NE.0)call Stop('hb_init: failed check')

  end subroutine hb_init
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hb_kill(this)
    type(hb),intent(inout)::this

    call Note('Begin hb_kill')
    if(associated(this%mass))nullify(this%mass)
    if(associated(this%charge))nullify(this%charge)
    if(associated(this%Qmin))nullify(this%Qmin)
    if(associated(this%Qmax))nullify(this%Qmax)
    if(associated(this%PBC))nullify(this%PBC)
    if(associated(this%Q))nullify(this%Q)
    if(associated(this%P))nullify(this%P)
    if(associated(this%F))nullify(this%F)
    if(associated(this%rmass))nullify(this%rmass)
    if(associated(this%mode))nullify(this%mode)
    if(associated(this%W))nullify(this%W)
    if(associated(this%EigenVec))nullify(this%EigenVec)
    this%initialized=.false.
  end subroutine hb_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hb_update(this)
    type(hb),intent(inout)::this
    
    call Note('Begin hb_update.')

  end subroutine hb_update
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hb_resample(this)
    type(hb),intent(inout)::this
    
    call Note('Begin hb_resample.')
    
  end subroutine hb_resample
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hb_display(this,msg)
    type(hb),intent(in)::this
    character*(*),intent(in),optional::msg
    integer(long)::idof,jdof

    call Note('Begin hb_display.')
    if(check(this).NE.0)then
       call warn('hb_display: failed check','displaying nothing.')
       return
    end if

    write(Dunit,*)'____________________________________________________'
    write(Dunit,*)'-------------------     hb     ---------------------'
    if(present(msg))write(Dunit,*)msg
    write(Dunit,*)'----------------------------------------------------'
    write(Dunit,*)'npart=',this%npart
    write(Dunit,*)'ave mass=',sum(this%mass)/real(this%npart)
    write(Dunit,*)'Total charge=',sum(this%charge)
    write(Dunit,*)'ndim=',this%ndim
    write(Dunit,*)'Qmin=',this%Qmin
    write(Dunit,*)'Qmax=',this%Qmax
    write(Dunit,*)'PBC=',this%PBC
    write(Dunit,*)'V=',this%V
    write(Dunit,*)'V(hbar/cm)=',this%V/invcm
    write(Dunit,*)'T=',this%T
    write(Dunit,*)'T(hbar/cm)=',this%T/invcm
    write(Dunit,*)'thermostat(K)=',this%temperature
    write(Dunit,*)'beta=',1._double/(kb*this%temperature)
    write(Dunit,*)'beta(ps/hbar)=',1._double/(kb*this%temperature*ps)
    write(Dunit,*)'kT(hbar/cm)=',(kb*this%temperature)/invcm
    write(Dunit,*)'===================================================='

  end subroutine hb_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hb_save(this,file)
    type(hb),intent(in)::this
    character*(*),intent(in)::file
    
    integer(long)::ipart,idim,idof
    integer(long)::unit

    call Note('Begin hb_save.')
    call Note('input file= '//file)

    if(check(this).NE.0)then
       call warn('hb_save: failed check.','not saving object.')
    else
       unit=newunit()
       open(unit,file=file)
       
       write(unit,*)'hb'
       write(unit,*)this%npart
       write(unit,*)(this%mass(ipart),ipart=1,this%npart)
       write(unit,*)(this%charge(ipart),ipart=1,this%npart)
       write(unit,*)this%ndim
       write(unit,*)(this%Qmin(idim),idim=1,this%ndim)
       write(unit,*)(this%Qmax(idim),idim=1,this%ndim)
       write(unit,*)(this%PBC(idim),idim=1,this%ndim)
       write(unit,*)(this%Q(idof),idof=1,this%ndof)
       write(unit,*)(this%P(idof),idof=1,this%ndof)
       write(unit,*)this%temperature

       close(unit)
    end if
  end subroutine hb_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short) function hb_check(this)
    type(hb),intent(in)::this
    
    integer(short)::ierr
    integer(long)::idof,ipart,idim,jdof
    real(double)::norm

    call Note('Checking hb.')

    hb_check=0
    
    !logical::initialized
    if(.not.this%initialized)then
       call Warn('hb_check: classical primitive not initialized.')
       hb_check=1
       return
    end if
    
    !integer(long)::npart
    if(this%npart.LE.0)then
       call Warn('hb_check: number of particles is less than 0.')
       hb_check=1
       return
    end if
    if(this%npart.NE.this%npart)then
       call Warn('hb_check: number of particles is not a number.')
       hb_check=1
       return
    end if
    if(this%npart.GE.huge(this%npart))then
       call Warn('hb_check: number of particles is too large.')
       hb_check=1
       return
    end if
    
    !real(double),dimension(:),pointer::mass
    if(.not.associated(this%mass))then
       call Warn('hb_check: particle mass memory not associated')
       hb_check=1
       return
    end if
    if(size(this%mass).ne.this%npart)then
       call Warn('hb_check: size of particle mass array not equal to number of particles')
       hb_check=1
       return
    end if
    ierr=0
    do ipart=1,this%npart
       if(this%mass(ipart).NE.this%mass(ipart))ierr=1
       if(abs(this%mass(ipart)).GE.huge(this%mass(ipart)))ierr=1
    end do
    if(ierr.NE.0)then
       call Warn('hb_check: particle masses have bad values')
       hb_check=1
       return
    end if

    !real(double),dimension(:),pointer::charge
    if(.not.associated(this%charge))then
       call Warn('hb_check: particle charge memory not associated')
       hb_check=1
       return
    end if
    if(size(this%charge).ne.this%npart)then
       call Warn('hb_check: size of particle charge array not equal to number of particles')
       hb_check=1
       return
    end if
    ierr=0
    do ipart=1,this%npart
       if(this%charge(ipart).NE.this%charge(ipart))ierr=1
       if(abs(this%charge(ipart)).GE.huge(this%charge(ipart)))ierr=1
    end do
    if(ierr.NE.0)then
       call Warn('hb_check: particle chargees have bad values')
       hb_check=1
       return
    end if
    
    !integer(byte)::ndim
    if(this%ndim.NE.this%ndim)then
       call Warn('hb_check: value of phase-space dimensions is not a number.')
       hb_check=1
       return
    end if
    if(this%ndim.LE.0)then
       call Warn('hb_check: number of phase-space dimensions is less than 0.')
       hb_check=1
       return
    end if
    if(this%ndim.GE.huge(this%ndim))then
       call Warn('hb_check: number of phase-space dimensions is too large.')
       hb_check=1
       return
    end if
    
    !integer(long)::ndof
    if(this%ndof.NE.this%ndof)then
       call Warn('hb_check: number of degrees of freedom not a number.')
       hb_check=1
       return
    end if
    if(this%ndof.LE.0_long)then
       call Warn('hb_check: number of degrees of freedom not a positive integer.')
       hb_check=1
       return
    end if
    if(this%ndof.GE.huge(this%ndof))then
       call Warn('hb_check: number of degrees of freedom too large.')
       hb_check=1
       return
    end if
        
    !real(double),dimension(:),pointer::rmass
    if(.not.associated(this%rmass))then
       call Warn('hb_check: rmass memory not associated')
       hb_check=1
       return
    end if
    if(size(this%rmass).ne.this%ndof)then
       call Warn('hb_check: size of rmass array not equal to number of degrees of freedom')
       hb_check=1
       return
    end if
    ierr=0
    do idof=1,this%ndof
       if(this%rmass(idof).NE.this%rmass(idof))ierr=1
       if(abs(this%rmass(idof)).GE.huge(this%rmass(idof)))ierr=1
    end do
    if(ierr.NE.0)then
       call Warn('hb_check: rmass has bad values')
       hb_check=1
       return
    end if
        
    !real(double),dimension(:),pointer::Qmin,Qmax
    ierr=0
    if(.not.associated(this%Qmin))ierr=1
    if(.not.associated(this%Qmax))ierr=1
    if(ierr.NE.0)then
       call Warn('hb_check: space domain association error.')
       hb_check=1
       return
    end if
    ierr=0
    if(size(this%Qmin).NE.this%ndim)ierr=1
    if(size(this%Qmax).NE.this%ndim)ierr=1
    if(ierr.NE.0)then
       call Warn('hb_check: space dimensions not equal to ndim.')
       hb_check=1
       return
    end if
    ierr=0
    if(any(this%Qmin.NE.this%Qmin))then
       call Warn('hb_check: Qmin has NAN values.')
       hb_check=1
       return
    end if
    if(any(abs(this%Qmin).GE.huge(this%Qmin)))then
       call Warn('hb_check: Qmin has huge values.')
       hb_check=1
       return
    end if
    if(any(this%Qmax.NE.this%Qmax))then
       call Warn('hb_check: Qmax has NAN values.')
       hb_check=1
       return
    end if
    if(any(abs(this%Qmax).GE.huge(this%Qmax)))then
       call Warn('hb_check: Qmax has huge values.')
       hb_check=1
       return
    end if
    
    !logical,dimension(:),pointer::PBC
    if(.not.associated(this%PBC))then
       call Warn('hb_check: Periodic Boundary Condition memory association error.')
       hb_check=1
       return
    end if
    if(size(this%PBC).NE.this%ndim)then
       call Warn('hb_check: Periodic Boundary Condition memory dimensions not equal to ndim.')
       hb_check=1
       return
    end if
    do idim=1,this%ndim
       if(this%PBC(idim).NEQV.this%PBC(idim))then
          call Warn('hb_check: Periodic Boundary Conditions have bad values.')
          hb_check=1
          return
       end if
    end do
    
    !real(double),dimension(:),pointer::Q,P,F
    ierr=0
    if(.not.associated(this%Q))ierr=1
    if(.not.associated(this%P))ierr=1
    if(.not.associated(this%F))ierr=1
    if(ierr.NE.0)then
       call Warn('hb_check: phase-space association error.')
       hb_check=1
       return
    end if
    ierr=0
    if(size(this%Q).NE.this%ndof)ierr=1
    if(size(this%P).NE.this%ndof)ierr=1
    if(size(this%F).NE.this%ndof)ierr=1
    if(ierr.NE.0)then
       call Warn('hb_check: phase-space memory size error.')
       hb_check=1
       return
    end if

    if(any(this%Q.NE.this%Q))then
       call Warn('hb_check: coordinates have NAN values.')
       hb_check=1
       return
    end if
    if(any(this%P.NE.this%P))then
       call Warn('hb_check: momenta have NAN values.')
       hb_check=1
       return
    end if
    if(any(this%F.NE.this%F))then
       call Warn('hb_check: forcefield has NAN values.')
       hb_check=1
       return
    end if

    if(any(this%Q.GE.huge(this%Q)))then
       call Warn('hb_check: coordinates have huge values.')
       hb_check=1
       return
    end if
    if(any(this%P.GE.huge(this%P)))then
       call Warn('hb_check: momenta have huge values.')
       hb_check=1
       return
    end if
    if(any(this%F.GE.huge(this%F)))then
       call Warn('hb_check: forcefield has huge values.')
       hb_check=1
       return
    end if
    
    !real(double)::V
    if(this%V.NE.this%V)then
       call Warn('hb_check: potential energy is not a number.')
       hb_check=1
       return
    end if
    if(this%V.GE.huge(this%V))then
       call Warn('hb_check: potential energy is huge.')
       hb_check=1
       return
    end if
    
    !real(double)::T
    if(this%T.NE.this%T)then
       call Warn('hb_check: kinetic energy is not a number.')
       hb_check=1
       return
    end if
    if(this%T.GE.huge(this%T))then
       call Warn('hb_check: kinetic energy is huge.')
       hb_check=1
       return
    end if
    
    !real(double)::temperature
    if(this%temperature.NE.this%temperature)then
       call Warn('hb_check: temperature is not a number.')
       hb_check=1
       return
    end if
    if(this%temperature.GE.huge(this%temperature))then
       call Warn('hb_check: temperature is huge.')
       hb_check=1
       return
    end if

    !normal modes
    if(.not.associated(this%mode))then
       call Warn('hb_check: Normal mode association error.')
       hb_check=1
       return
    end if
    if(size(this%mode).NE.this%ndof)then
       call Warn('hb_check: Normal mode memory size error.')
       hb_check=1
       return
    end if
    if(any(this%mode.NE.this%mode))then
       call Warn('hb_check: NAN normal mode values.')
       hb_check=1
       return
    end if
    if(any(abs(this%mode).GE.huge(this%mode)))then
       call Warn('hb_check: Huge normal mode values.')
       hb_check=1
       return
    end if

    !normal mode coordinates
    if(.not.associated(this%W))then
       call Warn('hb_check: Normal mode coodinate association error.')
       hb_check=1
       return
    end if
    if(size(this%W).NE.this%ndof)then
       call Warn('hb_check: Normal mode coordinate memory size error.')
       hb_check=1
       return
    end if
    if(any(this%W.NE.this%W))then
       call Warn('hb_check: NAN normal mode coordinate values.')
       hb_check=1
       return
    end if
    if(any(abs(this%W).GE.huge(this%W)))then
       call Warn('hb_check: Huge normal mode coordinate values.')
       hb_check=1
       return
    end if

    !eigen vectors
    if(.not.associated(this%EigenVec))then
       call Warn('hb_check: Eigen vector association error.')
       hb_check=1
       return
    end if
    if(size(this%EigenVec).NE.this%ndof**2)then
       call Warn('hb_check: Eigen vector memory size error.')
       hb_check=1
       return
    end if
    if(any(this%EigenVec.NE.this%EigenVec))then
       call Warn('hb_check: NAN eigen vector values.')
       hb_check=1
       return
    end if
    if(any(abs(this%EigenVec).GE.huge(this%EigenVec)))then
       call Warn('hb_check: Huge eigen vector values.')
       hb_check=1
       return
    end if
    do idof=1,this%ndof
       norm=sum(this%EigenVec(:,idof)**2)
       if(norm.NE.1_double)then
          call Warn('hb_check: Sum of eigen vector '//trim(int2str(idof))//' not equal to 1.')
          hb_check=1
          return
       end if
    end do
    
  end function hb_check
  
end module hb_class
