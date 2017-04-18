!-------------------------------
!>\brief
!! Transition dipoles with Gaussian site disorder
!!\details
!! Defines an inhomogeniously broadend set of exctionic states 
!! with site transition dipoles, dipole-dipole coupling, and Gaussian
!! site disorder.
!<--------------------------------  
module dipoles_class
  use type_kinds
  use ErrorLog
  use string
  use atomicunits
  use rand
  use hs_class
  use textgraphs
  use filemanager
  use outputdisplay
  implicit none
  private

  public::dipoles
  public::new, display, save, update, resample, kill, check

  !-----------------------------------------
  !>\brief
  !! derived quantum subsystem object
  !!\sa hs
  !<----------------------------------------
  type dipoles 
!>\sa Evar
!<----

!left off here

     logical::initialized=.false.
     type(hs)::hs !>quantum primitive
     type(hs)::hs0!>disorder-less quantum primitive
     integer(short)::ndim!>Dimensionality of transition dipoles and site center of mass
     real(double),dimension(:,:),pointer::dielectric!>site to site optical dielectric screening
     real(double)::Gvar!>groundstate energy variance
     real(double)::Gdisorder!>current realization of groundstate disorder
     real(double),dimension(:,:),pointer::Evar!>variance in diabatic hamiltonian
     real(double),dimension(:,:),pointer::Disorder!>current realization of diabatic hamiltonian disorder
     real(double),dimension(:,:),pointer::mu!>Site Transition dipole vector
     real(double),dimension(:,:),pointer::Q!>Site center of mass
  end type dipoles

  interface new
     module procedure dipoles_init
  end interface

  interface kill
     module procedure dipoles_kill
  end interface

  interface display
     module procedure dipoles_display
  end interface

  interface save
     module procedure dipoles_save
  end interface

  interface update
     module procedure dipoles_update
  end interface

  interface resample
     module procedure dipoles_resample
  end interface

  interface check
     module procedure dipoles_check
  end interface

contains
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine dipoles_init(this,file)
    type(dipoles),intent(inout)::this
    character*(*),intent(in),optional::file
    character(len=title)::filetype
    character(len=path)::infile

    integer::unit
    logical::usedefault,usedunit

    integer(long)::istate,jstate
    integer(short)::idim
    !real(double)::R
    !real(double),allocatable,dimension(:)::rhat



    call Note('Begin dipoles_init.')
    if(present(file))call Note('input file= '//file)

    !check if input file is present and valid
    if(present(file))then
       if(check(file).EQ.1)call stop('dipoles_init')
       unit=newunit()
       open(unit,file=file)
       read(unit,*)filetype
       filetype=adjustl(filetype)
       if(trim(filetype).NE.'dipoles')call Stop('dipoles_init: input file is not valid.')
    end if

    if(present(file))then
       read(unit,*)infile
       infile=adjustl(infile)
       call new(this%hs0,trim(infile))
       call new(this%hs,trim(infile))
    else
       call new(this%hs0)
       call save(this%hs0,file='temp.dipoles.hs')
       call new(this%hs,file='temp.dipoles.hs')
    end if

    if(present(file))then
       read(unit,*)this%Gvar
    else
       write(*,*)'Enter the variance in groundstate energy.'
       read(*,*)this%Gvar
    end if

    if(associated(this%dielectric))nullify(this%dielectric)
    allocate(this%dielectric(this%hs%nstate,this%hs%nstate))
    this%dielectric=1.0_double
    if(present(file))then
       read(unit,*)((this%dielectric(istate,jstate),jstate=1,this%hs%nstate),istate=1,this%hs%nstate)
    else
       usedefault=.true.
       write(*,*)'Use default optical dielectric (no dipole screening)? (T/F)'
       read(*,*)usedefault
       if(usedefault)then
          this%dielectric=1.0_double
       else
          write(*,*)'Enter optical dielectric constant:'
          do istate=1,this%hs%nstate
             do jstate=istate,this%hs%nstate
                write(*,*)'between states',istate,jstate
                read(*,*)this%dielectric(istate,jstate)
                this%dielectric(jstate,istate)=this%dielectric(istate,jstate)
             end do
          end do
       end if
    end if

    if(present(file))then
       read(unit,*)this%ndim
    else
       write(*,*)'Enter number of dipole dimensions.'
       read(*,*)this%ndim
    end if

    if(associated(this%mu))nullify(this%mu)
    allocate(this%mu(this%hs%nstate,this%ndim))
    this%mu=0.0_double
    if(present(file))then
       read(unit,*)((this%mu(istate,idim),idim=1,this%ndim),istate=1,this%hs%nstate)
    else
       usedefault=.true.
       write(*,*)'Use default site dipoles (unit Debye vector in principal dimension)? (T/F)'
       read(*,*)usedefault
       if(usedefault)then
          this%mu(:,1)=1.0_double*Debye
       else
          write(*,*)'Enter dipole component for.'
          do istate=1,this%hs%nstate
             do idim=1,this%ndim
                write(*,*)'dimension',idim,'in state',istate
                read(*,*)this%mu(istate,idim)
             end do
          end do
       end if
    end if

    if(associated(this%Q))nullify(this%Q)
    allocate(this%Q(this%hs%nstate,this%ndim))
    this%Q=0.0_double
    if(present(file))then
       read(unit,*)((this%Q(istate,idim),idim=1,this%ndim),istate=1,this%hs%nstate)
    else
       usedefault=.true.
       write(*,*)'Use default site locations (atomic unit increments along principal dimension)? (T/F)'
       read(*,*)usedefault
       if(usedefault)then
          this%Q(:,1)=1.0_double
       else
          write(*,*)'Enter dipole origin coordinates.'
          do istate=1,this%hs%nstate
             do idim=1,this%ndim
                write(*,*)'dimension',idim,'in state',istate
                read(*,*)this%Q(istate,idim)
             end do
          end do
       end if
    end if

    if(associated(this%Evar))nullify(this%Evar)
    allocate(this%Evar(this%hs%nstate,this%hs%nstate))
    this%Evar=0.0_double
    if(present(file))then
       read(unit,*)((this%Evar(istate,jstate),jstate=1,this%hs%nstate),istate=1,this%hs%nstate)
    else
       usedefault=.true.
       write(*,*)'Use default Diabatic energy variance (no variance)? (T/F)'
       read(*,*)usedefault
       if(usedefault)then
          this%Evar=0.0_double
       else
          write(*,*)'Enter Variance in Diabatic energy:'
          do istate=1,this%hs%nstate
             do jstate=istate,this%hs%nstate
                write(*,*)'for state',istate,jstate
                read(*,*)this%Evar(istate,jstate)
                this%Evar(jstate,istate)=this%Evar(istate,jstate)
             end do
          end do
       end if
    end if

    if(associated(this%Disorder))nullify(this%Disorder)
    allocate(this%Disorder(this%hs%nstate,this%hs%nstate))
    this%initialized=.true.

!!$    !compute dipole-dipole coupling
!!$    allocate(rhat(this%ndim))
!!$    do istate=1,this%hs%nstate
!!$       do jstate=istate,this%hs%nstate
!!$          Rhat=this%Q(istate,:)-this%Q(jstate,:)
!!$          R=sqrt(sum(rhat**2))
!!$          if(R.NE.0.0_double.and.this%dielectric(istate,jstate).NE.0.0_double)then
!!$             Rhat=Rhat/R
!!$             this%hs0%diabat(istate,jstate)=kC*(sum(this%mu(istate,:)*this%mu(jstate,:))&
!!$                  -(3.0_double*sum(this%mu(istate,:)*rhat)&
!!$                  *sum(this%mu(jstate,:)*rhat)))&
!!$                  /(this%dielectric(istate,jstate)*R*R*R)
!!$             this%hs0%diabat(jstate,istate)=this%hs0%diabat(istate,jstate)
!!$          end if
!!$       end do
!!$    end do
!!$    deallocate(rhat)

    this%Disorder=0.0_double
    if(present(file))then
       read(unit,*)this%Gdisorder
       read(unit,*)((this%Disorder(istate,jstate),jstate=1,this%hs%nstate),istate=1,this%hs%nstate)
       call update(this)
    else
       call resample(this)
    end if

    if(present(file))close(unit)

    if(check(this).EQ.1)call stop('dipoles_init, final object check')

  end subroutine dipoles_init
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine dipoles_display(this,msg)
    type(dipoles),intent(in)::this
    character*(*),intent(in),optional::msg
    integer(long)::istate,jstate

    real(double),allocatable::Q(:,:),Mu(:,:)
    real(double)::Rx(3,3),Ry(3,3),Rz(3,3)
    real(double)::Ax,Ay,Az

    call Note('Begin dipoles_display.')
    if(check(this).NE.0)then
       call warn('dipoles_display: failed check','displaying nothing.')
       return
    end if

    Ax=pi/4
    Ay=0
    Az=pi/4

    Rx=0_double
    Rx(1,1)=1_double
    Rx(2,2)=cos(Ax)
    Rx(2,3)=-sin(Ax)
    Rx(3,2)=sin(Ax)
    Rx(3,3)=cos(Ax)

    Rz=0_double
    Rz(3,3)=1_double
    Rz(1,1)=cos(Az)
    Rz(1,2)=-sin(Az)
    Rz(2,1)=sin(Az)
    Rz(2,2)=cos(Az)

    if(allocated(Q))deallocate(Q)
    allocate(Q(this%hs%nstate,this%ndim))
    Q=this%Q

    if(allocated(Mu))deallocate(Mu)
    allocate(Mu(this%hs%nstate,this%ndim))
    Mu=this%Mu

    do istate=1,this%hs%nstate
       Q(istate,:)=matmul(Rx,Q(istate,:))
       Mu(istate,:)=matmul(Rx,Mu(istate,:))
       Q(istate,:)=matmul(Rz,Q(istate,:))
       Mu(istate,:)=matmul(Rz,Mu(istate,:))
    end do


    write(Dunit,*)'____________________________________________________'
    write(Dunit,*)'-------------------  dipoles   ---------------------'
    if(present(msg))write(Dunit,*)msg
    write(Dunit,*)'----------------------------------------------------'
    write(Dunit,*)'Groundstate energy stdev (1/cm)=',sqrt(this%Gvar)/invcm
    write(Dunit,*)'ndim =',this%ndim
    call display(this%Q(:,1:2),this%Mu(:,1:2),msg='Dipoles x-y projection')
    call display(Q(:,1:2),Mu(:,1:2),msg='Rotated view (z points in x+y+z direction)')
    call display(this%dielectric-1._double,msg='Optical dielectric matrix')
    call display(this%Evar,msg='Variance in site energy')
    call display(this%Disorder,msg='Realization of Disorder')
    call display(this%hs%diabat,msg='Quantum subsystem Hamiltonain')
    call display(this%hs0%diabat,msg='Disorder-less quantum subsystem Hamiltonain')
    write(Dunit,*)'===================================================='

    if(allocated(Q))deallocate(Q)
    if(allocated(Mu))deallocate(Mu)
  end subroutine dipoles_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine dipoles_save(this,file)
    type(dipoles),intent(in)::this
    character*(*),intent(in)::file

    integer(short)::unit
    integer(long)::istate,jstate
    integer(short)::idim
    logical::usedunit      

    call note('Begin dipoles_save.')
    call Note('input file= '//file)
    if(check(this).NE.0)then
       call warn('dipoles_save: failed check.','not saving object.')
    else
       unit=newunit()
       open(unit,file=file)
       write(unit,*)'dipoles'
       write(unit,*)quote(file//'.hs')
       call save(this%hs0,file//'.hs')
       write(unit,*)this%Gvar
       write(unit,*)((this%dielectric(istate,jstate),jstate=1,this%hs%nstate),istate=1,this%hs%nstate)
       write(unit,*)this%ndim
       write(unit,*)((this%mu(istate,idim),idim=1,this%ndim),istate=1,this%hs%nstate)
       write(unit,*)((this%Q(istate,idim),idim=1,this%ndim),istate=1,this%hs%nstate)
       write(unit,*)((this%Evar(istate,jstate),jstate=1,this%hs%nstate),istate=1,this%hs%nstate)
       write(unit,*)this%Gdisorder
       write(unit,*)((this%Disorder(istate,jstate),jstate=1,this%hs%nstate),istate=1,this%hs%nstate)
       close(unit)
    end if
  end subroutine dipoles_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-
  subroutine dipoles_update(this)
    type(dipoles),intent(inout)::this
    real(double),allocatable,dimension(:)::rhat
    integer(long)::istate,jstate
    real(double)::R,DD

    call Note('Begin dipoles_update.')

    !compute dipole-dipole coupling
    if(allocated(rhat))deallocate(rhat)
    allocate(rhat(this%ndim))

    do istate=1,this%hs%nstate
       do jstate=istate,this%hs%nstate
          if(istate.EQ.jstate)then
             this%hs%diabat(istate,jstate)=this%hs0%diabat(istate,jstate)
          else
             DD=0.0_double
             Rhat=this%Q(istate,:)-this%Q(jstate,:)
             R=sqrt(sum(rhat**2))
             if(R.NE.0.0_double.and.this%dielectric(istate,jstate).NE.0.0_double)then
                Rhat=Rhat/R
                DD=kC*(sum(this%mu(istate,:)*this%mu(jstate,:))&
                     -(3.0_double*sum(this%mu(istate,:)*rhat)&
                     *sum(this%mu(jstate,:)*rhat)))&
                     /(this%dielectric(istate,jstate)*R*R*R)
             end if
             this%hs%diabat(istate,jstate)=this%hs0%diabat(istate,jstate)+DD
             this%hs%diabat(jstate,istate)=this%hs%diabat(istate,jstate)
          end if
       end do
    end do
    if(allocated(rhat))deallocate(rhat)
    
    !add static disorder
    this%hs%Eg=this%Gdisorder+this%hs0%Eg
    this%hs%diabat=this%hs%diabat+this%Disorder
    
    !recompute eigen states
    this%hs%EigenVec=cmplx(this%hs%diabat)
    call diagonalize(this%hs%nstate,this%hs%EigenVec,this%hs%EigenVal)

  end subroutine dipoles_update
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-
  subroutine dipoles_resample(this)
    type(dipoles),intent(inout)::this
    integer(long)::istate,jstate

    call Note('Begin dipoles_resample.')
    !call resample(this%hs)

    !resample static disorder
    this%Gdisorder=sqrt(this%Gvar)*gran()
    this%Disorder=0.0_double
    do istate=1,this%hs%nstate
       do jstate=1,this%hs%nstate
          this%Disorder(istate,jstate)=sqrt(this%Evar(istate,jstate))*gran()
       end do
    end do

    call update(this)    

  end subroutine dipoles_resample
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine dipoles_kill(this)
    type(dipoles),intent(inout)::this
 
    call note('Begin dipoles_kill.')

    call kill(this%hs0)
    call kill(this%hs)
    if(associated(this%dielectric))nullify(this%dielectric)
    if(associated(this%mu))nullify(this%mu)
    if(associated(this%Q))nullify(this%Q)
    if(associated(this%Evar))nullify(this%Evar)
    if(associated(this%Disorder))nullify(this%Disorder)
    this%initialized=.false.

  end subroutine dipoles_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short)function dipoles_check(this)
    type(dipoles),intent(in)::this

    integer(long)::istate,jstate
    integer(short)::idim
    integer(short)::ierr

    call Note('Checking dipoles.')

    dipoles_check=0

    !logical::initialized=.false.
    if(.not.this%initialized)then
       call Warn('dipoles_check: dipoles object not initialized.')
       dipoles_check=1
       return
    end if

    !type(hs)::hs0
    if(check(this%hs0).EQ.1)then
       call Warn('dipoles_check: hs0 failed check')
       dipoles_check=1
       return
    end if

    !type(hs)::hs
    if(check(this%hs).EQ.1)then
       call Warn('dipoles_check: hs failed check')
       dipoles_check=1
       return
    end if

    !integer::ndim
    if(this%ndim.NE.this%ndim)then
       call Warn('dipoles_check: ndim not a number.')
       dipoles_check=1
       return
    end if
    if(abs(this%ndim).GE.huge(this%ndim))then
       call Warn('dipoles_check: ndim is too big.')
       dipoles_check=1
       return
    end if
    if(this%ndim.LE.0)then
       call Warn('dipoles_check: ndim not a positive integer.')
       dipoles_check=1
       return
    end if

    !real(double)::Gvar
    if(this%Gvar.NE.this%Gvar)then
       call Warn('dipoles_check: Gvar not a number.')
       dipoles_check=1
       return
    end if
    if(abs(this%Gvar).GE.huge(this%Gvar))then
       call Warn('dipoles_check: Gvar is too big.')
       dipoles_check=1
       return
    end if

    !real(double)::Gdisorder
    if(this%Gdisorder.NE.this%Gdisorder)then
       call Warn('dipoles_check: Gdisorder not a number.')
       dipoles_check=1
       return
    end if
    if(abs(this%Gdisorder).GE.huge(this%Gdisorder))then
       call Warn('dipoles_check: Gdisorder is too big.')
       dipoles_check=1
       return
    end if

    !real(double),dimension(:,:),pointer::dielectric
    if(.not.associated(this%dielectric))then
       call Warn('dipoles_check: dielectric matrix memory not associated.')
       dipoles_check=1
       return
    end if
    if(size(this%dielectric).NE.this%hs%nstate**2)then
       call Warn('dipoles_check: number of dielectric matrix elements must equal number of quantum states squared.')
       dipoles_check=1
       return
    end if
    ierr=0
    do istate=1,this%hs%nstate
       do jstate=1, this%hs%nstate
          if(this%dielectric(istate,jstate).NE.this%dielectric(istate,jstate))ierr=1
          if(abs(this%dielectric(istate,jstate)).GE.huge(this%dielectric(istate,jstate)))ierr=1
       end do
    end do
    if(ierr.NE.0)then
       call Warn('dipoles_check: dielectric matrix elements have bad values.')
       dipoles_check=1
       return
    end if

    !real(double),dimension(:,:),pointer::Evar
    if(.not.associated(this%Evar))then
       call Warn('dipoles_check: Evar matrix memory not associated.')
       dipoles_check=1
       return
    end if
    if(size(this%Evar).NE.this%hs%nstate**2)then
       call Warn('dipoles_check: number of Evar matrix elements must equal number of quantum states squared.')
       dipoles_check=1
       return
    end if
    ierr=0
    do istate=1,this%hs%nstate
       do jstate=1, this%hs%nstate
          if(this%Evar(istate,jstate).NE.this%Evar(istate,jstate))ierr=1
          if(abs(this%Evar(istate,jstate)).GE.huge(this%Evar(istate,jstate)))ierr=1
       end do
    end do
    if(ierr.NE.0)then
       call Warn('dipoles_check: Evar matrix elements have bad values.')
       dipoles_check=1
       return
    end if

    !real(double),dimension(:,:),pointer::Disorder
    if(.not.associated(this%Disorder))then
       call Warn('dipoles_check: Disorder matrix memory not associated.')
       dipoles_check=1
       return
    end if
    if(size(this%Disorder).NE.this%hs%nstate**2)then
       call Warn('dipoles_check: number of Disorder matrix elements must equal number of quantum states squared.')
       dipoles_check=1
       return
    end if
    ierr=0
    do istate=1,this%hs%nstate
       do jstate=1, this%hs%nstate
          if(this%Disorder(istate,jstate).NE.this%Disorder(istate,jstate))ierr=1
          if(abs(this%Disorder(istate,jstate)).GE.huge(this%Disorder(istate,jstate)))ierr=1
       end do
    end do
    if(ierr.NE.0)then
       call Warn('dipoles_check: Disorder matrix elements have bad values.')
       dipoles_check=1
       return
    end if

    !real(double),dimension(:,:),pointer::mu
    if(.not.associated(this%mu))then
       call Warn('dipoles_check: Mu matrix memory not associated.')
       dipoles_check=1
       return
    end if
    if(size(this%mu).NE.this%hs%nstate*this%ndim)then
       call Warn('dipoles_check: Mu matrix must have dimensions quantum states by ndim .')
       dipoles_check=1
       return
    end if
    ierr=0
    do istate=1,this%hs%nstate
       do idim=1, this%ndim
          if(this%mu(istate,idim).NE.this%mu(istate,idim))ierr=1
          if(abs(this%mu(istate,idim)).GE.huge(this%mu(istate,idim)))ierr=1
       end do
    end do
    if(ierr.NE.0)then
       call Warn('dipoles_check: Mu matrix elements have bad values.')
       dipoles_check=1
       return
    end if

    !real(double),dimension(:,:),pointer::Q
    if(.not.associated(this%Q))then
       call Warn('dipoles_check: Q matrix memory not associated.')
       dipoles_check=1
       return
    end if
    if(size(this%Q).NE.this%hs%nstate*this%ndim)then
       call Warn('dipoles_check: Q matrix must have dimensions quantum states by ndim .')
       dipoles_check=1
       return
    end if
    ierr=0
    do istate=1,this%hs%nstate
       do idim=1, this%ndim
          if(this%Q(istate,idim).NE.this%Q(istate,idim))ierr=1
          if(abs(this%Q(istate,idim)).GE.huge(this%Q(istate,idim)))ierr=1
       end do
    end do
    if(ierr.NE.0)then
       call Warn('dipoles_check: Q matrix elements have bad values.')
       dipoles_check=1
       return
    end if

  end function dipoles_check

end module dipoles_class

