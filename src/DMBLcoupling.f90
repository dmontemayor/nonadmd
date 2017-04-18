!-----------------------------------------------------------------------------
!> \brief
!!Diagonal Multi-Bath Linear Coupling.
!! \details
!! This coupling scheme is a Caldeira-Leggett model, coupling the normal modes of the classical subsystem linearly to the diagonal states of the quantum subsystem. Assumes the Hamiltonian
!!        \f[  H=hs+hb+H_c+ct \f]
!! where
!!         \f[ hs=\frac{\hat p^2}{2m}+V(\hat q) \f]
!! ,projected on to a diabatic basis of hs (ie diabats attribute of hs).
!!        \f[ hb=V(Q) \f]
!! , projected onto the normal modes of hb (ie modes attribute of hb).
!! Each mode is assigned to one of \f$ N \f$ dissipative baths
!! with bath \f$ n \f$ having \f$ K_n \f$ modes.
!! Each bath is coupled uniquely to the diagonal basis states
!! which are assumed constant in \f$ \hat q \f$ such that, 
!!   \f[ \hat q_n=\sum_i f_{i,n} |i \rangle\langle i|. \f]
!!   \f[ \hat q_n^2=\sum_i f_{i,n}^2 |i \rangle\langle i|. \f]
!! The bilinear coupling term is written
!!   \f[ H_c=-\sum_n^N \sum_k^{K_n} c_{k,n} Q_{k,n} \hat q_n \f]
!! with counter term
!!   \f[ ct=\sum_n^N \sum_k^{K_n} \hat q_n^2 \frac{c_{k,n}^2}{2m_{k,n}\omega_{k,n}^2}. \f]
!! The spectral density of bath \f$ n \f$ is defined
!! \f[ J_n(\omega)=\frac{\pi}{\hbar}\sum_k^{K_n}\frac{c_{k,n}^2}{2m_{k,n}\omega_{k,n}}\delta_{\omega-\omega_{k,n}}.\f]
!! The re-organization energy of bath \f$ n \f$ is defined
!!   \f[ \lambda_n=\frac{2}{\pi} \int_0^\infty \frac{J_n(\omega)}{\omega} d\omega \f]
!! such that,
!!   \f[ c_{k,n}^2=\frac{\hbar \lambda_n m_{k,n} \omega_{k,n}^2}{K_n}\f]
!! and 
!!   \f[ ct=2 \hbar \sum_n^N \lambda_n \hat q_n^2 .\f]
!<----------------------------------------------------------------------------
!>\authors
!! Daniel Montemayor
!!
!!\date
!! Aug 2011, Sept 2011
!<----------------------------------------------------------------------------
!>\todo
!!Change from rmass (particle basis) to normal mode mass weighting.
!!\todo
!!Check that coupling and counter terms are treated as defined in details section.
!<-----------------------------------------------------------------------------
! !>\bug
! !<-----------------------------------------------------------------------------
!Change Log:
!
!=== v1.1 Sept 2011 ===
! - removed automatic this%f normalization
! 
! === v1.0 Aug 2011 ===
! + beta 1.0 compatable
! + defines multiple baths with fractions of the total normal modes
!       associated to each bath
! + defines weights(g) using spectral denisty def 
!       J(w)=pi/2*sum_k delta(w_k - w)*(g_k**2)/w_k
!       and reoganization energy defined 
!       lambda=1/pi * integral_{0}^{infinity} J(w)/w dw
! + read input reorganization energies 
! + maps multiple baths to input reorganization energies
! + added counterterm Hr=quadratic in quantum dof 
!    times sum_{normal modes k} (g_k/w_k)**2/(2m)
!-----------------------------------------------------------------------------

module DMBLcoupling_class
  use type_kinds
  use atomicunits
  use string
  use ErrorLog
  use filemanager
  use hc_class
  use quantum_class
  use classical_class
  use textgraphs
  implicit none
  private

  public::DMBLcoupling
  public::new,kill,update,resample,display,save,check

  type DMBLcoupling
     logical::initialized=.false.
     type(hc)::hc
     integer(long)::nbath
     real(double),dimension(:),pointer::lambda,X!,wc
     integer(long),dimension(:),pointer::map,bdof
     real(double),dimension(:,:),pointer::g
     real(double),dimension(:,:),pointer::f
  end type DMBLcoupling

  interface new
     module procedure DMBLcoupling_init
  end interface

  interface kill
     module procedure DMBLcoupling_kill
  end interface

  interface update
     module procedure DMBLcoupling_update
  end interface

  interface resample
     module procedure DMBLcoupling_resample
  end interface

  interface display
     module procedure DMBLcoupling_display
  end interface

  interface save
     module procedure DMBLcoupling_save
  end interface

  interface check
     module procedure DMBLcoupling_check
  end interface

contains
!======================================================================
  subroutine DMBLcoupling_init(this,qs,cs,file)
    type(DMBLcoupling),intent(inout)::this
!!$    type(hs),intent(inout),target::hs
!!$    type(hb),intent(inout),target::hb
    type(quantum),intent(inout),target::qs
    type(classical),intent(inout),target::cs
    character*(*),intent(in),optional::file      
    character(len=path)::filename
    character(len=title)::filetype
    integer(short)::ierr
    integer(long)::unit
    logical::usedefault,pass

    integer(long)::idof,istate,ibath,iden,bdof
    real(double)::norm


    call Note('Begin DMBLcoupling_init.')
    filename='unknown'
    filetype='unknown'

!!$    if(check(hs))call stop('DMBLcoupling_init, hs failed check')
!!$    if(check(hb))call stop('DMBLcoupling_init, hb failed check')
    if(check(qs).EQ.1)call stop('DMBLcoupling_init, qs failed check')
    if(check(cs).EQ.1)call stop('DMBLcoupling_init, cs failed check')

    if(present(file))call Note('input file= '//trim(file))


    if(present(file))then
       if(check(trim(file)).EQ.1)&
            call stop('DMBLcoupling_init, cannot find input file')
       
       unit=newunit()
       open(unit,file=file)
       read(unit,*)filetype
       filetype=adjustl(filetype)
       if(trim(filetype).EQ.'DMBLcoupling')then
          read(unit,*)filename
          filename=adjustl(filename)
!!$          call new(this%hc,hs,hb,trim(filename))
          call new(this%hc,qs,cs,trim(filename))
       else
          call Stop('DMBLcoupling_init: input file not valid type.','Check input file format for DMBLcouling object.')
       end if
    else
!!$       call new(this%hc,hs,hb)
       call new(this%hc,qs,cs)
    end if
    if(check(this%hc).EQ.1)call stop('DMBLcoupling_init: coupling primitive failed check')    

    this%nbath=0
    if(present(file))then
       read(unit,*)this%nbath
    else
       do while(this%nbath.LT.1)
          write(*,*)'Please enter the number of independent baths.'
          read(*,*)this%nbath
          if(this%nbath.LT.1)write(*,*)'Number of baths must be a positive integer, Try again!'
       end do
    end if

    if(associated(this%map))nullify(this%map)
    allocate(this%map(this%hc%hb%ndof),stat=ierr)
    if(ierr.NE.0)call stop('DMBLcoupling_init: failed to allocate memory for map.')
    this%map=0
    if(present(file))then
       read(unit,*)(this%map(idof),idof=1,this%hc%hb%ndof)
    else
       pass=.false.
       if(mod(this%hc%hb%ndof,this%nbath).EQ.0)pass=.true.
       if(pass)then
          write(*,*)'Use Default dof to nbath mapping? [T/F] (Divide dofs equally among nbaths.)'
          read(*,*)usedefault
          if(usedefault)then
             bdof=this%hc%hb%ndof/this%nbath
             do ibath=1,this%nbath
                do idof=1,bdof
                   this%map(idof+(ibath-1)*bdof)=ibath
                end do
             end do
          else
             do idof=1,this%hc%hb%ndof
                ibath=0
                do while(ibath.LE.0.or.ibath.GT.this%nbath)
                   write(*,*)'Associate dof '//trim(int2str(idof))//' with a bath [1:'//trim(int2str(this%nbath))//'].'
                   read(*,*)ibath
                   if(ibath.LE.0.or.ibath.GT.this%nbath)&
                        write(*,*)'Not a valid entry. Please enter an integer between(including) 1 and '//trim(int2str(this%nbath))//'. Try again.'
                end do
                this%map(idof)=ibath
             end do
          end if
       end if
    end if
    
    if(associated(this%bdof))nullify(this%bdof)
    allocate(this%bdof(this%nbath),stat=ierr)
    if(ierr.NE.0)call stop('DMBLcoupling_init: failed to allocate memory for bdof.')
    this%bdof=0
    do ibath=1,this%nbath
       do idof=1,this%hc%hb%ndof
          if(this%map(idof).EQ.ibath)this%bdof(ibath)=this%bdof(ibath)+1
       end do
    end do

    if(associated(this%lambda))nullify(this%lambda)
    allocate(this%lambda(this%nbath),stat=ierr)
    if(ierr.NE.0)call stop('DMBLcoupling_init: failed to allocate memory for lambda.')
    this%lambda=0._double
    if(present(file))then
       read(unit,*)(this%lambda(ibath),ibath=1,this%nbath)
    else
       do ibath=1,this%nbath
          write(*,*)'Please enter reorganization energy for bath '//trim(int2str(ibath))//'.'
          read(*,*)this%lambda(ibath)
       end do
    end if

!!$    if(associated(this%wc))nullify(this%wc)
!!$    allocate(this%wc(this%nbath),stat=ierr)
!!$    if(ierr.NE.0)call stop('DMBLcoupling_init: failed to allocate memory for wc.')
!!$    this%wc=0._double
!!$    if(present(file))then
!!$       read(unit,*)(this%wc(ibath),ibath=1,this%nbath)
!!$    else
!!$       do ibath=1,this%nbath
!!$          write(*,*)'Please enter Debey cutoff energy for bath '//trim(int2str(ibath))//'.'
!!$          read(*,*)this%wc(ibath)
!!$       end do
!!$    end if
    
    if(associated(this%g))nullify(this%g)
    allocate(this%g(this%hc%hb%ndof,this%nbath),stat=ierr)
    if(ierr.NE.0)call Stop('DMBLcoupling_init: mode-bath interaction memory allocation error.')
    
    if(associated(this%f))nullify(this%f)
    allocate(this%f(this%hc%hs%nstate,this%nbath),stat=ierr)
    if(ierr.NE.0)call Stop('DMBLcoupling_init: state-bath interaction memory allocation error.')
    this%f=0.0_double
    if(present(file))then
       read(unit,*)((this%f(istate,ibath),ibath=1,this%nbath),istate=1,this%hc%hs%nstate)
    else
       usedefault=.true.
       write(*,*)'Use default state-bath interaction? [T/F] (Bath n associated only to population state |n><n|. Nbath must = nstate!)'
       read(*,*)usedefault
       if(usedefault.and.this%nbath.NE.this%hc%hs%nstate)then
          call warn('DMBLcoupling_init: nstate must equal nbath.','Switching to manual state to bath mapping.')
          usedefault=.false.
       end if
       if(usedefault)then
          do istate=1,this%hc%hs%nstate
             ibath=istate
             this%f(istate,ibath)=1.0_double
          end do
       else
          write(*,*)'Enter state-bath interaction terms'
          do istate=1,this%hc%hs%nstate
             do ibath=1,this%nbath
                write(*,*)'Enter state '//trim(int2str(istate))//' to bath '//trim(int2str(ibath))//' relative coupling strength.'
                read(*,*)this%f(istate,ibath)
             end do
          end do
       end if
    end if
    !!normalize state-bath interaction terms
    !do ibath=1,this%nbath
    !   norm=sum(this%f(:,ibath)**2)
    !   if(norm.NE.0.)this%f(:,ibath)=this%f(:,ibath)/sqrt(norm)
    !end do

    if(associated(this%X))nullify(this%X)
    allocate(this%X(this%hc%hb%ndof),stat=ierr)
    if(ierr.NE.0)call Stop('DMBLcoupling_init: Normal Mode Coordinate memory allocation error.')
    this%X=0.0_double

    this%initialized=.true.
    call resample(this)
    if(check(this).EQ.1)call stop('DMBLcoupling_init: failed final object check')

  end subroutine DMBLcoupling_init
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine DMBLcoupling_kill(this)

    type(DMBLcoupling),intent(inout)::this

    call Note('Begin DMBLcoupling_kill.')
    call kill(this%hc)
    if(associated(this%lambda))nullify(this%lambda)
    !if(associated(this%wc))nullify(this%wc)
    if(associated(this%map))nullify(this%map)
    if(associated(this%g))nullify(this%g)
    if(associated(this%f))nullify(this%f)
    if(associated(this%X))nullify(this%X)
    this%initialized=.false.
  end subroutine DMBLcoupling_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine DMBLcoupling_update(this)
    type(DMBLcoupling),intent(inout)::this
    integer(long)::istate,idof,jdof,ibath,icount
    real(double)::wr,counterterm,tol

    call Note('Begin DMBLcoupling_update.')

    !call update(this%hc)

    this%X=0._double
    do idof=1,this%hc%hb%ndof
       this%X=sum(this%hc%hb%EigenVec(:,idof)*this%hc%hb%Q)
    end do

    this%hc%V=0.0_double
    this%hc%dV=0.0_double
    do ibath=1,this%nbath
       do idof=1,this%hc%hb%ndof
          if(this%map(idof).EQ.ibath)then
             counterterm=0._double
             if(this%hc%hb%canonical)&
                  counterterm=.5_double*this%hc%hb%rmass(idof)&
                  *(this%g(idof,ibath)/this%hc%hb%mode(idof))**2
             do istate=1,this%hc%hs%nstate
                this%hc%V(istate,istate)=this%hc%V(istate,istate)&
                     -this%X(idof)*this%g(idof,ibath)*this%f(istate,ibath)&
                     +2_double*this%lambda(ibath)*this%f(istate,ibath)!**2counterterm
                     !+this%f(istate,ibath)**2*counterterm
                this%hc%dV(idof,istate,istate)=&
                     -this%g(idof,ibath)*this%f(istate,ibath)
             end do
          end if
       end do
    end do
    
  end subroutine DMBLcoupling_update
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine DMBLcoupling_resample(this)
    type(DMBLcoupling),intent(inout)::this
    integer(long)::istate,idof,jdof,ibath,bdof
    real(double)::wc,wo,wm,wj,wr,Jw

    call Note('Begin DMBLcoupling_resample.')

    !call resample(this%hc)

    !resample normal mode oscillator strengths
    this%g=0._double
    do ibath=1,this%nbath
       wr=hbar*this%lambda(ibath)/real(this%bdof(ibath))
       do idof=1,this%hc%hb%ndof
          if(this%map(idof).EQ.ibath)this%g(idof,ibath)=&
               this%hc%hb%mode(idof)*sqrt(wr/this%hc%hb%rmass(idof))      
       end do
    end do

    call update(this)
    
  end subroutine DMBLcoupling_resample
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine DMBLcoupling_save(this,file)
    type(DMBLcoupling),intent(in)::this
    character*(*),intent(in)::file

    integer(short)::ierr
    integer(long)::unit
    logical::usedefault

    integer(long)::istate,idof,ibath

    call Note('Begin DMBLcoupling_save.')
    call Note('input file= '//file)
    if(check(this).EQ.1)call stop('DMBLcoupling_save')
    
    unit=newunit()
    open(unit,file=file)
    write(unit,*)'DMBLcoupling'
    write(unit,*)quote(file//'.hc')
    call save(this%hc,file//'.hc')
    write(unit,*)this%nbath
    write(unit,*)(this%map(idof),idof=1,this%hc%hb%ndof)
    write(unit,*)(this%lambda(ibath),ibath=1,this%nbath)
    !write(unit,*)(this%wc(ibath),ibath=1,this%nbath)
    write(unit,*)((this%f(istate,ibath),ibath=1,this%nbath),istate=1,this%hc%hs%nstate)
    close(unit)      

  end subroutine DMBLcoupling_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine DMBLcoupling_display(this,msg)
    type(DMBLcoupling),intent(in)::this
    character*(*),intent(in),optional::msg
    !character(len=comment)::mapstring
    integer(long)::idof,ibath,istate
    logical::mask(this%hc%hb%ndof)

    call Note('Begin DMBLcoupling_display.')
    if(check(this).EQ.0)then
       write(*,*)'################## DMBLcoupling ################'
       write(*,*)'nbath=',this%nbath
       write(*,*)'re-organization energy=',this%lambda
       !write(*,*)'Debey cutoff energy=',this%wc
       write(*,*)
       write(*,*)'mode to bath mapping'
       do ibath=1,this%nbath
          write(*,*)
          write(*,*)'bath '//trim(int2str(ibath))
          write(*,'(A7)',advance='no')'modes: '
          mask=.false.
          do idof=1,this%hc%hb%ndof
             if(this%map(idof).EQ.ibath)then
                mask(idof)=.true.
                write(*,'(A3,1X)',advance='no')trim(int2str(idof))
             end if
          end do
          write(*,*)
          call display(this%hc%hb%mode,this%g(:,ibath),mask=mask)
       end do
       write(*,*)
       write(*,*)'state(col)-bath(row) interaction terms'
       call display(transpose(this%f))
       write(*,*)
       write(*,*)'mode(col)-bath(row) interaction terms'
       call display(transpose(this%g))
       call display(this%hc,'DMBLcoupling primitive coupling term')
       write(*,*)'############### end DMBLcoupling ################'
    else
       call warn('DMBLcoupling_display: coupling object failed check.','displaying nothing.')
    end if
    
  end subroutine DMBLcoupling_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short) function DMBLcoupling_check(this)
    type(DMBLcoupling),intent(in)::this
    
    integer(long)::istate,idof,ibath
    integer(short)::ierr
    
    call Note('Checking DMBLcoupling.')

    DMBLcoupling_check=0
    
    !logical::initialized=.false.
    if(.not.this%initialized)then
       call Warn('DMBLcoupling_check: coupling primitive not initialized.')
       DMBLcoupling_check=1
       return
    end if
    !type(hc)::hc
    if(check(this%hc).NE.0)then
       call Warn('DMBLcoupling_check: quantum primitive failed check')
       DMBLcoupling_check=1
       return
    end if
    
    !integer(long)::nbath
    if(this%nbath.NE.this%nbath)then
       call Warn('DMBLcoupling_check: nbath not a number.')
       DMBLcoupling_check=1
       return
    end if
    if(abs(this%nbath).GE.huge(this%nbath))then
       call Warn('DMBLcoupling_check: nbath is too big.')
       DMBLcoupling_check=1
       return
    end if
    if(this%nbath.LE.0)then
       call Warn('DMBLcoupling_check: nbath not a positive integer.')
       DMBLcoupling_check=1
       return
    end if
    
    !real(double),dimension(:),pointer::lambda
    if(.not.associated(this%lambda))then
       call Warn('DMBLcoupling_check: lambda memory not allocated')
       DMBLcoupling_check=1
       return
    end if
    if(size(this%lambda).NE.this%nbath)then
       call Warn('DMBLcoupling_check: size of lambda not equal number of baths')
       DMBLcoupling_check=1
       return
    end if
    ierr=0
    do ibath=1,this%nbath
       if(this%lambda(ibath).NE.this%lambda(ibath))ierr=1 
       if(abs(this%lambda(ibath)).GE.huge(this%lambda(ibath)))ierr=1 
    end do
    if(ierr.EQ.1)then
       call Warn('DMBLcoupling_check: lambda has bad values.')
       DMBLcoupling_check=1
       return
    end if
    
    !real(double),dimension(:),pointer::X
    if(.not.associated(this%X))then
       call Warn('DMBLcoupling_check: X memory not allocated')
       DMBLcoupling_check=1
       return
    end if
    if(size(this%X).NE.this%hc%hb%ndof)then
       call Warn('DMBLcoupling_check: size of X not equal number of classcial dofs')
       DMBLcoupling_check=1
       return
    end if
    ierr=0
    do idof=1,this%hc%hb%ndof
       if(this%X(idof).NE.this%X(idof))ierr=1 
       if(abs(this%X(idof)).GE.huge(this%X(idof)))ierr=1 
    end do
    if(ierr.EQ.1)then
       call Warn('DMBLcoupling_check: X has bad values.')
       DMBLcoupling_check=1
       return
    end if
    
    !integer(long),dimension(:),pointer::map
    if(.not.associated(this%map))then
       call Warn('DMBLcoupling_check: map memory not allocated')
       DMBLcoupling_check=1
       return
    end if
    if(size(this%map).NE.this%hc%hb%ndof)then
       call Warn('DMBLcoupling_check: size of map not equal number of classcial dofs')
       DMBLcoupling_check=1
       return
    end if
    ierr=0
    do idof=1,this%hc%hb%ndof
       if(this%map(idof).NE.this%map(idof))ierr=1 
       if(abs(this%map(idof)).GE.huge(this%map(idof)))ierr=1 
       if(this%map(idof).GT.this%nbath)ierr=1 
       if(this%map(idof).GT.this%nbath)ierr=1 
       if(this%map(idof).LT.1)ierr=1 
    end do
    if(ierr.EQ.1)then
       call Warn('DMBLcoupling_check: map has bad values.')
       DMBLcoupling_check=1
       return
    end if
    
    !integer(long),dimension(:),pointer::bdof
    if(.not.associated(this%bdof))then
       call Warn('DMBLcoupling_check: bdof memory not allocated')
       DMBLcoupling_check=1
       return
    end if
    if(size(this%bdof).NE.this%nbath)then
       call Warn('DMBLcoupling_check: size of bdof array  not equal number of baths')
       DMBLcoupling_check=1
       return
    end if
    if(any(this%bdof.NE.this%bdof))then
       call Warn('DMBLcoupling_check: bdof has NAN values.')
       DMBLcoupling_check=1
       return
    end if
    if(any(abs(this%bdof).GE.huge(this%bdof)))then
       call Warn('DMBLcoupling_check: bdof has huge values.')
       DMBLcoupling_check=1
       return
    end if
    if(any(this%bdof.GT.this%hc%hb%ndof))then
       call Warn('DMBLcoupling_check: bdof is greater than total dofs.') 
       DMBLcoupling_check=1
       return
    end if
    if(any(this%bdof.LE.0))then
       call Warn('DMBLcoupling_check: bdof must be positive.')
       DMBLcoupling_check=1
       return
    end if
    
    !real(double),dimension(:,:),pointer::g
    if(.not.associated(this%g))then
       call Warn('DMBLcoupling_check: mode-bath memory not allocated')
       DMBLcoupling_check=1
       return
    end if
    if(size(this%g).NE.this%hc%hb%ndof*this%nbath)then
       call Warn('DMBLcoupling_check: size of mode-bath interactions not equal number of classcial dofs times nbath')
       DMBLcoupling_check=1
       return
    end if
    ierr=0
    do idof=1,this%hc%hb%ndof
       do ibath=1,this%nbath
          if(this%g(idof,ibath).NE.this%g(idof,ibath))ierr=1 
          if(abs(this%g(idof,ibath)).GE.huge(this%g(idof,ibath)))ierr=1
       end do
    end do
    if(ierr.EQ.1)then
       call Warn('DMBLcoupling_check: mode-bath interaction terms have bad values.')
       DMBLcoupling_check=1
       return
    end if
    !real(double),dimension(:,:),pointer::f
    if(.not.associated(this%f))then
       call Warn('DMBLcoupling_check: state-bath memory not allocated')
       DMBLcoupling_check=1
       return
    end if
    if(size(this%f).NE.this%hc%hs%nstate*this%nbath)then
       call Warn('DMBLcoupling_check: size of state-bath interactions not equal number of quatnum states times nbath')
       DMBLcoupling_check=1
       return
    end if
    ierr=0
    do istate=1,this%hc%hs%nstate
       do ibath=1,this%nbath
          if(this%f(istate,ibath).NE.this%f(istate,ibath))ierr=1 
          if(abs(this%f(istate,ibath)).GE.huge(this%f(istate,ibath)))ierr=1
       end do
    end do
    if(ierr.EQ.1)then
       call Warn('DMBLcoupling_check: stae-bath interaction terms have bad values.')
       DMBLcoupling_check=1
       return
    end if
    
  end function DMBLcoupling_check
  
end module DMBLcoupling_class

