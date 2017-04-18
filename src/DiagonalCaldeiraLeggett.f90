!-----------------------------------------------------------------------------
!> \brief
!! Diagonal Bi-linear Coupling.
!! \details
!! This coupling scheme follows the Caldeira-Leggett model, coupling the normal modes
!! of the classical subsystem linearly to the diagonal states of the quantum subsystem.
!! Assumes Hamiltonian of the form
!!        \f[  H=h_s+h_b+h_c+ct \f]
!! where the quantum subsystem Hamiltonian matrix elements  \f$ \langle\psi_i|hs|\psi_j\rangle \f$
!! are the matrix elements of the diabats attribute of the hs type.
!! The classical subsystem potential, here \f$ hb=V(Q) \f$
!! ,is written in terms of the normal modes of the classical subsystem \f$Q\f$.
!! These normal mode coordinates can be found in the W attribute of the hb type.
!! Each state is coupled to one of \f$ N \f$ normal modes.
!! Here we assume that the diabatic basis states \f$ \psi \f$
!! are eigenstates of some quantum subsystem operator \f$ \hat q \f$ 
!! such that \f[ \hat q |\psi_j\rangle =f_j|\psi_j\rangle \f]
!! and the eigenvalues \f$ f_j \f$ are supplied by the user.
!! The bilinear coupling term is written
!!   \f[ h_c=-\sum_n^N \sum_k^{K_n} g_{k,n} Q_{k} f_n |\psi_n\rangle\langle\psi_n|\f]
!! with counter term
!!   \f[ ct=\sum_n^N \sum_k^{K_n} f_n^2 \frac{g_{k,n}^2}{2m_{k,n}\omega_{k,n}^2}|\psi_n\rangle\langle\psi_n|. \f]
!! The spectral density associated with state \f$ n \f$ is defined
!! \f[ J_n(\omega)=\frac{\pi}{\hbar}\sum_k^{K_n}\frac{g_{k,n}^2}{2m_{k,n}\omega_{k,n}}\delta_{\omega-\omega_{k,n}}.\f]
!! The re-organization energy for each state is defined
!!   \f[ \lambda_n=\frac{1}{\pi} \int_0^\infty \frac{J_n(\omega)}{\omega} d\omega \f]
!! such that,
!!   \f[ g_{k,n}^2=\frac{2m_{k,n} \lambda_n m_{k,n} \omega_{k,n}^2}{K_n}\f]
!! and 
!!   \f[ ct=2 \hbar \sum_n^N \lambda_n f_n^2|\psi_n\rangle\langle\psi_n| .\f]
!! \authors
!! Daniel Montemayor
!!
!! \date
!! Aug 2011, Sept 2011, Oct 2012
!<-----------------------------------------------------------------------------
!Change Log:
!=== v0.1 Oct 2012 ===
! + introduced Wbasis logical variable (default false)
! + Assumes hb%W is mass weighted transform of Q coordinates
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

module DiagonalCaldeiraLeggett_class
  use type_kinds
  use atomicunits
  use string
  use ErrorLog
  use filemanager
  use hc_class
  use quantum_class
  use classical_class
  use textgraphs
  use outputdisplay
  implicit none
  private

  public::DiagonalCaldeiraLeggett
  public::new,kill,update,resample,display,save,check

  !> Derived quantum-classical coupling term. 
  type DiagonalCaldeiraLeggett
     !> True if coupling term has been properly initialized
     logical::initialized=.false.
     !> Inherited coupling primitive type
     type(hc)::hc
     !> \brief Reorganization energy
     !! \details nstate long array containing the vaules of
     !! the reorganization energy of each states dissipative bath.
     real(double),dimension(:),pointer::lambda
     !> Number of dissipative modes associated to each state. 
     integer(long),dimension(:),pointer::bdof
     !> \brief Coupling coefficient coupling quantum state to classical mode.
     !! \details The matrix g couples the
     !! native diabatic basis states of the quantum subsystem with the normal
     !! modes of the classical subsystem and has dimensions of ndofs by nstate. 
     real(double),dimension(:,:),pointer::g
     !> \brief Eigen values of the \f$ \hat q \f$ operator
     !! \details The native diabatic basis states are assumed eigen states
     !! of the \f$ \hat q \f$ operator with eigen values stored in the array f of dimension nstate.  
     real(double),dimension(:),pointer::f
     !! \brief Logically maps the classical modes to a quantum state.
     !! \details Has dimensions ndof by nstate.
     logical,dimension(:,:),pointer::gmap
     !! True if the classical subsystem Potential is assummed is diagonal, i.e. non-interacting particles.
     logical::Wbasis
  end type DiagonalCaldeiraLeggett

  !> Derived coupling type creator
  interface new
     module procedure DiagonalCaldeiraLeggett_init
  end interface

  !> Derived coupling type destroyer.
  interface kill
     module procedure DiagonalCaldeiraLeggett_kill
  end interface

  !> Recalculates the coupling primitive's attributes.
  interface update
     module procedure DiagonalCaldeiraLeggett_update
  end interface

  !> Reinitializes the coupling primitive.
  interface resample
     module procedure DiagonalCaldeiraLeggett_resample
  end interface

  !> Displays the current state of the derived coupling type.
  interface display
     module procedure DiagonalCaldeiraLeggett_display
  end interface

  !> Saves the current state of the derived coupling type to file.
  interface save
     module procedure DiagonalCaldeiraLeggett_save
  end interface

  !> Checks the derived coupling type.
  interface check
     module procedure DiagonalCaldeiraLeggett_check
  end interface

contains
  !======================================================================
  !> \brief
  !! Creates and initializes the DiagonalCaldeiraLegget coupling type.
  !! \param[inout] this is the derived coupling type to be initialized.
  !! \param[inout] qs is a general quantum subsystem to be coupled to a classical subsystem.
  !! \param[inout] cs is a general classical subsysterm to be coupled to the quantum subsystem.
  !! \param[in] file is a string containing the name of a previously saved DiagonalCaldeiraLeggett
  !!            used as an input file to initialize THIS.
  !! \remark If no input file provided the user must manually initialize THIS using stout.
  !!<====================================================================
  subroutine DiagonalCaldeiraLeggett_init(this,qs,cs,file)
    type(DiagonalCaldeiraLeggett),intent(inout)::this
    type(quantum),intent(inout),target::qs
    type(classical),intent(inout),target::cs
    character*(*),intent(in),optional::file      
    character(len=path)::filename
    character(len=title)::filetype
    integer(short)::ierr
    integer(long)::unit
    logical::usedefault,pass

    integer(long)::ndof,nstate,bdof
    integer(long)::idof,istate
    real(double)::norm


    call Note('Begin DiagonalCaldeiraLeggett_init.')
    filename='unknown'
    filetype='unknown'

    if(check(qs).EQ.1)call stop('DiagonalCaldeiraLeggett_init, qs failed check')
    if(check(cs).EQ.1)call stop('DiagonalCaldeiraLeggett_init, cs failed check')

    if(present(file))call Note('input file= '//trim(file))

    if(present(file))then
       if(check(trim(file)).EQ.1)&
            call stop('DiagonalCaldeiraLeggett_init, cannot find input file')
       
       unit=newunit()
       open(unit,file=file)
       read(unit,*)filetype
       filetype=adjustl(filetype)
       if(trim(filetype).EQ.'DiagonalCaldeiraLeggett')then
          read(unit,*)filename
          filename=adjustl(filename)
          call new(this%hc,qs,cs,trim(filename))
       else
          call Stop('DiagonalCaldeiraLeggett_init: input file not valid type.','Check input file format for DiagonalCaldeiraLeggett object.')
       end if
    else
       call new(this%hc,qs,cs)
    end if
    if(check(this%hc).EQ.1)call stop('DiagonalCaldeiraLeggett_init: coupling primitive failed check')    

    ndof=this%hc%hb%ndof
    nstate=this%hc%hs%nstate

    if(associated(this%gmap))nullify(this%gmap)
    allocate(this%gmap(ndof,nstate),stat=ierr)
    if(ierr.NE.0)call stop('DiagonalCaldeiraLeggett_init: failed to allocate memory for mode-state array gmap.')
    this%gmap=.false.
    if(present(file))then
       read(unit,*)((this%gmap(idof,istate),istate=1,nstate),idof=1,ndof)
    else
       pass=.false.
       if(mod(ndof,nstate).EQ.0)pass=.true.
       usedefault=.false.
       if(pass)then
          write(*,*)'Use default mode to state mapping? [T/F] (Divide modes equally among quantum states.)'
          read(*,*)usedefault
       end if
       if(usedefault)then
          bdof=ndof/nstate
          do istate=1,nstate
             do idof=1,bdof
                this%gmap(idof+(istate-1)*bdof,istate)=.true.
             end do
          end do
       else
          write(*,*)'OK, would you like all states coupled to all modes? T/F'
          read(*,*)pass
          if(pass)then
             this%gmap=.true.
          else
             do idof=1,ndof
                do istate=1,nstate
                   write(*,*)'Associate mode '//trim(int2str(idof))//' with state '//trim(int2str(istate))//'? T/F'
                   read(*,*)this%gmap(idof,istate)
                end do
             end do
          end if
       end if
    end if
    if(associated(this%bdof))nullify(this%bdof)
    allocate(this%bdof(nstate),stat=ierr)
    if(ierr.NE.0)call stop('DiagonalCaldeiraLeggett_init: failed to allocate memory for bdof.')
    this%bdof=0
    do istate=1,nstate
       do idof=1,ndof
          if(this%gmap(idof,istate))this%bdof(istate)=this%bdof(istate)+1
       end do
    end do

    if(associated(this%lambda))nullify(this%lambda)
    allocate(this%lambda(nstate),stat=ierr)
    if(ierr.NE.0)call stop('DiagonalCaldeiraLeggett_init: failed to allocate memory for lambda.')
    this%lambda=0._double
    if(present(file))then
       read(unit,*)(this%lambda(istate),istate=1,nstate)
    else
       do istate=1,nstate
          write(*,*)'Please enter reorganization energy for state '//trim(int2str(istate))//'.'
          read(*,*)this%lambda(istate)
       end do
    end if
    
    
    if(associated(this%g))nullify(this%g)
    allocate(this%g(this%hc%hb%ndof,nstate),stat=ierr)
    if(ierr.NE.0)call Stop('DiagonalCaldeiraLeggett_init: mode-state coefficient memory allocation error.')
    
    if(associated(this%f))nullify(this%f)
    allocate(this%f(this%hc%hs%nstate),stat=ierr)
    if(ierr.NE.0)call Stop('DiagonalCaldeiraLeggett_init: qhat operator eigenvalue allocation error.')
    this%f=0.0_double
    if(present(file))then
       read(unit,*)(this%f(istate),istate=1,nstate)
    else
       usedefault=.true.
       write(*,*)'Use default qhat operator eigenvalues? [all = 1.0] T/F '
       read(*,*)usedefault
       if(usedefault)then
          this%f=1.0_double
       else
          write(*,*)'Enter qhat operator eignvalues.'
          do istate=1,nstate
             write(*,*)'Enter eigenvalue '//trim(int2str(istate))
             read(*,*)this%f(istate)
          end do
       end if
    end if

    this%Wbasis=.false.
    this%initialized=.true.
    call resample(this)
    if(check(this).EQ.1)call stop('DiagonalCaldeiraLeggett_init: failed final object check')

  end subroutine DiagonalCaldeiraLeggett_init
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine DiagonalCaldeiraLeggett_kill(this)

    type(DiagonalCaldeiraLeggett),intent(inout)::this

    call Note('Begin DiagonalCaldeiraLeggett_kill.')
    call kill(this%hc)
    if(associated(this%lambda))nullify(this%lambda)
    if(associated(this%gmap))nullify(this%gmap)
    if(associated(this%bdof))nullify(this%bdof)
    if(associated(this%g))nullify(this%g)
    if(associated(this%f))nullify(this%f)
    this%initialized=.false.
  end subroutine DiagonalCaldeiraLeggett_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine DiagonalCaldeiraLeggett_update(this)
    type(DiagonalCaldeiraLeggett),intent(inout)::this

    integer(long)::nstate,ndof
    integer(long)::istate,idof
    real(double)::counterterm

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof

    call Note('Begin DiagonalCaldeiraLeggett_update.')

    !call update(this%hc)

    this%hc%V=0.0_double
    this%hc%dV=0.0_double
    do istate=1,nstate
       do idof=1,ndof
          if(this%gmap(idof,istate))then
             counterterm=0._double
             
             if(this%hc%hb%canonical)&
                  counterterm=.5_double*(this%g(idof,istate)*this%f(istate)&
                  /this%hc%hb%mode(idof))**2
             
             this%hc%V(istate,istate)=this%hc%V(istate,istate)&
                  -this%hc%hb%W(idof)*this%g(idof,istate)*this%f(istate)&
                  +counterterm
             this%hc%dV(idof,istate,istate)=&
                  -this%g(idof,istate)*this%f(istate)
          end if
       end do
       !change dV into Q basis
       if(.not.this%Wbasis)this%hc%dV(:,istate,istate)=matmul(transpose(this%hc%hb%EigenVec),this%hc%dV(:,istate,istate))
    end do
    
  end subroutine DiagonalCaldeiraLeggett_update
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine DiagonalCaldeiraLeggett_resample(this)
    type(DiagonalCaldeiraLeggett),intent(inout)::this
    integer(long)::nstate,ndof
    integer(long)::istate,idof
    call Note('Begin DiagonalCaldeiraLeggett_resample.')

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof

    !call resample(this%hc)

    !resample normal mode oscillator strengths
    this%g=0._double
    do istate=1,nstate
       do idof=1,ndof
          if(this%gmap(idof,istate))this%g(idof,istate)=this%hc%hb%mode(idof)&
               *sqrt(2._double*this%lambda(istate)/this%bdof(istate))   
       end do
    end do
    if(all(this%hc%hb%EigenVec.eq.iden(ndof)))this%Wbasis=.true.
    if(check(this).EQ.1)call stop('DiagonalCaldeiraLeggett_resample: object failed exit check.')
    call update(this)
    
  end subroutine DiagonalCaldeiraLeggett_resample
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine DiagonalCaldeiraLeggett_save(this,file)
    type(DiagonalCaldeiraLeggett),intent(in)::this
    character*(*),intent(in)::file

    integer(short)::ierr
    integer(long)::unit
    logical::usedefault

    integer(long)::nstate,ndof
    integer(long)::istate,idof

    call Note('Begin DiagonalCaldeiraLeggett_save.')
    call Note('input file= '//file)
    if(check(this).EQ.1)call stop('DiagonalCaldeiraLeggett_save')

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof
    
    unit=newunit()
    open(unit,file=file)
    write(unit,*)'DiagonalCaldeiraLeggett'
    write(unit,*)quote(file//'.hc')
    call save(this%hc,file//'.hc')
    write(unit,*)((this%gmap(idof,istate),istate=1,nstate),idof=1,ndof)
    write(unit,*)(this%lambda(istate),istate=1,nstate)
    write(unit,*)(this%f(istate),istate=1,this%hc%hs%nstate)
    close(unit)      

  end subroutine DiagonalCaldeiraLeggett_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine DiagonalCaldeiraLeggett_display(this,msg)
    type(DiagonalCaldeiraLeggett),intent(in)::this
    character*(*),intent(in),optional::msg

    integer(long)::ndof,nstate
    integer(long)::idof,istate
    logical::mask(this%hc%hb%ndof)

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof

    call Note('Begin DiagonalCaldeiraLeggett_display.')
    if(check(this).NE.0)then
       call warn('DiagonalCaldeiraLeggett_display: failed check','displaying nothing.')
       return
    end if
    write(Dunit,*)'____________________________________________________'
    write(Dunit,*)'-----------  DiagonalCaldeiraLeggett  --------------'
    write(Dunit,*)'re-organization energy by state=',this%lambda
    call display(this%lambda)
    write(Dunit,*)'mode to state mapping'
    do istate=1,nstate
       write(Dunit,*)
       write(Dunit,*)'state '//trim(int2str(istate))
       write(Dunit,'(A7)',advance='no')'modes: '
       mask=.false.
       do idof=1,ndof
          if(this%gmap(idof,istate))then
             mask(idof)=.true.
             write(Dunit,'(A3,1X)',advance='no')trim(int2str(idof))
          end if
       end do
       call display(this%hc%hb%mode,this%g(:,istate),mask=mask)
       write(Dunit,*)
    end do
    write(Dunit,*)'qhat operator eigenvalues:',this%f
    call display(this%f)
    call display(this%hc,'DiagonalCaldeiraLeggett primitive')
    write(Dunit,*)'===================================================='
    
  end subroutine DiagonalCaldeiraLeggett_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short) function DiagonalCaldeiraLeggett_check(this)
    type(DiagonalCaldeiraLeggett),intent(in)::this
    
    integer(long)::nstate,ndof
    integer(long)::istate,idof
    integer(short)::ierr
    
    call Note('Checking DiagonalCaldeiraLeggett.')

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof

    DiagonalCaldeiraLeggett_check=0
    
    !logical::initialized=.false.
    if(.not.this%initialized)then
       call Warn('DiagonalCaldeiraLeggett_check: coupling primitive not initialized.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    !type(hc)::hc
    if(check(this%hc).NE.0)then
       call Warn('DiagonalCaldeiraLeggett_check: quantum primitive failed check')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    
    !real(double),dimension(:),pointer::lambda
    if(.not.associated(this%lambda))then
       call Warn('DiagonalCaldeiraLeggett_check: lambda memory not allocated')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(size(this%lambda).NE.nstate)then
       call Warn('DiagonalCaldeiraLeggett_check: size of lambda not equal number of states')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(this%lambda.NE.this%lambda))then
       call Warn('DiagonalCaldeiraLeggett_check: lambda has NAN values.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(abs(this%lambda).GT.huge(this%lambda)))then
       call Warn('DiagonalCaldeiraLeggett_check: lambda has Huge values.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(abs(this%lambda).LT.epsilon(this%lambda)))then
       call Warn('DiagonalCaldeiraLeggett_check: lambda has tiny values.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    
    !integer(long),dimension(:,:),pointer::gmap
    if(.not.associated(this%gmap))then
       call Warn('DiagonalCaldeiraLeggett_check: gmap memory not allocated')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(size(this%gmap,1).NE.ndof)then
       call Warn('DiagonalCaldeiraLeggett_check: dimension 1 of gmap array not equal number of classcial dofs.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(size(this%gmap,2).NE.nstate)then
       call Warn('DiagonalCaldeiraLeggett_check: dimension 2 of gmap array not equal number of states.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(this%gmap.NE.this%gmap))then
       call Warn('DiagonalCaldeiraLeggett_check: gmap array has bad values.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if

    
    !integer(long),dimension(:),pointer::bdof
    if(.not.associated(this%bdof))then
       call Warn('DiagonalCaldeiraLeggett_check: bdof memory not allocated')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(size(this%bdof).NE.nstate)then
       call Warn('DiagonalCaldeiraLeggett_check: size of bdof array  not equal number of states')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(this%bdof.NE.this%bdof))then
       call Warn('DiagonalCaldeiraLeggett_check: bdof has NAN values.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(abs(this%bdof).GE.huge(this%bdof)))then
       call Warn('DiagonalCaldeiraLeggett_check: bdof has huge values.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(this%bdof.GT.this%hc%hb%ndof))then
       call Warn('DiagonalCaldeiraLeggett_check: bdof is greater than total dofs.') 
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(this%bdof.LE.0))then
       call Warn('DiagonalCaldeiraLeggett_check: bdof must be positive.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    
    !real(double),dimension(:,:),pointer::g
    if(.not.associated(this%g))then
       call Warn('DiagonalCaldeiraLeggett_check: mode-state memory not allocated')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(size(this%g,1).NE.ndof)then
       call Warn('DiagonalCaldeiraLeggett_check: size of mode-state coefficent array  in dimension 1 does not equal number of classcial dofs.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(size(this%g,2).NE.nstate)then
       call Warn('DiagonalCaldeiraLeggett_check: size of mode-state coefficent array in dimension 2 does not equal number of quatnum states.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(this%g.NE.this%g))then
       call Warn('DiagonalCaldeiraLeggett_check: mode-state coefficent array has NAN values.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(abs(this%g).GT.huge(this%g)))then
       call Warn('DiagonalCaldeiraLeggett_check: mode-state coefficent array has Huge values.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(this%g.LT.0._double))then
       call Warn('DiagonalCaldeiraLeggett_check: mode-state coefficent array has negative values.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if

    !real(double),dimension(:),pointer::f
    if(.not.associated(this%f))then
       call Warn('DiagonalCaldeiraLeggett_check: qhat eigenvalue array memory not allocated')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(size(this%f).NE.nstate)then
       call Warn('DiagonalCaldeiraLeggett_check: size of qhat eigenvalue array not equal number of quatnum states')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(this%f.NE.this%f))then
       call Warn('DiagonalCaldeiraLeggett_check: qhat eigenvalue array has NAN values.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    if(any(abs(this%f).GT.huge(this%f)))then
       call Warn('DiagonalCaldeiraLeggett_check: qhat eigenvalue array has Huge values.')
       DiagonalCaldeiraLeggett_check=1
       return
    end if
    
  end function DiagonalCaldeiraLeggett_check
  
end module DiagonalCaldeiraLeggett_class

