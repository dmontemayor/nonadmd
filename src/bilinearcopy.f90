!-----------------------------------------------------------------------------
!> \brief
!! Diagonal Bi-linear system-bath coupling.
!! \details
!! This coupling scheme follows the Caldeira-Leggett model, coupling the normal modes
!! of the classical subsystem linearly to the diagonal states of the quantum subsystem.
!! Assumes Hamiltonian of the form
!!        \f[  H=h_s+h_b+h_c+ct \f]
!! where the quantum subsystem Hamiltonian matrix elements  \f$ \langle\psi_i|hs|\psi_j\rangle \f$
!! are the matrix elements of \f$qs.hs.diabats\f$.
!! The classical subsystem \f$ hb \f$ potential energy
!! ,is written in terms of the normal modes of the classical subsystem \f$cs.hb.W\f$.
!! Each one of the \f$ N\f$ states is coupled to one or more of the \f$ K \f$ normal modes.
!! Here we assume that the diabatic basis states \f$ \psi \f$
!! are eigenstates of some quantum subsystem operator \f$ \hat q \f$ 
!! such that \f[ \hat q |\psi_j\rangle =f_j|\psi_j\rangle \f]
!! and the eigenvalues \f$ f_j \f$ are supplied by the user.
!! The bilinear coupling term is written
!!   \f[ h_c=-\sum_n^N \sum_k^{K_n} g_{k,n} Q_{k} f_n |\psi_n\rangle\langle\psi_n|\f]
!! with counter term
!!   \f[ ct=\sum_n^N \sum_k^{K_n} \frac{g_{k,n}^2}{2m_{k}\omega_{k}^2}f_n^2 |\psi_n\rangle\langle\psi_n|. \f]
!! The coefficients \f$ g_{k,n} \f$ are sampled from the spectral density \f$ J_n(\omega)\f$ 
!! describing the solvent response about the state \f$ n \f$. The spectral density is written
!! \f[ J_n(\omega)=\frac{\pi}{\hbar}\sum_k^{K_n}\frac{g_{k,n}^2}{2m_{k,n}\omega_{k,n}}\delta_{\omega-\omega_{k,n}}.\f]
!! The re-organization energy for each state is defined
!!   \f[ \lambda_n=\frac{1}{\pi} \int_0^\infty \frac{J_n(\omega)}{\omega} d\omega \f]
!! such that,
!!   \f[ g_{k,n}^2=\frac{2\lambda_n m_{k} \omega_{k}^2}{K_n}\f]
!! and 
!!   \f[ ct=\sum_n^N \lambda_n f_n^2|\psi_n\rangle\langle\psi_n| .\f]
!! NEED TO UPDATE counterterm now shifts classical normal modes on resample
!! \authors
!! Daniel Montemayor
!!
!! \date
!! Aug 2011, Sept 2011, Oct 2012 ,Feb 2013
!<-----------------------------------------------------------------------------
!Change Log:
!=== July 17 2013
! - counterm now shifts classical bath normal modes on resample instead of adding reorginzation energy
!=== Apr 5 2013
! + fix bug: this%g not properly initialize from file - init now resamples at end 
! + allow zero bdof
! + allow for zero lambda and restrict negative lambda
! + moving counter term away from hb and to here.
!
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

module bilinear_class
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

  public::bilinear
  public::new,kill,update,resample,display,save,check

  !> Derived quantum-classical coupling term. 
  type bilinear

     !> True if coupling term has been properly initialized
     logical::initialized=.false.

     !> Inherited coupling primitive type
     type(hc)::hc

     !> Array containing the vaules of the bath reorganization
     !! energy about each state and has dimensions (qs.hs.nstate).
     real(double),dimension(:),pointer::lambda

     !> Number of dissipative bath modes associated to each state has dimensions (qs.hs.nstate). 
     integer(long),dimension(:),pointer::bdof

     !> Coefficient coupling each quantum state to the normal modes of the classical subsystem.
     real(double),dimension(:,:),pointer::g

     !> Eigen values of the \f$ \hat q \f$ operator.
     !! The native diabatic basis states are assumed eigen states
     !! of the \f$ \hat q \f$ operator with eigen values stored in the array f of dimension (qs.hs.nstate).  
     real(double),dimension(:),pointer::f

     !! Logically maps the classical modes to the quantum states. It is true if the mode is coupled to the state.
     !! and has dimensions (cs.hb.ndof by qs.hs.nstate).
     logical,dimension(:,:),pointer::gmap

     !! True if the classical subsystem potential is assummed is diagonal, i.e. non-interacting particles.
     logical::Wbasis

     !! True if counter term is applied
     logical::CTused


  end type bilinear

  !> Creates the bilinear object
  interface new
     module procedure bilinear_init
  end interface

  !> Destroys the bilinear object.
  interface kill
     module procedure bilinear_kill
  end interface

  !> Recalculates bilinear object.
  interface update
     module procedure bilinear_update
  end interface

  !> Reinitializes the bilinear object.
  interface resample
     module procedure bilinear_resample
  end interface

  !> Displays the current state of the bilinear object.
  interface display
     module procedure bilinear_display
  end interface

  !> Saves the current state of the bilinear object to file.
  interface save
     module procedure bilinear_save
  end interface

  !> Checks the bilinear object.
  interface check
     module procedure bilinear_check
  end interface

contains

  !======================================================================
  !> \brief Creates and initializes the bilinear coupling object.
  !> \param this is the bilinear coupling object to be initialized.
  !> \param qs is a general quantum subsystem to be coupled to a classical subsystem.
  !> \param cs is a general classical subsysterm to be coupled to the quantum subsystem.
  !> \param[in] file is an optional string containing the name of a previously saved bilinear file.
  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine bilinear_init(this,qs,cs,file)
    type(bilinear),intent(inout)::this
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


    call Note('Begin bilinear_init.')

    filename='unknown'
    filetype='unknown'

    if(check(qs).EQ.1)call stop('bilinear_init, qs failed check')
    if(check(cs).EQ.1)call stop('bilinear_init, cs failed check')

    if(present(file))then
       if(check(trim(file)).EQ.1)&
            call stop('bilinear_init, cannot find input file')
       
       unit=newunit()
       open(unit,file=file)
       read(unit,*)filetype
       filetype=adjustl(filetype)
       if(trim(filetype).EQ.'bilinear')then
          read(unit,*)filename
          filename=adjustl(filename)
          call new(this%hc,qs,cs,trim(filename))
       else
          call Stop('bilinear_init: input file not valid type.','Check input file format for bilinear object.')
       end if
    else
       call new(this%hc,qs,cs)
    end if
    if(check(this%hc).EQ.1)call stop('bilinear_init: coupling primitive failed check')    


    ndof=this%hc%hb%ndof
    nstate=this%hc%hs%nstate

    if(associated(this%gmap))nullify(this%gmap)
    allocate(this%gmap(ndof,nstate),stat=ierr)
    if(ierr.NE.0)call stop('bilinear_init: failed to allocate memory for mode-state array gmap.')
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
    if(ierr.NE.0)call stop('bilinear_init: failed to allocate memory for bdof.')
    this%bdof=0
    do istate=1,nstate
       do idof=1,ndof
          if(this%gmap(idof,istate))this%bdof(istate)=this%bdof(istate)+1
       end do
    end do

    if(associated(this%lambda))nullify(this%lambda)
    allocate(this%lambda(nstate),stat=ierr)
    if(ierr.NE.0)call stop('bilinear_init: failed to allocate memory for lambda.')
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
    allocate(this%g(ndof,nstate),stat=ierr)
    if(ierr.NE.0)call Stop('bilinear_init: mode-state coefficient memory allocation error.')
    
    if(associated(this%f))nullify(this%f)
    allocate(this%f(nstate),stat=ierr)
    if(ierr.NE.0)call Stop('bilinear_init: qhat operator eigenvalue allocation error.')
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

    if(present(file))then
       read(unit,*)this%CTused
    else
       write(*,*)'Use counter term? T/F '
       read(*,*)this%CTused
    end if

    this%Wbasis=.false.
    if(all(this%hc%hb%EigenVec.eq.iden(ndof)))this%Wbasis=.true.

    this%initialized=.true.

    call resample(this)

    !if(present(file))then
    !   call update(this)
    !else
    !end if

    if(check(this).EQ.1)call stop('bilinear_init: failed final object check')

  end subroutine bilinear_init

  !======================================================================
  !> \brief Destroys the bilinear coupling object.
  !> \param this is the bilinear coupling object to be destroyed.
  !====================================================================
  subroutine bilinear_kill(this)

    type(bilinear),intent(inout)::this

    call Note('Begin bilinear_kill.')
    call kill(this%hc)
    if(associated(this%lambda))nullify(this%lambda)
    if(associated(this%gmap))nullify(this%gmap)
    if(associated(this%bdof))nullify(this%bdof)
    if(associated(this%g))nullify(this%g)
    if(associated(this%f))nullify(this%f)
    this%initialized=.false.
  end subroutine bilinear_kill

  !======================================================================
  !> \brief Computes the current state of bilinear coupling object.
  !> \param this is the bilinear coupling object to be updated.
  !======================================================================
  subroutine bilinear_update(this)
    type(bilinear),intent(inout)::this

    integer(long)::nstate,ndof
    integer(long)::istate,idof
    real(double)::counterterm

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof

    call Note('Begin bilinear_update.')

    !call update(this%hc)

    this%hc%V=0.0_double
    this%hc%dV=0.0_double
    do istate=1,nstate
       !counterterm=0._double
       !if(this%CTused)&
       !     counterterm=this%lambda(istate)*this%f(istate)*this%f(istate)

       do idof=1,ndof
          if(this%gmap(idof,istate))then
             
             this%hc%V(istate,istate)=this%hc%V(istate,istate)&
                  -this%hc%hb%W(idof)*this%g(idof,istate)*this%f(istate)!+counterterm
             this%hc%dV(idof,istate,istate)=&
                  -this%g(idof,istate)*this%f(istate)
          end if
       end do
       !change dV into Q basis
       if(.not.this%Wbasis)this%hc%dV(:,istate,istate)=matmul(transpose(this%hc%hb%EigenVec),this%hc%dV(:,istate,istate))
    end do
    
  end subroutine bilinear_update

  !======================================================================
  !> \brief Re-initiallizes the bilinear coupling object.
  !> \param this is the bilinear coupling object to be re-initialized.
  !======================================================================
  subroutine bilinear_resample(this)
    type(bilinear),intent(inout)::this
    integer(long)::nstate,ndof
    integer(long)::istate,idof
    call Note('Begin bilinear_resample.')

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

    !shift normal mode corrdinates if using counterterm
    if(this%CTused)then
       do idof=1,ndof
          do istate=1,nstate
             if(this%gmap(idof,istate))this%hc%hb%W(idof)=&
                  this%hc%hb%W(idof)-(this%g(idof,istate)&
                  /this%hc%hb%mode(idof)**2)*this%f(istate)
          end do
       end do
       !manually update partical coordinates from shifted normal mode coordinates
       if(.not.this%Wbasis)then
          this%hc%hb%Q=matmul(transpose(this%hc%hb%EigenVec),this%hc%hb%W)
       else
          this%hc%hb%Q=this%hc%hb%W
       end if
    end if


    call update(this)

    if(check(this).EQ.1)call stop('bilinear_resample: object failed exit check.')    
  end subroutine bilinear_resample

  !======================================================================
  !> \brief Saves the current state of the bilinear coupling object to file.
  !> \param[in] this is the bilinear coupling object to be updated.
  !> \param[in] file is a string containing the location of the save file.
  !======================================================================
  subroutine bilinear_save(this,file)
    type(bilinear),intent(in)::this
    character*(*),intent(in)::file

    integer(short)::ierr
    integer(long)::unit
    logical::usedefault

    integer(long)::nstate,ndof
    integer(long)::istate,idof

    call Note('Begin bilinear_save.')
    call Note('input file= '//file)
    if(check(this).EQ.1)call stop('bilinear_save')

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof
    
    unit=newunit()
    open(unit,file=file)
    write(unit,*)'bilinear'
    write(unit,*)quote(file//'.hc')
    call save(this%hc,file//'.hc')
    write(unit,*)((this%gmap(idof,istate),istate=1,nstate),idof=1,ndof)
    write(unit,*)(this%lambda(istate),istate=1,nstate)
    write(unit,*)(this%f(istate),istate=1,this%hc%hs%nstate)
    write(unit,*)this%CTused
    close(unit)      

  end subroutine bilinear_save

  !======================================================================
  !> \brief Displays the bilinear coupling object.
  !> \param[in] this is the bilinear coupling object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will exit after first failed check.
  !======================================================================
  subroutine bilinear_display(this,msg)
    type(bilinear),intent(in)::this
    character*(*),intent(in),optional::msg

    integer(long)::ndof,nstate
    integer(long)::idof,istate
    logical::mask(this%hc%hb%ndof)

    call Note('Begin bilinear_display.')
    if(check(this).NE.0)then
       call warn('bilinear_display: failed check','displaying nothing.')
       return
    end if

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof

    write(Dunit,*)'____________________________________________________'
    write(Dunit,*)'-------------------   bilinear   -------------------'
    if(present(msg))write(Dunit,*)msg
    write(Dunit,*)'----------------------------------------------------'
    write(Dunit,*)'re-organization energy by state (Eh)=',this%lambda
    write(Dunit,*)'re-organization energy by state (1/cm)=',this%lambda/invcm
    write(Dunit,*)'mode to state coupling'
    do istate=1,nstate
       write(Dunit,*)
       write(Dunit,*)'state '//trim(int2str(istate))
       write(Dunit,*)'[ mode, frequency(1/cm), coupling const (hbar/(cm*Angstrom)) ] '
       mask=.false.
       do idof=1,ndof
          if(this%gmap(idof,istate))then
             mask(idof)=.true.
             write(Dunit,*)'[ '//trim(int2str(idof))//', '&
                  //trim(float2str(this%hc%hb%mode(idof)/invcm))//', '&
                  //trim(float2str(this%g(idof,istate)*angstrom/invcm))//' ]'
          end if
       end do
       call display(this%hc%hb%mode/invcm,this%g(:,istate)*angstrom/invcm&
            ,mask=mask,msg='Coupling constant plot: X(1/cm),Y(hbar/(cm*Angstrom))')
    end do
    write(Dunit,*)'qhat operator eigenvalues:',this%f
    call display(this%hc,'bilinear primitive coupling term')
    write(Dunit,*)'===================================================='

  end subroutine bilinear_display

  !======================================================================
  !> \brief Checks the bilinear coupling object.
  !> \param[in] this is the bilinear coupling object to be updated.
  !> \param[in] msg is an optional string message to preface the displayed object.
  !> \param[in] file is an optional string containing the location of a file where the display output gets forwarded.
  !======================================================================
  integer(short) function bilinear_check(this)
    type(bilinear),intent(in)::this
    
    integer(long)::nstate,ndof
    integer(long)::istate,idof
    integer(short)::ierr
    
    call Note('Checking bilinear.')

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof

    bilinear_check=0
    
    !logical::initialized=.false.
    if(.not.this%initialized)then
       call Warn('bilinear_check: coupling primitive not initialized.')
       bilinear_check=1
       return
    end if
    !type(hc)::hc
    if(check(this%hc).NE.0)then
       call Warn('bilinear_check: quantum primitive failed check')
       bilinear_check=1
       return
    end if
    
    !real(double),dimension(:),pointer::lambda
    if(.not.associated(this%lambda))then
       call Warn('bilinear_check: lambda memory not allocated')
       bilinear_check=1
       return
    end if
    if(size(this%lambda).NE.nstate)then
       call Warn('bilinear_check: size of lambda not equal number of states')
       bilinear_check=1
       return
    end if
    if(any(this%lambda.NE.this%lambda))then
       call Warn('bilinear_check: lambda has NAN values.')
       bilinear_check=1
       return
    end if
    if(any(abs(this%lambda).GT.huge(this%lambda)))then
       call Warn('bilinear_check: lambda has Huge values.')
       bilinear_check=1
       return
    end if
    if(any(this%lambda.LT.0._double))then
       call Warn('bilinear_check: lambda has negative values.')
       bilinear_check=1
       return
    end if
    
    !integer(long),dimension(:,:),pointer::gmap
    if(.not.associated(this%gmap))then
       call Warn('bilinear_check: gmap memory not allocated')
       bilinear_check=1
       return
    end if
    if(size(this%gmap,1).NE.ndof)then
       call Warn('bilinear_check: dimension 1 of gmap array not equal number of classcial dofs.')
       bilinear_check=1
       return
    end if
    if(size(this%gmap,2).NE.nstate)then
       call Warn('bilinear_check: dimension 2 of gmap array not equal number of states.')
       bilinear_check=1
       return
    end if

    
    !integer(long),dimension(:),pointer::bdof
    if(.not.associated(this%bdof))then
       call Warn('bilinear_check: bdof memory not allocated')
       bilinear_check=1
       return
    end if
    if(size(this%bdof).NE.nstate)then
       call Warn('bilinear_check: size of bdof array  not equal number of states')
       bilinear_check=1
       return
    end if
    if(any(this%bdof.NE.this%bdof))then
       call Warn('bilinear_check: bdof has NAN values.')
       bilinear_check=1
       return
    end if
    if(any(abs(this%bdof).GE.huge(this%bdof)))then
       call Warn('bilinear_check: bdof has huge values.')
       bilinear_check=1
       return
    end if
    if(any(this%bdof.GT.this%hc%hb%ndof))then
       call Warn('bilinear_check: bdof is greater than total dofs.') 
       bilinear_check=1
       return
    end if
    if(any(this%bdof.LT.0))then
       call Warn('bilinear_check: bdof cannot be negative.')
       bilinear_check=1
       return
    end if
    
    !real(double),dimension(:,:),pointer::g
    if(.not.associated(this%g))then
       call Warn('bilinear_check: mode-state memory not allocated')
       bilinear_check=1
       return
    end if
    if(size(this%g,1).NE.ndof)then
       call Warn('bilinear_check: size of mode-state coefficent array  in dimension 1 does not equal number of classcial dofs.')
       bilinear_check=1
       return
    end if
    if(size(this%g,2).NE.nstate)then
       call Warn('bilinear_check: size of mode-state coefficent array in dimension 2 does not equal number of quatnum states.')
       bilinear_check=1
       return
    end if
    if(any(this%g.NE.this%g))then
       call Warn('bilinear_check: mode-state coefficent array has NAN values.')
       bilinear_check=1
       return
    end if
    if(any(abs(this%g).GT.huge(this%g)))then
       call Warn('bilinear_check: mode-state coefficent array has Huge values.')
       bilinear_check=1
       return
    end if
    if(any(this%g.LT.0._double))then
       call Warn('bilinear_check: mode-state coefficent array has negative values.')
       bilinear_check=1
       return
    end if

    !real(double),dimension(:),pointer::f
    if(.not.associated(this%f))then
       call Warn('bilinear_check: qhat eigenvalue array memory not allocated')
       bilinear_check=1
       return
    end if
    if(size(this%f).NE.nstate)then
       call Warn('bilinear_check: size of qhat eigenvalue array not equal number of quatnum states')
       bilinear_check=1
       return
    end if
    if(any(this%f.NE.this%f))then
       call Warn('bilinear_check: qhat eigenvalue array has NAN values.')
       bilinear_check=1
       return
    end if 
   if(any(abs(this%f).GT.huge(this%f)))then
       call Warn('bilinear_check: qhat eigenvalue array has Huge values.')
       bilinear_check=1
       return
    end if
    
  end function bilinear_check
  
end module bilinear_class

