!-----------------------------------------------------------------------------
!> \brief
!!Bi-linear Coupling.
!! \details
!! This coupling scheme is a Caldeira-Leggett model, coupling the normal modes of the classical subsystem linearly to the diagonal states of the quantum subsystem. Assumes the Hamiltonian
!!        \f[  H=hs+hb+H_c+ct \f]
!! where
!!         \f[ hs=\frac{\hat p^2}{2m}+V(\hat q) \f]
!! ,       \f[ \langle\psi_i|hs|\psi_j\rangle =\f]
!! is the diabats attribute of hs
!! , and       \f[ hb=V(Q) \f]
!!  projected onto the normal modes of hb (ie modes attribute of hb).
!! Each matrix element \f$ \langle\psi_n|hs|\psi_m\rangle \f$ is coupled to \f$ K \f$ modes such that
!! the bilinear coupling term is written
!!   \f[ H_c=-\sum_n^N \sum_m^N \sum_k^{K_{n,m}} c_{k,n,m}} Q_{k} q_{n,m} \f]
!! with counter term
!!   \f[ ct=\sum_n^N \sum_m^N \sum_k^{K_{n,m}} q^2_{n,m} \frac{c_{k,n,m}^2}{2m_k\omega_k^2}. \f].
!! Here \f$ q_{n,m} \f$ and \f$ q^2_{n,m} \f$ are quantum subsystem operator matrix elements
!!    \f[ \langle \psi_n|\hat q|\psi_m \rangle \f]
!! and
!!    \f[ \langle \psi_n|\hat q \hat q|\psi_m \rangle \f]
!! respectively.
!! The spectral density associated with matrix element \f$ \langle n|H|m\rangle \f$, abbreviated \f$ H_{n,m} \f$, is defined
!! \f[ J_{n,m}(\omega)=\frac{\pi}{\hbar}\sum_k^{K_{n,m}\frac{c_{k,n,m}^2}{2m_k \omega_k}\delta_{\omega-\omega_k}.\f]
!! The bath re-organization energy \f$ \lambda_{n,m} \f$ associated to \f$ H_{n,m} \f$ is defined in terms of the normal modes of the classical subsystem coupled to that matrix element.
!! \f$ \lambda_{n,m} \f$ is computed on the frequency domain from 0 to \f$ \omega_{max} \f$, where \f$ \omega_{max} \f$ is the largest frequency mode of classical subsystem associated with \f$ H_{n,m} \f$ and is written
!!   \f[ \lambda_{n,m}=\frac{2}{\pi} \int_0^{\omega_{max}} \frac{J_{n,m}(\omega)}{\omega} d\omega \f]
!! such that,
!!   \f[ c_{k,n,m}^2=\frac{\hbar \lambda_{n,m} m_k \omega_k^2}{K_{n,m}}\f]
!! and 
!!   \f[ ct=2 \hbar \sum_n^N \sum_m^N \lambda_{n,m} q^2_{n,m} .\f]
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

module CaldeiraLeggett_class
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

  public::CaldeiraLeggett
  public::new,kill,update,resample,display,save,check

  type CaldeiraLeggett
     logical::initialized=.false.
     type(hc)::hc
     real(double),dimension(:),pointer::X
     integer(long),dimension(:,:),pointer::bdof
     real(double),dimension(:,:),pointer::lambda
     complex(double),dimension(:,:),pointer::f,f2
     real(double),dimension(:,:,:),pointer::g
     logical,dimension(:,:,:),pointer::gmap
  end type CaldeiraLeggett

  interface new
     module procedure CaldeiraLeggett_init
  end interface

  interface kill
     module procedure CaldeiraLeggett_kill
  end interface

  interface update
     module procedure CaldeiraLeggett_update
  end interface

  interface resample
     module procedure CaldeiraLeggett_resample
  end interface

  interface display
     module procedure CaldeiraLeggett_display
  end interface

  interface save
     module procedure CaldeiraLeggett_save
  end interface

  interface check
     module procedure CaldeiraLeggett_check
  end interface

contains
!======================================================================
  subroutine CaldeiraLeggett_init(this,qs,cs,file)
    type(CaldeiraLeggett),intent(inout)::this
    type(quantum),intent(inout),target::qs
    type(classical),intent(inout),target::cs
    character*(*),intent(in),optional::file      
    character(len=path)::filename
    character(len=title)::filetype
    integer(short)::ierr
    integer(long)::unit
    logical::usedefault,pass

    integer(long)::ndof,nstate,bdof
    integer(long)::idof,istate,jstate
    real(double)::norm


    call Note('Begin CaldeiraLeggett_init.')
    filename='unknown'
    filetype='unknown'

    if(check(qs).EQ.1)call stop('CaldeiraLeggett_init, qs failed check')
    if(check(cs).EQ.1)call stop('CaldeiraLeggett_init, cs failed check')

    if(present(file))call Note('input file= '//trim(file))

    if(present(file))then
       if(check(trim(file)).EQ.1)&
            call stop('CaldeiraLeggett_init, cannot find input file')
       
       unit=newunit()
       open(unit,file=file)
       read(unit,*)filetype
       filetype=adjustl(filetype)
       if(trim(filetype).EQ.'CaldeiraLeggett')then
          read(unit,*)filename
          filename=adjustl(filename)
          call new(this%hc,qs,cs,trim(filename))
       else
          call Stop('CaldeiraLeggett_init: input file not valid type.','Check input file format for CaldeiraLeggett object.')
       end if
    else
       call new(this%hc,qs,cs)
    end if
    if(check(this%hc).EQ.1)call stop('CaldeiraLeggett_init: coupling primitive failed check')    


    ndof=this%hc%hb%ndof
    nstate=this%hc%hs%nstate

    if(associated(this%gmap))nullify(this%gmap)
    allocate(this%gmap(ndof,nstate,nstate),stat=ierr)
    if(ierr.NE.0)call stop('CaldeiraLeggett_init: failed to allocate memory for mode-state array gmap.')
    this%gmap=.false.
    if(present(file))then
       read(unit,*)(((this%gmap(idof,istate,jstate)&
            ,jstate=1,nstate),istate=1,nstate),idof=1,ndof)
    else
       pass=.false.
       if(mod(ndof,nstate).EQ.0)pass=.true.
       usedefault=.false.
       if(pass)then
          write(*,*)'Use default mode to state mapping? [T/F] (Divide modes equally among Hailtonian matrix elements.)'
          read(*,*)usedefault
       end if
       if(usedefault)then
          bdof=ndof/(nstate**2)
          do istate=1,nstate
             do jstate=1,nstate
                do idof=1,bdof
                   this%gmap(idof+(istate-1)*bdof*nstate+(jstate-1)*bdof&
                        ,istate,jstate)=.true.
                end do
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
                   do jstate=1,nstate
                      write(*,*)'Associate mode '//trim(int2str(idof))//' with matrix element '//trim(int2str(istate))//','//trim(int2str(jstate))//'? T/F'
                      read(*,*)this%gmap(idof,istate,jstate)
                   end do
                end do
             end do
          end if
       end if
    end if
    if(associated(this%bdof))nullify(this%bdof)
    allocate(this%bdof(nstate,nstate),stat=ierr)
    if(ierr.NE.0)call stop('CaldeiraLeggett_init: failed to allocate memory for bdof.')
    this%bdof=0
    do istate=1,nstate
       do istate=1,nstate
          do idof=1,ndof
             if(this%gmap(idof,istate,jstate))this%bdof(istate,jstate)=this%bdof(istate,jstate)+1
          end do
       end do
    end do
    
    if(associated(this%lambda))nullify(this%lambda)
    allocate(this%lambda(nstate,nstate),stat=ierr)
    if(ierr.NE.0)call stop('CaldeiraLeggett_init: failed to allocate memory for lambda.')
    this%lambda=0._double
    if(present(file))then
       read(unit,*)((this%lambda(istate,jstate),jstate=1,nstate),istate=1,nstate)
    else
       do istate=1,nstate
          do jstate=1,nstate
             write(*,*)'Please enter reorganization energy for state '&
                  //trim(int2str(istate))//','//trim(int2str(jstate))
             read(*,*)this%lambda(istate,jstate)
          end do
       end do
    end if
    
    if(associated(this%g))nullify(this%g)
    allocate(this%g(ndof,nstate,nstate),stat=ierr)
    if(ierr.NE.0)call Stop('CaldeiraLeggett_init: mode-state coefficient memory allocation error.')
    
    if(associated(this%f))nullify(this%f)
    allocate(this%f(nstate,nstate),stat=ierr)
    if(ierr.NE.0)call Stop('CaldeiraLeggett_init: quantum subsystem operator matrix allocation error.')
    this%f=0.0_double
    if(present(file))then
       read(unit,*)((this%f(istate,jstate),jstate=1,nstate),istate=1,nstate)
    else
       usedefault=.true.
       write(*,*)'Use default quantum subsystem operator matrix values? [all = 1.0 on diagonal] T/F '
       read(*,*)usedefault
       if(usedefault)then
          do istate=1,nstate
             this%f(istate,istate)=1.0_double
          end do
       else
          write(*,*)'Enter quantum subsystem operator matrix.'
          do istate=1,nstate
             do jstate=istate,nstate
                write(*,*)'Element '//trim(int2str(istate))&
                     //','//trim(int2str(jstate))
                read(*,*)this%f(istate,jstate)
                this%f(jstate,istate)=conjg(this%f(istate,jstate))
             end do
          end do
       end if
    end if

    if(associated(this%f2))nullify(this%f2)
    allocate(this%f2(nstate,nstate),stat=ierr)
    if(ierr.NE.0)call Stop('CaldeiraLeggett_init: quantum subsystem operator^2 matrix allocation error.')
    this%f2=matmul(this%f,this%f)

    if(associated(this%X))nullify(this%X)
    allocate(this%X(ndof),stat=ierr)
    if(ierr.NE.0)call Stop('CaldeiraLeggett_init: Normal mode coordinates memory allocation error.')
    this%X=0.0_double

    this%initialized=.true.
    call resample(this)
    if(check(this).EQ.1)call stop('CaldeiraLeggett_init: failed final object check')

  end subroutine CaldeiraLeggett_init
  !-=-=-=-=-=-=-=-=-=-==--=-=-=-=-=-=-=-=-=-=-=-
  subroutine CaldeiraLeggett_kill(this)

    type(CaldeiraLeggett),intent(inout)::this

    call Note('Begin CaldeiraLeggett_kill.')
    call kill(this%hc)
    if(associated(this%lambda))nullify(this%lambda)
    if(associated(this%gmap))nullify(this%gmap)
    if(associated(this%bdof))nullify(this%bdof)
    if(associated(this%g))nullify(this%g)
    if(associated(this%f))nullify(this%f)
    if(associated(this%f2))nullify(this%f2)
    if(associated(this%X))nullify(this%X)
    this%initialized=.false.
  end subroutine CaldeiraLeggett_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine CaldeiraLeggett_update(this)
    type(CaldeiraLeggett),intent(inout)::this

    integer(long)::nstate,ndof
    integer(long)::istate,jstate,idof
    complex(double)::counterterm

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof

    call Note('Begin CaldeiraLeggett_update.')

    !call update(this%hc)

    this%X=0._double
    do idof=1,ndof
       this%X=sum(this%hc%hb%EigenVec(:,idof)*this%hc%hb%Q)
    end do

    this%hc%V=0.0_double
    this%hc%dV=0.0_double
    do istate=1,nstate
       do jstate=1,nstate
          do idof=1,ndof
             if(this%gmap(idof,istate,jstate))then
                counterterm=0._double
                
                if(this%hc%hb%canonical)&
                     counterterm=.5_double*this%hc%hb%rmass(idof)&
                     *(this%g(idof,istate,jstate)/this%hc%hb%mode(idof))**2&
                     *this%f2(istate,jstate)
             
                this%hc%V(istate,jstate)=this%hc%V(istate,jstate)&
                  -this%X(idof)*this%g(idof,istate,jstate)&
                  *this%f(istate,jstate)+counterterm
                !+2_double*this%lambda(istate,jstate)*this%f(istate)!**2counterterm

                this%hc%dV(idof,istate,jstate)=&
                     -this%g(idof,istate,jstate)*this%f(istate,jstate)
             end if
          end do
       end do
    end do
    
  end subroutine CaldeiraLeggett_update
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine CaldeiraLeggett_resample(this)
    type(CaldeiraLeggett),intent(inout)::this
    integer(long)::nstate,ndof
    integer(long)::istate,jstate,idof
    real(double)::wr
    call Note('Begin CaldeiraLeggett_resample.')

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof

    !call resample(this%hc)

    !resample normal mode oscillator strengths
    this%g=0._double
    do istate=1,nstate
       do jstate=1,nstate
          wr=hbar*this%lambda(istate,jstate)/real(this%bdof(istate,jstate))
          do idof=1,ndof
             if(this%gmap(idof,istate,jstate))this%g(idof,istate,jstate)=&
                  this%hc%hb%mode(idof)*sqrt(wr/this%hc%hb%rmass(idof))      
          end do
       end do
    end do

    call update(this)
    
  end subroutine CaldeiraLeggett_resample
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine CaldeiraLeggett_save(this,file)
    type(CaldeiraLeggett),intent(in)::this
    character*(*),intent(in)::file

    integer(short)::ierr
    integer(long)::unit
    logical::usedefault

    integer(long)::nstate,ndof
    integer(long)::istate,jstate,idof

    call Note('Begin CaldeiraLeggett_save.')
    call Note('input file= '//file)
    if(check(this).EQ.1)call stop('CaldeiraLeggett_save')

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof
    
    unit=newunit()
    open(unit,file=file)
    write(unit,*)'CaldeiraLeggett'
    write(unit,*)quote(file//'.hc')
    call save(this%hc,file//'.hc')
    write(unit,*)((this%gmap(idof,istate,jstate)&
         ,jstate=1,nstate),istate=1,nstate),idof=1,ndof)
    write(unit,*)((this%lambda(istate,jstate),jstate=1,nstate),istate=1,nstate)
    write(unit,*)((this%f(istate,jstate),jstate=1,nstate),istate=1,nstate)
    close(unit)      

  end subroutine CaldeiraLeggett_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine CaldeiraLeggett_display(this,msg)
    type(CaldeiraLeggett),intent(in)::this
    character*(*),intent(in),optional::msg

    integer(long)::ndof,nstate
    integer(long)::idof,istate,jstate
    logical::mask(this%hc%hb%ndof)

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof

    call Note('Begin CaldeiraLeggett_display.')
    if(check(this).EQ.0)then
       write(*,*)'################## CaldeiraLeggett ################'
       write(*,*)'re-organization energy by state='!,this%lambda
       call display(this%lambda)
       write(*,*)

       write(*,*)'mode to state mapping'
       do istate=1,nstate
          do jstate=1,nstate
             write(*,*)
             write(*,*)'element '//trim(int2str(istate))&
                  //','//trim(int2str(istate))
             write(*,'(A7)',advance='no')'modes: '
             mask=.false.
             do idof=1,ndof
                if(this%gmap(idof,istate,jstate))then
                   mask(idof)=.true.
                   write(*,'(A3,1X)',advance='no')trim(int2str(idof))
                end if
             end do
             call display(this%hc%hb%mode,this%g(:,istate,jstate),mask=mask)
             write(*,*)
          end do
       end do
       write(*,*)
          
       write(*,*)'quantum subsystem operator matrix:'!,this%f
       call display(this%f)
       write(*,*)
       !write(*,*)'mode(col)-state(row) interaction terms'
       !call display(transpose(this%g(:,:,jstate))
       write(*,*)
       call display(this%hc,'CaldeiraLeggett primitive coupling term')
       write(*,*)'############### end CaldeiraLeggett ################'
    else
       call warn('CaldeiraLeggett_display: coupling object failed check.','displaying nothing.')
    end if
    
  end subroutine CaldeiraLeggett_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short) function CaldeiraLeggett_check(this)
    type(CaldeiraLeggett),intent(in)::this
    
    integer(long)::nstate,ndof
    integer(long)::istate,idof
    integer(short)::ierr
    
    call Note('Checking CaldeiraLeggett.')

    nstate=this%hc%hs%nstate
    ndof=this%hc%hb%ndof

    CaldeiraLeggett_check=0
    
    !logical::initialized=.false.
    if(.not.this%initialized)then
       call Warn('CaldeiraLeggett_check: coupling primitive not initialized.')
       CaldeiraLeggett_check=1
       return
    end if
    !type(hc)::hc
    if(check(this%hc).NE.0)then
       call Warn('CaldeiraLeggett_check: quantum primitive failed check')
       CaldeiraLeggett_check=1
       return
    end if
    
    !real(double),dimension(:),pointer::lambda
    if(.not.associated(this%lambda))then
       call Warn('CaldeiraLeggett_check: lambda memory not allocated')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%lambda,1).NE.nstate)then
       call Warn('CaldeiraLeggett_check: size of lambda dimension 1 not equal number of states')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%lambda,2).NE.nstate)then
       call Warn('CaldeiraLeggett_check: size of lambda dimension 2 not equal number of states')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(this%lambda.NE.this%lambda))then
       call Warn('CaldeiraLeggett_check: lambda has NAN values.')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(abs(this%lambda).GT.huge(this%lambda)))then
       call Warn('CaldeiraLeggett_check: lambda has Huge values.')
       CaldeiraLeggett_check=1
       return
    end if
    
    !real(double),dimension(:),pointer::X
    if(.not.associated(this%X))then
       call Warn('CaldeiraLeggett_check: X memory not allocated')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%X).NE.this%hc%hb%ndof)then
       call Warn('CaldeiraLeggett_check: size of X not equal number of classcial dofs')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(this%X.NE.this%X))then
       call Warn('CaldeiraLeggett_check: Normal mode coordinate array X has NAN values.')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(abs(this%X).GT.huge(this%X)))then
       call Warn('CaldeiraLeggett_check: Normal mode coordinate array X has Huge values.')
       CaldeiraLeggett_check=1
       return
    end if
    
    !integer(long),dimension(:,:),pointer::gmap
    if(.not.associated(this%gmap))then
       call Warn('CaldeiraLeggett_check: gmap memory not allocated')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%gmap,1).NE.ndof)then
       call Warn('CaldeiraLeggett_check: dimension 1 of gmap array not equal number of classcial dofs.')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%gmap,2).NE.nstate)then
       call Warn('CaldeiraLeggett_check: dimension 2 of gmap array not equal number of states.')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%gmap,3).NE.nstate)then
       call Warn('CaldeiraLeggett_check: dimension 3 of gmap array not equal number of states.')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(this%gmap.NE.this%gmap))then
       call Warn('CaldeiraLeggett_check: gmap array has bad values.')
       CaldeiraLeggett_check=1
       return
    end if

    
    !integer(long),dimension(:),pointer::bdof
    if(.not.associated(this%bdof))then
       call Warn('CaldeiraLeggett_check: bdof memory not allocated')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%bdof,1).NE.nstate)then
       call Warn('CaldeiraLeggett_check: size of bdof array dimension 1 not equal number of states')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%bdof,2).NE.nstate)then
       call Warn('CaldeiraLeggett_check: size of bdof array dimension 2 not equal number of states')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(this%bdof.NE.this%bdof))then
       call Warn('CaldeiraLeggett_check: bdof has NAN values.')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(abs(this%bdof).GE.huge(this%bdof)))then
       call Warn('CaldeiraLeggett_check: bdof has huge values.')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(this%bdof.GT.this%hc%hb%ndof))then
       call Warn('CaldeiraLeggett_check: bdof is greater than total dofs.') 
       CaldeiraLeggett_check=1
       return
    end if
    if(any(this%bdof.LE.0))then
       call Warn('CaldeiraLeggett_check: bdof must be positive.')
       CaldeiraLeggett_check=1
       return
    end if
    
    !real(double),dimension(:,:),pointer::g
    if(.not.associated(this%g))then
       call Warn('CaldeiraLeggett_check: mode-state memory not allocated')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%g,1).NE.ndof)then
       call Warn('CaldeiraLeggett_check: size of mode-state coefficent array  in dimension 1 does not equal number of classcial dofs.')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%g,2).NE.nstate)then
       call Warn('CaldeiraLeggett_check: size of mode-state coefficent array in dimension 2 does not equal number of quatnum states.')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%g,3).NE.nstate)then
       call Warn('CaldeiraLeggett_check: size of mode-state coefficent array in dimension 3 does not equal number of quatnum states.')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(this%g.NE.this%g))then
       call Warn('CaldeiraLeggett_check: mode-state coefficent array has NAN values.')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(abs(this%g).GT.huge(this%g)))then
       call Warn('CaldeiraLeggett_check: mode-state coefficent array has Huge values.')
       CaldeiraLeggett_check=1
       return
    end if

    !complex(double),dimension(:,:),pointer::f
    if(.not.associated(this%f))then
       call Warn('CaldeiraLeggett_check: quantum subsystem operator matrix memory not allocated')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%f,1).NE.nstate)then
       call Warn('CaldeiraLeggett_check: size of quantum subsystem operator matrix dimension 1 not equal number of quatnum states')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%f,2).NE.nstate)then
       call Warn('CaldeiraLeggett_check: size of quantum subsystem operator matrix dimension 2 not equal number of quatnum states')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(this%f.NE.this%f))then
       call Warn('CaldeiraLeggett_check: quantum subsystem operator matrix has NAN values.')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(abs(this%f).GT.huge(this%f)))then
       call Warn('CaldeiraLeggett_check: quantum subsystem operator matrix has Huge values.')
       CaldeiraLeggett_check=1
       return
    end if

    !complex(double),dimension(:,:),pointer::f2
    if(.not.associated(this%f2))then
       call Warn('CaldeiraLeggett_check: quantum subsystem operator^2 matrix memory not allocated')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%f2,1).NE.nstate)then
       call Warn('CaldeiraLeggett_check: size of quantum subsystem operator^2 matrix dimension 1 not equal number of quatnum states')
       CaldeiraLeggett_check=1
       return
    end if
    if(size(this%f2,2).NE.nstate)then
       call Warn('CaldeiraLeggett_check: size of quantum subsystem operator^2 matrix dimension 2 not equal number of quatnum states')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(this%f2.NE.this%f2))then
       call Warn('CaldeiraLeggett_check: quantum subsystem operator^2 matrix has NAN values.')
       CaldeiraLeggett_check=1
       return
    end if
    if(any(abs(this%f2).GT.huge(this%f2)))then
       call Warn('CaldeiraLeggett_check: quantum subsystem operator^2 matrix has Huge values.')
       CaldeiraLeggett_check=1
       return
    end if
    
  end function CaldeiraLeggett_check
  
end module CaldeiraLeggett_class

