!============================================================================
!> \brief
!! Quantum subsystem primitive class.
!!
!! \details
!! This class defines the quantum subsystem primitive type \f$hs/f$ inherited 
!! by all derived quantum subsystems. Methods include the standard interface
!! set: NEW, KILL, UPDATE, RESAMPLE, DISPLAY, SAVE, and CHECK.
!!
!! \todo
!! Make diabat complex
!!
!!\date
!! 12 Nov 2012
!!\authors 
!!Daniel Montemayor
!<============================================================================
!============================ C H A N G E   L O G ============================
!- corrected initial density input for non-default case, July 2012 Daniel M.
!=============================================================================
module hs_class
  use type_kinds
  use filemanager
  use ErrorLog
  use string
  use atomicunits
  use textgraphs
  use outputdisplay
  implicit none
  private

  public::hs
  public::new, kill, update, resample, display, save, check


!============================================================================
!> \brief
!! Quantum subsystem primitive type.
!<===========================================================================
  type hs
     !> true if hs type has been properly initialized
     logical::initialized=.false.
     !> defines the number of quantum states.
     integer(long) :: nstate 
     !>is the energy origin often defined to be the ground state energy.
     real(double)::Eg 
     !>is a vestigial parameter to be removed in later versions.
     real(double)::tol 
     !> is an NxN matrix containing the elements of the quantum 
     !!        subsystem hamiltonian in the native diabatic basis.
     real(double),dimension(:,:),pointer::diabat 
     !>is an NxN matrix containing along its columns the expansion
     !!        coefficients that transform the diabatic basis into the eigen state
     !!        basis of the quantum subsystem hamiltonian
     complex(double),dimension(:,:),pointer::EigenVec 
     !>is array of length N containing the eigen vaules of the
     !!        eigen states.
     real(double),dimension(:),pointer::EigenVal 
     !>is an NxN matrix containing the reduced density matrix of the
     !!        total Hamiltonian. It describes the current state of the quantum
     !!        subsystem in the diabatic basis
     complex(double),dimension(:,:),pointer::den 
  end type hs

  !> Creates hs type
  interface new
     module procedure hs_init
  end interface

  !> Destroys hs type
  interface kill
     module procedure hs_kill
  end interface

  !> \brief
  !!  Null method
  !! \details
  !! The hs type cannot update itself. It relies on a derived 
  !! quantum subsystem type to update its attributes.
  interface update
     module procedure hs_update
  end interface

  !> Null method
  !! \details
  !! The hs type cannot resample itself. It relies on a derived 
  !! quantum subsystem type to reinitialize its attributes.
  interface resample
     module procedure hs_resample
  end interface

  !> Prints out current state of the hs type
  interface display
     module procedure hs_display
  end interface

  !> Saves the current state of the hs type to file 
  interface save
     module procedure hs_save
  end interface

  !> \brief Check is a function that checks if the attributes of the hs type
  !! are within acceptable values
  !> \return
  !! The integer 0 if no attribute check fails. The first failed attribute
  !! check will cause this function to return the integer 1 and issue
  !! a warning stating the failed check. This warning may cause the program
  !! to stop according to the loglevel set by the user.
  interface check
     module procedure hs_check
  end interface

contains
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hs_init(this,file)
    type(hs),intent(inout)::this
    character*(*),intent(in),optional::file

    character(len=title)::filetype
    integer(short)::ierr
    integer(long)::istate,jstate,rval,imval,trace

    integer(long)::unit
    logical::usedefault

    call Note('Begin hs_init.')
    if(present(file))call Note('input file= '//file)

    !ratio of allowed imaginary component of diagonal density matrix element
    !   to density trace. 
    this%tol=1E-5


    if(present(file))then

       if(check(file))call stop('hs_init') 

       unit=newunit()
       open(unit,file=file)

       read(unit,*)filetype
       filetype=adjustl(filetype)
       if(trim(filetype).NE.'hs')call stop('hs_init: input file not valid.')
    end if

    if(present(file))then
       read(unit,*)this%nstate
    else
       this%nstate=0
       do while (this%nstate.LE.0)
          write(*,*)'Please enter the number of quantum states.'
          read(*,*)this%nstate
          if(this%nstate.LE.0)write(*,*)'Number of quantum states must be positive integer. Try again.'
       end do
    end if

    this%Eg=0.0_double
    if(present(file))then
       read(unit,*)this%Eg
    else
       write(*,*)'Please enter an energy origin.'
       read(*,*)this%Eg
    end if
    
    if(associated(this%diabat))nullify(this%diabat)
    allocate(this%diabat(this%nstate,this%nstate),stat=ierr)
    if(ierr.NE.0)call stop('hs_init: failed diabatic energy memory allocation.')
    this%diabat(:,:)=0.0
    if(present(file))then
       read(unit,*)((this%diabat(istate,jstate),jstate=1,this%nstate),istate=1,this%nstate)
    else
       write(*,*)'Use default diabatic energy settings [zero engergy for all states]? (T/F)'
       read(*,*)usedefault
       if(.not.usedefault)then
          do istate=1,this%nstate
             do jstate=istate,this%nstate
                if(istate.EQ.jstate)then
                   write(*,*)'Enter diabatic state ',istate,'energy.'
                else
                   write(*,*)'Enter diabatic state coupling energy for states',istate,jstate
                end if
                read(*,*)this%diabat(istate,jstate)
                this%diabat(jstate,istate)=this%diabat(istate,jstate)
             end do
          end do
       end if
    end if
    
    if(associated(this%den))nullify(this%den)
    allocate(this%den(this%nstate,this%nstate),stat=ierr)
    if(ierr.NE.0)call stop('hs_init: failed density matrix memory allocation.')

    this%den=(0._double,0._double)
    this%den(1,1)=(1._double,0._double)
    if(present(file))then
       read(unit,*)((this%den(istate,jstate),jstate=1,this%nstate),istate=1,this%nstate)
    else
       write(*,*)'Use default initial density settings [only state 1 populated]? (T/F)'
       read(*,*)usedefault
       if(.not.usedefault)then
          write(*,*)'Enter density for each element.'
          do istate=1,this%nstate
             write(*,*)'|',istate,'><',istate,'| (real number only)'
             read(*,*)rval
             this%den(istate,istate)=rval
             if(istate.LT.this%nstate)then
                do jstate=istate+1,this%nstate
                   write(*,*)'|',istate,'><',jstate,'|= (2 numbers: real followed by imaginary)'
                   read(*,*)rval,imval
                   this%den(istate,jstate)=rval+eye*imval
                   this%den(jstate,istate)=rval-eye*imval
                end do
             end if
          end do
       end if
    end if
!!$    trace=0_double
!!$    do istate=1,this%nstate
!!$       trace=trace+real(this%den(istate,istate))
!!$    end do
!!$    this%den=this%den/trace

    if(associated(this%EigenVec))nullify(this%EigenVec)
    allocate(this%EigenVec(this%nstate,this%nstate),stat=ierr)
    if(ierr.NE.0)call stop('hs_init: failed Eigen vector matrix memory allocation.')

    if(associated(this%EigenVal))nullify(this%EigenVal)
    allocate(this%EigenVal(this%nstate),stat=ierr)
    if(ierr.NE.0)call stop('hs_init: failed Eigen value array memory allocation.')

    this%EigenVec=this%diabat
    call diagonalize(this%nstate,this%EigenVec,this%EigenVal)

    if(present(file))close(unit)

    this%initialized=.true.
    if(check(this))call stop('hs_init, failed final object check')

  end subroutine hs_init
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hs_kill(this)
    type(hs),intent(inout)::this

    call Note('Begin hs_kill.')

    if(associated(this%diabat))nullify(this%diabat)
    if(associated(this%den))nullify(this%den)
    if(associated(this%EigenVec))nullify(this%EigenVec)
    if(associated(this%EigenVal))nullify(this%EigenVal)
    this%initialized=.false.

  end subroutine hs_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hs_update(this)
    type(hs),intent(inout)::this

    call Note('Begin hs_update.')

  end subroutine hs_update
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hs_resample(this)
    type(hs),intent(inout)::this
    
    call Note('Begin hs_resample.')

    this%EigenVec=this%diabat
    call diagonalize(this%nstate,this%EigenVec,this%EigenVal)

    call update(this)

  end subroutine hs_resample
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hs_display(this,msg)
    type(hs),intent(in)::this
    character*(*),intent(in),optional::msg
    integer(long)::istate,jstate
    integer(short)::ierr
    real(double)::diabat(this%nstate)

    call Note('Begin hs_display.')
    if(check(this).NE.0)then
       call warn('hs_display: failed check','displaying nothing.')
       return
    end if

    write(Dunit,*)'____________________________________________________'
    write(Dunit,*)'----------------------- hs -------------------------'
    if(present(msg))write(Dunit,*)msg
    write(Dunit,*)'----------------------------------------------------'
    write(Dunit,*)'Nstate=',this%nstate
    write(Dunit,*)'Energy origin=',this%Eg
    write(Dunit,*)'Energy origin(1/cm)=',this%Eg/invcm
    write(Dunit,*)'Diagonal diabatic energies (1/cm)=',(this%diabat(istate,istate)/invcm,istate=1,this%nstate)
    write(Dunit,*)'Adiabatic energies (1/cm)=',(this%EigenVal(jstate)/invcm,jstate=1,this%nstate)
    call display(this%diabat,msg='Diabatic Hamiltonian')
    call display(this%EigenVec,msg='Adiabatic Transformation Matrix')
    call display(this%den,msg='Density Matrix')
    write(Dunit,*)'===================================================='

  end subroutine hs_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine hs_save(this,file)
    type(hs),intent(in)::this
    character*(*),intent(in)::file
    integer(long)::n,istate,jstate
    integer(long)::unit      
    
    call note('Begin hs_save.')
    call note('input file= '//file)
    if(check(this))call stop('hs_save')

    unit=newunit()
    open(unit,file=file)
    n=this%nstate
    write(unit,*)'hs'
    write(unit,*)n
    write(unit,*)this%Eg
    write(unit,*)((this%diabat(istate,jstate),jstate=1,n),istate=1,n)
    write(unit,*)((this%den(istate,jstate),jstate=1,n),istate=1,n)
    close(unit)

  end subroutine hs_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short)function hs_check(this)
    type(hs),intent(in)::this

    integer(long)::istate,jstate
    real(double)::trace
    integer(short)::ierr

    call note('Checking hs.')
    hs_check=0

    !logical::initialized=.false.
    if(.not.this%initialized)then
       call Warn('hs_check: quantum primative not initialized.')
       hs_check=1
       return
    end if

    !integer(long) :: nstate
    if(this%nstate.NE.this%nstate)then
       call Warn('hs_check: nstate not a number.')
       hs_check=1
       return
    end if
    if(abs(this%nstate).GE.huge(this%nstate))then
       call Warn('hs_check: nstate is too big.')
       hs_check=1
       return
    end if
    if(this%nstate.LE.0)then
       call Warn('hs_check: nstate not a positive integer.')
       hs_check=1
       return
    end if

    !real(double)::Eg
    if(this%Eg.NE.this%Eg)then
       call Warn('hs_check: Energy origin is not a number.')
       hs_check=1
       return
    end if
    if(abs(this%Eg).GE.huge(this%Eg))then
       call Warn('hs_check: Energy origin absolute value is too big.')
       hs_check=1
       return
    end if

    !real(double),dimension(:,:),pointer::diabat
    if(.not.associated(this%diabat))then
       call Warn('hs_check: diabatic matrix memory not associated.')
       hs_check=1
       return
    end if
    if(size(this%diabat).NE.this%nstate**2)then
       call Warn('hs_check: number of diabatic matrix elements must equal number of quantum states squared.')
       hs_check=1
       return
    end if
    ierr=0
    do istate=1,this%nstate
       do jstate=1, this%nstate
          if(this%diabat(istate,jstate).NE.this%diabat(istate,jstate))ierr=1
          if(abs(this%diabat(istate,jstate)).GE.huge(this%diabat(istate,jstate)))ierr=1
       end do
    end do
    if(ierr.NE.0)then
       call Warn('hs_check: diabatic matrix elements have bad values.')
       hs_check=1
       return
    end if

    !complex(double),dimension(:,:),pointer::den
    if(.not.associated(this%den))then
       call Warn('hs_check: density matrix memory not associated.')
       hs_check=1
       return
    end if
    if(size(this%den).NE.this%nstate**2)then
       call Warn('hs_check: number of density matrix elements must equal number of quantum states squared.')
       hs_check=1
       return
    end if

    if(any(this%den.NE.this%den))then
       call Warn('hs_check: density matrix has NAN values')
       hs_check=1
       return
    end if

    if(any(abs(real(this%den)).GE.huge(real(this%den))))then
       call Warn('hs_check: real parts of density matrix elements have huge values')
       hs_check=1
       return
    end if


    if(any(abs(aimag(this%den)).GE.huge(aimag(this%den))))then
       call Warn('hs_check: imaginary parts of density matrix elements have huge values')
       hs_check=1
       return
    end if
        
    !complex(double),dimension(:,:),pointer::EigenVec
    if(.not.associated(this%EigenVec))then
       call Warn('hs_check: Eigen vector memory not associated.')
       hs_check=1
       return
    end if
    if(size(this%EigenVec).NE.this%nstate**2)then
       call Warn('hs_check: Eigen Vector elements must equal number of quantum states squared.')
       hs_check=1
       return
    end if
    ierr=0
    do istate=1,this%nstate
       do jstate=1,this%nstate
          if(this%EigenVec(istate,jstate).NE.this%EigenVec(istate,jstate))ierr=1
          if(abs(real(this%EigenVec(istate,jstate))).GE.huge(real(this%EigenVec(istate,jstate))))ierr=1
          if(abs(aimag(this%EigenVec(istate,jstate))).GE.huge(aimag(this%EigenVec(istate,jstate))))ierr=1
       end do
    end do
    if(ierr.NE.0)then
       call Warn('hs_check: Eigen Vector elements have bad values.')
       hs_check=1
       return
    end if

    !real(double),dimension(:),pointer::EigenVal
    if(.not.associated(this%EigenVal))then
       call Warn('hs_check: Eigen Value memory not associated.')
       hs_check=1
       return
    end if
    if(size(this%EigenVal).NE.this%nstate)then
       call Warn('hs_check: Eigen Value elements must equal number of quantum states.')
       hs_check=1
       return
    end if
    ierr=0
    do istate=1,this%nstate
       if(this%EigenVal(istate).NE.this%EigenVal(istate))ierr=1
       if(abs(this%EigenVal(istate)).GE.huge(this%EigenVal(istate)))ierr=1
    end do
    if(ierr.NE.0)then
       call Warn('hs_check: Eigen Value elements have bad values.')
       hs_check=1
       return
    end if
  end function hs_check

end module hs_class
