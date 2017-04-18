!>\breif
!! Asymmetric vibrating diatomic molecule.
!!\details
!! Describes the vibrational states of a Morse oscillating diatomic molecule. Solutions to the eigen state wave functions are projected onto the internuclear coordinate and are solved analytically using Laguerre polynomials.
!\todo
!! finish check method
!<---------------------------------------------------------------
module MorseOscillator_class
  use type_kinds
  use ErrorLog
  use string
  use atomicunits
  use rand
  use textgraphs
  use filemanager
  use hs_class 
  implicit none
  private

  public::MorseOscillator
  public::new, display, save, update, resample, kill, check

  type MorseOscillator
     logical::initialized=.false.
     type(hs)::hs

     integer(long)::npt
     real(double),dimension(:,:),pointer::wf
     real(double)::De,alpha,q0,mu,qmin,qmax

  end type MorseOscillator

  interface new
     module procedure MorseOscillator_init
  end interface

  interface kill
     module procedure MorseOscillator_kill
  end interface

  interface display
     module procedure MorseOscillator_display
  end interface

  interface save
     module procedure MorseOscillator_save
  end interface

  interface update
     module procedure MorseOscillator_update
  end interface

  interface resample
     module procedure MorseOscillator_resample
  end interface

  interface check
     module procedure MorseOscillator_check
  end interface

contains
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine MorseOscillator_init(this,file)
    type(MorseOscillator),intent(inout)::this
    character*(*),intent(in),optional::file
    character(len=title)::filetype
    character(len=path)::infile

    integer::unit
    logical::usedefault,usedunit

    integer::i,j

    call Note('Begin MorseOscillator_init.')
    if(present(file))call Note('input file= '//file)

    !check if input file is present and valid
    if(present(file))then

       !check if file is there
       if(check(file).EQ.1)&
            call stop('MorseOscillator_init: cannot find input file '//file)

       !assign a unique unit label
       unit=newunit()

       !open the file
       open(unit,file=file)

       !read the file type - should always be on the first line
       read(unit,*)filetype
       filetype=adjustl(filetype)

       !check if input file is the right kind of file
       if(trim(filetype).NE.'MorseOscillator')& 
            call Stop('MorseOscillator_init: input file is not valid.')
    end if

    ! prepare quantum primative
    if(present(file))then
       read(unit,*)infile
       infile=adjustl(infile)
       call new(this%hs,trim(infile))
    else
       call new(this%hs)
    end if

     if(present(file))then
        read(unit,*)this%npt
     else
        write(*,*)'Please enter npt.'
        read(*,*)this%npt
     end if

     if(present(file))then
        read(unit,*)this%De
     else
        write(*,*)'Please enter De.'
        read(*,*)this%De
     end if

     if(present(file))then
        read(unit,*)this%alpha
     else
        write(*,*)'Please enter alpha.'
        read(*,*)this%alpha
     end if

     if(present(file))then
        read(unit,*)this%q0
     else
        write(*,*)'Please enter q0.'
        read(*,*)this%q0
     end if

     if(present(file))then
        read(unit,*)this%mu
     else
        write(*,*)'Please enter mu.'
        read(*,*)this%mu
     end if

     if(present(file))then
        read(unit,*)this%qmin
     else
        write(*,*)'Please enter qmin.'
        read(*,*)this%qmin
     end if

     if(present(file))then
        read(unit,*)this%qmax
     else
        write(*,*)'Please enter qmax.'
        read(*,*)this%qmax
     end if

     if(associated(this%wf))nullify(this%wf)
     allocate(this%wf(this%npt,this%hs%nstate))
     
    !finished reading data - now close input file
    if(present(file))close(unit)

    !declare that initialization is complete
    this%initialized=.true.

    !update (or resample if brand new) this derived type
    call resample(this)

    !do a final check before we exit
    if(check(this).EQ.1)call stop('MorseOscillator_init: failed final check!')

  end subroutine MorseOscillator_init
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine MorseOscillator_display(this,msg)
    type(MorseOscillator),intent(in)::this
    character*(*),intent(in),optional::msg

    !notice these new declarations!
    integer(long)::istate,ipt
    real(double)::grid(this%npt),dq,L

    call Note('Begin MorseOscillator_display.')
    if(check(this).NE.0)then
       call warn('MorseOscillator_display: failed check','displaying nothing.')
    else
       
       !display the hs
       call display(this%hs,msg)

       !****    Display the derived type's attributes here if you want   ****!

       write(*,*)'npt = '//trim(int2str(this%npt))
       write(*,*)'De = '//trim(float2str(this%De))
       write(*,*)'alpha = '//trim(float2str(this%alpha))
       write(*,*)'equillibrium bond distance = '//trim(float2str(this%q0))
       write(*,*)'reduced mass = '//trim(float2str(this%mu))
       write(*,*)'bond domain = ['//trim(float2str(this%qmin))&
            //':'//trim(float2str(this%qmax))//']'

       L=this%qmax-this%qmin
       dq=L/real(this%npt,double)
       do ipt=1,this%npt
          grid(ipt)=(ipt-1)*dq-this%qmin
       end do
       do istate=1,this%hs%nstate
          call display(grid(:),this%wf(:,istate))
       end do

       !*********************************************************************!
       !=====================================================================!
       !*******        Example display scalar attribute 'var'      **********!
       !                                                                     !
       ! write(*,*)'VAR=',this%var                                           !
       !                                                                     !
       !*********************************************************************!
       !=====================================================================!
       !*******    Example display NxM matrix attribute 'matrix'   **********!
       !                                                                     !
       ! write(*,*)'MATRIX=',this%matrix   !a simple example                 !
       !                                                                     !
       ! write(*,*)'MATRIX='               !a better example:                !
       ! do i=1,N                          !write each row on a new line     !
       !    write(*,*)i,(this%matrix(i,j),j=1,M)                             !
       ! end do                                                              !
       !*********************************************************************!


    end if
  end subroutine MorseOscillator_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine MorseOscillator_save(this,file)
    type(MorseOscillator),intent(in)::this
    character*(*),intent(in)::file

    integer(short)::unit
    logical::usedunit      

    call note('Begin MorseOscillator_save.')
    call Note('input file= '//file)
    if(check(this).NE.0)then
       call warn('MorseOscillator_save: failed check.','not saving object.')
    else

       !assign a unique unit label
       unit=newunit()

       !open save file
       open(unit,file=file)

       !always write the data type on the first line
       write(unit,*)'MorseOscillator'

       !save the hs
       call save(this%hs,file//'.hs')

       !write the location of the hs
       write(unit,*)quote(file//'.hs')

       !******      Save below all the derived type's attributes       ******!
       !******         in the order the NEW command reads them         ******!

       write(unit,*)this%npt
       write(unit,*)this%De
       write(unit,*)this%alpha
       write(unit,*)this%q0
       write(unit,*)this%mu
       write(unit,*)this%qmin
       write(unit,*)this%qmax

       !*********************************************************************!
       !=====================================================================!
       !******      Example - Save an attribute called 'var  '    ***********!
       ! write(unit,*)this%var                                               !
       !*********************************************************************!
       !=====================================================================!
       !***  Example - Save an NxM matrix attribute called 'matrix'  ********!
       ! write(unit,*)((this%matrix(i,j),j=1,M),i=1,N)                       !
       !*********************************************************************!


       !finished saving all attributes - now close save file
       close(unit)
    end if
  end subroutine MorseOscillator_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-
  subroutine MorseOscillator_update(this)
    type(MorseOscillator),intent(inout)::this

    call Note('Begin MorseOscillator_update.')

    !update the hs
    call update(this%hs)
    
    !******   Recompute any attribute values that might have evolved   ******!


    !************************************************************************!
    !========================================================================!
    !*****    Example - attribute 'var' is always equall to the trace   *****!
    !                   of the hs's denisity                          !
    !                                                                        !
    ! this%var=0._double                                                     !
    ! do istate=1,this%hs%nstate                                             !
    !    this%var=this%var+this%hs%den(istate,istate)                        !
    ! end do                                                                 !
    !                                                                        !
    !************************************************************************!

    !recompute eigen states
    this%hs%EigenVec=cmplx(this%hs%diabat)
    call diagonalize(this%hs%nstate,this%hs%EigenVec,this%hs%EigenVal)

    !do a final check before we exit (usually left out for performance boost)
    !if(check(this).EQ.1)call stop('MorseOscillator_update: failed final check!')

  end subroutine MorseOscillator_update
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-
  subroutine MorseOscillator_resample(this)
    type(MorseOscillator),intent(inout)::this
    integer(long)::istate,jstate,ipt
    real(double)::dq,L

    call Note('Begin MorseOscillator_resample.')

    !resample the hs
    call resample(this%hs)

    !****  Resample any attributes to satisfy re-initiation conditions   ****! 

    !call function that recomputes wavefunctions
    L=this%qmax-this%qmin
    dq=L/real(this%npt,double)
    do istate=1,this%hs%nstate
       do ipt=1,this%npt
          this%wf(ipt,istate)=cos(((ipt-1)*dq+this%q0)*istate*twopi/L)
       end do
    end do


    !************************************************************************!
    !========================================================================!
    !********    Example - attribute 'var' is always initially a     ********!
    !                      Gaussian random number                            !
    !                                                                        !
    ! this%var=gran()                                                        !
    !                                                                        !
    !************************************************************************!

    !update now that we have changed some values
    call update(this)    


  !if (check(this))call warn('hey this has failed the check!')

  end subroutine MorseOscillator_resample
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine MorseOscillator_kill(this)
    type(MorseOscillator),intent(inout)::this
 
    call note('Begin MorseOscillator_kill.')

    !kill the hs
    call kill(this%hs)


    !*******************       Nullify all pointers    **********************!
 
     if(associated(this%wf))nullify(this%wf)

    !************************************************************************!
    !========================================================================!
    !******        Example - cleanup matrix attribute 'matrix'       ********!
    !                                                                        !
    ! if(associated(this%matrix))nullify(this%matrix)                        !
    !                                                                        !
    !************************************************************************!


    !Decalare this derived type un-initialized
    this%initialized=.false.

  end subroutine MorseOscillator_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short)function MorseOscillator_check(this)
    type(MorseOscillator),intent(in)::this

    call Note('Checking MorseOscillator.')

    !initiate with no problems found 
    MorseOscillator_check=0

    !check that this derived type is initialized
    if(.not.this%initialized)then
       call Warn('MorseOscillator_check: MorseOscillator object not initialized.')
       MorseOscillator_check=1
       return
    end if

    !check the hs
    if(check(this%hs).EQ.1)call stop('MorseOscillator_check: failed hs check!')


    !********    Check all attributes are within acceptable values    *******!

    if(this%npt.NE.this%npt)then
       MorseOscillator_check=1
       call Warn('MorseOscillator_check: npt is NAN.')
       return
    end if

    if(this%npt.GT.Huge(this%npt))then
       MorseOscillator_check=1
       call Warn('MorseOscillator_check: npt is Huge.')
       return
    end if

    if (this%npt.LT.1)then 
       MorseOscillator_check=1
       call Warn('MorseOscillator_check: npt cannot be less than 1.')
       return
    end if


    if(this%De.NE.this%De)then
       MorseOscillator_check=1
       call Warn('MorseOscillator_check: De is NAN.')
       return
    end if

    if(this%De.GE.huge(this%De))then
       MorseOscillator_check=1
       call Warn('MorseOscillator_check: De is Huge.')
       return
    end if

    if(this%De.LT.epsilon(this%De))then
       MorseOscillator_check=1
       call Warn('MorseOscillator_check: De is tiny.')
       return
    end if


    ! real(double),dimension(:,:),pointer::wf
    ! real(double)::De,alpha,q0,mu,qmin,qmax



    !**************************************************************************************!
    !======================================================================================!
    !**********     Example - check an integer attribute 'ndim'    ************************!
    !                                                                                      !
    ! !check if integer 'ndim' is NAN (not a number)                                       !
    ! if(this%ndim.NE.this%ndim)then                                                       !
    !    call Warn('MorseOscillator_check: ndim not a number.')                   !
    !    MorseOscillator_check=1                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    ! !check if 'ndim' is too big to fit in its memory                                     !
    ! if(abs(this%ndim).GE.huge(this%ndim))then                                            !
    !    call Warn('MorseOscillator_check: ndim is too big.')                     !
    !    MorseOscillator_check=1                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    ! !add a constrain that says 'ndim' can only be positive                               !
    ! if(this%ndim.LE.0)then                                                               !
    !    call Warn('MorseOscillator_check: ndim not a positive integer.')         !
    !    MorseOscillator_check=1                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    !**************************************************************************************!
    !======================================================================================!
    !**********    Example - check a real number attribute 'var'   ************************!
    !                                                                                      !
    ! !check if 'var' is not a number                                                      !
    ! if(this%var.NE.this%var)then                                                         !
    !    call Warn('MorseOscillator_check: var is not a number.')                 !
    !    MorseOscillator_check=1                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    ! !check if 'var' is too big to fit in its memory                                      !
    ! if(abs(this%var).GE.huge(this%var))then                                              !
    !    call Warn('MorseOscillator_check: var is too big.')                      !
    !    MorseOscillator_check=1                                                  !
    !   return                                                                             !
    ! end if                                                                               !
    !                                                                                      !
    ! !add a constrain that says 'var' can not be zero:                                    !
    ! !      'var' can not be smaller than the smallest computable value                   !
    ! if(abs(this%var).LE.epsilon(this%var))then                                           !
    !    call Warn('MorseOscillator_check: var is too small.')                    !
    !    MorseOscillator_check=1                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    !**************************************************************************************!
    !========================================================================!
    !*********    Example - check an NxM matrix attribute 'matrix'   ********!
    !                                                                        !
    ! !check that 'matrix' points to something                               !
    ! if(.not.associated(this%matrix))then                                   !
    !    call Warn('MorseOscillator_check: matrix memory not associated.')        !
    !    MorseOscillator_check=1                                                  !
    !    return                                                              !
    ! end if                                                                 !
    !                                                                        !
    ! !check that 'matrix' has the right dimensions                          !
    ! if(size(this%matrix).NE.N*M)then                                       !
    !    call Warn('MorseOscillator_check: number of matrix elements not = N*M.') !
    !    MorseOscillator_check=1                                                  !
    !    return                                                              !
    ! end if                                                                 !
    !                                                                        !
    ! !check for NAN values in the matrix                                    !
    ! if(any(this%matrix.NE.this%matrix))then                                !
    !    call Warn('MorseOscillator_check: matirx has NAN values.')               !
    !    MorseOscillator_check=1                                                  !
    !    return                                                              !
    ! end if                                                                 !
    !                                                                        !
    ! !check if any matrix element values are too big for thier memory       !
    ! if(any(abs(this%matirx).GT.huge(this%matirx)))then                     !
    !    call Warn('MorseOscillator_check: matrix has huge values.')              !
    !    mappingH_check=1                                                    !
    !    MorseOscillator_check=1                                                  !
    !    return                                                              !
    ! end if                                                                 !
    !                                                                        !
    !************************************************************************!

  end function MorseOscillator_check

end module MorseOscillator_class

