!>\brief
!! Derived Subsystem Template
!!\details
!! Use this template to help you start building a derived quantum, classical, or coupling subsystem. 
!! Replace any instance of \a primitive with \a hs, \a hb, or \a hc for your derived
!! quantum, classical or coupling subsystem respectively.
!! Do not remove the \a new, \a display, \a save, \a update, \a resample, \a kill, and \a check subroutines and functions.
!! These methods must be there for your derived subsystem to integrate properly with the application.
!! You can leave these methods blank if you dont want those particular funcitionalities, although, this is not advised.
!! Follow the examples commented in the source code to define the various attributes of your derived subsystem.
!! DO NOT OVERWRITE THIS FILE. Make a copy and title it \a myderivedsusbsystem.f90 where 'myderivedsubsystem' is your choice
!! one-word title for your derived subsystem.
!! Finally, Replace instances of the string 'TEMPLATE' in this template with the same
!! one-word title 'myderivedsubsystem' you chose.
!<------------------------------------------------------------------------
module TEMPLATE_class
  use type_kinds
  use ErrorLog
  use string
  use atomicunits
  use rand
  use textgraphs
  use filemanager
  use outputdisplay
  use primitive_class !<replace 'primitive' for hs, hb, or hc for your derived quantum, classical or coupling subsystem respectively 
  implicit none
  private

  public::TEMPLATE
  public::new, display, save, update, resample, kill, check

  type TEMPLATE
     logical::initialized=.false.
     type(primitive)::primitive
     !**********     Enter your derived type's attributes here     **********!



     !***********************************************************************!
     !=======================================================================!
     !***   Here are some example attributes your derived type may have   ***!
     ! type(primitive)::primitive2                     !an extra primitive   !
     ! integer(short)::label                           !a short integer      !
     ! integer(long)::ndim                             !a long integer       !
     ! real(double)::var                               !a real variable      !
     ! complex(double)::zed                            !a complex variable   !
     ! real(double),dimension(:,:),pointer::matrix     !a real matrix        !
     ! complex(double),dimension(:,:),pointer::Zmatrix !a complex matrix     !
     !***********************************************************************!
  end type TEMPLATE

  !> Creates the TEMPLATE object.
  interface new
     module procedure TEMPLATE_init
  end interface

  !> Destroys the TEMPLATE object.
  interface kill
     module procedure TEMPLATE_kill
  end interface

  !> Displays the current state of the TEMPLATE object.
  interface display
     module procedure TEMPLATE_display
  end interface

  !> Saves the current state of the TEMPLATE object.
  interface save
     module procedure TEMPLATE_save
  end interface

  !> Recaluclates the TEMPLATE object.
  interface update
     module procedure TEMPLATE_update
  end interface

  !> Reinitializes the TEMPLATE object.
  interface resample
     module procedure TEMPLATE_resample
  end interface

  !> Checks that the TEMPLATE object.
  interface check
     module procedure TEMPLATE_check
  end interface

contains

  !======================================================================
  !> \brief Creates and initializes the TEMPLATE object.
  !> \param this is the TEMPLATE object to be initialized.
  !> \param[in] file is an optional string containing the name of a previously saved TEMPLATE file.
  !> \remark If no input file is provided the user must manually initialize THIS using stout.
  !=====================================================================
  subroutine TEMPLATE_init(this,file)
    type(TEMPLATE),intent(inout)::this
    character*(*),intent(in),optional::file
    character(len=title)::filetype
    character(len=path)::infile

    integer::unit
    logical::usedefault,usedunit

    call Note('Begin TEMPLATE_init.')

    !check if input file is present and valid
    if(present(file))then

       !check if file is there
       if(check(file).EQ.1)call stop('TEMPLATE_init: cannot find input file '//file)

       !assign a unique unit label
       unit=newunit()

       !open the file
       open(unit,file=file)

       !read the file type - should always be on the first line
       read(unit,*)filetype
       filetype=adjustl(filetype)

       !check if input file is the right kind of file
       if(trim(filetype).NE.'TEMPLATE')& 
            call Stop('TEMPLATE_init: input file is not valid.')

    end if
    
    ! prepare primative
    if(present(file))then
       read(unit,*)infile
       infile=adjustl(infile)
       call new(this%primitive,trim(infile))
    else
       call new(this%primitive)
    end if

    !******     Initiate all your derived type's attributes below      ******!
    !************************************************************************!




    !************************************************************************!
    !========================================================================!
    !********           Example - Setup attribute 'var'            **********!
    !                                                                        !
    ! if(present(file))then                                                  !
    !    read(unit,*)this%var          !Read from input file if present      !
    ! else                                                                   !
    !    write(*,*)'Please enter VAR.' !manually enter information           !
    !    read(*,*)this%var                                                   !
    ! end if                                                                 !
    !                                                                        !
    !************************************************************************!
    !========================================================================!
    !********   Example - Setup a matrix attribute called 'matrix'   ********!
    !                                                                        !
    ! if(associated(this%matrix))nullify(this%matrix) !cleanup memory        !
    ! allocate(this%matirx(N,M))                      !create an NxM matrix  !
    !                                                                        !
    ! !Now fill that matrix with data from the input file or manually        !
    ! if(present(file))then                                                  !
    !    read(unit,*)((this%matrix(i,j),j=1,M),i=1,N) !read M loop first     !
    ! else                                                                   !
    !    write(*,*)'Enter value of matrix element:'                          !
    !    do i=1,N                                                            !
    !       do j=1,M                                                         !
    !          write(*,*)'index=',i,j      !prompt user the matrix index     !
    !          read(*,*)this%matrix(i,j)   !read data manually from keyboard !
    !       end do                                                           !
    !    end do                                                              !
    ! end if                                                                 !
    !                                                                        !
    !************************************************************************!


    !finished reading data - now close input file
    if(present(file))close(unit)

    !declare that initialization is complete
    this%initialized=.true.

    !update (or resample if brand new) this derived type
    if(present(file))then
       call update(this)
    else
       call resample(this)
    end if

    !do a final check before we exit
    if(check(file).EQ.1)call stop('TEMPLATE_init: failed final check!')

  end subroutine TEMPLATE_init

  !======================================================================
  !> \brief Destroys the TEMPLATE object.
  !> \param this is the TEMPLATE object to be destroyed.
  !====================================================================
  subroutine TEMPLATE_kill(this)
    type(TEMPLATE),intent(inout)::this
 
    call note('Begin TEMPLATE_kill.')

    !kill the primitive
    call kill(this%primitive)


    !*******************       Nullify all pointers    **********************!

 



    !************************************************************************!
    !========================================================================!
    !******        Example - cleanup matrix attribute 'matrix'       ********!
    !                                                                        !
    ! if(associated(this%matrix))nullify(this%matrix)                        !
    !                                                                        !
    !************************************************************************!


    !Decalare this derived type un-initialized
    this%initialized=.false.

  end subroutine TEMPLATE_kill

  !======================================================================
  !> \brief Computes the current state of TEMPLATE object.
  !> \param this is the TEMPLATE  object to be updated.
  !======================================================================
  subroutine TEMPLATE_update(this)
    type(TEMPLATE),intent(inout)::this

    call Note('Begin TEMPLATE_update.')
    
    !Primitives usually dont get updated
    !    call update(this%primitive)    


    !******   Recompute any attribute values that might have evolved   ******!





    !************************************************************************!
    !========================================================================!
    !*****    Example - attribute 'var' is always equall to the trace   *****!
    !                   of the primitive's denisity                          !
    !                                                                        !
    ! this%var=0._double                                                     !
    ! do istate=1,this%primitive%nstate                                             !
    !    this%var=this%var+this%primitive%den(istate,istate)                        !
    ! end do                                                                 !
    !                                                                        !
    !************************************************************************!

    !Usually leave out final check before we exit
    !if(check(this).EQ.1)call stop('TEMPLATE_update: failed final check!')

  end subroutine TEMPLATE_update

  !======================================================================
  !> \brief Re-initiallizes the TEMPLATE object.
  !> \param this is the TEMPLATE  object to be re-initialized.
  !======================================================================
  subroutine TEMPLATE_resample(this)
    type(TEMPLATE),intent(inout)::this
    integer(long)::istate,jstate

    call Note('Begin TEMPLATE_resample.')

    !resample the primitive
    call resample(this%primitive)

    !****  Resample any attributes to satisfy re-initiation conditions   ****! 





    !************************************************************************!
    !========================================================================!
    !********    Example - attribute 'var' is always initially a     ********!
    !                      Gaussian random number                            !
    !                                                                        !
    ! this%var=gran()                                                        !
    !                                                                        !
    !************************************************************************!

    !update now that we have changed some values and do a final check
    call update(this)
    if(check(this).EQ.1)call stop('TEMPLATE_resample: failed final check!')

  end subroutine TEMPLATE_resample

  !======================================================================
  !> \brief Saves the current state of the TEMPLATE object to file.
  !> \param[in] this is the TEMPLATE  object to be updated.
  !> \param[in] file is a string containing the location of the save file.
  !======================================================================
  subroutine TEMPLATE_save(this,file)
    type(TEMPLATE),intent(in)::this
    character*(*),intent(in)::file

    integer(short)::unit
    logical::usedunit      

    call note('Begin TEMPLATE_save.')
    call Note('input file= '//file)
    if(check(this).NE.0)then
       call warn('TEMPLATE_save: failed check.','not saving object.')
    else

       !assign a unique unit label
       unit=newunit()

       !open save file
       open(unit,file=file)

       !always write the data type on the first line
       write(unit,*)'TEMPLATE'

       !save the primitive
       call save(this%primitive,file//'.primitive')

       !write the location of the primitive
       write(unit,*)quote(file//'.primitive')

       !******      Save below all the derived type's attributes       ******!
       !******         in the order the NEW command reads them         ******!





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
  end subroutine TEMPLATE_save


  !======================================================================
  !> \brief Checks the TEMPLATE object.
  !> \param[in] this is the TEMPLATE object to be updated.
  !> \param[in] msg is an optional string message to preface the displayed object.
  !======================================================================
  subroutine TEMPLATE_display(this,msg)
    type(TEMPLATE),intent(in)::this
    character*(*),intent(in),optional::msg

    call Note('Begin TEMPLATE_display.')

    if(check(this).NE.0)then
       call warn('TEMPLATE_display: failed check','displaying nothing.')
       return
    end if

    write(Dunit,*)'____________________________________________________'
    write(Dunit,*)'-------------------   TEMPLATE   -------------------'
    if(present(msg))write(Dunit,*)msg
    write(Dunit,*)'____________________________________________________'
    
    
    !****    Display the derived type's attributes here if you want   ****!
    
    
    
    !*********************************************************************!
    !=====================================================================!
    !*******        Example display scalar attribute 'var'      **********!
    !                                                                     !
    ! write(Dunit,*)'VAR=',this%var                                       !
    !                                                                     !
    !*********************************************************************!
    !=====================================================================!
    !*******    Example display NxM matrix attribute 'matrix'   **********!
    !                                                                     !
    ! write(Dunit,*)'MATRIX=',this%matrix   !a simple example             !
    !                                                                     !
    ! write(Dunit,*)'MATRIX='               !a better example:            !
    ! do i=1,N                          !write each row on a new line     !
    !    write(Dunit,*)i,(this%matrix(i,j),j=1,M)                         !
    ! end do                                                              !
    !*********************************************************************!

    call display(this%primitive,msg='TEMPLATE primitive')
    write(Dunit,*)'===================================================='

  end subroutine TEMPLATE_display

  !======================================================================
  !> \brief Displays the TEMPLATE object.
  !> \param[in] this is the TEMPLATE object to be checked.
  !> \return Nothing if all checks pass or 1 and a warn for the first failed check.
  !> \remark Will exit after first failed check.
  !======================================================================
  integer(short)function TEMPLATE_check(this)
    type(TEMPLATE),intent(in)::this

    call Note('Checking TEMPLATE.')

    !initiate with no problems found 
    TEMPLATE_check=0

    !check that this derived type is initialized
    if(.not.this%initialized)then
       call Warn('TEMPLATE_check: TEMPLATE object not initialized.')
       TEMPLATE_check=1
       return
    end if

    !check the primitive
    if(check(this%primitive))call stop('TEMPLATE_check: failed primitive check!')


    !********    Check all attributes are within acceptable values    *******!





    !**************************************************************************************!
    !======================================================================================!
    !**********     Example - check an integer attribute 'ndim'    ************************!
    !                                                                                      !
    ! !check if integer 'ndim' is NAN (not a number)                                       !
    ! if(this%ndim.NE.this%ndim)then                                                       !
    !    call Warn('TEMPLATE_check: ndim not a number.')                                   !
    !    TEMPLATE_check=1                                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    ! !check if 'ndim' is too big to fit in its memory                                     !
    ! if(abs(this%ndim).GE.huge(this%ndim))then                                            !
    !    call Warn('TEMPLATE_check: ndim is too big.')                                     !
    !    TEMPLATE_check=1                                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    ! !add a constrain that says 'ndim' can only be positive                               !
    ! if(this%ndim.LE.0)then                                                               !
    !    call Warn('TEMPLATE_check: ndim not a positive integer.')                         !
    !    TEMPLATE_check=1                                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    !**************************************************************************************!
    !======================================================================================!
    !**********    Example - check a real number attribute 'var'   ************************!
    !                                                                                      !
    ! !check if 'var' is not a number                                                      !
    ! if(this%var.NE.this%var)then                                                         !
    !    call Warn('TEMPLATE_check: var is not a number.')                                 !
    !    TEMPLATE_check=1                                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    ! !check if 'var' is too big to fit in its memory                                      !
    ! if(abs(this%var).GE.huge(this%var))then                                              !
    !    call Warn('TEMPLATE_check: var is too big.')                                      !
    !    TEMPLATE_check=1                                                                  !
    !   return                                                                             !
    ! end if                                                                               !
    !                                                                                      !
    ! !add a constrain that says 'var' can not be zero:                                    !
    ! !      'var' can not be smaller than the smallest computable value                   !
    ! if(abs(this%var).LE.epsilon(this%var))then                                           !
    !    call Warn('TEMPLATE_check: var is too small.')                                    !
    !    TEMPLATE_check=1                                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    !**************************************************************************************!
    !======================================================================================!
    !*********          Example - check an NxM matrix attribute 'matrix'        ***********!
    !                                                                                      !
    ! !check that 'matrix' points to something                                             !
    ! if(.not.associated(this%matrix))then                                                 !
    !    call Warn('TEMPLATE_check: matrix memory not associated.')                        !
    !    TEMPLATE_check=1                                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    ! !check that 'matrix' has the right dimensions                                        !
    ! if(size(this%matrix).NE.N*M)then                                                     !
    !    call Warn('TEMPLATE_check: number of matrix elements not = N*M.')                 !
    !    TEMPLATE_check=1                                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    ! !check for NAN values in the matrix                                                  !
    ! if(any(this%matrix.NE.this%matrix))then                                              !
    !    call Warn('TEMPLATE_check: matirx has NAN values.')                               !
    !    TEMPLATE_check=1                                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    ! !check if any matrix element values are too big for thier memory                     !
    ! if(any(abs(this%matirx).GT.huge(this%matirx)))then                                   !
    !    call Warn('TEMPLATE_check: matrix has huge values.')                              !
    !    mappingH_check=1                                                                  !
    !    TEMPLATE_check=1                                                                  !
    !    return                                                                            !
    ! end if                                                                               !
    !                                                                                      !
    !**************************************************************************************!

  end function TEMPLATE_check

end module TEMPLATE_class

