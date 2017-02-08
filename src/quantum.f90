!---------------------------------------------------
!> \brief
!! Generic quantum Class
!! \details
!! Can assume any derived quantum type
!! \authors
!! Daniel Montemayor
!<---------------------------------------------------
module quantum_class
  use type_kinds
  use filemanager
  use ErrorLog
  use string
  use hs_class
  use dipoles_class
  implicit none
  private
  public::quantum
  public::new,kill,update,resample,display,save,check
  type quantum
     logical::initialized=.false.
     type(hs),pointer::hs
     type(hs)::primitive
     character(len=title)::type=''
     type(dipoles)::dipoles 
  end type quantum
  interface new
     module procedure quantum_new
  end interface
  interface kill
     module procedure quantum_kill
  end interface
  interface update
     module procedure quantum_update
  end interface
  interface resample
     module procedure quantum_resample
  end interface
  interface display
     module procedure quantum_display
  end interface
  interface save
     module procedure quantum_save
  end interface
  interface check
     module procedure quantum_check
  end interface
contains
  subroutine quantum_new(this,type,file)
    type(quantum),intent(inout),target::this    
    character*(*),intent(in),optional::type
    character*(*),intent(in),optional::file
    character(len=title)::filetype
    character(len=path)::infile
    integer(long)::unit
    call Note('Begin quantum_new.')
    filetype='unknown'
    infile='unknown'
    if (present(type))call note('input type= '//type)
    if (present(file))call note('input file= '//file)
    call quantum_kill(this)
    if(.not.(present(file).or.present(type)))&
         call stop('quantum_new: niether quantum input type nor input file was found.')
    if(present(type))this%type=type
    if(present(file))then
       if(present(type))then
          infile=file
          if(check(trim(infile)))call stop('quantum_new: input file not found.')
       else
          if(check(file))call stop('quantum_new: input file not found.')
          unit=newunit()
          open(unit,file=file)
          read(unit,*)filetype
          filetype=adjustl(filetype)
          if(trim(filetype).NE.'quantum')call stop('quantum_new: input file is not valid.')
          read(unit,*)this%type
          this%type=adjustl(this%type)
          call Note('type found in input file= '//trim(this%type))
          read(unit,*)infile
          infile=adjustl(infile)
          if(check(trim(infile)))call stop('quantum_new: input file not found.')
          close(unit)
       end if
    end if
    select case(trim(this%type))
    case('hs')
       if(trim(infile).NE.'unknown')then;call new(this%primitive,trim(infile));else;call new(this%primitive);end if
       this%hs => this%primitive
    case('dipoles')
       if(trim(infile).NE.'unknown')then;call new(this%dipoles,trim(infile));else;call new(this%dipoles);end if
       this%hs => this%dipoles%hs
    case default
       call stop('quantum_new: unknown quantum type.')
    end select
    this%initialized=.true.    
  end subroutine quantum_new
  subroutine quantum_kill(this)
    type(quantum),intent(inout)::this
    call Note('Begin quantum_kill.')
    nullify(this%hs)
    call kill(this%primitive)
    call kill(this%dipoles)
    this%type=''
    this%initialized=.false.
  end subroutine quantum_kill
  subroutine quantum_update(this)
    type(quantum),intent(inout)::this
    call Note('Begin quantum_update.')
    !if(check(this))call stop('quantum_update: failed object check.')
    select case(trim(this%type))
    case('hs');call update(this%primitive)
    case('dipoles');call update(this%dipoles)
    case default
       call stop('quantum_update: unknown quantum type.')
    end select
  end subroutine quantum_update
  subroutine quantum_resample(this)
    type(quantum),intent(inout)::this    
    call Note('Begin quantum_resample.')
    if(check(this))call stop('quantum_resample: failed object check.')
    select case(trim(this%type))
    case('hs');call resample(this%primitive)
    case('dipoles');call resample(this%dipoles)
    case default
       call stop('quantum_update: unknown quantum type.')
    end select
  end subroutine quantum_resample
  subroutine quantum_display(this,msg)
    type(quantum),intent(in)::this    
    character*(*),intent(in),optional::msg    
    logical::withmsg
    call Note('Begin quantum_display.')
    if(check(this))then
       call warn('quantum_display, failed object check','displaying nothing.')
       return
    end if
    withmsg=.false.
    if(present(msg))withmsg=.true.
    select case(trim(this%type))
    case('hs')
       if(withmsg)call display(this%hs,msg)
       if(.not.withmsg)call display(this%hs)
    case('dipoles')
       if(withmsg)call display(this%dipoles,msg)
       if(.not.withmsg)call display(this%dipoles)
    case default
       call stop('quantum_display: unknown quantum type.')
    end select
  end subroutine quantum_display
  subroutine quantum_save(this,file)
    type(quantum),intent(in)::this
    character*(*),intent(in)::file    
    integer(long)::unit
    call Note('Begin quantum_save.')
    call Note('input file= '//file)
    if(check(this))call stop('quantum_save: failed object check.')
    unit=newunit()
    open(unit,file=file)
    write(unit,*)'quantum'
    write(unit,*)trim(this%type)
    write(unit,*)quote(file//'.qs')
    select case(trim(this%type))
    case('hs');call save(this%primitive,file//'.qs')
    case('dipoles');call save(this%dipoles,file//'.qs')
    case default
       call stop('quantum_update: unknown quantum type.')
    end select
  end subroutine quantum_save
  integer(short)function quantum_check(this)
    type(quantum),intent(in)::this
    call Note('Checking quantum object.')
    quantum_check=0
    if(.not.this%initialized)then
       call Warn('quantum_check: quantum not initialized.')
       quantum_check=1
       return
    end if
    if(.not.associated(this%hs))then
       call Warn('quantum_check: quantum primitive not associated.')
       quantum_check=1
       return
    end if
    if(check(this%hs))then
       call Warn('quantum_check: quantum primitive failed check.')
       quantum_check=1
       return
    end if
    select case(trim(this%type))
    case('hs')
    case('dipoles')
    case default
       call stop('quantum_update: unknown quantum type: '//trim(this%type))
       quantum_check=1
       return
    end select
  end function quantum_check
end module quantum_class
