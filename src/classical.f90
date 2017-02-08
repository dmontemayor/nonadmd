!---------------------------------------------------
!> \brief
!! Generic classical Class
!! \details
!! Can assume any derived classical type
!! \authors
!! Daniel Montemayor
!<---------------------------------------------------
module classical_class
  use type_kinds
  use filemanager
  use ErrorLog
  use string
  use hb_class
  use harmonicbath_class
  implicit none
  private
  public::classical
  public::new,kill,update,resample,display,save,check
  type classical
     logical::initialized=.false.
     type(hb),pointer::hb
     type(hb)::primitive
     character(len=title)::type=''
     type(harmonicbath)::harmonicbath 
  end type classical
  interface new
     module procedure classical_new
  end interface
  interface kill
     module procedure classical_kill
  end interface
  interface update
     module procedure classical_update
  end interface
  interface resample
     module procedure classical_resample
  end interface
  interface display
     module procedure classical_display
  end interface
  interface save
     module procedure classical_save
  end interface
  interface check
     module procedure classical_check
  end interface
contains
  subroutine classical_new(this,type,file)
    type(classical),intent(inout),target::this    
    character*(*),intent(in),optional::type
    character*(*),intent(in),optional::file
    character(len=title)::filetype
    character(len=path)::infile
    integer(long)::unit
    call Note('Begin classical_new.')
    filetype='unknown'
    infile='unknown'
    if (present(type))call note('input type= '//type)
    if (present(file))call note('input file= '//file)
    call classical_kill(this)
    if(.not.(present(file).or.present(type)))&
         call stop('classical_new: niether classical input type nor input file was found.')
    if(present(type))this%type=type
    if(present(file))then
       if(present(type))then
          infile=file
          if(check(trim(infile)))call stop('classical_new: input file not found.')
       else
          if(check(file))call stop('classical_new: input file not found.')
          unit=newunit()
          open(unit,file=file)
          read(unit,*)filetype
          filetype=adjustl(filetype)
          if(trim(filetype).NE.'classical')call stop('classical_new: input file is not valid.')
          read(unit,*)this%type
          this%type=adjustl(this%type)
          call Note('type found in input file= '//trim(this%type))
          read(unit,*)infile
          infile=adjustl(infile)
          if(check(trim(infile)))call stop('classical_new: input file not found.')
          close(unit)
       end if
    end if
    select case(trim(this%type))
    case('hb')
       if(trim(infile).NE.'unknown')then;call new(this%primitive,trim(infile));else;call new(this%primitive);end if
       this%hb => this%primitive
    case('harmonicbath')
       if(trim(infile).NE.'unknown')then;call new(this%harmonicbath,trim(infile));else;call new(this%harmonicbath);end if
       this%hb => this%harmonicbath%hb
    case default
       call stop('classical_new: unknown classical type.')
    end select
    this%initialized=.true.    
  end subroutine classical_new
  subroutine classical_kill(this)
    type(classical),intent(inout)::this
    call Note('Begin classical_kill.')
    nullify(this%hb)
    call kill(this%primitive)
    call kill(this%harmonicbath)
    this%type=''
    this%initialized=.false.
  end subroutine classical_kill
  subroutine classical_update(this)
    type(classical),intent(inout)::this
    call Note('Begin classical_update.')
    !if(check(this))call stop('classical_update: failed object check.')
    select case(trim(this%type))
    case('hb');call update(this%primitive)
    case('harmonicbath');call update(this%harmonicbath)
    case default
       call stop('classical_update: unknown classical type.')
    end select
  end subroutine classical_update
  subroutine classical_resample(this)
    type(classical),intent(inout)::this    
    call Note('Begin classical_resample.')
    if(check(this))call stop('classical_resample: failed object check.')
    select case(trim(this%type))
    case('hb');call resample(this%primitive)
    case('harmonicbath');call resample(this%harmonicbath)
    case default
       call stop('classical_update: unknown classical type.')
    end select
  end subroutine classical_resample
  subroutine classical_display(this,msg)
    type(classical),intent(in)::this    
    character*(*),intent(in),optional::msg    
    logical::withmsg
    call Note('Begin classical_display.')
    if(check(this))then
       call warn('classical_display, failed object check','displaying nothing.')
       return
    end if
    withmsg=.false.
    if(present(msg))withmsg=.true.
    select case(trim(this%type))
    case('hb')
       if(withmsg)call display(this%hb,msg)
       if(.not.withmsg)call display(this%hb)
    case('harmonicbath')
       if(withmsg)call display(this%harmonicbath,msg)
       if(.not.withmsg)call display(this%harmonicbath)
    case default
       call stop('classical_display: unknown classical type.')
    end select
  end subroutine classical_display
  subroutine classical_save(this,file)
    type(classical),intent(in)::this
    character*(*),intent(in)::file    
    integer(long)::unit
    call Note('Begin classical_save.')
    call Note('input file= '//file)
    if(check(this))call stop('classical_save: failed object check.')
    unit=newunit()
    open(unit,file=file)
    write(unit,*)'classical'
    write(unit,*)trim(this%type)
    write(unit,*)quote(file//'.cs')
    select case(trim(this%type))
    case('hb');call save(this%primitive,file//'.cs')
    case('harmonicbath');call save(this%harmonicbath,file//'.cs')
    case default
       call stop('classical_update: unknown classical type.')
    end select
  end subroutine classical_save
  integer(short)function classical_check(this)
    type(classical),intent(in)::this
    call Note('Checking classical object.')
    classical_check=0
    if(.not.this%initialized)then
       call Warn('classical_check: classical not initialized.')
       classical_check=1
       return
    end if
    if(.not.associated(this%hb))then
       call Warn('classical_check: classical primitive not associated.')
       classical_check=1
       return
    end if
    if(check(this%hb))then
       call Warn('classical_check: classical primitive failed check.')
       classical_check=1
       return
    end if
    select case(trim(this%type))
    case('hb')
    case('harmonicbath')
    case default
       call stop('classical_update: unknown classical type: '//trim(this%type))
       classical_check=1
       return
    end select
  end function classical_check
end module classical_class
