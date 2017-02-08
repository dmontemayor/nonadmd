!---------------------------------------------------
!> \brief
!! Generic coupling Class
!! \details
!! Can assume any derived coupling type
!! \authors
!! Daniel Montemayor
!<---------------------------------------------------
module coupling_class
  use type_kinds
  use filemanager
  use ErrorLog
  use string
  use hc_class
  use bilinear_class
  use quantum_class
  use classical_class
  implicit none
  private
  public::coupling
  public::new,kill,update,resample,display,save,check
  type coupling
     logical::initialized=.false.
     type(hc),pointer::hc
     type(hc)::primitive
     character(len=title)::type=''
     type(bilinear)::bilinear 
  end type coupling
  interface new
     module procedure coupling_new
  end interface
  interface kill
     module procedure coupling_kill
  end interface
  interface update
     module procedure coupling_update
  end interface
  interface resample
     module procedure coupling_resample
  end interface
  interface display
     module procedure coupling_display
  end interface
  interface save
     module procedure coupling_save
  end interface
  interface check
     module procedure coupling_check
  end interface
contains
  subroutine coupling_new(this,qs,cs,type,file)
    type(coupling),intent(inout),target::this    
    type(quantum),intent(inout),target::qs    
    type(classical),intent(inout),target::cs    
    character*(*),intent(in),optional::type
    character*(*),intent(in),optional::file
    character(len=title)::filetype
    character(len=path)::infile
    integer(short)::unit
    call Note('Begin coupling_new.')
    filetype='unknown'
    infile='unknown'
    if (present(type))call note('input type= '//type)
    if (present(file))call note('input file= '//file)
    call coupling_kill(this)
    if(.not.(present(file).or.present(type)))&
         call stop('coupling_new: niether coupling input type nor input file was found.')
    if(present(type))this%type=type
    if(present(file))then
       if(present(type))then
          infile=file
          if(check(trim(infile)))call stop('coupling_new: input file not found.')
       else
          if(check(file))call stop('coupling_new: input file not found.')
          unit=newunit()
          open(unit,file=file)
          read(unit,*)filetype
          filetype=adjustl(filetype)
          if(trim(filetype).NE.'coupling')call stop('coupling_new: input file is not valid.')
          read(unit,*)this%type
          this%type=adjustl(this%type)
          call Note('type found in input file= '//trim(this%type))
          read(unit,*)infile
          infile=adjustl(infile)
          if(check(trim(infile)))call stop('coupling_new: input file not found.')
          close(unit)
       end if
    end if
    select case(trim(this%type))
    case('hc')
       if(trim(infile).NE.'unknown')then;call new(this%primitive,qs,cs,trim(infile));else;call new(this%primitive,qs,cs);end if
       this%hc => this%primitive
    case('bilinear')
       if(trim(infile).NE.'unknown')then;call new(this%bilinear,qs,cs,trim(infile));else;call new(this%bilinear,qs,cs);end if
       this%hc => this%bilinear%hc
    case default
       call stop('coupling_new: unknown coupling type.')
    end select
    this%initialized=.true.    
  end subroutine coupling_new
  subroutine coupling_kill(this)
    type(coupling),intent(inout)::this
    call Note('Begin coupling_kill.')
    nullify(this%hc)
    call kill(this%primitive)
    call kill(this%bilinear)
    this%type=''
    this%initialized=.false.
  end subroutine coupling_kill
  subroutine coupling_update(this)
    type(coupling),intent(inout)::this
    call Note('Begin coupling_update.')
    !if(check(this))call stop('coupling_update: failed object check.')
    select case(trim(this%type))
    case('hc');call update(this%primitive)
    case('bilinear');call update(this%bilinear)
    case default
       call stop('coupling_update: unknown coupling type.')
    end select
  end subroutine coupling_update
  subroutine coupling_resample(this)
    type(coupling),intent(inout)::this    
    call Note('Begin coupling_resample.')
    if(check(this))call stop('coupling_resample: failed object check.')
    select case(trim(this%type))
    case('hc');call resample(this%primitive)
    case('bilinear');call resample(this%bilinear)
    case default
       call stop('coupling_update: unknown coupling type.')
    end select
  end subroutine coupling_resample
  subroutine coupling_display(this,msg)
    type(coupling),intent(in)::this    
    character*(*),intent(in),optional::msg    
    logical::withmsg
    call Note('Begin coupling_display.')
    if(check(this))then
       call warn('coupling_display, failed object check','displaying nothing.')
       return
    end if
    withmsg=.false.
    if(present(msg))withmsg=.true.
    select case(trim(this%type))
    case('hc')
       if(withmsg)call display(this%hc,msg)
       if(.not.withmsg)call display(this%hc)
    case('bilinear')
       if(withmsg)call display(this%bilinear,msg)
       if(.not.withmsg)call display(this%bilinear)
    case default
       call stop('coupling_display: unknown coupling type.')
    end select
  end subroutine coupling_display
  subroutine coupling_save(this,file)
    type(coupling),intent(in)::this
    character*(*),intent(in)::file    
    integer(long)::unit
    call Note('Begin coupling_save.')
    call Note('input file= '//file)
    if(check(this))call stop('coupling_save: failed object check.')
    unit=newunit()
    open(unit,file=file)
    write(unit,*)'coupling'
    write(unit,*)trim(this%type)
    write(unit,*)quote(file//'.cp')
    select case(trim(this%type))
    case('hc');call save(this%primitive,file//'.cp')
    case('bilinear');call save(this%bilinear,file//'.cp')
    case default
       call stop('coupling_update: unknown coupling type.')
    end select
  end subroutine coupling_save
  integer(short)function coupling_check(this)
    type(coupling),intent(in)::this
    call Note('Checking coupling object.')
    coupling_check=0
    if(.not.this%initialized)then
       call Warn('coupling_check: coupling not initialized.')
       coupling_check=1
       return
    end if
    if(.not.associated(this%hc))then
       call Warn('coupling_check: coupling primitive not associated.')
       coupling_check=1
       return
    end if
    if(check(this%hc))then
       call Warn('coupling_check: coupling primitive failed check.')
       coupling_check=1
       return
    end if
    select case(trim(this%type))
    case('hc')
    case('bilinear')
    case default
       call stop('coupling_update: unknown coupling type: '//trim(this%type))
       coupling_check=1
       return
    end select
  end function coupling_check
end module coupling_class
