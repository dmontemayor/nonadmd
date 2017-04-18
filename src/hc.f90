!============================================================================
!>\brief
!! Coupling Term Primitive
!!\details
!!This class defines the coupling term kernal inherited by all derived
!!coupling terms.
!!\authors
!! Daniel Montemayor
!<============================================================================
  module hc_class
    use type_kinds
    use atomicunits
    use ErrorLog
    use filemanager
    use textgraphs
    use string
    use outputdisplay

    use hs_class
    use hb_class


    implicit none
    private

    public::hc
    public::new,kill,update,resample,display,save,check
    
    type hc
       logical::initialized=.false.
       type(hs),pointer::hs
       type(hb),pointer::hb
       complex(double),dimension(:,:),pointer::V
       complex(double),dimension(:,:,:),pointer::dV
    end type hc
    
    interface new
       module procedure hc_init
    end interface

    interface kill
       module procedure hc_kill
    end interface
    
    interface update
       module procedure hc_update
    end interface
    
    interface resample
       module procedure hc_resample
    end interface
    
    interface display
       module procedure hc_display
    end interface

    interface save
       module procedure hc_save
    end interface

    interface check
       module procedure hc_check
    end interface
       
  contains
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    subroutine hc_init(this,qs,cs,file)
      use quantum_class
      use classical_class

      type(hc),intent(inout)::this
      type(quantum),intent(inout),target::qs
      type(classical),intent(inout),target::cs
      character*(*),intent(in),optional::file
      
      integer(long)::istate,jstate
      integer(short)::ierr
      logical::usedefault

      call Note('Begin hc_init.')
            
      if(check(qs).NE.0)call Stop('hc_init: quantum object failed check')
      if(check(cs).NE.0)call Stop('hc_init: classical object failed check')

      if(associated(this%hs))nullify(this%hs)
      this%hs => qs%hs
      if(associated(this%hb))nullify(this%hb)
      this%hb => cs%hb

      if(associated(this%V))nullify(this%V)
      allocate(this%V(this%hs%nstate,this%hs%nstate),stat=ierr)
      if(ierr.NE.0)then
         write(*,*)'hc_init: potential V allocation error'
         stop
      end if
      this%V=0_double

      if(associated(this%dV))nullify(this%dV)
      allocate(this%dV(this%hb%ndof,this%hs%nstate,this%hs%nstate),stat=ierr)
      if(ierr.NE.0)then
         write(*,*)'hc_init: vector dV allocation error'
         stop
      end if
      this%dV=0_double

      this%initialized=.true.
      if(check(this).NE.0)call stop('hc_init: coupling primitive failed check.')
      
    end subroutine hc_init
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    subroutine hc_kill(this)
      type(hc),intent(inout)::this

      call Note('Begin hc_kill.')

      if(associated(this%hs))nullify(this%hs)
      if(associated(this%hb))nullify(this%hb)
      if(associated(this%V))nullify(this%V)
      if(associated(this%dV))nullify(this%dV)
      this%initialized=.false.
    end subroutine hc_kill
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    subroutine hc_update(this)
      type(hc),intent(inout)::this

      call Note('Begin hc_update.')

    end subroutine hc_update
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    subroutine hc_resample(this)
      type(hc),intent(inout)::this

      call Note('Begin hc_resample.')
      
    end subroutine hc_resample
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    subroutine hc_display(this,msg)
      type(hc),intent(in)::this
      character*(*),intent(in),optional::msg
      integer(long)::idof,istate,jstate
      
      call Note('Begin hc_display.')
      if(check(this).NE.0)then
         call warn('hc_display: failed check','displaying nothing.')
         return
      end if
      write(Dunit,*)'____________________________________________________'
      write(Dunit,*)'-----------------       hc       -------------------'
      if(present(msg))write(Dunit,*)msg
      write(Dunit,*)'----------------------------------------------------'
      call display(this%V/invcm,msg='Energy in 1/cm')
      write(Dunit,*)'Diagonal Values in 1/cm =',(this%V(istate,istate)/invcm,istate=1,this%hs%nstate)
      write(Dunit,*)'===================================================='
     
    end subroutine hc_display
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    subroutine hc_save(this,file)
      type(hc),intent(in)::this
      character*(*)::file
      integer(long)::unit

      call Note('Begin hc_save.')
      call Note('input file='//file)
      
      if(check(this).NE.0)then
         call warn('hc_save: failed check.','not saving object.')
      else
         unit=newunit()
         open(unit,file=file)
         
         write(unit,*)'hc'
         close(unit)
      end if

    end subroutine hc_save
    !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    integer(short) function hc_check(this)
      type(hc),intent(in)::this
      
      integer(long)::istate,jstate,idof
      integer(short)::ierr

      call Note('Checking hc.')

      hc_check=0

      !logical::initialized=.false.
      if(.not.this%initialized)then
         call Warn('hc_check: coupling primitive not initialized.')
         hc_check=1
         return
      end if

      !type(hs),pointer::hs
      if(check(this%hs).NE.0)then
         call Warn('hc_check: quantum primitive failed check')
         hc_check=1
         return
      end if

      !type(hb),pointer::hb
      if(check(this%hb).NE.0)then
         call Warn('hc_check: classical primitive failed check')
         hc_check=1
         return
      end if

      !real(double),dimension(:,:),pointer::V
      if(.not.associated(this%V))then
         call Warn('hc_check: potential V memory not allocated')
         hc_check=1
         return
      end if
      if(size(this%V).NE.this%hs%nstate**2)then
         call Warn('hc_check: number of potential V elements not equal number of quantum states squared')
         hc_check=1
         return
      end if
      ierr=0
      do istate=1,this%hs%nstate
         do jstate=1,this%hs%nstate
            if(this%V(istate,jstate).NE.this%V(istate,jstate))ierr=1 
            if(abs(this%V(istate,jstate)).GE.huge(real(this%V(istate,jstate))))ierr=1 
         end do
      end do
      if(ierr.EQ.1)then
         call Warn('hc_check: potential V has bad values.')
         hc_check=1
         return
      end if

      !real(double),dimension(:,:,:),pointer::dV
      if(.not.associated(this%dV))then
         call Warn('hc_check: vector dV memory not allocated')
         hc_check=1
      end if
      if(size(this%dV).NE.this%hs%nstate**2*this%hb%ndof)then
         call Warn('hc_check: number of potential V elements not equal number'&
              //' of quantum states squared times classical degrees of freedom')
         hc_check=1
         return
      end if
      ierr=0
      do idof=1,this%hb%ndof
         do istate=1,this%hs%nstate
            do jstate=1,this%hs%nstate
               if(this%dV(idof,istate,jstate).NE.this%dV(idof,istate,jstate))ierr=1 
               if(abs(this%dV(idof,istate,jstate)).GE.huge(real(this%dV(idof,istate,jstate))))ierr=1 
            end do
         end do
      end do
      if(ierr.EQ.1)then
         call Warn('hc_check: potential dV has bad values.')
         hc_check=1
         return
      end if
      
    end function hc_check

  end module hc_class
  
