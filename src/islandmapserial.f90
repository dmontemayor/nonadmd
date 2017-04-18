module islandmapserial_class
  use math
  use atomicunits
  use rand
  use quantum_class
  use classical_class
  use coupling_class
  implicit none
  Private

  Public::islandmapserial
  Public::new,save,run

  type islandmapserial
     logical::initialized=.false.
     real*8::runtime
     real*8::shorttime
     real*8::dtn
     real*8::dte
     integer::ntraj
     integer::nattempt
     integer::nbstep
     integer::nlstep
     integer::nstep
     integer::nbstepmax
     integer::nstepmax
     type(quantum),pointer::qs
     type(classical),pointer::cs
     type(coupling),pointer::cp
     real*8,dimension(:),pointer::q,p,qt,pt
     complex*16,dimension(:,:,:),pointer::rhor
  end type islandmapserial

  interface new
     module procedure new_islandmapserial
  end interface

  interface save
     module procedure save_islandmapserial
  end interface

  interface run
     module procedure run_islandmapserial
  end interface

contains
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine new_islandmapserial(this,qs,cs,cp,file)
    implicit none
    type(islandmapserial),intent(inout)::this
    type(quantum),intent(in),target::qs
    type(classical),intent(in),target::cs
    type(coupling),intent(in),target::cp
    character*(*),intent(in),optional::file

    integer::ierr,unit,nstate
    logical::usedefault,usedunit
    character(len=50)::filetype      

    if(present(file))then
       inquire(file=file,exist=usedefault)
       if(.not.usedefault)then
          write(*,*)'new_islandmapserial error cannot find inputfile '//file
          stop
       end if
    end if

    unit=1000
    usedunit=.true.
    do while (usedunit)
       unit=unit+1
       inquire(unit,opened=usedunit)
    end do

    if(present(file))then
       open(unit,file=file)
       read(unit,'(A50)')filetype
       filetype=adjustl(filetype)
       if(trim(filetype).NE.'islandmapserial')then
          write(*,*)'islandmapserial_init error: not a valid input file.'
          stop
       end if
       read(unit,*)this%ntraj
       read(unit,*)this%runtime
       read(unit,*)this%shorttime
       read(unit,*)this%dtn
       read(unit,*)this%dte
    else
       this%ntraj=0
       do while(this%ntraj.LT.1)
          write(*,*)'Enter total number of trajectories.'
          read(*,*)this%ntraj
          if(this%ntraj.LT.1)then
             write(*,*)'Number of trajectories must be positive. try again!'
          end if
       end do
       this%runtime=0.0
       do while(this%runtime.LE.0.0)
          write(*,*)'Enter total runtime.'
          read(*,*)this%runtime
          if(this%runtime.LE.0.0)then
             write(*,*)'runtime must be positive. try again!'
          end if
       end do
       this%shorttime=this%runtime*2.0_8
       do while(this%shorttime.GT.this%runtime)
          write(*,*)'Enter total short time approx.'
          read(*,*)this%shorttime
          if(this%shorttime.GT.this%runtime)then
             write(*,*)'short time must be < or = to runtime. try again!'
             write(*,*)'respective vaules:',this%shorttime,this%runtime
          end if
       end do
       this%dtn=this%shorttime
       do while(this%dtn.GE.this%shorttime)
          write(*,*)'Enter nuclear timestep.'
          read(*,*)this%dtn
          if(this%dtn.GE.this%shorttime)then
             write(*,*)'nuclear time step must be less than shorttime. try again'
          end if
       end do
       this%dte=this%dtn
       do while(this%dte.GE.this%dtn)
          write(*,*)'Enter electronic timestep.'
          read(*,*)this%dte
          if(this%dte.GE.this%dtn)then
             write(*,*)'electronic time step must be less than nuclear time step. try again.'
          end if
       end do
    end if
    this%nattempt=int(this%runtime/this%shorttime)
    this%nbstep=int(this%shorttime/this%dtn)
    this%nstep=this%nattempt*this%nbstep
    this%nlstep=int(this%dtn/this%dte)

    this%nbstepmax=100
    if(this%nbstepmax.GT.this%nbstep)this%nbstepmax=this%nbstep
    this%nstepmax=this%nbstepmax*this%nattempt

    if(.NOT.present(file))then
       write(*,*)'Warning nstepmax is hard coded: nstepmax=',this&
            &%nstepmax
       write(*,*)'nstep=',this%nstep
    end if

    if(.not.qs%initialized)then
       write(*,*)'islandmapserial_init error: quantum object not initialized.'
       stop
    end if
    if(associated(this%qs))nullify(this%qs)
    this%qs=>qs
    if(associated(this%rhor))nullify(this%rhor)
    allocate(this%rhor(this%qs%hs%nstate,this%qs%hs%nstate,0:this%nstepmax))
    if(associated(this%q))nullify(this%q)
    allocate(this%q(this%qs%hs%nstate))
    if(associated(this%p))nullify(this%p)
    allocate(this%p(this%qs%hs%nstate))
    if(associated(this%qt))nullify(this%qt)
    allocate(this%qt(this%qs%hs%nstate))
    if(associated(this%pt))nullify(this%pt)
    allocate(this%pt(this%qs%hs%nstate))

    if(.not.cs%initialized)then
       write(*,*)'islandmapserial_init error: classical object not initialized.'
       stop
    end if
    if(associated(this%cs))nullify(this%cs)
    this%cs=>cs

    if(.not.cp%initialized)then
       write(*,*)'islandmapserial_init error: coupling object not initialized.'
       stop
    end if
    if(associated(this%cp))nullify(this%cp)
    this%cp=>cp

    if(present(file))close(unit)

    this%initialized=.true.

  end subroutine new_islandmapserial
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine save_islandmapserial(this,file)
    implicit none
    type(islandmapserial),intent(in)::this
    character*(*),intent(in)::file

    integer::unit
    logical::usedunit    

    !write(*,*)'Warning islandmapserial_save: is corrupting the islandmapserial data type'

    if(.not.this%initialized)then
       write(*,*)'islandmapserial_save error: dynamics not initialized.'
       stop
    end if

    unit=1000
    usedunit=.true.
    do while (usedunit)
       unit=unit+1
       inquire(unit,opened=usedunit)
    end do

    open(unit,file=file)
    write(unit,'(A50)')'islandmapserial'
    write(unit,*)this%ntraj
    write(unit,*)this%runtime
    write(unit,*)this%shorttime
    write(unit,*)this%dtn
    write(unit,*)this%dte
    close(unit)  

  end subroutine save_islandmapserial
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine run_islandmapserial(this,rhortfile,tol,logfile)
    implicit none
    type(islandmapserial),intent(inout)::this
    character*(*),intent(in),optional::rhortfile,logfile
    real*8,intent(in),optional::tol

    integer::rhortunit
    logical::usedunit

    !real*8,parameter::tol=1e-5
    real*8,parameter:: x0i=1.0_8,p0i=1.0_8
    real*8,parameter::x0it=1.0_8,p0it=-1.0_8

    real*8,allocatable::qsjudge(:,:,:),psjudge(:,:,:)
    real*8,allocatable::qsin(:),psin(:),qsnew(:),psnew(:),fnon(:)
    complex*16,allocatable::bigA(:,:),bigAjudge(:,:,:,:,:)
    real*8,allocatable::PE(:,:,:),KE(:,:,:),CE(:,:,:)

    real*8:: dtn,dtn2,tol0,Etol
    real*8::wmc,wnorm,cume,rsq,norm
    real*8::w,wt,wn,wtn,sumfnon
    complex::aweight

    integer::nstate,ndof,nstep,ntraj

    integer::iatt,itraj,istep
    integer::istate_read,istatet_read
    integer::istate,istatet
    integer::fstate,fstatet
    integer::fs,fst
    integer::i,j

    integer::ntrajdump,idump
    real*8::dtdump

    real*8::gamma
    integer::maxstep,nadj


    tol0=1E-3
    if(present(tol))tol0=tol
    dtn=this%dtn
    dtn2=0.5*dtn*dtn
    nstate=this%qs%hs%nstate
    ndof=this%cs%hb%ndof

    ntrajdump=1!this%cs%hb%npt/100
    if(ntrajdump.EQ.0)ntrajdump=1
    dtdump=this%shorttime/this%nbstepmax
    ntraj=this%ntraj

    if(allocated(PE))deallocate(PE)
    allocate(PE(nstate,nstate,this%nbstepmax))
    if(allocated(KE))deallocate(KE)
    allocate(KE(nstate,nstate,this%nbstepmax))
    if(allocated(CE))deallocate(CE)
    allocate(CE(nstate,nstate,this%nbstepmax))


    if(allocated(bigA))deallocate(bigA)
    allocate(bigA(nstate,nstate))
    if(allocated(bigAjudge))deallocate(bigAjudge)
    allocate(bigAjudge(nstate,nstate,nstate,nstate,this%nstepmax))
    if(allocated(qsjudge))deallocate(qsjudge)
    allocate(qsjudge(ndof,nstate,nstate))
    if(allocated(psjudge))deallocate(psjudge)
    allocate(psjudge(ndof,nstate,nstate))
    if(allocated(qsin))deallocate(qsin)
    allocate(qsin(ndof))
    if(allocated(psin))deallocate(psin)
    allocate(psin(ndof))
    if(allocated(qsnew))deallocate(qsnew)
    allocate(qsnew(ndof))
    if(allocated(psnew))deallocate(psnew)
    allocate(psnew(ndof))
    if(allocated(fnon))deallocate(fnon)
    allocate(fnon(ndof))

    if(present(logfile))&
         open(1234,file=trim(logfile))
    !open(1234,file='debug/dyn.con')

    istate_read=1
    istatet_read=1

    if(present(logfile))&
         write(*,*)'WARNING istate_read and istatet_read are hard coded:', istate_read,istatet_read

    !initialize quantum subsystem
    Etol=5.0*invcm
    if(present(logfile))then
       write(*,*)'Warning islandmapserial will not consider weakly coupled final pairs of states'
       write(*,*)'with coupling strength less than ',Etol/invcm,'wavenumbers.'
    end if

    do fs=1,nstate
       do fst=1,nstate
          if(.not.(abs(this%qs%hs%diabat(fs,fst)).GT.Etol&
               .or.(fs.eq.fst))) this%qs%hs%diabat(fs,fst)=0.0_8
       end do
    end do



    !call display(this%qs%hs)
    !stop

    !initialize density matrix
    this%rhor=cmplx(0.,0.)
    !norm=0.0_8
    !do i=1,nstate
    !   do j=1,nstate
    !      this%rhor(i,j,0)=this%qs%hs%psi0(i)*this%qs%hs%psi0(j)
    !   end do
    !   norm=norm+this%rhor(i,i,0)
    !end do
    !this%rhor=this%rhor*ntraj

    do itraj=1,ntraj

       call resample(this%cs)

!!$       !initial Monte Carlo branching
!!$       wnorm=0.
!!$       do i=1,nstate
!!$          do j=1,nstate
!!$             wnorm=wnorm+sqrt(this%rhor(i,j,0)&
!!$                  *conjg(this%rhor(i,j,0)))
!!$          end do
!!$       end do
!!$       
!!$       w=ran0()*wnorm
!!$       cume=0.
!!$       
!!$       do i=1,nstate
!!$          do j=1,nstate
!!$             rsq=sqrt(this%rhor(i,j,0)&
!!$                  *conjg(this%rhor(i,j,0)))
!!$             cume=cume+rsq
!!$             if(w.le.cume) then
!!$                istate=i
!!$                istatet=j
!!$                goto 20
!!$             end if
!!$          end do
!!$       end do
!!$20     continue


       istatet = istatet_read
       istate = istate_read

       this%rhor(istate,istatet,0)=this%rhor(istate,istatet,0)+1.0_8


       !Initialize trajectory MC weight
       wmc=1._8

       !Initialize trajectory Phase weight
       aweight=(1._8,0.)

       !Compute force 
       call update(this%cs)
       call update(this%cp)


       !initialize mapping variables before time loop
       qsin=this%cs%hb%Q
       psin=this%cs%hb%P


       !Hop attempt loop
       do iatt=1,this%nattempt

          PE=0.0_8
          KE=0.0_8
          CE=0.0_8

          !Loop over possible final states
          do fs=1,nstate
             do fst=1,nstate
                if(.not.(abs(this%qs%hs%diabat(fs,fst)).GT.Etol&
                     .or.(fs.eq.fst)))then
                   !write(*,*)'skipping final state pair',fs,fst
                else
                   
                   !write(*,*)itraj,iatt,fs,fst,abs(this%qs%hs%diabat(fs,fst))
                   
                   this%q=0.
                   this%p=0.
                   this%qt=0.
                   this%pt=0.

                   if(iatt.eq.1) then
                      this%q(istate_read)=x0i
                      this%p(istate_read)=p0i
                      this%qt(istatet_read)=x0it
                      this%pt(istatet_read)=p0it

                      this%cs%hb%Q=qsin
                      this%cs%hb%P=psin
                   else
                      this%q(istate)=x0i
                      this%p(istate)=p0i
                      this%qt(istatet)=x0it
                      this%pt(istatet)=p0it

                      this%cs%hb%Q=qsnew
                      this%cs%hb%P=psnew
                   end if
                   !Compute New force 
                   call update(this%cs)
                   call update(this%cp)

                   !Calculate nonadiabatic Forces
                   wn=(this%q(fs)**2+this%p(fs)**2)
                   wtn=(this%qt(fst)**2+this%pt(fst)**2)
                   fnon=0.0_8
                   do i=1,nstate
                      w=0.0_8
                      wt=0.0_8
                      if(wn.GT.tol0)w=(this%q(fs)*this%q(i)&
                           +this%p(fs)*this%p(i))/wn
                      if(i.EQ.fs)w=1.0
                      !if(abs(w).GT.1E15)w=0.0
                      if(wtn.GT.tol0)wt=(this%qt(fst)*this%qt(i)&
                           +this%pt(fst)*this%pt(i))/wtn
                      if(i.EQ.fst)wt=1.0
                      !if(abs(wt).GT.1E15)wt=0.0
                      fnon=fnon-this%cp%hc%dV(:,i,fs)*w&
                           -this%cp%hc%dV(:,i,fst)*wt
                   end do
                   fnon=0.5_8*fnon+this%cs%hb%F
                   sumfnon=sum(fnon)
                   if(sumfnon.NE.sumfnon)then
                      write(*,*)'Warning NAN nondaiabatic force.'
                      !write(*,*)'Program will recover by assuming 0 force.'
                      write(*,*)'Program will stop.'
                      stop
                      fnon=0.0
                   end if

                   istep=0
                   if(fs.EQ.1.and.fst.EQ.1.and.itraj.EQ.1)then
                      call update(this%cs)
                      call update(this%cp)

                      if(present(logfile))then
                         write(1234,"(1000(ES18.10E2,1X))")(iatt-1)*this%shorttime+istep*dtn,this&
                              &%cs%hb%V,this%cs%hb%T&
                              ,this%cp%hc%V(fs&
                              &,fst),sum(this%q**2+this%p**2)&
                              &,sum(this%qt**2+this%pt**2)
                      end if
                   end if

                   !Main time stepping loop
                   idump=1
                   do istep=1,this%nbstep

                      !Advance nuclear positions with current forces
                      !first half Ver(let) step 
                      this%cs%hb%Q=this%cs%hb%Q&
                           +this%cs%hb%rmass*this%cs%hb%P*dtn&
                           +this%cs%hb%rmass*fnon*dtn2
                      this%cs%hb%P=this%cs%hb%P&
                           +0.5*fnon*dtn

                      !Compute New force 
                      call update(this%cs)
                      call update(this%cp)

                      !Advance mapping variables 
                      call mapverlet(this)

                      !Construct bigA
                      do i=1,nstate
                         do j=1,nstate
                            bigA(i,j)=0.25*(this%q(i)+eye*this%p(i))&
                                 *(this%qt(j)-eye*this%pt(j))*(x0i-eye*p0i)&
                                 *(x0it+eye*p0it)
                         end do
                      end do

                      !Calculate nonadiabatic Forces
                      wn=(this%q(fs)**2+this%p(fs)**2)
                      wtn=(this%qt(fst)**2+this%pt(fst)**2)
                      fnon=0.0_8
                      do i=1,nstate
                         w=0.0_8
                         wt=0.0_8
                         if(wn.GT.tol0)w=(this%q(fs)*this%q(i)&
                              +this%p(fs)*this%p(i))/wn
                         if(i.EQ.fs)w=1.0
                         !if(abs(w).GT.1E15)w=0.0
                         if(wtn.GT.tol0)wt=(this%qt(fst)*this%qt(i)&
                              +this%pt(fst)*this%pt(i))/wtn
                         if(i.EQ.fst)wt=1.0
                         !if(abs(wt).GT.1E15)wt=0.0
                         fnon=fnon-this%cp%hc%dV(:,i,fs)*w&
                              -this%cp%hc%dV(:,i,fst)*wt
                      end do
                      fnon=0.5_8*fnon+this%cs%hb%F
                      sumfnon=sum(fnon)
                      if(sumfnon.NE.sumfnon)then
                         write(*,*)'Warning NAN nondaiabatic force.'
                         !write(*,*)'Program will recover by assuming 0 force.'
                         write(*,*)'Program will stop.'
                         stop
                         fnon=0.0
                      end if

                      if(fs.EQ.1.and.fst.EQ.1.and.itraj.EQ.1)then
                         call update(this%cs)
                         call update(this%cp)
                         if(present(logfile))then
                            write(1234,"(1000(ES18.10E2,1X))")(iatt-1)*this%shorttime+istep*dtn,this&
                                 &%cs%hb%V,this%cs%hb%T&
                                 ,this%cp%hc%V(fs&
                                 &,fst),sum(this%q**2+this%p**2)&
                                 &,sum(this%qt**2+this%pt**2)
                         end if
                      end if
                      
                      !second half of (ver)Let
                      this%cs%hb%P=this%cs%hb%P+0.5*fnon*dtn


                      !record all bigA for reconstructing dynamics
                      ! during nbstep
                      !loop after states are decided
                      if(istep*dtn.GE.idump*dtdump)then

                         !write(*,*)itraj,iatt,fs,fst,istep

                         bigAjudge(:,:,fs,fst,idump)=bigA
                         PE(fs,fst,idump)=this%cs%hb%V
                         KE(fs,fst,idump)=this%cs%hb%T
                         CE(fs,fst,idump)=this%cp%hc%V(fs,fst)
                         !write(*,*)fs,fst,idump,real(bigA),aimag(bigA)
                         !write(bigAunit)fs,fst,idump,bigA
                         idump=idump+1
                      end if

                   end do!End of Nuclear Stepping loop


                   !save all current trajectory information
                   qsjudge(:,fs,fst)=this%cs%hb%Q
                   psjudge(:,fs,fst)=this%cs%hb%P

                end if!end Etol condition

             end do!end fst loop
          end do!end fs loop

          !Monte Carlo branching
          wnorm=0.
          do i=1,nstate
             do j=1,nstate
                wnorm=wnorm+sqrt(bigAjudge(i,j,i,j,this%nbstepmax)&
                     *conjg(bigAjudge(i,j,i,j,this%nbstepmax)))
             end do
          end do

          w=ran0()*wnorm
          cume=0.

          do i=1,nstate
             do j=1,nstate
                rsq=sqrt(bigAjudge(i,j,i,j,this%nbstepmax)&
                     *conjg(bigAjudge(i,j,i,j,this%nbstepmax)))
                cume=cume+rsq
                if(w.le.cume) then
                   fstate=i
                   fstatet=j
                   goto 10
                end if
             end do
          end do
10        continue

          qsnew=qsjudge(:,fstate,fstatet)
          psnew=psjudge(:,fstate,fstatet)
          bigA=bigAjudge(:,:,fstate,fstatet,this%nbstepmax)

          !reconstruct redmat up to current hopping time for recording
          do istep=1,this%nbstepmax
             !composite trajectory weight
             !accumulate average reduced density matrix
             if(istep.eq.this%nbstepmax) then
                wmc=wmc*wnorm
                !new calculation for the phase:
                aweight = aweight*bigA(fstate,fstatet)&
                     /sqrt(bigA(fstate,fstatet)*conjg(bigA(fstate,fstatet)))
                !accumulate average reduced density matrix
                this%rhor(fstate,fstatet,(iatt-1)*this%nbstepmax+istep)&
                     =this%rhor(fstate,fstatet,(iatt-1)*this%nbstepmax+istep)&
                     +wmc*aweight
             else
                do fs=1,nstate
                   do fst=1,nstate
                      this%rhor(fs,fst,(iatt-1)*this%nbstepmax+istep)&
                           =this%rhor(fs,fst,(iatt-1)*this%nbstepmax+istep)&
                           +wmc*aweight*bigAjudge(fs,fst,fstate,fstatet,istep)
                   end do
                end do
             end if

          end do!end reconstruction loop
          istate=fstate
          istatet=fstatet


!          rhortunit=1000
!          usedunit=.true.
!          do while (usedunit)
!             rhortunit=rhortunit+1
!             inquire(rhortunit,opened=usedunit)
!          end do
!          open(rhortunit,file=rhortfile//'.PE',POSITION='APPEND')
!          do istep=1,this%nbstepmax
!             write(rhortunit,*)((iatt-1)*this%nbstepmax+istep)*dtdump,itraj&
!                  ,PE(istate,istatet,istep),KE(istate,istatet,istep),CE(istate,istatet,istep),istate,istatet
!          end do
!          close(rhortunit)

       end do!end nattempt loop

       if(present(rhortfile))then
          if (mod(itraj,ntrajdump).EQ.0)then
             !write(*,*)'run_islandmapserial: saving intermediate reduced density in '//rhortfile
             rhortunit=1000
             usedunit=.true.
             do while (usedunit)
                rhortunit=rhortunit+1
                inquire(rhortunit,opened=usedunit)
             end do
             open(rhortunit,file=rhortfile)
             do istep=0,this%nstepmax
                write(rhortunit,"(1000(ES18.10E2,1X))")istep*dtdump,((real(this%rhor(i,j&
                     &,istep)),aimag(this%rhor(i,j&
                     &,istep)),j=1,nstate),i=1,nstate)
             end do
             close(rhortunit)
          end if
       end if

    end do!end ntraj loop

    close(1234)

    if(allocated(PE))deallocate(PE)
    if(allocated(KE))deallocate(KE)
    if(allocated(CE))deallocate(CE)


    if(allocated(bigA))deallocate(bigA)
    if(allocated(bigAjudge))deallocate(bigAjudge)
    if(allocated(qsjudge))deallocate(qsjudge)
    if(allocated(psjudge))deallocate(psjudge)
    if(allocated(qsin))deallocate(qsin)
    if(allocated(psin))deallocate(psin)
    if(allocated(qsnew))deallocate(qsnew)
    if(allocated(psnew))deallocate(psnew)
    if(allocated(fnon))deallocate(fnon)

    if(present(rhortfile))then

       rhortunit=1000
       usedunit=.true.
       do while (usedunit)
          rhortunit=rhortunit+1
          inquire(rhortunit,opened=usedunit)
       end do
       open(rhortunit,file=rhortfile)
       do istep=0,this%nstepmax
          write(rhortunit,"(1000(ES18.10E2,1X))")istep*dtdump,((real(this%rhor(i,j,istep))&
               &,aimag(this%rhor(i,j,istep)),j=1,nstate),i=1,nstate)
       end do
       close(rhortunit)

    end if

  end subroutine run_islandmapserial
  !------------------------------------------------------------------------

  subroutine mapverlet(this)
    implicit none
    type(islandmapserial),intent(inout)::this

    real*8,allocatable,dimension(:)::xrhs,prhs,xtrhs,ptrhs,force,forcet
    real*8,allocatable,dimension(:,:)::hel
    real*8 xsumdum,psumdum,dt
    integer i,n,m,nstate

    dt=this%dte
    nstate=this%qs%hs%nstate

    if(allocated(hel))deallocate(hel)
    if(allocated(xrhs))deallocate(xrhs)
    if(allocated(prhs))deallocate(prhs)
    if(allocated(xtrhs))deallocate(xtrhs)
    if(allocated(ptrhs))deallocate(ptrhs)
    if(allocated(force))deallocate(force)
    if(allocated(forcet))deallocate(forcet)
    allocate(hel(nstate,nstate))
    allocate(xrhs(nstate))
    allocate(prhs(nstate))
    allocate(xtrhs(nstate))
    allocate(ptrhs(nstate))
    allocate(force(nstate))
    allocate(forcet(nstate))

    hel=0.0_8
    do n=1,nstate
       do m=1,nstate
          hel(n,m)=hel(n,m)+this%cp%hc%V(n,m)+this%qs%hs%diabat(n,m)
       end do
       hel(n,n)=hel(n,n)+this%cs%hb%V
    end do

    do i=1,this%nlstep

       !     Generate current derivatives

       !     Forward

       do n=1,nstate
          xrhs(n)=hel(n,n)*this%p(n)
          prhs(n)=-hel(n,n)*this%q(n)
       end do

       do n=1,nstate
          xsumdum=0.
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                xsumdum=xsumdum+hel(n,m)*this%p(m)
                psumdum=psumdum+hel(n,m)*this%q(m)
             end if
          end do
          xrhs(n)=xrhs(n)+xsumdum
          prhs(n)=prhs(n)-psumdum
       end do

       !     Backward

       do n=1,nstate
          xtrhs(n)=hel(n,n)*this%pt(n)
          ptrhs(n)=-hel(n,n)*this%qt(n)
       end do

       do n=1,nstate
          xsumdum=0.
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                xsumdum=xsumdum+hel(n,m)*this%pt(m)
                psumdum=psumdum+hel(n,m)*this%qt(m)
             end if
          end do
          xtrhs(n)=xtrhs(n)+xsumdum
          ptrhs(n)=ptrhs(n)-psumdum
       end do

       !     Generate current second derivatives

       !     Forward

       do n=1,nstate
          force(n)=hel(n,n)*prhs(n)
       end do


       do n=1,nstate
          xsumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                xsumdum=xsumdum+hel(n,m)*prhs(m)
             end if
          end do
          force(n)=force(n)+xsumdum
       end do

       !     backward

       do n=1,nstate
          forcet(n)=hel(n,n)*ptrhs(n)
       end do


       do n=1,nstate
          xsumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                xsumdum=xsumdum+hel(n,m)*ptrhs(m)
             end if
          end do
          forcet(n)=forcet(n)+xsumdum
       end do

       !     Advance Ver step

       do n=1,nstate
          this%q(n)=this%q(n)+xrhs(n)*dt+0.5*force(n)*dt*dt
          this%p(n)=this%p(n)+0.5*prhs(n)*dt

          this%qt(n)=this%qt(n)+xtrhs(n)*dt+0.5*forcet(n)*dt*dt
          this%pt(n)=this%pt(n)+0.5*ptrhs(n)*dt
       end do

       !     Compute new first derivatives

       !     Forward

       do n=1,nstate
          prhs(n)=-hel(n,n)*this%q(n)
       end do

       do n=1,nstate
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                psumdum=psumdum+hel(n,m)*this%q(m)
             end if
          end do
          prhs(n)=prhs(n)-psumdum
       end do

       !     Backward

       do n=1,nstate
          ptrhs(n)=-hel(n,n)*this%qt(n)
       end do

       do n=1,nstate
          psumdum=0.
          do m=1,nstate
             if(m.ne.n) then
                psumdum=psumdum+hel(n,m)*this%qt(m)
             end if
          end do
          ptrhs(n)=ptrhs(n)-psumdum
       end do

       !    Advance let step

       do n=1,nstate
          this%p(n)=this%p(n)+0.5*prhs(n)*dt
          this%pt(n)=this%pt(n)+0.5*ptrhs(n)*dt
       end do

    end do

    if(allocated(hel))deallocate(hel)
    if(allocated(xrhs))deallocate(xrhs)
    if(allocated(prhs))deallocate(prhs)
    if(allocated(xtrhs))deallocate(xtrhs)
    if(allocated(ptrhs))deallocate(ptrhs)
    if(allocated(force))deallocate(force)
    if(allocated(forcet))deallocate(forcet)

    return
  end subroutine mapverlet


end module islandmapserial_class
