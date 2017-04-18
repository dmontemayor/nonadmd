!-----------------------------------------------------------------------------
!> \brief
!! Harmonic bath(s) of classical non-interacting particles.
!> \details
!! Assumes Hamiltonian of the form
!!        \f[  H=hs+H_b+hc \f]
!! with standard primative quantum and coupling subsystems hs and hc repsectively.
!! The classical subsystem is composed of \f$ N \f$ independent baths
!! each with \f$ K_n \f$ harmonic modes.
!! The classical subsystem hamiltonian takes the form
!!  \f[ H_b=\sum_n^N \sum_k^{K_n} \frac{P_{k,n}^2}{2m_{k,n}}+\frac{1}{2}m_{k,n} \omega_{k,n}^2 Q_{k,n}^2 .\f]
!! The spectral properties of these \f$ N \f$ independent baths are determined by
!!\f$ M \f$ input spectral denisties \f$ J_m(\omega) \f$.
!! More than one harmonic bath can be defined by the same input spectral density
!! , hence \f$ N>M \f$.
!! The oscillator frequencies \f$ \omega_{k,n} \f$ are sampled uniformly from the density of states 
!! of these input spectral densities. The density of states is written as follows 
!!\f[ \rho_m(\omega)=\frac{J_m(\omega)}{\omega} .\f]
!! The re-organization energy associated with input spectral density \f$ m \f$ is calculated
!!   \f[ \lambda_m=\frac{2}{\pi}\int_0^\infty\frac{J_m(\omega)}{\omega} d\omega \f]
!! and is computed at every NEW or RESAMPLE call. Because \f$H_b\f$ is quadratic in \f$Q\f$, the natural frequencies
!! \f$ \omega_{k,n} \f$ also define the normal modes of the classical subsystem primative hb. 
!!
!! \authors
!! Daniel Montemayor
!!
!! \date
!! Aug 2011
!!
!! \todo
!! Make stochasitcsample attribute editable on creation of oject.
!<-----------------------------------------------------------------------------
!Change Log:
!=== 12 Apr 2013 ===
! - fixed bug: map was not saving properly
! 
!===   v0.1 Oct 2012   ===
! + added stochasticsample logical variable (default false)
!
!===   v1.0 Aug 2011   ===
! + beta 1.0 compatable
! + defines multiple baths with equal number of dofs composed of one dimensional
!   non-interacting harmonic oscillators
! + maps multiple baths to input spectral densities
! + samples frequencies from input spectral densities and assigns them to the
!   normal modes of hb
!------------------------------------------------------------------------------!
module harmonicbath_class
  use type_kinds
  use ErrorLog
  use filemanager
  use string
  use math
  use atomicunits
  use hb_class
  use rand
  use textgraphs
  use outputdisplay
  implicit none
  private

  public::harmonicbath
  public::new,kill,update,resample,display,save,check

  type harmonicbath
     logical::initialized=.false.
     type(hb)::hb
     integer(long)::nbath,npt,nden
     real(double),dimension(:),pointer::wgrid
     real(double),dimension(:,:),pointer::Jw,intDOS
     real(double),dimension(:),pointer::lambda
     integer(long),dimension(:),pointer::map
     logical::stochasticsample
  end type harmonicbath

  interface new
     module procedure harmonicbath_init
  end interface

  interface kill
     module procedure harmonicbath_kill
  end interface

  interface display
     module procedure harmonicbath_display
  end interface

  interface save
     module procedure harmonicbath_save
  end interface

  interface update
     module procedure harmonicbath_update
  end interface

  interface resample
     module procedure harmonicbath_resample
  end interface

  interface check
     module procedure harmonicbath_check
  end interface

contains

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine harmonicbath_init(this,file)

    type(harmonicbath),intent(inout)::this
    character*(*),intent(in),optional::file
    character(len=path)::filename
    character(len=path)::hbfilename
    character(len=title)::filetype
    integer(short)::ierr
    integer(long)::unit
    logical::usedefault
    integer(long):: idof,jdof,i,ipt,ibath,bdof,idensity
    real(double)::dw,jbar,wbar,xscale,yscale
    real(double)::yr,slope,yint,y0

    real(double)::params(2),wo,wc,wmax,wj,cj,lo

    call Note('Begin harmonicbath_init.')
    filename='unknown'
    hbfilename='unknown'
    filetype='unknown'

    if(.not.present(file))then
       call Note('harmonicbath_init: input file is required.')
       write(*,*)'Please enter a specden file.'
       read(*,*)filename
       filename=adjustl(filename)
    else
       filename=file
    end if
    call Note('input file= '//trim(filename))
    if(check(trim(filename)).EQ.1)call stop('harmonicbath_init: input file '//trim(filename)//' not found')
    call note('Loading spectral density found in '//trim(filename))
    

    unit=newunit()
    open(unit,file=filename)
    read(unit,*)filetype
    filetype=adjustl(filetype)
    if(trim(filetype).EQ.'harmonicbath')then
       read(unit,*)hbfilename
       hbfilename=adjustl(hbfilename)
       call new(this%hb,trim(hbfilename))
       read(unit,*)this%npt
       read(unit,*)this%nden
    elseif(trim(filetype).EQ.'specden')then
       call new(this%hb)
       read(unit,*)this%nden
       if(this%nden.LE.0)call stop('harmonicbath_init: bad value in line 2 in '&
            //file,'checking that input file has correct format.')
       read(unit,*)xscale 
       read(unit,*)yscale
       this%npt=0
       do
          read(unit,*,END=222)
          this%npt=this%npt+1
       end do
222    rewind(unit)
       read(unit,*)!filetype
       read(unit,*)!nden
       read(unit,*)!xscale
       read(unit,*)!yscale
    else
       call Stop('harmonicbath_init: input file not valid type.')
    end if

    !allocate memory for specden and grid
    if(associated(this%wgrid))nullify(this%wgrid)
    allocate(this%wgrid(this%npt),stat=ierr)
    if(ierr.NE.0)&
         call stop('harmonicbath_init: failed to allocate memory for wgrid.')
    if(associated(this%Jw))nullify(this%Jw)
    allocate(this%Jw(this%npt,this%nden),stat=ierr)
    if(ierr.NE.0)&
         call stop('harmonicbath_init: failed to allocate memory for Jw.')
    this%wgrid=0._double
    this%Jw=0._double

    if(trim(filetype).EQ.'harmonicbath')then
       read(unit,*)(this%wgrid(ipt),ipt=1,this%npt)
       read(unit,*)((this%Jw(ipt,idensity),idensity=1,this%nden),ipt=1,this%npt)
    elseif(trim(filetype).EQ.'specden')then
       do ipt=1,this%npt
          read(unit,*) this%wgrid(ipt),(this%Jw(ipt,idensity),idensity=1,this%nden)
       end do
       !write(*,*) 'xscale,yscale',xscale,yscale
       this%wgrid=this%wgrid*xscale
       this%Jw=this%Jw*yscale
    end if
    
    if(trim(filetype).EQ.'harmonicbath')then
       read(unit,*)this%nbath
    elseif(trim(filetype).EQ.'specden')then
       this%nbath=0
       write(*,*)'Enter the number of isolated baths.'
       do while(this%nbath.LT.1)
          read(*,*)this%nbath
          if(this%nbath.LT.1)&
               write(*,*)'Number of baths must be a positve integer. Try again.'
       end do
    end if

    if(associated(this%map))nullify(this%map)
    allocate(this%map(this%nbath),stat=ierr)
    if(ierr.NE.0)&
         call stop('harmonicbath_init: failed to allocate memory for map.')

    if(trim(filetype).EQ.'harmonicbath')then
       read(unit,*)(this%map(ibath),ibath=1,this%nbath)
    elseif(trim(filetype).EQ.'specden')then
       this%map=0
       write(*,*)'Mapping '//trim(int2str(this%nden))&
            //' input spectral densities to '//trim(int2str(this%nbath))&
            //' independent baths.'
       if(this%nden.EQ.1)then
          this%map=1
       else
          do ibath=1,this%nbath
             i=0
             do while(i.LE.0.or.i.GT.this%nden)
                write(*,*)'Please associate bath '//trim(int2str(ibath))&
                     //' with a density [ options are from 1 to'//trim(int2str(this%nden))//'].'
                read(*,*)i
                if(i.LE.0.or.i.GT.this%nden)write(*,*)'Not a valid entry.'&
                     //' Please enter an integer between(including) 1 and '&
                     //trim(int2str(this%nden))//'. Try again.'
             end do
             this%map(ibath)=i
          end do
       end if
    end if

    this%stochasticsample=.false.
    if(trim(filetype).EQ.'harmonicbath')then
       read(unit,*)this%stochasticsample
    else
       write(*,*)'Use a stochasitc sample of the density of states. T/F'
       read(*,*)this%stochasticsample
    end if

    close(unit)

    if(associated(this%intDOS))nullify(this%intDOS)
    allocate(this%intDOS(this%npt,this%nden),stat=ierr)
    if(ierr.NE.0)&
         call stop('harmonicbath_init: failed to allocate memory for intDOS.')
    !integration by trapezoidal rule
    this%intDOS=0._double
    do idensity=1,this%nden
       y0=0._double
       do i=2,this%npt
          yr=0._double
          if(this%wgrid(i).LT.epsilon(this%wgrid))&
               call stop('harmonicbath_init: wgrid has tiny or negative values after first point.')
          yr=this%Jw(i,idensity)/this%wgrid(i)
          dw=this%wgrid(i)-this%wgrid(i-1)
          jbar=.5_double*(yr+y0)
          this%intDOS(i,idensity)=this%intDOS(i-1,idensity)+jbar*dw
          y0=yr
       end do
    end do

    !re-orginazation energy for each bath of bdof modes
    if(associated(this%lambda))nullify(this%lambda)
    allocate(this%lambda(this%nden),stat=ierr)
    if(ierr.NE.0)&
         call stop('harmonicbath_init: failed to allocate memory for lambda.')
    do idensity=1,this%nden    
       this%lambda(idensity)=maxval(this%intDOS(:,idensity))/pi
    end do

    !dofs per bath (total dofs / nbath)
    if(this%nbath.EQ.0)call Stop('harmonicbath_init:'&
         //' the number isolated baths is not a positive integer.')
    if(mod(this%hb%ndof,this%nbath).ne.0)call Stop('harmonicbath_init:'&
         //' the total number of classical dofs must divide equally'&
         //'among the isolated baths.')
    bdof=this%hb%ndof/this%nbath
    call Note('harmonicbath will divide the total number dofs equally'&
         //' among the isolated baths.')
    do ibath=0,this%nbath-1
       call Note('dofs '//trim(int2str(1+(ibath*bdof)))&
            //' through '//trim(int2str((1+ibath)*bdof))&
            //' will be associated with bath '//trim(int2str(1+ibath)))
    end do

    !assign static modes
    this%hb%mode=0.0_double
    do ibath=1,this%nbath
       idensity=this%map(ibath)
       do idof=1,bdof
          jdof=idof+(ibath-1)*bdof
          yr=pi*this%lambda(idensity)*idof/real(bdof,double)
          do i=2,this%npt
             dw=this%wgrid(i)-this%wgrid(i-1)
             if(yr.GT.this%intDOS(i-1,idensity).and.yr.LE.this%intDOS(i,idensity)+epsilon(yr))then
             !epsilon yr is used above to smooth out precision errors from reading from file.  
                if(dw.GT.epsilon(dw))then
                   !linearly extrapolate w to get inverse function
                   slope=(this%intDOS(i,idensity)-this%intDOS(i-1,idensity))/dw
                   yint=this%intDOS(i,idensity)-(slope*this%wgrid(i))
                   this%hb%mode(jdof)=(yr-yint)/slope
                else
                   this%hb%mode(jdof)=this%wgrid(i)
                end if
             end if
          end do
       end do
    end do

    !harmonic clasical subsystem, no work needed to compute normal modes 
    !so EigenVec is unit matrix
    this%hb%EigenVec=iden(this%hb%ndof)

    !W coodinates = Q coordinates beacuse EigenVec=iden
    this%hb%W=this%hb%Q

    this%initialized=.true.
    call resample(this)
    if(check(this).EQ.1)call stop('harmonicbath_init, failed exit object check')

  end subroutine harmonicbath_init
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine harmonicbath_kill(this)

    type(harmonicbath),intent(inout)::this

    call Note('Begin harmonicbath_kill.')
    call kill(this%hb)
    if(associated(this%wgrid))nullify(this%wgrid)
    if(associated(this%Jw))nullify(this%Jw)
    if(associated(this%intDOS))nullify(this%intDOS)
    if(associated(this%lambda))nullify(this%lambda)
    if(associated(this%map))nullify(this%map)
    this%initialized=.false.
  end subroutine harmonicbath_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine harmonicbath_display(this,msg)
    type(harmonicbath),intent(in)::this
    character*(*),intent(in),optional::msg
    integer(long)::idof,bdof,i,j,ibath,idensity
    real(double)::dw,w

    call Note('Begin harmonicbath_display.')
    if(check(this).NE.0)then
       call warn('harmonicbath_display: failed check','displaying nothing.')
       return
    end if
    bdof=this%hb%ndof/this%nbath

    write(Dunit,*)'____________________________________________________'
    write(Dunit,*)'-----------------  harmonicbath  -------------------'
    if(present(msg))write(Dunit,*)msg
    write(Dunit,*)'----------------------------------------------------'
    write(Dunit,*)'npt=',this%npt
    write(Dunit,*)'nden=',this%nden
    do idensity=1,this%nden
       write(Dunit,*)'Input spectral density '//trim(int2str(idensity))
       write(Dunit,*)'Re-organization energy(Eh)='//trim(float2str(this%lambda(idensity)))
       write(Dunit,*)'Re-organization energy(hbar/cm)='//trim(float2str(this%lambda(idensity)/invcm))
       call display(this%wgrid/invcm,this%Jw(:,idensity)/invcm,msg='Spectral Density in wavenumbers')
    end do
    write(Dunit,*)'nbath=',this%nbath
    do ibath=1,this%nbath
       i=1+(ibath-1)*bdof
       j=ibath*bdof
       write(Dunit,*)'Bath '//trim(int2str(ibath))//' definition.'
       write(Dunit,*)'Use input spectal density '//trim(int2str(this%map(ibath)))
       write(Dunit,*)'Classical subsystem DOFs ['//trim(int2str(i))//':'//trim(int2str(j))//']'
       call display(this%hb%Q(i:j)/angstrom,this%hb%P(i:j)*ps/angstrom,msg='Phase-space plot X in Angstrom, Y in Angstrom/ps')
       !call display(this%hb%mode(i:j)/invcm,msg='Normal modes (1/cm).')
       write(Dunit,*)'Normal modes in 1/cm :',this%hb%mode(i:j)/invcm
       write(Dunit,*)
    end do
    call display(this%hb,'harmonicbath primitive classical subsystem')
    write(Dunit,*)'===================================================='

  end subroutine harmonicbath_display
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine harmonicbath_save(this,file)
    type(harmonicbath),intent(in)::this
    character*(*),intent(in)::file

    integer(long)::unit
    integer(long)::ipt,idensity,ibath

    call Note('Begin harmonicbath_save.')
    call Note('input file= '//file)
    if(check(this).NE.0)then
       call warn('harmonicbath_save: failed check','not saving object.')
    else       
       unit=newunit()
       open(unit,file=file)
       write(unit,*)'harmonicbath'
       write(unit,*)quote(file//'.hb')
       call save(this%hb,file//'.hb')
       write(unit,*)this%npt
       write(unit,*)this%nden
       write(unit,*)(this%wgrid(ipt),ipt=1,this%npt)
       write(unit,*)((this%Jw(ipt,idensity),idensity=1,this%nden),ipt=1,this%npt)
       write(unit,*)this%nbath
       write(unit,*)(this%map(ibath),ibath=1,this%nbath)
       write(unit,*)this%stochasticsample
       close(unit)
    end if
  end subroutine harmonicbath_save
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine harmonicbath_update(this)
    type(harmonicbath),intent(inout)::this
    !integer(long)::ipt,idof

    call Note('Begin harmonicbath_update.')

    this%hb%F=-this%hb%mode**2*this%hb%Q/this%hb%rmass
    this%hb%V=0.5*sum((this%hb%mode*this%hb%Q)**2/this%hb%rmass)
    this%hb%T=0.5*sum(this%hb%rmass*this%hb%P**2)

    !normal modes and eigen vectors do not change (set in resample)
    this%hb%W=this%hb%Q

  end subroutine harmonicbath_update
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine harmonicbath_resample(this)

    type(harmonicbath),intent(inout)::this

    real(double)::sigp,sigq,Q0
    integer(long):: idof,jdof,ipt

    integer(long)::ibath,i,bdof,idensity
    real(double)::dw,wr,yr,slope,yint,beta,wbar,wbar0
    real(double)::wc,wo,ysum,ysum0,wj

    call Note('Begin harmonicbath_resample.')
    !call resample(this%hb)

    beta=1._double/(kb*this%hb%temperature)

    if(this%stochasticsample)then
       !assign random modes
       bdof=this%hb%ndof/this%nbath
       this%hb%mode=0.0_double
       do ibath=1,this%nbath
          idensity=this%map(ibath)
          do idof=1,bdof
             jdof=idof+(ibath-1)*bdof
             yr=pi*this%lambda(idensity)*ran0()
             do i=2,this%npt
                dw=this%wgrid(i)-this%wgrid(i-1)
                if(yr.GT.this%intDOS(i-1,idensity).and.yr.LE.this%intDOS(i,idensity)+epsilon(yr))then
                   !epsilon yr is used above to smooth out precision errors from reading from file.  
                   if(dw.GT.epsilon(dw))then
                      !linearly extrapolate w to get inverse function
                      slope=(this%intDOS(i,idensity)-this%intDOS(i-1,idensity))/dw
                      yint=this%intDOS(i,idensity)-(slope*this%wgrid(i))
                      this%hb%mode(jdof)=(yr-yint)/slope
                   else
                      this%hb%mode(jdof)=this%wgrid(i)
                   end if
                end if
             end do
          end do
       end do
    end if

    !Bonella and Coker, JCP 122, 194102 (2005)
    !EQ 53: thermal dist = exp(-[tanh(beta*wj/2)/wj]*[(Pj**2)/2+((wj*Qj)**2)/2])
    !EQ 54: canonical dist= exp(-[tanh(beta*wj/2)/wj]
    !                    *[(Pj**2)/2+{(wj**2)/2}*{Qj-(cj/wj**2)}**2])

    != = = = = = = = = = = = = u n i t s = = = = = = = = = = = = =
    !sigp=mass*sqrt(hbar*omega)/sqrt(2*tanh(.5*beta*hbar*omega))
    !sigq=mass*sqrt(hbar)/sqrt(2*omega*tanh(.5*beta*hbar*omega))
    != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    !resample phase-space point
    Q0=0.
    do idof=1,this%hb%ndof
       sigp=sqrt(this%hb%mode(idof)&
            /(2._double*tanh(0.5_double*beta*this%hb%mode(idof))))/this%hb%rmass(idof)
       sigq=sigp/this%hb%mode(idof)
       this%hb%Q(idof)=Q0+sigq*gran()
       this%hb%P(idof)=sigp*gran()
    end do
    this%hb%Qmin=minval(this%hb%Q)
    this%hb%Qmax=maxval(this%hb%Q)

    !Eigenvec = identity (usually computed in update)
    this%hb%EigenVec=iden(this%hb%ndof)

    call update(this)
    if(check(this).EQ.1)call stop('harmonicbath_resample, object failed exit check')
  end subroutine harmonicbath_resample
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short) function harmonicbath_check(this)
    type(harmonicbath),intent(in)::this
    
    integer(short)::ierr
    integer(long)::ipt,loc(1)
    
    call Note('Begin harmonicbath_check.')

    harmonicbath_check=0
    
    !logical::initialized
    if(.not.this%initialized)then
       call Warn('harmonicbath_check: harmonicbath object not initialized.')
       harmonicbath_check=1
       return
    end if

    !type(hb)::hb
    if(check(this%hb).NE.0)then
       call Warn('harmonicbath_check: classical primitive failed check.')
       harmonicbath_check=1
       return
    end if

    !integer(long)::nbath
    if(this%nbath.LE.0)then
       call Warn('harmonicbath_check: number of baths is less than 0.')
       harmonicbath_check=1
       return
    end if
    if(this%nbath.NE.this%nbath)then
       call Warn('harmonicbath_check: number of baths is not a number.')
       harmonicbath_check=1
       return
    end if
    if(this%nbath.GE.huge(this%nbath))then
       call Warn('harmonicbath_check: number of baths is too large.')
       harmonicbath_check=1
       return
    end if

    !integer(long)::npt
    if(this%npt.LE.0)then
       call Warn('harmonicbath_check: number of spectral points is negative.')
       harmonicbath_check=1
       return
    end if
    if(this%npt.NE.this%npt)then
       call Warn('harmonicbath_check: number of spectral points is NAN.')
       harmonicbath_check=1
       return
    end if
    if(this%npt.GE.huge(this%npt))then
       call Warn('harmonicbath_check: number of spectral points is too large.')
       harmonicbath_check=1
       return
    end if

    !integer(long)::,nden
    if(this%nden.NE.this%nden)then
       call Warn('harmonicbath_check: number of spectral densities is NAN.')
       harmonicbath_check=1
       return
    end if
    if(this%nden.LE.0)then
       call Warn('harmonicbath_check: number of spectral densities is negative.')
       harmonicbath_check=1
       return
    end if
    if(this%nden.GE.huge(this%nden))then
       call Warn('harmonicbath_check: number of spectral densities is huge.')
       harmonicbath_check=1
       return
    end if

    !real(double),dimension(:),pointer::wgrid
    if(.not.associated(this%wgrid))then
       call Warn('harmonicbath_check: wgrid memory not associated')
       harmonicbath_check=1
       return
    end if
    if(size(this%wgrid).ne.this%npt)then
       call Warn('harmonicbath_check: size of spectral grid (wgrid) array'&
            //' is not equal to number of spectral data points')
       harmonicbath_check=1
       return
    end if
    if(any(this%wgrid.NE.this%wgrid))then
       call Warn('harmonicbath_check: wgrid has NAN values')
       harmonicbath_check=1
       return
    end if
    if(any(abs(this%wgrid).GE.huge(this%wgrid)))then
       call Warn('harmonicbath_check: wgrid has huge values')
       harmonicbath_check=1
       return
    end if
    do ipt=this%npt,1,-1
       loc=maxloc(this%wgrid(1:ipt))
       if(loc(1).NE.ipt)then
          call Warn('harmonicbath_check: wgrid is unordered from min to max. val='&
               //trim(float2str(this%wgrid(ipt)))//' loc='//trim(int2str(ipt))&
               //' maxloc='//trim(int2str(loc(1))))
          harmonicbath_check=1
          return
       end if
    end do

    !real(double),dimension(:,:),pointer::Jw
    if(.not.associated(this%Jw))then
       call Warn('harmonicbath_check: Jw memory not associated')
       harmonicbath_check=1
       return
    end if
    if(size(this%Jw).ne.this%nden*(this%npt))then
       call Warn('harmonicbath_check: size of Jw matrix is not equal to'&
            //' the number of spectral data points per spectral density')
       harmonicbath_check=1
       return
    end if
    if(any(this%Jw.NE.this%Jw))then
       call Warn('harmonicbath_check: Jw has NAN values')
       harmonicbath_check=1
       return
    end if
    if(any(abs(this%Jw).GE.huge(this%Jw)))then
       call Warn('harmonicbath_check: Jw has huge values')
       harmonicbath_check=1
       return
    end if

    !real(double),dimension(:,:),pointer::intDOS
    if(.not.associated(this%intDOS))then
       call Warn('harmonicbath_check: intDOS memory not associated')
       harmonicbath_check=1
       return
    end if
    if(size(this%intDOS).ne.this%nden*(this%npt))then
       call Warn('harmonicbath_check: size of intDOS matrix is not equal to'&
            //' the number of spectral points per spectral density')
       harmonicbath_check=1
       return
    end if

    if(any(this%intDOS.NE.this%intDOS))then
       call Warn('harmonicbath_check: intDOS has NAN values')
       harmonicbath_check=1
       return
    end if
    if(any(abs(this%intDOS).GE.huge(this%intDOS)))then
       call Warn('harmonicbath_check: intDOS has huge values')
       harmonicbath_check=1
       return
    end if

    !real(double),dimension(:),pointer::lambda
    if(.not.associated(this%lambda))then
       call Warn('harmonicbath_check: lambda memory not associated')
       harmonicbath_check=1
       return
    end if
    if(size(this%lambda).ne.this%nden)then
       call Warn('harmonicbath_check: size of lambda array not equal to'&
            //' the number of spectral points')
       harmonicbath_check=1
       return
    end if
    if(any(this%lambda.NE.this%lambda))then
       call Warn('harmonicbath_check: lambda has NAN values.')
       harmonicbath_check=1
       return
    end if
    if(any(abs(this%lambda).GE.huge(this%lambda)))then
       call Warn('harmonicbath_check: lambda has huge values.')
       harmonicbath_check=1
       return
    end if
    if(any(this%lambda.LE.0._double))then
       call Warn('harmonicbath_check: lambda has non positive values.')
       harmonicbath_check=1
       return
    end if

    !integer(long),dimension(:),pointer::map
    if(.not.associated(this%map))then
       call Warn('harmonicbath_check: map memory not associated.')
       harmonicbath_check=1
       return
    end if
    if(size(this%map).ne.this%nbath)then
       call Warn('harmonicbath_check: size of map array not equal to number of'&
            //' the number of independent baths.')
       harmonicbath_check=1
       return
    end if
    if(any(this%map.NE.this%map))then
       call Warn('harmonicbath_check: map has NAN values')
       harmonicbath_check=1
       return
    end if
    if(any(abs(this%map).GE.huge(this%map)))then
       call Warn('harmonicbath_check: map has huge values')
       harmonicbath_check=1
       return
    end if
    
  end function harmonicbath_check

end module harmonicbath_class

