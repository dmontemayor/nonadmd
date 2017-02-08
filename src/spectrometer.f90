!-----------------------------------------------------------------------------!
!> \brief 
!! Spectrometric functionality 
!! \details
!! A spectrum type is defined here. The spectrum type is the standard output
!! of the associated methods with perform various spectrocopic mesurements of
!! the total Hamiltonian or subsystems.
!! \author
!! Daniel Montemayor
!! \date
!! May 2011
!<----------------------------------------------------------------------------!
!----------------     Change Log     -----------------------------------------!
! + ABspec does not change cs sampling method (wignerbath now obsolete)
! + ABspec saves intermediate files.
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
module spectrometer_class
  use type_kinds
  use ErrorLog
  use filemanager
  use MPIframework
  use string
  use rand
  use math
  use atomicunits
  use quantum_class
  use classical_class
  use coupling_class
  use textgraphs
  use molreader
  implicit none
  private 

  !> Standard spectrum type  
  type spectrum
     !> True if spectrum type is properly initialized
     logical::initialized=.false.
     !> number of points in the frequency domain
     integer(long)::npt
     !> minimum and maximum values of the domain
     real(double)::Emin,Emax
     !> Array of size npt containing the observed spectrum as a function of frequency.
     real(double),dimension(:),pointer::spectrum
  end type spectrum

  !> Machine limit of integer values. Used to exit out of non-converging loops.
  integer(long),parameter::killint=huge(killint)
  
  public::observe,display,check,kill
  public::spectrum
  
  !> Measures the Hamiltonian
  interface observe
     module procedure get_spectrum
  end interface

  !> Display the spectrum type 
  interface display
     module procedure display_spectrum
  end interface

  !> Function that checks the spectrum type
  interface check
     module procedure spectrum_check
  end interface

  !> Destroys the spectrum type
  interface kill
     module procedure spectrum_kill
  end interface
contains

  !--------------------------------------------------------------------------------
  !> \brief General call to compute verious spectra.
  !! \param[inout] this spectrum type, it is overwritten with the spectroscopic measurement on return.
  !! \param[inout] qs general quantum subsystem to be measured.
  !! \param[inout] cs general classical subsystem to be measured.
  !! \param[inout] cp general coupling term to be measured.
  !! \param[in] npt number of points in the frequency domain.
  !! \param[in] Emin minimum value of frequncy domain.
  !! \param[in] Emax maximum value of the frequency domain.
  !! \param[in] type string containing the type of the spectroscopic measurement. Options are 'AB' and 'ABPDB'.
  !! \param[in] samples minimum number of measurements to make before convergance protocols take effect. Guarantees good statistics.
  !! \param[in] tol \f$ R^2 \f$ difference between successive measurements. Used to determine spectrum convergance.
  !! \param[in] file name of file to output spectrum on return.
  !! \param[in] pdb name of input pdb file.
  !! \param[in] psf name of input psf file.
  !! \param[in] psfx name of psf file containing excited state charges.
  !! \param[in] atm_select atom type used to determine a subset of atoms.
  !! \param[in] dipole_to atom type specifing direction a dipole points to.
  !! \param[in] dipole_from atom type specifing direction a dipole points from.
  !! \param[in] dipole_center atom type specifing the origin of a dipole.
  !<-------------------------------------------------------------------------------
  subroutine get_spectrum(this,qs,cs,cp,Emin,Emax,npt,type,tol,samples,file,pdb,psf,psfx,atm_select,dipole_to,dipole_from,dipole_center)
    type(spectrum),intent(inout)::this
    type(quantum),intent(inout),target::qs
    type(classical),intent(inout),target,optional::cs
    type(coupling),intent(inout),target,optional::cp
    real(double),intent(in)::Emin,Emax
    integer(long),intent(in)::npt
    character(len=*),intent(in)::type

    integer(long),intent(in),optional::samples
    real(double),intent(in),optional::tol
    character*(*),intent(in),optional::file
    character*(*),intent(in),optional::pdb
    character*(*),intent(in),optional::psf
    character*(*),intent(in),optional::psfx
    type(atom),intent(in),optional::atm_select(:),dipole_to(:),dipole_from(:),dipole_center(:)

    integer(long)::nsample
    real(double)::tol0
    character(len=path)::outfile='AB.spectrum'

    nsample=100
    if(present(samples))nsample=samples
    tol0=1E-3
    if(present(tol))tol0=tol
    if(present(file))outfile=adjustl(file)

    !setup energy domain and allocate spectrum
    this%npt=npt
    this%Emax=Emax
    this%Emin=Emin
    if(associated(this%spectrum))nullify(this%spectrum)
    allocate(this%spectrum(0:this%npt-1))
    this%spectrum=0_double
    this%initialized=.true.
    if(check(this).NE.0)call Stop('spectrum: failed check')

    select case (trim(type))
    case('AB')
       if(present(cs).and.present(cp))then
          call get_ABspectrum(this,qs,cs,cp,tol0,nsample,file=trim(outfile))
       else
          call warn('get_spectrum: cannot compute AB spectrum, missing classical subsystem or coupling term.','exiting.')
       end if
    case('ABPDB')
       call get_ABPDBspectrum(this,qs,pdb=trim(pdb),psf=trim(psf),psfx=trim(psfx),file=trim(outfile),atm_select=atm_select,dipole_to=dipole_to,dipole_from=dipole_from,dipole_center=dipole_center)
    case default
       write(*,*)'Unknown spectrum type!'
    end select

  end subroutine get_spectrum
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !> \brief nullifies spectrum pointer
  !! \param [inout] this spectrum type
  !<=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine spectrum_kill(this)
    type(spectrum),intent(inout)::this

    call Note('Begin spectrum_kill')
    if(associated(this%spectrum))nullify(this%spectrum)
    this%initialized=.false.
  end subroutine spectrum_kill
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !> \breif Checks spectrum type
  !<=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  integer(short) function spectrum_check(this)
    type(spectrum),intent(in)::this
    

    call Note('Checking spectrum.')

    spectrum_check=0
    
    !logical::initialized
    if(.not.this%initialized)then
       call Warn('spectrum_check: spectrum not initialized.')
       spectrum_check=1
       return
    end if

    !integer(long)::npt
    if(this%npt.LE.0)then
       call Warn('spectrum_check: number of points is less than 0.')
       spectrum_check=1
       return
    end if
    if(this%npt.NE.this%npt)then
       call Warn('spectrum_check: number of points is not a number.')
       spectrum_check=1
       return
    end if
    if(this%npt.GE.huge(this%npt))then
       call Warn('spectrum_check: number of points is huge.')
       spectrum_check=1
       return
    end if

    !real(double)::Emin
    if(this%Emin.NE.this%Emin)then
       call Warn('spectrum_check: Emin is not a number.')
       spectrum_check=1
       return
    end if
    if(this%Emin.GE.huge(this%Emin))then
       call Warn('spectrum_check: Emin is huge.')
       spectrum_check=1
       return
    end if

    !real(double)::Emax
    if(this%Emax.NE.this%Emax)then
       call Warn('spectrum_check: Emax is not a number.')
       spectrum_check=1
       return
    end if
    if(this%Emax.GE.huge(this%Emax))then
       call Warn('spectrum_check: Emax is huge.')
       spectrum_check=1
       return
    end if

    !real(double),dimension(:),pointer::spectrum
    if(.not.associated(this%spectrum))then
       call Warn('spectrum_check: spectrum memory not associated')
       spectrum_check=1
       return
    end if
    if(size(this%spectrum).ne.this%npt)then
       call Warn('spectrum_check: size of spectrum array not equal to number of points')
       spectrum_check=1
       return
    end if
    if(any(this%spectrum.NE.this%spectrum))then
       call Warn('spectrum_check: spectrum has NAN values')
       spectrum_check=1
       return
    end if
    if(any(this%spectrum.GE.huge(this%spectrum)))then
       call Warn('spectrum_check: spectrum has huge values')
       spectrum_check=1
       return
    end if

  end function spectrum_check
  !-----------------------------------------------------------  
  subroutine display_spectrum(this)
    type(spectrum),intent(inout)::this

    real(double)::Xset(this%npt),Yset(this%npt)
    integer(long)::i
    real(double)::dw

    dw=(this%Emax-this%Emin)/real(this%npt)
    do i=1,this%npt
       Xset(i)=(this%Emin+(i-1)*dw)/invcm
       Yset(i)=this%spectrum(i)
    end do
    call display(Xset,Yset,nptx=line)
    write(*,*)'Frequency (x axis) in (1/cm)'

  end subroutine display_spectrum
  !-----------------------------------------------------------  

  !--------------------------------------------------------------------------------
  !> \brief Computes absorption spectrum using groundstate to adiabatic state transition dipole.
  !! \details Requires a quantum subsystem described in terms of transistion dipole moments.
  !! N points in spectrum with energy domain from Emin to Emax.
  !! Integrate over thermal realizations to compute spectrum.
  !! Recompute spectrum with twice as many thermal realizations.
  !! Keep doubling number of realizations until sum of difference in spectrum to 
  !! previous spectrum points squared is less than some tolerance tol.
  !! \param[inout] this spectrum type, it is overwritten with the spectroscopic measurement on return.
  !! \param[inout] qs general quantum subsystem to be measured.
  !! \param[inout] cs general classical subsystem to be measured.
  !! \param[inout] cp general coupling term to be measured.
  !! \param[in] samples minimum number of measurements to make before convergance protocols take effect. Guarantees good statistics.
  !! \param[in] tol \f$ R^2 \f$ difference between successive measurements. Used to determine spectrum convergance.
  !! \param[in] file name of file to output spectrum on return.
  !! \todo
  !! * Add method to check if quantum subsystem is supported.
  !! * Rewrite details section and add brief comment. 
  !< ---------------------------------------------------------------------------------
  subroutine get_ABspectrum(this,qs,cs,cp,tol,samples,file)
    type(spectrum),intent(inout)::this
    type(quantum),intent(inout),target::qs
    type(classical),intent(inout),optional::cs
    type(coupling),intent(inout),optional::cp

    integer(long),intent(in),optional::samples
    real(double),intent(in),optional::tol
    character*(*),intent(in),optional::file

    real(double)::Emin,Emax

    real(double),parameter::factor=1.20_double

    integer(long)::nstate,nsample,ndim
    integer(short)::unit
    logical::pass,cnvrg,stable

    real(double)::w,dw,tr,R2,norm,thisnorm,tolerance
    complex(double),allocatable::H(:,:)
    real(double),pointer::Mu(:,:)
    real(double),allocatable::Eng(:),MM(:,:),spectrum(:)
    integer(long)::i,idim,istate,jstate,kstate,isample,lastsample

    pass=.true.
    !check quantum subsystem
    if(check(qs).NE.0)pass=.false.
    nstate=qs%hs%nstate

    !check supported quantum subsystem types for transition dipole information 
    if(associated(Mu))nullify(Mu)
    select case(qs%type)
    case('dipoles')
       ndim=qs%dipoles%ndim
       allocate(Mu(nstate,ndim))
       Mu=>qs%dipoles%mu
    case default
       call warn('get_ABspectrum: quantum subsystem type not supported.')
       pass=.false.
    end select

    if(present(cs))then
       if(check(cs).NE.0)pass=.false.
    end if
    if(present(cp))then
       if(check(cp).NE.0)pass=.false.
    end if
    
    if(pass)then

       stable=.true.
       call note('Computing Absorption Spectrum.')

       !setup energy domain and allocate spectrum
       Emax=this%Emax
       Emin=this%Emin
       dw=(Emax-Emin)/real(this%npt)
       if(allocated(spectrum))deallocate(spectrum)
       allocate(spectrum(0:this%npt-1))

       tolerance=1E-3
       if(present(tol))tolerance=tol*tol
       if(tolerance.GT.huge(tolerance))then
          call warn('get_ABspectrum: convergence tolerance (opt param [tol]) is huge.'& 
               ,'setting default to tol='//trim(float2str(huge(tolerance)))//'.')
          tolerance=huge(tolerance)
       end if
       if(tolerance.LT.epsilon(tolerance))then
          call warn('get_ABspectrum: convergence tolerance (opt param [tol]) is tiny.'& 
               ,'setting default to tol='//trim(float2str(epsilon(tolerance)))//'.')
          tolerance=epsilon(tolerance)
       end if

       !allocate EigenVectors and EigenValues
       if(allocated(H))deallocate(H)
       allocate(H(nstate,nstate))
       if(allocated(Eng))deallocate(Eng)
       allocate(Eng(nstate))

       !allocate and compute Mu dot Mu 
       if(allocated(MM))deallocate(MM)
       allocate(MM(nstate,nstate))
       

       !sum over realizations
       nsample=10
       if(present(samples))nsample=samples
       lastsample=0
       cnvrg=.false.
       this%spectrum=0.0_double
       spectrum=0.0_double
       do while(.not.cnvrg)
          do isample=lastsample+1,nsample
             !call note('sample '//trim(int2str(isample)))
             call resample(qs)
             if(present(cs))call resample(cs)
             if(present(cp))call resample(cp)
             
             H=cmplx(qs%hs%diabat)
             if(present(cp))H=H+cmplx(cp%hc%V)
             call diagonalize(nstate,H,Eng)
             Eng=Eng-qs%hs%Eg


             MM=0.0_double
             do istate=1,nstate
                do jstate=istate,nstate
                   do idim=1,ndim
                      MM(istate,jstate)=MM(istate,jstate)+mu(istate,idim)*mu(jstate,idim)
                   end do
                   MM(jstate,istate)=MM(istate,jstate)
                end do
             end do


             !sum over adiabats
             do kstate=1,nstate
                !sum over diabats
                tr=0.0_double
                do istate=1,nstate
                   do jstate=1,nstate
                      tr=tr+real(H(istate,kstate)*conjg(H(jstate,kstate)))&
                           *MM(istate,jstate)
                   end do
                end do                
                tr=Eng(kstate)*tr


                !build up spectrum 
                do i=0,this%npt-1
                   w=Emin+i*dw
                   if(Eng(kstate).GT.w-dw/2_double.and.Eng(kstate).LE.w+dw/2_double)spectrum(i)=spectrum(i)+tr
                end do
             end do
          end do

          !decide to keep or add more realizations 
          norm=sum(spectrum)
          if(norm.EQ.0.)then
             call warn('get_ABspectrum: spectrum has zero valued integral.'&
                  ,'exiting without computing spectrum.')
             stable=.false.
          end if
          R2=0.0_double
          if(stable)R2=sum((spectrum/real(nsample)-this%spectrum)**2)/real(this%npt)
          if(R2.LE.tolerance)cnvrg=.true.
          !call Note('ABspectrum: samples='//trim(int2str(nsample)&
          !     //' R='//trim(float2str(sqrt(R2)))))

          if(nsample.GT.killint)then
             call warn('get_ABspectrum: spectrum failed to converge after '&
                  //trim(int2str(killint))//' realizations.'&
                  ,'Stop trying to converge and continue as is.')
             cnvrg=.true.
          end if
          
          !save spectrum as old spectrum and double number of samples
          this%spectrum=spectrum/real(nsample)
          lastsample=nsample
          nsample=int(nsample*factor)

          if(present(file))then
             call Note('Saving intermediate ABspec')
             unit=newunit()
             open(unit,file=file)
             norm=maxval(abs(this%spectrum))
             do i=0,this%npt-1
                w=Emin+i*dw
                write(unit,"(10(ES18.10E2,1X))")w,this%spectrum(i)/norm
             end do
             close(unit)
          end if

       end do
       nsample=lastsample
       !call Note('ABspectrum: samples='//trim(int2str(nsample)&
       !     //' R_squared ='//trim(float2str(R2))))
       
       if(.not.stable)spectrum=0.0_double !incase of runtime error
       
       if(present(file))then
          call note('input file= '//file)

          norm=maxval(abs(this%spectrum))

          unit=newunit()
          open(unit,file=file)
          do i=0,this%npt-1
             w=Emin+i*dw
             write(unit,"(10(ES18.10E2,1X))")w,this%spectrum(i)/norm
          end do
          close(unit)

       end if
    else
       call warn('get_ABspectrum: failed check.','not computing absorption spectrum.')
    end if
    
    if(allocated(H))deallocate(H)
    if(allocated(Eng))deallocate(Eng)
    if(allocated(MM))deallocate(MM)
    if(allocated(spectrum))deallocate(spectrum)
    if(associated(Mu))nullify(Mu)


  end subroutine get_ABspectrum
!--------------------------------------------------------
!!$ subroutine get_SpectralDensity(this,qs,cs,cp,Emin,Emax,N,tol,samples,file)
!!$    !------------------- ------------------------------------------------------
!!$    !*Computes Spectal density from cp%DMBLcoupling%g factors.
!!$    !*Requires DMBL type coupiling.
!!$    !   N points in spectrum with energy domain from Emin to Emax.
!!$    !*Integrate over thermal realizations to compute spectrum.
!!$    !*Recompute spectrum with twice as many thermal realizations.
!!$    !*Keep doubling number of realizations until sum of dif in spectrum to 
!!$    !    previous spectrum points squared is less than some tolerance tol.
!!$    !---------------------------
!!$    !TODO:
!!$    !add method to check if coupling subsystem is supported 
!!$    !---------------------------------------------------------------------------------
!!$    type(spectrum),intent(inout)::this
!!$    type(quantum),intent(inout),target::qs
!!$    type(classical),intent(inout),target::cs
!!$    type(coupling),intent(inout),target::cp
!!$    real(double),intent(in)::Emin,Emax
!!$    integer(long),intent(in),optional::N,samples
!!$    real(double),intent(in),optional::tol
!!$    character*(*),intent(in),optional::file
!!$
!!$    real(double),parameter::factor=1.20_double
!!$
!!$    integer(long)::nbath,nsample
!!$    integer(short)::unit
!!$    logical::pass,cnvrg,stable
!!$
!!$    real(double)::w,dw,tolerance,g,wj,rmass,R2,norm,Jw
!!$    integer(long)::lastsample,isample,i,idof,ibath
!!$    real(double),allocatable::spectrum(:)
!!$
!!$    pass=.true.
!!$
!!$    if(check(qs).NE.0)pass=.false.
!!$
!!$    if(check(cs).NE.0)pass=.false.
!!$
!!$    if(check(cp).NE.0)pass=.false.
!!$
!!$    !check for supported coupling types 
!!$    select case(cp%type)
!!$    case('DMBLcoupling')
!!$    case default
!!$       call warn('get_SpectralDensity: coupling type not supported.')
!!$       pass=.false.
!!$    end select
!!$    nbath=cp%DMBLcoupling%nbath
!!$
!!$
!!$    if(pass)then
!!$       stable=.true.
!!$       call note('Computing Spectral Density.')
!!$       
!!$       !setup energy domain and allocate spectrum
!!$       this%npt=100
!!$       if(present(N))this%npt=N
!!$       if(this%npt.LT.1)then
!!$          call warn('get_spectraldensity: Number of spectrum points (opt param [N]) is less than 1.'&
!!$               ,'using default 100 point spectrum.')
!!$          this%npt=100
!!$       end if
!!$
!!$       dw=(Emax-Emin)/real(this%npt)
!!$       if(associated(this%spectrum))nullify(this%spectrum)
!!$       allocate(this%spectrum(0:this%npt-1))
!!$       if(allocated(spectrum))deallocate(spectrum)
!!$       allocate(spectrum(0:this%npt-1))
!!$
!!$
!!$       tolerance=1E-3
!!$       if(present(tol))tolerance=tol*tol
!!$       if(tolerance.GT.huge(tolerance))then
!!$          call warn('get_SpectralDensity: convergence tolerance (opt param [tol]) is huge.'& 
!!$               ,'setting default to tol='//trim(float2str(huge(tolerance)))//'.')
!!$          tolerance=huge(tolerance)
!!$       end if
!!$       if(tolerance.LT.epsilon(tolerance))then
!!$          call warn('get_SpectralDensity: convergence tolerance (opt param [tol]) is tiny.'& 
!!$               ,'setting default to tol='//trim(float2str(epsilon(tolerance)))//'.')
!!$          tolerance=epsilon(tolerance)
!!$       end if
!!$
!!$       !sum over realizations
!!$       nsample=10
!!$       if(present(samples))nsample=samples
!!$       lastsample=0
!!$       cnvrg=.false.
!!$       this%spectrum=0.0_double
!!$       spectrum=0.0_double
!!$       do while(.not.cnvrg)
!!$          do isample=lastsample+1,nsample
!!$             !call note('sample '//trim(int2str(isample)))
!!$             call resample(qs)
!!$             call resample(cs)
!!$             call resample(cp)
!!$
!!$             !compute for each independent bath
!!$             do ibath=1,1!nbath
!!$                !sum over modes
!!$                do idof=1,cs%hb%ndof
!!$                   if(cp%DMBLcoupling%map(idof).EQ.ibath)then
!!$                      g=cp%DMBLcoupling%g(idof,ibath)
!!$                      wj=cs%hb%mode(idof)
!!$                      rmass=cs%hb%rmass(idof)
!!$                      !build up spectrum 
!!$                      do i=0,this%npt-1
!!$                         w=Emin+i*dw
!!$                         if(wj.GT.w-dw/2_double.and.wj.LE.w+dw/2_double)&
!!$                              spectrum(i)=spectrum(i)+rmass*g**2/wj
!!$                      end do
!!$                   end if
!!$                end do
!!$             end do
!!$          end do
!!$
!!$          !decide to keep or add more realizations 
!!$          norm=sum(spectrum)
!!$          if(norm.EQ.0.)then
!!$             call warn('get_SpectralDensity: spectrum has zero valued integral.'&
!!$                  ,'exiting without computing spectrum.')
!!$             stable=.false.
!!$          end if
!!$
!!$          R2=0.0_double
!!$          if(stable)R2=sum((spectrum/real(nsample)-this%spectrum)**2)/real(this%npt)
!!$          if(R2.LE.tolerance)cnvrg=.true.
!!$          call Note('SpectralDensity: samples='//trim(int2str(nsample)&
!!$               //' R='//trim(float2str(sqrt(R2)))))
!!$
!!$          if(nsample.GT.killint)then
!!$             call warn('get_SpectraDensity: spectrum failed to converge after '&
!!$                  //trim(int2str(killint))//' realizations.'&
!!$                  ,'Stop trying to converge and continue as is.')
!!$             cnvrg=.true.
!!$          end if
!!$          
!!$          !save spectrum as old spectrum and double number of samples
!!$          this%spectrum=spectrum/real(nsample)
!!$          lastsample=nsample
!!$          nsample=int(nsample*factor)
!!$
!!$          if(present(file))then
!!$             call Note('Saving intermediate SpectralDensity')
!!$             unit=newunit()
!!$             open(unit,file=file)
!!$             do i=0,this%npt-1
!!$                w=Emin+i*dw
!!$                write(unit,"(10(ES18.10E2,1X))")w,this%spectrum(i)*halfpi/dw
!!$             end do
!!$             close(unit)
!!$          end if
!!$       end do
!!$       nsample=lastsample
!!$       call Note('Spectral density: samples='//trim(int2str(nsample)&
!!$            //' R_squared ='//trim(float2str(R2))))
!!$       
!!$       if(.not.stable)spectrum=0.0_double !incase of runtime error
!!$       
!!$       
!!$    else
!!$       call warn('get_SpectralDensity: failed check.','not computing spectrum.')
!!$    end if
!!$    
!!$  end subroutine get_SpectralDensity
!--------------------------------------------------------

!!$  subroutine get_StimEM(this,qs,cs,cp,Emin,Emax,N,tol,file)
!!$    !--------------------------------------------------------------------------------
!!$    !*Computes change in Energy of the field per unit time due to stimulated
!!$    !emmission of initial excited state populations to groundstate transitions.
!!$    !*Requires a quantum subsystem described in terms of transistion dipole moments.
!!$    !N points in spectrum with energy domain from Emin to Emax.
!!$    !*Integrate over thermal realizations to compute spectrum.
!!$    !*Compute spectrum with twice as many thermal realizations.
!!$    !*Keep doubling number of realizations untill sum of difference
!!$    !in spectrum to previous spectrum points squared is less than
!!$    !some tolerance tol.
!!$    !---------------------------
!!$    !TODO:
!!$    !add method to check if quantum subsystem is supported 
!!$    !---------------------------------------------------------------------------------
!!$    type(StimEM),intent(inout)::this
!!$    type(quantum),intent(inout)::qs
!!$    type(classical),intent(inout),optional::cs
!!$    type(coupling),intent(inout),optional::cp
!!$    real(double),intent(in)::Emin,Emax
!!$    integer(long),intent(in),optional::N
!!$    real(double),intent(in),optional::tol
!!$    character*(*),intent(in),optional::file
!!$
!!$    integer(long)::nstate,nsample,ndim
!!$    integer(short)::unit
!!$    logical::usedunit,pass,cnvrg
!!$
!!$    real(double)::w,dw,tr,Mk,Pk,norm,thisnorm,tolerance,Eg
!!$    complex(double),allocatable::H(:,:)
!!$    real(double),allocatable::Eng(:),MM(:,:),Mu(:,:),spectrum(:)
!!$    integer(long)::i,idim,istate,jstate,kstate,isample,lastsample
!!$
!!$    pass=.true.
!!$    !check quantum subsystem is initialized
!!$    if(.not.qs%initialized)then
!!$       write(*,*)'Warning in get_StimEM: quantum subsystem not initialized'
!!$       write(*,*)'Exiting get_StimEM without doing anything'
!!$       pass=.false.
!!$    else
!!$       nstate=qs%hs%nstate
!!$    end if
!!$
!!$
!!$    !check if quantum subsystem type is supported
!!$    if(.not.associated(qs%disordereddipole_v1r0))then
!!$       write(*,*)'Warning in get_StimEM: quantum subsystem type not supported'
!!$       write(*,*)'Exiting get_StimEM without doing anything'
!!$       pass=.false.
!!$    else
!!$       ndim=qs%disordereddipole_v1r0%ndim
!!$       if(allocated(Mu))deallocate(Mu)
!!$       allocate(Mu(nstate,ndim))
!!$       Mu=qs%disordereddipole_v1r0%mu
!!$    end if
!!$
!!$    if(present(cs))then
!!$       if(.not.cs%initialized)then
!!$          write(*,*)'Warning in get_StimEM: classical subsystem not initialized'
!!$          write(*,*)'Exiting get_StimEM without doing anything'
!!$          pass=.false.
!!$       end if
!!$    end if
!!$    if(present(cp))then
!!$       if(.not.cp%initialized)then
!!$          write(*,*)'Warning in get_StimEM: coupling subsystem not initialized'
!!$          write(*,*)'Exiting get_StimEM without doing anything'
!!$          pass=.false.
!!$       end if
!!$    end if
!!$    
!!$    if(pass)then
!!$       
!!$       write(*,*)'Computing Stimulated Emission Spectrum.'
!!$
!!$       !setup energy domain and allocate spectrum
!!$       this%npt=100
!!$       if(present(N))this%npt=N
!!$       if(this%npt.LT.1)then
!!$          write(*,*)'Warning in get StimEM: Number of spectrum points (opt param [N]) is less than 1.' 
!!$          write(*,*)'Program will default to 100 spectrum points.'
!!$          this%npt=100
!!$       end if
!!$       dw=(Emax-Emin)/real(this%npt)
!!$       if(associated(this%spectrum))nullify(this%spectrum)
!!$       allocate(this%spectrum(0:this%npt-1))
!!$       if(allocated(spectrum))deallocate(spectrum)
!!$       allocate(spectrum(0:this%npt-1))
!!$       Eg=qs%hs%Eg
!!$       
!!$       tolerance=.01
!!$       if(present(tol))tolerance=tol
!!$       if(tolerance.GT.1E15)then
!!$          write(*,*)'Warning in get StimEM: convergence tolerance (opt param [tol]) is too large.' 
!!$          write(*,*)'Program will default to tol=1E15.'
!!$          tolerance=1E-15
!!$       end if
!!$       if(tolerance.LT.1E-15)then
!!$          write(*,*)'Warning in get StimEM: convergence tolerance (opt param [tol]) is too small.' 
!!$          write(*,*)'Program will default to tol=1E-15.'
!!$          tolerance=1E-15
!!$       end if
!!$
!!$       !allocate EigenVectors and EigenValues
!!$       if(allocated(H))deallocate(H)
!!$       allocate(H(nstate,nstate))
!!$       if(allocated(Eng))deallocate(Eng)
!!$       allocate(Eng(nstate))
!!$
!!$       !allocate and compute Mu dot Mu 
!!$       if(allocated(MM))deallocate(MM)
!!$       allocate(MM(nstate,nstate))
!!$       MM=0.0_double
!!$       do istate=1,nstate
!!$          do jstate=istate,nstate
!!$             do idim=1,ndim
!!$                MM(istate,jstate)=MM(istate,jstate)+mu(istate,idim)*mu(jstate,idim)
!!$             end do
!!$             MM(jstate,istate)=MM(istate,jstate)
!!$          end do
!!$       end do
!!$       
!!$
!!$       !sum over realizations
!!$       nsample=1
!!$       lastsample=0
!!$       cnvrg=.false.
!!$       this%spectrum=0.0_double
!!$       spectrum=0.0_double
!!$       do while(.not.cnvrg)
!!$          do isample=lastsample+1,nsample
!!$             !write(*,*)'sample',isample
!!$             call update(qs)
!!$             if(present(cs))call resample(cs)
!!$             if(present(cp))call update(cp)
!!$             H=cmplx(qs%hs%diabat)
!!$             if(present(cp))H=H+cmplx(cp%hc%V)
!!$             call diagonalize(nstate,H,Eng)
!!$
!!$             !sum over adiabats
!!$             do kstate=1,nstate
!!$
!!$                !compute dipole contribution
!!$                Mk=0.0_double
!!$                do istate=1,nstate
!!$                   do jstate=1,nstate
!!$                      !Mk=Mk+real(H(kstate,istate)*conjg(H(kstate,jstate)))&
!!$                      Mk=Mk+real(H(istate,kstate)*conjg(H(jstate,kstate)))&
!!$                           *MM(istate,jstate)
!!$                   end do
!!$                end do
!!$                
!!$                !compute population contribution
!!$                !Pk=real(sum(H(kstate,:)*qs%hs%psi0)*sum(conjg(H(kstate,:)*qs%hs%psi0)))
!!$                !write(*,*)kstate,Pk
!!$
!!$                Pk=0.0_double
!!$                do istate=1,nstate
!!$                   do jstate=1,nstate
!!$                      !Pk=Pk+real(H(kstate,istate)*conjg(H(kstate,jstate)))&
!!$                      Pk=Pk+real(H(istate,kstate)*conjg(H(jstate,kstate)))&
!!$                           *qs%hs%psi0(istate)*qs%hs%psi0(jstate)
!!$                   end do
!!$                end do
!!$
!!$                !write(*,*)kstate,Pk
!!$                !write(*,*)
!!$
!!$                !satisfy resonance condition
!!$                tr=Eng(kstate)*Pk*Mk
!!$                do i=0,this%npt-1
!!$                   w=(Emin-Eg)+i*dw
!!$                   if(Eng(kstate).GT.w.and.Eng(kstate).LE.w+dw)spectrum(i)=spectrum(i)+tr
!!$                end do
!!$
!!$             end do
!!$!stop
!!$          end do
!!$
!!$          !decide to keep or add more realizations 
!!$          if(nsample.GT.1)then!compare normalized spectrumtra
!!$             norm=sum(spectrum)
!!$             thisnorm=sum(this%spectrum)
!!$             !tr=sqrt(sum((spectrum/norm-this%spectrum/thisnorm)**2)/real(this%npt))
!!$             tr=sum(abs(spectrum/norm-this%spectrum/thisnorm))/real(this%npt)
!!$             if(tr.LE.tolerance)cnvrg=.true.
!!$
!!$             if(nsample.GT.killint)then
!!$                write(*,*)'Warning in get_StimEM: spectrum failed to converge after '&
!!$                     //trim(int2str(killint))//' realizations.'
!!$                write(*,*)'Programm stop trying to converge and continue as is.'
!!$                cnvrg=.true.
!!$             end if
!!$             write(*,*)'samples',nsample,tr,cnvrg
!!$
!!$          end if
!!$
!!$          !save spectrum as old spectrum and double number of samples
!!$          this%spectrum=spectrum
!!$          lastsample=nsample
!!$          nsample=nsample*2
!!$
!!$       end do
!!$       this%spectrum=this%spectrum/real(lastsample)
!!$       !this%spectrum=this%spectrum/sum(this%spectrum)
!!$              
!!$       if(present(file))then
!!$          unit=1000
!!$          usedunit=.true.
!!$          do while (usedunit)
!!$             unit=unit+1
!!$             inquire(unit,opened=usedunit)
!!$          end do
!!$          open(unit,file=file)
!!$          norm=maxval(abs(this%spectrum))
!!$          do i=0,this%npt-1
!!$             w=(i*dw+Emin)/invcm
!!$             write(unit,"(10(ES18.10E2,1X))")w,this%spectrum(i)/norm
!!$          end do
!!$          close(unit)
!!$
!!$          write(*,*) "Stimulated Emission spectrum output in "//file
!!$
!!$       end if
!!$
!!$    end if
!!$    
!!$    if(allocated(H))deallocate(H)
!!$    if(allocated(Eng))deallocate(Eng)
!!$    if(allocated(MM))deallocate(MM)
!!$    if(allocated(spectrum))deallocate(spectrum)
!!$    if(allocated(Mu))deallocate(Mu)
!!$
!!$  end subroutine get_StimEM
!!$    !--------------------------------------------------------------------------------
!!$  subroutine get_DOS(this,qs,cs,cp,file)
!!$    !--------------------------------------------------------------------------------
!!$    !*Computes Density of adiabatic state spectrum
!!$    !*Does not include Bath energy. Bath serves to modulate state energies through
!!$    ! coupling object. 
!!$    !---------------------------
!!$    !TODO:
!!$    ! double samples until converged
!!$    !---------------------------------------------------------------------------------
!!$    type(DOS),intent(inout)::this
!!$    type(quantum),intent(inout),optional::qs
!!$    type(classical),intent(inout),optional::cs
!!$    type(coupling),intent(inout),optional::cp
!!$    character*(*),intent(in),optional::file
!!$
!!$    integer(long)::nstate
!!$    integer(short)::unit
!!$    logical::usedunit,pass
!!$
!!$    integer(long),parameter::nsteps=40000,EQsteps=100
!!$    real(double)::Emin=0.0_double
!!$    real(double)::Emax=0.0_double
!!$    real(double)::w,dw,Ewn
!!$    complex(double),allocatable::H(:,:)
!!$    real(double),allocatable::Eng(:)
!!$    integer(long)::i,j,istate,jstate
!!$
!!$    this%npt=100
!!$    if(associated(this%spectrum))nullify(this%spectrum)
!!$    allocate(this%spectrum(0:this%npt-1))
!!$    pass=.true.
!!$
!!$    nstate=1
!!$    if(present(qs))then
!!$       if(.not.qs%initialized)then
!!$          write(*,*)'Warning in get_DOS: quantum subsystem not initialized'
!!$          write(*,*)'Exiting get_DOS without doing anything'
!!$          pass=.false.
!!$       else
!!$          nstate=qs%hs%nstate
!!$       end if
!!$    end if
!!$    if(present(cs))then
!!$       if(.not.cs%initialized)then
!!$          write(*,*)'Warning in get_DOS: classical subsystem not initialized'
!!$          write(*,*)'Exiting get_DOS without doing anything'
!!$          pass=.false.
!!$       end if
!!$    end if
!!$
!!$
!!$    if(present(cp))then
!!$       if(.not.cp%initialized)then
!!$          write(*,*)'Warning in get_DOS: coupling subsystem not initialized'
!!$          write(*,*)'Exiting get_DOS without doing anything'
!!$          pass=.false.
!!$       end if
!!$    end if
!!$    
!!$    if(pass)then
!!$
!!$       write(*,*)'Computing Density of States Spectrum.'
!!$       
!!$       if(allocated(H))deallocate(H)
!!$       allocate(H(nstate,nstate))
!!$       if(allocated(Eng))deallocate(Eng)
!!$       allocate(Eng(nstate))
!!$
!!$       write(*,*)'Initializing energy domain.'
!!$       do i=1,EQsteps !get energy domain
!!$          if(present(qs))call update(qs)
!!$          if(present(cs))call resample(cs)
!!$          if(present(cp))call update(cp)
!!$          H=0.0_double
!!$          if(present(qs))H=H+cmplx(qs%hs%diabat)
!!$          if(present(cp))H=H+cmplx(cp%hc%V)
!!$          call diagonalize(nstate,H,Eng)
!!$          if(i.EQ.1)Emin=minval(Eng)
!!$          if(i.EQ.1)Emax=maxval(Eng)
!!$          if(minval(Eng).LT.Emin)Emin=minval(Eng)
!!$          if(maxval(Eng).GT.Emax)Emax=maxval(Eng) 
!!$      end do
!!$
!!$       !add a bufffer to domain
!!$       dw=(Emax-Emin)/10.0_double
!!$       Emax=Emax+dw
!!$       Emin=Emin-dw
!!$       dw=(Emax-Emin)/real(this%npt)
!!$
!!$       this%spectrum=0.0_double
!!$       do i=1,Nsteps !build histogram
!!$          !write(*,*)'sample',i
!!$          if(present(qs))call update(qs)
!!$          if(present(cs))call resample(cs)
!!$          if(present(cp))call update(cp)
!!$          H=0.0_double
!!$          if(present(qs))H=H+cmplx(qs%hs%diabat)
!!$          if(present(cp))H=H+cmplx(cp%hc%V)
!!$          call diagonalize(nstate,H,Eng)
!!$          do j=0,this%npt-1
!!$             w=j*dw+Emin
!!$             do istate=1,nstate
!!$                if(w.LE.Eng(istate).and.w+dw.GT.Eng(istate))&
!!$                     this%spectrum(j)=this%spectrum(j)+1.0_double              
!!$             end do
!!$          end do
!!$       end do
!!$       this%spectrum=this%spectrum/sum(this%spectrum)
!!$       if(allocated(H))deallocate(H)
!!$       if(allocated(Eng))deallocate(Eng)
!!$  
!!$       if(present(file))then
!!$          unit=1000
!!$          usedunit=.true.
!!$          do while (usedunit)
!!$             unit=unit+1
!!$             inquire(unit,opened=usedunit)
!!$          end do
!!$          open(unit,file=file)
!!$          
!!$          do j=0,this%npt-1
!!$             w=j*dw+Emin
!!$             write(unit,"(2(ES18.10E2,1X))")w/invcm,this%spectrum(j)/maxval(this%spectrum)
!!$          end do
!!$          close(unit)
!!$       end if
!!$
!!$       write(*,*) "density of states spectrum output in "//file
!!$
!!$    end if
!!$    
!!$  end subroutine get_DOS
!!$!--------------------------------------------------------

  subroutine get_ABPDBspectrum(this,qs,file,pdb,psf,psfx,atm_select,dipole_to,dipole_from,dipole_center)
!--------------------------------------------------------------------------------
!> \brief
!! Computes Absorption spectrum using a pdb trajectory and a defined transition dipole.
!! \todo 
!! * Use psfx to define transition charges.
!! * Add params section.
!< ---------------------------------------------------------------------------------
    type(spectrum),intent(inout)::this
    type(quantum),intent(inout),target::qs
    character*(*),intent(in)::pdb
    character*(*),intent(in)::psf
    character*(*),intent(in)::psfx
    type(ATOM),intent(in)::atm_select(:),dipole_to(:),dipole_from(:),dipole_center(:)
    character*(*),intent(in),optional::file

    real(double)::Emin,Emax

    integer(long)::nstate,ndim
    integer(short)::unit
    logical::pass

    real(double)::w,dw,tr,R2,norm,thisnorm
    complex(double),allocatable::H(:,:)
    real(double),allocatable::Mu(:,:)
    real(double),allocatable::Eng(:),MM(:,:),spectrum(:)
    integer(long)::i,idim,istate,jstate,kstate,isample,lastsample

    type(ATOM)::atm
    type(psfATOM)::psfatm
    type(psfATOM),allocatable::psfxatm(:)
    
    integer(long)::atomto,atomfrom,atomcenter
    integer(long)::natom,npsfatom,npsfxatom,iatom,jatom,nmatch
    logical::match,match_res,match_resname,match_segname,match_name
    !integer(long)::nframe,nmatch,natom,iatom,jatom,katom,npsfatom,nchromatm

    logical,allocatable::mask(:)
    real(double),allocatable::XYZ(:,:),Qg(:),Qx(:)
    real*8::Dg,Dx,Dmag,Tmag


    pass=.true.
    !check quantum subsystem
    if(check(qs).NE.0)pass=.false.
    nstate=qs%hs%nstate
    ndim=3

    if(size(atm_select).ne.nstate)then
       call warn('ABPDBspectrum: size of atm_select must equal number of states.','exiting.')
       return
    end if
    if(size(dipole_to).ne.nstate)then
       call warn('ABPDBspectrum: size of dipole_to must equal number of states.','exiting.')
       return
    end if
    if(size(dipole_from).ne.nstate)then
       call warn('ABPDBspectrum: size of dipole_from must equal number of states.','exiting.')
       return
    end if
    if(size(dipole_center).ne.nstate)then
       call warn('ABPDBspectrum: size of dipole_center must equal number of states.','exiting.')
       return
    end if

    if(allocated(Mu))deallocate(Mu)
    allocate(Mu(nstate,3))




!!$    
!!$    
!!$    !count atoms
!!$    call open_pdb(pdb)
!!$    natom=0
!!$    do while(trim(Next_Record()).NE.'END')
!!$       if(trim(Next_Record()).EQ.'ATOM')then
!!$          call read(atm)
!!$          !call display(atm)
!!$          natom=natom+1
!!$       else
!!$          call note('Record '//Next_Record()//' will get skiped')
!!$          call advance_record
!!$       end if
!!$    end do
!!$    call note('PDB natom = '//trim(int2str(natom)))
!!$
!!$    !allocate arrays
!!$    call note('ABPDB: Allocating dynamic arrays.')
!!$    if(allocated(mask))deallocate(mask)
!!$    allocate(mask(natom))
!!$    if(allocated(XYZ))deallocate(XYZ)
!!$    allocate(XYZ(3,natom))
!!$    if(allocated(Qg))deallocate(Qg)
!!$    allocate(Qg(natom))
!!$    if(allocated(Qx))deallocate(Qx)
!!$    allocate(Qx(natom))
!!$
!!$    !load the ground state psf and skip to atoms section
!!$    call open_psf(psf)
!!$    call advance_psf(section='atoms',nrecord=npsfatom)
!!$
!!$    !check for the same number of atoms
!!$    if(npsfatom.NE.natom)then
!!$       call warn('ABPDBspectrum: Number of atoms in psf does not equall the number of atoms in the pdb.','exiting.')
!!$       pass=.false.
!!$    else
!!$       !read psf atoms
!!$       do iatom=1,npsfatom
!!$          call read(psfatm)
!!$          !call display(psfatm)
!!$          Qg(iatom)=psfatm%charge
!!$       end do
!!$    end if
!!$    !close ground state psf
!!$    call close_psf
!!$
!!$    !load the excited state psf and skip to atoms section
!!$    call open_psf(psfx)
!!$    call advance_psf(section='atoms',nrecord=npsfxatm)
!!$
!!$    !allocate excited state charges
!!$    call note('ABPDBspectrum: Allocating excited state charge dynamic array.')
!!$    if(allocated(psfxatm))deallocate(psfxatm)
!!$    allocate(psfxatm(npsfxatm))
!!$
!!$    !read psfx atoms
!!$    do iatom=1,npsfxatm
!!$       call read(psfxatm(iatom))
!!$       !call display(psfxatm(iatom))
!!$    end do
!!$    !close excited state psf
!!$    call close_psf
!!$
!!$    call note('Looking for atom matches!')
!!$    !rewind the pdb
!!$    call rewind_pdb
!!$    iatom=0
!!$    nmatch=0
!!$    mask=.false.
!!$    atomto=0
!!$    atomfrom=0
!!$    atomcenter=0
!!$    do while(trim(Next_Record()).NE.'END')
!!$       if(trim(Next_Record()).EQ.'ATOM')then
!!$          match_res=.False.
!!$          match_resname=.False.
!!$          match_segname=.False.
!!$        
!!$          call read(atm)
!!$          iatom=iatom+1
!!$        
!!$          if(adjustl(atm%resname).EQ.adjustl(atm_select%resname))match_resname=.True.
!!$          if(adjustl(atm%resSeq).EQ.adjustl(atm_select%resSeq))match_res=.True.
!!$          if(adjustl(atm%segName).EQ.adjustl(atm_select%segName))match_segname=.True.
!!$                
!!$          match=match_resname.and.match_res.and.match_segname
!!$        
!!$          if(match)then
!!$             nmatch=nmatch+1
!!$             mask(iatom)=.true.
!!$             !call display(atm) !< uncomment if you want to see the matched atoms
!!$           
!!$             if(adjustl(atm%name).EQ.adjustl(dipole_to%name))then
!!$                atomto=iatom
!!$                call display(atm)
!!$             end if
!!$             if(adjustl(atm%name).EQ.adjustl(dipole_from%name))then
!!$                atomfrom=iatom
!!$                call display(atm)
!!$             end if
!!$             if(adjustl(atm%name).EQ.adjustl(dipole_center%name))then
!!$                atomcenter=iatom
!!$                call display(atm)
!!$             end if
!!$
!!$             !assign excited state charges if available
!!$             do jatom=1,npsfxatm
!!$                match_name=.false.
!!$                if(adjustl(psfxatm(jatom)%name).EQ.adjustl(atm%name))&
!!$                     Qx(iatom)=psfxatm(jatom)%charge
!!$             end do
!!$          end if
!!$       else
!!$          call advance_record
!!$       end if
!!$    end do
!!$    call note('atom matches= '//trim(int2str(nmatch)))
!!$    if(atomto==0)then
!!$       call warn('ABPDB: Cannot find TO atom.','exiting')
!!$       pass=.false.
!!$    end if
!!$    if(atomfrom==0)then
!!$       call warn('cannot find FROM atom','exiting')
!!$       pass=.false.
!!$    end if
!!$    if(atomcenter==0)then
!!$       call warn('cannot find center atom','exiting')
!!$       pass=.false.
!!$    end if
!!$    call close_pdb
!!$
!!$
!!$    if(pass)then
!!$
!!$       call note('Computing Absorption Spectrum.')
!!$
!!$       !setup energy domain and allocate spectrum
!!$       Emax=this%Emax
!!$       Emin=this%Emin
!!$       dw=(Emax-Emin)/real(this%npt)
!!$       if(allocated(spectrum))deallocate(spectrum)
!!$       allocate(spectrum(0:this%npt-1))
!!$
!!$       !allocate EigenVectors and EigenValues
!!$       if(allocated(H))deallocate(H)
!!$       allocate(H(nstate,nstate))
!!$       if(allocated(Eng))deallocate(Eng)
!!$       allocate(Eng(nstate))
!!$
!!$       !allocate and compute Mu dot Mu 
!!$       if(allocated(MM))deallocate(MM)
!!$       allocate(MM(nstate,nstate))
!!$       
!!$
!!$       !sum over realizations
!!$       call open_pdb(pdb)
!!$       this%spectrum=0.0_double
!!$       spectrum=0.0_double
!!$       nframe=0
!!$       do while(len(trim(Next_Record())).NE.0)
!!$          iatom=0
!!$
!!$          !get current frame coordinates
!!$          call note('Processing Frame: '//trim(int2str(nframe)))
!!$          do while(trim(Next_Record()).NE.'END')
!!$             if(trim(Next_Record()).EQ.'ATOM')then
!!$                call read(atm)
!!$                iatom=iatom+1
!!$                XYZ(1,iatom)=atm%X
!!$                XYZ(2,iatom)=atm%Y
!!$                XYZ(3,iatom)=atm%Z   
!!$             else
!!$                call advance_record
!!$             end if
!!$          end do
!!$          if(trim(Next_Record()).EQ.'END')then
!!$             call advance_record
!!$             nframe=nframe+1
!!$             if(iatom.NE.natom)&
!!$                  call stop('ABPDB: number of atoms in frame '&
!!$                  //trim(int2str(nframe))//' not equal to natoms!','exiting.')
!!$          end if
!!$          if(iatom.NE.0)then
!!$             
!!$             !compute groundstate dipole moment magnitude
!!$             Dg=0._double
!!$             do jatom=1,natom!loop over atoms j (pigment mask)
!!$                if(mask(katom))then
!!$                   Dg=Dg+(Qg(katom)*(XYZ(:,atomcenter)-XYZ(:,katom)))**2
!!$                end if
!!$             end do
!!$
!!$left off here
!!$             !compute excitedstate dipole moment magnitude
!!$             Dg=0._double
!!$             do jatom=1,natom!loop over atoms j (pigment mask)
!!$                if(mask(katom))then
!!$                   Dg=Dg+(Qg(katom)*(XYZ(:,atomcenter)-XYZ(:,katom)))**2
!!$                end if
!!$             end do
!!$
!!$compute diference dmag=dx-dg
!!$
!!$             !accumulate sum delE
!!$             delE=0._double
!!$             !compute dipole vector
!!$             dm=XYZ(:,atomto)-XYZ(:,atomfrom)
!!$             dm=dm/sqrt(sum(dm**2))*Dmag
!!$             
!!$             !N loop over atoms k (not pigment mask)
!!$             do katom=1,natom
!!$                if(.not.mask(katom))then
!!$                   
!!$                   !compute dist vector
!!$                   R=XYZ(:,katom)-XYZ(:,atomcenter)
!!$                   R3=sqrt(sum(R**2))
!!$                   R3=R3*R3*R3
!!$                   if(R3.EQ.0)call stop('dipole center'//trim(int2str(atomcenter))&
!!$                        //' and env atom '//trim(int2str(katom))&
!!$                        //' sit on top of each other, distance = 0') 
!!$                   delE=delE+Qg(katom)*sum(dm*R)/R3
!!$                end if
!!$             end do
!!$             delE=delE*kc/eff
!!$             inquire(unit,opened=match)
!!$             if(.not.match)open(unit,file=trim(outfile),access='append')
!!$             write(unit,*)nframe,delE
!!$          end if
!!$
!!$
!!$          call resample(qs)
!!$
!!$
!!$             
!!$             H=cmplx(qs%hs%diabat)
!!$             if(present(cp))H=H+cmplx(cp%hc%V)
!!$             call diagonalize(nstate,H,Eng)
!!$             Eng=Eng-qs%hs%Eg
!!$
!!$
!!$             MM=0.0_double
!!$             do istate=1,nstate
!!$                do jstate=istate,nstate
!!$                   do idim=1,ndim
!!$                      MM(istate,jstate)=MM(istate,jstate)+mu(istate,idim)*mu(jstate,idim)
!!$                   end do
!!$                   MM(jstate,istate)=MM(istate,jstate)
!!$                end do
!!$             end do
!!$
!!$
!!$             !sum over adiabats
!!$             do kstate=1,nstate
!!$                !sum over diabats
!!$                tr=0.0_double
!!$                do istate=1,nstate
!!$                   do jstate=1,nstate
!!$                      tr=tr+real(H(istate,kstate)*conjg(H(jstate,kstate)))&
!!$                           *MM(istate,jstate)
!!$                   end do
!!$                end do                
!!$                tr=Eng(kstate)*tr
!!$
!!$
!!$                !build up spectrum 
!!$                do i=0,this%npt-1
!!$                   w=Emin+i*dw
!!$                   if(Eng(kstate).GT.w-dw/2_double.and.Eng(kstate).LE.w+dw/2_double)spectrum(i)=spectrum(i)+tr
!!$                end do
!!$             end do
!!$          end do
!!$
!!$          !decide to keep or add more realizations 
!!$          norm=sum(spectrum)
!!$          if(norm.EQ.0.)then
!!$             call warn('get_ABspectrum: spectrum has zero valued integral.'&
!!$                  ,'exiting without computing spectrum.')
!!$             stable=.false.
!!$          end if
!!$          R2=0.0_double
!!$          if(stable)R2=sum((spectrum/real(nsample)-this%spectrum)**2)/real(this%npt)
!!$          if(R2.LE.tolerance)cnvrg=.true.
!!$          call Note('ABspectrum: samples='//trim(int2str(nsample)&
!!$               //' R='//trim(float2str(sqrt(R2)))))
!!$
!!$          if(nsample.GT.killint)then
!!$             call warn('get_ABspectrum: spectrum failed to converge after '&
!!$                  //trim(int2str(killint))//' realizations.'&
!!$                  ,'Stop trying to converge and continue as is.')
!!$             cnvrg=.true.
!!$          end if
!!$          
!!$          !save spectrum as old spectrum and double number of samples
!!$          this%spectrum=spectrum/real(nsample)
!!$          lastsample=nsample
!!$          nsample=int(nsample*factor)
!!$
!!$          if(present(file))then
!!$             call Note('Saving intermediate ABspec')
!!$             unit=newunit()
!!$             open(unit,file=file)
!!$             norm=maxval(abs(this%spectrum))
!!$             do i=0,this%npt-1
!!$                w=Emin+i*dw
!!$                write(unit,"(10(ES18.10E2,1X))")w,this%spectrum(i)/norm
!!$             end do
!!$             close(unit)
!!$          end if
!!$
!!$       end do
!!$       nsample=lastsample
!!$       call Note('ABspectrum: samples='//trim(int2str(nsample)&
!!$            //' R_squared ='//trim(float2str(R2))))
!!$       
!!$       if(.not.stable)spectrum=0.0_double !incase of runtime error
!!$       
!!$       if(present(file))then
!!$          call note('input file= '//file)
!!$
!!$          norm=maxval(abs(this%spectrum))
!!$
!!$          unit=newunit()
!!$          open(unit,file=file)
!!$          do i=0,this%npt-1
!!$             w=Emin+i*dw
!!$             write(unit,"(10(ES18.10E2,1X))")w,this%spectrum(i)/norm
!!$          end do
!!$          close(unit)
!!$
!!$       end if
!!$    else
!!$       call warn('get_ABspectrum: failed check.','not computing absorption spectrum.')
!!$    end if
!!$
!!$
!!$
!!$
!!$
!!$       end do
!!$  close(unit)
!!$

  !deallocate arrays
  if(allocated(mask))deallocate(mask)
  if(allocated(XYZ))deallocate(XYZ)
  if(allocated(Qg))deallocate(Qg)



    
    if(allocated(H))deallocate(H)
    if(allocated(Eng))deallocate(Eng)
    if(allocated(MM))deallocate(MM)
    if(allocated(spectrum))deallocate(spectrum)

    if(allocated(Mu))deallocate(Mu)

  end subroutine get_ABPDBspectrum
!--------------------------------------------------------
  
end module Spectrometer_class
   
