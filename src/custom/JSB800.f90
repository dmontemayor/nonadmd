!.................................................................
!LH2 Dimer Model
!------------------------------
!Daniel Montemayor
!Feb 2011
!==============================
!Updates:
! ? No
! (+) Convert to subroutine that builds the dipole object
!------------v1.2 Aug 2011
!(*)Bias changed from 260 to 280 /cm to match monomer AB lineshape
!   [JPC B, 2007, 111, 6807-6814] 
!------------v1.1 March 2011
!(*)Added groundstate fluctuations
!(*)Added Kruegger cutoff distances
!(*)Expects optical dielectric matrix 
!(*)Intercomplex distance is now optional parameter
!(*)Added hs%Eg groundstate reference energy. [-850nm excitation]
!   [B850-BChl ave site energy is 0 so ground state is 850nm lower]
!(*)Gvar is 0 and Evar is written in terms of sigEn and sigED [absolute energies]
!   not sigEng and sigEDg [energy gap between excited state and groundstate]
!   assumes states have ground states of equall energy but not same ground state itself. 
!------------v1.0 Feb 2011
!==============================
!TODO:
!()put in dimer tilt parameter
!==============================
!Observations:
!(+)B800 ABspec lineshape is a bit broader than in the Volker spectrum. B800 bath friction may be too strong.
!
!(+)Is the JS B800-B850 optical dielectric due to breakdown of resonance dipole approx.
!==============================
!KEY:
!(+)Issue not yet addressed
!(-)Issue not valid or cannot be considered
!(?)Issue under investigation.
!(*)Issue addressed. No follow up needed.
!(&)Terminal Solution/Validation.
!(&+)Intermediate Solution/Validation. Follow up expected but not yet addressed.
!(&?)Intermediate Solution/Validation. Follow up under investigation. 
!(&*)Intermediate Solution/Validation. Follow up concluded. 
!(&-)Solution cannot be applied. 
!.................................................................

 
subroutine JSB800(file)
  use type_kinds
  use math
  use atomicunits
  use rand
  use ErrorLog
  use String
  use dipoles_class
  use filemanager
  implicit none


  type(dipoles)::qs
  character*(*),intent(in)::file
  integer(long),parameter::Rstate=9,shift=0
  integer(long)::map(Rstate)
  !                                                            !Krueger !JS
  real(double),parameter::trans          =78.0_double*angstrom !NA(78)  !NA(78)
  real(double),parameter::B800B800_cutoff=12.0_double*angstrom !32      !12
  real(double),parameter::B800B850_cutoff=12.0_double*angstrom !22      !12
  real(double),parameter::B850B850_cutoff=18.0_double*angstrom !18      !18

  integer(short)::unit,ierr
  logical::usedunit

  !-----------quantum subsystem--------
  integer(short),parameter::nD=9,nstate=6*nD,ndim=3
  real(double),parameter::nu=10.3_double*Deg,nuD=23.55_double*Deg

  real(double),parameter::muD=6.13_double*Debye
  !real(double),parameter::muD=sqrt(1.895_double*1E5*invcm*(angstrom**3)/kC)
  real(double),parameter::muA=MuD
  real(double),parameter::muB=MuD

  real(double),parameter::thetaD=90.0_double*Deg
  !real(double),parameter::thetaD=5.1_double*Deg

  !real(double),parameter::thetaA=5.1_double*Deg 
  real(double),parameter::thetaA=84.9_double*Deg !(from ref)

  real(double),parameter::thetaB=thetaA          !(from ref)
  real(double),parameter::phiD=244.57_double*Deg
  real(double),parameter::phiA=-112.5_double*Deg
  real(double),parameter::phiB=63.2_double*Deg

  real(double),parameter::RD=31.0_double*angstrom,ZD=16.6_double*angstrom
  real(double),parameter::RA=26.0_double*angstrom,RB=27.2_double*angstrom
 
  !real(double),parameter::Eg=-invcm*1E7/802.0_double !with counterterm
  !real(double),parameter::Eg=-invcm*1E7/812.0_double !without counterterm
  real(double),parameter::Eg=0.0
 
  !real(double),parameter::En=0.0_double*invcm
  real(double),parameter::En=invcm*1E7/812.0_double

  real(double),parameter::bias=260.0_double*invcm

  real(double),parameter::sigEg=40.0_double*invcm    !external site disorder

  real(double),parameter::sigEDg=54.0_double*invcm   !B800 internal site disorder (from ref)
  real(double),parameter::sigEng=246.78_double*invcm !B850 internal site disorder [sqrt(sigEn**2-sigEg**2)]

!!$  real(double),parameter::sigEn=sqrt(sigEng**2+sigEg**2) !250 invcm (from ref) =Total B850 site energy fluct
!!$  real(double),parameter::sigED=sqrt(sigEDg**2+sigEg**2) !from convolution  =Total B800 site energy fluct
  real(double)::sigEn
  real(double)::sigED

  !default values for optional parameters
  real(double)::intercmplxdist=76_double*angstrom
  real(double)::optdielect(nstate,nstate)=1.0_double

  real(double)::muhat(nstate,ndim),Q(nstate,ndim)
  real(double)::diabat(nstate,nstate)
  real(double)::Evar(nstate,nstate)
  integer(long)::istate,idim,n,jstate,m,index
  real(double)::phi,R,phi0,nu0,theta0
  integer(long)::tmap(nstate)
  logical::IDA

  !debug vars
  integer(long),parameter::nmcstep=1E7,nsamp=1E3,ntrainset=14
  integer(long)::i,j,mcstep,naccept,ncount,trainset
  integer(long)::ni,nj,ki,kj
  character(len=5),parameter::Head='(A51)'
  character(len=42),parameter::Body='(A8,A8,F8.3,4X,F6.3,1X,F8.2,1X,F9.2)'
  character(len=8)::donor,acceptor
  real(double),dimension(ntrainset)::KExp,Kt,Qt,VExp,Vt
  real(double)::kappa,kvart,kvar0,kvarmin,beta
  real(double)::V,Vvart,Vvar0,Vvarmin
  real(double),dimension(3)::rhat1,rhat2,Rhat
  real(double)::randtheta(3),dtheta
  real(double)::thetaDt,thetaAt,thetaBt
  real(double)::thetaDt0,thetaAt0,thetaBt0
  real(double)::thetaDtmin,thetaAtmin,thetaBtmin
  real(double)::randphi(3),dphi
  real(double)::phiDt,phiAt,phiBt
  real(double)::phiDt0,phiAt0,phiBt0
  real(double)::phiDtmin,phiAtmin,phiBtmin
  real(double)::randmu(3),dmu,muistate,mujstate
  real(double)::muDt,muAt,muBt
  real(double)::muDt0,muAt0,muBt0
  real(double)::muDtmin,muAtmin,muBtmin

  call note('Creating JSB800.')
  sigEn=sqrt(sigEng**2+sigEg**2) !250 invcm (from ref) =Total B850 site energy fluct
  sigED=sqrt(sigEDg**2+sigEg**2) !from convolution  =Total B800 site energy fluct

  !optional parameters
  intercmplxdist=trans
  tmap=0
  do istate=1,Rstate
     Map(istate)=istate+shift
     tmap(Map(istate))=istate
  end do

  !Krueger orientation factors
  Kexp(1)=-1.43_double
  Kexp(2)=1.11_double
  Kexp(3)=1.68_double
  Kexp(4)=-1.26_double
  Kexp(5)=-1.18_double
  Kexp(6)=1.05_double
  Kexp(7)=-0.87_double
  Kexp(8)=0.13_double
  Kexp(9)=0.79_double
  Kexp(10)=0.17_double
  Kexp(11)=-0.55_double
  Kexp(12)=1.30_double
  Kexp(13)=-1.33_double
  Kexp(14)=-1.04_double

  !Krueger Coulombic potentials
  Vexp(1)=-46._double
  Vexp(2)=213._double
  Vexp(3)=238._double
  Vexp(4)=-37._double
  Vexp(5)=-4._double
  Vexp(6)=7._double
  Vexp(7)=-13._double
  Vexp(8)=5._double
  Vexp(9)=27._double
  Vexp(10)=23._double
  Vexp(11)=-4._double
  Vexp(12)=5._double
  Vexp(13)=-26._double
  Vexp(14)=-3._double

  !begin montecarlo to best fit orientation factors
  thetaDt=thetaD
  thetaAt=thetaA
  thetaBt=thetaB
  dtheta=pi*1E-4
  phiDt=phiD
  phiAt=phiA
  phiBt=phiB
  dphi=twopi*1E-4
  muDt=muD
  muAt=muA
  muBt=muB
  dmu=1E-4
  beta=1.0_double
  naccept=0
  ncount=0

  !------make quantum subsystem primitive
  diabat=0.0
  !B800 energy
  do istate=1,nD
     diabat(istate,istate)=bias+En
  end do
  !B850 energy
  do istate=nD+1,3*nD
     diabat(istate,istate)=En
  end do

  !......monomer 1 with monomer 1 optical dielectric.......
  !B800 with B800
  do istate=1,nD
     do jstate=1,nD
        optdielect(istate,jstate)=1.0_double
     end do
  end do
  !B800 with B850
  do istate=1,nD
     do jstate=nD+1,3*nD
        optdielect(istate,jstate)=(MuD/(5.3_double*Debye))**2!1.343
     end do
     !write(*,*)optdielect(istate,:)
  end do
  !B850 with B850
  do istate=nD+1,3*nD
     do jstate=nD+1,3*nD
        optdielect(istate,jstate)=1.0_double
     end do
  end do
  !......monomer 1 with monomer 2 optical dielectric.......
  !B800 with B800
  do istate=1,nD
     do jstate=nD*3+1,nD*4
        optdielect(istate,jstate)=1E3
     end do
  end do
  !B800 with B850
  do istate=1,nD
     do jstate=nD*4+1,nstate
        optdielect(istate,jstate)=1E3
     end do
  end do
  !B850 with B800
  do istate=nD+1,3*nD
     do jstate=3*nD+1,4*nD
        optdielect(istate,jstate)=1E3
     end do
  end do
  !B850 with B850
  do istate=nD+1,3*nD
     do jstate=4*nD+1,nstate
        optdielect(istate,jstate)=1E3
     end do
  end do

  !copy values to ring 2
  do istate=1,nstate/2
     do jstate=istate,nstate/2
        diabat(istate+nstate/2,jstate+nstate/2)=diabat(istate,jstate)
        optdielect(istate+nstate/2,jstate+nstate/2)=optdielect(istate,jstate)
     end do
  end do
  !Make symmetric
  do istate=1,nstate
     do jstate=istate,nstate
        diabat(jstate,istate)=diabat(istate,jstate)
        optdielect(jstate,istate)=optdielect(istate,jstate)
     end do
  end do

  unit=newunit()
  open(unit,file=file//'.dipoles.hs')
  !hs
  write(unit,*)'hs'
  write(unit,*)Rstate
  write(unit,*)Eg
  write(unit,*)((diabat(Map(istate),Map(jstate))&
       ,jstate=1,Rstate),istate=1,Rstate)
  write(unit,*)1.0_double,(0.0_double,istate=2,Rstate*Rstate)
  close(unit)

  !------make dipole quantum subsystem
  do istate=1,nD
     phi0=twopi*(istate-1)/real(nD)
     phi=phiDt+nuD
     muhat(istate,1)=sin(thetaDt)*cos(phi0+phi)
     muhat(istate,2)=sin(thetaDt)*sin(phi0+phi)
     muhat(istate,3)=cos(thetaDt)
     !Note: muhat not a unit vector anymore
     muhat(istate,:)=muDt*muhat(istate,:)
     Q(istate,1)=RD*cos(phi0+nuD)
     Q(istate,2)=RD*sin(phi0+nuD)
     Q(istate,3)=ZD
  end do
  n=0
  do istate=nD+1,nstate/2
     if(mod(istate,2).EQ.0)then
        R=RA
        phi=phiAt-nu
        nu0=-nu
        muistate=muAt
        theta0=thetaAt
     else
        R=RB
        phi=phiBt+nu
        nu0=nu
        muistate=muBt
        theta0=thetaBt
     end if
     phi0=twopi*n/real(nD)
     muhat(istate,1)=sin(theta0)*cos(phi0+phi)
     muhat(istate,2)=sin(theta0)*sin(phi0+phi)
     muhat(istate,3)=cos(theta0)
     Q(istate,1)=R*cos(phi0+nu0)
     Q(istate,2)=R*sin(phi0+nu0)
     Q(istate,3)=0.0_double
     n=n+mod(istate,2)
     !Note: muhat not a unit vector anymore
     muhat(istate,:)=muistate*muhat(istate,:)
  end do
  !translate transition dipoles, coordinates and energies onto second LHII
  do istate=1,nstate/2
     muhat(istate+nstate/2,:)=muhat(istate,:)
     Q(istate+nstate/2,:)=Q(istate,:)
     Q(istate+nstate/2,1)=Q(istate+nstate/2,1)+intercmplxdist
  end do
  unit=newunit()
  open(unit,file=file//'.dipoles')
  write(unit,*)'dipoles'
  write(unit,*)quote(file//'.dipoles.hs')
  write(unit,*)0.0_double
  write(unit,*)((optdielect(Map(istate),Map(jstate)),jstate=1,Rstate)&
       ,istate=1,Rstate)
  write(unit,*)ndim
  write(unit,*)((muhat(Map(istate),idim),idim=1,ndim),istate=1,Rstate)
  write(unit,*)((Q(Map(istate),idim),idim=1,ndim),istate=1,Rstate)
  write(unit,*)((0.0_double,jstate=1,Rstate),istate=1,Rstate)
  write(unit,*)0.0_double
  write(unit,*)((0.0_double,jstate=1,Rstate),istate=1,Rstate)
  close(unit)

  !create disorderless qs and update
  call new(qs,file=file//'.dipoles')
  qs%Evar=0._double
  qs%Gvar=0._double
  call resample(qs)

  !make corrections to hs according to Krueger and save
  do istate=1,nD
     do jstate=1,nD
        !construct ring1
        if(tmap(istate).NE.0.and.tmap(jstate).NE.0)then
           Rhat=qs%Q(tmap(istate),:)-qs%Q(tmap(jstate),:)
           R=sqrt(sum(Rhat**2))
           if(R.LT.B800B800_cutoff)then
              index=0

              m=1 !sites D(n)&D(n+1)
              if(mod(istate+m,9).EQ.mod(jstate,9))index=13

              m=2 !sites D(n)&D(n+2)
              if(mod(istate+m,9).EQ.mod(jstate,9))index=14

              if(index.NE.0)&
                   qs%hs0%diabat(tmap(istate),tmap(jstate))=Vexp(index)*invcm&
                   -qs%hs%diabat(tmap(istate),tmap(jstate))

              !symmeterize ring 1
              qs%hs0%diabat(tmap(jstate),tmap(istate))=&
                   qs%hs0%diabat(tmap(istate),tmap(jstate))

              !copy ring1 to ring2
              if(tmap(istate+nstate/2).NE.0.and.tmap(jstate+nstate/2).NE.0)then
                 qs%hs0%diabat(tmap(istate+nstate/2),tmap(jstate+nstate/2))=&
                      qs%hs0%diabat(tmap(istate),tmap(jstate))
                 qs%hs0%diabat(tmap(jstate+nstate/2),tmap(istate+nstate/2))=&
                      qs%hs0%diabat(tmap(istate),tmap(jstate))
              end if
           end if
        end if
     end do
  end do

  do istate=1,nD
     do jstate=nD+1,3*nD
        !construct ring1
        if(tmap(istate).NE.0.and.tmap(jstate).NE.0)then
           Rhat=qs%Q(tmap(istate),:)-qs%Q(tmap(jstate),:)
           R=sqrt(sum(Rhat**2))
           if(R.LT.B800B850_cutoff)then
              index=0

              m=1 !sites D(n+1)&A(n)
              if(mod(istate,9).EQ.mod((jstate-nD-1)/2+1+m,nD).and.mod(jstate,2).EQ.0)&
                   index=5

              m=1 !sites D(n+1)&B(n)
              if(mod(istate,9).EQ.mod((jstate-nD-1)/2+1+m,nD).and.mod(jstate,2).EQ.1)&
                   index=6
              
              m=0 !sites D(n)&A(n)
              if(mod(istate,9).EQ.mod((jstate-nD-1)/2+1+m,nD).and.mod(jstate,2).EQ.0)&
                   index=7

              m=0 !sites D(n)&B(n)
              if(mod(istate,9).EQ.mod((jstate-nD-1)/2+1+m,nD).and.mod(jstate,2).EQ.1)&
                   index=8
              
              m=1 !sites D(n)&A(n+1)
              if(mod(istate+m,9).EQ.mod((jstate-nD-1)/2+1,nD).and.mod(jstate,2).EQ.0)&
                   index=9
              
              m=1 !sites D(n)&B(n+1)
              if(mod(istate,9).EQ.mod((jstate-nD-1)/2+1+m,nD).and.mod(jstate,2).EQ.1)&
                   index=10

              m=2 !sites D(n)&A(n+2)
              if(mod(istate+m,9).EQ.mod((jstate-nD-1)/2+1,nD).and.mod(jstate,2).EQ.0)&
                   index=11
              
              m=2 !sites D(n)&B(n+2)
              if(mod(istate,9).EQ.mod((jstate-nD-1)/2+1+m,nD).and.mod(jstate,2).EQ.1)&
                   index=12
              

              if(index.NE.0)&
                   qs%hs0%diabat(tmap(istate),tmap(jstate))=Vexp(index)*invcm&
                   -qs%hs%diabat(tmap(istate),tmap(jstate))

              !symmeterize ring 1
              qs%hs0%diabat(tmap(jstate),tmap(istate))=&
                   qs%hs0%diabat(tmap(istate),tmap(jstate))
              !copy ring1 to ring2
              if(tmap(istate+nstate/2).NE.0.and.tmap(jstate+nstate/2).NE.0)then
                 qs%hs0%diabat(tmap(istate+nstate/2),tmap(jstate+nstate/2))=&
                      qs%hs0%diabat(tmap(istate),tmap(jstate))
                 qs%hs0%diabat(tmap(jstate+nstate/2),tmap(istate+nstate/2))=&
                      qs%hs0%diabat(tmap(istate),tmap(jstate))
              end if
           end if
        end if
     end do
  end do

  do istate=nD+1,3*nD
     do jstate=nD+1,3*nD
        !construct ring1
        if(tmap(istate).NE.0.and.tmap(jstate).NE.0)then
           Rhat=qs%Q(tmap(istate),:)-qs%Q(tmap(jstate),:)
           R=sqrt(sum(Rhat**2))
           if(R.LT.B850B850_cutoff)then
              index=0

              m=1 !sites A(n+1)&A(n)
              if(mod((istate-nD-1)/2+1,nD).EQ.mod((jstate-nD-1)/2+1+m,nD)&
                   .and.mod(istate,2).EQ.0.and.mod(jstate,2).EQ.0)&
                   index=1

              m=1 !sites A(n+1)&B(n)
              if(mod((istate-nD-1)/2+1,nD).EQ.mod((jstate-nD-1)/2+1+m,nD)&
                   .and.mod(istate,2).EQ.0.and.mod(jstate,2).EQ.1)&
                   index=2

              m=0 !sites A(n)&B(n)
              if(mod((istate-nD-1)/2+1,nD).EQ.mod((jstate-nD-1)/2+1+m,nD)&
                   .and.mod(istate,2).EQ.0.and.mod(jstate,2).EQ.1)&
                   index=3
                         
              m=1 !sites B(n+1)&B(n)
              if(mod((istate-nD-1)/2+1,nD).EQ.mod((jstate-nD-1)/2+1+m,nD)&
                   .and.mod(istate,2).EQ.1.and.mod(jstate,2).EQ.1)&
                   index=4
              
              if(index.NE.0)&
                   qs%hs0%diabat(tmap(istate),tmap(jstate))=Vexp(index)*invcm&
                   -qs%hs%diabat(tmap(istate),tmap(jstate))

              !symmeterize ring 1
              qs%hs0%diabat(tmap(jstate),tmap(istate))=&
                   qs%hs0%diabat(tmap(istate),tmap(jstate))
              !copy ring1 to ring2
              if(tmap(istate+nstate/2).NE.0.and.tmap(jstate+nstate/2).NE.0)then
                 qs%hs0%diabat(tmap(istate+nstate/2),tmap(jstate+nstate/2))=&
                      qs%hs0%diabat(tmap(istate),tmap(jstate))
                 qs%hs0%diabat(tmap(jstate+nstate/2),tmap(istate+nstate/2))=&
                      qs%hs0%diabat(tmap(istate),tmap(jstate))
              end if
           end if
        end if
     end do
  end do

  !add disorder to qs, save and update
  !qs%Gvar=0.0_double !when using total disorder
  qs%Gvar=sigEg*sigEg !external disorder
  Evar=0.0_double
  do istate=1,nD
     !Evar(istate,istate)=sigED*sigED  !total disorder
     Evar(istate,istate)=sigEDg*sigEDg !internal disorder
  end do
  do istate=nD+1,3*nD
     !Evar(istate,istate)=sigEn*sigEn  !total disorder
     Evar(istate,istate)=sigEng*sigEng !internal disorder
  end do
  do istate=1,nstate/2
     Evar(istate+nstate/2,istate+nstate/2)=Evar(istate,istate)
  end do
  do istate=1,qs%hs%nstate
     do jstate=1,qs%hs%nstate
        qs%Evar(istate,jstate)=Evar(Map(istate),Map(jstate))
     end do
  end do
  call resample(qs)
  call save(qs,file=file//'.dipoles')

  !write quantum object
  unit=newunit()
  open(unit,file=file)
  write(unit,*)'quantum'
  write(unit,*)'dipoles'
  write(unit,*)quote(file//'.dipoles')
  close(unit)

  qs%Evar=0._double
  qs%Gvar=0._double
  call update(qs)

  unit=newunit()
  open(unit,file=file//'.log')
  write(unit,'(6X,3(A12,1X))')'B800','aB850','bB850'
  write(unit,'(A5,1X,3(ES12.4,1X))')'Mu'&
       ,muDt/debye,muAt/debye,muBt/debye
  write(unit,'(A5,1X,3(ES12.4,1X))')'Theta'&
       ,thetaDt/Deg,thetaAt/Deg,thetaBt/Deg
  write(unit,'(A5,1X,3(ES12.4,1X))')'Phi'&
       ,phiDt/Deg,phiAt/Deg,phiBt/Deg
  write(unit,*)
  write(unit,Head)'Donor Acceptor Separation Orientation   V    V(corrected)'
  write(unit,Head)'Trans  Trans   (Angstrom)   Factor    (1/cm) (1/cm)'
  Vvarmin=0.0
  do istate=1,Rstate-1
     select case(Map(istate))
     case(1)
        donor='LB800_1'
     case(2)
        donor='LB800_2'
     case(3)
        donor='LB800_3'
     case(4)
        donor='LB800_4'
     case(5)
        donor='LB800_5'
     case(6)
        donor='LB800_6'
     case(7)
        donor='LB800_7'
     case(8)
        donor='LB800_8'
     case(9)
        donor='LB800_9'
     case(10)
        donor='LaB850_1'
     case(11)
        donor='LbB850_1'
     case(12)
        donor='LaB850_2'
     case(13)
        donor='LbB850_2'
     case(14)
        donor='LaB850_3'
     case(15)
        donor='LbB850_3'
     case(16)
        donor='LaB850_4'
     case(17)
        donor='LbB850_4'
     case(18)
        donor='LaB850_5'
     case(19)
        donor='LbB850_5'
     case(20)
        donor='LaB850_6'
     case(21)
        donor='LbB850_6'
     case(22)
        donor='LaB850_7'
     case(23)
        donor='LbB850_7'
     case(24)
        donor='LaB850_8'
     case(25)
        donor='LbB850_8'
     case(26)
        donor='LaB850_9'
     case(27)
        donor='LbB850_9'
     case(28)
        donor='RB800_1'
     case(29)
        donor='RB800_2'
     case(30)
        donor='RB800_3'
     case(31)
        donor='RB800_4'
     case(32)
        donor='RB800_5'
     case(33)
        donor='RB800_6'
     case(34)
        donor='RB800_7'
     case(35)
        donor='RB800_8'
     case(36)
        donor='RB800_9'
     case(37)
        donor='RaB850_1'
     case(38)
        donor='RbB850_1'
     case(39)
        donor='RaB850_2'
     case(40)
        donor='RbB850_2'
     case(41)
        donor='RaB850_3'
     case(42)
        donor='RbB850_3'
     case(43)
        donor='RaB850_4'
     case(44)
        donor='RbB850_4'
     case(45)
        donor='RaB850_5'
     case(46)
        donor='RbB850_5'
     case(47)
        donor='RaB850_6'
     case(48)
        donor='RbB850_6'
     case(49)
        donor='RaB850_7'
     case(50)
        donor='RbB850_7'
     case(51)
        donor='RaB850_8'
     case(52)
        donor='RbB850_8'
     case(53)
        donor='RaB850_9'
     case(54)
        donor='RbB850_9'
     end select


     rhat1=qs%mu(istate,:)/sqrt(sum(qs%mu(istate,:)**2))

     do jstate=istate+1,Rstate
        select case(Map(jstate))
        case(1)
           acceptor='LB800_1'
        case(2)
           acceptor='LB800_2'
        case(3)
           acceptor='LB800_3'
        case(4)
           acceptor='LB800_4'
        case(5)
           acceptor='LB800_5'
        case(6)
           acceptor='LB800_6'
        case(7)
           acceptor='LB800_7'
        case(8)
           acceptor='LB800_8'
        case(9)
           acceptor='LB800_9'
        case(10)
           acceptor='LaB850_1'
        case(11)
           acceptor='LbB850_1'
        case(12)
           acceptor='LaB850_2'
        case(13)
           acceptor='LbB850_2'
        case(14)
           acceptor='LaB850_3'
        case(15)
           acceptor='LbB850_3'
        case(16)
           acceptor='LaB850_4'
        case(17)
           acceptor='LbB850_4'
        case(18)
           acceptor='LaB850_5'
        case(19)
           acceptor='LbB850_5'
        case(20)
           acceptor='LaB850_6'
        case(21)
           acceptor='LbB850_6'
        case(22)
           acceptor='LaB850_7'
        case(23)
           acceptor='LbB850_7'
        case(24)
           acceptor='LaB850_8'
        case(25)
           acceptor='LbB850_8'
        case(26)
           acceptor='LaB850_9'
        case(27)
           acceptor='LbB850_9'
        case(28)
           acceptor='RB800_1'
        case(29)
           acceptor='RB800_2'
        case(30)
           acceptor='RB800_3'
        case(31)
           acceptor='RB800_4'
        case(32)
           acceptor='RB800_5'
        case(33)
           acceptor='RB800_6'
        case(34)
           acceptor='RB800_7'
        case(35)
           acceptor='RB800_8'
        case(36)
           acceptor='RB800_9'
        case(37)
           acceptor='RaB850_1'
        case(38)
           acceptor='RbB850_1'
        case(39)
           acceptor='RaB850_2'
        case(40)
           acceptor='RbB850_2'
        case(41)
           acceptor='RaB850_3'
        case(42)
           acceptor='RbB850_3'
        case(43)
           acceptor='RaB850_4'
        case(44)
           acceptor='RbB850_4'
        case(45)
           acceptor='RaB850_5'
        case(46)
           acceptor='RbB850_5'
        case(47)
           acceptor='RaB850_6'
        case(48)
           acceptor='RbB850_6'
        case(49)
           acceptor='RaB850_7'
        case(50)
           acceptor='RbB850_7'
        case(51)
           acceptor='RaB850_8'
        case(52)
           acceptor='RbB850_8'
        case(53)
           acceptor='RaB850_9'
        case(54)
           acceptor='RbB850_9'
        end select

        rhat2=qs%mu(jstate,:)&
             /sqrt(sum(qs%mu(jstate,:)**2))

        Rhat=qs%Q(istate,:)-qs%Q(jstate,:)
        R=sqrt(sum(Rhat**2))
        kappa=0.0_double
        V=0.0_double
        Rhat=Rhat/R
        kappa=sum(rhat1*rhat2)&
             -3.0_double*sum(rhat1*rhat)*sum(rhat2*rhat)
        V=kC*kappa*sqrt(sum(qs%mu(istate,:)**2))&
             *sqrt(sum(qs%mu(jstate,:)**2))/(qs%dielectric(istate,jstate)*R*R*R)
        write(unit,Body)donor,acceptor,R/angstrom,kappa&
             ,V/invcm,qs%hs%diabat(istate,jstate)/invcm


        Vvart=R
        if(Vvarmin.EQ.0.)then
           Vvarmin=Vvart
           vvar0=vvart
        end if


        if(donor(1:1).NE.acceptor(1:1))then
           if(Vvart.LE.Vvarmin)then
              Vvarmin=Vvart
              call Note('JSB800: min pair ('//donor//','//acceptor//') distance in angstrom ='//trim(float2str(R/angstrom)))
           end if
           if(Vvart.GE.Vvar0)then
              Vvar0=Vvart
              call Note('JSB800: max pair ('//donor//','//acceptor//') distance in angstrom ='//trim(float2str(R/angstrom)))
           end if
        end if

     end do
  end do

  write(unit,*)'==================================================='
  write(unit,*)'JSB800 coordinates'
  write(unit,*)'  site       X                Y                Z'
  do istate=1,Rstate
     write(unit,'(I5,3(F16.8,1X))')istate,qs%Q(istate,:)
  end do
  close(unit)

  !stop

end subroutine JSB800
