!-----------------------------------------------------------------------------
!> \brief
!! Compute LH3 Spectral Density using Renger's method.
!<----------------------------------------------------------------------------
!> \details
!! Implements Eq5 in JPC B 2012, 116, 1164-1171. 
!! Input: multiframe PDB, ground state psf, excited state psf.
!! Output: The Qy transition site energy shift 
!! Program will use time dep coord and charge information. 
!! The program assumes the pdb and ground state psf files all have the same atom count
!! and atom ordering. The excited state psf may contain default charges for the chormophore only,
!! written either as the charge difference between ground and excited state(default) or the absolute
!! excited state charge. Edit the 'usedifQ' parameter to choose between the two.
!! Edit the code defining atm_select. This tells the program which atoms make 
!! up the pigment for which we are computing the site energy gap.
!! The excited state psf may contian default charges for a prototypical chromophre. These
!! charges may be used for different yet molecularly equivalent chromophores found in the pdb and ground state psf.
!<----------------------------------------------------------------------------
!>\authors
!! Daniel Montemayor
!!
!!\date
!! 21 Nov 2013 
!<----------------------------------------------------------------------------
program main
  use type_kinds
  use atomicunits
  use math
  use rand
  use string
  use MPIframework
  use filemanager
  use ErrorLog
  use molreader
  use wallclock
  use outputdisplay

  !hamiltonian  
  use quantum_class
  use classical_class
  use coupling_class

  !spectrometer
  use spectrometer_class

  !propagators
  use PLDM_class  

  implicit none
# include "config.h"

  type(quantum),target::qs
  type(classical),target::cs
  type(coupling),target::cp
  type(spectrum)::spec
  type(PLDM)::prop

  integer(long)::ierr

  !parameters used for fitting
  real(double)::param(2)
  real(double),external::f

  !-----variable pertinent to this experiment-----
  !params
  integer(long),parameter::Cn=9,nstate=Cn*3
  integer,parameter::nfile=1,fpf=400

  !output
  integer(long)::unit
  character(len=path)::outfile
  character(len=3)::index

  !input
  character(len=path)::pdb,psf,psfx

  type(ATOM)::atm,atm_select(nstate)
  type(psfATOM)::psfatm

  type(psfATOM),allocatable::psfxatm(:)

  logical::match,match_res,match_chain,match_resname,match_segname,match_name,chargeaveraged,usedifQ
  integer(long)::nframe,nmatch,natom,maxatomid,iatom,jatom,katom,npsfatom,nchromatm

  logical,allocatable::mask(:,:)
  real(double),allocatable::XYZ(:,:),Qg(:),Qx(:),delE(:)

  real(double)::delQ,R

  !book keeping
  integer::ifile,istate
  !------------------------------------------

!!$# ifdef MPI
!!$  call MPI_INIT(ierr)
!!$  call MPI_COMM_rank(MPI_COMM_World,myid,ierr)
!!$  call MPI_COMM_size(MPI_COMM_World,nproc,ierr)
!!$# endif

  !Developers can run system tests here
  !call tests

  call openLog(level=1_short,file='./runtime.log')
  call display('display.out')
  !--------------------   Body of experiment   ----------------------
  usedifQ=.true.

  psf='../gnu/share/examples/molreader/lh3mem.psf'
  psfx='../gnu/share/examples/molreader/RengerB3LYP_diff.psf'

  !assign atm_select
  atm_select%resname='BCL'
  do istate=1,nstate
     atm_select(istate)%resSeq=mod(istate-1,Cn)+1
     atm_select(istate)%segName='BCL'//trim(int2str(1+(istate-1)/Cn))
  end do

  !open output file
  unit=newunit()
  outfile='Cn'//trim(int2str(Cn))//'.out'
  open(unit,file=trim(outfile))

  !load the ground state psf
  call open_psf(trim(psf))

  !skip to atoms section
  call advance_psf(section='atoms',nrecord=npsfatom)
  natom=npsfatom

  !allocate arrays
  call Prompt('Allocating dynamic arrays for '&
       //trim(int2str(npsfatom))//' atoms.')
  if(allocated(mask))deallocate(mask)
  allocate(mask(natom,nstate))
  if(allocated(XYZ))deallocate(XYZ)
  allocate(XYZ(3,natom))
  if(allocated(Qg))deallocate(Qg)
  allocate(Qg(natom))
  if(allocated(Qx))deallocate(Qx)
  allocate(Qx(natom))
  if(allocated(delE))deallocate(delE)
  allocate(delE(natom))

  !read psf atoms
  do iatom=1,npsfatom
     call read(psfatm)
     !call display(psfatm)
     Qg(iatom)=psfatm%charge
  end do
  !close ground state psf
  call close_psf

  !load the excited state psf
  call open_psf(trim(psfx))

  !skip to atoms section
  call advance_psf(section='atoms',nrecord=nchromatm)

  !allocate excited state charges
  call Prompt('Allocating excited state charge dynamic array.')
  if(allocated(psfxatm))deallocate(psfxatm)
  allocate(psfxatm(nchromatm))

  !read psfx atoms
  do iatom=1,nchromatm
     call read(psfxatm(iatom))
     !call display(psfxatm(iatom))
  end do
  !close ground state psf
  call close_psf

  !loop over files
  do ifile=1,nfile
     if(ierr.LT.10)index='00'//trim(int2str(ifile))
     if(ierr.GE.10)index='0'//trim(int2str(ifile))

     index='001' !hard code for example purpose only

     pdb='../gnu/share/examples/molreader/lh3_'//index//'.pdb'
     if(check(pdb).EQ.1)call stop('cannot find file '//trim(pdb))

     !write(*,*)myid,ifile,atm_select%segName//'-'//trim(int2str(atm_select%resSeq))

     !open pdb
     call open_pdb(trim(pdb))

     call Prompt('Looking for atom matches!')
     !find matches
     natom=0
     nmatch=0
     mask=.false.
     do while(trim(Next_Record()).NE.'END')
        if(trim(Next_Record()).EQ.'ATOM')then

           call read(atm)
           natom=natom+1

           do istate=1,nstate

              match_res=.False.
              match_resname=.False.
              match_segname=.False.

              if(atm%resname.EQ.atm_select(istate)%resname)match_resname=.True.
              if(atm%resSeq.EQ.atm_select(istate)%resSeq)match_res=.True.
              if(atm%segName.EQ.atm_select(istate)%segName)match_segname=.True.

              match=match_resname.and.match_res.and.match_segname

              if(match)then
                 nmatch=nmatch+1
                 mask(natom,istate)=.true. 
                 !call display(atm) !< uncomment if you want to see the matched atoms

                 !assign excited state charges if available
                 do jatom=1,nchromatm
                    match_name=.false.
                    if(adjustl(psfxatm(jatom)%name).EQ.adjustl(atm%name))&
                         Qx(natom)=psfxatm(jatom)%charge
                 end do
              end if
           end do
        else
           call advance_record
        end if
     end do
     !check for the same number of atoms
     if(npsfatom.NE.natom)call stop('number of atoms in psf does not equall the number of atoms in the pdb')
     call Prompt('atom matches= '//trim(int2str(nmatch)))
     
     !-------begin analysis--------
     !rewind the pdb 
     call rewind_pdb

     !loop over frames
     nframe=0
     do while(len(trim(Next_Record())).NE.0)
        iatom=0
        !new frame
        call Prompt('Processing File_Frame: '//trim(int2str(ifile))//'_'//trim(int2str(nframe)))
        do while(trim(Next_Record()).NE.'END')
           if(trim(Next_Record()).EQ.'ATOM')then
              call read(atm)
              iatom=iatom+1
              XYZ(1,iatom)=atm%X
              XYZ(2,iatom)=atm%Y
              XYZ(3,iatom)=atm%Z

              !for debugging only
              !if(mask(iatom))then
              !   call display(atm)
              !   write(*,*)Qg(iatom),Qx(iatom)
              !   stop
              !end if

           else
              call advance_record
           end if

        end do
        if(trim(Next_Record()).EQ.'END')then
           call advance_record
           nframe=nframe+1
           if(iatom.NE.natom)call stop('number of atoms in frame '&
                //trim(int2str(nframe))//' not equal to natoms!')
        end if

        if(iatom.NE.0)then
           !accumulate sum delE
           delE=0._double
           do istate=1,nstate
              !N loop over atoms k (pigment mask)
              do katom=1,natom
                 if(mask(katom,istate))then
                    delQ=Qx(katom)-Qg(katom)
                    if(usedifQ)delQ=Qx(katom)
                    !M loop over atoms j (not pigment mask)
                    do jatom=1,natom
                       if(.not.mask(jatom,istate))then
                          !compute dist |Rij|
                          R=sqrt(sum((XYZ(:,katom)-XYZ(:,jatom))**2))
                          if(R.EQ.0)call stop('pigment atom '//trim(int2str(katom))&
                               //' and env atom '//trim(int2str(jatom))&
                               //' sit on top of each other, distance = 0') 
                          delE(istate)=delE(istate)+delQ*Qg(jatom)/R
                       end if
                    end do
                 end if
              end do
           end do
           delE=kc*delE
           inquire(unit,opened=match)
           if(.not.match)open(unit,file=trim(outfile),access='append')
           write(unit,'(I5,1X,'//trim(int2str(nstate))//'(F16.8,1X))')(ifile-1)*fpf+nframe&
                ,(delE(istate)/invcm,istate=1,nstate)
           !write(*,'(I5,1X,'//trim(int2str(nstate))//'(F16.8,1X))')(ifile-1)*fpf+nframe&
           !     ,(delE(istate)/invcm,istate=1,nstate)
           !stop
        end if
     end do
     call close_pdb
  end do
  close(unit)
!-----------------------------------------------------------------
call closeLog

end program
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!function to fit
FUNCTION F(param,n)
  use type_kinds
  use filemanager
  use errorlog
  
  INTEGER(long),intent(in)::n
  REAL(double),intent(inout):: param(n)
  INTEGER(long):: i,ierr,unit
  REAL(double):: F,xobs,yobs,ycal
  !character(len=path)::file='example.dat'
  F=0.0

!!$  if(check(trim(file)).NE.0)call stop('cannot find file '//trim(file))
!!$  unit=newunit()
!!$  open(unit,file=trim(file))
!!$  ierr=0
!!$  DO while(ierr.eq.0)
!!$     read(unit,*,iostat=ierr)xobs,yobs
!!$     if(ierr.EQ.0)then
!!$        ycal=param(1)*xobs*exp(-param(2)*xobs) !example fit to ohmic with debey cuttoff
!!$        F=F+(ycal-yobs)**2         !minimize error calulated as sum of the diference squared
!!$     end if
!!$  END DO
!!$  close(unit)  
  

END FUNCTION F
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
