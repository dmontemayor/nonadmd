!---------------------------------------------------
!>\brief
!! Molecule Building Utility
!!\details
!! Supports PDB - Atomic Coordinate Entry Format Version 3.3
!! plus some extra PDB entries typical of some MD packages. 
!!\authors
!!Daniel Montemayor
!!\date
!! 23 Feb 2012
!!\bug
!! psf file cannont have extra charaters in 'section' line. (go to advance_psf subrtouine)
!!\todo
!! depreciate next_record, and advance_record use instead next_pdb and advance_pdb
!<---------------------------------------------------
module molreader
  use type_kinds
  use atomicunits
  use filemanager
  use string
  use errorlog
  use wallclock
  implicit none
  private

  !Public types
  public::CRYST1
  public::MODEL,ATOM,psfATOM

  !public methods
  public::read,display,new
  public::next_Record,advance_Record
  public::open_pdb,rewind_pdb,close_pdb
  public::open_psf,rewind_psf,advance_psf,close_psf,next_psf

  !. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
  !file parameters
  integer(long)::pdbunit=1000,psfunit=2000,parsingunit=100
  character(len=path)::pdbfile='',psffile=''
  logical::pdbisopen,psfisopen
  !integer::natoms,nbonds,ndihedrals,nimpropers,ncrossterms

  !buffer parameters
  integer(long),parameter::REClen=6,Linelen=80,maxrecord=1000
  character(len=5)::lineFMT='(A80)'
  character(len=LineLen)::buffer(maxrecord),psfbuffer(maxrecord)
  integer(long)::buffer_cursor,buffer_count
  integer(long)::psfbuffer_cursor,psfbuffer_count
  logical::buffermask(maxrecord),psfbuffermask(maxrecord)
  logical::recording_parsingdata=.False.!>to save parsing timing data
  real(double)::parsing_rate,time0

  type CRYST1
     character(len=30)::FMT="(A6,3(F9.3),3(F7.2),1X,A11,I4)"
     character(len=REClen)::RECNAME
     real(double)::a,b,c
     real(double)::alpha,beta,gamma
     character(len=11)::sGroup
     integer(long)::z
  end type CRYST1

  type MODEL
     character(len=11)::FMT="(A6,4X,I4)"
     character(len=REClen)::RECNAME
     integer(long)::serial
  end type MODEL

  type ATOM
    character(len=62)::FMT="(A6,A5,1X,A4,A1,A4,A1,A4,A1,3X,3(F8.3),2(F6.2),6X,A4,A2,A2)"
     character(len=REClen)::RECNAME
     integer(long)::serial
     character(len=4)::name
     character::altLoc
     character(len=4)::resname 
     character::chainID
     integer(long)::resSeq
     character::iCode
     real(double)::x,y,z
     real(double)::occupancy
     real(double)::tempFactor
     character(len=4)::segName !<found from some MD packages
     character(len=2)::element 
     character(len=2)::charge
  end type ATOM

  type psfATOM
     character(len=58):: FMT="(I8,1X,A4,1X,I4,1X,A4,1X,A4,1X,A4,1X,F10.6,6X,F8.4,11X,I1)"
     integer(long)::serial
     character(len=4)::segName
     integer(long)::resSeq
     character(len=4)::resname
     character(len=4)::name
     character(len=4)::type
     real(double)::charge
     real(double)::mass
  end type PSFATOM

  interface advance_Record
     module procedure advance_buffer
  end interface

  interface read
     module procedure READ_CRYST1
     module procedure READ_MODEL
     module procedure READ_ATOM
     module procedure READ_psfATOM
  end interface

  interface display
     module procedure display_CRYST1
     module procedure display_MODEL
     module procedure display_ATOM
     module procedure display_psfATOM
  end interface

  interface new
     module procedure new_CRYST1
     module procedure new_MODEL
     module procedure new_ATOM
     module procedure new_psfATOM
  end interface

contains
!!\brief Retrun Next PDB Record name.
!!\details
!! This Function will return the next PDB record name from opened PDB file.
!! Does not advance the buffer! 
!<-------------------------------------------------------------
  function Next_Record()
    character(len=REClen)::Next_Record

    Next_Record=buffer(buffer_cursor)(1:6)

    return
  end function Next_Record
!----------------------------------------------------------------
!!\brief Reads and Returns CRYST1 Record type.
!!\details
!!The CRYST1 record presents the unit cell parameters, space group, and Z value. If the structure was not determined by crystallographic means, CRYST1 simply provides the unitary values, with an appropriate REMARK.
!!
!!Record Format
!!
!!COLUMNS       DATA  TYPE    FIELD          DEFINITION
!!-------------------------------------------------------------
!! 1 -  6       Record name   "CRYST1"
!! 7 - 15       Real(9.3)     a              a (Angstroms).
!!16 - 24       Real(9.3)     b              b (Angstroms).
!!25 - 33       Real(9.3)     c              c (Angstroms).
!!34 - 40       Real(7.2)     alpha          alpha (degrees).
!!41 - 47       Real(7.2)     beta           beta (degrees).
!!48 - 54       Real(7.2)     gamma          gamma (degrees).
!!55            Character       -            unused space.
!!56 - 66       LString(11)   sGroup         Space  group.
!!67 - 70       Integer       z              Z value.
!<-------------------------------------------------------------
  subroutine READ_CRYST1(this)
    type(CRYST1),intent(inout)::this
    integer(long)::ierr
    character(len=70)::line
    call new(this)

    line=buffer(buffer_cursor)(1:70)
 
    this%RecNAME=line(1:6)
    this%a =str2float(line(7:15))*Angstrom
    this%b =str2float(line(16:24))*Angstrom
    this%c =str2float(line(25:33))*Angstrom
    this%alpha =str2float(line(34:40))*Deg
    this%beta =str2float(line(41:47))*Deg
    this%gamma =str2float(line(48:54))*Deg
    this%sGroup =line(56:66)
    this%z =str2int(line(67:70))

    call advance_buffer

  end subroutine READ_CRYST1
!----------------------------------------------------------------
  subroutine display_CRYST1(this,unit)
    type(CRYST1),intent(inout)::this
    integer(long),optional::unit
    integer(long)::ierr,u

    u=6
    if(present(unit))u=unit
    write(u,this%FMT,iostat=ierr)this%RECNAME,this%a/Angstrom,this%b/Angstrom,this%c/Angstrom&
         ,this%alpha/Deg,this%beta/Deg,this%gamma/Deg,this%sGroup,this%z

    if(ierr)call warn('WRITE_CRYST1: could not write record!')
  end subroutine DISPLAY_CRYST1
!----------------------------------------------------------------
  subroutine new_CRYST1(this)
    type(CRYST1),intent(inout)::this
     this%RECNAME=''
     this%a=0_double
     this%b=0_double
     this%c=0_double
     this%alpha=0_double
     this%beta=0_double
     this%gamma=0_double
     this%sGroup=''
     this%z=0
   end subroutine new_CRYST1
!----------------------------------------------------------------
!>\brief Reads and returns MODEL Record type
!!\details
!!The MODEL record specifies the model serial number when multiple models of the same structure are presented in a single coordinate entry, as is often the case with structures determined by NMR.
!!
!!Record Format
!!
!!COLUMNS        DATA  TYPE    FIELD          DEFINITION
!!----------------------------------------------------------------
!! 1 -  6        Record name   "MODEL "
!! 7 - 10        LString(4)      -            unused space.
!!11 - 14        Integer       serial         Model serial number.
!<----------------------------------------------------------------
  subroutine READ_MODEL(this)
    type(MODEL),intent(inout)::this
    integer(long)::ierr
    character(len=14)::line
    call new(this)

    line=buffer(buffer_cursor)(1:14)

    this%RECNAME=line(1:6)
    this%serial=str2int(line(11:14))
    
    call advance_buffer

  end subroutine READ_MODEL
!----------------------------------------------------------------
  subroutine display_MODEL(this,unit)
    type(MODEL),intent(inout)::this
    integer(long),optional::unit
    integer(long)::ierr,u

    u=6
    if(present(unit))u=unit
    write(u,this%FMT,iostat=ierr)this%RECNAME,this%serial
    if(ierr)call warn('WRITE_MODEL: could not write record!')

  end subroutine DISPLAY_MODEL
!----------------------------------------------------------------
  subroutine new_Model(this)
    type(Model),intent(inout)::this
     this%RECNAME='MODEL'
     this%serial=0
   end subroutine new_Model
!----------------------------------------------------------------
!>\brief Reads and returns ATOM/HETATM/TER Record types
!!\details
!!The ATOM records present the atomic coordinates for standard amino acids and nucleotides.
!!They also present the occupancy and temperature factor for each atom.
!!Non-polymer chemical coordinates use the HETATM record type.
!!The element symbol is always present on each ATOM record; charge is optional.
!!
!!Changes in ATOM/HETATM records result from the standardization atom and residue nomenclature.
!!This nomenclature is described in the Chemical Component Dictionary (ftp://ftp.wwpdb.org/pub/pdb/data/monomers).
!!
!!Record Format
!!
!!COLUMNS        DATA  TYPE    FIELD        DEFINITION
!!-------------------------------------------------------------------------------------
!! 1 -  6        Record name   "ATOM  "
!! 7 - 11        Integer       serial       Atom  serial number.
!!12             Character       -          unused space.
!!13 - 16        Atom          name         Atom name.
!!17             Character     altLoc       Alternate location indicator.
!!18 - 20        Residue name  resName      Residue name.
!!21             Character       -          unused space.
!!22             Character     chainID      Chain identifier.
!!23 - 26        Integer       resSeq       Residue sequence number.
!!27             AChar         iCode        Code for insertion of residues.
!!28 - 30        LString(3)      -          unused space.
!!31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
!!39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
!!47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
!!55 - 60        Real(6.2)     occupancy    Occupancy.
!!61 - 66        Real(6.2)     tempFactor   Temperature  factor.
!!67 - 72        LString(6)      -          Unused section.
!!73 - 76        LString(4)    segName      Segment Name.
!!77 - 78        LString(2)    element      Element symbol, right-justified.
!!79 - 80        LString(2)    charge       Charge  on the atom.
!!=============================================================================
!!Non-polymer or other “non-standard” chemical coordinates,
!!such as water molecules or atoms presented in HET groups use the HETATM record type.
!!They also present the occupancy and temperature factor for each atom.
!!The ATOM records present the atomic coordinates for standard residues.
!!The element symbol is always present on each HETATM record; charge is optional.
!!
!!Changes in ATOM/HETATM records will require standardization in atom and residue nomenclature.
!!This nomenclature is described in the Chemical Component Dictionary, ftp://ftp.wwpdb.org/pub/pdb/data/monomers.
!!
!!Record Format
!!
!!COLUMNS       DATA  TYPE     FIELD         DEFINITION
!!-----------------------------------------------------------------------
!! 1 - 6        Record name    "HETATM"
!! 7 - 11       Integer        serial        Atom serial number.
!!13 - 16       Atom           name          Atom name.
!!17            Character      altLoc        Alternate location indicator.
!!18 - 20       Residue name   resName       Residue name.
!!22            Character      chainID       Chain identifier.
!!23 - 26       Integer        resSeq        Residue sequence number.
!!27            AChar          iCode         Code for insertion of residues.
!!31 - 38       Real(8.3)      x             Orthogonal coordinates for X.
!!39 - 46       Real(8.3)      y             Orthogonal coordinates for Y.
!!47 - 54       Real(8.3)      z             Orthogonal coordinates for Z.
!!55 - 60       Real(6.2)      occupancy     Occupancy.
!!61 - 66       Real(6.2)      tempFactor    Temperature factor.
!!77 - 78       LString(2)     element       Element symbol; right-justified.
!!79 - 80       LString(2)     charge        Charge on the atom.
!!=============================================================================
!!The TER record indicates the end of a list of ATOM/HETATM records for a chain.
!!
!!Record Format
!!
!!COLUMNS        DATA  TYPE    FIELD           DEFINITION
!!-------------------------------------------------------------------------
!! 1 -  6        Record name   "TER   "
!! 7 - 11        Integer       serial          Serial number.
!!12             Character       -          unused space.
!!13 - 16        LString(4)      -          unused space.
!!17             Character       -          unused space.
!!18 - 20        Residue name  resName      Residue name.
!!21             Character       -          unused space.
!!22             Character     chainID      Chain identifier.
!!23 - 26        Integer       resSeq       Residue sequence number.
!!27             AChar         iCode        Insertion code.
!!
!! Notice ATOM, HETATM, and TER can all be defined by the ATOM type
!<---------------------------------------------------------------
  subroutine READ_ATOM(this)
    type(ATOM),intent(inout)::this
    integer::ierr
    character(len=80)::line
    call new(this)
    

    line=buffer(buffer_cursor)(1:80)

    this%RECNAME=line(1:6)
    this%serial=str2int(line(7:11))
    this%name=line(13:16)
    this%altLoc=line(17:17)
    this%resNAME=line(18:20)
    this%chainID=line(22:22)
    this%resSeq=str2int(line(23:26))
    this%iCode=line(27:27)
    this%x=str2float(line(31:38))*Angstrom
    this%y=str2float(line(39:46))*Angstrom
    this%z=str2float(line(47:54))*Angstrom
    this%occupancy=str2float(line(55:60))
    this%tempFactor=str2float(line(61:66))
    this%segName=line(73:76)
    this%element=line(77:78)
    this%charge=line(79:80)

    call advance_buffer

    return
  end subroutine READ_ATOM
!----------------------------------------------------------------
  subroutine new_ATOM(this)
    type(ATOM),intent(inout)::this
    this%RECNAME='ATOM  '
    this%serial=0
    this%name='    '
    this%altLoc=' '
    this%resname='    '
    this%chainID=' '
    this%resSeq=0
    this%iCode=' '
    this%x=0_double
    this%y=0_double
    this%z=0_double
    this%occupancy=0_double
    this%tempFactor=0_double
    this%segName='    '
    this%element='  ' 
    this%charge='  '
    return
  end subroutine new_ATOM
!----------------------------------------------------------------
  subroutine display_ATOM(this,unit)
    type(ATOM),intent(inout)::this
    integer(long),optional::unit
    integer(long)::ierr,u
    character(len=5)::serial
    character(len=4)::resSeq

    u=6
    if(present(unit))u=unit

    if(this%serial.GT.99999)then
       serial=trim(int2hex(this%serial))
    else
       serial=trim(int2str(this%serial))
    end if
    serial=adjustr(serial)

    if(this%resSeq.GT.9999)then
       resSeq=trim(int2hex(this%resSeq))
    else
       resSeq=trim(int2str(this%resSeq))
    end if
    resSeq=adjustr(resSeq)

    write(u,this%FMT,iostat=ierr)this%RECNAME,serial,this%name&
         ,this%altLoc,this%resName,this%chainID,resSeq,this%iCode&
         ,this%x/Angstrom,this%y/Angstrom,this%z/Angstrom,this%occupancy,this%tempFactor&
         ,this%segName,this%element,this%charge

    if(ierr)call warn('WRITE_ATOM: could not write record!')
  end subroutine DISPLAY_ATOM
!----------------------------------------------------------------
!>\brief Initialize PDB file for parsing
!<---------------------------------------------------------------
  subroutine open_pdb(file)
    character*(*),intent(in)::file
    integer(long)::unit,ierr,init_time(8)
    logical::pass

    !check file exists
    if(check(file))then
       call warn('open_pdb: cannot find input file '//file,'not opening file.')
       pdbisopen=.false.
    end if

    !assign unit
    unit=newunit()
 
    !openfile
    open(unit,file=file)

    pass=.True.
    !check file is correct format
    if(.not.pass)then
       call warn('open_pdb: cannot read file '//file,'closing file.')
       pdbisopen=.false.
       close(unit)
       close(parsingunit)
    else
       pdbisopen=.true.
       
       !save unit and associated file
       pdbunit=unit
       pdbfile=file
       init_time=timearray()
       time0=(init_time(8)*1E-3+init_time(7)&
            +init_time(6)*60_double+init_time(5)*3600_double)
       parsingunit=newunit()
       if(recording_parsingdata)open(parsingunit,file='parsing.dat')

       !load buffer data
       buffer_count=0
       call refresh_Buffer
    end if

    if(pdbisopen)call NOTE('PDB file '//file//' successfully opened.')

    return
  end subroutine open_pdb
!----------------------------------------------------------------
  subroutine close_pdb
    if(.not.pdbisopen)then
       call NOTE('PDB file '//trim(pdbfile)//' is already closed.')
    else
       close(pdbunit)
       pdbisopen=.false.
       call NOTE('PDB file '//trim(pdbfile)//' successfully closed.')
    end if
  end subroutine close_pdb
!----------------------------------------------------------------
!>\brief Rewind an initialized PDB file for parsing
!<---------------------------------------------------------------
  subroutine rewind_pdb
    integer(long)::ierr
    logical::pass

    !check file is open
    if(.not.pdbisopen)call warn('rewind_pdb:'//trim(pdbfile)//' is not open.')

    if(pdbisopen)then
       !rewind file and load buffer data
       rewind(pdbunit,iostat=ierr)
       if(ierr)then
          call warn('rewind_pdb: cannot rewind file.')
       else
          buffer_count=0
          call refresh_Buffer
          call NOTE('PDB file '//trim(pdbfile)//' successfully rewound.')
       end if
    end if
    
  end subroutine rewind_pdb
!------------------------------------------------------------------
!>\details reads the next chunk of the pdb file.
  subroutine refresh_buffer
    integer(long)::irecord,ierr,curr_time(8),pdblinecount
    real(double)::elapsedtime,time1
    buffer=''
    buffermask=.false. !valid lines to read
    do irecord=1,maxrecord
       read(pdbunit,lineFMT,iostat=ierr)buffer(irecord)
       if(.not.ierr)buffermask(irecord)=.true.
    end do
    !close pdb file if no elements of the buffer are readable
    if(all(.not.buffermask))then
       call note('buffer is empty!')
       close(pdbunit)
       close(parsingunit)
       pdbisopen=.false.
       call note('closing pdb file '//trim(pdbfile))
    end if
    buffer_cursor=1
    buffer_count=buffer_count+1
    if(recording_parsingdata.and.buffer_count.GT.1)then
       curr_time=timearray()
       time1=(curr_time(8)*1E-3+curr_time(7)&
            +curr_time(6)*60_double+curr_time(5)*3600_double)
       elapsedtime=time1-time0
       pdblinecount=buffer_count*maxrecord
       parsing_rate=pdblinecount/elapsedtime
       call note('Refreashing buffer. Parsing rate in records/second= '&
            //float2str(parsing_rate))
       write(parsingunit,'(3(F16.9))')time1-time0,real(pdblinecount),parsing_rate
    end if

  end subroutine refresh_buffer
!----------------------
  subroutine advance_buffer

    buffer_cursor=buffer_cursor+1
    if(buffer_cursor.GT.maxrecord)call refresh_buffer
    do while(.not.buffermask(buffer_cursor).and.pdbisopen)
       buffer_cursor=buffer_cursor+1
       if(buffer_cursor.GT.maxrecord)call refresh_buffer
    end do
    
  end subroutine advance_buffer
!-----------------------------------------
!< \brief
!! Exctract PSF atom information
!! \authors
!! Daniel Montemayor
!! \date
!! 20 Feb 2012
!<------------------------------------------
  subroutine READ_psfatom(this)
    type(psfatom),intent(inout)::this
    integer(long)::ierr
    character(len=Linelen)::line
    call new(this)
    
    !          1         2         3         4         5         6         7
    !01234567890123456789012345678901234567890123456789012345678901234567890123456789
    !!IIIIIIIIXAAAAXIIIIXAAAAXAAAAXAAAAX+FF.ddddddXXXXXX+FF.ddddXXXXXXXXXXXI
    !  serial  sgnm rsid rsnm name type   charge          mass
    !character(len=58):: FMT="(I8,1X,A4,1X,I4,1X,A4,1X,A4,1X,A4,1X,F10.6,6X,F8.4,11X,I1)"<---this is the format
    line=psfbuffer(psfbuffer_cursor)(1:lineLen)
    
    this%serial=str2int(line(1:8))
    this%segname=line(10:13)
    this%resSeq=str2int(line(15:18))
    this%resname=line(20:23)
    this%name=line(25:28)
    this%type=line(30:33)
    this%charge=str2float(line(35:44))*e
    this%mass=str2float(line(51:58))*mp
    
    call advance_psfbuffer
    
  end subroutine READ_PSFATOM
  !----------------------
  subroutine display_psfatom(this,unit)
    type(psfatom),intent(in)::this
    integer(long),optional::unit
    integer(long)::ierr,u

    u=6
    if(present(unit))u=unit
    write(u,this%FMT,iostat=ierr)this%serial,this%segname&
         ,this%resSeq,this%resname,this%name,this%type&
         ,this%charge/e,this%mass/mp,0
    if(ierr)call warn('display: could not display psf atom!')

  end subroutine DISPLAY_PSFATOM
  !----------------------
  subroutine new_psfatom(this)
    type(psfatom),intent(inout)::this
    integer(long)::ierr
    
    this%serial=0
    this%segname=''
    this%resSeq=0
    this%resname=''
    this%name=''
    this%type=''
    this%charge=0_double
    this%mass=0_double
        
  end subroutine NEW_PSFATOM
!!\brief Retrun Next PDB Record name.
!!\details
!! This Function will return the next PDB record name from opened PDB file.
!! Does not advance the buffer! 
!<-------------------------------------------------------------
  function next_psf()
    character(len=Linelen)::next_psf

    next_psf=psfbuffer(psfbuffer_cursor)(10:Linelen)

    return
  end function next_psf
!----------------------------------------------------------------

  !----------------------
  subroutine advance_psfbuffer

    psfbuffer_cursor=psfbuffer_cursor+1
    if(psfbuffer_cursor.GT.maxrecord)call refresh_psfbuffer
    do while(.not.psfbuffermask(psfbuffer_cursor).and.psfisopen)
       psfbuffer_cursor=psfbuffer_cursor+1
       if(psfbuffer_cursor.GT.maxrecord)call refresh_psfbuffer
    end do
    
  end subroutine advance_psfbuffer
!------------------------------------------------------------------
!>\details reads the next chunk of the psf file.
  subroutine refresh_psfbuffer
    integer(long)::irecord,ierr,psflinecount
    psfbuffer=''
    psfbuffermask=.false. !valid lines to read
    do irecord=1,maxrecord
       read(psfunit,lineFMT,iostat=ierr)psfbuffer(irecord)
       if(.not.ierr)psfbuffermask(irecord)=.true.
    end do
    !close psf file if no elements of the psfbuffer are readable
    if(all(.not.psfbuffermask))then
       call note('psfbuffer is empty!')
       close(psfunit)
       psfisopen=.false.
       call note('closing psf file '//trim(psffile))
    end if
    psfbuffer_cursor=1
    psfbuffer_count=psfbuffer_count+1
  end subroutine refresh_psfbuffer
!----------------------------------------------------------------
!>\brief Initialize PSF file for parsing
!<---------------------------------------------------------------
  subroutine open_psf(file)
    character*(*),intent(in)::file
    character(len=3)::TITLE
    integer(long)::unit,ierr
    logical::pass

    !check file exists
    if(check(file))then
       call warn('open_psf: cannot find input file '//file,'not opening file.')
       psfisopen=.false.
    end if

    !assign unit
    unit=newunit()
 
    !openfile
    open(unit,file=file)

    pass=.True.
    !check file is correct format
    read(unit,'(A3)')TITLE
    if(TITLE.NE.'PSF')pass=.false.
    if(.not.pass)then
       call warn('open_psf: cannot read file '//file,'closing file.')
       psfisopen=.false.
       close(unit)
    else
       psfisopen=.true.
       
       !save unit and associated file
       psfunit=unit
       psffile=file

       !load buffer data
       psfbuffer_count=0
       call refresh_psfBuffer
    end if

    if(psfisopen)call NOTE('PSF file '//file//' successfully opened.')

    return
  end subroutine open_psf
!----------------------------------------------------------------
  subroutine close_psf
    if(.not.psfisopen)then
       call NOTE('PSF file '//trim(psffile)//' is already closed.')
    else
       close(psfunit)
       psfisopen=.false.
       call NOTE('PSF file '//trim(psffile)//' successfully closed.')
    end if
  end subroutine close_psf
!----------------------------------------------------------------
!>\brief Rewind an initialized PSF file for parsing
!<---------------------------------------------------------------
  subroutine rewind_psf
    integer(long)::ierr
    logical::pass

    !check file is open
    if(.not.psfisopen)call warn('rewind_psf:'//trim(psffile)//'is not open.')

    if(psfisopen)then
       !rewind file and load buffer data
       rewind(psfunit,iostat=ierr)
       if(ierr)then
          call warn('rewind_psf: cannot rewind file.')
       else
          psfbuffer_count=0
          call refresh_psfBuffer
          call NOTE('PSF file '//trim(psffile)//' successfully rewound.')
       end if
    end if
    
  end subroutine rewind_psf
!----------------------------------------------------------------
!>\brief
!!Advances opened psf file
!!\details
!!Advances one line or optionally to one of the following sections:
!! title, atoms, bonds, angles, dihedrals, impropers, or cross-terms
!!\bug
!!trim(str).ne.'section': means that nothing else can be written on section line,
!! the bug can be fixed by selecting the columns of the next_psf() string not the trim. 
!<---------------------------------------------------------------
  subroutine advance_psf(section,nrecord)
    character*(*),intent(in),optional::section
    integer(long),optional::nrecord
    integer(long)::N
    integer(long)::ititle

    if(psfisopen)then

       if(present(section))then
          select case(section)

          case('title')
             !skip to next section
             do while(trim(next_psf()).NE.'!NTITLE'.and.psfisopen)
                call advance_psfbuffer
             end do
             if(psfisopen)then
                N=str2int(psfbuffer(psfbuffer_cursor)(1:8))
                if(present(nrecord))nrecord=N
                call note('advanced psf to title section.') 
                call note('advance_psf found '//trim(int2str(N))//' remarks.') 
                call advance_psfbuffer
             else
                call warn('advance_psf did not find title section.','closing file.') 
                call close_psf
             end if

          case('atoms')
             !skip to next section
             do while(trim(next_psf()).NE.'!NATOM'.and.psfisopen)
                call advance_psfbuffer
             end do
             if(psfisopen)then
                N=str2int(psfbuffer(psfbuffer_cursor)(1:8))
                if(present(nrecord))nrecord=N
                call note('advanced psf to atoms section.') 
                call note('advance_psf found '//trim(int2str(N))//' atoms.') 
                call advance_psfbuffer
             else
                call warn('advance_psf did not find atoms section.','closing file.') 
                call close_psf
             end if

          case('bonds')
             !skip to next section
             do while(trim(next_psf()).NE.'!NBOND'.and.psfisopen)
                call advance_psfbuffer
             end do
             if(psfisopen)then
                N=str2int(psfbuffer(psfbuffer_cursor)(1:8))
                if(present(nrecord))nrecord=N
                call note('advanced psf to bonds section.') 
                call note('advance_psf found '//trim(int2str(N))//' bonds.') 
                call advance_psfbuffer
             else
                call warn('advance_psf did not find bonds section.','closing file.') 
                call close_psf
             end if

          case('angles')
             !skip to next section
             do while(trim(next_psf()).NE.'!NTHETA'.and.psfisopen)
                call advance_psfbuffer
             end do
             if(psfisopen)then
                N=str2int(psfbuffer(psfbuffer_cursor)(1:8))
                if(present(nrecord))nrecord=N
                call note('advanced psf to angles section.') 
                call note('advance_psf found '//trim(int2str(N))//' angles.') 
                call advance_psfbuffer
             else
                call warn('advance_psf did not find angles section.','closing file.') 
                call close_psf
             end if

          case('dihedrals')
             !skip to next section
             do while(trim(next_psf()).NE.'!NPHI'.and.psfisopen)
                call advance_psfbuffer
             end do
             if(psfisopen)then
                N=str2int(psfbuffer(psfbuffer_cursor)(1:8))
                if(present(nrecord))nrecord=N
                call note('advanced psf to dihedrals section.') 
                call note('advance_psf found '//trim(int2str(N))//' dihedrals.') 
                call advance_psfbuffer
             else
                call warn('advance_psf did not find dihedrals section.','closing file.') 
                call close_psf
             end if

          case('impropers')
             !skip to next section
             do while(trim(next_psf()).NE.'!NIMPHI'.and.psfisopen)
                call advance_psfbuffer
             end do
             if(psfisopen)then
                N=str2int(psfbuffer(psfbuffer_cursor)(1:8))
                if(present(nrecord))nrecord=N
                call note('advanced psf to impropers section.') 
                call note('advance_psf found '//trim(int2str(N))//' impropers.') 
                call advance_psfbuffer
             else
                call warn('advance_psf did not find impropers section.','closing file.') 
                call close_psf
             end if

          case('cross-terms')
             !skip to next section
             do while(trim(next_psf()).NE.'!NCRTERM'.and.psfisopen)
                call advance_psfbuffer
             end do
             if(psfisopen)then
                N=str2int(psfbuffer(psfbuffer_cursor)(1:8))
                if(present(nrecord))nrecord=N
                call note('advanced psf to cross-terms section.') 
                call note('advance_psf found '//trim(int2str(N))//' cross-terms.') 
                call advance_psfbuffer
             else
                call warn('advance_psf did not find cross-terms section.','closing file.') 
                call close_psf
             end if

          case default
             call warn('advance_psf: unknown section! Valid options are: title, atoms, bonds, angles, dihedrals, impropers, or cross-terms','ignorning advance_psf command.')
          end select
       else
          !advance buffer by one record
          call advance_psfbuffer
       end if
    else
       call warn('advance_psf: PSF file is not open.','ignorning advance_psf command.')
    end if
  end subroutine advance_psf
!!$  subroutine advance_psf(section,nrecord)
!!$    character*(*),intent(in),optional::section
!!$    integer(long),optional::nrecord
!!$    integer(long)::Ntitle,Natom
!!$    integer(long)::ititle
!!$
!!$    if(present(section))then
!!$       select case(section)
!!$       case('atoms')
!!$          call rewind_psf
!!$          !read PSF title
!!$          call advance_psfbuffer
!!$
!!$          !read blank line
!!$          call advance_psfbuffer
!!$
!!$          !read NTITLE
!!$          Ntitle=str2int(psfbuffer(psfbuffer_cursor)(1:8))
!!$          !read(psfbuffer(psfbuffer_cursor),*,iostat=ierr)Ntitle          
!!$          call prompt('advance_psf found '//trim(int2str(Ntitle))//' REMARKS') 
!!$          call advance_psfbuffer
!!$
!!$          !skip Ntitle lines
!!$          do ititle=1,Ntitle
!!$             call advance_psfbuffer
!!$          end do
!!$ 
!!$          !read blank line
!!$          call advance_psfbuffer
!!$
!!$          !read and prompt NATOM
!!$          Natom=str2int(psfbuffer(psfbuffer_cursor)(1:8))
!!$          !read(psfbuffer(psfbuffer_cursor),*)Natom
!!$          if(present(nrecord))nrecord=Natom
!!$          call prompt('advance_psf found '//trim(int2str(Natom))//' atoms') 
!!$          call prompt('advanced psf to atoms section.') 
!!$          call advance_psfbuffer
!!$
!!$       case('bonds')
!!$          call warn('advance_psf: cannot advance to bonds section! Please advance to section manually.','doing nothing.')
!!$       case('angles')
!!$          call warn('advance_psf: cannot advance to anngles section! Please advance to section manually.','doing nothing.')
!!$       case('dihedrals')
!!$          call warn('advance_psf: cannot advance to dihedrals section! Please advance to section manually.','doing nothing.')
!!$       case('impropers')
!!$          call warn('advance_psf: cannot advance to impropers section! Please advance to section manually.','doing nothing.')
!!$       case('cross-terms')
!!$          call warn('advance_psf: cannot advance to cross-terms section! Please advance to section manually.','doing nothing.')
!!$       case default
!!$          call warn('advance_psf: unknown section!','doing nothing.')
!!$       end select
!!$    else
!!$       !advance buffer by one record
!!$       call advance_psfbuffer
!!$    end if
!!$  end subroutine advance_psf

end module molreader

