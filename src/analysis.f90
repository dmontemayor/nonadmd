!==================================================================================================================================!
!                                                                                                                                  !
! THIS MODULE CONTAINS ALL SUBROUTINES NEEDED TO COMPUTE AUTO- AND CROSS-CORRELATION FUNCTIONS OF SITE ENERGIES, THEIR             !
! AVERAGE VALUE AND THEIR PROBABILITY DISTRIBUTION FUNCTION.                                                                       !
!                                                                                                                                  !
! COMMENTS:                                                                                                                        !
! 1 - The site energies must be contained on a file named results.out, where the step number and the energy are given in nearby    !
!     columns; for example: step(site 1), en(site 1), step(site 2), en(site 2) ....                                                !
! 2 - The file results.out must have a fixed length (nconfs < ndata), and cannot contain columns with different length. This is    !
!     an important warning because it could happen that different samples are computed for any chromophore. In this case the file  !
!     length should be equal to the minimum number of samples;                                                                     !
! 3 - Given the format of the input it is possible that some time steps are not included for some site; it is fundamental, though, !
!     that the step number is given in agreement with a single time step;                                                          !
! 4 - The units for input site energies and for output files are given in variables input_units and output_units;                  !
! 5 - it is assumed that the correlation time and the time step are in ps;                                                         !
! 6 - be careful to set the right temperature in the main program, otherwise the SD are wrong!                                     !
!                                                                                                                                  !
!                                                                                                                                  !
! Record of revisions:                                                                                                             !
!     Date        Programmer                                           Description of change                                       !
!   ========     ============      =============================================================================================== !
!   25/05/11      Marco Masia      Original code                                                                                   !
!   01/03/12      Marco Masia      Added the automatic calculation of the Spectral Density through FFT of site energies acf        !
!   05/03/12      Marco Masia      Added input and output units for more flexibility                                               !
!                                                                                                                                  !
!==================================================================================================================================!


MODULE corr_pdf_av_en 

!---------------------------------------------------- Begin Variables Declarations ------------------------------------------------!

implicit none

save

! variables and parameters
integer, parameter :: sgl = selected_real_kind(p = 6, r = 37)
integer, parameter :: dbl = selected_real_kind(p = 13, r = 200)

integer :: i, j, k, norigin, it, it0, nstored, istat, nconfs, npoint, lcor, ntorig, nsites, nhist
integer, parameter :: ndata = 60000, newtorig = 1, nmxpow = 5000 

real(kind=dbl) :: pi, idE, integral, diff, tcor, tstep, dE, temp, cnv
real(kind=dbl), parameter :: ev2wvn = 8065.5d0,            &! eV  ->  1/cm
                             wvn2ev = 1.0d0/8065.5d0        ! 1/cm  ->  eV 

! allocatable arrays
integer, allocatable, dimension(:,:) :: ntime0, step

real(kind=dbl), allocatable, dimension(:) :: ener, av_en, sd_en, freq, re_fft, im_fft  
real(kind=dbl), allocatable, dimension(:,:) :: ccf0, en, en0, en_pdf, acf, power
real(kind=dbl), allocatable, dimension(:,:,:) :: ccf, cnorm

character(len=2) :: input_units, output_units 
character(len=5), allocatable, dimension(:,:) :: corr_type


!------------------------------------------------------ End Variables Declarations ------------------------------------------------!


CONTAINS


!**********************************************************************************************************************************!
!                                                                                                                                  !
!                                                 B E G I N       S U B R O U T I N E S                                            !
!                                                                                                                                  !
!**********************************************************************************************************************************!


!================================================================!
!                                                                !
! THIS SUBROUTINE INITIALIZES VARIABLES AND ALLOCATES ALL ARRAYS !
!                                                                !
!================================================================!

subroutine init

! check the input variables and set them to default values if negative
if (nsites < 0) nsites = 1
if (nhist < 0) nhist = 200
if (tcor < 0.0d0) tcor = 4.0d0 
if (tstep < 0.0d0) tcor = 0.001d0 
if (dE < 0.0d0) dE = 20.0d0

! initialize variables
lcor = int( tcor/tstep ) + 1
ntorig = lcor + 5
idE = 1.0d0/dE

! allocate arrays
allocate ( ntime0(nsites,ntorig) )
allocate ( step(nsites,ndata) )
allocate ( ener(nsites) )
allocate ( av_en(nsites) )
allocate ( sd_en(nsites) )
allocate ( en(nsites,ndata) )
allocate ( en0(ntorig,nsites) )
allocate ( en_pdf(nsites,-nhist:nhist) )
allocate ( ccf(ntorig,nsites,nsites) )
allocate ( ccf0(nsites,nsites) )
allocate ( cnorm(ntorig,nsites,nsites) )
allocate ( corr_type(nsites,nsites) )
allocate ( acf(lcor,nsites) )
allocate ( power(nmxpow,nsites) )
allocate ( freq(nmxpow) )
allocate ( re_fft(nmxpow) )
allocate ( im_fft(nmxpow) )

! initialize arrays
ntime0 = 0
step = 0
ener = 0.0d0
av_en = 0.0d0
sd_en = 0.0d0
en = 0.0d0
en0 = 0.0d0
en_pdf = 0.0d0
ccf = 0.0d0
ccf0 = 0.0d0
cnorm = 0.0d0
acf = 0.0d0
power = 0.0d0
freq = 0.0d0
re_fft = 0.0d0
im_fft = 0.0d0

end subroutine init


!============================================================================!
!                                                                            !
! THIS SUBROUTINE READS THE DATA AND COMPUTES AVERAGE AND STANDARD DEVIATION !
!                                                                            !
!============================================================================!

subroutine read_data 

! change input_units from wavenumbers to ev if needed
if (input_units == 'wn') then 
  cnv = wvn2ev
 elseif (input_units == 'ev') then
  cnv = 1.0d0
endif

open (unit=10, file='results.out', status='old')

nconfs = 0
norigin = 0

DO_READ: do

  nconfs = nconfs + 1

  if (nconfs > ndata) then
    open (unit=20, file='analysis-error.out', status='unknown')
    write (20,'(a)') '#################    E  R  R  O  R    #################'
    write (20,'(a)') '#                                                     #'
    write (20,'(a)') '#                  nconfs > ndata                     #'
    write (20,'(a)') '#                                                     #'
    write (20,'(a)') '#  Please, increment ndata parameter in the source!   #'
    write (20,'(a)') '#                                                     #'
    write (20,'(a)') '##############   E N D   E  R  R  O  R   ##############'
    close (unit=20)
    write (6,*) ' ABNORMAL TERMINATION OF MODULE corr_pdf_av_en!'
    write (6,*) ' SEE FILE analysis-error.out FOR DETAILS.'
    stop
  endif

  read (10,*,IOSTAT=istat) (step(i,nconfs), en(i,nconfs), i = 1,nsites) 
  if (istat < 0) exit DO_READ

  en(:,nconfs) = en(:,nconfs)*cnv

  av_en(:) = av_en(:) + en(:,nconfs)
  sd_en(:) = sd_en(:) + en(:,nconfs)*en(:,nconfs)

enddo DO_READ

av_en = av_en/dble(nconfs - 1)
sd_en = dsqrt( sd_en/dble(nconfs - 1) - av_en*av_en )

close (unit=10)

end subroutine read_data


!===============================================================================================!
!                                                                                               !
! THIS SUBROUTINE COMPUTES THE CORRELATION FUNCTIONS AND THE PROBABILITY DISTRIBUTION FUNCTIONS !
!                                                                                               !
!===============================================================================================!

subroutine compute_cf_pdf

do i = 1,nconfs - 1

  ! compute probability distribution functions for site energies
  do j = 1,nsites
    diff = en(j,i) - av_en(j)
    npoint = nint(diff*idE)
    if (npoint > - nhist .and. npoint < nhist) en_pdf(j,npoint) =  en_pdf(j,npoint) + 1.0d0
  enddo
  
  ! store values at new origin 
  if (mod(i,newtorig) == 0) then
    norigin = norigin + 1
    it0 = mod(norigin - 1,ntorig) + 1
    ntime0(:,it0) = step(:,i) !- 1
    en0(it0,:) = en(:,i)
  endif

  ! compute the site energies auto- and cross-corrrelation functions
  nstored = min(norigin,ntorig)
  do it0 = 1,nstored
    do j = 1,nsites
      do k = 1,nsites
        it = step(k,i) - ntime0(j,it0) + 1
        if (it <= lcor) then
          ccf(it,j,k) = ccf(it,j,k) + (en(k,i) - av_en(k))*(en0(it0,j) - av_en(j)) 
          cnorm(it,j,k) = cnorm(it,j,k) + 1.0d0
        endif
      enddo
    enddo
  enddo

enddo

end subroutine compute_cf_pdf


!=======================================================!
!                                                       !
! THIS SUBROUTINE WRITES OUT AUTO-CORRELATION FUNCTIONS !
!                                                       !
!=======================================================!

subroutine wrt_acf 

open (unit=10, file='acf.dat', status='unknown')

! change output_units from eV to wavenumbers if needed
if (output_units == 'wn') then
  cnv = ev2wvn*ev2wvn
 elseif (output_units == 'ev') then
  cnv = 1.0d0
endif

!ccf0(:,:) = sign( ccf(1,:,:)/cnorm(1,:,:) , 1.0d0 )
ccf0(:,:) = 1.0d0 
do i = 1,lcor 
  do j = 1,nsites
    acf(i,j) = ccf(i,j,j) / ( ccf0(j,j)*cnorm(i,j,j) )
  enddo
  write (10,'(30es20.10)') (i - 1)*tstep, (acf(i,j)*cnv, j = 1,nsites)
enddo

close (unit=10)

end subroutine wrt_acf


!========================================================!
!                                                        !
! THIS SUBROUTINE WRITES OUT CROSS-CORRELATION FUNCTIONS !
!                                                        !
!========================================================!

subroutine wrt_ccf

open (unit=10, file='ccf.dat', status='unknown')

! set the header of the file
corr_type = '-----'
do j = 1,nsites
  do k = 1,nsites
    if (j == k) write(corr_type(j,k)(1:1),'(a1)') 'a'
    if (j /= k) write(corr_type(j,k)(1:1),'(a1)') 'x'
    write(corr_type(j,k)(3:3),'(i1)') j
    write(corr_type(j,k)(5:5),'(i1)') k
  enddo
enddo
write (10,'(a12,6x,30(a5,7x))') '#     tstep ', ( (corr_type(j,k), k = 1,nsites), j = 1,nsites)

! change output_units from eV to wavenumbers if needed
if (output_units == 'wn') then
  cnv = ev2wvn*ev2wvn
 elseif (output_units == 'ev') then
  cnv = 1.0d0
endif

! normalize ccf
!ccf0(:,:) = sign( ccf(1,:,:)/cnorm(1,:,:) , 1.0d0 )
ccf0(:,:) = 1.0d0 
do i = 1,lcor 
  write (10,'(30es20.10)') (i - 1)*tstep, ( (ccf(i,j,k)/(ccf0(j,k)*cnorm(i,j,k))*cnv, k = 1,nsites), j = 1,nsites)
enddo

close (unit=10)

end subroutine wrt_ccf


!============================================================================!
!                                                                            !
! THIS SUBROUTINE WRITES OUT SITE ENERGIES PROBABILITY DISTRIBUTION FUNCTION !
!                                                                            !
!============================================================================!

subroutine wrt_pdf

! normalize site energy probability distribution function and write them out
do i = 1,nsites
  integral = 0.0d0
  do j = - nhist, nhist - 2
      integral = integral + en_pdf(i,j) + 4.0d0*en_pdf(i,j + 1) + en_pdf(i,j + 2)
  enddo
  integral = integral*dE/3.0d0
  if (integral /= 0.0d0) then
    en_pdf(i,:) = en_pdf(i,:)/integral
   else
    en_pdf(i,:) = 0.0d0
  endif
enddo

! change output_units from eV to wavenumbers if needed
if (output_units == 'wn') then
  cnv = ev2wvn
 elseif (output_units == 'ev') then 
  cnv = 1.0d0
endif

! open output file and write results
open (unit=10, file='en_pdf.dat', status='unknown')
write (10,'(a)') '#  energy and relative pdf for each chromophore ' 
do j = - nhist + 1, nhist - 1
  ener(:) = (j - 1)*dE + av_en(:)
  write (10,'(99e16.8)') (ener(i)*cnv, en_pdf(i,j), i = 1,nsites)
enddo
close (unit=10)

end subroutine wrt_pdf


!=========================================================================!
!                                                                         !
! THIS SUBROUTINE WRITES OUT SITE ENERGIES AVERAGE AND STANDARD DEVIATION !
!                                                                         !
!=========================================================================!

subroutine wrt_av

! change output_units from eV to wavenumbers if needed
if (output_units == 'wn') then
  cnv = ev2wvn
 elseif (output_units == 'ev') then 
  cnv = 1.0d0
endif

! open output file and write results
open (unit=10, file='averages-en.out', status='unknown')
write (10,'(a)') '# site     av_en        sd_en'
do i = 1, nsites
  write (10,'(i4,2x,2f14.6)') i, av_en(i)*cnv, sd_en(i)*cnv
enddo
close (unit=10)

end subroutine wrt_av


!==================================================================================================================================!
!                                                                                                                                  !
!  THIS SUBROUTINE COMPUTES THE FAST FOURIER TRANSFORM OF THE SITE ENERGY ACF. THE OUTPUT IS THE SPECTRAL DENSITY.                 ! 
!                                                                                                                                  !
!==================================================================================================================================!

  subroutine sd 

!---------------------------------------------------- Begin Variables Declarations ------------------------------------------------!

  implicit none

  ! parameters
  real(kind=dbl), parameter :: tmax = 2.0d0, deltaw = 0.1d0
  real(kind=dbl), parameter :: light = 2.99792458d-2,      & ! speed of light in cm/ps 
                               planck = 4.135667516d-3,    & ! Planck constant in eV*ps 
                               betahbar = 7.638233025d0      ! product of h bar with inverse Boltzmann constant in K*ps

  ! declare local variables
  integer :: id, i, k, j, n1, n2

  real(kind=dbl) :: b1, b2, b3, b4, freqcnv, arg, deltat, dt, inst_time, f1, w, const, ihkt, hyp_tg  

  real(kind=dbl), allocatable, dimension(:) :: dip, bw, p, tm, f, time

!------------------------------------------------------ End Variables Declarations ------------------------------------------------!

  ! allocate arrays
  allocate ( dip(nstored) )
  allocate ( bw(nstored) )
  allocate ( p(nstored) )
  allocate ( tm(nstored) )
  allocate ( f(nstored) )
  allocate ( time(nstored) )

  ! assign values to useful variables
  freqcnv = wvn2ev/(2.0d0*pi*light)
  const = 4.0d0/planck
  ihkt = 0.5d0*betahbar/temp

  deltat = tstep 
  dt = tmax/dble(nstored - 1)

  !Definition of a Blackman window
  !b1 = 0.42D0
  !b2 =-0.50D0
  !b3 = 0.08D0
  b1 = 0.40217D0
  b2 =-0.49703D0
  b3 = 0.09392D0
  b4 =-0.00183D0
  n1 = nstored - 1
  n2 = n1 - 1
  do i = 1,nstored 
    time(i) = (i - 1)*tstep
    arg = pi*dble(i + n2)/dble(n1) 
    bw(i) = b1 + b2*dcos(arg) + b3*dcos(2.0d0*arg) + b4*dcos(3.0d0*arg)
  enddo

  ! Loop over all species
  do id = 1,nsites
    
    ! Copy autocorrelation function to a temporary array (be careful: their dimensions are different!)
    do i = 1,nstored
      dip(i) = acf(i,id)
    enddo
  
    ! If necessary, improves data point number by polynomial interpolation
    do k = 1,nstored 
      inst_time = (k - 1)*dt
      call polrea(time,dip,inst_time,f1,nstored)
      f(k) = f1
    enddo

    ! Multiplying by the Blackman window
    do i = 1,nstored 
      p(i) = f(i)*bw(i)
      tm(i) = (i - 1)*dt
    enddo
  
    ! Fourier transforming
    do j = 1,nmxpow
      w = dble(j - 1)*deltaw              ! frequency
      re_fft(j) = fourtr(p,tm,w,nstored)  ! real part of the transform
      im_fft(j) = foursi(p,tm,w,nstored)  ! immaginary part of the transform
      hyp_tg = dtanh(ihkt*w)          
      if (id == 1) freq(j) = w * freqcnv
      power(j,id) = const * hyp_tg * re_fft(j)
    enddo

  enddo

  ! change output_units from eV to wavenumbers if needed
  if (output_units == 'wn') then
    cnv = ev2wvn
   elseif (output_units == 'ev') then 
    cnv = 1.0d0
  endif

  ! write output
  open (unit=90, file='sd.dat', status='unknown')
  do j = 1,nmxpow
    write (90,'(25es20.10)') freq(j)*cnv, (power(j,id)*cnv, id = 1,nsites)
  enddo
  close (unit=90)

  ! deallocate arrays
  deallocate ( dip )
  deallocate ( bw )
  deallocate ( p )
  deallocate ( tm )
  deallocate ( f )
  deallocate ( time )

  end subroutine sd 


!==================================================================================================================================!
!                                                                                                                                  !
!  THIS FUNCTION COMPUTES THE REAL PART OF THE FAST FOURIER TRANSFORM.                                                             ! 
!                                                                                                                                  !
!==================================================================================================================================!

  real(kind=dbl) function fourtr(y,x,c,n)

!---------------------------------------------------- Begin Variables Declarations ------------------------------------------------!

  implicit none

  ! input variables
  integer, intent(in) :: n

  real(kind=dbl), intent(in) :: c
  real(kind=dbl), intent(in), dimension(n) :: x, y

  ! local variables
  integer :: i

  real(kind=dbl) :: temp1, temp2

!------------------------------------------------------ End Variables Declarations ------------------------------------------------!

  fourtr = 0.0d0
  temp1 = y(1)*dcos( x(1)*c )
  do i = 2,n
    temp2 = y(i)*dcos( x(i)*c )
    fourtr = fourtr + (temp2 + temp1)*( x(i) - x(i-1) )
    temp1 = temp2
  enddo
  fourtr = 0.5d0*fourtr
  
  end function fourtr


!==================================================================================================================================!
!                                                                                                                                  !
!  THIS FUNCTION COMPUTES THE IMAGINARY PART OF THE FAST FOURIER TRANSFORM.                                                        ! 
!                                                                                                                                  !
!==================================================================================================================================!

  real(kind=dbl) function foursi(y,x,c,n)

!---------------------------------------------------- Begin Variables Declarations ------------------------------------------------!

  implicit none

  ! input variables
  integer, intent(in) :: n

  real(kind=dbl), intent(in) :: c
  real(kind=dbl), intent(in), dimension(n) :: x, y

  ! local variables
  integer :: i

  real(kind=dbl) :: temp1, temp2

!------------------------------------------------------ End Variables Declarations ------------------------------------------------!

  foursi = 0.0d0
  temp1 = y(1)*dsin( x(1)*c )
  do i = 2,n
    temp2 = y(i)*dsin( x(i)*c )
    foursi = foursi + (temp2 + temp1)*( x(i) - x(i-1) )
    temp1 = temp2
  enddo
  foursi = 0.5d0*foursi
  
  end function foursi


!==================================================================================================================================!
!                                                                                                                                  !
!  THE FOLLOWING SUBROUTINES ARE USED TO IMPROVE THE INPUT FUNCTION THROUGH A POLYNOMIAL INTERPOLATION.                            ! 
!                                                                                                                                  !
!==================================================================================================================================!

!****************************************** Polynomial interpolation: subroutine 1 ************************************************!

  subroutine polrea(x,yr,xa,yyr,ndim)

!---------------------------------------------------- Begin Variables Declarations ------------------------------------------------!

  integer, intent(in) :: ndim

  real(kind=dbl), intent(in) :: xa, x(ndim), yr(ndim)
  real(kind=dbl), intent(out) :: yyr

  integer :: j, j1, j2, j3, n1, n2

  real(kind=dbl) :: dyr, x1(5), yr1(5) 

!------------------------------------------------------ End Variables Declarations ------------------------------------------------!
 
  call locate1(x,ndim,xa,j)

  n1 = ndim - 2
  n2 = ndim - 4
  if (j <= 2) then
    do j1 = 1,5
      x1(j1) = x(j1)
      yr1(j1) = yr(j1)
    enddo

   else if (j > n1) then
    do j1 = n2,ndim
      j3 = j1 - n2 + 1
      x1(j3) = x(j1)
      yr1(j3) = yr(j1)
    enddo

   else
    do j2 = 1,5
      x1(j2) = x(j - 2 + j2 - 1)
      yr1(j2) = yr(j - 2 + j2 - 1)
    enddo
  endif
 
  call polint1(x1,yr1,5,xa,yyr,dyr)
 
  end subroutine polrea


!****************************************** Polynomial interpolation: subroutine 2 ************************************************!

  subroutine locate1(xx,n,x,j)

!---------------------------------------------------- Begin Variables Declarations ------------------------------------------------!

  integer, intent(in) :: n
  integer, intent(out) :: j

  real(kind=dbl), intent(in) :: x, xx(n)

  integer :: jl, jm, ju

!------------------------------------------------------ End Variables Declarations ------------------------------------------------!

  jl = 0
  ju = n + 1
  do while (ju -jl > 1) 
    jm = (ju + jl)/2
    if ( (xx(n) > xx(1)) .eqv. (x > xx(jm)) ) then
      jl = jm
     else
      ju = jm
    endif
  enddo
     
  j = jl

  end subroutine locate1


!****************************************** Polynomial interpolation: subroutine 3 ************************************************!

  subroutine polint1(xa,ya,n,x,y,dy)
  
!---------------------------------------------------- Begin Variables Declarations ------------------------------------------------!

  integer, intent(in) :: n

  real(kind=dbl), intent(in) :: x, xa(n), ya(n)
  real(kind=dbl), intent(out) :: dy, y

  integer :: i, m, ns
  integer, parameter :: nmax = 100

  real(kind=dbl) :: den, dif, dift, ho, hp, w, c(nmax), d(nmax)

!------------------------------------------------------ End Variables Declarations ------------------------------------------------!

  ns = 1
  dif = dabs( x - xa(1) )
  do i = 1,n
     dift = dabs( x - xa(i) )
     if (dift < dif) then
       ns = i
       dif = dift
     endif
     c(i) = ya(i)
     d(i) = ya(i)
  enddo

  y = ya(ns)
  ns = ns - 1
  do m = 1,n - 1

    do i = 1,n - m
      ho = xa(i) - x
      hp = xa(i + m) - x
      w = c(i + 1) - d(i)
      den = ho - hp
      if (den == 0.0d0) pause 'failure in polint'
      den = w/den
      d(i) = hp*den
      c(i) = ho*den
    enddo

    if (2*ns < n - m)then
      dy = c(ns + 1)
     else
      dy = d(ns)
      ns = ns - 1
    endif

    y = y + dy

  enddo

  end subroutine polint1


!========================================!
!                                        !
! THIS SUBROUTINE DEALLOCATES ALL ARRAYS !
!                                        !
!========================================!

subroutine finish 

deallocate ( ntime0 )
deallocate ( step )
deallocate ( ener )
deallocate ( av_en )
deallocate ( sd_en )
deallocate ( ccf0 )
deallocate ( en )
deallocate ( en0 )
deallocate ( ccf )
deallocate ( cnorm )
deallocate ( en_pdf )
deallocate ( corr_type )
deallocate ( acf )
deallocate ( power )
deallocate ( freq )
deallocate ( re_fft )
deallocate ( im_fft )

end subroutine finish 


!**********************************************************************************************************************************!
!                                                                                                                                  !
!                                                   E N D       S U B R O U T I N E S                                              !
!                                                                                                                                  !
!**********************************************************************************************************************************!

END MODULE corr_pdf_av_en 


!**********************************************************************************************************************************!
!                                                                                                                                  !
!                                                  M  A  I  N     P  R  O  G  R  A  M                                              !
!                                                                                                                                  !
!**********************************************************************************************************************************!

PROGRAM analysis

USE corr_pdf_av_en

implicit none

! initial values for variables; if negative, use default (see subroutine init) 
nsites = 1                ! number of sites
nhist = 100               ! number of bins for site energies PDFs
tcor = 50.0d0             ! correlation time (ps)
tstep = 0.005d0           ! time step (ps)
dE = 2.0d-2               ! energy difference for binning site energies PDFs in eV
temp = 300.0d0            ! reference temperature
pi = 4.0d0*datan(1.0d0)   ! value of pi

! units in the input and output files
!input_units = 'ev'        
input_units = 'wn' 
!output_units = 'ev'      
output_units = 'wn'

! initialize arrays and variables
call init

! read site enrgies from results.out and compute averages
call read_data

! compute correlation functions and probability distribution functions
call compute_cf_pdf

! write auto-correlation functions
call wrt_acf

! write cross-correlation functions
call wrt_ccf

! write probability distribution functions
call wrt_pdf

! write site energies average and standard deviation
call wrt_av

! compute Spectral Density and write it to a file
call sd 

! deallocate arrays
call finish

END PROGRAM analysis

