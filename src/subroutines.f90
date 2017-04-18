!---------------------------------------------------
! Subroutines.f90 are a collection of very general subroutines,
!  functions, and wrappers used by several objects.
!
! Author: Daniel Montemayor
!---------------------------------------------------
subroutine diagonalize(N,A,W)
  implicit none

  integer,  intent(in)::N
  character,parameter::JOBZ='V',UPLO='U'
  complex*16:: WORK(2*N)
  real*8::RWORK(3*N-2)
  integer::INFO,i,j

  complex*16,intent(inout):: A(N,N)
  real*8,intent(out)::W(N)

  complex*16::B(N,N)

  write(*,*)'Diagonalize subroutine has been disabled...enable MKL package'
!!$  !solve problem H*psi=E*psi using intel MKL package
!!$
!!$  call zheev(JOBZ, UPLO, N, A, N&
!!$       , W, WORK, 2*N, RWORK, INFO)
!!$
!!$  do i=1,N
!!$     B(i,1)=W(N-i+1)
!!$  end do
!!$  W=B(:,1)
!!$
!!$  do i=1,N
!!$     do j=1,N
!!$        B(i,j)=A(i,N-j+1)
!!$     end do
!!$  end do
!!$  A=B

end subroutine diagonalize
!---------------------------------------------------
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

SUBROUTINE PRAXISFIT(param,n,verbosity)
  use type_kinds
  use string

  IMPLICIT NONE
  integer(long),intent(in)::n
  real(double),intent(inout)::param(n)
  integer(long),intent(in),optional::verbosity
  real(double)::t0,machep,h0,f,aux
  real(double)::praxis
  integer(long):: prin,i,idum
  external f,praxis
  
  t0=1.D-05
  machep=epsilon(param)
  h0=10.0
  prin=0
  if(present(verbosity))prin=verbosity
  aux=praxis(t0,machep,h0,n,prin,param,f,idum)
  !write(*,FMT='('//trim(int2str(n))//'(F16.8,1X))') param
  
END SUBROUTINE PRAXISFIT
