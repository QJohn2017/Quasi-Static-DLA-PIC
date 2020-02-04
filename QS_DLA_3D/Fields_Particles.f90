module global_variables
integer, parameter ::RP=8
real(RP), parameter :: pi=3.14159265358979d0
complex(RP), parameter :: ci=(0d0,1d0)
end module global_variables

module output_var
	USE global_variables
	real(RP),allocatable :: r2d(:,:),V2d(:,:),force2d(:,:),gama2d(:,:)
 	real(RP),allocatable :: psi2d(:,:),Er2d(:,:),Et2d(:,:),Bth2d(:,:)
  real(RP),allocatable :: den2d(:,:),den_e2d(:,:),jr2d(:,:),jz2d(:,:)
end module output_var

    module pulse_var
    use global_variables
    implicit none

    integer :: Types
    real(RP)::a0,Wist0,k0,FocalPlane,zcenter,FWHM
    complex(RP),allocatable ::Acold2d(:,:),Atemp2d(:,:),Ac2d(:,:)
    end module pulse_var


module TestPart_var
USE global_variables
integer :: NTP,NTP_i(5),TP_perZ(5),TP_perR(5),N_per_peri(5),N_per_dt,Tp_Dt_per_out,pushtype,Num_of_Beam   !NTP=TP_len*TP_perZ+TP_Rad*TP_perR
real(RP) :: TP_len(5),TP_Rad(5),TP_z0(5),GamaT(5)
logical :: TestParts,LaserEZ
real(RP),allocatable:: Err2d(:,:),Ezz2d(:,:),Bthh2d(:,:),&
Err2dold(:,:),Ezz2dold(:,:),Bthh2dold(:,:)
complex(RP),allocatable:: EL2d(:,:),BL2d(:,:),EL2dold(:,:),BL2dold(:,:),ELZ2dold(:,:),ELZ2d(:,:)
real(RP),allocatable:: XXT(:),YYT(:),ZZT(:),PXT(:),PYT(:),PZT(:)
real(RP),allocatable:: WX_L(:),WZ_L(:),WX_W(:),WY_W(:),WZ_W(:)
real(RP),allocatable:: XXT2d(:,:),YYT2d(:,:),ZZT2d(:,:),PXT2d(:,:),PYT2d(:,:),PZT2d(:,:)

! (1) Err full-z full-r , (2)Ez full-z and half-r. (3) Bth full-z full-r 
! (4) ELaser  half-z half-r (5) BL, half-z, half-r
! (4) EZLaser  Full-z Full-r

end module TestPart_var


module particles_var
USE global_variables
	real(RP),allocatable :: r(:),V(:),force(:),gama(:),mass(:)
end module particles_var


module fields_var
	USE global_variables
    real(RP),allocatable :: psi(:),Er(:),Et(:),Bth(:)
    real(RP),allocatable :: den(:),den_e(:),jr(:),jv(:),jz(:),Chi2d(:,:)
end module fields_var

module domain
	USE global_variables
    real(RP):: Rmax,Zmax,Tmax
    real(RP):: hr  ! meshsize
    real(RP):: dr,dz0,dt0
    integer :: Nr,Nz  ! grid number
    integer :: Part_per_cell,Push_Order,Np,NWake !particle number 
    real(RP):: Push_per_cell
    integer :: Dt_per_out
    real(RP):: RayLen
    real(Rp):: WakeLen
    real(Rp):: ramp_len,rampBegin
end module domain


  module InAndOut
  use pulse_var
  use domain
  use netcdf
  use TestPart_var
    implicit none 
    integer :: &
    ncidf,ncidp,ncidtp,ncids,ztime_dim, p_dim, &
    r_dim,dims(2),var_de,var_jr,var_dr,&
    var_bt,var_ps,var_jz,var_ez,var_er,&
    var_vr,var_rr,var_fr,var_ri,var_rf,&
    var_ga,var_AR,var_AI,var_BR,var_BI,var_ELR,var_ELI,&
    var_Tx,var_Ty,var_Tz,var_Tpx,var_Tpy,var_Tpz,&
    var_Wxl,var_Wxw,var_Wyw,var_Wzw,var_EZLR,var_EZLI,var_Wzl

    NAMELIST/Lasers/ a0,Wist0,k0,FocalPlane,zcenter,FWHM,Types
    NAMELIST/Simulation/ Rmax,Zmax,Nr,Nz,Part_per_cell,Push_per_cell,Push_Order,dt0,Tmax,Dt_per_out,ramp_len,rampBegin
    NAMELIST/TestParticles/ TestParts,LaserEZ,Num_of_Beam,&
    TP_len,TP_Rad,TP_z0,GamaT,TP_perZ,TP_perR,N_per_peri,N_per_dt,Tp_Dt_per_out,pushtype

    contains

    subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
    print *, trim(nf90_strerror(status))
    stop 2
    end if
    end subroutine check  
end module InAndOut


module Functions
USE global_variables
use domain
implicit none
contains



 SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!
!  -- LAPACK routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )
!     ..
!
!  Purpose =======
!  DGTSV  solves the equation
!     A*X = B,
!  where A is an n by n tridiagonal matrix, by Gaussian elimination with
!  partial pivoting.
! Note that the equation  A**T*X = B  may be solved by interchanging the
!  order of the arguments DU and DL.
!
!  Arguments
!  =========
! N       (input) INTEGER
!          The order of the matrix A.  N >= 0.!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, DL must contain the (n-1) sub-diagonal elements of
!          A.
!          On exit, DL is overwritten by the (n-2) elements of the
!          second super-diagonal of the upper triangular matrix U from
!          the LU factorization of A, in DL(1), ..., DL(n-2).!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, D must contain the diagonal elements of A.
!
!         On exit, D is overwritten by the n diagonal elements of U.
!
!  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, DU must contain the (n-1) super-diagonal elements
!          of A.
!
!          On exit, DU is overwritten by the (n-1) elements of the first
!          super-diagonal of U.
!
!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N by NRHS matrix of right hand side matrix B.
!          On exit, if INFO = 0, the N by NRHS solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = i, U(i,i) is exactly zero, and the solution
!               has not been computed.  The factorization has not been
!               completed unless i = N.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
!     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   FACT, TEMP
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Executable Statements ..
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
       !  CALL XERBLA( 'DGTSV ', -INFO )
         RETURN
      END IF
!
      IF( N.EQ.0 )   RETURN
!
      IF( NRHS.EQ.1 ) THEN
         DO 10 I = 1, N - 2
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
!
!              No row interchange required
!
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               END IF
               DL( I ) = ZERO
            ELSE
!
!              Interchange rows I and I+1
!
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            END IF
   10    CONTINUE
         IF( N.GT.1 ) THEN
            I = N - 1
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               END IF
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            END IF
         END IF
         IF( D( N ).EQ.ZERO ) THEN
            INFO = N
            RETURN
         END IF
      ELSE
         DO 40 I = 1, N - 2
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
!
!              No row interchange required
!
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 20 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   20             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               END IF
               DL( I ) = ZERO
            ELSE
!
!              Interchange rows I and I+1
!
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               DO 30 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   30          CONTINUE
            END IF
   40    CONTINUE
         IF( N.GT.1 ) THEN
            I = N - 1
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 50 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   50             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               END IF
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               DO 60 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   60          CONTINUE
            END IF
         END IF
         IF( D( N ).EQ.ZERO ) THEN
            INFO = N
            RETURN
         END IF
      END IF
!     Back solve with the matrix U from the factorization.
!
      IF( NRHS.LE.2 ) THEN
         J = 1
   70    CONTINUE
         B( N, J ) = B( N, J ) / D( N )
         IF( N.GT.1 )      B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
         DO 80 I = N - 2, 1, -1
            B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )*B( I+2, J ) ) / D( I )
   80    CONTINUE
         IF( J.LT.NRHS ) THEN
            J = J + 1
            GO TO 70
         END IF
      ELSE
         DO 100 J = 1, NRHS
            B( N, J ) = B( N, J ) / D( N )
            IF( N.GT.1 ) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) /D( N-1 )
            DO 90 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )*B( I+2, J ) ) / D( I )
   90       CONTINUE
  100    CONTINUE
      END IF
!
      RETURN
!
!     End of DGTSV
!
      END


    SUBROUTINE SGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      REAL(RP)::               B( LDB, * ), D( * ), DL( * ), DU( * )
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL(RP)::               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
      INTEGER            I, J
      REAL(RP)::               FACT, TEMP
      INTRINSIC          ABS, MAX
      EXTERNAL           XERBLA
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
       !  CALL XERBLA( 'SGTSV ', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 )   RETURN
      IF( NRHS.EQ.1 ) THEN
         DO 10 I = 1, N - 2
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               END IF
               DL( I ) = ZERO
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            END IF
   10    CONTINUE
         IF( N.GT.1 ) THEN
            I = N - 1
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  B( I+1, 1 ) = B( I+1, 1 ) - FACT*B( I, 1 )
               ELSE
                  INFO = I
                  RETURN
               END IF
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               TEMP = B( I, 1 )
               B( I, 1 ) = B( I+1, 1 )
               B( I+1, 1 ) = TEMP - FACT*B( I+1, 1 )
            END IF
         END IF
         IF( D( N ).EQ.ZERO ) THEN
            INFO = N
            RETURN
         END IF
      ELSE
         DO 40 I = 1, N - 2
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 20 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   20             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               END IF
               DL( I ) = ZERO
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DL( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DL( I )
               DU( I ) = TEMP
               DO 30 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   30          CONTINUE
            END IF
   40    CONTINUE
         IF( N.GT.1 ) THEN
            I = N - 1
            IF( ABS( D( I ) ).GE.ABS( DL( I ) ) ) THEN
               IF( D( I ).NE.ZERO ) THEN
                  FACT = DL( I ) / D( I )
                  D( I+1 ) = D( I+1 ) - FACT*DU( I )
                  DO 50 J = 1, NRHS
                     B( I+1, J ) = B( I+1, J ) - FACT*B( I, J )
   50             CONTINUE
               ELSE
                  INFO = I
                  RETURN
               END IF
            ELSE
               FACT = D( I ) / DL( I )
               D( I ) = DL( I )
               TEMP = D( I+1 )
               D( I+1 ) = DU( I ) - FACT*TEMP
               DU( I ) = TEMP
               DO 60 J = 1, NRHS
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - FACT*B( I+1, J )
   60          CONTINUE
            END IF
         END IF
         IF( D( N ).EQ.ZERO ) THEN
            INFO = N
            RETURN
         END IF
      END IF
      IF( NRHS.LE.2 ) THEN
         J = 1
   70    CONTINUE
         B( N, J ) = B( N, J ) / D( N )
         IF( N.GT.1 )     B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
         DO 80 I = N - 2, 1, -1
            B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )*B( I+2, J ) ) / D( I )
   80    CONTINUE
         IF( J.LT.NRHS ) THEN
            J = J + 1
            GO TO 70
         END IF
      ELSE
         DO 100 J = 1, NRHS
            B( N, J ) = B( N, J ) / D( N )
            IF( N.GT.1 ) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) /D( N-1 )
            DO 90 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DL( I )*B( I+2, J ) ) / D( I )
   90       CONTINUE
  100    CONTINUE
      END IF
      RETURN
!
!     End of SGTSV
!
      END



!  =====================================================================
      SUBROUTINE CGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
      INTEGER            INFO, LDB, N, NRHS
      COMPLEX(RP)::            B( LDB, * ), D( * ), DL( * ), DU( * )
      COMPLEX(RP)::            ZERO
      PARAMETER          ( ZERO = ( 0.0E+0, 0.0E+0 ) )
      INTEGER            J, K
      COMPLEX(RP)::            MULT, TEMP, ZDUM
      INTRINSIC          ABS, AIMAG, MAX, REAL
      EXTERNAL           XERBLA
      REAL(RP)::               CABS1
!     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
        ! CALL XERBLA( 'CGTSV ', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 )   RETURN
      DO 30 K = 1, N - 1
         IF( DL( K ).EQ.ZERO ) THEN
            IF( D( K ).EQ.ZERO ) THEN
               INFO = K
               RETURN
            END IF
         ELSE IF( CABS1( D( K ) ).GE.CABS1( DL( K ) ) ) THEN
            MULT = DL( K ) / D( K )
            D( K+1 ) = D( K+1 ) - MULT*DU( K )
            DO 10 J = 1, NRHS
               B( K+1, J ) = B( K+1, J ) - MULT*B( K, J )
   10       CONTINUE
            IF( K.LT.( N-1 ) )         DL( K ) = ZERO
         ELSE
            MULT = D( K ) / DL( K )
            D( K ) = DL( K )
            TEMP = D( K+1 )
            D( K+1 ) = DU( K ) - MULT*TEMP
            IF( K.LT.( N-1 ) ) THEN
               DL( K ) = DU( K+1 )
               DU( K+1 ) = -MULT*DL( K )
            END IF
            DU( K ) = TEMP
            DO 20 J = 1, NRHS
               TEMP = B( K, J )
               B( K, J ) = B( K+1, J )
               B( K+1, J ) = TEMP - MULT*B( K+1, J )
   20       CONTINUE
         END IF
   30 CONTINUE
      IF( D( N ).EQ.ZERO ) THEN
         INFO = N
         RETURN
      END IF
      DO 50 J = 1, NRHS
         B( N, J ) = B( N, J ) / D( N )
         IF( N.GT.1 ) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
         DO 40 K = N - 2, 1, -1
            B( K, J ) = ( B( K, J )-DU( K )*B( K+1, J )-DL( K )*B( K+2, J ) ) / D( K )
   40    CONTINUE
   50 CONTINUE
      RETURN
      END



    subroutine solve_tridiag(a,b,c,d,x,n)
    
!      |b1 c1       0    | |x1|  |d1|
!      |a2 b2 c2         | |x2|  |d2| 
!      |   a3 b3 .       | |. |= |. |
!      |      .  .  c_n-1| |. |  |. |
!      |0        an bn   | |xn|  |dn|
    
    implicit none
!	 a - sub-diagonal (means it is the diagonal below the main diagonal)
!	 b - the main diagonal
!	 c - sup-diagonal (means it is the diagonal above the main diagonal)
!	 d - right part
!	 x - the answer
!	 n - number of equations
   integer,intent(in) :: n
        real(8),dimension(n),intent(in) :: a,b,c,d
        real(8),dimension(n),intent(out) :: x
        real(8),dimension(n) :: cp,dp
        real(8) :: m
        integer i

! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         enddo
! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do

    end subroutine solve_tridiag
    

!     ************************************************
       real(RP) function FastestInDomain(V,r)
       implicit none
       real(RP) :: V(1:Np),r(1:Np)
       integer :: n
       FastestInDomain=0d0
       do n=1,Np
       if (V(n)>=FastestInDomain.and.r(n)<Rmax-hr/2) FastestInDomain=V(n)
       end do
       end function FastestInDomain



        subroutine calck1( arg, result1, jint )
        ! at jint=1 subroutine calculates K1(arg,result)
        !*****************************************************************************80
        !
        !! CALCK1 computes variouss K1 Bessel functions.
        !
        !  Discussion:
        !
        !    This routine computes modified Bessel functions of the second kind
        !    and order one, K1(X) and EXP(X)*K1(X), for real arguments X.
        !
        !    The main computation evaluates slightly modified forms of near
        !    minimax rational approximations generated by Russon and Blair,
        !    Chalk River (Atomic Energy of Canada Limited) Report AECL-3461,
        !    1969.
        ! 
        !  Licensing:
        !
        !    This code is distributed under the GNU LGPL license.
        !
        !  Modified:
        !
        !    03 April 2007
        !
        !  Author:
        !
        !    Original FORTRAN77 version by William Cody, Laura Stoltz.
        !    FORTRAN90 version by John Burkardt.
        !
        !  Parameters:
        !
        !    Input, real ( kind = 8 ) ARG, the argument.  XLEAST < ARG is
        !    always required.  If JINT = 1, then the argument must also be
        !    less than XMAX.
        !
        !    Output, real ( kind = 8 ) RESULT, the value of the function,
        !    which depends on the input value of JINT:
        !    1, RESULT = K1(x);
        !    2, RESULT = exp(x) * K1(x);
        !
        !    Input, integer ( kind = 4 ) JINT, chooses the function to be computed.
        !    1, K1(x);
        !    2, exp(x) * K1(x);
        implicit none
        real(kind=8):: arg
        real(kind=8):: f(5)
        real(kind=8):: g(3)
        integer:: i
        integer:: jint
        real(kind=8):: p(5)
        real(kind=8):: pp(11)
        real(kind=8):: q(3)
        real(kind=8):: qq(9)
        real(kind=8):: result1
        real(kind=8):: sumf
        real(kind=8):: sumg
        real(kind=8):: sump
        real(kind=8):: sumq
        real(kind=8):: x
        real(kind=8):: xinf
        real(kind=8):: xmax
        real(kind=8):: xleast
        real(kind=8):: xsmall
        real(kind=8):: xx
        !  Machine-dependent constants
        data xleast /2.23d-308/
        data xsmall /1.11d-16/
        data xinf /1.79d+308/
        data xmax /705.343d+0/
        !  Coefficients for  XLEAST <=  ARG  <= 1.0
        data   p/ 4.8127070456878442310d-1, 9.9991373567429309922d+1,&
           &  7.1885382604084798576d+3, 1.7733324035147015630d+5,7.1938920065420586101d+5/
        data   q/-2.8143915754538725829d+2, 3.7264298672067697862d+4, -2.2149374878243304548d+6/
        data   f/-2.2795590826955002390d-1,-5.3103913335180275253d+1,-4.5051623763436087023d+3&
           &,-1.4758069205414222471d+5, -1.3531161492785421328d+6/
        data   g/-3.0507151578787595807d+2, 4.3117653211351080007d+4,-2.7062322985570842656d+6/
        !
        !  Coefficients for  1.0 < ARG
        data  pp/ 6.4257745859173138767d-2, 7.5584584631176030810d+0,&
        &               1.3182609918569941308d+2, 8.1094256146537402173d+2,&
        &               2.3123742209168871550d+3, 3.4540675585544584407d+3,&
        &               2.8590657697910288226d+3, 1.3319486433183221990d+3,&
        &               3.4122953486801312910d+2, 4.4137176114230414036d+1,&
        &               2.2196792496874548962d+0/
        data  qq/ 3.6001069306861518855d+1, 3.3031020088765390854d+2,&
        &               1.2082692316002348638d+3, 2.1181000487171943810d+3,& 
        &               1.9448440788918006154d+3, 9.6929165726802648634d+2,&
        &               2.5951223655579051357d+2, 3.4552228452758912848d+1,&
        &               1.7710478032601086579d+0/
     
        x = arg
        !  Error return for ARG < XLEAST.
         if ( x < xleast ) then  
        result1 = xinf
         !  XLEAST <= ARG <= 1.0.
        else if ( x <= 1.0D+00 ) then     
        if ( x < xsmall ) then
        !
        !  Return for small ARG.
        !
            result1 = 1.0D+00 / x
                   else       
            xx = x * x    
            sump = ((((p(1)* xx + p(2) )* xx + p(3) )* xx + p(4) )* xx + p(5) )* xx + q(3)
            sumq = (( xx + q(1) ) * xx + q(2) )  * xx + q(3)
            sumf = ((( f(1) * xx + f(2) )  * xx + f(3) ) * xx + f(4) ) * xx + f(5)
            sumg = (( xx + g(1) ) * xx + g(2) ) * xx + g(3)
            result1 = ( xx * log ( x ) * sumf / sumg + sump / sumq ) / x
            if ( jint == 2 ) then
                result1 = result1 * exp ( x )
            end if   
        end if  
         else if ( jint == 1 .AND. xmax < x ) then
            !
        !  Error return for XMAX < ARG.
        result1 = 0.0D+00 
        else
        !
        !  1.0 < ARG.
        xx = 1.0D+00 / x  
        sump = pp(1)
        do i = 2, 11
            sump = sump * xx + pp(i)
        end do     
        sumq = xx
        do i = 1, 8
            sumq = ( sumq + qq(i) ) * xx
        end do
        sumq = sumq + qq(9)     
        result1 = sump / sumq / sqrt ( x )    
        if ( jint == 1 ) then
            result1 = result1 * exp ( -x )
        end if       
        end if
        return
        end subroutine calck1


end module Functions

