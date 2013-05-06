! Marsaglia & Tsang generator for random normals & random exponentials.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integers are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integers in Fortran.

! Latest version - 1 January 2001

! Parallel version - October 2006

! This version has been customised for parallel processing use,
! specifically with OpenMP.  Each thread uses its own pseudo-random
! sequence. (Gib Bogle)
! Note: thread index kpar is 0-based
!--------------------------------------------------------------------------

MODULE Par_Zig_mod

   IMPLICIT NONE

   PRIVATE

   INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )
   REAL(DP), PARAMETER  ::  m1=2147483648.0_DP,   m2=2147483648.0_DP,      &
                            half=0.5_DP
   REAL(DP)             ::  dn0=3.442619855899_DP, tn0=3.442619855899_DP,    &
                            vn=0.00991256303526217_DP,                     &
                            q,                    de0=7.697117470131487_DP, &
                            te0=7.697117470131487_DP,                       &
                            ve=0.003949659822581572_DP
!   INTEGER,  SAVE       ::  iz, jz, jsr=123456789, kn(0:127),              &
!                            ke(0:255), hz
!   REAL(DP), SAVE       ::  wn(0:127), fn(0:127), we(0:255), fe(0:255)
!   LOGICAL,  SAVE       ::  initialized=.FALSE.

	integer, save :: par_n = 0, par_step
	integer, allocatable, save ::  par_jsr(:), par_kn(:,:), par_ke(:,:)
	real(DP), allocatable, save :: par_wn(:,:), par_fn(:,:), par_we(:,:), par_fe(:,:)

   PUBLIC  :: par_zigset, par_zigfree, par_shr3, par_uni, par_rnor, par_rexp, par_test

CONTAINS

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
SUBROUTINE par_zigset( npar, par_jsrseed, grainsize)

   INTEGER, INTENT(IN)  :: npar, grainsize, par_jsrseed(0:*)

   INTEGER  :: i, kpar
	REAL(DP) dn, tn, de, te

	par_n = npar
	par_step = grainsize
! First we need to allocate all the non-volatile arrays with the size npar
	allocate(par_jsr(0:npar*par_step))
	allocate(par_kn(0:127,0:npar-1))
	allocate(par_ke(0:255,0:npar-1))
	allocate(par_wn(0:127,0:npar-1))
	allocate(par_fn(0:127,0:npar-1))
	allocate(par_we(0:255,0:npar-1))
	allocate(par_fe(0:255,0:npar-1))

! Now treat each instance separately
do kpar = 0,npar-1
   !  Set the seed
   par_jsr(kpar*par_step) = par_jsrseed(kpar)

   !  Tables for RNOR
   dn = dn0
   tn = tn0
   q = vn*EXP(half*dn*dn)
   par_kn(0,kpar) = (dn/q)*m1
   par_kn(1,kpar) = 0
   par_wn(0,kpar) = q/m1
   par_wn(127,kpar) = dn/m1
   par_fn(0,kpar) = 1.0_DP
   par_fn(127,kpar) = EXP( -half*dn*dn )
   DO  i = 126, 1, -1
      dn = SQRT( -2.0_DP * LOG( vn/dn + EXP( -half*dn*dn ) ) )		! dn
      par_kn(i+1,kpar) = (dn/tn)*m1
      tn = dn														! tn
      par_fn(i,kpar) = EXP(-half*dn*dn)
      par_wn(i,kpar) = dn/m1
   END DO

   !  Tables for REXP
   de = de0
   te = te0
   q = ve*EXP( de )
   par_ke(0,kpar) = (de/q)*m2
   par_ke(1,kpar) = 0
   par_we(0,kpar) = q/m2
   par_we(255,kpar) = de/m2
   par_fe(0,kpar) = 1.0_DP
   par_fe(255,kpar) = EXP( -de )
   DO  i = 254, 1, -1
      de = -LOG( ve/de + EXP( -de ) )								! de
      par_ke(i+1,kpar) = m2 * (de/te)
      te = de														! te
      par_fe(i,kpar) = EXP( -de )
      par_we(i,kpar) = de/m2
   END DO
enddo
RETURN
END SUBROUTINE par_zigset

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine par_zigfree
deallocate(par_jsr)
deallocate(par_kn)
deallocate(par_ke)
deallocate(par_wn)
deallocate(par_fn)
deallocate(par_we)
deallocate(par_fe)
end subroutine

!--------------------------------------------------------------------------
!  Generate random 32-bit integers
!--------------------------------------------------------------------------
FUNCTION par_shr3(kpar) RESULT( ival )
   INTEGER  ::  ival, kpar
   integer :: jz, jsr

   jsr = par_jsr(kpar*par_step)
   jz = jsr
   jsr = IEOR( jsr, ISHFT( jsr,  13 ) )
   jsr = IEOR( jsr, ISHFT( jsr, -17 ) )
   jsr = IEOR( jsr, ISHFT( jsr,   5 ) )
   par_jsr(kpar*par_step) = jsr
   ival = jz + jsr
   RETURN
END FUNCTION par_shr3

!--------------------------------------------------------------------------
!  Generate uniformly distributed random numbers, sequence kpar
!--------------------------------------------------------------------------
FUNCTION par_uni(kpar) RESULT( fn_val )
	integer :: kpar
   REAL(DP)  ::  fn_val

	if (kpar >= par_n) then
		write(*,*) 'thread number: ',kpar,' exceeds initialized max: ',par_n-1
		write(21,*) 'thread number: ',kpar,' exceeds initialized max: ',par_n-1
		stop
	endif
   fn_val = half + 0.2328306e-9_DP * par_shr3(kpar)
   RETURN
END FUNCTION par_uni

FUNCTION par_test(kpar) RESULT( fn_val )
	integer :: kpar
	integer :: fn_val
	if (kpar >= par_n) then
		write(*,*) 'thread number: ',kpar,' exceeds initialized max: ',par_n-1
		stop
	endif
    fn_val = kpar
    return
end function

!--------------------------------------------------------------------------
!  Generate random normals, sequence kpar
!--------------------------------------------------------------------------
FUNCTION par_rnor(kpar) RESULT( fn_val )
   REAL(DP)             ::  fn_val
   integer :: kpar

   REAL(DP), PARAMETER  ::  r = 3.442620_DP
   REAL(DP)             ::  x, y
   integer :: iz, hz

!   IF( .NOT. initialized ) CALL zigset( jsr )
	if (kpar >= par_n) then
		write(*,*) 'thread number exceeds initialized max: ',kpar,par_n-1
		stop
	endif
   hz = par_shr3(kpar)
   iz = IAND( hz, 127 )
   IF( ABS( hz ) < par_kn(iz,kpar) ) THEN
      fn_val = hz * par_wn(iz,kpar)
   ELSE
      DO
         IF( iz == 0 ) THEN
            DO
               x = -0.2904764_DP* LOG( par_uni(kpar) )
               y = -LOG( par_uni(kpar) )
               IF( y+y >= x*x ) EXIT
            END DO
            fn_val = r+x
            IF( hz <= 0 ) fn_val = -fn_val
            RETURN
         END IF
         x = hz * par_wn(iz,kpar)
         IF( par_fn(iz,kpar) + par_uni(kpar)*(par_fn(iz-1,kpar)-par_fn(iz,kpar)) < EXP(-half*x*x) ) THEN
            fn_val = x
            RETURN
         END IF
         hz = par_shr3(kpar)
         iz = IAND( hz, 127 )
         IF( ABS( hz ) < par_kn(iz,kpar) ) THEN
            fn_val = hz * par_wn(iz,kpar)
            RETURN
         END IF
      END DO
   END IF
   RETURN
END FUNCTION par_rnor

!--------------------------------------------------------------------------
!  Generate random exponentials, sequence kpar
!--------------------------------------------------------------------------
FUNCTION par_rexp(kpar) RESULT( fn_val )
   REAL(DP)  ::  fn_val
   integer :: kpar

   REAL(DP)  ::  x
   integer :: iz, jz

!   IF( .NOT. initialized ) CALL Zigset( jsr )
	if (kpar >= par_n) then
		write(*,*) 'thread number exceeds initialized max: ',kpar,par_n-1
		stop
	endif
   jz = par_shr3(kpar)
   iz = IAND( jz, 255 )
   IF( ABS( jz ) < par_ke(iz,kpar) ) THEN
      fn_val = ABS(jz) * par_we(iz,kpar)
      RETURN
   END IF
   DO
      IF( iz == 0 ) THEN
         fn_val = 7.69711 - LOG(par_uni(kpar) )
         RETURN
      END IF
      x = ABS( jz ) * par_we(iz,kpar)
      IF( par_fe(iz,kpar) + par_uni(kpar)*(par_fe(iz-1,kpar) - par_fe(iz,kpar)) < EXP( -x ) ) THEN
         fn_val = x
         RETURN
      END IF
      jz = par_shr3(kpar)
      iz = IAND( jz, 255 )
      IF( ABS( jz ) < par_ke(iz,kpar) ) THEN
         fn_val = ABS( jz ) * par_we(iz,kpar)
         RETURN
      END IF
   END DO
   RETURN
END FUNCTION par_rexp

END MODULE par_zig_mod



