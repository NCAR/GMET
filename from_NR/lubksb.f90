	SUBROUTINE lubksb(a,indx,b)
	! http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/lu_f90.txt
	!  ******************************************************************
	!  * Solves the set of N linear equations A . X = B.  Here A is     *
	!  * input, not as the matrix A but rather as its LU decomposition, *
	!  * determined by the routine LUDCMP. INDX is input as the permuta-*
	!  * tion vector returned by LUDCMP. B is input as the right-hand   *
	!  * side vector B, and returns with the solution vector X. A, N and*
	!  * INDX are not modified by this routine and can be used for suc- *
	!  * cessive calls with different right-hand sides. This routine is *
	!  * also efficient for plain matrix inversion.                     *
	!  ******************************************************************
	USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
	INTEGER(I4B) :: i,n,ii,ll
	REAL(SP) :: summ
	n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
	ii=0
	do i=1,n
		ll=indx(i)
		summ=b(ll)
		b(ll)=b(i)
		if (ii /= 0) then
			summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
		else if (summ /= 0.0) then
			ii=i
		end if
		b(i)=summ
	end do
	do i=n,1,-1
		b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
	end do
	END SUBROUTINE lubksb
