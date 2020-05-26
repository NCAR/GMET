	SUBROUTINE ludcmp(a,indx,d)
	! http://jean-pierre.moreau.pagesperso-orange.fr/Fortran/lu_f90.txt
	!  ***************************************************************
	!  * Given an N x N matrix A, this routine replaces it by the LU *
	!  * decomposition of a rowwise permutation of itself. A and N   *
	!  * are input. INDX is an output vector which records the row   *
	!  * permutation effected by the partial pivoting; D is output   *
	!  * as -1 or 1, depending on whether the number of row inter-   *
	!  * changes was even or odd, respectively. This routine is used *
	!  * in combination with LUBKSB to solve linear equations or to  *
	!  * invert a matrix. Return code is 1, if matrix is singular.   *
	!  ***************************************************************
	USE nrtype; USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
	REAL(SP), INTENT(OUT) :: d
	REAL(SP), DIMENSION(size(a,1)) :: vv
	REAL(SP), PARAMETER :: TINY=1.0e-20_sp
	INTEGER(I4B) :: j,n,imax
	! determine whether the three numbers are equal. if not, stop
	n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
	d=1.0
	vv=maxval(abs(a),dim=2)
	if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
	vv=1.0_sp/vv
	do j=1,n
		imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
		if (j /= imax) then
			call swap(a(imax,:),a(j,:))
			d=-d
			vv(imax)=vv(j)
		end if
		indx(j)=imax
		if (a(j,j) == 0.0) a(j,j)=TINY
		a(j+1:n,j)=a(j+1:n,j)/a(j,j)
		a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
	end do
	END SUBROUTINE ludcmp
