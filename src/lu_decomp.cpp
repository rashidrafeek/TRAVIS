/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Martin Brehm.

    ---------------------------------------------------------------------------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*****************************************************************************/


#include "config.h"

#include "lu_decomp.h"
#include <math.h>
#include "tools.h"


const char *GetRevisionInfo_lu_decomp(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_lu_decomp() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


#define TINY 1.0e-20 // A small number.


void ludcmp(double **a, int n, int *indx, double *d)
/* Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
indx[1..n] is an output vector that records the row permutation effected by the partial
pivoting; d is output as +/-1 depending on whether the number of row interchanges was even
or odd, respectively. This routine is used in combination with lubksb to solve linear equations
or invert a matrix. */
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv; // vv stores the implicit scaling of each row.


	imax = -1; // Avoid compiler warning

	vv=new double[n+1];
	*d=1.0; // No row interchanges yet.
	for (i=1;i<=n;i++) { // Loop over rows to get the implicit scaling information.
		big=0.0; 
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) { eprintf("Singular matrix in routine ludcmp."); abort(); }
		// No nonzero largest element.
		vv[i]=1.0/big; // Save the scaling.
	}
	for (j=1;j<=n;j++) { // This is the loop over columns of Crout's method.
		for (i=1;i<j;i++) { // This is equation (2.3.12) except for i = j.
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0; // Initialize for the search for largest pivot element.
		for (i=j;i<=n;i++) { // This is i = j of equation (2.3.12) and i = j+1. . .N of equation (2.3.13).
			sum=a[i][j]; 
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				// Is the figure of merit for the pivot better than the best so far?
				big=dum;
				imax=i;
			}
		}
		if (j != imax) { // Do we need to interchange rows?
			for (k=1;k<=n;k++) { // Yes, do so...
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d); // ...and change the parity of d.
			vv[imax]=vv[j]; // Also interchange the scale factor.
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		// If the pivot element is zero the matrix is singular (at least to the precision of the
		// algorithm). For some applications on singular matrices, it is desirable to substitute
		// TINY for zero.
		if (j != n) { // Now, finally, divide by the pivot element.
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	} // Go back for the next column in the reduction.
	delete[] vv;
}





void lubksb(double **a, int n, int *indx, double b[])
/* Solves the set of n linear equations A x X = B. Here a[1..n][1..n] is input, not as the matrix
A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
B, and returns with the solution vector X. a, n, and indx are not modified by this routine
and can be left in place for successive calls with different right-hand sides b. This routine takes
into account the possibility that b will begin with many zero elements, so it is efficient for use
in matrix inversion. */
{
	int i,ii=0,ip,j;
	double sum;


	for (i=1;i<=n;i++) {
		// When ii is set to a positive value, it will become the
		// index of the first nonvanishing element of b. We now
		// do the forward substitution, equation (2.3.6). The
		// only new wrinkle is to unscramble the permutation
		// as we go.
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i; // A nonzero element was encountered, so from now on we will have to do the sums in the loop above.
			b[i]=sum;
	}
	for (i=n;i>=1;i--) { // Now we do the backsubstitution, equation (2.3.7).
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i]; // Store a component of the solution vector X.
	} // All done!
}

