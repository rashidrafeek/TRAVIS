/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

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


// This must always be the first include directive
#include "config.h"

#include "diagonal.h"
#include <math.h>


const char *GetRevisionInfo_diagonal(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_diagonal() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


bool Diagonalize(const CxDMatrixMN &inmat, CxDMatrixMN *evec, CxDoubleArray *eval, bool verbose) {

	UNUSED(verbose);

	int n, z, z2;
	double **a, *d, *e;

	if (inmat.GetRows() != inmat.GetCols()) {
		eprintf("Diagonalize(): Error: rows != cols (%d, %d).\n",inmat.GetRows(),inmat.GetCols());
		return false;
	}

	n = inmat.GetRows();

	d = new double[n+1];
	e = new double[n+1];
	a = new double*[n+1];

	for (z=0;z<n+1;z++)
		a[z] = new double[n+1];

	for (z=0;z<n;z++)
		for (z2=0;z2<n;z2++)
			a[z+1][z2+1] = inmat(z,z2);

	tred2(a,n,d,e,verbose);	

	tqli(d,e,n,a,verbose);

	if (evec != NULL) {
		*evec = CxDMatrixMN(n,n);
		for (z=0;z<n;z++)
			for (z2=0;z2<n;z2++)
				(*evec)(z,z2) = a[z+1][z2+1];
	}

	if (eval != NULL) {
		eval->SetSize(n);
		for (z=0;z<n;z++)
			(*eval)[z] = d[z+1];
	}

	for (z=0;z<n+1;z++)
		delete[] a[z];
	delete[] a;
	delete[] e;
	delete[] d;

	return true;
}


static double sqrarg;

#define NRSQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


double pythag(double a, double b)
// Computes (a^2 + b^2)^1/2 without destructive underflow or overflow.
{
	double absa,absb;

	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb)
		return absa*sqrt(1.0+NRSQR(absb/absa));
	else
		return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+NRSQR(absa/absb)));
}


void tred2(double **a, int n, double d[], double e[], bool verbose)
/* Householder reduction of a real, symmetric matrix a[1..n][1..n]. On output, a is replaced
by the orthogonal matrix Q effecting the transformation. d[1..n] returns the diagonal elements
of the tridiagonal matrix, and e[1..n] the off-diagonal elements, with e[1]=0. Several
statements, as noted in comments, can be omitted if only eigenvalues are to be found, in which
case a contains no useful information on output. Otherwise they are to be included. */
{
	int l,k,j,i;
	double scale,hh,h,g,f;
	double tfs;

	if (verbose) {
		mprintf("      Householder  ");
		mprintf(WHITE,"[");
	}

	tfs = n/60.0;

	for (i=n;i>=2;i--) {
		if (verbose)
			if (fmod(i,tfs) < 1.0)
				mprintf(WHITE,"#");
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0) // Skip transformation.
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale; // Use scaled a's for transformation.
					h += a[i][k]*a[i][k]; // Form sigma in h.
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g; // Now h is equation (11.2.4).
				a[i][l]=f-g; // Store u in the ith row of a.
				f=0.0;
				for (j=1;j<=l;j++) {
					/* Next statement can be omitted if eigenvectors not wanted */
					a[j][i]=a[i][j]/h; // Store u/H in ith column of a.
					g=0.0; // Form an element of A x u in g.
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h; //Form element of p in temporarily unused element of e.
					f += e[j]*a[i][j];
				}
				hh=f/(h+h); // Form K, equation (11.2.11).
				for (j=1;j<=l;j++) { // Form q and store in e overwriting p.
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++) // Reduce a, equation (11.2.13).
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
	}
	if (verbose)
		mprintf(WHITE,"]\n");
	/* Next statement can be omitted if eigenvectors not wanted */
	d[1]=0.0;
	e[1]=0.0;
	if (verbose) {
		mprintf("      Eigenvectors ");
		mprintf(WHITE,"[");
	}
	/* Contents of this loop can be omitted if eigenvectors not
	wanted except for statement d[i]=a[i][i]; */
	for (i=1;i<=n;i++) { // Begin accumulation of transformation matrices.
		if (verbose)
			if (fmod(i,tfs) < 1.0)
				mprintf(WHITE,"#");
		l=i-1; 
		if (d[i]) { // This block skipped when i=1.
			for (j=1;j<=l;j++) {
				g=0.0;
				for (k=1;k<=l;k++) // Use u and u/H stored in a to form P x Q.
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i]; // This statement remains.
		a[i][i]=1.0; // Reset row and column of a to identity matrix for next iteration.
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0; 
	}
	if (verbose)
		mprintf(WHITE,"]\n");
}




void tqli(double d[], double e[], int n, double **z, bool verbose)
/* QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real, symmetric,
tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2 11.2. On
input, d[1..n] contains the diagonal elements of the tridiagonal matrix. On output, it returns
the eigenvalues. The vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrix,
with e[1] arbitrary. On output e is destroyed. When finding only the eigenvalues, several lines
may be omitted, as noted in the comments. If the eigenvectors of a tridiagonal matrix are desired,
the matrix z[1..n][1..n] is input as the identity matrix. If the eigenvectors of a matrix
that has been reduced by tred2 are required, then z is input as the matrix output by tred2.
In either case, the kth column of z returns the normalized eigenvector corresponding to d[k]. */
{
	double pythag(double a, double b);
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;
	double tfs;

	if (verbose) { 
		mprintf("      QL           ");
		mprintf(WHITE,"[");
	}

	tfs = n/60.0;

	for (i=2;i<=n;i++) e[i-1]=e[i]; // Convenient to renumber the elements of e.
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		if (verbose)
			if (fmod(l,tfs) < 1.0)
				mprintf(WHITE,"#");
		iter=0;
		do {
			for (m=l;m<=n-1;m++) { // Look for a single small subdiagonal element to split the matrix.
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) eprintf("\nError: Too many iterations in tqli.\n");
				g=(d[l+1]-d[l])/(2.0*e[l]); // Form shift.
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g)); // This is dm - ks.
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) { // A plane rotation as in the original QL, followed by Givens rotations to restore tridiagonal form.
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) { // Recover from underflow.
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					/* Next loop can be omitted if eigenvectors not wanted */
					for (k=1;k<=n;k++) { // Form eigenvectors.
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
	if (verbose)
		mprintf(WHITE,"]\n");
}


void TestDiagonalization() {

	CxDMatrixMN inmat, evec;
	CxDoubleArray eval;
	int n = 6, z;

	inmat = CxDMatrixMN(n,n);

	/* Mathematica
	{{302.700231596574, 170.443219393374, -139.038947072181, -46.245819859365, -103.181902861864, 29.010383308527},
	{170.443219393374, 260.815132417353, -12.935985300894, 44.502937022346, -94.698662374553, -43.419516133519},
	{-139.038947072181, -12.935985300894, 217.193731596619, 128.387185615870, 43.079225195861, -69.443726271268},
	{-46.245819859365, 44.502937022346, 128.387185615870, 209.996615121875, 74.217503978145, -96.339813050717},
	{-103.181902861864, -94.698662374553, 43.079225195861, 74.217503978145, 208.245347398560, -102.388238316016},
	{29.010383308527, -43.419516133519, -69.443726271268, -96.339813050717, -102.388238316016, 138.236759834129}} */

	inmat(0,0) = 302.700231596574;
	inmat(0,1) = 170.443219393374;
	inmat(0,2) = -139.038947072181;
	inmat(0,3) = -46.245819859365;
	inmat(0,4) = -103.181902861864;
	inmat(0,5) = 29.010383308527;
	inmat(1,0) = 170.443219393374;
	inmat(1,1) = 260.815132417353;
	inmat(1,2) = -12.935985300894;
	inmat(1,3) = 44.502937022346;
	inmat(1,4) = -94.698662374553;
	inmat(1,5) = -43.419516133519;
	inmat(2,0) = -139.038947072181;
	inmat(2,1) = -12.935985300894;
	inmat(2,2) = 217.193731596619;
	inmat(2,3) = 128.387185615870;
	inmat(2,4) = 43.079225195861;
	inmat(2,5) = -69.443726271268;
	inmat(3,0) = -46.245819859365;
	inmat(3,1) = 44.502937022346;
	inmat(3,2) = 128.387185615870;
	inmat(3,3) = 209.996615121875;
	inmat(3,4) = 74.217503978145;
	inmat(3,5) = -96.339813050717;
	inmat(4,0) = -103.181902861864;
	inmat(4,1) = -94.698662374553;
	inmat(4,2) = 43.079225195861;
	inmat(4,3) = 74.217503978145;
	inmat(4,4) = 208.245347398560;
	inmat(4,5) = -102.388238316016;
	inmat(5,0) = 29.010383308527;
	inmat(5,1) = -43.419516133519;
	inmat(5,2) = -69.443726271268;
	inmat(5,3) = -96.339813050717;
	inmat(5,4) = -102.388238316016;
	inmat(5,5) = 138.236759834129;

	mprintf("Input Matrix:\n");

	inmat.Dump();

	mprintf("\nDiagonalizing...\n");

	Diagonalize(inmat,&evec,&eval,true);

	mprintf("Eigenvalues:\n");

	for (z=0;z<n;z++)
		printf("%f\n",eval[z]);

	mprintf("\nEigenvectors:\n");
	evec.Dump();

	/* Mathematica Eigensystem[] result
	{{601.9106394832, 387.2145035075, 197.14653642685, 79.13859090556, 49.15230921061, 22.625238431395},
	{{0.6156803019077, 0.3703696018678, -0.4207499291255, -0.30015545190547, -0.4102706883168, 0.21981421245112},
	{0.26646953855962, 0.6223460655665, 0.31137391907636, 0.5401190235001, 0.03868016185508, -0.3892312008883},
	{0.31753278850618, -0.12959558740343, -0.5100908695359, -0.015818987016767, 0.7069115299446, -0.3495869923871},
	{-0.4125610205684, 0.4687247567437, 0.003398913549180, -0.6510228490342, 0.04697698892909, -0.4290004273383},
	{-0.5276787728848, 0.28596095292004, -0.6514921326311, 0.4264753464293, -0.09330630273689, 0.15732819178134},
	{-0.021288401548652, 0.3965351841040, 0.20339501710323, -0.11052571518852, 0.5652831083617, 0.6849643731083}}} */

	return;
}

