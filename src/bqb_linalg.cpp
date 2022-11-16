/***********************************************************************************

    LibBQB - File Format and Compression Algorithms for Trajectories of
             Volumetric Data and Atom Positions

    https://brehm-research.de/bqb

    Free software, licensed under GNU LGPL v3

    Copyright (c) Martin Brehm and Martin Thomas,
                  Martin Luther University Halle-Wittenberg, Germany,
                  2016 - 2022.

    Please cite:  M. Brehm, M. Thomas: "An Efficient Lossless Compression Algorithm
                  for Trajectories of Atom Positions and Volumetric Data",
                  J. Chem. Inf. Model. 2018, 58 (10), pp 2092-2107.

    --------------------------------------------------------------------------------

    LibBQB is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

***********************************************************************************/



// This must always be the first include directive
#include "bqb_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bqb_interface.h"
#include "bqb_linalg.h"
#include "bqb_tools.h"


const char *GetRevisionInfo_bqb_linalg(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_linalg() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


//#define MIN(x,y) ( (x) < (y) ? (x) : (y) )
//#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SIGNF(a, b) ((b) >= 0.0 ? fabsf(a) : -fabsf(a))
//#define SQR(x) ((x)*(x))


/* 
 * svdcomp - SVD decomposition routine. 
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is 
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/

 
static double PYTHAG(double a, double b)
{
    double at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}


/* m - rows, n - cols
   a[zm][zn] = a[zm*n + zn]; */
int CBQBLinAlg::BQB_ComputeSVD_Flat(double *a, int m, int n, double *w, double *v)
{
    int flag, i, its, j, jj, k, l = 0, nm = 0;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        m_IF.eprintf("ComputeSVD(): #rows must be > #cols \n");
        return(0);
    }
  
    rv1 = (double *)malloc((unsigned int) n*sizeof(double));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++) 
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m) 
        {
            for (k = i; k < m; k++) 
                scale += fabs(a[k*n+i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k*n+i] = a[k*n+i]/scale;
                    s += a[k*n+i] * a[k*n+i];
                }
                f = a[i*n+i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i*n+i] = f - g;
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += a[k*n+i] * a[k*n+j];
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k*n+j] += f * a[k*n+i];
                    }
                }
                for (k = i; k < m; k++) 
                    a[k*n+i] = a[k*n+i]*scale;
            }
        }
        w[i] = scale * g;
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs(a[i*n+k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i*n+k] = a[i*n+k]/scale;
                    s += a[i*n+k] * a[i*n+k];
                }
                f = a[i*n+l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i*n+l] = f - g;
                for (k = l; k < n; k++) 
                    rv1[k] = a[i*n+k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += a[j*n+k] * a[i*n+k];
                        for (k = l; k < n; k++) 
                            a[j*n+k] += s * rv1[k];
                    }
                }
                for (k = l; k < n; k++) 
                    a[i*n+k] = a[i*n+k]*scale;
            }
        }
        anorm = MAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }
  
    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        if (i < n - 1) 
        {
            if (g) 
            {
                for (j = l; j < n; j++)
                    v[j*n+i] = (a[i*n+j] / a[i*n+l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += a[i*n+k] * v[k*n+j];
                    for (k = l; k < n; k++) 
                        v[k*n+j] += s * v[k*n+i];
                }
            }
            for (j = l; j < n; j++) 
                v[i*n+j] = v[j*n+i] = 0.0;
        }
        v[i*n+i] = 1.0;
        g = rv1[i];
        l = i;
    }
  
    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--) 
    {
        l = i + 1;
        g = w[i];
        if (i < n - 1) 
            for (j = l; j < n; j++) 
                a[i*n+j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += a[k*n+i] * a[k*n+j];
                    f = (s / a[i*n+i]) * g;
                    for (k = i; k < m; k++) 
                        a[k*n+j] += f * a[k*n+i];
                }
            }
            for (j = i; j < m; j++) 
                a[j*n+i] = a[j*n+i]*g;
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j*n+i] = 0.0;
        }
        ++a[i*n+i];
    }

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--) 
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++) 
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--) 
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm) 
                {
                    flag = 0;
                    break;
                }
                if (fabs(w[nm]) + anorm == anorm) 
                    break;
            }
            if (flag) 
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++) 
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm) 
                    {
                        g = w[i];
                        h = PYTHAG(f, g);
                        w[i] = h; 
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++) 
                        {
                            y = a[j*n+nm];
                            z = a[j*n+i];
                            a[j*n+nm] = y * c + z * s;
                            a[j*n+i] = z * c - y * s;
                        }
                    }
                }
            }
            z = w[k];
            if (l == k) 
            {                  /* convergence */
                if (z < 0.0) 
                {              /* make singular value nonnegative */
                    w[k] = -z;
                    for (j = 0; j < n; j++) 
                        v[j*n+k] = (-v[j*n+k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                m_IF.eprintf("ComputeSVD(): No convergence after 30,000! iterations \n");
                return(0);
            }
    
            /* shift from bottom 2 x 2 minor */
            x = w[l];
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
          
            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++) 
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++) 
                {
                    x = v[jj*n+j];
                    z = v[jj*n+i];
                    v[jj*n+j] = x * c + z * s;
                    v[jj*n+i] = z * c - x * s;
                }
                z = PYTHAG(f, h);
                w[j] = z;
                if (z) 
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++) 
                {
                    y = a[jj*n+j];
                    z = a[jj*n+i];
                    a[jj*n+j] = y * c + z * s;
                    a[jj*n+i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    free((void*) rv1);
    return(1);
}


void CBQBLinAlg::BQB_ComputePseudoInverse(int m, double *m_in, double *m_out)
{
	double *mat_a, *mat_b, *mat_c, tf;
	int z1, z2, z3;
//	int z;

	mat_a = new double[m*m];
	mat_b = new double[m];
	mat_c = new double[m*m];

	memcpy(mat_a,m_in,m*m*sizeof(double));

	BQB_ComputeSVD_Flat(mat_a,m,m,mat_b,mat_c);

/*	m_IF.printf("*** U Matrix:\n");
	m_IF.printf("{ ");
	for (z=0;z<m;z++)
	{
		m_IF.printf("{ ");
		for (z2=0;z2<m;z2++)
		{
			m_IF.printf("%10G",mat_a[z*m+z2]);
			if (z2 < m-1)
				m_IF.printf(", ");
		}
		m_IF.printf(" }");
		if (z < m-1)
			m_IF.printf(",\n");
	}
	m_IF.printf(" }\n\n");

	m_IF.printf("*** D Matrix:\n");
	m_IF.printf("{ ");
	for (z=0;z<m;z++)
	{
		m_IF.printf("{ ");
		for (z2=0;z2<m;z2++)
		{
			if (z == z2)
				m_IF.printf("%10G",mat_b[z]);
					else m_IF.printf("%10G",0.0);
			if (z2 < m-1)
				m_IF.printf(", ");
		}
		m_IF.printf(" }");
		if (z < m-1)
			m_IF.printf(",\n");
	}
	m_IF.printf(" }\n\n");

	m_IF.printf("*** V Matrix:\n");
	m_IF.printf("{ ");
	for (z=0;z<m;z++)
	{
		m_IF.printf("{ ");
		for (z2=0;z2<m;z2++)
		{
			m_IF.printf("%10G",mat_c[z*m+z2]);
			if (z2 < m-1)
				m_IF.printf(", ");
		}
		m_IF.printf(" }");
		if (z < m-1)
			m_IF.printf(",\n");
	}
	m_IF.printf(" }\n\n");*/

	for (z1=0;z1<m;z1++)
		if (mat_b[z1] != 0)
			mat_b[z1] = 1.0 / mat_b[z1];

	for (z1=0;z1<m;z1++)
	{
		for (z2=0;z2<m;z2++)
		{
			tf = 0;

			for (z3=0;z3<m;z3++)
				tf += mat_c[z1*m+z3] * mat_b[z3] * mat_a[z2*m+z3];

			m_out[z1*m+z2] = tf;
		}
	}

/*	m_IF.printf("*** Pseudoinverse:\n");
	m_IF.printf("{ ");
	for (z=0;z<m;z++)
	{
		m_IF.printf("{ ");
		for (z2=0;z2<m;z2++)
		{
			m_IF.printf("%10G",m_out[z*m+z2]);
			if (z2 < m-1)
				m_IF.printf(", ");
		}
		m_IF.printf(" }");
		if (z < m-1)
			m_IF.printf(",\n");
	}
	m_IF.printf(" }\n\n");*/

	delete[] mat_a;
	delete[] mat_b;
	delete[] mat_c;
}


CBQBDMatrixMN CBQBLinAlg::BQB_ComputePseudoInverse(CBQBDMatrixMN input) {

	double *mat_a, *mat_b, *mat_c, tf;
	int z1, z2, z3;
	int m, n;
	CBQBDMatrixMN mout;


	m = input.GetRows();
	n = input.GetCols();

	mat_a = new double[m*n];
	mat_b = new double[n];
	mat_c = new double[n*n];

	for (z1=0;z1<m;z1++)
		for (z2=0;z2<n;z2++)
			mat_a[z1*n+z2] = input(z1,z2);

	BQB_ComputeSVD_Flat(mat_a,m,n,mat_b,mat_c);

//	m_IF.printf("@@@");
	for (z1=0;z1<n;z1++) {
//		m_IF.printf(" %.3E",mat_b[z1]);
		if (mat_b[z1] != 0) {
			if (fabs(mat_b[z1]) > 1.0E-13)
				mat_b[z1] = 1.0 / mat_b[z1];
			else
				mat_b[z1] = 0;
		}
	}
//	m_IF.printf(" @@@\n");

	mout = CBQBDMatrixMN(n,m);

	for (z1=0;z1<n;z1++) {

		for (z2=0;z2<m;z2++) {

			tf = 0;

			for (z3=0;z3<n;z3++)
				tf += mat_c[z1*n+z3] * mat_b[z3] * mat_a[z2*n+z3];

			mout(z1,z2) = tf;
		}
	}

	delete[] mat_a;
	delete[] mat_b;
	delete[] mat_c;

	return mout;
}


void CBQBLinAlg::BQB_TestPseudoInverse()
{
	double *mat_in, *mat_out;
	int mat_M, z, z2;

	mat_M = 5;

	mat_in = new double[mat_M*mat_M];
	mat_out = new double[mat_M*mat_M];

	for (z=0;z<mat_M;z++)
	{
//		mat_a[z] = new double[mat_M];
//		mat_c[z] = new double[mat_M];
		for (z2=0;z2<mat_M;z2++)
			mat_in[z*mat_M+z2] = ((rand()%20001)-10000)/1000.0;
	}

	m_IF.printf("*** Input Matrix:\n");
	m_IF.printf("{ ");
	for (z=0;z<mat_M;z++)
	{
		m_IF.printf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			m_IF.printf("%10G",mat_in[z*mat_M+z2]);
			if (z2 < mat_M-1)
				m_IF.printf(", ");
		}
		m_IF.printf(" }");
		if (z < mat_M-1)
			m_IF.printf(",\n");
	}
	m_IF.printf(" }\n\n");

	BQB_ComputePseudoInverse(mat_M,mat_in,mat_out);

	m_IF.printf("*** Moore-Penrose Pseudoinverse:\n");
	m_IF.printf("{ ");
	for (z=0;z<mat_M;z++)
	{
		m_IF.printf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			m_IF.printf("%10G",mat_out[z*mat_M+z2]);
			if (z2 < mat_M-1)
				m_IF.printf(", ");
		}
		m_IF.printf(" }");
		if (z < mat_M-1)
			m_IF.printf(",\n");
	}
	m_IF.printf(" }\n\n");

	delete[] mat_in;
	delete[] mat_out;
}


void CBQBLinAlg::BQB_TestSVD()
{
	// SVD-Test

	double *mat_a, *mat_b, *mat_c;
	int mat_M, z, z2;

	mat_M = 5;

	mat_a = new double[mat_M*mat_M];
	mat_b = new double[mat_M];
	mat_c = new double[mat_M*mat_M];

	for (z=0;z<mat_M;z++)
	{
//		mat_a[z] = new double[mat_M];
//		mat_c[z] = new double[mat_M];
		for (z2=0;z2<mat_M;z2++)
			mat_a[z*mat_M+z2] = ((rand()%20001)-10000)/1000.0;
	}

	m_IF.printf("*** Input Matrix:\n");
	m_IF.printf("{ ");
	for (z=0;z<mat_M;z++)
	{
		m_IF.printf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			m_IF.printf("%10G",mat_a[z*mat_M+z2]);
			if (z2 < mat_M-1)
				m_IF.printf(", ");
		}
		m_IF.printf(" }");
		if (z < mat_M-1)
			m_IF.printf(",\n");
	}
	m_IF.printf(" }\n\n");

	BQB_ComputeSVD_Flat(mat_a,mat_M,mat_M,mat_b,mat_c);

	m_IF.printf("*** U Matrix:\n");
	m_IF.printf("{ ");
	for (z=0;z<mat_M;z++)
	{
		m_IF.printf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			m_IF.printf("%10G",mat_a[z*mat_M+z2]);
			if (z2 < mat_M-1)
				m_IF.printf(", ");
		}
		m_IF.printf(" }");
		if (z < mat_M-1)
			m_IF.printf(",\n");
	}
	m_IF.printf(" }\n\n");

	m_IF.printf("*** D Matrix:\n");
	m_IF.printf("{ ");
	for (z=0;z<mat_M;z++)
	{
		m_IF.printf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			if (z == z2)
				m_IF.printf("%10G",mat_b[z]);
					else m_IF.printf("%10G",0.0);
			if (z2 < mat_M-1)
				m_IF.printf(", ");
		}
		m_IF.printf(" }");
		if (z < mat_M-1)
			m_IF.printf(",\n");
	}
	m_IF.printf(" }\n\n");

	m_IF.printf("*** V Matrix:\n");
	m_IF.printf("{ ");
	for (z=0;z<mat_M;z++)
	{
		m_IF.printf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			m_IF.printf("%10G",mat_c[z*mat_M+z2]);
			if (z2 < mat_M-1)
				m_IF.printf(", ");
		}
		m_IF.printf(" }");
		if (z < mat_M-1)
			m_IF.printf(",\n");
	}
	m_IF.printf(" }\n\n");
}

