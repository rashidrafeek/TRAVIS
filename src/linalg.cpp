/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

    This file was written by Martin Brehm.

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

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

#include "linalg.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tools.h"
#include "qr_fact.h"
#include "lu_decomp.h"
//#include "defs_and_types.h"


const char *GetRevisionInfo_linalg(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_linalg() {
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


int ComputeSVD(double **a, int m, int n, double *w, double **v)
{
    int flag, i, its, j, jj, k, l = 0, nm = 0;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        eprintf("ComputeSVD(): #rows must be > #cols \n");
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
                scale += fabs(a[k][i]);
            if (scale) 
            {
                for (k = i; k < m; k++) 
                {
                    a[k][i] = a[k][i]/scale;
                    s += a[k][i] * a[k][i];
                }
                f = a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = f - g;
                if (i != n - 1) 
                {
                    for (j = l; j < n; j++) 
                    {
                        for (s = 0.0, k = i; k < m; k++) 
                            s += a[k][i] * a[k][j];
                        f = s / h;
                        for (k = i; k < m; k++) 
                            a[k][j] += f * a[k][i];
                    }
                }
                for (k = i; k < m; k++) 
                    a[k][i] = a[k][i]*scale;
            }
        }
        w[i] = scale * g;
    
        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1) 
        {
            for (k = l; k < n; k++) 
                scale += fabs(a[i][k]);
            if (scale) 
            {
                for (k = l; k < n; k++) 
                {
                    a[i][k] = a[i][k]/scale;
                    s += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = f - g;
                for (k = l; k < n; k++) 
                    rv1[k] = a[i][k] / h;
                if (i != m - 1) 
                {
                    for (j = l; j < m; j++) 
                    {
                        for (s = 0.0, k = l; k < n; k++) 
                            s += a[j][k] * a[i][k];
                        for (k = l; k < n; k++) 
                            a[j][k] += s * rv1[k];
                    }
                }
                for (k = l; k < n; k++) 
                    a[i][k] = a[i][k]*scale;
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
                    v[j][i] = (a[i][j] / a[i][l]) / g;
                    /* double division to avoid underflow */
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < n; k++) 
                        s += a[i][k] * v[k][j];
                    for (k = l; k < n; k++) 
                        v[k][j] += s * v[k][i];
                }
            }
            for (j = l; j < n; j++) 
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
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
                a[i][j] = 0.0;
        if (g) 
        {
            g = 1.0 / g;
            if (i != n - 1) 
            {
                for (j = l; j < n; j++) 
                {
                    for (s = 0.0, k = l; k < m; k++) 
                        s += a[k][i] * a[k][j];
                    f = (s / a[i][i]) * g;
                    for (k = i; k < m; k++) 
                        a[k][j] += f * a[k][i];
                }
            }
            for (j = i; j < m; j++) 
                a[j][i] = a[j][i]*g;
        }
        else 
        {
            for (j = i; j < m; j++) 
                a[j][i] = 0.0;
        }
        ++a[i][i];
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
                            y = a[j][nm];
                            z = a[j][i];
                            a[j][nm] = y * c + z * s;
                            a[j][i] = z * c - y * s;
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
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                eprintf("ComputeSVD(): No convergence after 30,000! iterations \n");
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
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
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
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y * c + z * s;
                    a[jj][i] = z * c - y * s;
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


/* m - rows, n - cols
   a[zm][zn] = a[zm*n + zn]; */
int ComputeSVD_Flat(double *a, int m, int n, double *w, double *v)
{
    int flag, i, its, j, jj, k, l = 0, nm = 0;
    double c, f, h, s, x, y, z;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < n) 
    {
        eprintf("ComputeSVD(): #rows must be > #cols \n");
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
                eprintf("ComputeSVD(): No convergence after 30,000! iterations \n");
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


void ComputePseudoInverse(int m, double *m_in, double *m_out)
{
	double *mat_a, *mat_b, *mat_c, tf;
	int z1, z2, z3;
//	int z;

	mat_a = new double[m*m];
	mat_b = new double[m];
	mat_c = new double[m*m];

	memcpy(mat_a,m_in,m*m*sizeof(double));

	ComputeSVD_Flat(mat_a,m,m,mat_b,mat_c);

/*	mprintf("*** U Matrix:\n");
	mprintf("{ ");
	for (z=0;z<m;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<m;z2++)
		{
			mprintf("%10G",mat_a[z*m+z2]);
			if (z2 < m-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < m-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	mprintf("*** D Matrix:\n");
	mprintf("{ ");
	for (z=0;z<m;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<m;z2++)
		{
			if (z == z2)
				mprintf("%10G",mat_b[z]);
					else mprintf("%10G",0.0);
			if (z2 < m-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < m-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	mprintf("*** V Matrix:\n");
	mprintf("{ ");
	for (z=0;z<m;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<m;z2++)
		{
			mprintf("%10G",mat_c[z*m+z2]);
			if (z2 < m-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < m-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");*/

	for (z1=0;z1<m;z1++)
		if (fabs(mat_b[z1]) > 1.0E-13)
			mat_b[z1] = 1.0 / mat_b[z1];
		else
			mat_b[z1] = 0;

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

/*	mprintf("*** Pseudoinverse:\n");
	mprintf("{ ");
	for (z=0;z<m;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<m;z2++)
		{
			mprintf("%10G",m_out[z*m+z2]);
			if (z2 < m-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < m-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");*/

	delete[] mat_a;
	delete[] mat_b;
	delete[] mat_c;
}


CxDMatrixMN ComputePseudoInverse(CxDMatrixMN input) {

	double *mat_a, *mat_b, *mat_c, tf;
	int z1, z2, z3;
	int m, n;
	CxDMatrixMN mout;


	m = input.GetRows();
	n = input.GetCols();

	mat_a = new double[m*n];
	mat_b = new double[n];
	mat_c = new double[n*n];

	for (z1=0;z1<m;z1++)
		for (z2=0;z2<n;z2++)
			mat_a[z1*n+z2] = input(z1,z2);

	ComputeSVD_Flat(mat_a,m,n,mat_b,mat_c);

//	mprintf("@@@");
	for (z1=0;z1<n;z1++) {
//		mprintf(" %.3E",mat_b[z1]);
//		if (mat_b[z1] != 0) {
			if (fabs(mat_b[z1]) > 1.0E-13)
				mat_b[z1] = 1.0 / mat_b[z1];
			else
				mat_b[z1] = 0;
//		}
	}
//	mprintf(" @@@\n");

	mout = CxDMatrixMN(n,m);

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


void TestPseudoInverse()
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

	mprintf("*** Input Matrix:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			mprintf("%10G",mat_in[z*mat_M+z2]);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	ComputePseudoInverse(mat_M,mat_in,mat_out);

	mprintf("*** Moore-Penrose Pseudoinverse:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			mprintf("%10G",mat_out[z*mat_M+z2]);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	delete[] mat_in;
	delete[] mat_out;
}


void TestSVD()
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

	mprintf("*** Input Matrix:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			mprintf("%10G",mat_a[z*mat_M+z2]);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	ComputeSVD_Flat(mat_a,mat_M,mat_M,mat_b,mat_c);

	mprintf("*** U Matrix:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			mprintf("%10G",mat_a[z*mat_M+z2]);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	mprintf("*** D Matrix:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			if (z == z2)
				mprintf("%10G",mat_b[z]);
					else mprintf("%10G",0.0);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");

	mprintf("*** V Matrix:\n");
	mprintf("{ ");
	for (z=0;z<mat_M;z++)
	{
		mprintf("{ ");
		for (z2=0;z2<mat_M;z2++)
		{
			mprintf("%10G",mat_c[z*mat_M+z2]);
			if (z2 < mat_M-1)
				mprintf(", ");
		}
		mprintf(" }");
		if (z < mat_M-1)
			mprintf(",\n");
	}
	mprintf(" }\n\n");
}


void Solve_LeastSquares_QR(double *a, double *b, double *x, int m, int n)
{
	double **ta;
	int z;

	ta = new double*[m];
	for (z=0;z<m;z++)
	{
		ta[z] = new double[n];
		memcpy(ta[z],&a[z*n],n*sizeof(double));
	}

	QR_least_squares(ta, b, x, m, n);

	for (z=0;z<m;z++)
		delete[] ta[z];
	delete[] ta;
}


void TestLeastSquaresQR()
{
#define M 6
#define N 3

	double A[M*N], b[M], x[N];
	int z, z2;


	A[0*N+0] = 1;   A[0*N+1] = 2;   A[0*N+2] = 3;    b[0] = 3; 
	A[1*N+0] = 4;   A[1*N+1] = 5;   A[1*N+2] = 6;    b[1] = 9; 
	A[2*N+0] = 7;   A[2*N+1] = 8;   A[2*N+2] = 9;    b[2] = 15; 
	A[3*N+0] = 10;  A[3*N+1] = 11;  A[3*N+2] = 12;   b[3] = 22; 
	A[4*N+0] = 13;  A[4*N+1] = 14;  A[4*N+2] = 15;   b[4] = 27; 
	A[5*N+0] = 16;  A[5*N+1] = 17;  A[5*N+2] = -5;   b[5] = 33;

	mprintf("Testing Least-Squares Solver for system:\n\n");

	for (z=0;z<M;z++)
	{
		mprintf(" |");
		for (z2=0;z2<N;z2++)
		{
			mprintf(" %7.3f ",A[z2+z*N]);
		}
		mprintf("|");
		if (z == M/2)
			mprintf(" * ");
				else mprintf("   ");
		if ((z > (M-N)/2) && (z-(M-N)/2 <= N))
			mprintf("| x%d |",(z-(M-N)/2));
				else mprintf("      ");
		if (z == M/2)
			mprintf(" = ");
				else mprintf("   ");
		mprintf("| %7.3f |\n",b[z]);
	}
	mprintf("\n");

	mprintf("Solving...\n");
	Solve_LeastSquares_QR(A,b,x,M,N);

	mprintf("\nResult:\n\n");
	for (z=0;z<N;z++)
		mprintf("  x%d = %G\n",z+1,x[z]);
	mprintf("\n");
}


void SolveLinearLU(CxDMatrixMN *a, CxDVectorN *b, CxDVectorN *x) {

	double **ta;
	double *tb, /*  *tx,  */ d;
	int *indx;
	int z, z2, n;

	if ((a->GetCols() != b->GetDim()) || (a->GetRows() != b->GetDim())) {
		eprintf("SolveLinearLU(): Error: Dimension mismatch (matrix %d x %d, vector %d).\n",a->GetCols(),a->GetRows(),b->GetDim());
		abort();
	}

	n = a->GetCols();

	ta = new double*[n+1];
	tb = new double[n+1];
	//tx = new double[n+1];
	indx = new int[n+1];
	for (z=0;z<n;z++) {
		ta[z+1] = new double[n+1];
		for (z2=0;z2<n;z2++)
			ta[z+1][z2+1] = (*a)(z,z2);
		tb[z+1] = (*b)[z];
	}

	ludcmp(ta,n,indx,&d);

/*	mprintf("Matrix:\n");
	for (z=0;z<n;z++) {
		for (z2=0;z2<n;z2++)
			mprintf("  %12.8f",ta[z+1][z2+1]);
		mprintf("\n");
	}
	mprintf("\nIndex:\n");
	for (z=0;z<n;z++)
		mprintf("  %d",indx[z+1]);
	mprintf("\n");*/

	lubksb(ta,n,indx,tb);

	if (x->GetDim() != n)
		*x = CxDVectorN(n);

	for (z=0;z<n;z++)
		(*x)[z] = tb[z+1];

	delete[] indx;
	delete[] tb;
	//delete[] tx;
	for (z=0;z<n;z++)
		delete[] ta[z+1];
	delete[] ta;
}


void TestLinearLU() {

/*	Mathematica Input:

	mat = {{2.23987, 8.40283, -6.59509, -3.27052, 0.367372},
	{-7.68671, 7.21389,  2.78113, -3.51783, 3.63235},
	{-9.17321, 6.71374, -2.48653, 6.97584,  4.37098},
	{-5.90448, 8.52263, -3.77812, -8.7871,  1.80073},
	{7.15649, 0.948323, -6.86932, 9.37041, -5.43305}};

	vec = {7.64096, 0.579755, 1.44277, 7.44359, 7.72779};
	
	LinearSolve[mat, vec] gives:
	
	{-0.289021, 0.540435, -0.569921, -0.133297, -1.21805}  */


	CxDMatrixMN mat;
	CxDVectorN b, x;


	mat = CxDMatrixMN(5,5);

	mat(0,0) = 2.23987;
	mat(0,1) = 8.40283;
	mat(0,2) = -6.59509;
	mat(0,3) = -3.27052;
	mat(0,4) = 0.367372;
	mat(1,0) = -7.68671;
	mat(1,1) = 7.21389;
	mat(1,2) = 2.78113;
	mat(1,3) = -3.51783;
	mat(1,4) = 3.63235;
	mat(2,0) = -9.17321;
	mat(2,1) = 6.71374;
	mat(2,2) = -2.48653;
	mat(2,3) = 6.97584;
	mat(2,4) = 4.37098;
	mat(3,0) = -5.90448;
	mat(3,1) = 8.52263;
	mat(3,2) = -3.77812;
	mat(3,3) = -8.7871;
	mat(3,4) = 1.80073;
	mat(4,0) = 7.15649;
	mat(4,1) = 0.948323;
	mat(4,2) = -6.86932;
	mat(4,3) = 9.37041;
	mat(4,4) = -5.43305;

	b = CxDVectorN(5);

	b[0] = 7.64096;
	b[1] = 0.579755;
	b[2] = 1.44277;
	b[3] = 7.44359;
	b[4] = 7.72779;

	SolveLinearLU(&mat,&b,&x);

	x.Dump();
}


