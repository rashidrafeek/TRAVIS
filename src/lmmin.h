/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

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



/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lmmin.h
 *
 * Contents: Public interface to the Levenberg-Marquardt core implementation.
 *
 * Author:   Joachim Wuttke 2004-2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 */
 
#ifndef LMMIN_H
#define LMMIN_H

//#ifdef __cplusplus
//extern "C" {
//#endif


// This must always be the first include directive
#include "config.h"

#include <float.h>


/** Compact high-level interface. **/

/* Collection of control (input) parameters. */
typedef struct {
    double ftol;      /* relative error desired in the sum of squares. */
    double xtol;      /* relative error between last two approximations. */
    double gtol;      /* orthogonality desired between fvec and its derivs. */
    double epsilon;   /* step used to calculate the jacobian. */
    double stepbound; /* initial bound to steps in the outer loop. */
    int maxcall;      /* maximum number of iterations. */
    int scale_diag;   /* UNDOCUMENTED, TESTWISE automatical diag rescaling? */
    int printflags;   /* OR'ed to produce more noise */
} lm_control_struct;

/* Collection of status (output) parameters. */
typedef struct {
    double fnorm;     /* norm of the residue vector fvec. */
    int nfev;	      /* actual number of iterations. */
    int info;	      /* status of minimization. */
} lm_status_struct;

/* Recommended control parameter settings. */
extern const lm_control_struct lm_control_double;
extern const lm_control_struct lm_control_float;

/* Standard monitoring routine. */
void lm_printout_std( int n_par, const double *par, int m_dat,
                      const void *data, const double *fvec,
                      int printflags, int iflag, int iter, int nfev);

/* Refined calculation of Eucledian norm, typically used in printout routine. */
double lm_enorm( int, const double * );

/* The actual minimization. */
void lmmin( int n_par, double *par, int m_dat, const void *data, 
            void (*evaluate) (const double *par, int m_dat, const void *data,
                              double *fvec, int *info),
            const lm_control_struct *control, lm_status_struct *status,
            void (*printout) (int n_par, const double *par, int m_dat,
                              const void *data, const double *fvec,
                              int printflags, int iflag, int iter, int nfev) );


/** Legacy low-level interface. **/

/* Alternative to lm_minimize, allowing full control, and read-out
   of auxiliary arrays. For usage, see implementation of lm_minimize. */
void lm_lmdif( int m, int n, double *x, double *fvec, double ftol,
	       double xtol, double gtol, int maxfev, double epsfcn,
	       double *diag, int mode, double factor, int *info, int *nfev,
	       double *fjac, int *ipvt, double *qtf, double *wa1,
	       double *wa2, double *wa3, double *wa4,
               void (*evaluate) (const double *par, int m_dat, const void *data,
                                 double *fvec, int *info),
               void (*printout) (int n_par, const double *par, int m_dat,
                                 const void *data, const double *fvec,
                                 int printflags, int iflag, int iter, int nfev),
               int printflags, const void *data );

extern const char *lm_infmsg[];
extern const char *lm_shortmsg[];


/*****************************************************************************/
/*  set numeric constants                                                    */
/*****************************************************************************/

/* machine-dependent constants from float.h */
#define LM_MACHEP     DBL_EPSILON   /* resolution of arithmetic */
#define LM_DWARF      DBL_MIN       /* smallest nonzero number */
#define LM_SQRT_DWARF sqrt(DBL_MIN) /* square should not underflow */
#define LM_SQRT_GIANT sqrt(DBL_MAX) /* square should not overflow */
#define LM_USERTOL    30*LM_MACHEP  /* users are recommended to require this */

/* If the above values do not work, the following seem good for an x86:
 LM_MACHEP     .555e-16
 LM_DWARF      9.9e-324	
 LM_SQRT_DWARF 1.e-160   
 LM_SQRT_GIANT 1.e150 
 LM_USERTOL   1.e-14*/

 /*  The following values should work on any machine:
#define LM_MACHEP     1.2e-16
#define LM_DWARF      1.0e-38
#define LM_SQRT_DWARF 3.834e-20
#define LM_SQRT_GIANT 1.304e19
#define LM_USERTOL    1.e-14*/


/*const lm_control_struct lm_control_double = {
    LM_USERTOL, LM_USERTOL, LM_USERTOL, LM_USERTOL, 1.0, 2000, 1, 0 };

const lm_control_struct lm_control_float = {
    1.e-7, 1.e-7, 1.e-7, 1.e-7, 1.0, 100, 0, 0 };*/



//#ifdef __cplusplus
//}
//#endif

#endif /* LMMIN_H */
