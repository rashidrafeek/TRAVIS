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



/******************************************************************************/
/*

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 December 20080

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.
*/	 
		
#ifndef ZIGGURAT_H
#define ZIGGURAT_H


// This must always be the first include directive
#include "config.h"



double r4_exp ( unsigned long int *jsr, int ke[256], double fe[256], 
  double we[256] );
void r4_exp_setup ( int ke[256], double fe[256], double we[256] );
double r4_nor ( unsigned long int *jsr, int kn[128], double fn[128], 
  double wn[128] );
void r4_nor_setup ( int kn[128], double fn[128], double wn[128] );
double r4_uni ( unsigned long int *jsr );
unsigned long int shr3 ( unsigned long int *jsr );
void timestamp ( void );


#endif
