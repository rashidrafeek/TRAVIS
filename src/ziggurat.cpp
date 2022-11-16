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


// This must always be the first include directive
#include "config.h"

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

#include "ziggurat.h"
#include "tools.h"


const char *GetRevisionInfo_ziggurat(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_ziggurat() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



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



/******************************************************************************/

double r4_exp ( unsigned long int *jsr, int ke[256], double fe[256], 
  double we[256] )

/******************************************************************************/
/*
  Purpose:

    R4_EXP returns an exponentially distributed single precision real value.

  Discussion:

    The underlying algorithm is the ziggurat method.

    Before the first call to this function, the user must call R4_EXP_SETUP
    to determine the values of KE, FE and WE.

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

  Parameters:

    Input/output, unsigned long int *JSR, the seed.

    Input, int KE[256], data computed by R4_EXP_SETUP.

    Input, double FE[256], WE[256], data computed by R4_EXP_SETUP.

    Output, double R4_EXP, an exponentially distributed random value.
*/
{
  int iz;
  int jz;
  double value;
  double x;

  jz = shr3 ( jsr );
  iz = ( jz & 255 );

  if ( abs ( jz  ) < ke[iz] )
  {
    value = ( double ) ( abs ( jz ) ) * we[iz];
  }
  else
  {
    for ( ; ; )
    {
      if ( iz == 0 )
      {
        value = 7.69711 - log ( r4_uni ( jsr ) );
        break;
      }

      x = ( double ) ( abs ( jz ) ) * we[iz];

      if ( fe[iz] + r4_uni ( jsr ) * ( fe[iz-1] - fe[iz] ) < exp ( - x ) )
      {
        value = x;
        break;
      }

      jz = shr3 ( jsr );
      iz = ( jz & 255 );

      if ( abs ( jz ) < ke[iz] )
      {
        value = ( double ) ( abs ( jz ) ) * we[iz];
        break;
      }
    }
  }
  return value;
}
/******************************************************************************/

void r4_exp_setup ( int ke[256], double fe[256], double we[256] )

/******************************************************************************/
/*
  Purpose:

    R4_EXP_SETUP sets data needed by R4_EXP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 December 2008

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

  Parameters:

    Output, int KE[256], data needed by R4_EXP.

    Output, double FE[256], WE[256], data needed by R4_EXP.
*/
{
  double de = 7.697117470131487;
  int i;
  const double m2 = 2147483648.0;
  double q;
  double te = 7.697117470131487;
  const double ve = 3.949659822581572E-03;

  q = ve / exp ( - de );

  ke[0] = ( int ) ( ( de / q ) * m2 );
  ke[1] = 0;

  we[0] = ( double ) ( q / m2 );
  we[255] = ( double ) ( de / m2 );

  fe[0] = 1.0;
  fe[255] = ( double ) ( exp ( - de ) );

  for ( i = 254; 1 <= i; i-- )
  {
    de = - log ( ve / de + exp ( - de ) );
    ke[i+1] = ( int ) ( ( de / te ) * m2 );
    te = de;
    fe[i] = ( double ) ( exp ( - de ) );
    we[i] = ( double ) ( de / m2 );
  }
  return;
}
/******************************************************************************/

double r4_nor ( unsigned long int *jsr, int kn[128], double fn[128], 
  double wn[128] )

/******************************************************************************/
/*
  Purpose:

    R4_NOR returns a normally distributed single precision real value.

  Discussion:

    The value returned is generated from a distribution with mean 0 and 
    variance 1.

    The underlying algorithm is the ziggurat method.

    Before the first call to this function, the user must call R4_NOR_SETUP
    to determine the values of KN, FN and WN.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 December 2008

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

  Parameters:

    Input/output, unsigned long int *JSR, the seed.

    Input, int KN[128], data computed by R4_NOR_SETUP.

    Input, double FN[128], WN[128], data computed by R4_NOR_SETUP.

    Output, double R4_NOR, a normally distributed random value.
*/
{
  int hz;
  int iz;
  const double r = 3.442620;
  double value;
  double x;
//  double y;

  hz = shr3 ( jsr );
  iz = ( hz & 127 );

  if ( abs ( hz ) < kn[iz] )
  {
    value = ( double ) ( hz ) * wn[iz];
  }
  else
  {
    for ( ; ; )
    {
      if ( iz == 0 )
      {
//        for ( ; ; )
//        {
          x = - 0.2904764 * log ( r4_uni ( jsr ) );
     //     y = - log ( r4_uni ( jsr ) );
     //     if ( x * x <= y + y );  // Martin B. doesn't understand this
     //     {
//            break;
     //     }
//        }

        if ( hz <= 0 )
        {
          value = - r - x;
        }
        else
        {
          value = + r + x;
        }
        break;
      }

      x = ( double ) ( hz ) * wn[iz];

      if ( fn[iz] + r4_uni ( jsr ) * ( fn[iz-1] - fn[iz] ) < exp ( - 0.5 * x * x ) )
      {
        value = x;
        break;
      }

      hz = shr3 ( jsr );
      iz = ( hz & 127 );

      if ( abs ( hz ) < kn[iz] )
      {
        value = ( double ) ( hz ) * wn[iz];
        break;
      }
    }
  }

  return value;
}
/******************************************************************************/

void r4_nor_setup ( int kn[128], double fn[128], double wn[128] )

/******************************************************************************/
/*
  Purpose:

    R4_NOR_SETUP sets data needed by R4_NOR.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 December 2008

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

  Parameters:

    Output, int KN[128], data needed by R4_NOR.

    Output, double FN[128], WN[128], data needed by R4_NOR.
*/
{
  double dn = 3.442619855899;
  int i;
  const double m1 = 2147483648.0;
  double q;
  double tn = 3.442619855899;
  const double vn = 9.91256303526217E-03;

  q = vn / exp ( - 0.5 * dn * dn );

  kn[0] = ( int ) ( ( dn / q ) * m1 );
  kn[1] = 0;

  wn[0] = ( double ) ( q / m1 );
  wn[127] = ( double ) ( dn / m1 );

  fn[0] = 1.0;
  fn[127] = ( double ) ( exp ( - 0.5 * dn * dn ) );

  for ( i = 126; 1 <= i; i-- )
  {
    dn = sqrt ( - 2.0 * log ( vn / dn + exp ( - 0.5 * dn * dn ) ) );
    kn[i+1] = ( int ) ( ( dn / tn ) * m1 );
    tn = dn;
    fn[i] = ( double ) ( exp ( - 0.5 * dn * dn ) );
    wn[i] = ( double ) ( dn / m1 );
  }
  return;
}
/******************************************************************************/

double r4_uni ( unsigned long int *jsr )

/******************************************************************************/
/*
  Purpose:

    R4_UNI returns a uniformly distributed real value.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 December 2008

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

  Parameters:

    Input/output, unsigned long int *JSR, the seed.

    Output, double R4_UNI, a uniformly distributed random value in
    the range [0,1].
*/
{
  unsigned long int jsr_input;
  double value;

  jsr_input = *jsr;

  *jsr = ( *jsr ^ ( *jsr <<   13 ) );
  *jsr = ( *jsr ^ ( *jsr >>   17 ) );
  *jsr = ( *jsr ^ ( *jsr <<    5 ) );

  value = (double)fmod ( 0.5 + ( double ) ( jsr_input + *jsr ) / 65536.0 / 65536.0, 1.0 );

  return value;
}
/******************************************************************************/

unsigned long int shr3 ( unsigned long int *jsr )

/******************************************************************************/
/*
  Purpose:

    SHR3 evaluates the SHR3 generator for integers.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 December 2008

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

  Parameters:

    Input/output, unsigned long int *JSR, the seed, which is updated 
    on each call.

    Output, unsigned long int SHR3, the new value.
*/
{
  unsigned long int value;

  value = *jsr;

  *jsr = ( *jsr ^ ( *jsr <<   13 ) );
  *jsr = ( *jsr ^ ( *jsr >>   17 ) );
  *jsr = ( *jsr ^ ( *jsr <<    5 ) );

  value = value + *jsr;

  return value;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  UNUSED(len); // Avoid "unused parameter" warning

  mprintf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
