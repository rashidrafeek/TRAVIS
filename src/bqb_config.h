/***********************************************************************************

    LibBQB - File Format and Compression Algorithms for Trajectories of
             Volumetric Data and Atom Positions

    https://brehm-research.de/bqb

    Free software, licensed under GNU LGPL v3

    Copyright (c) Martin Brehm and Martin Thomas,
                  Martin Luther University Halle-Wittenberg, Germany,
                  2016 - 2021.

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



// This file is always included before any other include file.


#ifndef BQB_CONFIG_H
#define BQB_CONFIG_H


// Build this library together with TRAVIS?
#define BQB_INSIDE_TRAVIS


// Options for stand-alone build
#ifndef BQB_INSIDE_TRAVIS


	// Last change to this version of the source code
	#define BQB_SOURCE_VERSION "Dec 25 2020"


	// Use faster hard-coded string to floating point conversion (atof).
	// This speeds up reading CUBE files by a factor of > 2, and also other file formats.
	// Seems to give identical results to system atof(), but comment out if problems occur.
	#define BQB_FAST_ATOF




	#define _FILE_OFFSET_BITS 64


	#ifdef _MSC_VER
	  #pragma warning(disable:4786) // Warning "Debug Info truncated to 255 characters"
	  #pragma warning(disable:4702) // Warning "Unreachable Code"
	  #define _CRT_SECURE_NO_WARNINGS
	  #define _CRT_SECURE_NO_DEPRECATE
	#endif

	// Revision Information

	#include <stdio.h>

	#define GET_REVISION_INFO(X,LEN)  { if (LEN == 0) sprintf( X, "%s", "\"" __FILE__ "\"" ); else sprintf( X, "File %-*s compiled at %s %s, source version %s.", LEN, "\"" __FILE__ "\"", __DATE__,  __TIME__, BQB_SOURCE_VERSION ); }

	#define GET_SOURCE_VERSION(X)  sprintf( X, "%s", BQB_SOURCE_VERSION )

#else // BQB_INSIDE_TRAVIS

	// Include the TRAVIS main config file
	#include "config.h"

	#define BQB_SOURCE_VERSION SOURCE_VERSION

	#ifdef FAST_ATOF
		#define BQB_FAST_ATOF
	#endif

#endif


const char *GetRevisionInfo_bqb_alphabet(unsigned int len);
const char *GetSourceVersion_bqb_alphabet();
const char *GetRevisionInfo_bqb_bitset(unsigned int len);
const char *GetSourceVersion_bqb_bitset();
const char *GetRevisionInfo_bqb_crc(unsigned int len);
const char *GetSourceVersion_bqb_crc();
const char *GetRevisionInfo_bqb_cubeframe(unsigned int len);
const char *GetSourceVersion_bqb_cubeframe();
const char *GetRevisionInfo_bqb_driver(unsigned int len);
const char *GetSourceVersion_bqb_driver();
const char *GetRevisionInfo_bqb_engine(unsigned int len);
const char *GetSourceVersion_bqb_engine();
const char *GetRevisionInfo_bqb_extrapolator(unsigned int len);
const char *GetSourceVersion_bqb_extrapolator();
const char *GetRevisionInfo_bqb_format(unsigned int len);
const char *GetSourceVersion_bqb_format();
const char *GetRevisionInfo_bqb_hilbert(unsigned int len);
const char *GetSourceVersion_bqb_hilbert();
const char *GetRevisionInfo_bqb_hufftree(unsigned int len);
const char *GetSourceVersion_bqb_hufftree();
const char *GetRevisionInfo_bqb_integerengine(unsigned int len);
const char *GetSourceVersion_bqb_integerengine();
const char *GetRevisionInfo_bqb_interface(unsigned int len);
const char *GetSourceVersion_bqb_interface();
const char *GetRevisionInfo_bqb_largeinteger(unsigned int len);
const char *GetSourceVersion_bqb_largeinteger();
const char *GetRevisionInfo_bqb_linalg(unsigned int len);
const char *GetSourceVersion_bqb_linalg();
const char *GetRevisionInfo_bqb_math(unsigned int len);
const char *GetSourceVersion_bqb_math();
const char *GetRevisionInfo_bqb_parmset(unsigned int len);
const char *GetSourceVersion_bqb_parmset();
const char *GetRevisionInfo_bqb_tools(unsigned int len);
const char *GetSourceVersion_bqb_tools();
const char *GetRevisionInfo_bqbtool(unsigned int len);
const char *GetSourceVersion_bqbtool();



#endif


