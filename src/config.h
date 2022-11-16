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


// This file is always included before any other include file.


#ifndef CONFIG_H
#define CONFIG_H


// Last change to this version of the source code
#define SOURCE_VERSION "Jul 29 2022"

/* Please uncomment / comment out the flags you want to use / not to use. */

// You have to change this according to your target operating system.
// For a generic platform-independent build, comment out both lines.
//#define TARGET_WINDOWS    // Tested with Microsoft Visual Studio
#define TARGET_LINUX      // Tested with GNU C++ compiler on GNU/Linux

// Enable to switch on minimization of charge variance via two-step Voronoi tessellation.
// Requires the "libtravis" sub-directory and a compiler with C++14 standard.
#define NEW_CHARGEVAR

//>>>OMP
// Activate to use OpenMPI
//#define USE_OMP
//<<<OMP

// Use color for screen output?
#define USE_COLOR

// Maximum number of bonds any atom may form
#define MAX_BONDS 8

// Activate if you have OpenBabel in the system search path
//#define USE_OPENBABEL

// Use the FFTW library? Otherwise, the built-in KISS-FFT routine is used (default)
//#define USE_FFTW

// Uncomment the following line for an official release
//#define RELEASE_VERSION "1.14.0"

// Activate to enable reading NetCDF-trajectories (AMBER9+)
// It is necessary to have the NetCDF-libraries installed (e.g. libnetcdff-dev and libnetcdff-dev)
// compile with: make netcdf
//#define USE_NETCDF

// Bound checking in dynamic arrays?
#define DEBUG_ARRAYS


// Handle volumetric data in single or double precision? Active=double, inactive=single
#define VORI_DOUBLE


// Use faster hard-coded string to floating point conversion (atof).
// This speeds up reading CUBE files by a factor of > 2, and also other file formats.
// Seems to give identical results to system atof(), but comment out if problems occur.
#define FAST_ATOF


#ifdef VORI_DOUBLE
	#define VORI_FLOAT double
#else
	#define VORI_FLOAT float
#endif


/***************************************************************************/
// Some "hardcore" debug flags
// Warning: They may drastically decrease performance
/***************************************************************************/


// Activate debug backtrace?
//#define DEBUG_BACKTRACE
//#define DEBUG_EXTENDED_BACKTRACE


//#define DEBUG_DATABASE
//#define DEBUG_CSTRING
//#define DEBUG_MATULTRA
//#define DEBUG_CVECTOR3
//#define DEBUG_CDVECTOR3


//#define DEBUG_COBARRAY
//#define DEBUG_CPTRARRAY
//#define DEBUG_CWORDARRAY
//#define DEBUG_CLONGARRAY
//#define DEBUG_CINTARRAY
//#define DEBUG_CFLOATARRAY
//#define DEBUG_CDOUBLEARRAY
//#define DEBUG_CVEC3ARRAY
//#define DEBUG_CDVEC3ARRAY


#define _FILE_OFFSET_BITS 64


#ifdef TARGET_WINDOWS
	#pragma warning(disable:4503) // Warning "Name length exceeded"
	#pragma warning(disable:4786) // Warning "Debug Info truncated to 255 characters"
	#pragma warning(disable:4702) // Warning "Unreachable Code"
	#define _CRT_SECURE_NO_WARNINGS
	#define _CRT_SECURE_NO_DEPRECATE
#endif


#ifndef __PRETTY_FUNCTION__
	#ifdef __func__
		#define __PRETTY_FUNCTION__ __func__
	#elif defined(__FUNCTION__)
		#define __PRETTY_FUNCTION__ __FUNCTION__
	#else
		#define __PRETTY_FUNCTION__ ""
	#endif
#endif


// Revision Information

#include <stdio.h>


#define GET_REVISION_INFO(X,LEN)  { if (LEN == 0) sprintf( X, "%s", "\"" __FILE__ "\"" ); else sprintf( X, "File %-*s compiled at %s %s, source version %s.", LEN, "\"" __FILE__ "\"", __DATE__,  __TIME__, SOURCE_VERSION ); }

#define GET_SOURCE_VERSION(X)  sprintf( X, "%s", SOURCE_VERSION )


const char *GetRevisionInfo_2df(unsigned int len);
const char *GetSourceVersion_2df();
const char *GetRevisionInfo_acf(unsigned int len);
const char *GetSourceVersion_acf();
const char *GetRevisionInfo_aggrtopo(unsigned int len);
const char *GetSourceVersion_aggrtopo();
const char *GetRevisionInfo_analysisgroup(unsigned int len);
const char *GetSourceVersion_analysisgroup();
const char *GetRevisionInfo_asciiart(unsigned int len);
const char *GetSourceVersion_asciiart();
const char *GetRevisionInfo_atomgroup(unsigned int len);
const char *GetSourceVersion_atomgroup();
const char *GetRevisionInfo_backtrace(unsigned int len);
const char *GetSourceVersion_backtrace();
const char *GetRevisionInfo_base64(unsigned int len);
const char *GetSourceVersion_base64();
const char *GetRevisionInfo_bicgstab(unsigned int len);
const char *GetSourceVersion_bicgstab();
const char *GetRevisionInfo_bintools(unsigned int len);
const char *GetSourceVersion_bintools();
const char *GetRevisionInfo_bintree(unsigned int len);
const char *GetSourceVersion_bintree();

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


const char *GetRevisionInfo_chargevar2(unsigned int len);
const char *GetSourceVersion_chargevar2();
const char *GetRevisionInfo_chiral(unsigned int len);
const char *GetSourceVersion_chiral();
const char *GetRevisionInfo_cluster(unsigned int len);
const char *GetSourceVersion_cluster();
const char *GetRevisionInfo_contactmatrix(unsigned int len);
const char *GetSourceVersion_contactmatrix();
const char *GetRevisionInfo_cubetool(unsigned int len);
const char *GetSourceVersion_cubetool();
const char *GetRevisionInfo_dacf(unsigned int len);
const char *GetSourceVersion_dacf();
const char *GetRevisionInfo_database(unsigned int len);
const char *GetSourceVersion_database();
const char *GetRevisionInfo_df(unsigned int len);
const char *GetSourceVersion_df();
const char *GetRevisionInfo_diagonal(unsigned int len);
const char *GetSourceVersion_diagonal();
const char *GetRevisionInfo_domain(unsigned int len);
const char *GetSourceVersion_domain();
const char *GetRevisionInfo_dpol(unsigned int len);
const char *GetSourceVersion_dpol();
const char *GetRevisionInfo_elementdata(unsigned int len);
const char *GetSourceVersion_elementdata();
const char *GetRevisionInfo_fdf(unsigned int len);
const char *GetSourceVersion_fdf();
const char *GetRevisionInfo_fft(unsigned int len);
const char *GetSourceVersion_fft();
const char *GetRevisionInfo_fixplproj(unsigned int len);
const char *GetSourceVersion_fixplproj();
const char *GetRevisionInfo_gather(unsigned int len);
const char *GetSourceVersion_gather();
const char *GetRevisionInfo_geodens(unsigned int len);
const char *GetSourceVersion_geodens();
const char *GetRevisionInfo_globalvar(unsigned int len);
const char *GetSourceVersion_globalvar();
const char *GetRevisionInfo_grace(unsigned int len);
const char *GetSourceVersion_grace();
const char *GetRevisionInfo_hbond(unsigned int len);
const char *GetSourceVersion_hbond();
const char *GetRevisionInfo_interface(unsigned int len);
const char *GetSourceVersion_interface();
const char *GetRevisionInfo_internalcoord(unsigned int len);
const char *GetSourceVersion_internalcoord();
const char *GetRevisionInfo_ionpair(unsigned int len);
const char *GetSourceVersion_ionpair();
const char *GetRevisionInfo_ir(unsigned int len);
const char *GetSourceVersion_ir();
const char *GetRevisionInfo_kiss_fft(unsigned int len);
const char *GetSourceVersion_kiss_fft();
const char *GetRevisionInfo_largeinteger(unsigned int len);
const char *GetSourceVersion_largeinteger();
const char *GetRevisionInfo_linalg(unsigned int len);
const char *GetSourceVersion_linalg();
const char *GetRevisionInfo_lmmin(unsigned int len);
const char *GetSourceVersion_lmmin();
const char *GetRevisionInfo_lmwrapper(unsigned int len);
const char *GetSourceVersion_lmwrapper();
const char *GetRevisionInfo_lu_decomp(unsigned int len);
const char *GetSourceVersion_lu_decomp();
const char *GetRevisionInfo_luzar(unsigned int len);
const char *GetSourceVersion_luzar();
const char *GetRevisionInfo_maintools(unsigned int len);
const char *GetSourceVersion_maintools();
const char *GetRevisionInfo_matrixplot(unsigned int len);
const char *GetSourceVersion_matrixplot();
const char *GetRevisionInfo_moltools(unsigned int len);
const char *GetSourceVersion_moltools();
const char *GetRevisionInfo_nbexchange(unsigned int len);
const char *GetSourceVersion_nbexchange();
const char *GetRevisionInfo_nbsearch(unsigned int len);
const char *GetSourceVersion_nbsearch();
const char *GetRevisionInfo_normalcoordinate(unsigned int len);
const char *GetSourceVersion_normalcoordinate();
const char *GetRevisionInfo_normalmode(unsigned int len);
const char *GetSourceVersion_normalmode();
const char *GetRevisionInfo_order(unsigned int len);
const char *GetSourceVersion_order();
const char *GetRevisionInfo_order_chain(unsigned int len);
const char *GetSourceVersion_order_chain();
const char *GetRevisionInfo_order_vector(unsigned int len);
const char *GetSourceVersion_order_vector();
const char *GetRevisionInfo_pdf(unsigned int len);
const char *GetSourceVersion_pdf();
const char *GetRevisionInfo_plproj(unsigned int len);
const char *GetSourceVersion_plproj();
const char *GetRevisionInfo_posdomain(unsigned int len);
const char *GetSourceVersion_posdomain();
const char *GetRevisionInfo_qr_fact(unsigned int len);
const char *GetSourceVersion_qr_fact();
const char *GetRevisionInfo_raman(unsigned int len);
const char *GetSourceVersion_raman();
const char *GetRevisionInfo_random(unsigned int len);
const char *GetSourceVersion_random();
const char *GetRevisionInfo_reflux(unsigned int len);
const char *GetSourceVersion_reflux();
const char *GetRevisionInfo_region(unsigned int len);
const char *GetSourceVersion_region();
const char *GetRevisionInfo_reordyn(unsigned int len);
const char *GetSourceVersion_reordyn();
const char *GetRevisionInfo_revsdf(unsigned int len);
const char *GetSourceVersion_revsdf();
const char *GetRevisionInfo_roa(unsigned int len);
const char *GetSourceVersion_roa();
const char *GetRevisionInfo_sankey(unsigned int len);
const char *GetSourceVersion_sankey();
const char *GetRevisionInfo_sdfmap(unsigned int len);
const char *GetSourceVersion_sdfmap();
const char *GetRevisionInfo_spectrum(unsigned int len);
const char *GetSourceVersion_spectrum();
const char *GetRevisionInfo_statistics(unsigned int len);
const char *GetSourceVersion_statistics();
const char *GetRevisionInfo_structurefactor(unsigned int len);
const char *GetSourceVersion_structurefactor();
const char *GetRevisionInfo_svgwriter(unsigned int len);
const char *GetSourceVersion_svgwriter();
const char *GetRevisionInfo_tddf(unsigned int len);
const char *GetSourceVersion_tddf();
const char *GetRevisionInfo_tensor(unsigned int len);
const char *GetSourceVersion_tensor();
const char *GetRevisionInfo_tetrapak(unsigned int len);
const char *GetSourceVersion_tetrapak();
const char *GetRevisionInfo_timestep(unsigned int len);
const char *GetSourceVersion_timestep();
const char *GetRevisionInfo_tools(unsigned int len);
const char *GetSourceVersion_tools();
const char *GetRevisionInfo_travis(unsigned int len);
const char *GetSourceVersion_travis();
const char *GetRevisionInfo_void(unsigned int len);
const char *GetSourceVersion_void();
const char *GetRevisionInfo_v_base(unsigned int len);
const char *GetSourceVersion_v_base();
const char *GetRevisionInfo_v_base_wl(unsigned int len);
const char *GetSourceVersion_v_base_wl();
const char *GetRevisionInfo_v_cell(unsigned int len);
const char *GetSourceVersion_v_cell();
const char *GetRevisionInfo_v_common(unsigned int len);
const char *GetSourceVersion_v_common();
const char *GetRevisionInfo_v_compute(unsigned int len);
const char *GetSourceVersion_v_compute();
const char *GetRevisionInfo_v_container(unsigned int len);
const char *GetSourceVersion_v_container();
const char *GetRevisionInfo_v_container_prd(unsigned int len);
const char *GetSourceVersion_v_container_prd();
const char *GetRevisionInfo_v_pre_container(unsigned int len);
const char *GetSourceVersion_v_pre_container();
const char *GetRevisionInfo_v_unitcell(unsigned int len);
const char *GetSourceVersion_v_unitcell();
const char *GetRevisionInfo_v_wall(unsigned int len);
const char *GetSourceVersion_v_wall();
const char *GetRevisionInfo_vorowrapper(unsigned int len);
const char *GetSourceVersion_vorowrapper();
const char *GetRevisionInfo_xbytearray(unsigned int len);
const char *GetSourceVersion_xbytearray();
const char *GetRevisionInfo_xdf(unsigned int len);
const char *GetSourceVersion_xdf();
const char *GetRevisionInfo_xdmatrix3(unsigned int len);
const char *GetSourceVersion_xdmatrix3();
const char *GetRevisionInfo_xdmatrixmn(unsigned int len);
const char *GetSourceVersion_xdmatrixmn();
const char *GetRevisionInfo_xdoublearray(unsigned int len);
const char *GetSourceVersion_xdoublearray();
const char *GetRevisionInfo_xdvec3array(unsigned int len);
const char *GetSourceVersion_xdvec3array();
const char *GetRevisionInfo_xdvector3(unsigned int len);
const char *GetSourceVersion_xdvector3();
const char *GetRevisionInfo_xdvectorn(unsigned int len);
const char *GetSourceVersion_xdvectorn();
const char *GetRevisionInfo_xintarray(unsigned int len);
const char *GetSourceVersion_xintarray();
const char *GetRevisionInfo_xlongarray(unsigned int len);
const char *GetSourceVersion_xlongarray();
const char *GetRevisionInfo_xmemfile(unsigned int len);
const char *GetSourceVersion_xmemfile();
const char *GetRevisionInfo_xobarray(unsigned int len);
const char *GetSourceVersion_xobarray();
const char *GetRevisionInfo_xptrarray(unsigned int len);
const char *GetSourceVersion_xptrarray();
const char *GetRevisionInfo_string(unsigned int len);
const char *GetSourceVersion_string();
const char *GetRevisionInfo_xwordarray(unsigned int len);
const char *GetSourceVersion_xwordarray();
const char *GetRevisionInfo_ziggurat(unsigned int len);
const char *GetSourceVersion_ziggurat();





#endif

