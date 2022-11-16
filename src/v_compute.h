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


// Voro++, a 3D cell-based Voronoi library
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

/** \file v_compute.hh
 * \brief Header file for the voro_compute template and related classes. */


#ifndef VOROPP_V_COMPUTE_HH
#define VOROPP_V_COMPUTE_HH


 // This must always be the first include directive
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include "v_config.h"
#include "v_worklist.h"
#include "v_cell.h"


class container_periodic_poly;


/** \brief Structure for holding information about a particle.
 *
 * This small structure holds information about a single particle, and is used
 * by several of the routines in the voro_compute template for passing
 * information by reference between functions. */
struct particle_record {
	/** The index of the block that the particle is within. */
	int ijk;
	/** The number of particle within its block. */
	int l;
	/** The x-index of the block. */
	int di;
	/** The y-index of the block. */
	int dj;
	/** The z-index of the block. */
	int dk;
};

/** \brief Template for carrying out Voronoi cell computations. */
class voro_compute {
	public:
		/** A reference to the container class on which to carry out*/
		container_periodic_poly &con;
		/** The size of an internal computational block in the x
		 * direction. */
		const double boxx;
		/** The size of an internal computational block in the y
		 * direction. */
		const double boxy;
		/** The size of an internal computational block in the z
		 * direction. */
		const double boxz;
		/** The inverse box length in the x direction, set to
		 * nx/(bx-ax). */
		const double xsp;
		/** The inverse box length in the y direction, set to
		 * ny/(by-ay). */
		const double ysp;
		/** The inverse box length in the z direction, set to
		 * nz/(bz-az). */
		const double zsp;
		/** The number of boxes in the x direction for the searching mask. */
		const int hx;
		/** The number of boxes in the y direction for the searching mask. */
		const int hy;
		/** The number of boxes in the z direction for the searching mask. */
		const int hz;
		/** A constant, set to the value of hx multiplied by hy, which
		 * is used in the routines which step through mask boxes in
		 * sequence. */
		const int hxy;
		/** A constant, set to the value of hx*hy*hz, which is used in
		 * the routines which step through mask boxes in sequence. */
		const int hxyz;
		/** The number of floating point entries to store for each
		 * particle. */
		const int ps;
		/** This array holds the numerical IDs of each particle in each
		 * computational box. */
		int **id;
		/** A two dimensional array holding particle positions. For the
		 * derived container_poly class, this also holds particle
		 * radii. */
		double **p;
		/** An array holding the number of particles within each
		 * computational box of the container. */
		int *co;
		voro_compute(container_periodic_poly &con_,int hx_,int hy_,int hz_);
		/** The class destructor frees the dynamically allocated memory
		 * for the mask and queue. */
		~voro_compute() {
			delete [] qu;
			delete [] mask;
		}
		bool compute_cell(voronoicell_neighbor &c,int ijk,int s,int ci,int cj,int ck);
		void find_voronoi_cell(double x,double y,double z,int ci,int cj,int ck,int ijk,particle_record &w,double &mrs);
	private:
		/** A constant set to boxx*boxx+boxy*boxy+boxz*boxz, which is
		 * frequently used in the computation. */
		const double bxsq;
		/** This sets the current value being used to mark tested blocks
		 * in the mask. */
		unsigned int mv;
		/** The current size of the search list. */
		int qu_size;
		/** A pointer to the array of worklists. */
		const unsigned int *wl;
		/** An pointer to the array holding the minimum distances
		 * associated with the worklists. */
		double *mrad;
		/** This array is used during the cell computation to determine
		 * which blocks have been considered. */
		unsigned int *mask;
		/** An array is used to store the queue of blocks to test
		 * during the Voronoi cell computation. */
		int *qu;
		/** A pointer to the end of the queue array, used to determine
		 * when the queue is full. */
		int *qu_l;
		bool corner_test(voronoicell_neighbor &c,double xl,double yl,double zl,double xh,double yh,double zh);
		inline bool edge_x_test(voronoicell_neighbor &c,double x0,double yl,double zl,double x1,double yh,double zh);
		inline bool edge_y_test(voronoicell_neighbor &c,double xl,double y0,double zl,double xh,double y1,double zh);
		inline bool edge_z_test(voronoicell_neighbor &c,double xl,double yl,double z0,double xh,double yh,double z1);
		inline bool face_x_test(voronoicell_neighbor &c,double xl,double y0,double z0,double y1,double z1);
		inline bool face_y_test(voronoicell_neighbor &c,double x0,double yl,double z0,double x1,double z1);
		inline bool face_z_test(voronoicell_neighbor &c,double x0,double y0,double zl,double x1,double y1);
		bool compute_min_max_radius(int di,int dj,int dk,double fx,double fy,double fz,double gx,double gy,double gz,double& crs,double mrs);
		bool compute_min_radius(int di,int dj,int dk,double fx,double fy,double fz,double mrs);
		inline void add_to_mask(int ei,int ej,int ek,int *&qu_e);
		inline void scan_bits_mask_add(unsigned int q,unsigned int *mijk,int ei,int ej,int ek,int *&qu_e);
		inline void scan_all(int ijk,double x,double y,double z,int di,int dj,int dk,particle_record &w,double &mrs);
		void add_list_memory(int*& qu_s,int*& qu_e);
		/** Resets the mask in cases where the mask counter wraps
		 * around. */
		inline void reset_mask() {
			for(unsigned int *mp = (mask);mp<mask+hxyz;mp++) *mp=0;
		}
};

#endif

