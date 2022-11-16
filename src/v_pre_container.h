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

/** \file pre_container.hh
 * \brief Header file for the pre_container and related classes. */


#ifndef VOROPP_PRE_CONTAINER_HH
#define VOROPP_PRE_CONTAINER_HH


 // This must always be the first include directive
#include "config.h"

#include "v_c_loops.h"
#include "v_container.h"


/** \brief A class for storing an arbitrary number of particles, prior to setting
 * up a container geometry.
 *
 * The pre_container_base class can dynamically import and store an arbitrary
 * number of particles. Once the particles have been read in, an appropriate
 * container class can be set up with the optimal grid size, and the particles
 * can be transferred.
 *
 * The pre_container_base class is not intended for direct use, but forms the
 * base of the pre_container and pre_container_poly classes, that add routines
 * depending on whether particle radii need to be tracked or not. */
class pre_container_base {
	public:
		/** The minimum x coordinate of the container. */
		const double ax;
		/** The maximum x coordinate of the container. */
		const double bx;
		/** The minimum y coordinate of the container. */
		const double ay;
		/** The maximum y coordinate of the container. */
		const double by;
		/** The minimum z coordinate of the container. */
		const double az;
		/** The maximum z coordinate of the container. */
		const double bz;
		/** A boolean value that determines if the x coordinate in
		 * periodic or not. */
		const bool xperiodic;
		/** A boolean value that determines if the y coordinate in
		 * periodic or not. */
		const bool yperiodic;
		/** A boolean value that determines if the z coordinate in
		 * periodic or not. */
		const bool zperiodic;
		void guess_optimal(int &nx,int &ny,int &nz);
		pre_container_base(double ax_,double bx_,double ay_,double by_,double az_,double bz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int ps_);
		~pre_container_base();
		/** Calculates and returns the total number of particles stored
		 * within the class.
		 * \return The number of particles. */
		inline int total_particles() {
			return (int)((end_id-pre_id)*pre_container_chunk_size+(ch_id-*end_id));
		}
	protected:
		/** The number of doubles associated with a single particle
		 * (three for the standard container, four when radius
		 * information is stored). */
		const int ps;
		void new_chunk();
		void extend_chunk_index();
		/** The size of the chunk index. */
		int index_sz;
		/** A pointer to the chunk index to store the integer particle
		 * IDs. */
		int **pre_id;
		/** A pointer to the last allocated integer ID chunk. */
		int **end_id;
		/** A pointer to the end of the integer ID chunk index, used to
		 * determine when the chunk index is full. */
		int **l_id;
		/** A pointer to the next available slot on the current
		 * particle ID chunk. */
		int *ch_id;
		/** A pointer to the end of the current integer chunk. */
		int *e_id;
		/** A pointer to the chunk index to store the floating point
		 * information associated with particles. */
		double **pre_p;
		/** A pointer to the last allocated chunk of floating point
		 * information. */
		double **end_p;
		/** A pointer to the next available slot on the current
		 * floating point chunk. */
		double *ch_p;
};


/** \brief A class for storing an arbitrary number of particles with radius
 * information, prior to setting up a container geometry.
 *
 * The pre_container_poly class is an extension of the pre_container_base class
 * for cases when particle radius information is available. */
class pre_container_poly : public pre_container_base {
	public:
		/** The class constructor sets up the geometry of container,
		 * initializing the minimum and maximum coordinates in each
		 * direction.
		 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
		 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
		 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
		 * \param[in] (xperiodic_,yperiodic_,zperiodic_ ) flags setting whether the
		 *                                                container is periodic in
		 *                                                each coordinate direction. */
		pre_container_poly(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
				bool xperiodic_,bool yperiodic_,bool zperiodic_)
			: pre_container_base(ax_,bx_,ay_,by_,az_,bz_,xperiodic_,yperiodic_,zperiodic_,4) {};
		void put(int n,double x,double y,double z,double r);
		void import(FILE *fp=stdin);
		/** Imports particles from a file.
		 * \param[in] filename the name of the file to read from. */
		inline void import(const char* filename) {
			FILE *fp=safe_fopen(filename,"r");
			import(fp);
			fclose(fp);
		}
		void setup(container_periodic_poly &con);
//		void setup(particle_order &vo,container_periodic_poly &con);
};


#endif


