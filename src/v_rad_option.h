/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

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

/** \file rad_option.hh
 * \brief Header file for the classes encapsulating functionality for the
 * regular and radical Voronoi tessellations. */


#ifndef VOROPP_RAD_OPTION_HH
#define VOROPP_RAD_OPTION_HH


 // This must always be the first include directive
#include "config.h"


/** \brief Class containing all of the routines that are specific to computing 
 * the regular Voronoi tessellation.
 *
 * The container and container_periodic classes are derived from this class,
 * and during the Voronoi cell computation, these routines are used to create
 * the regular Voronoi tessellation. */
class radius_mono {
	protected:
		/** This is called prior to computing a Voronoi cell for a
		 * given particle to initialize any required constants.
		 * \param[in] ijk the block that the particle is within.
		 * \param[in] s the index of the particle within the block. */
		inline void r_init(int ijk,int s) { UNUSED(ijk); UNUSED(s); /* Suppress "unused parameter" warning */ }
		/** Sets a required constant to be used when carrying out a
		 * plane bounds check. */
		inline void r_prime(double rv) { UNUSED(rv); /* Suppress "unused parameter" warning */ }
		/** Carries out a radius bounds check.
		 * \param[in] crs the radius squared to be tested.
		 * \param[in] mrs the current maximum distance to a Voronoi
		 *                vertex multiplied by two.
		 * \return True if particles at this radius could not possibly
		 * cut the cell, false otherwise. */
		inline bool r_ctest(double crs,double mrs) {return crs>mrs;}
		/** Scales a plane displacement during a plane bounds check.
		 * \param[in] lrs the plane displacement.
		 * \return The scaled value. */
		inline double r_cutoff(double lrs) {return lrs;}
		/** Adds the maximum radius squared to a given value.
		 * \param[in] rs the value to consider.
		 * \return The value with the radius squared added. */
		inline double r_max_add(double rs) {return rs;}
		/** Subtracts the radius squared of a particle from a given 
		 * value.
		 * \param[in] rs the value to consider.
		 * \param[in] ijk the block that the particle is within.
		 * \param[in] q the index of the particle within the block. 
		 * \return The value with the radius squared subtracted. */
		inline double r_current_sub(double rs,int ijk,int q) { UNUSED(ijk); UNUSED(q); /* Suppress "unused parameter" warning */ return rs;}
		/** Scales a plane displacement prior to use in the plane cutting
		 * algorithm.
		 * \param[in] rs the initial plane displacement.
		 * \param[in] ijk the block that the particle is within.
		 * \param[in] q the index of the particle within the block.
		 * \return The scaled plane displacement. */ 
		inline double r_scale(double rs,int ijk,int q) { UNUSED(ijk); UNUSED(q); /* Suppress "unused parameter" warning */ return rs;}
		/** Scales a plane displacement prior to use in the plane
		 * cutting algorithm, and also checks if it could possibly cut
		 * the cell.
		 * \param[in,out] rs the plane displacement to be scaled.
		 * \param[in] mrs the current maximum distance to a Voronoi
		 *                vertex multiplied by two.
		 * \param[in] ijk the block that the particle is within.
		 * \param[in] q the index of the particle within the block.
		 * \return True if the cell could possibly cut the cell, false
		 * otherwise. */		
		inline bool r_scale_check(double &rs,double mrs,int ijk,int q) { UNUSED(ijk); UNUSED(q); /* Suppress "unused parameter" warning */ return rs<mrs;}
};

/**  \brief Class containing all of the routines that are specific to computing 
 * the radical Voronoi tessellation.
 *
 * The container_poly and container_periodic_poly classes are derived from this
 * class, and during the Voronoi cell computation, these routines are used to
 * create the radical Voronoi tessellation. */
class radius_poly {
	public:
		/** A two-dimensional array holding particle positions and radii. */			
		double **ppr;
		/** The current maximum radius of any particle, used to
		 * determine when to cut off the radical Voronoi computation.
		 * */
		double max_radius;
		/** The class constructor sets the maximum particle radius to
		 * be zero. */
		radius_poly() : max_radius(0) {}
	protected:
		/** This is called prior to computing a Voronoi cell for a
		 * given particle to initialize any required constants.
		 * \param[in] ijk the block that the particle is within.
		 * \param[in] s the index of the particle within the block. */
		// Patched on 12.04.2021 to fix vectorization issue
		inline void r_init(int ijk,int s) {
			volatile double rsq = max_radius*max_radius;
			r_rad=ppr[ijk][4*s+3]*ppr[ijk][4*s+3];
			r_mul = r_rad - rsq;
		}
		/** Sets a required constant to be used when carrying out a
		 * plane bounds check. */
		inline void r_prime(double rv) {r_val=1+r_mul/rv;}
		/** Carries out a radius bounds check.
		 * \param[in] crs the radius squared to be tested.
		 * \param[in] mrs the current maximum distance to a Voronoi
		 *                vertex multiplied by two.
		 * \return True if particles at this radius could not possibly
		 * cut the cell, false otherwise. */		
		inline bool r_ctest(double crs,double mrs) {return crs+r_mul>sqrt(mrs*crs);}
		/** Scales a plane displacement during a plane bounds check.
		 * \param[in] lrs the plane displacement.
		 * \return The scaled value. */		
		inline double r_cutoff(double lrs) {return lrs*r_val;}
		/** Adds the maximum radius squared to a given value.
		 * \param[in] rs the value to consider.
		 * \return The value with the radius squared added. */		
		inline double r_max_add(double rs) {return rs+max_radius*max_radius;}
		/** Subtracts the radius squared of a particle from a given 
		 * value.
		 * \param[in] rs the value to consider.
		 * \param[in] ijk the block that the particle is within.
		 * \param[in] q the index of the particle within the block. 
		 * \return The value with the radius squared subtracted. */
		inline double r_current_sub(double rs,int ijk,int q) {
			return rs-ppr[ijk][4*q+3]*ppr[ijk][4*q+3];
		}
		/** Scales a plane displacement prior to use in the plane cutting
		 * algorithm.
		 * \param[in] rs the initial plane displacement.
		 * \param[in] ijk the block that the particle is within.
		 * \param[in] q the index of the particle within the block.
		 * \return The scaled plane displacement. */ 
		inline double r_scale(double rs,int ijk,int q) {
			return rs+r_rad-ppr[ijk][4*q+3]*ppr[ijk][4*q+3];
		}
		/** Scales a plane displacement prior to use in the plane
		 * cutting algorithm, and also checks if it could possibly cut
		 * the cell.
		 * \param[in,out] rs the plane displacement to be scaled.
		 * \param[in] mrs the current maximum distance to a Voronoi
		 *                vertex multiplied by two.
		 * \param[in] ijk the block that the particle is within.
		 * \param[in] q the index of the particle within the block.
		 * \return True if the cell could possibly cut the cell, false
		 * otherwise. */
		inline bool r_scale_check(double &rs,double mrs,int ijk,int q) {
			double trs=rs;
			rs+=r_rad-ppr[ijk][4*q+3]*ppr[ijk][4*q+3];
			return rs<sqrt(mrs*trs);
		}
	private:
		double r_rad,r_mul,r_val;
};

#endif


