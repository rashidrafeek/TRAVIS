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

/** \file unitcell.hh
 * \brief Header file for the unitcell class. */


#ifndef VOROPP_UNITCELL_HH
#define VOROPP_UNITCELL_HH


// This must always be the first include directive
#include "config.h"

#include <vector>

#include "v_config.h"
#include "v_cell.h"


/** \brief Class for computation of the unit Voronoi cell associated with
 * a 3D non-rectangular periodic domain. */
class unitcell {
	public:
		/** The x coordinate of the first vector defining the periodic
		 * domain. */
		const double bx;
		/** The x coordinate of the second vector defining the periodic
		 * domain. */
		const double bxy;
		/** The y coordinate of the second vector defining the periodic
		 * domain. */
		const double by;
		/** The x coordinate of the third vector defining the periodic
		 * domain. */
		const double bxz;
		/** The y coordinate of the third vector defining the periodic
		 * domain. */
		const double byz;
		/** The z coordinate of the third vector defining the periodic
		 * domain. */
		const double bz;
		/** The computed unit Voronoi cell corresponding the given
		 * 3D non-rectangular periodic domain geometry. */
		voronoicell_neighbor unit_voro;
		unitcell(double bx_,double bxy_,double by_,double bxz_,double byz_,double bz_);
		/** Draws an outline of the domain in Gnuplot format.
		 * \param[in] filename the filename to write to. */
		inline void draw_domain_gnuplot(const char* filename) {
			FILE *fp = (safe_fopen(filename,"w"));
			draw_domain_gnuplot(fp);
			fclose(fp);
		}
		void draw_domain_gnuplot(FILE *fp=stdout);
		/** Draws an outline of the domain in Gnuplot format.
		 * \param[in] filename the filename to write to. */
		inline void draw_domain_pov(const char* filename) {
			FILE *fp = (safe_fopen(filename,"w"));
			draw_domain_pov(fp);
			fclose(fp);
		}
		void draw_domain_pov(FILE *fp=stdout);
		bool intersects_image(double dx,double dy,double dz,double &vol);
		void images(std::vector<int> &vi,std::vector<double> &vd);
	protected:
		/** The maximum y-coordinate that could possibly cut the
		 * computed unit Voronoi cell. */
		double max_uv_y;
		/** The maximum z-coordinate that could possibly cut the
		 * computed unit Voronoi cell. */
		double max_uv_z;
	private:
		inline void unit_voro_apply(int i,int j,int k);
		bool unit_voro_intersect(int l);
		inline bool unit_voro_test(int i,int j,int k);
};

#endif


