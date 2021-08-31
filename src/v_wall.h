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

/** \file wall.hh
 * \brief Header file for the derived wall classes. */


#ifndef VOROPP_WALL_HH
#define VOROPP_WALL_HH


 // This must always be the first include directive
#include "config.h"

#include "v_cell.h"
#include "v_container.h"


/** \brief A class representing a spherical wall object.
 *
 * This class represents a spherical wall object. */
struct wall_sphere : public wall {
	public:
		/** Constructs a spherical wall object.
		 * \param[in] w_id_ an ID number to associate with the wall for
		 *		    neighbor tracking.
		 * \param[in] (xc_,yc_,zc_) a position vector for the sphere's
		 * 			    center.
		 * \param[in] rc_ the radius of the sphere. */
		wall_sphere(double xc_,double yc_,double zc_,double rc_,int w_id_=-99)
			: w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), rc(rc_) {}
		bool point_inside(double x,double y,double z);
		bool cut_cell_base(voronoicell_neighbor &c,double x,double y,double z);
		bool cut_cell(voronoicell_neighbor &c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
	private:
		const int w_id;
		const double xc,yc,zc,rc;
};

/** \brief A class representing a plane wall object.
 *
 * This class represents a single plane wall object. */
struct wall_plane : public wall {
	public:
		/** Constructs a plane wall object.
		 * \param[in] (xc_,yc_,zc_) a normal vector to the plane.
		 * \param[in] ac_ a displacement along the normal vector.
		 * \param[in] w_id_ an ID number to associate with the wall for
		 *		    neighbor tracking. */
		wall_plane(double xc_,double yc_,double zc_,double ac_,int w_id_=-99)
			: w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), ac(ac_) {}
		bool point_inside(double x,double y,double z);
		bool cut_cell_base(voronoicell_neighbor &c,double x,double y,double z);
		bool cut_cell(voronoicell_neighbor &c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
	private:
		const int w_id;
		const double xc,yc,zc,ac;
};

/** \brief A class representing a cylindrical wall object.
 *
 * This class represents a open cylinder wall object. */
struct wall_cylinder : public wall {
	public:
		/** Constructs a cylinder wall object.
		 * \param[in] (xc_,yc_,zc_) a point on the axis of the
		 *			    cylinder.
		 * \param[in] (xa_,ya_,za_) a vector pointing along the
		 *			    direction of the cylinder.
		 * \param[in] rc_ the radius of the cylinder
		 * \param[in] w_id_ an ID number to associate with the wall for
		 *		    neighbor tracking. */
		wall_cylinder(double xc_,double yc_,double zc_,double xa_,double ya_,double za_,double rc_,int w_id_=-99)
			: w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), xa(xa_), ya(ya_), za(za_),
			asi(1/(xa_*xa_+ya_*ya_+za_*za_)), rc(rc_) {}
		bool point_inside(double x,double y,double z);
		bool cut_cell_base(voronoicell_neighbor &c,double x,double y,double z);
		bool cut_cell(voronoicell_neighbor &c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
	private:
		const int w_id;
		const double xc,yc,zc,xa,ya,za,asi,rc;
};


/** \brief A class representing a conical wall object.
 *
 * This class represents a cone wall object. */
struct wall_cone : public wall {
	public:
		/** Constructs a cone wall object.
		 * \param[in] (xc_,yc_,zc_) the apex of the cone.
		 * \param[in] (xa_,ya_,za_) a vector pointing along the axis of
		 *			    the cone.
		 * \param[in] ang the angle (in radians) of the cone, measured
		 *		  from the axis.
		 * \param[in] w_id_ an ID number to associate with the wall for
		 *		    neighbor tracking. */
		wall_cone(double xc_,double yc_,double zc_,double xa_,double ya_,double za_,double ang,int w_id_=-99)
			: w_id(w_id_), xc(xc_), yc(yc_), zc(zc_), xa(xa_), ya(ya_), za(za_),
			asi(1/(xa_*xa_+ya_*ya_+za_*za_)),
			gra(tan(ang)), sang(sin(ang)), cang(cos(ang)) {}
		bool point_inside(double x,double y,double z);
		bool cut_cell_base(voronoicell_neighbor &c,double x,double y,double z);
		bool cut_cell(voronoicell_neighbor &c,double x,double y,double z) {return cut_cell_base(c,x,y,z);}
	private:
		const int w_id;
		const double xc,yc,zc,xa,ya,za,asi,gra,sang,cang;
};

#endif

