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

/** \file wall.cc
 * \brief Function implementations for the derived wall classes. */


 // This must always be the first include directive
#include "config.h"

#include "v_wall.h"


const char *GetRevisionInfo_v_wall(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_v_wall() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


/** Tests to see whether a point is inside the sphere wall object.
 * \param[in,out] (x,y,z) the vector to test.
 * \return True if the point is inside, false if the point is outside. */
bool wall_sphere::point_inside(double x,double y,double z) {
	return (x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc)<rc*rc;
}

/** Cuts a cell by the sphere wall object. The spherical wall is approximated by
 * a single plane applied at the point on the sphere which is closest to the center
 * of the cell. This works well for particle arrangements that are packed against
 * the wall, but loses accuracy for sparse particle distributions.
 * \param[in,out] c the Voronoi cell to be cut.
 * \param[in] (x,y,z) the location of the Voronoi cell.
 * \return True if the cell still exists, false if the cell is deleted. */
bool wall_sphere::cut_cell_base(voronoicell_neighbor &c,double x,double y,double z) {
	double xd=x-xc,yd=y-yc,zd=z-zc,dq;
	dq=xd*xd+yd*yd+zd*zd;
	if (dq>1e-5) {
		dq=2*(sqrt(dq)*rc-dq);
		return c.nplane(xd,yd,zd,dq,w_id);
	}
	return true;
}

/** Tests to see whether a point is inside the plane wall object.
 * \param[in] (x,y,z) the vector to test.
 * \return True if the point is inside, false if the point is outside. */
bool wall_plane::point_inside(double x,double y,double z) {
	return x*xc+y*yc+z*zc<ac;
}

/** Cuts a cell by the plane wall object.
 * \param[in,out] c the Voronoi cell to be cut.
 * \param[in] (x,y,z) the location of the Voronoi cell.
 * \return True if the cell still exists, false if the cell is deleted. */
bool wall_plane::cut_cell_base(voronoicell_neighbor &c,double x,double y,double z) {
	double dq=2*(ac-x*xc-y*yc-z*zc);
	return c.nplane(xc,yc,zc,dq,w_id);
}

/** Tests to see whether a point is inside the cylindrical wall object.
 * \param[in] (x,y,z) the vector to test.
 * \return True if the point is inside, false if the point is outside. */
bool wall_cylinder::point_inside(double x,double y,double z) {
	double xd=x-xc,yd=y-yc,zd=z-zc;
	double pa=(xd*xa+yd*ya+zd*za)*asi;
	xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
	return xd*xd+yd*yd+zd*zd<rc*rc;
}

/** Cuts a cell by the cylindrical wall object. The cylindrical wall is
 * approximated by a single plane applied at the point on the cylinder which is
 * closest to the center of the cell. This works well for particle arrangements
 * that are packed against the wall, but loses accuracy for sparse particle
 * distributions.
 * \param[in,out] c the Voronoi cell to be cut.
 * \param[in] (x,y,z) the location of the Voronoi cell.
 * \return True if the cell still exists, false if the cell is deleted. */
bool wall_cylinder::cut_cell_base(voronoicell_neighbor &c,double x,double y,double z) {
	double xd=x-xc,yd=y-yc,zd=z-zc;
	double pa=(xd*xa+yd*ya+zd*za)*asi;
	xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
	pa=xd*xd+yd*yd+zd*zd;
	if(pa>1e-5) {
		pa=2*(sqrt(pa)*rc-pa);
		return c.nplane(xd,yd,zd,pa,w_id);
	}
	return true;
}

/** Tests to see whether a point is inside the cone wall object.
 * \param[in] (x,y,z) the vector to test.
 * \return True if the point is inside, false if the point is outside. */
bool wall_cone::point_inside(double x,double y,double z) {
	double xd=x-xc,yd=y-yc,zd=z-zc;
	double pa=(xd*xa+yd*ya+zd*za)*asi;
	xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
	pa*=gra;
	if (pa<0) return false;
	pa*=pa;
	return xd*xd+yd*yd+zd*zd<pa;
}

/** Cuts a cell by the cone wall object. The conical wall is
 * approximated by a single plane applied at the point on the cone which is
 * closest to the center of the cell. This works well for particle arrangements
 * that are packed against the wall, but loses accuracy for sparse particle
 * distributions.
 * \param[in,out] c the Voronoi cell to be cut.
 * \param[in] (x,y,z) the location of the Voronoi cell.
 * \return True if the cell still exists, false if the cell is deleted. */
bool wall_cone::cut_cell_base(voronoicell_neighbor &c,double x,double y,double z) {
	double xd=x-xc,yd=y-yc,zd=z-zc;
	double xf,yf,zf,imoda;
	double pa=(xd*xa+yd*ya+zd*za)*asi;
	xd-=xa*pa;yd-=ya*pa;zd-=za*pa;
	pa=xd*xd+yd*yd+zd*zd;
	if(pa>1e-5) {
		pa=1/sqrt(pa);
		imoda=sqrt(asi);
		xf=-sang*imoda*xa+cang*pa*xd;
		yf=-sang*imoda*ya+cang*pa*yd;
		zf=-sang*imoda*za+cang*pa*zd;
		pa=2*(xf*(xc-x)+yf*(yc-y)+zf*(zc-z));
		return c.nplane(xf,yf,zf,pa,w_id);
	}
	return true;
}


