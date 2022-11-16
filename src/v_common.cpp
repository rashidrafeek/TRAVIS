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

/** \file common.cc
 * \brief Implementations of the small helper functions. */


 // This must always be the first include directive
#include "config.h"

#include "v_common.h"


using namespace std;


const char *GetRevisionInfo_v_common(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_v_common() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


/** \brief Prints a vector of integers.
 *
 * Prints a vector of integers.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
void voro_print_vector(vector<int> &v,FILE *fp) {
	int k=0,s=(int)v.size();
	while(k+4<s) {
		fprintf(fp,"%d %d %d %d ",v[k],v[k+1],v[k+2],v[k+3]);
		k+=4;
	}
	if(k+3<=s) {
		if(k+4==s) fprintf(fp,"%d %d %d %d",v[k],v[k+1],v[k+2],v[k+3]);
		else fprintf(fp,"%d %d %d",v[k],v[k+1],v[k+2]);
	} else {
		if(k+2==s) fprintf(fp,"%d %d",v[k],v[k+1]);
		else fprintf(fp,"%d",v[k]);
	}
}

/** \brief Prints a vector of doubles.
 *
 * Prints a vector of doubles.
 * \param[in] v the vector to print.
 * \param[in] fp the file stream to print to. */
void voro_print_vector(vector<double> &v,FILE *fp) {
	int k=0,s=(int)v.size();
	while(k+4<s) {
		fprintf(fp,"%g %g %g %g ",v[k],v[k+1],v[k+2],v[k+3]);
		k+=4;
	}
	if(k+3<=s) {
		if(k+4==s) fprintf(fp,"%g %g %g %g",v[k],v[k+1],v[k+2],v[k+3]);
		else fprintf(fp,"%g %g %g",v[k],v[k+1],v[k+2]);
	} else {
		if(k+2==s) fprintf(fp,"%g %g",v[k],v[k+1]);
		else fprintf(fp,"%g",v[k]);
	}
}

/** \brief Prints a vector a face vertex information.
 *
 * Prints a vector of face vertex information. A value is read, which
 * corresponds to the number of vertices in the next face. The routine reads
 * this number of values and prints them as a bracked list. This is repeated
 * until the end of the vector is reached.
 * \param[in] v the vector to interpret and print.
 * \param[in] fp the file stream to print to. */
void voro_print_face_vertices(vector<int> &v,FILE *fp) {
	int j,k=0,l;
	if(v.size()>0) {
		l=v[k++];
		if(l<=1) {
			if(l==1) fprintf(fp,"(%d)",v[k++]);
			else fprintf( fp, "()");
		} else {
			j=k+l;
			fprintf(fp,"(%d",v[k++]);
			while(k<j) fprintf(fp,",%d",v[k++]);
			fprintf( fp, ")");
		}
		while((unsigned int) k<v.size()) {
			l=v[k++];
			if(l<=1) {
				if(l==1) fprintf(fp," (%d)",v[k++]);
				else fprintf( fp," ()");
			} else {
				j=k+l;
				fprintf(fp," (%d",v[k++]);
				while(k<j) fprintf(fp,",%d",v[k++]);
				fprintf( fp,")");
			}
		}
	}
}

