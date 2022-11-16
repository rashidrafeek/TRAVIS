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

/** \file pre_container.cc
 * \brief Function implementations for the pre_container and related classes.
 */


 // This must always be the first include directive
#include "config.h"

#include <stdio.h>
#include <math.h>

#include "v_config.h"
#include "v_pre_container.h"
#include "v_container_prd.h"


const char *GetRevisionInfo_v_pre_container(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_v_pre_container() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


class container_periodic_poly;


/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction. It allocates an initial
 * chunk into which to store particle information.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
 * \param[in] (xperiodic_,yperiodic_,zperiodic_ ) flags setting whether the
 *                                                container is periodic in each
 *                                                coordinate direction.
 * \param[in] ps_ the number of floating point entries to store for each
 *                particle. */
pre_container_base::pre_container_base(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
	bool xperiodic_,bool yperiodic_,bool zperiodic_,int ps_) :
	ax(ax_), bx(bx_), ay(ay_), by(by_), az(az_), bz(bz_),
	xperiodic(xperiodic_), yperiodic(yperiodic_), zperiodic(zperiodic_), ps(ps_),
	index_sz(init_chunk_size), pre_id(new int*[index_sz]), end_id(pre_id),
	pre_p(new double*[index_sz]), end_p(pre_p) {
		ch_id=*end_id=new int[pre_container_chunk_size];
		l_id=end_id+index_sz;e_id=ch_id+pre_container_chunk_size;
		ch_p=*end_p=new double[ps*pre_container_chunk_size];
}

/** The destructor frees the dynamically allocated memory. */
pre_container_base::~pre_container_base() {
	delete [] *end_p;
	delete [] *end_id;
	while (end_id!=pre_id) {
		end_p--;
		delete [] *end_p;
		end_id--;
		delete [] *end_id;
	}
	delete [] pre_p;
	delete [] pre_id;
}

/** Makes a guess at the optimal grid of blocks to use, computing in
 * a way that
 * \param[out] (nx,ny,nz) the number of blocks to use. */
void pre_container_base::guess_optimal(int &nx,int &ny,int &nz) {
	double dx = (bx-ax),dy = (by-ay),dz = (bz-az);
	double ilscale = (mypow(total_particles()/(optimal_particles*dx*dy*dz),1/3.0));
	nx=int(dx*ilscale+1);
	ny=int(dy*ilscale+1);
	nz=int(dz*ilscale+1);
}


/** Stores a particle ID and position, allocating a new memory chunk if necessary.
 * \param[in] n the numerical ID of the inserted particle.
 * \param[in] (x,y,z) the position vector of the inserted particle.
 * \param[in] r the radius of the particle. */
void pre_container_poly::put(int n,double x,double y,double z,double r) {
	if((xperiodic||(x>=ax&&x<=bx))&&(yperiodic||(y>=ay&&y<=by))&&(zperiodic||(z>=az&&z<=bz))) {
		if(ch_id==e_id) new_chunk();
		*(ch_id++)=n;
		*(ch_p++)=x;*(ch_p++)=y;*(ch_p++)=z;*(ch_p++)=r;
	}
#if VOROPP_REPORT_OUT_OF_BOUNDS ==1
	else fprintf(stderr,"Out of bounds: (x,y,z)=(%g,%g,%g)\n",x,y,z);
#endif
}


/** Transfers the particles stored within the class to a container_poly class.
 * \param[in] con the container_poly class to transfer to. */
void pre_container_poly::setup(container_periodic_poly &con) {
	int **c_id=pre_id,*idp,*ide,n;
	double **c_p=pre_p,*pp,x,y,z,r;
	while(c_id<end_id) {
		idp=*(c_id++);ide=idp+pre_container_chunk_size;
		pp=*(c_p++);
		while(idp<ide) {
			n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);r=*(pp++);
			con.put(n,x,y,z,r);
		}
	}
	idp=*c_id;
	pp=*c_p;
	while(idp<ch_id) {
		n=*(idp++);x=*(pp++);y=*(pp++);z=*(pp++);r=*(pp++);
		con.put(n,x,y,z,r);
	}
}


/** Import a list of particles from an open file stream, also storing the order
 * of that the particles are read. Entries of four numbers (Particle ID, x
 * position, y position, z position) are searched for. If the file cannot be
 * successfully read, then the routine causes a fatal error.
 * \param[in] fp the file handle to read from. */
void pre_container_poly::import(FILE *fp) {
	int i,j;
	double x,y,z,r;
	while((j=fscanf(fp,"%d %lg %lg %lg %lg",&i,&x,&y,&z,&r))==5) put(i,x,y,z,r);
	if(j!=EOF) voro_fatal_error("File import error",VOROPP_FILE_ERROR);
}

/** Allocates a new chunk of memory for storing particles. */
void pre_container_base::new_chunk() {
	end_id++;end_p++;
	if(end_id==l_id) extend_chunk_index();
	ch_id=*end_id=new int[pre_container_chunk_size];
	e_id=ch_id+pre_container_chunk_size;
	ch_p=*end_p=new double[ps*pre_container_chunk_size];
}

/** Extends the index of chunks. */
void pre_container_base::extend_chunk_index() {
	index_sz<<=1;
	if(index_sz>max_chunk_size)
		voro_fatal_error("Absolute memory limit on chunk index reached",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=2
	fprintf(stderr,"Pre-container chunk index scaled up to %d\n",index_sz);
#endif
	int **n_id=new int*[index_sz],**p_id=n_id,**c_id=pre_id;
	double **n_p=new double*[index_sz],**p_p=n_p,**c_p=pre_p;
	while(c_id<end_id) {
		*(p_id++)=*(c_id++);
		*(p_p++)=*(c_p++);
	}
	delete [] pre_id;pre_id=n_id;end_id=p_id;l_id=pre_id+index_sz;
	delete [] pre_p;pre_p=n_p;end_p=p_p;
}

