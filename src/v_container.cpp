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

/** \file container.cc
 * \brief Function implementations for the container and related classes. */


 // This must always be the first include directive
#include "config.h"

#include "v_container.h"
#include "tools.h"


const char *GetRevisionInfo_v_container(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_v_container() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


/** The class constructor sets up the geometry of container, initializing the
 * minimum and maximum coordinates in each direction, and setting whether each
 * direction is periodic or not. It divides the container into a rectangular
 * grid of blocks, and allocates memory for each of these for storing particle
 * positions and IDs.
 * \param[in] (ax_,bx_) the minimum and maximum x coordinates.
 * \param[in] (ay_,by_) the minimum and maximum y coordinates.
 * \param[in] (az_,bz_) the minimum and maximum z coordinates.
 * \param[in] (nx_,ny_,nz_) the number of grid blocks in each of the three
 *			    coordinate directions.
 * \param[in] (xperiodic_,yperiodic_,zperiodic_) flags setting whether the
 *                                               container is periodic in each
 *                                               coordinate direction.
 * \param[in] init_mem the initial memory allocation for each block.
 * \param[in] ps_ the number of floating point entries to store for each
 *                particle. */
container_base::container_base(double ax_,double bx_,double ay_,double by_,double az_,double bz_,
		int nx_,int ny_,int nz_,bool xperiodic_,bool yperiodic_,bool zperiodic_,int init_mem,int ps_)
	: voro_base(nx_,ny_,nz_,(bx_-ax_)/nx_,(by_-ay_)/ny_,(bz_-az_)/nz_),
	ax(ax_), bx(bx_), ay(ay_), by(by_), az(az_), bz(bz_),
	xperiodic(xperiodic_), yperiodic(yperiodic_), zperiodic(zperiodic_),
	id(new int*[nxyz]), p(new double*[nxyz]), co(new int[nxyz]), mem(new int[nxyz]), ps(ps_) {
	int l;
	for(l=0;l<nxyz;l++) co[l]=0;
	for(l=0;l<nxyz;l++) mem[l]=init_mem;
	for(l=0;l<nxyz;l++) id[l]=new int[init_mem];
	for(l=0;l<nxyz;l++) p[l]=new double[ps*init_mem];
}

/** The container destructor frees the dynamically allocated memory. */
container_base::~container_base() {
	int l;
	for(l=0;l<nxyz;l++) delete [] p[l];
	for(l=0;l<nxyz;l++) delete [] id[l];
	delete [] id;
	delete [] p;
	delete [] co;
	delete [] mem;
}



/** This routine takes a particle position vector, tries to remap it into the
 * primary domain. If successful, it computes the region into which it can be
 * stored and checks that there is enough memory within this region to store
 * it.
 * \param[out] ijk the region index.
 * \param[in,out] (x,y,z) the particle position, remapped into the primary
 *                        domain if necessary.
 * \return True if the particle can be successfully placed into the container,
 * false otherwise. */
inline bool container_base::put_locate_block(int &ijk,double &x,double &y,double &z) {
	if(put_remap(ijk,x,y,z)) {
		if(co[ijk]==mem[ijk]) add_particle_memory(ijk);
		return true;
	}
#if VOROPP_REPORT_OUT_OF_BOUNDS ==1
	fprintf(stderr,"Out of bounds: (x,y,z)=(%g,%g,%g)\n",x,y,z);
#endif
	return false;
}

/** Takes a particle position vector and computes the region index into which
 * it should be stored. If the container is periodic, then the routine also
 * maps the particle position to ensure it is in the primary domain. If the
 * container is not periodic, the routine bails out.
 * \param[out] ijk the region index.
 * \param[in,out] (x,y,z) the particle position, remapped into the primary
 *                        domain if necessary.
 * \return True if the particle can be successfully placed into the container,
 * false otherwise. */
inline bool container_base::put_remap(int &ijk,double &x,double &y,double &z) {
	int l;

	ijk=step_int((x-ax)*xsp);
	if(xperiodic) {l=step_mod(ijk,nx);x+=boxx*(l-ijk);ijk=l;}
	else if(ijk<0||ijk>=nx) return false;

	int j=step_int((y-ay)*ysp);
	if(yperiodic) {l=step_mod(j,ny);y+=boxy*(l-j);j=l;}
	else if(j<0||j>=ny) return false;

	int k=step_int((z-az)*zsp);
	if(zperiodic) {l=step_mod(k,nz);z+=boxz*(l-k);k=l;}
	else if(k<0||k>=nz) return false;

	ijk+=nx*j+nxy*k;
	return true;
}

/** Takes a position vector and attempts to remap it into the primary domain.
 * \param[out] (ai,aj,ak) the periodic image displacement that the vector is in,
 *                       with (0,0,0) corresponding to the primary domain.
 * \param[out] (ci,cj,ck) the index of the block that the position vector is
 *                        within, once it has been remapped.
 * \param[in,out] (x,y,z) the position vector to consider, which is remapped
 *                        into the primary domain during the routine.
 * \param[out] ijk the block index that the vector is within.
 * \return True if the particle is within the container or can be remapped into
 * it, false if it lies outside of the container bounds. */
inline bool container_base::remap(int &ai,int &aj,int &ak,int &ci,int &cj,int &ck,double &x,double &y,double &z,int &ijk) {
	ci=step_int((x-ax)*xsp);
	if(ci<0||ci>=nx) {
		if(xperiodic) {ai=step_div(ci,nx);x-=ai*(bx-ax);ci-=ai*nx;}
		else return false;
	} else ai=0;

	cj=step_int((y-ay)*ysp);
	if(cj<0||cj>=ny) {
		if(yperiodic) {aj=step_div(cj,ny);y-=aj*(by-ay);cj-=aj*ny;}
		else return false;
	} else aj=0;

	ck=step_int((z-az)*zsp);
	if(ck<0||ck>=nz) {
		if(zperiodic) {ak=step_div(ck,nz);z-=ak*(bz-az);ck-=ak*nz;}
		else return false;
	} else ak=0;

	ijk=ci+nx*cj+nxy*ck;
	return true;
}



/** Increase memory for a particular region.
 * \param[in] i the index of the region to reallocate. */
void container_base::add_particle_memory(int i) {
	int l,nmem=mem[i]<<1;

	// Carry out a check on the memory allocation size, and
	// print a status message if requested
	if(nmem>max_particle_memory)
		voro_fatal_error("Absolute maximum memory allocation exceeded",VOROPP_MEMORY_ERROR);
#if VOROPP_VERBOSE >=3
	fprintf(stderr,"Particle memory in region %d scaled up to %d\n",i,nmem);
#endif

	// Allocate new memory and copy in the contents of the old arrays
	int *idp=new int[nmem];
	for(l=0;l<co[i];l++) idp[l]=id[i][l];
	double *pp=new double[ps*nmem];
	for(l=0;l<ps*co[i];l++) pp[l]=p[i][l];

	// Update pointers and delete old arrays
	mem[i]=nmem;
	delete [] id[i];id[i]=idp;
	delete [] p[i];p[i]=pp;
}



/** Outputs the a list of all the container regions along with the number of
 * particles stored within each. */
void container_base::region_count() {
	int i,j,k,*cop=co;
	for(k=0;k<nz;k++) for(j=0;j<ny;j++) for(i=0;i<nx;i++)
		mprintf("Region (%d,%d,%d): %d particles\n",i,j,k,*(cop++));
}




/** This function tests to see if a given vector lies within the container
 * bounds and any walls.
 * \param[in] (x,y,z) the position vector to be tested.
 * \return True if the point is inside the container, false if the point is
 *         outside. */
bool container_base::point_inside(double x,double y,double z) {
	if(x<ax||x>bx||y<ay||y>by||z<az||z>bz) return false;
	return point_inside_walls(x,y,z);
}

/** Draws an outline of the domain in gnuplot format.
 * \param[in] fp the file handle to write to. */
void container_base::draw_domain_gnuplot(FILE *fp) {
	fprintf(fp,"%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n",ax,ay,az,bx,ay,az,bx,by,az,ax,by,az);
	fprintf(fp,"%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n",ax,by,bz,bx,by,bz,bx,ay,bz,ax,ay,bz);
	fprintf(fp,"%g %g %g\n\n%g %g %g\n%g %g %g\n\n",ax,by,bz,ax,ay,az,ax,ay,bz);
	fprintf(fp,"%g %g %g\n%g %g %g\n\n%g %g %g\n%g %g %g\n\n",bx,ay,az,bx,ay,bz,bx,by,az,bx,by,bz);
}

/** Draws an outline of the domain in POV-Ray format.
 * \param[in] fp the file handle to write to. */
void container_base::draw_domain_pov(FILE *fp) {
	fprintf(fp,"cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
		   "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",ax,ay,az,bx,ay,az,ax,by,az,bx,by,az);
	fprintf(fp,"cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
		   "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",ax,by,bz,bx,by,bz,ax,ay,bz,bx,ay,bz);
	fprintf(fp,"cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
		   "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",ax,ay,az,ax,by,az,bx,ay,az,bx,by,az);
	fprintf(fp,"cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
		   "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",bx,ay,bz,bx,by,bz,ax,ay,bz,ax,by,bz);
	fprintf(fp,"cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
		   "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",ax,ay,az,ax,ay,bz,bx,ay,az,bx,ay,bz);
	fprintf(fp,"cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n"
		   "cylinder{<%g,%g,%g>,<%g,%g,%g>,rr}\n",bx,by,az,bx,by,bz,ax,by,az,ax,by,bz);
	fprintf(fp,"sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n"
		   "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n",ax,ay,az,bx,ay,az,ax,by,az,bx,by,az);
	fprintf(fp,"sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n"
		   "sphere{<%g,%g,%g>,rr}\nsphere{<%g,%g,%g>,rr}\n",ax,ay,bz,bx,ay,bz,ax,by,bz,bx,by,bz);
}


/** The wall_list constructor sets up an array of pointers to wall classes. */
wall_list::wall_list() : walls(new wall*[init_wall_size]), wep(walls), wel(walls+init_wall_size),
	current_wall_size(init_wall_size) {}

/** The wall_list destructor frees the array of pointers to the wall classes.
 */
wall_list::~wall_list() {
	delete [] walls;
}

/** Adds all of the walls on another wall_list to this class.
 * \param[in] wl a reference to the wall class. */
void wall_list::add_wall(wall_list &wl) {
	for(wall **wp=wl.walls;wp<wl.wep;wp++) add_wall(*wp);
}

/** Deallocates all of the wall classes pointed to by the wall_list. */
void wall_list::deallocate() {
	for(wall **wp=walls;wp<wep;wp++) delete *wp;
}

/** Increases the memory allocation for the walls array. */
void wall_list::increase_wall_memory() {
	current_wall_size<<=1;
	if(current_wall_size>max_wall_size)
		voro_fatal_error("Wall memory allocation exceeded absolute maximum",VOROPP_MEMORY_ERROR);
	wall **nwalls=new wall*[current_wall_size],**nwp=nwalls,**wp=walls;
	while(wp<wep) *(nwp++)=*(wp++);
	delete [] walls;
	walls=nwalls;wel=walls+current_wall_size;wep=nwp;
}

