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


// This must always be the first include directive
#include "config.h"

#include "posdomain.h"


const char *GetRevisionInfo_posdomain(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_posdomain() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



bool CPosDomainEngine::Create(const CxDVec3Array &atoms, const CxDMatrix3 &cell, double mindist) {

	return Create(atoms,atoms.GetSize(),cell,mindist);
}


bool CPosDomainEngine::Create(const CxDVec3Array &atoms, int count, const CxDMatrix3 &cell, double mindist) {

	CxDVector3 tv;
	CxDMatrix3 invcell;
	int z, dx, dy, dz, d, ti;
	double tf;


	//mprintf("\n>>> Init >>>\n");

	invcell = cell;
	if (!invcell.Invert()) {
		eprintf("CPosDomainEngine::Create(): Error: Encountered singular cell matrix (cell volume is zero).\n");
		abort();
	}

	tf = sqrt( cell(0,0)*cell(0,0) + cell(0,1)*cell(0,1) + cell(0,2)*cell(0,2) );
	z = (int)floor(tf/mindist);
	if (z < 1)
		z = 1;
	if (z > 64)
		z = 64;
	//m_iRes[0] = AskUnsignedInteger("    Domain X grid: [%d] ",z,z);
	m_iRes[0] = z;

	tf = sqrt( cell(1,0)*cell(1,0) + cell(1,1)*cell(1,1) + cell(1,2)*cell(1,2) );
	z = (int)floor(tf/mindist);
	if (z < 1)
		z = 1;
	if (z > 64)
		z = 64;
	//m_iRes[1] = AskUnsignedInteger("    Domain Y grid: [%d] ",z,z);
	m_iRes[1] = z;

	tf = sqrt( cell(2,0)*cell(2,0) + cell(2,1)*cell(2,1) + cell(2,2)*cell(2,2) );
	z = (int)floor(tf/mindist);
	if (z < 1)
		z = 1;
	if (z > 64)
		z = 64;
	//m_iRes[2] = AskUnsignedInteger("    Domain Z grid: [%d] ",z,z);
	m_iRes[2] = z;

	if ((m_iRes[0] < 3) || (m_iRes[1] < 3) || (m_iRes[2] < 3)) {

		m_bTrivial = true;

		m_iaTrivialList.resize(atoms.GetSize());
		for (z=0;z<atoms.GetSize();z++)
			m_iaTrivialList[z] = z;

		return true;
	}

	if (m_iRes[0]*m_iRes[1]*m_iRes[2] > (int)m_oaDomains.size()) {
		ti = (int)m_oaDomains.size();
		m_oaDomains.resize(m_iRes[0]*m_iRes[1]*m_iRes[2]);
		for (z=ti;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
			m_oaDomains[z] = new CPosDomain();
	}

	for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
		m_oaDomains[z]->m_iaAtoms.clear();

	m_iaAtomDomain.resize(count);

	for (z=0;z<count;z++) {
		tv = invcell * atoms[z];
		tv[0] = fmod(tv[0],1.0);
		if (tv[0] < 0)
			tv[0] += 1.0;
		tv[1] = fmod(tv[1],1.0);
		if (tv[1] < 0)
			tv[1] += 1.0;
		tv[2] = fmod(tv[2],1.0);
		if (tv[2] < 0)
			tv[2] += 1.0;

		dx = (int)(tv[0] * m_iRes[0]);
		dy = (int)(tv[1] * m_iRes[1]);
		dz = (int)(tv[2] * m_iRes[2]);
		if (dx >= m_iRes[0])
			dx = m_iRes[0]-1;
		if (dy >= m_iRes[1])
			dy = m_iRes[1]-1;
		if (dz >= m_iRes[2])
			dz = m_iRes[2]-1;

		d = dx+dy*m_iRes[0]+dz*m_iRes[0]*m_iRes[1];

		m_oaDomains[d]->m_iaAtoms.push_back(z);
		m_iaAtomDomain[z] = d;
	}

	//mprintf("\n>>> Init Done >>>\n");

	return true;
}


bool CPosDomainEngine::CreateTrivial(const CxDVec3Array &atoms) {

	int z;


	m_bTrivial = true;

	m_iaTrivialList.resize(atoms.GetSize());
	for (z=0;z<atoms.GetSize();z++)
		m_iaTrivialList[z] = z;

	return true;
}


void CPosDomainEngine::FindNeighbors(int atom, std::vector<int> &nb) const {

	int d, d2, ox, oy, oz, rx, ry, rz, nx, ny, nz;


	//mprintf("\n*** Neighbors of %d ***\n",atom);

	if (m_bTrivial) {
		nb.insert( nb.end(), m_iaTrivialList.begin(), m_iaTrivialList.end() );
		return;
	}

	// Guaranteed resolution >= 3 in all directions
	
	d = m_iaAtomDomain[atom];
	ox = d % m_iRes[0];
	oy = (d / m_iRes[0]) % m_iRes[1];
	oz = (d / m_iRes[0] / m_iRes[1]);

	//mprintf("  is in Domain %d, which is %d %d %d.\n",d,ox,oy,oz);

	nb.clear();

	for (rz=oz-1;rz<=oz+1;rz++) {

		if (rz < 0)
			nz = rz+m_iRes[2];
		else if (rz >= m_iRes[2])
			nz = rz-m_iRes[2];
		else
			nz = rz;

		for (ry=oy-1;ry<=oy+1;ry++) {

			//mprintf("@ ry=%d (%d..%d)\n",ry,oy-1,oy+1);

			if (ry < 0)
				ny = ry+m_iRes[1];
			else if (ry >= m_iRes[1])
				ny = ry-m_iRes[1];
			else
				ny = ry;

			for (rx=ox-1;rx<=ox+1;rx++) {

				if (rx < 0)
					nx = rx+m_iRes[0];
				else if (rx >= m_iRes[0])
					nx = rx-m_iRes[0];
				else
					nx = rx;

				d2 = nx+ny*m_iRes[0]+nz*m_iRes[0]*m_iRes[1];

				//mprintf("    %2d %2d %2d / %4d: %d atoms.\n",nx,ny,nz,d2,(int)m_oaDomains[d2]->m_iaAtoms.size());

				nb.insert( nb.end(), m_oaDomains[d2]->m_iaAtoms.begin(), m_oaDomains[d2]->m_iaAtoms.end() );
			}
		}
	}
}


void CPosDomainEngine::GetDomainAndNeighbors(int domain, std::vector<int> &dom, std::vector<int> &nbh) const {

	int d2, ox, oy, oz, rx, ry, rz, nx, ny, nz;


	if (m_bTrivial) {
		dom.assign( m_iaTrivialList.begin(), m_iaTrivialList.end() );
		nbh.assign( m_iaTrivialList.begin(), m_iaTrivialList.end() );
		return;
	}

	dom.assign( m_oaDomains[domain]->m_iaAtoms.begin(), m_oaDomains[domain]->m_iaAtoms.end() );

	ox = domain % m_iRes[0];
	oy = (domain / m_iRes[0]) % m_iRes[1];
	oz = (domain / m_iRes[0] / m_iRes[1]);

	nbh.clear();

	for (rz=oz-1;rz<=oz+1;rz++) {

		if (rz < 0)
			nz = rz+m_iRes[2];
		else if (rz >= m_iRes[2])
			nz = rz-m_iRes[2];
		else
			nz = rz;

		for (ry=oy-1;ry<=oy+1;ry++) {

			if (ry < 0)
				ny = ry+m_iRes[1];
			else if (ry >= m_iRes[1])
				ny = ry-m_iRes[1];
			else
				ny = ry;

			for (rx=ox-1;rx<=ox+1;rx++) {

				if (rx < 0)
					nx = rx+m_iRes[0];
				else if (rx >= m_iRes[0])
					nx = rx-m_iRes[0];
				else
					nx = rx;

				d2 = nx+ny*m_iRes[0]+nz*m_iRes[0]*m_iRes[1];

				nbh.insert( nbh.end(), m_oaDomains[d2]->m_iaAtoms.begin(), m_oaDomains[d2]->m_iaAtoms.end() );
			}
		}
	}
}




