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


#ifndef POSDOMAIN_H
#define POSDOMAIN_H


// This must always be the first include directive
#include "config.h"

#include "tools.h"
#include "xdvec3array.h"
#include "xdmatrix3.h"
#include <vector>


class CPosDomain {
public:

	std::vector<int> m_iaAtoms;
};


class CPosDomainEngine {
public:


	CPosDomainEngine() {

		m_iRes[0] = 0;
		m_iRes[1] = 0;
		m_iRes[2] = 0;
		m_bTrivial = false;
	}


	~CPosDomainEngine() {

		int z;

		for (z=0;z<(int)m_oaDomains.size();z++)
			delete m_oaDomains[z];
		m_oaDomains.clear();
	}


	bool Create(const CxDVec3Array &atoms, const CxDMatrix3 &cell, double mindist);


	bool Create(const CxDVec3Array &atoms, int count, const CxDMatrix3 &cell, double mindist);


	bool CreateTrivial(const CxDVec3Array &atoms);


	void FindNeighbors(int atom, std::vector<int> &nb) const;


	void GetDomainAndNeighbors(int domain, std::vector<int> &dom, std::vector<int> &nbh) const;


	std::vector<CPosDomain*> m_oaDomains;
	std::vector<int> m_iaAtomDomain;
	std::vector<int> m_iaTrivialList;

	bool m_bTrivial;

	int m_iRes[3];
};


#endif


