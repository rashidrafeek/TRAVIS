/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

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


#ifndef REACTIVE_H
#define REACTIVE_H


// This must always be the first include directive
#include "config.h"

#include <vector>
#include "xdvector3.h"
#include "xdvec3array.h"
#include "xdmatrix3.h"


class CReactiveRingBuffer {
public:

	CReactiveRingBuffer();

	~CReactiveRingBuffer();

	void Init(int depth, int atomcount);

	void PushFrame(const CxDVec3Array &va);

	void PushCellMatrix(const CxDMatrix3 &m);

	void PushCellVector(double x, double y, double z);


	const CxDVector3 & Get(int atom, int depth) const {
		int t = m_iCurrentPos - depth;
		if (t < 0)
			t += m_iDepth;
		return m_vaaCoords[atom][t];
	}


	int m_iCurrentPos;
	int m_iDepth;
	int m_iAtomCount;

	std::vector<std::vector<CxDVector3> > m_vaaCoords;

	std::vector<CxDMatrix3> m_maCellMatrix;

	std::vector<CxDVector3> m_vaCellVector;
};


class CReactiveEngine {
public:

	CReactiveEngine();

	~CReactiveEngine();

	void Parse();

	bool m_bFirst;

	int m_iBondPersist;
};


#endif


