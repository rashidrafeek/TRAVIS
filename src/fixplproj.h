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


#ifndef FIXPLPROJ_H
#define FIXPLPROJ_H


// This must always be the first include directive
#include "config.h"

#include "tools.h"
#include "timestep.h"
#include "2df.h"
#include <vector>



class CFixedPlProjObservation {
public:


	CFixedPlProjObservation() : m_p2DF(NULL), m_pPointPlaneDF(NULL) {
	}


	~CFixedPlProjObservation() {

		if (m_p2DF != NULL) {
			delete m_p2DF;
			m_p2DF = NULL;
		}
		if (m_pPointPlaneDF != NULL) {
			delete m_pPointPlaneDF;
			m_pPointPlaneDF = NULL;
		}
	}


	bool Parse(int i);

	bool m_bDifference;

	double m_fLayerWidth;

	std::vector<int> m_iaAtomList;

	C2DF *m_p2DF;

	CDF *m_pPointPlaneDF;

	CxString m_sName;

	CxDVector3 m_vPlaneBase;
	CxDVector3 m_vPlaneVec1;
	CxDVector3 m_vPlaneVec2;
	CxDVector3 m_vPlaneNormal;
};


class CFixedPlProj {
public:

	~CFixedPlProj() {
		int z;
		for (z=0;z<(int)m_oaObservations.size();z++)
			delete m_oaObservations[z];
		m_oaObservations.clear();
	}

	bool Parse();

	void Process(CTimeStep *ts);

	void Finish();

	std::vector<CFixedPlProjObservation*> m_oaObservations;
};



#endif


