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


#ifndef AGGRTOPO_H
#define AGGRTOPO_H


// This must always be the first include directive
#include "config.h"

#include <vector>
#include <string>



class CTimeStep;
class CAtomGroup;
class CDF;



class CAggrTopoGroup {
public:

	bool Parse();

	int m_iObservation;
	std::vector<CAtomGroup*> m_oaAtoms;
	std::vector<int> m_iaAtoms;
	std::string m_sLongName;
	std::string m_sShortName;
	std::vector<CDF*> m_oaRDFs;
	std::vector<CAggrTopoGroup*> m_oaObsGroups;
	std::vector<double> m_faAggrValues;
	bool m_bDonor;
};



class CAggrTopoEngine {
public:


	CAggrTopoEngine() : m_pNeighborMatrix(NULL) {
	}


	~CAggrTopoEngine() {
		if (m_pNeighborMatrix != NULL) {
			delete[] m_pNeighborMatrix;
			m_pNeighborMatrix = NULL;
		}
	}


	bool Parse();

	void Initialize();

	void BuildNeighborMatrix();

	void ProcessStep(CTimeStep *ts);

	void Finish();

	std::vector<CAggrTopoGroup*> m_oaGroups;
	std::vector<CAggrTopoGroup*> m_oaDonorGroups;
	std::vector<CAggrTopoGroup*> m_oaAcceptorGroups;

	std::vector<int> m_iaIntra;

	//std::vector<CDF*> m_oaRDFs;

	double m_fMaxDist;

	int m_iResolution;

	unsigned char *m_pNeighborMatrix;

	bool m_bMinMax;

	double m_fCutoffDist;

	//CPosDomainEngine *m_pDomainEngine;
};



#endif





