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


#ifndef CONTACTMATRIX_H
#define CONTACTMATRIX_H


// This must always be the first include directive
#include "config.h"

#include <vector>
#include "xstring.h"
#include "df.h"
#include "posdomain.h"


class CTimeStep;
class CAtomGroup;


class CContactMatrixObservation {
public:

	CContactMatrixObservation();

	~CContactMatrixObservation();

	bool Parse();

	void Initialize();

	void Finish();

	void BuildNeighborMatrices();

	bool m_bInter;
	bool m_bIntra;

	bool m_bJoinEquivalent;

	bool m_bWriteRDFs;

	bool m_bExclude12;
	bool m_bExclude13;
	bool m_bExclude14;

	double m_fMaxDist;

	int m_iResolution;

	double m_fHeightCut;
	double m_fDistRange;

	int m_iRows;
	int m_iCols;

	int m_iObservation;

	std::vector<CAtomGroup*> m_oaRows;
	std::vector<CAtomGroup*> m_oaColumns;

	std::vector<int> m_iaRowIndex;
	std::vector<int> m_iaColumnIndex;

	std::vector<int> m_iaRowMolecule;
	std::vector<int> m_iaColumnMolecule;

	std::vector<int> m_iaRowAtomCount;
	std::vector<int> m_iaColumnAtomCount;

	std::vector<int> m_iaRowRealElement;
	std::vector<int> m_iaColumnRealElement;

	std::vector<CxString> m_saRowLabels;
	std::vector<CxString> m_saColumnLabels;

	CDF *m_oaMatrix;

	std::vector<int> m_iaFirstMaximum;
	std::vector<int> m_iaFirstMinimum;

	std::vector<std::vector<int> > m_iaaNeighborMatrices;
	std::vector<int> m_iaLocalAtomIndex;
};


class CContactMatrix {
public:

	CContactMatrix();

	~CContactMatrix();

	bool Parse();

	void Initialize();

	void ProcessStep(CTimeStep *ts);

	void Finish();

	std::vector<CContactMatrixObservation*> m_oaObservations;

	double m_fMaxDist;

	CPosDomainEngine *m_pDomainEngine;
};


#endif



