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


#ifndef ORDER_H
#define ORDER_H


// This must always be the first include directive
#include "config.h"

#include "timestep.h"
#include <vector>



class COrderEngine;


class COrderAnalysis {
public:

	COrderAnalysis(COrderEngine *parent) : m_pParent(parent) {
	}

	virtual bool Parse(CTimeStep *ts) = 0;

	virtual void ProcessStep(CTimeStep *ts) = 0;

	virtual void Finish() = 0;

	int m_iNumber;

	COrderEngine *m_pParent;

private:
	COrderAnalysis(); // "delete" the default constructor by making it private
};


class COrderChain {
public:

	int m_iMolecule;
	int m_iCarbonIndex;

	std::vector<int> m_iaAtoms;
};


class COrderEngine {
public:

	bool Parse(CTimeStep *ts);

	void ProcessStep(CTimeStep *ts);

	void Finish();

	void DetectChains();

	void REC_DetectChains(int depth, int pos, std::vector<int> &touched);

	std::vector<COrderAnalysis*> m_oaAnalyses;

	std::vector<std::vector<COrderChain*> > m_oaChains;

	int m_iChainsTotal;

	CSingleMolecule *m_pTempSM;
	std::vector<std::string> m_saChainAllowedSub;
};



#endif


