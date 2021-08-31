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


#ifndef TDDF_H
#define TDDF_H


// This must always be the first include directive
#include "config.h"

#include <vector>
#include "timestep.h"
#include "2df.h"
#include "atomgroup.h"



class CTDDFObservationBase {
public:


	explicit CTDDFObservationBase(int idx) : m_pHistogram(NULL), m_iIndex(idx) {
	}


	virtual ~CTDDFObservationBase() {
		if (m_pHistogram != NULL) {
			delete m_pHistogram;
			m_pHistogram = NULL;
		}
	}


	virtual bool Parse() = 0;

	virtual bool Initialize() = 0;

	virtual void ProcessStep(const CTimeStep *ts) = 0;

	virtual void Finish() = 0;

protected:
	C2DF *m_pHistogram;
	int m_iIndex;
};



class CTDDFObservationRDF : public CTDDFObservationBase {
public:


	explicit CTDDFObservationRDF(int idx) : CTDDFObservationBase(idx) {
	}


	~CTDDFObservationRDF() {
	}


	bool Parse();

	bool Initialize();

	void ProcessStep(const CTimeStep *ts);

	void Finish();

private:
	bool m_bInter;
	int m_iRMIndex;
	int m_iOMIndex;
	CAtomGroup m_oRMAG;
	CAtomGroup m_oOMAG;

	double m_fMinDist;
	double m_fMaxDist;

	int m_iResX;
	int m_iResY;

	double m_fOverlap;

	bool m_bCorrectRadial;
	bool m_bUniform;
	bool m_bNS;

	CxString m_sName;

	std::vector<int> m_iaSliceCounter;
};



class CTDDFEngine {
public:

	bool Parse();

	bool Initialize();

	void ProcessStep(const CTimeStep *ts);

	void Finish();

private:
	std::vector<CTDDFObservationBase*> m_oaObservations;
};




#endif



