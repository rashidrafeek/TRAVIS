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


#ifndef DACF_H
#define DACF_H


// This must always be the first include directive
#include "config.h"

#include "xobject.h"
#include "xobarray.h"
#include "xdoublearray.h"
#include "nbsearch.h"
#include "moltools.h"
#include "nbexchange.h"


class CDACFSub : public CxObject
{
public:
	void ProcessPair(const CxIntArray &intervals);
	void Create(CConditionGroup *c);
	void Parse();
	CDACFSub();
	~CDACFSub();

	bool m_bNewMode;
	bool m_bBorderMode;
	bool m_bDistTrace;
	bool m_bIntermittend;
	double m_fIntGap;
//	bool m_bIntTravisStyle;
	bool m_bCorrectEq;
	bool m_bDebugOut;

	bool m_bIntNew;

	double m_fEqCounter;

	int m_iRefMol;
	int m_iShowMol;

	char *m_sName;
	void BuildName(const char *n);

	CConditionGroup *m_pCondition;
	CxObArray *m_oaAggregates;
	CxObArray *m_oaNbExPairs;

	C2DF *m_pDLDisp;
	CDF *m_pDACF;
	CDF *m_pDLDF;
	CDF *m_pNDF;
	CDF *m_pDDisp;
	CAF *m_pPairMSD;

	CDF **m_pNbExDF;
	CConditionGroup **m_pNbExConditions;

	CxIntArray *m_piaIntervals;
};


class CDACF : public CxObject
{
public:
	void CreateSubDACFStack(char *s);
	void CreateGridFit2DF(C2DF *df, int degree, bool intermittend);
//	void CreateGridFitDF(CDF *df, int degree, bool intermittend);
	void CalcGridFitParms();
	void UpdateNbEx(int rm, CDACFSub *dacfsub);
	void UpdateDACFSub(int rm, CTimeStep *t, CDACFSub *dacfsub);
	void FinishDACFSub(CTimeStep *t, CDACFSub *dacfsub);
	bool m_bRemoveMaxVel;
	double m_fMaxVel;
	CDACF();
	~CDACF();
	void Parse();
	CxDVector3 CalcCenter(CxDVec3Array *v, int i1, int i2);

	CxObArray m_oaSubDACFs;
	CConditionGroup *m_pCondition;
	int m_iFirstMol;
	int m_iSecondMol;
	char *m_sName;

	CAtomGroup *m_pCenterAtoms1;
	CAtomGroup *m_pCenterAtoms2;
	CxDoubleArray m_faWeight1;
	CxDoubleArray m_faWeight2;
	double m_fWeightMol1;

	int m_iLifetimeRes;
	int m_iDACFRes;
	int m_iDisplacementRes;
	double m_fLargestDisplacement;
	double m_fLargestLifetime;
	
	bool m_bDACFGrid;
	int m_iGridMode;
	bool m_bGridCon;
	bool m_bGridInt;
	bool m_bGridIntTravisStyle;
	double m_fGridIntGap;

	CNbExchange *m_pNbExchange;
	bool m_bFitDACF;
	int m_iFitDegreeMin;
	int m_iFitDegreeMax;

	double *m_pFitRMin;
	double *m_pFitRAvg;
	double *m_pFitRMax;

	int m_iGridX;
	int m_iGridY;

	double m_fGridXMin;
	double m_fGridXMax;
	double m_fGridYMin;
	double m_fGridYMax;

	bool m_bLifetimeSpectrum;
	CxObArray m_oaLTSpectra;
};

#endif
