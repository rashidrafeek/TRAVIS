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


#ifndef CHARGEVAR2_H
#define CHARGEVAR2_H


// This must always be the first include directive
#include "config.h"


#ifdef NEW_CHARGEVAR


#include <array>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>

#include "xobject.h"
#include "xobarray.h"
#include "xdoublearray.h"
#include "xdvectorn.h"
#include "xintarray.h"
#include "xptrarray.h"
#include "bqb.h"

#include "libtravis/calcvoronoi.h"
#include "libtravis/calcvoronoiintegrals.h"
#include "libtravis/cellinfo.h"
#include "libtravis/constants.h"
#include "libtravis/filereadercube.h"
#include "libtravis/moleculeset.h"
#include "libtravis/recognizemolecules.h"
#include "libtravis/snapshot.h"
#include "libtravis/trajectory.h"
#include "libtravis/voronoitessellation.h"



class CChargeVar2_Snapshot {
public:

	CChargeVar2_Snapshot() : m_pCellInfo(NULL), m_pElectronDensity(NULL) {
	}

	~CChargeVar2_Snapshot() {
		if (m_pCellInfo != NULL) {
			delete m_pCellInfo;
			m_pCellInfo = NULL;
		}
		if (m_pElectronDensity != NULL) {
			delete m_pElectronDensity;
			m_pElectronDensity = NULL;
		}
	}

	std::vector<Libtravis::Travis::Vector3> m_vaAtomPositions;
	Libtravis::Travis::CellInfo *m_pCellInfo;
	Libtravis::Travis::DataGrid3D<double> *m_pElectronDensity;
};



class CChargeVar2 {

public:
	CChargeVar2();
	~CChargeVar2();

	void Parse();

	void PerformAnalysis();

	bool CalculateCharges();

	void CalculateMinMaxVar();

	double CalculateObjective(CxDVectorN *vec);

	bool LineSearch(CxDVectorN *pos, CxDVectorN *dir, double maxlen, double *result);

	double REC_GoldenSectionSearch(double a, double b, double c, double fa, double fb, double fc, CxDVectorN *pos, CxDVectorN *dir, int depth);

	void ComputeGradient(CxDVectorN *pos, CxDVectorN *grad);

	void WriteRadiiRestart();

	bool ReadRestartRadii( bool mol );

	int m_iDimensions;

	CxIntArray m_iaAtomVectorAssignment;

	CxDVectorN m_vPosition;
	CxDVectorN m_vPositionI1;
	CxDVectorN m_vPositionI2;

	CxPtrArray m_saComponentNames;
	CxIntArray m_iaComponentElements;

	CxObArray m_oaAtomCharges; // Contains CxDoubleArrays

	CxDoubleArray m_faMolChargeVar;
	CxDoubleArray m_faMolChargeMean;
	CxDoubleArray m_faMolChargeMin;
	CxDoubleArray m_faMolChargeMax;

	CxDoubleArray m_faAtomChargeMean;
	CxDoubleArray m_faAtomChargeMin;
	CxDoubleArray m_faAtomChargeMax;
	CxDoubleArray m_faAtomChargeVar;

	double m_fGradientEps;
	double m_fLSEps;
	double m_fVarEps;
	double m_fInitialStepLength;
	int m_iGradientSamples;
	bool m_bGradientOutlier;
//	bool m_bMolVar;
	int m_iVarMode; // 0 - Atomwise Charge Variance, 1 - Molecule-wise Charge Variance, 2 - mix both
	double m_fVarMixing; // Percentage of Molecule-wise Charge Variance in total objective
	bool m_bConjugateGradient;
	bool m_bNoOpt;
	bool m_bOutputEveryStep;

	bool m_bOptimize;
	bool m_bTwoStep;
	bool m_bOneStepTargetAtomic;
	bool m_bWrapCoords;
	bool m_bVerboseTiming;

	int m_iInitialRadii1;
	int m_iInitialRadii2;

	int m_iGradCounter;
	int m_iLSCounter;

	CxDoubleArray m_faLSX;
	CxDoubleArray m_faLSY;

	bool m_bElement;

	int m_iLoadFrames;
	int m_iStride;

	std::vector<double> m_faRadii1;
	std::vector<double> m_faRadii2;

	bool m_bSecPass;

	bool m_bPerformFirstPass;
	bool m_bPerformSecondPass;
	bool m_bUseThreads;
	int m_iThreads;

	std::vector<CChargeVar2_Snapshot*> m_oaSnapshots;
	std::vector<std::vector<std::size_t> > m_iaaMaskCellIds;

	std::array<unsigned int, 3> m_iaRefining; 

	bool m_bWeightGradient;

	CBQBFile *m_pBQBFile;

	std::mutex m_oLock;
	std::vector<std::thread*> m_oaThreads;

	void ThreadMain( int id, int nfirst, int nlast );

	bool loadData( const char *filename, int steps, int stride );

	std::vector<double> calcCharges( const std::vector<double>& radii1, const std::vector<double>& radii2, std::size_t frame, int threadid );


};



#endif


#endif



