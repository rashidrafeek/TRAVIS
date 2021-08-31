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


#ifndef ROA_H
#define ROA_H


// This must always be the first include directive
#include "config.h"

#include <vector>
#include "tools.h"
#include "timestep.h"
#include "tensor.h"
#include "bqb.h"


class CCrossCorrelation;


#define ROA_FIELD_NONE  0
#define ROA_FIELD_EXP   1
#define ROA_FIELD_EXN   2
#define ROA_FIELD_EYP   3
#define ROA_FIELD_EYN   4
#define ROA_FIELD_EZP   5
#define ROA_FIELD_EZN   6

#define ROA_MODE_COMBINED  0
#define ROA_MODE_GATHER    1
#define ROA_MODE_ANALYZE   2

#define ROA_SPECTRUM_IR     1
#define ROA_SPECTRUM_RAMAN  2
#define ROA_SPECTRUM_VCD    3
#define ROA_SPECTRUM_SFG    4
#define ROA_SPECTRUM_ROA    5




class CSmoothener {
public:
	CSmoothener() : m_pForward(NULL), m_pInverse(NULL), m_iLength(-1), m_iInternalLength(-1), m_fWaveNumber(0) {
	}

	~CSmoothener() {
		if (m_pForward != NULL) {
			delete m_pForward;
			m_pForward = NULL;
		}
		if (m_pInverse != NULL) {
			delete m_pInverse;
			m_pInverse = NULL;
		}
	}

	void Init(int length, double wavenumber);

	void Smoothen(const std::vector<double> &inp, std::vector<double> &outp);
	void Smoothen(const std::vector<CxDVector3> &inp, std::vector<CxDVector3> &outp);
	void Smoothen(const std::vector<CxDMatrix3> &inp, std::vector<CxDMatrix3> &outp);

	void Smoothen(std::vector<double> &data) {
		Smoothen(data,data);
	}

	void Smoothen(std::vector<CxDVector3> &data) {
		Smoothen(data,data);
	}

	void Smoothen(std::vector<CxDMatrix3> &data) {
		Smoothen(data,data);
	}

	CFFT *m_pForward;
	CFFT *m_pInverse;
	int m_iLength;
	int m_iInternalLength;
	double m_fWaveNumber;
};






class CROAObservation {
public:
	CROAObservation();
	CROAObservation(const CROAObservation &t);
	~CROAObservation();


	CxString m_sName;
	CxString m_sMolName;
	CxString m_sTypeName;

	std::vector<std::vector<int> > m_iaMolSelection;

	std::vector<double> m_faACF;
	std::vector<double> m_faACF2;
	std::vector<double> m_faACF3;
	std::vector<double> m_faACF4;
	std::vector<double> m_faSpectrum;
	std::vector<double> m_faSpectrum2;
	std::vector<double> m_faSpectrum3;

	std::vector<double> m_faACF_CTP;
	std::vector<double> m_faACF2_CTP;
	std::vector<double> m_faACF3_CTP;
	std::vector<double> m_faACF4_CTP;
	std::vector<double> m_faACF_CTN;
	std::vector<double> m_faACF2_CTN;
	std::vector<double> m_faACF3_CTN;
	std::vector<double> m_faACF4_CTN;

	bool m_bSingleACFs;
	std::vector<std::vector<double> > m_faaSingleACFs;
	std::vector<std::vector<double> > m_faaSingleACFs_CTP;
	std::vector<std::vector<double> > m_faaSingleACFs_CTN;

	std::vector<double> m_faWFN;

	bool m_bMolPattern;


	int m_iType;
	int m_iDepth;
	int m_iWindowFunction;
	int m_iWindowFunctionParameter;
	int m_iZeroPadding;
	double m_fSpecResolution;
	int m_iSpecLength;
	bool m_bFiniteDifferenceCorrection;
	bool m_bSaveACF;
	int m_iQuantumCorrection;
	double m_fQuantumCorrectionTemperature;
	bool m_bCorrectFrequency;
	bool m_bCrossCorrelation;
	double m_fLaserFreq;
	double m_fTemperature;

	bool m_bCorrectTemperature;
	double m_fCorrectTemperature;
	bool m_bUseCommutator;

};



class CROAMolecule {
public:
	CROAMolecule() {
	}

	~CROAMolecule() {
	}

	std::vector<double> m_faCharge;           // Unit: e
	std::vector<CxDVector3> m_vaElDip;        // Unit: Debye
	std::vector<CxDMatrix3> m_maElQuad;       // Unit: Debye*pm
	std::vector<CxDVector3> m_vaElCurr;       // Unit: MB/pm
	std::vector<CxDVector3> m_vaMagDip;       // Unit: MB

	std::vector<CxDMatrix3> m_maPolElDip;     // Unit: e*pm^2/V
	std::vector<CDTensor3> m_taPolElQuad;     // Unit: e*pm^3/V
	std::vector<CxDMatrix3> m_maPolMagDip;    // Unit: e*pm^3/(fs*V)

	std::vector<double> m_faDCharge;          // Unit: e/fs
	std::vector<CxDVector3> m_vaDElDip;       // Unit: Debye/fs
	std::vector<CxDMatrix3> m_maDElQuad;      // Unit: Debye*pm/fs
	std::vector<CxDVector3> m_vaDElCurr;      // Unit: MB/(pm*fs)
	std::vector<CxDVector3> m_vaDMagDip;      // Unit: MB/fs

	std::vector<CxDMatrix3> m_maDPolElDip;    // Unit: e*pm^2/(V*fs)
	std::vector<CDTensor3> m_taDPolElQuad;    // Unit: e*pm^3/(V*fs)
	std::vector<CxDMatrix3> m_maDPolMagDip;   // Unit: e*pm^3/(V*fs^2)
};



class CROATrajectory {
public:
	CROATrajectory();
	~CROATrajectory();
	bool OpenFile();

	CTimeStep* GetTimeStep(int i);
	bool ReadTimeStep();
	bool SkipTimeStep();
	void CalcVelocities();

	void Rewind();

	void SetHistory(int i);

	std::string m_sFileName;
	std::vector<CTimeStep*> m_oaTSHistory;
	int m_iTSPos;
	std::vector<CROAMolecule*> m_oaMolecules;
	std::vector<CxDVec3Array> m_vaAtomVelocities;
	std::vector<CxDVec3Array> m_vaCOMCoord;

	FILE *m_pFile;
	bool m_bBQB;
	bool m_bVoronoi;
	CBQBFile *m_pBQB;

	int m_iHistory;
};



class CROAAtom {
public:
	std::vector<CxDVector3>  m_vaPos;
	std::vector<CxDVector3>  m_vaVel;

	std::vector<double> m_faCharge;           // Unit: e
	std::vector<CxDVector3> m_vaElDip;        // Unit: Debye
	std::vector<CxDMatrix3> m_maElQuad;       // Unit: Debye*pm
	std::vector<CxDVector3> m_vaElCurr;       // Unit: MB/pm
	std::vector<CxDVector3> m_vaMagDip;       // Unit: MB

	std::vector<CxDMatrix3> m_maPolElDip;     // Unit: e*pm^2/V
	std::vector<CDTensor3> m_taPolElQuad;     // Unit: e*pm^3/V
	std::vector<CxDMatrix3> m_maPolMagDip;    // Unit: e*pm^3/(fs*V)

	std::vector<double> m_faDCharge;          // Unit: e/fs
	std::vector<CxDVector3> m_vaDElDip;       // Unit: Debye/fs
	std::vector<CxDMatrix3> m_maDElQuad;      // Unit: Debye*pm/fs
	std::vector<CxDVector3> m_vaDElCurr;      // Unit: MB/(pm*fs)
	std::vector<CxDVector3> m_vaDMagDip;      // Unit: MB/fs

	std::vector<CxDMatrix3> m_maDPolElDip;    // Unit: e*pm^2/(V*fs)
	std::vector<CDTensor3> m_taDPolElQuad;    // Unit: e*pm^3/(V*fs)
	std::vector<CxDMatrix3> m_maDPolMagDip;   // Unit: e*pm^3/(V*fs^2)
};



class CROAEngine {
public:
	CROAEngine();
	~CROAEngine();

	bool Parse();
	void MainLoop();
	void MainLoop_Combined();
	void MainLoop_Gather();
	void MainLoop_Analyze();
	bool ReadStep(bool verbose);
	bool SkipStep();
	void Finish();
	void CalcVelocities();

	void ParseObservations();

	bool CheckTrajColumns(CROATrajectory *tr);

	void ComputeACF(CROAObservation *obs);

	void ComputeACFPair(CROAObservation *obs, CROAMolecule *mol1, CROAMolecule *mol2, std::vector<double> *wfn=NULL);

	void ComputeSpectrum(CROAObservation *obs);
	void ComputeSpectrum(CROAObservation *obs, const std::vector<double> &acf, std::vector<double> &spectrum, const char *name, bool im, bool ghack=false);
	void WriteSpectrum(CROAObservation *obs, const std::vector<double> &spectrum, const char *name);
	void WriteSpectrumSFG(CROAObservation *obs, const std::vector<double> &spectrumre, const std::vector<double> &spectrumim, const char *name);

	CTimeStep* GetROAStep(int traj, int step);


	bool PrepareVori(CTimeStep *ts, bool dipole, bool quadrupole, bool magmom);
	bool PrepareCurrentDensity(CTimeStep *ts, CTimeStep *tsref, bool first);
	bool CalculateCurrentDensity(CTimeStep *ts1, CTimeStep *ts2, CTimeStep *ts3, bool reset);

	bool m_bAniso;
	bool m_bCentral;
	double m_fFieldStrength; // In atomic units: Eh / (e * a0) = 5.14220652E11 V/m

	bool m_bSmoothData;

	bool m_bReplaceOutliers;
	bool m_bReplaceOutliersElectric;
	bool m_bReplaceOutliersMagnetic;
	double m_fReplaceOutliersThreshold;

	std::vector<int> m_iaOutliers[8];
	double m_fSmoothWaveNumber;

	bool m_bDumpMolecularProps;
	bool m_bDumpMol1Props;
	bool m_bWriteACFs;

	bool m_bGatherWriteCSV;
	bool m_bGatherCheck;

	bool m_bReverseTraj;

	int m_iMode;
	bool m_bPattern3;


	bool m_bIR;
	bool m_bRaman;
	bool m_bVCD;
	bool m_bROA;
	bool m_bSFG;
	bool m_bMagMom;
	bool m_bPola;

	int m_iStepsProcessed;

	std::vector<CROATrajectory*> m_oaTrajectories;
	std::vector<CROAObservation*> m_oaObservations;

	std::vector<CROAAtom> m_oaAtoms;

	std::vector<int> m_iaTrajKinds;
	std::vector<int> m_iaInvTrajKinds;

	C3DF<VORI_FLOAT> m_pDensity;
	C3DF<VORI_FLOAT> m_pDensityGrad[3];

	double *m_pPDESolution;

	FILE *m_fPDEInfo;

	bool m_bPDEFastMode;

	std::vector<FILE*> m_fMolIntegralFiles;

	CCrossCorrelation *m_pCrossCorr;
	CxDoubleArray m_faInput1;
	CxDoubleArray m_faInput2;
	CxDoubleArray m_faOutput;

	bool m_bPDERestart;
};


class CFitParabola {
public:

	// Finds a best-fit parabola through four given points
	void Init4(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

	double Evaluate(double x);

private:
	double m_fA;
	double m_fB;
	double m_fC;
};



#endif



