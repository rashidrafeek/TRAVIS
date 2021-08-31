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


#ifndef VOROWRAPPER_H
#define VOROWRAPPER_H


// This must always be the first include directive
#include "config.h"

#include "tools.h"
#include "xobject.h"
#include "voro++.h"
#include "xobarray.h"
#include "df.h"
#include "xdmatrix3.h"
#include "xwordarray.h"
#include "xdoublearray.h"
#include "xintarray.h"
#include "xstring.h"
#include "xlongarray.h"
#include <string>


class CTimeStep;
class CVoroWrapper;
class CAtomGroup;



class CVoroNbMatrix {
public:

	bool Parse( int i );

	void Finish();

	void BuildNeighborMatrices();

	bool m_bJoinEquivalent;

	bool m_bInter;
	bool m_bIntra;

	bool m_bExclude12;
	bool m_bExclude13;
	bool m_bExclude14;

	std::vector<CAtomGroup*> m_oaRows;
	std::vector<CAtomGroup*> m_oaColumns;

	std::vector<CxString> m_saRowLabels;
	std::vector<CxString> m_saColumnLabels;

	std::vector<int> m_iaRowIndex;
	std::vector<int> m_iaColumnIndex;

	std::vector<int> m_iaRowMolecule;
	std::vector<int> m_iaColumnMolecule;

	std::vector<int> m_iaRowAtomCount;
	std::vector<int> m_iaColumnAtomCount;

	std::vector<int> m_iaRowRealElement;
	std::vector<int> m_iaColumnRealElement;

	std::vector<double> m_iaTempMatrix;
	std::vector<double> m_faMatrix;
	std::vector<double> m_faCounter;

	std::vector<std::vector<int> > m_iaaNeighborMatrices;
	std::vector<int> m_iaLocalAtomIndex;

	int m_iIndex;
	int m_iRows;
	int m_iCols;
};



class CVoroAtom : public CxObject
{
public:
	void Dump(const char *s);
	void InitAnalyses();
	CVoroAtom();
	~CVoroAtom();

	CVoroWrapper *m_pParent;

	int m_iMolecule;
	int m_iAtomType;
	int m_iRealAtomType;
	int m_iAtom;

	int *m_pNbhTempMol;

	CDF *m_pVolume;
	CDF *m_pSurfaceArea;
	CDF *m_pExposedSurface;
	CDF *m_pExposedSurfacePerc;
	CDF *m_pAVRatio;
	CDF *m_pFaces;
	CDF *m_pFaceOrders;
	CDF *m_pMaxRadius;
	CDF **m_pNbhAtoms;
	CDF **m_pNbhMolecules;
	CDF **m_pNbhDistAtoms;
	CDF **m_pNbhDistMolecules;
};


class CVoroMolecule : public CxObject
{
public:
	void Add(CVoroMolecule *m);
	void Dump(const char *s);
	void InitAnalyses();
	CVoroMolecule();
	~CVoroMolecule();

	CVoroWrapper *m_pParent;

	int m_iMolecule;
	int m_iSingleMol;
	int m_iCenterOffset;

	int *m_pNbhTempMol;
	int *m_pNbhTempAtoms;

	double tempvol;
	double tempsurf;
	int tempfaces;
	double *m_pNbhTempArea;

	CDF *m_pVolume;
	CDF *m_pSurfaceArea;
	CDF *m_pAVRatio;
	CDF *m_pFaces;
	CDF **m_pNbhAtoms;
	CDF **m_pNbhMolecules;
	CDF **m_pNbhDistAtoms;
	CDF **m_pNbhDistMolecules;
	CDF **m_pNbhAreaMolecules;
};




class CVoroWrapper : public CxObject
{
public:
	bool m_bIncludeWritten;
	bool m_bWritePOVMovie;
	CxDVector3 m_vPOVFaceBleach;
	double m_fPOVFaceBleach;
	double m_fPOVNbHighlightFac;
	double m_fPOVRayTrans;
	double m_fPOVNbTrans;
	double m_fPOVFaceTrans;
	double m_fPOVEdgeTrans;
	CxDoubleArray m_faVoroMetric;
	bool m_bVoroMetricCDF;
	bool m_bVoroMetric;
	C2DF *m_pVoroMetricCDF2;
	C2DF *m_pVoroMetricCDF;
	int m_iPOVThreads;
	int m_iPOVThreadPos;
	FILE *m_fPOVCT;
	double m_fPOVAtomCenter;
	double m_fPOVScale;
	bool m_bPOVMolecules;
	double m_fPOVExplode;
	double m_fPOVClip;
	CxDVector3 m_vPOVBoxClipLow;
	CxDVector3 m_vPOVBoxClipHigh;
	int m_iPOVBoxTotalSteps;
	int m_iPOVBoxTotalFrames;
	int m_iPOVCurrentChoreographyPart;


	double m_fPOV_FPS;
	CxObArray m_oaBoxChoreography;
	double m_fPOVBoxOuterOpacity;
	double m_fPOVZoom;
	int m_iPOVFrameCounter;
	char m_sPOVText[256];
	double m_fPOVAngle;
	CxLongArray m_iaPOVBoxAtomColors;
	FILE **m_fPOVBoxScript;
	bool m_bPOVBox;
	void FinishSurfCover(const char *s);
	void ProcessSurfCover(CTimeStep *ts);
//	void WriteXYZCell(CTimeStep *ts, int mol, int smol);
//	void WriteXYZCell_Start(CTimeStep *ts, const char *s, int mol, int smol);
	void WriteMoleculeInfo(const char *s);
	void WriteAtomInfo(const char *s);
	void WriteNbAtomDistMatrix(const char *s);
	void WriteNbAtomCountMatrix(const char *s);
	void WriteNbMoleculeDistMatrix(const char *s);
	void WriteNbMoleculeAreaMatrix(const char *s);
	void WriteNbMoleculeCountMatrix(const char *s);
	double m_fBoxDens;
	void Finish();
	void Dump(const char *s, CTimeStep *ts);
	void DumpPlastic(FILE *s, CTimeStep *ts);
	void Parse();
	void Init();
	void Init2();
	void Build(CTimeStep *ts);
	CVoroWrapper();
	~CVoroWrapper();
	int m_iBlocksX;
	int m_iBlocksY;
	int m_iBlocksZ;
	container_periodic_poly *m_pContainer;
	CxObArray m_oaVoroAtoms;
	CxObArray m_oaVoroMolecules;
	long *m_pAssignAtoms;
	long *m_pAssignMolecules;
	long *m_pAssignMoleculeTypes;
	double m_fMaxVol;
	double m_fMaxSurf;
	bool *m_pAtomTouched;
	bool m_bVoroStat;

	bool m_bVoroTimeSeries;
	std::vector<int> m_iaVoroTS;
	std::string m_sVoroTSHead;
	FILE *m_fTSVolume;
	FILE *m_fTSSurface;
	FILE *m_fTSFaces;
	FILE *m_fTSAVRatio;
	FILE *m_fTSMaxRadius;

	CxDoubleArray m_faPOVBoxAtomVisible;
	bool m_bWritePOV;
	int m_iPOVMol;
	int m_iPOVSM;
	bool m_bPOVEdges;
	bool m_bPOVFaces;
//	double m_fPOVFaceOpac;
	bool m_bPOVFaceColor;
	bool m_bPOVAtoms;
	bool m_bPOVRot;
	CxDMatrix3 m_mPOVMat;
	bool m_bPOVVertices;
	CxString m_sPOVExe;
	CxDoubleArray m_faPOVFaceColor;
	CxWordArray m_waPOVFaceColorMol;
	CxWordArray m_waPOVFaceColorElem;
	CxWordArray m_waPOVFaceColorAtom;
	FILE *m_fPOVScript;
	CxDVector3 m_vPOVRotInc;
	CxDVector3 m_vPOVRotPos;
	bool m_bPOVAtomGrey;
	double m_fPOVCameraDist;
	int m_iPOVResX;
	int m_iPOVResY;
	bool m_bPOVFaceColorRef;
	bool m_bPOVDrawNeighbors;
	bool m_bPOVNbColorFace;
	bool m_bPOVHighlightNbAtoms;

	bool m_bVoroNbMatrix;
	std::vector<CVoroNbMatrix*> m_oaNbMatrices;

	bool m_bSurfCover;
	int m_iSurfCoverMol;
//	int m_iSurfCoverSM;
	CxWordArray m_waSurfCoverSM;
	CxObArray m_oaSurfCoverData;
	CxWordArray m_waSurfCoverMol;
	CxWordArray m_waSurfCoverElem;
	CxWordArray m_waSurfCoverAtom;
};


#endif


