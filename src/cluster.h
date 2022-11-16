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


#ifndef CLUSTER_H
#define CLUSTER_H


// This must always be the first include directive
#include "config.h"

#include <math.h>
#include <stdio.h>
#include "xobarray.h"
#include "xdvec3array.h"
#include "xdoublearray.h"
#include "xlongarray.h"
#include "df.h"
#include "timestep.h"
#include "random.h"


CxDVector3 FoldVector(CxDVector3 v);
extern CxObArray g_oaSingleMolecules;
extern CxObArray g_oaMolecules;
extern CxIntArray g_laAtomSMIndex;


class CClusterAnalysis;


class CClusterNode : public CxObject
{
public:
	CClusterNode();
	~CClusterNode();

	double m_fPosX;
	int m_iParent;
//	CxLongArray m_laChildren;
	int m_iChildren[2];
	CxDVec3Array m_vaAtoms;
	CxIntArray m_iaAtoms;
	double m_fDistance;
	double m_fNNDistance;
	int m_iNextNeighbor;
	int m_iIndex;
	int m_iMonomers;
	CClusterNode *m_pParent;
};


class CExtendedCluster : public CxObject
{
public:
	bool IsHydrogenBond(int i1, int i2);
	int CountDifferentAtomCodesUndir();
	int CountDifferentAtomCodes();
	void BuildAtomCodesUndir();
	void BuildAtomCodes();
	void CreateFrom(CClusterNode *n, CTimeStep *ts, CClusterAnalysis *ca, bool wrap);
	CExtendedCluster();
	~CExtendedCluster();

	int m_iClusterIndex;
	int m_iIndex;
	int m_iInterBonds;
	int m_iMonomers;
	int m_iCount;
	int m_iCountUndir;
	double m_fSignificance;
	double m_fSignificanceUndir;
	CxIntArray m_iaAtoms;
	CxObArray m_oaBonds;
	CxObArray m_oaRelBonds;
	CxIntArray m_iaMonomerSM;
	CxObArray m_oaRelBondsUndir;
	CxDoubleArray m_faAtomCodes;
	CxDoubleArray m_faTempAtomCodes;
	CxDoubleArray m_faAtomCodesUndir;
	CxDoubleArray m_faTempAtomCodesUndir;
	CxDoubleArray m_faCountAC;
};


class CClusterAnalysis : public CxObject
{
public:
	void REC_DumpAgr(CClusterNode *n, CGrace *g, int ec, unsigned long color);
	unsigned int POVColorFunction(int i);
	void WriteClusterXYZ(CClusterNode *n, CTimeStep *ts, FILE *a);
	void RenderStepPOV(CTimeStep *ts, const char *s);
	void REC_IsConnected(CExtendedCluster *ec, int i);
	bool IsConnected(CExtendedCluster *ec);
	void RenderFormula(const char *s);
	void DumpAtomDOT(CExtendedCluster *ec, const char *s);
	void DumpSchematicDirectedDOT(CExtendedCluster *ec, const char *s);
	void DumpSchematicUndirectedDOT(CExtendedCluster *ec, const char *s);
	bool AddToClusterTopo(CExtendedCluster *ec, int *i);
	bool AddToClusterTopoUndir(CExtendedCluster *ec, int *i);
	void WriteOutput(const char *multibuf);
	void BuildClusterDistribution();
	void BuildDistCache(CTimeStep *ts);
	void BinPoly();
	void CleanUp();
	void Process(CTimeStep *ts);
	void Create();
	void BinDistances(CTimeStep *ts);
	void REC_DumpPoints(CxDoubleArray *fa, CClusterNode *n);
	void REC_AvgX(CClusterNode *n);
	void REC_SortX(double pos, int depth, CClusterNode *n);
	void DumpAgr(const char *s);
	void REC_DumpDot(FILE *a, CClusterNode *n);
	void DumpDot(const char *s);
	void FindNearestNeighbor(CClusterNode *n);
	void TraceNeighbors();
	void BuildTree();
	CClusterAnalysis();
	~CClusterAnalysis();
	void Parse();
	void AddCluster(int i);
	void AddParticle(double x, double y, double z);
	void AddParticle(int i);


	double Distance(CxDVector3 *a, CxDVector3 *b)
	{
		return FoldVector(*a - *b).GetLength();
	}


	double Distance_Cache(int a, int b)
	{
		return m_pDistCache[a*m_iAtomListSize+b];
	}


	int MassCenter(int i)
	{
		return ((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[g_laAtomSMIndex[i]])->m_oaAtomOffset[((CSingleMolecule*)g_oaSingleMolecules[g_laAtomSMIndex[i]])->m_oaAtomOffset.GetSize()-1])->GetAt(1);
	}


	double m_fCSDFThresh;
	bool m_bUseBondCutoffs;
	bool m_bCOMTrick;
	bool m_bPOVDiagram;
	int m_iPOVMonomersMin;
	int m_iPOVMonomersMax;
	double m_fPOVAngleIncrement;
	double m_fPOVAngle;
	int m_iPOVColCount;
	int m_iPOVFrameCounter;
	CxObArray m_oaPOVSMColorHistory;
	FILE *m_fPOVScript;
	CxObArray m_oaPOVExtendedCluster;
	CxIntArray m_iaPOVSMColor;
	CxIntArray m_iaPOVSMCluster;
	bool m_bPOVTrajectory;
	bool m_bClusterTopoTrajectories;
	CxIntArray m_iaConnectedBuffer;
	bool m_bSelectBonds;
	CxIntArray m_iaAtomBondFrom;
	CxIntArray m_iaAtomBondTo;
	int m_iClusterTopoTries;
	CxIntArray m_iaClusterTopoSMTemp;
	CxObArray m_oaClusterTopo;
	CxObArray m_oaClusterTopoUndir;
	bool m_bClusterTopoAllBonds;
	bool m_bClusterTopoDrawAtom;
	bool m_bClusterTopoDrawDirected;
	bool m_bClusterTopoDrawUndirected;
	double m_fClusterTopoHBThres;
	double m_fClusterTopoDispThres;
	CxIntArray m_iaClusterTopoMolColor;
	bool m_bClusterTopo;
	int m_iClusterTopoMax;
	int m_iCounter;
	int m_iMolecule;
	CxObArray m_oaAtomGroups;
	int m_iResX;
	int m_iResY;
	double m_fMaxDist;
	int m_iRes;
	bool m_bAnim;

	double *m_pDistCache;
	int *m_pAtomIndex;
	CxIntArray m_iaAtomList;
	int m_iAtomListSize;
	bool m_bDistCache;
	bool m_b2DPlots;
	bool m_bDiffPlot;
	double m_fDiffInterval;

	CDF *m_pClusterDistanceDF;
	CDF *m_pPolymerDF;
	CDF *m_pClusterCountDF;
	CDF *m_pClusterDistributionDF;
	CDF *m_pClusterDistributionDF_Mod;
//	CDF *m_pClusterDistribution2DF;
	CDF *m_pClusterSizeDF;

	C2DF *m_pClusterDistance2D;
	C2DF *m_pClusterDistribution2D;
	C2DF *m_pClusterSize2D;
	C2DF *m_pClusterCount2D;
	C2DF *m_pClusterDiff2D;

	CxObArray *m_poaClusterPoly;
	int m_iMonomers;
	FILE *m_fAnim;
	CClusterNode *m_pTop;
	CxObArray m_oaClusters;
	CxObArray m_oaBaseClusters;
	CxObArray m_oaTopClusters;

	CxDoubleArray m_faMinDist;
	CxDoubleArray m_faMaxDist;
	CxDoubleArray m_faHetNumber;
	CxDoubleArray m_faIntHetNumber;
	CxDoubleArray m_faIntHetNumber2;
	CxDoubleArray m_faCenter;

	int m_i2DResX;
	int m_i2DResY;
	int m_i2DStride;

	bool m_bHetMeasure;
	bool m_bIdealGas;

	bool m_bCutClusters;
	int m_iCutClusterMaxSize;
	int *m_pCutClusterCounter;
	double *m_pCutClusterDifference;
	FILE **m_pCutClusterFiles;

	CRandom *m_pRandom;

};

#endif


