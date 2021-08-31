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


#ifndef TIMESTEP_H
#define TIMESTEP_H


// This must always be the first include directive
#include "config.h"

#include <math.h>
#include <stdio.h>
#include "xobject.h"
#include "xobarray.h"
#include "xdoublearray.h"
#include "xwordarray.h"
#include "xdvector3.h"
#include "xdmatrix3.h"
#include "nbsearch.h"
#include "backtrace.h"
#include "statistics.h"
#include "xbytearray.h"
#include "atomgroup.h"
#include "moltools.h"
#include "internalcoord.h"
#include "xptrarray.h"
#include "xquaternion.h"
#include "xdvec3array.h"
#include "3df.h"
#include "bqb.h"
#include "posdomain.h"


class CTimeStep : public CxObject
{
public:
	void DumpDipoles();
	void WritePOV(const char *s);
	void FoldAtomsPositive();
	double FoldedDistance(int i1, int i2);
	void CenterCOM();
	int REC_UniteMolecules(CSingleMolecule *sm, int i0, int depth);
	bool ReadCellVector(FILE *a);
	void WriteMol2(FILE *a);
	void CalcMinMax();

	bool SkipTimestep(FILE *a);

	bool ReadTimestep(FILE *a, bool needinfo);
	bool ReadTimestep(CxMemFile *file);

	bool ReadMol2(FILE *a, bool needinfo);
	bool SkipMol2(FILE *a);

	bool ReadPDB(FILE *a, bool needinfo, CxDVec3Array *v);
	bool SkipPDB(FILE *a);

	bool ReadLAMMPS(FILE *a, bool needinfo);
	bool SkipLAMMPS(FILE *a);
	bool AnalyzeLAMMPSColumns( const char *s );

	bool ReadDLPOLY(FILE *a, bool needinfo);
	bool SkipDLPOLY(FILE *a);

	bool ReadAmber(FILE *a, bool needinfo);
	bool SkipAmber(FILE *a);

	bool ReadXYZ(FILE *a, bool needinfo, CxDVec3Array *v);
	bool SkipXYZ(FILE *a);

	bool ReadGRO(FILE *a, bool needinfo);
	bool SkipGRO(FILE *a);

	bool ReadDCDHeader(FILE *a);
	bool ReadDCD(FILE *a, bool needinfo);
	bool SkipDCDHeader(FILE *a);
	bool SkipDCD(FILE *a);

	bool ReadVoronoi(FILE *a, bool needinfo);
	bool SkipVoronoi(FILE *a);

	bool ReadCube(FILE *a, bool needinfo);
	bool ReadCube(CxMemFile *file);
	bool SkipCube(FILE *a);

	bool ReadBQB(CBQBFile *source, bool needinfo);
	bool ReadBQBCubeStep(bool needinfo);
	bool ReadBQBTrajStep(const std::vector<unsigned char> *data, bool needinfo);
	bool SkipBQB();
	bool SkipBQB(CBQBFile *source);
	bool FetchBQBCompTrajFrame(bool needinfo);

	void DoubleBox();
	void DoubleBoxVelocity();
	void DoubleBoxForce();

	CTimeStep();
	~CTimeStep();

	bool ScanWannier(bool verbose);
	void CalcMagneticDipoles();
	void CalcDipoles(bool verbose);
	void CalcPolarizabilities();
	long ExtractNumber(int i);
	int GetCommentNumberCount();
//	double MolDist(CSingleMolecule *ref, CSingleMolecule *sm2, CNbSearch *nb);
//	bool CheckIfOrdered();
	bool ScanMolecules();
	void CalcCenters();
	void CalcCenterVelocities();
	void RECURSION_ScanMolecules(int i, CxByteArray *ta, CSingleMolecule *sm, int depth, int *stack, unsigned long bmask, bool w);
	void PrintMegaTree();
	void PrintMatrix(bool onlyfirst, bool onlybind);
	void RECURSION_MegaTree(int i, char* ta, int depth, unsigned long bmask, bool w, int *stack);
	bool BondRange(int i1, int i2, double *f);
	bool MirrorBond(int i1, int i2);
	void FoldMolecules();
	void FoldAtoms();
	void UniteMolecules(bool verbose);
	void CenterPos(const CxDVector3 &vec);
//	void Transform(const CxDMatrix3 &mat);
	void Transform(const CxDMatrix3 &mat);
	void Transform(const CxQuaternion &q);
	void SubVelocities(const CxDVector3 &vec);
	bool ReadTimestepVel(FILE *a);
	bool ReadTimestepForce(FILE *a);
	void WriteTimestep(FILE *a);
	void WriteTimestepNb(FILE *a, CNbSet *nbs, int singlemol = -1);
	void ExportSingleMolecule_PDB(CSingleMolecule *sm, const char *s);
//	void WriteTimestepNb(int refmol, FILE *a);
//	void GatherNbDiagram(int refmol, CNbSearch *nb);
	void AddAtoms();
	void BuildAtomIndex();
	void CopyFrom(CTimeStep *t);

//	void ScanNeighborhood(int fixmol, int refmol, CNbSearch *nb, CNbSearch *prev);
/*	void ScanAngle(int fixmol, int refmol, CCondition *co, CNbSearch *prev);
	void ScanNewNeighborhood(int fixmol, int refmol, CNbSearch *nb, CNbSearch *prev);*/

	bool m_bAbortRing;

	CxDVector3 m_vMin, m_vMax;

	CxDVec3Array m_vaCoords;
	CxDVec3Array m_vaCoords_Unfolded;
	CxDVec3Array m_vaCoords_Original;
	CxDVec3Array m_vaVelocities;
	CxDVec3Array m_vaForces;
	unsigned long m_iConnectedAtoms;
	unsigned long m_iGesAtomCount;
	unsigned long m_iGesAtomModulo;
	char *m_pComment;
//	char *m_pLabels;
	CxPtrArray m_paLabels;
	CxPtrArray m_paMol2Types;

	CxDoubleArray m_faVolume;                 //            Unit: pm^3
	CxDoubleArray m_faCharge;                 //            Unit: e
	CxDoubleArray m_faCoreCharge;             //            Unit: e
	CxDVec3Array m_dipoleMoments;             // Ref: Atom, Unit: e*pm
	CxDVec3Array m_chargeCenters;             // Ref: Atom, Unit: pm
	CxDVec3Array m_totalCurrents;             //            Unit: MB/pm
	CxDVec3Array m_magneticDipoleMoments;     // Ref: Atom, Unit: MB

	std::vector<CxDMatrix3> m_maQuadTensor;   // Ref: Atom, Unit: e*nm^2  (!!!)

//	int m_iSizeBytes;
	C3DF<VORI_FLOAT> *m_pVolumetricData;
	C3DF<VORI_FLOAT> *m_pVolumetricDataTimeDev;
	CxDoubleArray *m_pCurrentDensity;

	CPosDomainEngine *m_pPosDomainEngine;
	
private:

	static bool m_bFirst;


};

#endif
