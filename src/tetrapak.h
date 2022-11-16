/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Martin Brehm and Martin Thomas.

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


#ifndef TETRAPAK_H
#define TETRAPAK_H


// This must always be the first include directive
#include "config.h"

#include "xobject.h"
#include "xobarray.h"
#include "xdvec3array.h"
#include "timestep.h"
#include "3df.h"


#define VORI_EPS 1.0E-11


extern int glv_iRayNonConv;
extern int glv_iRayFail;
extern bool glv_bJitter;
extern int glv_iJitterSeed;
extern double glv_fJitterAmplitude;



void DumpMolecularProps();



class CTetraPak : public CxObject
{
public:
	CTetraPak();
	~CTetraPak();
	void Parse();
	bool ParseSilent(CTimeStep *ts);
	bool BuildVoronoi(CTimeStep *ts, int verbose, bool sanity, bool sanityall=false);
	void BuildVoronoiBuffer(CTimeStep *ts);
// 	VORI_FLOAT Integrate(C3DF<VORI_FLOAT> *df, double mi[2], double ma[2]);
	VORI_FLOAT Integrate_Refine(C3DF<VORI_FLOAT> *df, double mi[2], double ma[2]);
	void Integrate_Verbose(C3DF<VORI_FLOAT> *df, double mi[2], double ma[2]);
	void ProcessStep(CTimeStep *ts, int verbose);
	
	CxObArray m_oaFaces;
	CxObArray m_oaVoroBuffer;
	C3DF<VORI_FLOAT> *m_p3DF;
	int m_fRayCount;
	int m_fRayHisto[4];
	int *m_hitCount;
	int m_iFaces;
	CxDoubleArray m_faCharge;
	CxDoubleArray m_faVolume;
	CxDoubleArray m_faCoreCharge;
	bool m_bVoronoiCharges;
	FILE *m_fCubePipe;
	bool m_bInterpolation;
	int m_iInterpolationFactor;
//	bool m_bInterpolationLinear;
	int m_iInterpolationFactorCube;


	double m_totCharge;
	double m_totChargeElec;
	double m_totChargeCore;
	double m_totDip[3];
	double m_totDipElec[3];
	double m_totDipCore[3];
	double m_totQuad[9];
	double m_totQuadElec[9];
	double m_totQuadCore[9];
	double m_totCurrent[3];
	double m_totCurrentElec[3];
	double m_totCurrentCore[3];
	double m_totMagMom[3];
	double m_totMagMomElec[3];
	double m_totMagMomCore[3];
	double m_qRef[3];

	bool m_bSaveAtomIntegrals;

private:
	void integrateCell(C3DF<VORI_FLOAT> *electronDensity, CxDoubleArray *currentDensity, double mi[2], double ma[2], double refPoint[3], double *charge, CxDVector3 *dipoleMoment, CxDMatrix3 *quadrupoleMoment, CxDVector3 *totalCurrent, CxDVector3 *magneticMoment, bool sanity = false/*, int debug = 0*/);
	
// 	CxDVector3 integrateMoment(C3DF<VORI_FLOAT> *df, double mi[2], double ma[2], double refPoint[3]);
// 	CxDVector3 integrateTotalCurrent(C3DF<VORI_FLOAT> *df, CxFloatArray *currentDensity, double mi[2], double ma[2]);
// 	CxDVector3 integrateMagneticMoment(C3DF<VORI_FLOAT> *df, CxFloatArray *currentDensity, double mi[2], double ma[2], double refPoint[3]);
	
	CxDVec3Array m_moments;
	CxDVec3Array m_chargeCenters;
	std::vector<CxDMatrix3> m_maQuadTensor;
	CxDVec3Array m_totalCurrent;
	CxDVec3Array m_magneticMoments;
	bool m_saveTotalIntegrals;
	FILE *m_totalIntegralFile;

	FILE *m_fAtomIntegralFile;
// 	bool m_saveAtomCharges;
// 	FILE *m_atomChargeFile;
// 	bool m_saveMolCharges;
// 	FILE *m_molChargeFile;

};


class CTetraFace : public CxObject
{
public:
	CTetraFace()
	{
	}


	~CTetraFace()
	{
	}


	bool XRay_Hit(double py, double pz, double &px/*, bool verbose=false*/);
// 	{
// 		double tf, a, b, ty, tz;
// 		int z;
// 
// 		#define lvec1x (m_vSpan1[0])
// 		#define lvec1y (m_vSpan1[1])
// 		#define lvec1z (m_vSpan1[2])
// 		#define lvec2x (m_vSpan2[0])
// 		#define lvec2y (m_vSpan2[1])
// 		#define lvec2z (m_vSpan2[2])
// 
// 		// Hack from 09.10.2014: Ray cannot hit walls which are almost parallel...
// //		if (fabs(m_vNormal[0]) < 1.0E-14)
// //			return false;
// 
// 		ty = py - m_vCenter[1];
// 		tz = pz - m_vCenter[2];
// 		if (verbose)
// 			mprintf(GREEN, "%.20g %.20g\n", ty, tz);
// 		
// 		tf = -lvec1z*lvec2y + lvec1y*lvec2z;
// 
// 		bool run;
// 		do {
// 			run = false;
// 			
// 			a = (-tz*lvec2y + ty*lvec2z)/tf;
// 			b = ( tz*lvec1y - ty*lvec1z)/tf;
// 			
// // 			if (verbose)
// // 				mprintf(GREEN, "%.20g %.20g\n", a, b);
// 
// #ifdef TARGET_WINDOWS
// 			if (!_finite(a) || !_finite(b))
// 				return false;
// #else
// 			if (!isfinite(a) || !isfinite(b))
// 				return false;
// #endif
// 
// 			for (z=0;z<m_vaEdges.GetSize();z++) {
// 				double area = m_vaEdges[z][0]*a + m_vaEdges[z][1]*b + m_vaEdges[z][2];
// 				double dist = 2.0 * area / m_faEdgesLength[z];
// // 				if (verbose)
// // 					mprintf(GREEN, "%d %.20g\n", z, dist);
// 				if (fabs(dist) < VORI_EPS) {
// // 					if (verbose)
// // 						mprintf(GREEN, "repeat %d\n", z);
// 					ty += VORI_EPS;
// 					tz += VORI_EPS;
// 					run = true;
// 					break;
// 				}
// 				if (area * m_vaEdges[z][2] < 0)
// // 				if ((m_vaEdges[z][0]*a + m_vaEdges[z][1]*b + m_vaEdges[z][2])*m_vaEdges[z][2] < 0)
// 					return false;
// 			}
// 		} while (run);
// 
// 		px = a*lvec1x + b*lvec2x + m_vCenter[0];
// 
// /*
// #ifdef TARGET_WINDOWS
// 		if (!_finite(px))
// #else
// 		if (!isfinite(px))
// #endif
// 		{
// 			mprintf("\nXRay_Hit(py=%.10f, pz=%.10f) returned %.10f:\n",py,pz,px);
// 			mprintf("  ty=%.20G  tz=%.20G\n",ty,tz);
// 			mprintf("  m_vSpan1: %f  %f  %f\n",m_vSpan1[0],m_vSpan1[1],m_vSpan1[2]);
// 			mprintf("  m_vSpan2: %f  %f  %f\n",m_vSpan2[0],m_vSpan2[1],m_vSpan2[2]);
// 			mprintf("  m_vNormal: %f  %f  %f\n",m_vNormal[0],m_vNormal[1],m_vNormal[2]);
// 			mprintf("  CrossP: %10G\n",CrossP(m_vSpan1,m_vSpan2).GetLength());
// 			mprintf("  m_vCenter: %f  %f  %f\n",m_vCenter[0],m_vCenter[1],m_vCenter[2]);
// 			mprintf("  tf=%.20G,  a=%f,  b=%f.\n",tf,a,b);
// 		}*/
// 
// /*		if (verbose)
// 		{
// 			mprintf("\nXRay_Hit(py=%.10f, pz=%.10f) returned %.10f:\n",py,pz,px);
// 			mprintf("  ty=%.20G  tz=%.20G\n",ty,tz);
// 			mprintf("  m_vSpan1: %f  %f  %f\n",m_vSpan1[0],m_vSpan1[1],m_vSpan1[2]);
// 			mprintf("  m_vSpan2: %f  %f  %f\n",m_vSpan2[0],m_vSpan2[1],m_vSpan2[2]);
// 			mprintf("  m_vNormal: %g  %g  %g\n",m_vNormal[0],m_vNormal[1],m_vNormal[2]);
// 			mprintf("  CrossP: %10G\n",CrossP(m_vSpan1,m_vSpan2).GetLength());
// 			mprintf("  m_vCenter: %f  %f  %f\n",m_vCenter[0],m_vCenter[1],m_vCenter[2]);
// 			mprintf("  tf=%.20G,  a=%f,  b=%f.\n",tf,a,b);
// 		}*/
// 
// 		return true;
// 	}


	CxDVector3 m_vCenter;
	CxDVector3 m_vSpan1;
	CxDVector3 m_vSpan2;
	CxDVector3 m_vNormal;

	CxDVec3Array m_vaEdges;
	CxDoubleArray m_faEdgesLength;

	double m_fMin[2];
	double m_fMax[2];
};


class CTetraCellBuffer : public CxObject
{
public:
	CTetraCellBuffer()
	{
	}

	~CTetraCellBuffer()
	{
	}

	int m_iFaces;
	CxDVector3 m_vCenter;
	CxObArray m_oaFaces;
};


class CTetraFaceBuffer : public CxObject
{
public:
	CTetraFaceBuffer()
	{
	}

	~CTetraFaceBuffer()
	{
	}

	CxDVec3Array m_vaVertices;
};

void printCubeMemFileName(char *name, int length, bool nextFile);

#endif
