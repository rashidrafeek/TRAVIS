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


// This must always be the first include directive
#include "config.h"


//#define USE_CUBETOOL_FFTW


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>


#include "tools.h"
#include "2df.h"
#include "globalvar.h"


#ifdef USE_CUBETOOL_FFTW
	#include "fftw3.h"
#endif



const char *GetRevisionInfo_cubetool(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_cubetool() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



#define LEN_BOHR2ANGSTROM (0.52917721092)




class CAtomPos {
public:
	int m_iOrd;
	double m_fPos[3];
};



class CCube {
public:

	CCube() : m_pBin(NULL) {
		m_iRes[0] = 0;
		m_iRes[1] = 0;
		m_iRes[2] = 0;
		m_fStrideA[0] = 0;
		m_fStrideA[1] = 0;
		m_fStrideA[2] = 0;
		m_fStrideB[0] = 0;
		m_fStrideB[1] = 0;
		m_fStrideB[2] = 0;
		m_fStrideC[0] = 0;
		m_fStrideC[1] = 0;
		m_fStrideC[2] = 0;
		m_fCenter[0] = 0;
		m_fCenter[1] = 0;
		m_fCenter[2] = 0;
		m_iAtomCount = 0;
	}


	~CCube() {
		int z;
		if (m_pBin != NULL) {
			delete[] m_pBin;
			m_pBin = NULL;
		}
		for (z=0;z<(int)m_oaAtoms.size();z++)
			delete m_oaAtoms[z];
		m_oaAtoms.clear();
	}


	void CopyFrom( const CCube *cube ) {
		int z;
		CAtomPos *ap;
		if (m_pBin != NULL)
			delete[] m_pBin;
		for (z=0;z<(int)m_oaAtoms.size();z++)
			delete m_oaAtoms[z];
		m_oaAtoms.clear();
		m_iRes[0] = cube->m_iRes[0];
		m_iRes[1] = cube->m_iRes[1];
		m_iRes[2] = cube->m_iRes[2];
		m_fStrideA[0] = cube->m_fStrideA[0];
		m_fStrideA[1] = cube->m_fStrideA[1];
		m_fStrideA[2] = cube->m_fStrideA[2];
		m_fStrideB[0] = cube->m_fStrideB[0];
		m_fStrideB[1] = cube->m_fStrideB[1];
		m_fStrideB[2] = cube->m_fStrideB[2];
		m_fStrideC[0] = cube->m_fStrideC[0];
		m_fStrideC[1] = cube->m_fStrideC[1];
		m_fStrideC[2] = cube->m_fStrideC[2];
		m_fCenter[0] = cube->m_fCenter[0];
		m_fCenter[1] = cube->m_fCenter[1];
		m_fCenter[2] = cube->m_fCenter[2];
		m_iAtomCount = cube->m_iAtomCount;
		m_pBin = new double[m_iRes[0]*m_iRes[1]*m_iRes[2]];
		for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
			m_pBin[z] = cube->m_pBin[z];
		for (z=0;z<(int)cube->m_oaAtoms.size();z++) {
			ap = new CAtomPos();
			ap->m_iOrd = cube->m_oaAtoms[z]->m_iOrd;
			ap->m_fPos[0] = cube->m_oaAtoms[z]->m_fPos[0];
			ap->m_fPos[1] = cube->m_oaAtoms[z]->m_fPos[1];
			ap->m_fPos[2] = cube->m_oaAtoms[z]->m_fPos[2];
			m_oaAtoms.push_back(ap);
		}
	}


	bool ReadCube(const char *s);

	bool WriteCube(const char *s);

	bool ImportBinary( const char *s );

	bool ImportBinary4D( const char *s, int totalframes, int frame );

	bool AnalyzeBinary4D( const char *s, int totalframes );

	bool AddCube( CCube *tcu );

	bool SubtractCube( CCube *tcu );

	void Add( double f );

	void LimitRange( double mi, double ma );

	void Multiply( double f );

	void AddLayer( int ix, int iy, int iz, double val );

	void AddGaussian( double h, double px, double py, double pz, double s );

	int m_iAtomCount;
	int m_iRes[3];

	std::vector<CAtomPos*> m_oaAtoms;

	double *m_pBin;

	double m_fCenter[3];

	double m_fStrideA[3];
	double m_fStrideB[3];
	double m_fStrideC[3];
};




CCube *g_pCubeWork;
CCube *g_pVectorField[3];
CCube *g_pCurl[3];
CCube *g_pDivergence;
int g_iCutNormal;
int g_iCutSliceFrom;
int g_iCutSliceTo;
bool g_bShowAtoms;
double g_fAtomRadius;
double g_fVectorScale;
double g_fVectorHeadSize;
double g_fValueRangeMin, g_fValueRangeMax;
double g_fVectorRange;
bool g_bCustomValueRange;
int g_iVectorStride;
int g_iPlaneX1, g_iPlaneY1, g_iPlaneX2, g_iPlaneY2;
int g_iGPInterpolation;
int g_iContourLines;
double g_fAspectRatio;
bool g_bShowTicks;
bool g_bDrawMesh;
int g_iPlotPixel;
double g_fPlotExpo;





void Convolution3D( bool fft, bool debug, CCube *pin, CCube *pker, CCube *pout ) {

	CCube *tempin, *tempker, *tempout;
	int ix, iy, iz, ix2, iy2, iz2, dx, dy, dz;
	int tresx, tresy, tresz, tresxy, kerrad, kshx, kshy, kshz;
	double tf;


	if (debug)
		mprintf("*** Convolution3D ***\n");

	kerrad = (MAX3( pker->m_iRes[0], pker->m_iRes[1], pker->m_iRes[2] ) - 1) / 2;
	tresx = pin->m_iRes[0] + kerrad;
	tresy = pin->m_iRes[1] + kerrad;
	tresz = pin->m_iRes[2] + kerrad;
	tresxy = tresx * tresy;

	if (debug) {
		mprintf("    kerradius = %d\n",kerrad);
		mprintf("    tresx = %d, tresy = %d, tresz = %d\n",tresx,tresy,tresz);

		mprintf("    Writing input CUBE to \"conv_in.cube\"...\n");
		pin->WriteCube("conv_in.cube");

		mprintf("    Writing kernel CUBE to \"conv_ker.cube\"...\n");
		pker->WriteCube("conv_ker.cube");
	}

	tempin = new CCube();
	tempin->m_fStrideA[0] = pin->m_fStrideA[0];
	tempin->m_fStrideB[1] = pin->m_fStrideB[1];
	tempin->m_fStrideC[2] = pin->m_fStrideC[2];
	tempin->m_fCenter[0] = -tresx/2*pin->m_fStrideA[0];
	tempin->m_fCenter[1] = -tresy/2*pin->m_fStrideB[1];
	tempin->m_fCenter[2] = -tresz/2*pin->m_fStrideC[2];
	tempin->m_iRes[0] = tresx;
	tempin->m_iRes[1] = tresy;
	tempin->m_iRes[2] = tresz;
	tempin->m_pBin = new double[tresx*tresy*tresz];
	for (ix=0;ix<tresx*tresy*tresz;ix++)
		tempin->m_pBin[ix] = 0;

	tempker = new CCube();
	tempker->m_fStrideA[0] = pin->m_fStrideA[0];
	tempker->m_fStrideB[1] = pin->m_fStrideB[1];
	tempker->m_fStrideC[2] = pin->m_fStrideC[2];
	tempker->m_fCenter[0] = -tresx/2*pin->m_fStrideA[0];
	tempker->m_fCenter[1] = -tresy/2*pin->m_fStrideB[1];
	tempker->m_fCenter[2] = -tresz/2*pin->m_fStrideC[2];
	tempker->m_iRes[0] = tresx;
	tempker->m_iRes[1] = tresy;
	tempker->m_iRes[2] = tresz;
	tempker->m_pBin = new double[tresx*tresy*tresz];
	for (ix=0;ix<tresx*tresy*tresz;ix++)
		tempker->m_pBin[ix] = 0;

	tempout = new CCube();
	tempout->m_fStrideA[0] = pin->m_fStrideA[0];
	tempout->m_fStrideB[1] = pin->m_fStrideB[1];
	tempout->m_fStrideC[2] = pin->m_fStrideC[2];
	tempout->m_fCenter[0] = -tresx/2*pin->m_fStrideA[0];
	tempout->m_fCenter[1] = -tresy/2*pin->m_fStrideB[1];
	tempout->m_fCenter[2] = -tresz/2*pin->m_fStrideC[2];
	tempout->m_iRes[0] = tresx;
	tempout->m_iRes[1] = tresy;
	tempout->m_iRes[2] = tresz;
	tempout->m_pBin = new double[tresx*tresy*tresz];
	for (ix=0;ix<tresx*tresy*tresz;ix++)
		tempout->m_pBin[ix] = 0;

	for (iz=0;iz<pin->m_iRes[2];iz++)
		for (iy=0;iy<pin->m_iRes[1];iy++)
			for (ix=0;ix<pin->m_iRes[0];ix++)
				tempin->m_pBin[(iz+kerrad/2)*tresy*tresx+(iy+kerrad/2)*tresx+(ix+kerrad/2)] = pin->m_pBin[iz*pin->m_iRes[1]*pin->m_iRes[0]+iy*pin->m_iRes[0]+ix];

	kshx = (tresx/2) - (pker->m_iRes[0]-1)/2;
	kshy = (tresy/2) - (pker->m_iRes[1]-1)/2;
	kshz = (tresz/2) - (pker->m_iRes[2]-1)/2;
	for (iz=0;iz<pker->m_iRes[2];iz++)
		for (iy=0;iy<pker->m_iRes[1];iy++)
			for (ix=0;ix<pker->m_iRes[0];ix++)
				tempker->m_pBin[(iz+kshz)*tresy*tresx+(iy+kshy)*tresx+(ix+kshx)] = pker->m_pBin[iz*pker->m_iRes[1]*pker->m_iRes[0]+iy*pker->m_iRes[0]+ix];

	if (debug) {
		mprintf("    Writing temp input CUBE to \"conv_temp_in.cube\"...\n");
		tempin->WriteCube("conv_temp_in.cube");

		mprintf("    Writing temp kernel CUBE to \"conv_temp_ker.cube\"...\n");
		tempker->WriteCube("conv_temp_ker.cube");
	}

	if (fft) {

		#ifdef USE_CUBETOOL_FFTW

			double tfr, tfi, fac;

			if (debug) {
				mprintf("    FFT-based convolution...\n");
				mprintf("      Preparation...\n");
			}

			fftw_complex *pFTBufWork1;
			fftw_complex *pFTBufWork2;
			fftw_complex *pFTBufInOut;

			fftw_plan plComplexForward1;
			fftw_plan plComplexForward2;
			fftw_plan plComplexBackward;

			pFTBufWork1 = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * tresx * tresy * tresz );
			pFTBufWork2 = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * tresx * tresy * tresz );
			pFTBufInOut = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * tresx * tresy * tresz );

			plComplexForward1 = fftw_plan_dft_3d( tresx, tresy, tresz, pFTBufInOut, pFTBufWork1, FFTW_FORWARD, FFTW_ESTIMATE );

			plComplexForward2 = fftw_plan_dft_3d( tresx, tresy, tresz, pFTBufInOut, pFTBufWork2, FFTW_FORWARD, FFTW_ESTIMATE );

			plComplexBackward = fftw_plan_dft_3d( tresx, tresy, tresz, pFTBufWork1, pFTBufInOut, FFTW_BACKWARD, FFTW_ESTIMATE );

			if (debug)
				mprintf("      Packing 1...\n");

			for (iz=0;iz<tresz;iz++) {
				for (iy=0;iy<tresy;iy++) {
					for (ix=0;ix<tresx;ix++) {
						pFTBufInOut[ix*tresy*tresz+iy*tresz+iz][0] = tempin->m_pBin[iz*tresx*tresy+iy*tresx+ix];
						pFTBufInOut[ix*tresy*tresz+iy*tresz+iz][1] = 0;
					}
				}
			}

			if (debug)
				mprintf("      Forward Transform 1...\n");

			fftw_execute( plComplexForward1 );

			if (debug)
				mprintf("      Packing 2...\n");

			for (iz=0;iz<tresz;iz++) {
				for (iy=0;iy<tresy;iy++) {
					for (ix=0;ix<tresx;ix++) {
						pFTBufInOut[ix*tresy*tresz+iy*tresz+iz][0] = tempker->m_pBin[iz*tresx*tresy+iy*tresx+ix];
						pFTBufInOut[ix*tresy*tresz+iy*tresz+iz][1] = 0;
					}
				}
			}

			if (debug)
				mprintf("      Forward Transform 2...\n");

			fftw_execute( plComplexForward2 );

			if (debug)
				mprintf("      Complex Multiplication...\n");

			for (ix=0;ix<tresx*tresy*tresz;ix++) {
				tfr = pFTBufWork1[ix][0]*pFTBufWork2[ix][0] - pFTBufWork1[ix][1]*pFTBufWork2[ix][1];
				tfi = pFTBufWork1[ix][0]*pFTBufWork2[ix][1] + pFTBufWork1[ix][1]*pFTBufWork2[ix][0];
				pFTBufWork1[ix][0] = tfr;
				pFTBufWork1[ix][1] = tfi;
			}

			if (debug)
				mprintf("      Backward Transform...\n");

			fftw_execute( plComplexBackward );

			if (debug)
				mprintf("      Unacking...\n");

			fac = 1.0/((double)tresx*tresy*tresz);
			for (iz=0;iz<tresz;iz++) {
				iz2 = iz + tresz/2;
				if (iz2 >= tresz)
					iz2 -= tresz;
				for (iy=0;iy<tresy;iy++) {
					iy2 = iy + tresy/2;
					if (iy2 >= tresy)
						iy2 -= tresy;
					for (ix=0;ix<tresx;ix++) {
						ix2 = ix + tresx/2;
						if (ix2 >= tresx)
							ix2 -= tresx;
						tempout->m_pBin[iz*tresx*tresy+iy*tresx+ix] = pFTBufInOut[ix2*tresy*tresz+iy2*tresz+iz2][0] * fac;
					}
				}
			}

			if (debug)
				mprintf("      Deallocation...\n");

			fftw_destroy_plan( plComplexForward1 );
			fftw_destroy_plan( plComplexForward2 );
			fftw_destroy_plan( plComplexBackward );

			fftw_free( pFTBufWork1 );
			fftw_free( pFTBufWork2 );
			fftw_free( pFTBufInOut );

			if (debug)
				mprintf("    Done.\n");

		#else

			eprintf("Error: Cubetool not compiled with USE_CUBETOOL_FFTW.\n");
			abort();

		#endif

	} else {

		if (debug) {
			mprintf("    Loop-based convolution...\n");
			mprintf(WHITE,"      [");
		}

		for (iz=0;iz<pin->m_iRes[2];iz++) {
			for (iy=0;iy<pin->m_iRes[1];iy++) {
				for (ix=0;ix<pin->m_iRes[0];ix++) {

					tf = 0;

					for (iz2=-kerrad;iz2<=kerrad;iz2++) {
						dz = iz+kerrad-iz2;
						for (iy2=-kerrad;iy2<=kerrad;iy2++) {
							dy = iy+kerrad-iy2;
							for (ix2=-kerrad;ix2<=kerrad;ix2++) {
								dx = ix+kerrad-ix2;
								tf += tempin->m_pBin[dz*tresxy+dy*tresx+dx] * tempker->m_pBin[(iz2+tresz/2)*tresxy+(iy2+tresy/2)*tresx+(ix2+tresx/2)];
							}
						}
					}

					tempout->m_pBin[(iz+kerrad)*tresy*tresx+(iy+kerrad)*tresx+(ix+kerrad)] = tf;
				}

				if (debug)
					if (fmod(iz*pin->m_iRes[1]+iy,((double)pin->m_iRes[1]*pin->m_iRes[2])/60.0) < 1.0)
						mprintf(WHITE,"#");
			}
		}

		if (debug)
			mprintf(WHITE,"] Done.\n");
	}

	if (debug) {
		mprintf("    Writing temp output CUBE to \"conv_temp_out.cube\"...\n");
		tempout->WriteCube("conv_temp_out.cube");
	}

	for (iz=0;iz<pin->m_iRes[2];iz++)
		for (iy=0;iy<pin->m_iRes[1];iy++)
			for (ix=0;ix<pin->m_iRes[0];ix++)
				pout->m_pBin[iz*pin->m_iRes[1]*pin->m_iRes[0]+iy*pin->m_iRes[0]+ix] = tempout->m_pBin[(iz+kerrad/2)*tresy*tresx+(iy+kerrad/2)*tresx+(ix+kerrad/2)];

	if (debug) {
		mprintf("    Writing output CUBE to \"conv_out.cube\"...\n");
		pout->WriteCube("conv_out.cube");
	}

	delete tempin;
	delete tempker;
	delete tempout;

	if (debug)
		mprintf("*** Convolution3D Finished ***\n");
}



bool CCube::AddCube( CCube *tcu ) {

	int z;

	if ((tcu->m_iRes[0] != m_iRes[0]) || (tcu->m_iRes[1] != m_iRes[1]) || (tcu->m_iRes[2] != m_iRes[2])) {
		eprintf("CCube::AddCube(): Error: Resolution mismatch (%d x %d x %d vs. %d x %d x %d).\n",m_iRes[0],m_iRes[1],m_iRes[2],tcu->m_iRes[0],tcu->m_iRes[1],tcu->m_iRes[2]);
		return false;
	}

	for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
		m_pBin[z] += tcu->m_pBin[z];

	return true;
}



bool CCube::SubtractCube( CCube *tcu ) {

	int z;

	if ((tcu->m_iRes[0] != m_iRes[0]) || (tcu->m_iRes[1] != m_iRes[1]) || (tcu->m_iRes[2] != m_iRes[2])) {
		eprintf("CCube::SubtractCube(): Error: Resolution mismatch (%d x %d x %d vs. %d x %d x %d).\n",m_iRes[0],m_iRes[1],m_iRes[2],tcu->m_iRes[0],tcu->m_iRes[1],tcu->m_iRes[2]);
		return false;
	}

	for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
		m_pBin[z] -= tcu->m_pBin[z];

	return true;
}



void CCube::Add( double f ) {

	int z;

	for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
		m_pBin[z] += f;
}



void CCube::LimitRange( double mi, double ma ) {

	int z;

	for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++) {
		if (m_pBin[z] < mi)
			m_pBin[z] = mi;
		if (m_pBin[z] > ma)
			m_pBin[z] = ma;
	}
}



void CCube::Multiply( double f ) {

	int z;

	for (z=0;z<m_iRes[0]*m_iRes[1]*m_iRes[2];z++)
		m_pBin[z] *= f;
}



void CCube::AddLayer( int ix, int iy, int iz, double val ) {

	int nrx, nry, nrz, zx, zy, zz;
	double *tbin;

	mprintf("        Adding layer to working buffer...\n");

	nrx = m_iRes[0]+2*ix;
	nry = m_iRes[1]+2*iy;
	nrz = m_iRes[2]+2*iz;

	tbin = new double[nrx*nry*nrz];

	for (zz=0;zz<nrx*nry*nrz;zz++)
		tbin[zz] = val;

	for (zz=0;zz<m_iRes[2];zz++)
		for (zy=0;zy<m_iRes[1];zy++)
			for (zx=0;zx<m_iRes[0];zx++)
				tbin[(zz+iz)*nry*nrx+(zy+iy)*nrx+zx+ix] = m_pBin[zz*m_iRes[1]*m_iRes[0]+zy*m_iRes[0]+zx];

	m_iRes[0] = nrx;
	m_iRes[1] = nry;
	m_iRes[2] = nrz;
	delete[] m_pBin;
	m_pBin = tbin;
}



void CCube::AddGaussian( double h, double px, double py, double pz, double s ) {

	int ix, iy, iz;
	unsigned long idx;
	double tf, xsq, ysq, zsq, ssq;


	mprintf("        Adding Gaussian to cube...\n");
	mprintf(WHITE,"          [");

	ssq = -1.0 / (s * s);
	idx = 0;
	for (iz=0;iz<m_iRes[2];iz++) {

		tf = (m_fCenter[2] + iz*m_fStrideC[2]) - pz;
		if (tf > m_iRes[2]*m_fStrideC[2]/2.0)
			tf -= m_iRes[2]*m_fStrideC[2];
		if (tf < -m_iRes[2]*m_fStrideC[2]/2.0)
			tf += m_iRes[2]*m_fStrideC[2];
		zsq = SQR(tf);

		for (iy=0;iy<m_iRes[1];iy++) {

			tf = (m_fCenter[1] + iy*m_fStrideB[1]) - py;
			if (tf > m_iRes[1]*m_fStrideB[1]/2.0)
				tf -= m_iRes[1]*m_fStrideB[1];
			if (tf < -m_iRes[1]*m_fStrideB[1]/2.0)
				tf += m_iRes[1]*m_fStrideB[1];
			ysq = SQR(tf);

			for (ix=0;ix<m_iRes[0];ix++) {

				tf = (m_fCenter[0] + ix*m_fStrideA[0]) - px;
				if (tf > m_iRes[0]*m_fStrideA[0]/2.0)
					tf -= m_iRes[0]*m_fStrideA[0];
				if (tf < -m_iRes[0]*m_fStrideA[0]/2.0)
					tf += m_iRes[0]*m_fStrideA[0];
				xsq = SQR(tf);

				m_pBin[idx] += h * exp( (xsq + ysq + zsq) * ssq );

				idx++;
			}

			if (fmod(iz*m_iRes[1]+iy,(m_iRes[1]*m_iRes[2])/60.0) < 1.0)
				mprintf(WHITE,"#");
		}
	}

	mprintf(WHITE,"]\n");
}



bool CCube::ImportBinary(const char *s) {

	FILE *a;
	int ix, iy, iz;
	double tf, tmi, tma;


	a = fopen(s,"rb");
	if (a == NULL) {
		eprintf("Error: Could not open \"%s\" for reading.\n",s);
		return false;
	}

	mprintf("        Importing binary data...\n");

	mprintf(WHITE,"          [");

	tmi = 1.0e300;
	tma = -1.0e300;
	for (ix=0;ix<m_iRes[0];ix++) {
		for (iy=0;iy<m_iRes[1];iy++) {
			for (iz=0;iz<m_iRes[2];iz++) {

				(void)!fread(&tf,8,1,a);
				if (feof(a)) {
					eprintf("Error: Unexpected end of file at X=%d/%d, Y=%d/%d, Z=%d/%d.\n",ix+1,m_iRes[0],iy+1,m_iRes[1],iz+1,m_iRes[2]);
					fclose(a);
					return false;
				}

				if (tf > tma)
					tma = tf;
				if (tf < tmi)
					tmi = tf;

				m_pBin[iz*m_iRes[0]*m_iRes[1]+iy*m_iRes[0]+ix] = tf;
	//			if ((ix < 4) && (iy < 4) && (iz < 4))
	//				mprintf("( %3d | %3d | %3d ): %.10e\n",ix+1,iy+1,iz+1,tf);
			}

			if (fmod(iy*m_iRes[0]+ix,(m_iRes[0]*m_iRes[1])/60.0) < 1.0)
				mprintf(WHITE,"#");
		}
	}

	mprintf(WHITE,"]\n");

	mprintf("       Read %s from file.\n",FormatBytes((double)m_iRes[0]*m_iRes[1]*m_iRes[2]*8) );
	mprintf("       Min value: %.10e\n",tmi);
	mprintf("       Max value: %.10e\n",tma);

	fclose(a);

	return true;
}



bool CCube::ImportBinary4D(const char *s, int totalframes, int frame) {

	FILE *a;
	int ix, iy, iz;
	double tf, tmi, tma, *tbuf;


	a = fopen(s,"rb");
	if (a == NULL) {
		eprintf("Error: Could not open \"%s\" for reading.\n",s);
		return false;
	}

	mprintf("        Importing binary data...\n");

	mprintf(WHITE,"          [");

	tbuf = new double[totalframes];

	tmi = 1.0e300;
	tma = -1.0e300;
	for (ix=0;ix<m_iRes[0];ix++) {
		for (iy=0;iy<m_iRes[1];iy++) {
			for (iz=0;iz<m_iRes[2];iz++) {

				(void)!fread(tbuf,8*totalframes,1,a);
				if (feof(a)) {
					eprintf("Error: Unexpected end of file at X=%d/%d, Y=%d/%d, Z=%d/%d.\n",ix+1,m_iRes[0],iy+1,m_iRes[1],iz+1,m_iRes[2]);
					fclose(a);
					return false;
				}

				tf = tbuf[frame];

				if (tf > tma)
					tma = tf;
				if (tf < tmi)
					tmi = tf;

				m_pBin[iz*m_iRes[0]*m_iRes[1]+iy*m_iRes[0]+ix] = tf;
	//			if ((ix < 4) && (iy < 4) && (iz < 4))
	//				mprintf("( %3d | %3d | %3d ): %.10e\n",ix+1,iy+1,iz+1,tf);
			}

			if (fmod(ix*m_iRes[1]+iy,(m_iRes[0]*m_iRes[1])/60.0) < 1.0)
				mprintf(WHITE,"#");
		}
	}

	mprintf(WHITE,"]\n");
	delete[] tbuf;

	mprintf("       Read %s from file.\n",FormatBytes((double)m_iRes[0]*m_iRes[1]*m_iRes[2]*8) );
	mprintf("       Min value: %.10e\n",tmi);
	mprintf("       Max value: %.10e\n",tma);

	fclose(a);

	return true;
}



bool CCube::AnalyzeBinary4D(const char *s, int totalframes) {

	FILE *a;
	int ix, iy, iz, in;
	double *tbuf;
	std::vector<double> vmi, vma, vavg;


	a = fopen(s,"rb");
	if (a == NULL) {
		eprintf("Error: Could not open \"%s\" for reading.\n",s);
		return false;
	}

	mprintf("        Analyzing binary data...\n");

	vmi.resize(totalframes);
	vma.resize(totalframes);
	vavg.resize(totalframes);

	for (ix=0;ix<totalframes;ix++) {
		vmi[ix] = 1.0e300;
		vma[ix] = -1.0e300;
		vavg[ix] = 0;
	}

	mprintf(WHITE,"          [");

	tbuf = new double[totalframes];

	for (ix=0;ix<m_iRes[0];ix++) {
		for (iy=0;iy<m_iRes[1];iy++) {
			for (iz=0;iz<m_iRes[2];iz++) {

				(void)!fread(tbuf,8*totalframes,1,a);
				if (feof(a)) {
					eprintf("Error: Unexpected end of file at X=%d/%d, Y=%d/%d, Z=%d/%d.\n",ix+1,m_iRes[0],iy+1,m_iRes[1],iz+1,m_iRes[2]);
					fclose(a);
					return false;
				}

				for (in=0;in<totalframes;in++) {
					if (tbuf[in] < vmi[in])
						vmi[in] = tbuf[in];
					if (tbuf[in] > vma[in])
						vma[in] = tbuf[in];
					vavg[in] += tbuf[in];
				}
			}

			if (fmod(ix*m_iRes[1]+iy,(m_iRes[0]*m_iRes[1])/60.0) < 1.0)
				mprintf(WHITE,"#");
		}
	}

	mprintf(WHITE,"]\n");
	delete[] tbuf;

	for (ix=0;ix<totalframes;ix++)
		vavg[ix] /= (double)m_iRes[0]*m_iRes[1]*m_iRes[2];

	mprintf("       Read %s from file.\n",FormatBytes((double)m_iRes[0]*m_iRes[1]*m_iRes[2]*8) );

	fclose(a);

	mprintf(WHITE,"       *** Results ***\n");
	for (ix=0;ix<totalframes;ix++)
		mprintf("        Frame %4d:  Min %20.10E,   Max %20.10E,   Avg %20.10E\n",ix+1,vmi[ix],vma[ix],vavg[ix]);
	mprintf("\n");

	return true;
}



bool CCube::ReadCube(const char *s) {

	FILE *a;
	int z, ix, iy, iz;
	char buf[256], *p, *q;
	double tf;
	double mi, ma, sum;
	CAtomPos *ap;


	mprintf(WHITE,"    ------------------------------------------------------\n");
	mprintf("        CCube::ReadCube(): Reading cube file \"%s\"...\n",s);
	a = fopen(s,"rt");
	if (a == NULL) {
		eprintf("CCube::ReadCube(): Error: Could not open \"%s\" for reading.\n",s);
		return false;
	}

	(void)!fgets(buf,256,a);
	p = &buf[strlen(buf)-1];
	while ((*p == '\r') || (*p == '\n'))
		p--;
	p++;
	*p = 0;
	mprintf("        Comment line 1:  \"%s\"\n",buf);
	(void)!fgets(buf,256,a);
	p = &buf[strlen(buf)-1];
	while ((*p == '\r') || (*p == '\n'))
		p--;
	p++;
	*p = 0;
	mprintf("        Comment line 2:  \"%s\"\n",buf);


	(void)!fgets(buf,256,a);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_iAtomCount = atoi(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_fCenter[0] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_fCenter[1] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != '\r') && (*q != '\n') && (*q != 0))
		q++;
	*q = 0;
	m_fCenter[2] = atof(p);


	(void)!fgets(buf,256,a);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_iRes[0] = atoi(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_fStrideA[0] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_fStrideA[1] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != '\r') && (*q != '\n') && (*q != 0))
		q++;
	*q = 0;
	m_fStrideA[2] = atof(p);


	(void)!fgets(buf,256,a);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_iRes[1] = atoi(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_fStrideB[0] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_fStrideB[1] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != '\r') && (*q != '\n') && (*q != 0))
		q++;
	*q = 0;
	m_fStrideB[2] = atof(p);


	(void)!fgets(buf,256,a);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_iRes[2] = atoi(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_fStrideC[0] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while (*q != ' ')
		q++;
	*q = 0;
	m_fStrideC[1] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != '\r') && (*q != '\n') && (*q != 0))
		q++;
	*q = 0;
	m_fStrideC[2] = atof(p);

	for (z=0;z<3;z++) {
		m_fCenter[z]  *= LEN_BOHR2ANGSTROM;
		m_fStrideA[z] *= LEN_BOHR2ANGSTROM;
		m_fStrideB[z] *= LEN_BOHR2ANGSTROM;
		m_fStrideC[z] *= LEN_BOHR2ANGSTROM;
	}


	mprintf("        Resolution:  %d x %d x %d\n",m_iRes[0],m_iRes[1],m_iRes[2]);
	mprintf("        Center:      %10.6f  %10.6f  %10.6f\n",m_fCenter[0],m_fCenter[1],m_fCenter[2]);
	mprintf("        Stride A:    %10.6f  %10.6f  %10.6f\n",m_fStrideA[0],m_fStrideA[1],m_fStrideA[2]);
	mprintf("        Stride B:    %10.6f  %10.6f  %10.6f\n",m_fStrideB[0],m_fStrideB[1],m_fStrideB[2]);
	mprintf("        Stride C:    %10.6f  %10.6f  %10.6f\n",m_fStrideC[0],m_fStrideC[1],m_fStrideC[2]);
	mprintf("        Vector A:    %10.6f  %10.6f  %10.6f\n",m_fStrideA[0]*m_iRes[0],m_fStrideA[1]*m_iRes[0],m_fStrideA[2]*m_iRes[0]);
	mprintf("        Vector B:    %10.6f  %10.6f  %10.6f\n",m_fStrideB[0]*m_iRes[1],m_fStrideB[1]*m_iRes[1],m_fStrideB[2]*m_iRes[1]);
	mprintf("        Vector C:    %10.6f  %10.6f  %10.6f\n",m_fStrideC[0]*m_iRes[2],m_fStrideC[1]*m_iRes[2],m_fStrideC[2]*m_iRes[2]);
	mprintf("        Atom count:  %d\n",m_iAtomCount);

	mprintf("        Reading atom positions...\n");

	for (z=0;z<m_iAtomCount;z++) {

		ap = new CAtomPos();

		(void)!fgets(buf,256,a);
		p = buf;
		while (*p == ' ')
			p++;
		q = p;
		while (*q != ' ')
			q++;
		*q = 0;
		ap->m_iOrd = atoi(p);
		p = q+1;
		while (*p == ' ')
			p++;
		q = p;
		while (*q != ' ')
			q++;
		*q = 0;
		p = q+1;
		while (*p == ' ')
			p++;
		q = p;
		while (*q != ' ')
			q++;
		*q = 0;
		ap->m_fPos[0] = atof(p) * LEN_BOHR2ANGSTROM;
		p = q+1;
		while (*p == ' ')
			p++;
		q = p;
		while (*q != ' ')
			q++;
		*q = 0;
		ap->m_fPos[1] = atof(p) * LEN_BOHR2ANGSTROM;
		p = q+1;
		while (*p == ' ')
			p++;
		q = p;
		while ((*q != ' ') && (*q != '\r') && (*q != '\n') && (*q != 0))
			q++;
		*q = 0;
		ap->m_fPos[2] = atof(p) * LEN_BOHR2ANGSTROM;

		m_oaAtoms.push_back(ap);
	}

	if (feof(a)) {
		eprintf("CCube::ReadCube(): Error: Unexpected end of cube file.\n");
		return false;
	}



	m_pBin = new double[m_iRes[0]*m_iRes[1]*m_iRes[2]];

	mprintf("        Reading volumetric data...\n");
	mprintf(WHITE,"          [");

	ix = 0;
	iy = 0;
	iz = 0;

	mi = 1.0e30;
	ma = -1.0e30;
	sum = 0;

	while (!feof(a)) {
_read:
		(void)!fgets(buf,256,a);
		if (feof(a))
			break;
		p = buf;

_next:
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != '\n') && (*p != 0))
			p++;
		if ((p-q) < 5)
			goto _read;

		*p = 0;

		#ifdef FAST_ATOF
			tf = fast_atof(q);
		#else
			tf = atof(q);
		#endif

		if (tf < mi)
			mi = tf;
		if (tf > ma)
			ma = tf;
		sum += tf;

		m_pBin[iz*m_iRes[0]*m_iRes[1]+iy*m_iRes[0]+ix] = tf;

		iz++;
		if (iz >= m_iRes[2]) {
			iz = 0;
			iy++;
			if (iy >= m_iRes[1]) {
				iy = 0;
				ix++;
				if (fmod(ix,m_iRes[0]/60.0) < 1.0)
					mprintf(WHITE,"#");
			}
		}

		if (ix == m_iRes[0])
			break;

		p++;
		goto _next;
	}

	if (feof(a)) {
		eprintf("CCube::ReadCube(): Error: Unexpected end of cube file stream.\n");
		return false;
	}

	mprintf(WHITE,"]\n");

	fclose(a);

	mprintf("        Finished reading cube file.\n");
	mprintf("        The value range is %.6f .. %.6f\n",mi,ma);
	mprintf("        Sum = %.10E, Integral = %.10E.\n",sum,sum*m_fStrideA[0]*m_fStrideB[1]*m_fStrideC[2]/LEN_BOHR2ANGSTROM/LEN_BOHR2ANGSTROM/LEN_BOHR2ANGSTROM);
	mprintf(WHITE,"    ------------------------------------------------------\n");

	return true;
}



bool CCube::WriteCube(const char *s) {

	FILE *a;
	int z, zi, ix, iy, iz, ti;
	double tfs, tf, tmi, tma;


	a = fopen(s,"wt");
	if (a == NULL) {
		eprintf("CCube::WriteCube(): Error: Could not open \"%s\" for writing.\n",s);
		return false;
	}

	fprintf(a,"Gaussian Cube File\n");
	fprintf(a,"Written by CubeTool, (c) Martin Brehm, 2020 - 2022\n");


	ti = (int)m_oaAtoms.size();
	if (ti == 0)
		ti = 1;

	fprintf(
		a,
		"%10d  %10.6f  %10.6f  %10.6f\n",
		ti,
		m_fCenter[0] / LEN_BOHR2ANGSTROM,
		m_fCenter[1] / LEN_BOHR2ANGSTROM,
		m_fCenter[2] / LEN_BOHR2ANGSTROM
	);

	fprintf(
		a,
		"%10d  %10.6f  %10.6f  %10.6f\n",
		m_iRes[0],
		m_fStrideA[0] / LEN_BOHR2ANGSTROM,
		m_fStrideA[1] / LEN_BOHR2ANGSTROM,
		m_fStrideA[2] / LEN_BOHR2ANGSTROM
	);

	fprintf(
		a,
		"%10d  %10.6f  %10.6f  %10.6f\n",
		m_iRes[1],
		m_fStrideB[0] / LEN_BOHR2ANGSTROM,
		m_fStrideB[1] / LEN_BOHR2ANGSTROM,
		m_fStrideB[2] / LEN_BOHR2ANGSTROM
	);

	fprintf(
		a,
		"%10d  %10.6f  %10.6f  %10.6f\n",
		m_iRes[2],
		m_fStrideC[0] / LEN_BOHR2ANGSTROM,
		m_fStrideC[1] / LEN_BOHR2ANGSTROM,
		m_fStrideC[2] / LEN_BOHR2ANGSTROM
	);


	for (z=0;z<(int)m_oaAtoms.size();z++)
		fprintf(
			a,
			"%10d  0.00  %10.6f  %10.6f  %10.6f\n",
			m_oaAtoms[z]->m_iOrd,
			m_oaAtoms[z]->m_fPos[0] / LEN_BOHR2ANGSTROM,
			m_oaAtoms[z]->m_fPos[1] / LEN_BOHR2ANGSTROM,
			m_oaAtoms[z]->m_fPos[2] / LEN_BOHR2ANGSTROM
		);

	if (m_oaAtoms.size() == 0)
		fprintf(a,"0 0 0 0 0\n");


	mprintf(WHITE,"      [");

	tfs = (m_iRes[0]*m_iRes[1]*m_iRes[2])/60.0;
	z = 0;
	zi = 0;
	tmi = 1.0e300;
	tma = -1.0e300;
	for (ix=0;ix<m_iRes[0];ix++) {
		for (iy=0;iy<m_iRes[1];iy++) {
			for (iz=0;iz<m_iRes[2];iz++) {
				tf = m_pBin[iz*m_iRes[0]*m_iRes[1]+iy*m_iRes[0]+ix];
				// Limit to two-digit exponents
				if ((tf > 0) && (tf < 1.0e-98))
					tf = 1.0e-98;
				if ((tf < 0) && (tf > -1.0e-98))
					tf = -1.0e-98;
				if (tf > 1.0e98)
					tf = 1.0e98;
				if (tf < -1.0e98)
					tf = -1.0e98;
				if (tf > tma)
					tma = tf;
				if (tf < tmi)
					tmi = tf;
				fprintf(
					a,
					"  %11.5E",
					tf
				);
				z++;
				zi++;
				if (zi == 6) {
					fprintf(a,"\n");
					zi = 0;
				}
				if (fmod(z,tfs) < 1.0)
					mprintf(WHITE,"#");
			}
			if (zi != 0) {
				fprintf(a,"\n");
				zi = 0;
			}
		}
	}

	fclose(a);

	mprintf(WHITE,"] Done.\n");
	mprintf("    Min value: %.10e\n",tmi);
	mprintf("    Max value: %.10e\n",tma);

	return true;
}



bool SORT_AtomsX(CAtomPos *a1, CAtomPos *a2) {
	return (a1->m_fPos[0] > a2->m_fPos[0]);
}



bool SORT_AtomsY(CAtomPos *a1, CAtomPos *a2) {
	return (a1->m_fPos[1] > a2->m_fPos[1]);
}



bool SORT_AtomsZ(CAtomPos *a1, CAtomPos *a2) {
	return (a1->m_fPos[2] > a2->m_fPos[2]);
}



bool WriteSlice( CCube *voldat, int normal, int slice1, int slice2, const char *s ) {

	C2DF *df;
	int z, z2, ix, iy, iz;
	double rad, cr, cg, cb, px, py, fac, mi, ma;


	mprintf(WHITE,"    ------------------------------------------------------\n");
	mprintf("      WriteSlice(): Creating cut plane.\n");
	df = new C2DF();

	mi = 1.0e30;
	ma = -1.0e30;

	switch( normal ) {
		case 1: // X
			if (g_iPlaneX1 == -1) {
				g_iPlaneX1 = 0;
				g_iPlaneY1 = 0;
				g_iPlaneX2 = voldat->m_iRes[1] - 1;
				g_iPlaneY2 = voldat->m_iRes[2] - 1;
			}
			df->m_iRes[0] = g_iPlaneX2 - g_iPlaneX1 + 1;
			df->m_iRes[1] = g_iPlaneY2 - g_iPlaneY1 + 1;
			df->m_fMinVal[0] = floor(100.0*voldat->m_fStrideB[1]*(g_iPlaneX1-voldat->m_iRes[1]/2)+0.5);
			df->m_fMaxVal[0] = floor(100.0*voldat->m_fStrideB[1]*(g_iPlaneX2-voldat->m_iRes[1]/2)+0.5);
			df->m_fMinVal[1] = floor(100.0*voldat->m_fStrideC[2]*(g_iPlaneY1-voldat->m_iRes[2]/2)+0.5);
			df->m_fMaxVal[1] = floor(100.0*voldat->m_fStrideC[2]*(g_iPlaneY2-voldat->m_iRes[2]/2)+0.5);
			if (g_bShowTicks) {
				df->SetLabelX("Y Position / pm");
				df->SetLabelY("Z Position / pm");
				df->SetLabelZ("Value");
			}
			df->Create();
			for (ix=slice1;ix<=slice2;ix++)
				for (iy=g_iPlaneX1;iy<=g_iPlaneX2;iy++)
					for (iz=g_iPlaneY1;iz<=g_iPlaneY2;iz++)
						df->m_pBin[(iz-g_iPlaneY1)*df->m_iRes[0]+iy-g_iPlaneX1] += voldat->m_pBin[iz*voldat->m_iRes[0]*voldat->m_iRes[1]+iy*voldat->m_iRes[0]+ix];
			for (z=0;z<df->m_iRes[0]*df->m_iRes[1];z++) {
			 	df->m_pBin[z] /= (double)(slice2-slice1+1);
				if (g_bCustomValueRange) {
					if (df->m_pBin[z] > g_fValueRangeMax)
						df->m_pBin[z] = g_fValueRangeMax;
					if (df->m_pBin[z] < g_fValueRangeMin)
						df->m_pBin[z] = g_fValueRangeMin;
				}
				if (df->m_pBin[z] < mi)
					mi = df->m_pBin[z];
				if (df->m_pBin[z] > ma)
					ma = df->m_pBin[z];
			}
			break;

		case 2: // Y
			if (g_iPlaneX1 == -1) {
				g_iPlaneX1 = 0;
				g_iPlaneY1 = 0;
				g_iPlaneX2 = voldat->m_iRes[0] - 1;
				g_iPlaneY2 = voldat->m_iRes[2] - 1;
			}
			df->m_iRes[0] = g_iPlaneX2 - g_iPlaneX1 + 1;
			df->m_iRes[1] = g_iPlaneY2 - g_iPlaneY1 + 1;
			df->m_fMinVal[0] = floor(100.0*voldat->m_fStrideA[0]*(g_iPlaneX1-voldat->m_iRes[0]/2)+0.5);
			df->m_fMaxVal[0] = floor(100.0*voldat->m_fStrideA[0]*(g_iPlaneX2-voldat->m_iRes[0]/2)+0.5);
			df->m_fMinVal[1] = floor(100.0*voldat->m_fStrideC[2]*(g_iPlaneY1-voldat->m_iRes[2]/2)+0.5);
			df->m_fMaxVal[1] = floor(100.0*voldat->m_fStrideC[2]*(g_iPlaneY2-voldat->m_iRes[2]/2)+0.5);
			if (g_bShowTicks) {
				df->SetLabelX("X Position / pm");
				df->SetLabelY("Z Position / pm");
				df->SetLabelZ("Value");
			}
			df->Create();
			for (iy=slice1;iy<=slice2;iy++)
				for (ix=g_iPlaneX1;ix<=g_iPlaneX2;ix++)
					for (iz=g_iPlaneY1;iz<=g_iPlaneY2;iz++)
						df->m_pBin[(iz-g_iPlaneY1)*df->m_iRes[0]+ix-g_iPlaneX1] += voldat->m_pBin[iz*voldat->m_iRes[0]*voldat->m_iRes[1]+iy*voldat->m_iRes[0]+ix];
			for (z=0;z<df->m_iRes[0]*df->m_iRes[1];z++) {
			 	df->m_pBin[z] /= (double)(slice2-slice1+1);
				if (g_bCustomValueRange) {
					if (df->m_pBin[z] > g_fValueRangeMax)
						df->m_pBin[z] = g_fValueRangeMax;
					if (df->m_pBin[z] < g_fValueRangeMin)
						df->m_pBin[z] = g_fValueRangeMin;
				}
				if (df->m_pBin[z] < mi)
					mi = df->m_pBin[z];
				if (df->m_pBin[z] > ma)
					ma = df->m_pBin[z];
			}
			break;

		case 3: // Z
			if (g_iPlaneX1 == -1) {
				g_iPlaneX1 = 0;
				g_iPlaneY1 = 0;
				g_iPlaneX2 = voldat->m_iRes[0] - 1;
				g_iPlaneY2 = voldat->m_iRes[1] - 1;
			}
			df->m_iRes[0] = g_iPlaneX2 - g_iPlaneX1 + 1;
			df->m_iRes[1] = g_iPlaneY2 - g_iPlaneY1 + 1;
			df->m_fMinVal[0] = floor(100.0*voldat->m_fStrideA[0]*(g_iPlaneX1-voldat->m_iRes[0]/2)+0.5);
			df->m_fMaxVal[0] = floor(100.0*voldat->m_fStrideA[0]*(g_iPlaneX2-(voldat->m_iRes[0]-1)/2)+0.5);
			df->m_fMinVal[1] = floor(100.0*voldat->m_fStrideB[1]*(g_iPlaneY1-voldat->m_iRes[1]/2)+0.5);
			df->m_fMaxVal[1] = floor(100.0*voldat->m_fStrideB[1]*(g_iPlaneY2-(voldat->m_iRes[1]-1)/2)+0.5);
			if (g_bShowTicks) {
				df->SetLabelX("X Position / pm");
				df->SetLabelY("Y Position / pm");
				df->SetLabelZ("Value");
			}
			df->Create();
			for (iz=slice1;iz<=slice2;iz++)
				for (ix=g_iPlaneX1;ix<=g_iPlaneX2;ix++)
					for (iy=g_iPlaneY1;iy<=g_iPlaneY2;iy++)
						df->m_pBin[(iy-g_iPlaneY1)*df->m_iRes[0]+ix-g_iPlaneX1] += voldat->m_pBin[iz*voldat->m_iRes[0]*voldat->m_iRes[1]+iy*voldat->m_iRes[0]+ix];
			for (z=0;z<df->m_iRes[0]*df->m_iRes[1];z++) {
			 	df->m_pBin[z] /= (double)(slice2-slice1+1);
				if (g_bCustomValueRange) {
					if (df->m_pBin[z] > g_fValueRangeMax)
						df->m_pBin[z] = g_fValueRangeMax;
					if (df->m_pBin[z] < g_fValueRangeMin)
						df->m_pBin[z] = g_fValueRangeMin;
				}
				if (df->m_pBin[z] < mi)
					mi = df->m_pBin[z];
				if (df->m_pBin[z] > ma)
					ma = df->m_pBin[z];
			}
			break;

		default:
			eprintf("WriteSlice): Error: Unknown normal setting %d.\n",normal);
			return false;
	}

	mprintf("      The value range in the cut plane is %.6f .. %.6f\n",mi,ma);

	if (mi < 0) {
		mprintf("      There are negative values, switching on plus-minus color scale.\n");
		df->m_iColorScale = 4;
	}


	if (g_bShowAtoms) {

		switch( normal ) {
			case 1: // X
				std::stable_sort( voldat->m_oaAtoms.begin(), voldat->m_oaAtoms.end(), SORT_AtomsX );
				break;
			case 2: // Y
				std::stable_sort( voldat->m_oaAtoms.begin(), voldat->m_oaAtoms.end(), SORT_AtomsY );
				break;
			case 3: // Z
				std::stable_sort( voldat->m_oaAtoms.begin(), voldat->m_oaAtoms.end(), SORT_AtomsZ );
				break;
		}


		for (z=0;z<(int)voldat->m_oaAtoms.size();z++) {

			px = 0;
			py = 0;

			switch( normal ) {
				case 1: // X
					px = (voldat->m_oaAtoms[z]->m_fPos[1] - voldat->m_fCenter[1]) * 100.0 - voldat->m_fStrideB[1]*voldat->m_iRes[1]*50.0;
					py = (voldat->m_oaAtoms[z]->m_fPos[2] - voldat->m_fCenter[2]) * 100.0 - voldat->m_fStrideC[2]*voldat->m_iRes[2]*50.0;
					break;
				case 2: // Y
					px = (voldat->m_oaAtoms[z]->m_fPos[0] - voldat->m_fCenter[0]) * 100.0 - voldat->m_fStrideA[0]*voldat->m_iRes[0]*50.0;
					py = (voldat->m_oaAtoms[z]->m_fPos[2] - voldat->m_fCenter[2]) * 100.0 - voldat->m_fStrideC[2]*voldat->m_iRes[2]*50.0;
					break;
				case 3: // Z
					px = (voldat->m_oaAtoms[z]->m_fPos[0] - voldat->m_fCenter[0]) * 100.0 - voldat->m_fStrideA[0]*voldat->m_iRes[0]*50.0;
					py = (voldat->m_oaAtoms[z]->m_fPos[1] - voldat->m_fCenter[1]) * 100.0 - voldat->m_fStrideB[1]*voldat->m_iRes[1]*50.0;
					break;
			}

			for (z2=0;z2<g_oaElements.GetSize();z2++)
				if (((voldat->m_oaAtoms[z]->m_iOrd == 1) && (strcmp(((CElement*)g_oaElements[z2])->m_sLabel,"H") == 0)) ||
					((voldat->m_oaAtoms[z]->m_iOrd > 1) && (((CElement*)g_oaElements[z2])->m_iOrd == voldat->m_oaAtoms[z]->m_iOrd))) {
					rad = ((CElement*)g_oaElements[z2])->m_fRadius * 0.6;
					cr = ((CElement*)g_oaElements[z2])->m_iColorR / 255.0;
					cg = ((CElement*)g_oaElements[z2])->m_iColorG / 255.0;
					cb = ((CElement*)g_oaElements[z2])->m_iColorB / 255.0;
					goto _elfound;
				}

			rad = 25.0;
			cr = 1.0;
			cg = 0.0;
			cb = 1.0;

_elfound:
			if (((px - rad*g_fAtomRadius) < df->m_fMinVal[0]) || ((px + rad*g_fAtomRadius) > df->m_fMaxVal[0]) || ((py - rad*g_fAtomRadius) < df->m_fMinVal[1]) || ((py + rad*g_fAtomRadius) > df->m_fMaxVal[1]))
				continue;

			df->AddCircle(
				px,
				py,
				rad * g_fAtomRadius,
				cr,
				cg,
				cb
			);
		}
	}


	if (g_pVectorField[0] != NULL) {

		df->m_bVectorField = true;

		mi = -1.0e30;
		ma = -1.0e30;

		switch( normal ) {

			case 1: // X
				df->m_iVectorRes[0] = df->m_iRes[0] / g_iVectorStride;
				df->m_iVectorRes[1] = df->m_iRes[1] / g_iVectorStride;
				df->m_faVectorData.resize(df->m_iVectorRes[0]*df->m_iVectorRes[1]*2);
				for (z=0;z<df->m_iVectorRes[0]*df->m_iVectorRes[1]*2;z++)
					df->m_faVectorData[z] = 0;
				for (ix=slice1;ix<=slice2;ix++) {
					for (iy=0;iy<df->m_iVectorRes[0];iy++) {
						for (iz=0;iz<df->m_iVectorRes[1];iz++) {
							df->m_faVectorData[(iz*df->m_iVectorRes[0]+iy)*2]   += g_pVectorField[1]->m_pBin[(iz*g_iVectorStride+g_iPlaneY1)*voldat->m_iRes[0]*voldat->m_iRes[1]+(iy*g_iVectorStride+g_iPlaneX1)*voldat->m_iRes[0]+ix];
							df->m_faVectorData[(iz*df->m_iVectorRes[0]+iy)*2+1] += g_pVectorField[2]->m_pBin[(iz*g_iVectorStride+g_iPlaneY1)*voldat->m_iRes[0]*voldat->m_iRes[1]+(iy*g_iVectorStride+g_iPlaneX1)*voldat->m_iRes[0]+ix];
						}
					}
				}
				for (z=0;z<df->m_iVectorRes[0]*df->m_iVectorRes[1]*2;z++)
					if (df->m_faVectorData[z] > ma)
						ma = df->m_faVectorData[z];
				if (g_fVectorRange > 0)
					fac = g_fVectorScale * 0.5 / g_fVectorRange;
				else
					fac = g_fVectorScale * 0.5 / ma;
				for (z=0;z<df->m_iVectorRes[0]*df->m_iVectorRes[1]*2;z++) {
					df->m_faVectorData[z] *= fac;
					if (df->m_faVectorData[z] > mi)
						mi = df->m_faVectorData[z];
				}
				break;

			case 2: // Y
				df->m_iVectorRes[0] = df->m_iRes[0] / g_iVectorStride;
				df->m_iVectorRes[1] = df->m_iRes[1] / g_iVectorStride;
				df->m_faVectorData.resize(df->m_iVectorRes[0]*df->m_iVectorRes[1]*2);
				for (z=0;z<df->m_iVectorRes[0]*df->m_iVectorRes[1]*2;z++)
					df->m_faVectorData[z] = 0;
				for (iy=slice1;iy<=slice2;iy++) {
					for (ix=0;ix<df->m_iVectorRes[0];ix++) {
						for (iz=0;iz<df->m_iVectorRes[1];iz++) {
							df->m_faVectorData[(iz*df->m_iVectorRes[0]+ix)*2]   += g_pVectorField[0]->m_pBin[(iz*g_iVectorStride+g_iPlaneY1)*voldat->m_iRes[0]*voldat->m_iRes[1]+iy*voldat->m_iRes[0]+ix*g_iVectorStride+g_iPlaneX1];
							df->m_faVectorData[(iz*df->m_iVectorRes[0]+ix)*2+1] += g_pVectorField[2]->m_pBin[(iz*g_iVectorStride+g_iPlaneY1)*voldat->m_iRes[0]*voldat->m_iRes[1]+iy*voldat->m_iRes[0]+ix*g_iVectorStride+g_iPlaneX1];
						}
					}
				}
				for (z=0;z<df->m_iVectorRes[0]*df->m_iVectorRes[1]*2;z++)
					if (df->m_faVectorData[z] > ma)
						ma = df->m_faVectorData[z];
				if (g_fVectorRange > 0)
					fac = g_fVectorScale * 0.5 / g_fVectorRange;
				else
					fac = g_fVectorScale * 0.5 / ma;
				for (z=0;z<df->m_iVectorRes[0]*df->m_iVectorRes[1]*2;z++) {
					df->m_faVectorData[z] *= fac;
					if (df->m_faVectorData[z] > mi)
						mi = df->m_faVectorData[z];
				}
				break;

			case 3: // Z
				df->m_iVectorRes[0] = df->m_iRes[0] / g_iVectorStride;
				df->m_iVectorRes[1] = df->m_iRes[1] / g_iVectorStride;
				df->m_faVectorData.resize(df->m_iVectorRes[0]*df->m_iVectorRes[1]*2);
				for (z=0;z<df->m_iVectorRes[0]*df->m_iVectorRes[1]*2;z++)
					df->m_faVectorData[z] = 0;
				for (iz=slice1;iz<=slice2;iz++) {
					for (ix=0;ix<df->m_iVectorRes[0];ix++) {
						for (iy=0;iy<df->m_iVectorRes[1];iy++) {
							df->m_faVectorData[(iy*df->m_iVectorRes[0]+ix)*2]   += g_pVectorField[0]->m_pBin[iz*voldat->m_iRes[0]*voldat->m_iRes[1]+(iy*g_iVectorStride+g_iPlaneY1)*voldat->m_iRes[0]+ix*g_iVectorStride+g_iPlaneX1];
							df->m_faVectorData[(iy*df->m_iVectorRes[0]+ix)*2+1] += g_pVectorField[1]->m_pBin[iz*voldat->m_iRes[0]*voldat->m_iRes[1]+(iy*g_iVectorStride+g_iPlaneY1)*voldat->m_iRes[0]+ix*g_iVectorStride+g_iPlaneX1];
						}
					}
				}
				for (z=0;z<df->m_iVectorRes[0]*df->m_iVectorRes[1]*2;z++)
					if (df->m_faVectorData[z] > ma)
						ma = df->m_faVectorData[z];
				if (g_fVectorRange > 0)
					fac = g_fVectorScale * 0.5 / g_fVectorRange;
				else
					fac = g_fVectorScale * 0.5 / ma;
				for (z=0;z<df->m_iVectorRes[0]*df->m_iVectorRes[1]*2;z++) {
					df->m_faVectorData[z] *= fac;
					if (df->m_faVectorData[z] > mi)
						mi = df->m_faVectorData[z];
				}
				break;
		}

		mprintf("      The maximum projected vector length in the cut plane is %.6f\n",ma);
		mprintf("      The resulting maximum vector length after rescaling is %.6f\n",mi);

		df->m_fVectorHeadSize = g_fVectorHeadSize * 20.0;
	}

	if (g_iContourLines > 0) {
		df->m_bContourLines = true;
		df->m_iContourLines = g_iContourLines;
	} else
		df->m_bContourLines = false;

	df->m_iGPInterpolation = g_iGPInterpolation;
	df->m_bDrawMesh = g_bDrawMesh;
	df->m_bShowXTicks = g_bShowTicks;
	df->m_bShowYTicks = g_bShowTicks;
	df->m_iPlotPixel = g_iPlotPixel;
	df->m_fAspectRatio = g_fAspectRatio;
	df->m_fPlotExp = g_fPlotExpo;

	if (g_bCustomValueRange) {
		df->m_fCustomValueMin = g_fValueRangeMin;
		df->m_fCustomValueMax = g_fValueRangeMax;
	}

	mprintf("      Writing Gnuplot input data...\n");
	df->WriteGnuplotInput("",s,"",g_bCustomValueRange);

	mprintf("      Writing Mathematica notebook...\n");
	df->WriteMathematicaNb("",s,".nb",g_bCustomValueRange);

	delete df;

	mprintf("      WriteSlice(): Creating cut plane done.\n");
	mprintf(WHITE,"    ------------------------------------------------------\n");

	return true;
}



void ComputeCurl( CCube *pinx, CCube *piny, CCube *pinz, CCube *poutx, CCube *pouty, CCube *poutz ) {

	int z, ix, iy, iz, idx;
	int resx, resy, resz, resxy;


	resx = pinx->m_iRes[0];
	resy = pinx->m_iRes[1];
	resz = pinx->m_iRes[2];
	resxy = resx*resy;

	for (z=0;z<resx*resy*resz;z++) {
		poutx->m_pBin[z] = 0;
		pouty->m_pBin[z] = 0;
		poutz->m_pBin[z] = 0;
	}

	mprintf(WHITE,"        [");

	for (iz=1;iz<resz-1;iz++) {
		for (iy=1;iy<resy-1;iy++) {
			for (ix=1;ix<resx-1;ix++) {

				idx = ix + iy*resx + iz*resxy;

				poutx->m_pBin[idx] =  (pinz->m_pBin[ix+(iy+1)*resx+iz*resxy] - pinz->m_pBin[ix+(iy-1)*resx+iz*resxy]) / 2.0
				                     -(piny->m_pBin[ix+iy*resx+(iz+1)*resxy] - piny->m_pBin[ix+iy*resx+(iz-1)*resxy]) / 2.0;

				pouty->m_pBin[idx] =  (pinx->m_pBin[ix+iy*resx+(iz+1)*resxy] - pinx->m_pBin[ix+iy*resx+(iz-1)*resxy]) / 2.0
				                     -(pinz->m_pBin[(ix+1)+iy*resx+iz*resxy] - pinz->m_pBin[(ix-1)+iy*resx+iz*resxy]) / 2.0;

				poutz->m_pBin[idx] =  (piny->m_pBin[(ix+1)+iy*resx+iz*resxy] - piny->m_pBin[(ix-1)+iy*resx+iz*resxy]) / 2.0
				                     -(pinx->m_pBin[ix+(iy+1)*resx+iz*resxy] - pinx->m_pBin[ix+(iy-1)*resx+iz*resxy]) / 2.0;
			}
			if (fmod(iz*resy+iy,(double)resy*resz/60.0) < 1.0)
				mprintf(WHITE,"#");
		}
	}
	mprintf(WHITE,"] Done.\n");
}



void ComputeCurl() {

	int z;


	for (z=0;z<3;z++) {
		if (g_pCurl[z] != NULL)
			delete g_pCurl[z];
		g_pCurl[z] = new CCube();
		g_pCurl[z]->CopyFrom(g_pVectorField[z]);
	}

	ComputeCurl( g_pVectorField[0], g_pVectorField[1], g_pVectorField[2], g_pCurl[0], g_pCurl[1], g_pCurl[2] );

	mprintf("      Writing curl to \"curl_x.cube\", \"curl_y.cube\", and \"curl_z.cube\"...\n");
	g_pCurl[0]->WriteCube( "curl_x.cube" );
	g_pCurl[1]->WriteCube( "curl_y.cube" );
	g_pCurl[2]->WriteCube( "curl_z.cube" );
}



void ComputeDivergence( CCube *pinx, CCube *piny, CCube *pinz, CCube *pout ) {

	int z, ix, iy, iz, idx;
	int resx, resy, resz, resxy;


	resx = pinx->m_iRes[0];
	resy = pinx->m_iRes[1];
	resz = pinx->m_iRes[2];
	resxy = resx*resy;

	for (z=0;z<resx*resy*resz;z++)
		pout->m_pBin[z] = 0;

	mprintf(WHITE,"        [");

	for (iz=1;iz<resz-1;iz++) {
		for (iy=1;iy<resy-1;iy++) {
			for (ix=1;ix<resx-1;ix++) {

				idx = ix + iy*resx + iz*resxy;

				pout->m_pBin[idx] =   (pinx->m_pBin[(ix+1)+iy*resx+iz*resxy] - pinx->m_pBin[(ix-1)+iy*resx+iz*resxy]) / 2.0
				                    + (piny->m_pBin[ix+(iy+1)*resx+iz*resxy] - piny->m_pBin[ix+(iy-1)*resx+iz*resxy]) / 2.0
				                    + (pinz->m_pBin[ix+iy*resx+(iz+1)*resxy] - pinz->m_pBin[ix+iy*resx+(iz-1)*resxy]) / 2.0;

			}
			if (fmod(iz*resy+iy,(double)resy*resz/60.0) < 1.0)
				mprintf(WHITE,"#");
		}
	}
	mprintf(WHITE,"] Done.\n");
}



void ComputeDivergence() {

	if (g_pDivergence != NULL)
		delete g_pDivergence;
	g_pDivergence = new CCube();
	g_pDivergence->CopyFrom(g_pVectorField[0]);

	ComputeDivergence( g_pVectorField[0], g_pVectorField[1], g_pVectorField[2], g_pDivergence );

	mprintf("      Writing divergence to \"div.cube\"...\n");
	g_pDivergence->WriteCube( "div.cube" );
}



void ComputeGradient( CCube *pin, CCube *poutx, CCube *pouty, CCube *poutz ) {

	int z, ix, iy, iz, idx;
	int resx, resy, resz, resxy;


	resx = pin->m_iRes[0];
	resy = pin->m_iRes[1];
	resz = pin->m_iRes[2];
	resxy = resx*resy;

	for (z=0;z<resx*resy*resz;z++) {
		poutx->m_pBin[z] = 0;
		pouty->m_pBin[z] = 0;
		poutz->m_pBin[z] = 0;
	}

	mprintf(WHITE,"        [");

	for (iz=1;iz<resz-1;iz++) {
		for (iy=1;iy<resy-1;iy++) {
			for (ix=1;ix<resx-1;ix++) {

				idx = ix + iy*resx + iz*resxy;

				poutx->m_pBin[idx] = (pin->m_pBin[(ix+1)+iy*resx+iz*resxy] - pin->m_pBin[(ix-1)+iy*resx+iz*resxy]) / 2.0;
				pouty->m_pBin[idx] = (pin->m_pBin[ix+(iy+1)*resx+iz*resxy] - pin->m_pBin[ix+(iy-1)*resx+iz*resxy]) / 2.0;
				poutz->m_pBin[idx] = (pin->m_pBin[ix+iy*resx+(iz+1)*resxy] - pin->m_pBin[ix+iy*resx+(iz-1)*resxy]) / 2.0;
			}
			if (fmod(iz*resy+iy,(double)resy*resz/60.0) < 1.0)
				mprintf(WHITE,"#");
		}
	}
	mprintf(WHITE,"] Done.\n");
}



void CreateAbs() {

	int z;

	if (g_pCubeWork != NULL)
		delete g_pCubeWork;

	g_pCubeWork = new CCube();
	g_pCubeWork->CopyFrom( g_pVectorField[0] );

	for (z=0;z<g_pCubeWork->m_iRes[0]*g_pCubeWork->m_iRes[1]*g_pCubeWork->m_iRes[2];z++)
		g_pCubeWork->m_pBin[z] = sqrt( SQR(g_pVectorField[0]->m_pBin[z]) + SQR(g_pVectorField[1]->m_pBin[z]) + SQR(g_pVectorField[2]->m_pBin[z]) );
}



void CropCube( CCube *cufrom, CCube *cuto, int x1, int y1, int z1, int x2, int y2, int z2 ) {

	int ix, iy, iz;


	if ((x1 >= cufrom->m_iRes[0]) || (x2 >= cufrom->m_iRes[0]) || (x1 > x2)) {
		eprintf("CropCube(): Error in X range. 0 <= %d <= %d < %d is not fulfilled.\n",x1,x2,cufrom->m_iRes[0]);
		abort();
	}

	if ((y1 >= cufrom->m_iRes[1]) || (y2 >= cufrom->m_iRes[1]) || (y1 > y2)) {
		eprintf("CropCube(): Error in Y range. 0 <= %d <= %d < %d is not fulfilled.\n",y1,y2,cufrom->m_iRes[1]);
		abort();
	}

	if ((z1 >= cufrom->m_iRes[2]) || (z2 >= cufrom->m_iRes[2]) || (z1 > z2)) {
		eprintf("CropCube(): Error in Z range. 0 <= %d <= %d < %d is not fulfilled.\n",z1,z2,cufrom->m_iRes[2]);
		abort();
	}

	cuto->CopyFrom( cufrom );

	cuto->m_iRes[0] = x2-x1+1;
	cuto->m_iRes[1] = y2-y1+1;
	cuto->m_iRes[2] = z2-z1+1;

	if (cuto->m_pBin != NULL)
		delete[] cuto->m_pBin;

	cuto->m_pBin = new double[cuto->m_iRes[0]*cuto->m_iRes[1]*cuto->m_iRes[2]];

	mprintf(WHITE,"      [");
	for (iz=0;iz<cuto->m_iRes[2];iz++) {
		for (iy=0;iy<cuto->m_iRes[1];iy++) {

			for (ix=0;ix<cuto->m_iRes[0];ix++)
				cuto->m_pBin[iz*cuto->m_iRes[0]*cuto->m_iRes[1]+iy*cuto->m_iRes[0]+ix] = cufrom->m_pBin[(iz+z1)*cufrom->m_iRes[0]*cufrom->m_iRes[1]+(iy+y1)*cufrom->m_iRes[0]+ix+x1];

			if (fmod(iz*cuto->m_iRes[1]+iy,(double)cuto->m_iRes[1]*cuto->m_iRes[2]/60.0) < 1.0)
				mprintf(WHITE,"#");
		}
	}
	mprintf(WHITE,"] Done.\n");
}



void EnvelopeCosSq( CCube *pcu, double inner ) {

	int ix, iy, iz;
	double tf, fx, fy, fz;


	mprintf(WHITE,"      [");

	for (iz=0;iz<pcu->m_iRes[2];iz++) {

		tf = (double)iz/(pcu->m_iRes[2]-1)*2.0;
		if (tf > 1.0)
			tf = 2.0 - tf;
		tf /= 1.0 - inner;
		if (tf < 1.0)
			fz = pow2( sin( tf*Pi/2.0 ) );
		else
			fz = 1.0;

		for (iy=0;iy<pcu->m_iRes[1];iy++) {

			tf = (double)iy/(pcu->m_iRes[1]-1)*2.0;
			if (tf > 1.0)
				tf = 2.0 - tf;
			tf /= 1.0 - inner;
			if (tf < 1.0)
				fy = pow2( sin( tf*Pi/2.0 ) );
			else
				fy = 1.0;

			for (ix=0;ix<pcu->m_iRes[0];ix++) {

				tf = (double)ix/(pcu->m_iRes[0]-1)*2.0;
				if (tf > 1.0)
					tf = 2.0 - tf;
				tf /= 1.0 - inner;
				if (tf < 1.0)
					fx = pow2( sin( tf*Pi/2.0 ) );
				else
					fx = 1.0;

				pcu->m_pBin[iz*pcu->m_iRes[0]*pcu->m_iRes[1]+iy*pcu->m_iRes[0]+ix] *= fx*fy*fz;
			}

			if (fmod(iz*pcu->m_iRes[1]+iy,(double)pcu->m_iRes[1]*pcu->m_iRes[2]/60.0) < 1.0)
				mprintf(WHITE,"#");
		}
	}

	mprintf(WHITE,"] Done.\n");
}



void Helmholtz_Decomp(int kerrad, double eps) {

	// See  https://en.wikipedia.org/wiki/Helmholtz_decomposition

	int ix, iy, iz;
	double rxs, rys, rzs, tama, tdma;
	CCube *pcPot, *pcAX, *pcAY, *pcAZ, *pcKer, *pcRAX, *pcRAY, *pcRAZ, *pcGPotX, *pcGPotY, *pcGPotZ, *pcOutX, *pcOutY, *pcOutZ;


	mprintf("*** Computing Helmholtz decomposition of vector field (experimental) ***\n\n");

	mprintf("    Grid %d x %d x %d, Kernel radius %d\n",g_pVectorField[0]->m_iRes[0],g_pVectorField[0]->m_iRes[1],g_pVectorField[2]->m_iRes[0],kerrad);

	mprintf("    Computing curl...\n");

	ComputeCurl();

	mprintf("    Computing divergence...\n");

	ComputeDivergence();

	mprintf("    Allocating...\n");

	pcPot = new CCube();
	pcPot->CopyFrom( g_pDivergence );

	pcAX = new CCube();
	pcAX->CopyFrom( g_pDivergence );

	pcAY = new CCube();
	pcAY->CopyFrom( g_pDivergence );

	pcAZ = new CCube();
	pcAZ->CopyFrom( g_pDivergence );

	pcRAX = new CCube();
	pcRAX->CopyFrom( g_pDivergence );

	pcRAY = new CCube();
	pcRAY->CopyFrom( g_pDivergence );

	pcRAZ = new CCube();
	pcRAZ->CopyFrom( g_pDivergence );

	pcGPotX = new CCube();
	pcGPotX->CopyFrom( g_pDivergence );

	pcGPotY = new CCube();
	pcGPotY->CopyFrom( g_pDivergence );

	pcGPotZ = new CCube();
	pcGPotZ->CopyFrom( g_pDivergence );

	pcOutX = new CCube();
	pcOutX->CopyFrom( g_pDivergence );

	pcOutY = new CCube();
	pcOutY->CopyFrom( g_pDivergence );

	pcOutZ = new CCube();
	pcOutZ->CopyFrom( g_pDivergence );

	pcKer = new CCube();
	pcKer->m_iRes[0] = 2*kerrad+1;
	pcKer->m_iRes[1] = 2*kerrad+1;
	pcKer->m_iRes[2] = 2*kerrad+1;
	pcKer->m_fStrideA[0] = g_pDivergence->m_fStrideA[0];
	pcKer->m_fStrideB[1] = g_pDivergence->m_fStrideB[1];
	pcKer->m_fStrideC[2] = g_pDivergence->m_fStrideC[2];
	pcKer->m_fCenter[0] = -kerrad*pcKer->m_fStrideA[0];
	pcKer->m_fCenter[1] = -kerrad*pcKer->m_fStrideB[1];
	pcKer->m_fCenter[2] = -kerrad*pcKer->m_fStrideC[2];
	pcKer->m_pBin = new double[pcKer->m_iRes[0]*pcKer->m_iRes[1]*pcKer->m_iRes[2]];

	mprintf("    Forming kernel...\n");

	for (iz=0;iz<pcKer->m_iRes[2];iz++) {
		rzs = pow2(iz-kerrad);
		for (iy=0;iy<pcKer->m_iRes[1];iy++) {
			rys = pow2(iy-kerrad);
			for (ix=0;ix<pcKer->m_iRes[0];ix++) {
				rxs = pow2(ix-kerrad);
				pcKer->m_pBin[iz*pcKer->m_iRes[0]*pcKer->m_iRes[1]+iy*pcKer->m_iRes[0]+ix] = 1.0 / (sqrt( rxs + rys + rzs ) + eps );
			}
		}
	}

	mprintf("    Computing Potential...\n");

	Convolution3D( true, false, g_pDivergence, pcKer, pcPot );
	pcPot->Multiply( -1.0 / (4.0 * Pi) );

	mprintf("    Computing A_X...\n");

	Convolution3D( true, false, g_pCurl[0], pcKer, pcAX );
	pcAX->Multiply( 1.0 / (4.0 * Pi) );

	mprintf("    Computing A_Y...\n");

	Convolution3D( true, false, g_pCurl[1], pcKer, pcAY );
	pcAY->Multiply( 1.0 / (4.0 * Pi) );

	mprintf("    Computing A_Z...\n");

	Convolution3D( true, false, g_pCurl[2], pcKer, pcAZ );
	pcAZ->Multiply( 1.0 / (4.0 * Pi) );

	mprintf("    Computing rot(A)...\n");

	ComputeCurl( pcAX, pcAY, pcAZ, pcRAX, pcRAY, pcRAZ );

	mprintf("    Computing grad(Pot)...\n");

	ComputeGradient( pcPot, pcGPotX, pcGPotY, pcGPotZ );

	mprintf("    Computing output...\n");

	for (iz=0;iz<g_pDivergence->m_iRes[0]*g_pDivergence->m_iRes[1]*g_pDivergence->m_iRes[2];iz++) {
		pcOutX->m_pBin[iz] = pcGPotX->m_pBin[iz] + pcRAX->m_pBin[iz];
		pcOutY->m_pBin[iz] = pcGPotY->m_pBin[iz] + pcRAY->m_pBin[iz];
		pcOutZ->m_pBin[iz] = pcGPotZ->m_pBin[iz] + pcRAZ->m_pBin[iz];
	}

	mprintf("    Checking...\n");

	tama = 0;
	tdma = 0;
	for (iz=0;iz<g_pDivergence->m_iRes[0]*g_pDivergence->m_iRes[1]*g_pDivergence->m_iRes[2];iz++) {

		if (fabs(g_pVectorField[0]->m_pBin[iz]) > tama)
			tama = fabs(g_pVectorField[0]->m_pBin[iz]);
		if (fabs(g_pVectorField[0]->m_pBin[iz] - pcOutX->m_pBin[iz]) > tdma)
			tdma = fabs(g_pVectorField[0]->m_pBin[iz] - pcOutX->m_pBin[iz]);

		if (fabs(g_pVectorField[1]->m_pBin[iz]) > tama)
			tama = fabs(g_pVectorField[1]->m_pBin[iz]);
		if (fabs(g_pVectorField[1]->m_pBin[iz] - pcOutY->m_pBin[iz]) > tdma)
			tdma = fabs(g_pVectorField[1]->m_pBin[iz] - pcOutY->m_pBin[iz]);

		if (fabs(g_pVectorField[2]->m_pBin[iz]) > tama)
			tama = fabs(g_pVectorField[2]->m_pBin[iz]);
		if (fabs(g_pVectorField[2]->m_pBin[iz] - pcOutZ->m_pBin[iz]) > tdma)
			tdma = fabs(g_pVectorField[2]->m_pBin[iz] - pcOutZ->m_pBin[iz]);
	}

	mprintf("      Largest absolute vector component: %20.10E\n",tama);
	mprintf("      Largest absolute deviation:        %20.10E\n",tdma);

	mprintf("    Writing all outputs...\n");
	g_pVectorField[0]->WriteCube("helmholtz_input_x.cube");
	g_pVectorField[1]->WriteCube("helmholtz_input_y.cube");
	g_pVectorField[2]->WriteCube("helmholtz_input_z.cube");
	pcPot->WriteCube("helmholtz_pot.cube");
	pcGPotX->WriteCube("helmholtz_grad_pot_x.cube");
	pcGPotY->WriteCube("helmholtz_grad_pot_y.cube");
	pcGPotZ->WriteCube("helmholtz_grad_pot_z.cube");
	pcAX->WriteCube("helmholtz_a_x.cube");
	pcAY->WriteCube("helmholtz_a_y.cube");
	pcAZ->WriteCube("helmholtz_a_z.cube");
	pcRAX->WriteCube("helmholtz_rot_a_x.cube");
	pcRAY->WriteCube("helmholtz_rot_a_y.cube");
	pcRAZ->WriteCube("helmholtz_rot_a_z.cube");
	pcOutX->WriteCube("helmholtz_output_x.cube");
	pcOutY->WriteCube("helmholtz_output_y.cube");
	pcOutZ->WriteCube("helmholtz_output_z.cube");
	pcKer->WriteCube("helmholtz_kernel.cube");

	mprintf("    Deallocating...\n");

	delete pcKer;
	delete pcPot;
	delete pcAX;
	delete pcAY;
	delete pcAZ;
	delete pcRAX;
	delete pcRAY;
	delete pcRAZ;
	delete pcGPotX;
	delete pcGPotY;
	delete pcGPotZ;
	delete pcOutX;
	delete pcOutY;
	delete pcOutZ;

	mprintf("*** All done.\n\n");
}



void Helmholtz_Decomp_Old(double eps) {

	int ix, iy, iz, ix2, iy2, iz2, rxs, rys, rzs;
	double tf;
	CCube *pcGPotX, *pcGPotY, *pcGPotZ;


	mprintf("*** Old Helmholtz decomposition (experimental) ***\n\n");

	mprintf("    Computing divergence...\n");

	ComputeDivergence();

	if (g_pCubeWork != NULL)
		delete g_pCubeWork;
	g_pCubeWork = new CCube();
	g_pCubeWork->CopyFrom( g_pDivergence );


//	mprintf(WHITE,"      [");

	mprintf("\n  Computing Helmholtz decomposition...\n");

	for (iz=0;iz<g_pDivergence->m_iRes[2];iz++) {

		for (iy=0;iy<g_pDivergence->m_iRes[1];iy++) {

			for (ix=0;ix<g_pDivergence->m_iRes[0];ix++) {

				tf = 0;

				for (iz2=0;iz2<g_pDivergence->m_iRes[2];iz2++) {

					rzs = SQR( iz2-iz );

					for (iy2=0;iy2<g_pDivergence->m_iRes[1];iy2++) {

						rys = SQR( iy2-iy );

						for (ix2=0;ix2<g_pDivergence->m_iRes[0];ix2++) {

							rxs = SQR( ix2-ix );

							tf += g_pDivergence->m_pBin[iz2*g_pDivergence->m_iRes[1]*g_pDivergence->m_iRes[0]+iy2*g_pDivergence->m_iRes[0]+ix2] / (sqrt( rzs + rys + rxs ) + eps);
						}

					}
				}

				tf *= -1.0 / (4.0 * Pi);

				g_pCubeWork->m_pBin[iz*g_pDivergence->m_iRes[1]*g_pDivergence->m_iRes[0]+iy*g_pDivergence->m_iRes[0]+ix] = tf;
			}

			mprintf("      %6d / %6d = %8.4f %% ...\n",iz*g_pDivergence->m_iRes[1]+iy,g_pDivergence->m_iRes[1]*g_pDivergence->m_iRes[2],((double)iz*g_pDivergence->m_iRes[1]+iy)/(g_pDivergence->m_iRes[1]*g_pDivergence->m_iRes[2])*100.0);

	//		if (fmod(iz*g_pDivergence->m_iRes[1]+iy,(double)g_pDivergence->m_iRes[1]*g_pDivergence->m_iRes[2]/60.0) < 1.0)
	//			mprintf(WHITE,"#");
		}
	}

//	mprintf(WHITE,"] Done.\n");


	mprintf("\n  Writing resulting potential to \"helmholtz_pot.cube\" ...\n");

	g_pCubeWork->WriteCube( "helmholtz_pot.cube" );

	mprintf("    Computing grad(Pot)...\n");

	pcGPotX = new CCube();
	pcGPotX->CopyFrom( g_pDivergence );

	pcGPotY = new CCube();
	pcGPotY->CopyFrom( g_pDivergence );

	pcGPotZ = new CCube();
	pcGPotZ->CopyFrom( g_pDivergence );

	ComputeGradient( g_pCubeWork, pcGPotX, pcGPotY, pcGPotZ );

	mprintf("\n  Writing resulting fields to \"helmholtz_field{x,y,z}.cube\" ...\n");

	pcGPotX->WriteCube( "helmholtz_fieldx.cube" );
	pcGPotY->WriteCube( "helmholtz_fieldy.cube" );
	pcGPotZ->WriteCube( "helmholtz_fieldz.cube" );

	delete pcGPotX;
	delete pcGPotY;
	delete pcGPotZ;

	mprintf("\n  All done.\n\n");
}



void VectorFieldFromGradient( CCube *tcu ) {

	int z, ix, iy, iz;


	mprintf("    Constructing vector field from gradient of volumetric data...\n" );

	for (z=0;z<3;z++) {
		if (g_pVectorField[z] != NULL)
			delete g_pVectorField[z];
		g_pVectorField[z] = new CCube();
		g_pVectorField[z]->CopyFrom( tcu );
		for (iz=0;iz<tcu->m_iRes[0]*tcu->m_iRes[1]*tcu->m_iRes[2];iz++)
			g_pVectorField[z]->m_pBin[iz] = 0;
	}

	mprintf(WHITE,"      [");

	for (iz=1;iz<tcu->m_iRes[2]-1;iz++) {

		for (iy=1;iy<tcu->m_iRes[1]-1;iy++) {

			for (ix=1;ix<tcu->m_iRes[0]-1;ix++) {

				g_pVectorField[0]->m_pBin[iz*tcu->m_iRes[0]*tcu->m_iRes[1]+iy*tcu->m_iRes[0]+ix] = (tcu->m_pBin[iz*tcu->m_iRes[0]*tcu->m_iRes[1]+iy*tcu->m_iRes[0]+(ix+1)] - tcu->m_pBin[iz*tcu->m_iRes[0]*tcu->m_iRes[1]+iy*tcu->m_iRes[0]+(ix-1)]) / 2.0;
				g_pVectorField[1]->m_pBin[iz*tcu->m_iRes[0]*tcu->m_iRes[1]+iy*tcu->m_iRes[0]+ix] = (tcu->m_pBin[iz*tcu->m_iRes[0]*tcu->m_iRes[1]+(iy+1)*tcu->m_iRes[0]+ix] - tcu->m_pBin[iz*tcu->m_iRes[0]*tcu->m_iRes[1]+(iy-1)*tcu->m_iRes[0]+ix]) / 2.0;
				g_pVectorField[2]->m_pBin[iz*tcu->m_iRes[0]*tcu->m_iRes[1]+iy*tcu->m_iRes[0]+ix] = (tcu->m_pBin[(iz+1)*tcu->m_iRes[0]*tcu->m_iRes[1]+iy*tcu->m_iRes[0]+ix] - tcu->m_pBin[(iz-1)*tcu->m_iRes[0]*tcu->m_iRes[1]+iy*tcu->m_iRes[0]+ix]) / 2.0;
			}

			if (fmod(iz*tcu->m_iRes[1]+iy,(double)tcu->m_iRes[1]*tcu->m_iRes[2]/60.0) < 1.0)
				mprintf(WHITE,"#");
		}
	}

	mprintf(WHITE,"] Done.\n");
}



void Usage() {

	mprintf("\n");
	mprintf(WHITE,"Command Line Usage:\n\n");

	mprintf("    -readcube <input.cube>                   Reads a cube file into the working buffer.\n");
	mprintf("    -newcube <rx> <ry> <rz> <sx> <sy> <sz>   Creates a new working buffer with resolution (rx,ry,rz) and extent (sx,sy,sz).\n");
	mprintf("    -setoffset <ox> <oy> <oz>                Sets the grid offset, i.e. the position of voxel (0|0|0) in space, to (ox|oy|oz).\n");
	mprintf("    -addcube <input.cube>                    Adds a cube file to the working buffer.\n");
	mprintf("    -subtractcube <input.cube>               Subtracts a cube file from the working buffer.\n");
	mprintf("    -writecube <output.cube>                 Writes the contents of the working buffer to a cube file.\n");

	mprintf("\n");
	mprintf("    -vectorfield <vx> <vy> <vz>              Reads three cube files with x, y, and z component of a vector field.\n");
	mprintf("    -vectorfromgrad                          Computes a vector field as the gradient of the working buffer.\n");
	mprintf("    -pushtovectorfield {x|y|z}               Copies the working buffer to component x, y, or z of the vector field.\n");
	mprintf("    -pullfromvectorfield {x|y|z}             Copies component x, y, or z of the vector field into the working buffer.\n");
	mprintf("    -vectorabs                               Fills the working buffer with the absolute values of the vector field.\n");

	mprintf("\n");
	mprintf("    -add <const>                             Adds a constant to the volumetric data in the working buffer.\n");
	mprintf("    -multiply <const>                        Multiplies the volumetric data in the working buffer with a constant.\n");
	mprintf("    -addgauss <h> <x> <y> <z> <s>            Adds a Gaussian hill with height h at (x,y,z) with width s^2 to the working buffer.\n");
	mprintf("    -addlayer <dx> <dy> <dz> <v>             Adds a layer of voxel width (dx,dy,dz) with value v at each side of the volumetric data.\n");
	mprintf("    -limitrange <min> <max>                  Limits the range of the volumetric data to [min,max] by clipping.\n");
	mprintf("    -crop <x1> <y1> <z1> <x2> <y2> <z2>      Crops the working buffer and the vector field to the grid point range (x1,y1,z1) to (x2,y2,z2).\n");
	mprintf("    -resample <type> <rx> <ry> <rz>          Resamples the working buffer and the vector field to a new resolution (rx,ry,rz). type=\"linear\".\n");
	mprintf("    -envelope <type> <inner>                 Multiplies the working buffer and vector field with an envelope. type=\"cossq\". inner=fraction of inner grid points (0-1).\n");

	mprintf("\n");
	mprintf("    -normal{x|y|z}                           Defines a cut plane with normal vector along x / y / z direction.\n");
	mprintf("    -planesection <x1> <x2> <y1> <y2>        Uses only the given range of points from the volumetric data for the cut plane.\n");
	mprintf("    -valuerange <min> <max>                  Use the specified value range for the cut plane plot (automatically determined by default).\n");
	mprintf("    -slice <n>                               Takes the n-th slice along the normal vector from the CUBE file.\n");
	mprintf("    -slicerange <n1> <n2>                    Averages from n1-th to n2-th slice along the normal vector.\n");
	mprintf("    -showatoms <r>                           Displays the atoms in the cut plane plot with radius r (default 1.0).\n");
	mprintf("    -vectorrange <max>                       Use the specified max. vector length for the cut plane plot (automatically determined by default).\n");
	mprintf("    -vectorscale <f>                         Scales vectors of the vector field with factor f for the cut plane (default 1.0).\n");
	mprintf("    -vectorstride <n>                        Shows only every n-th vector of the vector field in the cut plane (default 1).\n");
	mprintf("    -vectorhead <f>                          Scales the head size of the vector arrows in the cut plane (default 1.0).\n");
	mprintf("    -nogrid                                  Disables the grid in the cut plane plot.\n");
	mprintf("    -noticks                                 Disables ticks and tick labels in the cut plane plot.\n");
	mprintf("    -nointerpolate                           Disables interpolation between the values in the cut plane plot.\n");
	mprintf("    -ncontour <n>                            Sets the number of contour lines to n. A value of 0 disables contour lines.\n");
	mprintf("    -aspectratio <f>                         Sets the aspect ratio to f=height/width.\n");
	mprintf("    -plotexpo <e>                            Sets the plot exponent to e. A value of 1 means a linear plot (default).\n");
	mprintf("    -plotpixel <x>                           Sets the pixel width of the plot to x.\n");
	mprintf("    -writeplane <output>                     Writes the contents of the previously defined cut plane.\n");

	mprintf("\n");
}



void TestConvolution() {

	int resx, resy, resz, kernelsize;
	int ix, iy, iz, c1, c2;
	CCube *pin, *pker, *pout;
	double px, py, pz, tf, tsum;


	resx = 64;
	resy = 64;
	resz = 64;
	kernelsize = 24;
	c1 = 16;
	c2 = 48;

	pin = new CCube();
	pin->m_iRes[0] = resx;
	pin->m_iRes[1] = resy;
	pin->m_iRes[2] = resz;
	pin->m_fStrideA[0] = 20.0/resx;
	pin->m_fStrideB[1] = 20.0/resy;
	pin->m_fStrideC[2] = 20.0/resz;
	pin->m_fCenter[0] = -10.0;
	pin->m_fCenter[1] = -10.0;
	pin->m_fCenter[2] = -10.0;
	pin->m_pBin = new double[resx*resy*resz];

	pout = new CCube();
	pout->CopyFrom( pin );

	pker = new CCube();
	pker->m_iRes[0] = 2*kernelsize+1;
	pker->m_iRes[1] = 2*kernelsize+1;
	pker->m_iRes[2] = 2*kernelsize+1;
	pker->m_fStrideA[0] = 20.0/resx;
	pker->m_fStrideB[1] = 20.0/resy;
	pker->m_fStrideC[2] = 20.0/resz;
	pker->m_fCenter[0] = -kernelsize*20.0/resx;
	pker->m_fCenter[1] = -kernelsize*20.0/resy;
	pker->m_fCenter[2] = -kernelsize*20.0/resz;
	pker->m_pBin = new double[pker->m_iRes[0]*pker->m_iRes[1]*pker->m_iRes[2]];

	for (iz=0;iz<resz;iz++) {
		for (iy=0;iy<resy;iy++) {
			for (ix=0;ix<resx;ix++) {
				if ((iz >= c1) && (iz <= c2) && (iy >= c1) && (iy <= c2) && (ix >= c1) && (ix <= c2))
					pin->m_pBin[iz*resy*resx+iy*resx+ix] = 1.0;
				else
					pin->m_pBin[iz*resy*resx+iy*resx+ix] = 0;
			}
		}
	}

	tsum = 0;
	for (iz=0;iz<2*kernelsize+1;iz++) {
		pz = pow2(((double)iz-kernelsize)/kernelsize);
		for (iy=0;iy<2*kernelsize+1;iy++) {
			py = pow2(((double)iy-kernelsize)/kernelsize);
			for (ix=0;ix<2*kernelsize+1;ix++) {
				px = pow2(((double)ix-kernelsize)/kernelsize);
				tf = exp( -(px+py+pz)*5.0 );
				tsum += tf;
				pker->m_pBin[iz*pker->m_iRes[0]*pker->m_iRes[1]+iy*pker->m_iRes[0]+ix] = tf;
			}
		}
	}

	// Normalize kernel integral to 1
	for (ix=0;ix<pker->m_iRes[0]*pker->m_iRes[1]*pker->m_iRes[2];ix++)
		pker->m_pBin[ix] /= tsum;

	Convolution3D( true, true, pin, pker, pout );
}



class CElectrostaticObject {
public:

	CElectrostaticObject( CCube *c ) {
		m_pPresence = new CCube();
		m_pPresence->CopyFrom(c);
		for (int z=0;z<c->m_iRes[0]*c->m_iRes[1]*c->m_iRes[2];z++)
			m_pPresence->m_pBin[z] = 0;
	}

	~CElectrostaticObject() {
		if (m_pPresence != NULL) {
			delete m_pPresence;
			m_pPresence = NULL;
		}
	}

	void AddSphere( double posx, double posy, double posz, double radius );

	void AddBox( double px1, double py1, double pz1, double px2, double py2, double pz2 );

	double m_fEpsilon;
	bool m_bFixPotential;
	double m_fPotential;

	CCube *m_pPresence;
};



void CElectrostaticObject::AddSphere( double posx, double posy, double posz, double radius ) {

	int ix, iy, iz;
	double px, py, pz, eps, tf, tf2;


	eps = m_pPresence->m_fStrideA[0] * sqrt(2.0);

	for (iz=0;iz<m_pPresence->m_iRes[2];iz++) {

		pz = m_pPresence->m_fCenter[2] + iz*m_pPresence->m_fStrideC[2];

		for (iy=0;iy<m_pPresence->m_iRes[1];iy++) {

			py = m_pPresence->m_fCenter[1] + iy*m_pPresence->m_fStrideB[1];

			for (ix=0;ix<m_pPresence->m_iRes[0];ix++) {

				px = m_pPresence->m_fCenter[0] + ix*m_pPresence->m_fStrideA[0];

				tf = sqrt(pow2(px-posx)+pow2(py-posy)+pow2(pz-posz));

				if (tf < radius+eps/2.0) {

					if (tf > radius-eps/2.0)
						tf2 = 1.0 - (tf - radius + eps/2.0)/eps;
					else
						tf2 = 1.0;

					m_pPresence->m_pBin[iz*m_pPresence->m_iRes[1]*m_pPresence->m_iRes[0]+iy*m_pPresence->m_iRes[0]+ix] = tf2;
				}
			}
		}
	}
}



void CElectrostaticObject::AddBox( double px1, double py1, double pz1, double px2, double py2, double pz2 ) {

	int ix, iy, iz;
	double px, py, pz, vx, vy, vz;


	for (iz=0;iz<m_pPresence->m_iRes[2];iz++) {

		pz = m_pPresence->m_fCenter[2] + iz*m_pPresence->m_fStrideC[2];
		
		if ((pz > pz1-m_pPresence->m_fStrideC[2]/2.0) && (pz <= pz1+m_pPresence->m_fStrideC[2]/2.0))
			vz = ((pz - (pz1-m_pPresence->m_fStrideC[2]/2.0)) / m_pPresence->m_fStrideC[2]);
		else if ((pz > pz1+m_pPresence->m_fStrideC[2]/2.0) && (pz <= pz2-m_pPresence->m_fStrideC[2]/2.0))
			vz = 1.0;
		else if ((pz > pz2-m_pPresence->m_fStrideC[2]/2.0) && (pz <= pz2+m_pPresence->m_fStrideC[2]/2.0))
			vz = 1.0 - ((pz - (pz2-m_pPresence->m_fStrideC[2]/2.0)) / m_pPresence->m_fStrideC[2]);
		else
			vz = 0;

		for (iy=0;iy<m_pPresence->m_iRes[1];iy++) {

			py = m_pPresence->m_fCenter[1] + iy*m_pPresence->m_fStrideB[1];

			if ((py > py1-m_pPresence->m_fStrideB[1]/2.0) && (py <= py1+m_pPresence->m_fStrideB[1]/2.0))
				vy = ((py - (py1-m_pPresence->m_fStrideB[1]/2.0)) / m_pPresence->m_fStrideB[1]);
			else if ((py > py1+m_pPresence->m_fStrideB[1]/2.0) && (py <= py2-m_pPresence->m_fStrideB[1]/2.0))
				vy = 1.0;
			else if ((py > py2-m_pPresence->m_fStrideB[1]/2.0) && (py <= py2+m_pPresence->m_fStrideB[1]/2.0))
				vy = 1.0 - ((py - (py2-m_pPresence->m_fStrideB[1]/2.0)) / m_pPresence->m_fStrideB[1]);
			else
				vy = 0;

			for (ix=0;ix<m_pPresence->m_iRes[0];ix++) {

				px = m_pPresence->m_fCenter[0] + ix*m_pPresence->m_fStrideA[0];

				if ((px > px1-m_pPresence->m_fStrideA[0]/2.0) && (px <= px1+m_pPresence->m_fStrideA[0]/2.0))
					vx = ((px - (px1-m_pPresence->m_fStrideA[0]/2.0)) / m_pPresence->m_fStrideA[0]);
				else if ((px > px1+m_pPresence->m_fStrideA[0]/2.0) && (px <= px2-m_pPresence->m_fStrideA[0]/2.0))
					vx = 1.0;
				else if ((px > px2-m_pPresence->m_fStrideA[0]/2.0) && (px <= px2+m_pPresence->m_fStrideA[0]/2.0))
					vx = 1.0 - ((px - (px2-m_pPresence->m_fStrideA[0]/2.0)) / m_pPresence->m_fStrideA[0]);
				else
					vx = 0;

				m_pPresence->m_pBin[iz*m_pPresence->m_iRes[1]*m_pPresence->m_iRes[0]+iy*m_pPresence->m_iRes[0]+ix] = vx * vy * vz;
			}
		}
	}
}



#ifdef USE_CUBETOOL_FFTW



class CElectrostaticSolver {
public:


	CElectrostaticSolver() : m_pPotential(NULL), m_pCharge(NULL), m_pFTBufWork(NULL), m_pFTBufInOut(NULL),
		m_plComplexForward(NULL), m_plComplexBackward(NULL) {
	}


	~CElectrostaticSolver() {
		if (m_pPotential != NULL) {
			delete m_pPotential;
			m_pPotential = NULL;
		}
		if (m_pCharge != NULL) {
			delete m_pCharge;
			m_pCharge = NULL;
		}
		if (m_pFTBufWork != NULL) {
			fftw_free( m_pFTBufWork );
			m_pFTBufWork = NULL;
		}
		if (m_pFTBufInOut != NULL) {
			fftw_free( m_pFTBufInOut );
			m_pFTBufInOut = NULL;
		}
		if (m_plComplexForward != NULL) {
			fftw_destroy_plan( m_plComplexForward );
			m_plComplexForward = NULL;
		}
		if (m_plComplexBackward != NULL) {
			fftw_destroy_plan( m_plComplexBackward );
			m_plComplexBackward = NULL;
		}
	}


	void InitSolver( int rx, int ry, int rz, double sx, double sy, double sz, double ox, double oy, double oz );


	void AddObject( double pot ); 


	void SolvePoisson( CCube *density, CCube *potential );


	CCube *m_pPotential;
	CCube *m_pPotentialGradX;
	CCube *m_pPotentialGradY;
	CCube *m_pPotentialGradZ;
	CCube *m_pCharge;
	CCube *m_pEpsilon;
	CCube *m_pEpsilonLog;
	CCube *m_pEpsilonLogGradX;
	CCube *m_pEpsilonLogGradY;
	CCube *m_pEpsilonLogGradZ;


	fftw_complex *m_pFTBufWork;
	double       *m_pFTBufInOut;

	fftw_plan m_plComplexForward;
	fftw_plan m_plComplexBackward;

	std::vector<CElectrostaticObject*> m_oaObjects;
};



void CElectrostaticSolver::AddObject( double pot ) {

	CElectrostaticObject *obj;

	obj = new CElectrostaticObject( m_pPotential );
	obj->m_fPotential = pot;

	m_oaObjects.push_back(obj);
}



void CElectrostaticSolver::InitSolver( int rx, int ry, int rz, double sx, double sy, double sz, double ox, double oy, double oz ) {

	int z;


	m_pPotential = new CCube();
	m_pPotential->m_iRes[0] = rx;
	m_pPotential->m_iRes[1] = ry;
	m_pPotential->m_iRes[2] = rz;
	m_pPotential->m_fCenter[0] = ox;
	m_pPotential->m_fCenter[1] = oy;
	m_pPotential->m_fCenter[2] = oz;
	m_pPotential->m_fStrideA[0] = sx/rx;
	m_pPotential->m_fStrideB[1] = sy/ry;
	m_pPotential->m_fStrideC[2] = sz/rz;
	m_pPotential->m_pBin = new double[rx*ry*rz];
	for (z=0;z<rx*ry*rz;z++)
		m_pPotential->m_pBin[z] = 0;

	m_pPotentialGradX = new CCube();
	m_pPotentialGradX->CopyFrom( m_pPotential );

	m_pPotentialGradY = new CCube();
	m_pPotentialGradY->CopyFrom( m_pPotential );

	m_pPotentialGradZ = new CCube();
	m_pPotentialGradZ->CopyFrom( m_pPotential );

	m_pCharge = new CCube();
	m_pCharge->CopyFrom( m_pPotential );

	m_pEpsilon = new CCube();
	m_pEpsilon->CopyFrom( m_pPotential );

	m_pEpsilonLog = new CCube();
	m_pEpsilonLog->CopyFrom( m_pPotential );

	m_pEpsilonLogGradX = new CCube();
	m_pEpsilonLogGradX->CopyFrom( m_pPotential );

	m_pEpsilonLogGradY = new CCube();
	m_pEpsilonLogGradY->CopyFrom( m_pPotential );

	m_pEpsilonLogGradZ = new CCube();
	m_pEpsilonLogGradZ->CopyFrom( m_pPotential );

	m_pFTBufWork  = (fftw_complex*)fftw_malloc( sizeof(fftw_complex) * m_pPotential->m_iRes[0] * m_pPotential->m_iRes[1] * (m_pPotential->m_iRes[2]/2+1) );
	m_pFTBufInOut =       (double*)fftw_malloc( sizeof(double)       * m_pPotential->m_iRes[0] * m_pPotential->m_iRes[1] *  m_pPotential->m_iRes[2] );

	m_plComplexForward  = fftw_plan_dft_r2c_3d( m_pPotential->m_iRes[0], m_pPotential->m_iRes[1], m_pPotential->m_iRes[2], m_pFTBufInOut, m_pFTBufWork, FFTW_ESTIMATE );
	m_plComplexBackward = fftw_plan_dft_c2r_3d( m_pPotential->m_iRes[0], m_pPotential->m_iRes[1], m_pPotential->m_iRes[2], m_pFTBufWork, m_pFTBufInOut, FFTW_ESTIMATE );

	for (z=0;z<m_pPotential->m_iRes[0]*m_pPotential->m_iRes[1]*m_pPotential->m_iRes[2];z++) {
		m_pPotential->m_pBin[z] = 0;
		m_pCharge->m_pBin[z] = 0;
		m_pEpsilon->m_pBin[z] = 1.0;
	}
}



void CElectrostaticSolver::SolvePoisson( CCube *density, CCube *potential ) {

	int i, z, ix, iy, iz;
	double wzr, wzi, wyr, wyi, wxr, wxi;
	double rwzr, rwzi, rwyr, rwyi, rwxr, rwxi;
	double denomr, denomi, tf;


	for (z=0;z<m_pPotential->m_iRes[0]*m_pPotential->m_iRes[1]*m_pPotential->m_iRes[2];z++)
		m_pFTBufInOut[z] = density->m_pBin[z];

	fftw_execute( m_plComplexForward );

	i = 0;
	for (iz=0;iz<m_pPotential->m_iRes[2];iz++) {

		wzr = cos(2.0*Pi*iz/(double)m_pPotential->m_iRes[2]);
		wzi = sin(2.0*Pi*iz/(double)m_pPotential->m_iRes[2]);
		rwzr =  wzr / (wzr*wzr + wzi*wzi);
		rwzi = -wzi / (wzr*wzr + wzi*wzi);

		for (iy=0;iy<m_pPotential->m_iRes[1];iy++) {

			wyr = cos(2.0*Pi*iy/(double)m_pPotential->m_iRes[1]);
			wyi = sin(2.0*Pi*iy/(double)m_pPotential->m_iRes[1]);
			rwyr =  wyr / (wyr*wyr + wyi*wyi);
			rwyi = -wyi / (wyr*wyr + wyi*wyi);

			for (ix=0;ix<m_pPotential->m_iRes[0]/2+1;ix++) {

				wxr = cos(2.0*Pi*ix/(double)m_pPotential->m_iRes[0]);
				wxi = sin(2.0*Pi*ix/(double)m_pPotential->m_iRes[0]);
				rwxr =  wxr / (wxr*wxr + wxi*wxi);
				rwxi = -wxi / (wxr*wxr + wxi*wxi);

				denomr = 6.0 - wxr - wyr - wzr - rwxr - rwyr - rwzr;
				denomi =     - wxi - wyi - wzi - rwxi - rwyi - rwzi;

				if ((denomr != 0) || (denomi != 0)) {

					tf = (m_pFTBufWork[i][0]*denomr + m_pFTBufWork[i][1]*denomi) / (denomr*denomr + denomi*denomi);
					m_pFTBufWork[i][1] = (m_pFTBufWork[i][1]*denomr - m_pFTBufWork[i][0]*denomi) / (denomr*denomr + denomi*denomi);
					m_pFTBufWork[i][0] = tf;

				} else {

					m_pFTBufWork[i][0] = 0;
					m_pFTBufWork[i][1] = 0;
				}
				i++;
			}
		}
	}

	fftw_execute( m_plComplexBackward );

	for (z=0;z<m_pPotential->m_iRes[0]*m_pPotential->m_iRes[1]*m_pPotential->m_iRes[2];z++)
		potential->m_pBin[z] = m_pFTBufInOut[z] / (8.0*m_pPotential->m_iRes[0]*m_pPotential->m_iRes[1]*m_pPotential->m_iRes[2]);
}



void TestPoisson() {

	// https://my.ece.utah.edu/~ece6340/LECTURES/Feb1/Nagel%202012%20-%20Solving%20the%20Generalized%20Poisson%20Equation%20using%20FDM.pdf
	// https://de.wikipedia.org/wiki/SOR-Verfahren
	// https://de.wikipedia.org/wiki/Gau%C3%9F-Seidel-Verfahren
	// !! https://arxiv.org/pdf/1509.00680.pdf

	CElectrostaticSolver *esolv;
	CElectrostaticObject *obj;
	CCube *pot, *chg;
	int z, z2, i;
	double devi, tf, eta, fmi, fma;
	char buf[256];


	mprintf("*** Poisson Test ***\n");

	eta = 0.5;

	esolv = new CElectrostaticSolver();

	esolv->InitSolver( 128, 128, 128, 20.0, 20.0, 20.0, -10.0, -10.0, -10.0 );

	esolv->AddObject( 2.0 );
	esolv->m_oaObjects.back()->m_bFixPotential = true;
	esolv->m_oaObjects.back()->m_fEpsilon = 100000.0;
	esolv->m_oaObjects.back()->AddSphere( -5.0, 0.0, 0.0, 3.0 );

	esolv->AddObject( -2.0 );
	esolv->m_oaObjects.back()->m_bFixPotential = true;
	esolv->m_oaObjects.back()->m_fEpsilon = 100000.0;
	esolv->m_oaObjects.back()->AddSphere( 5.0, 0.0, 0.0, 2.0 );

	esolv->AddObject( 0.0 );
	esolv->m_oaObjects.back()->m_bFixPotential = false;
	esolv->m_oaObjects.back()->m_fEpsilon = 100.0;
	esolv->m_oaObjects.back()->AddBox( -1.5, -1.0, 0.0, 1.5, 1.0, 2.0 );

	i = 0;

	pot = esolv->m_pPotential;
	chg = esolv->m_pCharge;

	for (z2=0;z2<(int)esolv->m_oaObjects.size();z2++) {
		mprintf("    Object %d:\n",z2+1);
		mprintf("      Epsilon: %E\n",esolv->m_oaObjects[z2]->m_fEpsilon);
		fmi = 1.0e30;
		fma = -1.0e30;
		for (z=0;z<chg->m_iRes[0]*chg->m_iRes[1]*chg->m_iRes[2];z++) {
			if (esolv->m_oaObjects[z2]->m_pPresence->m_pBin[z] > fma)
				fma = esolv->m_oaObjects[z2]->m_pPresence->m_pBin[z];
			if (esolv->m_oaObjects[z2]->m_pPresence->m_pBin[z] < fmi)
				fmi = esolv->m_oaObjects[z2]->m_pPresence->m_pBin[z];
		}
		mprintf("      Presence: %E .. %E\n",fmi,fma);
	}


	for (z2=0;z2<(int)esolv->m_oaObjects.size();z2++)
		for (z=0;z<chg->m_iRes[0]*chg->m_iRes[1]*chg->m_iRes[2];z++)
			if (esolv->m_oaObjects[z2]->m_pPresence->m_pBin[z] != 0)
				esolv->m_pEpsilon->m_pBin[z] = MAX(esolv->m_pEpsilon->m_pBin[z],esolv->m_oaObjects[z2]->m_pPresence->m_pBin[z] * esolv->m_oaObjects[z2]->m_fEpsilon);

	fmi = 1.0e30;
	fma = -1.0e30;
	for (z=0;z<chg->m_iRes[0]*chg->m_iRes[1]*chg->m_iRes[2];z++) {
		if (esolv->m_pEpsilon->m_pBin[z] > fma)
			fma = esolv->m_pEpsilon->m_pBin[z];
		if (esolv->m_pEpsilon->m_pBin[z] < fmi)
			fmi = esolv->m_pEpsilon->m_pBin[z];
	}
	mprintf("    Epsilon values within %E .. %E\n",fmi,fma);

	for (z=0;z<chg->m_iRes[0]*chg->m_iRes[1]*chg->m_iRes[2];z++)
		esolv->m_pEpsilonLog->m_pBin[z] = log( esolv->m_pEpsilon->m_pBin[z] );

	fmi = 1.0e30;
	fma = -1.0e30;
	for (z=0;z<chg->m_iRes[0]*chg->m_iRes[1]*chg->m_iRes[2];z++) {
		if (esolv->m_pEpsilonLog->m_pBin[z] > fma)
			fma = esolv->m_pEpsilonLog->m_pBin[z];
		if (esolv->m_pEpsilonLog->m_pBin[z] < fmi)
			fmi = esolv->m_pEpsilonLog->m_pBin[z];
	}
	mprintf("    Log(Epsilon) values within %E .. %E\n",fmi,fma);

	ComputeGradient( esolv->m_pEpsilonLog, esolv->m_pEpsilonLogGradX, esolv->m_pEpsilonLogGradY, esolv->m_pEpsilonLogGradZ );

	while (true) {

		mprintf("    Iteration %d ...\n",i+1);

		mprintf("      Computing gradient of potential...\n");

		ComputeGradient( pot, esolv->m_pPotentialGradX, esolv->m_pPotentialGradY, esolv->m_pPotentialGradZ );

		mprintf("      Computing polarization charge...\n");

		fmi = 1.0e30;
		fma = -1.0e30;
		for (z=0;z<chg->m_iRes[0]*chg->m_iRes[1]*chg->m_iRes[2];z++) {
			tf = (1.0-eta)*chg->m_pBin[z] + eta*(esolv->m_pEpsilonLogGradX->m_pBin[z]*esolv->m_pPotentialGradX->m_pBin[z]+esolv->m_pEpsilonLogGradY->m_pBin[z]*esolv->m_pPotentialGradY->m_pBin[z]+esolv->m_pEpsilonLogGradZ->m_pBin[z]*esolv->m_pPotentialGradZ->m_pBin[z]);
			if (tf > fma)
				fma = tf;
			if (tf < fmi)
				fmi = tf;
			chg->m_pBin[z] = tf;
		}

		mprintf("      Charge values are within %E .. %E\n",fmi,fma);

		mprintf("      Enforcing potential on metallic objects...\n");

		devi = 0;

		for (z2=0;z2<(int)esolv->m_oaObjects.size();z2++) {

			obj = esolv->m_oaObjects[z2];

			if (!obj->m_bFixPotential)
				continue;

			for (z=0;z<chg->m_iRes[0]*chg->m_iRes[1]*chg->m_iRes[2];z++) {

				if (obj->m_pPresence->m_pBin[z] == 1.0) {
					tf = fabs(pot->m_pBin[z]-obj->m_fPotential);
					if (tf > devi)
						devi = tf;
				}

				if (obj->m_pPresence->m_pBin[z] != 0.0)
					chg->m_pBin[z] += 0.1 * obj->m_pPresence->m_pBin[z] * (obj->m_fPotential - pot->m_pBin[z]);
			}
		}

		mprintf("      Max. deviation: %20.10E\n",devi);

		mprintf("      Solving Poisson equation...\n");

		esolv->SolvePoisson( chg, pot );

		fmi = 1.0e30;
		fma = -1.0e30;
		for (z=0;z<chg->m_iRes[0]*chg->m_iRes[1]*chg->m_iRes[2];z++) {
			if (pot->m_pBin[z] > fma)
				fma = pot->m_pBin[z];
			if (pot->m_pBin[z] < fmi)
				fmi = pot->m_pBin[z];
		}

		mprintf("      Potential values are within %E .. %E\n",fmi,fma);

		if (((i%20) == 0) && (i != 0)) {
			mprintf("      Writing current state to cube files:\n");
			sprintf(buf,"poisson_pot_%04d.cube",i);
			mprintf("        Writing \"%s\"...\n",buf);
			pot->WriteCube(buf);
			sprintf(buf,"poisson_charge_%04d.cube",i);
			mprintf("        Writing \"%s\"...\n",buf);
			chg->WriteCube(buf);
		}

		i++;
	}

	delete esolv;

	mprintf("*** Poisson Test Finished ***\n");



/*	CCube *pob1, *pob2, *poball, *ppot, *ppotfix, *pres;
	int ix, iy, iz, iter;
	double px, py, pz, tf, tf2, eps, tfma;
	double rad1, pos1x, rad2, pos2x, vlt1, vlt2, omega;
	char buf[256];


	mprintf("*** Poisson Test ***\n");

	rad1 = 3.0;
	pos1x = -5.0;
	rad2 = 2.0;
	pos2x = 4.0;
	vlt1 = 2.0;
	vlt2 = 0;

	ppot = new CCube();
	ppot->m_iRes[0] = 128;
	ppot->m_iRes[1] = 128;
	ppot->m_iRes[2] = 128;
	ppot->m_fStrideA[0] = 20.0 / ppot->m_iRes[0];
	ppot->m_fStrideB[1] = 20.0 / ppot->m_iRes[1];
	ppot->m_fStrideC[2] = 20.0 / ppot->m_iRes[2];
	ppot->m_fCenter[0] = -10.0;
	ppot->m_fCenter[1] = -10.0;
	ppot->m_fCenter[2] = -10.0;
	ppot->m_pBin = new double[ppot->m_iRes[0]*ppot->m_iRes[1]*ppot->m_iRes[2]];

	for (iz=0;iz<ppot->m_iRes[0]*ppot->m_iRes[1]*ppot->m_iRes[2];iz++)
		ppot->m_pBin[iz] = 0;

	pob1 = new CCube();
	pob1->CopyFrom(ppot);

	pob2 = new CCube();
	pob2->CopyFrom(ppot);

	poball = new CCube();
	poball->CopyFrom(ppot);

	ppotfix = new CCube();
	ppotfix->CopyFrom(ppot);

	pres = new CCube();
	pres->CopyFrom(ppot);

	eps = ppot->m_fStrideA[0] * sqrt(2.0);

	for (iz=0;iz<pob1->m_iRes[2];iz++) {

		pz = pob1->m_fCenter[2] + iz*pob1->m_fStrideC[2];

		for (iy=0;iy<pob1->m_iRes[1];iy++) {

			py = pob1->m_fCenter[1] + iy*pob1->m_fStrideB[1];

			for (ix=0;ix<pob1->m_iRes[0];ix++) {

				px = pob1->m_fCenter[0] + ix*pob1->m_fStrideA[0];

				// Ob1

				tf = sqrt(pow2(px-pos1x)+pow2(py)+pow2(pz));

				if (tf < rad1+eps/2.0) {

					if (tf > rad1-eps/2.0)
						tf2 = 1.0 - (tf - rad1 + eps/2.0)/eps;
					else
						tf2 = 1.0;

					pob1->m_pBin[iz*pob1->m_iRes[1]*pob1->m_iRes[0]+iy*pob1->m_iRes[0]+ix] = tf2;
				}

				// Ob2

				tf = sqrt(pow2(px-pos2x)+pow2(py)+pow2(pz));

				if (tf < rad2+eps/2.0) {

					if (tf > rad2-eps/2.0)
						tf2 = 1.0 - (tf - rad2+ eps/2.0)/eps;
					else
						tf2 = 1.0;

					pob2->m_pBin[iz*pob1->m_iRes[1]*pob1->m_iRes[0]+iy*pob1->m_iRes[0]+ix] = tf2;
				}
			}
		}
	}

	for (iz=0;iz<pob1->m_iRes[0]*pob1->m_iRes[1]*pob1->m_iRes[2];iz++) {

		if (pob1->m_pBin[iz] != 0)
			ppotfix->m_pBin[iz] = vlt1;
		if (pob2->m_pBin[iz] != 0)
			ppotfix->m_pBin[iz] = vlt2;
		poball->m_pBin[iz] = pob1->m_pBin[iz] + pob2->m_pBin[iz];
		ppot->m_pBin[iz] = ppotfix->m_pBin[iz] * poball->m_pBin[iz];
	}

	pob1->WriteCube("poisson_ob1.cube");
	pob2->WriteCube("poisson_ob2.cube");
	poball->WriteCube("poisson_oball.cube");
	ppotfix->WriteCube("poisson_potfix.cube");
	ppot->WriteCube("poisson_pot.cube");

	omega = 0.2;
	iter = 0;

	while (true) {

		mprintf("    Iteration %d...\n",iter+1);

		mprintf("      Computing residues...\n");

		tfma = 0;
		for (iz=1;iz<pob1->m_iRes[2]-1;iz++) {
			for (iy=1;iy<pob1->m_iRes[1]-1;iy++) {
				for (ix=1;ix<pob1->m_iRes[0]-1;ix++) {
					tf =  ppot->m_pBin[iz*ppot->m_iRes[1]*ppot->m_iRes[0]+iy*ppot->m_iRes[0]+(ix+1)]
						+ ppot->m_pBin[iz*ppot->m_iRes[1]*ppot->m_iRes[0]+iy*ppot->m_iRes[0]+(ix-1)]
						+ ppot->m_pBin[iz*ppot->m_iRes[1]*ppot->m_iRes[0]+(iy+1)*ppot->m_iRes[0]+ix]
						+ ppot->m_pBin[iz*ppot->m_iRes[1]*ppot->m_iRes[0]+(iy-1)*ppot->m_iRes[0]+ix]
						+ ppot->m_pBin[(iz+1)*ppot->m_iRes[1]*ppot->m_iRes[0]+iy*ppot->m_iRes[0]+ix]
						+ ppot->m_pBin[(iz-1)*ppot->m_iRes[1]*ppot->m_iRes[0]+iy*ppot->m_iRes[0]+ix]
						- 6.0 * ppot->m_pBin[iz*ppot->m_iRes[1]*ppot->m_iRes[0]+iy*ppot->m_iRes[0]+ix];
					if (fabs(tf) > 10.0)
						mprintf("@ X=%d, Y=%d, Z=%d, tf=%E\n",ix,iy,iz,tf);
					pres->m_pBin[iz*ppot->m_iRes[1]*ppot->m_iRes[0]+iy*ppot->m_iRes[0]+ix] = tf;
					if (fabs(tf) > tfma)
						tfma = fabs(tf);
				}
			}
		}

		mprintf("      Max deviation: %20.10E\n",tfma);

		mprintf("      Correcting...\n");

		for (iz=0;iz<pob1->m_iRes[0]*pob1->m_iRes[1]*pob1->m_iRes[2];iz++)
			ppot->m_pBin[iz] = poball->m_pBin[iz] * ppotfix->m_pBin[iz] + (1.0 - poball->m_pBin[iz]) * (ppot->m_pBin[iz] + omega * pres->m_pBin[iz]);

		sprintf(buf,"poisson_pot_iter%04d.cube",iter+1);

		mprintf("      Writing current potential to \"%s\"...\n",buf);

		ppot->WriteCube(buf);

		iter++;
	}

	mprintf("*** Poisson Test Finished ***\n");*/
}



#endif



int cubetool_main( int argc, const char *argv[] ) {

	int z, z2, ti, ti2;
	int tix, tiy, tiz, tix2, tiy2, tiz2;
	CCube *tcu;
	double tf, tf2, tfx, tfy, tfz;


	mprintf("\n");
	mprintf(WHITE,"    #############################################\n");
	mprintf(WHITE,"    ####                                     ####\n");
	mprintf(WHITE,"    ####           *** CubeTool ***          ####\n");
	mprintf(WHITE,"    ####    (c) Martin Brehm, 2020 - 2022    ####\n");
	mprintf(WHITE,"    ####      https://brehm-research.de/     ####\n");
	mprintf(WHITE,"    ####                                     ####\n");
	mprintf(WHITE,"    #############################################\n");
	mprintf("\n");

	mprintf("    A tool to manipulate and process Gaussian Cube files.\n\n");

	if (argc < 2) {

		Usage();
		goto _ende;
	}


	g_pCubeWork = NULL;
	g_iCutNormal = -1;
	g_iCutSliceFrom = -1;
	g_iCutSliceTo = -1;
	g_pVectorField[0] = NULL;
	g_pVectorField[1] = NULL;
	g_pVectorField[2] = NULL;
	tcu = NULL;
	g_bShowAtoms = false;
	g_fAtomRadius = 1.0;
	g_fVectorScale = 1.0;
	g_fVectorHeadSize = 1.0;
	g_iVectorStride = 1;
	g_iPlaneX1 = -1;
	g_iPlaneY1 = -1;
	g_iPlaneX2 = -1;
	g_iPlaneY2 = -1;
	g_fValueRangeMin = 0;
	g_fValueRangeMax = 0;;
	g_fVectorRange = -1.0;
	g_bCustomValueRange = false;
	g_iGPInterpolation = 5;
	g_iContourLines = 30;
	g_fAspectRatio = 1.0;
	g_bShowTicks = true;
	g_bDrawMesh = true;
	g_iPlotPixel = 1024;
	g_fPlotExpo = 1.0;


	for (z=1;z<argc;z++) {

		mprintf(YELLOW,"\nNext command line argument: \"%s\"\n",argv[z]);

		/*****************************************************************************************************************/
		if (strcmp(argv[z],"-readcube") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			mprintf("    Trying to read cube file \"%s\" ...\n",argv[z+1]);
			g_pCubeWork = new CCube();
			if (!g_pCubeWork->ReadCube( argv[z+1] )) {
				eprintf("Reading failed.\n");
				goto _ende;
			}

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-testconv") == 0) {
		/*****************************************************************************************************************/

			TestConvolution();

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-testpoisson") == 0) {
		/*****************************************************************************************************************/

			#ifdef USE_CUBETOOL_FFTW
				TestPoisson();
			#else
				eprintf("Error: This requires to switch on USE_CUBETOOL_FFTW.\n");
				goto _ende;
			#endif

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-newcube") == 0) {
		/*****************************************************************************************************************/

			if (z+6 >= argc) {
				mprintf("    Error: This command requires six additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidUnsigned(argv[z+1])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			tix = atoi(argv[z+1]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidUnsigned(argv[z+2])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			tiy = atoi(argv[z+2]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+3]);
			if (!IsValidUnsigned(argv[z+3])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			tiz = atoi(argv[z+3]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+4]);
			if (!IsValidFloat(argv[z+4])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tfx = atof(argv[z+4]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+5]);
			if (!IsValidFloat(argv[z+5])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tfy = atof(argv[z+5]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+6]);
			if (!IsValidFloat(argv[z+5])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tfz = atof(argv[z+6]);

			mprintf("    Creating empty working buffer with resolution ( %d, %d, %d ) and extent ( %f, %f, %f ).\n",tix,tiy,tiz,tfx,tfy,tfz);

			if (g_pCubeWork != NULL)
				delete g_pCubeWork;

			g_pCubeWork = new CCube();

			g_pCubeWork->m_iRes[0] = tix;
			g_pCubeWork->m_iRes[1] = tiy;
			g_pCubeWork->m_iRes[2] = tiz;
			g_pCubeWork->m_fStrideA[0] = tfx/tix;
			g_pCubeWork->m_fStrideA[1] = 0;
			g_pCubeWork->m_fStrideA[2] = 0;
			g_pCubeWork->m_fStrideB[0] = 0;
			g_pCubeWork->m_fStrideB[1] = tfy/tiy;
			g_pCubeWork->m_fStrideB[2] = 0;
			g_pCubeWork->m_fStrideC[0] = 0;
			g_pCubeWork->m_fStrideC[1] = 0;
			g_pCubeWork->m_fStrideC[2] = tfz/tiz;
			g_pCubeWork->m_fCenter[0] = 0;
			g_pCubeWork->m_fCenter[1] = 0;
			g_pCubeWork->m_fCenter[2] = 0;
			g_pCubeWork->m_iAtomCount = 0;
			g_pCubeWork->m_pBin = new double[tix*tiy*tiz];
			for (ti=0;ti<tix*tiy*tiz;ti++)
				g_pCubeWork->m_pBin[z] = 0;

			z += 6;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-importbinary") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			mprintf("    Trying to read binary FP64 file \"%s\" ...\n",argv[z+1]);
			if (!g_pCubeWork->ImportBinary( argv[z+1] )) {
				eprintf("Import failed.\n");
				goto _ende;
			}

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-importbinary4d") == 0) {
		/*****************************************************************************************************************/

			if (z+3 >= argc) {
				eprintf("    Error: This command requires three additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidUnsigned(argv[z+2])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			ti = atoi(argv[z+2]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+3]);
			if (!IsValidUnsigned(argv[z+3])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			ti2 = atoi(argv[z+3]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			mprintf("    Trying to read frame %d/%d from binary FP64 file \"%s\" ...\n",ti,ti2,argv[z+1]);
			if (!g_pCubeWork->ImportBinary4D( argv[z+1], ti2, ti-1 )) {
				eprintf("Import failed.\n");
				goto _ende;
			}

			z += 3;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-analyzebinary4d") == 0) {
		/*****************************************************************************************************************/

			if (z+2 >= argc) {
				eprintf("    Error: This command requires two additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidUnsigned(argv[z+2])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			ti = atoi(argv[z+2]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			mprintf("    Trying to analyze %d frames from binary FP64 file \"%s\" ...\n",ti,argv[z+1]);
			if (!g_pCubeWork->AnalyzeBinary4D( argv[z+1], ti )) {
				eprintf("Analysis failed.\n");
				goto _ende;
			}

			z += 2;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-setoffset") == 0) {
		/*****************************************************************************************************************/

			if (z+3 >= argc) {
				eprintf("    Error: This command requires three additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tfx = atoi(argv[z+1]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidFloat(argv[z+2])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tfy = atoi(argv[z+2]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+3]);
			if (!IsValidFloat(argv[z+3])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tfz = atoi(argv[z+3]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			mprintf("    Setting Cube offset to ( %f | %f | %f ) ...\n",tfx,tfy,tfz);
			g_pCubeWork->m_fCenter[0] = tfx;
			g_pCubeWork->m_fCenter[1] = tfy;
			g_pCubeWork->m_fCenter[2] = tfz;

			z += 3;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-addcube") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				mprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}
			mprintf("    Trying to read cube file \"%s\" ...\n",argv[z+1]);
			tcu = new CCube();
			if (!tcu->ReadCube( argv[z+1] )) {
				eprintf("Reading failed.\n");
				goto _ende;
			}
			mprintf("    Adding volumetric data to the working buffer...\n");
			if (!g_pCubeWork->AddCube( tcu )) {
				eprintf("Failed.\n");
				goto _ende;
			}
			delete tcu;
			tcu = NULL;

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-subtractcube") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			mprintf("    Trying to read cube file \"%s\" ...\n",argv[z+1]);
			tcu = new CCube();
			if (!tcu->ReadCube( argv[z+1] )) {
				eprintf("Reading failed.\n");
				goto _ende;
			}
			mprintf("    Subtracting volumetric data from the working buffer...\n");
			if (!g_pCubeWork->SubtractCube( tcu )) {
				eprintf("Failed.\n");
				goto _ende;
			}
			delete tcu;
			tcu = NULL;

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-add") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tf = atof(argv[z+1]);

			mprintf("    Adding constant %f to working buffer...\n",tf);
			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}
			g_pCubeWork->Add( tf );

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-limitrange") == 0) {
		/*****************************************************************************************************************/

			if (z+2 >= argc) {
				eprintf("    Error: This command requires two additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tf = atof(argv[z+1]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidFloat(argv[z+2])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tf2 = atof(argv[z+2]);

			mprintf("    Limiting data range to [ %.6E, %.6E ]...\n",tf,tf2);
			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}
			g_pCubeWork->LimitRange( tf, tf2 );

			z += 2;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-multiply") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tf = atof(argv[z+1]);

			mprintf("    Multiplying working buffer with constant %f ...\n",tf);
			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}
			g_pCubeWork->Multiply( tf );

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-addlayer") == 0) {
		/*****************************************************************************************************************/

			if (z+4 >= argc) {
				eprintf("    Error: This command requires four additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidUnsigned(argv[z+1])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			tix = atoi(argv[z+1]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidUnsigned(argv[z+2])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			tiy = atoi(argv[z+2]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+3]);
			if (!IsValidUnsigned(argv[z+3])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			tiz = atoi(argv[z+3]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+4]);
			if (!IsValidFloat(argv[z+4])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tf = atof(argv[z+4]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			mprintf("    Adding layer of width (%d,%d,%d) and value %f ...\n",tix,tiy,tiz,tf);

			g_pCubeWork->AddLayer( tix, tiy, tiz, tf );

			z += 4;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-addgauss") == 0) {
		/*****************************************************************************************************************/

			if (z+5 >= argc) {
				eprintf("    Error: This command requires five additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tf = atof(argv[z+1]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidFloat(argv[z+2])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tfx = atof(argv[z+2]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+3]);
			if (!IsValidFloat(argv[z+3])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tfy = atof(argv[z+3]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+4]);
			if (!IsValidFloat(argv[z+4])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tfz = atof(argv[z+4]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+5]);
			if (!IsValidFloat(argv[z+5])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tf2 = atof(argv[z+5]);

			mprintf("    Adding Gaussian (height=%f, x=%f, y=%f, z=%f, s=%f) to working buffer ...\n",tf,tfx,tfy,tfz,tf2);
			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}
			g_pCubeWork->AddGaussian( tf, tfx, tfy, tfz, tf2 );

			z += 5;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-crop") == 0) {
		/*****************************************************************************************************************/

			if (z+6 >= argc) {
				eprintf("    Error: This command requires six additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidUnsigned(argv[z+1])) {
				eprintf("Error: Invalid argument, expected an integer number >= 0.\n");
				goto _ende;
			}
			tix = atoi(argv[z+1]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidUnsigned(argv[z+2])) {
				eprintf("Error: Invalid argument, expected an integer number >= 0.\n");
				goto _ende;
			}
			tiy = atoi(argv[z+2]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+3]);
			if (!IsValidUnsigned(argv[z+3])) {
				eprintf("Error: Invalid argument, expected an integer number >= 0.\n");
				goto _ende;
			}
			tiz = atoi(argv[z+3]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+4]);
			if (!IsValidUnsigned(argv[z+4])) {
				eprintf("Error: Invalid argument, expected an integer number >= 0.\n");
				goto _ende;
			}
			tix2 = atoi(argv[z+4]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+5]);
			if (!IsValidUnsigned(argv[z+5])) {
				eprintf("Error: Invalid argument, expected an integer number >= 0.\n");
				goto _ende;
			}
			tiy2 = atoi(argv[z+5]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+6]);
			if (!IsValidUnsigned(argv[z+6])) {
				eprintf("Error: Invalid argument, expected an integer number >= 0.\n");
				goto _ende;
			}
			tiz2 = atoi(argv[z+6]);

			if (g_pCubeWork != NULL) {
				mprintf("    Cropping working buffer to (%d,%d,%d) - (%d,%d,%d) ...\n",tix,tiy,tiz,tix2,tiy2,tiz2);
				tcu = new CCube();
				CropCube( g_pCubeWork, tcu, tix, tiy, tiz, tix2, tiy2, tiz2 );
				delete g_pCubeWork;
				g_pCubeWork = tcu;
			}

			if (g_pVectorField[0] != NULL) {
				mprintf("    Cropping vector field X component to (%d,%d,%d) - (%d,%d,%d) ...\n",tix,tiy,tiz,tix2,tiy2,tiz2);
				tcu = new CCube();
				CropCube( g_pVectorField[0], tcu, tix, tiy, tiz, tix2, tiy2, tiz2 );
				delete g_pVectorField[0];
				g_pVectorField[0] = tcu;
			}

			if (g_pVectorField[1] != NULL) {
				mprintf("    Cropping vector field Y component to (%d,%d,%d) - (%d,%d,%d) ...\n",tix,tiy,tiz,tix2,tiy2,tiz2);
				tcu = new CCube();
				CropCube( g_pVectorField[1], tcu, tix, tiy, tiz, tix2, tiy2, tiz2 );
				delete g_pVectorField[1];
				g_pVectorField[1] = tcu;
			}

			if (g_pVectorField[2] != NULL) {
				mprintf("    Cropping vector field Z component to (%d,%d,%d) - (%d,%d,%d) ...\n",tix,tiy,tiz,tix2,tiy2,tiz2);
				tcu = new CCube();
				CropCube( g_pVectorField[2], tcu, tix, tiy, tiz, tix2, tiy2, tiz2 );
				delete g_pVectorField[2];
				g_pVectorField[2] = tcu;
			}

			z += 6;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-envelope") == 0) {
		/*****************************************************************************************************************/

			if (z+2 >= argc) {
				eprintf("    Error: This command requires two additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			if (mystricmp(argv[z+1],"cossq") != 0) {
				eprintf("Only \"cossq\" type of envelope currently implemented.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidFloat(argv[z+2])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tf = atof(argv[z+2]);

			if ((tf < 0) || (tf >= 1.0)) {
				eprintf("Error: Inner point fraction (%f) outside of valid range (0 <= i < 1.0).\n",tf);
				goto _ende;
			}

			if (g_pCubeWork != NULL) {
				mprintf("      Applying Cos^2 envelope with fraction %f of inner points to working buffer...\n",tf);
				EnvelopeCosSq( g_pCubeWork, tf );
			}

			if (g_pVectorField[0] != NULL) {
				mprintf("      Applying Cos^2 envelope with fraction %f of inner points to X component of vector field...\n",tf);
				EnvelopeCosSq( g_pVectorField[0], tf );
			}

			if (g_pVectorField[1] != NULL) {
				mprintf("      Applying Cos^2 envelope with fraction %f of inner points to Y component of vector field...\n",tf);
				EnvelopeCosSq( g_pVectorField[1], tf );
			}

			if (g_pVectorField[2] != NULL) {
				mprintf("      Applying Cos^2 envelope with fraction %f of inner points to Z component of vector field...\n",tf);
				EnvelopeCosSq( g_pVectorField[2], tf );
			}

			z += 2;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-helmholtz_decomp") == 0) {
		/*****************************************************************************************************************/

			if (z+2 >= argc) {
				eprintf("    Error: This command requires two additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidUnsigned(argv[z+1])) {
				eprintf("Error: Invalid argument, expected an unsigned integer.\n");
				goto _ende;
			}
			ti = atoi(argv[z+1]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidFloat(argv[z+2])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tf = atof(argv[z+2]);

			if (g_pVectorField[0] == NULL) {
				eprintf("Error: Component X of vector field is empty.\n");
				goto _ende;
			}

			if (g_pVectorField[1] == NULL) {
				eprintf("Error: Component Y of vector field is empty.\n");
				goto _ende;
			}

			if (g_pVectorField[2] == NULL) {
				eprintf("Error: Component Z of vector field is empty.\n");
				goto _ende;
			}

			Helmholtz_Decomp( ti, tf );

			z += 2;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-helmholtz_old") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			tf = atof(argv[z+1]);

			if (g_pVectorField[0] == NULL) {
				eprintf("Error: Component X of vector field is empty.\n");
				goto _ende;
			}

			if (g_pVectorField[1] == NULL) {
				eprintf("Error: Component Y of vector field is empty.\n");
				goto _ende;
			}

			if (g_pVectorField[2] == NULL) {
				eprintf("Error: Component Z of vector field is empty.\n");
				goto _ende;
			}

			Helmholtz_Decomp_Old( tf );

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-writecube") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			mprintf("    Writing working buffer to cube file \"%s\" ...\n",argv[z+1]);
			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}
			g_pCubeWork->WriteCube( argv[z+1] );

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-vectorfield") == 0) {
		/*****************************************************************************************************************/

			if (z+3 >= argc) {
				eprintf("    Error: This command requires three additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			mprintf("    Consuming next argument: \"%s\"\n",argv[z+3]);

			mprintf("    Trying to read X component of vector field from cube file \"%s\" ...\n",argv[z+1]);
			g_pVectorField[0] = new CCube();
			if (!g_pVectorField[0]->ReadCube( argv[z+1] )) {
				eprintf("Reading failed.\n");
				goto _ende;
			}

			mprintf("    Trying to read Y component of vector field from cube file \"%s\" ...\n",argv[z+2]);
			g_pVectorField[1] = new CCube();
			if (!g_pVectorField[1]->ReadCube( argv[z+2] )) {
				eprintf("Reading failed.\n");
				goto _ende;
			}

			mprintf("    Trying to read Z component of vector field from cube file \"%s\" ...\n",argv[z+3]);
			g_pVectorField[2] = new CCube();
			if (!g_pVectorField[2]->ReadCube( argv[z+3] )) {
				eprintf("Reading failed.\n");
				goto _ende;
			}

			tf = 0;
			for (z2=0;z2<g_pVectorField[0]->m_iRes[0]*g_pVectorField[0]->m_iRes[1]*g_pVectorField[0]->m_iRes[2];z2++) {
				tf2 = SQR(g_pVectorField[0]->m_pBin[z2]) + SQR(g_pVectorField[1]->m_pBin[z2]) + SQR(g_pVectorField[2]->m_pBin[z2]);
				if (tf2 > tf)
					tf = tf2;
			}
			tf = sqrt(tf);

			mprintf("    The maximum vector length in the vector field is %.6f\n",tf);

			z += 3;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-pushtovectorfield") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (mystricmp(argv[z+1],"x") == 0) {
				mprintf("      Copying working buffer to X component of vector field.\n");
				if (g_pVectorField[0] != NULL)
					delete g_pVectorField[0];
				g_pVectorField[0] = new CCube();
				g_pVectorField[0]->CopyFrom( g_pCubeWork );
			} else if (mystricmp(argv[z+1],"y") == 0) {
				mprintf("      Copying working buffer to Y component of vector field.\n");
				if (g_pVectorField[1] != NULL)
					delete g_pVectorField[1];
				g_pVectorField[1] = new CCube();
				g_pVectorField[1]->CopyFrom( g_pCubeWork );
			} else if (mystricmp(argv[z+1],"z") == 0) {
				mprintf("      Copying working buffer to Z component of vector field.\n");
				if (g_pVectorField[2] != NULL)
					delete g_pVectorField[2];
				g_pVectorField[2] = new CCube();
				g_pVectorField[2]->CopyFrom( g_pCubeWork );
			} else {
				eprintf("Invalid argument: \"%s\". Allowed are \"x\", \"y\", and \"z\".\n",argv[z+1]);
				goto _ende;
			}

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-vectorfromgradient") == 0) {
		/*****************************************************************************************************************/

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			VectorFieldFromGradient( g_pCubeWork );

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-pullfromvectorfield") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			if (mystricmp(argv[z+1],"x") == 0) {
				mprintf("      Copying X component of vector field to working buffer.\n");
				if (g_pVectorField[0] == NULL) {
					eprintf("Error: The X component of the vector field is empty.\n");
					goto _ende;
				}
				if (g_pCubeWork != NULL)
					delete g_pCubeWork;
				g_pCubeWork = new CCube();
				g_pCubeWork->CopyFrom( g_pVectorField[0] );
			} else if (mystricmp(argv[z+1],"y") == 0) {
				mprintf("      Copying Y component of vector field to working buffer.\n");
				if (g_pVectorField[1] == NULL) {
					eprintf("Error: The Y component of the vector field is empty.\n");
					goto _ende;
				}
				if (g_pCubeWork != NULL)
					delete g_pCubeWork;
				g_pCubeWork = new CCube();
				g_pCubeWork->CopyFrom( g_pVectorField[1] );
			} else if (mystricmp(argv[z+1],"z") == 0) {
				mprintf("      Copying Z component of vector field to working buffer.\n");
				if (g_pVectorField[2] == NULL) {
					eprintf("Error: The Z component of the vector field is empty.\n");
					goto _ende;
				}
				if (g_pCubeWork != NULL)
					delete g_pCubeWork;
				g_pCubeWork = new CCube();
				g_pCubeWork->CopyFrom( g_pVectorField[2] );
			} else {
				eprintf("Invalid argument: \"%s\". Allowed are \"x\", \"y\", and \"z\".\n",argv[z+1]);
				goto _ende;
			}

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-computecurl") == 0) {
		/*****************************************************************************************************************/

			if (g_pVectorField[0] == NULL) {
				eprintf("Error: No vector field loaded.\n");
				goto _ende;
			}
			mprintf("    Computing curl of vector field...\n");
			ComputeCurl();

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-computedivergence") == 0) {
		/*****************************************************************************************************************/

			if (g_pVectorField[0] == NULL) {
				eprintf("Error: No vector field loaded.\n");
				goto _ende;
			}
			mprintf("    Computing divergence of vector field...\n");
			ComputeDivergence();

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-vectorabs") == 0) {
		/*****************************************************************************************************************/

			if (g_pVectorField[0] == NULL) {
				eprintf("Error: No vector field loaded.\n");
				goto _ende;
			}
			mprintf("    Creating working buffer from absolute value of vector field...\n");
			CreateAbs();

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-normalx") == 0) {
		/*****************************************************************************************************************/

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}
			mprintf("    Defining cut plane normal to X direction.\n");
			g_iCutNormal = 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-normaly") == 0) {
		/*****************************************************************************************************************/

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}
			mprintf("    Defining cut plane normal to Y direction.\n");
			g_iCutNormal = 2;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-normalz") == 0) {
		/*****************************************************************************************************************/

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}
			mprintf("    Defining cut plane normal to Z direction.\n");
			g_iCutNormal = 3;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-planesection") == 0) {
		/*****************************************************************************************************************/

			if (z+4 >= argc) {
				eprintf("    Error: This command requires four additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidUnsigned(argv[z+1])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidUnsigned(argv[z+2])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			mprintf("    Consuming next argument: \"%s\"\n",argv[z+3]);
			if (!IsValidUnsigned(argv[z+3])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			mprintf("    Consuming next argument: \"%s\"\n",argv[z+4]);
			if (!IsValidUnsigned(argv[z+4])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}

			ti = atoi(argv[z+1]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			g_iPlaneX1 = atoi(argv[z+1]);
			g_iPlaneX2 = atoi(argv[z+2]);
			g_iPlaneY1 = atoi(argv[z+3]);
			g_iPlaneY2 = atoi(argv[z+4]);

			switch(g_iCutNormal) {

				case 1: // X
					if ((g_iPlaneX1	< 0) || (g_iPlaneX1 >= g_pCubeWork->m_iRes[1])) {
						eprintf("Error: X1 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[1],g_iPlaneX1);
						goto _ende;
					}
					if ((g_iPlaneY1	< 0) || (g_iPlaneY1 >= g_pCubeWork->m_iRes[2])) {
						eprintf("Error: Y1 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[2],g_iPlaneY1);
						goto _ende;
					}
					if ((g_iPlaneX2	< 0) || (g_iPlaneX2 >= g_pCubeWork->m_iRes[1])) {
						eprintf("Error: X2 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[1],g_iPlaneX2);
						goto _ende;
					}
					if ((g_iPlaneY2	< 0) || (g_iPlaneY2 >= g_pCubeWork->m_iRes[2])) {
						eprintf("Error: Y2 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[2],g_iPlaneY2);
						goto _ende;
					}
					if (g_iPlaneX1 >= g_iPlaneX2) {
						eprintf("Error: X1 >= X2.\n");
						goto _ende;
					}
					if (g_iPlaneY1 >= g_iPlaneY2) {
						eprintf("Error: Y1 >= Y2.\n");
						goto _ende;
					}
					break;

				case 2: // Y
					if ((g_iPlaneX1	< 0) || (g_iPlaneX1 >= g_pCubeWork->m_iRes[0])) {
						eprintf("Error: X1 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[0],g_iPlaneX1);
						goto _ende;
					}
					if ((g_iPlaneY1	< 0) || (g_iPlaneY1 >= g_pCubeWork->m_iRes[2])) {
						eprintf("Error: Y1 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[2],g_iPlaneY1);
						goto _ende;
					}
					if ((g_iPlaneX2	< 0) || (g_iPlaneX2 >= g_pCubeWork->m_iRes[0])) {
						eprintf("Error: X2 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[0],g_iPlaneX2);
						goto _ende;
					}
					if ((g_iPlaneY2	< 0) || (g_iPlaneY2 >= g_pCubeWork->m_iRes[2])) {
						eprintf("Error: Y2 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[2],g_iPlaneY2);
						goto _ende;
					}
					if (g_iPlaneX1 >= g_iPlaneX2) {
						eprintf("Error: X1 >= X2.\n");
						goto _ende;
					}
					if (g_iPlaneY1 >= g_iPlaneY2) {
						eprintf("Error: Y1 >= Y2.\n");
						goto _ende;
					}
					break;

				case 3: // Z
					if ((g_iPlaneX1	< 0) || (g_iPlaneX1 >= g_pCubeWork->m_iRes[0])) {
						eprintf("Error: X1 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[0],g_iPlaneX1);
						goto _ende;
					}
					if ((g_iPlaneY1	< 0) || (g_iPlaneY1 >= g_pCubeWork->m_iRes[1])) {
						eprintf("Error: Y1 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[1],g_iPlaneY1);
						goto _ende;
					}
					if ((g_iPlaneX2	< 0) || (g_iPlaneX2 >= g_pCubeWork->m_iRes[0])) {
						eprintf("Error: X2 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[0],g_iPlaneX2);
						goto _ende;
					}
					if ((g_iPlaneY2	< 0) || (g_iPlaneY2 >= g_pCubeWork->m_iRes[1])) {
						eprintf("Error: Y2 argument out of range (allowed 0 .. %d, found %d).\n",g_pCubeWork->m_iRes[1],g_iPlaneY2);
						goto _ende;
					}
					if (g_iPlaneX1 >= g_iPlaneX2) {
						eprintf("Error: X1 >= X2.\n");
						goto _ende;
					}
					if (g_iPlaneY1 >= g_iPlaneY2) {
						eprintf("Error: Y1 >= Y2.\n");
						goto _ende;
					}
					break;
			}

			printf("    Using { %d .. %d } x { %d .. %d } as the section of the cut plane.\n",g_iPlaneX1,g_iPlaneX2,g_iPlaneY1,g_iPlaneY2);

			z += 4;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-slice") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidUnsigned(argv[z+1])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			ti = atoi(argv[z+1]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			switch(g_iCutNormal) {
				case 1: // X
					if ((ti < 0) || (ti >= g_pCubeWork->m_iRes[0])) {
						eprintf("Error: Slice index out of range (found %d, allowed 0 .. %d).\n",ti,g_pCubeWork->m_iRes[0]-1);
						goto _ende;
					}
					mprintf("    Taking slice %d from range 0 .. %d in X direction.\n",ti,g_pCubeWork->m_iRes[0]-1);
					g_iCutSliceFrom = ti;
					g_iCutSliceTo = ti;
					break;
				case 2: // Y
					if ((ti < 0) || (ti >= g_pCubeWork->m_iRes[1])) {
						eprintf("Error: Slice index out of range (found %d, allowed 0 .. %d).\n",ti,g_pCubeWork->m_iRes[1]-1);
						goto _ende;
					}
					mprintf("    Taking slice %d from range 0 .. %d in Y direction.\n",ti,g_pCubeWork->m_iRes[1]-1);
					g_iCutSliceFrom = ti;
					g_iCutSliceTo = ti;
					break;
				case 3: // Z
					if ((ti < 0) || (ti >= g_pCubeWork->m_iRes[2])) {
						eprintf("Error: Slice index out of range (found %d, allowed 0 .. %d).\n",ti,g_pCubeWork->m_iRes[2]-1);
						goto _ende;
					}
					mprintf("    Taking slice %d from range 0 .. %d in Z direction.\n",ti,g_pCubeWork->m_iRes[2]-1);
					g_iCutSliceFrom = ti;
					g_iCutSliceTo = ti;
					break;
			}

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-slicerange") == 0) {
		/*****************************************************************************************************************/

			if (z+2 >= argc) {
				eprintf("    Error: This command requires two additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidUnsigned(argv[z+1])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			ti = atoi(argv[z+1]);

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidUnsigned(argv[z+2])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}
			ti2 = atoi(argv[z+2]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			if (ti > ti2) {
				eprintf("Error: First index needs to be smaller than second index.\n");
				goto _ende;
			}

			switch(g_iCutNormal) {
				case 1: // X
					if ((ti < 0) || (ti >= g_pCubeWork->m_iRes[0])) {
						eprintf("Error: Slice index 1 out of range (found %d, allowed 0 .. %d).\n",ti,g_pCubeWork->m_iRes[0]-1);
						goto _ende;
					}
					if ((ti2 < 0) || (ti2 >= g_pCubeWork->m_iRes[0])) {
						eprintf("Error: Slice index 2 out of range (found %d, allowed 0 .. %d).\n",ti2,g_pCubeWork->m_iRes[0]-1);
						goto _ende;
					}
					mprintf("    Taking slice range %d .. %d from total range 0 .. %d in X direction.\n",ti,ti2,g_pCubeWork->m_iRes[0]-1);
					g_iCutSliceFrom = ti;
					g_iCutSliceTo = ti2;
					break;
				case 2: // Y
					if ((ti < 0) || (ti >= g_pCubeWork->m_iRes[1])) {
						eprintf("Error: Slice index 1 out of range (found %d, allowed 0 .. %d).\n",ti,g_pCubeWork->m_iRes[1]-1);
						goto _ende;
					}
					if ((ti2 < 0) || (ti2 >= g_pCubeWork->m_iRes[1])) {
						eprintf("Error: Slice index 2 out of range (found %d, allowed 0 .. %d).\n",ti2,g_pCubeWork->m_iRes[1]-1);
						goto _ende;
					}
					mprintf("    Taking slice range %d .. %d from total range 0 .. %d in Y direction.\n",ti,ti2,g_pCubeWork->m_iRes[1]-1);
					g_iCutSliceFrom = ti;
					g_iCutSliceTo = ti2;
					break;
				case 3: // Z
					if ((ti < 0) || (ti >= g_pCubeWork->m_iRes[2])) {
						eprintf("Error: Slice index 1 out of range (found %d, allowed 0 .. %d).\n",ti,g_pCubeWork->m_iRes[2]-1);
						goto _ende;
					}
					if ((ti2 < 0) || (ti2 >= g_pCubeWork->m_iRes[2])) {
						eprintf("Error: Slice index 2 out of range (found %d, allowed 0 .. %d).\n",ti2,g_pCubeWork->m_iRes[2]-1);
						goto _ende;
					}
					mprintf("    Taking slice range %d .. %d from total range 0 .. %d in Z direction.\n",ti,ti2,g_pCubeWork->m_iRes[2]-1);
					g_iCutSliceFrom = ti;
					g_iCutSliceTo = ti2;
					break;
			}

			z += 2;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-nogrid") == 0) {
		/*****************************************************************************************************************/

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			g_bDrawMesh = false;
			mprintf("    Disabling the grid in the cut plane plot.\n");

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-noticks") == 0) {
		/*****************************************************************************************************************/

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			g_bShowTicks = false;
			mprintf("    Disabling the axis ticks and tick labels in the cut plane plot.\n");

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-nointerpolate") == 0) {
		/*****************************************************************************************************************/

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			g_iGPInterpolation = 1;
			mprintf("    Disabling the interpolation between the data points in the cut plane plot.\n");

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-ncontour") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidUnsigned(argv[z+1])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			g_iContourLines = atoi(argv[z+1]);
			if (g_iContourLines > 0)
				mprintf("    Drawing %d contour lines in the cut plane plot.\n",g_iContourLines);
			else
				mprintf("    Disabling contour lines in the cut plane plot.\n");

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-aspectratio") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			g_fAspectRatio = atof(argv[z+1]);
			mprintf("    Setting the aspect ratio of the cut plane plot to %f.\n",g_fAspectRatio);

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-plotexpo") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			g_fPlotExpo = atof(argv[z+1]);
			mprintf("    Setting the plot exponent of the cut plane plot to %f.\n",g_fPlotExpo);

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-plotpixel") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidUnsigned(argv[z+1])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			g_iPlotPixel = atoi(argv[z+1]);
			mprintf("    Setting the width of the cut plane plot to %d pixel.\n",g_iPlotPixel);

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-showatoms") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			g_fAtomRadius = atof(argv[z+1]);
			g_bShowAtoms = true;
			mprintf("    Displaying the atoms in the cut plane plot with radius %.3f.\n",g_fAtomRadius);

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-valuerange") == 0) {
		/*****************************************************************************************************************/

			if (z+2 >= argc) {
				eprintf("    Error: This command requires two additional arguments.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}
			mprintf("    Consuming next argument: \"%s\"\n",argv[z+2]);
			if (!IsValidFloat(argv[z+2])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			g_bCustomValueRange = true;
			g_fValueRangeMin = atof(argv[z+1]);
			g_fValueRangeMax = atof(argv[z+2]);
			mprintf("    Using a value range of %.5f .. %.5f for the cut plane plot.\n",g_fValueRangeMin,g_fValueRangeMax);

			z += 2;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-vectorrange") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			if (g_pVectorField[0] == NULL) {
				eprintf("Error: No vector field loaded.\n");
				goto _ende;
			}

			g_fVectorRange = atof(argv[z+1]);
			mprintf("    Using a max. vector length of %.5f for the cut plane plot.\n",g_fVectorRange);

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-vectorscale") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			if (g_pVectorField[0] == NULL) {
				eprintf("Error: No vector field loaded.\n");
				goto _ende;
			}

			g_fVectorScale = atof(argv[z+1]);
			mprintf("    Scaling the vector field in the cut plane with factor %.3f.\n",g_fVectorScale);

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-vectorhead") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidFloat(argv[z+1])) {
				eprintf("Error: Invalid argument, expected a real number.\n");
				goto _ende;
			}

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			if (g_pVectorField[0] == NULL) {
				eprintf("Error: No vector field loaded.\n");
				goto _ende;
			}

			g_fVectorHeadSize = atof(argv[z+1]);
			mprintf("    Scaling the vector arrow size in the cut plane with factor %.3f.\n",g_fVectorHeadSize);

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-vectorstride") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);
			if (!IsValidUnsigned(argv[z+1])) {
				eprintf("Error: Invalid argument, expected unsigned integer.\n");
				goto _ende;
			}

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			if (g_pVectorField[0] == NULL) {
				eprintf("Error: No vector field loaded.\n");
				goto _ende;
			}

			g_iVectorStride = atoi(argv[z+1]);
			mprintf("    Using a stride of %d to show the vector field in the cut plane.\n",g_iVectorStride);

			z += 1;

		/*****************************************************************************************************************/
		} else if (strcmp(argv[z],"-writeplane") == 0) {
		/*****************************************************************************************************************/

			if (z+1 >= argc) {
				eprintf("    Error: This command requires an additional argument.\n");
				goto _ende;
			}

			mprintf("    Consuming next argument: \"%s\"\n",argv[z+1]);

			if (g_pCubeWork == NULL) {
				eprintf("Error: The working buffer is empty.\n");
				goto _ende;
			}

			if (g_iCutNormal == -1) {
				eprintf("Error: No cut plane defined.\n");
				goto _ende;
			}

			if (g_iCutSliceFrom == -1) {
				eprintf("Error: No slice / slice range defined.\n");
				goto _ende;
			}

			mprintf("    Writing cut plane to \"%s\" ...\n",argv[z+1]);

			if (!WriteSlice( g_pCubeWork, g_iCutNormal, g_iCutSliceFrom, g_iCutSliceTo, argv[z+1] )) {
				eprintf("Failed.\n");
				goto _ende;
			}

			z += 1;

		/*****************************************************************************************************************/
		} else {
		/*****************************************************************************************************************/

			eprintf("Error: Unrecognized command line argument: \"%s\".\n",argv[z]);
			Usage();
			goto _ende;
		}
		/*****************************************************************************************************************/
	}

	mprintf("Finished all tasks.\n");

	if (g_pCubeWork != NULL) {
		delete g_pCubeWork;
		g_pCubeWork = NULL;
	}

_ende:
	mprintf(WHITE,"\nCubeTool leaving.\n\n");

	return 0;
}


