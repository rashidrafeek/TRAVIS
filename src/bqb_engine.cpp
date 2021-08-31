/***********************************************************************************

    LibBQB - File Format and Compression Algorithms for Trajectories of
             Volumetric Data and Atom Positions

    https://brehm-research.de/bqb

    Free software, licensed under GNU LGPL v3

    Copyright (c) Martin Brehm and Martin Thomas,
                  Martin Luther University Halle-Wittenberg, Germany,
                  2016 - 2021.

    Please cite:  M. Brehm, M. Thomas: "An Efficient Lossless Compression Algorithm
                  for Trajectories of Atom Positions and Volumetric Data",
                  J. Chem. Inf. Model. 2018, 58 (10), pp 2092-2107.

    --------------------------------------------------------------------------------

    LibBQB is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

***********************************************************************************/



// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"
#include <math.h>
#include <algorithm>
#include <time.h>
#include "bqb_engine.h"
#include "bqb_hilbert.h"
#include "bqb_integerengine.h"
#include "bqb_bitset.h"
#include "bqb_crc.h"
#include "bqb_format.h"
#include "bqb_math.h"


const char *GetRevisionInfo_bqb_engine(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_engine() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



#define AR_EPS (0.0001)
#define CR_EPS (0.0001)




CBQBEngine::CBQBEngine(CBQBInterface &i) :
 		m_pExtrapolator(NULL),
		m_pExtrapolatorCorr(NULL),
		m_pExtrapolatorXYZ(NULL),
		m_pReadCacheXYZ(NULL),
		m_pReadCacheCube(NULL),
		m_bDeleteOutFrames(true),
		m_IF(i)
	{
}


CBQBEngine::~CBQBEngine() {

	int z;

	if (m_bDeleteOutFrames) {

		for (z=0;z<(int)m_oaOutputAtomBuf.size();z++) {
			if (m_oaOutputAtomBuf[z] != NULL) {
				delete m_oaOutputAtomBuf[z];
				m_oaOutputAtomBuf[z] = NULL;
			}
		}
		for (z=0;z<(int)m_oaOutputCubeBuf.size();z++) {
			if (m_oaOutputCubeBuf[z] != NULL) {
				delete m_oaOutputCubeBuf[z];
				m_oaOutputCubeBuf[z] = NULL;
			}
		}
		for (z=0;z<(int)m_oaOutput2AtomBuf.size();z++) {
			if (m_oaOutput2AtomBuf[z] != NULL) {
				delete m_oaOutput2AtomBuf[z];
				m_oaOutput2AtomBuf[z] = NULL;
			}
		}
		for (z=0;z<(int)m_oaOutput2CubeBuf.size();z++) {
			if (m_oaOutput2CubeBuf[z] != NULL) {
				delete m_oaOutput2CubeBuf[z];
				m_oaOutput2CubeBuf[z] = NULL;
			}
		}
	}

	if (m_pExtrapolator != NULL) {
		delete m_pExtrapolator;
		m_pExtrapolator = NULL;
	}
	if (m_pExtrapolatorCorr != NULL) {
		delete m_pExtrapolatorCorr;
		m_pExtrapolatorCorr = NULL;
	}
	if (m_pExtrapolatorXYZ != NULL) {
		delete m_pExtrapolatorXYZ;
		m_pExtrapolatorXYZ = NULL;
	}
	if (m_pReadCacheXYZ != NULL) {
		delete m_pReadCacheXYZ;
		m_pReadCacheXYZ = NULL;
	}
	if (m_pReadCacheCube != NULL) {
		delete m_pReadCacheCube;
		m_pReadCacheCube = NULL;
	}
}


void CBQBEngine::Reset() {

	int z;

	if (!m_bDeleteOutFrames) {

		for (z=0;z<(int)m_oaOutputAtomBuf.size();z++)
			m_oaOutputAtomBuf[z] = NULL;
		for (z=0;z<(int)m_oaOutputCubeBuf.size();z++)
			m_oaOutputCubeBuf[z] = NULL;
		for (z=0;z<(int)m_oaOutput2AtomBuf.size();z++)
			m_oaOutput2AtomBuf[z] = NULL;
		for (z=0;z<(int)m_oaOutput2CubeBuf.size();z++)
			m_oaOutput2CubeBuf[z] = NULL;

	} else {

		for (z=0;z<(int)m_oaOutputAtomBuf.size();z++) {
			if (m_oaOutputAtomBuf[z] != NULL) {
				delete m_oaOutputAtomBuf[z];
				m_oaOutputAtomBuf[z] = NULL;
			}
		}
		for (z=0;z<(int)m_oaOutputCubeBuf.size();z++) {
			if (m_oaOutputCubeBuf[z] != NULL) {
				delete m_oaOutputCubeBuf[z];
				m_oaOutputCubeBuf[z] = NULL;
			}
		}
		for (z=0;z<(int)m_oaOutput2AtomBuf.size();z++) {
			if (m_oaOutput2AtomBuf[z] != NULL) {
				delete m_oaOutput2AtomBuf[z];
				m_oaOutput2AtomBuf[z] = NULL;
			}
		}
		for (z=0;z<(int)m_oaOutput2CubeBuf.size();z++) {
			if (m_oaOutput2CubeBuf[z] != NULL) {
				delete m_oaOutput2CubeBuf[z];
				m_oaOutput2CubeBuf[z] = NULL;
			}
		}
	}

	if (m_pExtrapolator != NULL)
		m_pExtrapolator->Reset();

	if (m_pExtrapolatorCorr != NULL)
		m_pExtrapolatorCorr->Reset();

	if (m_pExtrapolatorXYZ != NULL)
		m_pExtrapolatorXYZ->Reset();

	if (m_pReadCacheXYZ != NULL) {
		delete m_pReadCacheXYZ;
		m_pReadCacheXYZ = NULL;
	}
	if (m_pReadCacheCube != NULL) {
		delete m_pReadCacheCube;
		m_pReadCacheCube = NULL;
	}
}


int CompareHilbertPair( const void *arg1, const void *arg2 ) {
	
	if (((THilbertPair*)arg1)->m_iH > ((THilbertPair*)arg2)->m_iH)
		return 1;
	else if (((THilbertPair*)arg1)->m_iH < ((THilbertPair*)arg2)->m_iH)
		return -1;
	else
		return 0;
}


void CBQBEngine::BuildHilbertIdx(int resx, int resy, int resz) {

	int b, bt;
	int resxyz, resyz;
	bitmask_t mi;
	unsigned long crd[3], /*crd2[3],*/ tul;
	FILE *a;
	unsigned char *uc;
	THilbertPair *hlist;
	CBQBTools bqbtools(m_IF);


	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("\n");
		m_IF.printf("Creating Hilbert curve index table...\n");
	}
	else if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("      Creating Hilbert curve index table (%dx%dx%d)...\n",resx,resy,resz);

	b = (int)(ceil(mylog2(resx))+0.5);
	bt = (int)(ceil(mylog2(resy))+0.5);
	if (bt > b)
		b = bt;
	bt = (int)(ceil(mylog2(resz))+0.5);
	if (bt > b)
		b = bt;
	mi = ((unsigned long)2)<<(3*b-1);

	if (b > 10) {
		m_IF.eprintf("CBQBEngine::BuildHilbertIdx(): Error: Cube file resolution is too large (%d x %d x %d).\n",
			resx,resy,resz);
		m_IF.printf("\n");
		m_IF.printf("Resolutions >= 1024 are not supported by Hilbert curve indexing.\n");
		m_IF.printf("You can disable Hilbert curve indexing by specifying \"-chilbert no\".\n");
		abort();
	}

	resxyz = resx*resy*resz;
	resyz = resy*resz;

	m_iaHilbertIdx.resize(resxyz);

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("Resolution %d | %d | %d --> Bits %d\n",resx,resy,resz,b);
		m_IF.printf("Requires %lu indices, %.3f MiB of RAM.\n",mi,resxyz*sizeof(int)/1024.0/1024.0);
		m_IF.printf("Running...\n");
	}
	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.printf("        [");
		m_IF.FlushLog();
	}

	uc = new unsigned char[resxyz];
	memset(uc,0,resxyz);

/*
	int ti;
	bitmask_t index, coords[3];
  
	ti = 0;
	for (index=0;index<mi;index++) {
		if ((index%100) == 0)
			if (fmod(index/100,mi/6000.0) < 1.0) {
				m_IF.printf("#");
				fflush(stdout);
			}
		hilbert_i2c(3,b,index,coords);
		if (((long)coords[2] >= resx) || ((long)coords[1] >= resy) || ((long)coords[0] >= resz))
			continue;
		m_iaHilbertIdx[ti] = coords[2]*resyz+coords[1]*resz+coords[0];
		uc[m_iaHilbertIdx[ti]]++;
		ti++;
	}*/

	hlist = new THilbertPair[resxyz];

	for (crd[0]=0;(long)crd[0]<resz;crd[0]++) {
		for (crd[1]=0;(long)crd[1]<resy;crd[1]++) {
			if (fmod(crd[0]*resy+crd[1],resz*resy/60.0) < 1.0) {
				if (m_IF.IsPL(BQB_PL_STANDARD)) {
					m_IF.printf("#");
					m_IF.FlushLog();
				}
			}
			for (crd[2]=0;(long)crd[2]<resx;crd[2]++) {
				tul = crd[0]+crd[1]*resz+crd[2]*resyz;
				hlist[tul].m_iN = tul;
				hlist[tul].m_iH = hilbert_c2i(3,b,crd);
	//			hilbert_i2c(3,b,hlist[tul].m_iH,crd2);
	//			if ((crd[0] != crd2[0]) || (crd[1] != crd2[1]) || (crd[2] != crd2[2]))
	//				m_IF.eprintf("Mismatch.\n");
			}
		}
	}

	qsort(hlist,resxyz,sizeof(THilbertPair),CompareHilbertPair);

	for (b=0;b<resxyz;b++) {
		m_iaHilbertIdx[b] = hlist[b].m_iN;
		uc[m_iaHilbertIdx[b]]++;
	}

	delete[] hlist;

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("]\n");

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("Sanity check...\n");

	for (b=0;b<resxyz;b++) {
		if (uc[b] > 1) {
			m_IF.eprintf("CBQBEngine::BuildHilbertIdx(): Error: Element %d touched multiple times (%u).\n",
				b,uc[b]);
			m_IF.eprintf("Hilbert curve index creation failed.\n");
			m_IF.printf("\n");
			m_IF.printf("You can disable Hilbert curve indexing by specifying \"-chilbert no\".\n");
			abort();
		}
		if (uc[b] == 0) {
			m_IF.eprintf("CBQBEngine::BuildHilbertIdx(): Error: Element %d never touched.\n",b);
			m_IF.eprintf("Hilbert curve index creation failed.\n");
			m_IF.printf("\n");
			m_IF.printf("You can disable Hilbert curve indexing by specifying \"-chilbert no\".\n");
			abort();
		}
	}
	delete[] uc;

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("Finished. %lu entries.\n",(unsigned long)m_iaHilbertIdx.size());
		m_IF.printf("Dumping Hilbert index to file \"hilbert.csv\"...\n");
		a = bqbtools.BQBOpenFileWrite("hilbert.csv",true);
		m_IF.printf("  [");
		for (b=0;b<resxyz;b++) {
			if ((b%100) == 0)
				if (fmod(b/100,resxyz/6000.0) < 1.0) {
					m_IF.printf("#");
					m_IF.FlushLog();
				}
			fprintf(a,"%d;  %d\n",b,m_iaHilbertIdx[b]);
		}
		m_IF.printf("]\n");
		m_IF.printf("\n");
		fclose(a);
	}
}


void CBQBEngine::ExportCellInfo(CBQBOldCellInfo *ci, CBQBBitSet *bs) {

	if (ci->IsCubic()) {

		bs->WriteBit(1);
		bs->WriteBits( (unsigned long)(ci->m_mCell(0,0)*10000.0+0.5), 32 );

	} else {
		
		bs->WriteBit(0);

		if (ci->IsOrthorhombic()) {

			bs->WriteBit(1);
			bs->WriteBits( (unsigned long)(ci->m_mCell(0,0)*10000.0+0.5), 32 );
			bs->WriteBits( (unsigned long)(ci->m_mCell(1,1)*10000.0+0.5), 32 );
			bs->WriteBits( (unsigned long)(ci->m_mCell(2,2)*10000.0+0.5), 32 );

		} else {

			bs->WriteBit(0);
			bs->WriteSignedBits( (long)(ci->m_mCell(0,0)*10000.0+0.5), 32 );
			bs->WriteSignedBits( (long)(ci->m_mCell(0,1)*10000.0+0.5), 32 );
			bs->WriteSignedBits( (long)(ci->m_mCell(0,2)*10000.0+0.5), 32 );
			bs->WriteSignedBits( (long)(ci->m_mCell(1,0)*10000.0+0.5), 32 );
			bs->WriteSignedBits( (long)(ci->m_mCell(1,1)*10000.0+0.5), 32 );
			bs->WriteSignedBits( (long)(ci->m_mCell(1,2)*10000.0+0.5), 32 );
			bs->WriteSignedBits( (long)(ci->m_mCell(2,0)*10000.0+0.5), 32 );
			bs->WriteSignedBits( (long)(ci->m_mCell(2,1)*10000.0+0.5), 32 );
			bs->WriteSignedBits( (long)(ci->m_mCell(2,2)*10000.0+0.5), 32 );
		}
	}
}


bool CBQBEngine::ImportCellInfo(CBQBOldCellInfo *ci, CBQBBitSet *bs) {

	double tfxx, tfxy, tfxz, tfyx, tfyy, tfyz, tfzx, tfzy, tfzz;

	if (bs->ReadBit()) { // Cubic

		tfxx = ((double)bs->ReadBitsInteger(32))/10000.0;
		ci->SetCellInfo( tfxx, tfxx, tfxx );

	} else { // Not Cubic

		if (bs->ReadBit()) { // Orthorhombic

			tfxx = ((double)bs->ReadBitsInteger(32))/10000.0;
			tfyy = ((double)bs->ReadBitsInteger(32))/10000.0;
			tfzz = ((double)bs->ReadBitsInteger(32))/10000.0;
			ci->SetCellInfo( tfxx, tfyy, tfzz );

		} else { // Non-orthorhombic

			tfxx = ((double)bs->ReadBitsSignedInteger(32))/10000.0;
			tfxy = ((double)bs->ReadBitsSignedInteger(32))/10000.0;
			tfxz = ((double)bs->ReadBitsSignedInteger(32))/10000.0;
			tfyx = ((double)bs->ReadBitsSignedInteger(32))/10000.0;
			tfyy = ((double)bs->ReadBitsSignedInteger(32))/10000.0;
			tfyz = ((double)bs->ReadBitsSignedInteger(32))/10000.0;
			tfzx = ((double)bs->ReadBitsSignedInteger(32))/10000.0;
			tfzy = ((double)bs->ReadBitsSignedInteger(32))/10000.0;
			tfzz = ((double)bs->ReadBitsSignedInteger(32))/10000.0;
			ci->SetCellInfo( tfxx, tfxy, tfxz, tfyx, tfyy, tfyz, tfzx, tfzy, tfzz );
		}
	}

	return true;
}


//#define DEBUG_VOL_ELEMENT_FROM 0
//#define DEBUG_VOL_ELEMENT_TO 2


bool CBQBEngine::CubeToIntegerArray(
		std::vector<int> &outp,
		int order,
		CBQBParameterSet_Volumetric *parm,
		int &histused,
		bool skipcomp,
		double *resi
	) {

	int z, z2, ti, tiov, tiuf, tiuf2, zx, zy, zz, zi, ii, ex, res, mantis, split;
	const CBQBCubeFrame *cfr;
	double tf=0, tf2, tf3, tfa2;
	int maxsymb;
	bool uf, of;


	if (parm == NULL) {
		m_IF.eprintf("CBQBEngine::CubeToIntegerArray(): Error: parm == NULL.\n");
		return false;
	}

/*	CDF df;
	if (resi != NULL) {
		df.m_iResolution = 250;
		df.m_fMinVal = 0.0;
		df.m_fMaxVal = 8.0;
		df.Create();
	}*/

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("** CubeToIntegerArray order=%d ***\n",order);

	split = ((unsigned long)1) << parm->GetSplit();

	cfr = m_pReadCacheCube->GetFrameHistory(0);

	if (!skipcomp) // Reserve some space (lower boundary estimate)
		outp.reserve(outp.size()+cfr->m_iResXYZ);
	
	maxsymb = pow10i(parm->GetSigni())-1;

	if (m_pExtrapolator != NULL) {
		if (order+1 < m_pExtrapolator->m_iCubeFrameVisibleCount)
			m_pExtrapolator->m_iCubeFrameVisibleCount = order+1;
		histused = MIN( m_pExtrapolator->m_iTRange, m_pExtrapolator->m_iCubeFrameVisibleCount ) - 1;
	} else
		histused = order;

	if (m_pExtrapolator != NULL) {

		if (m_IF.IsPL(BQB_PL_VERBOSE)) {
			m_IF.printf("Extrapolator running...\n");
			m_IF.printf("    [");
		}

		if (m_pExtrapolatorCorr != NULL) {

			m_faTempPred.resize(cfr->m_iResXYZ);
			m_faTempPred2.resize(cfr->m_iResXYZ);

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				tf = cfr->m_iRes[0] / 60.0;

			zi = 0;
			for (zx=0;zx<cfr->m_iRes[0];zx++) {

				if (m_IF.IsPL(BQB_PL_VERBOSE)) {
					if (fmod(zx,tf) < 1.0) {
						m_IF.printf("#");
						m_IF.FlushLog();
					}
				}

				for (zy=0;zy<cfr->m_iRes[1];zy++) {

					if ((zx >= m_pExtrapolatorCorr->m_iSOffset[0]) && 
						(zx <= m_pExtrapolatorCorr->m_iTempVal[0]) &&
						(zy >= m_pExtrapolatorCorr->m_iSOffset[1]) && 
						(zy <= m_pExtrapolatorCorr->m_iTempVal[1])) {

						for (zz=0;zz<m_pExtrapolatorCorr->m_iSOffset[2];zz++) {
							tf2 = m_pExtrapolator->ExtrapolatePred(zi);
							//m_faTempPred[zi] = ((double)cfr->m_iaMantis[zi]*pow10(cfr->m_iaExpo[zi])) - tf2;
							m_faTempPred[zi] = cfr->m_faBin[zi] - tf2;
							m_faTempPred2[zi] = tf2 + m_pExtrapolatorCorr->ExtrapolateCorr(m_faTempPred,zx,zy,zz,zi);
							zi++;
						}

						for (;zz<=m_pExtrapolatorCorr->m_iTempVal[2];zz++) {
							tf2 = m_pExtrapolator->ExtrapolatePred(zi);
							//m_faTempPred[zi] = ((double)cfr->m_iaMantis[zi]*pow10(cfr->m_iaExpo[zi])) - tf2;
							m_faTempPred[zi] = cfr->m_faBin[zi] - tf2;
							m_faTempPred2[zi] = tf2 + m_pExtrapolatorCorr->ExtrapolateKnownCorr(m_faTempPred,0,zi);
							zi++;
						}

						for (;zz<cfr->m_iRes[2];zz++) {
							tf2 = m_pExtrapolator->ExtrapolatePred(zi);
							//m_faTempPred[zi] = ((double)cfr->m_iaMantis[zi]*pow10(cfr->m_iaExpo[zi])) - tf2;
							m_faTempPred[zi] = cfr->m_faBin[zi] - tf2;
							m_faTempPred2[zi] = tf2 + m_pExtrapolatorCorr->ExtrapolateCorr(m_faTempPred,zx,zy,zz,zi);
							zi++;
						}

					} else {

						for (zz=0;zz<cfr->m_iRes[2];zz++) {
							tf2 = m_pExtrapolator->ExtrapolatePred(zi);
							//m_faTempPred[zi] = ((double)cfr->m_iaMantis[zi]*pow10(cfr->m_iaExpo[zi])) - tf2;
							m_faTempPred[zi] = cfr->m_faBin[zi] - tf2;
							m_faTempPred2[zi] = tf2 + m_pExtrapolatorCorr->ExtrapolateCorr(m_faTempPred,zx,zy,zz,zi);
							zi++;
						}
					}
				}
			}

		} else {

			m_faTempPred2.resize(cfr->m_iResXYZ);

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				tf = cfr->m_iRes[0] / 60.0;

			zi = 0;
			for (zx=0;zx<cfr->m_iRes[0];zx++) {

				if (m_IF.IsPL(BQB_PL_VERBOSE)) {
					if (fmod(zx,tf) < 1.0) {
						m_IF.printf("#");
						m_IF.FlushLog();
					}
				}

				for (zy=0;zy<cfr->m_iRes[1];zy++) {
					for (zz=0;zz<cfr->m_iRes[2];zz++) {
						m_faTempPred2[zi] = m_pExtrapolator->Extrapolate(zx,zy,zz,zi);
						zi++;
					}
				}
			}
		}

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("] Done.\n");

	} else { // No Extrapolator

		m_faTempPred.resize(cfr->m_iResXYZ);
		m_faTempPred2.resize(cfr->m_iResXYZ);
		for (z=0;z<cfr->m_iResXYZ;z++)
			m_faTempPred[z] = 0;

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Polynomial extrapolation...\n");

		for (z=1;z<order+1;z++) {
			cfr = m_pReadCacheCube->GetFrameHistory(z);
			tf = 1.0;
			for (z2=1;z2<=order;z2++) { 
				if (z2 == z)
					continue;
				tf *= z2 / ((double)z2 - z);
			}
			for (z2=0;z2<cfr->m_iResXYZ;z2++)
				m_faTempPred[z2] += tf * cfr->m_iaMantis[z2] * pow10(cfr->m_iaExpo[z2]);
		}

		cfr = m_pReadCacheCube->GetFrameHistory(0);

		if (parm->GetNbhFac() != 0) {

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Neighborhood extrapolation...\n");

			zi = 0;
			for (zx=0;zx<cfr->m_iRes[0];zx++) {
				for (zy=0;zy<cfr->m_iRes[1];zy++) {
					for (zz=0;zz<cfr->m_iRes[2];zz++) {

						tf2 = 0;
						tf3 = 0;

						if (zz > 0) {
							tf2++;
							tf3 += ((double)cfr->m_iaMantis[zi-1]*pow10(cfr->m_iaExpo[zi-1])) - m_faTempPred[zi-1];
						}

						if (zy > 0) {
							tf2++;
							tf3 += ((double)cfr->m_iaMantis[zi-cfr->m_iRes[2]]*pow10(cfr->m_iaExpo[zi-cfr->m_iRes[2]])) - m_faTempPred[zi-cfr->m_iRes[2]];
						}

						if (zx > 0) {
							tf2++;
							tf3 += ((double)cfr->m_iaMantis[zi-cfr->m_iResYZ]*pow10(cfr->m_iaExpo[zi-cfr->m_iResYZ])) - m_faTempPred[zi-cfr->m_iResYZ];
						}

						if (tf2 != 0)
							tf3 /= tf2 * parm->GetNbhFac();

						tf3 += m_faTempPred[zi];

						m_faTempPred2[zi] = tf3;

						zi++;
					}
				}
			}
		} else
			m_faTempPred2.assign(m_faTempPred.begin(),m_faTempPred.end());

	} // End No Extrapolator

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("Discretizing and serializing symbols...\n");

	tiov = 0;
	tiuf = 0;
	tiuf2 = 0;
	tfa2 = 0;
	for (zi=0;zi<cfr->m_iResXYZ;zi++) {

		if (parm->GetHilbert())
			ii = m_iaHilbertIdx[zi];
		else
			ii = zi;

		if (m_faTempPred2[ii] == 0)
			ex = -parm->GetEps();
		else
			ex = (int)floor(log10(fabs(m_faTempPred2[ii])))-parm->GetSigni()+1;

		if (ex < -parm->GetEps())
			ex = -parm->GetEps();

		uf = false;
		of = false;

_again:
		mantis = (int)floor(cfr->m_faBin[ii]*pow10(-ex)+0.5);
		// Bug (01.07.2019): If ex=-12 (as it happens by default, because GetEps() == 12),
		//                   result is larger than maximum int, so mantis is any (negative) garbage...

		if (floor(m_faTempPred2[ii]*pow10(-ex)+0.5+CR_EPS) != floor(m_faTempPred2[ii]*pow10(-ex)+0.5-CR_EPS)) {
			ti = (int)floor(m_faTempPred2[ii]*pow10(-ex));
	//		m_IF.printf("Sanitized rounding of %f.\n",m_faTempPred2[ii]*pow10(-ex));
		} else
			ti = (int)floor(m_faTempPred2[ii]*pow10(-ex)+0.5);

		res = mantis - ti;

//		if (mantis != 0)
//			m_IF.printf("@@ mantis=%d, ti=%d, res=%d, fabs(res)=%E, sqrt(abs(res))=%E\n",mantis,ti,res,fabs((double)res),sqrt(fabs((double)res)));

		if (resi != NULL) {
//			if (m_faTempPred2[ii] != 0) // Fix MB 01.07.2019
			// Fix MB 02.12.2020
				*resi += sqrt(fabs((double)res));

		//	if (sqrt(bqbabs(res)) != sqrt(bqbabs(res))) {
		//	if (!_finite(*resi)) {
		//		printf("@XX mantis=%d, ti=%d, res=%d, abs(res)=%d, sqrt(abs(res))=%f\n",mantis,ti,res,bqbabs(res),sqrt(bqbabs(res)));
		//		printf("ex=%d, m_faBin=%f, temppred=%f\n",ex,cfr->m_faBin[ii],m_faTempPred2[ii]);
		//	}
			//*resi += mypow(fabs(res),m_fEntropyExpoCube);
			//if (res != 0)
			//	df.AddToBin(log10(fabs(res)));
		}

		if (skipcomp)
			continue;

		if (ExpMantisEqual(cfr->m_iaExpo[ii],cfr->m_iaMantis[ii],ex,mantis) && (abs(res) <= maxsymb)) {
			if (uf)
				outp.push_back(C_UNDERFLOW);
			if (of)
				outp.push_back(C_OVERFLOW);
			if (abs(res) >= split) {
				#ifdef DEBUG_VOL_ELEMENT_FROM
					if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
						m_IF.printf("@> %d: Equal with split (%d): %d = %d & %d\n",ii,split,res,res/split,abs(res)%split);
				#endif
				outp.push_back(C_SPLIT);
				outp.push_back(res/split);
				if (res < 0)
					outp.push_back(-(abs(res)%split));
				else
					outp.push_back(res%split);
			} else {
				outp.push_back(res);
				#ifdef DEBUG_VOL_ELEMENT_FROM
					if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
						m_IF.printf("@> %d: Equal without split: %d\n",ii,res);
				#endif
			}
		} else if ((ex == cfr->m_iaExpo[ii]+1) && (ex > -parm->GetEps())) {
			tiuf2++;
			ex--;
	//		tia.push_back(maxsymb+2);
			uf = true;
			#ifdef DEBUG_VOL_ELEMENT_FROM
				if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
					m_IF.printf("@> %d: Soft Underflow 1.\n",ii);
			#endif
			goto _again;
		} else if (ex == cfr->m_iaExpo[ii]-1) {
			tiuf2++;
			ex++;
			of = true;
			#ifdef DEBUG_VOL_ELEMENT_FROM
				if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
					m_IF.printf("@> %d: Soft Underflow 2.\n",ii);
			#endif
			goto _again;
		} else if (ex > cfr->m_iaExpo[ii]) {
			tiuf++;
			outp.push_back(C_FULLNUMBER);
			outp.push_back(cfr->m_iaExpo[ii]);
			outp.push_back(cfr->m_iaMantis[ii]/split);
			if (cfr->m_iaMantis[ii] < 0)
				outp.push_back(-(abs(cfr->m_iaMantis[ii])%split));
			else
				outp.push_back(cfr->m_iaMantis[ii]%split);
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Underflow %d\n",ex-cfr->m_iaExpo[ii]);
			#ifdef DEBUG_VOL_ELEMENT_FROM
				if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO) && (!verbose))
					m_IF.printf("@> %d: Underflow %d\n",ii,ex-cfr->m_iaExpo[ii]);
			#endif
		} else if (abs(res) > maxsymb) {
			tiov++;
			outp.push_back(C_FULLNUMBER);
			#ifdef DEBUG_VOL_ELEMENT_FROM
				if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
					m_IF.printf("@> %d: @0: %d\n",ii,outp[outp.size()-1]);
			#endif
			outp.push_back(cfr->m_iaExpo[ii]);
			outp.push_back(cfr->m_iaMantis[ii]/split);
			if (cfr->m_iaMantis[ii] < 0)
				outp.push_back(-(abs(cfr->m_iaMantis[ii])%split));
			else
				outp.push_back(cfr->m_iaMantis[ii]%split);
			#ifdef DEBUG_VOL_ELEMENT_FROM
				if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
					m_IF.printf("@> %d: Full Number: %d %d %d\n",
						ii,cfr->m_iaExpo[ii],cfr->m_iaMantis[ii]/split,abs(cfr->m_iaMantis[ii])%split);
			#endif
		} else {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("UE: idx=%d, val=%.10f, pred=%.10f, cexpo=%d, pexpo=%d, cmantis=%d, pmantis=%d, res=%d.\n",
					ii,cfr->m_faBin[ii],m_faTempPred2[ii],cfr->m_iaExpo[ii],ex,cfr->m_iaMantis[ii],mantis,res);
			outp.push_back(C_FULLNUMBER);
			outp.push_back(cfr->m_iaExpo[ii]);
			outp.push_back(cfr->m_iaMantis[ii]/split);
			if (cfr->m_iaMantis[ii] < 0)
				outp.push_back(-(abs(cfr->m_iaMantis[ii])%split));
			else
				outp.push_back(cfr->m_iaMantis[ii]%split);
		}

		#ifdef DEBUG_VOL_ELEMENT_FROM
			if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
				m_IF.printf("@> UE: idx=%6d, val=%.10f, pred=%.22f, cexpo=%d, pexpo=%d, cmantis=%d, pmantis=%d, res=%d.\n",
					ii,cfr->m_faBin[ii],m_faTempPred2[ii],cfr->m_iaExpo[ii],ex,cfr->m_iaMantis[ii],mantis,res);
		#endif

		tfa2 += res;
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		tfa2 /= cfr->m_iResXYZ;
		m_IF.printf("      tfa2 = %8.6f\n",tfa2);
		m_IF.printf("      Overflow: %d, Underflow: %d, Soft Underflow: %d\n",tiov,tiuf,tiuf2);
		m_IF.printf("      Produced %lu symbols.\n",(unsigned long)outp.size());
	}

/*	if (resi != 0) {
		df.NormBinIntegral(1000.0);
		CxString buf;
		buf.sprintf("histo%06d.csv",gc_iHistogramCounter++);
		df.Write("",(const char*)buf,"",false);
	}*/

	return true;
}


bool CBQBEngine::IntegerArrayToCube(
		std::vector<int> &inp,
		int order,
		CBQBParameterSet_Volumetric *parm,
		int ver,
		bool second
	) {

	int z, z2, zx, zy, zz, zi, ii, ex, split;
	CBQBCubeFrame *ofr, *cfr;
	double tf, tf2, tf3;
	std::vector<int> tresi, texpo, tmantis;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("** IntegerArrayToCube order=%d ***\n",order);
	
	if (ver >= 1)
		split = ((unsigned long)1) << parm->GetSplit();
	else
		split = parm->GetSplit(); // Backward compatibility with bug

	if (second)
		ofr = GetOutput2CubeFrame(0);
	else
		ofr = GetOutputCubeFrame(0);

	m_faTempPred.resize(ofr->m_iResXYZ);

	if (m_pExtrapolator == NULL) {

		for (z=0;z<ofr->m_iResXYZ;z++)
			m_faTempPred[z] = 0;

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Polynomial extrapolation...\n");

		for (z=1;z<order+1;z++) {
			if (second)
				cfr = GetOutput2CubeFrame(z);
			else
				cfr = GetOutputCubeFrame(z);
			if (cfr == NULL) {
				m_IF.eprintf("CBQBEngine::IntegerArrayToCube(): Error: Frame is of temporal order %d, step history is too short (%d).\n",
					order,z);
				abort();
			}
			tf = 1.0;
			for (z2=1;z2<=order;z2++) { 
				if (z2 == z)
					continue;
				tf *= z2 / ((double)z2 - z);
			}
			for (z2=0;z2<ofr->m_iResXYZ;z2++)
				m_faTempPred[z2] += tf * cfr->m_iaMantis[z2] * pow10(cfr->m_iaExpo[z2]);
		}

	} else {

		if (order+1 > m_pExtrapolator->m_iCubeFrameVisibleCount) {
			m_IF.eprintf("CBQBEngine::IntegerArrayToCube(): Error: Frame is of temporal order %d, step history is too short (%d).\n",
				order,m_pExtrapolator->m_iCubeFrameVisibleCount);
			abort();
		}

		m_pExtrapolator->m_iCubeFrameVisibleCount = order+1;

		if (m_pExtrapolatorCorr != NULL)
			for (z2=0;z2<ofr->m_iResXYZ;z2++)
				m_faTempPred[z2] = m_pExtrapolator->ExtrapolatePred(z2);
	}

	tresi.resize(ofr->m_iResXYZ);
	texpo.resize(ofr->m_iResXYZ);
	tmantis.resize(ofr->m_iResXYZ);

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("Decoding input...\n");

	zi = 0;
	z = 0;
	while (zi < (int)inp.size()) {
		if (parm->GetHilbert())
			ii = m_iaHilbertIdx[z];
		else
			ii = z;
		if (inp[zi] == C_FULLNUMBER) { // Full Number
			texpo[ii] = inp[++zi];
			tmantis[ii] = inp[++zi]*split;
			tmantis[ii] += inp[zi+1];
	//		if (inp[zi+1] < 0)
	//			tmantis[ii] = -tmantis[ii];
			#ifdef DEBUG_VOL_ELEMENT_FROM
				if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
					m_IF.printf("@< %d: Full Number: %d %d %d --> %d %d\n",
						ii,inp[zi-1],inp[zi],inp[zi+1],texpo[ii],tmantis[ii]);
			#endif
			zi++;
		} else if (inp[zi] == C_UNDERFLOW) { // Soft Underflow
			if (inp[++zi] == C_SPLIT) {
				#ifdef DEBUG_VOL_ELEMENT_FROM
					if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
						m_IF.printf("@< %d: Soft Underflow with split: %d & %d = %d\n",
							ii,inp[zi+1],inp[zi+2],inp[zi+1]*split+inp[zi+2]);
				#endif
				tresi[ii] = inp[++zi]*split;
				tresi[ii] += inp[zi+1];
	//			if (inp[zi+1] < 0)
	//				tresi[ii] = -tresi[ii];
				zi++;
				#ifdef DEBUG_VOL_ELEMENT_FROM
					if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
						m_IF.printf("@< %d:  --> %d\n",ii,tresi[ii]);
				#endif
			} else {
				#ifdef DEBUG_VOL_ELEMENT_FROM
					if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
						m_IF.printf("@< %d: Soft Underflow without split.\n",ii);
				#endif
				tresi[ii] = inp[zi];
			}
			texpo[ii] = -101;
		} else if (inp[zi] == C_OVERFLOW) { // Soft Overflow
			if (inp[++zi] == C_SPLIT) {
				#ifdef DEBUG_VOL_ELEMENT_FROM
					if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
						m_IF.printf("@< %d: Soft Underflow with split: %d & %d = %d\n",
							ii,inp[zi+1],inp[zi+2],inp[zi+1]*split+inp[zi+2]);
				#endif
				tresi[ii] = inp[++zi]*split;
				tresi[ii] += inp[zi+1];
	//			if (inp[zi+1] < 0)
	//				tresi[ii] = -tresi[ii];
				zi++;
				#ifdef DEBUG_VOL_ELEMENT_FROM
					if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
						m_IF.printf("@< %d:  --> %d\n",ii,tresi[ii]);
				#endif
			} else {
				#ifdef DEBUG_VOL_ELEMENT_FROM
					if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
						m_IF.printf("@< %d: Soft Underflow without split.\n",ii);
				#endif
				tresi[ii] = inp[zi];
			}
			texpo[ii] = -102;
		} else if (inp[zi] == C_SPLIT) { // Splitted Residue
			tresi[ii] = inp[++zi]*split;
			tresi[ii] += inp[zi+1];
	//		if (inp[zi+1] < 0)
	//			tresi[ii] = -tresi[ii];
			zi++;
			texpo[ii] = -100;
			#ifdef DEBUG_VOL_ELEMENT_FROM
				if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
					m_IF.printf("@< %d: Splitted residue.\n",ii);
			#endif
		} else { // Standard Residue
			tresi[ii] = inp[zi];
			texpo[ii] = -100;
			#ifdef DEBUG_VOL_ELEMENT_FROM
				if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
					m_IF.printf("@< %d: Standard residue.\n",ii);
			#endif
		}
		#ifdef DEBUG_VOL_ELEMENT_FROM
			if ((ii >= DEBUG_VOL_ELEMENT_FROM) && (ii <= DEBUG_VOL_ELEMENT_TO))
				m_IF.printf("@< %d: resi=%d, expo=%d, mantis=%d\n",ii,tresi[ii],texpo[ii],tmantis[ii]);
		#endif
		zi++;
		z++;
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("Read %d symbols and %d grid points.\n",zi,z);
		m_IF.printf("Neighborhood extrapolation...\n");
	}

	zi = 0;
	for (zx=0;zx<ofr->m_iRes[0];zx++) {
		for (zy=0;zy<ofr->m_iRes[1];zy++) {
			for (zz=0;zz<ofr->m_iRes[2];zz++) {

				if (texpo[zi] > -100) {
					ofr->m_iaMantis[zi] = tmantis[zi];
					ofr->m_iaExpo[zi] = (char)texpo[zi];
					ofr->m_faBin[zi] = ofr->m_iaMantis[zi] * pow10(ofr->m_iaExpo[zi]);
					if (m_pExtrapolatorCorr != NULL) {
						#ifdef DEBUG_VOL_ELEMENT_FROM
							if ((zi >= DEBUG_VOL_ELEMENT_FROM) && (zi <= DEBUG_VOL_ELEMENT_TO))
								m_IF.printf("@< %d: X-Updating TempPred %.10G with true %.10G: %.10G\n",
									zi,m_faTempPred[zi],ofr->m_faBin[zi],ofr->m_faBin[zi]-m_faTempPred[zi]);
						#endif
						m_faTempPred[zi] = ofr->m_faBin[zi] - m_faTempPred[zi];
					}
					zi++;
					continue;
				}

				if (m_pExtrapolator != NULL) {

					if (m_pExtrapolatorCorr != NULL) {
						tf3 = m_faTempPred[zi] + m_pExtrapolatorCorr->ExtrapolateCorr(m_faTempPred,zx,zy,zz,zi);
						#ifdef DEBUG_VOL_ELEMENT_FROM
							if ((zi >= DEBUG_VOL_ELEMENT_FROM) && (zi <= DEBUG_VOL_ELEMENT_TO))
								m_IF.printf("@< %d: Pred %.10G, Corr %.10G --> %.10G\n",
									zi,m_faTempPred[zi],tf3-m_faTempPred[zi],tf3);
						#endif
					} else
						tf3 = m_pExtrapolator->Extrapolate(zx,zy,zz,zi);

				} else { // Classical Neighborhood correction

					tf2 = 0;
					tf3 = 0;

					if (parm->GetNbhFac() != 0) {

						if (zz > 0) {
							tf2++;
							tf3 += ((double)ofr->m_iaMantis[zi-1]*pow10(ofr->m_iaExpo[zi-1])) - m_faTempPred[zi-1];
						}

						if (zy > 0) {
							tf2++;
							tf3 += ((double)ofr->m_iaMantis[zi-ofr->m_iRes[2]]*pow10(ofr->m_iaExpo[zi-ofr->m_iRes[2]])) - m_faTempPred[zi-ofr->m_iRes[2]];
						}

						if (zx > 0) {
							tf2++;
							tf3 += ((double)ofr->m_iaMantis[zi-ofr->m_iResYZ]*pow10(ofr->m_iaExpo[zi-ofr->m_iResYZ])) - m_faTempPred[zi-ofr->m_iResYZ];
						}

						if (tf2 != 0)
							tf3 /= tf2 * parm->GetNbhFac();
					}

					tf3 += m_faTempPred[zi];
				}


				if (tf3 == 0)
					ex = -parm->GetEps();
				else
					ex = (int)floor(log10(fabs(tf3)))-parm->GetSigni()+1;

				if (ex < -parm->GetEps())
					ex = -parm->GetEps();

				if (texpo[zi] == -101)
					ex--;
				else if (texpo[zi] == -102)
					ex++;

				if (floor(tf3*pow10(-ex)+0.5+CR_EPS) != floor(tf3*pow10(-ex)+0.5-CR_EPS)) {
					ofr->m_iaMantis[zi] = (int)floor(tf3*pow10(-ex));
	//				m_IF.printf("Sanitized rounding of %f.\n",tf3*pow10(-ex));
				} else
					ofr->m_iaMantis[zi] = (int)floor(tf3*pow10(-ex)+0.5);

				ofr->m_iaMantis[zi] += tresi[zi];
				ofr->m_iaExpo[zi] = (char)ex;
				//while (ofr->m_iaMantis[zi] >= pow10i(signi)) {
				while (abs(ofr->m_iaMantis[zi]) >= pow10i(parm->GetSigni())) {
					#ifdef DEBUG_VOL_ELEMENT_FROM
						if ((zi >= DEBUG_VOL_ELEMENT_FROM) && (zi <= DEBUG_VOL_ELEMENT_TO))
							m_IF.printf("@< %d: A %d %d %d\n",zi,ofr->m_iaMantis[zi],ofr->m_iaExpo[zi],parm->GetSigni());
					#endif
					ofr->m_iaMantis[zi] /= 10;
					ofr->m_iaExpo[zi]++;
				}
				//while ((ofr->m_iaMantis[zi] < pow10i(signi-1)) && (ofr->m_iaExpo[zi] > -eps)) {
				while ((abs(ofr->m_iaMantis[zi]) < pow10i(parm->GetSigni()-1)) && (ofr->m_iaExpo[zi] > -parm->GetEps())) {
					#ifdef DEBUG_VOL_ELEMENT_FROM
						if ((zi >= DEBUG_VOL_ELEMENT_FROM) && (zi <= DEBUG_VOL_ELEMENT_TO))
							m_IF.printf("@< %d: B %d %d %d\n",zi,ofr->m_iaMantis[zi],ofr->m_iaExpo[zi],parm->GetSigni());
					#endif
					ofr->m_iaMantis[zi] *= 10;
					ofr->m_iaExpo[zi]--;
				}
				ofr->m_faBin[zi] = ofr->m_iaMantis[zi] * pow10(ofr->m_iaExpo[zi]);

				if (m_pExtrapolatorCorr != NULL) {
					#ifdef DEBUG_VOL_ELEMENT_FROM
						if ((zi >= DEBUG_VOL_ELEMENT_FROM) && (zi <= DEBUG_VOL_ELEMENT_TO))
							m_IF.printf("@< %d: Updating TempPred %.10G with true %.10G: %.10G\n",
								zi,m_faTempPred[zi],ofr->m_faBin[zi],ofr->m_faBin[zi]-m_faTempPred[zi]);
					#endif
					m_faTempPred[zi] = ofr->m_faBin[zi] - m_faTempPred[zi];
				}

				#ifdef DEBUG_VOL_ELEMENT_FROM
					if ((zi >= DEBUG_VOL_ELEMENT_FROM) && (zi <= DEBUG_VOL_ELEMENT_TO))
						m_IF.printf("@< %d: Pred=%.22f, expo=%d, predmantis=%d, resi=%d, mantis=%d, result=%.10f\n",
							zi,tf3,ex,(int)floor(tf3*pow10(-ex)+0.5),tresi[zi],ofr->m_iaMantis[zi],ofr->m_faBin[zi]);
				#endif

				zi++;
			}
		}
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("All done.\n");

	return true;
}


void CBQBEngine::ExportCubeHeader(CBQBBitSet *bs, int order) {

	const CBQBCubeFrame *cfr;
	int i;

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("Writing cube header...\n");

	if (order != 0) {
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Cube header already written before.\n");
		bs->WriteBit(0);
		return;
	} else
		bs->WriteBit(1);

	i = bs->GetLength();

	cfr = m_pReadCacheCube->GetFrameHistory(0);

	bs->WriteBits(cfr->m_iRes[0],10);
	bs->WriteBits(cfr->m_iRes[1],10);
	bs->WriteBits(cfr->m_iRes[2],10);

	bs->WriteSignedBits(cfr->m_iCenter[0],32);
	bs->WriteSignedBits(cfr->m_iCenter[1],32);
	bs->WriteSignedBits(cfr->m_iCenter[2],32);

	bs->WriteSignedBits(cfr->m_iStrideA[0],32);
	bs->WriteSignedBits(cfr->m_iStrideA[1],32);
	bs->WriteSignedBits(cfr->m_iStrideA[2],32);

	bs->WriteSignedBits(cfr->m_iStrideB[0],32);
	bs->WriteSignedBits(cfr->m_iStrideB[1],32);
	bs->WriteSignedBits(cfr->m_iStrideB[2],32);

	bs->WriteSignedBits(cfr->m_iStrideC[0],32);
	bs->WriteSignedBits(cfr->m_iStrideC[1],32);
	bs->WriteSignedBits(cfr->m_iStrideC[2],32);

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("%d bytes written.\n",(bs->GetLength()-i)/8);
}


void CBQBEngine::ImportCubeHeader(CBQBBitSet *bs, bool second) {

	CBQBCubeFrame *ofr, *ofr2;
	int z, i;

	if (second)
		ofr = GetOutput2CubeFrame(0);
	else
		ofr = GetOutputCubeFrame(0);

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("Reading cube header...\n");

	i = bs->GetReadPos();

	if (!bs->ReadBit()) { // Header already found in previous frame
		if (second)
			ofr2 = GetOutput2CubeFrame(1);
		else
			ofr2 = GetOutputCubeFrame(1);
		if (ofr2 == NULL) {
			m_IF.printf("CEngine::ImportCubeHeader(): Error: First frame does not contain header.\n");
			return;
		}
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Cube header already read before.\n");
		for (z=0;z<3;z++) {
			ofr->m_iRes[z] = ofr2->m_iRes[z];
			ofr->m_iCenter[z] = ofr2->m_iCenter[z];
			ofr->m_iStrideA[z] = ofr2->m_iStrideA[z];
			ofr->m_iStrideB[z] = ofr2->m_iStrideB[z];
			ofr->m_iStrideC[z] = ofr2->m_iStrideC[z];
			ofr->m_fCenter[z] = ofr2->m_fCenter[z];
			ofr->m_fStrideA[z] = ofr2->m_fStrideA[z];
			ofr->m_fStrideB[z] = ofr2->m_fStrideB[z];
			ofr->m_fStrideC[z] = ofr2->m_fStrideC[z];
		}
	} else { // Header here

		for (z=0;z<3;z++)
			ofr->m_iRes[z] = bs->ReadBitsInteger(10);

		for (z=0;z<3;z++) {
			ofr->m_iCenter[z] = bs->ReadBitsSignedInteger(32);
			ofr->m_fCenter[z] = FixedToFloat(ofr->m_iCenter[z],6);
		}

		for (z=0;z<3;z++) {
			ofr->m_iStrideA[z] = bs->ReadBitsSignedInteger(32);
			ofr->m_fStrideA[z] = FixedToFloat(ofr->m_iStrideA[z],6);
		}

		for (z=0;z<3;z++) {
			ofr->m_iStrideB[z] = bs->ReadBitsSignedInteger(32);
			ofr->m_fStrideB[z] = FixedToFloat(ofr->m_iStrideB[z],6);
		}

		for (z=0;z<3;z++) {
			ofr->m_iStrideC[z] = bs->ReadBitsSignedInteger(32);
			ofr->m_fStrideC[z] = FixedToFloat(ofr->m_iStrideC[z],6);
		}
	}

	ofr->m_fMinVal[0] = ofr->m_fCenter[0];
	ofr->m_fMaxVal[0] = ofr->m_fCenter[0] + ofr->m_fStrideA[0] * ofr->m_iRes[0];
	ofr->m_fMinVal[1] = ofr->m_fCenter[1];
	ofr->m_fMaxVal[1] = ofr->m_fCenter[1] + ofr->m_fStrideB[1] * ofr->m_iRes[1];
	ofr->m_fMinVal[2] = ofr->m_fCenter[2];
	ofr->m_fMaxVal[2] = ofr->m_fCenter[2] + ofr->m_fStrideC[2] * ofr->m_iRes[2];

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("Resolution: %d x %d x %d\n",ofr->m_iRes[0],ofr->m_iRes[1],ofr->m_iRes[2]);
		m_IF.printf("Center: %f | %f | %f\n",ofr->m_fCenter[0],ofr->m_fCenter[1],ofr->m_fCenter[2]);
		m_IF.printf("Stride A: %f | %f | %f\n",ofr->m_fStrideA[0],ofr->m_fStrideA[1],ofr->m_fStrideA[2]);
		m_IF.printf("Stride B: %f | %f | %f\n",ofr->m_fStrideB[0],ofr->m_fStrideB[1],ofr->m_fStrideB[2]);
		m_IF.printf("Stride C: %f | %f | %f\n",ofr->m_fStrideC[0],ofr->m_fStrideC[1],ofr->m_fStrideC[2]);
		m_IF.printf("Range: X %f - %f, Y %f - %f, Z %f - %f\n",
			ofr->m_fMinVal[0],ofr->m_fMaxVal[0],ofr->m_fMinVal[1],ofr->m_fMaxVal[1],ofr->m_fMinVal[2],ofr->m_fMaxVal[2]);
		m_IF.printf("%d bytes read.\n",(bs->GetReadPos()-i)/8);
	}

	ofr->m_iResXY = ofr->m_iRes[0] * ofr->m_iRes[1];
	ofr->m_iResYZ = ofr->m_iRes[1] * ofr->m_iRes[2];
	ofr->m_iResXYZ = ofr->m_iRes[0] * ofr->m_iRes[1] * ofr->m_iRes[2];

	ofr->m_faBin.resize(ofr->m_iResXYZ);
	ofr->m_iaExpo.resize(ofr->m_iResXYZ);
	ofr->m_iaMantis.resize(ofr->m_iResXYZ);
}


//#define DEBUG_POS_ELEMENT 0


void CBQBEngine::AtomsToIntegerArray(
		std::vector<int> &outp,
		int order,
		CBQBParameterSet_Position *parm,
		int &histused,
		bool skipcomp,
		double *resi
	) {

	int z, z2, i, zc, zi;
	long mantis, ti, res;
	int split;
	double tf;
	const CBQBAtomSet *cfr;
	std::vector<double> tpred;
	double compred;


	if (parm == NULL) {
		m_IF.eprintf("CBQBEngine::AtomsToIntegerArray(): Error: parm == NULL.\n");
		abort();
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("\n");
		m_IF.printf("*** AtomsToIntegerArray order=%d, signi=%d, esplit=%d, split=%lu ***\n",
			order,parm->GetPrecision(),parm->GetSplit(),((unsigned long)1)<<parm->GetSplit());
	}

	i = (int)outp.size();

	if (m_pReadCacheCube != NULL)
		cfr = m_pReadCacheCube->GetFrameHistory(0)->m_pAtoms;
	else
		cfr = m_pReadCacheXYZ->GetFrameHistory(0);

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("Processing %lu atoms...\n",(unsigned long)cfr->m_oaAtoms.size());

	split = ((unsigned long)1)<<parm->GetSplit();

	tpred.resize(cfr->m_oaAtoms.size());

	if (m_pExtrapolatorXYZ != NULL) {
		if (order+1 < m_pExtrapolatorXYZ->m_iCubeFrameVisibleCount)
			m_pExtrapolatorXYZ->m_iCubeFrameVisibleCount = order+1;
		histused = MIN( m_pExtrapolatorXYZ->m_iTRange, m_pExtrapolatorXYZ->m_iCubeFrameVisibleCount ) - 1;
	} else
		histused = order;

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">> histused %d\n",histused);

	for (zc=0;zc<3;zc++) {

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("*** Now coordinate %d/3...\n",zc+1);

		if (m_pExtrapolatorXYZ != NULL) {

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Extrapolator...\n");

			for (z=0;z<(int)cfr->m_oaAtoms.size();z++)
				tpred[z] = m_pExtrapolatorXYZ->ExtrapolateXYZ(m_iaAtomSort[z],zc);

			compred = m_pExtrapolatorXYZ->ExtrapolateXYZ((int)cfr->m_oaAtoms.size(),zc);

		} else {
	
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Polynomial extrapolation...\n");

			for (z=0;z<(int)cfr->m_oaAtoms.size();z++)
				tpred[z] = 0;
			compred = 0;

			for (z=1;z<order+1;z++) {

				if (m_pReadCacheCube != NULL)
					cfr = m_pReadCacheCube->GetFrameHistory(z)->m_pAtoms;
				else
					cfr = m_pReadCacheXYZ->GetFrameHistory(z);

				tf = 1.0;
				for (z2=1;z2<=order;z2++) { 
					if (z2 == z)
						continue;
					tf *= z2 / ((double)z2 - z);
				}
				for (z2=0;z2<(int)cfr->m_oaAtoms.size();z2++)
					tpred[z2] += tf * cfr->m_oaAtoms[m_iaAtomSort[z2]]->m_fRelCoord[zc];

				compred += tf * cfr->m_faCOM[zc];
			}
		}

		if (m_pReadCacheCube != NULL)
			cfr = m_pReadCacheCube->GetFrameHistory(0)->m_pAtoms;
		else
			cfr = m_pReadCacheXYZ->GetFrameHistory(0);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Discretizing and serializing symbols...\n");

		// COM component
		mantis = cfr->m_iaCOM[zc];
		if (floor(compred*pow10(parm->GetPrecision())+0.5+AR_EPS) != floor(compred*pow10(parm->GetPrecision())+0.5-AR_EPS)) {
			ti = (int)floor(compred*pow10(parm->GetPrecision()));
//			m_IF.printf("1A) Sanitized rounding of %f.\n",compred*pow10(asigni));
		} else
			ti = (int)floor(compred*pow10(parm->GetPrecision())+0.5);

		res = mantis - ti;

		#ifdef DEBUG_POS_ELEMENT
			m_IF.printf("@> CoM zc=%d, true %.6f, pred %.6f, mantis %ld, ti %ld, res %ld.\n",
				zc,cfr->m_faCOM[zc],compred,mantis,ti,res);
		#endif

		if (resi != NULL)
			*resi += bqbpow2((double)res);

		if (!skipcomp) {

			if (bqbabs(res) >= split) {
				outp.push_back(C_SPLIT);
				outp.push_back(res/split);
				if (res < 0)
					outp.push_back(-(bqbabs(res)%split));
				else
					outp.push_back(res%split);
			} else {
				outp.push_back(res);
			}
		}

		for (zi=0;zi<(int)cfr->m_oaAtoms.size();zi++) {

			mantis = cfr->m_oaAtoms[m_iaAtomSort[zi]]->m_iRelCoord[zc];

			if (floor(tpred[zi]*pow10(parm->GetPrecision())+0.5+AR_EPS) != floor(tpred[zi]*pow10(parm->GetPrecision())+0.5-AR_EPS)) {
				ti = (int)floor(tpred[zi]*pow10(parm->GetPrecision()));
//				m_IF.printf("1B) %3d Sanitized rounding of %f.\n",zi,tpred[zi]*pow10(asigni));
			} else
				ti = (int)floor(tpred[zi]*pow10(parm->GetPrecision())+0.5);

			res = mantis - ti;

			#ifdef DEBUG_POS_ELEMENT
				if (m_iaAtomSort[zi] == DEBUG_POS_ELEMENT)
					m_IF.printf("@> A%d zc=%d, true %.6f, pred %.6f, mantis %ld, ti %ld, res %ld.\n",
						m_iaAtomSort[zi],zc,cfr->m_oaAtoms[m_iaAtomSort[zi]]->m_fRelCoord[zc],tpred[zi],mantis,ti,res);
			#endif

			if (resi != NULL)
				*resi += bqbpow2((double)res);

			if (skipcomp)
				continue;

			if (bqbabs(res) >= split) {
				outp.push_back(C_SPLIT);
				outp.push_back(res/split);
				if (res < 0) {
					outp.push_back(-(bqbabs(res)%split));
					#ifdef DEBUG_POS_ELEMENT
						if (m_iaAtomSort[zi] == DEBUG_POS_ELEMENT)
							m_IF.printf("@> A%d zc=%d, Negative Split(%d): %ld | %ld = %ld.\n",
								m_iaAtomSort[zi],zc,split,res/split,-(abs(res)%split),res);
					#endif
				} else {
					outp.push_back(res%split);
					#ifdef DEBUG_POS_ELEMENT
						if (m_iaAtomSort[zi] == DEBUG_POS_ELEMENT)
							m_IF.printf("@> A%d zc=%d, Positive Split(%d): %ld | %ld = %ld.\n",
								m_iaAtomSort[zi],zc,split,res/split,res%split,res);
					#endif
				}

			} else {
				outp.push_back(res);
				#ifdef DEBUG_POS_ELEMENT
					if (m_iaAtomSort[zi] == DEBUG_POS_ELEMENT)
						m_IF.printf("@> A%d zc=%d, No Split: %ld.\n",
							m_iaAtomSort[zi],zc,res);
				#endif
			}
		}
	}

	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		if ((int)outp.size()-i > 20) {
			m_IF.printf("First 20 symbols (from %lu):\n",(unsigned long)outp.size());
			for (z=0;z<20;z++)
				m_IF.printf("%2d: %12d\n",z,outp[z+i]);
			m_IF.printf("\n");
		}
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("*** AtomsToIntegerArray finished, %ld symbols written ***\n",((long)outp.size())-i);
		m_IF.printf("\n");
	}
}


void CBQBEngine::IntegerArrayToAtoms(
		const std::vector<int> &inp,
		int order,
		CBQBParameterSet_Position *parm,
		int ver,
		bool second
	) {

	int z, z2, i=0, zc, zi, ii;
	std::vector<double> tpred;
	double compred;
	long ti;
	int split;
	double tf;
	CBQBAtomSet *ofr, *ofr0;


	ii = 0;

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("\n");
		m_IF.printf("*** IntegerArrayToAtoms ***\n");
	}

	i = ii;

	if (ver >= 1)
		split = ((unsigned long)1)<<parm->GetSplit();
	else
		split = parm->GetSplit(); // Backward compatibility with bug

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("Order=%d, asigni=%d, esplit=%d, split=%d\n",
			order,parm->GetPrecision(),parm->GetSplit(),split);

	if (second)
		ofr0 = GetOutput2AtomFrame(0);
	else
		ofr0 = GetOutputAtomFrame(0);

	tpred.resize(ofr0->m_oaAtoms.size());

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("First 20 symbols:\n");
		for (z=0;z<20;z++)
			m_IF.printf("%2d: %12d\n",z,inp[z]);
		m_IF.printf("\n");
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		if (m_pExtrapolatorXYZ != NULL)
			m_IF.printf("<< histused %d, visco %d, trange %d\n",
				order,m_pExtrapolatorXYZ->m_iCubeFrameVisibleCount,parm->GetExtraTRange());

	for (zc=0;zc<3;zc++) {

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("*** Now coordinate %d/3...\n",zc+1);
	
		if (m_pExtrapolatorXYZ != NULL) {

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Extrapolator...\n");

			if (order+1 > m_pExtrapolatorXYZ->m_iCubeFrameVisibleCount) {
				m_IF.eprintf("CBQBEngine::IntegerArrayToAtoms(): Error: Frame is of temporal order %d, step history is too short (%d).\n",
					order,m_pExtrapolatorXYZ->m_iCubeFrameVisibleCount);
				abort();
			}

			m_pExtrapolatorXYZ->m_iCubeFrameVisibleCount = order+1;

			for (z=0;z<(int)ofr0->m_oaAtoms.size();z++)
				tpred[z] = m_pExtrapolatorXYZ->ExtrapolateXYZ(m_iaAtomSort[z],zc);

			compred = m_pExtrapolatorXYZ->ExtrapolateXYZ((int)ofr0->m_oaAtoms.size(),zc);

		} else {

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Polynomial extrapolation...\n");

			for (z=0;z<(int)ofr0->m_oaAtoms.size();z++)
				tpred[z] = 0;
			compred = 0;

			for (z=1;z<order+1;z++) {
				if (second)
					ofr = GetOutput2AtomFrame(z);
				else
					ofr = GetOutputAtomFrame(z);
				if (ofr == NULL) {
					m_IF.eprintf("CBQBEngine::IntegerArrayToAtoms(): Error: Frame is of order %d, step history is too short (%d).\n",
						order,z);
					abort();
				}
				tf = 1.0;
				for (z2=1;z2<=order;z2++) { 
					if (z2 == z)
						continue;
					tf *= z2 / ((double)z2 - z);
				}
				for (z2=0;z2<(int)ofr0->m_oaAtoms.size();z2++)
					tpred[z2] += tf * ofr->m_oaAtoms[m_iaAtomSort[z2]]->m_fRelCoord[zc];
				compred += tf * ofr->m_faCOM[zc];
			}
		}

	//	if (second)
	//		ofr = GetOutput2AtomFrame(0);
	//	else
	//		ofr = GetOutputAtomFrame(0);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Decoding symbols...\n");

		// COM component
		if (floor(compred*pow10(parm->GetPrecision())+0.5+AR_EPS) != floor(compred*pow10(parm->GetPrecision())+0.5-AR_EPS)) {
			ti = (int)floor(compred*pow10(parm->GetPrecision()));
//			m_IF.printf("2A) Sanitized rounding of %f.\n",compred*pow10(asigni));
		} else
			ti = (int)floor(compred*pow10(parm->GetPrecision())+0.5);
		if (inp[ii] == C_SPLIT) {
			ii++;
			ofr0->m_iaCOM[zc] = ti + inp[ii++] * split;
			ofr0->m_iaCOM[zc] += inp[ii++];
		} else
			ofr0->m_iaCOM[zc] = ti + inp[ii++];
		ofr0->m_faCOM[zc] = FixedToFloat(ofr0->m_iaCOM[zc],parm->GetPrecision());

		#ifdef DEBUG_POS_ELEMENT
			m_IF.printf("@< CoM zc=%d, pred %.6f, mantis %d, ti %ld, true %.6f.\n",
				zc,compred,ofr0->m_iaCOM[zc],ti,ofr0->m_faCOM[zc]);
		#endif

		for (zi=0;zi<(int)ofr0->m_oaAtoms.size();zi++) {

			if (floor(tpred[zi]*pow10(parm->GetPrecision())+0.5+AR_EPS) != floor(tpred[zi]*pow10(parm->GetPrecision())+0.5-AR_EPS)) {
				ti = (int)floor(tpred[zi]*pow10(parm->GetPrecision()));
//				m_IF.printf("2B) %3d Sanitized rounding of %f.\n",zi,tpred[zi]*pow10(asigni));
			} else
				ti = (int)floor(tpred[zi]*pow10(parm->GetPrecision())+0.5);

			if (inp[ii] == C_SPLIT) {
				ii++;
				ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_iRelCoord[zc] = ti + inp[ii++] * split;
				ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_iRelCoord[zc] += inp[ii++];
				#ifdef DEBUG_POS_ELEMENT
					if (m_iaAtomSort[zi] == DEBUG_POS_ELEMENT)
						m_IF.printf("@< A%d ii=%d Split(%d): %ld + %d*%d + %d = %ld.\n",
							m_iaAtomSort[zi],ii-3,split,ti,inp[ii-2],split,inp[ii-1],
							ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_iRelCoord[zc]);
				#endif
			} else {
				ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_iRelCoord[zc] = ti + inp[ii++];
				#ifdef DEBUG_POS_ELEMENT
					if (m_iaAtomSort[zi] == DEBUG_POS_ELEMENT)
						m_IF.printf("@< A%d ii=%d No Split: %ld + %d = %ld.\n",
							m_iaAtomSort[zi],ii-3,ti,inp[ii-1],ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_iRelCoord[zc]);
				#endif
			}

			ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_iCoord[zc] = ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_iRelCoord[zc] + ofr0->m_iaCOM[zc];

			ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_fRelCoord[zc] = FixedToFloat(ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_iRelCoord[zc],parm->GetPrecision());
			ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_fCoord[zc] = FixedToFloat(ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_iCoord[zc],parm->GetPrecision());

			#ifdef DEBUG_POS_ELEMENT
				if (m_iaAtomSort[zi] == DEBUG_POS_ELEMENT)
					m_IF.printf("@< A%d zc=%d, pred %.6f, mantis %ld, ti %ld, true %.6f.\n",
						m_iaAtomSort[zi],zc,tpred[zi],ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_iRelCoord[zc],
						ti,ofr0->m_oaAtoms[m_iaAtomSort[zi]]->m_fRelCoord[zc]);
			#endif
		}
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("*** IntegerArrayToAtoms finished, %d symbols read ***\n",ii-i);
		m_IF.printf("\n");
	}
}



void CBQBEngine::CompressString(const char *s, CBQBBitSet *bs, CBQBStatistics *stat) {

	CBQBIntegerEngine ie(m_IF);
	std::vector<int> ia;
	int z, mi, ma, bits, co;
	CBQBBitSet bs2(m_IF);
	unsigned long tpos;


	tpos = bs->GetLength();

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CompressString >>>\n");

	for (z=0;z<(int)strlen(s);z++)
		ia.push_back(s[z]);

	mi = 255;
	ma = 0;
	for (z=0;z<(int)ia.size();z++) {
		if (ia[z] > ma)
			ma = ia[z];
		if (ia[z] < mi)
			mi = ia[z];
	}
	if (ia.size() == 0) {
		mi = 0;
		ma = 0;
	}
	bits = (int)ceil(mylog2(ma-mi+1));
	co = (int)ceil((bits*ia.size()+4.0)/8.0)+2;

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("CompressString: Trying conventional storage, length %lu (%lu), range %d - %d, %d bits.\n",
			(unsigned long)strlen(s),(unsigned long)ia.size(),mi,ma,bits);

	if (stat != NULL)
		stat->PushStatistics();

	if (ia.size() > 1)
		ie.CompressSingle(
			ia,
			&bs2,
			true,
			true,
			true,
			50,
			1,
			false,
			true,
			10,
			stat
		);

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("CompressString: Source %lu Bytes, IntegerEngine %d Bytes, Conventional %d Bytes.\n",
			(unsigned long)ia.size(),bs2.GetByteLength(),co);

	if ((ia.size() <= 1) || ((co <= bs2.GetByteLength()) && (ia.size() < 256))) {
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("--> Using conventional storage.\n");
		bs->WriteBit(0);
		bs->WriteBits((unsigned long)strlen(s),8);
		bs->WriteBits(mi,8);
		bs->WriteBits(bits,4);
		for (z=0;z<(int)ia.size();z++)
			bs->WriteBits(ia[z]-mi,bits);
	} else {
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("--> Using IntegerEngine storage.\n");
		bs->WriteBit(1);
		bs->WriteBits(&bs2);
	}

	if (stat != NULL) {
		stat->PopStatistics();
		stat->m_oStat.m_lOverhead += bs->GetLength() - tpos;
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("<<< CompressString <<<\n");
}


void CBQBEngine::DecompressString(CBQBBitSet *bs, char **s) {

	CBQBIntegerEngine ie(m_IF);
	std::vector<int> ia;
	int z, mi, c, bits;

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> DecompressString >>>\n");

	if (!bs->ReadBit()) {

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Classical storage.\n");

		c = bs->ReadBitsInteger(8);
		mi = bs->ReadBitsInteger(8);
		bits = bs->ReadBitsInteger(4);

		for (z=0;z<c;z++)
			ia.push_back(mi+bs->ReadBitsInteger(bits));

	} else {

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("IntegerEngine was used.\n");

		ie.DecompressSingle(
			bs,
			ia,
			NULL
		);
	}

	*s = new char[ia.size()+1];

	for (z=0;z<(int)ia.size();z++)
		(*s)[z] = (char)ia[z];

	(*s)[z] = 0;

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("<<< DecompressString <<<\n");
}


bool CBQBEngine::CompressCubeFrame(
		CBQBBitSet *bs,
		int corder,
		CBQBParameterSet_Volumetric *parm,
		int &histused,
		bool staticinfo,
		CBQBStatistics *stat,
		bool skipcomp,
		double *resi,
		std::vector<std::vector<int> > *tciaa
	) {

	CBQBIntegerEngine ie(m_IF);
	const CBQBCubeFrame *cfr;
	unsigned long tpos;


	histused = 0;
	
	if (parm == NULL) {
		m_IF.eprintf("CBQBEngine::CompressCubeFrame(): Error: parm == NULL.\n");
		return false;
	}

	if (resi != NULL)
		*resi = 0.0;

	m_iaTempCompressCube.clear();

	cfr = m_pReadCacheCube->GetFrameHistory(0);

	if (parm->GetHilbert() && ((int)m_iaHilbertIdx.size() != cfr->m_iResXYZ))
		BuildHilbertIdx(cfr->m_iRes[0],cfr->m_iRes[1],cfr->m_iRes[2]/*,verbose*/);

	if (tciaa != NULL) {
		tciaa->push_back(std::vector<int>());
		CubeToIntegerArray(
			tciaa->back(),
			corder,
			parm,
			histused,
			skipcomp,
			resi
		);
	} else
		CubeToIntegerArray(
			m_iaTempCompressCube,
			corder,
			parm,
			histused,
			skipcomp,
			resi
		);

	if (skipcomp || (tciaa != NULL))
		return true;

	tpos = bs->GetLength();

	// Current Version 2

	// History Frames Used
	bs->WriteBits(histused,4);

	if (staticinfo) {

		bs->WriteBit(1);

		// Compression Settings Version Info
		bs->WriteBits(0,3);

		// Compression Settings
		parm->ExportBits(bs);

		ExportCubeHeader(bs,corder/*,verbose*/);

		// Reserved for future use, must be 0
		bs->WriteBit(0);

	} else
		bs->WriteBit(0);

	// Reserved for future use, must be 0
	bs->WriteBit(0);


/*
	// Legacy Version 0 and 1
	bs->WriteBits(corder,4);
	bs->WriteBits(parm->GetSigni(),4);
	bs->WriteBits(parm->GetEps(),6);
	bs->WriteBits(parm->GetSplit(),6);
	bs->WriteBit(parm->GetHilbert()?1:0);

	// Additional fields in version 1
	if (m_pExtrapolator != NULL) {
		bs->WriteBit(1);
		ExportExtrapolatorSettings(bs,true,verbose);
	} else {
		bs->WriteBit(0);
		bs->WriteBits((unsigned long)floor(parm->GetNbhFac()*1000.0+0.5),11);
	}

	bs->WriteBit(0); // Reserved for future use, must be 0
	// End Additional fields
*/

	if (stat != NULL)
		stat->m_oStat.m_lOverhead += bs->GetLength() - tpos;

	ie.Compress(
		m_iaTempCompressCube,
		bs,
		parm->GetCompressorParameterSet(),
		false,
		stat
	);

	return true;
}


bool CBQBEngine::DecompressCubeFrame(
		CBQBBitSet *bs,
		int ver,
		CBQBParameterSet_Volumetric *parm,
		int *histused,
		bool second,
		bool check
	) {

	CBQBIntegerEngine ie(m_IF);
	int tco, csver;
	CBQBCubeFrame *cfr;
	bool statinfo;


	statinfo = false;

	if (ver > BQB_FRAMETYPE_COMPCUBE_VERSION) {
		m_IF.eprintf("CBQBEngine::DecompressCubeFrame(): Unknown Cube Frame Version: %d.\n",ver);
		m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
		m_IF.printf("\n");
		return false;
	}

	m_iaTempCompressCube.clear();

	if (second)
		cfr = GetOutput2CubeFrame(0);
	else
		cfr = GetOutputCubeFrame(0);

	if (ver >= 2) { // Current Version 2

		// History Frames Used
		tco = bs->ReadBitsInteger(4);

		if (histused != NULL)
			*histused = tco;

		if (bs->ReadBit()) { // Static info present

			statinfo = true;

			if (!check)
				m_IF.printf("      Reading volumetric trajectory static information and compression settings...\n");

			// Compression Settings Version Info
			csver = bs->ReadBitsInteger(3);

			// Compression Settings
			parm->ImportBits(bs,csver);

			ImportCubeHeader(bs,second);

			// Reserved for future use, must be 0
			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressCubeFrame(): Reserved bit A expected to be 0, but was 1.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}

		} else {

			if (second)
				cfr->CopyHeader( GetOutput2CubeFrame(1) );
			else
				cfr->CopyHeader( GetOutputCubeFrame(1) );
		}

		// Reserved for future use, must be 0
		if (bs->ReadBit()) {
			m_IF.eprintf("CBQBEngine::DecompressCubeFrame(): Reserved bit B expected to be 0, but was 1.\n");
			m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
			m_IF.printf("\n");
			return false;
		}

	} else { // Legacy Versions 0 and 1

		tco = bs->ReadBitsInteger(4);
		parm->SetSigni(bs->ReadBitsInteger(4));
		parm->SetEps(bs->ReadBitsInteger(6));
		parm->SetSplit(bs->ReadBitsInteger(6));
		parm->SetHilbert(bs->ReadBit());

		if (histused != NULL)
			*histused = tco;

		// Additional fields in version 1
		if (ver == 1) {
			if (bs->ReadBit()) {
				ImportExtrapolatorSettings(bs,true);
			} else {
				parm->SetNbhFac(bs->ReadBitsInteger(11) / 1000.0);
				parm->SetUseExtra(false);
			}

			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressCubeFrame(): Reserved bit expected to be 0, but was 1.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}
		} else {
			parm->SetNbhFac(1.075);
			parm->SetUseExtra(false);
		}
		// End Additional fields

		if (m_pExtrapolator != NULL)
			if (tco > m_pExtrapolator->m_iTRange-1)
				tco = m_pExtrapolator->m_iTRange-1;
	}


	cfr->m_iEps = parm->GetEps();
	cfr->m_iSigni = parm->GetSigni();


	if (ver < 2)
		ImportCubeHeader(bs,second);


	if (!check && parm->GetUseExtra() && (ver >= 2) && statinfo) {

		if (parm->GetExtraPredCorr()) {

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("      Initializing Volumetric Predictor Extrapolator (%d/%d)...\n",
					parm->GetExtraTRange(),parm->GetExtraTOrder());

			if (m_pExtrapolator != NULL)
				delete m_pExtrapolator;
			m_pExtrapolator = new CBQBExtrapolator(m_IF);
			m_pExtrapolator->Initialize(
				cfr->m_iRes[0], // resx,
				cfr->m_iRes[1], // resy,
				cfr->m_iRes[2], // resz,
				1,              // srangex,
				1,              // srangey,
				1,              // srangez,
				parm->GetExtraTRange(),   // trange,
				0,              // sorder,
				parm->GetExtraTOrder(),   // torder,
				0,              // offsetx,
				0,              // offsety,
				0,              // offsetz,
				false,          // crosss,
				false,          // crosst,
				false,          // wrap,
				false,          // crossranges,
				false,          // crossranget
				0,              // distexpo
				parm->GetExtraTimeExpo(), // timeexpo
				false
			);

			if (m_IF.IsPL(BQB_PL_STANDARD)) {
				if (parm->IsExtraSRangeXYZEqual())
					m_IF.printf("      Initializing Volumetric Corrector Extrapolator (%d/%d)...\n",
						parm->GetExtraSRangeX(),parm->GetExtraSOrder());
				else
					m_IF.printf("      Initializing Volumetric Corrector Extrapolator (%d|%d|%d/%d)...\n",
						parm->GetExtraSRangeX(),parm->GetExtraSRangeY(),parm->GetExtraSRangeZ(),parm->GetExtraSOrder());
			}

			if (m_pExtrapolatorCorr != NULL)
				delete m_pExtrapolatorCorr;
			m_pExtrapolatorCorr = new CBQBExtrapolator(m_IF);
			m_pExtrapolatorCorr->Initialize(
				cfr->m_iRes[0],     // resx,
				cfr->m_iRes[1],     // resy,
				cfr->m_iRes[2],     // resz,
				parm->GetExtraSRangeX(),      // srangex,
				parm->GetExtraSRangeY(),      // srangey,
				parm->GetExtraSRangeZ(),      // srangez,
				1,               // trange,
				parm->GetExtraSOrder(),       // sorder,
				0,               // torder,
				parm->GetExtraOffsetX(),      // offsetx,
				parm->GetExtraOffsetY(),      // offsety,
				parm->GetExtraOffsetZ(),      // offsetz,
				parm->GetExtraCrossS(),       // crosss,
				false,           // crosst,
				parm->GetExtraWrap(),         // wrap,
				parm->GetExtraCrossRangeS(),  // crossranges,
				false,           // crossranget
				parm->GetExtraDistExpo(),     // distexpo
				0,                  // timeexpo
				false
			);

		} else {

			if (m_IF.IsPL(BQB_PL_STANDARD)) {

				if (parm->IsExtraSRangeXYZEqual())
					m_IF.printf("      Initializing Volumetric Extrapolator (%d/%d;%d/%d)...\n",
						parm->GetExtraSRangeX(),
						parm->GetExtraSOrder(),
						parm->GetExtraTRange(),
						parm->GetExtraTOrder()
					);
				else
					m_IF.printf("      Initializing Volumetric Extrapolator (%d|%d|%d/%d;%d|%d)...\n",
						parm->GetExtraSRangeX(),
						parm->GetExtraSRangeY(),
						parm->GetExtraSRangeZ(),
						parm->GetExtraSOrder(),
						parm->GetExtraTRange(),
						parm->GetExtraTOrder()
					);
			}

			if (m_pExtrapolator != NULL)
				delete m_pExtrapolator;
			m_pExtrapolator = new CBQBExtrapolator(m_IF);
			m_pExtrapolator->Initialize(
				cfr->m_iRes[0],     // resx,
				cfr->m_iRes[1],     // resy,
				cfr->m_iRes[2],     // resz,
				parm->GetExtraSRangeX(),      // srangex,
				parm->GetExtraSRangeY(),      // srangey,
				parm->GetExtraSRangeZ(),      // srangez,
				parm->GetExtraTRange(),       // trange,
				parm->GetExtraSOrder(),       // sorder,
				parm->GetExtraTOrder(),       // torder,
				parm->GetExtraOffsetX(),      // offsetx,
				parm->GetExtraOffsetY(),      // offsety,
				parm->GetExtraOffsetZ(),      // offsetz,
				parm->GetExtraCrossS(),       // crosss,
				parm->GetExtraCrossT(),       // crosst,
				parm->GetExtraWrap(),         // wrap,
				parm->GetExtraCrossRangeS(),  // crossranges,
				parm->GetExtraCrossRangeT(),  // crossranget
				parm->GetExtraDistExpo(),     // distexpo
				parm->GetExtraTimeExpo(),     // timeexpo
				false
			);
		}
	}


	if (parm->GetHilbert() && ((int)m_iaHilbertIdx.size() != cfr->m_iResXYZ))
		BuildHilbertIdx(cfr->m_iRes[0],cfr->m_iRes[1],cfr->m_iRes[2]);

	if (!check && (m_pExtrapolator != NULL))
		m_pExtrapolator->PushCubeFrame(cfr);

	ie.Decompress(
		bs,
		m_iaTempCompressCube,
		parm->GetCompressorParameterSet()
	);

	IntegerArrayToCube(
		m_iaTempCompressCube,
		tco,
		parm,
		ver,
		second
	);

	return true;
}


void CBQBEngine::ExportAtomLabels(const CBQBAtomSet *as, CBQBBitSet *bs, CBQBStatistics *stat) {

	//std::vector<int> ia;
	char tch[5001];
	std::string tst;
	int z, z2, zt, rem;


	if ((as->m_oaAtoms.size()%5000) == 0)
		zt = (int)as->m_oaAtoms.size()/5000;
	else
		zt = (int)as->m_oaAtoms.size()/5000 + 1;

	if (zt > 255) {
		m_IF.eprintf("CBQBEngine::ExportAtomLabels(): Error: More than 1'275'000 atoms are not supported (%d).\n",
			zt*5000);
		abort();
	}

	bs->WriteBits(zt,8);

	for (z=0;z<zt;z++) {
	
		rem = (int)as->m_oaAtoms.size() - z*5000;
		if (rem > 5000)
			rem = 5000;

		tst = "";
		for (z2=0;z2<rem;z2++) {
			tch[z2] = (char)as->m_oaAtoms[z*5000+z2]->m_sLabel.length();
			tst += as->m_oaAtoms[z*5000+z2]->m_sLabel;
			if (as->m_oaAtoms[z*5000+z2]->m_sLabel.length() > 127) {
				m_IF.eprintf("CBQBEngine::ExportAtomLabels(): Error: Labels longer than 127 characters are not supported (%lu).\n",
					(unsigned long)as->m_oaAtoms[z*5000+z2]->m_sLabel.length());
				abort();
			}
		}
		tch[rem] = 0;

		CompressString(tch,bs,stat);
		CompressString(tst.c_str(),bs,stat);
	}
}


void CBQBEngine::ImportAtomLabels(CBQBAtomSet *as, CBQBBitSet *bs) {

	int z, z2, zt, i;
	char *tch, *tch2, *p, buf[128];
	unsigned char tuc;


	zt = bs->ReadBitsInteger(8);

	i = 0;
	for (z=0;z<zt;z++) {

		DecompressString(bs,&tch);
		DecompressString(bs,&tch2);

		z2 = 0;
		p = tch2;
		while (tch[z2] != 0) {
			if (tch[z2] < 0) {
				m_IF.eprintf("CBQBEngine::ImportAtomLabels(): Error: Label length < 0 (%d).\n",
					tch[z2]);
				abort();
			}
			tuc = (unsigned char)tch[z2];
			memcpy(buf,p,tuc);
			buf[tuc] = 0;
			as->m_oaAtoms[i]->m_sLabel = buf;
			p += tuc;
			z2++;
			i++;
		}

		delete[] tch;
		delete[] tch2;
	}

	if (i != (int)as->m_oaAtoms.size()) {
		m_IF.eprintf("CBQBEngine::ImportAtomLabels(): Error: Expected %lu atom labels, but found %d (%d).\n",
			(unsigned long)as->m_oaAtoms.size(),i,zt);
		abort();
	}
}


bool CBQBEngine::CompressAtomFrame(
		CBQBBitSet *bs,
		int order,
		CBQBParameterSet_Position *parm,
		int &histused,
		bool staticinfo,
		bool comment,
		CBQBStatistics *stat,
		bool skipcomp,
		double *resi,
		std::vector<std::vector<int> > *tciaa
	) {

	CBQBIntegerEngine ie(m_IF);
	//std::string ts;
	int bits, z;
	const CBQBAtomSet *as;
	unsigned long tpos;


	histused = 0;

	if (parm == NULL) {
		m_IF.eprintf("CBQBEngine::CompressAtomFrame(): Error: parm == NULL.\n");
		return false;
	}

	if (resi != NULL)
		*resi = 0.0;

	m_iaTempCompressXYZ.clear();

	if (skipcomp) {
		AtomsToIntegerArray(
			m_iaTempCompressXYZ,
			order,
			parm,
			histused,
			skipcomp,
			resi
		);
		return true;
	}

	if (m_pReadCacheCube != NULL)
		as = m_pReadCacheCube->GetFrameHistory(0)->m_pAtoms;
	else
		as = m_pReadCacheXYZ->GetFrameHistory(0);


	if (parm->GetSortAtom())
		BuildAtomSort(as);
	else
		BuildIdentityAtomSort(as);


	if (tciaa != NULL) {
		tciaa->push_back(std::vector<int>());
		AtomsToIntegerArray(
			tciaa->back(),
			order,
			parm,
			histused,
			skipcomp,
			resi
		);
	} else
		AtomsToIntegerArray(
			m_iaTempCompressXYZ,
			order,
			parm,
			histused,
			skipcomp,
			resi
		);


	tpos = bs->GetLength();

	if (stat != NULL)
		stat->PushStatistics();

	// New in Version 2: Atom Count always written
	bits = (int)ceil(mylog2((double)as->m_oaAtoms.size()+1));
	bs->WriteBits(bits,5);
	bs->WriteBits((unsigned long)as->m_oaAtoms.size(),bits);

	if (staticinfo) {

		// Compression settings stored
		bs->WriteBit(1);

		// 3 Bits Compression Setting Version Info (Currently Version 0)
		bs->WriteBits(0,3);

		// Compression Settings
		parm->ExportBits(bs);

	} else
		bs->WriteBit(0); // No compression settings stored


	if (staticinfo) {

		bs->WriteBit(1); // Static trajectory information present

		if (as->m_bOrd) {
			bs->WriteBit(1); // Ord numbers
			for (z=0;z<(int)as->m_oaAtoms.size();z++)
				bs->WriteBits(as->m_oaAtoms[z]->m_iOrd,8);
		} else
			bs->WriteBit(0); // No Ord numbers
		
		if (as->m_bLabels) {
			bs->WriteBit(1); // Atom label strings
			ExportAtomLabels(as,bs,stat);
		} else
			bs->WriteBit(0); // No Atom label strings

		// Additional fields since version 1
		if (parm->GetSortAtom())
			bs->WriteBit(1);
		else
			bs->WriteBit(0);
		// End Additional fields since version 1

		// Additional Fields in version 2

		// Trajectory Comment
		bs->WriteBit(0);

		// Static Cell Info
		bs->WriteBit(0);

		// Static Topology
		bs->WriteBit(0);

		// Integrator Properties
		bs->WriteBit(0);

		// Atom Masses
		bs->WriteBit(0);

		// Fixed atomic partial charges
		bs->WriteBit(0);

		// Atom Radii
		bs->WriteBit(0);

		// Additional per-atom Datasets: Header
		bs->WriteBit(0);

		// Reserved for future use, must be 0
		bs->WriteBit(0);

		// End Additional Fields in version 2

	} else
		bs->WriteBit(0); // No static trajectory information present

	// Comment
	if ((as->m_sComment != NULL) && comment) {
		bs->WriteBit(1);
		CompressString(as->m_sComment,bs,stat);
	} else
		bs->WriteBit(0);

	// Additional Fields in version 2

	// History Frames Used
	bs->WriteBits(histused,4);

	// Simulation Time
	bs->WriteBit(0);

	// Thermodynamic Properties
	bs->WriteBit(0);

	// Dynamic Cell Info
	bs->WriteBit(0);

	// Dynamic Topology
	bs->WriteBit(0);

	// Voronoi Tessellation
	bs->WriteBit(0);

	// Additional per-atom Datasets: Data
	bs->WriteBit(0);

	// Reserved for future use, must be 0
	bs->WriteBit(0);

	// End Additional Fields in version 2

/*
	// Legacy Version 0 / 1
	bs->WriteBits(order,4);
	bs->WriteBits(parm->GetPrecision(),4);
	bs->WriteBits(parm->GetSplit(),6);

	// Additional fields in version 1
	if (m_pExtrapolatorXYZ != NULL) {
		bs->WriteBit(1);
		ExportExtrapolatorSettings(bs,false,verbose);
	} else
		bs->WriteBit(0);

	bs->WriteBit(0); // Reserved for future use, must be 0
	// End Additional fields
*/


	if (stat != NULL) {
		stat->PopStatistics();
		stat->m_oStat.m_lOverhead += bs->GetLength() - tpos;
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("CBQBEngine::CompressAtomFrame(): Will compress %lu symbols.\n",(unsigned long)m_iaTempCompressXYZ.size());

	if (tciaa == NULL)
		ie.Compress(
			m_iaTempCompressXYZ,
			bs,
			parm->GetCompressorParameterSet(),
			false,
			stat
		);

	return true;
}


bool CBQBEngine::DecompressAtomStep(
		const std::vector<unsigned char> *inp,
		int ver,
		CBQBParameterSet_Position *parm,
		int *histused,
		bool second,
		bool check
	) {

	CBQBBitSet bs(m_IF);
	CBQBAtomSet *afr;
	bool b;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBEngine::DecompressAtomStep() >>>\n");

	afr = new CBQBAtomSet(m_IF);
	PushOutputAtomFrame(afr);

	bs.m_iaData.assign(inp->begin(),inp->end());

	b =  DecompressAtomFrame(
		&bs,
		ver,
		parm,
		histused,
		second,
		check
	);

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("<<< CBQBEngine::DecompressAtomStep() <<<\n");

	return b;
}


bool CBQBEngine::DecompressAtomFrame(
		CBQBBitSet *bs,
		int ver,
		CBQBParameterSet_Position *parm,
		int *histused,
		bool second,
		bool check
	) {

	CBQBIntegerEngine ie(m_IF);
	int z, bits, atc, csver;
	int tao;
	CBQBAtom *at;
	CBQBAtomSet *as, *as2;
	char *tch, *p, *q;


	atc = 0;

	if (ver > BQB_FRAMETYPE_COMPTRAJ_VERSION) {
		m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Unknown Atom Frame Version: %d.\n",ver);
		m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
		m_IF.printf("\n");
		return false;
	}

	m_iaTempCompressXYZ.clear();

	if (second)
		as = GetOutput2AtomFrame(0);
	else
		as = GetOutputAtomFrame(0);


	if (ver >= 2) {

		// New in Version 2: Atom count is always written
		bits = bs->ReadBitsInteger(5);
		atc = bs->ReadBitsInteger(bits);


		if (bs->ReadBit()) { // Compression settings stored

			//m_IF.printf("Reading compression settings...\n");
			if (!check)
				m_IF.printf("      Reading position trajectory compression settings...\n");

			// Compression Settings Version Info
			csver = bs->ReadBitsInteger(3);

			// Compression Settings
			parm->ImportBits(bs,csver);

			//m_IF.printf("%s",parm->ToString(4).c_str());

			if (!check && parm->GetUseExtra()) {
				m_IF.printf("      Initializing Extrapolator (%d/%d)...\n",
						parm->GetExtraTRange(),
						parm->GetExtraTOrder()
				);
				if (m_pExtrapolatorXYZ != NULL)
					delete m_pExtrapolatorXYZ;
				m_pExtrapolatorXYZ = new CBQBExtrapolator(m_IF);
				m_pExtrapolatorXYZ->InitializeXYZ(parm,false);
			}
		}
	}


	if (bs->ReadBit()) { // Static trajectory information stored

		if (!check)
			m_IF.printf("      Reading position trajectory static information...\n");

		if (ver < 2) {
			bits = bs->ReadBitsInteger(5);
			atc = bs->ReadBitsInteger(bits);
		}

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Will read %d atoms.\n",atc);

		for (z=0;z<atc;z++) {
			at = new CBQBAtom();
			as->m_oaAtoms.push_back(at);
		}

		if (bs->ReadBit()) { // Read Ord Numbers
			as->m_bOrd = true;
			for (z=0;z<atc;z++)
				as->m_oaAtoms[z]->m_iOrd = bs->ReadBitsInteger(8);
		}

		if (bs->ReadBit()) { // Read Atom label strings

			as->m_bLabels = true;

			if (ver >= 1) { // New method

				ImportAtomLabels(as,bs);
				if (!as->m_bOrd)
					for (z=0;z<atc;z++)
						as->m_oaAtoms[z]->m_iOrd = GetAtomOrd(as->m_oaAtoms[z]->m_sLabel.c_str());

			} else { // Old method

				DecompressString(bs,&tch);
				if (m_IF.IsPL(BQB_PL_VERBOSE))
					m_IF.printf("Decompressed string to %lu characters.\n",(unsigned long)strlen(tch));
				p = tch;
				for (z=0;z<atc;z++) {
					q = strchr(p,' ');
					if (q != NULL)
						*q = 0;
					as->m_oaAtoms[z]->m_sLabel = p;
					if (!as->m_bOrd)
						as->m_oaAtoms[z]->m_iOrd = GetAtomOrd(p);
					p = q+1;
				}
				delete[] tch;
			}
		}

		// Additional fields since version 1
		if (ver >= 1) {
			if (bs->ReadBit()) // Use Atom Sorting
				BuildAtomSort(as);
			else
				BuildIdentityAtomSort(as);
		} else
			BuildAtomSort(as);
		// End Additional fields since version 1

		if (ver >= 2) { // Additional fields in version 2

			// Trajectory Comment
			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Trajectory comment found in BQB file, but not yet implemented.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}

			// Static Cell Info
			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Static cell info found in BQB file, but not yet implemented.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}

			// Static Topology
			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Static topology found in BQB file, but not yet implemented.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}

			// Integrator Properties
			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Integrator properties found in BQB file, but not yet implemented.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}

			// Atom Masses
			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Atom masses found in BQB file, but not yet implemented.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}

			// Fixed atomic partial charges
			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Fixed atomic partial charges found in BQB file, but not yet implemented.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}

			// Atom Radii
			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Atom radii found in BQB file, but not yet implemented.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}

			// Additional per-atom Datasets: Header
			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Additional per-atom dataset header found in BQB file, but not yet implemented.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}

			// Reserved for future use, must be 0
			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Reserved bit A expected to be 0, but was 1.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}

		} // End additional fields in version 2

	} else { // No static trajectory information stored

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Atom type information already read before.\n");
		if (second)
			as2 = GetOutput2AtomFrame(1);
		else
			as2 = GetOutputAtomFrame(1);
		if (as2 == NULL) {
			m_IF.printf("CEngine::DecompressAtomFrame(): Error: First frame does not contain atom type information.\n");
			return false;
		}
		for (z=0;z<(int)as2->m_oaAtoms.size();z++) {
			at = new CBQBAtom();
			at->m_iOrd = as2->m_oaAtoms[z]->m_iOrd;
			at->m_sLabel = as2->m_oaAtoms[z]->m_sLabel;
			as->m_oaAtoms.push_back(at);
		}
		as->m_bLabels = as2->m_bLabels;
		as->m_bOrd = as2->m_bOrd;
	}

	if (bs->ReadBit()) // Comment was stored
		DecompressString(
			bs,
			&as->m_sComment
		);
	else
		as->m_sComment = NULL;


	if (ver == 2) { // Current Version 2

		// History Frames Used
		tao = bs->ReadBitsInteger(4);

		if (histused != NULL)
			*histused = tao;

		// Simulation Time
		if (bs->ReadBit()) {
			m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Simulation time found in BQB file, but not yet implemented.\n");
			m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
			m_IF.printf("\n");
			return false;
		}

		// Thermodynamic Properties
		if (bs->ReadBit()) {
			m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Thermodynamic properties found in BQB file, but not yet implemented.\n");
			m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
			m_IF.printf("\n");
			return false;
		}

		// Dynamic Cell Info
		if (bs->ReadBit()) {
			m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Dynamic cell info found in BQB file, but not yet implemented.\n");
			m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
			m_IF.printf("\n");
			return false;
		}

		// Dynamic Topology
		if (bs->ReadBit()) {
			m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Dynamic topology found in BQB file, but not yet implemented.\n");
			m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
			m_IF.printf("\n");
			return false;
		}

		// Voronoi Tessellation
		if (bs->ReadBit()) {
			m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Voronoi tessellation found in BQB file, but not yet implemented.\n");
			m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
			m_IF.printf("\n");
			return false;
		}

		// Additional per-atom Datasets: Data
		if (bs->ReadBit()) {
			m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Additional per-atom dataset data found in BQB file, but not yet implemented.\n");
			m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
			m_IF.printf("\n");
			return false;
		}

		// Reserved for future use, must be 0
		if (bs->ReadBit()) {
			m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Reserved bit B expected to be 0, but was 1.\n");
			m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
			m_IF.printf("\n");
			return false;
		}

	} else { // Legacy versions 0 and 1

		tao = bs->ReadBitsInteger(4);
		parm->SetPrecision(bs->ReadBitsInteger(4));
		parm->SetSplit(bs->ReadBitsInteger(6));

		if (histused != NULL)
			*histused = tao;

		// Additional fields in version 1
		if (ver == 1) {
			if (bs->ReadBit())
				ImportExtrapolatorSettings(bs,false);
			else
				parm->SetUseExtra(false);
			if (bs->ReadBit()) {
				m_IF.eprintf("CBQBEngine::DecompressAtomFrame(): Reserved bit expected to be 0, but was 1.\n");
				m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
				m_IF.printf("\n");
				return false;
			}
		} else
			parm->SetUseExtra(false);
		// End Additional fields

		if (m_pExtrapolatorXYZ != NULL)
			if (tao > m_pExtrapolatorXYZ->m_iTRange-1)
				tao = m_pExtrapolatorXYZ->m_iTRange-1;
	}


	as->m_iSigni = parm->GetPrecision();

	if (!check && (m_pExtrapolatorXYZ != NULL))
		m_pExtrapolatorXYZ->PushAtomFrame(as);

	ie.Decompress(
		bs,
		m_iaTempCompressXYZ,
		parm->GetCompressorParameterSet()
	);

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("CBQBEngine::DecompressAtomFrame(): Decompressed %lu symbols.\n",(unsigned long)m_iaTempCompressXYZ.size());

	IntegerArrayToAtoms(
		m_iaTempCompressXYZ,
		tao,
		parm,
		ver,
		second
	);

	return true;
}


bool CBQBEngine::PrepareDecompressCube() {

	int z;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBEngine::PrepareDecompressCube() >>>\n");

	for (z=0;z<(int)m_oaOutputCubeBuf.size();z++)
		if (m_oaOutputCubeBuf[z] != NULL)
			delete m_oaOutputCubeBuf[z];

	for (z=0;z<(int)m_oaOutputAtomBuf.size();z++)
		if (m_oaOutputAtomBuf[z] != NULL)
			delete m_oaOutputAtomBuf[z];

	m_oaOutputCubeBuf.resize(10);
	for (z=0;z<10;z++)
		m_oaOutputCubeBuf[z] = NULL;
	m_iOutputCubeBufPos = 0;

	m_oaOutputAtomBuf.resize(10);
	for (z=0;z<10;z++)
		m_oaOutputAtomBuf[z] = NULL;
	m_iOutputAtomBufPos = 0;

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("<<< CBQBEngine::PrepareDecompressCube() <<<\n");

	return true;
}


bool CBQBEngine::DecompressCubeStep(
		const std::vector<unsigned char> *inp,
		int ver,
		CBQBParameterSet_PosAndVol *parm,
		int *ahistused,
		int *chistused,
		bool skipvol
	) {

	CBQBCubeFrame *cfr;
	int al, cl, ahi, chi;
	CBQBBitSet bsat(m_IF), bscu(m_IF), bshe(m_IF);


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBEngine::DecompressCubeStep() >>>\n");

	bshe.Clear();
	bshe.m_iaData.assign(inp->begin(),inp->begin()+8);

	al = bshe.ReadBitsInteger(32);
	cl = bshe.ReadBitsInteger(32);

	cfr = new CBQBCubeFrame(m_IF);
	PushOutputCubeFrame(cfr);

	cfr->m_pAtoms = new CBQBAtomSet(m_IF);
	PushOutputAtomFrame(cfr->m_pAtoms);

	bsat.Clear();
	bsat.m_iaData.assign(inp->begin()+8,inp->begin()+8+al);

	bscu.Clear();
	bscu.m_iaData.assign(inp->begin()+8+al,inp->begin()+8+al+cl);

	ahi = 0;

	if (!DecompressAtomFrame(
			&bsat,
			ver,
			parm->GetPositionParameterSet(),
			&ahi,
			false,
			false
		))
		return false;

	if (ahistused != NULL)
		*ahistused = ahi;

	if (!skipvol) {

		chi = 0;

		if (!DecompressCubeFrame(
				&bscu,
				ver,
				parm->GetVolumetricParameterSet(),
				&chi,
				false,
				false
			))
			return false;

		if (chistused != NULL)
			*chistused = chi;
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("<<< CBQBEngine::DecompressCubeStep() <<<\n");

	return true;
}


void CBQBEngine::PushOutputCubeFrame(CBQBCubeFrame *cfr) {

	m_iOutputCubeBufPos++;
	if (m_iOutputCubeBufPos == (int)m_oaOutputCubeBuf.size())
		m_iOutputCubeBufPos = 0;

	if (m_bDeleteOutFrames)
		if (m_oaOutputCubeBuf[m_iOutputCubeBufPos] != NULL)
			delete m_oaOutputCubeBuf[m_iOutputCubeBufPos];
	m_oaOutputCubeBuf[m_iOutputCubeBufPos] = cfr;
}


void CBQBEngine::PushOutput2CubeFrame(CBQBCubeFrame *cfr) {

	m_iOutput2CubeBufPos++;
	if (m_iOutput2CubeBufPos == (int)m_oaOutput2CubeBuf.size())
		m_iOutput2CubeBufPos = 0;

	if (m_bDeleteOutFrames)
		if (m_oaOutput2CubeBuf[m_iOutput2CubeBufPos] != NULL)
			delete m_oaOutput2CubeBuf[m_iOutput2CubeBufPos];
	m_oaOutput2CubeBuf[m_iOutput2CubeBufPos] = cfr;
}


CBQBCubeFrame* CBQBEngine::GetOutputCubeFrame(int depth) {

	int ti;

	if (depth >= (int)m_oaOutputCubeBuf.size()) {
		m_IF.eprintf("CBQBEngine::GetOutputCubeFrame(): Out of bounds (%d/%d).\n",
			depth,(int)m_oaOutputCubeBuf.size());
		abort();
	}

	ti = m_iOutputCubeBufPos - depth;
	if (ti < 0)
		ti += (int)m_oaOutputCubeBuf.size();

	return m_oaOutputCubeBuf[ti];
}


CBQBCubeFrame* CBQBEngine::GetOutput2CubeFrame(int depth) {

	int ti;

	if (depth >= (int)m_oaOutput2CubeBuf.size()) {
		m_IF.eprintf("CBQBEngine::GetOutput2CubeFrame(): Out of bounds (%d/%d).\n",
			depth,(int)m_oaOutput2CubeBuf.size());
		abort();
	}

	ti = m_iOutput2CubeBufPos - depth;
	if (ti < 0)
		ti += (int)m_oaOutput2CubeBuf.size();

	return m_oaOutput2CubeBuf[ti];
}


void CBQBEngine::PushOutputAtomFrame(CBQBAtomSet *cfr) {

	m_iOutputAtomBufPos++;
	if (m_iOutputAtomBufPos == (int)m_oaOutputAtomBuf.size())
		m_iOutputAtomBufPos = 0;

	if (m_bDeleteOutFrames)
		if (m_oaOutputAtomBuf[m_iOutputAtomBufPos] != NULL)
			delete m_oaOutputAtomBuf[m_iOutputAtomBufPos];
	m_oaOutputAtomBuf[m_iOutputAtomBufPos] = cfr;
}


void CBQBEngine::PushOutput2AtomFrame(CBQBAtomSet *cfr) {

	m_iOutput2AtomBufPos++;
	if (m_iOutput2AtomBufPos == (int)m_oaOutput2AtomBuf.size())
		m_iOutput2AtomBufPos = 0;

	if (m_bDeleteOutFrames)
		if (m_oaOutput2AtomBuf[m_iOutput2AtomBufPos] != NULL)
			delete m_oaOutput2AtomBuf[m_iOutput2AtomBufPos];
	m_oaOutput2AtomBuf[m_iOutput2AtomBufPos] = cfr;
}


CBQBAtomSet* CBQBEngine::GetOutputAtomFrame(int depth) {

	int ti;

	if (depth >= (int)m_oaOutputAtomBuf.size()) {
		m_IF.eprintf("CBQBEngine::GetOutputAtomFrame(): Out of bounds (%d/%d).\n",
			depth,(int)m_oaOutputAtomBuf.size());
		abort();
	}

	ti = m_iOutputAtomBufPos - depth;
	if (ti < 0)
		ti += (int)m_oaOutputAtomBuf.size();

	return m_oaOutputAtomBuf[ti];
}


CBQBAtomSet* CBQBEngine::GetOutput2AtomFrame(int depth) {

	int ti;

	if (depth >= (int)m_oaOutput2AtomBuf.size()) {
		m_IF.eprintf("CBQBEngine::GetOutput2AtomFrame(): Out of bounds (%d/%d).\n",
			depth,(int)m_oaOutput2AtomBuf.size());
		abort();
	}

	ti = m_iOutput2AtomBufPos - depth;
	if (ti < 0)
		ti += (int)m_oaOutput2AtomBuf.size();

	return m_oaOutput2AtomBuf[ti];
}


void CBQBEngine::BuildAtomSort(const CBQBAtomSet *as) {

	int z;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("CBQBEngine::BuildAtomSort() running...\n");

	m_iaAtomSort.resize(as->m_oaAtoms.size());
	m_iaAtomOrd.resize(as->m_oaAtoms.size());
	for (z=0;z<(int)m_iaAtomSort.size();z++) {
		m_iaAtomSort[z] = z;
		m_iaAtomOrd[z] = as->m_oaAtoms[z]->m_iOrd;
	}

	std::stable_sort(
		m_iaAtomSort.begin(),
		m_iaAtomSort.end(),
		SSortAtomOrd(*this)
	);
}


void CBQBEngine::BuildIdentityAtomSort(const CBQBAtomSet *as) {

	int z;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("CBQBEngine::BuildIdentityAtomSort() running...\n");

	m_iaAtomSort.resize(as->m_oaAtoms.size());
	m_iaAtomOrd.resize(as->m_oaAtoms.size());
	for (z=0;z<(int)m_iaAtomSort.size();z++) {
		m_iaAtomSort[z] = z;
		m_iaAtomOrd[z] = as->m_oaAtoms[z]->m_iOrd;
	}
}


void CBQBEngine::ExportExtrapolatorSettings(CBQBBitSet *bs, bool volumetric) {

	if (volumetric) {

		// Use Volumetric Extrapolator?
		if (m_pExtrapolator == NULL) {
			bs->WriteBit(0);
			return;
		}
		bs->WriteBit(1);

		// ResX
		if ((m_pExtrapolator->m_iRes[0] < 1) || (m_pExtrapolator->m_iRes[0] > 1023)) {
			m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.ResX out of range (%d, allowed 1..1023).\n",
				m_pExtrapolator->m_iRes[0]);
			abort();
		}
		bs->WriteBits(m_pExtrapolator->m_iRes[0],10);

		// ResY
		if ((m_pExtrapolator->m_iRes[1] < 1) || (m_pExtrapolator->m_iRes[1] > 1023)) {
			m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.ResY out of range (%d, allowed 1..1023).\n",
				m_pExtrapolator->m_iRes[1]);
			abort();
		}
		bs->WriteBits(m_pExtrapolator->m_iRes[1],10);

		// ResZ
		if ((m_pExtrapolator->m_iRes[2] < 1) || (m_pExtrapolator->m_iRes[2] > 1023)) {
			m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.ResZ out of range (%d, allowed 1..1023).\n",
				m_pExtrapolator->m_iRes[2]);
			abort();
		}
		bs->WriteBits(m_pExtrapolator->m_iRes[2],10);

		// PredCorr
		if (m_pExtrapolatorCorr != NULL)
			bs->WriteBit(1);
		else
			bs->WriteBit(0);

		if (m_pExtrapolatorCorr != NULL) {

			// SRangeX
			if ((m_pExtrapolatorCorr->m_iSRange[0] < 1) || (m_pExtrapolatorCorr->m_iSRange[0] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorCorr.SRangeX out of range (%d, allowed 1..15).\n",
					m_pExtrapolatorCorr->m_iSRange[0]);
				abort();
			}
			bs->WriteBits(m_pExtrapolatorCorr->m_iSRange[0],4);

			// SRangeY
			if ((m_pExtrapolatorCorr->m_iSRange[1] < 1) || (m_pExtrapolatorCorr->m_iSRange[1] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorCorr.SRangeY out of range (%d, allowed 1..15).\n",
					m_pExtrapolatorCorr->m_iSRange[1]);
				abort();
			}
			bs->WriteBits(m_pExtrapolatorCorr->m_iSRange[1],4);

			// SRangeZ
			if ((m_pExtrapolatorCorr->m_iSRange[2] < 1) || (m_pExtrapolatorCorr->m_iSRange[2] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorCorr.SRangeZ out of range (%d, allowed 1..15).\n",
					m_pExtrapolatorCorr->m_iSRange[2]);
				abort();
			}
			bs->WriteBits(m_pExtrapolatorCorr->m_iSRange[2],4);

			// SOrder
			if ((m_pExtrapolatorCorr->m_iSOrder < 0) || (m_pExtrapolatorCorr->m_iSOrder > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorCorr.SOrder out of range (%d, allowed 0..15).\n",
					m_pExtrapolatorCorr->m_iSOrder);
				abort();
			}
			bs->WriteBits(m_pExtrapolatorCorr->m_iSOrder,4);

			// OffsetX
			if ((m_pExtrapolatorCorr->m_iSOffset[0] < 0) || (m_pExtrapolatorCorr->m_iSOffset[0] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorCorr.SOffsetX out of range (%d, allowed 0..15).\n",
					m_pExtrapolatorCorr->m_iSOffset[0]);
				abort();
			}
			if (m_pExtrapolatorCorr->m_iSOffset[0] >= m_pExtrapolatorCorr->m_iSRange[0]) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorCorr.SOffsetX (%d) >= SRangeX (%d.\n",
					m_pExtrapolatorCorr->m_iSOffset[0],m_pExtrapolatorCorr->m_iSRange[0]);
				abort();
			}
			bs->WriteBits(m_pExtrapolatorCorr->m_iSOffset[0],4);

			// OffsetY
			if ((m_pExtrapolatorCorr->m_iSOffset[1] < 0) || (m_pExtrapolatorCorr->m_iSOffset[1] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorCorr.SOffsetY out of range (%d, allowed 0..15).\n",
					m_pExtrapolatorCorr->m_iSOffset[1]);
				abort();
			}
			if (m_pExtrapolatorCorr->m_iSOffset[1] >= m_pExtrapolatorCorr->m_iSRange[1]) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorCorr.SOffsetY (%d) >= SRangeY (%d.\n",
					m_pExtrapolatorCorr->m_iSOffset[1],m_pExtrapolatorCorr->m_iSRange[1]);
				abort();
			}
			bs->WriteBits(m_pExtrapolatorCorr->m_iSOffset[1],4);

			// OffsetZ
			if ((m_pExtrapolatorCorr->m_iSOffset[2] < 0) || (m_pExtrapolatorCorr->m_iSOffset[2] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorCorr.SOffsetZ out of range (%d, allowed 0..15).\n",
					m_pExtrapolatorCorr->m_iSOffset[2]);
				abort();
			}
			if (m_pExtrapolatorCorr->m_iSOffset[2] >= m_pExtrapolatorCorr->m_iSRange[2]) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorCorr.SOffsetZ (%d) >= SRangeZ (%d.\n",
					m_pExtrapolatorCorr->m_iSOffset[2],m_pExtrapolatorCorr->m_iSRange[2]);
				abort();
			}
			bs->WriteBits(m_pExtrapolatorCorr->m_iSOffset[2],4);

			// CrossS
			if (m_pExtrapolatorCorr->m_bCrossS)
				bs->WriteBit(1);
			else
				bs->WriteBit(0);

			// CrossRangeS
			if (m_pExtrapolatorCorr->m_bCrossRangeS)
				bs->WriteBit(1);
			else
				bs->WriteBit(0);

			// Wrap
			if (m_pExtrapolatorCorr->m_bWrap)
				bs->WriteBit(1);
			else
				bs->WriteBit(0);

			// DistExpo
			if ((m_pExtrapolatorCorr->m_fDistExpo < 0) || (m_pExtrapolatorCorr->m_fDistExpo > 16.0)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorCorr.DistExpo out of range (%.3f, allowed 0..16).\n",
					m_pExtrapolatorCorr->m_fDistExpo);
				abort();
			}
			bs->WriteBits((unsigned int)floor(m_pExtrapolatorCorr->m_fDistExpo*1000.0+0.5),14);

			// CorrFactor
			bs->WriteBits(0,11);

		} else {

			// SRangeX
			if ((m_pExtrapolator->m_iSRange[0] < 1) || (m_pExtrapolator->m_iSRange[0] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.SRangeX out of range (%d, allowed 1..15).\n",
					m_pExtrapolator->m_iSRange[0]);
				abort();
			}
			bs->WriteBits(m_pExtrapolator->m_iSRange[0],4);

			// SRangeY
			if ((m_pExtrapolator->m_iSRange[1] < 1) || (m_pExtrapolator->m_iSRange[1] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.SRangeY out of range (%d, allowed 1..15).\n",
					m_pExtrapolator->m_iSRange[1]);
				abort();
			}
			bs->WriteBits(m_pExtrapolator->m_iSRange[1],4);

			// SRangeZ
			if ((m_pExtrapolator->m_iSRange[2] < 1) || (m_pExtrapolator->m_iSRange[2] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.SRangeZ out of range (%d, allowed 1..15).\n",
					m_pExtrapolator->m_iSRange[2]);
				abort();
			}
			bs->WriteBits(m_pExtrapolator->m_iSRange[2],4);

			// SOrder
			if ((m_pExtrapolator->m_iSOrder < 0) || (m_pExtrapolator->m_iSOrder > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.SOrder out of range (%d, allowed 0..15).\n",
					m_pExtrapolator->m_iSOrder);
				abort();
			}
			bs->WriteBits(m_pExtrapolator->m_iSOrder,4);

			// OffsetX
			if ((m_pExtrapolator->m_iSOffset[0] < 0) || (m_pExtrapolator->m_iSOffset[0] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.SOffsetX out of range (%d, allowed 0..15).\n",
					m_pExtrapolator->m_iSOffset[0]);
				abort();
			}
			if (m_pExtrapolator->m_iSOffset[0] >= m_pExtrapolator->m_iSRange[0]) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.SOffsetX (%d) >= SRangeX (%d.\n",
					m_pExtrapolator->m_iSOffset[0],m_pExtrapolator->m_iSRange[0]);
				abort();
			}
			bs->WriteBits(m_pExtrapolator->m_iSOffset[0],4);

			// OffsetY
			if ((m_pExtrapolator->m_iSOffset[1] < 0) || (m_pExtrapolator->m_iSOffset[1] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.SOffsetY out of range (%d, allowed 0..15).\n",
					m_pExtrapolator->m_iSOffset[1]);
				abort();
			}
			if (m_pExtrapolator->m_iSOffset[1] >= m_pExtrapolator->m_iSRange[1]) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.SOffsetY (%d) >= SRangeY (%d.\n",
					m_pExtrapolator->m_iSOffset[1],m_pExtrapolator->m_iSRange[1]);
				abort();
			}
			bs->WriteBits(m_pExtrapolator->m_iSOffset[1],4);

			// OffsetZ
			if ((m_pExtrapolator->m_iSOffset[2] < 0) || (m_pExtrapolator->m_iSOffset[2] > 15)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.SOffsetZ out of range (%d, allowed 0..15).\n",
					m_pExtrapolator->m_iSOffset[2]);
				abort();
			}
			if (m_pExtrapolator->m_iSOffset[2] >= m_pExtrapolator->m_iSRange[2]) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.SOffsetZ (%d) >= SRangeZ (%d.\n",
					m_pExtrapolator->m_iSOffset[2],m_pExtrapolator->m_iSRange[2]);
				abort();
			}
			bs->WriteBits(m_pExtrapolator->m_iSOffset[2],4);

			// CrossS
			if (m_pExtrapolator->m_bCrossS)
				bs->WriteBit(1);
			else
				bs->WriteBit(0);

			// CrossRangeS
			if (m_pExtrapolator->m_bCrossRangeS)
				bs->WriteBit(1);
			else
				bs->WriteBit(0);

			// Wrap
			if (m_pExtrapolator->m_bWrap)
				bs->WriteBit(1);
			else
				bs->WriteBit(0);

			// DistExpo
			if ((m_pExtrapolator->m_fDistExpo < 0) || (m_pExtrapolator->m_fDistExpo > 16.0)) {
				m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.DistExpo out of range (%.3f, allowed 0..16).\n",
					m_pExtrapolator->m_fDistExpo);
				abort();
			}
			bs->WriteBits((unsigned int)floor(m_pExtrapolator->m_fDistExpo*1000.0+0.5),14);

			// CorrFactor
			bs->WriteBits(0,11);
		}

		// TRange
		if ((m_pExtrapolator->m_iTRange < 0) || (m_pExtrapolator->m_iTRange > 15)) {
			m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.TRange out of range (%d, allowed 0..15).\n",
				m_pExtrapolator->m_iTRange);
			abort();
		}
		bs->WriteBits(m_pExtrapolator->m_iTRange,4);

		// TOrder
		if ((m_pExtrapolator->m_iTOrder < 0) || (m_pExtrapolator->m_iTOrder > 15)) {
			m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.TOrder out of range (%d, allowed 0..15).\n",
				m_pExtrapolator->m_iTOrder);
			abort();
		}
		bs->WriteBits(m_pExtrapolator->m_iTOrder,4);

		// CrossT
		if (m_pExtrapolator->m_bCrossT)
			bs->WriteBit(1);
		else
			bs->WriteBit(0);

		// CrossRangeT
		if (m_pExtrapolator->m_bCrossRangeT)
			bs->WriteBit(1);
		else
			bs->WriteBit(0);

		// TimeExpo
		if ((m_pExtrapolator->m_fTimeExpo < 0) || (m_pExtrapolator->m_fTimeExpo > 16.0)) {
			m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolator.TimeExpo out of range (%.3f, allowed 0..16).\n",
				m_pExtrapolator->m_fTimeExpo);
			abort();
		}
		bs->WriteBits((unsigned int)floor(m_pExtrapolator->m_fTimeExpo*1000.0+0.5),14);

		// Reserved
		bs->WriteBit(0);

	} else { // Position

		// Use Position Extrapolator?
		if (m_pExtrapolatorXYZ == NULL) {
			bs->WriteBit(0);
			return;
		}
		bs->WriteBit(1);

		// TRange
		if ((m_pExtrapolatorXYZ->m_iTRange < 0) || (m_pExtrapolatorXYZ->m_iTRange > 15)) {
			m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorXYZ.TRange out of range (%d, allowed 0..15).\n",
				m_pExtrapolatorXYZ->m_iTRange);
			abort();
		}
		bs->WriteBits(m_pExtrapolatorXYZ->m_iTRange,4);

		// TOrder
		if ((m_pExtrapolatorXYZ->m_iTOrder < 0) || (m_pExtrapolatorXYZ->m_iTOrder > 15)) {
			m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorXYZ.TOrder out of range (%d, allowed 0..15).\n",
				m_pExtrapolatorXYZ->m_iTOrder);
			abort();
		}
		bs->WriteBits(m_pExtrapolatorXYZ->m_iTOrder,4);

		// TimeExpo
		if ((m_pExtrapolatorXYZ->m_fTimeExpo < 0) || (m_pExtrapolatorXYZ->m_fTimeExpo > 16.0)) {
			m_IF.eprintf("CBQBEngine::ExportExtrapolatorSettings(): Error: m_pExtrapolatorXYZ.TimeExpo out of range (%.3f, allowed 0..16).\n",
				m_pExtrapolatorXYZ->m_fTimeExpo);
			abort();
		}
		bs->WriteBits((unsigned int)floor(m_pExtrapolatorXYZ->m_fTimeExpo*1000.0+0.5),14);

		// Reserved
		bs->WriteBit(0);
	}
}


void CBQBEngine::ImportExtrapolatorSettings(CBQBBitSet *bs, bool volumetric) {

	int srangex, srangey, srangez, trange, sorder, torder, offsetx, offsety, offsetz, postrange, postorder, resx, resy, resz;
	bool crosss, crosst, wrap, crossranges, crossranget;
	double distexpo, timeexpo, postimeexpo;
	bool predcorr;


	if (volumetric) {

		if (!bs->ReadBit()) {
			if (m_pExtrapolator != NULL) {
				m_IF.eprintf("CBQBEngine::ImportExtrapolatorSettings(): Warning: No volumetric extrapolator specified in input frame, but m_pExtrapolator != NULL.\n");
				delete m_pExtrapolator;
				m_pExtrapolator = NULL;
			}
			if (m_pExtrapolatorCorr != NULL) {
				m_IF.eprintf("CBQBEngine::ImportExtrapolatorSettings(): Warning: No volumetric extrapolator specified in input frame, but m_pExtrapolatorCorr != NULL.\n");
				delete m_pExtrapolatorCorr;
				m_pExtrapolatorCorr = NULL;
			}
			return;
		}

		// ResX
		resx = bs->ReadBitsInteger(10);

		// ResY
		resy = bs->ReadBitsInteger(10);

		// ResZ
		resz = bs->ReadBitsInteger(10);

		// PredCorr
		predcorr = bs->ReadBit();

		// SRangeX
		srangex = bs->ReadBitsInteger(4);

		// SRangeY
		srangey = bs->ReadBitsInteger(4);

		// SRangeZ
		srangez = bs->ReadBitsInteger(4);

		// SOrder
		sorder = ((int)bs->ReadBitsInteger(4));

		// OffsetX
		offsetx = bs->ReadBitsInteger(4);

		// OffsetY
		offsety = bs->ReadBitsInteger(4);

		// OffsetZ
		offsetz = bs->ReadBitsInteger(4);

		// CrossS
		crosss = bs->ReadBit();

		// CrossRangeS
		crossranges = bs->ReadBit();

		// Wrap
		wrap = bs->ReadBit();

		// DistExpo
		distexpo = bs->ReadBitsInteger(14) / 1000.0;

		// CorrFactor
		bs->ReadBitsInteger(11);

		// TRange
		trange = bs->ReadBitsInteger(4);

		// TOrder
		torder = ((int)bs->ReadBitsInteger(4));

		// CrossT
		crosst = bs->ReadBit();

		// CrossRangeT
		crossranget = bs->ReadBit();

		// TimeExpo
		timeexpo = bs->ReadBitsInteger(14) / 1000.0;

		// Reserved
		if (bs->ReadBit()) {
			m_IF.eprintf("CBQBEngine::ImportExtrapolatorSettings(): Error: Reserved bit after volumetric extrapolator is 1.\n");
			m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
			m_IF.printf("\n");
			abort();
		}


		if (m_pExtrapolator != NULL) {

			if (predcorr) {

				if (m_pExtrapolatorCorr == NULL)
					goto _volcheckfail;

				// SRangeX
				if (m_pExtrapolatorCorr->m_iSRange[0] != srangex)
					goto _volcheckfail;

				// SRangeY
				if (m_pExtrapolatorCorr->m_iSRange[1] != srangey)
					goto _volcheckfail;

				// SRangeZ
				if (m_pExtrapolatorCorr->m_iSRange[2] != srangez)
					goto _volcheckfail;

				// SOrder
				if (m_pExtrapolatorCorr->m_iSOrder != sorder)
					goto _volcheckfail;

				// OffsetX
				if (m_pExtrapolatorCorr->m_iSOffset[0] != offsetx)
					goto _volcheckfail;

				// OffsetY
				if (m_pExtrapolatorCorr->m_iSOffset[1] != offsety)
					goto _volcheckfail;

				// OffsetZ
				if (m_pExtrapolatorCorr->m_iSOffset[2] != offsetz)
					goto _volcheckfail;

				// CrossS
				if (m_pExtrapolatorCorr->m_bCrossS != crosss)
					goto _volcheckfail;

				// CrossRangeS
				if (m_pExtrapolatorCorr->m_bCrossRangeS != crossranges)
					goto _volcheckfail;

				// Wrap
				if (m_pExtrapolatorCorr->m_bWrap != wrap)
					goto _volcheckfail;

				// DistExpo
				if (fabs(m_pExtrapolatorCorr->m_fDistExpo-distexpo) > 1.0e-6)
					goto _volcheckfail;

			} else {

				if (m_pExtrapolatorCorr != NULL)
					goto _volcheckfail;

				// SRangeX
				if (m_pExtrapolator->m_iSRange[0] != srangex)
					goto _volcheckfail;

				// SRangeY
				if (m_pExtrapolator->m_iSRange[1] != srangey)
					goto _volcheckfail;

				// SRangeZ
				if (m_pExtrapolator->m_iSRange[2] != srangez)
					goto _volcheckfail;

				// SOrder
				if (m_pExtrapolator->m_iSOrder != sorder)
					goto _volcheckfail;

				// OffsetX
				if (m_pExtrapolator->m_iSOffset[0] != offsetx)
					goto _volcheckfail;

				// OffsetY
				if (m_pExtrapolator->m_iSOffset[1] != offsety)
					goto _volcheckfail;

				// OffsetZ
				if (m_pExtrapolator->m_iSOffset[2] != offsetz)
					goto _volcheckfail;

				// CrossS
				if (m_pExtrapolator->m_bCrossS != crosss)
					goto _volcheckfail;

				// CrossRangeS
				if (m_pExtrapolator->m_bCrossRangeS != crossranges)
					goto _volcheckfail;

				// Wrap
				if (m_pExtrapolator->m_bWrap != wrap)
					goto _volcheckfail;

				// DistExpo
				if (fabs(m_pExtrapolator->m_fDistExpo-distexpo) > 1.0e-6)
					goto _volcheckfail;
			}

			// ResX
			if (m_pExtrapolator->m_iRes[0] != resx)
				goto _volcheckfail;

			// ResY
			if (m_pExtrapolator->m_iRes[1] != resy)
				goto _volcheckfail;

			// ResZ
			if (m_pExtrapolator->m_iRes[2] != resz)
				goto _volcheckfail;

			// TRange
			if (m_pExtrapolator->m_iTRange != trange)
				goto _volcheckfail;

			// TOrder
			if (m_pExtrapolator->m_iTOrder != torder)
				goto _volcheckfail;

			// CrossT
			if (m_pExtrapolator->m_bCrossT != crosst)
				goto _volcheckfail;

			// CrossRangeT
			if (m_pExtrapolator->m_bCrossRangeT != crossranget)
				goto _volcheckfail;

			// TimeExpo
			if (fabs(m_pExtrapolator->m_fTimeExpo-timeexpo) > 1.0e-6)
				goto _volcheckfail;

			return;

_volcheckfail:
			m_IF.eprintf("CBQBEngine::ImportExtrapolatorSettings(): Warning: Reloading volumetric extrapolator with different settings.\n");
			delete m_pExtrapolator;
			m_pExtrapolator = NULL;
		}

		if (m_pExtrapolatorCorr != NULL) {
			delete m_pExtrapolatorCorr;
			m_pExtrapolatorCorr = NULL;
		}

		if (predcorr) {

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("      Initializing Volumetric Predictor Extrapolator (%d|%d)...\n",trange,torder);
			m_pExtrapolator = new CBQBExtrapolator(m_IF);
			m_pExtrapolator->Initialize(
				resx,     // resx,
				resy,     // resy,
				resz,     // resz,
				1,        // srangex,
				1,        // srangey,
				1,        // srangez,
				trange,   // trange,
				0,        // sorder,
				torder,   // torder,
				0,        // offsetx,
				0,        // offsety,
				0,        // offsetz,
				false,    // crosss,
				false,    // crosst,
				false,    // wrap,
				false,    // crossranges,
				false,    // crossranget
				0,        // distexpo
				timeexpo, // timeexpo
				false
			);

			if (m_IF.IsPL(BQB_PL_STANDARD)) {
				if ((srangex == srangey) && (srangex == srangez))
					m_IF.printf("      Initializing Volumetric Corrector Extrapolator (%d/%d)...\n",srangex,sorder);
				else
					m_IF.printf("      Initializing Volumetric Corrector Extrapolator (%d|%d|%d/%d)...\n",
						srangex,srangey,srangez,sorder);
			}

			m_pExtrapolatorCorr = new CBQBExtrapolator(m_IF);
			m_pExtrapolatorCorr->Initialize(
				resx,        // resx,
				resy,        // resy,
				resz,        // resz,
				srangex,     // srangex,
				srangey,     // srangey,
				srangez,     // srangez,
				1,           // trange,
				sorder,      // sorder,
				0,           // torder,
				offsetx,     // offsetx,
				offsety,     // offsety,
				offsetz,     // offsetz,
				crosss,      // crosss,
				false,       // crosst,
				wrap,        // wrap,
				crossranges, // crossranges,
				false,       // crossranget
				distexpo,    // distexpo
				0,           // timeexpo
				false
			);

		} else {

			if (m_IF.IsPL(BQB_PL_STANDARD)) {
				if ((srangex == srangey) && (srangex == srangez))
					m_IF.printf("      Initializing Volumetric Extrapolator (%d/%d;%d/%d)...\n",
						srangex,sorder,trange,torder);
				else
					m_IF.printf("      Initializing Volumetric Extrapolator (%d|%d|%d/%d;%d/%d)...\n",
						srangex,srangey,srangez,sorder,trange,torder);
			}

			m_pExtrapolator = new CBQBExtrapolator(m_IF);
			m_pExtrapolator->Initialize(
				resx,        // resx,
				resy,        // resy,
				resz,        // resz,
				srangex,     // srangex,
				srangey,     // srangey,
				srangez,     // srangez,
				trange,      // trange,
				sorder,      // sorder,
				torder,      // torder,
				offsetx,     // offsetx,
				offsety,     // offsety,
				offsetz,     // offsetz,
				crosss,      // crosss,
				crosst,      // crosst,
				wrap,        // wrap,
				crossranges, // crossranges,
				crossranget, // crossranget
				distexpo,    // distexpo
				timeexpo,    // timeexpo
				false
			);
		}
	
	} else {

		// Use Position Extrapolator?
		if (!bs->ReadBit()) {
			if (m_pExtrapolatorXYZ != NULL) {
				m_IF.eprintf("CBQBEngine::ImportExtrapolatorSettings(): Warning: No position extrapolator specified in input frame, but m_pExtrapolatorXYZ != NULL.\n");
				delete m_pExtrapolatorXYZ;
				m_pExtrapolatorXYZ = NULL;
			}
			return;
		}

		// TRange
		postrange = bs->ReadBitsInteger(4);

		// TOrder
		postorder = ((int)bs->ReadBitsInteger(4));

		// TimeExpo
		postimeexpo = bs->ReadBitsInteger(14) / 1000.0;

		// Reserved
		if (bs->ReadBit()) {
			m_IF.eprintf("CBQBEngine::ImportExtrapolatorSettings(): Error: Reserved bit after position extrapolator is 1.\n");
			m_IF.printf("Either the BQB file is damaged, or was written with a more recent software version.\n");
			m_IF.printf("\n");
			abort();
		}

		if (m_pExtrapolatorXYZ != NULL) {

			// TRange
			if (m_pExtrapolatorXYZ->m_iTRange != postrange)
				goto _poscheckfail;

			// TOrder
			if (m_pExtrapolatorXYZ->m_iTOrder != postorder)
				goto _poscheckfail;

			// TimeExpo
			if (fabs(m_pExtrapolatorXYZ->m_fTimeExpo-postimeexpo) > 1.0e-6)
				goto _poscheckfail;

			return;

_poscheckfail:
			m_IF.eprintf("CBQBEngine::ImportExtrapolatorSettings(): Warning: Reloading position extrapolator with different settings.\n");
			delete m_pExtrapolatorXYZ;
			m_pExtrapolatorXYZ = NULL;
		}

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("      Initializing Extrapolator (%d|%d)...\n",postrange,postorder);

		m_pExtrapolatorXYZ = new CBQBExtrapolator(m_IF);
		m_pExtrapolatorXYZ->InitializeXYZ(postrange,postorder,postimeexpo,false);
	}
}


// Renders an image of a 3D Hilbert Curve in POVRay with Ambient Occlusion through Radiosity
void CBQBEngine::ExportPOVRayScene(int degree, const char *s) {

	FILE *a;
	unsigned long z, res, crd[3], oldcrd[3];
	double dres;
	CBQBTools bqbtools(m_IF);


	crd[0] = 0;
	crd[1] = 0;
	crd[2] = 0;

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Exporting POVRay scene of 3D Hilbert curve to \"%s\"...\n",s);

	a = bqbtools.BQBOpenFileWrite(s,true);

	res = pow2i(degree);
	dres = (double)res;

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.printf("    Resolution is %lu x %lu x %lu.\n",res,res,res);
		m_IF.printf("    Render with POVRay using \"+W1280 +H960 +J +A0.1 +AM2\" for high quality.\n");
	}

	fprintf(a,"global_settings {\n");
	fprintf(a,"    assumed_gamma 2.2\n");
	fprintf(a,"    adc_bailout 1/255\n");
	fprintf(a,"    charset utf8\n");
	fprintf(a,"    // 1 unit in coordinate system = 1 cm\n");
	fprintf(a,"    // To calculate SSLT for different materials\n");
	fprintf(a,"    mm_per_unit 10\n");
	fprintf(a,"    // Only with +q9 or higher\n");
	fprintf(a,"    subsurface {\n");
	fprintf(a,"        // Number of sample rays for diffuse (first value) and single\n");
	fprintf(a,"        // (second value) scattering\n");
	fprintf(a,"        // Lower values -> better performance\n");
	fprintf(a,"        samples 50, 50\n");
	fprintf(a,"        // Use subsurface scattering on radiosity photons\n");
	fprintf(a,"        radiosity on\n");
	fprintf(a,"    }\n");
	fprintf(a,"\n");
	fprintf(a,"    radiosity {\n");
	fprintf(a,"        gray_threshold 0.0\n");
	fprintf(a,"        brightness 1\n");
	fprintf(a,"        // not smaller than 0.0039\n");
	fprintf(a,"        adc_bailout 0.005\n");
	fprintf(a,"        // not lower than 2 for isosurfaces, not bigger than 3\n");
	fprintf(a,"        recursion_limit 3\n");
	fprintf(a,"        // lowers error bound during pretrace, removes blotches (default 0.5)\n");
	fprintf(a,"        low_error_factor 0.5\n");
	fprintf(a,"        // Fraction of screen to follow bounced rays (default 0.015)\n");
	fprintf(a,"        // Too small ? rendering gets slow; too high ? no natural darkening\n");
	fprintf(a,"        // of crevices as rays get reused\n");
	fprintf(a,"        minimum_reuse 0.015\n");
	fprintf(a,"        maximum_reuse 0.2\n");
	fprintf(a,"        // Min and maxsize of block in mosaic preview during pretrace\n");
	fprintf(a,"        // default 0.08\n");
	fprintf(a,"        pretrace_start 0.08\n");
	fprintf(a,"        // not smaller than 0.02 (default 0.04)\n");
	fprintf(a,"        pretrace_end 0.003 // Use 0.04 for testing\n");
	fprintf(a,"        count 1200 // Use 100 or even 20 for testing\n");
	fprintf(a,"        // Max number of radiosity values to be blended together (default\n");
	fprintf(a,"        // 4, max 20)\n");
	fprintf(a,"        // Smaller than 4 ? patchiness\n");
	fprintf(a,"        nearest_count 20\n");
	fprintf(a,"        // lower error bound ? more artifacts (requires higher count), but\n");
	fprintf(a,"        // also higher quality\n");
	fprintf(a,"        error_bound 0.2 // Use 0.4 for testing\n");
	fprintf(a,"    }\n");
	fprintf(a,"}\n");
	fprintf(a,"\n");
	fprintf(a,"#declare rad = 0.23;\n");
	fprintf(a,"\n");
	fprintf(a,"#declare rot_angle = 0.5;\n");
	fprintf(a,"#declare camera_dist = 800.0;\n");
	fprintf(a,"#declare camera_x = 10.0 * camera_dist / 1566.0 * sin(rot_angle);\n");
	fprintf(a,"#declare camera_y = 4.0 * camera_dist / 1566.0;\n");
	fprintf(a,"#declare camera_z = 10.0 * camera_dist /1566.0 * cos(rot_angle);\n");
	fprintf(a,"\n");
	fprintf(a,"#macro Hilbert_Sphere( pos )\n");
	fprintf(a,"  sphere {\n");
	fprintf(a,"    pos, rad / 16.0\n");
	fprintf(a,"    pigment { rgb < 0.8, 0.8, 0.8 > }\n");
	fprintf(a,"    finish {\n");
	fprintf(a,"      reflection 0\n");
	fprintf(a,"      specular 0.0\n");
	fprintf(a,"      ambient 0.0\n");
	fprintf(a,"      diffuse 0.6\n");
	fprintf(a,"      roughness 0.12\n");
	fprintf(a,"    }\n");
	fprintf(a,"    no_shadow\n");
	fprintf(a,"  }\n");
	fprintf(a,"#end\n");
	fprintf(a,"\n");
	fprintf(a,"#macro Hilbert_Cylinder( pos1, pos2 )\n");
	fprintf(a,"  cylinder {\n");
	fprintf(a,"    pos1, pos2, rad / 16.0\n");
	fprintf(a,"    open\n");
	fprintf(a,"    pigment { rgb < 0.8, 0.8, 0.8 > }\n");
	fprintf(a,"    finish {\n");
	fprintf(a,"      reflection 0\n");
	fprintf(a,"      specular 0.0\n");
	fprintf(a,"      ambient 0.0\n");
	fprintf(a,"      diffuse 0.6\n");
	fprintf(a,"      roughness 0.12\n");
	fprintf(a,"    }\n");
	fprintf(a,"    no_shadow\n");
	fprintf(a,"  }\n");
	fprintf(a,"#end\n");
	fprintf(a,"\n");
	fprintf(a,"camera {\n");
	fprintf(a,"	location < camera_x, camera_y, camera_z >\n");
	fprintf(a,"	sky y\n");
	fprintf(a,"	right -0.3*x*image_width/image_height\n");
	fprintf(a,"	up 0.3*y\n");
	fprintf(a,"	look_at < 0.5, 0.42, 0.5 >\n");
	fprintf(a,"}\n");
	fprintf(a,"\n");
	fprintf(a,"// Solid colored Background\n");
	fprintf(a,"background { rgb < 1.0, 1.0, 1.0 > }\n");
	fprintf(a,"\n");
	fprintf(a,"/***** Color Gradient Background ****/  \n");
	fprintf(a,"/*\n");
	fprintf(a,"sky_sphere {\n");
	fprintf(a,"  pigment {\n");
	fprintf(a,"    gradient y\n");
	fprintf(a,"    color_map {\n");
	fprintf(a,"      [ 0 color rgb < 0.9, 0.9, 0.9 > ]\n");
	fprintf(a,"      [ 1 color rgb < 0.85, 0.75, 1.0 > ]\n");
	fprintf(a,"    }\n");
	fprintf(a,"    scale 1\n");
	fprintf(a,"    translate -0.05\n");
	fprintf(a,"  }\n");
	fprintf(a,"}\n");
	fprintf(a,"\n");
	fprintf(a,"// Two light sources\n");
	fprintf(a,"light_source { < 21.54*sin(rot_angle-0.7297), 20, 21.54*cos(rot_angle-0.7297) > color rgb 0.8 }\n");
	fprintf(a,"light_source { < 32.02*sin(rot_angle+0.547), 12, 32.02*cos(rot_angle+0.547) > color rgb 0.5 }      */\n");
	fprintf(a,"\n");
	fprintf(a,"// Weak highlight within the sky sphere\n");
	fprintf(a,"light_source {\n");
	fprintf(a,"    500\n");
	fprintf(a,"    // Light intensity depends on color brightness\n");
	fprintf(a,"    rgb .5\n");
	fprintf(a,"    // Cylindrical and spot lights can be area lights as well\n");
	fprintf(a,"    area_light 100*x, 100*y, 100, 100\n");
	fprintf(a,"    // Automatically adapt number of light sources\n");
	fprintf(a,"    adaptive 1\n");
	fprintf(a,"    // Add random positional noise to light source positions\n");
	fprintf(a,"    jitter\n");
	fprintf(a,"    // Use circular area light\n");
	fprintf(a,"    circular\n");
	fprintf(a,"    // Always normalize the light area with regard to the ray\n");
	fprintf(a,"    orient\n");
	fprintf(a,"}\n");
	fprintf(a,"// invisible pseudo-sky sphere\n");
	fprintf(a,"sphere {\n");
	fprintf(a,"    0, 1\n");
	fprintf(a,"    pigment {\n");
	fprintf(a,"        color rgb 1\n");
	fprintf(a,"    }\n");
	fprintf(a,"    finish {\n");
	fprintf(a,"        emission .5\n");
	fprintf(a,"    }\n");
	fprintf(a,"    // make sure the sphere is big enough to cover the entire scene\n");
	fprintf(a,"    scale 1000\n");
	fprintf(a,"    hollow\n");
	fprintf(a,"    // make sphere invisible for direct rays\n");
	fprintf(a,"    no_image\n");
	fprintf(a,"}\n");

	fprintf(a,"// Here the actual Hilbert curve data follows\n\n");

	for (z=0;z<res*res*res;z++) {
		oldcrd[0] = crd[0];
		oldcrd[1] = crd[1];
		oldcrd[2] = crd[2];
		hilbert_i2c(3,degree,z,crd);
		fprintf(a,"Hilbert_Sphere( < %f, %f, %f > )\n",crd[0]/dres,crd[1]/dres,crd[2]/dres);
		if (z > 0) {
			fprintf(a,"Hilbert_Cylinder( < %f, %f, %f >, < %f, %f, %f > )\n",
				oldcrd[0]/dres,oldcrd[1]/dres,oldcrd[2]/dres,crd[0]/dres,crd[1]/dres,crd[2]/dres);
		}
	}
	fprintf(a,"\n");


	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Done.\n");

	fclose(a);
}


void CBQBEngine::ExportXYZScene(int degree, const char *s) {

	unsigned long z, res, crd[3];
	FILE *a;
	CBQBTools bqbtools(m_IF);


	res = pow2i(degree);

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Exporting XYZ scene of 3D Hilbert curve to \"%s\"...\n",s);

	a = bqbtools.BQBOpenFileWrite(s,true);

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("    Resolution is %lu x %lu x %lu.\n",res,res,res);

	fprintf(a,"%lu\nHilbert Scene\n",res*res*res);
	for (z=0;z<res*res*res;z++)
		fprintf(a,"C  %.1f  0.0  0.0\n",(double)z*1.4);

	fprintf(a,"%lu\nHilbert Scene\n",res*res*res);
	for (z=0;z<res*res*res;z++) {
		hilbert_i2c(3,degree,z,crd);
		fprintf(a,"C  %lu.0  %lu.0  %lu.0\n",crd[0],crd[1],crd[2]);
	}

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Done.\n");
}



