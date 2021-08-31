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

#include "bqb_driver.h"
#include "bqb_format.h"
#include "bqb_integerengine.h"
#include "bqb_bitset.h"
#include <algorithm>
#include <time.h>


const char *GetRevisionInfo_bqb_driver(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_driver() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



CBQBDriver::CBQBDriver(CBQBInterface &i) :
		m_ParmPosAndVol(i),
		m_ParmPos(i),
		m_iHistogramCounter(0),
		m_iOptIncludeFirst(-1),
		m_IF(i)
	{
	m_pEngine = new CBQBEngine(m_IF);
}


CBQBDriver::~CBQBDriver() {

	if (m_pEngine != NULL) {
		delete m_pEngine;
		m_pEngine = NULL;
	}
}


bool CBQBDriver::CompressCube(
		const char *inp,
		const char *outp,
		const char *ref,
		int start,
		int steps,
		int stride,
		int keyfreq,
		CBQBParameterSet_PosAndVol *parm,
		int optimize,
		int optsteps,
		bool onlyopt,
		bool comment,
		bool compare,
		bool dummyread
		) {

	const CBQBCubeFrame *cfr;
	CBQBCubeFrame *cfr2;
	int z, z2, i, k, co, ao, ft, lcsize, lasize, morder, histused, thu, vmo, pmo;
	double tb, tb2, tf, tf2, tot;
	long insize;
	//std::vector<int> ia, ia2;
	CBQBBitSet bsat(m_IF), bscu(m_IF), bshe(m_IF), bstemp(m_IF);
	bool err, to;
	unsigned long t0, rxyz=1;
	CBQBFile bf(m_IF);
	CBQBStatistics stat(m_IF);


	if (parm == NULL) {
		m_IF.eprintf("CBQBDriver::CompressCube(): Error: parm == NULL.\n");
		return false;
	}

	if (parm->GetVolUseExtra())
		if (parm->GetVolExtraTRange() > parm->GetVolOrder())
			parm->SetVolOrder(parm->GetVolExtraTRange());

	if (parm->GetPosUseExtra())
		if (parm->GetPosExtraTRange() > parm->GetPosOrder())
			parm->SetPosOrder(parm->GetPosExtraTRange());

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Opening cube file \"%s\" ...\n",inp);

	m_pEngine->m_pReadCacheCube = new CBQBReadCubeCache(m_IF);
	m_pEngine->m_pReadCacheCube->SetReadParameters(parm->GetVolEps(),parm->GetVolSigni(),parm->GetPosPrecision());
	if (!m_pEngine->m_pReadCacheCube->FileOpenRead(inp,ref)) {
		m_IF.eprintf("Error: Could not open file for reading.\n");
		return false;
	}

	if (optimize != 0) {

		if (MAX(parm->GetPosOrder(),parm->GetVolOrder())+1 > optsteps)
			m_pEngine->m_pReadCacheCube->SetHistoryDepth(MAX(parm->GetPosOrder(),parm->GetVolOrder())+1);
		else
			m_pEngine->m_pReadCacheCube->SetHistoryDepth(optsteps);

		m_IF.printf("\n");
		m_pEngine->m_pReadCacheCube->CacheSteps(optsteps,true);
		m_IF.printf("\n");

		if (m_pEngine->m_pReadCacheCube->GetCachedSteps() == 0) {
			m_IF.eprintf("Error: Zero frames read. Aborting.\n");
			return false;
		}

		if (m_iOptIncludeFirst == -1) {
			if ((m_pEngine->m_pReadCacheCube->GetCachedSteps() < optsteps) || (optsteps < 10)) {
				if (m_IF.IsPL(BQB_PL_STANDARD)) {
					m_IF.printf("        Low step count (%d), switching on -optstart.\n",m_pEngine->m_pReadCacheCube->GetCachedSteps());
					m_IF.printf("\n");
				}
				m_iOptIncludeFirst = 1;
			} else {
				if (m_IF.IsPL(BQB_PL_STANDARD)) {
					m_IF.printf("        Large step count (%d), switching off -optstart.\n",m_pEngine->m_pReadCacheCube->GetCachedSteps());
					m_IF.printf("\n");
				}
				m_iOptIncludeFirst = 0;
			}
		}

		m_pEngine->m_pReadCacheCube->SetStartLock();
	
		if (!OptimizeXYZParameters(
				optimize, // int     level,
				parm->GetPositionParameterSet(),
				comment, // bool    comment, 
				true // bool    cubepart
		)) {
			m_IF.eprintf("Error: Position data parameter optimization failed.\n");
			delete m_pEngine->m_pReadCacheCube;
			m_pEngine->m_pReadCacheCube = NULL;
			return false;
		}

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("\n");
		
		if (!OptimizeCubeParameters(
				optimize, // int     level,
				parm->GetVolumetricParameterSet()
		)) {
			m_IF.eprintf("Error: Volumetric data parameter optimization failed.\n");
			delete m_pEngine->m_pReadCacheCube;
			m_pEngine->m_pReadCacheCube = NULL;
			return false;
		}

		if (m_IF.IsPL(BQB_PL_STANDARD)) {
			
			m_IF.bprintf("\n");
			m_IF.bprintf("    The optimal parameters for this trajectory are:\n");
			m_IF.bprintf("\n");

			m_IF.printf("%s",parm->ToString(4).c_str());
			m_IF.printf("\n");

			m_IF.bprintf("    Parameter key:");
			m_IF.printf("  %s\n",parm->ToKey().c_str());
			m_IF.printf("\n");
		}

		if (onlyopt) {
			m_pEngine->m_pReadCacheCube->CloseFile();
			delete m_pEngine->m_pReadCacheCube;
			m_pEngine->m_pReadCacheCube = NULL;
			return true;
		}

		m_pEngine->m_pReadCacheCube->RewindReadPos();

		m_pEngine->m_pReadCacheCube->LiftStartLock();

		if (parm->GetVolUseExtra())
			if (parm->GetVolExtraTRange() > parm->GetVolOrder())
				parm->SetVolOrder(parm->GetVolExtraTRange());

		if (parm->GetPosUseExtra())
			if (parm->GetPosExtraTRange() > parm->GetPosOrder())
				parm->SetPosOrder(parm->GetPosExtraTRange());

	} else {

		m_pEngine->m_pReadCacheCube->SetHistoryDepth(MAX(parm->GetPosOrder(),parm->GetVolOrder())+1);
	}

	t0 = (unsigned long)time(NULL);

	stat.m_oStat.reset();

	if (!dummyread) {
		if (outp != NULL) {
			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("Opening compressed file \"%s\" ...\n",outp);
			if (!bf.OpenWriteAppend(outp)) {
				m_IF.eprintf("Error: Could not open file for writing.\n");
				m_pEngine->m_pReadCacheCube->CloseFile();
				delete m_pEngine->m_pReadCacheCube;
				m_pEngine->m_pReadCacheCube = NULL;
				return false;
			}
		}
	} else {
		ref = NULL;
		outp = NULL;
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("Dummy Read Mode (no compression / output is performed).\n");
	}

	m_pEngine->m_oaOutputCubeBuf.resize(parm->GetVolOrder()+1);
	for (z=0;z<parm->GetVolOrder()+1;z++)
		m_pEngine->m_oaOutputCubeBuf[z] = NULL;
	m_pEngine->m_iOutputCubeBufPos = 0;

	m_pEngine->m_oaOutputAtomBuf.resize(parm->GetPosOrder()+1);
	for (z=0;z<parm->GetPosOrder()+1;z++)
		m_pEngine->m_oaOutputAtomBuf[z] = NULL;
	m_pEngine->m_iOutputAtomBufPos = 0;

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Starting process...\n");

	if (start != 0) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("Fast-forwarding to step %d: ",start+1);
		for (z=0;z<start;z++) {
			if (!m_pEngine->m_pReadCacheCube->SkipOneStep()) {
				m_pEngine->m_pReadCacheCube->CloseFile();
				delete m_pEngine->m_pReadCacheCube;
				m_pEngine->m_pReadCacheCube = NULL;
				return false;
			}
			if (m_IF.IsPL(BQB_PL_STANDARD)) {
				m_IF.printf(".");
				m_IF.FlushLog();
			}
		}
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf(" Done.\n");
	}

	if (parm->GetVolUseExtra())
		vmo = parm->GetVolExtraTRange();
	else
		vmo = parm->GetVolOptOrder();

	if (parm->GetPosUseExtra())
		pmo = parm->GetPosExtraTRange();
	else
		pmo = parm->GetPosOptOrder();

	morder = vmo;
	if (pmo > morder)
		morder = pmo;

	stat.ResetStatistics();

	err = false;
	i = 0;
	k = 0;
	tb = 0;
	tb2 = 0;
	tf = 0;
	tf2 = 0;
	to = false;
	if (parm->GetVolOptOrder() && !parm->GetVolUseExtra())
		lcsize = 1000000000;
	else
		lcsize = -1;
	if (parm->GetPosOptOrder() && !parm->GetPosUseExtra())
		lasize = 1000000000;
	else
		lasize = -1;


	while (true) {

		if (m_IF.IsPL(BQB_PL_VERBOSE)) {
			m_IF.printf("\n");
			m_IF.printf("\n");
		}

		if (m_IF.IsPL(BQB_PL_STANDARD)) {
			if (i != 0)
				m_IF.printf("    Step %6d ... (avg %5.1f s/step)\n",i+1,((double)(time(NULL)-t0))/i);
			else
				m_IF.printf("    Step %6d ...\n",i+1);
		}

		cfr = m_pEngine->m_pReadCacheCube->GetNextFrame();

		if (cfr == NULL) {
			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("      Reached end of input trajectory.\n");
			break;
		}

		if (i == 0) {
			rxyz = cfr->m_iResXYZ;

			if (parm->GetVolUseExtra()) {

				if (parm->GetVolExtraPredCorr()) {

					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("      Initializing Volumetric Predictor Extrapolator (%d/%d)...\n",
							parm->GetVolExtraTRange(),parm->GetVolExtraTOrder());

					m_pEngine->m_pExtrapolator = new CBQBExtrapolator(m_IF);
					m_pEngine->m_pExtrapolator->Initialize(
						cfr->m_iRes[0], // resx,
						cfr->m_iRes[1], // resy,
						cfr->m_iRes[2], // resz,
						1,              // srangex,
						1,              // srangey,
						1,              // srangez,
						parm->GetVolExtraTRange(),   // trange,
						0,              // sorder,
						parm->GetVolExtraTOrder(),   // torder,
						0,              // offsetx,
						0,              // offsety,
						0,              // offsetz,
						false,          // crosss,
						false,          // crosst,
						false,          // wrap,
						false,          // crossranges,
						false,          // crossranget
						0,              // distexpo
						parm->GetVolExtraTimeExpo(), // timeexpo
						false
					);

					if (m_IF.IsPL(BQB_PL_STANDARD)) {
						if (parm->IsVolExtraSRangeXYZEqual())
							m_IF.printf("      Initializing Volumetric Corrector Extrapolator (%d/%d)...\n",
								parm->GetVolExtraSRangeX(),parm->GetVolExtraSOrder());
						else
							m_IF.printf("      Initializing Volumetric Corrector Extrapolator (%d|%d|%d/%d)...\n",
								parm->GetVolExtraSRangeX(),parm->GetVolExtraSRangeY(),parm->GetVolExtraSRangeZ(),parm->GetVolExtraSOrder());
					}

					m_pEngine->m_pExtrapolatorCorr = new CBQBExtrapolator(m_IF);
					m_pEngine->m_pExtrapolatorCorr->Initialize(
						cfr->m_iRes[0],     // resx,
						cfr->m_iRes[1],     // resy,
						cfr->m_iRes[2],     // resz,
						parm->GetVolExtraSRangeX(),      // srangex,
						parm->GetVolExtraSRangeY(),      // srangey,
						parm->GetVolExtraSRangeZ(),      // srangez,
						1,                  // trange,
						parm->GetVolExtraSOrder(),       // sorder,
						0,                  // torder,
						parm->GetVolExtraOffsetX(),      // offsetx,
						parm->GetVolExtraOffsetY(),      // offsety,
						parm->GetVolExtraOffsetZ(),      // offsetz,
						parm->GetVolExtraCrossS(),       // crosss,
						false,              // crosst,
						parm->GetVolExtraWrap(),         // wrap,
						parm->GetVolExtraCrossRangeS(),  // crossranges,
						false,              // crossranget
						parm->GetVolExtraDistExpo(),     // distexpo
						0,                  // timeexpo
						false
					);

				} else {

					if (m_IF.IsPL(BQB_PL_STANDARD)) {

						if (parm->IsVolExtraSRangeXYZEqual())
							m_IF.printf("      Initializing Volumetric Extrapolator (%d/%d;%d/%d)...\n",
								parm->GetVolExtraSRangeX(),
								parm->GetVolExtraSOrder(),
								parm->GetVolExtraTRange(),
								parm->GetVolExtraTOrder()
							);
						else
							m_IF.printf("      Initializing Volumetric Extrapolator (%d|%d|%d/%d;%d|%d)...\n",
								parm->GetVolExtraSRangeX(),
								parm->GetVolExtraSRangeY(),
								parm->GetVolExtraSRangeZ(),
								parm->GetVolExtraSOrder(),
								parm->GetVolExtraTRange(),
								parm->GetVolExtraTOrder()
							);
					}

					m_pEngine->m_pExtrapolator = new CBQBExtrapolator(m_IF);
					m_pEngine->m_pExtrapolator->Initialize(
						cfr->m_iRes[0],     // resx,
						cfr->m_iRes[1],     // resy,
						cfr->m_iRes[2],     // resz,
						parm->GetVolExtraSRangeX(),      // srangex,
						parm->GetVolExtraSRangeY(),      // srangey,
						parm->GetVolExtraSRangeZ(),      // srangez,
						parm->GetVolExtraTRange(),       // trange,
						parm->GetVolExtraSOrder(),       // sorder,
						parm->GetVolExtraTOrder(),       // torder,
						parm->GetVolExtraOffsetX(),      // offsetx,
						parm->GetVolExtraOffsetY(),      // offsety,
						parm->GetVolExtraOffsetZ(),      // offsetz,
						parm->GetVolExtraCrossS(),       // crosss,
						parm->GetVolExtraCrossT(),       // crosst,
						parm->GetVolExtraWrap(),         // wrap,
						parm->GetVolExtraCrossRangeS(),  // crossranges,
						parm->GetVolExtraCrossRangeT(),  // crossranget
						parm->GetVolExtraDistExpo(),     // distexpo
						parm->GetVolExtraTimeExpo(),     // timeexpo
						false
					);
				}
			}

			if (parm->GetPosUseExtra()) {

				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("      Initializing Position Extrapolator (%d/%d)...\n",
						parm->GetPosExtraTRange(),parm->GetPosExtraTOrder());

				m_pEngine->m_pExtrapolatorXYZ = new CBQBExtrapolator(m_IF);
				m_pEngine->m_pExtrapolatorXYZ->InitializeXYZ(
					parm->GetPosExtraTRange(),
					parm->GetPosExtraTOrder(),
					parm->GetPosExtraTimeExpo(),
					false
				);
			}

		}

		for (z=0;z<stride-1;z++) {
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Skipping frame...\n");
			if (!m_pEngine->m_pReadCacheCube->SkipOneStep()) {
				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("      Reached end of input trajectory.\n");
				break;
			}
		}

		if (dummyread) // Just benchmark the CUBE reading speed
			goto _dummy;

		if (m_pEngine->m_pExtrapolator != NULL)
			m_pEngine->m_pExtrapolator->PushCubeFrame(cfr);

		if (m_pEngine->m_pExtrapolatorXYZ != NULL)
			m_pEngine->m_pExtrapolatorXYZ->PushAtomFrame(cfr->m_pAtoms);

		if ((i != 0) && (keyfreq != 0) && ((i % keyfreq) == 0)) {
			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("      Enforcing key frame.\n");
			k = 0; // Enforce this frame to be a keyframe
		}

_again:
		if (i < vmo)
			co = i;
		else if (k < vmo)
			co = k;
		else
			co = vmo;

		if (i < pmo)
			ao = i;
		else if (k < pmo)
			ao = k;
		else
			ao = pmo;

_aagain:
		bsat.Clear();
		stat.PushStatistics();

		histused = 0;

		m_pEngine->CompressAtomFrame(
			&bsat,
			ao,
			parm->GetPositionParameterSet(),
			histused,
			(i==0), // Store static trajectory information?
			comment,
			&stat
		);

		if (m_IF.IsPL(BQB_PL_STANDARD)) {
			if (m_pEngine->m_pExtrapolatorXYZ == NULL)
				m_IF.printf("      Atoms: Order %d,    output size %9.3f KiB.\n",ao,bsat.GetByteLength()/1024.0);
			else
				m_IF.printf("      Atoms: History %2d, output size %9.3f KiB.\n",histused,bsat.GetByteLength()/1024.0);
		}

		if (lasize >= 0) {
			if (to) {
				to = false;
				if (bsat.GetByteLength() <= bstemp.GetByteLength()) {
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("      Size is indeed smaller with lower order. Limiting atom order to %d.\n",ao);
					stat.PopDiffStatistics();
					parm->SetPosOrder(ao);
					pmo = ao;
					lasize = -1;
				} else {
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("      Size was smaller with higher order. Not limiting.\n");
					stat.PopStatistics();
					bsat = bstemp;
					lasize = bsat.GetByteLength();
				}
			} else {
				if ((ao != 0) && ((double)bsat.GetByteLength() >= (double)lasize*0.97)) {
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("      Size did not decrease any further with order %d. Trying order %d...\n",ao,ao-1);
					bstemp = CBQBBitSet(bsat);
					to = true;
					ao--;
					goto _aagain;
				} else if (ao == pmo)
					lasize = -1;
				else
					lasize = bsat.GetByteLength();
			}
		}

		stat.PopIgnoreStatistics();

_cagain:
		bscu.Clear();
		stat.PushStatistics();

		m_pEngine->CompressCubeFrame(
			&bscu,
			co,
			parm->GetVolumetricParameterSet(),
			thu,
			(i==0), // Store static info?
			&stat
		);

		if (thu > histused)
			histused = thu;

		if (m_IF.IsPL(BQB_PL_STANDARD)) {
			if (m_pEngine->m_pExtrapolator == NULL)
				m_IF.printf("      Cube:  Order %d,    output size %9.3f KiB.\n",co,bscu.GetByteLength()/1024.0);
			else
				m_IF.printf("      Cube:  History %2d, output size %9.3f KiB.\n",thu,bscu.GetByteLength()/1024.0);
		}

		if (lcsize >= 0) {
			if (to) {
				to = false;
				if (bscu.GetByteLength() <= bstemp.GetByteLength()) {
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("      Size is indeed smaller with lower order. Limiting cube order to %d.\n",co);
					stat.PopDiffStatistics();
					parm->SetVolOrder(co);
					vmo = co;
					lcsize = -1;
				} else {
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("      Size was smaller with higher order. Not limiting.\n");
					stat.PopStatistics();
					bscu = CBQBBitSet(bstemp);
					lcsize = bscu.GetByteLength();
				}
			} else {
				if ((co != 0) && ((double)bscu.GetByteLength() >= (double)lcsize*0.97)) {
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("      Size did not decrease any further with order %d. Trying order %d...\n",co,co-1);
					bstemp = CBQBBitSet(bscu);
					to = true;
					co--;
					goto _cagain;
				} else if (co == vmo)
					lcsize = -1;
				else
					lcsize = bscu.GetByteLength();
			}
		}

		stat.PopIgnoreStatistics();

		if (MAX(pmo,vmo) < morder)
			morder = MAX(pmo,vmo);

		if (compare) {

			cfr2 = new CBQBCubeFrame(m_IF);
			m_pEngine->PushOutputCubeFrame(cfr2);
			cfr2->m_pAtoms = new CBQBAtomSet(m_IF);
			m_pEngine->PushOutputAtomFrame(cfr2->m_pAtoms);

			m_pEngine->DecompressAtomFrame(
				&bsat,
				BQB_FRAMETYPE_COMPTRAJ_VERSION,
				parm->GetPositionParameterSet(),
				NULL, // histused
				false,
				true
			);

			m_pEngine->DecompressCubeFrame(
				&bscu,
				BQB_FRAMETYPE_COMPCUBE_VERSION,
				parm->GetVolumetricParameterSet(),
				NULL, // histused
				false,
				true
			);

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Comparing input and output...\n");

			for (z=0;z<(int)cfr2->m_pAtoms->m_oaAtoms.size();z++)
				for (z2=0;z2<3;z2++)
					if (cfr->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2] != cfr2->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2]) {
						m_IF.eprintf("        Error in atom coordinate %d[%d]: %.6f (%ld) != %.6f (%ld)\n",
							z,z2,cfr->m_pAtoms->m_oaAtoms[z]->m_fCoord[z2],cfr->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2],
							cfr2->m_pAtoms->m_oaAtoms[z]->m_fCoord[z2],cfr2->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2]);
						err = true;
						goto _skerr;
					}

			for (z=0;z<cfr->m_iResXYZ;z++)
				if (!ExpMantisEqual(cfr->m_iaExpo[z],cfr->m_iaMantis[z],cfr2->m_iaExpo[z],cfr2->m_iaMantis[z])) {
					m_IF.eprintf("        Error in volumetric data element %7d: %.10G vs %.10G\n",
						z,cfr->m_faBin[z],cfr2->m_faBin[z]);
					err = true;
					goto _skerr;
				}
_skerr:
			if (err) {
				if (k != 0) {
					m_IF.eprintf("Errors occurred. Compressing frame again with history zero.\n");
					err = false;
					k = 0;
					goto _again;
				} else {
					m_IF.eprintf("Errors occurred despite of zero history. Aborting.\n");
					goto _end;
				}
			}

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Done.\n");
		}

		bshe.Clear();

		bshe.WriteBits(bsat.GetByteLength(),32);
		bshe.WriteBits(bscu.GetByteLength(),32);

		tb += bshe.GetByteLength();
		tb += bsat.GetByteLength();
		tb += bscu.GetByteLength();
		tf++;

		stat.m_oStat.m_lOverhead += bshe.GetLength();

		if ((bsat.GetLength()%8) != 0)
			stat.m_oStat.m_lOverhead += 8 - (bsat.GetLength()%8);

		if ((bscu.GetLength()%8) != 0)
			stat.m_oStat.m_lOverhead += 8 - (bscu.GetLength()%8);

		if (i >= morder) {
			tb2 += bshe.GetByteLength();
			tb2 += bsat.GetByteLength();
			tb2 += bscu.GetByteLength();
			tf2++;
		}

		if (outp != NULL) {

			if (i == 0)
				ft = BQB_FRAMETYPE_COMPCUBESTART;
			else if (k == 0)
				ft = BQB_FRAMETYPE_COMPCUBEKEY;
			else
				ft = BQB_FRAMETYPE_COMPCUBE;

			bf.CreateShortFrame(
				ft,
				BQB_FRAMETYPE_COMPCUBE_VERSION,
				i+1
			);

			bf.PushPayload(bshe.m_iaData);
			bf.PushPayload(bsat.m_iaData);
			bf.PushPayload(bscu.m_iaData);
			bf.FinalizeFrame(&stat);
		}

_dummy:
		i++;
		k++;
		if (i == steps)
			break;
	}

	if (outp != NULL)
		bf.WriteIndexFrame(true,&stat);

_end:
	if (outp != NULL)
		bf.Close();

	insize = m_pEngine->m_pReadCacheCube->ftell();

	m_pEngine->m_pReadCacheCube->CloseFile();
	delete m_pEngine->m_pReadCacheCube;
	m_pEngine->m_pReadCacheCube = NULL;

	if (m_pEngine->m_pExtrapolator != NULL) {
		delete m_pEngine->m_pExtrapolator;
		m_pEngine->m_pExtrapolator = NULL;
	}

	if (m_pEngine->m_pExtrapolatorCorr != NULL) {
		delete m_pEngine->m_pExtrapolatorCorr;
		m_pEngine->m_pExtrapolatorCorr = NULL;
	}

	if (m_pEngine->m_pExtrapolatorXYZ != NULL) {
		delete m_pEngine->m_pExtrapolatorXYZ;
		m_pEngine->m_pExtrapolatorXYZ = NULL;
	}

	if (err) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.eprintf("Errors occurred while compressing the cube file.\n");
		return false;
	} else if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.printf("Finished compressing the cube file.\n");
		m_IF.printf("\n");
		m_IF.bprintf("Size contributions:\n");
		tot = double(stat.m_oStat.m_lOverhead+stat.m_oStat.m_lAlphabet+stat.m_oStat.m_lHuffmanTables+stat.m_oStat.m_lTableSwitch+stat.m_oStat.m_lHuffmanData);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) overhead.\n",
			double(stat.m_oStat.m_lOverhead)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lOverhead>>3).string(),
			double(stat.m_oStat.m_lOverhead)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) alphabet data.\n",
			double(stat.m_oStat.m_lAlphabet)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lAlphabet>>3).string(),
			double(stat.m_oStat.m_lAlphabet)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) Huffman tables.\n",
			double(stat.m_oStat.m_lHuffmanTables)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lHuffmanTables>>3).string(),
			double(stat.m_oStat.m_lHuffmanTables)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) table switching.\n",
			double(stat.m_oStat.m_lTableSwitch)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lTableSwitch>>3).string(),
			double(stat.m_oStat.m_lTableSwitch)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) payload data.\n",
			double(stat.m_oStat.m_lHuffmanData)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lHuffmanData>>3).string(),
			double(stat.m_oStat.m_lHuffmanData)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) in total.\n",
			tot/1024.0/1024.0/8.0,
			((stat.m_oStat.m_lOverhead+stat.m_oStat.m_lAlphabet+stat.m_oStat.m_lHuffmanTables+stat.m_oStat.m_lTableSwitch+stat.m_oStat.m_lHuffmanData)>>3).string(),
			100.0
		);
		m_IF.printf("\n");
		if (tf > 0) {
			m_IF.bprintf("\n");
			m_IF.bprintf("  * Totals:\n");
			m_IF.printf("      %9.3f KiB per frame on average.\n",tb/1024.0/tf);
			m_IF.printf("      %9.3f bits per bin entry on average.\n",tb/tf/rxyz*8.0);
			m_IF.printf("      Compression ratio of %.3f : 1\n",(double)insize/tb);
		}
		if (tf2 > 0) {
			m_IF.bprintf("\n");
			m_IF.bprintf("  * Starting from step %d:\n",morder+1);
			m_IF.printf("      %9.3f KiB per frame on average.\n",tb2/1024.0/tf2);
			m_IF.printf("      %9.3f bits per bin entry on average.\n",(tb2/tf2)/rxyz*8.0);
			m_IF.printf("      Compression ratio of %.3f : 1\n",((double)insize/tf*tf2)/tb2);
		}
		m_IF.printf("\n");
	}

	return true;
}


double CBQBDriver::CompressCubeBenchmark(
		CBQBParameterSet_Volumetric *parm,
		bool realsize,
		double *presi,
		std::vector<std::vector<int> > *tciaa
	) {

	const CBQBCubeFrame *cfr;
	int i, co, histused, mo;
	double tb, tf, resi, tresi;
	CBQBBitSet bscu(m_IF);
	unsigned long t0, rxyz=1;


	if (parm == NULL) {
		m_IF.eprintf("CBQBDriver::CompressCube(): Error: parm == NULL.\n");
		return false;
	}

	resi = 0;

	t0 = (unsigned long)time(NULL);

	if (parm->GetUseExtra())
		if (parm->GetExtraTRange() > parm->GetOrder())
			parm->SetOrder(parm->GetExtraTRange());

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("Starting process...\n");

	m_pEngine->m_pReadCacheCube->RewindReadPos();

	if (parm->GetUseExtra())
		mo = parm->GetExtraTRange();
	else
		mo = parm->GetOrder();

	i = 0;
	tb = 0;
	tf = 0;
	while (true) {

		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			m_IF.printf("\n");
			m_IF.printf("\n");
		}

		if (m_IF.IsPL(BQB_PL_VERBOSE)) {
			if (i != 0)
				m_IF.printf("    Step %6d ... (avg %5.1f s/step)\n",
					i+1,((double)(time(NULL)-t0))/i);
			else
				m_IF.printf("    Step %6d ...\n",i+1);
		}

		cfr = m_pEngine->m_pReadCacheCube->GetNextFrame();

		if (cfr == NULL)
			break;

		if (i == 0) {

			rxyz = cfr->m_iResXYZ;

			if (parm->GetUseExtra()) {

				if (parm->GetExtraPredCorr()) {

					if (m_IF.IsPL(BQB_PL_VERBOSE))
						m_IF.printf("      Initializing Predictor Extrapolator...\n");

					m_pEngine->m_pExtrapolator = new CBQBExtrapolator(m_IF);
					m_pEngine->m_pExtrapolator->Initialize(
						cfr->m_iRes[0],  // resx,
						cfr->m_iRes[1],  // resy,
						cfr->m_iRes[2],  // resz,
						1,               // srangex,
						1,               // srangey,
						1,               // srangez,
						parm->GetExtraTRange(),    // trange,
						0,               // sorder,
						parm->GetExtraTOrder(),    // torder,
						0,               // offsetx,
						0,               // offsety,
						0,               // offsetz,
						false,           // crosss,
						false,           // crosst,
						false,           // wrap,
						false,           // crossranges,
						false,           // crossranget
						0,               // distexpo
						parm->GetExtraTimeExpo(),  // timeexpo
						!m_IF.IsPL(BQB_PL_VERBOSE)
					);

					if (m_IF.IsPL(BQB_PL_VERBOSE))
						m_IF.printf("      Initializing Corrector Extrapolator...\n");

					m_pEngine->m_pExtrapolatorCorr = new CBQBExtrapolator(m_IF);
					m_pEngine->m_pExtrapolatorCorr->Initialize(
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
						!m_IF.IsPL(BQB_PL_VERBOSE)
					);

				} else {

					if (m_IF.IsPL(BQB_PL_VERBOSE))
						m_IF.printf("      Initializing Extrapolator...\n");

					m_pEngine->m_pExtrapolator = new CBQBExtrapolator(m_IF);
					m_pEngine->m_pExtrapolator->Initialize(
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
						!m_IF.IsPL(BQB_PL_VERBOSE)
					);
				}
			}
		}

		if (m_pEngine->m_pExtrapolator != NULL)
			m_pEngine->m_pExtrapolator->PushCubeFrame(cfr);

		if (i < mo)
			co = i;
		else
			co = mo;

		if ((i+1 > mo) || (m_iOptIncludeFirst > 0)) {

			bscu.Clear();

			m_pEngine->CompressCubeFrame(
				&bscu,
				co,
				parm,
				histused,
				(i==0), // Store static info
				NULL,
				!realsize,
				&tresi,
				tciaa
			);

			//printf("@ i=%d: %f + %f --> %f\n",i,resi,tresi,resi+tresi);

			resi += tresi;

			tf++;

			if (realsize) {
				tb += bscu.GetByteLength();
/*				if (m_fResiCorrFrame != NULL) {
					mfprintf(m_fResiCorrFrame,"%.3f;  %.3f; %d\n", bscu.GetByteLength()/1024.0, mypow( tresi / rxyz, 1.0/m_fEntropyExpoCube ), m_iHistogramCounter-1 );
					fflush(m_fResiCorrFrame);
				}*/
			}
		}

		i++;
		if (m_pEngine->m_pReadCacheCube->GetCachedSteps() == 0)
			break;
	}

	if (m_pEngine->m_pExtrapolator != NULL) {
		delete m_pEngine->m_pExtrapolator;
		m_pEngine->m_pExtrapolator = NULL;
	}

	if (m_pEngine->m_pExtrapolatorCorr != NULL) {
		delete m_pEngine->m_pExtrapolatorCorr;
		m_pEngine->m_pExtrapolatorCorr = NULL;
	}

	if (m_pEngine->m_pExtrapolatorXYZ != NULL) {
		delete m_pEngine->m_pExtrapolatorXYZ;
		m_pEngine->m_pExtrapolatorXYZ = NULL;
	}

//	printf("@ tf=%f, resi=%f, frac=%f, pow2=%f\n",tf,resi, resi / (tf * rxyz), bqbpow2( resi / (tf * rxyz) ) );

	resi = bqbpow2( resi / (tf * rxyz) );

	if (presi != NULL)
		*presi = resi;

	if (realsize) {
/*		if (m_fResiCorrTraj != NULL) {
			mfprintf(m_fResiCorrTraj,"%.4f;  %.4f\n", tb / tf / 1024.0, resi );
			fflush(m_fResiCorrTraj);
		}*/
		return tb / tf / 1024.0;
	} else
		return resi;
}


double CBQBDriver::OptimizeCube_ObjectiveTExpo(double texpo) {

	double entropy, resi;

	m_pVolParm->SetExtraTimeExpo(texpo);

	entropy = CompressCubeBenchmark(
		m_pVolParm,
		m_bOptCubeRealSize,     // bool realsize
		&resi,                   // double *presi
		NULL                     // std::vector<std::vector<int> > *tciaa
	);

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		if (m_bOptCubeRealSize)
			m_IF.printf("        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f (%10.3f kiB)\n",
				m_pVolParm->GetExtraSRangeX(),
				m_pVolParm->GetExtraSOrder(),
				m_pVolParm->GetExtraTRange(),
				m_pVolParm->GetExtraTOrder(),
				m_pVolParm->GetExtraDistExpo(),
				texpo,
				resi,
				entropy
			);
		else
			m_IF.printf("        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f\n",
				m_pVolParm->GetExtraSRangeX(),
				m_pVolParm->GetExtraSOrder(),
				m_pVolParm->GetExtraTRange(),
				m_pVolParm->GetExtraTOrder(),
				m_pVolParm->GetExtraDistExpo(),
				texpo,
				entropy
			);
	}

	return entropy;
}


double CBQBDriver::OptimizeCube_ObjectiveSExpo(double sexpo) {

	double entropy, resi;

	m_pVolParm->SetExtraDistExpo(sexpo);

	entropy = CompressCubeBenchmark(
		m_pVolParm,
		m_bOptCubeRealSize,     // bool realsize
		&resi,                   // double *presi
		NULL                     // std::vector<std::vector<int> > *tciaa
	);

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		if (m_bOptCubeRealSize)
			m_IF.printf("        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f (%10.3f kiB)\n",
				m_pVolParm->GetExtraSRangeX(),
				m_pVolParm->GetExtraSOrder(),
				m_pVolParm->GetExtraTRange(),
				m_pVolParm->GetExtraTOrder(),
				sexpo,
				m_pVolParm->GetExtraTimeExpo(),
				resi,
				entropy
			);
		else
			m_IF.printf("        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f\n",
				m_pVolParm->GetExtraSRangeX(),
				m_pVolParm->GetExtraSOrder(),
				m_pVolParm->GetExtraTRange(),
				m_pVolParm->GetExtraTOrder(),
				sexpo,
				m_pVolParm->GetExtraTimeExpo(),
				entropy
			);
	}

	return entropy;
}


class VARPARAMS{
public:
	VARPARAMS(int srx, int sry, int srz, int so, int ofx, int ofy, int ofz, int tr, int to, double de, double te)
		: srangex(srx), srangey(sry), srangez(srz), sorder(so), offsetx(ofx), offsety(ofy), offsetz(ofz),
		  trange(tr), torder(to), distexp(de), timeexp(te) { }
	bool operator < (const VARPARAMS& v) const { UNUSED(v); return true; }
	int srangex;
	int srangey;
	int srangez;
	int sorder;
	int offsetx;
	int offsety;
	int offsetz;
	int trange;
	int torder;
	double distexp;
	double timeexp;
};


bool CBQBDriver::OptimizeCubeParameters(
		int     level,
		CBQBParameterSet_Volumetric *parm
	) {

	int steps, z, z2, z3, ti, zsr,zso, btr, bto, bsrx, bsry, bsrz, bsofx, bsofy, bsofz, bso, ca, bi;
	int zsrx, zsry, zsrz, zsofx, zsofy, zsofz;
	int tia[12];
	double tf, tflast, btv, resi, tfa[12], tfa2[12], fac, ode;
	bool realsize;
	unsigned long t0;
	const CBQBCubeFrame *cfr;
	char buf[256];
	std::vector<std::pair<double,std::string> > rlist;
	std::vector<std::pair<double,VARPARAMS> > rlist2, rlist3;
	std::vector<std::vector<int> > tciaa;
	std::vector<int> iatc, iabl;
	CBQBBitSet bs(m_IF);
	CBQBIntegerEngine *tie;


	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.bprintf("    >>> Volumetric Trajectory Parameter Optimization >>>\n");

/*	if (level == 3) {
		m_fResiCorrFrame = OpenFileWrite("corrframecube.csv",true);
		m_fResiCorrTraj = OpenFileWrite("corrtrajcube.csv",true);
	}*/

	if (level >= 3)
		realsize = true;
	else
		realsize = false;

	t0 = (unsigned long)time(NULL);

	m_pEngine->m_pReadCacheCube->RewindReadPos();

	cfr = m_pEngine->m_pReadCacheCube->GetNextFrame();

	if (parm->GetHilbert())
		m_pEngine->BuildHilbertIdx(cfr->m_iRes[0],cfr->m_iRes[1],cfr->m_iRes[2]/*,verbose*/);

	m_pEngine->m_pReadCacheCube->RewindReadPos();

	steps = m_pEngine->m_pReadCacheCube->GetCachedSteps();

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.printf("        Optimization level %d\n",level);
		m_IF.printf("        Have %d cached steps.\n",steps);

		m_IF.printf("        SR  SO  TR  TO  DistExp  TimeExp        Residuum\n");
	}

	tflast = 0;
	btv = 1.0e30;
	tf = 0;
	btr = -1;
	bto = -1;
	bsrx = -1;
	bsry = -1;
	bsrz = -1;
	bsofx = -1;
	bsofy = -1;
	bsofz = -1;
	bso = -1;

	parm->SetUseExtra(true);
	//cextratimeexpo = 1.0;
	//cextradistexpo = 3.0;
	parm->SetExtraSOrder(0);
	parm->SetExtraSRange(1);
	parm->SetExtraOffset(0);

	if (level == 2)
		parm->SetTableCount(1);

	/********************* Separate Code Block for Level 1 *************************************************/
	if (level == 1) {

		rlist2.clear();

		switch(steps) {
			case 1: //                                                       RX RY RZ SO OX OY OZ TR TO  SE   TE
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 5, 3, 3, 3, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 6, 3, 3, 3, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 5, 7, 7, 3, 2, 2, 2, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 5, 5, 5, 2, 2, 2, 2, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 3, 3, 3, 2, 1, 1, 1, 1, 0, 3.5, 1.5 ) ) );
				break;
			case 2: //                                                       RX RY RZ SO OX OY OZ TR TO  SE   TE
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 2, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 5, 3, 3, 3, 2, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 5, 3, 3, 3, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 5, 5, 5, 3, 2, 2, 2, 2, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 3, 3, 3, 2, 1, 1, 1, 2, 0, 3.5, 1.5 ) ) );
				break;
			case 3: //                                                       RX RY RZ SO OX OY OZ TR TO  SE   TE
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 3, 1, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 2, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 5, 3, 3, 3, 2, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 5, 3, 3, 3, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 5, 5, 5, 3, 2, 2, 2, 3, 1, 3.5, 1.5 ) ) );
				break;
			case 4: //                                                       RX RY RZ SO OX OY OZ TR TO  SE   TE
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 4, 1, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 2, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 5, 3, 3, 3, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 5, 3, 3, 3, 2, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 5, 5, 5, 3, 2, 2, 2, 4, 1, 3.5, 1.5 ) ) );
				break;
			default: //                                                      RX RY RZ SO OX OY OZ TR TO  SE   TE
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 5, 1, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 2, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 4, 3, 3, 3, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 5, 3, 3, 3, 1, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 7, 7, 7, 5, 3, 3, 3, 2, 0, 3.5, 1.5 ) ) );
				rlist2.push_back( std::pair<double, VARPARAMS>( 0, VARPARAMS( 5, 5, 5, 3, 2, 2, 2, 4, 1, 3.5, 1.5 ) ) );
				break;
		}

		for (z=0;z<(int)rlist2.size();z++) {

			parm->SetExtraSRangeX(rlist2[z].second.srangex);
			parm->SetExtraSRangeY(rlist2[z].second.srangey);
			parm->SetExtraSRangeZ(rlist2[z].second.srangez);
			parm->SetExtraTRange(rlist2[z].second.trange);
			parm->SetExtraSOrder(rlist2[z].second.sorder);
			parm->SetExtraTOrder(rlist2[z].second.torder);
			parm->SetExtraOffsetX(rlist2[z].second.offsetx);
			parm->SetExtraOffsetY(rlist2[z].second.offsety);
			parm->SetExtraOffsetZ(rlist2[z].second.offsetz);
			parm->SetExtraDistExpo(rlist2[z].second.distexp);
			parm->SetExtraTimeExpo(rlist2[z].second.timeexp);

			tf = CompressCubeBenchmark(
				parm,
				false,                     // bool realsize
				&resi,                     // double *presi
				NULL                       // std::vector<std::vector<int> > *tciaa
			);

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f",
					rlist2[z].second.srangex,rlist2[z].second.sorder,rlist2[z].second.trange,
					rlist2[z].second.torder,rlist2[z].second.distexp,rlist2[z].second.timeexp,tf);

			if (tf < btv) {
				btv = tf;
				btr = z;
				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("  <---");
			}

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("\n");
		}

		parm->SetExtraSRangeX(rlist2[btr].second.srangex);
		parm->SetExtraSRangeY(rlist2[btr].second.srangey);
		parm->SetExtraSRangeZ(rlist2[btr].second.srangez);
		parm->SetExtraTRange(rlist2[btr].second.trange);
		parm->SetExtraSOrder(rlist2[btr].second.sorder);
		parm->SetExtraTOrder(rlist2[btr].second.torder);
		parm->SetExtraOffsetX(rlist2[btr].second.offsetx);
		parm->SetExtraOffsetY(rlist2[btr].second.offsety);
		parm->SetExtraOffsetZ(rlist2[btr].second.offsetz);
		parm->SetExtraDistExpo(rlist2[btr].second.distexp);
		parm->SetExtraTimeExpo(rlist2[btr].second.timeexp);

		tflast = parm->GetExtraTimeExpo();
		ode = parm->GetExtraDistExpo();

		// Fix MB 11.04.2020
		if (parm->GetExtraTRange() > parm->GetExtraTOrder()+2) {
//		if (parm->GetExtraTRange() != parm->GetExtraTOrder()+2) {

			m_IF.printf( "\n" );

			tfa[0] = 0.5;
			tfa[1] = 1.0;

			for (z=0;z<2;z++) {

				parm->SetExtraTimeExpo(tfa[z]);

				tfa2[z] = CompressCubeBenchmark(
					parm,
					false,              // bool realsize
					&resi,              // double *presi
					NULL                // std::vector<std::vector<int> > *tciaa
				);

				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f\n",
						parm->GetExtraSRangeX(),
						parm->GetExtraSOrder(),
						parm->GetExtraTRange(),
						parm->GetExtraTOrder(),
						parm->GetExtraDistExpo(),
						tfa[z],
						tfa2[z]
					);
			}

			if ((tfa2[0] < tfa2[1]) && (tfa2[0] < btv)) {
				tflast = tfa[0];
				btv = tfa2[0];
			} else if ((tfa2[1] < tfa2[0]) && (tfa2[1] < btv)) {
				tflast = tfa[1];
				btv = tfa2[1];
			}
		}

		parm->SetExtraTimeExpo(tflast);

		tfa[0] = 4.5;
		tfa[1] = 6.0;

		m_IF.printf( "\n" );

		for (z=0;z<2;z++) {

			parm->SetExtraDistExpo(tfa[z]);

			tfa2[z] = CompressCubeBenchmark(
				parm,
				false,              // bool realsize
				&resi,              // double *presi
				NULL                // std::vector<std::vector<int> > *tciaa
			);

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f\n",
					parm->GetExtraSRangeX(),
					parm->GetExtraSOrder(),
					parm->GetExtraTRange(),
					parm->GetExtraTOrder(),
					tfa[z],
					parm->GetExtraTimeExpo(),
					tfa2[z]
				);
		}

		//printf("@ btv=%f, tflast=%f, tfa2[0]=%f, tfa2[1]=%f\n",btv,tflast,tfa2[0],tfa2[1]);

		parm->SetExtraDistExpo(ode);

		if ((tfa2[0] < tfa2[1]) && (tfa2[0] < btv))
			parm->SetExtraDistExpo(tfa[0]);
		else if ((tfa2[1] < tfa2[0]) && (tfa2[1] < btv))
			parm->SetExtraDistExpo(tfa[1]);

		goto _end;
	}
	/********************* End of Separate Code Block for Level 1 ******************************************/

	for (z=0;z<3;z++)
		tfa[z] = 1.0e30;

	for (z=0;z<3;z++) {

		if (z+1 > steps)
			break;

		parm->SetExtraSRange(1);
		parm->SetExtraTRange(z+1);
		parm->SetExtraSOrder(0);
		parm->SetExtraTOrder(MAX(z-1,0));
		parm->SetExtraOffset(0);

		tfa[z] = CompressCubeBenchmark(
			parm,
			realsize,           // bool realsize
			&resi,              // double *presi
			NULL                // std::vector<std::vector<int> > *tciaa
		);

		if (m_IF.IsPL(BQB_PL_STANDARD)) {
			if (realsize)
				sprintf(buf,"        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f (%10.3f kiB)",
					1,
					0,
					z+1,
					MAX(z-1,0),
					parm->GetExtraDistExpo(),
					parm->GetExtraTimeExpo(),
					resi,
					tfa[z]
				);
			else
				sprintf(buf,"        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f",
					1,
					0,
					z+1,
					MAX(z-1,0),
					parm->GetExtraDistExpo(),
					parm->GetExtraTimeExpo(),
					tfa[z]
				);

			m_IF.printf("%s",(const char*)buf);
		} else
			buf[0] = 0;

		if (tfa[z] < btv) {
			btv = tfa[z];
			bto = MAX(z-1,0);
			btr = z+1;
			bsrx = 1;
			bsry = 1;
			bsrz = 1;
			bso = 0;
			bsofx = 0;
			bsofy = 0;
			bsofz = 0;
			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("  <---");
		}

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("\n");

		rlist.push_back(std::pair<double,std::string>(tfa[z],(const char*)buf));

		rlist2.push_back(std::pair<double,VARPARAMS>(tfa[z],VARPARAMS(
			1,
			1,
			1,
			0,
			0,
			0,
			0,
			z+1,
			MAX(z-1,0),
			parm->GetExtraDistExpo(),
			parm->GetExtraTimeExpo()
		)));

	}

	if ((tfa[2]+0.001 < tfa[0]) && (tfa[2]+0.001 < tfa[1])) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("        ---> Time-continuous case.\n");
		ca = 0;
	} else if ((tfa[1]+0.001 < tfa[0]) && (tfa[1]+0.001 < tfa[2])) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("        ---> Semi-time-continuous case.\n");
		ca = 1;
	} else {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("        ---> Time-discontinuous case.\n");
		ca = 2;
	}

	if (level == 4) {

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("        SRX  SRY  SRZ  OFX  OFY  OFZ  SO  TR  TO  DistExp  TimeExp        Residuum\n");

		bi = 0;

		for (zsrx=1;zsrx<10;zsrx++) {

			for (zsry=1;zsry<10;zsry++) {

				for (zsrz=1;zsrz<10;zsrz++) {

					if (MAX3(zsrx,zsry,zsrz) - MIN3(zsrx,zsry,zsrz) > 2)
						continue;

					for (zsofx=0;zsofx<zsrx;zsofx++) {

						for (zsofy=0;zsofy<zsry;zsofy++) {

							for (zsofz=0;zsofz<zsrz;zsofz++) {

								if (MAX3(zsrx,zsry,zsrz) <= 3) {
									if (MAX3(abs(zsofx-(zsrx/2)),abs(zsofy-(zsry/2)),abs(zsofz-(zsrz/2))) > 0)
										continue;
								} else if (MAX3(zsrx,zsry,zsrz) <= 5) {
									if (MAX3(abs(zsofx-(zsrx/2)),abs(zsofy-(zsry/2)),abs(zsofz-(zsrz/2))) > 1)
										continue;
								} else {
									if (MAX3(abs(zsofx-(zsrx/2)),abs(zsofy-(zsry/2)),abs(zsofz-(zsrz/2))) > 2)
										continue;
								}

								for (zso=0;zso<MIN(MAX3(zsrx,zsry,zsrz),9);zso++) {

									if (m_IF.IsPL(BQB_PL_STANDARD))
										if (bi > 1)
											m_IF.printf("\n");

									bi = 0;

									if (((zsrx >= 5) || (zsry >= 5) || (zsrz >= 5)) && (zso == 0))
										continue;

									if (((zsrx >= 7) || (zsry >= 7) || (zsrz >= 7)) && (zso < 2))
										continue;

									if (((zsrx == 9) || (zsry == 9) || (zsrz == 9)) && (zso < 3))
										continue;

									for (z=1;z<MIN(steps+1,7);z++) { // trange

										for (z2=z-2;z2>=-1;z2--) { // torder

											if ((z > 1) && (z2 == -1))
												continue;
											if ((z > 3) && (z2 == 0))
												continue;

											// The three cases from above
											if ((zsrx == 1) && (zsry == 1) && (zsrz == 1) && (z == 1) && (z2 == -1))
												continue;
											if ((zsrx == 1) && (zsry == 1) && (zsrz == 1) && (z == 2) && (z2 == 0))
												continue;
											if ((zsrx == 1) && (zsry == 1) && (zsrz == 1) && (z == 3) && (z2 == 1))
												continue;

											if (ca == 0) {

												if ((z2 < 0) || ((z2 < 1) && (z > 3)))
													continue;

											} else if (ca == 1) {

												if ((z2 < 0) || (z2 > 1))
													continue;

											} else {

												if (z2 > -1)
													continue;
											}

											tflast = tf;

											parm->SetExtraSRangeX(zsrx);
											parm->SetExtraSRangeY(zsry);
											parm->SetExtraSRangeZ(zsrz);
											parm->SetExtraTRange(z);
											parm->SetExtraSOrder(zso);
											parm->SetExtraTOrder(z2);
											parm->SetExtraOffsetX(zsofx);
											parm->SetExtraOffsetY(zsofy);
											parm->SetExtraOffsetZ(zsofz);

											tf = CompressCubeBenchmark(
												parm,
												realsize,           // bool realsize
												&resi,              // double *presi
												NULL                // std::vector<std::vector<int> > *tciaa
											);

											if (m_IF.IsPL(BQB_PL_STANDARD)) {
												if (realsize)
													sprintf(buf,"        %3d  %3d  %3d  %3d  %3d  %3d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f (%10.3f kiB)",
														zsrx,
														zsry,
														zsrz,
														zsofx,
														zsofy,
														zsofz,
														zso,
														z,
														z2,
														parm->GetExtraDistExpo(),
														parm->GetExtraTimeExpo(),
														resi,
														tf
													);
												else
													sprintf(buf,"        %3d  %3d  %3d  %3d  %3d  %3d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f",
														zsrx,
														zsry,
														zsrz,
														zsofx,
														zsofy,
														zsofz,
														zso,
														z,
														z2,
														parm->GetExtraDistExpo(),
														parm->GetExtraTimeExpo(),
														tf
													);

												m_IF.printf("%s",(const char*)buf);
											} else
												buf[0] = 0;

											bi++;

											rlist.push_back(std::pair<double,std::string>(tf,(const char*)buf));

											rlist2.push_back(std::pair<double,VARPARAMS>(tf,VARPARAMS(
												zsrx,
												zsry,
												zsrz,
												zso,
												zsofx,
												zsofy,
												zsofz,
												z,
												z2,
												parm->GetExtraDistExpo(),
												parm->GetExtraTimeExpo()
											)));

											if (tf < btv) {
												btv = tf;
												bto = z2;
												btr = z;
												bsrx = zsrx;
												bsry = zsry;
												bsrz = zsrz;
												bsofx = zsofx;
												bsofy = zsofy;
												bsofz = zsofz;
												bso = zso;
												if (m_IF.IsPL(BQB_PL_STANDARD))
													m_IF.printf("  <---");
											}

											if (m_IF.IsPL(BQB_PL_STANDARD))
												m_IF.printf("\n");

											if ((z2 < z-2) && (tflast < tf))
												break;
										}
									}
								}
							}
						}
					}
				}
			}
		}


	} else {


		bi = 0;

		for (zsr=1;zsr<10;zsr+=2) {

			if ((zsr == 9) && (level < 3))
				break;

			for (zso=0;zso<MIN(zsr,9);zso++) {

				if (m_IF.IsPL(BQB_PL_STANDARD))
					if (bi > 1)
						m_IF.printf("\n");

				bi = 0;

				if ((zsr >= 5) && (zso == 0))
					continue;

				if ((zsr >= 7) && (zso < 2))
					continue;

				if ((zsr == 9) && (zso < 3))
					continue;

				for (z=1;z<MIN(steps+1,7);z++) { // trange

					if ((level != 3) && (z > 2)) {
						if (zso > 4)
							continue;
						if ((zsr == 5) && (zso == 4))
							continue;
					}

					for (z2=z-2;z2>=-1;z2--) { // torder

						if ((z > 1) && (z2 == -1))
							continue;
						if ((z > 3) && (z2 == 0))
							continue;

						// The three cases from above
						if ((zsr == 1) && (z == 1) && (z2 == -1))
							continue;
						if ((zsr == 1) && (z == 2) && (z2 == 0))
							continue;
						if ((zsr == 1) && (z == 3) && (z2 == 1))
							continue;

						if (ca == 0) {

							if ((z2 < 0) || ((z2 < 1) && (z > 3)))
								continue;

						} else if (ca == 1) {

							if ((z2 < 0) || (z2 > 1))
								continue;

						} else {

							if (z2 > -1)
								continue;
						}

				//		printf("@@ srange=%d, sorder=%d, trange=%d, torder=%d, offset=%d\n",zsr,zso,z,z2,zsr/2);

						tflast = tf;

						// Bugfix MB 02.12.2020
						//parm->SetExtraSRangeX(zsr);
						parm->SetExtraSRange(zsr);
						parm->SetExtraTRange(z);
						parm->SetExtraSOrder(zso);
						parm->SetExtraTOrder(z2);
						parm->SetExtraOffset(zsr/2);

						tf = CompressCubeBenchmark(
							parm,
							realsize,           // bool realsize
							&resi,              // double *presi
							NULL                // std::vector<std::vector<int> > *tciaa
						);

						if (m_IF.IsPL(BQB_PL_STANDARD)) {
							if (realsize)
								sprintf(buf,"        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f (%10.3f kiB)",
									zsr,
									zso,
									z,
									z2,
									parm->GetExtraDistExpo(),
									parm->GetExtraTimeExpo(),
									resi,
									tf
								);
							else
								sprintf(buf,"        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f",
									zsr,
									zso,
									z,
									z2,
									parm->GetExtraDistExpo(),
									parm->GetExtraTimeExpo(),
									tf
								);

							m_IF.printf("%s",(const char*)buf);
						} else
							buf[0] = 0;

						bi++;

						rlist.push_back(std::pair<double,std::string>(tf,(const char*)buf));

						rlist2.push_back(std::pair<double,VARPARAMS>(tf,VARPARAMS(
							zsr,
							zsr,
							zsr,
							zso,
							zsr/2,
							zsr/2,
							zsr/2,
							z,
							z2,
							parm->GetExtraDistExpo(),
							parm->GetExtraTimeExpo()
						)));

						if (tf < btv) {
							btv = tf;
							bto = z2;
							btr = z;
							bsrx = zsr;
							bsry = zsr;
							bsrz = zsr;
							bsofx = zsr/2;
							bsofy = zsr/2;
							bsofz = zsr/2;
							bso = zso;
							if (m_IF.IsPL(BQB_PL_STANDARD))
								m_IF.printf("  <---"); 
						}

						if (m_IF.IsPL(BQB_PL_STANDARD))
							m_IF.printf("\n");

						if ((z2 < z-2) && (tflast < tf))
							break;
					}
				}
			}
		}
	}

	std::sort(rlist.begin(),rlist.end());

/*	FILE *a;
	a = OpenFileWrite("opticube_list.txt",true);
	mfprintf(a,"        sr  so  tr  to  distexp  timeexp       residuum\n");
	for (z=0;z<(int)rlist.size();z++)
		mfprintf(a,"%s\n",rlist[z].second.c_str());
	fclose(a);*/

	if (level == 2) {

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("        Evaluating best 10:\n");

		std::sort(rlist2.begin(),rlist2.end());

		for (z=0;z<10;z++)
			rlist3.push_back(rlist2[z]);

		//int tia[4];

		for (z=0;z<4;z++)
			tia[z] = 0;

		for (z=0;z<(int)rlist3.size();z++)
			tia[rlist3[z].second.srangex/2]++;

		// Make sure that at least 2 sets from each srange are included
		for (z=0;z<4;z++) {
			if ((tia[z] < 1) || ((z > 0) && (tia[z] < 2))) {
				for (z2=10;z2<(int)rlist2.size();z2++) {
					if ((rlist2[z2].second.srangex/2) == z) {
						rlist3.push_back(rlist2[z2]);
						tia[z]++;
						if ((tia[z] >= 2) || ((z == 0) && (tia[z] == 1)))
							break;
					}
				}
			}
		}

		// Crop the total number back to 10
		if (rlist3.size() > 10) {
			for (z=(int)rlist3.size()-1;z>=0;z--) {
				if ((tia[rlist3[z].second.srangex/2] > 2) || ((rlist3[z].second.srangex == 1) && (tia[rlist3[z].second.srangex/2] > 1))) {
					tia[rlist3[z].second.srangex/2]--;
					rlist3.erase(rlist3.begin()+z);
				}
				if (rlist3.size() == 10)
					break;
			}
		}

		std::sort(rlist3.begin(),rlist3.end());

		btr = -1;
		btv = 1.0e30;
		for (z=0;z<(int)rlist3.size();z++) {

			parm->SetExtraSRangeX(rlist3[z].second.srangex);
			parm->SetExtraSRangeY(rlist3[z].second.srangey);
			parm->SetExtraSRangeZ(rlist3[z].second.srangez);
			parm->SetExtraTRange(rlist3[z].second.trange);
			parm->SetExtraSOrder(rlist3[z].second.sorder);
			parm->SetExtraTOrder(rlist3[z].second.torder);
			parm->SetExtraOffsetX(rlist3[z].second.offsetx);
			parm->SetExtraOffsetY(rlist3[z].second.offsety);
			parm->SetExtraOffsetZ(rlist3[z].second.offsetz);
			parm->SetExtraDistExpo(rlist3[z].second.distexp);
			parm->SetExtraTimeExpo(rlist3[z].second.timeexp);

			tf = CompressCubeBenchmark(
				parm,
				true,                      // bool realsize
				&resi,                     // double *presi
				NULL                       // std::vector<std::vector<int> > *tciaa
			);

			if (m_IF.IsPL(BQB_PL_STANDARD)) {
				sprintf(buf,"        %2d  %2d  %2d  %2d   %6.3f   %6.3f  %14.4f (%10.3f kiB)",
					rlist3[z].second.srangex,rlist3[z].second.sorder,rlist3[z].second.trange,
					rlist3[z].second.torder,rlist3[z].second.distexp,rlist3[z].second.timeexp,resi,tf);

				m_IF.printf("%s",(const char*)buf);
			} else
				buf[0] = 0;

			if (tf < btv) {
				btv = tf;
				btr = z;
				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("  <---");
			}

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("\n");
		}

		parm->SetExtraSRangeX(rlist3[btr].second.srangex);
		parm->SetExtraSRangeY(rlist3[btr].second.srangey);
		parm->SetExtraSRangeZ(rlist3[btr].second.srangez);
		parm->SetExtraTRange(rlist3[btr].second.trange);
		parm->SetExtraSOrder(rlist3[btr].second.sorder);
		parm->SetExtraTOrder(rlist3[btr].second.torder);
		parm->SetExtraOffsetX(rlist3[btr].second.offsetx);
		parm->SetExtraOffsetY(rlist3[btr].second.offsety);
		parm->SetExtraOffsetZ(rlist3[btr].second.offsetz);

	} else {

		parm->SetExtraSRangeX(bsrx);
		parm->SetExtraSRangeY(bsry);
		parm->SetExtraSRangeZ(bsrz);
		parm->SetExtraTRange(btr);
		parm->SetExtraSOrder(bso);
		parm->SetExtraTOrder(bto);
		parm->SetExtraOffsetX(bsofx);
		parm->SetExtraOffsetY(bsofy);
		parm->SetExtraOffsetZ(bsofz);
	}

	if (level >= 2) {

		m_pVolParm = parm;

		m_bOptCubeRealSize = realsize;

		// Fix MB 11.04.2020
		if (parm->GetExtraTOrder()+2 < parm->GetExtraTRange()) {
		//if (parm->GetExtraTOrder()+2 != parm->GetExtraTRange()) {

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("        Optimizing time exponent:\n");

			tf = GoldenSectionSearch(4.0, (level>=3)?0.02:0.2, &CBQBDriver::OptimizeCube_ObjectiveTExpo);

			if (level >= 3)
 				parm->SetExtraTimeExpo( floor(tf*1000.0+0.5)/1000.0 );
			else
				parm->SetExtraTimeExpo( floor(tf*100.0+0.5)/100.0 );

		} else
			parm->SetExtraTimeExpo( 0 );

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("        Optimizing distance exponent:\n");

		tf = GoldenSectionSearch(8.0, (level>=3)?0.02:0.2, &CBQBDriver::OptimizeCube_ObjectiveSExpo);

		if (level >= 3)
			parm->SetExtraDistExpo( floor(tf*1000.0+0.5)/1000.0 );
		else
			parm->SetExtraDistExpo( floor(tf*100.0+0.5)/100.0 );


		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("        Optimizing Compressor settings:\n");

		tciaa.clear();

		tf = CompressCubeBenchmark(
			parm,
			true,               // bool realsize
			NULL,               // double *presi
			&tciaa              // std::vector<std::vector<int> > *tciaa
		);

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("        TC  BL   BW  MTF  Frame Size [kiB]\n");

		if (level >= 4) {

			iatc.push_back(1);
			iatc.push_back(2);
			iatc.push_back(3);
			iatc.push_back(4);
			iatc.push_back(6);
			iatc.push_back(8);
			iatc.push_back(10);
			iatc.push_back(12);
			iatc.push_back(15);
			iatc.push_back(18);
			iatc.push_back(21);
			iatc.push_back(24);

			iabl.push_back(10);
			iabl.push_back(12);
			iabl.push_back(14);
			iabl.push_back(16);
			iabl.push_back(20);
			iabl.push_back(24);
			iabl.push_back(28);
			iabl.push_back(32);
			iabl.push_back(36);
			iabl.push_back(40);
			iabl.push_back(48);
			iabl.push_back(56);
			iabl.push_back(64);

			btv = 1.0e30;
			btr = -1;

			for (z=0;z<(int)iatc.size();z++) {

				for (z2=0;z2<(int)iabl.size();z2++) {

					tf = 0;
					for (z3=0;z3<(int)tciaa.size();z3++) {

						bs.Clear();
						tie = new CBQBIntegerEngine(m_IF);

						parm->SetBlockLength(iabl[z2]);
						parm->SetTableCount(iatc[z]);

						tie->Compress(
							tciaa[z3],
							&bs,
							parm->GetCompressorParameterSet(),
							false,
							NULL
						);

						delete tie;
						tf += bs.GetByteLength();
					}

					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("        %2d  %2d  %3s  %3s  %10.4f",
							iatc[z],
							iabl[z2],
							parm->GetBW()?"yes":"no",
							parm->GetMTF()?"yes":"no",
							tf/tciaa.size()/1024.0
						);

					if (tf < btv) {
						btv = tf;
						parm->SetTableCount(iatc[z]);
						parm->SetBlockLength(iabl[z2]);
						if (m_IF.IsPL(BQB_PL_STANDARD))
							m_IF.printf("  <---");
					}
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("\n");

					if (iatc[z] == 1)
						break; // One table: Block length does not matter
				}
			}

		} else {
			
			if (level >= 3) {

				ti = 11;
				tia[0] = 28;
				tia[1] = 24;
				tia[2] = 20;
				tia[3] = 16;
				tia[4] = 12;
				tia[5] = 8;
				tia[6] = 6;
				tia[7] = 4;
				tia[8] = 3;
				tia[9] = 2;
				tia[10] = 1;
			} else {
				ti = 4;
				tia[0] = 10;
				tia[1] = 6;
				tia[2] = 3;
				tia[3] = 1;
			}

			btv = 1.0e30;
			btr = -1;
			for (z=0;z<ti;z++) {

				tfa[z] = 0;
				for (z2=0;z2<(int)tciaa.size();z2++) {

					bs.Clear();
					tie = new CBQBIntegerEngine(m_IF);

					parm->SetTableCount(tia[z]);

					tie->Compress(
						tciaa[z2],
						&bs,
						parm->GetCompressorParameterSet(),
						false,
						NULL
					);

					delete tie;
					tfa[z] += bs.GetByteLength();
				}

				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("        %2d  %2d  %3s  %3s  %10.3f",
						tia[z],
						parm->GetBlockLength(),
						parm->GetBW()?"yes":"no",
						parm->GetMTF()?"yes":"no",
						tfa[z]/1024.0/tciaa.size()
					);

				if (tfa[z] < btv) {
					btv = tfa[z];
					btr = z;
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("  <---");
				}
				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("\n");

				if ((z > 0) && (tfa[z] > tflast))
					break;

				tflast = tfa[z];
			}
			parm->SetTableCount(tia[btr]);

			if (parm->GetTableCount() != 1) {

				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("\n");

				if (level >= 3) {
					ti = 10;
					tia[0] = 10;
					tia[1] = 12;
					tia[2] = 14;
					tia[3] = 16;
					tia[4] = 20;
					tia[5] = 24;
					tia[6] = 28;
					tia[7] = 32;
					tia[8] = 36;
					tia[9] = 40;
				} else {
					ti = 4;
					tia[0] = 15;
					tia[1] = 20;
					tia[2] = 30;
					tia[3] = 40;
				}

				btv = 1.0e30;
				btr = -1;
				for (z=0;z<ti;z++) {

					tfa[z] = 0;
					for (z2=0;z2<(int)tciaa.size();z2++) {

						bs.Clear();
						tie = new CBQBIntegerEngine(m_IF);

						parm->SetBlockLength(tia[z]);

						tie->Compress(
							tciaa[z2],
							&bs,
							parm->GetCompressorParameterSet(),
							false,
							NULL
						);

						delete tie;
						tfa[z] += bs.GetByteLength();
					}

					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("        %2d  %2d  %3s  %3s  %10.3f",
							parm->GetTableCount(),
							tia[z],
							parm->GetBW()?"yes":"no",
							parm->GetMTF()?"yes":"no",
							tfa[z]/1024.0/tciaa.size()
						);

					if (tfa[z] < btv) {
						btv = tfa[z];
						btr = z;
						if (m_IF.IsPL(BQB_PL_STANDARD))
							m_IF.printf("  <---");
					}
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("\n");

					if ((z > 0) && (tfa[z] > tflast))
						break;

					tflast = tfa[z];
				}
				parm->SetBlockLength(tia[btr]);
			}
		}

		if (level >= 3) {

			for (z=0;z<3;z++) {

				switch(z) {
					case 0: parm->SetBW(true);  parm->SetMTF(true); break;
					case 1: parm->SetBW(true);  parm->SetMTF(false); break;
					case 2: parm->SetBW(false); parm->SetMTF(false); break;
				}
				tfa[z] = 0;
				for (z2=0;z2<(int)tciaa.size();z2++) {

					bs.Clear();
					tie = new CBQBIntegerEngine(m_IF);

					tie->Compress(
						tciaa[z2],
						&bs,
						parm->GetCompressorParameterSet(),
						false,
						NULL
					);

					delete tie;
					tfa[z] += bs.GetByteLength();
				}

				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("        %2d  %2d  %3s  %3s  %10.3f\n",
						parm->GetTableCount(),
						parm->GetBlockLength(),
						parm->GetBW()?"yes":"no",
						parm->GetMTF()?"yes":"no",
						tfa[z]/1024.0/tciaa.size()
					);

			}
			fac = 1.0;
			if ((tfa[0] <= tfa[1]) && (tfa[0] <= tfa[2]/fac)) {
				parm->SetBW(true);
				parm->SetMTF(true);
			} else if ((tfa[1] <= tfa[0]) && (tfa[1] <= tfa[2]/fac)) {
				parm->SetBW(true);
				parm->SetMTF(false);
			} else {
				parm->SetBW(false);
				parm->SetMTF(false);
			}
		}
	}
_end:

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		BQBFormatTime((unsigned long)time(NULL)-t0,buf);
		m_IF.printf("        Optimization took %s\n",(const char*)buf);
	}

/*	if (m_fResiCorrFrame != NULL) {
		fclose(m_fResiCorrFrame);
		m_fResiCorrFrame = NULL;
	}

	if (m_fResiCorrTraj != NULL) {
		fclose(m_fResiCorrTraj);
		m_fResiCorrTraj = NULL;
	}*/

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.bprintf("    <<< Volumetric Trajectory Parameter Optimization Done <<<\n");

	return true;
}


bool CBQBDriver::DecompressCube(
		const char *inp,
		const char *outp,
		const char *readref,
		int steps,
		int stride,
		CBQBParameterSet_PosAndVol *parm
	) {

	UNUSED(stride);

	CBQBCubeFrame *cfr, *cfr2;
	int z, i, al, cl, ty, ver;
	int ahistused, chistused;
	CBQBBitSet bsat(m_IF), bscu(m_IF), bshe(m_IF);
	FILE *b, *fref;
	CBQBFile bf(m_IF);


	if (parm == NULL)
		parm = &m_ParmPosAndVol;

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Opening archive file \"%s\" ...\n",inp);
	if (!bf.OpenRead(inp)) {
		m_IF.eprintf("Error: Could not open file for reading.\n");
		return false;
	}

	if (outp != NULL) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("Opening output cube file \"%s\" ...\n",outp);
		b = fopen(outp,"wb");
		if (b == NULL) {
			m_IF.eprintf("Error: Could not open file for writing.\n");
			return false;
		}
	} else {
		b = NULL;
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("No output file specified; no output will be written.\n");
	}

	fref = NULL;
	if (readref != NULL) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("Opening reference cube file \"%s\" ...\n",readref);
		fref = fopen(readref,"rb");
		if (fref == NULL) {
			m_IF.eprintf("Error: Could not open file for reading.\n");
			if (outp != NULL)
				fclose(b);
			return false;
		}
	}

	m_pEngine->m_oaOutputCubeBuf.resize(7);
	for (z=0;z<7;z++)
		m_pEngine->m_oaOutputCubeBuf[z] = NULL;
	m_pEngine->m_iOutputCubeBufPos = 0;

	m_pEngine->m_oaOutputAtomBuf.resize(11);
	for (z=0;z<11;z++)
		m_pEngine->m_oaOutputAtomBuf[z] = NULL;
	m_pEngine->m_iOutputAtomBufPos = 0;

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Starting process...\n");

	i = 0;
	while (true) {

		if (m_IF.IsPL(BQB_PL_VERBOSE)) {
			m_IF.printf("\n");
			m_IF.printf("\n");
		}

		if (!bf.ReadFrame())
			break;

		ty = bf.GetFrameType();
		ver = bf.GetFrameTypeVersion();

		if ((ty == BQB_FRAMETYPE_IDX) || (ty == BQB_FRAMETYPE_COMPIDX))
			continue;

		if ((ty == BQB_FRAMETYPE_COMPCUBE) || (ty == BQB_FRAMETYPE_COMPCUBESTART) || (ty == BQB_FRAMETYPE_COMPCUBEKEY)) {
			if (ver > BQB_FRAMETYPE_COMPCUBE_VERSION) {
				m_IF.eprintf("Error: Unexpected frame type version (%d v%d).\n",ty,ver);
				m_IF.eprintf("       Probably the BQB file is more recent than the libbqb version.\n");
				return false;
			}
		} else {
			m_IF.eprintf("Error: Unexpected frame type (%d v%d).\n",ty,ver);
			return false;
		}

		bshe.Clear();
		bshe.m_iaData.assign(bf.GetFramePayload()->begin(),bf.GetFramePayload()->begin()+8);

		if (m_IF.IsPL(BQB_PL_STANDARD)) {
			if (i == 0)
				m_IF.printf("Reading first time step from compressed trajectory...\n");
			else
				m_IF.printf("    Step %d ...\n",i+1);
		}

		al = bshe.ReadBitsInteger(32);
		cl = bshe.ReadBitsInteger(32);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Atom block %d bytes, Cube block %d bytes.\n",al,cl);

		cfr = new CBQBCubeFrame(m_IF);
		m_pEngine->PushOutputCubeFrame(cfr);

		cfr->m_pAtoms = new CBQBAtomSet(m_IF);
		m_pEngine->PushOutputAtomFrame(cfr->m_pAtoms);

		bsat.Clear();
		bsat.m_iaData.assign(bf.GetFramePayload()->begin()+8,bf.GetFramePayload()->begin()+8+al);

		bscu.Clear();
		bscu.m_iaData.assign(bf.GetFramePayload()->begin()+8+al,bf.GetFramePayload()->begin()+8+al+cl);

		if (!m_pEngine->DecompressAtomFrame(
				&bsat,
				ver,
				parm->GetPositionParameterSet(),
				&ahistused,
				false,
				false
			))
			return false;

		if (!m_pEngine->DecompressCubeFrame(
				&bscu,
				ver,
				parm->GetVolumetricParameterSet(),
				&chistused,
				false,
				false
			))
			return false;

		if (m_IF.IsPL(BQB_PL_STANDARD) && (i == 0)) {
			m_IF.bprintf("\n");
			m_IF.bprintf("The following parameters have been used for compressing this file:\n");
			m_IF.bprintf("\n");
			m_IF.printf("%s\n",parm->ToString(4).c_str());
			m_IF.bprintf("Parameter Key:");
			m_IF.printf(" %s\n",parm->ToKey().c_str());
			m_IF.bprintf("\n");
		}

		if (readref != NULL) {

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("      Comparing to reference...\n");

			cfr2 = new CBQBCubeFrame(m_IF);

			if (!cfr2->ReadFrame(
					fref,
					cfr->m_iEps,
					cfr->m_iSigni,
					cfr->m_pAtoms->m_iSigni
				))
				break;

			for (z=0;z<cfr->m_iResXYZ;z++)
				if (!ExpMantisEqual(cfr->m_iaExpo[z],cfr->m_iaMantis[z],cfr2->m_iaExpo[z],cfr2->m_iaMantis[z])) {
					m_IF.eprintf("        Error %7d: %.10G vs %.10G\n",z,cfr->m_faBin[z],cfr2->m_faBin[z]);
					break;
				}

			delete cfr2;

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Done.\n");
		}

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("      %9.3f KiB compressed data unpacked (history: position %2d, volumetric %2d).\n",
				(al+cl)/1024.0, ahistused, chistused );

		if (outp != NULL)
			cfr->WriteFrame(b);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Done.\n");

		i++;
		if (i == steps)
			break;
	}

	bf.Close();

	if (outp != NULL)
		fclose(b);

	if (readref != NULL)
		fclose(fref);

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.printf("Finished decompressing the cube file.\n");
		m_IF.printf("\n");
	}

	return true;
}



bool CBQBDriver::CompressXYZ(
		const char *inp,
		const char *outp,
		const char *ref,
		int start,
		int steps,
		int stride,
		int keyfreq,
		CBQBParameterSet_Position *parm,
		bool comment,
		bool compare,
		int optimize,
		int optsteps,
		bool onlyopt
	) {

	const CBQBAtomSet *as;
	CBQBAtomSet *as2;
	int z, z2, i, k, o, ft, lasize, morder, atc, histused, mo;
	long insize;
	double tb, tb2, tf, tf2, tot;
	CBQBBitSet bsat(m_IF), bstemp(m_IF);
	bool err, to;
	CBQBFile bf(m_IF);
	CBQBStatistics stat(m_IF);


	atc = 0;

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Opening XYZ file \"%s\" ...\n",inp);

	m_pEngine->m_pReadCacheXYZ = new CBQBReadXYZCache(m_IF);
	m_pEngine->m_pReadCacheXYZ->SetReadParameters( parm->GetPrecision() );
	if (!m_pEngine->m_pReadCacheXYZ->FileOpenRead(inp,ref)) {
		m_IF.eprintf("Error: Could not open file for reading.\n");
		return false;
	}

	if (parm->GetUseExtra() && (parm->GetExtraTRange() > parm->GetOrder()))
		parm->SetOrder(parm->GetExtraTRange());

	if (optimize != 0) {

		if (parm->GetOrder()+1 > optsteps)
			m_pEngine->m_pReadCacheXYZ->SetHistoryDepth(parm->GetOrder()+1);
		else
			m_pEngine->m_pReadCacheXYZ->SetHistoryDepth(optsteps);

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("\n");
		m_pEngine->m_pReadCacheXYZ->CacheSteps(optsteps,true);
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("\n");

		if (m_pEngine->m_pReadCacheXYZ->GetCachedSteps() == 0) {
			m_IF.eprintf("Error: Zero frames read. Aborting.\n");
			return false;
		}

		if (m_iOptIncludeFirst == -1) {
			if ((m_pEngine->m_pReadCacheXYZ->GetCachedSteps() < optsteps) || (optsteps < 10)) {
				if (m_IF.IsPL(BQB_PL_STANDARD)) {
					m_IF.printf("        Low step count (%d), switching on -optstart.\n",m_pEngine->m_pReadCacheXYZ->GetCachedSteps());
					m_IF.printf("\n");
				}
				m_iOptIncludeFirst = 1;
			} else {
				if (m_IF.IsPL(BQB_PL_STANDARD)) {
					m_IF.printf("        Large step count (%d), switching off -optstart.\n",m_pEngine->m_pReadCacheXYZ->GetCachedSteps());
					m_IF.printf("\n");
				}
				m_iOptIncludeFirst = 0;
			}
		}

		m_pEngine->m_pReadCacheXYZ->SetStartLock();

		if (!OptimizeXYZParameters(
				optimize,       // int     level,
				parm,
				comment,        // bool    comment, 
				false           // bool    cubepart
		)) {
			m_IF.eprintf("Error: Parameter optimization failed.\n");
			delete m_pEngine->m_pReadCacheXYZ;
			m_pEngine->m_pReadCacheXYZ = NULL;
			return false;
		}

		if (m_IF.IsPL(BQB_PL_STANDARD)) {

			m_IF.printf("\n");
			m_IF.bprintf("    The optimal parameters for this trajectory are:\n");
			m_IF.printf("\n");

			m_IF.printf("%s",parm->ToString(4).c_str());
			m_IF.printf("\n");

			m_IF.bprintf("    Parameter key:");
			m_IF.printf("  %s\n",parm->ToKey().c_str());
			m_IF.printf("\n");
		}

		if (onlyopt) {
			m_pEngine->m_pReadCacheXYZ->CloseFile();
			delete m_pEngine->m_pReadCacheXYZ;
			m_pEngine->m_pReadCacheXYZ = NULL;
			return true;
		}

		m_pEngine->m_pReadCacheXYZ->RewindReadPos();

		m_pEngine->m_pReadCacheXYZ->LiftStartLock();

	} else {

		m_pEngine->m_pReadCacheXYZ->SetHistoryDepth(parm->GetOrder()+1);
	}

	stat.m_oStat.reset();

	if (outp != NULL) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("Opening compressed output file \"%s\" ...\n",outp);
		if (!bf.OpenWriteAppend(outp)) {
			m_IF.eprintf("Error: Could not open file for writing.\n");
			m_pEngine->m_pReadCacheXYZ->CloseFile();
			delete m_pEngine->m_pReadCacheXYZ;
			m_pEngine->m_pReadCacheXYZ = NULL;
			return false;
		}
	}

	m_pEngine->m_oaOutputAtomBuf.resize(parm->GetOrder()+1);
	for (z=0;z<parm->GetOrder()+1;z++)
		m_pEngine->m_oaOutputAtomBuf[z] = NULL;
	m_pEngine->m_iOutputAtomBufPos = 0;

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Starting process...\n");

	if (start != 0) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("Fast-forwarding to step %d ...\n",start+1);
		for (z=0;z<start;z++)
			if (!m_pEngine->m_pReadCacheXYZ->SkipOneStep()) {
				delete m_pEngine->m_pReadCacheXYZ;
				m_pEngine->m_pReadCacheXYZ = NULL;
				return false;
			}
	}

	if (parm->GetUseExtra()) {
		// Fix MB 11.04.2020
//		mo = parm->GetExtraTOrder();
		mo = parm->GetExtraTRange();
	} else
		mo = parm->GetOrder();

	morder = mo;

	i = 0;
	k = 0;
	tb = 0;
	tb2 = 0;
	tf = 0;
	tf2 = 0;
	to = false;

	if (parm->GetOptOrder())
		lasize = 1000000000;
	else
		lasize = -1;

	err = false;

	while (true) {

		if (m_IF.IsPL(BQB_PL_VERBOSE)) {
			m_IF.printf("\n");
			m_IF.printf("\n");
		}

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("    Step %d ...\n",i+1);

		as = m_pEngine->m_pReadCacheXYZ->GetNextFrame();

		if (i == 0)
			atc = (int)as->m_oaAtoms.size();

		if (as == NULL) {
			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("      Reached end of input trajectory.\n");
			break;
		}

		for (z=0;z<stride-1;z++) {
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Skipping frame...\n");
			if (!m_pEngine->m_pReadCacheXYZ->SkipOneStep()) {
				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("      Reached end of input trajectory.\n");
				break;
			}
		}

		if (i == 0) {
			if (parm->GetUseExtra()) {
				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("      Initializing Extrapolator (%d/%d)...\n",
						parm->GetExtraTRange(),
						parm->GetExtraTOrder()
					);
				m_pEngine->m_pExtrapolatorXYZ = new CBQBExtrapolator(m_IF);
				//m_IF.SetPrintLevel(BQB_PL_DEBUG);
				m_pEngine->m_pExtrapolatorXYZ->InitializeXYZ(
					parm->GetExtraTRange(),
					parm->GetExtraTOrder(),
					parm->GetExtraTimeExpo(),
					false
				);
				//m_IF.SetPrintLevel(BQB_PL_STANDARD);
			}
		}

		if (m_pEngine->m_pExtrapolatorXYZ != NULL)
			m_pEngine->m_pExtrapolatorXYZ->PushAtomFrame(as);

		if ((i != 0) && (keyfreq != 0) && ((i % keyfreq) == 0)) {
			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("      Enforcing key frame.\n");
			k = 0; // Enforce this frame to be a keyframe
		}

_again:
		if (i < mo)
			o = i;
		else if (k < mo)
			o = k;
		else
			o = mo;

_aagain:
		bsat.Clear();
		stat.PushStatistics();

		m_pEngine->CompressAtomFrame(
			&bsat,
			o,
			parm,
			histused,
			(i==0), // Store static trajectory information?
			comment,
			&stat
		);

		if (m_IF.IsPL(BQB_PL_STANDARD)) {
			if (parm->GetUseExtra())
				m_IF.printf("      History %2d, output size %9.3f KiB.\n",histused,bsat.GetByteLength()/1024.0);
			else
				m_IF.printf("      Order %d, output size %9.3f KiB.\n",o,bsat.GetByteLength()/1024.0);
		}

		if (lasize >= 0) { // Only happens if optorder == true

			if (!parm->GetUseExtra()) {

				if (to) {
					to = false;
					if (bsat.GetByteLength() <= bstemp.GetByteLength()) {
						if (m_IF.IsPL(BQB_PL_STANDARD))
							m_IF.printf("    Size is indeed smaller with lower order. Limiting atom order to %d.\n",o);
						stat.PopDiffStatistics();
						parm->SetOrder(o);
						mo = o;
						lasize = -1;
					} else {
						if (m_IF.IsPL(BQB_PL_STANDARD))
							m_IF.printf("    Size was smaller with higher order. Not limiting.\n");
						stat.PopStatistics();
						bsat = CBQBBitSet(bstemp);
						lasize = bsat.GetByteLength();
					}
				} else {
					if ((o != 0) && ((double)bsat.GetByteLength() >= (double)lasize*0.97)) {
						if (m_IF.IsPL(BQB_PL_STANDARD))
							m_IF.printf("    Size did not decrease any further with order %d. Trying order %d...\n",o,o-1);
						bstemp = CBQBBitSet(bsat);
						to = true;
						o--;
						goto _aagain;
					} else if (o == mo)
						lasize = -1;
					else
						lasize = bsat.GetByteLength();
				}
			}
		}

		stat.PopIgnoreStatistics();

		tb += bsat.GetByteLength();
		tf++;

		if ((bsat.GetLength()%8) != 0)
			stat.m_oStat.m_lOverhead += 8 - (bsat.GetLength()%8);

		if (i >= morder) {
			tb2 += bsat.GetByteLength();
			tf2++;
		}

		if (compare) {

			as2 = new CBQBAtomSet(m_IF);
			m_pEngine->PushOutputAtomFrame(as2);

			m_pEngine->DecompressAtomFrame(
				&bsat,
				BQB_FRAMETYPE_COMPTRAJ_VERSION,
				parm,
				NULL, // histused
				false,
				true
			);

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Comparing input and output...\n");

			for (z=0;z<(int)as2->m_oaAtoms.size();z++) {
				for (z2=0;z2<3;z2++) {
					if (as->m_oaAtoms[z]->m_iCoord[z2] != as2->m_oaAtoms[z]->m_iCoord[z2]) {
						m_IF.eprintf("        Error %d[%d]: %.6f != %.6f\n",
							z,
							z2,
							as->m_oaAtoms[z]->m_fCoord[z2],
							as2->m_oaAtoms[z]->m_fCoord[z2]
						);
						err = true;
						goto _skerr;
					}
				}
			}

_skerr:
			if (err) {
				if (k != 0) {
					m_IF.eprintf("Errors occurred. Compressing frame again with history zero.\n");
					err = false;
					k = 0;
					goto _again;
				} else {
					m_IF.eprintf("Errors occurred despite of zero history. Aborting.\n");
					goto _end;
				}
			}

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Done.\n");
		}

		if (outp != NULL) {
			if (i == 0)
				ft = BQB_FRAMETYPE_COMPTRAJSTART;
			else if (k == 0)
				ft = BQB_FRAMETYPE_COMPTRAJKEY;
			else
				ft = BQB_FRAMETYPE_COMPTRAJ;

			bf.CreateShortFrame(
				ft,
				BQB_FRAMETYPE_COMPTRAJ_VERSION,
				i+1
			);
			bf.PushPayload(bsat.m_iaData);
			bf.FinalizeFrame(&stat);
		}

		i++;
		k++;
		if (i == steps)
			break;
	}
_end:

	if (outp != NULL) {
		bf.WriteIndexFrame(true,&stat);
		bf.Close();
	}

	insize = m_pEngine->m_pReadCacheXYZ->ftell();

	m_pEngine->m_pReadCacheXYZ->CloseFile();

	delete m_pEngine->m_pReadCacheXYZ;
	m_pEngine->m_pReadCacheXYZ = NULL;

	if (m_pEngine->m_pExtrapolatorXYZ != NULL) {
		delete m_pEngine->m_pExtrapolatorXYZ;
		m_pEngine->m_pExtrapolatorXYZ = NULL;
	}

	if (err) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.eprintf("Errors occurred while compressing the XYZ file.\n");
		return false;
	} else if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.printf("Finished compressing the XYZ file.\n");
		m_IF.printf("\n");
		m_IF.bprintf("Size contributions:\n");
		tot = double(stat.m_oStat.m_lOverhead+stat.m_oStat.m_lAlphabet+stat.m_oStat.m_lHuffmanTables+stat.m_oStat.m_lTableSwitch+stat.m_oStat.m_lHuffmanData);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) overhead.\n",
			double(stat.m_oStat.m_lOverhead)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lOverhead>>3).string(),
			double(stat.m_oStat.m_lOverhead)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) alphabet data.\n",
			double(stat.m_oStat.m_lAlphabet)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lAlphabet>>3).string(),
			double(stat.m_oStat.m_lAlphabet)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) Huffman tables.\n",
			double(stat.m_oStat.m_lHuffmanTables)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lHuffmanTables>>3).string(),
			double(stat.m_oStat.m_lHuffmanTables)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) table switching.\n",
			double(stat.m_oStat.m_lTableSwitch)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lTableSwitch>>3).string(),
			double(stat.m_oStat.m_lTableSwitch)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) payload data.\n",
			double(stat.m_oStat.m_lHuffmanData)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lHuffmanData>>3).string(),
			double(stat.m_oStat.m_lHuffmanData)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) in total.\n",
			tot/1024.0/1024.0/8.0,
			((stat.m_oStat.m_lOverhead+stat.m_oStat.m_lAlphabet+stat.m_oStat.m_lHuffmanTables+stat.m_oStat.m_lTableSwitch+stat.m_oStat.m_lHuffmanData)>>3).string(),
			100.0
		);
		if (tf > 0) {
			m_IF.bprintf("\n");
			m_IF.bprintf("  * Totals:\n");
			m_IF.printf("      %9.3f KiB per frame on average.\n",tb/1024.0/tf);
			m_IF.printf("      %9.3f bits per coordinate on average.\n",tb/tf/atc/3.0*8.0);
			m_IF.printf("      Compression ratio of %.3f : 1\n",(double)insize/tb);
		}
		if (tf2 > 0) {
			m_IF.bprintf("\n");
			m_IF.bprintf("  * Starting from step %d:\n",morder+1);
			m_IF.printf("      %9.3f KiB per frame on average.\n",tb2/1024.0/tf2);
			m_IF.printf("      %9.3f bits per coordinate on average.\n",tb2/tf2/atc/3.0*8.0);
			m_IF.printf("      Compression ratio of %.3f : 1\n",((double)insize/tf*tf2)/tb2);
		}
		m_IF.printf("\n");
	}

	return true;
}


double CBQBDriver::CompressXYZBenchmark(
		CBQBParameterSet_Position *parm,
		bool comment, 
		bool realsize,
		double *presi,
		std::vector<std::vector<int> > *tciaa
	) {

	const CBQBAtomSet *as;
	int i, o, histused, mo;
	double tb, tf, resi, tresi;
	CBQBBitSet bsat(m_IF);


	resi = 0;

	if (parm->GetUseExtra() && (parm->GetExtraTRange() > parm->GetOrder()))
		parm->SetOrder(parm->GetExtraTRange());

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("Starting process...\n");

	i = 0;
	tb = 0;
	tf = 0;

	if (m_pEngine->m_pReadCacheCube != NULL)
		m_pEngine->m_pReadCacheCube->RewindReadPos();
	else
		m_pEngine->m_pReadCacheXYZ->RewindReadPos();

	if (parm->GetUseExtra())
		mo = parm->GetExtraTRange();
	else
		mo = parm->GetOrder();


	while (true) {

		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			m_IF.printf("\n");
			m_IF.printf("\n");
		}

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Step %d ...\n",i+1);

		if (m_pEngine->m_pReadCacheCube != NULL)
			as = m_pEngine->m_pReadCacheCube->GetNextFrame()->m_pAtoms;
		else
			as = m_pEngine->m_pReadCacheXYZ->GetNextFrame();

		if (as == NULL)
			break;

		if (i == 0) {

			if (!realsize)
				m_pEngine->BuildIdentityAtomSort(as);

			if (parm->GetUseExtra()) {
				if (m_IF.IsPL(BQB_PL_VERBOSE))
					m_IF.printf("      Initializing Extrapolator (%d|%d)...\n",
					parm->GetExtraTRange(),
					parm->GetExtraTOrder()
				);
				m_pEngine->m_pExtrapolatorXYZ = new CBQBExtrapolator(m_IF);
				m_pEngine->m_pExtrapolatorXYZ->InitializeXYZ(
					parm->GetExtraTRange(),
					parm->GetExtraTOrder(),
					parm->GetExtraTimeExpo(),
					!m_IF.IsPL(BQB_PL_VERBOSE)
				);
			}
		}

		if (m_pEngine->m_pExtrapolatorXYZ != NULL) 
			m_pEngine->m_pExtrapolatorXYZ->PushAtomFrame(as);

		if (i < mo)
			o = i;
		else
			o = mo;

		if ((i+1 >= mo) || (m_iOptIncludeFirst > 0)) {

			bsat.Clear();

			m_pEngine->CompressAtomFrame(
				&bsat,
				o,
				parm,
				histused,
				(i==0), // Store static trajectory information?
				comment,
				NULL,
				!realsize,
				&tresi,
				tciaa
			);

			resi += tresi;

			tf++;

			if (realsize) {
				tb += bsat.GetByteLength();
/*				if (m_fResiCorrFrame != NULL) {
					mfprintf(m_fResiCorrFrame,"%d;  %.3f\n", bsat.GetByteLength(), mypow( tresi / ((as->m_oaAtoms.size() + 1) * 3.0), 1.0/m_fEntropyExpoXYZ ) );
					fflush(m_fResiCorrFrame);
				}*/
			}
		}

		i++;

		if (m_pEngine->m_pReadCacheCube != NULL) {
			if (m_pEngine->m_pReadCacheCube->GetCachedSteps() == 0)
				break;
		} else {
			if (m_pEngine->m_pReadCacheXYZ->GetCachedSteps() == 0)
				break;
		}
	}

	if (m_pEngine->m_pExtrapolatorXYZ != NULL) {
		delete m_pEngine->m_pExtrapolatorXYZ;
		m_pEngine->m_pExtrapolatorXYZ = NULL;
	}

	resi = sqrt( resi / (tf * (as->m_oaAtoms.size() + 1) * 3.0) );

	if (presi != NULL)
		*presi = resi;

	if (realsize) {
/*		if (m_fResiCorrTraj != NULL) {
			mfprintf(m_fResiCorrTraj,"%.1f;  %.4f\n",tb / tf,resi);
			fflush(m_fResiCorrTraj);
		}*/
		return tb / tf;
	} else
		return resi;
}


double CBQBDriver::GoldenSectionSearch(double maxl, double eps, double (CBQBDriver::*fobj)(double)) {

	double a, c, fa, fc;


	a = 0;
	fa = (this->*fobj)(a);

	c = maxl;
	fc = (this->*fobj)(c);

	return REC_GoldenSectionSearch(a,c,c,fa,fc,fc,fobj,eps,0);
}


double CBQBDriver::REC_GoldenSectionSearch(double a, double b, double c, double fa, double fb, double fc, double (CBQBDriver::*fobj)(double), double eps, int depth) {

	double x, fx;


	if ((c - b) > (b - a))
		x = b + (2 - (1 + sqrt(5)) / 2) * (c - b);
	else
		x = b - (2 - (1 + sqrt(5)) / 2) * (b - a);

	if (fabs(c - a) < eps)
		return (c + a) / 2.0; 

	fx = (this->*fobj)(x);

	if (fx < fb) {
		if ((c - b) > (b - a))
			return REC_GoldenSectionSearch(b, x, c, fb, fx, fc, fobj, eps, depth+1);
		else
			return REC_GoldenSectionSearch(a, x, b, fa, fx, fb, fobj, eps, depth+1);
	} else {
		if ((c - b) > (b - a))
			return REC_GoldenSectionSearch(a, b, x, fa, fb, fx, fobj, eps, depth+1);
		else
			return REC_GoldenSectionSearch(x, b, c, fx, fb, fc, fobj, eps, depth+1);
	}
}


double CBQBDriver::OptimizeXYZ_ObjectiveTExpo(double texpo) {

	double entropy, resi;

	m_pPosParm->SetExtraTimeExpo(texpo);

	entropy = CompressXYZBenchmark(
		m_pPosParm,
		m_bOptXYZComment,   // bool comment, 
		m_bOptXYZRealSize,  // bool realsize
		&resi,               // double *presi
		NULL                 // std::vector<std::vector<int> > *tciaa
	);

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		if (m_bOptXYZRealSize)
			m_IF.printf("        %5d  %5d  %8.3f  %14.5f (%10.4f kiB)\n",
				m_pPosParm->GetExtraTRange(),
				m_pPosParm->GetExtraTOrder(),
				texpo,
				resi,
				entropy/1024.0
			);
		else
			m_IF.printf("        %5d  %5d  %8.3f  %14.5f\n",
				m_pPosParm->GetExtraTRange(),
				m_pPosParm->GetExtraTOrder(),
				texpo,
				entropy
			);
	}

	return entropy;
}


bool CBQBDriver::OptimizeXYZParameters(
		int     level,
		CBQBParameterSet_Position *parm,
		bool    comment, 
		bool    cubepart
	) {

	int z, z2, z3, smr, smo, steps, btr, btr2;
	int tia[5];
	double smv, tf, tflast, resi, btv, tfa[5], fac;
	bool realsize;
	unsigned long t0;
	char buf[256];
	std::vector<std::pair<double,std::string> > rlist;
	std::vector<std::vector<int> > tciaa;
	std::vector<int> iabl, iatc;
	CBQBBitSet bs(m_IF);
	CBQBIntegerEngine *tie;


	UNUSED(cubepart);
	tflast = 0;

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.bprintf("    >>> Position Trajectory Parameter Optimization >>>\n");

/*	if (level == 3) {
		m_fResiCorrFrame = OpenFileWrite("corrframexyz.csv",true);
		m_fResiCorrTraj = OpenFileWrite("corrtrajxyz.csv",true);
	}*/

	t0 = (unsigned long)time(NULL);

	if (level >= 3)
		realsize = true;
	else
		realsize = false;

	if (m_pEngine->m_pReadCacheCube != NULL) {
		m_pEngine->m_pReadCacheCube->RewindReadPos();
		steps = m_pEngine->m_pReadCacheCube->GetCachedSteps();
	} else {
		m_pEngine->m_pReadCacheXYZ->RewindReadPos();
		steps = m_pEngine->m_pReadCacheXYZ->GetCachedSteps();
	}

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.printf("        Optimization level %d\n",level);
		if ((level == 1) && !cubepart)
			m_IF.printf("          Please note: Level 1 and 2 are identical for position trajectories.\n");
		m_IF.printf("        Have %d cached steps.\n",steps);

		m_IF.printf("        Range  Order  Exponent        Residuum\n");
	}

	smv = 1.0e30;
	tf = 0;
	smr = -1;
	smo = -1;
	parm->SetUseExtra(true);
	//extratimeexpo = 1.0;

	for (z=1;z<MIN(steps+1,10);z++) { // trange

		for (z2=MAX(0,z-2);z2>=0;z2--) { // torder

	//		if ((z > 1) && (z2 == -1))
	//			continue;

			tflast = tf;

			parm->SetExtraTRange(z);
			parm->SetExtraTOrder(z2);

			tf = CompressXYZBenchmark(
				parm,
				comment,        // bool comment, 
				realsize,       // bool realsize
				&resi,          // double *presi
				NULL            // std::vector<std::vector<int> > *tciaa
			);

			if (m_IF.IsPL(BQB_PL_STANDARD)) {
				if (realsize)
					sprintf(buf,"        %5d  %5d  %8.3f  %14.5f (%10.4f kiB)",
						z,
						z2,
						parm->GetExtraTimeExpo(),
						resi,
						tf/1024.0
					);
				else
					sprintf(buf,"        %5d  %5d  %8.3f  %14.5f",
						z,
						z2,
						parm->GetExtraTimeExpo(),
						tf
					);

				m_IF.printf("%s",(const char*)buf);
			} else
				buf[0] = 0;

			rlist.push_back(std::pair<double,std::string>(tf,(const char*)buf));

			if (tf < smv) {
				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("  <---");
				smv = tf;
				smo = z2;
				smr = z;
			}

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("\n");

			if ((z2 < z-2) && (tflast < tf))
				break;
		}
	}

	std::sort(rlist.begin(),rlist.end());

/*	FILE *a;
	a = OpenFileWrite("optixyz_list.txt",true);
	mfprintf(a,"        Range  Order  Exponent       Entropy\n");
	for (z=0;z<(int)rlist.size();z++)
		mfprintf(a,"%s\n",rlist[z].second.c_str());
	fclose(a);*/

	parm->SetExtraTRange(smr);
	parm->SetExtraTOrder(smo);

	// Fix MB 11.04.2020
	if (parm->GetExtraTOrder()+2 < parm->GetExtraTRange()) {
//	if (parm->GetExtraTOrder()+2 != parm->GetExtraTRange()) {

		m_pPosParm = parm;

		m_bOptXYZRealSize = realsize;

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("\n");

		tf = GoldenSectionSearch(7.0, (level>=3)?0.02:0.1, &CBQBDriver::OptimizeXYZ_ObjectiveTExpo);

		if (level >= 3)
			parm->SetExtraTimeExpo( floor(tf*1000.0+0.5)/1000.0 );
		else
			parm->SetExtraTimeExpo( floor(tf*100.0+0.5)/100.0 );

	} else
		parm->SetExtraTimeExpo( 0 );

	if (level >= 2) {

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("        Optimizing Compressor settings:\n");

		tciaa.clear();

		tf = CompressXYZBenchmark(
			parm,
			comment,        // bool comment, 
			true,           // bool realsize
			NULL,           // double *presi
			&tciaa          // std::vector<std::vector<int> > *tciaa
		);

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("        TC   BL    BW   MTF   Frame Size [kiB]\n");

		if (level >= 4) {

			iatc.push_back(1);
			iatc.push_back(2);
			iatc.push_back(3);
			iatc.push_back(4);
			iatc.push_back(6);
			iatc.push_back(8);
		/*	iatc.push_back(10);
			iatc.push_back(12);
			iatc.push_back(15);
			iatc.push_back(18);
			iatc.push_back(21);
			iatc.push_back(24);*/

			iabl.push_back(10);
			iabl.push_back(12);
			iabl.push_back(14);
			iabl.push_back(16);
			iabl.push_back(20);
			iabl.push_back(24);
			iabl.push_back(28);
			iabl.push_back(32);
			iabl.push_back(36);
			iabl.push_back(40);
			iabl.push_back(48);
			iabl.push_back(56);
			iabl.push_back(64);

			btv = 1.0e30;
			btr = -1;
			btr2 = -1;

			for (z=0;z<(int)iatc.size();z++) {

				for (z2=0;z2<(int)iabl.size();z2++) {

					parm->SetBlockLength(iabl[z2]);
					parm->SetTableCount(iatc[z]);

					tf = 0;
					for (z3=0;z3<(int)tciaa.size();z3++) {

						bs.Clear();
						tie = new CBQBIntegerEngine(m_IF);

						tie->Compress(
							tciaa[z3],
							&bs,
							parm->GetCompressorParameterSet(),
							false,
							NULL
						);

						delete tie;
						tf += bs.GetByteLength();
					}

					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("        %2d   %2d   %3s   %3s   %10.4f",
							iatc[z],
							iabl[z2],
							parm->GetBW()?"yes":"no",
							parm->GetMTF()?"yes":"no",
							tf/tciaa.size()/1024.0
						);

					if (tf < btv) {
						btv = tf;
						if (m_IF.IsPL(BQB_PL_STANDARD))
							m_IF.printf("  <---");
						btr = z;
						btr2 = z2;
					}
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("\n");

					if (iatc[z] == 1)
						break; // One table: Block length does not matter
				}
			}

			parm->SetTableCount(iatc[btr]);
			parm->SetBlockLength(iabl[btr2]);

		} else {

			tia[0] = 1;
			tia[1] = 2;
			tia[2] = 3;
			tia[3] = 6;
			tia[4] = 12;

			btv = 1.0e30;
			btr = -1;
			for (z=0;z<5;z++) {

				tfa[z] = 0;
				for (z2=0;z2<(int)tciaa.size();z2++) {

					bs.Clear();
					tie = new CBQBIntegerEngine(m_IF);

					parm->SetTableCount(tia[z]);

					tie->Compress(
						tciaa[z2],
						&bs,
						parm->GetCompressorParameterSet(),
						false,
						NULL
					);

					delete tie;
					tfa[z] += bs.GetByteLength();
				}

				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("        %2d   %2d   %3s   %3s   %10.4f",
						tia[z],
						parm->GetBlockLength(),
						parm->GetBW()?"yes":"no",
						parm->GetMTF()?"yes":"no",
						tfa[z]/tciaa.size()/1024.0
					);

				if (tfa[z] < btv) {
					btv = tfa[z];
					btr = z;
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("  <---");
				}
				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("\n");

				if ((z > 0) && (tfa[z] > tflast))
					break;

				tflast = tfa[z];
			}
			parm->SetTableCount(tia[btr]);

			if (parm->GetTableCount() != 1) {

				tia[0] = 12;
				tia[1] = 16;
				tia[2] = 24;
				tia[3] = 32;
				tia[4] = 40;

				btv = 1.0e30;
				btr = -1;
				for (z=0;z<5;z++) {

					tfa[z] = 0;
					for (z2=0;z2<(int)tciaa.size();z2++) {

						bs.Clear();
						tie = new CBQBIntegerEngine(m_IF);

						parm->SetBlockLength(tia[z]);

						tie->Compress(
							tciaa[z2],
							&bs,
							parm->GetCompressorParameterSet(),
							false,
							NULL
						);

						delete tie;
						tfa[z] += bs.GetByteLength();
					}

					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("        %2d   %2d   %3s   %3s   %10.4f",
							parm->GetTableCount(),
							tia[z],
							parm->GetBW()?"yes":"no",
							parm->GetMTF()?"yes":"no",
							tfa[z]/tciaa.size()/1024.0
						);

					if (tfa[z] < btv) {
						btv = tfa[z];
						btr = z;
						if (m_IF.IsPL(BQB_PL_STANDARD))
							m_IF.printf("  <---");
					}
					if (m_IF.IsPL(BQB_PL_STANDARD))
						m_IF.printf("\n");

					if ((z > 0) && (tfa[z] > tflast))
						break;

					tflast = tfa[z];
				}
				parm->SetBlockLength(tia[btr]);
			}
		}

		if (level >= 3) {

			for (z=0;z<3;z++) {

				switch(z) {
					case 0: parm->SetBW(true); parm->SetMTF(true); break;
					case 1: parm->SetBW(true); parm->SetMTF(false); break;
					case 2: parm->SetBW(false); parm->SetMTF(false); break;
				}
				tfa[z] = 0;
				for (z2=0;z2<(int)tciaa.size();z2++) {

					bs.Clear();
					tie = new CBQBIntegerEngine(m_IF);

					tie->Compress(
						tciaa[z2],
						&bs,
						parm->GetCompressorParameterSet(),
						false,
						NULL
					);

					delete tie;
					tfa[z] += bs.GetByteLength();
				}

				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("        %2d   %2d   %3s   %3s   %10.4f\n",
						parm->GetTableCount(),
						parm->GetBlockLength(),
						parm->GetBW()?"yes":"no",
						parm->GetMTF()?"yes":"no",
						tfa[z]/tciaa.size()/1024.0
					);

			}
			fac = 1.0;
			if ((tfa[0] <= tfa[1]) && (tfa[0] <= tfa[2]/fac)) {
				parm->SetBW(true);
				parm->SetMTF(true);
			} else if ((tfa[1] <= tfa[0]) && (tfa[1] <= tfa[2]/fac)) {
				parm->SetBW(true);
				parm->SetMTF(false);
			} else {
				parm->SetBW(false);
				parm->SetMTF(false);
			}
		}
	}
	
	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		BQBFormatTime((unsigned long)time(NULL)-t0,buf);
		m_IF.printf("        Optimization took %s\n",(const char*)buf);
	}

/*	if (m_fResiCorrFrame != NULL) {
		fclose(m_fResiCorrFrame);
		m_fResiCorrFrame = NULL;
	}

	if (m_fResiCorrTraj != NULL) {
		fclose(m_fResiCorrTraj);
		m_fResiCorrTraj = NULL;
	}*/

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.bprintf("    <<< Position Trajectory Parameter Optimization Done <<<\n");

	return true;
}


bool CBQBDriver::DecompressXYZ(
		const char *inp,
		const char *outp,
		const char *readref,
		int steps,
		int stride,
		CBQBParameterSet_Position *parm
	) {

	UNUSED(stride);

	CBQBAtomSet *as, *as2;
	int z, z2, i, ty, ver, histused;
	CBQBBitSet bsat(m_IF);
	FILE *b, *fref;
	CBQBFile bf(m_IF);


	if (parm == NULL)
		parm = &m_ParmPos;

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Opening archive file \"%s\" ...\n",inp);
	if (!bf.OpenRead(inp)) {
		m_IF.eprintf("Error: Could not open file for reading.\n");
		return false;
	}

	b = NULL;
	if (outp != NULL) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("Opening output XYZ file \"%s\" ...\n",outp);
		b = fopen(outp,"wb");
		if (b == NULL) {
			m_IF.eprintf("Error: Could not open file for writing.\n");
			return false;
		}
	}

	fref = NULL;
	if (readref != NULL) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("Opening reference XYZ file \"%s\" ...\n",readref);
		fref = fopen(readref,"rb");
		if (fref == NULL) {
			m_IF.eprintf("Error: Could not open file for reading.\n");
			fclose(b);
			return false;
		}
	}

	m_pEngine->m_oaOutputAtomBuf.resize(11);
	for (z=0;z<11;z++)
		m_pEngine->m_oaOutputAtomBuf[z] = NULL;
	m_pEngine->m_iOutputAtomBufPos = 0;

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("Starting process...\n");

	i = 0;
	while (true) {

		if (m_IF.IsPL(BQB_PL_VERBOSE)) {
			m_IF.printf("\n");
			m_IF.printf("\n");
		}

		if (!bf.ReadFrame())
			break;

		ty = bf.GetFrameType();
		ver = bf.GetFrameTypeVersion();

		if ((ty == BQB_FRAMETYPE_IDX) || (ty == BQB_FRAMETYPE_COMPIDX))
			continue;

		if ((ty == BQB_FRAMETYPE_COMPTRAJ) || (ty == BQB_FRAMETYPE_COMPTRAJSTART) || (ty == BQB_FRAMETYPE_COMPTRAJKEY)) {
			if (ver > BQB_FRAMETYPE_COMPTRAJ_VERSION) {
				m_IF.eprintf("Error: Unexpected frame type version (%d v%d).\n",ty,ver);
				m_IF.eprintf("       Probably the BQB file is more recent than the libbqb version.\n");
				return false;
			}
		} else {
			m_IF.eprintf("Error: Unexpected frame type (%d v%d).\n",ty,ver);
			return false;
		}

		if (m_IF.IsPL(BQB_PL_STANDARD)) {
			if (i == 0)
				m_IF.printf("Reading first time step from compressed trajectory...\n");
			else
				m_IF.printf("    Step %d ...\n",i+1);
		}

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Atom block %lu bytes.\n",(unsigned long)bf.GetFramePayload()->size());

		as = new CBQBAtomSet(m_IF);
		m_pEngine->PushOutputAtomFrame(as);

		bsat.Clear();
		bsat.m_iaData.assign(bf.GetFramePayload()->begin(),bf.GetFramePayload()->end());

		m_pEngine->DecompressAtomFrame(
			&bsat,
			ver,
			parm,
			&histused,
			false,
			false
		);

		if (m_IF.IsPL(BQB_PL_STANDARD) && (i == 0)) {
			m_IF.bprintf("\n");
			m_IF.bprintf("The following parameters have been used for compressing this file:\n");
			m_IF.bprintf("\n");
			m_IF.printf("%s\n",parm->ToString(4).c_str());
			m_IF.bprintf("Parameter Key:");
			m_IF.printf(" %s\n",parm->ToKey().c_str());
			m_IF.bprintf("\n");
		}

		if (readref != NULL) {

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("      Comparing to reference...\n");

			as2 = new CBQBAtomSet(m_IF);

			if (!as2->ReadXYZ(fref,as->m_iSigni,NULL))
				break;

			for (z=0;z<(int)as2->m_oaAtoms.size();z++) {
				for (z2=0;z2<3;z2++)
					if (as->m_oaAtoms[z]->m_iCoord[z2] != as2->m_oaAtoms[z]->m_iCoord[z2])
						m_IF.eprintf("        Error %d[%d]: %.6f != %.6f\n",
							z,z2,as->m_oaAtoms[z]->m_fCoord[z2],as2->m_oaAtoms[z]->m_fCoord[z2]);
			}

			delete as2;

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("Done.\n");
		}

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("      %8.3f KiB compressed data unpacked (history %2d).\n",
				bf.GetFramePayload()->size()/1024.0,histused);

		if (b != NULL)
			as->WriteXYZ(b,as->m_iSigni);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("Done.\n");

		i++;
		if (i == steps)
			break;
	}

	bf.Close();

	if (b != NULL)
		fclose(b);

	if (readref != NULL)
		fclose(fref);

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.printf("Finished decompressing the XYZ file.\n");
		m_IF.printf("\n");
	}

	return true;
}


bool CBQBDriver::CompressFile(
		const char *inp,
		const char *outp,
		const CBQBParameterSet_Compressor *parm
	) {

	FILE *a;
	unsigned char buf[4096];
	double tot;
	std::vector<int> ia;
	CBQBBitSet bs(m_IF);
	CBQBIntegerEngine ie(m_IF);
	CBQBFile bf(m_IF);
	int i, z;
	CBQBStatistics stat(m_IF);


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CEngine::CompressFile >>>\n");

	stat.m_oStat.reset();

	a = fopen(inp,"rb");
	if (a == NULL) {
		m_IF.eprintf("Error: Could not open file \"%s\" for reading.\n",inp);
		return false;
	}

	if (outp != NULL) {
		if (!bf.OpenWriteReplace(outp)) {
			m_IF.eprintf("Error: Could not open file \"%s\" for writing.\n",outp);
			fclose(a);
			return false;
		}
	}

	while (!feof(a)) {
		i = (int)fread(buf,1,4096,a);
		for (z=0;z<i;z++)
			ia.push_back(buf[z]);
		if (i < 4096)
			break;
	}

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("%lu Bytes read from input file.\n",(unsigned long)ia.size());

	if (!ie.Compress(
			ia,
			&bs,
			parm,
			true,
			&stat
	)) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.eprintf("Compressor returned an error.\n");
		return false;
	}

	fclose(a);

	if (outp != NULL) {
		bf.CreateShortFrame(BQB_FRAMETYPE_COMPFILE,0,0);
		bf.PushPayload(bs.m_iaData);
		bf.FinalizeFrame(&stat);
		bf.Close();
	}

	if (m_IF.IsPL(BQB_PL_STANDARD)) {

		m_IF.printf("Compressed %lu to %d bytes (%.3f:1, %7.3f%%, %.3f bits/byte).\n",
			(unsigned long)ia.size(),
			bs.GetByteLength(),
			(double)ia.size()/bs.GetByteLength(),
			(double)bs.GetByteLength()/ia.size()*100.0,
			(double)bs.GetLength()/ia.size()
		);
		m_IF.printf("\n");

		tot = double(stat.m_oStat.m_lOverhead+stat.m_oStat.m_lAlphabet+stat.m_oStat.m_lHuffmanTables+stat.m_oStat.m_lTableSwitch+stat.m_oStat.m_lHuffmanData);

		m_IF.bprintf("Size contributions:\n");
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) overhead.\n",
			double(stat.m_oStat.m_lOverhead)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lOverhead>>3).string(),
			double(stat.m_oStat.m_lOverhead)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) alphabet data.\n",
			double(stat.m_oStat.m_lAlphabet)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lAlphabet>>3).string(),
			double(stat.m_oStat.m_lAlphabet)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) Huffman tables.\n",
			double(stat.m_oStat.m_lHuffmanTables)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lHuffmanTables>>3).string(),
			double(stat.m_oStat.m_lHuffmanTables)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) table switching.\n",
			double(stat.m_oStat.m_lTableSwitch)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lTableSwitch>>3).string(),
			double(stat.m_oStat.m_lTableSwitch)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) payload data.\n",
			double(stat.m_oStat.m_lHuffmanData)/1024.0/1024.0/8.0,
			(stat.m_oStat.m_lHuffmanData>>3).string(),
			double(stat.m_oStat.m_lHuffmanData)/tot*100.0
		);
		m_IF.printf("%10.3f MiB (%12s Bytes, %7.3f %%) in total.\n",
			tot/1024.0/1024.0/8.0,
			((stat.m_oStat.m_lOverhead+stat.m_oStat.m_lAlphabet+stat.m_oStat.m_lHuffmanTables+stat.m_oStat.m_lTableSwitch+stat.m_oStat.m_lHuffmanData)>>3).string(),
			100.0
		);
		m_IF.printf("\n");
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("<<< CEngine::CompressFile <<<\n");

	return true;
}


bool CBQBDriver::DecompressFile(const char *inp, const char *outp, const char *ref) {

	UNUSED(ref);

	FILE *b;
	unsigned char buf[4096];
	std::vector<int> ia;
	CBQBBitSet bs(m_IF);
	CBQBIntegerEngine ie(m_IF);
	int z, z2;
	CBQBFile bf(m_IF);
	CBQBParameterSet_Compressor parm(m_IF);


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CEngine::DecompressFile >>>\n");

	if (!bf.OpenRead(inp)) {
		m_IF.eprintf("Error: Could not open file \"%s\" for reading.\n",inp);
		return false;
	}

	b = fopen(outp,"wb");
	if (b == NULL) {
		m_IF.eprintf("Error: Could not open file \"%s\" for writing.\n",outp);
		return false;
	}

	bf.ReadFrame();

	if ((bf.GetFrameType() != BQB_FRAMETYPE_COMPFILE) || (bf.GetFrameTypeVersion() != 0)) {
		m_IF.eprintf("Error: Unexpected frame type (expected %d v0, found %d v%d).\n",
			BQB_FRAMETYPE_COMPFILE,bf.GetFrameType(),bf.GetFrameTypeVersion());
		fclose(b);
		return false;
	}

	bs.m_iaData.assign(bf.GetFramePayload()->begin(),bf.GetFramePayload()->end());

	bf.Close();

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("%lu bytes of payload read from input file.\n",(unsigned long)bs.m_iaData.size());

	if (!ie.Decompress(
			&bs,
			ia,
			&parm
		)) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.eprintf("IntegerEngine returned an error.\n");
		fclose(b);
		return false;
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("Decompressed to %lu Bytes.\n",(unsigned long)ia.size());

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.bprintf("\n");
		m_IF.bprintf("The following parameters have been used for compressing this file:\n");
		m_IF.bprintf("\n");
		m_IF.printf("%s\n",parm.ToString(4).c_str());
		m_IF.bprintf("Parameter Key:");
		m_IF.printf(" %s\n",parm.ToKey().c_str());
		m_IF.bprintf("\n");
	}

	for (z=0;z<(int)ia.size()/4096;z++) {
		for (z2=0;z2<4096;z2++)
			buf[z2] = (unsigned char)ia[z*4096+z2];
		(void)!fwrite(buf,4096,1,b);
	}

	if ((ia.size()%4096) != 0) {
		for (z2=0;z2<(int)ia.size()%4096;z2++)
			buf[z2] = (unsigned char)ia[z*4096+z2];
		(void)!fwrite(buf,ia.size()%4096,1,b);
	}

	fclose(b);

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("<<< CEngine::DecompressFile <<<\n");

	return true;
}



bool CBQBDriver::MergeBQB(const char *outfile, const std::vector<const char*> &infile) {

	int z, i, zi, zo;
	CBQBFile *pout;
	std::vector<CBQBFile*> pin;
	CBQBStatistics stat(m_IF);


	if (m_IF.IsPL(BQB_PL_STANDARD)) {

		m_IF.printf("\n");
		m_IF.printf("    *****************************\n");
		m_IF.printf("    ***   Merging BQB Files   ***\n");
		m_IF.printf("    *****************************\n");
		m_IF.printf("\n");

		m_IF.printf("    Will merge the input files\n");
		for (z=0;z<(int)infile.size();z++)
			m_IF.printf("      [%d] %s\n",z+1,infile[z]);
		m_IF.printf("    into the output file %s.\n",outfile);
		m_IF.printf("\n");
	}

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("    Opening input files...\n");
	pin.resize(infile.size());
	for (z=0;z<(int)infile.size();z++) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("      [%d] %s... ",z+1,infile[z]);
		pin[z] = new CBQBFile(m_IF);
		if (!pin[z]->OpenRead(infile[z])) {
			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.eprintf("Error.\n");
			return false;
		}
		if (m_IF.IsPL(BQB_PL_STANDARD)) {
			if ((i = pin[z]->GetTotalFrameCount()) != -1)
				m_IF.printf("Ok, %d frames (index).",i);
			else
				m_IF.printf("Ok, no index.");
			m_IF.printf("\n");
		}
	}

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.printf("\n");
		m_IF.printf("    Opening output file...\n");
	}
	pout = new CBQBFile(m_IF);
	if (!pout->OpenWriteReplace(outfile)) {
		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.eprintf("Error.\n");
		return false;
	}

	zo = 0;
	for (z=0;z<(int)pin.size();z++) {

		if (m_IF.IsPL(BQB_PL_STANDARD))
			m_IF.printf("    Processing input file %d.\n",z+1);
		zi = 0;

		while (true) {
			if (!pin[z]->ReadFrame()) {
				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("    Reached end of input file %d after reading %d frames.\n",
						z+1,zi);
				break;
			}
			zi++;

			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf("      Input file %d frame %6d: Type %2d.%d, ID %6d, Size %7.2f KiB",
					z+1,
					zi,
					pin[z]->GetFrameType(),
					pin[z]->GetFrameTypeVersion(),
					pin[z]->GetFrameID(),
					pin[z]->GetFramePayload()->size()/1024.0
				);

			if ((pin[z]->GetFrameType() == 2) || (pin[z]->GetFrameType() == 3)) {
				if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf(" --> Skipping index frame.\n");
				continue;
			}
			if (m_IF.IsPL(BQB_PL_STANDARD))
				m_IF.printf(" --> Output frame %6d.\n",zo);
			pout->CreateShortFrame(
				pin[z]->GetFrameType(),
				pin[z]->GetFrameTypeVersion(),
				pin[z]->GetFrameID()
			);
			pout->PushPayload(*pin[z]->GetFramePayload());
			pout->FinalizeFrame(&stat);
			zo++;
		}
	}

	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.printf("    Merged %d frames from %lu input files into output file.\n",
			zo,(unsigned long)infile.size());
		m_IF.printf("\n");
	}

	if (m_IF.IsPL(BQB_PL_STANDARD))
		m_IF.printf("    Writing new index frame for output file...\n");
	pout->WriteIndexFrame(true,&stat);
	if (m_IF.IsPL(BQB_PL_STANDARD)) {
		m_IF.printf("    All done.\n");
		m_IF.printf("\n");
	}
	pout->Close();
	delete pout;
	for (z=0;z<(int)pin.size();z++) {
		pin[z]->Close();
		delete pin[z];
	}

	return true;
}


/*bool CBQBEngine::SplitBQB(const char *infile, const char *outbase, int splitlength, int steps, bool verbose) {

	int z, z2, i, k, ofc, ofc2, al, cl;
	int ao, co, aorder, corder;
	int lcsize, lasize, ty, ver;
	bool err, to;
	CBQBFile bqin, *bqout;
	CBQBCubeFrame *cfr, *cfr2;
	CBQBBitSet bshe, bsat, bscu, bstemp;
	CxString buf;


	lcsize = 1000000000; // To avoid "use of uninitialized variable" warning
	lasize = 1000000000;
	ofc = 0;
	ofc2 = 0;

	m_IF.printf("\n");
	m_IF.printf(WHITE,"    ******************************\n");
	m_IF.printf(WHITE,"    ***   Splitting BQB File   ***\n");
	m_IF.printf(WHITE,"    ******************************\n\n");

	SetBarbecubeVerbose(verbose);

	if (!bqin.OpenRead(infile)) {
		m_IF.eprintf("CBQBEngine::SplitBQB(): Error: Could not open BQB file \"%s\" for reading.\n",infile);
		return false;
	}

	m_oaOutputCubeBuf.resize(9);
	m_oaOutput2CubeBuf.resize(9);
	m_oaInputCubeBuf.resize(9);
	m_oaInputAtomBuf.resize(9);
	m_oaOutputAtomBuf.resize(9);
	m_oaOutput2AtomBuf.resize(9);
	for (z=0;z<9;z++) {
		m_oaInputCubeBuf[z] = NULL;
		m_oaOutputCubeBuf[z] = NULL;
		m_oaOutput2CubeBuf[z] = NULL;
		m_oaInputAtomBuf[z] = NULL;
		m_oaOutputAtomBuf[z] = NULL;
		m_oaOutput2AtomBuf[z] = NULL;
	}
	m_iInputCubeBufPos = 0;
	m_iOutputCubeBufPos = 0;
	m_iOutput2CubeBufPos = 0;
	m_iInputAtomBufPos = 0;
	m_iOutputAtomBufPos = 0;
	m_iOutput2AtomBufPos = 0;

	i = 0;
	k = 0;
	aorder = 8;
	corder = 8;
	bqout = NULL;
	to = false;

	m_IF.printf("\n");
	while (true) {

		m_IF.printf("    Reading input frame %d...\n",i+1);
		if (!bqin.ReadFrame()) {
			m_IF.printf("Could not read further frame from input file.\n");
			break;
		}

		ty = bqin.GetFrameType();
		ver = bqin.GetFrameTypeVersion();

		if ((ty == 2) || (ty == 3)) {
			m_IF.printf("      Skipping index frame.\n");
			goto _skip;
		}
		if ((ty != 8) && (ty != 9)) {
			m_IF.eprintf("CBQBEngine::SplitBQB(): Error: Splitting only implemented for compressed cube frames yet.\n");
			return false;
		}

		if (bqout == NULL) {
			buf.sprintf("%s%03d.bqb",outbase,k+1);
			m_IF.printf("Creating next output file \"%s\".\n",(const char*)buf);
			bqout = new CBQBFile();
			if (!bqout->OpenWriteReplace((const char*)buf)) {
				m_IF.eprintf("CBQBEngine::SplitBQB(): Error: Could not open BQB file \"%s\" for writing.\n",(const char*)buf);
				return false;
			}
			lcsize = 1000000000;
			lasize = 1000000000;
			aorder = 8;
			corder = 8;
			ofc = 0;
			ofc2 = 0;
		}

		m_IF.printf("      Decompressing...\n");

		bshe.Clear();
		bshe.m_iaData.assign(bqin.GetFramePayload()->begin(),bqin.GetFramePayload()->begin()+8);

		al = bshe.ReadBitsInteger(32);
		cl = bshe.ReadBitsInteger(32);

		cfr = new CBQBCubeFrame();
		PushOutputCubeFrame(cfr);
		PushInputCubeFrame_NoDelete(cfr);

		cfr->m_pAtoms = new CBQBAtomSet();
		PushOutputAtomFrame(cfr->m_pAtoms);
		PushInputAtomFrame_NoDelete(cfr->m_pAtoms);

		bsat.Clear();
		bsat.m_iaData.assign(bqin.GetFramePayload()->begin()+8,bqin.GetFramePayload()->begin()+8+al);

		bscu.Clear();
		bscu.m_iaData.assign(bqin.GetFramePayload()->begin()+8+al,bqin.GetFramePayload()->begin()+8+al+cl);

		if (!DecompressAtomFrame(&bsat,ver,verbose,false))
			return false;

		if (!DecompressCubeFrame(&bscu,ver,verbose,false))
			return false;

		if (ofc >= 9) {

			m_IF.printf("      Pass-through writing output file %d, frame %d...\n",k+1,ofc+1);

			bqout->CreateShortFrame((ofc==0)?BQB_FRAMETYPE_COMPCUBESTART:BQB_FRAMETYPE_COMPCUBE,0,i+1);
			bqout->PushPayload(*bqin.GetFramePayload());
			bqout->FinalizeFrame();

		} else {

			m_IF.printf("      Re-compressing...\n");

_again:
			if (ofc2 < corder)
				co = ofc2;
			else if (ofc < corder)
				co = ofc;
			else
				co = corder;

			if (ofc2 < aorder)
				ao = ofc2;
			else if (ofc < aorder)
				ao = ofc;
			else
				ao = aorder;

_aagain:
			bsat.Clear();
			CompressAtomFrame(
				&bsat,
				ao,
				6,         // aprecision
				31,        // asplit
				40,        // ablock
				4,         // atables
				false,     // optatables
				true,      // acoderun
				true,      // abw
				true,      // amtf
				(ofc2>=2), // apreopt
				10,        // amaxiter
				true,      // asortatom
				(ofc==0),  // atominfo
				true,      // keepcomment
				0,         // amaxchunk
				verbose
				);

			m_IF.printf("        Atoms: Order %d, output size %9.3f KiB.\n",ao,bsat.GetByteLength()/1024.0);

			if (lasize >= 0) {
				if (to) {
					to = false;
					if (bsat.GetByteLength() <= bstemp.GetByteLength()) {
						m_IF.printf("        Size is indeed smaller with lower order. Limiting atom order to %d.\n",ao);
						aorder = ao;
						lasize = -1;
					} else {
						m_IF.printf("        Size was smaller with higher order. Not limiting.\n");
						bsat = CBQBBitSet(bstemp);
						lasize = bsat.GetByteLength();
					}
				} else {
					if ((ao != 0) && ((double)bsat.GetByteLength() >= (double)lasize*0.97)) {
						m_IF.printf("        Size did not decrease any further with order %d. Trying order %d...\n",ao,ao-1);
						bstemp = CBQBBitSet(bsat);
						to = true;
						ao--;
						goto _aagain;
					} else if (ao == aorder)
						lasize = -1;
					else
						lasize = bsat.GetByteLength();
				}
			}

_cagain:
			bscu.Clear();
			CompressCubeFrame(
				&bscu,
				co,
				12,       // ceps
				5,        // csigni
				31,       // csplit
				40,       // cblock
				6,        // ctables
				false,    // copttables
				true,     // chilbert
				1.075,    // nbhfac
				true,     // ccoderun
				false,    // cbw
				false,    // cmtf
				false,    // cpreopt
				10,       // cmaxiter
				(ofc==0), // atominfo
				0,        // cmaxchunk
				verbose
				);

			m_IF.printf("        Cube:  Order %d, output size %9.3f KiB.\n",co,bscu.GetByteLength()/1024.0);

			if (lcsize >= 0) {
				if (to) {
					to = false;
					if (bscu.GetByteLength() <= bstemp.GetByteLength()) {
						m_IF.printf("        Size is indeed smaller with lower order. Limiting cube order to %d.\n",co);
						corder = co;
						lcsize = -1;
					} else {
						m_IF.printf("        Size was smaller with higher order. Not limiting.\n");
						bscu = CBQBBitSet(bstemp);
						lcsize = bscu.GetByteLength();
					}
				} else {
					if ((co != 0) && ((double)bscu.GetByteLength() >= (double)lcsize*0.97)) {
						m_IF.printf("        Size did not decrease any further with order %d. Trying order %d...\n",co,co-1);
						bstemp = CBQBBitSet(bscu);
						to = true;
						co--;
						goto _cagain;
					} else if (co == corder)
						lcsize = -1;
					else
						lcsize = bscu.GetByteLength();
				}
			}

			m_IF.printf("        Comparing input and output...\n");

			cfr2 = new CBQBCubeFrame();
			PushOutput2CubeFrame(cfr2);
			cfr2->m_pAtoms = new CBQBAtomSet();
			PushOutput2AtomFrame(cfr2->m_pAtoms);

			DecompressAtomFrame(&bsat,1,verbose,true);

			DecompressCubeFrame(&bscu,1,verbose,true);

			err = false;
			for (z=0;z<(int)cfr2->m_pAtoms->m_oaAtoms.size();z++)
				for (z2=0;z2<3;z2++)
					if (cfr->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2] != cfr2->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2]) {
						m_IF.eprintf("        Error in atom coordinate %d[%d]: %.6f (%ld) != %.6f (%ld)\n",z,z2,cfr->m_pAtoms->m_oaAtoms[z]->m_fCoord[z2],cfr->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2],cfr2->m_pAtoms->m_oaAtoms[z]->m_fCoord[z2],cfr2->m_pAtoms->m_oaAtoms[z]->m_iCoord[z2]);
						err = true;
					}

			for (z=0;z<cfr->m_iResXYZ;z++)
				if (!ExpMantisEqual(cfr->m_iaExpo[z],cfr->m_iaMantis[z],cfr2->m_iaExpo[z],cfr2->m_iaMantis[z])) {
					m_IF.eprintf("        Error in volumetric data element %7d: %.10G vs %.10G\n",z,cfr->m_faBin[z],cfr2->m_faBin[z]);
					err = true;
				}

			if (err) {
				if (ofc2 != 0) {
					m_IF.eprintf("Errors occurred. Compressing frame again with order zero.\n");
					err = false;
					ofc2 = 0;
					goto _again;
				} else {
					m_IF.eprintf("Errors occurred despite of order zero. Aborting.\n");
					return false;
				}
			}

			bshe.Clear();
			bshe.WriteBits(bsat.GetByteLength(),32);
			bshe.WriteBits(bscu.GetByteLength(),32);

			m_IF.printf("      Writing output file %d, frame %d...\n",k+1,ofc+1);

			bqout->CreateShortFrame((ofc2==0)?BQB_FRAMETYPE_COMPCUBESTART:BQB_FRAMETYPE_COMPCUBE,0,i+1);
			bqout->PushPayload(bshe.m_iaData);
			bqout->PushPayload(bsat.m_iaData);
			bqout->PushPayload(bscu.m_iaData);
			bqout->FinalizeFrame();
		}

_skip:
		i++;
		ofc++;
		ofc2++;

		if ((i % splitlength) == 0) {
			m_IF.printf("Closing output file \"%s\".\n",(const char*)buf);
			bqout->WriteIndexFrame(true);
			bqout->Close();
			delete bqout;
			bqout = NULL;
			k++;
		}

		if (i == steps) {
			m_IF.printf("Reached step limit.\n\n");
			break;
		}
	}

	if (bqout != NULL) {
		m_IF.printf("Closing output file \"%s\".\n\n",(const char*)buf);
		bqout->WriteIndexFrame(true);
		bqout->Close();
		delete bqout;
		bqout = NULL;
		k++;
	}

	m_IF.printf("Processed %d frames, wrote %d output files.\n\n",i,k);

	bqin.Close();

	// Those are duplicates of the output buffers and deleted there
	m_oaInputAtomBuf.clear();
	m_oaInputCubeBuf.clear();

	m_IF.printf("All done. Leaving.\n");

	return true;
}*/


