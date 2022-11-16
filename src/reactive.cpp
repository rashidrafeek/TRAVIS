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

#include "reactive.h"
#include "tools.h"
#include "globalvar.h"
#include "maintools.h"
#include "timestep.h"
#include "posdomain.h"


const char *GetRevisionInfo_reactive(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_reactive() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



CReactiveRingBuffer::CReactiveRingBuffer() {

	m_iCurrentPos = -1;
	m_iAtomCount = 0;
	m_iDepth = 0;
}


CReactiveRingBuffer::~CReactiveRingBuffer() {

}


void CReactiveRingBuffer::Init(int depth, int atomcount) {

	int z;

	m_iAtomCount = atomcount;
	m_iDepth = depth;

	m_vaaCoords.resize(atomcount);

	for (z=0;z<atomcount;z++)
		m_vaaCoords[z].resize(depth);

	if (g_bNPT) {
		if (g_bBoxNonOrtho)
			m_maCellMatrix.resize(depth);
		else
			m_vaCellVector.resize(depth);
	}

	m_iCurrentPos = 0;
}


void CReactiveRingBuffer::PushFrame(const CxDVec3Array &va) {

	int z;

	m_iCurrentPos++;
	if (m_iCurrentPos >= m_iDepth)
		m_iCurrentPos = 0;

	for (z=0;z<m_iAtomCount;z++)
		m_vaaCoords[z][m_iCurrentPos] = va[z];
}


void CReactiveRingBuffer::PushCellMatrix(const CxDMatrix3 &m) {

	m_maCellMatrix[m_iCurrentPos] = m;
}


void CReactiveRingBuffer::PushCellVector(double x, double y, double z) {

	m_vaCellVector[m_iCurrentPos] = CxDVector3(x,y,z);
}



/*****************************************************************************************************************/



CReactiveEngine::CReactiveEngine() {

	m_bFirst = true;
	m_iBondPersist = 0;
}


CReactiveEngine::~CReactiveEngine() {

}


void CReactiveEngine::Parse() {

	CTimeStep ts;
	CPosDomainEngine de;
	CReactiveRingBuffer rb;
	int step, z;
	double maxbond;


	mprintf("\n");
	mprintf(WHITE,"    ######################################\n");
	mprintf(WHITE,"    ####    Reactive Molecule Scan    ####\n");
	mprintf(WHITE,"    ######################################\n");
	mprintf("\n");

	if (m_bFirst) {

		m_bFirst = false;

		mprintf("    You specified that bonds might form or break inside your trajectory.\n");
		mprintf("    TRAVIS will now scan the trajectory file to detect these events, and list all occurring species.\n\n");

		mprintf("    To avoid fast fluctuations of bonds, a bond formation / breaking is only adopted if it persists\n");
		mprintf("    for a given number of trajectory frames without interruption.\n\n");

		m_iBondPersist = AskUnsignedInteger("    Enter number of frames that bond changes need to persist: [200] ",200);

		mprintf("\n");
	}

	maxbond = 0;
	for (z=0;z<g_oaAtoms.GetSize()-1;z++)
		if (((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius > maxbond)
			maxbond = ((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius;
	maxbond *= 2 * g_fBondFactor;


	if (!OpenInputTrajectory()) {
		eprintf("Error: Could not open position trajectory for reading.\n");
		abort();
	}
	
	if (g_bNPT && (g_sNPTFile[0] != 0)) {
		g_fNPTFile = fopen((const  char*)g_sNPTFile,"rt");
		if (g_fNPTFile == NULL) {
			eprintf("Error: Could not open cell vector file \"%s\" for reading.\n",(const  char*)g_sNPTFile);
			abort();
		}
	}

	rb.Init(2*m_iBondPersist,g_iGesAtomCount);

	step = 0;
	while (true) {

		if (!ts.ReadTimestep(g_fPos,true)) {

			if (InputTrajectoryEOF())
				mprintf("End of trajectory file reached.\n");
			else
				eprintf("Error while reading trajectory frame.\n");
			break;
		}

		rb.PushFrame(ts.m_vaCoords);

		if (g_bNPT) {

			if (g_sNPTFile[0] != 0)
				ts.ReadCellVector(g_fNPTFile);

			if (g_bBoxNonOrtho)
				rb.PushCellMatrix(g_mBoxFromOrtho);
			else
				rb.PushCellVector(g_fBoxX,g_fBoxY,g_fBoxZ);
		}

		if (g_bPeriodicX && g_bPeriodicY && g_bPeriodicZ)
			de.Create( ts.m_vaCoords, g_mBoxFromOrtho, maxbond );
		else
			de.CreateTrivial( ts.m_vaCoords );

		step++;
	}

	if ((g_bNPT) && (g_sNPTFile[0] != 0))
		fclose(g_fNPTFile);

	CloseInputTrajectory();

	mprintf("\n");
	mprintf(WHITE,"    ##############################################\n");
	mprintf(WHITE,"    ####    Reactive Molecule Scan Finished   ####\n");
	mprintf(WHITE,"    ##############################################\n");
	mprintf("\n");
}
