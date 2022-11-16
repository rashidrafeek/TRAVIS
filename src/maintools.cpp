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


// This must always be the first include directive
#include "config.h"

#include "maintools.h"
#include "conversion.h"
#include "globalvar.h"

#include "interface.h"
#include "bicgstab.h"

#ifdef TARGET_LINUX
#include <unistd.h>
#endif

#ifdef TARGET_WINDOWS
#include <windows.h>
#endif

//>>>OMP
#ifdef USE_OMP
#include <omp.h>
#endif
//<<<OMP



const char *GetRevisionInfo_maintools(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_maintools() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



bool OpenInputTrajectory() {

	if (g_bVerbose)
		mprintf(">>> OpenInputTrajectory() >>>\n");

	if (g_iTrajFormat == 7) {

		if (g_pBQBFile == NULL) {

			g_pBQBFile = new CBQBFile(*g_pBQBInterface);
			//g_pBQBEngine = new CBQBEngine();
			//g_pBQBEngine->PrepareDecompressCube(false);
			if (g_bVerbose) {
				mprintf("\nSwitching on BQB verbosity.\n\n");
				g_pBQBInterface->SetPrintLevel(BQB_PL_DEBUG);
			}
			//SetBQBVerbose(g_bVerbose);
			if (!g_pBQBFile->OpenRead(g_sInputTraj)) {
				eprintf("Error: Could not open BQB file \"%s\".\n",g_sInputTraj);
				return false;
			}

		} else {

			if (!g_pBQBFile->Rewind())
				return false;

		}

	} else {

		g_fPos = fopen(g_sInputTraj,"rb"); 
		if (g_fPos == NULL) {
			eprintf("\nError: Could not open position trajectory \"%s\" for reading.\n",g_sInputTraj);
			return false;
		}
	}

	if (g_bVerbose)
		mprintf("<<< OpenInputTrajectory() <<<\n");

	return true;
}


bool CloseInputTrajectory() {

	if (g_bVerbose)
		mprintf(">>> CloseInputTrajectory() >>>\n");

	if (g_iTrajFormat == 7) {

	} else {

		if (g_fPos == NULL) {
			eprintf("\nCloseInputTrajectory(): Error: File pointer is NULL.\n");
			return false;
		}

		fclose(g_fPos);
		g_fPos = NULL;
	}

	if (g_bVerbose)
		mprintf("<<< CloseInputTrajectory() <<<\n");

	return true;
}


bool InputTrajectoryEOF() {

#ifdef USE_NETCDF
	if (g_iTrajFormat == 6 && (int)g_iSteps >= g_iTrajSteps)
		return true;
#endif
	if (g_iTrajFormat == 7)
		return g_pBQBFile->IsEOF();

	if (g_fPos == NULL) {
		eprintf("\nInputTrajectoryEOF(): Error: File pointer is NULL.\n");
		return true;
	}

	return (feof(g_fPos) != 0);
}


bool SeekInputTrajectory(int step) {

	if (g_bVerbose)
		mprintf(">>> SeekInputTrajectory() >>>\n");

	if (g_iTrajFormat == 7) {

		eprintf("Error: SeekInputTrajectory() not yet implemented for BQB files.\n");
		abort();

	} else {

		if (g_fPos == NULL) {
			eprintf("\nSeekInputTrajectory(): Error: File pointer is NULL.\n");
			return false;
		}

		if (step == g_iBeginStep) {
			fseek(g_fPos,g_iFastForwardPos,SEEK_SET);
		} else {
			eprintf("\nSeekInputTrajectory(): Error: Seek only supported for g_iBeginStep.\n");
			return false;
		}
	}

	if (g_bVerbose)
		mprintf("<<< SeekInputTrajectory() <<<\n");

	return true;
}


CAnalysisGroup* AddAnalysisGroup(const char *name)
{
	BTIN;

	CAnalysisGroup *g;
	try { g = new CAnalysisGroup(); } catch(...) { g = NULL; }
	if (g == NULL) NewException((double)sizeof(CAnalysisGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { g->m_sGroupName = new char[strlen(name)+1]; } catch(...) { g->m_sGroupName = NULL; }
	if (g->m_sGroupName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(g->m_sGroupName,name);
	g_oaAnalysisGroups.Add(g);
	BTOUT;
	return g;
}


void AddAnalysis(CAnalysisGroup* g, const char *name, const char *abbrev)
{
	BTIN;

	CAnalysis *a;
	int z;


	for (z=0;z<g_oaAnalyses.GetSize();z++) {
		a = (CAnalysis*)g_oaAnalyses[z];
		if (mystricmp(abbrev,a->m_sAbbrev) == 0) {
			eprintf("\nAddAnalysis(): Internal Error: Found duplicate analysis abbreviation:\n");
			eprintf("\"%s\" - \"%s\"\n",a->m_sAbbrev,a->m_sName);
			eprintf("\"%s\" - \"%s\"\n",abbrev,name);
			abort();
		}
	}

	try { a = new CAnalysis(); } catch(...) { a = NULL; }
	if (a == NULL) NewException((double)sizeof(CAnalysis),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { a->m_sName = new char[strlen(name)+1]; } catch(...) { a->m_sName = NULL; }
	if (a->m_sName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { a->m_sAbbrev = new char[strlen(abbrev)+1]; } catch(...) { a->m_sAbbrev = NULL; }
	if (a->m_sAbbrev == NULL) NewException((double)(strlen(abbrev)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(a->m_sName,name);
	strcpy(a->m_sAbbrev,abbrev);
	g_oaAnalyses.Add(a);
	g->m_oaAnalyses.Add(a);
	BTOUT;
}


void InitAnalyses()
{
	BTIN;
	CAnalysisGroup *g;

	/**********************************************************************/
	g = AddAnalysisGroup("Static (time independent) Functions");

	AddAnalysis(g,"Combined Distribution Function","cdf");
	AddAnalysis(g,"Radial Distribution Function","rdf");
	AddAnalysis(g,"Angular Distribution Function","adf");
	AddAnalysis(g,"Dihedral Distribution Function","ddf");
	AddAnalysis(g,"Point-Plane Distance Distribution","pldf");
	AddAnalysis(g,"Point-Line Distance Distribution","lidf");

	AddAnalysis(g,"Plane Projection Distribution","plproj");
	AddAnalysis(g,"Fixed Plane Density Profile","dprof");
	AddAnalysis(g,"Fixed Plane Projection","fpp");
	AddAnalysis(g,"Density Distribution Function","dens");
	AddAnalysis(g,"Spatial Distribution Function","sdf");
	AddAnalysis(g,"Pseudo SDF (only 2 ref. atoms)","psdf");
	AddAnalysis(g,"Dipole Distribution Function","dip");
	AddAnalysis(g,"Evaluate Structural Condition","cond");
	AddAnalysis(g,"Compute Structure Factor","sfac");
	AddAnalysis(g,"Aggregation Topology Analysis", "agtopo");
	AddAnalysis(g,"Sankey diagram plot", "sankey");
	AddAnalysis(g,"Contact Matrix / Connection Matrix", "cmat");
	AddAnalysis(g,"Particle Density inside Geometric Object", "geodens");



	/**********************************************************************/
	g = AddAnalysisGroup("Dynamic (time dependent) Functions");

	AddAnalysis(g,"Velocity Distribution Function","vdf");
	AddAnalysis(g,"Force Distribution Function","fdf");
	AddAnalysis(g,"Mean Square Displacement / Diffusion Coefficients","msd");
	AddAnalysis(g,"Velocity Autocorrelation Functions","acf");
	AddAnalysis(g,"Vector Reorientation Dynamics","rdyn");
	AddAnalysis(g,"Van Hove Correlation Function","vhcf");
	AddAnalysis(g,"Aggregation Functions (DACF, DLDF, DDisp)","aggr");
	AddAnalysis(g,"Hydrogen Bond Dynamics", "hbond");
	AddAnalysis(g,"Diffusion Corrected Ion Pair Dynamics", "ionpair");
	AddAnalysis(g,"Time Dependent Distribution Functions","tddf");



	/**********************************************************************/
	g = AddAnalysisGroup("Spectroscopic Functions");

	AddAnalysis(g, "New Spectroscopy Module (IR, Raman, VCD, ROA)", "spec");
	AddAnalysis(g,"Calculate Power Spectrum","power");
	AddAnalysis(g, "Normal Coordinate Analysis", "nc");
	AddAnalysis(g,"(old) Calculate IR Spectrum","ir");
	AddAnalysis(g,"(old) Calculate Raman Spectrum","raman");
	AddAnalysis(g, "(old) Calculate VCD Spectrum", "vcd");
	AddAnalysis(g, "Calculate Resonance Raman Spectrum", "rera");
	AddAnalysis(g, "Save Dipole Restart File", "drst");
	AddAnalysis(g, "Save Magnetic Moment Restart File", "mrst");
	AddAnalysis(g, "Set up Polarizability Calculation", "pol");

	/**********************************************************************/
	g = AddAnalysisGroup("Miscellaneous Functions");

	AddAnalysis(g,"Save Trajectory of RM Environment / TDO Plot","env");
	AddAnalysis(g,"Save Processed Trajectory","proc");
	AddAnalysis(g,"Cut Clusters","cut");
	AddAnalysis(g,"Region-specific Analysis","reg");
	AddAnalysis(g,"Chirality Analysis", "chi");
	AddAnalysis(g,"Transform to Eckart Frame", "eck");
	AddAnalysis(g,"Sort Wannier Centers", "swan");
	AddAnalysis(g,"Basic Voronoi Analysis","voro");
	AddAnalysis(g,"Voronoi Integration Functions","vori");
	AddAnalysis(g,"Voronoi Charge Optimization","vorochg");
	AddAnalysis(g,"Domain Analysis","doma");
	AddAnalysis(g,"Order Parameters","order");
	AddAnalysis(g,"Hole Theory based Analysis","hole");
	AddAnalysis(g,"Cluster Analysis","cla");

	/**********************************************************************/
	BTOUT;
}


void DumpAnalyses()
{
	BTIN;
	int z, z2;
	CAnalysis *a;
	CAnalysisGroup *g;

	for (z=0;z<g_oaAnalysisGroups.GetSize();z++)
	{
		g = (CAnalysisGroup*)g_oaAnalysisGroups[z];
		mprintf(YELLOW," *** %s\n",g->m_sGroupName);
		for (z2=0;z2<g->m_oaAnalyses.GetSize();z2++)
		{
			a = (CAnalysis*)g->m_oaAnalyses[z2];
			mprintf(WHITE," %-9s",a->m_sAbbrev);
			mprintf("- %s\n",a->m_sName);
		}
		mprintf("\n");
	}
	BTOUT;
}


/*void UniteNb()
{
	BTIN;
	int z, z2, z3, z4;
	CNbSearch *nb;
	CxIntArray *w, *w2;

	if (g_pNbAll != NULL)
		delete g_pNbAll;
	g_pNbAll = new CNbSearch();
	g_pNbAll->Create();

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		w2 = (CxIntArray*)g_pNbAll->m_oaScanNeighbors[z];
		for (z2=0;z2<g_oaNbSearches.GetSize();z2++)
		{
			nb = (CNbSearch*)g_oaNbSearches[z2];
			w = (CxIntArray*)nb->m_oaScanNeighbors[z];
			for (z3=0;z3<nb->m_waScanNeighborCount[z];z3++)
			{
				if (!g_bKeepNbCount)
					for (z4=0;z4<g_pNbAll->m_waScanNeighborCount[z];z4++)
					{
						if (w->GetAt(z3) == w2->GetAt(z4))
							goto _unb_found;
					}
				w2->Add(w->GetAt(z3));
				g_pNbAll->m_waScanNeighborCount[z]++;
_unb_found:;
			}
		}
	}
	BTOUT; 
}*/


bool ParseAtom(const char *s, int refmol, int &ty, int &rty, int &atom)
{
	BTIN;
	CMolecule *mol;
	char buf[64];
	const char *r, *t;
	int a, b;

	mol = (CMolecule*)g_oaMolecules[refmol];
	t = s;
	while ((*t == ' ') /* && (*t != 0) */ )
		t++;
	r = t;
	while ((!isdigit(*r)) && (*r != 0))
		r++;
	if (r-t == 0)
	{
		eprintf("ParseAtom(): Missing Atom Label for Atom.\n");
		BTOUT;
		return false;
	}
	if (r-t >= 64)
	{
		eprintf("Internal Error in ParseAtom(): A, %ld >= 64.\n",(long)(r-t));
		return false;
	}
	memcpy(buf,t,r-t);
	buf[r-t] = 0;
//	printf("Label \"%s\"\n",buf);
	a = mol->FindAtomInMol(buf);
	if (a == -1)
	{
		eprintf("ParseAtom(): Atom Type \"%s\" not found in Molecule.\n",buf);
		BTOUT;
		return false;
	}
	b = atoi(r);
	if (b != 0)
		b--;
	if (b >= ((CMolecule*)g_oaMolecules[refmol])->m_waAtomCount[a])
	{
		eprintf("ParseAtom(): The Molecule only has %d %s-Atoms.\n",((CMolecule*)g_oaMolecules[refmol])->m_waAtomCount[a],buf);
		BTOUT;
		return false;
	}
//	printf("ParseRefSystem: Atom 1 ist in Molekuel %d Typ %d Nummer %d.\n",refmol,a,b);
	ty = a;
	rty = ((CMolecule*)g_oaMolecules[refmol])->m_baAtomIndex[a];
	atom = b;
	BTOUT;
	return true;
}


bool ParseRefSystem(int refmol, const char *s, int points)
{
	BTIN;
	char buf[256];
	char *p, *q;

	if (points == 1)
	{
		BTOUT;
		return ParseAtom(s,refmol,g_iFixAtomType[0],g_iFixRealAtomType[0],g_iFixAtom[0]);
	}
	strcpy(buf,s);
	if (points == 2)
	{
		p = strchr(buf,',');
		if (p == NULL)
		{
			eprintf("ParseRefSystem(): No comma found.\n");
			BTOUT;
			return false;
		}
		*p = 0;
		p++;
		if  (!ParseAtom(buf,refmol,g_iFixAtomType[0],g_iFixRealAtomType[0],g_iFixAtom[0]))
		{
			BTOUT;
			return false;
		}
		BTOUT;
		return ParseAtom(p,refmol,g_iFixAtomType[1],g_iFixRealAtomType[1],g_iFixAtom[1]);
	}
	if (points == 3)
	{
		p = strchr(buf,',');
		if (p == NULL)
		{
			eprintf("ParseRefSystem(): No comma found.\n");
			BTOUT;
			return false;
		}
		*p = 0;
		p++;
		q = strchr(p,',');
		if (q == NULL)
		{
			eprintf("ParseRefSystem(): No second comma found.\n");
			BTOUT;
			return false;
		}
		*q = 0;
		q++;
		if  (!ParseAtom(buf,refmol,g_iFixAtomType[0],g_iFixRealAtomType[0],g_iFixAtom[0]))
		{
			BTOUT;
			return false;
		}
		if  (!ParseAtom(p,refmol,g_iFixAtomType[1],g_iFixRealAtomType[1],g_iFixAtom[1]))
		{
			BTOUT;
			return false;
		}
		BTOUT;
		return ParseAtom(q,refmol,g_iFixAtomType[2],g_iFixRealAtomType[2],g_iFixAtom[2]);
	}
	BTOUT;
	return false;
}


CTimeStep* GetTimeStep(int i)
{
	BTIN;
	if (g_iCurrentTimeStep != -1)
	{
		if (g_bUseVelocities || g_bUseForces)
			i++;
		if (i >= g_iStepHistory) {
			eprintf("GetTimeStep(): Error! Requested depth (%d) >= step history (%d).\n",i,g_iStepHistory);
			BTOUT;
			return NULL;
		}
		if (g_iCurrentTimeStep-i+g_iStepHistory < 0)
		{
			eprintf("GetTimeStep(): Error! i=%d, g_iCurrentTimeStep=%d, g_iStepHistory=%d.\n",i,g_iCurrentTimeStep,g_iStepHistory);
			BTOUT;
			return NULL;
		}
		if (g_iCurrentTimeStep-i >= 0)
		{
			BTOUT;
			return (CTimeStep*)g_oaTimeSteps[g_iCurrentTimeStep-i];
		} else
		{
			BTOUT;
			return (CTimeStep*)g_oaTimeSteps[g_iCurrentTimeStep-i+g_iStepHistory];
		}
	} else 
	{
		BTOUT;
		return NULL;
	}
}


CTimeStep** GetTimeStepAddress(int i)
{
	BTIN;
	if (g_iCurrentTimeStep != -1)
	{
		if (g_bUseVelocities || g_bUseForces)
			i++;
		if (g_iCurrentTimeStep-i >= 0)
		{
			BTOUT;
			return (CTimeStep**)&g_oaTimeSteps[g_iCurrentTimeStep-i];
		} else
		{
			BTOUT;
			return (CTimeStep**)&g_oaTimeSteps[g_iCurrentTimeStep-i+g_iStepHistory];
		}
	} else
	{
		BTOUT;
		return NULL;
	}
}


void CalcVelocities()
{
	BTIN;
	int z, z2;
	unsigned long index;
	double v, imv, ilmv;
	CxDVector3 dist;
	CMolecule *m;
	CSingleMolecule *sm;
	CVirtualAtom *vi;

	g_fLMidVel = 0;
	imv = g_fMaxVel;
	ilmv = g_fLMaxVel;
	GetTimeStep(0)->m_vaVelocities.SetSize(g_iGesVirtAtomCount);
// 	for (z=0;z<g_iGesVirtAtomCount;z++)
	for (z=0;z<g_iGesAtomCount;z++)
	{
//		if (z == 0)
//			mprintf("\n(%f - %f) / (2 * %f) * 100000 = %f.",GetTimeStep(1)->m_vaCoords[z][0],GetTimeStep(-1)->m_vaCoords[z][0],g_fTimestepLength,(GetTimeStep(1)->m_vaCoords[z][0] - GetTimeStep(-1)->m_vaCoords[z][0]) / (2.0f*g_fTimestepLength) * 1000.0f);
		
		dist = GetTimeStep(-1)->m_vaCoords[z] - GetTimeStep(1)->m_vaCoords[z];
// 		mprintf(RED, "d: %d ( %12G | %12G | %12G )\n", g_waAtomElement[z], dist[0], dist[1], dist[2]);
		if(g_bPeriodicX) {
			while(dist[0] > g_fBoxX / 2.0) dist[0] -= g_fBoxX;
			while(dist[0] < -g_fBoxX / 2.0) dist[0] += g_fBoxX;
		}
		if(g_bPeriodicY) {
			while(dist[1] > g_fBoxY / 2.0) dist[1] -= g_fBoxY;
			while(dist[1] < -g_fBoxY / 2.0) dist[1] += g_fBoxY;
		}
		if(g_bPeriodicZ) {
			while(dist[2] > g_fBoxZ / 2.0) dist[2] -= g_fBoxZ;
			while(dist[2] < -g_fBoxZ / 2.0) dist[2] += g_fBoxZ;
		}
		GetTimeStep(0)->m_vaVelocities[z] = dist * 1000.0 / 2.0 / g_fTimestepLength / g_iStride;
// 		mprintf(RED, "v: ( %12G | %12G | %12G )\n", GetTimeStep(0)->m_vaVelocities[z][0], GetTimeStep(0)->m_vaVelocities[z][1], GetTimeStep(0)->m_vaVelocities[z][2]);
		
// 		GetTimeStep(0)->m_vaVelocities[z] = (GetTimeStep(1)->m_vaCoords[z] - GetTimeStep(-1)->m_vaCoords[z]) / (2.0f*g_fTimestepLength) * 1000.0f;
		v = GetTimeStep(0)->m_vaVelocities[z].GetLength();
		g_fLMidVel += v;
		if (v > imv)
			imv = v;
		if (v > ilmv)
			ilmv = v;
	}
	for (z = 0; z < g_oaVirtualAtoms.GetSize(); z++) {
		vi = (CVirtualAtom*)g_oaVirtualAtoms[z];
		m = (CMolecule*)g_oaMolecules[vi->m_iMolecule];
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++) {
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];

			index = ((CxIntArray*)sm->m_oaAtomOffset[sm->m_baAtomIndex.GetSize()-1])->GetAt(vi->m_iMolVirtAtom);
			if (((CVirtualAtom *)g_oaVirtualAtoms[z])->m_iMode == 2) {
				dist = CxDVector3(0.0, 0.0, 0.0);
			} else {
				dist = GetTimeStep(-1)->m_vaCoords[index] - GetTimeStep(1)->m_vaCoords[index];
		// 		mprintf(RED, "d: %d ( %12G | %12G | %12G )\n", g_waAtomElement[index], dist[0], dist[1], dist[2]);
				if(g_bPeriodicX) {
					while(dist[0] > g_fBoxX / 2.0) dist[0] -= g_fBoxX;
					while(dist[0] < -g_fBoxX / 2.0) dist[0] += g_fBoxX;
				}
				if(g_bPeriodicY) {
					while(dist[1] > g_fBoxY / 2.0) dist[1] -= g_fBoxY;
					while(dist[1] < -g_fBoxY / 2.0) dist[1] += g_fBoxY;
				}
				if(g_bPeriodicZ) {
					while(dist[2] > g_fBoxZ / 2.0) dist[2] -= g_fBoxZ;
					while(dist[2] < -g_fBoxZ / 2.0) dist[2] += g_fBoxZ;
				}
			}
			GetTimeStep(0)->m_vaVelocities[index] = dist * 1000.0 / 2.0 / g_fTimestepLength / g_iStride;
	//		mprintf("\nVA[%lu/%d]: (%10g|%10g|%10g)",index,g_oaVirtualAtoms.GetSize(),GetTimeStep(0)->m_vaVelocities[index][0],GetTimeStep(0)->m_vaVelocities[index][1],GetTimeStep(0)->m_vaVelocities[index][2]);
	// 		mprintf(RED, "v: ( %12G | %12G | %12G )\n", GetTimeStep(0)->m_vaVelocities[index][0], GetTimeStep(0)->m_vaVelocities[index][1], GetTimeStep(0)->m_vaVelocities[index][2]);
			v = GetTimeStep(0)->m_vaVelocities[index].GetLength();
			g_fLMidVel += v;
			if (v > imv)
				imv = v;
			if (v > ilmv)
				ilmv = v;
		}
	}
	g_fLMidVel /= g_iGesVirtAtomCount;
	g_fLMaxVel = ilmv;
	if (imv < g_fUnsteadyLimit)
	{
		g_fMaxVel = imv;
	} else if (!g_bWarnUnsteady)
	{
		g_bWarnUnsteady = true;
		mprintf("Warning: Discontinuity in time step %d (maximum velocity = %f > %f pm/ps).\n",(int)g_iSteps,imv,g_fUnsteadyLimit);
	}
	BTOUT; 
}

void CalcVolumetricDataTimeDev() {
	if (GetTimeStep(0)->m_pVolumetricDataTimeDev == NULL)
		return;
	if (GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin == NULL) {
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] = GetTimeStep(0)->m_pVolumetricData->m_iRes[0];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] = GetTimeStep(0)->m_pVolumetricData->m_iRes[1];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] = GetTimeStep(0)->m_pVolumetricData->m_iRes[2];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMinVal[0] = GetTimeStep(0)->m_pVolumetricData->m_fMinVal[0];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMaxVal[0] = GetTimeStep(0)->m_pVolumetricData->m_fMaxVal[0];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMinVal[1] = GetTimeStep(0)->m_pVolumetricData->m_fMinVal[1];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMaxVal[1] = GetTimeStep(0)->m_pVolumetricData->m_fMaxVal[1];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMinVal[2] = GetTimeStep(0)->m_pVolumetricData->m_fMinVal[2];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->m_fMaxVal[2] = GetTimeStep(0)->m_pVolumetricData->m_fMaxVal[2];
		GetTimeStep(0)->m_pVolumetricDataTimeDev->Create();
	}
	
	int i;
	double integral = 0.0;
	int timedevSize = 0;
	for (i = 0; i < GetTimeStep(0)->m_pVolumetricData->m_iRes[0] * GetTimeStep(0)->m_pVolumetricData->m_iRes[1] * GetTimeStep(0)->m_pVolumetricData->m_iRes[2]; i++) {
		if (fabs(GetTimeStep(0)->m_pVolumetricData->m_pBin[i]) > -1.0e-20) {
			GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i] = (GetTimeStep(-1)->m_pVolumetricData->m_pBin[i] - GetTimeStep(1)->m_pVolumetricData->m_pBin[i]) / (2.0 * g_fTimestepLength);
			integral += GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i];
			timedevSize++;
		} else {
			GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i] = 0.0;
		}
	}
// 	mprintf(GREEN, "%g %g\n", integral, integral / timedevSize);
	for (i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2]; i++) {
		if (fabs(GetTimeStep(0)->m_pVolumetricData->m_pBin[i]) > -1.0e-20) {
			GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i] -= integral / timedevSize;
		}
	}
// 	integral = 0.0;
// 	for (i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2]; i++) {
// 		integral += GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i];
// 	}
// 	mprintf(GREEN, "%g\n", integral);
	
// 	mprintf(GREEN, "%#20.10g %#20.10g %#20.10g %#20.10g\n", GetTimeStep(-1)->m_pVolumetricData->m_pBin[10], GetTimeStep(0)->m_pVolumetricData->m_pBin[10], GetTimeStep(1)->m_pVolumetricData->m_pBin[10], GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[10]);
	
// 		FILE *cubeFile = fopen("test.cube", "w");
// 		if (cubeFile == NULL) {
// 			printf("Could not open output file!\n");
// 			abort();
// 		}
// 		
// 		fprintf(cubeFile, "\n\n");
// 		fprintf(cubeFile, "%5lu %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_iGesAtomCount, 0.0f, 0.0f, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0], g_fCubeXStep, 0.0f, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1], 0.0f, g_fCubeYStep, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2], 0.0f, 0.0f, g_fCubeZStep);
// 		for (int i = 0; i < (int)GetTimeStep(0)->m_iGesAtomCount; i++) {
// 			fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f %12.6f\n", 1, 0.0f, GetTimeStep(0)->m_vaCoords[i][0] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][1] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][2] / 0.529177249f / 100.0f);
// 		}
// 		for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 			for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 				for (int k = 0; k < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6; k++) {
// 					for (int l = 0; l < 6; l++) {
// 						fprintf(cubeFile, "%13.5E", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + (k * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]]);
// 					}
// 					fprintf(cubeFile, "\n");
// 				}
// 				if (GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6 != 0) {
// 					for (int l = 0; l < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6; l++) {
// 						fprintf(cubeFile, "%13.5E", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + ((GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6) * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]]);
// 					}
// 					fprintf(cubeFile, "\n");
// 				}
// 			}
// 		}
// 		
// 		fclose(cubeFile);
}

void CalcCurrentDensity() {
	int res[3];
	res[0] = GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0];
	res[1] = GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1];
	res[2] = GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2];
	
	C3DF<VORI_FLOAT> density;
	density.CopyFrom(GetTimeStep(0)->m_pVolumetricData);
	
	C3DF<VORI_FLOAT> densityGradient[3];
	int i;
	for (i = 0; i < 3; i++) {
		densityGradient[i].m_iRes[0] = density.m_iRes[0];
		densityGradient[i].m_iRes[1] = density.m_iRes[1];
		densityGradient[i].m_iRes[2] = density.m_iRes[2];
		densityGradient[i].m_fMinVal[0] = density.m_fMinVal[0];
		densityGradient[i].m_fMaxVal[0] = density.m_fMaxVal[0];
		densityGradient[i].m_fMinVal[1] = density.m_fMinVal[1];
		densityGradient[i].m_fMaxVal[1] = density.m_fMaxVal[1];
		densityGradient[i].m_fMinVal[2] = density.m_fMinVal[2];
		densityGradient[i].m_fMaxVal[2] = density.m_fMaxVal[2];
		densityGradient[i].Create();
	}
	
	for (i = 0; i < res[2]; i++) {
		int j;
		for (j = 0; j < res[1]; j++) {
			int k;
			for (k = 1; k < res[0] - 1; k++) {
				densityGradient[0].m_pBin[i * res[0] * res[1] + j * res[0] + k] = (density.m_pBin[i * res[0] * res[1] + j * res[0] + k + 1] - density.m_pBin[i * res[0] * res[1] + j * res[0] + k - 1]) / (2.0 * g_fCubeXStep);
			}
			densityGradient[0].m_pBin[i * res[0] * res[1] + j * res[0]] = (density.m_pBin[i * res[0] * res[1] + j * res[0] + 1] - density.m_pBin[i * res[0] * res[1] + j * res[0] + res[0] - 1]) / (2.0 * g_fCubeXStep);
			densityGradient[0].m_pBin[i * res[0] * res[1] + j * res[0] + res[0] - 1] = (density.m_pBin[i * res[0] * res[1] + j * res[0]] - density.m_pBin[i * res[0] * res[1] + j * res[0] + res[0] - 2]) / (2.0 * g_fCubeXStep);
		}
	}
	
	for (i = 0; i < res[2]; i++) {
		int j;
		for (j = 1; j < res[1] - 1; j++) {
			int k;
			for (k = 0; k < res[0]; k++) {
				densityGradient[1].m_pBin[i * res[0] * res[1] + j * res[0] + k] = (density.m_pBin[i * res[0] * res[1] + (j + 1) * res[0] + k] - density.m_pBin[i * res[0] * res[1] + (j - 1) * res[0] + k]) / (2.0 * g_fCubeYStep);
			}
		}
		int k;
		for (k = 0; k < res[0]; k++) {
			densityGradient[1].m_pBin[i * res[0] * res[1] + k] = (density.m_pBin[i * res[0] * res[1] + res[0] + k] - density.m_pBin[i * res[0] * res[1] + (res[1] - 1) * res[0] + k]) / (2.0 * g_fCubeYStep);
			densityGradient[1].m_pBin[i * res[0] * res[1] + (res[1] - 1) * res[0] + k] = (density.m_pBin[i * res[0] * res[1] + k] - density.m_pBin[i * res[0] * res[1] + (res[1] - 2) * res[0] + k]) / (2.0 * g_fCubeYStep);
		}
	}
	
	for (i = 1; i < res[2] - 1; i++) {
		int j;
		for (j = 0; j < res[1]; j++) {
			int k;
			for (k = 0; k < res[0]; k++) {
				densityGradient[2].m_pBin[i * res[0] * res[1] + j * res[0] + k] = (density.m_pBin[(i + 1) * res[0] * res[1] + j * res[0] + k] - density.m_pBin[(i - 1) * res[0] * res[1] + j * res[0] + k]) / (2.0 * g_fCubeZStep);
			}
		}
	}
	int j;
	for (j = 0; j < res[1]; j++) {
		int k;
		for (k = 0; k < res[0]; k++) {
			densityGradient[2].m_pBin[j * res[0] + k] = (density.m_pBin[res[0] * res[1] + j * res[0] + k] - density.m_pBin[(res[2] - 1) * res[0] * res[1] + j * res[0] + k]) / (2.0 * g_fCubeZStep);
			densityGradient[2].m_pBin[(res[2] - 1) * res[0] * res[1] + j * res[0] + k] = (density.m_pBin[j * res[0] + k] - density.m_pBin[(res[2] - 2) * res[0] * res[1] + j * res[0] + k]) / (2.0 * g_fCubeZStep);
		}
	}
// 	mprintf(GREEN, "%#.10g %#.10g %#.10g %#.10g\n", density.m_pBin[1942857], densityGradient[0].m_pBin[1942857], densityGradient[1].m_pBin[1942857], densityGradient[2].m_pBin[1942857]);
	
// 	FILE *cubeFile2 = fopen("test2.cube", "w");
// 	if (cubeFile2 == NULL) {
// 		printf("Could not open output file!\n");
// 		abort();
// 	}
// 	
// 	fprintf(cubeFile2, "\n\n");
// 	fprintf(cubeFile2, "%5lu %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_iGesAtomCount, 0.0f, 0.0f, 0.0f);
// 	fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0], g_fCubeXStep, 0.0f, 0.0f);
// 	fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1], 0.0f, g_fCubeYStep, 0.0f);
// 	fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2], 0.0f, 0.0f, g_fCubeZStep);
// 	for (int i = 0; i < (int)GetTimeStep(0)->m_iGesAtomCount; i++) {
// 		fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f %12.6f\n", 1, 0.0f, GetTimeStep(0)->m_vaCoords[i][0] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][1] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][2] / 0.529177249f / 100.0f);
// 	}
// 	for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 		for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 			for (int k = 0; k < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6; k++) {
// 				for (int l = 0; l < 6; l++) {
// 					fprintf(cubeFile2, "%13.5E", densityGradient[2].m_pBin[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + (k * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]]);
// 				}
// 				fprintf(cubeFile2, "\n");
// 			}
// 			if (GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6 != 0) {
// 				for (int l = 0; l < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6; l++) {
// 					fprintf(cubeFile2, "%13.5E", densityGradient[2].m_pBin[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + ((GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6) * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]]);
// 				}
// 				fprintf(cubeFile2, "\n");
// 			}
// 		}
// 	}
// 	
// 	fclose(cubeFile2);
	
// 	CSparseMatrix testMatrix;
// 	testMatrix.setTest();
// 	double solution[9];
// 	for (int i = 0; i < 9; i++)
// 		solution[i] = 0.0;
// 	double rhs[9];
// 	rhs[0] = 1.0;
// 	rhs[1] = 2.0;
// 	rhs[2] = 1.0;
// 	rhs[3] = 2.0;
// 	rhs[4] = 1.0;
// 	rhs[5] = 2.0;
// 	rhs[6] = 1.0;
// 	rhs[7] = 2.0;
// 	rhs[8] = 1.0;
// 	
// 	CCurrentPDESolver::bicgstabl(4, &testMatrix, solution, rhs, 100, 1e-10, stdout);
// 	
// 	return;
	
	for (i = 0; i < res[0] * res[1] * res[2]; i++) {
		density.m_pBin[i] += g_fBackgroundDensity;
	}
	
	CSparseMatrix pdeMatrix;
	CCurrentPDEDiscretizer::discretize(&pdeMatrix, density, densityGradient[0], densityGradient[1], densityGradient[2]);
	
	static double *pdeSolution = NULL;
	if (pdeSolution == NULL) {
		try { pdeSolution = new double[res[0] * res[1] * res[2]]; } catch(...) { pdeSolution = NULL; }
		if (pdeSolution == NULL) NewException((double)res[0] * res[1] * res[2] * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		memset(pdeSolution, 0, res[0] * res[1] * res[2] * sizeof(double));
	}
	
	static bool first = true;
	static double thresh;
	if (first) {
		memset(pdeSolution, 0, res[0] * res[1] * res[2] * sizeof(double));
		thresh = g_fPDEConvThresh * CCurrentPDESolver::calcResidual(&pdeMatrix, pdeSolution, GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin);
		first = false;
	}
	
	if (g_fPDESolverInfoFile != NULL) {
		fprintf(g_fPDESolverInfoFile, "Step %lu\n", g_iSteps - 1);
	}
	double thresh2 = thresh;
	if (!CCurrentPDESolver::bicgstabl(4, &pdeMatrix, pdeSolution, GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin, g_iPDEMaxIter, &thresh2, g_fPDESolverInfoFile)) {
		memset(pdeSolution, 0, res[0] * res[1] * res[2] * sizeof(double));
		if (g_fPDESolverInfoFile != NULL) {
			fprintf(g_fPDESolverInfoFile, "Resetting solution guess\n");
		}
		thresh2 = thresh;
		if (!CCurrentPDESolver::bicgstabl(4, &pdeMatrix, pdeSolution, GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin, g_iPDEMaxIter, &thresh2, g_fPDESolverInfoFile)) {
			memset(pdeSolution, 0, res[0] * res[1] * res[2] * sizeof(double));
			if (g_fPDESolverInfoFile != NULL) {
				fprintf(g_fPDESolverInfoFile, "Trying different threshold\n");
			}
			thresh2 *= 1.05;
			if (!CCurrentPDESolver::bicgstabl(4, &pdeMatrix, pdeSolution, GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin, g_iPDEMaxIter, &thresh2, g_fPDESolverInfoFile)) {
				memset(pdeSolution, 0, res[0] * res[1] * res[2] * sizeof(double));
				if (g_fPDESolverInfoFile != NULL) {
					fprintf(g_fPDESolverInfoFile, "Severe convergence problem! Resetting solution guess\n");
				}
			}
		}
	}
	if (g_fPDESolverInfoFile != NULL) {
		fprintf(g_fPDESolverInfoFile, "\n");
	}
	
	GetTimeStep(0)->m_pCurrentDensity->SetSize(3 * res[0] * res[1] * res[2]);
	
// 	fftwf_complex *fft_data;
// 	fft_data = fftwf_alloc_complex(res[0] * res[1] * res[2]);
// 	fftwf_plan plan1, plan2;
// 	plan1 = fftwf_plan_dft_3d(res[2], res[1], res[0], fft_data, fft_data, FFTW_FORWARD, FFTW_MEASURE);
// 	plan2 = fftwf_plan_dft_3d(res[2], res[1], res[0], fft_data, fft_data, FFTW_BACKWARD, FFTW_MEASURE);
// 	
// 	int i;
// 	for (i = 0; i < res[0] * res[1] * res[2]; i++) {
// 		fft_data[i][0] = GetTimeStep(0)->m_pVolumetricDataTimeDev->m_pBin[i];
// 		fft_data[i][1] = 0.0f;
// 	}
// 	fftwf_execute(plan1);
// 	for (i = 0; i < res[2]; i++) {
// 		int j;
// 		for (j = 0; j < res[1]; j++) {
// 			int k;
// 			for (k = 0; k < res[0]; k++) {
// 				if (i == 0 && j == 0 && k == 0) {
// 					fft_data[i * res[1] * res[0] + j * res[0] + k][0] = 0.0f;
// 					fft_data[i * res[1] * res[0] + j * res[0] + k][1] = 0.0f;
// 					continue;
// 				}
// 				double f = (2.0f * cos(2.0f * (double)Pi / res[2] * i) - 2.0) / g_fCubeXStep / g_fCubeXStep + (2.0f * cos(2.0f * (double)Pi / res[1] * j) - 2.0f) / g_fCubeYStep / g_fCubeYStep + (2.0f * cos(2.0f * (double)Pi / res[0] * k) - 2.0f) / g_fCubeZStep / g_fCubeZStep;
// 				fft_data[i * res[1] * res[0] + j * res[0] + k][0] /= f;
// 				fft_data[i * res[1] * res[0] + j * res[0] + k][1] /= f;
// 			}
// 		}
// 	}
// 	fftwf_execute(plan2);
// 	for (i = 0; i < res[0] * res[1] * res[2]; i++) {
// 		fft_data[i][0] /= res[0] * res[1] * res[2];
// 		fft_data[i][1] /= res[0] * res[1] * res[2];
// 	}
	
// 		FILE *cubeFile2 = fopen("test2.cube", "w");
// 		if (cubeFile2 == NULL) {
// 			printf("Could not open output file!\n");
// 			abort();
// 		}
// 		
// 		fprintf(cubeFile2, "\n\n");
// 		fprintf(cubeFile2, "%5lu %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_iGesAtomCount, 0.0f, 0.0f, 0.0f);
// 		fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0], g_fCubeXStep, 0.0f, 0.0f);
// 		fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1], 0.0f, g_fCubeYStep, 0.0f);
// 		fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2], 0.0f, 0.0f, g_fCubeZStep);
// 		for (int i = 0; i < (int)GetTimeStep(0)->m_iGesAtomCount; i++) {
// 			fprintf(cubeFile2, "%5d %12.6f %12.6f %12.6f %12.6f\n", 1, 0.0f, GetTimeStep(0)->m_vaCoords[i][0] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][1] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][2] / 0.529177249f / 100.0f);
// 		}
// 		for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 			for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 				for (int k = 0; k < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6; k++) {
// 					for (int l = 0; l < 6; l++) {
// 						fprintf(cubeFile2, "%13.5E", fft_data[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + (k * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]][0]);
// 					}
// 					fprintf(cubeFile2, "\n");
// 				}
// 				if (GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6 != 0) {
// 					for (int l = 0; l < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6; l++) {
// 						fprintf(cubeFile2, "%13.5E", fft_data[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + ((GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6) * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]][0]);
// 					}
// 					fprintf(cubeFile2, "\n");
// 				}
// 			}
// 		}
// 		
// 		fclose(cubeFile2);
// 	
// 		FILE *testFile2 = fopen("test2.dat", "w");
// 		if (testFile2 == NULL) {
// 			printf("Could not open output file!\n");
// 			abort();
// 		}
// 		
// 		for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 			for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 				fprintf(testFile2, "%.6f %.6f %.14f\n", i * g_fCubeXStep, j * g_fCubeYStep, fft_data[i + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] + GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 2 * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]][0]);
// 			}
// 			fprintf(testFile2, "\n");
// 		}
// 		
// 		fclose(testFile2);
		
// 	for (i = 0; i < res[2]; i++) {
// 		int j;
// 		for (j = 0; j < res[1]; j++) {
// 			int k;
// 			for (k = 1; k < res[0] - 1; k++) {
// 				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3) = (fft_data[i * res[1] * res[0] + j * res[0] + k + 1][0] - fft_data[i * res[1] * res[0] + j * res[0] + k - 1][0]) / g_fCubeXStep;
// 			}
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3) = (fft_data[i * res[1] * res[0] + j * res[0] + 1][0] - fft_data[i * res[1] * res[0] + j * res[0] + res[0] - 1][0]) / g_fCubeXStep;
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + (res[0] - 1) * 3) = (fft_data[i * res[1] * res[0] + j * res[0]][0] - fft_data[i * res[1] * res[0] + j * res[0] + res[0] - 2][0]) / g_fCubeXStep;
// 		}
// 	}
// 	for (i = 0; i < res[2]; i++) {
// 		int j;
// 		for (j = 1; j < res[1] - 1; j++) {
// 			int k;
// 			for (k = 0; k < res[0]; k++) {
// 				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 1) = (fft_data[i * res[1] * res[0] + (j + 1) * res[0] + k][0] - fft_data[i * res[1] * res[0] + (j - 1) * res[0] + k][0]) / g_fCubeYStep;
// 			}
// 		}
// 		int k;
// 		for (k = 0; k < res[0]; k++) {
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + k * 3 + 1) = (fft_data[i * res[1] * res[0] + res[0] + k][0] - fft_data[i * res[1] * res[0] + (res[1] - 1) * res[0] + k][0]) / g_fCubeYStep;
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + (res[1] - 1) * res[0] * 3 + k * 3 + 1) = (fft_data[i * res[1] * res[0] + k][0] - fft_data[i * res[1] * res[0] + (res[1] - 2) * res[0] + k][0]) / g_fCubeYStep;
// 		}
// 	}
// 	for (i = 1; i < res[2] - 1; i++) {
// 		int j;
// 		for (j = 0; j < res[1]; j++) {
// 			int k;
// 			for (k = 0; k < res[0]; k++) {
// 				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 2) = (fft_data[(i + 1) * res[1] * res[0] + j * res[0] + k][0] - fft_data[(i - 1) * res[1] * res[0] + j * res[0] + k][0]) / g_fCubeZStep;
// 			}
// 		}
// 	}
// 	int j;
// 	for (j = 0; j < res[1]; j++) {
// 		int k;
// 		for (k = 0; k < res[0]; k++) {
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt(j * res[0] * 3 + k * 3 + 2) = (fft_data[res[1] * res[0] + j * res[0] + k][0] - fft_data[(res[2] - 1) * res[1] * res[0] + j * res[0] + k][0]) / g_fCubeZStep;
// 			GetTimeStep(0)->m_pCurrentDensity->GetAt((res[2] - 1) * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 2) = (fft_data[j * res[0] + k][0] - fft_data[(res[2] - 2) * res[1] * res[0] + j * res[0] + k][0]) / g_fCubeZStep;
// 		}
// 	}
// 	
// 	fftwf_destroy_plan(plan1);
// 	fftwf_destroy_plan(plan2);
// 	fftwf_free(fft_data);
	
	for (i = 0; i < res[2]; i++) {
		for (j = 0; j < res[1]; j++) {
			int k;
			for (k = 1; k < res[0] - 1; k++) {
				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3) = (pdeSolution[i * res[1] * res[0] + j * res[0] + k + 1] - pdeSolution[i * res[1] * res[0] + j * res[0] + k - 1]) / (2.0 * g_fCubeXStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + k];
			}
			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3) = (pdeSolution[i * res[1] * res[0] + j * res[0] + 1] - pdeSolution[i * res[1] * res[0] + j * res[0] + res[0] - 1]) / (2.0 * g_fCubeXStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0]];
			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + (res[0] - 1) * 3) = (pdeSolution[i * res[1] * res[0] + j * res[0]] - pdeSolution[i * res[1] * res[0] + j * res[0] + res[0] - 2]) / (2.0 * g_fCubeXStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + res[0] - 1];
		}
	}
	for (i = 0; i < res[2]; i++) {
		for (j = 1; j < res[1] - 1; j++) {
			int k;
			for (k = 0; k < res[0]; k++) {
				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 1) = (pdeSolution[i * res[1] * res[0] + (j + 1) * res[0] + k] - pdeSolution[i * res[1] * res[0] + (j - 1) * res[0] + k]) / (2.0 * g_fCubeYStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + k];
			}
		}
		int k;
		for (k = 0; k < res[0]; k++) {
			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + k * 3 + 1) = (pdeSolution[i * res[1] * res[0] + res[0] + k] - pdeSolution[i * res[1] * res[0] + (res[1] - 1) * res[0] + k]) / (2.0 * g_fCubeYStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + k];
			GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + (res[1] - 1) * res[0] * 3 + k * 3 + 1) = (pdeSolution[i * res[1] * res[0] + k] - pdeSolution[i * res[1] * res[0] + (res[1] - 2) * res[0] + k]) / (2.0 * g_fCubeYStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + (res[1] - 1) * res[0] + k];
		}
	}
	for (i = 1; i < res[2] - 1; i++) {
		for (j = 0; j < res[1]; j++) {
			int k;
			for (k = 0; k < res[0]; k++) {
				GetTimeStep(0)->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 2) = (pdeSolution[(i + 1) * res[1] * res[0] + j * res[0] + k] - pdeSolution[(i - 1) * res[1] * res[0] + j * res[0] + k]) / (2.0 * g_fCubeZStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + k];
			}
		}
	}
	for (j = 0; j < res[1]; j++) {
		int k;
		for (k = 0; k < res[0]; k++) {
			GetTimeStep(0)->m_pCurrentDensity->GetAt(j * res[0] * 3 + k * 3 + 2) = (pdeSolution[res[1] * res[0] + j * res[0] + k] - pdeSolution[(res[2] - 1) * res[1] * res[0] + j * res[0] + k]) / (2.0 * g_fCubeZStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[j * res[0] + k];
			GetTimeStep(0)->m_pCurrentDensity->GetAt((res[2] - 1) * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 2) = (pdeSolution[j * res[0] + k] - pdeSolution[(res[2] - 2) * res[1] * res[0] + j * res[0] + k]) / (2.0 * g_fCubeZStep) * GetTimeStep(0)->m_pVolumetricData->m_pBin[(res[2] - 1) * res[1] * res[0] + j * res[0] + k];
		}
	}
// 	mprintf(GREEN, "%#.10g\n", GetTimeStep(0)->m_pCurrentDensity->GetAt(0));
	
// 		FILE *cubeFile = fopen("test.cube", "w");
// 		if (cubeFile == NULL) {
// 			printf("Could not open output file!\n");
// 			abort();
// 		}
// 		
// 		fprintf(cubeFile, "\n\n");
// 		fprintf(cubeFile, "%5lu %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_iGesAtomCount, 0.0f, 0.0f, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0], g_fCubeXStep, 0.0f, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1], 0.0f, g_fCubeYStep, 0.0f);
// 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2], 0.0f, 0.0f, g_fCubeZStep);
// 		for (int i = 0; i < (int)GetTimeStep(0)->m_iGesAtomCount; i++) {
// 			fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f %12.6f\n", 1, 0.0f, GetTimeStep(0)->m_vaCoords[i][0] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][1] / 0.529177249f / 100.0f, GetTimeStep(0)->m_vaCoords[i][2] / 0.529177249f / 100.0f);
// 		}
// 		double sum = 0.0;
// 		for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 			for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 				for (int k = 0; k < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6; k++) {
// 					for (int l = 0; l < 6; l++) {
// 						fprintf(cubeFile, "%13.5E", GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + (k * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * 3));
// 						sum += GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + (k * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * 3);
// 					}
// 					fprintf(cubeFile, "\n");
// 				}
// 				if (GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6 != 0) {
// 					for (int l = 0; l < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] % 6; l++) {
// 						fprintf(cubeFile, "%13.5E", GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + ((GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6) * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * 3));
// 						sum += GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + ((GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 6) * 6 + l) * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * 3);
// 					}
// 					fprintf(cubeFile, "\n");
// 				}
// 			}
// 		}
// 		mprintf(GREEN, "%g\n", sum);
// 		
// 		fclose(cubeFile);
// 		
// 		FILE *testFile = fopen("test.dat", "w");
// 		if (testFile == NULL) {
// 			printf("Could not open output file!\n");
// 			abort();
// 		}
// 		
// 		for (int i = 0; i < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0]; i++) {
// 			for (int j = 0; j < GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1]; j++) {
// 				fprintf(testFile, "%.6f %.6f %.14f %.14f\n", i * g_fCubeXStep, j * g_fCubeYStep, GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 2 * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3), GetTimeStep(0)->m_pCurrentDensity->GetAt(i * 3 + j * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[2] / 2 * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[1] * GetTimeStep(0)->m_pVolumetricDataTimeDev->m_iRes[0] * 3 + 1));
// 			}
// 			fprintf(testFile, "\n");
// 		}
// 		
// 		fclose(testFile);
}

void CalcForces()
{
	BTIN;
	int z;
	double f, imf, ilmf;

	g_fLMidForce = 0;
	imf = g_fMaxForce;
	ilmf = g_fLMaxForce;
	GetTimeStep(0)->m_vaForces.SetSize(g_iGesVirtAtomCount);
	for (z=0;z<g_iGesVirtAtomCount;z++)
	{
		GetTimeStep(0)->m_vaForces[z] = (GetTimeStep(-1)->m_vaCoords[z] - 2*GetTimeStep(0)->m_vaCoords[z] + GetTimeStep(1)->m_vaCoords[z])/(double)pow2(g_fTimestepLength/1000.0);
		f = GetTimeStep(0)->m_vaForces[z].GetLength();
// 		mprintf(GREEN, "%g\n", f);
		g_fLMidForce += f;
		if (f > imf)
			imf = f;
		if (f > ilmf)
			ilmf = f;
	}
	g_fLMidForce /= g_iGesVirtAtomCount;
	g_fLMaxForce = ilmf;
	if (imf < g_fUnsteadyLimit)
	{
		g_fMaxForce = imf;
	} else if (!g_bWarnUnsteady)
	{
		g_bWarnUnsteady = true;
		mprintf("Warning: Discontinuity in time step %d (maximum acceleration = %f > %f pm/ps^2).\n",(int)g_iSteps,imf,g_fUnsteadyLimit);
	}
	BTOUT;
}


/*double AtomMass(char *s)
{
	BTIN;
	int z;
	for (z=0;z<g_oaElements.GetSize();z++)
		if (mystricmp(s,((CElement*)g_oaElements[z])->m_sLabel)==0)
		{
			BTOUT;
			return ((CElement*)g_oaElements[z])->m_fMass;
		}
	eprintf("Atom \"%s\" not found.\n",s);
	BTOUT;
	return -1.0f;
}

int AtomOrd(char *s)
{
	BTIN;
	int z;
	for (z=0;z<g_oaElements.GetSize();z++)
		if (mystricmp(s,((CElement*)g_oaElements[z])->m_sLabel)==0)
		{
			BTOUT;
			return ((CElement*)g_oaElements[z])->m_iOrd;
		}
	eprintf("Atom \"%s\" not found.\n",s);
	BTOUT;
	return -1;
}

double AtomRadius(char *s)
{
	BTIN;
	int z;
	for (z=0;z<g_oaElements.GetSize();z++)
		if (mystricmp(s,((CElement*)g_oaElements[z])->m_sLabel)==0)
		{
			BTOUT;
			return ((CElement*)g_oaElements[z])->m_fRadius;
		}
	mprintf("No Atom Radius found for Atom \"%s\".\n",s);
	BTOUT;
	return 0;
}*/


CElement* FindElement(const char *s, bool quiet)
{
	BTIN;
	int z;

	for (z=0;z<g_oaElements.GetSize();z++)
	{
		if (mystricmp(s,((CElement*)g_oaElements[z])->m_sLabel)==0)
		{
			BTOUT;
			return (CElement*)g_oaElements[z];
		}
	}
	if (!quiet)
		eprintf("    No element data found for atom \"%s\".\n",s);
	g_bUnknownElements = true;
	BTOUT;
	return NULL;
}


double GuessBoxSize()
{
	BTIN;
	int z;
	double m, gm;

	gm = 0;
	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		if (z == g_iVirtAtomType)
			continue;
		m = ((CAtom*)g_oaAtoms[z])->m_pElement->m_fMass;
/*		if (m < 0) // Dieses Atom ist nicht in der Liste
		{
			if (z != g_iVirtAtomType)
				mprintf("Atom \"%s\" nicht in der Liste.\n",((CAtom*)g_oaAtoms[z])->m_sName);
			continue;
		}*/
		gm += m * ((CAtom*)g_oaAtoms[z])->m_iCount;
	}
	m = (double)mypow(gm/6.02214086*10,0.33333333) * 100.0;  // In pm
	BTOUT;
	return m;
}


void strtolower(char *s)
{
	BTIN;
	char *p;
	p = s;
	while (*p != 0)
	{
		if ((*p >= 'A') && (*p <= 'Z'))
			*p += 32;
		p++;
	}
	BTOUT;
}


void SortAtoms()
{
	BTIN;
	int z, z2;
	int i;
	CAtom *a;

	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		i = z;
		for (z2=z+1;z2<g_oaAtoms.GetSize();z2++)
			if (strcmp(((CAtom*)g_oaAtoms[i])->m_sName,((CAtom*)g_oaAtoms[z2])->m_sName) > 0)
				i = z2;
		if (i != z)
		{
			a = (CAtom*)g_oaAtoms[z];
			g_oaAtoms[z] = g_oaAtoms[i];
			g_oaAtoms[i] = a;
		}
	}
	for (z=0;z<g_oaAtoms.GetSize();z++)
		((CAtom*)g_oaAtoms[z])->m_iIndex = z;
	BTOUT;
}


void SortElementsLabel()
{
	BTIN;
	int z, z2;
	int i;
	CElement *e;

	for (z=0;z<g_oaElements.GetSize();z++)
	{
		i = z;
		for (z2=z+1;z2<g_oaElements.GetSize();z2++)
			if (strcmp(((CElement*)g_oaElements[i])->m_sLabel,((CElement*)g_oaElements[z2])->m_sLabel) > 0)
				i = z2;
		if (i != z)
		{
			e = (CElement*)g_oaElements[z];
			g_oaElements[z] = g_oaElements[i];
			g_oaElements[i] = e;
		}
	}
	BTOUT;
}


void SortElementsMass()
{
	BTIN;
	int z, z2;
	int i;
	CElement *e;

	for (z=0;z<g_oaElements.GetSize();z++)
	{
		i = z;
		for (z2=z+1;z2<g_oaElements.GetSize();z2++)
			if (((CElement*)g_oaElements[i])->m_fMass > ((CElement*)g_oaElements[z2])->m_fMass)
				i = z2;
		if (i != z)
		{
			e = (CElement*)g_oaElements[z];
			g_oaElements[z] = g_oaElements[i];
			g_oaElements[i] = e;
		}
	}
	BTOUT;
}


bool SetAnalysis(const char *s)
{
	BTIN;

	if (mystricmp(s, "fpp") == 0) {
		g_bFixedPlProj = true;
		return true;
	}
	if (mystricmp(s, "chdf") == 0) {
		g_bCHDF = true;
		return true;
	}
	if(mystricmp(s, "chi") == 0) {
		g_bChiral = true;
		return true;
	}
	if(mystricmp(s, "eck") == 0) {
		g_bEckartTransform = true;
		return true;
	}
	if (mystricmp(s, "pol") == 0) {
		g_bSetUpPolarizabilityCalc = true;
		return true;
	}
	if(mystricmp(s, "swan") == 0) {
		g_bSortWannier = true;
		return true;
	}
	if (mystricmp(s, "order") == 0) {
		g_bOrder = true;
		return true;
	}
	if (mystricmp(s, "agtopo") == 0) {
		g_bAggrTopo = true;
		return true;
	}
	if (mystricmp(s, "cmat") == 0) {
		g_bContactMatrix = true;
		return true;
	}
	if (mystricmp(s, "geodens") == 0) {
		g_bGeoDens = true;
		return true;
	}
	if (mystricmp(s, "tddf") == 0) {
		g_bTDDF = true;
		return true;
	}
	if (mystricmp(s, "sankey") == 0) {
		mprintf("\n");
		mprintf("    To create a Sankey diagram, please perform an aggregation topology analysis first.\n");
		mprintf("    Choose \"agtopo\" in the main function menu to do so. As a result, you will obtain\n");
		mprintf("    an \"aggrtopo.dat\" file. If you have that file, start TRAVIS again with the parameters\n\n");
		mprintf("    \"travis -sankey aggrtopo.dat\".\n\n");
		mprintf("    Leaving now.\n\n");
		exit(0);
	}

	if (mystricmp(s, "rera") == 0) {
		g_bReRa = true;
		return true;
	}
	if(mystricmp(s, "vcd") == 0) {
		g_bVCD = true;
		return true;
	}
	if (mystricmp(s, "spec") == 0) {
		g_bROA = true;
		return true;
	}
	if (mystricmp(s, "drst") == 0) {
		g_bDipoleRestart = true;
		return true;
	}
	if (mystricmp(s, "mrst") == 0) {
		g_bMagneticDipoleRestart = true;
		return true;
	}

	if(mystricmp(s, "nc") == 0) {
		g_bNormalCoordinate = true;
		return true;
	}
	if (mystricmp(s,"dprof")==0)
	{
/*		if ((!g_bPeriodicX) || (!g_bPeriodicY) || (!g_bPeriodicZ))
		{
			eprintf("\n    Error: Requires XYZ-periodic box!\n\n");
			BTOUT;
			return false;
		}*/
		g_bPDF = true;
		return true;
	}

	if (mystricmp(s,"vori")==0)
	{
		g_bTegri = true;
		BTOUT;
		return true;
	}
	
	if (mystricmp(s,"doma")==0)
	{
		g_bDomA = true;
		BTOUT;
		return true;
	}

	if (mystricmp(s,"plproj")==0)
	{
		g_bPlProj = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"pldf")==0)
	{
		g_bPlDF = true;
		BTOUT;
		return true;
	}


	if (mystricmp(s,"lidf")==0)
	{
		g_bLiDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"power")==0)
	{
// 		g_bACF = true;
// 		g_bPowerSpec = true;
		g_bPower = true;
		BTOUT;
		return true;
	}


	if (mystricmp(s,"ir")==0)
	{
// 		g_bRDyn = true;
// 		g_bIRSpec = true;
		g_bIR = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"reg")==0)
	{
		g_bRegionAnalysis = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"raman")==0)
	{
		g_bRaman = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"sfac")==0)
	{
		g_bSFac = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"dens")==0)
	{
		if ((!g_bPeriodicX) || (!g_bPeriodicY) || (!g_bPeriodicZ))
		{
			eprintf("\n    Error: Requires XYZ-periodic box!\n\n");
			BTOUT;
			return false;
		}
		g_bDens = true;
		BTOUT;
		return true;
	}


	if (mystricmp(s,"cond")==0)
	{
		g_bCond = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"sdf")==0)
	{
		g_bSDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"psdf")==0)
	{
		g_bRevSDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"rdf")==0)
	{
		g_bRDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"vdf")==0)
	{
		g_bVDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"fdf")==0)
	{
		g_bFDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"adf")==0)
	{
		g_bADF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"ddf")==0)
	{
		g_bDDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"cdf")==0)
	{
		g_bCDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"msd")==0)
	{
		g_bMSD = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"aggr")==0)
	{
		g_bAggregation = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"hbond") == 0) 
	{
		g_bHBond = true;
		return true;
	}
	
	if (mystricmp(s,"ionpair") == 0) 
	{
		g_bIonPair = true;
		return true;
	}

	if (mystricmp(s,"hole") == 0) 
	{
	        g_bHole = true;
	        return true;
	}

	
	if (mystricmp(s,"nbex")==0)
	{
		g_bNbExchange = true;
		BTOUT;
		return true;
	}
/*	if (mystricmp(s,"ddisp")==0)
	{
		g_bDDisp = true;
		g_bAggregation = true;
		BTOUT;
		return true;
	}*/
	if (mystricmp(s,"acf")==0)
	{
		//g_bACF = true;
		g_bVACFNew = true;
		BTOUT;
		return true;
	}

	if (mystricmp(s,"cla")==0)
	{
		g_bClusterAnalysis = true;
		BTOUT;
		return true;
	}


	if (mystricmp(s,"nbx")==0)
	{
		g_bNbExchange = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"rdyn")==0)
	{
		g_bRDyn = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"cut")==0)
	{
		g_bCutCluster = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"avg")==0)
	{
		g_bAvg = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"env")==0)
	{
		g_bSaveRefEnv = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"proc")==0)
	{
		g_bSaveJustTraj = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"vhcf")==0)
	{
		g_bVHDF = true;
		BTOUT;
		return true;
	}
	if (mystricmp(s,"nbh")==0)
	{
		g_bNbAnalysis = true;
		BTOUT;
		return true;
	}
/*	if (mystricmp(s,"vfd")==0)
	{
		g_bSaveVelForce = true;
		BTOUT;
		return true;
	}*/
	if (mystricmp(s,"dip")==0)
	{
		g_bDipDF = true;
		g_bDipole = true;
		BTOUT;
		return true;
	}

	if (mystricmp(s,"voro")==0)
	{
		if ((!g_bPeriodicX) || (!g_bPeriodicY) || (!g_bPeriodicZ))
		{
			eprintf("\n    Error: Basic Voronoi Analysis currently requires XYZ-periodic box.\n\n");
			BTOUT;
			return false;
		}
		g_bVoro = true;
		BTOUT;
		return true;
	}

	if (mystricmp(s,"vorochg")==0)
	{
		#ifdef NEW_CHARGEVAR
			g_bMinimizeChargeVar2 = true;
		#else
			eprintf("Error: TRAVIS was not compiled with support for charge variance minimization.\n");
			eprintf("       You need to switch on the NEW_CHARGEVAR preprocessor definition in \"config.h\"\n");
			eprintf("       and supply the libtravis subdirectory.\n\n");
			return false;
		#endif
		return true;
	}


	BTOUT;
	return false;
}


bool ParseFunctions(const char *s)
{
	BTIN;
	const char *cp;
	char *p, *q;
	char buf[256];
	int z;
	CAnalysis *a;



	g_bReRa = false;
	g_bMinimizeChargeVar2 = false;
	g_bAmberCellAsked = false;
	g_bAmberCellInfo = false;
	g_bTDDF = false;
	g_bGeoDens = false;
	g_bContactMatrix = false;
	g_bAggrTopo = false;
	g_bOrder = false;
	g_bFixedPlProj = false;
	g_bROA = false;
	g_bSFac = false;
	g_bVoid = false;
	g_bTegri = false;
	g_bVoro = false;
	g_bDomA = false;
	g_bPlProj = false;
	g_bPlDF = false;
	g_bLiDF = false;
	g_bDens = false;
	g_bCond = false;
	g_bSDF = false;
	g_bRDF = false;
	g_bVDF = false;
	g_bFDF = false;
	g_bADF = false;
	g_bAvg = false;
	g_bSaveVelForce = false;
	g_bSaveRefEnv = false;
	g_bRefEnvCenter = false;
	g_bSaveJustTraj = false;
	g_bVHDF = false;
	g_bCutCluster = false;
	g_bNbAnalysis = false;
	g_bVACF = false;
	g_bDipACF = false;
	g_bDDF = false;
	g_bMSD = false;
	g_bDipole = false;
	g_bDipDF = false;
	g_bCDF = false;
	g_bAggregation = false;
	g_bDDisp = false; 
	g_bDACF = false; 
	g_bDLDF = false; 
	g_bRDyn = false;
	g_bNbExchange = false;
	g_bRaman = false;
	g_bIRSpec = false;
	g_bPowerSpec = false;
	g_bRegionAnalysis = false;
	g_bVACFNew = false;

	cp = s;
	q = buf;

	while (*cp != 0)
	{
		if (*cp == ' ')
			cp++;
		else {
			*q = *cp;
			cp++;
			q++;
		}
	}
	*q = 0;

	p = buf;
	
	while (true)
	{
		q = strchr(p,',');
		if (q != NULL)
			*q = 0;

		for (z=0;z<g_oaAnalyses.GetSize();z++)
		{
			a = (CAnalysis*)g_oaAnalyses[z];
			if (mystricmp(p,a->m_sAbbrev)==0)
			{
				mprintf(WHITE," %-8s",a->m_sAbbrev);
				mprintf("- %s\n",a->m_sName);
				if (!SetAnalysis(p))
				{
					eprintf("Function cannot be applied to this system: \"%s\"\n",p);
					return false;
				}
				goto _done;
			}
		}
		eprintf("Unknown Function: \"%s\"\n",p);
		return false;
_done:
		if ((q == NULL) || (*(q+1) == 0))
		{
			BTOUT;
			return true;
		}
		p = q+1;
	}
}


bool ParsePeriodic(const char *s)
{
	const char *p;

	BTIN;
	if (strlen(s) == 0)
	{
		g_bPeriodicX = true;
		g_bPeriodicY = true;
		g_bPeriodicZ = true;
		g_bPeriodic = true;
		BTOUT; 
		return true;
	}
	g_bPeriodicX = false;
	g_bPeriodicY = false;
	g_bPeriodicZ = false;

	p = s;
	while (*p != 0)
	{
		if ((*p == 'x') || (*p == 'X'))
			g_bPeriodicX = true;
		else if ((*p == 'y') || (*p == 'Y'))
			g_bPeriodicY = true;
		else if ((*p == 'z') || (*p == 'Z'))
			g_bPeriodicZ = true;
		else if (*p != '0')
		{
			eprintf("Invalid input.\n");
			BTOUT;
			return false;
		}
		p++;
	}
	g_bPeriodic = g_bPeriodicX || g_bPeriodicY || g_bPeriodicZ;
	BTOUT;
	return true;
}


void WriteHeader()
{
	struct tm *today;
	time_t ltime;
	char buf[64];
	unsigned long l, *lp;
	unsigned char *ca, *cb, *cc, *cd;
	bool b;

	BTIN;
	mprintf("\n");  /* http://patorjk.com/software/taag/ "Big Money-SW" */

/*	mprintf(" ________          ______   __     __  __           \n");
	mprintf("/        |        /      \\ /  |   /  |/  |          \n");
	mprintf("########/______  /######  |## |   ## |##/   _______ \n");
	mprintf("   ## | /      \\ ## |__## |## |   ## |/  | /       |\n");
	mprintf("   ## |/######  |##    ## |##  \\ /##/ ## |/#######/ \n");
	mprintf("   ## |## |  ##/ ######## | ##  /##/  ## |##      \\ \n");
	mprintf("   ## |## |      ## |  ## |  ## ##/   ## | ######  |\n");
	mprintf("   ## |## |      ## |  ## |   ###/    ## |/     ##/ \n");
	mprintf("   ##/ ##/       ##/   ##/     #/     ##/ #######/  \n");*/

/*	mprintf("   ________                             __\n");
	mprintf("  /        |                           /  |\n");
	mprintf("  ########/______   ______   __     __ ##/   _______\n");
	mprintf("     ## | /      \\ /      \\ /  \\   /  |/  | /       |\n");
	mprintf("     ## |/######  |######  |##  \\ /##/ ## |/#######/\n");
	mprintf("     ## |## |  ##/ /    ## | ##  /##/  ## |##      \\\n");
	mprintf("     ## |## |     /####### |  ## ##/   ## | ######  |\n");
	mprintf("     ## |## |     ##    ## |   ###/    ## |/     ##/\n");
	mprintf("     ##/ ##/       #######/     #/     ##/ #######/\n");*/

	mprintf(YELLOW,"  ________                                 __\n");
	mprintf(YELLOW," /        |                               /  |\n");
	mprintf(YELLOW," ########/ ______    ______    __     __  ##/    _______\n");
	mprintf(YELLOW,"    ## |  /      \\  /      \\  /  \\   /  | /  |  /       |\n");
	mprintf(YELLOW,"    ## | /######  | ######  | ##  \\ /##/  ## | /#######/\n");
	mprintf(YELLOW,"    ## | ## |  ##/  /    ## |  ##  /##/   ## | ##      \\\n");
	mprintf(YELLOW,"    ## | ## |      /####### |   ## ##/    ## |  ######  |\n");
	mprintf(YELLOW,"    ## | ## |      ##    ## |    ###/     ## | /     ##/\n");
	mprintf(YELLOW,"    ##/  ##/        #######/      #/      ##/  #######/\n");

	mprintf(WHITE,"\n    TRajectory Analyzer and VISualizer  -  Open-source free software under GNU GPL v3\n\n");
//	mprintf("");
	mprintf("    Copyright (c) Martin Brehm      (2009-2022), University of Halle (Saale)\n");
	mprintf("                  Martin Thomas     (2012-2022)\n");
	mprintf("                  Sascha Gehrke     (2016-2022), University of Bonn\n");
	mprintf("                  Barbara Kirchner  (2009-2022), University of Bonn\n");
	mprintf("\n");
	mprintf(YELLOW,"    http://www.travis-analyzer.de\n\n");
//	mprintf("     Open-source freeware; Licensed under the GNU General Public License v3.\n\n");
	mprintf("    Please cite: ");
	mprintf(WHITE,"J. Chem. Phys. 2020, 152 (16), 164105.         ");
	mprintf("(DOI 10.1063/5.0005078 )\n");
	mprintf(WHITE,"                 J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  ");
	mprintf("(DOI 10.1021/ci200217w )\n\n");
	mprintf("    There is absolutely no warranty on any results obtained from TRAVIS.\n\n");


	time(&ltime);
	today = localtime(&ltime);
	strcpy(buf,asctime(today));
	buf[strlen(buf)-1] = 0;
	mprintf(WHITE," #  ");
	if (g_sHostName != NULL)
		mprintf("Running on %s at %s",g_sHostName,buf);
			else mprintf("Running at %s",buf);

	#ifdef TARGET_LINUX
		mprintf(" (PID %d)\n",getpid());
	#elif defined(TARGET_WINDOWS)
		mprintf(" (PID %d)\n",GetCurrentProcessId());
	#else
		mprintf("\n");
	#endif

	if (g_sWorkingDir != NULL)
	{
		mprintf(WHITE," #  ");
		mprintf("Running in %s\n",g_sWorkingDir);
	}

/*	mprintf(WHITE,"  #  ");
	mprintf("Source code version: ");
	mprintf("%s.\n",SOURCE_VERSION);

	#ifdef RELEASE_VERSION
		mprintf(WHITE,"  #  ");
		mprintf("Release version: ");
		mprintf("%s.\n",RELEASE_VERSION);
	#endif

	mprintf(WHITE,"  #  ");
	mprintf("Compiled at ");
	mprintf("%s %s.\n",__DATE__,__TIME__);*/

	mprintf(WHITE," #  ");
	mprintf("Version: ");
	mprintf("%s",SOURCE_VERSION);

	#ifdef RELEASE_VERSION
		mprintf(",  Release ");
		mprintf("%s",RELEASE_VERSION);
	#endif

	mprintf(", built at ");
	mprintf("%s, %s",__DATE__,__TIME__);

	#if defined(__VERSION) || defined(__GNUC__) || defined(_MSC_VER) || defined(__clang__) || defined(__llvm__) || defined(__INTEL_COMPILER)
		mprintf(", compiler");
		b = false;
		#ifdef __VERSION__
			mprintf(" \"%s\"",__VERSION__);
			b = true;
		#endif
		#ifdef __llvm__
			if (b)
				mprintf(",");
			mprintf(" LLVM");
			b = true;
		#endif
		#ifdef __INTEL_COMPILER
				mprintf(" Intel Compiler \"%s\"", TOSTRING(__INTEL_COMPILER) );
		#elif defined( __clang__ )
			if (b)
				mprintf(",");
			#if defined(__clang_major__) && defined(__clang_minor__) && defined(__clang_patchlevel__)
				mprintf(" CLANG %d.%d.%d",__clang_major__,__clang_minor__,__clang_patchlevel__);
			#elif defined(__clang_version__)
				mprintf(" CLANG \"%s\"",__clang_version__);
			#else
				mprintf(" CLANG");
			#endif
			b = true;
		#elif defined( __GNUC__ )
			if (b)
				mprintf(",");
			#ifdef __MINGW32__
				mprintf(" MinGW GCC %d.%d.%d",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__);
			#else
				mprintf(" GCC %d.%d.%d",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__);
			#endif
			b = true;
		#endif
		#ifdef _MSC_VER
			if (b)
				mprintf(",");
			mprintf(" MSVC++ ");
			// See  https://docs.microsoft.com/de-de/cpp/preprocessor/predefined-macros?view=vs-2019
			switch(_MSC_VER) {
				case 1932: mprintf("17.2, VS 2022"); break;
				case 1931: mprintf("17.1, VS 2022"); break;
				case 1930: mprintf("17.0, VS 2022"); break;
				case 1929: mprintf("16.10/11, VS 2019"); break;
				case 1928: mprintf("16.8/9, VS 2019"); break;
				case 1927: mprintf("16.7, VS 2019"); break;
				case 1926: mprintf("16.6, VS 2019"); break;
				case 1925: mprintf("16.5, VS 2019"); break;
				case 1924: mprintf("16.4, VS 2019"); break;
				case 1923: mprintf("16.3, VS 2019"); break;
				case 1922: mprintf("16.2, VS 2019"); break;
				case 1921: mprintf("16.1, VS 2019"); break;
				case 1920: mprintf("16.0, VS 2019"); break;
				case 1916: mprintf("15.9, VS 2017"); break;
				case 1915: mprintf("15.8, VS 2017"); break;
				case 1914: mprintf("15.7, VS 2017"); break;
				case 1913: mprintf("15.6, VS 2017"); break;
				case 1912: mprintf("15.5, VS 2017"); break;
				case 1911: mprintf("15.3, VS 2017"); break;
				case 1910: mprintf("15.1, VS 2017"); break;
				case 1900: mprintf("14.0, VS 2015"); break;
				case 1800: mprintf("12.0, VS 2013"); break;
				case 1700: mprintf("11.0, VS 2012"); break;
				case 1600: mprintf("10.0, VS 2010"); break;
				case 1500: mprintf("9.0, VS 2008"); break;
				case 1400: mprintf("8.0, VS 2005"); break;
				case 1310: mprintf("7.1, VS 2003"); break;
				case 1300: mprintf("7.0"); break;
				case 1200: mprintf("6.0"); break;
				case 1100: mprintf("5.0"); break;
				default: mprintf("%d",_MSC_VER); break;
			}
			b = true;
		#endif
	#endif

	mprintf("\n");

	#ifdef TARGET_WINDOWS
		mprintf(WHITE," #  ");
		mprintf("Target platform: Windows");
	#elif defined(TARGET_LINUX)
		mprintf(WHITE," #  ");
		mprintf("Target platform: Linux");
	#else
		mprintf(WHITE," #  ");
		mprintf("Target platform: Generic");
	#endif

	#define STRINGIFY(x) #x
	#define TOSTRING(x) STRINGIFY(x)

	#ifdef __cplusplus
		mprintf( ", __cplusplus=%s", TOSTRING(__cplusplus) );
	#endif

	mprintf(", Compile flags: ");

	#ifdef NEW_CHARGEVAR
		mprintf("NEW_CHARGEVAR ");
	#endif

	#ifdef USE_OMP
		mprintf("USE_OMP ");
	#endif

	#ifdef DEBUG_BACKTRACE
		mprintf("DEBUG_BACKTRACE ");
	#endif

	#ifdef DEBUG_EXTENDED_BACKTRACE
		mprintf("DEBUG_EXTENDED_BACKTRACE ");
	#endif

	#ifdef USE_FFTW
		mprintf("USE_FFTW ");
	#endif

	#ifdef DEBUG_ARRAYS
		mprintf("DEBUG_ARRAYS ");
	#endif

	#ifdef DEBUG_COBARRAY
		mprintf("DEBUG_COBARRAY ");
	#endif

	#ifdef DEBUG_CSTRING
		mprintf("DEBUG_CSTRING ");
	#endif

	#ifdef DEBUG_CVEC3ARRAY
		mprintf("DEBUG_CVEC3ARRAY ");
	#endif

	#ifdef DEBUG_DATABASE
		mprintf("DEBUG_DATABASE ");
	#endif

	mprintf("\n");


	#ifdef _OPENMP
		#ifdef USE_OMP
			mprintf(WHITE," #  ");
			mprintf( "OpenMP is active, max_threads=%d\n", omp_get_max_threads() );
		#else
			mprintf(WHITE," #  ");
			mprintf( RED, "Compiled with OpenMP, but USE_OMP not switched on in \"config.h\"!\n" );
		#endif
	#else
		#ifdef USE_OMP
			mprintf(WHITE," #  ");
			mprintf( RED, "USE_OMP is switched on, but code not compiled with -fopenmp!\n" );
		#endif
	#endif


	l = 0xA0B0C0D0;
	lp = &l;
	ca = (unsigned char*)lp;
	cb = ca+1;
	cc = ca+2;
	cd = ca+3;
	mprintf(WHITE," #  ");
	mprintf("Machine: ");
	
	// Source:  https://stackoverflow.com/questions/152016/detecting-cpu-architecture-compile-time
    #if defined(__x86_64__) || defined(_M_X64)
		mprintf( "x86_64" );
    #elif defined(i386) || defined(__i386__) || defined(__i386) || defined(_M_IX86)
		mprintf( "i386" );
    #elif defined(__ARM_ARCH_2__)
		mprintf( "ARM2" );
    #elif defined(__ARM_ARCH_3__) || defined(__ARM_ARCH_3M__)
		mprintf( "ARM3" );
    #elif defined(__ARM_ARCH_4T__) || defined(__TARGET_ARM_4T)
		mprintf( "ARM4T" );
    #elif defined(__ARM_ARCH_5_) || defined(__ARM_ARCH_5E_)
		mprintf( "ARM5"
    #elif defined(__ARM_ARCH_6T2_) || defined(__ARM_ARCH_6T2_)
		mprintf( "ARM6T2" );
    #elif defined(__ARM_ARCH_6__) || defined(__ARM_ARCH_6J__) || defined(__ARM_ARCH_6K__) || defined(__ARM_ARCH_6Z__) || defined(__ARM_ARCH_6ZK__)
		mprintf( "ARM6" );
    #elif defined(__ARM_ARCH_7__) || defined(__ARM_ARCH_7A__) || defined(__ARM_ARCH_7R__) || defined(__ARM_ARCH_7M__) || defined(__ARM_ARCH_7S__)
		mprintf( "ARM7" );
    #elif defined(__ARM_ARCH_7A__) || defined(__ARM_ARCH_7R__) || defined(__ARM_ARCH_7M__) || defined(__ARM_ARCH_7S__)
		mprintf( "ARM7A" );
    #elif defined(__ARM_ARCH_7R__) || defined(__ARM_ARCH_7M__) || defined(__ARM_ARCH_7S__)
		mprintf( "ARM7R" );
    #elif defined(__ARM_ARCH_7M__)
		mprintf( "ARM7M" );
    #elif defined(__ARM_ARCH_7S__)
		mprintf( "ARM7S" );
    #elif defined(__aarch64__) || defined(_M_ARM64)
		mprintf( "ARM64" );
    #elif defined(mips) || defined(__mips__) || defined(__mips)
		mprintf( "MIPS" );
    #elif defined(__sh__)
		mprintf( "SuperH" );
    #elif defined(__powerpc) || defined(__powerpc__) || defined(__powerpc64__) || defined(__POWERPC__) || defined(__ppc__) || defined(__PPC__) || defined(_ARCH_PPC)
		mprintf( "PowerPC" );
    #elif defined(__PPC64__) || defined(__ppc64__) || defined(_ARCH_PPC64)
		mprintf( "PowerPC64" );
    #elif defined(__sparc__) || defined(__sparc)
		mprintf( "Sparc" );
    #elif defined(__m68k__)
		mprintf( "M68K" );
    #else
		mprintf( "Unknown" );
    #endif

	#ifdef __CHAR_UNSIGNED__
		mprintf( ", __CHAR_UNSIGNED__" );
	#endif

	mprintf(", char=%lub, int=%lub, long=%lub, addr=%lub, 0xA0B0C0D0=%02X,%02X,%02X,%02X.\n",(unsigned long)sizeof(char),(unsigned long)sizeof(int),(unsigned long)sizeof(long),(unsigned long)sizeof(void*),*ca,*cb,*cc,*cd);

	mprintf(WHITE," #  ");
	if (g_sHomeDir != NULL)
		mprintf("Home: %s,  ",g_sHomeDir);
	mprintf("Executable: %s\n",g_sExeName);

	mprintf(WHITE," #  ");
	if ((!IsTTY(stdin)) || (g_sInputFile != NULL))
	{
		if (g_sInputFile != NULL)
			mprintf("Input from %s, ",g_sInputFile);
				else mprintf("Input is redirected,  ");
	} else
		mprintf("Input from terminal,  ");

	if (IsTTY(stdout))
		mprintf("Output to terminal\n");
			else mprintf("Output is redirected\n");

	if (g_bStreamInput)
	{
		mprintf(WHITE," #  ");
		mprintf("Input trajectory treated as stream\n");
	}

	mprintf("\n");
	mprintf("    >>> Please use a color scheme with dark background or specify \"-nocolor\"! <<<\n\n");
	BTOUT;
}


void CommandLineHelp()
{
	BTIN;
	mprintf(WHITE,"    List of supported command line options:\n\n");
	mprintf("      -p <file>       Load position data from specified trajectory file.\n");
	mprintf("                      Format may be *.xyz, *.pdb, *.lmp (LAMMPS), HISTORY (DLPOLY), POSCAR/XDATCAR (VASP),\n");
	mprintf("                                    *.gro, *.dcd, or *.prmtop/*.mdcrd (Amber).\n");
	mprintf("                      The bqb format (*.bqb, *.btr, *.emp, *.blist) as well as *.voronoi are also supported.\n");
	mprintf("      -vel <file>     Read atom velocities from a file in addition to the position trajectory.\n");
	mprintf("                      Currently, only .xyz format is supported for velocity data.\n");
	mprintf("      -i <file>       Read input from specified text file.\n");
	mprintf("      -c <file>       Read and execute control file (experimental).\n");
	mprintf("      cubetool        Execute the CubeTool for modifying Gaussian Cube files.\n");
	mprintf("      -sankey <file>  Create Sankey diagrams (file name is optional).\n");
	mprintf("      -ramanfrompola  Compute Raman spectra from existing polarizability time series.\n");
	mprintf("      (de-)compress   Start built-in bqbtool (compress trajectories to BQB format).\n");
	mprintf("      check           Check BQB file integrity.\n");
//	mprintf("      compare         Compare atom coordinates in BQB file to reference trajectory.\n");
//	mprintf("      split           Split BQB file into parts.\n");
//	mprintf("      merge           Merge several BQB files into one.\n");
	mprintf("\n");
	mprintf("      -config <file>  Load the specified configuration file.\n");
	mprintf("      -stream         Treat input trajectory as a stream (e.g. named pipe): No fseek, etc.\n");
	mprintf("      -showconf       Show a tree structure of the configuration file.\n");
	mprintf("      -writeconf      Write the default configuration file, including all defines values.\n\n");
	mprintf("      -verbose        Show detailed information about what's going on.\n");
	mprintf("      -nocolor        Execute TRAVIS in monochrome mode (suitable for white background).\n");
	mprintf("      -dimcolor       Use dim instead of bright colors.\n\n");
	mprintf("      -credits        Display a list of persons who contributed to TRAVIS.\n");
	mprintf("      -help, -?       Shows this help.\n");
	mprintf("\n");
	mprintf("    If only one argument is specified, it is assumed to be the name of a trajectory file.\n");
	mprintf("    If no argument is specified at all, TRAVIS asks for the trajectory file to open.\n");
	BTOUT;
}


bool ParseArgs(int argc, const char *argv[])
{
	BTIN;
	const char *p;
	char SModeFlag[64];
	int z;
//	bool namegiven;


	UnBase64((unsigned char *)SModeFlag,(const unsigned char *)"LXNjaGlzcw==",12);

	g_bAsciiArt = false;
	z = 1;

//	g_sInputTraj[0] = 0;

//	g_sInputVel[0] = 0;
	g_sInputVel = "";

//	g_sInputForce[0] = 0;
	g_sInputForce = "";

//	g_sInputCtrl[0] = 0;
	g_sInputCtrl = "";

//	namegiven = false;

	while (z < argc)
	{
		if ((memcmp(argv[z],"-?",2)==0) || (memcmp(argv[z],"--?",3)==0) || (memcmp(argv[z],"-help",5)==0) || (memcmp(argv[z],"--help",6)==0))
		{
			CommandLineHelp();
			BTOUT;
			return false;
		}

		if (mystricmp(argv[z],"-p")==0)
		{
//			mprintf("Dateiname: \"%s\"\n",argv[z+1]);
			z++;

			try { g_sInputTraj = new char[strlen(argv[z])+1]; } catch(...) { g_sInputTraj = NULL; }
			if (g_sInputTraj == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			strcpy(g_sInputTraj,argv[z]);
			goto _argend;
		}

		if (mystricmp(argv[z],"-vel")==0)
		{
			z++;
			g_sInputVel.strcpy(argv[z]);
			goto _argend;
		}

		if (mystricmp(argv[z],"-conf")==0)
		{
//			mprintf("Dateiname: \"%s\"\n",argv[z+1]);
			z++;

			try { g_sConfFile = new char[strlen(argv[z])+1]; } catch(...) { g_sConfFile = NULL; }
			if (g_sConfFile == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			strcpy(g_sConfFile,argv[z]);
			goto _argend;
		}

		if (mystricmp(argv[z],"-i")==0)
		{
//			mprintf("Dateiname: \"%s\"\n",argv[z+1]);
			z++;

/*			try { g_sInputFile = new char[strlen(argv[z])+1]; } catch(...) { g_sInputFile = NULL; }
			if (g_sInputFile == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(g_sInputFile,argv[z]);*/
			goto _argend;
		}

		if (mystricmp(argv[z],"-c")==0) {

			z++;
			goto _argend;
		}

		if (mystricmp(argv[z],"-react")==0) {

			goto _argend;
		}

		if (mystricmp(argv[z],"-verbose")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-credits")==0)
		{
			goto _argend;
		}


		if (mystricmp(argv[z],"-showconf")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-stream")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-writeconf")==0)
		{
			goto _argend;
		}

/*		if (mystricmp(argv[z],"-f")==0)
		{
			z++;
			strcpy(g_sInputForce,argv[z]);
			goto _argend;
		}*/

		if ((mystricmp(argv[z],"-nocolor")==0) || (mystricmp(argv[z],"--nocolor")==0))
		{
			goto _argend;
		}

		if ((mystricmp(argv[z],"-dimcolor")==0) || (mystricmp(argv[z],"--dimcolor")==0))
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-a")==0)
		{
//			mprintf("AA ist an.\n");
			g_bAsciiArt = true;
			goto _argend;
		}

		if (mystricmp(argv[z],"-lsd")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-readdipole")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],"-sax")==0)
		{
			goto _argend;
		}

		if (mystricmp(argv[z],SModeFlag)==0)
		{
			g_bSMode = true;
			goto _argend;
		}

		if ((argc > 2) || (argv[z][0] == '-'))
		{
			eprintf("Unknown parameter: \"%s\".\n\n",argv[z]);
			CommandLineHelp();
			BTOUT;
			return false;
		}

		p = strrchr(argv[z],'.');
		if (p != NULL)
		{
			p++;
			if ((mystricmp(p,"xyz")==0) || (mystricmp(p,"cube")==0) || (mystricmp(p,"pdb")==0) || (mystricmp(p,"mol2")==0) ||
				(mystricmp(p,"lmp")==0) || (mystricmp(p,"HISTORY")==0) || (mystricmp(p,"prmtop")==0) || (mystricmp(p,"mdcrd")==0) ||
				(mystricmp(p,"gro")==0) || (mystricmp(p,"dcd")==0) || (mystricmp(p,"voronoi")==0) || 
				(mystricmp(p,"bqb")==0) || (mystricmp(p,"bbq")==0) || (mystricmp(p,"btr")==0) || (mystricmp(p,"blist")==0) || 
				(mystricmp(p,"emp")==0))	{

				try { g_sInputTraj = new char[strlen(argv[z])+1]; } catch(...) { g_sInputTraj = NULL; }
				if (g_sInputTraj == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				strcpy(g_sInputTraj,argv[z]);

			} else
			{
				eprintf("Unknown parameter: \"%s\".\n\n",argv[z]);
				CommandLineHelp();
				BTOUT;
				return false;
			}

//			try { g_sInputTraj = new char[strlen(argv[z])+1]; } catch(...) { g_sInputTraj = NULL; }
//			if (g_sInputTraj == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

//			strcpy(g_sInputTraj,argv[z]);
		} else
		{

			if (mystricmp(&argv[z][strlen(argv[z])-7],"HISTORY")==0)
			{
				try { g_sInputTraj = new char[strlen(argv[z])+1]; } catch(...) { g_sInputTraj = NULL; }
				if (g_sInputTraj == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				strcpy(g_sInputTraj,argv[z]);
			} else
			{
				eprintf("Unknown parameter: \"%s\".\n\n",argv[z]);
				CommandLineHelp();
				BTOUT;
				return false;
			}

		}
		BTOUT;
		return true;

_argend:
		z++;
	}
	BTOUT;
	return true;
}


void ParsePassiveArgs(int argc, const char *argv[])
{
	BTIN;
	int z;

	g_bGlobalPsycho = false;

	z = 1;
	while (z < argc)
	{
		if (mystricmp(argv[z],"-credits")==0)
			g_bShowCredits = true;

		if (mystricmp(argv[z],"-lsd")==0)
			g_bGlobalPsycho = true;

		if (mystricmp(argv[z],"-readdipole")==0)
			g_bDipolGrimme = true;

		if (mystricmp(argv[z],"-sax")==0)
			g_bSaxonize = true;

		if (mystricmp(argv[z],"-verbose")==0)
			g_bVerbose = true;

		if ((mystricmp(argv[z],"-ramanfrompola")==0) || (mystricmp(argv[z],"-ramanfrompolarizability")==0))
			g_bRamanFromPolarizability = true;


		if (mystricmp(argv[z],"-sankey")==0) {
			g_bSankey = true;
			if (z+1 < argc) {
				strcpy( g_sSankeyFile, argv[z+1] );
				z++;
			} else
				g_sSankeyFile[0] = 0;
		}

		if (mystricmp(argv[z],"-showconf")==0)
			g_bShowConf = true;

		if (mystricmp(argv[z],"-stream")==0)
			g_bStreamInput = true;


		if (mystricmp(argv[z],"-vel")==0)
			g_bReadVelocity = true;

		if (mystricmp(argv[z],"-writeconf")==0)
			g_bWriteConf = true;

		if ((mystricmp(argv[z],"-nocolor")==0) || (mystricmp(argv[z],"--nocolor")==0))
			g_bNoColor = true;

		if ((mystricmp(argv[z],"-dimcolor")==0) || (mystricmp(argv[z],"--dimcolor")==0))
			g_iColorIntensity = 2;

		if (mystricmp(argv[z],"-i")==0)
		{
			z++;

			try { g_sInputFile = new char[strlen(argv[z])+1]; } catch(...) { g_sInputFile = NULL; }
			if (g_sInputFile == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(g_sInputFile,argv[z]);
		}

		if (mystricmp(argv[z],"-c")==0) {

			z++;

			try { g_sControlFile = new char[strlen(argv[z])+1]; } catch(...) { g_sControlFile = NULL; }
			if (g_sControlFile == NULL) NewException((double)(strlen(argv[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(g_sControlFile,argv[z]);

			g_bControlRun = true;
		}

		z++;
	}
	BTOUT;
}


void CreateDatabaseDefaults()
{
	g_pDatabase->AddString("/GLOBAL/TRAVIS_VERSION",SOURCE_VERSION);
	g_pDatabase->AddBool("/GLOBAL/SHOWCREDITS",false);

	g_pDatabase->AddBool("/PLOT2D/FORMATS/WRITE_MATHEMATICA",true);
	g_pDatabase->AddBool("/PLOT2D/FORMATS/WRITE_GNUPLOT",true);
	g_pDatabase->AddBool("/PLOT2D/FORMATS/WRITE_TRIPLES",true);
	g_pDatabase->AddBool("/PLOT2D/FORMATS/WRITE_MATRIX",true);

	g_pDatabase->AddInt("/PLOT2D/DEFAULTS/BIN_RES",100);
	g_pDatabase->AddInt("/PLOT2D/DEFAULTS/IMAGE_RES",1000);
	g_pDatabase->AddFloat("/PLOT2D/DEFAULTS/PLOT_EXP",0.5);
	g_pDatabase->AddInt("/PLOT2D/DEFAULTS/CONTOUR_LINES",30);

	g_pDatabase->AddBool("/PLOT3D/FORMATS/WRITE_CUBE",true);
	g_pDatabase->AddBool("/PLOT3D/FORMATS/WRITE_PLT",true);
	g_pDatabase->AddInt("/PLOT3D/DEFAULTS/BIN_RES",100);

/*	g_pDatabase->AddString("/BLA/BLUBB/PLOEPP/STRING1","String 1 Content");
	g_pDatabase->AddInt("/BLA/BLUBB/PLOEPP/INT1",123456);
	g_pDatabase->AddFloat("/BLA/BLUBB/PLOEPP/FLOAT1",1.23456);
	g_pDatabase->AddBool("/BLA/BLUBB/PLOEPP/BOOL1",true);*/
	//g_pDatabase->AddInt("/BLA/BLUBB/PLOEPP/STRING1",17);
}


void LoadSettings()
{
//	char buf[256];
	CxString buf;
	char sep;
	FILE *a;

#ifdef TARGET_WINDOWS
	sep = '\\';
#else
	sep = '/';
#endif

	if (g_sConfFile != NULL)
	{
		mprintf(WHITE,"    Custom configuration file specified: %s\n",g_sConfFile);
		if (FileExist(g_sConfFile))
		{
			mprintf("    Loading configuration from %s ...\n",g_sConfFile);
			g_pDatabase->ParseInputFile(g_sConfFile);
			return;
		} else
		{
			eprintf("LoadSettings(): Error: File does not exist or cannot be read.\n");
			abort();
		}
	}

//	sprintf(buf,"%s%ctravis.conf",g_sWorkingDir,sep);
	buf.sprintf("%s%ctravis.conf",g_sWorkingDir,sep);
	if (FileExist(buf))
	{
		mprintf("    Loading configuration from %s ...\n",(const char*)buf);
		g_pDatabase->ParseInputFile(buf);
		return;
	}

//	sprintf(buf,"%s%c.travis.conf",g_sWorkingDir,sep);
	buf.sprintf("%s%c.travis.conf",g_sWorkingDir,sep);
	if (FileExist(buf))
	{
		mprintf("    Loading configuration from %s ...\n",(const char*)buf);
		g_pDatabase->ParseInputFile(buf);
		return;
	}

//	sprintf(buf,"%s%ctravis.conf",g_sHomeDir,sep);
	buf.sprintf("%s%ctravis.conf",g_sHomeDir,sep);
	if (FileExist(buf))
	{
		mprintf("    Loading configuration from %s ...\n",(const char*)buf);
		g_pDatabase->ParseInputFile(buf);
		return;
	}

//	sprintf(buf,"%s%c.travis.conf",g_sHomeDir,sep);
	buf.sprintf("%s%c.travis.conf",g_sHomeDir,sep);
	if (FileExist(buf))
	{
		mprintf("    Loading configuration from %s ...\n",(const char*)buf);
		g_pDatabase->ParseInputFile(buf);
		return;
	}

	mprintf(WHITE,"    No configuration file found.\n");
	if (g_sHomeDir == NULL)
	{
		eprintf("\n    Could not detect user home directory, writing config file to current directory.\n\n");
//		sprintf(buf,".travis.conf");
		buf.sprintf(".travis.conf");
	} else
//		sprintf(buf,"%s%c.travis.conf",g_sHomeDir,sep);
		buf.sprintf("%s%c.travis.conf",g_sHomeDir,sep);
	mprintf("    Writing default configuration to %s ...\n",(const char*)buf);
	a = fopen(buf,"wt");
	if (a == NULL)
	{
		eprintf("    Cannot open %s for writing.\n",(const char*)buf);
		return;
	}
	fclose(a);

	g_pDatabase->WriteOutputFile(buf);
}


void InitDatabase()
{
//	char buf[256];
	CxString buf;
	char sep;
	FILE *a;

#ifdef TARGET_WINDOWS
	sep = '\\';
#else
	sep = '/';
#endif

	try { g_pDatabase = new CDatabase(); } catch(...) { g_pDatabase = NULL; }
	if (g_pDatabase == NULL) NewException((double)sizeof(CDatabase),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	CreateDatabaseDefaults();

	/******* Interface *********/
	Interface_DefaultConf();

	LoadSettings();
	
//	mprintf("\"%s\" vs \"%s\"\n",g_pDatabase->GetString("/GLOBAL/TRAVIS_VERSION"),SOURCE_VERSION);

	if (strcmp(g_pDatabase->GetString("/GLOBAL/TRAVIS_VERSION"),SOURCE_VERSION) != 0)
	{
		g_pDatabase->SetString("/GLOBAL/TRAVIS_VERSION",SOURCE_VERSION);
//		sprintf(buf,"%s%c.travis.conf",g_sHomeDir,sep);
		buf.sprintf("%s%c.travis.conf",g_sHomeDir,sep);
		mprintf(WHITE,"      TRAVIS version has changed.");
		if (!g_bWriteConf)
			mprintf(" You should write a new configuration file (specify \"-writeconf\").\n");
	}

	if (g_bWriteConf)
	{
		if (g_sHomeDir == NULL)
		{
			eprintf("\n    Could not detect user's home directory, writing config file to current directory.\n\n");
//			sprintf(buf,".travis.conf");
			buf.sprintf(".travis.conf");
		} else
//			sprintf(buf,"%s%c.travis.conf",g_sHomeDir,sep);
			buf.sprintf("%s%c.travis.conf",g_sHomeDir,sep);

		mprintf("    Writing new default configuration to %s ...\n",(const char*)buf);
		a = fopen(buf,"wt");
		if (a == NULL)
		{
			eprintf("    Error: Cannot open %s for writing.\n",(const char*)buf);
			goto _writeend;
		}
		fclose(a);
		g_pDatabase->WriteOutputFile(buf);
	}
_writeend:

	if (g_bShowConf)
	{
		mprintf(WHITE,"\n  Output of Database Tree:\n\n");
		g_pDatabase->DumpTree();
	}

	if (g_pDatabase->GetBool("/GLOBAL/SHOWCREDITS"))
		g_bShowCredits = true;

	mprintf("\n");
}


void RECURSION_BuildCDF(CObservation *o, int channel, int om, CxDoubleArray **data, double *result)
{
	BXIN;
	int z/*, z2*/;

	if (channel == g_iCDFChannels)
	{
		o->m_pCDF->AddToBin(result);
		BXOUT;
		return;
	}
/*	if (o->m_pCDF->m_bChannelAll[channel])
	{
		for (z=0;z<o->m_iShowMolCount;z++)
			for (z2=0;z2<data[channel][z].GetSize();z2++)
			{
				result[channel] = data[channel][z][z2];
				RECURSION_BuildCDF(o,channel+1,om,data,result);
			}
	} else*/
	{
		for (z=0;z<data[channel][om].GetSize();z++)
		{
			result[channel] = data[channel][om][z];
			RECURSION_BuildCDF(o,channel+1,om,data,result);
		}
	}
	BXOUT;
}


CVirtualAtom* AddVirtualAtom(int mol)
{
	BTIN;
	int z2;
	CMolecule *m;
	CVirtualAtom *va;
	CSingleMolecule *sm;
	CxIntArray *la;

	m = (CMolecule*)g_oaMolecules[mol];

	try { va = new CVirtualAtom(); } catch(...) { va = NULL; }
	if (va == NULL) NewException((double)sizeof(CVirtualAtom),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	g_oaVirtualAtoms.Add(va);
	//va->m_iIndex = g_iGesVirtAtomCount;
	va->m_iMolecule = (unsigned short)mol;
	va->m_iMolVirtAtom = (unsigned char)m->m_laVirtualAtoms.GetSize();
	sprintf(va->m_sLabel,"#%d",va->m_iMolVirtAtom+1);
	if (m->m_laVirtualAtoms.GetSize()==0)
	{
		m->m_baAtomIndex.Add(g_iVirtAtomType);
		m->m_waAtomCount.Add(1);
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			sm->m_baAtomIndex.Add(g_iVirtAtomType);

			try { la = new CxIntArray("AddVirtualAtom():sm->m_oaAtomOffset[]"); } catch(...) { la = NULL; }
			if (la == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			sm->m_oaAtomOffset.Add(la);
		}
	} else
	{
		m->m_waAtomCount[m->m_baAtomIndex.GetSize()-1]++;
//		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
//		{
//			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
//			sm->m_iAtomCount[sm->m_iElements-1]++;
//		}
	}
	for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
	{
		sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
		((CxIntArray*)sm->m_oaAtomOffset[sm->m_baAtomIndex.GetSize()-1])->Add(g_iGesVirtAtomCount);
		g_iGesVirtAtomCount++;
	}
	m->m_iAtomGes++;
	m->m_laVirtualAtoms.Add((unsigned short)g_oaVirtualAtoms.GetSize()-1);
	BTOUT;
	return va;
}


void AddElement(const char *s, int ord, double mass, double radius, double vdw)
{
	BTIN;
	CElement *e;

	try { e = new CElement(); } catch(...) { e = NULL; }
	if (e == NULL) NewException((double)sizeof(CElement),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { e->m_sLabel = new char[strlen(s)+1]; } catch(...) { e->m_sLabel = NULL; }
	if (e->m_sLabel == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(e->m_sLabel,s);
	e->m_iOrd = ord;
	e->m_fMass = mass;
	e->m_fRadius = radius;
	e->m_fVdWRadius = vdw;
	e->m_fCoherentNCS = 0;
	g_oaElements.Add(e);
	BTOUT;
}


void AddElement(const char *s, int ord, double mass, double radius, double vdw, double ncs)
{
	BTIN;
	CElement *e;

	try { e = new CElement(); } catch(...) { e = NULL; }
	if (e == NULL) NewException((double)sizeof(CElement),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { e->m_sLabel = new char[strlen(s)+1]; } catch(...) { e->m_sLabel = NULL; }
	if (e->m_sLabel == NULL) NewException((double)(strlen(s)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(e->m_sLabel,s);
	e->m_iOrd = ord;
	e->m_fMass = mass;
	e->m_fRadius = radius;
	e->m_fVdWRadius = vdw;
	e->m_fCoherentNCS = ncs;
	g_oaElements.Add(e);
	BTOUT;
}


void SetElementColor(const char *s, unsigned char r, unsigned char g, unsigned char b, unsigned char br, unsigned char bb, unsigned char bg)
{
	BTIN;
	CElement *e;
	e = FindElement(s,true);
	if (e == NULL)
	{
		eprintf("SetElementColor(): Element \"%s\" not found.\n",s);
		return;
	}
	e->m_iColorR = r;
	e->m_iColorG = g;
	e->m_iColorB = b;
	e->m_iColorBleachedR = br;
	e->m_iColorBleachedG = bg;
	e->m_iColorBleachedB = bb;
	BTOUT;
}


void RemoveAllElements()
{
	BTIN;
	int z;
	for (z=0;z<g_oaElements.GetSize();z++)
		delete (CElement*)g_oaElements[z];
	g_oaElements.RemoveAll();
	BTOUT;
}


void RemoveAllAtoms()
{
	BTIN;
	int z;
	for (z=0;z<g_oaAtoms.GetSize();z++)
		delete (CAtom*)g_oaAtoms[z];
	g_oaAtoms.RemoveAll();
	BTOUT;
}


void RemoveAllAnalyses()
{
	BTIN;
	int z;
	for (z=0;z<g_oaAnalyses.GetSize();z++)
		delete (CAnalysis*)g_oaAnalyses[z];
	g_oaAnalyses.RemoveAll();
	for (z=0;z<g_oaAnalysisGroups.GetSize();z++)
		delete (CAnalysisGroup*)g_oaAnalysisGroups[z];
	g_oaAnalysisGroups.RemoveAll();
	BTOUT;
}


void RemoveAllMolecules()
{
	BTIN;
	int z;

	for (z=0;z<g_oaMolecules.GetSize();z++)
		delete (CMolecule*)g_oaMolecules[z];
	g_oaMolecules.RemoveAll();
	for (z=0;z<g_oaSingleMolecules.GetSize();z++)
		delete (CSingleMolecule*)g_oaSingleMolecules[z];
	g_oaSingleMolecules.RemoveAll();
	BTOUT;
}


void RemoveAllObservations()
{
	BTIN;
	int z;

	for (z=0;z<g_oaObserv.GetSize();z++)
		delete (CObservation*)g_oaObserv[z];
	g_oaObserv.RemoveAll();
	BTOUT;
}


void GetTravisPath()
{
	BTIN;
	char *p, *q, *tmp, *env;

	if (FileExist(g_sExeName))
		return;
	if ((strchr(g_sExeName,'/')!=NULL) || (strchr(g_sExeName,'\\')!=NULL))
		return;

	tmp = getenv("PATH");
	if (tmp == NULL)
	{
		BTOUT;
		return;
	}

	try { env = new char[strlen(tmp)+1]; } catch(...) { env = NULL; }
	if (env == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(env,tmp);
	p = env;
	while (true)
	{
#ifdef TARGET_WINDOWS
		q = strchr(p,';');
#else
		q = strchr(p,':');
#endif
		if (q != NULL)
			*q = 0;

		try { tmp = new char[strlen(p)+strlen(g_sExeName)+2]; } catch(...) { tmp = NULL; }
		if (tmp == NULL) NewException((double)(strlen(p)+strlen(g_sExeName)+2)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		strcpy(tmp,p);
#ifdef TARGET_WINDOWS
		strcat(tmp,"\\");
#else
		strcat(tmp,"/");
#endif
		strcat(tmp,g_sExeName);
		if (FileExist(tmp))
		{
			delete[] g_sExeName;

			try { g_sExeName = new char[strlen(tmp)+1]; } catch(...) { g_sExeName = NULL; }
			if (g_sExeName == NULL) NewException((double)(strlen(tmp)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(g_sExeName,tmp);
			delete[] env;
			BTOUT;
			return;
		}
		delete[] tmp;
		if (q != NULL)
			p = q+1;
				else break;
	}
	delete[] env;
	BTOUT;
}


void ReorderAtoms(int molecule)
{ // Written @ Ballmer Peak
	int z, z2, z3;
	CMolecule *mol;
	CSingleMolecule *sm;
	CxIntArray *twa;

	mol = (CMolecule*)g_oaMolecules[molecule];

	twa = NULL;

	for (z=0;z<mol->m_laSingleMolIndex.GetSize();z++)
	{
		sm = (CSingleMolecule*)g_oaSingleMolecules[mol->m_laSingleMolIndex[z]];
//		sm->m_oaTempAtomOffset.SetSize(mol->m_baAtomIndex.GetSize());

		for (z2=0;z2<mol->m_baAtomIndex.GetSize();z2++)
		{
/*			if (sm->m_oaTempAtomOffset[z2] != NULL)
				delete sm->m_oaTempAtomOffset[z2];
			sm->m_oaTempAtomOffset[z2] = new CxIntArray();*/
			if (twa != NULL)
				delete twa;

			try { twa = new CxIntArray("ReorderAtoms():twa"); } catch(...) { twa = NULL; }
			if (twa == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			twa->SetSize(((CxIntArray*)mol->m_oaNewNumbers[z2])->GetSize());
			for (z3=0;z3<((CxIntArray*)mol->m_oaNewNumbers[z2])->GetSize();z3++)
			{
//				ti = ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3);
//				(*((CxIntArray*)sm->m_oaAtomOffset[z2]))[z3] = ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(((CxIntArray*)mol->m_oaNewNumbers[z2])->GetAt(z3));
				(*twa)[z3] = ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(((CxIntArray*)mol->m_oaNewNumbers[z2])->GetAt(z3));
//				mprintf("  (%s%d wird zu %s%d; %d -> %d)\n",((CAtom*)g_oaAtoms[mol->m_baAtomIndex[z2]])->m_sName,z3+1,((CAtom*)g_oaAtoms[mol->m_baAtomIndex[z2]])->m_sName,((CxIntArray*)mol->m_oaNewNumbers[z2])->GetAt(z3)+1,((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3),((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(((CxIntArray*)mol->m_oaNewNumbers[z2])->GetAt(z3)));
			}
			for (z3=0;z3<((CxIntArray*)mol->m_oaNewNumbers[z2])->GetSize();z3++)
				(*((CxIntArray*)sm->m_oaAtomOffset[z2]))[z3] = (*twa)[z3];
		}
	}
	delete twa;
}


void ReorderLikeInput()
{ // Written @ Ballmer Peak
	int z, z2, z3, z4, z5, ti, tl;
	CMolecule *mol;
	CSingleMolecule *sm;

	mprintf("\n");
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		mol = (CMolecule*)g_oaMolecules[z];
		mprintf("Reordering atoms in %s...",mol->m_sName);
		mol->m_oaNewNumbers.SetSize(mol->m_baAtomIndex.GetSize());
		sm = (CSingleMolecule*)g_oaSingleMolecules[mol->m_laSingleMolIndex[0]];
		for (z2=0;z2<mol->m_baAtomIndex.GetSize();z2++)
		{
			if (mol->m_oaNewNumbers[z2] != NULL)
				delete mol->m_oaNewNumbers[z2];

			try { mol->m_oaNewNumbers[z2] = new CxIntArray("ReorderLikeInput():mol->m_oaNewNumbers[z2]"); } catch(...) { mol->m_oaNewNumbers[z2] = NULL; }
			if (mol->m_oaNewNumbers[z2] == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			// F***ing Stack sort is incompatible with Ballmer Peak ^^
			for (z3=0;z3<mol->m_waAtomCount[z2];z3++)
			{
				// Sortiere nach kleinstem Offset
				tl = 9999999;
				ti = -1;
				for (z4=0;z4<mol->m_waAtomCount[z2];z4++)
				{
					for (z5=0;z5<z3;z5++)
						if (((CxIntArray*)mol->m_oaNewNumbers[z2])->GetAt(z5) == z4)
							goto _nextreord; // Der ist schon einsortiert. Keyswapping hier offensichtlich nicht moeglich
					if (((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z4) < tl)
					{
						tl = ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z4);
						ti = z4;
					}
_nextreord:;
				}
				((CxIntArray*)mol->m_oaNewNumbers[z2])->Add(ti);
			}
		}
		ReorderAtoms(z);
		mprintf("Done.\n");
	}
	mprintf("\n");
}


unsigned long GraceColor(int z, double bleach)
{
	int r, g, b;

	z = z % 14;
	switch(z)
	{
		case 0: r = 0; g = 0; b = 0; break;
		case 1: r = 0; g = 0; b = 255; break;
		case 2: r = 255; g = 0; b = 0; break;
		case 3: r = 0; g = 255; b = 0; break;
		case 4: r = 255; g = 0; b = 255; break;
		case 5: r = 0; g = 255; b = 255; break;
		case 6: r = 255; g = 255; b = 0; break;
		case 7: r = 128; g = 128; b = 128; break;
		case 8: r = 0; g = 0; b = 128; break;
		case 9: r = 128; g = 0; b = 0; break;
		case 10: r = 0; g = 128; b = 0; break;
		case 11: r = 128; g = 0; b = 128; break;
		case 12: r = 0; g = 128; b = 128; break;
		case 13: r = 128; g = 255; b = 0; break;
		default: r = 0; g = 0; b = 0;
	}
	r = (int)((1.0-bleach)*r + bleach*255.0);
	g = (int)((1.0-bleach)*g + bleach*255.0);
	b = (int)((1.0-bleach)*b + bleach*255.0);
	if (r > 255)
		r = 255;
	if (g > 255)
		g = 255;
	if (b > 255)
		b = 255;
	return r + g*0x100 + b*0x10000;
}


unsigned long CalcFFTSize(unsigned long i, bool silent)
{
	unsigned long t, r, r2, r3, r5, t2, t3, t5;
	int i2, i3, i5, m2, m3, m5;
	bool b=false;

//	mprintf("*** CalcFFTSize %d ***\n",i);

	m2 = (int)ceil(log((double)i) / log(2.0));
	m3 = (int)ceil(log((double)i) / log(3.0));
	m5 = (int)ceil(log((double)i) / log(5.0));

	r = 2147483647;
	r2 = m2;
	r3 = m3;
	r5 = m5;
//	mprintf("m2=%d, m3=%d, m5=%d, r=%d.\n",m2,m3,m5,r);
	t2 = 1;
	for (i2=0;i2<=m2;i2++)
	{
		t3 = 1;
		for (i3=0;i3<=m3;i3++)
		{
			t5 = 1;
			for (i5=0;i5<=m5;i5++)
			{
				t = t2 * t3 * t5;
	//			mprintf("%d * %d * %d = %d\n",t2,t3,t5,t);
				if ((t >= i) && (t < r))
				{
					r = t;
					r2 = i2;
					r3 = i3;
					r5 = i5;
//					mprintf("2^%d + 3^%d +5^%d = %d.\n",r2,r3,r5,r);
				}
				t5 *= 5;
			}
			t3 *= 3;
		}
		t2 *= 2;
	}
	if (!silent)
	{
		mprintf("\n    CalcFFTSize(): %lu = ",i);
		if (r2 != 0)
		{
			mprintf("2^%lu",r2);
			b = true;
		}
		if (r3 != 0)
		{
			if (b)
				mprintf(" *");
			mprintf(" 3^%lu",r3);
			b = true;
		}
		if (r5 != 0)
		{
			if (b)
				mprintf(" *");
			mprintf(" 5^%lu",r5);
		}
		if ((i - r) == 0)
			mprintf(". All prime factors within 2, 3, 5. Fine.\n\n");
				else mprintf(" - %lu. Prime factors other than 2, 3, 5 not allowed.\n",r-i);
	}
//	mprintf("*** CalcFFTSize done ***\n");
	return r;
}


void BuildAtomIndices()
{
	int z, z2, z3, z4, z5, ti, ti2;
	CMolecule *m;
	CSingleMolecule *sm;

	g_waAtomMolIndex.SetSize(g_iGesVirtAtomCount);
	g_waAtomMolUID.SetSize(g_iGesVirtAtomCount);
	g_laAtomSMLocalIndex.SetSize(g_iGesVirtAtomCount);
	g_waAtomMolNumber.SetSize(g_iGesVirtAtomCount);
	g_waAtomElement.SetSize(g_iGesVirtAtomCount);
	g_waAtomRealElement.SetSize(g_iGesVirtAtomCount);
	//g_faAtomCode.SetSize(g_iGesAtomCount);
	g_liAtomCode.resize(g_iGesAtomCount);
	g_faVdWRadius.SetSize(g_iGesAtomCount);

	for (z=0;z<g_iGesVirtAtomCount;z++)
	{
		g_waAtomMolIndex[z] = 60000;
		g_waAtomRealElement[z] = 60000;
	}

	ti2 = 0;
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		if (m->m_bPseudo)
			continue;
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			for (z3=0;z3<sm->m_baAtomIndex.GetSize();z3++)
			{
				for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
				{
					ti = ((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4);
					if (g_waAtomMolIndex[ti] != 60000)
						eprintf("BuildAtomIndices(): Atom in more than one molecule! z=%d, z2=%d, z3=%d, z4=%d.\n",z,z2,z3,z4);
					g_waAtomMolIndex[ti] = (unsigned short)z;
					g_waAtomMolUID[ti] = (unsigned short)ti2;
					g_laAtomSMLocalIndex[ti] = z2;
					g_waAtomElement[ti] = (unsigned short)z3;
					g_waAtomMolNumber[ti] = (unsigned short)z4;
					g_waAtomRealElement[ti] = sm->m_baAtomIndex[z3];

					if (sm->m_baAtomIndex[z3] == g_iVirtAtomType)
						continue;

					g_faVdWRadius[ti] = ((CAtom*)g_oaAtoms[sm->m_baAtomIndex[z3]])->m_pElement->m_fVdWRadius;

					for (z5=0;z5<sm->m_oaMolAtoms.GetSize();z5++)
					{
						if (((CMolAtom*)sm->m_oaMolAtoms[z5])->m_iOffset == ti)
						{
							//g_faAtomCode[ti] = ((CMolAtom*)sm->m_oaMolAtoms[z5])->m_fAtomCode;
							g_liAtomCode[ti] = ((CMolAtom*)sm->m_oaMolAtoms[z5])->m_liAtomCode;
							goto _found;
						}
					}
					eprintf("BuildAtomIndices(): Atom Code not found (z=%d, z2=%d, z3=%d, z4=%d).\n",z,z2,z3,z4);
_found:;
				}
			}
			ti2++;
		}
	}
	g_baAtomPassedCondition.SetSize(g_iGesVirtAtomCount);
	for (z=0;z<g_iGesVirtAtomCount;z++)
		g_baAtomPassedCondition[z] = 0;
}


bool DetermineTrajFormat()
{
	char *p;


	// DLPoly trajectories have no file extension
	if (strlen(g_sInputTraj) >= 7)
	{
		if (mystricmp(&g_sInputTraj[strlen(g_sInputTraj)-7],"HISTORY") == 0)
		{
			g_iTrajFormat = 4;
			return true;
		}
	}

	// VASP trajectories have no file extension
	if (strlen(g_sInputTraj) >= 6) {
		if (mystricmp(&g_sInputTraj[strlen(g_sInputTraj)-6],"POSCAR") == 0) {
			g_iTrajFormat = 11;
			return true;
		}
	}
	if (strlen(g_sInputTraj) >= 7) {
		if (mystricmp(&g_sInputTraj[strlen(g_sInputTraj)-7],"XDATCAR") == 0) {
			g_iTrajFormat = 11;
			return true;
		}
	}

	// All other file formats require an extension
	p = strrchr(g_sInputTraj,'.');
	if (p == NULL)
		goto _unk;
	p++;
	if (mystricmp(p,"xyz")==0)
	{
		g_iTrajFormat = 0;
		return true;
	}
	if (mystricmp(p,"pdb")==0)
	{
		g_iTrajFormat = 1;
		return true;
	}
	if (mystricmp(p,"mol2")==0)
	{
		g_iTrajFormat = 2;
		return true;
	}
	if (mystricmp(p,"lmp")==0)
	{
		g_iTrajFormat = 3;
		return true;
	}
	if (mystricmp(p,"cube")==0)
	{
		g_iTrajFormat = 5;
		g_bVolumetricData = true;
		return true;
	}
	if ((mystricmp(p,"prmtop")==0) || (mystricmp(p,"mdcrd")==0))
	{
		g_iTrajFormat = 6;
		*p = 0; // Delete file extension
		return true;
	}
	if ( (mystricmp(p,"parm7")==0) || (mystricmp(p,"nc")==0) )
	{
#ifdef USE_NETCDF
		g_iTrajFormat = 6;
		*p = 0; // Delete file extension
		return true;
#else
		mprintf(WHITE,"The file extension determines the file to be of new Amber type.\n");
		mprintf(WHITE,"This is only supported if TRAVIS is compiled with activated\n");
		mprintf(WHITE,"NetCDF support. Please, see config.h for further information.\n\n");
		return false;
#endif
	}
	if ((mystricmp(p,"bqb")==0) || (mystricmp(p,"bbq")==0) || (mystricmp(p,"btr")==0) || (mystricmp(p,"emp")==0) || (mystricmp(p,"blist")==0))
	{
		g_iTrajFormat = 7;
		return true;
	}
	if (mystricmp(p,"gro")==0)
	{
		g_iTrajFormat = 8;
		return true;
	}
	if (mystricmp(p,"voronoi")==0)
	{
		g_iTrajFormat = 9;
		return true;
	}
	if (mystricmp(p,"dcd")==0)
	{
		g_iTrajFormat = 10;
		return true;
	}
_unk:
	if (p == NULL)
		eprintf("Could not determine trajectory file format: No file extension found.\n\n");
	else
		eprintf("Could not determine trajectory file format from extension (\"%s\").\n\n",p-1);
	mprintf(WHITE,"The following formats are supported:\n");
	mprintf("  - XYZ trajectories (extension .xyz)\n");
	mprintf("  - PDB trajectories (extension .pdb)\n");
	mprintf("  - LAMMPS trajectories \"dump custom element xu yu zu\" (extension .lmp)\n");
	mprintf("  - DLPOLY trajectories (file name \"HISTORY\", no extension)\n");
	mprintf("  - VASP trajectories (file name \"POSCAR\" or \"XDATCAR\", no extension)\n");
	mprintf("  - Gromacs text format trajectories (extension .gro)\n");
	mprintf("  - DCD binary trajectories (extension .dcd)\n");
#ifdef USE_NETCDF
	mprintf("  - Amber trajectories (extension .parm7 / .nc)\n");
#else
	mprintf("  - Amber trajectories (extension .prmtop / .mdcrd)\n");
#endif
	mprintf("  - Cube file trajectories (extension .cube)\n");
	mprintf("  - CP2k Voronoi text trajectories (extension .voronoi)\n");
	mprintf("  - BQB format (extensions .bqb, .btr, .emp, and .blist)\n");
	mprintf("\n");
	return false;
}


const char* GetFileExtension(const char *s) {

	const char *p;

	if (s == NULL)
		return NULL;

	p = strrchr(s,'.');
	if (p == NULL)
		return NULL;

	return p+1;
}


void PrintSMode()
{
	const char *stext[16];
	char buf[256];
	int z;

	stext[0]  = "ICAgU1NTU1NTU1NTU1NTU1NTICAgICAgICAgICAgICAgICAgaGhoaGhoICAgICAgICAgICAgICAgaWlpICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgISEh";
	stext[1]  = "IFNTOjo6Ojo6Ojo6Ojo6Ojo6UyAgICAgICAgICAgICAgICAgaDo6OjpoICAgICAgICAgICAgICBpOjo6aSAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAhITohIQ==";
	stext[2]  = "Uzo6Ojo6U1NTU1NTOjo6Ojo6UyAgICAgICAgICAgICAgICAgaDo6OjpoICAgICAgICAgICAgICAgaWlpICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAhOjo6IQ==";
	stext[3]  = "Uzo6Ojo6UyAgICAgU1NTU1NTUyAgICAgICAgICAgICAgICAgIGg6OjpoICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAhOjo6IQ==";
	stext[4]  = "Uzo6Ojo6UyAgICAgICAgICAgICAgICAgY2NjY2NjY2NjY2NjIGg6OjpoIGhoaGhoICAgICAgIGlpaWlpaSAgICAgIHNzc3Nzc3NzICAgICAgICBzc3Nzc3NzcyAgICAhOjo6IQ==";
	stext[5]  = "Uzo6Ojo6UyAgICAgICAgICAgICAgIGNjOjo6Ojo6Ojo6OjpjIGg6OjpoaDo6Ojo6aGggICAgIGk6Ojo6aSAgICBzczo6Ojo6Ojo6cyAgICAgc3M6Ojo6Ojo6OnMgICAhOjo6IQ==";
	stext[6]  = "IFM6Ojo6U1NTUyAgICAgICAgICAgYzo6Ojo6Ojo6Ojo6OjpjIGg6Ojo6Ojo6Ojo6OjpoaCAgICBpOjo6aSAgc3M6Ojo6Ojo6Ojo6OnMgIHNzOjo6Ojo6Ojo6OjpzICAhOjo6IQ==";
	stext[7]  = "ICBTUzo6Ojo6OlNTU1NTICAgICBjOjo6OjpjY2NjY2M6OjpjIGg6Ojo6OjpoaGg6Ojo6OmggICBpOjo6aSAgczo6Ojo6c3Nzczo6OjpzIHM6Ojo6OnNzc3M6Ojo6cyAhOjo6IQ==";
	stext[8]  = "ICAgIFNTUzo6Ojo6Ojo6U1MgICBjOjo6OmMgICAgIGNjY2NjIGg6Ojo6OmggICBoOjo6OjpoICBpOjo6aSAgIHM6Ojo6cyAgIHNzc3MgICBzOjo6OnMgICBzc3NzICAhOjo6IQ==";
	stext[9]  = "ICAgICAgIFNTU1NTUzo6OjpTICBjOjo6YyAgICAgICAgICAgIGg6Ojo6aCAgICAgaDo6OjpoICBpOjo6aSAgICAgczo6Ojo6cyAgICAgICAgIHM6Ojo6OnMgICAgICAhOjo6IQ==";
	stext[10] = "ICAgICAgICAgICAgUzo6Ojo6UyBjOjo6YyAgICAgICAgICAgIGg6Ojo6aCAgICAgaDo6OjpoICBpOjo6aSAgICAgICBzOjo6OjpzICAgICAgICAgczo6Ojo6cyAgICAhITohIQ==";
	stext[11] = "ICAgICAgICAgICAgUzo6Ojo6UyBjOjo6OmMgICAgIGNjY2NjIGg6Ojo6aCAgICAgaDo6OjpoICBpOjo6aSAgc3NzcyAgICBzOjo6OnMgIHNzc3MgICAgczo6OjpzICAgISEh";
	stext[12] = "U1NTU1NTUyAgICAgUzo6Ojo6UyBjOjo6OjpjY2NjY2M6OjpjIGg6Ojo6aCAgICAgaDo6OjpoIGk6Ojo6Omkgczo6Ojpzc3NzOjo6OjpzIHM6Ojo6c3Nzczo6Ojo6cw==";
	stext[13] = "Uzo6Ojo6OlNTU1NTUzo6Ojo6UyAgYzo6Ojo6Ojo6Ojo6OjpjIGg6Ojo6aCAgICAgaDo6OjpoIGk6Ojo6Omkgczo6Ojo6Ojo6Ojo6OnMgIHM6Ojo6Ojo6Ojo6OjpzICAgISEh";
	stext[14] = "Uzo6Ojo6Ojo6Ojo6Ojo6OlNTICAgIGNjOjo6Ojo6Ojo6OjpjIGg6Ojo6aCAgICAgaDo6OjpoIGk6Ojo6OmkgIHM6Ojo6Ojo6OjpzcyAgICBzOjo6Ojo6Ojo6c3MgICAhITohIQ==";
	stext[15] = "IFNTU1NTU1NTU1NTU1NTUyAgICAgICAgY2NjY2NjY2NjY2NjIGhoaGhoaCAgICAgaGhoaGhoIGlpaWlpaWkgICBzc3Nzc3Nzc3MgICAgICAgc3Nzc3Nzc3NzICAgICAgISEh";

	for (z=0;z<16;z++)
	{
		UnBase64((unsigned char*)buf,(const unsigned char*)stext[z],(int)strlen(stext[z]));
		mprintf(YELLOW,"%s\n",buf);
	}
	mprintf("\n");
}


void PrintBMode() {

	const char *stext[15];
	char buf[256];
	int z;

	stext[0]   = "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF8sIC4u";
	stext[1]   = "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICdfKF8wbylvKCku";
	stext[2]   = "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAsbyhfXyksXyksKF9fKW8=";
	stext[3]   = "ICAgICAgICAgXy4oLSkuXyAgICAgICAgICAgICAgICAgICAgbyhfLC1vKF8gKSgoXylvTylf";
	stext[4]   = "ICAgICAgLicgICAgICAgICAnLiAgICAgICAgICAgICAgICAgLk8oX18pbyxfXylvKF8pT29fKQ==";
	stext[5]   = "ICAgICAvICAgICAgICAgICAgIFwgICAgICAgICAgICAuLS0tLXwgICB8ICAgfCAgIHxfKTA=";
	stext[6]   = "ICAgICB8Jy0uLi5fX18uLi4tJ3wgICAgICAgICAgIC8gIC4tLXwgICB8ICAgfCAgIHwsXyk=";
	stext[7]   = "ICAgICAgXCAgICAnPScgICAgLyAgICAgICAgICAgfCAgLyAgIHwgICB8ICAgfCAgIHxvKF8p";
	stext[8]   = "ICAgICAgIGAnLl9fX19fLidgICAgICAgICAgICAgfCAgfCAgIHwgICB8ICAgfCAgIHxfL2Ap";
	stext[9]   = "ICAgICAgICAvICAgfCAgIFwgICAgICAgICAgICAgfCAgfCAgIHwgICB8ICAgfCAgIHxPXyk=";
	stext[10]  = "ICAgICAgIC8uLS0nfCctLS5cICAgICAgICAgICAgfCAgXCAgIHwgICB8ICAgfCAgIHw=";
	stext[11]  = "ICAgIFtdLyctLl9ffF9fLi0nXFtdICAgICAgICAgIFwgICctLXwgICB8ICAgfCAgIHw=";
	stext[12]  = "ICAgICAgICAgICAgfCAgICAgICAgICAgICAgICAgICAnLS0tLXwgICB8ICAgfCAgIHw=";
	stext[13]  = "ICAgICAgICAgICBbXSAgICAgICAgICAgICAgICAgICAgICAgIFwgICBcICAgLyAgIC8=";
	stext[14]  = "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBgIiIiIiIiIiIiYA==";

	mprintf("\n");
	for (z=0;z<15;z++)
	{
		UnBase64((unsigned char*)buf,(const unsigned char*)stext[z],(int)strlen(stext[z]));
		mprintf(YELLOW,"    %s\n",buf);
	}
//	mprintf("\n");
}


void WriteCredits()
{
	bool tb;


	mprintf("\n");
	if (g_bShowCredits)
	{
		WriteCredits_Long();
	} else
	{
		mprintf(YELLOW,"    Note:");
		mprintf(" To show a list of all persons who contributed to TRAVIS,\n");
		mprintf("          please add \"");
		mprintf(WHITE,"-credits");
		mprintf("\" to your command line arguments, or set the\n");
		mprintf("          variable \"SHOWCREDITS\" to \"TRUE\" in your travis.conf file.\n");
	}
	mprintf("\n");
	mprintf("    Source code from other projects used in TRAVIS:\n");
	mprintf("      - lmfit     from "); mprintf(WHITE,"Joachim Wuttke\n");
	mprintf("      - kiss_fft  from "); mprintf(WHITE,"Mark Borgerding\n");
	mprintf("      - voro++    from "); mprintf(WHITE,"Chris Rycroft\n");
	mprintf("\n");
	mprintf(WHITE,"    http://www.travis-analyzer.de\n\n");
	mprintf(RED,"    Please cite all of the following articles for the analyses you have used:\n\n");

	bool citeall = false;

	mprintf(YELLOW,"  * For TRAVIS in general:\n\n");

	mprintf(WHITE,"    \"TRAVIS - A Free Analyzer for Trajectories from Molecular Simulation\",\n");
    mprintf("    M. Brehm, M. Thomas, S. Gehrke, B. Kirchner; J. Chem. Phys. 2020, 152 (16), 164105.   (DOI 10.1063/5.0005078 )\n\n");

	mprintf(WHITE,"    \"TRAVIS - A Free Analyzer and Visualizer for Monte Carlo and Molecular Dynamics Trajectories\",\n");
	mprintf("    M. Brehm, B. Kirchner; J. Chem. Inf. Model. 2011, 51 (8), pp 2007-2023.   (DOI 10.1021/ci200217w )\n\n");

	if (g_bUseBQB || citeall) {
		mprintf(YELLOW,"  * For the compressed BQB file format:\n\n");
		mprintf(WHITE,"    \"An Efficient Lossless Compression Algorithm for Trajectories of Atom Positions and Volumetric Data\",\n");
		mprintf("    M. Brehm, M. Thomas; J. Chem. Inf. Model. 2018, 58 (10), pp 2092-2107.   (DOI 10.1021/acs.jcim.8b00501 )\n\n");
	}

	if (g_bIRSpec || g_bPowerSpec || g_bRaman || g_bPower || g_bIR || g_bROA || citeall)
	{
		mprintf(YELLOW,"  * For computing spectra from AIMD simulations:\n\n");
		mprintf(WHITE,"    \"Computing Vibrational Spectra from ab initio Molecular Dynamics\",\n");
		mprintf("    M. Thomas, M. Brehm, R. Fligg, P. Voehringer, B. Kirchner; Phys. Chem. Chem. Phys. 2013, 15, pp 6608-6622.\n    (DOI 10.1039/C3CP44302G )\n\n");
	}

	if (g_bTegri || g_bROA || g_bMinimizeChargeVar2 || g_bVoro || g_bDomA || g_bHole || citeall) {
		if (g_bMinimizeChargeVar2)
			mprintf(YELLOW,"  * For Voronoi charges via charge variance minimization:\n\n");
		else
			mprintf(YELLOW,"  * For the radical Voronoi tessellation:\n\n");
		mprintf(WHITE,"    \"Optimized Atomic Partial Charges and Radii Defined by Radical Voronoi Tessellation of Bulk Phase Simulations\",\n");
		mprintf("    M. Brehm, M. Thomas; Molecules 2021, 26 (7), 1875.   (DOI 10.3390/molecules26071875 )\n\n");
	}

	if (g_bTegri || g_bROA || g_bMinimizeChargeVar2 || citeall) {
		mprintf(YELLOW,"  * For Voronoi integration of molecular electromagnetic moments:\n\n");
		mprintf(WHITE,"    \"Voronoi Dipole Moments for the Simulation of Bulk Phase Vibrational Spectra\",\n");
		mprintf("    M. Thomas, M. Brehm, B. Kirchner; Phys. Chem. Chem. Phys. 2015, 17, pp 3207-3213.   (DOI 10.1039/C4CP05272B )\n\n");
	}


	if (g_bNormalCoordinate || g_bEckartTransform || citeall) {
		mprintf(YELLOW,"  * For bulk phase normal mode analysis:\n\n");
		mprintf(WHITE,"    \"Simulating the Vibrational Spectra of Ionic Liquid Systems: 1-Ethyl-3-methylimidazolium Acetate and its Mixtures\",\n");
		mprintf("    M. Thomas, M. Brehm, O. Holloczki, Z. Kelemen, L. Nyulaszi, T. Pasinszki, B. Kirchner;\n    J. Chem. Phys. 2014, 141, 024510.   (DOI 10.1063/1.4887082 )\n\n");
	}

	if (g_bVoro || g_bDomA || g_bHole || citeall) {
		mprintf(YELLOW,"  * For Voronoi analysis and Voronoi-based domain analysis:\n\n");
		mprintf(WHITE,"    \"Domain Analysis in Nanostructured Liquids: A Post-Molecular Dynamics Study at the Example of Ionic Liquids\",\n");
		mprintf("    M. Brehm, H. Weber, M. Thomas, O. Holloczki, B. Kirchner; ChemPhysChem 2015, 16, pp 3271-3277.\n    (DOI 10.1002/cphc.201500471 )\n\n");
	}

	if (g_bSFac || citeall) {
		mprintf(YELLOW,"  * For structure factor computations:\n\n");
		mprintf(WHITE,"    \"Triphilic Ionic-Liquid Mixtures: Fluorinated and Non-fluorinated Aprotic Ionic-Liquid Mixtures\",\n");
		mprintf("    O. Holloczki, M. Macchiagodena, H. Weber, M. Thomas, M. Brehm, A. Stark, O. Russina, A. Triolo, B. Kirchner;\n    ChemPhysChem 2015, 16, pp 3325-3333.   (DOI 10.1002/cphc.201500473 )\n\n");
	}

	if (g_bKuehneLong) {
		mprintf(YELLOW,"  * For long-range RDF computations:\n\n");
		mprintf(WHITE,"    \"Optimal Calculation of the Pair Correlation Function for an Orthorhombic System\",\n");
		mprintf("    K. A. F. Roehrig, T. D. Kuehne; Phys. Rev. E 2013, 87, 045301.   (DOI 10.1103/PhysRevE.87.045301 )\n\n");
	}

	if (g_bClusterAnalysis) {
		mprintf(YELLOW,"  * For the hierarchical cluster analysis:\n\n");
		mprintf(WHITE,"    \"Cluster Analysis in Liquids: A Novel Tool in TRAVIS\",\n");
		mprintf("    T. Froembgen, J. Blasius, V. Alizadeh, A. Chaumont, M. Brehm, B. Kirchner; submitted.\n\n");
	}

	tb = false;
	if (g_bVCD || g_bCubeTimeDev || citeall)
		tb = true;
	if (g_bROA)
		if (g_pROAEngine->m_bVCD || g_pROAEngine->m_bROA || g_pROAEngine->m_bMagMom)
			tb = true;
	if (tb) {
		mprintf(YELLOW,"  * For magnetic moment calculation from BOMD:\n\n");
		mprintf(WHITE,"    \"Classical Magnetic Dipole Moments for the Simulation of Vibrational Circular Dichroism by\n      ab initio Molecular Dynamics\",\n");
		mprintf("    M. Thomas, B. Kirchner; J. Phys. Chem. Lett. 2016, 7, pp 509-513.   (DOI 10.1021/acs.jpclett.5b02752 )\n\n");
	}

	tb = false;
	if (citeall)
		tb = true;
	if (g_bROA)
		if (g_pROAEngine->m_bROA)
			tb = true;
	if (tb) {
		mprintf(YELLOW,"  * For computing Raman optical activity spectra:\n\n");
		mprintf(WHITE,"    \"Computing Bulk Phase Raman Optical Activity Spectra from ab initio Molecular Dynamics Simulations\",\n");
		mprintf("    M. Brehm, M. Thomas; J. Phys. Chem. Lett. 2017, 8 (14), pp 3409-3414.   (DOI 10.1021/acs.jpclett.7b01616 )\n\n");
	}

	if (g_bHBond || g_bIonPair || citeall) {
		mprintf(YELLOW,"  * For the reactive flux analysis:\n\n");
		mprintf(WHITE,"    \"Structure and Lifetimes in Ionic Liquids and their Mixtures\",\n");
		mprintf("    S. Gehrke, M. von Domaros, R. Clark, O. Holloczki, M. Brehm, T. Welton, A. Luzar, B. Kirchner;\n    Faraday discuss. 2018, 206, pp 219-245.   (DOI 10.1039/C7FD00166E )\n\n");
		mprintf(WHITE,"    \"Robustness of the Hydrogen Bond and Ion Pair Dynamics in Ionic Liquids to Dierent Parameters from the Reactive Flux Method\",\n");
		mprintf("    S. Gehrke, B. Kirchner; J. Chem. Eng. Data 2020, 65 (3), pp  1146-1158.   (DOI 10.1021/acs.jced.9b00529 )\n\n");
	}

	if (g_bHole || citeall) {
		mprintf(YELLOW,"  * For the void analysis:\n\n");
		mprintf(WHITE,"    \"Understanding the Fluidity of Condensed Phase Systems in Terms of Voids - Novel Algorithm, Implementation and Application\",\n");
		mprintf("    S. Gehrke, R. Macchieraldo, B. Kirchner;\n    Phys. Chem. Chem. Phys. 2019, 21 (9), pp 4988-4997.   (DOI 10.1039/C8CP07120A )\n\n");
	}

	if (g_bOrder || citeall) {
		mprintf(YELLOW,"  * For computing order parameters:\n\n");
		mprintf(WHITE,"    \"Influence of Small Fluorophilic and Lipophilic Organic Molecules on Dipalmitoylphosphatidylcholine Bilayers\",\n");
		mprintf("    M. Brehm, G. Saddiq, T. Watermann, D. Sebastiani; J. Phys. Chem. B 2017, 121 (35), 8311-8321.\n    (DOI 10.1021/acs.jpcb.7b06520 )\n\n");
	}
}


void WriteCredits_Long()
{
	mprintf(YELLOW,"  >>> Credits <<<\n");
	mprintf(YELLOW,"\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Martin Brehm              "); mprintf("Main Developer (2009 - now)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Martin Thomas             "); mprintf("Developer (2012 - now)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Sascha Gehrke             "); mprintf("Developer (2016 - now)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Barbara Kirchner          "); mprintf("Group Leader and idea\n\n");

	mprintf(YELLOW,"    > "); mprintf(WHITE,"Marc Bruessel             "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Philipp di Dio            "); mprintf("Atom Parameters and Math Support\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Christian Dressler        "); mprintf("Bug Finding and Discussion\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Michael von Domaros       "); mprintf("New Ideas and Linux Repository updating\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Dzmitry Firaha            "); mprintf("Testing and new Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Dorothea Golze            "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Oldamur Holloczki         "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Daniela Kerle             "); mprintf("Name Finding\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Miriam Kohagen            "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Matthias Krack            "); mprintf("Testing and valuable suggestions\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Marina Macchiagodena      "); mprintf("Testing and new Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Friedrich \"Fred\" Malberg  "); mprintf("Leading Bug Finder :-)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Christopher Peschel       "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Eliane Roos               "); mprintf("Force field development and testing\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Matthias Schoeppke        "); mprintf("Testing and new Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Daniele Tedesco           "); mprintf("Testing and new Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Jens Thar                 "); mprintf("Fruitful Discussion and Scientific Input\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Henry Weber               "); mprintf("Testing, Debugging and creative Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Moritz Weiss              "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Stefan Zahn               "); mprintf("Testing and Debugging\n");

/*	mprintf(YELLOW,"  >>> Credits <<<\n");
	mprintf(YELLOW,"\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Martin Brehm "); mprintf("("); mprintf(RED,"*"); mprintf(")    "); mprintf("Main Developer (2009 - now)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Martin Thomas "); mprintf("("); mprintf(RED,"*"); mprintf(")   "); mprintf("Main Developer (2012 - now)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Barbara Kirchner    "); mprintf("Group Leader\n\n");

	mprintf(YELLOW,"    > "); mprintf(WHITE,"Marc Bruessel       "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Philipp di Dio "); mprintf("("); mprintf(RED,"*"); mprintf(")  "); mprintf("Atom Parameters and Math Support\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Dorothea Golze      "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Daniela Kerle       "); mprintf("Name Finding\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Miriam Kohagen      "); mprintf("Testing and Debugging\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Fred Malberg        "); mprintf("Leading Bug Finder :-)\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Matthias Schoeppke  "); mprintf("Testing and new Ideas\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Jens Thar "); mprintf("("); mprintf(RED,"*"); mprintf(")       Fruitful Discussion and Scientific Input\n");
	mprintf(YELLOW,"    > "); mprintf(WHITE,"Henry Weber "); mprintf("("); mprintf(RED,"*"); mprintf(")     ");mprintf("Testing, Debugging and creative Ideas\n");
	mprintf("\n");
	mprintf("    ("); mprintf(RED,"*"); mprintf(")"); mprintf(" Office 42 - \"answer to all questions\".\n");*/
}


CAutoCorrelation::CAutoCorrelation()
{
	m_iInput = 0;
	m_iDepth = 0;
	m_pFFT = NULL;
	m_pFFT2 = NULL;
}


CAutoCorrelation::~CAutoCorrelation()
{
}


void CAutoCorrelation::Init(int input, int depth, bool fft)
{
	m_iInput = input;
	m_iDepth = depth;
	if (fft)
	{
		m_bFFT = true;
		m_iFFTSize = CalcFFTSize(input,true);

		try { m_pFFT = new CFFT(); } catch(...) { m_pFFT = NULL; }
		if (m_pFFT == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pFFT->PrepareFFT_C2C(2*m_iFFTSize);

		try { m_pFFT2 = new CFFT(); } catch(...) { m_pFFT2 = NULL; }
		if (m_pFFT2 == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pFFT2->PrepareInverseFFT_C2C(2*m_iFFTSize);
	} else
	{
		m_bFFT = false;
	}
}


void CAutoCorrelation::AutoCorrelate(CxDoubleArray *inp, CxDoubleArray *outp)
{
	int z, z2;
	double tf;

	outp->SetSize(m_iDepth);

	if (m_bFFT)
	{
		for (z=0;z<m_iInput;z++)
		{
			m_pFFT->m_pInput[z*2] = (*inp)[z];
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		for (z=m_iInput;z<2*m_iFFTSize;z++)
		{
			m_pFFT->m_pInput[z*2] = 0;
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		m_pFFT->DoFFT();

		for (z=0;z<m_iFFTSize*2;z++)
		{
			m_pFFT2->m_pInput[z*2] = m_pFFT->m_pOutput[z*2]*m_pFFT->m_pOutput[z*2] + m_pFFT->m_pOutput[z*2+1]*m_pFFT->m_pOutput[z*2+1];
			m_pFFT2->m_pInput[z*2+1] = 0;
		}
		m_pFFT2->DoFFT();

		for (z=0;z<m_iDepth;z++)
			(*outp)[z] = m_pFFT2->m_pOutput[2*z] / m_iFFTSize / 2.0 / ((double)m_iInput - z);
	} else
	{
		for (z=0;z<m_iDepth;z++) // Tau
		{
			tf = 0;
			for (z2=0;z2<m_iInput-z;z2++)
				tf += (*inp)[z2] * (*inp)[z2+z];
			(*outp)[z] = tf / (double)(m_iInput-z);
		}
	}
}


void CAutoCorrelation::AutoCorrelate(std::vector<double> &inp, std::vector<double> &outp) {

	int z, z2;
	double tf;

	outp.resize(m_iDepth);

	if (m_bFFT) {

		for (z=0;z<m_iInput;z++) {
			m_pFFT->m_pInput[z*2] = inp[z];
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		for (z=m_iInput;z<2*m_iFFTSize;z++) {
			m_pFFT->m_pInput[z*2] = 0;
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		m_pFFT->DoFFT();

		for (z=0;z<m_iFFTSize*2;z++) {
			m_pFFT2->m_pInput[z*2] = m_pFFT->m_pOutput[z*2]*m_pFFT->m_pOutput[z*2] + m_pFFT->m_pOutput[z*2+1]*m_pFFT->m_pOutput[z*2+1];
			m_pFFT2->m_pInput[z*2+1] = 0;
		}
		m_pFFT2->DoFFT();

		for (z=0;z<m_iDepth;z++)
			outp[z] = m_pFFT2->m_pOutput[2*z] / m_iFFTSize / 2.0 / ((double)m_iInput - z);

	} else {

		for (z=0;z<m_iDepth;z++) { // Tau
			tf = 0;
			for (z2=0;z2<m_iInput-z;z2++)
				tf += inp[z2] * inp[z2+z];
			outp[z] = tf / (double)(m_iInput-z);
		}
	}
}


void CAutoCorrelation::AutoCorrelateSqrt(CxDoubleArray *inp, CxDoubleArray *outp) {

	int z;

	outp->SetSize(m_iDepth);

	if (m_bFFT) {

		for (z=0;z<m_iInput;z++) {
			m_pFFT->m_pInput[z*2] = (*inp)[z];
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		for (z=m_iInput;z<2*m_iFFTSize;z++) {
			m_pFFT->m_pInput[z*2] = 0;
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		m_pFFT->DoFFT();

		for (z=0;z<m_iFFTSize*2;z++) {
			m_pFFT2->m_pInput[z*2] = sqrt(m_pFFT->m_pOutput[z*2]*m_pFFT->m_pOutput[z*2] + m_pFFT->m_pOutput[z*2+1]*m_pFFT->m_pOutput[z*2+1]);
//			m_pFFT2->m_pInput[z*2] = fabs(m_pFFT->m_pOutput[z*2]*m_pFFT->m_pOutput[z*2] + m_pFFT->m_pOutput[z*2+1]*m_pFFT->m_pOutput[z*2+1]);
			m_pFFT2->m_pInput[z*2+1] = 0;
		}
		m_pFFT2->DoFFT();

		for (z=0;z<m_iDepth;z++)
			(*outp)[z] = m_pFFT2->m_pOutput[2*z]; // / m_iFFTSize / 2.0 / ((double)m_iInput - z);

	} else {

		eprintf("CAutoCorrelation::AutoCorrelateSqrt(): Error: Only implemented for FFT.\n");
		abort();
	}
}


/*
void FormatTime(unsigned long eta, char *buf)
{
	char tbuf[256], tbuf2[256];

	if ((eta/60) > 0)
		sprintf(tbuf,"%02lus",eta%60);
			else sprintf(tbuf,"%2lus",eta%60);
	eta /= 60;
	if (eta > 0)
	{
		strcpy(tbuf2,tbuf);
		if ((eta/60) > 0)
			sprintf(tbuf,"%02lum",eta%60);
				else sprintf(tbuf,"%2lum",eta%60);
		strcat(tbuf,tbuf2);
		eta /= 60;
		if (eta > 0)
		{
			strcpy(tbuf2,tbuf);
			if ((eta/60) > 0)
				sprintf(tbuf,"%02luh",eta%24);
					else sprintf(tbuf,"%2luh",eta%24);
			strcat(tbuf,tbuf2);
			eta /= 24;
		}
		if (eta > 0)
		{
			strcpy(tbuf2,tbuf);
			if ((eta/60) > 0)
				sprintf(tbuf,"%02lud",eta);
					else sprintf(tbuf,"%2lud",eta);
			strcat(tbuf,tbuf2);
		}
	}
	strcpy(buf,tbuf);
}
*/


void FormatTime(unsigned long eta, CxString *buf)
{
//	char tbuf[256], tbuf2[256];
	CxString tbuf, tbuf2;

	if (buf == NULL)
		return;

	if ((eta/60) > 0)
//		sprintf(tbuf,"%02lus",eta%60);
		tbuf.sprintf("%02lus",eta%60);
	else
//		sprintf(tbuf,"%2lus",eta%60);
		tbuf.sprintf("%2lus",eta%60);

	eta /= 60;
	if (eta > 0)
	{
//		strcpy(tbuf2,tbuf);
		tbuf2.strcpy(tbuf);

		if ((eta/60) > 0)
//			sprintf(tbuf,"%02lum",eta%60);
			tbuf.sprintf("%02lum",eta%60);
		else
//			sprintf(tbuf,"%2lum",eta%60);
			tbuf.sprintf("%2lum",eta%60);

//		strcat(tbuf,tbuf2);
		tbuf.strcat(tbuf2);

		eta /= 60;
		if (eta > 0)
		{
//			strcpy(tbuf2,tbuf);
			tbuf2.strcpy(tbuf);

			if ((eta/60) > 0)
//				sprintf(tbuf,"%02luh",eta%24);
				tbuf.sprintf("%02luh",eta%24);
			else
//				sprintf(tbuf,"%2luh",eta%24);
				tbuf.sprintf("%2luh",eta%24);

//			strcat(tbuf,tbuf2);
			tbuf.strcat(tbuf2);

			eta /= 24;
		}
		if (eta > 0)
		{
//			strcpy(tbuf2,tbuf);
			tbuf2.strcpy(tbuf);

			if ((eta/60) > 0)
//				sprintf(tbuf,"%02lud",eta);
				tbuf.sprintf("%02lud",eta);
			else
//				sprintf(tbuf,"%2lud",eta);
				tbuf.sprintf("%2lud",eta);

//			strcat(tbuf,tbuf2);
			tbuf.strcat(tbuf2);
		}
	}

//	strcpy(buf,tbuf);
	buf->strcpy(tbuf);
}


void RenderStructFormulas(int tries, bool allsm)
{
	int z, z2, z3, z0, ti;
	bool nohyd, hex, onlyh;
	double tf;
	FILE *a;
//	char buf[256];
	CxString buf;
	CMolecule *m;
	CSingleMolecule *sm;
	CMolAtom *ma1, *ma2;

	mprintf("\n");

	hex = false;

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];

		if ((m->m_iAtomGesNoVirt > 150) || (m->m_bPolymer)) {
			mprintf("\n");
			mprintf(YELLOW,"    Warning: ");
			mprintf("Molecule %d (%s) is rather large. Rendering of structural formula may take a lot of time.\n\n",z+1,(const char*)m->m_sName);
			if (!AskYesNo("    Render this structural formula anyways (y/n)? [no] ",false))
				continue;
			onlyh = AskYesNo("    Show only non-hydrogen atoms in this structural formula (y/n)? [yes] ",true);
			mprintf("\n");
		} else
			onlyh = false;

		if (m->m_iAtomGesNoVirt == 1) {
			mprintf(WHITE,"    * Molecule %d: %s - skipping (only one atom).\n",z+1,m->m_sName);
			continue;
		}

		for (z0=0;z0<(allsm?m->m_laSingleMolIndex.GetSize():1);z0++) {

			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z0]];
			nohyd = onlyh;

	_again:
			if (allsm) {
				if (nohyd) {
					buf.sprintf("mol%d_%s_sm%d_no_H.dot",z+1,m->m_sName,z0+1);
					mprintf(WHITE,"    * Molecule %d: %s, SM %d (without hydrogen atoms)\n",z+1,m->m_sName,z0+1);
				} else {
					buf.sprintf("mol%d_%s_sm%d.dot",z+1,m->m_sName,z0+1);
					mprintf(WHITE,"    * Molecule %d: %s, SM %d\n",z+1,m->m_sName,z0+1);
				}
			} else {
				if (nohyd) {
					buf.sprintf("mol%d_%s_no_H.dot",z+1,m->m_sName);
					mprintf(WHITE,"    * Molecule %d: %s (without hydrogen atoms)\n",z+1,m->m_sName);
				} else {
					buf.sprintf("mol%d_%s.dot",z+1,m->m_sName);
					mprintf(WHITE,"    * Molecule %d: %s\n",z+1,m->m_sName);
				}
			}
			mprintf("      Writing dot file %s ...\n",(const char*)buf);
			a = OpenFileWrite(buf,true);
			mfprintf(a,"graph molecule {\n");
			mfprintf(a,"  graph [pack=true,splines=true,overlap=false];\n");
			mfprintf(a,"  node [shape=none,fontsize=16,fontname=\"Arial\",margin=0,fixedsize=true,height=0.28];\n");
			mfprintf(a,"  edge [style=bold,len=0.70];\n");

			for (z2=0;z2<sm->m_oaMolAtoms.GetSize();z2++)
			{
				ma1 = (CMolAtom*)sm->m_oaMolAtoms[z2];
				if (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,"H")==0) {
					if (nohyd)
						continue;
					else
						hex = true;
				}
				mfprintf(a,"  %s%d [label=\"%s%d\",width=%.2f];\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,ma1->m_iNumber+1,(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,ma1->m_iNumber+1,(strlen((const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName)+(((ma1->m_iNumber+1)>9)?2:1))*0.17);
			}

			for (z2=0;z2<sm->m_oaMolAtoms.GetSize();z2++)
			{
				ma1 = (CMolAtom*)sm->m_oaMolAtoms[z2];
				if (nohyd && (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,"H")==0))
					continue;
				for (z3=0;z3<ma1->m_oaBonds.GetSize();z3++)
				{
					ma2 = (CMolAtom*)ma1->m_oaBonds[z3];
					if (ma2->m_iMolAtomNumber < ma1->m_iMolAtomNumber)
						continue;
					if (nohyd && (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma2->m_iType]])->m_sName,"H")==0))
						continue;
					// Calculate edge weight: C-C is most important
					if ((mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,"C")==0) && (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma2->m_iType]])->m_sName,"C")==0))
					{
						ti = 10000;
						//tf = 2.0;
					} else 
					{
						ti = ((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_pElement->m_iOrd * ((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma2->m_iType]])->m_pElement->m_iOrd;
						//tf = 1.0;
					}
					if (nohyd)
						ti = 40;
					tf = 0.5 + (((ti>80)?80:ti)/80.0*2.0);
					mfprintf(a,"  %s%d -- %s%d [weight=%d, penwidth=%f];\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma1->m_iType]])->m_sName,ma1->m_iNumber+1,(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[ma2->m_iType]])->m_sName,ma2->m_iNumber+1,ti,tf);
				}
			}
			mfprintf(a,"}\n");
			fclose(a);

			//if (nohyd)
			//	buf.sprintf("mol%d_%s_no_H",z+1,m->m_sName);
			//else
			//	buf.sprintf("mol%d_%s",z+1,m->m_sName);

			if (allsm) {
				if (nohyd)
					buf.sprintf("mol%d_%s_sm%d_no_H",z+1,m->m_sName,z0+1);
				else
					buf.sprintf("mol%d_%s_sm%d",z+1,m->m_sName,z0+1);
			} else {
				if (nohyd)
					buf.sprintf("mol%d_%s_no_H",z+1,m->m_sName);
				else
					buf.sprintf("mol%d_%s",z+1,m->m_sName);
			}

			RenderFormula(buf,tries);
			mprintf("\n");

			if ((!nohyd) && (m->m_iAtomGes > 30) && hex)
			{
				nohyd = true;
				goto _again;
			}
		}
	}
	mprintf("    If the command above worked, you can now view the PNG files that have been created.\n\n");
}


void RenderFormula(const char *s, int tries)
{
	int z, zm;
	char buf[256], buf2[256], *p, *q;
	double tf, mi, ma, av;
	FILE *a;
	
//	mprintf("%d Tries:\n",tries);
	zm = -1;
	av = 0;
	ma = 0;
	mi = 1e30;
	mprintf("      Command: dot %s.dot -Tpng -o%s.png -Kneato -Gepsilon=0.000001 -Gnslimit=5000 -Gmclimit=5000 -v -Gstart=%d\n",s,s,rand());
	mprintf("      Optimizing (%d tries): ",tries);
	mprintf(WHITE,"[");
	for (z=0;z<tries;z++)
	{
		mprintf(WHITE,"#");
		sprintf(buf,"dot %s.dot -Tpng -o%s.%d.png -Kneato -Gepsilon=0.000001 -Gnslimit=5000 -Gmclimit=5000 -v -Gstart=%d > dot%d.log 2>&1",s,s,z,rand(),z);

//		mprintf("    (%s)\n",buf);

		(void)!system(buf);

		sprintf(buf,"dot%d.log",z);
		a = fopen(buf,"rt");
		if (a == NULL)
		{
			eprintf("\nRenderFormula(): Error opening %s for reading.\n",buf);
			eprintf("It seems that GraphViz is not installed or in system path (try command \"dot\" in command line).\n\n");
			return;
		}
		while (!feof(a))
		{
			(void)!fgets(buf,256,a);
			if (strstr(buf,"final") != NULL)
			{
				p = buf;
	//			mprintf("    buf=\"%s\".\n",buf);
				while (((*p < '0') || (*p > '9')) && (*p != 0))
					p++;
	//			mprintf("    p=\"%s\".\n",p);
				if (*p == 0)
					continue;
				q = p;
//				mprintf("    *p = %d.\n",*p);
				while (((*p >= '0') && (*p <= '9')) || (*p == '.'))
				{
	//				mprintf("    *p = %d --> p++;\n",*p);
					p++;
				}
	//			mprintf("    p2=\"%s\".\n",p);
				*p = 0;
		//		mprintf("    q=\"%s\".\n",q);
				tf = atof(q);
//				mprintf("    tf = %g.\n",tf);
				if (tf > ma)
					ma = tf;
				if (tf < mi)
				{
					mi = tf;
					zm = z;
				}
				av += tf;
				goto _done;
			}
		}
		eprintf("\nError: Unexpected GraphViz output. See dot%d.log.\n",z);
		eprintf("It seems that GraphViz is not working properly (try command \"dot\" in command line).\n\n");
		return;
_done:
		fclose(a);
		sprintf(buf,"dot%d.log",z);
		remove(buf);
	}
	mprintf(WHITE,"]\n");
	av /= tries;
	mprintf("      Quality statistics: Min %g, Max %g, Avg %g.\n",mi,ma,av);
	mprintf("      Using best result to produce output file.\n");
	for (z=0;z<tries;z++)
	{
		if (z == zm)
			continue;
		sprintf(buf,"%s.%d.png",s,z);
		if (remove(buf) != 0)
			eprintf("      Error removing \"%s\".\n",buf);
	}
	sprintf(buf,"%s.%d.png",s,zm);
	sprintf(buf2,"%s.png",s);
	remove(buf2);
	if (rename(buf,buf2) != 0)
		eprintf("      Error renaming \"%s\" to \"%s\".\n",buf,buf2);
}


void parseCoreCharges() {
	static bool coreChargesDefined = false;
	if (coreChargesDefined)
		return;
	
	mprintf(WHITE, "\n  > Core Charge Definition >\n\n");
		mprintf("    Please mind pseudopotentials while entering core charges.\n\n");
		int i;
		for (i = 0; i < g_oaAtoms.GetSize(); i++) {
			if (i == g_iWannierAtomType)
				continue;
			if (i == g_iVirtAtomType)
				continue;
			double def = 0.0;
			if (((CAtom *)g_oaAtoms[i])->m_fDefCharge == 0) {
				if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "B") == 0) def = 3.0;
				else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "C") == 0) def = 4.0;
				else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "N") == 0) def = 5.0;
				else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "O") == 0) def = 6.0;
				else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "F") == 0) def = 7.0;
				else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Na") == 0) def = 9.0;
				else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Si") == 0) def = 4.0;
				else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "P") == 0) def = 5.0;
				else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "S") == 0) def = 6.0;
				else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Cl") == 0) def = 7.0;
				else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Br") == 0) def = 7.0;
				else if (mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "I") == 0) def = 7.0;
				else def = (double)((CAtom *)g_oaAtoms[i])->m_pElement->m_iOrd;
			} else
				def = ((CAtom *)g_oaAtoms[i])->m_fDefCharge;
			((CAtom *)g_oaAtoms[i])->m_fCharge = AskFloat("    Enter core charge for atom type %s: [%.1f] ", def, (const char*)((CAtom*)g_oaAtoms[i])->m_sName, def);
		}
	mprintf(WHITE, "\n  < End Core Charge Definition <\n");
	
	coreChargesDefined = true;
}


bool setupWannier() {
	if(!g_bWannier) {
		mprintf(WHITE, "\n  > Wannier Centers >\n\n");
		g_bWannier = true;
		int watom = -1;
		int i;
		for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
			if(((CAtom *)g_oaAtoms[i])->m_bExclude) {
				watom = i;
				break;
			}
		}
		if(watom == -1)
			eprintf("    Atom label of Wannier centers not found.\n\n");
		else
			mprintf("    Atom type %s is excluded from the system, probably these are the Wannier centers.\n\n", (const char*)((CAtom *)g_oaAtoms[watom])->m_sName);
		
		bool ok = false;
		while(!ok) {
			CxString buf, buf2;
			
			buf.sprintf("    Which atom label do the wannier centers have (");
			for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
				buf2.sprintf("%s", (const char*)((CAtom *)g_oaAtoms[i])->m_sName);
				buf.strcat(buf2);
				if(i < g_oaAtoms.GetSize() - 2) {
					buf2.sprintf(", ");
					buf.strcat(buf2);
				}
			}
			buf2.sprintf(")? ");
			buf.strcat(buf2);
			if(watom != -1) {
				buf2.sprintf("[%s] ", (const char*)((CAtom *)g_oaAtoms[watom])->m_sName);
				buf.strcat(buf2);
			}
			
			if(watom == -1)
				AskString_ND("%s", &buf2, (const char*)buf);
			else
				AskString("%s", &buf2, ((CAtom *)g_oaAtoms[watom])->m_sName, (const char*)buf);
			
			for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
				if(mystricmp(buf2, ((CAtom *)g_oaAtoms[i])->m_sName) == 0) {
					g_iWannierAtomType = i;
					ok = true;
					break;
				}
			}
			if(!ok) {
				eprintf("    Invalid input.\n");
				inpprintf("! Invalid input.\n");
			}
		}
		
		g_fWannierCharge = fabs(AskFloat("\n    Enter the negative charge of the Wannier centers (without sign): [2.0] ", 2.0));
		parseCoreCharges();
		mprintf("\n");
		if(!g_TimeStep.ScanWannier(true)) {
			eprintf("\n    Setup of Wannier centers failed.\n");
			return false;
		}
		mprintf(WHITE, "\n  < End Wannier Centers <\n");
	}
	return true;
}


static void setupFluctuatingCharges()
{
	if ((!g_bReadChargesFrom4thXYZ) && (!g_bReadLAMMPSCharges))
	{
		if (g_bXYZ4thCol)
			g_bReadChargesFrom4thXYZ = true;

		if (g_bLAMMPSCharge)
			g_bReadLAMMPSCharges = true;

		mprintf(WHITE, "\n  > Fluctuating Atomic Partial Charges >\n");
		
		FILE *a;
		int z, z2, z3, z4;
		double tf;
		CMolecule *m;
		CSingleMolecule *sm;
		mprintf("\n    Trying to read charges from first time step...\n");
		a = fopen(g_sInputTraj,"rb");
		if (a == NULL)
		{
			eprintf("setupFluctuatingCharges(): Error: Could not open \"%s\" for reading.\n",g_sInputTraj);
			abort();
		}
		
		if (!g_TimeStep.ReadTimestep(a,true))
		{
			eprintf("setupFluctuatingCharges(): Error: Could not read first time step of \"%s\".\n",g_sInputTraj);
			abort();
		}
		
		fclose(a);
		
		mprintf("\n    Uniting molecules which have been broken by wrapping...\n");
		g_TimeStep.UniteMolecules(true);

		g_TimeStep.CalcCenters();
		
		mprintf("\n    Found the following atomic charges:\n\n");
		
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			m = (CMolecule*)g_oaMolecules[z];
			m->m_bChargesAssigned = true;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
					continue;
				for (z3=0;z3<m->m_waAtomCount[z2];z3++)
				{
					for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
					{
						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
						mprintf("      %s[%d] %s%d:  %8.6f\n",m->m_sName,z4+1,(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)]);
					}
				}
			}
		}
		mprintf("\n    Total charges of the molecules:\n\n");
		
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			m = (CMolecule*)g_oaMolecules[z];
			for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
				tf = 0;
				for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
				{
					if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
						continue;
					
					for (z3=0;z3<m->m_waAtomCount[z2];z3++)
						tf += g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)];
				}
				mprintf("      %s[%d]:  Total charge is %8.6f.\n",m->m_sName,z4+1,tf);
				if (z4 == 0)
					m->m_fCharge = tf;
			}
		}
		mprintf(WHITE, "\n  < End Fluctuating Atomic Partial Charges <\n");
	}
}


static void setupDipoleRestartFile() {
	if (!g_bLoadDipoleRestart) {
		mprintf(WHITE, "\n  > Dipole Restart File >\n\n");
		g_bLoadDipoleRestart = true;
		while (true) {
			CxString filename;
			AskString("    Enter name of dipole restart file: [dipole.restart] ", &filename, "dipole.restart");
			g_fDipoleRestartFile = fopen((const char *)filename, "r");
			if (g_fDipoleRestartFile == NULL) {
				eprintf("    Could not open \"%s\": %s.\n\n", (const char *)filename, strerror(errno));
				continue;
			}
			int numMolecules;
			(void)!fread(&numMolecules, sizeof(int), 1, g_fDipoleRestartFile);
			if (numMolecules != g_oaSingleMolecules.GetSize()) {
				eprintf("    This dipole restart file was written for a different number of molecules.\n\n");
				continue;
			}
			break;
		}
		mprintf(WHITE, "\n  < End Dipole Restart File <\n");
	}
}


static void setupMagneticDipoleRestartFile() {
	if (!g_bLoadMagneticDipoleRestart) {
		mprintf(WHITE, "\n  > Magnetic Moment Restart File >\n\n");
		g_bLoadMagneticDipoleRestart = true;
		while (true) {
			CxString filename;
			AskString("    Enter name of magnetic moment restart file: [magnetic.restart] ", &filename, "magnetic.restart");
			g_fMagneticDipoleRestartFile = fopen((const char *)filename, "r");
			if (g_fMagneticDipoleRestartFile == NULL) {
				eprintf("    Could not open \"%s\": %s.\n\n", (const char *)filename, strerror(errno));
				continue;
			}
			int numMolecules;
			(void)!fread(&numMolecules, sizeof(int), 1, g_fMagneticDipoleRestartFile);
			if (numMolecules != g_oaSingleMolecules.GetSize()) {
				eprintf("    This magnetic moment restart file was written for a different number of molecules.\n\n");
				continue;
			}
			break;
		}
		mprintf(WHITE, "\n  < End Magnetic Moment Restart File <\n");
	}
}


void ParseDipole()
{
	if (g_bDipoleDefined)
		return;
	
	g_bDipole = true;
	
	if (g_bTegri && (g_pTetraPak == NULL)) {
		try { g_pVoroWrapper = new CVoroWrapper(); } catch(...) { g_pVoroWrapper = NULL; }
		if (g_pVoroWrapper == NULL) NewException((double)sizeof(CVoroWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { g_pTetraPak = new CTetraPak(); } catch(...) { g_pTetraPak = NULL; }
		if (g_pTetraPak == NULL) NewException((double)sizeof(CTetraPak),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		g_pTetraPak->Parse();
	}
	
	mprintf(WHITE, "\n>>> Dipole Definition >>>\n\n");
	
	if(g_bAdvanced2) {
		g_bDipoleRefFixed = AskYesNo("    Use a box-fixed reference point for dipole calculation (y/n)? [no] ", false);
		mprintf("\n\n");
	} else {
		g_bDipoleRefFixed = false;
	}
	
	mprintf("    There are the following possibilities to provide dipole moments:\n");
	mprintf("    (1) Calculate dipole moments from Wannier centers\n        (These need to be atoms in the trajectory)\n");
	mprintf("    (2) Read dipole moments from an external file\n        (This file is expected to contain one line per timestep\n         with the three dipole vector components separated by spaces)\n");
	mprintf("    (3) Read dipole moments from the trajectory file\n        (These need to be specified in the comment line)\n");
	mprintf("    (4) Calculate dipole moments from fixed atomic partial charges\n");
	mprintf("    (5) Calculate dipole moments from fluctuating atomic partial charges\n        (These need to be in the fourth column of the trajectory)\n");
	mprintf("    (6) Calculate dipole moments from Voronoi partial charges\n        (This requires volumetric electron density data)\n");
	mprintf("    (7) Calculate Voronoi dipole moments\n        (This requires volumetric electron density data)\n");
	mprintf("    (8) Load a dipole restart file\n");
	mprintf("\n");
	
	mprintf("    The following dipole modes can be set for all molecules at once: 1, 5, 6, 7, 8\n");
	
	while (true) {
		int globalMode = AskRangeInteger("    Which dipole mode to set for all molecules at once? [none] ", 0, 8, 0);
		
		if (globalMode == 0) {
			int i;
			for (i = 0; i < g_oaMolecules.GetSize(); i++)
				((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 0;
			break;
		} else if (globalMode == 1) {
			if (setupWannier()) {
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 1;
				break;
			}
			continue;
		} else if (globalMode == 2) {
			eprintf("\n    This mode cannot be set for all molecules at once.\n\n");
			continue;
		} else if (globalMode == 3) {
			eprintf("\n    This mode cannot be set for all molecules at once.\n\n");
			continue;
		} else if (globalMode == 4) {
			eprintf("\n    This mode cannot be set for all molecules at once.\n\n");
			continue;
		} else if (globalMode == 5) {
			if (g_bXYZ4thCol || g_bLAMMPSCharge) {
				setupFluctuatingCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 5;
				break;
			}
			eprintf("\n    The trajectory does not have a fourth column.\n\n");
			continue;
		} else if (globalMode == 6) {
			if (!g_bVolumetricData) {
				eprintf("\n    This requires volumetric electron density data in each step of the trajectory (.cube or .bqb).\n\n");
				continue;
			}
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				parseCoreCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++) {
					int rty;
					ParseAtom("#2", i, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterIndex);
					((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 6;
				}
				break;
			}
			eprintf("\n    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n\n");
			continue;
		} else if (globalMode == 7) {
			if (!g_bVolumetricData) {
				eprintf("\n    This requires volumetric electron density data in each step of the trajectory (.cube or .bqb).\n\n");
				continue;
			}
			if (g_bTegri) {
				if (!g_pTetraPak->m_bVoronoiCharges) {
					eprintf("\n    Error: You need to activate \"Voronoi integration\" (question above).\n\n");
					continue;
				}
				g_bVoroIntegrateCharge = true;
				g_bVoroIntegrateDipoleMoment = true;
				mprintf("\n");

				g_bVoroIntegrateQuadrupoleMoment = AskYesNo("    Also activate quadrupole integration (y/n)? [no] ",false);

				parseCoreCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++) {
					int rty;
					ParseAtom("#2", i, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterIndex);
					((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 7;
				}
				break;
			}
			eprintf("\n    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n\n");
			continue;
		} else if (globalMode == 8) {
			setupDipoleRestartFile();
			int i;
			for (i = 0; i < g_oaMolecules.GetSize(); i++)
				((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 8;
			break;
		}
		eprintf("\n    Error: Invalid input.\n\n");
	}
	
	while (true) {
		mprintf("\n    The following dipole modes are set up:\n\n");
		int longest = 0;
		int i;
		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			int length = (int)strlen(((CMolecule *)g_oaMolecules[i])->m_sName);
			if (length > longest)
				longest = length;
		}
		CxString buf;
		buf.sprintf("      %%%ds - %%d", longest);
		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			mprintf(buf, ((CMolecule *)g_oaMolecules[i])->m_sName, ((CMolecule *)g_oaMolecules[i])->m_iDipoleMode);
			if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 0)
				mprintf(" (none)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 1)
				mprintf(" (Wannier centers)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 2)
				mprintf(" (external file)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 3)
				mprintf(" (trajectory comment line)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 4)
				mprintf(" (fixed atomic partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 5)
				mprintf(" (fluctuating atomic partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 6)
				mprintf(" (Voronoi partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 7)
				mprintf(" (Voronoi dipole moments)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 8)
				mprintf(" (restart file)\n");
			else
				mprintf("\n");
		}
		mprintf("\n");
		
		if (!AskYesNo("    Change dipole mode for a molecule (y/n)? [no] ", false))
			break;
		
		int mol;
		if (g_oaMolecules.GetSize() > 1) {
			buf.sprintf("    Change dipole mode for which molecule (");
			for (i = 0; i < g_oaMolecules.GetSize(); i++) {
				CxString buf2;
				buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i + 1);
				buf.strcat(buf2);
				if(i < g_oaMolecules.GetSize() - 1) {
					buf.strcat(", ");
				}
			}
			buf.strcat(")? ");
			
			mol = AskRangeInteger_ND("%s", 1, g_oaMolecules.GetSize(), (const char*)buf) - 1;
		} else {
			mol = 0;
		}
		
		int mode = AskRangeInteger_ND("    Which dipole mode to set for molecule %s? ", 1, 8, ((CMolecule *)g_oaMolecules[mol])->m_sName);
		
		if (mode == 0) {
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 0;
		} else if (mode == 1) {
			if (setupWannier()) {
				((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 1;
			}
		} else if (mode == 2) {
			FILE *dipoleFile = NULL;
			while (true) {
				if(dipoleFile != NULL) {
					fclose(dipoleFile);
					dipoleFile = NULL;
				}
				AskString_ND("    Name of file for dipole moments: ", &buf);
				mprintf("\n    Trying to read first line...\n\n");
				dipoleFile = fopen(buf, "r");
				if(dipoleFile == NULL) {
					eprintf("    Could not open \"%s\".\n\n", (const char*)buf);
					continue;
				}
				if (buf.fgets(1024, dipoleFile) == NULL) {
					eprintf("    Could not read first line.\n\n");
					continue;
				}
				double dipole[3];
				if(sscanf((const char *)buf, "%lf %lf %lf", &dipole[0], &dipole[1], &dipole[2]) < 3) {
					eprintf("    Could not find three real numbers in first line.\n\n");
					continue;
				}
				mprintf("    Found the following dipole moment in the first line:\n");
				mprintf("      ( %8G | %8G | %8G ) Debye\n", dipole[0], dipole[1], dipole[2]);
				rewind(dipoleFile);
				break;
			}
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 2;
			((CMolecule *)g_oaMolecules[mol])->m_pDipoleFile = dipoleFile;
			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
				eprintf("\n    Warning: There are several molecules %s. The same dipole moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
		} else if (mode == 3) {
			CxString comment;
			comment.strcpy(g_TimeStep.m_pComment);
			mprintf("    The comment line is \"%s\".\n", (const char *)comment);
			buf = "";
			const char delim[] = "\t ";
			char *tok = strtok(comment.GetWritePointer(), delim);
			int num = 0;
			while(tok != NULL) {
				double f;
				if(sscanf(tok, "%lf", &f) == 1) {
					CxString buf2;
					num++;
					if(num > 1) {
						buf2.sprintf(", %s (%d)", tok, num);
					} else {
						buf2.sprintf("%s (%d)", tok, num);
					}
					buf.strcat(buf2);
				}
				tok = strtok(NULL, delim);
			}
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[0] = AskRangeInteger("    Choose the X dipole component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[1] = AskRangeInteger("    Choose the Y dipole component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[2] = AskRangeInteger("    Choose the Z dipole component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 3;
			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
				eprintf("\n    Warning: There are several molecules %s. The same dipole moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
		} else if (mode == 4) {
			CMolecule *m = (CMolecule*)g_oaMolecules[mol];
			CSingleMolecule *sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
			m->m_bChargesAssigned = true;
			int z2, z3, ti;
			CMolAtom *ma;
			CxDoubleArray *tfa;
			bool b;
			double tf;
			//double td;
			LargeInteger li;
			CxString buf2;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				mprintf("\n");
				ma = NULL;
				//td = 1.0e50;
				li = 1;
				for (z3=0;z3<50;z3++)
					li *= 10;
				
				try { tfa = new CxDoubleArray("ParseDipole():tfa"); } catch(...) { tfa = NULL; }
				if (tfa == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				m->m_oaCharges.Add(tfa);
				_cnm:
				buf = "";
				buf2 = "";
				b = false;
				ti = 0;
				for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
				{
					if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_iType != z2)
						continue;
					
					ma = (CMolAtom*)sm->m_oaMolAtoms[z3];
					
					//if (ma->m_fAtomCode > td)
					if (ma->m_liAtomCode > li)
					{
						continue;
					}
					if (b)
					{
						//if (ma->m_fAtomCode < td)
						if (ma->m_liAtomCode < li)
						{
							continue;
						}
						ti++;
						buf2.sprintf(", %s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
						buf.strcat(buf2);
						
					} else
					{
						//if (ma->m_fAtomCode < td)
						if (ma->m_liAtomCode < li)
						{
							li = ma->m_liAtomCode;
							//td = ma->m_fAtomCode;
						}
						
						buf.sprintf("    Enter partial atomic charge on %s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
						b = true;
						ti++;
					}
				}
				if (b)
				{
					buf.strcat(": [0.0] ");
					tf = AskFloat("%s",0, (const char*)buf);
					for (z3=0;z3<ti;z3++)
						tfa->Add(tf);
					//td -= 1.0;
					li -= 1;
					goto _cnm;
				}
			}
			
			tf = 0;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				for (z3=0;z3<((CxDoubleArray*)m->m_oaCharges[z2])->GetSize();z3++)
				{
					mprintf("    %s%d: %7.4f\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,((CxDoubleArray*)m->m_oaCharges[z2])->GetAt(z3));
					tf += ((CxDoubleArray*)m->m_oaCharges[z2])->GetAt(z3);
				}
			}
			m->m_fCharge = tf;
			mprintf(WHITE,"\n    Molecule %s: Total charge is %.4f.\n\n",m->m_sName,m->m_fCharge);
			
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 4;
		} else if (mode == 5) {
			if (g_bXYZ4thCol || g_bLAMMPSCharge) {
				setupFluctuatingCharges();
				((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 5;
			}
			eprintf("    The trajectory does not have a fourth column.\n");
		} else if (mode == 6) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				parseCoreCharges();
				int rty;
				ParseAtom("#2", mol, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterIndex);
				((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 6;
			} else {
				eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			}
		} else if (mode == 7) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				g_bVoroIntegrateDipoleMoment = true;
				parseCoreCharges();
				int rty;
				ParseAtom("#2", mol, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterIndex);
				((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 7;
			} else {
				eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			}
		} else if (mode == 8) {
			setupDipoleRestartFile();
			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 8;
		} else {
			eprintf("    This is impossible.\n");
		}
	}
	
	if (!g_bDipoleRefFixed) {
		int i;
		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			if (fabs(((CMolecule *)g_oaMolecules[i])->m_fCharge) > 0.001) {
				mprintf("\n");
				mprintf("    Molecule %s is charged (%.4f). Input of a dipole reference point is required.\n",((CMolecule *)g_oaMolecules[i])->m_sName, ((CMolecule *)g_oaMolecules[i])->m_fCharge);
				int rty;
				CxString buf;
				do {
					AskString("    Please enter dipole reference point: [#2] ", &buf, "#2");
				} while (!ParseAtom(buf, i, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterIndex));
			}
		}
	}
	mprintf("\n");
	
// 	CxString buf, buf2, comment;
// 
// //	const int BUF_SIZE = 1024;
// 	
// 	if (g_bDipoleDefined)
// 		return;
// 
// 	
// 	if(g_bAdvanced2) {
// 		g_bDipoleRefFixed = AskYesNo("    Use a box-fixed reference point for dipole calculation (y/n)? [no] ", false);
// 		mprintf("\n");
// 	} else {
// 		g_bDipoleRefFixed = false;
// 	}
// 	
// 	mprintf(WHITE, "\n>>> Dipole Definition >>>\n\n");
// 	mprintf("    There are different ways to provide dipole moments in TRAVIS:\n");
// 	mprintf("    (1) Calculate dipole moments from Wannier centers\n        (These need to be atoms in the trajectory)\n");
// 	mprintf("    (2) Read dipole moments from an external file\n        (This file is expected to contain one line per timestep\n         with the three dipole vector components separated by spaces)\n");
// 	mprintf("    (3) Read dipole moments from the trajectory file\n        (These need to be specified in the comment line)\n");
// 	mprintf("    (4) Calculate dipole moments from fixed atomic partial charges\n");
// 	mprintf("    (5) Calculate dipole moments from fluctuating atomic partial charges\n        (These need to be in the fourth column of the trajectory)\n");
// 	mprintf("    (8) Load dipole moments from a dipole restart file\n");
// 	mprintf("\n");
// 	
// 	while(true) {
// //		char buf[BUF_SIZE];
// //		char buf2[BUF_SIZE];
// //		size_t remaining = BUF_SIZE;
// 		int mol;
// 		if(g_oaMolecules.GetSize() > 1) {
// /*#ifdef TARGET_LINUX
// 			remaining -= snprintf(buf, remaining, "    Add dipole information for which molecule (");
// #else
// 			remaining -= sprintf(buf, "    Add dipole information for which molecule (");
// #endif*/
// 
// 			buf.sprintf("    Add dipole information for which molecule (");
// 
// 			int i;
// 			for(i = 0; i < g_oaMolecules.GetSize(); i++) {
// 
// /*				if(remaining < 1)
// 					break;
// #ifdef TARGET_LINUX
// 				size_t length = snprintf(buf2, remaining, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
// #else
// 				size_t length = sprintf(buf2, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
// #endif
// 				strncat(buf, buf2, remaining - 1);
// 				remaining -= length;*/
// 
// 				buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
// 				buf.strcat(buf2);
// 
// 				if(i < g_oaMolecules.GetSize() - 1) {
// /*#ifdef TARGET_LINUX
// 					length = snprintf(buf2, remaining, ", ");
// #else
// 					length = sprintf(buf2, ", ");
// #endif
// 					strncat(buf, buf2, remaining - 1);
// 					remaining -= length;*/
// 
// 					buf2.sprintf(", ");
// 					buf.strcat(buf2);
// 				}
// 			}
// //			strncat(buf, ")? ", remaining - 1);
// 			buf.strcat(")? ");
// 
// 			mol = AskRangeInteger_ND(buf, 1, g_oaMolecules.GetSize()) - 1;
// 		} else {
// 			mol = 0;
// 		}
// 		
// 		int mode;
// 
// 		{
// 			if(g_bXYZ4thCol)
// 				mode = AskRangeInteger("\n    Use Wannier centers (1), external file (2), comment line (3), fixed charges (4), or fluctuating charges (5) for molecule %s? [1] ", 1, 5, 1, ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 			else
// 				mode = AskRangeInteger("\n    Use Wannier centers (1), external file (2), comment line (3), or fixed charges (4) for molecule %s? [1] ", 1, 4, 1, ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		}
// 		mprintf("\n");
// 		
// 		if(mode == 1) {
// 			if(!g_bWannier) {
// 				mprintf(WHITE, "  > Wannier Centers >\n\n");
// 				g_bWannier = true;
// 				int watom = -1;
// 				int i;
// 				for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
// 					if(((CAtom *)g_oaAtoms[i])->m_bExclude) {
// 						watom = i;
// 						break;
// 					}
// 				}
// 				if(watom == -1)
// 					eprintf("    Atom label of Wannier centers not found.\n\n");
// 				else
// 					mprintf("    Atom type %s is excluded from the system, probably these are the Wannier centers.\n\n", ((CAtom *)g_oaAtoms[watom])->m_sName);
// 				
// 				bool ok = false;
// 				while(!ok) {
// /*					remaining = BUF_SIZE;
// #ifdef TARGET_LINUX
// 					remaining -= snprintf(buf, remaining, "    Which atom label do the wannier centers have (");
// #else
// 					remaining -= sprintf(buf, "    Which atom label do the wannier centers have (");
// #endif*/
// 
// 					buf.sprintf("    Which atom label do the wannier centers have (");
// 
// 					for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
// 
// /*						if(remaining < 1)
// 							break;
// #ifdef TARGET_LINUX
// 						size_t length = snprintf(buf2, remaining, "%s", ((CAtom *)g_oaAtoms[i])->m_sName);
// #else
// 						size_t length = sprintf(buf2, "%s", ((CAtom *)g_oaAtoms[i])->m_sName);
// #endif
// 						strncat(buf, buf2, remaining - 1);
// 						remaining -= length;*/
// 
// 						buf2.sprintf("%s", ((CAtom *)g_oaAtoms[i])->m_sName);
// 						buf.strcat(buf2);
// 
// 						if(i < g_oaAtoms.GetSize() - 2) {
// /*#ifdef TARGET_LINUX
// 							length = snprintf(buf2, remaining, ", ");
// #else
// 							length = sprintf(buf2, ", ");
// #endif
// 							strncat(buf, buf2, remaining - 1);
// 							remaining -= length;*/
// 
// 							buf2.sprintf(", ");
// 							buf.strcat(buf2);
// 						}
// 					}
// 
// /*#ifdef TARGET_LINUX
// 					size_t length = snprintf(buf2, remaining, ")? ");
// #else
// 					size_t length = sprintf(buf2, ")? ");
// #endif
// 					strncat(buf, buf2, remaining - 1);
// 					remaining -= length;*/
// 
// 					buf2.sprintf(")? ");
// 					buf.strcat(buf2);
// 
// 					if(watom != -1) {
// /*#ifdef TARGET_LINUX
// 						snprintf(buf2, remaining, "[%s] ", ((CAtom *)g_oaAtoms[watom])->m_sName);
// #else
// 						sprintf(buf2, "[%s] ", ((CAtom *)g_oaAtoms[watom])->m_sName);
// #endif
// 						strncat(buf, buf2, remaining - 1);*/
// 
// 						buf2.sprintf("[%s] ", ((CAtom *)g_oaAtoms[watom])->m_sName);
// 						buf.strcat(buf2);
// 					}
// 					
// 					if(watom == -1)
// 						AskString_ND(buf, &buf2);
// 					else
// 						AskString(buf, &buf2, ((CAtom *)g_oaAtoms[watom])->m_sName);
// 					
// 					for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
// 						if(mystricmp(buf2, ((CAtom *)g_oaAtoms[i])->m_sName) == 0) {
// 							g_iWannierAtomType = i;
// 							ok = true;
// 							break;
// 						}
// 					}
// 					if(!ok) {
// 						eprintf("    Wrong input.\n");
// 						inpprintf("! Wrong input.\n");
// 					}
// 				}
// 				
// 				g_fWannierCharge = fabs(AskFloat("\n    Enter the negative charge of the Wannier centers (without sign): [2.0] ", 2.0f));
// 				for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
// 					if(i == g_iWannierAtomType)
// 						continue;
// 					if(i == g_iVirtAtomType)
// 						continue;
// 					double def;
// 					if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "B") == 0) def = 3.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "C") == 0) def = 4.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "N") == 0) def = 5.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "O") == 0) def = 6.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "F") == 0) def = 7.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Si") == 0) def = 4.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "P") == 0) def = 5.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "S") == 0) def = 6.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Cl") == 0) def = 7.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "Br") == 0) def = 7.0f;
// 					else if(mystricmp(((CAtom *)g_oaAtoms[i])->m_sName, "I") == 0) def = 7.0f;
// 					else def = (double)((CAtom *)g_oaAtoms[i])->m_pElement->m_iOrd;
// 					((CAtom *)g_oaAtoms[i])->m_fCharge = AskFloat("    Enter core charge (mind pseudopotentials!) for atom type %s: [%.1f] ", def, ((CAtom*)g_oaAtoms[i])->m_sName, def);
// 				}
// 				mprintf("\n");
// 				if(!g_TimeStep.ScanWannier(true)) {
// 					eprintf("\n    No dipole information for molecule %s added.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 				}
// 				mprintf(WHITE, "\n  < End Wannier Centers <\n\n");
// 			}
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 1;
// 			mprintf("    Set dipole information for %s to Wannier centers.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		} else if(mode == 2) {
// 			FILE *dipoleFile = NULL;
// 			while(true) {
// 				if(dipoleFile != NULL) {
// 					fclose(dipoleFile);
// 					dipoleFile = NULL;
// 				}
// 				AskString_ND("    Name of file for dipole moments: ", &buf);
// 				mprintf("\n    Trying to read first line...\n\n");
// 				dipoleFile = fopen(buf, "r");
// 				if(dipoleFile == NULL) {
// 					eprintf("    Could not open \"%s\".\n\n", (const char*)buf);
// 					continue;
// 				}
// //				if(fgets(buf, BUF_SIZE, dipoleFile) == NULL) {
// 				if (buf.fgets(1024, dipoleFile) == NULL) {
// 					eprintf("    Could not read first line.\n\n");
// 					continue;
// 				}
// 				double dipole[3];
// 				if(sscanf(buf, "%f %f %f", &dipole[0], &dipole[1], &dipole[2]) < 3) {
// 					eprintf("    Could not find three real numbers in first line.\n\n");
// 					continue;
// 				}
// 				mprintf("    Found the following dipole moment in the first line:\n");
// 				mprintf("      ( %8G | %8G | %8G ) Debye\n", dipole[0], dipole[1], dipole[2]);
// 				rewind(dipoleFile);
// 				break;
// 			}
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 2;
// 			mprintf("\n    Set dipole information for %s to external file.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 			((CMolecule *)g_oaMolecules[mol])->m_pDipoleFile = dipoleFile;
// 			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
// 				eprintf("\n    Warning: There are several molecules %s. The same dipole moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		} else if(mode == 3) {
// //			char comment[BUF_SIZE];
// //			strncpy(comment, g_TimeStep.m_pComment, BUF_SIZE);
// //			comment[BUF_SIZE-1] = '\0';
// 			comment.strcpy(g_TimeStep.m_pComment);
// 			mprintf("    The comment line is \"%s\".\n", (const char*)comment);
// //			char buf[BUF_SIZE];
// //			buf[0] = '\0';
// //			size_t remaining = BUF_SIZE;
// 			buf.sprintf("");
// 			const char delim[] = "\t ";
// 			char *tok = strtok(comment.GetWritePointer(), delim);
// 			int num = 0;
// 			while(tok != NULL) {
// 				double f;
// 				if(sscanf(tok, "%f", &f) == 1) {
// 					num++;
// //					char buf2[BUF_SIZE];
// //					size_t length;
// 					if(num > 1) {
// /*#ifdef TARGET_LINUX
// 						length = snprintf(buf2, BUF_SIZE, ", %s (%d)", tok, num);
// #else
// 						length = sprintf(buf2, ", %s (%d)", tok, num);
// #endif*/
// 						buf2.sprintf(", %s (%d)", tok, num);
// 					} else {
// /*#ifdef TARGET_LINUX
// 						length = snprintf(buf2, BUF_SIZE, "%s (%d)", tok, num);
// #else
// 						length = sprintf(buf2, "%s (%d)", tok, num);
// #endif*/
// 						buf2.sprintf("%s (%d)", tok, num);
// 					}
// //					strncat(buf, buf2, remaining - 1);
// 					buf.strcat(buf2);
// /*					remaining -= length;
// 					if(remaining < 1)
// 						break;*/
// 				}
// 				tok = strtok(NULL, delim);
// 			}
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[0] = AskRangeInteger("    Choose the X dipole component: %s [1] ", 1, num, 1, (const char*)buf) - 1;
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[1] = AskRangeInteger("    Choose the Y dipole component: %s [1] ", 1, num, 1, (const char*)buf) - 1;
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleCommentIndex[2] = AskRangeInteger("    Choose the Z dipole component: %s [1] ", 1, num, 1, (const char*)buf) - 1;
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 3;
// 			mprintf("\n    Set dipole information for %s to comment line.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
// 				eprintf("\n    Warning: There are several molecules %s. The same dipole moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		} else if(mode == 4) {
// 			CMolecule *m = (CMolecule*)g_oaMolecules[mol];
// 			CSingleMolecule *sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
// 			m->m_bChargesAssigned = true;
// 			int z2, z3, ti;
// 			CMolAtom *ma;
// 			CxFloatArray *tfa;
// 			bool b;
// 			double tf;
// 			double td;
// 			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 			{
// 				mprintf("\n");
// //					mprintf("%s Anfang.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 				ma = NULL;
// 				td = 1.0e50;
// 
// 				try { tfa = new CxFloatArray("ParseDipole():tfa"); } catch(...) { tfa = NULL; }
// 				if (tfa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
// 							
// 				m->m_oaCharges.Add(tfa);
// _cnm:
// //				buf[0] = 0;
// //				buf2[0] = 0;
// 				buf.sprintf("");
// 				buf2.sprintf("");
// 				b = false;
// 				ti = 0;
// //					mprintf("D.\n");
// 				for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
// 				{
// //						mprintf("z3=%d.\n",z3);
// 					if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_iType != z2)
// 						continue;
// 
// 					ma = (CMolAtom*)sm->m_oaMolAtoms[z3];
// 
// //						mprintf("AtomType %s, z3=%d, td=%f, ac=%f.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3,td,ma->m_fAtomCode);
// 
// 					if (ma->m_fAtomCode > td)
// 					{
// //							mprintf("  AtomCode = %f > %f = td. Skip.\n",ma->m_fAtomCode,td);
// 						continue;
// 					}
// 					if (b)
// 					{
// 						if (ma->m_fAtomCode < td)
// 						{
// //								mprintf("  AtomCode = %f < %f = td. Skip.\n",ma->m_fAtomCode,td);
// 							continue;
// 						}
// 						ti++;
// //							mprintf("  Taking %s%d into account.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 
// //						sprintf(buf2,", %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// //						strcat(buf,buf2);
// 
// 						buf2.sprintf(", %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						buf.strcat(buf2);
// 
// 					} else
// 					{
// 						if (ma->m_fAtomCode < td)
// 						{
// //								mprintf("  Changing td from %f to %f.\n",td,ma->m_fAtomCode);
// //								mprintf("  Taking %s%d into account.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 							td = ma->m_fAtomCode;
// 						}
// 
// //						sprintf(buf,"    Enter partial atomic charge on %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						buf.sprintf("    Enter partial atomic charge on %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						b = true;
// 						ti++;
// 					}
// 				}
// //					mprintf("%s Mitte.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 				if (b)
// 				{
// //					strcat(buf,": [0.0] ");
// 					buf.strcat(": [0.0] ");
// 					tf = AskFloat(buf,0);
// 					for (z3=0;z3<ti;z3++)
// 						tfa->Add(tf);
// 					td -= 1.0;
// 					goto _cnm;
// 				}
// //					mprintf("%s Ende.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 			}
// 
// //			mprintf("\n");
// 			tf = 0;
// 			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 			{
// 				for (z3=0;z3<((CxFloatArray*)m->m_oaCharges[z2])->GetSize();z3++)
// 				{
// 					mprintf("    %s%d: %7.4f\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3));
// 					tf += ((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3);
// 				}
// 			}
// 			m->m_fCharge = tf;
// 			mprintf(WHITE,"\n    Molecule %s: Total charge is %.4f.\n\n",m->m_sName,m->m_fCharge);
// 			
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 4;
// 			mprintf("\n    Set dipole information for %s to fixed atomic partial charges.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		} else if(mode == 5) {
// 			if(!g_bReadChargesFrom4thXYZ) {
// 				mprintf(WHITE, "  > Fluctuating Atomic Partial Charges >\n");
// 				g_bReadChargesFrom4thXYZ = true;
// 				
// 				FILE *a;
// 				int z, z2, z3, z4;
// 				double tf;
// 				CMolecule *m;
// 				CSingleMolecule *sm;
// 				mprintf("\n    Trying to read charges from first time step...\n");
// 				a = fopen(g_sInputTraj,"rb");
// 				if (a == NULL)
// 				{
// 					eprintf("    Could not open \"%s\" for reading.\n",g_sInputTraj);
// 					abort();
// 				}
// 
// 				if (!g_TimeStep.ReadTimestep(a,true))
// 				{
// 					eprintf("    Could not read first time step of \"%s\".\n",g_sInputTraj);
// 					abort();
// 				}
// 
// 				fclose(a);
// 
// 				mprintf("\n    Uniting molecules which have been broken by wrapping...\n");
// 				g_TimeStep.UniteMolecules(true);
// 
// 				mprintf("\n    Found the following atomic charges:\n\n");
// 
// 				for (z=0;z<g_oaMolecules.GetSize();z++)
// 				{
// 					m = (CMolecule*)g_oaMolecules[z];
// 					m->m_bChargesAssigned = true;
// 					for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 					{
// 						if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
// 							continue;
// 						for (z3=0;z3<m->m_waAtomCount[z2];z3++)
// 						{
// 							for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
// 							{
// 								sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
// 								mprintf("      %s[%d] %s%d:  %8.6f\n",m->m_sName,z4+1,((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)]);
// 							}
// 						}
// 					}
// 				}
// 				mprintf("\n    Total charges of the molecules:\n\n");
// 
// 				for (z=0;z<g_oaMolecules.GetSize();z++)
// 				{
// 					m = (CMolecule*)g_oaMolecules[z];
// 					for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
// 					{
// 						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
// 						tf = 0;
// 						for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 						{
// 							if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
// 								continue;
// 
// 							for (z3=0;z3<m->m_waAtomCount[z2];z3++)
// 								tf += g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)];
// 						}
// 						mprintf("      %s[%d]:  Total charge is %8.6f.\n",m->m_sName,z4+1,tf);
// 						if (z4 == 0)
// 							m->m_fCharge = tf;
// 					}
// 				}
// 				mprintf(WHITE, "\n  < End Fluctuating Atomic Partial Charges <\n\n");
// 			}
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 5;
// 			mprintf("    Set dipole information for %s to fluctuating atomic partial charges.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 		} else if(mode == 6) {
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 6;
// 			mprintf("    Set dipole information for %s to Voronoi partial charges.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 			unsigned char rty;
// 			ParseAtom("#2", mol, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterIndex);
// 		} else if (mode == 7) {
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 7;
// 			mprintf("    Set dipole information for %s to Voronoi dipole moments.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
// 			unsigned char rty;
// 			ParseAtom("#2", mol, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterIndex);
// 		} else if (mode == 8) {
// 			((CMolecule *)g_oaMolecules[mol])->m_iDipoleMode = 8;
// // 			g_bDipoleRestartFile = fopen("dipole.restart", "r");
// 			mprintf("    Using dipole restart file.\n");
// 		} else {
// 			eprintf("Internal error.\n");
// 			abort();
// 		}
// 		
// 		if(!g_bDipoleRefFixed) {
// 			if(fabs(((CMolecule *)g_oaMolecules[mol])->m_fCharge) > 0.001f) {
// 				mprintf("\n");
// 				mprintf("    The molecule is charged (%.4f). Input of a dipole reference point is required.\n", ((CMolecule *)g_oaMolecules[mol])->m_fCharge);
// 				unsigned char rty;
// 				do {
// 					AskString("    Please enter dipole reference point: [#2] ", &buf, "#2");
// 				} while(!ParseAtom(buf, mol, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[mol])->m_iDipoleCenterIndex));
// 			}
// 		}
// 		
// 		if(!(g_oaMolecules.GetSize() > 1) || !AskYesNo("\n    Add dipole information for another molecule (y/n)? [yes] ", true))
// 			break;
// 		mprintf("\n");
// 	}
// 	mprintf("\n");
	
	size_t nameLength = 0;
	int i;
	for(i = 0; i < g_oaMolecules.GetSize(); i++) {
		CMolecule *m = (CMolecule *)g_oaMolecules[i];
		if(m->m_iDipoleMode == 1 || m->m_iDipoleMode == 4 || m->m_iDipoleMode == 5)
			if(strlen(m->m_sName) > nameLength)
				nameLength = strlen(m->m_sName);
	}
	
	if(!g_bDipoleRefFixed) {
		for(i = 0; i < g_oaMolecules.GetSize(); i++) {
			CMolecule *m = (CMolecule *)g_oaMolecules[i];
			if(m->m_iDipoleMode == 1 || m->m_iDipoleMode == 4 || m->m_iDipoleMode == 5) {
				mprintf(WHITE, "  - %s:", m->m_sName);
				size_t j;
				for(j = strlen(m->m_sName); j <= nameLength; j++)
					mprintf(" ");
				mprintf(" Dipole reference point is %s%d.\n", (const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[m->m_iDipoleCenterType]])->m_sName, m->m_iDipoleCenterIndex + 1);
			}
		}
	}
	
	mprintf("\n    Calculating dipole moments from first time step...\n\n");
	if (g_bTegri) {
		g_pTetraPak->ProcessStep(&g_TimeStep, 2);
		mprintf("\n");
	}
//	for (i=0;i<50;i++)
//		mprintf("@ %f %f %f\n",g_TimeStep.m_vaCoords[i][0],g_TimeStep.m_vaCoords[i][1],g_TimeStep.m_vaCoords[i][2]);
	g_TimeStep.CalcDipoles(true);
	for(i = 0; i < g_oaMolecules.GetSize(); i++)
		if(((CMolecule *)g_oaMolecules[i])->m_iDipoleMode == 2)
			rewind(((CMolecule *)g_oaMolecules[i])->m_pDipoleFile);
	if (g_bLoadDipoleRestart) {
		fseek(g_fDipoleRestartFile, sizeof(int), SEEK_SET);
	}
	
	for(i = 0; i < g_oaMolecules.GetSize(); i++) {
		CMolecule *m = (CMolecule *)g_oaMolecules[i];
		if(m->m_iDipoleMode != 0) {
			double min = 1.0e6, max = 0.0, ave = 0.0;
			int j;
			for(j = 0; j < m->m_laSingleMolIndex.GetSize(); j++) {
				double val = ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[j]])->m_vDipole.GetLength();
				if(val > max)
					max = val;
				if(val < min)
					min = val;
				ave += val;
			}
			ave /= (double)m->m_laSingleMolIndex.GetSize();
			mprintf(WHITE, "  - %s:", m->m_sName);
			size_t k;
			for(k = strlen(m->m_sName); k <= nameLength; k++)
				mprintf(" ");
			mprintf(" Min. %12G, Max. %12G, Avg. %12G Debye.\n", min, max, ave);
			mprintf("        Cartesian dipole vector of first molecule is ( %8G | %8G | %8G ) Debye.\n\n", ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[0], ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[1], ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[2]);
		}
	}
	mprintf("\n");
	
	if (g_bAdvanced2)
	{
		int z, z2, ti;
		CxIntArray *tia;
//		char buf[256];
		CxString buf;
		if (AskYesNo("    Write out Cartesian dipole vectors of molecules in each step (y/n)? [no] ",false))
		{
			g_bDumpDipoleVector = true;
			mprintf("\n");
			ti = 0;
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				if (AskYesNo("      Consider %s molecules for writing out dipole vectors (y/n)? [yes] ",true,((CMolecule*)g_oaMolecules[z])->m_sName))
				{
					try { tia = new CxIntArray("ParseDipole():tia"); } catch(...) { tia = NULL; }
					if (tia == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					g_oaDumpDipoleVector.Add(tia);
					AskString("        Which representants to consider (e.g. 1,3-6,8)? [all] ",&buf,"*");
					if (buf[0] == '*')
					{
						for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
							tia->Add(z2+1);
					} else ParseIntList(buf,tia);
					for (z2=0;z2<tia->GetSize();z2++)
						tia->GetAt(z2)--;
					ti += tia->GetSize();
					mprintf("\n");
				} else g_oaDumpDipoleVector.Add(NULL);
			}
			g_iDumpDipoleSMCount = ti;
			g_bDumpDipoleAbs = AskYesNo("    Write out absolute value as 4th column per molecule (y/n)? [yes] ",true);
			g_bDumpDipoleXYZ = AskYesNo("    Write out XYZ trajectory to visualize dipole vectors (y/n)? [no] ",false);

			if (g_bDumpDipoleXYZ)
			{
				g_fDumpDipoleScale = AskFloat("    Please enter dipole scale in trajectory (unit: Angstrom/Debye): [1.0] ",1.0);
				g_iDumpDipoleXYZAtoms = 0;
				for (z=0;z<g_oaMolecules.GetSize();z++)
				{
					if (g_oaDumpDipoleVector[z] == NULL)
						continue;
					g_iDumpDipoleXYZAtoms += ((CMolecule*)g_oaMolecules[z])->m_iAtomGesNoVirt * ((CxIntArray*)g_oaDumpDipoleVector[z])->GetSize();
				}
				g_iDumpDipoleXYZAtoms += g_iDumpDipoleSMCount*2;
			}
		} else g_bDumpDipoleVector = false;
		mprintf("\n");
	}

	if (g_bTegri && (g_iTrajFormat != 5) && (g_iTrajFormat != 7))
	{
		mprintf("    Rewinding density cube file to beginning...\n\n");
		fseek(g_pTetraPak->m_fCubePipe,0,SEEK_SET);
	}
	
	mprintf(WHITE, "<<< End of Dipole Definition <<<\n\n\n");
	
	g_bDipoleDefined = true;
	
//----------------- OLD -----------------------
// 	
// 	int z, z2, z3, z4, ti;
// 	bool b;
// 	char buf[256], buf2[256];
// 	double td, td2, td3;
// 	double tf;
// 	unsigned char rty;
// 	CMolecule *m;
// 	CSingleMolecule *sm;
// 	CMolAtom *ma;
// 	CxFloatArray *tfa;
// 	CxIntArray *tia;
// 	FILE *a;
// 
// 	if (g_bDipoleDefined)
// 		return;
// 
// _dipbeg:
// 	mprintf(WHITE,"\n>>> Dipole definition >>>\n\n");
// 	mprintf("    TRAVIS can calculate dipole moments either from wannier centers (you need\n");
// 	mprintf("    to have those in the trajectory then) or from fixed atomic partial charges.\n\n");
// 	g_bWannier = AskYesNo("    Obtain dipole vectors from wannier centers (y/n)? [yes] ",true);
// 	mprintf("\n");
// 	if (g_bWannier)
// 	{
// _wantype:
// 		ti = -1;
// 		for (z=0;z<g_oaAtoms.GetSize()-1;z++)
// 		{
// 			if (((CAtom*)g_oaAtoms[z])->m_bExclude)
// 			{
// 				ti = z;
// 				break;
// 			}
// 		}
// 
// 		if (ti == -1)
// 			eprintf("    Atom label of wannier centers not found.\n\n");
// 				else mprintf("    Atom type %s is excluded from system, probably the wannier centers.\n\n",((CAtom*)g_oaAtoms[ti])->m_sName);
// 		
// 		sprintf(buf,"    Which atom label do the wannier centers have (");
// 
// 		for (z=0;z<g_oaAtoms.GetSize()-1;z++)
// 		{
// 			if (z < (int)g_oaAtoms.GetSize()-2)
// 			{
// 				sprintf(buf2,"%s, ",((CAtom*)g_oaAtoms[z])->m_sName);
// 				strcat(buf,buf2);
// 			} else 
// 			{
// 				if (ti == -1)
// 					sprintf(buf2,"%s)? ",((CAtom*)g_oaAtoms[z])->m_sName);
// 						else sprintf(buf2,"%s)? [%s] ",((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[ti])->m_sName);
// 				strcat(buf,buf2);
// 			}
// 		}
// 
// 		if (ti == -1)
// 			AskString_ND(buf,buf2);
// 				else AskString(buf,buf2,((CAtom*)g_oaAtoms[ti])->m_sName);
// 
// 		for (z=0;z<g_oaAtoms.GetSize()-1;z++)
// 		{
// 			if (mystricmp(buf2,((CAtom*)g_oaAtoms[z])->m_sName)==0)
// 			{
// 				g_iWannierAtomType = z;
// 				goto _wandone;
// 			}
// 		}
// 		eprintf("    Wrong input.\n");
// 		inpprintf("! Wrong input.\n");
// 		goto _wantype;
// _wandone:
// 		g_fWannierCharge = (double)fabs(AskFloat("    Enter the negative charge of the wannier centers (without sign): [2.0] ",2.0f));
// 		mprintf("\n");
// 		for (z=0;z<g_oaAtoms.GetSize()-1;z++)
// 		{
// 			if (z == g_iWannierAtomType)
// 				continue;
// 			if (z == g_iVirtAtomType)
// 				continue;
// 
// 			// Pseudopotential core charges for frequently used atoms
// 			if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"C") == 0)
// 				z2 = 4;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"N") == 0)
// 				z2 = 5;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"O") == 0)
// 				z2 = 6;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"F") == 0)
// 				z2 = 7;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"P") == 0)
// 				z2 = 5;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"S") == 0)
// 				z2 = 6;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"Cl") == 0)
// 				z2 = 7;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"B") == 0)
// 				z2 = 3;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"Br") == 0)
// 				z2 = 7;
// 			else if (mystricmp(((CAtom*)g_oaAtoms[z])->m_sName,"I") == 0)
// 				z2 = 7;
// 			else z2 = ((CAtom*)g_oaAtoms[z])->m_pElement->m_iOrd;
// 
// 			((CAtom*)g_oaAtoms[z])->m_fCharge = AskFloat("    Enter core charge (mind pseudopotentials!) for atom type %s: [%d.0] ",(double)z2,((CAtom*)g_oaAtoms[z])->m_sName,z2);
// 		}
// 		mprintf("\n");
// 		if (!g_TimeStep.ScanWannier(true))
// 			goto _dipbeg;
// /*		if (g_bAdvanced2)
// 		{
// 			g_bUnwrapWannier = AskYesNo("    Write out trajectory with unwrapped molecules and wannier centers (y/n)? [no] ",false);
// 		} else g_bUnwrapWannier = false;*/
// 	} else // NOT WANNIER
// 	{
// 		if (g_bBetaFeatures && g_bXYZ4thCol)
// 		{
// 			if (AskYesNo("    Read atomic charges from 4th numeric column of XYZ file (y/n)? [yes] ",true))
// 			{
// 				g_bReadChargesFrom4thXYZ = true;
// 
// 				mprintf("\n    Trying to read charges from first time step...\n");
// 				a = fopen(g_sInputTraj,"rb");
// 				if (a == NULL)
// 				{
// 					eprintf("    Could not open \"%s\" for reading.\n",g_sInputTraj);
// 					abort();
// 				}
// 
// 				if (!g_TimeStep.ReadTimestep(a,true))
// 				{
// 					eprintf("    Could not read first time step of \"%s\".\n",g_sInputTraj);
// 					abort();
// 				}
// 
// 				fclose(a);
// 
// 				mprintf("\n    Uniting molecules which have been broken by wrapping...\n");
// 				g_TimeStep.UniteMolecules(true);
// 
// 				mprintf("\n    Found the following atomic charges:\n\n");
// 
// 				for (z=0;z<g_oaMolecules.GetSize();z++)
// 				{
// 					m = (CMolecule*)g_oaMolecules[z];
// 					m->m_bChargesAssigned = true;
// 					for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 					{
// 						if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
// 							continue;
// 						for (z3=0;z3<m->m_waAtomCount[z2];z3++)
// 						{
// 							for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
// 							{
// 								sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
// 								mprintf("      %s[%d] %s%d:  %8.6f\n",m->m_sName,z4+1,((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)]);
// 							}
// 						}
// 					}
// 				}
// 				mprintf("\n    Total charges of the molecules:\n\n");
// 
// 				for (z=0;z<g_oaMolecules.GetSize();z++)
// 				{
// 					m = (CMolecule*)g_oaMolecules[z];
// 					for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
// 					{
// 						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
// 						tf = 0;
// 						for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 						{
// 							if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
// 								continue;
// 
// 							for (z3=0;z3<m->m_waAtomCount[z2];z3++)
// 								tf += g_TimeStep.m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3)];
// 						}
// 						mprintf("      %s[%d]:  Total charge is %8.6f.\n",m->m_sName,z4+1,tf);
// 						if (z4 == 0)
// 							m->m_fCharge = tf;
// 					}
// 				}
// 				mprintf("\n");
// 
// 				goto _chargedone;
// 			}
// 		}
// 		for (z=0;z<g_oaMolecules.GetSize();z++)
// 		{
// 			m = (CMolecule*)g_oaMolecules[z];
// 			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
// 			mprintf(YELLOW,"\n  * Molecule %s\n\n",m->m_sName);
// 
// 			m->m_bChargesAssigned = AskYesNo("    Enter partial atomic charges for this molecule (y/n)? [yes] ",true);
// 
// 			if (!m->m_bChargesAssigned)
// 			{
// 				mprintf("\n    Dipole moment will be always zero for this molecule.\n");
// 				goto _nocharge;
// 			}
// 
// 			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 			{
// 				mprintf("\n");
// //					mprintf("%s Anfang.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 				ma = NULL;
// 				td = 1.0e50;
// 
// 				try { tfa = new CxFloatArray("ParseDipole():tfa"); } catch(...) { tfa = NULL; }
// 				if (tfa == NULL) NewException((double)sizeof(CxFloatArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
// 							
// 				m->m_oaCharges.Add(tfa);
// _cnm:
// 				buf[0] = 0;
// 				buf2[0] = 0;
// 				b = false;
// 				ti = 0;
// //					mprintf("D.\n");
// 				for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
// 				{
// //						mprintf("z3=%d.\n",z3);
// 					if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_iType != z2)
// 						continue;
// 
// 					ma = (CMolAtom*)sm->m_oaMolAtoms[z3];
// 
// //						mprintf("AtomType %s, z3=%d, td=%f, ac=%f.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3,td,ma->m_fAtomCode);
// 
// 					if (ma->m_fAtomCode > td)
// 					{
// //							mprintf("  AtomCode = %f > %f = td. Skip.\n",ma->m_fAtomCode,td);
// 						continue;
// 					}
// 					if (b)
// 					{
// 						if (ma->m_fAtomCode < td)
// 						{
// //								mprintf("  AtomCode = %f < %f = td. Skip.\n",ma->m_fAtomCode,td);
// 							continue;
// 						}
// 						ti++;
// //							mprintf("  Taking %s%d into account.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						sprintf(buf2,", %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						strcat(buf,buf2);
// 					} else
// 					{
// 						if (ma->m_fAtomCode < td)
// 						{
// //								mprintf("  Changing td from %f to %f.\n",td,ma->m_fAtomCode);
// //								mprintf("  Taking %s%d into account.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 							td = ma->m_fAtomCode;
// 						}
// 						sprintf(buf,"    Enter partial atomic charge on %s%d",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
// 						b = true;
// 						ti++;
// 					}
// 				}
// //					mprintf("%s Mitte.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 				if (b)
// 				{
// 					strcat(buf,": [0.0] ");
// 					tf = AskFloat(buf,0);
// 					for (z3=0;z3<ti;z3++)
// 						tfa->Add(tf);
// 					td -= 1.0;
// 					goto _cnm;
// 				}
// //					mprintf("%s Ende.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
// 			}
// 
// //			mprintf("\n");
// 			tf = 0;
// 			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
// 			{
// 				for (z3=0;z3<((CxFloatArray*)m->m_oaCharges[z2])->GetSize();z3++)
// 				{
// 					mprintf("    %s%d: %7.4f\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3));
// 					tf += ((CxFloatArray*)m->m_oaCharges[z2])->GetAt(z3);
// 				}
// 			}
// 			m->m_fCharge = tf;
// 			mprintf(WHITE,"\n    Molecule %s: Total charge is %.4f.\n\n",m->m_sName,m->m_fCharge);
// _nocharge:;
// 		}
// 	}
// _chargedone:
// 	if (g_bAdvanced2)
// 	{
// 		mprintf("\n");
// 		g_bDipoleRefFixed = AskYesNo("    Use a box-fixed reference point for dipole calculation (y/n)? [no] ",false);
// 		if (g_bDipoleRefFixed)
// 			mprintf("\n    Warning: Absolute dipole moments of charged molecules will be erroneous.\n    Only the derivatives will be useful.\n\n");
// 	} else g_bDipoleRefFixed = false;
// 
// 	if (!g_bDipoleRefFixed)
// 	{
// 		for (z=0;z<g_oaMolecules.GetSize();z++)
// 		{
// 			m = (CMolecule*)g_oaMolecules[z];
// 			if ((!g_bWannier) && (!m->m_bChargesAssigned))
// 				continue;
// 
// 			if (fabs(m->m_fCharge) > 0.001)
// 			{
// 				mprintf("\n");
// 				mprintf("    Molecule %s is charged (%.4f). Input of a dipole reference point is required.\n",m->m_sName,m->m_fCharge);
// 	_diprefagain:
// 				AskString("    Please enter dipole reference point for %s: [#1] ",buf,"#1",m->m_sName);
// 				if (!ParseAtom(buf,z,m->m_iDipoleCenterType,rty,m->m_iDipoleCenterIndex))
// 					goto _diprefagain;
// 			}
// 		}
// 		mprintf("\n");
// 	}
// 
// 	ti = 0;
// 	for (z=0;z<g_oaMolecules.GetSize();z++)
// 	{
// 		m = (CMolecule*)g_oaMolecules[z];
// 		if ((!g_bWannier) && (!m->m_bChargesAssigned))
// 			continue;
// 		if (((int)strlen(m->m_sName)) > ti)
// 			ti = strlen(m->m_sName);
// 	}
// 
// 	if (!g_bDipoleRefFixed)
// 	{
// 		for (z=0;z<g_oaMolecules.GetSize();z++)
// 		{
// 			m = (CMolecule*)g_oaMolecules[z];
// 			if ((!g_bWannier) && (!m->m_bChargesAssigned))
// 				continue;
// 			mprintf(WHITE,"  - %s:",m->m_sName);
// 			for (z2=strlen(m->m_sName);z2<=ti;z2++)
// 				mprintf(" ");
// 			mprintf(" Dipole reference point is %s%d.\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[m->m_iDipoleCenterType]])->m_sName,m->m_iDipoleCenterIndex+1);
// 		}
// 		mprintf("\n");
// 	}
// 
// 	mprintf("    Calculating dipole moments from first time step...\n\n");
// 	g_TimeStep.CalcDipoles();
// 
// 	sm = NULL;
// 	for (z=0;z<g_oaMolecules.GetSize();z++)
// 	{
// 		m = (CMolecule*)g_oaMolecules[z];
// 		if ((!g_bWannier) && (!m->m_bChargesAssigned))
// 			continue;
// 		td = 0;
// 		td2 = 1.0e10;
// 		td3 = 0;
// 		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
// 		{
// 			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
// 			tf = sm->m_vDipole.GetLength();
// 			if (tf > td)
// 				td = tf;
// 			if (tf < td2)
// 				td2 = tf;
// 			td3 += tf;
// 		}
// 		td3 /= (double)m->m_laSingleMolIndex.GetSize();
// 		mprintf(WHITE,"  - %s:",m->m_sName);
// 		for (z2=strlen(m->m_sName);z2<=ti;z2++)
// 			mprintf(" ");
// 		mprintf(" Min. %12G, Max. %12G, Avg. %12G Debye.\n",td2,td,td3);
// 		if (m->m_laSingleMolIndex.GetSize() == 1)
// 		{
// 			for (z2=strlen(m->m_sName);z2<=ti;z2++)
// 				mprintf(" ");
// 			mprintf("    Cartesian dipole vector is ( %8G | %8G | %8G ).\n\n",sm->m_vDipole[0],sm->m_vDipole[1],sm->m_vDipole[2]);
// 		}
// 	}
// 
// 	mprintf("\n");
// 	if (g_bAdvanced2)
// 	{
// 		if (AskYesNo("    Write out cartesian dipole vectors of molecules in each step (y/n)? [no] ",false))
// 		{
// 			g_bDumpDipoleVector = true;
// 			mprintf("\n");
// 			ti = 0;
// 			for (z=0;z<g_oaMolecules.GetSize();z++)
// 			{
// 				if (AskYesNo("      Consider %s molecules for writing out dipole vectors (y/n)? [yes] ",true,((CMolecule*)g_oaMolecules[z])->m_sName))
// 				{
// 					try { tia = new CxIntArray("ParseDipole():tia"); } catch(...) { tia = NULL; }
// 					if (tia == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
// 					
// 					g_oaDumpDipoleVector.Add(tia);
// 					AskString("        Which representants to consider (e.g. 1,3-6,8)? [all] ",buf,"*");
// 					if (buf[0] == '*')
// 					{
// 						for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
// 							tia->Add(z2+1);
// 					} else ParseIntList(buf,tia);
// 					for (z2=0;z2<tia->GetSize();z2++)
// 						tia->GetAt(z2)--;
// 					ti += tia->GetSize();
// 					mprintf("\n");
// 				} else g_oaDumpDipoleVector.Add(NULL);
// 			}
// 			g_iDumpDipoleSMCount = ti;
// 			g_bDumpDipoleAbs = AskYesNo("    Write out absolute value as 4th column per molecule (y/n)? [yes] ",true);
// 			g_bDumpDipoleXYZ = AskYesNo("    Write out XYZ trajectory to visualize dipole vectors (y/n)? [no] ",false);
// 
// 			if (g_bDumpDipoleXYZ)
// 			{
// 				g_fDumpDipoleScale = AskFloat("    Please enter dipole scale in trajectory (unit: Angstrom/Debye): [1.0] ",1.0f);
// 				g_iDumpDipoleXYZAtoms = 0;
// 				for (z=0;z<g_oaMolecules.GetSize();z++)
// 				{
// 					if (g_oaDumpDipoleVector[z] == NULL)
// 						continue;
// 					g_iDumpDipoleXYZAtoms += ((CMolecule*)g_oaMolecules[z])->m_iAtomGesNoVirt * ((CxIntArray*)g_oaDumpDipoleVector[z])->GetSize();
// 				}
// 				g_iDumpDipoleXYZAtoms += g_iDumpDipoleSMCount*2;
// 			}
// 		} else g_bDumpDipoleVector = false;
// 	}
// 
// 	mprintf(WHITE,"\n<<< End of Dipole definition <<<\n\n");
// 	
// 	g_bDipoleDefined = true;
}

void parseMagneticDipole() {
	static bool magneticDipoleDefined = false;
	if (magneticDipoleDefined)
		return;
	
	g_bUseVelocities = true;
	
	if (g_bTegri && (g_pTetraPak == NULL)) {
		try { g_pVoroWrapper = new CVoroWrapper(); } catch(...) { g_pVoroWrapper = NULL; }
		if (g_pVoroWrapper == NULL) NewException((double)sizeof(CVoroWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { g_pTetraPak = new CTetraPak(); } catch(...) { g_pTetraPak = NULL; }
		if (g_pTetraPak == NULL) NewException((double)sizeof(CTetraPak),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		g_pTetraPak->Parse();
	}
		
	mprintf(WHITE, "\n>>> Magnetic Moment Definition >>>\n\n");
	
	mprintf("    There are the following possibilities to provide magnetic moments:\n");
	mprintf("    (1) Calculate magnetic moments from Wannier centers\n        (These need to be atoms in the trajectory and have to be sorted, use \"swan\" in the main menu)\n");
	mprintf("    (2) Read magnetic moments from an external file\n        (This file is expected to contain one line per timestep\n         with the three vector components separated by spaces)\n");
	mprintf("    (3) Read magnetic moments from the trajectory file\n        (These need to be specified in the comment line)\n");
	mprintf("    (4) Calculate magnetic moments from fixed atomic partial charges\n");
	mprintf("    (5) Calculate magnetic moments from fluctuating atomic partial charges\n        (These need to be in the fourth column of the trajectory)\n");
	mprintf("    (6) Calculate magnetic moments from Voronoi partial charges\n        (This needs electron density data)\n");
	mprintf("    (7) Calculate Voronoi magnetic moments\n        (This needs electron density data)\n");
	mprintf("    (8) Load a magnetic moment restart file\n");
	mprintf("\n");
	
	mprintf("    The following magnetic moment modes can be set for all molecules at once: 1, 5, 6, 7, 8\n");
	
	while (true) {
		int globalMode = AskRangeInteger("    Which magnetic moment mode to set for all molecules at once? [none] ", 0, 8, 0);
		
		if (globalMode == 0) {
			int i;
			for (i = 0; i < g_oaMolecules.GetSize(); i++)
				((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 0;
			break;
		} else if (globalMode == 1) {
			if (setupWannier()) {
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 1;
				break;
			}
			continue;
		} else if (globalMode == 2) {
			eprintf("    This mode cannot be set for all molecules at once.\n");
			continue;
		} else if (globalMode == 3) {
			eprintf("    This mode cannot be set for all molecules at once.\n");
			continue;
		} else if (globalMode == 4) {
			eprintf("    This mode cannot be set for all molecules at once.\n");
			continue;
		} else if (globalMode == 5) {
			if (g_bXYZ4thCol || g_bLAMMPSCharge) {
				setupFluctuatingCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 5;
				break;
			}
			eprintf("    The trajectory does not have a fourth column.\n");
			continue;
		} else if (globalMode == 6) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				parseCoreCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 6;
				break;
			}
			eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			continue;
		} else if (globalMode == 7) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				g_bVoroIntegrateDipoleMoment = true;
				g_bVoroIntegrateTotalCurrent = true;
				g_bVoroIntegrateMagneticMoment = true;
				g_bCubeTimeDev = true;
				parseCoreCharges();
				int i;
				for (i = 0; i < g_oaMolecules.GetSize(); i++)
					((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 7;
				break;
			}
			eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			continue;
		} else if (globalMode == 8) {
			setupMagneticDipoleRestartFile();
			int i;
			for (i = 0; i < g_oaMolecules.GetSize(); i++)
				((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode = 8;
			break;
		}
		eprintf("    This is impossible.\n");
	}
	
	int i;
	for (i = 0; i < g_oaMolecules.GetSize(); i++) {
		int rty;
		ParseAtom("#2", i, ((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleCenterIndex);
	}
		
	while (true) {
		mprintf("\n    The following magnetic moment modes are set up:\n\n");
		int longest = 0;
		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			int length = (int)strlen(((CMolecule *)g_oaMolecules[i])->m_sName);
			if (length > longest)
				longest = length;
		}
		CxString buf;
		buf.sprintf("      %%%ds - %%d", longest);
		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			mprintf(buf, ((CMolecule *)g_oaMolecules[i])->m_sName, ((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode);
			if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 0)
				mprintf(" (none)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 1)
				mprintf(" (Wannier centers)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 2)
				mprintf(" (external file)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 3)
				mprintf(" (trajectory comment line)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 4)
				mprintf(" (fixed atomic partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 5)
				mprintf(" (fluctuating atomic partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 6)
				mprintf(" (Voronoi partial charges)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 7)
				mprintf(" (Voronoi magnetic moments)\n");
			else if (((CMolecule *)g_oaMolecules[i])->m_iMagneticDipoleMode == 8)
				mprintf(" (restart file)\n");
			else
				mprintf("\n");
		}
		mprintf("\n");
		
		if (!AskYesNo("    Change magnetic moment mode for a molecule (y/n)? [no] ", false))
			break;
		
		int mol;
		if (g_oaMolecules.GetSize() > 1) {
			buf.sprintf("    Change dipole mode for which molecule (");
			for (i = 0; i < g_oaMolecules.GetSize(); i++) {
				CxString buf2;
				buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i + 1);
				buf.strcat(buf2);
				if(i < g_oaMolecules.GetSize() - 1) {
					buf.strcat(", ");
				}
			}
			buf.strcat(")? ");
			
			mol = AskRangeInteger_ND("%s", 1, g_oaMolecules.GetSize(), (const char*)buf) - 1;
		} else {
			mol = 0;
		}
		
		int mode = AskRangeInteger_ND("    Which magnetic moment mode to set for molecule %s? ", 1, 8, ((CMolecule *)g_oaMolecules[mol])->m_sName);
		
		if (mode == 0) {
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 0;
		} else if (mode == 1) {
			if (setupWannier()) {
				((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 1;
			}
		} else if (mode == 2) {
			FILE *dipoleFile = NULL;
			while (true) {
				if(dipoleFile != NULL) {
					fclose(dipoleFile);
					dipoleFile = NULL;
				}
				AskString_ND("    Name of file for magnetic moments: ", &buf);
				mprintf("\n    Trying to read first line...\n\n");
				dipoleFile = fopen(buf, "r");
				if(dipoleFile == NULL) {
					eprintf("    Could not open \"%s\".\n\n", (const char*)buf);
					continue;
				}
				if (buf.fgets(1024, dipoleFile) == NULL) {
					eprintf("    Could not read first line.\n\n");
					continue;
				}
				double dipole[3];
				if(sscanf((const char *)buf, "%lf %lf %lf", &dipole[0], &dipole[1], &dipole[2]) < 3) {
					eprintf("    Could not find three real numbers in first line.\n\n");
					continue;
				}
				mprintf("    Found the following magnetic moment in the first line:\n");
				mprintf("      ( %8G | %8G | %8G ) Bohr magneton.\n", dipole[0], dipole[1], dipole[2]);
				rewind(dipoleFile);
				break;
			}
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 2;
			((CMolecule *)g_oaMolecules[mol])->m_pMagneticDipoleFile = dipoleFile;
			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
				eprintf("\n    Warning: There are several molecules %s. The same magnetic moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
		} else if (mode == 3) {
			CxString comment;
			comment.strcpy(g_TimeStep.m_pComment);
			mprintf("    The comment line is \"%s\".\n", (const char *)comment);
			buf = "";
			const char delim[] = "\t ";
			char *tok = strtok(comment.GetWritePointer(), delim);
			int num = 0;
			while(tok != NULL) {
				double f;
				if(sscanf(tok, "%lf", &f) == 1) {
					CxString buf2;
					num++;
					if(num > 1) {
						buf2.sprintf(", %s (%d)", tok, num);
					} else {
						buf2.sprintf("%s (%d)", tok, num);
					}
					buf.strcat(buf2);
				}
				tok = strtok(NULL, delim);
			}
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleCommentIndex[0] = AskRangeInteger("    Choose the X magnetic moment component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleCommentIndex[1] = AskRangeInteger("    Choose the Y magnetic moment component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleCommentIndex[2] = AskRangeInteger("    Choose the Z magnetic moment component: %s [1] ", 1, num, 1, (const char *)buf) - 1;
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 3;
			if(((CMolecule *)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize() > 1)
				eprintf("\n    Warning: There are several molecules %s. The same dipole moment will be used for all of them.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
		} else if (mode == 4) {
			CMolecule *m = (CMolecule*)g_oaMolecules[mol];
			CSingleMolecule *sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
			m->m_bChargesAssigned = true;
			int z2, z3, ti;
			CMolAtom *ma;
			CxDoubleArray *tfa;
			bool b;
			double tf;
			//double td;
			LargeInteger li;
			CxString buf2;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				mprintf("\n");
				ma = NULL;
				//td = 1.0e50;
				li = 1;
				for (z3=0;z3<50;z3++)
					li *= 10;
				
				try { tfa = new CxDoubleArray("ParseDipole():tfa"); } catch(...) { tfa = NULL; }
				if (tfa == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				m->m_oaCharges.Add(tfa);
				_cnm:
				buf = "";
				buf2 = "";
				b = false;
				ti = 0;
				for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
				{
					if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_iType != z2)
						continue;
					
					ma = (CMolAtom*)sm->m_oaMolAtoms[z3];
					
					//if (ma->m_fAtomCode > td)
					if (ma->m_liAtomCode > li)
					{
						continue;
					}
					if (b)
					{
						//if (ma->m_fAtomCode < td)
						if (ma->m_liAtomCode < li)
						{
							continue;
						}
						ti++;
						buf2.sprintf(", %s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
						buf.strcat(buf2);
						
					} else
					{
						//if (ma->m_fAtomCode < td)
						if (ma->m_liAtomCode < li)
						{
							li = ma->m_liAtomCode;
							//td = ma->m_fAtomCode;
						}
						
						buf.sprintf("    Enter partial atomic charge on %s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
						b = true;
						ti++;
					}
				}
				if (b)
				{
					buf.strcat(": [0.0] ");
					tf = AskFloat("%s",0, (const char*)buf);
					for (z3=0;z3<ti;z3++)
						tfa->Add(tf);
					//td -= 1.0;
					li -= 1;
					goto _cnm;
				}
			}
			
			tf = 0;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				for (z3=0;z3<((CxDoubleArray*)m->m_oaCharges[z2])->GetSize();z3++)
				{
					mprintf("    %s%d: %7.4f\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1,((CxDoubleArray*)m->m_oaCharges[z2])->GetAt(z3));
					tf += ((CxDoubleArray*)m->m_oaCharges[z2])->GetAt(z3);
				}
			}
			m->m_fCharge = tf;
			mprintf(WHITE,"\n    Molecule %s: Total charge is %.4f.\n\n",m->m_sName,m->m_fCharge);
			
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 4;
		} else if (mode == 5) {
			if (g_bXYZ4thCol || g_bLAMMPSCharge) {
				setupFluctuatingCharges();
				((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 5;
			}
			eprintf("    The trajectory does not have a fourth column.\n");
		} else if (mode == 6) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				parseCoreCharges();
				((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 6;
			} else {
				eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			}
		} else if (mode == 7) {
			if (g_bTegri) {
				g_bVoroIntegrateCharge = true;
				g_bVoroIntegrateDipoleMoment = true;
				g_bVoroIntegrateTotalCurrent = true;
				g_bVoroIntegrateMagneticMoment = true;
				g_bCubeTimeDev = true;
				parseCoreCharges();
				((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 7;
			} else {
				eprintf("    Voronoi integration needs to be active. Use \"vori\" in the main menu.\n");
			}
		} else if (mode == 8) {
			setupMagneticDipoleRestartFile();
			((CMolecule *)g_oaMolecules[mol])->m_iMagneticDipoleMode = 8;
		} else {
			eprintf("    This is impossible.\n");
		}
	}
	mprintf("\n");
	
	if (g_bCubeTimeDev) {
		g_fBackgroundDensity = AskFloat("    Background density to improve PDE solver convergence in atomic units [1e-6] ", 1.0e-6);
		g_fPDEConvThresh = AskFloat("    Relative convergence threshold for PDE solver [1e-2] ", 0.01);
		g_iPDEMaxIter = AskInteger("    Maximum number of iterations in PDE solver [20] ", 20);
		mprintf("\n");
	}
	
	mprintf(WHITE, "<<< End of Magnetic Moment Definition <<<\n\n");
	
	magneticDipoleDefined = true;
}


void DipolGrimme(const char *s)
{
	FILE *a;
	char buf[256], *p, *q;
	const char *sep = " ,;\t";
	int i, ti, z, z3;
	CxDoubleArray *fa, *ptfa2, *ptfa3;
	double fx, fy, fz, tf;
	CFFT tfft;
	CAutoCorrelation *ac;
	CReorDyn *m_pRDyn;

	mprintf("\n    This analysis reads dipole vectors from a text file\n");
	mprintf("    and calculates the dipole vector ACF / IR spectrum.\n\n");
	mprintf("    The text file needs to contain the three Cartesian components\n");
	mprintf("    of the dipole vector (three numbers per line); one line equals one time step.\n\n");

	if (s == NULL)
	{
		eprintf("Please specify dipole text file as first command line parameter.\n\n");
		return;
	}

	mprintf("    If you are not sure what to enter, use the default values (just press return).\n\n");

	g_fTimestepLength = AskFloat("    Enter the physical time distance between successive trajectory frames in fs: [0.5] ",0.5);
	mprintf("\n");

	try { m_pRDyn = new CReorDyn(); } catch(...) { m_pRDyn = NULL; }
	if (m_pRDyn == NULL) NewException((double)sizeof(CReorDyn),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
//	m_pRDyn->m_sName = "dipole";

	try { m_pRDyn->m_pRDyn = new CDF(); } catch(...) { m_pRDyn->m_pRDyn = NULL; }
	if (m_pRDyn->m_pRDyn == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pRDyn->m_pRDyn->m_bLeft = true;
	mprintf(WHITE,"\n*** Vector Reorientation Dynamics ***\n\n");

	m_pRDyn->m_iDepth = AskUnsignedInteger("    Enter the correlation depth of the vector ACF (in trajectory frames): [4096] ",4096);

	mprintf("\n    This corresponds to a spectral resolution of %.4f cm^-1.\n",33356.41/g_fTimestepLength/2.0/m_pRDyn->m_iDepth);

	ti = CalcFFTSize(m_pRDyn->m_iDepth,false);
	if (m_pRDyn->m_iDepth != ti)
	{
		mprintf(WHITE,"\n    The next \"fast\" size for FFT is %d. Using this instead of %d as depth.\n",tfft.NextFastSize(m_pRDyn->m_iDepth),m_pRDyn->m_iDepth);
		m_pRDyn->m_iDepth = tfft.NextFastSize(m_pRDyn->m_iDepth);
	}

	try { m_pRDyn->m_pACF = new CACF(); } catch(...) { m_pRDyn->m_pACF = NULL; }
	if (m_pRDyn->m_pACF == NULL) NewException((double)sizeof(CACF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	m_pRDyn->m_pACF->m_iSize = m_pRDyn->m_iDepth;
	m_pRDyn->m_pACF->m_bSpectrum = true;

	m_pRDyn->m_pACF->m_bWindowFunction = AskYesNo("    Apply window function (Cos^2) to Autocorrelation function (y/n)? [yes] ",true);

	tf = 33356.41 / g_fTimestepLength / 2.0;
	mprintf("\n    A time step length of %.1f fs allows a spectral range up to %.1f cm^-1.\n\n",g_fTimestepLength,tf);
	m_pRDyn->m_pACF->m_fSpecWaveNumber = AskRangeFloat("    Calculate spectrum up to which wave number (cm^-1)? [%.1f cm^-1] ",0,tf,(tf<5000.0)?tf:5000.0,(tf<5000.0)?tf:5000.0);
	m_pRDyn->m_pACF->m_iMirror = AskUnsignedInteger("    No mirroring (0), short-time enhancing (1) or long-time enhancing (2)? [1] ",1);
	m_pRDyn->m_pACF->m_iZeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ",m_pRDyn->m_iDepth*3,m_pRDyn->m_iDepth*3);
	m_pRDyn->m_pACF->m_iZeroPadding0 = m_pRDyn->m_pACF->m_iZeroPadding;

	ti = CalcFFTSize(m_pRDyn->m_iDepth+m_pRDyn->m_pACF->m_iZeroPadding,false);
	if (m_pRDyn->m_iDepth+m_pRDyn->m_pACF->m_iZeroPadding != ti)
	{
		mprintf(WHITE,"\n    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n",ti,ti-m_pRDyn->m_iDepth);
		m_pRDyn->m_pACF->m_iZeroPadding = ti-m_pRDyn->m_iDepth;
	}

	mprintf("    Zero padding increases the spectral resolution to %.4f cm^-1.\n\n",33356.41/g_fTimestepLength/2.0/(m_pRDyn->m_iDepth+m_pRDyn->m_pACF->m_iZeroPadding));

	m_pRDyn->m_pACF->m_bACF_DB = AskYesNo("    Convert intensity axis of spectrum to decibel (y/n)? [no] ",false);
	m_pRDyn->m_pACF->Create();

	m_pRDyn->m_pRDyn->m_fMinVal = 0;
	m_pRDyn->m_pRDyn->m_fMaxVal = m_pRDyn->m_iDepth * g_fTimestepLength / 1000.0;
	m_pRDyn->m_pRDyn->m_iResolution = m_pRDyn->m_iDepth;
	m_pRDyn->m_pRDyn->SetLabelX("Tau / ps");
	m_pRDyn->m_pRDyn->SetLabelY("Vector autocorrelation");
	m_pRDyn->m_pRDyn->Create();

	mprintf("    Trying to open dipole file \"%s\"...\n",s);

	a = fopen(s,"rt");

	if (a == NULL)
	{
		eprintf("Error: Could not open file for reading.\n\n");
		return;
	}

	mprintf("\nReading dipole data:\n");

	try { fa = new CxDoubleArray("DipolGrimme():fa"); } catch(...) { fa = NULL; }
	if (fa == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	i = 0;
	while (!feof(a))
	{
		(void)!fgets(buf,256,a);
		if (feof(a))
		{
			mprintf("\n\nEnd of dipole file reached.");
			break;
		}
		buf[strlen(buf)-1] = 0;

		p = buf;

		while (strchr(sep,*p) != NULL)
			p++;
		q = p;
		while ((strchr(sep,*p) == NULL) && (*p != 0))
			p++;
		*p = 0;
		p++;
		fx = atof(q);

		while (strchr(sep,*p) != NULL)
			p++;
		q = p;
		while ((strchr(sep,*p) == NULL) && (*p != 0))
			p++;
		*p = 0;
		p++;
		fy = atof(q);

		while (strchr(sep,*p) != NULL)
			p++;
		q = p;
		while ((strchr(sep,*p) == NULL) && (*p != 0))
			p++;
		*p = 0;
		fz = atof(q);

//		mprintf("\n( %f | %f | %f)",fx,fy,fz);

		fa->Add(fx);
		fa->Add(fy);
		fa->Add(fz);
	
		if ((i % 50) == 0)
			mprintf("\n%6d ",i);
		mprintf(".");

		i++;
	}
	mprintf("\n\n%d time steps read.\n\n",i);

	fclose(a);

	mprintf("    Autocorrelating cached vectors via FFT...\n");

	try { ptfa2 = new CxDoubleArray("DipolGrimme():ptfa2"); } catch(...) { ptfa2 = NULL; }
	if (ptfa2 == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	ptfa2->SetSize(i);

	try { ptfa3 = new CxDoubleArray("DipolGrimme():ptfa3"); } catch(...) { ptfa3 = NULL; }
	if (ptfa3 == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	ptfa3->SetSize(i);

	try { ac = new CAutoCorrelation(); } catch(...) { ac = NULL; }
	if (ac == NULL) NewException((double)sizeof(CAutoCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	ac->Init(i,m_pRDyn->m_iDepth,true);

	/* X */
	for (z3=0;z3<(int)i;z3++)
		(*ptfa2)[z3] = (*fa)[z3*3];

	ac->AutoCorrelate(ptfa2,ptfa3);

	for (z3=0;z3<(int)m_pRDyn->m_iDepth;z3++)
		m_pRDyn->m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3]);

	/* Y */
	for (z3=0;z3<(int)i;z3++)
		(*ptfa2)[z3] = (*fa)[z3*3+1];

	ac->AutoCorrelate(ptfa2,ptfa3);

	for (z3=0;z3<(int)m_pRDyn->m_iDepth;z3++)
		m_pRDyn->m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3]);

	/* Z */
	for (z3=0;z3<(int)i;z3++)
		(*ptfa2)[z3] = (*fa)[z3*3+2];

	ac->AutoCorrelate(ptfa2,ptfa3);

	for (z3=0;z3<(int)m_pRDyn->m_iDepth;z3++)
	{
		m_pRDyn->m_pRDyn->AddToBin_Int(z3,(*ptfa3)[z3]);
		m_pRDyn->m_pRDyn->m_fBinEntries += (double)(i-z3) - 3.0;
	}

	delete ac;
	delete ptfa2;
	delete ptfa3;

	mprintf("    Finished.\n\n");
	mprintf("    %.0f bin entries.\n",m_pRDyn->m_pRDyn->m_fBinEntries);
	mprintf("    Starting value is %f.\n",m_pRDyn->m_pRDyn->m_pBin[0]);
	m_pRDyn->m_pRDyn->MultiplyBin(1.0/m_pRDyn->m_pRDyn->m_pBin[0]);
	sprintf(buf,"rdyn_dipole.csv");
	mprintf("    Saving result as %s ...\n",buf);
	m_pRDyn->m_pRDyn->Write("",buf,"",false);

	mprintf("\n    Creating reorientation spectrum:\n");

	for (z=0;z<m_pRDyn->m_iDepth;z++)
 		m_pRDyn->m_pACF->m_pData[z] = m_pRDyn->m_pRDyn->m_pBin[z];

	if (m_pRDyn->m_pACF->m_iMirror != 0)
	{
		mprintf("    Mirroring ACF...\n");
		m_pRDyn->m_pACF->Mirror(m_pRDyn->m_pACF->m_iMirror);
		sprintf(buf,"rdyn_dipole.mirrored.csv");
		mprintf("    Saving mirrored ACF as %s ...\n",buf);
		m_pRDyn->m_pACF->WriteACF("",buf,"");
	}

	if (m_pRDyn->m_pACF->m_bWindowFunction)
	{
		mprintf("    Applying window function to ACF...\n");
		m_pRDyn->m_pACF->Window();
		sprintf(buf,"rdyn_dipole.windowed.csv");
		mprintf("    Saving windowed ACF as %s ...\n",buf);
		m_pRDyn->m_pACF->WriteACF("",buf,"");
	}

	mprintf("    Performing Fourier transformation...\n");

	try { g_pFFT = new CFFT(); } catch(...) { g_pFFT = NULL; }
	if (g_pFFT == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	g_pFFT->PrepareFFT_C2C(m_pRDyn->m_pACF->m_iSize+m_pRDyn->m_pACF->m_iZeroPadding);

	m_pRDyn->m_pACF->Transform(g_pFFT);
	delete g_pFFT;
	m_pRDyn->m_pACF->m_pSpectrum->SetMaxRWL(1E7/299.792/g_fTimestepLength);
	if (m_pRDyn->m_pACF->m_bACF_DB)
	{
		mprintf("    Normalizing spectrum to decibel...\n");
		m_pRDyn->m_pACF->m_pSpectrum->MakeDB();
	}

	sprintf(buf,"spectrum_dipole.csv");
	mprintf("    Saving spectrum as %s ...\n",buf);
	m_pRDyn->m_pACF->m_pSpectrum->Write("",buf,"");

	mprintf("\n\n*** Finished ***\n\n");

	WriteCredits();
}


void InitGlobalVars()
{
	g_bUseBQB = false;
	g_bFlipCellVectorX = false;
	g_bFlipCellVectorY = false;
	g_bFlipCellVectorZ = false;
	g_iRemoveCOMFixAtom = -1;
	g_bReadLAMMPSCharges = false;
	g_bLAMMPSCharge = false;
	g_bLAMMPSForce = false;
	g_iVoroPrintLevel = 0;
	g_bSDFVoro = false;
	g_bDoubleBox = false;
	g_bWriteOrtho = false;
	g_bBoxNonOrtho = false;
	g_bFoundNonOrtho = false;
	g_bWriteInputOrder = false;
	g_bShowCredits = false;
	g_bPairMSD = false;
	g_bStreamInput = false;
	g_bReactive = false;
	g_sConfFile = NULL;
	g_bShowConf = false;
	g_bWriteConf = false;
	g_bCheckWrite = true;
	g_bSaxonize = false;
	g_bUnknownElements = false;
	g_iStartTime = 0;
	g_bSMode = false;
	g_bNoColor = false;
	g_iColorIntensity = 1;
	g_bNPT = true; // Nur anfaenglich
//	g_sNPTFile[0] = 0;
	g_sNPTFile = "";
	g_fNPTFile = NULL;
	g_fBoxX = 0;
	g_fBoxY = 0;
	g_fBoxZ = 0;
	g_bVerbose = false;
	g_bClusterAnalysis = false;
	g_bTDDF = false;
	g_pTDDFEngine = NULL;
	g_bGeoDens = false;
	g_pGeoDensEngine = NULL;
	g_pAggrTopoEngine = NULL;
	g_bAggrTopo = false;
	g_bSankey = false;
	g_sSankeyFile[0] = 0;
	g_bNeedMoleculeWrap = false;
	g_bDipoleDefined = false;
	g_bDipolGrimme = false;
	g_bRamanFromPolarizability = false;
	g_fTimestepLength = 0;
	g_sInputTraj = NULL;
	g_sInputFile = NULL;
	g_bKeepOriginalCoords = false;
	g_bAbortAnalysis = true;
	g_pTempTimestep = NULL;
	g_pTempTimestepSnap = NULL;
	g_bMiddleAvg = true;
	g_bVFHisto = false;
	g_bUseVelocities = false;
	g_bUseForces = false;
	g_bWriteAtomwise = false;
	g_fBondFactor = 1.15;
	g_fVelPercentage = 0.95;
	g_fForcePercentage = 0.95;
	g_fPos = NULL;
	g_fVel = NULL;
	g_fForce = NULL;
	g_iRefSystemDim = 0;
	g_iFixMol = -1;
	g_fInputFile = NULL;
	g_fInput = NULL;
	g_fSaveCondFile = NULL;
	g_bSaveCondSnapshot = false;
	g_bScanVelocities = false;
	g_bSilentProgress = false;
	g_bCreateRevSDF = false;
	g_bDeriv = false;
	g_iStrideDetect = -1;
	g_bGlobalIR = false;
	g_pGlobalIR = NULL;
	g_bLMFitSilent = false;
	g_iLMMaxIter = 200;
	g_bXYZ4thCol = false;
	g_bXYZComment3Numbers = false;
	g_bXYZComment6Numbers = false;
	g_bXYZComment9Numbers = false;
	g_bReadChargesFrom4thXYZ = false;
	g_bEnvWriteDetailedInfo = false;
	g_bEnvSortNb = false;
	g_bEnvDisableSortNb = false;
	g_bSDFMap = false;
	g_iSDFMapSmoothGrade = 0;
	g_bVoroSilent = false;
	g_pOrderEngine = NULL;
	g_bMinimizeChargeVar2 = false;
	g_pContactMatrix = NULL;
	#ifdef NEW_CHARGEVAR
		g_pChargeVar2 = NULL;
	#endif
	g_bReRa = false;
	g_bQuadrupoleKeepTrace = false;
	g_bVoronoiMoments = true;
	g_iWannierAtomType = -1;
	g_iRefMolNum = 0;
	g_bProcAlternativeLabels = false;
	g_bProcDumpAwkScript = false;
	g_bProcWriteComments = false;
	g_bProcCellComment = false;
	g_bInterWarning = false;
	g_fMolIntegralFile = NULL;
	g_bVolumetricData = false;
	g_bElMagProperties = false;
	g_bReadVelocity = false;
	g_fVelocityConversion = 0;
	g_bStrideParsed = false;
	g_pFixedPlProj = NULL;
	g_sControlFile = NULL;
	g_bControlRun = false;
	g_pControlDB = NULL;
	g_iControlRT = -1;
	g_pReactEngine = NULL;
	g_pBQBInterface = NULL;
	g_bOrder = false;
	g_pDCDLabelsTS = NULL;
	g_bDCDSkipHeader = false;
	g_bDCDFirstFrame = false;

	
	g_bTDDF = false;
	g_bContactMatrix = false;
	g_bFixedPlProj = false;
	g_bROA = false;
	g_bSFac = false;
	g_bVoid = false;
	g_bTegri = false;
	g_bVoro = false;
	g_bDomA = false;
	g_bPlProj = false;
	g_bPlDF = false;
	g_bLiDF = false;
	g_bDens = false;
	g_bCond = false;
	g_bSDF = false;
	g_bRDF = false;
	g_bVDF = false;
	g_bFDF = false;
	g_bADF = false;
	g_bAvg = false;
	g_bSaveVelForce = false;
	g_bSaveRefEnv = false;
	g_bRefEnvCenter = false;
	g_bSaveJustTraj = false;
	g_bVHDF = false;
	g_bCutCluster = false;
	g_bNbAnalysis = false;
	g_bVACF = false;
	g_bDipACF = false;
	g_bDDF = false;
	g_bMSD = false;
	g_bDipole = false;
	g_bDipDF = false;
	g_bCDF = false;
	g_bAggregation = false;
	g_bDDisp = false; 
	g_bDACF = false; 
	g_bDLDF = false; 
	g_bRDyn = false;
	g_bNbExchange = false;
	g_bRaman = false;
	g_bIRSpec = false;
	g_bPowerSpec = false;
	g_bRegionAnalysis = false;
	g_bVACFNew = false;
	g_bMolRecNoDomain = false;
	g_bKuehneLong = false;
	g_bSaveTrajPattern = false;
	g_iSaveTrajPatternPos = 0;
	g_iCellVectorFileColumns = -1;
	g_iDCDAtomCount = 0;
}


CCrossCorrelation::CCrossCorrelation()
{
	m_iInput = 0;
	m_iDepth = 0;
	m_pFFT = NULL;
	m_pFFT2 = NULL;
}


CCrossCorrelation::~CCrossCorrelation()
{
}


void CCrossCorrelation::Init(int input, int depth, bool fft)
{
	m_iInput = input;
	m_iDepth = depth;
	if (fft)
	{
		m_bFFT = true;
		m_iFFTSize = CalcFFTSize(input,true);

		try { m_pFFT = new CFFT(); } catch(...) { m_pFFT = NULL; }
		if (m_pFFT == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pFFT->PrepareFFT_C2C(2*m_iFFTSize);

		try { m_pFFT2 = new CFFT(); } catch(...) { m_pFFT2 = NULL; }
		if (m_pFFT2 == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pFFT2->PrepareFFT_C2C(2*m_iFFTSize);

		try { m_pFFTback = new CFFT(); } catch(...) { m_pFFTback = NULL; }
		if (m_pFFTback == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m_pFFTback->PrepareInverseFFT_C2C(2*m_iFFTSize);
	} else
	{
		m_bFFT = false;
	}
}


void CCrossCorrelation::CrossCorrelate(CxDoubleArray *inp1, CxDoubleArray *inp2, CxDoubleArray *outp)
{
	int z, z2;
	double tf;

	outp->SetSize(m_iDepth);

	if (m_bFFT)
	{
		for (z=0;z<m_iInput;z++)
		{
			m_pFFT->m_pInput[z*2] = (*inp1)[z];
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		for (z=m_iInput;z<2*m_iFFTSize;z++)
		{
			m_pFFT->m_pInput[z*2] = 0;
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		m_pFFT->DoFFT();

		for (z=0;z<m_iInput;z++)
		{
			m_pFFT2->m_pInput[z*2] = (*inp2)[z];
			m_pFFT2->m_pInput[z*2+1] = 0;
		}
		for (z=m_iInput;z<2*m_iFFTSize;z++)
		{
			m_pFFT2->m_pInput[z*2] = 0;
			m_pFFT2->m_pInput[z*2+1] = 0;
		}
		m_pFFT2->DoFFT();

		for (z=0;z<m_iFFTSize*2;z++)
		{
			// a1*a2 + b1*b2
			m_pFFTback->m_pInput[z*2]   = (m_pFFT->m_pOutput[z*2]*m_pFFT2->m_pOutput[z*2]    + m_pFFT->m_pOutput[z*2+1]*m_pFFT2->m_pOutput[z*2+1]);
			// a2*b1 - a1*b2
			m_pFFTback->m_pInput[z*2+1] = (-m_pFFT2->m_pOutput[z*2]*m_pFFT->m_pOutput[z*2+1] + m_pFFT->m_pOutput[z*2]*m_pFFT2->m_pOutput[z*2+1]);
		}
		m_pFFTback->DoFFT();

		for (z=0;z<m_iDepth;z++)
			(*outp)[z] = m_pFFTback->m_pOutput[2*z] / m_iFFTSize / 2.0 / ((double)m_iInput - z);
	} else
	{
		for (z=0;z<m_iDepth;z++) // Tau
		{
			tf = 0;
			for (z2=0;z2<m_iInput-z;z2++)
				tf += (*inp1)[z2] * (*inp2)[z2+z];
			(*outp)[z] = tf / (double)(m_iInput-z);
		}
	}
}


void CCrossCorrelation::CrossCorrelate(std::vector<double> &inp1, std::vector<double> &inp2, std::vector<double> &outp) {

	int z, z2;
	double tf;

	outp.resize(m_iDepth);

	if (m_bFFT) {

		for (z=0;z<m_iInput;z++) {
			m_pFFT->m_pInput[z*2] = inp1[z];
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		for (z=m_iInput;z<2*m_iFFTSize;z++) {
			m_pFFT->m_pInput[z*2] = 0;
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		m_pFFT->DoFFT();

		for (z=0;z<m_iInput;z++) {
			m_pFFT2->m_pInput[z*2] = inp2[z];
			m_pFFT2->m_pInput[z*2+1] = 0;
		}
		for (z=m_iInput;z<2*m_iFFTSize;z++) {
			m_pFFT2->m_pInput[z*2] = 0;
			m_pFFT2->m_pInput[z*2+1] = 0;
		}
		m_pFFT2->DoFFT();

		for (z=0;z<m_iFFTSize*2;z++) {
			// a1*a2 + b1*b2
			m_pFFTback->m_pInput[z*2]   = (m_pFFT->m_pOutput[z*2]*m_pFFT2->m_pOutput[z*2]    + m_pFFT->m_pOutput[z*2+1]*m_pFFT2->m_pOutput[z*2+1]);
			// a2*b1 - a1*b2
			m_pFFTback->m_pInput[z*2+1] = (-m_pFFT2->m_pOutput[z*2]*m_pFFT->m_pOutput[z*2+1] + m_pFFT->m_pOutput[z*2]*m_pFFT2->m_pOutput[z*2+1]);
		}
		m_pFFTback->DoFFT();

		for (z=0;z<m_iDepth;z++)
			outp[z] = m_pFFTback->m_pOutput[2*z] / m_iFFTSize / 2.0 / ((double)m_iInput - z);

	} else {

		for (z=0;z<m_iDepth;z++) { // Tau
			tf = 0;
			for (z2=0;z2<m_iInput-z;z2++)
				tf += inp1[z2] * inp2[z2+z];
			outp[z] = tf / (double)(m_iInput-z);
		}
	}
}


void CCrossCorrelation::CrossCorrelateSymmetric(CxDoubleArray *inp1, CxDoubleArray *inp2, CxDoubleArray *outp) {

	int z, z2;
	double tf;

	outp->SetSize(m_iDepth);

	if (m_bFFT)
	{
		for (z=0;z<m_iInput;z++)
		{
			m_pFFT->m_pInput[z*2] = (*inp1)[z];
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		for (z=m_iInput;z<2*m_iFFTSize;z++)
		{
			m_pFFT->m_pInput[z*2] = 0;
			m_pFFT->m_pInput[z*2+1] = 0;
		}
		m_pFFT->DoFFT();

		for (z=0;z<m_iInput;z++)
		{
			m_pFFT2->m_pInput[z*2] = (*inp2)[z];
			m_pFFT2->m_pInput[z*2+1] = 0;
		}
		for (z=m_iInput;z<2*m_iFFTSize;z++)
		{
			m_pFFT2->m_pInput[z*2] = 0;
			m_pFFT2->m_pInput[z*2+1] = 0;
		}
		m_pFFT2->DoFFT();

		for (z=0;z<m_iFFTSize*2;z++)
		{
			// a1*a2 + b1*b2
			m_pFFTback->m_pInput[z*2]   = (m_pFFT->m_pOutput[z*2]*m_pFFT2->m_pOutput[z*2]    + m_pFFT->m_pOutput[z*2+1]*m_pFFT2->m_pOutput[z*2+1]);
			// a2*b1 - a1*b2
			m_pFFTback->m_pInput[z*2+1] = (-m_pFFT2->m_pOutput[z*2]*m_pFFT->m_pOutput[z*2+1] + m_pFFT->m_pOutput[z*2]*m_pFFT2->m_pOutput[z*2+1]);
		}
		m_pFFTback->DoFFT();

		for (z=0;z<m_iDepth;z++)
			(*outp)[z] = m_pFFTback->m_pOutput[2*z] / m_iFFTSize / 2.0 / ((double)m_iInput - z);

		for (z=0;z<m_iFFTSize*2;z++)
		{
			// a2*a1 + b2*b1
			m_pFFTback->m_pInput[z*2]   = (m_pFFT2->m_pOutput[z*2]*m_pFFT->m_pOutput[z*2]    + m_pFFT2->m_pOutput[z*2+1]*m_pFFT->m_pOutput[z*2+1]);
			// a1*b2 - a2*b1
			m_pFFTback->m_pInput[z*2+1] = (-m_pFFT->m_pOutput[z*2]*m_pFFT2->m_pOutput[z*2+1] + m_pFFT2->m_pOutput[z*2]*m_pFFT->m_pOutput[z*2+1]);
		}
		m_pFFTback->DoFFT();

		for (z=0;z<m_iDepth;z++)
			(*outp)[z] -= m_pFFTback->m_pOutput[2*z] / m_iFFTSize / 2.0 / ((double)m_iInput - z);
	} else
	{
		for (z=0;z<m_iDepth;z++) // Tau
		{
			tf = 0;
			for (z2=0;z2<m_iInput-z;z2++)
				tf += (*inp1)[z2] * (*inp2)[z2+z];
			(*outp)[z] = tf / (double)(m_iInput-z);
			tf = 0;
			for (z2=0;z2<m_iInput-z;z2++)
				tf += (*inp2)[z2] * (*inp1)[z2+z];
			(*outp)[z] -= tf / (double)(m_iInput-z);
		}
	}
}


void ParseVoronoiRadii()
{
	int i, z, z2, z3, z4;
	double *f, tf;
	CMolecule *m;

	if (g_faVoronoiRadii.GetSize() != 0)
		return; // Already parsed

	mprintf("\n    When performing the Voronoi decomposition, radii may be assigned to the atoms (\"radical Voronoi tesselation\").\n\n");

	i = AskRangeInteger("    Do not assign radii (0), use covalent radii (1), use Van-der-Waals radii (2), or specify other radii (3)? [2] ",0,3,2);

	if (i == 0)
	{
		g_faVoronoiRadii.SetSize(g_iGesAtomCount);

		for (z=0;z<g_iGesAtomCount;z++)
			g_faVoronoiRadii[z] = 0.5;
	}

	if (i == 1)
	{
		mprintf("\n    Using the following atom radii:\n");
		for (z=0;z<g_oaAtoms.GetSize();z++)
		{
			if (z == g_iVirtAtomType)
				continue;

			mprintf("      %2s: %5.1f pm.\n",(const char*)((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius);
			
			if (((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius < 1.0)
				mprintf(RED, "Warning: The radius of %s is very small.\n", (const char*)((CAtom*)g_oaAtoms[z])->m_sName);
		}

		g_faVoronoiRadii.SetSize(g_iGesAtomCount);

		for (z=0;z<g_iGesAtomCount;z++)
			g_faVoronoiRadii[z] = ((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_fRadius;
	}

	if (i == 2)
	{
		mprintf("\n    Using the following atom radii:\n");
		for (z=0;z<g_oaAtoms.GetSize();z++)
		{
			if (z == g_iVirtAtomType)
				continue;

			mprintf("      - %2s: %5.1f pm.\n",(const char*)((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[z])->m_pElement->m_fVdWRadius);
			
			if (((CAtom*)g_oaAtoms[z])->m_pElement->m_fVdWRadius < 1.0)
				mprintf(RED, "Warning: The radius of %s is very small.\n", (const char*)((CAtom*)g_oaAtoms[z])->m_sName);
		}

		g_faVoronoiRadii.SetSize(g_iGesAtomCount);

		for (z=0;z<g_iGesAtomCount;z++)
			g_faVoronoiRadii[z] = ((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_fVdWRadius;
	}

	if (i == 3)
	{
		mprintf("    Enter 0 as radius to exclude atoms from the Voronoi decomposition.\n\n");

		g_faVoronoiRadii.SetSize(g_iGesAtomCount);

		if (AskYesNo("    Assign atom radii per element type (y) or per atom (n)? [yes] ",true))
		{
			mprintf("\n");
			f = new double[g_oaAtoms.GetSize()];

			for (z=0;z<g_oaAtoms.GetSize();z++)
			{
				if (z == g_iVirtAtomType)
					continue;
		
				f[z] = AskFloat_ND("      Which radius to use for element type %s (in pm)? ",(const char*)((CAtom*)g_oaAtoms[z])->m_sName);
			}

			for (z=0;z<g_iGesAtomCount;z++)
				g_faVoronoiRadii[z] = f[g_waAtomRealElement[z]];

			delete[] f;
		} else
		{
			mprintf("\n");
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				m = (CMolecule*)g_oaMolecules[z];
				mprintf("    * Molecule %s\n\n",m->m_sName);
				for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
				{
					if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
						continue;
					for (z3=0;z3<m->m_waAtomCount[z2];z3++)
					{
						tf = AskFloat_ND("      Radius for %2s%-2d (in pm): ",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1);
						for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
							g_faVoronoiRadii[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]])->m_oaAtomOffset[z2])->GetAt(z3)] = tf;
					}
				}
				mprintf("\n");
			}
		}
	}

	if (i == 0)
		mprintf("\n    Not using Voronoi radii.\n\n");
	else
		mprintf("\n    Voronoi radii defined.\n\n");


	mprintf("    The Voronoi integration can have issues with highly symmetric structures.\n");
	mprintf("    It is recommended to add a tiny amount of random noise (\"jitter\") to the\n");
	mprintf("    atom coordinates to break the symmetry. The jitter will be switched on\n");
	mprintf("    by default. Switch it off for backward compatibility with earlier TRAVIS versions.\n\n");

	glv_bJitter = AskYesNo("    Use position jitter to mitigate symmetry issues (y/n)? [yes] ",true);

	if (glv_bJitter) {

		if (g_bAdvanced2) {
			// The jitter amplitude variable is in Angstrom internally...
			glv_fJitterAmplitude = AskFloat("    Enter jitter amplitude (in pm): [0.1] ",0.1) * 0.01;
			glv_iJitterSeed = AskUnsignedInteger("    Enter jitter random seed: [0] ",0);
		} else {
			mprintf("\n");
			mprintf("    Using a jitter amplitude of 0.1 pm and a random seed of 0.\n");
			mprintf("    Switch on the advanced mode to modify these parameters.\n");
		}

		mprintf("\n");
		mprintf("    Jitter switched on.\n");

	} else {

		mprintf("\n");
		mprintf("    Jitter switched off.\n");
	}

	mprintf("\n");
}


void DumpNonOrthoCellData()
{
	double tf;
	CxDVector3 veca, vecb, vecc;

	veca = CxDVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
	vecb = CxDVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
	vecc = CxDVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

	mprintf(WHITE,"    *** Non-orthorhombic cell geometry ***\n");
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *"); mprintf("   Cell vector A: ( %10.4f | %10.4f | %10.4f ), Length %10.3f pm.\n",g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2),veca.GetLength());
	mprintf(WHITE,"    *"); mprintf("   Cell vector B: ( %10.4f | %10.4f | %10.4f ), Length %10.3f pm.\n",g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2),vecb.GetLength());
	mprintf(WHITE,"    *"); mprintf("   Cell vector C: ( %10.4f | %10.4f | %10.4f ), Length %10.3f pm.\n",g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2),vecc.GetLength());
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *"); mprintf("   Cell angle Alpha (between B and C): %8.3f degree.\n",g_fBoxAngleA);
	mprintf(WHITE,"    *"); mprintf("   Cell angle Beta  (between A and C): %8.3f degree.\n",g_fBoxAngleB);
	mprintf(WHITE,"    *"); mprintf("   Cell angle Gamma (between A and B): %8.3f degree.\n",g_fBoxAngleC);
	mprintf(WHITE,"    *\n");

	tf = DotP(CrossP(veca,vecb),vecc)/1000000.0;
	mprintf(WHITE,"    *"); mprintf("   Cell volume:  %.3f Angstrom^3 = %.6f nm^3.\n",tf,tf/1000.0);
	g_fBoxVolume = tf;
	mprintf(WHITE,"    *"); mprintf("   Cell density: %.6f g/cm^3.\n",pow3(GuessBoxSize())/tf/1000000.0);
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *"); mprintf("   Orthogonal bounding box of the cell:\n");
	mprintf(WHITE,"    *"); mprintf("     dX = %10.4f pm,  dY = %10.4f pm,  dZ = %10.4f pm.\n",g_fBoxX,g_fBoxY,g_fBoxZ);
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *"); mprintf("   Minimal periodic diameter (smallest distance between atom and its periodic image):\n");
	mprintf(WHITE,"    *"); mprintf("     Min. Diam. = %10.4f pm.\n",g_fBoxMinDiam);
	mprintf(WHITE,"    *"); mprintf("     ( dA = %10.4f pm,  dB = %10.4f pm,  dC = %10.4f pm )\n",g_fBoxMinDiamA,g_fBoxMinDiamB,g_fBoxMinDiamC);
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *"); mprintf("   Determinant of the forward transformation matrix is  %12G.\n",g_mBoxFromOrtho.Det());
	mprintf(WHITE,"    *"); mprintf("   Determinant of the backward transformation matrix is %12G.\n",g_mBoxToOrtho.Det());
	mprintf(WHITE,"    *\n");
	mprintf(WHITE,"    *** End of non-orthorhombic cell geometry ***\n\n");
}


void ExtractXYZCellGeometry6(const char *s)
{
	int z;
	char buf[64];
	const char *p, *q;
	double data[6];
	CxDVector3 veca, vecb, vecc;

	z = 0;
	p = s;
	while (*p != 0)
	{
//		mprintf("# (A) z=%d \"%s\".\n",z,p);
		while (strchr("0123456789+-.Ee",*p) == NULL)
			p++;
		q = p;
//		mprintf("# (B) z=%d \"%s\".\n",z,p);
		while ((strchr("0123456789+-.Ee",*q) != NULL) && (*q != 0))
			q++;
		if (*(q-1) == 0)
			q--;
//		mprintf("# (C) z=%d q=\"%s\" p=\"%s\".\n",z,q,p);
		if (q > p)
		{
//			mprintf("# q>p, q-p = %d\n",q-p);
//			mprintf("# p=\"%s\".\n",p);
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
//			mprintf("# buf = \"%s\".\n",buf);
			if (IsValidFloat(buf))
			{
				data[z] = atof(buf);
				z++;
			}
		}// else mprintf("# q <= p\n");
//		mprintf("# z=%d.\n",z);
		if (*q == 0)
			break;
//		mprintf("# X.\n");
		if (z == 6)
			break;
//		mprintf("# Y.\n");
		p = q;
	}
	if (z < 6)
	{
		eprintf("\nError: ExtractXYZCellGeometry6(): Comment line of XYZ trajectory does not contain 6 numbers.\n");
		eprintf("  (Comment line is \"%s\").   ",s);
		return;
	}

	g_fBoxX = data[0] * 100.0; // Angstrom --> pm
	g_fBoxY = data[1] * 100.0;
	g_fBoxZ = data[2] * 100.0;
	g_fBoxAngleA = data[3];
	g_fBoxAngleB = data[4];
	g_fBoxAngleC = data[5];
//	mprintf("    Computing vectors...\n");

	g_mBoxFromOrtho(0,0) = g_fBoxX;
	g_mBoxFromOrtho(0,1) = 0;
	g_mBoxFromOrtho(0,2) = 0;

	g_mBoxFromOrtho(1,0) = g_fBoxY*cos(g_fBoxAngleC*Pi/180.0);
	g_mBoxFromOrtho(1,1) = g_fBoxY*sin(g_fBoxAngleC*Pi/180.0);
	g_mBoxFromOrtho(1,2) = 0;

	g_mBoxFromOrtho(2,0) = g_fBoxZ*cos(g_fBoxAngleB*Pi/180.0);
	g_mBoxFromOrtho(2,1) = (-(g_mBoxFromOrtho(1,0)*g_fBoxZ*cos(g_fBoxAngleB*Pi/180.0)) + g_fBoxY*g_fBoxZ*cos(g_fBoxAngleA*Pi/180.0))/g_mBoxFromOrtho(1,1);
	g_mBoxFromOrtho(2,2) = sqrt(-((pow2(g_fBoxZ)*(pow2(g_mBoxFromOrtho(1,1))*(-1 + pow2(cos(g_fBoxAngleB*Pi/180.0))) + pow2(g_mBoxFromOrtho(1,0)*cos(g_fBoxAngleB*Pi/180.0) - g_fBoxY*cos(g_fBoxAngleA*Pi/180.0))))/pow2(g_mBoxFromOrtho(1,1))));

	if (g_bDoubleBox)
	{
		g_mBoxFromOrtho(0,0) *= g_iDoubleBoxX;
		g_mBoxFromOrtho(1,0) *= g_iDoubleBoxX;
		g_mBoxFromOrtho(2,0) *= g_iDoubleBoxX;

		g_mBoxFromOrtho(0,1) *= g_iDoubleBoxY;
		g_mBoxFromOrtho(1,1) *= g_iDoubleBoxY;
		g_mBoxFromOrtho(2,1) *= g_iDoubleBoxY;

		g_mBoxFromOrtho(0,2) *= g_iDoubleBoxZ;
		g_mBoxFromOrtho(1,2) *= g_iDoubleBoxZ;
		g_mBoxFromOrtho(2,2) *= g_iDoubleBoxZ;
	}

//	mprintf("    Recalculating angles from vectors...\n");

	veca = CxDVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
	vecb = CxDVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
	vecc = CxDVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

	g_fBoxAngleA = acos(DotP(vecb,vecc) / vecb.GetLength() / vecc.GetLength()) * 180.0 / Pi;
	g_fBoxAngleB = acos(DotP(veca,vecc) / veca.GetLength() / vecc.GetLength()) * 180.0 / Pi;
	g_fBoxAngleC = acos(DotP(veca,vecb) / veca.GetLength() / vecb.GetLength()) * 180.0 / Pi;

//	mprintf("    Computing inverse of transformation matrix...\n\n");
	g_mBoxToOrtho = CxDMatrix3(g_mBoxFromOrtho);
	if (!g_mBoxToOrtho.Invert()) {
		eprintf("ExtractXYZCellGeometry6(): Error: Encountered singular cell matrix (cell volume is zero).\n");
		abort();
	}

	// Orthogonal bounding box
	g_fBoxX = fabs(g_mBoxFromOrtho(0,0)) + fabs(g_mBoxFromOrtho(1,0)) + fabs(g_mBoxFromOrtho(2,0));
	g_fBoxY = fabs(g_mBoxFromOrtho(0,1)) + fabs(g_mBoxFromOrtho(1,1)) + fabs(g_mBoxFromOrtho(2,1));
	g_fBoxZ = fabs(g_mBoxFromOrtho(0,2)) + fabs(g_mBoxFromOrtho(1,2)) + fabs(g_mBoxFromOrtho(2,2));

	veca = CxDVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
	vecb = CxDVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
	vecc = CxDVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

	// Minimal diameters
	g_fBoxMinDiamA = fabs(DotP(veca,Normalize(CrossP(vecb,vecc))));
	g_fBoxMinDiamB = fabs(DotP(vecb,Normalize(CrossP(veca,vecc))));
	g_fBoxMinDiamC = fabs(DotP(vecc,Normalize(CrossP(veca,vecb))));
	g_fBoxMinDiam = MIN3(g_fBoxMinDiamA,g_fBoxMinDiamB,g_fBoxMinDiamC);
}



void ExtractXYZCellGeometry9(const char *s) {

	int z;
	char buf[64];
	const char *p, *q;
	double data[9];
	CxDVector3 veca, vecb, vecc;


	z = 0;
	p = s;
	while (*p != 0) {

		while (strchr("0123456789+-.Ee",*p) == NULL)
			p++;
		q = p;

		while ((strchr("0123456789+-.Ee",*q) != NULL) && (*q != 0))
			q++;
		if (*(q-1) == 0)
			q--;

		if (q > p) {
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
			if (IsValidFloat(buf)) {
				data[z] = atof(buf);
				z++;
			}
		}
		if (*q == 0)
			break;
		if (z == 9)
			break;
		p = q;
	}
	if (z < 9) {
		eprintf("\nError: ExtractXYZCellGeometry9(): Comment line of trajectory does not contain 9 numbers.\n");
		eprintf("  (Comment line is \"%s\").   ",s);
		return;
	}

	g_mBoxFromOrtho(0,0) = data[0] * 100.0; // Angstrom --> pm
	g_mBoxFromOrtho(0,1) = data[1] * 100.0;
	g_mBoxFromOrtho(0,2) = data[2] * 100.0;

	g_mBoxFromOrtho(1,0) = data[3] * 100.0;
	g_mBoxFromOrtho(1,1) = data[4] * 100.0;
	g_mBoxFromOrtho(1,2) = data[5] * 100.0;

	g_mBoxFromOrtho(2,0) = data[6] * 100.0;
	g_mBoxFromOrtho(2,1) = data[7] * 100.0;
	g_mBoxFromOrtho(2,2) = data[8] * 100.0;

	if (g_bDoubleBox) {

		g_mBoxFromOrtho(0,0) *= g_iDoubleBoxX;
		g_mBoxFromOrtho(1,0) *= g_iDoubleBoxX;
		g_mBoxFromOrtho(2,0) *= g_iDoubleBoxX;

		g_mBoxFromOrtho(0,1) *= g_iDoubleBoxY;
		g_mBoxFromOrtho(1,1) *= g_iDoubleBoxY;
		g_mBoxFromOrtho(2,1) *= g_iDoubleBoxY;

		g_mBoxFromOrtho(0,2) *= g_iDoubleBoxZ;
		g_mBoxFromOrtho(1,2) *= g_iDoubleBoxZ;
		g_mBoxFromOrtho(2,2) *= g_iDoubleBoxZ;
	}

	veca = CxDVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
	vecb = CxDVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
	vecc = CxDVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

	g_fBoxAngleA = acos(DotP(vecb,vecc) / vecb.GetLength() / vecc.GetLength()) * 180.0 / Pi;
	g_fBoxAngleB = acos(DotP(veca,vecc) / veca.GetLength() / vecc.GetLength()) * 180.0 / Pi;
	g_fBoxAngleC = acos(DotP(veca,vecb) / veca.GetLength() / vecb.GetLength()) * 180.0 / Pi;

	g_mBoxToOrtho = CxDMatrix3(g_mBoxFromOrtho);
	if (!g_mBoxToOrtho.Invert()) {
		eprintf("ExtractXYZCellGeometry9(): Error: Encountered singular cell matrix (cell volume is zero).\n");
		abort();
	}

	// Orthogonal bounding box
	g_fBoxX = fabs(g_mBoxFromOrtho(0,0)) + fabs(g_mBoxFromOrtho(1,0)) + fabs(g_mBoxFromOrtho(2,0));
	g_fBoxY = fabs(g_mBoxFromOrtho(0,1)) + fabs(g_mBoxFromOrtho(1,1)) + fabs(g_mBoxFromOrtho(2,1));
	g_fBoxZ = fabs(g_mBoxFromOrtho(0,2)) + fabs(g_mBoxFromOrtho(1,2)) + fabs(g_mBoxFromOrtho(2,2));

	veca = CxDVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
	vecb = CxDVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
	vecc = CxDVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

	// Minimal diameters
	g_fBoxMinDiamA = fabs(DotP(veca,Normalize(CrossP(vecb,vecc))));
	g_fBoxMinDiamB = fabs(DotP(vecb,Normalize(CrossP(veca,vecc))));
	g_fBoxMinDiamC = fabs(DotP(vecc,Normalize(CrossP(veca,vecb))));
	g_fBoxMinDiam = MIN3(g_fBoxMinDiamA,g_fBoxMinDiamB,g_fBoxMinDiamC);
}



void ExtractXYZCellGeometry3(const char *s)
{
	int z;
	char buf[64];
	const char *p, *q;
	double data[6];
	CxDVector3 veca, vecb, vecc;

	if (strspn(s,"0123456789+-.Ee \r\n\t") != strlen(s)){
		eprintf("\nError: ExtractXYZCellGeometry3(): Comment line of XYZ trajectory contains non-numeric characters.\n");
		eprintf("  (Comment line is \"%s\").    ",s);
		return;
	}

	z = 0;
	p = s;
	while (*p != 0)
	{
//		mprintf("# (A) z=%d \"%s\".\n",z,p);
		while (strchr("0123456789+-.Ee",*p) == NULL)
			p++;
		q = p;
//		mprintf("# (B) z=%d \"%s\".\n",z,p);
		while ((strchr("0123456789+-.Ee",*q) != NULL) && (*q != 0))
			q++;
		if (*(q-1) == 0)
			q--;
//		mprintf("# (C) z=%d q=\"%s\" p=\"%s\".\n",z,q,p);
		if (q > p)
		{
//			mprintf("# q>p, q-p = %d\n",q-p);
//			mprintf("# p=\"%s\".\n",p);
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
//			mprintf("# buf = \"%s\".\n",buf);
			if (IsValidFloat(buf))
			{
				data[z] = atof(buf);
				z++;
			}
		}// else mprintf("# q <= p\n");
//		mprintf("# z=%d.\n",z);
		if (*q == 0)
			break;
//		mprintf("# X.\n");
		if (z == 3)
			break;
//		mprintf("# Y.\n");
		p = q;
	}
	if (z < 3)
	{
		eprintf("\nError: ExtractXYZCellGeometry3(): Comment line of XYZ trajectory does not contain 3 numbers.\n");
		eprintf("  (Comment line is \"%s\").   ",s);
		return;
	}

	if (g_bDoubleBox) {
		g_fBoxX *= g_iDoubleBoxX;
		g_fBoxY *= g_iDoubleBoxY;
		g_fBoxZ *= g_iDoubleBoxZ;
	}

	g_fBoxX = data[0] * 100.0; // Angstrom --> pm
	g_fBoxY = data[1] * 100.0;
	g_fBoxZ = data[2] * 100.0;

	g_mBoxFromOrtho(0,0) = g_fBoxX;
	g_mBoxFromOrtho(0,1) = 0;
	g_mBoxFromOrtho(0,2) = 0;

	g_mBoxFromOrtho(1,0) = 0;
	g_mBoxFromOrtho(1,1) = g_fBoxY;
	g_mBoxFromOrtho(1,2) = 0;

	g_mBoxFromOrtho(2,0) = 0;
	g_mBoxFromOrtho(2,1) = 0;
	g_mBoxFromOrtho(2,2) = g_fBoxZ;

	g_fBoxAngleA = 90.0;
	g_fBoxAngleB = 90.0;
	g_fBoxAngleC = 90.0;

	g_fBoxMinDiam = MIN3(g_fBoxX,g_fBoxY,g_fBoxZ);

	g_mBoxToOrtho = CxDMatrix3(g_mBoxFromOrtho);
	if (!g_mBoxToOrtho.Invert()) {
		eprintf("ExtractXYZCellGeometry3(): Error: Encountered singular cell matrix (cell volume is zero).\n");
		abort();
	}

/*	if (g_bDoubleBox)
	{
		g_mBoxFromOrtho(0,0) *= g_iDoubleBoxX;
		g_mBoxFromOrtho(1,0) *= g_iDoubleBoxX;
		g_mBoxFromOrtho(2,0) *= g_iDoubleBoxX;

		g_mBoxFromOrtho(0,1) *= g_iDoubleBoxY;
		g_mBoxFromOrtho(1,1) *= g_iDoubleBoxY;
		g_mBoxFromOrtho(2,1) *= g_iDoubleBoxY;

		g_mBoxFromOrtho(0,2) *= g_iDoubleBoxZ;
		g_mBoxFromOrtho(1,2) *= g_iDoubleBoxZ;
		g_mBoxFromOrtho(2,2) *= g_iDoubleBoxZ;
	}*/
}


bool IsElementMetal(const char *s) {

	// 2nd Period
	if (mystricmp(s,"Li") == 0) return true;

	// 3rd Period
	if (mystricmp(s,"Na") == 0) return true;
	if (mystricmp(s,"Mg") == 0) return true;
	if (mystricmp(s,"Al") == 0) return true;

	// 4th Period
	if (mystricmp(s, "K") == 0) return true;
	if (mystricmp(s,"Ca") == 0) return true;
	if (mystricmp(s,"Sc") == 0) return true;
	if (mystricmp(s,"Ti") == 0) return true;
	if (mystricmp(s, "V") == 0) return true;
	if (mystricmp(s,"Cr") == 0) return true;
	if (mystricmp(s,"Mn") == 0) return true;
	if (mystricmp(s,"Fe") == 0) return true;
	if (mystricmp(s,"Co") == 0) return true;
	if (mystricmp(s,"Ni") == 0) return true;
	if (mystricmp(s,"Cu") == 0) return true;
	if (mystricmp(s,"Zn") == 0) return true;
	if (mystricmp(s,"Ga") == 0) return true;
	if (mystricmp(s,"Ge") == 0) return true;

	// 5th Period
	if (mystricmp(s,"Rb") == 0) return true;
	if (mystricmp(s,"Sr") == 0) return true;
	if (mystricmp(s, "Y") == 0) return true;
	if (mystricmp(s,"Zr") == 0) return true;
	if (mystricmp(s,"Nb") == 0) return true;
	if (mystricmp(s,"Mo") == 0) return true;
	if (mystricmp(s,"Tc") == 0) return true;
	if (mystricmp(s,"Ru") == 0) return true;
	if (mystricmp(s,"Rh") == 0) return true;
	if (mystricmp(s,"Pd") == 0) return true;
	if (mystricmp(s,"Ag") == 0) return true;
	if (mystricmp(s,"Cd") == 0) return true;
	if (mystricmp(s,"In") == 0) return true;
	if (mystricmp(s,"Sn") == 0) return true;
	if (mystricmp(s,"Sb") == 0) return true;

	// 6th Period
	if (mystricmp(s,"Cs") == 0) return true;
	if (mystricmp(s,"Ba") == 0) return true;
	if (mystricmp(s,"La") == 0) return true;
	if (mystricmp(s,"Ce") == 0) return true;
	if (mystricmp(s,"Pr") == 0) return true;
	if (mystricmp(s,"Nd") == 0) return true;
	if (mystricmp(s,"Pm") == 0) return true;
	if (mystricmp(s,"Sm") == 0) return true;
	if (mystricmp(s,"Eu") == 0) return true;
	if (mystricmp(s,"Gd") == 0) return true;
	if (mystricmp(s,"Tb") == 0) return true;
	if (mystricmp(s,"Dy") == 0) return true;
	if (mystricmp(s,"Ho") == 0) return true;
	if (mystricmp(s,"Er") == 0) return true;
	if (mystricmp(s,"Tm") == 0) return true;
	if (mystricmp(s,"Yb") == 0) return true;
	if (mystricmp(s,"Lu") == 0) return true;
	if (mystricmp(s,"Hf") == 0) return true;
	if (mystricmp(s,"Ta") == 0) return true;
	if (mystricmp(s, "W") == 0) return true;
	if (mystricmp(s,"Re") == 0) return true;
	if (mystricmp(s,"Os") == 0) return true;
	if (mystricmp(s,"Ir") == 0) return true;
	if (mystricmp(s,"Pt") == 0) return true;
	if (mystricmp(s,"Au") == 0) return true;
	if (mystricmp(s,"Hg") == 0) return true;
	if (mystricmp(s,"Tl") == 0) return true;
	if (mystricmp(s,"Pb") == 0) return true;
	if (mystricmp(s,"Bi") == 0) return true;

	// 7th Period
	if (mystricmp(s,"Fr") == 0) return true;
	if (mystricmp(s,"Ra") == 0) return true;
	if (mystricmp(s,"Ac") == 0) return true;
	if (mystricmp(s,"Th") == 0) return true;
	if (mystricmp(s,"Pa") == 0) return true;
	if (mystricmp(s, "U") == 0) return true;
	if (mystricmp(s,"Np") == 0) return true;
	if (mystricmp(s,"Pu") == 0) return true;
	if (mystricmp(s,"Am") == 0) return true;
	if (mystricmp(s,"Cm") == 0) return true;
	if (mystricmp(s,"Bk") == 0) return true;
	if (mystricmp(s,"Cf") == 0) return true;
	if (mystricmp(s,"Es") == 0) return true;
	if (mystricmp(s,"Fm") == 0) return true;
	if (mystricmp(s,"Md") == 0) return true;
	if (mystricmp(s,"No") == 0) return true;
	if (mystricmp(s,"Lr") == 0) return true;

	return false;
}


bool IsElementNobleGas(const char *s) {

	if (mystricmp(s,"He") == 0) return true;
	if (mystricmp(s,"Ne") == 0) return true;
	if (mystricmp(s,"Ar") == 0) return true;
	if (mystricmp(s,"Kr") == 0) return true;
	if (mystricmp(s,"Xe") == 0) return true;
	if (mystricmp(s,"Rn") == 0) return true;

	return false;
}


void ParseCorrectWavenumber() {

	if (g_fIntegratorTimestep > 0)
		return;

	g_fIntegratorTimestep = AskFloat("    Please enter the timestep of the Verlet integrator used in the simulation: [%.2f] ",g_fTimestepLength,g_fTimestepLength);
}


double CorrectWavenumber(double w) {

	return sqrt( 1.0 - cos( 2.0 * 100.0 * g_fIntegratorTimestep * 1.0e-15 * Pi * CONST_SPEED_OF_LIGHT * w ) ) / ( sqrt(2.0) * 100.0 * Pi * CONST_SPEED_OF_LIGHT * g_fIntegratorTimestep * 1.0e-15 );
}





void RamanFromPolarizability(int argc, const char *argv[]) {

	FILE *a;
	char buf[1024], *p, *q;
	int z, z2, i, firstcol, skiphead, mode, fdir;
	int tia[4], tib[4], tic[4], tca[4], tcb[4], tcc[4], tialast[4], tiblast[4], ticlast[4];
	double tf, conv, field, tfs;
	bool ifberry;
	std::vector<CxDVector3> rawdip[4];
	CxDVector3 tvec;
	CROAEngine *roa;
	CROATrajectory *roatraj;
	CROAMolecule *roamol;
	CxDMatrix3 mat=CxDMatrix3(0), berry, invberry;
	CxString tabuf[4];


	mprintf("\n");
	mprintf(WHITE,"    ################################################\n");
	mprintf(WHITE,"    ####    Raman from Polarizability Module    ####\n");
	mprintf(WHITE,"    ####     (c) Martin Brehm, 2020 - 2022      ####\n");
	mprintf(WHITE,"    ####       https://brehm-research.de        ####\n");
	mprintf(WHITE,"    ################################################\n\n");
	mprintf("\n");


	g_bAdvanced2 = true;

	roa = new CROAEngine();

	roa->m_iMode = ROA_MODE_ANALYZE;
	roa->m_bMagMom = false;
	roa->m_bPola = true;

	roa->m_bDumpMolecularProps = false;
	roa->m_bDumpMol1Props = true;
	roa->m_bWriteACFs = true;
	roa->m_bSmoothData = false;
	roa->m_bReplaceOutliers = false;
	roa->m_bReplaceOutliersElectric = false;
	roa->m_bReplaceOutliersMagnetic = false;
	g_bQuadrupoleKeepTrace = false;
	roa->m_bReverseTraj = false;


	roatraj = new CROATrajectory();

	roa->m_oaTrajectories.push_back(roatraj);

	roamol = new CROAMolecule();

	roatraj->m_oaMolecules.push_back(roamol);


	mprintf("    The following operating modes are supported:\n\n");
	mprintf("    1.) Read total electric polarizability tensor (9 entries) from text file (one line per frame)\n");
	mprintf("    2.) Read total electric dipole vector with and without external electric field from 4 text files.\n");
	mprintf("\n");

	mode = AskRangeInteger( "    Which operating mode to select (1-2)? [1] ",1,2,1);

	mprintf("\n");

	if (mode == 1) { // Read polarizability tensor

		if (argc < 2) {
			eprintf("Error: Name of polarizability text file is expected as command line argument.\n\n");
			abort();
		}

		mprintf("    The entries of the polarizability tensor are expected as consecutive columns in a text file\n");
		mprintf("    (order XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ), where each line corresponds to one simulation step.\n");
		mprintf("    The columns can be separated by space, tabulators, and semicolons (or a mixture of those).\n\n");

		firstcol = AskUnsignedInteger("    Which is the first column for the tensor elements (XX)? [3] ",3) - 1;
		skiphead = AskUnsignedInteger("    How many lines to skip at the beginning of the file (header)? [1] ",1);

		mprintf("\n");

		conv = AskFloat("    Enter the factor to convert the polarizabilities from the text file to units of e*pm^2/V: [1.0] ",1.0);

		mprintf("\n");
		mprintf("    Trying to read data from \"%s\" ...\n",argv[1]);

		a = fopen( argv[1], "rt" );
		if (a == NULL) {
			eprintf("Error: Could not open file for reading.\n\n");
			abort();
		}

		i = 0;

		for (z=0;z<skiphead;z++)
			(void)!fgets(buf,1023,a);

		while (!feof(a)) {

			(void)!fgets(buf,1023,a);
			if (feof(a))
				break;

			p = buf;

			for (z=0;z<firstcol;z++) {
				while ((*p == ' ') || (*p == '\t') || (*p == ';'))
					p++;
				q = p;
				while ((*q != ' ') && (*q != '\t') && (*q != ';') && (*q != '\r') && (*q != '\n') && (*q != 0))
					q++;
				if (*q == 0) {
					eprintf("Error: Unexpected line end after column %d in line %d.\n\n",z+1,i+1);
					abort();
				}
				p = q+1;
			}

			for (z=0;z<9;z++) {
				while ((*p == ' ') || (*p == '\t') || (*p == ';'))
					p++;
				q = p;
				while ((*q != ' ') && (*q != '\t') && (*q != ';') && (*q != '\r') && (*q != '\n') && (*q != 0))
					q++;
				if ((*q == 0) && (z != 8)) {
					eprintf("Error: Unexpected line end after column %d in line %d.\n\n",z+1+firstcol,i+1);
					abort();
				}
				*q = 0;

				switch(z) {
					case 0:
						tf = atof(p);
						mat(0,0) = tf * conv;
						break;

					case 1:
						tf = atof(p);
						mat(1,0) = tf * conv;
						break;

					case 2:
						tf = atof(p);
						mat(2,0) = tf * conv;
						break;

					case 3:
						tf = atof(p);
						mat(0,1) = tf * conv;
						break;

					case 4:
						tf = atof(p);
						mat(1,1) = tf * conv;
						break;

					case 5:
						tf = atof(p);
						mat(2,1) = tf * conv;
						break;

					case 6:
						tf = atof(p);
						mat(0,2) = tf * conv;
						break;

					case 7:
						tf = atof(p);
						mat(1,2) = tf * conv;
						break;

					case 8:
						tf = atof(p);
						mat(2,2) = tf * conv;
						break;
				}

				p = q+1;
			}

			i++;

			roamol->m_maPolElDip.push_back(mat);
		}

		mprintf("    Read %d lines from file.\n",i);

		fclose(a);

		g_iMaxStep = i;
		g_iTrajSteps = i;
		g_iStride = 1;

	} else if (mode == 2) { // Read dipole vectors with and without field

		mprintf("    You need to supply four text files (one without external electric field, one for\n");
		mprintf("    each Cartesian field direction X, Y, Z) with dipole vectors of the total system.\n\n");

		mprintf("    The entries of the dipole vector are expected as consecutive columns in a text file\n");
		mprintf("    (order X, Y, Z), where each line corresponds to one simulation step.\n");
		mprintf("    The columns can be separated by space, tabulators, and semicolons (or a mixture of those).\n\n");

		firstcol = AskUnsignedInteger("    Which is the first column for the vector elements (X)? [1] ",1) - 1;
		skiphead = AskUnsignedInteger("    How many lines to skip at the beginning of the file (header)? [0] ",0);

		mprintf("\n");

		conv = AskFloat("    Enter the factor to convert the dipole vectors from the text files to units of Debye: [1.0] ",1.0);

		mprintf("\n");

		field = AskFloat("    Enter the field strength of the external field in units of a.u.: [0.005] ",0.005);

		field *= 0.514; // a.u. to V/pm

		mprintf("      This corresponds to %.7f V/pm.\n",field);

		g_iMaxStep = -1;

		for (fdir=0;fdir<4;fdir++) {

			mprintf("\n");

			switch(fdir) {
				case 0:
					AskString_ND( "    Enter the name of the text file with the field-free dipole vectors: ", &tabuf[fdir] );
					break;
				case 1:
					AskString_ND( "    Enter the name of the text file with the dipole vectors at field direction X: ", &tabuf[fdir] );
					break;
				case 2:
					AskString_ND( "    Enter the name of the text file with the dipole vectors at field direction Y: ", &tabuf[fdir] );
					break;
				case 3:
					AskString_ND( "    Enter the name of the text file with the dipole vectors at field direction Z: ", &tabuf[fdir] );
					break;
			}

			mprintf("\n");
			mprintf("      Trying to read data from \"%s\" ...\n",(const char*)tabuf[fdir]);

			a = fopen( (const char*)tabuf[fdir], "rt" );
			if (a == NULL) {
				eprintf("Error: Could not open file for reading.\n\n");
				abort();
			}

			i = 0;

			for (z=0;z<skiphead;z++)
				(void)!fgets(buf,1023,a);

			while (!feof(a)) {

				(void)!fgets(buf,1023,a);
				if (feof(a))
					break;

				p = buf;

				for (z=0;z<firstcol;z++) {
					while ((*p == ' ') || (*p == '\t') || (*p == ';'))
						p++;
					q = p;
					while ((*q != ' ') && (*q != '\t') && (*q != ';') && (*q != '\r') && (*q != '\n') && (*q != 0))
						q++;
					if (*q == 0) {
						eprintf("Error: Unexpected line end after column %d in line %d.\n\n",z+1,i+1);
						abort();
					}
					p = q+1;
				}

				for (z=0;z<3;z++) {
					while ((*p == ' ') || (*p == '\t') || (*p == ';'))
						p++;
					q = p;
					while ((*q != ' ') && (*q != '\t') && (*q != ';') && (*q != '\r') && (*q != '\n') && (*q != 0))
						q++;
					if ((*q == 0) && (z != 8)) {
						eprintf("Error: Unexpected line end after column %d in line %d.\n\n",z+1+firstcol,i+1);
						abort();
					}
					*q = 0;

					switch(z) {
						case 0:
							rawdip[fdir].push_back( CxDVector3(0,0,0) );
							rawdip[fdir].back()[0] = atof(p) * conv;
							break;

						case 1:
							rawdip[fdir].back()[1] = atof(p) * conv;
							break;

						case 2:
							rawdip[fdir].back()[2] = atof(p) * conv;
							break;
					}

					p = q+1;
				}

				i++;
			}

			mprintf("      Read %d lines from file.\n",i);

			if ((g_iMaxStep != -1) && (i != g_iMaxStep)) {
				eprintf("Error: Line count of %d differs from line count of field-free trajectory (%d).\n",i,g_iMaxStep);
				abort();
			}

			fclose(a);

			g_iMaxStep = i;
			g_iTrajSteps = i;
			g_iStride = 1;
		}

		mprintf("\n");

		mprintf("    If the dipole vectors are from Berry phase formalism, they are only defined up to a modulus.\n\n");

		ifberry = AskYesNo("    Are the dipole vectors from Berry phase formalism (y/n)? [yes] ",true);

		if (ifberry) {

			mprintf("\n");
			mprintf("    You now need to input the modulus. Look in the log file of your CP2k calculation\n");
			mprintf("    for the following section, and enter the nine numbers A_X to C_Z as marked below.\n\n");

			mprintf("    \"  Dipole vectors are based on the periodic (Berry phase) operator.        \"\n");
			mprintf("    \"  They are defined modulo integer multiples of the cell matrix [Debye].   \"\n");
			mprintf("    \"  [X] [   A_X   B_X   C_X   ] [i]                                         \"\n");
			mprintf("    \"  [Y]=[   A_Y   B_Y   C_Y   ]*[j]                                         \"\n");
			mprintf("    \"  [Z] [   A_Z   B_Z   C_Z   ] [k]                                         \"\n");
		
			mprintf("\n");
			berry(0,0) = AskFloat("      Enter component A_X: [0.0] ",0) * conv;
			berry(0,1) = AskFloat("      Enter component A_Y: [0.0] ",0) * conv;
			berry(0,2) = AskFloat("      Enter component A_Z: [0.0] ",0) * conv;
			mprintf("\n");
			berry(1,0) = AskFloat("      Enter component B_X: [0.0] ",0) * conv;
			berry(1,1) = AskFloat("      Enter component B_Y: [0.0] ",0) * conv;
			berry(1,2) = AskFloat("      Enter component B_Z: [0.0] ",0) * conv;
			mprintf("\n");
			berry(2,0) = AskFloat("      Enter component C_X: [0.0] ",0) * conv;
			berry(2,1) = AskFloat("      Enter component C_Y: [0.0] ",0) * conv;
			berry(2,2) = AskFloat("      Enter component C_Z: [0.0] ",0) * conv;
			mprintf("\n");

			mprintf("    Inverting matrix...\n");

			invberry = berry;
			if (!invberry.Invert()) {
				eprintf("Error: Could not invert Berry phase matrix!\n");
				abort();
			}

			mprintf("\n");
			mprintf("    Inverse Berry phase matrix (Debye^-1):\n");
			mprintf("      ( %16.10f %16.10f %16.10f )\n",invberry(0,0),invberry(1,0),invberry(2,0));
			mprintf("      ( %16.10f %16.10f %16.10f )\n",invberry(0,1),invberry(1,1),invberry(2,1));
			mprintf("      ( %16.10f %16.10f %16.10f )\n",invberry(0,2),invberry(1,2),invberry(2,2));
			mprintf("\n");

			mprintf("    Removing jumps from Berry phase formalism...\n");

			mprintf(WHITE,"      [");
			tfs = (double)g_iMaxStep / 60.0;

			// Bring dipole vector as close to the origin as possible in the first step
			for (z2=0;z2<4;z2++) {

				if (z2 == 0)
					tvec = invberry * rawdip[z2][0];
				else
					tvec = invberry * (rawdip[z2][0] - rawdip[0][0]);

				//mprintf("@ z2=%d\n",z2);
				//mprintf("@ rawdip = %16.10f %16.10f %16.10f\n",rawdip[z2][0][0],rawdip[z2][0][1],rawdip[z2][0][2]);
				//mprintf("@ tvec = %16.10f %16.10f %16.10f\n",tvec[0],tvec[1],tvec[2]);
				//mprintf("@ changeA = %16.10f %16.10f %16.10f\n",floor(tvec[0]+0.5) * berry(0,0),floor(tvec[0]+0.5) * berry(0,1),floor(tvec[0]+0.5) * berry(0,2));
				//mprintf("@ changeB = %16.10f %16.10f %16.10f\n",floor(tvec[0]+0.5) * berry(1,0),floor(tvec[0]+0.5) * berry(1,1),floor(tvec[0]+0.5) * berry(1,2));
				//mprintf("@ changeC = %16.10f %16.10f %16.10f\n",floor(tvec[0]+0.5) * berry(2,0),floor(tvec[0]+0.5) * berry(2,1),floor(tvec[0]+0.5) * berry(2,2));

				tia[z2] = (int)floor(tvec[0]+0.5);
				tib[z2] = (int)floor(tvec[1]+0.5);
				tic[z2] = (int)floor(tvec[2]+0.5);

				if (tia[z2] != 0) {
					rawdip[z2][0][0] -= tia[z2] * berry(0,0);
					rawdip[z2][0][1] -= tia[z2] * berry(0,1);
					rawdip[z2][0][2] -= tia[z2] * berry(0,2);
				}
				if (tib[z2] != 0) {
					rawdip[z2][0][0] -= tib[z2] * berry(1,0);
					rawdip[z2][0][1] -= tib[z2] * berry(1,1);
					rawdip[z2][0][2] -= tib[z2] * berry(1,2);
				}
				if (tic[z2] != 0) {
					rawdip[z2][0][0] -= tic[z2] * berry(2,0);
					rawdip[z2][0][1] -= tic[z2] * berry(2,1);
					rawdip[z2][0][2] -= tic[z2] * berry(2,2);
				}
				tialast[z2] = tia[z2];
				tiblast[z2] = tib[z2];
				ticlast[z2] = tic[z2];
			}

			for (z2=0;z2<4;z2++) {
				tca[z2] = 0;
				tcb[z2] = 0;
				tcc[z2] = 0;
			}

			// Remove jumps
			for (z=1;z<g_iMaxStep;z++) {

				for (z2=0;z2<4;z2++) {

					tvec = invberry * (rawdip[z2][z] - rawdip[z2][z-1]);

					tia[z2] = (int)floor(tvec[0]+0.5);
					tib[z2] = (int)floor(tvec[1]+0.5);
					tic[z2] = (int)floor(tvec[2]+0.5);

					if (tia[z2] != 0) {
						rawdip[z2][z][0] -= tia[z2] * berry(0,0);
						rawdip[z2][z][1] -= tia[z2] * berry(0,1);
						rawdip[z2][z][2] -= tia[z2] * berry(0,2);
						if (tia[z2] != tialast[z2])
							tca[z2]++;
					}
					if (tib[z2] != 0) {
						rawdip[z2][z][0] -= tib[z2] * berry(1,0);
						rawdip[z2][z][1] -= tib[z2] * berry(1,1);
						rawdip[z2][z][2] -= tib[z2] * berry(1,2);
						if (tib[z2] != tiblast[z2])
							tcb[z2]++;
					}
					if (tic[z2] != 0) {
						rawdip[z2][z][0] -= tic[z2] * berry(2,0);
						rawdip[z2][z][1] -= tic[z2] * berry(2,1);
						rawdip[z2][z][2] -= tic[z2] * berry(2,2);
						if (tic[z2] != ticlast[z2])
							tcc[z2]++;
					}

					tialast[z2] = tia[z2];
					tiblast[z2] = tib[z2];
					ticlast[z2] = tic[z2];
				}

				if (fmod(z,tfs) < 1.0)
					mprintf(WHITE,"#");
			}
			mprintf(WHITE,"] Done.\n");

			mprintf("\n");

			mprintf("    Number of removed jumps:\n\n");
			mprintf("      Trajectory   |   Along Cell Vector A   |   Along Cell Vector B   |   Along Cell Vector C\n");
			mprintf("     --------------|-------------------------|-------------------------|-----------------------\n");
			mprintf("      No Field     |                %6d   |                %6d   |                %6d\n",tca[0],tcb[0],tcc[0]);
			mprintf("      Field X      |                %6d   |                %6d   |                %6d\n",tca[0],tcb[0],tcc[0]);
			mprintf("      Field Y      |                %6d   |                %6d   |                %6d\n",tca[0],tcb[0],tcc[0]);
			mprintf("      Field Z      |                %6d   |                %6d   |                %6d\n",tca[0],tcb[0],tcc[0]);
			mprintf("\n");
		}

		mprintf("    Computing polarizabilities from finite differences...\n");
		mprintf("    Writing all data to \"result_polarizability.csv\" ...\n");

		a = OpenFileWrite("result_polarizability.csv",true);

		fprintf(a,"#Step;");
		fprintf(a,"  NoField DipX/Debye;  NoField DipY/Debye;  NoField DipZ/Debye;");
		fprintf(a,"  FieldX DipX/Debye;  FieldX DipY/Debye;  FieldX DipZ/Debye;");
		fprintf(a,"  FieldY DipX/Debye;  FieldY DipY/Debye;  FieldY DipZ/Debye;");
		fprintf(a,"  FieldZ DipX/Debye;  FieldZ DipY/Debye;  FieldZ DipZ/Debye;");
		fprintf(a,"  Polarizability_XX / e*pm^2/V;  Polarizability_XY / e*pm^2/V;  Polarizability_XZ / e*pm^2/V;");
		fprintf(a,"  Polarizability_YX / e*pm^2/V;  Polarizability_YY / e*pm^2/V;  Polarizability_YZ / e*pm^2/V;");
		fprintf(a,"  Polarizability_ZX / e*pm^2/V;  Polarizability_ZY / e*pm^2/V;  Polarizability_ZZ / e*pm^2/V\n");

		mprintf(WHITE,"      [");
		tfs = (double)g_iMaxStep / 60.0;
		for (z=0;z<g_iMaxStep;z++) {

			// 20.81943 - Debye to e*pm
			mat(0,0) = (rawdip[0][z][0] - rawdip[1][z][0]) * 20.81943 / field;
			mat(0,1) = (rawdip[0][z][1] - rawdip[1][z][1]) * 20.81943 / field;
			mat(0,2) = (rawdip[0][z][2] - rawdip[1][z][2]) * 20.81943 / field;

			mat(1,0) = (rawdip[0][z][0] - rawdip[2][z][0]) * 20.81943 / field;
			mat(1,1) = (rawdip[0][z][1] - rawdip[2][z][1]) * 20.81943 / field;
			mat(1,2) = (rawdip[0][z][2] - rawdip[2][z][2]) * 20.81943 / field;

			mat(2,0) = (rawdip[0][z][0] - rawdip[3][z][0]) * 20.81943 / field;
			mat(2,1) = (rawdip[0][z][1] - rawdip[3][z][1]) * 20.81943 / field;
			mat(2,2) = (rawdip[0][z][2] - rawdip[3][z][2]) * 20.81943 / field;

			roamol->m_maPolElDip.push_back(mat);

			fprintf(a,"%d;",z+1);
			fprintf(a,"  %.10E;  %.10E;  %.10E;",rawdip[0][z][0],rawdip[0][z][1],rawdip[0][z][2]);
			fprintf(a,"  %.10E;  %.10E;  %.10E;",rawdip[1][z][0],rawdip[1][z][1],rawdip[1][z][2]);
			fprintf(a,"  %.10E;  %.10E;  %.10E;",rawdip[2][z][0],rawdip[2][z][1],rawdip[2][z][2]);
			fprintf(a,"  %.10E;  %.10E;  %.10E;",rawdip[3][z][0],rawdip[3][z][1],rawdip[3][z][2]);
			fprintf(a,"  %.10E;  %.10E;  %.10E;",mat(0,0),mat(1,0),mat(2,0));
			fprintf(a,"  %.10E;  %.10E;  %.10E;",mat(0,1),mat(1,1),mat(2,1));
			fprintf(a,"  %.10E;  %.10E;  %.10E\n",mat(0,2),mat(1,2),mat(2,2));

			if (fmod(z,tfs) < 1.0)
				mprintf(WHITE,"#");
		}
		mprintf(WHITE,"] Done.\n");

		fclose(a);
	}


	mprintf("\n");

	g_fTimestepLength = AskFloat("    Enter the physical time distance between the lines in the text file (in fs): [4.0] ",4.0);

	CMolecule *m;
	CSingleMolecule *sm;

	m = new CMolecule();
	g_oaMolecules.Add(m);
	m->m_sName = new char[strlen("import")+1];
	strcpy(m->m_sName,"import");
	m->m_laSingleMolIndex.Add(0);

	sm = new CSingleMolecule();
	g_oaSingleMolecules.Add(sm);

	mprintf("\n");

	roa->ParseObservations();

	roa->m_iStepsProcessed = g_iTrajSteps;

	for (z=0;z<g_iTrajSteps;z++) {
		roamol->m_faCharge.push_back( 0 );
		roamol->m_vaElDip.push_back( CxDVector3(0,0,0) );
		roamol->m_maElQuad.push_back( CxDMatrix3(0,0,0,0,0,0,0,0,0) );
		roamol->m_taPolElQuad.push_back( CDTensor3(3,3,3) );
	}

	roa->Finish();

	mprintf(WHITE,"\n### Raman from Polarizability Module Leaving ###\n\n");
}



void WriteRevisionInfo() {

	std::vector<const char*> sa;
	unsigned int z, len;

	mprintf("Revision information of the source code files:\n\n");

	len = 0;

_begin:
	sa.push_back( GetRevisionInfo_2df( len ) );
	sa.push_back( GetRevisionInfo_acf( len ) );
	sa.push_back( GetRevisionInfo_aggrtopo( len ) );
	sa.push_back( GetRevisionInfo_analysisgroup( len ) );
	sa.push_back( GetRevisionInfo_asciiart( len ) );
	sa.push_back( GetRevisionInfo_atomgroup( len ) );
	sa.push_back( GetRevisionInfo_backtrace( len ) );
	sa.push_back( GetRevisionInfo_base64( len ) );
	sa.push_back( GetRevisionInfo_bicgstab( len ) );
	sa.push_back( GetRevisionInfo_bintools( len ) );
	sa.push_back( GetRevisionInfo_bintree( len ) );

	sa.push_back( GetRevisionInfo_bqb_alphabet( len ) );
	sa.push_back( GetRevisionInfo_bqb_bitset( len ) );
	sa.push_back( GetRevisionInfo_bqb_crc( len ) );
	sa.push_back( GetRevisionInfo_bqb_cubeframe( len ) );
	sa.push_back( GetRevisionInfo_bqb_driver( len ) );
	sa.push_back( GetRevisionInfo_bqb_engine( len ) );
	sa.push_back( GetRevisionInfo_bqb_extrapolator( len ) );
	sa.push_back( GetRevisionInfo_bqb_format( len ) );
	sa.push_back( GetRevisionInfo_bqb_hilbert( len ) );
	sa.push_back( GetRevisionInfo_bqb_hufftree( len ) );
	sa.push_back( GetRevisionInfo_bqb_integerengine( len ) );
	sa.push_back( GetRevisionInfo_bqb_interface( len ) );
	sa.push_back( GetRevisionInfo_bqb_largeinteger( len ) );
	sa.push_back( GetRevisionInfo_bqb_linalg( len ) );
	sa.push_back( GetRevisionInfo_bqb_math( len ) );
	sa.push_back( GetRevisionInfo_bqb_parmset( len ) );
	sa.push_back( GetRevisionInfo_bqb_tools( len ) );
	sa.push_back( GetRevisionInfo_bqbtool( len ) );

	sa.push_back( GetRevisionInfo_chargevar2( len ) );
	sa.push_back( GetRevisionInfo_chiral( len ) );
	sa.push_back( GetRevisionInfo_cluster( len ) );
	sa.push_back( GetRevisionInfo_contactmatrix( len ) );
	sa.push_back( GetRevisionInfo_cubetool( len ) );
	sa.push_back( GetRevisionInfo_dacf( len ) );
	sa.push_back( GetRevisionInfo_database( len ) );
	sa.push_back( GetRevisionInfo_df( len ) );
	sa.push_back( GetRevisionInfo_diagonal( len ) );
	sa.push_back( GetRevisionInfo_domain( len ) );
	sa.push_back( GetRevisionInfo_dpol( len ) );
	sa.push_back( GetRevisionInfo_elementdata( len ) );
	sa.push_back( GetRevisionInfo_fdf( len ) );
	sa.push_back( GetRevisionInfo_fft( len ) );
	sa.push_back( GetRevisionInfo_fixplproj( len ) );
	sa.push_back( GetRevisionInfo_gather( len ) );
	sa.push_back( GetRevisionInfo_geodens( len ) );
	sa.push_back( GetRevisionInfo_globalvar( len ) );
	sa.push_back( GetRevisionInfo_grace( len ) );
	sa.push_back( GetRevisionInfo_hbond( len ) );
	sa.push_back( GetRevisionInfo_interface( len ) );
	sa.push_back( GetRevisionInfo_internalcoord( len ) );
	sa.push_back( GetRevisionInfo_ionpair( len ) );
	sa.push_back( GetRevisionInfo_ir( len ) );
	sa.push_back( GetRevisionInfo_kiss_fft( len ) );
	sa.push_back( GetRevisionInfo_largeinteger( len ) );
	sa.push_back( GetRevisionInfo_linalg( len ) );
	sa.push_back( GetRevisionInfo_lmmin( len ) );
	sa.push_back( GetRevisionInfo_lmwrapper( len ) );
	sa.push_back( GetRevisionInfo_lu_decomp( len ) );
	sa.push_back( GetRevisionInfo_luzar( len ) );
	sa.push_back( GetRevisionInfo_maintools( len ) );
	sa.push_back( GetRevisionInfo_matrixplot( len ) );
	sa.push_back( GetRevisionInfo_moltools( len ) );
	sa.push_back( GetRevisionInfo_nbexchange( len ) );
	sa.push_back( GetRevisionInfo_nbsearch( len ) );
	sa.push_back( GetRevisionInfo_normalcoordinate( len ) );
	sa.push_back( GetRevisionInfo_normalmode( len ) );
	sa.push_back( GetRevisionInfo_order( len ) );
	sa.push_back( GetRevisionInfo_order_chain( len ) );
	sa.push_back( GetRevisionInfo_order_vector( len ) );
	sa.push_back( GetRevisionInfo_pdf( len ) );
	sa.push_back( GetRevisionInfo_plproj( len ) );
	sa.push_back( GetRevisionInfo_posdomain( len ) );
	sa.push_back( GetRevisionInfo_qr_fact( len ) );
	sa.push_back( GetRevisionInfo_raman( len ) );
	sa.push_back( GetRevisionInfo_random( len ) );
	sa.push_back( GetRevisionInfo_reflux( len ) );
	sa.push_back( GetRevisionInfo_region( len ) );
	sa.push_back( GetRevisionInfo_reordyn( len ) );
	sa.push_back( GetRevisionInfo_revsdf( len ) );
	sa.push_back( GetRevisionInfo_roa( len ) );
	sa.push_back( GetRevisionInfo_sankey( len ) );
	sa.push_back( GetRevisionInfo_sdfmap( len ) );
	sa.push_back( GetRevisionInfo_spectrum( len ) );
	sa.push_back( GetRevisionInfo_statistics( len ) );
	sa.push_back( GetRevisionInfo_structurefactor( len ) );
	sa.push_back( GetRevisionInfo_svgwriter( len ) );
	sa.push_back( GetRevisionInfo_tddf( len ) );
	sa.push_back( GetRevisionInfo_tensor( len ) );
	sa.push_back( GetRevisionInfo_tetrapak( len ) );
	sa.push_back( GetRevisionInfo_timestep( len ) );
	sa.push_back( GetRevisionInfo_tools( len ) );
	sa.push_back( GetRevisionInfo_travis( len ) );
	sa.push_back( GetRevisionInfo_void( len ) );
	sa.push_back( GetRevisionInfo_v_base( len ) );
	sa.push_back( GetRevisionInfo_v_base_wl( len ) );
	sa.push_back( GetRevisionInfo_v_cell( len ) );
	sa.push_back( GetRevisionInfo_v_common( len ) );
	sa.push_back( GetRevisionInfo_v_compute( len ) );
	sa.push_back( GetRevisionInfo_v_container( len ) );
	sa.push_back( GetRevisionInfo_v_container_prd( len ) );
	sa.push_back( GetRevisionInfo_v_pre_container( len ) );
	sa.push_back( GetRevisionInfo_v_unitcell( len ) );
	sa.push_back( GetRevisionInfo_v_wall( len ) );
	sa.push_back( GetRevisionInfo_vorowrapper( len ) );
	sa.push_back( GetRevisionInfo_xbytearray( len ) );
	sa.push_back( GetRevisionInfo_xdf( len ) );
	sa.push_back( GetRevisionInfo_xdmatrix3( len ) );
	sa.push_back( GetRevisionInfo_xdmatrixmn( len ) );
	sa.push_back( GetRevisionInfo_xdoublearray( len ) );
	sa.push_back( GetRevisionInfo_xdvec3array( len ) );
	sa.push_back( GetRevisionInfo_xdvector3( len ) );
	sa.push_back( GetRevisionInfo_xdvectorn( len ) );
	sa.push_back( GetRevisionInfo_xintarray( len ) );
	sa.push_back( GetRevisionInfo_xlongarray( len ) );
	sa.push_back( GetRevisionInfo_xmemfile( len ) );
	sa.push_back( GetRevisionInfo_xobarray( len ) );
	sa.push_back( GetRevisionInfo_xptrarray( len ) );
	sa.push_back( GetRevisionInfo_string( len ) );
	sa.push_back( GetRevisionInfo_xwordarray( len ) );
	sa.push_back( GetRevisionInfo_ziggurat( len ) );



	if (len == 0) {
		for (z=0;z<sa.size();z++)
			if (strlen(sa[z]) > len)
				len = (unsigned int)strlen(sa[z]);
		sa.clear();
		goto _begin;
	}

	for (z=0;z<sa.size();z++)
		mprintf( "  %s\n", sa[z] );

	mprintf( "\n");
}


void CheckSourceVersion() {

	const char *p;
	bool b;


	b = false;

	p = GetSourceVersion_2df();

	if (strcmp( p, GetSourceVersion_acf() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_aggrtopo() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_analysisgroup() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_asciiart() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_atomgroup() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_backtrace() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_base64() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bicgstab() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bintools() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bintree() ) != 0) b = true;

	if (strcmp( p, GetSourceVersion_bqb_alphabet() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_bitset() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_crc() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_cubeframe() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_driver() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_engine() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_extrapolator() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_format() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_hilbert() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_hufftree() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_integerengine() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_interface() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_largeinteger() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_linalg() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_math() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_parmset() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_tools() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqbtool() ) != 0) b = true;

	if (strcmp( p, GetSourceVersion_chargevar2() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_chiral() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_cluster() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_contactmatrix() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_cubetool() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_dacf() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_database() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_df() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_diagonal() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_domain() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_dpol() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_elementdata() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_fdf() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_fft() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_fixplproj() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_gather() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_geodens() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_globalvar() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_grace() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_hbond() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_interface() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_internalcoord() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_ionpair() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_ir() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_kiss_fft() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_largeinteger() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_linalg() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_lmmin() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_lmwrapper() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_lu_decomp() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_luzar() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_maintools() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_matrixplot() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_moltools() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_nbexchange() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_nbsearch() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_normalcoordinate() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_normalmode() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_order() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_order_chain() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_order_vector() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_pdf() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_plproj() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_posdomain() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_qr_fact() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_raman() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_random() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_reflux() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_region() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_reordyn() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_revsdf() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_roa() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_sankey() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_sdfmap() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_spectrum() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_statistics() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_structurefactor() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_svgwriter() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_tddf() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_tensor() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_tetrapak() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_timestep() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_tools() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_travis() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_void() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_v_base() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_v_base_wl() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_v_cell() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_v_common() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_v_compute() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_v_container() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_v_container_prd() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_v_pre_container() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_v_unitcell() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_v_wall() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_vorowrapper() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xbytearray() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xdf() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xdmatrix3() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xdmatrixmn() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xdoublearray() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xdvec3array() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xdvector3() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xdvectorn() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xintarray() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xlongarray() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xmemfile() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xobarray() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xptrarray() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_string() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_xwordarray() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_ziggurat() ) != 0) b = true;



	if (b) {

		mprintf("\n");
		eprintf("!!!!    Warning    !!!!\n");
		mprintf("The object files do not belong to the same source code version.\n");
		mprintf("Something went wrong with compiling. Expect problems and crashes.\n");
		mprintf("\n");

		WriteRevisionInfo();

		mprintf("\n");
		eprintf("!!!!    Warning    !!!!\n");
		mprintf("The object files do not belong to the same source code version.\n");
		mprintf("Something went wrong with compiling. Expect problems and crashes.\n");
		mprintf("\n");
	}
}


