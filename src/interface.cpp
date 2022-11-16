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

#include "interface.h"

#include "chiral.h"
#include "fdf.h"
#include "ir.h"
#include "normalcoordinate.h"
#include "normalmode.h"
#include "pdf.h"
#include "raman.h"
#include "region.h"
#include "structurefactor.h"
#include "hbond.h"
#include "ionpair.h"
#include "void.h"


const char *GetRevisionInfo_interface(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_interface() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


void Interface_DefaultConf() {
	/*	g_pDatabase->AddString("/BLA/BLUBB/PLOEPP/STRING1","String 1 Content");
	g_pDatabase->AddInt("/BLA/BLUBB/PLOEPP/INT1",123456);
	g_pDatabase->AddFloat("/BLA/BLUBB/PLOEPP/FLOAT1",1.23456);
	g_pDatabase->AddBool("/BLA/BLUBB/PLOEPP/BOOL1",true);*/
}

bool Interface_BeforeAnalysis() {
	if(g_bRegionAnalysis)
		if(!gatherRegionAnalysis())
			return false;
	if(g_bRaman)
		if(!gatherRaman())
			return false;
	if(g_bPDF)
		if(!gatherPDF())
			return false;
	if(g_bNormalCoordinate)
		if(!gatherNormalCoordinate())
			return false;
	if(g_bSFac)
		if(!gatherStructureFactor())
			return false;
	if(g_bChiral)
		if(!gatherChiral())
			return false;
	if(g_bSortWannier)
		if(!gatherSortWannier())
			return false;
	if(g_bEckartTransform)
		if(!gatherEckartTransform())
			return false;
	if(g_bPower || g_bVACFNew)
		if(!gatherPowerSpectrum())
			return false;
	if(g_bIR)
		if(!gatherIR())
			return false;
	if(g_bVCD)
		if(!gatherVCD())
			return false;

	
			
	if (g_bDipoleRestart)
		if (!gatherDipoleRestart())
			return false;
	if (g_bMagneticDipoleRestart)
		if (!gatherMagneticDipoleRestart())
			return false;
	
	if (g_bSetUpPolarizabilityCalc)
		if (!gatherPolarizabilityCalc())
			return false;
		
	if (g_bFDF)
		if (!gatherFDF())
			return false;
		
	if (g_bHBond)
		if (!gatherHBond())
			return false;
	
	if (g_bIonPair)
		if (!gatherIonPair())
			return false;


	if (g_bHole)
		if (!gatherVoid())
			return false;
	
	return true;
}

bool Interface_BeforeAnalysis2() {
	
	return true;
}

bool Interface_Initialization() {
	if(g_bRaman)
		if(!initializeRaman())
			return false;
	if(g_bPDF)
		if(!initializePDF())
			return false;
	if(g_bNormalCoordinate)
		if(!initializeNormalCoordinate())
			return false;
	if(g_bSFac)
		if(!initializeStructureFactor())
			return false;
	if(g_bChiral)
		if(!initializeChiral())
			return false;
	if(g_bSortWannier)
		if(!initializeSortWannier())
			return false;
	if(g_bEckartTransform)
		if(!initializeEckartTransform())
			return false;
	if(g_bPower || g_bVACFNew)
		if(!initializePowerSpectrum())
			return false;
	if(g_bIR)
		if(!initializeIR())
			return false;
	if(g_bVCD)
		if(!initializeVCD())
			return false;

	
	
	
	if (g_bDipoleRestart)
		if (!initializeDipoleRestart())
			return false;
	if (g_bMagneticDipoleRestart)
		if (!initializeMagneticDipoleRestart())
			return false;
	
	if (g_bSetUpPolarizabilityCalc)
		if (!initializePolarizabilityCalc())
			return false;

	if (g_bCubeTimeDev)
		g_fPDESolverInfoFile = OpenFileWrite("pdesolver_info", false);
		
	if (g_bFDF)
		if (!initializeFDF())
			return false;
		
	if (g_bHBond)
		if (!initializeHBond())
			return false;

	if (g_bIonPair)
		if (!initializeIonPair())
			return false;

	
	if (g_bHole)
		if (!initializeVoid())
			return false;

	return true;
}

void Interface_ProcessStep(CTimeStep *ts) {
	if(g_bRegionAnalysis)
		processRegionAnalysis(ts);
	if(g_bRaman)
		processRaman(ts);
	if(g_bPDF)
		processPDF(ts);
	if(g_bNormalCoordinate)
		processNormalCoordinate(ts);
	if(g_bSFac)
		processStructureFactor(ts);
	if(g_bChiral)
		processChiral(ts);
	if(g_bSortWannier)
		processSortWannier(ts);
	if(g_bEckartTransform)
		processEckartTransform(ts);
	if(g_bPower || g_bVACFNew)
		processPowerSpectrum(ts);
	if(g_bIR)
		processIR(ts);
	if(g_bVCD)
		processVCD(ts);

		
		
	
	if (g_bDipoleRestart)
		processDipoleRestart(ts);
	if (g_bMagneticDipoleRestart)
		processMagneticDipoleRestart(ts);
	
	if (g_bSetUpPolarizabilityCalc)
		processPolarizabilityCalc(ts);
	
	if (g_bFDF)
		processFDF(ts);
	
	if (g_bHBond)
		processHBond(ts);

	if (g_bIonPair)
		processIonPair(ts);


	if (g_bHole)
		processVoid(ts);
}

void Interface_AfterAnalysis() {
	if(g_bRegionAnalysis)
		finalizeRegionAnalysis();
	if(g_bRaman)
		finalizeRaman();
	if(g_bPDF)
		finalizePDF();
	if(g_bNormalCoordinate)
		finalizeNormalCoordinate();
	if(g_bSFac)
		finalizeStructureFactor();
	if(g_bChiral)
		finalizeChiral();
	if(g_bSortWannier)
		finalizeSortWannier();
	if(g_bEckartTransform)
		finalizeEckartTransform();
	if(g_bPower || g_bVACFNew)
		finalizePowerSpectrum();
	if(g_bIR)
		finalizeIR();
	if(g_bVCD)
		finalizeVCD();
		
		
		
	
	if (g_bDipoleRestart)
		finalizeDipoleRestart();
	if (g_bMagneticDipoleRestart)
		finalizeMagneticDipoleRestart();
	
	if (g_bSetUpPolarizabilityCalc)
		finalizePolarizabilityCalc();
	
	if (g_bCubeTimeDev)
		fclose(g_fPDESolverInfoFile);
		
	if (g_bFDF)
		finalizeFDF();
	
	if (g_bHBond)
		finalizeHBond();
	
	if (g_bIonPair)
		finalizeIonPair();
	
	
	if (g_bHole)
		finalizeVoid();
}

void Interface_DecomposeModes(int n, CxObArray *ccr_matrix) {
	normalModeAnalysis(n, ccr_matrix);
}

/************************************************ 
 
  Output: 
 
  mprintf( [Color,] ... ) 
 
#define GREY 8 
#define BLUE 9 
#define GREEN 10 
#define CYAN 11 
#define RED 12 
#define PINK 13 
#define YELLOW 14 
#define WHITE 15 
 
 
  A list of useful global variables: 
 
  g_bAdvanced2 
 
  g_fTimestepLength 
 
  g_iTrajSteps  // -1 means no information 
 
  g_fBoxX, g_fBoxY, g_fBoxZ  -  The cell vector in pm 
 
  g_oaMolecules              -  Array of all molecule kinds (members of type CMolecule) 
   
  g_oaSingleMolecules        -  Array of all molecules (members of type CSingleMolecule) 
 
 
 
		for (z=0;z<g_oaMolecules.GetSize();z++) 
		{ 
			m = (CMolecule*)g_oaMolecules[z]; 
			if (m->m_bPseudo) 
				continue; 
			for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++) 
			{ 
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]]; 
				for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++) 
				{ 
					if ((!g_bSaveVirtAtoms) && (m->m_baAtomIndex[z3] == g_iVirtAtomType)) 
						continue; 
					for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++) 
						mfprintf(a,"  %s  %8.5f  %8.5f  %8.5f\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][0]/100.0f,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][1]/100.0f,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][2]/100.0f); 
				} 
			} 
		} 
 
 
*************************************************/
