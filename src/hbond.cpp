/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Sascha Gehrke.

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

#include "hbond.h"
#include "luzar.h"
#include <math.h>

#include "globalvar.h"
#include "maintools.h"


const char *GetRevisionInfo_hbond(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_hbond() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


static CxObArray g_HBondObserv;


bool gatherHBond()
{
	CHBondObservation *o;

	try { o = new CHBondObservation(); } catch(...) { o = NULL; }
	if (o == NULL) NewException((double)sizeof(CHBondObservation),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	g_HBondObserv.Add(o);

	return true;
}


bool initializeHBond() 
{
	int i;

	for(i = 0; i < g_HBondObserv.GetSize(); i++) 
	{
		((CHBondObservation *)g_HBondObserv[i])->initialize();
	}
	return true;
}


bool processHBond(CTimeStep* ts) 
{
	int i;
	for(i = 0; i < g_HBondObserv.GetSize(); i++) 
	{
		((CHBondObservation *)g_HBondObserv[i])->process(ts);
	}
	return true;
}


bool finalizeHBond() 
{
	int i;
	for(i = 0; i < g_HBondObserv.GetSize(); i++) 
	{
		((CHBondObservation *)g_HBondObserv[i])->finalize();
	}
	return true;
}


CHBondObservation::CHBondObservation()
{
	int z, z2;
	CxString buf, buf2;

	try { m_laAccep = new CxIntArray(); } catch(...) { m_laAccep = NULL; }
	if (m_laAccep == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_laHydro = new CxIntArray(); } catch(...) { m_laHydro = NULL; }
	if (m_laHydro == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_laDonor = new CxIntArray(); } catch(...) { m_laDonor = NULL; }
	if (m_laDonor == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if (g_oaMolecules.GetSize() > 1)
	{
		buf.sprintf("\n    Which of the molecules should be the donor molecule (");
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
			buf.strcat(buf2);
			if (z < g_oaMolecules.GetSize()-1)
				buf.strcat(", ");
		}
		buf.strcat(")? ");
		m_iShowMol = AskRangeInteger_ND("%s",1,g_oaMolecules.GetSize(),(const char*)buf) - 1;
		mprintf(WHITE,"\n    %s is the donor molecule.\n\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);

		buf.sprintf("\n    Which of the molecules should be the acceptor molecule (");
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
			buf.strcat(buf2);
			if (z < g_oaMolecules.GetSize()-1)
				buf.strcat(", ");
		}
		buf.strcat(")? ");
		m_iShowMol2 = AskRangeInteger_ND("%s",1,g_oaMolecules.GetSize(),(const char*)buf) - 1;
		mprintf(WHITE,"\n    %s is the acceptor molecule.\n\n",((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName);

	} else 
	{
		m_iShowMol = 0; 	
		m_iShowMol2 = 0;		
	}

	mprintf(WHITE,"    Please note:");
	mprintf(" You have to enter the hydrogen bond donor atoms below.\n");
	mprintf("                 In this analysis, the donors are not the hydrogen atoms themselves,\n");
	mprintf("                 but the atoms they are bound to (e.g., O, N, C). Enter those as donors.\n");
	mprintf("                 The relevant hydrogen atoms will then be automatically detected.\n\n");

	mprintf("    Which atom(s) from %s are involved as hydrogen bond donor (e.g. \"C1,O3-5,N\")?  ",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	inpprintf("! Which atom(s) from %s are involved as hydrogen bond donor (e.g. \"C1,O3-5,N\")?  \n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	myget(&buf);
	while ((strlen(buf) == 0) || (!checkAtomChoice((CMolecule*)g_oaMolecules[m_iShowMol],buf,((CxIntArray*)m_laDonor))))
	{
		mprintf("    Please enter donating atom!  ");
		myget(&buf);
	}

	mprintf("\n    Added %d donor atoms in %d %s molecules:",m_laDonor->GetSize(),((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(),((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	for (z=0;z<m_laDonor->GetSize();z++) {
		if (g_waAtomMolIndex[m_laDonor->GetAt(z)] != m_iShowMol) {
			eprintf("\nError: Donor atom %s[%d] %s%d is not from donor molecule %s!\n",((CMolecule*)g_oaMolecules[g_waAtomMolIndex[m_laDonor->GetAt(z)]])->m_sName,g_laAtomSMLocalIndex[m_laDonor->GetAt(z)]+1,(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[m_laDonor->GetAt(z)]])->m_sName,g_waAtomMolNumber[m_laDonor->GetAt(z)]+1,((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
			abort();
		}
		if (g_laAtomSMLocalIndex[m_laDonor->GetAt(z)] == 0)
			mprintf(" %s%d",(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[m_laDonor->GetAt(z)]])->m_sName,g_waAtomMolNumber[m_laDonor->GetAt(z)]+1);
	}
	for (z=0;z<m_laDonor->GetSize();z++) {
		for (z2=z+1;z2<m_laDonor->GetSize();z2++) {
			if (m_laDonor->GetAt(z) == m_laDonor->GetAt(z2)) {
				eprintf("\nError: Donor atom %s[%d] %s%d is listed more than one time!\n",((CMolecule*)g_oaMolecules[g_waAtomMolIndex[m_laDonor->GetAt(z)]])->m_sName,g_laAtomSMLocalIndex[m_laDonor->GetAt(z)]+1,(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[m_laDonor->GetAt(z)]])->m_sName,g_waAtomMolNumber[m_laDonor->GetAt(z)]+1);
				abort();
			}
		}
	}
	mprintf("\n\n");

	mprintf("    Which atom(s) from %s are involved as hydrogen bond acceptor (e.g. \"C1,O3-5,N\")?  ",((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName);
	inpprintf("! Which atom(s) from %s are involved as hydrogen bond acceptor (e.g. \"C1,O3-5,N\")?  \n",((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName);
	myget(&buf);
	while ((strlen(buf) == 0) || (!checkAtomChoice((CMolecule*)g_oaMolecules[m_iShowMol2],buf,((CxIntArray*)m_laAccep))))
	{
		mprintf("    Please enter accepting atom!  ");
		myget(&buf);
	}

	mprintf("\n    Added %d acceptor atoms in %d %s molecules:",m_laAccep->GetSize(),((CMolecule*)g_oaMolecules[m_iShowMol2])->m_laSingleMolIndex.GetSize(),((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName);
	for (z=0;z<m_laAccep->GetSize();z++) {
		if (g_waAtomMolIndex[m_laAccep->GetAt(z)] != m_iShowMol2) {
			eprintf("\nError: Acceptor atom %s[%d] %s%d is not from acceptor molecule %s!\n",((CMolecule*)g_oaMolecules[g_waAtomMolIndex[m_laAccep->GetAt(z)]])->m_sName,g_laAtomSMLocalIndex[m_laAccep->GetAt(z)]+1,(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[m_laAccep->GetAt(z)]])->m_sName,g_waAtomMolNumber[m_laAccep->GetAt(z)]+1,((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName);
			abort();
		}
		if (g_laAtomSMLocalIndex[m_laAccep->GetAt(z)] == 0)
			mprintf(" %s%d",(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[m_laAccep->GetAt(z)]])->m_sName,g_waAtomMolNumber[m_laAccep->GetAt(z)]+1);
	}
	for (z=0;z<m_laAccep->GetSize();z++) {
		for (z2=z+1;z2<m_laAccep->GetSize();z2++) {
			if (m_laAccep->GetAt(z) == m_laAccep->GetAt(z2)) {
				eprintf("\nError: Acceptor atom %s[%d] %s%d is listed more than one time!\n",((CMolecule*)g_oaMolecules[g_waAtomMolIndex[m_laAccep->GetAt(z)]])->m_sName,g_laAtomSMLocalIndex[m_laAccep->GetAt(z)]+1,(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[m_laAccep->GetAt(z)]])->m_sName,g_waAtomMolNumber[m_laAccep->GetAt(z)]+1);
				abort();
			}
		}
	}
	mprintf("\n\n");

	BuildHydroArray();

	m_bLuzar = AskYesNo("    Analyse bond per donor-acceptor pair (y) or bond per hydrogen atom (n)? [no] ",false);

	if ( !m_bLuzar )	
	{
		BuildDonorArray();

		if ( m_laHydro->GetSize() != m_laDonor->GetSize() )
		{
			eprintf("Number of hydrogens is not equal number of donor atoms. Something is wrong. Aborting ...\n");
			abort();
		}
	}
	else
	{
		eprintf("Sorry! The \"donor-acceptor pair\" mode is currently not working properly. Aborting ...\n");
		abort();
	}

	m_frAD = 0;
	m_frAH = 0;
	m_fwinkel = 0;

	m_bHysteresis = AskYesNo("    Use hysteresis for the hydrogen bond criterion (y/n)? [no]", false);	

	while(m_frAD == 0 || m_frAH == 0)
	{
		m_frAH = AskFloat("    Maximum distance between hydrogen atom and accepting atom for bond formation in pm? [245] ",245.0);
		if ( m_bHysteresis )
		{
			m_frAH_out = AskFloat("    Maximum distance between hydrogen atom and accepting atom for bond breaking in pm? [%f] ",m_frAH,m_frAH);
		}		
		else
		{
			m_frAH_out = m_frAH;
		}
		
		m_frAD = AskFloat("    Maximum distance between donating atom and accepting atom for bond formation in pm? [%d] ",int(m_frAH + 105),int(m_frAH + 105));
		if ( m_bHysteresis )
		{
			m_frAD_out = AskFloat("    Maximum distance between donating atom and accepting atom for bond formation in pm? [%f] ",m_frAD,m_frAD);
		}
		else
		{
			m_frAD_out = m_frAD;
		}

		if(m_frAD < m_frAH)
		{
			if(!AskYesNo("    The distance between acceptor and donor is smaller than the distance between acceptor and hydrogen. This is very unusual. Keep the values anyway (y/n)? [yes] ",true))
			{
				m_frAD = 0;
				m_frAH = 0;
			}
		}
	}
	while( m_fwinkel == 0 )
	{
		m_fwinkel = cos((M_PI/180)*AskRangeFloat("    Maximum angle between donor-acceptor and donor-hydrogen for bond formation? [30] ", 0.0, 180.0, 30.0));
		if ( m_bHysteresis )
		{
			m_fwinkel_out = cos((M_PI/180)*AskRangeFloat("    Maximum angle between donor-acceptor and donor-hydrogen for bond breaking? [%f] ", 0.0, 180.0, acos(m_fwinkel) * 180 / M_PI, acos(m_fwinkel) * 180 / M_PI));
		}
		else
		{
			m_fwinkel_out = m_fwinkel;
		}
	}
	
	g_fTimestepLength = AskFloat("    Enter the physical time distance between successive trajectory frames in fs: [0.5] ",0.5f);
}


CHBondObservation::~CHBondObservation()
{
}


void CHBondObservation::initialize()
{
	int n,m,m_start = 0;

	try { m_laLastNeigh = new CxIntArray(); } catch(...) { m_laLastNeigh = NULL; }
	if (m_laLastNeigh == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_laLastBond = new CxIntArray(); } catch(...) { m_laLastBond = NULL; }
	if (m_laLastBond == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_oaNeigh = new CxObArray(); } catch(...) { m_oaNeigh = NULL; }
	if (m_oaNeigh == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_oaBond = new CxObArray(); } catch(...) { m_oaBond = NULL; }
	if (m_oaBond == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_oaPartners = new CxObArray(); } catch(...) { m_oaPartners = NULL; }
	if (m_oaPartners == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if ( g_iMaxStep != -1 )
		m_iDepth = g_iMaxStep * 3 / 4 / g_iStride;
	else
		m_iDepth = g_iTrajSteps * 3 / 4 / g_iStride;

	m_iDepth = AskUnsignedInteger("    Enter the correlation depth of the ACF (in trajectory frames): [%d] ", m_iDepth, m_iDepth);
	if ( m_iDepth > g_iTrajSteps )
	{
		mprintf("    Correlation Depth is exceeding number of total timesteps. Truncating ....\n");
		m_iDepth = g_iTrajSteps;
	}

	m_iTrunc = m_iDepth;

	while ( m_iTrunc > (int)(m_iDepth/10) )
	{
		m_iTrunc = (int)(50000/(g_fTimestepLength*g_iStride));
		if ( m_iTrunc > m_iDepth/10 )
			m_iTrunc = (int)(m_iDepth/10);	

		m_iTrunc = AskUnsignedInteger("    Enter the truncation for the fitting procedure (in trajectory frames): [%d] ", m_iTrunc, m_iTrunc);

		if ( m_iTrunc > m_iDepth/10 )
			mprintf("    Truncation must be at least 10 times smaller than correlation depth!\n");
	}

	if ( !m_bLuzar )			// dynamics of hydrogen bonds
	{
		for ( n=0; n < m_laAccep->GetSize(); n++)
		{
			for (m=0; m < m_laDonor->GetSize(); m++)
			{
				if ( m_laDonor->GetAt(m) != m_laAccep->GetAt(n) )	// Donor could be Accep ... processing would not make any sense ...
				{
					CxIntArray* IncludeME;
					CxIntArray* IncludeME2;

					try { IncludeME = new CxIntArray(); } catch(...) { IncludeME = NULL; }
					if (IncludeME == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					try { IncludeME2 = new CxIntArray(); } catch(...) { IncludeME2 = NULL; }
					if (IncludeME2 == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					m_laLastNeigh->Add(0);
					m_laLastBond->Add(0);
					m_oaNeigh->Add(IncludeME);
					m_oaBond->Add(IncludeME2);
				}
			}	
		}
	}
	else					// dynamics of hydrogen bonded molecules 
	{
		if ( m_laAccep->Contains(m_laDonor->GetAt(0)) ) 	// Donor and Acceptor same atomtype?
		{
			m_bSameAtom = true;
		}
		else
		{
			m_bSameAtom = false;
		}
		for ( n=0; n < m_laAccep->GetSize(); n++ )		
		{
			if ( m_bSameAtom )
			{
				m_start = n + 1;			// If same atomtype: No doubles!
			}
			for ( m = m_start; m < m_laDonor->GetSize(); m++ )
			{
				CxIntArray* IncludeME;
				CxIntArray* IncludeME2;
				CxIntArray* IncludeME3;

				try { IncludeME = new CxIntArray(); } catch(...) { IncludeME = NULL; }
				if (IncludeME == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				try { IncludeME2 = new CxIntArray(); } catch(...) { IncludeME2 = NULL; }
				if (IncludeME2 == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				try { IncludeME3 = new CxIntArray(); } catch(...) { IncludeME3 = NULL; }
				if (IncludeME3 == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				m_laLastNeigh->Add(0);
				m_laLastBond->Add(0);
				m_oaNeigh->Add(IncludeME);
				m_oaBond->Add(IncludeME2);

				FindHydro(IncludeME3,m_laAccep->GetAt(n));
				FindHydro(IncludeME3,m_laDonor->GetAt(m));
				m_oaPartners->Add(IncludeME3);
			}
		}
	}
}


void CHBondObservation::FindHydro(CxIntArray* IncludeME3, int donor)
{
	int m,i,j;

	for (m=0; m < ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(); m++)		// Walk over all SingleMolecules
	{
		for (i=0; i < ((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetAt(m)])->m_oaMolAtoms.GetSize(); i++)	// Walk over all Atoms in Molecule
		{
			if ( ((CMolAtom*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetAt(m)])->m_oaMolAtoms[i])->m_iOffset == donor ) // Search Atom
			{
				for (j=0; j < ((CMolAtom*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetAt(m)])->m_oaMolAtoms[i])->m_oaBonds.GetSize(); j++) // Walk over all bondpartners
				{
					if ( ((CAtom*)g_oaAtoms[((CMolAtom*)((CMolAtom*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetAt(m)])->m_oaMolAtoms[i])->m_oaBonds[j])->m_iType])->m_pElement->m_iOrd == 1 )
					{
						IncludeME3->Add(((CMolAtom*)((CMolAtom*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetAt(m)])->m_oaMolAtoms[i])->m_oaBonds[j])->m_iOffset);
					}
				}
				m = ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
				break;
			}
		}
	}
}


void CHBondObservation::process(CTimeStep* ts)
{	
	int n,m,k,i=0,m_start=0;

	if ( !m_bLuzar )		// Dynamics of hydrogen bonds: Check for all combinations if conditions are fulfilled
	{
		for (n=0; n < m_laAccep->GetSize(); n++)
		{
			for (m=0; m < m_laDonor->GetSize(); m++)
			{
				if ( m_laDonor->GetAt(m) != m_laAccep->GetAt(n) )	// Donor could be Accep ... processing would not make any sense ...
				{
					if ( checkCOND(ts->m_vaCoords[m_laAccep->GetAt(n)], ts->m_vaCoords[m_laHydro->GetAt(m)], ts->m_vaCoords[m_laDonor->GetAt(m)], i) )
					{
						mprintf(RED,"Something is wrong ... aborting ...");
						abort();
					}
					i++;
				}
			}
		}
	}
	else				// Dynamics of hydrogen bonded molecules: Check every pair only once
	{
		for ( n=0; n < m_laAccep->GetSize(); n++ )
		{
			if ( m_bSameAtom )
			{
				m_start = n + 1;						// If same atomtype: No doubles!
			}
			for ( m = m_start; m < m_laDonor->GetSize(); m++ )
			{
				for ( k=0; k < ((CxIntArray*)m_oaPartners->GetAt(i))->GetSize(); k++ )	// Check condition for every possible hydrogenatom
				{
					if ( checkCOND(ts->m_vaCoords[m_laAccep->GetAt(n)], ts->m_vaCoords[((CxIntArray*)m_oaPartners->GetAt(i))->GetAt(k)], ts->m_vaCoords[m_laDonor->GetAt(m)], i) )
					{
						if ( k == (((CxIntArray*)m_oaPartners->GetAt(i))->GetSize()-1) )
						{
							m_laLastBond->SetAt(i,0);						  // Kein H-Bond
							((CxIntArray*)m_oaBond->GetAt(i))->Add(g_iSteps);
						}
						continue;
					}
					if ( m_laLastNeigh->GetAt(i) == 0 || m_laLastBond->GetAt(i) == 1 ) // No Neighbor or hbond fulfilled
					{
						break;
					}
				} 
				i++;
			}
		}
	}
}


void CHBondObservation::finalize()	// reconstruct functions -> autocorrelate -> sum results -> normalize -> fit
{
	long n,i,interval;
	double BondCor, DiffCor;
	CxString buf,buf2;
	unsigned long t0,eta;

	CxDoubleArray* NeighFunction;	// Input
	CxDoubleArray* BondFunction;
	
	CxDoubleArray* temp1Function;	// Temp
	CxDoubleArray* temp2Function;

	CLuzarCorrelation* luzar;

	try { BondFunction = new CxDoubleArray(); } catch(...) { BondFunction = NULL; }
	if (BondFunction == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { NeighFunction = new CxDoubleArray(); } catch(...) { NeighFunction = NULL; }
	if (NeighFunction == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { BondDynamic = new CxDoubleArray(); } catch(...) { BondDynamic = NULL; }
	if (BondDynamic == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { DiffusionDynamic = new CxDoubleArray(); } catch(...) { DiffusionDynamic = NULL; }
	if (DiffusionDynamic == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { k_t = new CxDoubleArray(); } catch(...) { k_t = NULL; }
	if (k_t == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { a = new CxDoubleArray(); } catch(...) { a = NULL; }
	if (a == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { b = new CxDoubleArray(); } catch(...) { b = NULL; }
	if (b == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { temp1Function = new CxDoubleArray(); } catch(...) { temp1Function = NULL; }
	if (temp1Function == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { temp2Function = new CxDoubleArray(); } catch(...) { temp2Function = NULL; }
	if (temp2Function == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for ( i=0; i < m_iDepth; i++ )
	{
		BondDynamic->Add(0.0);
		DiffusionDynamic->Add(0.0);
	}

	try { luzar = new CLuzarCorrelation(); } catch(...) { luzar = NULL; }
	if (luzar == NULL) NewException((double)sizeof(CLuzarCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	NeighFunction->SetSize(g_iSteps / g_iStride);
	BondFunction->SetSize(g_iSteps / g_iStride);

	luzar->Init(g_iSteps / g_iStride,m_iDepth); 	// Inititalize Luzar-Correlation

	mprintf("    Finishing %6i Functions ... \nFunc %6li ",m_laLastNeigh->GetSize(),(long)0);
	g_iDotCounter = 0;
	t0 = (unsigned long)time(NULL);
	interval = 2;

	for ( n=0; n < m_laLastNeigh->GetSize(); n++)		// Walk over all pairs
	{
		if ( n % interval == 0 )
		{
			if ( g_iDotCounter >= 50 )
			{
				if ( interval == 2 && n > 100 )
				{
					interval = 20;
				}
				else if ( interval == 20 && n > 2000 )
				{
					interval = 200;
				}
				g_iDotCounter = 0;
				if ((time(NULL) - t0) > 5)
				{
					eta = (unsigned long)(((double)time(NULL) - t0) / n * ((double)MAX((long)0,(long)m_laLastNeigh->GetSize() - (long)n)));
					FormatTime(eta,&buf);
					mprintf(" ETA %s",(const char*)buf);
				}
				mprintf("\nFunc %6li ",n);
			}
			g_iDotCounter++;
			mprintf(".");
		}

		if ( ((CxIntArray*)m_oaNeigh->GetAt(n))->GetSize() > 0 )
		{
			RestoreFunction(((CxIntArray*)m_oaNeigh->GetAt(n)),((CxDoubleArray*)NeighFunction));	// Restore function of pair n as NeighFunction and BondFunction
			RestoreFunction(((CxIntArray*)m_oaBond->GetAt(n)),((CxDoubleArray*)BondFunction));	

			luzar->Correlate(BondFunction,NeighFunction,temp1Function,temp2Function);

			for ( i=0; i < m_iDepth; i++ )		// Sum results to Output-array
			{
				(*BondDynamic)[i] += (*temp1Function)[i];
				(*DiffusionDynamic)[i] += (*temp2Function)[i];
			}
		}
	}
	mprintf("\n\n");

	if ( (*BondDynamic)[0] == 0 || (*DiffusionDynamic)[0] == 0 )
	{
		eprintf("First correlation value is zero ... something is wrong ...\n");
		abort();
	}

	BondCor = (double)(*BondDynamic)[0];
	DiffCor = (double)(*DiffusionDynamic)[0];

	for ( i=0; i < m_iDepth; i++ )			// Normalize
	{
		(*BondDynamic)[i] /= BondCor;
		(*DiffusionDynamic)[i] /= DiffCor;
		(*DiffusionDynamic)[i] -= (*BondDynamic)[i];
	}

	luzar->Fit(BondDynamic, DiffusionDynamic, k_t, a, b, m_iTrunc);

	WriteOutput();
}


void CHBondObservation::RestoreFunction(CxIntArray* input, CxDoubleArray* output)
{
	unsigned int n;
	bool eins = false;

	for ( n=0; n < (g_iSteps / g_iStride); n++ )		// Walk over all processed steps
	{
		if ( input->Contains((n+1)*g_iStride) )
		{
			eins = !eins;
		}
		if (!eins)
		{
			output->SetAt(n,0.0);
		}
		else
		{
			output->SetAt(n,1.0);
		}
	}
}


bool CHBondObservation::includeAtom(int atom, int number, int mol, CxIntArray* target)	// atom: element number, number: atomsort number, mol m_iShowMol, target: in which array to include
{
	CxIntArray includeME;

	for ( int n = 0; n < ((CMolecule*)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize(); n++)
	{

		includeME.CopyFrom((CxIntArray*)(((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[mol])->m_laSingleMolIndex[n]])->m_oaAtomOffset[atom]));
		if (number >= includeME.GetSize()) {
			eprintf("CHBondObservation::includeAtom(): Internal error (%d/%d).\n",number,includeME.GetSize());
			abort();
		}

		target->Add(includeME[number]);
	}
	
	return true;
}


bool CHBondObservation::checkAtomChoice(CMolecule *mol, const char *s, CxIntArray* m_laInclude)
{
	const char *p, *q;
	char buf[32];
	const char *allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890,-#_ ";
	int inclu, la = -1, i = -1, i2 = -1, atom = -1; 
	bool m = false, fin = false;

	p = s;

	while (*p != 0)
	{
		fin = false;
		while (*p == ' ')
			p++;
		if (strchr(allowed,*p) == NULL)
		{
			eprintf("Error: Character \"%c\" not allowed.\n",*p);
			return false;
		}
		q = p;
		if (*q == '#' || *q == '*')
		{
			eprintf("No virtual atoms allowed!");
			return false;
		}
		if (isalpha(*q))
		{
			if (m)
			{
				eprintf("Error: Only digit allowed after \"-\".\n");
				return false;
			}
			while (isalpha(*q))
				q++;
			if (q-p >= 32)
			{
				eprintf("Internal Error (%ld >= 32).\n",(long)(q-p));
				return false;
			}
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
			atom = mol->FindAtomInMol(buf); 	// atom holds number of element 
			if (atom == -1)
			{
				eprintf("Error: Atom \"%s\" not in the molecule.\n",buf);
				return false;
			} 
		} 
		else if (isdigit(*q))
		{
			if (atom == -1)
				atom = la;
			if (atom == -1)
			{
				eprintf("Error: Number in the beginning not possible.\n");
				return false;
			}
			while (isdigit(*q)) 
				q++;
			if ((*q != '-') && (*q != ',') && (*q != ' ') && (*q != 0))
			{
				eprintf("Error: Only \",\" or \"-\" may follow after a number.\n");
				return false;
			}
			if (q-p >= 32)
			{
				eprintf("Internal Error (%ld >= 32).\n",(long)(q-p));
				return false;
			}
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
			if (atoi(buf)-1 >= mol->m_waAtomCount[atom])
			{
				eprintf("Error: Only %d %s atoms in the molecule (requested: %d)\n",mol->m_waAtomCount[atom],(const char*)((CAtom*)g_oaAtoms[mol->m_baAtomIndex[atom]])->m_sName,atoi(buf));
				return false;
			}
			if (m)	// Range: Include all atoms from i to i2				case: X1-3	A
			{
				i2 = atoi(buf)-1;
				if (i2 < i) 	// If range from high to low -> wrong input
				{
					eprintf("Error: Invalid atom range, %d < %d.\n",i2+1,i+1);
					return false;
				}
				else
				{
					//>>>MB170422
					//for (inclu = i; inclu < i2; inclu++)
					for (inclu = i; inclu <= i2; inclu++)
					//<<<MB170422
					{
						includeAtom(atom, inclu, mol->m_iIndex, m_laInclude);	 //hier
					}
					fin = true;
				}
			} 
			else	
			{
				i = atoi(buf)-1;
			}
		} 
		else if (*q == '-')
		{
			if (i == -1)
			{
				eprintf("Error: \"-\" without preceding number.\n");
				return false;
			} 
			m = true;
			q++;
		} 
		else if (*q == ',')   
		{
			if (atom == -1)
			{
				eprintf("Error: Comma without atom.\n");
				return false;
			}
			if ( i == -1 )	// include all atoms of kind				case: X,	C
			{
				for (inclu = 0; inclu < mol->m_waAtomCount[atom]; inclu++ )
				{
					includeAtom(atom, inclu, mol->m_iIndex, m_laInclude);	 //hier
				}
			}
			else if (!m)  // include single atom							case: X1,	B
			{
				includeAtom(atom, i, mol->m_iIndex, m_laInclude);	// hier
			}
			la = atom;
			m = false;
			i = -1;
			i2 = -1;
			atom = -1;
			q++;
			fin = true;
		}
		p = q;
	}
	if ( !fin )
	{
		if ( i == -1 )	// case: X
		{
			for (inclu = 0; inclu < mol->m_waAtomCount[atom]; inclu++ )
			{
				includeAtom(atom, inclu, mol->m_iIndex, m_laInclude);	 //hier
			}
		}
		else			// case: X1
		{
			includeAtom(atom, i, mol->m_iIndex, m_laInclude); // hier
		}
	}

	return true;
}


void CHBondObservation::BuildDonorArray()
{
	int n,m,i;

	m_laDonor->RemoveAll();

	for (n=0; n < m_laHydro->GetSize();n++) 	// Walk over all hydrogens
	{
		for (m=0; m < ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(); m++)	// Walk over all SingleMolecules
		{
			for (i=0; i < ((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetAt(m)])->m_oaMolAtoms.GetSize(); i++)	// Walk over all Atoms in Molecule
			{
				if ( ((CMolAtom*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetAt(m)])->m_oaMolAtoms[i])->m_iOffset == m_laHydro->GetAt(n) ) // Search Atom
				{
					if ( ((CMolAtom*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetAt(m)])->m_oaMolAtoms[i])->m_oaBonds.GetSize() > 1 )
					{
						eprintf("Something is wrong: Hydrogen has more than one bond! Aborting ...\n");
						abort();
					}
					m_laDonor->Add(((CMolAtom*)((CMolAtom*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetAt(m)])->m_oaMolAtoms[i])->m_oaBonds[0])->m_iOffset);
					m = ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
					break;
				}
			}
		}
	}
}


void CHBondObservation::BuildHydroArray()
{
	int n,m,i,j;
	bool tb;
	CMolecule *mol;
	CSingleMolecule *sm;
	CMolAtom *ma;

	mol = (CMolecule*)g_oaMolecules[m_iShowMol];

	for (n=0; n < m_laDonor->GetSize();n++) 		// Walk over all donors
	{
		tb = false;
		for (m=0; m < mol->m_laSingleMolIndex.GetSize(); m++)		// Walk over all SingleMolecules
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[mol->m_laSingleMolIndex[m]];
			for (i=0; i < sm->m_oaMolAtoms.GetSize(); i++)	// Walk over all Atoms in Molecule
			{
				ma = (CMolAtom*)sm->m_oaMolAtoms[i];
				if ( ma->m_iOffset == m_laDonor->GetAt(n) ) // Search Atom
				{
					for (j=0; j < ma->m_oaBonds.GetSize(); j++) // Walk over all bondpartners
					{
						if ( ((CAtom*)g_oaAtoms[mol->m_baAtomIndex[((CMolAtom*)ma->m_oaBonds[j])->m_iType]])->m_pElement->m_iOrd == 1 )
						{
							m_laHydro->Add(((CMolAtom*)ma->m_oaBonds[j])->m_iOffset);
							tb = true;
						}
					}
					if (!tb)
					{
						mprintf("Warning: Did not find hydrogen atoms on donor atom %s%d[%d].\n",(const char*)((CAtom*)g_oaAtoms[mol->m_baAtomIndex[ma->m_iType]])->m_sName,g_waAtomMolNumber[((CMolAtom*)((CSingleMolecule*)g_oaSingleMolecules[m])->m_oaMolAtoms[i])->m_iOffset]+1,m+1);
					}
					m = ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
					break;
				}
			}
		}
	}

	mprintf("    BuildHydroArray(): Found %d hydrogen atoms on %d donor atoms.\n\n",m_laHydro->GetSize(),m_laDonor->GetSize());

	if ((m_laHydro->GetSize() == 0) && (m_laDonor->GetSize() != 0)) {
		mprintf("\n");
		mprintf("    Probably you entered the hydrogen atoms as donor atoms. Instead, you have to\n");
		mprintf("    enter the bond partners of the hydrogen atoms as donor atoms (see note above).\n");
		mprintf("    Aborting execution.\n\n");
		abort();
	}
}


bool CHBondObservation::checkCOND(CxDVector3 m_vAccep, CxDVector3 m_vHydro, CxDVector3 m_vDonor, int i)
{
	CxDVector3 m_vrAD, m_vrAH, m_vrHD;
	double rAD, rAH, rHD, winkel;
	double rAD_crit, rAH_crit, winkel_crit;

	m_vrAD = m_vAccep - m_vDonor;

	m_vrAD = FoldVector(m_vrAD);

	rAD = m_vrAD.GetLength();

	if ( m_bHysteresis && m_laLastNeigh->GetAt(i) == 1 ) 	// pair is neighbor -> out_crit if hysteresis
	{
		rAD_crit = m_frAD_out;			
	}
	else
	{
		rAD_crit = m_frAD;
	}

	if ( m_bHysteresis && m_laLastBond->GetAt(i) == 1 )	// pair is hbonded -> out_crit if hysteresis
	{
		rAH_crit = m_frAH_out;
		winkel_crit = m_fwinkel_out;
	}
	else
	{
		rAH_crit = m_frAH;
		winkel_crit = m_fwinkel;
	}

	if ( rAD < rAD_crit ) 		// distance acceptor-donor smaller reference?
	{
		if ( m_laLastNeigh->GetAt(i) == 0 )
		{
			m_laLastNeigh->SetAt(i,1);				// Nachbar
			((CxIntArray*)m_oaNeigh->GetAt(i))->Add(g_iSteps);
		}
		m_vrHD = m_vDonor - m_vHydro;

		m_vrHD = FoldVector(m_vrHD);

		rHD = m_vrHD.GetLength();
		m_vrAH = m_vAccep - m_vHydro;

		m_vrAH = FoldVector(m_vrAH);

		rAH = m_vrAH.GetLength();

		if ( rAH < rHD )		// "acceptor" is in actually donor
		{
			winkel = (rAD*rAD + rAH*rAH - rHD*rHD) / (2 * rAD * rAH);
		}
		else				// "donor" is de facto donor
		{
			winkel = (rAD*rAD + rHD*rHD - rAH*rAH) / (2 * rAD * rHD);
		}
			
		if ( rAH < rAH_crit && rHD < rAH_crit ) 	// distance acceptor-hydrogen smaller reference? Additional: distance donor-hydrogen smaller (for LUZAR)
		{
			if ( winkel > winkel_crit )	// cos(winkel) bigger than reference (or in other words: angle smaller?)
			{
				if ( m_laLastBond->GetAt(i) == 0 )
				{
					m_laLastBond->SetAt(i,1);						  // H-Bond
					((CxIntArray*)m_oaBond->GetAt(i))->Add(g_iSteps);
				}
				return false;
			}
		}
		if ( m_laLastBond->GetAt(i) == 1 )
		{
			if ( m_bLuzar )
			{
				return true;
			}
			m_laLastBond->SetAt(i,0);						  // Kein H-Bond
			((CxIntArray*)m_oaBond->GetAt(i))->Add(g_iSteps);
		}
		return false;
	}
	else
	{
		if ( m_laLastNeigh->GetAt(i) == 1 )
		{
			m_laLastNeigh->SetAt(i,0);						   // Kein Nachbar
			((CxIntArray*)m_oaNeigh->GetAt(i))->Add(g_iSteps);
		}
		if ( m_laLastBond->GetAt(i) == 1 )
		{
			m_laLastBond->SetAt(i,0);						  // Kein H-Bond ... geht nicht, wenn kein Nachbar ... echt nicht .... geht gar nicht ...
			((CxIntArray*)m_oaBond->GetAt(i))->Add(g_iSteps);
		}
		return false;
	}
	return true;	// This point should never be reached ... it's not a nice point ...
}


void CHBondObservation::WriteOutput()
{
	CGrace* gc;
	int n,m;

	mprintf("    Writing \"Correlations.agr\"...\n");
	
	try { gc = new CGrace(); } catch(...) { gc = NULL; }
	if (gc == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	gc->SetRangeX(0,m_iDepth*g_iStride*g_fTimestepLength/1000);
	gc->SetRangeY(0,1);
	gc->MakeTicks();
	gc->SetLabelX("tau / ps");
	gc->SetLabelY("f(tau)");
	gc->CurrentGraph()->m_bInvertXAxis = false;
	gc->CurrentGraph()->m_bLegend = true;

	gc->AddDataset();
	gc->SetDatasetName(0, "c(t)");
	gc->AddDataset();
	gc->SetDatasetName(1, "n(t)");
	for ( n = 0; n < BondDynamic->GetSize(); n++)
	{
		gc->AddXYTupel(0,n*g_fTimestepLength*g_iStride/1000,BondDynamic->GetAt(n));
		gc->AddXYTupel(1,n*g_fTimestepLength*g_iStride/1000,DiffusionDynamic->GetAt(n));
	}

	gc->WriteAgr("Correlations.agr",false);
	gc->WriteCSV("Correlations.csv");
	delete gc;

	mprintf("    Writing \"Fitted.agr\"...\n");

	CalcMinMax();

	try { gc = new CGrace(); } catch(...) { gc = NULL; }
	if (gc == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	gc->SetRangeX(0,m_dGraceMax);
	gc->SetRangeY(0,m_dGraceMax);
	gc->MakeTicks();
	gc->SetLabelX("kf c(t) - kb n(t) [1/ps]");
	gc->SetLabelY("-dc/dt [1/ps]");
	gc->CurrentGraph()->m_bInvertXAxis = false;
	gc->CurrentGraph()->m_bLegend = true;
		
	gc->AddDataset();
	gc->SetDatasetName(0, "Data");

	for ( m=0; m < k_t->GetSize(); m++)
		gc->AddXYTupel(0,(*a)[0] * (*BondDynamic)[m] - (*b)[0] * (*DiffusionDynamic)[m], -(*k_t)[m]);

	gc->WriteCSV("Fitted.csv");
	gc->AddDataset();
	gc->SetDatasetName(1, "Unity");
	gc->AddXYTupel(1,0,0);
	gc->AddXYTupel(1,m_dGraceMax,m_dGraceMax);
	gc->WriteAgr("Fitted.agr",false);
	delete gc;
}


void CHBondObservation::CalcMinMax()
{
	double temp;
	int n;	

	m_dGraceMin =  1E20;
	m_dGraceMax = -1E20;

	for ( n = m_iTrunc / 5; n < BondDynamic->GetSize(); n++)
	{
		temp = ((*a)[0] * (*BondDynamic)[n] - (*b)[0] * (*DiffusionDynamic)[n]);
		if ( m_dGraceMin >= temp )
			m_dGraceMin = temp;
		if ( m_dGraceMax <= temp )
			m_dGraceMax = temp;
	}

	for ( n = m_iTrunc / 5; n < k_t->GetSize(); n++)
	{
		temp = -k_t->GetAt(n);
		if ( m_dGraceMin >= temp )
			m_dGraceMin = temp;
		if ( m_dGraceMax <= temp )
			m_dGraceMax = temp;
	} 
}
