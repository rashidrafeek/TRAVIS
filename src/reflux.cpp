/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

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

#include "reflux.h"
#include "luzar.h"
#include "math.h"

#include "globalvar.h"
#include "maintools.h"


const char *GetRevisionInfo_reflux(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_reflux() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


static CxObArray g_ReFluxObserv;

bool gatherReFlux()
{
	CReFluxObservation *o;

	try { o = new CReFluxObservation(); } catch(...) { o = NULL; }
	if (o == NULL) NewException((double)sizeof(CReFluxObservation),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	g_ReFluxObserv.Add(o);

	return true;
}

bool initializeReFlux() 
{
	int i;

	for(i = 0; i < g_ReFluxObserv.GetSize(); i++) 
	{
		((CReFluxObservation *)g_ReFluxObserv[i])->initialize();
	}
	return true;
}

bool processReFlux(CTimeStep* ts) 
{
	int i;
	for(i = 0; i < g_ReFluxObserv.GetSize(); i++) 
	{
		((CReFluxObservation *)g_ReFluxObserv[i])->process(ts);
	}
	return true;
}

bool finalizeReFlux() 
{
	int i;
	for(i = 0; i < g_ReFluxObserv.GetSize(); i++) 
	{
		((CReFluxObservation *)g_ReFluxObserv[i])->finalize();
	}
	return true;
}

CReFluxObservation::CReFluxObservation()
{
	int i_temp, mol = 0;
	bool b_temp;
	double d_temp;
	int i,j,z;
	CxString buf, buf2, buf_temp;
	CReFluxCond* temp;

	m_iShowMol3 = -1;
	mprintf("    The generelized reactive flux analysis is organized as followed:\n");
	mprintf("      First two types of molecules are defined: the reference molecule RM\n");
	mprintf("                                                the observed molecule OM\n");
	mprintf("      All possible combinations of these molecules are the total number of pairs.\n");
	mprintf("      The neighboring of each pair is defined by a single distance criterion \n");
	mprintf("      defined by a maximum distance between one atom of RM and one of OM.\n");
	mprintf("      For the aggregation a third molecule type may be include: The TM.\n");
	if (g_oaMolecules.GetSize() > 1)
	{
		buf.sprintf("\n    Which of the molecules should be the reference molecule (");
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
			buf.strcat(buf2);
			if (z < g_oaMolecules.GetSize()-1)
				buf.strcat(", ");
		}
		buf.strcat(")? ");
		m_iShowMol = AskRangeInteger_ND("%s",1,g_oaMolecules.GetSize(),(const char*)buf) - 1;

		buf.sprintf("\n    Which of the molecules should be the observed molecule (");
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
			buf.strcat(buf2);
			if (z < g_oaMolecules.GetSize()-1)
				buf.strcat(", ");
		}
		buf.strcat(")? ");
		m_iShowMol2 = AskRangeInteger_ND("%s",1,g_oaMolecules.GetSize(),(const char*)buf) - 1;

	} else 
	{
		m_iShowMol = 0; 	
		m_iShowMol2 = 0;		
	}
	mprintf(WHITE,"\n    %s is the reference molecule.\n\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	mprintf(WHITE,"\n    %s is the observed molecule.\n\n",((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName);
	for ( i=0; i < ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(); i++)
	{
		for ( j=0; j < ((CMolecule*)g_oaMolecules[m_iShowMol2])->m_laSingleMolIndex.GetSize(); j++)
		{
			try { temp = new CReFluxCond(); } catch(...) { temp = NULL; }
			if (temp == NULL) NewException((double)sizeof(CReFluxCond),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			m_oaCond.Add(temp);
		}
	}
	mprintf("\n    Analyzing %d pairs of molecules!\n", m_oaCond.GetSize());

	mprintf("    Which atom from %s defines the neighboring condition (e.g. \"C1,O3,#2\" )?",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	inpprintf("! Which atom from %s defines the neighboring condition (e.g. \"C1,O3,#2\" )?\n",((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName);
	myget(&buf);
	while ((strlen(buf) == 0) || (!checkAtomChoice(buf,1,1)))
	{
		mprintf("    Please enter single atom!\n");
		myget(&buf);
	}

	mprintf("    Which atom from %s defines the neighboring condition (e.g. \"C1,O3,#2\" )?",((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName);
	inpprintf("! Which atom from %s defines the neighboring condition (e.g. \"C1,O3,#2\" )?\n",((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName);
	myget(&buf);
	while ((strlen(buf) == 0) || (!checkAtomChoice(buf,2,2)))
	{
		mprintf("    Please enter atom!\n");
		myget(&buf);
	}
	d_temp = AskFloat("    Enter maximum distance for neighbor criterion [300] ",300.0);
	for ( i=0; i < m_oaCond.GetSize(); i++ )
		((CReFluxCond*)m_oaCond[i])->m_dNeighDist = d_temp;

	if ( AskYesNo("    Perform three body analysis? [no]", false) )
	{
		if (g_oaMolecules.GetSize() > 1)
		{
			buf.sprintf("\n    Which of the molecules should be the third molecule (");
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
				buf.strcat(buf2);
				if (z < g_oaMolecules.GetSize()-1)
					buf.strcat(", ");
			}
			buf.strcat(")? ");
			m_iShowMol3 = AskRangeInteger_ND("%s",1,g_oaMolecules.GetSize(),(const char*)buf) - 1;
		}
		else
			m_iShowMol3 = 0;
		mprintf(WHITE,"\n    %s is the third molecule.\n\n",((CMolecule*)g_oaMolecules[m_iShowMol3])->m_sName);
	} 
	if ( m_iShowMol3 != -1 )
	{
		mprintf("    For a three body analysis a distance criterion containing at least one atom of the TM is needed!\n");
		mprintf("    Add which atom(s) from %s to distance condition (e.g. \"C1,O3-5,#2\" *=all)?",((CMolecule*)g_oaMolecules[m_iShowMol3])->m_sName);
		inpprintf("! Add which atom(s) from %s to distance condition (e.g. \"C1,O3-5,#2\" *=all)?\n",((CMolecule*)g_oaMolecules[m_iShowMol3])->m_sName);
		while ((strlen(buf) == 0) || (!checkAtomChoice(buf,3,3)))
		{
			mprintf("    Please enter atom(s)!\n");
			myget(&buf);
		}
		i_temp = AskRangeInteger_ND("\n    Take second atom(s) from which molecule (%s(RM)=1, %s(OM)=2, %s(TM)=3 )?",
				1, 3, 
				((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,
				((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName,
				((CMolecule*)g_oaMolecules[m_iShowMol3])->m_sName );
		if 	( i_temp == 1 )
			mol = m_iShowMol;
		else if	( i_temp == 2 )
			mol = m_iShowMol2;
		else if	( i_temp == 3 )
			mol = m_iShowMol3;
		mprintf("    Add which atom(s) from %s to distance condition (e.g. \"C1,O3-5,#2\" *=all)?",((CMolecule*)g_oaMolecules[mol])->m_sName);
		inpprintf("! Add which atom(s) from %s to distance condition (e.g. \"C1,O3-5,#2\" *=all)?\n",((CMolecule*)g_oaMolecules[mol])->m_sName);
		myget(&buf);
		while ((strlen(buf) == 0) || (!checkAtomChoice(buf,i_temp,4)))
		{
			mprintf("    Please enter atom(s)!\n");
			myget(&buf);
		}
		d_temp = AskFloat("    Enter maximum distance for distance criterion (in pm) [300] ",300.0);
		for ( i=0; i < m_oaCond.GetSize(); i++ )
			for ( z=0; z < ((CReFluxCond*)m_oaCond[i])->m_oaDistConds.GetSize(); z++ )
				if ( ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[i])->m_oaDistConds[z])->m_dAggrMinAngle == -1.0 )
					((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[i])->m_oaDistConds[z])->m_dAggrMinAngle = d_temp;
	}	 
	while ( AskYesNo("    Define a(nother) distance criterion for the aggregation? [no]", false) )
	{
		if ( m_iShowMol3 != -1 )
			i_temp = AskRangeInteger_ND("\n    Take first atom(s) from which molecule (%s(RM)=1, %s(OM)=2, %s(TM)=3)?",
				1, 3, 
				((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,
				((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName,
				((CMolecule*)g_oaMolecules[m_iShowMol3])->m_sName );
		else
			i_temp = AskRangeInteger_ND("\n    Take first atom(s) from which molecule (%s(RM)=1, %s(OM)=2)?",
				1, 2, 
				((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,
				((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName );
		if 	( i_temp == 1 )
			mol = m_iShowMol;
		else if	( i_temp == 2 )
			mol = m_iShowMol2;
		else if	( i_temp == 3 )
			mol = m_iShowMol3;
		mprintf("    Add which atom(s) from %s to distance condition (e.g. \"C1,O3-5,#2\" *=all)?",((CMolecule*)g_oaMolecules[mol])->m_sName);
		inpprintf("! Add which atom(s) from %s to distance condition (e.g. \"C1,O3-5,#2\" *=all)?\n",((CMolecule*)g_oaMolecules[mol])->m_sName);
		myget(&buf);
		while ((strlen(buf) == 0) || (!checkAtomChoice(buf,i_temp,3)))
		{
			mprintf("    Please enter atom(s)!\n");
			myget(&buf);
		}
		if ( m_iShowMol3 != -1 )
			i_temp = AskRangeInteger_ND("\n    Take second atom(s) from which molecule (%s(RM)=1, %s(OM)=2, %s(TM)=3)?",
				1, 3, 
				((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,
				((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName,
				((CMolecule*)g_oaMolecules[m_iShowMol3])->m_sName );
		else
			i_temp = AskRangeInteger_ND("\n    Take second atom(s) from which molecule (%s(RM)=1, %s(OM)=2)?",
				1, 2, 
				((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,
				((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName );
		if 	( i_temp == 1 )
			mol = m_iShowMol;
		else if	( i_temp == 2 )
			mol = m_iShowMol2;
		else if	( i_temp == 3 )
			mol = m_iShowMol3;
		mprintf("    Add which atom(s) from %s to distance condition (e.g. \"C1,O3-5,#2\" *=all)?",((CMolecule*)g_oaMolecules[mol])->m_sName);
		inpprintf("! Add which atom(s) from %s to distance condition (e.g. \"C1,O3-5,#2\" *=all)?\n",((CMolecule*)g_oaMolecules[mol])->m_sName);
		myget(&buf);
		while ((strlen(buf) == 0) || (!checkAtomChoice(buf,i_temp,4)))
		{
			mprintf("    Please enter atom(s)!\n");
			myget(&buf);
		}
		d_temp = AskFloat("    Enter maximum distance for distance criterion (in pm) [300] ",300.0);
		for ( i=0; i < m_oaCond.GetSize(); i++ )
			for ( z=0; z < ((CReFluxCond*)m_oaCond[i])->m_oaDistConds.GetSize(); z++ )
				if ( ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[i])->m_oaDistConds[z])->m_dAggrMinAngle == -1.0 )
					((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[i])->m_oaDistConds[z])->m_dAggrMinAngle = d_temp;
	}	
	SortDist();

	while ( AskYesNo("\n    Define an(other) angle criterion for the aggregation? [no]", false) )
	{
		mprintf("    Two vectors need to be defined.\n");
		mprintf("    Each vector can be defined by two points or perpendicular to a plane by three points.\n");
		mprintf("    Each point can be an atom from the reference molecule (RM), the observed molecule (OM),\n");
		mprintf("    or the third molecule (TM).\n");

		for ( z=5; z<11; z++ )
		{
			if ( z == 5 )
				mprintf("\n    First vector:\n");
			else
				mprintf("\n    Second vector:\n");
			b_temp = !AskYesNo("    Define vector by two points (y) or perpendicular to a plane (n) ? [yes]", true);

			if ( m_iShowMol3 != -1 )
				i_temp = AskRangeInteger_ND("\n    Take base atom(s) from which molecule (%s(RM)=1, %s(OM)=2, %s(TM)=3)?",
					1, 3, 
					((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,
					((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName,
					((CMolecule*)g_oaMolecules[m_iShowMol3])->m_sName );
			else
				i_temp = AskRangeInteger_ND("\n    Take base atom(s) from which molecule (%s(RM)=1, %s(OM)=2)?",
					1, 2, 
					((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,
					((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName );
			if 	( i_temp == 1 )
				mol = m_iShowMol;
			else if	( i_temp == 2 )
				mol = m_iShowMol2;
			else if	( i_temp == 3 )
				mol = m_iShowMol3;
			mprintf("    Use which atom(s) from %s as base atom(s) (e.g. \"C1,O3-5,#2\" *=all)?",((CMolecule*)g_oaMolecules[mol])->m_sName);
			inpprintf("! Use which atom(s) from %s as base atom(s) (e.g. \"C1,O3-5,#2\" *=all)?\n",((CMolecule*)g_oaMolecules[mol])->m_sName);
			myget(&buf);
			while ((strlen(buf) == 0) || (!checkAtomChoice(buf, i_temp, z++)))
			{
				mprintf("    Please enter atom(s)!\n");
				myget(&buf);
			}
			
			if ( m_iShowMol3 != -1 )
				i_temp = AskRangeInteger_ND("\n    Take second atom(s) from which molecule (%s(RM)=1, %s(OM)=2, %s(TM)=3)?",
					1, 3, 
					((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,
					((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName,
					((CMolecule*)g_oaMolecules[m_iShowMol3])->m_sName );
			else
				i_temp = AskRangeInteger_ND("\n    Take second atom(s) from which molecule (%s(RM)=1, %s(OM)=2)?",
					1, 2, 
					((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,
					((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName );
			if 	( i_temp == 1 )
				mol = m_iShowMol;
			else if	( i_temp == 2 )
				mol = m_iShowMol2;
			else if	( i_temp == 3 )
				mol = m_iShowMol3;
			mprintf("    Use which atom(s) from %s as second atom(s) (e.g. \"C1,O3-5,#2\" *=all)?",((CMolecule*)g_oaMolecules[mol])->m_sName);
			inpprintf("! Use which atom(s) from %s as second atom(s) (e.g. \"C1,O3-5,#2\" *=all)?\n",((CMolecule*)g_oaMolecules[mol])->m_sName);
			myget(&buf);
			while ((strlen(buf) == 0) || (!checkAtomChoice(buf, i_temp, z++)))
			{
				mprintf("    Please enter atom(s)!\n");
				myget(&buf);
			}
			if ( b_temp )
			{
				if ( m_iShowMol3 != -1 )
					i_temp = AskRangeInteger_ND("\n    Take third atom(s) from which molecule (%s(RM)=1, %s(OM)=2, %s(TM)=3)?",
						1, 3, 
						((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,
						((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName,
						((CMolecule*)g_oaMolecules[m_iShowMol3])->m_sName );
				else
					i_temp = AskRangeInteger_ND("\n    Take third atom(s) from which molecule (%s(RM)=1, %s(OM)=2)?",
						1, 2, 
						((CMolecule*)g_oaMolecules[m_iShowMol])->m_sName,
						((CMolecule*)g_oaMolecules[m_iShowMol2])->m_sName );
				if 	( i_temp == 1 )
					mol = m_iShowMol;
				else if	( i_temp == 2 )
					mol = m_iShowMol2;
				else if	( i_temp == 3 )
					mol = m_iShowMol3;
				mprintf("    Use which atom(s) from %s as third atom(s) (e.g. \"C1,O3-5,#2\" *=all)?",((CMolecule*)g_oaMolecules[mol])->m_sName);
				inpprintf("! Use which atom(s) from %s as third atom(s) (e.g. \"C1,O3-5,#2\" *=all)?\n",((CMolecule*)g_oaMolecules[mol])->m_sName);
				myget(&buf);
				while ((strlen(buf) == 0) || (!checkAtomChoice(buf, i_temp, z)))
				{
					mprintf("    Please enter atom(s)!\n");
					myget(&buf);
				}
			}
			else
				for ( i=0; i < m_oaCond.GetSize(); i++ )
					for ( j=0; j < ((CReFluxCond*)m_oaCond[i])->m_oaAngleConds.GetSize(); j++ )
						((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[i])->m_oaAngleConds[j])->m_iaTypes.Add(-1);
		}
		d_temp = cos((M_PI/180)*AskRangeFloat("    Enter minimum angle for aggregation criterion [0] ",0.0,179.0,0.0));
		for ( i=0; i < m_oaCond.GetSize(); i++ )
			for ( j=0; j < ((CReFluxCond*)m_oaCond[i])->m_oaAngleConds.GetSize(); j++ )
				if ( ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[i])->m_oaAngleConds[j])->m_dAggrMinAngle == -1 )
					((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[i])->m_oaAngleConds[j])->m_dAggrMinAngle = d_temp;
		d_temp = cos((M_PI/180)*AskRangeFloat("    Enter maximum angle for aggregation criterion [180] ",d_temp,180.0,180.0));
		for ( i=0; i < m_oaCond.GetSize(); i++ )
			for ( j=0; j < ((CReFluxCond*)m_oaCond[i])->m_oaAngleConds.GetSize(); j++ )
				if ( ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[i])->m_oaAngleConds[j])->m_dAggrMaxAngle == -1 )
					((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[i])->m_oaAngleConds[j])->m_dAggrMaxAngle = d_temp;
	}
	SortAngle();
	SortFinal();

	if ( ((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize() == 0 ) 
	{
		eprintf("Error: No aggregation criterion left?\n");
		abort();
	}

	z = m_oaCond.GetSize();
	j = ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
	if ( m_iShowMol == m_iShowMol2 )
		for ( i=0; i < j; i++ )
			m_oaCond.RemoveAt(z-(i*(j+1))-1,1);

	mprintf("\n    Analysing %d pairs with %d aggregation criteria each ...\n",m_oaCond.GetSize(),((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize());
	g_fTimestepLength = AskFloat("    Enter the physical time distance between successive trajectory frames in fs: [0.5] ",0.5f);
}

CReFluxObservation::~CReFluxObservation()
{
}

void CReFluxObservation::initialize()
{
	CxIntArray* temp;
	int i;

	for ( i=0; i < m_oaCond.GetSize(); i++ )
	{
		try { temp = new CxIntArray(); } catch(...) { temp = NULL; }
		if (temp == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		m_oaNeigh.Add(temp);

		try { temp = new CxIntArray(); } catch(...) { temp = NULL; }
		if (temp == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		m_oaAggr.Add(temp);
	}

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
}

void CReFluxObservation::process(CTimeStep* ts)
{
	int i;

	for ( i=0; i < m_oaCond.GetSize(); i++ )
		((CReFluxCond*)m_oaCond[i])->check(ts,((CxIntArray*)m_oaNeigh[i]),((CxIntArray*)m_oaAggr[i]));
}

void CReFluxObservation::finalize()	// reconstruct functions -> autocorrelate -> sum results -> normalize -> fit
{
	long n,i,interval;
	double AggrCor, DiffCor;
	CxString buf,buf2;
	unsigned long t0,eta;

	CxDoubleArray NeighFunction;	// Input
	CxDoubleArray AggrFunction;
	
	CxDoubleArray temp1Function;	// Temp
	CxDoubleArray temp2Function;

	CLuzarCorrelation luzar;

	for ( i=0; i < m_iDepth; i++ )
	{
		AggrDynamic.Add(0.0);
		DiffusionDynamic.Add(0.0);
	}

	NeighFunction.SetSize(g_iSteps / g_iStride);
	AggrFunction.SetSize(g_iSteps / g_iStride);

	luzar.Init(g_iSteps / g_iStride,m_iDepth); 	// Inititalize Luzar-Correlation

	mprintf("    Finishing %6i Functions ... \nFunc %6li ",m_oaNeigh.GetSize(),(long)0);
	g_iDotCounter = 0;
	t0 = (unsigned long)time(NULL);
	interval = 2;

	for ( n=0; n < m_oaNeigh.GetSize(); n++)		// Walk over all pairs
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
					eta = (unsigned long)(((double)time(NULL) - t0) / n * ((double)MAX((long)0,(long)m_oaNeigh.GetSize() - (long)n)));
					FormatTime(eta,&buf);
					mprintf(" ETA %s",(const char*)buf);
				}
				mprintf("\nFunc %6li ",n);
			}
			g_iDotCounter++;
			mprintf(".");
		}

		if ( ((CxIntArray*)m_oaNeigh[n])->GetSize() > 0 )
		{
			RestoreFunction(((CxIntArray*)m_oaNeigh[n]),&NeighFunction);	// Restore function of pair n as NeighFunction and AggrFunction
			RestoreFunction(((CxIntArray*)m_oaAggr[n]),&AggrFunction);	

			luzar.Correlate(&AggrFunction,&NeighFunction,&temp1Function,&temp2Function);

			for ( i=0; i < m_iDepth; i++ )		// Sum results to Output-array
			{
				AggrDynamic[i] += temp1Function[i];
				DiffusionDynamic[i] += temp2Function[i];
			}
		}
	}
	mprintf("\n\n");

	if ( AggrDynamic[0] == 0 || DiffusionDynamic[0] == 0 )
	{
		eprintf("First correlation value is zero ... something is wrong ...\n");
		abort();
	}

	AggrCor = (double)AggrDynamic[0];
	DiffCor = (double)DiffusionDynamic[0];

	for ( i=0; i < m_iDepth; i++ )			// Normalize
	{
		AggrDynamic[i] /= AggrCor;
		DiffusionDynamic[i] /= DiffCor;
		DiffusionDynamic[i] -= AggrDynamic[i];
	}

	luzar.Fit(&AggrDynamic, &DiffusionDynamic, &k_t, &a, &b, m_iTrunc);

	WriteOutput();
}


void CReFluxObservation::RestoreFunction(CxIntArray* input, CxDoubleArray* output)
{
	unsigned int n;
	bool eins = false;

	for ( n=0; n < (g_iSteps / g_iStride); n++ )		// Walk over all processed steps
	{
		if ( input->Contains((n+1)*g_iStride) )
			eins = !eins;
		if (!eins)
			output->SetAt(n,0.0);
		else
			output->SetAt(n,1.0);
	}
}

bool CReFluxObservation::includeAtom(int atom, int number, int moltype, int target)	// atom: element number, number: atomsort number, mol m_iShowMol, target: in which array to include
{
	int n, m, i, smol = -1, count = -1, mol;
	CxIntArray includeME, temp;

	if 	( moltype == 1 )
		mol = m_iShowMol;
	else if ( moltype == 2 )
		mol = m_iShowMol2;
	else if ( moltype == 3 )
		mol = m_iShowMol3;
	else
	{
		eprintf("CReFluxObservation::includeAtom Error: Unknown moltype %d!",moltype);
		abort();
	}

	for ( n=0; n < ((CMolecule*)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize(); n++)
	{
		for ( m=0; m < ((CMolecule*)g_oaMolecules[m_iShowMol2])->m_laSingleMolIndex.GetSize(); m++)
		{
			count++;
			if      ( moltype == 1 )		// Ref 
				smol = n;
			if	( moltype == 2 )		// Obs 
				smol = m;
			
			if ( smol != -1 )
			{
				includeME.CopyFrom((CxIntArray*)(((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[mol])->m_laSingleMolIndex[smol]])->m_oaAtomOffset[atom]));
						
				if (number >= includeME.GetSize()) 
				{
					eprintf("CReFluxObservation::includeAtom(): Internal error (%d/%d).\n",number,includeME.GetSize());
					abort();
				}
				((CReFluxCond*)m_oaCond[count])->Add(moltype, target,includeME[number]);
			}
			else
			{
				for ( i=0; i < ((CMolecule*)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize(); i++ )
				{
					temp.CopyFrom((CxIntArray*)(((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[mol])->m_laSingleMolIndex[smol]])->m_oaAtomOffset[atom]));

					if (number >= temp.GetSize()) 
					{
						eprintf("CReFluxObservation::includeAtom(): Internal error (%d/%d).\n",number,temp.GetSize());
						abort();
					}
					includeME.Add(temp[number]);
				}
				((CReFluxCond*)m_oaCond[count])->AddArray(moltype, target,&includeME);
			}	
		}
	}
	return true;
}

bool CReFluxObservation::checkAtomChoice(const char *s, int moltype, int target)
{
	const char *p, *q;
	char buf[32];
	const char *allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890-,#*_ ";
	int inclu, la = -1, i = -1, i2 = -1, atom = -1, n, mol; 
	bool m = false, fin = false, single = false;

	if 	( moltype == 1 )
		mol = m_iShowMol;
	else if ( moltype == 2 )
		mol = m_iShowMol2;
	else if ( moltype == 3 )
		mol = m_iShowMol3;
	else
	{
		eprintf("CReFluxObservation::includeAtom Error: Unknown moltype %d!",moltype);
		abort();
	}

	if ( target == 1 || target == 2 )
		single = true;

	p = s;

	while (*p != 0)
	{
		fin = false;
		while (*p == ' ')
			p++;
		if ( single && ( *p == ',' || *p == '*' || *p == '-' ) )
		{
			eprintf("Error: Must be single atom!\n");
			return false;
		}
		if (strchr(allowed,*p) == NULL)
		{
			eprintf("Error: Character \"%c\" not allowed.\n",*p);
			return false;
		}
		if ( *p == '*' )
		{
			for ( n=0; n < ((CMolecule*)g_oaMolecules[mol])->m_baAtomIndex.GetSize()-1; n++ )
				for (inclu = 0; inclu < ((CMolecule*)g_oaMolecules[mol])->m_waAtomCount[n]; inclu++ )
					includeAtom(n, inclu, moltype, target);	 //hier
			return true;
		}
		q = p;
		if (isalpha(*q) || *q == '#')
		{
			if (m)
			{
				eprintf("Error: Only digit allowed after \"-\".\n");
				return false;
			}
			while ( isalpha(*q) || *q == '#' )
				q++;
			if (q-p >= 32)
			{
				eprintf("Internal Error (%ld >= 32).\n",(long)(q-p));
				return false;
			}
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
			atom = ((CMolecule*)g_oaMolecules[mol])->FindAtomInMol(buf); 	// atom holds number of element 
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
			if (atoi(buf)-1 >= ((CMolecule*)g_oaMolecules[mol])->m_waAtomCount[atom])
			{
				eprintf("Error: Only %d %s atoms in the molecule (requested: %d)\n",((CMolecule*)g_oaMolecules[mol])->m_waAtomCount[atom],(const char*)((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[mol])->m_baAtomIndex[atom]])->m_sName,atoi(buf));
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
					for (inclu = i; inclu <= i2; inclu++)
					{
						includeAtom(atom, inclu, moltype, target);	 //hier
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
				for (inclu = 0; inclu < ((CMolecule*)g_oaMolecules[mol])->m_waAtomCount[atom]; inclu++ )
					includeAtom(atom, inclu, moltype, target);	 //hier
			else if (!m)  // include single atom							case: X1,	B
				includeAtom(atom, i, moltype, target);	// hier
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
			for (inclu = 0; inclu < ((CMolecule*)g_oaMolecules[mol])->m_waAtomCount[atom]; inclu++ )
			{
				includeAtom(atom, inclu, moltype, target);	 //hier
			}
		}
		else			// case: X1
		{
			includeAtom(atom, i, moltype, target); // hier
		}
	}

	return true;
}

void CReFluxObservation::WriteOutput()
{
	CGrace gc, gc2;
	int n,m;

	mprintf("    Writing \"Correlations.agr\"...\n");
	
	gc.SetRangeX(0,m_iDepth*g_iStride*g_fTimestepLength/1000);
	gc.SetRangeY(0,1);
	gc.MakeTicks();
	gc.SetLabelX("t / ps");
	gc.SetLabelY("f(t)");
	gc.CurrentGraph()->m_bInvertXAxis = false;
	gc.CurrentGraph()->m_bLegend = true;

	gc.AddDataset();
	gc.SetDatasetName(0, "c(t)");
	gc.AddDataset();
	gc.SetDatasetName(1, "n(t)");
	for ( n = 0; n < AggrDynamic.GetSize(); n++)
	{
		gc.AddXYTupel(0,n*g_fTimestepLength*g_iStride/1000,AggrDynamic[n]);
		gc.AddXYTupel(1,n*g_fTimestepLength*g_iStride/1000,DiffusionDynamic[n]);
	}

	gc.WriteAgr("Correlations.agr",false);
	gc.WriteCSV("Correlations.csv");

	mprintf("    Writing \"Fitted.agr\"...\n");

	CalcMinMax();

	gc2.SetRangeX(0,m_dGraceMax);
	gc2.SetRangeY(0,m_dGraceMax);
	gc2.MakeTicks();
	gc2.SetLabelX("kf c(t) - kb n(t) [1/ps]");
	gc2.SetLabelY("-dc/dt [1/ps]");
	gc2.CurrentGraph()->m_bInvertXAxis = false;
	gc2.CurrentGraph()->m_bLegend = true;
		
	gc2.AddDataset();
	gc2.SetDatasetName(0, "Data");

	for ( m=0; m < k_t.GetSize(); m++)
		gc2.AddXYTupel(0,a[0] * AggrDynamic[m] - b[0] * DiffusionDynamic[m], - k_t[m]);

	gc2.WriteCSV("Fitted.csv");
	gc2.AddDataset();
	gc2.SetDatasetName(1, "Unity");
	gc2.AddXYTupel(1,0,0);
	gc2.AddXYTupel(1,m_dGraceMax,m_dGraceMax);
	gc2.WriteAgr("Fitted.agr",false);
}

void CReFluxObservation::CalcMinMax()
{
	double temp;
	int n;	

	m_dGraceMin =  1E20;
	m_dGraceMax = -1E20;

	for ( n = m_iTrunc / 5; n < AggrDynamic.GetSize(); n++)
	{
		temp = (a[0] * AggrDynamic[n] - b[0] * DiffusionDynamic[n]);
		if ( m_dGraceMin >= temp )
			m_dGraceMin = temp;
		if ( m_dGraceMax <= temp )
			m_dGraceMax = temp;
	}

	for ( n = m_iTrunc / 5; n < k_t.GetSize(); n++)
	{
		temp = -k_t.GetAt(n);
		if ( m_dGraceMin >= temp )
			m_dGraceMin = temp;
		if ( m_dGraceMax <= temp )
			m_dGraceMax = temp;
	} 
}

CReFluxCond::CReFluxCond()
{
}

CReFluxCond::~CReFluxCond()
{
}

void CReFluxCond::check(CTimeStep* ts, CxIntArray* na, CxIntArray* aa)
{
	int i;
	CxDVector3 Vec;

	if ( m_iRefNeigh == m_iObsNeigh )
		return;
	Vec = FoldVector( ts->m_vaCoords[m_iRefNeigh] - ts->m_vaCoords[m_iObsNeigh] );
	if ( Vec.GetLength() > m_dNeighDist )
	{
		if ( na->GetSize() % 2 == 1 )
			na->Add(g_iSteps);
		if ( aa->GetSize() % 2 == 1 )
			aa->Add(g_iSteps);
		return;
	}
	if ( na->GetSize() % 2 == 0 )
		na->Add(g_iSteps);
	
	for ( i=0; i < m_oaAggrConds.GetSize(); i++ )
	{
		if ( ((CReFluxSubCond*)m_oaAggrConds[i])->check(ts) )
		{
			if ( aa->GetSize() % 2 == 0 )
				aa->Add(g_iSteps);
			return;
		}
	}
	if ( aa->GetSize() % 2 == 1 )
		aa->Add(g_iSteps);
}

void CReFluxCond::Add(int moltype, int target, int number)
{
	CReFluxSingleCond* sc;
	int i, i_temp, j;

	if 	( target == 1 )		// Neighbor 1
		m_iRefNeigh = number;
	else if ( target == 2 )		// Neighbor 2
		m_iObsNeigh = number;
	else if ( target == 3 || target == 5 )		// Dist 1 or Angle 1
	{
		try { sc = new CReFluxSingleCond(); } catch(...) { sc = NULL; }
                if (sc == NULL) NewException((double)sizeof(CReFluxSingleCond),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		if ( target == 3 )
			sc->m_iRefDist = number;
		if ( target == 5 )
			sc->m_iRefVec = number;
		sc->m_iaTypes.Add(moltype);

		if ( target == 3 )
			m_oaDistConds.Add(sc);
		if ( target == 5 )
			m_oaAngleConds.Add(sc);
	}
	else if ( target == 4 || ( target > 5 && target < 11 ) )		// Dist 2 or Angle 2
	{
		j = -1;
		i_temp = m_oaAngleConds.GetSize();
		if ( target == 4 )
			i_temp = m_oaDistConds.GetSize();
		
		for ( i=0; i < i_temp; i++ )
		{
			if ( target == 4 )
				sc = (CReFluxSingleCond*)m_oaDistConds[i];
			else
				sc = (CReFluxSingleCond*)m_oaAngleConds[i];
			if ( sc->m_dAggrMinAngle != -1 )
				continue;
			if ( j == -1 )
				j = i;
			if ( target == 4 && i != j && sc->m_iRefDist == ((CReFluxSingleCond*)m_oaDistConds[j])->m_iRefDist )
				break;
			if (       (target ==  4 && sc->m_iObsVec  != -1)
				|| (target ==  6 && sc->m_iRefVec2 != -1)
				|| (target ==  7 && sc->m_iRefVec3 != -1)
				|| (target ==  8 && sc->m_iObsVec  != -1)
				|| (target ==  9 && sc->m_iObsVec2 != -1)
				|| (target == 10 && sc->m_iObsVec3 != -1) )
			{
				try { sc = new CReFluxSingleCond(); } catch(...) { sc = NULL; }
                		if (sc == NULL) NewException((double)sizeof(CReFluxSingleCond),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				if ( target == 4 )
				{
					sc->Copy( (CReFluxSingleCond*)m_oaDistConds[i] );
					m_oaDistConds.Add(sc);
				}
				else
				{
					sc->Copy( (CReFluxSingleCond*)m_oaAngleConds[i] );
					m_oaAngleConds.Add(sc);
				}
				sc->m_iaTypes.RemoveAt(sc->m_iaTypes.GetSize()-1,1);
			}
			if 	( target ==  4 || target == 8 )
				sc->m_iObsVec = number;
			else if ( target ==  6 )
				sc->m_iRefVec2 = number;
			else if ( target ==  7 )
				sc->m_iRefVec3 = number;
			else if ( target ==  9 )
				sc->m_iObsVec2 = number;
			else if ( target == 10 )
				sc->m_iObsVec3 = number;
			sc->m_iaTypes.Add(moltype);
		}
	}
	else
	{
		eprintf("CReFluxCond::Add(): Internal error (target %d).\n",target);
		abort();
	}
}

void CReFluxCond::AddArray(int moltype, int target, CxIntArray* addme)
{
	int i;

	if	( target == 1 || target == 2 )
	{
		eprintf("CReFluxCond::AddArray(): Internal error (target %d).\n",target);
		abort();
	}
	for ( i=0; i < addme->GetSize(); i++ )
		Add(moltype, target, (*addme)[i]);
}

bool CReFluxObservation::GetArray(const char* s, CxIntArray* temp, int max)
{
	const char *p, *q;
	char buf[32];

	temp->RemoveAll();

	p = s;

	while ( *p != 0 )
	{
		while ( *p == ' ' || *p == ',' )
			p++;
		q = p;
		while ( isdigit(*q) )
			q++;
		if ( q == p )
		{
			mprintf("\n    Please enter comma separated numbers!\n");
			return false;
		}
		memcpy(buf,p,q-p);
		buf[q-p] = 0;
		if ( atoi(buf) > max )
		{
			mprintf("\n    There is no option %d!\n", atoi(buf));
			return false;
		}
		temp->Add(atoi(buf)-1);
		p = q;
	}	

	return true;
}

bool CReFluxObservation::GetDoubleArray( const char* s, CxIntArray* temp1, CxIntArray* temp2, int max1, int max2 )
{
	const char *p, *q;
	char buf[32];
	bool minus = false;

	temp1->RemoveAll();
	temp2->RemoveAll();

	p = s;

	while ( *p != 0 )
	{
		while ( *p == ' ' || *p == ',' )
			p++;
		if ( *p == '-' )
		{
			minus = true;
			p++;
		}
		q = p;
		while ( isdigit(*q) )
			q++;
		if ( q == p )
		{
			mprintf("\n    Please enter: x1,x2-y1,y2,y3 !\n");
			return false;
		}
		memcpy(buf,p,q-p);
		buf[q-p] = 0;
		if ( ( atoi(buf) > max1 && !minus ) || ( atoi(buf) > max2 && minus ) )
		{
			if ( !minus )
				mprintf("\n    There is no option %d for distances!\n", atoi(buf));
			else
				mprintf("\n    There is no option %d for distances!\n", atoi(buf));
			return false;
		}
		if ( !minus )
			temp1->Add(atoi(buf)-1);
		else
			temp2->Add(atoi(buf)-1);
		p = q;
	}	

	return true;
}

CxString CReFluxObservation::FindType(CReFluxSingleCond* sc, int k)
{
	CxString out;

	if 	( sc->m_iaTypes[k] == 1 )
		out = ("RM");
	else if ( sc->m_iaTypes[k] == 2 )
		out = ("OM");
	else if ( sc->m_iaTypes[k] == 3 )
		out = ("TM");

	return out;
}

void CReFluxObservation::DumpDist()
{
	int i;
	int a1, a2;
	CxString out, temp;

	for ( i=0; i < ((CReFluxCond*)m_oaCond[0])->m_oaDistConds.GetSize(); i++ )
	{
		a1 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaDistConds[i])->m_iRefDist;
		a2 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaDistConds[i])->m_iObsVec;

		out.Format("%7d: %s%d (", i+1,
			(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[a1]])->m_sName,
			g_waAtomMolNumber[a1]+1 );
		out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaDistConds[i], 0);
		temp.Format(") - %s%d (",
			(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[a2]])->m_sName,
			g_waAtomMolNumber[a2]+1 );
		out += temp;
		out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaDistConds[i], 1);
		out += (")\n");
		out.Dump();
	}
}

void CReFluxObservation::DumpAngle()
{
	int i;
	int a1, a2, a3, b1, b2, b3;
	CxString out, temp;

	for ( i=0; i < ((CReFluxCond*)m_oaCond[0])->m_oaAngleConds.GetSize(); i++ )
	{
		a1 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i])->m_iRefVec;
		a2 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i])->m_iRefVec2;
		a3 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i])->m_iRefVec3;
		b1 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i])->m_iObsVec;
		b2 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i])->m_iObsVec2;
		b3 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i])->m_iObsVec3;

		out.Format("%7d: ", i+1);
		if ( a3 == -1 )
			out += ("         vec[");
		else
			out += ("plane[");
		temp.Format("%s%d(",
			(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[a1]])->m_sName,
			g_waAtomMolNumber[a1]+1 );
		out += temp;
		out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i], 0);
		temp.Format(")-%s%d(",
			(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[a2]])->m_sName,
			g_waAtomMolNumber[a2]+1 );
		out += temp;
		out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i], 1);
		if ( a3 != -1 )
		{
			temp.Format(")-%s%d(",
				(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[a3]])->m_sName,
				g_waAtomMolNumber[a3]+1 );
			out += temp;
			out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i], 2);
		}
		if ( b3 == -1 )
			out += (")] - vec[");
		else
			out += (")] - plane[");
		temp.Format("%s%d(",
			(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[b1]])->m_sName,
			g_waAtomMolNumber[b1]+1 );
		out += temp;
		out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i], 3);
		temp.Format(")-%s%d(",
			(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[b2]])->m_sName,
			g_waAtomMolNumber[b2]+1 );
		out += temp;
		out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i], 4);
		if ( b3 != -1 )
		{
			temp.Format(")-%s%d(",
				(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[b3]])->m_sName,
				g_waAtomMolNumber[b3]+1 );
			out += temp;
			out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i], 5);
		}
		out += (")]\n");
		out.Dump();
	}
}

void CReFluxObservation::DumpAll()
{
	int i, j;
	int a1, a2, a3, b1, b2, b3;
	int i_temp, j_max, j_temp;
	CxString out, temp;

	for ( i=0; i < ((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize(); i++ )
	{
		j_max = ((CReFluxSubCond*)((CReFluxCond*)m_oaCond[0])->m_oaAggrConds[i])->m_iaDistConds.GetSize();
		for ( j=0; j < j_max; j++ )
		{
			i_temp = ((CReFluxSubCond*)((CReFluxCond*)m_oaCond[0])->m_oaAggrConds[i])->m_iaDistConds[j];
			a1 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaDistConds[i_temp])->m_iRefDist;
			a2 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaDistConds[i_temp])->m_iObsVec;
			if ( j == 0 )
				out.Format("%7d:   ", i+1);
			else
				out = ("         + ");
			temp.Format("Distance: %s%d (",
				(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[a1]])->m_sName,
				g_waAtomMolNumber[a1]+1 );
			out += temp;
			out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaDistConds[i], 0);
			temp.Format(") - %s%d (",
				(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[a2]])->m_sName,
				g_waAtomMolNumber[a2]+1 );
			out += temp;
			out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaDistConds[i], 1);
			temp = (")\n");
			out += temp;
			out.Dump();	
		}
		j_temp = j_max;
		for ( j_max += ((CReFluxSubCond*)((CReFluxCond*)m_oaCond[0])->m_oaAggrConds[i])->m_iaAngleConds.GetSize(); j < j_max; j++ )
		{
			i_temp = ((CReFluxSubCond*)((CReFluxCond*)m_oaCond[0])->m_oaAggrConds[i])->m_iaAngleConds[j-j_temp];
			a1 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp])->m_iRefVec;
			a2 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp])->m_iRefVec2;
			a3 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp])->m_iRefVec3;
			b1 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp])->m_iObsVec;
			b2 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp])->m_iObsVec2;
			b3 = ((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp])->m_iObsVec3;
			if ( j == 0 )
				out.Format("%7d:   ", i+1);
			else
				out = ("         + ");
			out += ("Angle: ");
			if ( a3 == -1 )
				out += ("         vec[");
			else
				out += ("plane[");
			temp.Format("%s%d(",
				(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[a1]])->m_sName,
				g_waAtomMolNumber[a1]+1 );
			out += temp;
			out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp], 0);
			temp.Format(")-%s%d(",
				(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[a2]])->m_sName,
				g_waAtomMolNumber[a2]+1 );
			out += temp;
			out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp], 1);
			if ( a3 != -1 )
			{
				temp.Format(")-%s%d(",
					(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[a3]])->m_sName,
					g_waAtomMolNumber[a3]+1 );
				out += temp;
				out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp], 2);
			}
			if ( b3 == -1 )
				out += (")] - vec[");
			else
				out += (")] - plane[");
			temp.Format("%s%d(",
				(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[b1]])->m_sName,
				g_waAtomMolNumber[b1]+1 );
			out += temp;
			out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp], 3);
			temp.Format(")-%s%d(",
				(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[b2]])->m_sName,
				g_waAtomMolNumber[b2]+1 );
			out += temp;
			out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp], 4);
			if ( b3 != -1 )
			{
				temp.Format(")-%s%d(",
					(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[b3]])->m_sName,
					g_waAtomMolNumber[b3]+1 );
				out += temp;
				out += FindType((CReFluxSingleCond*)((CReFluxCond*)m_oaCond[0])->m_oaAngleConds[i_temp], 5);
			}
			out += (")]\n");
			out.Dump();	
		}
	}
}

void CReFluxObservation::SortDist()
{
	if ( ((CReFluxCond*)m_oaCond[0])->m_oaDistConds.GetSize() == 0 )
		return;

	int i, j;
	CxString buf;
	CxIntArray temp;
	CReFluxSubCond* subcond;

	mprintf("\n    Found the following distance couples:\n");
	DumpDist();
	while ( true )
	{
		mprintf("\n    Which couples should be ignored (comma separated numbers) [none]?");
		inpprintf("! Which couples should be ignored (comma separated numbers) [none]?\n");
		myget(&buf);
		if ( strlen(buf) != 0 ) 
			if ( !GetArray(buf, &temp, ((CReFluxCond*)m_oaCond[0])->m_oaDistConds.GetSize()) )
				continue;
		break;
	}
	if ( temp.GetSize() != 0 )
		for ( i=((CReFluxCond*)m_oaCond[0])->m_oaDistConds.GetSize()-1; i != -1; i-- )
			if ( temp.Contains(i) )
				for ( j=0; j < m_oaCond.GetSize(); j++ )
					((CReFluxCond*)m_oaCond[j])->m_oaDistConds.RemoveAt(i,1);
	mprintf("\n    Left are the following distance couples:\n");
	DumpDist();
	mprintf("\n    It is possible to define groups of criteria which have to be fulfilled simultaneously.\n");
	mprintf("    You have now the chance to define one or more groups.\n");
	while ( true )
	{
		mprintf("\n    Which criteria should be paired (comma separated numbers) [end]?");
		inpprintf("! Which criteria should be paired (comma separated numbers) [end]?\n");
		myget(&buf);
		if ( strlen(buf) == 0 )
			break;
		if ( !GetArray(buf, &temp, ((CReFluxCond*)m_oaCond[0])->m_oaDistConds.GetSize()) )
			continue;
		
		for ( i=0; i < m_oaCond.GetSize(); i++ )
		{
			try { subcond = new CReFluxSubCond(); } catch(...) { subcond = NULL; }
			if (subcond == NULL) NewException((double)sizeof(CReFluxSubCond),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			subcond->m_iaDistConds.CopyFrom(&temp);
			subcond->m_oaDistConds = &((CReFluxCond*)m_oaCond[i])->m_oaDistConds;
			subcond->m_oaAngleConds = &((CReFluxCond*)m_oaCond[i])->m_oaAngleConds;

			((CReFluxCond*)m_oaCond[i])->m_oaAggrConds.Add(subcond);
		}
	}
	for ( j=0; j < ((CReFluxCond*)m_oaCond[0])->m_oaDistConds.GetSize(); j++ )
	{
		for ( i=0; i < ((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize(); i++ )
		{
			if ( ((CReFluxSubCond*)((CReFluxCond*)m_oaCond[0])->m_oaAggrConds[i])->m_iaDistConds.Contains(j) )
			{
				i = -1;
				break;
			}
		}
		if ( i == -1 )
			continue;
		for ( i=0; i < m_oaCond.GetSize(); i++ )
		{
			try { subcond = new CReFluxSubCond(); } catch(...) { subcond = NULL; }
			if (subcond == NULL) NewException((double)sizeof(CReFluxSubCond),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			subcond->m_iaDistConds.Add(j);
			subcond->m_oaDistConds = &((CReFluxCond*)m_oaCond[i])->m_oaDistConds;
			subcond->m_oaAngleConds = &((CReFluxCond*)m_oaCond[i])->m_oaAngleConds;

			((CReFluxCond*)m_oaCond[i])->m_oaAggrConds.Add(subcond);
		}
	}
}

void CReFluxObservation::SortAngle()
{
	int i, j, n, m;
	CxString buf;
	CxIntArray temp, temp2;
	CReFluxSubCond* subcond;

	if ( ((CReFluxCond*)m_oaCond[0])->m_oaAngleConds.GetSize() == 0 )
		return;
	
	mprintf("\n    Found the following angle definitions:\n");
	DumpAngle();
	while ( true )
	{
		mprintf("\n    Which definitions should be ignored (comma separated numbers) [none]?");
		inpprintf("! Which definitions should be ignored (comma separated numbers) [none]?\n");
		myget(&buf);
		if ( strlen(buf) != 0 ) 
			if ( !GetArray(buf, &temp, ((CReFluxCond*)m_oaCond[0])->m_oaAngleConds.GetSize()) )
				continue;
		break;
	}
	if ( temp.GetSize() != 0 )
		for ( i=((CReFluxCond*)m_oaCond[0])->m_oaAngleConds.GetSize()-1; i != -1; i-- )
			if ( temp.Contains(i) )
				for ( j=0; j < m_oaCond.GetSize(); j++ )
					((CReFluxCond*)m_oaCond[j])->m_oaAngleConds.RemoveAt(i,1);
	mprintf("\n    Left are the following angle definitions ...\n");
	DumpAngle();
	if ( ((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize() > 0 )
	{
		mprintf("\n    ... and the following distance criteria:\n");
		DumpAll();
	}
	mprintf("\n    It is possible to define groups of criteria which have to be fulfilled simultaneously.\n");
	mprintf("    You have now the chance to define one or more groups.\n");
	if ( ((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize() > 0 )
	{
		mprintf("    Please define the groups like this: x1,x2-y1,y2,y3 (x=distances, y=angles).\n");
		mprintf("    Or, if a group only contains angles: -y1,y2\n");
	}
	else
		mprintf("    Please define the groups one by one as comma separated numbers.\n");
	mprintf("    Just press enter to end.\n");
	n = ((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize();
	m = ((CReFluxCond*)m_oaCond[0])->m_oaAngleConds.GetSize();
	while ( true )
	{
		inpprintf("! Which criteria should be paired [end]?\n");
		myget(&buf);
		if ( strlen(buf) == 0 )
			break;
		if ( n > 0 )
			if ( !GetDoubleArray( buf, &temp, &temp2, n, m ) )
				continue;
		if ( m == 0 )
			if ( !GetArray( buf, &temp2, m ) )
				continue;
		for ( i=0; i < m_oaCond.GetSize(); i++ )
		{
			try { subcond = new CReFluxSubCond(); } catch(...) { subcond = NULL; }
			if (subcond == NULL) NewException((double)sizeof(CReFluxSubCond),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
			for ( j=0; j < temp.GetSize(); j++ )
				subcond->m_iaDistConds.Append( &((CReFluxSubCond*)((CReFluxCond*)m_oaCond[i])->m_oaAggrConds[temp[j]])->m_iaDistConds );
			subcond->m_iaAngleConds.CopyFrom(&temp2);
			subcond->m_oaDistConds = &((CReFluxCond*)m_oaCond[i])->m_oaDistConds;
			subcond->m_oaAngleConds = &((CReFluxCond*)m_oaCond[i])->m_oaAngleConds;
			((CReFluxCond*)m_oaCond[i])->m_oaAggrConds.Add(subcond);
		}
	}
	for ( j=0; j < n; j++ )
	{
		for ( i=n; i < ((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize(); i++ )
		{
			if ( ((CReFluxSubCond*)((CReFluxCond*)m_oaCond[0])->m_oaAggrConds[i])->m_iaDistConds.Contains(((CReFluxSubCond*)((CReFluxCond*)m_oaCond[0])->m_oaAggrConds[j])->m_iaDistConds[0]) )
			{
				i = -1;
				break;
			}
		}
		if ( i != -1 )
		{
			for ( i=0; i < m_oaCond.GetSize(); i++ )
			{
				try { subcond = new CReFluxSubCond(); } catch(...) { subcond = NULL; }
				if (subcond == NULL) NewException((double)sizeof(CReFluxSubCond),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				subcond->m_iaDistConds.CopyFrom( &((CReFluxSubCond*)((CReFluxCond*)m_oaCond[i])->m_oaAggrConds[j])->m_iaDistConds );
				subcond->m_oaDistConds = &((CReFluxCond*)m_oaCond[i])->m_oaDistConds;
				subcond->m_oaAngleConds = &((CReFluxCond*)m_oaCond[i])->m_oaAngleConds;
				((CReFluxCond*)m_oaCond[i])->m_oaAggrConds.Add(subcond);
			}
		}
	}
	for ( j=0; j < m; j++ )
	{
		for ( i=n; i < ((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize(); i++ )
		{
			if ( ((CReFluxSubCond*)((CReFluxCond*)m_oaCond[0])->m_oaAggrConds[i])->m_iaAngleConds.Contains(j) )
			{
				i = -1;	
				break;
			}
		}
		if ( i != -1 )
		{
			for ( i=0; i < m_oaCond.GetSize(); i++ )
			{
				try { subcond = new CReFluxSubCond(); } catch(...) { subcond = NULL; }
				if (subcond == NULL) NewException((double)sizeof(CReFluxSubCond),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				subcond->m_iaAngleConds.Add(j);
				subcond->m_oaDistConds = &((CReFluxCond*)m_oaCond[i])->m_oaDistConds;
				subcond->m_oaAngleConds = &((CReFluxCond*)m_oaCond[i])->m_oaAngleConds;
				((CReFluxCond*)m_oaCond[i])->m_oaAggrConds.Add(subcond);
			}
		}
	}
	for ( i=0; i < m_oaCond.GetSize(); i++ )
		((CReFluxCond*)m_oaCond[i])->m_oaAggrConds.RemoveAt(0,n);	
}

void CReFluxObservation::SortFinal()
{
	int i, j;
	CxString buf;
	CxIntArray temp;

	if ( ((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize() == 0 )
		return;

	mprintf("\n    Your choice result in the following aggregation criteria:\n");
	DumpAll();

	mprintf("\n    If some of these criteria are incorrect it is now possible to discard them.");
	while ( true )
	{
		mprintf("\n    Which criteria should be discarded (comma separated numbers) [none]?");
		inpprintf("! Which criteria should be discarded (comma separated numbers) [none]?\n");
		myget(&buf);
		if ( strlen(buf) != 0 ) 
			if ( !GetArray(buf, &temp, ((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize()) )
				continue;
		break;
	}
	if ( temp.GetSize() != 0 )
		for ( i=((CReFluxCond*)m_oaCond[0])->m_oaAggrConds.GetSize()-1; i != -1; i-- )
			if ( temp.Contains(i) )
				for ( j=0; j < m_oaCond.GetSize(); j++ )
					((CReFluxCond*)m_oaCond[j])->m_oaAggrConds.RemoveAt(i,1);
}

CReFluxSubCond::CReFluxSubCond()
{
}

bool CReFluxSubCond::check(CTimeStep* ts)
{
	int i, a1, a2, a3;
	CxDVector3 Vec, Vec2, Vec3;
	double cos;

	for ( i=0; i < m_iaDistConds.GetSize(); i++ )
	{
		a1 = ((CReFluxSingleCond*)(*m_oaDistConds)[i])->m_iRefDist;
		a2 = ((CReFluxSingleCond*)(*m_oaDistConds)[i])->m_iObsVec;
		Vec = FoldVector( ts->m_vaCoords[a1] - ts->m_vaCoords[a2] );
		if ( Vec.GetLength() > ((CReFluxSingleCond*)(*m_oaDistConds)[i])->m_dAggrMinAngle )
			return false;
	}
	
	for ( i=0; i < m_iaAngleConds.GetSize(); i++ )
	{
		a1 = ((CReFluxSingleCond*)(*m_oaAngleConds)[i])->m_iRefVec;
		a2 = ((CReFluxSingleCond*)(*m_oaAngleConds)[i])->m_iRefVec2;
		a3 = ((CReFluxSingleCond*)(*m_oaAngleConds)[i])->m_iRefVec3;
		Vec = FoldVector( ts->m_vaCoords[a1] - ts->m_vaCoords[a2] );
		if ( a3 != -1 )
		{
			Vec3 = FoldVector( ts->m_vaCoords[a3] - ts->m_vaCoords[a2] );
			Vec = CrossP(Vec, Vec3);
		}
		Vec.Normalize();
		a1 = ((CReFluxSingleCond*)(*m_oaAngleConds)[i])->m_iObsVec;
		a2 = ((CReFluxSingleCond*)(*m_oaAngleConds)[i])->m_iObsVec2;
		a3 = ((CReFluxSingleCond*)(*m_oaAngleConds)[i])->m_iObsVec3;
		Vec2 = FoldVector( ts->m_vaCoords[a1] - ts->m_vaCoords[a2] );
		if ( a3 != -1 )
		{
			Vec3 = FoldVector( ts->m_vaCoords[a3] - ts->m_vaCoords[a2] );
			Vec2 = CrossP(Vec2, Vec3);
		}
		Vec2.Normalize();
		cos = DotP( Vec, Vec2 );
		if ( ( cos > ((CReFluxSingleCond*)(*m_oaAngleConds)[i])->m_dAggrMinAngle )
			|| ( cos < ((CReFluxSingleCond*)(*m_oaAngleConds)[i])->m_dAggrMaxAngle ) )
			return false;
	}

	return true;
}

CReFluxSubCond::~CReFluxSubCond()
{
}

CReFluxSingleCond::CReFluxSingleCond()
{
	m_iRefDist = -1;		// Ref for Dist
	m_iRefVec  = -1;		// Reference vectors base for angle crit
	m_iRefVec2 = -1;		// Reference vectors tipp
	m_iRefVec3 = -1;		// only used if Vectors are defined by plane
	m_iObsVec  = -1;		// Observed vectors base for angle crit of Obs for Dist
	m_iObsVec2 = -1;		// Observed vectors tipp
	m_iObsVec3 = -1;		// only used if Vectors are defined by plane
	m_dAggrMinAngle = -1.0;		// Angle for AngleCrit of Dist for Dist
	m_dAggrMaxAngle = -1.0;		// Angle for AngleCrit
}

CReFluxSingleCond::~CReFluxSingleCond()
{
}

void CReFluxSingleCond::Copy(CReFluxSingleCond* sc)
{
	m_iRefDist = sc->m_iRefDist;		
	m_iRefVec  = sc->m_iRefVec;		
	m_iRefVec2 = sc->m_iRefVec2;		
	m_iRefVec3 = sc->m_iRefVec3;		
	m_iObsVec  = sc->m_iObsVec;		
	m_iObsVec2 = sc->m_iObsVec2;		
	m_iObsVec3 = sc->m_iObsVec3;		
	m_dAggrMinAngle = sc->m_dAggrMinAngle;		
	m_dAggrMaxAngle = sc->m_dAggrMaxAngle;	
	m_iaTypes.CopyFrom(&sc->m_iaTypes);
}
