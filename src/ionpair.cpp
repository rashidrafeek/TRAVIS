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

#include "ionpair.h"
#include "luzar.h"
#include <math.h>

#include "globalvar.h"
#include "maintools.h"


const char *GetRevisionInfo_ionpair(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_ionpair() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


static CxObArray g_IonPairObserv;


bool gatherIonPair()
{
	CIonPairObservation *o;

	try { o = new CIonPairObservation(); } catch(...) { o = NULL; }
	if (o == NULL) NewException((double)sizeof(CIonPairObservation),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	g_IonPairObserv.Add(o);

	return true;
}


bool initializeIonPair() 
{
    int i;

    for(i = 0; i < g_IonPairObserv.GetSize(); i++) 
	{
        ((CIonPairObservation *)g_IonPairObserv[i])->initialize();
    }
    return true;
}


bool processIonPair(CTimeStep* ts) 
{
    int i;
 
	for(i = 0; i < g_IonPairObserv.GetSize(); i++) 
	{
        ((CIonPairObservation *)g_IonPairObserv[i])->process(ts);
    }
    return true;
}


bool finalizeIonPair() 
{
    int i;

    for(i = 0; i < g_IonPairObserv.GetSize(); i++) 
	{
		((CIonPairObservation *)g_IonPairObserv[i])->finalize();
    }
    return true;
}


CIonPairObservation::CIonPairObservation()
{
	int z,n;
	double d_temp;
	CxString buf, buf2;
	CxObArray* oa_temp;

	try { m_oaAnion = new CxObArray(); } catch(...) { m_oaAnion = NULL; }
        if (m_oaAnion == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { m_oaCation = new CxObArray(); } catch(...) { m_oaCation = NULL; }
        if (m_oaCation == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	mprintf("\n    First TRAVIS will gather which molecules are to be handled as cations and anions. \n");
	if (g_oaMolecules.GetSize() > 1)
	{
       		buf.sprintf("\n    Which of the molecules are cations (");
        	for (z=0;z<g_oaMolecules.GetSize();z++)
        	{
        	    buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
        	    buf.strcat(buf2);
        	    if (z < g_oaMolecules.GetSize()-1)
        	        buf.strcat(", ");
        	}
        	buf.strcat(")? [enter if finished]");
		z = -1;
		while ( z != 0 )
		{
        		z = AskRangeInteger("%s",0,g_oaMolecules.GetSize(),0,(const char*)buf);
			if ( ia_cation.Contains(z-1) )
			{
				mprintf(WHITE,"\n    %s is already chosen as cation.\n\n",((CMolecule*)g_oaMolecules[z-1])->m_sName);
			}
			else if ( z != 0 )
			{
				ia_cation.Add(z-1);
				mprintf(WHITE,"\n    Adding %s to the cations.\n\n",((CMolecule*)g_oaMolecules[z-1])->m_sName);
			}
			else if ( ia_cation.GetSize() == 0 )
			{
				mprintf(WHITE,"\n    Need at least one cation for analysis!\n\n");
				z = -1;
			}
		}

        	buf.sprintf("\n    Which of the molecules are anions (");
        	for (z=0;z<g_oaMolecules.GetSize();z++)
        	{
        	    buf2.sprintf("%s=%d",((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
        	    buf.strcat(buf2);
        	    if (z < g_oaMolecules.GetSize()-1)
        	        buf.strcat(", ");
        	}
        	buf.strcat(")? [enter if finished]");
		z = -1;
		while ( z != 0 )
		{
        		z = AskRangeInteger("%s",0,g_oaMolecules.GetSize(),0,(const char*)buf);
			if ( ia_cation.Contains(z-1) )
			{
				mprintf(WHITE,"\n    %s is already chosen as cation.\n",((CMolecule*)g_oaMolecules[z-1])->m_sName);
				if ( !AskYesNo("    Choose as anion as well (This may cause some weird effects!) ? [yes]", true) )
				{
					continue;
				}
			}
			if ( ia_anion.Contains(z-1) )
			{
				mprintf(WHITE,"\n    %s is already chosen as anion.\n\n",((CMolecule*)g_oaMolecules[z-1])->m_sName);
			}
			else if ( z != 0 )
			{
				ia_anion.Add(z-1);
				mprintf(WHITE,"\n    Adding %s to the anions.\n\n",((CMolecule*)g_oaMolecules[z-1])->m_sName);
			}
			else if ( ia_anion.GetSize() == 0 )
			{
				mprintf(WHITE,"\n    Need at least one anion for analysis!\n\n");
				z = -1;
			}
		}
	} 
	else 
	{
		mprintf("\n    There is only one molecule, which is handled as cation AND as anion. \n");
	 	ia_anion.Add(0);
		ia_cation.Add(0);		
	}

	mprintf("\n    Now TRAVIS needs to know which (virtual) atom is the center of the respective ions. ");
	mprintf("\n    If the analysis is performed with more than one cation/anion it may be necessary to ");
	mprintf("\n    correct the size of the smaller ions.\n");
	for ( z=0; z < ia_cation.GetSize(); z++ )
	{
		try { oa_temp = new CxObArray(); } catch(...) { oa_temp = NULL; }
        	if (oa_temp == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		while ( oa_temp->GetSize() == 0 )
		{
	        	mprintf("    Which atom from %s is the cationic center? [#2] ",((CMolecule*)g_oaMolecules[ia_cation[z]])->m_sName);
	        	inpprintf("! Which atom from %s is the cationic center? [#2] \n",((CMolecule*)g_oaMolecules[ia_cation[z]])->m_sName);
	        	myget(&buf);
			if (buf.GetLength() == 0)
			{
				buf.sprintf("#2");
			}
        		checkAtomChoice((CMolecule*)g_oaMolecules[ia_cation[z]],buf,oa_temp);
		}

		if ( ia_cation.GetSize() > 1 )
		{
			d_temp = AskFloat("    Please enter size correction in pm! [0.0]",0.0 );
			
			if ( d_temp != 0.0 )
			{
				for ( n=0; n < oa_temp->GetSize(); n++ )
				{
					((CIon*)(*oa_temp)[n])->m_dCorr = d_temp;
				}
			}
		}

		m_oaCation->Add(oa_temp);		
	}

	for ( z=0; z < ia_anion.GetSize(); z++ )
	{
		try { oa_temp = new CxObArray(); } catch(...) { oa_temp = NULL; }
        	if (oa_temp == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		while ( oa_temp->GetSize() == 0 )
		{
	        	mprintf("    Which atom from %s is the anionic center? [#2] ",((CMolecule*)g_oaMolecules[ia_anion[z]])->m_sName);
	        	inpprintf("! Which atom from %s is the anionic center? [#2] \n",((CMolecule*)g_oaMolecules[ia_anion[z]])->m_sName);
	        	myget(&buf);
			if (buf.GetLength() == 0)
			{
				buf.sprintf("#2");
			}
        		checkAtomChoice((CMolecule*)g_oaMolecules[ia_anion[z]],buf,oa_temp); 
		}

		if ( ia_anion.GetSize() > 1 )
		{
			d_temp = AskFloat("    Please enter size correction in pm! [0.0]",0.0 );
			
			if ( d_temp != 0.0 )
			{
				for ( n=0; n < oa_temp->GetSize(); n++ )
				{
					((CIon*)oa_temp->GetAt(n))->m_dCorr = d_temp;
				}
			}
		}
		
		m_oaAnion->Add(oa_temp);		
	}

	m_fDist = AskFloat("    Maximum distance between cation and anion in pm (for this distance no corrections are applied!)? [350] ",350.0);
	g_fTimestepLength = AskFloat("    Enter the physical time distance between successive trajectory frames in fs: [0.5] ",0.5f);

}


CIonPairObservation::~CIonPairObservation()
{
}


void CIonPairObservation::initialize()
{
	int n_c, n_a, m_c, m_a;

	try { m_oaLastNeigh = new CxObArray(); } catch(...) { m_oaLastNeigh = NULL; }
	if (m_oaLastNeigh == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_laLastPair = new CxIntArray(); } catch(...) { m_laLastPair = NULL; }
	if (m_laLastPair == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_oaNeigh = new CxObArray(); } catch(...) { m_oaNeigh = NULL; }
	if (m_oaNeigh == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_oaPair = new CxObArray(); } catch(...) { m_oaPair = NULL; }
	if (m_oaPair == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	if ( g_iMaxStep != -1 )
	{
	    m_iDepth = g_iMaxStep * 3 / 4 / g_iStride;
	}
	else
	{
	    m_iDepth = g_iTrajSteps * 3 / 4 / g_iStride;
	}
	m_iDepth = AskUnsignedInteger("    Enter the correlation depth of the ACF (in trajectory frames): [%d] ", m_iDepth, m_iDepth);
	if ( m_iDepth > g_iTrajSteps )
	{
	    mprintf("    Correlation Depth is exceeding number of total timesteps. Truncating ....\n");
	    m_iDepth = g_iTrajSteps;
	}

	m_iTotalCats = 0;
	m_iTotalAns = 0;	
			
	for ( n_c=0; n_c < m_oaCation->GetSize(); n_c++)		// Walks over all Cation Types
	{
		for ( m_c=0; m_c < ((CxObArray*)(*m_oaCation)[n_c])->GetSize(); m_c++ )	// Walks over all cations of type
		{
			m_iTotalCats++;
			m_laLastPair->Add(-1);	

	        	CxIntArray* IncludeME3;

			try { IncludeME3 = new CxIntArray(); } catch(...) { IncludeME3 = NULL; }
			if (IncludeME3 == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			for ( n_a=0; n_a < m_oaAnion->GetSize(); n_a++)		// Walks over all Anion Types
			{
				for ( m_a=0; m_a < ((CxObArray*)(*m_oaAnion)[n_a])->GetSize(); m_a++ )	// Walks over all anions of type
				{
					if ( m_iTotalCats == 1 )
					{
						m_iTotalAns++;
					}

					CxIntArray* IncludeME;
            				CxIntArray* IncludeME2;

					try { IncludeME = new CxIntArray(); } catch(...) { IncludeME = NULL; }
					if (IncludeME == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					try { IncludeME2 = new CxIntArray(); } catch(...) { IncludeME2 = NULL; }
					if (IncludeME2 == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					m_oaNeigh->Add(IncludeME);
					m_oaPair->Add(IncludeME2);
		
					IncludeME3->Add(0);
				}
			}
			m_oaLastNeigh->Add(IncludeME3);
		}
	}

	for ( n_a=0; n_a < m_oaAnion->GetSize(); n_a++)		// Walks over all Anion Types
	{
		for ( m_a=0; m_a < ((CxObArray*)(*m_oaAnion)[n_a])->GetSize(); m_a++ )	// Walks over all anions of type
		{
			m_laLastPair->Add(-1);	

	        	CxIntArray* IncludeME3;

			try { IncludeME3 = new CxIntArray(); } catch(...) { IncludeME3 = NULL; }
			if (IncludeME3 == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			for ( n_c=0; n_c < m_oaCation->GetSize(); n_c++)		// Walks over all Cation Types
			{
				for ( m_c=0; m_c < ((CxObArray*)(*m_oaCation)[n_c])->GetSize(); m_c++ )	// Walks over all cations of type
				{
					CxIntArray* IncludeME;
            				CxIntArray* IncludeME2;

					try { IncludeME = new CxIntArray(); } catch(...) { IncludeME = NULL; }
					if (IncludeME == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					try { IncludeME2 = new CxIntArray(); } catch(...) { IncludeME2 = NULL; }
					if (IncludeME2 == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

					m_oaNeigh->Add(IncludeME);
					m_oaPair->Add(IncludeME2);

					IncludeME3->Add(0);
				}
			}
			m_oaLastNeigh->Add(IncludeME3);
		}
	}

	mprintf("    Found in total %d cations and %d anions!\n",m_iTotalCats, m_iTotalAns);

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


void CIonPairObservation::process(CTimeStep* ts)
{
	int Cur_Pair = -1;
	double Min_Pair;
	int n_c, n_a, m_c, m_a;
	int m, n = -1;
	CxDVector3 Vec;
	CIon *cat, *an;

	for ( n_c=0; n_c < m_oaCation->GetSize(); n_c++)		// Walks over all Cation Types
	{
		for ( m_c=0; m_c < ((CxObArray*)(*m_oaCation)[n_c])->GetSize(); m_c++ )	// Walks over all cations of type
		{
			n++;
			cat = ((CIon*)(*(CxObArray*)(*m_oaCation)[n_c])[m_c]);

			Min_Pair = 2 * g_fBoxX;
			m = -1;

			for ( n_a=0; n_a < m_oaAnion->GetSize(); n_a++)		// Walks over all Anion Types
			{
				for ( m_a=0; m_a < ((CxObArray*)(*m_oaAnion)[n_a])->GetSize(); m_a++ )	// Walks over all anions of type
				{
					m++;
					an = ((CIon*)(*(CxObArray*)(*m_oaAnion)[n_a])[m_a]);

					if ( cat->m_iOffset == an->m_iOffset )
					{
						continue;
					}
			
					Vec = ts->m_vaCoords[cat->m_iOffset] - ts->m_vaCoords[an->m_iOffset];
					Vec = FoldVector(Vec);
					
					if ( Vec.GetLength() <= m_fDist ) 	// In cutoff range -> Neighbor
					{
						if ( Vec.GetLength() < Min_Pair )
						{
							Min_Pair = Vec.GetLength();
							Cur_Pair = m;
						}

						if ( (*(CxIntArray*)(*m_oaLastNeigh)[n])[m] == 0 ) 	// Was not a neighbor ... hmmmm ....
						{
							((CxIntArray*)(*m_oaNeigh)[n*m_iTotalAns+m])->Add(g_iSteps);
							(*(CxIntArray*)(*m_oaLastNeigh)[n])[m] = 1;
						}
					}
					if ( Vec.GetLength() > m_fDist ) 	// Out of cutoff range -> No Neighbor
					{
						if ( (*(CxIntArray*)(*m_oaLastNeigh)[n])[m] == 1 ) 	// Was a neighbor ... hmmmm ....
						{
							((CxIntArray*)(*m_oaNeigh)[n*m_iTotalAns+m])->Add(g_iSteps);
							(*(CxIntArray*)(*m_oaLastNeigh)[n])[m] = 0;
						}
					}
				} // m_a end
			} // n_a end	

			if ( Cur_Pair == -1 )
			{
				eprintf("\n     Failure! There was no next counterion found for at least one of the cations.\n");
				eprintf("     This is mostly caused by a to large maximum distance.\n");
				abort();
			}

			if ( Cur_Pair != (*m_laLastPair)[n] )
			{
				if ( (*m_laLastPair)[n] != -1 )
				{
					((CxIntArray*)(*m_oaPair)[ ( n * m_iTotalAns ) + (*m_laLastPair)[n] ])->Add(g_iSteps);
				}
				((CxIntArray*)(*m_oaPair)[ ( n * m_iTotalAns ) + Cur_Pair ])->Add(g_iSteps);
				(*m_laLastPair)[n] = Cur_Pair;
			}
		} // m_c end
	} // n_c end

	m = -1;
	for ( n_a=0; n_a < m_oaAnion->GetSize(); n_a++)		// Walks over all Anion Types
	{
		for ( m_a=0; m_a < ((CxObArray*)(*m_oaAnion)[n_a])->GetSize(); m_a++ )	// Walks over all anions of type
		{
			m++;
			an = ((CIon*)(*(CxObArray*)(*m_oaAnion)[n_a])[m_a]);

			Min_Pair = 2 * g_fBoxX;

			n = -1;
			for ( n_c=0; n_c < m_oaCation->GetSize(); n_c++)		// Walks over all Cation Types
			{
				for ( m_c=0; m_c < ((CxObArray*)(*m_oaCation)[n_c])->GetSize(); m_c++ )	// Walks over all cations of type
				{
					n++;
					cat = ((CIon*)(*(CxObArray*)(*m_oaCation)[n_c])[m_c]);

					if ( cat->m_iOffset == an->m_iOffset )
					{
						continue;
					}
			
					Vec = ts->m_vaCoords[cat->m_iOffset] - ts->m_vaCoords[an->m_iOffset];
					Vec = FoldVector(Vec);

					if ( Vec.GetLength() <= m_fDist ) 	// In cutoff range -> Neighbor
					{
						if ( Vec.GetLength() < Min_Pair )
						{
							Min_Pair = Vec.GetLength();
							Cur_Pair = n;
						}

						if ( (*(CxIntArray*)(*m_oaLastNeigh)[m+m_iTotalCats])[n] == 0 ) 	// Was not a neighbor ... hmmmm ....
						{
							((CxIntArray*)(*m_oaNeigh)[m_iTotalCats*m_iTotalAns+m*m_iTotalCats+n])->Add(g_iSteps);
							(*(CxIntArray*)(*m_oaLastNeigh)[m+m_iTotalCats])[n] = 1;
						}
					}
					if ( Vec.GetLength() > m_fDist ) 	// Out of cutoff range -> No Neighbor
					{
						if ( (*(CxIntArray*)(*m_oaLastNeigh)[m+m_iTotalCats])[n] == 1 ) 	// Was a neighbor ... hmmmm ....
						{
							((CxIntArray*)(*m_oaNeigh)[m_iTotalCats*m_iTotalAns+m*m_iTotalCats+n])->Add(g_iSteps);
							(*(CxIntArray*)(*m_oaLastNeigh)[m+m_iTotalCats])[n] = 0;
						}
					}
				} // m_a end
			} // n_a end	
		
			if ( Cur_Pair == -1 )
			{
				eprintf("\n     Failure! There was no next counterion found for at least one of the anions.\n");
				eprintf("     This is mostly caused by a to small maximum distance.\n");
				abort();
			}

			if ( Cur_Pair != (*m_laLastPair)[m+m_iTotalCats] )
			{
				if ( (*m_laLastPair)[m+m_iTotalCats] != -1 )
				{
					((CxIntArray*)(*m_oaPair)[ m_iTotalCats * m_iTotalAns + m * m_iTotalCats + (*m_laLastPair)[m+m_iTotalCats] ])->Add(g_iSteps);
				}
				((CxIntArray*)(*m_oaPair)[ m_iTotalCats * m_iTotalAns + m * m_iTotalCats + Cur_Pair ])->Add(g_iSteps);
				(*m_laLastPair)[m+m_iTotalCats] = Cur_Pair;
			}
		} // m_a end
	} // n_a end
}

void CIonPairObservation::finalize()	// reconstruct functions -> autocorrelate -> sum results -> normalize -> fit
{
	long n,m,i,interval, count, z;
	double PairCor, DiffCor;
	CxString buf,buf2;
	unsigned long t0,eta;

	CxDoubleArray* NeighFunction;	// Input
	CxDoubleArray* PairFunction;
	
	CxDoubleArray* temp1Function;	// Temp
	CxDoubleArray* temp2Function;
	CxDoubleArray* temp3Function;

	CLuzarCorrelation* luzar;

	try { PairFunction = new CxDoubleArray(); } catch(...) { PairFunction = NULL; }
	if (PairFunction == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { NeighFunction = new CxDoubleArray(); } catch(...) { NeighFunction = NULL; }
	if (NeighFunction == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { PairDynamic = new CxObArray(); } catch(...) { PairDynamic = NULL; }
	if (PairDynamic == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { DiffusionDynamic = new CxObArray(); } catch(...) { DiffusionDynamic = NULL; }
	if (DiffusionDynamic == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
 
	try { k_t = new CxObArray(); } catch(...) { k_t = NULL; }
	if (k_t == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { a = new CxObArray(); } catch(...) { a = NULL; }
	if (a == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { b = new CxObArray(); } catch(...) { b = NULL; }
	if (b == NULL) NewException((double)sizeof(CxObArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for ( n=0; n < (m_oaCation->GetSize()+m_oaAnion->GetSize()); n++)                // Walks over all Cation Types
	{
		try { temp1Function = new CxDoubleArray(); } catch(...) { temp1Function = NULL; }
		if (temp1Function == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { temp2Function = new CxDoubleArray(); } catch(...) { temp2Function = NULL; }
		if (temp2Function == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		for ( i=0; i < m_iDepth; i++ )
		{
			temp1Function->Add(0.0);
			temp2Function->Add(0.0);
		}
	
		PairDynamic->Add(temp1Function);
		DiffusionDynamic->Add(temp2Function);
	
		try { temp1Function = new CxDoubleArray(); } catch(...) { temp1Function = NULL; }
		if (temp1Function == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { temp2Function = new CxDoubleArray(); } catch(...) { temp2Function = NULL; }
		if (temp2Function == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		try { temp3Function = new CxDoubleArray(); } catch(...) { temp3Function = NULL; }
		if (temp3Function == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		k_t->Add(temp1Function);
		a->Add(temp2Function);
		b->Add(temp3Function);

	}

	try { temp1Function = new CxDoubleArray(); } catch(...) { temp1Function = NULL; }
	if (temp1Function == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { temp2Function = new CxDoubleArray(); } catch(...) { temp2Function = NULL; }
	if (temp2Function == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { luzar = new CLuzarCorrelation(); } catch(...) { luzar = NULL; }
	if (luzar == NULL) NewException((double)sizeof(CLuzarCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	NeighFunction->SetSize(g_iSteps / g_iStride);
	PairFunction->SetSize(g_iSteps / g_iStride);

	for ( i=0; i < (int)(g_iSteps/g_iStride); i++ ) 
	{
		(*NeighFunction)[i] = 0;
		(*PairFunction)[i] = 0;
	}

	luzar->Init(g_iSteps / g_iStride,m_iDepth);		// Inititalize Luzar-Correlation

	mprintf("    Finishing %6i Functions ... \nFunc %6li ",m_oaNeigh->GetSize(),(long)0);
	g_iDotCounter = 0;
	t0 = (unsigned long)time(NULL);
	interval = 2;
	count = -1;

	for ( n=0; n < (m_oaCation->GetSize() + m_oaAnion->GetSize()); n++)		// Walks over all ion Types
	{
		if ( n < m_oaCation->GetSize() )
			z = (((CxObArray*)(*m_oaCation)[n])->GetSize() * m_iTotalAns) ;
		else
			z = (((CxObArray*)(*m_oaAnion)[n-m_oaCation->GetSize()])->GetSize() * m_iTotalCats);
			
		for ( m=0; m < z; m++ )	// Walks over all ions of type
		{
			count++;
			if ( count % interval == 0 )
			{
				if ( g_iDotCounter >= 50 )
				{
					if ( interval == 2 && count > 100 )
					{
						interval = 20;
					}
					else if ( interval == 20 && count > 2000 )
					{
						interval = 200;
					}
				    g_iDotCounter = 0;
			        if ((time(NULL) - t0) > 5)
			        {
			                eta = (unsigned long)(((double)time(NULL) - t0) / count * ((double)MAX((long)0,(long)m_oaNeigh->GetSize() - (long)count)));
			                FormatTime(eta,&buf);
			                mprintf(" ETA %s",(const char*)buf);
			        }
			        mprintf("\nFunc %6li ",count);
				}
				g_iDotCounter++;
				mprintf(".");
			}

			if ( ((CxIntArray*)(*m_oaPair)[count])->GetSize() > 0 )
			{
				RestoreFunction(((CxIntArray*)(*m_oaNeigh)[count]),((CxDoubleArray*)NeighFunction));	// Restore function of pair n as NeighFunction and BondFunction
				RestoreFunction(((CxIntArray*)(*m_oaPair)[count]),((CxDoubleArray*)PairFunction));	

				luzar->Correlate(PairFunction,NeighFunction,temp1Function,temp2Function);
				
	        		for ( i=0; i < m_iDepth; i++ )		// Sum results to Output-array
	        		{
	        			(*(CxDoubleArray*)(*PairDynamic)[n])[i] += (*temp1Function)[i];
	        			(*(CxDoubleArray*)(*DiffusionDynamic)[n])[i] += (*temp2Function)[i];
	        		}
			}
		}	// end of individuals

		if ( (*(CxDoubleArray*)(*PairDynamic)[n])[0] == 0 || (*(CxDoubleArray*)(*DiffusionDynamic)[n])[0] == 0 )
		{
			eprintf("\n    First correlation value is zero ... something is wrong ...\n");
			abort();
		}

		PairCor = (double)(*(CxDoubleArray*)(*PairDynamic)[n])[0];
		DiffCor = (double)(*(CxDoubleArray*)(*DiffusionDynamic)[n])[0];

		for ( i=0; i < m_iDepth; i++ )			// Normalize
		{
			(*(CxDoubleArray*)(*PairDynamic)[n])[i] /= PairCor;
			(*(CxDoubleArray*)(*DiffusionDynamic)[n])[i] /= DiffCor;
			(*(CxDoubleArray*)(*DiffusionDynamic)[n])[i] -= (*(CxDoubleArray*)(*PairDynamic)[n])[i];
		}
	}	// end of types
	mprintf("\n\n");

	for ( n=0; n < (m_oaCation->GetSize() + m_oaAnion->GetSize()); n++)		// Walks over all Cation Types
	{
		if ( n < m_oaCation->GetSize() )
			z = ia_cation[n];
		else
			z = ia_anion[n-m_oaCation->GetSize()];
			
		mprintf("Fitting %s ...\n",((CMolecule*)g_oaMolecules[z])->m_sName);
		luzar->Fit(((CxDoubleArray*)(*PairDynamic)[n]), ((CxDoubleArray*)(*DiffusionDynamic)[n]), ((CxDoubleArray*)(*k_t)[n]), ((CxDoubleArray*)(*a)[n]), ((CxDoubleArray*)(*b)[n]), m_iTrunc);
		WriteOutput(((CxDoubleArray*)(*PairDynamic)[n]), ((CxDoubleArray*)(*DiffusionDynamic)[n]), ((CxDoubleArray*)(*a)[n]), ((CxDoubleArray*)(*b)[n]), ((CxDoubleArray*)(*k_t)[n]), ((CMolecule*)g_oaMolecules[z])->m_sName);
	}

}


void CIonPairObservation::RestoreFunction(CxIntArray* input, CxDoubleArray* output)
{
	unsigned int n;
	bool eins = false;

        for ( n=0; n < (g_iSteps / g_iStride); n++ )            // Walk over all processed steps
	{
                if ( input->Contains((n+1)*g_iStride) )
			eins = !eins;
		if (!eins)
			output->SetAt(n,0.0);
		else
			output->SetAt(n,1.0);
	}
}


bool CIonPairObservation::includeAtom(int atom, int number, int mol, CxObArray* target)  // atom: element number, number: atomsort number, mol m_iShowMol, target: in which array to include
{
	int n;
	CIon* ion;

	for ( n=0; n < ((CMolecule*)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize(); n++)
	{
		try { ion = new CIon(); } catch(...) { ion = NULL; }
		if (ion == NULL) NewException((double)sizeof(CIon),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		ion->m_iOffset = ((CxIntArray*)(((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[mol])->m_laSingleMolIndex[n]])->m_oaAtomOffset[atom]))->GetAt(number);

		target->Add(ion);
	}

	return true;
}

bool CIonPairObservation::checkAtomChoice(CMolecule *mol, const char *s, CxObArray* m_oaInclude)
{
	const char *p, *q;
	char buf[32];
	const char *allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890,# ";
	int i = -1, atom = -1;	

	p = s;

       	while (*p == ' ')
        	p++;
	if (strchr(allowed,*p) == NULL)
	{
		eprintf("Error: Character \"%c\" not allowed.\n",*p);
		return false;
	}
	q = p;
	if (isalpha(*q) || (*q == '#'))
	{
		while (isalpha(*q) || (*q == '#'))
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
	p = q;
	if (isdigit(*q))
	{
		if (atom == -1)
		{
		    eprintf("Error: Number in the beginning not possible.\n");
		    return false;
		}
		while (isdigit(*q)) 
			q++;
		if (*q != 0) 
		{
			eprintf("Error: Please type only one atom.\n");
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
	}
	i = atoi(buf)-1;

	if ( i == -1 )	
	{
		eprintf("Error: Please type only one atom.\n");
		return false;
	}
	else 	
		includeAtom(atom, i, mol->m_iIndex, m_oaInclude); // hier
	
	return true;
}


void CIonPairObservation::WriteOutput(CxDoubleArray* pairs, CxDoubleArray* diffs, CxDoubleArray* as, CxDoubleArray* bs, CxDoubleArray* k_ts, std::string name)
{
	CGrace* gc;
	int n,m;
	char buf[256];

	mprintf("    Writing \"%s_Correlations.agr\"...\n", name.c_str());
	
	try { gc = new CGrace(); } catch(...) { gc = NULL; }
	if (gc == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	gc->SetRangeX(0,m_iDepth*g_fTimestepLength/1000*g_iStride);
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
	for ( n = 0; n < pairs->GetSize(); n++)
	{
		gc->AddXYTupel(0,n*g_fTimestepLength/1000*g_iStride,(*pairs)[n]);
		gc->AddXYTupel(1,n*g_fTimestepLength/1000*g_iStride,(*diffs)[n]);
	}

	sprintf(buf,"%s_Correlations.agr",name.c_str());
	gc->WriteAgr(buf,false);
	sprintf(buf,"%s_Correlations.csv",name.c_str());
	gc->WriteCSV(buf);
	delete gc;

	mprintf("    Writing \"%s_Fitted.agr\"...\n", name.c_str());

	CalcMinMax(pairs, diffs, as, bs, k_ts);

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

	for ( m=0; m < k_ts->GetSize(); m++)
		gc->AddXYTupel(0,(*as)[0] * (*pairs)[m] - (*bs)[0] * (*diffs)[m], -(*k_ts)[m]);

	sprintf(buf,"%s_Fitted.csv",name.c_str());
	gc->WriteCSV(buf);
	gc->AddDataset();
	gc->SetDatasetName(1, "Unity");
	gc->AddXYTupel(1,m_dGraceMin,m_dGraceMin);
	gc->AddXYTupel(1,m_dGraceMax,m_dGraceMax);
	sprintf(buf,"%s_Fitted.agr",name.c_str());
	gc->WriteAgr(buf,false);
	delete gc;
}

void CIonPairObservation::CalcMinMax(CxDoubleArray* pairs, CxDoubleArray* diffs, CxDoubleArray* as, CxDoubleArray* bs, CxDoubleArray* k_ts)
{
	double temp;
	int n;
	m_dGraceMin =  1E10;
	m_dGraceMax = -1E10;

	for ( n = m_iTrunc; n < PairDynamic->GetSize(); n++)
	{
		temp = ((*as)[0] * (*pairs)[n] - (*bs)[0] * (*diffs)[n]);
		if ( m_dGraceMin >= temp )
			m_dGraceMin = temp;
		if ( m_dGraceMax <= temp )
		    m_dGraceMax = temp;
	}

	for ( n = m_iTrunc; n < k_ts->GetSize(); n++)
	{
		temp = -k_ts->GetAt(n);
        	if ( m_dGraceMin >= temp )
        	    m_dGraceMin = temp;
        	if ( m_dGraceMax <= temp )
        	    m_dGraceMax = temp;
	} 
}

void CIonPairObservation::Initiate_Offsets(int mol, CxObArray* array)
{
	int n;

	for ( n=0; n<((CMolecule*)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize(); n++ )
	{
		CxIntArray* temp;

		try { temp = new CxIntArray(); } catch(...) { temp = NULL; }
        	if (temp == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		array->Add(temp);
	}
}

CIon::CIon()
{
	m_iOffset = -1;
	m_dCorr   = -1.0;
}

CIon::~CIon()
{
}

