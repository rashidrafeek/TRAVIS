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

#include <math.h>
#include "void.h"

#include "globalvar.h"
#include "maintools.h"


const char *GetRevisionInfo_void(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_void() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


using namespace std;


#define CUT_OFF 1E-5
#define SHIFT_START 1.0

static CxObArray g_VoidObserv;

bool gatherVoid()
{
	CVoidObservation *o;

	try { o = new CVoidObservation(); } catch(...) { o = NULL; }
	if (o == NULL) NewException((double)sizeof(CVoidObservation),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	g_VoidObserv.Add(o);

	return true;
}


bool initializeVoid() 
{
	int i;
	
	for(i = 0; i < g_VoidObserv.GetSize(); i++) 
		((CVoidObservation *)g_VoidObserv[i])->initialize();

	return true;
}

bool processVoid(CTimeStep* ts) 
{
	int i;
	
	for(i = 0; i < g_VoidObserv.GetSize(); i++) 
		((CVoidObservation *)g_VoidObserv[i])->process(ts);

	return true;
}

bool finalizeVoid() 
{
	int i;
	
	for(i = 0; i < g_VoidObserv.GetSize(); i++) 
	{
		((CVoidObservation *)g_VoidObserv[i])->finalize();
		delete g_VoidObserv[i];
	}

	return true;
}

CVoidObservation::CVoidObservation()
{
//	time_tot = omp_get_wtime();
	time_dom   = 0.0;
	time_voro1 = 0.0;
	time_voro2 = 0.0;
	time_grow  = 0.0;
	time_crit  = 0.0;

	char temp[20];
	int i, j;
	int bin1 = 0, bin2 = 0;
	double d1 = 0.0, d2 = 0.0, save_d1;
	CxString buf, buf2;
	CxIntArray* ia;
	double cur_rad;
	bool found;
	CDF  *m_pVoidSphDF;
	CDF  *m_pVoidDomDF;
	CDF  *m_pVoidAVDF;
	C2DF *m_pVoidAV2DF;
	CDF  *m_pHoleDF;


/*	for ( i=0; i < g_oaAtoms.GetSize(); i++ )
		if ( ((CAtom*)g_oaAtoms[i])->m_pElement->m_iOrd == 1 )
			m_bNoH = AskYesNo("    Exclude hydrogen atoms (y/n)? [yes] ",true);
*/
	m_bSurfTens = false;
	m_bViscosity = false;
	m_bNoH = false;
	m_bFull = true;
	m_bHoleDF = false;


	if ( m_bFull )
	{
		try { ia = new CxIntArray(); } catch(...) { ia = NULL; }
		if ( ia == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		for ( i=0; i < g_iGesAtomCount; i++ )
		{
			found = false;
			for ( j=0; j < (int)m_oaDomains.size(); j++ )
			{
				if ( m_oaDomains[j]->Contains(i) )
				{
					found = true;
					break;
				}
			}
			if ( !found )
				ia->Add(i);
		}
		m_oaDomains.push_back(ia);
	}	

	m_dMinRadius = AskFloat("    Please enter minimum spherical void radius to consider (in pm): [10]", 10.0);
	m_dCutOff = AskFloat("    Please enter the minimum distance between two void sphere centers (in pm): [%d]",m_dMinRadius, int(m_dMinRadius));

	m_bVoidSphDF = AskYesNo("    Do you want to create a spherical void distribution function (y/n)? [yes] ",true);

	if (m_bVoidSphDF)
	{
		d1 = AskFloat("    Enter the maximal void radius (in pm): [250.0] ",250.0);
		bin1 = AskUnsignedInteger("    Enter the binning resolution: [%d] ", 100, 100);
	}
	save_d1 = d1;
	m_dAV = AskFloat("    Which isoperimetric factor is needed to detect a void domain as spheric hole? [0.5]",0.5);
	d2 = AskFloat("    Enter the maximal representative hole radius (in pm): [500.0] ",500.0);
	bin2 = AskUnsignedInteger("    Enter the binning resolution of the representative radii distribution: [%d] ", 50, 50);

	for ( i=0; i < (int)m_oaDomains.size()+1; i++ )
	{
		try { m_pVoidSphDF = new CDF(); } catch(...) { m_pVoidSphDF = NULL; }
		if (m_pVoidSphDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		if (m_bVoidSphDF)
		{
			m_pVoidSphDF->m_fMinVal = m_dMinRadius;
			m_pVoidSphDF->m_fMaxVal = d1;
			m_pVoidSphDF->m_iResolution = bin1;
			m_pVoidSphDF->Create();
			m_pVoidSphDF->SetLabelX("Void radius / pm");
			m_pVoidSphDF->SetLabelY("Occurrence");
			mprintf("\n");
		}
		else
			m_pVoidSphDF->m_fMaxVal = 250;
		m_oaVoidSphDF.push_back(m_pVoidSphDF);

		try { m_pHoleDF = new CDF(); } catch(...) { m_pHoleDF = NULL; }
		if (m_pHoleDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		m_pHoleDF->m_fMinVal = m_dMinRadius;
		m_pHoleDF->m_fMaxVal = d2;
		m_pHoleDF->m_iResolution = bin2;
		m_pHoleDF->Create();
		m_pHoleDF->SetLabelX("Representative Hole radius / pm");
		m_pHoleDF->SetLabelY("Occurrence");
		mprintf("\n");

		m_oaHoleDF.push_back(m_pHoleDF);
		if ( m_oaDomains.size() == 1 )
			break;
	}


	if ( !m_bHoleDF )
	{
		m_bVoidDomDF = AskYesNo("    Do you want to create a void domain distribution function (y/n)? [yes] ",true);

		if (m_bVoidDomDF)
		{
			d1 = AskFloat("    Enter the maximal void volume (in A**3): [%.1f] ",0.5*pow(g_fBoxX/100,3),0.5*pow(g_fBoxX/100,3));
			bin1 = AskUnsignedInteger("    Enter the binning resolution for the volume distribution: [%d] ",100,100);
			bin2 = AskUnsignedInteger("    Enter the binning resolution of the AV distribution: [%d] ",100,100);

			for ( i=0; i < (int)m_oaDomains.size()+1; i++ )
			{
				try { m_pVoidAVDF = new CDF(); } catch(...) { m_pVoidAVDF = NULL; }
				if (m_pVoidAVDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				try { m_pVoidDomDF = new CDF(); } catch(...) { m_pVoidDomDF = NULL; }
				if (m_pVoidDomDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				try { m_pVoidAV2DF = new C2DF(); } catch(...) { m_pVoidAV2DF = NULL; }
				if (m_pVoidAV2DF == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				m_pVoidAVDF->m_fMinVal = 0;
				m_pVoidAVDF->m_fMaxVal = 1;
				m_pVoidAVDF->m_iResolution = bin2;
				m_oaVoidAVDF.push_back(m_pVoidAVDF);

				m_pVoidAV2DF->m_fMinVal[0] = 0;
				m_pVoidAV2DF->m_fMaxVal[0] = 1;
				m_pVoidAV2DF->m_iRes[0] = bin2;

				// Change MB 16.05.2020
				// This produced an "use of possibly undefined variable 'm_pVoidSphDF' error.
				//m_pVoidAV2DF->m_fMinVal[1] = m_pVoidSphDF->m_fMinVal;
				//m_pVoidAV2DF->m_fMaxVal[1] = m_pVoidSphDF->m_fMaxVal;
				m_pVoidAV2DF->m_fMinVal[1] = m_dMinRadius;
				m_pVoidAV2DF->m_fMaxVal[1] = save_d1;

				m_oaVoidAV2DF.push_back(m_pVoidAV2DF);

				m_pVoidDomDF->m_fMinVal = 0;
				m_pVoidDomDF->m_fMaxVal = d1;
				m_pVoidDomDF->m_iResolution = bin1;
				m_oaVoidDomDF.push_back(m_pVoidDomDF);

				if ( m_oaDomains.size() == 1 )
					break;
			}	
		}
	}
	else
		m_bVoidDomDF = false;
	m_iInterDom = 1;
	if ( m_bHoleDF || ( m_bVoidDomDF && AskYesNo("    Perform Domain decomposition analysis (y/n)? [yes] ",true) ) )
		m_iInterDom = AskUnsignedInteger("    Use how many intervals? [24] ",24);
	if ( m_iInterDom > 50 )
		mprintf(YELLOW,"\n    Warning! Too large interval numbers may cause memory problems!\n");
	if ( m_bVoidDomDF )
	{
		for ( i=0; i < (int)m_oaDomains.size()+1; i++ )
		{
			if ( m_iInterDom == 1 )
			{
				m_oaVoidAVDF[i]->Create();
				m_oaVoidDomDF[i]->Create();
			}
			else
			{
				m_oaVoidAVDF[i]->CreateMulti(m_iInterDom);
				m_oaVoidDomDF[i]->CreateMulti(m_iInterDom);
				m_oaVoidAV2DF[i]->m_iRes[1] = m_iInterDom;
				m_oaVoidAV2DF[i]->Create();

				cur_rad = m_dMinRadius;
				for ( j=0; j<m_iInterDom; j++ )
				{
					// Change MB 16.05.2020
					// This produced an "use of possibly undefined variable 'm_pVoidSphDF' error.
					//cur_rad += (m_pVoidSphDF->m_fMaxVal - m_dMinRadius) / m_iInterDom;
					cur_rad += (save_d1 - m_dMinRadius) / m_iInterDom;
					sprintf(temp, "%10.1f",cur_rad);
					m_oaVoidAVDF[i]->SetLabelMulti(j,temp);
					m_oaVoidDomDF[i]->SetLabelMulti(j,temp);
				}
			}
			
			m_oaVoidAV2DF[i]->SetLabelX("AV factor");
			m_oaVoidAV2DF[i]->SetLabelY("Threshold");
			m_oaVoidAV2DF[i]->SetLabelZ("Occurrence");

			m_oaVoidAVDF[i]->SetLabelX("AV factor");
			m_oaVoidAVDF[i]->SetLabelY("Occurrence");
			m_oaVoidDomDF[i]->SetLabelX("Void volume");
			m_oaVoidDomDF[i]->SetLabelY("Occurrence");

			if ( m_oaDomains.size() == 1 )
				break;
		}
	}
	if ( !m_bHoleDF )
		m_bVoidACF = AskYesNo("    Do you want to analyse the void fluidity (y/n)? [yes] ",true);
	else
		m_bVoidACF = false;

	if ( m_bVoidACF )
	{
		m_iSpotGrid = AskUnsignedInteger("    How many grid spots per dimension should be used? [25]", 25);
		if ( m_iInterDom > 1 && m_iSpotGrid > 1 )
			mprintf(YELLOW,"\n    Warning! Large interval numbers in combination with a high number of grid points may cause memory problems!\n");
		m_iDepth = g_iTrajSteps * 3 / 4;

		m_iDepth = AskUnsignedInteger("    Enter the correlation depth of the ACF (in trajectory frames): [%d] ", m_iDepth, m_iDepth);

		g_fTimestepLength = AskFloat("    Enter the physical time distance between successive trajectory frames in fs: [0.5] ",0.5);
	}

	if ( !m_bHoleDF )
		m_bMakePics = AskYesNo("    Do you want to create xyz files of the void centers for the first simulation step (y/n)? [no] ",false);
	else
		m_bMakePics = false;

	try { g_pVoroWrapper = new CVoroWrapper(); } catch(...) { g_pVoroWrapper = NULL; }
	if (g_pVoroWrapper == NULL) NewException((double)sizeof(CVoroWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mprintf("    Initializing Voronoi tesselation...\n");
	g_pVoroWrapper->Init();
	if ( AskYesNo("    Scale all radii by a constant factor (y/n)? [no]", false) )
	{
		m_dScale = AskFloat("    Enter factor [1.000]",1.0);
		mprintf("\n    Voronoi radii redefined.\n");
	}
	else
		m_dScale = 1.0;
	mprintf("*** Voro: Box density is %f particles / Angstrom^3.\n",g_pVoroWrapper->m_fBoxDens);
	mprintf("*** Voro: Using %d x %d x %d blocks.\n",g_pVoroWrapper->m_iBlocksX,g_pVoroWrapper->m_iBlocksY,g_pVoroWrapper->m_iBlocksZ);
}

CVoidObservation::~CVoidObservation()
{
}

void CVoidObservation::initialize()
{
	int i,n;
	CxIntArray* ia;
	double max;

	g_bVoroSilent = true;
	m_dCutOff *= m_dCutOff;
	m_iPics = 0;
	m_oaVolume.resize(m_iInterDom);

	if ( m_bVoidACF )
	{
		m_oaACFSpots.resize(m_iInterDom);
		for ( i=0; i < m_iInterDom; i++ )
		{
			for ( n=0; n < pow(m_iSpotGrid,3); n++ )
			{
				try { ia = new CxIntArray(); } catch(...) { ia = NULL; }
				if (ia == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				m_oaACFSpots[i].push_back(ia);
			}
		}
		if ( g_iMaxStep != -1 )
			n = g_iMaxStep;
		else
			n = g_iTrajSteps;
	
		if ( m_iDepth > n / g_iStride )
		{
			m_iDepth = int(n / g_iStride);
			mprintf("    Correlation Depth is exceeding number of total timesteps. Truncating to %d ....\n",m_iDepth);
		}
	}

	max = 0.0;
	for ( n=0; n < g_faVoronoiRadii.GetSize(); n++ )
		if ( g_faVoronoiRadii[n] * m_dScale > max )
			max = g_faVoronoiRadii[n] * m_dScale;
	m_dMinBox = m_oaVoidSphDF[0]->m_fMaxVal + max;
}

void CVoidObservation::process(CTimeStep* ts)
{
//	double time_temp;

	int n, m, i, j, k, i2=0, i3;
	container_periodic_poly *con;
	voronoicell_neighbor c;
	CxDVector3 Vec, Vec2;	
	vector< CxDVec3Array > all_v;
	vector< CxIntArray >   all_i;
	double *pp, cur_rad;
	int q, ijk, id = 0, faces, boxes[3], xyz[3], nos=0;
	vector< vector<CSphere*> > temp_voids;
	vector<CSphere*> grows;
	vector<CSphere*> allatoms;
	CSphere* atom;
	vector<int> nb;
	vector<double> vert, fa;
	CDomain *domain;
	int nmi;
	bool partofvoid, found;


	boxes[0] = (int)ceil(g_fBoxX / m_dMinBox);
	m_dBox[0] = g_fBoxX / boxes[0];
	boxes[1] = (int)ceil(g_fBoxY / m_dMinBox);
	m_dBox[1] = g_fBoxY / boxes[1];
	boxes[2] = (int)ceil(g_fBoxZ / m_dMinBox);
	m_dBox[2] = g_fBoxZ / boxes[2];

	m_iPics++;

	for ( n=0; n < (int)spheres.size(); n++ )
		spheres[n].clear();
	spheres.clear();
	spheres.resize(boxes[0]*boxes[1]*boxes[2]);
	all_v.resize(boxes[0]*boxes[1]*boxes[2]);
	all_i.resize(boxes[0]*boxes[1]*boxes[2]);
	temp_voids.resize(boxes[0]*boxes[1]*boxes[2]);
	for ( n=0; n < (int)VoidKey.size(); n++ )
		VoidKey[n].clear();
	VoidKey.clear();
	VoidKey.resize(boxes[0]*boxes[1]*boxes[2]);

	try { con = new container_periodic_poly(g_fBoxX/1000.0,0,g_fBoxY/1000.0,0,0,g_fBoxZ/1000.0,g_pVoroWrapper->m_iBlocksX,g_pVoroWrapper->m_iBlocksY,g_pVoroWrapper->m_iBlocksZ,g_iVoroMemory); } catch(...) { con = NULL; }
	if (con == NULL) NewException((double)sizeof(container_periodic_poly),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for ( n=0; n < int(ts->m_iGesAtomCount); n++ )
	{
		try { atom = new CSphere(); } catch(...) { atom = NULL; }
		if ( atom == NULL ) NewException((double)sizeof(CSphere),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		Vec = FoldVector(ts->m_vaCoords[n]);
		atom->dom = n;
		atom->m_vCenter = Vec;
		atom->m_dRadius = g_faVoronoiRadii[n] * m_dScale;
		allatoms.push_back(atom);

		if ( m_bNoH && ((CAtom*)g_oaAtoms[((CMolAtom*)g_oaMolAtoms[n])->m_iType])->m_pElement->m_iOrd == 1 )
			continue;

		Vec = FoldVectorPositive(Vec);
		for ( i=0; i<3; i++ )
			xyz[i] = (int)floor(Vec[i]/m_dBox[i]);	
		k = (xyz[0]*boxes[1]*boxes[2])+(xyz[1]*boxes[2])+xyz[2];
		spheres[k].push_back(atom);
		nos++;
		con->put(n,Vec[0]/1000.0,Vec[1]/1000.0,Vec[2]/1000.0,atom->m_dRadius/1000.0);
	}
//	time_temp = omp_get_wtime();
	c_loop_all_periodic vl1(*con);

	id = -1;
	if ( vl1.start() )
	{
		do
		{
			if ( con->compute_cell(c,vl1) )
			{
				id++;
				found = false;
				if ( m_bFull )
					found = true;
				else
					for ( i=0; i < (int)m_oaDomains.size(); i++ )
						if ( m_oaDomains[i]->Contains(id) )
							found = true;
				if ( !found )
					continue;
				ijk = vl1.ijk;
				q = vl1.q;
				pp = con->p[ijk] + con->ps * q;
				c.vertices(pp[0],pp[1],pp[2],vert);

				for ( n=0; n < int(vert.size())/3; n++ )
				{
					for ( m=0; m<3; m++ )
						Vec[m] = (vert[n*3+m]) * 1000;
					Vec = FoldVector(Vec);

					Vec2 = FoldVectorPositive(Vec);
					for ( i=0; i<3; i++ )
						xyz[i] = (int)floor(Vec2[i]/m_dBox[i]);	
					k = (xyz[0]*boxes[1]*boxes[2])+(xyz[1]*boxes[2])+xyz[2];
					all_v[k].Add(Vec);
					all_i[k].Add(all_i[k].GetSize());
				}
			}	
		} while ( vl1.inc() );
	}

	delete con;
//	time_voro1 += omp_get_wtime() - time_temp;

	m_oaVoids.clear();

	for ( k=0; k < (int)all_i.size(); k++ )
	{
		nb.clear();
		FindXYZ(k,&xyz[0],&boxes[0]);
		
		for ( i=0; i < 3; i++ )
			xyz[i] -= 2;

		for ( i=1; i < 4; i++ )
		{
			xyz[0]++;
			if ( xyz[0] < 0 )
				xyz[0] += boxes[0];
			if ( xyz[0] == boxes[0] )
				xyz[0] -= boxes[0];
			for ( i2=1; i2 < 4; i2++ )
			{
				xyz[1]++;
				if ( xyz[1] < 0 )
					xyz[1] += boxes[1];
				if ( xyz[1] == boxes[1] )
					xyz[1] -= boxes[1];
				for ( i3=1; i3 < 4; i3++ )
				{
					xyz[2]++;
					if ( xyz[2] < 0 )
						xyz[2] += boxes[2];
					if ( xyz[2] == boxes[2] )
						xyz[2] -= boxes[2];
					j = (xyz[0]*boxes[1]*boxes[2])+(xyz[1]*boxes[2])+xyz[2];
					nb.push_back(j);
					if ( i3 == boxes[2] )
						break;
				}
				if ( boxes[2] > 2 )
					xyz[2] -= 3;
				else
					xyz[2] -= boxes[2];
				if ( i2 == boxes[1] )
					break;
			}
			if ( boxes[1] > 2 )
				xyz[1] -= 3;
			else
				xyz[1] -= boxes[1];
			if ( i == boxes[0] )
				break;
		}
		for ( j=0; j < all_i[k].GetSize(); j++ )
		{
			if ( all_i[k][j] == -1 )
				continue;
			for ( i=0; i < (int)nb.size(); i++ )
			{
				for ( n=0; n < all_i[nb[i]].GetSize(); n++ )
				{
					if ( all_i[nb[i]][n] == -1 || (nb[i] == k && n == j))
						continue;
					Vec = FoldVector( all_v[nb[i]][all_i[nb[i]][n]] - all_v[k][all_i[k][j]] );
					if ( Vec.GetLengthSqr() < m_dCutOff )
						all_i[nb[i]][n] = -1;
				}
			}
		}
		grows.clear();
		for ( i=0; i < (int)nb.size(); i++ )
			grows.insert(grows.end(), spheres[nb[i]].begin(), spheres[nb[i]].end()); 
		for ( j=0; j < all_i[k].GetSize(); j++ )
		{
			if ( all_i[k][j] == -1 )
				continue;
			try { atom = new CSphere(); } catch(...) { atom = NULL; }
			if ( atom == NULL ) NewException((double)sizeof(CSphere),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			atom->m_vCenter = all_v[k][all_i[k][j]];
			if ( !atom->GrowIn(&grows,m_dMinRadius,m_oaVoidSphDF[0]->m_fMaxVal) || atom->m_dRadius < m_dMinRadius)
				all_i[k][j] = -1;
			else
				temp_voids[k].push_back(atom);
		}
		for ( j=0; j < (int)temp_voids[k].size(); j++ )
		{
//			time_temp = omp_get_wtime();
			for ( i=0; i < (int)nb.size(); i++ )
			{
				for ( i2=0; i2 < (int)temp_voids[nb[i]].size(); i2++ )
				{
					if ( (nb[i] == k && i2 == j) || temp_voids[nb[i]][i2]->m_dRadius < 0 ) 
						continue;
					Vec = FoldVector( temp_voids[nb[i]][i2]->m_vCenter - temp_voids[k][j]->m_vCenter );
					if ( Vec.GetLengthSqr() < m_dCutOff )
					{
						i2 = -1;
						break;
					}
				}
				if ( i2 == -1 )
					break;
			}
	
//HIER
			if ( i2 != -1 )
			{
				Vec = FoldVectorPositive(temp_voids[k][j]->m_vCenter);
				for ( i=0; i<3; i++ )
					xyz[i] = (int)floor(Vec[i]/m_dBox[i]);	
				VoidKey[(xyz[0]*boxes[1]*boxes[2])+(xyz[1]*boxes[2])+xyz[2]].push_back((int)m_oaVoids.size());
				m_oaVoids.push_back(temp_voids[k][j]);
			}
			else
				temp_voids[k][j]->m_dRadius = -1;
//			time_crit += omp_get_wtime() - time_temp;
		}
	}
	if ( m_bMakePics )
	{
		FILE* outpic;
		char name[256];
		int nov;
	
		outpic = OpenFileWrite("ref.xyz",true);
		fprintf(outpic,"%lu\n\n",(unsigned long)allatoms.size());
		for ( n=0; n < (int)allatoms.size(); n++ )
			fprintf(outpic,"%s %f %f %f\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[n]])->m_pElement->m_sLabel,
								allatoms[n]->m_vCenter[0]/100,
								allatoms[n]->m_vCenter[1]/100,
								allatoms[n]->m_vCenter[2]/100);
	
		cur_rad = m_dMinRadius;
	
		for ( k=0; k < m_iInterDom; k++ )
		{
			cur_rad += (m_oaVoidSphDF[0]->m_fMaxVal - m_dMinRadius) / m_iInterDom;
	
			nov = 0;
			for ( n=0; n < (int)m_oaVoids.size(); n++ )
				if ( m_oaVoids[n]->m_dRadius > cur_rad )
					nov++;
			if ( nov == 0 )
				continue;
			sprintf(name,"%.3f.xyz",cur_rad);
			outpic = OpenFileWrite(name,true);
			fprintf(outpic,"%d\n\n",nov);
			for ( n=0; n < (int)m_oaVoids.size(); n++ )
				if ( m_oaVoids[n]->m_dRadius > cur_rad )
					fprintf(outpic,"Ne %f %f %f\n",m_oaVoids[n]->m_vCenter[0]/100,
								       m_oaVoids[n]->m_vCenter[1]/100,
								       m_oaVoids[n]->m_vCenter[2]/100);
			fclose(outpic);
		}
		m_bMakePics = false;
	}

	for ( m=0; m < (int)m_oaVoids.size(); m++ )
		for ( n=0; n < (int)m_oaDomains.size(); n++ )
			if ( m_oaDomains[n]->Contains( m_oaVoids[m]->m_iaTouch[0] ) 
				&& m_oaDomains[n]->Contains( m_oaVoids[m]->m_iaTouch[1] )
				&& m_oaDomains[n]->Contains( m_oaVoids[m]->m_iaTouch[2] )
				&& m_oaDomains[n]->Contains( m_oaVoids[m]->m_iaTouch[3] ))
					m_oaVoids[m]->dom = n;
	for ( m=0; m < (int)m_oaVoids.size(); m++ )
		if ( m_oaVoids[m]->dom == -1 )
			m_oaVoids[m]->dom = (int)m_oaDomains.size();

	if ( m_bVoidSphDF )
		for ( i=0; i < (int)m_oaDomains.size()+1; i++ )
			for ( n=0; n < (int)m_oaVoids.size(); n++ )
				if ( m_oaVoids[n]->dom == i )
					m_oaVoidSphDF[i]->AddToBin( m_oaVoids[n]->m_dRadius );

	cur_rad = m_dMinRadius;

// Ab hier ist ein wenig Kaese ...
//	for ( j=0; j < (int)m_oaDomains.size()+1; j++ )
//	{
	for ( i=0; i < m_iInterDom; i++ )
	{
		cur_rad += (m_oaVoidSphDF[0]->m_fMaxVal - m_dMinRadius) / m_iInterDom;

		try { con = new container_periodic_poly(g_fBoxX/1000.0,0,g_fBoxY/1000.0,0,0,g_fBoxZ/1000.0,g_pVoroWrapper->m_iBlocksX,g_pVoroWrapper->m_iBlocksY,g_pVoroWrapper->m_iBlocksZ,g_iVoroMemory); } catch(...) { con = NULL; }
		if (con == NULL) NewException((double)sizeof(container_periodic_poly),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		m=0;
		for ( k=0; k < (int)spheres.size(); k++ )
		{
			for ( n=0; n < (int)spheres[k].size(); n++ )
			{
				con->put(m,spheres[k][n]->m_vCenter[0]/1000,spheres[k][n]->m_vCenter[1]/1000,spheres[k][n]->m_vCenter[2]/1000,spheres[k][n]->m_dRadius/1000);
				m++;
			}
		}
		for ( n=0; n < (int)m_oaVoids.size(); n++ )
		{
			if ( m_oaVoids[n]->m_dRadius >= cur_rad && m_oaVoids[n]->m_bActive )
			{
				con->put(m,m_oaVoids[n]->m_vCenter[0]/1000,m_oaVoids[n]->m_vCenter[1]/1000,m_oaVoids[n]->m_vCenter[2]/1000,m_oaVoids[n]->m_dRadius/1000);
				m++;
			}
		}
	
//		time_temp = omp_get_wtime();
		m_oaVoroVoids.clear();

		c_loop_all_periodic vl2(*con);
		voronoicell_neighbor c2;

		if (vl2.start()) 
		{
			do
			{
				ijk = vl2.ijk;
				q = vl2.q;
				id = con->id[ijk][q];
				if ( id < nos )
					continue;
				if (con->compute_cell(c2,vl2))
				{
					c2.neighbors(nb);
					c2.face_areas(fa);
					faces = c2.number_of_faces();

					try { domain = new CDomain(); } catch(...) { domain = NULL; }
					if ( domain == NULL ) NewException((double)sizeof(CDomain),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					domain->m_iaCells.Add(id);
					domain->m_fSurfaceArea = 0;
					domain->m_iFaces = 0;

					for ( n=0; n < faces; n++ )
					{
						if ( nb[n] > nos )
							domain->m_iaNeighbors.Add(nb[n]);
						else
							domain->m_fSurfaceArea += fa[n]*100.0;
					}
					domain->m_fVolume = c2.volume()*1000.0;
					m_oaVoroVoids.push_back(domain);
					m_oaVolume[i] += domain->m_fVolume; //HIERHER	
				}
			} while ( vl2.inc() );
		}
////////////////////////HIER/////////////////////////////////////////HIER

		delete con;
//		time_voro2 += omp_get_wtime() - time_temp;

		if ( Start_Domainize() )
		{
			for ( n=0; n < (int)m_oaDomVoids.size(); n++ )
			{
				j = m_oaVoids[m_oaDomVoids[n]->m_iaCells[0]-nos]->dom;
				for ( q=1; q < m_oaDomVoids[n]->m_iaCells.GetSize(); q++ )
					if ( j != m_oaVoids[m_oaDomVoids[n]->m_iaCells[q]-nos]->dom )
					{
						j = (int)m_oaDomains.size();
						break;
					}
				if ( m_oaDomVoids[n]->m_fAVRatio >= m_dAV )
					m_oaHoleDF[j]->AddToBin( pow( m_oaDomVoids[n]->m_fVolume * 0.238732414637843, 1/3.0 )*100.0 );
				if ( m_bVoidDomDF )
				{
					if ( m_iInterDom == 1 )
					{
						m_oaVoidAVDF[j]->AddToBin(  m_oaDomVoids[n]->m_fAVRatio );
						m_oaVoidDomDF[j]->AddToBin( m_oaDomVoids[n]->m_fVolume );
					}
					else
					{
						m_oaVoidAVDF[j]->AddToBin_Multi( i, m_oaDomVoids[n]->m_fAVRatio );
						m_oaVoidAV2DF[j]->AddToBin_IntY(    m_oaDomVoids[n]->m_fAVRatio, m_iInterDom-1-i );
						m_oaVoidDomDF[j]->AddToBin_Multi( i, m_oaDomVoids[n]->m_fVolume );
					}
				}
			}
			m_oaDomVoids.clear();
		}
		for ( n=0; n < (int)m_oaVoroVoids.size(); n++ )
			if ( m_oaVoroVoids[n] != NULL )
				delete m_oaVoroVoids[n];
		m_oaVoroVoids.clear();
	}

	if ( m_bVoidACF )
	{
	        cur_rad = m_dMinRadius;
	        for ( k=0; k < m_iInterDom; k++ )
       		{
			cur_rad += (m_oaVoidSphDF[0]->m_fMaxVal - m_dMinRadius) / m_iInterDom;

			for ( n=0; n < m_iSpotGrid; n++ )
			{
				for ( m=0; m < m_iSpotGrid; m++ )
				{
					for ( i=0; i < m_iSpotGrid; i++ )
					{
						nmi = n * m_iSpotGrid * m_iSpotGrid + m * m_iSpotGrid + i;
						Vec = CxDVector3(g_fBoxX/m_iSpotGrid*n,g_fBoxY/m_iSpotGrid*m,g_fBoxZ/m_iSpotGrid*i);

						partofvoid = true;
	
						nb.clear();
						xyz[0] = (int)floor(boxes[0]/m_iSpotGrid*n); 
						xyz[1] = (int)floor(boxes[1]/m_iSpotGrid*m); 
						xyz[2] = (int)floor(boxes[2]/m_iSpotGrid*i); 
						
						for ( j=0; j < 3; j++ )
							xyz[j] -= 2;
				
						for ( j=1; j < 4; j++ )
						{
							xyz[0]++;
							if ( xyz[0] < 0 )
								xyz[0] += boxes[0];
							if ( xyz[0] == boxes[0] )
								xyz[0] -= boxes[0];
							for ( i2=1; i2 < 4; i2++ )
							{
								xyz[1]++;
								if ( xyz[1] < 0 )
									xyz[1] += boxes[1];
								if ( xyz[1] == boxes[1] )
									xyz[1] -= boxes[1];
								for ( i3=1; i3 < 4; i3++ )
								{
									xyz[2]++;
									if ( xyz[2] < 0 )
										xyz[2] += boxes[2];
									if ( xyz[2] == boxes[2] )
										xyz[2] -= boxes[2];
									nb.push_back( (xyz[0]*boxes[1]*boxes[2])+(xyz[1]*boxes[2])+xyz[2] );
									if ( i3 == boxes[2] )
										break;
								}
								if ( boxes[2] > 2 )
									xyz[2] -= 3;
								else
									xyz[2] -= boxes[2];
								if ( i2 == boxes[1] )
									break;
							}
							if ( boxes[1] > 2 )
								xyz[1] -= 3;
							else
								xyz[1] -= boxes[1];
							if ( j == boxes[0] )
								break;
						}
						partofvoid = CheckIfVoid(nb,&Vec,cur_rad);

						if ( ( m_oaACFSpots[k][nmi]->GetSize()%2 == 0 && partofvoid ) || ( m_oaACFSpots[k][nmi]->GetSize()%2 == 1 && !partofvoid ) )
							m_oaACFSpots[k][nmi]->Add(m_iPics);
					}
				}
			}
		}
	}

}

void CVoidObservation::finalize()
{
	double cur_rad, norm, temp_d;
	vector< vector<double> > space_filling;
	int i, n, m, FFTSize = m_iPics;
	unsigned long j;
	CFFT FFT_in, FFT_out;
	CxDoubleArray results;
	CGrace gc, gc2, *gc3;
	bool zero;
	char temp[30];


	mprintf(WHITE,"\n    Producing the following results:\n");
	if ( m_bVoidSphDF )
	{
		mprintf("      void_radii*agr    - Distribution of the void spheres\n");
		mprintf("                          found during the fitting procedure.\n");
		for ( j=0; j < m_oaVoidSphDF.size(); j++ )
		{
			if ( m_oaVoidSphDF.size() == 1 )
				sprintf(temp, ".");
			else	
				sprintf(temp, "_Domain_%lu.",j+1);
			m_oaVoidSphDF[j]->NormBinIntegral(1);
			m_oaVoidSphDF[j]->Write("void_radii",temp,"csv",false);
			m_oaVoidSphDF[j]->WriteAgr("void_radii",temp,"agr","",false);
		}
	}

	if ( m_bVoidDomDF )
	{
		mprintf("      void_AV*agr       - Occurrence of AV factors of the \n");
		mprintf("                          void domains formed of the combined\n");
		mprintf("                          voronoi cells of the vois spheres\n");
		if ( m_iInterDom == 1 )
		{
			for ( j=0; j < m_oaVoidAVDF.size(); j++ )
			{
				if ( m_oaVoidAVDF.size() == 1 )
					sprintf(temp, ".");
				else	
					sprintf(temp, "_Domain_%lu.",j+1);
				m_oaVoidAVDF[j]->NormBinIntegral(1);
				m_oaVoidAVDF[j]->Write("void_AV",temp,"csv",false);
				m_oaVoidAVDF[j]->WriteAgr("void_AV",temp,"agr","",false);
			}
/*			for ( j=0; j < m_oaVoidDomDF.size(); j++ )
			{
				if ( m_oaVoidDomDF.size() == 1 )
					sprintf(temp, ".");
				else	
					sprintf(temp, "_Domain_%lu.",j+1);
				m_oaVoidDomDF[j]->Write("void_dom",temp,"csv",false);
				m_oaVoidDomDF[j]->WriteAgr("void_dom",temp,"agr","",false);
			}*/
		}
		else
		{
			mprintf("      void_AV_ratio*agr - Ratio of the number of void domains\n");
			mprintf("                          with an AV factor larger %f (\"spheric\")\n",m_dAV); 
			mprintf("                          to the number of those with an AV factor\n");
			mprintf("                          smaller than %f (\"non-spheric\")\n",m_dAV); 
			mprintf("                          in dependence of the void sphere threshold\n");

			for ( j=0; j < m_oaVoidAVDF.size(); j++ )
			{
				if ( m_oaVoidAVDF.size() == 1 )
					sprintf(temp, ".");
				else	
					sprintf(temp, "_Domain_%lu.",j+1);
				m_oaVoidAVDF[j]->NormAllBin(1);
				m_oaVoidAVDF[j]->WriteMulti("void_AV",temp,"csv");
				m_oaVoidAVDF[j]->WriteMultiAgr("void_AV",temp,"agr","",false);
			}
			
/*			for ( j=0; j < m_oaVoidAV2DF.size(); j++ )
			{
				if ( m_oaVoidAV2DF.size() == 1 )
					sprintf(temp, "dom");
				else	
					sprintf(temp, "dom_Domain_%lu",j+1);
				m_oaVoidAV2DF[j]->Log();
				m_oaVoidAV2DF[j]->WriteGnuplotInput("void_AV_",temp,"",false);
			}

			for ( j=0; j < m_oaVoidDomDF.size(); j++ )
			{
				if ( m_oaVoidDomDF.size() == 1 )
					sprintf(temp, ".");
				else	
					sprintf(temp, "_Domain_%lu.",j+1);
				m_oaVoidDomDF[j]->WriteMulti("void_dom",temp,"csv");
				m_oaVoidDomDF[j]->WriteMultiAgr("void_dom",temp,"agr","",false);
			}*/

			for ( j=0; j < m_oaVoidSphDF.size(); j++ )
			{
				try { gc3 = new CGrace(); } catch(...) { gc3 = NULL; }
				if ( gc3 == NULL ) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);

				gc3->SetRangeX(m_dMinRadius,m_oaVoidSphDF[j]->m_fMaxVal);
				gc3->SetRangeY(0,1);
				gc3->MakeTicks();
				gc3->SetLabelX("void radii threshold / pm");
				gc3->SetLabelY("ratio of domain AV > 0.5");
				gc3->CurrentGraph()->m_bLegend = false;
				gc3->CurrentGraph()->m_bInvertXAxis = false;
				gc3->AddDataset();
				
				cur_rad = m_dMinRadius;
				for ( i=0; i < m_iInterDom; i++ )
				{
					cur_rad += (m_oaVoidSphDF[j]->m_fMaxVal - m_dMinRadius) / m_iInterDom ;

					norm = 0.0;
					for ( n=0; n < m_oaVoidAVDF[j]->m_iResolution / 2; n++ )
						norm += m_oaVoidAVDF[j]->m_pMultiBin[i][n] / m_oaVoidAVDF[j]->m_iResolution;

					if ( !myisnan(norm) )
						gc3->AddXYTupel(0,cur_rad,1-norm);
				}
				if ( m_oaVoidSphDF.size() == 1 )
					sprintf(temp, "void_AV_ratio.agr");
				else	
					sprintf(temp, "void_AV_ratio_Domain_%lu.agr",j+1);
				gc3->WriteAgr(temp,false);
				if ( m_oaVoidSphDF.size() == 1 )
					sprintf(temp, "void_AV_ratio.csv");
				else	
					sprintf(temp, "void_AV_ratio_Domain_%lu.csv",j+1);
				gc3->WriteCSV(temp);
				delete gc3;
			}
		}
	}
/*
	space_filling.resize( m_oaVoidSphDF.size() );
	for ( j=0; j < m_oaVoidSphDF.size(); j++ )
	{
		space_filling[j].resize( m_iInterDom );
		if ( m_oaVoidSphDF.size() != 1 )
			mprintf(GREEN,"\n     Effective space filling in Domain %ld:",j+1);
		cur_rad = m_dMinRadius;
		for ( i=0; i < m_iInterDom; i++ )
		{
			cur_rad += (m_oaVoidSphDF[j]->m_fMaxVal - m_dMinRadius) / m_iInterDom ;
			m_oaVolume[j][i] /= m_iPics;
			space_filling[j][i] = 1 - m_oaVolume[j][i] / pow( g_fBoxX/100, 3 );
			mprintf(WHITE,"\n       With void sphere radius > %.1f:",cur_rad);
			mprintf(WHITE,"\n          Average effective void volume   = %10.2f A**3",m_oaVolume[j][i]);
			mprintf(WHITE,"\n          Average effective space filling = %10.2f %% \n",100*space_filling[j][i]);
		}
	}
*/
	if ( m_bVoidACF )
	{
		mprintf("      void_fluidity*agr - Auto-correlation functions of the void grid\n");
		mprintf("                          native and corrected by the assembly average\n"); 
		cur_rad = m_dMinRadius;
		for ( i=0; i < m_iInterDom; i++ )
		{
			cur_rad += (m_oaVoidSphDF[0]->m_fMaxVal - m_dMinRadius) / m_iInterDom ;
			while ( true )
			{
				n = FFTSize;
				while ( n%2 == 0 )
					n /= 2;
				while ( n%3 == 0 )
					n /= 3;
				while ( n%5 == 0 )
					n /= 5;
				if ( n == 1 )
					break;
				FFTSize++;
			}
	
			FFT_in.PrepareFFT_C2C(2*FFTSize);		
			FFT_out.PrepareInverseFFT_C2C(2*FFTSize);		
			results.RemoveAll();
			for ( n=0; n < FFTSize/2; n++ )
				results.Add(0.0);
	
			for ( n=0; n < (int)m_oaACFSpots[i].size(); n++ )
			{
				zero = true;
				
				for ( m=0; m < m_iPics; m++ )
				{
					if ( m_oaACFSpots[i][n]->Contains(m+1) )
						zero = !zero;
					if ( zero )
						FFT_in.m_pInput[m*2] = 0.0;
					else
						FFT_in.m_pInput[m*2] = 1.0;
					FFT_in.m_pInput[m*2+1] = 0.0;
				}
				for ( m=m_iPics; m < 2*FFTSize; m++ )
				{
					FFT_in.m_pInput[m*2] = 0.0;
					FFT_in.m_pInput[m*2+1] = 0.0;
				}
				FFT_in.DoFFT();
	
				for ( m=0; m < 2*FFTSize; m++ )
				{
					FFT_out.m_pInput[m*2] = FFT_in.m_pOutput[m*2]*FFT_in.m_pOutput[m*2] + FFT_in.m_pOutput[m*2+1]*FFT_in.m_pOutput[m*2+1];
					FFT_out.m_pInput[m*2+1] = 0.0;
				}
				FFT_out.DoFFT();
	
				for ( m=0; m < FFTSize/2; m++ )
					results[m] += FFT_out.m_pOutput[2*m] / ( m_iPics - m);
			}
			norm = results[0];
			for ( m=0; m < FFTSize/2; m++ )
				results[m] /= norm;
			sprintf(temp, "%10.1f",cur_rad);

			gc.SetRangeX(0,m_iPics/2*g_fTimestepLength*g_iStride);
			gc.SetRangeY(0,1);
			gc.MakeTicks();
			if ( g_fTimestepLength*g_iStride > 100 )
				gc.SetLabelX("t / ps");
			else
				gc.SetLabelX("t / fs");
			gc.SetLabelY("c(t)");
			gc.CurrentGraph()->m_bLegend = true;
			gc.CurrentGraph()->m_bInvertXAxis = false;
			gc.AddDataset();
			gc.SetDatasetName(i,temp);
			for ( n=0; n < results.GetSize(); n++ )
			{
				if ( g_fTimestepLength*g_iStride > 100 )
					gc.AddXYTupel(i,n*g_fTimestepLength*g_iStride/1000,results[n]);
				else
					gc.AddXYTupel(i,n*g_fTimestepLength*g_iStride,results[n]);
			}
			temp_d = m_oaVolume[i] / m_iPics / (g_fBoxX/100*g_fBoxY/100*g_fBoxZ/100); //HIERHER
			for ( m=0; m < FFTSize/2; m++ )
				results[m] -= temp_d;	
			norm = results[0];
			for ( m=0; m < FFTSize/2; m++ )
				results[m] /= norm;
			
			gc2.SetRangeX(0,m_iPics/2);
			gc2.SetRangeY(0,1);
			gc2.MakeTicks();
			if ( g_fTimestepLength*g_iStride > 100 )
				gc2.SetLabelX("t / ps");
			else
				gc2.SetLabelX("t / fs");
			gc2.SetLabelY("c(t)");
			gc2.CurrentGraph()->m_bLegend = true;
			gc2.CurrentGraph()->m_bInvertXAxis = false;
			gc2.AddDataset();
			gc2.SetDatasetName(i,temp);
			for ( n=0; n < results.GetSize(); n++ )
			{
				if ( g_fTimestepLength*g_iStride > 100 )
					gc2.AddXYTupel(i,n*g_fTimestepLength*g_iStride/1000,results[n]);
				else
					gc2.AddXYTupel(i,n*g_fTimestepLength*g_iStride,results[n]);
			}
		}
		gc.WriteAgr("void_fluidity.agr",false);
		gc.WriteCSV("void_fluidity.csv");
		gc2.WriteAgr("void_fluidity_cor.agr",false);
		gc2.WriteCSV("void_fluidity_cor.csv");
	}

	mprintf("      hole_radii*agr    - Distribution of the hole radii. A \"hole\" is defined\n");
	mprintf("                          as the radius of a sphere with the same volume as a\n"); 
	mprintf("                          void domain with a spheric shape (AV factor > %f)\n",m_dAV); 
	for ( j=0; j < m_oaHoleDF.size(); j++ )
	{
		if ( m_oaHoleDF.size() == 1 )
			sprintf(temp, ".");
		else	
			sprintf(temp, "_Domain_%ld.",j+1);
		m_oaHoleDF[j]->NormBinIntegral(1);
		m_oaHoleDF[j]->Write("hole_radii",temp,"csv",false);
		m_oaHoleDF[j]->WriteAgr("hole_radii",temp,"agr","",false);
	}

/*
	time_tot = omp_get_wtime() - time_tot;

	mprintf(RED,"\nTimings:\n  total: %7.3f s\n  voro1: %7.3f s (%7.1f%%)\n  voro2: %7.3f s (%7.1f%%)\n    dom: %7.3f s (%7.1f%%)\n   grow: %7.3f s (%7.1f%%)\n   crit: %7.3f s (%7.1f%%)\n",
		time_tot,
		time_voro1,time_voro1/time_tot*100.0,
		time_voro2,time_voro2/time_tot*100.0,
		time_dom,time_dom/time_tot*100.0,
		time_grow,time_grow/time_tot*100.0,
		time_crit,time_crit/time_tot*100.0
		);
*/
}

bool CVoidObservation::Start_Domainize()
{
//	double time_temp = omp_get_wtime();

	int n, next = 0;
	CxIntArray* stack;
	CDomain* domain, *bdomain;

	if ( m_oaVoroVoids.size() == 0 )
		return false;

	for ( n=0; n < (int)m_oaDomVoids.size(); n++ )
		if ( m_oaDomVoids[n] != NULL )
			delete m_oaDomVoids[n];
	m_oaDomVoids.clear();

	try { stack = new CxIntArray(); } catch(...) { stack = NULL; }
	if ( stack == NULL ) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	for ( n=0; n < (int)m_oaVoroVoids.size(); n++ )
	{
		stack->Add(m_oaVoroVoids[n]->m_iaCells[0]);
		m_oaVoroVoids[n]->m_bActive = true;
	}

	while ( next != -1 )
	{
		bdomain = m_oaVoroVoids[next];
		bdomain->m_bActive = false;

		for ( n=0; n < bdomain->m_iaNeighbors.GetSize(); n++ )
		{
			domain = m_oaVoroVoids[ stack->GetPosition( bdomain->m_iaNeighbors[n] ) ];
			if ( !domain->m_bActive )
				continue;
			domain->m_bActive = false;
			if ( !Domainize(stack, bdomain, domain ) )
				return false;
		}	

// MAGIC NUMBER:
		bdomain->m_fAVRatio = 10.6347231054330961 * bdomain->m_fVolume / mypow(bdomain->m_fSurfaceArea,1.5);
		m_oaDomVoids.push_back(bdomain);
		
		next = -1;
		for ( n=0; n < (int)m_oaVoroVoids.size(); n++ )
			if ( m_oaVoroVoids[n]->m_bActive )
			{
				next = stack->GetPosition(m_oaVoroVoids[n]->m_iaCells[0]);
				break;
			}
	}

	delete stack;

//	time_dom += omp_get_wtime() - time_temp;
	return true;
}

bool CVoidObservation::Domainize(CxIntArray* stack, CDomain* bdom, CDomain* dom)
{
	int n;
	CDomain* domain;

	for ( n=0; n < dom->m_iaNeighbors.GetSize(); n++ )
	{
		domain = m_oaVoroVoids[ stack->GetPosition( dom->m_iaNeighbors[n] ) ];
		if ( !domain->m_bActive )
			continue;
		domain->m_bActive = false;
		if ( !Domainize(stack, bdom, domain ) )
			return false;
	}
	bdom->Assimilate(dom);

	return true;
}

bool CVoidObservation::CheckIfVoid(vector<int> nb, CxDVector3* point, double treshold)
{
	unsigned int n, m;
	CxDVector3 Vec;	
	double dist_atom = 999999.9, dist_void= 999999.9;
	double temp;

	for ( n=0; n < nb.size(); n++ )
	{
		for ( m=0; m < spheres[nb[n]].size(); m++ )
		{
			Vec = FoldVector( *point - spheres[nb[n]][m]->m_vCenter );
			if ( Vec.GetLength() < spheres[nb[n]][m]->m_dRadius )
				return false;
			temp = Vec.GetLengthSqr() - spheres[nb[n]][m]->m_dRadius * spheres[nb[n]][m]->m_dRadius;
			if ( temp < dist_atom )
				dist_atom = temp;
		}	

		for ( m=0; m < VoidKey[nb[n]].size(); m++ )
		{
			if ( m_oaVoids[VoidKey[nb[n]][m]]->m_dRadius < treshold )
				continue;
			Vec = FoldVector( *point - m_oaVoids[VoidKey[nb[n]][m]]->m_vCenter );
			if ( Vec.GetLength() < m_oaVoids[VoidKey[nb[n]][m]]->m_dRadius )
				return true;
			temp = Vec.GetLengthSqr() - m_oaVoids[VoidKey[nb[n]][m]]->m_dRadius * m_oaVoids[VoidKey[nb[n]][m]]->m_dRadius;
			if ( temp < dist_void )
				dist_void = temp;
		}
	}
	if ( dist_atom < dist_void )
		return false;
	return true;
}

CSphere::CSphere()
{
	int n;

	m_bActive = true;
	m_dRadius = 0.0;
	dom = -1;
	for ( n=0; n<3; n++ )
		m_vCenter[n] = -999.9;
}

void CSphere::FindCenter(CxDVector3* a,CxDVector3* b,CxDVector3* c,CxDVector3* d)
{
	double A[4][5];
	double faktor;
	int n, m, i;

	for ( n=0; n<4; n++ )
		A[n][0] =   1.0;
	for ( n=1; n<4; n++ )
	{
		A[0][n] = (*a)[n-1];
		A[1][n] = (*b)[n-1];
		A[2][n] = (*c)[n-1];
		A[3][n] = (*d)[n-1];
	}
	A[0][4] = -(a->GetLengthSqr());
	A[1][4] = -(b->GetLengthSqr());
	A[2][4] = -(c->GetLengthSqr());
	A[3][4] = -(d->GetLengthSqr());

	for ( m=0; m<4; m++ )
	{
		for ( i=0; i<4; i++ )
		{
			if ( m != i )
			{
				faktor = A[i][m] / A[m][m];
				for ( n=1; n<5; n++ )
					A[i][n] -= A[m][n] * faktor;
			}
			else
			{
				faktor = A[i][m];
				for ( n=1; n<5; n++ )
					A[i][n] /= faktor;
			}
		}
	}
	
	for ( i=0; i<3; i++ )
		m_vCenter[i] = -A[i+1][4]/2;
	m_dRadius = sqrt(A[1][4]*A[1][4]/4 + A[2][4]*A[2][4]/4 + A[3][4]*A[3][4]/4 - A[0][4]);
}

bool CSphere::GrowIn(vector<CSphere*> *spheres, double MinRadius, double cutoff)//, double *time_grow)
{
//	double time_temp = omp_get_wtime();

	if ( spheres->size() < 4 )
	{
		eprintf("     Trying to grow into a \"hole\" of less then 4 atoms!\n");	
		return false;
	}

	double shift = SHIFT_START;
	int max_scf = 50, n;
	unsigned int m, i, j;
	int sign = 0;
	double dist;
	CxDVector3 Vec, Vec1(0,0,0), Vec2(0,0,0), Vec3(0,0,0);
	CxDVec3Array points;
	CSphere Last;
	CxIntArray m_iaNL;
	int touch[4];

	m_iaTouch.resize(4,-1);
	m_dRadius = g_fBoxX * 2;

	for ( n=0; n < (int)spheres->size(); n++ )
	{
		Vec = FoldVector(m_vCenter - (*spheres)[n]->m_vCenter);
		dist = Vec.GetLength() - (*spheres)[n]->m_dRadius;

		if ( dist < 0 )
			return false;

		if ( dist < cutoff + cutoff )
			m_iaNL.Add(n);

		if ( dist < m_dRadius )
		{
			m_dRadius = dist;
			touch[0] = (*spheres)[n]->dom;
			m_iaTouch[0] = n;
		}
	}
	Vec = FoldVector(m_vCenter - (*spheres)[m_iaTouch[0]]->m_vCenter);
	Vec.Normalize();

	m = 1;
	while ( m < 4 )
	{
		m_vCenter += shift * Vec;
		Vec1 = FoldVector(m_vCenter - (*spheres)[m_iaTouch[0]]->m_vCenter);
		if ( m == 1 )
			m_dRadius += shift;
		else
			m_dRadius = Vec1.GetLength() - (*spheres)[m_iaTouch[0]]->m_dRadius;
		if ( m_dRadius < cutoff )
			j = m_iaNL.GetSize();
		else
			j = (unsigned int)spheres->size();
		for ( i=0; i<j; i++ )
		{
			if ( m_dRadius < cutoff )
				n = m_iaNL[i];
			else
				n = i;
			if ( n != m_iaTouch[0] && n != m_iaTouch[1] && n != m_iaTouch[2] )
			{
				Vec2 = FoldVector(m_vCenter - (*spheres)[n]->m_vCenter);
				dist = Vec2.GetLength() - (*spheres)[n]->m_dRadius;
				if ( dist <= m_dRadius )
				{
					if ( m == 2 && m_dRadius < 0 )
						return false;
					if ( shift > 0.01 )
					{
						m_vCenter -= shift * Vec;
						m_dRadius -= shift;
						shift *= 0.2;
						break;
					}
				 	touch[m] = (*spheres)[n]->dom;
				 	m_iaTouch[m] = n;
					m++;
					shift = SHIFT_START;
					break;
				}
			}
		}
		if ( m > 1 )
		{
			Vec2 = FoldVector(m_vCenter - (*spheres)[m_iaTouch[1]]->m_vCenter);
			Vec1.Normalize();
			Vec2.Normalize();
		}
		if ( m == 2 )
			Vec = Vec1 + Vec2;
		if ( m == 3 )
		{
			Vec3 = FoldVector(m_vCenter - (*spheres)[m_iaTouch[2]]->m_vCenter);
			Vec3.Normalize();
			Vec1 = FoldVector((*spheres)[m_iaTouch[0]]->m_vCenter + Vec1 * (*spheres)[m_iaTouch[0]]->m_dRadius );
			Vec2 = FoldVector((*spheres)[m_iaTouch[1]]->m_vCenter + Vec2 * (*spheres)[m_iaTouch[1]]->m_dRadius );
			Vec  = FoldVector((*spheres)[m_iaTouch[2]]->m_vCenter + Vec3 * (*spheres)[m_iaTouch[2]]->m_dRadius );
			Vec  = CrossP(FoldVector(Vec - Vec1), FoldVector(Vec - Vec2));
			if ( sign == 0 )
			{
				if ( DotP(Vec,Vec3) < 0 )
					sign = -1;
				else
					sign = 1;
			}
			Vec *= sign;
		}
		if ( m > 1 )
			Vec.Normalize();
	}

	if ( m_dRadius < MinRadius/2 )
		return false;

	while ( max_scf > 0 )
	{
		max_scf--;

		points.RemoveAll();

		Last.m_dRadius = m_dRadius;
		Last.m_vCenter = m_vCenter;
	
		for ( m=0; m<4; m++ )		// Collecting 4 points (which are on the spheres in a line between m_vCenter and the m_vCenter of the sphere)
		{
			Vec = FoldVector(m_vCenter - (*spheres)[m_iaTouch[m]]->m_vCenter);
			dist = Vec.GetLength() - (*spheres)[m_iaTouch[m]]->m_dRadius;
			Vec.Normalize();
			points.Add(m_vCenter - dist * Vec);
		}
		
		FindCenter(&points[0],&points[1],&points[2],&points[3]);

		Vec = FoldVector(m_vCenter - Last.m_vCenter);
		if ( fabs(Last.m_dRadius - m_dRadius) < CUT_OFF && Vec.GetLengthSqr() < CUT_OFF )
		{
			for ( m=0; m<4; m++ )
				m_iaTouch[m] = touch[m];
//			(*time_grow) += omp_get_wtime() - time_temp;
			return true;
		}	
	}
//	(*time_grow) += omp_get_wtime() - time_temp;
	return false;
}

CSphere::~CSphere()
{
}

bool CVoidObservation::checkAtomChoice(const char *s, int mol, CxIntArray* ia)
{
	const char *p, *q;
	char buf[256];
	const char *allowed = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890-,_ ";
	int inclu, la = -1, i = -1, i2 = -1, atom = -1;
	bool m = false, fin = false;

	p = s;

	while (*p != 0)
	{
		fin = false;
		while (*p == ' ')
			p++;
		if ( strchr(allowed,*p) == NULL )
		{
			eprintf("Error: Character \"%c\" not allowed.\n",*p);
			return false;
		}
		q = p;
		if (isalpha(*q) || *q == '#')
		{
			if (m)
			{
				eprintf("Error: Only digit allowed after \"-\".\n");
				return false;
			}
			while ( isalpha(*q) )
				q++;
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
			atom = ((CMolecule*)g_oaMolecules[mol])->FindAtomInMol(buf);
			if (atom == -1)
			{
				eprintf("Error: Atom \"%s\" not in the molecule.\n",buf);
				return false;
			}
		}
		else if ( isdigit(*q) )
		{
			if (atom == -1)
				atom = la;
			if (atom == -1)
			{
				eprintf("Error: Number in the beginning not possible.\n");
				return false;
			}
			while ( isdigit(*q) )
				q++;
			if ((*q != '-') && (*q != ',') && (*q != ' ') && (*q != 0))
			{
				eprintf("Error: Only \",\" or \"-\" may follow after a number.\n");
				return false;
			}
			memcpy(buf,p,q-p);
			buf[q-p] = 0;
			if ( atoi(buf)-1 >= ((CMolecule*)g_oaMolecules[mol])->m_waAtomCount[atom] )
			{
				eprintf("Error: Only %d %s atoms in the molecule (requested: %d)\n",((CMolecule*)g_oaMolecules[mol])->m_waAtomCount[atom],(const char*)((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[mol])->m_baAtomIndex[atom]])->m_sName,atoi(buf));
				return false;
			}
			if ( m )
			{
				i2 = atoi(buf)-1;
				if ( i2 < i )
				{
					eprintf("Error: Invalid atom range, %d < %d.\n",i2+1,i+1);
					return false;
				}
				else
				{
					for ( inclu=i; inclu <= i2; inclu++ )
						Add(atom, inclu, mol, ia);
					fin = true;
				}
			}
			else
				i = atoi(buf)-1;
		}
		else if (*q == '-')
		{
			if ( i == -1 )
			{
				eprintf("Error: \"-\" without preceding number.\n");
				return false;
			}
			m = true;
			q++;
		}
		else if ( *q == ',' )
		{
			if ( atom == -1 )
			{
				eprintf("Error: Comma without atom.\n");
				return false;
			}
			if ( i == -1 )
				for ( inclu=0; inclu < ((CMolecule*)g_oaMolecules[mol])->m_waAtomCount[atom]; inclu++ )
					Add(atom, inclu, mol, ia);
			else if ( !m )
				Add(atom, i, mol, ia);
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
		if ( i == -1 )
			for ( inclu=0; inclu < ((CMolecule*)g_oaMolecules[mol])->m_waAtomCount[atom]; inclu++ )
				Add(atom, inclu, mol, ia);
		else
			Add(atom, i, mol, ia);
	}

	return true;
}

void CVoidObservation::AddAll(int mol, CxIntArray* ia)
{
	int i,j;

	for ( j=0; j < ((CMolecule*)g_oaMolecules[mol])->m_waAtomCount.GetSize()-1; j++ )
		for ( i=0; i < ((CMolecule*)g_oaMolecules[mol])->m_waAtomCount[j]; i++ )
			Add(j,i,mol,ia);
}

void CVoidObservation::Add(int atom, int include, int mol, CxIntArray* ia)
{
	CxIntArray includeME;
	int i;

	for ( i=0; i < ((CMolecule*)g_oaMolecules[mol])->m_laSingleMolIndex.GetSize(); i++ )
	{
		includeME.CopyFrom( (CxIntArray*)(((CSingleMolecule*)g_oaSingleMolecules[ ((CMolecule*)g_oaMolecules[mol])->m_laSingleMolIndex[i] ])->m_oaAtomOffset[atom]) );
		ia->Add(includeME[include]);
	}
}


void CVoidObservation::FindXYZ(int n, int *xyz, int* boxes)
{
	int i;

	for ( i=0; i<3; i++ )
		xyz[i] = 0;

	while ( n >= boxes[2]*boxes[1] )
	{
		n -= boxes[2]*boxes[1];
		xyz[0]++;
	}
	while ( n >= boxes[2] )
	{
		n -= boxes[2];
		xyz[1]++;
	}
	xyz[2] = n;
}



