/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

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

#include "timestep.h"
#include "travis.h"
#include "tools.h"
#include "maintools.h"
#include "conversion.h"
#ifdef USE_NETCDF
#include <netcdf.h>
#include <vector>
#endif


const char *GetRevisionInfo_timestep(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_timestep() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


bool CTimeStep::m_bFirst = true;


void CTimeStep::CalcCenters()
{
	BTIN;
	int z, z2, z3, z4, c, k;
	CxDVector3 tv;
	CMolecule *m;
	CVirtualAtom *v;
	CSingleMolecule *sm;

	tv = CxDVector3(0);

//	mprintf("** CalcCenters **");

	if (m_vaCoords.GetSize() < g_iGesVirtAtomCount)
		m_vaCoords.SetSize(g_iGesVirtAtomCount);

//	mprintf("CalcCenters(): Size is now %d.\n",m_vaCoords.GetSize());

	if (g_bUseVelocities)
		if (m_vaVelocities.GetSize() < g_iGesVirtAtomCount)
			m_vaVelocities.SetSize(g_iGesVirtAtomCount);

	if (g_bUseForces)
		if (m_vaForces.GetSize() < g_iGesVirtAtomCount)
			m_vaForces.SetSize(g_iGesVirtAtomCount);

	for (z=0;z<g_oaVirtualAtoms.GetSize();z++)
	{
		v = (CVirtualAtom*)g_oaVirtualAtoms[z];
		m = (CMolecule*)g_oaMolecules[v->m_iMolecule];

		// Find virtual atom type in this molecule
		for (k = 0; k < m->m_baAtomIndex.GetSize(); k++)
			if (m->m_baAtomIndex[k] == g_iVirtAtomType)
				break;
		
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			if (v->m_iMode == 0) // Mittel aus Atompositionen
			{
				tv = CxDVector3(0);
				c = 0;
				for (z3=0;z3<v->m_oCenterAtoms.m_baAtomType.GetSize();z3++)
				{
					for (z4=0;z4<((CxIntArray*)v->m_oCenterAtoms.m_oaAtoms[z3])->GetSize();z4++)
					{
//						mprintf("{%d:(%G|%G|%G) }\n",((CxIntArray*)sm->m_oaAtomOffset[v->m_oCenterAtoms.m_baAtomType[z3]])->GetAt(((CxIntArray*)v->m_oCenterAtoms.m_oaAtoms[z3])->GetAt(z4)),m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[v->m_oCenterAtoms.m_baAtomType[z3]])->GetAt(((CxIntArray*)v->m_oCenterAtoms.m_oaAtoms[z3])->GetAt(z4))][0],m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[v->m_oCenterAtoms.m_baAtomType[z3]])->GetAt(((CxIntArray*)v->m_oCenterAtoms.m_oaAtoms[z3])->GetAt(z4))][1],m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[v->m_oCenterAtoms.m_baAtomType[z3]])->GetAt(((CxIntArray*)v->m_oCenterAtoms.m_oaAtoms[z3])->GetAt(z4))][2]);
						tv += m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[v->m_oCenterAtoms.m_baAtomType[z3]])->GetAt(((CxIntArray*)v->m_oCenterAtoms.m_oaAtoms[z3])->GetAt(z4))] * v->m_faWeight[c];
						c++;
					}
				}
				tv /= v->m_fGesWeight;
			} else if (v->m_iMode == 1) // Abstand, Winkel, Dihedralwinkel
			{
				tv = PointFromRAD(m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[v->m_iAtomType[0]])->GetAt(v->m_iAtom[0])],m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[v->m_iAtomType[1]])->GetAt(v->m_iAtom[1])],m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[v->m_iAtomType[2]])->GetAt(v->m_iAtom[2])],v->m_fValues[0],v->m_fValues[1],v->m_fValues[2]);
			} else if (v->m_iMode == 2) // Dipolvektor
			{
/*				if (!g_bDipole)
				{
					eprintf("Cannot use dipole vectors.\n");
					BTOUT;
					return;
				}
				tv = sm->m_vDipole + m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[sm->m_baAtomIndex.GetSize()-1])->GetAt(0)];*/

				//tv = CxDVector3(0,0,0); // Dipolmoment wird erst spaeter ausgerechnet

				// Fix MB Mar 2020
				// #3 is equivalent to #2 (equals zero dipole vector)
				tv = m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[k])->GetAt(1)];

			} else if (v->m_iMode == 3) // Geschwindigkeitsvektor
			{
			} else if (v->m_iMode == 4) // Kraftvektor
			{
			}
//			mprintf("%d:(%G|%G|%G), ",((CxIntArray*)sm->m_oaAtomOffset[sm->m_baAtomIndex.GetSize()-1])->GetAt(v->m_iMolVirtAtom),tv[0],tv[1],tv[2]);
			m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[sm->m_baAtomIndex.GetSize()-1])->GetAt(v->m_iMolVirtAtom)] = tv;
		}
	}
	BTOUT;
}


void CTimeStep::CalcCenterVelocities()
{
	int z, z2, z3, z4, c;
	CxDVector3 tv;
	CMolecule *m;
	CVirtualAtom *v;
	CSingleMolecule *sm;

	if (g_bUseVelocities)
		if (m_vaVelocities.GetSize() < g_iGesVirtAtomCount)
			m_vaVelocities.SetSize(g_iGesVirtAtomCount);

	for (z=0;z<g_oaVirtualAtoms.GetSize();z++) {
		v = (CVirtualAtom*)g_oaVirtualAtoms[z];
		m = (CMolecule*)g_oaMolecules[v->m_iMolecule];
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++) {
	//		mprintf("\nCenter %d.%d: ",z+1,z2+1);
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			if (v->m_iMode == 0) { // Mittel aus Atompositionen
				tv = CxDVector3(0);
				c = 0;
				for (z3=0;z3<v->m_oCenterAtoms.m_baAtomType.GetSize();z3++) {
					for (z4=0;z4<((CxIntArray*)v->m_oCenterAtoms.m_oaAtoms[z3])->GetSize();z4++) {

	/*					mprintf(" + %5g * (%10g|%10g|%10g)",
							v->m_faWeight[c],
							m_vaVelocities[((CxIntArray*)sm->m_oaAtomOffset[v->m_oCenterAtoms.m_baAtomType[z3]])->GetAt(((CxIntArray*)v->m_oCenterAtoms.m_oaAtoms[z3])->GetAt(z4))][0],
							m_vaVelocities[((CxIntArray*)sm->m_oaAtomOffset[v->m_oCenterAtoms.m_baAtomType[z3]])->GetAt(((CxIntArray*)v->m_oCenterAtoms.m_oaAtoms[z3])->GetAt(z4))][1],
							m_vaVelocities[((CxIntArray*)sm->m_oaAtomOffset[v->m_oCenterAtoms.m_baAtomType[z3]])->GetAt(((CxIntArray*)v->m_oCenterAtoms.m_oaAtoms[z3])->GetAt(z4))][2]
							);*/

						tv += m_vaVelocities[((CxIntArray*)sm->m_oaAtomOffset[v->m_oCenterAtoms.m_baAtomType[z3]])->GetAt(((CxIntArray*)v->m_oCenterAtoms.m_oaAtoms[z3])->GetAt(z4))] * v->m_faWeight[c];
						c++;
					}
				}
				tv /= v->m_fGesWeight;

	/*			mprintf("\n");
				mprintf("  --> (%10g|%10g|%10g)",tv[0],tv[1],tv[2]);*/

				m_vaVelocities[((CxIntArray*)sm->m_oaAtomOffset[sm->m_baAtomIndex.GetSize()-1])->GetAt(v->m_iMolVirtAtom)] = tv;
			} else	
				m_vaVelocities[((CxIntArray*)sm->m_oaAtomOffset[sm->m_baAtomIndex.GetSize()-1])->GetAt(v->m_iMolVirtAtom)] = CxDVector3(0);
		}
	}
}


bool CTimeStep::BondRange(int i1, int i2, double *f)
{
	BTIN;
	int a, b;
	double d;
	CxDVector3 v;

	if (g_iScanMolStep == -1)
		return false;

	v = FoldVector(m_vaCoords[i1] - m_vaCoords[i2]);
	d = v.GetLength();

	if (f != NULL)
		*f = d;

	a = g_baAtomIndex[i1];
	b = g_baAtomIndex[i2];

	if ((a >= g_oaAtoms.GetSize()) || (b >= g_oaAtoms.GetSize()))
	{
		BTOUT; 
		return false;
	}

	if (((CAtom*)g_oaAtoms[a])->m_pElement->m_fRadius == 0)
	{
		BTOUT;
		return false;
	}
	
	if (((CAtom*)g_oaAtoms[b])->m_pElement->m_fRadius == 0)
	{
		BTOUT;
		return false;
	}

	if (d < (((CAtom*)g_oaAtoms[a])->m_pElement->m_fRadius+((CAtom*)g_oaAtoms[b])->m_pElement->m_fRadius)*g_fBondFactor)
	{
//		printf("    \"%s\" br=%f, \"%s\" br=%f Hat Nachbarn %s%d im Abstand von %.3f.\n",m_pLabels[i1],c1,m_pLabels[i2],c2,m_pLabels[i2],i2+1,d);
		BTOUT; 
		return true;
	} else 
	{
		BTOUT;
		return false;
	}
}


bool CTimeStep::MirrorBond(int i1, int i2)
{
	BTIN;
	bool changed;
	int n;

	changed = false;

	if (g_bBoxNonOrtho)
	{
		CxDVector3 w1, w2;

		w1 = g_mBoxToOrtho * m_vaCoords[i1];
		w2 = g_mBoxToOrtho * m_vaCoords[i2];

		n = PeriodicImage1D( w2[0]-w1[0], -0.5, 0.5 );
		if (n != 0) {
			w2[0] -= n;
			changed = true;
		}

		n = PeriodicImage1D( w2[1]-w1[1], -0.5, 0.5 );
		if (n != 0) {
			w2[1] -= n;
			changed = true;
		}

		n = PeriodicImage1D( w2[2]-w1[2], -0.5, 0.5 );
		if (n != 0) {
			w2[2] -= n;
			changed = true;
		}

/*		while (w1[0]-w2[0] > 0.5)
		{
			w2[0] += 1.0;
			changed = true;
		}
		while (w1[0]-w2[0] <= -0.5)
		{
			w2[0] -= 1.0;
			changed = true;
		}
		while (w1[1]-w2[1] > 0.5)
		{
			w2[1] += 1.0;
			changed = true;
		}
		while (w1[1]-w2[1] <= -0.5)
		{
			w2[1] -= 1.0;
			changed = true;
		}
		while (w1[2]-w2[2] > 0.5)
		{
			w2[2] += 1.0;
			changed = true;
		}
		while (w1[2]-w2[2] <= -0.5)
		{
			w2[2] -= 1.0;
			changed = true;
		}*/

		if (changed)
			m_vaCoords[i2] = g_mBoxFromOrtho * w2;

	} else {

		if (g_bPeriodicX) {

			n = PeriodicImage1D( m_vaCoords[i2][0]-m_vaCoords[i1][0], -g_fBoxX/2.0, g_fBoxX/2.0 );
			if (n != 0) {
				m_vaCoords[i2][0] -= n*g_fBoxX;
				changed = true;
			}
		}

		if (g_bPeriodicY) {

			n = PeriodicImage1D( m_vaCoords[i2][1]-m_vaCoords[i1][1], -g_fBoxY/2.0, g_fBoxY/2.0 );
			if (n != 0) {
				m_vaCoords[i2][1] -= n*g_fBoxY;
				changed = true;
			}
		}

		if (g_bPeriodicZ) {

			n = PeriodicImage1D( m_vaCoords[i2][2]-m_vaCoords[i1][2], -g_fBoxZ/2.0, g_fBoxZ/2.0 );
			if (n != 0) {
				m_vaCoords[i2][2] -= n*g_fBoxZ;
				changed = true;
			}
		}

/*		if (g_bPeriodicX)
		{
			while (m_vaCoords[i1][0]-m_vaCoords[i2][0] > g_fBoxX/2)
			{
				m_vaCoords[i2][0] += g_fBoxX;
				changed = true;
			}
			while (m_vaCoords[i2][0]-m_vaCoords[i1][0] > g_fBoxX/2)
			{
				m_vaCoords[i2][0] -= g_fBoxX;
				changed = true;
			}
		}

		if (g_bPeriodicY)
		{
			while (m_vaCoords[i1][1]-m_vaCoords[i2][1] > g_fBoxY/2)
			{
				m_vaCoords[i2][1] += g_fBoxY;
				changed = true;
			}
			while (m_vaCoords[i2][1]-m_vaCoords[i1][1] > g_fBoxY/2)
			{
				m_vaCoords[i2][1] -= g_fBoxY;
				changed = true;
			}
		}

		if (g_bPeriodicZ)
		{
			while (m_vaCoords[i1][2]-m_vaCoords[i2][2] > g_fBoxZ/2)
			{
				m_vaCoords[i2][2] += g_fBoxZ;
				changed = true;
			}
			while (m_vaCoords[i2][2]-m_vaCoords[i1][2] > g_fBoxZ/2)
			{
				m_vaCoords[i2][2] -= g_fBoxZ;
				changed = true; 
			}
		}*/
	}

	if (changed)
	{
		BTOUT; 
		return true;
	} 
	BTOUT; 
	return false;
}


void CTimeStep::RECURSION_ScanMolecules(int i, CxByteArray *ta, CSingleMolecule *sm, int depth, int *stack, unsigned long bmask, bool w)
{
	BTIN;
	int z, z2, z3, z0;
	int nblist[64], nbs;
	double f=0;
	bool tb;
	CxIntArray *wa;	
	CMolAtom *ma;
	std::vector<int> vnb;


//	mprintf("  Rekursion fuer Atom %d.\n",i+1);

//	mprintf("Depth=%d, i=%d.\n",depth,i);
	if (g_bVerbose)
	{
		mprintf("  ");
		for (z=1;z<depth;z++)
		{
			if ((bmask & (int)mypow(2.0,z)) != 0)
				mprintf(WHITE,"|  ");
					else mprintf("   ");
		}
		if (depth != 0)
		{
			if (w)
				mprintf(WHITE,"|--");
					else mprintf(WHITE,"`--");
		}
		mprintf(CYAN,"%s",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[i]])->m_sName);
		mprintf("(%d)",i+1);
	}

	stack[depth] = i;
	g_laAtomSMIndex[i] = (unsigned long)g_oaSingleMolecules.GetSize()-1;
	(*ta)[i] = 1;
//	mprintf("  sm->m_baAtomIndex.GetSize()=%d\n",sm->m_baAtomIndex.GetSize());
	for (z=0;z<sm->m_baAtomIndex.GetSize();z++)
	{
		if (sm->m_baAtomIndex[z] == g_baAtomIndex[i])
		{
			try { ma = new CMolAtom(); } catch(...) { ma = NULL; }
			if (ma == NULL) NewException((double)sizeof(CMolAtom),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			g_oaMolAtoms[i] = ma;

			ma->m_iOffset = i;
			ma->m_iType = z;
			ma->m_iNumber = -1;
			sm->m_oaMolAtoms.Add(ma);
			goto _ok;
		}
	}

	try { ma = new CMolAtom(); } catch(...) { ma = NULL; }
	if (ma == NULL) NewException((double)sizeof(CMolAtom),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	g_oaMolAtoms[i] = ma;
	
	ma->m_iOffset = i;
	ma->m_iType = sm->m_baAtomIndex.GetSize();
	ma->m_iNumber = -1;
	sm->m_oaMolAtoms.Add(ma);

	sm->m_baAtomIndex.Add(g_baAtomIndex[i]);
_ok:
	m_iConnectedAtoms++;

	if (!g_bVerbose)
		if ((m_iConnectedAtoms % m_iGesAtomModulo)==0)
			mprintf(WHITE,"#");

	nbs = 0;

//	for (z=0;z<(long)m_iGesAtomCount;z++) // Schon mal alle Nachbarn raussuchen

	m_pPosDomainEngine->FindNeighbors(i,vnb);

	for (z0=0;z0<(int)vnb.size();z0++)
	{
		z = vnb[z0];
		if (z == i)
			continue;
		if (((CAtom*)g_oaAtoms[g_baAtomIndex[z]])->m_bExclude)
			continue;
		tb = false;
		for (z3=0;z3<g_laBondAddList.GetSize()/2;z3++) {
			if (((g_laBondAddList[z3*2] == i) && (g_laBondAddList[z3*2+1] == z)) ||
				((g_laBondAddList[z3*2] == z) && (g_laBondAddList[z3*2+1] == i))) {
				tb = true;
				break;
			}
		}
		if (BondRange(i,z,&f) || tb)
		{
			if ((*ta)[z] != 0)
			{
				if ((depth > 0) && (z != stack[depth-1]))
				{
					if (sm->m_oaRings.GetSize() < 100)
					{
						if (g_bVerbose)
						{
							mprintf(GREEN,"  <-- Ring closure: ");
							mprintf("%s(%d)",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[stack[depth]]])->m_sName,stack[depth]+1);
						}

						try { wa = new CxIntArray("CTimeStep::RECURSION_ScanMolecules():wa"); } catch(...) { wa = NULL; }
						if (wa == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						sm->m_oaRings.Add(wa);
						wa->Add(stack[depth]);
						z2 = depth-1;
						while ((stack[z2] != z) && (z2 >= 0))
						{
							wa->Add(stack[z2]);
							if (g_bVerbose)
								mprintf(" - %s(%d)",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[stack[z2]]])->m_sName,stack[z2]+1);
							z2--;
						}
						wa->Add(z);
						if (g_bVerbose)
							mprintf(" - %s(%d)",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[z]])->m_sName,z+1);
					} else if (!m_bAbortRing)
					{
						m_bAbortRing = true;
						eprintf("\n### More than 100 rings, aborting ring system scan for molecule %d.",sm->m_iIndex+1);
					}
				}
				continue;
			}
/*			for (z2=0;z2<g_waBondBlackList.GetSize();z2+=2)
				if (((g_waBondBlackList[z2] == i) && (g_waBondBlackList[z2+1] == z)) || 
					((g_waBondBlackList[z2] == z) && (g_waBondBlackList[z2+1] == i)))
				{
					if (g_bVerbose)
					{
						mprintf("\n  ");
						for (z2=1;z2<depth+1;z2++)
						{
							if ((bmask & (int)pow(2,z2)) != 0)
								mprintf(WHITE,"|  ");
							else mprintf("   ");
						}
						if (z+1 < nbs)
							mprintf(WHITE,"|--");
								else mprintf("`--");
						mprintf("%s(%d)",((CAtom*)g_oaAtoms[g_baAtomIndex[z]])->m_sName,z+1);
						mprintf(GREEN,"  <-- This bond has been broken, skipping.\n");
					}
					g_iBondBlackListUsed++;
					goto _nextnb;
				}*/
			if (f < 50.0)
			{
				if (g_iCloseAtomCounter < 25)
					eprintf("\n### The atoms %s(%d) and %s(%d) are VERY close to each other. (Distance=%.4f pm)",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[i]])->m_sName,i+1,(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[z]])->m_sName,z+1,f);
				else if (g_iCloseAtomCounter == 25)
					eprintf("\n### Suppressing further \"atoms close to each other\" warnings.");

				g_iCloseAtomCounter++;
			}
			nblist[nbs] = z;
			nbs++;
			if (nbs >= MAX_BONDS)
			{
				eprintf("\n### Atom %s(%d) has more than %d bonds. Ignoring further bonds (defined via MAX_BONDS in config.h).",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[i]])->m_sName,i+1,MAX_BONDS);
				goto _nbdone;
			}
		}
	}
_nbdone:
	if (g_bVerbose)
		mprintf("\n");
	for (z=0;z<nbs;z++) // Fuer das aktuelle Atom z2 alle Nachbarn durchgehen
	{
		for (z2=0;z2<g_laBondBlackList.GetSize();z2+=2)
		{
			if (((g_laBondBlackList[z2] == i) && (g_laBondBlackList[z2+1] == nblist[z])) || 
				((g_laBondBlackList[z2] == nblist[z]) && (g_laBondBlackList[z2+1] == i)))
			{
				if (g_bVerbose)
				{
					mprintf("  ");
					for (z2=1;z2<depth+1;z2++)
					{
						if ((bmask & (int)mypow(2.0,z2)) != 0)
							mprintf(WHITE,"|  ");
						else mprintf("   ");
					}
					if (z+1 < nbs)
						mprintf(WHITE,"|--");
							else mprintf(WHITE,"`--");
					mprintf(CYAN,"%s",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[nblist[z]]])->m_sName);
					mprintf("(%d)",nblist[z]+1);
					mprintf(GREEN,"  <-- This bond shall be broken, skipping.\n");
				}
				g_iBondBlackListUsed++;
				goto _nextnb;
			}
		}

		sm->m_laBonds.Add(i);
		sm->m_laBonds.Add(nblist[z]);

		if ((*ta)[nblist[z]] == 0) // Der Nachbar ist noch immer frei
		{
			if (z+1 == nbs)
				bmask -= (int)mypow(2.0,depth+1);
			RECURSION_ScanMolecules(nblist[z],ta,sm,depth+1,stack,bmask,(z+1==nbs)?false:true);
		} else // Ringschluss
		{
			if (g_bVerbose)
			{
				mprintf("  ");
				for (z2=1;z2<depth+1;z2++)
				{
					if ((bmask & (int)mypow(2.0,z2)) != 0)
						mprintf(WHITE,"|  ");
					else mprintf("   ");
				}
				if (z+1 < nbs)
					mprintf(WHITE,"|--");
						else mprintf(WHITE,"`--");
				mprintf(CYAN,"%s",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[nblist[z]]])->m_sName);
				mprintf("(%d)",nblist[z]+1);
				mprintf(GREEN,"  <-- Ring closure\n");
			}
		}
_nextnb:;
	}
	BTOUT; 
}


void CTimeStep::BuildAtomIndex()
{
	BTIN;
	int z, z2;
	CAtom *at;
//	char buf[64];
	CxString buf;

	g_baAtomIndex.SetSize(m_iGesAtomCount);
	for (z=0;z<(long)m_iGesAtomCount;z++)
	{
//		strcpy(buf,(char*)m_paLabels[z]);
		buf.strcpy((char*)m_paLabels[z]);
		ReplaceDigits(&buf);
		for (z2=0;z2<g_oaAtoms.GetSize();z2++)
		{
			if (mystricmp(buf,((CAtom*)g_oaAtoms[z2])->m_sName)==0)
			{
				at = (CAtom*)g_oaAtoms[z2];
	//			mprintf("Atom %d: Type %d (%s).\n",z,z2,at->m_sName);
				while (at->m_pMergedTo != NULL)
				{
					at = at->m_pMergedTo;
	//				mprintf("  ...was merged to %s.\n",at->m_sName);
				}
	//			mprintf("  Index is %d.\n",at->m_iIndex);
				g_baAtomIndex[z] = (unsigned char)at->m_iIndex;
				goto _e;
			}
		}
		g_baAtomIndex[z] = 255;
		eprintf("Error: CTimeStep::BuildAtomIndex(): Atom type \"%s\" not known.\n",(const char*)buf);
		_e:;
	}
	BTOUT; 
}


bool CTimeStep::ScanMolecules()
{
	BTIN;
	int z, z2, z3, z4, z5, i, i2, ti, ti2, i1;
	CxByteArray ta;
	CSingleMolecule *sm, *sm2;
	CxIntArray *wa, *waz2, *waz3, *waneu;
	CMolecule *m;
	CMolAtom *ma, *ma2, *ma3;
	CMolBond *bond, *bond2;
	CMolBondGroup *bg;
	CMolAngle *angle, *angle2;
	CMolAngleGroup *ag;
	CAtom *at;
	double f, tfs;
	bool b, topo;
	int *stack;

	BuildAtomIndex();
//	CalcRadii();

	try { stack = new int[g_iGesAtomCount]; } catch(...) { stack = NULL; }
	if (stack == NULL) NewException((double)g_iGesAtomCount*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	b = false;
	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		if (z == g_iVirtAtomType)
			continue;
		at = (CAtom*)g_oaAtoms[z];
		if (at->m_pMergedTo != NULL)
			continue;
		if (at->m_pElement->m_fRadius == 0)
		{
			if (!b)
			{
				b = true;
				mprintf("\n");
			}
			//at->m_bExclude = AskYesNo("    Atom type \"%s\" has bond radius 0. Exclude it from the system (y/n)? [yes] ",true,(const char*)at->m_sName);
//			at->m_bExclude = false;
	//		if (!at->m_bExclude)
	//			at->m_pElement->m_fRadius = AskFloat_ND("    Please enter covalent bond radius for %s (in pm): ",at->m_sName);
		}
	}

	ta.SetSize(g_iGesAtomCount);
	g_laAtomSMIndex.SetSize(g_iGesAtomCount);
	g_oaMolAtoms.SetSize(g_iGesAtomCount);
	for (z=0;z<g_iGesAtomCount;z++)
	{
		ta[z] = 0;
		g_laAtomSMIndex[z] = -1;
		g_oaMolAtoms[z] = NULL;
	}
	g_iBondBlackListUsed = 0;
	g_oaSingleMolecules.RemoveAll();
	g_oaMolecules.RemoveAll();
	m_iConnectedAtoms = 0;
	m_iGesAtomModulo = (m_iGesAtomCount / 60)+1;
	if (g_iScanMolStep == -1)
	{
		mprintf(WHITE,"\n    Skipping molecule recognition, using atoms as molecules...\n\n");
		for (z=0;z<(long)m_iGesAtomCount;z++)
		{
			if (((CAtom*)g_oaAtoms[g_baAtomIndex[z]])->m_bExclude)
				continue;
			
			try { sm = new CSingleMolecule(); } catch(...) { sm = NULL; }
			if (sm == NULL) NewException((double)sizeof(CSingleMolecule),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			sm->m_iIndex = g_oaSingleMolecules.GetSize();

			g_oaSingleMolecules.Add(sm);

			g_laAtomSMIndex[z] = (unsigned long)g_oaSingleMolecules.GetSize()-1;

			try { ma = new CMolAtom(); } catch(...) { ma = NULL; }
			if (ma == NULL) NewException((double)sizeof(CMolAtom),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			ma->m_iOffset = z;
			ma->m_iType = sm->m_baAtomIndex.GetSize();
			ma->m_iNumber = -1;
			sm->m_oaMolAtoms.Add(ma);
			sm->m_baAtomIndex.Add(g_baAtomIndex[z]);
		}
		mprintf("    %d molecules found.\n\n",g_oaSingleMolecules.GetSize());
		topo = false;
	} else
	{
	//	fflush(stdout);

		m_pPosDomainEngine = new CPosDomainEngine();

		f = 0;
		for (z=0;z<g_oaAtoms.GetSize()-1;z++)
			if (((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius > f)
				f = ((CAtom*)g_oaAtoms[z])->m_pElement->m_fRadius;
		f *= 2 * g_fBondFactor;

		if (g_bPeriodicX && g_bPeriodicY && g_bPeriodicZ && !g_bMolRecNoDomain)
			m_pPosDomainEngine->Create( m_vaCoords, g_mBoxFromOrtho, f );
		else
			m_pPosDomainEngine->CreateTrivial( m_vaCoords );

		if (g_bAdvanced1 && !m_pPosDomainEngine->m_bTrivial) {

			if (!AskYesNo("    Use spatial domain decomposition for faster mol. rec. (disable if there are problems) (y/n)? [yes] ",true)) {

				delete m_pPosDomainEngine;
				m_pPosDomainEngine = new CPosDomainEngine();
				m_pPosDomainEngine->CreateTrivial( m_vaCoords );
			}
		}

		if (g_bAdvanced1)
			topo = AskYesNo("    Use refined topological atom ordering (can take very long for large molecules) (y/n)? [yes] ",true);
		else
			topo = true;

		if (m_pPosDomainEngine->m_bTrivial)
			mprintf("\n    Not using spatial domain engine (min. diameter = %.0f pm).\n",f);
		else
			mprintf("\n    Initialized spatial domain engine with %d x %d x %d cells, min. diameter = %.0f pm.\n",m_pPosDomainEngine->m_iRes[0],m_pPosDomainEngine->m_iRes[1],m_pPosDomainEngine->m_iRes[2],f);


		if (g_bVerbose)
		{
			mprintf(WHITE,"\n    Molecule recognition...\n\n");
			mprintf(WHITE,">>> Output of the molecule tree >>>\n");
		} else
			mprintf("\n    Molecule recognition [");

		for (z=0;z<(long)m_iGesAtomCount;z++)
		{
			if (ta[z] != 0) // Dieses Atom wurde bereits in irgendein Molekuel eingebaut
			{
	//			mprintf("> Atom %d: Schon vergeben.\n",z+1);
				continue;
			}
	//		mprintf("# Atom %d: Starte Rekursion.\n",z+1);

			if (((CAtom*)g_oaAtoms[g_baAtomIndex[z]])->m_bExclude)
				continue;
			
			try { sm = new CSingleMolecule(); } catch(...) { sm = NULL; }
			if (sm == NULL) NewException((double)sizeof(CSingleMolecule),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			sm->m_iIndex = g_oaSingleMolecules.GetSize();

			g_oaSingleMolecules.Add(sm);

	//		sm->m_iAtomGes = 0;
	/*		for (z2=0;z2<16;z2++)
				sm->m_iAtomCount[z2] = 0;
			sm->m_iElements = 0;*/
			if (g_bVerbose)
				mprintf(YELLOW,"\nThe next molecule starts with %s(%d):\n",(const char*)m_paLabels[z],z+1);
			
			m_bAbortRing = false;

			RECURSION_ScanMolecules(z,&ta,sm,0,stack,0xFFFFFFFF,true);

	//		printf("%d Atome in diesem Molekuel.\n",g_pSingleMolecules[g_oaSingleMolecules.GetSize()].AtomGes);

	//		g_oaSingleMolecules.GetSize()++;
		}

		delete m_pPosDomainEngine;

		if (g_bVerbose)
		{
			mprintf(WHITE,"\n<<< Output of the molecule tree <<<\n\n");
			mprintf("%d molecules found.\n\n",g_oaSingleMolecules.GetSize());
		} else mprintf("]\n\n    %d molecules found.\n\n",g_oaSingleMolecules.GetSize());
	}


	if (g_iCloseAtomCounter != 0)
	{
		mprintf(RED,"\n*** Warning: ");
		mprintf("Some of the atoms were found to be very close to each other (distance < 50 pm).\n");
		mprintf("             This usually indicates a problem, e.g. an incorrect cell vector.\n\n");
		if (!AskYesNo("    Continue with molecule recognition anyways (y/n)? [yes] ",true))
			return false;
		mprintf("\n");
	}

	mprintf("    Sorting atom types...\n");
	SortSingleMolAtomTypes();

	delete[] stack;

/*	mprintf("\n");
	for (z=0;z<g_oaSingleMolecules.GetSize();z++)
	{
		mprintf("%d ",z+1);
		((CSingleMolecule*)g_oaSingleMolecules[z])->Dump();
	}

	mprintf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");*/

/*	mprintf("\n");
	for (z=0;z<g_oaSingleMolecules.GetSize();z++)
	{
		mprintf("%d ",z+1);
		((CSingleMolecule*)g_oaSingleMolecules[z])->Dump();
	}
	mprintf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");*/
	if (g_iBondBlackListUsed != 0)
		mprintf("    %d bonds have been broken.\n",g_iBondBlackListUsed);
	
/*	if (g_oaSingleMolecules.GetSize() >= 4096)
	{
		mprintf("\n>>> Mehr als 4096 Molekuele! Dies wird nicht unterstuetzt. <<<\n\n");
//		g_oaSingleMolecules.GetSize() = 4096;
	}*/

	mprintf("    Setting up bond lists...\n");
	// Die Bindungen in m_oaMolAtoms aufbauen
	for (z=0;z<g_oaSingleMolecules.GetSize();z++) 
	{
		sm = (CSingleMolecule*)g_oaSingleMolecules[z];
		for (z2=0;z2<sm->m_laBonds.GetSize()/2;z2++)
		{
			i = sm->m_laBonds[z2*2];
			i2 = sm->m_laBonds[z2*2+1];
			ti = -1;
			for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
			{
				if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_iOffset == i)
				{
					ti = z3;
					break;
				}
			}
			if (ti == -1)
			{
				eprintf("CTimeStep::ScanMolecules(): Atom 1 of bond (%d) not found.\n",i);
				return false;
			}
			ti2 = -1;
			for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
			{
				if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_iOffset == i2)
				{
					ti2 = z3;
					break;
				}
			}
			if (ti2 == -1)
			{
				eprintf("CTimeStep::ScanMolecules(): Atom 2 of bond (%d) not found.\n",i2);
				return false;
			}
			((CMolAtom*)sm->m_oaMolAtoms[ti])->m_oaBonds.Add((CMolAtom*)sm->m_oaMolAtoms[ti2]);
			((CMolAtom*)sm->m_oaMolAtoms[ti2])->m_oaBonds.Add((CMolAtom*)sm->m_oaMolAtoms[ti]);
		}
	}

	// Die AtomCodes berechnen
	mprintf("    Building atom codes...\n");
	for (z=0;z<g_oaSingleMolecules.GetSize();z++) 
	{
//		mprintf("Molecule %d: ",z+1);
		sm = (CSingleMolecule*)g_oaSingleMolecules[z];
		if (g_bVerbose)
			mprintf(YELLOW,"\n*** Singlemolecule %d ***\n\n",z+1);
		sm->BuildAtomCodes();
	}

	if (topo) {
		mprintf("    Creating refined topological atom order...\n");
		mprintf(WHITE,"      [");
		tfs = g_oaSingleMolecules.GetSize() / 60.0;
		for (z=0;z<g_oaSingleMolecules.GetSize();z++) 
		{
			if (fmod(z,tfs) < 1.0)
				mprintf(WHITE,"#");
			sm = (CSingleMolecule*)g_oaSingleMolecules[z];
			sm->BuildAtomOrder(false);
	/*		for (z2=0;z2<sm->m_baAtomIndex.GetSize();z2++)
			{
				try { wa = new CxIntArray("CTimeStep::ScanMolecules():wa"); } catch(...) { wa = NULL; }
				if (wa == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				sm->m_oaAtomOffset.Add(wa);
				for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
				{
					ma = (CMolAtom*)sm->m_oaMolAtoms[z3];
					if (ma->m_iType != z2)
						continue;
					ma->m_iNumber = wa->GetSize();
					wa->Add(ma->m_iOffset);
				}
			}*/
		}
		mprintf(WHITE,"] Done.\n");
	} else {
		mprintf("    Creating non-refined topological atom order...\n");
		for (z=0;z<g_oaSingleMolecules.GetSize();z++) {
			sm = (CSingleMolecule*)g_oaSingleMolecules[z];
			sm->BuildAtomOrder(true);
		}
	}

	mprintf("    Creating bond list...\n");
	for (z=0;z<g_oaSingleMolecules.GetSize();z++) 
	{
		sm = (CSingleMolecule*)g_oaSingleMolecules[z];
//		mprintf("%d:\n",z+1);
		for (z2=0;z2<sm->m_oaMolAtoms.GetSize();z2++)
		{
			ma = (CMolAtom*)sm->m_oaMolAtoms[z2];
			for (z3=0;z3<ma->m_oaBonds.GetSize();z3++)
			{
				ma2 = (CMolAtom*)ma->m_oaBonds[z3];
				if (ma2->m_iMolAtomNumber < z2 /* == ma->m_iMolAtomNumber*/)
					continue;

		//		mprintf("  %s%2d <--> %s%2d\n",((CAtom*)g_oaAtoms[sm->m_baAtomIndex[ma->m_iType]])->m_sName,ma->m_iNumber+1,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[ma2->m_iType]])->m_sName,ma2->m_iNumber+1);

				try { bond = new CMolBond(); } catch(...) { bond = NULL; }
				if (bond == NULL) NewException((double)sizeof(CMolBond),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				bond->m_iAtom[0] = ma->m_iNumber;
				bond->m_iAtom[1] = ma2->m_iNumber;
				//mprintf("@ %d-%d\n",ma->m_iNumber,ma2->m_iNumber);
				bond->m_iAtomType[0] = ma->m_iType;
				bond->m_iAtomType[1] = ma2->m_iType;
				bond->m_iMolAtom[0] = z2;
				bond->m_iMolAtom[1] = ma2->m_iMolAtomNumber;
				bond->m_iAtomOffset[0] = ma->m_iOffset;
				bond->m_iAtomOffset[1] = ma2->m_iOffset;
				sm->m_oaBonds.Add(bond);
			}
		}
/*		if (z == 0)
			mprintf("  %d Bonds added.\n",oatemp.GetSize());*/
		for (z2=0;z2<sm->m_oaBonds.GetSize();z2++)
		{
			bond = (CMolBond*)sm->m_oaBonds[z2];
			for (z3=0;z3<sm->m_oaBondGroups.GetSize();z3++)
			{
				bg = (CMolBondGroup*)sm->m_oaBondGroups[z3];
				for (z4=0;z4<bg->m_oaBonds.GetSize();z4++)
				{
					bond2 = (CMolBond*)bg->m_oaBonds[z4];
/*					if (((((CMolAtom*)sm->m_oaMolAtoms[bond->m_iMolAtom[0]])->m_fAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[bond2->m_iMolAtom[0]])->m_fAtomCode) &&
						 (((CMolAtom*)sm->m_oaMolAtoms[bond->m_iMolAtom[1]])->m_fAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[bond2->m_iMolAtom[1]])->m_fAtomCode)) ||
						((((CMolAtom*)sm->m_oaMolAtoms[bond->m_iMolAtom[0]])->m_fAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[bond2->m_iMolAtom[1]])->m_fAtomCode) &&
						 (((CMolAtom*)sm->m_oaMolAtoms[bond->m_iMolAtom[1]])->m_fAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[bond2->m_iMolAtom[0]])->m_fAtomCode)))	
*/					if (((((CMolAtom*)sm->m_oaMolAtoms[bond->m_iMolAtom[0]])->m_liAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[bond2->m_iMolAtom[0]])->m_liAtomCode) &&
						 (((CMolAtom*)sm->m_oaMolAtoms[bond->m_iMolAtom[1]])->m_liAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[bond2->m_iMolAtom[1]])->m_liAtomCode)) ||
						((((CMolAtom*)sm->m_oaMolAtoms[bond->m_iMolAtom[0]])->m_liAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[bond2->m_iMolAtom[1]])->m_liAtomCode) &&
						 (((CMolAtom*)sm->m_oaMolAtoms[bond->m_iMolAtom[1]])->m_liAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[bond2->m_iMolAtom[0]])->m_liAtomCode)))	
					{
						bg->m_oaBonds.Add(bond);
						goto _bonddone;
					}
				}
			}

			try { bg = new CMolBondGroup(); } catch(...) { bg = NULL; }
			if (bg == NULL) NewException((double)sizeof(CMolBondGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			bg->m_oaBonds.Add(bond);
			sm->m_oaBondGroups.Add(bg);
_bonddone:;
		}
/*		if (z == 0)
		{
			mprintf("  %d Bond groups found:\n",sm->m_oaBondGroups.GetSize());
			for (z2=0;z2<sm->m_oaBondGroups.GetSize();z2++)
			{
				mprintf("  - ");
				bg = (CMolBondGroup*)sm->m_oaBondGroups[z2];
				for (z3=0;z3<bg->m_oaBonds.GetSize();z3++)
				{
					bond = (CMolBond*)bg->m_oaBonds[z3];
					mprintf("%s%2d <--> %s%2d",((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[0]]])->m_sName,bond->m_iAtom[0]+1,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[1]]])->m_sName,bond->m_iAtom[1]+1);
					if (z3 < bg->m_oaBonds.GetSize()-1)
						mprintf(", ");
				}
				mprintf("\n");
			}
		}*/
	}

	mprintf("    Creating angle list...\n");
	for (z=0;z<g_oaSingleMolecules.GetSize();z++) 
	{
		sm = (CSingleMolecule*)g_oaSingleMolecules[z];
//		mprintf("%d:\n",z+1);
		for (z2=0;z2<sm->m_oaMolAtoms.GetSize();z2++)
		{
			ma = (CMolAtom*)sm->m_oaMolAtoms[z2];
			for (z3=0;z3<ma->m_oaBonds.GetSize();z3++)
			{
				ma2 = (CMolAtom*)ma->m_oaBonds[z3];
				for (z4=z3+1;z4<ma->m_oaBonds.GetSize();z4++)
				{
					ma3 = (CMolAtom*)ma->m_oaBonds[z4];
			//		mprintf("  (%s%2d, %s%2d, %s%2d)\n",((CAtom*)g_oaAtoms[sm->m_baAtomIndex[ma2->m_iType]])->m_sName,ma2->m_iNumber+1,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[ma->m_iType]])->m_sName,ma->m_iNumber+1,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[ma3->m_iType]])->m_sName,ma3->m_iNumber+1);

					try { angle = new CMolAngle(); } catch(...) { angle = NULL; }
					if (angle == NULL) NewException((double)sizeof(CMolAngle),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					angle->m_iAtom[0] = ma2->m_iNumber;
					angle->m_iAtom[1] = ma->m_iNumber;
					angle->m_iAtom[2] = ma3->m_iNumber;
					angle->m_iAtomType[0] = ma2->m_iType;
					angle->m_iAtomType[1] = ma->m_iType;
					angle->m_iAtomType[2] = ma3->m_iType;
					angle->m_iMolAtom[0] = ma2->m_iMolAtomNumber;
					angle->m_iMolAtom[1] = z2;
					angle->m_iMolAtom[2] = ma3->m_iMolAtomNumber;
					angle->m_iAtomOffset[0] = ma2->m_iOffset;
					angle->m_iAtomOffset[1] = ma->m_iOffset;
					angle->m_iAtomOffset[2] = ma3->m_iOffset;
					sm->m_oaAngles.Add(angle);
				}
			}
		}
/*		if (z ==  0)
			mprintf("  %d Angles added.\n",oatemp.GetSize());*/
		if (sm->m_oaAngles.GetSize() > 1000)
			mprintf(WHITE,"      [");
		tfs = sm->m_oaAngles.GetSize() / 50.0;
		for (z2=0;z2<sm->m_oaAngles.GetSize();z2++)
		{
			if (sm->m_oaAngles.GetSize() > 1000)
				if (fmod(z2,tfs) < 1)
					mprintf(WHITE,"#");
			angle = (CMolAngle*)sm->m_oaAngles[z2];
			for (z3=0;z3<sm->m_oaAngleGroups.GetSize();z3++)
			{
				ag = (CMolAngleGroup*)sm->m_oaAngleGroups[z3];
				for (z4=0;z4<ag->m_oaAngles.GetSize();z4++)
				{
					angle2 = (CMolAngle*)ag->m_oaAngles[z4];
/*					if (( ((CMolAtom*)sm->m_oaMolAtoms[angle->m_iMolAtom[1]])->m_fAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[angle2->m_iMolAtom[1]])->m_fAtomCode) &&
						(((( ((CMolAtom*)sm->m_oaMolAtoms[angle->m_iMolAtom[0]])->m_fAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[angle2->m_iMolAtom[0]])->m_fAtomCode) &&
						( ((CMolAtom*)sm->m_oaMolAtoms[angle->m_iMolAtom[2]])->m_fAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[angle2->m_iMolAtom[2]])->m_fAtomCode))) ||
						(( ((CMolAtom*)sm->m_oaMolAtoms[angle->m_iMolAtom[0]])->m_fAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[angle2->m_iMolAtom[2]])->m_fAtomCode) &&
						( ((CMolAtom*)sm->m_oaMolAtoms[angle->m_iMolAtom[2]])->m_fAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[angle2->m_iMolAtom[0]])->m_fAtomCode))))	
*/					if (( ((CMolAtom*)sm->m_oaMolAtoms[angle->m_iMolAtom[1]])->m_liAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[angle2->m_iMolAtom[1]])->m_liAtomCode) &&
						(((( ((CMolAtom*)sm->m_oaMolAtoms[angle->m_iMolAtom[0]])->m_liAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[angle2->m_iMolAtom[0]])->m_liAtomCode) &&
						( ((CMolAtom*)sm->m_oaMolAtoms[angle->m_iMolAtom[2]])->m_liAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[angle2->m_iMolAtom[2]])->m_liAtomCode))) ||
						(( ((CMolAtom*)sm->m_oaMolAtoms[angle->m_iMolAtom[0]])->m_liAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[angle2->m_iMolAtom[2]])->m_liAtomCode) &&
						( ((CMolAtom*)sm->m_oaMolAtoms[angle->m_iMolAtom[2]])->m_liAtomCode == ((CMolAtom*)sm->m_oaMolAtoms[angle2->m_iMolAtom[0]])->m_liAtomCode))))	
					{
						ag->m_oaAngles.Add(angle);
						goto _angledone;
					}
				}
			}

			try { ag = new CMolAngleGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CMolAngleGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			ag->m_oaAngles.Add(angle);
			sm->m_oaAngleGroups.Add(ag);
_angledone:;
		}
		if (sm->m_oaAngles.GetSize() > 1000)
			mprintf(WHITE,"]\n");
/*		if (z == 0)
		{
			mprintf("  %d Angle groups found:\n",sm->m_oaAngleGroups.GetSize());
			for (z2=0;z2<sm->m_oaAngleGroups.GetSize();z2++)
			{
				mprintf("  - ");
				ag = (CMolAngleGroup*)sm->m_oaAngleGroups[z2];
				for (z3=0;z3<ag->m_oaAngles.GetSize();z3++)
				{
					angle = (CMolAngle*)ag->m_oaAngles[z3];
					mprintf("  (%s%2d, %s%2d, %s%2d)",((CAtom*)g_oaAtoms[sm->m_baAtomIndex[angle->m_iAtomType[0]]])->m_sName,angle->m_iAtom[0]+1,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[angle->m_iAtomType[1]]])->m_sName,angle->m_iAtom[1]+1,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[angle->m_iAtomType[2]]])->m_sName,angle->m_iAtom[2]+1);
					if (z3 < ag->m_oaAngles.GetSize()-1)
						mprintf(", ");
				}
				mprintf("\n");
			}
		}*/
	}

	if (g_iaaMoleculeTypeSplit.size() != 0)
		mprintf("    Grouping together equivalent molecules (%lu type splits)...\n",(unsigned long)g_iaaMoleculeTypeSplit.size());
	else
		mprintf("    Grouping together equivalent molecules...\n");

	for (z=0;z<g_oaSingleMolecules.GetSize();z++) // Alle Molekuelvertreter durchgehen
	{
//		mprintf("Pruefe SM %d...\n",z+1);
		sm = (CSingleMolecule*)g_oaSingleMolecules[z];

		// Split Mechanism
		ti = -1;
		for (z2=0;z2<(int)g_iaaMoleculeTypeSplit.size();z2++)
			for (z3=0;z3<(int)g_iaaMoleculeTypeSplit[z2].size();z3++)
				if (z == g_iaaMoleculeTypeSplit[z2][z3]) {
					ti = z2;
					//mprintf("@SM %d is on list %d.\n",z,ti);
					goto _split;
				}
_split:
		// End Split Mechanism

		for (z2=0;z2<g_oaMolecules.GetSize();z2++) // Alle schon vorhandenen Molekueltemplates durchgehen
		{
			m = (CMolecule*)g_oaMolecules[z2];
			sm2 = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];

			// Split Mechanism
			ti2 = m->m_laSingleMolIndex[0];
			if (ti != -1) { // SingleMolecule z is on split list ti
				for (z3=0;z3<(int)g_iaaMoleculeTypeSplit[ti].size();z3++)
					if (g_iaaMoleculeTypeSplit[ti][z3] == ti2)
						goto _nosplit; // SingleMolecule ti2 is on the same list --> may merge
				//mprintf("  @SM %d is not on that list --> skipping MT %d.\n",ti2,z2);
				goto _cont; // ti2 is not on that list --> may not merge
_nosplit:;
				//mprintf("  @SM %d is also on that list --> Taking MT %d.\n",ti2,z2);
			} else { // SingleMolecule z is NOT on any split list
				for (z3=0;z3<(int)g_iaaMoleculeTypeSplit.size();z3++)
					for (z4=0;z4<(int)g_iaaMoleculeTypeSplit[z3].size();z4++)
						if (ti2 == g_iaaMoleculeTypeSplit[z3][z4])
							goto _cont; // But ti2 is on some list --> may not merge
			}
			// End Split Mechanism

/*			if (m->m_baAtomIndex.GetSize() != sm->m_baAtomIndex.GetSize())
				goto _cont; // Nich die Gleiche Anzahl verschiedener Elemente -> Keine Uebereinstimmung
			for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
			{
				if (m->m_baAtomIndex[z3] != sm->m_baAtomIndex[z3])
					goto _cont; // Elemente an Position z3 unterscheiden sich -> keine Uebereinstimmung
				if (m->m_waAtomCount[z3] != ((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize())
					goto _cont; // Anzahl an Atomen dieser Elementsorte unterscheidet sich -> Keine Uebereinstimmung*/
//				mprintf("  Stimmt mit M%d ueberein.\n",z2);
//				((CMolecule*)g_oaMolecules[z2])->Count++;
//				g_pSingleMolecules[z].m_iMoleculeOffset = z2;
//			}
			if (sm->m_oaMolAtoms.GetSize() != sm2->m_oaMolAtoms.GetSize())
				goto _cont;
			for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
//				if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_fAtomCode != ((CMolAtom*)sm2->m_oaMolAtoms[z3])->m_fAtomCode)
				if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_liAtomCode != ((CMolAtom*)sm2->m_oaMolAtoms[z3])->m_liAtomCode)
					goto _cont;
			sm->m_iMolSMIndex = m->m_laSingleMolIndex.GetSize();
			m->m_laSingleMolIndex.Add(z);
//			sm->m_iMolType = z2;
			goto _fert;
			_cont:;
		}

		try { m = new CMolecule(); } catch(...) { m = NULL; }
		if (m == NULL) NewException((double)sizeof(CMolecule),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		g_oaMolecules.Add(m);
//		m->Count = 1; // Neues Molekueltemplate hinzufuegen
//		m->Elements = ((CSingleMolecule*)g_oaSingleMolecules[z])->m_iElements;
		m->m_iAtomGes = 0;
		f = 0;
		for (z2=0;z2<sm->m_oaAtomOffset.GetSize();z2++)
		{
			m->m_iAtomGes += ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetSize();
			f += ((CAtom*)g_oaAtoms[sm->m_baAtomIndex[z2]])->m_pElement->m_fMass * ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetSize();
		}
		m->m_iAtomGesNoVirt = m->m_iAtomGes;
		m->m_fMass = f;
		((CSingleMolecule*)g_oaSingleMolecules[z])->m_iMolType = g_oaMolecules.GetSize()-1;
		for (z2=0;z2<((CSingleMolecule*)g_oaSingleMolecules[z])->m_baAtomIndex.GetSize();z2++)
		{
			m->m_baAtomIndex.Add(((CSingleMolecule*)g_oaSingleMolecules[z])->m_baAtomIndex[z2]);
			m->m_waAtomCount.Add((unsigned short)((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[z])->m_oaAtomOffset[z2])->GetSize());
//			m->AtomCount[z2] = ((CSingleMolecule*)g_oaSingleMolecules[z])->m_iAtomCount[z2];
		}
		sm->m_iMolSMIndex = m->m_laSingleMolIndex.GetSize();
		m->m_laSingleMolIndex.Add(z);
//		g_pSingleMolecules[z].m_iMoleculeOffset = g_oaMolecules.GetSize();
		m->BuildName();

		// Standardmaessig geom. Zentrum als Dipolzentrum
		m->m_iDipoleCenterType = m->m_baAtomIndex.GetSize();
		m->m_iDipoleCenterIndex = 0;

//		mprintf("### Erzeuge neues Molekuel %d \"%s\".\n",g_oaMolecules.GetSize(),m->Name);
//		g_oaMolecules.GetSize()++;
		_fert:;
	}
	mprintf("    Found %d unique molecule types.\n",g_oaMolecules.GetSize());

	// Sort by molecular mass
	mprintf("    Sorting molecule types by molecular mass...\n");
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
//		mprintf("@ z=%d\n",z);
		f = -1.0;
		z3 = -1;
		for (z2=z;z2<g_oaMolecules.GetSize();z2++)
		{
//			mprintf("@   z2=%d\n",z2);
			if (((CMolecule*)g_oaMolecules[z2])->m_fMass > f)
			{
//				mprintf("@     %f > %f --> z3 = %d\n",((CMolecule*)g_oaMolecules[z2])->m_fMass,f,z2);
				f = ((CMolecule*)g_oaMolecules[z2])->m_fMass;
				z3 = z2;
			}
		}
		if (z3 == -1) {
			eprintf("CTimeStep::ScanMolecules(): Internal error.\n");
			abort();
		}
		if (f > 0)
		{
//			mprintf("@ Swapping %d with %d.\n",z,z3);
			m = (CMolecule*)g_oaMolecules[z3];
			g_oaMolecules[z3] = g_oaMolecules[z];
			g_oaMolecules[z] = m;
		}
	}

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		((CMolecule*)g_oaMolecules[z])->m_iIndex = z;
		for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
			((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex[z2]])->m_iMolType = z;
	}


/*	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex[0]];
		mprintf("Molecule %d Ring systems:\n",z+1);
		for (z2=0;z2<sm->m_oaRings.GetSize();z2++)
		{
			wa = (CxIntArray*)sm->m_oaRings[z2];
			mprintf("  - Ring %d: ",z2+1);
			for (z3=0;z3<wa->GetSize();z3++)
			{
				mprintf("%d",(*wa)[z3]);
				if (z3+1 < wa->GetSize())
					mprintf(", ");
			}
			mprintf("\n");
		}
		mprintf("\n");
	}*/

	ti = 0;
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex[0]];
		ti += sm->m_oaRings.GetSize();
	}

	mprintf("    Found %d rings.\n",ti);
	if (ti != 0)
	{
		mprintf("    Refining ring systems...\n");

		try { wa = new CxIntArray("CTimeStep::ScanMolecules():wa"); } catch(...) { wa = NULL; }
		if (wa == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		// Process ring systems
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			m = (CMolecule*)g_oaMolecules[z];
			sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex[0]];
			if (sm->m_oaRings.GetSize() > 100)
			{
				eprintf("      More than 100 rings in molecule %d, skipping refinement. Ring system may be improper.\n",z+1);
				continue;
			}
			ti2 = 0;
_again:
			ti2++;
			if (ti2 > 1000)
			{
				eprintf("      Too many iterations in molecule %d, aborting refinement. Ring system may be improper.\n",z+1);
				continue;
			}
//			mprintf("### again\n");
			for (z2=0;z2<sm->m_oaRings.GetSize();z2++)
			{
				waz2 = (CxIntArray*)sm->m_oaRings[z2];
				for (z3=z2+1;z3<sm->m_oaRings.GetSize();z3++)
				{
					waz3 = (CxIntArray*)sm->m_oaRings[z3];

//					mprintf("*A* Checking if %d(%d) contains %d(%d)...\n",z2+1,waz2->GetSize(),z3+1,waz3->GetSize());
					// Check if waz2 contains waz3
					wa->SetSize(waz2->GetSize());
					for (z4=0;z4<waz2->GetSize();z4++)
						(*wa)[z4] = 0;
					for (z4=0;z4<waz3->GetSize();z4++)
					{
						for (z5=0;z5<waz2->GetSize();z5++)
						{
							if ((*waz3)[z4] == (*waz2)[z5])
							{
								(*wa)[z5] = 1;
								goto _found2;
							}
						}
						goto _notfound2;
_found2:;
					}
//					mprintf("  Yes!\n");

					try { waneu = new CxIntArray("CTimeStep::ScanMolecules():waneu"); } catch(...) { waneu = NULL; }
					if (waneu == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					for (z4=0;z4<waz2->GetSize();z4++)
					{
						if (z4 == 0)
						{
							if (((*wa)[z4] == 0) || ((*wa)[z4+1] == 0) || ((*wa)[waz2->GetSize()-1] == 0))
								waneu->Add((*waz2)[z4]);
						} else if (z4 == waz2->GetSize()-1)
						{
							if (((*wa)[z4] == 0) || ((*wa)[z4-1] == 0) || ((*wa)[0] == 0))
								waneu->Add((*waz2)[z4]);
						} else
						{
							if (((*wa)[z4] == 0) || ((*wa)[z4+1] == 0) || ((*wa)[z4-1] == 0))
								waneu->Add((*waz2)[z4]);
						}
					}
//					mprintf("   Took %d from %d atoms from waz2.\n",waneu->GetSize(),waz2->GetSize());
					i1 = waz2->GetSize();
					delete waz2;
					sm->m_oaRings[z2] = waneu;
					if (waneu->GetSize() != i1)
						goto _again;
_notfound2:

//					mprintf("*B* Checking if %d contains %d...\n",z3+1,z2+1);
					// Check if waz3 contains waz2
					wa->SetSize(waz3->GetSize());
					for (z4=0;z4<waz3->GetSize();z4++)
						(*wa)[z4] = 0;
					for (z4=0;z4<waz2->GetSize();z4++)
					{
						for (z5=0;z5<waz3->GetSize();z5++)
						{
							if ((*waz2)[z4] == (*waz3)[z5])
							{
								(*wa)[z5] = 1;
								goto _found3;
							}
						}
						goto _notfound3;
_found3:;
					}
//					mprintf("  Yes!\n");

					try { waneu = new CxIntArray("CTimeStep::ScanMolecules():waneu"); } catch(...) { waneu = NULL; }
					if (waneu == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					for (z4=0;z4<waz3->GetSize();z4++)
					{
						if (z4 == 0)
						{
							if (((*wa)[z4] == 0) || ((*wa)[z4+1] == 0) || ((*wa)[waz3->GetSize()-1] == 0))
								waneu->Add((*waz3)[z4]);
						} else if (z4 == waz3->GetSize()-1)
						{
							if (((*wa)[z4] == 0) || ((*wa)[z4-1] == 0) || ((*wa)[0] == 0))
								waneu->Add((*waz3)[z4]);
						} else
						{
							if (((*wa)[z4] == 0) || ((*wa)[z4+1] == 0) || ((*wa)[z4-1] == 0))
								waneu->Add((*waz3)[z4]);
						}
					}
//					mprintf("   Took %d from %d atoms from waz3.\n",waneu->GetSize(),waz3->GetSize());
					i1 = waz3->GetSize();
					delete waz3;
					sm->m_oaRings[z3] = waneu;
					if (waneu->GetSize() != i1)
						goto _again;
_notfound3:;
				}
			}
		}
		delete wa;

		mprintf("    Assigning ring systems to molecule types...\n");
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			m = (CMolecule*)g_oaMolecules[z];
			sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex[0]];
			for (z2=0;z2<sm->m_oaRings.GetSize();z2++)
			{
				wa = (CxIntArray*)sm->m_oaRings[z2];

				try { waz2 = new CxIntArray("CTimeStep::ScanMolecules():waz2"); } catch(...) { waz2 = NULL; }
				if (waz2 == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				try { waz3 = new CxIntArray("CTimeStep::ScanMolecules():waz3"); } catch(...) { waz3 = NULL; }
				if (waz3 == NULL) NewException((double)sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				m->m_oaRingAtomTypes.Add(waz2);
				m->m_oaRingAtoms.Add(waz3);
				for (z3=0;z3<wa->GetSize();z3++)
				{
					for (z4=0;z4<sm->m_baAtomIndex.GetSize();z4++)
					{
						for (z5=0;z5<((CxIntArray*)sm->m_oaAtomOffset[z4])->GetSize();z5++)
						{
							if (((CxIntArray*)sm->m_oaAtomOffset[z4])->GetAt(z5) == wa->GetAt(z3))
							{
								waz2->Add(z4);
								waz3->Add(z5);
								goto _done;
							}
						}
					}
_done:;
				}
			}
		}

		mprintf("    Sorting rings by size...\n");
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			m = (CMolecule*)g_oaMolecules[z];

			for (z2=0;z2<m->m_oaRingAtoms.GetSize();z2++)
			{
				i1 = 0;
				ti = -1;
				for (z3=z2;z3<m->m_oaRingAtoms.GetSize();z3++)
				{
					if (((CxIntArray*)m->m_oaRingAtoms[z3])->GetSize() > i1)
					{
						i1 = ((CxIntArray*)m->m_oaRingAtoms[z3])->GetSize();
						ti = z3;
					}
				}
				if ((ti != -1) && (ti != z2))
				{
					wa = (CxIntArray*)m->m_oaRingAtoms[z2];
					m->m_oaRingAtoms[z2] = m->m_oaRingAtoms[ti];
					m->m_oaRingAtoms[ti] = wa;
					wa = (CxIntArray*)m->m_oaRingAtomTypes[z2];
					m->m_oaRingAtomTypes[z2] = m->m_oaRingAtomTypes[ti];
					m->m_oaRingAtomTypes[ti] = wa;
				}
			}
		}
	}

/*	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex[0]];
		mprintf("Molecule %d Ring systems:\n",z+1);
		for (z2=0;z2<sm->m_oaRings.GetSize();z2++)
		{
			wa = (CxIntArray*)sm->m_oaRings[z2];
			mprintf("  - Ring %d: ",z2+1);
			for (z3=0;z3<wa->GetSize();z3++)
			{
				mprintf("%d",(*wa)[z3]);
				if (z3+1 < wa->GetSize())
					mprintf(", ");
			}
			mprintf("\n");
		}
		mprintf("\n");
	}*/

/*	mprintf("\n%d Molekuelsorten gefunden.\n",g_oaMolecules.GetSize());
	for (z=0;z<g_oaMolecules.GetSize();z++)
	mprintf("Molekuel \"%s\": %d Atome gesamt, %d mal in der Simulation.\n",g_pMolecules[z].Name,g_pMolecules[z].AtomGes,g_pMolecules[z].Count);
	*/
/*	mprintf("\n");
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		mprintf("%d ",z+1);
		((CMolecule*)g_oaMolecules[z])->Dump();
	}*/

	mprintf("    Molecule recognition finished.\n");
	BTOUT; 
	return true;
}


void CTimeStep::RECURSION_MegaTree(int i, char *ta, int depth, unsigned long bmask, bool w, int *stack)
{
	BTIN;
	int z, z2;
	int nblist[20], nbs;

	stack[depth] = i;
	ta[i] = 1;
	for (z=1;z<depth;z++)
	{
		if ((bmask & (int)mypow(2.0,z)) != 0)
			mprintf(WHITE,"|  ");
				else mprintf("   ");
	}
	if (depth != 0)
	{
		if (w)
			mprintf(WHITE,"|--");
				else mprintf(WHITE,"`--");
	}
	mprintf("%s(%d)",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[i]])->m_sName,i+1);
	nbs = 0;
	for (z=0;z<(long)m_iGesAtomCount;z++) // Schon mal alle Nachbarn raussuchen
	{
		if (z == i)
			continue;
		if (BondRange(i,z,NULL))
		{
			if (ta[z] != 0)
			{
				if ((depth > 0) && (z != stack[depth-1]))
				{
					mprintf(GREEN,"  <-- Ring closure: ");
					mprintf("%s(%d)",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[stack[depth]]])->m_sName,stack[depth]+1);
					z2 = depth-1;
					while ((stack[z2] != z) && (z2 >= 0))
					{
						mprintf(" - %s(%d)",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[stack[z2]]])->m_sName,stack[z2]+1);
						z2--;
					}
					mprintf(" - %s(%d)",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[z]])->m_sName,z+1);
//					mprintf("Ringschluss! [d-1=%d,d=%d,i=%d] %s(%d)\n",stack[depth-1]+1,stack[depth]+1,i+1,((CAtom*)g_oaAtoms[g_baAtomIndex[z]])->m_sName,z+1);
				}
				continue;
			}
			nblist[nbs] = z;
			nbs++;
		}
	}
	mprintf("\n");
	for (z=0;z<nbs;z++) // Fuer das aktuelle Atom z2 alle Nachbarn finden
	{
		if (ta[nblist[z]] != 0) // Der Nachbar wurde uns weggeschnappt -> Ringschluss
		{
			for (z2=1;z2<depth+1;z2++)
			{
				if ((bmask & (int)mypow(2.0,z2)) != 0)
					mprintf(WHITE,"|  ");
						else mprintf("   ");
			}
			if (z+1 < nbs)
				mprintf(WHITE,"|--");
					else mprintf("`--");
			mprintf("%s(%d)",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[nblist[z]]])->m_sName,nblist[z]+1);
			mprintf(GREEN,"  <-- Ring closure\n");
		} else
		{
			if (z+1 == nbs)
				bmask -= (int)mypow(2.0,depth+1);
			RECURSION_MegaTree(nblist[z],ta,/*sm,*/depth+1,bmask,(z+1==nbs)?false:true,stack);
		}
	}
	BTOUT; 
}


void CTimeStep::PrintMegaTree()
{
	BTIN;
	int z;
//	char ta[16384];
//	int stack[16384];
	char *ta;
	int *stack;

	try { ta = new char[m_iGesAtomCount]; } catch(...) { ta = NULL; }
	if (ta == NULL) NewException((double)m_iGesAtomCount*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { stack = new int[m_iGesAtomCount]; } catch(...) { stack = NULL; }
	if (stack == NULL) NewException((double)m_iGesAtomCount*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z=0;z<(int)m_iGesAtomCount;z++)
		ta[z] = 0;

	for (z=0;z<(int)m_iGesAtomCount;z++)
	{
		if (ta[z] != 0) // Dieses Atom wurde bereits in irgendein Molekuel eingebaut
			continue;
		mprintf(YELLOW,"\nThe next molecule starts with %s(%d):\n",(const char*)m_paLabels[z],z+1);
		RECURSION_MegaTree(z,ta,0,0xFFFFFFFF,true,stack);
	}

	delete[] stack;
	delete[] ta;

	mprintf("\n");
	BTOUT; 
}


void CTimeStep::PrintMatrix(bool onlyfirst, bool onlybind)
{
	BTIN;
	int z, z2, z3, z4, z5, z6, z7, z8, ti, ti2, hc;
	double tf;
	bool b, c, noh, hex, log;
	CxDVector3 vec1, vec2;
	CSingleMolecule *sm;
	CMolecule *m;
	CMolAtom *ma, *ma0;
	CxIntArray *wa, *wat;


	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		for (z6=0;z6<(onlyfirst?1:m->m_laSingleMolIndex.GetSize());z6++)
		{
			mprintf(YELLOW,"\n*** Molecule %d: %s; Representant %d; Distances in pm ***\n",z+1,m->m_sName,z6+1);
			if ((m->m_baAtomIndex.GetSize() == 1) && (m->m_waAtomCount[0] == 1))
			{
				mprintf("[only 1 Atom]\n");
				continue;
			}
			hex = false;
			noh = false;
			hc = 0;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++) {
				if (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,"H") == 0) {
					hex = true;
					hc += m->m_waAtomCount[z2];
				}
			}
_matagain:
			log = false;
			if (noh) {
				if (m->m_iAtomGes-hc > 55) {
					mprintf(WHITE,"\n    This molecule has many atoms (even without H). Output would flood the screen. Only writing to travis.log file.\n");
					log = true;
				}
			} else {
				if (m->m_iAtomGes > 55) {
					mprintf(WHITE,"\n    This molecule has many atoms. Output would flood the screen. Only writing to travis.log file.\n");
					log = true;
				}
			}
			if (!log && (m->m_oaRingAtoms.GetSize() != 0))
				mprintf(GREEN,"    (ring bonds are shown in green)\n");
			if (log)
				mprintf_nos("\n   ");
			else
				mprintf("\n   ");
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				if (noh) // Skip H atoms
					if (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,"H") == 0)
						continue;

				for (z3=0;z3<m->m_waAtomCount[z2];z3++)
				{
					if (((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName[1] == 0) {
						if (log)
							mprintf_nos(" ");
						else
							mprintf(" ");
					}
					if (log)
						mprintf_nos("%s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1);	
					else
						mprintf(WHITE,"%s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1);	
					if (z3+1 < 10) {
						if (log)
							mprintf_nos(" ");
						else
							mprintf(" ");
					}
				}
			}
			if (log)
				mprintf_nos("\n");
			else
				mprintf("\n");
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				if (noh) // Skip H atoms
					if (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,"H") == 0)
						continue;
				for (z3=0;z3<m->m_waAtomCount[z2];z3++)
				{
					if (log)
						mprintf_nos("%s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1);
					else
						mprintf(WHITE,"%s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1);
					if (((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName[1] == 0) {
						if (log)
							mprintf_nos(" ");
						else
							mprintf(" ");
					}
					if (z3+1 < 10) {
						if (log)
							mprintf_nos(" ");
						else
							mprintf(" ");
					}
					for (z4=0;z4<m->m_baAtomIndex.GetSize();z4++)
					{
						if (noh) // Skip H atoms
							if (mystricmp(((CAtom*)g_oaAtoms[m->m_baAtomIndex[z4]])->m_sName,"H") == 0)
								continue;
						for (z5=0;z5<m->m_waAtomCount[z4];z5++)
						{
							if ((z2 == z4) && (z3 == z5))
							{
								if (log)
									mprintf_nos("*** ");
								else
									mprintf(BLUE,"*** ");
								continue;
							}
							ti = ((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z6]])->m_oaAtomOffset[z2])->GetAt(z3);
							//vec1 = g_TimeStep.m_vaCoords[ti];
							ti2 = ((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z6]])->m_oaAtomOffset[z4])->GetAt(z5);
							//vec2 = g_TimeStep.m_vaCoords[ti2];
							//vec1 -= vec2;
							vec1 = FoldVector(g_TimeStep.m_vaCoords[ti] - g_TimeStep.m_vaCoords[ti2]);
							tf = vec1.GetLength();
							if ((g_TimeStep.BondRange(ti,ti2,NULL)) || (!onlybind))
							{
								for (z8=0;z8<m->m_oaRingAtomTypes.GetSize();z8++)
								{
									wat = (CxIntArray*)m->m_oaRingAtomTypes[z8];
									wa = (CxIntArray*)m->m_oaRingAtoms[z8];
									for (z7=0;z7<wa->GetSize();z7++)
									{
										if (((*wat)[z7] == z2) && ((*wa)[z7] == z3))
										{
											if (z7 == 0)
											{
												if (((*wat)[z7+1] == z4) && ((*wa)[z7+1] == z5))
													goto _green;
												if (((*wat)[wa->GetSize()-1] == z4) && ((*wa)[wa->GetSize()-1] == z5))
													goto _green;
											} else if (z7 == wa->GetSize()-1)
											{
												if (((*wat)[0] == z4) && ((*wa)[0] == z5))
													goto _green;
												if (((*wat)[z7-1] == z4) && ((*wa)[z7-1] == z5))
													goto _green;
											} else
											{
												if (((*wat)[z7+1] == z4) && ((*wa)[z7+1] == z5))
													goto _green;
												if (((*wat)[z7-1] == z4) && ((*wa)[z7-1] == z5))
													goto _green;
											}
										}
									}
								}
								if (log)
									mprintf_nos("%3.0f ",tf);
								else
									mprintf("%3.0f ",tf);
								continue;
_green:
								if (log)
									mprintf_nos("%3.0f ",tf);
								else
									mprintf(GREEN,"%3.0f ",tf);
							} else {
								if (log)
									mprintf_nos(" -  ");
								else
									mprintf(" -  ");
							}
						}
					}
					if (log)
						mprintf_nos("\n");
					else
						mprintf("\n");
				}
			}
			if (hex && (hc > 10) && (!noh) && (m->m_iAtomGes > 30))
			{
				mprintf(WHITE,"\n    This was a very large molecule.\n    Printing the matrix again without H atoms.\n");
				noh = true;
				goto _matagain;
			}
//			if (((!onlyfirst) && ((z2 < ((CMolecule*)g_oaMolecules[z])->Elements-1) || (g_oaMolecules.GetSize()-1))) || (z < g_oaMolecules.GetSize()-1))
//				mprintf(">>>");
		}

		if (m->m_oaRingAtoms.GetSize() != 0)
		{
			mprintf("\n");
			for (z2=0;z2<m->m_oaRingAtoms.GetSize();z2++)
			{
				mprintf("    %d-ring: ",((CxIntArray*)m->m_oaRingAtoms[z2])->GetSize());
				if (((CxIntArray*)m->m_oaRingAtoms[z2])->GetSize() > 20)
				{
					for (z3=0;z3<5;z3++)
					{
						mprintf("%s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[((CxIntArray*)m->m_oaRingAtomTypes[z2])->GetAt(z3)]])->m_sName,((CxIntArray*)m->m_oaRingAtoms[z2])->GetAt(z3)+1);
						if (z3+1 < ((CxIntArray*)m->m_oaRingAtoms[z2])->GetSize())
							mprintf(" - ");
					}
					mprintf(RED,"...");
					mprintf(" - ");
					for (z3=((CxIntArray*)m->m_oaRingAtoms[z2])->GetSize()-5;z3<((CxIntArray*)m->m_oaRingAtoms[z2])->GetSize();z3++)
					{
						mprintf("%s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[((CxIntArray*)m->m_oaRingAtomTypes[z2])->GetAt(z3)]])->m_sName,((CxIntArray*)m->m_oaRingAtoms[z2])->GetAt(z3)+1);
						if (z3+1 < ((CxIntArray*)m->m_oaRingAtoms[z2])->GetSize())
							mprintf(" - ");
					}
				} else
				{
					for (z3=0;z3<((CxIntArray*)m->m_oaRingAtoms[z2])->GetSize();z3++)
					{
						mprintf("%s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[((CxIntArray*)m->m_oaRingAtomTypes[z2])->GetAt(z3)]])->m_sName,((CxIntArray*)m->m_oaRingAtoms[z2])->GetAt(z3)+1);
						if (z3+1 < ((CxIntArray*)m->m_oaRingAtoms[z2])->GetSize())
							mprintf(" - ");
					}
				}
				mprintf("\n");
				if (z2 >= 99)
				{
					mprintf(RED,"\n    Only showing the first 100 rings.\n");
					break;
				}
			}
		}

		c = false;
		sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
		for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
		{
			b = false;
			ma = NULL;
			//for (z3=0;z3<sm->m_oaMolAtoms.GetSize();z3++)
			for (z3=0;z3<sm->m_oaMolAtomsSorted.GetSize();z3++)
			{
				//if (((CMolAtom*)sm->m_oaMolAtoms[z3])->m_iType != z2)
				if (((CMolAtom*)sm->m_oaMolAtomsSorted[z3])->m_iType != z2)
					continue;
				ma0 = ma;
				//ma = (CMolAtom*)sm->m_oaMolAtoms[z3];
				ma = (CMolAtom*)sm->m_oaMolAtomsSorted[z3];
				if (ma0 == NULL)
					continue;
				z4 = ma0->m_iNumber;
//				if (ma->m_fAtomCode == ma0->m_fAtomCode)
				if (ma->m_liAtomCode == ma0->m_liAtomCode)
				{
					if (b)
					{
						mprintf(", %s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
					} else
					{
						if (!c)
						{
							c = true;
							mprintf("\n");
						}
						mprintf("    Atoms %s%d, %s%d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z4+1,(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,ma->m_iNumber+1);
						b = true;
					}
				} else if (b)
				{
					mprintf(" are equivalent.\n");
					b = false;
				}
			}
			if (b)
				mprintf(" are equivalent.\n");
		}
	} // Tabelle Ende
	BTOUT; 
}


/*void CTimeStep::Transform(const CxDMatrix3 &mat)
{
	BTIN;
	int z;
//	CxDVector3 v;

	for (z=0;z<g_iGesVirtAtomCount;z++)
	{
//		if (z < g_iGesAtomCount)
//			mprintf("@ %2d: ( %g | %g | %g ) --> ",z,m_vaCoords[z][0],m_vaCoords[z][1],m_vaCoords[z][2]);
		m_vaCoords[z] = mat * m_vaCoords[z];
//		if (z < g_iGesAtomCount)
//			mprintf("( %g | %g | %g )\n",m_vaCoords[z][0],m_vaCoords[z][1],m_vaCoords[z][2]);
	}

	if (g_bUseVelocities)
		for (z=0;z<g_iGesVirtAtomCount;z++)
			m_vaVelocities[z] = mat * m_vaVelocities[z];

	if (g_bUseForces)
		for (z=0;z<g_iGesVirtAtomCount;z++)
			m_vaForces[z] = mat * m_vaForces[z];
	BTOUT; 
}*/


void CTimeStep::Transform(const CxDMatrix3 &mat)
{
	BTIN;
	int z;

	for (z=0;z<g_iGesVirtAtomCount;z++)
		m_vaCoords[z] = mat * m_vaCoords[z];

	if (g_bUseVelocities)
		for (z=0;z<g_iGesVirtAtomCount;z++)
			m_vaVelocities[z] = mat * m_vaVelocities[z];

	if (g_bUseForces)
		for (z=0;z<g_iGesVirtAtomCount;z++)
			m_vaForces[z] = mat * m_vaForces[z];
	BTOUT; 
}


void CTimeStep::Transform(const CxQuaternion &q)
{
	BTIN;
	int z;

	for (z=0;z<g_iGesVirtAtomCount;z++)
		m_vaCoords[z] = q.Transform(m_vaCoords[z]);

	if (g_bUseVelocities)
		for (z=0;z<g_iGesVirtAtomCount;z++)
			m_vaVelocities[z] = q.Transform(m_vaVelocities[z]);

	if (g_bUseForces)
		for (z=0;z<g_iGesVirtAtomCount;z++)
			m_vaForces[z] = q.Transform(m_vaForces[z]);
	BTOUT; 
}


void CTimeStep::SubVelocities(const CxDVector3 &vec)
{
	BTIN;
	int z;

	for (z=0;z<g_iGesVirtAtomCount;z++)
		m_vaVelocities[z] -= vec;
	BTOUT; 
}


void CTimeStep::FoldMolecules()
{
	BTIN;
	int z, z2, z3, z4, n;
	CxDVector3 v, vc;
	CMolecule *m;
	CSingleMolecule *sm;

	if (!g_bPeriodic)
		return;

//	mprintf("*** Fold ***\n");
	
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		if (m->m_bPseudo)
			continue;

		if (m->m_bPolymer) // Polymers need to be wrapped atom-wise
		{
			for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
				for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
				{
					for (z4=0;z4<m->m_waAtomCount[z3];z4++)
					{
						if (g_bBoxNonOrtho) {

							v = g_mBoxToOrtho * m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)];

							n = PeriodicImage1D( v[0], -0.5, 0.5 );
							if (n != 0)
								v[0] -= n;

							n = PeriodicImage1D( v[1], -0.5, 0.5 );
							if (n != 0)
								v[1] -= n;

							n = PeriodicImage1D( v[2], -0.5, 0.5 );
							if (n != 0)
								v[2] -= n;

			/*				while (v[0] > 0.5)
								v[0] -= 1.0;
							while (v[0] <= -0.5)
								v[0] += 1.0;
							while (v[1] > 0.5)
								v[1] -= 1.0;
							while (v[1] <= -0.5)
								v[1] += 1.0;
							while (v[2] > 0.5)
								v[2] -= 1.0;
							while (v[2] <= -0.5)
								v[2] += 1.0;*/

							m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)] = g_mBoxFromOrtho * v;

						} else {

							v = m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)];

							if (g_bPeriodicX) {
								n = PeriodicImage1D( v[0], -g_fBoxX/2.0, g_fBoxX/2.0 );
								if (n != 0)
									v[0] -= n*g_fBoxX;
							}

							if (g_bPeriodicY) {
								n = PeriodicImage1D( v[1], -g_fBoxY/2.0, g_fBoxY/2.0 );
								if (n != 0)
									v[1] -= n*g_fBoxY;
							}

							if (g_bPeriodicZ) {
								n = PeriodicImage1D( v[2], -g_fBoxZ/2.0, g_fBoxZ/2.0 );
								if (n != 0)
									v[2] -= n*g_fBoxZ;
							}

			/*				if (g_bPeriodicX)
							{
								while (v[0] < -g_fBoxX/2) { v[0] += g_fBoxX; vc[0] += g_fBoxX; }
								while (v[0] >= g_fBoxX/2) { v[0] -= g_fBoxX; vc[0] -= g_fBoxX; }
							}
							if (g_bPeriodicY)
							{
								while (v[1] < -g_fBoxY/2) { v[1] += g_fBoxY; vc[1] += g_fBoxY; }
								while (v[1] >= g_fBoxY/2) { v[1] -= g_fBoxY; vc[1] -= g_fBoxY; }
							}
							if (g_bPeriodicZ)
							{
								while (v[2] < -g_fBoxZ/2) { v[2] += g_fBoxZ; vc[2] += g_fBoxZ; }
								while (v[2] >= g_fBoxZ/2) { v[2] -= g_fBoxZ; vc[2] -= g_fBoxZ; }
							}*/

							m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)] = v;
						}
					}
				}
			}
		} else
		{
			for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];

	//			mprintf("\n%d, %d",m->m_baAtomIndex.GetSize(),((CxIntArray*)sm->m_oaAtomOffset[m->m_baAtomIndex.GetSize()-1])->GetAt(1));

				// Massenzentrum

	/*			mprintf("  .Mol ");
				v.Dump();
				mprintf("\n");			

	*/			
				vc = CxDVector3(0,0,0);

				if (g_bBoxNonOrtho) {

					v = g_mBoxToOrtho *  m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[m->m_baAtomIndex.GetSize()-1])->GetAt(1)];

					n = PeriodicImage1D( v[0], -0.5, 0.5 );
					if (n != 0)
						vc[0] -= n;

					n = PeriodicImage1D( v[1], -0.5, 0.5 );
					if (n != 0)
						vc[1] -= n;

					n = PeriodicImage1D( v[2], -0.5, 0.5 );
					if (n != 0)
						vc[2] -= n;

	/*				while (v[0] > 0.5)
					{
						v[0] -= 1.0;
						vc[0] -= 1.0;
					}
					while (v[0] <= -0.5)
					{
						v[0] += 1.0;
						vc[0] += 1.0;
					}
					while (v[1] > 0.5)
					{
						v[1] -= 1.0;
						vc[1] -= 1.0;
					}
					while (v[1] <= -0.5)
					{
						v[1] += 1.0;
						vc[1] += 1.0;
					}
					while (v[2] > 0.5)
					{
						v[2] -= 1.0;
						vc[2] -= 1.0;
					}
					while (v[2] <= -0.5)
					{
						v[2] += 1.0;
						vc[2] += 1.0;
					}*/

					vc = g_mBoxFromOrtho * vc;

				} else {

					v = m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[m->m_baAtomIndex.GetSize()-1])->GetAt(1)];

					if (g_bPeriodicX) {
						n = PeriodicImage1D( v[0], -g_fBoxX/2.0, g_fBoxX/2.0 );
						if (n != 0)
							vc[0] -= n*g_fBoxX;
					}

					if (g_bPeriodicY) {
						n = PeriodicImage1D( v[1], -g_fBoxY/2.0, g_fBoxY/2.0 );
						if (n != 0)
							vc[1] -= n*g_fBoxY;
					}

					if (g_bPeriodicZ) {
						n = PeriodicImage1D( v[2], -g_fBoxZ/2.0, g_fBoxZ/2.0 );
						if (n != 0)
							vc[2] -= n*g_fBoxZ;
					}

		/*			if (g_bPeriodicX)
					{
						while (v[0] < -g_fBoxX/2) { v[0] += g_fBoxX; vc[0] += g_fBoxX; }
						while (v[0] >= g_fBoxX/2) { v[0] -= g_fBoxX; vc[0] -= g_fBoxX; }
					}
					if (g_bPeriodicY)
					{
						while (v[1] < -g_fBoxY/2) { v[1] += g_fBoxY; vc[1] += g_fBoxY; }
						while (v[1] >= g_fBoxY/2) { v[1] -= g_fBoxY; vc[1] -= g_fBoxY; }
					}
					if (g_bPeriodicZ)
					{
						while (v[2] < -g_fBoxZ/2) { v[2] += g_fBoxZ; vc[2] += g_fBoxZ; }
						while (v[2] >= g_fBoxZ/2) { v[2] -= g_fBoxZ; vc[2] -= g_fBoxZ; }
					}*/
				}
				
				if ((vc[0] == 0) && (vc[1] == 0) && (vc[2] == 0))
					 continue;
				
				for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
					for (z4=0;z4<m->m_waAtomCount[z3];z4++)
						 m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)] += vc;

				if (g_bWannier)
					for (z3=0;z3<sm->m_laWannier.GetSize();z3++)
						 m_vaCoords[sm->m_laWannier[z3]] += vc;
			}
		}
	}
//	mprintf("   Fold Done\n");
	BTOUT;
}


void CTimeStep::FoldAtoms()
{
	BTIN;
	CxDVector3 v;
	int z, n;


	if (!g_bPeriodic)
		return;

	for (z=0;z<g_iGesVirtAtomCount;z++) {

		if (g_bBoxNonOrtho) {

			v = g_mBoxToOrtho *  m_vaCoords[z];

			n = PeriodicImage1D( v[0], -0.5, 0.5 );
			if (n != 0)
				v[0] -= n;

			n = PeriodicImage1D( v[1], -0.5, 0.5 );
			if (n != 0)
				v[1] -= n;

			n = PeriodicImage1D( v[2], -0.5, 0.5 );
			if (n != 0)
				v[2] -= n;

/*			while (v[0] > 0.5)
				v[0] -= 1.0;
			while (v[0] <= -0.5)
				v[0] += 1.0;
			while (v[1] > 0.5)
				v[1] -= 1.0;
			while (v[1] <= -0.5)
				v[1] += 1.0;
			while (v[2] > 0.5)
				v[2] -= 1.0;
			while (v[2] <= -0.5)
				v[2] += 1.0;*/

			m_vaCoords[z] = g_mBoxFromOrtho * v;

		} else {

			if (g_bPeriodicX) {
				n = PeriodicImage1D( m_vaCoords[z][0], -g_fBoxX/2.0, g_fBoxX/2.0 );
				if (n != 0)
					m_vaCoords[z][0] -= n*g_fBoxX;
			}

			if (g_bPeriodicY) {
				n = PeriodicImage1D( m_vaCoords[z][1], -g_fBoxY/2.0, g_fBoxY/2.0 );
				if (n != 0)
					m_vaCoords[z][1] -= n*g_fBoxY;
			}

			if (g_bPeriodicZ) {
				n = PeriodicImage1D( m_vaCoords[z][2], -g_fBoxZ/2.0, g_fBoxZ/2.0 );
				if (n != 0)
					m_vaCoords[z][2] -= n*g_fBoxZ;
			}

/*			if (g_bPeriodicX)
			{
				while (m_vaCoords[z][0] < -g_fBoxX/2) m_vaCoords[z][0] += g_fBoxX;
				while (m_vaCoords[z][0] >= g_fBoxX/2) m_vaCoords[z][0] -= g_fBoxX;
			}

			if (g_bPeriodicY)
			{
				while (m_vaCoords[z][1] < -g_fBoxY/2) m_vaCoords[z][1] += g_fBoxY;
				while (m_vaCoords[z][1] >= g_fBoxY/2) m_vaCoords[z][1] -= g_fBoxY;
			}

			if (g_bPeriodicZ)
			{
				while (m_vaCoords[z][2] < -g_fBoxZ/2) m_vaCoords[z][2] += g_fBoxZ;
				while (m_vaCoords[z][2] >= g_fBoxZ/2) m_vaCoords[z][2] -= g_fBoxZ;
			}*/
		}
	}
	BTOUT; 
}


void CTimeStep::CenterPos(const CxDVector3 &vec)
{
	BTIN;
	int z/*, z2, z3, z4*/;

/*	 for (z=0;z<g_oaMolecules.GetSize();z++)
		 for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->Count;z2++)
			 for (z3=0;z3<((CMolecule*)g_oaMolecules[z])->Elements;z3++)
				 for (z4=0;z4<((CMolecule*)g_oaMolecules[z])->AtomCount[z3];z4++)
					 m_vaCoords[((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[z])->SingleMolIndex[z2]])->m_iAtomOffset[z3][z4]] -= vec;*/

	for (z=0;z<g_iGesVirtAtomCount;z++)
		m_vaCoords[z] -= vec;
	BTOUT; 
}


void CTimeStep::AddAtoms()
{
	BTIN;
	char *q;
	int z;
//	bool labeledatoms;

//	labeledatoms = false;
	for (z=0;z<(long)m_iGesAtomCount;z++)
	{
		q = (char*)m_paLabels[z];
		q += strlen(q)-1;
/*		while ((*q >= '0') && (*q <= '9'))
		{
			if (!labeledatoms)
			{
				mprintf(">>>\n>>> The atoms in the input file are numbered.\n>>> Ignoring this, using own numbers ;-)\n>>>\n");
				labeledatoms = true;
			}
			*q = 0;
			q--;
		}*/
//		printf("AddAtoms: \"%s\"\n",q);
		xAddAtom((char*)m_paLabels[z],false);
	}
	BTOUT; 
}


void CTimeStep::WriteTimestep(FILE *a)
{
	BTIN;
	int z, z2, z3, z4, z0;
	CMolecule *m;
	CSingleMolecule *sm;

	if (g_bSaveVirtAtoms)
		mfprintf(a,"  %d\n",g_iGesVirtAtomCount);
			else mfprintf(a,"  %d\n",g_iGesAtomCount);
	if (m_pComment != NULL)
		mfprintf(a,"%s\n",m_pComment);
			else mfprintf(a,"\n");

	if (g_bWriteAtomwise)
	{
		for (z0=0;z0<g_oaAtoms.GetSize();z0++)
		{
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				m = (CMolecule*)g_oaMolecules[z];
				if (m->m_bPseudo)
					continue;
				for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
				{
					if (m->m_baAtomIndex[z3] != z0)
						continue;
					if ((!g_bSaveVirtAtoms) && (m->m_baAtomIndex[z3] == g_iVirtAtomType))
						continue;
					for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
					{
						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
						for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
							mfprintf(a,"  %s  %8.5f  %8.5f  %8.5f\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][0]/100.0,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][1]/100.0,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][2]/100.0);
					}
				}
			}
		}
	} else
	{
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
						mfprintf(a,"  %s  %8.5f  %8.5f  %8.5f\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][0]/100.0,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][1]/100.0,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][2]/100.0);
				}
			}
		}
	}
			
	BTOUT; 
}


void CTimeStep::ExportSingleMolecule_PDB(CSingleMolecule *sm, const char *s) {

	FILE *a;
	CMolecule *m;
	std::vector<std::vector<int> > iaa;
	int z, z2, z3, z4, z5, ti, ti2, o, o2;
	CMolBond *mb;


	a = OpenFileWrite(s,true);

	mfprintf(a,"COMPND    UNNAMED\n");
	mfprintf(a,"AUTHOR    GENERATED BY TRAVIS\n");

	m = (CMolecule*)g_oaMolecules[sm->m_iMolType];

	iaa.resize(m->m_iAtomGesNoVirt);

	ti = 0;
	for (z=0;z<m->m_baAtomIndex.GetSize();z++) {

		if (m->m_baAtomIndex[z] == g_iVirtAtomType)
			continue;

		for (z2=0;z2<((CxIntArray*)sm->m_oaAtomOffset[z])->GetSize();z2++) {

			o = ((CxIntArray*)sm->m_oaAtomOffset[z])->GetAt(z2);

			mfprintf(
				a,
				"HETATM%5d %4s LIG     1    %8.3f%8.3f%8.3f  1.00  0.00          %2s  \n",
				ti+1,
				(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z]])->m_sName,
				m_vaCoords[o][0]/100.0,
				m_vaCoords[o][1]/100.0,
				m_vaCoords[o][2]/100.0,
				(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z]])->m_sName
			);

			ti2 = 0;
			for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++) {

				if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
					continue;

				for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++) {

					o2 = ((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4);

					for (z5=0;z5<sm->m_oaBonds.GetSize();z5++) {

						mb = (CMolBond*)sm->m_oaBonds[z5];

						if (((mb->m_iAtomOffset[0] == o) && (mb->m_iAtomOffset[1] == o2)) ||
							((mb->m_iAtomOffset[1] == o) && (mb->m_iAtomOffset[0] == o2))) {
							iaa[ti].push_back(ti2);
							break;
						}
					}

					ti2++;
				}
			}

			ti++;
		}
	}

	for (z=0;z<m->m_iAtomGesNoVirt;z++) {
		ti2 = 0;
_next:
		mfprintf(a,"CONECT%5d",z+1);
		for (z2=0;z2<4;z2++) {
			if (ti2 < (int)iaa[z].size())
				mfprintf(a,"%5d",iaa[z][ti2]+1);
			else
				mfprintf(a,"     ");
			ti2++;
		}
		mfprintf(a,"                                       \n");
		if (ti2 < (int)iaa[z].size())
			goto _next;
	}

	mfprintf(a,"END\n");

	fclose(a);
}

	
void CTimeStep::WriteTimestepNb(FILE *a, CNbSet *nbs, int singlemol)
{
	BTIN;
	int z, z2, z3, z4, z0, n, ti;
	CMolecule *m;
	CSingleMolecule *sm;
	CConditionGroup *cg;
	CConditionSubGroup *cs;
	CNbSearch *nb;

	n = 0;

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		if (nbs->m_oaConditionGroups[z] == NULL)
			continue;
		m = (CMolecule*)g_oaMolecules[z];
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			if (((CConditionGroup*)nbs->m_oaConditionGroups[z])->Contains(z2))
			{
				n += m->m_iAtomGes;
				if (!g_bSaveVirtAtoms)
					n -= m->m_laVirtualAtoms.GetSize();
			}
		}
	}

	mfprintf(a,"  %d\n",n);

	mfprintf(a,"# Step %d",((int)g_iSteps));

	if (singlemol >= 0)
	{
		sm = (CSingleMolecule*)g_oaSingleMolecules[singlemol];
		m = (CMolecule*)g_oaMolecules[sm->m_iMolType];
		if (!g_bSaveVirtAtoms)
		{
			if (g_bWriteAtomwise)
				mfprintf(a,", RM=%s[%d]",m->m_sName,sm->m_iMolSMIndex+1);
					else mfprintf(a,", RM=%s[%d] (%d atoms)",m->m_sName,sm->m_iMolSMIndex+1,m->m_iAtomGes-m->m_laVirtualAtoms.GetSize());
		} else
		{
			if (g_bWriteAtomwise)
				mfprintf(a,", RM=%s[%d]",m->m_sName,sm->m_iMolSMIndex+1);
					else mfprintf(a,", RM=%s[%d] (%d atoms)",m->m_sName,sm->m_iMolSMIndex+1,m->m_iAtomGes);
		}
	}

	if (g_bEnvWriteDetailedInfo)
	{
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			if (nbs->m_oaConditionGroups[z] == NULL)
				continue;
			m = (CMolecule*)g_oaMolecules[z];

			cg = (CConditionGroup*)nbs->m_oaConditionGroups[z];

			if (g_bEnvSortNb)
			{
				nb = NULL;
				if (cg->m_oaConditionSubGroups.GetSize() == 1)
				{
					cs = (CConditionSubGroup*)cg->m_oaConditionSubGroups[0];
					if (cs->m_oaConditions.GetSize() == 1)
						nb = (CNbSearch*)cs->m_oaConditions[0];
				}
				if (nb == NULL)
					goto _nosortinfo;

				if (nb->m_iNbCountMin <= -1)
					goto _nosortinfo;

				for (z2=nb->m_iNbCountMin;z2<=nb->m_iNbCountMax;z2++)
				{
					ti = nb->m_pNbSort[z2].m_iOM;
					if (!g_bSaveVirtAtoms)
					{
						if (m->m_laSingleMolIndex[ti] != singlemol)
						{
							if (g_bWriteAtomwise)
								mfprintf(a,", %s[%d] d=%.3fpm",m->m_sName,z2+1,nb->m_pNbSort[z2].m_fMinDist);
									else mfprintf(a,", %s[%d] d=%.3fpm (%d atoms)",m->m_sName,ti+1,nb->m_pNbSort[z2].m_fMinDist,m->m_iAtomGes-m->m_laVirtualAtoms.GetSize());
						}
					} else
					{
						if (m->m_laSingleMolIndex[ti] != singlemol)
						{
							if (g_bWriteAtomwise)
								mfprintf(a,", %s[%d] d=%.3fpm",m->m_sName,z2+1,nb->m_pNbSort[z2].m_fMinDist);
									else mfprintf(a,", %s[%d] d=%.3fpm (%d atoms)",m->m_sName,ti+1,nb->m_pNbSort[z2].m_fMinDist,m->m_iAtomGes);
						}
					}
				}
			} else
			{
_nosortinfo:
				for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
				{
					if (cg->Contains(z2))
					{
						if (!g_bSaveVirtAtoms)
						{
							if (m->m_laSingleMolIndex[z2] != singlemol)
							{
								if (g_bWriteAtomwise)
									mfprintf(a,", %s[%d]",m->m_sName,z2+1);
										else mfprintf(a,", %s[%d] (%d atoms)",m->m_sName,z2+1,m->m_iAtomGes-m->m_laVirtualAtoms.GetSize());
							}
						} else
						{
							if (m->m_laSingleMolIndex[z2] != singlemol)
							{
								if (g_bWriteAtomwise)
									mfprintf(a,", %s[%d]",m->m_sName,z2+1);
										else mfprintf(a,", %s[%d] (%d atoms)",m->m_sName,z2+1,m->m_iAtomGes);
							}
						}
					}
				}
			}
		}
	}

	mfprintf(a,"\n");

	if (g_bWriteAtomwise)
	{
		for (z0=0;z0<g_oaAtoms.GetSize();z0++)
		{
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				if (nbs->m_oaConditionGroups[z] == NULL)
					continue;
				m = (CMolecule*)g_oaMolecules[z];
				for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
				{
					if (m->m_baAtomIndex[z3] != z0)
						continue;
					if ((!g_bSaveVirtAtoms) && (m->m_baAtomIndex[z3] == g_iVirtAtomType))
						continue;
					for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
					{
						if (((CConditionGroup*)nbs->m_oaConditionGroups[z])->Contains(z2))
						{
							sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
							for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
								mfprintf(a,"  %s  %8.5f  %8.5f  %8.5f\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][0]/100.0f,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][1]/100.0f,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][2]/100.0f);
						}
					}
				}
			}
		}
	} else
	{
		if (singlemol >= 0)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[singlemol];
			m = (CMolecule*)g_oaMolecules[sm->m_iMolType];
			for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
			{
				if ((!g_bSaveVirtAtoms) && (m->m_baAtomIndex[z3] == g_iVirtAtomType))
					continue;
				for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
					mfprintf(a,"  %s  %8.5f  %8.5f  %8.5f\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][0]/100.0f,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][1]/100.0f,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][2]/100.0f);
			}
		}

		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			if (nbs->m_oaConditionGroups[z] == NULL)
				continue;
			m = (CMolecule*)g_oaMolecules[z];
			cg = (CConditionGroup*)nbs->m_oaConditionGroups[z];

			if (g_bEnvSortNb)
			{
				nb = NULL;
				if (cg->m_oaConditionSubGroups.GetSize() == 1)
				{
					cs = (CConditionSubGroup*)cg->m_oaConditionSubGroups[0];
					if (cs->m_oaConditions.GetSize() == 1)
						nb = (CNbSearch*)cs->m_oaConditions[0];
				}
				if (nb == NULL)
					goto _nosort;

				if (nb->m_iNbCountMin <= -1)
					goto _nosort;

				for (z2=nb->m_iNbCountMin;z2<=nb->m_iNbCountMax;z2++)
				{
					ti = nb->m_pNbSort[z2].m_iOM;

					if (m->m_laSingleMolIndex[ti] == singlemol)
						continue;

					sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[ti]];
					for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
					{
						if ((!g_bSaveVirtAtoms) && (m->m_baAtomIndex[z3] == g_iVirtAtomType))
							continue;
						for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
							mfprintf(a,"  %s  %8.5f  %8.5f  %8.5f\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][0]/100.0f,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][1]/100.0f,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][2]/100.0f);
					}
				}
			} else
			{
_nosort:
				for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
				{
					if (((CConditionGroup*)nbs->m_oaConditionGroups[z])->Contains(z2))
					{
						if (m->m_laSingleMolIndex[z2] == singlemol)
							continue;
						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
						for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
						{
							if ((!g_bSaveVirtAtoms) && (m->m_baAtomIndex[z3] == g_iVirtAtomType))
								continue;
							for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
								mfprintf(a,"  %s  %8.5f  %8.5f  %8.5f\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][0]/100.0f,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][1]/100.0f,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][2]/100.0f);
						}
					}
				}
			}
		}
	}
			
	BTOUT; 
}


/*void CTimeStep::WriteTimestepNb(int refmol, FILE *a)
{
	BTIN;
	int z, c, z2, z3, z4, z0;
	int ti;
	CMolecule *m, *m2;

	m = (CMolecule*)g_oaMolecules[g_iFixMol];

	c = 0;
	if (g_bRefEnvVirt)
	{
		for (z=0;z<g_oaMolecules.GetSize();z++)
			c += g_pNbAll->m_waScanNeighborCount[z]*((CMolecule*)g_oaMolecules[z])->m_iAtomGes;
		if (g_bSaveRefWithEnv)
			c += m->m_iAtomGes;
	} else
	{
		for (z=0;z<g_oaMolecules.GetSize();z++)
			c += g_pNbAll->m_waScanNeighborCount[z]*(((CMolecule*)g_oaMolecules[z])->m_iAtomGes-((CMolecule*)g_oaMolecules[z])->m_waVirtualAtoms.GetSize());
		if (g_bSaveRefWithEnv)
			c += m->m_iAtomGes - m->m_waVirtualAtoms.GetSize();
	}

	mfprintf(a,"  %d\n\n",c);

	if (g_bRefEnvAtomwise)
	{
		for (z0=0;z0<g_oaAtoms.GetSize();z0++)
		{
			if (g_bSaveRefWithEnv)
			{
				for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
				{
					if (m->m_baAtomIndex[z3] != z0)
						continue;
					if (!g_bRefEnvVirt)
						if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
							continue;
					for (z4=0;z4<m->m_waAtomCount[z3];z4++)
					{
						ti = ((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[refmol]])->m_oaAtomOffset[z3])->GetAt(z4);
						mfprintf(a,"  %s  %f  %f  %f\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,m_vaCoords[ti][0]/100.0f,m_vaCoords[ti][1]/100.0f,m_vaCoords[ti][2]/100.0f);
					}
				}
			}
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				m2 = (CMolecule*)g_oaMolecules[z];
				for (z2=0;z2<g_pNbAll->m_waScanNeighborCount[z];z2++)
				{
					for (z3=0;z3<m2->m_baAtomIndex.GetSize();z3++)
					{
						if (m2->m_baAtomIndex[z3] != z0)
							continue;
						if (!g_bRefEnvVirt)
							if (m2->m_baAtomIndex[z3] == g_iVirtAtomType)
								continue;
						for (z4=0;z4<m2->m_waAtomCount[z3];z4++)
						{
							ti = ((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[m2->m_laSingleMolIndex[((CxIntArray*)g_pNbAll->m_oaScanNeighbors[z])->GetAt(z2)]])->m_oaAtomOffset[z3])->GetAt(z4);
							mfprintf(a,"  %s  %f  %f  %f\n",((CAtom*)g_oaAtoms[m2->m_baAtomIndex[z3]])->m_sName,m_vaCoords[ti][0]/100.0f,m_vaCoords[ti][1]/100.0f,m_vaCoords[ti][2]/100.0f);
						}
					}
				}
			}
		}
	} else
	{
		if (g_bSaveRefWithEnv)
		{
			for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
			{
				if (!g_bRefEnvVirt)
					if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
						continue;
				for (z4=0;z4<m->m_waAtomCount[z3];z4++)
				{
					ti = ((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[refmol]])->m_oaAtomOffset[z3])->GetAt(z4);
					mfprintf(a,"  %s  %f  %f  %f\n",((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,m_vaCoords[ti][0]/100.0,m_vaCoords[ti][1]/100.0,m_vaCoords[ti][2]/100.0);
				}
			}
		}
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			m2 = (CMolecule*)g_oaMolecules[z];
			for (z2=0;z2<g_pNbAll->m_waScanNeighborCount[z];z2++)
			{
				for (z3=0;z3<m2->m_baAtomIndex.GetSize();z3++)
				{
					if (!g_bRefEnvVirt)
						if (m2->m_baAtomIndex[z3] == g_iVirtAtomType)
							continue;
					for (z4=0;z4<m2->m_waAtomCount[z3];z4++)
					{
						ti = ((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[m2->m_laSingleMolIndex[((CxIntArray*)g_pNbAll->m_oaScanNeighbors[z])->GetAt(z2)]])->m_oaAtomOffset[z3])->GetAt(z4);
						mfprintf(a,"  %s  %f  %f  %f\n",((CAtom*)g_oaAtoms[m2->m_baAtomIndex[z3]])->m_sName,m_vaCoords[ti][0]/100.0,m_vaCoords[ti][1]/100.0,m_vaCoords[ti][2]/100.0);
					}
				}
			}
		}
	}
	BTOUT; 
}*/


/*double CTimeStep::MolDist(CSingleMolecule *ref, CSingleMolecule *sm2, CNbSearch *nb)
{
	BXIN;
	int z3, z4;
	double d;
	CxDVector3 vec;

	d = 99999.0f;
	for (z3=0;z3<nb->m_waRefElements.GetSize();z3++)
		for (z4=0;z4<((CxIntArray*)nb->m_oaNbElements[sm2->m_iMolType])->GetSize();z4++)
		{
			vec = m_vaCoords[((CxIntArray*)ref->m_oaAtomOffset[nb->m_waRefElements[z3]])->GetAt(nb->m_waRefAtoms[z3])] - m_vaCoords[((CxIntArray*)sm2->m_oaAtomOffset[((CxIntArray*)nb->m_oaNbElements[sm2->m_iMolType])->GetAt(z4)])->GetAt(((CxIntArray*)nb->m_oaNbAtoms[sm2->m_iMolType])->GetAt(z4))];
			if (g_bFold)
			{
				while (vec[0] >= g_fBoxX/2) vec[0] -= g_fBoxX;
				while (vec[0] < -g_fBoxX/2) vec[0] += g_fBoxX;
				while (vec[1] >= g_fBoxY/2) vec[1] -= g_fBoxY;
				while (vec[1] < -g_fBoxY/2) vec[1] += g_fBoxY;
				while (vec[2] >= g_fBoxZ/2) vec[2] -= g_fBoxZ;
				while (vec[2] < -g_fBoxZ/2) vec[2] += g_fBoxZ;
			}
			if (vec.GetLength() < d)
				d = vec.GetLength();
		}
	BXOUT;
	if (d < 90000.0f)
		return d;
			else return -1.0f;
}*/


/*void CTimeStep::ScanNeighborhood(int fixmol, int refmol, CNbSearch *nb, CNbSearch *prev)
{
	BTIN;
	double *del, d;
	int *best, b;
	int z, z2, z3, c;
	CMolecule *m, *m2;
	CSingleMolecule *sm, *sm2;
	CxDVector3 vec;

	m = (CMolecule*)g_oaMolecules[fixmol];
	sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[refmol]];
	del = new double[g_oaSingleMolecules.GetSize()];
	best = new int[g_oaSingleMolecules.GetSize()];

//	printf("Suche Nachbarschaft von Molekuel %d...\n",refmol+1);

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		if (((!nb->m_bDistMode) && (nb->m_waMolCount[z] == 0)) || (nb->m_bDistMode && (nb->m_faMolDist[z] == 0)))
		{
			nb->m_waScanNeighborCount[z] = 0;
			continue;
		}
		m2 = (CMolecule*)g_oaMolecules[z];
//		printf("*** Molekuel %d: %s\n",z+1,m2->Name);
		for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
		{
			del[z2] = 99999.0f;
			if ((z == fixmol) && (z2 == refmol))
				continue;
			if (prev != NULL)
				if (!prev->Contains(z,z2))
					continue;
//			printf("  - Vertreter %d\n",z2+1);
			sm2 = (CSingleMolecule*)g_oaSingleMolecules[m2->m_laSingleMolIndex[z2]];
			del[z2] = MolDist(sm,sm2,nb);
//			printf("  - Finaler Abstand: %f\n",del[z2]);
		}
		if (nb->m_bDistMode)
		{
//			printf("**DistMode**\n");
//			mprintf("Die Nachbarn: ");
			c = 0;
			for (z3=0;z3<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z3++)
				if ((del[z3] <= nb->m_faMolDist[z]) && (del[z3] >= nb->m_faMolMinDist[z]))
				{
//					mprintf("%d (%f), ",z3+1,del[z3]);
					best[c] = z3;
					c++;
				}
//			mprintf("\n");
//			printf("  # Die naechsten Nachbarn: ");
			for (z2=0;z2<c;z2++)
			{
//				printf("%d, ",best[z2]+1);
				nb->m_iNeighbourCount++;
				for (z3=0;z3<nb->m_waScanNeighborCount[z];z3++)
					if (((CxIntArray*)nb->m_oaScanNeighbors[z])->GetAt(z3) == best[z2])
					{
						((CxIntArray*)nb->m_oaScanNeighborHits[z])->GetAt(z3)++;
						goto _enddist;
					}
				((CxIntArray*)nb->m_oaScanNeighbors[z])->Add(best[z2]);
				((CxIntArray*)nb->m_oaScanNeighborHits[z])->Add(1);
				nb->m_waScanNeighborCount[z]++;
_enddist:;
			}
//			printf("\n");
		} else
		{
//			printf("**CountMode**\n");
			for (z2=0;z2<nb->m_waMolCount[z];z2++)
			{
				d = 999999.0f;
				b = -1;
				for (z3=0;z3<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z3++)
				{
					if (del[z3] < d)
					{
						d = del[z3];
						b = z3;
					}
				}
				best[z2] = b;
				del[b] = 1000000.0f;
			}
//			printf("  # Die naechsten Nachbarn: ");
			for (z2=nb->m_waMolCountMin[z]-1;z2<nb->m_waMolCount[z];z2++)
			{
//				printf("%d, ",best[z2]);
				nb->m_iNeighbourCount++;
				for (z3=0;z3<nb->m_waScanNeighborCount[z];z3++)
					if (((CxIntArray*)nb->m_oaScanNeighbors[z])->GetAt(z3) == best[z2])
					{
						((CxIntArray*)nb->m_oaScanNeighborHits[z])->GetAt(z3)++;
						goto _end;
					}
				((CxIntArray*)nb->m_oaScanNeighbors[z])->Add(best[z2]);
				((CxIntArray*)nb->m_oaScanNeighborHits[z])->Add(1);
				nb->m_waScanNeighborCount[z]++;
_end:;
			}
//			printf("\n");
		}
	}
	delete del;
	delete best;
	BTOUT; 
}*/


/*void CTimeStep::ScanAngle(int fixmol, int refmol, CCondition *co, CNbSearch *prev)
{
	BTIN;
	int z2, z3, z4;
	double tf;
	CMolecule *m, *m2;
	CSingleMolecule *sm, *sm2;
	CxDVector3 vec0, vec1, vec2, vec3, vec4, vec5;
	CxIntArray tempwa;
	CNbSearch *nb;

//	mprintf("*** ScanAngle ***\n");
//	mprintf("  FixMol = %d, RefMol = %d\n",fixmol,refmol);
	m = (CMolecule*)g_oaMolecules[fixmol];
//	mprintf("  m->m_laSingleMolIndex[refmol] = %d\n",m->m_laSingleMolIndex[refmol]);

	nb = co->m_pTempNbS;
	sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[refmol]];

//	printf("Suche Nachbarschaft von Molekuel %d...\n",refmol+1);

//	mprintf("  SecondMol = %d\n",co->m_iSecondMol);
	m2 = (CMolecule*)g_oaMolecules[co->m_iSecondMol];

//	printf("*** Molekuel %d: %s\n",z+1,m2->Name);
	for (z2=0;z2<m2->m_laSingleMolIndex.GetSize();z2++)
	{
		if ((co->m_iSecondMol == fixmol) && (z2 == refmol))
			continue;
		if (prev != NULL)
			if (!prev->Contains(co->m_iSecondMol,z2))
				continue;
//		printf("  - Vertreter %d\n",z2+1);
		sm2 = (CSingleMolecule*)g_oaSingleMolecules[m2->m_laSingleMolIndex[z2]];

		co->m_pADF->BuildAtomList(sm,sm2,NULL,&tempwa);

		for (z4=0;z4<tempwa.GetSize();z4+=6)
		{
			if (co->m_pADF->m_bOrtho[0])
			{
				vec0 = m_vaCoords[tempwa[z4]];
				vec2 = m_vaCoords[tempwa[z4+1]];
				vec3 = m_vaCoords[tempwa[z4+2]];
				vec1 = CrossP(vec2-vec0,vec3-vec0);
			} else
			{
				vec0 = m_vaCoords[tempwa[z4]];
				vec2 = m_vaCoords[tempwa[z4+1]];
				vec1 = vec2-vec0;
			}
			if (co->m_pADF->m_bOrtho[1])
			{
				vec4 = m_vaCoords[tempwa[z4+3]];
				vec3 = m_vaCoords[tempwa[z4+4]];
				vec5 = m_vaCoords[tempwa[z4+5]];
				vec2 = CrossP(vec3-vec4,vec5-vec4);
			} else
			{
				vec4 = m_vaCoords[tempwa[z4+3]];
				vec3 = m_vaCoords[tempwa[z4+4]];
				vec2 = vec3-vec4;
			}

			tf = Angle_Deg(vec1,vec2);

  if ((tf >= co->m_pADF->m_fMinAngle) && (tf <= co->m_pADF->m_fMaxAngle))
			{
				nb->m_iNeighbourCount++;
				for (z3=0;z3<nb->m_waScanNeighborCount[co->m_iSecondMol];z3++)
					if (((CxIntArray*)nb->m_oaScanNeighbors[co->m_iSecondMol])->GetAt(z3) == z2)
					{
						((CxIntArray*)nb->m_oaScanNeighborHits[co->m_iSecondMol])->GetAt(z3)++;
						goto _endang;
					}
				((CxIntArray*)nb->m_oaScanNeighbors[co->m_iSecondMol])->Add(z2);
				((CxIntArray*)nb->m_oaScanNeighborHits[co->m_iSecondMol])->Add(1);
				nb->m_waScanNeighborCount[co->m_iSecondMol]++;
_endang:;		break;
			}
		}
	}
	BTOUT; 
}*/

/*void CTimeStep::GatherNbDiagram(int refmol, CNbSearch *nb)
{
	BTIN;
	double *del, d;
	int *best, b;
	bool *done;
	int z, z2, z3, z4, c, z0;
	CMolecule *m, *m2;
	CSingleMolecule *sm, *sm2;
	CxDVector3 vec;
//	FILE *a;

	m = (CMolecule*)g_oaMolecules[g_iFixMol];
	sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[refmol]];
	del = new double[g_oaSingleMolecules.GetSize()];
	best = new int[g_oaSingleMolecules.GetSize()];
	done = new bool[g_oaSingleMolecules.GetSize()];

//	printf("Suche Nachbarschaft von Molekuel %d...\n",refmol+1);

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		if (nb->m_bDistMode && (nb->m_faMolDist[z] == 0))
		{
			nb->m_waScanNeighborCount[z] = 0;
			continue;
		}
		m2 = (CMolecule*)g_oaMolecules[z];
//		printf("*** Molekuel %d: %s\n",z+1,m2->Name);
//		for (z2=0;z2<=((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
//			del[z2] = 99999.0f;
		for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
		{
			del[z2] = 99999.0f;
			if ((z == g_iFixMol) && (z2 == refmol))
				continue;
//			printf("  - Vertreter %d\n",z2+1);
			sm2 = (CSingleMolecule*)g_oaSingleMolecules[m2->m_laSingleMolIndex[z2]];
			for (z3=0;z3<nb->m_waRefElements.GetSize();z3++)
				for (z4=0;z4<((CxIntArray*)nb->m_oaNbElements[z])->GetSize();z4++)
				{
					vec = m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[nb->m_waRefElements[z3]])->GetAt(nb->m_waRefAtoms[z3])] - m_vaCoords[((CxIntArray*)sm2->m_oaAtomOffset[((CxIntArray*)nb->m_oaNbElements[z])->GetAt(z4)])->GetAt(((CxIntArray*)nb->m_oaNbAtoms[z])->GetAt(z4))];
					if (g_bFold)
					{
						while (vec[0] >= g_fBoxX/2) vec[0] -= g_fBoxX;
						while (vec[0] < -g_fBoxX/2) vec[0] += g_fBoxX;
						while (vec[1] >= g_fBoxY/2) vec[1] -= g_fBoxY;
						while (vec[1] < -g_fBoxY/2) vec[1] += g_fBoxY;
						while (vec[2] >= g_fBoxZ/2) vec[2] -= g_fBoxZ;
						while (vec[2] < -g_fBoxZ/2) vec[2] += g_fBoxZ;
					}
					d = vec.GetLength();
//					printf("    = m1 %d (%d|%d); m2 %d (%d|%d); Dist = %f\n",z3+1,m->m_waNbElements[z3]+1,m->m_waNbAtoms[z3]+1,z4+1,m2->m_waNbElements[z4]+1,m2->m_waNbAtoms[z4]+1,d);
					if (d < del[z2])
						del[z2] = d;
				}
//			printf("  - Finaler Abstand: %f\n",del[z2]);
		}

		for (z0=0;z0<nb->m_pAF->m_iResolution;z0++)
		{
			c = 0;
			for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
				done[z2] = false;
			while (true)
			{
				d = 999999.0f;
				b = -1;
				for (z3=0;z3<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z3++)
				{
					if (done[z3])
						continue;
					if (del[z3] < d)
					{
						d = del[z3];
						b = z3;
					}
				}
				if (d > nb->m_pAF->m_fMaxVal*z0/nb->m_pAF->m_iResolution)
					break;
				best[c] = b;
				done[b] = true;
				c++;
//				if (c >= ((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize())
//					break;
			}
			nb->m_pAF->AddToBin(nb->m_pAF->m_fMaxVal*z0/nb->m_pAF->m_iResolution,(double)c);
//			if (g_iSteps == 1)
//				mfprintf(a,"%f;%f\n",nb->m_pAF->m_fMaxVal*z0/nb->m_pAF->m_iResolution,(double)c);
		}
//		if (g_iSteps == 1)
//			fclose(a);
	}
	delete del;
	delete best;
	BTOUT; 
}*/


int CTimeStep::REC_UniteMolecules(CSingleMolecule *sm, int i0, int depth)
{
	int z, z2, n;
	CMolAtom *m;
	FILE *a;


	if (depth > 1000) {

		eprintf("CTimeStep::REC_UniteMolecules(): Error: depth > 1000.\n");
		eprintf("Unwrapping of periodic bonds failed.\n\n");
		mprintf("Is your molecule periodic across the cell boundaries? In this case, you need to\n");
		mprintf("answer the \"polymer\" question after the molecule recognition with \"yes\".\n\n");
		mprintf("Writing current trajectory snapshot to \"unwrap_failed.xyz\" for debugging purposes...\n");

		a = OpenFileWrite( "unwrap_failed.xyz", true );

		mfprintf( a, "%lu\nComment\n", m_iGesAtomCount );

		for (z=0;z<(int)m_iGesAtomCount;z++)
			mfprintf(
				a,
				"%-4s  %16.10f  %16.10f  %16.10f\n",
				((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_sLabel,
				m_vaCoords_Original[z][0]/100.0,
				m_vaCoords_Original[z][1]/100.0,
				m_vaCoords_Original[z][2]/100.0
			);

		fclose(a);

		mprintf("Done. Aborting execution.\n\n");

		abort();
	}


	n = 0;
	m = (CMolAtom*)sm->m_oaMolAtoms[i0];

/*	if (g_bVerbose)
	{
		mprintf("# ");
		for (z2=0;z2<depth;z2++)
			mprintf(". ");
		mprintf(">>> REC_UniteMolecules MolAtom=%d, Offset=%d.\n",i0,m->m_iOffset);
	}*/

	g_pUniteTemp[m->m_iOffset] = true;

	for (z=0;z<m->m_oaBonds.GetSize();z++)
	{
		if (!g_pUniteTemp[((CMolAtom*)m->m_oaBonds[z])->m_iOffset])
		{
			if (MirrorBond(m->m_iOffset,((CMolAtom*)m->m_oaBonds[z])->m_iOffset))
			{
				if (g_bVerbose)
				{
					mprintf("# ");
					for (z2=0;z2<depth;z2++)
						mprintf(". ");
					mprintf("    Bond %s(%d) <--> %s(%d) unwrapped.\n",(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[m->m_iOffset]])->m_sName,m->m_iOffset+1,(const char*)((CAtom*)g_oaAtoms[g_baAtomIndex[((CMolAtom*)m->m_oaBonds[z])->m_iOffset]])->m_sName,((CMolAtom*)m->m_oaBonds[z])->m_iOffset+1);
				}
				n++;
			}
			n += REC_UniteMolecules(sm,((CMolAtom*)m->m_oaBonds[z])->m_iMolAtomNumber,depth+1);
		}
	}
/*	if (g_bVerbose)
	{
		mprintf("# ");
		for (z2=0;z2<depth;z2++)
			mprintf(". ");
		mprintf("<<< REC_UniteMolecules MolAtom=%d, Offset=%d.\n",i0,m->m_iOffset);
	}*/
	return n;
}


void CTimeStep::UniteMolecules(bool verbose)
{
	BTIN;
	int z, z2, n;
	CMolecule *m;
	CSingleMolecule *sm;

	for (z=0;z<g_iGesAtomCount;z++)
		g_pUniteTemp[z] = false;
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];

		if (m->m_bPseudo)
			continue;

		if (m->m_bPolymer)
		{
			if (verbose)
				mprintf("      Skipping molecule %s: Is a polymer of infinite extent.\n",m->m_sName);

			continue;
		}

		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];

//			if (g_bVerbose)
//				mprintf("  # UniteMolecules molecule %s (%d)\n",m->m_sName,z2+1);

			n = REC_UniteMolecules(sm,0,0);

			if ((n != 0) && verbose)
				mprintf("      - Molecule %s[%d] united, %d bonds unwrapped.\n",m->m_sName,z2+1,n);

			if ((n != 0) && (!verbose) && g_bVerbose)
				mprintf("\n  # UniteMolecules: Molecule %s[%d] united, %d bonds unwrapped.",m->m_sName,z2+1,n);
	
		}
	}
	BTOUT; 
}



bool CTimeStep::ReadTimestep(FILE *a, bool needinfo) {

	int z;


	switch(g_iTrajFormat) {

		case 0:
			if (!ReadXYZ(a,needinfo,&m_vaCoords))
				return false;
			break;

		case 1:
			if (!ReadPDB(a,needinfo,&m_vaCoords))
				return false;
			break;

		case 2:
			if (!ReadMol2(a,needinfo))
				return false;
			break;

		case 3:
			if (!ReadLAMMPS(a,needinfo))
				return false;
			break;

		case 4:
			if (!ReadDLPOLY(a,needinfo))
				return false;
			break;

		case 5:
			if (!ReadCube(a,needinfo))
				return false;
			break;

		case 6:
			if (!ReadAmber(a,needinfo))
				return false;
			break;

		case 7:
			if (!ReadBQB(g_pBQBFile,needinfo))
				return false;
			break;

		case 8:
			if (!ReadGRO(a,needinfo))
				return false;
			break;

		case 9:
			if (!ReadVoronoi(a,needinfo))
				return false;
			break;

		case 10:
			if (!ReadDCD(a,needinfo))
				return false;
			break;
	}

	for (z=0;z<(int)m_iGesAtomCount;z++) {
		if ((m_vaCoords[z][0] < -1.0e6) || (m_vaCoords[z][0] > 1.0e6) || (m_vaCoords[z][1] < -1.0e6) || (m_vaCoords[z][1] > 1.0e6) || (m_vaCoords[z][2] < -1.0e6) || (m_vaCoords[z][2] > 1.0e6)) {
			eprintf("\nWarning:");
			mprintf(" The position of atom %d in the trajectory is outside of the range -10^6 .. 10^6 pm:\n",z+1);
			mprintf("           ( %.10G | %.10G | %.10G ) pm\n",m_vaCoords[z][0],m_vaCoords[z][1],m_vaCoords[z][2]);
			mprintf("         This likely means that your simulation \"exploded\". Suppressing this warning for further atoms.\n");
			break;
		}
	}

	if (g_bDoubleBox)
		DoubleBox();

	BTOUT; 
	return true;
}



bool CTimeStep::ReadTimestep(CxMemFile *file) {
	if (!ReadCube(file))
		return false;
	if (g_bDoubleBox)
		DoubleBox();
	return true;
}


bool CTimeStep::SkipTimestep(FILE *a)
{
	BTIN;
	switch(g_iTrajFormat)
	{
		case 0:
			if (!SkipXYZ(a))
				return false;
			break;

		case 1:
			if (!SkipPDB(a))
				return false;
			break;

		case 2:
			if (!SkipMol2(a))
				return false;
			break;

		case 3:
			if (!SkipLAMMPS(a))
				return false;
			break;

		case 4:
			if (!SkipDLPOLY(a))
				return false;
			break;

		case 5:
			if (!SkipCube(a))
				return false;
			break;

		case 6:
			if (!SkipAmber(a))
				return false;
			break;

		case 7:
			if (!SkipBQB())
				return false;
			break;

		case 8:
			if (!SkipGRO(a))
				return false;
			break;

		case 9:
			if (!SkipVoronoi(a))
				return false;
			break;

		case 10:
			if (!SkipDCD(a))
				return false;
			break;
	}
	BTOUT; 
	return true;
}


bool CTimeStep::ReadTimestepVel(FILE *a)
{
	BTIN;
	char buf[256], *p, *q;
	int z, tc;

	buf[0] = 0;
	(void)!fgets(buf,256,a);
	if (feof(a))
	{
		BTOUT; 
		return false;
	}
	if (strlen(buf) > 0)
		buf[strlen(buf)-1] = 0;
	tc = atoi(buf);
	if (tc == 0)
	{
		BTOUT; 
		return false;
	}
	m_vaVelocities.SetSize(tc);
	buf[0] = 0;
	(void)!fgets(buf,256,a); // Zeitschritt - egal hier
	if (strlen(buf) > 0)
		buf[strlen(buf)-1] = 0;
	for (z=0;z<tc;z++) // den ersten Zeitschritt einlesen
	{
		buf[0] = 0;
		(void)!fgets(buf,256,a);
		if (feof(a))
		{
			BTOUT; 
			return false;
		}
		if (strlen(buf) > 0)
			buf[strlen(buf)-1] = 0;
		q = buf;
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p == NULL)
		{
			eprintf("CTimeStep::ReadTimestepVel(): Error 1. %d, \"%s\"\n",z+1,buf);
			continue;
		}
		while (isdigit(*(p-1)) && (p > buf))
			p--;
		if (p == buf)
		{
			eprintf("CTimeStep::ReadTimestepVel(): No Atom laben found. %d, \"%s\"\n",z+1,buf);
			continue;
		}
		*p = 0;
		p++;
		q = strchr(p,' ');
		if (q == NULL)
		{
			eprintf("CTimeStep::ReadTimestepVel(): Error 2. %d, \"%s\"\n",z+1,p);
			continue;
		}
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p == NULL)
		{
			eprintf("CTimeStep::ReadTimestepVel(): Error 3. %d, \"%s\"\n",z+1,q);
			continue;
		}
		*p = 0;
		// Convert to pm/ps = m/s as used in TRAVIS
		m_vaVelocities[z][0] = atof(q) * g_fVelocityConversion;
		q = p+1;
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p == NULL)
		{
			eprintf("CTimeStep::ReadTimestepVel(): Error 4. %d, \"%s\"\n",z+1,q);
			continue;
		}
		*p = 0;
		m_vaVelocities[z][1] = atof(q) * g_fVelocityConversion;
		q = p+1;
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p != NULL)
			*p = 0;
		m_vaVelocities[z][2] = atof(q) * g_fVelocityConversion;
	}
	if (g_bDoubleBox)
		DoubleBoxVelocity();
	BTOUT; 
	return true;
}


bool CTimeStep::ReadTimestepForce(FILE *a)
{
	BTIN;
	char buf[256], *p, *q;
	int z, tc;

	buf[0] = 0;
	(void)!fgets(buf,256,a);
	if (feof(a))
	{
		BTOUT; 
		return false;
	}
	if (strlen(buf) > 0)
		buf[strlen(buf)-1] = 0;
	tc = atoi(buf);
	if (tc == 0)
	{
		BTOUT; 
		return false;
	}
	m_vaForces.SetSize(tc);
	buf[0] = 0;
	(void)!fgets(buf,256,a); // Zeitschritt - egal hier
	if (strlen(buf) > 0)
		buf[strlen(buf)-1] = 0;
	for (z=0;z<tc;z++) // den ersten Zeitschritt einlesen
	{
		buf[0] = 0;
		(void)!fgets(buf,256,a);
		if (feof(a))
		{
			BTOUT; 
			return false;
		}
		if (strlen(buf) > 0)
			buf[strlen(buf)-1] = 0;
		q = buf;
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p == NULL)
		{
			eprintf("CTimeStep::ReadTimestepForce(): Error 1. %d, \"%s\"\n",z+1,buf);
			continue;
		}
		while (isdigit(*(p-1)) && (p > buf))
			p--;
		if (p == buf)
		{
			eprintf("CTimeStep::ReadTimestepForce(): No Atom laben found. %d, \"%s\"\n",z+1,buf);
			continue;
		}
		*p = 0;
		p++;
		q = strchr(p,' ');
		if (q == NULL)
		{
			eprintf("CTimeStep::ReadTimestepForce(): Error 2. %d, \"%s\"\n",z+1,p);
			continue;
		}
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p == NULL)
		{
			eprintf("CTimeStep::ReadTimestepForce(): Error 3. %d, \"%s\"\n",z+1,q);
			continue;
		}
		*p = 0;
		m_vaForces[z][0] = atof(q);
		q = p+1;
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p == NULL)
		{
			eprintf("CTimeStep::ReadTimestepForce(): Error 4. %d, \"%s\"\n",z+1,q);
			continue;
		}
		*p = 0;
		m_vaForces[z][1] = atof(q);
		q = p+1;
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p != NULL)
			*p = 0;
		m_vaForces[z][2] = atof(q);
	}
	if (g_bDoubleBox)
		DoubleBoxForce();
	BTOUT; 
	return true;
}


bool CTimeStep::SkipXYZ(FILE *a)
{
	BTIN;
	char buf[256];
	int z, tc;

//	mprintf("*** Skip Anfang ***\n");
	buf[0] = 0;
	(void)fgets_bin(buf,256,a);
	if (feof(a))
	{
		BTOUT; 
		return false;
	}
	if (strlen(buf)==0)
	{
		BTOUT; 
		return false;
	}
	buf[strlen(buf)-1] = 0;
	tc = atoi(buf);
//	mprintf("SkipA: \"%s\".\n",buf);
	if (tc == 0)
	{
		BTOUT; 
		return false;
	}
//	buf[0] = 0;
	(void)fgets_bin(buf,256,a); // Zeitschritt - egal hier
//	mprintf("SkipB: \"%s\".\n",buf);
	for (z=0;z<tc;z++) // den ersten Zeitschritt einlesen
	{
//		buf[0] = 0;
		//fgets_bin(buf,256,a);
		(void)!fgets(buf,256,a);
//		mprintf("SkipC%d: \"%s\".\n",z,buf);
	}
	if (feof(a))
	{
		BTOUT; 
		return false;
	}
//	mprintf("*** Skip Ende ***\n");
	BTOUT; 
	return true;
}


void CTimeStep::CopyFrom(CTimeStep *t)
{
	BTIN;
	int z;
	char *p;
	CxDVector3 v;

	if (t == NULL) {
		eprintf("CTimeStep::CopyFrom(): Error: Argument == NULL.\n");
		abort();
	}

	m_iGesAtomCount = t->m_iGesAtomCount;
	m_vaCoords.CopyFrom(&t->m_vaCoords);
	if (g_bKeepUnfoldedCoords)
		m_vaCoords_Unfolded.CopyFrom(&t->m_vaCoords_Unfolded);
	m_vaForces.CopyFrom(&t->m_vaForces);
	m_vaVelocities.CopyFrom(&t->m_vaVelocities);
	if (t->m_paLabels.GetSize() != 0)
	{
		for (z=0;z<m_paLabels.GetSize();z++)
			delete[] (char*)m_paLabels[z];
		m_paLabels.RemoveAll();
		for (z=0;z<t->m_paLabels.GetSize();z++)
		{
			try { p = new char[strlen((char*)t->m_paLabels[z])+1]; } catch(...) { p = NULL; }
			if (p == NULL) NewException((double)(strlen((char*)t->m_paLabels[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(p,(char*)t->m_paLabels[z]);
			m_paLabels.Add(p);
		}
	}
	if (t->m_pComment != NULL)
	{
		if (m_pComment == NULL)
		{
			try { m_pComment = new char[256]; } catch(...) { m_pComment = NULL; }
			if (m_pComment == NULL) NewException((double)256*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}
		strcpy(m_pComment,t->m_pComment);
	}
	if (t->m_pVolumetricData != NULL)
	{
		if (m_pVolumetricData != NULL)
			delete m_pVolumetricData;
		m_pVolumetricData = new C3DF<VORI_FLOAT>();
		m_pVolumetricData->CopyFrom(t->m_pVolumetricData);
	}
	if (t->m_pVolumetricDataTimeDev != NULL) {
		if (m_pVolumetricDataTimeDev != NULL)
			delete m_pVolumetricDataTimeDev;
		m_pVolumetricDataTimeDev = new C3DF<VORI_FLOAT>();
		m_pVolumetricDataTimeDev->CopyFrom(t->m_pVolumetricDataTimeDev);
	}
	if (t->m_pCurrentDensity != NULL) {
		if (m_pCurrentDensity != NULL)
			delete m_pCurrentDensity;
		m_pCurrentDensity = new CxDoubleArray();
		m_pCurrentDensity->CopyFrom(t->m_pCurrentDensity);
	}
	BTOUT; 
}

	  
long CTimeStep::ExtractNumber(int i)
{
	BXIN;
	char *p, *q, buf[20];
	int z;
	long l;

	if (m_pComment == NULL)
		return -1;
	p = m_pComment;
	for (z=0;z<i;z++)
	{
		while ((!isdigit(*p)) && (*p != 0))
			p++;
		if (*p == 0)
		{
			BXOUT;
			return -1;
		}
		q = p+1;
		while ((isdigit(*q)) && (*q != 0))
			q++;
		if (*q == 0)
		{
			BXOUT;
			return -1;
		}
		p = q+1;
	}
	while ((!isdigit(*p)) && (*p != 0))
		p++;
	if (*p == 0)
	{
		BXOUT;
		return -1;
	}
	q = p;
	while ((isdigit(*q)) && (*q != 0))
		q++;
	if (q == p)
	{
		BXOUT;
		return -1;
	}
	if (q-p >= 256)
	{
		eprintf("Internal Error in ExtractNumber(): %ld >= 256.\n",(long)(q-p));
		return 0;
	}
	memcpy(buf,p,q-p);
	buf[q-p] = 0;
	l = atoi(buf);
	BXOUT;
	return l;
}


int CTimeStep::GetCommentNumberCount()
{
	BXIN;
	char *p, *q;
	int z;

	if (m_pComment == NULL)
		return 0;
	p = m_pComment;
	z = 0;
	while (true)
	{
		while ((!isdigit(*p)) && (*p != 0))
			p++;
		if (*p == 0)
		{
//			mprintf("GetCommentNumberCount 1: %d (p=\"%s\", q=\"%s\")\n",z,p,q);
			BXOUT;
			return z;
		}
		q = p;
		while ((isdigit(*q)) && (*q != 0))
			q++;
		if (q == p)
		{
//			mprintf("GetCommentNumberCount 2: %d (p=\"%s\", q=\"%s\")\n",z,p,q);
			BXOUT;
			return z;
		}
		z++;
		if (*q == 0)
		{
//			mprintf("GetCommentNumberCount 3: %d (p=\"%s\", q=\"%s\")\n",z,p,q);
			BXOUT;
			return z;
		}
		p = q+1;
	}
	return 0; // Never happens
}


bool CTimeStep::ScanWannier(bool verbose)
{
	BTIN;
	int z, z2, z3;
	double td, d, dx, dy, dz;
	CMolecule *m;

	if (g_bVerbose)
	{
		mprintf(WHITE,"\n*** ScanWannier ***\n\n");
		verbose = true;
	}

	for (z=0;z<g_oaSingleMolecules.GetSize();z++)
		((CSingleMolecule*)g_oaSingleMolecules[z])->m_laWannier.RemoveAll_KeepSize();

	for (z=0;z<g_iGesAtomCount;z++)
	{
		if (g_baAtomIndex[z] != g_iWannierAtomType)
			continue;
		d = 9999.0;
		z3 = -1;
		for (z2=0;z2<g_iGesAtomCount;z2++)
		{
			if (g_baAtomIndex[z2] == g_iWannierAtomType)
				continue;
			dx = m_vaCoords[z][0]-m_vaCoords[z2][0];
			dy = m_vaCoords[z][1]-m_vaCoords[z2][1];
			dz = m_vaCoords[z][2]-m_vaCoords[z2][2];
			
			if (g_bBoxNonOrtho)
			{
				CxDVector3 v = g_mBoxToOrtho * CxDVector3(dx, dy, dz);
				while (v[0] > 0.5)
					v[0] -= 1.0;
				while (v[0] <= -0.5)
					v[0] += 1.0;
				while (v[1] > 0.5)
					v[1] -= 1.0;
				while (v[1] <= -0.5)
					v[1] += 1.0;
				while (v[2] > 0.5)
					v[2] -= 1.0;
				while (v[2] <= -0.5)
					v[2] += 1.0;
				v = g_mBoxFromOrtho * v;
				dx = v[0];
				dy = v[1];
				dz = v[2];
			} else {
				if (g_bPeriodicX)
				{
					while (dx > g_fBoxX/2.0) dx -= g_fBoxX;
					while (dx < -g_fBoxX/2.0) dx += g_fBoxX;
				}
				if (g_bPeriodicY)
				{
					while (dy > g_fBoxY/2.0) dy -= g_fBoxY;
					while (dy < -g_fBoxY/2.0) dy += g_fBoxY;
				}
				if (g_bPeriodicZ)
				{
					while (dz > g_fBoxZ/2.0) dz -= g_fBoxZ;
					while (dz < -g_fBoxZ/2.0) dz += g_fBoxZ;
				}
			}
			
			td = sqrt(dx*dx+dy*dy+dz*dz);
			if (td < d)
			{
				d = td;
				z3 = z2;
			}
		}
		if (z3 == -1) {
			eprintf("CTimeStep::ScanWannier(): Internal error.\n");
			abort();
		}
		
		if (g_bBoxNonOrtho) {
			CxDVector3 v = g_mBoxToOrtho * (m_vaCoords[z] - m_vaCoords[z3]);
			while (v[0] > 0.5)
				v[0] -= 1.0;
			while (v[0] <= -0.5)
				v[0] += 1.0;
			while (v[1] > 0.5)
				v[1] -= 1.0;
			while (v[1] <= -0.5)
				v[1] += 1.0;
			while (v[2] > 0.5)
				v[2] -= 1.0;
			while (v[2] <= -0.5)
				v[2] += 1.0;
			v = g_mBoxFromOrtho * v;
			m_vaCoords[z] = v + m_vaCoords[z3];
		} else {
			if (g_bPeriodicX)
			{
				while (m_vaCoords[z][0]-m_vaCoords[z3][0] > g_fBoxX/2.0) m_vaCoords[z][0] -= g_fBoxX;
				while (m_vaCoords[z][0]-m_vaCoords[z3][0] < -g_fBoxX/2.0) m_vaCoords[z][0] += g_fBoxX;
			}
			if (g_bPeriodicY)
			{
				while (m_vaCoords[z][1]-m_vaCoords[z3][1] > g_fBoxY/2.0) m_vaCoords[z][1] -= g_fBoxY;
				while (m_vaCoords[z][1]-m_vaCoords[z3][1] < -g_fBoxY/2.0) m_vaCoords[z][1] += g_fBoxY;
			}
			if (g_bPeriodicZ)
			{
				while (m_vaCoords[z][2]-m_vaCoords[z3][2] > g_fBoxZ/2.0) m_vaCoords[z][2] -= g_fBoxZ;
				while (m_vaCoords[z][2]-m_vaCoords[z3][2] < -g_fBoxZ/2.0) m_vaCoords[z][2] += g_fBoxZ;
			}
		}
		
		z2 = g_laAtomSMIndex[z3];
		if (d > 105.0)
		{
			eprintf("Step %lu: Wannier center at offset %d too far away from any atom (closest atom is %s[%d] %s%d, %.0f pm).\n",g_iSteps,z+1,((CMolecule*)g_oaMolecules[((CSingleMolecule*)g_oaSingleMolecules[z2])->m_iMolType])->m_sName,((CSingleMolecule*)g_oaSingleMolecules[z2])->m_iMolSMIndex+1,(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[z3]])->m_sName,g_waAtomMolNumber[z3]+1,d);
		} else if (verbose)
			mprintf("  - Wannier center %d belongs to %s[%d] %s%d (%.0f pm).\n",z+1,((CMolecule*)g_oaMolecules[((CSingleMolecule*)g_oaSingleMolecules[z2])->m_iMolType])->m_sName,((CSingleMolecule*)g_oaSingleMolecules[z2])->m_iMolSMIndex+1,(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[z3]])->m_sName,g_waAtomMolNumber[z3]+1,d);
/*		if (z2 == -1)
		{
			eprintf("Wannier center %d: Atom %d does not belong to any molecule.\n",z+1,z3+1);
			continue;
		}*/
//		mprintf("Wannier Center %d hat Abstand %.2f zu Atom %d in SM %d.\n",z+1,d,z3+1,z2+1);
		((CSingleMolecule*)g_oaSingleMolecules[z2])->m_laWannier.Add(z);
	}
	if (verbose)
		mprintf("\n\n");
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		if (m->m_iWannierCount == 0 && !m->m_bPseudo)
		{
			m->m_iWannierCount = ((CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_laWannier.GetSize();
			td = 0;
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++)
			{
				if (m->m_baAtomIndex[z2] == g_iWannierAtomType)
					continue;
				if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
					continue;
				td += ((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_fCharge * m->m_waAtomCount[z2];
			}
			m->m_fCharge = td-m->m_iWannierCount*g_fWannierCharge;
			mprintf("  - Molecule %s contains %d wannier centers. Total charge is %.2f - %.2f = %.2f\n",m->m_sName,m->m_iWannierCount,td,m->m_iWannierCount*g_fWannierCharge,m->m_fCharge);
			if (m->m_fCharge > 5.0)
			{
				eprintf("\n    This molecular charge seems to be too high.\n\n");
				if (!AskYesNo("    Do you want to continue (y/n)? [yes] ",true))
					return false;
			}
		}
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			if (m->m_iWannierCount != ((CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]])->m_laWannier.GetSize())
				eprintf("Step %lu: Molecule %s[%d] contains %d instead of %d wannier centers.\n",g_iSteps,m->m_sName,z2+1,((CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]])->m_laWannier.GetSize(),m->m_iWannierCount);
		}
	}
	BTOUT;
	return true;
}


void CTimeStep::CalcMagneticDipoles() {

	const int BUF_SIZE = 1024;
	int i;


	for (i = 0; i < g_oaSingleMolecules.GetSize(); i++) {
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
		if (g_bLoadMagneticDipoleRestart) {
			int j;
			for (j = 0; j < 3; j++) {
				double val;
				if (fread(&val, sizeof(double), 1, g_fMagneticDipoleRestartFile) < 1) {
					eprintf("Could not read from magnetic moment restart file. Setting value to 0.0.\n");
					val = 0.0;
				}
				sm->m_vMagneticDipole[j] = val;
			}
		} else {
			sm->m_vMagneticDipole = CxDVector3(0.0, 0.0, 0.0);
			sm->m_vCurrent = CxDVector3(0.0, 0.0, 0.0);
		}
	}
	
	for (i = 0; i < g_oaMolecules.GetSize(); i++) {
		CMolecule *m = (CMolecule *)g_oaMolecules[i];
		int j;
		for (j = 0; j < m->m_laSingleMolIndex.GetSize(); j++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[j]];
			if (m->m_iMagneticDipoleMode == 0) {
				sm->m_vMagneticDipole = CxDVector3(0.0, 0.0, 0.0);
			} else if (m->m_iMagneticDipoleMode == 1) {
				CxDVector3 ref(0.0, 0.0, 0.0);
				ref = m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];
				int k;
				for (k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
					if (m->m_baAtomIndex[k] == g_iVirtAtomType)
						continue;
					if (m->m_baAtomIndex[k] == g_iWannierAtomType)
						continue;
					CAtom *a = (CAtom *)g_oaAtoms[m->m_baAtomIndex[k]];
					int l;
					for (l = 0; l < m->m_waAtomCount[k]; l++) {
						sm->m_vMagneticDipole += 0.5 * a->m_fCharge * CrossP(m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref, m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)]);
// 						mprintf(RED, "c2: ( %14.10G | %14.10G | %14.10G )\n", (m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref)[0], (m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref)[1], (m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref)[2]);
// 						mprintf(RED, "v2: ( %14.10G | %14.10G | %14.10G )\n", m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][0], m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][1], m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][2]);
					}
				}
				for (k = 0; k < sm->m_laWannier.GetSize(); k++) {
					sm->m_vMagneticDipole -= g_fWannierCharge * CrossP(m_vaCoords[sm->m_laWannier[k]] - ref, m_vaVelocities[sm->m_laWannier[k]]);
// 					mprintf(RED, "c2: ( %14.10G | %14.10G | %14.10G )\n", (m_vaCoords[sm->m_laWannier[k]] - ref)[0], (m_vaCoords[sm->m_laWannier[k]] - ref)[1], (m_vaCoords[sm->m_laWannier[k]] - ref)[2]);
// 					mprintf(RED, "v2: ( %14.10G | %14.10G | %14.10G )\n", m_vaVelocities[sm->m_laWannier[k]][0], m_vaVelocities[sm->m_laWannier[k]][1], m_vaVelocities[sm->m_laWannier[k]][2]);
				}
				sm->m_vMagneticDipole *= MAG_EPMMS2MB;
// 				mprintf(RED, "M: ( %14.10G | %14.10G | %14.10G )\n", sm->m_vMagneticDipole[0], sm->m_vMagneticDipole[1], sm->m_vMagneticDipole[2]);
			} else if (m->m_iMagneticDipoleMode == 2) {
				if (j == 0) {
					char buf[BUF_SIZE];
					double dipole[3];
					if (fgets(buf, BUF_SIZE, m->m_pMagneticDipoleFile) == NULL || sscanf(buf, "%lf %lf %lf", &dipole[0], &dipole[1], &dipole[2]) < 3)
						eprintf("Could not read from dipole file. Setting magnetic moment to ( 0.0 | 0.0 | 0.0 ).\n");
					else
						sm->m_vMagneticDipole = CxDVector3(dipole[0], dipole[1], dipole[2]);
				} else {
					sm->m_vMagneticDipole = ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vMagneticDipole;
				}
			} else if (m->m_iMagneticDipoleMode == 3) {
				if (j == 0) {
					char comment[BUF_SIZE];
					strncpy(comment, m_pComment, BUF_SIZE-1 );
					comment[BUF_SIZE-1] = '\0';
					const char delim[] = "\t ";
					char *tok = strtok(comment, delim);
					int num = 0;
					bool ok[3] = { false, false, false };
					while (tok != NULL) {
						double f;
						if (sscanf(tok, "%lf", &f) == 1) {
							if (num == m->m_iMagneticDipoleCommentIndex[0]) {
								sm->m_vMagneticDipole[0] = f;
								ok[0] = true;
							}
							if (num == m->m_iMagneticDipoleCommentIndex[1]) {
								sm->m_vMagneticDipole[1] = f;
								ok[1] = true;
							}
							if (num == m->m_iMagneticDipoleCommentIndex[2]) {
								sm->m_vMagneticDipole[2] = f;
								ok[2] = true;
							}
							num++;
						}
						tok = strtok(NULL, delim);
					}
					if (!(ok[0] && ok[1] && ok[2])) {
						eprintf("Could not read from comment line. Setting magnetic moment to ( 0.0 | 0.0 | 0.0 ).\n");
						sm->m_vMagneticDipole = CxDVector3(0.0, 0.0, 0.0);
					}
				} else {
					sm->m_vMagneticDipole = ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vMagneticDipole;
				}
			} else if (m->m_iMagneticDipoleMode == 4) {
				CxDVector3 ref(0.0, 0.0, 0.0);
				ref = m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];
				int k;
				for (k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
					if (m->m_baAtomIndex[k] == g_iVirtAtomType)
						continue;
					int l;
					for (l = 0; l < m->m_waAtomCount[k]; l++)  {
						sm->m_vMagneticDipole += ((CxDoubleArray *)m->m_oaCharges[k])->GetAt(l) * CrossP(m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref, m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)]);
					}
					sm->m_vMagneticDipole *= MAG_EPMMS2MB;
				}
			} else if (m->m_iMagneticDipoleMode == 5) {
				CxDVector3 ref(0.0, 0.0, 0.0);
				ref = m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];
				int k;
				for (k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
					if (m->m_baAtomIndex[k] == g_iVirtAtomType)
						continue;
					int l;
					for (l = 0; l < m->m_waAtomCount[k]; l++) {
						sm->m_vMagneticDipole += m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] * CrossP(m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref, m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)]);
					}
					sm->m_vMagneticDipole *= MAG_EPMMS2MB;
				}
			} else if (m->m_iMagneticDipoleMode == 6) {
				CxDVector3 ref(0.0, 0.0, 0.0);
				ref = m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];
				int k;
				for (k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
					if (m->m_baAtomIndex[k] == g_iVirtAtomType)
						continue;
					int l;
					for (l = 0; l < m->m_waAtomCount[k]; l++) {
						sm->m_vMagneticDipole += m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] * CrossP(m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref, m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)]);
					}
					sm->m_vMagneticDipole *= MAG_EPMMS2MB;
				}
			} else if (m->m_iMagneticDipoleMode == 7) {

				CxDVector3 ref(0.0, 0.0, 0.0);

				if (!g_bDipoleRefFixed)
					ref = m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];

				int k;
// 				double totalCharge = 0.0f;
				for (k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
					if (m->m_baAtomIndex[k] == g_iVirtAtomType)
						continue;
					CAtom *a = (CAtom *)g_oaAtoms[m->m_baAtomIndex[k]];
					int l;
					for (l = 0; l < m->m_waAtomCount[k]; l++) {

			/*			if (j == 1)
							mprintf("\n#%d.%d: += (%10g|%10g|%10g) + 0.5*CrossP[ (%10g|%10g|%10g)-(%10g|%10g|%10g), (%10g|%10g|%10g) ] + 0.5*%10g*CrossP[ (%10g|%10g|%10g)-(%10g|%10g|%10g), (%10g|%10g|%10g) ]",
								k+1,
								l+1,
								m_magneticDipoleMoments[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][0],
								m_magneticDipoleMoments[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][1],
								m_magneticDipoleMoments[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][2],
								m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][0],
								m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][1],
								m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][2],
								ref[0],
								ref[1],
								ref[2],
								m_totalCurrents[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][0],
								m_totalCurrents[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][1],
								m_totalCurrents[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][2],
								a->m_fCharge,
								m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][0],
								m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][1],
								m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][2],
								ref[0],
								ref[1],
								ref[2],
								m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][0],
								m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][1],
								m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][2]);*/

						sm->m_vMagneticDipole += m_magneticDipoleMoments[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] + 0.5 * CrossP(m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref, m_totalCurrents[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)]) + 0.5 * a->m_fCharge * CrossP(m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref, m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)]) * MAG_EPMMS2MB;
						sm->m_vCurrent += m_totalCurrents[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)];
// 						totalCharge += m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)];
// 						sm->m_vMagneticDipole += 0.5f * a->m_fCharge * CrossP(m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref, m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)]) * 1.727599e-8f;
// 						sm->m_vMagneticDipole += m_magneticDipoleMoments[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - 0.5f * CrossP(m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref, m_totalCurrents[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)]);
// 						sm->m_vMagneticDipole += m_magneticDipoleMoments[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)];
					}
				}

				if (!g_bDipoleRefFixed)
					sm->m_vMagneticDipole -= 0.5 * CrossP(sm->m_vDipole, m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)]) * MAG_EPMMS2MB / DIP_EPM2DEBYE;

	/*			if (j == 1) {
					mprintf("\n-= CrossP[ (%10g|%10g|%10g), (%10g|%10g|%10g) ]",
						sm->m_vDipole[0],
						sm->m_vDipole[1],
						sm->m_vDipole[2],
						m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)][0],
						m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)][1],
						m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)][2]);

					mprintf("\n--> (%10g|%10g|%10g).",sm->m_vMagneticDipole[0],sm->m_vMagneticDipole[1],sm->m_vMagneticDipole[2]);
				}*/

// 				sm->m_vMagneticDipole += 0.5f * totalCharge * CrossP(ref, m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)]) * (double)MAG_EPMMS2MB;
// 				mprintf(GREEN, "%f %f %f\n", sm->m_vDipole[0], sm->m_vDipole[1], sm->m_vDipole[2]);
// 				mprintf(GREEN, "%f %f %f\n", m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)][0], m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)][1], m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)][2]);
// 				mprintf(GREEN, "%f %f %f\n", (0.5f * CrossP(sm->m_vDipole, m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)]) * 3.596761e-7f)[0], (0.5f * CrossP(sm->m_vDipole, m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)]) * 3.596761e-7f)[1], (0.5f * CrossP(sm->m_vDipole, m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)]) * 3.596761e-7f)[2]);

			} else if (m->m_iMagneticDipoleMode == 8) {

				eprintf("CTimeStep::CalcMagneticDipoles(): Mode 8 not implemented.\n");
				abort();

			} else {

				eprintf("CTimeStep::CalcMagneticDipoles(): Internal error.\n");
				abort();

			}
		}
	}
}


void CTimeStep::CalcDipoles(bool verbose)
{
	const int BUF_SIZE = 1024;
	int i;
	int z, z2;
	double tf;


	UNUSED(verbose);

	for (i = 0; i < g_oaSingleMolecules.GetSize(); i++) {
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
		if (g_bLoadDipoleRestart) {
			int j;
			for (j = 0; j < 3; j++) {
				double val;
				if (fread(&val, sizeof(double), 1, g_fDipoleRestartFile) < 1) {
					eprintf("Could not read from dipole restart file. Setting value to 0.0.\n");
					val = 0.0;
				}
				sm->m_vDipole[j] = val;
			}
		} else {
			sm->m_fCharge = 0;
			sm->m_vDipole = CxDVector3(0.0, 0.0, 0.0);
			if (g_bVoroIntegrateQuadrupoleMoment)
				for (z=0;z<9;z++)
					sm->m_mQuadrupole[z] = 0.0;
		}
	}
	
// 	if (g_bTegri)
// 	{
// 		if (g_pTetraPak->m_bVoronoiCharges)
// 			g_pTetraPak->ProcessStep(this,verbose);
// 
// 		if (verbose)
// 			mprintf("\n");
// 	}
	
	for(i = 0; i < g_oaMolecules.GetSize(); i++) {
		CMolecule *m = (CMolecule *)g_oaMolecules[i];
		int j;
		for(j = 0; j < m->m_laSingleMolIndex.GetSize(); j++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[j]];
			CxDVector3 ref(0.0, 0.0, 0.0);
			if (!g_bDipoleRefFixed)
				ref = m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];
			if(m->m_iDipoleMode == 0) {
				sm->m_vDipole = CxDVector3(0.0, 0.0, 0.0);
			} else if(m->m_iDipoleMode == 1) {
				int k;
				for(k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
					if(m->m_baAtomIndex[k] == g_iVirtAtomType)
						continue;
					if(m->m_baAtomIndex[k] == g_iWannierAtomType)
						continue;
					CAtom *a = (CAtom *)g_oaAtoms[m->m_baAtomIndex[k]];
					int l;
					for(l = 0; l < m->m_waAtomCount[k]; l++) {
						sm->m_vDipole += a->m_fCharge * (m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref);
					}
				}
				for(k = 0; k < sm->m_laWannier.GetSize(); k++) {
					sm->m_vDipole -= g_fWannierCharge * (m_vaCoords[sm->m_laWannier[k]] - ref);
				}
				sm->m_vDipole *= DIP_EPM2DEBYE;
			} else if(m->m_iDipoleMode == 2) {
				if(j == 0) {
					char buf[BUF_SIZE];
					double dipole[3];
					if(fgets(buf, BUF_SIZE, m->m_pDipoleFile) == NULL || sscanf(buf, "%lf %lf %lf", &dipole[0], &dipole[1], &dipole[2]) < 3)
						eprintf("Could not read from dipole file. Setting dipole moment to ( 0.0 | 0.0 | 0.0 ).\n");
					else
						sm->m_vDipole = CxDVector3(dipole[0], dipole[1], dipole[2]);
				} else {
					sm->m_vDipole = ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole;
				}
			} else if(m->m_iDipoleMode == 3) {
				if(j == 0) {
					char comment[BUF_SIZE];
					strncpy(comment, m_pComment, BUF_SIZE-1 );
					comment[BUF_SIZE-1] = '\0';
					const char delim[] = "\t ";
					char *tok = strtok(comment, delim);
					int num = 0;
					bool ok[3] = { false, false, false };
					while(tok != NULL) {
						double f;
						if(sscanf(tok, "%lf", &f) == 1) {
							if(num == m->m_iDipoleCommentIndex[0]) {
								sm->m_vDipole[0] = f;
								ok[0] = true;
							}
							if(num == m->m_iDipoleCommentIndex[1]) {
								sm->m_vDipole[1] = f;
								ok[1] = true;
							}
							if(num == m->m_iDipoleCommentIndex[2]) {
								sm->m_vDipole[2] = f;
								ok[2] = true;
							}
							num++;
						}
						tok = strtok(NULL, delim);
					}
					if(!(ok[0] && ok[1] && ok[2])) {
						eprintf("Could not read from comment line. Setting dipole moment to ( 0.0 | 0.0 | 0.0 ).\n");
						sm->m_vDipole = CxDVector3(0.0, 0.0, 0.0);
					}
				} else {
					sm->m_vDipole = ((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole;
				}
			} else if(m->m_iDipoleMode == 4) {
				int k;
				for(k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
					if(m->m_baAtomIndex[k] == g_iVirtAtomType)
						continue;
					int l;
					for(l = 0; l < m->m_waAtomCount[k]; l++) {
						sm->m_vDipole += ((CxDoubleArray *)m->m_oaCharges[k])->GetAt(l) * (m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref);
					}
				}
				sm->m_vDipole *= DIP_EPM2DEBYE;
			} else if(m->m_iDipoleMode == 5) {
				int k;
				for(k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
					if(m->m_baAtomIndex[k] == g_iVirtAtomType)
						continue;
					int l;
					for(l = 0; l < m->m_waAtomCount[k]; l++) {
						sm->m_vDipole += m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] * (m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref);
					}
				}
				sm->m_vDipole *= DIP_EPM2DEBYE;
			} else if(m->m_iDipoleMode == 6) {
				int k;
				for(k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
					if(m->m_baAtomIndex[k] == g_iVirtAtomType)
						continue;
					int l;
					for(l = 0; l < m->m_waAtomCount[k]; l++) {
						sm->m_vDipole += m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] * (m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref);
					}
				}
				sm->m_vDipole *= DIP_EPM2DEBYE;
			} else if(m->m_iDipoleMode == 7) {
				int k;
				for(k = 0; k < m->m_baAtomIndex.GetSize(); k++) {
					if(m->m_baAtomIndex[k] == g_iVirtAtomType)
						continue;
					int l;
					for(l = 0; l < m->m_waAtomCount[k]; l++) {
						int ti = ((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l);
						sm->m_fCharge += m_faCharge[ti];
						sm->m_vDipole += m_dipoleMoments[ti] + m_faCharge[ti] * (m_vaCoords[ti] - ref);

						//mprintf("@ [%d|%d] %f + %f * (%f - %f)\n",k,l,m_dipoleMoments[ti][0],m_faCharge[ti],m_vaCoords[ti][0],ref[0]);

						if (g_bVoroIntegrateQuadrupoleMoment) {
							for (z=0;z<3;z++)
								for (z2=0;z2<3;z2++)

									sm->m_mQuadrupole(z,z2) +=
										  (m_maQuadTensor[ti](z,z2) * 1.0E6) // e*nm^2 to e*pm^2
										+ (m_dipoleMoments[ti][z]  * (m_vaCoords[ti][z2] - ref[z2]) )
										+ (m_dipoleMoments[ti][z2] * (m_vaCoords[ti][z]  - ref[z] ) )
										+ (m_faCharge[ti] * (m_vaCoords[ti][z] - ref[z]) * (m_vaCoords[ti][z2] - ref[z2]) );

					/*				sm->m_mQuadrupole(z,z2) +=
										  (m_maQuadTensor[ti](z,z2) * 1.0E6) // e*nm^2 to e*pm^2
										- (m_dipoleMoments[ti][z]  * (m_vaCoords[ti][z2] - ref[z2]) )
										- (m_dipoleMoments[ti][z2] * (m_vaCoords[ti][z]  - ref[z] ) )
										+ (m_faCharge[ti] * (m_vaCoords[ti][z] - ref[z]) * (m_vaCoords[ti][z2] - ref[z2]) ); */
						}

// 						mprintf(GREEN, "Dip: %12.6f %12.6f %12.6f\n", m_dipoleMoments[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][0], m_dipoleMoments[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][1], m_dipoleMoments[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][2]);
// 						mprintf(GREEN, "Crd: %12.6f %12.6f %12.6f\n", m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][0], m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][1], m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][2]);
// 						mprintf(GREEN, "Chg: %12.6f %12.6f %12.6f %12.6f\n", m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)], m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] * (m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref)[0], m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] * (m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref)[1], m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] * (m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)] - ref)[2]);
// 						mprintf(GREEN, "%.20g %.20g %.20g %.20g\n", m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)], ref[0], ref[1], ref[2]);
// 						mprintf(GREEN, "%g %g %g\n", m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][0], m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][1], m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[k])->GetAt(l)][2]);
					}
				}

				sm->m_vDipole *= DIP_EPM2DEBYE;

				if (g_bVoroIntegrateQuadrupoleMoment) {

					sm->m_mQuadrupole *= DIP_EPM2DEBYE; // e*pm^2 to Debye*pm

					if (!g_bQuadrupoleKeepTrace) {
						tf = 0.5 * (sm->m_mQuadrupole(0,0) + sm->m_mQuadrupole(1,1) + sm->m_mQuadrupole(2,2) );
						sm->m_mQuadrupole *= 1.5;
						sm->m_mQuadrupole(0,0) -= tf;
						sm->m_mQuadrupole(1,1) -= tf;
						sm->m_mQuadrupole(2,2) -= tf;
					}
				}

// 				mprintf(GREEN, "%g %g %g\n", sm->m_vDipole[0], sm->m_vDipole[1], sm->m_vDipole[2]);
			} else if(m->m_iDipoleMode == 8) {

				eprintf("CTimeStep::CalcDipoles(): Dipole mode 8 not implemented.\n");
				abort();

			} else {
				eprintf("CTimeStep::CalcDipoles(): Internal error.\n");
				abort();
			}

			// Last virtual atom of each molecule is defined as tip of molecular diople vector starting in #2
			int k;
			for(k = 0; k < m->m_baAtomIndex.GetSize(); k++)
				if(m->m_baAtomIndex[k] == g_iVirtAtomType)
					break;
			m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[k])->GetAt(m->m_laVirtualAtoms.GetSize()-1)] = 100.0*sm->m_vDipole + m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[k])->GetAt(1)];
		}
	}
	
//----------------- OLD -----------------------
// 	BTIN;
// 	int z, z2, z3, z4;
// 	CMolecule *m;
// 	CSingleMolecule *sm;
// 	CAtom *a;
// 	CxDVector3 dc;
// 
// //	mprintf("\n*** CalcDipoles ***");
// 
// 	for (z=0;z<g_oaMolecules.GetSize();z++)
// 	{
// 		m = (CMolecule*)g_oaMolecules[z];
// 		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
// 		{
// 			if (g_bVerbose)
// 				mprintf("\nCalcDipoles %s (%d):\n",m->m_sName,z2+1);
// 
// 			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
// 
// 			sm->m_vDipole = CxDVector3(0,0,0);
// 
// 			if ((!g_bWannier) && (!m->m_bChargesAssigned))
// 				continue;
// 
// 			if (g_bDipoleRefFixed)
// 				dc = CxDVector3(0,0,0);
// 					else dc = m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];
// 
// 			if (g_bVerbose)
// 				mprintf("  Ref. point is ( %f | %f | %f )\n",dc[0],dc[1],dc[2]);
// 
// 			for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
// 			{
// 				if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
// 					continue;
// 
// 				if (g_bWannier)
// 					if (m->m_baAtomIndex[z3] == g_iWannierAtomType)
// 						continue;
// 
// 				a = (CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]];
// 				for (z4=0;z4<m->m_waAtomCount[z3];z4++)
// 				{
// //					mprintf("  %s%d %.2f x ",a->m_sName,z4+1,a->m_fCharge);
// //					(m_vaCoords[sm->m_iAtomOffset[z3][z4]]-m_vaCoords[sm->m_iAtomOffset[m->Elements-1][0]]).Dump();
// //					mprintf("\n");
// //					mprintf("  %s%d %.2f (%d)\n",a->m_sName,z4+1,a->m_fCharge,((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4));
// 
// 					if (g_bWannier)
// 					{
// 						sm->m_vDipole += a->m_fCharge * (m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)] - dc);
// 
// 						if (g_bVerbose)
// 							mprintf("  %.2f: ( %f | %f | %f )\n",a->m_fCharge,m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][0],m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][1],m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][2]);
// 					} else
// 					{
// 						if (g_bReadChargesFrom4thXYZ)
// 						{
// 				//			mprintf("Moep. z3=%d, z4=%d, c=%f.\n",z3,z4,m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)]);
// 							sm->m_vDipole += m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)] * (m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)] - dc);
// 							if (g_bVerbose)
// 								mprintf("  %.2f: ( %f | %f | %f )\n",m_faCharge[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)],m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][0],m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][1],m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][2]);
// 						} else
// 						{
// 							sm->m_vDipole += ((CxFloatArray*)m->m_oaCharges[z3])->GetAt(z4) * (m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)] - dc);
// 							if (g_bVerbose)
// 								mprintf("  %.2f: ( %f | %f | %f )\n",((CxFloatArray*)m->m_oaCharges[z3])->GetAt(z4),m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][0],m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][1],m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)][2]);
// 						}
// 					}
// 				}
// 			}
// 			for (z3=0;z3<sm->m_laWannier.GetSize();z3++)
// 			{
// //				mprintf("  - %.2f x ",g_fWannierCharge);
// //				(m_vaCoords[sm->m_waWannier[z3]] - m_vaCoords[sm->m_iAtomOffset[m->Elements-1][0]]).Dump();
// //				mprintf("\n");
// 				sm->m_vDipole -= g_fWannierCharge * (m_vaCoords[sm->m_laWannier[z3]] - dc);
// 
// 				if (g_bVerbose)
// 					mprintf("  %.2f: ( %f | %f | %f )\n",-g_fWannierCharge,m_vaCoords[sm->m_laWannier[z3]][0],m_vaCoords[sm->m_laWannier[z3]][1],m_vaCoords[sm->m_laWannier[z3]][2]);
// 			}
// //			mprintf("  = ");
// //			sm->m_vDipole.Dump();
// //			mprintf("\n");
// 			sm->m_vDipole *= 0.048008f; // Conversion e*pm --> Debye
// 
// 			if (g_bVerbose)
// 				mprintf("  Result: %f.\n",sm->m_vDipole.GetLength());
// 
// //			mprintf("Molecule %s - Dipole %.3f\n",m->Name,sm->m_vDipole.GetLength());
// 		}
// 	}
// 	BTOUT;
}


void CTimeStep::CalcPolarizabilities() {

	int i;
	int j;

	if (g_iPolarizabilityMode == 1) {
		eprintf("CTimeStep::CalcPolarizabilities(): This is not implemented yet.\n");
		abort();
	} else if (g_iPolarizabilityMode == 2) {
		char fieldChar[4] = "xyz";
		for (i = 0; i < 3; i++) {
			if (g_iPolarizabilityConf[i] > 0) {
				CTimeStep tempTs1;
				CxDVec3Array tempDipole1;
				tempDipole1.SetMaxSize(g_oaSingleMolecules.GetSize());
				if (!tempTs1.ReadXYZ(g_fPolarizabilityFile[2 * i], false, &tempTs1.m_vaCoords)) {
					eprintf("CTimeStep::CalcPolarizabilities(): Error while reading Wannier centers with field along positive %c axis.\n", fieldChar[i]);
					abort();
				}
				if (g_bDoubleBox)
					tempTs1.DoubleBox();
				if (!g_bSaveCoordsUnchanged) {
					tempTs1.UniteMolecules(false);
					if (g_bRemoveCOM)
						tempTs1.CenterCOM();
				}
				tempTs1.CalcCenters();
				tempTs1.ScanWannier(false);
				for (j = 0; j < g_oaMolecules.GetSize(); j++) {
					CMolecule *m = (CMolecule *)g_oaMolecules[j];
					int k;
					for (k = 0; k < m->m_laSingleMolIndex.GetSize(); k++) {
						CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[k]];
						CxDVector3 ref(0.0, 0.0, 0.0);
						if (!g_bDipoleRefFixed) {
							ref = tempTs1.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];
						} else {
							CxDVector3 v = tempTs1.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_baAtomIndex.GetSize()-1])->GetAt(1)];
							int shift[3] = { 0, 0, 0 };
							if (g_bPeriodicX) {
								while (v[0] < -g_fBoxX / 2.0) { v[0] += g_fBoxX; shift[0]--; }
								while (v[0] > g_fBoxX / 2.0) { v[0] -= g_fBoxX; shift[0]++; }
							}
							if (g_bPeriodicY) {
								while (v[1] < -g_fBoxY / 2.0) { v[1] += g_fBoxY; shift[1]--; }
								while (v[1] > g_fBoxY / 2.0) { v[1] -= g_fBoxY; shift[1]++; }
							}
							if (g_bPeriodicZ) {
								while (v[2] < -g_fBoxZ / 2.0) { v[2] += g_fBoxZ; shift[2]--; }
								while (v[2] > g_fBoxZ / 2.0) { v[2] -= g_fBoxZ; shift[2]++; }
							}
							ref[0] = g_fBoxX * shift[0];
							ref[1] = g_fBoxY * shift[1];
							ref[2] = g_fBoxZ * shift[2];
						}
						CxDVector3 dip(0.0, 0.0, 0.0);
						int l;
						for (l = 0; l < m->m_baAtomIndex.GetSize(); l++) {
							if (m->m_baAtomIndex[l] == g_iVirtAtomType)
								continue;
							if (m->m_baAtomIndex[l] == g_iWannierAtomType)
								continue;
							CAtom *a = (CAtom *)g_oaAtoms[m->m_baAtomIndex[l]];
							int n;
							for (n = 0; n < m->m_waAtomCount[l]; n++) {
								dip += a->m_fCharge * (tempTs1.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[l])->GetAt(n)] - ref);
							}
						}
						for (l = 0; l < sm->m_laWannier.GetSize(); l++) {
							dip -= g_fWannierCharge * (tempTs1.m_vaCoords[sm->m_laWannier[l]] - ref);
						}
						dip *= DIP_EPM2DEBYE;
						tempDipole1.Add(dip);
					}
				}
				if (g_iPolarizabilityConf[i] > 1) {
					CTimeStep tempTs2;
					if (!tempTs2.ReadXYZ(g_fPolarizabilityFile[2 * i + 1], false, &tempTs2.m_vaCoords)) {
						eprintf("CTimeStep::CalcPolarizabilities(): Error while reading Wannier centers with field along negative %c axis.\n", fieldChar[i]);
						abort();
					}
					if (g_bDoubleBox)
						tempTs2.DoubleBox();
					if (!g_bSaveCoordsUnchanged) {
						tempTs2.UniteMolecules(false);
						if (g_bRemoveCOM)
							tempTs2.CenterCOM();
					}
					tempTs2.CalcCenters();
					tempTs2.ScanWannier(false);
					int z = 0;
					for (j = 0; j < g_oaMolecules.GetSize(); j++) {
						CMolecule *m = (CMolecule *)g_oaMolecules[j];
						int k;
						for (k = 0; k < m->m_laSingleMolIndex.GetSize(); k++) {
							CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[k]];
							CxDVector3 ref(0.0, 0.0, 0.0);
							if (!g_bDipoleRefFixed) {
								ref = tempTs2.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];
							} else {
								CxDVector3 v = tempTs2.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_baAtomIndex.GetSize()-1])->GetAt(1)];
								int shift[3] = { 0, 0, 0 };
								if (g_bPeriodicX) {
									while (v[0] < -g_fBoxX / 2.0) { v[0] += g_fBoxX; shift[0]--; }
									while (v[0] > g_fBoxX / 2.0) { v[0] -= g_fBoxX; shift[0]++; }
								}
								if (g_bPeriodicY) {
									while (v[1] < -g_fBoxY / 2.0) { v[1] += g_fBoxY; shift[1]--; }
									while (v[1] > g_fBoxY / 2.0) { v[1] -= g_fBoxY; shift[1]++; }
								}
								if (g_bPeriodicZ) {
									while (v[2] < -g_fBoxZ / 2.0) { v[2] += g_fBoxZ; shift[2]--; }
									while (v[2] > g_fBoxZ / 2.0) { v[2] -= g_fBoxZ; shift[2]++; }
								}
								ref[0] = g_fBoxX * shift[0];
								ref[1] = g_fBoxY * shift[1];
								ref[2] = g_fBoxZ * shift[2];
							}
							CxDVector3 dip(0.0, 0.0, 0.0);
							int l;
							for (l = 0; l < m->m_baAtomIndex.GetSize(); l++) {
								if (m->m_baAtomIndex[l] == g_iVirtAtomType)
									continue;
								if (m->m_baAtomIndex[l] == g_iWannierAtomType)
									continue;
								CAtom *a = (CAtom *)g_oaAtoms[m->m_baAtomIndex[l]];
								int n;
								for (n = 0; n < m->m_waAtomCount[l]; n++) {
									dip += a->m_fCharge * (tempTs2.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[l])->GetAt(n)] - ref);
								}
							}
							for (l = 0; l < sm->m_laWannier.GetSize(); l++) {
								dip -= g_fWannierCharge * (tempTs2.m_vaCoords[sm->m_laWannier[l]] - ref);
							}
							dip *= DIP_EPM2DEBYE;
							dip = tempDipole1[z] - dip;
							for (l = 0; l < 3; l++) {
								sm->m_polarizability[3 * l + i] = -1.0 * 0.5 * dip[l] / g_fPolarizabilityFieldStrength * 0.393430; // Multiply by -1.0f to work with CP2K field definition
							}
							z++;
						}
					}
				} else {
					int z = 0;
					for (j = 0; j < g_oaMolecules.GetSize(); j++) {
						CMolecule *m = (CMolecule *)g_oaMolecules[j];
						int k;
						for (k = 0; k < m->m_laSingleMolIndex.GetSize(); k++) {
							CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[k]];
							int l;
							CxDVector3 dip = tempDipole1[z] - sm->m_vDipole;
							for (l = 0; l < 3; l++) {
								sm->m_polarizability[3 * l + i] = -1.0 * dip[l] / g_fPolarizabilityFieldStrength * 0.393430;
							}
							z++;
						}
					}
				}
			} else {
				for (j = 0; j < g_oaMolecules.GetSize(); j++) {
					CMolecule *m = (CMolecule *)g_oaMolecules[j];
					int k;
					for (k = 0; k < m->m_laSingleMolIndex.GetSize(); k++) {
						CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[k]];
						int l;
						for (l = 0; l < 3; l++) {
							sm->m_polarizability[3 * l + i] = 0.0;
						}
					}
				}
			}
		}
	} else if (g_iPolarizabilityMode == 3) {
		char fieldChar[4] = "xyz";
		for (i = 0; i < 3; i++) {
			if (g_iPolarizabilityConf[i] > 0) {
				CTimeStep tempTs1;
				CxDVec3Array tempDipole1;
				tempDipole1.SetMaxSize(g_oaSingleMolecules.GetSize());
				if (!tempTs1.ReadCube(g_fPolarizabilityFile[2 * i], false)) {
					eprintf("CTimeStep::CalcPolarizabilities(): Error while reading electron density with field along positive %c axis.\n", fieldChar[i]);
					abort();
				}
				if (g_bDoubleBox)
					tempTs1.DoubleBox();
				if (!g_bSaveCoordsUnchanged) {
					tempTs1.UniteMolecules(false);
					if (g_bRemoveCOM)
						tempTs1.CenterCOM();
				}
				tempTs1.CalcCenters();
				g_pTetraPak->ProcessStep(&tempTs1, false);
				for (j = 0; j < g_oaMolecules.GetSize(); j++) {
					CMolecule *m = (CMolecule *)g_oaMolecules[j];
					int k;
					for (k = 0; k < m->m_laSingleMolIndex.GetSize(); k++) {
						CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[k]];
						CxDVector3 ref(0.0, 0.0, 0.0);
						if (!g_bDipoleRefFixed) {
							ref = tempTs1.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];
						} else {
							CxDVector3 v = tempTs1.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_baAtomIndex.GetSize()-1])->GetAt(1)];
							int shift[3] = { 0, 0, 0 };
							if (g_bPeriodicX) {
								while (v[0] < -g_fBoxX / 2.0) { v[0] += g_fBoxX; shift[0]--; }
								while (v[0] > g_fBoxX / 2.0) { v[0] -= g_fBoxX; shift[0]++; }
							}
							if (g_bPeriodicY) {
								while (v[1] < -g_fBoxY / 2.0) { v[1] += g_fBoxY; shift[1]--; }
								while (v[1] > g_fBoxY / 2.0) { v[1] -= g_fBoxY; shift[1]++; }
							}
							if (g_bPeriodicZ) {
								while (v[2] < -g_fBoxZ / 2.0) { v[2] += g_fBoxZ; shift[2]--; }
								while (v[2] > g_fBoxZ / 2.0) { v[2] -= g_fBoxZ; shift[2]++; }
							}
							ref[0] = g_fBoxX * shift[0];
							ref[1] = g_fBoxY * shift[1];
							ref[2] = g_fBoxZ * shift[2];
						}
						CxDVector3 dip(0.0, 0.0, 0.0);
						int l;
						for (l = 0; l < m->m_baAtomIndex.GetSize(); l++) {
							if (m->m_baAtomIndex[l] == g_iVirtAtomType)
								continue;
							int n;
							for (n = 0; n < m->m_waAtomCount[l]; n++) {
								dip += tempTs1.m_dipoleMoments[((CxIntArray *)sm->m_oaAtomOffset[l])->GetAt(n)] + tempTs1.m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[l])->GetAt(n)] * (tempTs1.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[l])->GetAt(n)] - ref);
							}
						}
						dip *= DIP_EPM2DEBYE;
						tempDipole1.Add(dip);
					}
				}
				if (g_iPolarizabilityConf[i] > 1) {
					CTimeStep tempTs2;
					if (!tempTs2.ReadCube(g_fPolarizabilityFile[2 * i + 1], false)) {
						eprintf("CTimeStep::CalcPolarizabilities(): Error while reading electron density with field along negative %c axis.\n", fieldChar[i]);
						abort();
					}
					if (g_bDoubleBox)
						tempTs2.DoubleBox();
					if (!g_bSaveCoordsUnchanged) {
						tempTs2.UniteMolecules(false);
						if (g_bRemoveCOM)
							tempTs2.CenterCOM();
					}
					tempTs2.CalcCenters();
					g_pTetraPak->ProcessStep(&tempTs2, false);
					int z = 0;
					for (j = 0; j < g_oaMolecules.GetSize(); j++) {
						CMolecule *m = (CMolecule *)g_oaMolecules[j];
						int k;
						for (k = 0; k < m->m_laSingleMolIndex.GetSize(); k++) {
							CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[k]];
							CxDVector3 ref(0.0, 0.0, 0.0);
							if (!g_bDipoleRefFixed) {
								ref = tempTs2.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];
							} else {
								CxDVector3 v = tempTs2.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[m->m_baAtomIndex.GetSize()-1])->GetAt(1)];
								int shift[3] = { 0, 0, 0 };
								if (g_bPeriodicX) {
									while (v[0] < -g_fBoxX / 2.0) { v[0] += g_fBoxX; shift[0]--; }
									while (v[0] > g_fBoxX / 2.0) { v[0] -= g_fBoxX; shift[0]++; }
								}
								if (g_bPeriodicY) {
									while (v[1] < -g_fBoxY / 2.0) { v[1] += g_fBoxY; shift[1]--; }
									while (v[1] > g_fBoxY / 2.0) { v[1] -= g_fBoxY; shift[1]++; }
								}
								if (g_bPeriodicZ) {
									while (v[2] < -g_fBoxZ / 2.0) { v[2] += g_fBoxZ; shift[2]--; }
									while (v[2] > g_fBoxZ / 2.0) { v[2] -= g_fBoxZ; shift[2]++; }
								}
								ref[0] = g_fBoxX * shift[0];
								ref[1] = g_fBoxY * shift[1];
								ref[2] = g_fBoxZ * shift[2];
							}
							CxDVector3 dip(0.0, 0.0, 0.0);
							int l;
							for (l = 0; l < m->m_baAtomIndex.GetSize(); l++) {
								if (m->m_baAtomIndex[l] == g_iVirtAtomType)
									continue;
								int n;
								for (n = 0; n < m->m_waAtomCount[l]; n++) {
									dip += tempTs2.m_dipoleMoments[((CxIntArray *)sm->m_oaAtomOffset[l])->GetAt(n)] + tempTs2.m_faCharge[((CxIntArray *)sm->m_oaAtomOffset[l])->GetAt(n)] * (tempTs2.m_vaCoords[((CxIntArray *)sm->m_oaAtomOffset[l])->GetAt(n)] - ref);
								}
							}
							dip *= DIP_EPM2DEBYE;
							dip = tempDipole1[z] - dip;
							for (l = 0; l < 3; l++) {
								sm->m_polarizability[3 * l + i] = -1.0 * 0.5 * dip[l] / g_fPolarizabilityFieldStrength * 0.393430; // Multiply by -1.0f to work with CP2K field definition
							}
							z++;
						}
					}
				} else {
					int z = 0;
					for (j = 0; j < g_oaMolecules.GetSize(); j++) {
						CMolecule *m = (CMolecule *)g_oaMolecules[j];
						int k;
						for (k = 0; k < m->m_laSingleMolIndex.GetSize(); k++) {
							CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[k]];
							int l;
							CxDVector3 dip = tempDipole1[z] - sm->m_vDipole;
							for (l = 0; l < 3; l++) {
								sm->m_polarizability[3 * l + i] = -1.0 * dip[l] / g_fPolarizabilityFieldStrength * 0.393430;
							}
							z++;
						}
					}
				}
			} else {
				for (j = 0; j < g_oaMolecules.GetSize(); j++) {
					CMolecule *m = (CMolecule *)g_oaMolecules[j];
					int k;
					for (k = 0; k < m->m_laSingleMolIndex.GetSize(); k++) {
						CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[k]];
						int l;
						for (l = 0; l < 3; l++) {
							sm->m_polarizability[3 * l + i] = 0.0;
						}
					}
				}
			}
		}
	} else if (g_iPolarizabilityMode == 4) {
		char fieldChar[4] = "xyz";
		for (i = 0; i < 3; i++) {
			if (g_iPolarizabilityConf[i] > 0) {
				CxDVec3Array tempDipole1;
				tempDipole1.SetMaxSize(g_oaSingleMolecules.GetSize());
				for (j = 0; j < g_oaSingleMolecules.GetSize(); j++) {
					CxDVector3 dip;
					int k;
					for (k = 0; k < 3; k++) {
						double val;
						if (fread(&val, sizeof(double), 1, g_fPolarizabilityFile[2 * i]) < 1) {
							eprintf("CTimeStep::CalcPolarizabilities(): Error while reading dipole restart file with field along positive %c axis.\n", fieldChar[i]);
							abort();
						}
						dip[k] = val;
					}
					tempDipole1.Add(dip);
				}
				if (g_iPolarizabilityConf[i] > 1) {
					for (j = 0; j < g_oaSingleMolecules.GetSize(); j++) {
						CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[j];
						CxDVector3 dip;
						int k;
						for (k = 0; k < 3; k++) {
							double val;
							if (fread(&val, sizeof(double), 1, g_fPolarizabilityFile[2 * i + 1]) < 1) {
								eprintf("CTimeStep::CalcPolarizabilities(): Error while reading dipole restart file with field along negative %c axis.\n", fieldChar[i]);
								abort();
							}
							dip[k] = val;
						}
						dip = tempDipole1[j] - dip;
						for (k = 0; k < 3; k++) {
							sm->m_polarizability[3 * k + i] = -1.0 * 0.5 * dip[k] / g_fPolarizabilityFieldStrength * 0.393430; // Multiply by -1.0f to work with CP2K field definition
						}
					}
				} else {
					for (j = 0; j < g_oaSingleMolecules.GetSize(); j++) {
						CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[j];
						CxDVector3 dip;
						dip = tempDipole1[j] - sm->m_vDipole;
						int k;
						for (k = 0; k < 3; k++) {
							sm->m_polarizability[3 * k + i] = -1.0 * dip[k] / g_fPolarizabilityFieldStrength * 0.393430;
						}
					}
				}
			} else {
				for (j = 0; j < g_oaMolecules.GetSize(); j++) {
					CMolecule *m = (CMolecule *)g_oaMolecules[j];
					int k;
					for (k = 0; k < m->m_laSingleMolIndex.GetSize(); k++) {
						CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[k]];
						int l;
						for (l = 0; l < 3; l++) {
							sm->m_polarizability[3 * l + i] = 0.0;
						}
					}
				}
			}
		}
	} else {
		eprintf("CTimeStep::CalcPolarizabilities(): Internal error.\n");
		abort();
	}
}

void CTimeStep::DoubleBox()
{
	int px, py, pz, z;
	char *p;
	CxDVector3 vec;

	if (m_vaCoords.GetSize() < g_iGesAtomCount)
		m_vaCoords.SetSize(g_iGesAtomCount);

	if (m_paLabels.GetSize() != 0)
		m_paLabels.SetSize(g_iGesAtomCount);

/*	if (g_bBoxNonOrtho)
	{
		for (z=0;z<g_iGesAtomCount/g_iDoubleBoxFactor;z++)
		{
			vec = g_mBoxToOrtho * m_vaCoords[z];

			while (vec[0] >= 1.0/g_iDoubleBoxX)
				vec[0] -= 1.0/g_iDoubleBoxX;
			while (vec[0] < 0)
				vec[0] += 1.0/g_iDoubleBoxX;
			while (vec[1] >= 1.0/g_iDoubleBoxY)
				vec[1] -= 1.0/g_iDoubleBoxY;
			while (vec[1] < 0)
				vec[1] += 1.0/g_iDoubleBoxY;
			while (vec[2] >= 1.0/g_iDoubleBoxZ)
				vec[2] -= 1.0/g_iDoubleBoxZ;
			while (vec[2] < 0)
				vec[2] += 1.0/g_iDoubleBoxZ;

			m_vaCoords[z] = g_mBoxFromOrtho * vec;
		}
	}*/

	for (pz=0;pz<g_iDoubleBoxZ;pz++)
	{
		for (py=0;py<g_iDoubleBoxY;py++)
		{
			for (px=0;px<g_iDoubleBoxX;px++)
			{
				if ((px == 0) && (py == 0) && (pz == 0))
					continue;

				for (z=0;z<g_iGesAtomCount/g_iDoubleBoxFactor;z++)
				{
					if (g_bBoxNonOrtho)
					{
//						mprintf("### px=%d py=%d pz=%d z=%d\n",px,py,pz,z);
//						mprintf("    vec ( %f | %f | %f )\n",m_vaCoords[z][0],m_vaCoords[z][1],m_vaCoords[z][2]);
						vec = g_mBoxToOrtho * m_vaCoords[z];
//						mprintf("    wird zu ( %f | %f | %f )\n",vec[0],vec[1],vec[2]);

/*						while (vec[0] >= 1.0/g_iDoubleBoxX)
							vec[0] -= 1.0/g_iDoubleBoxX;
						while (vec[0] < 0)
							vec[0] += 1.0/g_iDoubleBoxX;
						while (vec[1] >= 1.0/g_iDoubleBoxY)
							vec[1] -= 1.0/g_iDoubleBoxY;
						while (vec[1] < 0)
							vec[1] += 1.0/g_iDoubleBoxY;
						while (vec[2] >= 1.0/g_iDoubleBoxZ)
							vec[2] -= 1.0/g_iDoubleBoxZ;
						while (vec[2] < 0)
							vec[2] += 1.0/g_iDoubleBoxZ;*/

//						mprintf("    gefaltet ( %f | %f | %f )\n",vec[0],vec[1],vec[2]);

						vec += CxDVector3((double)px/g_iDoubleBoxX,(double)py/g_iDoubleBoxY,(double)pz/g_iDoubleBoxZ);

//						mprintf("    geshiftet ( %f | %f | %f )\n",vec[0],vec[1],vec[2]);

						m_vaCoords[(pz*g_iDoubleBoxX*g_iDoubleBoxY+py*g_iDoubleBoxX+px)*g_iGesAtomCount/g_iDoubleBoxFactor+z] = g_mBoxFromOrtho * vec;

//						mprintf("    Ergebnis ( %f | %f | %f )\n",m_vaCoords[(pz*g_iDoubleBoxX*g_iDoubleBoxY+py*g_iDoubleBoxX+px)*g_iGesAtomCount/g_iDoubleBoxFactor+z][0],m_vaCoords[(pz*g_iDoubleBoxX*g_iDoubleBoxY+py*g_iDoubleBoxX+px)*g_iGesAtomCount/g_iDoubleBoxFactor+z][1],m_vaCoords[(pz*g_iDoubleBoxX*g_iDoubleBoxY+py*g_iDoubleBoxX+px)*g_iGesAtomCount/g_iDoubleBoxFactor+z][2]);
					} else
					{
						m_vaCoords[(pz*g_iDoubleBoxX*g_iDoubleBoxY+py*g_iDoubleBoxX+px)*g_iGesAtomCount/g_iDoubleBoxFactor+z] = m_vaCoords[z] + CxDVector3(px*g_fBoxX/g_iDoubleBoxX,py*g_fBoxY/g_iDoubleBoxY,pz*g_fBoxZ/g_iDoubleBoxZ);
					}

					if (m_paLabels.GetSize() != 0)
					{
						try { p = new char[strlen((char*)m_paLabels[z])+1]; } catch(...) { p = NULL; }
						if (p == NULL) NewException((double)(strlen((char*)m_paLabels[z])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						strcpy(p,(char*)m_paLabels[z]);
						m_paLabels[(pz*g_iDoubleBoxX*g_iDoubleBoxY+py*g_iDoubleBoxX+px)*g_iGesAtomCount/g_iDoubleBoxFactor+z] = p;
					}
				}
			}
		}
	}
	m_iGesAtomCount = g_iGesAtomCount;
}


void CTimeStep::DoubleBoxVelocity()
{
	int px, py, pz, z;

	if (m_vaVelocities.GetSize() < g_iGesAtomCount)
		m_vaVelocities.SetSize(g_iGesAtomCount);

	for (pz=0;pz<g_iDoubleBoxZ;pz++)
	{
		for (py=0;py<g_iDoubleBoxY;py++)
		{
			for (px=0;px<g_iDoubleBoxX;px++)
			{
				if ((px == 0) && (py == 0) && (pz == 0))
					continue;
				for (z=0;z<g_iGesAtomCount/g_iDoubleBoxFactor;z++)
					m_vaVelocities[(pz*g_iDoubleBoxX*g_iDoubleBoxY+py*g_iDoubleBoxX+px)*g_iGesAtomCount/g_iDoubleBoxFactor+z] = m_vaVelocities[z];
			}
		}
	}
}


void CTimeStep::DoubleBoxForce()
{
	int px, py, pz, z;

	if (m_vaForces.GetSize() < g_iGesAtomCount)
		m_vaForces.SetSize(g_iGesAtomCount);

	for (pz=0;pz<g_iDoubleBoxZ;pz++)
	{
		for (py=0;py<g_iDoubleBoxY;py++)
		{
			for (px=0;px<g_iDoubleBoxX;px++)
			{
				if ((px == 0) && (py == 0) && (pz == 0))
					continue;
				for (z=0;z<g_iGesAtomCount/g_iDoubleBoxFactor;z++)
					m_vaForces[(pz*g_iDoubleBoxX*g_iDoubleBoxY+py*g_iDoubleBoxX+px)*g_iGesAtomCount/g_iDoubleBoxFactor+z] = m_vaForces[z];
			}
		}
	}
}


bool CTimeStep::ReadXYZ(FILE *a, bool needinfo, CxDVec3Array *v)
{
	BTIN;
	char buf[256], obuf[256],  *p, *q, *r;
	int z, /*i,*/ j;
	const char *separators = " \t";

//	mprintf("*** Read Anfang.\n");
	j = 0;
//	m_iSizeBytes = 0;
_readagain:
	buf[0] = 0;
	(void)fgets_bin(buf,256,a);
	if (feof(a)) return false;
//	m_iSizeBytes += strlen(buf);
	if (strlen(buf) > 0)
		buf[strlen(buf)-1] = 0;
//	mprintf("Read1: \"%s\".\n",buf);
	goto _firsttry;
_again:
	eprintf("\nTrajectory file seems to be damaged. Searching for next time step...");
_firsttry:
	if (feof(a))
	{
//		mprintf("CTimeStep::ReadTimestep(): Unexpected End of File (1).\n"); 
		BTOUT; 
		return false;
	}
	if (strchr(buf,'.') != NULL)
	{
		if (j < 10)
			eprintf("x");
//			mprintf("x: %d - \"%s\"\n",j,buf);
		j++;
		goto _readagain;
	}
	if (!IsValidInteger(buf)) {
		eprintf("\nCTimeStep::ReadXYZ(): Error: Atom count field is not a valid integer (\"%s\"). Your trajectory file is damaged.",buf);
		m_iGesAtomCount = 0;
		goto _readagain;
	}
	m_iGesAtomCount = atoi(buf);
	if (atoi(buf) < 0)
	{
		eprintf("\nCTimeStep::ReadXYZ(): Error: Atom count = %d < 0. \"%s\"",atoi(buf),buf);
		m_iGesAtomCount = 0;
		return false;
	}
	if (m_iGesAtomCount == 0)
	{
		eprintf("\nCTimeStep::ReadXYZ(): Error: Atom count = 0 (\"%s\"). Your trajectory file is damaged.",buf);
		goto _readagain;
	}
	if (m_pComment == NULL)
	{
		try { m_pComment = new char[256]; } catch(...) { m_pComment = NULL; }
		if (m_pComment == NULL) NewException((double)256*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}
	m_pComment[0] = 0;
	(void)fgets_bin(m_pComment,255,a);
	if (feof(a))
		return false;
	if (strlen(m_pComment) >= 254)
		eprintf("\nCTimeStep::ReadXYZ(): Warning: Comment line longer than 254 characters.  ");
	if (needinfo)
	{
		for (z=0;z<m_paLabels.GetSize();z++)
			delete[] (char*)m_paLabels[z];
		m_paLabels.RemoveAll();
		if (g_bDoubleBox)
			m_paLabels.SetSize(m_iGesAtomCount*g_iDoubleBoxFactor);
		else
			m_paLabels.SetSize(m_iGesAtomCount);
	}
	if (v->GetSize() < (long)m_iGesAtomCount)
		v->SetSize(m_iGesAtomCount);

	if (g_bKeepOriginalCoords)
	{
		if (m_vaCoords_Original.GetSize() < (long)m_iGesAtomCount)
			m_vaCoords_Original.SetSize(m_iGesAtomCount);
	}

//	if (strlen(m_pComment) > 0)
//		m_pComment[strlen(m_pComment)-1] = 0;

	p = m_pComment;
	while (*p != 0)
	{
		if ((*p == 10) || (*p == 13))
			*p = 0;
		p++;
	}
//	mprintf("# Comment: \"%s\".\n",m_pComment);

	if (g_bNPT && g_bXYZComment3Numbers)
		ExtractXYZCellGeometry3(m_pComment);

	if (g_bNPT && g_bXYZComment6Numbers)
		ExtractXYZCellGeometry6(m_pComment);

	if (g_bNPT && g_bXYZComment9Numbers)
		ExtractXYZCellGeometry9(m_pComment);

	if (needinfo)
	{
		z = 0;
		p = m_pComment;
		while (*p != 0)
		{
//			mprintf("# (A) \"%s\".\n",p);

//			while (strchr(" \t,",*p) != NULL)
//				p++;

			while (strchr("0123456789+-.Ee",*p) == NULL)
				p++;

			q = p;
//			mprintf("# (B) \"%s\".\n",p);
//			while (strchr("0123456789+-.Ee",*q) != NULL)
			while ((strchr("0123456789+-.Ee",*q) != NULL) && (*q != 0))
				q++;
			if (*(q-1) == 0)
				q--;
//			mprintf("# (C) \"%s\".\n",q);
			if (q > p)
			{
				memcpy(buf,p,q-p);
				buf[q-p] = 0;
				if (atof(buf) != 0)
				{
	//				mprintf("\n## \"%s\" --> %f.",buf,atof(buf));
					z++;
				}
			}
			if (*q == 0)
				break;
			p = q;
	//		mprintf("\nNow p=%d (%c).",*p,*p);
		}
		if (z == 6)
			g_bXYZComment6Numbers = true;
		else if (z == 9)
			g_bXYZComment9Numbers = true;
		else if (z == 3) {
			if (strspn(m_pComment,"0123456789+-.Ee \r\n\t") == strlen(m_pComment))
				g_bXYZComment3Numbers = true;
		}
	}

//	mprintf("Comment: \"%s\".\n",m_pComment);
//	m_iSizeBytes += strlen(m_pComment);

	if (g_bReadChargesFrom4thXYZ)
		m_faCharge.SetSize(m_iGesAtomCount);

	for (z=0;z<(long)m_iGesAtomCount;z++) // den ersten Zeitschritt einlesen
	{
		buf[0] = 0;
		(void)fgets_bin(buf,256,a);
		if (feof(a))
		{
			eprintf("\nCTimeStep::ReadXYZ(): Unexpected end of file (2). \"%s\"\n",buf);
			BTOUT; 
			return false;
		}
//		m_iSizeBytes += strlen(buf);
		buf[strlen(buf)-1] = 0;
//		mprintf("  %d: \"%s\".\n",z,buf);
		strcpy(obuf,buf);
//		i = 0;

		p = &buf[0];
		while (strchr(separators,*p) != NULL)
			p++;
		q = p+1;
		while ((strchr(separators,*q) == NULL) && (*q != 0))
			q++;
		if (*q == 0)
		{
			eprintf("\nCTimeStep::ReadXYZ(): %d: Incomplete line (1): \"%s\"\n",z+1,obuf);
			BTOUT; 
			return false;
		}
		*q = 0;

		if (needinfo)
		{
	//		if (strlen(p) > 7)
	//		{
	//			eprintf("\nCTimeStep::ReadXYZ(): \"%s\" - Maximum length of atom labels is 7 chars; truncating.\n",p);
	//			p[7] = 0;
	//		}

			try { r = new char[strlen(p)+1]; } catch(...) { r = NULL; }
			if (r == NULL) NewException((double)(strlen(p)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(r,p);
			m_paLabels[z] = r;
		}

		q++;
		while (strchr(separators,*q) != NULL)
			q++;
		p = q;
		while ((strchr(separators,*p) == NULL) && (*p != 0))
			p++;
		if (*p == 0)
		{
			eprintf("\nCTimeStep::ReadXYZ(): %d: Incomplete line (2): \"%s\"",z+1,obuf);
			goto _again;
		}
		*p = 0;

		#ifdef FAST_ATOF
			(*v)[z][0] = fast_atof(q) * 100.0;
		#else
			(*v)[z][0] = atof(q) * 100.0;
		#endif

	//	if (g_bKeepOriginalCoords)
	//		m_vaCoords_Original[z][0] = (*v)[z][0];

		q = p+1;
		while (strchr(separators,*q) != NULL)
			q++;
		p = q;
		while ((strchr(separators,*p) == NULL) && (*p != 0))
			p++;
		if (*p == 0)
		{
			eprintf("\nCTimeStep::ReadXYZ(): %d: Incomplete line (3) \"%s\"",z+1,obuf);
			goto _again;
		}
		*p = 0;

		#ifdef FAST_ATOF
			(*v)[z][1] = fast_atof(q) * 100.0;
		#else
			(*v)[z][1] = atof(q) * 100.0;
		#endif

	//	if (g_bKeepOriginalCoords)
	//		m_vaCoords_Original[z][1] = (*v)[z][1];

		q = p+1;
		while (strchr(separators,*q) != NULL)
			q++;
		p = q;
		while ((strchr(separators,*p) == NULL) && (*p != 0))
			p++;
		if (g_bReadChargesFrom4thXYZ)
		{
			if (*p == 0)
			{
				eprintf("\nCTimeStep::ReadXYZ(): %d: Incomplete line (4) \"%s\"",z+1,obuf);
				goto _again;
			}
			*p = 0;
		} else
		{
			if (*p != 0)
			{
				*p = 0;
				if (needinfo && (z==0))
					p++;
			}
		}

		#ifdef FAST_ATOF
			(*v)[z][2] = fast_atof(q) * 100.0;
		#else
			(*v)[z][2] = atof(q) * 100.0;
		#endif

//		mprintf("\n Have %f  %f  %f.",(*v)[z][0],(*v)[z][1],(*v)[z][2]);

	//	if (g_bKeepOriginalCoords)
	//		m_vaCoords_Original[z][2] = (*v)[z][2];

		if (g_bKeepOriginalCoords)
			m_vaCoords_Original[z] = (*v)[z];

		if (g_bReadChargesFrom4thXYZ)
		{
			q = p+1;
			while (strchr(separators,*q) != NULL)
				q++;
			p = q;
			while ((strchr(separators,*p) == NULL) && (*p != 0))
				p++;
			*p = 0;

			#ifdef FAST_ATOF
				m_faCharge[z] = fast_atof(q);
			#else
				m_faCharge[z] = atof(q);
			#endif

		} else if (needinfo && (z==0))
		{
			if (*p == 0)
				goto _no4;
			q = p+1;
			if (*q == 0)
				goto _no4;
			while (strchr(separators,*q) != NULL)
				q++;
			p = q;
			while ((strchr(separators,*p) == NULL) && (*p != 0))
				p++;
			if ((p-q) > 2)
			{
				for (;q<p;q++)
					if (!(((*q >= '0') && (*q <= '9')) || (*q == '-') || (*q == '.')))
						goto _no4;
				g_bXYZ4thCol = true;
_no4:;
			}
		}

/*		q = strchr(buf,'.');
		if (q != NULL)
		{
			i++;
			q = strchr(q+1,'.');
		}
		if (q != NULL)
		{
			i++;
			q = strchr(q+1,'.');
		}
		if (q==NULL)
		{
			eprintf("\nCTimeStep::ReadXYZ(): Error - only %d/3 dots. %d, \"%s\"",i,z+1,buf);
			goto _again;
		}
		q = buf;
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p == NULL)
		{
			eprintf("\nCTimeStep::ReadXYZ(): Error 1. %d, \"%s\"",z+1,buf);
			goto _again;
		}
		while (isdigit(*(p-1)) && (p > buf))
			p--;
		if (p == buf)
		{
			eprintf("\nCTimeStep::ReadXYZ(): No Atom label found. %d, \"%s\"",z+1,buf);
			goto _again;
		}
		*p = 0;
		if (needinfo)
		{
			r = new char[strlen(q)+1];
			strcpy(r,q);
			m_paLabels[z] = r;
		}
		q = p+1;
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p == NULL)
		{
			eprintf("\nCTimeStep::ReadXYZ(): Error 3. %d, \"%s\"",z+1,q);
			goto _again;
		}
		*p = 0;
		(*v)[z][0] = (double)atof(q) * 100.0f;
		q = p+1;
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p == NULL)
		{
			eprintf("\nCTimeStep::ReadXYZ(): Error 4. %d, \"%s\"",z+1,q);
			goto _again;
		}
		*p = 0;
		(*v)[z][1] = (double)atof(q) * 100.0f;
		q = p+1;
		while (*q == ' ')
			q++;
		p = strchr(q,' ');
		if (p != NULL)
			*p = 0;
		(*v)[z][2] = (double)atof(q) * 100.0f;*/


	}
//	mprintf("*** Read Ende.\n");
	BTOUT; 
	return true;
}


bool CTimeStep::ReadCube(FILE *a, bool needinfo)
{
	int rx, ry, rz, ac, z, z2, k, i;
	char buf[256], *p, *r;
	const char *q;
	double fx, fy, fz, tf;

	(void)!fgets(buf,256,a);
	(void)!fgets(buf,256,a);

	if (feof(a))
	{
		eprintf("\nCTimeStep::ReadCube(): Error: End of file reached while reading cube file.\n");
		return false;
	}

	(void)!fgets(buf,256,a);
	if (feof(a)) { eprintf("CTimeStep::ReadCube(): Error: Unexpected end of file reached while reading cube file.\n"); return false; }
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;

	ac = atoi(buf);
	m_iGesAtomCount = ac;

	if (ac <= 0)
	{
		eprintf("CTimeStep::ReadCube(): Error: Encountered atom count of %d while reading cube file.\n",ac);
		return false;
	}

	(void)!fgets(buf,256,a);
	if (feof(a)) { eprintf("CTimeStep::ReadCube(): Error: Unexpected end of file reached while reading cube file.\n"); return false; }
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;
	rx = atoi(buf);
	g_iCubeRes[0] = rx;
	
	//static bool first = true;
	// replaced by m_bFirst
	
	for (k = 0; k < 3; k++) {
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		*p = 0;
		
		tf = atof(q);
		if (g_bFlipCellVectorX && (k == 0) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorY && (k == 1) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorZ && (k == 2) && (tf != 0))
			tf *= -1.0;
		if (m_bFirst) {
			g_fCubeXVector[k] = tf;
		} else {
			g_fCubeXVector[k] = g_mCubeCell(0, k) / (rx * g_iCubeXStride - g_iCubeXMismatch) / LEN_AU2PM;
			if (fabs(g_fCubeXVector[k] - tf / g_iCubeXStride) > 0.1) {
				eprintf("CTimeStep::ReadCube(): Error: Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeXVector[k] * LEN_AU2PM, tf / g_iCubeXStride * LEN_AU2PM);
				return false;
			}
		}
	}
	g_fCubeXStep = g_fCubeXVector[0];
	
// 		p++;
// 		while (*p == ' ')
// 			p++;
// 		q = p;
// 		while ((*p != ' ') && (*p != 0))
// 			p++;
// 	// 	if (*p == 0)
// 	// 	{
// 	// 		eprintf("Incomplete line while reading grid data from cube file.\n");
// 	// 		return false;
// 	// 	}
// 		*p = 0;
// 		
// 		if (first) {
// 			g_fCubeXStep = atof(q);
// 		} else {
// 			g_fCubeXStep = g_fBoxX / (rx * g_iCubeXStride - g_iCubeXMismatch) / LEN_AU2PM;
// 			if (fabs(g_fCubeXStep - atof(q) / g_iCubeXStride) > 0.1f) {
// 				eprintf("Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeXStep * LEN_AU2PM, atof(q) / g_iCubeXStride * LEN_AU2PM);
// 				return false;
// 			}
// 		}
	
	(void)!fgets(buf,256,a);
	if (feof(a)) { eprintf("CTimeStep::ReadCube(): Error: Unexpected end of file reached while reading cube file.\n"); return false; }
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;
	ry = atoi(buf);
	g_iCubeRes[1] = ry;

	for (k = 0; k < 3; k++) {
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		*p = 0;
		
		tf = atof(q);
		if (g_bFlipCellVectorX && (k == 0) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorY && (k == 1) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorZ && (k == 2) && (tf != 0))
			tf *= -1.0;
		if (m_bFirst) {
			g_fCubeYVector[k] = tf;
		} else {
			g_fCubeYVector[k] = g_mCubeCell(1, k) / (ry * g_iCubeYStride - g_iCubeYMismatch) / LEN_AU2PM;
			if (fabs(g_fCubeYVector[k] - tf / g_iCubeYStride) > 0.1) {
				eprintf("CTimeStep::ReadCube(): Error: Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeYVector[k] * LEN_AU2PM, tf / g_iCubeYStride * LEN_AU2PM);
				return false;
			}
		}
	}
	g_fCubeYStep = g_fCubeYVector[1];
	
// 		for (z = 0; z < 2; z++) {
// 			p++;
// 			while (*p == ' ')
// 				p++;
// 			q = p;
// 			while ((*p != ' ') && (*p != 0))
// 				p++;
// 	// 		if (*p == 0)
// 	// 		{
// 	// 			eprintf("Incomplete line while reading grid data from cube file.\n");
// 	// 			return false;
// 	// 		}
// 			*p = 0;
// 		}
// 		
// 		if (first) {
// 			g_fCubeYStep = atof(q);
// 		} else {
// 			g_fCubeYStep = g_fBoxY / (ry * g_iCubeYStride - g_iCubeYMismatch) / LEN_AU2PM;
// 			if (fabs(g_fCubeYStep - atof(q) / g_iCubeYStride) > 0.1f) {
// 				eprintf("Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeYStep * LEN_AU2PM, atof(q) / g_iCubeYStride * LEN_AU2PM);
// 				return false;
// 			}
// 		}
	
	(void)!fgets(buf,256,a);
	if (feof(a)) { eprintf("CTimeStep::ReadCube(): Error: Unexpected end of file reached while reading cube file.\n"); return false; }
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;
	rz = atoi(buf);
	g_iCubeRes[2] = rz;

	for (k = 0; k < 3; k++) {
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		*p = 0;
		
		tf = atof(q);
		if (g_bFlipCellVectorX && (k == 0) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorY && (k == 1) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorZ && (k == 2) && (tf != 0))
			tf *= -1.0;
		if (m_bFirst) {
			g_fCubeZVector[k] = tf;
		} else {
			g_fCubeZVector[k] = g_mCubeCell(2, k) / (rz * g_iCubeZStride - g_iCubeZMismatch) / LEN_AU2PM;
			if (fabs(g_fCubeZVector[k] - tf / g_iCubeZStride) > 0.1) {
				eprintf("CTimeStep::ReadCube(): Error: Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeZVector[k] * LEN_AU2PM, tf / g_iCubeZStride * LEN_AU2PM);
				return false;
			}
		}
	}
	g_fCubeZStep = g_fCubeZVector[2];
	m_bFirst = false;
	
// 		for (z = 0; z < 3; z++) {
// 			p++;
// 			while (*p == ' ')
// 				p++;
// 			q = p;
// 			while ((*p != ' ') && (*p != 0))
// 				p++;
// 	// 		if (*p == 0)
// 	// 		{
// 	// 			eprintf("Incomplete line while reading grid data from cube file.\n");
// 	// 			return false;
// 	// 		}
// 			*p = 0;
// 		}
// 		
// 		if (first) {
// 			g_fCubeZStep = atof(q);
// 			first = false;
// 		} else {
// 			g_fCubeZStep = g_fBoxZ / (rz * g_iCubeZStride - g_iCubeZMismatch) / LEN_AU2PM;
// 			if (fabs(g_fCubeZStep - atof(q) / g_iCubeZStride) > 0.1f) {
// 				eprintf("Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeZStep * LEN_AU2PM, atof(q) / g_iCubeZStride * LEN_AU2PM);
// 				return false;
// 			}
// 		}
	
	if (needinfo)
	{
		for (z=0;z<m_paLabels.GetSize();z++)
			delete[] (char*)m_paLabels[z];
		m_paLabels.RemoveAll();
		if (g_bDoubleBox)
			m_paLabels.SetSize(m_iGesAtomCount*g_iDoubleBoxFactor);
		else
			m_paLabels.SetSize(m_iGesAtomCount);
	}

	m_vaCoords.RemoveAll_KeepSize();
	CxIntArray atomNumbers;
	atomNumbers.SetSize(ac);
	for (z=0;z<ac;z++)
	{
		(void)!fgets(buf,256,a);
		if (feof(a)) { eprintf("CTimeStep::ReadCube(): Error: Unexpected end of file reached while reading cube file.\n"); return false; }
		buf[strlen(buf)-1] = 0;

		p = buf;
		while (*p == ' ')
			p++;
		while (*p != ' ')
			p++;
		*p = 0;

		atomNumbers[z] = atoi(buf);
		if (needinfo)
		{
			i = atoi(buf);

			for (z2=0;z2<g_oaElements.GetSize();z2++)
			{
				if (((CElement*)g_oaElements[z2])->m_iOrd == i)
				{
					q = ((CElement*)g_oaElements[z2])->m_sLabel;
					goto _found;
				}
			}
			eprintf("CTimeStep::ReadCube(): Error: Encountered unknown atom type %d while reading cube file.\n",i);
			return false;

_found:
			try { r = new char[strlen(q)+1]; } catch(...) { r = NULL; }
			if (r == NULL) NewException((double)(strlen(q)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(r,q);
			m_paLabels[z] = r;
		}

		p++; // Column is always zero
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		if (*p == 0)
		{
			eprintf("CTimeStep::ReadCube(): Error: Incomplete line while reading coordinates from cube file (A).\n");
			return false;
		}
		*p = 0;

		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		if (*p == 0)
		{
			eprintf("CTimeStep::ReadCube(): Error: Incomplete line while reading coordinates from cube file (B).\n");
			return false;
		}
		*p = 0;
		fx = atof(q);

		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		if (*p == 0)
		{
			eprintf("CTimeStep::ReadCube(): Error: Incomplete line while reading coordinates from cube file (C).\n");
			return false;
		}
		fy = atof(q);

		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		*p = 0;
		fz = atof(q);

		m_vaCoords.Add(CxDVector3(fx*LEN_AU2PM,fy*LEN_AU2PM,fz*LEN_AU2PM));
	}

	if (feof(a)) { eprintf("CTimeStep::ReadCube(): Error: Unexpected end of file reached while reading cube file.\n"); return false; }

	if (m_pVolumetricData == NULL)
	{
		m_pVolumetricData = new C3DF<VORI_FLOAT>();
	}
	
	if (m_pVolumetricDataTimeDev == NULL && g_bCubeTimeDev) {
		m_pVolumetricDataTimeDev = new C3DF<VORI_FLOAT>();
	}
	if (m_pCurrentDensity == NULL && g_bCubeTimeDev) {
		m_pCurrentDensity = new CxDoubleArray();
	}

	if (g_bVerbose) {
		mprintf("\n");
		mprintf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		mprintf(">>>   CTimeStep::ReadCube()   >>>\n");
		mprintf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		mprintf("    Resolution %d x %d x %d\n",g_iCubeRes[0],g_iCubeRes[1],g_iCubeRes[2]);
		mprintf("    Stride %d x %d x %d\n",g_iCubeXStride,g_iCubeYStride,g_iCubeZStride);
		mprintf("    Mismatch %d x %d x %d\n",g_iCubeXMismatch,g_iCubeYMismatch,g_iCubeZMismatch);
		mprintf("    Resulting Resolution %d x %d x %d\n",rx * g_iCubeXStride - g_iCubeXMismatch,ry * g_iCubeYStride - g_iCubeYMismatch,rz * g_iCubeZStride - g_iCubeZMismatch);
		mprintf("    Cell Vector A  ( %12.6f | %12.6f | %12.6f ) pm\n",
			g_fCubeXVector[0] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch),
			g_fCubeXVector[1] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch),
			g_fCubeXVector[2] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch));
		mprintf("    Cell Vector B  ( %12.6f | %12.6f | %12.6f ) pm\n",
			g_fCubeYVector[0] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch),
			g_fCubeYVector[1] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch),
			g_fCubeYVector[2] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch));
		mprintf("    Cell Vector C  ( %12.6f | %12.6f | %12.6f ) pm\n",
			g_fCubeZVector[0] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch),
			g_fCubeZVector[1] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch),
			g_fCubeZVector[2] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch));
		mprintf("    First 10 atoms:\n");
		for (z=0;z<((ac<10)?ac:10);z++)
			mprintf("      %3d  %12.6f  %12.6f  %12.6f\n",atomNumbers[z],m_vaCoords[z][0],m_vaCoords[z][1],m_vaCoords[z][2]);
	}

	static bool second = false;
	if (m_pVolumetricData->m_pBin == NULL)
	{
		m_pVolumetricData->m_iRes[0] = rx * g_iCubeXStride - g_iCubeXMismatch;
		m_pVolumetricData->m_iRes[1] = ry * g_iCubeYStride - g_iCubeYMismatch;
		m_pVolumetricData->m_iRes[2] = rz * g_iCubeZStride - g_iCubeZMismatch;
		
		m_pVolumetricData->m_fMinVal[0] = 0.0;
		m_pVolumetricData->m_fMinVal[1] = 0.0;
		m_pVolumetricData->m_fMinVal[2] = 0.0;
		g_mCubeCell(0, 0) = g_fCubeXVector[0] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
		g_mCubeCell(0, 1) = g_fCubeXVector[1] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
		g_mCubeCell(0, 2) = g_fCubeXVector[2] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
		g_mCubeCell(1, 0) = g_fCubeYVector[0] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
		g_mCubeCell(1, 1) = g_fCubeYVector[1] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
		g_mCubeCell(1, 2) = g_fCubeYVector[2] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
		g_mCubeCell(2, 0) = g_fCubeZVector[0] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
		g_mCubeCell(2, 1) = g_fCubeZVector[1] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
		g_mCubeCell(2, 2) = g_fCubeZVector[2] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
		m_pVolumetricData->m_fMaxVal[0] = sqrt(g_mCubeCell(0, 0) * g_mCubeCell(0, 0) + g_mCubeCell(0, 1) * g_mCubeCell(0, 1) + g_mCubeCell(0, 2) * g_mCubeCell(0, 2));
		m_pVolumetricData->m_fMaxVal[1] = sqrt(g_mCubeCell(1, 0) * g_mCubeCell(1, 0) + g_mCubeCell(1, 1) * g_mCubeCell(1, 1) + g_mCubeCell(1, 2) * g_mCubeCell(1, 2));
		m_pVolumetricData->m_fMaxVal[2] = sqrt(g_mCubeCell(2, 0) * g_mCubeCell(2, 0) + g_mCubeCell(2, 1) * g_mCubeCell(2, 1) + g_mCubeCell(2, 2) * g_mCubeCell(2, 2));
// 	// 		m_pVolumetricData->m_fMaxVal[0] = g_fBoxX;
// 			m_pVolumetricData->m_fMaxVal[0] = g_fCubeXStep * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
// 	// 		m_pVolumetricData->m_fMaxVal[1] = g_fBoxY;
// 			m_pVolumetricData->m_fMaxVal[1] = g_fCubeYStep * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
// 	// 		m_pVolumetricData->m_fMaxVal[2] = g_fBoxZ;
// 			m_pVolumetricData->m_fMaxVal[2] = g_fCubeZStep * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
		
		m_pVolumetricData->Create();
		
		if (m_pVolumetricDataTimeDev != NULL) {
// 			m_pVolumetricDataTimeDev->m_iRes[0] = rx * g_iCubeXStride - g_iCubeXMismatch;
// 			m_pVolumetricDataTimeDev->m_iRes[1] = ry * g_iCubeYStride - g_iCubeYMismatch;
// 			m_pVolumetricDataTimeDev->m_iRes[2] = rz * g_iCubeZStride - g_iCubeZMismatch;
// 			m_pVolumetricDataTimeDev->m_fMinVal[0] = 0;
// 			m_pVolumetricDataTimeDev->m_fMaxVal[0] = g_fCubeXStep * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
// 			m_pVolumetricDataTimeDev->m_fMinVal[1] = 0;
// 			m_pVolumetricDataTimeDev->m_fMaxVal[1] = g_fCubeYStep * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
// 			m_pVolumetricDataTimeDev->m_fMinVal[2] = 0;
// 			m_pVolumetricDataTimeDev->m_fMaxVal[2] = g_fCubeZStep * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
			m_pVolumetricDataTimeDev->m_iRes[0] = m_pVolumetricData->m_iRes[0];
			m_pVolumetricDataTimeDev->m_iRes[1] = m_pVolumetricData->m_iRes[1];
			m_pVolumetricDataTimeDev->m_iRes[2] = m_pVolumetricData->m_iRes[2];
			m_pVolumetricDataTimeDev->m_fMinVal[0] = m_pVolumetricData->m_fMinVal[0];
			m_pVolumetricDataTimeDev->m_fMinVal[1] = m_pVolumetricData->m_fMinVal[1];
			m_pVolumetricDataTimeDev->m_fMinVal[2] = m_pVolumetricData->m_fMinVal[2];
			m_pVolumetricDataTimeDev->m_fMaxVal[0] = m_pVolumetricData->m_fMaxVal[0];
			m_pVolumetricDataTimeDev->m_fMaxVal[1] = m_pVolumetricData->m_fMaxVal[1];
			m_pVolumetricDataTimeDev->m_fMaxVal[2] = m_pVolumetricData->m_fMaxVal[2];
			
			m_pVolumetricDataTimeDev->Create();
		}
		if (m_pCurrentDensity != NULL) {
			m_pCurrentDensity->SetSize(3 * m_pVolumetricDataTimeDev->m_iRes[0] * m_pVolumetricDataTimeDev->m_iRes[1] * m_pVolumetricDataTimeDev->m_iRes[2]);
		}
		
		second = true;
	} else
	{
		if ((rx * g_iCubeXStride - g_iCubeXMismatch != m_pVolumetricData->m_iRes[0]) || (ry * g_iCubeYStride - g_iCubeYMismatch!= m_pVolumetricData->m_iRes[1]) || (rz * g_iCubeZStride - g_iCubeZMismatch != m_pVolumetricData->m_iRes[2]))
		{
			if (second) {
				second = false;
				
				m_pVolumetricData->m_iRes[0] = rx * g_iCubeXStride - g_iCubeXMismatch;
				m_pVolumetricData->m_iRes[1] = ry * g_iCubeYStride - g_iCubeYMismatch;
				m_pVolumetricData->m_iRes[2] = rz * g_iCubeZStride - g_iCubeZMismatch;
				
				m_pVolumetricData->m_fMinVal[0] = 0.0;
				m_pVolumetricData->m_fMinVal[1] = 0.0;
				m_pVolumetricData->m_fMinVal[2] = 0.0;
				g_mCubeCell(0, 0) = g_fCubeXVector[0] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
				g_mCubeCell(0, 1) = g_fCubeXVector[1] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
				g_mCubeCell(0, 2) = g_fCubeXVector[2] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
				g_mCubeCell(1, 0) = g_fCubeYVector[0] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
				g_mCubeCell(1, 1) = g_fCubeYVector[1] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
				g_mCubeCell(1, 2) = g_fCubeYVector[2] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
				g_mCubeCell(2, 0) = g_fCubeZVector[0] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
				g_mCubeCell(2, 1) = g_fCubeZVector[1] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
				g_mCubeCell(2, 2) = g_fCubeZVector[2] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
				m_pVolumetricData->m_fMaxVal[0] = sqrt(g_mCubeCell(0, 0) * g_mCubeCell(0, 0) + g_mCubeCell(0, 1) * g_mCubeCell(0, 1) + g_mCubeCell(0, 2) * g_mCubeCell(0, 2));
				m_pVolumetricData->m_fMaxVal[1] = sqrt(g_mCubeCell(1, 0) * g_mCubeCell(1, 0) + g_mCubeCell(1, 1) * g_mCubeCell(1, 1) + g_mCubeCell(1, 2) * g_mCubeCell(1, 2));
				m_pVolumetricData->m_fMaxVal[2] = sqrt(g_mCubeCell(2, 0) * g_mCubeCell(2, 0) + g_mCubeCell(2, 1) * g_mCubeCell(2, 1) + g_mCubeCell(2, 2) * g_mCubeCell(2, 2));
// 					// 		m_pVolumetricData->m_fMaxVal[0] = g_fBoxX;
// 					m_pVolumetricData->m_fMaxVal[0] = g_fCubeXStep * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
// 					// 		m_pVolumetricData->m_fMaxVal[1] = g_fBoxY;
// 					m_pVolumetricData->m_fMaxVal[1] = g_fCubeYStep * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
// 					// 		m_pVolumetricData->m_fMaxVal[2] = g_fBoxZ;
// 					m_pVolumetricData->m_fMaxVal[2] = g_fCubeZStep * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
				
				m_pVolumetricData->Create();
				
				if (m_pVolumetricDataTimeDev != NULL) {
// 					m_pVolumetricDataTimeDev->m_iRes[0] = rx * g_iCubeXStride - g_iCubeXMismatch;
// 					m_pVolumetricDataTimeDev->m_iRes[1] = ry * g_iCubeYStride - g_iCubeYMismatch;
// 					m_pVolumetricDataTimeDev->m_iRes[2] = rz * g_iCubeZStride - g_iCubeZMismatch;
// 					m_pVolumetricDataTimeDev->m_fMinVal[0] = 0;
// 					m_pVolumetricDataTimeDev->m_fMaxVal[0] = g_fCubeXStep * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
// 					m_pVolumetricDataTimeDev->m_fMinVal[1] = 0;
// 					m_pVolumetricDataTimeDev->m_fMaxVal[1] = g_fCubeYStep * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
// 					m_pVolumetricDataTimeDev->m_fMinVal[2] = 0;
// 					m_pVolumetricDataTimeDev->m_fMaxVal[2] = g_fCubeZStep * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
					m_pVolumetricDataTimeDev->m_iRes[0] = m_pVolumetricData->m_iRes[0];
					m_pVolumetricDataTimeDev->m_iRes[1] = m_pVolumetricData->m_iRes[1];
					m_pVolumetricDataTimeDev->m_iRes[2] = m_pVolumetricData->m_iRes[2];
					m_pVolumetricDataTimeDev->m_fMinVal[0] = m_pVolumetricData->m_fMinVal[0];
					m_pVolumetricDataTimeDev->m_fMinVal[1] = m_pVolumetricData->m_fMinVal[1];
					m_pVolumetricDataTimeDev->m_fMinVal[2] = m_pVolumetricData->m_fMinVal[2];
					m_pVolumetricDataTimeDev->m_fMaxVal[0] = m_pVolumetricData->m_fMaxVal[0];
					m_pVolumetricDataTimeDev->m_fMaxVal[1] = m_pVolumetricData->m_fMaxVal[1];
					m_pVolumetricDataTimeDev->m_fMaxVal[2] = m_pVolumetricData->m_fMaxVal[2];

					m_pVolumetricDataTimeDev->Create();
				}
				if (m_pCurrentDensity != NULL) {
					m_pCurrentDensity->SetSize(3 * m_pVolumetricDataTimeDev->m_iRes[0] * m_pVolumetricDataTimeDev->m_iRes[1] * m_pVolumetricDataTimeDev->m_iRes[2]);
				}
			} else {
				eprintf("\nCTimeStep::ReadCube(): Error: Cube file dimension mismatch (%d-%d, %d-%d, %d-%d).\n",rx * g_iCubeXStride - g_iCubeXMismatch,m_pVolumetricData->m_iRes[0],ry * g_iCubeYStride - g_iCubeYMismatch,m_pVolumetricData->m_iRes[1],rz * g_iCubeZStride - g_iCubeZMismatch,m_pVolumetricData->m_iRes[2]);
				return false;
			}
		}
		for (z=0;z<m_pVolumetricData->m_iRes[0]*m_pVolumetricData->m_iRes[1]*m_pVolumetricData->m_iRes[2];z++)
			m_pVolumetricData->m_pBin[z] = 0;
	}

	m_pVolumetricData->ReadCubeData(a,false);

	if (g_bVerbose) {
		mprintf("    First 100 volumetric data elements:");
		for (z=0;z<((m_pVolumetricData->m_iRes[0]*m_pVolumetricData->m_iRes[1]*m_pVolumetricData->m_iRes[2]<100)?m_pVolumetricData->m_iRes[0]*m_pVolumetricData->m_iRes[1]*m_pVolumetricData->m_iRes[2]:100);z++) {
			if ((z%10) == 0)
				mprintf("\n   ");
			mprintf("   %10G",m_pVolumetricData->m_pBin[z]);
		}
		mprintf("\n");
		mprintf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");
	}
	
	
/* 	FILE *cubeFile = fopen("test.cube", "w");
 	if (cubeFile == NULL) {
 		printf("Could not open output file!\n");
 		return 1;
 	}
 	
 	fprintf(cubeFile, "\n\n");
 	fprintf(cubeFile, "%5lu %12.6f %12.6f %12.6f\n", m_iGesAtomCount, 0.0f, 0.0f, 0.0f);
 	fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", m_pVolumetricData->m_iRes[0], g_fCubeXStep, 0.0f, 0.0f);
 	fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", m_pVolumetricData->m_iRes[1], 0.0f, g_fCubeYStep, 0.0f);
 	fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", m_pVolumetricData->m_iRes[2], 0.0f, 0.0f, g_fCubeZStep);
 	for (int i = 0; i < (int)m_iGesAtomCount; i++) {
 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f %12.6f\n", 1, 0.0f, m_vaCoords[i][0] / 0.529177249f / 100.0f, m_vaCoords[i][1] / 0.529177249f / 100.0f, m_vaCoords[i][2] / 0.529177249f / 100.0f);
 	}
 	for (int i = 0; i < m_pVolumetricData->m_iRes[0]; i++) {
 		for (int j = 0; j < m_pVolumetricData->m_iRes[1]; j++) {
 			for (int k = 0; k < m_pVolumetricData->m_iRes[2] / 6; k++) {
 				for (int l = 0; l < 6; l++) {
 					fprintf(cubeFile, "%13.5E", m_pVolumetricData->m_pBin[i + j * m_pVolumetricData->m_iRes[0] + (k * 6 + l) * m_pVolumetricData->m_iRes[0] * m_pVolumetricData->m_iRes[1]]);
 				}
 				fprintf(cubeFile, "\n");
 			}
 			if (m_pVolumetricData->m_iRes[2] % 6 != 0) {
 				for (int l = 0; l < m_pVolumetricData->m_iRes[2] % 6; l++) {
 					fprintf(cubeFile, "%13.5E", m_pVolumetricData->m_pBin[i + j * m_pVolumetricData->m_iRes[0] + ((m_pVolumetricData->m_iRes[2] / 6) * 6 + l) * m_pVolumetricData->m_iRes[0] * m_pVolumetricData->m_iRes[1]]);
 				}
 				fprintf(cubeFile, "\n");
 			}
 		}
 	}
 	
 	fclose(cubeFile);*/
	
	return true;
}


bool CTimeStep::ReadBQBCubeStep(bool needinfo)
{
	int rx, ry, rz, ac, z, z2;
	char *r, i;
	const char *q;
	double fx, fy, fz, tf;
	CBQBCubeFrame *cfr;
	CBQBAtomSet *afr;


	if (g_bVerbose)
		mprintf(">>> CTimeStep::ReadBQBCubeStep() >>>\n");

	cfr = g_pBQBFile->GetEngine()->GetOutputCubeFrame(0);
	afr = g_pBQBFile->GetEngine()->GetOutputAtomFrame(0);

	ac = (int)afr->m_oaAtoms.size();
	m_iGesAtomCount = ac;

	if (ac <= 0) {
		eprintf("CTimeStep::ReadBQBCubeStep(): Error: Encountered atom count of %d while reading cube file.\n",ac);
		return false;
	}

	rx = cfr->m_iRes[0];
	g_iCubeRes[0] = rx;
	
	//static bool first = true;
	// replaced by m_bFirst
	
	int k;
	for (k = 0; k < 3; k++) {
		tf = cfr->m_fStrideA[k];
		if (g_bFlipCellVectorX && (k == 0) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorY && (k == 1) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorZ && (k == 2) && (tf != 0))
			tf *= -1.0;
		if (m_bFirst)
			g_fCubeXVector[k] = tf;
		else {
			g_fCubeXVector[k] = g_mCubeCell(0, k) / (rx * g_iCubeXStride - g_iCubeXMismatch) / LEN_AU2PM;
			if (fabs(g_fCubeXVector[k] - tf / g_iCubeXStride) > 0.1) {
				eprintf("CTimeStep::ReadBQBCubeStep(): Error: Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeXVector[k] * LEN_AU2PM, tf / g_iCubeXStride * LEN_AU2PM);
				return false;
			}
		}
	}
	g_fCubeXStep = g_fCubeXVector[0];
	
	ry = cfr->m_iRes[1];
	g_iCubeRes[1] = ry;

	for (k = 0; k < 3; k++) {
		tf = cfr->m_fStrideB[k];
		if (g_bFlipCellVectorX && (k == 0) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorY && (k == 1) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorZ && (k == 2) && (tf != 0))
			tf *= -1.0;
		if (m_bFirst)
			g_fCubeYVector[k] = tf;
		else {
			g_fCubeYVector[k] = g_mCubeCell(1, k) / (ry * g_iCubeYStride - g_iCubeYMismatch) / LEN_AU2PM;
			if (fabs(g_fCubeYVector[k] - tf / g_iCubeYStride) > 0.1) {
				eprintf("CTimeStep::ReadBQBCubeStep(): Error: Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeYVector[k] * LEN_AU2PM, tf / g_iCubeYStride * LEN_AU2PM);
				return false;
			}
		}
	}
	g_fCubeYStep = g_fCubeYVector[1];
	
	rz = cfr->m_iRes[2];
	g_iCubeRes[2] = rz;

	for (k = 0; k < 3; k++) {
		tf = cfr->m_fStrideC[k];
		if (g_bFlipCellVectorX && (k == 0) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorY && (k == 1) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorZ && (k == 2) && (tf != 0))
			tf *= -1.0;
		if (m_bFirst)
			g_fCubeZVector[k] = tf;
		else {
			g_fCubeZVector[k] = g_mCubeCell(2, k) / (rz * g_iCubeZStride - g_iCubeZMismatch) / LEN_AU2PM;
			if (fabs(g_fCubeZVector[k] - tf / g_iCubeZStride) > 0.1) {
				eprintf("CTimeStep::ReadBQBCubeStep(): Error: Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeZVector[k] * LEN_AU2PM, tf / g_iCubeZStride * LEN_AU2PM);
				return false;
			}
		}
	}
	g_fCubeZStep = g_fCubeZVector[2];
	m_bFirst = false;
	
	if (needinfo)
	{
		for (z=0;z<m_paLabels.GetSize();z++)
			delete[] (char*)m_paLabels[z];
		m_paLabels.RemoveAll();
		if (g_bDoubleBox)
			m_paLabels.SetSize(m_iGesAtomCount*g_iDoubleBoxFactor);
		else
			m_paLabels.SetSize(m_iGesAtomCount);
	}

	m_vaCoords.RemoveAll_KeepSize();
	CxIntArray atomNumbers;
	atomNumbers.SetSize(ac);
	for (z=0;z<ac;z++)
	{
		atomNumbers[z] = afr->m_oaAtoms[z]->m_iOrd;
		if (needinfo)
		{
			i = (char)afr->m_oaAtoms[z]->m_iOrd;

			for (z2=0;z2<g_oaElements.GetSize();z2++)
			{
				if (((CElement*)g_oaElements[z2])->m_iOrd == i)
				{
					q = ((CElement*)g_oaElements[z2])->m_sLabel;
					goto _found;
				}
			}
			eprintf("CTimeStep::ReadBQBCubeStep(): Error: Encountered unknown atom type %d while reading cube file.\n",i);
			return false;

_found:
			try { r = new char[strlen(q)+1]; } catch(...) { r = NULL; }
			if (r == NULL) NewException((double)(strlen(q)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(r,q);
			m_paLabels[z] = r;
		}

		fx = afr->m_oaAtoms[z]->m_fCoord[0];
		fy = afr->m_oaAtoms[z]->m_fCoord[1];
		fz = afr->m_oaAtoms[z]->m_fCoord[2];

	//	if (z < 10)
	//		mprintf("@%d: %f | %f | %f --> %f | %f | %f\n",z,fx,fy,fz,fx*LEN_AU2PM,fy*LEN_AU2PM,fz*LEN_AU2PM);

		m_vaCoords.Add(CxDVector3(fx*LEN_AU2PM,fy*LEN_AU2PM,fz*LEN_AU2PM));
	}

	if (m_pVolumetricData == NULL)
		m_pVolumetricData = new C3DF<VORI_FLOAT>();
	
	if (m_pVolumetricDataTimeDev == NULL && g_bCubeTimeDev)
		m_pVolumetricDataTimeDev = new C3DF<VORI_FLOAT>();

	if (m_pCurrentDensity == NULL && g_bCubeTimeDev)
		m_pCurrentDensity = new CxDoubleArray();

	if (g_bVerbose) {
		mprintf("\n");
		mprintf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		mprintf(">>>   CTimeStep::ReadBQBCubeStep()   >>>\n");
		mprintf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		mprintf("    Resolution %d x %d x %d\n",g_iCubeRes[0],g_iCubeRes[1],g_iCubeRes[2]);
		mprintf("    Stride %d x %d x %d\n",g_iCubeXStride,g_iCubeYStride,g_iCubeZStride);
		mprintf("    Mismatch %d x %d x %d\n",g_iCubeXMismatch,g_iCubeYMismatch,g_iCubeZMismatch);
		mprintf("    Resulting Resolution %d x %d x %d\n",rx * g_iCubeXStride - g_iCubeXMismatch,ry * g_iCubeYStride - g_iCubeYMismatch,rz * g_iCubeZStride - g_iCubeZMismatch);
		mprintf("    Cell Vector A  ( %12.6f | %12.6f | %12.6f ) pm\n",
			g_fCubeXVector[0] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch),
			g_fCubeXVector[1] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch),
			g_fCubeXVector[2] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch));
		mprintf("    Cell Vector B  ( %12.6f | %12.6f | %12.6f ) pm\n",
			g_fCubeYVector[0] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch),
			g_fCubeYVector[1] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch),
			g_fCubeYVector[2] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch));
		mprintf("    Cell Vector C  ( %12.6f | %12.6f | %12.6f ) pm\n",
			g_fCubeZVector[0] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch),
			g_fCubeZVector[1] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch),
			g_fCubeZVector[2] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch));
		mprintf("    First 10 atoms:\n");
		for (z=0;z<((ac<10)?ac:10);z++)
			mprintf("      %3d  %12.6f  %12.6f  %12.6f\n",atomNumbers[z],m_vaCoords[z][0],m_vaCoords[z][1],m_vaCoords[z][2]);
	}

	static bool second = false;
	if (m_pVolumetricData->m_pBin == NULL)
	{
		m_pVolumetricData->m_iRes[0] = rx * g_iCubeXStride - g_iCubeXMismatch;
		m_pVolumetricData->m_iRes[1] = ry * g_iCubeYStride - g_iCubeYMismatch;
		m_pVolumetricData->m_iRes[2] = rz * g_iCubeZStride - g_iCubeZMismatch;
		
		m_pVolumetricData->m_fMinVal[0] = 0.0;
		m_pVolumetricData->m_fMinVal[1] = 0.0;
		m_pVolumetricData->m_fMinVal[2] = 0.0;
		g_mCubeCell(0, 0) = g_fCubeXVector[0] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
		g_mCubeCell(0, 1) = g_fCubeXVector[1] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
		g_mCubeCell(0, 2) = g_fCubeXVector[2] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
		g_mCubeCell(1, 0) = g_fCubeYVector[0] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
		g_mCubeCell(1, 1) = g_fCubeYVector[1] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
		g_mCubeCell(1, 2) = g_fCubeYVector[2] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
		g_mCubeCell(2, 0) = g_fCubeZVector[0] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
		g_mCubeCell(2, 1) = g_fCubeZVector[1] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
		g_mCubeCell(2, 2) = g_fCubeZVector[2] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
		m_pVolumetricData->m_fMaxVal[0] = sqrt(g_mCubeCell(0, 0) * g_mCubeCell(0, 0) + g_mCubeCell(0, 1) * g_mCubeCell(0, 1) + g_mCubeCell(0, 2) * g_mCubeCell(0, 2));
		m_pVolumetricData->m_fMaxVal[1] = sqrt(g_mCubeCell(1, 0) * g_mCubeCell(1, 0) + g_mCubeCell(1, 1) * g_mCubeCell(1, 1) + g_mCubeCell(1, 2) * g_mCubeCell(1, 2));
		m_pVolumetricData->m_fMaxVal[2] = sqrt(g_mCubeCell(2, 0) * g_mCubeCell(2, 0) + g_mCubeCell(2, 1) * g_mCubeCell(2, 1) + g_mCubeCell(2, 2) * g_mCubeCell(2, 2));
		
		m_pVolumetricData->Create();
		
		if (m_pVolumetricDataTimeDev != NULL) {
			m_pVolumetricDataTimeDev->m_iRes[0] = m_pVolumetricData->m_iRes[0];
			m_pVolumetricDataTimeDev->m_iRes[1] = m_pVolumetricData->m_iRes[1];
			m_pVolumetricDataTimeDev->m_iRes[2] = m_pVolumetricData->m_iRes[2];
			m_pVolumetricDataTimeDev->m_fMinVal[0] = m_pVolumetricData->m_fMinVal[0];
			m_pVolumetricDataTimeDev->m_fMinVal[1] = m_pVolumetricData->m_fMinVal[1];
			m_pVolumetricDataTimeDev->m_fMinVal[2] = m_pVolumetricData->m_fMinVal[2];
			m_pVolumetricDataTimeDev->m_fMaxVal[0] = m_pVolumetricData->m_fMaxVal[0];
			m_pVolumetricDataTimeDev->m_fMaxVal[1] = m_pVolumetricData->m_fMaxVal[1];
			m_pVolumetricDataTimeDev->m_fMaxVal[2] = m_pVolumetricData->m_fMaxVal[2];
			
			m_pVolumetricDataTimeDev->Create();
		}
		if (m_pCurrentDensity != NULL) {
			m_pCurrentDensity->SetSize(3 * m_pVolumetricDataTimeDev->m_iRes[0] * m_pVolumetricDataTimeDev->m_iRes[1] * m_pVolumetricDataTimeDev->m_iRes[2]);
		}
		
		second = true;
	} else
	{
		if ((rx * g_iCubeXStride - g_iCubeXMismatch != m_pVolumetricData->m_iRes[0]) || (ry * g_iCubeYStride - g_iCubeYMismatch!= m_pVolumetricData->m_iRes[1]) || (rz * g_iCubeZStride - g_iCubeZMismatch != m_pVolumetricData->m_iRes[2]))
		{
			if (second) {
				second = false;
				
				m_pVolumetricData->m_iRes[0] = rx * g_iCubeXStride - g_iCubeXMismatch;
				m_pVolumetricData->m_iRes[1] = ry * g_iCubeYStride - g_iCubeYMismatch;
				m_pVolumetricData->m_iRes[2] = rz * g_iCubeZStride - g_iCubeZMismatch;
				
				m_pVolumetricData->m_fMinVal[0] = 0.0;
				m_pVolumetricData->m_fMinVal[1] = 0.0;
				m_pVolumetricData->m_fMinVal[2] = 0.0;
				g_mCubeCell(0, 0) = g_fCubeXVector[0] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
				g_mCubeCell(0, 1) = g_fCubeXVector[1] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
				g_mCubeCell(0, 2) = g_fCubeXVector[2] * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
				g_mCubeCell(1, 0) = g_fCubeYVector[0] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
				g_mCubeCell(1, 1) = g_fCubeYVector[1] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
				g_mCubeCell(1, 2) = g_fCubeYVector[2] * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
				g_mCubeCell(2, 0) = g_fCubeZVector[0] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
				g_mCubeCell(2, 1) = g_fCubeZVector[1] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
				g_mCubeCell(2, 2) = g_fCubeZVector[2] * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
				m_pVolumetricData->m_fMaxVal[0] = sqrt(g_mCubeCell(0, 0) * g_mCubeCell(0, 0) + g_mCubeCell(0, 1) * g_mCubeCell(0, 1) + g_mCubeCell(0, 2) * g_mCubeCell(0, 2));
				m_pVolumetricData->m_fMaxVal[1] = sqrt(g_mCubeCell(1, 0) * g_mCubeCell(1, 0) + g_mCubeCell(1, 1) * g_mCubeCell(1, 1) + g_mCubeCell(1, 2) * g_mCubeCell(1, 2));
				m_pVolumetricData->m_fMaxVal[2] = sqrt(g_mCubeCell(2, 0) * g_mCubeCell(2, 0) + g_mCubeCell(2, 1) * g_mCubeCell(2, 1) + g_mCubeCell(2, 2) * g_mCubeCell(2, 2));
				
				m_pVolumetricData->Create();
				
				if (m_pVolumetricDataTimeDev != NULL) {
					m_pVolumetricDataTimeDev->m_iRes[0] = m_pVolumetricData->m_iRes[0];
					m_pVolumetricDataTimeDev->m_iRes[1] = m_pVolumetricData->m_iRes[1];
					m_pVolumetricDataTimeDev->m_iRes[2] = m_pVolumetricData->m_iRes[2];
					m_pVolumetricDataTimeDev->m_fMinVal[0] = m_pVolumetricData->m_fMinVal[0];
					m_pVolumetricDataTimeDev->m_fMinVal[1] = m_pVolumetricData->m_fMinVal[1];
					m_pVolumetricDataTimeDev->m_fMinVal[2] = m_pVolumetricData->m_fMinVal[2];
					m_pVolumetricDataTimeDev->m_fMaxVal[0] = m_pVolumetricData->m_fMaxVal[0];
					m_pVolumetricDataTimeDev->m_fMaxVal[1] = m_pVolumetricData->m_fMaxVal[1];
					m_pVolumetricDataTimeDev->m_fMaxVal[2] = m_pVolumetricData->m_fMaxVal[2];

					m_pVolumetricDataTimeDev->Create();
				}
				if (m_pCurrentDensity != NULL) {
					m_pCurrentDensity->SetSize(3 * m_pVolumetricDataTimeDev->m_iRes[0] * m_pVolumetricDataTimeDev->m_iRes[1] * m_pVolumetricDataTimeDev->m_iRes[2]);
				}
			} else {
				eprintf("\nCTimeStep::ReadBQBCubeStep(): Error: Cube file dimension mismatch (%d-%d, %d-%d, %d-%d).\n",rx * g_iCubeXStride - g_iCubeXMismatch,m_pVolumetricData->m_iRes[0],ry * g_iCubeYStride - g_iCubeYMismatch,m_pVolumetricData->m_iRes[1],rz * g_iCubeZStride - g_iCubeZMismatch,m_pVolumetricData->m_iRes[2]);
				return false;
			}
		}
		for (z=0;z<m_pVolumetricData->m_iRes[0]*m_pVolumetricData->m_iRes[1]*m_pVolumetricData->m_iRes[2];z++)
			m_pVolumetricData->m_pBin[z] = 0;
	}

	m_pVolumetricData->ReadCubeData(cfr,g_bVerbose);

	if (g_bVerbose) {
		mprintf("    First 100 volumetric data elements:");
		for (z=0;z<((m_pVolumetricData->m_iRes[0]*m_pVolumetricData->m_iRes[1]*m_pVolumetricData->m_iRes[2]<100)?m_pVolumetricData->m_iRes[0]*m_pVolumetricData->m_iRes[1]*m_pVolumetricData->m_iRes[2]:100);z++) {
			if ((z%10) == 0)
				mprintf("\n   ");
			mprintf("   %10G",m_pVolumetricData->m_pBin[z]);
		}
		mprintf("\n");
		mprintf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");
	}

//	for (z=0;z<10;z++)
//		mprintf("@%d: %10G\n",z,m_pVolumetricData->m_pBin[z]);
	
	
/* 	FILE *cubeFile = fopen("test.cube", "w");
 	if (cubeFile == NULL) {
 		printf("Could not open output file!\n");
 		return 1;
 	}
 	
 	fprintf(cubeFile, "\n\n");
 	fprintf(cubeFile, "%5lu %12.6f %12.6f %12.6f\n", m_iGesAtomCount, 0.0f, 0.0f, 0.0f);
 	fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", m_pVolumetricData->m_iRes[0], g_fCubeXStep, 0.0f, 0.0f);
 	fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", m_pVolumetricData->m_iRes[1], 0.0f, g_fCubeYStep, 0.0f);
 	fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", m_pVolumetricData->m_iRes[2], 0.0f, 0.0f, g_fCubeZStep);
 	for (int i = 0; i < (int)m_iGesAtomCount; i++) {
 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f %12.6f\n", 1, 0.0f, m_vaCoords[i][0] / 0.529177249f / 100.0f, m_vaCoords[i][1] / 0.529177249f / 100.0f, m_vaCoords[i][2] / 0.529177249f / 100.0f);
 	}
 	for (int i = 0; i < m_pVolumetricData->m_iRes[0]; i++) {
 		for (int j = 0; j < m_pVolumetricData->m_iRes[1]; j++) {
 			for (int k = 0; k < m_pVolumetricData->m_iRes[2] / 6; k++) {
 				for (int l = 0; l < 6; l++) {
 					fprintf(cubeFile, "%13.5E", m_pVolumetricData->m_pBin[i + j * m_pVolumetricData->m_iRes[0] + (k * 6 + l) * m_pVolumetricData->m_iRes[0] * m_pVolumetricData->m_iRes[1]]);
 				}
 				fprintf(cubeFile, "\n");
 			}
 			if (m_pVolumetricData->m_iRes[2] % 6 != 0) {
 				for (int l = 0; l < m_pVolumetricData->m_iRes[2] % 6; l++) {
 					fprintf(cubeFile, "%13.5E", m_pVolumetricData->m_pBin[i + j * m_pVolumetricData->m_iRes[0] + ((m_pVolumetricData->m_iRes[2] / 6) * 6 + l) * m_pVolumetricData->m_iRes[0] * m_pVolumetricData->m_iRes[1]]);
 				}
 				fprintf(cubeFile, "\n");
 			}
 		}
 	}
 	
 	fclose(cubeFile);*/

	if (g_bVerbose)
		mprintf("<<< CTimeStep::ReadBQBCubeStep() <<<\n");
	
	return true;
}


bool CTimeStep::ReadBQBTrajStep(const std::vector<unsigned char> *data, bool needinfo) {

	int z;
	char *p;
	CBQBTrajectoryFrame btf(*g_pBQBInterface);
	CBQBTrajectoryFrameColumn *btfc, *btfcx, *btfcy, *btfcz, *btfcxx, *btfcxy, *btfcxz, *btfcyx, *btfcyy, *btfcyz, *btfczx, *btfczy, *btfczz, *velx, *vely, *velz;
	CxDVector3 veca, vecb, vecc;


	if (g_bVerbose)
		mprintf(">>> CTimeStep::ReadBQBTrajStep() >>>\n");

	if (!btf.ReadFrame(data)) {
		eprintf("CTimeStep::ReadBQBTrajStep(): Error: Could not read trajectory frame.\n");
		return false;
	}

	m_iGesAtomCount = btf.m_iAtomCount;

	/*if (m_iGesAtomCount < 0) {
		eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Atom count = %d < 0.\n",m_iGesAtomCount);
		m_iGesAtomCount = 0;
		return false;
	}*/
	if (m_iGesAtomCount == 0) {
		eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Atom count = 0.\n");
		return false;
	}
	if ((btfcx = btf.GetColumn("PosX")) == NULL) {
		eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"PosX\" does not exist in input trajectory.\n");
		return false;
	}
	if ((btfcx->m_iType != BQB_TYPE_FLOAT) && (btfcx->m_iType != BQB_TYPE_DOUBLE)) {
		eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"PosX\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcx->m_iType);
		return false;
	}
	if ((btfcy = btf.GetColumn("PosY")) == NULL) {
		eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"PosY\" does not exist in input trajectory.\n");
		return false;
	}
	if ((btfcy->m_iType != BQB_TYPE_FLOAT) && (btfcy->m_iType != BQB_TYPE_DOUBLE)) {
		eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"PosY\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcy->m_iType);
		return false;
	}
	if ((btfcz = btf.GetColumn("PosZ")) == NULL) {
		eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"PosZ\" does not exist in input trajectory.\n");
		return false;
	}
	if ((btfcz->m_iType != BQB_TYPE_FLOAT) && (btfcz->m_iType != BQB_TYPE_DOUBLE)) {
		eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"PosZ\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcz->m_iType);
		return false;
	}

	if (m_vaCoords.GetSize() < (long)m_iGesAtomCount)
		m_vaCoords.SetSize(m_iGesAtomCount);
	if (g_bKeepOriginalCoords)
		if (m_vaCoords_Original.GetSize() < (long)m_iGesAtomCount)
			m_vaCoords_Original.SetSize(m_iGesAtomCount);

	for (z=0;z<(int)m_iGesAtomCount;z++) {
		m_vaCoords[z] = CxDVector3(btfcx->m_aReal[z],btfcy->m_aReal[z],btfcz->m_aReal[z]);
		if (g_bKeepOriginalCoords)
			m_vaCoords_Original[z] = m_vaCoords[z];
	}

	if (needinfo) {
		if ((btfc = btf.GetColumn("Label")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"Label\" does not exist in input trajectory.\n");
			return false;
		}
		if (btfc->m_iType != BQB_TYPE_STRING) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"Label\" is of type %d, which is not STRING.\n",btfc->m_iType);
			return false;
		}
		for (z=0;z<m_paLabels.GetSize();z++)
			delete[] (char*)m_paLabels[z];
		m_paLabels.RemoveAll();
		if (g_bDoubleBox)
			m_paLabels.SetSize(m_iGesAtomCount*g_iDoubleBoxFactor);
		else
			m_paLabels.SetSize(m_iGesAtomCount);
		for (z=0;z<(int)m_iGesAtomCount;z++) {
			try { p = new char[btfc->m_aString[z].length()+1]; } catch(...) { p = NULL; }
			if (p == NULL) NewException((double)(btfc->m_aString[z].length()+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			strcpy(p,(const char*)btfc->m_aString[z].c_str());
			m_paLabels[z] = p;
		}
	}

	if (g_bNPT) {

		if (btf.m_pCellMatrix != NULL) {

			if ((btf.m_pCellMatrix->GetAt(0,1) != 0) || (btf.m_pCellMatrix->GetAt(0,2) != 0) || (btf.m_pCellMatrix->GetAt(1,2) != 0) || (btf.m_pCellMatrix->GetAt(1,0) != 0) || (btf.m_pCellMatrix->GetAt(2,0) != 0) || (btf.m_pCellMatrix->GetAt(2,1) != 0)) {

				g_bFoundNonOrtho = true;

				g_mBoxFromOrtho = *btf.m_pCellMatrix;

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

				// Backward transformation matrix
				g_mBoxToOrtho = CxDMatrix3(g_mBoxFromOrtho);
				if (!g_mBoxToOrtho.Invert()) {
					eprintf("CTimeStep::ReadBQBTrajStep(): Error: Encountered singular cell matrix (cell volume is zero).\n");
					abort();
				}

				veca = CxDVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
				vecb = CxDVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
				vecc = CxDVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

				// Cell angles
				g_fBoxAngleA = acos(DotP(vecb,vecc) / vecb.GetLength() / vecc.GetLength()) * 180.0 / Pi;
				g_fBoxAngleB = acos(DotP(veca,vecc) / veca.GetLength() / vecc.GetLength()) * 180.0 / Pi;
				g_fBoxAngleC = acos(DotP(veca,vecb) / veca.GetLength() / vecb.GetLength()) * 180.0 / Pi;

				// Orthogonal bounding box
				g_fBoxX = fabs(g_mBoxFromOrtho(0,0)) + fabs(g_mBoxFromOrtho(1,0)) + fabs(g_mBoxFromOrtho(2,0));
				g_fBoxY = fabs(g_mBoxFromOrtho(0,1)) + fabs(g_mBoxFromOrtho(1,1)) + fabs(g_mBoxFromOrtho(2,1));
				g_fBoxZ = fabs(g_mBoxFromOrtho(0,2)) + fabs(g_mBoxFromOrtho(1,2)) + fabs(g_mBoxFromOrtho(2,2));

				// Minimal diameters
				g_fBoxMinDiamA = fabs(DotP(veca,Normalize(CrossP(vecb,vecc))));
				g_fBoxMinDiamB = fabs(DotP(vecb,Normalize(CrossP(veca,vecc))));
				g_fBoxMinDiamC = fabs(DotP(vecc,Normalize(CrossP(veca,vecb))));
				g_fBoxMinDiam = MIN3(g_fBoxMinDiamA,g_fBoxMinDiamB,g_fBoxMinDiamC);

			} else {

				g_fBoxX = btf.m_pCellMatrix->GetAt(0,0);
				g_fBoxY = btf.m_pCellMatrix->GetAt(1,1);
				g_fBoxZ = btf.m_pCellMatrix->GetAt(2,2);
				g_fBoxMinDiam = MIN3(g_fBoxX,g_fBoxY,g_fBoxZ);
			}
		}
	}

/*	if (g_bUseVelocities && ((btfcx = btf.GetColumn("VelX")) != NULL) && ((btfcy = btf.GetColumn("VelY")) != NULL) && ((btfcz = btf.GetColumn("VelZ")) != NULL)) {
		if ((btfcx->m_iType != BQB_TYPE_FLOAT) && (btfcx->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"VelX\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcx->m_iType);
			return false;
		}
		if ((btfcy->m_iType != BQB_TYPE_FLOAT) && (btfcy->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"VelY\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcy->m_iType);
			return false;
		}
		if ((btfcz->m_iType != BQB_TYPE_FLOAT) && (btfcz->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"VelZ\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcz->m_iType);
			return false;
		}
		if (m_vaVelocities.GetSize() < (int)m_iGesAtomCount)
			m_vaVelocities.SetSize(m_iGesAtomCount);
		for (z=0;z<(int)m_iGesAtomCount;z++)
			m_vaVelocities[z] = CxDVector3(btfcx->m_aReal[z],btfcy->m_aReal[z],btfcz->m_aReal[z]);
	}*/

	if (g_bVoroIntegrateCharge) {
		if ((btfc = btf.GetColumn("EChg")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EChg\" does not exist in input trajectory, but g_bVoroIntegrateCharge==true.\n");
			return false;
		}
		if ((btfc->m_iType != BQB_TYPE_FLOAT) && (btfc->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EChg\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfc->m_iType);
			return false;
		}
		if (m_faCharge.GetSize() < (int)m_iGesAtomCount)
			m_faCharge.SetSize(m_iGesAtomCount);
		for (z=0;z<(int)m_iGesAtomCount;z++)
			m_faCharge[z] = btfc->m_aReal[z];

		if (btf.GetColumn("ECCg") != NULL) { // Core charge included
			if ((btfc = btf.GetColumn("ECCg")) == NULL) {
				eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"ECCg\" does not exist in input trajectory.\n");
				return false;
			}
			if ((btfc->m_iType != BQB_TYPE_FLOAT) && (btfc->m_iType != BQB_TYPE_DOUBLE)) {
				eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"ECCg\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfc->m_iType);
				return false;
			}
			if (m_faCoreCharge.GetSize() < (int)m_iGesAtomCount)
				m_faCoreCharge.SetSize(m_iGesAtomCount);
			for (z=0;z<(int)m_iGesAtomCount;z++)
				m_faCoreCharge[z] = btfc->m_aReal[z];
		}
	}

	if (g_bVoroIntegrateDipoleMoment) {
		if ((btfcx = btf.GetColumn("EDipX")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EDipX\" does not exist in input trajectory, but g_bVoroIntegrateDipoleMoment==true.\n");
			return false;
		}
		if ((btfcx->m_iType != BQB_TYPE_FLOAT) && (btfcx->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EDipX\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcx->m_iType);
			return false;
		}
		if ((btfcy = btf.GetColumn("EDipY")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EDipY\" does not exist in input trajectory, but g_bVoroIntegrateDipoleMoment==true.\n");
			return false;
		}
		if ((btfcy->m_iType != BQB_TYPE_FLOAT) && (btfcy->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EDipY\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcy->m_iType);
			return false;
		}
		if ((btfcz = btf.GetColumn("EDipZ")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EDipZ\" does not exist in input trajectory, but g_bVoroIntegrateDipoleMoment==true.\n");
			return false;
		}
		if ((btfcz->m_iType != BQB_TYPE_FLOAT) && (btfcz->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EDipZ\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcz->m_iType);
			return false;
		}
		if (m_dipoleMoments.GetSize() < (int)m_iGesAtomCount)
			m_dipoleMoments.SetSize(m_iGesAtomCount);
		for (z=0;z<(int)m_iGesAtomCount;z++)
			m_dipoleMoments[z] = CxDVector3(btfcx->m_aReal[z],btfcy->m_aReal[z],btfcz->m_aReal[z]);

		if (btf.GetColumn("EChZX") != NULL) { // Charge centers are included
			if ((btfcx = btf.GetColumn("EChZX")) == NULL) {
				eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EChZX\" does not exist in input trajectory.\n");
				return false;
			}
			if ((btfcx->m_iType != BQB_TYPE_FLOAT) && (btfcx->m_iType != BQB_TYPE_DOUBLE)) {
				eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EChZX\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcx->m_iType);
				return false;
			}
			if ((btfcy = btf.GetColumn("EChZY")) == NULL) {
				eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EChZY\" does not exist in input trajectory.\n");
				return false;
			}
			if ((btfcy->m_iType != BQB_TYPE_FLOAT) && (btfcy->m_iType != BQB_TYPE_DOUBLE)) {
				eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EChZY\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcy->m_iType);
				return false;
			}
			if ((btfcz = btf.GetColumn("EChZZ")) == NULL) {
				eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EChZZ\" does not exist in input trajectory.\n");
				return false;
			}
			if ((btfcz->m_iType != BQB_TYPE_FLOAT) && (btfcz->m_iType != BQB_TYPE_DOUBLE)) {
				eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EChZZ\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcz->m_iType);
				return false;
			}
			if (m_chargeCenters.GetSize() < (int)m_iGesAtomCount)
				m_chargeCenters.SetSize(m_iGesAtomCount);
			for (z=0;z<(int)m_iGesAtomCount;z++)
				m_chargeCenters[z] = CxDVector3(btfcx->m_aReal[z],btfcy->m_aReal[z],btfcz->m_aReal[z]);
		}
	}

	if (g_bVoroIntegrateQuadrupoleMoment) {
		if ((btfcxx = btf.GetColumn("EQXX")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQXX\" does not exist in input trajectory, but g_bVoroIntegrateQuadrupoleMoment==true.\n");
			return false;
		}
		if ((btfcxx->m_iType != BQB_TYPE_FLOAT) && (btfcxx->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQXX\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcxx->m_iType);
			return false;
		}
		if ((btfcxy = btf.GetColumn("EQXY")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQXY\" does not exist in input trajectory, but g_bVoroIntegrateQuadrupoleMoment==true.\n");
			return false;
		}
		if ((btfcxy->m_iType != BQB_TYPE_FLOAT) && (btfcxy->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQXY\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcxy->m_iType);
			return false;
		}
		if ((btfcxz = btf.GetColumn("EQXZ")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQXZ\" does not exist in input trajectory, but g_bVoroIntegrateQuadrupoleMoment==true.\n");
			return false;
		}
		if ((btfcxz->m_iType != BQB_TYPE_FLOAT) && (btfcxz->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQXZ\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcxz->m_iType);
			return false;
		}
		if ((btfcyx = btf.GetColumn("EQYX")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQYX\" does not exist in input trajectory, but g_bVoroIntegrateQuadrupoleMoment==true.\n");
			return false;
		}
		if ((btfcyx->m_iType != BQB_TYPE_FLOAT) && (btfcyx->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQYX\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcyx->m_iType);
			return false;
		}
		if ((btfcyy = btf.GetColumn("EQYY")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQYY\" does not exist in input trajectory, but g_bVoroIntegrateQuadrupoleMoment==true.\n");
			return false;
		}
		if ((btfcyy->m_iType != BQB_TYPE_FLOAT) && (btfcyy->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQYY\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcyy->m_iType);
			return false;
		}
		if ((btfcyz = btf.GetColumn("EQYZ")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQYZ\" does not exist in input trajectory, but g_bVoroIntegrateQuadrupoleMoment==true.\n");
			return false;
		}
		if ((btfcyz->m_iType != BQB_TYPE_FLOAT) && (btfcyz->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQYZ\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcyz->m_iType);
			return false;
		}
		if ((btfczx = btf.GetColumn("EQZX")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQZX\" does not exist in input trajectory, but g_bVoroIntegrateQuadrupoleMoment==true.\n");
			return false;
		}
		if ((btfczx->m_iType != BQB_TYPE_FLOAT) && (btfczx->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQZX\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfczx->m_iType);
			return false;
		}
		if ((btfczy = btf.GetColumn("EQZY")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQZY\" does not exist in input trajectory, but g_bVoroIntegrateQuadrupoleMoment==true.\n");
			return false;
		}
		if ((btfczy->m_iType != BQB_TYPE_FLOAT) && (btfczy->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQZY\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfczy->m_iType);
			return false;
		}
		if ((btfczz = btf.GetColumn("EQZZ")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQZZ\" does not exist in input trajectory, but g_bVoroIntegrateQuadrupoleMoment==true.\n");
			return false;
		}
		if ((btfczz->m_iType != BQB_TYPE_FLOAT) && (btfczz->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"EQZZ\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfczz->m_iType);
			return false;
		}

		if (m_maQuadTensor.size() < m_iGesAtomCount)
			m_maQuadTensor.resize(m_iGesAtomCount);
		for (z=0;z<(int)m_iGesAtomCount;z++)
			m_maQuadTensor[z] = CxDMatrix3(btfcxx->m_aReal[z],btfcxy->m_aReal[z],btfcxz->m_aReal[z],btfcyx->m_aReal[z],btfcyy->m_aReal[z],btfcyz->m_aReal[z],btfczx->m_aReal[z],btfczy->m_aReal[z],btfczz->m_aReal[z]);
	}

	if (g_bVoroIntegrateTotalCurrent) {
		if ((btfcx = btf.GetColumn("ECurX")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"ECurX\" does not exist in input trajectory, but g_bVoroIntegrateDipoleMoment==true.\n");
			return false;
		}
		if ((btfcx->m_iType != BQB_TYPE_FLOAT) && (btfcx->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"ECurX\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcx->m_iType);
			return false;
		}
		if ((btfcy = btf.GetColumn("ECurY")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"ECurY\" does not exist in input trajectory, but g_bVoroIntegrateDipoleMoment==true.\n");
			return false;
		}
		if ((btfcy->m_iType != BQB_TYPE_FLOAT) && (btfcy->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"ECurY\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcy->m_iType);
			return false;
		}
		if ((btfcz = btf.GetColumn("ECurZ")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"ECurZ\" does not exist in input trajectory, but g_bVoroIntegrateDipoleMoment==true.\n");
			return false;
		}
		if ((btfcz->m_iType != BQB_TYPE_FLOAT) && (btfcz->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"ECurZ\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcz->m_iType);
			return false;
		}
		if (m_totalCurrents.GetSize() < (int)m_iGesAtomCount)
			m_totalCurrents.SetSize(m_iGesAtomCount);
		for (z=0;z<(int)m_iGesAtomCount;z++)
			m_totalCurrents[z] = CxDVector3(btfcx->m_aReal[z],btfcy->m_aReal[z],btfcz->m_aReal[z]);
	}

	if (g_bVoroIntegrateMagneticMoment) {
		if ((btfcx = btf.GetColumn("MDipX")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"MDipX\" does not exist in input trajectory, but g_bVoroIntegrateDipoleMoment==true.\n");
			return false;
		}
		if ((btfcx->m_iType != BQB_TYPE_FLOAT) && (btfcx->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"MDipX\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcx->m_iType);
			return false;
		}
		if ((btfcy = btf.GetColumn("MDipY")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"MDipY\" does not exist in input trajectory, but g_bVoroIntegrateDipoleMoment==true.\n");
			return false;
		}
		if ((btfcy->m_iType != BQB_TYPE_FLOAT) && (btfcy->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"MDipY\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcy->m_iType);
			return false;
		}
		if ((btfcz = btf.GetColumn("MDipZ")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"MDipZ\" does not exist in input trajectory, but g_bVoroIntegrateDipoleMoment==true.\n");
			return false;
		}
		if ((btfcz->m_iType != BQB_TYPE_FLOAT) && (btfcz->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"MDipZ\" is of type %d, which is neither FLOAT nor DOUBLE.\n",btfcz->m_iType);
			return false;
		}
		if (m_magneticDipoleMoments.GetSize() < (int)m_iGesAtomCount)
			m_magneticDipoleMoments.SetSize(m_iGesAtomCount);
		for (z=0;z<(int)m_iGesAtomCount;z++)
			m_magneticDipoleMoments[z] = CxDVector3(btfcx->m_aReal[z],btfcy->m_aReal[z],btfcz->m_aReal[z]);
	}

	if (g_bUseVelocities) {
		if ((velx = btf.GetColumn("VelX")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"VelX\" does not exist in input trajectory, but g_bUseVelocities==true.\n");
			return false;
		}
		if ((velx->m_iType != BQB_TYPE_FLOAT) && (velx->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"VelX\" is of type %d, which is neither FLOAT nor DOUBLE.\n",velx->m_iType);
			return false;
		}
		if ((vely = btf.GetColumn("VelY")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"VelY\" does not exist in input trajectory, but g_bUseVelocities==true.\n");
			return false;
		}
		if ((vely->m_iType != BQB_TYPE_FLOAT) && (vely->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"VelY\" is of type %d, which is neither FLOAT nor DOUBLE.\n",vely->m_iType);
			return false;
		}
		if ((velz = btf.GetColumn("VelZ")) == NULL) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"VelZ\" does not exist in input trajectory, but g_bUseVelocities==true.\n");
			return false;
		}
		if ((velz->m_iType != BQB_TYPE_FLOAT) && (velz->m_iType != BQB_TYPE_DOUBLE)) {
			eprintf("\nCTimeStep::ReadBQBTrajStep(): Error: Column \"VelZ\" is of type %d, which is neither FLOAT nor DOUBLE.\n",velz->m_iType);
			return false;
		}
		if (m_vaVelocities.GetSize() < (int)m_iGesAtomCount)
			m_vaVelocities.SetSize(m_iGesAtomCount);
		for (z=0;z<(int)m_iGesAtomCount;z++)
			m_vaVelocities[z] = CxDVector3(velx->m_aReal[z],vely->m_aReal[z],velz->m_aReal[z]);
	}

	if (g_bVerbose)
		mprintf("<<< CTimeStep::ReadBQBTrajStep() <<<\n");
	
	return true;
}


bool CTimeStep::ReadBQB(CBQBFile *source, bool needinfo) {

	int ty, ver;


	if (g_bVerbose)
		mprintf(">>> CTimeStep::ReadBQB() >>>\n");

	if (source == NULL) {
		eprintf("CTimeStep::ReadBQB(): Error: source == NULL.\n");
		abort();
	}

_again:
	if (!source->ReadFrame()) {
		if (g_bVerbose)
			mprintf("CTimeStep::ReadBQB(): Error while reading BQB frame.\n");
		return false;
	}

	ty = source->GetFrameType();
	ver = source->GetFrameTypeVersion();

	if (g_bVerbose)
		mprintf("CTimeStep::ReadBQB(): Frame type is %d v%d.\n",ty,ver);

	if (((ty == BQB_FRAMETYPE_COMPCUBESTART) || 
		(ty == BQB_FRAMETYPE_COMPCUBE) ||
		(ty == BQB_FRAMETYPE_COMPCUBEKEY)) 
		&& (ver < 3)) {

		if (g_bVerbose)
			mprintf("  --> Processing compressed cube frame.\n");

		if (!g_pBQBFile->GetEngine()->DecompressCubeStep(
				source->GetFramePayload(),
				ver,
				&source->m_oParmSetCube,
				NULL, // ahistused
				NULL, // chistused
				false // skipvol
			))
			return false;

		if (!ReadBQBCubeStep(needinfo))
			return false;

	} else if (((ty == BQB_FRAMETYPE_COMPTRAJSTART) || 
		(ty == BQB_FRAMETYPE_COMPTRAJ) ||
		(ty == BQB_FRAMETYPE_COMPTRAJKEY)) 
		&& (ver < 3)) {

		if (g_bVerbose)
			mprintf("  --> Processing compressed trajectory frame.\n");

		if (!g_pBQBFile->GetEngine()->DecompressAtomStep(
				source->GetFramePayload(),
				ver,
				&source->m_oParmSetPos,
				NULL, // histused
				false, // second
				false // check
			))
			return false;

		if (!FetchBQBCompTrajFrame(needinfo))
			return false;

	} else if ((ty == 4) && (ver == 0)) {

		if (g_bVerbose)
			mprintf("  --> Processing uncompressed trajectory frame.\n");

		if (!ReadBQBTrajStep(source->GetFramePayload(),needinfo))
			return false;

	} else if ((ty == 2) || (ty == 3)) {

		if (g_bVerbose)
			mprintf("CTimeStep::ReadBQB(): Skipping index frame.\n");
		goto _again;

	} else {

		eprintf("CTimeStep::ReadBQB(): Error: Frame type %d v%d not yet implemented.\n",
			source->GetFrameType(),source->GetFrameTypeVersion());
		return false;
	}

	if (g_bVerbose)
		mprintf("<<< CTimeStep::ReadBQB() <<<\n");

	return true;
}


bool CTimeStep::FetchBQBCompTrajFrame(bool needinfo) {

	CBQBAtomSet *as;
	CBQBAtom *at;
	int z;
	char *p;


	if (m_pComment == NULL) {
		try { m_pComment = new char[256]; } catch(...) { m_pComment = NULL; }
		if (m_pComment == NULL) NewException((double)256*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}
	m_pComment[0] = 0;

	as = g_pBQBFile->GetEngine()->GetOutputAtomFrame(0);

	strcpy(m_pComment,as->m_sComment);

	m_iGesAtomCount = (unsigned long)as->m_oaAtoms.size();
	
	if (m_iGesAtomCount == 0) {
		eprintf("\nCTimeStep::FetchBQBCompTrajFrame(): Error: Atom count = 0.");
		return false;
	}

	if (needinfo) {
		for (z=0;z<m_paLabels.GetSize();z++)
			delete[] (char*)m_paLabels[z];
		m_paLabels.RemoveAll();
		if (g_bDoubleBox)
			m_paLabels.SetSize(m_iGesAtomCount*g_iDoubleBoxFactor);
		else
			m_paLabels.SetSize(m_iGesAtomCount);
	}
	if (m_vaCoords.GetSize() < (long)m_iGesAtomCount)
		m_vaCoords.SetSize(m_iGesAtomCount);

	if (g_bKeepOriginalCoords)
		if (m_vaCoords_Original.GetSize() < (long)m_iGesAtomCount)
			m_vaCoords_Original.SetSize(m_iGesAtomCount);

	for (z=0;z<(int)m_iGesAtomCount;z++) {

		at = as->m_oaAtoms[z];

		if (needinfo) {
			try { p = new char[at->m_sLabel.length()+1]; } catch(...) { p = NULL; }
			if (p == NULL) NewException((double)(at->m_sLabel.length()+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			strcpy(p,at->m_sLabel.c_str());
			m_paLabels[z] = p;
		}

		m_vaCoords[z][0] = at->m_fCoord[0] * 100.0;
		m_vaCoords[z][1] = at->m_fCoord[1] * 100.0;
		m_vaCoords[z][2] = at->m_fCoord[2] * 100.0;

		if (g_bKeepOriginalCoords)
			m_vaCoords_Original[z] = m_vaCoords[z];
	}

	return true;
}


bool CTimeStep::SkipBQB() {

	int ty, ver;


	if (g_bVerbose)
		mprintf(">>> CTimeStep::SkipBQB() >>>\n");

	if (g_pBQBFile == NULL) {
		eprintf("CTimeStep::SkipBQB(): Error: g_pBQBFile == NULL.\n");
		abort();
	}

_again:
	if (!g_pBQBFile->ReadFrame()) {
		if (g_bVerbose)
			mprintf("CTimeStep::SkipBQB(): Error while reading BQB frame.\n");
		return false;
	}

	ty = g_pBQBFile->GetFrameType();
	ver = g_pBQBFile->GetFrameTypeVersion();

	if (g_bVerbose)
		mprintf("CTimeStep::SkipBQB(): Frame type is %d v%d.\n",ty,ver);

	if (((ty == BQB_FRAMETYPE_COMPCUBESTART) || 
		(ty == BQB_FRAMETYPE_COMPCUBE) ||
		(ty == BQB_FRAMETYPE_COMPCUBEKEY)) 
		&& (ver < 3)) {

		if (g_bVerbose)
			mprintf("  --> Skipping compressed cube frame.\n");

		if (!g_pBQBFile->GetEngine()->DecompressCubeStep(
				g_pBQBFile->GetFramePayload(),
				ver,
				&g_pBQBFile->m_oParmSetCube,
				NULL, // ahistused
				NULL, // chistused
				false // skipvol
			))
			return false;

	} else if (((ty == BQB_FRAMETYPE_COMPTRAJSTART) || 
		(ty == BQB_FRAMETYPE_COMPTRAJ) ||
		(ty == BQB_FRAMETYPE_COMPTRAJKEY)) 
		&& (ver < 3)) {

		if (g_bVerbose)
			mprintf("  --> Skipping compressed trajectory frame.\n");

		if (!g_pBQBFile->GetEngine()->DecompressAtomStep(
				g_pBQBFile->GetFramePayload(),
				ver,
				&g_pBQBFile->m_oParmSetPos,
				NULL, // histused
				false, // second
				false // check
			))
			return false;

	} else if ((ty == 4) && (ver == 0)) {

		if (g_bVerbose)
			mprintf("  --> Skipping uncompressed trajectory frame.\n");

	} else if ((ty == 2) || (ty == 3)) {

		if (g_bVerbose)
			mprintf("  Skipping index frame.\n");

		goto _again;

	} else {

		eprintf("CTimeStep::SkipBQB(): Error: Frame type %d v%d not yet implemented.\n",ty,ver);
		return false;
	}

	if (g_bVerbose)
		mprintf("<<< CTimeStep::SkipBQB() <<<\n");

	return true;
}


bool CTimeStep::SkipBQB(CBQBFile *source) {

	int ty, ver;


	if (g_bVerbose)
		mprintf(">>> CTimeStep::SkipBQB() >>>\n");

	if (source == NULL) {
		eprintf("CTimeStep::SkipBQB(): Error: source == NULL.\n");
		abort();
	}

_again:
	if (!source->ReadFrame()) {
		if (g_bVerbose)
			mprintf("CTimeStep::SkipBQB(): Error while reading BQB frame.\n");
		return false;
	}

	ty = source->GetFrameType();
	ver = source->GetFrameTypeVersion();

	if (g_bVerbose)
		mprintf("CTimeStep::SkipBQB(): Frame type is %d v%d.\n",ty,ver);

	if (((ty == BQB_FRAMETYPE_COMPCUBESTART) || 
		(ty == BQB_FRAMETYPE_COMPCUBE) ||
		(ty == BQB_FRAMETYPE_COMPCUBEKEY)) 
		&& (ver < 3)) {

		if (g_bVerbose)
			mprintf("  --> Skipping compressed cube frame.\n");

		if (!g_pBQBFile->GetEngine()->DecompressCubeStep(
				source->GetFramePayload(),
				ver,
				&source->m_oParmSetCube,
				NULL, // ahistused
				NULL, // chistused
				false // skipvol
			))
			return false;

	} else if (((ty == BQB_FRAMETYPE_COMPTRAJSTART) || 
		(ty == BQB_FRAMETYPE_COMPTRAJ) ||
		(ty == BQB_FRAMETYPE_COMPTRAJKEY)) 
		&& (ver < 3)) {

		if (g_bVerbose)
			mprintf("  --> Skipping compressed trajectory frame.\n");

		if (!g_pBQBFile->GetEngine()->DecompressAtomStep(
				source->GetFramePayload(),
				ver,
				&source->m_oParmSetPos,
				NULL, // histused
				false, // second
				false // check
			))
			return false;

	} else if ((ty == 4) && (ver == 0)) {

		if (g_bVerbose)
			mprintf("  --> Skipping uncompressed trajectory frame.\n");

	} else if ((ty == 2) || (ty == 3)) {

		if (g_bVerbose)
			mprintf("  Skipping index frame.\n");

		goto _again;

	} else {

		eprintf("CTimeStep::SkipBQB(): Error: Frame type %d v%d not yet implemented.\n",ty,ver);
		return false;
	}

	if (g_bVerbose)
		mprintf("<<< CTimeStep::SkipBQB() <<<\n");

	return true;
}


bool CTimeStep::ReadCube(CxMemFile *file)
{
	int rx, ry, rz, ac, z, k;
	char buf[256], *p;
	const char *q;
	double fx, fy, fz, tf;
	
	file->Seek(0);
	(void)!file->fgets(buf, 256);
	(void)!file->fgets(buf, 256);
	
	if (file->Eof())
	{
		eprintf("\nEnd of file reached while reading cube file.\n");
		return false;
	}
	
	(void)!file->fgets(buf, 256);
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;
	
	if (file->Eof())
	{
		eprintf("Unexpected end of file reached while reading cube file.\n");
		return false;
	}
	
	ac = atoi(buf);
	m_iGesAtomCount = ac;
	
	if (ac <= 0)
	{
		eprintf("Encountered atom count of %d while reading cube file.\n",ac);
		return false;
	}
	
	(void)!file->fgets(buf, 256);
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;
	rx = atoi(buf);
	
/*	p++;
	while (*p == ' ')
		p++;
	q = p;
	while ((*p != ' ') && (*p != 0))
		p++;
	*p = 0;
	g_fCubeXStep = atof(q) / g_iCubeXStride;*/

	for (k = 0; k < 3; k++) {
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		*p = 0;
		
		tf = atof(q);
		if (g_bFlipCellVectorX && (k == 0) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorY && (k == 1) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorZ && (k == 2) && (tf != 0))
			tf *= -1.0;
		if (m_bFirst) {
			g_fCubeXVector[k] = tf;
		} else {
			g_fCubeXVector[k] = g_mCubeCell(0, k) / (rx * g_iCubeXStride - g_iCubeXMismatch) / LEN_AU2PM;
			if (fabs(g_fCubeXVector[k] - tf / g_iCubeXStride) > 0.1) {
				eprintf("CTimeStep::ReadCube(): Error: Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeXVector[k] * LEN_AU2PM, tf / g_iCubeXStride * LEN_AU2PM);
				return false;
			}
		}
	}
	g_fCubeXStep = g_fCubeXVector[0];

	
	(void)!file->fgets(buf, 256);
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;
	ry = atoi(buf);
	
/*	for (z = 0; z < 2; z++) {
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		*p = 0;
	}
	g_fCubeYStep = atof(q) / g_iCubeYStride;*/

	for (k = 0; k < 3; k++) {
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		*p = 0;
		
		tf = atof(q);
		if (g_bFlipCellVectorX && (k == 0) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorY && (k == 1) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorZ && (k == 2) && (tf != 0))
			tf *= -1.0;
		if (m_bFirst) {
			g_fCubeYVector[k] = tf;
		} else {
			g_fCubeYVector[k] = g_mCubeCell(1, k) / (ry * g_iCubeYStride - g_iCubeYMismatch) / LEN_AU2PM;
			if (fabs(g_fCubeYVector[k] - tf / g_iCubeYStride) > 0.1) {
				eprintf("CTimeStep::ReadCube(): Error: Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeYVector[k] * LEN_AU2PM, tf / g_iCubeYStride * LEN_AU2PM);
				return false;
			}
		}
	}
	g_fCubeYStep = g_fCubeYVector[1];

	
	(void)!file->fgets(buf, 256);
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;
	rz = atoi(buf);
	
/*	for (z = 0; z < 3; z++) {
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		*p = 0;
	}
	g_fCubeZStep = atof(q) / g_iCubeZStride;*/

	for (k = 0; k < 3; k++) {
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		*p = 0;
		
		tf = atof(q);
		if (g_bFlipCellVectorX && (k == 0) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorY && (k == 1) && (tf != 0))
			tf *= -1.0;
		if (g_bFlipCellVectorZ && (k == 2) && (tf != 0))
			tf *= -1.0;
		if (m_bFirst) {
			g_fCubeZVector[k] = tf;
		} else {
			g_fCubeZVector[k] = g_mCubeCell(2, k) / (rz * g_iCubeZStride - g_iCubeZMismatch) / LEN_AU2PM;
			if (fabs(g_fCubeZVector[k] - tf / g_iCubeZStride) > 0.1) {
				eprintf("CTimeStep::ReadCube(): Error: Grid point distance in the cube file changed from %.3f pm to %.3f pm!\n", g_fCubeZVector[k] * LEN_AU2PM, tf / g_iCubeZStride * LEN_AU2PM);
				return false;
			}
		}
	}
	g_fCubeZStep = g_fCubeZVector[2];

	
	m_vaCoords.RemoveAll_KeepSize();
	CxIntArray atomNumbers;
	atomNumbers.SetSize(ac);
	for (z=0;z<ac;z++)
	{
		(void)!file->fgets(buf, 256);
		buf[strlen(buf)-1] = 0;
		
		if (file->Eof())
		{
			eprintf("Unexpected end of cube file while reading coordinates (%d/%d).\n",z+1,ac);
			return false;
		}
		
		p = buf;
		while (*p == ' ')
			p++;
		while (*p != ' ')
			p++;
		*p = 0;
		
		atomNumbers[z] = atoi(buf);
		
		p++; // Column is always zero
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		if (*p == 0)
		{
			eprintf("Incomplete line while reading coordinates from cube file (A).\n");
			return false;
		}
		*p = 0;
		
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		if (*p == 0)
		{
			eprintf("Incomplete line while reading coordinates from cube file (B).\n");
			return false;
		}
		*p = 0;
		fx = atof(q);
		
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		if (*p == 0)
		{
			eprintf("Incomplete line while reading coordinates from cube file (C).\n");
			return false;
		}
		fy = atof(q);
		
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		*p = 0;
		fz = atof(q);
		
		m_vaCoords.Add(CxDVector3(fx*LEN_AU2PM,fy*LEN_AU2PM,fz*LEN_AU2PM));
	}
	
	if (file->Eof())
	{
		eprintf("Unexpected end of file while reading cube file.\n");
		return false;
	}
	
	if (m_pVolumetricData == NULL)
	{
		m_pVolumetricData = new C3DF<VORI_FLOAT>();
	}
	
	if (m_pVolumetricDataTimeDev == NULL && g_bCubeTimeDev) {
		m_pVolumetricDataTimeDev = new C3DF<VORI_FLOAT>();
	}
	if (m_pCurrentDensity == NULL && g_bCubeTimeDev) {
		m_pCurrentDensity = new CxDoubleArray();
	}
	
	static bool second = false;
	if (m_pVolumetricData->m_pBin == NULL)
	{
		m_pVolumetricData->m_iRes[0] = rx * g_iCubeXStride - g_iCubeXMismatch;
		m_pVolumetricData->m_iRes[1] = ry * g_iCubeYStride - g_iCubeYMismatch;
		m_pVolumetricData->m_iRes[2] = rz * g_iCubeZStride - g_iCubeZMismatch;
		
		m_pVolumetricData->m_fMinVal[0] = 0;
		m_pVolumetricData->m_fMaxVal[0] = g_fCubeXStep * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
		m_pVolumetricData->m_fMinVal[1] = 0;
		m_pVolumetricData->m_fMaxVal[1] = g_fCubeYStep * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
		m_pVolumetricData->m_fMinVal[2] = 0;
		m_pVolumetricData->m_fMaxVal[2] = g_fCubeZStep * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
		
		m_pVolumetricData->Create();
		
		if (m_pVolumetricDataTimeDev != NULL) {
			m_pVolumetricDataTimeDev->m_iRes[0] = rx * g_iCubeXStride - g_iCubeXMismatch;
			m_pVolumetricDataTimeDev->m_iRes[1] = ry * g_iCubeYStride - g_iCubeYMismatch;
			m_pVolumetricDataTimeDev->m_iRes[2] = rz * g_iCubeZStride - g_iCubeZMismatch;
			m_pVolumetricDataTimeDev->m_fMinVal[0] = 0;
			m_pVolumetricDataTimeDev->m_fMaxVal[0] = g_fCubeXStep * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
			m_pVolumetricDataTimeDev->m_fMinVal[1] = 0;
			m_pVolumetricDataTimeDev->m_fMaxVal[1] = g_fCubeYStep * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
			m_pVolumetricDataTimeDev->m_fMinVal[2] = 0;
			m_pVolumetricDataTimeDev->m_fMaxVal[2] = g_fCubeZStep * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
			
			m_pVolumetricDataTimeDev->Create();
		}
		if (m_pCurrentDensity != NULL) {
			m_pCurrentDensity->SetSize(3 * m_pVolumetricDataTimeDev->m_iRes[0] * m_pVolumetricDataTimeDev->m_iRes[1] * m_pVolumetricDataTimeDev->m_iRes[2]);
		}
		
		second = true;
	} else
	{
		if ((rx * g_iCubeXStride - g_iCubeXMismatch != m_pVolumetricData->m_iRes[0]) || (ry * g_iCubeYStride - g_iCubeYMismatch!= m_pVolumetricData->m_iRes[1]) || (rz * g_iCubeZStride - g_iCubeZMismatch != m_pVolumetricData->m_iRes[2]))
		{
			if (second) {
				second = false;
				
				m_pVolumetricData->m_iRes[0] = rx * g_iCubeXStride - g_iCubeXMismatch;
				m_pVolumetricData->m_iRes[1] = ry * g_iCubeYStride - g_iCubeYMismatch;
				m_pVolumetricData->m_iRes[2] = rz * g_iCubeZStride - g_iCubeZMismatch;
				
				m_pVolumetricData->m_fMinVal[0] = 0;
				m_pVolumetricData->m_fMaxVal[0] = g_fCubeXStep * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
				m_pVolumetricData->m_fMinVal[1] = 0;
				m_pVolumetricData->m_fMaxVal[1] = g_fCubeYStep * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
				m_pVolumetricData->m_fMinVal[2] = 0;
				m_pVolumetricData->m_fMaxVal[2] = g_fCubeZStep * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
				
				m_pVolumetricData->Create();
				if (m_pVolumetricDataTimeDev != NULL) {
					m_pVolumetricDataTimeDev->m_iRes[0] = rx * g_iCubeXStride - g_iCubeXMismatch;
					m_pVolumetricDataTimeDev->m_iRes[1] = ry * g_iCubeYStride - g_iCubeYMismatch;
					m_pVolumetricDataTimeDev->m_iRes[2] = rz * g_iCubeZStride - g_iCubeZMismatch;
					m_pVolumetricDataTimeDev->m_fMinVal[0] = 0;
					m_pVolumetricDataTimeDev->m_fMaxVal[0] = g_fCubeXStep * LEN_AU2PM * (rx * g_iCubeXStride - g_iCubeXMismatch);
					m_pVolumetricDataTimeDev->m_fMinVal[1] = 0;
					m_pVolumetricDataTimeDev->m_fMaxVal[1] = g_fCubeYStep * LEN_AU2PM * (ry * g_iCubeYStride - g_iCubeYMismatch);
					m_pVolumetricDataTimeDev->m_fMinVal[2] = 0;
					m_pVolumetricDataTimeDev->m_fMaxVal[2] = g_fCubeZStep * LEN_AU2PM * (rz * g_iCubeZStride - g_iCubeZMismatch);
					
					m_pVolumetricDataTimeDev->Create();
				}
				if (m_pCurrentDensity != NULL) {
					m_pCurrentDensity->SetSize(3 * m_pVolumetricDataTimeDev->m_iRes[0] * m_pVolumetricDataTimeDev->m_iRes[1] * m_pVolumetricDataTimeDev->m_iRes[2]);
				}
			} else {
				eprintf("\nCube file dimension mismatch (%d-%d, %d-%d, %d-%d).\n",rx * g_iCubeXStride - g_iCubeXMismatch,m_pVolumetricData->m_iRes[0],ry * g_iCubeYStride - g_iCubeYMismatch,m_pVolumetricData->m_iRes[1],rz * g_iCubeZStride - g_iCubeZMismatch,m_pVolumetricData->m_iRes[2]);
				return false;
			}
		}
		for (z=0;z<m_pVolumetricData->m_iRes[0]*m_pVolumetricData->m_iRes[1]*m_pVolumetricData->m_iRes[2];z++)
			m_pVolumetricData->m_pBin[z] = 0;
	}
	
	m_pVolumetricData->ReadCubeData(file, false);
	
	
/* 	FILE *cubeFile = fopen("test.cube", "w");
 	if (cubeFile == NULL) {
 		printf("Could not open output file!\n");
 		return 1;
 	}
 	
 	fprintf(cubeFile, "\n\n");
 	fprintf(cubeFile, "%5lu %12.6f %12.6f %12.6f\n", m_iGesAtomCount, 0.0f, 0.0f, 0.0f);
 	fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", m_pVolumetricData->m_iRes[0], g_fCubeXStep, 0.0f, 0.0f);
 	fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", m_pVolumetricData->m_iRes[1], 0.0f, g_fCubeYStep, 0.0f);
 	fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f\n", m_pVolumetricData->m_iRes[2], 0.0f, 0.0f, g_fCubeZStep);
 	for (int i = 0; i < (int)m_iGesAtomCount; i++) {
 		fprintf(cubeFile, "%5d %12.6f %12.6f %12.6f %12.6f\n", 1, 0.0f, m_vaCoords[i][0] / 0.529177249f / 100.0f, m_vaCoords[i][1] / 0.529177249f / 100.0f, m_vaCoords[i][2] / 0.529177249f / 100.0f);
 	}
 	for (int i = 0; i < m_pVolumetricData->m_iRes[0]; i++) {
 		for (int j = 0; j < m_pVolumetricData->m_iRes[1]; j++) {
 			for (int k = 0; k < m_pVolumetricData->m_iRes[2] / 6; k++) {
 				for (int l = 0; l < 6; l++) {
 					fprintf(cubeFile, "%13.5E", m_pVolumetricData->m_pBin[i + j * m_pVolumetricData->m_iRes[0] + (k * 6 + l) * m_pVolumetricData->m_iRes[0] * m_pVolumetricData->m_iRes[1]]);
 				}
 				fprintf(cubeFile, "\n");
 			}
 			if (m_pVolumetricData->m_iRes[2] % 6 != 0) {
 				for (int l = 0; l < m_pVolumetricData->m_iRes[2] % 6; l++) {
 					fprintf(cubeFile, "%13.5E", m_pVolumetricData->m_pBin[i + j * m_pVolumetricData->m_iRes[0] + ((m_pVolumetricData->m_iRes[2] / 6) * 6 + l) * m_pVolumetricData->m_iRes[0] * m_pVolumetricData->m_iRes[1]]);
 				}
 				fprintf(cubeFile, "\n");
 			}
 		}
 	}
 	
 	fclose(cubeFile);*/
	
	return true;
}


bool CTimeStep::SkipCube(FILE *a)
{
	int rx, ry, rz, ac, z, ix, iy, iz;
	char buf[256], *p, *q;

	(void)!fgets(buf,256,a);
	(void)!fgets(buf,256,a);

	(void)!fgets(buf,256,a);
	if (feof(a)) return false;
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;

	ac = atoi(buf);

	(void)!fgets(buf,256,a);
	if (feof(a)) { eprintf("\nCTimeStep::SkipCube(): Unexpected end of file while skipping in cube file.\n"); return false; }
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;
	rx = atoi(buf);

	(void)!fgets(buf,256,a);
	if (feof(a)) { eprintf("\nCTimeStep::SkipCube(): Unexpected end of file while skipping in cube file.\n"); return false; }
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;
	ry = atoi(buf);

	(void)!fgets(buf,256,a);
	if (feof(a)) { eprintf("\nCTimeStep::SkipCube(): Unexpected end of file while skipping in cube file.\n"); return false; }
	p = buf;
	while (*p == ' ')
		p++;
	while (*p != ' ')
		p++;
	*p = 0;
	rz = atoi(buf);

	for (z=0;z<ac;z++)
		(void)!fgets(buf,256,a);

	if (feof(a)) { eprintf("\nCTimeStep::SkipCube(): Unexpected end of file while skipping in cube file.\n"); return false; }

	ix = 0;
	iy = 0;
	iz = 0;

	while (!feof(a))
	{
_read:
		(void)!fgets(buf,256,a);
		if (feof(a)) { eprintf("\nCTimeStep::SkipCube(): Unexpected end of file while skipping in cube file.\n"); return false; }
		p = buf;

_next:
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != '\n') && (*p != 0))
			p++;
		if ((p-q) < 8)
			goto _read;

		*p = 0;

		iz++;
		if (iz >= rz)
		{
			iz = 0;
			iy++;
			if (iy >= ry)
			{
				iy = 0;
				ix++;
			}
		}

		if (ix == rx)
			return true;

		p++;
		goto _next;
	}
	
	return false;
}



bool CTimeStep::ReadPDB(FILE *a, bool needinfo, CxDVec3Array *v)
{
	int i;
	static char buf[256], obuf[256], buf2[64];
	char *p, *q, *r;
	double x, y, z;
	bool b;


	v->RemoveAll_KeepSize();
	for (i=0;i<m_paLabels.GetSize();i++)
		delete[] (char*)m_paLabels[i];
	m_paLabels.RemoveAll();
	if (g_bKeepOriginalCoords)
		m_vaCoords_Original.RemoveAll_KeepSize();
	while (true)
	{
		(void)fgets_bin(buf,256,a);
		if (feof(a)) {
			//return false;
			eprintf("\nWarning: Unexpected end of PDB file.\n");
			goto _ende;
		}
		if (strlen(buf) == 0)
			continue;
		buf[strlen(buf)-1] = 0;
		strcpy(obuf,buf);
		if (g_bNPT)
		{
			if (strstr(buf,"CRYST1") != 0) // Boxlaenge
			{
	//			mprintf(GREY,"\"%s\".\n",buf);
				p = &buf[6];
				while (*p == ' ')
					p++;
				q = p;
				while ((*q != ' ') && (*q != 0))
					q++;
				if (*q == 0)
				{
					eprintf("Error 5 reading PDB line: \"%s\".\n",obuf);
					return false;
				}
				*q = 0;
				g_fBoxX = atof(p)*100.0;
				p = q+1;

				while (*p == ' ')
					p++;
				q = p;
				while ((*q != ' ') && (*q != 0))
					q++;
				if (*q == 0)
				{
					eprintf("Error 6 reading PDB line: \"%s\".\n",obuf);
					return false;
				}
				*q = 0;
				g_fBoxY = atof(p)*100.0;
				p = q+1;

				while (*p == ' ')
					p++;
				q = p;
				while ((*q != ' ') && (*q != 0))
					q++;
				*q = 0;
				g_fBoxZ = atof(p)*100.0;
//				mprintf(GREY,"--> %f %f %f\n",g_fBoxX,g_fBoxY,g_fBoxZ);

				if (g_bDoubleBox) {
					g_fBoxX *= g_iDoubleBoxX;
					g_fBoxY *= g_iDoubleBoxY;
					g_fBoxZ *= g_iDoubleBoxZ;
				}

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

				g_mBoxToOrtho = CxDMatrix3(g_mBoxFromOrtho);
				if (!g_mBoxToOrtho.Invert()) {
					eprintf("CTimeStep::ReadPDB(): Error: Encountered singular cell matrix (cell volume is zero).\n");
					abort();
				}
			}
		}

		if (strstr(buf,"END") == buf)
			break;

		if ((strstr(buf,"ATOM")==0) && (strstr(buf,"HETATM")==0))
			continue;

		p = &buf[7];
		while (!(isalpha(*p) || (*p == '_')) && (*p != 0))
			p++;
		if (*p == 0)
		{
			eprintf("Error 4 reading PDB line: \"%s\".\n",obuf);
			return false;
		}
		q = p;
		while (isalpha(*q) || (*q == '_'))
			q++;
		if (needinfo)
		{
			*q = 0;
	//		if (strlen(p) > 7)
	//		{
	//			eprintf("\nCTimeStep::ReadPDB(): \"%s\" - maximum length for atom labels is 7 chars; truncating.\n",p);
	//			p[7] = 0;
	//		}

			try { r = new char[strlen(p)+1]; } catch(...) { r = NULL; }
			if (r == NULL) NewException((double)(strlen(p)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(r,p);
			m_paLabels.Add(r);
		}
		p = q+1;

		while ((*p != '.') && (*p != 0))
			p++;

		while (*p != ' ')
			p--;

		p++;

		q = p;
		b = false;
		while (*q != 0)
		{
			if (b)
			{
				if ((*q == ' ') || (*q == '-'))
					break;
			}
			if (*q == 0)
				break;
			if ((*q == '-') || (*q == '.') || ((*q >= '0') && (*q <= '9')))
				b = true;
			q++;
		}
//		q = strchr(p,' ');
		if (*q == 0)
		{
			eprintf("\nCTimeStep::ReadPDB(): Error 1. \"%s\".",obuf);
			return false;
		}
		memcpy(buf2,p,q-p);
		buf2[q-p] = 0;
	//	x = (double)atof(buf2) * 100.0;
		#ifdef FAST_ATOF
			x = fast_atof(buf2) * 100.0;
		#else
			x = atof(buf2) * 100.0;
		#endif
//		mprintf("\n  \"%s\".",buf2);

		p = q;
		b = false;
		while (*q != 0)
		{
			if (b)
			{
				if ((*q == ' ') || (*q == '-'))
					break;
			}
			if (*q == 0)
				break;
			if ((*q == '-') || (*q == '.') || ((*q >= '0') && (*q <= '9')))
				b = true;
			q++;
		}
//		q = strchr(p,' ');
		if (*q == 0)
		{
			eprintf("\nCTimeStep::ReadPDB(): Error 2. \"%s\".",obuf);
			return false;
		}
		memcpy(buf2,p,q-p);
		buf2[q-p] = 0;
		//y = (double)atof(buf2) * 100.0;
		#ifdef FAST_ATOF
			y = fast_atof(buf2) * 100.0;
		#else
			y = atof(buf2) * 100.0;
		#endif
//		mprintf("\n  \"%s\".",buf2);

		//z = (double)atof(q) * 100.0;
		#ifdef FAST_ATOF
			z = fast_atof(q) * 100.0;
		#else
			z = atof(q) * 100.0;
		#endif
//		mprintf("\n  \"%s\".",q);

//		mprintf("\n%f, %f, %f",x,y,z);

		m_vaCoords.Add(CxDVector3(x,y,z));

		if (g_bKeepOriginalCoords)
			m_vaCoords_Original.Add(CxDVector3(x,y,z));
	}
_ende:
	m_iGesAtomCount = m_vaCoords.GetSize();

	return true;
}


bool CTimeStep::SkipPDB(FILE *a)
{
	char buf[256];
//	bool b;

	while (true)
	{
		(void)fgets_bin(buf,256,a);
		if (feof(a))
			return false;
		if (strlen(buf) == 0)
			continue;
		buf[strlen(buf)-1] = 0;
		if (strstr(buf,"END") == buf)
			break;
	}
	return true;
}



bool CTimeStep::ReadGRO(FILE *a, bool needinfo) {

	static char buf[256], buf2[64];
	int i, z;
	char *p, *q, *r;
	double fx, fy, fz;
	std::vector<double> cell;
	CxDVector3 veca, vecb, vecc;


	m_vaCoords.RemoveAll_KeepSize();
	for (i=0;i<m_paLabels.GetSize();i++)
		delete[] (char*)m_paLabels[i];
	m_paLabels.RemoveAll();
	if (g_bKeepOriginalCoords)
		m_vaCoords_Original.RemoveAll_KeepSize();

	if (m_pComment == NULL) {
		try { m_pComment = new char[256]; } catch(...) { m_pComment = NULL; }
		if (m_pComment == NULL) NewException((double)256*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}
	m_pComment[0] = 0;
	(void)fgets_bin(m_pComment,255,a);
	if (feof(a))
		return false;
	if (strlen(m_pComment) >= 254)
		eprintf("\nCTimeStep::ReadGRO(): Warning: Comment line longer than 254 characters.  ");
	if (strlen(m_pComment) != 0)
		m_pComment[strlen(m_pComment)-1] = 0;

	(void)fgets_bin(buf,255,a);
	m_iGesAtomCount = atoi(buf);
	if (atoi(buf) < 0) {
		eprintf("\nCTimeStep::ReadGRO(): Error: Atom count = %d < 0. \"%s\"",atoi(buf),buf);
		m_iGesAtomCount = 0;
		return false;
	}
	if (m_iGesAtomCount == 0) {
		eprintf("\nCTimeStep::ReadGRO(): Error: Atom count = 0 (\"%s\"). Your trajectory file is damaged.",buf);
		return false;
	}

	if (needinfo) {
		if (g_bDoubleBox)
			m_paLabels.SetSize(m_iGesAtomCount*g_iDoubleBoxFactor);
		else
			m_paLabels.SetSize(m_iGesAtomCount);
	}
	if (m_vaCoords.GetSize() < (long)m_iGesAtomCount)
		m_vaCoords.SetSize(m_iGesAtomCount);

	if (g_bKeepOriginalCoords)
		if (m_vaCoords_Original.GetSize() < (long)m_iGesAtomCount)
			m_vaCoords_Original.SetSize(m_iGesAtomCount);

	for (z=0;z<(int)m_iGesAtomCount;z++) {

		(void)fgets_bin(buf,255,a);

		if (feof(a)) {
			eprintf("\nCTimeStep::ReadGRO(): Error: Unexpected end of file while reading atom %d/%lu.",z+1,m_iGesAtomCount);
			return false;
		}

		// Field 1: Residue number (%5d)
		p = buf;
		memcpy(buf2,p,5);
		buf2[5] = 0;

		// Field 2: Residue name (%5s)
		p += 5;
		memcpy(buf2,p,5);
		buf2[5] = 0;

		// Field 3: Atom name (%5s)
		p += 5;
		memcpy(buf2,p,5);
		buf2[5] = 0;

		if (needinfo) {
			q = buf2;
			while (*q == ' ')
				q++;
			try { r = new char[strlen(q)+1]; } catch(...) { r = NULL; }
			if (r == NULL) NewException((double)(strlen(q)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			strcpy(r,q);
			m_paLabels[z] = r;
		}

		// Field 4: Atom number (%5d)
		p += 5;
		memcpy(buf2,p,5);
		buf2[5] = 0;

		// Field 5: Atom position X in nm (%8.3f)
		p += 5;
		memcpy(buf2,p,8);
		buf2[8] = 0;

		#ifdef FAST_ATOF
			fx = fast_atof(buf2) * 1000.0;
		#else
			fx = atof(buf2) * 1000.0;
		#endif

		// Field 6: Atom position Y in nm (%8.3f)
		p += 8;
		memcpy(buf2,p,8);
		buf2[8] = 0;

		#ifdef FAST_ATOF
			fy = fast_atof(buf2) * 1000.0;
		#else
			fy = atof(buf2) * 1000.0;
		#endif

		// Field 7: Atom position Z in nm (%8.3f)
		p += 8;
		memcpy(buf2,p,8);
		buf2[8] = 0;

		#ifdef FAST_ATOF
			fz = fast_atof(buf2) * 1000.0;
		#else
			fz = atof(buf2) * 1000.0;
		#endif

		m_vaCoords[z] = CxDVector3(fx,fy,fz);

		if (g_bKeepOriginalCoords)
			m_vaCoords_Original[z] = CxDVector3(fx,fy,fz);
	}

	// Cell vector
	(void)fgets_bin(buf,255,a);

	p = buf;
	while (true) {
		while ((strchr(" \t\r\n",*p) != NULL) && (*p != 0))
			p++;
		if (*p == 0)
			break;
		q = p;
		while ((strchr(" \t\r\n",*q) == NULL) && (*q != 0))
			q++;
		memcpy(buf2,p,q-p);
		buf2[q-p] = 0;
		cell.push_back( atof(buf2) );
		if (*q == 0)
			break;
		p = q;
	}

	if ((cell.size() != 3) && (cell.size() != 9)) {
		eprintf("CTimeStep::ReadGRO(): Error: Expected three or nine cell vector entries, found %lu.\n",(unsigned long)cell.size());
		return false;
	}

	if ((cell.size() == 9) && ((cell[3] != 0) || (cell[4] != 0) || (cell[5] != 0) || (cell[6] != 0) || (cell[7] != 0) || (cell[8] != 0))) {

		g_bFoundNonOrtho = true;

		g_mBoxFromOrtho(0,0) = cell[0] * 1000.0;
		g_mBoxFromOrtho(0,1) = cell[3] * 1000.0;
		g_mBoxFromOrtho(0,2) = cell[4] * 1000.0;

		g_mBoxFromOrtho(1,0) = cell[5] * 1000.0;
		g_mBoxFromOrtho(1,1) = cell[1] * 1000.0;
		g_mBoxFromOrtho(1,2) = cell[6] * 1000.0;

		g_mBoxFromOrtho(2,0) = cell[7] * 1000.0;
		g_mBoxFromOrtho(2,1) = cell[8] * 1000.0;
		g_mBoxFromOrtho(2,2) = cell[2] * 1000.0;

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

		// Backward transformation matrix
		g_mBoxToOrtho = CxDMatrix3(g_mBoxFromOrtho);
		if (!g_mBoxToOrtho.Invert()) {
			eprintf("CTimeStep::ReadGRO(): Error: Encountered singular cell matrix (cell volume is zero).\n");
			abort();
		}

		veca = CxDVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
		vecb = CxDVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
		vecc = CxDVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

		// Cell angles
		g_fBoxAngleA = acos(DotP(vecb,vecc) / vecb.GetLength() / vecc.GetLength()) * 180.0 / Pi;
		g_fBoxAngleB = acos(DotP(veca,vecc) / veca.GetLength() / vecc.GetLength()) * 180.0 / Pi;
		g_fBoxAngleC = acos(DotP(veca,vecb) / veca.GetLength() / vecb.GetLength()) * 180.0 / Pi;

		// Orthogonal bounding box
		g_fBoxX = fabs(g_mBoxFromOrtho(0,0)) + fabs(g_mBoxFromOrtho(1,0)) + fabs(g_mBoxFromOrtho(2,0));
		g_fBoxY = fabs(g_mBoxFromOrtho(0,1)) + fabs(g_mBoxFromOrtho(1,1)) + fabs(g_mBoxFromOrtho(2,1));
		g_fBoxZ = fabs(g_mBoxFromOrtho(0,2)) + fabs(g_mBoxFromOrtho(1,2)) + fabs(g_mBoxFromOrtho(2,2));

		// Minimal diameters
		g_fBoxMinDiamA = fabs(DotP(veca,Normalize(CrossP(vecb,vecc))));
		g_fBoxMinDiamB = fabs(DotP(vecb,Normalize(CrossP(veca,vecc))));
		g_fBoxMinDiamC = fabs(DotP(vecc,Normalize(CrossP(veca,vecb))));
		g_fBoxMinDiam = MIN3(g_fBoxMinDiamA,g_fBoxMinDiamB,g_fBoxMinDiamC);

	} else {

		g_fBoxX = cell[0] * 1000.0;
		g_fBoxY = cell[1] * 1000.0;
		g_fBoxZ = cell[2] * 1000.0;

		if (g_bDoubleBox) {
			g_fBoxX *= g_iDoubleBoxX;
			g_fBoxY *= g_iDoubleBoxY;
			g_fBoxZ *= g_iDoubleBoxZ;
		}

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
			eprintf("CTimeStep::ReadGRO(): Error: Encountered singular cell matrix (cell volume is zero).\n");
			abort();
		}
	}

	return true;
}



bool CTimeStep::SkipGRO(FILE *a) {

	static char buf[256];
	int z, i;

	(void)fgets_bin(buf,255,a);
	(void)fgets_bin(buf,255,a);
	i = atoi(buf);
	for (z=0;z<i;z++)
		(void)fgets_bin(buf,255,a);
	(void)fgets_bin(buf,255,a);

	if (feof(a))
		return false;

	return true;
}



bool CTimeStep::ReadDCDHeader(FILE *a) {

	char head[4];
	char buf[256];
	int z, ti, nframes, nstart, nremark, natoms;
	float dt;


//	mprintf(">>> ReadDCDHeader >>>\n");

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 84) {
		eprintf("CTimeStep::ReadDCDHeader(): Format error (1): Expected 84, found %d.\n",ti);
		return false;
	}
//	mprintf("@ 84=%d\n",ti);

	(void)!fread(head,4,1,a);
//	mprintf("@ head=\"%c%c%c%c\"\n",head[0],head[1],head[2],head[3]);
	if (memcmp(head,"CORD",4) != 0) {
		eprintf("CTimeStep::ReadDCDHeader(): Error: Found \"%c%c%c%c\" instead of \"CORD\".\n",head[0],head[1],head[2],head[3]);
		return false;
	}

	nframes = 0;
	(void)!fread(&nframes,4,1,a);
//	mprintf("@ nframes=%d\n",nframes);

	nstart = 0;
	(void)!fread(&nstart,4,1,a);
//	mprintf("@ nstart=%d\n",nstart);

	ti = 0;
	(void)!fread(&ti,4,1,a);
//	mprintf("@ stride=%d\n",ti);

	ti = 0;
	(void)!fread(&ti,4,1,a);
//	mprintf("@ %d\n",ti);

	ti = 0;
	(void)!fread(&ti,4,1,a);
//	mprintf("@ %d\n",ti);

	for (z=0;z<4;z++)
		(void)!fread(&ti,4,1,a);

	(void)!fread(&dt,4,1,a);
//	mprintf("@ dt=%f\n",dt);

	ti = 0;
	(void)!fread(&ti,4,1,a);
//	mprintf("@ haveunitcell=%d\n",ti);

	for (z=0;z<8;z++)
		(void)!fread(&ti,4,1,a);

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 24) {
		eprintf("CTimeStep::ReadDCDHeader(): Format error (2): Expected 24, found %d.\n",ti);
		return false;
	}
//	mprintf("@ 24=%d\n",ti);

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 84) {
		eprintf("CTimeStep::ReadDCDHeader(): Format error (3): Expected 84, found %d.\n",ti);
		return false;
	}
//	mprintf("@ 84=%d\n",ti);

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 164) {
		eprintf("CTimeStep::ReadDCDHeader(): Format error (4): Expected 164, found %d.\n",ti);
		return false;
	}
//	mprintf("@ 164=%d\n",ti);

	nremark = 0;
	(void)!fread(&nremark,4,1,a);
//	mprintf("@ nremark=%d\n",nremark);

	g_saDCDRemarkLines.resize( nremark );

	for (z=0;z<nremark;z++) {
		(void)!fread(buf,80,1,a);
		buf[80] = 0;
		g_saDCDRemarkLines[z] = (const char*)buf;
//		mprintf("@ remark%d = \"%s\"\n",z+1,buf);
	}

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 164) {
		eprintf("CTimeStep::ReadDCDHeader(): Format error (5): Expected 164, found %d.\n",ti);
		return false;
	}
//	mprintf("@ 164=%d\n",ti);

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 4) {
		eprintf("CTimeStep::ReadDCDHeader(): Format error (6): Expected 4, found %d.\n",ti);
		return false;
	}
//	mprintf("@ 4=%d\n",ti);

	natoms = 0;
	(void)!fread(&natoms,4,1,a);
	m_iGesAtomCount = natoms;
//	mprintf("@ natoms=%d\n",natoms);

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 4) {
		eprintf("CTimeStep::ReadDCDHeader(): Format error (7): Expected 4, found %d.\n",ti);
		return false;
	}

	if (feof(a)) {
		eprintf("CTimeStep::ReadDCDHeader(): Error: Unexpected end of file while reading header.\n");
		return false;
	}

//	mprintf("@ 4=%d\n",ti);

//	mprintf("<<< ReadDCDHeader <<<\n");

	return true;
}



bool CTimeStep::SkipDCDHeader(FILE *a) {

	char head[4];
	char buf[256];
	int z, ti, nframes, nstart, nremark, natoms;
	float dt;


	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 84) {
		eprintf("CTimeStep::SkipDCDHeader(): Format error (1): Expected 84, found %d.\n",ti);
		return false;
	}

	(void)!fread(head,4,1,a);
	if (memcmp(head,"CORD",4) != 0) {
		eprintf("CTimeStep::SkipDCDHeader(): Error: Found \"%c%c%c%c\" instead of \"CORD\".\n",head[0],head[1],head[2],head[3]);
		return false;
	}

	nframes = 0;
	(void)!fread(&nframes,4,1,a);

	nstart = 0;
	(void)!fread(&nstart,4,1,a);

	ti = 0;
	(void)!fread(&ti,4,1,a);

	ti = 0;
	(void)!fread(&ti,4,1,a);

	ti = 0;
	(void)!fread(&ti,4,1,a);

	for (z=0;z<4;z++)
		(void)!fread(&ti,4,1,a);

	(void)!fread(&dt,4,1,a);

	ti = 0;
	(void)!fread(&ti,4,1,a);

	for (z=0;z<8;z++)
		(void)!fread(&ti,4,1,a);

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 24) {
		eprintf("CTimeStep::SkipDCDHeader(): Format error (2): Expected 24, found %d.\n",ti);
		return false;
	}

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 84) {
		eprintf("CTimeStep::SkipDCDHeader(): Format error (3): Expected 84, found %d.\n",ti);
		return false;
	}

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 164) {
		eprintf("CTimeStep::SkipDCDHeader(): Format error (4): Expected 164, found %d.\n",ti);
		return false;
	}

	nremark = 0;
	(void)!fread(&nremark,4,1,a);

	for (z=0;z<nremark;z++) {
		(void)!fread(buf,80,1,a);
		buf[80] = 0;
	}

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 164) {
		eprintf("CTimeStep::SkipDCDHeader(): Format error (5): Expected 164, found %d.\n",ti);
		return false;
	}

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 4) {
		eprintf("CTimeStep::SkipDCDHeader(): Format error (6): Expected 4, found %d.\n",ti);
		return false;
	}

	natoms = 0;
	(void)!fread(&natoms,4,1,a);
	//m_iGesAtomCount = natoms;

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 4) {
		eprintf("CTimeStep::SkipDCDHeader(): Format error (7): Expected 4, found %d.\n",ti);
		return false;
	}

	if (feof(a)) {
		eprintf("CTimeStep::SkipDCDHeader(): Error: Unexpected end of file while skipping header.\n");
		return false;
	}

	return true;
}



bool CTimeStep::ReadDCD(FILE *a, bool needinfo) {

	int ti, z;
	double fa, fb, fc, alpha, beta, gamma;
	float tf;
	CxDVector3 veca, vecb, vecc;
	CxString buf;
	FILE *b;


//	mprintf(">>> ReadDCD >>>\n");

	m_paLabels.RemoveAll();

	if (g_bDCDFirstFrame) {
		if (!ReadDCDHeader( a ))
			return false;
		g_bDCDFirstFrame = false;

		mprintf("\n");
		mprintf("    You are opening a DCD file. DCD files do not contain any atom labels, but TRAVIS requires them.\n");
		mprintf("    You need to supply a single .xyz snapshot with the correct atom labels (coordinates are ignored).\n\n");

_again:
		AskString_ND("    Enter the file name of the .xyz frame with the atom labels: ",&buf);

		b = fopen((const char*)buf,"rb");
		if (b == NULL) {
			eprintf("Error: Could not open file \"%s\" for reading.\n",(const char*)buf);
			goto _again;
		}

		CTimeStep tempts;

		g_pDCDLabelsTS = new CTimeStep();

		if (!g_pDCDLabelsTS->ReadXYZ( b, true, &g_pDCDLabelsTS->m_vaCoords )) {
			eprintf("Error: Not a valid XYZ file.\n");
			fclose(b);
			delete g_pDCDLabelsTS;
			goto _again;
		}

		if (g_pDCDLabelsTS->m_iGesAtomCount != m_iGesAtomCount) {
			eprintf("Error: Atom count mismatch (%lu != %lu).\n",g_pDCDLabelsTS->m_iGesAtomCount,m_iGesAtomCount);
			fclose(b);
			delete g_pDCDLabelsTS;
			goto _again;
		}

		fclose(b);

	} else
		m_iGesAtomCount = g_iGesAtomCount;

	if (g_bDCDSkipHeader) {
		if (!SkipDCDHeader( a ))
			return false;
		g_bDCDSkipHeader = false;
	}

	if (needinfo) {
		if (needinfo) {
			if (g_bDoubleBox)
				m_paLabels.SetSize(m_iGesAtomCount*g_iDoubleBoxFactor);
			else
				m_paLabels.SetSize(m_iGesAtomCount);
		}
		for (z=0;z<(int)m_iGesAtomCount;z++)
			m_paLabels[z] = g_pDCDLabelsTS->m_paLabels[z];
	}


	if (m_vaCoords.GetSize() < (long)m_iGesAtomCount)
		m_vaCoords.SetSize(m_iGesAtomCount);

	if (g_bKeepOriginalCoords)
		if (m_vaCoords_Original.GetSize() < (long)m_iGesAtomCount)
			m_vaCoords_Original.SetSize(m_iGesAtomCount);

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (feof(a))
		return false;
	if (ti != 48) {
		eprintf("CTimeStep::ReadDCD(): Format error (1): Expected 48, found %d.\n",ti);
		return false;
	}
//	mprintf("@ 48=%d\n",ti);

	(void)!fread(&fa,8,1,a);
	(void)!fread(&gamma,8,1,a);
	(void)!fread(&fb,8,1,a);
	(void)!fread(&beta,8,1,a);
	(void)!fread(&alpha,8,1,a);
	(void)!fread(&fc,8,1,a);

	if (g_bNPT) {

		if ((fabs(fa) < 0.01) || (fabs(fb) < 0.01) || (fabs(fc) < 0.01)) // Zero cell vector --> No Data
			goto _nocell;

		if ((fabs(alpha-90.0) > 0.001) || (fabs(beta-90.0) > 0.001) || (fabs(gamma-90.0) > 0.001)) {

			g_bFoundNonOrtho = true;

			g_fBoxX = fa * 100.0; // Angstrom --> pm
			g_fBoxY = fb * 100.0;
			g_fBoxZ = fc * 100.0;
			g_fBoxAngleA = alpha;
			g_fBoxAngleB = beta;
			g_fBoxAngleC = gamma;

			g_mBoxFromOrtho(0,0) = g_fBoxX;
			g_mBoxFromOrtho(0,1) = 0;
			g_mBoxFromOrtho(0,2) = 0;

			g_mBoxFromOrtho(1,0) = g_fBoxY*cos(g_fBoxAngleC*Pi/180.0);
			g_mBoxFromOrtho(1,1) = g_fBoxY*sin(g_fBoxAngleC*Pi/180.0);
			g_mBoxFromOrtho(1,2) = 0;

			g_mBoxFromOrtho(2,0) = g_fBoxZ*cos(g_fBoxAngleB*Pi/180.0);
			g_mBoxFromOrtho(2,1) = (-(g_mBoxFromOrtho(1,0)*g_fBoxZ*cos(g_fBoxAngleB*Pi/180.0)) + g_fBoxY*g_fBoxZ*cos(g_fBoxAngleA*Pi/180.0))/g_mBoxFromOrtho(1,1);
			g_mBoxFromOrtho(2,2) = sqrt(-((pow2(g_fBoxZ)*(pow2(g_mBoxFromOrtho(1,1))*(-1 + pow2(cos(g_fBoxAngleB*Pi/180.0))) + pow2(g_mBoxFromOrtho(1,0)*cos(g_fBoxAngleB*Pi/180.0) - g_fBoxY*cos(g_fBoxAngleA*Pi/180.0))))/pow2(g_mBoxFromOrtho(1,1))));

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
				eprintf("CTimeStep::ReadDCD(): Error: Encountered singular cell matrix (cell volume is zero).\n");
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

		} else {

			g_fBoxX = fa*100.0;
			g_fBoxY = fb*100.0;
			g_fBoxZ = fc*100.0;

			if (g_bDoubleBox) {
				g_fBoxX *= g_iDoubleBoxX;
				g_fBoxY *= g_iDoubleBoxY;
				g_fBoxZ *= g_iDoubleBoxZ;
			}

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
				eprintf("CTimeStep::ReadDCD(): Error: Encountered singular cell matrix (cell volume is zero).\n");
				abort();
			}
		}
	}

_nocell:
	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 48) {
		eprintf("CTimeStep::ReadDCD(): Format error (2): Expected 48, found %d.\n",ti);
		return false;
	}
//	mprintf("@ 48=%d\n",ti);

//	mprintf("@ Cell a=%f, b=%f, c=%f, alpha=%f, beta=%f, gamma=%f\n",fa,fb,fc,alpha,beta,gamma);

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::ReadDCD(): Format error (3): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}
//	mprintf("@ nx=%d\n",ti);
	for (z=0;z<(int)m_iGesAtomCount;z++) {
		(void)!fread(&tf,4,1,a);
		m_vaCoords[z][0] = tf * 100.0;
		//mprintf("@  X%3d: %16.8f\n",z+1,tf);
	}
	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::ReadDCD(): Format error (4): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}
//	mprintf("@ nx=%d\n",ti);

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::ReadDCD(): Format error (5): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}
//	mprintf("@ nx=%d\n",ti);
	for (z=0;z<(int)m_iGesAtomCount;z++) {
		(void)!fread(&tf,4,1,a);
		m_vaCoords[z][1] = tf * 100.0;
		//mprintf("@  Y%3d: %16.8f\n",z+1,tf);
	}
	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::ReadDCD(): Format error (6): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}
//	mprintf("@ nx=%d\n",ti);

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::ReadDCD(): Format error (7): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}
//	mprintf("@ nx=%d\n",ti);
	for (z=0;z<(int)m_iGesAtomCount;z++) {
		(void)!fread(&tf,4,1,a);
		m_vaCoords[z][2] = tf * 100.0;
		//mprintf("@  Z%3d: %16.8f\n",z+1,tf);
	}
	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::ReadDCD(): Format error (8): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}
//	mprintf("@ nx=%d\n",ti);

	if (feof(a)) {
		eprintf("CTimeStep::ReadDCD(): Error: Unexpected end of file while reading frame.\n");
		return false;
	}

	if (g_bKeepOriginalCoords)
		for (z=0;z<(int)m_iGesAtomCount;z++)
			m_vaCoords_Original[z] = m_vaCoords[z];

//	mprintf("<<< ReadDCD <<<\n");

	return true;
}



bool CTimeStep::SkipDCD(FILE *a) {

	int ti, z;
	double fa, fb, fc, alpha, beta, gamma;
	float tf;


	if (g_bDCDFirstFrame) {
		if (!ReadDCDHeader( a ))
			return false;
		g_bDCDFirstFrame = false;
	}

	if (g_bDCDSkipHeader) {
		if (!SkipDCDHeader( a ))
			return false;
		g_bDCDSkipHeader = false;
	}

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (feof(a))
		return false;
	if (ti != 48) {
		eprintf("CTimeStep::SkipDCD(): Format error (1): Expected 48, found %d.\n",ti);
		return false;
	}

	(void)!fread(&fa,8,1,a);
	(void)!fread(&gamma,8,1,a);
	(void)!fread(&fb,8,1,a);
	(void)!fread(&beta,8,1,a);
	(void)!fread(&alpha,8,1,a);
	(void)!fread(&fc,8,1,a);

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != 48) {
		eprintf("CTimeStep::SkipDCD(): Format error (2): Expected 48, found %d.\n",ti);
		return false;
	}

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::SkipDCD(): Format error (3): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}
	for (z=0;z<(int)m_iGesAtomCount;z++)
		(void)!fread(&tf,4,1,a);
	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::SkipDCD(): Format error (4): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::SkipDCD(): Format error (5): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}
	for (z=0;z<(int)m_iGesAtomCount;z++)
		(void)!fread(&tf,4,1,a);
	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::SkipDCD(): Format error (6): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}

	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::SkipDCD(): Format error (7): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}
	for (z=0;z<(int)m_iGesAtomCount;z++)
		(void)!fread(&tf,4,1,a);
	ti = 0;
	(void)!fread(&ti,4,1,a);
	if (ti != (int)m_iGesAtomCount*4) {
		eprintf("CTimeStep::SkipDCD(): Format error (8): Expected %lu, found %d.\n",m_iGesAtomCount*4,ti);
		return false;
	}

	if (feof(a)) {
		eprintf("CTimeStep::SkipDCD(): Error: Unexpected end of file while skipping frame.\n");
		return false;
	}

	return true;
}



bool CTimeStep::ReadVoronoi(FILE *a, bool needinfo) {

	static char buf[512];
	int z, z2, z3, ti;
	char *p, *q, *r;
	const char *tl;
	double tf;


	(void)fgets_bin(buf,511,a);
	if (feof(a))
		return false;
	if (!IsValidInteger(buf)) {
		eprintf("\nCTimeStep::ReadVoronoi(): Error: The atom number field is not a valid integer (\"%s\").\n",buf);
		eprintf("Maybe the .voronoi trajectory was written with CP2k 8.1. The format has changed in May 2021 - sorry.\n");
		eprintf("Contact the developers if you need to convert your files to the new format.\n");
		m_iGesAtomCount = 0;
		return false;
	}
	m_iGesAtomCount = atoi(buf);
	if (atoi(buf) < 0) {
		eprintf("\nCTimeStep::ReadVoronoi(): Error: Atom count = %d < 0. \"%s\"",atoi(buf),buf);
		m_iGesAtomCount = 0;
		return false;
	}
	if (m_iGesAtomCount == 0) {
		eprintf("\nCTimeStep::ReadVoronoi(): Error: Atom count = 0 (\"%s\"). Your trajectory file is damaged.",buf);
		return false;
	}

	m_vaCoords.RemoveAll_KeepSize();
	for (z=0;z<m_paLabels.GetSize();z++)
		delete[] (char*)m_paLabels[z];
	m_paLabels.RemoveAll();
	if (g_bKeepOriginalCoords)
		m_vaCoords_Original.RemoveAll_KeepSize();

	if (m_pComment == NULL) {
		try { m_pComment = new char[256]; } catch(...) { m_pComment = NULL; }
		if (m_pComment == NULL) NewException((double)256*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	}
	m_pComment[0] = 0;
	// First comment line
	(void)fgets_bin(m_pComment,511,a);
	if (feof(a))
		return false;
	if (strlen(m_pComment) >= 511)
		eprintf("\nCTimeStep::ReadVoronoi(): Warning: Comment line longer than 511 characters.  ");
	if (strlen(m_pComment) != 0)
		m_pComment[strlen(m_pComment)-1] = 0;

	// Second comment line: Cell
	(void)fgets_bin(buf,511,a);
	if (strlen(buf) >= 511)
		eprintf("\nCTimeStep::ReadVoronoi(): Warning: Comment line longer than 511 characters.  ");
	if (g_bNPT) {
		p = buf;
		if (*p != '#') {
			eprintf("\nCTimeStep::ReadVoronoi(): Error: Cell geometry line is in invalid format (1).\n");
			mprintf("  \"%s\"\n",buf);
			return false;
		}
		p++;
		while (*p == ' ')
			p++;
		q = p;
		while ((*q != ' ') && (*q != 0))
			q++;
		if (q-p != 4) {
			eprintf("\nCTimeStep::ReadVoronoi(): Error: Cell geometry line is in invalid format (2).\n");
			mprintf("  \"%s\"\n",buf);
			return false;
		}
		if (memcmp(p,"Cell",4) != 0) {
			eprintf("\nCTimeStep::ReadVoronoi(): Error: Cell geometry line is in invalid format (3).\n");
			mprintf("  \"%s\"\n",buf);
			return false;
		}
		p = q;
		g_bFoundNonOrtho = true;
		ExtractXYZCellGeometry9(p);
	}

	// Third comment line
	(void)fgets_bin(buf,511,a);
	if (strlen(buf) >= 511)
		eprintf("\nCTimeStep::ReadVoronoi(): Warning: Comment line longer than 511 characters.  ");

	if (needinfo) {
		if (g_bDoubleBox)
			m_paLabels.SetSize(m_iGesAtomCount*g_iDoubleBoxFactor);
		else
			m_paLabels.SetSize(m_iGesAtomCount);
	}

	if (m_vaCoords.GetSize() < (long)m_iGesAtomCount)
		m_vaCoords.SetSize(m_iGesAtomCount);

	if (g_bKeepOriginalCoords)
		if (m_vaCoords_Original.GetSize() < (long)m_iGesAtomCount)
			m_vaCoords_Original.SetSize(m_iGesAtomCount);

	if (g_faVoronoiRadii.GetSize() < (int)m_iGesAtomCount)
		g_faVoronoiRadii.SetSize(m_iGesAtomCount);

	if (m_faVolume.GetSize() < (int)m_iGesAtomCount)
		m_faVolume.SetSize(m_iGesAtomCount);

	if (m_chargeCenters.GetSize() < (int)m_iGesAtomCount)
		m_chargeCenters.SetSize(m_iGesAtomCount);

	if (m_maQuadTensor.size() < m_iGesAtomCount)
		m_maQuadTensor.resize(m_iGesAtomCount);

	if (m_dipoleMoments.GetSize() < (int)m_iGesAtomCount)
		m_dipoleMoments.SetSize(m_iGesAtomCount);

	if (m_faCharge.GetSize() < (int)m_iGesAtomCount)
		m_faCharge.SetSize(m_iGesAtomCount);

	if (m_faCoreCharge.GetSize() < (int)m_iGesAtomCount)
		m_faCoreCharge.SetSize(m_iGesAtomCount);


	for (z=0;z<(int)m_iGesAtomCount;z++) {

		(void)fgets_bin(buf,511,a);
		if (strlen(buf) >= 511)
			eprintf("\nCTimeStep::ReadVoronoi(): Warning: Data line longer than 511 characters.  ");

		if (feof(a)) {
			eprintf("\nCTimeStep::ReadVoronoi(): Error: Unexpected end of file while reading atom %d/%lu.",z+1,m_iGesAtomCount);
			return false;
		}

		p = buf;
		for (z2=0;z2<24;z2++) { // Loop over columns

			while ((*p == ' ') || (*p == '\t') || (*p == ';'))
				p++;
			q = p;
			while ((*q != ' ') && (*q != '\t') && (*q != ';') && (*q != '\r') && (*q != '\n') && (*q != 0))
				q++;
			if (z2 < 23) {
				if ((*q == 0) || (*q == '\r') || (*q == '\n')) {
					eprintf("\nCTimeStep::ReadVoronoi(): Error: Incomplete line while reading atom %d/%lu:",z+1,m_iGesAtomCount);
					eprintf("  \"%s\"\n",buf);
					return false;
				}
			}
			*q = 0;

			switch(z2) {

				case 0: // Column 1: ID
					if (atoi(p) != z+1) {
						eprintf("\nCTimeStep::ReadVoronoi(): Error: Unexpected atom ID while reading atom %d/%lu: %d",z+1,m_iGesAtomCount,atoi(p));
						return false;
					}
					break;

				case 1: // Column 2: Ord
					if (needinfo) {
						ti = atoi(p);
						for (z3=0;z3<g_oaElements.GetSize();z3++) {
							if (((CElement*)g_oaElements[z3])->m_iOrd == ti) {
								tl = ((CElement*)g_oaElements[z3])->m_sLabel;
								goto _lfound;
							}
						}
						tl = "X"; // Unknown element
_lfound:
						try { r = new char[strlen(tl)+1]; } catch(...) { r = NULL; }
						if (r == NULL) NewException((double)(strlen(tl)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						strcpy(r,tl);
						m_paLabels[z] = r;
					}
					break;

				case 2: // Column 3: Radius / Angstrom
					#ifdef FAST_ATOF
						tf = fast_atof(p) * 100.0;
					#else
						tf = atof(p) * 100.0;
					#endif
					g_faVoronoiRadii[z] = tf;
					break;

				case 3: // Column 4: PosX / Angstrom
					#ifdef FAST_ATOF
						tf = fast_atof(p) * 100.0;
					#else
						tf = atof(p) * 100.0;
					#endif
					m_vaCoords[z][0] = tf;
					if (g_bKeepOriginalCoords)
						m_vaCoords_Original[z][0] = tf;
					break;

				case 4: // Column 5: PosY / Angstrom
					#ifdef FAST_ATOF
						tf = fast_atof(p) * 100.0;
					#else
						tf = atof(p) * 100.0;
					#endif
					m_vaCoords[z][1] = tf;
					if (g_bKeepOriginalCoords)
						m_vaCoords_Original[z][1] = tf;
					break;

				case 5: // Column 6: PosZ / Angstrom
					#ifdef FAST_ATOF
						tf = fast_atof(p) * 100.0;
					#else
						tf = atof(p) * 100.0;
					#endif
					m_vaCoords[z][2] = tf;
					if (g_bKeepOriginalCoords)
						m_vaCoords_Original[z][2] = tf;
					break;

				case 6: // Column 7: Volume / Angstrom^3
					#ifdef FAST_ATOF
						tf = fast_atof(p) * 1.0e6;
					#else
						tf = atof(p) * 1.0e6;
					#endif
					m_faVolume[z] = tf;
					break;

				case 7: // Column 8: Z(eff) / e
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_faCoreCharge[z] = tf;
					break;

				case 8: // Column 9: Charge / e
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_faCharge[z] = tf;
					break;

				case 9: // Column 10: Dipole X / e*pm
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_dipoleMoments[z][0] = tf;
					break;

				case 10: // Column 11: Dipole Y / e*pm
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_dipoleMoments[z][1] = tf;
					break;

				case 11: // Column 12: Dipole Z / e*pm
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_dipoleMoments[z][2] = tf;
					break;

				case 12: // Column 13: Charge Center X / Angstrom
					#ifdef FAST_ATOF
						tf = fast_atof(p) * 100.0;
					#else
						tf = atof(p) * 100.0;
					#endif
					m_chargeCenters[z][0] = tf;
					break;

				case 13: // Column 14: Charge Center Y / Angstrom
					#ifdef FAST_ATOF
						tf = fast_atof(p) * 100.0;
					#else
						tf = atof(p) * 100.0;
					#endif
					m_chargeCenters[z][1] = tf;
					break;

				case 14: // Column 15: Charge Center Z / Angstrom
					#ifdef FAST_ATOF
						tf = fast_atof(p) * 100.0;
					#else
						tf = atof(p) * 100.0;
					#endif
					m_chargeCenters[z][2] = tf;
					break;

				case 15: // Column 16: Quadrupole Tensor XX / e*nm^2
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_maQuadTensor[z][0] = tf;
					break;

				case 16: // Column 17: Quadrupole Tensor XY / e*nm^2
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_maQuadTensor[z][1] = tf;
					break;

				case 17: // Column 18: Quadrupole Tensor XZ / e*nm^2
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_maQuadTensor[z][2] = tf;
					break;

				case 18: // Column 19: Quadrupole Tensor YX / e*nm^2
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_maQuadTensor[z][3] = tf;
					break;

				case 19: // Column 20: Quadrupole Tensor YY / e*nm^2
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_maQuadTensor[z][4] = tf;
					break;

				case 20: // Column 21: Quadrupole Tensor YZ / e*nm^2
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_maQuadTensor[z][5] = tf;
					break;

				case 21: // Column 22: Quadrupole Tensor ZX / e*nm^2
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_maQuadTensor[z][6] = tf;
					break;

				case 22: // Column 23: Quadrupole Tensor ZY / e*nm^2
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_maQuadTensor[z][7] = tf;
					break;

				case 23: // Column 24: Quadrupole Tensor ZZ / e*nm^2
					#ifdef FAST_ATOF
						tf = fast_atof(p);
					#else
						tf = atof(p);
					#endif
					m_maQuadTensor[z][8] = tf;
					break;
			}

			p = q+1;
		}
	}

	return true;
}



bool CTimeStep::SkipVoronoi(FILE *a) {

	static char buf[256];
	int z, i;

	(void)fgets_bin(buf,255,a);
	i = atoi(buf);
	// Three comment lines
	(void)fgets_bin(buf,255,a);
	(void)fgets_bin(buf,255,a);
	(void)fgets_bin(buf,255,a);
	for (z=0;z<i;z++)
		(void)fgets_bin(buf,255,a);

	if (feof(a))
		return false;

	return true;
}



bool CTimeStep::AnalyzeLAMMPSColumns( const char *s ) {

	const char *p, *q;
	static char buf[256];
	int i, ti, z;


	if (g_bLAMMPSColumnsAnalyzed) { // Already analyzed before

		if (strcmp(s,g_sLAMMPSColumnString) != 0) {
			eprintf("CTimeStep::ReadLAMMPS(): Error: LAMMPS trajectory column definition has changed within the trajectory.\n");
			eprintf("Original definition: \"%s\"\n",g_sLAMMPSColumnString);
			eprintf("Current definition:  \"%s\"\n",s);
			return false;
		}

		return true;
	}

	g_bLAMMPSColumnsAnalyzed = true;

	if (g_sLAMMPSColumnString != NULL)
		delete[] g_sLAMMPSColumnString;
	g_sLAMMPSColumnString = new char[strlen(s)+1];
	strcpy(g_sLAMMPSColumnString,s);

	p = s;
	i = 0;
	g_iaLAMMPSOffsetToColumn.clear();
	g_iaLAMMPSColumnToOffset.resize( LAMMPS_FIELD_MAX+1 );
	for (z=0;z<=LAMMPS_FIELD_MAX;z++)
		g_iaLAMMPSColumnToOffset[z] = -1;

	mprintf("    Analyzing LAMMPS trajectory column definition:\n");

	while (true) {

		while ((*p == ' ') || (*p == '\t') || (*p == '\r') || (*p == '\n'))
			p++;
		if (*p == 0)
			break;
		q = p;
		while ((*q != ' ') && (*q != '\t') && (*q != '\r') && (*q != '\n') && (*q != 0))
			q++;

		memcpy(buf,p,q-p);
		buf[q-p] = 0;

		if (i == 0) {
			if (strcmp(buf,"ITEM:") != 0) {
				eprintf("CTimeStep::ReadLAMMPS(): Internal error: Missing \"ITEM:\".\n");
				return false;
			}
		} else if (i == 1) {
			if (strcmp(buf,"ATOMS") != 0) {
				eprintf("CTimeStep::ReadLAMMPS(): Internal error: Missing \"ATOMS\".\n");
				return false;
			}
		} else {

			if (mystricmp(buf,"id") == 0)
				ti = LAMMPS_FIELD_ID;
			else if (mystricmp(buf,"element") == 0)
				ti = LAMMPS_FIELD_ELEMENT;
			else if (mystricmp(buf,"xu") == 0)
				ti = LAMMPS_FIELD_XU;
			else if (mystricmp(buf,"yu") == 0)
				ti = LAMMPS_FIELD_YU;
			else if (mystricmp(buf,"zu") == 0)
				ti = LAMMPS_FIELD_ZU;
			else if (mystricmp(buf,"x") == 0)
				ti = LAMMPS_FIELD_X;
			else if (mystricmp(buf,"y") == 0)
				ti = LAMMPS_FIELD_Y;
			else if (mystricmp(buf,"z") == 0)
				ti = LAMMPS_FIELD_Z;
			else if (mystricmp(buf,"vx") == 0)
				ti = LAMMPS_FIELD_VX;
			else if (mystricmp(buf,"vy") == 0)
				ti = LAMMPS_FIELD_VY;
			else if (mystricmp(buf,"vz") == 0)
				ti = LAMMPS_FIELD_VZ;
			else if (mystricmp(buf,"fx") == 0)
				ti = LAMMPS_FIELD_FX;
			else if (mystricmp(buf,"fy") == 0)
				ti = LAMMPS_FIELD_FY;
			else if (mystricmp(buf,"fz") == 0)
				ti = LAMMPS_FIELD_FZ;
			else if (mystricmp(buf,"q") == 0)
				ti = LAMMPS_FIELD_Q;
			else
				ti = -1; // Unknown field: ignore

			g_iaLAMMPSOffsetToColumn.push_back( ti );
			if (ti != -1)
				g_iaLAMMPSColumnToOffset[ti] = i-2;

			mprintf("      - Field %2d: %-8s",i-1,buf);
			if (ti == -1)
				mprintf(" (ignored)");
			mprintf("\n");
		}

		i++;

		if (*q == 0)
			break;
		p = q;
	}

	g_iLAMMPSColumnCount = i-2;
	mprintf("    Found %d columns in total.\n",g_iLAMMPSColumnCount);

	if (g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_ELEMENT ] == -1) {
		mprintf("\n");
		eprintf("    CTimeStep::ReadLAMMPS(): Error: The \"element\" field is missing in the trajectory.\n");
		mprintf("    You should use a LAMMPS dump style similar to the following:\n");
		mprintf("      \"dump custom id element xu yu zu\"\n\n");
		return false;
	}

	if (g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_ID ] == -1) {
		mprintf("\n");
		eprintf("    CTimeStep::ReadLAMMPS(): Warning: The \"id\" field is not present in the trajectory.\n");
		mprintf("    This is fine as long as you use \"dump_modify sort id\" to write your trajectory.\n");
		mprintf("    Otherwise, the atom ordering will be random, and TRAVIS has no chance to detect it,\n");
		mprintf("    which will lead to *wrong* analysis results!\n\n");
		mprintf("    If in doubt, you should use a LAMMPS dump style similar to the following:\n");
		mprintf("      \"dump custom id element xu yu zu\"\n\n");
		if (!AskYesNo("    I have acknowledged this note (y/n). [no] ",false))
			return false;
	}

	if ((g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_XU ] == -1) || (g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_YU ] == -1) || (g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_ZU ] == -1)) {

		if ((g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_X ] == -1) || (g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_Y ] == -1) || (g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_Z ] == -1)) {

			mprintf("\n");
			eprintf("    CTimeStep::ReadLAMMPS(): Error: The trajectory does not contain atom positions.\n");
			mprintf("    You should use a LAMMPS dump style similar to the following:\n");
			mprintf("      \"dump custom id element xu yu zu\"\n\n");
			return false;

		} else {

			mprintf("\n");
			eprintf("    Warning: Your LAMMPS trajectory contains wrapped coordinates.\n");
			eprintf("             Dynamical analyses (like MSD, ...) will yield erroneous results!\n");
			eprintf("             Better use \"dump custom id element xu yu zu\" in future to write your trajectory.\n\n");

			g_bLAMMPSWrappedCoords = true;
		}

	} else
		g_bLAMMPSWrappedCoords = false;


	if ((g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_FX ] != -1) && (g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_FY ] != -1) && (g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_FZ ] != -1))
		g_bLAMMPSForce = true;

	if (g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_FX ] != -1)
		g_bLAMMPSCharge = true;

	return true;
}



bool CTimeStep::ReadLAMMPS(FILE *a, bool needinfo) {

	int i, ti, id;
	char buf[256], obuf[256], *p, *q, *tlabel;
	double x, y, z, fx, fy, fz, vx, vy, vz, fq;
	double lowbound[3], highbound[3], skew[3];
	bool last;
	CxDVector3 veca, vecb, vecc;


	for (i=0;i<m_paLabels.GetSize();i++)
		if (m_paLabels[i] != 0)
			delete[] (char*)m_paLabels[i];
	m_paLabels.RemoveAll();
	m_iGesAtomCount = 0;

	while (true) {

		(void)fgets_bin(buf,256,a);
		if (feof(a))
			return false;

		if (strlen(buf) == 0)
			continue;

		buf[strlen(buf)-1] = 0;
		strcpy(obuf,buf);

		if (g_bNPT) {

			if (strstr(buf,"ITEM: BOX BOUNDS") != 0) { // Boxlaenge

				if (strstr(buf,"xy xz yz") != 0) {

					g_bFoundNonOrtho = true;

					for (i=0;i<3;i++) {

						(void)fgets_bin(buf,256,a);
						if (feof(a))
							return false;

						p = &buf[0];
						q = p;
						while ((*q != ' ') && (*q != 0))
							q++;
						if (*q == 0) {
							eprintf("\nCTimeStep::ReadLAMMPS(): Incomplete line: \"%s\"\n",obuf);
							return false;
						}
						*q = 0;
						lowbound[i] = atof(p) * 100.0;

						p = q+1;
						while (*p == ' ')
							p++;
						q = p;
						while ((*q != ' ') && (*q != 0))
							q++;
						if (*q == 0) {
							eprintf("\nCTimeStep::ReadLAMMPS(): Incomplete line: \"%s\"\n",obuf);
							return false;
						}
						*q = 0;
						highbound[i] = atof(p) * 100.0;

						p = q+1;
						while (*p == ' ')
							p++;
						q = p;
						skew[i] = atof(q) * 100.0;
					}

					// This is the weird way how LAMMPS saves the cell geometry ^^
					// See http://lammps.sandia.gov/doc/Section_howto.html, Section "Triclinic (non-orthogonal) simulation boxes"
					lowbound[0]  -= MIN4(0.0,skew[0],skew[1],skew[0]+skew[1]);
					highbound[0] -= MAX4(0.0,skew[0],skew[1],skew[0]+skew[1]);
					lowbound[1]  -= MIN(0.0,skew[2]);
					highbound[1] -= MAX(0.0,skew[2]);

					g_mBoxFromOrtho(0,0) = highbound[0] - lowbound[0];
					g_mBoxFromOrtho(0,1) = 0;
					g_mBoxFromOrtho(0,2) = 0;

					g_mBoxFromOrtho(1,0) = skew[0];
					g_mBoxFromOrtho(1,1) = highbound[1] - lowbound[1];
					g_mBoxFromOrtho(1,2) = 0;

					g_mBoxFromOrtho(2,0) = skew[1];
					g_mBoxFromOrtho(2,1) = skew[2];
					g_mBoxFromOrtho(2,2) = highbound[2] - lowbound[2];

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

					// Backward transformation matrix
					g_mBoxToOrtho = CxDMatrix3(g_mBoxFromOrtho);
					if (!g_mBoxToOrtho.Invert()) {
						eprintf("CTimeStep::ReadLAMMPS(): Error: Encountered singular cell matrix (cell volume is zero).\n");
						abort();
					}

					veca = CxDVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
					vecb = CxDVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
					vecc = CxDVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

					// Cell angles
					g_fBoxAngleA = acos(DotP(vecb,vecc) / vecb.GetLength() / vecc.GetLength()) * 180.0 / Pi;
					g_fBoxAngleB = acos(DotP(veca,vecc) / veca.GetLength() / vecc.GetLength()) * 180.0 / Pi;
					g_fBoxAngleC = acos(DotP(veca,vecb) / veca.GetLength() / vecb.GetLength()) * 180.0 / Pi;

					// Orthogonal bounding box
					g_fBoxX = fabs(g_mBoxFromOrtho(0,0)) + fabs(g_mBoxFromOrtho(1,0)) + fabs(g_mBoxFromOrtho(2,0));
					g_fBoxY = fabs(g_mBoxFromOrtho(0,1)) + fabs(g_mBoxFromOrtho(1,1)) + fabs(g_mBoxFromOrtho(2,1));
					g_fBoxZ = fabs(g_mBoxFromOrtho(0,2)) + fabs(g_mBoxFromOrtho(1,2)) + fabs(g_mBoxFromOrtho(2,2));

					// Minimal diameters
					g_fBoxMinDiamA = fabs(DotP(veca,Normalize(CrossP(vecb,vecc))));
					g_fBoxMinDiamB = fabs(DotP(vecb,Normalize(CrossP(veca,vecc))));
					g_fBoxMinDiamC = fabs(DotP(vecc,Normalize(CrossP(veca,vecb))));
					g_fBoxMinDiam = MIN3(g_fBoxMinDiamA,g_fBoxMinDiamB,g_fBoxMinDiamC);

				} else {

					for (i=0;i<3;i++) {

						(void)fgets_bin(buf,256,a);
						if (feof(a))
							return false;

						p = &buf[0];
						while (strchr(" ",*p) != NULL)
							p++;
						q = p+1;
						while ((strchr(" ",*q) == NULL) && (*q != 0))
							q++;
						if (*q == 0) {
							eprintf("\nCTimeStep::ReadLAMMPS(): Incomplete line: \"%s\"\n",obuf);
							return false;
						}
						*q = 0;
						x = atof(p);

						p = q+1;
						while (strchr(" ",*p) != NULL)
							p++;

						y = atof(p);

						switch(i) {
							case 0: g_fBoxX = (y-x)*100.0; break;
							case 1: g_fBoxY = (y-x)*100.0; break;
							case 2: g_fBoxZ = (y-x)*100.0; break;
						}
					}

					if (g_bDoubleBox) {
						g_fBoxX *= g_iDoubleBoxX;
						g_fBoxY *= g_iDoubleBoxY;
						g_fBoxZ *= g_iDoubleBoxZ;
					}

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
						eprintf("CTimeStep::ReadLAMMPS(): Error: Encountered singular cell matrix (cell volume is zero).\n");
						abort();
					}
				}
			}
		}

		if (strstr(buf,"ITEM: NUMBER OF ATOMS") != 0) {
			(void)fgets_bin(buf,256,a);
			if (feof(a))
				return false;
			buf[strlen(buf)-1] = 0;
			if (!IsValidInteger(buf)) {
				eprintf("CTimeStep::ReadLAMMPS(): Error: \"ITEM: NUMBER OF ATOMS\" is not a valid integer (\"%s\"). Your trajectory file is damaged.\n",buf);
				m_iGesAtomCount = 0;
				return false;
			}
			m_iGesAtomCount = atoi(buf);
			continue;
		}

		if (strstr(buf,"ITEM: ATOMS") != 0) {

			if (!AnalyzeLAMMPSColumns(buf))
				return false;

			if (m_iGesAtomCount == 0) {
				eprintf("CTimeStep::ReadLAMMPS(): \"ITEM: ATOMS\" before \"ITEM: NUMBER OF ATOMS\".\n");
				return false;
			}

			if (needinfo) {
				if (g_bDoubleBox)
					m_paLabels.SetSize(m_iGesAtomCount*g_iDoubleBoxFactor);
				else
					m_paLabels.SetSize(m_iGesAtomCount);
			}

			if (g_bReadLAMMPSCharges)
				if (m_faCharge.GetSize() < (int)m_iGesAtomCount)
					m_faCharge.SetSize(m_iGesAtomCount);

			if (m_vaCoords.GetSize() < (long)m_iGesAtomCount)
				m_vaCoords.SetSize(m_iGesAtomCount);

			if (g_bLAMMPSForce)
				if (m_vaForces.GetSize() < (long)m_iGesAtomCount)
					m_vaForces.SetSize(m_iGesAtomCount);

			if (g_bLAMMPSVelocity)
				if (m_vaVelocities.GetSize() < (long)m_iGesAtomCount)
					m_vaVelocities.SetSize(m_iGesAtomCount);

			if (g_bKeepOriginalCoords)
				if (m_vaCoords_Original.GetSize() < (long)m_iGesAtomCount)
					m_vaCoords_Original.SetSize(m_iGesAtomCount);

			for (i=0;i<(int)m_iGesAtomCount;i++) {

				(void)fgets_bin(buf,256,a);
				if (feof(a))
					return false;

				ti = 0;
				p = buf;
				last = false;
				x = 0;
				y = 0;
				z = 0;
				vx = 0;
				vy = 0;
				vz = 0;
				fx = 0;
				fy = 0;
				fz = 0;
				fq = 0;
				tlabel = NULL;
				id = 0;

				while (true) {

					while ((*p == ' ') || (*p == '\t') || (*p == '\r') || (*p == '\n'))
						p++;
					if (*p == 0)
						break;
					q = p;
					while ((*q != ' ') && (*q != '\t') && (*q != '\r') && (*q != '\n') && (*q != 0))
						q++;
					if (*q == 0)
						last = true;
					*q = 0;

					if (g_iaLAMMPSOffsetToColumn[ti] != -1) {

						switch( g_iaLAMMPSOffsetToColumn[ti] ) {

							case LAMMPS_FIELD_ID:
								id = atoi( p ) - 1;
								if (id < 0) {
									eprintf("\nCTimeStep::ReadLAMMPS(): Error: Invalid ID value (\"%s\") for atom %d.\n",p,i+1);
									return false;
								}
								break;

							case LAMMPS_FIELD_ELEMENT:
								if (needinfo) {
									try { tlabel = new char[strlen(p)+1]; } catch(...) { tlabel = NULL; }
									if (tlabel == NULL) NewException((double)(strlen(p)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
									strcpy(tlabel,p);
								}
								break;

							case LAMMPS_FIELD_XU:
								if (!g_bLAMMPSWrappedCoords) {
									#ifdef FAST_ATOF
										x = fast_atof(p) * 100.0; // Angstrom --> pm
									#else
										x = atof(p) * 100.0; // Angstrom --> pm
									#endif
								}
								break;

							case LAMMPS_FIELD_YU:
								if (!g_bLAMMPSWrappedCoords) {
									#ifdef FAST_ATOF
										y = fast_atof(p) * 100.0; // Angstrom --> pm
									#else
										y = atof(p) * 100.0; // Angstrom --> pm
									#endif
								}
								break;

							case LAMMPS_FIELD_ZU:
								if (!g_bLAMMPSWrappedCoords) {
									#ifdef FAST_ATOF
										z = fast_atof(p) * 100.0; // Angstrom --> pm
									#else
										z = atof(p) * 100.0; // Angstrom --> pm
									#endif
								}
								break;

							case LAMMPS_FIELD_X:
								if (g_bLAMMPSWrappedCoords) {
									#ifdef FAST_ATOF
										x = fast_atof(p) * 100.0; // Angstrom --> pm
									#else
										x = atof(p) * 100.0; // Angstrom --> pm
									#endif
								}
								break;

							case LAMMPS_FIELD_Y:
								if (g_bLAMMPSWrappedCoords) {
									#ifdef FAST_ATOF
										y = fast_atof(p) * 100.0; // Angstrom --> pm
									#else
										y = atof(p) * 100.0; // Angstrom --> pm
									#endif
								}
								break;

							case LAMMPS_FIELD_Z:
								if (g_bLAMMPSWrappedCoords) {
									#ifdef FAST_ATOF
										z = fast_atof(p) * 100.0; // Angstrom --> pm
									#else
										z = atof(p) * 100.0; // Angstrom --> pm
									#endif
								}
								break;

							case LAMMPS_FIELD_VX:
								if (g_bLAMMPSVelocity) {
									#ifdef FAST_ATOF
										vx = fast_atof(p) * 1.0e5; // Angstrom/fs --> m/s
									#else
										vx = atof(p) * 1.0e5; // Angstrom/fs --> m/s
									#endif
								}
								break;

							case LAMMPS_FIELD_VY:
								if (g_bLAMMPSVelocity) {
									#ifdef FAST_ATOF
										vy = fast_atof(p) * 1.0e5; // Angstrom/fs --> m/s
									#else
										vy = atof(p) * 1.0e5; // Angstrom/fs --> m/s
									#endif
								}
								break;

							case LAMMPS_FIELD_VZ:
								if (g_bLAMMPSVelocity) {
									#ifdef FAST_ATOF
										vz = fast_atof(p) * 1.0e5; // Angstrom/fs --> m/s
									#else
										vz = atof(p) * 1.0e5; // Angstrom/fs --> m/s
									#endif
								}
								break;

							case LAMMPS_FIELD_FX:
								if (g_bLAMMPSForce) {
									#ifdef FAST_ATOF
										fx = fast_atof(p);
									#else
										fx = atof(p);
									#endif
								}
								break;

							case LAMMPS_FIELD_FY:
								if (g_bLAMMPSForce) {
									#ifdef FAST_ATOF
										fy = fast_atof(p);
									#else
										fy = atof(p);
									#endif
								}
								break;

							case LAMMPS_FIELD_FZ:
								if (g_bLAMMPSForce) {
									#ifdef FAST_ATOF
										fz = fast_atof(p);
									#else
										fz = atof(p);
									#endif
								}
								break;

							case LAMMPS_FIELD_Q:
								if (g_bReadLAMMPSCharges) {
									#ifdef FAST_ATOF
										fq = fast_atof(p);
									#else
										fq = atof(p);
									#endif
								}
								break;

						}
					}

					ti++;

					if (last)
						break;
					p = q+1;
				}

				if (ti != g_iLAMMPSColumnCount) {
					eprintf("\nCTimeStep::ReadLAMMPS(): %d: Incomplete line - found only %d / %d columns.\n",i+1,ti,g_iLAMMPSColumnCount);
					mprintf("  \"%s\"\n",buf);
					return false;
				}

				if (g_iaLAMMPSColumnToOffset[ LAMMPS_FIELD_ID ] == -1)
					id = i;

				m_vaCoords[id] = CxDVector3(x,y,z);

				if (g_bKeepOriginalCoords)
					m_vaCoords_Original[id] = CxDVector3(x,y,z);

				if (needinfo)
					m_paLabels[id] = tlabel;

				if (g_bLAMMPSVelocity)
					m_vaVelocities[id] = CxDVector3(vx,vy,vz);

				if (g_bLAMMPSForce)
					m_vaForces[id] = CxDVector3(fx,fy,fz);

				if (g_bReadLAMMPSCharges)
					m_faCharge[id] = fq;

			}
			break;
		}
	}
	return true;
}


bool CTimeStep::SkipLAMMPS(FILE *a)
{
	int i, j;
	char buf[256];

	j = 0;
	while (true)
	{
		(void)fgets_bin(buf,256,a);
		if (feof(a))
			return false;
		if (strlen(buf) == 0)
			continue;
		buf[strlen(buf)-1] = 0;
		if (strstr(buf,"ITEM: NUMBER OF ATOMS") != 0)
		{
			(void)fgets_bin(buf,256,a);
			if (feof(a))
				return false;
			buf[strlen(buf)-1] = 0;
			j = atoi(buf);
			continue;
		}
		if (strstr(buf,"ITEM: ATOMS") != 0)
		{
			if (j == 0)
			{
				eprintf("CTimeStep::SkipLAMMPS():  \"ITEM: ATOMS\" before \"ITEM: NUMBER OF ATOMS\".\n");
				return false;
			}
			for (i=0;i<j;i++)
			{
				(void)fgets_bin(buf,256,a);
				if (feof(a))
					return false;
			}
			break;
		}
	}
	return true;
}


bool CTimeStep::ReadDLPOLY(FILE *a, bool needinfo)
{
	int i;
	char buf[256], obuf[256], *p, *q, *r;
	double x, y, z;

	m_vaCoords.RemoveAll_KeepSize();
	for (i=0;i<m_paLabels.GetSize();i++)
		delete[] (char*)m_paLabels[i];
	m_paLabels.RemoveAll();
	m_iGesAtomCount = 0;
	if (g_bKeepOriginalCoords)
		m_vaCoords_Original.RemoveAll_KeepSize();
	while (true)
	{
		(void)fgets_bin(buf,256,a);
		if (feof(a))
			return false;
		if (strlen(buf) == 0)
			continue;
		buf[strlen(buf)-1] = 0;
		strcpy(obuf,buf);

		if (strstr(buf,"timestep") != 0)
		{
			p = &buf[0];
			while (*p == ' ')
				p++;
			q = p;
			while ((*q != ' ') && (*q != 0))
				q++;
			if (*q == 0)
			{
				eprintf("\nCTimeStep::ReadDLPOLY(): Incomplete line A: \"%s\"\n",obuf);
				return false;
			}
			*q = 0;
			p=q+1;

			while (*p == ' ')
				p++;
			q = p;
			while ((*q != ' ') && (*q != 0))
				q++;
			if (*q == 0)
			{
				eprintf("\nCTimeStep::ReadDLPOLY(): Incomplete line B: \"%s\"\n",obuf);
				return false;
			}
			*q = 0;
			p=q+1;

			while (*p == ' ')
				p++;
			q = p;
			while ((*q != ' ') && (*q != 0))
				q++;
			if (*q == 0)
			{
				eprintf("\nCTimeStep::ReadDLPOLY(): Incomplete line C: \"%s\"\n",obuf);
				return false;
			}
			*q = 0;
			m_iGesAtomCount = atoi(p);

			if (g_bNPT)
			{
				for (i=0;i<3;i++)
				{
					(void)fgets_bin(buf,256,a);
					if (feof(a))
						return false;

					p = &buf[0];

					while (*p == ' ')
						p++;
					q = p;
					while ((*q != ' ') && (*q != 0))
						q++;
					if (*q == 0)
					{
						eprintf("\nCTimeStep::ReadDLPOLY(): Incomplete line D: \"%s\"\n",obuf);
						return false;
					}
					*q = 0;
					x = atof(p);
					p = q+1;

					while (*p == ' ')
						p++;
					q = p;
					while ((*q != ' ') && (*q != 0))
						q++;
					if (*q == 0)
					{
						eprintf("\nCTimeStep::ReadDLPOLY(): Incomplete line E: \"%s\"\n",obuf);
						return false;
					}
					*q = 0;
					y = atof(p);
					p = q+1;

					while (*p == ' ')
						p++;
					q = p;
					while ((*q != ' ') && (*q != 0))
						q++;
					*q = 0;
					z = atof(p);

					switch(i)
					{
						case 0:
							if ((y != 0) || (z != 0))
							{
								eprintf("\nCTimeStep::ReadDLPOLY(): X: Only orthorhombic cells are supported.\n");
								return false;
							}
							g_fBoxX = x*100.0;
							break;

						case 1:
							if ((x != 0) || (z != 0))
							{
								eprintf("\nCTimeStep::ReadDLPOLY(): Y: Only orthorhombic cells are supported.\n");
								return false;
							}
							g_fBoxY = y*100.0;
							break;

						case 2:
							if ((x != 0) || (y != 0))
							{
								eprintf("\nCTimeStep::ReadDLPOLY(): Z: Only orthorhombic cells are supported.\n");
								return false;
							}
							g_fBoxZ = z*100.0;
							break;

					}
				}

				if (g_bDoubleBox) {
					g_fBoxX *= g_iDoubleBoxX;
					g_fBoxY *= g_iDoubleBoxY;
					g_fBoxZ *= g_iDoubleBoxZ;
				}

			} else
			{
				(void)fgets_bin(buf,256,a);
				(void)fgets_bin(buf,256,a);
				(void)fgets_bin(buf,256,a);
			}

			if (m_iGesAtomCount == 0)
			{
				eprintf("CTimeStep::ReadDLPOLY(): Error: Atom count is 0.\n");
				return false;
			}

			if (needinfo)
			{
				if (g_bDoubleBox)
					m_paLabels.SetSize(m_iGesAtomCount*g_iDoubleBoxFactor);
						else m_paLabels.SetSize(m_iGesAtomCount);
			}

			for (i=0;i<(int)m_iGesAtomCount;i++)
			{
_readagain:
				(void)fgets_bin(buf,256,a);
				if (feof(a))
					return false;

				p = &buf[0];
				while (*p == ' ')
					p++;
				q = p;

				// Numer at beginning of line: Likely velocities. Skip that line
				if ((*p == '-') || (*p == '.') || ((*p >= '0') && (*p <= '9')))
					goto _readagain;

				while ((*q != ' ') && (*q != 0))
					q++;
				if (*q == 0)
				{
					eprintf("\nCTimeStep::ReadDLPOLY(): %d: Incomplete line F: \"%s\"\n",i+1,obuf);
					return false;
				}
				*q = 0;

				if (needinfo)
				{
		//			if (strlen(p) > 7)
		//			{
		//				eprintf("\nCTimeStep::ReadDLPOLY(): \"%s\" - maximum length for atom labels is 7 chars; truncating.\n",p);
		//				p[7] = 0;
		//			}

					try { r = new char[strlen(p)+1]; } catch(...) { r = NULL; }
					if (r == NULL) NewException((double)(strlen(p)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					strcpy(r,p);
					m_paLabels[i] = r;
				}

				(void)fgets_bin(buf,256,a);
				if (feof(a))
					return false;

				p = &buf[0];

				while (*p == ' ')
					p++;
				q = p;
				while ((*q != ' ') && (*q != 0))
					q++;
				if (*q == 0)
				{
					eprintf("\nCTimeStep::ReadDLPOLY(): %d: Incomplete line G: \"%s\"\n",i+1,obuf);
					return false;
				}
				*q = 0;
				//x = atof(p)*100.0;
				#ifdef FAST_ATOF
					x = fast_atof(p)*100.0;
				#else
					x = atof(p)*100.0;
				#endif
				p = q+1;

				while (*p == ' ')
					p++;
				q = p;
				while ((*q != ' ') && (*q != 0))
					q++;
				if (*q == 0)
				{
					eprintf("\nCTimeStep::ReadDLPOLY(): %d: Incomplete line H: \"%s\"\n",i+1,obuf);
					return false;
				}
				*q = 0;
				//y = atof(p)*100.0;
				#ifdef FAST_ATOF
					y = fast_atof(p)*100.0;
				#else
					y = atof(p)*100.0;
				#endif
				p = q+1;

				while (*p == ' ')
					p++;
				q = p;
				while ((*q != ' ') && (*q != 0))
					q++;
				*q = 0;
				//z = atof(p)*100.0;
				#ifdef FAST_ATOF
					z = fast_atof(p)*100.0;
				#else
					z = atof(p)*100.0;
				#endif

				m_vaCoords.Add(CxDVector3(x,y,z));

				if (g_bKeepOriginalCoords)
					m_vaCoords_Original.Add(CxDVector3(x,y,z));
			}
			break;
		}
	}
	return true;
}


bool CTimeStep::SkipDLPOLY(FILE *a)
{
	int i, j;
	char buf[256], obuf[256], *p, *q;

	while (true)
	{
		(void)fgets_bin(buf,256,a);
		if (feof(a))
			return false;
		if (strlen(buf) == 0)
			continue;
		buf[strlen(buf)-1] = 0;
		strcpy(obuf,buf);

		if (strstr(buf,"timestep") != 0)
		{
			p = &buf[0];
			while (*p == ' ')
				p++;
			q = p;
			while ((*q != ' ') && (*q != 0))
				q++;
			if (*q == 0)
			{
				eprintf("\nCTimeStep::ReadDLPOLY(): Incomplete line A: \"%s\"\n",obuf);
				return false;
			}
			*q = 0;
			p=q+1;

			while (*p == ' ')
				p++;
			q = p;
			while ((*q != ' ') && (*q != 0))
				q++;
			if (*q == 0)
			{
				eprintf("\nCTimeStep::ReadDLPOLY(): Incomplete line B: \"%s\"\n",obuf);
				return false;
			}
			*q = 0;
			p=q+1;

			while (*p == ' ')
				p++;
			q = p;
			while ((*q != ' ') && (*q != 0))
				q++;
			if (*q == 0)
			{
				eprintf("\nCTimeStep::ReadDLPOLY(): Incomplete line C: \"%s\"\n",obuf);
				return false;
			}
			*q = 0;
			j = atoi(p);

			(void)fgets_bin(buf,256,a);
			(void)fgets_bin(buf,256,a);
			(void)fgets_bin(buf,256,a);

			if (j == 0)
			{
				eprintf("CTimeStep::ReadDLPOLY(): Error: Atom count is 0.\n");
				return false;
			}

			for (i=0;i<j;i++)
			{
				(void)fgets_bin(buf,256,a);
				(void)fgets_bin(buf,256,a);
				if (feof(a))
					return false;
			}
			break;
		}
	}
	return true;
}


bool CTimeStep::ReadAmber(FILE *a, bool needinfo)
{
	FILE *b;
	int z;
	char buf[1024], *p, *q, *r;
	double tf;

	if (needinfo)
	{
		b = fopen(g_sAmberParmFile,"rb");
		if (b == NULL)
		{
			eprintf("ReadAmber(): Error opening file \"%s\".\n",(const char*)g_sAmberParmFile);
			return false;
		}
		do {
			(void)fgets_bin(buf,1024,b);
			if (feof(b))
				break;
		} while (strstr(buf,"%FLAG ATOM_NAME") == NULL);
		if (feof(b))
		{
			eprintf("ReadAmber(): Did not find \"%cFLAG ATOM_NAME\" in .prmtop file.\n",'%');
			fclose(b);
			return false;
		}
		(void)fgets_bin(buf,1024,b); // "%FORMAT" line

		m_paLabels.RemoveAll_KeepSize();
		while (true)
		{
			(void)fgets_bin(buf,1024,b);
			buf[strlen(buf)-1] = 0;
			if (strchr(buf,'%') != NULL)
				break;
			if (feof(b))
			{
				eprintf("ReadAmber(): Unexpected end of file while reading .prmtop file.\n");
				fclose(b);
				return false;
			}
			p = buf;

			while (true)
			{
				while ((*p == ' ') || ((*p >= '0') && (*p <= '9')))
					p++;
				if (*p == 0)
					break;
				q = p;
				while ((*p != ' ') && (*p != 0) && !((*p >= '0') && (*p <= '9')))
					p++;
				if (*p == 0)
				{
					r = p;
					while (isdigit(*(r-1)))
						r--;
					*r = 0;

					r = new char[strlen(q)+1];
					strcpy(r,q);
					m_paLabels.Add(r);
					break;
				} else
				{
					r = p;
					while (isdigit(*(r-1)))
						r--;
					*r = 0;

					*p = 0;
					p++;
					r = new char[strlen(q)+1];
					strcpy(r,q);
					m_paLabels.Add(r);
				}
			}
		}
		fclose(b);
	/*	mprintf("Found %d atoms:\n",m_paLabels.GetSize());
		for (z=0;z<m_paLabels.GetSize();z++)
			mprintf("  %4d: \"%s\".\n",z,m_paLabels[z]);*/
		m_iGesAtomCount = m_paLabels.GetSize();
	} else
		m_iGesAtomCount = g_iGesAtomCount;

	//mprintf("@%d\n",m_iGesAtomCount);

	m_vaCoords.RemoveAll_KeepSize();

#ifdef USE_NETCDF
	// START new Amber Version
// Reopening File and Sanity Checks:
	int ncid;
	nc_open(g_sInputTraj, NC_NOWRITE, &ncid); 

// Reading NPT-Data:
	size_t start[3] = {(size_t)g_iSteps, 0, 0};
	size_t count[3] = {1, 3, 3};

	std::vector<float> xyz;
	if (g_bNPT)
	{
		xyz.resize(6);
	
		nc_get_vara_float(ncid, 5, start, count, &xyz[0]);	// cell length
		nc_get_vara_float(ncid, 6, start, count, &xyz[3]);	// cell vector
	
		g_fBoxX = xyz[0]*100;
		g_fBoxY = xyz[1]*100;
		g_fBoxZ = xyz[2]*100;
		g_fBoxMinDiam = MIN3(g_fBoxX,g_fBoxY,g_fBoxZ);	
		g_fBoxAngleA = xyz[3];
		g_fBoxAngleB = xyz[4];
		g_fBoxAngleC = xyz[5];
	}

// Reading Coordinates:
	count[1] = m_iGesAtomCount;
	long unsigned int n;

	xyz.clear();
	xyz.resize(m_iGesAtomCount*3);
	nc_get_vara_float(ncid, 2, start, count, &xyz[0]);

	for ( n=0; n<m_iGesAtomCount; n++ )
		m_vaCoords.Add(CxDVector3(xyz[n*3]*100,xyz[n*3+1]*100,xyz[n*3+2]*100));

	return true;
	// END new Amber Version
#endif

	z = 0;
	if (ftell(a) == 0)
		(void)fgets_bin(buf,1024,a); // Comment line only in the first time step
	while (true)
	{
		(void)fgets_bin(buf,1024,a);
		buf[strlen(buf)-1] = 0;
		if (feof(a))
		{
			eprintf("\nEnd of trajectory file reached.\n");
			return false;
		}
		p = buf;

		while (true)
		{
			while (*p == ' ')
				p++;
			if (*p == 0)
				break;
			q = p;
			while ((*p != ' ') && (*p != 0))
				p++;
			if (*p == 0)
			{
				//tf = atof(q)*100.0;
				#ifdef FAST_ATOF
					tf = fast_atof(q)*100.0;
				#else
					tf = atof(q)*100.0;
				#endif
				switch(z)
				{
					case 0: m_vaCoords.Add(CxDVector3(tf,0,0)); break;            // X
					case 1: m_vaCoords[m_vaCoords.GetSize()-1][1] = tf; break;    // Y
					case 2: m_vaCoords[m_vaCoords.GetSize()-1][2] = tf; break;    // Z
				}
				z++;
				if (z > 2)
					z = 0;
				break;
			} else
			{
				*p = 0;
				p++;
				//tf = atof(q)*100.0;
				#ifdef FAST_ATOF
					tf = fast_atof(q)*100.0;
				#else
					tf = atof(q)*100.0;
				#endif
				switch(z)
				{
					case 0: m_vaCoords.Add(CxDVector3(tf,0,0)); break;            // X
					case 1: m_vaCoords[m_vaCoords.GetSize()-1][1] = tf; break;    // Y
					case 2: m_vaCoords[m_vaCoords.GetSize()-1][2] = tf; break;    // Z
				}
				z++;
				if (z > 2)
					z = 0;
			}
			if ((m_vaCoords.GetSize() >= (int)m_iGesAtomCount) && (z == 0))
				break;
		}
		//mprintf("@ %d/%d ...\n",m_vaCoords.GetSize(),m_iGesAtomCount);
		if ((m_vaCoords.GetSize() >= (int)m_iGesAtomCount) && (z == 0))
			break;
	}
/*	mprintf("Found %d atoms:\n",m_paLabels.GetSize());
	for (z=0;z<m_paLabels.GetSize();z++)
		mprintf("  %4d: \"%s\"  %8f  %8f  %8f.\n",z,m_paLabels[z],m_vaCoords[z][0],m_vaCoords[z][1],m_vaCoords[z][2]);*/

	if (!g_bAmberCellAsked) {

		mprintf("\n");
		g_bAmberCellInfo = AskYesNo("    Does your Amber .mdcrd file contain cell information (y/n)? [no] ",false);
	}

	if (g_bAmberCellInfo) {

		(void)fgets_bin(buf,1024,a); // Box size
		buf[strlen(buf)-1] = 0;
		if (g_bNPT)
		{
			p = buf;

			while (*p == ' ')
				p++;
			q = p;
			while ((*p != ' ') && (*p != 0))
				p++;
			if (*p == 0)
			{
				eprintf("ReadAmber(): Unexpected end of line while reading box length.\n");
				return false;
			}
			*p = 0;
			g_fBoxX = atof(q)*100.0;
			p++;

			while (*p == ' ')
				p++;
			q = p;
			while ((*p != ' ') && (*p != 0))
				p++;
			if (*p == 0)
			{
				eprintf("ReadAmber(): Unexpected end of line while reading box length.\n");
				return false;
			}
			*p = 0;
			g_fBoxY = atof(q)*100.0;
			p++;

			while (*p == ' ')
				p++;
			q = p;
			while ((*p != ' ') && (*p != 0))
				p++;
			*p = 0;
			g_fBoxZ = atof(q)*100.0;
		}
	}

	if (g_bAmberCellInfo && !g_bAmberCellAsked) {
		mprintf("\n");
		mprintf("    Assuming that the following numbers are the cell size:\n");
		mprintf("      X=%12.6f pm,  Y=%12.6f pm,  Z=%12.6f pm\n\n",g_fBoxX,g_fBoxY,g_fBoxZ);
		if (!AskYesNo("    Are this the correct simulation cell dimensions (y/n)? [yes] ",true)) {
			mprintf("\n");
			eprintf("    Please restart the analysis and specify that your trajectory does not contain cell information!\n\n");
			abort();
		}
		g_bAmberCellAsked = true;
	}

	if (feof(a))
	{
		mprintf("\nReadAmber(): End of trajectory reached.\n");
		return false;
	}

	return true;
}


bool CTimeStep::SkipAmber(FILE *a)
{
#ifdef USE_NETCDF
	return true;
#endif
	int z, i;
	char buf[1024], *p;

	m_iGesAtomCount = g_iGesAtomCount;

	if (ftell(a) == 0)
		(void)fgets_bin(buf,1024,a); // Comment line only in the first time step
	z = 0;
	i = 0;
	while (true)
	{
		(void)fgets_bin(buf,1024,a);
		buf[strlen(buf)-1] = 0;
		if (feof(a))
		{
			eprintf("ReadAmber(): Unexpected end of file while reading .mdcrd file.\n");
			return false;
		}
		p = buf;

		while (true)
		{
			while (*p == ' ')
				p++;
			if (*p == 0)
				break;
			while ((*p != ' ') && (*p != 0))
				p++;
			if (*p == 0)
			{
				if (z == 0)
					i++;
				z++;
				if (z > 2)
					z = 0;
				break;
			} else
			{
				*p = 0;
				p++;
				if (z == 0)
					i++;
				z++;
				if (z > 2)
					z = 0;
			}
		}
		if ((i == (int)m_iGesAtomCount) && (z == 0))
			break;
	}
	(void)fgets_bin(buf,1024,a); // Box size

	if (feof(a))
	{
		mprintf("ReadAmber(): End of trajectory reached.\n");
		return false;
	}

	return true;
}


bool CTimeStep::ReadMol2(FILE *a, bool needinfo)
{
	int i;
	char buf[256], obuf[256], *p, *q, *r;
	double x, y, z;

	(void)fgets_bin(buf,256,a);
	(void)fgets_bin(buf,256,a);
	(void)fgets_bin(buf,256,a);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	*q = 0;
	m_iGesAtomCount = atoi(p);
	
	if (needinfo)
	{
		m_vaCoords.RemoveAll_KeepSize();
		for (i=0;i<m_paLabels.GetSize();i++)
			delete[] (char*)m_paLabels[i];
		m_paLabels.RemoveAll();
		m_paLabels.SetSize(m_iGesAtomCount);
		for (i=0;i<m_paMol2Types.GetSize();i++)
			delete[] (char*)m_paMol2Types[i];
		m_paMol2Types.RemoveAll();
		m_paMol2Types.SetSize(m_iGesAtomCount);
	}

/*	fgets_bin(buf,256,a);
	fgets_bin(buf,256,a);
	fgets_bin(buf,256,a);
	fgets_bin(buf,256,a);
	fgets_bin(buf,256,a);*/

	do {
		(void)fgets_bin(buf,256,a);
		if (feof(a))
		{
			eprintf("ReadMol2(): Error: Unexpected end of file.\n");
			return false;
		}
	} while (strstr(buf,"ATOM") == NULL);

	for (i=0;i<(long)m_iGesAtomCount;i++)
	{
		(void)fgets_bin(buf,256,a);
		if (feof(a))
			return false;
		buf[strlen(buf)-1] = 0;
		strcpy(obuf,buf);
//		mprintf("  \"%s\"\n",buf);
		p = buf;
		while (*p == ' ')
			p++;
		while ((*p != ' ') && (*p != 0))
			p++;
		while (*p == ' ')
			p++;
		q = p;
		while (isalpha(*q) || (*q == '_'))
			q++;
		if (needinfo)
		{
			*q = 0;
	//		if (strlen(p) > 7)
	//		{
	//			eprintf("\nCTimeStep::ReadMol2(): \"%s\" - maximum length for atom labels is 7 chars; truncating.\n",p);
	//			p[7] = 0;
	//		}

			try { r = new char[strlen(p)+1]; } catch(...) { r = NULL; }
			if (r == NULL) NewException((double)(strlen(p)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(r,p);
			m_paLabels[i] = r;
		}
		p = q+1;
		while ((*p != ' ') && (*p != 0))
			p++;
		while (*p == ' ')
			p++;

		q = strchr(p,' ');
		if (q == NULL)
		{
			eprintf("\nCTimeStep::ReadMol2(): Error 1. %d. \"%s\".",i+1,obuf);
			return false;
		}
		*q = 0;
		x = atof(p) * 100.0;
		p = q+1;
		while (*p == ' ')
			p++;

		q = strchr(p,' ');
		if (q == NULL)
		{
			eprintf("\nCTimeStep::ReadMol2(): Error 2. %d. \"%s\".",i+1,obuf);
			return false;
		}
		*q = 0;
		y = atof(p) * 100.0;
		p = q+1;
		while (*p == ' ')
			p++;

		q = strchr(p,' ');
		if (q == NULL)
		{
			eprintf("\nCTimeStep::ReadMol2(): Error 3. %d. \"%s\".",i+1,obuf);
			return false;
		}
		*q = 0;
		z = atof(p) * 100.0;
		m_vaCoords.Add(CxDVector3(x,y,z));

		p = q+1;
		while (*p == ' ')
			p++;
		q = p;
		while ((*q != ' ') && (*q != 0))
			q++;
		if (needinfo)
		{
			*q = 0;

			try { r = new char[strlen(p)+1]; } catch(...) { r = NULL; }
			if (r == NULL) NewException((double)(strlen(p)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(r,p);
			m_paMol2Types[i] = r;
		}
/*		mprintf("  X=%f, Y=%f, Z=%f",x,y,z);
		if (needinfo)
			mprintf(", A=%s, B=%s",(char*)m_paLabels[i],(char*)m_paMol2Types[i]);
		mprintf("\n");*/
	}
	m_iGesAtomCount = m_vaCoords.GetSize();

	return true;
}


bool CTimeStep::SkipMol2(FILE *a)
{
	int i;
	char buf[256];

	(void)fgets_bin(buf,256,a);
	(void)fgets_bin(buf,256,a);
	(void)fgets_bin(buf,256,a);
	(void)fgets_bin(buf,256,a);
	(void)fgets_bin(buf,256,a);
	(void)fgets_bin(buf,256,a);
	(void)fgets_bin(buf,256,a);
	(void)fgets_bin(buf,256,a);

	for (i=0;i<(long)m_iGesAtomCount;i++)
	{
		(void)fgets_bin(buf,256,a);
		if (feof(a))
			return false;
	}

	return true;
}


void CTimeStep::CalcMinMax()
{
	int z;

	m_vMin[0] = 1E20;
	m_vMin[1] = 1E20;
	m_vMin[2] = 1E20;
	m_vMax[0] = -1E20;
	m_vMax[1] = -1E20;
	m_vMax[2] = -1E20;

	for (z=0;z<(long)m_iGesAtomCount;z++)
	{
		if (m_vMin[0] > m_vaCoords[z][0])
			m_vMin[0] = m_vaCoords[z][0];
		if (m_vMin[1] > m_vaCoords[z][1])
			m_vMin[1] = m_vaCoords[z][1];
		if (m_vMin[2] > m_vaCoords[z][2])
			m_vMin[2] = m_vaCoords[z][2];
		if (m_vMax[0] < m_vaCoords[z][0])
			m_vMax[0] = m_vaCoords[z][0];
		if (m_vMax[1] < m_vaCoords[z][1])
			m_vMax[1] = m_vaCoords[z][1];
		if (m_vMax[2] < m_vaCoords[z][2])
			m_vMax[2] = m_vaCoords[z][2];
	}
}


void CTimeStep::WriteMol2(FILE *a)
{
	int z, z2, z3, z4, ti, c, mc;
	CMolecule *m;
	CSingleMolecule *sm;
	int bonds;
	int *tpi;
	const char *cc;

	bonds = 0;
	for (z=0;z<g_oaSingleMolecules.GetSize();z++)
		bonds += ((CSingleMolecule*)g_oaSingleMolecules[z])->m_oaBonds.GetSize();

	try { tpi = new int[g_iGesAtomCount]; } catch(...) { tpi = NULL; }
	if (tpi == NULL) NewException((double)g_iGesAtomCount*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	mfprintf(a," @<TRIPOS>MOLECULE\n");
	mfprintf(a,"MOL\n");
	mfprintf(a,"    %d    %d     %d  0  0\n",g_iGesAtomCount,bonds,g_oaSingleMolecules.GetSize());
	mfprintf(a," SMALL\n");
	mfprintf(a,"resp\n\n\n");
	mfprintf(a," @<TRIPOS>ATOM\n");
	c = 0;
	mc = 0;
	if (m_faCharge.GetSize() != g_iGesAtomCount)
	{
		m_faCharge.SetSize(g_iGesAtomCount);
		for (z=0;z<g_iGesAtomCount;z++)
			m_faCharge[z] = 0;
	}
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			mc++;
			for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
			{
				if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
					continue;
				for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
				{
					ti = ((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4);
					if (m_paMol2Types.GetSize() != 0)
						cc = (char*)m_paMol2Types[ti];
							else cc = (const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName;
					tpi[ti] = c++;
					mfprintf(a,"  %6d  %2s  % 11.4f  % 11.4f  % 11.4f  %2s  %4d  MOL  % 8.4f\n",c,(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,m_vaCoords[ti][0]/100.0,m_vaCoords[ti][1]/100.0,m_vaCoords[ti][2]/100.0,cc,mc,m_faCharge[ti]);
				}
			}
		}
	}
	mfprintf(a," @<TRIPOS>BOND\n");
	c = 0;
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			for (z3=0;z3<sm->m_oaBonds.GetSize();z3++)
				mfprintf(a,"  %6d  %6d  %6d  1\n",++c,tpi[((CMolBond*)sm->m_oaBonds[z3])->m_iAtomOffset[0]]+1,tpi[((CMolBond*)sm->m_oaBonds[z3])->m_iAtomOffset[1]]+1);
		}
	}
	mfprintf(a," @<TRIPOS>SUBSTRUCTURE\n");
	c = 0;
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			mfprintf(a,"  %4d  MOL  %4d  TEMP              0 ****  ****    0 ROOT\n",c+1,c+1);
			c++;
		}
	}
	delete[] tpi;
}


bool CTimeStep::ReadCellVector(FILE *a)
{
	char buf[256], obuf[256];
	char *p, *q;
	int i;
	double tfa[9];
	bool last;
	CxDVector3 veca, vecb, vecc;


	(void)!fgets(buf,256,a);
	if (feof(a)) {
		eprintf("\nCTimeStep::ReadCellVector(): End of file.\n");
		eprintf("Your cell vector text file is too short.\n");
		return false;
	}
	if (strlen(buf) == 0) {
		eprintf("\nCTimeStep::ReadCellVector(): Empty line.\n");
		return false;
	}
	buf[strlen(buf)-1] = 0;
	strcpy(obuf,buf);
//	mprintf(GREY,"\nReadCellVector: \"%s\".\n",buf);
	p = buf;

	i = 0;
	last = false;
	while (true) {

		while ((!isnumeric(*p)) && (*p != 0))
			p++;
		if (*p == 0)
			break;
		q = p;
		while (isnumeric(*q))
			q++;
		if (*q == 0)
			last = true;
		*q = 0;
		tfa[i] = atof(p) * 100.0;

		i++;
		if (i == 9)
			break;
		if (last)
			break;
		p = q+1;
	}

	if (g_iCellVectorFileColumns != -1) {

		if (i != g_iCellVectorFileColumns) {
			eprintf("\nCTimeStep::ReadCellVector(): Error: Found %d numbers in line, expected %d numbers.\n",i,g_iCellVectorFileColumns);
			eprintf("Your cell vector file is damaged.\n");
			return false;
		}

	} else {

		if ((i != 3) && (i != 9)) {
			eprintf("\nCTimeStep::ReadCellVector(): Error: Found %d numbers in first line. Only 3 or 9 numbers are allowed.\n",i);
			return false;
		}
		g_iCellVectorFileColumns = i;
		mprintf("\n    Found %d numeric columns in first line of cell vector file.\n",g_iCellVectorFileColumns);
	}

	if (g_iCellVectorFileColumns == 3) {

		g_bFoundNonOrtho = false;
		g_bBoxNonOrtho = false;

		g_fBoxX = tfa[0];
		g_fBoxY = tfa[1];
		g_fBoxZ = tfa[2];

		if (g_bDoubleBox) {
			g_fBoxX *= g_iDoubleBoxX;
			g_fBoxY *= g_iDoubleBoxY;
			g_fBoxZ *= g_iDoubleBoxZ;
		}

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
			eprintf("CTimeStep::ReadCellVector(): Error: Encountered singular cell matrix (cell volume is zero).\n");
			abort();
		}

	} else if (g_iCellVectorFileColumns == 9) {

		g_bFoundNonOrtho = true;
		g_bBoxNonOrtho = true;

		g_mBoxFromOrtho(0,0) = tfa[0];
		g_mBoxFromOrtho(0,1) = tfa[1];
		g_mBoxFromOrtho(0,2) = tfa[2];

		g_mBoxFromOrtho(1,0) = tfa[3];
		g_mBoxFromOrtho(1,1) = tfa[4];
		g_mBoxFromOrtho(1,2) = tfa[5];

		g_mBoxFromOrtho(2,0) = tfa[6];
		g_mBoxFromOrtho(2,1) = tfa[7];
		g_mBoxFromOrtho(2,2) = tfa[8];

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

		// Backward transformation matrix
		g_mBoxToOrtho = CxDMatrix3(g_mBoxFromOrtho);
		if (!g_mBoxToOrtho.Invert()) {
			eprintf("CTimeStep::ReadCellVector(): Error: Encountered singular cell matrix (cell volume is zero).\n");
			abort();
		}

		veca = CxDVector3(g_mBoxFromOrtho(0,0),g_mBoxFromOrtho(0,1),g_mBoxFromOrtho(0,2));
		vecb = CxDVector3(g_mBoxFromOrtho(1,0),g_mBoxFromOrtho(1,1),g_mBoxFromOrtho(1,2));
		vecc = CxDVector3(g_mBoxFromOrtho(2,0),g_mBoxFromOrtho(2,1),g_mBoxFromOrtho(2,2));

		// Cell angles
		g_fBoxAngleA = acos(DotP(vecb,vecc) / vecb.GetLength() / vecc.GetLength()) * 180.0 / Pi;
		g_fBoxAngleB = acos(DotP(veca,vecc) / veca.GetLength() / vecc.GetLength()) * 180.0 / Pi;
		g_fBoxAngleC = acos(DotP(veca,vecb) / veca.GetLength() / vecb.GetLength()) * 180.0 / Pi;

		// Orthogonal bounding box
		g_fBoxX = fabs(g_mBoxFromOrtho(0,0)) + fabs(g_mBoxFromOrtho(1,0)) + fabs(g_mBoxFromOrtho(2,0));
		g_fBoxY = fabs(g_mBoxFromOrtho(0,1)) + fabs(g_mBoxFromOrtho(1,1)) + fabs(g_mBoxFromOrtho(2,1));
		g_fBoxZ = fabs(g_mBoxFromOrtho(0,2)) + fabs(g_mBoxFromOrtho(1,2)) + fabs(g_mBoxFromOrtho(2,2));

		// Minimal diameters
		g_fBoxMinDiamA = fabs(DotP(veca,Normalize(CrossP(vecb,vecc))));
		g_fBoxMinDiamB = fabs(DotP(vecb,Normalize(CrossP(veca,vecc))));
		g_fBoxMinDiamC = fabs(DotP(vecc,Normalize(CrossP(veca,vecb))));
		g_fBoxMinDiam = MIN3(g_fBoxMinDiamA,g_fBoxMinDiamB,g_fBoxMinDiamC);

	} else {

		eprintf("CTimeStep::ReadCellVector(): Internal error.\n");
		abort();
	}

	return true;
}


void CTimeStep::CenterCOM()
{
	CxDVector3 vc;
	double m;
	int z;

	if (g_iRemoveCOMFixAtom != -1) {

		CalcCenters();
		vc = m_vaCoords[g_iRemoveCOMFixAtom];

	} else {

		vc = CxDVector3(0,0,0);
		m = 0;
		for (z=0;z<g_iGesAtomCount;z++)
		{
			if (g_bWannier)
				if (g_baAtomIndex[z] == g_iWannierAtomType)
					continue;
			vc += m_vaCoords[z] * ((CAtom*)g_oaAtoms[g_baAtomIndex[z]])->m_pElement->m_fMass;
			m += ((CAtom*)g_oaAtoms[g_baAtomIndex[z]])->m_pElement->m_fMass;
		}
		vc /= m;

	}

	for (z=0;z<g_iGesAtomCount;z++)
		m_vaCoords[z] -= vc;
}


double CTimeStep::FoldedDistance(int i1, int i2)
{
	CxDVector3 t;

	t = m_vaCoords[i2]-m_vaCoords[i1];

	if (g_bPeriodicX)
	{
		while (t[0] < -g_fBoxX/2) t[0] += g_fBoxX;
		while (t[0] >= g_fBoxX/2) t[0] -= g_fBoxX;
	}

	if (g_bPeriodicY)
	{
		while (t[1] < -g_fBoxY/2) t[1] += g_fBoxY;
		while (t[1] >= g_fBoxY/2) t[1] -= g_fBoxY;
	}

	if (g_bPeriodicZ)
	{
		while (t[2] < -g_fBoxZ/2) t[2] += g_fBoxZ;
		while (t[2] >= g_fBoxZ/2) t[2] -= g_fBoxZ;
	}

	return t.GetLength();
}


void CTimeStep::FoldAtomsPositive()
{
	BTIN;
	int z, n;


	if (!g_bPeriodic)
		return;

	for (z=0;z<g_iGesVirtAtomCount;z++)
	{
		if (g_bBoxNonOrtho) {

			CxDVector3 coord = g_mBoxToOrtho * m_vaCoords[z];

			n = PeriodicImage1D( coord[0], 0, 1.0 );
			if (n != 0)
				coord[0] -= n;

			n = PeriodicImage1D( coord[1], 0, 1.0 );
			if (n != 0)
				coord[1] -= n;

			n = PeriodicImage1D( coord[2], 0, 1.0 );
			if (n != 0)
				coord[2] -= n;

/*			if (g_bPeriodicX) {
				while (coord[0] < 0.0) coord[0] += 1.0;
				while (coord[0] >= 1.0) coord[0] -= 1.0;
			}
			if (g_bPeriodicY) {
				while (coord[1] < 0.0) coord[1] += 1.0f;
				while (coord[1] >= 1.0) coord[1] -= 1.0f;
			}
			if (g_bPeriodicZ) {
				while (coord[2] < 0.0) coord[2] += 1.0;
				while (coord[2] >= 1.0) coord[2] -= 1.0;
			}*/

			m_vaCoords[z] = g_mBoxFromOrtho * coord;

		} else {

			if (g_bPeriodicX) {
				n = PeriodicImage1D( m_vaCoords[z][0], 0, g_fBoxX );
				if (n != 0)
					m_vaCoords[z][0] -= n*g_fBoxX;
			}

			if (g_bPeriodicY) {
				n = PeriodicImage1D( m_vaCoords[z][1], 0, g_fBoxY );
				if (n != 0)
					m_vaCoords[z][1] -= n*g_fBoxY;
			}

			if (g_bPeriodicZ) {
				n = PeriodicImage1D( m_vaCoords[z][2], 0, g_fBoxZ );
				if (n != 0)
					m_vaCoords[z][2] -= n*g_fBoxZ;
			}

/*			if (g_bPeriodicX)
			{
				while (m_vaCoords[z][0] < 0) m_vaCoords[z][0] += g_fBoxX;
				while (m_vaCoords[z][0] >= g_fBoxX) m_vaCoords[z][0] -= g_fBoxX;
			}

			if (g_bPeriodicX)
			{
				while (m_vaCoords[z][1] < 0) m_vaCoords[z][1] += g_fBoxY;
				while (m_vaCoords[z][1] >= g_fBoxY) m_vaCoords[z][1] -= g_fBoxY;
			}

			if (g_bPeriodicX)
			{
				while (m_vaCoords[z][2] < 0) m_vaCoords[z][2] += g_fBoxZ;
				while (m_vaCoords[z][2] >= g_fBoxZ) m_vaCoords[z][2] -= g_fBoxZ;
			}*/
		}
	}
	BTOUT; 
}


void CTimeStep::WritePOV(const char *s)
{
	FILE *b;
	CMolecule *m;
	CSingleMolecule *sm;
	CElement *el, *el2;
	int z, z2, z3, z4, o, o2;
	CxDVector3 vec1, vec2, vec3, vec1b, vec2b, vec3b, vecA, vecB, vecC, vecD;
	CMolBond *mb;
	CxDVector3 cam;
//	double cr, cg, cb;

	b = OpenFileWrite(s,true);

	mfprintf(b,"// Written by TRAVIS\n");
	mfprintf(b,"// See http://www.travis-analyzer.de\n\n");
	mfprintf(b,"#version 3.6;\n");
	mfprintf(b,"\n");

	mfprintf(b,"/**** Atoms ****/\n");

	mfprintf(b,"#declare atom_draw          = true;\n");
	mfprintf(b,"#declare atom_r             = 0.65;\n");
	mfprintf(b,"#declare atom_specular      = 0.7;\n");
	mfprintf(b,"#declare atom_reflect       = 0;\n");
	mfprintf(b,"#declare atom_ambient       = 0.2;\n");
	mfprintf(b,"#declare atom_diffuse       = 0.7;\n");
	mfprintf(b,"//#declare atom_color         = < 1.0, 1.0, 1.0, 0, 0 >;\n");
	mfprintf(b,"//#declare atom_trans         = 0.7;\n");

	mfprintf(b,"\n");

	mfprintf(b,"#declare atom_draw_halo1      = true;\n");
	mfprintf(b,"#declare atom_r_halo1         = 0.008;\n");
	mfprintf(b,"#declare atom_d_halo1         = 0.0125;\n");
	mfprintf(b,"#declare atom_color_halo1     = < 0, 0, 0, 0, 0 >;\n");
	mfprintf(b,"#declare atom_specular_halo1  = 0;\n");
	mfprintf(b,"#declare atom_reflect_halo1   = 0;\n");
	mfprintf(b,"#declare atom_ambient_halo1   = 1.0;\n");
	mfprintf(b,"#declare atom_diffuse_halo1   = 0;\n");

	mfprintf(b,"\n");

	mfprintf(b,"#declare atom_draw_halo2      = true;\n");
	mfprintf(b,"#declare atom_r_halo2         = 0.003;\n");
	mfprintf(b,"#declare atom_d_halo2         = 0.01875;\n");
	mfprintf(b,"#declare atom_color_halo2     = < 1, 1, 1, 0, 0 >;\n");
	mfprintf(b,"#declare atom_specular_halo2  = 0;\n");
	mfprintf(b,"#declare atom_reflect_halo2   = 0;\n");
	mfprintf(b,"#declare atom_ambient_halo2   = 1.0;\n");
	mfprintf(b,"#declare atom_diffuse_halo2   = 0;\n");

	mfprintf(b,"\n");

	mfprintf(b,"#declare atom_draw_halo3      = true;\n");
	mfprintf(b,"#declare atom_r_halo3         = 0.003;\n");
	mfprintf(b,"#declare atom_d_halo3         = 0.025;\n");
	mfprintf(b,"#declare atom_color_halo3     = < 1, 1, 1, 0, 0.5 >;\n");
	mfprintf(b,"#declare atom_specular_halo3  = 0;\n");
	mfprintf(b,"#declare atom_reflect_halo3   = 0;\n");
	mfprintf(b,"#declare atom_ambient_halo3   = 1.0;\n");
	mfprintf(b,"#declare atom_diffuse_halo3   = 0;\n");

	mfprintf(b,"\n");

	mfprintf(b,"/**** Bonds ****/\n");

	mfprintf(b,"#declare bond_draw        = true;\n");
	mfprintf(b,"#declare bond_r           = 0.015;\n");
	mfprintf(b,"#declare bond_specular    = 0.7;\n");
	mfprintf(b,"#declare bond_reflect     = 0;\n");
	mfprintf(b,"#declare bond_ambient     = 0.2;\n");
	mfprintf(b,"#declare bond_diffuse     = 0.7;\n");

	mfprintf(b,"\n");

	mfprintf(b,"#declare bond_draw_halo1      = true;\n");
	mfprintf(b,"#declare bond_r_halo1         = 0.008;\n");
	mfprintf(b,"#declare bond_d_halo1         = 0;\n");
	mfprintf(b,"#declare bond_color_halo1     = < 0, 0, 0, 0, 0 >;\n");
	mfprintf(b,"#declare bond_specular_halo1  = 0;\n");
	mfprintf(b,"#declare bond_reflect_halo1   = 0;\n");
	mfprintf(b,"#declare bond_ambient_halo1   = 1.0;\n");
	mfprintf(b,"#declare bond_diffuse_halo1   = 0;\n");

	mfprintf(b,"\n");

	mfprintf(b,"#declare bond_draw_halo2      = true;\n");
	mfprintf(b,"#declare bond_r_halo2         = 0.003;\n");
	mfprintf(b,"#declare bond_d_halo2         = 0.0286;\n");
	mfprintf(b,"#declare bond_color_halo2     = < 1, 1, 1, 0, 0 >;\n");
	mfprintf(b,"#declare bond_specular_halo2  = 0;\n");
	mfprintf(b,"#declare bond_reflect_halo2   = 0;\n");
	mfprintf(b,"#declare bond_ambient_halo2   = 1.0;\n");
	mfprintf(b,"#declare bond_diffuse_halo2   = 0;\n");

	mfprintf(b,"\n");

	mfprintf(b,"#declare bond_draw_halo3      = true;\n");
	mfprintf(b,"#declare bond_r_halo3         = 0.003;\n");
	mfprintf(b,"#declare bond_d_halo3         = 0.0333;\n");
	mfprintf(b,"#declare bond_color_halo3     = < 1, 1, 1, 0, 0.5 >;\n");
	mfprintf(b,"#declare bond_specular_halo3  = 0;\n");
	mfprintf(b,"#declare bond_reflect_halo3   = 0;\n");
	mfprintf(b,"#declare bond_ambient_halo3   = 1.0;\n");
	mfprintf(b,"#declare bond_diffuse_halo3   = 0;\n");

	mfprintf(b,"\n");

	mfprintf(b,"#declare bond_draw_stub       = true;\n");
	mfprintf(b,"#declare bond_r_stub          = 0.004;\n");
	mfprintf(b,"#declare bond_l_stub          = 0.004;\n");
	mfprintf(b,"#declare bond_color_stub      = < 0, 0, 0, 0, 0 >;\n");
	mfprintf(b,"#declare bond_specular_stub   = 0;\n");
	mfprintf(b,"#declare bond_reflect_stub    = 0;\n");
	mfprintf(b,"#declare bond_ambient_stub    = 1.0;\n");
	mfprintf(b,"#declare bond_diffuse_stub    = 0;\n");

	mfprintf(b,"\n");

	mfprintf(b,"/**** Element Colors ****/\n");
	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		if (z == g_iVirtAtomType)
			continue;
		el = ((CAtom*)g_oaAtoms[z])->m_pElement;
		mfprintf(b,"#declare elem_%s_col   = < %f, %f, %f, 0, 0 >;\n",el->m_sLabel,el->m_iColorR/255.0,el->m_iColorG/255.0,el->m_iColorB/255.0);
	}

	mfprintf(b,"\n/**** Element Radii ****/\n");
	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		if (z == g_iVirtAtomType)
			continue;
		el = ((CAtom*)g_oaAtoms[z])->m_pElement;
		mfprintf(b,"#declare elem_%s_r  = %f;\n",el->m_sLabel,el->m_fRadius/1000.0);
	}
	mfprintf(b,"\n");

	cam[0] = 0;
	cam[1] = 0;
	cam[2] = 40.0;

	mfprintf(b,"global_settings {\n");
	mfprintf(b,"  assumed_gamma 1\n");
	mfprintf(b,"/*  radiosity {\n");
	mfprintf(b,"    pretrace_start 0.08\n");
	mfprintf(b,"    pretrace_end   0.04\n");
	mfprintf(b,"    count 100\n\n");
	mfprintf(b,"    nearest_count 5\n");
	mfprintf(b,"    error_bound 0.4\n");
	mfprintf(b,"    recursion_limit 1\n\n");
	mfprintf(b,"    low_error_factor .5\n");
	mfprintf(b,"    gray_threshold 0.0\n");
	mfprintf(b,"    minimum_reuse 0.015\n");
	mfprintf(b,"    brightness 1\n\n");
	mfprintf(b,"    adc_bailout 0.01/2\n");
	mfprintf(b,"  }*/\n");
	mfprintf(b,"}\n");

	mfprintf(b,"\ncamera {\n");
	mfprintf(b,"	location <%f, %f, %f>\n",cam[0],cam[1],cam[2]);
	mfprintf(b,"	sky y\n");
	mfprintf(b,"	right -0.06*x*image_width/image_height\n");
	mfprintf(b,"	up 0.06*y\n");
	mfprintf(b,"	look_at <0, 0, 0>\n");
	mfprintf(b,"}\n");
	mfprintf(b,"\n");

	mfprintf(b,"// Solid background\n");
	mfprintf(b,"background { rgb < 0.15, 0.1, 0.3 > }\n");
	mfprintf(b,"\n");

	mfprintf(b,"// Gradient background\n");
	mfprintf(b,"sky_sphere {\n");
	mfprintf(b,"  pigment {\n");
	mfprintf(b,"    gradient y\n");
	mfprintf(b,"    color_map {\n");
	mfprintf(b,"      [ 0 color rgb < 0.05, 0.05, 0.05 > ]\n");
	mfprintf(b,"      [ 1 color rgb < 0.20, 0.16, 0.50 > ]\n");
	mfprintf(b,"    }\n");
	mfprintf(b,"    scale 0.1\n");
	mfprintf(b,"    translate -0.05\n");
	mfprintf(b,"  }\n");
	mfprintf(b,"}\n\n");

	mfprintf(b,"/**** Invisible, only for Radiosity ****/\n");
	mfprintf(b,"sphere {\n");
	mfprintf(b,"  <0, 0, 0>, 1\n");
	mfprintf(b,"  texture {\n");
	mfprintf(b,"    pigment {color rgb < 1.0, 1.0, 1.0 > }\n");
	mfprintf(b,"    finish { diffuse 0 ambient 1 }\n");
	mfprintf(b,"  }\n");
	mfprintf(b,"  hollow on\n");
	mfprintf(b,"  no_shadow\n");
	mfprintf(b,"  no_image\n");
	mfprintf(b,"  scale 30000\n");
	mfprintf(b,"}\n\n");

	mfprintf(b,"light_source { < -8, 20, 20 > color rgb 0.8 }\n");
	mfprintf(b,"//light_source { < 25, 12, 20 > color rgb 0.5 }\n\n");

	mfprintf(b,"#macro m_atom_color(col)\n");
	mfprintf(b,"  #ifdef(atom_color)\n");
	mfprintf(b,"    atom_color\n");
	mfprintf(b,"  #else #if (defined(atom_trans))\n");
	mfprintf(b,"    col + < 0, 0, 0, 0, atom_trans >\n");
	mfprintf(b,"  #else\n");
	mfprintf(b,"    col\n");
	mfprintf(b,"  #end #end\n");
	mfprintf(b,"#end\n");

	mfprintf(b,"\nunion {\n");

	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		m = (CMolecule*)g_oaMolecules[z];

		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];

			for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
			{
				if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
					continue;

				el = ((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_pElement;

				for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
				{
					o = ((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4);
					vec1 = m_vaCoords[o];
					vec1 /= 1000.0;

					mfprintf(b,"#if (atom_draw)\n");
					mfprintf(b,"  sphere { <%g, %g, %g>, elem_%s_r*atom_r\n",vec1[0],vec1[1],vec1[2],el->m_sLabel);
					mfprintf(b,"    pigment { rgbft m_atom_color(elem_%s_col) } finish { reflection atom_reflect specular atom_specular ambient atom_ambient diffuse atom_diffuse } }\n",el->m_sLabel);
					mfprintf(b,"#end\n");

					vec2 = cam - vec1;
					vec2.Normalize();

	//				if (shadow)
					{
						mfprintf(b,"#if (atom_draw_halo1)\n");
						mfprintf(b,"  disc { < %g - (%g) * atom_d_halo1, %g - (%g) * atom_d_halo1, %g - (%g) * atom_d_halo1 >,\n",vec1[0],vec2[0],vec1[1],vec2[1],vec1[2],vec2[2]);
						mfprintf(b,"    < %g, %g, %g >, (elem_%s_r * atom_r) + atom_r_halo1, elem_%s_r * atom_r\n",vec2[0],vec2[1],vec2[2],el->m_sLabel,el->m_sLabel);
						mfprintf(b,"    pigment { rgbft atom_color_halo1 } finish { reflection atom_reflect_halo1 specular atom_specular_halo1 ambient atom_ambient_halo1 diffuse atom_diffuse_halo1 } no_reflection no_radiosity }\n");
						mfprintf(b,"#end\n");
					}

	//				if (halo)
					{
						mfprintf(b,"#if (atom_draw_halo2)\n");
						mfprintf(b,"  disc { < %g - (%g) * atom_d_halo2, %g - (%g) * atom_d_halo2, %g - (%g) * atom_d_halo2 >,\n",vec1[0],vec2[0],vec1[1],vec2[1],vec1[2],vec2[2]);
						mfprintf(b,"    < %g, %g, %g >, (elem_%s_r * atom_r) + atom_r_halo1 + atom_r_halo2, elem_%s_r * atom_r\n",vec2[0],vec2[1],vec2[2],el->m_sLabel,el->m_sLabel);
						mfprintf(b,"    pigment { rgbft atom_color_halo2 } finish { reflection atom_reflect_halo2 specular atom_specular_halo2 ambient atom_ambient_halo2 diffuse atom_diffuse_halo2 } no_reflection no_radiosity }\n");
						mfprintf(b,"#end\n");

						mfprintf(b,"#if (atom_draw_halo3)\n");
						mfprintf(b,"  disc { < %g - (%g) * atom_d_halo3, %g - (%g) * atom_d_halo3, %g - (%g) * atom_d_halo3 >,\n",vec1[0],vec2[0],vec1[1],vec2[1],vec1[2],vec2[2]);
						mfprintf(b,"    < %g, %g, %g >, (elem_%s_r * atom_r) + atom_r_halo1 + atom_r_halo2 + atom_r_halo3, elem_%s_r * atom_r\n",vec2[0],vec2[1],vec2[2],el->m_sLabel,el->m_sLabel);
						mfprintf(b,"    pigment { rgbft atom_color_halo3 } finish { reflection atom_reflect_halo3 specular atom_specular_halo3 ambient atom_ambient_halo3 diffuse atom_diffuse_halo3 } no_reflection no_radiosity }\n");
						mfprintf(b,"#end\n");
					}
				}
			}

			for (z3=0;z3<sm->m_oaBonds.GetSize();z3++)
			{
				mb = (CMolBond*)sm->m_oaBonds[z3];
				o = mb->m_iAtomOffset[0];
				o2 = mb->m_iAtomOffset[1];

				el = ((CAtom*)g_oaAtoms[g_waAtomRealElement[o]])->m_pElement;
				el2 = ((CAtom*)g_oaAtoms[g_waAtomRealElement[o2]])->m_pElement;

				vec1 = m_vaCoords[o];
				vec1 /= 1000.0;

				vec2 = m_vaCoords[o2];
				vec2 /= 1000.0;

				if ((vec1-vec2).GetLength() > 0.3)
					continue;

				vec3 = (vec1/el->m_fRadius + vec2/el2->m_fRadius) / (1.0/el->m_fRadius+1.0/el2->m_fRadius);

				vec3b = vec2 - vec1;
				vec3b.Normalize();

				mfprintf(b,"#local vec1c = < %g + (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ),\n",vec1[0],vec3b[0],el->m_sLabel,el->m_sLabel,el->m_sLabel,el->m_sLabel);
				mfprintf(b,"  %g + (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ),\n",vec1[1],vec3b[1],el->m_sLabel,el->m_sLabel,el->m_sLabel,el->m_sLabel);
				mfprintf(b,"  %g + (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ) >;\n",vec1[2],vec3b[2],el->m_sLabel,el->m_sLabel,el->m_sLabel,el->m_sLabel);
				mfprintf(b,"#local vec2c = < %g - (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ),\n",vec2[0],vec3b[0],el2->m_sLabel,el2->m_sLabel,el2->m_sLabel,el2->m_sLabel);
				mfprintf(b,"  %g - (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ),\n",vec2[1],vec3b[1],el2->m_sLabel,el2->m_sLabel,el2->m_sLabel,el2->m_sLabel);
				mfprintf(b,"  %g - (%g) * ( elem_%s_r*atom_r - ( elem_%s_r*atom_r - 0.5*sqrt(4*elem_%s_r*elem_%s_r*atom_r*atom_r - 4*bond_r*bond_r ) ) ) >;\n",vec2[2],vec3b[2],el2->m_sLabel,el2->m_sLabel,el2->m_sLabel,el2->m_sLabel);

				mfprintf(b,"#if (bond_draw)\n");
				mfprintf(b,"  cylinder { < %g, %g, %g >, vec1c, bond_r open\n",vec3[0],vec3[1],vec3[2]);
				mfprintf(b,"    pigment { rgbft m_atom_color(elem_%s_col) } finish { reflection bond_reflect specular bond_specular ambient bond_ambient diffuse bond_diffuse } }\n",el->m_sLabel);
				mfprintf(b,"  cylinder { < %g, %g, %g >, vec2c, bond_r open\n",vec3[0],vec3[1],vec3[2]);
				mfprintf(b,"    pigment { rgbft m_atom_color(elem_%s_col) } finish { reflection bond_reflect specular bond_specular ambient bond_ambient diffuse bond_diffuse } }\n",el2->m_sLabel);
				mfprintf(b,"#end\n");

				vec3 = cam - (vec1 + vec2) / 2.0;
				vec3.Normalize();

				vec2b = vec2-vec1;
				vec1b = CrossP(vec3,vec2b);
				vec1b.Normalize();

				mfprintf(b,"#local vec3  = < %g, %g, %g >;\n",vec3[0],vec3[1],vec3[2]);
				mfprintf(b,"#local vec1b = < %g, %g, %g >;\n",vec1b[0],vec1b[1],vec1b[2]);

//				if (shadow)
				{
					mfprintf(b,"#if (bond_draw_halo1)\n");
					mfprintf(b,"  #local vecA = vec1c + vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecB = vec2c + vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecC = vec2c + vec1b * (bond_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecD = vec1c + vec1b * (bond_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  triangle { vecA, vecB, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
					mfprintf(b,"  triangle { vecA, vecD, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
					mfprintf(b,"  #local vecA = vec1c - vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecB = vec2c - vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecC = vec2c - vec1b * (bond_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  #local vecD = vec1c - vec1b * (bond_r) - vec3 * bond_d_halo1;\n");
					mfprintf(b,"  triangle { vecA, vecB, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
					mfprintf(b,"  triangle { vecA, vecD, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo1 } finish { reflection bond_reflect_halo1 specular bond_specular_halo1 ambient bond_ambient_halo1 diffuse bond_diffuse_halo1 } no_reflection no_radiosity }\n");
					mfprintf(b,"#end\n");
				}

//				if (halo)
				{
					mfprintf(b,"#if (bond_draw_halo2)\n");
					mfprintf(b,"  #local vecA = vec1c + vec1b * (bond_r + bond_r_halo1 + bond_r_halo2) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecB = vec2c + vec1b * (bond_r + bond_r_halo1 + bond_r_halo2) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecC = vec2c + vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecD = vec1c + vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  triangle { vecA, vecB, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo2 } finish { reflection bond_reflect_halo2 specular bond_specular_halo2 ambient bond_ambient_halo2 diffuse bond_diffuse_halo2 } no_reflection no_radiosity }\n");
					mfprintf(b,"  triangle { vecA, vecD, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo2 } finish { reflection bond_reflect_halo2 specular bond_specular_halo2 ambient bond_ambient_halo2 diffuse bond_diffuse_halo2 } no_reflection no_radiosity }\n");
					mfprintf(b,"  #local vecA = vec1c - vec1b * (bond_r + bond_r_halo1 + bond_r_halo2) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecB = vec2c - vec1b * (bond_r + bond_r_halo1 + bond_r_halo2) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecC = vec2c - vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecD = vec1c - vec1b * (bond_r + bond_r_halo1) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  triangle { vecA, vecB, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo2 } finish { reflection bond_reflect_halo2 specular bond_specular_halo2 ambient bond_ambient_halo2 diffuse bond_diffuse_halo2 } no_reflection no_radiosity }\n");
					mfprintf(b,"  triangle { vecA, vecD, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo2 } finish { reflection bond_reflect_halo2 specular bond_specular_halo2 ambient bond_ambient_halo2 diffuse bond_diffuse_halo2 } no_reflection no_radiosity }\n");
					mfprintf(b,"#end\n");

					mfprintf(b,"#if (bond_draw_halo3)\n");
					mfprintf(b,"  #local vecA = vec1c + vec1b * (bond_r + bond_r_halo1 + bond_r_halo2 + bond_r_halo3) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecB = vec2c + vec1b * (bond_r + bond_r_halo1 + bond_r_halo2 + bond_r_halo3) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecC = vec2c + vec1b * (bond_r + bond_r_halo1 + bond_r_halo2) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecD = vec1c + vec1b * (bond_r + bond_r_halo1 + bond_r_halo2) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  triangle { vecA, vecB, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo3 } finish { reflection bond_reflect_halo3 specular bond_specular_halo3 ambient bond_ambient_halo3 diffuse bond_diffuse_halo3 } no_reflection no_radiosity }\n");
					mfprintf(b,"  triangle { vecA, vecD, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo3 } finish { reflection bond_reflect_halo3 specular bond_specular_halo3 ambient bond_ambient_halo3 diffuse bond_diffuse_halo3 } no_reflection no_radiosity }\n");
					mfprintf(b,"  #local vecA = vec1c - vec1b * (bond_r + bond_r_halo1 + bond_r_halo2 + bond_r_halo3) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecB = vec2c - vec1b * (bond_r + bond_r_halo1 + bond_r_halo2 + bond_r_halo3) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecC = vec2c - vec1b * (bond_r + bond_r_halo1 + bond_r_halo2) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  #local vecD = vec1c - vec1b * (bond_r + bond_r_halo1 + bond_r_halo2) - vec3 * bond_d_halo2;\n");
					mfprintf(b,"  triangle { vecA, vecB, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo3 } finish { reflection bond_reflect_halo3 specular bond_specular_halo3 ambient bond_ambient_halo3 diffuse bond_diffuse_halo3 } no_reflection no_radiosity }\n");
					mfprintf(b,"  triangle { vecA, vecD, vecC\n");
					mfprintf(b,"    pigment { rgbft bond_color_halo3 } finish { reflection bond_reflect_halo3 specular bond_specular_halo3 ambient bond_ambient_halo3 diffuse bond_diffuse_halo3 } no_reflection no_radiosity }\n");
					mfprintf(b,"#end\n");
				}

//				if (shadow)
				{
					vec3 = vec2 - vec1;
					vec3.Normalize();

					mfprintf(b,"#if (bond_draw_stub)\n");
					mfprintf(b,"  #local vec3  = < %g, %g, %g >;\n",vec3[0],vec3[1],vec3[2]);
					mfprintf(b,"  difference {\n");
					mfprintf(b,"    cylinder { vec1c, vec1c + vec3 * bond_l_stub, bond_r + bond_r_stub }\n");
					mfprintf(b,"    cylinder { vec1c - vec3 * (bond_l_stub+0.001), vec1c + vec3 * (bond_l_stub+0.001), bond_r }\n");
					mfprintf(b,"    pigment { rgbft bond_color_stub } finish { reflection bond_reflect_stub specular bond_specular_stub ambient bond_ambient_stub diffuse bond_diffuse_stub } }\n");
					mfprintf(b,"  difference {\n");
					mfprintf(b,"    cylinder { vec2c, vec2c - vec3 * bond_l_stub, bond_r + bond_r_stub }\n");
					mfprintf(b,"    cylinder { vec2c + vec3 * (bond_l_stub+0.001), vec2c - vec3 * (bond_l_stub+0.001), bond_r }\n");
					mfprintf(b,"    pigment { rgbft bond_color_stub } finish { reflection bond_reflect_stub specular bond_specular_stub ambient bond_ambient_stub diffuse bond_diffuse_stub } }\n");
					mfprintf(b,"#end\n");
				}
			}
		}
	}

	mfprintf(b,"\n  no_shadow\n}\n\n");

	fclose(b);
}


void CTimeStep::DumpDipoles()
{
	int z, z2, z3, z4, ti;
	CMolecule *m;
	CSingleMolecule *sm;
	CxDVector3 dc;

	if (!(g_bDipole && g_bDumpDipoleVector))
		return;

	fprintf(g_fDumpDipole,"%lu",g_iSteps);
	for (z=0;z<g_oaMolecules.GetSize();z++)
	{
		if (g_oaDumpDipoleVector[z] == NULL)
			continue;
		m = (CMolecule*)g_oaMolecules[z];
		for (z2=0;z2<((CxIntArray*)g_oaDumpDipoleVector[z])->GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[((CxIntArray*)g_oaDumpDipoleVector[z])->GetAt(z2)]];
			fprintf(g_fDumpDipole,";  %.10G;  %.10G;  %.10G",sm->m_vDipole[0],sm->m_vDipole[1],sm->m_vDipole[2]);
			if (g_bDumpDipoleAbs)
				fprintf(g_fDumpDipole,";  %.10G",sm->m_vDipole.GetLength());
		}
	}
	fprintf(g_fDumpDipole,"\n");

	if (g_bDumpDipoleXYZ)
	{
		fprintf(g_fDumpDipoleXYZ,"%d\n\n",g_iDumpDipoleXYZAtoms);
		fprintf(g_fDumpDipolePURE,"%d\n\n",g_iDumpDipoleSMCount);

		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			if (g_oaDumpDipoleVector[z] == NULL)
				continue;
			m = (CMolecule*)g_oaMolecules[z];
			for (z2=0;z2<((CxIntArray*)g_oaDumpDipoleVector[z])->GetSize();z2++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[((CxIntArray*)g_oaDumpDipoleVector[z])->GetAt(z2)]];

				for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
				{
					if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
						continue;
					for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
					{
						ti = ((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4);
						fprintf(g_fDumpDipoleXYZ,"%s  %12f  %12f  %12f\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,m_vaCoords[ti][0]/100.0,m_vaCoords[ti][1]/100.0,m_vaCoords[ti][2]/100.0);
					}
				}
			}
		}

		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			if (g_oaDumpDipoleVector[z] == NULL)
				continue;
			m = (CMolecule*)g_oaMolecules[z];
			for (z2=0;z2<((CxIntArray*)g_oaDumpDipoleVector[z])->GetSize();z2++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[((CxIntArray*)g_oaDumpDipoleVector[z])->GetAt(z2)]];
				dc = m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[m->m_iDipoleCenterType])->GetAt(m->m_iDipoleCenterIndex)];
				fprintf(g_fDumpDipoleXYZ,"B  %12f  %12f  %12f\n",dc[0]/100.0,dc[1]/100.0,dc[2]/100.0);
				fprintf(g_fDumpDipoleXYZ,"B  %12f  %12f  %12f\n",dc[0]/100.0+sm->m_vDipole[0]*g_fDumpDipoleScale,dc[1]/100.0+sm->m_vDipole[1]*g_fDumpDipoleScale,dc[2]/100.0+sm->m_vDipole[2]*g_fDumpDipoleScale);

				fprintf(g_fDumpDipolePURE,"%s  %12f  %12f  %12f\n",(const char*)((CMolecule*)g_oaMolecules[sm->m_iMolType])->m_sName,sm->m_vDipole[0]*g_fDumpDipoleScale,sm->m_vDipole[1]*g_fDumpDipoleScale,sm->m_vDipole[2]*g_fDumpDipoleScale);
			}
		}
	}
}


CTimeStep::CTimeStep()
{
	m_pComment = NULL;
	m_pVolumetricData = NULL;
	m_pVolumetricDataTimeDev = NULL;
	m_pCurrentDensity = NULL;
	m_pPosDomainEngine = NULL;

	m_vaCoords.SetName("CTimeStep::m_vaCoords");
	m_vaCoords_Unfolded.SetName("CTimeStep::m_vaCoords_Unfolded");
	m_vaCoords_Original.SetName("CTimeStep::m_vaCoords_Original");
	m_vaVelocities.SetName("CTimeStep::m_vaVelocities");
	m_vaForces.SetName("CTimeStep::m_vaForces");
	m_paLabels.SetName("CTimeStep::m_paLabels");
	m_paMol2Types.SetName("CTimeStep::m_paMol2Types");
	m_faCharge.SetName("CTimeStep::m_faCharge");
}


CTimeStep::~CTimeStep()
{
/*	if (m_pLabels != NULL)
		delete[] m_pLabels;*/
	if (m_pComment != NULL)
	{
		delete[] m_pComment;
		m_pComment = NULL;
	}
	if (m_pVolumetricData != NULL)
	{
		delete m_pVolumetricData;
		m_pVolumetricData = NULL;
	}
	if (m_pVolumetricDataTimeDev != NULL) {
		delete m_pVolumetricDataTimeDev;
		m_pVolumetricDataTimeDev = NULL;
	}
	if (m_pCurrentDensity != NULL) {
		delete m_pCurrentDensity;
		m_pCurrentDensity = NULL;
	}
}

