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

#include "geodens.h"
#include "globalvar.h"
#include "maintools.h"


const char *GetRevisionInfo_geodens(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_geodens() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



CGeoDensObservation::CGeoDensObservation() {

	m_iObservation = -1;
	m_iType = GEODENS_TYPE_UNKNOWN;
	m_fTimeDev = NULL;
}


bool CGeoDensObservation::Parse() {

	int i, z, z2, z3, z4;
	CMolecule *m;
	CSingleMolecule *sm;
	CAtomGroup ag;
	CxString buf;
	std::vector<int> ml;


	mprintf(WHITE,"\n      >>> Geometric Density Observation %d >>>\n\n",m_iObservation+1);

	buf.sprintf("geodens_obs%d",m_iObservation+1);
	m_sName = (const char*)buf;

	mprintf("      The following geometric shapes are currently implemented:\n\n");
	mprintf("        1) Sphere\n");
	mprintf("        2) Cylinder\n");
	mprintf("\n");

_tagain:
	i = AskRangeInteger_ND("      Which type to use (1-2)? ",1,2);

	switch(i) {
		case 1:
			m_iType = GEODENS_TYPE_SPHERE;
			break;
		case 2:
			m_iType = GEODENS_TYPE_CYLINDER;
			break;
		default:
			eprintf("Invalid input.\n");
			goto _tagain;
	}

	mprintf("\n");

	if (m_iType == GEODENS_TYPE_SPHERE) {

		if (!ParseSphere()) {
			eprintf("Error in sphere definition.\n");
			return false;
		}

	} else if (m_iType == GEODENS_TYPE_CYLINDER) {

		if (!ParseCylinder()) {
			eprintf("Error in cylinder definition.\n");
			return false;
		}
	}

	mprintf(YELLOW,"      *** Observed Atom Definition ***\n\n");

	m_iaObserved.clear();

	m_sName += "_O";

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];

		if (!AskYesNo("      Observe atoms from molecule %d (%s) (y/n)? [no] ",false,z+1,m->m_sName))
			continue;

		if (m->m_iAtomGesNoVirt == 1) {
			mprintf("        %s is only one atom, there is no choice.\n",m->m_sName);
			ag.Reset();
			ag.m_pMolecule = m;
			ag.AddAtom(0,0,false);
			ag.SortAtoms();
			ag.BuildName();
		} else {
_as:
			ag.Reset();
			AskString("        Which atom(s) to take from %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ", &buf, "#2", m->m_sName );
			if (!ag.ParseAtoms(m,buf))
				goto _as;
		}

		m_sName += "_";
		m_sName += m->m_sName;
		if (m->m_iAtomGesNoVirt > 1) {
			m_sName += "_";
			m_sName += ag.m_sName;
		}

_ms:
		ml.clear();

		if (m->m_laSingleMolIndex.GetSize() > 1) {

			AskString("        Which molecules of %s to use (1-%d, e.g. \"1,3-7\")? [all] ",&buf,"",m->m_sName,m->m_laSingleMolIndex.GetSize());

			if (buf.GetLength() == 0) {
				for (z2=1;z2<=m->m_laSingleMolIndex.GetSize();z2++)
					ml.push_back(z2);
			} else if (!ParseIntList((const char*)buf,ml)) {
				eprintf("Invalid input.\n");
				goto _ms;
			}

			for (z2=0;z2<(int)ml.size();z2++) {
				if ((ml[z2] < 1) || (ml[z2] > m->m_laSingleMolIndex.GetSize())) {
					eprintf("Error: Molecule index %d out of range (allowed: 1 ... %d).\n",ml[z2],m->m_laSingleMolIndex.GetSize());
					goto _ms;
				}
				ml[z2]--;
			}

		} else
			ml.push_back(0);


		for (z2=0;z2<ag.m_oaAtoms.GetSize();z2++) {

			for (z3=0;z3<((CxIntArray*)ag.m_oaAtoms[z2])->GetSize();z3++) {

				for (z4=0;z4<(int)ml.size();z4++) {

					sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[ml[z4]]];

					m_iaObserved.push_back(((CxIntArray*)sm->m_oaAtomOffset[ag.m_baAtomType[z2]])->GetAt(((CxIntArray*)ag.m_oaAtoms[z2])->GetAt(z3)));
				}
			}
		}
	}

	mprintf("\n      Observing %lu atoms.\n\n",(unsigned long)m_iaObserved.size());

	m_bTimeDev = AskYesNo("      Save temporal development of this observation (y/n)? [yes] ",true);

	mprintf("\n      Name of this observation is \"%s\".\n\n",m_sName.c_str());


	mprintf(WHITE,"      <<< Geometric Density Observation %d <<<\n\n",m_iObservation+1);

	return true;
}


bool CGeoDensObservation::ParseSphere() {

	eprintf("Sphere: Not yet implemented.\n");

	return false;
}


bool CGeoDensObservation::ParseCylinder() {

	std::vector<int> set[2], ml;
	int z0, z, z2, z3, z4;
	CMolecule *m;
	CSingleMolecule *sm;
	CAtomGroup ag;
	CxString buf;
	bool tb;


	mprintf(YELLOW,"      *** Cylinder Definition ***\n\n");

	mprintf("      Lines are spawned between atoms from set A and set B.\n");
	mprintf("      Around these lines, cylinders with fixed radius are created.\n\n");

	m_sName += "_cylinder";

	for (z0=0;z0<2;z0++) {

		mprintf(WHITE,"      *** Definition of Atom Set %c ***\n\n",'A'+z0);

		buf.sprintf("_%c",'A'+z0);
		m_sName += (const char*)buf;

		set[z0].clear();

		for (z=0;z<g_oaMolecules.GetSize();z++) {

			m = (CMolecule*)g_oaMolecules[z];

			if (!AskYesNo("      Use atoms from molecule %d (%s) for this set (y/n)? [no] ",false,z+1,m->m_sName))
				continue;

			if (m->m_iAtomGesNoVirt == 1) {
				mprintf("        %s is only one atom, there is no choice.\n",m->m_sName);
				ag.Reset();
				ag.m_pMolecule = m;
				ag.AddAtom(0,0,false);
				ag.SortAtoms();
				ag.BuildName();
			} else {
_as:
				ag.Reset();
				AskString("        Which atom(s) to take from %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ", &buf, "#2", m->m_sName );
				if (!ag.ParseAtoms(m,buf))
					goto _as;
			}

			m_sName += "_";
			m_sName += m->m_sName;
			if (m->m_iAtomGesNoVirt > 1) {
				m_sName += "_";
				m_sName += ag.m_sName;
			}

_ms:
			ml.clear();

			if (m->m_laSingleMolIndex.GetSize() > 1) {

				AskString("        Which molecules of %s to use (1-%d, e.g. \"1,3-7\")? [all] ",&buf,"",m->m_sName,m->m_laSingleMolIndex.GetSize());

				if (buf.GetLength() == 0) {
					for (z2=1;z2<=m->m_laSingleMolIndex.GetSize();z2++)
						ml.push_back(z2);
				} else if (!ParseIntList((const char*)buf,ml)) {
					eprintf("Invalid input.\n");
					goto _ms;
				}

				for (z2=0;z2<(int)ml.size();z2++) {
					if ((ml[z2] < 1) || (ml[z2] > m->m_laSingleMolIndex.GetSize())) {
						eprintf("Error: Molecule index %d out of range (allowed: 1 ... %d).\n",ml[z2],m->m_laSingleMolIndex.GetSize());
						goto _ms;
					}
					ml[z2]--;
				}

			} else
				ml.push_back(0);


			for (z2=0;z2<ag.m_oaAtoms.GetSize();z2++) {

				for (z3=0;z3<((CxIntArray*)ag.m_oaAtoms[z2])->GetSize();z3++) {

					for (z4=0;z4<(int)ml.size();z4++) {

						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[ml[z4]]];

						set[z0].push_back(((CxIntArray*)sm->m_oaAtomOffset[ag.m_baAtomType[z2]])->GetAt(((CxIntArray*)ag.m_oaAtoms[z2])->GetAt(z3)));
					}
				}
			}
		}

		mprintf("\n      Set %c contains %lu atoms.\n\n",'A'+z0,(unsigned long)set[z0].size());
		if (set[z0].size() == 0) {
			eprintf("Error: Empty set is not allowed.\n");
			return false;
		}
	}

	mprintf("      Set definitions completed.\n\n");

	m_iaTuples.clear();

	if (AskYesNo("      Automatically create combinations (y) or enter combinations per hand (n)? [yes] ",true)) {

		mprintf("\n");
		mprintf("      Combinations of identical atoms (i.e., empty cylinders) will be excluded.\n\n");

		tb = AskYesNo("      Exclude combinations where both atoms are from the same molecule (y/n)? [yes] ",true);

		for (z=0;z<(int)set[0].size();z++) {

			for (z2=0;z2<(int)set[1].size();z2++) {

				if (set[0][z] == set[1][z2])
					continue;

				if ((g_waAtomMolIndex[set[0][z]] == g_waAtomMolIndex[set[1][z2]]) &&
					(g_laAtomSMLocalIndex[set[0][z]] == g_laAtomSMLocalIndex[set[1][z2]]) &&
					tb)
					continue;

				m_iaTuples.push_back(set[0][z]);
				m_iaTuples.push_back(set[1][z2]);
			}
		}

	} else {

		for (z0=0;z0<2;z0++) {
			mprintf("\n");
			mprintf("      Set %c contains the following %lu atoms:\n",'A'+z0,(unsigned long)set[0].size());
			for (z=0;z<(int)set[z0].size();z++)
				mprintf("        %03d) %s[%d] %s%d\n",
					z+1,
					((CMolecule*)g_oaMolecules[g_waAtomMolIndex[set[z0][z]]])->m_sName,
					g_laAtomSMLocalIndex[set[z0][z]]+1,
					(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[set[z0][z]]])->m_sName,
					g_waAtomMolNumber[set[z0][z]]+1
				);
		}
		mprintf("\n");

		while (true) {
_pairagain:
			AskString("      Enter comma-separated pair to add (e.g. \"1,3\"): [done] ",&buf,"");
			if (buf.GetLength() == 0)
				break;
			ml.clear();
			if (!ParseIntList(buf,ml)) {
				eprintf("Invalid input: Not a comma-separated list of integers.\n");
				goto _pairagain;
			}
			if (ml.size() != 2) {
				eprintf("Invalid input: Expected two comma-separated integers, found %lu.\n",(unsigned long)ml.size());
				goto _pairagain;
			}
			for (z0=0;z0<2;z0++) {
				if ((ml[z0] < 1) || (ml[z0] > (int)set[z0].size())) {
					eprintf("Invalid input: Index %d out of range (allowed 1 ... %lu, found %d).\n",z0+1,(unsigned long)set[z0].size(),ml[z0]);
					goto _pairagain;
				}
			}
			for (z0=0;z0<2;z0++)
				m_iaTuples.push_back(set[z0][ml[z0]-1]);
		}
	}

	mprintf("\n      %lu cylinders have been defined.\n\n",(unsigned long)m_iaTuples.size()/2);

//	for (z=0;z<(int)m_iaTuples.size()/2;z++)
//		mprintf("%d : %d\n",m_iaTuples[z*2],m_iaTuples[z*2+1]);

	if (m_iaTuples.size() == 0) {
		eprintf("At least one cylinder is required.\n");
		return false;
	}

	m_faParameters.push_back( AskFloat("      Enter radius of cylinders (in pm): [300] ",300) );

	mprintf("\n");
	mprintf(YELLOW,"      *** End of Cylinder Definition ***\n\n");

	return true;
}


void CGeoDensObservation::Initialize() {

	CxString buf;
	int z;

	if (m_bTimeDev) {
		buf.sprintf("%s_timedev.csv",m_sName.c_str());
		m_fTimeDev = OpenFileWrite((const char*)buf,true);
		mfprintf(m_fTimeDev,"# Step;  Uniform Density (nm^-3)");
		switch(m_iType) {
			case GEODENS_TYPE_CYLINDER:
				for (z=0;z<(int)m_iaTuples.size()/2;z++)
					mfprintf(m_fTimeDev,"; Cylinder %d Volume (pm^3); Count; Density (nm^-3); Ratio",z+1);
				mfprintf(m_fTimeDev,"\n");
				break;
			default:
				eprintf("CGeoDensObservation::Initialize(): Internal error: m_iType == %d.\n",m_iType);
				abort();
		}
	}
}


void CGeoDensObservation::ProcessStep(CTimeStep *ts) {

	int z, z2, count;
	double tvol, cvol, dp, r, len;
	CxDVector3 veca, vecb, vecab, veco;


	tvol = g_mBoxFromOrtho.Det(); // Cell volume in pm^3

	if (m_bTimeDev)
		mfprintf(m_fTimeDev,"%lu;  %f",g_iSteps,(double)m_iaObserved.size()/tvol*1.0E9);

	if (m_iType == GEODENS_TYPE_CYLINDER) {

		for (z=0;z<(int)m_iaTuples.size()/2;z++) {

			veca = ts->m_vaCoords[m_iaTuples[z*2]];
			vecb = ts->m_vaCoords[m_iaTuples[z*2+1]];
			vecab = FoldVector(vecb-veca);

			len = vecab.GetLength(); // Cylinder length

			cvol = Pi * pow2(m_faParameters[0]) * len; // Cylinder volume in pm^3

			vecab.Normalize();

			count = 0;

			for (z2=0;z2<(int)m_iaObserved.size();z2++) {

				veco = FoldVector( ts->m_vaCoords[m_iaObserved[z2]] - veca );

				dp = DotP( vecab, veco );

				if ((dp < 0) || (dp > len))
					continue;

				r = (veco - vecab*dp).GetLength();

				if (r > m_faParameters[0])
					continue;

				count++;
			}

			if (m_bTimeDev)
				mfprintf(m_fTimeDev,"; %f; %d; %f; %f",cvol,count,(double)count/cvol*1.0E9,((double)count/cvol)/((double)m_iaObserved.size()/tvol));
		}

	} else {

		eprintf("CGeoDensObservation::ProcessStep(): Internal error: m_iType == %d.\n",m_iType);
		abort();
	}

	if (m_bTimeDev)
		mfprintf(m_fTimeDev,"\n");
}


void CGeoDensObservation::Finish() {

	mprintf(WHITE,"    *** Observation %d ***\n",m_iObservation+1);

	if (m_bTimeDev)
		fclose(m_fTimeDev);
}



/*****************************************************************************************************************/


	
bool CGeoDensEngine::Parse() {

	CGeoDensObservation *o;


	mprintf(WHITE,"\n    >>> Geometric Density Analysis >>>\n");

	while (true) {

		o = new CGeoDensObservation();
		o->m_iObservation = (int)m_oaObservations.size();
		m_oaObservations.push_back(o);

		if (!o->Parse())
			return false;

		if (!AskYesNo("    Add another geometric density observation (y/n)? [no] ",false))
			break;
	}

	mprintf(WHITE,"\n    <<< Geometric Density Analysis <<<\n\n");

	return true;
}


void CGeoDensEngine::Initialize() {

	int z;

	mprintf("  Initializing Geometric Density Observations...\n");

	for (z=0;z<(int)m_oaObservations.size();z++) {
		mprintf("    Observation %d...\n",z+1);
		m_oaObservations[z]->Initialize();
	}
}


void CGeoDensEngine::ProcessStep(CTimeStep *ts) {

	int z;

	for (z=0;z<(int)m_oaObservations.size();z++)
		m_oaObservations[z]->ProcessStep(ts);
}


void CGeoDensEngine::Finish() {

	int z;

	mprintf(WHITE,"\n    Finishing Geometric Density Observations...\n");
	for (z=0;z<(int)m_oaObservations.size();z++)
		m_oaObservations[z]->Finish();
}


