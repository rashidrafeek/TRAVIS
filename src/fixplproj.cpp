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

#include "fixplproj.h"
#include "globalvar.h"


const char *GetRevisionInfo_fixplproj(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_fixplproj() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



bool CFixedPlProjObservation::Parse(int i) {

	int z, z2, z3, z4, ti;
	std::vector<int> tia;
	CxString buf;
	CMolecule *m;
	CSingleMolecule *sm;
	CAtomGroup ag, ag2;


	mprintf(WHITE,"\n    >>> Observation %d >>>\n\n",i+1);

	try { m_p2DF = new C2DF(); } catch(...) { m_p2DF = NULL; }
	if (m_p2DF == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	m_vPlaneBase[0] = AskFloat("      Enter X coordinate of base point (in pm): [0] ",0);
	m_vPlaneBase[1] = AskFloat("      Enter Y coordinate of base point (in pm): [0] ",0);
	m_vPlaneBase[2] = AskFloat("      Enter Z coordinate of base point (in pm): [0] ",0);

	mprintf("\n");

_nvagain:
	m_vPlaneVec1[0] = AskFloat("      Enter X coordinate of first spanning vector: [1.0] ",1.0);
	m_vPlaneVec1[1] = AskFloat("      Enter Y coordinate of first spanning vector: [0]   ",0);
	m_vPlaneVec1[2] = AskFloat("      Enter Z coordinate of first spanning vector: [0]   ",0);
	if (m_vPlaneVec1.GetLength() == 0) {
		eprintf("\n      Error: The length of the spanning vector is zero.\n\n");
		goto _nvagain;
	}

	mprintf("\n");

_nvagain2:
	m_vPlaneVec2[0] = AskFloat("      Enter X coordinate of second spanning vector: [0]   ",0);
	m_vPlaneVec2[1] = AskFloat("      Enter Y coordinate of second spanning vector: [1.0] ",1.0);
	m_vPlaneVec2[2] = AskFloat("      Enter Z coordinate of second spanning vector: [0]   ",0);
	if (m_vPlaneVec2.GetLength() == 0) {
		eprintf("\n      Error: The length of the spanning vector is zero.\n\n");
		goto _nvagain2;
	}

	mprintf("\n      Normalizing spanning vectors...\n");
	m_vPlaneVec1.Normalize();
	m_vPlaneVec2.Normalize();
	m_vPlaneNormal = CrossP(m_vPlaneVec1,m_vPlaneVec2);
	m_vPlaneNormal.Normalize();
	mprintf("        Vector 1:     ( %9.6f | %9.6f | %9.6f )\n",m_vPlaneVec1[0],m_vPlaneVec1[1],m_vPlaneVec1[2]);
	mprintf("        Vector 2:     ( %9.6f | %9.6f | %9.6f )\n",m_vPlaneVec2[0],m_vPlaneVec2[1],m_vPlaneVec2[2]);
	mprintf("        Plane Normal: ( %9.6f | %9.6f | %9.6f )\n\n",m_vPlaneNormal[0],m_vPlaneNormal[1],m_vPlaneNormal[2]);

	if (m_vPlaneNormal.GetLength() == 0) {
		eprintf("      Error: Both spanning vectors are linear. No plane defined.\n\n");
		return false;
	}

	m_bDifference = AskYesNo("      Observe absolute atom position (n) or difference between 2 atoms in a molecule (y)? [yes] ",true);

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];

		mprintf("\n");
		if (AskYesNo("      Observe some atoms from %s molecules (y/n)? [no] ",false,m->m_sName)) {

			tia.clear();
			if (AskYesNo("      Take into account all %s molecules (y) or only some of them (n)? [yes] ",true,m->m_sName)) {

				for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
					tia.push_back(z2);

			} else {

				mprintf("\n");
_smagain:
				AskString_ND("      Which of the %d %s molecules to take into account (e.g. 3-7,9): ",&buf,m->m_laSingleMolIndex.GetSize(),m->m_sName);
				if (!ParseIntList((const char*)buf,tia)) {
					eprintf("\n      Invalid input.\n\n");
					goto _smagain;
				}
				mprintf("      Observing %lu molecules of %s.\n\n",(unsigned long)tia.size(),m->m_sName);
			}

			ti = (int)m_iaAtomList.size();
			if (m_bDifference) {

				mprintf("      You can now enter several pairs of atoms between which the distance vector is processed.\n");

				while (true) {
					mprintf("\n");
_atagain2:
					AskString("      Enter first atom of pair from %s (e.g. C3): [done] ",&buf,"",m->m_sName);
					if (buf.GetLength() == 0)
						break;
					ag.Reset();
					if (!ag.ParseAtoms(m,(const char*)buf)) {
						eprintf("\n      Invalid input.\n\n");
						goto _atagain2;
					}
					if (ag.m_iAtomGes != 1) {
						eprintf("\n      Please enter only one atom at a time.\n\n");
						goto _atagain2;
					}

_atagain3:
					AskString("      Enter second atom of pair from %s (e.g. C3): [done] ",&buf,"",m->m_sName);
					if (buf.GetLength() == 0)
						break;
					ag2.Reset();
					if (!ag2.ParseAtoms(m,(const char*)buf)) {
						eprintf("\n      Invalid input.\n\n");
						goto _atagain3;
					}
					if (ag2.m_iAtomGes != 1) {
						eprintf("\n      Please enter only one atom at a time.\n\n");
						goto _atagain3;
					}

					for (z2=0;z2<(int)tia.size();z2++) {

						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[tia[z2]]];
						m_iaAtomList.push_back(((CxIntArray*)sm->m_oaAtomOffset[ag.m_baAtomType[0]])->GetAt(((CxIntArray*)ag.m_oaAtoms[0])->GetAt(0)));
						m_iaAtomList.push_back(((CxIntArray*)sm->m_oaAtomOffset[ag2.m_baAtomType[0]])->GetAt(((CxIntArray*)ag2.m_oaAtoms[0])->GetAt(0)));
					}

				}

				mprintf("\n      Added %d pairs of atoms for molecule %s.\n",((int)m_iaAtomList.size()-ti)/2,m->m_sName);

			} else {

_atagain:
				AskString("      Which atoms in %s to observe (e.g. C3)? [#2] ",&buf,"#2",m->m_sName);
				ag.Reset();
				if (!ag.ParseAtoms(m,(const char*)buf)) {
					eprintf("\n      Invalid input.\n\n");
					goto _atagain;
				}

				for (z2=0;z2<(int)tia.size();z2++) {

					sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[tia[z2]]];

					for (z3=0;z3<ag.m_oaAtoms.GetSize();z3++) {

						for (z4=0;z4<((CxIntArray*)ag.m_oaAtoms[z3])->GetSize();z4++)
							m_iaAtomList.push_back(((CxIntArray*)sm->m_oaAtomOffset[ag.m_baAtomType[z3]])->GetAt(((CxIntArray*)ag.m_oaAtoms[z3])->GetAt(z4)));
					}
				}

				mprintf("\n      Added %d atoms for molecule %s.\n",(int)m_iaAtomList.size()-ti,m->m_sName);
			}
		}
	}
	mprintf("\n");

	mprintf("      You can choose to take into account only atoms within a certain maximum distance\n");
	mprintf("      from the plane. For atom pairs, the center is used to compute the plane distance.\n\n");

	m_fLayerWidth = AskFloat("      Enter maximum distance from plane for atoms to take into account (in pm): [infinite] ",0);

	mprintf("\n");
	m_p2DF->m_fMinVal[0] = AskFloat("      Enter minimum value in X direction (in pm): [-1000.0] ",-1000.0);
	m_p2DF->m_fMaxVal[0] = AskFloat("      Enter maximum value in X direction (in pm): [ 1000.0] ",1000.0);
	m_p2DF->m_fMinVal[1] = AskFloat("      Enter minimum value in Y direction (in pm): [-1000.0] ",-1000.0);
	m_p2DF->m_fMaxVal[1] = AskFloat("      Enter maximum value in Y direction (in pm): [ 1000.0] ",1000.0);

	mprintf("\n");
	m_p2DF->m_iRes[0] = AskUnsignedInteger("      Enter resolution in X direction: [100] ",100);
	m_p2DF->m_iRes[1] = AskUnsignedInteger("      Enter resolution in Y direction: [100] ",100);

	m_p2DF->Create();

	m_p2DF->SetLabelX("Plane Vector 1 / pm");
	m_p2DF->SetLabelY("Plane Vector 2 / pm");
	m_p2DF->SetLabelZ("Occurrence");

	m_sName.sprintf("fixplproj_obs%d",i+1);

	if (m_bDifference)
		m_sName.strcat("_diff");
	else
		m_sName.strcat("_abs");

	if (m_fLayerWidth > 0) {
		m_pPointPlaneDF = new CDF();
		m_pPointPlaneDF->m_fMinVal = -MAX3(g_fBoxX,g_fBoxY,g_fBoxZ);
		m_pPointPlaneDF->m_fMaxVal = MAX3(g_fBoxX,g_fBoxY,g_fBoxZ);
		m_pPointPlaneDF->m_iResolution = 300;
		m_pPointPlaneDF->Create();
		m_pPointPlaneDF->SetLabelX("Distance from plane / pm");
		m_pPointPlaneDF->SetLabelY("Occurrence");
	}
	
	mprintf(WHITE,"\n    <<< Observation %d <<<\n\n",i+1);

	return true;
}


bool CFixedPlProj::Parse() {

	CFixedPlProjObservation *obs;


	mprintf(WHITE,"\n>>> Fixed Plane Projection Function >>>\n\n");

	mprintf("    The plane will be fixed in terms of simulation cell coordinates.\n");
	mprintf("    The plane will be defined via a base point and two spanning vectors.\n\n");

	do {

		obs = new CFixedPlProjObservation();
		if (!obs->Parse((int)m_oaObservations.size()))
			return false;
		m_oaObservations.push_back(obs);

	} while (AskYesNo("    Add another observation (y/n)? [no] ",false));

	mprintf("\n    Added %lu observations.\n",(unsigned long)m_oaObservations.size());

	mprintf(WHITE,"\n<<< Fixed Plane Projection Function <<<\n\n");

	return true;
}


void CFixedPlProj::Process(CTimeStep *ts) {

	int z, z2;
	double dist, ca, cb;
	CFixedPlProjObservation *obs;
	CxDVector3 avg, diff;
	CTimeStep t2;

	
	t2.CopyFrom(ts);
	t2.FoldMolecules();
	
	for (z=0;z<(int)m_oaObservations.size();z++) {

		obs = m_oaObservations[z];

		if (obs->m_bDifference) {

			for (z2=0;z2<(int)obs->m_iaAtomList.size()/2;z2++) {

				if (obs->m_fLayerWidth > 0) {

					avg = (t2.m_vaCoords[obs->m_iaAtomList[z2*2+1]] + t2.m_vaCoords[obs->m_iaAtomList[z2*2]]) / 2.0;

					dist = DotP(obs->m_vPlaneNormal,avg-obs->m_vPlaneBase);

					obs->m_pPointPlaneDF->AddToBin(dist);

					if (fabs(dist) > obs->m_fLayerWidth)
						continue;
				}

				diff = t2.m_vaCoords[obs->m_iaAtomList[z2*2+1]] - t2.m_vaCoords[obs->m_iaAtomList[z2*2]];

				PointRootCoefficients(obs->m_vPlaneVec1,obs->m_vPlaneVec2,diff,ca,cb);

				obs->m_p2DF->AddToBin(ca,cb);
			}

		} else {

			for (z2=0;z2<(int)obs->m_iaAtomList.size();z2++) {

				if (obs->m_fLayerWidth > 0) {

					dist = DotP(obs->m_vPlaneNormal,t2.m_vaCoords[obs->m_iaAtomList[z2]]-obs->m_vPlaneBase);

					obs->m_pPointPlaneDF->AddToBin(dist);

					if (fabs(dist) > obs->m_fLayerWidth)
						continue;
				}

				PointRootCoefficients(obs->m_vPlaneVec1,obs->m_vPlaneVec2,t2.m_vaCoords[obs->m_iaAtomList[z2]],ca,cb);

				obs->m_p2DF->AddToBin(ca,cb);
			}
		}
	}
}


void CFixedPlProj::Finish() {

	int z;
	CFixedPlProjObservation *obs;


	mprintf(WHITE,"\n    >>> Finishing Fixed Plane Projection Function >>>\n");

	for (z=0;z<(int)m_oaObservations.size();z++) {

		obs = m_oaObservations[z];

		mprintf(WHITE,"\n    *** Observation %d ***\n",z+1);

		mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",obs->m_p2DF->m_fBinEntries,obs->m_p2DF->m_fSkipEntries,ZeroDivide(obs->m_p2DF->m_fSkipEntries,obs->m_p2DF->m_fBinEntries+obs->m_p2DF->m_fSkipEntries)*100.0);
		if (obs->m_p2DF->m_fBinEntries == 0) {
			eprintf("    There were no bin entries. Check your function definition. Skipping this observation.\n\n");
			goto _skip;
		}

		obs->m_p2DF->CalcMaxEntry();
		mprintf("    Raw data range from %.0f to %.0f hits.\n",obs->m_p2DF->m_fMinEntry,obs->m_p2DF->m_fMaxEntry);
		mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",obs->m_p2DF->m_fMaxEntry,obs->m_p2DF->m_fEps);
		if (obs->m_p2DF->m_fEps > 1.0E-4) {
			eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
			mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
		}

		mprintf("    Normalizing integral value to %.2f.\n",1000000.0);
		obs->m_p2DF->NormalizeBinIntegral(1000000.0);

		obs->m_p2DF->CalcMaxEntry();
		mprintf("    Resulting data range from %.6f to %.6f.\n",obs->m_p2DF->m_fMinEntry,obs->m_p2DF->m_fMaxEntry);

		mprintf("    Saving Plane Projection DF triples as \"%s_triples.csv\"...\n",(const char*)obs->m_sName);
		obs->m_p2DF->Write((const char*)obs->m_sName,"","_triples.csv");
		mprintf("    Saving Plane Projection DF matrix as \"%s_matrix.csv\"...\n",(const char*)obs->m_sName);
		obs->m_p2DF->WriteCSV((const char*)obs->m_sName,"","_matrix.csv");
		mprintf("    Saving Plane Projection DF Mathematica Notebook as \"%s.nb\"...\n",(const char*)obs->m_sName);
		obs->m_p2DF->WriteMathematicaNb((const char*)obs->m_sName,"",".nb",false);
		mprintf("    Saving Plane Projection DF Gnuplot Input as \"%s.gp\"...\n",(const char*)obs->m_sName);
		obs->m_p2DF->WriteGnuplotInput((const char*)obs->m_sName,"","",false);

		if (obs->m_fLayerWidth > 0) {
			mprintf("    Saving plane distance distribution as \"%s_planedistance.csv\"...\n",(const char*)obs->m_sName);
			obs->m_pPointPlaneDF->Write((const char*)obs->m_sName,"","_planedistance.csv",true);
		}
_skip:;
	}

	mprintf(WHITE,"\n    <<< Finishing Fixed Plane Projection Function done <<<\n\n");
}


