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

#include "order.h"
#include "order_vector.h"
#include "globalvar.h"
#include "maintools.h"


const char *GetRevisionInfo_order_vector(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_order_vector() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



bool COrderAnalysis_Vector::Parse(CTimeStep *ts) {

	COrderAnalysis_VectorObs *obs;
	CxString buf;
	CAtomGroup ag;
	CMolecule *m;
	bool chain, b;
	char *p, *q, cbuf[16];
	int z, z2, z3, z4, ti, ia, it;
	std::vector<int> tia;
	std::vector<std::string> tsa;
	CSingleMolecule *sm;
	CMolBond *mb;
	COrderChain *oc;


	UNUSED(ts);

	mprintf(WHITE,"\n    >>> Analysis %d: Vector-Based Order Parameters >>>\n\n",m_iNumber+1);

	m_vFixedVector[0] = AskFloat("      Enter X component of fixed vector: [0.0] ",0);
	m_vFixedVector[1] = AskFloat("      Enter Y component of fixed vector: [0.0] ",0);
	m_vFixedVector[2] = AskFloat("      Enter Z component of fixed vector: [1.0] ",1.0);

	mprintf("\n");

	if (m_pParent->m_iChainsTotal != 0) {
		chain = AskYesNo("      Use C-H bonds along an alkyl chain as vectors (y/n)? [no] ",false);
		mprintf("\n");
	} else
		chain = false;

	if (chain) {

		AskString("      Form vectors to which elements as carbon bond partners? [H,D,F] ",&buf,"H,D,F");

		p = buf.GetWritePointer();

		while (*p != 0) {
			while (*p == ' ')
				p++;
			q = p;
			while ((*q != 0) && (*q != ',') && (*q != ' '))
				q++;
			memcpy(cbuf,p,q-p);
			cbuf[q-p] = 0;
			tsa.push_back((const char*)cbuf);
			while ((*q == ' ') || (*q == ','))
				q++;
			if (*q == 0)
				break;
			p = q;
		}

		if (g_oaMolecules.GetSize() == 1) {
			mprintf("      Using molecule 1.\n");
			ti = 0;
		} else
			ti = AskRangeInteger_ND("      Define vectors in which molecule (1-%d)? ",1,g_oaMolecules.GetSize(),g_oaMolecules.GetSize()) - 1;
		if (m_pParent->m_oaChains[ti].size() == 0) {
			eprintf("\n      Error: No alkyl chains have been detected in molecule %d.\n",ti+1);
			abort();
		}
		m = (CMolecule*)g_oaMolecules[ti];
_ilagain:
		tia.clear();
		AskString_ND("      Which chain(s) to use in molecule %d (comma separated, see list above)? ",&buf,ti+1);
		if (!ParseIntList((const char*)buf,tia)) {
			eprintf("\n      Invalid input.\n\n");
			goto _ilagain;
		}
		for (z=0;z<(int)tia.size();z++) {
			if ((tia[z] < 1) || (tia[z] > (int)m_pParent->m_oaChains[ti].size())) {
				eprintf("\n      Chain number %d outside of allowed range (%d .. %lu).\n\n",tia[z],1,(unsigned long)m_pParent->m_oaChains[ti].size());
				goto _ilagain;
			}
			tia[z]--;
			if (z > 0) {
				if (m_pParent->m_oaChains[ti][tia[0]]->m_iaAtoms.size() != m_pParent->m_oaChains[ti][tia[z]]->m_iaAtoms.size()) {
					eprintf("\n      Chains have different lengths (%d: %lu vs. %d: %lu).\n\n",tia[0]+1,(unsigned long)m_pParent->m_oaChains[ti][tia[0]]->m_iaAtoms.size(),tia[z]+1,(unsigned long)m_pParent->m_oaChains[ti][tia[z]]->m_iaAtoms.size());
					goto _ilagain;
				}
			}
		}
		mprintf("\n");
		mprintf("      Using %lu chains of length %lu:\n\n",(unsigned long)tia.size(),(unsigned long)m_pParent->m_oaChains[ti][tia[0]]->m_iaAtoms.size());

		sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
		for (z=0;z<(int)m_pParent->m_oaChains[ti][tia[0]]->m_iaAtoms.size();z++) {
			mprintf("        Vector %2d: ",z+1);
			obs = new COrderAnalysis_VectorObs();
			m_oaObservations.push_back(obs);
			obs->m_sName = m->m_sName;
			obs->m_iMolecule = ti;
			b = false;
			for (z2=0;z2<(int)tia.size();z2++) {
				oc = m_pParent->m_oaChains[ti][tia[z2]];
				for (z3=0;z3<sm->m_oaBonds.GetSize();z3++) {
					mb = (CMolBond*)sm->m_oaBonds[z3];
					if ((mb->m_iAtomType[0] == oc->m_iCarbonIndex) && (mb->m_iAtom[0] == oc->m_iaAtoms[z])) {
						ia = mb->m_iAtom[1];
						it = mb->m_iAtomType[1];
					} else if ((mb->m_iAtomType[1] == oc->m_iCarbonIndex) && (mb->m_iAtom[1] == oc->m_iaAtoms[z])) {
						ia = mb->m_iAtom[0];
						it = mb->m_iAtomType[0];
					} else
						continue;
					for (z4=0;z4<(int)tsa.size();z4++)
						if (mystricmp(tsa[z4].c_str(),((CAtom*)g_oaAtoms[m->m_baAtomIndex[it]])->m_sName) == 0)
							goto _elok;
					continue;
_elok:
					obs->m_iaAtomTypes.push_back(oc->m_iCarbonIndex);
					obs->m_iaAtomTypes.push_back(it);
					obs->m_iaAtoms.push_back(oc->m_iaAtoms[z]);
					obs->m_iaAtoms.push_back(ia);
					buf.sprintf("_C%d-%s%d",oc->m_iaAtoms[z]+1,(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[it]])->m_sName,ia+1);
					obs->m_sName += buf;
					if (b)
						mprintf(", ");
					b = true;
					mprintf("C%d --> %s%d",oc->m_iaAtoms[z]+1,(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[it]])->m_sName,ia+1);
				}

			}
			mprintf("\n");
		}
		mprintf("\n");
		goto _defdone;
	}


	while (true) {

		mprintf(WHITE,"      *** Observation %lu ***\n\n",(unsigned long)m_oaObservations.size());

		obs = new COrderAnalysis_VectorObs();

		m_oaObservations.push_back(obs);

		if (g_oaMolecules.GetSize() == 1) {
			mprintf("        Defining vector in molecule 1.\n");
			obs->m_iMolecule = 0;
		} else
			obs->m_iMolecule = AskRangeInteger_ND("        Define vector in which molecule (1-%d)? ",1,g_oaMolecules.GetSize(),g_oaMolecules.GetSize()) - 1;

		mprintf("\n");

		m = (CMolecule*)g_oaMolecules[obs->m_iMolecule];

		obs->m_sName = m->m_sName;

_again1:
		AskString_ND("        Which atom to use as vector base point (e.g. C3)? ",&buf);
		ag.Reset();
		if (!ag.ParseAtoms(m,(const char*)buf)) {
			eprintf("\n      Invalid input.\n\n");
			goto _again1;
		}
		if (ag.m_iAtomGes != 1) {
			eprintf("\n      Please enter only one atom at a time.\n\n");
			goto _again1;
		}
		obs->m_iaAtomTypes.push_back(ag.m_baAtomType[0]);
		obs->m_iaAtoms.push_back(((CxIntArray*)ag.m_oaAtoms[0])->GetAt(0));
		buf.sprintf("_%s%d-",(const char*)((CAtom*)g_oaAtoms[ag.m_baRealAtomType[0]])->m_sName,((CxIntArray*)ag.m_oaAtoms[0])->GetAt(0)+1);
		obs->m_sName += buf;

_again2:
		AskString_ND("        Which atom to use as vector tip point (e.g. C3)? ",&buf);
		ag.Reset();
		if (!ag.ParseAtoms(m,(const char*)buf)) {
			eprintf("\n      Invalid input.\n\n");
			goto _again2;
		}
		if (ag.m_iAtomGes != 1) {
			eprintf("\n      Please enter only one atom at a time.\n\n");
			goto _again2;
		}
		obs->m_iaAtomTypes.push_back(ag.m_baAtomType[0]);
		obs->m_iaAtoms.push_back(((CxIntArray*)ag.m_oaAtoms[0])->GetAt(0));
		buf.sprintf("%s%d",(const char*)((CAtom*)g_oaAtoms[ag.m_baRealAtomType[0]])->m_sName,((CxIntArray*)ag.m_oaAtoms[0])->GetAt(0)+1);
		obs->m_sName += buf;

		mprintf("\n");

		if (AskYesNo("        Add another vector to this observation (y/n)? [no] ",false))
			goto _again1;
		
		if (!AskYesNo("        Add another observation (y/n)? [no] ",false))
			break;

		mprintf("\n");
	}
	mprintf("\n");

_defdone:
	m_bHistograms = AskYesNo("      Create histograms of each vector order parameter (y) or only average values (n)? [yes] ",true);

	if (m_bHistograms) {
		ti = AskUnsignedInteger("      Which resolution to use for the histograms? [100] ",100);
		for (z=0;z<(int)m_oaObservations.size();z++) {
			obs = m_oaObservations[z];
			obs->m_pHistogram = new CDF();
			obs->m_pHistogram->m_fMinVal = 0.0;
			obs->m_pHistogram->m_fMaxVal = 180.0;
			obs->m_pHistogram->m_iResolution = ti;
			obs->m_pHistogram->Create();
			obs->m_pHistogram->SetLabelX("Angle");
			obs->m_pHistogram->SetLabelY("Occurrence");
		}
		mprintf("\n");
	}

	for (z=0;z<(int)m_oaObservations.size();z++) {
		obs = m_oaObservations[z];
		obs->m_fAverage = 0;
		obs->m_fAvgCounter = 0;
	}

	if (m_oaObservations.size() > 1)
		m_bCombined = AskYesNo("      Create an overview datafile with the average values of all order parameters (y/n)? [yes] ",true);
	else
		m_bCombined = false;

	mprintf(WHITE,"\n    <<< Analysis %d: Vector-Based Order Parameters <<<\n\n",m_iNumber+1);

	return true;
}


void COrderAnalysis_Vector::ProcessStep(CTimeStep *ts) {

	int z, z2, z3;
	COrderAnalysis_VectorObs *obs;
	CxDVector3 vec;
	CMolecule *m;
	CSingleMolecule *sm;
	double tf, tf2;


	for (z=0;z<(int)m_oaObservations.size();z++) {

		obs = m_oaObservations[z];
		m = (CMolecule*)g_oaMolecules[obs->m_iMolecule];

		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++) {

			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];

			for (z3=0;z3<(int)obs->m_iaAtoms.size()/2;z3++) {

				vec = ts->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[obs->m_iaAtomTypes[z3*2+1]])->GetAt(obs->m_iaAtoms[z3*2+1])] - 
					ts->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[obs->m_iaAtomTypes[z3*2]])->GetAt(obs->m_iaAtoms[z3*2])];

				vec = FoldVector(vec);

				vec.Normalize();

				tf = DotP(vec,m_vFixedVector);
				tf2 = 1.5*tf*tf - 0.5;

				obs->m_fAverage += tf2;
				obs->m_fAvgCounter++;

				if (m_bHistograms)
					obs->m_pHistogram->AddToBin(acos(tf)*180.0/Pi);
			}
		}
	}
}


void COrderAnalysis_Vector::Finish() {

	int z /*, z2*/;
	//double tf , tf2;
	//std::vector<double> tfa;
	COrderAnalysis_VectorObs *obs;
	FILE *a;
	CxString buf;


	mprintf(WHITE,"\n    *** Vector Order Parameter Analysis %d ***\n\n",m_iNumber+1);

	for (z=0;z<(int)m_oaObservations.size();z++) {

		obs = m_oaObservations[z];
		mprintf("      Observation %d:\n",z+1);
		obs->m_fAverage /= obs->m_fAvgCounter;
		mprintf("        Average value: %.6f\n",obs->m_fAverage);
		if (m_bHistograms) {
			buf.sprintf("order_vector%d_histogram%d_%s.csv",m_iNumber+1,z+1,(const char*)obs->m_sName);
			mprintf("        Writing histogram to %s ...\n",(const char*)buf);
			obs->m_pHistogram->AngleCorrect();
			obs->m_pHistogram->NormBinSum(100.0);
			obs->m_pHistogram->Integrate(false,1.0);
			obs->m_pHistogram->Write((const char*)buf,"","",true);
	/*		tfa.resize(obs->m_pHistogram->m_iResolution);
			for (z2=0;z2<obs->m_pHistogram->m_iResolution;z2++) {
				tf = -0.5 + 1.5*(z2+0.5)/obs->m_pHistogram->m_iResolution;
				tf2 = sin(acos(sqrt((tf+0.5)/1.5)));
				tfa[z2] = obs->m_pHistogram->m_pBin[z2] / tf2 / tf2;
			}
			tf = 0;
			for (z2=0;z2<obs->m_pHistogram->m_iResolution;z2++)
				tf += tfa[z2];
			for (z2=0;z2<obs->m_pHistogram->m_iResolution;z2++)
				tfa[z2] *= 100.0 / tf;
			obs->m_pHistogram->NormBinSum(100.0);
			obs->m_pHistogram->Integrate(false,1.0);
			a = OpenFileWrite((const char*)buf,true);
			mfprintf(a,"# Order Parameter;  Occurrence;  Integral;  Corrected Occurrence;  Corrected Integral\n");
			tf = 0;
			tf2 = 0;
			for (z2=0;z2<obs->m_pHistogram->m_iResolution;z2++) {
				tf += obs->m_pHistogram->m_pBin[z2];
				tf2 += tfa[z2];
				mfprintf(a,"%f;  %f;  %f;  %f;  %f\n",
					-0.5 + 1.5*(z2+0.5)/obs->m_pHistogram->m_iResolution,
					obs->m_pHistogram->m_pBin[z2],
					tf,
					tfa[z2],
					tf2);
			}
			fclose(a);*/
		}
		mprintf("\n");
	}

	if (m_bCombined) {
		buf.sprintf("order_vector%d_combined.csv",m_iNumber+1);
		mprintf("      Creating combined data file %s ...\n",(const char*)buf);
		a = OpenFileWrite((const char*)buf,true);
		mfprintf(a,"# Vector Set;  Average Order Parameter\n");
		for (z=0;z<(int)m_oaObservations.size();z++)
			mfprintf(a,"%d;  %.6f\n",z+1,m_oaObservations[z]->m_fAverage);
		fclose(a);
	}

}


