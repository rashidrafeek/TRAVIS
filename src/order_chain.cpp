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
#include "order_chain.h"
#include "globalvar.h"
#include "maintools.h"


const char *GetRevisionInfo_order_chain(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_order_chain() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



bool COrderAnalysis_Chain::Parse(CTimeStep *ts) {

	int z;
	//CMolecule *m;
	//CSingleMolecule *sm;
	CxString buf;


	UNUSED(ts);

	mprintf(WHITE,"\n    >>> Analysis %d: Alkyl Chain Order Parameters >>>\n\n",m_iNumber+1);

	if (g_oaMolecules.GetSize() == 1) {
		mprintf("      Using molecule 1.\n");
		m_iMolecule = 0;
	} else
		m_iMolecule = AskRangeInteger_ND("      Define vectors in which molecule (1-%d)? ",1,g_oaMolecules.GetSize(),g_oaMolecules.GetSize()) - 1;
	if (m_pParent->m_oaChains[m_iMolecule].size() == 0) {
		eprintf("\n      Error: No alkyl chains have been detected in molecule %d.\n",m_iMolecule+1);
		abort();
	}
	//m = (CMolecule*)g_oaMolecules[m_iMolecule];
_ilagain:
	m_iaChains.clear();
	AskString_ND("      Which chain(s) to use in molecule %d (comma separated, see list above)? ",&buf,m_iMolecule+1);
	if (!ParseIntList((const char*)buf,m_iaChains)) {
		eprintf("\n      Invalid input.\n\n");
		goto _ilagain;
	}
	for (z=0;z<(int)m_iaChains.size();z++) {
		if ((m_iaChains[z] < 1) || (m_iaChains[z] > (int)m_pParent->m_oaChains[m_iMolecule].size())) {
			eprintf("\n      Chain number %d outside of allowed range (%d .. %lu).\n\n",m_iaChains[z],1,(unsigned long)m_pParent->m_oaChains[m_iMolecule].size());
			goto _ilagain;
		}
		m_iaChains[z]--;
		if (z > 0) {
			if (m_pParent->m_oaChains[m_iMolecule][m_iaChains[0]]->m_iaAtoms.size() != m_pParent->m_oaChains[m_iMolecule][m_iaChains[z]]->m_iaAtoms.size()) {
				eprintf("\n      Chains have different lengths (%d: %lu vs. %d: %lu).\n\n",m_iaChains[0]+1,(unsigned long)m_pParent->m_oaChains[m_iMolecule][m_iaChains[0]]->m_iaAtoms.size(),m_iaChains[z]+1,(unsigned long)m_pParent->m_oaChains[m_iMolecule][m_iaChains[z]]->m_iaAtoms.size());
				goto _ilagain;
			}
		} else {
			if (m_pParent->m_oaChains[m_iMolecule][m_iaChains[0]]->m_iaAtoms.size() < 4) {
				eprintf("\n      Chains need to have at least length 4 (found: %lu).\n\n",(unsigned long)m_pParent->m_oaChains[m_iMolecule][m_iaChains[0]]->m_iaAtoms.size());
				goto _ilagain;
			}
		}
	}
	mprintf("\n");
	m_iDihedrals = (int)m_pParent->m_oaChains[m_iMolecule][m_iaChains[0]]->m_iaAtoms.size()-3;
	mprintf("      Using %lu chains of length %lu, therefore %d dihedrals.\n\n",(unsigned long)m_iaChains.size(),(unsigned long)m_pParent->m_oaChains[m_iMolecule][m_iaChains[0]]->m_iaAtoms.size(),m_iDihedrals);

	m_bSaveHistograms = AskYesNo("      Save histogram for each group of dihedrals (y/n)? [yes] ",true);

	m_faGauche.resize(m_iDihedrals);
	m_faTotal.resize(m_iDihedrals);
	m_faEcliptic1.resize(m_iDihedrals);
	m_faEcliptic2.resize(m_iDihedrals);

	for (z=0;z<m_iDihedrals;z++) {
		m_faGauche[z] = 0;
		m_faTotal[z] = 0;
		m_faEcliptic1[z] = 0;
		m_faEcliptic2[z] = 0;
	}

	m_fGaucheCount = 0;
	m_fGaucheSqSum = 0;
	m_fCountTotal = 0;

	m_pGaucheCountHistogram = new CDF();
	m_pGaucheCountHistogram->m_fMinVal = 0;
	m_pGaucheCountHistogram->m_fMaxVal = (double)m_iDihedrals+1.0;
	m_pGaucheCountHistogram->m_iResolution = (int)m_iDihedrals+1;
	m_pGaucheCountHistogram->Create();
	m_pGaucheCountHistogram->SetLabelX("Number of Gauche Dihedrals in Chain");
	m_pGaucheCountHistogram->SetLabelY("Occurrence");
	m_pGaucheCountHistogram->m_bLeft = true;

	if (m_bSaveHistograms) {
		m_oaHistograms.resize(m_iDihedrals);
		for (z=0;z<m_iDihedrals;z++) {
			m_oaHistograms[z] = new CDF();
			m_oaHistograms[z]->m_fMinVal = 0;
			m_oaHistograms[z]->m_fMaxVal = 180.0;
			m_oaHistograms[z]->m_iResolution = 100;
			m_oaHistograms[z]->Create();
			m_oaHistograms[z]->SetLabelX("Dihedral Angle");
			m_oaHistograms[z]->SetLabelY("Occurrence");
		}
	}
	
	mprintf(WHITE,"\n    <<< Analysis %d: Alkyl Chain Order Parameters <<<\n\n",m_iNumber+1);

	return true;
}


void COrderAnalysis_Chain::ProcessStep(CTimeStep *ts) {

	int z, z2, z3, i;
	CMolecule *m;
	CSingleMolecule *sm;
	COrderChain *oc;
	double tf;


	m = (CMolecule*)g_oaMolecules[m_iMolecule];

	for (z=0;z<m->m_laSingleMolIndex.GetSize();z++) {

		sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z]];

		for (z2=0;z2<(int)m_iaChains.size();z2++) {

			oc = m_pParent->m_oaChains[m_iMolecule][m_iaChains[z2]];

			i = 0;

			for (z3=0;z3<m_iDihedrals;z3++) {

				tf = Dihedral(
					FoldVector(ts->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[oc->m_iCarbonIndex])->GetAt(oc->m_iaAtoms[z3])]
					- ts->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[oc->m_iCarbonIndex])->GetAt(oc->m_iaAtoms[z3+1])]),
					FoldVector(ts->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[oc->m_iCarbonIndex])->GetAt(oc->m_iaAtoms[z3+3])]
					- ts->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[oc->m_iCarbonIndex])->GetAt(oc->m_iaAtoms[z3+2])]),
					FoldVector(ts->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[oc->m_iCarbonIndex])->GetAt(oc->m_iaAtoms[z3+1])]
					- ts->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[oc->m_iCarbonIndex])->GetAt(oc->m_iaAtoms[z3+2])]),
					true
					);

				if (m_bSaveHistograms)
					m_oaHistograms[z3]->AddToBin(tf);

				if (tf < 120.0) {
					i++;
					m_faGauche[z3]++;
				}

				if (tf < 20.0)
					m_faEcliptic1[z3]++;

				if ((tf > 100.0) && (tf < 140.0))
					m_faEcliptic2[z3]++;

				m_faTotal[z3]++;
			}

			m_fGaucheCount += (double)i/m_iDihedrals;
			m_fGaucheSqSum += (double)i*i/m_iDihedrals/m_iDihedrals;
			m_fCountTotal++;
			m_pGaucheCountHistogram->AddToBin((double)i);
		}
	}
}


void COrderAnalysis_Chain::Finish() {

	int z;
	double mean, sd;
	CxString buf;
	FILE *a;


	mprintf(WHITE,"\n    *** Chain Order Parameter Analysis %d ***\n\n",m_iNumber+1);

	for (z=0;z<m_iDihedrals;z++) {
		mprintf("      Dihedral %2d:  %5.2f%c Gauche,  %5.2f%c Ecliptic 1,  %5.2f%c Ecliptic 2,  %.0f hits.\n",z+1,m_faGauche[z]/m_faTotal[z]*100.0,'%',m_faEcliptic1[z]/m_faTotal[z]*100.0,'%',m_faEcliptic2[z]/m_faTotal[z]*100.0,'%',m_faTotal[z]);
		if (m_bSaveHistograms) {
			buf.sprintf("order_chain%d_histogram_%d.csv",m_iNumber+1,z+1);
			mprintf("        Saving histogram to %s ...\n",(const char*)buf);
			m_oaHistograms[z]->NormBinSum(100.0);
			m_oaHistograms[z]->Write((const char*)buf,"","",false);
		}
	}

	mprintf("\n");
	buf.sprintf("order_chain%d_results.csv",m_iNumber+1);
	mprintf("      Saving results to %s ...\n",(const char*)buf);
	a = OpenFileWrite((const char*)buf,true);
	mfprintf(a,"# Dihedral;  Gauche percentage;  Anti percentage;  Ecliptic 1 percentage;  Ecliptic 2 percentage\n");

	for (z=0;z<m_iDihedrals;z++)
		mfprintf(a,"%d;  %.6f;  %.6f;  %.6f;  %.6f\n",
			z+1,
			m_faGauche[z]/m_faTotal[z]*100.0,
			(1.0-m_faGauche[z]/m_faTotal[z])*100.0,
			m_faEcliptic1[z]/m_faTotal[z]*100.0,
			m_faEcliptic2[z]/m_faTotal[z]*100.0
			);

	fclose(a);
	mprintf("\n");

	mean = m_fGaucheCount/m_fCountTotal;
	sd = sqrt(m_fGaucheSqSum/m_fCountTotal - mean*mean);
	mprintf("      Average Gauche count per chain is %.3f%c +/- %.3f%c.\n",mean*100.0,'%',100.0*sd,'%');
	mprintf("        (Uncertainty is the standard deviation)\n\n");
	buf.sprintf("order_chain%d_gauchecount.csv",m_iNumber+1);
	mprintf("      Saving Gauche count histogram to %s ...\n",(const char*)buf);
	m_pGaucheCountHistogram->NormBinSum(100.0);
	m_pGaucheCountHistogram->Write((const char*)buf,"","",false);
	mprintf("\n");
}



