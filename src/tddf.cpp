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

#include "tddf.h"
#include "globalvar.h"
#include "maintools.h"


const char *GetRevisionInfo_tddf(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_tddf() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



bool CTDDFObservationRDF::Parse() {

	int z;
	CxString buf, buf2;
	

	mprintf( WHITE, "    *** TDDF Observation %d: RDF ***\n\n", m_iIndex+1 );

	m_bInter = AskYesNo("      Shall this observation be intermolecular (y) or intramolecular (n)? [yes] ",true);

	if (!m_bInter) {
		eprintf("Intramolecular RDFs not yet supported in TDDF module.\n");
		abort();
	}

	mprintf("\n");

	for (z=0;z<g_oaMolecules.GetSize();z++) {
		buf2.sprintf("%s=%d",(const char*)((CMolecule*)g_oaMolecules[z])->m_sName,z+1);
		buf.strcat(buf2);
		if (z < g_oaMolecules.GetSize()-1)
			buf.strcat(", ");
	}

	m_iRMIndex = AskRangeInteger_ND("      Which of the molecules should be the reference molecule (%s)? ",1,g_oaMolecules.GetSize(),(const char*)buf) - 1;

	m_iOMIndex = AskRangeInteger_ND("      Which of the molecules should be the observed molecule (%s)? ",1,g_oaMolecules.GetSize(),(const char*)buf) - 1;

	mprintf("\n");

	if (((CMolecule*)g_oaMolecules[m_iRMIndex])->m_iAtomGesNoVirt == 1) {
		mprintf("      %s is only one atom, there is no choice.\n",(const char*)((CMolecule*)g_oaMolecules[m_iRMIndex])->m_sName);
		m_oRMAG.Reset();
		m_oRMAG.m_pMolecule = (CMolecule*)g_oaMolecules[m_iRMIndex];
		m_oRMAG.AddAtom(0,0,false);
		m_oRMAG.SortAtoms();
		m_oRMAG.BuildName();
	} else {
_rma:
		m_oRMAG.Reset();
		AskString("      Which atom(s) to take from RM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ", &buf, "#2", (const char*)((CMolecule*)g_oaMolecules[m_iRMIndex])->m_sName );
		if (!m_oRMAG.ParseAtoms((CMolecule*)g_oaMolecules[m_iRMIndex],buf))
			goto _rma;
	}

	if (((CMolecule*)g_oaMolecules[m_iOMIndex])->m_iAtomGesNoVirt == 1) {
		mprintf("      %s is only one atom, there is no choice.\n",(const char*)((CMolecule*)g_oaMolecules[m_iOMIndex])->m_sName);
		m_oOMAG.Reset();
		m_oOMAG.m_pMolecule = (CMolecule*)g_oaMolecules[m_iOMIndex];
		m_oOMAG.AddAtom(0,0,false);
		m_oOMAG.SortAtoms();
		m_oOMAG.BuildName();
	} else {
_oma:
		m_oOMAG.Reset();
		AskString("      Which atom(s) to take from OM %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [#2] ", &buf, "#2", (const char*)((CMolecule*)g_oaMolecules[m_iOMIndex])->m_sName );
		if (!m_oOMAG.ParseAtoms((CMolecule*)g_oaMolecules[m_iOMIndex],buf))
			goto _oma;
	}

	mprintf("\n");

	m_fMinDist = AskFloat("      Enter minimal RDF distance (in pm): [0] ",0);
	m_fMaxDist = AskFloat("      Enter maximal RDF distance (in pm): [%d] ",(double)HalfBox(),HalfBox());

	mprintf("\n");

	m_bCorrectRadial = AskYesNo("      Correct radial distribution of this RDF (y/n)? [yes] ",true);

	if (m_bCorrectRadial)
		m_bUniform = AskYesNo("      Plot this RDF relative to uniform density (y) or as particle density (n)? [yes] ",true);
	else {
		m_bUniform = false;
		mprintf("\n      Uniform density is only possible with radial distribution correction.\n");
	}

	mprintf("\n");

	m_iResX = AskUnsignedInteger("      Enter resolution in distance dimension: [100] ",100);
	m_iResY = AskUnsignedInteger("      Enter resolution in temporal dimension: [100] ",100);

	mprintf("\n");

	if (g_iTrajSteps * g_fTimestepLength >= 3.0e6)
		m_bNS = AskYesNo("      Use nanoseconds (y) or picoseconds (n) for temporal axis? [yes] ",true);
	else
		m_bNS = AskYesNo("      Use nanoseconds (y) or picoseconds (n) for temporal axis? [no] ",false);

	m_fOverlap = AskFloat("      Enter overlap factor (0=no overlap, 1=one point in two slices, etc.): [2.0] ",2.0);

	mprintf( WHITE, "\n    *** TDDF Observation %d done ***\n\n", m_iIndex+1 );

	return true;
}


bool CTDDFObservationRDF::Initialize() {

	int z;


	mprintf("    Observation %d...\n",m_iIndex+1);

	m_pHistogram = new C2DF();

	m_pHistogram->m_iRes[0] = m_iResX;
	m_pHistogram->m_iRes[1] = m_iResY;
	m_pHistogram->m_fMinVal[0] = m_fMinDist;
	m_pHistogram->m_fMaxVal[0] = m_fMaxDist;
	m_pHistogram->m_fMinVal[1] = 0;
	if (m_bNS)
		m_pHistogram->m_fMaxVal[1] = g_iTrajSteps * g_fTimestepLength / 1000000.0;
	else
		m_pHistogram->m_fMaxVal[1] = g_iTrajSteps * g_fTimestepLength / 1000.0;
	m_pHistogram->Create();
	m_pHistogram->SetLabelX("Distance / pm");
	if (m_bNS)
		m_pHistogram->SetLabelY("Simulation Time / ns");
	else
		m_pHistogram->SetLabelY("Simulation Time / ps");
	m_pHistogram->SetLabelZ("g(r)");

	m_pHistogram->m_fPlotExp = 1.0;

	m_iaSliceCounter.resize(m_iResY);
	for (z=0;z<m_iResY;z++)
		m_iaSliceCounter[z] = 0;

	m_sName = ((CMolecule*)g_oaMolecules[m_iRMIndex])->m_sName;
	if (((CMolecule*)g_oaMolecules[m_iRMIndex])->m_iAtomGesNoVirt > 1) {
		m_sName += "_";
		m_sName += m_oRMAG.m_sName;
	}
	m_sName += "_";
	m_sName += ((CMolecule*)g_oaMolecules[m_iOMIndex])->m_sName;
	if (((CMolecule*)g_oaMolecules[m_iOMIndex])->m_iAtomGesNoVirt > 1) {
		m_sName += "_";
		m_sName += m_oOMAG.m_sName;
	}

	return true;
}


void CTDDFObservationRDF::ProcessStep(const CTimeStep *ts) {

	const CMolecule *rm, *om;
	const CSingleMolecule *rsm, *osm;
	int z1s, z2s, z1t, z2t, z1a, z2a, zt, ztmin, ztmax;
	double dist, tf;
	const CxIntArray *a1, *a2;


	tf = (((double)g_iSteps) / g_iMaxStep) * m_iResY;
	ztmin = (int)floor( tf - m_fOverlap/2.0 );
	ztmax = (int)floor( tf + m_fOverlap/2.0 );
	if (ztmin < 0)
		ztmin = 0;
	if (ztmax >= m_iResY)
		ztmax = m_iResY-1;

	for (zt=ztmin;zt<=ztmax;zt++)
		m_iaSliceCounter[zt]++;

	rm = (CMolecule*)g_oaMolecules[m_iRMIndex];
	om = (CMolecule*)g_oaMolecules[m_iOMIndex];

	for (z1s=0;z1s<rm->m_laSingleMolIndex.GetSize();z1s++) {

		rsm = (CSingleMolecule*)g_oaSingleMolecules[ rm->m_laSingleMolIndex[z1s] ];

		for (z2s=0;z2s<om->m_laSingleMolIndex.GetSize();z2s++) {

			if ((m_iRMIndex == m_iOMIndex) && (z1s == z2s))
				continue;

			osm = (CSingleMolecule*)g_oaSingleMolecules[ om->m_laSingleMolIndex[z2s] ];

			for (z1t=0;z1t<m_oRMAG.m_baAtomType.GetSize();z1t++) {

				a1 = (CxIntArray*)m_oRMAG.m_oaAtoms[z1t];

				for (z1a=0;z1a<a1->GetSize();z1a++) {

					for (z2t=0;z2t<m_oOMAG.m_baAtomType.GetSize();z2t++) {

						a2 = (CxIntArray*)m_oOMAG.m_oaAtoms[z2t];

						for (z2a=0;z2a<a2->GetSize();z2a++) {

							dist = FoldedLength(
								ts->m_vaCoords[ ((CxIntArray*)rsm->m_oaAtomOffset[m_oRMAG.m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)) ] -
								ts->m_vaCoords[ ((CxIntArray*)osm->m_oaAtomOffset[m_oOMAG.m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)) ]
							);

							for (zt=ztmin;zt<=ztmax;zt++)
								m_pHistogram->AddToBin_IntY( dist, zt );
						}
					}
				}
			}
		}
	}
}


void CTDDFObservationRDF::Finish() {

	CxString buf;
	int z, z2, ti;


	mprintf(WHITE,"    *** Observation %d ***\n",m_iIndex+1);

	mprintf("      %.0f bin entries, %.0f out of bin range (%.2f percent).\n",m_pHistogram->m_fBinEntries,m_pHistogram->m_fSkipEntries,ZeroDivide(m_pHistogram->m_fSkipEntries,m_pHistogram->m_fBinEntries+m_pHistogram->m_fSkipEntries)*100.0);
	m_pHistogram->CalcMaxEntry();
	mprintf("      Max. bin entry: %.6E  -->  EPS: %.6E\n",m_pHistogram->m_fMaxEntry,m_pHistogram->m_fEps);
	if (m_pHistogram->m_fEps > 1.0E-4) {
		eprintf("\n      Warning: Very large bin entries - probably loss of accuracy occurred.\n");
		mprintf("      Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
	}

	for (z=0;z<m_iResY;z++)
		for (z2=0;z2<m_iResX;z2++)
			if (m_iaSliceCounter[z] != 0)
				m_pHistogram->m_pBin[z*m_iResX+z2] /= (double)m_iaSliceCounter[z];

	if (m_bCorrectRadial) {

		mprintf("      Correcting radial distribution...\n");
		m_pHistogram->CorrectRadialDist(0);

		if (m_bUniform) {

			mprintf("      Normalizing RDF to uniform density...\n");

			if (m_iRMIndex == m_iOMIndex)
				ti = ((CMolecule*)g_oaMolecules[m_iOMIndex])->m_laSingleMolIndex.GetSize() - 1;
			else
				ti = ((CMolecule*)g_oaMolecules[m_iOMIndex])->m_laSingleMolIndex.GetSize();

			if (g_bBoxNonOrtho)
				m_pHistogram->MultiplyBin( g_fBoxVolume*1000000.0 / (4.0/3.0*Pi) /
					((CMolecule*)g_oaMolecules[m_iRMIndex])->m_laSingleMolIndex.GetSize() /
					m_oOMAG.m_iAtomGes / ti / m_oRMAG.m_iAtomGes );
			else
				m_pHistogram->MultiplyBin( g_fBoxX*g_fBoxY*g_fBoxZ / (4.0/3.0*Pi) /
					((CMolecule*)g_oaMolecules[m_iRMIndex])->m_laSingleMolIndex.GetSize() /
					m_oOMAG.m_iAtomGes / ti / m_oRMAG.m_iAtomGes );

		} else {

			mprintf("      Normalizing RDF to particle density (nm^-3)...\n");
			m_pHistogram->MultiplyBin( 1e9 / (4.0/3.0*Pi) /
				((CMolecule*)g_oaMolecules[m_iRMIndex])->m_laSingleMolIndex.GetSize() / m_oRMAG.m_iAtomGes );
		}

	} else {

		mprintf("      Normalizing bin integral to 10^6...\n");
		m_pHistogram->NormalizeBinIntegral( 1.0e6 );
	}

	if (g_bDoubleBox)
		m_pHistogram->MultiplyBin(g_iDoubleBoxFactor);

	buf.sprintf("tddf_obs%d_%s.nb",m_iIndex+1,(const char*)m_sName);
	mprintf("      Writing Mathematica Notebook to \"%s\"...\n",(const char*)buf);
	m_pHistogram->WriteMathematicaNb( "", (const char*)buf, "", false );
	buf.sprintf("tddf_obs%d_%s",m_iIndex+1,(const char*)m_sName);
	mprintf("      Writing Gnuplot Input to \"%s.gp\"...\n",(const char*)buf);
	m_pHistogram->WriteGnuplotInput( "", (const char*)buf, "", false );
}


bool CTDDFEngine::Parse() {

	mprintf( WHITE, "\n" );
	mprintf( WHITE, "    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
	mprintf( WHITE, "    >>>>>>      Time Dependent Distribution Functions      >>>>>>\n");
	mprintf( WHITE, "    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
	mprintf( WHITE, "\n" );

	g_bStrideParsed = true;
	g_iBeginStep = AskUnsignedInteger("    In which trajectory frame to start processing the trajectory? [1] ",1) - 1;
	if (g_iBeginStep == -1)
		g_iBeginStep = 0;
	if (g_iTrajSteps > 0)
		g_iMaxStep = AskUnsignedInteger("    How many trajectory frames to read (from this position on)? [%d] ",g_iTrajSteps-g_iBeginStep,g_iTrajSteps-g_iBeginStep);
	else
		g_iMaxStep = AskUnsignedInteger_ND("    How many trajectory frames to read (from this position on)? ");
	g_iStride = AskUnsignedInteger("    Use every n-th read trajectory frame for the analysis: [1] ",1);

	mprintf("\n");

	if (g_iTrajSteps * g_fTimestepLength >= 1.0e6)
		mprintf("    This corresponds to a physical time of %.3f ns.\n",g_iTrajSteps * g_fTimestepLength / 1.0e6 );
	else if (g_iTrajSteps * g_fTimestepLength >= 1.0e3)
		mprintf("    This corresponds to a physical time of %.3f ps.\n",g_iTrajSteps * g_fTimestepLength / 1.0e3 );
	else
		mprintf("    This corresponds to a physical time of %.1f fs.\n",g_iTrajSteps * g_fTimestepLength );

	do {

		mprintf("\n    Currently, only RDFs are supported.\n\n");

		m_oaObservations.push_back( new CTDDFObservationRDF((int)m_oaObservations.size()) );

		if (!m_oaObservations.back()->Parse())
			return false;

	} while (AskYesNo("    Add another TDDF observation (y/n)? [no] ",false));

	mprintf( WHITE, "\n" );
	mprintf( WHITE, "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
	mprintf( WHITE, "    <<<<<<      Time Dependent Distribution Functions      <<<<<<\n");
	mprintf( WHITE, "    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
	mprintf( WHITE, "\n" );

	return true;
}


bool CTDDFEngine::Initialize() {

	int z;

	mprintf("  Initializing TDDF Observations...\n");

	for (z=0;z<(int)m_oaObservations.size();z++)
		if (!m_oaObservations[z]->Initialize())
			return false;

	return true;
}


void CTDDFEngine::ProcessStep(const CTimeStep *ts) {

	int z;

	for (z=0;z<(int)m_oaObservations.size();z++)
		m_oaObservations[z]->ProcessStep(ts);
}


void CTDDFEngine::Finish() {

	int z;

	mprintf(WHITE,"\n    Finishing TDDF Observations...\n");

	for (z=0;z<(int)m_oaObservations.size();z++)
		m_oaObservations[z]->Finish();
}







