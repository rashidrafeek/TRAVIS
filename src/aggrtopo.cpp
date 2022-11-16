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

#include "aggrtopo.h"
#include "tools.h"
#include "timestep.h"
#include "globalvar.h"
#include "df.h"
#include "maintools.h"



const char *GetRevisionInfo_aggrtopo(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}



const char *GetSourceVersion_aggrtopo() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



bool CAggrTopoGroup::Parse() {

	int z, z2, z3, z4;
	CAtomGroup *ag;
	CMolecule *m;
	CSingleMolecule *sm;
	CxString buf, buflong, bufshort;
	bool b;


	mprintf(WHITE,"\n    >>> Aggregation Topology Group %d >>>\n\n",m_iObservation+1);

	m_bDonor = AskYesNo("      Is this a donor (y) or acceptor (n) group? [yes] ",true);

	mprintf("\n");

	buflong = "";
	bufshort = "";
	b = false;

	m_oaAtoms.resize( g_oaMolecules.GetSize() );

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];

		if (AskYesNo("      Use atoms from %s for this group (y/n)? [no] ",false,m->m_sName)) {

			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);
_rowagain:
			AskString("        Which atoms from %s to use (e.g. \"C1,C5-7\", *=all)? [all] ",&buf,"*",m->m_sName);
			if (!ag->ParseAtoms(m,(const char*)buf)) {
				eprintf("Invalid input.\n");
				goto _rowagain;
			}
			m_oaAtoms[z] = ag;

			for (z2=0;z2<ag->m_oaAtoms.GetSize();z2++) {
				for (z3=0;z3<((CxIntArray*)ag->m_oaAtoms[z2])->GetSize();z3++) {
					for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++) {
						sm = (CSingleMolecule*)g_oaSingleMolecules[ m->m_laSingleMolIndex[z4] ];
						m_iaAtoms.push_back( ((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z3)) );
					}
				}
			}

			if (b) {
				bufshort += "_";
				buflong += ", ";
			} else
				b = true;

			bufshort += m->m_sName;
			bufshort += "_";
			bufshort += ag->m_sName;

			buflong += ag->m_sName;
			buflong += " from ";
			buflong += m->m_sName;

		} else
			m_oaAtoms[z] = NULL;
	}

	m_sShortName = (const char*)bufshort;
	m_sLongName = (const char*)buflong;

	mprintf("\n      This yields %lu atoms.\n\n",(unsigned long)m_iaAtoms.size());

	mprintf(WHITE,"    <<< Aggregation Topology Group %d <<<\n\n",m_iObservation+1);

	return true;
}



bool CAggrTopoEngine::Parse() {

	CAggrTopoGroup *g;
	int z;


	mprintf(WHITE,"\n    >>> Aggregation Topology Analysis >>>\n\n");

	m_fMaxDist = AskFloat("    Enter maximal observation distance (in pm): [500] ",500.0);
	m_iResolution = AskUnsignedInteger("    Enter RDF binning resolution: [150] ",150);

	mprintf("\n");
	m_bMinMax = AskYesNo("    Use first RDF minimum as cutoff (y), or enter fixed geometric criterion (n)? [n] ",false);

	if (!m_bMinMax) {

		mprintf("\n");
		m_fCutoffDist = AskFloat("    Enter cutoff distance for interactions (in pm): [300] ",300.0);
	}

	//m_bDomDec = AskYesNo("    Use domain decomposition to speed up the analysis (y/n)? [yes] ",true);

	while (true) {

		try { g = new CAggrTopoGroup(); } catch(...) { g = NULL; }
		if (g == NULL) NewException((double)sizeof(CAggrTopoGroup),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		g->m_iObservation = (int)m_oaGroups.size();
		m_oaGroups.push_back(g);

		if (!g->Parse())
			return false;

		if (!AskYesNo("    Add another aggregation topology group (y/n)? [no] ",false))
			break;
	}

	for (z=0;z<(int)m_oaGroups.size();z++)
		if (m_oaGroups[z]->m_bDonor)
			m_oaDonorGroups.push_back( m_oaGroups[z] );
		else
			m_oaAcceptorGroups.push_back( m_oaGroups[z] );

	mprintf("\n    Have %lu donor groups and %lu acceptor groups:\n\n", (unsigned long)m_oaDonorGroups.size(), (unsigned long)m_oaAcceptorGroups.size() );

	for (z=0;z<(int)m_oaGroups.size();z++)
		mprintf(
			"      Group %2d %s %5lu atoms:  %s.\n",
			z+1,
			m_oaGroups[z]->m_bDonor?"(donor),   ":"(acceptor),",
			(unsigned long)m_oaGroups[z]->m_iaAtoms.size(),
			(const char*)m_oaGroups[z]->m_sLongName.c_str()
		);

	mprintf("\n");
	mprintf("    This yields %lu combinations for interactions.\n", (unsigned long)m_oaDonorGroups.size() * (unsigned long)m_oaAcceptorGroups.size() );

	mprintf("\n");
	mprintf("    You can now decide for each molecule type if you want to include intramolecular interactions.\n");
	mprintf("\n");
	mprintf("    Note: Interactions between atoms with less than 4 covalent bonds distance (i.e., 1-2, 1-3, and 1-4\n");
	mprintf("          interactions) are automatically excluded, even if you answer \"yes\" below.\n");
	mprintf("\n");

	m_iaIntra.resize(g_oaMolecules.GetSize());
	for (z=0;z<g_oaMolecules.GetSize();z++)
		m_iaIntra[z] = AskYesNo("    Include intramolecular interactions in molecule type %d (%s) (y/n)? [no] ",false,z+1,((CMolecule*)g_oaMolecules[z])->m_sName) ? 1 : 0;

	mprintf(WHITE,"\n    <<< Aggregation Topology Analysis <<<\n\n");

	return true;
}



void CAggrTopoEngine::Initialize() {

	int zd, za;
	CAggrTopoGroup *gd, *ga;
	CDF *df;

	mprintf("  Initializing Aggregation Topology Analysis...\n");

	BuildNeighborMatrix();

	for (zd=0;zd<(int)m_oaDonorGroups.size();zd++) {

		gd = m_oaDonorGroups[zd];

		gd->m_oaRDFs.resize( m_oaAcceptorGroups.size() );
		gd->m_oaObsGroups.resize( m_oaAcceptorGroups.size() );
		
		for (za=0;za<(int)m_oaAcceptorGroups.size();za++) {
			try { df = new CDF(); } catch(...) { df = NULL; }
			if (df == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			df->m_iResolution = m_iResolution;
			df->m_fMinVal = 0;
			df->m_fMaxVal = m_fMaxDist;
			df->SetLabelX("Distance / pm");
			df->SetLabelY("g(r)");
			df->Create();
			gd->m_oaRDFs[za] = df;
			gd->m_oaObsGroups[za] = m_oaAcceptorGroups[za];
		}
	}


	for (za=0;za<(int)m_oaAcceptorGroups.size();za++) {

		ga = m_oaAcceptorGroups[za];

		ga->m_oaRDFs.resize( m_oaDonorGroups.size() );
		ga->m_oaObsGroups.resize( m_oaDonorGroups.size() );
		
		for (zd=0;zd<(int)m_oaDonorGroups.size();zd++) {
			try { df = new CDF(); } catch(...) { df = NULL; }
			if (df == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			df->m_iResolution = m_iResolution;
			df->m_fMinVal = 0;
			df->m_fMaxVal = m_fMaxDist;
			df->SetLabelX("Distance / pm");
			df->SetLabelY("g(r)");
			df->Create();
			ga->m_oaRDFs[zd] = df;
			ga->m_oaObsGroups[zd] = m_oaDonorGroups[zd];
		}
	}
}



void CAggrTopoEngine::ProcessStep(CTimeStep *ts) {

	int zd, za, zd2, za2, ia, id;
	double tf;
	CAggrTopoGroup *gd, *ga;
	CDF *dfd, *dfa;


	for (zd=0;zd<(int)m_oaDonorGroups.size();zd++) {

		gd = m_oaDonorGroups[zd];

		for (za=0;za<(int)m_oaAcceptorGroups.size();za++) {

			ga = m_oaAcceptorGroups[za];

			dfd = gd->m_oaRDFs[za];
			dfa = ga->m_oaRDFs[zd];

			for (zd2=0;zd2<(int)gd->m_iaAtoms.size();zd2++) {
				id = gd->m_iaAtoms[zd2];
				for (za2=0;za2<(int)ga->m_iaAtoms.size();za2++) {
					ia = ga->m_iaAtoms[za2];
					if (g_laAtomSMIndex[id] == g_laAtomSMIndex[ia]) {
						if (m_iaIntra[ g_waAtomMolIndex[id] ] == 0)
							continue;
						if (m_pNeighborMatrix[id*g_iGesAtomCount+ia] <= 4)
							continue;
					}
					tf = FoldedLength( ts->m_vaCoords[ id ] - ts->m_vaCoords[ ia ] );
					dfd->AddToBin( tf );
					dfa->AddToBin( tf );
				}
			}
		}
	}
}



void CAggrTopoEngine::BuildNeighborMatrix() {

	int z, z2, z3, z4, z5, z2s;
	int o1, o2, o3, o4;
	CMolAtom *ma;
	CMolecule *m;
	CSingleMolecule *sm;


	mprintf("    Neighbor matrix occupies %s of RAM.\n",FormatBytes((double)g_iGesAtomCount*g_iGesAtomCount*sizeof(unsigned char)));

	try { m_pNeighborMatrix = new unsigned char[g_iGesAtomCount*g_iGesAtomCount]; } catch(...) { m_pNeighborMatrix = NULL; }
	if (m_pNeighborMatrix == NULL) NewException((double)sizeof(unsigned char)*g_iGesAtomCount*g_iGesAtomCount,__FILE__,__LINE__,__PRETTY_FUNCTION__);

	for (z2=0;z2<g_iGesAtomCount;z2++)
		for (z3=0;z3<g_iGesAtomCount;z3++)
			if (z2 == z3)
				m_pNeighborMatrix[z2*g_iGesAtomCount+z3] = 0;
			else
				m_pNeighborMatrix[z2*g_iGesAtomCount+z3] = 9;

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];

		for (z2s=0;z2s<m->m_laSingleMolIndex.GetSize();z2s++) {

			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2s]];

			for (z2=0;z2<sm->m_oaMolAtoms.GetSize();z2++) {

				ma = (CMolAtom*)sm->m_oaMolAtoms[z2];
				o1 = ma->m_iOffset;

				for (z3=0;z3<sm->m_oaBonds.GetSize();z3++) {

					if (((CMolBond*)sm->m_oaBonds[z3])->m_iAtomOffset[0] == o1)
						o2 = ((CMolBond*)sm->m_oaBonds[z3])->m_iAtomOffset[1];
					else if (((CMolBond*)sm->m_oaBonds[z3])->m_iAtomOffset[1] == o1)
						o2 = ((CMolBond*)sm->m_oaBonds[z3])->m_iAtomOffset[0];
					else
						continue;

					if (m_pNeighborMatrix[o1*g_iGesAtomCount+o2] > 1)
						m_pNeighborMatrix[o1*g_iGesAtomCount+o2] = 1;

					for (z4=0;z4<sm->m_oaBonds.GetSize();z4++) {

						if (((CMolBond*)sm->m_oaBonds[z4])->m_iAtomOffset[0] == o2)
							o3 = ((CMolBond*)sm->m_oaBonds[z4])->m_iAtomOffset[1];
						else if (((CMolBond*)sm->m_oaBonds[z4])->m_iAtomOffset[1] == o2)
							o3 = ((CMolBond*)sm->m_oaBonds[z4])->m_iAtomOffset[0];
						else
							continue;

						if (o3 == o1)
							continue;

						if (m_pNeighborMatrix[o1*g_iGesAtomCount+o3] > 2)
							m_pNeighborMatrix[o1*g_iGesAtomCount+o3] = 2;

						for (z5=0;z5<sm->m_oaBonds.GetSize();z5++) {

							if (((CMolBond*)sm->m_oaBonds[z5])->m_iAtomOffset[0] == o3)
								o4 = ((CMolBond*)sm->m_oaBonds[z5])->m_iAtomOffset[1];
							else if (((CMolBond*)sm->m_oaBonds[z5])->m_iAtomOffset[1] == o3)
								o4 = ((CMolBond*)sm->m_oaBonds[z5])->m_iAtomOffset[0];
							else
								continue;

							if ((o4 == o1) || (o4 == o2))
								continue;

							if (m_pNeighborMatrix[o1*g_iGesAtomCount+o4] > 3)
								m_pNeighborMatrix[o1*g_iGesAtomCount+o4] = 3;
						}
					}
				}
			}
		}
	}
}



void CAggrTopoEngine::Finish() {

	int z, z2, z3;
	double tf;
	FILE *a;
	CAggrTopoGroup *g1, *g2;
	CxString buf;
	CDF *df;


	mprintf(WHITE,"\n    Finishing Aggregation Topology Analysis...\n");

	for (z=0;z<(int)m_oaGroups.size();z++) {

		g1 = m_oaGroups[z];

		g1->m_faAggrValues.resize( g1->m_oaObsGroups.size() );

		for (z2=0;z2<(int)g1->m_oaObsGroups.size();z2++) {

			g2 = g1->m_oaObsGroups[z2];
			df = g1->m_oaRDFs[z2];

			mprintf("\n");
			mprintf(WHITE,"      * RDF:  %s  ---  %s\n",g1->m_sLongName.c_str(),g2->m_sLongName.c_str());
			mprintf("        %.0f bin entries, %.0f out of bin range (%.2f percent).\n",df->m_fBinEntries,df->m_fSkipEntries,ZeroDivide(df->m_fSkipEntries,df->m_fBinEntries+df->m_fSkipEntries)*100.0);

			df->CalcMinMax();
			mprintf("        Max. bin entry: %.6E  -->  EPS: %.6E\n",df->m_fMaxEntry,df->m_fEps);
			if (df->m_fEps > 1.0E-4) {
				eprintf("\n        Warning: Very large bin entries - probably loss of accuracy occurred.\n");
				mprintf("        Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
			}
			df->CalcMeanSD();
			mprintf("        Mean value: %10G pm    Standard deviation: %10G pm\n",df->m_fMean,df->m_fSD);
			mprintf("        Min. value: %10G pm    Max.value:          %10G pm\n",df->m_fMinInput,df->m_fMaxInput);

			df->MultiplyBin( 1.0 / g_iSteps * g_iStride );

			df->CorrectRadialDist();

			if (g_bBoxNonOrtho) {
				df->Integrate( true, 1.0 / (double)g1->m_iaAtoms.size() );
				df->MultiplyBin( g_fBoxVolume*1000000.0 / (4.0/3.0*Pi) / (double)g1->m_iaAtoms.size() / (double)g2->m_iaAtoms.size() );
			} else {
				df->Integrate( true, 1.0 / (double)g1->m_iaAtoms.size() );
				df->MultiplyBin( g_fBoxX*g_fBoxY*g_fBoxZ / (4.0/3.0*Pi) / (double)g1->m_iaAtoms.size() / (double)g2->m_iaAtoms.size() );
			}

			if (g_bDoubleBox) {
				df->MultiplyBin( g_iDoubleBoxFactor );
				df->MultiplyIntegral( g_iDoubleBoxFactor );
			}

			df->FindFirstMinMax(false);
			if (df->m_faMaxPos.size() != 0) {
				mprintf(WHITE,"        >>> Peak Data >>>\n");
				for (z3=0;z3<(int)df->m_faMaxPos.size();z3++) {
					mprintf("          %d. maximum at %7.2f pm, height %7.4f\n",z3+1,df->m_faMaxPos[z3],df->m_faMaxVal[z3]);
					if ((int)df->m_faMinPos.size() > z3) {
						mprintf("          %d. minimum at %7.2f pm, height %7.4f",z3+1,df->m_faMinPos[z3],df->m_faMinVal[z3]);
						if ((int)df->m_faShellInt.size() > z3)
							mprintf(", %7.4f particles in %d. solvation shell.",df->m_faShellInt[z3],z3+1);
						mprintf("\n");
					}
				}
				mprintf(WHITE,"        <<< Peak Data <<<\n");
			}

			buf.sprintf("aggrtopo_rdf_%s___%s.csv",g1->m_sShortName.c_str(),g2->m_sShortName.c_str());
			mprintf("        Saving RDF as \"%s\"...\n",(const char*)buf);
			df->Write("",buf,"",true);

			if (m_bMinMax) {

				if (df->m_faShellInt.size() > 0)
					g1->m_faAggrValues[z2] = df->m_faShellInt[0];
				else
					g1->m_faAggrValues[z2] = 0;

				mprintf("        Using first shell integral as result: %10.6f\n",g1->m_faAggrValues[z2]);

			} else {

				g1->m_faAggrValues[z2] = df->m_pIntegral[ (int)((m_fCutoffDist-df->m_fMinVal)*df->m_fFac) ];

				mprintf("        Using integral up to %.2f pm as result: %10.6f\n",m_fCutoffDist,g1->m_faAggrValues[z2]);
			}
		}
	}
	mprintf("\n\n");
	mprintf(WHITE,"    *** Result Matrix ***\n\n");
	mprintf("      (rows = reference groups, columns = observed groups)\n\n");

	mprintf("              |");
	for (z=0;z<(int)m_oaGroups.size();z++)
		mprintf("    Group %-2d ",z+1);
	mprintf("\n");
	mprintf("    ----------|");
	for (z=0;z<(int)m_oaGroups.size();z++)
		mprintf("-------------");
	mprintf("---\n");

	for (z=0;z<(int)m_oaGroups.size();z++) {
		g1 = m_oaGroups[z];
		mprintf("     Group %-2d |",z+1);
		for (z2=0;z2<(int)m_oaGroups.size();z2++) {
			for (z3=0;z3<(int)g1->m_oaObsGroups.size();z3++) {
				if (g1->m_oaObsGroups[z3] == m_oaGroups[z2]) {
					mprintf("   %8.5f  ",g1->m_faAggrValues[z3]);
					goto _found;
				}
			}
			mprintf("        ---  ");
_found:;
		}
		mprintf("\n");
	}

	mprintf("\n\n");

	for (z=0;z<(int)m_oaGroups.size();z++)
		mprintf(
			"      Group %2d %s %5lu atoms:  %s.\n",
			z+1,
			m_oaGroups[z]->m_bDonor?"(donor),   ":"(acceptor),",
			(unsigned long)m_oaGroups[z]->m_iaAtoms.size(),
			(const char*)m_oaGroups[z]->m_sLongName.c_str()
		);

	mprintf("\n");

	mprintf("    Writing text data file to \"aggrtopo.dat\"...\n");

	a = OpenFileWrite( "aggrtopo.dat", true );

	struct tm *today;
	time_t ltime;
	char cbuf[256];

	time(&ltime);
	today = localtime(&ltime);
	strcpy(cbuf,asctime(today));
	cbuf[strlen(cbuf)-1] = 0;


	mfprintf(a,"# Data file from aggregation topology analysis\n");
	mfprintf(a,"# Written by TRAVIS at %s\n",cbuf);
	mfprintf(a,"# See http://www.travis-analyzer.de/\n");
	mfprintf(a,"# Number of groups\n");
	mfprintf(a,"  %lu\n",(unsigned long)m_oaGroups.size());
	for (z=0;z<(int)m_oaGroups.size();z++) {
		g1 = m_oaGroups[z];
		mfprintf(a,"# Group %d Description\n",z+1);
		mfprintf(a,"  %s%s\n",g1->m_bDonor?"Donor, ":"Acceptor, ",g1->m_sLongName.c_str());
		mfprintf(a,"# Outgoing connections to Groups 1-%d, one per line\n",(int)m_oaGroups.size());
		tf = 0;
		for (z2=0;z2<(int)m_oaGroups.size();z2++) {
			for (z3=0;z3<(int)g1->m_oaObsGroups.size();z3++) {
				if (g1->m_oaObsGroups[z3] == m_oaGroups[z2]) {
					tf += g1->m_faAggrValues[z3];
					goto _found2a;
				}
			}
_found2a:;
		}
		mfprintf(a,"# Outgoing sum:  %14.10f\n",tf);
		for (z2=0;z2<(int)m_oaGroups.size();z2++) {
			for (z3=0;z3<(int)g1->m_oaObsGroups.size();z3++) {
				if (g1->m_oaObsGroups[z3] == m_oaGroups[z2]) {
					mfprintf(a,"%14.10f\n",g1->m_faAggrValues[z3]);
					goto _found2;
				}
			}
			mfprintf(a,"  0\n");
_found2:;
		}
	}
	mfprintf(a,"# End of data file\n");
	fclose(a);

	mprintf("    To create a Sankey diagram from the result file, use the following command:\n\n");
	mprintf("      travis -sankey aggrtopo.dat\n");

	mprintf("\n");
}




