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

#include "contactmatrix.h"
#include "tools.h"
#include "globalvar.h"
#include "matrixplot.h"
#include "maintools.h"


const char *GetRevisionInfo_contactmatrix(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_contactmatrix() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



CContactMatrixObservation::CContactMatrixObservation() {

	m_oaMatrix = NULL;
	m_iRows = 0;
	m_iCols = 0;
}


CContactMatrixObservation::~CContactMatrixObservation() {
}


bool CContactMatrixObservation::Parse() {

	int z, z2, z3, z4;
	CAtomGroup *ag;
	CMolecule *m;
	CSingleMolecule *sm;
	CMolAtom *ma, *ma0;
	CxString buf, ts;


	mprintf(WHITE,"\n      >>> Contact Matrix Observation %d >>>\n\n",m_iObservation+1);

	m_bInter = AskYesNo("      Take into account intermolecular contacts (y/n)? [yes] ",true);

	if (m_bInter)
		m_bIntra = AskYesNo("      Take into account intramolecular contacts (y/n)? [no] ",false);
	else
		m_bIntra = AskYesNo("      Take into account intramolecular contacts (y/n)? [yes] ",true);

	if (m_bIntra) {
		m_bExclude12 = AskYesNo("        Exclude 1,2-bonded contacts (y/n)? [yes] ",true);
		m_bExclude13 = AskYesNo("        Exclude 1,3-bonded contacts (y/n)? [yes] ",true);
		m_bExclude14 = AskYesNo("        Exclude 1,4-bonded contacts (y/n)? [yes] ",true);
	}

	mprintf("\n");

	m_fMaxDist = AskFloat("      Enter maximum contact distance (in pm): [350] ",350.0);

	m_iResolution = AskUnsignedInteger("      Enter distance histogram resolution: [100] ",100);

	mprintf("\n");

	m_oaRows.resize(g_oaMolecules.GetSize());
	m_oaColumns.resize(g_oaMolecules.GetSize());

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];

		if (AskYesNo("      Use atoms from molecule type %s for the rows (y/n)? [yes] ",true,m->m_sName)) {

			ag = new CAtomGroup();
_rowagain:
			AskString("        Which atoms from %s to use (e.g. \"C1,C5-7\", *=all)? [all] ",&buf,"*",m->m_sName);
			if (!ag->ParseAtoms(m,(const char*)buf)) {
				eprintf("Invalid input.\n");
				goto _rowagain;
			}
			m_oaRows[z] = ag;
		} else
			m_oaRows[z] = NULL;
	}

	mprintf("\n");

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];

		if (AskYesNo("      Use atoms from molecule type %s for the columns (y/n)? [yes] ",true,m->m_sName)) {

			ag = new CAtomGroup();
_colagain:
			AskString("        Which atoms from %s to use (e.g. \"C1,C5-7\", *=all)? [all] ",&buf,"*",m->m_sName);
			if (!ag->ParseAtoms(m,(const char*)buf)) {
				eprintf("Invalid input.\n");
				goto _colagain;
			}
			m_oaColumns[z] = ag;
		} else
			m_oaColumns[z] = NULL;
	}

	mprintf("\n");

	m_bJoinEquivalent = AskYesNo("      Merge topologically equivalent atoms (y/n)? [yes] ",true);

	m_bWriteRDFs = AskYesNo("      Write all RDFs to a subdirectory (y/n)? [no] ",false);

	m_fDistRange = AskFloat("      Enter maximum plot distance for first maximum (in pm): [250] ",250.0);

	m_fHeightCut = AskFloat("      Enter maximum height cutoff for combined plot (in nm^-3, -1 = no cutoff): [5.0] ",5.0);
	if (m_fHeightCut <= 0)
		m_fHeightCut = 1.0e30;

	mprintf("\n");

	m_iaColumnIndex.resize(g_iGesVirtAtomCount);
	for (z=0;z<g_iGesVirtAtomCount;z++)
		m_iaColumnIndex[z] = -1;
	m_iCols = 0;

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];

		if (m_oaColumns[z] == NULL)
			continue;

		for (z3=0;z3<m_oaColumns[z]->m_oaAtoms.GetSize();z3++) {

			ma0 = NULL;

			for (z4=0;z4<((CxIntArray*)m_oaColumns[z]->m_oaAtoms[z3])->GetSize();z4++) {

				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
				ma = NULL;
				for (z2=0;z2<sm->m_oaMolAtoms.GetSize();z2++) {
					if ((((CMolAtom*)sm->m_oaMolAtoms[z2])->m_iType == m_oaColumns[z]->m_baAtomType[z3]) &&
						(((CMolAtom*)sm->m_oaMolAtoms[z2])->m_iNumber == ((CxIntArray*)m_oaColumns[z]->m_oaAtoms[z3])->GetAt(z4))) {
						ma = (CMolAtom*)sm->m_oaMolAtoms[z2];
						break;
					}
				}

				if (z4 > 0) {
					if ((ma0 != NULL) && (ma != NULL) && m_bJoinEquivalent)
						if (ma->m_liAtomCode == ma0->m_liAtomCode)
							goto _merge;
					m_iCols++;
_merge:;
				}

				if ((int)m_saColumnLabels.size() <= m_iCols) {
					m_saColumnLabels.push_back(CxString());
					m_iaColumnMolecule.push_back(z);
					m_iaColumnAtomCount.push_back(1);
					m_iaColumnRealElement.push_back(m_oaColumns[z]->m_baRealAtomType[z3]);
				}

				if (m_saColumnLabels[m_iCols].GetLength() != 0) {
					ts.Format(",%d",((CxIntArray*)m_oaColumns[z]->m_oaAtoms[z3])->GetAt(z4)+1);
					m_saColumnLabels[m_iCols] += ts;
					m_iaColumnAtomCount.back()++;
				} else
					m_saColumnLabels[m_iCols].Format("%s%d",(const char*)((CAtom*)g_oaAtoms[m_oaColumns[z]->m_baRealAtomType[z3]])->m_sName,((CxIntArray*)m_oaColumns[z]->m_oaAtoms[z3])->GetAt(z4)+1);
					//m_saColumnLabels[m_iCols].Format("%d.%s%d",z+1,(const char*)((CAtom*)g_oaAtoms[m_oaColumns[z]->m_baRealAtomType[z3]])->m_sName,((CxIntArray*)m_oaColumns[z]->m_oaAtoms[z3])->GetAt(z4)+1);

				ma0 = ma;

				for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++) {
					sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
					m_iaColumnIndex[((CxIntArray*)sm->m_oaAtomOffset[m_oaColumns[z]->m_baAtomType[z3]])->GetAt(((CxIntArray*)m_oaColumns[z]->m_oaAtoms[z3])->GetAt(z4))] = m_iCols;
				}
			}
			if (z4 != 0)
				m_iCols++;
		}
	}

	m_iaRowIndex.resize(g_iGesVirtAtomCount);
	for (z=0;z<g_iGesVirtAtomCount;z++)
		m_iaRowIndex[z] = -1;
	m_iRows = 0;

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];

		if (m_oaRows[z] == NULL)
			continue;

		for (z3=0;z3<m_oaRows[z]->m_oaAtoms.GetSize();z3++) {

			ma0 = NULL;

			for (z4=0;z4<((CxIntArray*)m_oaRows[z]->m_oaAtoms[z3])->GetSize();z4++) {

				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
				ma = NULL;
				for (z2=0;z2<sm->m_oaMolAtoms.GetSize();z2++) {
					if ((((CMolAtom*)sm->m_oaMolAtoms[z2])->m_iType == m_oaRows[z]->m_baAtomType[z3]) &&
						(((CMolAtom*)sm->m_oaMolAtoms[z2])->m_iNumber == ((CxIntArray*)m_oaRows[z]->m_oaAtoms[z3])->GetAt(z4))) {
						ma = (CMolAtom*)sm->m_oaMolAtoms[z2];
						break;
					}
				}

				if (z4 > 0) {
					if ((ma0 != NULL) && (ma != NULL) && m_bJoinEquivalent)
						if (ma->m_liAtomCode == ma0->m_liAtomCode)
							goto _merge2;
					m_iRows++;
_merge2:;
				}

				if ((int)m_saRowLabels.size() <= m_iRows) {
					m_saRowLabels.push_back(CxString());
					m_iaRowMolecule.push_back(z);
					m_iaRowAtomCount.push_back(1);
					m_iaRowRealElement.push_back(m_oaRows[z]->m_baRealAtomType[z3]);
				}

				if (m_saRowLabels[m_iRows].GetLength() != 0) {
					ts.Format(",%d",((CxIntArray*)m_oaRows[z]->m_oaAtoms[z3])->GetAt(z4)+1);
					m_saRowLabels[m_iRows] += ts;
					m_iaRowAtomCount.back()++;
				} else
					m_saRowLabels[m_iRows].Format("%s%d",(const char*)((CAtom*)g_oaAtoms[m_oaRows[z]->m_baRealAtomType[z3]])->m_sName,((CxIntArray*)m_oaRows[z]->m_oaAtoms[z3])->GetAt(z4)+1);
					//m_saRowLabels[m_iRows].Format("%d.%s%d",z+1,(const char*)((CAtom*)g_oaAtoms[m_oaRows[z]->m_baRealAtomType[z3]])->m_sName,((CxIntArray*)m_oaRows[z]->m_oaAtoms[z3])->GetAt(z4)+1);

				ma0 = ma;

				for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++) {
					sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
					m_iaRowIndex[((CxIntArray*)sm->m_oaAtomOffset[m_oaRows[z]->m_baAtomType[z3]])->GetAt(((CxIntArray*)m_oaRows[z]->m_oaAtoms[z3])->GetAt(z4))] = m_iRows;
				}
			}
			if (z4 != 0)
				m_iRows++;
		}
	}

	mprintf("      Contact matrix has %d rows and %d columns.\n\n",m_iRows,m_iCols);

	for (z=0;z<m_iRows;z++)
		mprintf("        Row %3d (molecule %d, %2d atoms): \"%s\"\n",z+1,m_iaRowMolecule[z]+1,m_iaRowAtomCount[z],(const char*)m_saRowLabels[z]);
	mprintf("\n");
	for (z=0;z<m_iCols;z++)
		mprintf("        Col %3d (molecule %d, %2d atoms): \"%s\"\n",z+1,m_iaColumnMolecule[z]+1,m_iaColumnAtomCount[z],(const char*)m_saColumnLabels[z]);
	mprintf("\n");

	if (m_bIntra) {
		mprintf("      Building bond neighbor matrices...\n");
		BuildNeighborMatrices();
	}

	mprintf(WHITE,"      <<< Contact Matrix Observation %d <<<\n\n",m_iObservation+1);

	return true;
}


void CContactMatrixObservation::BuildNeighborMatrices() {

	int z, z2, z3, z4, z5, z6;
	int ti, ti1, ti2, ti3, ti4;
	int o1, o2, o3, o4;
	CMolecule *m;
	CSingleMolecule *sm;
	CMolAtom *ma;


	m_iaLocalAtomIndex.resize(g_iGesAtomCount);

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];
		ti = 0;

		for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++) {

			if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
				continue;

			for (z4=0;z4<m->m_waAtomCount[z3];z4++) {

				for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++) {

					sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
					m_iaLocalAtomIndex[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4)] = ti;
				}

				ti++;
			}
		}

		if (ti != m->m_iAtomGesNoVirt)
			eprintf("CContactMatrixObservation::BuildNeighborMatrices(): Error: Count mismatch for molecule type %d: %d vs %d\n",z+1,ti,m->m_iAtomGesNoVirt);
	}

	m_iaaNeighborMatrices.resize(g_oaMolecules.GetSize());

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];
		sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];

		m_iaaNeighborMatrices[z].resize(m->m_iAtomGesNoVirt*m->m_iAtomGesNoVirt);
		for (z2=0;z2<m->m_iAtomGesNoVirt;z2++)
			for (z3=0;z3<m->m_iAtomGesNoVirt;z3++)
				if (z2 == z3)
					m_iaaNeighborMatrices[z][z2*m->m_iAtomGesNoVirt+z3] = 0;
				else
					m_iaaNeighborMatrices[z][z2*m->m_iAtomGesNoVirt+z3] = 9;

		for (z2=0;z2<sm->m_oaMolAtoms.GetSize();z2++) {

			ma = (CMolAtom*)sm->m_oaMolAtoms[z2];
			o1 = ma->m_iOffset;
			ti1 = m_iaLocalAtomIndex[o1];

			for (z3=0;z3<sm->m_oaBonds.GetSize();z3++) {

				if (((CMolBond*)sm->m_oaBonds[z3])->m_iAtomOffset[0] == o1)
					o2 = ((CMolBond*)sm->m_oaBonds[z3])->m_iAtomOffset[1];
				else if (((CMolBond*)sm->m_oaBonds[z3])->m_iAtomOffset[1] == o1)
					o2 = ((CMolBond*)sm->m_oaBonds[z3])->m_iAtomOffset[0];
				else
					continue;

				ti2 = m_iaLocalAtomIndex[o2];

				if (m_iaaNeighborMatrices[z][ti1*m->m_iAtomGesNoVirt+ti2] > 1)
					m_iaaNeighborMatrices[z][ti1*m->m_iAtomGesNoVirt+ti2] = 1;

				for (z4=0;z4<sm->m_oaBonds.GetSize();z4++) {

					if (((CMolBond*)sm->m_oaBonds[z4])->m_iAtomOffset[0] == o2)
						o3 = ((CMolBond*)sm->m_oaBonds[z4])->m_iAtomOffset[1];
					else if (((CMolBond*)sm->m_oaBonds[z4])->m_iAtomOffset[1] == o2)
						o3 = ((CMolBond*)sm->m_oaBonds[z4])->m_iAtomOffset[0];
					else
						continue;

					if (o3 == o1)
						continue;

					ti3 = m_iaLocalAtomIndex[o3];

					if (m_iaaNeighborMatrices[z][ti1*m->m_iAtomGesNoVirt+ti3] > 2)
						m_iaaNeighborMatrices[z][ti1*m->m_iAtomGesNoVirt+ti3] = 2;

					for (z5=0;z5<sm->m_oaBonds.GetSize();z5++) {

						if (((CMolBond*)sm->m_oaBonds[z5])->m_iAtomOffset[0] == o3)
							o4 = ((CMolBond*)sm->m_oaBonds[z5])->m_iAtomOffset[1];
						else if (((CMolBond*)sm->m_oaBonds[z5])->m_iAtomOffset[1] == o3)
							o4 = ((CMolBond*)sm->m_oaBonds[z5])->m_iAtomOffset[0];
						else
							continue;

						if ((o4 == o1) || (o4 == o2))
							continue;

						ti4 = m_iaLocalAtomIndex[o4];

						if (m_iaaNeighborMatrices[z][ti1*m->m_iAtomGesNoVirt+ti4] > 3)
							m_iaaNeighborMatrices[z][ti1*m->m_iAtomGesNoVirt+ti4] = 3;
					}
				}
			}
		}
	}

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];
		ti = 0;

		mprintf(WHITE,"\n    *** Molecule %d ***\n\n",z+1);

		mprintf("     ");

		for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++) {

			if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
				continue;

			for (z4=0;z4<m->m_waAtomCount[z3];z4++) 
				mprintf(" %2s%-2d",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,z4+1);
		}

		mprintf("\n");

		for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++) {

			if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
				continue;

			for (z4=0;z4<m->m_waAtomCount[z3];z4++) {

				mprintf("%2s%-2d ",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,z4+1);

				ti2 = 0;
				for (z5=0;z5<m->m_baAtomIndex.GetSize();z5++) {

					if (m->m_baAtomIndex[z5] == g_iVirtAtomType)
						continue;

					for (z6=0;z6<m->m_waAtomCount[z5];z6++) {
						if (m_iaaNeighborMatrices[z][ti*m->m_iAtomGesNoVirt+ti2] <= 3)
							mprintf("  %d  ",m_iaaNeighborMatrices[z][ti*m->m_iAtomGesNoVirt+ti2]);
						else
							mprintf("  -  ");
						ti2++;
					}
				}
				mprintf("\n");
				ti++;
			}
		}
	}
	mprintf("\n");
}


void CContactMatrixObservation::Initialize() {

	int z;

	m_iaFirstMaximum.resize(m_iRows*m_iCols);
	m_iaFirstMinimum.resize(m_iRows*m_iCols);
	m_oaMatrix = new CDF[m_iRows*m_iCols];

	for (z=0;z<m_iRows*m_iCols;z++) {

		m_oaMatrix[z].m_iResolution = m_iResolution;
		m_oaMatrix[z].m_fMinVal = 0;
		m_oaMatrix[z].m_fMaxVal = m_fMaxDist;
		m_oaMatrix[z].Create();
		m_oaMatrix[z].SetLabelX("Distance / pm");
		m_oaMatrix[z].SetLabelY("g(r) / nm^-3");

		m_iaFirstMaximum[z] = 0;
		m_iaFirstMinimum[z] = 0;
	}
}


void CContactMatrixObservation::Finish() {

	CMatrixPlot mp;
	CxString buf;
	CDF *df;
	int zr, zc, z;
	double tf;
	char dirsep;
	FILE *a;


	mprintf(WHITE,"    *** Observation %d ***\n",m_iObservation+1);

	mprintf("        Normalizing RDFs...\n");

	for (zr=0;zr<m_iRows;zr++) {

		for (zc=0;zc<m_iCols;zc++) {

			df = &m_oaMatrix[zr*m_iCols+zc];

			df->MultiplyBin(1.0 / g_iSteps * g_iStride);

			df->CorrectRadialDist();

			df->Integrate(true,
				1.0 /
				(double)((CMolecule*)g_oaMolecules[m_iaRowMolecule[zr]])->m_laSingleMolIndex.GetSize() / 
				m_iaRowAtomCount[zr]
			);

			if (((CMolecule*)g_oaMolecules[m_iaRowMolecule[zr]])->m_laSingleMolIndex.GetSize()*m_iaRowAtomCount[zr] >
				((CMolecule*)g_oaMolecules[m_iaColumnMolecule[zc]])->m_laSingleMolIndex.GetSize()*m_iaColumnAtomCount[zc])
				tf = ((CMolecule*)g_oaMolecules[m_iaColumnMolecule[zc]])->m_laSingleMolIndex.GetSize()*m_iaColumnAtomCount[zc];
			else
				tf = ((CMolecule*)g_oaMolecules[m_iaRowMolecule[zr]])->m_laSingleMolIndex.GetSize()*m_iaRowAtomCount[zr];

		/*	df->MultiplyBin(
				tvol /
				((CMolecule*)g_oaMolecules[m_iaRowMolecule[zr]])->m_laSingleMolIndex.GetSize() /
				m_iaColumnAtomCount[zc] /
				((CMolecule*)g_oaMolecules[m_iaColumnMolecule[zc]])->m_laSingleMolIndex.GetSize() /
				m_iaRowAtomCount[zr]
			);*/

			// Values in nm^-3
			df->MultiplyBin( 1e9 / (4.0/3.0*Pi) / tf );

			if (g_bDoubleBox) {
				df->MultiplyBin(g_iDoubleBoxFactor);
				df->MultiplyIntegral(g_iDoubleBoxFactor);
			}

			tf = 0;
			for (z=0;z<m_iResolution;z++) {
				if (df->m_pBin[z] > tf) {
					tf = df->m_pBin[z];
					m_iaFirstMaximum[zr*m_iCols+zc] = z;
				}
				if (tf == 0)
					continue;
				// 0.2 nm^-3 is the threshold for the first maximum
				if ((tf > 0.2) && ((df->m_pBin[z] / tf) < 0.8)) {
					m_iaFirstMinimum[zr*m_iCols+zc] = z;
					break;
				}
			}
			if (m_iaFirstMinimum[zr*m_iCols+zc] == 0)
				m_iaFirstMaximum[zr*m_iCols+zc] = 0;
		}
	}

	if (m_bWriteRDFs) {

		#ifdef TARGET_WINDOWS
			dirsep = '\\';
		#else
			dirsep = '/';
		#endif

		mprintf("        Creating subdirectory \"cmat_rdf_%d\"..\n",m_iObservation+1);

		buf.sprintf("mkdir cmat_rdf_%d",m_iObservation+1);
		(void)!system((const char*)buf);

		mprintf("        Writing %d RDFs to subdirectory \"cmat_rdf_%d\"..\n",m_iRows*m_iCols,m_iObservation+1);

		for (zr=0;zr<m_iRows;zr++) {

			for (zc=0;zc<m_iCols;zc++) {

				df = &m_oaMatrix[zr*m_iCols+zc];

				buf.sprintf("cmat_rdf_%d%crdf_r%03d_c%03d_%s_%s_%s_%s.csv",
					m_iObservation+1,
					dirsep,
					zr+1,
					zc+1,
					(const char*)((CMolecule*)g_oaMolecules[m_iaRowMolecule[zr]])->m_sName,
					(const char*)m_saRowLabels[zr],
					(const char*)((CMolecule*)g_oaMolecules[m_iaColumnMolecule[zc]])->m_sName,
					(const char*)m_saColumnLabels[zc]
				);

				df->Write("",buf,"",true);
			}
		}
	}

	mprintf("        Preparing Matrix Plot...\n");

	mp.Init(m_iRows,m_iCols);

	for (z=0;z<m_iRows;z++) {
		mp.SetRowLabel(z,(const char*)m_saRowLabels[z]);
		if (z+1 < m_iRows) {
			if (m_iaRowMolecule[z+1] != m_iaRowMolecule[z])
				mp.AddRowCategory(z);
			else if (m_iaRowRealElement[z+1] != m_iaRowRealElement[z])
				mp.AddRowSubCategory(z);
		}
	}

	for (z=0;z<m_iCols;z++) {
		mp.SetColumnLabel(z,(const char*)m_saColumnLabels[z]);
		if (z+1 < m_iCols) {
			if (m_iaColumnMolecule[z+1] != m_iaColumnMolecule[z])
				mp.AddColumnCategory(z);
			else if (m_iaColumnRealElement[z+1] != m_iaColumnRealElement[z])
				mp.AddColumnSubCategory(z);
		}
	}

	buf.Format("contactmatrix_%d_maxheight.csv",m_iObservation+1);
	mprintf("        Writing \"%s\"...\n",(const char*)buf);
	a = OpenFileWrite((const char*)buf,true);
	for (z=0;z<m_iCols;z++)
		mfprintf(a,"; %s",(const char*)m_saColumnLabels[z]);
	mfprintf(a,"\n");
	for (zr=0;zr<m_iRows;zr++) {
		mfprintf(a,"%s",(const char*)m_saRowLabels[zr]);
		for (zc=0;zc<m_iCols;zc++) {
			tf = m_oaMatrix[zr*m_iCols+zc].m_pBin[m_iaFirstMaximum[zr*m_iCols+zc]];
			if (m_iaFirstMaximum[zr*m_iCols+zc] != 0)
				mp.m_faBin[zr*m_iCols+zc] = tf;
			else
				mp.m_iaActive[zr*m_iCols+zc] = 0;
			mfprintf(a,"; %.3f",tf);
		}
		mfprintf(a,"\n");
	}
	fclose(a);

	buf.Format("contactmatrix_%d_maxheight.svg",m_iObservation+1);
	mprintf("        Writing \"%s\"...\n",(const char*)buf);
	mp.SetValueLabel("First Maximum Height / nm^-3");
	mp.SetRangeZero(true);
	mp.SetPlotExponent(0.5);
	mp.WriteSVGPlot((const char*)buf);

	buf.Format("contactmatrix_%d_maxdist.csv",m_iObservation+1);
	mprintf("        Writing \"%s\"...\n",(const char*)buf);
	a = OpenFileWrite((const char*)buf,true);
	for (z=0;z<m_iCols;z++)
		mfprintf(a,"; %s",(const char*)m_saColumnLabels[z]);
	mfprintf(a,"\n");
	for (zr=0;zr<m_iRows;zr++) {
		mfprintf(a,"%s",(const char*)m_saRowLabels[zr]);
		for (zc=0;zc<m_iCols;zc++) {
			tf = (m_iaFirstMaximum[zr*m_iCols+zc]+0.5) * m_fMaxDist / m_iResolution;
			if (m_iaFirstMaximum[zr*m_iCols+zc] != 0)
				mp.m_faBin[zr*m_iCols+zc] = tf;
			else
				mp.m_iaActive[zr*m_iCols+zc] = 0;
			mfprintf(a,"; %.3f",tf);
		}
		mfprintf(a,"\n");
	}
	fclose(a);

	buf.Format("contactmatrix_%d_maxdist.svg",m_iObservation+1);
	mprintf("        Writing \"%s\"...\n",(const char*)buf);
	mp.SetValueLabel("First Maximum Distance / pm");
	mp.SetRangeZero(false);
	mp.SetPlotExponent(1.0);
	mp.SetInvertColorScale(true);
	mp.WriteSVGPlot((const char*)buf);

	for (zr=0;zr<m_iRows;zr++) {
		for (zc=0;zc<m_iCols;zc++) {
			if (m_iaFirstMaximum[zr*m_iCols+zc] != 0) {
				mp.m_iaActive[zr*m_iCols+zc] = 1;
				mp.m_faBin[zr*m_iCols+zc]  = (m_iaFirstMaximum[zr*m_iCols+zc]+0.5) * m_fMaxDist / m_iResolution;
				mp.m_faBin2[zr*m_iCols+zc] = m_oaMatrix[zr*m_iCols+zc].m_pBin[m_iaFirstMaximum[zr*m_iCols+zc]];
			} else
				mp.m_iaActive[zr*m_iCols+zc] = 0;
		}
	}
	buf.Format("contactmatrix_%d_combined.svg",m_iObservation+1);
	mprintf("        Writing \"%s\"...\n",(const char*)buf);
	mp.Set2DMode(true);
	mp.SetValueLabel("First Maximum Distance / pm");
	mp.SetValueLabel2("First Maximum Height / nm^-3");
	mp.SetRangeZero(false);
	mp.SetRangeZero2(true);
	mp.SetRangeCut(m_fDistRange);
	mp.SetRangeCut2(m_fHeightCut);
	mp.SetPlotExponent(1.0);
	mp.SetPlotExponent2(1.0);
	mp.SetInvertColorScale(true);
	mp.SetInvertColorScale2(false);
	mp.WriteSVGPlot((const char*)buf);
}


CContactMatrix::CContactMatrix() {

	m_pDomainEngine = NULL;
}


CContactMatrix::~CContactMatrix() {
}


bool CContactMatrix::Parse() {

	CContactMatrixObservation *o;


	mprintf(WHITE,"\n    >>> Contact Matrix Analysis >>>\n");

	m_fMaxDist = 0;

	while (true) {

		o = new CContactMatrixObservation();
		o->m_iObservation = (int)m_oaObservations.size();
		m_oaObservations.push_back(o);

		if (!o->Parse())
			return false;

		if (o->m_fMaxDist > m_fMaxDist)
			m_fMaxDist = o->m_fMaxDist;

		if (!AskYesNo("    Add another contact matrix observation (y/n)? [no] ",false))
			break;
	}

	mprintf("\n    The maximum contact distance is %.3f pm.\n",m_fMaxDist);

	mprintf(WHITE,"\n    <<< Contact Matrix Analysis <<<\n\n");

	return true;
}


void CContactMatrix::Initialize() {

	int z;

	mprintf("  Initializing Contact Matrix Observations...\n");

	m_pDomainEngine = new CPosDomainEngine();

	for (z=0;z<(int)m_oaObservations.size();z++) {
		mprintf("    Observation %d...\n",z+1);
		m_oaObservations[z]->Initialize();
	}
}


void CContactMatrix::ProcessStep(CTimeStep *ts) {

	int z, z2, z3, z4, o1, o2, ti;
	CContactMatrixObservation *o;
	std::vector<int> dom, nbh;
	double tf;


	m_pDomainEngine->Create(ts->m_vaCoords,g_iGesAtomCount,g_mBoxFromOrtho,m_fMaxDist);

	if (m_pDomainEngine->m_bTrivial) {

	} else {

		for (z=0;z<m_pDomainEngine->m_iRes[0]*m_pDomainEngine->m_iRes[1]*m_pDomainEngine->m_iRes[2];z++) {

			m_pDomainEngine->GetDomainAndNeighbors(z,dom,nbh);

			for (z2=0;z2<(int)dom.size();z2++) {

				o1 = dom[z2];

				for (z3=0;z3<(int)nbh.size();z3++) {

					o2 = nbh[z3];

					if (o1 >= o2)
						continue;

					for (z4=0;z4<(int)m_oaObservations.size();z4++) {

						o = m_oaObservations[z4];

						if (g_laAtomSMIndex[o1] == g_laAtomSMIndex[o2]) {
							if (!o->m_bIntra)
								continue;
							ti = o->m_iaaNeighborMatrices[g_waAtomMolIndex[o1]][o->m_iaLocalAtomIndex[o1]*((CMolecule*)g_oaMolecules[g_waAtomMolIndex[o1]])->m_iAtomGesNoVirt+o->m_iaLocalAtomIndex[o2]];
							if (o->m_bExclude12 && (ti == 1))
								continue;
							if (o->m_bExclude13 && (ti == 2))
								continue;
							if (o->m_bExclude14 && (ti == 3))
								continue;
						} else {
							if (!o->m_bInter)
								continue;
							if (o->m_bIntra && (g_waAtomMolIndex[o1] == g_waAtomMolIndex[o2]))
								continue;
						}

						if ((o->m_iaRowIndex[o1] != -1) && (o->m_iaColumnIndex[o2] != -1)) {

							tf = FoldedLength(ts->m_vaCoords[o2] - ts->m_vaCoords[o1]);
							o->m_oaMatrix[o->m_iaRowIndex[o1]*o->m_iCols+o->m_iaColumnIndex[o2]].AddToBin(tf);

							if ((o->m_iaRowIndex[o2] != -1) && (o->m_iaColumnIndex[o1] != -1))
								o->m_oaMatrix[o->m_iaRowIndex[o2]*o->m_iCols+o->m_iaColumnIndex[o1]].AddToBin(tf);

						} else if ((o->m_iaRowIndex[o2] != -1) && (o->m_iaColumnIndex[o1] != -1)) {

							tf = FoldedLength(ts->m_vaCoords[o2] - ts->m_vaCoords[o1]);
							o->m_oaMatrix[o->m_iaRowIndex[o2]*o->m_iCols+o->m_iaColumnIndex[o1]].AddToBin(tf);
						}
					}
				}
			}
		}
	}
}


void CContactMatrix::Finish() {

	int z;

	mprintf(WHITE,"\n    Finishing Contact Matrix Observations...\n");
	for (z=0;z<(int)m_oaObservations.size();z++)
		m_oaObservations[z]->Finish();
}





