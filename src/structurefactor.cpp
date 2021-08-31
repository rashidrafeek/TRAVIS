/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Martin Thomas.

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

#include "structurefactor.h"

#include "globalvar.h"
#include "maintools.h"
#include "timestep.h"
#include "xobarray.h"
//>>>OMP
#ifdef USE_OMP
#include <omp.h>
#endif
//<<<OMP

const char *GetRevisionInfo_structurefactor(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_structurefactor() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


#define BUF_SIZE 4096


static CxObArray g_isotopes;
static CStructureFactor *g_structureFactor;
static CxObArray g_sFacObserv;


static void createIsotopeList() {
	// Neutron factors from http://www.ncnr.nist.gov/resources/n-lengths/
	// X-Ray factors from http://www.ruppweb.org/new_comp/scattering_factors.htm
	CIsotope *isotope;
	try { isotope = new CIsotope("H", -3.7390, 0.493, 0.323, 0.140, 0.041, 10.511, 26.126, 3.142, 57.800, 0.003 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("1H", -3.7406, 0.493, 0.323, 0.140, 0.041, 10.511, 26.126, 3.142, 57.800, 0.003 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("D", 6.671, 0.493, 0.323, 0.140, 0.041, 10.511, 26.126, 3.142, 57.800, 0.003 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("2H", 6.671, 0.493, 0.323, 0.140, 0.041, 10.511, 26.126, 3.142, 57.800, 0.003 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("3H", 4.792, 0.493, 0.323, 0.140, 0.041, 10.511, 26.126, 3.142, 57.800, 0.003 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("Li", -1.90, 1.128, 0.751, 0.618, 0.465, 3.955, 1.052, 85.391, 168.261, 0.038 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("6Li", 2.00, 1.128, 0.751, 0.618, 0.465, 3.955, 1.052, 85.391, 168.261, 0.038 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("7Li", -2.22, 1.128, 0.751, 0.618, 0.465, 3.955, 1.052, 85.391, 168.261, 0.038 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("B", 5.30, 2.055, 1.333, 1.098, 0.707, 23.219, 1.021, 60.350, 0.140, -0.193 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("10B", -0.1, 2.055, 1.333, 1.098, 0.707, 23.219, 1.021, 60.350, 0.140, -0.193 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("11B", 6.65, 2.055, 1.333, 1.098, 0.707, 23.219, 1.021, 60.350, 0.140, -0.193 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("C", 6.6460, 2.310, 1.020, 1.589, 0.865, 20.844, 10.208, 0.569, 51.651, 0.216 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("12C", 6.6511, 2.310, 1.020, 1.589, 0.865, 20.844, 10.208, 0.569, 51.651, 0.216 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("13C", 6.19, 2.310, 1.020, 1.589, 0.865, 20.844, 10.208, 0.569, 51.651, 0.216 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("N", 9.36, 12.213, 3.132, 2.013, 1.166, 0.006, 9.893, 28.997, 0.583, -11.529 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("14N", 9.37, 12.213, 3.132, 2.013, 1.166, 0.006, 9.893, 28.997, 0.583, -11.529 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("15N", 6.44, 12.213, 3.132, 2.013, 1.166, 0.006, 9.893, 28.997, 0.583, -11.529 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("O", 5.803, 3.049, 2.287, 1.546, 0.867, 13.277, 5.701, 0.324, 32.909, 0.251 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("16O", 5.803, 3.049, 2.287, 1.546, 0.867, 13.277, 5.701, 0.324, 32.909, 0.251 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("17O", 5.78, 3.049, 2.287, 1.546, 0.867, 13.277, 5.701, 0.324, 32.909, 0.251 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("18O", 5.84, 3.049, 2.287, 1.546, 0.867, 13.277, 5.701, 0.324, 32.909, 0.251 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("F", 5.654, 3.539, 2.641, 1.517, 1.024, 10.283, 4.294, 0.262, 26.148, 0.278 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("P", 5.13, 6.435, 4.179, 1.780, 1.491, 1.907, 27.157, 0.526, 68.164, 1.115 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("S", 2.847, 6.905, 5.203, 1.438, 1.586, 1.468, 22.215, 0.254, 56.172, 0.867 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("32S", 2.804, 6.905, 5.203, 1.438, 1.586, 1.468, 22.215, 0.254, 56.172, 0.867 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("33S", 4.74, 6.905, 5.203, 1.438, 1.586, 1.468, 22.215, 0.254, 56.172, 0.867 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("34S", 3.48, 6.905, 5.203, 1.438, 1.586, 1.468, 22.215, 0.254, 56.172, 0.867 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("36S", 3.1, 6.905, 5.203, 1.438, 1.586, 1.468, 22.215, 0.254, 56.172, 0.867 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("Cl", 9.577, 11.460, 7.196, 6.256, 1.645, 0.010, 1.166, 18.519, 47.778, -9.557 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("35Cl", 11.65, 11.460, 7.196, 6.256, 1.645, 0.010, 1.166, 18.519, 47.778, -9.557 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("37Cl", 3.08, 11.460, 7.196, 6.256, 1.645, 0.010, 1.166, 18.519, 47.778, -9.557 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("Br", 6.795, 17.179, 5.236, 5.638, 3.985, 2.172, 16.580, 0.261, 41.433, -2.956 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("79Br", 6.80, 17.179, 5.236, 5.638, 3.985, 2.172, 16.580, 0.261, 41.433, -2.956 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
	try { isotope = new CIsotope("81Br", 6.79, 17.179, 5.236, 5.638, 3.985, 2.172, 16.580, 0.261, 41.433, -2.956 ); } catch(...) { isotope = NULL; }
	if(isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_isotopes.Add(isotope);
}


static void deleteIsotopeList() {
	int i;
	for(i = 0; i < g_isotopes.GetSize(); i++)
		delete (CIsotope *)g_isotopes[i];
}


struct StructureFactorAtomKind: public CxObject
{
	int atomType;
	CxObArray isotopeList;
};


struct StructureFactorMolecule: public CxObject
{
	int moleculeType;
	CxObArray atomKinds;
};


CIsotope::CIsotope(const char *label, double neutronFactor, double cma1, double cma2, double cma3, double cma4, double cmb1, double cmb2, double cmb3, double cmb4, double cmc) {
	try { _label = new char[strlen(label) + 1]; } catch(...) { _label = NULL; }
	if(_label == NULL) NewException((double)sizeof(char) * (strlen(label) + 1), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	strcpy(_label, label);
	_neutronFactor = neutronFactor;
	_cma[0] = cma1;
	_cma[1] = cma2;
	_cma[2] = cma3;
	_cma[3] = cma4;
	_cmb[0] = cmb1;
	_cmb[1] = cmb2;
	_cmb[2] = cmb3;
	_cmb[3] = cmb4;
	_cmc = cmc;
}


CIsotope::~CIsotope() {
	delete[] _label;
}


double CIsotope::xrayFactor(double q) {
	double x = q * 100.0 / 4.0 / Pi;
	double factor = 0.0;
	int i;
	for(i = 0; i < 4; i++) {
		factor += _cma[i] * exp(-_cmb[i] * x * x);
	}
	factor += _cmc;
	return factor;
}


CStructureFactorGroup::CStructureFactorGroup(bool global, CxObArray &isotopeAssignList) {

	int i, j, k, l;


	if (global) {

		for (i = 0; i < g_oaMolecules.GetSize(); i++) {
			CAtomGroup *ag;
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			if(!ag->ParseAtoms((CMolecule *)g_oaMolecules[i], "*")) {
				eprintf("CStructureFactorGroup::CStructureFactorGroup(): Internal error.\n");
				abort();
			}
			m_atomGroupList.Add(ag);
		}
		m_sepInterIntra = AskYesNo("    Separate intramolecular and intermolecular contributions in the total structure factor (y/n)? [no] ", false);
		mprintf("\n");

	} else {

		CxString buf, buf2;

		while (true) {
			buf.sprintf("      Take atoms from which molecule (");
			for (i = 0; i < g_oaMolecules.GetSize(); i++) {
				buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
				buf.strcat(buf2);
				if (i < g_oaMolecules.GetSize() - 1) {
					buf.strcat(", ");
				}
			}
			buf.strcat(")? ");
			
			int mol = AskRangeInteger_ND("%s", 1, g_oaMolecules.GetSize(),(const char*)buf) - 1;
			
			CAtomGroup *ag;
			try { ag = new CAtomGroup(); } catch(...) { ag = NULL; }
			if (ag == NULL) NewException((double)sizeof(CAtomGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			
			while (true) {
				if (((CMolecule *)g_oaMolecules[mol])->m_iAtomGesNoVirt == 1) {
					mprintf("      %s is only one atom, there is no choice.\n", ((CMolecule *)g_oaMolecules[mol])->m_sName);
					ag->Reset();
					ag->m_pMolecule = (CMolecule *)g_oaMolecules[mol];
					ag->AddAtom(0, 0, false);
					ag->SortAtoms();
					ag->BuildName();
				} else {
					AskString("      Which atom(s) to take from molecule %s (e.g. \"C1,C3-5,H\", \"*\"=all)? [*] ", &buf, "*", ((CMolecule*)g_oaMolecules[mol])->m_sName);
					if(!ag->ParseAtoms((CMolecule *)g_oaMolecules[mol], buf)) {
						continue;
					}
					bool containsVirt = false;
					for (i = 0; i < ag->m_baAtomType.GetSize(); i++) {
						if (ag->m_baAtomType[i] == g_iVirtAtomType) {
							containsVirt = true;
							break;
						}
					}
					if (containsVirt) {
						eprintf("The selection must not contain virtual atoms.\n");
						continue;
					}
				}
				break;
			}
			m_atomGroupList.Add(ag);
			
			if (!AskYesNo("\n      Add atoms from another molecule (y/n)? [no] ", false))
				break;
			mprintf("\n");
		}
		mprintf("\n");
		
		m_sepInterIntra = AskYesNo("      Separate intramolecular and intermolecular contributions (y/n)? [no] ", false);
		mprintf("\n");
	}
	
	if (global) {
		m_name.sprintf("total");
	} else {
		CxString temp;
		for (i = 0; i < m_atomGroupList.GetSize(); i++) {
			if (i > 0)
				temp.strcat("__");
			temp.strcat(((CAtomGroup *)m_atomGroupList[i])->m_pMolecule->m_sName);
			if (((CAtomGroup *)m_atomGroupList[i])->m_pMolecule->m_iAtomGes > 1) {
				temp.strcat("_");
				temp.strcat(((CAtomGroup *)m_atomGroupList[i])->m_sName);
			}
		}
		
		m_name.sprintf("[%s]", (const char *)temp);
	}
	
	for (i = 0; i < m_atomGroupList.GetSize(); i++) {
		CAtomGroup *ag = (CAtomGroup *)m_atomGroupList[i];
		CMolecule *mol = ag->m_pMolecule;
		for (j = 0; j < ag->m_baAtomType.GetSize(); j++) {
			CxIntArray *a = (CxIntArray *)ag->m_oaAtoms[j];
			for (k = 0; k < a->GetSize(); k++) {
				CIsotope *isotope = (CIsotope *)((StructureFactorAtomKind *)((StructureFactorMolecule *)isotopeAssignList.GetAt(mol->m_iIndex))->atomKinds[ag->m_baAtomType[j]])->isotopeList[k];
				int isotopeIndex = -1;
				for (l = 0; l < m_isotopeTypeList.GetSize(); l++) {
					if((CIsotope *)m_isotopeTypeList[l] == isotope) {
						isotopeIndex = l;
						break;
					}
				}
				if (isotopeIndex == -1) {
					isotopeIndex = m_isotopeTypeList.GetSize();
					m_isotopeTypeList.Add(isotope);
					m_isotopeTypeCount.Add(0);
				}
				for (l = 0; l < mol->m_laSingleMolIndex.GetSize(); l++) {
					CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[mol->m_laSingleMolIndex[l]];
					m_atomIndexList.Add(((CxIntArray *)sm->m_oaAtomOffset[ag->m_baAtomType[j]])->GetAt(a->GetAt(k)));
					m_singleMolList.Add(mol->m_laSingleMolIndex[l]);
					m_isotopeList.Add(isotopeIndex);
					m_isotopeTypeCount[isotopeIndex]++;
				}
			}
		}
	}
	
	m_isotopeTypeTotalCount.SetSize(m_isotopeTypeCount.GetSize());
	for (i = 0; i < m_isotopeTypeTotalCount.GetSize(); i++)
		m_isotopeTypeTotalCount[i] = 0;
	
	for (i = 0; i < isotopeAssignList.GetSize(); i++) {
		StructureFactorMolecule *smol = (StructureFactorMolecule *)isotopeAssignList[i];
		for (j = 0; j < smol->atomKinds.GetSize(); j++) {
			StructureFactorAtomKind *satom = (StructureFactorAtomKind *)smol->atomKinds[j];
			for (k = 0; k < satom->isotopeList.GetSize(); k++) {
				CIsotope *isotope = (CIsotope *)satom->isotopeList[k];
				int isotopeIndex = -1;
				for (l = 0; l < m_isotopeTypeList.GetSize(); l++) {
					if ((CIsotope *)m_isotopeTypeList[l] == isotope) {
						isotopeIndex = l;
						break;
					}
				}
				if (isotopeIndex != -1) {
					m_isotopeTypeTotalCount[isotopeIndex] += ((CMolecule *)g_oaMolecules[i])->m_laSingleMolIndex.GetSize();
				}
			}
		}
	}
	
	m_global = global;
}


CStructureFactorGroup::~CStructureFactorGroup() {
	int i;
	for (i = 0; i < m_atomGroupList.GetSize(); i++)
		delete (CAtomGroup *)m_atomGroupList[i];
	for (i = 0; i < m_rdfList.GetSize(); i++)
		delete (CDF *)m_rdfList[i];
	for (i = 0; i < m_rdfIntraList.GetSize(); i++)
		delete (CDF *)m_rdfIntraList[i];
	for (i = 0; i < m_rdfInterList.GetSize(); i++)
		delete (CDF *)m_rdfInterList[i];
}


void CStructureFactorGroup::initialize(double rdfMax, int rdfRes) {
	int i;
	mprintf("      Creating %d RDFs\n", m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2);
	for (i = 0; i < m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2; i++) {
		CDF *df;
		try { df = new CDF(); } catch(...) { df = NULL; }
		if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		df->m_fMinVal = 0.0;
		df->m_fMaxVal = rdfMax;
		df->m_iResolution = rdfRes;
		df->SetLabelX("Distance / pm");
		df->SetLabelY("g(r)");
		df->Create();
		m_rdfList.Add(df);
	}
	
	if (m_sepInterIntra) {
		mprintf("      Creating %d intramolecular RDFs\n", m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2);
		for (i = 0; i < m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2; i++) {
			CDF *df;
			try { df = new CDF(); } catch(...) { df = NULL; }
			if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			df->m_fMinVal = 0.0;
			df->m_fMaxVal = rdfMax;
			df->m_iResolution = rdfRes;
			df->SetLabelX("Distance / pm");
			df->SetLabelY("g(r)");
			df->Create();
			m_rdfIntraList.Add(df);
		}
		mprintf("      Creating %d intermolecular RDFs\n", m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2);
		for (i = 0; i < m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2; i++) {
			CDF *df;
			try { df = new CDF(); } catch(...) { df = NULL; }
			if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			df->m_fMinVal = 0.0;
			df->m_fMaxVal = rdfMax;
			df->m_iResolution = rdfRes;
			df->SetLabelX("Distance / pm");
			df->SetLabelY("g(r)");
			df->Create();
			m_rdfInterList.Add(df);
		}
	}
//>>>OMP
#ifdef USE_OMP
	m_tempLists.resize(omp_get_max_threads());
	for ( int j=0; j < omp_get_max_threads(); j++ )
	{
		m_tempLists[j].resize(m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2);
		for (i = 0; i < m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2; i++) 
		{
			CDF *df;
			try { df = new CDF(); } catch(...) { df = NULL; }
			if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			df->m_fMinVal = 0.0;
			df->m_fMaxVal = rdfMax;
			df->m_iResolution = rdfRes;
			df->SetLabelX("Distance / pm");
			df->SetLabelY("g(r)");
			df->Create();
			m_tempLists[j][i] = df;
		}
	}
	if (m_sepInterIntra) 
	{
		m_tempInterLists.resize(omp_get_max_threads());
		m_tempIntraLists.resize(omp_get_max_threads());
		for ( int j=0; j < omp_get_max_threads(); j++ )
		{
			m_tempInterLists[j].resize(m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2);
			m_tempIntraLists[j].resize(m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2);
			for (i = 0; i < m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2; i++) 
			{
				CDF *df;
				try { df = new CDF(); } catch(...) { df = NULL; }
				if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				df->m_fMinVal = 0.0;
				df->m_fMaxVal = rdfMax;
				df->m_iResolution = rdfRes;
				df->SetLabelX("Distance / pm");
				df->SetLabelY("g(r)");
				df->Create();
				m_tempInterLists[j][i] = df;
			}
			for (i = 0; i < m_isotopeTypeList.GetSize() * (m_isotopeTypeList.GetSize() + 1) / 2; i++) 
			{
				CDF *df;
				try { df = new CDF(); } catch(...) { df = NULL; }
				if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				df->m_fMinVal = 0.0;
				df->m_fMaxVal = rdfMax;
				df->m_iResolution = rdfRes;
				df->SetLabelX("Distance / pm");
				df->SetLabelY("g(r)");
				df->Create();
				m_tempIntraLists[j][i] = df;
			}
		}
	}
#endif
//<<<OMP
}


void CStructureFactorGroup::process(CTimeStep *ts) {
	int i, j, a, b, index;
	double dist;
//>>>OMP
#ifdef USE_OMP
#pragma omp parallel for private(i,j,a,b,index,dist)
#endif
//<<<OMP
	for (i = 0; i < m_atomIndexList.GetSize(); i++) {
		//if ((i%100) == 0)
		//	mprintf("\n@Group i=%d",i);
		for (j = i+1; j < m_atomIndexList.GetSize(); j++) {
			dist = FoldedLength(ts->m_vaCoords[m_atomIndexList[i]] - ts->m_vaCoords[m_atomIndexList[j]]);
			a = m_isotopeList[i];
			b = m_isotopeList[j];
			if (a < b)
				index = (m_isotopeTypeList.GetSize()-1)*a - a*(a-1)/2 + b;
			else
				index = (m_isotopeTypeList.GetSize()-1)*b - b*(b-1)/2 + a;
//>>>OMP
#ifdef USE_OMP
			m_tempLists[omp_get_thread_num()][index]->AddToBin(dist);
			if (m_sepInterIntra) 
			{
				if (m_singleMolList[i] == m_singleMolList[j])
					m_tempIntraLists[omp_get_thread_num()][index]->AddToBin(dist);
				else
					m_tempInterLists[omp_get_thread_num()][index]->AddToBin(dist);
			}
#else
//<<<OMP
			((CDF *)m_rdfList[index])->AddToBin(dist);
			if (m_sepInterIntra) {
				if (m_singleMolList[i] == m_singleMolList[j])
					((CDF *)m_rdfIntraList[index])->AddToBin(dist);
				else
					((CDF *)m_rdfInterList[index])->AddToBin(dist);
			}
//>>>OMP
#endif
//<<<OMP
		}
	}
}


void CStructureFactorGroup::finalize(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch) {
	int i, j, k, l;
	CxString filename;

//>>>OMP
#ifdef USE_OMP
	for ( i=0; i < (int)m_tempLists.size(); i++ )
		for ( j=0; j < (int)m_tempLists[i].size(); j++ )
			((CDF *)m_rdfList[j])->AddFrom(m_tempLists[i][j]);
	for ( i=0; i < (int)m_tempInterLists.size(); i++ )
		for ( j=0; j < (int)m_tempInterLists[i].size(); j++ )
			((CDF *)m_rdfInterList[j])->AddFrom(m_tempInterLists[i][j]);
	for ( i=0; i < (int)m_tempIntraLists.size(); i++ )
		for ( j=0; j < (int)m_tempIntraLists[i].size(); j++ )
			((CDF *)m_rdfIntraList[j])->AddFrom(m_tempIntraLists[i][j]);
#endif
//<<<OMP
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_isotopeTypeList.GetSize(); i++) {
		for (j = i; j < m_isotopeTypeList.GetSize(); j++) {
			mprintf("      Processing partial structure factor %s-%s\n", ((CIsotope *)m_isotopeTypeList[i])->label(), ((CIsotope *)m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfList[(m_isotopeTypeList.GetSize()-1)*i - i*(i-1)/2 + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if(i == j)
				df->MultiplyBin(2.0);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0 / (4.0/3.0 * Pi) / ((double)m_isotopeTypeTotalCount[i]) / ((double)m_isotopeTypeTotalCount[j]));
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0/3.0 * Pi) / ((double)m_isotopeTypeTotalCount[i]) / ((double)m_isotopeTypeTotalCount[j]));
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
// 				df->Integrate(true, (4.0f/3.0f * Pi) * m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j] / (g_fBoxX * g_fBoxY * g_fBoxZ) / m_isotopeTypeCount[i]);
				filename.sprintf("sfac_%s%s_rdf_%s_%s.csv", m_global ? "" : "self_", (const char *)m_name, (const char *)(((CIsotope *)m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_isotopeTypeList[j])->label()));
				mprintf("      Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			double shift = ((double)m_isotopeTypeCount[i]) * ((double)m_isotopeTypeCount[j]) / (((double)m_isotopeTypeTotalCount[i]) * ((double)m_isotopeTypeTotalCount[j]));

	/*		mprintf("@ m_isotopeTypeCount[i] = %d\n",m_isotopeTypeCount[i]);
			mprintf("@ m_isotopeTypeCount[j] = %d\n",m_isotopeTypeCount[j]);
			mprintf("@ m_isotopeTypeTotalCount[i] = %d\n",m_isotopeTypeTotalCount[i]);
			mprintf("@ m_isotopeTypeTotalCount[j] = %d\n",m_isotopeTypeTotalCount[j]);
			mprintf("@ shift = %.10G\n",shift);*/

			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				double f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					double r = (0.5 + l) / df->m_fFac;
//>>> SG 27.08.2019
					if ( m_bLorch )
						f += r * (df->m_pBin[l] - shift) / q * sin(r * q) * sin(M_PI * r / df->m_iResolution) * df->m_iResolution / (M_PI * r);
					else
//<<< SG 27.08.2019
					f += r * (df->m_pBin[l] - shift) / q * sin(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0 : 2.0) * ((double)m_isotopeTypeTotalCount[i]) * ((double)m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount);
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0 : 2.0) * ((double)m_isotopeTypeTotalCount[i]) * ((double)m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0 : 2.0) * ((double)m_isotopeTypeTotalCount[i]) * ((double)m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0);
				else
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_%s%s_sfac_%s_%s.csv", m_global ? "" : "self_", (const char *)m_name, (const char *)((CIsotope *)m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_isotopeTypeList[j])->label());
				mprintf("      Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}


void CStructureFactorGroup::finalizeIntra(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch) {
	int i, j, k, l;
	CxString filename;
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_isotopeTypeList.GetSize(); i++) {
		for (j = i; j < m_isotopeTypeList.GetSize(); j++) {
			mprintf("      Processing intramolecular partial structure factor %s-%s\n", ((CIsotope *)m_isotopeTypeList[i])->label(), ((CIsotope *)m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfIntraList[(m_isotopeTypeList.GetSize()-1)*i - i*(i-1)/2 + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if(i == j)
				df->MultiplyBin(2.0);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0 / (4.0/3.0 * Pi) / ((double)m_isotopeTypeTotalCount[i]) / ((double)m_isotopeTypeTotalCount[j]));
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0/3.0 * Pi) / ((double)m_isotopeTypeTotalCount[i]) / ((double)m_isotopeTypeTotalCount[j]));
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
				filename.sprintf("sfac_self_%s_rdf_%s_%s_intra.csv", (const char *)m_name, (const char *)(((CIsotope *)m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_isotopeTypeList[j])->label()));
				mprintf("      Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			double shift = 0.0;
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				double f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					double r = (0.5 + l) / df->m_fFac;
//>>> SG 27.08.2019
					if ( m_bLorch )
						f += r * (df->m_pBin[l] - shift) / q * sin(r * q) * sin(M_PI * r / df->m_iResolution) * df->m_iResolution / (M_PI * r);
					else
//<<< SG 27.08.2019
					f += r * (df->m_pBin[l] - shift) / q * sin(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0 : 2.0) * ((double)m_isotopeTypeTotalCount[i]) * ((double)m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount);
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0 : 2.0) * ((double)m_isotopeTypeTotalCount[i]) * ((double)m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0 : 2.0) * ((double)m_isotopeTypeTotalCount[i]) * ((double)m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0);
				else
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_self_%s_sfac_%s_%s_intra.csv", (const char *)m_name, (const char *)((CIsotope *)m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_isotopeTypeList[j])->label());
				mprintf("      Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}


void CStructureFactorGroup::finalizeInter(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch) {
	int i, j, k, l;
	CxString filename;
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_isotopeTypeList.GetSize(); i++) {
		for (j = i; j < m_isotopeTypeList.GetSize(); j++) {
			mprintf("      Processing intermolecular partial structure factor %s-%s\n", ((CIsotope *)m_isotopeTypeList[i])->label(), ((CIsotope *)m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfInterList[(m_isotopeTypeList.GetSize()-1)*i - i*(i-1)/2 + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if(i == j)
				df->MultiplyBin(2.0);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0 / (4.0/3.0 * Pi) / ((double)m_isotopeTypeTotalCount[i]) / ((double)m_isotopeTypeTotalCount[j]));
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0/3.0 * Pi) / ((double)m_isotopeTypeTotalCount[i]) / ((double)m_isotopeTypeTotalCount[j]));
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
				filename.sprintf("sfac_self_%s_rdf_%s_%s_inter.csv", (const char *)m_name, (const char *)(((CIsotope *)m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_isotopeTypeList[j])->label()));
				mprintf("      Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			//double shift = (double)m_isotopeTypeCount[i] * m_isotopeTypeCount[j] / (m_isotopeTypeTotalCount[i] * m_isotopeTypeTotalCount[j]);
			double shift = ((double)m_isotopeTypeCount[i]) * ((double)m_isotopeTypeCount[j]) / (((double)m_isotopeTypeTotalCount[i]) * ((double)m_isotopeTypeTotalCount[j]));
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				double f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					double r = (0.5 + l) / df->m_fFac;
//>>> SG 27.08.2019
					if ( m_bLorch )
						f += r * (df->m_pBin[l] - shift) / q * sin(r * q) * sin(M_PI * r / df->m_iResolution) * df->m_iResolution / (M_PI * r);
					else
//<<< SG 27.08.2019
					f += r * (df->m_pBin[l] - shift) / q * sin(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0 : 2.0) * ((double)m_isotopeTypeTotalCount[i]) * ((double)m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount);
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0 : 2.0) * ((double)m_isotopeTypeTotalCount[i]) * ((double)m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * ((i == j) ? 1.0 : 2.0) * ((double)m_isotopeTypeTotalCount[i]) * ((double)m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0);
				else
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_self_%s_sfac_%s_%s_inter.csv", (const char *)m_name, (const char *)((CIsotope *)m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_isotopeTypeList[j])->label());
				mprintf("      Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}


CStructureFactorCross::CStructureFactorCross(CStructureFactorGroup *sfacGroup1, CStructureFactorGroup *sfacGroup2) {
	m_sfacGroup1 = sfacGroup1;
	m_sfacGroup2 = sfacGroup2;
	m_name.sprintf("%s_%s", (const char *)sfacGroup1->getName(), (const char *)sfacGroup2->getName());
	m_sepInterIntra = false;
}


CStructureFactorCross::~CStructureFactorCross() {
	int i;
	for (i = 0; i < m_rdfList.GetSize(); i++)
		delete (CDF *)m_rdfList[i];
	for (i = 0; i < m_rdfIntraList.GetSize(); i++)
		delete (CDF *)m_rdfIntraList[i];
	for (i = 0; i < m_rdfInterList.GetSize(); i++)
		delete (CDF *)m_rdfInterList[i];
}


void CStructureFactorCross::initialize(double rdfMax, int rdfRes) {
	int i;
	mprintf("      Creating %d RDFs\n", m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize());
	for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize(); i++) {
		CDF *df;
		try { df = new CDF(); } catch(...) { df = NULL; }
		if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		df->m_fMinVal = 0.0;
		df->m_fMaxVal = rdfMax;
		df->m_iResolution = rdfRes;
		df->SetLabelX("Distance / pm");
		df->SetLabelY("g(r)");
		df->Create();
		m_rdfList.Add(df);
	}
	
	if (m_sepInterIntra) {
		mprintf("      Creating %d intramolecular RDFs\n", m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize());
		for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize(); i++) {
			CDF *df;
			try { df = new CDF(); } catch(...) { df = NULL; }
			if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			df->m_fMinVal = 0.0;
			df->m_fMaxVal = rdfMax;
			df->m_iResolution = rdfRes;
			df->SetLabelX("Distance / pm");
			df->SetLabelY("g(r)");
			df->Create();
			m_rdfIntraList.Add(df);
		}
		mprintf("      Creating %d intermolecular RDFs\n", m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize());
		for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize() * m_sfacGroup2->m_isotopeTypeList.GetSize(); i++) {
			CDF *df;
			try { df = new CDF(); } catch(...) { df = NULL; }
			if (df == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			df->m_fMinVal = 0.0;
			df->m_fMaxVal = rdfMax;
			df->m_iResolution = rdfRes;
			df->SetLabelX("Distance / pm");
			df->SetLabelY("g(r)");
			df->Create();
			m_rdfInterList.Add(df);
		}
	}
}


void CStructureFactorCross::process(CTimeStep* ts) {
	int i, j;
	for(i = 0; i < m_sfacGroup1->m_atomIndexList.GetSize(); i++) {
		//if ((i%100) == 0)
		//	mprintf("\n@Cross i=%d",i);
		for(j = 0; j < m_sfacGroup2->m_atomIndexList.GetSize(); j++) {
			double dist = FoldedLength(ts->m_vaCoords[m_sfacGroup1->m_atomIndexList[i]] - ts->m_vaCoords[m_sfacGroup2->m_atomIndexList[j]]);
			int a = m_sfacGroup1->m_isotopeList[i];
			int b = m_sfacGroup2->m_isotopeList[j];
			((CDF *)m_rdfList[a * m_sfacGroup2->m_isotopeTypeList.GetSize() + b])->AddToBin(dist);
			if (m_sepInterIntra) {
				if (m_sfacGroup1->m_singleMolList[i] == m_sfacGroup2->m_singleMolList[j])
					((CDF *)m_rdfIntraList[a * m_sfacGroup2->m_isotopeTypeList.GetSize() + b])->AddToBin(dist);
				else
					((CDF *)m_rdfInterList[a * m_sfacGroup2->m_isotopeTypeList.GetSize() + b])->AddToBin(dist);
			}
		}
	}
}


void CStructureFactorCross::finalize(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch) {
	int i, j, k, l;
	CxString filename;
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize(); i++) {
		for (j = 0; j < m_sfacGroup2->m_isotopeTypeList.GetSize(); j++) {
			mprintf("      Processing partial structure factor %s-%s\n", ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfList[i * m_sfacGroup2->m_isotopeTypeList.GetSize() + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				df->MultiplyBin(2.0);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0 / (4.0/3.0 * Pi) / ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) / ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]));
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0/3.0 * Pi) / ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) / ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]));
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
				filename.sprintf("sfac_cross_%s_rdf_%s_%s.csv", (const char *)m_name, (const char *)(((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label()));
				mprintf("      Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			double shift = ((double)m_sfacGroup1->m_isotopeTypeCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeCount[j]) / (((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]));
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				shift *= 2.0;
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				double f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					double r = (0.5 + l) / df->m_fFac;
//>>> SG 27.08.2019
                                        if ( m_bLorch )
                                                f += r * (df->m_pBin[l] - shift) / q * sin(r * q) * sin(M_PI * r / df->m_iResolution) * df->m_iResolution / (M_PI * r);
                                        else
//<<< SG 27.08.2019
					f += r * (df->m_pBin[l] - shift) / q * sin(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				double f = (CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j] ? 1.0 : 2.0;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * f * ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount);
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * f * ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * f * ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0);
				else
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_cross_%s_sfac_%s_%s.csv", (const char *)m_name, (const char *)((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
				mprintf("      Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}


void CStructureFactorCross::finalizeIntra(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch) {
	int i, j, k, l;
	CxString filename;
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize(); i++) {
		for (j = 0; j < m_sfacGroup2->m_isotopeTypeList.GetSize(); j++) {
			mprintf("      Processing intramolecular partial structure factor %s-%s\n", ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfIntraList[i * m_sfacGroup2->m_isotopeTypeList.GetSize() + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				df->MultiplyBin(2.0);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0 / (4.0/3.0 * Pi) / ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) / ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]));
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0/3.0 * Pi) / ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) / ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]));
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
				filename.sprintf("sfac_cross_%s_rdf_%s_%s_intra.csv", (const char *)m_name, (const char *)(((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label()));
				mprintf("      Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			double shift = 0.0;
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				shift *= 2.0;
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				double f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					double r = (0.5 + l) / df->m_fFac;
//>>> SG 27.08.2019
                                        if ( m_bLorch )
                                                f += r * (df->m_pBin[l] - shift) / q * sin(r * q) * sin(M_PI * r / df->m_iResolution) * df->m_iResolution / (M_PI * r);
                                        else
//<<< SG 27.08.2019
					f += r * (df->m_pBin[l] - shift) / q * sin(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				double f = (CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j] ? 1.0 : 2.0;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * f * ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount);
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * f * ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * f * ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0);
				else
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_cross_%s_sfac_%s_%s_intra.csv", (const char *)m_name, (const char *)((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
				mprintf("      Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}


void CStructureFactorCross::finalizeInter(CDF *sfac, CDF *xray, CDF *neutron, bool saveIntermediate, bool m_bLorch) {
	int i, j, k, l;
	CxString filename;
	
	CDF *sfacTemp;
	try { sfacTemp = new CDF(); } catch(...) { sfacTemp = NULL; }
	if(sfacTemp == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfacTemp->m_fMinVal = sfac->m_fMinVal;
	sfacTemp->m_fMaxVal = sfac->m_fMaxVal;
	sfacTemp->m_iResolution = sfac->m_iResolution;
	sfacTemp->SetLabelX("Wave vector modulus / nm^-1");
	sfacTemp->SetLabelY("Intensity");
	sfacTemp->Create();
	
	for (i = 0; i < m_sfacGroup1->m_isotopeTypeList.GetSize(); i++) {
		for (j = 0; j < m_sfacGroup2->m_isotopeTypeList.GetSize(); j++) {
			mprintf("      Processing intermolecular partial structure factor %s-%s\n", ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
			CDF *df = (CDF *)m_rdfInterList[i * m_sfacGroup2->m_isotopeTypeList.GetSize() + j];
			
			df->MultiplyBin(1.0 / g_iSteps * g_iStride);
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				df->MultiplyBin(2.0);
			df->CorrectRadialDist();
			if (g_bBoxNonOrtho)
				df->MultiplyBin(g_fBoxVolume * 1000000.0 / (4.0/3.0 * Pi) / ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) / ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]));
			else
				df->MultiplyBin(g_fBoxX * g_fBoxY * g_fBoxZ / (4.0/3.0 * Pi) / ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) / ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]));
			if (g_bDoubleBox)
				df->MultiplyBin(g_iDoubleBoxFactor);
			
			if (saveIntermediate) {
				filename.sprintf("sfac_cross_%s_rdf_%s_%s_inter.csv", (const char *)m_name, (const char *)(((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label()), (const char *)(((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label()));
				mprintf("      Writing RDF to %s...\n", (const char *)filename);
				df->Write("", filename, "", false);
			}
			
			double shift = ((double)m_sfacGroup1->m_isotopeTypeCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeCount[j]) / (((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]));
			if ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])
				shift *= 2.0;
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				double f = 0;
				for (l = 0; l < df->m_iResolution; l++) {
					double r = (0.5 + l) / df->m_fFac;
//>>> SG 27.08.2019
					if ( m_bLorch )
						f += r * (df->m_pBin[l] - shift) / q * sin(r * q) * sin(M_PI * r / df->m_iResolution) * df->m_iResolution / (M_PI * r);
					else
//<<< SG 27.08.2019
					f += r * (df->m_pBin[l] - shift) / q * sin(r * q);
				}
				sfacTemp->m_pBin[k] = f / df->m_fFac;
			}
			
			for (k = 0; k < sfacTemp->m_iResolution; k++) {
				double q = (sfacTemp->m_fMinVal + (0.5 + k) / sfacTemp->m_fFac) / 1000.0;
				double f = (CIsotope *)m_sfacGroup1->m_isotopeTypeList[i] == (CIsotope *)m_sfacGroup2->m_isotopeTypeList[j] ? 1.0 : 2.0;
				sfac->m_pBin[k] += sfacTemp->m_pBin[k] * f * ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount);
				neutron->m_pBin[k] += sfacTemp->m_pBin[k] * f * ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->neutronFactor() * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->neutronFactor();
				xray->m_pBin[k] += sfacTemp->m_pBin[k] * f * ((double)m_sfacGroup1->m_isotopeTypeTotalCount[i]) * ((double)m_sfacGroup2->m_isotopeTypeTotalCount[j]) / ((double)g_iGesAtomCount) / ((double)g_iGesAtomCount) * ((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->xrayFactor(q) * ((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->xrayFactor(q);
			}
			
			if (saveIntermediate) {
				if (g_bBoxNonOrtho)
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / g_fBoxVolume / 1000000.0);
				else
					sfacTemp->MultiplyBin(4.0 * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ));
				filename.sprintf("sfac_cross_%s_sfac_%s_%s_inter.csv", (const char *)m_name, (const char *)((CIsotope *)m_sfacGroup1->m_isotopeTypeList[i])->label(), (const char *)((CIsotope *)m_sfacGroup2->m_isotopeTypeList[j])->label());
				mprintf("      Writing Structure Factor to %s...\n", (const char *)filename);
				sfacTemp->Write("", filename, "", false);
			}
			
			sfacTemp->ZeroBin();
		}
	}
	
	delete sfacTemp;
}


CStructureFactor::CStructureFactor() {
	int i, j, k;
	CxString buf, buf2, temp;
	
	CxObArray isotopeAssignList;
	for (i = 0; i < g_oaMolecules.GetSize(); i++) {
		CMolecule *mol = (CMolecule *)g_oaMolecules[i];
		StructureFactorMolecule *smol;
		try { smol = new StructureFactorMolecule; } catch(...) { smol = NULL; }
		if (smol == NULL) NewException((double)sizeof(StructureFactorMolecule), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		smol->moleculeType = i;
		for (j = 0; j < mol->m_baAtomIndex.GetSize(); j++) {
			if (mol->m_baAtomIndex[j] == g_iVirtAtomType)
				continue;
			StructureFactorAtomKind *satom;
			try { satom = new StructureFactorAtomKind; } catch(...) { satom = NULL; }
			if (satom == NULL) NewException((double)sizeof(StructureFactorAtomKind), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			satom->atomType = mol->m_baAtomIndex[j];
			CIsotope *isotope = NULL;
			for (k = 0; k < g_isotopes.GetSize(); k++) {
				if (mystricmp(((CAtom *)g_oaAtoms[mol->m_baAtomIndex[j]])->m_sName, ((CIsotope *)g_isotopes[k])->label()) == 0) {
					isotope = (CIsotope *)g_isotopes[k];
					break;
				}
			}
			if (isotope == NULL) {
				mprintf(RED, "Warning: no isotope data for \"%s\" found. Setting all scattering factors to 0.\n", (const char*)((CAtom *)g_oaAtoms[mol->m_baAtomIndex[j]])->m_sName);
				try { isotope = new CIsotope(((CAtom *)g_oaAtoms[mol->m_baAtomIndex[j]])->m_sName, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0); } catch(...) { isotope = NULL; }
				if (isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				g_isotopes.Add(isotope);
			}
			for (k = 0; k < mol->m_waAtomCount[j]; k++) {
				satom->isotopeList.Add(isotope);
			}
			smol->atomKinds.Add(satom);
		}
		isotopeAssignList.Add(smol);
	}
	
	if (g_bAdvanced2) {
		if (!AskYesNo("    Use standard atom data (y) or specify isotopes (n)? [yes] ", true)) {
			while (true) {
				buf.sprintf("\n    Change isotopes in which molecule (");
				for(i = 0; i < g_oaMolecules.GetSize(); i++) {
					buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
					buf.strcat(buf2);
					if(i < g_oaMolecules.GetSize() - 1)
						buf.strcat(", ");
				}
				buf.strcat(")? ");
				int mol = AskRangeInteger_ND("%s", 1, g_oaMolecules.GetSize(),(const char*)buf) - 1;
				CMolecule *m = (CMolecule *)g_oaMolecules[mol];
				
				while (true) {
					mprintf("\n    The following isotopes are set up:\n");
					
					for (i = 0; i < m->m_baAtomIndex.GetSize(); i++) {
						if (m->m_baAtomIndex[i] == g_iVirtAtomType)
							continue;
						for (j = 0; j < m->m_waAtomCount[i]; j++) {
							mprintf("      %s%d: %-4s", (const char*)((CAtom *)g_oaAtoms[m->m_baAtomIndex[i]])->m_sName, j+1, ((CIsotope *)((StructureFactorAtomKind *)((StructureFactorMolecule *)isotopeAssignList[mol])->atomKinds[i])->isotopeList[j])->label());
						}
					}
					mprintf("\n\n");
					
					int ca, cb, cc;
					ca = 0;
					cc = 0;
					do {
						AskString("    Change isotope for which atom? [done] ", &buf, "");
						if (buf.GetLength() == 0)
							break;
					} while (!ParseAtom(buf, mol, ca, cb, cc));
					if (buf.GetLength() == 0)
						break;
					
					while (true) {
						AskString_ND("    Which isotope to set for %s (e.g. \"13C\")? ", &buf2, (const char *)buf);
						int isotopeIndex = -1;
						for (i = 0; i < g_isotopes.GetSize(); i++) {
							if(mystricmp(((CIsotope *)g_isotopes[i])->label(), buf2) == 0) {
								isotopeIndex = i;
								break;
							}
						}
						if (isotopeIndex != -1) {
							((StructureFactorAtomKind *)((StructureFactorMolecule *)isotopeAssignList[mol])->atomKinds[ca])->isotopeList[cc] = (CIsotope *)g_isotopes[isotopeIndex];
							break;
						}
						mprintf("\n    No isotope data for \"%s\" found.\n", (const char *)buf2);
						if (AskYesNo("    Enter data for \"%s\" (y/n)? [yes] ", true, (const char *)buf2)) {
							double nf = AskFloat_ND("    Enter neutron scattering factor: ");
							double cma1 = AskFloat_ND("    Enter Cromer-Mann coefficient a1: ");
							double cma2 = AskFloat_ND("    Enter Cromer-Mann coefficient a2: ");
							double cma3 = AskFloat_ND("    Enter Cromer-Mann coefficient a3: ");
							double cma4 = AskFloat_ND("    Enter Cromer-Mann coefficient a4: ");
							double cmb1 = AskFloat_ND("    Enter Cromer-Mann coefficient b1: ");
							double cmb2 = AskFloat_ND("    Enter Cromer-Mann coefficient b2: ");
							double cmb3 = AskFloat_ND("    Enter Cromer-Mann coefficient b3: ");
							double cmb4 = AskFloat_ND("    Enter Cromer-Mann coefficient b4: ");
							double cmc = AskFloat_ND("    Enter Cromer-Mann coefficient c: ");
							CIsotope *isotope;
							try { isotope = new CIsotope(buf2, nf, cma1, cma2, cma3, cma4, cmb1, cmb2, cmb3, cmb4, cmc); } catch(...) { isotope = NULL; }
							if (isotope == NULL) NewException((double)sizeof(CIsotope), __FILE__, __LINE__, __PRETTY_FUNCTION__);
							g_isotopes.Add(isotope);
							((StructureFactorAtomKind *)((StructureFactorMolecule *)isotopeAssignList[mol])->atomKinds[ca])->isotopeList[cc] = isotope;
							break;
						} else {
							mprintf("\n");
						}
					}
				}
				
				if (!AskYesNo("\n    Change isotopes in another molecule (y/n)? [no] ", false))
					break;
			}
			mprintf("\n");
		} else {
			mprintf("\n");
		}
		
		if (AskYesNo("    Change the predefined scattering factors (y/n)? [no] ", false)) {
			CxObArray isotopeTypeList;
			for (i = 0; i < isotopeAssignList.GetSize(); i++) {
				StructureFactorMolecule *smol = (StructureFactorMolecule *)isotopeAssignList[i];
				for (j = 0; j < smol->atomKinds.GetSize(); j++) {
					StructureFactorAtomKind *satom = (StructureFactorAtomKind *)smol->atomKinds[j];
					for (k = 0; k < satom->isotopeList.GetSize(); k++) {
						CIsotope *isotope = (CIsotope *)satom->isotopeList[k];
						bool found = false;
						int l;
						for (l = 0; l < isotopeTypeList.GetSize(); l++) {
							if ((CIsotope *)isotopeTypeList[l] == isotope) {
								found = true;
								break;
							}
						}
						if (!found) {
							isotopeTypeList.Add(isotope);
						}
					}
				}
			}
			
			while (true) {
				mprintf("\n    The following scattering factors are set up:\n\n");
				mprintf("    Name    neutron       a1       a2       a3       a4       b1       b2       b3       b4        c\n");
				for (i = 0; i < isotopeTypeList.GetSize(); i++) {
					CIsotope *isotope = (CIsotope *)isotopeTypeList[i];
					mprintf("    %-6s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", isotope->_label, isotope->_neutronFactor, isotope->_cma[0], isotope->_cma[1], isotope->_cma[2], isotope->_cma[3], isotope->_cmb[0], isotope->_cmb[1], isotope->_cmb[2], isotope->_cmb[3], isotope->_cmc);
				}
				mprintf("\n");
				CIsotope *isotope = NULL;
				while (true) {
					AskString("    Change scattering factors of which isotope? [done] ", &buf, "");
					if (buf.GetLength() == 0)
						break;
					for (i = 0; i < isotopeTypeList.GetSize(); i++) {
						if (mystricmp(((CIsotope *)isotopeTypeList[i])->label(), buf) == 0) {
							isotope = (CIsotope *)isotopeTypeList[i];
							break;
						}
					}
					if (isotope != NULL)
						break;
					eprintf("The system does not contain isotope \"%s\".\n", (const char *)buf);
				}
				if (buf.GetLength() == 0)
					break;
				isotope->_neutronFactor = AskFloat_ND("    Enter neutron scattering factor: ");
				isotope->_cma[0] = AskFloat_ND("    Enter Cromer-Mann coefficient a1: ");
				isotope->_cma[1] = AskFloat_ND("    Enter Cromer-Mann coefficient a2: ");
				isotope->_cma[2] = AskFloat_ND("    Enter Cromer-Mann coefficient a3: ");
				isotope->_cma[3] = AskFloat_ND("    Enter Cromer-Mann coefficient a4: ");
				isotope->_cmb[0] = AskFloat_ND("    Enter Cromer-Mann coefficient b1: ");
				isotope->_cmb[1] = AskFloat_ND("    Enter Cromer-Mann coefficient b2: ");
				isotope->_cmb[2] = AskFloat_ND("    Enter Cromer-Mann coefficient b3: ");
				isotope->_cmb[3] = AskFloat_ND("    Enter Cromer-Mann coefficient b4: ");
				isotope->_cmc = AskFloat_ND("    Enter Cromer-Mann coefficient c: ");
			}
			mprintf("\n");
		} else {
			mprintf("\n");
		}
	}

	mprintf("\n    The following scattering factors are set up:\n\n");
	mprintf("    Name    neutron       a1       a2       a3       a4       b1       b2       b3       b4        c\n");
	CxObArray tlist;
	for (i = 0; i < isotopeAssignList.GetSize(); i++) {
		StructureFactorMolecule *smol = (StructureFactorMolecule *)isotopeAssignList[i];
		for (j = 0; j < smol->atomKinds.GetSize(); j++) {
			StructureFactorAtomKind *satom = (StructureFactorAtomKind *)smol->atomKinds[j];
			for (k = 0; k < satom->isotopeList.GetSize(); k++) {
				CIsotope *isotope = (CIsotope *)satom->isotopeList[k];
				bool tfound = false;
				int tl;
				for (tl = 0; tl < tlist.GetSize(); tl++) {
					if ((CIsotope *)tlist[tl] == isotope) {
						tfound = true;
						break;
					}
				}
				if (!tfound) {
					mprintf("    %-6s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", isotope->_label, isotope->_neutronFactor, isotope->_cma[0], isotope->_cma[1], isotope->_cma[2], isotope->_cma[3], isotope->_cmb[0], isotope->_cmb[1], isotope->_cmb[2], isotope->_cmb[3], isotope->_cmc);
					tlist.Add(isotope);
				}
			}
		}
	}
	mprintf("\n");

	if (!g_bAdvanced2)
		mprintf("\n    To specify isotopes or change the scattering cross sections, use the advanced mode.\n\n");
	
	try { m_globalSFac = new CStructureFactorGroup(true, isotopeAssignList); } catch(...) { m_globalSFac = NULL; }
	if (m_globalSFac == NULL) NewException((double)sizeof(CStructureFactorGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	if (AskYesNo("    Compute partial structure factors of specific atom groups (y/n)? [no] ", false)) {
		while (true) {
			mprintf(WHITE, "\n    > Group %d >\n\n", m_sFacGroups.GetSize() + 1);
			CStructureFactorGroup *sfac;
			try { sfac = new CStructureFactorGroup(false, isotopeAssignList); } catch(...) { sfac = NULL; }
			if (sfac == NULL) NewException((double)sizeof(CStructureFactorGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			m_sFacGroups.Add(sfac);
			AskString("      Which label to use for this atom group? [%s] ",&buf,(const char*)sfac->getName(),(const char*)sfac->getName());
			sfac->getName() = (const char*)buf;
			mprintf("\n");
			mprintf(WHITE, "    < End of Group %d <\n\n", m_sFacGroups.GetSize());
			if (!AskYesNo("    Add another group (y/n)? [no] ", false))
				break;
		}
	}
	mprintf("\n");
	
	if (m_sFacGroups.GetSize() > 1) {
		if (AskYesNo("    Compute all cross terms of the defined atom groups (y) or select certain combinations (n)? [yes] ", true)) {
			mprintf("\n");
			int count = 0;
			for (i = 0; i < m_sFacGroups.GetSize(); i++) {
				for (j = i + 1; j < m_sFacGroups.GetSize(); j++) {
					CStructureFactorCross *sfac;
					try { sfac = new CStructureFactorCross((CStructureFactorGroup *)m_sFacGroups[i], (CStructureFactorGroup *)m_sFacGroups[j]); } catch(...) { sfac = NULL; }
					if (sfac == NULL) NewException((double)sizeof(CStructureFactorGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
					m_sFacCrosses.Add(sfac);
					mprintf("    Cross term %d: Group %d - Group %d\n", ++count, i + 1, j + 1);
				}
			}
		} else {
			mprintf("    The following groups are defined:\n\n");
			for (i = 0; i < m_sFacGroups.GetSize(); i++) {
				mprintf("    %3d - %s\n", i + 1, (const char *)((CStructureFactorGroup *)m_sFacGroups[i])->getName());
			}
			mprintf("\n    Enter each combination as a comma-separated pair of numbers (e.g. 1,1)\n\n");
			CxWordArray wa;
			while (true) {
				wa.RemoveAll_KeepSize();
				AskString("    Enter combination (return=finished): ", &buf, "");
				if (buf.GetLength() == 0)
					break;
				ParseIntList(buf, &wa);
				if (wa.GetSize() != 2) {
					eprintf("Invalid input. Please enter exactly two numbers separated by a comma.\n");
					continue;
				}
				if ((wa[0] < 1) || (wa[1] < 1)) {
					eprintf("Invalid input. Please enter positive numbers.\n");
					continue;
				}
				if ((wa[0] > m_sFacGroups.GetSize()) || (wa[1] > m_sFacGroups.GetSize())) {
					eprintf("Invalid input. There are only %d groups.\n", m_sFacGroups.GetSize());
					continue;
				}
				if (wa[0] == wa[1]) {
					eprintf("Invalid input. Please enter two different numbers.\n");
					continue;
				}
				CStructureFactorCross *sfac;
				try { sfac = new CStructureFactorCross((CStructureFactorGroup *)m_sFacGroups[wa[0] - 1], (CStructureFactorGroup *)m_sFacGroups[wa[1] - 1]); } catch(...) { sfac = NULL; }
				if (sfac == NULL) NewException((double)sizeof(CStructureFactorGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
				m_sFacCrosses.Add(sfac);
			}
		}
		mprintf("\n");
	}
	
	if (m_sFacCrosses.GetSize() > 0) {
		if (AskYesNo("    Separate intramolecular and intermolecular contributions in the cross terms (y/n)? [no] ", false))
			for (i = 0; i < m_sFacCrosses.GetSize(); i++)
				((CStructureFactorCross *)m_sFacCrosses[i])->setSepInterIntra(true);
		mprintf("\n");
	}

	m_fNormalizationFactor = 1.0;
	
	m_rdfMax = AskFloat("    Enter maximum RDF distance to observe (pm): [%d] ", (double)HalfBox(), HalfBox());
	m_rdfRes = AskUnsignedInteger("    Enter RDF binning resolution: [%d] ", 2 * HalfBox(), 2 * HalfBox());
	m_sfacMax = AskFloat("    Enter maximum wave vector modulus (nm^-1): [200] ", 200.0);
	m_sfacRes = AskUnsignedInteger("    Enter Structure Factor resolution: [2000] ", 2000);
	
	mprintf("\n    The following normalization factors are available:\n");
	mprintf("    (1) [Sum_i x_i*f_i(q)]^2\n");
	mprintf("    (2) Sum_i x_i*[f_i(q)]^2\n");
	mprintf("    (3) manual normalization factor\n");
	mprintf("    (4) 1 (no normalization)\n\n");

	mprintf("    They will all be written to the output files as consecutive columns.\n\n");

	m_fNormalizationFactor = AskFloat("    Enter manual normalization factor for (3): [%.5f] ",0.00714,0.00714);
	
	if (g_bAdvanced2) {
		mprintf("\n");
		m_saveIntermediate = AskYesNo("    Save all intermediate data (y/n)? [no] ", false);
	} else
		m_saveIntermediate = false;

	mprintf("\n");
	
	m_sharpening = AskYesNo("    Apply sharpening factor (y/n)? [no] ", false);
	if (m_sharpening) {
		while (true) {
			AskString_ND("    Which isotope to use as sharpening atom (e.g. \"N\" or \"14N\")? ", &buf);
			int isotopeIndex = -1;
			for (i = 0; i < g_isotopes.GetSize(); i++) {
				if(mystricmp(((CIsotope *)g_isotopes[i])->label(), buf) == 0) {
					isotopeIndex = i;
					break;
				}
			}
			if(isotopeIndex != -1) {
				m_sharpeningIsotope = (CIsotope *)g_isotopes[isotopeIndex];
				break;
			}
			eprintf("Isotope data for \"%s\" not found.\n", (const char *)buf);
		}
	} else {
		m_sharpeningIsotope = NULL;
	}
	mprintf("\n");

//>>> SG 27.08.2019
	m_bLorch = AskYesNo("    Apply Lorch type window function (y/n)? [no] ", false);
//<<< SG 27.08.2019
}


CStructureFactor::~CStructureFactor() {
	int i;
	delete m_globalSFac;
	for (i = 0; i < m_sFacGroups.GetSize(); i++)
		delete (CStructureFactorGroup *)m_sFacGroups[i];
	for (i = 0; i < m_sFacCrosses.GetSize(); i++)
		delete (CStructureFactorCross *)m_sFacCrosses[i];
}


void CStructureFactor::initialize() {
	int i;
	mprintf("  Initializing structure factor\n");
	mprintf(WHITE, "\n    > Total >\n\n");
	m_globalSFac->initialize(m_rdfMax, m_rdfRes);
	mprintf(WHITE, "\n    < End of Total <\n");
	for (i = 0; i < m_sFacGroups.GetSize(); i++) {
		mprintf(WHITE, "\n    > Group %d >\n\n", i + 1);
		((CStructureFactorGroup *)m_sFacGroups[i])->initialize(m_rdfMax, m_rdfRes);
		mprintf(WHITE, "\n    < End of Group %d <\n", i + 1);
	}
	for (i = 0; i < m_sFacCrosses.GetSize(); i++) {
		mprintf(WHITE, "\n    > Cross Term %d >\n\n", i + 1);
		((CStructureFactorCross *)m_sFacCrosses[i])->initialize(m_rdfMax, m_rdfRes);
		mprintf(WHITE, "\n    < End of Cross Term %d <\n", i + 1);
	}
	#ifdef _OPENMP
		#ifdef USE_OMP
			mprintf( "\n  OpenMP is active, max_threads=%d\n\n", omp_get_max_threads() );
		#else
			mprintf( RED, "\n  Compiled with OpenMP, but USE_OMP not switched on in \"config.h\"!\n\n" );
		#endif
	#else
		#ifdef USE_OMP
			mprintf( RED, "\n  USE_OMP is switched on, but code not compiled with -fopenmp!\n\n" );
		#endif
	#endif
}


void CStructureFactor::process(CTimeStep *ts) {
	int i;
	m_globalSFac->process(ts);
	for (i = 0; i < m_sFacGroups.GetSize(); i++) {
		((CStructureFactorGroup *)m_sFacGroups[i])->process(ts);
	}
	for (i = 0; i < m_sFacCrosses.GetSize(); i++) {
		((CStructureFactorCross *)m_sFacCrosses[i])->process(ts);
	}
}



void CStructureFactor::HandleGroup( const char *name, CStructureFactorGroup *sfgroup, CStructureFactorCross *sfcross, int mode ) {

	int i, z, z2, z3;
	CDF *sfac, *xray, *neutron, *dfwork;
	double *workbin, sf, xf, nf, q;
	CxString modesuffix, filename, clabel;


	try { sfac = new CDF(); } catch(...) { sfac = NULL; }
	if (sfac == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	sfac->m_fMinVal = 2.0 * Pi / m_rdfMax * 1000.0;
	sfac->m_fMaxVal = m_sfacMax;
	sfac->m_iResolution = m_sfacRes;
	sfac->SetLabelX("Wave vector modulus / nm^-1");
	sfac->SetLabelY("Intensity (normalization 1)");
	sfac->Create();

	try { xray = new CDF(); } catch(...) { xray = NULL; }
	if (xray == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	xray->m_fMinVal = 2.0 * Pi / m_rdfMax * 1000.0;
	xray->m_fMaxVal = m_sfacMax;
	xray->m_iResolution = m_sfacRes;
	xray->SetLabelX("Wave vector modulus / nm^-1");
	xray->SetLabelY("Intensity (normalization 1)");
	xray->Create();

	try { neutron = new CDF(); } catch(...) { neutron = NULL; }
	if (neutron == NULL) NewException((double)sizeof(CDF), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	neutron->m_fMinVal = 2.0 * Pi / m_rdfMax * 1000.0;
	neutron->m_fMaxVal = m_sfacMax;
	neutron->m_iResolution = m_sfacRes;
	neutron->SetLabelX("Wave vector modulus / nm^-1");
	neutron->SetLabelY("Intensity (normalization 1)");
	neutron->Create();
	
	switch(mode) {

		case 0: // Standard
			modesuffix = "";
			if (sfgroup != NULL)
				sfgroup->finalize( sfac, xray, neutron, m_saveIntermediate, m_bLorch );
			else if (sfcross != NULL)
				sfcross->finalize( sfac, xray, neutron, m_saveIntermediate, m_bLorch );
			else {
				eprintf("CStructureFactor::HandleGroup(): Internal error.\n");
				abort();
			}
			break;

		case 1: // Intra
			modesuffix.Format("_intra");
			if (sfgroup != NULL)
				sfgroup->finalizeIntra( sfac, xray, neutron, m_saveIntermediate, m_bLorch );
			else if (sfcross != NULL)
				sfcross->finalizeIntra( sfac, xray, neutron, m_saveIntermediate, m_bLorch );
			else {
				eprintf("CStructureFactor::HandleGroup(): Internal error.\n");
				abort();
			}
			break;

		case 2: // Inter
			modesuffix.Format("_inter");
			if (sfgroup != NULL)
				sfgroup->finalizeInter( sfac, xray, neutron, m_saveIntermediate, m_bLorch );
			else if (sfcross != NULL)
				sfcross->finalizeInter( sfac, xray, neutron, m_saveIntermediate, m_bLorch );
			else {
				eprintf("CStructureFactor::HandleGroup(): Internal error.\n");
				abort();
			}
			break;
	}

	for (i=0;i<3;i++) {

		switch(i) {
			case 0:
				dfwork = sfac;
				break;
			case 1:
				dfwork = xray;
				break;
			case 2:
				dfwork = neutron;
				break;
			default:
				eprintf("CStructureFactor::HandleGroup(): Internal error.\n");
				abort();
		}

		if (g_bBoxNonOrtho) {
			double f = 4.0 * Pi * g_iGesAtomCount / (g_fBoxVolume * 1000000.0);
			dfwork->MultiplyBin(f);
		} else {
			double f = 4.0 * Pi * g_iGesAtomCount / (g_fBoxX * g_fBoxY * g_fBoxZ);
			dfwork->MultiplyBin(f);
		}

		dfwork->m_iAdditionalSets = 3;

		try { dfwork->m_pAdditionalSets = new double*[dfwork->m_iAdditionalSets]; } catch(...) { dfwork->m_pAdditionalSets = NULL; }
		if (dfwork->m_pAdditionalSets == NULL) NewException((double)dfwork->m_iAdditionalSets*sizeof(double*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		try { dfwork->m_pAdditionalSetLabels = new char*[dfwork->m_iAdditionalSets]; } catch(...) { dfwork->m_pAdditionalSetLabels = NULL; }
		if (dfwork->m_pAdditionalSetLabels == NULL) NewException((double)dfwork->m_iAdditionalSets*sizeof(char*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z=0;z<dfwork->m_iAdditionalSets;z++) {

			try { dfwork->m_pAdditionalSets[z] = new double[dfwork->m_iResolution]; } catch(...) { dfwork->m_pAdditionalSets[z] = NULL; }
			if (dfwork->m_pAdditionalSets[z] == NULL) NewException((double)dfwork->m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			for (z2=0;z2<dfwork->m_iResolution;z2++)
				dfwork->m_pAdditionalSets[z][z2] = dfwork->m_pBin[z2];

			switch(z) {
				case 0:
					clabel.sprintf("Intensity (normalization 2)");
					break;
				case 1:
					clabel.sprintf("Intensity (manual normalization by %.8f)",m_fNormalizationFactor);
					break;
				case 2:
					clabel.sprintf("Intensity (no normalization)");
					break;
			}

			try { dfwork->m_pAdditionalSetLabels[z] = new char[strlen((const char*)clabel)+1]; } catch(...) { dfwork->m_pAdditionalSetLabels[z] = NULL; }
			if (dfwork->m_pAdditionalSetLabels[z] == NULL) NewException((double)strlen((const char*)clabel)+1,__FILE__,__LINE__,__PRETTY_FUNCTION__);
			strcpy( dfwork->m_pAdditionalSetLabels[z], (const char*)clabel );
		}

		for (z=0;z<4;z++) {

			if (z == 0)
				workbin = dfwork->m_pBin;
			else
				workbin = dfwork->m_pAdditionalSets[z-1];

			switch(z) {

				case 0: // Normalization factor 1
					switch(i) {
						case 0: // SFac
							sf = 0.0;
							for (z2=0;z2<m_globalSFac->getIsotopeTypeList().GetSize();z2++)
								sf += ((double)m_globalSFac->getIsotopeTypeTotalCount().GetAt(z2)) / ((double)g_iGesAtomCount);
							for (z2=0;z2<dfwork->m_iResolution;z2++)
								workbin[z2] *= 1.0 / (sf * sf);
							break;

						case 1: // Xray
							for (z2=0;z2<dfwork->m_iResolution;z2++) {
								q = (dfwork->m_fMinVal + (0.5 + z2) / dfwork->m_fFac) / 1000.0;
								xf = 0.0;
								for (z3=0;z3<m_globalSFac->getIsotopeTypeList().GetSize();z3++)
									xf += ((double)m_globalSFac->getIsotopeTypeTotalCount().GetAt(z3)) / ((double)g_iGesAtomCount) * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(z3))->xrayFactor(q);
								workbin[z2] *= 1.0 / (xf * xf);
							}
							break;
					
						case 2: // Neutron
							nf = 0.0;
							for (z2=0;z2<m_globalSFac->getIsotopeTypeList().GetSize();z2++)
								nf += ((double)m_globalSFac->getIsotopeTypeTotalCount().GetAt(z2)) / ((double)g_iGesAtomCount) * ((CIsotope *)m_globalSFac->getIsotopeTypeList().GetAt(z2))->neutronFactor();
							for (z2=0;z2<dfwork->m_iResolution;z2++)
								workbin[z2] *= 1.0 / (nf * nf);
							break;
					}
					break;

				case 2: // Normalization factor 2
					switch(i) {
						case 0: // SFac
							sf = 0.0;
							for (z2=0;z2<m_globalSFac->getIsotopeTypeList().GetSize();z2++)
								sf += ((double)m_globalSFac->getIsotopeTypeTotalCount().GetAt(z2)) / ((double)g_iGesAtomCount);
							for (z2=0;z2<dfwork->m_iResolution;z2++)
								workbin[z2] *= 1.0 / sf;
							break;

						case 1: // Xray
							for (z2=0;z2<dfwork->m_iResolution;z2++) {
								q = (dfwork->m_fMinVal + (0.5 + z2) / dfwork->m_fFac) / 1000.0;
								xf = 0.0;
								for (z3=0;z3<m_globalSFac->getIsotopeTypeList().GetSize();z3++)
									xf += ((double)m_globalSFac->getIsotopeTypeTotalCount().GetAt(z3)) / ((double)g_iGesAtomCount) * ((CIsotope*)m_globalSFac->getIsotopeTypeList().GetAt(z3))->xrayFactor(q) * ((CIsotope*)m_globalSFac->getIsotopeTypeList().GetAt(z3))->xrayFactor(q);
								workbin[z2] *= 1.0 / xf;
							}
							break;
					
						case 2: // Neutron
							nf = 0.0;
							for (z2=0;z2<m_globalSFac->getIsotopeTypeList().GetSize();z2++)
								nf += ((double)m_globalSFac->getIsotopeTypeTotalCount().GetAt(z2)) / ((double)g_iGesAtomCount) * ((CIsotope*)m_globalSFac->getIsotopeTypeList().GetAt(z2))->neutronFactor() * ((CIsotope*)m_globalSFac->getIsotopeTypeList().GetAt(z2))->neutronFactor();
							for (z2=0;z2<dfwork->m_iResolution;z2++)
								workbin[z2] *= 1.0 / nf;
							break;
					}
					break;

				case 3: // Manual Normalization
					for (z2=0;z2<dfwork->m_iResolution;z2++)
						workbin[z2] *= m_fNormalizationFactor;
					break;

				case 4: // No Normalization
					// Do nothing
					break;
			}
		}
	}

	mprintf("\n");
	filename.sprintf( "%s_unweighted%s.csv", name, (const char*)modesuffix );
	mprintf( "      Writing unweighted structure factor to %s ...\n", (const char *)filename );
	sfac->Write( "", filename, "", false );
	filename.sprintf( "%s_xray%s.csv", name, (const char*)modesuffix );
	mprintf( "      Writing X-ray scattering function to %s ...\n", (const char *)filename );
	xray->Write( "", filename, "", false );
	filename.sprintf( "%s_neutron%s.csv", name, (const char*)modesuffix );
	mprintf( "      Writing neutron scattering function to %s ...\n", (const char *)filename );
	neutron->Write( "", filename, "", false );
	
	if (m_sharpening) {

		for (i=0;i<xray->m_iResolution;i++) {
			q = (xray->m_fMinVal + (0.5 + i) / xray->m_fFac) / 1000.0;
			nf = q * exp(-0.01 * q * q);
			xf = q * exp(-0.01 * q * q) * m_sharpeningIsotope->xrayFactor(0.0) * m_sharpeningIsotope->xrayFactor(0.0) / m_sharpeningIsotope->xrayFactor(q) / m_sharpeningIsotope->xrayFactor(q);
			neutron->m_pBin[i] *= nf;
			xray->m_pBin[i] *= xf;
			for (z=0;z<3;z++) {
				neutron->m_pAdditionalSets[z][i] *= nf;
				xray->m_pAdditionalSets[z][i] *= xf;
			}
		}
		
		filename.sprintf( "%s_neutron%s_sharpened.csv", name,(const char*)modesuffix );
		mprintf( "      Writing sharpened neutron scattering function to %s ...\n", (const char *)filename );
		neutron->Write( "", filename, "", false );
		filename.sprintf( "%s_xray%s_sharpened.csv", name, (const char*)modesuffix );
		mprintf( "      Writing sharpened X-ray scattering function to %s ...\n", (const char *)filename );
		xray->Write( "", filename, "", false );
	}

	mprintf("\n");

	delete sfac;
	delete xray;
	delete neutron;
}



void CStructureFactor::finalize() {

	int i;
	CxString name;


	mprintf("  Finalizing structure factor\n");
	
	mprintf(WHITE, "\n    > Total >\n\n");
	
	HandleGroup( "sfac_total", m_globalSFac, NULL, 0 );

	if (m_globalSFac->sepInterIntra()) {
		HandleGroup( "sfac_total", m_globalSFac, NULL, 1 );
		HandleGroup( "sfac_total", m_globalSFac, NULL, 2 );
	}
		
	mprintf(WHITE, "    < End of Total <\n");
	
	for (i=0;i<m_sFacGroups.GetSize();i++) {

		mprintf(WHITE, "\n    > Group %d >\n\n", i + 1);

		name.sprintf( "sfac_self_%s", (const char*)((CStructureFactorGroup*)m_sFacGroups[i])->getName() );

		HandleGroup( (const char*)name, (CStructureFactorGroup*)m_sFacGroups[i], NULL, 0 );

		if (((CStructureFactorGroup*)m_sFacGroups[i])->sepInterIntra()) {
			HandleGroup( (const char*)name, (CStructureFactorGroup*)m_sFacGroups[i], NULL, 1 );
			HandleGroup( (const char*)name, (CStructureFactorGroup*)m_sFacGroups[i], NULL, 2 );
		}

		mprintf(WHITE, "    < End of Group %d <\n", i + 1);
	}

	for (i=0;i<m_sFacCrosses.GetSize();i++) {

		mprintf(WHITE, "\n    > Cross Term %d >\n\n", i + 1);

		name.sprintf( "sfac_cross_%s", (const char*)((CStructureFactorCross*)m_sFacCrosses[i])->getName() );

		HandleGroup( (const char*)name, NULL, (CStructureFactorCross*)m_sFacCrosses[i], 0 );

		if (((CStructureFactorGroup*)m_sFacGroups[i])->sepInterIntra()) {
			HandleGroup( (const char*)name, NULL, (CStructureFactorCross*)m_sFacCrosses[i], 1 );
			HandleGroup( (const char*)name, NULL, (CStructureFactorCross*)m_sFacCrosses[i], 2 );
		}

		mprintf(WHITE, "    < End of Cross Term %d <\n", i + 1);
	}
}


bool gatherStructureFactor() {
	createIsotopeList();
	
	try { g_structureFactor = new CStructureFactor(); } catch(...) { g_structureFactor = NULL; }
	if (g_structureFactor == NULL) NewException((double)sizeof(CStructureFactor), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	return true;
}


bool initializeStructureFactor() {
	g_structureFactor->initialize();
	return true;
}


void processStructureFactor(CTimeStep *ts) {
	g_structureFactor->process(ts);
}


void finalizeStructureFactor() {
	g_structureFactor->finalize();
	delete g_structureFactor;
	deleteIsotopeList();
}

