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

#include "roa.h"
#include "xstring.h"
#include "globalvar.h"
#include "maintools.h"
#include "bicgstab.h"
#include "conversion.h"


const char *GetRevisionInfo_roa(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_roa() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



void CSmoothener::Init(int length, double wavenumber) {

	m_iLength = length;
	m_fWaveNumber = wavenumber;

	m_pForward = new CFFT();
	// The constant "40" controls the padding between original data set and mirrored data set
	m_iInternalLength = m_pForward->NextFastSize(m_iLength+40);
	m_pForward->PrepareFFT_C2C(2*m_iInternalLength);
	m_pInverse = new CFFT();
	m_pInverse->PrepareInverseFFT_C2C(2*m_iInternalLength);
}


void CSmoothener::Smoothen(const std::vector<double> &inp, std::vector<double> &outp) {

	int z;


	for (z=0;z<m_iLength;z++) {
		m_pForward->m_pInput[2*z]   = inp[z];
		m_pForward->m_pInput[2*z+1] = 0.0;
	}

	for (z=m_iLength;z<m_iInternalLength;z++) {
		m_pForward->m_pInput[2*z]   = inp[m_iLength-1];
		m_pForward->m_pInput[2*z+1] = 0.0;
	}

	// Symmetrize by mirroring
	for (z=0;z<m_iLength;z++) {
		m_pForward->m_pInput[2*(z+m_iInternalLength)] = inp[m_iLength-z-1];
		m_pForward->m_pInput[2*(z+m_iInternalLength)+1] = 0.0;
	}

	for (z=m_iLength;z<m_iInternalLength;z++) {
		m_pForward->m_pInput[2*(z+m_iInternalLength)]   = inp[0];
		m_pForward->m_pInput[2*(z+m_iInternalLength)+1] = 0.0;
	}

	m_pForward->DoFFT();

	for (z=0;z<2*m_iInternalLength;z++) {
		if (z < 2*m_iInternalLength*(m_fWaveNumber / (33356.41 / g_fTimestepLength / g_iStride / 2.0))) {
			m_pInverse->m_pInput[2*z]   = m_pForward->m_pOutput[2*z];
			m_pInverse->m_pInput[2*z+1] = m_pForward->m_pOutput[2*z+1];
		} else {
			m_pInverse->m_pInput[2*z]   = 0.0;
			m_pInverse->m_pInput[2*z+1] = 0.0;
		}
	}

	m_pInverse->m_pInput[0] /= 2.0*m_iInternalLength;
	m_pInverse->m_pInput[1] /= 2.0*m_iInternalLength;
	for (z=1;z<m_iInternalLength;z++) {
		m_pInverse->m_pInput[2*z]   /= (double)m_iInternalLength;
		m_pInverse->m_pInput[2*z+1] /= (double)m_iInternalLength;
	}

	m_pInverse->DoFFT();

	outp.resize(m_iLength);

	// Borders: Linearly fade to original data, as smoothing might fail on borders
	for (z=0;z<MIN(20,m_iLength);z++)
		outp[z] = (double)z/19.0*m_pInverse->m_pOutput[2*z] + (1.0-(double)z/19.0)*inp[z];
	for (z=MIN(20,m_iLength);z<m_iLength-20;z++)
		outp[z] = m_pInverse->m_pOutput[2*z];
	for (z=MAX(1,m_iLength-20);z<m_iLength;z++)
		outp[z] = (double)(m_iLength-z-1)/19.0*m_pInverse->m_pOutput[2*z] + (1.0-(double)(m_iLength-z-1)/19.0)*inp[z];
}


void CSmoothener::Smoothen(const std::vector<CxDVector3> &inp, std::vector<CxDVector3> &outp) {

	int z, z2;
	std::vector<double> tv;


	tv.resize(m_iLength);
	for (z=0;z<3;z++) {
		for (z2=0;z2<m_iLength;z2++)
			tv[z2] = inp[z2][z];
		Smoothen(tv);
		for (z2=0;z2<m_iLength;z2++)
			outp[z2][z] = tv[z2];
	}
}


void CSmoothener::Smoothen(const std::vector<CxDMatrix3> &inp, std::vector<CxDMatrix3> &outp) {

	int z, z2;
	std::vector<double> tv;


	tv.resize(m_iLength);
	for (z=0;z<9;z++) {
		for (z2=0;z2<m_iLength;z2++)
			tv[z2] = inp[z2][z];
		Smoothen(tv);
		for (z2=0;z2<m_iLength;z2++)
			outp[z2][z] = tv[z2];
	}
}

	
CROAEngine::CROAEngine() {

	int z;


	m_iaTrajKinds.resize(7);
	for (z=0;z<7;z++)
		m_iaTrajKinds[z] = -1;
	m_pPDESolution = NULL;
	m_bIR = false;
	m_bSFG = false;
	m_bRaman = false;
	m_bVCD = false;
	m_bROA = false;
	m_bMagMom = false;
	m_bPola = false;
	m_fPDEInfo = NULL;
	m_bPDEFastMode = false;
	m_bPattern3 = false;
}


CROAEngine::~CROAEngine() {

	int z;


	for (z=0;z<(int)m_oaTrajectories.size();z++)
		delete m_oaTrajectories[z];

	for (z=0;z<(int)m_oaObservations.size();z++)
		delete m_oaObservations[z];

	if (m_pPDESolution != NULL) {
		delete[] m_pPDESolution;
		m_pPDESolution = NULL;
	}
}


CROAObservation::CROAObservation() {

	m_bMolPattern = false;
}


CROAObservation::CROAObservation(const CROAObservation &t) {

	int z;


	m_iType = t.m_iType;
	m_iDepth = t.m_iDepth;
	m_iWindowFunction = t.m_iWindowFunction;
	m_iWindowFunctionParameter = t.m_iWindowFunctionParameter;
	m_iZeroPadding = t.m_iZeroPadding;
	m_fSpecResolution = t.m_fSpecResolution;
	m_iSpecLength = t.m_iSpecLength;
	m_bFiniteDifferenceCorrection = t.m_bFiniteDifferenceCorrection;
	m_bSaveACF = t.m_bSaveACF;
	m_iQuantumCorrection = t.m_iQuantumCorrection;
	m_fQuantumCorrectionTemperature = t.m_fQuantumCorrectionTemperature;
	m_bCorrectFrequency = t.m_bCorrectFrequency;
	m_sName = t.m_sName;
	m_sMolName = t.m_sMolName;
	m_sTypeName = t.m_sTypeName;
	m_bCrossCorrelation = t.m_bCrossCorrelation;
	m_fLaserFreq = t.m_fLaserFreq;
	m_fTemperature = t.m_fTemperature;
	m_bCorrectTemperature = t.m_bCorrectTemperature;
	m_fCorrectTemperature = t.m_fCorrectTemperature;
	m_bUseCommutator = t.m_bUseCommutator;
	m_bSingleACFs = t.m_bSingleACFs;

	m_iaMolSelection.resize(t.m_iaMolSelection.size());
	for (z=0;z<(int)t.m_iaMolSelection.size();z++)
		m_iaMolSelection[z].assign(t.m_iaMolSelection[z].begin(),t.m_iaMolSelection[z].end());
}


CROAObservation::~CROAObservation() {
}






void CROAEngine::ParseObservations() {

	int iobs, i, z, z2, ti;
	CROAObservation *obs, *obs2;
	CxString buf;
	char buf2[64];
	char *p, *q;
	bool ir, raman, vcd, sfg, roa;
	CxIntArray tia;


	iobs = 0;
	while (true) {

		mprintf(WHITE,"\n    >>> Spectroscopy Observation %d >>>\n\n",iobs+1);

		obs = new CROAObservation();

		if (g_bRamanFromPolarizability) {
			mprintf("    Computing Raman spectrum from external polarizability data.\n\n");
			ir = false;
			raman = true;
			vcd = false;
			sfg = false;
			roa = false;
			m_bRaman = true;
			obs->m_iaMolSelection.resize(1);
			obs->m_iaMolSelection[0].push_back(0);
			goto _ramanfrompola;
		}

		mprintf("    The following types of spectra can be computed:\n");
		mprintf("      IR    - Infrared Spectrum\n");
		mprintf("      Raman - Raman Spectrum\n");
		mprintf("      VCD   - Vibrational Circular Dichroism Spectrum\n");
		mprintf("      ROA   - Raman Optical Activity Spectrum\n");
		mprintf("\n");

_again:
		ir = false;
		raman = false;
		vcd = false;
		sfg = false;
		roa = false;
		AskString_ND("    Spectra to compute for this observation (comma separated): ",&buf);

		p = buf.GetWritePointer();
		while (true) {
			while ((*p == ' ') || (*p == ','))
				p++;
			if (p == 0)
				break;
			q = p;
			while ((*q != ' ') && (*q != ',') && (*q != 0))
				q++;
			memcpy(buf2,p,q-p);
			buf2[q-p] = 0;
			if (mystricmp(buf2,"IR") == 0)
				ir = true;
			else if (mystricmp(buf2,"Raman") == 0)
				raman = true;
			else if (mystricmp(buf2,"VCD") == 0)
				vcd = true;
			else if (mystricmp(buf2,"ROA") == 0)
				roa = true;
			else {
				eprintf("Error: Unknown input: \"%s\".\n\n",buf2);
				goto _again;
			}
			if (*q == 0)
				break;
			p = q+1;
		}

		if (!ir && !raman && !vcd && !sfg && !roa) {
			eprintf("Error: At least one type of spectrum needs to be selected.\n\n");
			goto _again;
		}
		mprintf("\n");

		if (ir)
			m_bIR = true;
		if (raman)
			m_bRaman = true;
		if (vcd)
			m_bVCD = true;
		if (roa)
			m_bROA = true;

		{
			obs->m_iaMolSelection.resize(g_oaMolecules.GetSize());
			for (z=0;z<g_oaMolecules.GetSize();z++) {
				if (AskYesNo("    Observe molecules of type %s for this spectrum (y/n)? [yes] ",
					true,((CMolecule*)g_oaMolecules[z])->m_sName)) {
					mprintf("\n    Enter \"#\" to supply a string pattern for active molecules.\n\n");
_again2:
					AskString("      Which molecules of type %s to observe (e.g. 1,5-7)? [all] ",
						&buf,"*",((CMolecule*)g_oaMolecules[z])->m_sName);
					if (mystricmp((const char*)buf,"*") == 0) {
						for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
							obs->m_iaMolSelection[z].push_back(z2);
					} else if (mystricmp((const char*)buf,"#") == 0) {
						mprintf("\n");
						mprintf("      Please enter a string which consists of %d characters \"0\" or \"1\", where\n",((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize());
						mprintf("      \"0\" means to leave this molecule out, and \"1\" means to include it.\n\n");
						AskString_ND("      Enter pattern: ",&buf);
						if (buf.GetLength() != ((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize()) {
							eprintf("Error: String does not have a length of %d.\n",((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize());
							goto _again2;
						}
						for (z2=0;z2<buf.GetLength();z2++) {
							if ((buf[z2] != '0') && (buf[z2] != '1')) {
								eprintf("Error: Found invalid character \"%c\" - only \"0\" and \"1\" are allowed.\n",buf[z2]);
								goto _again2;
							}
						}
						for (z2=0;z2<buf.GetLength();z2++)
							if (buf[z2] == '1')
								obs->m_iaMolSelection[z].push_back(z2);
						obs->m_bMolPattern = true;
					} else {
						tia.RemoveAll();
						if (!ParseIntList((const char*)buf,&tia))
							goto _again2;
						for (z2=0;z2<tia.GetSize();z2++) {
							if ((tia[z2] < 1) || (tia[z2] > ((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize())) {
								eprintf("Error: Molecule index out of range (found %d, allowed 1..%d).\n",tia[z2],((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize());
								goto _again2;
							}
						}
						for (z2=0;z2<tia.GetSize();z2++)
							obs->m_iaMolSelection[z].push_back(tia[z2]-1);
					}
					mprintf("      Will include %lu/%d molecules of this type.\n",(unsigned long)obs->m_iaMolSelection[z].size(),((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize());
				}
			}
			mprintf("\n");
		}

_ramanfrompola:
		
_resagain:
		if (g_iTrajSteps != -1) {
			obs->m_iDepth = (int)(0.5 * g_iTrajSteps/g_iStride);
            if ((int)floor(0.5 + (g_fTimestepLength*g_iStride) / 0.5) >= 1) // MB 181227: Avoid division by zero
                ti = 4096 / (int)floor(0.5 + (g_fTimestepLength*g_iStride) / 0.5);
            else
                ti = 4096;
			//ti = 4096 / (int)floor(0.5 + (g_fTimestepLength*g_iStride) / 0.5);
			if (obs->m_iDepth > ti)
				obs->m_iDepth = ti;
			if (obs->m_iDepth < 1)
				obs->m_iDepth = 1;
			obs->m_iDepth = AskUnsignedInteger("    Enter the correlation depth of the correlation functions (in trajectory frames): [%d] ",
				obs->m_iDepth, obs->m_iDepth);
		} else
			obs->m_iDepth = AskUnsignedInteger("    Enter the correlation depth of the correlation functions (in trajectory frames): [512] ",
				512);

		if (obs->m_iDepth < 1) {
			eprintf("Error: Invalid input.\n");
			goto _resagain;
		}

		if ((g_iTrajSteps != -1) && (obs->m_iDepth >= (int)(0.6 * g_iTrajSteps/g_iStride))) {
			mprintf("\n");
			mprintf(RED,"    Warning: ");
			mprintf("The resolution is larger than half of the processed trajectory step count.\n");
			mprintf("             This leads to noisy and inaccurate spectra.\n\n");
			if (!AskYesNo("    Use this resolution anyway (y/n)? [no] ",false)) {
				mprintf("\n");
				goto _resagain;
			}
		}

		i = CalcFFTSize(obs->m_iDepth, false);
		if (obs->m_iDepth != i) {
			mprintf(WHITE,"    The next \"fast\" size for FFT is %d. Using this value instead of %d as depth.\n",
				i, obs->m_iDepth);
			obs->m_iDepth = i;
		}
		
		if(g_bAdvanced2) {
			obs->m_iWindowFunction = AskRangeInteger("    Window function: cos^2(a*t) (1), exp(-t/a) (2), exp(-(t/a)^2) (3) [1] ",
				1, 3, 1);
			if (obs->m_iWindowFunction == 1) {
				mprintf("\n    The parameter \"a\" is chosen according to the correlation depth.\n\n");
				obs->m_iWindowFunctionParameter = 0;
			} else if (obs->m_iWindowFunction == 2) {
				obs->m_iWindowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in trajectory frames): [%d] ",
					obs->m_iDepth / 4, obs->m_iDepth / 4);
			} else if (obs->m_iWindowFunction == 3) {
				obs->m_iWindowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in trajectory frames): [%d] ",
					obs->m_iDepth / 2, obs->m_iDepth / 2);
			} else {
				eprintf("Unknown input.\n");
				abort();
			}
		} else {
			obs->m_iWindowFunction = 1;
			obs->m_iWindowFunctionParameter = 0;
			mprintf("    Using Hann window function (cos^2).\n");
		}
		
		if (g_bAdvanced2) {

			obs->m_iZeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ",
				obs->m_iDepth * 3, obs->m_iDepth * 3);
			i = CalcFFTSize(obs->m_iDepth + obs->m_iZeroPadding, false);
			if (obs->m_iDepth + obs->m_iZeroPadding != i) {
				mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n",
					i, i-obs->m_iDepth);
				obs->m_iZeroPadding = i-obs->m_iDepth;
			}

		} else {

			obs->m_iZeroPadding = obs->m_iDepth * 3;
			i = CalcFFTSize(obs->m_iDepth + obs->m_iZeroPadding, false);
			if (obs->m_iDepth + obs->m_iZeroPadding != i) {
				mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n",
					i, i-obs->m_iDepth);
				obs->m_iZeroPadding = i-obs->m_iDepth;
			}
			mprintf("    Inserting %d zeros for zero padding.\n",obs->m_iZeroPadding);
		}
		
		double possibleRange = 33356.41 / g_fTimestepLength / g_iStride / 2.0;

		obs->m_fSpecResolution = possibleRange / (
			(obs->m_iDepth + obs->m_iZeroPadding));

		mprintf("    This results in a spectral resolution of %.2f cm^-1 per data point.\n",
			obs->m_fSpecResolution);
		mprintf("\n    A time step length of %.2f fs at stride %d allows a spectral range up to %.2f cm^-1.\n",
			g_fTimestepLength, g_iStride, possibleRange);
		double specLimit = AskRangeFloat("\n    Calculate spectrum up to which wave number (cm^-1)? [%.2f cm^-1] ",
			0, possibleRange, ((possibleRange*0.9) < 5000.0) ? (possibleRange*0.9) : 5000.0, ((possibleRange*0.9) < 5000.0) ? (possibleRange*0.9) : 5000.0 );
		obs->m_iSpecLength = (int)(specLimit / obs->m_fSpecResolution);
		mprintf("\n    This corresponds to %d data points.\n",obs->m_iSpecLength);
		mprintf("\n");
		
		if (g_bAdvanced2) {
			obs->m_bFiniteDifferenceCorrection = AskYesNo("    Apply finite difference correction (y/n)? [yes] ", true);
			mprintf("\n");
		} else {
			mprintf("    Automatically applying finite difference correction (sinus cardinalis).\n");
			mprintf("    Use the \"advanced mode\" to switch it off.\n\n");
			obs->m_bFiniteDifferenceCorrection = true;
		}
		
		if (g_bAdvanced2) {
			obs->m_bSaveACF = AskYesNo("    Save also correlation function(s) (y/n)? [no] ", false);
			if (sfg)
				obs->m_bSingleACFs = AskYesNo("    Save all SFG correlation component spectra separately (y/n)? [no] ", false);
			else
				obs->m_bSingleACFs = false;
			mprintf("\n");
		} else {
			obs->m_bSaveACF = false;
			obs->m_bSingleACFs = false;
		}

		if (g_bAdvanced2)
			obs->m_bUseCommutator = AskYesNo("    Use \"commutator trick\" for cross-correlations (y/n)? [yes] ",true);
		else {
			mprintf("\n    Using the \"commutator trick\" for enhanced sampling due to time reversibility.\n");
			mprintf("    Use the \"advanced mode\" to switch it off.\n\n");
			obs->m_bUseCommutator = true;
		}
		

		obs->m_bCorrectFrequency = AskYesNo("    Correct frequency shift of the Verlet integrator (y/n)? [no] ",false);
		if (obs->m_bCorrectFrequency) {
			ParseCorrectWavenumber();
			mprintf("\n    Due to the frequency correction, the spectral range is now up to %.2f cm^-1.\n\n",
				CorrectWavenumber(specLimit) );
		}

		{
			if (g_bAdvanced2 && !g_bRamanFromPolarizability)
				obs->m_bCrossCorrelation = AskYesNo("    Include intermolecular cross-correlations (y/n)? [no] ",false);
			else
				obs->m_bCrossCorrelation = false;

		}

		obs->m_sMolName = "";
		for (z=0;z<(int)obs->m_iaMolSelection.size();z++) {
			if (obs->m_iaMolSelection[z].size() == 0)
				continue;
			obs->m_sMolName.strcat("_");
			obs->m_sMolName.strcat(((CMolecule*)g_oaMolecules[z])->m_sName);
			if ((int)obs->m_iaMolSelection[z].size() != ((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize()) {
				if (!obs->m_bMolPattern) {
					for (z2=0;z2<(int)obs->m_iaMolSelection[z].size();z2++) {
						buf.sprintf("_%d",obs->m_iaMolSelection[z][z2]+1);
						obs->m_sMolName.strcat((const char*)buf);
					}
				} else
					obs->m_sMolName.strcat("_pattern");
			}
		}

		if (raman || roa) {
			obs->m_fLaserFreq = AskFloat("    Calculate scattering cross section for which laser wave number (cm^-1)? [20000.0] ",
				20000.0);
			obs->m_fTemperature = AskFloat("    Calculate scattering cross section for which temperature (K)? [350.0] ", 350.0);
		}

		obs->m_bCorrectTemperature = AskYesNo("    Correct spectrum for a certain simulation temperature (y/n)? [yes] ",true);

		if (obs->m_bCorrectTemperature) {
			if (raman || roa)
				obs->m_fCorrectTemperature = AskFloat("    Enter the simulation temperature (K): [%.1f] ",
					obs->m_fTemperature,obs->m_fTemperature);
			else
				obs->m_fCorrectTemperature = AskFloat("    Enter the simulation temperature (K): [350.0] ",350.0);
		}

		mprintf("\n    Will compute the following spectra:\n");
		if (ir) {
			mprintf("      - IR\n");
			obs2 = new CROAObservation(*obs);
			obs2->m_iType = ROA_SPECTRUM_IR;
			obs2->m_sTypeName = "ir";
			obs2->m_sName.Format("%s%s",(const char*)obs2->m_sTypeName,(const char*)obs2->m_sMolName);
			if (obs2->m_bCrossCorrelation)
				obs2->m_sName.strcat("_cc");
			m_oaObservations.push_back(obs2);
		}
		if (raman) {
			mprintf("      - Raman\n");
			obs2 = new CROAObservation(*obs);
			obs2->m_iType = ROA_SPECTRUM_RAMAN;
			obs2->m_sTypeName = "raman";
			obs2->m_sName.Format("%s%s",(const char*)obs2->m_sTypeName,(const char*)obs2->m_sMolName);
			if (obs2->m_bCrossCorrelation)
				obs2->m_sName.strcat("_cc");
			m_oaObservations.push_back(obs2);
		}
		if (vcd) {
			mprintf("      - VCD\n");
			obs2 = new CROAObservation(*obs);
			obs2->m_iType = ROA_SPECTRUM_VCD;
			obs2->m_sTypeName = "vcd";
			obs2->m_sName.Format("%s%s",(const char*)obs2->m_sTypeName,(const char*)obs2->m_sMolName);
			if (obs2->m_bCrossCorrelation)
				obs2->m_sName.strcat("_cc");
			m_oaObservations.push_back(obs2);
		}
		if (roa) {
			mprintf("      - ROA\n");
			obs2 = new CROAObservation(*obs);
			obs2->m_iType = ROA_SPECTRUM_ROA;
			obs2->m_sTypeName = "roa";
			obs2->m_sName.Format("%s%s",(const char*)obs2->m_sTypeName,(const char*)obs2->m_sMolName);
			if (obs2->m_bCrossCorrelation)
				obs2->m_sName.strcat("_cc");
			m_oaObservations.push_back(obs2);
		}

		delete obs;

		mprintf(WHITE,"\n    <<< Spectroscopy Observation %d <<<\n\n",iobs+1);


		iobs++;
		if (!AskYesNo("    Add another observation (y/n)? [no] ",false))
			break;
	}

	mprintf("\n    Added %lu observations.\n\n",(unsigned long)m_oaObservations.size());
}


CROATrajectory::CROATrajectory() {

	m_iTSPos = 0;
	m_iHistory = 0;
	m_pFile = NULL;
	m_bBQB = false;
	m_bVoronoi = false;
	m_pBQB = NULL;
}


CROATrajectory::~CROATrajectory() {

	int z;

	for (z=0;z<(int)m_oaTSHistory.size();z++)
		if (m_oaTSHistory[z] != NULL)
			delete m_oaTSHistory[z];

	for (z=0;z<(int)m_oaMolecules.size();z++)
		delete m_oaMolecules[z];

	if (m_pFile != NULL) {
		fclose(m_pFile);
		m_pFile = NULL;
	}

	if (m_pBQB != NULL) {
		delete m_pBQB;
		m_pBQB = NULL;
	}
}


void CROATrajectory::SetHistory(int i) {

	int z;

	m_iHistory = i;
	for (z=0;z<(int)m_oaTSHistory.size();z++)
		if (m_oaTSHistory[z] != NULL)
			delete m_oaTSHistory[z];
	m_oaTSHistory.clear();
	for (z=0;z<i;z++)
		m_oaTSHistory.push_back(new CTimeStep);
}


bool CROATrajectory::OpenFile() {

	if (!FileExist(m_sFileName.c_str())) {
		eprintf("\n    Error: File does not exist or could not be read: \"%s\".\n\n",
			m_sFileName.c_str());
		return false;
	}

	if ((mystricmp(GetFileExtension(m_sFileName.c_str()),"bqb") == 0) || 
		(mystricmp(GetFileExtension(m_sFileName.c_str()),"bbq") == 0) || 
		(mystricmp(GetFileExtension(m_sFileName.c_str()),"emp") == 0) || 
		(mystricmp(GetFileExtension(m_sFileName.c_str()),"blist") == 0)) {

		if (m_pBQB != NULL)
			delete m_pBQB;
		m_pBQB = new CBQBFile(*g_pBQBInterface);
		m_bBQB = true;

		if (!m_pBQB->OpenRead(m_sFileName)) {
			eprintf("\n    Error: Not a valid BQB file.\n\n");
			return false;
		}

		mprintf("\n      Is a valid BQB file, ");
		if (m_pBQB->GetTotalFrameCount() != -1)
			mprintf("%d frames (with index)",m_pBQB->GetTotalFrameCount());
		else
			mprintf("no index");
		mprintf(".\n\n");

	} else if (mystricmp(GetFileExtension(m_sFileName.c_str()),"cube") == 0) {

		m_pFile = fopen(m_sFileName.c_str(),"rb");

	} else if (mystricmp(GetFileExtension(m_sFileName.c_str()),"voronoi") == 0) {

		m_pFile = fopen(m_sFileName.c_str(),"rb");
		m_bVoronoi = true;

	} else {

		eprintf("\n    Error: Not a valid file extension: \"%s\".\n\n",
			GetFileExtension(m_sFileName.c_str()));
		return false;

	}

	return true;
}


bool CROATrajectory::ReadTimeStep() {

	if (m_iHistory == 0) {
		eprintf("CROATrajectory::ReadTimeStep(): Error: m_iHistory == 0.\n");
		abort();
	}

	m_iTSPos++;
	if (m_iTSPos >= m_iHistory)
		m_iTSPos = 0;

	if (m_bBQB) {
		if (!m_oaTSHistory[m_iTSPos]->ReadBQB(m_pBQB,true))
			return false;
	} else if (m_bVoronoi) {
		if (!m_oaTSHistory[m_iTSPos]->ReadVoronoi(m_pFile,true))
			return false;
	} else {
		if (!m_oaTSHistory[m_iTSPos]->ReadCube(m_pFile,true))
			return false;
	}

	m_oaTSHistory[m_iTSPos]->UniteMolecules(false);
	m_oaTSHistory[m_iTSPos]->CalcCenters();

	return true;
}


bool CROATrajectory::SkipTimeStep() {

	CTimeStep ts;

	if (m_bBQB) {
		if (!ts.SkipBQB(m_pBQB))
			return false;
	} else if (m_bVoronoi) {
		if (!ts.SkipVoronoi(m_pFile))
			return false;
	} else {
		if (!ts.SkipCube(m_pFile))
			return false;
	}

	return true;
}


void CROATrajectory::Rewind() {

	if (m_bBQB) {
		if (m_pBQB == NULL) {
			eprintf("CROATrajectory::Rewind(): Error: BQB trajectory not initialized.\n");
			abort();
		}
		if (!m_pBQB->Rewind()) {
			eprintf("CROATrajectory::Rewind(): Error rewinding BQB trajectory.\n");
			abort();
		}
	} else if (m_bVoronoi) {
		if (m_pFile != NULL)
			fclose(m_pFile);
		m_pFile = fopen(m_sFileName.c_str(),"rb");
		if (m_pFile == NULL) {
			eprintf("CROATrajectory::Rewind(): Error rewinding Voronoi trajectory.\n");
			abort();
		}
	} else {
		if (m_pFile != NULL)
			fclose(m_pFile);
		m_pFile = fopen(m_sFileName.c_str(),"rb");
		if (m_pFile == NULL) {
			eprintf("CROATrajectory::Rewind(): Error rewinding CUBE trajectory.\n");
			abort();
		}
	}
}


void CROATrajectory::CalcVelocities() {

	int z, z2;
	CMolecule *m;
	CSingleMolecule *sm;
	CVirtualAtom *vi;
	unsigned long index;
	CxDVector3 dist;

	this->GetTimeStep(0)->m_vaVelocities.SetSize(g_iGesVirtAtomCount);
	for (z=0;z<g_iGesAtomCount;z++)
	{
		dist = this->GetTimeStep(1)->m_vaCoords[z] - this->GetTimeStep(-1)->m_vaCoords[z];
		if(g_bPeriodicX) {
			while(dist[0] > g_fBoxX / 2.0) dist[0] -= g_fBoxX;
			while(dist[0] < -g_fBoxX / 2.0) dist[0] += g_fBoxX;
		}
		if(g_bPeriodicY) {
			while(dist[1] > g_fBoxY / 2.0) dist[1] -= g_fBoxY;
			while(dist[1] < -g_fBoxY / 2.0) dist[1] += g_fBoxY;
		}
		if(g_bPeriodicZ) {
			while(dist[2] > g_fBoxZ / 2.0) dist[2] -= g_fBoxZ;
			while(dist[2] < -g_fBoxZ / 2.0) dist[2] += g_fBoxZ;
		}
		this->GetTimeStep(0)->m_vaVelocities[z] = dist * 1000.0 / 2.0 / g_fTimestepLength;
	}
	for (z = 0; z < g_oaVirtualAtoms.GetSize(); z++) {
		vi = (CVirtualAtom*)g_oaVirtualAtoms[z];
		m = (CMolecule*)g_oaMolecules[vi->m_iMolecule];
		for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++) {
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
			index = ((CxIntArray*)sm->m_oaAtomOffset[sm->m_baAtomIndex.GetSize()-1])->GetAt(vi->m_iMolVirtAtom);
			if (((CVirtualAtom *)g_oaVirtualAtoms[z])->m_iMode == 2) {
				dist = CxDVector3(0.0, 0.0, 0.0);
			} else {
				dist = GetTimeStep(-1)->m_vaCoords[index] - GetTimeStep(1)->m_vaCoords[index];
				if(g_bPeriodicX) {
					while(dist[0] > g_fBoxX / 2.0) dist[0] -= g_fBoxX;
					while(dist[0] < -g_fBoxX / 2.0) dist[0] += g_fBoxX;
				}
				if(g_bPeriodicY) {
					while(dist[1] > g_fBoxY / 2.0) dist[1] -= g_fBoxY;
					while(dist[1] < -g_fBoxY / 2.0) dist[1] += g_fBoxY;
				}
				if(g_bPeriodicZ) {
					while(dist[2] > g_fBoxZ / 2.0) dist[2] -= g_fBoxZ;
					while(dist[2] < -g_fBoxZ / 2.0) dist[2] += g_fBoxZ;
				}
			}
			this->GetTimeStep(0)->m_vaVelocities[index] = dist * 1000.0 / 2.0 / g_fTimestepLength;
		}
	}
}


CTimeStep* CROATrajectory::GetTimeStep(int i) {

	int t;

	if (m_iHistory == 0) {
		eprintf("CROATrajectory::GetTimeStep(): Error: m_iHistory == 0.\n");
		abort();
	}

	if ((i < -1) || (i > 1)) {
		eprintf("CROATrajectory::GetTimeStep(): Error: %d.\n",i);
		abort();
	}

	t = (m_iTSPos+i-1) % m_iHistory;
	if (t < 0)
		t += m_iHistory;

	return m_oaTSHistory[t];
}


bool CROAEngine::CheckTrajColumns(CROATrajectory *tr) {

	CBQBTrajectoryFrame btf(*g_pBQBInterface);
	bool b;


	mprintf("    Checking data fields of trajectory...\n");

	if (tr->m_bVoronoi) {
		if (m_bMagMom) {
			eprintf("Error: .voronoi trajectory does not contain magnetic moments required for VCD / ROA spectra.\n");
			return false;
		}
		mprintf("      Ok (assuming that .voronoi trajectory is valid).\n");
		return true;
	}

	if (!tr->m_bBQB) {
		eprintf("Error: Trajectory is neither in EMP nor in Voronoi format.\n");
		return false;
	}

	if (!tr->m_pBQB->ReadFrame()) {
		eprintf("Error while reading first EMP frame.\n");
		return false;
	}

	if ((tr->m_pBQB->GetFrameType() != 4) || (tr->m_pBQB->GetFrameTypeVersion() != 0)) {
		eprintf("Error: Wrong BQB frame type (expected: 4v0, found: %dv%d).\n",
			tr->m_pBQB->GetFrameType(),tr->m_pBQB->GetFrameTypeVersion());
		eprintf("You need to supply the \"properties.emp\" file which was created in a gathering run before.\n");
		return false;
	}

	if (!btf.ReadFrame(tr->m_pBQB->GetFramePayload())) {
		eprintf("Error: Could not read trajectory frame from BQB payload.\n");
		return false;
	}

	b = false;
	mprintf("      Contains: ");
	if (btf.GetColumn("Label") != NULL) { if (b) mprintf(", "); mprintf("Labels"); b=true; }
	if (btf.GetColumn("PosX") != NULL) { if (b) mprintf(", "); mprintf("Positions"); b=true; }
	if (btf.GetColumn("VelX") != NULL) { if (b) mprintf(", "); mprintf("Velocities"); b=true; }
	if (btf.GetColumn("EChg") != NULL) { if (b) mprintf(", "); mprintf("El.Charges"); b=true; }
	if (btf.GetColumn("EDipX") != NULL) { if (b) mprintf(", "); mprintf("El.Dipoles"); b=true; }
	if (btf.GetColumn("EQXX") != NULL) { if (b) mprintf(", "); mprintf("El.Quadrupoles"); b=true; }
	if (btf.GetColumn("ECurX") != NULL) { if (b) mprintf(", "); mprintf("El.Currents"); b=true; }
	if (btf.GetColumn("MDipX") != NULL) { if (b) mprintf(", "); mprintf("Mag.Dipoles"); }
	mprintf("\n\n");


	if (btf.GetColumn("Label") == NULL) {
		eprintf("Error: Data file does not contain atom labels.\n");
		return false;
	}
	if (btf.GetColumn("PosX") == NULL) {
		eprintf("Error: Data file does not contain atom positions.\n");
		return false;
	}
	if (btf.GetColumn("EChg") == NULL) {
		eprintf("Error: Data file does not contain electric charges.\n");
		return false;
	}
	if (btf.GetColumn("EDipX") == NULL) {
		eprintf("Error: Data file does not contain electric dipoles.\n");
		return false;
	}

	if (m_bROA) {
		if (btf.GetColumn("EQXX") == NULL) {
			eprintf("Error: Data file does not contain electric quadrupoles required by ROA spectra.\n");
			return false;
		}
	}

	if (m_bMagMom) {
		if (btf.GetColumn("VelX") == NULL) {
			eprintf("Error: Data file does not contain velocities required by VCD and ROA spectra.\n");
			eprintf("       Switch on the magnetic dipole calculation in the gathering run to produce the data file.\n");
			return false;
		}
		if (btf.GetColumn("MDipX") == NULL) {
			eprintf("Error: Data file does not contain magnetic dipoles required by VCD and ROA spectra.\n");
			eprintf("       Switch on the magnetic dipole calculation in the gathering run to produce the data file.\n");
			return false;
		}
	}

	return true;
}


bool CROAEngine::Parse() {

	int i, z, z2, steps, ti;
	CROATrajectory *tr;
	CxString buf;
	CxDVec3Array va;
	int rty;
	CROAMolecule *rm;


	mprintf(WHITE,"\n    >>> New Spectroscopy Module >>>\n\n");

	mprintf("    This module computes vibrational spectra (IR/Raman/VCD/ROA) from trajectories of\n");
	mprintf("    volumetric electron density data. To use Wannier centers instead, use one of the\n");
	mprintf("    old spectroscopy modules in the main function menu (enter \"ir\", \"raman\", or \"vcd\").\n");
	mprintf("\n");

	if (g_bTegri) {
		eprintf("Error: Don't switch on Voronoi integration in the main function menu.\n");
		eprintf("       It will be switched on automatically.\n\n");
		return false;
	}

	mprintf("    There are two different operating modes available:\n\n");

	mprintf(YELLOW,"    (*) One-pass Mode:\n");
	mprintf("      Read all volumetric data trajectories, solve the current PDE, integrate the properties\n");
	mprintf("      and compute the desired spectra in one single pass. Easy to use, but takes a long time,\n");
	mprintf("      because only one single CPU core is utilized for the process.\n");

	mprintf("\n");

	mprintf(YELLOW,"    (*) Two-pass Mode:\n");
	mprintf("      Solve the current PDE for small parts of each trajectory individually in the first pass\n");
	mprintf("      (\"gathering run\"), allowing to parallelize this process on many CPU cores. In the second\n");
	mprintf("      pass (\"analyzing run\"), all data files (.emp) obtained in the first-pass runs are used to\n");
	mprintf("      yield the desired spectra.\n\n");

	if (AskYesNo("    Use two-pass mode (y/n)? [yes] ",true)) {
		if (AskYesNo("    Is this a gathering run (y) or the analyzing run (n)? [%s] ",
			!g_bElMagProperties,g_bElMagProperties?"no":"yes"))
			m_iMode = ROA_MODE_GATHER;
		else
			m_iMode = ROA_MODE_ANALYZE;
	} else
		m_iMode = ROA_MODE_COMBINED;

	if ((m_iMode == ROA_MODE_COMBINED) || (m_iMode == ROA_MODE_GATHER)) {

		if (!g_bVolumetricData) {

			mprintf("\n");
			eprintf("    Error: This requires volumetric electron density data in each frame of the trajectory.\n");
			eprintf("           The file formats .cube and .bqb are suitable for this.\n\n");
			return false;
		}
	}

	if (m_iMode == ROA_MODE_ANALYZE) {
		mprintf("\n");
		mprintf("    Make sure to adapt the time step length if you processed only every n-th frame while gathering.\n");
		mprintf("    If you, e.g., processed every 8-th frame from a trajectory with 0.5 fs time step, enter \"4.0\" below.\n");
	}

	m_bPattern3 = AskYesNo("    Are the frames in a non-uniform pattern in blocks of three (y/n)? [no] ",false);

	mprintf("\n");

	if (m_bPattern3)
		g_fTimestepLength = AskFloat("    Enter the physical time stride for the 3-frame blocks in fs: [4.0] ",4.0);
	else
		g_fTimestepLength = AskFloat("    Enter the physical time distance between successive trajectory frames in fs: [0.5] ",0.5);

	if (!g_bStrideParsed) {
		mprintf("\n");
		g_bStrideParsed = true;
		g_iBeginStep = AskUnsignedInteger("    In which trajectory frame to start processing the trajectory? [1] ",1) - 1;
		if (g_iBeginStep == -1)
			g_iBeginStep = 0;
		g_iMaxStep = AskUnsignedInteger("    How many trajectory frames to read (from this position on)? [all] ",0);
		if (g_iMaxStep == 0)
			g_iMaxStep = -1; // Alle verwenden
		if (g_iMaxStep > 0)
			g_iTrajSteps = g_iMaxStep;
		if (((m_iMode == ROA_MODE_GATHER) || (m_iMode == ROA_MODE_COMBINED)) && (!m_bPattern3)) {
			mprintf("\n    In most cases, it is sufficient to process frames only every 2 fs or even 4 fs for good spectra.\n\n");
	//		ti = (int)floor(4.01/g_fTimestepLength);
			ti = 1;
			g_iStride = AskUnsignedInteger("    Use every n-th read trajectory frame for the analysis: [%d] ",ti,ti);
		} else
			g_iStride = AskUnsignedInteger("    Use every n-th read trajectory frame for the analysis: [1] ",1);
	}


	if ((m_iMode == ROA_MODE_COMBINED) || (m_iMode == ROA_MODE_ANALYZE)) {


		ParseObservations();

		if (m_bVCD || m_bROA) {
			m_bMagMom = true;
			mprintf("    Magnetic dipole moments will be required.\n");
		}

		if (m_bRaman || m_bROA || m_bSFG) {
			m_bPola = true;
			mprintf("    Polarizabilities will be required.\n");
		}

		if (m_bMagMom || m_bPola)
			mprintf("\n");

		if (m_bPola) {

		/*	m_bAniso = AskYesNo("    Compute polarizability for all three spatial directions (y) or only for one (n)? [yes] ",true);

			if (!m_bAniso) {
				eprintf("\nNot yet implemented.\n\n");
				return false;
			}*/
			m_bAniso = true;

			m_bCentral = AskYesNo("    Use central finite differences (y) or one-sided differences (n) for polarizabilities? [no] ",false);

			m_fFieldStrength = AskFloat("    Enter electric field strength that was used for field trajectories (in a.u.): [5.0E-3] ",5.0E-3);
		}

		if (g_bAdvanced2) {

			m_bDumpMolecularProps = AskYesNo("    Write all molecular electronic/magnetic properties to file (y/n)? [no] ",false);

			m_bDumpMol1Props = AskYesNo("    Write electronic/magnetic properties of molecule 1 to file (y/n)? [no] ",false);

			m_bWriteACFs = AskYesNo("    Write autocorrelation functions of all quantities (y/n)? [no] ",false);

			m_bSmoothData = AskYesNo("    Smoothen properties before computing temporal derivatives (y/n)? [no] ",false);
			if (m_bSmoothData)
				m_fSmoothWaveNumber = AskFloat("    Smoothing: Remove all frequencies above which wave number (in cm^-1)? [5000] ",5000.0);

			m_bReplaceOutliers = AskYesNo("    Replace outliers in time series (y/n)? [no] ",false);

			if (m_bReplaceOutliers) {
				m_bReplaceOutliersElectric = AskYesNo("      Replace electric/polarizability outliers (y/n)? [yes] ",true);
				m_bReplaceOutliersMagnetic = AskYesNo("      Replace magnetic outliers (y/n)? [no] ",false);
				m_fReplaceOutliersThreshold = AskFloat("      Enter outlier threshold: [5.0] ",5.0);
			} else {
				m_bReplaceOutliersElectric = false;
				m_bReplaceOutliersMagnetic = false;
			}

			g_bQuadrupoleKeepTrace = AskYesNo("    Keep trace of quadrupole tensor (y/n)? [no] ",false);

			m_bReverseTraj = AskYesNo("    Reverse direction of input trajectory (y/n)? [no] ",false);

			g_bDipoleRefFixed = AskYesNo("    Use a common box-fixed reference point for electric / magnetic multipoles (y/n)? [no] ", false);

		} else {

			m_bDumpMolecularProps = false;
			m_bDumpMol1Props = false;
			m_bWriteACFs = false;
			m_bSmoothData = false;
			m_bReplaceOutliers = false;
			m_bReplaceOutliersElectric = false;
			m_bReplaceOutliersMagnetic = false;
			g_bQuadrupoleKeepTrace = false;
			m_bReverseTraj = false;
		}
	}

	if (m_iMode == ROA_MODE_GATHER) {

		mprintf("\n");
		m_bMagMom = AskYesNo("    Compute magnetic moments (time consuming, required for VCD and ROA) (y/n)? [yes] ",true);
		mprintf("\n");
	}

	parseCoreCharges();

	if ((m_iMode == ROA_MODE_COMBINED) || (m_iMode == ROA_MODE_GATHER)) {

		if (m_bMagMom) {

			if (g_fTimestepLength > 0.501) {
				mprintf("\n");
				eprintf("    Warning: ");
				mprintf(             "The finite differences approach for the magnetic moments only works accurately\n");
				mprintf("             for small time steps <= 0.5 fs. Your trajectory has a time step of %.2f fs.\n\n",g_fTimestepLength);
				if (!AskYesNo("    Continue anyway (y/n)? [no] ",false))
					return false;
			}

			if (g_iStride != 1) {

				mprintf("\n");
				eprintf("    Warning: ");
				mprintf(             "You have chosen to process only every %d-",g_iStride);
				if (g_iStride == 2)
					mprintf("nd");
				else if (g_iStride == 3)
					mprintf("rd");
				else
					mprintf("th");
				mprintf(" step of the input trajectory.\n");
				mprintf("             The finite differences approach for the magnetic moments requires a time step\n");
				mprintf("             not larger than 0.5 fs. Therefore, TRAVIS will use the following scheme:\n\n");

				mprintf("      Read time steps 1 - 3.\n");
				mprintf("      Compute magnetic dipole in step 2 by finite differences from step 1 and 3.\n");

				if (g_iStride > 4) {

					mprintf("      Skip time steps 4 - %d.\n",g_iStride);
					mprintf("      Read time steps %d - %d.\n",g_iStride+1,g_iStride+3);

				} else if (g_iStride == 4) {

					mprintf("      Skip time step 4.\n");
					mprintf("      Read time steps 5 - 7.\n");

				} else if (g_iStride == 3) {

					mprintf("      Read time steps 4 - 6.\n");

				} else if (g_iStride == 2) {

					mprintf("      Read time steps 4 and 5.\n");

				}

				mprintf("      Compute magnetic dipole in step %d by finite differences from step %d and %d.\n",
					g_iStride+2,g_iStride+1,g_iStride+3);
				mprintf("      ...\n\n");

				if (!AskYesNo("    Acknowledged (y/n): [yes] ",true))
					return false;
			}
		}

		ParseVoronoiRadii();

		if (m_bMagMom) {

			if (g_bAdvanced2) {

				g_fBackgroundDensity = AskFloat("    Background density to improve PDE solver convergence (in a.u.): [1e-6] ", 1.0e-6);
				g_iPDEMaxIter = AskInteger("    Maximum number of iterations in PDE solver: [25] ", 25);
				m_bPDERestart = AskYesNo("    Restart PDE solver from solution of last step (y/n)? [yes] ",true);
				m_bPDEFastMode = AskYesNo("    Use PDE fast mode (activate unless there are problems) (y/n)? [yes] ",true);
				if (AskYesNo("    Write verbose info on PDE solving to text file (y/n)? [no] ",false))
					m_fPDEInfo = OpenFileWrite("pde_verbose.txt",true);

				mprintf("\n");

			} else {

				g_fBackgroundDensity = 1.0e-6;
				g_iPDEMaxIter = 25;
				m_bPDERestart = true;
				m_bPDEFastMode = true;

				mprintf("    Using a background density of %.2e and a maximum iteration count of %d in PDE solver.\n",
					g_fBackgroundDensity,g_iPDEMaxIter);
				mprintf("    Restarting PDE solver from solution of last step.\n\n");
			}

			mprintf("    ROA spectra require a threshold of <= 0.001, but for VCD, normally 0.01 is enough (faster).\n\n");

			if ((m_iMode == ROA_MODE_COMBINED) && !m_bROA)
				g_fPDEConvThresh = AskFloat("    Relative convergence threshold for PDE solver: [0.01] ", 0.01);
			else
				g_fPDEConvThresh = AskFloat("    Relative convergence threshold for PDE solver: [0.001] ", 0.001);
		}
	}

	for (z=0;z<g_oaMolecules.GetSize();z++) {
		((CMolecule*)g_oaMolecules[z])->m_iDipoleMode = 7;
		if (m_bMagMom)
			((CMolecule*)g_oaMolecules[z])->m_iMagneticDipoleMode = 7;
		ParseAtom("#2", z, ((CMolecule *)g_oaMolecules[z])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[z])->m_iDipoleCenterIndex);
	}

	if (m_iMode == ROA_MODE_GATHER) {

		mprintf("\n");
		m_bGatherWriteCSV = AskYesNo("    Write CSV text file with atomic properties in addition to EMP file (y/n)? [no] ",false);

		if (m_bGatherWriteCSV) {
			mprintf("\n");
			mprintf(WHITE,"    Here is a list of the columns in the \"properties.csv\" file:\n\n");

			mprintf("      \"Step\"        - Number of the time step.\n");
			mprintf("      \"AtomID\"      - Number of the atom in the system (in input trajectory order).\n");
			mprintf("      \"Label\"       - Element label of the atom.\n");
			mprintf("\n");

			mprintf("      \"PosX\"        - X position of the atom core (in pm).\n");
			mprintf("      \"PosY\"        - Y position of the atom core (in pm).\n");
			mprintf("      \"PosZ\"        - Z position of the atom core (in pm).\n");
			mprintf("\n");

			mprintf("      \"Vol\"         - Volume of the atom's Voronoi cell (in pm^3).\n");
			mprintf("\n");

			if (m_bMagMom) {
				mprintf("      \"VelX\"        - X velocity of the atom core (in pm/ps).\n");
				mprintf("      \"VelY\"        - Y velocity of the atom core (in pm/ps).\n");
				mprintf("      \"VelZ\"        - Z velocity of the atom core (in pm/ps).\n");
				mprintf("\n");
			}

			mprintf("      \"Charge\"      - Total electric charge of the atom (core + electrons) (in e).\n");
			mprintf("      \"Core Charge\" - Electric charge of the atom core (in e).\n");
			mprintf("\n");

			mprintf("      \"ElDipX\"      - X component of the atom's electric dipole vector (in e*pm).\n");
			mprintf("      \"ElDipY\"      - Y component of the atom's electric dipole vector (in e*pm).\n");
			mprintf("      \"ElDipZ\"      - Z component of the atom's electric dipole vector (in e*pm).\n");
			mprintf("\n");

			mprintf("      \"ElChZX\"      - X coordinate of the electronic charge center relative to the atom core position (in pm).\n");
			mprintf("      \"ElChZY\"      - Y coordinate of the electronic charge center relative to the atom core position (in pm).\n");
			mprintf("      \"ElChZZ\"      - Z coordinate of the electronic charge center relative to the atom core position (in pm).\n");
			mprintf("\n");

			mprintf("      \"ElQuadXX\"    - XX component of the atom's electric quadrupole tensor (in e*nm^2).\n");
			mprintf("      \"ElQuadXY\"    - XY component of the atom's electric quadrupole tensor (in e*nm^2).\n");
			mprintf("      \"ElQuadXZ\"    - XZ component of the atom's electric quadrupole tensor (in e*nm^2).\n");
			mprintf("      \"ElQuadYX\"    - YX component of the atom's electric quadrupole tensor (in e*nm^2).\n");
			mprintf("      \"ElQuadYY\"    - YY component of the atom's electric quadrupole tensor (in e*nm^2).\n");
			mprintf("      \"ElQuadYZ\"    - YZ component of the atom's electric quadrupole tensor (in e*nm^2).\n");
			mprintf("      \"ElQuadZX\"    - ZX component of the atom's electric quadrupole tensor (in e*nm^2).\n");
			mprintf("      \"ElQuadZY\"    - ZY component of the atom's electric quadrupole tensor (in e*nm^2).\n");
			mprintf("      \"ElQuadZZ\"    - ZZ component of the atom's electric quadrupole tensor (in e*nm^2).\n");
			mprintf("\n");

			if (m_bMagMom) {
				mprintf("      \"ElCurrX\"     - X component of the atom's electric current (in Bohr magnetons/pm).\n");
				mprintf("      \"ElCurrY\"     - Y component of the atom's electric current (in Bohr magnetons/pm).\n");
				mprintf("      \"ElCurrZ\"     - Z component of the atom's electric current (in Bohr magnetons/pm).\n");
				mprintf("\n");

				mprintf("      \"MagDipX\"     - X component of the atom's magnetic dipole vector (in Bohr magnetons).\n");
				mprintf("      \"MagDipY\"     - Y component of the atom's magnetic dipole vector (in Bohr magnetons).\n");
				mprintf("      \"MagDipZ\"     - Z component of the atom's magnetic dipole vector (in Bohr magnetons).\n");
				mprintf("\n");
			}

			mprintf("      All quantities which have a reference point (dipole, quadrupole, etc.)\n");
			mprintf("      are printed relative to the atom core position.\n");

			mprintf("\n\n");

			(void)AskYesNo("    Acknowledged: [yes] ",true);
			mprintf("\n");
		}

		if (g_bAdvanced2)
			m_bGatherCheck = AskYesNo("    Extract EMP file to CSV file (for checking) at the end of the run (y/n)? [no] ",false);
		else
			m_bGatherCheck = false;

		mprintf("\n    Opening input trajectory...\n");
		tr = new CROATrajectory;
		if (m_bMagMom)
			tr->SetHistory(3);
		else
			tr->SetHistory(1);
		m_iaTrajKinds[0] = (int)m_iaInvTrajKinds.size();
		m_iaInvTrajKinds.push_back(0);
		m_oaTrajectories.push_back(tr);
		tr->m_sFileName = (const char*)g_sInputTraj;
		if (!tr->OpenFile()) {
			eprintf("Error: Could not re-open input trajectory \"%s\".\n",(const char*)g_sInputTraj);
			abort();
		}

/*		if (g_pBQBEngine != NULL) {
			if (g_pBQBEngine->m_pExtrapolator != NULL)
				g_pBQBEngine->m_pExtrapolator->Clear();
			if (g_pBQBEngine->m_pExtrapolatorCorr != NULL)
				g_pBQBEngine->m_pExtrapolatorCorr->Clear();
			if (g_pBQBEngine->m_pExtrapolatorXYZ != NULL)
				g_pBQBEngine->m_pExtrapolatorXYZ->Clear();
		}*/

		if (tr->m_bBQB) {
			steps = tr->m_pBQB->GetTotalFrameCount();
			mprintf("    Trajectoriy contains %d steps.\n",steps);
		} else
			steps = -1;

		mprintf("\n    Reading first time step from trajectory...\n");

		for (z=0;z<(int)m_oaTrajectories.size();z++) {
			mprintf("      %s ...\n",m_oaTrajectories[z]->m_sFileName.c_str());
			if (!m_oaTrajectories[z]->ReadTimeStep()) {
				eprintf("Error: Could not read time step.\n");
				abort();
			}
		}

		mprintf("    Setting up Voronoi integration...\n");
		if (!PrepareVori(m_oaTrajectories[0]->GetTimeStep(1),true,true,true)) {
			eprintf("Error while setting up Voronoi integration.\n");
			return false;
		}
		
		mprintf("\n    Computing dipole moment of first molecule from first time step of trajectory...\n");

		if (m_bMagMom) {
			PrepareCurrentDensity( m_oaTrajectories[0]->GetTimeStep(1),  m_oaTrajectories[0]->GetTimeStep(1), true  );
			PrepareCurrentDensity( m_oaTrajectories[0]->GetTimeStep(0),  m_oaTrajectories[0]->GetTimeStep(1), false );
			PrepareCurrentDensity( m_oaTrajectories[0]->GetTimeStep(-1), m_oaTrajectories[0]->GetTimeStep(1), false );
		}

		g_bVoroIntegrateCharge = true;
		g_bVoroIntegrateDipoleMoment = true;
		g_bVoroIntegrateQuadrupoleMoment = true;
		g_bVoroIntegrateTotalCurrent = false;
		g_bVoroIntegrateMagneticMoment = false;
		g_bCubeTimeDev = false;

		g_pTetraPak->ProcessStep(m_oaTrajectories[0]->GetTimeStep(1),true);
		m_oaTrajectories[0]->GetTimeStep(1)->CalcDipoles(false);

		mprintf("\n");

		if (m_bMagMom && !m_bPattern3) {
			mprintf(WHITE,"    Important note:\n");
			mprintf("    If you perform multiple gathering runs on different intervals of a volumetric trajectory file,\n");
			mprintf("    make sure that the parts have an overlap of exactly 2 frames, because you will lose 2 frames\n");
			mprintf("    at the beginning of the trajectory due to the finite differences in the magnetic moment PDE.\n\n");
		}
	}


	if (m_iMode == ROA_MODE_ANALYZE) {

		mprintf("\n");
		if (m_bPola) {
			mprintf("    You now need to supply the data files (.emp) which have previously been generated\n");
			mprintf("    in gathering runs.\n\n");

			mprintf("    Assuming that TRAVIS was started with the field-free trajectory.\n\n");
		}

		tr = new CROATrajectory;
		tr->SetHistory(1);
		m_iaTrajKinds[0] = (int)m_iaInvTrajKinds.size();
		m_iaInvTrajKinds.push_back(0);
		m_oaTrajectories.push_back(tr);
/*	_again0:
		AskString_ND("    Enter data file name without field: ",&buf);
		tr->m_sFileName = (const char*)buf;*/
		tr->m_sFileName = (const char*)g_sInputTraj;

		mprintf("    Trying to open data file \"%s\"...\n",(const char*)g_sInputTraj);

		if ((mystricmp(GetFileExtension((const char*)g_sInputTraj),"bqb") != 0) &&
			(mystricmp(GetFileExtension((const char*)g_sInputTraj),"emp") != 0) &&
			(mystricmp(GetFileExtension((const char*)g_sInputTraj),"blist") != 0) &&
			(mystricmp(GetFileExtension((const char*)g_sInputTraj),"voronoi") != 0)) {
			eprintf("Error: File extension of data files should be .bqb, .emp, .blist, or .voronoi\n");
			if (m_bPola)
				eprintf("Please start TRAVIS with the field-free data file generated in a gathering run before.\n");
			else
				eprintf("Please start TRAVIS with the data file generated in a gathering run before.\n");
			abort();
		//	goto _again0;
		}
		if (!tr->OpenFile()) {
			eprintf("Error while opening data file.\n");
			if (m_bPola)
				eprintf("Please start TRAVIS with the field-free data file generated in a gathering run before.\n");
			else
				eprintf("Please start TRAVIS with the data file generated in a gathering run before.\n");
			abort();
		//	goto _again0;
		}

		if (tr->m_bBQB)
			steps = tr->m_pBQB->GetTotalFrameCount();
		else
			steps = -1;

		if (!CheckTrajColumns(tr)) {
			eprintf("\nSomething was wrong with the supplied data file (see above). Aborting.\n");
			abort();
		}

		if (m_bPola) {

			mprintf("    Polarizability calculation was requested. You need to supply the previously generated\n");
			mprintf("    data files from the trajectories with external electric field.\n\n");

			for (z=0;z<(m_bAniso?3:1);z++) {
				for (z2=0;z2<(m_bCentral?2:1);z2++) {
					tr = new CROATrajectory;
					tr->SetHistory(1);
					i = z*2+z2+1;
					m_iaTrajKinds[i] = (int)m_iaInvTrajKinds.size();
					m_iaInvTrajKinds.push_back(i);
					m_oaTrajectories.push_back(tr);
		_again:
					AskString_ND("    Enter data file name for field %c%s: ",
						&buf,'X'+z,m_bCentral?((z2==0)?" positive":" negative"):"");
					tr->m_sFileName = (const char*)buf;
					if ((mystricmp(GetFileExtension((const char*)buf),"bqb") != 0) &&
						(mystricmp(GetFileExtension((const char*)buf),"emp") != 0) &&
						(mystricmp(GetFileExtension((const char*)buf),"blist") != 0)) {
						eprintf("Error: File extension of data files is .bqb, .emp, or .blist\n\n");
						goto _again;
					}
					if (!tr->OpenFile())
						goto _again;
					if (tr->m_bBQB)
						if ((steps != -1) && (tr->m_pBQB->GetTotalFrameCount() != -1) && (tr->m_pBQB->GetTotalFrameCount() != steps))
							eprintf("Warning: Step count (%d) differs from first data file (%d).\n",
								tr->m_pBQB->GetTotalFrameCount(),steps);
					if (!CheckTrajColumns(tr)) {
						eprintf("\nSomething was wrong with the supplied data file (see above). Aborting.\n");
						abort();
					}
				}
			}
		}

	/*	if (steps != -1)
			mprintf("\n    First trajectory contains %d steps. Good.\n",steps);
		else 
			mprintf("\n    No index, could not determine step count.\n"); */

		mprintf("    Rewinding trajectories...\n");

		for (z=0;z<(int)m_oaTrajectories.size();z++)
			m_oaTrajectories[z]->Rewind();

	/*	if (g_pBQBEngine != NULL) {
			if (g_pBQBEngine->m_pExtrapolator != NULL)
				g_pBQBEngine->m_pExtrapolator->Clear();
			if (g_pBQBEngine->m_pExtrapolatorCorr != NULL)
				g_pBQBEngine->m_pExtrapolatorCorr->Clear();
			if (g_pBQBEngine->m_pExtrapolatorXYZ != NULL)
				g_pBQBEngine->m_pExtrapolatorXYZ->Clear();
		}*/

		mprintf("\n    Reading first time step from %d trajectories...\n",(int)m_oaTrajectories.size());

		g_bVoroIntegrateCharge = true;
		g_bVoroIntegrateDipoleMoment = true;
		g_bVoroIntegrateQuadrupoleMoment = true;
		if (m_bMagMom) {
			g_bVoroIntegrateTotalCurrent = true;
			g_bVoroIntegrateMagneticMoment = true;
			g_bCubeTimeDev = true;
			g_bUseVelocities = true;
		} else {
			g_bVoroIntegrateTotalCurrent = false;
			g_bVoroIntegrateMagneticMoment = false;
			g_bCubeTimeDev = false;
			g_bUseVelocities = false;
		}

		for (z=0;z<(int)m_oaTrajectories.size();z++) {
			mprintf("      %s ...\n",m_oaTrajectories[z]->m_sFileName.c_str());
			if (!m_oaTrajectories[z]->ReadTimeStep()) {
				eprintf("Error: Could not read time step.\n");
				abort();
			}
			for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++) {
				rm = new CROAMolecule();
				m_oaTrajectories[z]->m_oaMolecules.push_back(rm);
			}
		}

		for (z=0;z<g_oaMolecules.GetSize();z++) {
			((CMolecule*)g_oaMolecules[z])->m_iDipoleMode = 7;
			if (m_bMagMom)
				((CMolecule*)g_oaMolecules[z])->m_iMagneticDipoleMode = 7;
		}

		mprintf("\n    Computing dipole moment of first molecule from first time step of each trajectory...\n");

		mprintf("\n");
		for (z=0;z<(int)m_oaTrajectories.size();z++) {
			mprintf("      Trajectory %d: ",z+1);
			m_oaTrajectories[z]->GetTimeStep(1)->CalcDipoles(false);
			CMolecule *m = (CMolecule *)g_oaMolecules[0];
			mprintf("( %12.8f | %12.8f | %12.8f ) Debye.\n",
				((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[0],
				((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[1],
				((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[2]
			);
			va.Add(((CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole);
		}
		mprintf("\n");

		if (m_bPola && m_bCentral) {

			mprintf("    Checking symmetry of finite differences (should be all around 0.5)...\n");
			for (z=0;z<(m_bAniso?3:1);z++)
				mprintf("      Field +/- %c: ( %10.8f | %10.8f | %10.8f )\n",
					'X'+z,
					(va[0][0]-va[z*2+2][0])/(va[z*2+1][0]-va[z*2+2][0]),
					(va[0][1]-va[z*2+2][1])/(va[z*2+1][1]-va[z*2+2][1]),
					(va[0][2]-va[z*2+2][2])/(va[z*2+1][2]-va[z*2+2][2])
					);

			mprintf("\n");
		}
	}


	if (m_iMode == ROA_MODE_COMBINED) {

		mprintf("\n");
		if (m_bPola)
			mprintf("    Assuming that TRAVIS was started with the field-free trajectory.\n\n");
		tr = new CROATrajectory;
		if (m_bMagMom)
			tr->SetHistory(3);
		else
			tr->SetHistory(1);
		m_iaTrajKinds[0] = (int)m_iaInvTrajKinds.size();
		m_iaInvTrajKinds.push_back(0);
		m_oaTrajectories.push_back(tr);
/*	_again0b:
		AskString_ND("    Enter trajectory file name without field: ",&buf);
		tr->m_sFileName = (const char*)buf;*/
		tr->m_sFileName = (const char*)g_sInputTraj;

		mprintf("    Trying to open volumetric trajectory file \"%s\"...\n",(const char*)g_sInputTraj);

		if (!tr->OpenFile()) {
			eprintf("Error while opening volumetric trajectory file.\n");
			if (m_bPola)
				eprintf("Please start TRAVIS with the field-free volumetric trajectory file.\n");
			else
				eprintf("Please start TRAVIS with the volumetric trajectory file.\n");
			abort();
		//	goto _again0b;
		}

		if (tr->m_bBQB)
			steps = tr->m_pBQB->GetTotalFrameCount();
		else
			steps = -1;

		if (m_bPola) {

			mprintf("\n");
			mprintf("    Polarizability calculation was requested. You need to supply the\n");
			mprintf("    volumetric data trajectories with external electric field.\n\n");

			for (z=0;z<(m_bAniso?3:1);z++) {
				for (z2=0;z2<(m_bCentral?2:1);z2++) {
					tr = new CROATrajectory;
					if (m_bMagMom)
						tr->SetHistory(3);
					else
						tr->SetHistory(1);
					i = z*2+z2+1;
					m_iaTrajKinds[i] = (int)m_iaInvTrajKinds.size();
					m_iaInvTrajKinds.push_back(i);
					m_oaTrajectories.push_back(tr);
		_againb:
					AskString_ND("    Enter trajectory file name for field %c%s: ",
						&buf,'X'+z,m_bCentral?((z2==0)?" positive":" negative"):"");
					tr->m_sFileName = (const char*)buf;
					if (!tr->OpenFile())
						goto _againb;
					if (tr->m_bBQB)
						if ((steps != -1) && (tr->m_pBQB->GetTotalFrameCount() != -1) && (tr->m_pBQB->GetTotalFrameCount() != steps)) {
							eprintf("Warning: Step count (%d) differs from first trajectory (%d).\n",
								tr->m_pBQB->GetTotalFrameCount(),steps);
							//return false;
					}
				}
			}
		}

	/*	if (steps != -1)
			mprintf("\n    First trajectory contains %d steps. Good.\n",steps);
		else 
			mprintf("\n    No index, could not determine step count.\n"); */

		mprintf("\n    Reading first time step from %d trajectories...\n",(int)m_oaTrajectories.size());

		for (z=0;z<(int)m_oaTrajectories.size();z++) {

			mprintf("      %s ... ",m_oaTrajectories[z]->m_sFileName.c_str());
			if (!m_oaTrajectories[z]->ReadTimeStep()) {
				eprintf("Error: Could not read time step.\n");
				abort();
			}

			mprintf("%d x %d x %d bins, %9.4f x %9.4f x %9.4f pm.\n",
				m_oaTrajectories[z]->GetTimeStep(1)->m_pVolumetricData->m_iRes[0],
				m_oaTrajectories[z]->GetTimeStep(1)->m_pVolumetricData->m_iRes[1],
				m_oaTrajectories[z]->GetTimeStep(1)->m_pVolumetricData->m_iRes[2],
				m_oaTrajectories[z]->GetTimeStep(1)->m_pVolumetricData->m_fMaxVal[0],
				m_oaTrajectories[z]->GetTimeStep(1)->m_pVolumetricData->m_fMaxVal[1],
				m_oaTrajectories[z]->GetTimeStep(1)->m_pVolumetricData->m_fMaxVal[2]);

			for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++) {
				rm = new CROAMolecule();
				m_oaTrajectories[z]->m_oaMolecules.push_back(rm);
			}
		}

		for (z=0;z<g_oaMolecules.GetSize();z++) {
			((CMolecule*)g_oaMolecules[z])->m_iDipoleMode = 7;
			if (m_bMagMom)
				((CMolecule*)g_oaMolecules[z])->m_iMagneticDipoleMode = 7;
		}

		g_bVoroIntegrateCharge = true;
		g_bVoroIntegrateDipoleMoment = true;
		g_bVoroIntegrateQuadrupoleMoment = true;
		g_bVoroIntegrateTotalCurrent = false;
		g_bVoroIntegrateMagneticMoment = false;
		g_bCubeTimeDev = false;

		mprintf("    Setting up Voronoi integration...\n");
		if (!PrepareVori(m_oaTrajectories[0]->GetTimeStep(1),true,true,true)) {
			eprintf("Error while setting up Voronoi integration.\n");
			return false;
		}
		
		mprintf("\n    Computing dipole moment of first molecule from first time step of each trajectory...\n");

		mprintf("\n");
		for (z=0;z<(int)m_oaTrajectories.size();z++) {
			mprintf("    Trajectory %d: ",z+1);
			g_pTetraPak->ProcessStep(m_oaTrajectories[z]->GetTimeStep(1),false);
			m_oaTrajectories[z]->GetTimeStep(1)->CalcDipoles(false);
			CMolecule *m = (CMolecule *)g_oaMolecules[0];
			mprintf("( %12.8f | %12.8f | %12.8f ) Debye.\n",
				((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[0],
				((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[1],
				((CSingleMolecule *)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole[2]
			);
			va.Add(((CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]])->m_vDipole);
		}
		mprintf("\n");

		if (m_bPola && m_bCentral) {

			mprintf("    Checking symmetry of finite differences (should be all around 0.5)...\n");
			for (z=0;z<(m_bAniso?3:1);z++)
				mprintf("      Field +/- %c: ( %10.8f | %10.8f | %10.8f )\n",
					'X'+z,
					(va[0][0]-va[z*2+2][0])/(va[z*2+1][0]-va[z*2+2][0]),
					(va[0][1]-va[z*2+2][1])/(va[z*2+1][1]-va[z*2+2][1]),
					(va[0][2]-va[z*2+2][2])/(va[z*2+1][2]-va[z*2+2][2])
					);

			mprintf("\n");
		}

		g_bVoroIntegrateCharge = true;
		g_bVoroIntegrateDipoleMoment = true;
		g_bVoroIntegrateQuadrupoleMoment = true;
		if (m_bMagMom) {
			g_bVoroIntegrateTotalCurrent = true;
			g_bVoroIntegrateMagneticMoment = true;
			g_bCubeTimeDev = true;
			g_bUseVelocities = true;
		} else {
			g_bVoroIntegrateTotalCurrent = false;
			g_bVoroIntegrateMagneticMoment = false;
			g_bCubeTimeDev = false;
			g_bUseVelocities = false;
		}
	}

	mprintf(WHITE,"\n    <<< New Spectroscopy Module <<<\n\n");

	return true;
}


void CROAEngine::MainLoop() {

	if (m_iMode == ROA_MODE_COMBINED)
		MainLoop_Combined();
	else if (m_iMode == ROA_MODE_GATHER)
		MainLoop_Gather();
	else if (m_iMode == ROA_MODE_ANALYZE)
		MainLoop_Analyze();
}


void CROAEngine::MainLoop_Combined() {

	int z, z2, tsteps;
	unsigned long eta, t0;
	double tfs;
	CROATrajectory *tr;
	CTimeStep *ts, tts;
	CxString buf;
	CROAMolecule *rm;
	CSingleMolecule *sm;


	mprintf(WHITE,"\n    >>> Processing Trajectories - Combined Mode >>>\n\n");

	for (z=0;z<(int)m_oaTrajectories.size();z++) {

		tr = m_oaTrajectories[z];

		mprintf(WHITE,"*** Processing trajectory %d / %d ***\n\n",
			z+1,(int)m_oaTrajectories.size());

		mprintf("    Rewinding trajectory...\n");

		tr->Rewind();

	/*	if (g_pBQBEngine != NULL) {
			if (g_pBQBEngine->m_pExtrapolator != NULL)
				g_pBQBEngine->m_pExtrapolator->Clear();
			if (g_pBQBEngine->m_pExtrapolatorCorr != NULL)
				g_pBQBEngine->m_pExtrapolatorCorr->Clear();
			if (g_pBQBEngine->m_pExtrapolatorXYZ != NULL)
				g_pBQBEngine->m_pExtrapolatorXYZ->Clear();
		}*/

		if (m_bMagMom)
			tr->m_vaAtomVelocities.resize(g_iGesAtomCount);

		if (g_iBeginStep != 0) {

			mprintf("    Fast-forwarding to step %d...\n",g_iBeginStep+1);
			if (tr->m_bBQB) 
			{
				if (!tr->m_pBQB->SeekFrame(g_iBeginStep)) {
					eprintf("      Error while fast-forwarding trajectory (BQB).\n");
					return;
				}

			} else {

				mprintf(WHITE,"      [");
				tfs = g_iBeginStep/60.0;
				for (z2=0;z2<g_iBeginStep;z2++) {
					if (!tts.SkipCube(m_oaTrajectories[0]->m_pFile)) {
						eprintf("      Error while fast-forwarding trajectory (CUBE).\n");
						return;
					}
					if (fmod(z2,tfs) < 1)
						mprintf(WHITE,"#");
				}
				mprintf(WHITE,"]");
				mprintf("\n");
			}
		}

		tsteps = 0;
		if (m_bMagMom && !m_bPattern3) {

			mprintf("\n");
			mprintf("    Reading first 2 time steps for finite differences...\n");

			tsteps += 2;

			for (z2=0;z2<2;z2++)
				if (!tr->ReadTimeStep())
					return;

			if (z == 0)
				PrepareCurrentDensity( tr->GetTimeStep(1),  tr->GetTimeStep(1), true  );
			else
				PrepareCurrentDensity( tr->GetTimeStep(1),  tr->GetTimeStep(1), false );
			PrepareCurrentDensity( tr->GetTimeStep(0),  tr->GetTimeStep(1), false );
			PrepareCurrentDensity( tr->GetTimeStep(-1), tr->GetTimeStep(1), false );
		}

		mprintf("\n");
		mprintf("    Starting run...\n\n");

		t0 = (unsigned long)time(NULL);
		g_iSteps = 0;

		while (true) {

			mprintf("    Step %5lu",g_iSteps+1);
			if ((g_iSteps != 0) && (g_iTrajSteps > 0) && ((time(NULL) - t0) > 5)) {
				eta = (unsigned long)(((double)time(NULL) - t0) / tsteps * ((double)MAX(long(0),g_iTrajSteps - ((long)tsteps))));
				FormatTime(eta,&buf);
				mprintf(": ETA %s",(const char*)buf);
			}
			mprintf("\n");

			if ((g_iStride != 1) && (g_iSteps != 0)) {

				tsteps += g_iStride-1;
				mprintf("      Skipping steps ");
				mprintf("[");

				if (g_iStride > 4) {

					for (z2=0;z2<g_iStride-3;z2++) {
						if (!tr->SkipTimeStep()) {
							mprintf("\n    Reached end of file.\n");
							goto _end;
						}
						mprintf("S");
					}

					for (z2=0;z2<2;z2++) {
						if (!tr->ReadTimeStep()) {
							mprintf("\n    Reached end of file.\n");
							goto _end;
						}
						mprintf("R");
					}

				} else if (g_iStride == 4) {

					if (!tr->SkipTimeStep()) {
						mprintf("\n    Reached end of file.\n");
						goto _end;
					}
					mprintf("S");

					for (z2=0;z2<2;z2++) {
						if (!tr->ReadTimeStep()) {
							mprintf("\n    Reached end of file.\n");
							goto _end;
						}
						mprintf("R");
					}

				} else if (g_iStride == 3) {

					for (z2=0;z2<2;z2++) {
						if (!tr->ReadTimeStep()) {
							mprintf("\n    Reached end of file.\n");
							goto _end;
						}
						mprintf("R");
					}

				} else if (g_iStride == 2) {

					if (!tr->ReadTimeStep()) {
						mprintf("\n    Reached end of file.\n");
						goto _end;
					}
					mprintf("R");
				}
				mprintf("]\n");
			}

			if ((g_iMaxStep > 0) && (tsteps >= g_iMaxStep)) {
				mprintf("\nReached step limit.\n");
				break;
			}

			if (m_bPattern3) {

				mprintf("      Reading block of 3 steps ");
				mprintf("[");
				for (z2=0;z2<3;z2++) {
					if (!ReadStep(false)) {
						mprintf("\n    Reached end of file.\n");
						goto _end;
					}
					mprintf("R");
					tsteps++;
				}
				mprintf("]\n");

				if (m_bMagMom && (g_iSteps == 0)) {

					if (z == 0)
						PrepareCurrentDensity( tr->GetTimeStep(1),  tr->GetTimeStep(1), true  );
					else
						PrepareCurrentDensity( tr->GetTimeStep(1),  tr->GetTimeStep(1), false );
					PrepareCurrentDensity( tr->GetTimeStep(0),  tr->GetTimeStep(1), false );
					PrepareCurrentDensity( tr->GetTimeStep(-1), tr->GetTimeStep(1), false );
				}

			} else {

				mprintf("      Reading step...\n");
				if (!ReadStep(true)) {
					mprintf("\nReached end of file.\n");
					break;
				}
				tsteps++;
			}

			if (m_bMagMom) {
				tr->CalcVelocities();
				mprintf("      Current ");
				CalculateCurrentDensity(tr->GetTimeStep(-1),tr->GetTimeStep(0),tr->GetTimeStep(1),g_iSteps==0);
				mprintf("\n");
			}

			ts = tr->GetTimeStep(0);

			mprintf("      Integrating...\n");
			g_pTetraPak->ProcessStep(ts,false);

			if (m_bMagMom)
				ts->CalcCenterVelocities();

			ts->CalcDipoles(false);

			if (m_bMagMom)
				ts->CalcMagneticDipoles();

			for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++) {
				rm = tr->m_oaMolecules[z2];
				sm = (CSingleMolecule*)g_oaSingleMolecules[z2];

				rm->m_faCharge.push_back(sm->m_fCharge);          // Unit: e
				rm->m_vaElDip.push_back(sm->m_vDipole);           // Unit: Debye
				rm->m_maElQuad.push_back(sm->m_mQuadrupole);      // Unit: Debye*pm

				if (m_bMagMom) {
					rm->m_vaElCurr.push_back(sm->m_vCurrent);         // Unit: MB/pm
					rm->m_vaMagDip.push_back(sm->m_vMagneticDipole);  // Unit: MB
				}
			}

			if (m_bMagMom)
				for (z2=0;z2<g_iGesAtomCount;z2++)
					tr->m_vaAtomVelocities[z2].Add(ts->m_vaVelocities[z2]);

			g_iSteps++;

			if ((g_iMaxStep > 0) && (tsteps >= g_iMaxStep)) {
				mprintf("\n    Reached step limit.\n");
				break;
			}
		}
_end:
		if (g_iStride != 1)
			mprintf("\n    Processed %lu steps out of %d input steps.\n\n",g_iSteps,tsteps);
		else
			mprintf("\n    Processed %lu steps.\n\n",g_iSteps);

		for (z2=0;z2<(int)tr->m_oaTSHistory.size();z2++)
			delete tr->m_oaTSHistory[z2];
		tr->m_oaTSHistory.clear();

		m_iStepsProcessed = (int)tr->m_oaMolecules[0]->m_faCharge.size();
	}

	mprintf(WHITE,"    Finished processing %d trajectories:\n\n",(int)m_oaTrajectories.size());

	for (z=0;z<(int)m_oaTrajectories.size();z++)
		mprintf("      Trajectory %d: %d steps.\n",z+1,(int)m_oaTrajectories[z]->m_oaMolecules[0]->m_faCharge.size());

	for (z=0;z<(int)m_oaTrajectories.size();z++) {
		if ((int)m_oaTrajectories[z]->m_oaMolecules[0]->m_faCharge.size() != m_iStepsProcessed) {
			eprintf("\nError: Step count of trajectories does not match.\n");
			abort();
		}
	}

	mprintf(WHITE,"\n    <<< Processing Trajectories - Combined Mode <<<\n\n");

	Finish();
}


void CROAEngine::MainLoop_Analyze() {

	int z, z2, z3, tsteps;
	unsigned long eta, t0;
	CxString buf;
	CROATrajectory *tr;
	CTimeStep *ts;
	CROAMolecule *rm;
	CSingleMolecule *sm;
	FILE *a;


	mprintf(WHITE,"\n    >>> Processing Trajectories - Analyze >>>\n\n");

	if (m_bDumpMolecularProps) {
		for (z=0;z<(int)m_oaTrajectories.size();z++) {
			buf.sprintf("bqb_molecular_integrals_traj%d.csv",z+1);
			a = OpenFileWrite((const char*)buf, true);
			mfprintf(a,"#Step; Molecule; SingleMolecule; Charge; ElDipX; ElDipY; ElDipZ; ElQuadXX; ElQuadXY; ElQuadXZ; ElQuadYX; ElQuadYY; ElQuadYZ; ElQuadZX; ElQuadZY; ElQuadZZ");
			if (m_bMagMom)
				mfprintf(a,"; ElCurrX; ElCurrY; ElCurrZ; MagDipX; MagDipY; MagDipZ");
			mfprintf(a,"\n");
			fflush(a);
			m_fMolIntegralFiles.push_back(a);
		}
	}

	mprintf("    Rewinding trajectories...\n\n");

	for (z=0;z<(int)m_oaTrajectories.size();z++)
		m_oaTrajectories[z]->Rewind();

/*	if (g_pBQBEngine != NULL) {
		if (g_pBQBEngine->m_pExtrapolator != NULL)
			g_pBQBEngine->m_pExtrapolator->Clear();
		if (g_pBQBEngine->m_pExtrapolatorCorr != NULL)
			g_pBQBEngine->m_pExtrapolatorCorr->Clear();
		if (g_pBQBEngine->m_pExtrapolatorXYZ != NULL)
			g_pBQBEngine->m_pExtrapolatorXYZ->Clear();
	}*/

	if (g_iBeginStep != 0) {
		mprintf("    Fast-forwarding to step %d...\n",g_iBeginStep+1);
		for (z=0;z<(int)m_oaTrajectories.size();z++) {
			if (!m_oaTrajectories[z]->m_pBQB->SeekFrame(g_iBeginStep)) {
				eprintf("      Error while fast-forwarding trajectory %d (BQB).\n",z+1);
				return;
			}
		}
	}

	if (m_bMagMom)
		m_oaTrajectories[0]->m_vaAtomVelocities.resize(g_iGesAtomCount);

	m_oaTrajectories[0]->m_vaCOMCoord.resize(g_oaSingleMolecules.GetSize());

	g_bVoroIntegrateCharge = true;
	g_bVoroIntegrateDipoleMoment = true;
	g_bVoroIntegrateQuadrupoleMoment = true;
	if (m_bMagMom) {
		g_bVoroIntegrateTotalCurrent = true;
		g_bVoroIntegrateMagneticMoment = true;
		g_bCubeTimeDev = true;
		g_bUseVelocities = true;
	} else {
		g_bVoroIntegrateTotalCurrent = false;
		g_bVoroIntegrateMagneticMoment = false;
		g_bCubeTimeDev = false;
		g_bUseVelocities = false;
	}

	t0 = (unsigned long)time(NULL);
	g_iSteps = 0;
	tsteps = 0;

	while (true) {

		if ((g_iSteps % 10) == 0) {
			if ((g_iSteps % 500) == 0) {
				if (g_iSteps != 0) {
					if ((g_iTrajSteps > 0) && ((time(NULL) - t0) > 5)) {
						eta = (unsigned long)(((double)time(NULL) - t0) / tsteps * ((double)MAX(long(0),g_iTrajSteps - ((long)tsteps))));
						FormatTime(eta,&buf);
						mprintf(" ETA %s\nStep %5lu ",(const char*)buf,g_iSteps);
					} else
						mprintf("\nStep %5lu ",g_iSteps);
				} else
					mprintf("Step %5d ",1);
			} else
				mprintf(".");
		}

		if ((g_iStride != 1) && (g_iSteps != 0)) {
			tsteps += g_iStride-1;
			for (z=0;z<g_iStride-1;z++) {
				if (!SkipStep()) {
					mprintf("\nReached end of file.\n");
					goto _end;
				}
			}
		}

		if (!ReadStep(false)) {
			mprintf("\nReached end of file.\n");
			break;
		}

		tsteps++;


			for (z=0;z<(int)m_oaTrajectories.size();z++) {

				tr = m_oaTrajectories[z];
				ts = tr->GetTimeStep(1);

				if (m_bMagMom)
					ts->CalcCenterVelocities();

				ts->CalcDipoles(false);

				if (m_bMagMom)
					ts->CalcMagneticDipoles();

				if (m_bDumpMolecularProps) {
					g_fMolIntegralFile = m_fMolIntegralFiles[z];
					DumpMolecularProps();
				}

				for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++) {
					rm = tr->m_oaMolecules[z2];
					sm = (CSingleMolecule*)g_oaSingleMolecules[z2];

					rm->m_faCharge.push_back(sm->m_fCharge);          // Unit: e
					rm->m_vaElDip.push_back(sm->m_vDipole);           // Unit: Debye
					rm->m_maElQuad.push_back(sm->m_mQuadrupole);      // Unit: Debye*pm

					if (m_bMagMom) {
						rm->m_vaElCurr.push_back(sm->m_vCurrent);         // Unit: MB/pm
						rm->m_vaMagDip.push_back(sm->m_vMagneticDipole);  // Unit: MB
					}
				}
			}

		ts = m_oaTrajectories[0]->GetTimeStep(1);

		for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++) {
			sm = (CSingleMolecule*)g_oaSingleMolecules[z2];
			m_oaTrajectories[0]->m_vaCOMCoord[z2].Add(ts->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[sm->m_baAtomIndex.GetSize()-1])->GetAt(1)]);
		}

		if (m_bMagMom)
			for (z2=0;z2<g_iGesAtomCount;z2++)
				m_oaTrajectories[0]->m_vaAtomVelocities[z2].Add(ts->m_vaVelocities[z2]);

		g_iSteps++;
		if ((g_iMaxStep > 0) && (tsteps >= g_iMaxStep)) {
			mprintf("\nReached step limit.\n");
			break;
		}
	}
_end:
	if (g_iStride != 1)
		mprintf("\nProcessed %lu steps out of %d input steps.\n\n",g_iSteps,tsteps);
	else
		mprintf("\nProcessed %lu steps.\n\n",g_iSteps);

	m_iStepsProcessed = g_iSteps;

	if (m_bDumpMolecularProps) {
		mprintf("  Molecular properties have been written to the following files:\n");
		for (z=0;z<(int)m_oaTrajectories.size();z++) {
			mprintf("    bqb_molecular_integrals_traj%d.csv\n",z+1);
			fclose(m_fMolIntegralFiles[z]);
		}
		m_fMolIntegralFiles.clear();
	}

	if (m_bReverseTraj) {
		mprintf("\n    Reversing direction of input trajectory...\n");
		for (z=0;z<(int)m_oaTrajectories.size();z++) {
			tr = m_oaTrajectories[z];
			for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++) {
				rm = tr->m_oaMolecules[z2];
				std::reverse(rm->m_faCharge.begin(),rm->m_faCharge.end());
				std::reverse(rm->m_vaElDip.begin(),rm->m_vaElDip.end());
				std::reverse(rm->m_maElQuad.begin(),rm->m_maElQuad.end());
				if (m_bMagMom) {
					std::reverse(rm->m_vaElCurr.begin(),rm->m_vaElCurr.end());
					std::reverse(rm->m_vaMagDip.begin(),rm->m_vaMagDip.end());
					for (z3=0;z3<(int)rm->m_vaElCurr.size();z3++)
						rm->m_vaElCurr[z3] *= -1.0;
					for (z3=0;z3<(int)rm->m_vaMagDip.size();z3++)
						rm->m_vaMagDip[z3] *= -1.0;
				}
			}
		}
	}

	mprintf(WHITE,"\n    <<< Processing Trajectories - Analyze <<<\n\n");

		Finish();
}


void CROAEngine::MainLoop_Gather() {

	int z, z2, z3, tsteps;
	unsigned long eta, t0;
	double tfs;
	CxString buf;
	FILE *a;
	CBQBFile *bqwrite;
	CBQBTrajectoryFrame *btf;
	CBQBTrajectoryFrameColumn *bco;
	CTimeStep *ts, tts;
	std::vector<unsigned char> ia;


	a = NULL;

	mprintf(WHITE,"\n    >>> Processing Trajectory - Gather >>>\n\n");

	mprintf("    Rewinding trajectory...\n");

	m_oaTrajectories[0]->Rewind();

/*	if (g_pBQBEngine != NULL) {
		if (g_pBQBEngine->m_pExtrapolator != NULL)
			g_pBQBEngine->m_pExtrapolator->Clear();
		if (g_pBQBEngine->m_pExtrapolatorCorr != NULL)
			g_pBQBEngine->m_pExtrapolatorCorr->Clear();
		if (g_pBQBEngine->m_pExtrapolatorXYZ != NULL)
			g_pBQBEngine->m_pExtrapolatorXYZ->Clear();
	}*/

	if (g_iBeginStep != 0) {
		mprintf("    Fast-forwarding to step %d...\n",g_iBeginStep+1);
		if (m_oaTrajectories[0]->m_bBQB) {
			if (!m_oaTrajectories[0]->m_pBQB->SeekFrame(g_iBeginStep)) {
				eprintf("      Error while fast-forwarding trajectory (BQB).\n");
				return;
			}
		} else {
			mprintf(WHITE,"      [");
			tfs = g_iBeginStep/60.0;
			for (z2=0;z2<g_iBeginStep;z2++) {
				if (!tts.SkipCube(m_oaTrajectories[0]->m_pFile)) {
					eprintf("      Error while fast-forwarding trajectory (CUBE).\n");
					return;
				}
				if (fmod(z2,tfs) < 1)
					mprintf(WHITE,"#");
			}
			mprintf(WHITE,"]");
			mprintf("\n");
		}
	}

	tsteps = 0;

	if (m_bMagMom) {
		
		if (!m_bPattern3) {

			mprintf("    Reading first 2 time steps for finite differences...\n");

			tsteps += 2;

			for (z=0;z<2;z++)
				if (!ReadStep(true))
					return;

			mprintf("\n");
		}

		if (g_iTrajSteps != -1)
			mprintf("    Starting run (%d remaining steps)...\n\n",g_iTrajSteps-2);
		else
			mprintf("    Starting run...\n\n");

	} else {

		if (g_iTrajSteps != -1)
			mprintf("    Starting run (%d steps)...\n\n",g_iTrajSteps);
		else
			mprintf("    Starting run...\n\n");
	}

	bqwrite = new CBQBFile(*g_pBQBInterface);
	bqwrite->OpenWriteReplace("properties.emp");

	btf = new CBQBTrajectoryFrame(*g_pBQBInterface);
	btf->m_iAtomCount = g_iGesAtomCount;
	btf->AddColumn(BQB_TYPE_STRING,"Label");
	btf->AddColumn(BQB_TYPE_DOUBLE,"PosX");    // Unit: pm
	btf->AddColumn(BQB_TYPE_DOUBLE,"PosY");    // Unit: pm
	btf->AddColumn(BQB_TYPE_DOUBLE,"PosZ");    // Unit: pm

	btf->AddColumn(BQB_TYPE_DOUBLE,"Vol");     // Unit: pm^3

	if (m_bMagMom) {
		btf->AddColumn(BQB_TYPE_DOUBLE,"VelX");    // Unit: pm/ps
		btf->AddColumn(BQB_TYPE_DOUBLE,"VelY");    // Unit: pm/ps
		btf->AddColumn(BQB_TYPE_DOUBLE,"VelZ");    // Unit: pm/ps
	}

	btf->AddColumn(BQB_TYPE_DOUBLE,"EChg");    // Unit: e
	btf->AddColumn(BQB_TYPE_DOUBLE,"ECCg");    // Unit: e
	btf->AddColumn(BQB_TYPE_DOUBLE,"EDipX");   // Unit: e*pm
	btf->AddColumn(BQB_TYPE_DOUBLE,"EDipY");   // Unit: e*pm
	btf->AddColumn(BQB_TYPE_DOUBLE,"EDipZ");   // Unit: e*pm
	btf->AddColumn(BQB_TYPE_DOUBLE,"EChZX");   // Unit: pm
	btf->AddColumn(BQB_TYPE_DOUBLE,"EChZY");   // Unit: pm
	btf->AddColumn(BQB_TYPE_DOUBLE,"EChZZ");   // Unit: pm
	btf->AddColumn(BQB_TYPE_DOUBLE,"EQXX");    // Unit: e*nm^2
	btf->AddColumn(BQB_TYPE_DOUBLE,"EQXY");    // Unit: e*nm^2
	btf->AddColumn(BQB_TYPE_DOUBLE,"EQXZ");    // Unit: e*nm^2
	btf->AddColumn(BQB_TYPE_DOUBLE,"EQYX");    // Unit: e*nm^2
	btf->AddColumn(BQB_TYPE_DOUBLE,"EQYY");    // Unit: e*nm^2
	btf->AddColumn(BQB_TYPE_DOUBLE,"EQYZ");    // Unit: e*nm^2
	btf->AddColumn(BQB_TYPE_DOUBLE,"EQZX");    // Unit: e*nm^2
	btf->AddColumn(BQB_TYPE_DOUBLE,"EQZY");    // Unit: e*nm^2
	btf->AddColumn(BQB_TYPE_DOUBLE,"EQZZ");    // Unit: e*nm^2

	if (m_bMagMom) {
		btf->AddColumn(BQB_TYPE_DOUBLE,"ECurX");   // Unit: MB/pm
		btf->AddColumn(BQB_TYPE_DOUBLE,"ECurY");   // Unit: MB/pm
		btf->AddColumn(BQB_TYPE_DOUBLE,"ECurZ");   // Unit: MB/pm
		btf->AddColumn(BQB_TYPE_DOUBLE,"MDipX");   // Unit: MB
		btf->AddColumn(BQB_TYPE_DOUBLE,"MDipY");   // Unit: MB
		btf->AddColumn(BQB_TYPE_DOUBLE,"MDipZ");   // Unit: MB
	}

	btf->m_oaColumns[0]->m_aString.resize(g_iGesAtomCount);
	for (z=1;z<(int)btf->m_oaColumns.size();z++)
		btf->m_oaColumns[z]->m_aReal.resize(g_iGesAtomCount);

	btf->m_pCellMatrix = new CxDMatrix3();
	for (z=0;z<9;z++)
		(*btf->m_pCellMatrix)[z] = g_mBoxFromOrtho[z];

	for (z=0;z<g_iGesAtomCount;z++)
		btf->m_oaColumns[0]->m_aString[z] = (const char*)m_oaTrajectories[0]->GetTimeStep(1)->m_paLabels[z];

	if (m_bGatherWriteCSV) {

		a = OpenFileWrite("properties.csv",true);

		mfprintf(a,"#Step");
		mfprintf(a,"; AtomID");
		mfprintf(a,"; Label");
		mfprintf(a,"; PosX; PosY; PosZ");

		mfprintf(a,"; Vol");

		if (m_bMagMom)
			mfprintf(a,"; VelX; VelY; VelZ");

		mfprintf(a,"; Charge; CoreCharge");
		mfprintf(a,"; ElDipX; ElDipY; ElDipZ");
		mfprintf(a,"; ElChZX; ElChZY; ElChZZ");
		mfprintf(a,"; ElQuadXX; ElQuadXY; ElQuadXZ; ElQuadYX; ElQuadYY; ElQuadYZ; ElQuadZX; ElQuadZY; ElQuadZZ");

		if (m_bMagMom) {
			mfprintf(a,"; ElCurrX; ElCurrY; ElCurrZ");
			mfprintf(a,"; MagDipX; MagDipY; MagDipZ");
		}
		mfprintf(a,"\n");
	}

	g_bVoroIntegrateCharge = true;
	g_bVoroIntegrateDipoleMoment = true;
	g_bVoroIntegrateQuadrupoleMoment = true;
	if (m_bMagMom) {
		g_bVoroIntegrateTotalCurrent = true;
		g_bVoroIntegrateMagneticMoment = true;
		g_bCubeTimeDev = true;
	} else {
		g_bVoroIntegrateTotalCurrent = false;
		g_bVoroIntegrateMagneticMoment = false;
		g_bCubeTimeDev = false;
	}

	t0 = (unsigned long)time(NULL);
	g_iSteps = 0;
	while (true) {

		mprintf("    Step %5lu",g_iSteps+1);
		if ((g_iSteps != 0) && (g_iTrajSteps > 0)) {
			eta = (unsigned long)(((double)time(NULL) - t0) / tsteps * ((double)MAX(long(0),g_iTrajSteps - tsteps)));
			FormatTime(eta,&buf);
			mprintf(": ETA %s",(const char*)buf);
		}
		mprintf("\n");

		if ((g_iStride != 1) && (g_iSteps != 0)) {

			mprintf("      Skipping steps ");
			mprintf("[");

			tsteps += g_iStride-1;

			if (g_iStride > 4) {

				for (z2=0;z2<g_iStride-3;z2++) {
					if (!SkipStep()) {
						mprintf("\n    Reached end of file.\n");
						goto _end;
					}
					mprintf("S");
				}

				for (z2=0;z2<2;z2++) {
					if (!ReadStep(false)) {
						mprintf("\n    Reached end of file.\n");
						goto _end;
					}
					mprintf("R");
				}

			} else if (g_iStride == 4) {

				if (!SkipStep()) {
					mprintf("\n    Reached end of file.\n");
					goto _end;
				}
				mprintf("S");

				for (z2=0;z2<2;z2++) {
					if (!ReadStep(false)) {
						mprintf("\n    Reached end of file.\n");
						goto _end;
					}
					mprintf("R");
				}

			} else if (g_iStride == 3) {

				for (z2=0;z2<2;z2++) {
					if (!ReadStep(false)) {
						mprintf("\n    Reached end of file.\n");
						goto _end;
					}
					mprintf("R");
				}

			} else if (g_iStride == 2) {

				if (!ReadStep(false)) {
					mprintf("\n    Reached end of file.\n");
					goto _end;
				}
				mprintf("R");
			}
			mprintf("]\n");
		}

		if ((g_iMaxStep > 0) && (tsteps >= g_iMaxStep)) {
			mprintf("\nReached step limit.\n");
			break;
		}

		if (m_bPattern3) {

			mprintf("      Reading block of 3 steps ");
			mprintf("[");
			for (z2=0;z2<3;z2++) {
				if (!ReadStep(false)) {
					mprintf("\n    Reached end of file.\n");
					goto _end;
				}
				mprintf("R");
				tsteps++;
			}
			mprintf("]\n");

		} else {

			if (!ReadStep(true)) {
				mprintf("\nReached end of file.\n");
				break;
			}
			tsteps++;
		}


		if (m_bMagMom) {

			this->CalcVelocities();

			mprintf("      Current ");
			CalculateCurrentDensity(m_oaTrajectories[0]->GetTimeStep(-1),m_oaTrajectories[0]->GetTimeStep(0),m_oaTrajectories[0]->GetTimeStep(1),false);
			mprintf("\n");
		}

		ts = m_oaTrajectories[0]->GetTimeStep(0);

		mprintf("      Integrating...\n");
		g_pTetraPak->ProcessStep(ts,false);

		if (m_bGatherWriteCSV) {
			for (z=0;z<g_iGesAtomCount;z++) {
				mfprintf(a,"%6lu",g_iSteps+1);
				mfprintf(a,"; %4d",z+1);
				mfprintf(a,"; %-3s",(const char*)ts->m_paLabels[z]);
				mfprintf(a,"; %17.10f; %17.10f; %17.10f",
					ts->m_vaCoords[z][0],
					ts->m_vaCoords[z][1],
					ts->m_vaCoords[z][2]
				);
				mfprintf(a,"; %17.10f",ts->m_faVolume[z]);

				if (m_bMagMom)
					mfprintf(a,"; %17.10f; %17.10f; %17.10f",
						ts->m_vaVelocities[z][0],
						ts->m_vaVelocities[z][1],
						ts->m_vaVelocities[z][2]
					);

				mfprintf(a,"; %17.10f",ts->m_faCharge[z]);
				mfprintf(a,"; %17.10f",ts->m_faCoreCharge[z]);
				mfprintf(a,"; %17.10f; %17.10f; %17.10f",
					ts->m_dipoleMoments[z][0],
					ts->m_dipoleMoments[z][1],
					ts->m_dipoleMoments[z][2]
				);
				mfprintf(a,"; %17.10f; %17.10f; %17.10f",
					ts->m_chargeCenters[z][0],
					ts->m_chargeCenters[z][1],
					ts->m_chargeCenters[z][2]
				);
				mfprintf(a,"; %17.10f; %17.10f; %17.10f; %17.10f; %17.10f; %17.10f; %17.10f; %17.10f; %17.10f",
					ts->m_maQuadTensor[z][0],
					ts->m_maQuadTensor[z][1],
					ts->m_maQuadTensor[z][2],
					ts->m_maQuadTensor[z][3],
					ts->m_maQuadTensor[z][4],
					ts->m_maQuadTensor[z][5],
					ts->m_maQuadTensor[z][6],
					ts->m_maQuadTensor[z][7],
					ts->m_maQuadTensor[z][8]
				);
				if (m_bMagMom) {
					mfprintf(a,"; %17.10f; %17.10f; %17.10f",
						ts->m_totalCurrents[z][0],
						ts->m_totalCurrents[z][1],
						ts->m_totalCurrents[z][2]
					);
					mfprintf(a,"; %17.10f; %17.10f; %17.10f",
						ts->m_magneticDipoleMoments[z][0],
						ts->m_magneticDipoleMoments[z][1],
						ts->m_magneticDipoleMoments[z][2]
					);
				}
				mfprintf(a,"\n");
			}
			fflush(a);
		}

		bco = btf->GetColumn("PosX");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_vaCoords[z][0];
		bco = btf->GetColumn("PosY");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_vaCoords[z][1];
		bco = btf->GetColumn("PosZ");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_vaCoords[z][2];

		bco = btf->GetColumn("Vol");   for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_faVolume[z];

		if (m_bMagMom) {
			bco = btf->GetColumn("VelX");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_vaVelocities[z][0];
			bco = btf->GetColumn("VelY");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_vaVelocities[z][1];
			bco = btf->GetColumn("VelZ");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_vaVelocities[z][2];
		}

		bco = btf->GetColumn("EChg");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_faCharge[z];
		bco = btf->GetColumn("ECCg");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_faCoreCharge[z];
		bco = btf->GetColumn("EDipX"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_dipoleMoments[z][0];
		bco = btf->GetColumn("EDipY"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_dipoleMoments[z][1];
		bco = btf->GetColumn("EDipZ"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_dipoleMoments[z][2];
		bco = btf->GetColumn("EChZX"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_chargeCenters[z][0];
		bco = btf->GetColumn("EChZY"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_chargeCenters[z][1];
		bco = btf->GetColumn("EChZZ"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_chargeCenters[z][2];
		bco = btf->GetColumn("EQXX");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_maQuadTensor[z][0];
		bco = btf->GetColumn("EQXY");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_maQuadTensor[z][1];
		bco = btf->GetColumn("EQXZ");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_maQuadTensor[z][2];
		bco = btf->GetColumn("EQYX");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_maQuadTensor[z][3];
		bco = btf->GetColumn("EQYY");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_maQuadTensor[z][4];
		bco = btf->GetColumn("EQYZ");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_maQuadTensor[z][5];
		bco = btf->GetColumn("EQZX");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_maQuadTensor[z][6];
		bco = btf->GetColumn("EQZY");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_maQuadTensor[z][7];
		bco = btf->GetColumn("EQZZ");  for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_maQuadTensor[z][8];
		if (m_bMagMom) {
			bco = btf->GetColumn("ECurX"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_totalCurrents[z][0];
			bco = btf->GetColumn("ECurY"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_totalCurrents[z][1];
			bco = btf->GetColumn("ECurZ"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_totalCurrents[z][2];
			bco = btf->GetColumn("MDipX"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_magneticDipoleMoments[z][0];
			bco = btf->GetColumn("MDipY"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_magneticDipoleMoments[z][1];
			bco = btf->GetColumn("MDipZ"); for (z=0;z<g_iGesAtomCount;z++) bco->m_aReal[z] = ts->m_magneticDipoleMoments[z][2];
		}

		ia.clear();
		btf->WriteFrame(&ia);

		bqwrite->CreateShortFrame(BQB_FRAMETYPE_TRAJ,0,g_iSteps+1);
		bqwrite->PushPayload(ia);
		bqwrite->FinalizeFrame(NULL);

		g_iSteps++;

		if ((g_iMaxStep > 0) && (tsteps >= g_iMaxStep)) {
			mprintf("\nReached step limit.\n");
			break;
		}
	}
_end:
	if (g_iStride != 1)
		mprintf("\nProcessed %lu steps out of %d input steps.\n\n",g_iSteps,tsteps);
	else
		mprintf("\nProcessed %lu steps.\n\n",g_iSteps);

	mprintf("Writing BQB index...\n");
	bqwrite->WriteIndexFrame(true,NULL);
	bqwrite->Close();
	delete bqwrite;
	delete btf;

	mprintf("Finished writing BQB file.\n\n");

	if (m_bGatherWriteCSV)
		fclose(a);

	if (m_bGatherCheck) {

		mprintf("Now re-opening BQB file to check.\n");

		bqwrite = new CBQBFile(*g_pBQBInterface);
		bqwrite->OpenRead("properties.emp");

		btf = new CBQBTrajectoryFrame(*g_pBQBInterface);

		a = OpenFileWrite("properties_out.csv",true);

		z = 0;
		while (bqwrite->ReadFrame()) {
			if (bqwrite->GetFrameType() != BQB_FRAMETYPE_TRAJ)
				continue;
			mprintf("Frame %d...\n",z+1);
			btf->ReadFrame(bqwrite->GetFramePayload());
			if (z == 0) {
				mfprintf(a,"#Step; AtomID");
				for (z2=0;z2<(int)btf->m_oaColumns.size();z2++)
					mfprintf(a,"; %s",btf->m_oaColumns[z2]->m_sLabel.c_str());
				mfprintf(a,"\n");
			}
			for (z2=0;z2<(int)btf->m_iAtomCount;z2++) {
				mfprintf(a,"%6d; %4d",z+1,z2+1);
				mfprintf(a,"; %-3s",(const char*)btf->m_oaColumns[0]->m_aString[z2].c_str());
				for (z3=1;z3<(int)btf->m_oaColumns.size();z3++)
					mfprintf(a,"; %17.10f",btf->m_oaColumns[z3]->m_aReal[z2]);
				mfprintf(a,"\n");
			}
			fflush(a);
			z++;
		}

		fclose(a);

		bqwrite->Close();
		delete bqwrite;
		delete btf;
	}

	mprintf(WHITE,"\n    <<< Processing Trajectory - Gather <<<\n\n");
}


CTimeStep* CROAEngine::GetROAStep(int traj, int step) {

	if (m_iaTrajKinds[traj] == -1) {
		eprintf("CROAEngine::GetROAStep(): Error: Trajectory %d not present.\n",traj);
		abort();
	}

	return m_oaTrajectories[m_iaTrajKinds[traj]]->GetTimeStep(step);
}


bool CROAEngine::ReadStep(bool verbose) {

	int z;

	if (verbose)
		mprintf("      Reading [");
	for (z=0;z<(int)m_oaTrajectories.size();z++) {
		if (!m_oaTrajectories[z]->ReadTimeStep()) {
			mprintf("\nCROAEngine::ReadStep(): Could not read further time step from trajectory %d (%s).\n",
				z+1,m_oaTrajectories[z]->m_sFileName.c_str());
			return false;
		}
		if (verbose)
			mprintf("#");
	}
	if (verbose)
		mprintf("]\n");

	return true;
}


bool CROAEngine::SkipStep() {

	int z;

	for (z=0;z<(int)m_oaTrajectories.size();z++) {
		if (!m_oaTrajectories[z]->SkipTimeStep()) {
			mprintf("\nCROAEngine::SkipStep(): Could not skip further time step in trajectory %d (%s).\n",
				z+1,m_oaTrajectories[z]->m_sFileName.c_str());
			return false;
		}
	}

	return true;
}


void CROAEngine::CalcVelocities() {

	int z;

	for (z=0;z<(int)m_oaTrajectories.size();z++)
		m_oaTrajectories[z]->CalcVelocities();
}


void CROAEngine::Finish() {

	int z, z2, z3, o, ic;
	bool b, b2;
	double f, tfs, cfac;
	CROATrajectory *tn, *txp, *txn, *typ, *tyn, *tzp, *tzn;
	CROAMolecule *rn, *rxp, *rxn, *ryp, *ryn, *rzp, *rzn;
	FILE *a;
	CSmoothener smo;
	CAutoCorrelation *ac;
	CxDoubleArray dain, daout;
	CxDoubleArray daacf[64];
	CxDoubleArray daspec[64];
	CxDoubleArray *tda;
	CxString buf, buf2;
	CFFT *fft;
	std::vector<double> temp;


	mprintf(WHITE,"\n    >>> Finishing Analysis >>>\n\n");

	if (g_bRamanFromPolarizability) {
		tn=NULL; txp=NULL; txn=NULL; typ=NULL; tyn=NULL; tzp=NULL; tzn=NULL;
		rn=NULL; rxp=NULL; rxn=NULL; ryp=NULL; ryn=NULL; rzp=NULL; rzn=NULL;
		b = false;
		tda = NULL;
		cfac = 0;
		ic = 0;
		a = NULL;
		tn = m_oaTrajectories[0];
		goto _ramanfrompola;
	}

	if (m_fPDEInfo != NULL)
		fclose(m_fPDEInfo);

	b = false;
	tda = NULL;
	cfac = 0;
	ic = 0;
	a = NULL;

	tn=NULL; txp=NULL; txn=NULL; typ=NULL; tyn=NULL; tzp=NULL; tzn=NULL;
	rn=NULL; rxp=NULL; rxn=NULL; ryp=NULL; ryn=NULL; rzp=NULL; rzn=NULL;

	if (m_bPola) {
		if (m_bCentral && m_bAniso) {
			tn = m_oaTrajectories[0];
			txp = m_oaTrajectories[1];
			txn = m_oaTrajectories[2];
			typ = m_oaTrajectories[3];
			tyn = m_oaTrajectories[4];
			tzp = m_oaTrajectories[5];
			tzn = m_oaTrajectories[6];
		} else if (!m_bCentral && m_bAniso) {
			tn = m_oaTrajectories[0];
			txp = m_oaTrajectories[1];
			typ = m_oaTrajectories[2];
			tzp = m_oaTrajectories[3];
		} else {
			eprintf("Error: Not yet implemented.\n");
			return;
		}
	} else
		tn = m_oaTrajectories[0];

	if (m_bMagMom) {

		mprintf("    Computing power spectrum of total system...\n");
		ac = new CAutoCorrelation();
		ac->Init(m_iStepsProcessed-2,m_oaObservations[0]->m_iDepth,g_bACFFFT);
		dain.SetSize(m_iStepsProcessed-2);
		daout.SetSize(m_oaObservations[0]->m_iDepth);

		fft = new CFFT();

		fft->PrepareFFT_C2C(2*(m_oaObservations[0]->m_iDepth+m_oaObservations[0]->m_iZeroPadding));
		temp.resize(2*(m_oaObservations[0]->m_iDepth+m_oaObservations[0]->m_iZeroPadding));

		daacf[0].SetSize(m_oaObservations[0]->m_iDepth);
		for (z2=0;z2<m_oaObservations[0]->m_iDepth;z2++)
			daacf[0][z2] = 0;

		tfs = g_iGesAtomCount / 60.0;
		mprintf(WHITE,"      [");
		for (z=0;z<g_iGesAtomCount;z++) {
			if (fmod(z,tfs) < 1.0)
				mprintf(WHITE,"#");
			for (z3=0;z3<3;z3++) {
				for (z2=0;z2<m_iStepsProcessed-2;z2++)
					dain[z2] = tn->m_vaAtomVelocities[z][z2][z3];
				ac->AutoCorrelate(&dain,&daout);
				for (z2=0;z2<m_oaObservations[0]->m_iDepth;z2++)
					daacf[0][z2] += daout[z2] * ((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_fMass;
			}
		}
		mprintf(WHITE,"]\n");

		for (z3=0;z3<m_oaObservations[0]->m_iDepth;z3++)
			daacf[0][z3] /= (double)3*g_iGesAtomCount;

		for (z3=0;z3<m_oaObservations[0]->m_iDepth;z3++)
			temp[z3] = daacf[0][z3];

		for (z3=m_oaObservations[0]->m_iDepth;z3<m_oaObservations[0]->m_iDepth+m_oaObservations[0]->m_iZeroPadding;z3++)
			temp[z3] = 0.0;

		for (z3=0;z3<m_oaObservations[0]->m_iDepth;z3++)
			temp[z3] *= pow2(cos((double)z3 / (m_oaObservations[0]->m_iDepth - 1) / 2.0 * Pi));

		o = m_oaObservations[0]->m_iDepth+m_oaObservations[0]->m_iZeroPadding;
		for (z3=1;z3<o;z3++)
			temp[o+z3] = temp[o-z3];
		temp[o] = 0.0;

		for (z3=0;z3<2*(m_oaObservations[0]->m_iDepth+m_oaObservations[0]->m_iZeroPadding);z3++) {
			fft->m_pInput[2*z3] = temp[z3];
			fft->m_pInput[2*z3+1] = 0.0;
		}

		fft->DoFFT();

		daspec[0].SetSize(m_oaObservations[0]->m_iSpecLength);

		for (z3=0;z3<m_oaObservations[0]->m_iSpecLength;z3++)
			daspec[0][z3] = fft->m_pOutput[2*z3];

		for(z3 = 0; z3 < m_oaObservations[0]->m_iSpecLength; z3++)
			daspec[0][z3] *= 7.211349e-9 * g_fTimestepLength * g_iStride; // Output in K*cm

		delete ac;
		delete fft;

		mprintf("    Writing power spectrum to \"power_spectrum.csv\"...\n");

		a = OpenFileWrite("power_spectrum.csv",true);
		mfprintf(a,"#Wavenumber / cm^-1; Power Spectrum / K cm; Integral / K\n");
		f = 0;
		for (z=0;z<m_oaObservations[0]->m_iSpecLength;z++) {
			f += daspec[0][z] * m_oaObservations[0]->m_fSpecResolution;
			mfprintf(a,"%.2f; %.10G; %.10G\n",m_oaObservations[0]->m_fSpecResolution*z,daspec[0][z],f);
		}
		fclose(a);

		if (m_bWriteACFs || m_oaObservations[0]->m_bSaveACF) {
			mprintf("    Writing velocity autocorrelation function to \"power_vacf.csv\"...\n");
			a = OpenFileWrite("power_vacf.csv",true);
			mfprintf(a,"#Depth / fs; VACF\n");
			for (z=0;z<m_oaObservations[0]->m_iDepth;z++)
				mfprintf(a,"%.2f; %.10G\n",z*g_fTimestepLength*g_iStride,daacf[0][z]);
			fclose(a);
		}

		mprintf("\n    Assuming 3n = %d degrees of freedom, the average system temperature is %.2f K.\n\n",
			3 * g_iGesAtomCount, f );
	}


	if (m_bPola) {

	_again:
		if (!b && m_bDumpMolecularProps && m_bCentral) {
			mprintf("    Writing results to \"central.csv\"...\n");
			a = OpenFileWrite("central.csv",true);
			mfprintf(a,"#Step");
			for (z=0;z<g_oaSingleMolecules.GetSize();z++) {
				mfprintf(a,"; SM%d_PolElDipXX; SM%d_PolElDipXY; SM%d_PolElDipXZ; SM%d_PolElDipYX; SM%d_PolElDipYY; SM%d_PolElDipYZ; SM%d_PolElDipZX; SM%d_PolElDipZY; SM%d_PolElDipZZ",
					z+1,z+1,z+1,z+1,z+1,z+1,z+1,z+1,z+1);
				if (m_bMagMom)
					mfprintf(a,"; SM%d_PolMagDipXX; SM%d_PolMagDipXY; SM%d_PolMagDipXZ; SM%d_PolMagDipYX; SM%d_PolMagDipYY; SM%d_PolMagDipYZ; SM%d_PolMagDipZX; SM%d_PolMagDipZY; SM%d_PolMagDipZZ",
						z+1,z+1,z+1,z+1,z+1,z+1,z+1,z+1,z+1);
			}
			mfprintf(a,"\n");
			if (m_bMagMom)
				ic = 18;
			else
				ic = 9;

			tda = new CxDoubleArray[g_oaSingleMolecules.GetSize()*ic];
			for (z=0;z<g_oaSingleMolecules.GetSize()*ic;z++)
				tda[z].SetSize(m_iStepsProcessed);
		}
		mprintf("    Computing polarizabilities for %d molecules and %d steps...\n",
			g_oaSingleMolecules.GetSize(),m_iStepsProcessed);
		tfs = g_oaSingleMolecules.GetSize() / 60.0;
		if (g_oaSingleMolecules.GetSize() > 1)
			mprintf(WHITE,"      [");
		for (z=0;z<g_oaSingleMolecules.GetSize();z++) {

			if (g_oaSingleMolecules.GetSize() > 1)
				if (fmod(z,tfs) < 1.0)
					mprintf(WHITE,"#");

			if (m_bCentral && m_bAniso) {
				rn = tn->m_oaMolecules[z];
				rxp = txp->m_oaMolecules[z];
				rxn = txn->m_oaMolecules[z];
				ryp = typ->m_oaMolecules[z];
				ryn = tyn->m_oaMolecules[z];
				rzp = tzp->m_oaMolecules[z];
				rzn = tzn->m_oaMolecules[z];
				cfac = 2.0;
			} else if (!m_bCentral && m_bAniso) {
				rn = tn->m_oaMolecules[z];
				rxp = txp->m_oaMolecules[z];
				rxn = tn->m_oaMolecules[z];
				ryp = typ->m_oaMolecules[z];
				ryn = tn->m_oaMolecules[z];
				rzp = tzp->m_oaMolecules[z];
				rzn = tn->m_oaMolecules[z];
				cfac = 1.0;
			}

			rn->m_maPolElDip.resize(m_iStepsProcessed);
			rn->m_taPolElQuad.resize(m_iStepsProcessed,CDTensor3(3,3,3));
			if (m_bMagMom)
				rn->m_maPolMagDip.resize(m_iStepsProcessed);

			for (z2=0;z2<m_iStepsProcessed;z2++) {

				// Electric Dipole-Electric Dipole Polarizability
				rn->m_maPolElDip[z2](0,0) = (rxp->m_vaElDip[z2][0] - rxn->m_vaElDip[z2][0]) / (cfac * m_fFieldStrength);
				rn->m_maPolElDip[z2](0,1) = (rxp->m_vaElDip[z2][1] - rxn->m_vaElDip[z2][1]) / (cfac * m_fFieldStrength);
				rn->m_maPolElDip[z2](0,2) = (rxp->m_vaElDip[z2][2] - rxn->m_vaElDip[z2][2]) / (cfac * m_fFieldStrength);

				rn->m_maPolElDip[z2](1,0) = (ryp->m_vaElDip[z2][0] - ryn->m_vaElDip[z2][0]) / (cfac * m_fFieldStrength);
				rn->m_maPolElDip[z2](1,1) = (ryp->m_vaElDip[z2][1] - ryn->m_vaElDip[z2][1]) / (cfac * m_fFieldStrength);
				rn->m_maPolElDip[z2](1,2) = (ryp->m_vaElDip[z2][2] - ryn->m_vaElDip[z2][2]) / (cfac * m_fFieldStrength);

				rn->m_maPolElDip[z2](2,0) = (rzp->m_vaElDip[z2][0] - rzn->m_vaElDip[z2][0]) / (cfac * m_fFieldStrength);
				rn->m_maPolElDip[z2](2,1) = (rzp->m_vaElDip[z2][1] - rzn->m_vaElDip[z2][1]) / (cfac * m_fFieldStrength);
				rn->m_maPolElDip[z2](2,2) = (rzp->m_vaElDip[z2][2] - rzn->m_vaElDip[z2][2]) / (cfac * m_fFieldStrength);

				if (!b && m_bDumpMolecularProps && m_bCentral) {
					tda[z*ic  ][z2] = (rn->m_vaElDip[z2][0] - rxn->m_vaElDip[z2][0]) / (rxp->m_vaElDip[z2][0] - rxn->m_vaElDip[z2][0]);
					tda[z*ic+1][z2] = (rn->m_vaElDip[z2][1] - rxn->m_vaElDip[z2][1]) / (rxp->m_vaElDip[z2][1] - rxn->m_vaElDip[z2][1]);
					tda[z*ic+2][z2] = (rn->m_vaElDip[z2][2] - rxn->m_vaElDip[z2][2]) / (rxp->m_vaElDip[z2][2] - rxn->m_vaElDip[z2][2]);
					tda[z*ic+3][z2] = (rn->m_vaElDip[z2][0] - ryn->m_vaElDip[z2][0]) / (ryp->m_vaElDip[z2][0] - ryn->m_vaElDip[z2][0]);
					tda[z*ic+4][z2] = (rn->m_vaElDip[z2][1] - ryn->m_vaElDip[z2][1]) / (ryp->m_vaElDip[z2][1] - ryn->m_vaElDip[z2][1]);
					tda[z*ic+5][z2] = (rn->m_vaElDip[z2][2] - ryn->m_vaElDip[z2][2]) / (ryp->m_vaElDip[z2][2] - ryn->m_vaElDip[z2][2]);
					tda[z*ic+6][z2] = (rn->m_vaElDip[z2][0] - rzn->m_vaElDip[z2][0]) / (rzp->m_vaElDip[z2][0] - rzn->m_vaElDip[z2][0]);
					tda[z*ic+7][z2] = (rn->m_vaElDip[z2][1] - rzn->m_vaElDip[z2][1]) / (rzp->m_vaElDip[z2][1] - rzn->m_vaElDip[z2][1]);
					tda[z*ic+8][z2] = (rn->m_vaElDip[z2][2] - rzn->m_vaElDip[z2][2]) / (rzp->m_vaElDip[z2][2] - rzn->m_vaElDip[z2][2]);
				}

				for (z3=0;z3<9;z3++)
					rn->m_maPolElDip[z2][z3] *= 1.0 / DIP_EPM2DEBYE / EFIELD_AU2VPM; // Conversion Debye*e*a0/Eh to e*pm^2/V


				// Electric Dipole-Electric Quadrupole Polarizability
				rn->m_taPolElQuad[z2](0,0,0) = (rxp->m_maElQuad[z2](0,0) - rxn->m_maElQuad[z2](0,0)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](0,0,1) = (rxp->m_maElQuad[z2](0,1) - rxn->m_maElQuad[z2](0,1)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](0,0,2) = (rxp->m_maElQuad[z2](0,2) - rxn->m_maElQuad[z2](0,2)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](0,1,0) = (rxp->m_maElQuad[z2](1,0) - rxn->m_maElQuad[z2](1,0)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](0,1,1) = (rxp->m_maElQuad[z2](1,1) - rxn->m_maElQuad[z2](1,1)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](0,1,2) = (rxp->m_maElQuad[z2](1,2) - rxn->m_maElQuad[z2](1,2)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](0,2,0) = (rxp->m_maElQuad[z2](2,0) - rxn->m_maElQuad[z2](2,0)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](0,2,1) = (rxp->m_maElQuad[z2](2,1) - rxn->m_maElQuad[z2](2,1)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](0,2,2) = (rxp->m_maElQuad[z2](2,2) - rxn->m_maElQuad[z2](2,2)) / (cfac * m_fFieldStrength);

				rn->m_taPolElQuad[z2](1,0,0) = (ryp->m_maElQuad[z2](0,0) - ryn->m_maElQuad[z2](0,0)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](1,0,1) = (ryp->m_maElQuad[z2](0,1) - ryn->m_maElQuad[z2](0,1)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](1,0,2) = (ryp->m_maElQuad[z2](0,2) - ryn->m_maElQuad[z2](0,2)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](1,1,0) = (ryp->m_maElQuad[z2](1,0) - ryn->m_maElQuad[z2](1,0)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](1,1,1) = (ryp->m_maElQuad[z2](1,1) - ryn->m_maElQuad[z2](1,1)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](1,1,2) = (ryp->m_maElQuad[z2](1,2) - ryn->m_maElQuad[z2](1,2)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](1,2,0) = (ryp->m_maElQuad[z2](2,0) - ryn->m_maElQuad[z2](2,0)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](1,2,1) = (ryp->m_maElQuad[z2](2,1) - ryn->m_maElQuad[z2](2,1)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](1,2,2) = (ryp->m_maElQuad[z2](2,2) - ryn->m_maElQuad[z2](2,2)) / (cfac * m_fFieldStrength);

				rn->m_taPolElQuad[z2](2,0,0) = (rzp->m_maElQuad[z2](0,0) - rzn->m_maElQuad[z2](0,0)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](2,0,1) = (rzp->m_maElQuad[z2](0,1) - rzn->m_maElQuad[z2](0,1)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](2,0,2) = (rzp->m_maElQuad[z2](0,2) - rzn->m_maElQuad[z2](0,2)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](2,1,0) = (rzp->m_maElQuad[z2](1,0) - rzn->m_maElQuad[z2](1,0)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](2,1,1) = (rzp->m_maElQuad[z2](1,1) - rzn->m_maElQuad[z2](1,1)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](2,1,2) = (rzp->m_maElQuad[z2](1,2) - rzn->m_maElQuad[z2](1,2)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](2,2,0) = (rzp->m_maElQuad[z2](2,0) - rzn->m_maElQuad[z2](2,0)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](2,2,1) = (rzp->m_maElQuad[z2](2,1) - rzn->m_maElQuad[z2](2,1)) / (cfac * m_fFieldStrength);
				rn->m_taPolElQuad[z2](2,2,2) = (rzp->m_maElQuad[z2](2,2) - rzn->m_maElQuad[z2](2,2)) / (cfac * m_fFieldStrength);

				for (z3=0;z3<27;z3++)
					rn->m_taPolElQuad[z2][z3] *= 1.0 / DIP_EPM2DEBYE / EFIELD_AU2VPM; // Conversion Debye*pm*e*a0/Eh to e*pm^3/V


				if (m_bMagMom) {

					// Electric Dipole-Magnetic Dipole Polarizability
					rn->m_maPolMagDip[z2](0,0) = (rxp->m_vaMagDip[z2][0] - rxn->m_vaMagDip[z2][0]) / (cfac * m_fFieldStrength);
					rn->m_maPolMagDip[z2](0,1) = (rxp->m_vaMagDip[z2][1] - rxn->m_vaMagDip[z2][1]) / (cfac * m_fFieldStrength);
					rn->m_maPolMagDip[z2](0,2) = (rxp->m_vaMagDip[z2][2] - rxn->m_vaMagDip[z2][2]) / (cfac * m_fFieldStrength);

					rn->m_maPolMagDip[z2](1,0) = (ryp->m_vaMagDip[z2][0] - ryn->m_vaMagDip[z2][0]) / (cfac * m_fFieldStrength);
					rn->m_maPolMagDip[z2](1,1) = (ryp->m_vaMagDip[z2][1] - ryn->m_vaMagDip[z2][1]) / (cfac * m_fFieldStrength);
					rn->m_maPolMagDip[z2](1,2) = (ryp->m_vaMagDip[z2][2] - ryn->m_vaMagDip[z2][2]) / (cfac * m_fFieldStrength);

					rn->m_maPolMagDip[z2](2,0) = (rzp->m_vaMagDip[z2][0] - rzn->m_vaMagDip[z2][0]) / (cfac * m_fFieldStrength);
					rn->m_maPolMagDip[z2](2,1) = (rzp->m_vaMagDip[z2][1] - rzn->m_vaMagDip[z2][1]) / (cfac * m_fFieldStrength);
					rn->m_maPolMagDip[z2](2,2) = (rzp->m_vaMagDip[z2][2] - rzn->m_vaMagDip[z2][2]) / (cfac * m_fFieldStrength);

					if (!b && m_bDumpMolecularProps && m_bCentral) {
						tda[z*ic+ 9][z2] = (rn->m_vaMagDip[z2][0] - rxn->m_vaMagDip[z2][0]) / (rxp->m_vaMagDip[z2][0] - rxn->m_vaMagDip[z2][0]);
						tda[z*ic+10][z2] = (rn->m_vaMagDip[z2][1] - rxn->m_vaMagDip[z2][1]) / (rxp->m_vaMagDip[z2][1] - rxn->m_vaMagDip[z2][1]);
						tda[z*ic+11][z2] = (rn->m_vaMagDip[z2][2] - rxn->m_vaMagDip[z2][2]) / (rxp->m_vaMagDip[z2][2] - rxn->m_vaMagDip[z2][2]);
						tda[z*ic+12][z2] = (rn->m_vaMagDip[z2][0] - ryn->m_vaMagDip[z2][0]) / (ryp->m_vaMagDip[z2][0] - ryn->m_vaMagDip[z2][0]);
						tda[z*ic+13][z2] = (rn->m_vaMagDip[z2][1] - ryn->m_vaMagDip[z2][1]) / (ryp->m_vaMagDip[z2][1] - ryn->m_vaMagDip[z2][1]);
						tda[z*ic+14][z2] = (rn->m_vaMagDip[z2][2] - ryn->m_vaMagDip[z2][2]) / (ryp->m_vaMagDip[z2][2] - ryn->m_vaMagDip[z2][2]);
						tda[z*ic+15][z2] = (rn->m_vaMagDip[z2][0] - rzn->m_vaMagDip[z2][0]) / (rzp->m_vaMagDip[z2][0] - rzn->m_vaMagDip[z2][0]);
						tda[z*ic+16][z2] = (rn->m_vaMagDip[z2][1] - rzn->m_vaMagDip[z2][1]) / (rzp->m_vaMagDip[z2][1] - rzn->m_vaMagDip[z2][1]);
						tda[z*ic+17][z2] = (rn->m_vaMagDip[z2][2] - rzn->m_vaMagDip[z2][2]) / (rzp->m_vaMagDip[z2][2] - rzn->m_vaMagDip[z2][2]);
					}

					for (z3=0;z3<9;z3++)
						rn->m_maPolMagDip[z2][z3] *= 0.001 / MAG_EPMMS2MB / EFIELD_AU2VPM; // Conversion MB*e*a0/Eh to e*pm^3/(fs*V)
				}
			}
		}
		if (g_oaSingleMolecules.GetSize() > 1)
			mprintf(WHITE,"]\n");

		if (!b && m_bDumpMolecularProps && m_bCentral) {
			for (z=0;z<m_iStepsProcessed;z++) {
				mfprintf(a,"%d",z+1);
				for (z2=0;z2<g_oaSingleMolecules.GetSize()*ic;z2++)
					mfprintf(a,"; %.10G",tda[z2][z]);
				mfprintf(a,"\n");
			}
			fclose(a);
			delete[] tda;
		}

		if (!b && m_bReplaceOutliers) {

			double tfmi, tfma;
			double tfout;
			bool tb, tb2;
			CFitParabola para;

			tfout = m_fReplaceOutliersThreshold;
			tb2 = false;

	_outagain:
			mprintf("    Detecting outliers...\n");
			for (z=0;z<8;z++) {
				m_iaOutliers[z].resize(m_iStepsProcessed);
				for (z2=0;z2<m_iStepsProcessed;z2++)
					m_iaOutliers[z][z2] = 0;
			}
			tfs = g_oaSingleMolecules.GetSize() / 60.0;
			if (g_oaSingleMolecules.GetSize() > 1)
				mprintf(WHITE,"      [");
			for (z=0;z<g_oaSingleMolecules.GetSize();z++) {

				if (g_oaSingleMolecules.GetSize() > 1)
					if (fmod(z,tfs) < 1.0)
						mprintf(WHITE,"#");

				rn = tn->m_oaMolecules[z];

				for (z2=0;z2<m_iStepsProcessed-4;z2++) {
					tfmi = MIN4(rn->m_faCharge[z2],rn->m_faCharge[z2+1],rn->m_faCharge[z2+3],rn->m_faCharge[z2+4]);
					tfma = MAX4(rn->m_faCharge[z2],rn->m_faCharge[z2+1],rn->m_faCharge[z2+3],rn->m_faCharge[z2+4]);
					if (fabs((tfma-tfmi) / (tfma+tfmi)) < 0.001) {
							tfma += 0.001 * fabs(tfma);
							tfmi -= 0.001 * fabs(tfmi);
					}
					if ((rn->m_faCharge[z2+2] < tfmi-tfout*fabs(tfma-tfmi)) || (rn->m_faCharge[z2+2] > tfma+tfout*fabs(tfma-tfmi))) {
						mprintf("Charge %3d %6d: min=%f, max=%f, tfout=%f, lb=%f, up=%f, value=%f\n",z,z2,tfmi,tfma,tfout,tfmi-tfout*fabs(tfma-tfmi),tfma+tfout*fabs(tfma-tfmi),rn->m_faCharge[z2+2]);
						m_iaOutliers[0][z2]++;
					}
				}

				for (z2=0;z2<m_iStepsProcessed-4;z2++) {
					tb = false;
					for (z3=0;z3<3;z3++) {
						tfmi = MIN4(rn->m_vaElDip[z2][z3],rn->m_vaElDip[z2+1][z3],rn->m_vaElDip[z2+3][z3],rn->m_vaElDip[z2+4][z3]);
						tfma = MAX4(rn->m_vaElDip[z2][z3],rn->m_vaElDip[z2+1][z3],rn->m_vaElDip[z2+3][z3],rn->m_vaElDip[z2+4][z3]);
						if (fabs((tfma-tfmi) / (tfma+tfmi)) < 0.001) {
								tfma += 0.001 * fabs(tfma);
								tfmi -= 0.001 * fabs(tfmi);
						}
						if ((rn->m_vaElDip[z2+2][z3] < tfmi-tfout*fabs(tfma-tfmi)) || (rn->m_vaElDip[z2+2][z3] > tfma+tfout*fabs(tfma-tfmi))) {
							mprintf("ElDip %3d %6d: min=%f, max=%f, tfout=%f, lb=%f, up=%f, value=%f\n",z,z2,tfmi,tfma,tfout,tfmi-tfout*fabs(tfma-tfmi),tfma+tfout*fabs(tfma-tfmi),rn->m_vaElDip[z2+2][z3]);
							tb = true;
						}
					}
					if (tb)
						m_iaOutliers[1][z2]++;
				}

				for (z2=0;z2<m_iStepsProcessed-4;z2++) {
					tb = false;
					for (z3=0;z3<9;z3++) {
						tfmi = MIN4(rn->m_maElQuad[z2][z3],rn->m_maElQuad[z2+1][z3],rn->m_maElQuad[z2+3][z3],rn->m_maElQuad[z2+4][z3]);
						tfma = MAX4(rn->m_maElQuad[z2][z3],rn->m_maElQuad[z2+1][z3],rn->m_maElQuad[z2+3][z3],rn->m_maElQuad[z2+4][z3]);
						if (fabs((tfma-tfmi) / (tfma+tfmi)) < 0.001) {
								tfma += 0.001 * fabs(tfma);
								tfmi -= 0.001 * fabs(tfmi);
						}
						if ((rn->m_maElQuad[z2+2][z3] < tfmi-tfout*fabs(tfma-tfmi)) || (rn->m_maElQuad[z2+2][z3] > tfma+tfout*fabs(tfma-tfmi))) {
							mprintf("ElQuad %3d %6d: min=%f, max=%f, tfout=%f, lb=%f, up=%f, value=%f\n",z,z2,tfmi,tfma,tfout,tfmi-tfout*fabs(tfma-tfmi),tfma+tfout*fabs(tfma-tfmi),rn->m_maElQuad[z2+2][z3]);
							tb = true;
						}
					}
					if (tb)
						m_iaOutliers[2][z2]++;
				}

				for (z2=0;z2<m_iStepsProcessed-4;z2++) {
					tb = false;
					for (z3=0;z3<3;z3++) {
						tfmi = MIN4(rn->m_vaElCurr[z2][z3],rn->m_vaElCurr[z2+1][z3],rn->m_vaElCurr[z2+3][z3],rn->m_vaElCurr[z2+4][z3]);
						tfma = MAX4(rn->m_vaElCurr[z2][z3],rn->m_vaElCurr[z2+1][z3],rn->m_vaElCurr[z2+3][z3],rn->m_vaElCurr[z2+4][z3]);
						if (fabs((tfma-tfmi) / (tfma+tfmi)) < 0.001) {
								tfma += 0.001 * fabs(tfma);
								tfmi -= 0.001 * fabs(tfmi);
						}
						if ((rn->m_vaElCurr[z2+2][z3] < tfmi-2.0*tfout*fabs(tfma-tfmi)) || (rn->m_vaElCurr[z2+2][z3] > tfma+2.0*tfout*fabs(tfma-tfmi))) {
							mprintf("ElCurr %3d %6d: min=%f, max=%f, tfout=%f, lb=%f, up=%f, value=%f\n",z,z2,tfmi,tfma,tfout,tfmi-tfout*fabs(tfma-tfmi),tfma+tfout*fabs(tfma-tfmi),rn->m_vaElCurr[z2+2][z3]);
							tb = true;
						}
					}
					if (tb)
						m_iaOutliers[3][z2]++;
				}

				for (z2=0;z2<m_iStepsProcessed-4;z2++) {
					tb = false;
					for (z3=0;z3<3;z3++) {
						tfmi = MIN4(rn->m_vaMagDip[z2][z3],rn->m_vaMagDip[z2+1][z3],rn->m_vaMagDip[z2+3][z3],rn->m_vaMagDip[z2+4][z3]);
						tfma = MAX4(rn->m_vaMagDip[z2][z3],rn->m_vaMagDip[z2+1][z3],rn->m_vaMagDip[z2+3][z3],rn->m_vaMagDip[z2+4][z3]);
						if (fabs((tfma-tfmi) / (tfma+tfmi)) < 0.001) {
								tfma += 0.001 * fabs(tfma);
								tfmi -= 0.001 * fabs(tfmi);
						}
						if ((rn->m_vaMagDip[z2+2][z3] < tfmi-2.0*tfout*fabs(tfma-tfmi)) || (rn->m_vaMagDip[z2+2][z3] > tfma+2.0*tfout*fabs(tfma-tfmi))) {
							mprintf("MagDip %3d %6d: min=%f, max=%f, tfout=%f, lb=%f, up=%f, value=%f\n",z,z2,tfmi,tfma,tfout,tfmi-tfout*fabs(tfma-tfmi),tfma+tfout*fabs(tfma-tfmi),rn->m_vaMagDip[z2+2][z3]);
							tb = true;
						}
					}
					if (tb)
						m_iaOutliers[4][z2]++;
				}

				for (z2=0;z2<m_iStepsProcessed-4;z2++) {
					tb = false;
					for (z3=0;z3<9;z3++) {
						tfmi = MIN4(rn->m_maPolElDip[z2][z3],rn->m_maPolElDip[z2+1][z3],rn->m_maPolElDip[z2+3][z3],rn->m_maPolElDip[z2+4][z3]);
						tfma = MAX4(rn->m_maPolElDip[z2][z3],rn->m_maPolElDip[z2+1][z3],rn->m_maPolElDip[z2+3][z3],rn->m_maPolElDip[z2+4][z3]);
						if (fabs((tfma-tfmi) / (tfma+tfmi)) < 0.001) {
								tfma += 0.001 * fabs(tfma);
								tfmi -= 0.001 * fabs(tfmi);
						}
						if ((rn->m_maPolElDip[z2+2][z3] < tfmi-tfout*fabs(tfma-tfmi)) || (rn->m_maPolElDip[z2+2][z3] > tfma+tfout*fabs(tfma-tfmi))) {
							mprintf("PolElDip %3d %6d: min=%f, max=%f, tfout=%f, lb=%f, up=%f, value=%f\n",z,z2,tfmi,tfma,tfout,tfmi-tfout*fabs(tfma-tfmi),tfma+tfout*fabs(tfma-tfmi),rn->m_maPolElDip[z2+2][z3]);
							tb = true;
						}
					}
					if (tb)
						m_iaOutliers[5][z2]++;
				}

				for (z2=0;z2<m_iStepsProcessed-4;z2++) {
					tb = false;
					for (z3=0;z3<27;z3++) {
						tfmi = MIN4(rn->m_taPolElQuad[z2][z3],rn->m_taPolElQuad[z2+1][z3],rn->m_taPolElQuad[z2+3][z3],rn->m_taPolElQuad[z2+4][z3]);
						tfma = MAX4(rn->m_taPolElQuad[z2][z3],rn->m_taPolElQuad[z2+1][z3],rn->m_taPolElQuad[z2+3][z3],rn->m_taPolElQuad[z2+4][z3]);
						if (fabs((tfma-tfmi) / (tfma+tfmi)) < 0.001) {
								tfma += 0.001 * fabs(tfma);
								tfmi -= 0.001 * fabs(tfmi);
						}
						if ((rn->m_taPolElQuad[z2+2][z3] < tfmi-tfout*fabs(tfma-tfmi)) || (rn->m_taPolElQuad[z2+2][z3] > tfma+tfout*fabs(tfma-tfmi))) {
							mprintf("PolElQuad %3d %6d: min=%f, max=%f, tfout=%f, lb=%f, up=%f, value=%f\n",z,z2,tfmi,tfma,tfout,tfmi-tfout*fabs(tfma-tfmi),tfma+tfout*fabs(tfma-tfmi),rn->m_taPolElQuad[z2+2][z3]);
							tb = true;
						}
					}
					if (tb)
						m_iaOutliers[6][z2]++;
				}

				for (z2=0;z2<m_iStepsProcessed-4;z2++) {
					tb = false;
					for (z3=0;z3<3;z3++) {
						tfmi = MIN4(rn->m_maPolMagDip[z2][z3],rn->m_maPolMagDip[z2+1][z3],rn->m_maPolMagDip[z2+3][z3],rn->m_maPolMagDip[z2+4][z3]);
						tfma = MAX4(rn->m_maPolMagDip[z2][z3],rn->m_maPolMagDip[z2+1][z3],rn->m_maPolMagDip[z2+3][z3],rn->m_maPolMagDip[z2+4][z3]);
						if (fabs((tfma-tfmi) / (tfma+tfmi)) < 0.001) {
								tfma += 0.001 * fabs(tfma);
								tfmi -= 0.001 * fabs(tfmi);
						}
						if ((rn->m_maPolMagDip[z2+2][z3] < tfmi-2.0*tfout*(tfma-tfmi)) || (rn->m_maPolMagDip[z2+2][z3] > tfma+2.0*tfout*(tfma-tfmi))) {
							mprintf("Charge %3d %6d: min=%f, max=%f, tfout=%f, lb=%f, up=%f, value=%f\n",z,z2,tfmi,tfma,tfout,tfmi-tfout*(tfma-tfmi),tfma+tfout*(tfma-tfmi),rn->m_faCharge[z2+2]);
							tb = true;
						}
					}
					if (tb)
						m_iaOutliers[7][z2]++;
				}
			}
			if (g_oaSingleMolecules.GetSize() > 1)
				mprintf(WHITE,"]\n");

	/*		a = OpenFileWrite("outliers.csv",true);
			fprintf(a,"#Step;  Charge;  ElDip;  ElQuad;  Curr;  MagDip;  PolElDip;  PolElQuad;  PolMagDip\n");
			for (z=0;z<m_iStepsProcessed-4;z++) {
				fprintf(a,"%d",z+3);
				for (z2=0;z2<8;z2++)
					fprintf(a,";  %d",m_iaOutliers[z2][z]);
				fprintf(a,"\n");
			}
			fclose(a);*/

	_outagain2:

			tb = false;
			for (z=0;z<m_iStepsProcessed-4;z++) {

				if (tb2) {

					if ((m_iaOutliers[3][z] < MAX(1,g_oaSingleMolecules.GetSize()/4)) &&
						(m_iaOutliers[4][z] < MAX(1,g_oaSingleMolecules.GetSize()/4)) &&
						(m_iaOutliers[7][z] < MAX(1,g_oaSingleMolecules.GetSize()/4)))
						continue;

					mprintf("      Found magnetic outlier in step %d.\n",z+3);

				} else {

					if ((m_iaOutliers[0][z] < MAX(1,g_oaSingleMolecules.GetSize()/4)) &&
						(m_iaOutliers[1][z] < MAX(1,g_oaSingleMolecules.GetSize()/4)) &&
						(m_iaOutliers[2][z] < MAX(1,g_oaSingleMolecules.GetSize()/4)) &&
						(m_iaOutliers[5][z] < MAX(1,g_oaSingleMolecules.GetSize()/4)) &&
						(m_iaOutliers[6][z] < MAX(1,g_oaSingleMolecules.GetSize()/4)))
						continue;

					mprintf("      Found outlier in step %d.\n",z+3);
					tb = true;
				}


				for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++) {

					rn = tn->m_oaMolecules[z2];

					if (m_bReplaceOutliersElectric) {
						if (!tb2) {

							para.Init4(-2.0,rn->m_faCharge[z],-1.0,rn->m_faCharge[z+1],1.0,rn->m_faCharge[z+3],2.0,rn->m_faCharge[z+4]);
							rn->m_faCharge[z+2] = para.Evaluate(0);

							for (z3=0;z3<3;z3++) {
								para.Init4(-2.0,rn->m_vaElDip[z][z3],-1.0,rn->m_vaElDip[z+1][z3],1.0,rn->m_vaElDip[z+3][z3],2.0,rn->m_vaElDip[z+4][z3]);
								rn->m_vaElDip[z+2][z3] = para.Evaluate(0);
							}

							for (z3=0;z3<9;z3++) {
								para.Init4(-2.0,rn->m_maElQuad[z][z3],-1.0,rn->m_maElQuad[z+1][z3],1.0,rn->m_maElQuad[z+3][z3],2.0,rn->m_maElQuad[z+4][z3]);
								rn->m_maElQuad[z+2][z3] = para.Evaluate(0);
							}

							for (z3=0;z3<9;z3++) {
								para.Init4(-2.0,rn->m_maPolElDip[z][z3],-1.0,rn->m_maPolElDip[z+1][z3],1.0,rn->m_maPolElDip[z+3][z3],2.0,rn->m_maPolElDip[z+4][z3]);
								rn->m_maPolElDip[z+2][z3] = para.Evaluate(0);
							}

							for (z3=0;z3<27;z3++) {
								para.Init4(-2.0,rn->m_taPolElQuad[z][z3],-1.0,rn->m_taPolElQuad[z+1][z3],1.0,rn->m_taPolElQuad[z+3][z3],2.0,rn->m_taPolElQuad[z+4][z3]);
								rn->m_taPolElQuad[z+2][z3] = para.Evaluate(0);
							}
						}
					}

					if (m_bReplaceOutliersMagnetic) {
						if ((z > 0) && (z < m_iStepsProcessed-6)) {

							for (z3=0;z3<3;z3++) {
								para.Init4(-3.0,rn->m_vaMagDip[z-1][z3],-2.0,rn->m_vaMagDip[z][z3],3.0,rn->m_vaMagDip[z+5][z3],4.0,rn->m_vaMagDip[z+6][z3]);
								rn->m_vaMagDip[z+1][z3] = para.Evaluate(-1.0);
								rn->m_vaMagDip[z+2][z3] = para.Evaluate(0);
								rn->m_vaMagDip[z+3][z3] = para.Evaluate(1.0);
								rn->m_vaMagDip[z+4][z3] = para.Evaluate(2.0);
							}

							for (z3=0;z3<3;z3++) {
								para.Init4(-3.0,rn->m_vaElCurr[z-1][z3],-2.0,rn->m_vaElCurr[z][z3],3.0,rn->m_vaElCurr[z+5][z3],4.0,rn->m_vaElCurr[z+6][z3]);
								rn->m_vaElCurr[z+1][z3] = para.Evaluate(-1.0);
								rn->m_vaElCurr[z+2][z3] = para.Evaluate(0);
								rn->m_vaElCurr[z+3][z3] = para.Evaluate(1.0);
								rn->m_vaElCurr[z+4][z3] = para.Evaluate(2.0);
							}

							for (z3=0;z3<9;z3++) {
								para.Init4(-3.0,rn->m_maPolMagDip[z-1][z3],-2.0,rn->m_maPolMagDip[z][z3],3.0,rn->m_maPolMagDip[z+5][z3],4.0,rn->m_maPolMagDip[z+6][z3]);
								rn->m_maPolMagDip[z+1][z3] = para.Evaluate(-1.0);
								rn->m_maPolMagDip[z+2][z3] = para.Evaluate(0);
								rn->m_maPolMagDip[z+3][z3] = para.Evaluate(1.0);
								rn->m_maPolMagDip[z+4][z3] = para.Evaluate(2.0);
							}
						}
					}
				}

				z += 2;
			}

			if (tb) {
				mprintf("    Found some outliers, repeating detection.\n");
				goto _outagain;
			}

			if (!tb2 && m_bReplaceOutliersMagnetic) {
				mprintf("    Now scanning for magnetic outliers...\n");
				tb2 = true;
				goto _outagain2;
			}
			mprintf("    Outlier correction finished.\n");
		}
	}

_ramanfrompola:

	if (m_bDumpMol1Props) {
		if (b) {
			mprintf("    Writing results to \"prop_mol1_smooth.csv\"...\n");
			a = OpenFileWrite("prop_mol1_smooth.csv",true);
		} else {
			mprintf("    Writing results to \"prop_mol1.csv\"...\n");
			a = OpenFileWrite("prop_mol1.csv",true);
		}

		mfprintf(a,"#Step");
		mfprintf(a,"; Charge; ElDipX; ElDipY; ElDipZ; ElQuadXX; ElQuadXY; ElQuadXZ; ElQuadYX; ElQuadYY; ElQuadYZ; ElQuadZX; ElQuadZY; ElQuadZZ");
		if (m_bMagMom)
			mfprintf(a,"; ElCurrX; ElCurrY; ElCurrZ; MagDipX; MagDipY; MagDipZ");
		if (m_bPola) {
			mfprintf(a,"; PolElDipXX; PolElDipXY; PolElDipXZ; PolElDipYX; PolElDipYY; PolElDipYZ; PolElDipZX; PolElDipZY; PolElDipZZ");
			mfprintf(a,"; PolElQuadXXX; PolElQuadXXY; PolElQuadXXZ; PolElQuadXYX; PolElQuadXYY; PolElQuadXYZ; PolElQuadXZX; PolElQuadXZY; PolElQuadXZZ; PolElQuadYXX; PolElQuadYXY; PolElQuadYXZ; PolElQuadYYX; PolElQuadYYY; PolElQuadYYZ; PolElQuadYZX; PolElQuadYZY; PolElQuadYZZ; PolElQuadZXX; PolElQuadZXY; PolElQuadZXZ; PolElQuadZYX; PolElQuadZYY; PolElQuadZYZ; PolElQuadZZX; PolElQuadZZY; PolElQuadZZZ");
		}
		if (m_bMagMom && m_bPola)
			mfprintf(a,"; PolMagDipXX; PolMagDipXY; PolMagDipXZ; PolMagDipYX; PolMagDipYY; PolMagDipYZ; PolMagDipZX; PolMagDipZY; PolMagDipZZ");
		mfprintf(a,"\n");

		rn = tn->m_oaMolecules[0];

		tfs = m_iStepsProcessed / 60.0;
		mprintf(WHITE,"      [");
		for (z=0;z<m_iStepsProcessed;z++) {

			if (fmod(z,tfs) < 1.0)
				mprintf(WHITE,"#");

			mfprintf(a,"%d",z+1);

			mfprintf(a,"; %.10G",rn->m_faCharge[z]);

			mfprintf(a,"; %.10G; %.10G; %.10G",
				rn->m_vaElDip[z][0],
				rn->m_vaElDip[z][1],
				rn->m_vaElDip[z][2]
			);

			mfprintf(a,"; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G",
				rn->m_maElQuad[z](0,0),rn->m_maElQuad[z](0,1),rn->m_maElQuad[z](0,2),
				rn->m_maElQuad[z](1,0),rn->m_maElQuad[z](1,1),rn->m_maElQuad[z](1,2),
				rn->m_maElQuad[z](2,0),rn->m_maElQuad[z](2,1),rn->m_maElQuad[z](2,2)
			);

			if (m_bMagMom) {

				mfprintf(a,"; %.10G; %.10G; %.10G",
					rn->m_vaElCurr[z][0],
					rn->m_vaElCurr[z][1],
					rn->m_vaElCurr[z][2]
				);

				mfprintf(a,"; %.10G; %.10G; %.10G",
					rn->m_vaMagDip[z][0],
					rn->m_vaMagDip[z][1],
					rn->m_vaMagDip[z][2]
				);
			}

			if (m_bPola) {

				mfprintf(a,"; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G",
					rn->m_maPolElDip[z](0,0),
					rn->m_maPolElDip[z](0,1),
					rn->m_maPolElDip[z](0,2),
					rn->m_maPolElDip[z](1,0),
					rn->m_maPolElDip[z](1,1),
					rn->m_maPolElDip[z](1,2),
					rn->m_maPolElDip[z](2,0),
					rn->m_maPolElDip[z](2,1),
					rn->m_maPolElDip[z](2,2)
				);

				mfprintf(a,"; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G",
					rn->m_taPolElQuad[z](0,0,0),
					rn->m_taPolElQuad[z](0,0,1),
					rn->m_taPolElQuad[z](0,0,2),
					rn->m_taPolElQuad[z](0,1,0),
					rn->m_taPolElQuad[z](0,1,1),
					rn->m_taPolElQuad[z](0,1,2),
					rn->m_taPolElQuad[z](0,2,0),
					rn->m_taPolElQuad[z](0,2,1),
					rn->m_taPolElQuad[z](0,2,2),
					rn->m_taPolElQuad[z](1,0,0),
					rn->m_taPolElQuad[z](1,0,1),
					rn->m_taPolElQuad[z](1,0,2),
					rn->m_taPolElQuad[z](1,1,0),
					rn->m_taPolElQuad[z](1,1,1),
					rn->m_taPolElQuad[z](1,1,2),
					rn->m_taPolElQuad[z](1,2,0),
					rn->m_taPolElQuad[z](1,2,1),
					rn->m_taPolElQuad[z](1,2,2),
					rn->m_taPolElQuad[z](2,0,0),
					rn->m_taPolElQuad[z](2,0,1),
					rn->m_taPolElQuad[z](2,0,2),
					rn->m_taPolElQuad[z](2,1,0),
					rn->m_taPolElQuad[z](2,1,1),
					rn->m_taPolElQuad[z](2,1,2),
					rn->m_taPolElQuad[z](2,2,0),
					rn->m_taPolElQuad[z](2,2,1),
					rn->m_taPolElQuad[z](2,2,2)
				);
			}

			if (m_bMagMom && m_bPola)
				mfprintf(a,"; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G",
					rn->m_maPolMagDip[z](0,0),
					rn->m_maPolMagDip[z](0,1),
					rn->m_maPolMagDip[z](0,2),
					rn->m_maPolMagDip[z](1,0),
					rn->m_maPolMagDip[z](1,1),
					rn->m_maPolMagDip[z](1,2),
					rn->m_maPolMagDip[z](2,0),
					rn->m_maPolMagDip[z](2,1),
					rn->m_maPolMagDip[z](2,2)
				);

			mfprintf(a,"\n");
			//fflush(a);
		}
		mprintf(WHITE,"]\n");

		fclose(a);
	}


	mprintf("    Deriving all properties...\n");

	if (m_iStepsProcessed < 3) {
		eprintf("Error: Less than 3 frames processed, cannot compute derivative.\n\n");
		abort();
	}

	tfs = g_oaSingleMolecules.GetSize() / 60.0;
	if (g_oaSingleMolecules.GetSize() > 1)
		mprintf(WHITE,"      [");
	for (z=0;z<g_oaSingleMolecules.GetSize();z++) {

		if (g_oaSingleMolecules.GetSize() > 1)
			if (fmod(z,tfs) < 1.0)
				mprintf(WHITE,"#");

		rn = tn->m_oaMolecules[z];

		rn->m_faDCharge.resize(m_iStepsProcessed-2);
		rn->m_vaDElDip.resize(m_iStepsProcessed-2);
		rn->m_maDElQuad.resize(m_iStepsProcessed-2);
		if (m_bMagMom) {
			rn->m_vaDElCurr.resize(m_iStepsProcessed-2);
			rn->m_vaDMagDip.resize(m_iStepsProcessed-2);
		}
		if (m_bPola) {
			rn->m_maDPolElDip.resize(m_iStepsProcessed-2);
			rn->m_taDPolElQuad.resize(m_iStepsProcessed-2,CDTensor3(3,3,3));
		}
		if (m_bMagMom && m_bPola)
			rn->m_maDPolMagDip.resize(m_iStepsProcessed-2);

		for (z2=0;z2<m_iStepsProcessed-2;z2++) {

			// Unit: e/fs
			rn->m_faDCharge[z2] = (rn->m_faCharge[z2+2] - rn->m_faCharge[z2]) / (2.0 * g_fTimestepLength * g_iStride);

			// Unit: Debye/fs
			for (z3=0;z3<3;z3++)
				rn->m_vaDElDip[z2][z3] = (rn->m_vaElDip[z2+2][z3] - rn->m_vaElDip[z2][z3]) / (2.0 * g_fTimestepLength * g_iStride);

			// Unit: Debye*pm/fs
			for (z3=0;z3<9;z3++)
				rn->m_maDElQuad[z2][z3] = (rn->m_maElQuad[z2+2][z3] - rn->m_maElQuad[z2][z3]) / (2.0 * g_fTimestepLength * g_iStride);

			if (m_bMagMom) {

				// Unit: MB/(pm*fs)
				for (z3=0;z3<3;z3++)
					rn->m_vaDElCurr[z2][z3] = (rn->m_vaElCurr[z2+2][z3] - rn->m_vaElCurr[z2][z3]) / (2.0 * g_fTimestepLength * g_iStride);

				// Unit: MB/fs
				for (z3=0;z3<3;z3++)
					rn->m_vaDMagDip[z2][z3] = (rn->m_vaMagDip[z2+2][z3] - rn->m_vaMagDip[z2][z3]) / (2.0 * g_fTimestepLength * g_iStride);
			}

			if (m_bPola) {

				// Unit: e*pm^2/(V*fs)
				for (z3=0;z3<9;z3++)
					rn->m_maDPolElDip[z2][z3] = (rn->m_maPolElDip[z2+2][z3] - rn->m_maPolElDip[z2][z3]) / (2.0 * g_fTimestepLength * g_iStride);

				// Unit: e*pm^3/(V*fs)
				for (z3=0;z3<27;z3++)
					rn->m_taDPolElQuad[z2][z3] = (rn->m_taPolElQuad[z2+2][z3] - rn->m_taPolElQuad[z2][z3]) / (2.0 * g_fTimestepLength * g_iStride);
			}

			if (m_bMagMom && m_bPola) {

				// Unit: e*pm^3/(V*fs^2)
				for (z3=0;z3<9;z3++)
					rn->m_maDPolMagDip[z2][z3] = (rn->m_maPolMagDip[z2+2][z3] - rn->m_maPolMagDip[z2][z3]) / (2.0 * g_fTimestepLength * g_iStride);
			}
		}
	}

	if (g_oaSingleMolecules.GetSize() > 1)
		mprintf(WHITE,"]\n");

	if (m_bDumpMol1Props) {

		if (b) {
			mprintf("    Writing results to \"prop_mol1_deriv_smooth.csv\"...\n");
			a = OpenFileWrite("prop_mol1_deriv_smooth.csv",true);
		} else {
			mprintf("    Writing results to \"prop_mol1_deriv.csv\"...\n");
			a = OpenFileWrite("prop_mol1_deriv.csv",true);
		}

		mfprintf(a,"#Step");
		mfprintf(a,"; DCharge; DElDipX; DElDipY; DElDipZ; DElQuadXX; DElQuadXY; DElQuadXZ; DElQuadYX; DElQuadYY; DElQuadYZ; DElQuadZX; DElQuadZY; DElQuadZZ");
		if (m_bMagMom)
			mfprintf(a,"; DElCurrX; DElCurrY; DElCurrZ; DMagDipX; DMagDipY; DMagDipZ");
		if (m_bPola) {
			mfprintf(a,"; DPolElDipXX; DPolElDipXY; DPolElDipXZ; DPolElDipYX; DPolElDipYY; DPolElDipYZ; DPolElDipZX; DPolElDipZY; DPolElDipZZ");
			mfprintf(a,"; DPolElQuadXXX; DPolElQuadXXY; DPolElQuadXXZ; DPolElQuadXYX; DPolElQuadXYY; DPolElQuadXYZ; DPolElQuadXZX; DPolElQuadXZY; DPolElQuadXZZ; DPolElQuadYXX; DPolElQuadYXY; DPolElQuadYXZ; DPolElQuadYYX; DPolElQuadYYY; DPolElQuadYYZ; DPolElQuadYZX; DPolElQuadYZY; DPolElQuadYZZ; DPolElQuadZXX; DPolElQuadZXY; DPolElQuadZXZ; DPolElQuadZYX; DPolElQuadZYY; DPolElQuadZYZ; DPolElQuadZZX; DPolElQuadZZY; DPolElQuadZZZ");
		}
		if (m_bMagMom && m_bPola)
			mfprintf(a,"; DPolMagDipXX; DPolMagDipXY; DPolMagDipXZ; DPolMagDipYX; DPolMagDipYY; DPolMagDipYZ; DPolMagDipZX; DPolMagDipZY; DPolMagDipZZ");
		mfprintf(a,"\n");

		rn = tn->m_oaMolecules[0];

		tfs = (m_iStepsProcessed-2) / 60.0;
		mprintf(WHITE,"      [");
		for (z=0;z<m_iStepsProcessed-2;z++) {

			if (fmod(z,tfs) < 1.0)
				mprintf(WHITE,"#");

			mfprintf(a,"%d",z+2);

			mfprintf(a,"; %.10G",rn->m_faDCharge[z]);

			mfprintf(a,"; %.10G; %.10G; %.10G",
				rn->m_vaDElDip[z][0],
				rn->m_vaDElDip[z][1],
				rn->m_vaDElDip[z][2]
			);

			mfprintf(a,"; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G",
				rn->m_maDElQuad[z](0,0),
				rn->m_maDElQuad[z](0,1),
				rn->m_maDElQuad[z](0,2),
				rn->m_maDElQuad[z](1,0),
				rn->m_maDElQuad[z](1,1),
				rn->m_maDElQuad[z](1,2),
				rn->m_maDElQuad[z](2,0),
				rn->m_maDElQuad[z](2,1),
				rn->m_maDElQuad[z](2,2)
			);

			if (m_bMagMom) {

				mfprintf(a,"; %.10G; %.10G; %.10G",
					rn->m_vaDElCurr[z][0],
					rn->m_vaDElCurr[z][1],
					rn->m_vaDElCurr[z][2]
				);

				mfprintf(a,"; %.10G; %.10G; %.10G",
					rn->m_vaDMagDip[z][0],
					rn->m_vaDMagDip[z][1],
					rn->m_vaDMagDip[z][2]
				);
			}

			if (m_bPola) {

				mfprintf(a,"; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G",
					rn->m_maDPolElDip[z](0,0),
					rn->m_maDPolElDip[z](0,1),
					rn->m_maDPolElDip[z](0,2),
					rn->m_maDPolElDip[z](1,0),
					rn->m_maDPolElDip[z](1,1),
					rn->m_maDPolElDip[z](1,2),
					rn->m_maDPolElDip[z](2,0),
					rn->m_maDPolElDip[z](2,1),
					rn->m_maDPolElDip[z](2,2)
				);

				mfprintf(a,"; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G",
					rn->m_taDPolElQuad[z](0,0,0),
					rn->m_taDPolElQuad[z](0,0,1),
					rn->m_taDPolElQuad[z](0,0,2),
					rn->m_taDPolElQuad[z](0,1,0),
					rn->m_taDPolElQuad[z](0,1,1),
					rn->m_taDPolElQuad[z](0,1,2),
					rn->m_taDPolElQuad[z](0,2,0),
					rn->m_taDPolElQuad[z](0,2,1),
					rn->m_taDPolElQuad[z](0,2,2),
					rn->m_taDPolElQuad[z](1,0,0),
					rn->m_taDPolElQuad[z](1,0,1),
					rn->m_taDPolElQuad[z](1,0,2),
					rn->m_taDPolElQuad[z](1,1,0),
					rn->m_taDPolElQuad[z](1,1,1),
					rn->m_taDPolElQuad[z](1,1,2),
					rn->m_taDPolElQuad[z](1,2,0),
					rn->m_taDPolElQuad[z](1,2,1),
					rn->m_taDPolElQuad[z](1,2,2),
					rn->m_taDPolElQuad[z](2,0,0),
					rn->m_taDPolElQuad[z](2,0,1),
					rn->m_taDPolElQuad[z](2,0,2),
					rn->m_taDPolElQuad[z](2,1,0),
					rn->m_taDPolElQuad[z](2,1,1),
					rn->m_taDPolElQuad[z](2,1,2),
					rn->m_taDPolElQuad[z](2,2,0),
					rn->m_taDPolElQuad[z](2,2,1),
					rn->m_taDPolElQuad[z](2,2,2)
				);
			}

			if (m_bMagMom && m_bPola) {

				mfprintf(a,"; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G; %.10G",
					rn->m_maDPolMagDip[z](0,0),
					rn->m_maDPolMagDip[z](0,1),
					rn->m_maDPolMagDip[z](0,2),
					rn->m_maDPolMagDip[z](1,0),
					rn->m_maDPolMagDip[z](1,1),
					rn->m_maDPolMagDip[z](1,2),
					rn->m_maDPolMagDip[z](2,0),
					rn->m_maDPolMagDip[z](2,1),
					rn->m_maDPolMagDip[z](2,2)
				);
			}

			mfprintf(a,"\n");
			//fflush(a);
		}
		mprintf(WHITE,"]\n");

		fclose(a);
	}

	if (m_bWriteACFs) {

		b2 = false;

_normalized:
		buf.sprintf("acfs");
		if (b)
			buf.strcat("_smooth");
		if (b2)
			buf.strcat("_normalized");
		buf.strcat(".csv");
		mprintf("    Writing autocorrelation functions to \"%s\"...\n",(const char*)buf);
		a = OpenFileWrite((const char*)buf,true);

		ac = new CAutoCorrelation();
		ac->Init(m_iStepsProcessed-2,m_oaObservations[0]->m_iDepth,g_bACFFFT);
		dain.SetSize(m_iStepsProcessed-2);
		daout.SetSize(m_oaObservations[0]->m_iDepth);

		fft = new CFFT();
		fft->PrepareFFT_C2C(2*(m_oaObservations[0]->m_iDepth+m_oaObservations[0]->m_iZeroPadding));
		temp.resize(2*(m_oaObservations[0]->m_iDepth+m_oaObservations[0]->m_iZeroPadding));

		buf2 = "";

		tfs = 64 / 60.0;
		mprintf(WHITE,"      [");
		for (z=0;z<64;z++) {

			if (fmod(z,tfs) < 1.0)
				mprintf(WHITE,"#");

			if (!m_bMagMom)
				if (((z >= 13) && (z <= 18)) || ((z >= 55) && (z <= 63)))
					continue;

			if (!m_bPola)
				if ((z >= 19) && (z <= 63))
					continue;

			daacf[z].SetSize(m_oaObservations[0]->m_iDepth);
			for (z2=0;z2<m_oaObservations[0]->m_iDepth;z2++)
				daacf[z][z2] = 0;

			for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++) {

				rn = tn->m_oaMolecules[z2];

				switch(z) {
					case  0: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_faDCharge[z3];        if (z2 == 0) buf2 += "; DCharge"; break;

					case  1: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_vaDElDip[z3][0];      if (z2 == 0) buf2 += "; DElDipX"; break;
					case  2: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_vaDElDip[z3][1];      if (z2 == 0) buf2 += "; DElDipY"; break;
					case  3: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_vaDElDip[z3][2];      if (z2 == 0) buf2 += "; DElDipZ"; break;

					case  4: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDElQuad[z3][0];     if (z2 == 0) buf2 += "; DElQuadXX"; break;
					case  5: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDElQuad[z3][1];     if (z2 == 0) buf2 += "; DElQuadXY"; break;
					case  6: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDElQuad[z3][2];     if (z2 == 0) buf2 += "; DElQuadXZ"; break;
					case  7: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDElQuad[z3][3];     if (z2 == 0) buf2 += "; DElQuadYX"; break;
					case  8: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDElQuad[z3][4];     if (z2 == 0) buf2 += "; DElQuadYY"; break;
					case  9: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDElQuad[z3][5];     if (z2 == 0) buf2 += "; DElQuadYZ"; break;
					case 10: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDElQuad[z3][6];     if (z2 == 0) buf2 += "; DElQuadZX"; break;
					case 11: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDElQuad[z3][7];     if (z2 == 0) buf2 += "; DElQuadZY"; break;
					case 12: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDElQuad[z3][8];     if (z2 == 0) buf2 += "; DElQuadZZ"; break;

					case 13: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_vaDElCurr[z3][0];     if (z2 == 0) buf2 += "; DElCurrX"; break;
					case 14: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_vaDElCurr[z3][1];     if (z2 == 0) buf2 += "; DElCurrY"; break;
					case 15: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_vaDElCurr[z3][2];     if (z2 == 0) buf2 += "; DElCurrZ"; break;

					case 16: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_vaDMagDip[z3][0];     if (z2 == 0) buf2 += "; DMagDipX"; break;
					case 17: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_vaDMagDip[z3][1];     if (z2 == 0) buf2 += "; DMagDipY"; break;
					case 18: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_vaDMagDip[z3][2];     if (z2 == 0) buf2 += "; DMagDipZ"; break;

					case 19: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolElDip[z3][0];   if (z2 == 0) buf2 += "; DPolElDipXX"; break;
					case 20: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolElDip[z3][1];   if (z2 == 0) buf2 += "; DPolElDipXY"; break;
					case 21: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolElDip[z3][2];   if (z2 == 0) buf2 += "; DPolElDipXZ"; break;
					case 22: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolElDip[z3][3];   if (z2 == 0) buf2 += "; DPolElDipYX"; break;
					case 23: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolElDip[z3][4];   if (z2 == 0) buf2 += "; DPolElDipYY"; break;
					case 24: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolElDip[z3][5];   if (z2 == 0) buf2 += "; DPolElDipYZ"; break;
					case 25: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolElDip[z3][6];   if (z2 == 0) buf2 += "; DPolElDipZX"; break;
					case 26: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolElDip[z3][7];   if (z2 == 0) buf2 += "; DPolElDipZY"; break;
					case 27: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolElDip[z3][8];   if (z2 == 0) buf2 += "; DPolElDipZZ"; break;

					case 28: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][0];  if (z2 == 0) buf2 += "; DPolElQuadXXX"; break;
					case 29: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][1];  if (z2 == 0) buf2 += "; DPolElQuadXXY"; break;
					case 30: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][2];  if (z2 == 0) buf2 += "; DPolElQuadXXZ"; break;
					case 31: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][3];  if (z2 == 0) buf2 += "; DPolElQuadXYX"; break;
					case 32: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][4];  if (z2 == 0) buf2 += "; DPolElQuadXYY"; break;
					case 33: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][5];  if (z2 == 0) buf2 += "; DPolElQuadXYZ"; break;
					case 34: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][6];  if (z2 == 0) buf2 += "; DPolElQuadXZX"; break;
					case 35: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][7];  if (z2 == 0) buf2 += "; DPolElQuadXZY"; break;
					case 36: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][8];  if (z2 == 0) buf2 += "; DPolElQuadXZZ"; break;
					case 37: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][9];  if (z2 == 0) buf2 += "; DPolElQuadYXX"; break;
					case 38: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][10]; if (z2 == 0) buf2 += "; DPolElQuadYXY"; break;
					case 39: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][11]; if (z2 == 0) buf2 += "; DPolElQuadYXZ"; break;
					case 40: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][12]; if (z2 == 0) buf2 += "; DPolElQuadYYX"; break;
					case 41: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][13]; if (z2 == 0) buf2 += "; DPolElQuadYYY"; break;
					case 42: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][14]; if (z2 == 0) buf2 += "; DPolElQuadYYZ"; break;
					case 43: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][15]; if (z2 == 0) buf2 += "; DPolElQuadYZX"; break;
					case 44: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][16]; if (z2 == 0) buf2 += "; DPolElQuadYZY"; break;
					case 45: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][17]; if (z2 == 0) buf2 += "; DPolElQuadYZZ"; break;
					case 46: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][18]; if (z2 == 0) buf2 += "; DPolElQuadZXX"; break;
					case 47: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][19]; if (z2 == 0) buf2 += "; DPolElQuadZXY"; break;
					case 48: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][20]; if (z2 == 0) buf2 += "; DPolElQuadZXZ"; break;
					case 49: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][21]; if (z2 == 0) buf2 += "; DPolElQuadZYX"; break;
					case 50: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][22]; if (z2 == 0) buf2 += "; DPolElQuadZYY"; break;
					case 51: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][23]; if (z2 == 0) buf2 += "; DPolElQuadZYZ"; break;
					case 52: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][24]; if (z2 == 0) buf2 += "; DPolElQuadZZX"; break;
					case 53: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][25]; if (z2 == 0) buf2 += "; DPolElQuadZZY"; break;
					case 54: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_taDPolElQuad[z3][26]; if (z2 == 0) buf2 += "; DPolElQuadZZZ"; break;

					case 55: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolMagDip[z3][0];  if (z2 == 0) buf2 += "; DPolMagDipXX"; break;
					case 56: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolMagDip[z3][1];  if (z2 == 0) buf2 += "; DPolMagDipXY"; break;
					case 57: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolMagDip[z3][2];  if (z2 == 0) buf2 += "; DPolMagDipXZ"; break;
					case 58: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolMagDip[z3][3];  if (z2 == 0) buf2 += "; DPolMagDipYX"; break;
					case 59: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolMagDip[z3][4];  if (z2 == 0) buf2 += "; DPolMagDipYY"; break;
					case 60: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolMagDip[z3][5];  if (z2 == 0) buf2 += "; DPolMagDipYZ"; break;
					case 61: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolMagDip[z3][6];  if (z2 == 0) buf2 += "; DPolMagDipZX"; break;
					case 62: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolMagDip[z3][7];  if (z2 == 0) buf2 += "; DPolMagDipZY"; break;
					case 63: for (z3=0;z3<m_iStepsProcessed-2;z3++) dain[z3] = rn->m_maDPolMagDip[z3][8];  if (z2 == 0) buf2 += "; DPolMagDipZZ"; break;
				}

				ac->AutoCorrelate(&dain,&daout);
				for (z3=0;z3<m_oaObservations[0]->m_iDepth;z3++)
					daacf[z][z3] += daout[z3];
			}

			for (z3=0;z3<m_oaObservations[0]->m_iDepth;z3++)
				daacf[z][z3] /= (double)g_oaSingleMolecules.GetSize();

			if (b2) {
				f = daacf[z][0];
				for (z3=0;z3<m_oaObservations[0]->m_iDepth;z3++)
					daacf[z][z3] /= f;
			}

			for (z3=0;z3<m_oaObservations[0]->m_iDepth;z3++)
				temp[z3] = daacf[z][z3];

			for (z3=m_oaObservations[0]->m_iDepth;z3<m_oaObservations[0]->m_iDepth+m_oaObservations[0]->m_iZeroPadding;z3++)
				temp[z3] = 0.0;

			for (z3=0;z3<m_oaObservations[0]->m_iDepth;z3++)
				temp[z3] *= pow2(cos((double)z3 / (m_oaObservations[0]->m_iDepth - 1) / 2.0 * Pi));

			o = m_oaObservations[0]->m_iDepth+m_oaObservations[0]->m_iZeroPadding;
			for (z3=1;z3<o;z3++)
				temp[o+z3] = temp[o-z3];
			temp[o] = 0.0;

			for (z3=0;z3<2*(m_oaObservations[0]->m_iDepth+m_oaObservations[0]->m_iZeroPadding);z3++) {
				fft->m_pInput[2*z3] = temp[z3];
				fft->m_pInput[2*z3+1] = 0.0;
			}
			fft->DoFFT();
	
			daspec[z].SetSize(m_oaObservations[0]->m_iSpecLength);

			for (z3=0;z3<m_oaObservations[0]->m_iSpecLength;z3++)
				daspec[z][z3] = fft->m_pOutput[2*z3];

			f = m_oaObservations[0]->m_fSpecResolution * g_fTimestepLength * g_iStride * 1.883652e-4;
			for (z3=1;z3<m_oaObservations[0]->m_iSpecLength;z3++)
				daspec[z][z3] *= pow2(f * z3 / sin(f * z3)); // Divide by sinc function to correct finite difference derivation
		}
		delete ac;
		delete fft;

		mprintf(WHITE,"]\n");

		mfprintf(a,"#Depth/fs%s\n",(const char*)buf2);

		for (z=0;z<m_oaObservations[0]->m_iDepth;z++) {
			mfprintf(a,"%.2f",z*g_fTimestepLength*g_iStride);
			for (z2=0;z2<64;z2++) {
				if (!m_bMagMom)
					if (((z2 >= 13) && (z2 <= 18)) || ((z2 >= 55) && (z2 <= 63)))
						continue;
				if (!m_bPola)
					if ((z2 >= 19) && (z2 <= 63))
						continue;
				mfprintf(a,"; %.10G",daacf[z2][z]);
			}
			mfprintf(a,"\n");
		}

		fclose(a);

		buf.sprintf("acf_spectra");
		if (b)
			buf.strcat("_smooth");
		if (b2)
			buf.strcat("_normalized");
		buf.strcat(".csv");
		mprintf("    Writing autocorrelation spectra to \"%s\"...\n",(const char*)buf);
		a = OpenFileWrite((const char*)buf,true);
		mfprintf(a,"#Wavenumber%s\n",(const char*)buf2);

		for (z=0;z<m_oaObservations[0]->m_iSpecLength;z++) {
			mfprintf(a,"%.2f",m_oaObservations[0]->m_fSpecResolution*z);
			for (z2=0;z2<64;z2++) {
				if (!m_bMagMom)
					if (((z2 >= 13) && (z2 <= 18)) || ((z2 >= 55) && (z2 <= 63)))
						continue;
				if (!m_bPola)
					if ((z2 >= 19) && (z2 <= 63))
						continue;
				mfprintf(a,"; %.10G",daspec[z2][z]);
			}
			mfprintf(a,"\n");
		}

		fclose(a);

		if (!b2) {
			b2 = true;
			goto _normalized;
		}
	}

	if (m_bSmoothData && !b) {
		mprintf("\n    Smoothing properties...\n\n");
		b = true;
		smo.Init(m_iStepsProcessed,m_fSmoothWaveNumber);
		for (z=0;z<(int)m_oaTrajectories.size();z++) {
			for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++) {
				smo.Smoothen(m_oaTrajectories[z]->m_oaMolecules[z2]->m_faCharge);
				smo.Smoothen(m_oaTrajectories[z]->m_oaMolecules[z2]->m_vaElDip);
				smo.Smoothen(m_oaTrajectories[z]->m_oaMolecules[z2]->m_maElQuad);
				if (m_bMagMom) {
					smo.Smoothen(m_oaTrajectories[z]->m_oaMolecules[z2]->m_vaElCurr);
					smo.Smoothen(m_oaTrajectories[z]->m_oaMolecules[z2]->m_vaMagDip);
				}
			}
		}
		goto _again;
	}

	for (z=0;z<(int)m_oaObservations.size();z++) {
		mprintf(WHITE,"\n    * Finishing Observation %d: %s\n",
			z+1,(const char*)m_oaObservations[z]->m_sName);
		if (m_iStepsProcessed-2 <= m_oaObservations[z]->m_iDepth) {
			mprintf(RED,"      Warning: ");
			mprintf("The ACF depth (%d) is larger than the number of data points (%d).\n",
				m_oaObservations[z]->m_iDepth,m_iStepsProcessed-2);
			mprintf("               The results of this observation will be meaningless.\n");
		}
		mprintf("      Computing ACF...\n");
		ComputeACF(m_oaObservations[z]);
		mprintf("      Computing Spectrum...\n");
		ComputeSpectrum(m_oaObservations[z]);
	}

	mprintf(WHITE,"\n    <<< Finished Analysis <<<\n");
}




void CROAEngine::ComputeACFPair(CROAObservation *obs, CROAMolecule *mol1, CROAMolecule *mol2, std::vector<double> *wfn) {

	int zc, z, zi, zj, zk;
	//FILE *a;


	if (wfn != NULL) {

		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = (*wfn)[z];
			m_faInput2[z] = (*wfn)[z];
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faWFN[z] += m_faOutput[z];
	}

	/*********************************************************************************************************************/
	if (obs->m_iType == ROA_SPECTRUM_IR) {

		for (zc=0;zc<3;zc++) {
			if (wfn != NULL) {
				for (z=0;z<m_iStepsProcessed-2;z++) {
					m_faInput1[z] = mol1->m_vaDElDip[z][zc] * (*wfn)[z];
					m_faInput2[z] = mol2->m_vaDElDip[z][zc] * (*wfn)[z];
				}
			} else {
				for (z=0;z<m_iStepsProcessed-2;z++) {
					m_faInput1[z] = mol1->m_vaDElDip[z][zc];
					m_faInput2[z] = mol2->m_vaDElDip[z][zc];
				}
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++)
				obs->m_faACF[z] += m_faOutput[z];
		}

	/*********************************************************************************************************************/
	} else if (obs->m_iType == ROA_SPECTRUM_RAMAN) {

		// Isotropic Part
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = (mol1->m_maDPolElDip[z](0,0) + mol1->m_maDPolElDip[z](1,1) + mol1->m_maDPolElDip[z](2,2)) / 3.0;
			m_faInput2[z] = (mol2->m_maDPolElDip[z](0,0) + mol2->m_maDPolElDip[z](1,1) + mol2->m_maDPolElDip[z](2,2)) / 3.0;
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF[z] += m_faOutput[z];
		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF_CTN[z] += 0.5 * m_faOutput[z];
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF_CTN[z] -= 0.5 * m_faOutput[z];
			}
		}

		// First Anisotropic Part
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](0,0) - mol1->m_maDPolElDip[z](1,1);
			m_faInput2[z] = mol2->m_maDPolElDip[z](0,0) - mol2->m_maDPolElDip[z](1,1);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] / 2.0;
		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z];
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z];
			}
		}

		// Second Anisotropic Part
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](1,1) - mol1->m_maDPolElDip[z](2,2);
			m_faInput2[z] = mol2->m_maDPolElDip[z](1,1) - mol2->m_maDPolElDip[z](2,2);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] / 2.0;
		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z];
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z];
			}
		}

		// Third Anisotropic Part
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](2,2) - mol1->m_maDPolElDip[z](0,0);
			m_faInput2[z] = mol2->m_maDPolElDip[z](2,2) - mol2->m_maDPolElDip[z](0,0);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] / 2.0;
		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z];
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z];
			}
		}

		// Fourth Anisotropic Part
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](0,1);
			m_faInput2[z] = mol2->m_maDPolElDip[z](0,1);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] * 3.0;
		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z];
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z];
			}
		}

		// Fifth Anisotropic Part
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](1,2);
			m_faInput2[z] = mol2->m_maDPolElDip[z](1,2);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] * 3.0;
		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z];
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z];
			}
		}

		// Sixth Anisotropic Part
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](2,0);
			m_faInput2[z] = mol2->m_maDPolElDip[z](2,0);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] * 3.0;
		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z];
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z];
			}
		}

	/*********************************************************************************************************************/
	} else if (obs->m_iType == ROA_SPECTRUM_VCD) {

		for (zc=0;zc<3;zc++) {
			for (z=0;z<m_iStepsProcessed-2;z++) {
				m_faInput1[z] = mol1->m_vaDElDip[z][zc];
				m_faInput2[z] = mol2->m_vaDMagDip[z][zc];
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++)
				obs->m_faACF[z] += m_faOutput[z];

			if (obs->m_bUseCommutator) {
				for (z=0;z<obs->m_iDepth;z++) {
					obs->m_faACF_CTP[z] += 0.5 * m_faOutput[z];
					obs->m_faACF_CTN[z] += 0.5 * m_faOutput[z];
				}
				m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
				for (z=0;z<obs->m_iDepth;z++) {
					obs->m_faACF_CTP[z] += 0.5 * m_faOutput[z];
					obs->m_faACF_CTN[z] -= 0.5 * m_faOutput[z];
				}
			}
		}

	/*********************************************************************************************************************/
	} else if (obs->m_iType == ROA_SPECTRUM_SFG) {

		for (zi=0;zi<3;zi++){
			for (zj=0;zj<3;zj++){
				for (zk=0;zk<3;zk++){
					for (z=0;z<m_iStepsProcessed-2;z++) {
						m_faInput1[z] = mol1->m_maDPolElDip[z](zi,zj);
						m_faInput2[z] = mol2->m_vaDElDip[z][zk];
					}
					m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
					for (z=0;z<obs->m_iDepth;z++)
						obs->m_faACF[z] += m_faOutput[z];
					if (obs->m_bSingleACFs)
						for (z=0;z<obs->m_iDepth;z++)
							obs->m_faaSingleACFs[zi*9+zj*3+zk][z] += m_faOutput[z];
					if (obs->m_bUseCommutator) {
						for (z=0;z<obs->m_iDepth;z++) {
							obs->m_faACF_CTP[z] += 0.5 * m_faOutput[z];
							obs->m_faACF_CTN[z] += 0.5 * m_faOutput[z];
						}
						if (obs->m_bSingleACFs) {
							for (z=0;z<obs->m_iDepth;z++) {
								obs->m_faaSingleACFs_CTP[zi*9+zj*3+zk][z] += 0.5 * m_faOutput[z];
								obs->m_faaSingleACFs_CTN[zi*9+zj*3+zk][z] += 0.5 * m_faOutput[z];
							}
						}
						m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
						for (z=0;z<obs->m_iDepth;z++) {
							obs->m_faACF_CTP[z] += 0.5 * m_faOutput[z];
							obs->m_faACF_CTN[z] -= 0.5 * m_faOutput[z];
						}
						if (obs->m_bSingleACFs) {
							for (z=0;z<obs->m_iDepth;z++) {
								obs->m_faaSingleACFs_CTP[zi*9+zj*3+zk][z] += 0.5 * m_faOutput[z];
								obs->m_faaSingleACFs_CTN[zi*9+zj*3+zk][z] -= 0.5 * m_faOutput[z];
							}
						}
					}
				}
			}
		}

	/*********************************************************************************************************************/
	} else if (obs->m_iType == ROA_SPECTRUM_ROA) {

		// aG
		for (z=0;z<m_iStepsProcessed-2;z++) {
			// Unit: e*pm^2/(V*fs)
			m_faInput1[z] = (mol1->m_maDPolElDip[z](0,0) + mol1->m_maDPolElDip[z](1,1) + mol1->m_maDPolElDip[z](2,2)) / 3.0;
			// Unit: e*pm^3/(V*fs^2)
			m_faInput2[z] = -(mol2->m_maDPolMagDip[z](0,0) + mol2->m_maDPolMagDip[z](1,1) + mol2->m_maDPolMagDip[z](1,1)) / 3.0;
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF[z] += m_faOutput[z]; // Unit: e^2*pm^5/(V^2*fs^3)

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF_CTN[z] += 0.5 * m_faOutput[z];
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF_CTP[z] += 0.5 * m_faOutput[z];
				obs->m_faACF_CTN[z] -= 0.5 * m_faOutput[z];
			}
		}


		// GammaG - Part 1
		for (z=0;z<m_iStepsProcessed-2;z++) {
			// Unit: e*pm^2/(V*fs)
			m_faInput1[z] = mol1->m_maDPolElDip[z](0,0) - mol1->m_maDPolElDip[z](1,1);
			// Unit: e*pm^3/(V*fs^2)
			m_faInput2[z] = -(mol2->m_maDPolMagDip[z](0,0) - mol2->m_maDPolMagDip[z](1,1));
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] / 2.0; // Unit: e^2*pm^5/(V^2*fs^3)

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] / 2.0;
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z] / 2.0;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] / 2.0;
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z] / 2.0;
			}
		}

		// GammaG - Part 2
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](1,1) - mol1->m_maDPolElDip[z](2,2);
			m_faInput2[z] = -(mol2->m_maDPolMagDip[z](1,1) - mol2->m_maDPolMagDip[z](2,2));
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] / 2.0;

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] / 2.0;
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z] / 2.0;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] / 2.0;
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z] / 2.0;
			}
		}


		// GammaG - Part 3
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](2,2) - mol1->m_maDPolElDip[z](0,0);
			m_faInput2[z] = -(mol2->m_maDPolMagDip[z](2,2) - mol2->m_maDPolMagDip[z](0,0));
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] / 2.0;

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] / 2.0;
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z] / 2.0;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] / 2.0;
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z] / 2.0;
			}
		}

		// GammaG - Part 4
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](0,1);
			m_faInput2[z] = -(mol2->m_maDPolMagDip[z](0,1) + mol2->m_maDPolMagDip[z](1,0));
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] * 1.5;

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] * 1.5;
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z] * 1.5;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] * 1.5;
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z] * 1.5;
			}
		}

		// GammaG - Part 5
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](1,2);
			m_faInput2[z] = -(mol2->m_maDPolMagDip[z](1,2) + mol2->m_maDPolMagDip[z](2,1));
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] * 1.5;

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] * 1.5;
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z] * 1.5;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] * 1.5;
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z] * 1.5;
			}
		}

		// GammaG - Part 6
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](2,0);
			m_faInput2[z] = -(mol2->m_maDPolMagDip[z](2,0) + mol2->m_maDPolMagDip[z](0,2));
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF2[z] += m_faOutput[z] * 1.5;

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] * 1.5;
				obs->m_faACF2_CTN[z] += 0.5 * m_faOutput[z] * 1.5;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF2_CTP[z] += 0.5 * m_faOutput[z] * 1.5;
				obs->m_faACF2_CTN[z] -= 0.5 * m_faOutput[z] * 1.5;
			}
		}


		// GammaA - Part 1
		for (z=0;z<m_iStepsProcessed-2;z++) {
			// Unit: e*pm^2/(V*fs)
			m_faInput1[z] = mol1->m_maDPolElDip[z](1,1) - mol1->m_maDPolElDip[z](0,0);
			// Unit: e*pm^3/(V*fs)
			m_faInput2[z] = mol2->m_taDPolElQuad[z](2,0,1);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			// Unit: e^2*pm^5/(V^2*fs^3)
			obs->m_faACF3[z] += m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] -= 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
		}

		// GammaA - Part 2
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](0,0) - mol1->m_maDPolElDip[z](2,2);
			m_faInput2[z] = mol2->m_taDPolElQuad[z](1,2,0);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF3[z] += m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] -= 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
		}

		// GammaA - Part 3
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](2,2) - mol1->m_maDPolElDip[z](1,1);
			m_faInput2[z] = mol2->m_taDPolElQuad[z](0,1,2);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF3[z] += m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] -= 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
		}

		// GammaA - Part 4
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](0,1);
			m_faInput2[z] = mol2->m_taDPolElQuad[z](1,1,2) - mol2->m_taDPolElQuad[z](2,1,1) + mol2->m_taDPolElQuad[z](2,0,0) - mol2->m_taDPolElQuad[z](0,0,2);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF3[z] += m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] -= 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
		}

		// GammaA - Part 5
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](1,2);
			m_faInput2[z] = mol2->m_taDPolElQuad[z](2,2,0) - mol2->m_taDPolElQuad[z](0,2,2) + mol2->m_taDPolElQuad[z](0,1,1) - mol2->m_taDPolElQuad[z](1,1,0);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF3[z] += m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] -= 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
		}

		// GammaA - Part 6
		for (z=0;z<m_iStepsProcessed-2;z++) {
			m_faInput1[z] = mol1->m_maDPolElDip[z](2,0);
			m_faInput2[z] = mol2->m_taDPolElQuad[z](1,2,2) - mol2->m_taDPolElQuad[z](2,2,1) + mol2->m_taDPolElQuad[z](0,0,1) - mol2->m_taDPolElQuad[z](1,0,0);
		}
		m_pCrossCorr->CrossCorrelate(&m_faInput1,&m_faInput2,&m_faOutput);
		for (z=0;z<obs->m_iDepth;z++)
			obs->m_faACF3[z] += m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
			m_pCrossCorr->CrossCorrelate(&m_faInput2,&m_faInput1,&m_faOutput);
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF3_CTP[z] += 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
				obs->m_faACF3_CTN[z] -= 0.5 * m_faOutput[z] / 2.0 * 1.0E-15 * obs->m_fLaserFreq * 100.0 * CONST_SPEED_OF_LIGHT * 2.0 * Pi;
			}
		}

	/*********************************************************************************************************************/
	} else {
		eprintf("CROAEngine::ComputeACFPair(): Error: Unknown spectrum type %d.\n",obs->m_iType);
		abort();
	}
}


//int g_iTC=0;




void CROAEngine::ComputeACF(CROAObservation *obs) {

	int z, z2, z3, z4, imol, imol2, acc;
	double pc, tfs;
	CROATrajectory *tr;
	CxString buf;
	std::vector<double> wfn;



	try { m_pCrossCorr = new CCrossCorrelation(); } catch(...) { m_pCrossCorr = NULL; }
	if (m_pCrossCorr == NULL) NewException((double)sizeof(CCrossCorrelation), __FILE__, __LINE__, __PRETTY_FUNCTION__);

	m_pCrossCorr->Init(m_iStepsProcessed-2, obs->m_iDepth, g_bACFFFT);
	
	m_faInput1.SetSize(m_iStepsProcessed-2);
	m_faInput2.SetSize(m_iStepsProcessed-2);
	m_faOutput.SetSize(obs->m_iDepth);

	obs->m_faACF.resize(obs->m_iDepth);
	obs->m_faACF2.resize(obs->m_iDepth);
	obs->m_faACF3.resize(obs->m_iDepth);
	obs->m_faACF4.resize(obs->m_iDepth);
	for (z=0;z<obs->m_iDepth;z++) {
		obs->m_faACF[z] = 0;
		obs->m_faACF2[z] = 0;
		obs->m_faACF3[z] = 0;
		obs->m_faACF4[z] = 0;
	}

	if (obs->m_bSingleACFs && m_bSFG) {
		obs->m_faaSingleACFs.resize(27);
		for (z=0;z<27;z++) {
			obs->m_faaSingleACFs[z].resize(obs->m_iDepth);
			for (z2=0;z2<obs->m_iDepth;z2++)
				obs->m_faaSingleACFs[z][z2] = 0;
		}
	}

	if (obs->m_bUseCommutator) {
		obs->m_faACF_CTP.resize(obs->m_iDepth);
		obs->m_faACF2_CTP.resize(obs->m_iDepth);
		obs->m_faACF3_CTP.resize(obs->m_iDepth);
		obs->m_faACF4_CTP.resize(obs->m_iDepth);
		obs->m_faACF_CTN.resize(obs->m_iDepth);
		obs->m_faACF2_CTN.resize(obs->m_iDepth);
		obs->m_faACF3_CTN.resize(obs->m_iDepth);
		obs->m_faACF4_CTN.resize(obs->m_iDepth);
		for (z=0;z<obs->m_iDepth;z++) {
			obs->m_faACF_CTP[z] = 0;
			obs->m_faACF2_CTP[z] = 0;
			obs->m_faACF3_CTP[z] = 0;
			obs->m_faACF4_CTP[z] = 0;
			obs->m_faACF_CTN[z] = 0;
			obs->m_faACF2_CTN[z] = 0;
			obs->m_faACF3_CTN[z] = 0;
			obs->m_faACF4_CTN[z] = 0;
		}
		if (obs->m_bSingleACFs && m_bSFG) {
			obs->m_faaSingleACFs_CTP.resize(27);
			obs->m_faaSingleACFs_CTN.resize(27);
			for (z=0;z<27;z++) {
				obs->m_faaSingleACFs_CTP[z].resize(obs->m_iDepth);
				obs->m_faaSingleACFs_CTN[z].resize(obs->m_iDepth);
				for (z2=0;z2<obs->m_iDepth;z2++) {
					obs->m_faaSingleACFs_CTP[z][z2] = 0;
					obs->m_faaSingleACFs_CTN[z][z2] = 0;
				}
			}
		}
	}

	tr = m_oaTrajectories[0];

	tfs = 0;
	acc = 0;
	for (z=0;z<g_oaMolecules.GetSize();z++) {
		for (z2=0;z2<(int)obs->m_iaMolSelection[z].size();z2++) {
			tfs++;
			acc++;
			if (obs->m_bCrossCorrelation) {
				for (z3=0;z3<g_oaMolecules.GetSize();z3++) {
					for (z4=0;z4<(int)obs->m_iaMolSelection[z3].size();z4++) {
						if ((z == z3) && (z2 == z4))
							continue;
						tfs++;
					}
				}
			}
		}
	}
	tfs /= 60.0;


	pc = 0;
	mprintf(WHITE,"        [");
	for (z=0;z<g_oaMolecules.GetSize();z++) {

		for (z2=0;z2<(int)obs->m_iaMolSelection[z].size();z2++) {

			imol = ((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex[ obs->m_iaMolSelection[z][z2] ];

			ComputeACFPair(obs,tr->m_oaMolecules[imol],tr->m_oaMolecules[imol]);

			pc++;
			if (fmod(pc,tfs) < 1.0)
				mprintf(WHITE,"#");

			if (obs->m_bCrossCorrelation) {

					for (z3=0;z3<g_oaMolecules.GetSize();z3++) {
						for (z4=0;z4<(int)obs->m_iaMolSelection[z3].size();z4++) {
							if ((z == z3) && (z2 == z4))
								continue;
							imol2 = ((CMolecule*)g_oaMolecules[z3])->m_laSingleMolIndex[ obs->m_iaMolSelection[z3][z4] ];
							ComputeACFPair(obs,tr->m_oaMolecules[imol],tr->m_oaMolecules[imol2]);
							pc++;
							if (fmod(pc,tfs) < 1.0)
								mprintf(WHITE,"#");
						}
					}
			}
		}
	}
	mprintf(WHITE,"]\n");

	delete m_pCrossCorr;

		for (z=0;z<obs->m_iDepth;z++) {
			obs->m_faACF[z] /= pc;
			obs->m_faACF2[z] /= pc;
			obs->m_faACF3[z] /= pc;
			obs->m_faACF4[z] /= pc;
		}

		if (obs->m_bUseCommutator) {
			for (z=0;z<obs->m_iDepth;z++) {
				obs->m_faACF_CTP[z] /= pc;
				obs->m_faACF2_CTP[z] /= pc;
				obs->m_faACF3_CTP[z] /= pc;
				obs->m_faACF4_CTP[z] /= pc;
				obs->m_faACF_CTN[z] /= pc;
				obs->m_faACF2_CTN[z] /= pc;
				obs->m_faACF3_CTN[z] /= pc;
				obs->m_faACF4_CTN[z] /= pc;
			}
		}
}




void CROAEngine::ComputeSpectrum(CROAObservation *obs, const std::vector<double> &acf, std::vector<double> &spectrum, const char *name, bool im, bool ghack) {

	int z, o;
	double f, val;
	CFFT *fft;
	std::vector<double> temp;
	FILE *a;
	CxString buf, unit;


	if (obs->m_bSaveACF && (name != NULL)) {
		buf.sprintf("acf_%s%s%s",(const char*)obs->m_sTypeName,name,(const char*)obs->m_sMolName);
		if (obs->m_bCrossCorrelation)
			buf.strcat("_cc");
		if (m_bSmoothData)
			buf.strcat("_smooth");
		buf.strcat(".csv");
		mprintf("      Writing ACF to \"%s\"...\n",(const char*)buf);
		a = OpenFileWrite((const char*)buf,true);
		mfprintf(a,"#Depth/fs;  ACF\n");
		for (z=0;z<obs->m_iDepth;z++)
			mfprintf(a,"%.2f;  %.10G\n",(double)z*g_fTimestepLength*g_iStride,acf[z]);
		fclose(a);
	}

	temp.resize( 2*
		(obs->m_iDepth+obs->m_iZeroPadding) );

	for (z=0;z<obs->m_iDepth;z++)
		temp[z] = acf[z];


	for (z=obs->m_iDepth;z<
		(obs->m_iDepth+obs->m_iZeroPadding);z++)
		temp[z] = 0.0;

	{
		if (obs->m_iWindowFunction == 1) {
			for (z=0;z<obs->m_iDepth;z++)
				temp[z] *= pow2(cos((double)z / (obs->m_iDepth - 1) / 2.0 * Pi));
		} else if (obs->m_iWindowFunction == 2) {
			for (z=0;z<obs->m_iDepth;z++)
				temp[z] *= exp(-(double)z / obs->m_iWindowFunctionParameter);
		} else if (obs->m_iWindowFunction == 3) {
			for (z=0;z<obs->m_iDepth;z++)
				temp[z] *= exp(-(double)z * z / obs->m_iWindowFunctionParameter / obs->m_iWindowFunctionParameter);
		} else if (obs->m_iWindowFunction != 0) {
			eprintf("CROAEngine::ComputeSpectrum(): Error: Unknown window function: %d.\n",obs->m_iWindowFunction);
			abort();
		}
	}

	if (im) {
		o = obs->m_iDepth+obs->m_iZeroPadding;
		for (z=1;z<o;z++)
			temp[o+z] = -temp[o-z];
		temp[o] = 0.0;
	} else {
		o = obs->m_iDepth+obs->m_iZeroPadding;
		for (z=1;z<o;z++)
			temp[o+z] = temp[o-z];
		temp[o] = 0.0;
	}


	if (obs->m_bSaveACF && (name != NULL)) {
		buf.sprintf("acf_%s%s_windowed%s",(const char*)obs->m_sTypeName,name,(const char*)obs->m_sMolName);
		if (obs->m_bCrossCorrelation)
			buf.strcat("_cc");
		if (m_bSmoothData)
			buf.strcat("_smooth");
		buf.strcat(".csv");
		mprintf("      Writing windowed/mirrored ACF to \"%s\"...\n",(const char*)buf);
		a = OpenFileWrite((const char*)buf,true);
		mfprintf(a,"#Depth/fs;  ACF\n");
		for (z=0;z<2*
			(obs->m_iDepth+obs->m_iZeroPadding);z++)
			mfprintf(a,"%.2f;  %.10G\n",(double)z*g_fTimestepLength*g_iStride,temp[z]);
		fclose(a);
	}

	try { fft = new CFFT(); } catch(...) { fft = NULL; }
	if (fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);

	fft->PrepareFFT_C2C(2*
		(obs->m_iDepth+obs->m_iZeroPadding));

	for (z=0;z<2*
		(obs->m_iDepth+obs->m_iZeroPadding);z++) {
		fft->m_pInput[2*z] = temp[z];
		fft->m_pInput[2*z+1] = 0.0;
	}
	fft->DoFFT();
	
	spectrum.resize(obs->m_iSpecLength);

	unit = "";
	switch(obs->m_iType) {

		case ROA_SPECTRUM_IR:
			for (z=0;z<obs->m_iSpecLength;z++) {
				if (im)
					val = fft->m_pOutput[2*z+1];
				else
					val = fft->m_pOutput[2*z];
				spectrum[z] = 3047.2310 * val * g_fTimestepLength * g_iStride;
			}
			if (obs->m_bCorrectTemperature)
				unit.sprintf(";  Spectrum at %.2f K (cm*km*mol^-1);  Integral (km*mol^-1)",obs->m_fCorrectTemperature);
			else
				unit.sprintf(";  Spectrum (K*cm*km*mol^-1);  Integral (K*km*mol^-1)");
			break;

		case ROA_SPECTRUM_RAMAN:
			for (z=1;z<obs->m_iSpecLength;z++) {
				f = z * obs->m_fSpecResolution;
				if (im)
					val = fft->m_pOutput[2*z+1];
				else
					val = fft->m_pOutput[2*z];
				spectrum[z] = 1.0E40 * 100.0 * 2.0 * CONST_PLANCK / (8.0 * CONST_BOLTZMANN * CONST_EPSILON0 * CONST_EPSILON0) * 1.0E6 * CONST_ELEMENTARY_CHARGE * CONST_ELEMENTARY_CHARGE * 1.0E-48 / 1.0E-15 * pow4(obs->m_fLaserFreq - f) / f / (1.0 - exp(-1.438777 * f / obs->m_fTemperature)) * val * g_fTimestepLength * g_iStride; // Output in 1e-30*m^2*cm
			}
			spectrum[0] = 0.0;
			if (obs->m_bCorrectTemperature)
				unit.sprintf(";  Spectrum at %.2f K (10^-40*m^2*cm);  Integral (10^-40*m^2)",obs->m_fCorrectTemperature);
			else
				unit.sprintf(";  Spectrum (10^-40*K*m^2*cm);  Integral (10^-40*K*m^2)");
			break;

		case ROA_SPECTRUM_ROA:
			for (z=1;z<obs->m_iSpecLength;z++) {
				f = z * obs->m_fSpecResolution;
				if (im)
					val = fft->m_pOutput[2*z+1];
				else
					val = fft->m_pOutput[2*z];
				spectrum[z] = 1.0/90.0 * 1.0E40 * 100.0 * 2.0 * CONST_PLANCK / (8.0 * CONST_BOLTZMANN * CONST_EPSILON0 * CONST_EPSILON0) * 1.0E6 * CONST_ELEMENTARY_CHARGE * CONST_ELEMENTARY_CHARGE / CONST_SPEED_OF_LIGHT * 1.0E3 * 1.0E-48 / 1.0E-15 * pow4(obs->m_fLaserFreq - f) / f / (1.0 - exp(-1.438777 * f / obs->m_fTemperature)) * val * g_fTimestepLength * g_iStride;
				if (ghack)
					spectrum[z] *= obs->m_fLaserFreq / f;
			}
			spectrum[0] = 0.0;
			if (obs->m_bCorrectTemperature)
				unit.sprintf(";  Spectrum at %.2f K (10^-40*m^2*cm);  Integral (10^-40*m^2)",obs->m_fCorrectTemperature);
			else
				unit.sprintf(";  Spectrum (10^-40*K*m^2*cm);  Integral (10^-40*K*m^2)");
			break;

		case ROA_SPECTRUM_VCD:
			for (z=0;z<obs->m_iSpecLength;z++) {
				if (im)
					val = fft->m_pOutput[2*z+1];
				else
					val = fft->m_pOutput[2*z];
				spectrum[z] = 2.0 * 28.260058 * val * g_fTimestepLength * g_iStride;
			}
			if (obs->m_bCorrectTemperature)
				unit.sprintf(";  Spectrum at %.2f K (cm*km*mol^-1);  Integral (km*mol^-1)",obs->m_fCorrectTemperature);
			else
				unit.sprintf(";  Spectrum (K*cm*km*mol^-1);  Integral (K*km*mol^-1)");
			break;

		case ROA_SPECTRUM_SFG:
			for (z=1;z<obs->m_iSpecLength;z++) {
				f = z * obs->m_fSpecResolution;
				if (im)
					val = fft->m_pOutput[2*z+1];
				else
					val = fft->m_pOutput[2*z];
				spectrum[z] = val / f * g_fTimestepLength * g_iStride; 
			}
			spectrum[0] = 0.0;
			if (obs->m_bCorrectTemperature)
				unit.sprintf(";  Spectrum at %.2f K;  Integral",obs->m_fCorrectTemperature);
			else
				unit.sprintf(";  Spectrum;  Integral");
			break;

		default:
			for (z=0;z<obs->m_iSpecLength;z++) {
				if (im)
					val = fft->m_pOutput[2*z+1];
				else
					val = fft->m_pOutput[2*z];
				spectrum[z] = val * g_fTimestepLength * g_iStride;
			}
			unit.sprintf(";  Spectrum;  Integral");
			break;
	}

	if (obs->m_bFiniteDifferenceCorrection) {
		f = obs->m_fSpecResolution * g_fTimestepLength * g_iStride * 1.883652e-4;
		for (z=1;z<obs->m_iSpecLength;z++)
			spectrum[z] *= pow2(f * z / sin(f * z)); // Divide by sinc function to correct finite difference derivation
	}
	

	if (obs->m_bCorrectTemperature) {
		if (name != NULL) // Silent if no file output
			mprintf("      Correcting spectrum for %.2f K simulation temperature...\n",obs->m_fCorrectTemperature);
		for (z=0;z<obs->m_iSpecLength;z++)
			spectrum[z] /= obs->m_fCorrectTemperature;
	}

	if (name != NULL) {

		buf.sprintf("spectrum_%s%s%s",(const char*)obs->m_sTypeName,name,(const char*)obs->m_sMolName);
		if (obs->m_bCrossCorrelation)
			buf.strcat("_cc");
		if (m_bSmoothData)
			buf.strcat("_smooth");
		buf.strcat(".csv");
		mprintf("      Writing spectrum to \"%s\"...\n",(const char*)buf);
		a = OpenFileWrite((const char*)buf,true);

		if (obs->m_bCorrectFrequency)
			mfprintf(a, "#Corrected Wavenumber (cm^-1)%s\n",(const char*)unit);
		else
			mfprintf(a, "#Wavenumber (cm^-1)%s\n",(const char*)unit);

		f = 0.0;
		for (z=0;z<obs->m_iSpecLength;z++) {

			f += spectrum[z] * obs->m_fSpecResolution;

			if (obs->m_bCorrectFrequency)
				mfprintf(a, "%.2f; %.8G; %.14G\n", CorrectWavenumber(obs->m_fSpecResolution * z), spectrum[z], f);
			else
				mfprintf(a, "%.2f; %.8G; %.14G\n", obs->m_fSpecResolution * z, spectrum[z], f);
		}

		fclose(a);
	}
}


void CROAEngine::WriteSpectrum(CROAObservation *obs, const std::vector<double> &spectrum, const char *name) {

	int z;
	double f;
	FILE *a;
	CxString buf, unit;


	buf.sprintf("spectrum_%s%s%s",(const char*)obs->m_sTypeName,name,(const char*)obs->m_sMolName);
	if (obs->m_bCrossCorrelation)
		buf.strcat("_cc");
	if (m_bSmoothData)
		buf.strcat("_smooth");
	buf.strcat(".csv");
	mprintf("      Writing spectrum to \"%s\"...\n",(const char*)buf);
	a = OpenFileWrite((const char*)buf,true);

	unit = "";
	switch(obs->m_iType) {

		case ROA_SPECTRUM_IR:
			if (obs->m_bCorrectTemperature)
				unit.sprintf(";  Spectrum at %.2f K (cm*km*mol^-1);  Integral (km*mol^-1)",obs->m_fCorrectTemperature);
			else
				unit.sprintf(";  Spectrum (K*cm*km*mol^-1);  Integral (K*km*mol^-1)");
			break;

		case ROA_SPECTRUM_RAMAN:
			if (obs->m_bCorrectTemperature)
				unit.sprintf(";  Spectrum at %.2f K (10^-40*m^2*cm);  Integral (10^-40*m^2)",obs->m_fCorrectTemperature);
			else
				unit.sprintf(";  Spectrum (10^-40*K*m^2*cm);  Integral (10^-40*K*m^2)");
			break;

		case ROA_SPECTRUM_ROA:
			if (obs->m_bCorrectTemperature)
				unit.sprintf(";  Spectrum at %.2f K (10^-40*m^2*cm);  Integral (10^-40*m^2)",obs->m_fCorrectTemperature);
			else
				unit.sprintf(";  Spectrum (10^-40*K*m^2*cm);  Integral (10^-40*K*m^2)");
			break;

		case ROA_SPECTRUM_VCD:
			if (obs->m_bCorrectTemperature)
				unit.sprintf(";  Spectrum at %.2f K (cm*km*mol^-1);  Integral (km*mol^-1)",obs->m_fCorrectTemperature);
			else
				unit.sprintf(";  Spectrum (K*cm*km*mol^-1);  Integral (K*km*mol^-1)");
			break;

		case ROA_SPECTRUM_SFG:
			if (obs->m_bCorrectTemperature)
				unit.sprintf(";  Spectrum at %.2f K;  Integral",obs->m_fCorrectTemperature);
			else
				unit.sprintf(";  Spectrum;  Integral");
			break;

		default:
			unit.sprintf(";  Spectrum;  Integral");
			break;
	}

	if (obs->m_bCorrectFrequency)
		mfprintf(a, "#Corrected Wavenumber (cm^-1)%s\n",(const char*)unit);
	else
		mfprintf(a, "#Wavenumber (cm^-1)%s\n",(const char*)unit);

	f = 0.0;
	for (z=0;z<obs->m_iSpecLength;z++) {

		f += spectrum[z] * obs->m_fSpecResolution;

		if (obs->m_bCorrectFrequency)
			mfprintf(a, "%.2f; %.8G; %.14G\n", CorrectWavenumber(obs->m_fSpecResolution * z), spectrum[z], f);
		else
			mfprintf(a, "%.2f; %.8G; %.14G\n", obs->m_fSpecResolution * z, spectrum[z], f);
	}

	fclose(a);
}


void CROAEngine::WriteSpectrumSFG(CROAObservation *obs, const std::vector<double> &spectrumre, const std::vector<double> &spectrumim, const char *name) {

	int z;
	FILE *a;
	CxString buf;


	buf.sprintf("spectrum_%s%s%s",(const char*)obs->m_sTypeName,name,(const char*)obs->m_sMolName);
	if (obs->m_bCrossCorrelation)
		buf.strcat("_cc");
	if (m_bSmoothData)
		buf.strcat("_smooth");
	buf.strcat(".csv");
	mprintf("      Writing spectrum to \"%s\"...\n",(const char*)buf);
	a = OpenFileWrite((const char*)buf,true);

	if (obs->m_bCorrectFrequency)
		mfprintf(a, "#Corrected Wavenumber (cm^-1)");
	else
		mfprintf(a, "#Wavenumber (cm^-1)");

	mfprintf(a,";  Spectrum Re;  Spectrum Im;  Spectrum Re^2;  Spectrum Im^2;  Spectrum Abs^2\n");

	for (z=0;z<obs->m_iSpecLength;z++) {

		if (obs->m_bCorrectFrequency)
			mfprintf(a, "%.2f", CorrectWavenumber(obs->m_fSpecResolution * z) );
		else
			mfprintf(a, "%.2f", obs->m_fSpecResolution * z );

		mfprintf(a, ";  %.8G;  %.8G;  %.8G;  %.8G;  %.8G\n", spectrumre[z], spectrumim[z], SQR(spectrumre[z]), SQR(spectrumim[z]), SQR(spectrumre[z])+SQR(spectrumim[z]) );
	}

	fclose(a);
}


void CROAEngine::ComputeSpectrum(CROAObservation *obs) {

	int z;
	std::vector<double> spectrum;
	CxString buf;


	/*********************************************************************************************************************/
	if (obs->m_iType == ROA_SPECTRUM_IR) {

		ComputeSpectrum(obs,obs->m_faACF,obs->m_faSpectrum,"",false);

	/*********************************************************************************************************************/
	} else if (obs->m_iType == ROA_SPECTRUM_RAMAN) {

		if (!obs->m_bUseCommutator || g_bAdvanced2) {
			ComputeSpectrum(obs,obs->m_faACF,obs->m_faSpectrum,"_iso",false);
			ComputeSpectrum(obs,obs->m_faACF2,obs->m_faSpectrum2,"_aniso",false);

			spectrum.resize(obs->m_iSpecLength);

			for (z=0;z<obs->m_iSpecLength;z++)
				spectrum[z] = obs->m_faSpectrum[z] + 4.0 / 45.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_para");

			for (z=0;z<obs->m_iSpecLength;z++)
				spectrum[z] = obs->m_faSpectrum2[z] / 15.0;
			WriteSpectrum(obs,spectrum,"_ortho");

			for (z=0;z<obs->m_iSpecLength;z++)
				spectrum[z] = obs->m_faSpectrum[z] + 7.0 / 45.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_unpol");

			for (z=0;z<obs->m_iSpecLength;z++)
				if ((obs->m_faSpectrum[z] + 4.0 / 45.0 * obs->m_faSpectrum2[z]) != 0)
					spectrum[z] = (obs->m_faSpectrum2[z] / 15.0) / (obs->m_faSpectrum[z] + 4.0 / 45.0 * obs->m_faSpectrum2[z]);
				else
					spectrum[z] = 0;
			WriteSpectrum(obs,spectrum,"_depol_ratio");
		}

		if (obs->m_bUseCommutator) {
			ComputeSpectrum(obs,obs->m_faACF_CTP,obs->m_faSpectrum,"_iso_ct",false);
			ComputeSpectrum(obs,obs->m_faACF2_CTP,obs->m_faSpectrum2,"_aniso_ct",false);

			spectrum.resize(obs->m_iSpecLength);

			for (z=0;z<obs->m_iSpecLength;z++)
				spectrum[z] = obs->m_faSpectrum[z] + 4.0 / 45.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_para_ct");

			for (z=0;z<obs->m_iSpecLength;z++)
				spectrum[z] = obs->m_faSpectrum2[z] / 15.0;
			WriteSpectrum(obs,spectrum,"_ortho_ct");

			for (z=0;z<obs->m_iSpecLength;z++)
				spectrum[z] = obs->m_faSpectrum[z] + 7.0 / 45.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_unpol_ct");

			for (z=0;z<obs->m_iSpecLength;z++)
				if ((obs->m_faSpectrum[z] + 4.0 / 45.0 * obs->m_faSpectrum2[z]) > 0)
					spectrum[z] = (obs->m_faSpectrum2[z] / 15.0) / (obs->m_faSpectrum[z] + 4.0 / 45.0 * obs->m_faSpectrum2[z]);
				else
					spectrum[z] = 0;
			WriteSpectrum(obs,spectrum,"_depol_ratio_ct");
		}

	/*********************************************************************************************************************/
	} else if (obs->m_iType == ROA_SPECTRUM_VCD) {

		if (!obs->m_bUseCommutator || g_bAdvanced2)
			ComputeSpectrum(obs,obs->m_faACF,obs->m_faSpectrum,"",true);

		if (obs->m_bUseCommutator) {
	//		ComputeSpectrum(obs,obs->m_faACF_CTP,obs->m_faSpectrum,"_ct_p",false);
	//		ComputeSpectrum(obs,obs->m_faACF_CTP,obs->m_faSpectrum,"_ct_p_im",true);
	//		ComputeSpectrum(obs,obs->m_faACF_CTN,obs->m_faSpectrum,"_ct_n_re",false);
			ComputeSpectrum(obs,obs->m_faACF_CTN,obs->m_faSpectrum,"_ct",true);
		}

	/*********************************************************************************************************************/
	} else if (obs->m_iType == ROA_SPECTRUM_SFG) {

		ComputeSpectrum(obs,obs->m_faACF,obs->m_faSpectrum,NULL,false); // Re
		ComputeSpectrum(obs,obs->m_faACF,obs->m_faSpectrum2,NULL,true); // Im
		WriteSpectrumSFG(obs,obs->m_faSpectrum,obs->m_faSpectrum2,"");

		if (obs->m_bSingleACFs) {
			for (z=0;z<27;z++) {
				buf.sprintf("_single_%c%c_%c",'x'+(z/9),'x'+((z/3)%3),'x'+(z%3));
				ComputeSpectrum(obs,obs->m_faaSingleACFs[z],obs->m_faSpectrum,NULL,false); // Re
				ComputeSpectrum(obs,obs->m_faaSingleACFs[z],obs->m_faSpectrum2,NULL,true); // Im
				WriteSpectrumSFG(obs,obs->m_faSpectrum,obs->m_faSpectrum2,(const char*)buf);
			}
		}

		if (obs->m_bUseCommutator) {
			ComputeSpectrum(obs,obs->m_faACF_CTP,obs->m_faSpectrum,NULL,false); // Re
			ComputeSpectrum(obs,obs->m_faACF_CTP,obs->m_faSpectrum2,NULL,true); // Im
			WriteSpectrumSFG(obs,obs->m_faSpectrum,obs->m_faSpectrum2,"_ct");

			if (obs->m_bSingleACFs) {
				for (z=0;z<27;z++) {
					buf.sprintf("_single_%c%c_%c_ct",'x'+(z/9),'x'+((z/3)%3),'x'+(z%3));
					ComputeSpectrum(obs,obs->m_faaSingleACFs_CTP[z],obs->m_faSpectrum,NULL,false); // Re
					ComputeSpectrum(obs,obs->m_faaSingleACFs_CTP[z],obs->m_faSpectrum2,NULL,true); // Im
					WriteSpectrumSFG(obs,obs->m_faSpectrum,obs->m_faSpectrum2,(const char*)buf);
				}
			}
		}

	/*********************************************************************************************************************/
	} else if (obs->m_iType == ROA_SPECTRUM_ROA) {

		if (!obs->m_bUseCommutator || g_bAdvanced2) {

			ComputeSpectrum(obs,obs->m_faACF,obs->m_faSpectrum,"_ag",true,true);
			ComputeSpectrum(obs,obs->m_faACF2,obs->m_faSpectrum2,"_gamma_g",true,true);
			ComputeSpectrum(obs,obs->m_faACF3,obs->m_faSpectrum3,"_gamma_a",false);

			spectrum.resize(obs->m_iSpecLength);

			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 180.0 * obs->m_faSpectrum[z] + 28.0 * obs->m_faSpectrum2[z] +  4.0 * obs->m_faSpectrum3[z];
			WriteSpectrum(obs,spectrum,"_90deg_ortho");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] =                                24.0 * obs->m_faSpectrum2[z] -  8.0 * obs->m_faSpectrum3[z];
			WriteSpectrum(obs,spectrum,"_90deg_para");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 180.0 * obs->m_faSpectrum[z] + 52.0 * obs->m_faSpectrum2[z] -  4.0 * obs->m_faSpectrum3[z];
			WriteSpectrum(obs,spectrum,"_90deg_unpol");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 720.0 * obs->m_faSpectrum[z] + 16.0 * obs->m_faSpectrum2[z] - 16.0 * obs->m_faSpectrum3[z];
			WriteSpectrum(obs,spectrum,"_forward");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] =                                96.0 * obs->m_faSpectrum2[z] + 32.0 * obs->m_faSpectrum3[z];
			WriteSpectrum(obs,spectrum,"_backward");

	/*		for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 180.0 * obs->m_faSpectrum[z] + 28.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_90deg_ortho_noa");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] =                                24.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_90deg_para_noa");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 180.0 * obs->m_faSpectrum[z] + 52.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_90deg_unpol_noa");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 720.0 * obs->m_faSpectrum[z] + 16.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_forward_noa");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] =                                96.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_backward_noa");

			ComputeSpectrum(obs,obs->m_faACF,obs->m_faSpectrum,"_ag_re",false,true);
			ComputeSpectrum(obs,obs->m_faACF2,obs->m_faSpectrum2,"_gamma_g_re",false,true);
			ComputeSpectrum(obs,obs->m_faACF3,obs->m_faSpectrum3,"_gamma_a_im",true);*/
		}

		if (obs->m_bUseCommutator) {

			ComputeSpectrum(obs,obs->m_faACF_CTN,obs->m_faSpectrum,"_ag_ct",true,true);
			ComputeSpectrum(obs,obs->m_faACF2_CTN,obs->m_faSpectrum2,"_gamma_g_ct",true,true);
			ComputeSpectrum(obs,obs->m_faACF3_CTP,obs->m_faSpectrum3,"_gamma_a_ct",false);

			spectrum.resize(obs->m_iSpecLength);

			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 180.0 * obs->m_faSpectrum[z] + 28.0 * obs->m_faSpectrum2[z] +  4.0 * obs->m_faSpectrum3[z];
			WriteSpectrum(obs,spectrum,"_90deg_ortho_ct");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] =                                24.0 * obs->m_faSpectrum2[z] -  8.0 * obs->m_faSpectrum3[z];
			WriteSpectrum(obs,spectrum,"_90deg_para_ct");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 180.0 * obs->m_faSpectrum[z] + 52.0 * obs->m_faSpectrum2[z] -  4.0 * obs->m_faSpectrum3[z];
			WriteSpectrum(obs,spectrum,"_90deg_unpol_ct");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 720.0 * obs->m_faSpectrum[z] + 16.0 * obs->m_faSpectrum2[z] - 16.0 * obs->m_faSpectrum3[z];
			WriteSpectrum(obs,spectrum,"_forward_ct");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] =                                96.0 * obs->m_faSpectrum2[z] + 32.0 * obs->m_faSpectrum3[z];
			WriteSpectrum(obs,spectrum,"_backward_ct");

	/*		for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 180.0 * obs->m_faSpectrum[z] + 28.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_90deg_ortho_noa_ct");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] =                                24.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_90deg_para_noa_ct");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 180.0 * obs->m_faSpectrum[z] + 52.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_90deg_unpol_noa_ct");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] = 720.0 * obs->m_faSpectrum[z] + 16.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_forward_noa_ct");
			for (z=0;z<obs->m_iSpecLength;z++) spectrum[z] =                                96.0 * obs->m_faSpectrum2[z];
			WriteSpectrum(obs,spectrum,"_backward_noa_ct");

			ComputeSpectrum(obs,obs->m_faACF_CTP,obs->m_faSpectrum,"_ag_ct_p",false,true);
			ComputeSpectrum(obs,obs->m_faACF2_CTP,obs->m_faSpectrum2,"_gamma_g_ct_p",false,true);
			ComputeSpectrum(obs,obs->m_faACF3_CTN,obs->m_faSpectrum3,"_gamma_a_ct_n",true);*/
		}

	/*********************************************************************************************************************/
	} else {

		eprintf("CROAEngine::ComputeSpectrum(): Error: Unknown spectrum type %d.\n",obs->m_iType);
		abort();

	}
}


bool CROAEngine::PrepareVori(CTimeStep *ts, bool dipole, bool quadrupole, bool magmom) {

	try { g_pVoroWrapper = new CVoroWrapper(); } catch(...) { g_pVoroWrapper = NULL; }
	if (g_pVoroWrapper == NULL) NewException((double)sizeof(CVoroWrapper),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { g_pTetraPak = new CTetraPak(); } catch(...) { g_pTetraPak = NULL; }
	if (g_pTetraPak == NULL) NewException((double)sizeof(CTetraPak),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	g_bTegri = true;
	g_iTrajFormat = 5; // Cube format

	g_bVoroIntegrateCharge = true;
	g_bVoroIntegrateDipoleMoment = dipole;
	g_bVoroIntegrateQuadrupoleMoment = quadrupole;
	g_bVoroIntegrateTotalCurrent = magmom;
	g_bVoroIntegrateMagneticMoment = magmom;
	g_bCubeTimeDev = magmom;

	g_pTetraPak->ParseSilent(ts);

	for (int i = 0; i < g_oaMolecules.GetSize(); i++) {
		int rty;
		ParseAtom("#2", i, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterType, rty, ((CMolecule *)g_oaMolecules[i])->m_iDipoleCenterIndex);
		((CMolecule *)g_oaMolecules[i])->m_iDipoleMode = 7;
	}

	return true;
}


bool CROAEngine::PrepareCurrentDensity(CTimeStep *ts, CTimeStep *tsref, bool first) {

	if (first) {

		m_pDensity.m_iRes[0] = tsref->m_pVolumetricData->m_iRes[0];
		m_pDensity.m_iRes[1] = tsref->m_pVolumetricData->m_iRes[1];
		m_pDensity.m_iRes[2] = tsref->m_pVolumetricData->m_iRes[2];
		m_pDensity.m_fMinVal[0] = tsref->m_pVolumetricData->m_fMinVal[0];
		m_pDensity.m_fMaxVal[0] = tsref->m_pVolumetricData->m_fMaxVal[0];
		m_pDensity.m_fMinVal[1] = tsref->m_pVolumetricData->m_fMinVal[1];
		m_pDensity.m_fMaxVal[1] = tsref->m_pVolumetricData->m_fMaxVal[1];
		m_pDensity.m_fMinVal[2] = tsref->m_pVolumetricData->m_fMinVal[2];
		m_pDensity.m_fMaxVal[2] = tsref->m_pVolumetricData->m_fMaxVal[2];
		m_pDensity.Create();

		m_pDensityGrad[0].m_iRes[0] = tsref->m_pVolumetricData->m_iRes[0];
		m_pDensityGrad[0].m_iRes[1] = tsref->m_pVolumetricData->m_iRes[1];
		m_pDensityGrad[0].m_iRes[2] = tsref->m_pVolumetricData->m_iRes[2];
		m_pDensityGrad[0].m_fMinVal[0] = tsref->m_pVolumetricData->m_fMinVal[0];
		m_pDensityGrad[0].m_fMaxVal[0] = tsref->m_pVolumetricData->m_fMaxVal[0];
		m_pDensityGrad[0].m_fMinVal[1] = tsref->m_pVolumetricData->m_fMinVal[1];
		m_pDensityGrad[0].m_fMaxVal[1] = tsref->m_pVolumetricData->m_fMaxVal[1];
		m_pDensityGrad[0].m_fMinVal[2] = tsref->m_pVolumetricData->m_fMinVal[2];
		m_pDensityGrad[0].m_fMaxVal[2] = tsref->m_pVolumetricData->m_fMaxVal[2];
		m_pDensityGrad[0].Create();

		m_pDensityGrad[1].m_iRes[0] = tsref->m_pVolumetricData->m_iRes[0];
		m_pDensityGrad[1].m_iRes[1] = tsref->m_pVolumetricData->m_iRes[1];
		m_pDensityGrad[1].m_iRes[2] = tsref->m_pVolumetricData->m_iRes[2];
		m_pDensityGrad[1].m_fMinVal[0] = tsref->m_pVolumetricData->m_fMinVal[0];
		m_pDensityGrad[1].m_fMaxVal[0] = tsref->m_pVolumetricData->m_fMaxVal[0];
		m_pDensityGrad[1].m_fMinVal[1] = tsref->m_pVolumetricData->m_fMinVal[1];
		m_pDensityGrad[1].m_fMaxVal[1] = tsref->m_pVolumetricData->m_fMaxVal[1];
		m_pDensityGrad[1].m_fMinVal[2] = tsref->m_pVolumetricData->m_fMinVal[2];
		m_pDensityGrad[1].m_fMaxVal[2] = tsref->m_pVolumetricData->m_fMaxVal[2];
		m_pDensityGrad[1].Create();

		m_pDensityGrad[2].m_iRes[0] = tsref->m_pVolumetricData->m_iRes[0];
		m_pDensityGrad[2].m_iRes[1] = tsref->m_pVolumetricData->m_iRes[1];
		m_pDensityGrad[2].m_iRes[2] = tsref->m_pVolumetricData->m_iRes[2];
		m_pDensityGrad[2].m_fMinVal[0] = tsref->m_pVolumetricData->m_fMinVal[0];
		m_pDensityGrad[2].m_fMaxVal[0] = tsref->m_pVolumetricData->m_fMaxVal[0];
		m_pDensityGrad[2].m_fMinVal[1] = tsref->m_pVolumetricData->m_fMinVal[1];
		m_pDensityGrad[2].m_fMaxVal[1] = tsref->m_pVolumetricData->m_fMaxVal[1];
		m_pDensityGrad[2].m_fMinVal[2] = tsref->m_pVolumetricData->m_fMinVal[2];
		m_pDensityGrad[2].m_fMaxVal[2] = tsref->m_pVolumetricData->m_fMaxVal[2];
		m_pDensityGrad[2].Create();
	}

	ts->m_pVolumetricDataTimeDev = new C3DF<VORI_FLOAT>;
	ts->m_pVolumetricDataTimeDev->m_iRes[0] = tsref->m_pVolumetricData->m_iRes[0];
	ts->m_pVolumetricDataTimeDev->m_iRes[1] = tsref->m_pVolumetricData->m_iRes[1];
	ts->m_pVolumetricDataTimeDev->m_iRes[2] = tsref->m_pVolumetricData->m_iRes[2];
	ts->m_pVolumetricDataTimeDev->m_fMinVal[0] = tsref->m_pVolumetricData->m_fMinVal[0];
	ts->m_pVolumetricDataTimeDev->m_fMaxVal[0] = tsref->m_pVolumetricData->m_fMaxVal[0];
	ts->m_pVolumetricDataTimeDev->m_fMinVal[1] = tsref->m_pVolumetricData->m_fMinVal[1];
	ts->m_pVolumetricDataTimeDev->m_fMaxVal[1] = tsref->m_pVolumetricData->m_fMaxVal[1];
	ts->m_pVolumetricDataTimeDev->m_fMinVal[2] = tsref->m_pVolumetricData->m_fMinVal[2];
	ts->m_pVolumetricDataTimeDev->m_fMaxVal[2] = tsref->m_pVolumetricData->m_fMaxVal[2];
	ts->m_pVolumetricDataTimeDev->Create();

	ts->m_pCurrentDensity = new CxDoubleArray();

	return true;
}


bool CROAEngine::CalculateCurrentDensity(CTimeStep *ts1, CTimeStep *ts2, CTimeStep *ts3, bool reset) {

	int i, j, k;
	double integral = 0.0;
	int timedevSize = 0;
	int res[3];
	bool c=true;


	for (i = 0; i < ts2->m_pVolumetricData->m_iRes[0] * ts2->m_pVolumetricData->m_iRes[1] * ts2->m_pVolumetricData->m_iRes[2]; i++) {
		if (fabs(ts2->m_pVolumetricData->m_pBin[i]) > -1.0e-20) {
			ts2->m_pVolumetricDataTimeDev->m_pBin[i] = (ts3->m_pVolumetricData->m_pBin[i] - ts1->m_pVolumetricData->m_pBin[i]) / (2.0 * g_fTimestepLength);
			integral += ts2->m_pVolumetricDataTimeDev->m_pBin[i];
			timedevSize++;
		} else
			ts2->m_pVolumetricDataTimeDev->m_pBin[i] = 0.0;
	}
	for (i = 0; i < ts2->m_pVolumetricDataTimeDev->m_iRes[0] * ts2->m_pVolumetricDataTimeDev->m_iRes[1] * ts2->m_pVolumetricDataTimeDev->m_iRes[2]; i++)
		if (fabs(ts2->m_pVolumetricData->m_pBin[i]) > -1.0e-20)
			ts2->m_pVolumetricDataTimeDev->m_pBin[i] -= integral / timedevSize;

	res[0] = ts2->m_pVolumetricDataTimeDev->m_iRes[0];
	res[1] = ts2->m_pVolumetricDataTimeDev->m_iRes[1];
	res[2] = ts2->m_pVolumetricDataTimeDev->m_iRes[2];
	
	m_pDensity.CopyFrom(ts2->m_pVolumetricData);
	
	for (i = 0; i < res[2]; i++) {
		for (j = 0; j < res[1]; j++) {
			for (k = 1; k < res[0] - 1; k++)
				m_pDensityGrad[0].m_pBin[i * res[0] * res[1] + j * res[0] + k] = (m_pDensity.m_pBin[i * res[0] * res[1] + j * res[0] + k + 1] - m_pDensity.m_pBin[i * res[0] * res[1] + j * res[0] + k - 1]) / (2.0 * g_fCubeXStep);
			m_pDensityGrad[0].m_pBin[i * res[0] * res[1] + j * res[0]] = (m_pDensity.m_pBin[i * res[0] * res[1] + j * res[0] + 1] - m_pDensity.m_pBin[i * res[0] * res[1] + j * res[0] + res[0] - 1]) / (2.0 * g_fCubeXStep);
			m_pDensityGrad[0].m_pBin[i * res[0] * res[1] + j * res[0] + res[0] - 1] = (m_pDensity.m_pBin[i * res[0] * res[1] + j * res[0]] - m_pDensity.m_pBin[i * res[0] * res[1] + j * res[0] + res[0] - 2]) / (2.0 * g_fCubeXStep);
		}
	}
	
	for (i = 0; i < res[2]; i++) {
		for (j = 1; j < res[1] - 1; j++)
			for (k = 0; k < res[0]; k++)
				m_pDensityGrad[1].m_pBin[i * res[0] * res[1] + j * res[0] + k] = (m_pDensity.m_pBin[i * res[0] * res[1] + (j + 1) * res[0] + k] - m_pDensity.m_pBin[i * res[0] * res[1] + (j - 1) * res[0] + k]) / (2.0 * g_fCubeYStep);
		for (k = 0; k < res[0]; k++) {
			m_pDensityGrad[1].m_pBin[i * res[0] * res[1] + k] = (m_pDensity.m_pBin[i * res[0] * res[1] + res[0] + k] - m_pDensity.m_pBin[i * res[0] * res[1] + (res[1] - 1) * res[0] + k]) / (2.0 * g_fCubeYStep);
			m_pDensityGrad[1].m_pBin[i * res[0] * res[1] + (res[1] - 1) * res[0] + k] = (m_pDensity.m_pBin[i * res[0] * res[1] + k] - m_pDensity.m_pBin[i * res[0] * res[1] + (res[1] - 2) * res[0] + k]) / (2.0 * g_fCubeYStep);
		}
	}
	
	for (i = 1; i < res[2] - 1; i++)
		for (j = 0; j < res[1]; j++)
			for (k = 0; k < res[0]; k++)
				m_pDensityGrad[2].m_pBin[i * res[0] * res[1] + j * res[0] + k] = (m_pDensity.m_pBin[(i + 1) * res[0] * res[1] + j * res[0] + k] - m_pDensity.m_pBin[(i - 1) * res[0] * res[1] + j * res[0] + k]) / (2.0 * g_fCubeZStep);
	for (j = 0; j < res[1]; j++) {
		for (k = 0; k < res[0]; k++) {
			m_pDensityGrad[2].m_pBin[j * res[0] + k] = (m_pDensity.m_pBin[res[0] * res[1] + j * res[0] + k] - m_pDensity.m_pBin[(res[2] - 1) * res[0] * res[1] + j * res[0] + k]) / (2.0 * g_fCubeZStep);
			m_pDensityGrad[2].m_pBin[(res[2] - 1) * res[0] * res[1] + j * res[0] + k] = (m_pDensity.m_pBin[j * res[0] + k] - m_pDensity.m_pBin[(res[2] - 2) * res[0] * res[1] + j * res[0] + k]) / (2.0 * g_fCubeZStep);
		}
	}
	
	for (i = 0; i < res[0] * res[1] * res[2]; i++)
		m_pDensity.m_pBin[i] += g_fBackgroundDensity;
	

	CSparseMatrix pdeMatrix;
	CCurrentPDEDiscretizer::discretize(&pdeMatrix, m_pDensity, m_pDensityGrad[0], m_pDensityGrad[1], m_pDensityGrad[2]);
	
	static double thresh;

	if (m_pPDESolution == NULL) {

		try { m_pPDESolution = new double[res[0] * res[1] * res[2]]; } catch(...) { m_pPDESolution = NULL; }
		if (m_pPDESolution == NULL) NewException((double)res[0] * res[1] * res[2] * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);

		memset(m_pPDESolution, 0, res[0] * res[1] * res[2] * sizeof(double));
		thresh = g_fPDEConvThresh * CCurrentPDESolver::calcResidual(&pdeMatrix, m_pPDESolution, ts2->m_pVolumetricDataTimeDev->m_pBin);

	} else if (reset) {

		memset(m_pPDESolution, 0, res[0] * res[1] * res[2] * sizeof(double));
		thresh = g_fPDEConvThresh * CCurrentPDESolver::calcResidual(&pdeMatrix, m_pPDESolution, ts2->m_pVolumetricDataTimeDev->m_pBin);
	}
	
	//memset(pdeSolution, 0, res[0] * res[1] * res[2] * sizeof(double));
	
	//mprintf("solution: %g\ntimedev: %g\ndensity: %g\ngradx: %g\ngrady: %g\ngradz: %g\n",pdeSolution[0],ts2->m_pVolumetricDataTimeDev->m_pBin[0],m_pDensity.m_pBin[0],m_pDensityGrad[0].m_pBin[0],m_pDensityGrad[1].m_pBin[0],m_pDensityGrad[2].m_pBin[0]);

	double thresh2 = thresh;

	mprintf("[");

	if (!m_bPDERestart)
		goto _norestart;

	mprintf("R");

	if (m_fPDEInfo != NULL)
		fprintf(m_fPDEInfo,"Starting PDE solver...\n");


	if (!CCurrentPDESolver::bicgstabl(4, &pdeMatrix, m_pPDESolution, ts2->m_pVolumetricDataTimeDev->m_pBin, g_iPDEMaxIter, &thresh2, m_fPDEInfo, m_bPDEFastMode )) {

		mprintf("|");
_norestart:
		c = true;
		memset(m_pPDESolution, 0, res[0] * res[1] * res[2] * sizeof(double));
		mprintf("Z");

		if (m_fPDEInfo != NULL)
			fprintf(m_fPDEInfo,"Starting PDE solver with zero initial guess...\n");

		thresh2 = thresh;

		if (!CCurrentPDESolver::bicgstabl(4, &pdeMatrix, m_pPDESolution, ts2->m_pVolumetricDataTimeDev->m_pBin, g_iPDEMaxIter, &thresh2, m_fPDEInfo, m_bPDEFastMode )) {

			c = false;

			if (!m_bPDEFastMode) {

				memset(m_pPDESolution, 0, res[0] * res[1] * res[2] * sizeof(double));
	//			if (verbose)
	//				mprintf("@PDE Trying different threshold\n");


				mprintf("|T");
				thresh2 *= 1.05;
			//	mprintf("\nSetting threshold to %f.\n",thresh2);

				if (m_fPDEInfo != NULL)
					fprintf(m_fPDEInfo,"Starting PDE solver with reduced threshold (%f)...\n",thresh2);

				CCurrentPDESolver::bicgstabl(4, &pdeMatrix, m_pPDESolution, ts2->m_pVolumetricDataTimeDev->m_pBin, g_iPDEMaxIter, &thresh2, m_fPDEInfo, m_bPDEFastMode );
			}
		}
	}

	mprintf("]");

	if (c)
		mprintf("  Converged.");
	else
		mprintf("  Not converged (%.2f).",thresh2/thresh);
	
	ts2->m_pCurrentDensity->SetSize(3 * res[0] * res[1] * res[2]);
	
	for (i = 0; i < res[2]; i++) {
		for (j = 0; j < res[1]; j++) {
			for (k = 1; k < res[0] - 1; k++)
				ts2->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3) = (m_pPDESolution[i * res[1] * res[0] + j * res[0] + k + 1] - m_pPDESolution[i * res[1] * res[0] + j * res[0] + k - 1]) / (2.0 * g_fCubeXStep) * ts2->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + k];
			ts2->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3) = (m_pPDESolution[i * res[1] * res[0] + j * res[0] + 1] - m_pPDESolution[i * res[1] * res[0] + j * res[0] + res[0] - 1]) / (2.0 * g_fCubeXStep) * ts2->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0]];
			ts2->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + (res[0] - 1) * 3) = (m_pPDESolution[i * res[1] * res[0] + j * res[0]] - m_pPDESolution[i * res[1] * res[0] + j * res[0] + res[0] - 2]) / (2.0 * g_fCubeXStep) * ts2->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + res[0] - 1];
		}
	}

	for (i = 0; i < res[2]; i++) {
		for (j = 1; j < res[1] - 1; j++)
			for (k = 0; k < res[0]; k++)
				ts2->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 1) = (m_pPDESolution[i * res[1] * res[0] + (j + 1) * res[0] + k] - m_pPDESolution[i * res[1] * res[0] + (j - 1) * res[0] + k]) / (2.0 * g_fCubeYStep) * ts2->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + k];
		for (k = 0; k < res[0]; k++) {
			ts2->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + k * 3 + 1) = (m_pPDESolution[i * res[1] * res[0] + res[0] + k] - m_pPDESolution[i * res[1] * res[0] + (res[1] - 1) * res[0] + k]) / (2.0 * g_fCubeYStep) * ts2->m_pVolumetricData->m_pBin[i * res[1] * res[0] + k];
			ts2->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + (res[1] - 1) * res[0] * 3 + k * 3 + 1) = (m_pPDESolution[i * res[1] * res[0] + k] - m_pPDESolution[i * res[1] * res[0] + (res[1] - 2) * res[0] + k]) / (2.0 * g_fCubeYStep) * ts2->m_pVolumetricData->m_pBin[i * res[1] * res[0] + (res[1] - 1) * res[0] + k];
		}
	}

	for (i = 1; i < res[2] - 1; i++)
		for (j = 0; j < res[1]; j++)
			for (k = 0; k < res[0]; k++)
				ts2->m_pCurrentDensity->GetAt(i * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 2) = (m_pPDESolution[(i + 1) * res[1] * res[0] + j * res[0] + k] - m_pPDESolution[(i - 1) * res[1] * res[0] + j * res[0] + k]) / (2.0 * g_fCubeZStep) * ts2->m_pVolumetricData->m_pBin[i * res[1] * res[0] + j * res[0] + k];
	for (j = 0; j < res[1]; j++) {
		for (k = 0; k < res[0]; k++) {
			ts2->m_pCurrentDensity->GetAt(j * res[0] * 3 + k * 3 + 2) = (m_pPDESolution[res[1] * res[0] + j * res[0] + k] - m_pPDESolution[(res[2] - 1) * res[1] * res[0] + j * res[0] + k]) / (2.0 * g_fCubeZStep) * ts2->m_pVolumetricData->m_pBin[j * res[0] + k];
			ts2->m_pCurrentDensity->GetAt((res[2] - 1) * res[1] * res[0] * 3 + j * res[0] * 3 + k * 3 + 2) = (m_pPDESolution[j * res[0] + k] - m_pPDESolution[(res[2] - 2) * res[1] * res[0] + j * res[0] + k]) / (2.0 * g_fCubeZStep) * ts2->m_pVolumetricData->m_pBin[(res[2] - 1) * res[1] * res[0] + j * res[0] + k];
		}
	}

	return true;
}




void CFitParabola::Init4(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {

	CxDMatrix3 a, t;
	double sx1, sx2, sx3, sx4, sy1, sy2, sy3, d;

	sx1 = x1          + x2          + x3          + x4;
	sx2 = x1*x1       + x2*x2       + x3*x3       + x4*x4;
	sx3 = x1*x1*x1    + x2*x2*x2    + x3*x3*x3    + x4*x4*x4;
	sx4 = x1*x1*x1*x1 + x2*x2*x2*x2 + x3*x3*x3*x3 + x4*x4*x4*x4;

	sy1 = y1          + y2          + y3          + y4;
	sy2 = x1*y1       + x2*y2       + x3*y3       + x4*y4;
	sy3 = x1*x1*y1    + x2*x2*y2    + x3*x3*y3    + x4*x4*y4;

	a(0,0) = 4.0;
	a(0,1) = sx1;
	a(0,2) = sx2;
	a(1,0) = sx1;
	a(1,1) = sx2;
	a(1,2) = sx3;
	a(2,0) = sx2;
	a(2,1) = sx3;
	a(2,2) = sx4;

	d = a.Det();

	if (d == 0) {
		eprintf("CParabola::Init4(): Error: Det == 0.\n");
		m_fA = 0;
		m_fB = 0;
		m_fC = 0;
		return;
	}

	t = a;
	t(0,0) = sy1;
	t(1,0) = sy2;
	t(2,0) = sy3;
	m_fA = t.Det() / d;

	t = a;
	t(0,1) = sy1;
	t(1,1) = sy2;
	t(2,1) = sy3;
	m_fB = t.Det() / d;

	t = a;
	t(0,2) = sy1;
	t(1,2) = sy2;
	t(2,2) = sy3;
	m_fC = t.Det() / d;

	//mprintf("%.10G + %.10G * x + %.10G * x^2\n",m_fA,m_fB,m_fC);
}


double CFitParabola::Evaluate(double x) {

	return m_fA + x * m_fB + x * x * m_fC;
}




