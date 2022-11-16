/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

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

#include "ir.h"

#include "globalvar.h"
#include "maintools.h"
#include "timestep.h"
#include "xstring.h"
#include "linalg.h"

#include <math.h>


const char *GetRevisionInfo_ir(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_ir() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


#define BUF_SIZE 4096


static CxObArray g_PowerObserv;
// static CxFloatArray g_TemperatureCache;


CPowerObservation::CPowerObservation(bool global) {

	int i, size;
	CxString buf, buf2;

	
	try { _atoms = new CAtomGroup(); } catch(...) { _atoms = NULL; }
	if (_atoms == NULL) NewException((double)sizeof(CAtomGroup), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	
	if (global) {
		m_iShowMol = -1;
		m_iShowMolCount = g_oaSingleMolecules.GetSize();
		_name = new char[7];
		sprintf(_name, "global");
	} else {
//		char buf[BUF_SIZE];
//		char buf2[BUF_SIZE];
//		size_t remaining = BUF_SIZE;
		if (g_oaMolecules.GetSize() > 1) {
/*#ifdef TARGET_LINUX
			remaining -= snprintf(buf, remaining, "    Which molecule should be observed (");
#else
			remaining -= sprintf(buf, "    Which molecule should be observed (");
#endif*/
			buf.sprintf( "    Which molecule should be observed (");
			for (i = 0; i < g_oaMolecules.GetSize(); i++) {

/*				if(remaining < 1)
					break;
#ifdef TARGET_LINUX
				size_t length = snprintf(buf2, remaining, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#else
				size_t length = sprintf(buf2, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#endif
				strncat(buf, buf2, remaining - 1);
				remaining -= length;*/

				buf2.sprintf("%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
				buf.strcat(buf2);

				if (i < g_oaMolecules.GetSize() - 1) {
/*#ifdef TARGET_LINUX
					length = snprintf(buf2, remaining, ", ");
#else
					length = sprintf(buf2, ", ");
#endif
					strncat(buf, buf2, remaining - 1);
					remaining -= length;*/

					buf2.sprintf(", ");
					buf.strcat(buf2);
				}
			}
//			strncat(buf, ")? ", remaining - 1);
			buf.strcat(")? ");
			m_iShowMol = AskRangeInteger_ND("%s", 1, g_oaMolecules.GetSize(),(const char*)buf) - 1;
		} else {
			m_iShowMol = 0;
			mprintf("    Observing molecule %s.\n", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
		}
		
		while(true) {
			mprintf("    Which atom(s) to observe (e.g. \"C1,C3-5,H\", \"*\"=all)? [*] ");
			inpprintf("! Which atom(s) to observe (e.g. \"C1,C3-5,H\", \"*\"=all)? [*]\n");
//			char buf[BUF_SIZE];
			myget(&buf);
			if(strlen(buf) == 0) {
				if(!_atoms->ParseAtoms((CMolecule *)g_oaMolecules[m_iShowMol], "*")) {
					eprintf("CPowerObservation::CPowerObservation(): Internal error 1.\n");
					continue;
				}
			} else if(!_atoms->ParseAtoms((CMolecule *)g_oaMolecules[m_iShowMol], buf)) {
				continue;
			}
			break;
		}
		mprintf("\n    Observing %d atoms of molecule %s.\n", _atoms->m_iAtomGes, ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
		m_iShowMolCount = ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
		_name = new char[strlen(((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName) + strlen(_atoms->m_sName) + 4];
		sprintf(_name, "[%s_%s]", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName, _atoms->m_sName);
		mprintf("\n");
	}
	
	if(g_iTrajSteps != -1) {
		_correlationDepth = (int)(0.75 * g_iTrajSteps);
		if (g_bPower) {
			if (_correlationDepth > 4096)
				_correlationDepth = 4096;
			if ((g_fTimestepLength > 1.0) && (_correlationDepth > 2048))
				_correlationDepth = 2048;
			if ((g_fTimestepLength > 2.0) && (_correlationDepth > 1024))
				_correlationDepth = 1024;
		}
		_correlationDepth = AskUnsignedInteger("    Enter the correlation depth of the ACF (in trajectory frames): [%d] ", _correlationDepth, _correlationDepth);
	} else {
		_correlationDepth = AskUnsignedInteger("    Enter the correlation depth of the ACF (in trajectory frames): [256] ", 256);
	}
	size = CalcFFTSize(_correlationDepth, false);
	if (_correlationDepth != size) {
		mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using this instead of %d as depth.\n\n", size, _correlationDepth);
		_correlationDepth = size;
	}

	m_iFinalDepth = _correlationDepth;

	if (g_bPower) {
		if (g_bAdvanced2) {
			m_bSplitCart = AskYesNo("    Write Cartesian contributions (x, y, z, xy, xz, yz) of this power spectrum (y/n)? [no] ",false);
			mprintf("\n");
		} else
			m_bSplitCart = false;
	} else {
		m_bSplitCart = AskYesNo("    Write Cartesian contributions (x, y, z, xy, xz, yz) of this autocorrelation (y/n)? [no] ",false);
		mprintf("\n");
	}
	
	if(g_bAdvanced2) {
		if (g_bPower)
			_windowFunction = AskRangeInteger("    Window function: None (0), cos^2(a*t) (1), exp(-t/a) (2), exp(-(t/a)^2) (3) [1] ", 0, 3, 1);
		else
			_windowFunction = AskRangeInteger("    Window function: None (0), cos^2(a*t) (1), exp(-t/a) (2), exp(-(t/a)^2) (3) [0] ", 0, 3, 0);
		if(_windowFunction == 1) {
			mprintf("    The parameter \"a\" is chosen according to the correlation depth.\n");
			_windowFunctionParameter = 0;
		} else if(_windowFunction == 2) {
			_windowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in trajectory frames): [%d] ", m_iFinalDepth / 4, m_iFinalDepth / 4);
		} else if(_windowFunction == 3) {
			_windowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in trajectory frames): [%d] ", m_iFinalDepth / 2, m_iFinalDepth / 2);
		} else if (_windowFunction != 0) {
			eprintf("CPowerObservation::CPowerObservation(): Internal error 2.\n");
			abort();
		}
		mprintf("\n");
	} else {
		if (g_bPower)
			_windowFunction = 1;
		else
			_windowFunction = 0;
		_windowFunctionParameter = 0;
	}
	
	if (g_bPower)
		_massWeighting = AskYesNo("    Weight autocorrelation functions by atomic mass (y/n)? [yes] ", true);
	else
		_massWeighting = AskYesNo("    Weight autocorrelation functions by atomic mass (y/n)? [no] ", false);
	mprintf("\n");
		
	if (g_bPower) {
		if(g_bAdvanced2) {
			_zeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ", m_iFinalDepth * 3, m_iFinalDepth * 3);
			size = CalcFFTSize(m_iFinalDepth + _zeroPadding, false);
			if(m_iFinalDepth + _zeroPadding != size) {
				mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-m_iFinalDepth);
				_zeroPadding = size - m_iFinalDepth;
			}
		} else {
			_zeroPadding = m_iFinalDepth * 3;
			size = CalcFFTSize(m_iFinalDepth + _zeroPadding, false);
			mprintf("    Using cos^2 window function; inserting %d zeros for zero padding.\n",_zeroPadding);
			if(m_iFinalDepth + _zeroPadding != size) {
				mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-m_iFinalDepth);
				_zeroPadding = size - m_iFinalDepth;
			}
		}
		
		double possibleRange = 33356.41 / g_fTimestepLength / 2.0;
		_specResolution = possibleRange / (m_iFinalDepth + _zeroPadding);
		mprintf("    This results in a spectral resolution of %.2f cm^-1.\n", _specResolution);
		mprintf("\n    A time step length of %.2f fs allows a spectral range up to %.2f cm^-1.\n", g_fTimestepLength, possibleRange);
		double specLimit = AskRangeFloat("\n    Calculate spectrum up to which wave number (cm^-1)? [%.2f cm^-1] ", 0, possibleRange, (possibleRange < 5000.0) ? possibleRange : 5000.0, (possibleRange < 5000.0) ? possibleRange : 5000.0);
		_specSize = (int)(specLimit / _specResolution);
		mprintf("\n");
		
		if(g_bAdvanced2) {
			_saveACF = AskYesNo("    Save autocorrelation function (y/n)? [no] ", false);
			mprintf("\n");
		} else
			_saveACF = false;

		_correctfreq = AskYesNo("    Correct frequency shift of the Verlet integrator (y/n)? [no] ",false);
		if (_correctfreq) {
			ParseCorrectWavenumber();
			mprintf("\n    Due to the frequency correction, the spectral range is now up to %.2f cm^-1.\n",CorrectWavenumber(specLimit));
		}
		mprintf("\n");


	} else {

		_saveACF = true;
		_correctfreq = false;
	}
}


CPowerObservation::~CPowerObservation() {
	delete[] _name;
	delete _atoms;
}


void CPowerObservation::initialize() {
	int i, j, k;
	int n;
	if(g_iTrajSteps != -1)
		n = (int)(1.1 * g_iTrajSteps / g_iStride);
	else
		n = 10000;
	
	if(m_iShowMol == -1) {
		mprintf("    Velocity cache: Trying to allocate %s of memory...\n", FormatBytes((double)g_iGesAtomCount * n * sizeof(CxDVector3)));
		for(i = 0; i < g_iGesAtomCount; i++) {
			CxDVec3Array *a;
			try { a = new CxDVec3Array(); } catch(...) { a = NULL; }
			if(a == NULL) NewException((double)sizeof(CxDVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			a->SetMaxSize(n);
			a->SetGrow(n / 10);
			_velocityCache.Add(a);
		}
	} else {
		mprintf("    Velocity cache: Trying to allocate %s of memory...\n", FormatBytes((double)m_iShowMolCount * _atoms->m_iAtomGes * n * sizeof(CxDVector3)));
		for(i = 0; i < m_iShowMolCount * _atoms->m_iAtomGes; i++) {
			CxDVec3Array *a;
			try { a = new CxDVec3Array(); } catch(...) { a = NULL; }
			if(a == NULL) NewException((double)sizeof(CxDVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
			a->SetMaxSize(n);
			a->SetGrow(n / 10);
			_velocityCache.Add(a);
		}
	}
	
	if(m_iShowMol == -1) {
		_masses.SetSize(g_iGesAtomCount);
		if(_massWeighting) {
			for(i = 0; i < g_iGesAtomCount; i++) {
				_masses[i] = ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_pElement->m_fMass;
				if (_masses[i] == 0) {
					eprintf("\nWarning: Zero mass for element %s, using 1 for mass weighting.",(const char*)((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_sName);
					_masses[i] = 1.0;
				}
			}
		} else {
			for(i = 0; i < g_iGesAtomCount; i++) {
				_masses[i] = 1.0;
			}
		}
	} else {
		_masses.SetSize(m_iShowMolCount * _atoms->m_iAtomGes);
		if(_massWeighting) {
			for(i = 0; i < m_iShowMolCount; i++) {
				n = 0;
				for(j = 0; j < _atoms->m_baRealAtomType.GetSize(); j++) {
					if ((i == 0) && (((CAtom *)g_oaAtoms[_atoms->m_baRealAtomType[j]])->m_pElement->m_fMass == 0))
						eprintf("\nWarning: Zero mass for element %s, using 1 for mass weighting.",(const char*)((CAtom *)g_oaAtoms[_atoms->m_baRealAtomType[j]])->m_sName);
					CxIntArray *a = (CxIntArray *)_atoms->m_oaAtoms[j];
					for(k = 0; k < a->GetSize(); k++) {
						_masses[i * _atoms->m_iAtomGes + n] = ((CAtom *)g_oaAtoms[_atoms->m_baRealAtomType[j]])->m_pElement->m_fMass;
						if (_masses[i * _atoms->m_iAtomGes + n] == 0)
							_masses[i * _atoms->m_iAtomGes + n] = 1.0;
						n++;
					}
				}
			}
		} else {
			for(i = 0; i < m_iShowMolCount * _atoms->m_iAtomGes; i++) {
				_masses[i] = 1.0;
			}
		}
	}
}


void CPowerObservation::process(CTimeStep *ts) {
	int i, j, k;
	if(m_iShowMol == -1) {
		for(i = 0; i < g_iGesAtomCount; i++) {
			((CxDVec3Array *)_velocityCache[i])->Add(ts->m_vaVelocities[i]);
		}
	} else {
		for(i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
			int n = 0;
			for(j = 0; j < _atoms->m_baAtomType.GetSize(); j++) {
				CxIntArray *a = (CxIntArray *)_atoms->m_oaAtoms[j];
				for(k = 0; k < a->GetSize(); k++) {
					((CxDVec3Array *)_velocityCache[i * _atoms->m_iAtomGes + n])->Add(ts->m_vaVelocities[((CxIntArray *)sm->m_oaAtomOffset[_atoms->m_baAtomType[j]])->GetAt(a->GetAt(k))]);
					n++;
				}
			}
		}
	}
}


void CPowerObservation::finalize() {

	int i, j, k, l, n, oldSize, split, cc=0;
	double step, possibleRange, integral, fac;
	CFFT *fft;
	CxDoubleArray *acf, *temp, *temp2, *spectrum;
	CAutoCorrelation *ac;
	FILE *specFile, *acfFile;
	char filename[BUF_SIZE];
	const char *contrib=NULL;


	split = 0;

	n = ((CxDVec3Array *)_velocityCache[0])->GetSize();
	if (n < _correlationDepth) {
		eprintf("\nError: Autocorrelation depth is %d, but only %d timesteps evaluated.\n",_correlationDepth,n);
		eprintf("       Reduce depth or increase trajectory length.\n\n");
		abort();
	}

	if ((g_iStride != 1) && g_bPower) {
		mprintf("\n    A trajectory stride of %d was entered, correcting spectral resolution...\n",g_iStride);
		possibleRange = 33356.41 / g_fTimestepLength / g_iStride / 2.0;
		_specResolution = possibleRange / (_correlationDepth + _zeroPadding);
		mprintf("    New spectral resolution is %.2f cm^-1.\n\n", _specResolution);
		_specSize *= g_iStride;
		if (_specSize > _correlationDepth+_zeroPadding)
			_specSize = _correlationDepth+_zeroPadding;
	}

	if(m_iShowMol == -1)
		step = (double)g_iGesAtomCount / 50.0;
	else
		step = (double)m_iShowMolCount * _atoms->m_iAtomGes / 50.0;
	
_beginsplit:
	switch(split) {
		case 0:
			cc = 3;
			contrib = "";
			if (m_bSplitCart)
				mprintf(WHITE,"  *** All Contributions ***\n");
			break;
		case 1:
			cc = 1;
			contrib = "_x";
			if (m_bSplitCart)
				mprintf(WHITE,"\n  *** X Projection ***\n");
			break;
		case 2:
			cc = 1;
			contrib = "_y";
			if (m_bSplitCart)
				mprintf(WHITE,"\n  *** Y Projection ***\n");
			break;
		case 3:
			cc = 1;
			contrib = "_z";
			if (m_bSplitCart)
				mprintf(WHITE,"\n  *** Z Projection ***\n");
			break;
		case 4:
			cc = 2;
			contrib = "_xy";
			if (m_bSplitCart)
				mprintf(WHITE,"\n  *** XY Projection ***\n");
			break;
		case 5:
			cc = 2;
			contrib = "_xz";
			if (m_bSplitCart)
				mprintf(WHITE,"\n  *** XZ Projection ***\n");
			break;
		case 6:
			cc = 2;
			contrib = "_yz";
			if (m_bSplitCart)
				mprintf(WHITE,"\n  *** YZ Projection ***\n");
			break;
	}

	mprintf("    Calculating autocorrelation...\n");
	try { acf = new CxDoubleArray(); } catch(...) { acf = NULL; }
	if(acf == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	acf->SetSize(_correlationDepth);
	for(i = 0; i < _correlationDepth; i++)
		acf->GetAt(i) = 0.0;
	try { ac = new CAutoCorrelation(); } catch(...) { ac = NULL; }
	if(ac == NULL) NewException((double)sizeof(CAutoCorrelation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	ac->Init(n, _correlationDepth, g_bACFFFT);
	
	try { temp = new CxDoubleArray(); } catch(...) { temp = NULL; }
	if(temp == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp->SetSize(n);
	try { temp2 = new CxDoubleArray(); } catch(...) { temp2 = NULL; }
	if(temp2 == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp2->SetSize(_correlationDepth);
	
	mprintf(WHITE, "      [");
	if (m_iShowMol == -1) {
		for (i = 0; i < g_iGesAtomCount; i++) {
			if (fmod(i, step) < 1.0)
				mprintf(WHITE, "#");
			for (j = 0; j < 3; j++) {
				switch(split) {
					case 0: break;
					case 1: if (j != 0) continue; break;
					case 2: if (j != 1) continue; break;
					case 3: if (j != 2) continue; break;
					case 4: if (j == 2) continue; break;
					case 5: if (j == 1) continue; break;
					case 6: if (j == 0) continue; break;
				}
				for(k = 0; k < n; k++)
					temp->GetAt(k) = ((CxDVec3Array *)_velocityCache[i])->GetAt(k)[j];
				ac->AutoCorrelate(temp, temp2);
				for(k = 0; k < _correlationDepth; k++)
					acf->GetAt(k) += temp2->GetAt(k) * _masses[i];
			}
		}
	} else {
		for(i = 0; i < m_iShowMolCount; i++) {
			for(j = 0; j < _atoms->m_iAtomGes; j++) {
				if(fmod(i * _atoms->m_iAtomGes + j, step) < 1.0)
					mprintf(WHITE, "#");
				for(k = 0; k < 3; k++) {
					switch(split) {
						case 0: break;
						case 1: if (k != 0) continue; break;
						case 2: if (k != 1) continue; break;
						case 3: if (k != 2) continue; break;
						case 4: if (k == 2) continue; break;
						case 5: if (k == 1) continue; break;
						case 6: if (k == 0) continue; break;
					}
					for(l = 0; l < n; l++)
						temp->GetAt(l) = ((CxDVec3Array *)_velocityCache[i * _atoms->m_iAtomGes + j])->GetAt(l)[k];
					ac->AutoCorrelate(temp, temp2);
					for(l = 0; l < _correlationDepth; l++)
						acf->GetAt(l) += temp2->GetAt(l) * _masses[i * _atoms->m_iAtomGes + j];
				}
			}
		}
	}
	mprintf(WHITE, "]\n");
	
	if (m_iShowMol != -1)
		for (i = 0; i < _correlationDepth; i++)
			acf->GetAt(i) /= (double)m_iShowMolCount;
	
	delete ac;
	delete temp2;
	
	if (_saveACF) {
		if (g_bPower) {
			#ifdef TARGET_LINUX
				snprintf(filename, BUF_SIZE, "power_acf%s_%s.csv", contrib, _name);
			#else
				sprintf(filename, "power_acf%s_%s.csv", contrib, _name);
			#endif
		} else {
			#ifdef TARGET_LINUX
				snprintf(filename, BUF_SIZE, "vacf%s_%s.csv", contrib, _name);
			#else
				sprintf(filename, "vacf%s_%s.csv", contrib, _name);
			#endif
		}
		mprintf("    Saving autocorrelation function as \"%s\"...\n", filename);
		acfFile = OpenFileWrite(filename, false);
		integral = 0;
		fprintf(acfFile, "#Correlation Depth / fs;  Velocity ACF / pm^2 ps^-2;  Integral / pm^2 ps^-1\n");
		for(i = 0; i < _correlationDepth; i++) {
			integral += acf->GetAt(i) * g_fTimestepLength / 1000.0;
			fprintf(acfFile, "%.2f; %.8G; %.8G\n", i * g_fTimestepLength, acf->GetAt(i), integral);
		}
		fclose(acfFile);
	}

	
	temp->CopyFrom(acf);
	
	if (_windowFunction == 1) {
		for(i = 0; i < temp->GetSize(); i++)
			temp->GetAt(i) *= pow2(cos((double)i / (temp->GetSize() - 1) / 2.0 * Pi));
	} else if(_windowFunction == 2) {
		for(i = 0; i < temp->GetSize(); i++)
			temp->GetAt(i) *= exp(-(double)i / _windowFunctionParameter);
	} else if(_windowFunction == 3) {
		for(i = 0; i < temp->GetSize(); i++)
			temp->GetAt(i) *= exp(-(double)i * i / _windowFunctionParameter / _windowFunctionParameter);
	} else if(_windowFunction != 0) {
		eprintf("Unknown window function.\n");
		abort();
	}
	
	if (_saveACF && (_windowFunction != 0)) {
		if (g_bPower) {
			#ifdef TARGET_LINUX
				snprintf(filename, BUF_SIZE, "power_acf%s_windowed_%s.csv", contrib, _name);
			#else
				sprintf(filename, "power_acf%s_windowed_%s.csv", contrib, _name);
			#endif
		} else {
			#ifdef TARGET_LINUX
				snprintf(filename, BUF_SIZE, "vacf%s_windowed_%s.csv", contrib, _name);
			#else
				sprintf(filename, "vacf%s_windowed_%s.csv", contrib, _name);
			#endif
		}
		mprintf("    Saving windowed autocorrelation function as \"%s\"...\n", filename);
		acfFile = OpenFileWrite(filename, false);
		fprintf(acfFile, "#Correlation Depth / fs; Windowed Velocity ACF / pm^2 ps^-2\n");
		for(i = 0; i < m_iFinalDepth; i++)
			fprintf(acfFile, "%.2f; %.8G\n", i * g_fTimestepLength * g_iStride, temp->GetAt(i));
		fclose(acfFile);
	}
	
	if (g_bPower) {


		if (_zeroPadding > 0)
			for(i = 0; i < _zeroPadding; i++)
				temp->Add(0.0);
		
		oldSize = temp->GetSize();
		temp->SetSize(2 * oldSize);
		for (i = 1; i < oldSize; i++)
			temp->GetAt(oldSize + i) = temp->GetAt(oldSize - i);
		temp->GetAt(oldSize) = 0.0;
		
		mprintf("    Fourier transforming autocorrelation function...\n");
		try { fft = new CFFT(); } catch(...) { fft = NULL; }
		if (fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		fft->PrepareFFT_C2C(temp->GetSize());
		for (i = 0; i < temp->GetSize(); i++) {
			fft->m_pInput[2*i] = temp->GetAt(i);
			fft->m_pInput[2*i+1] = 0.0;
		}
		fft->DoFFT();
		
		try { spectrum = new CxDoubleArray(); } catch(...) { spectrum = NULL; }
		if (spectrum == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		spectrum->SetSize(_specSize);
		if (m_iShowMol == -1)
			fac = 7.211349e-9 * g_fTimestepLength * g_iStride / (double)g_iGesAtomCount / (double)cc;
		else
			fac = 7.211349e-9 * g_fTimestepLength * g_iStride / (double)_atoms->m_iAtomGes / (double)cc;
		for (i = 0; i < _specSize; i++)
			// Fix MB 02.05.2021
			spectrum->GetAt(i) = fac * fft->m_pOutput[2*i]; // Output in K*cm
			//spectrum->GetAt(i) = 7.211349e-9 * fft->m_pOutput[2*i] * g_fTimestepLength * g_iStride; // Output in K*cm
		
		#ifdef TARGET_LINUX
			snprintf(filename, BUF_SIZE, "power_spectrum%s_%s.csv", contrib, _name);
		#else
			sprintf(filename, "power_spectrum%s_%s.csv", contrib, _name);
		#endif
			mprintf("    Saving spectrum as \"%s\"...\n", filename);
		specFile = OpenFileWrite(filename, false);

		if (_correctfreq)
			fprintf(specFile, "#Corrected Wavenumber (cm^-1); Spectrum (K*cm); Integral (K)\n");
		else
			fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (K*cm); Integral (K)\n");

		integral = 0.0;
		for (i = 0; i < _specSize; i++) {

			integral += (double)spectrum->GetAt(i) * _specResolution;

			if (_correctfreq)
				fprintf(specFile, "%.2f; %.8G; %.14G\n", CorrectWavenumber(_specResolution * i), spectrum->GetAt(i), integral);
			else
				fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum->GetAt(i), integral);
		}
		fclose(specFile);
		
		mprintf("\n");

		if (m_iShowMol == -1)
			mprintf("    Assuming %dn = %d degrees of freedom, the average temperature is %.3f K.\n", cc, cc * g_iGesAtomCount, integral );
		else
			mprintf("    Assuming %dn = %d degrees of freedom, the average temperature is %.3f K.\n", cc, cc * _atoms->m_iAtomGes, integral );
		
		delete fft;
		delete spectrum;
	}

	delete acf;
	delete temp;

	if (m_bSplitCart && (split < 6)) {
		split++;
		goto _beginsplit;
	}

	for(i = 0; i < _velocityCache.GetSize(); i++)
		delete (CxDVec3Array *)_velocityCache[i];
}
	

bool gatherPowerSpectrum() {

	bool tb;
	CPowerObservation *obs;
	const char *fname, *fnameu;


	if (g_bPower) {
		fname = "power spectrum";
		fnameu = "Power Spectrum";
	} else {
		fname = "velocity autocorrelation";
		fnameu = "Velocity Autocorrelation";
	}

	g_bUseVelocities = true;
	
	if (AskYesNo("\n    Compute global %s of whole system (y/n) [yes] ", true, fname)) {

		mprintf(YELLOW, "\n>>> Global %s >>>\n\n",fnameu);
		
		try { obs = new CPowerObservation(true); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CPowerObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_PowerObserv.Add(obs);
		
		mprintf(YELLOW, "<<< End of Global %s <<<\n\n",fnameu);
		tb = true;
	} else
		tb = false;

	if (tb)
		if (!AskYesNo("    Compute %s of certain atoms/molecules (y/n)? [no] ",false, fname))
			return true;
	
	while(true) {
		mprintf(YELLOW, "\n>>> %s %d >>>\n\n", fnameu, g_PowerObserv.GetSize() + 1);
		
		try { obs = new CPowerObservation(); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CPowerObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_PowerObserv.Add(obs);
		
		mprintf(YELLOW, "<<< End of %s %d <<<\n\n", fnameu, g_PowerObserv.GetSize());
		
		if(!AskYesNo("    Add another %s observation (y/n)? [no] ", false, fname))
			break;
		mprintf("\n");
	}
	
	return true;
}


bool initializePowerSpectrum() {
	int i;
	for(i = 0; i < g_PowerObserv.GetSize(); i++) {
		((CPowerObservation *)g_PowerObserv[i])->initialize();
	}
	
// 	int n;
// 	if(g_iTrajSteps != -1)
// 		n = (int)(1.1 * g_iTrajSteps / g_iStride);
// 	else
// 		n = 10000;
// 	
// 	mprintf("    Temperature cache: Trying to allocate %s of memory...\n", FormatBytes((double)g_iGesAtomCount * n * sizeof(double)));
// 	g_TemperatureCache.SetMaxSize(n);
// 	g_TemperatureCache.SetGrow(n / 10);
	
	return true;
}


void processPowerSpectrum(CTimeStep* ts) {
	int i;
	for(i = 0; i < g_PowerObserv.GetSize(); i++) {
		((CPowerObservation *)g_PowerObserv[i])->process(ts);
	}
	
// 	double temperature = 0.0f;
// 	for (i = 0; i < g_iGesAtomCount; i++) {
// 		temperature += ts->m_vaVelocities[i].GetLengthSqr() * ((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_pElement->m_fMass;
// 	}
// 	temperature /= (3.0f * g_iGesAtomCount * 8314.4621f);
// 	g_TemperatureCache.Add(temperature);
}


void finalizePowerSpectrum() {

	int i;
	const char *fnameu;


	if (g_bPower)
		fnameu = "Power Spectrum";
	else
		fnameu = "Velocity Autocorrelation";

	for(i = 0; i < g_PowerObserv.GetSize(); i++) {
		mprintf(YELLOW, "\n>>> %s %d >>>\n\n", fnameu, i+1);
		((CPowerObservation *)g_PowerObserv[i])->finalize();
		delete (CPowerObservation *)g_PowerObserv[i];
		mprintf(YELLOW, "\n<<< End of %s %d <<<\n\n", fnameu, i+1);
	}
	
// 	CxString filename;
// 	filename.sprintf("temperature.csv");
// 	mprintf("    Saving temperature as %s...\n", (const char *)filename);
// 	FILE *tempFile = OpenFileWrite(filename, false);
// 	fprintf(tempFile, "#Time (fs); Temperature (K)\n");
// 	double average = 0.0;
// 	for(i = 0; i < g_TemperatureCache.GetSize(); i++) {
// 		average += (double)g_TemperatureCache[i];
// 		fprintf(tempFile, "%.2f; %.8G\n", g_fTimestepLength * i, g_TemperatureCache[i]);
// 	}
// 	fclose(tempFile);
// 	average /= g_TemperatureCache.GetSize();
// 	
// 	mprintf("\n    The average temperature is %.2f K\n", average);
}


static CxObArray g_IRObserv;


CIRObservation::CIRObservation(bool global) {
	int i;
	if(global) {
		m_iShowMol = -1;
		m_iShowMolCount = g_oaSingleMolecules.GetSize();
		_name = new char[7];
		sprintf(_name, "global");
	} else {
		char buf[BUF_SIZE];
		char buf2[BUF_SIZE];
		size_t remaining = BUF_SIZE;
		if(g_oaMolecules.GetSize() > 1) {
#ifdef TARGET_LINUX
			remaining -= snprintf(buf, remaining, "    Which molecule should be observed (");
#else
			remaining -= sprintf(buf, "    Which molecule should be observed (");
#endif
			for(i = 0; i < g_oaMolecules.GetSize(); i++) {
				if(remaining < 1)
					break;
#ifdef TARGET_LINUX
				size_t length = snprintf(buf2, remaining, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#else
				size_t length = sprintf(buf2, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#endif
				strncat(buf, buf2, remaining - 1);
				remaining -= length;
				if(i < g_oaMolecules.GetSize() - 1) {
#ifdef TARGET_LINUX
					length = snprintf(buf2, remaining, ", ");
#else
					length = sprintf(buf2, ", ");
#endif
					strncat(buf, buf2, remaining - 1);
					remaining -= length;
				}
			}
			strncat(buf, ")? ", remaining - 1);
			m_iShowMol = AskRangeInteger_ND("%s", 1, g_oaMolecules.GetSize(),(const char*)buf) - 1;
		} else {
			m_iShowMol = 0;
			mprintf("    Observing molecule %s.\n", ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
		}
		m_iShowMolCount = ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
		_name = new char[strlen(((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName) + 1];
		strcpy(_name, ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
		mprintf("\n");
	}
	
	if(g_iTrajSteps != -1) {
		_correlationDepth = (int)(0.75 * g_iTrajSteps);
		if(_correlationDepth > 4096)
			_correlationDepth = 4096;
		if(g_fTimestepLength > 1.0)
			_correlationDepth = 2048;
		if(g_fTimestepLength > 2.0)
			_correlationDepth = 1024;
		_correlationDepth = AskUnsignedInteger("    Enter the correlation depth of the ACF (in trajectory frames): [%d] ", _correlationDepth, _correlationDepth);
	} else {
		_correlationDepth = AskUnsignedInteger("    Enter the correlation depth of the ACF (in trajectory frames): [256] ", 256);
	}
	int size = CalcFFTSize(_correlationDepth, false);
	if(_correlationDepth != size) {
		mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using this instead of %d as depth.\n", size, _correlationDepth);
		_correlationDepth = size;
	}
	
	if(g_bAdvanced2) {
		_windowFunction = AskRangeInteger("    Window function: cos^2(a*t) (1), exp(-t/a) (2), exp(-(t/a)^2) (3) [1] ", 1, 3, 1);
		if(_windowFunction == 1) {
			mprintf("    The parameter \"a\" is chosen according to the correlation depth.\n");
			_windowFunctionParameter = 0;
		} else if(_windowFunction == 2) {
			_windowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in trajectory frames): [%d] ", _correlationDepth / 4, _correlationDepth / 4);
		} else if(_windowFunction == 3) {
			_windowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in trajectory frames): [%d] ", _correlationDepth / 2, _correlationDepth / 2);
		} else {
			eprintf("This is impossible.\n");
			abort();
		}
	} else {
		_windowFunction = 1;
		_windowFunctionParameter = 0;
		mprintf("    Using cos^2 window function.\n");
	}
	
	if(g_bAdvanced2) {
		_zeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ", _correlationDepth * 3, _correlationDepth * 3);
		size = CalcFFTSize(_correlationDepth + _zeroPadding, false);
		if(_correlationDepth + _zeroPadding != size) {
			mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-_correlationDepth);
			_zeroPadding = size-_correlationDepth;
		}
	} else {
		_zeroPadding = _correlationDepth * 3;
		size = CalcFFTSize(_correlationDepth + _zeroPadding, false);
		if(_correlationDepth + _zeroPadding != size) {
			mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-_correlationDepth);
			_zeroPadding = size-_correlationDepth;
		}
		mprintf("    Inserting %d zeros for zero padding.\n",_zeroPadding);
	}
	
	double possibleRange = 33356.41 / g_fTimestepLength / 2.0;
	_specResolution = possibleRange / (_correlationDepth + _zeroPadding);
	mprintf("    This results in a spectral resolution of %.2f cm^-1.\n", _specResolution);
	mprintf("\n    A time step length of %.2f fs allows a spectral range up to %.2f cm^-1.\n", g_fTimestepLength, possibleRange);
	double specLimit = AskRangeFloat("\n    Calculate spectrum up to which wave number (cm^-1)? [%.2f cm^-1] ", 0, possibleRange, (possibleRange < 5000.0) ? possibleRange : 5000.0, (possibleRange < 5000.0) ? possibleRange : 5000.0);
	_specSize = (int)(specLimit / _specResolution);
	mprintf("\n");
	
	if(g_bAdvanced2) {
		_finiteDifferenceCorrection = AskYesNo("    Apply finite difference correction (y/n)? [yes] ", true);
		mprintf("\n");
	} else {
		_finiteDifferenceCorrection = true;
	}
	
	if(g_bAdvanced2) {
		_saveACF = AskYesNo("    Save autocorrelation function (y/n)? [no] ", false);
		mprintf("\n");
	} else {
		_saveACF = false;
	}
	
	if(g_bAdvanced2 && m_iShowMol == -1) {
		_includeCross = AskYesNo("    Include also cross-correlations (y/n)? [no] ", false);
		mprintf("\n");
	} else {
		_includeCross = false;
	}
	

	_correctfreq = AskYesNo("    Correct frequency shift of the Verlet integrator (y/n)? [no] ",false);
	if (_correctfreq) {
		ParseCorrectWavenumber();
		mprintf("\n    Due to the frequency correction, the spectral range is now up to %.2f cm^-1.\n\n",CorrectWavenumber(specLimit));
	}
}

CIRObservation::~CIRObservation() {
	delete[] _name;
}

void CIRObservation::initialize() {
	int n;
	if(g_iTrajSteps != -1)
		n = (int)(1.1 * g_iTrajSteps / g_iStride);
	else
		n = 10000;
	
	mprintf("    Moment cache: Trying to allocate %s of memory...\n", FormatBytes((double)m_iShowMolCount * n * sizeof(CxDVector3)));
	int i;
	for(i = 0; i < m_iShowMolCount; i++) {
		CxDVec3Array *a;
		try { a = new CxDVec3Array(); } catch(...) { a = NULL; }
		if(a == NULL) NewException((double)sizeof(CxDVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		a->SetMaxSize(n);
		a->SetGrow(n / 10);
		_dipoleCache.Add(a);
	}
}

void CIRObservation::process() {
	int i;
	if(m_iShowMol == -1) {
		for(i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
			((CxDVec3Array *)_dipoleCache[i])->Add(sm->m_vDipole);
		}
	} else {
		for(i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
			((CxDVec3Array *)_dipoleCache[i])->Add(sm->m_vDipole);
		}
	}
}

void CIRObservation::finalize() {
	int n = ((CxDVec3Array *)_dipoleCache[0])->GetSize() - 2;
	double step = (double)m_iShowMolCount / 20.0;

	if (n < _correlationDepth)
	{
		eprintf("\nError: Autocorrelation depth is %d, but only %d timesteps evaluated.\n",_correlationDepth,n);
		eprintf("       Reduce depth or increase trajectory length.\n\n");
		abort();
	}

	mprintf("    Deriving dipole moments...\n");
	mprintf(WHITE, "     [");
	int i, j, k, l;
	for(i = 0; i < m_iShowMolCount; i++) {
		if(fmod(i, step) < 1.0)
			mprintf(WHITE, "#");
		for(j = 0; j < n; j++) {
			((CxDVec3Array *)_dipoleCache[i])->GetAt(j) = 0.5 * (((CxDVec3Array *)_dipoleCache[i])->GetAt(j+2) - ((CxDVec3Array *)_dipoleCache[i])->GetAt(j)) / g_fTimestepLength;
		}
	}
	mprintf(WHITE, "]\n");
	
	mprintf("    Calculating autocorrelation...\n");
	CxDoubleArray *acf;
	try { acf = new CxDoubleArray(); } catch(...) { acf = NULL; }
	if(acf == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	acf->SetSize(_correlationDepth);
	for(i = 0; i < _correlationDepth; i++) {
		acf->GetAt(i) = 0.0;
	}
	CAutoCorrelation *ac;
	try { ac = new CAutoCorrelation(); } catch(...) { ac = NULL; }
	if(ac == NULL) NewException((double)sizeof(CAutoCorrelation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	ac->Init(n, _correlationDepth, g_bACFFFT);
	
	CxDoubleArray *temp;
	try { temp = new CxDoubleArray(); } catch(...) { temp = NULL; }
	if(temp == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp->SetSize(n);
	CxDoubleArray *temp2;
	try { temp2 = new CxDoubleArray(); } catch(...) { temp2 = NULL; }
	if(temp2 == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp2->SetSize(n);
	
	mprintf(WHITE, "     [");
	for(i = 0; i < m_iShowMolCount; i++) {
		if(fmod(i, step) < 1.0)
			mprintf(WHITE, "#");
		for(j = 0; j < 3; j++) {
			for(k = 0; k < n; k++) {
				temp->GetAt(k) = ((CxDVec3Array *)_dipoleCache[i])->GetAt(k)[j];
			}
			ac->AutoCorrelate(temp, temp2);
			for(k = 0; k < _correlationDepth; k++) {
				acf->GetAt(k) += temp2->GetAt(k);
			}
		}
	}
	mprintf(WHITE, "]\n");
// 	if(m_iShowMol == -1) {
// 		for(i = 0; i < _correlationDepth; i++) {
// 			acf->GetAt(i) /= 3.0f;
// 		}
// 	} else {
// 		for(i = 0; i < _correlationDepth; i++) {
// 			acf->GetAt(i) /= 3.0f * m_iShowMolCount;
// 		}
// 	}
// The 3.0f is included in the final normalization
	if(m_iShowMol != -1) {
		for(i = 0; i < _correlationDepth; i++) {
			acf->GetAt(i) /= (double)m_iShowMolCount;
		}
	}
	
	delete ac;
	
	CxDoubleArray *ccf = NULL;
	if(_includeCross) {
		mprintf("    Calculating cross-correlation...\n");
		try { ccf = new CxDoubleArray(); } catch(...) { ccf = NULL; }
		if(ccf == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		ccf->SetSize(_correlationDepth);
		for(i = 0; i < _correlationDepth; i++) {
			ccf->GetAt(i) = 0.0;
		}
		CCrossCorrelation *cc;
		try { cc = new CCrossCorrelation(); } catch(...) { cc = NULL; }
		if(cc == NULL) NewException((double)sizeof(CCrossCorrelation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		cc->Init(n, _correlationDepth, g_bACFFFT);
		CxDoubleArray *temp3;
		try { temp3 = new CxDoubleArray(); } catch(...) { temp3 = NULL; }
		if(temp3 == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		temp3->SetSize(n);
		temp->SetSize(n);
		temp2->SetSize(n);
		
		mprintf(WHITE, "     [");
		step = (double)m_iShowMolCount * (m_iShowMolCount - 1) / 20.0;
		int c = 0;
		for(i = 0; i < m_iShowMolCount; i++) {
			for(j = 0; j < m_iShowMolCount; j++) {
				if (i == j)
					continue;
				if(fmod((double)c++, step) < 1.0)
					mprintf(WHITE, "#");
				for(k = 0; k < 3; k++) {
					for(l = 0; l < n; l++) {
						temp->GetAt(l) = ((CxDVec3Array *)_dipoleCache[i])->GetAt(l)[k];
						temp2->GetAt(l) = ((CxDVec3Array *)_dipoleCache[j])->GetAt(l)[k];
					}
					cc->CrossCorrelate(temp, temp2, temp3);
					for(l = 0; l < _correlationDepth; l++) {
						ccf->GetAt(l) += temp3->GetAt(l);
					}
				}
			}
		}
		mprintf(WHITE, "]\n");
// 		for(i = 0; i < _correlationDepth; i++) {
// 			ccf->GetAt(i) /= 3.0f;
// 		}
		
		delete cc;
		delete temp3;
	}
	
	delete temp2;
	
	if(_saveACF) {
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "ir_acf_%s.csv", _name);
#else
		sprintf(filename, "ir_acf_%s.csv", _name);
#endif
		mprintf("    Saving autocorrelation function as %s...\n", filename);
		FILE *acfFile = OpenFileWrite(filename, false);
		for(i = 0; i < _correlationDepth; i++) {
			fprintf(acfFile, "%.2f; %.10G\n", i * g_fTimestepLength, acf->GetAt(i));
		}
		fclose(acfFile);
	}
	
	temp->CopyFrom(acf);
	
	if(_windowFunction == 1) {
		for(i = 0; i < temp->GetSize(); i++) {
			temp->GetAt(i) *= pow2(cos((double)i / (temp->GetSize() - 1) / 2.0 * Pi));
		}
	} else if(_windowFunction == 2) {
		for(i = 0; i < temp->GetSize(); i++) {
			temp->GetAt(i) *= exp(-(double)i / _windowFunctionParameter);
		}
	} else if(_windowFunction == 3) {
		for(i = 0; i < temp->GetSize(); i++) {
			temp->GetAt(i) *= exp(-(double)i * i / _windowFunctionParameter / _windowFunctionParameter);
		}
	} else if(_windowFunction != 0) {
		eprintf("Unknown window function.\n");
		abort();
	}
	
	if(_saveACF) {
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "ir_acf_windowed_%s.csv", _name);
#else
		sprintf(filename, "ir_acf_windowed_%s.csv", _name);
#endif
		mprintf("    Saving windowed autocorrelation function as %s...\n", filename);
		FILE *acfFile = OpenFileWrite(filename, false);
		for(i = 0; i < _correlationDepth; i++) {
			fprintf(acfFile, "%.2f; %.10G\n", i * g_fTimestepLength, temp->GetAt(i));
		}
		fclose(acfFile);
	}
	
	if(_zeroPadding > 0) {
		for(i = 0; i < _zeroPadding; i++) {
			temp->Add(0.0);
		}
	}
	
	int oldSize = temp->GetSize();
	temp->SetSize(2 * oldSize);
	for(i = 1; i < oldSize; i++) {
		temp->GetAt(oldSize + i) = temp->GetAt(oldSize - i);
	}
	temp->GetAt(oldSize) = 0.0;
	
	mprintf("    Fourier transforming autocorrelation function...\n");
	CFFT *fft;
	try { fft = new CFFT(); } catch(...) { fft = NULL; }
	if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	fft->PrepareFFT_C2C(temp->GetSize());
	for(i = 0; i < temp->GetSize(); i++) {
		fft->m_pInput[2*i] = temp->GetAt(i);
		fft->m_pInput[2*i+1] = 0.0;
	}
	fft->DoFFT();
	
	CxDoubleArray *spectrum;
	try { spectrum = new CxDoubleArray(); } catch(...) { spectrum = NULL; }
	if(spectrum == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	spectrum->SetSize(_specSize);
	for(i = 0; i < _specSize; i++) {
		spectrum->GetAt(i) = 3047.2310 * fft->m_pOutput[2*i] * g_fTimestepLength; // Output in K*cm*km/mol
	}
	
	if(_finiteDifferenceCorrection) {
		double f = _specResolution * g_fTimestepLength * 1.883652e-4;
		for(i = 1; i < _specSize; i++) {
			spectrum->GetAt(i) *= pow2(f * i / sin(f * i)); // Divide by sinc function to correct finite difference derivation
		}
	}
	
	
	if(_includeCross) {
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "ir_spectrum_auto_%s.csv", _name);
#else
		sprintf(filename, "ir_spectrum_auto_%s.csv", _name);
#endif
		mprintf("    Saving autocorrelation spectrum as %s...\n", filename);
		FILE *specFile = OpenFileWrite(filename, false);

		if (_correctfreq)
			fprintf(specFile, "#Corrected Wavenumber (cm^-1); Spectrum (K*cm*km*mol^-1); Integral (K*km*mol^-1)\n");
		else
			fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (K*cm*km*mol^-1); Integral (K*km*mol^-1)\n");

		double integral = 0.0;
		for(i = 0; i < _specSize; i++) {

			integral += (double)spectrum->GetAt(i) * _specResolution;

			if (_correctfreq)
				fprintf(specFile, "%.2f; %.8G; %.14G\n", CorrectWavenumber(_specResolution * i), spectrum->GetAt(i), integral);
			else
				fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum->GetAt(i), integral);
		}
		fclose(specFile);
		
		if(_saveACF) {
			//char filename[BUF_SIZE];
#ifdef TARGET_LINUX
			snprintf(filename, BUF_SIZE, "ir_ccf_%s.csv", _name);
#else
			sprintf(filename, "ir_ccf_%s.csv", _name);
#endif
			mprintf("    Saving cross-correlation function as %s...\n", filename);
			FILE *ccfFile = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(ccfFile, "%.2f; %.10G\n", i * g_fTimestepLength, ccf->GetAt(i));
			}
			fclose(ccfFile);
		}
		
		temp->CopyFrom(ccf);
		
		if(_windowFunction == 1) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= pow2(cos((double)i / temp->GetSize() / 2.0 * Pi));
			}
		} else if(_windowFunction == 2) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= exp(-(double)i / _windowFunctionParameter);
			}
		} else if(_windowFunction == 3) {
			for(i = 0; i < temp->GetSize(); i++) {
				temp->GetAt(i) *= exp(-(double)i * i / _windowFunctionParameter / _windowFunctionParameter);
			}
		} else if(_windowFunction != 0) {
			eprintf("Unknown window function.\n");
			abort();
		}
		
		if(_saveACF) {
			//char filename[BUF_SIZE];
#ifdef TARGET_LINUX
			snprintf(filename, BUF_SIZE, "ir_ccf_windowed_%s.csv", _name);
#else
			sprintf(filename, "ir_ccf_windowed_%s.csv", _name);
#endif
			mprintf("    Saving windowed cross-correlation function as %s...\n", filename);
			FILE *ccfFile = OpenFileWrite(filename, false);
			for(i = 0; i < _correlationDepth; i++) {
				fprintf(ccfFile, "%.2f; %.10G\n", i * g_fTimestepLength, ccf->GetAt(i));
			}
			fclose(ccfFile);
		}
		
		if(_zeroPadding > 0) {
			for(i = 0; i < _zeroPadding; i++) {
				temp->Add(0.0);
			}
		}
		
		//int oldSize = temp->GetSize();
		oldSize = temp->GetSize();
		temp->SetSize(2 * oldSize);
		for(i = 1; i < oldSize; i++) {
			temp->GetAt(oldSize + i) = temp->GetAt(oldSize - i);
		}
		temp->GetAt(oldSize) = 0.0;
		
		mprintf("    Fourier transforming cross-correlation function...\n");
		fft->PrepareFFT_C2C(temp->GetSize());
		for(i = 0; i < temp->GetSize(); i++) {
			fft->m_pInput[2*i] = temp->GetAt(i);
			fft->m_pInput[2*i+1] = 0.0;
		}
		fft->DoFFT();
		
		CxDoubleArray *ccSpectrum;
		try { ccSpectrum = new CxDoubleArray(); } catch(...) { ccSpectrum = NULL; }
		if(ccSpectrum == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		ccSpectrum->SetSize(_specSize);
		for(i = 0; i < _specSize; i++) {
			ccSpectrum->GetAt(i) = 3047.2310 * fft->m_pOutput[2*i] * g_fTimestepLength;
		}
		
		
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "ir_spectrum_cross_%s.csv", _name);
#else
		sprintf(filename, "ir_spectrum_cross_%s.csv", _name);
#endif
		mprintf("    Saving cross-correlation spectrum as %s...\n", filename);
		specFile = OpenFileWrite(filename, false);

		if (_correctfreq)
			fprintf(specFile, "#Corrected Wavenumber (cm^-1); Spectrum (K*cm*km*mol^-1); Integral (K*km*mol^-1)\n");
		else
			fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (K*cm*km*mol^-1); Integral (K*km*mol^-1)\n");

		integral = 0.0;
		for(i = 0; i < _specSize; i++) {

			integral += ccSpectrum->GetAt(i) * _specResolution;

			if (_correctfreq)
				fprintf(specFile, "%.2f; %.8G; %.14G\n", CorrectWavenumber(_specResolution * i), spectrum->GetAt(i), integral);
			else
				fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, ccSpectrum->GetAt(i), integral);
		}
		fclose(specFile);
		
		for(i = 0; i < _specSize; i++) {
			spectrum->GetAt(i) += ccSpectrum->GetAt(i);
		}
		
		delete ccSpectrum;
	}
	
	char filename[BUF_SIZE];
#ifdef TARGET_LINUX
	snprintf(filename, BUF_SIZE, "ir_spectrum_%s.csv", _name);
#else
	sprintf(filename, "ir_spectrum_%s.csv", _name);
#endif
	mprintf("    Saving spectrum as %s...\n", filename);
	FILE *specFile = OpenFileWrite(filename, false);

	if (_correctfreq)
		fprintf(specFile, "#Corrected Wavenumber (cm^-1); Spectrum (K*cm*km*mol^-1); Integral (K*km*mol^-1)\n");
	else
		fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (K*cm*km*mol^-1); Integral (K*km*mol^-1)\n");

	double integral = 0.0;
	for(i = 0; i < _specSize; i++) {

		integral += spectrum->GetAt(i) * _specResolution;

		if (_correctfreq)
			fprintf(specFile, "%.2f; %.8G; %.14G\n", CorrectWavenumber(_specResolution * i), spectrum->GetAt(i), integral);
		else
			fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum->GetAt(i), integral);
	}
	fclose(specFile);
	
	delete acf;
	if(_includeCross) {
		delete ccf;
	}
	delete temp;
	delete fft;
	delete spectrum;
	
	for(i = 0; i < _dipoleCache.GetSize(); i++)
		delete (CxDVec3Array *)_dipoleCache[i];
}


bool gatherIR() {
	g_bDipole = true;
	ParseDipole();
	
	while(true) {
		mprintf(YELLOW, "\n>>> IR Observation %d >>>\n\n", g_IRObserv.GetSize() + 1);
		
		CIRObservation *obs;
		try { obs = new CIRObservation(); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CIRObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_IRObserv.Add(obs);
	
		mprintf("\n");	
		mprintf(YELLOW, "<<< End of IR Observation %d <<<\n\n", g_IRObserv.GetSize());
		
		if(!AskYesNo("    Add another observation (y/n)? [no] ", false))
			break;
		mprintf("\n");
	}
	
	if(AskYesNo("\n    Compute IR spectrum of whole system (y/n) [no] ", false)) {
		mprintf(YELLOW, "\n>>> Global IR Observation >>>\n\n");
		
		CIRObservation *obs;
		try { obs = new CIRObservation(true); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CIRObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_IRObserv.Add(obs);
		
		mprintf(YELLOW, "<<< End of Global IR Observation <<<\n\n");
	}
	
	return true;
}

bool initializeIR() {
	int i;
	for(i = 0; i < g_IRObserv.GetSize(); i++) {
		((CIRObservation *)g_IRObserv[i])->initialize();
	}
	return true;
}

void processIR(CTimeStep* ts) {
	UNUSED(ts);
	int i;
	for(i = 0; i < g_IRObserv.GetSize(); i++) {
		((CIRObservation *)g_IRObserv[i])->process();
	}
}

void finalizeIR() {
	int i;
	for(i = 0; i < g_IRObserv.GetSize(); i++) {
		mprintf(YELLOW, "\n>>> IR Observation %d >>>\n\n", i+1);
		((CIRObservation *)g_IRObserv[i])->finalize();
		delete (CIRObservation *)g_IRObserv[i];
		mprintf(YELLOW, "\n<<< End of IR Observation %d <<<\n\n", i+1);
	}
}

static CxObArray g_DipoleRestartObserv;

CDipoleRestartObservation::CDipoleRestartObservation() 
{
	int i, j, include1, include2;
	FILE* temp;
	CSingleMolecule* sm;
	CxString filename, filename2;
	const char *p, *q;
	char buf[256];
	bool repeat;

	if ( AskYesNo("    Write dipole vectors into human readable files (yes/no)? [no]",false) )
	{
		(void)!system("mkdir dipoles");

		if ( g_oaSingleMolecules.GetSize() > 1000 )
			mprintf(RED,"    Number of Molecules does not allow to print dipole vectors of all molecules.");
		else if ( AskYesNo("    for all molecules (yes/no)? [yes]",true) )
			for ( i=0; i < g_oaSingleMolecules.GetSize(); i++ )
				m_iaMols.Add(i);
		if ( m_iaMols.GetSize() == 0 )
		{
			for ( i=0; i < g_oaMolecules.GetSize(); i++ )
			{
				if ( AskYesNo("    Add representatives of Molecule %d %s (yes/no)? [no]",false,i+1,((CMolecule*)g_oaMolecules[i])->m_sName) )
				{
					repeat = true;
					while ( repeat )
					{
						repeat = false;
						mprintf("    Which representatives of type %s should be added (e.g. \"1,3-5\")? [all] ",((CMolecule*)g_oaMolecules[i])->m_sName);
						inpprintf("! Which representatives of type %s should be added (e.g. \"1,3-5\")? [all] \n",((CMolecule*)g_oaMolecules[i])->m_sName);
						myget(&filename);
						if ( filename.GetLength() == 0 )
							for ( j=0; j < ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex.GetSize(); j++ )
								m_iaMols.Add( ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex[j] );
						else
						{
							p = filename;
							while ( *p != 0 )
							{
								while ( *p == ' ' )
									p++;
								if ( !isdigit(*p) )
								{
									eprintf("Error: Invalid input!\n");
									repeat = true;
									break;    
								}
								q = p;
								while ( isdigit(*q) )
									q++;
								memcpy(buf,p,q-p);
								buf[q-p] = 0;
								include1 = atoi(buf)-1;
								if ( include1 > ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex.GetSize() )
								{
									eprintf("Error: Requested %d of only %d molecules of type %s\n",include1+1,((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex.GetSize()+1,((CMolecule*)g_oaMolecules[i])->m_sName);
									repeat = true;
									break;
								}
								m_iaMols.Add( ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex[include1] );
								if ( *q == '-' )
								{
									q++;
									if ( !isdigit(*p) )
									{
										eprintf("Error: Invalid input!\n");
										repeat = true;
										break;    
									}
									p = q;
									while ( isdigit(*p) )
										p++;
									memcpy(buf,q,p-q);
									buf[p-q] = 0;
									include2 = atoi(buf)-1;
									if ( include1 > ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex.GetSize() )
									{
										eprintf("Error: Requested %d of only %d molecules of type %s\n",include2+1,((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex.GetSize()+1,((CMolecule*)g_oaMolecules[i])->m_sName);
										repeat = true;
										break;
									}
									for ( include1++; include1 <= include2; include1++ )
										m_iaMols.Add( ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex[include1] );
									q = p;
								}
								if ( *q == ',' )
									q++;
								p = q;
							}
						}
					}
				}
			}
		}
		if ( m_iaMols.GetSize() == 0 )
		{
			mprintf(RED,"\n    No molecule selected!\n    Proceeding without writing dipoles to file ...");
			(void)!system("rm -r dipoles");
		}

		for ( i=0; i < m_iaMols.GetSize(); i++ )
		{
			sm = (CSingleMolecule *)g_oaSingleMolecules[m_iaMols[i]];

			try { temp = new FILE(); } catch(...) { temp = NULL; }
			if(temp == NULL) NewException((double)sizeof(FILE), __FILE__, __LINE__, __PRETTY_FUNCTION__);

			#ifdef TARGET_WINDOWS
				filename.sprintf("dipoles\\");
			#else
				filename.sprintf("dipoles/");
			#endif
			filename2.sprintf("%s_%d",((CMolecule*)g_oaMolecules[sm->m_iMolType])->m_sName,sm->m_iMolSMIndex+1);
			filename += filename2;
			temp = OpenFileWrite(filename,true);
			mfprintf(temp," # Dipole Vector: <x>  <y>  <z>  <norm> \n");
			m_oaFiles.Add(temp);
		}
	}
}

CDipoleRestartObservation::~CDipoleRestartObservation() 
{
	int i;

	for ( i=0; i < m_oaFiles.GetSize(); i++ )
		if ( (FILE*)m_oaFiles[i] != NULL )
			fclose((FILE*)m_oaFiles[i]);
}

static FILE *g_dipoleRestartFile = NULL;

bool gatherDipoleRestart() {
	CDipoleRestartObservation *obs;
	try { obs = new CDipoleRestartObservation(); } catch(...) { obs = NULL; }
	if(obs == NULL) NewException((double)sizeof(CDipoleRestartObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_DipoleRestartObserv.Add(obs);

	g_bDipole = true;
	ParseDipole();

	mprintf(YELLOW, "\n>>> Dipole Restart File >>>\n\n");
	
	CxString filename;
	AskString("    Enter name of dipole restart file to write: [dipole.restart] ", &filename, "dipole.restart");
	
	g_dipoleRestartFile = OpenFileWrite((const char *)filename, false);
	int numMolecules = g_oaSingleMolecules.GetSize();
	(void)!fwrite(&numMolecules, sizeof(int), 1, g_dipoleRestartFile);
	
	mprintf(YELLOW, "\n<<< End of Dipole Restart File <<<\n\n");
	return true;
}

bool initializeDipoleRestart() {
	return true;
}

void processDipoleRestart(CTimeStep *ts) {
	UNUSED(ts);
	int i, n;

	for (i = 0; i < g_oaSingleMolecules.GetSize(); i++) {
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
		int j;
		for (j = 0; j < 3; j++) {
			double val = sm->m_vDipole[j];
			(void)!fwrite(&val, sizeof(double), 1, g_dipoleRestartFile);
		}
		for ( j=0; j < g_DipoleRestartObserv.GetSize(); j++ )
		{
			n = ((CDipoleRestartObservation*)g_DipoleRestartObserv[j])->m_iaMols.GetPosition(i);
			if ( n != -1 )
				mfprintf((FILE*)((CDipoleRestartObservation*)g_DipoleRestartObserv[j])->m_oaFiles[n],"%15.10f %15.10f %15.10f %16.10f\n",sm->m_vDipole[0],sm->m_vDipole[1],sm->m_vDipole[2],sm->m_vDipole.GetLength());
		}
	}
	fflush(g_dipoleRestartFile);
}

void finalizeDipoleRestart() {
	int i, j;

	if (g_dipoleRestartFile != NULL)
		fclose(g_dipoleRestartFile);
	for ( j=0; j < g_DipoleRestartObserv.GetSize(); j++ )
		for ( i=0; i < ((CDipoleRestartObservation*)g_DipoleRestartObserv[j])->m_oaFiles.GetSize(); i++ )
			if ( (FILE*)((CDipoleRestartObservation*)g_DipoleRestartObserv[j])->m_oaFiles[i] != NULL )
				fclose((FILE*)((CDipoleRestartObservation*)g_DipoleRestartObserv[j])->m_oaFiles[i]);
}

static CxObArray g_VCDObserv;

CVCDObservation::CVCDObservation(bool global) {
	int i;
	if(global) {
		m_iShowMol = -1;
		m_iShowMolCount = g_oaSingleMolecules.GetSize();
		_name = new char[7];
		sprintf(_name, "global");
	} else {
		char buf[BUF_SIZE];
		char buf2[BUF_SIZE];
		size_t remaining = BUF_SIZE;
		if(g_oaMolecules.GetSize() > 1) {
#ifdef TARGET_LINUX
			remaining -= snprintf(buf, remaining, "    Which molecule should be observed (");
#else
			remaining -= sprintf(buf, "    Which molecule should be observed (");
#endif
			for(i = 0; i < g_oaMolecules.GetSize(); i++) {
				if(remaining < 1)
					break;
#ifdef TARGET_LINUX
				size_t length = snprintf(buf2, remaining, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#else
				size_t length = sprintf(buf2, "%s=%d", ((CMolecule *)g_oaMolecules[i])->m_sName, i+1);
#endif
				strncat(buf, buf2, remaining - 1);
				remaining -= length;
				if(i < g_oaMolecules.GetSize() - 1) {
#ifdef TARGET_LINUX
					length = snprintf(buf2, remaining, ", ");
#else
					length = sprintf(buf2, ", ");
#endif
					strncat(buf, buf2, remaining - 1);
					remaining -= length;
				}
			}
			strncat(buf, ")? ", remaining - 1);
			m_iShowMol = AskRangeInteger_ND("%s", 1, g_oaMolecules.GetSize(),(const char*)buf) - 1;
		} else {
			m_iShowMol = 0;
		}
		m_iShowMolCount = ((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex.GetSize();
		_name = new char[strlen(((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName) + 1];
		strcpy(_name, ((CMolecule *)g_oaMolecules[m_iShowMol])->m_sName);
		mprintf("\n");
	}
	
	if(g_iTrajSteps != -1) {
		_correlationDepth = (int)(0.75 * g_iTrajSteps);
		if(_correlationDepth > 4096)
			_correlationDepth = 4096;
		if(g_fTimestepLength > 1.0)
			_correlationDepth = 2048;
		if(g_fTimestepLength > 2.0)
			_correlationDepth = 1024;
		_correlationDepth = AskUnsignedInteger("    Enter the correlation depth of the ACF (in trajectory frames): [%d] ", _correlationDepth, _correlationDepth);
	} else {
		_correlationDepth = AskUnsignedInteger("    Enter the correlation depth of the ACF (in trajectory frames): [256] ", 256);
	}
	int size = CalcFFTSize(_correlationDepth, false);
	if(_correlationDepth != size) {
		mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using this instead of %d as depth.\n", size, _correlationDepth);
		_correlationDepth = size;
	}
	
	if(g_bAdvanced2) {
		_windowFunction = AskRangeInteger("    Window function: cos^2(a*t) (1), exp(-t/a) (2), exp(-(t/a)^2) (3) [1] ", 1, 3, 1);
		if(_windowFunction == 1) {
			mprintf("    The parameter \"a\" is chosen according to the correlation depth.\n");
			_windowFunctionParameter = 0;
		} else if(_windowFunction == 2) {
			_windowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in trajectory frames): [%d] ", _correlationDepth / 4, _correlationDepth / 4);
		} else if(_windowFunction == 3) {
			_windowFunctionParameter = AskUnsignedInteger("    Parameter \"a\" (in trajectory frames): [%d] ", _correlationDepth / 2, _correlationDepth / 2);
		} else {
			eprintf("This is impossible.\n");
			abort();
		}
	} else {
		_windowFunction = 1;
		_windowFunctionParameter = 0;
		mprintf("    Using cos^2 window function.\n");
	}
	
	if(g_bAdvanced2) {
		_zeroPadding = AskUnsignedInteger("    Zero Padding: How many zeros to insert? [%d] ", _correlationDepth * 3, _correlationDepth * 3);
		size = CalcFFTSize(_correlationDepth + _zeroPadding, false);
		if(_correlationDepth + _zeroPadding != size) {
			mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-_correlationDepth);
			_zeroPadding = size-_correlationDepth;
		}
	} else {
		_zeroPadding = _correlationDepth * 3;
		size = CalcFFTSize(_correlationDepth + _zeroPadding, false);
		if(_correlationDepth + _zeroPadding != size) {
			mprintf(WHITE, "    The next \"fast\" size for FFT is %d. Using %d zeros for zero padding.\n", size, size-_correlationDepth);
			_zeroPadding = size-_correlationDepth;
		}
		mprintf("    Inserting %d zeros for zero padding.\n",_zeroPadding);
	}
	
	double possibleRange = 33356.41 / g_fTimestepLength / 2.0;
	_specResolution = possibleRange / (_correlationDepth + _zeroPadding);
	mprintf("    This results in a spectral resolution of %.2f cm^-1.\n", _specResolution);
	mprintf("\n    A time step length of %.2f fs allows a spectral range up to %.2f cm^-1.\n", g_fTimestepLength, possibleRange);
	double specLimit = AskRangeFloat("\n    Calculate spectrum up to which wave number (cm^-1)? [%.2f cm^-1] ", 0, possibleRange, (possibleRange < 5000.0) ? possibleRange : 5000.0, (possibleRange < 5000.0) ? possibleRange : 5000.0);
	_specSize = (int)(specLimit / _specResolution);
	mprintf("\n");
	
	if(g_bAdvanced2) {
		_smoothWidth = AskUnsignedInteger("    Window width of median filter to smooth magnetic dipole moments: [0] ", 0);
		mprintf("\n");
	} else {
		_smoothWidth = 0;
	}
	
	if(g_bAdvanced2) {
		_saveMoments = AskYesNo("    Save electric and magnetic moments (y/n)? [no] ", false);
	} else {
		_saveMoments = false;
	}
	
	if(g_bAdvanced2) {
		_saveACF = AskYesNo("    Save cross-correlation function (y/n)? [no] ", false);
		mprintf("\n");
	} else {
		_saveACF = false;
	}
}

CVCDObservation::~CVCDObservation() {
	delete[] _name;
}

void CVCDObservation::initialize() {
	int n;
	if(g_iTrajSteps != -1)
		n = (int)(1.1 * g_iTrajSteps / g_iStride);
	else
		n = 10000;
	
	mprintf("    Moment cache: Trying to allocate %s of memory...\n", FormatBytes((double)m_iShowMolCount * 2.0 * n * sizeof(CxDVector3)));
	int i;
	for(i = 0; i < m_iShowMolCount; i++) {
		CxDVec3Array *a, *b;
		try { a = new CxDVec3Array(); } catch(...) { a = NULL; }
		if(a == NULL) NewException((double)sizeof(CxDVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		try { b = new CxDVec3Array(); } catch(...) { b = NULL; }
		if(b == NULL) NewException((double)sizeof(CxDVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		a->SetMaxSize(n);
		a->SetGrow(n / 10);
		b->SetMaxSize(n);
		b->SetGrow(n / 10);
		_electricDipoleCache.Add(a);
		_magneticDipoleCache.Add(b);
	}
}

void CVCDObservation::process() {
	int i;
	if (m_iShowMol == -1) {
		for (i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
			((CxDVec3Array *)_electricDipoleCache[i])->Add(sm->m_vDipole);
			((CxDVec3Array *)_magneticDipoleCache[i])->Add(sm->m_vMagneticDipole);
		}
	} else {
		for (i = 0; i < m_iShowMolCount; i++) {
			CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[((CMolecule *)g_oaMolecules[m_iShowMol])->m_laSingleMolIndex[i]];
			((CxDVec3Array *)_electricDipoleCache[i])->Add(sm->m_vDipole);
			((CxDVec3Array *)_magneticDipoleCache[i])->Add(sm->m_vMagneticDipole);
		}
	}
}

static int compare_double(const void *f1, const void *f2) {
	if(*(double *)f1 < *(double *)f2)
		return -1;
	if(*(double *)f1 > *(double *)f2)
		return 1;
	return 0;
}

void CVCDObservation::finalize() {
	int i, j, k, l;
	if(_smoothWidth > 0) {
		mprintf("    Smoothing magnetic dipole moments...\n");
		CxDVec3Array *tempMoments;
		try { tempMoments = new CxDVec3Array(); } catch(...) { tempMoments = NULL; }
		if(tempMoments == NULL) NewException((double)sizeof(CxDVec3Array), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		double *sortArray;
		try { sortArray = new double[2*_smoothWidth+1]; } catch(...) { sortArray = NULL; }
		if(sortArray == NULL) NewException((double)_smoothWidth * sizeof(double), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		for(i = 0; i < m_iShowMolCount; i++) {
			tempMoments->SetSize(((CxDVec3Array *)_magneticDipoleCache[i])->GetSize());
			for(j = 0; j < tempMoments->GetSize(); j++) {
				for(k = 0; k < 3; k++) {
					int index = 0;
					for(l = j - _smoothWidth; l <= j + _smoothWidth; l++) {
						if(l < 0)
							continue;
						if(l >= tempMoments->GetSize())
							break;
						sortArray[index] = ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(l)[k];
						index++;
					}
					qsort(sortArray, index, sizeof(double), compare_double);
					if(index % 2 == 0) {
						tempMoments->GetAt(j)[k] = (sortArray[index / 2] + sortArray[index / 2 - 1]) / 2.0;
					} else {
						tempMoments->GetAt(j)[k] = sortArray[index / 2];
					}
				}
			}
			
// 			for(j = 0; j < 9; j++) {
// 				tempMoments->GetAt(j) = ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j);
// 			}
// 			for(j = 9; j < tempMoments->GetSize() - 9; j++) {
// 				int k;
// 				for(k = 0; k < 3; k++) {
// 					double array[19];
// 					int l;
// 					for(l = j - 9; l <= j + 9; l++) {
// 						array[l - j + 9] = ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(l)[k];
// 					}
// 					qsort(array, 19, sizeof(double), compare_double);
// 					tempMoments->GetAt(j)[k] = array[9];
// 				}
// 			}
// 			for(j = tempMoments->GetSize() - 9; j < tempMoments->GetSize(); j++) {
// 				tempMoments->GetAt(j) = ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j);
// 			}
// 		for(j = 0; j < tempMoments->GetSize(); j++) {
// 			tempMoments->GetAt(j) = CxDVector3(0.0f, 0.0f, 0.0f);
// 			int k;
// 			int count = 0;
// 			for(k = j - 10; k <= j + 10; k++) {
// 				if(k < 0)
// 					continue;
// 				if(k >= tempMoments->GetSize())
// 					break;
// 				tempMoments->GetAt(j) += ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(k);
// 				count++;
// 			}
// 			tempMoments->GetAt(j) /= (double)count;
// 		}
// 		for(j = 0; j < 12; j++) {
// 			tempMoments->GetAt(j) = ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j);
// 		}
// 		for(j = 12; j < tempMoments->GetSize() - 12; j++) {
// // 			tempMoments->GetAt(j) = (-3.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-2) + 12.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-1) + 17.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j) + 12.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+1) - 3.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+2)) / 35.0f;
// // 			tempMoments->GetAt(j) = (-21.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-4) + 14.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-3) + 39.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-2) + 54.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-1) + 59.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j) + 54.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+1) + 39.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+2) + 14.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+3) - 21.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+4)) / 231.0f;
// tempMoments->GetAt(j) = (-253.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-12) - 138.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-11) - 33.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-10) + 62.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-9) + 147.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-8) + 222.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-7) + 287.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-6) + 322.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-5) + 387.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-4) + 422.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-3) + 447.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-2) + 462.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j-1) + 467.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j) + 462.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+1) + 447.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+2) +
// 422.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+3) + 387.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+4) + 322.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+5) + 287.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+6) + 222.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+7) + 147.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+8) + 62.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+9) - 33.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+10) - 138.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+11) - 253.0f * ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+12)) / 5175.0f;
// 		}
// 		for(j = tempMoments->GetSize() - 12; j < tempMoments->GetSize(); j++) {
// 			tempMoments->GetAt(j) = ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j);
// 		}
			((CxDVec3Array *)_magneticDipoleCache[i])->CopyFrom(tempMoments);
		}
		delete tempMoments;
		delete[] sortArray;
	}
	
	if(_saveMoments) {
		mprintf("    Saving electric and magnetic dipole moments...\n");
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "vcd_moments_%s.csv", _name);
#else
		sprintf(filename, "vcd_moments_%s.csv", _name);
#endif
		FILE *momentsFile = OpenFileWrite(filename, false);
		for(i = 0; i < ((CxDVec3Array *)_magneticDipoleCache[0])->GetSize(); i++) {
			fprintf(momentsFile, "%.2f;", g_fTimestepLength * (i + 1));
			for(j = 0; j < m_iShowMolCount; j++) {
				fprintf(momentsFile, " %.10G; %.10G; %.10G; %.10G; %.10G; %.10G;", ((CxDVec3Array *)_electricDipoleCache[j])->GetAt(i)[0], ((CxDVec3Array *)_electricDipoleCache[j])->GetAt(i)[1], ((CxDVec3Array *)_electricDipoleCache[j])->GetAt(i)[2], ((CxDVec3Array *)_magneticDipoleCache[j])->GetAt(i)[0], ((CxDVec3Array *)_magneticDipoleCache[j])->GetAt(i)[1], ((CxDVec3Array *)_magneticDipoleCache[j])->GetAt(i)[2]);
			}
			fprintf(momentsFile, "\n");
		}
		fclose(momentsFile);
	}
	
	mprintf("    Deriving electric dipole moments...\n");
	int n = ((CxDVec3Array *)_electricDipoleCache[0])->GetSize() - 2;
	double step = (double)m_iShowMolCount / 20.0;
	mprintf(WHITE, "     [");
	for(i = 0; i < m_iShowMolCount; i++) {
		if(fmod(i, step) < 1.0)
			mprintf(WHITE, "#");
		for(j = 0; j < n; j++) {
			((CxDVec3Array *)_electricDipoleCache[i])->GetAt(j) = 0.5 * (((CxDVec3Array *)_electricDipoleCache[i])->GetAt(j+2) - ((CxDVec3Array *)_electricDipoleCache[i])->GetAt(j)) / g_fTimestepLength;
		}
	}
	mprintf(WHITE, "]\n");
	mprintf("    Deriving magnetic dipole moments...\n");
	mprintf(WHITE, "     [");
	for(i = 0; i < m_iShowMolCount; i++) {
		if(fmod(i, step) < 1.0)
			mprintf(WHITE, "#");
		for(j = 0; j < n; j++) {
			((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j) = 0.5 * (((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j+2) - ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(j)) / g_fTimestepLength;
		}
	}
	mprintf(WHITE, "]\n");
	
	if(_saveMoments) {
		mprintf("    Saving electric and magnetic dipole moment derivatives...\n");
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "vcd_moments_deriv_%s.csv", _name);
#else
		sprintf(filename, "vcd_moments_deriv_%s.csv", _name);
#endif
		FILE *momentsFile = OpenFileWrite(filename, false);
		for(i = 0; i < n; i++) {
			fprintf(momentsFile, "%.2f;", g_fTimestepLength * (i + 2));
			for(j = 0; j < m_iShowMolCount; j++) {
				fprintf(momentsFile, " %.10G; %.10G; %.10G; %.10G; %.10G; %.10G;", ((CxDVec3Array *)_electricDipoleCache[j])->GetAt(i)[0], ((CxDVec3Array *)_electricDipoleCache[j])->GetAt(i)[1], ((CxDVec3Array *)_electricDipoleCache[j])->GetAt(i)[2], ((CxDVec3Array *)_magneticDipoleCache[j])->GetAt(i)[0], ((CxDVec3Array *)_magneticDipoleCache[j])->GetAt(i)[1], ((CxDVec3Array *)_magneticDipoleCache[j])->GetAt(i)[2]);
			}
			fprintf(momentsFile, "\n");
		}
		fclose(momentsFile);
	}
	
	mprintf("    Calculating cross-correlation...\n");
	if ((n < 2) || (_correlationDepth < 2)) {
		mprintf("Less than 2 steps, skipping correlation.\n");
		for(i = 0; i < _electricDipoleCache.GetSize(); i++) {
			delete (CxDVec3Array *)_electricDipoleCache[i];
			delete (CxDVec3Array *)_magneticDipoleCache[i];
		}
		return;
	}
	CxDoubleArray *ccf;
	try { ccf = new CxDoubleArray(); } catch(...) { ccf = NULL; }
	if(ccf == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	ccf->SetSize(_correlationDepth);
	for(i = 0; i < _correlationDepth; i++) {
		ccf->GetAt(i) = 0.0;
	}
	CCrossCorrelation *cc;
	try { cc = new CCrossCorrelation(); } catch(...) { cc = NULL; }
	if(cc == NULL) NewException((double)sizeof(CCrossCorrelation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	cc->Init(n, _correlationDepth, g_bACFFFT);
	
	CxDoubleArray *temp;
	try { temp = new CxDoubleArray(); } catch(...) { temp = NULL; }
	if(temp == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp->SetSize(n);
	CxDoubleArray *temp2;
	try { temp2 = new CxDoubleArray(); } catch(...) { temp2 = NULL; }
	if(temp2 == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp2->SetSize(n);
	CxDoubleArray *temp3;
	try { temp3 = new CxDoubleArray(); } catch(...) { temp3 = NULL; }
	if(temp3 == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	temp3->SetSize(n);
	
	mprintf(WHITE, "     [");
	for(i = 0; i < m_iShowMolCount; i++) {
		if(fmod(i, step) < 1.0)
			mprintf(WHITE, "#");
		for(j = 0; j < 3; j++) {
			for(k = 0; k < n; k++) {
				temp->GetAt(k) = ((CxDVec3Array *)_electricDipoleCache[i])->GetAt(k)[j];
				temp2->GetAt(k) = ((CxDVec3Array *)_magneticDipoleCache[i])->GetAt(k)[j];
			}
			cc->CrossCorrelate(temp, temp2, temp3);
			for(k = 0; k < _correlationDepth; k++) {
				ccf->GetAt(k) += temp3->GetAt(k);
			}
			cc->CrossCorrelate(temp2, temp, temp3);
			for(k = 0; k < _correlationDepth; k++) {
				ccf->GetAt(k) -= temp3->GetAt(k);
			}
		}
	}
	mprintf(WHITE, "]\n");
// The 3.0f is included in the final normalization
// 	for(i = 0; i < _correlationDepth; i++) {
// 		ccf->GetAt(i) /= 3.0f * g_iSteps * m_iShowMolCount;
// 	}
	if(m_iShowMol != -1) {
		for(i = 0; i < _correlationDepth; i++) {
			ccf->GetAt(i) /= (double)m_iShowMolCount;
		}
	}
	
	delete cc;
	delete temp2;
	delete temp3;
	
	if(_saveACF) {
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "vcd_ccf_%s.csv", _name);
#else
		sprintf(filename, "vcd_ccf_%s.csv", _name);
#endif
		mprintf("    Saving cross-correlation function as %s...\n", filename);
		FILE *ccfFile = OpenFileWrite(filename, false);
		for(i = 0; i < _correlationDepth; i++) {
			fprintf(ccfFile, "%.2f; %.10G\n", i * g_fTimestepLength, ccf->GetAt(i));
		}
		fclose(ccfFile);
	}
	
	temp->CopyFrom(ccf);
	
	if(_windowFunction == 1) {
		for(i = 0; i < temp->GetSize(); i++) {
			temp->GetAt(i) *= pow2(cos((double)i / (temp->GetSize() - 1) / 2.0 * Pi));
		}
	} else if(_windowFunction == 2) {
		for(i = 0; i < temp->GetSize(); i++) {
			temp->GetAt(i) *= exp(-(double)i / _windowFunctionParameter);
		}
	} else if(_windowFunction == 3) {
		for(i = 0; i < temp->GetSize(); i++) {
			temp->GetAt(i) *= exp(-(double)i * i / _windowFunctionParameter / _windowFunctionParameter);
		}
	} else if(_windowFunction != 0) {
		eprintf("Unknown window function.\n");
		abort();
	}
	
	if(_saveACF) {
		char filename[BUF_SIZE];
#ifdef TARGET_LINUX
		snprintf(filename, BUF_SIZE, "vcd_ccf_windowed_%s.csv", _name);
#else
		sprintf(filename, "vcd_ccf_windowed_%s.csv", _name);
#endif
		mprintf("    Saving windowed cross-correlation function as %s...\n", filename);
		FILE *ccfFile = OpenFileWrite(filename, false);
		for(i = 0; i < _correlationDepth; i++) {
			fprintf(ccfFile, "%.2f; %.10G\n", i * g_fTimestepLength, temp->GetAt(i));
		}
		fclose(ccfFile);
	}
	
	if(_zeroPadding > 0) {
		for(i = 0; i < _zeroPadding; i++) {
			temp->Add(0.0);
		}
	}
	
	int oldSize = temp->GetSize();
	temp->SetSize(2 * oldSize);
	for(i = 1; i < oldSize; i++) {
		temp->GetAt(oldSize + i) = -temp->GetAt(oldSize - i);
	}
	temp->GetAt(oldSize) = 0.0;
	
	mprintf("    Fourier transforming cross-correlation function...\n");
	CFFT *fft;
	try { fft = new CFFT(); } catch(...) { fft = NULL; }
	if(fft == NULL) NewException((double)sizeof(CFFT), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	fft->PrepareFFT_C2C(temp->GetSize());
	for(i = 0; i < temp->GetSize(); i++) {
		fft->m_pInput[2*i] = temp->GetAt(i);
		fft->m_pInput[2*i+1] = 0.0;
	}
	fft->DoFFT();

	CxDoubleArray *spectrum;
	try { spectrum = new CxDoubleArray(); } catch(...) { spectrum = NULL; }
	if(spectrum == NULL) NewException((double)sizeof(CxDoubleArray), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	spectrum->SetSize(_specSize);
	for(i = 0; i < _specSize; i++) {
		spectrum->GetAt(i) = 28.260058 * fft->m_pOutput[2*i+1] * g_fTimestepLength; // Output in K*cm*km/mol
	}
	
	if(_finiteDifferenceCorrection) {
		double f = _specResolution * g_fTimestepLength * 1.883652e-4;
		for(i = 1; i < _specSize; i++) {
			spectrum->GetAt(i) *= pow2(f * i / sin(f * i)); // Divide by sinc function to correct finite difference derivation
		}
	}
	
	char filename[BUF_SIZE];
#ifdef TARGET_LINUX
	snprintf(filename, BUF_SIZE, "vcd_spectrum_%s.csv", _name);
#else
	sprintf(filename, "vcd_spectrum_%s.csv", _name);
#endif
	mprintf("    Saving spectrum as %s...\n", filename);
	FILE *specFile = OpenFileWrite(filename, false);
	fprintf(specFile, "#Wavenumber (cm^-1); Spectrum (K*cm*km*mol^-1); Integral (K*km*mol^-1)\n");
	double integral = 0.0;
	for(i = 0; i < _specSize; i++) {
		integral += spectrum->GetAt(i) * _specResolution;
		fprintf(specFile, "%.2f; %.8G; %.14G\n", _specResolution * i, spectrum->GetAt(i), integral);
	}
	fclose(specFile);
	
	delete ccf;
	delete temp;
	delete fft;
	delete spectrum;
	
	for(i = 0; i < _electricDipoleCache.GetSize(); i++) {
		delete (CxDVec3Array *)_electricDipoleCache[i];
		delete (CxDVec3Array *)_magneticDipoleCache[i];
	}
}

bool gatherVCD() {
	parseMagneticDipole();
	ParseDipole();
	
	while(true) {
		mprintf(YELLOW, "\n>>> VCD Observation %d >>>\n\n", g_VCDObserv.GetSize() + 1);
		
		CVCDObservation *obs;
		try { obs = new CVCDObservation(); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CVCDObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_VCDObserv.Add(obs);
		
		mprintf(YELLOW, "<<< End of VCD Observation %d <<<\n\n", g_VCDObserv.GetSize());
		
		if(!AskYesNo("    Add another observation (y/n)? [no] ", false))
			break;
		mprintf("\n");
	}
	
	if(AskYesNo("    Compute VCD spectrum of whole system (y/n) [no] ", false)) {
		mprintf(YELLOW, "\n>>> Global VCD Observation >>>\n\n");
		
		CVCDObservation *obs;
		try { obs = new CVCDObservation(true); } catch(...) { obs = NULL; }
		if(obs == NULL) NewException((double)sizeof(CVCDObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
		g_VCDObserv.Add(obs);
		
		mprintf(YELLOW, "<<< End of Global VCD Observation <<<\n\n");
	}
	
	return true;
}

bool initializeVCD() {
	int i;
	for(i = 0; i < g_VCDObserv.GetSize(); i++) {
		((CVCDObservation *)g_VCDObserv[i])->initialize();
	}
	return true;
}

void processVCD(CTimeStep *ts) {
	UNUSED(ts);
// 	ts->CalcMagneticDipoles();
	int i;
	for(i = 0; i < g_VCDObserv.GetSize(); i++) {
		((CVCDObservation *)g_VCDObserv[i])->process();
	}
}

void finalizeVCD() {
	int i;
	for(i = 0; i < g_VCDObserv.GetSize(); i++) {
		((CVCDObservation *)g_VCDObserv[i])->finalize();
	}
}

static CxObArray g_MagneticRestartObserv;

CMagneticRestartObservation::CMagneticRestartObservation() 
{
	int i, j, include1, include2;
	FILE* temp;
	CSingleMolecule* sm;
	CxString filename, filename2;
	const char *p, *q;
	char buf[256];
	bool repeat;

	if ( AskYesNo("    Write magnetic dipole vectors into human readable files (yes/no)? [no]",false) )
	{
		(void)!system("mkdir magnetics");

		if ( g_oaSingleMolecules.GetSize() > 1000 )
			mprintf(RED,"    Number of Molecules does not allow to print magnetic dipole vectors of all molecules.");
		else if ( AskYesNo("    for all molecules (yes/no)? [yes]",true) )
			for ( i=0; i < g_oaSingleMolecules.GetSize(); i++ )
				m_iaMols.Add(i);
		if ( m_iaMols.GetSize() == 0 )
		{
			for ( i=0; i < g_oaMolecules.GetSize(); i++ )
			{
				if ( AskYesNo("    Add representatives of Molecule %d %s (yes/no)? [no]",false,i+1,((CMolecule*)g_oaMolecules[i])->m_sName) )
				{
					repeat = true;
					while ( repeat )
					{
						repeat = false;
						mprintf("    Which representatives of type %s should be added (e.g. \"1,3-5\")? [all] ",((CMolecule*)g_oaMolecules[i])->m_sName);
						inpprintf("! Which representatives of type %s should be added (e.g. \"1,3-5\")? [all] \n",((CMolecule*)g_oaMolecules[i])->m_sName);
						myget(&filename);
						if ( filename.GetLength() == 0 )
							for ( j=0; j < ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex.GetSize(); j++ )
								m_iaMols.Add( ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex[j] );
						else
						{
							p = filename;
							while ( *p != 0 )
							{
								while ( *p == ' ' )
									p++;
								if ( !isdigit(*p) )
								{
									eprintf("Error: Invalid input!\n");
									repeat = true;
									break;    
								}
								q = p;
								while ( isdigit(*q) )
									q++;
								memcpy(buf,p,q-p);
								buf[q-p] = 0;
								include1 = atoi(buf)-1;
								if ( include1 > ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex.GetSize() )
								{
									eprintf("Error: Requested %d of only %d molecules of type %s\n",include1+1,((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex.GetSize()+1,((CMolecule*)g_oaMolecules[i])->m_sName);
									repeat = true;
									break;
								}
								m_iaMols.Add( ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex[include1] );
								if ( *q == '-' )
								{
									q++;
									if ( !isdigit(*p) )
									{
										eprintf("Error: Invalid input!\n");
										repeat = true;
										break;    
									}
									p = q;
									while ( isdigit(*p) )
										p++;
									memcpy(buf,q,p-q);
									buf[p-q] = 0;
									include2 = atoi(buf)-1;
									if ( include1 > ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex.GetSize() )
									{
										eprintf("Error: Requested %d of only %d molecules of type %s\n",include2+1,((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex.GetSize()+1,((CMolecule*)g_oaMolecules[i])->m_sName);
										repeat = true;
										break;
									}
									for ( include1++; include1 <= include2; include1++ )
										m_iaMols.Add( ((CMolecule*)g_oaMolecules[i])->m_laSingleMolIndex[include1] );
									q = p;
								}
								if ( *q == ',' )
									q++;
								p = q;
							}
						}
					}
				}
			}
		}
		if ( m_iaMols.GetSize() == 0 )
		{
			mprintf(RED,"\n    No molecule selected!\n    Proceeding without writing magnetic dipoles to file ...");
			(void)!system("rm -r magnetics");
		}

		for ( i=0; i < m_iaMols.GetSize(); i++ )
		{
			sm = (CSingleMolecule *)g_oaSingleMolecules[m_iaMols[i]];

			try { temp = new FILE(); } catch(...) { temp = NULL; }
			if(temp == NULL) NewException((double)sizeof(FILE), __FILE__, __LINE__, __PRETTY_FUNCTION__);

			#ifdef TARGET_WINDOWS
				filename.sprintf("magnetics\\");
			#else
				filename.sprintf("magnetics/");
			#endif
			filename2.sprintf("%s_%d",((CMolecule*)g_oaMolecules[sm->m_iMolType])->m_sName,sm->m_iMolSMIndex+1);
			filename += filename2;
			temp = OpenFileWrite(filename,true);
			mfprintf(temp," # Magnetic Dipole Vector: <x>  <y>  <z>  <norm> \n");
			m_oaFiles.Add(temp);
		}
	}
}

CMagneticRestartObservation::~CMagneticRestartObservation() 
{
	int i;

	for ( i=0; i < m_oaFiles.GetSize(); i++ )
		if ( (FILE*)m_oaFiles[i] != NULL )
			fclose((FILE*)m_oaFiles[i]);
}

static FILE *g_magneticDipoleRestartFile = NULL;

bool gatherMagneticDipoleRestart() {
	CMagneticRestartObservation *obs;
	try { obs = new CMagneticRestartObservation(); } catch(...) { obs = NULL; }
	if(obs == NULL) NewException((double)sizeof(CMagneticRestartObservation), __FILE__, __LINE__, __PRETTY_FUNCTION__);
	g_MagneticRestartObserv.Add(obs);

	parseMagneticDipole();
	ParseDipole();

	mprintf(YELLOW, "\n>>> Magnetic Moment Restart File >>>\n\n");
	
	CxString filename;
	AskString("    Enter name of magnetic moment restart file to write: [magnetic.restart] ", &filename, "magnetic.restart");
	
	g_magneticDipoleRestartFile = OpenFileWrite((const char *)filename, false);
	int numMolecules = g_oaSingleMolecules.GetSize();
	(void)!fwrite(&numMolecules, sizeof(int), 1, g_magneticDipoleRestartFile);
	
	mprintf(YELLOW, "\n<<< End of Magnetic Moment Restart File <<<\n\n");
	return true;
}

bool initializeMagneticDipoleRestart() {
	return true;
}

void processMagneticDipoleRestart(CTimeStep *ts) {
	UNUSED(ts);
	int i, n;

	for (i = 0; i < g_oaSingleMolecules.GetSize(); i++) {
		CSingleMolecule *sm = (CSingleMolecule *)g_oaSingleMolecules[i];
		int j;
		for (j = 0; j < 3; j++) {
			double val = sm->m_vMagneticDipole[j];
			(void)!fwrite(&val, sizeof(double), 1, g_magneticDipoleRestartFile);
		}
		for ( j=0; j < g_MagneticRestartObserv.GetSize(); j++ )
		{
			n = ((CMagneticRestartObservation*)g_MagneticRestartObserv[j])->m_iaMols.GetPosition(i);
			if ( n != -1 )
				mfprintf((FILE*)((CMagneticRestartObservation*)g_MagneticRestartObserv[j])->m_oaFiles[n],"%15.10f %15.10f %15.10f %16.10f\n",sm->m_vMagneticDipole[0],sm->m_vMagneticDipole[1],sm->m_vMagneticDipole[2],sm->m_vMagneticDipole.GetLength());
		}		
	}
	fflush(g_magneticDipoleRestartFile);
}

void finalizeMagneticDipoleRestart() {
	int i, j;

	if (g_magneticDipoleRestartFile != NULL)
		fclose(g_magneticDipoleRestartFile);
	for ( j=0; j < g_MagneticRestartObserv.GetSize(); j++ )
		for ( i=0; i < ((CDipoleRestartObservation*)g_MagneticRestartObserv[j])->m_oaFiles.GetSize(); i++ )
			if ( (FILE*)((CDipoleRestartObservation*)g_MagneticRestartObserv[j])->m_oaFiles[i] != NULL )
				fclose((FILE*)((CDipoleRestartObservation*)g_MagneticRestartObserv[j])->m_oaFiles[i]);
}

static FILE *g_sortWannierTrajFile;
static int g_sortWannierCount;
static CxDVec3Array g_sortWannierHistory;
static CxIntArray g_sortWannierHistoryIndex;
static CxDoubleArray g_sortWannierDistanceMatrix;
static CxIntArray g_sortWannierKMStarCol;
static CxIntArray g_sortWannierKMPrimeRow;
static CxDoubleArray g_sortWannierKMMinCol;
static CxDoubleArray g_sortWannierKMMinRow;

bool gatherSortWannier() {
	g_bKeepOriginalCoords = true;
	g_bWannier = true;
	int watom = -1;
	int i;
	for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
		if(((CAtom *)g_oaAtoms[i])->m_bExclude) {
			watom = i;
			break;
		}
	}
	if(watom == -1)
		eprintf("    Atom label of Wannier centers not found.\n\n");
	else
		mprintf("    Atom type %s is excluded from the system, probably these are the Wannier centers.\n\n", (const char*)((CAtom *)g_oaAtoms[watom])->m_sName);
	
	bool ok = false;
	char buf[BUF_SIZE];
	char buf2[BUF_SIZE];
	size_t remaining = BUF_SIZE;
	while(!ok) {
#ifdef TARGET_LINUX
		remaining -= snprintf(buf, remaining, "    Which atom label do the wannier centers have (");
#else
		remaining -= sprintf(buf, "    Which atom label do the wannier centers have (");
#endif
		for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
			if(remaining < 1)
				break;
#ifdef TARGET_LINUX
			size_t length = snprintf(buf2, remaining, "%s", (const char*)((CAtom *)g_oaAtoms[i])->m_sName);
#else
			size_t length = sprintf(buf2, "%s", (const char*)((CAtom *)g_oaAtoms[i])->m_sName);
#endif
			strncat(buf, buf2, remaining - 1);
			remaining -= length;
			if(i < g_oaAtoms.GetSize() - 2) {
#ifdef TARGET_LINUX
				length = snprintf(buf2, remaining, ", ");
#else
				length = sprintf(buf2, ", ");
#endif
				strncat(buf, buf2, remaining - 1);
				remaining -= length;
			}
		}
#ifdef TARGET_LINUX
		size_t length = snprintf(buf2, remaining, ")? ");
#else
		size_t length = sprintf(buf2, ")? ");
#endif
		strncat(buf, buf2, remaining - 1);
		remaining -= length;
		if(watom != -1) {
#ifdef TARGET_LINUX
			snprintf(buf2, remaining, "[%s] ", (const char*)((CAtom *)g_oaAtoms[watom])->m_sName);
#else
			sprintf(buf2, "[%s] ", (const char*)((CAtom *)g_oaAtoms[watom])->m_sName);
#endif
			strncat(buf, buf2, remaining - 1);
		}
		
		CxString buf3;
		if(watom == -1)
			AskString_ND("%s", &buf3, (const char*)buf);
		else
			AskString("%s", &buf3, ((CAtom *)g_oaAtoms[watom])->m_sName, (const char*)buf);
		
		for(i = 0; i < g_oaAtoms.GetSize() - 1; i++) {
			if(mystricmp((const char *)buf3, ((CAtom *)g_oaAtoms[i])->m_sName) == 0) {
				g_iWannierAtomType = i;
				ok = true;
				break;
			}
		}
		if(!ok) {
			eprintf("    Invalid input.\n");
			inpprintf("! Invalid input.\n");
		}
	}
	
	g_sortWannierCount = 0;
	for(i = 0; i < g_iGesAtomCount; i++) {
		if(g_baAtomIndex[i] == g_iWannierAtomType)
			g_sortWannierCount++;
	}
	mprintf("\n    There are %d Wannier centers in the first step.\n", g_sortWannierCount);
	
	return true;
}

bool initializeSortWannier() {
	char filename[BUF_SIZE];
#ifdef TARGET_LINUX
	char buf[BUF_SIZE];
	strncpy(buf, g_sInputTraj, BUF_SIZE-1);
	buf[BUF_SIZE-1] = 0;
	char *p = strrchr(buf, '/');
	if(p == NULL)
		p = buf;
	else
		p++;
	size_t s = strcspn(p, ".");
	if(s > BUF_SIZE - 8)
		s = BUF_SIZE - 8;
	strncpy(filename, p, s);
	filename[s] = 0;
	strcat(filename, "_out.xyz");
#else
	sprintf(filename, "out.xyz");
#endif
	mprintf("    Saving processed trajectory as %s\n\n", filename);
	g_sortWannierTrajFile = OpenFileWrite(filename, false);
	
	mprintf("    Coordinate history: Trying to allocate %s of memory...\n", FormatBytes((double)(g_sortWannierCount * sizeof(int) + g_sortWannierCount * sizeof(CxDVector3))));
	g_sortWannierHistory.SetSize(g_sortWannierCount);
	g_sortWannierHistoryIndex.SetSize(g_sortWannierCount);
	
	mprintf("    Hungarian algorithm: Trying to allocate %s of memory...\n", FormatBytes((double)(g_sortWannierCount * g_sortWannierCount * sizeof(double) + 2 * g_sortWannierCount * sizeof(double) + 2 * g_sortWannierCount * sizeof(int))));
	g_sortWannierDistanceMatrix.SetSize(g_sortWannierCount * g_sortWannierCount);
	g_sortWannierKMStarCol.SetSize(g_sortWannierCount);
	g_sortWannierKMPrimeRow.SetSize(g_sortWannierCount);
	g_sortWannierKMMinCol.SetSize(g_sortWannierCount);
	g_sortWannierKMMinRow.SetSize(g_sortWannierCount);
	
	return true;
}

void processSortWannier(CTimeStep *ts) {
	static bool first = true;
	if(first) {
		first = false;
		fprintf(g_sortWannierTrajFile, "%d\n", g_iGesAtomCount);
		if(ts->m_pComment != NULL)
			fprintf(g_sortWannierTrajFile, "%s\n", ts->m_pComment);
		else
			fprintf(g_sortWannierTrajFile, "\n");
		
		int i;
		int pos = 0;
		for(i = 0; i < g_iGesAtomCount; i++) {
			if(g_baAtomIndex[i] == g_iWannierAtomType) {
				g_sortWannierHistory[pos] = ts->m_vaCoords[i];
				g_sortWannierHistoryIndex[pos] = i;
				pos++;
				fprintf(g_sortWannierTrajFile, "%4s %14.10f %14.10f %14.10f\n", (const char*)((CAtom *)g_oaAtoms[g_iWannierAtomType])->m_sName, ts->m_vaCoords_Original[i][0] / 100.0f, ts->m_vaCoords_Original[i][1] / 100.0f, ts->m_vaCoords_Original[i][2] / 100.0f);
			} else {
				fprintf(g_sortWannierTrajFile, "%4s %14.10f %14.10f %14.10f\n", (const char*)((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_sName, ts->m_vaCoords_Original[i][0] / 100.0f, ts->m_vaCoords_Original[i][1] / 100.0f, ts->m_vaCoords_Original[i][2] / 100.0f);
			}
		}
	} else {
		CxIntArray sortIndex;
		sortIndex.SetSize(g_sortWannierHistory.GetSize());
		
		// Hungarian algorithm as described on http://de.wikipedia.org/wiki/Ungarische_Methode, 07.02.14
		const double thresh = 1.0e-5;
		int i, j;
		for(i = 0; i < g_sortWannierCount; i++) {
			for(j = 0; j < g_sortWannierCount; j++) {
				CxDVector3 dist = g_sortWannierHistory[i] - ts->m_vaCoords[g_sortWannierHistoryIndex[j]];
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
				g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] = dist.GetLengthSqr();
			}
		}
		for(i = 0; i < g_sortWannierCount; i++) {
			g_sortWannierKMStarCol[i] = -1;
		}
		
		for(j = 0; j < g_sortWannierCount; j++) {
			g_sortWannierKMMinCol[j] = g_sortWannierDistanceMatrix[j];
			for(i = 0; i < g_sortWannierCount; i++) {
				if(g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] < g_sortWannierKMMinCol[j]) {
					g_sortWannierKMMinCol[j] = g_sortWannierDistanceMatrix[i*g_sortWannierCount+j];
				}
			}
		}
		for(i = 0; i < g_sortWannierCount; i++) {
			g_sortWannierKMMinRow[i] = g_sortWannierDistanceMatrix[i*g_sortWannierCount] - g_sortWannierKMMinCol[0];
			for(j = 0; j < g_sortWannierCount; j++) {
				if(g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinCol[j] < g_sortWannierKMMinRow[i]) {
					g_sortWannierKMMinRow[i] = g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinCol[j];
				}
			}
		}
		
		while(true) {
			for(i = 0; i < g_sortWannierCount; i++) {
				g_sortWannierKMPrimeRow[i] = -1;
			}
			
			int count = 0;
			for(i = 0; i < g_sortWannierCount; i++) {
				bool star = false;
				for(j = 0; j < g_sortWannierCount; j++) {
					if(g_sortWannierKMStarCol[j] == i) {
						star = true;
						count++;
						break;
					}
				}
				if(star) {
					continue;
				}
				for(j = 0; j < g_sortWannierCount; j++) {
					if(g_sortWannierKMStarCol[j] == -1 && g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinRow[i] - g_sortWannierKMMinCol[j] < thresh) {
						g_sortWannierKMStarCol[j] = i;
						count++;
						break;
					}
				}
			}
			if(count == g_sortWannierCount) {
				break;
			}
			
			int row;
			while(true) {
				double min = 1.0e5;
				for(i = 0; i < g_sortWannierCount; i++) {
					if(g_sortWannierKMPrimeRow[i] != -1) {
						continue;
					}
					for(j = 0; j < g_sortWannierCount; j++) {
						if((g_sortWannierKMStarCol[j] == -1 || g_sortWannierKMPrimeRow[g_sortWannierKMStarCol[j]] != -1) && g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinRow[i] - g_sortWannierKMMinCol[j] < min) {
							min = g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinRow[i] - g_sortWannierKMMinCol[j];
						}
					}
				}
				
				for(i = 0; i < g_sortWannierCount; i++) {
					if(g_sortWannierKMStarCol[i] == -1 || g_sortWannierKMPrimeRow[g_sortWannierKMStarCol[i]] != -1) {
						g_sortWannierKMMinCol[i] += min;
					}
					if(g_sortWannierKMPrimeRow[i] != -1) {
						g_sortWannierKMMinRow[i] -= min;
					}
				}
				
				row = -1;
				for(i = 0; i < g_sortWannierCount; i++) {
					if(g_sortWannierKMPrimeRow[i] != -1) {
						continue;
					}
					bool found = false;
					for(j = 0; j < g_sortWannierCount; j++) {
						if((g_sortWannierKMStarCol[j] == -1 || g_sortWannierKMPrimeRow[g_sortWannierKMStarCol[j]] != -1) && g_sortWannierDistanceMatrix[i*g_sortWannierCount+j] - g_sortWannierKMMinRow[i] - g_sortWannierKMMinCol[j] < thresh) {
							g_sortWannierKMPrimeRow[i] = j;
							row = i;
							found = true;
							break;
						}
					}
					if(found) {
						break;
					}
				}
				
				bool star = false;
				for(j = 0; j < g_sortWannierCount; j++) {
					if(g_sortWannierKMStarCol[j] == row) {
						star = true;
						break;
					}
				}
				if(!star) {
					break;
				}
			}
			
			while(true) {
				if(g_sortWannierKMStarCol[g_sortWannierKMPrimeRow[row]] == -1) {
					break;
				}
				int tmp = g_sortWannierKMStarCol[g_sortWannierKMPrimeRow[row]];
				g_sortWannierKMStarCol[g_sortWannierKMPrimeRow[row]] = row;
				row = tmp;
			}
			g_sortWannierKMStarCol[g_sortWannierKMPrimeRow[row]] = row;
		}
		
		for(i = 0; i < g_sortWannierHistoryIndex.GetSize(); i++)
			sortIndex[g_sortWannierKMStarCol[i]] = g_sortWannierHistoryIndex[i];
		
		fprintf(g_sortWannierTrajFile, "%d\n", g_iGesAtomCount);
		if(ts->m_pComment != NULL)
			fprintf(g_sortWannierTrajFile, "%s\n", ts->m_pComment);
		else
			fprintf(g_sortWannierTrajFile, "\n");
		
		int pos = 0;
		for(i = 0; i < g_iGesAtomCount; i++) {
			if(g_baAtomIndex[i] == g_iWannierAtomType) {
				g_sortWannierHistory[pos] = ts->m_vaCoords[sortIndex[pos]];
				fprintf(g_sortWannierTrajFile, "%4s %14.10f %14.10f %14.10f\n", (const char*)((CAtom *)g_oaAtoms[g_iWannierAtomType])->m_sName, ts->m_vaCoords_Original[sortIndex[pos]][0] / 100.0f, ts->m_vaCoords_Original[sortIndex[pos]][1] / 100.0f, ts->m_vaCoords_Original[sortIndex[pos]][2] / 100.0f);
				pos++;
			} else {
				fprintf(g_sortWannierTrajFile, "%4s %14.10f %14.10f %14.10f\n", (const char*)((CAtom *)g_oaAtoms[g_waAtomRealElement[i]])->m_sName, ts->m_vaCoords_Original[i][0] / 100.0f, ts->m_vaCoords_Original[i][1] / 100.0f, ts->m_vaCoords_Original[i][2] / 100.0f);
			}
		}
	}
}

void finalizeSortWannier() {
	fclose(g_sortWannierTrajFile);
}
