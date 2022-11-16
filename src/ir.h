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


#ifndef IR_H
#define IR_H


// This must always be the first include directive
#include "config.h"

#include "moltools.h"
#include "xptrarray.h"

class CTimeStep;

class CPowerObservation: public CObservation
{
public:
	explicit CPowerObservation(bool global = false);
	~CPowerObservation();
	
	void initialize();
	void process(CTimeStep *ts);
	void finalize();
	
private:
	char *_name;
	CAtomGroup *_atoms;
	int _correlationDepth;
	int _windowFunction;
	int _windowFunctionParameter;
	int _zeroPadding;
	int _specSize;
	double _specResolution;
	bool _saveACF;
	bool _massWeighting;
	bool m_bSplitCart;
	int m_iFinalDepth;
	bool _correctfreq;

	
	CxObArray _velocityCache;
	CxDoubleArray _masses;
};

bool gatherPowerSpectrum();
bool initializePowerSpectrum();
void processPowerSpectrum(CTimeStep *ts);
void finalizePowerSpectrum();

class CIRObservation: public CObservation
{
public:
	explicit CIRObservation(bool global = false);
	~CIRObservation();
	
	void initialize();
	void process();
	void finalize();
	
private:
	char *_name;
	int _correlationDepth;
	int _windowFunction;
	int _windowFunctionParameter;
	int _zeroPadding;
	int _specSize;
	double _specResolution;
	bool _saveACF;
	bool _includeCross;
	bool _finiteDifferenceCorrection;
	bool _correctfreq;
	
	CxObArray _dipoleCache;
};

bool gatherIR();
bool initializeIR();
void processIR(CTimeStep *ts);
void finalizeIR();

class CDipoleRestartObservation: public CObservation
{
public:
	CDipoleRestartObservation();
	~CDipoleRestartObservation();

	CxPtrArray m_oaFiles;
	CxIntArray m_iaMols;
private:
};

bool gatherDipoleRestart();
bool initializeDipoleRestart();
void processDipoleRestart(CTimeStep *ts);
void finalizeDipoleRestart();

class CVCDObservation: public CObservation
{
public:
	explicit CVCDObservation(bool global = false);
	~CVCDObservation();
	
	void initialize();
	void process();
	void finalize();
	
private:
	char *_name;
	int _correlationDepth;
	int _windowFunction;
	int _windowFunctionParameter;
	int _zeroPadding;
	int _specSize;
	double _specResolution;
	bool _saveMoments;
	bool _saveACF;
	bool _finiteDifferenceCorrection;
	
	int _smoothWidth;
	
	CxObArray _electricDipoleCache;
	CxObArray _magneticDipoleCache;
};

bool gatherVCD();
bool initializeVCD();
void processVCD(CTimeStep *ts);
void finalizeVCD();

class CMagneticRestartObservation: public CObservation
{
public:
	CMagneticRestartObservation();
	~CMagneticRestartObservation();

	CxPtrArray m_oaFiles;
	CxIntArray m_iaMols;
private:
};

bool gatherMagneticDipoleRestart();
bool initializeMagneticDipoleRestart();
void processMagneticDipoleRestart(CTimeStep *ts);
void finalizeMagneticDipoleRestart();

bool gatherSortWannier();
bool initializeSortWannier();
void processSortWannier(CTimeStep *ts);
void finalizeSortWannier();

#endif
