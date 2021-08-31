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


#ifndef ACF_H
#define ACF_H


// This must always be the first include directive
#include "config.h"

#include "xobject.h"
#include "backtrace.h"
#include "spectrum.h"
#include "fft.h"
#include "atomgroup.h"
#include "xintarray.h"


class CSingleMolecule;


class CACF : public CxObject  
{
public:
	bool m_bExcludeR0;
//	bool m_b2ndDerivative;
	int m_iDerivative;
	bool m_bDerivative;
	void Mirror(int i);
	int m_iShowAtomGes;
	void BuildAtomList(CSingleMolecule *obs, CxIntArray *vec);
	CAtomGroup m_oAtoms;
	int m_iShowMol;
	void BuildName();
	CACF();
	~CACF();
	void Window();
	void WriteACF(const char *pre, const char *s, const char *post);
	void Normalize();
	void NormalizeCached();
	void Multiply(double f);
	void MultiplyCached(double f);
	void Create();
	void Transform(CFFT *fft);
	void Parse();

	char* m_sName;
	double *m_pData;
	unsigned long *m_pCounter;
	int m_iSize;
	CSpectrum *m_pSpectrum;
//	int m_iStride;
	bool m_bWindowFunction;
	bool m_bSpectrum;
	bool m_bMassWeight;
	double m_fSpecWaveNumber;
	int m_iZeroPadding;
	int m_iZeroPadding0;
	int m_iMirror;
	bool m_bACF_DB;
	CxObArray m_oaCache;
	bool m_bDecomposeModes; 
	int m_iParticles;       
	CxObArray m_oaCCRMatrix;

	//bool m_bSplitCart;

};

#endif 
