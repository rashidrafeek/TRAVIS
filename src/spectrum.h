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


#ifndef SPECTRUM_H
#define SPECTRUM_H


// This must always be the first include directive
#include "config.h"

#include "xobject.h"
#include "tools.h"
#include "backtrace.h"


class CSpectrum : public CxObject  
{
public:
	void SetIntegral(double f);
	void SetMaxValue(double f);
	void Multiply(double f);
	void MakeDB();
	void SetMaxRWL(double f);
	void Write(const char *pre, const char *s, const char *post);
	void FromComplex(double *f);
	void Create(int i);
	double *m_pComplex;
	double *m_pData;
	double m_fMaxRWL;
	long m_iSize;
	double m_fWaveNumber;
	CSpectrum();
	virtual ~CSpectrum();

};

#endif 
