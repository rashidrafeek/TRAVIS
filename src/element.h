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


#ifndef ELEMENT_H
#define ELEMENT_H


// This must always be the first include directive
#include "config.h"

#include "xobject.h"


class CElement : public CxObject  
{
public:
	CElement()
	{
		m_sLabel = NULL;
		m_iColorR = 0;
		m_iColorG = 255;
		m_iColorB = 255;
		m_iColorBleachedR = 170;
		m_iColorBleachedG = 170;
		m_iColorBleachedB = 170;
		m_fMass = 0;
		m_fRadius = 0;
		m_iOrd = 0;
		m_fVdWRadius = 0;
	}

	virtual ~CElement()
	{
		if (m_sLabel != NULL)
		{
			delete[] m_sLabel;
			m_sLabel = NULL;
		}
	}

	void CopyData(CElement *e)
	{
		m_fMass = e->m_fMass;
		m_iOrd = e->m_iOrd;
		m_fRadius = e->m_fRadius;
		m_iColorR = e->m_iColorR;
		m_iColorG = e->m_iColorG;
		m_iColorB = e->m_iColorB;
		m_iColorBleachedR = e->m_iColorBleachedR;
		m_iColorBleachedG = e->m_iColorBleachedG;
		m_iColorBleachedB = e->m_iColorBleachedB;
	}

	char *m_sLabel;
	double m_fMass;
	double m_fRadius;
	double m_fVdWRadius;
	int m_iOrd;
	double m_fCoherentNCS;

	unsigned char m_iColorR;
	unsigned char m_iColorG;
	unsigned char m_iColorB;
	unsigned char m_iColorBleachedR;
	unsigned char m_iColorBleachedG;
	unsigned char m_iColorBleachedB;

	unsigned char ColorR(double bleached)
	{
		return (unsigned char)(((double)m_iColorBleachedR-m_iColorR)*bleached+m_iColorR);
	}

	unsigned char ColorG(double bleached)
	{
		return (unsigned char)(((double)m_iColorBleachedG-m_iColorG)*bleached+m_iColorG);
	}

	unsigned char ColorB(double bleached)
	{
		return (unsigned char)(((double)m_iColorBleachedB-m_iColorB)*bleached+m_iColorB);
	}

};

#endif 
