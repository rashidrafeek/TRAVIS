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


#ifndef INTERNALCOORD_H
#define INTERNALCOORD_H


// This must always be the first include directive
#include "config.h"

#include "xobject.h"
#include "xobarray.h"
#include "xdoublearray.h"


class CMolBond : public CxObject
{
public:
	CMolBond();
	~CMolBond();
	int m_iAtomType[2];
	int m_iAtom[2];
	int m_iMolAtom[2];
	int m_iAtomOffset[2];
	CxDoubleArray m_faData;
};


class CMolBondGroup : public CxObject
{
public:
	CMolBondGroup();
	~CMolBondGroup();
	CxObArray m_oaBonds;
};


class CMolAngle : public CxObject
{
public:
	CMolAngle();
	~CMolAngle();
	int m_iAtomType[3];
	int m_iAtom[3];
	int m_iMolAtom[3];
	int m_iAtomOffset[3];
	CxDoubleArray m_faData;
};


class CMolAngleGroup : public CxObject
{
public:
	CMolAngleGroup();
	~CMolAngleGroup();
	CxObArray m_oaAngles;
};



#endif

