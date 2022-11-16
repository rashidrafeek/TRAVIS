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


#ifndef REVSDF_H
#define REVSDF_H


// This must always be the first include directive
#include "config.h"

#include "xobject.h"
#include "2df.h"
#include "xdvec3array.h"
#include "atomgroup.h"
#include "3df.h"


class CSingleMolecule;


class CRevSDF : public CxObject
{
public:
	CRevSDF();
	~CRevSDF();

	void BuildAtomList(CSingleMolecule *ref, CSingleMolecule *obs, CxIntArray *vec);
	void BuildName();
	void Parse();
	void CreateRevSDF();

	bool m_bIntra;
	int m_iRefOrSec;
	char *m_sName;
	double m_fRadius;
	int m_iResolution;
	bool m_bMirrorY;
	bool m_bMirrorBond;
	int m_iShowMol;
	CxDVec3Array *m_vaData; 
	C2DF *m_p2DF;
	C3DF<double> *m_pRevSDF;
	CAtomGroup m_oAtoms;
	double m_fParticleDensity;
	int m_iShowAtomGes;

	double m_fSecondAtomPosX;
	double m_fSecondAtomCount;

	bool m_bCorrectAngle;
	bool m_bCorrectRadial;
	bool m_bCreateRevSDF;
	int m_iRevSDFRes;
	bool m_bDrawAtoms;
	int m_iHistogramRes;
	CxByteArray *m_baDataEnabled; 
};


#endif
