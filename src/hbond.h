/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Sascha Gehrke.

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


#ifndef HBOND_H
#define HBOND_H


// This must always be the first include directive
#include "config.h"

#include "moltools.h"


bool gatherHBond();
bool initializeHBond();
bool processHBond(CTimeStep* ts);
bool finalizeHBond();


class CHBondObservation: public CObservation
{
public:
	CHBondObservation();
	~CHBondObservation();

	CxIntArray *m_laAccep;	 // Enthaelt m_iOffset fuer jedes Acceptoratom
	CxIntArray *m_laHydro;	 // 		   "                Wasserstoffatom
	CxIntArray *m_laDonor;	 // 		   "                Donoratom

	CxObArray *m_oaNeigh;
	CxObArray *m_oaBond;
	CxObArray *m_oaPartners;	// Contains relevant hydrogens for each pair

	CxIntArray *m_laLastNeigh;		// Contains the last states 
	CxIntArray *m_laLastBond;		// Contains the last states

	CxDoubleArray* BondDynamic; // Final Datasets after processing
	CxDoubleArray* DiffusionDynamic;
	CxDoubleArray* k_t; 	// Numerical derivation of BondDynamic
	CxDoubleArray* a;		
	CxDoubleArray* b;



	bool m_bSameAtom;
	bool m_bLuzar;
	bool m_bHysteresis;
	int m_iDepth, m_iTrunc;
	double m_frAD, m_frAH, m_fwinkel;
	double m_frAD_out, m_frAH_out, m_fwinkel_out;
	double m_dGraceMin, m_dGraceMax;

	void initialize();
	void process(CTimeStep* ts);
	void finalize();

	bool checkCOND(CxDVector3,CxDVector3,CxDVector3,int);
	bool checkAtomChoice(CMolecule*, const char*, CxIntArray*);
	bool includeAtom(int, int, int, CxIntArray*);
	void BuildDonorArray();
	void BuildHydroArray();
	void FindHydro(CxIntArray*, int);
	void RestoreFunction(CxIntArray*, CxDoubleArray*);
//	void CheckBoundary(CxDVector3*);
	void CalcMinMax();

	void WriteOutput();
	
private:

};


#endif

