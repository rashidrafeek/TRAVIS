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


#ifndef IONPAIR_H
#define IONPAIR_H


// This must always be the first include directive
#include "config.h"

#include "moltools.h"


bool gatherIonPair();
bool initializeIonPair();
bool processIonPair(CTimeStep* ts);
bool finalizeIonPair();


class CIonPairObservation: public CObservation
{
public:
	CIonPairObservation();
	~CIonPairObservation();

	CxIntArray ia_anion, ia_cation;	// Contains the type numbers of the ion types
	CxObArray *m_oaAnion;		// Contains an ObArray for each anion. This array contains the CIons 
	CxObArray *m_oaCation; 		//	
	
	CxObArray *m_oaNeigh;		// Arrays are two dimensional: First, all cations, second all anions
	CxObArray *m_oaPair;		// The ions are thereby sorted in an "A1,A2,A3,B1,B2,C1..." style

	CxObArray *m_oaLastNeigh;	// Contains the last neighborstates
	CxIntArray *m_laLastPair;	// Contains the last pair partner (there must be one in every step)

	CxObArray* PairDynamic; 	// Final Datasets after processing
	CxObArray* DiffusionDynamic;
	CxObArray* k_t; 		// Numerical derivation of BondDynamic
	CxObArray* a;		
	CxObArray* b;

	double m_fDist;			// Distance for the caging

	int m_iDepth, m_iTrunc;
	int m_iTotalCats, m_iTotalAns;
	double m_dGraceMin, m_dGraceMax;

	void initialize();
	void process(CTimeStep* ts);
	void finalize();

	void Initiate_Offsets(int, CxObArray*);
	bool checkAtomChoice(CMolecule*, const char*, CxObArray*);
	bool includeAtom(int, int, int, CxObArray*);

	void RestoreFunction(CxIntArray*, CxDoubleArray*);
	void CalcMinMax(CxDoubleArray*, CxDoubleArray*, CxDoubleArray*, CxDoubleArray*, CxDoubleArray*);

	void WriteOutput(CxDoubleArray*, CxDoubleArray*, CxDoubleArray*, CxDoubleArray*, CxDoubleArray*,std::string);

	
private:

};


class CIon: public CxObject
{
public:
	CIon();
	~CIon();
	
	int m_iOffset;
	double m_dCorr;

private:

};

#endif

