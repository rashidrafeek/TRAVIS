/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

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


#ifndef REFLUX_H
#define REFLUX_H


// This must always be the first include directive
#include "config.h"

#include "moltools.h"


bool gatherReFlux();
bool initializeReFlux();
bool processReFlux(CTimeStep* ts);
bool finalizeReFlux();

class CReFluxSingleCond: public CxObject
{
public:
	CReFluxSingleCond();
	~CReFluxSingleCond();

	CxIntArray m_iaTypes;		// Molecule types
	int m_iRefDist;			// Ref for Dist
	int m_iRefVec;			// Reference vectors base for angle crit
	int m_iRefVec2;			// Reference vectors tipp
	int m_iRefVec3;			// only used if Vectors are defined by plane
	int m_iObsVec;			// Observed vectors base for angle crit of Obs for Dist
	int m_iObsVec2;			// Observed vectors tipp
	int m_iObsVec3;			// only used if Vectors are defined by plane
	double m_dAggrMinAngle;		// Angle for AngleCrit of Dist for Dist
	double m_dAggrMaxAngle;		// Angle for AngleCrit

	void Copy(CReFluxSingleCond*);
};

class CReFluxSubCond: public CxObject
{
public:
	CReFluxSubCond();
	~CReFluxSubCond();

	CxIntArray m_iaDistConds;	// all of these must be fulfilled simultaneously!
	CxIntArray m_iaAngleConds;	// all of these must be fulfilled simultaneously!
	CxObArray* m_oaDistConds;
	CxObArray* m_oaAngleConds;

	bool check(CTimeStep*);		// checks fulfillment
};

class CReFluxCond: public CxObject	// Condition for each pair
{
public:
	CReFluxCond();
	~CReFluxCond();

	int m_iRefNeigh; 		// Reference for neighbor criterion
	int m_iObsNeigh; 		// Observed for neighbor criterion
	double m_dNeighDist;		// Distance for Neighbor

	CxObArray m_oaAggrConds;	// at least one of these must be fulfilled
	CxObArray m_oaDistConds;	// single conds (the AggrConds link on them)
	CxObArray m_oaAngleConds;	// single conds (the AggrConds link on them)

	void check(CTimeStep*,CxIntArray*,CxIntArray*);	// checks conditions for fulfillment
	void Add(int, int, int);
	void AddArray(int, int, CxIntArray*);
};

class CReFluxObservation: public CObservation
{
public:
	CReFluxObservation();
	~CReFluxObservation();

	CxObArray m_oaCond;		// Contains the conditions for the possible pairs

	CxObArray m_oaNeigh;		// Contains the switching timesteps 
	CxObArray m_oaAggr;

	CxDoubleArray AggrDynamic; // Final Datasets after processing
	CxDoubleArray DiffusionDynamic;
	CxDoubleArray k_t; 	// Numerical derivation of BondDynamic
	CxDoubleArray a;		
	CxDoubleArray b;

	int m_iShowMol3;
	bool m_bRefVecPlane;		// defined by plane?
	bool m_bObsVecPlane;		// defined by plane?
	CxIntArray m_iaAngleMols;
	CxIntArray tempRef, tempObs, tempRef2, tempObs2, tempRef3, tempObs3;
	int m_iDepth, m_iTrunc;
	double m_dNeighDist, m_dAggrDist;
	double m_dGraceMin, m_dGraceMax;

	void initialize();
	void process(CTimeStep* ts);
	void finalize();

	bool checkAtomChoice(const char*i, int, int);
	bool includeAtom(int, int, int, int);
	bool GetArray(const char*, CxIntArray*, int);
	bool GetDoubleArray( const char*, CxIntArray*, CxIntArray*, int, int );
	CxString FindType(CReFluxSingleCond*, int);
	void DumpDist();				// Dumps the Conds	
	void DumpAngle();				// Dumps the Conds	
	void DumpAll();					// Dumps the Conds	
	void SortDist();				// Sorts the Conds 
	void SortAngle();				// Sorts the Conds 
	void SortFinal();				// Sorts the Conds 
	void RestoreFunction(CxIntArray*, CxDoubleArray*);
	void CalcMinMax();
	void WriteOutput();
	
private:

};


#endif

