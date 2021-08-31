/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

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


#ifndef NORMALCOORDINATE_H
#define NORMALCOORDINATE_H


// This must always be the first include directive
#include "config.h"

#include "moltools.h"
#include "xdoublearray.h"
#include "xintarray.h"
#include "xobject.h"
#include "xobarray.h"
#include "xdvector3.h"


class CTimeStep;
class CxByteArray;
class CxDVec3Array;

class CNormalCoordinateObservation;


class CReferenceStructure: public CxObject
{
public:
	CReferenceStructure(int showMol, const char *basename, bool calcIR);
	~CReferenceStructure();
	
	void initialize(const char *basename);
	void processCoordinates(CxDVec3Array &coord, CxDVector3 &dipole, int showMol);
	void calcVelocities(int showMol);
	double calcMinimumDistance(int showMol, bool useInternals, CxObArray &internals);
	void nextStep();
	double getMinimumDistance(int showMol, int step);
	void setGlobalProbability(int showMol, int step, double prob);
	void finalize(const char *basename);
	
private:
	CTimeStep *_refTimestep;
	CSingleMolecule *_singleMol;
	int _showMolCount;
	int _atomCount;
	int _permutationCount;
	int _numSteps;
	
	bool _writeTransformedTrajectories;
	FILE **_transformedTrajectoryFiles;
	
	CxObArray _permutations;
	CxDVector3 _centroid;
	CxObArray _centroidCoords;
	CxDoubleArray _weights;
	double _weightSum;
	char _filename[128];
	
	CxObArray _distanceTimedev;
	CxObArray _coordHistory;
	int _historyIndex;
	bool _calcVelocity;
	CxObArray _velocityCache;
	CxObArray _globalProb;
	
	bool _useInternals;
	CxObArray _internals;
	double _gaussianWidth;
	
	int _correlationDepth;
	int _windowFunction;
	int _windowFunctionParameter;
	double _specResolution;
	int _specSize;
	int _zeroPadding;
	
	double _convergenceThreshold;
	int _maxIterations;
	
	bool _calcIR;
	CxObArray _dipoleHistory;
	CxObArray _dipoleDerivativeCache;
	
	bool recognizeMolecule(CTimeStep *ts, int showMol, const char *basename);
	void recognizeMoleculeRecursion(unsigned int index, bool *used, int depth, unsigned int *stack, CTimeStep *ts, CxByteArray &baAtomIndex);
	bool recognizeMoleculeBondRange(CTimeStep *ts, int i1, int i2, double *distance, CxByteArray &baAtomIndex);
	void askPermutations(CxIntArray *actions);
	void createPermutationsRecursion(int start, CxIntArray *permutation, CxIntArray *permutationActions);
	
	double offDiagonalNorm(CxObArray &matrix);
	double findRotationAngle(CxObArray &matrix, int i, int j);
	void calcIntegrals(CxObArray &matrix, CxDoubleArray *integrals, CxDoubleArray *centers);
};

class CNormalCoordinateObservation: public CObservation
{
public:
	CNormalCoordinateObservation();
	~CNormalCoordinateObservation();
	
	void initialize();
	void process(CTimeStep *ts);
	void finalize();
	
private:
	CxObArray _refStructures;
	int _refCount;
	int _numSteps;
	double _gaussianWidth;
	//int _correlationDepth;
	
	CxObArray _distanceTimedev;
	bool _useInternals;
	CxObArray _internals;
	
	bool _calcIR;
	
// 	bool _writeTransformedTrajectories;
// 	FILE **_transformedTrajectoryFiles;
// 	CxObArray _distanceTimedev;
// 	CxObArray _coordHistory;
// 	int _historyIndex;
// 	bool _calcVelocity;
// 	CxObArray _velocityCache;
};

bool gatherNormalCoordinate();
bool initializeNormalCoordinate();
void processNormalCoordinate(CTimeStep *ts);
void finalizeNormalCoordinate();

class CEckartReferenceStructure: public CxObject
{
public:
	CEckartReferenceStructure();
	~CEckartReferenceStructure();
	
	void process(CTimeStep *ts);
	const double *rotationMatrix() { return _rotationMatrix; }
	const CxDVector3 &translationVector() { return _translationVector; }
	
private:
	CTimeStep *_refTimestep;
	CSingleMolecule *_singleMol;
	char _filename[128];
	int _singleMolIndex;
	CAtomGroup *_mapAtoms;
	CxDVector3 _refCentroid;
	CxDVec3Array _refCentroidCoord;
	CxDoubleArray _weights;
	double _weightSum;
	
	double _rotationMatrix[9];
	CxDVector3 _translationVector;
	
	bool recognizeMolecule(CTimeStep *ts);
	void recognizeMoleculeRecursion(unsigned int index, bool *used, int depth, unsigned int *stack, CTimeStep *ts, CxByteArray &baAtomIndex);
	bool recognizeMoleculeBondRange(CTimeStep *ts, int i1, int i2, double *distance, CxByteArray &baAtomIndex);
};

bool gatherEckartTransform();
bool initializeEckartTransform();
void processEckartTransform(CTimeStep *ts);
void finalizeEckartTransform();

#endif
