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


#ifndef BICGSTAB_H
#define BICGSTAB_H


// This must always be the first include directive
#include "config.h"

#include "3df.h"
#include <stdio.h>



class CxDoubleArray;



class CSparseMatrix {
public:
	friend class CCurrentPDEDiscretizer;
	friend class CCurrentPDESolver;
	friend class CILUPreconditioner;
	
	CSparseMatrix();
	CSparseMatrix(const CSparseMatrix &matrix);
	~CSparseMatrix();
	
	void setTest();
	void matrixVectorProduct(const double *vect, double *result) const;
	
private:
	int m_size;
	double *m_val;
	int *m_colInd;
	int *m_rowPtr;
};



class CCurrentPDEDiscretizer {
public:
	static void discretize(CSparseMatrix *pdeMatrix, const C3DF<VORI_FLOAT> &electronDensity, const C3DF<VORI_FLOAT> &electronDensityGradientX, const C3DF<VORI_FLOAT> &electronDensityGradientY, const C3DF<VORI_FLOAT> &electronDensityGradientZ);
};



class CCurrentPDESolver {
public:
	static bool bicgstabl(const int l, const CSparseMatrix *pdeMatrix, double *solution, const double *rightHandSide, const int maxIter, double *thresh, FILE *infoFile, bool fastmode=false);
	static double calcResidual(const CSparseMatrix *pdeMatrix, const double *solution, const double *rightHandSide);
	
private:
	static void luFactorization(double *matrix, int *perm, int size);
	static void forwardBackward(double *luMatrix, int *perm, int size, double *sol);
};



class CILUPreconditioner {
public:
	static void precondition(CSparseMatrix *matrix);
	static void solve(CSparseMatrix *matrix, double *vect);
};



#endif
