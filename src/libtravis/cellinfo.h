/*****************************************************************************

    libTRAVIS - Class library for trajectory analysis and visualization

    Copyright (C) 2015-2021 Martin Thomas

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


#ifndef LT_CELLINFO_H
#define LT_CELLINFO_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "matrix.h"
#include "vector.h"

namespace Libtravis {

namespace Travis {

class CellInfo
{
public:
	CellInfo() {}
	CellInfo(const MatrixS3 &cellMatrix, bool periodicA = true, bool periodicB = true, bool periodicC = true) : m_cellMatrix(cellMatrix), m_inverseCellMatrix(m_cellMatrix.inverseMatrix()), m_volume(m_cellMatrix.determinant()), m_periodicA(periodicA), m_periodicB(periodicB), m_periodicC(periodicC) {}
	CellInfo(const Vector3 &vectorA, const Vector3 &vectorB, const Vector3 &vectorC, bool periodicA = true, bool periodicB = true, bool periodicC = true) : CellInfo({vectorA[0], vectorB[0], vectorC[0], vectorA[1], vectorB[1], vectorC[1], vectorA[2], vectorB[2], vectorC[2]}, periodicA, periodicB, periodicC) {}
	CellInfo(double orthoA, double orthoB, double orthoC, bool periodicA = true, bool periodicB = true, bool periodicC = true) : CellInfo({orthoA, 0.0, 0.0, 0.0, orthoB, 0.0, 0.0, 0.0, orthoC}, periodicA, periodicB, periodicC) {}
	CellInfo(double cubicA, bool periodicA = true, bool periodicB = true, bool periodicC = true) : CellInfo({cubicA, 0.0, 0.0, 0.0, cubicA, 0.0, 0.0, 0.0, cubicA}, periodicA, periodicB, periodicC) {}
	
	void setOrigin(const Vector3 &origin) { m_origin = origin; }
	const Vector3 &getOrigin() const { return m_origin; }
	double getVolume() const { return m_volume; }
	
	const MatrixS3 &getCellMatrix() const { return m_cellMatrix; }
	const MatrixS3 &getInverseCellMatrix() const { return m_inverseCellMatrix; }
	Vector3 getVectorA() const { return Vector3(m_cellMatrix[0], m_cellMatrix[3], m_cellMatrix[6]); }
	Vector3 getVectorB() const { return Vector3(m_cellMatrix[1], m_cellMatrix[4], m_cellMatrix[7]); }
	Vector3 getVectorC() const { return Vector3(m_cellMatrix[2], m_cellMatrix[5], m_cellMatrix[8]); }
	bool isPeriodicA() const { return m_periodicA; }
	bool isPeriodicB() const { return m_periodicB; }
	bool isPeriodicC() const { return m_periodicC; }
	bool isFullyPeriodic() const { return m_periodicA && m_periodicB && m_periodicC; }
	bool isOrthorhombic() const { return m_cellMatrix[1] == 0.0 && m_cellMatrix[2] == 0.0 && m_cellMatrix[3] == 0.0 && m_cellMatrix[5] == 0.0 && m_cellMatrix[6] == 0.0 && m_cellMatrix[7] == 0.0; }
	
	Vector3 foldVector(const Vector3 &vector) const;
	Vector3 vectorCartesianToFractional(const Vector3 &vector) const { return m_inverseCellMatrix * (vector - m_origin); }
	Vector3 vectorFractionalToCartesian(const Vector3 &vector) const { return m_cellMatrix * vector + m_origin; }
	
private:
	MatrixS3 m_cellMatrix = MatrixS3::unityMatrix();
	MatrixS3 m_inverseCellMatrix = MatrixS3::unityMatrix();
	Vector3 m_origin{0.0, 0.0, 0.0};
	double m_volume = 1.0;
	bool m_periodicA = false;
	bool m_periodicB = false;
	bool m_periodicC = false;
};

}

}

#endif

#endif


