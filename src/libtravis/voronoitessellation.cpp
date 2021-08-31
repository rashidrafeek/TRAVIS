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


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "voronoitessellation.h"

#include "cellinfo.h"
#include "matrix.h"
#include "vector.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Libtravis {

namespace Travis {

std::vector<std::size_t> VoronoiTessellation::Cell::getFaceIdsByNeighborId(std::size_t neighborId) const {
	std::vector<std::size_t> faceIds;
	for (auto &&f: m_faceIds)
		if (m_voronoiTessellation->m_faces[f].getCellIdFirst() == neighborId || m_voronoiTessellation->m_faces[f].getCellIdSecond() == neighborId)
			faceIds.push_back(f);
	return faceIds;
}

void VoronoiTessellation::Face::createLocalCoordinates(const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, double epsilon) {
	const std::size_t size = m_vertexIds.size();
	if (size < 2)
		throw std::runtime_error("VoronoiTessellation::Face: Cannot create local coordinates with less than two vertices");
	
	Vector3 centroid(0.0, 0.0, 0.0);
	std::vector<Vector2> vertexCoordinates(size);
	for (std::size_t i = 0; i < size; ++i) {
		Vector3 vec = cellInfo.vectorCartesianToFractional(m_voronoiTessellation->m_vertices[m_vertexIds[i]]);
		vertexCoordinates[i][0] = vec[0] * gridSizes[0];
		vertexCoordinates[i][1] = vec[1] * gridSizes[1];
		if (vertexCoordinates[i][0] < m_localCoordinates.minimumA)
			m_localCoordinates.minimumA = vertexCoordinates[i][0];
		if (vertexCoordinates[i][0] > m_localCoordinates.maximumA)
			m_localCoordinates.maximumA = vertexCoordinates[i][0];
		if (vertexCoordinates[i][1] < m_localCoordinates.minimumB)
			m_localCoordinates.minimumB = vertexCoordinates[i][1];
		if (vertexCoordinates[i][1] > m_localCoordinates.maximumB)
			m_localCoordinates.maximumB = vertexCoordinates[i][1];
		m_localCoordinates.localCentroid += vertexCoordinates[i];
		centroid += m_voronoiTessellation->m_vertices[m_vertexIds[i]];
	}
	m_localCoordinates.minimumA -= epsilon;
	m_localCoordinates.maximumA += epsilon;
	m_localCoordinates.minimumB -= epsilon;
	m_localCoordinates.maximumB += epsilon;
	centroid /= size;
	m_localCoordinates.localCentroid /= size;
	for (std::size_t i = 0; i < size; ++i)
		vertexCoordinates[i] -= m_localCoordinates.localCentroid;
	
	Vector3 normal = m_voronoiTessellation->m_vertices[m_vertexIds[0]] - centroid;
	normal = normal.cross(m_voronoiTessellation->m_vertices[m_vertexIds[1]] - centroid);
	normal.normalize();
	m_localCoordinates.normalCentroid = normal.dot(centroid);
	m_localCoordinates.normal[0] = normal.dot(cellInfo.getVectorA() / gridSizes[0]);
	m_localCoordinates.normal[1] = normal.dot(cellInfo.getVectorB() / gridSizes[1]);
	m_localCoordinates.normal[2] = normal.dot(cellInfo.getVectorC() / gridSizes[2]);
	
	m_localCoordinates.edgeVectors.resize(size);
	m_localCoordinates.triangleAreas.resize(size);
	for (std::size_t i = 0; i < size - 1; ++i) {
		m_localCoordinates.edgeVectors[i] = vertexCoordinates[i + 1] - vertexCoordinates[i];
		m_localCoordinates.triangleAreas[i] = vertexCoordinates[i][0] * vertexCoordinates[i + 1][1] - vertexCoordinates[i + 1][0] * vertexCoordinates[i][1];
	}
	m_localCoordinates.edgeVectors[size - 1] = vertexCoordinates[0] - vertexCoordinates[size - 1];
	m_localCoordinates.triangleAreas[size - 1] = vertexCoordinates[size - 1][0] * vertexCoordinates[0][1] - vertexCoordinates[0][0] * vertexCoordinates[size - 1][1];
	
	if (m_directNeighbor) {
		double cell1 = normal.dot(m_voronoiTessellation->m_cells[m_cellIdFirst].getCentroid() - centroid);
		double cell2 = normal.dot(m_voronoiTessellation->m_cells[m_cellIdSecond].getCentroid() - centroid);
		if (std::signbit(cell1) != std::signbit(cell2)) {
			if (std::signbit(cell1) == std::signbit(m_localCoordinates.normal[2])) {
				m_localCoordinates.cellIdForward = m_cellIdFirst;
				m_localCoordinates.cellIdBackward = m_cellIdSecond;
			} else {
				m_localCoordinates.cellIdForward = m_cellIdSecond;
				m_localCoordinates.cellIdBackward = m_cellIdFirst;
			}
		} else {
			throw std::runtime_error("VoronoiTessellation::Face: Direct neighbors on same side of face");
		}
	} else {
		double cell1 = normal.dot(m_voronoiTessellation->m_cells[m_cellIdFirst].getCentroid() - centroid);
		if (std::signbit(cell1) == std::signbit(m_localCoordinates.normal[2])) {
			m_localCoordinates.cellIdForward = m_cellIdFirst;
			m_localCoordinates.cellIdBackward = std::numeric_limits<std::size_t>::max();
		} else {
			m_localCoordinates.cellIdForward = std::numeric_limits<std::size_t>::max();
			m_localCoordinates.cellIdBackward = m_cellIdFirst;
		}
	}
}

void VoronoiTessellation::Face::setRayCandidates(std::vector<std::vector<std::size_t>> &candidateList, long int minimumA, long int maximumA, long int minimumB, long int maximumB, bool refining) const {
	const long int locMinA = std::min(static_cast<long int>(refining ? std::floor(m_localCoordinates.minimumA) : std::ceil(m_localCoordinates.minimumA)), maximumA);
	const long int locMaxA = std::max(static_cast<long int>(std::floor(m_localCoordinates.maximumA)), minimumA);
	const long int locMinB = std::min(static_cast<long int>(refining ? std::floor(m_localCoordinates.minimumB) : std::ceil(m_localCoordinates.minimumB)), maximumB);
	const long int locMaxB = std::max(static_cast<long int>(std::floor(m_localCoordinates.maximumB)), minimumB);
	const std::size_t iMax = static_cast<std::size_t>(locMaxA - minimumA);
	const std::size_t jMax = static_cast<std::size_t>(locMaxB - minimumB);
	const std::size_t sizeB = static_cast<std::size_t>(maximumB - minimumB + 1);
	for (std::size_t i = static_cast<std::size_t>(locMinA - minimumA); i <= iMax; ++i) {
		const std::size_t indexA = i * sizeB;
		for (std::size_t j = static_cast<std::size_t>(locMinB - minimumB); j <= jMax; ++j) {
			candidateList[indexA + j].push_back(m_id);
		}
	}
}

bool VoronoiTessellation::Face::calculateIntersectionRayC(double coordA, double coordB, double &coordC, double epsilon) const {
	if (std::fabs(m_localCoordinates.normal[2]) < epsilon)
		return false;
	
	bool tooClose = true;
	double localIntersection[2];
	unsigned int loopCount = 0;
	const std::size_t numTriangles = m_localCoordinates.triangleAreas.size();
	while (tooClose) {
		tooClose = false;
		localIntersection[0] = coordA - m_localCoordinates.localCentroid[0];
		localIntersection[1] = coordB - m_localCoordinates.localCentroid[1];
		double area;
		std::size_t i;
		for (i = 0; i < numTriangles; ++i) {
			area = localIntersection[1] * m_localCoordinates.edgeVectors[i][0] - localIntersection[0] * m_localCoordinates.edgeVectors[i][1] + m_localCoordinates.triangleAreas[i];
			if (std::fabs(area) < epsilon) {
				tooClose = true;
				break;
			}
			if (std::signbit(area) != std::signbit(m_localCoordinates.triangleAreas[i]))
				return false;
		}
		if (tooClose) {
			if (m_localCoordinates.edgeVectors[i][1] > 0) {
				coordA += epsilon * m_localCoordinates.edgeVectors[i][1];
				coordB -= epsilon * m_localCoordinates.edgeVectors[i][0];
			} else {
				coordA -= epsilon * m_localCoordinates.edgeVectors[i][1];
				coordB += epsilon * m_localCoordinates.edgeVectors[i][0];
			}
		}
		++loopCount;
		if (loopCount > 100)
			throw std::runtime_error("VoronoiTessellation::Face: Loop count limit exceeded");
	}
	
	coordC = (m_localCoordinates.normalCentroid - coordA * m_localCoordinates.normal[0] - coordB * m_localCoordinates.normal[1]) / m_localCoordinates.normal[2];
	if (!std::isfinite(coordC))
		return false;
	
	return true;
}

}

}

#endif


