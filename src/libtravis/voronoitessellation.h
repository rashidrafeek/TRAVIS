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


#ifndef LT_VORONOITESSELLATION_H
#define LT_VORONOITESSELLATION_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "matrix.h"
#include "vector.h"

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace Libtravis {

namespace Travis {

class CellInfo;

class VoronoiTessellation
{
public:
	friend class CalcVoronoi;
	friend class CalcVoronoiIntegrals;
	
	class Face
	{
	public:
		friend class CalcVoronoi;
		friend class VoronoiTessellation;
		
		std::size_t getId() const { return m_id; }
		std::size_t getNumVertices() const { return m_vertexIds.size(); }
		std::size_t getVertexId(std::size_t localId) const { return m_vertexIds.at(localId); }
		const Vector3 &getVertex(std::size_t localId) const { return m_voronoiTessellation->m_vertices[getVertexId(localId)]; }
		std::size_t getCellIdFirst() const { return m_cellIdFirst; }
		std::size_t getCellIdSecond() const { return m_cellIdSecond; }
		bool isDirectNeighbor() const { return m_directNeighbor; }
		
		void createLocalCoordinates(const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, double epsilon);
		void setRayCandidates(std::vector<std::vector<std::size_t>> &candidateList, long int minimumA, long int maximumA, long int minimumB, long int maximumB, bool refining) const;
		bool checkIntersectionRange(double coordA, double coordB) const { return coordA >= m_localCoordinates.minimumA && coordA <= m_localCoordinates.maximumA && coordB >= m_localCoordinates.minimumB && coordB <= m_localCoordinates.maximumB; }
		bool calculateIntersectionRayC(double coordA, double coordB, double &coordC, double epsilon) const;
		
		std::size_t getCellIdForward() const { return m_localCoordinates.cellIdForward; }
		std::size_t getCellIdBackward() const { return m_localCoordinates.cellIdBackward; }
		
	private:
		struct LocalCoordinates {
			double minimumA = std::numeric_limits<double>::max();
			double maximumA = std::numeric_limits<double>::lowest();
			double minimumB = std::numeric_limits<double>::max();
			double maximumB = std::numeric_limits<double>::lowest();
			Vector2 localCentroid{0.0, 0.0};
			Vector3 normal{0.0, 0.0, 0.0};
			double normalCentroid = 0.0;
			std::vector<Vector2> edgeVectors;
			std::vector<double> triangleAreas;
			std::size_t cellIdForward = std::numeric_limits<std::size_t>::max();
			std::size_t cellIdBackward = std::numeric_limits<std::size_t>::max();
		};
		
		VoronoiTessellation *m_voronoiTessellation;
		std::size_t m_id;
		std::vector<std::size_t> m_vertexIds;
		std::size_t m_cellIdFirst = 0;
		std::size_t m_cellIdSecond = 0;
		bool m_directNeighbor = false;
		double m_area = 0.0;
		LocalCoordinates m_localCoordinates;
		
		Face(VoronoiTessellation *voronoiTessellation, std::size_t id) : m_voronoiTessellation(voronoiTessellation), m_id(id) {}
	};
	
	class Cell
	{
	public:
		friend class CalcVoronoi;
		friend class CalcVoronoiIntegrals;
		friend class VoronoiTessellation;
		
		struct Integrals {
			double charge = 0.0;
			Vector3 electricDipole{0.0, 0.0, 0.0};
			MatrixS3 electricQuadrupole{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
			Vector3 totalCurrent{0.0, 0.0, 0.0};
			Vector3 magneticDipole{0.0, 0.0, 0.0};
		};
		
		std::size_t getId() const { return m_id; }
		std::size_t getNumFaces() const { return m_faceIds.size(); }
		std::size_t getFaceId(std::size_t localId) const { return m_faceIds.at(localId); }
		const Face &getFace(std::size_t localId) const { return m_voronoiTessellation->m_faces[getFaceId(localId)]; }
		const Vector3 &getCentroid() const { return m_centroid; }
		const Vector3 &getAtomPosition() const { return m_atomPosition; }
		double getAtomRadius() const { return m_atomRadius; }
		
		const Integrals &getIntegrals() const { return m_integrals; }
		void multiplyIntegrals(double factor) {
			m_integrals.charge *= factor;
			m_integrals.electricDipole *= factor;
			m_integrals.electricQuadrupole *= factor;
			m_integrals.totalCurrent *= factor;
			m_integrals.magneticDipole *= factor;
		}
		
		std::vector<std::size_t> getFaceIdsByNeighborId(std::size_t neighborId) const;
		
	private:
		VoronoiTessellation *m_voronoiTessellation;
		std::size_t m_id;
		Vector3 m_centroid;
		Vector3 m_atomPosition;
		double m_atomRadius;
		std::vector<std::size_t> m_faceIds;
		double m_volume;
		Integrals m_integrals;
		
		Cell(VoronoiTessellation *voronoiTessellation, std::size_t id) : m_voronoiTessellation(voronoiTessellation), m_id(id) {}
	};
	
	VoronoiTessellation() = default;
	VoronoiTessellation(const VoronoiTessellation &other) : m_vertices(other.m_vertices), m_faces(other.m_faces), m_cells(other.m_cells) {
		p_updatePointers();
	}
	VoronoiTessellation(VoronoiTessellation &&other) noexcept {
		swap(*this, other);
	}
	VoronoiTessellation &operator=(VoronoiTessellation other) noexcept {
		swap(*this, other);
		return *this;
	}
	friend void swap(VoronoiTessellation &first, VoronoiTessellation &second) noexcept {
		using std::swap;
		swap(first.m_vertices, second.m_vertices);
		swap(first.m_faces, second.m_faces);
		swap(first.m_cells, second.m_cells);
		first.p_updatePointers();
		second.p_updatePointers();
	}
	
	const std::vector<Vector3> &getVertices() const { return m_vertices; }
	const Vector3 &getVertex(std::size_t id) const { return m_vertices.at(id); }
	const std::vector<Face> &getFaces() const { return m_faces; }
	const Face &getFace(std::size_t id) const { return m_faces.at(id); }
	const std::vector<Cell> &getCells() const { return m_cells; }
	const Cell &getCell(std::size_t id) const { return m_cells.at(id); }
	
	void createFaceLocalCoordinates(const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, double epsilon) {
		for (auto &&f: m_faces)
			f.createLocalCoordinates(cellInfo, gridSizes, epsilon);
	}
	void multiplyCellIntegrals(double factor) {
		for (auto &&c: m_cells)
			c.multiplyIntegrals(factor);
	}
	
private:
	std::vector<Vector3> m_vertices;
	std::vector<Face> m_faces;
	std::vector<Cell> m_cells;
	
	void p_updatePointers() noexcept {
		for (auto &&f: m_faces)
			f.m_voronoiTessellation = this;
		for (auto &&c: m_cells)
			c.m_voronoiTessellation = this;
	}
};

}

}

#endif

#endif


