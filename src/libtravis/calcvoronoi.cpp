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


#include "calcvoronoi.h"

#include "cellinfo.h"
#include "voronoitessellation.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#ifdef __clang__
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wextra-semi"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#pragma clang diagnostic ignored "-Wshadow-field"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#endif

#include "voro++/voro++.hh"

#pragma GCC diagnostic pop

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace Libtravis {

namespace Travis {

VoronoiTessellation CalcVoronoi::calcVoronoi(const std::vector<Vector3> &atomPositions, const std::vector<double> &atomRadii, const CellInfo &cellInfo) const {
	if (atomPositions.size() >= std::numeric_limits<int>::max())
		throw std::runtime_error("CalcVoronoi: cannot handle more than " + std::to_string(std::numeric_limits<int>::max()) + " atoms");
	
	VoronoiTessellation tessellation;
	
	if (cellInfo.isFullyPeriodic()) {
		if (cellInfo.getCellMatrix()[3] != 0.0 || cellInfo.getCellMatrix()[6] != 0.0)
			throw std::runtime_error("CalcVoronoi: first cell vector has to be aligned with x axis");
		if (cellInfo.getCellMatrix()[7] != 0.0)
			throw std::runtime_error("CalcVoronoi: second cell vector has to be aligned with xy plane");
		
		int nblock = static_cast<int>(std::cbrt(atomPositions.size() / 5));
		if (nblock < 1)
			nblock = 1;
//		printf("$$ %f %f %f %f %f %f %f %f %f\n",cellInfo.getCellMatrix()[0],cellInfo.getCellMatrix()[1],cellInfo.getCellMatrix()[2],cellInfo.getCellMatrix()[3],cellInfo.getCellMatrix()[4],cellInfo.getCellMatrix()[5],cellInfo.getCellMatrix()[6],cellInfo.getCellMatrix()[7],cellInfo.getCellMatrix()[8]);
		Libtravis::voro::container_periodic_poly container(cellInfo.getCellMatrix()[0], cellInfo.getCellMatrix()[1], cellInfo.getCellMatrix()[4], cellInfo.getCellMatrix()[2], cellInfo.getCellMatrix()[5], cellInfo.getCellMatrix()[8], nblock, nblock, nblock, 8);
		Libtravis::voro::particle_order order(static_cast<int>(atomPositions.size()));
		for (std::size_t i = 0; i < atomPositions.size(); ++i) {
//			if (i < 10)
//					printf("$$ %f %f %f %f\n",atomPositions[i][0], atomPositions[i][1], atomPositions[i][2], atomRadii.at(i));
			container.put(order, static_cast<int>(i), atomPositions[i][0], atomPositions[i][1], atomPositions[i][2], atomRadii.at(i));
		}
		
		Libtravis::voro::c_loop_order_periodic cloop(container, order);
		Libtravis::voro::voronoicell_neighbor cell;
		std::vector<int> neighbors;
		std::vector<int> faceVertices;
		std::vector<int> vertexOrders;
		std::vector<double> vertices;
		Vector3 position;
		Vector3 centroid;
		double radius;
		int pid = -1;
		if (!cloop.start())
			throw std::runtime_error("CalcVoronoi: no elements in loop");
		do {
			if (!container.compute_cell(cell, cloop)) {
				++pid;
				tessellation.m_cells.emplace_back(VoronoiTessellation::Cell(&tessellation, static_cast<std::size_t>(pid)));
				tessellation.m_cells.back().m_centroid = atomPositions.at(pid);
				tessellation.m_cells.back().m_atomPosition = atomPositions.at(pid);
				tessellation.m_cells.back().m_atomRadius = atomRadii.at(pid);
				continue;
			}
			double px, py, pz, pr, cx, cy, cz;
			cloop.pos(pid, px, py, pz, pr);
			position[0] = px;
			position[1] = py;
			position[2] = pz;
			radius = pr;
			cell.vertices(px, py, pz, vertices);
			cell.vertex_orders(vertexOrders);
			cell.face_vertices(faceVertices);
			cell.neighbors(neighbors);
			cell.centroid(cx, cy, cz);
//			if (pid < 10)
//					printf("$$ %f\n",cell.volume());
			centroid[0] = cx + px;
			centroid[1] = cy + py;
			centroid[2] = cz + pz;
			
			std::vector<std::size_t> vertexIds(vertices.size() / 3, std::numeric_limits<std::size_t>::max());
			std::vector<std::size_t> neighborFaceIds(neighbors.size(), std::numeric_limits<std::size_t>::max());
			std::size_t vertIndex = 0;
			for (std::size_t i = 0; i < neighbors.size(); ++i) {
				int numVertices = faceVertices[vertIndex];
				++vertIndex;
				// Check if neighbor is known
				if (neighbors[i] < pid) {
					// Get face candidates
					std::vector<std::size_t> faces = tessellation.m_cells[static_cast<std::size_t>(neighbors[i])].getFaceIdsByNeighborId(static_cast<std::size_t>(pid));
					std::size_t vertIndexStart = vertIndex;
					for (std::size_t j = 0; j < faces.size(); ++j) {
						const VoronoiTessellation::Face &face = tessellation.m_faces[faces[j]];
						// Cannot be neighbor if number of vertices is different
						if (face.getNumVertices() != static_cast<std::size_t>(numVertices))
							continue;
						vertIndex = vertIndexStart;
						neighborFaceIds[i] = faces[j];
						for (int k = 0; k < numVertices; ++k) {
							// Check if vertex is known
							if (vertexIds[static_cast<std::size_t>(faceVertices[vertIndex])] == std::numeric_limits<std::size_t>::max()) {
								const Vector3 vertex(vertices[3 * static_cast<std::size_t>(faceVertices[vertIndex])], vertices[3 * static_cast<std::size_t>(faceVertices[vertIndex]) + 1], vertices[3 * static_cast<std::size_t>(faceVertices[vertIndex]) + 2]);
								std::size_t index = std::numeric_limits<std::size_t>::max();
								double dist = std::numeric_limits<double>::max();
								for (std::size_t l = 0; l < face.getNumVertices(); ++l) {
									Vector3 vertex1 = face.getVertex(l);
									vertex1 -= vertex;
									double dist1 = vertex1.norm();
									if (dist1 < dist) {
										index = l;
										dist = dist1;
									}
								}
								if (dist < m_vertexIdentityThreshold) {
									vertexIds[static_cast<std::size_t>(faceVertices[vertIndex])] = face.getVertexId(index);
								} else {
									neighborFaceIds[i] = std::numeric_limits<std::size_t>::max();
								}
							} else {
								for (std::size_t l = 0; l < face.getNumVertices(); ++l) {
									if (neighborFaceIds[i] == std::numeric_limits<std::size_t>::max())
										break;
									if (std::find(face.m_vertexIds.cbegin(), face.m_vertexIds.cend(), vertexIds[static_cast<std::size_t>(faceVertices[vertIndex])]) == face.m_vertexIds.cend()) {
										neighborFaceIds[i] = std::numeric_limits<std::size_t>::max();
									}
								}
							}
							++vertIndex;
						}
						if (neighborFaceIds[i] != std::numeric_limits<std::size_t>::max())
							break;
					}
				} else {
					vertIndex += static_cast<std::size_t>(numVertices);
				}
			}
			
			for (std::size_t i = 0; i < vertices.size() / 3; ++i) {
				if (vertexIds[i] == std::numeric_limits<std::size_t>::max()) {
					vertexIds[i] = tessellation.m_vertices.size();
					tessellation.m_vertices.emplace_back(vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
				}
			}
			
			vertIndex = 0;
			tessellation.m_cells.emplace_back(VoronoiTessellation::Cell(&tessellation, static_cast<std::size_t>(pid)));
			tessellation.m_cells.back().m_centroid = centroid;
			tessellation.m_cells.back().m_atomPosition = position;
			tessellation.m_cells.back().m_atomRadius = radius;
			tessellation.m_cells.back().m_faceIds.resize(neighbors.size());
			for (std::size_t i = 0; i < neighbors.size(); ++i) {
				int numVertices = faceVertices[vertIndex];
				++vertIndex;
				if (neighborFaceIds[i] != std::numeric_limits<std::size_t>::max()) {
					vertIndex += static_cast<std::size_t>(numVertices);
					tessellation.m_cells.back().m_faceIds[i] = neighborFaceIds[i];
					tessellation.m_faces[neighborFaceIds[i]].m_directNeighbor = true;
				} else {
					tessellation.m_faces.emplace_back(VoronoiTessellation::Face(&tessellation, tessellation.m_faces.size()));
					tessellation.m_faces.back().m_vertexIds.resize(static_cast<std::size_t>(numVertices));
					for (std::size_t j = 0; j < static_cast<std::size_t>(numVertices); ++j) {
						tessellation.m_faces.back().m_vertexIds[j] = vertexIds[static_cast<std::size_t>(faceVertices[vertIndex])];
						++vertIndex;
					}
					tessellation.m_faces.back().m_cellIdFirst = static_cast<std::size_t>(pid);
					tessellation.m_faces.back().m_cellIdSecond = static_cast<std::size_t>(neighbors[i]);
					tessellation.m_cells.back().m_faceIds[i] = tessellation.m_faces.size() - 1;
				}
			}
		} while (cloop.inc());
	} else {
		throw std::runtime_error("CalcVoronoi: Voronoi tessellation for nonperiodic cells is not implemented");
	}
	
	return tessellation;
}

}

}

#endif


