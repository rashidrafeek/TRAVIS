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


#include "calcvoronoiintegrals.h"

#include "cellinfo.h"
#include "datagrid.h"
#include "vector.h"
#include "voronoitessellation.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace Libtravis {

namespace Travis {

bool CalcVoronoiIntegrals::IntegrationMask::RaySegment::isInRange(std::size_t indexC, unsigned int refIndexC) const {
	if (indexC == m_segmentC.first) {
		if (indexC == m_segmentC.second)
			return refIndexC >= m_refSegmentC.first && refIndexC <= m_refSegmentC.second;
		return refIndexC >= m_refSegmentC.first;
	}
	if (indexC == m_segmentC.second)
		return refIndexC <= m_refSegmentC.second;
	return indexC > m_segmentC.first && indexC < m_segmentC.second;
}

bool CalcVoronoiIntegrals::IntegrationMask::RaySegmentation::isInRange(std::size_t indexC) const {
	for (auto &&s: m_segments)
		if (s.get().isInRange(indexC))
			return true;
	return false;
}

bool CalcVoronoiIntegrals::IntegrationMask::RaySegmentation::isInRange(std::size_t indexC, unsigned int refIndexC) const {
	for (auto &&s: m_segments)
		if (s.get().isInRange(indexC, refIndexC))
			return true;
	return false;
}

CalcVoronoiIntegrals::IntegrationMask::RaySegmentation CalcVoronoiIntegrals::IntegrationMask::getSegmentation(std::size_t indexA, std::size_t indexB) const {
	RaySegmentation segmentation;
	for (auto &&s: m_raySegments)
		if (s.m_indexA == indexA && s.m_indexB == indexB)
			segmentation.m_segments.emplace_back(s);
	return segmentation;
}

CalcVoronoiIntegrals::IntegrationMask::RaySegmentation CalcVoronoiIntegrals::IntegrationMask::getSegmentation(std::size_t indexA, unsigned int refIndexA, std::size_t indexB, unsigned int refIndexB) const {
	RaySegmentation segmentation;
	for (auto &&s: m_raySegments)
		if (s.m_indexA == indexA && s.m_refIndexA == refIndexA && s.m_indexB == indexB && s.m_refIndexB == refIndexB)
			segmentation.m_segments.emplace_back(s);
	return segmentation;
}

void CalcVoronoiIntegrals::IntegrationMask::addIndexLast(std::size_t indexA, std::size_t indexB, std::size_t indexC) {
	if (!m_raySegments.empty() && m_raySegments.back().m_indexA == indexA && m_raySegments.back().m_indexB == indexB && m_raySegments.back().m_segmentC.second + 1 == indexC)
		++m_raySegments.back().m_segmentC.second;
	else
		m_raySegments.emplace_back(indexA, indexB, indexC);
}

void CalcVoronoiIntegrals::IntegrationMask::addIndexLast(std::size_t indexA, unsigned int refIndexA, std::size_t indexB, unsigned int refIndexB, std::size_t indexC, unsigned int refIndexC) {
	if (!m_raySegments.empty() && m_raySegments.back().m_indexA == indexA && m_raySegments.back().m_refIndexA == refIndexA && m_raySegments.back().m_indexB == indexB && m_raySegments.back().m_refIndexB == refIndexB) {
		if (refIndexC == 0) {
			if (m_raySegments.back().m_segmentC.second + 1 == indexC && m_raySegments.back().m_refSegmentC.second + 1 == m_refining[2]) {
				++m_raySegments.back().m_segmentC.second;
				m_raySegments.back().m_refSegmentC.second = 0;
			} else {
				m_raySegments.emplace_back(indexA, refIndexA, indexB, refIndexB, indexC, refIndexC);
			}
		} else {
			if (m_raySegments.back().m_segmentC.second == indexC && m_raySegments.back().m_refSegmentC.second + 1 == refIndexC) {
				++m_raySegments.back().m_refSegmentC.second;
			} else {
				m_raySegments.emplace_back(indexA, refIndexA, indexB, refIndexB, indexC, refIndexC);
			}
		}
	} else {
		m_raySegments.emplace_back(indexA, refIndexA, indexB, refIndexB, indexC, refIndexC);
	}
}

void CalcVoronoiIntegrals::IntegrationConfig::setCalcElectricMoments(const DataGrid3D<double> &electronDensity, bool calcCharge, bool calcDipole, bool calcQuadrupole) {
	m_electronDensity = &electronDensity;
	m_flags[Charge] = calcCharge;
	m_flags[ElectricDipole] = calcDipole;
	m_flags[ElectricQuadrupole] = calcQuadrupole;
}

void CalcVoronoiIntegrals::IntegrationConfig::setCalcMagneticMoments(const DataGrid3D<Vector3> &electricCurrentDensity, bool calcCurrent, bool calcDipole) {
	m_electricCurrentDensity = &electricCurrentDensity;
	m_flags[Current] = calcCurrent;
	m_flags[MagneticDipole] = calcDipole;
}

void CalcVoronoiIntegrals::IntegrationConfig::setCalcMaskSet(std::vector<IntegrationMask> &integrationMasks, const std::vector<std::vector<std::size_t>> &maskCellIds) {
	m_integrationMasks = &integrationMasks;
	m_maskCellIds = &maskCellIds;
	m_flags[Mask] = true;
}

void CalcVoronoiIntegrals::IntegrationConfig::setCalcSanityCheck(DataGrid3D<unsigned int> &sanityData) {
	m_sanityData = &sanityData;
	m_flags[Sanity] = true;
}

void CalcVoronoiIntegrals::IntegrationConfig::setCalcSanityCheckRefined(DataGrid3D<double> &sanityData) {
	m_sanityDataRefined = &sanityData;
	m_flags[SanityRefined] = true;
}

void CalcVoronoiIntegrals::IntegrationConfig::setIntegrationMask(const IntegrationMask &integrationMask) {
	m_masked = true;
	m_integrationMask = &integrationMask;
}

void CalcVoronoiIntegrals::calcIntegrals(VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, const IntegrationConfig &config) const {
	if (config.m_masked && config.m_refining != config.m_integrationMask->getRefining())
		throw std::runtime_error("CalcVoronoiIntegrals: Mismatch in refining");
	
	double ed;
	Vector3 cd;
	std::array<std::size_t, 3> mi;
	std::vector<std::size_t> cellGroups;
	if (config.m_flags[IntegrationConfig::Mask]) {
		const std::vector<std::vector<std::size_t>> &maskCellIds = *config.m_maskCellIds;
		cellGroups.resize(tessellation.m_cells.size(), maskCellIds.size());
		for (std::size_t i = 0; i < maskCellIds.size(); ++i) {
			for (std::size_t j: maskCellIds[i]) {
				cellGroups.at(j) = i;
			}
		}
		config.m_integrationMasks->clear();
		config.m_integrationMasks->resize(maskCellIds.size(), IntegrationMask(config.m_refining));
	}
	
	std::vector<IntegrandType> integrands;
	if (config.m_flags[IntegrationConfig::Charge])
		integrands.emplace_back([&ed, &tessellation] (std::size_t cellId, const Vector3 &) { tessellation.m_cells[cellId].m_integrals.charge += ed; });
	if (config.m_flags[IntegrationConfig::ElectricDipole])
		integrands.emplace_back([&ed, &tessellation] (std::size_t cellId, const Vector3 &position) {
			const Vector3 pos = position - tessellation.m_cells[cellId].m_atomPosition;
			tessellation.m_cells[cellId].m_integrals.electricDipole += ed * pos;
		});
	if (config.m_flags[IntegrationConfig::ElectricQuadrupole])
		integrands.emplace_back([&ed, &tessellation] (std::size_t cellId, const Vector3 &position) {
			const Vector3 pos = position - tessellation.m_cells[cellId].m_atomPosition;
			const double posn = pos.dot(pos);
			tessellation.m_cells[cellId].m_integrals.electricQuadrupole[0] += ed * (1.5 * pos[0] * pos[0] - 0.5 * posn);
			tessellation.m_cells[cellId].m_integrals.electricQuadrupole[1] += ed * 1.5 * pos[0] * pos[1];
			tessellation.m_cells[cellId].m_integrals.electricQuadrupole[2] += ed * 1.5 * pos[0] * pos[2];
			tessellation.m_cells[cellId].m_integrals.electricQuadrupole[3] += ed * 1.5 * pos[1] * pos[0];
			tessellation.m_cells[cellId].m_integrals.electricQuadrupole[4] += ed * (1.5 * pos[1] * pos[1] - 0.5 * posn);
			tessellation.m_cells[cellId].m_integrals.electricQuadrupole[5] += ed * 1.5 * pos[1] * pos[2];
			tessellation.m_cells[cellId].m_integrals.electricQuadrupole[6] += ed * 1.5 * pos[2] * pos[0];
			tessellation.m_cells[cellId].m_integrals.electricQuadrupole[7] += ed * 1.5 * pos[2] * pos[1];
			tessellation.m_cells[cellId].m_integrals.electricQuadrupole[8] += ed * (1.5 * pos[2] * pos[2] - 0.5 * posn);
		});
	if (config.m_flags[IntegrationConfig::Current])
		integrands.emplace_back([&cd, &tessellation] (std::size_t cellId, const Vector3 &) { tessellation.m_cells[cellId].m_integrals.totalCurrent += cd; });
	if (config.m_flags[IntegrationConfig::MagneticDipole])
		integrands.emplace_back([&cd, &tessellation] (std::size_t cellId, const Vector3 &position) {
			const Vector3 pos = position - tessellation.m_cells[cellId].m_atomPosition;
			tessellation.m_cells[cellId].m_integrals.magneticDipole += 0.5 * pos.cross(cd);
		});
	
	if (config.m_refining[0] == 1 && config.m_refining[1] == 1 && config.m_refining[2] == 1) {
		if (config.m_flags[IntegrationConfig::Mask])
			integrands.emplace_back([&config, &cellGroups, &mi] (std::size_t cellId, const Vector3 &) {
				if (cellGroups[cellId] != std::numeric_limits<std::size_t>::max())
					(*config.m_integrationMasks)[cellGroups[cellId]].addIndexLast(mi[0], mi[1], mi[2]);
			});
		
		std::vector<PreparatorType> preparators;
		if (config.m_flags[IntegrationConfig::Charge] || config.m_flags[IntegrationConfig::ElectricDipole] || config.m_flags[IntegrationConfig::ElectricQuadrupole]) {
			if (config.m_electronDensity->getSizes() != gridSizes)
				throw std::runtime_error("CalcVoronoiIntegrals: Grid size mismatch in electron density");
			preparators.emplace_back([&config, &ed] (std::size_t index, const std::array<std::size_t, 3> &) { ed = (*config.m_electronDensity)[index]; });
		}
		if (config.m_flags[IntegrationConfig::Current] || config.m_flags[IntegrationConfig::MagneticDipole]) {
			if (config.m_electricCurrentDensity->getSizes() != gridSizes)
				throw std::runtime_error("CalcVoronoiIntegrals: Grid size mismatch in current density");
			preparators.emplace_back([&config, &cd] (std::size_t index, const std::array<std::size_t, 3> &) { cd = (*config.m_electricCurrentDensity)[index]; });
		}
		if (config.m_flags[IntegrationConfig::Mask])
			preparators.emplace_back([&mi] (std::size_t, const std::array<std::size_t, 3> &index) { mi = index; });
		if (config.m_flags[IntegrationConfig::Sanity]) {
			if (config.m_sanityData->getSizes() != gridSizes)
				throw std::runtime_error("CalcVoronoiIntegrals: Grid size mismatch in sanity data");
			preparators.emplace_back([&config] (std::size_t index, const::std::array<std::size_t, 3> &) { (*config.m_sanityData)[index] += 1; });
		}
		
		if (config.m_masked)
			p_performIntegrationMasked(tessellation, cellInfo, gridSizes, preparators, integrands, *config.m_integrationMask);
		else
			p_performIntegration(tessellation, cellInfo, gridSizes, preparators, integrands);
	} else {
		std::array<unsigned int, 3> mri;
		if (config.m_flags[IntegrationConfig::Mask])
			integrands.emplace_back([&config, &cellGroups, &mi, &mri] (std::size_t cellId, const Vector3 &) {
				if (cellGroups[cellId] != std::numeric_limits<std::size_t>::max())
					(*config.m_integrationMasks)[cellGroups[cellId]].addIndexLast(mi[0], mri[0], mi[1], mri[1], mi[2], mri[2]);
			});
		
		std::vector<PreparatorRefinedType> preparators;
		if (config.m_flags[IntegrationConfig::Charge] || config.m_flags[IntegrationConfig::ElectricDipole] || config.m_flags[IntegrationConfig::ElectricQuadrupole]) {
			if (config.m_electronDensity->getSizes() != gridSizes)
				throw std::runtime_error("CalcVoronoiIntegrals: Grid size mismatch in electron density");
			preparators.emplace_back([&config, &ed] (const std::array<std::size_t, 8> &indices, const std::array<double, 8> &fractions, unsigned int refiningFactor, const std::array<std::size_t, 3> &, const std::array<unsigned int, 3> &) {
				ed = 0.0;
				for (std::size_t i = 0; i < 8; ++i)
					ed += fractions[i] * (*config.m_electronDensity)[indices[i]];
				ed /= refiningFactor;
			});
		}
		if (config.m_flags[IntegrationConfig::Current] || config.m_flags[IntegrationConfig::MagneticDipole]) {
			if (config.m_electricCurrentDensity->getSizes() != gridSizes)
				throw std::runtime_error("CalcVoronoiIntegrals: Grid size mismatch in current density");
			preparators.emplace_back([&config, &cd] (const std::array<std::size_t, 8> &indices, const std::array<double, 8> &fractions, unsigned int refiningFactor, const std::array<std::size_t, 3> &, const std::array<unsigned int, 3> &) {
				cd = Vector3(0.0, 0.0, 0.0);
				for (std::size_t i = 0; i < 8; ++i)
					cd += fractions[i] * (*config.m_electricCurrentDensity)[indices[i]];
				cd /= refiningFactor;
			});
		}
		if (config.m_flags[IntegrationConfig::Mask])
			preparators.emplace_back([&mi, &mri] (const std::array<std::size_t, 8> &, const std::array<double, 8> &, unsigned int, const std::array<std::size_t, 3> &index, const std::array<unsigned int, 3> &refIndex) { mi = index; mri = refIndex; });
		if (config.m_flags[IntegrationConfig::SanityRefined]) {
			if (config.m_sanityDataRefined->getSizes() != gridSizes)
				throw std::runtime_error("CalcVoronoiIntegrals: Grid size mismatch in sanity data");
			preparators.emplace_back([&config] (const std::array<std::size_t, 8> &indices, const std::array<double, 8> &fractions, unsigned int refiningFactor, const std::array<std::size_t, 3> &, const std::array<unsigned int, 3> &) {
				for (std::size_t i = 0; i < 8; ++i)
					(*config.m_sanityDataRefined)[indices[i]] += fractions[i] / refiningFactor;
			});
		}
		
		if (config.m_masked)
			p_performIntegrationRefinedMasked(tessellation, cellInfo, gridSizes, config.m_refining, preparators, integrands, *config.m_integrationMask);
		else
			p_performIntegrationRefined(tessellation, cellInfo, gridSizes, config.m_refining, preparators, integrands);
	}
	
	const double volumeElement = cellInfo.getVolume() / (gridSizes[0] * gridSizes[1] * gridSizes[2]);
	tessellation.multiplyCellIntegrals(volumeElement);
}

void CalcVoronoiIntegrals::p_calculateNextStep() {
	if (!m_cellInfoAccessor)
		throw std::runtime_error("CalcVoronoiIntegrals: source for cell info not set");
	if (!m_voronoiTessellationAccessor)
		throw std::runtime_error("CalcVoronoiIntegrals: source for Voronoi tessellation not set");

	IntegrationConfig config;
	std::array<std::size_t, 3> gridSizes;
	config.setRefining(m_refining);
	if (m_calculateCharge || m_calculateElectricDipole || m_calculateElectricQuadrupole) {
		if (!m_electronDensityAccessor)
			throw std::runtime_error("CalcVoronoiIntegrals: source for electron density not set");
		gridSizes = m_electronDensityAccessor().getSizes();
		config.setCalcElectricMoments(m_electronDensityAccessor(), m_calculateCharge, m_calculateElectricDipole, m_calculateElectricQuadrupole);
	}
	if (m_calculateCurrent || m_calculateMagneticDipole) {
		if (!m_electricCurrentDensityAccessor)
			throw std::runtime_error("CalcVoronoiIntegrals: source for electric current density not set");
		gridSizes = m_electricCurrentDensityAccessor().getSizes();
		config.setCalcMagneticMoments(m_electricCurrentDensityAccessor(), m_calculateCurrent, m_calculateMagneticDipole);
	}
	if (m_calculateMaskSet) {
		if (!m_maskSetTargetAccessor)
			throw std::runtime_error("CalcVoronoiIntegrals: target for mask set not set");
		if (!m_maskCellIdAccessor)
			throw std::runtime_error("CalcVoronoiIntegrals: source for mask cell IDs not set");
		config.setCalcMaskSet(m_maskSetTargetAccessor(), m_maskCellIdAccessor());
	}
	calcIntegrals(m_voronoiTessellationAccessor(), m_cellInfoAccessor(), gridSizes, config);
}

void CalcVoronoiIntegrals::p_performIntegration(VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, const std::vector<PreparatorType> &preparators, const std::vector<IntegrandType> &integrands) const {
	tessellation.createFaceLocalCoordinates(cellInfo, gridSizes, m_epsilon);
	
	long int minPositionA, maxPositionA, minPositionB, maxPositionB;
	std::tie(minPositionA, maxPositionA, minPositionB, maxPositionB) = p_integrationCalcBoundsAB(tessellation, cellInfo, gridSizes, false);
	
	std::vector<std::vector<std::size_t>> faceCandidates = p_integrationGetFaceCandidates(tessellation, minPositionA, maxPositionA, minPositionB, maxPositionB, false);
	const std::size_t sizeB = static_cast<std::size_t>(maxPositionB - minPositionB + 1);
	
	const std::size_t gridSizes12 = gridSizes[1] * gridSizes[2];
	if (cellInfo.isOrthorhombic()) {
		const std::array<double, 3> incrementalCellVectors = { cellInfo.getVectorA()[0] / gridSizes[0], cellInfo.getVectorB()[1] / gridSizes[1], cellInfo.getVectorC()[2] / gridSizes[2] };
		for (long int i = minPositionA; i <= maxPositionA; ++i) {
			const std::size_t indA = DataGrid3D<double>::foldToSize(i, gridSizes[0]);
			const std::size_t indexA = indA * gridSizes12;
			const std::size_t candA = static_cast<std::size_t>(i - minPositionA) * sizeB;
			const double positionVectorA = i * incrementalCellVectors[0];
			for (long int j = minPositionB; j <= maxPositionB; ++j) {
				const std::vector<std::size_t> &cand = faceCandidates[candA + static_cast<std::size_t>(j - minPositionB)];
				if (cand.empty())
					continue;
				
				const std::size_t indB = DataGrid3D<double>::foldToSize(j, gridSizes[1]);
				const std::size_t indexB = indexA + indB * gridSizes[2];
				const double positionVectorB = j * incrementalCellVectors[1];
				
				const std::vector<std::pair<double, std::size_t>> facesHit = p_integrationCalcIntersectionsRayC(tessellation, static_cast<double>(i), static_cast<double>(j), cand);
				if (facesHit.empty())
					continue;
				
				long int positionC, gridShift;
				std::size_t gridPositionC;
				std::tie(positionC, gridShift, gridPositionC) = p_integrationInitLoopC(cellInfo.isPeriodicC(), gridSizes[2], facesHit.front().first);
				
				auto faceIt = facesHit.cbegin();
				std::size_t cellId = std::numeric_limits<std::size_t>::max();
				p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
				
				while (faceIt != facesHit.cend()) {
					if (cellId != std::numeric_limits<std::size_t>::max()) {
						const std::size_t indexC = indexB + gridPositionC;
						const Vector3 positionVector(positionVectorA, positionVectorB, positionC * incrementalCellVectors[2]);
						for (auto &&preparator: preparators)
							preparator(indexC, {indA, indB, gridPositionC});
						for (auto &&integrand: integrands)
							integrand(cellId, positionVector);
					}
					if (!p_integrationUpdateGridPosition(cellInfo.isPeriodicC(), gridSizes[2], positionC, gridShift, gridPositionC))
						break;
					p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
				}
			}
		}
	} else {
		const std::array<Vector3, 3> incrementalCellVectors = { cellInfo.getVectorA() / gridSizes[0], cellInfo.getVectorB() / gridSizes[1], cellInfo.getVectorC() / gridSizes[2] };
		for (long int i = minPositionA; i <= maxPositionA; ++i) {
			const std::size_t indA = DataGrid3D<double>::foldToSize(i, gridSizes[0]);
			const std::size_t indexA = indA * gridSizes12;
			const std::size_t candA = static_cast<std::size_t>(i - minPositionA) * sizeB;
			const Vector3 positionVectorA = i * incrementalCellVectors[0];
			for (long int j = minPositionB; j <= maxPositionB; ++j) {
				const std::vector<std::size_t> &cand = faceCandidates[candA + static_cast<std::size_t>(j - minPositionB)];
				if (cand.empty())
					continue;
				
				const std::size_t indB = DataGrid3D<double>::foldToSize(j, gridSizes[1]);
				const std::size_t indexB = indexA + indB * gridSizes[2];
				const Vector3 positionVectorB = positionVectorA + j * incrementalCellVectors[1];
				
				const std::vector<std::pair<double, std::size_t>> facesHit = p_integrationCalcIntersectionsRayC(tessellation, static_cast<double>(i), static_cast<double>(j), cand);
				if (facesHit.empty())
					continue;
				
				long int positionC, gridShift;
				std::size_t gridPositionC;
				std::tie(positionC, gridShift, gridPositionC) = p_integrationInitLoopC(cellInfo.isPeriodicC(), gridSizes[2], facesHit.front().first);
				
				auto faceIt = facesHit.cbegin();
				std::size_t cellId = std::numeric_limits<std::size_t>::max();
				p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
				
				while (faceIt != facesHit.cend()) {
					if (cellId != std::numeric_limits<std::size_t>::max()) {
						const std::size_t indexC = indexB + gridPositionC;
						const Vector3 positionVector = positionVectorB + positionC * incrementalCellVectors[2];
						for (auto &&preparator: preparators)
							preparator(indexC, {indA, indB, gridPositionC});
						for (auto &&integrand: integrands)
							integrand(cellId, positionVector);
					}
					if (!p_integrationUpdateGridPosition(cellInfo.isPeriodicC(), gridSizes[2], positionC, gridShift, gridPositionC))
						break;
					p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
				}
			}
		}
	}
}

void CalcVoronoiIntegrals::p_performIntegrationRefined(VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, const std::array<unsigned int, 3> &refining, const std::vector<PreparatorRefinedType> &preparators, const std::vector<IntegrandType> &integrands) const {
	tessellation.createFaceLocalCoordinates(cellInfo, gridSizes, m_epsilon);
	
	long int minPositionA, maxPositionA, minPositionB, maxPositionB;
	std::tie(minPositionA, maxPositionA, minPositionB, maxPositionB) = p_integrationCalcBoundsAB(tessellation, cellInfo, gridSizes, true);
	
	std::vector<std::vector<std::size_t>> faceCandidates = p_integrationGetFaceCandidates(tessellation, minPositionA, maxPositionA, minPositionB, maxPositionB, true);
	const std::size_t sizeB = static_cast<std::size_t>(maxPositionB - minPositionB + 1);
	
	const std::size_t gridSizes12 = gridSizes[1] * gridSizes[2];
	const unsigned int refiningFactor = refining[0] * refining[1] * refining[2];
	if (cellInfo.isOrthorhombic()) {
		const std::array<double, 3> incrementalCellVectors = { cellInfo.getVectorA()[0] / gridSizes[0], cellInfo.getVectorB()[1] / gridSizes[1], cellInfo.getVectorC()[2] / gridSizes[2] };
		for (long int i = minPositionA; i <= maxPositionA; ++i) {
			const std::size_t indA1 = DataGrid3D<double>::foldToSize(i, gridSizes[0]);
			const std::size_t indA2 = DataGrid3D<double>::foldToSize(i + 1, gridSizes[0]);
			const std::size_t indexA1 = indA1 * gridSizes12;
			const std::size_t indexA2 = indA2 * gridSizes12;
			const std::size_t candA = static_cast<std::size_t>(i - minPositionA) * sizeB;
			for (unsigned int j = 0; j < refining[0]; ++j) {
				const double fracA2 = static_cast<double>(j) / refining[0];
				const double fracA1 = 1.0 - fracA2;
				const double positionVectorA = (i + static_cast<double>(j) / refining[0]) * incrementalCellVectors[0];
				for (long int k = minPositionB; k <= maxPositionB; ++k) {
					const std::vector<std::size_t> &cand = faceCandidates[candA + static_cast<std::size_t>(k - minPositionB)];
					if (cand.empty())
						continue;
					
					const std::size_t indB1 = DataGrid3D<double>::foldToSize(k, gridSizes[1]);
					const std::size_t indB2 = DataGrid3D<double>::foldToSize(k + 1, gridSizes[1]);
					const std::size_t indexB1 = indB1 * gridSizes[2];
					const std::size_t indexB2 = indB2 * gridSizes[2];
					const std::array<std::size_t, 4> indexB = { indexA1 + indexB1, indexA1 + indexB2, indexA2 + indexB1, indexA2 + indexB2 };
					for (unsigned int l = 0; l < refining[1]; ++l) {
						const double fracB2 = static_cast<double>(l) / refining[1];
						const double fracB1 = 1.0 - fracB2;
						const std::array<double, 4> fracB = { fracA1 * fracB1, fracA1 * fracB2, fracA2 * fracB1, fracA2 * fracB2 };
						const double positionVectorB = (k + static_cast<double>(l) / refining[1]) * incrementalCellVectors[1];
						
						const std::vector<std::pair<double, std::size_t>> facesHit = p_integrationCalcIntersectionsRayC(tessellation, static_cast<double>(j) / refining[0] + i, static_cast<double>(l) / refining[1] + k, cand);
						if (facesHit.empty())
							continue;
						
						double positionC;
						long int gridShift;
						std::size_t gridPositionC;
						unsigned int refiningPositionC;
						std::tie(positionC, gridShift, gridPositionC, refiningPositionC) = p_integrationRefinedInitLoopC(cellInfo.isPeriodicC(), gridSizes[2], refining[2], facesHit.front().first);
						
						auto faceIt = facesHit.cbegin();
						std::size_t cellId = std::numeric_limits<std::size_t>::max();
						p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
						
						std::size_t indexC2 = (gridPositionC + 1) % gridSizes[2];
						std::array<std::size_t, 8> indexC = { indexB[0] + gridPositionC, indexB[0] + indexC2, indexB[1] + gridPositionC, indexB[1] + indexC2, indexB[2] + gridPositionC, indexB[2] + indexC2, indexB[3] + gridPositionC, indexB[3] + indexC2 };
						while (faceIt != facesHit.cend()) {
							if (cellId != std::numeric_limits<std::size_t>::max()) {
								const double fracC2 = static_cast<double>(refiningPositionC) / refining[2];
								const double fracC1 = 1.0 - fracC2;
								const std::array<double, 8> fracC = {fracB[0] * fracC1, fracB[0] * fracC2, fracB[1] * fracC1, fracB[1] * fracC2, fracB[2] * fracC1, fracB[2] * fracC2, fracB[3] * fracC1, fracB[3] * fracC2};
								const Vector3 positionVector(positionVectorA, positionVectorB, positionC * incrementalCellVectors[2]);
								for (auto &&preparator: preparators)
									preparator(indexC, fracC, refiningFactor, {indA1, indB1, gridPositionC}, {j, l, refiningPositionC});
								for (auto &&integrand: integrands)
									integrand(cellId, positionVector);
							}
							if (!p_integrationRefinedUpdateGridPosition(cellInfo.isPeriodicC(), gridSizes[2], refining[2], indexB, positionC, gridShift, gridPositionC, refiningPositionC, indexC2, indexC))
								break;
							p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
						}
					}
				}
			}
		}
	} else {
		const std::array<Vector3, 3> incrementalCellVectors = { cellInfo.getVectorA() / gridSizes[0], cellInfo.getVectorB() / gridSizes[1], cellInfo.getVectorC() / gridSizes[2] };
		for (long int i = minPositionA; i <= maxPositionA; ++i) {
			const std::size_t indA1 = DataGrid3D<double>::foldToSize(i, gridSizes[0]);
			const std::size_t indA2 = DataGrid3D<double>::foldToSize(i + 1, gridSizes[0]);
			const std::size_t indexA1 = indA1 * gridSizes12;
			const std::size_t indexA2 = indA2 * gridSizes12;
			const std::size_t candA = static_cast<std::size_t>(i - minPositionA) * sizeB;
			for (unsigned int j = 0; j < refining[0]; ++j) {
				const double fracA2 = static_cast<double>(j) / refining[0];
				const double fracA1 = 1.0 - fracA2;
				const Vector3 positionVectorA = (i + static_cast<double>(j) / refining[0]) * incrementalCellVectors[0];
				for (long int k = minPositionB; k <= maxPositionB; ++k) {
					const std::vector<std::size_t> &cand = faceCandidates[candA + static_cast<std::size_t>(k - minPositionB)];
					if (cand.empty())
						continue;
					
					const std::size_t indB1 = DataGrid3D<double>::foldToSize(k, gridSizes[1]);
					const std::size_t indB2 = DataGrid3D<double>::foldToSize(k + 1, gridSizes[1]);
					const std::size_t indexB1 = indB1 * gridSizes[2];
					const std::size_t indexB2 = indB2 * gridSizes[2];
					const std::array<std::size_t, 4> indexB = { indexA1 + indexB1, indexA1 + indexB2, indexA2 + indexB1, indexA2 + indexB2 };
					for (unsigned int l = 0; l < refining[1]; ++l) {
						const double fracB2 = static_cast<double>(l) / refining[1];
						const double fracB1 = 1.0 - fracB2;
						const std::array<double, 4> fracB = { fracA1 * fracB1, fracA1 * fracB2, fracA2 * fracB1, fracA2 * fracB2 };
						const Vector3 positionVectorB = positionVectorA + (k + static_cast<double>(l) / refining[1]) * incrementalCellVectors[1];
						
						const std::vector<std::pair<double, std::size_t>> facesHit = p_integrationCalcIntersectionsRayC(tessellation, static_cast<double>(j) / refining[0] + i, static_cast<double>(l) / refining[1] + k, cand);
						if (facesHit.empty())
							continue;
						
						double positionC;
						long int gridShift;
						std::size_t gridPositionC;
						unsigned int refiningPositionC;
						std::tie(positionC, gridShift, gridPositionC, refiningPositionC) = p_integrationRefinedInitLoopC(cellInfo.isPeriodicC(), gridSizes[2], refining[2], facesHit.front().first);
						
						auto faceIt = facesHit.cbegin();
						std::size_t cellId = std::numeric_limits<std::size_t>::max();
						p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
						
						std::size_t indexC2 = (gridPositionC + 1) % gridSizes[2];
						std::array<std::size_t, 8> indexC = { indexB[0] + gridPositionC, indexB[0] + indexC2, indexB[1] + gridPositionC, indexB[1] + indexC2, indexB[2] + gridPositionC, indexB[2] + indexC2, indexB[3] + gridPositionC, indexB[3] + indexC2 };
						while (faceIt != facesHit.cend()) {
							if (cellId != std::numeric_limits<std::size_t>::max()) {
								const double fracC2 = static_cast<double>(refiningPositionC) / refining[2];
								const double fracC1 = 1.0 - fracC2;
								const std::array<double, 8> fracC = {fracB[0] * fracC1, fracB[0] * fracC2, fracB[1] * fracC1, fracB[1] * fracC2, fracB[2] * fracC1, fracB[2] * fracC2, fracB[3] * fracC1, fracB[3] * fracC2};
								const Vector3 positionVector = positionVectorB + positionC * incrementalCellVectors[2];
								for (auto &&preparator: preparators)
									preparator(indexC, fracC, refiningFactor, {indA1, indB1, gridPositionC}, {j, l, refiningPositionC});
								for (auto &&integrand: integrands)
									integrand(cellId, positionVector);
							}
							if (!p_integrationRefinedUpdateGridPosition(cellInfo.isPeriodicC(), gridSizes[2], refining[2], indexB, positionC, gridShift, gridPositionC, refiningPositionC, indexC2, indexC))
								break;
							p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
						}
					}
				}
			}
		}
	}
}

void CalcVoronoiIntegrals::p_performIntegrationMasked(VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, const std::vector<PreparatorType> &preparators, const std::vector<IntegrandType> &integrands, const IntegrationMask &integrationMask) const {
	tessellation.createFaceLocalCoordinates(cellInfo, gridSizes, m_epsilon);
	
	long int minPositionA, maxPositionA, minPositionB, maxPositionB;
	std::tie(minPositionA, maxPositionA, minPositionB, maxPositionB) = p_integrationCalcBoundsAB(tessellation, cellInfo, gridSizes, false);
	
	std::vector<std::vector<std::size_t>> faceCandidates = p_integrationGetFaceCandidates(tessellation, minPositionA, maxPositionA, minPositionB, maxPositionB, false);
	const std::size_t sizeB = static_cast<std::size_t>(maxPositionB - minPositionB + 1);
	
	const std::size_t gridSizes12 = gridSizes[1] * gridSizes[2];
	if (cellInfo.isOrthorhombic()) {
		const std::array<double, 3> incrementalCellVectors = { cellInfo.getVectorA()[0] / gridSizes[0], cellInfo.getVectorB()[1] / gridSizes[1], cellInfo.getVectorC()[2] / gridSizes[2] };
		for (long int i = minPositionA; i <= maxPositionA; ++i) {
			const std::size_t indA = DataGrid3D<double>::foldToSize(i, gridSizes[0]);
			const std::size_t indexA = indA * gridSizes12;
			const std::size_t candA = static_cast<std::size_t>(i - minPositionA) * sizeB;
			const double positionVectorA = i * incrementalCellVectors[0];
			for (long int j = minPositionB; j <= maxPositionB; ++j) {
				const std::vector<std::size_t> &cand = faceCandidates[candA + static_cast<std::size_t>(j - minPositionB)];
				if (cand.empty())
					continue;
				
				const std::size_t indB = DataGrid3D<double>::foldToSize(j, gridSizes[1]);
				const IntegrationMask::RaySegmentation &integrationRanges = integrationMask.getSegmentation(indA, indB);
				if (integrationRanges.isEmpty())
					continue;
				
				const std::size_t indexB = indexA + indB * gridSizes[2];
				const double positionVectorB = j * incrementalCellVectors[1];
				const std::vector<std::pair<double, std::size_t>> facesHit = p_integrationCalcIntersectionsRayC(tessellation, static_cast<double>(i), static_cast<double>(j), cand);
				if (facesHit.empty())
					continue;
				
				long int positionC, gridShift;
				std::size_t gridPositionC;
				std::tie(positionC, gridShift, gridPositionC) = p_integrationInitLoopC(cellInfo.isPeriodicC(), gridSizes[2], facesHit.front().first);
				
				auto faceIt = facesHit.cbegin();
				std::size_t cellId = std::numeric_limits<std::size_t>::max();
				p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
				
				while (faceIt != facesHit.cend()) {
					if (cellId != std::numeric_limits<std::size_t>::max() && integrationRanges.isInRange(gridPositionC)) {
						const std::size_t indexC = indexB + gridPositionC;
						const Vector3 positionVector(positionVectorA, positionVectorB, positionC * incrementalCellVectors[2]);
						for (auto &&preparator: preparators)
							preparator(indexC, {indA, indB, gridPositionC});
						for (auto &&integrand: integrands)
							integrand(cellId, positionVector);
					}
					if (!p_integrationUpdateGridPosition(cellInfo.isPeriodicC(), gridSizes[2], positionC, gridShift, gridPositionC))
						break;
					p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
				}
			}
		}
	} else {
		const std::array<Vector3, 3> incrementalCellVectors = { cellInfo.getVectorA() / gridSizes[0], cellInfo.getVectorB() / gridSizes[1], cellInfo.getVectorC() / gridSizes[2] };
		for (long int i = minPositionA; i <= maxPositionA; ++i) {
			const std::size_t indA = DataGrid3D<double>::foldToSize(i, gridSizes[0]);
			const std::size_t indexA = indA * gridSizes12;
			const std::size_t candA = static_cast<std::size_t>(i - minPositionA) * sizeB;
			const Vector3 positionVectorA = i * incrementalCellVectors[0];
			for (long int j = minPositionB; j <= maxPositionB; ++j) {
				const std::vector<std::size_t> &cand = faceCandidates[candA + static_cast<std::size_t>(j - minPositionB)];
				if (cand.empty())
					continue;
				
				const std::size_t indB = DataGrid3D<double>::foldToSize(j, gridSizes[1]);
				const IntegrationMask::RaySegmentation &integrationRanges = integrationMask.getSegmentation(indA, indB);
				if (integrationRanges.isEmpty())
					continue;
				
				const std::size_t indexB = indexA + indB * gridSizes[2];
				const Vector3 positionVectorB = positionVectorA + j * incrementalCellVectors[1];
				const std::vector<std::pair<double, std::size_t>> facesHit = p_integrationCalcIntersectionsRayC(tessellation, static_cast<double>(i), static_cast<double>(j), cand);
				if (facesHit.empty())
					continue;
				
				long int positionC, gridShift;
				std::size_t gridPositionC;
				std::tie(positionC, gridShift, gridPositionC) = p_integrationInitLoopC(cellInfo.isPeriodicC(), gridSizes[2], facesHit.front().first);
				
				auto faceIt = facesHit.cbegin();
				std::size_t cellId = std::numeric_limits<std::size_t>::max();
				p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
				
				while (faceIt != facesHit.cend()) {
					if (cellId != std::numeric_limits<std::size_t>::max() && integrationRanges.isInRange(gridPositionC)) {
						const std::size_t indexC = indexB + gridPositionC;
						const Vector3 positionVector = positionVectorB + positionC * incrementalCellVectors[2];
						for (auto &&preparator: preparators)
							preparator(indexC, {indA, indB, gridPositionC});
						for (auto &&integrand: integrands)
							integrand(cellId, positionVector);
					}
					if (!p_integrationUpdateGridPosition(cellInfo.isPeriodicC(), gridSizes[2], positionC, gridShift, gridPositionC))
						break;
					p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
				}
			}
		}
	}
}

void CalcVoronoiIntegrals::p_performIntegrationRefinedMasked(VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, const std::array<unsigned int, 3> &refining, const std::vector<PreparatorRefinedType> &preparators, const std::vector<IntegrandType> &integrands, const IntegrationMask &integrationMask) const {
	tessellation.createFaceLocalCoordinates(cellInfo, gridSizes, m_epsilon);
	
	long int minPositionA, maxPositionA, minPositionB, maxPositionB;
	std::tie(minPositionA, maxPositionA, minPositionB, maxPositionB) = p_integrationCalcBoundsAB(tessellation, cellInfo, gridSizes, true);
	
	std::vector<std::vector<std::size_t>> faceCandidates = p_integrationGetFaceCandidates(tessellation, minPositionA, maxPositionA, minPositionB, maxPositionB, true);
	const std::size_t sizeB = static_cast<std::size_t>(maxPositionB - minPositionB + 1);
	
	const std::size_t gridSizes12 = gridSizes[1] * gridSizes[2];
	const unsigned int refiningFactor = refining[0] * refining[1] * refining[2];
	if (cellInfo.isOrthorhombic()) {
		const std::array<double, 3> incrementalCellVectors = { cellInfo.getVectorA()[0] / gridSizes[0], cellInfo.getVectorB()[1] / gridSizes[1], cellInfo.getVectorC()[2] / gridSizes[2] };
		for (long int i = minPositionA; i <= maxPositionA; ++i) {
			const std::size_t indA1 = DataGrid3D<double>::foldToSize(i, gridSizes[0]);
			const std::size_t indA2 = DataGrid3D<double>::foldToSize(i + 1, gridSizes[0]);
			const std::size_t indexA1 = indA1 * gridSizes12;
			const std::size_t indexA2 = indA2 * gridSizes12;
			const std::size_t candA = static_cast<std::size_t>(i - minPositionA) * sizeB;
			for (unsigned int j = 0; j < refining[0]; ++j) {
				const double fracA2 = static_cast<double>(j) / refining[0];
				const double fracA1 = 1.0 - fracA2;
				const double positionVectorA = (i + static_cast<double>(j) / refining[0]) * incrementalCellVectors[0];
				for (long int k = minPositionB; k <= maxPositionB; ++k) {
					const std::vector<std::size_t> &cand = faceCandidates[candA + static_cast<std::size_t>(k - minPositionB)];
					if (cand.empty())
						continue;
					
					const std::size_t indB1 = DataGrid3D<double>::foldToSize(k, gridSizes[1]);
					const std::size_t indB2 = DataGrid3D<double>::foldToSize(k + 1, gridSizes[1]);
					const std::size_t indexB1 = indB1 * gridSizes[2];
					const std::size_t indexB2 = indB2 * gridSizes[2];
					const std::array<std::size_t, 4> indexB = { indexA1 + indexB1, indexA1 + indexB2, indexA2 + indexB1, indexA2 + indexB2 };
					for (unsigned int l = 0; l < refining[1]; ++l) {
						const IntegrationMask::RaySegmentation &integrationRanges = integrationMask.getSegmentation(indA1, j, indB1, l);
						if (integrationRanges.isEmpty())
							continue;
						
						const double fracB2 = static_cast<double>(l) / refining[1];
						const double fracB1 = 1.0 - fracB2;
						const std::array<double, 4> fracB = { fracA1 * fracB1, fracA1 * fracB2, fracA2 * fracB1, fracA2 * fracB2 };
						const double positionVectorB = (k + static_cast<double>(l) / refining[1]) * incrementalCellVectors[1];
						
						const std::vector<std::pair<double, std::size_t>> facesHit = p_integrationCalcIntersectionsRayC(tessellation, static_cast<double>(j) / refining[0] + i, static_cast<double>(l) / refining[1] + k, cand);
						if (facesHit.empty())
							continue;
						
						double positionC;
						long int gridShift;
						std::size_t gridPositionC;
						unsigned int refiningPositionC;
						std::tie(positionC, gridShift, gridPositionC, refiningPositionC) = p_integrationRefinedInitLoopC(cellInfo.isPeriodicC(), gridSizes[2], refining[2], facesHit.front().first);
						
						auto faceIt = facesHit.cbegin();
						std::size_t cellId = std::numeric_limits<std::size_t>::max();
						p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
						
						std::size_t indexC2 = (gridPositionC + 1) % gridSizes[2];
						std::array<std::size_t, 8> indexC = { indexB[0] + gridPositionC, indexB[0] + indexC2, indexB[1] + gridPositionC, indexB[1] + indexC2, indexB[2] + gridPositionC, indexB[2] + indexC2, indexB[3] + gridPositionC, indexB[3] + indexC2 };
						while (faceIt != facesHit.cend()) {
							if (cellId != std::numeric_limits<std::size_t>::max() && integrationRanges.isInRange(gridPositionC, refiningPositionC)) {
								const double fracC2 = static_cast<double>(refiningPositionC) / refining[2];
								const double fracC1 = 1.0 - fracC2;
								const std::array<double, 8> fracC = {fracB[0] * fracC1, fracB[0] * fracC2, fracB[1] * fracC1, fracB[1] * fracC2, fracB[2] * fracC1, fracB[2] * fracC2, fracB[3] * fracC1, fracB[3] * fracC2};
								const Vector3 positionVector(positionVectorA, positionVectorB, positionC * incrementalCellVectors[2]);
								for (auto &&preparator: preparators)
									preparator(indexC, fracC, refiningFactor, {indA1, indB1, gridPositionC}, {j, l, refiningPositionC});
								for (auto &&integrand: integrands)
									integrand(cellId, positionVector);
							}
							if (!p_integrationRefinedUpdateGridPosition(cellInfo.isPeriodicC(), gridSizes[2], refining[2], indexB, positionC, gridShift, gridPositionC, refiningPositionC, indexC2, indexC))
								break;
							p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
						}
					}
				}
			}
		}
	} else {
		const std::array<Vector3, 3> incrementalCellVectors = { cellInfo.getVectorA() / gridSizes[0], cellInfo.getVectorB() / gridSizes[1], cellInfo.getVectorC() / gridSizes[2] };
		for (long int i = minPositionA; i <= maxPositionA; ++i) {
			const std::size_t indA1 = DataGrid3D<double>::foldToSize(i, gridSizes[0]);
			const std::size_t indA2 = DataGrid3D<double>::foldToSize(i + 1, gridSizes[0]);
			const std::size_t indexA1 = indA1 * gridSizes12;
			const std::size_t indexA2 = indA2 * gridSizes12;
			const std::size_t candA = static_cast<std::size_t>(i - minPositionA) * sizeB;
			for (unsigned int j = 0; j < refining[0]; ++j) {
				const double fracA2 = static_cast<double>(j) / refining[0];
				const double fracA1 = 1.0 - fracA2;
				const Vector3 positionVectorA = (i + static_cast<double>(j) / refining[0]) * incrementalCellVectors[0];
				for (long int k = minPositionB; k <= maxPositionB; ++k) {
					const std::vector<std::size_t> &cand = faceCandidates[candA + static_cast<std::size_t>(k - minPositionB)];
					if (cand.empty())
						continue;
					
					const std::size_t indB1 = DataGrid3D<double>::foldToSize(k, gridSizes[1]);
					const std::size_t indB2 = DataGrid3D<double>::foldToSize(k + 1, gridSizes[1]);
					const std::size_t indexB1 = indB1 * gridSizes[2];
					const std::size_t indexB2 = indB2 * gridSizes[2];
					const std::array<std::size_t, 4> indexB = { indexA1 + indexB1, indexA1 + indexB2, indexA2 + indexB1, indexA2 + indexB2 };
					for (unsigned int l = 0; l < refining[1]; ++l) {
						const IntegrationMask::RaySegmentation &integrationRanges = integrationMask.getSegmentation(indA1, j, indB1, l);
						if (integrationRanges.isEmpty())
							continue;
						
						const double fracB2 = static_cast<double>(l) / refining[1];
						const double fracB1 = 1.0 - fracB2;
						const std::array<double, 4> fracB = { fracA1 * fracB1, fracA1 * fracB2, fracA2 * fracB1, fracA2 * fracB2 };
						const Vector3 positionVectorB = positionVectorA + (k + static_cast<double>(l) / refining[1]) * incrementalCellVectors[1];
						
						const std::vector<std::pair<double, std::size_t>> facesHit = p_integrationCalcIntersectionsRayC(tessellation, static_cast<double>(j) / refining[0] + i, static_cast<double>(l) / refining[1] + k, cand);
						if (facesHit.empty())
							continue;
						
						double positionC;
						long int gridShift;
						std::size_t gridPositionC;
						unsigned int refiningPositionC;
						std::tie(positionC, gridShift, gridPositionC, refiningPositionC) = p_integrationRefinedInitLoopC(cellInfo.isPeriodicC(), gridSizes[2], refining[2], facesHit.front().first);
						
						auto faceIt = facesHit.cbegin();
						std::size_t cellId = std::numeric_limits<std::size_t>::max();
						p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
						
						std::size_t indexC2 = (gridPositionC + 1) % gridSizes[2];
						std::array<std::size_t, 8> indexC = { indexB[0] + gridPositionC, indexB[0] + indexC2, indexB[1] + gridPositionC, indexB[1] + indexC2, indexB[2] + gridPositionC, indexB[2] + indexC2, indexB[3] + gridPositionC, indexB[3] + indexC2 };
						while (faceIt != facesHit.cend()) {
							if (cellId != std::numeric_limits<std::size_t>::max() && integrationRanges.isInRange(gridPositionC, refiningPositionC)) {
								const double fracC2 = static_cast<double>(refiningPositionC) / refining[2];
								const double fracC1 = 1.0 - fracC2;
								const std::array<double, 8> fracC = {fracB[0] * fracC1, fracB[0] * fracC2, fracB[1] * fracC1, fracB[1] * fracC2, fracB[2] * fracC1, fracB[2] * fracC2, fracB[3] * fracC1, fracB[3] * fracC2};
								const Vector3 positionVector = positionVectorB + positionC * incrementalCellVectors[2];
								for (auto &&preparator: preparators)
									preparator(indexC, fracC, refiningFactor, {indA1, indB1, gridPositionC}, {j, l, refiningPositionC});
								for (auto &&integrand: integrands)
									integrand(cellId, positionVector);
							}
							if (!p_integrationRefinedUpdateGridPosition(cellInfo.isPeriodicC(), gridSizes[2], refining[2], indexB, positionC, gridShift, gridPositionC, refiningPositionC, indexC2, indexC))
								break;
							p_integrationUpdateFaceIterator(tessellation, facesHit, positionC, faceIt, cellId);
						}
					}
				}
			}
		}
	}
}

std::tuple<long int, long int, long int, long int> CalcVoronoiIntegrals::p_integrationCalcBoundsAB(const VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, bool refining) const {
	double minA = std::numeric_limits<double>::max();
	double maxA = std::numeric_limits<double>::lowest();
	double minB = std::numeric_limits<double>::max();
	double maxB = std::numeric_limits<double>::lowest();
	for (auto v: tessellation.m_vertices) {
		v = cellInfo.vectorCartesianToFractional(v);
		if (v[0] < minA)
			minA = v[0];
		if (v[0] > maxA)
			maxA = v[0];
		if (v[1] < minB)
			minB = v[1];
		if (v[1] > maxB)
			maxB = v[1];
	}
	minA -= m_epsilon;
	minA *= gridSizes[0];
	maxA += m_epsilon;
	maxA *= gridSizes[0];
	minB -= m_epsilon;
	minB *= gridSizes[1];
	maxB += m_epsilon;
	maxB *= gridSizes[1];
	long int minPositionA = static_cast<long int>(refining ? std::floor(minA) : std::ceil(minA));
	long int maxPositionA = static_cast<long int>(std::floor(maxA));
	long int minPositionB = static_cast<long int>(refining ? std::floor(minB) : std::ceil(minB));
	long int maxPositionB = static_cast<long int>(std::floor(maxB));
	if (!cellInfo.isPeriodicA()) {
		if (minPositionA < 0)
			minPositionA = 0;
		else if (static_cast<std::size_t>(minPositionA) >= gridSizes[0])
			minPositionA = static_cast<long int>(gridSizes[0]) - 1;
		if (maxPositionA < 0)
			maxPositionA = 0;
		else if (static_cast<std::size_t>(maxPositionA) >= gridSizes[0])
			maxPositionA = static_cast<long int>(gridSizes[0]) - 1;
	}
	if (!cellInfo.isPeriodicB()) {
		if (minPositionB < 0)
			minPositionB = 0;
		else if (static_cast<std::size_t>(minPositionB) >= gridSizes[1])
			minPositionB = static_cast<long int>(gridSizes[1]) - 1;
		if (maxPositionB < 0)
			maxPositionB = 0;
		else if (static_cast<std::size_t>(maxPositionB) >= gridSizes[1])
			maxPositionB = static_cast<long int>(gridSizes[1]) - 1;
	}
	return std::make_tuple(minPositionA, maxPositionA, minPositionB, maxPositionB);
}

std::vector<std::vector<std::size_t>> CalcVoronoiIntegrals::p_integrationGetFaceCandidates(VoronoiTessellation &tessellation, long int minPositionA, long int maxPositionA, long int minPositionB, long int maxPositionB, bool refining) const {
	const std::size_t numRays = static_cast<std::size_t>(maxPositionA - minPositionA + 1) * static_cast<std::size_t>(maxPositionB - minPositionB + 1);
	std::vector<std::vector<std::size_t>> faceCandidates(numRays);
	for (auto &&f: tessellation.m_faces)
		f.setRayCandidates(faceCandidates, minPositionA, maxPositionA, minPositionB, maxPositionB, refining);
	return faceCandidates;
}

std::vector<std::pair<double, std::size_t>> CalcVoronoiIntegrals::p_integrationCalcIntersectionsRayC(const VoronoiTessellation &tessellation, double coordA, double coordB, const std::vector<std::size_t> &faceCandidates) const {
	std::vector<std::pair<double, std::size_t>> facesHit;
	for (auto &&f: faceCandidates) {
		const VoronoiTessellation::Face &face = tessellation.m_faces[f];
		double coordC;
		if (face.checkIntersectionRange(coordA, coordB)) {
			if (face.calculateIntersectionRayC(coordA, coordB, coordC, m_epsilon)) {
				if (coordC - std::floor(coordC) < m_epsilon)
					coordC -= m_epsilon;
				facesHit.emplace_back(coordC, f);
			}
		}
	}
	if (!facesHit.empty()) {
		std::sort(facesHit.begin(), facesHit.end());
		if (tessellation.m_faces[facesHit.front().second].getCellIdBackward() != std::numeric_limits<std::size_t>::max())
			throw std::runtime_error("CalcVoronoiIntegrals: Ray does not begin in the void");
		if (tessellation.m_faces[facesHit.back().second].getCellIdForward() != std::numeric_limits<std::size_t>::max())
			throw std::runtime_error("CalcVoronoiIntegrals: Ray does not end in the void");
	}
	return facesHit;
}

std::tuple<long int, long int, std::size_t> CalcVoronoiIntegrals::p_integrationInitLoopC(bool periodicC, std::size_t gridSizeC, double firstHit) {
	long int positionC = static_cast<long int>(std::ceil(firstHit));
	long int gridShift = 0;
	if (periodicC)
		gridShift = static_cast<long int>(std::floor(static_cast<double>(positionC) / gridSizeC)) * static_cast<long int>(gridSizeC);
	else if (positionC < 0)
		positionC = 0;
	if (gridShift > positionC)
		throw std::runtime_error("CalcVoronoiIntegrals: Invalid grid shift");
	std::size_t gridPositionC = static_cast<std::size_t>(positionC - gridShift);
	if (gridPositionC >= gridSizeC)
		throw std::runtime_error("CalcVoronoiIntegrals: Grid position out of range");
	return std::make_tuple(positionC, gridShift, gridPositionC);
}

std::tuple<double, long int, std::size_t, unsigned int> CalcVoronoiIntegrals::p_integrationRefinedInitLoopC(bool periodicC, std::size_t gridSizeC, unsigned int refiningC, double firstHit) {
	double positionC = std::ceil(firstHit * refiningC) / refiningC;
	long int gridShift = 0;
	if (periodicC)
		gridShift = static_cast<long int>(std::floor(positionC / gridSizeC)) * static_cast<long int>(gridSizeC);
	else if (positionC < 0.0)
		positionC = 0.0;
	if (gridShift > positionC)
		throw std::runtime_error("CalcVoronoiIntegrals: Invalid grid shift");
	std::size_t gridPositionC = static_cast<std::size_t>(std::floor(positionC - gridShift));
	if (gridPositionC >= gridSizeC)
		throw std::runtime_error("CalcVoronoiIntegrals: Grid position out of range");
	unsigned int refiningPositionC = static_cast<unsigned int>(std::floor((positionC - gridShift - gridPositionC) * refiningC));
	if (refiningPositionC >= refiningC)
		throw std::runtime_error("CalcVoronoiIntegrals: Refining position out of range");
	return std::make_tuple(positionC, gridShift, gridPositionC, refiningPositionC);
}

void CalcVoronoiIntegrals::p_integrationUpdateFaceIterator(const VoronoiTessellation &tessellation, const std::vector<std::pair<double, std::size_t>> &facesHit, double positionC, std::vector<std::pair<double, std::size_t>>::const_iterator &faceIt, std::size_t &cellId) {
	while (positionC >= faceIt->first) {
		cellId = tessellation.m_faces[faceIt->second].getCellIdForward();
		++faceIt;
		if (faceIt == facesHit.cend())
			break;
		if (tessellation.m_faces[faceIt->second].getCellIdBackward() != cellId)
			throw std::runtime_error("CalcVoronoiIntegrals: Cell mismatch between subsequent faces");
	}
}

bool CalcVoronoiIntegrals::p_integrationUpdateGridPosition(bool periodicC, std::size_t gridSizeC, long int &positionC, long int &gridShift, std::size_t &gridPositionC) {
	++gridPositionC;
	if (gridPositionC >= gridSizeC) {
		if (periodicC) {
			gridPositionC = 0;
			gridShift += static_cast<long int>(gridSizeC);
		} else {
			return false;
		}
	}
	positionC = static_cast<long int>(gridPositionC) + gridShift;
	return true;
}

bool CalcVoronoiIntegrals::p_integrationRefinedUpdateGridPosition(bool periodicC, std::size_t gridSizeC, unsigned int refiningC, const std::array<std::size_t, 4> &indexB, double &positionC, long int &gridShift, std::size_t &gridPositionC, unsigned int &refiningPositionC, std::size_t &indexC2, std::array<std::size_t, 8> &indexC) {
	++refiningPositionC;
	if (refiningPositionC >= refiningC) {
		refiningPositionC = 0;
		++gridPositionC;
		if (gridPositionC >= gridSizeC) {
			if (periodicC) {
				gridPositionC = 0;
				gridShift += static_cast<long int>(gridSizeC);
			} else {
				return false;
			}
		}
		indexC2 = (gridPositionC + 1) % gridSizeC;
		indexC = { indexB[0] + gridPositionC, indexB[0] + indexC2, indexB[1] + gridPositionC, indexB[1] + indexC2, indexB[2] + gridPositionC, indexB[2] + indexC2, indexB[3] + gridPositionC, indexB[3] + indexC2 };
	}
	positionC = static_cast<double>(gridPositionC) + gridShift + static_cast<double>(refiningPositionC) / refiningC;
	return true;
}

}

}

#endif


