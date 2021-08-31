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


#ifndef LT_CALCVORONOIINTEGRALS_H
#define LT_CALCVORONOIINTEGRALS_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "datacalculator.h"
#include "datagrid.h"
#include "vector.h"

#include <array>
#include <bitset>
#include <cstddef>
#include <functional>
#include <tuple>
#include <utility>
#include <vector>

namespace Libtravis {

namespace Travis {

class CellInfo;
class VoronoiTessellation;

class CalcVoronoiIntegrals : public DataCalculator
{
public:
	class IntegrationMask
	{
	public:
		class RaySegment
		{
			friend class IntegrationMask;
			
		public:
			RaySegment(std::size_t indexA, std::size_t indexB, std::size_t indexC) : m_indexA(indexA), m_indexB(indexB), m_segmentC(indexC, indexC) {}
			RaySegment(std::size_t indexA, unsigned int refIndexA, std::size_t indexB, unsigned int refIndexB, std::size_t indexC, unsigned int refIndexC) : m_indexA(indexA), m_refIndexA(refIndexA), m_indexB(indexB), m_refIndexB(refIndexB), m_segmentC(indexC, indexC), m_refSegmentC(refIndexC, refIndexC) {}
			
			bool isInRange(std::size_t indexC) const { return indexC >= m_segmentC.first && indexC <= m_segmentC.second; }
			bool isInRange(std::size_t indexC, unsigned int refIndexC) const;
			
		private:
			std::size_t m_indexA;
			unsigned int m_refIndexA = 0;
			std::size_t m_indexB;
			unsigned int m_refIndexB = 0;
			std::pair<std::size_t, std::size_t> m_segmentC;
			std::pair<unsigned int, unsigned int> m_refSegmentC{0, 0};
		};
		
		class RaySegmentation
		{
			friend class IntegrationMask;
			
		public:
			bool isEmpty() const { return m_segments.empty(); }
			bool isInRange(std::size_t indexC) const;
			bool isInRange(std::size_t indexC, unsigned int refIndexC) const;
			
		private:
			std::vector<std::reference_wrapper<const RaySegment>> m_segments;
		};
		
		IntegrationMask() = default;
		IntegrationMask(const std::array<unsigned int, 3> &refining) : m_refining(refining) {}
		
		const std::array<unsigned int, 3> &getRefining() const { return m_refining; }
		
		RaySegmentation getSegmentation(std::size_t indexA, std::size_t indexB) const;
		RaySegmentation getSegmentation(std::size_t indexA, unsigned int refIndexA, std::size_t indexB, unsigned int refIndexB) const;
		
		void addIndexLast(std::size_t indexA, std::size_t indexB, std::size_t indexC);
		void addIndexLast(std::size_t indexA, unsigned int refIndexA, std::size_t indexB, unsigned int refIndexB, std::size_t indexC, unsigned int refIndexC);
		
	private:
		std::vector<RaySegment> m_raySegments;
		std::array<unsigned int, 3> m_refining{1, 1, 1};
	};
	
	class IntegrationConfig
	{
		friend class CalcVoronoiIntegrals;
		
	public:
		void setCalcElectricMoments(const DataGrid3D<double> &electronDensity, bool calcCharge, bool calcDipole, bool calcQuadrupole);
		void setCalcMagneticMoments(const DataGrid3D<Vector3> &electricCurrentDensity, bool calcCurrent, bool calcDipole);
		void setCalcMaskSet(std::vector<IntegrationMask> &integrationMasks, const std::vector<std::vector<std::size_t>> &maskCellIds);
		void setCalcSanityCheck(DataGrid3D<unsigned int> &sanityData);
		void setCalcSanityCheckRefined(DataGrid3D<double> &sanityData);
		void setIntegrationMask(const IntegrationMask &integrationMask);
		void setRefining(const std::array<unsigned int, 3> &refining) { m_refining = refining; }
		
	private:
		enum Flags { Charge = 0, ElectricDipole = 1, ElectricQuadrupole = 2, Current = 3, MagneticDipole = 4, Mask = 5, Sanity = 6, SanityRefined = 7 };
		std::bitset<8> m_flags;
		const DataGrid3D<double> *m_electronDensity = nullptr;
		const DataGrid3D<Vector3> *m_electricCurrentDensity = nullptr;
		DataGrid3D<unsigned int> *m_sanityData = nullptr;
		DataGrid3D<double> *m_sanityDataRefined = nullptr;
		std::vector<IntegrationMask> *m_integrationMasks = nullptr;
		const std::vector<std::vector<std::size_t>> *m_maskCellIds = nullptr;
		
		bool m_masked = false;
		const IntegrationMask *m_integrationMask = nullptr;
		std::array<unsigned int, 3> m_refining{1, 1, 1};
	};
	
	void setCellInfoAccessor(const std::function<CellInfo &()> &cellInfoAccessor) { m_cellInfoAccessor = cellInfoAccessor; }
	void setVoronoiTessellationAccessor(const std::function<VoronoiTessellation &()> &voronoiTessellationAccessor) { m_voronoiTessellationAccessor = voronoiTessellationAccessor; }
	void setElectronDensityAccessor(const std::function<DataGrid3D<double> &()> &electronDensityAccessor) { m_electronDensityAccessor = electronDensityAccessor; }
	void setElectricCurrentDensityAccessor(const std::function<DataGrid3D<Vector3> &()> &electricCurrentDensityAccessor) { m_electricCurrentDensityAccessor = electricCurrentDensityAccessor; }
	void setMaskSetTargetAccessor(const std::function<std::vector<IntegrationMask> &()> &maskSetTargetAccessor) { m_maskSetTargetAccessor = maskSetTargetAccessor; }
	void setMaskCellIdAccessor(const std::function<std::vector<std::vector<std::size_t>> &()> &maskCellIdAccessor) { m_maskCellIdAccessor = maskCellIdAccessor; }
	
	void setCalculateElectricMoments(bool calculateCharge, bool calculateElectricDipole, bool calculateElectricQuadrupole) { m_calculateCharge = calculateCharge; m_calculateElectricDipole = calculateElectricDipole; m_calculateElectricQuadrupole = calculateElectricQuadrupole; }
	void setCalculateMagneticMoments(bool calculateCurrent, bool calculateMagneticDipole) { m_calculateCurrent = calculateCurrent; m_calculateMagneticDipole = calculateMagneticDipole; }
	void setCalculateMaskSet(bool calculateMaskSet) { m_calculateMaskSet = calculateMaskSet; }
	void setRefining(const std::array<unsigned int, 3> &refining) { m_refining = refining; }
	
	void setEpsilon(double epsilon) { m_epsilon = epsilon; }
	
	void calcIntegrals(VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, const IntegrationConfig &config) const;
	
private:
	using IntegrandType = std::function<void(std::size_t, const Vector3 &)>;
	using PreparatorType = std::function<void(std::size_t, const std::array<std::size_t, 3> &)>;
	using PreparatorRefinedType = std::function<void(const std::array<std::size_t, 8> &, const std::array<double, 8> &, unsigned int, const std::array<std::size_t, 3> &, const std::array<unsigned int, 3> &)>;
	
	std::function<CellInfo &()> m_cellInfoAccessor;
	std::function<VoronoiTessellation &()> m_voronoiTessellationAccessor;
	std::function<DataGrid3D<double> &()> m_electronDensityAccessor;
	std::function<DataGrid3D<Vector3> &()> m_electricCurrentDensityAccessor;
	std::function<std::vector<IntegrationMask> &()> m_maskSetTargetAccessor;
	std::function<std::vector<std::vector<std::size_t>> &()> m_maskCellIdAccessor;
	
	bool m_calculateCharge = false;
	bool m_calculateElectricDipole = false;
	bool m_calculateElectricQuadrupole = false;
	bool m_calculateCurrent = false;
	bool m_calculateMagneticDipole = false;
	bool m_calculateMaskSet = false;
	std::array<unsigned int, 3> m_refining{1, 1, 1};
	
	double m_epsilon = 1e-10;
	
	virtual void p_calculateNextStep() override;
	
	void p_performIntegration(VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, const std::vector<PreparatorType> &preparators, const std::vector<IntegrandType> &integrands) const;
	void p_performIntegrationRefined(VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, const std::array<unsigned int, 3> &refining, const std::vector<PreparatorRefinedType> &preparators, const std::vector<IntegrandType> &integrands) const;
	void p_performIntegrationMasked(VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, const std::vector<PreparatorType> &preparators, const std::vector<IntegrandType> &integrands, const IntegrationMask &integrationMask) const;
	void p_performIntegrationRefinedMasked(VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, const std::array<unsigned int, 3> &refining, const std::vector<PreparatorRefinedType> &preparators, const std::vector<IntegrandType> &integrands, const IntegrationMask &integrationMask) const;
	
	std::tuple<long int, long int, long int, long int> p_integrationCalcBoundsAB(const VoronoiTessellation &tessellation, const CellInfo &cellInfo, const std::array<std::size_t, 3> &gridSizes, bool refining) const;
	std::vector<std::vector<std::size_t>> p_integrationGetFaceCandidates(VoronoiTessellation &tessellation, long int minPositionA, long int maxPositionA, long int minPositionB, long int maxPositionB, bool refining) const;
	std::vector<std::pair<double, std::size_t>> p_integrationCalcIntersectionsRayC(const VoronoiTessellation &tessellation, double coordA, double coordB, const std::vector<std::size_t> &faceCandidates) const;
	static std::tuple<long int, long int, std::size_t> p_integrationInitLoopC(bool periodicC, std::size_t gridSizeC, double firstHit);
	static std::tuple<double, long int, std::size_t, unsigned int> p_integrationRefinedInitLoopC(bool periodicC, std::size_t gridSizeC, unsigned int refiningC, double firstHit);
	static void p_integrationUpdateFaceIterator(const VoronoiTessellation &tessellation, const std::vector<std::pair<double, std::size_t>> &facesHit, double positionC, std::vector<std::pair<double, std::size_t>>::const_iterator &faceIt, std::size_t &cellId);
	static bool p_integrationUpdateGridPosition(bool periodicC, std::size_t gridSizeC, long int &positionC, long int &gridShift, std::size_t &gridPositionC);
	static bool p_integrationRefinedUpdateGridPosition(bool periodicC, std::size_t gridSizeC, unsigned int refiningC, const std::array<std::size_t, 4> &indexB, double &positionC, long int &gridShift, std::size_t &gridPositionC, unsigned int &refiningPositionC, std::size_t &indexC2, std::array<std::size_t, 8> &indexC);
};

}

}

#endif

#endif


