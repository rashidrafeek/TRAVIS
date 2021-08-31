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


#ifndef LT_CALCVORONOI_H
#define LT_CALCVORONOI_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "cellinfo.h"
#include "datacalculator.h"
#include "vector.h"
#include "voronoitessellation.h"

#include <functional>
#include <stdexcept>
#include <vector>

namespace Libtravis {

namespace Travis {

class CalcVoronoi : public DataCalculator
{
public:
	void setAtomPositionsAccessor(const std::function<std::vector<Vector3> &()> &atomPositionsAccessor) { m_atomPositionsAccessor = atomPositionsAccessor; }
	void setAtomRadiiAccessor(const std::function<std::vector<double> &()> &atomRadiiAccessor) { m_atomRadiiAccessor = atomRadiiAccessor; }
	void setCellInfoAccessor(const std::function<CellInfo &()> &cellInfoAccessor) { m_cellInfoAccessor = cellInfoAccessor; }
	void setVoronoiTessellationAccessor(const std::function<VoronoiTessellation &()> &voronoiTessellationAccessor) { m_voronoiTessellationAccessor = voronoiTessellationAccessor; }
	
	void setVertexIdentityThreshold(double threshold) { m_vertexIdentityThreshold = threshold; }
	
	VoronoiTessellation calcVoronoi(const std::vector<Vector3> &atomPositions, const std::vector<double> &atomRadii, const CellInfo &cellInfo = CellInfo()) const;
	
private:
	std::function<std::vector<Vector3> &()> m_atomPositionsAccessor;
	std::function<std::vector<double> &()> m_atomRadiiAccessor;
	std::function<CellInfo &()> m_cellInfoAccessor;
	std::function<VoronoiTessellation &()> m_voronoiTessellationAccessor;
	
	double m_vertexIdentityThreshold = 1e-10;

	virtual void p_calculateNextStep() override {
		if (!m_atomPositionsAccessor)
			throw std::runtime_error("CalcVoronoi: source for atom positions not set");
		if (!m_atomRadiiAccessor)
			throw std::runtime_error("CalcVoronoi: source for atom radii not set");
		if (!m_voronoiTessellationAccessor)
			throw std::runtime_error("CalcVoronoi: target for Voronoi tessellation not set");
		
		if (!m_cellInfoAccessor)
			m_voronoiTessellationAccessor() = calcVoronoi(m_atomPositionsAccessor(), m_atomRadiiAccessor());
		else
			m_voronoiTessellationAccessor() = calcVoronoi(m_atomPositionsAccessor(), m_atomRadiiAccessor(), m_cellInfoAccessor());
	}
};

}

}

#endif

#endif


