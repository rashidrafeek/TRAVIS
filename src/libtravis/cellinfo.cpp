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


#include "cellinfo.h"

namespace Libtravis {

namespace Travis {

Vector3 CellInfo::foldVector(const Vector3 &vector) const {
	if (!m_periodicA && !m_periodicB && !m_periodicC)
		return vector;
	
	Vector3 v = m_inverseCellMatrix * vector;
	if (m_periodicA) {
		while (v[0] > 0.5)
			v[0] -= 1.0;
		while (v[0] <= -0.5)
			v[0] += 1.0;
	}
	if (m_periodicB) {
		while (v[1] > 0.5)
			v[1] -= 1.0;
		while (v[1] <= -0.5)
			v[1] += 1.0;
	}
	if (m_periodicC) {
		while (v[2] > 0.5)
			v[2] -= 1.0;
		while (v[2] <= -0.5)
			v[2] += 1.0;
	}
	return m_cellMatrix * v;
}

}

}

#endif


