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


#ifndef LT_SNAPSHOT_H
#define LT_SNAPSHOT_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "datagrid.h"
#include "vector.h"

#include <cstddef>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace Libtravis {

namespace Travis {

struct SnapshotDataTypeBase {};

struct AtomNames : SnapshotDataTypeBase { std::vector<std::string> data; };
struct AtomicNumbers : SnapshotDataTypeBase { std::vector<unsigned int> data; };
struct AtomRadii : SnapshotDataTypeBase { std::vector<double> data; };
struct AtomPositions : SnapshotDataTypeBase { std::vector<Vector3> data; };
struct AtomVelocities : SnapshotDataTypeBase { std::vector<Vector3> data; };
struct ElectronDensity : SnapshotDataTypeBase { DataGrid3D<double> data; };
struct ElectricCurrentDensity : SnapshotDataTypeBase { DataGrid3D<Vector3> data; };

template<class... SnapshotData> class Trajectory;

template<class... SnapshotData>
class Snapshot
{
	friend class Trajectory<SnapshotData...>;
	
public:
	std::size_t getId() const { return m_id; }
	
	template<class DataType>
	std::enable_if_t<std::is_base_of<SnapshotDataTypeBase, DataType>::value, const decltype(std::declval<DataType>().data) &> getData() const { return std::get<DataType>(m_snapshotData).data; }
	template<class DataType>
	std::enable_if_t<!std::is_base_of<SnapshotDataTypeBase, DataType>::value, const DataType &> getData() const { return std::get<DataType>(m_snapshotData); }
	
private:
	std::size_t m_id;
	std::tuple<SnapshotData...> m_snapshotData;
};

}

}

#endif

#endif


