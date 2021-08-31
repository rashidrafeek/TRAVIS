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


#ifndef LT_DATAGRID_H
#define LT_DATAGRID_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "vector.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace Libtravis {

namespace Travis {

template<class T, unsigned int Dim>
class DataGridBase
{
	static_assert(Dim > 0, "DataGrid: Dimension must not be 0");
	
	friend class FileReaderCube;
	
public:
	DataGridBase() {}
	DataGridBase(const std::array<std::size_t, Dim> &sizes) {
		setSizes(sizes);
	}
	
	T &at(const std::array<std::size_t, Dim> &index) {
		if (index[0] >= m_sizes[0])
			throw std::out_of_range("DataGrid: Range error");
		std::size_t totalIndex = index[0];
		for (unsigned int i = 1; i < Dim; ++i) {
			if (index[i] >= m_sizes[i])
				throw std::out_of_range("DataGrid: Range error");
			totalIndex *= m_sizes[i];
			totalIndex += index[i];
		}
		return m_data[totalIndex];
	}
	const T &at(const std::array<std::size_t, Dim> &index) const {
		if (index[0] >= m_sizes[0])
			throw std::out_of_range("DataGrid: Range error");
		std::size_t totalIndex = index[0];
		for (unsigned int i = 1; i < Dim; ++i) {
			if (index[i] >= m_sizes[i])
				throw std::out_of_range("DataGrid: Range error");
			totalIndex *= m_sizes[i];
			totalIndex += index[i];
		}
		return m_data[totalIndex];
	}
	T &at(std::size_t index) {
		return m_data.at(index);
	}
	constexpr const T &at(std::size_t index) const {
		return m_data.at(index);
	}
	T &operator[](std::size_t index) {
		return m_data[index];
	}
	constexpr const T &operator[](std::size_t index) const {
		return m_data[index];
	}
	
	const std::array<std::size_t, Dim> &getSizes() const {
		return m_sizes;
	}
	void setSizes(const std::array<std::size_t, Dim> &sizes) {
		m_sizes = sizes;
		std::size_t totalSize = m_sizes[0];
		for (unsigned int i = 1; i < Dim; ++i)
			totalSize *= m_sizes[i];
		m_data.resize(totalSize);
	}
	
	void fill(const T &value) {
		std::fill(m_data.begin(), m_data.end(), value);
	}
	
	static std::size_t foldToSize(long int index, std::size_t size) {
		if (index >= 0)
			return static_cast<std::size_t>(index) % size;
		else
			return (size - static_cast<std::size_t>(-index) % size) % size;
	}
	std::size_t foldIndex(long int index, unsigned int dim) const {
		return foldToSize(index, m_sizes[dim]);
	}
	
private:
	std::array<std::size_t, Dim> m_sizes;
	std::vector<T> m_data;
};

template<class T, unsigned int Dim>
class DataGrid : public DataGridBase<T, Dim>
{
public:
	DataGrid() : DataGridBase<T, Dim>() {}
	DataGrid(const std::array<std::size_t, Dim> &sizes) : DataGridBase<T, Dim>(sizes) {}
};

template<class T>
class DataGrid<T, 3> : public DataGridBase<T, 3>
{
public:
	DataGrid() : DataGridBase<T, 3>() {}
	DataGrid(const std::array<std::size_t, 3> &sizes) : DataGridBase<T, 3>(sizes) {}
	DataGrid(std::size_t size1, std::size_t size2, std::size_t size3) : DataGridBase<T, 3>({size1, size2, size3}) {}
	
	using DataGridBase<T, 3>::at;
	T &at(std::size_t index1, std::size_t index2, std::size_t index3) {
		return DataGridBase<T, 3>::at({index1, index2, index3});
	}
	const T &at(std::size_t index1, std::size_t index2, std::size_t index3) const {
		return DataGridBase<T, 3>::at({index1, index2, index3});
	}
	
	void setSizes(std::size_t size1, std::size_t size2, std::size_t size3) {
		DataGridBase<T, 3>::setSizes({size1, size2, size3});
	}
};

template<class T>
using DataGrid3D = DataGrid<T, 3>;

}

}

#endif

#endif


