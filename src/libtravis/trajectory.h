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


#ifndef LT_TRAJECTORY_H
#define LT_TRAJECTORY_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "datacalculator.h"
#include "datasource.h"
#include "snapshot.h"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace Libtravis {

namespace Travis {

template<class... SnapshotData>
class Trajectory
{
public:
	class DataSourceEnded : public std::runtime_error
	{
	public:
		DataSourceEnded() : std::runtime_error("Trajectory: Data source ended before requested snapshot could be reached") {}
	};
	
	class Iterator
	{
	public:
		Iterator(Trajectory<SnapshotData...> &trajectory, bool endOfTrajectory) : m_trajectory(trajectory), m_endOfTrajectory(endOfTrajectory) {
			if (!m_endOfTrajectory)
				if (!m_trajectory.hasSnapshot(m_snapshotId))
					m_endOfTrajectory = true;
		}
		
		const Snapshot<SnapshotData...> &operator*() const { return m_trajectory.getSnapshot(m_snapshotId); }
		const Snapshot<SnapshotData...> *operator->() const { return &m_trajectory.getSnapshot(m_snapshotId); }
		Iterator &operator++() {
			++m_snapshotId;
			if (!m_trajectory.hasSnapshot(m_snapshotId))
				m_endOfTrajectory = true;
			return *this;
		}
		Iterator &operator+=(std::size_t increment) {
			m_snapshotId += increment;
			if (!m_trajectory.hasSnapshot(m_snapshotId))
				m_endOfTrajectory = true;
			return *this;
		}
		Iterator operator+(std::size_t increment) const {
			Iterator it(*this);
			it += increment;
			return it;
		}
		
		bool operator==(const Iterator &other) const {
			if (m_endOfTrajectory)
				return other.m_endOfTrajectory;
			else
				return (m_snapshotId == other.m_snapshotId) && !other.m_endOfTrajectory;
		}
		bool operator!=(const Iterator &other) const {
			return !(*this == other);
		}
		
	private:
		Trajectory<SnapshotData...> &m_trajectory;
		bool m_endOfTrajectory;
		std::size_t m_snapshotId = 0;
	};
	
	using SnapshotType = Snapshot<SnapshotData...>;
	
	Iterator begin() { return Iterator(*this, false); }
	Iterator end() { return Iterator(*this, true); }
	
	template<class DataType, long int shift = 0>
	std::enable_if_t<std::is_base_of<SnapshotDataTypeBase, DataType>::value, std::function<decltype(std::declval<DataType>().data) &()>> getAccessor() {
		return [this] () -> decltype(std::declval<DataType>().data) & { return std::get<DataType>(m_snapshotBuffer[(m_bufferAccessorIndex + m_snapshotBuffer.size() + static_cast<std::size_t>(shift)) % m_snapshotBuffer.size()].m_snapshotData).data; };
	}
	template<class DataType, long int shift = 0>
	std::enable_if_t<!std::is_base_of<SnapshotDataTypeBase, DataType>::value, std::function<DataType &()>> getAccessor() {
		return [this] () -> DataType & { return std::get<DataType>(m_snapshotBuffer[(m_bufferAccessorIndex + m_snapshotBuffer.size() + static_cast<std::size_t>(shift)) % m_snapshotBuffer.size()].m_snapshotData); };
	}
	
	const Snapshot<SnapshotData...> &getSnapshot(std::size_t snapshotId) {
		if (m_needReset || snapshotId + 1 != m_nextSnapshotId) {
			if (!p_prepareSnapshot(snapshotId)) {
				m_needReset = true;
				throw DataSourceEnded();
			}
		}
		return m_snapshotBuffer[m_bufferIndex];
	}
	const Snapshot<SnapshotData...> &operator[](std::size_t snapshotId) { return getSnapshot(snapshotId); }
	
	bool hasSnapshot(std::size_t snapshotId) {
		if (snapshotId < m_nextSnapshotId)
			return true;
		if (!p_prepareSnapshot(snapshotId)) {
			m_needReset = true;
			return false;
		}
		return true;
	}
	
	void reset() {
		m_nextSnapshotId = 0;
		for (auto &&s: m_dataSources)
			s->rewind();
		
		m_bufferAccessorIndex = 0;
		for (std::size_t i = 0; i < m_snapshotBuffer.size() - 1; ++i) {
			for (auto &&s: m_dataSources)
				if (!s->readNextStep())
					throw DataSourceEnded();
			++m_bufferAccessorIndex;
		}
		m_bufferIndex = (2 * m_snapshotBuffer.size() - m_bufferSizeForward - 2) % m_snapshotBuffer.size();
		m_needReset = false;
	}
	
	template<class T, class... Args>
	T &addDataSource(Args &&... args) {
		m_needReset = true;
		std::unique_ptr<T> ptr = std::make_unique<T>(std::forward<Args>(args)...);
		T &ref = *ptr;
		m_dataSources.push_back(std::move(ptr));
		return ref;
	}
	template<class T, class... Args>
	T &addDataCalculator(Args &&... args) {
		m_needReset = true;
		std::unique_ptr<T> ptr = std::make_unique<T>(std::forward<Args>(args)...);
		T &ref = *ptr;
		m_dataCalculators.push_back(std::move(ptr));
		return ref;
	}
	
	void setSnapshotBufferSize(std::size_t forward, std::size_t backward) {
		m_bufferSizeForward = forward;
		m_snapshotBuffer.resize(forward + backward + 1);
		m_needReset = true;
	}
	
private:
	std::size_t m_nextSnapshotId = 0;
	std::vector<Snapshot<SnapshotData...>> m_snapshotBuffer{1};
	std::size_t m_bufferIndex = 0;
	std::size_t m_bufferAccessorIndex = 0;
	std::size_t m_bufferSizeForward = 0;
	bool m_needReset = false;
	
	std::vector<std::unique_ptr<DataSource>> m_dataSources;
	std::vector<std::unique_ptr<DataCalculator>> m_dataCalculators;
	
	bool p_prepareSnapshot(std::size_t snapshotId) {
		if (m_needReset || snapshotId + 1 < m_nextSnapshotId)
			reset();
		
		if (snapshotId + 1 > m_nextSnapshotId) {
			if (snapshotId > m_nextSnapshotId + m_snapshotBuffer.size() - 1) {
				std::size_t numSteps = snapshotId + 1 - m_nextSnapshotId - m_snapshotBuffer.size();
				m_nextSnapshotId += numSteps;
				for (auto &&s: m_dataSources)
					if (!s->skipSteps(numSteps))
						return false;
			}
			
			while (snapshotId > m_nextSnapshotId) {
				m_bufferIndex = (m_bufferIndex + 1) % m_snapshotBuffer.size();
				m_bufferAccessorIndex = (m_bufferIndex + m_bufferSizeForward) % m_snapshotBuffer.size();
				for (auto &&s: m_dataSources)
					if (!s->readNextStep())
						return false;
				m_snapshotBuffer[m_bufferIndex].m_id = m_nextSnapshotId;
				++m_nextSnapshotId;
			}
			
			m_bufferIndex = (m_bufferIndex + 1) % m_snapshotBuffer.size();
			m_bufferAccessorIndex = (m_bufferIndex + m_bufferSizeForward) % m_snapshotBuffer.size();
			for (auto &&s: m_dataSources)
				if (!s->readNextStep())
					return false;
			m_bufferAccessorIndex = m_bufferIndex;
			for (auto &&c: m_dataCalculators)
				c->calculateNextStep();
			
			m_snapshotBuffer[m_bufferIndex].m_id = m_nextSnapshotId;
			++m_nextSnapshotId;
		}
		
		return true;
	}
};

}

}

#endif

#endif


