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


#ifndef LT_VECTOR_H
#define LT_VECTOR_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include <array>
#include <cmath>
#include <cstddef>

namespace Libtravis {

namespace Travis {

template<class T, std::size_t N, template<class, std::size_t> class Derived>
class VectorBase
{
	static_assert(N > 0, "Vector: Dimension must not be 0");
	
public:
	VectorBase() {}
	VectorBase(const std::array<T, N> &data) : m_data(data) {}
	
	T &operator[](std::size_t pos) {
		return m_data[pos];
	}
	constexpr const T &operator[](std::size_t pos) const {
		return m_data[pos];
	}
	T &at(std::size_t pos) {
		return m_data.at(pos);
	}
	constexpr const T &at(std::size_t pos) const {
		return m_data.at(pos);
	}
	
	Derived<T, N> &operator+=(const Derived<T, N> &other) {
		for (std::size_t i = 0; i < N; ++i)
			m_data[i] += other.m_data[i];
		return *static_cast<Derived<T, N> *>(this);
	}
	Derived<T, N> &operator-=(const Derived<T, N> &other) {
		for (std::size_t i = 0; i < N; ++i)
			m_data[i] -= other.m_data[i];
		return *static_cast<Derived<T, N> *>(this);
	}
	Derived<T, N> &operator*=(const T &scalar) {
		for (std::size_t i = 0; i < N; ++i)
			m_data[i] *= scalar;
		return *static_cast<Derived<T, N> *>(this);
	}
	Derived<T, N> &operator/=(const T &scalar) {
		for (std::size_t i = 0; i < N; ++i)
			m_data[i] /= scalar;
		return *static_cast<Derived<T, N> *>(this);
	}
	
	Derived<T, N> operator+() const {
		return Derived<T, N>(m_data);
	}
	Derived<T, N> operator-() const {
		Derived<T, N> v;
		for (std::size_t i = 0; i < N; ++i)
			v.m_data[i] = -m_data[i];
		return v;
	}
	Derived<T, N> operator+(const Derived<T, N> &other) const {
		Derived<T, N> v;
		for (std::size_t i = 0; i < N; ++i)
			v.m_data[i] = m_data[i] + other.m_data[i];
		return v;
	}
	Derived<T, N> operator-(const Derived<T, N> &other) const {
		Derived<T, N> v;
		for (std::size_t i = 0; i < N; ++i)
			v.m_data[i] = m_data[i] - other.m_data[i];
		return v;
	}
	friend Derived<T, N> operator*(const T &scalar, Derived<T, N> vector) {
		return vector *= scalar;
	}
	friend Derived<T, N> operator*(Derived<T, N> vector, const T &scalar) {
		return vector *= scalar;
	}
	Derived<T, N> operator/(const T &scalar) const {
		Derived<T, N> v;
		for (std::size_t i = 0; i < N; ++i)
			v.m_data[i] = m_data[i] / scalar;
		return v;
	}
	T dot(const Derived<T, N> &other) const {
		T val = m_data[0] * other.m_data[0];
		for (std::size_t i = 1; i < N; ++i)
			val += m_data[i] * other.m_data[i];
		return val;
	}
	
	T norm() const {
		return std::sqrt(dot(*static_cast<const Derived<T, N> *>(this)));
	}
	void normalize() {
		*this /= norm();
	}
	
protected:
	std::array<T, N> m_data;
};

template<class T, std::size_t N>
class Vector : public VectorBase<T, N, Vector>
{
public:
	Vector() {}
	Vector(const std::array<T, N> &data) : VectorBase<T, N, Travis::Vector>(data) {}
};

template<class T>
class Vector<T, 2> : public VectorBase<T, 2, Vector>
{
public:
	Vector() {}
	Vector(const std::array<T, 2> &data) : VectorBase<T, 2, Travis::Vector>(data) {}
	Vector(const T &x, const T &y) : VectorBase<T, 2, Travis::Vector>({x, y}) {}
	
	T &x() {
		return this->m_data[0];
	}
	constexpr const T &x() const {
		return this->m_data[0];
	}
	T &y() {
		return this->m_data[1];
	}
	constexpr const T &y() const {
		return this->m_data[1];
	}
};

using Vector2 = Vector<double, 2>;

template<class T>
class Vector<T, 3> : public VectorBase<T, 3, Vector>
{
public:
	Vector() {}
	Vector(const std::array<T, 3> &data) : VectorBase<T, 3, Travis::Vector>(data) {}
	Vector(const T &x, const T &y, const T &z) : VectorBase<T, 3, Travis::Vector>({x, y, z}) {}
	
	T &x() {
		return this->m_data[0];
	}
	constexpr const T &x() const {
		return this->m_data[0];
	}
	T &y() {
		return this->m_data[1];
	}
	constexpr const T &y() const {
		return this->m_data[1];
	}
	T &z() {
		return this->m_data[2];
	}
	constexpr const T &z() const {
		return this->m_data[2];
	}
	
	Vector<T, 3> cross(const Vector<T, 3> &other) const {
		return Vector<T, 3>(this->m_data[1] * other.m_data[2] - this->m_data[2] * other.m_data[1], this->m_data[2] * other.m_data[0] - this->m_data[0] * other.m_data[2], this->m_data[0] * other.m_data[1] - this->m_data[1] * other.m_data[0]);
	}
};

using Vector3 = Vector<double, 3>;

}

}

#endif

#endif


