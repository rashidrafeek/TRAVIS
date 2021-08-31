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


#ifndef LT_MATRIX_H
#define LT_MATRIX_H


#include "../config.h"


#ifdef NEW_CHARGEVAR


#include "vector.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace Libtravis {

namespace Travis {

template<class T, std::size_t M, std::size_t N, template<class, std::size_t, std::size_t> class Derived>
class MatrixBase
{
	static_assert(M > 0 && N > 0, "Matrix: Size must not be 0");
	
public:
	MatrixBase() {}
	MatrixBase(const std::array<T, M * N> &data) : m_data(data) {}
	
	T &operator[](std::size_t pos) {
		return m_data[pos];
	}
	constexpr const T&operator[](std::size_t pos) const {
		return m_data[pos];
	}
	T &at(std::size_t pos) {
		return m_data.at(pos);
	}
	constexpr const T&at(std::size_t pos) const {
		return m_data.at(pos);
	}
	T &at(std::size_t row, std::size_t col) {
		if (row > M || col > N)
			throw std::out_of_range("Matrix: Range error");
		return m_data[N * row + col];
	}
	constexpr const T&at(std::size_t row, std::size_t col) const {
		if (row > M || col > N)
			throw std::out_of_range("Matrix: Range error");
		return m_data[N * row + col];
	}
	
	Derived<T, M, N> &operator+=(const Derived<T, M, N> &other) {
		for (std::size_t i = 0; i < M * N; ++i)
			m_data[i] += other.m_data[i];
		return *static_cast<Derived<T, M, N> *>(this);
	}
	Derived<T, M, N> &operator-=(const Derived<T, M, N> &other) {
		for (std::size_t i = 0; i < M * N; ++i)
			m_data[i] -= other.m_data[i];
		return *static_cast<Derived<T, M, N> *>(this);
	}
	Derived<T, M, N> &operator*=(const T &scalar) {
		for (std::size_t i = 0; i < M * N; ++i)
			m_data[i] *= scalar;
		return *static_cast<Derived<T, M, N> *>(this);
	}
	Derived<T, M, N> &operator/=(const T &scalar) {
		for (std::size_t i = 0; i < M * N; ++i)
			m_data[i] /= scalar;
		return *static_cast<Derived<T, M, N> *>(this);
	}
	
	Derived<T, M, N> operator+() const {
		return Derived<T, M, N>(m_data);
	}
	Derived<T, M, N> operator-() const {
		Derived<T, M, N> m;
		for (std::size_t i = 0; i < M * N; ++i)
			m.m_data[i] = -m_data[i];
		return m;
	}
	Derived<T, M, N> operator+(const Derived<T, M, N> &other) const {
		Derived<T, M, N> m;
		for (std::size_t i = 0; i < M * N; ++i)
			m.m_data[i] = m_data[i] + other.m_data[i];
		return m;
	}
	Derived<T, M, N> operator-(const Derived<T, M, N> &other) const {
		Derived<T, M, N> m;
		for (std::size_t i = 0; i < M * N; ++i)
			m.m_data[i] = m_data[i] - other.m_data[i];
		return m;
	}
	friend Derived<T, M, N> operator*(const T &scalar, Derived<T, M, N> matrix) {
		return matrix *= scalar;
	}
	friend Derived<T, M, N> operator*(Derived<T, M, N> matrix, const T &scalar) {
		return matrix *= scalar;
	}
	Vector<T, M> operator*(const Vector<T, N> &vector) const {
		Vector<T, M> v;
		for (std::size_t i = 0; i < M; ++i) {
			v[i] = m_data[N * i] * vector[0];
			for (std::size_t j = 1; j < N; ++j)
				v[i] += m_data[N * i + j] * vector[j];
		}
		return v;
	}
	template<std::size_t O>
	Derived<T, M, O> operator*(const Derived<T, N, O> &other) const {
		Derived<T, M, O> m;
		for (std::size_t i = 0; i < M; ++i) {
			for (std::size_t j = 0; j < O; ++j) {
				m[O * i + j] = m_data[N * i] * other[j];
				for (std::size_t k = 1; k < N; ++k) {
					m[O * i + j] += m_data[N * i + k] * other[O * k + j];
				}
			}
		}
		return m;
	}
	Derived<T, M, N> operator/(const T &scalar) const {
		Derived<T, M, N> m;
		for (std::size_t i = 0; i < M * N; ++i)
			m.m_data[i] = m_data[i] / scalar;
		return m;
	}
	
protected:
	std::array<T, M * N> m_data;
};

template<class T, std::size_t M, std::size_t N>
class Matrix : public MatrixBase<T, M, N, Matrix>
{
public:
	Matrix() {}
	Matrix(const std::array<T, M * N> &data) : MatrixBase<T, M, N, Travis::Matrix>(data) {}
};

template<class T>
class Matrix<T, 2, 2> : public MatrixBase<T, 2, 2, Matrix>
{
public:
	Matrix() {}
	Matrix(const std::array<T, 4> &data) : MatrixBase<T, 2, 2, Travis::Matrix>(data) {}
	Matrix(const T &xx, const T &xy, const T &yx, const T &yy) : MatrixBase<T, 2, 2, Travis::Matrix>({xx, xy, yx, yy}) {}
	
	T determinant() const {
		return this->m_data[0] * this->m_data[3] - this->m_data[1] * this->m_data[2];
	}
	Matrix<T, 2, 2> inverseMatrix() const {
		Matrix<T, 2, 2> m;
		T det = determinant();
		m.m_data[0] = this->m_data[3] / det;
		m.m_data[1] = -this->m_data[1] / det;
		m.m_data[2] = -this->m_data[2] / det;
		m.m_data[3] = this->m_data[0] / det;
		return m;
	}
	void invert() {
		*this = inverseMatrix();
	}
};

using MatrixS2 = Matrix<double, 2, 2>;

template<class T>
class Matrix<T, 3, 3> : public MatrixBase<T, 3, 3, Matrix>
{
public:
	Matrix() {}
	Matrix(const std::array<T, 9> &data) : MatrixBase<T, 3, 3, Travis::Matrix>(data) {}
	Matrix(const T &xx, const T &xy, const T &xz, const T &yx, const T &yy, const T &yz, const T &zx, const T &zy, const T &zz) : MatrixBase<T, 3, 3, Travis::Matrix>({xx, xy, xz, yx, yy, yz, zx, zy, zz}) {}
	
	T determinant() const {
		return this->m_data[0] * (this->m_data[4] * this->m_data[8] - this->m_data[5] * this->m_data[7]) + this->m_data[1] * (this->m_data[5] * this->m_data[6] - this->m_data[3] * this->m_data[8]) + this->m_data[2] * (this->m_data[3] * this->m_data[7] - this->m_data[4] * this->m_data[6]);
	}
	Matrix<T, 3, 3> inverseMatrix() const {
		Matrix<T, 3, 3> m;
		T det = determinant();
		m.m_data[0] = (this->m_data[4] * this->m_data[8] - this->m_data[5] * this->m_data[7]) / det;
		m.m_data[1] = (this->m_data[2] * this->m_data[7] - this->m_data[1] * this->m_data[8]) / det;
		m.m_data[2] = (this->m_data[1] * this->m_data[5] - this->m_data[2] * this->m_data[4]) / det;
		m.m_data[3] = (this->m_data[5] * this->m_data[6] - this->m_data[3] * this->m_data[8]) / det;
		m.m_data[4] = (this->m_data[0] * this->m_data[8] - this->m_data[2] * this->m_data[6]) / det;
		m.m_data[5] = (this->m_data[2] * this->m_data[3] - this->m_data[0] * this->m_data[5]) / det;
		m.m_data[6] = (this->m_data[3] * this->m_data[7] - this->m_data[4] * this->m_data[6]) / det;
		m.m_data[7] = (this->m_data[1] * this->m_data[6] - this->m_data[0] * this->m_data[7]) / det;
		m.m_data[8] = (this->m_data[0] * this->m_data[4] - this->m_data[1] * this->m_data[3]) / det;
		return m;
	}
	void invert() {
		*this = inverseMatrix();
	}
	
	static Matrix<T, 3, 3> unityMatrix() {
		return Matrix<T, 3, 3>(T(1), T(0), T(0), T(0), T(1), T(0), T(0), T(0), T(1));
	}
	static Matrix<T, 3, 3> rotationMatrix(Vector<T, 3> axis, const T &angle) {
		axis.normalize();
		Matrix<T, 3, 3> m;
		T s = std::sin(angle);
		T c = std::cos(angle);
		m.m_data[0] = axis[0] * axis[0] * (T(1.0) - c) + c;
		m.m_data[1] = axis[0] * axis[1] * (T(1.0) - c) - axis[2] * s;
		m.m_data[2] = axis[0] * axis[2] * (T(1.0) - c) + axis[1] * s;
		m.m_data[3] = axis[0] * axis[1] * (T(1.0) - c) + axis[2] * s;
		m.m_data[4] = axis[1] * axis[1] * (T(1.0) - c) + c;
		m.m_data[5] = axis[1] * axis[2] * (T(1.0) - c) - axis[0] * s;
		m.m_data[6] = axis[0] * axis[2] * (T(1.0) - c) - axis[1] * s;
		m.m_data[7] = axis[1] * axis[2] * (T(1.0) - c) + axis[0] * s;
		m.m_data[8] = axis[2] * axis[2] * (T(1.0) - c) + c;
		return m;
	}
};

using MatrixS3 = Matrix<double, 3, 3>;

}

}

#endif

#endif


