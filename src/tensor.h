/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Martin Brehm.

    ---------------------------------------------------------------------------

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


#ifndef TENSOR_H
#define TENSOR_H


// This must always be the first include directive
#include "config.h"

#include <vector>
#include <algorithm>
#include <functional>
#include "tools.h"


#define TENSOR_BOUND_CHECK


/*****************************************************************************
***   CDTensor1   ************************************************************
*****************************************************************************/

class CDTensor1 {
public:

	CDTensor1() {
	}

	explicit CDTensor1(int i) {
		m_pData.resize(i,0);
	}

	~CDTensor1() {
	}

	CDTensor1(const CDTensor1 &t) : m_pData(t.m_pData) {
	}

	CDTensor1& operator=(const CDTensor1& t) {
		m_pData = t.m_pData;
		return *this;
	}

	double& operator[] (int i) {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= (int)m_pData.size())) {
				eprintf("CDTensor1 &operator [] boundary error (%d/%lu).\n",i,(unsigned long)m_pData.size());
				abort();
			}
		#endif
		return m_pData[i];
	}

	double operator[] (int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= (int)m_pData.size())) {
				eprintf("CDTensor1 &operator [] boundary error (%d/%lu).\n",i,(unsigned long)m_pData.size());
				abort();
			}
		#endif
		return m_pData[i];
	}

	double& operator() (int i) {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= (int)m_pData.size())) {
				eprintf("CDTensor1 &operator () boundary error (%d/%lu).\n",i,(unsigned long)m_pData.size());
				abort();
			}
		#endif
		return m_pData[i];
	}

	double operator() (int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= (int)m_pData.size())) {
				eprintf("CDTensor1 operator () boundary error (%d/%lu).\n",i,(unsigned long)m_pData.size());
				abort();
			}
		#endif
		return m_pData[i];
	}

	CDTensor1& operator += (const CDTensor1 &t) {
		#ifdef TENSOR_BOUND_CHECK
			if (!DimensionsMatch(t)) {
				eprintf("CDTensor1 &operator += dimension mismatch.\n");
				abort();
			}
		#endif
		std::transform(m_pData.begin(),m_pData.end(),t.m_pData.begin(),m_pData.begin(),std::plus<double>());
		return *this;
	}

	CDTensor1& operator -= (const CDTensor1 &t) {
		#ifdef TENSOR_BOUND_CHECK
			if (!DimensionsMatch(t)) {
				eprintf("CDTensor1 &operator -= dimension mismatch.\n");
				abort();
			}
		#endif
		std::transform(m_pData.begin(),m_pData.end(),t.m_pData.begin(),m_pData.begin(),std::minus<double>());
		return *this;
	}

	CDTensor1& operator *= (double f) {
		//std::transform(m_pData.begin(),m_pData.end(),m_pData.begin(),std::bind1st(std::multiplies<double>(), f));
		for (unsigned long z=0;z<m_pData.size();z++)
			m_pData[z] *= f;
		return *this;
	}

	CDTensor1& operator /= (double f) {
		//std::transform(m_pData.begin(),m_pData.end(),m_pData.begin(),std::bind1st(std::divides<double>(), f));
		for (unsigned long z=0;z<m_pData.size();z++)
			m_pData[z] /= f;
		return *this;
	}

	void setZero() {
		std::fill(m_pData.begin(),m_pData.end(),0);
	}

	int getDimension() const {
		return (int)m_pData.size();
	}

	bool DimensionsMatch(const CDTensor1 &t) const {
		if (m_pData.size() != t.m_pData.size()) {
			eprintf("Dimension mismatch %lu <--> %lu.\n",(unsigned long)m_pData.size(),(unsigned long)t.m_pData.size());
			return false;
		}
		return true;
	}

private:
	std::vector<double> m_pData;
};


inline CDTensor1 operator+( CDTensor1 lhs, CDTensor1 const & rhs ) {
   lhs += rhs;
   return lhs;
}


inline CDTensor1 operator-( CDTensor1 lhs, CDTensor1 const & rhs ) {
   lhs -= rhs;
   return lhs;
}


inline CDTensor1 operator*( CDTensor1 lhs, double f ) {
   lhs *= f;
   return lhs;
}


inline CDTensor1 operator/( CDTensor1 lhs, double f ) {
   lhs /= f;
   return lhs;
}


inline CDTensor1 operator*( double f, CDTensor1 lhs ) {
   lhs *= f;
   return lhs;
}


inline CDTensor1 operator/( double f, CDTensor1 lhs ) {
   lhs /= f;
   return lhs;
}



/*****************************************************************************
***   CDTensor2   ************************************************************
*****************************************************************************/

class CDTensor2 {
public:

	CDTensor2() {
		m_iDimensions[0] = 0;
		m_iDimensions[1] = 0;
	}

	CDTensor2(int i, int j) {
		m_iDimensions[0] = i;
		m_iDimensions[1] = j;
		m_pData.resize(i*j,0);
	}

	~CDTensor2() {
	}

	CDTensor2(const CDTensor2 &t) : m_pData(t.m_pData) {
		m_iDimensions[0] = t.m_iDimensions[0];
		m_iDimensions[1] = t.m_iDimensions[1];
	}

	CDTensor2& operator=(const CDTensor2& t) {
		m_pData = t.m_pData;
		m_iDimensions[0] = t.m_iDimensions[0];
		m_iDimensions[1] = t.m_iDimensions[1];
		return *this;
	}

	double& operator[] (int i) {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[0]*m_iDimensions[1])) {
				eprintf("CDTensor2 &operator [] boundary error (%d/%d).\n",i,m_iDimensions[0]*m_iDimensions[1]);
				abort();
			}
		#endif
		return m_pData[i];
	}

	double operator[] (int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[0]*m_iDimensions[1])) {
				eprintf("CDTensor2 &operator [] boundary error (%d/%d).\n",i,m_iDimensions[0]*m_iDimensions[1]);
				abort();
			}
		#endif
		return m_pData[i];
	}

	double& operator() (int i, int j) {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[1])) {
				eprintf("CDTensor2 &operator () boundary error (%d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[1]);
				abort();
			}
		#endif
		return m_pData[i*m_iDimensions[1]+j];
	}

	double operator() (int i, int j) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[1])) {
				eprintf("CDTensor2 operator () boundary error (%d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[1]);
				abort();
			}
		#endif
		return m_pData[i*m_iDimensions[1]+j];
	}

	CDTensor2& operator += (const CDTensor2 &t) {
		#ifdef TENSOR_BOUND_CHECK
			if (!DimensionsMatch(t)) {
				eprintf("CDTensor2 &operator += dimension mismatch.\n");
				abort();
			}
		#endif
		std::transform(m_pData.begin(),m_pData.end(),t.m_pData.begin(),m_pData.begin(),std::plus<double>());
		return *this;
	}

	CDTensor2& operator -= (const CDTensor2 &t) {
		#ifdef TENSOR_BOUND_CHECK
			if (!DimensionsMatch(t)) {
				eprintf("CDTensor2 &operator -= dimension mismatch.\n");
				abort();
			}
		#endif
		std::transform(m_pData.begin(),m_pData.end(),t.m_pData.begin(),m_pData.begin(),std::minus<double>());
		return *this;
	}

	CDTensor2& operator *= (double f) {
		//std::transform(m_pData.begin(),m_pData.end(),m_pData.begin(),std::bind1st(std::multiplies<double>(), f));
		for (unsigned long z=0;z<m_pData.size();z++)
			m_pData[z] *= f;
		return *this;
	}

	CDTensor2& operator /= (double f) {
		//std::transform(m_pData.begin(),m_pData.end(),m_pData.begin(),std::bind1st(std::divides<double>(), f));
		for (unsigned long z=0;z<m_pData.size();z++)
			m_pData[z] /= f;
		return *this;
	}

	CDTensor1 Reduce1(int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[0])) {
				eprintf("CDTensor2 Reduce1() boundary error (%d/%d).\n",i,m_iDimensions[0]);
				abort();
			}
		#endif
		CDTensor1 t(m_iDimensions[1]);
		for (int z=0;z<m_iDimensions[1];z++)
			t(z) = (*this)(i,z);
		return t;
	}

	CDTensor1 Reduce2(int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[1])) {
				eprintf("CDTensor2 Reduce2() boundary error (%d/%d).\n",i,m_iDimensions[1]);
				abort();
			}
		#endif
		CDTensor1 t(m_iDimensions[0]);
		for (int z=0;z<m_iDimensions[0];z++)
			t(z) = (*this)(z,i);
		return t;
	}

	void setZero() {
		std::fill(m_pData.begin(),m_pData.end(),0);
	}

	int getDimension(int i) const {
		return m_iDimensions[i];
	}

	bool DimensionsMatch(const CDTensor2 &t) const {
		if ((m_iDimensions[0] != t.m_iDimensions[0]) ||
			(m_iDimensions[1] != t.m_iDimensions[1])) {
			eprintf("Dimension mismatch ( %d, %d ) <--> ( %d, %d ).\n",m_iDimensions[0],m_iDimensions[1],t.m_iDimensions[0],t.m_iDimensions[1]);
			return false;
		}
		return true;
	}

private:
	std::vector<double> m_pData;
	int m_iDimensions[2];
};


inline CDTensor2 operator+( CDTensor2 lhs, CDTensor2 const & rhs ) {
   lhs += rhs;
   return lhs;
}


inline CDTensor2 operator-( CDTensor2 lhs, CDTensor2 const & rhs ) {
   lhs -= rhs;
   return lhs;
}


inline CDTensor2 operator*( CDTensor2 lhs, double f ) {
   lhs *= f;
   return lhs;
}


inline CDTensor2 operator/( CDTensor2 lhs, double f ) {
   lhs /= f;
   return lhs;
}


inline CDTensor2 operator*( double f, CDTensor2 lhs ) {
   lhs *= f;
   return lhs;
}


inline CDTensor2 operator/( double f, CDTensor2 lhs ) {
   lhs /= f;
   return lhs;
}



/*****************************************************************************
***   CDTensor3   ************************************************************
*****************************************************************************/

class CDTensor3 {
public:

	CDTensor3() {
		m_iDimensions[0] = 0;
		m_iDimensions[1] = 0;
		m_iDimensions[2] = 0;
	}

	CDTensor3(int i, int j, int k) {
		m_iDimensions[0] = i;
		m_iDimensions[1] = j;
		m_iDimensions[2] = k;
		m_pData.resize(i*j*k,0);
	}

	~CDTensor3() {
	}

	CDTensor3(const CDTensor3 &t) : m_pData(t.m_pData) {
		m_iDimensions[0] = t.m_iDimensions[0];
		m_iDimensions[1] = t.m_iDimensions[1];
		m_iDimensions[2] = t.m_iDimensions[2];
	}

	CDTensor3& operator=(const CDTensor3& t) {
		m_pData = t.m_pData;
		m_iDimensions[0] = t.m_iDimensions[0];
		m_iDimensions[1] = t.m_iDimensions[1];
		m_iDimensions[2] = t.m_iDimensions[2];
		return *this;
	}

	double& operator[] (int i) {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[0]*m_iDimensions[1]*m_iDimensions[2])) {
				eprintf("CDTensor3 &operator [] boundary error (%d/%d).\n",i,m_iDimensions[0]*m_iDimensions[1]*m_iDimensions[2]);
				abort();
			}
		#endif
		return m_pData[i];
	}

	double operator[] (int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[0]*m_iDimensions[1]*m_iDimensions[2])) {
				eprintf("CDTensor3 &operator [] boundary error (%d/%d).\n",i,m_iDimensions[0]*m_iDimensions[1]*m_iDimensions[2]);
				abort();
			}
		#endif
		return m_pData[i];
	}

	double& operator() (int i, int j, int k) {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (k < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[1]) || (k >= m_iDimensions[2])) {
				eprintf("CDTensor3 &operator () boundary error (%d/%d, %d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[1],k,m_iDimensions[2]);
				abort();
			}
		#endif
		return m_pData[i*m_iDimensions[1]*m_iDimensions[2]+j*m_iDimensions[2]+k];
	}

	double operator() (int i, int j, int k) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (k < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[1]) || (k >= m_iDimensions[2])) {
				eprintf("CDTensor3 operator () boundary error (%d/%d, %d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[1],k,m_iDimensions[2]);
				abort();
			}
		#endif
		return m_pData[i*m_iDimensions[1]*m_iDimensions[2]+j*m_iDimensions[2]+k];
	}

	CDTensor3& operator += (const CDTensor3 &t) {
		#ifdef TENSOR_BOUND_CHECK
			if (!DimensionsMatch(t)) {
				eprintf("CDTensor3 &operator += dimension mismatch.\n");
				abort();
			}
		#endif
		std::transform(m_pData.begin(),m_pData.end(),t.m_pData.begin(),m_pData.begin(),std::plus<double>());
		return *this;
	}

	CDTensor3& operator -= (const CDTensor3 &t) {
		#ifdef TENSOR_BOUND_CHECK
			if (!DimensionsMatch(t)) {
				eprintf("CDTensor3 &operator -= dimension mismatch.\n");
				abort();
			}
		#endif
		std::transform(m_pData.begin(),m_pData.end(),t.m_pData.begin(),m_pData.begin(),std::minus<double>());
		return *this;
	}

	CDTensor3& operator *= (double f) {
		//std::transform(m_pData.begin(),m_pData.end(),m_pData.begin(),std::bind1st(std::multiplies<double>(), f));
		for (unsigned long z=0;z<m_pData.size();z++)
			m_pData[z] *= f;
		return *this;
	}

	CDTensor3& operator /= (double f) {
		//std::transform(m_pData.begin(),m_pData.end(),m_pData.begin(),std::bind1st(std::divides<double>(), f));
		for (unsigned long z=0;z<m_pData.size();z++)
			m_pData[z] /= f;
		return *this;
	}

	CDTensor2 Reduce1(int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[0])) {
				eprintf("CDTensor3 Reduce1() boundary error (%d/%d).\n",i,m_iDimensions[0]);
				abort();
			}
		#endif
		CDTensor2 t(m_iDimensions[1],m_iDimensions[2]);
		for (int z=0;z<m_iDimensions[1];z++)
			for (int z2=0;z2<m_iDimensions[2];z2++)
				t(z,z2) = (*this)(i,z,z2);
		return t;
	}

	CDTensor2 Reduce2(int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[1])) {
				eprintf("CDTensor3 Reduce2() boundary error (%d/%d).\n",i,m_iDimensions[1]);
				abort();
			}
		#endif
		CDTensor2 t(m_iDimensions[0],m_iDimensions[2]);
		for (int z=0;z<m_iDimensions[0];z++)
			for (int z2=0;z2<m_iDimensions[2];z2++)
				t(z,z2) = (*this)(z,i,z2);
		return t;
	}

	CDTensor2 Reduce3(int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[2])) {
				eprintf("CDTensor3 Reduce3() boundary error (%d/%d).\n",i,m_iDimensions[2]);
				abort();
			}
		#endif
		CDTensor2 t(m_iDimensions[0],m_iDimensions[1]);
		for (int z=0;z<m_iDimensions[0];z++)
			for (int z2=0;z2<m_iDimensions[1];z2++)
				t(z,z2) = (*this)(z,z2,i);
		return t;
	}

	CDTensor1 Reduce12(int i, int j) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[1])) {
				eprintf("CDTensor3 Reduce12() boundary error (%d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[1]);
				abort();
			}
		#endif
		CDTensor1 t(m_iDimensions[2]);
		for (int z=0;z<m_iDimensions[2];z++)
			t(z) = (*this)(i,j,z);
		return t;
	}

	CDTensor1 Reduce13(int i, int j) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[2])) {
				eprintf("CDTensor3 Reduce13() boundary error (%d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[2]);
				abort();
			}
		#endif
		CDTensor1 t(m_iDimensions[1]);
		for (int z=0;z<m_iDimensions[1];z++)
			t(z) = (*this)(i,z,j);
		return t;
	}

	CDTensor1 Reduce23(int i, int j) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (i >= m_iDimensions[1]) || (j >= m_iDimensions[2])) {
				eprintf("CDTensor3 Reduce23() boundary error (%d/%d, %d/%d).\n",i,m_iDimensions[1],j,m_iDimensions[2]);
				abort();
			}
		#endif
		CDTensor1 t(m_iDimensions[0]);
		for (int z=0;z<m_iDimensions[0];z++)
			t(z) = (*this)(z,i,j);
		return t;
	}

	void setZero() {
		std::fill(m_pData.begin(),m_pData.end(),0);
	}

	int getDimension(int i) const {
		return m_iDimensions[i];
	}

	bool DimensionsMatch(const CDTensor3 &t) const {
		if ((m_iDimensions[0] != t.m_iDimensions[0]) ||
			(m_iDimensions[1] != t.m_iDimensions[1]) ||
			(m_iDimensions[2] != t.m_iDimensions[2])) {
			eprintf("Dimension mismatch ( %d, %d, %d ) <--> ( %d, %d, %d ).\n",m_iDimensions[0],m_iDimensions[1],m_iDimensions[2],t.m_iDimensions[0],t.m_iDimensions[1],t.m_iDimensions[2]);
			return false;
		}
		return true;
	}

private:
	std::vector<double> m_pData;
	int m_iDimensions[3];
};


inline CDTensor3 operator+( CDTensor3 lhs, CDTensor3 const & rhs ) {
   lhs += rhs;
   return lhs;
}


inline CDTensor3 operator-( CDTensor3 lhs, CDTensor3 const & rhs ) {
   lhs -= rhs;
   return lhs;
}


inline CDTensor3 operator*( CDTensor3 lhs, double f ) {
   lhs *= f;
   return lhs;
}


inline CDTensor3 operator/( CDTensor3 lhs, double f ) {
   lhs /= f;
   return lhs;
}


inline CDTensor3 operator*( double f, CDTensor3 lhs ) {
   lhs *= f;
   return lhs;
}


inline CDTensor3 operator/( double f, CDTensor3 lhs ) {
   lhs /= f;
   return lhs;
}



/*****************************************************************************
***   CDTensor4   ************************************************************
*****************************************************************************/

class CDTensor4 {
public:

	CDTensor4() {
		m_iDimensions[0] = 0;
		m_iDimensions[1] = 0;
		m_iDimensions[2] = 0;
		m_iDimensions[3] = 0;
	}

	CDTensor4(int i, int j, int k, int l) {
		m_iDimensions[0] = i;
		m_iDimensions[1] = j;
		m_iDimensions[2] = k;
		m_iDimensions[3] = l;
		m_pData.resize(i*j*k*l,0);
	}

	~CDTensor4() {
	}

	CDTensor4(const CDTensor4 &t) : m_pData(t.m_pData) {
		m_iDimensions[0] = t.m_iDimensions[0];
		m_iDimensions[1] = t.m_iDimensions[1];
		m_iDimensions[2] = t.m_iDimensions[2];
		m_iDimensions[3] = t.m_iDimensions[3];
	}

	CDTensor4& operator=(const CDTensor4& t) {
		m_pData = t.m_pData;
		m_iDimensions[0] = t.m_iDimensions[0];
		m_iDimensions[1] = t.m_iDimensions[1];
		m_iDimensions[2] = t.m_iDimensions[2];
		m_iDimensions[3] = t.m_iDimensions[3];
		return *this;
	}

	double& operator[] (int i) {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[0]*m_iDimensions[1]*m_iDimensions[2]*m_iDimensions[3])) {
				eprintf("CDTensor4 &operator [] boundary error (%d/%d).\n",i,m_iDimensions[0]*m_iDimensions[1]*m_iDimensions[2]*m_iDimensions[3]);
				abort();
			}
		#endif
		return m_pData[i];
	}

	double operator[] (int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[0]*m_iDimensions[1]*m_iDimensions[2]*m_iDimensions[3])) {
				eprintf("CDTensor4 &operator [] boundary error (%d/%d).\n",i,m_iDimensions[0]*m_iDimensions[1]*m_iDimensions[2]*m_iDimensions[3]);
				abort();
			}
		#endif
		return m_pData[i];
	}

	double& operator() (int i, int j, int k, int l) {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (k < 0) || (l < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[1]) || (k >= m_iDimensions[2]) || (l >= m_iDimensions[3])) {
				eprintf("CDTensor4 &operator () boundary error (%d/%d, %d/%d, %d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[1],k,m_iDimensions[2],l,m_iDimensions[3]);
				abort();
			}
		#endif
		return m_pData[i*m_iDimensions[1]*m_iDimensions[2]*m_iDimensions[3]+j*m_iDimensions[2]*m_iDimensions[3]+k*m_iDimensions[3]+l];
	}

	double operator() (int i, int j, int k, int l) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (k < 0) || (l < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[1]) || (k >= m_iDimensions[2]) || (l >= m_iDimensions[3])) {
				eprintf("CDTensor4 operator () boundary error (%d/%d, %d/%d, %d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[1],k,m_iDimensions[2],l,m_iDimensions[3]);
				abort();
			}
		#endif
		return m_pData[i*m_iDimensions[1]*m_iDimensions[2]*m_iDimensions[3]+j*m_iDimensions[2]*m_iDimensions[3]+k*m_iDimensions[3]+l];
	}

	CDTensor4& operator += (const CDTensor4 &t) {
		#ifdef TENSOR_BOUND_CHECK
			if (!DimensionsMatch(t)) {
				eprintf("CDTensor4 &operator += dimension mismatch.\n");
				abort();
			}
		#endif
		std::transform(m_pData.begin(),m_pData.end(),t.m_pData.begin(),m_pData.begin(),std::plus<double>());
		return *this;
	}

	CDTensor4& operator -= (const CDTensor4 &t) {
		#ifdef TENSOR_BOUND_CHECK
			if (!DimensionsMatch(t)) {
				eprintf("CDTensor4 &operator -= dimension mismatch.\n");
				abort();
			}
		#endif
		std::transform(m_pData.begin(),m_pData.end(),t.m_pData.begin(),m_pData.begin(),std::minus<double>());
		return *this;
	}

	CDTensor4& operator *= (double f) {
		//std::transform(m_pData.begin(),m_pData.end(),m_pData.begin(),std::bind1st(std::multiplies<double>(), f));
		for (unsigned long z=0;z<m_pData.size();z++)
			m_pData[z] *= f;
		return *this;
	}

	CDTensor4& operator /= (double f) {
		//std::transform(m_pData.begin(),m_pData.end(),m_pData.begin(),std::bind1st(std::divides<double>(), f));
		for (unsigned long z=0;z<m_pData.size();z++)
			m_pData[z] /= f;
		return *this;
	}

	CDTensor3 Reduce1(int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[0])) {
				eprintf("CDTensor4 Reduce1() boundary error (%d/%d).\n",i,m_iDimensions[0]);
				abort();
			}
		#endif
		CDTensor3 t(m_iDimensions[1],m_iDimensions[2],m_iDimensions[3]);
		for (int z=0;z<m_iDimensions[1];z++)
			for (int z2=0;z2<m_iDimensions[2];z2++)
				for (int z3=0;z3<m_iDimensions[3];z3++)
					t(z,z2,z3) = (*this)(i,z,z2,z3);
		return t;
	}

	CDTensor3 Reduce2(int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[1])) {
				eprintf("CDTensor4 Reduce2() boundary error (%d/%d).\n",i,m_iDimensions[1]);
				abort();
			}
		#endif
		CDTensor3 t(m_iDimensions[0],m_iDimensions[2],m_iDimensions[3]);
		for (int z=0;z<m_iDimensions[0];z++)
			for (int z2=0;z2<m_iDimensions[2];z2++)
				for (int z3=0;z3<m_iDimensions[3];z3++)
					t(z,z2,z3) = (*this)(z,i,z2,z3);
		return t;
	}

	CDTensor3 Reduce3(int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[2])) {
				eprintf("CDTensor4 Reduce3() boundary error (%d/%d).\n",i,m_iDimensions[2]);
				abort();
			}
		#endif
		CDTensor3 t(m_iDimensions[0],m_iDimensions[1],m_iDimensions[3]);
		for (int z=0;z<m_iDimensions[0];z++)
			for (int z2=0;z2<m_iDimensions[1];z2++)
				for (int z3=0;z3<m_iDimensions[3];z3++)
					t(z,z2,z3) = (*this)(z,z2,i,z3);
		return t;
	}

	CDTensor3 Reduce4(int i) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (i >= m_iDimensions[3])) {
				eprintf("CDTensor4 Reduce4() boundary error (%d/%d).\n",i,m_iDimensions[3]);
				abort();
			}
		#endif
		CDTensor3 t(m_iDimensions[0],m_iDimensions[1],m_iDimensions[2]);
		for (int z=0;z<m_iDimensions[0];z++)
			for (int z2=0;z2<m_iDimensions[1];z2++)
				for (int z3=0;z3<m_iDimensions[2];z3++)
					t(z,z2,z3) = (*this)(z,z2,z3,i);
		return t;
	}

	CDTensor2 Reduce12(int i, int j) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[1])) {
				eprintf("CDTensor4 Reduce12() boundary error (%d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[1]);
				abort();
			}
		#endif
		CDTensor2 t(m_iDimensions[2],m_iDimensions[3]);
		for (int z=0;z<m_iDimensions[2];z++)
			for (int z2=0;z2<m_iDimensions[3];z2++)
				t(z,z2) = (*this)(i,j,z,z2);
		return t;
	}

	CDTensor2 Reduce13(int i, int j) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[2])) {
				eprintf("CDTensor4 Reduce13() boundary error (%d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[2]);
				abort();
			}
		#endif
		CDTensor2 t(m_iDimensions[1],m_iDimensions[3]);
		for (int z=0;z<m_iDimensions[1];z++)
			for (int z2=0;z2<m_iDimensions[3];z2++)
				t(z,z2) = (*this)(i,z,j,z2);
		return t;
	}

	CDTensor2 Reduce14(int i, int j) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[3])) {
				eprintf("CDTensor4 Reduce14() boundary error (%d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[3]);
				abort();
			}
		#endif
		CDTensor2 t(m_iDimensions[1],m_iDimensions[2]);
		for (int z=0;z<m_iDimensions[1];z++)
			for (int z2=0;z2<m_iDimensions[2];z2++)
				t(z,z2) = (*this)(i,z,z2,j);
		return t;
	}

	CDTensor2 Reduce23(int i, int j) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (i >= m_iDimensions[1]) || (j >= m_iDimensions[2])) {
				eprintf("CDTensor4 Reduce23() boundary error (%d/%d, %d/%d).\n",i,m_iDimensions[1],j,m_iDimensions[2]);
				abort();
			}
		#endif
		CDTensor2 t(m_iDimensions[0],m_iDimensions[3]);
		for (int z=0;z<m_iDimensions[0];z++)
			for (int z2=0;z2<m_iDimensions[3];z2++)
				t(z,z2) = (*this)(z,i,j,z2);
		return t;
	}

	CDTensor2 Reduce24(int i, int j) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (i >= m_iDimensions[1]) || (j >= m_iDimensions[3])) {
				eprintf("CDTensor4 Reduce24() boundary error (%d/%d, %d/%d).\n",i,m_iDimensions[1],j,m_iDimensions[3]);
				abort();
			}
		#endif
		CDTensor2 t(m_iDimensions[0],m_iDimensions[2]);
		for (int z=0;z<m_iDimensions[0];z++)
			for (int z2=0;z2<m_iDimensions[2];z2++)
				t(z,z2) = (*this)(z,i,z2,j);
		return t;
	}

	CDTensor2 Reduce34(int i, int j) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (i >= m_iDimensions[2]) || (j >= m_iDimensions[3])) {
				eprintf("CDTensor4 Reduce34() boundary error (%d/%d, %d/%d).\n",i,m_iDimensions[2],j,m_iDimensions[3]);
				abort();
			}
		#endif
		CDTensor2 t(m_iDimensions[0],m_iDimensions[1]);
		for (int z=0;z<m_iDimensions[0];z++)
			for (int z2=0;z2<m_iDimensions[1];z2++)
				t(z,z2) = (*this)(z,z2,i,j);
		return t;
	}

	CDTensor1 Reduce123(int i, int j, int k) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (k < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[1]) || (k >= m_iDimensions[2])) {
				eprintf("CDTensor4 Reduce123() boundary error (%d/%d, %d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[1],k,m_iDimensions[2]);
				abort();
			}
		#endif
		CDTensor1 t(m_iDimensions[3]);
		for (int z=0;z<m_iDimensions[3];z++)
			t(z) = (*this)(i,j,k,z);
		return t;
	}

	CDTensor1 Reduce124(int i, int j, int k) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (k < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[1]) || (k >= m_iDimensions[3])) {
				eprintf("CDTensor4 Reduce124() boundary error (%d/%d, %d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[1],k,m_iDimensions[3]);
				abort();
			}
		#endif
		CDTensor1 t(m_iDimensions[2]);
		for (int z=0;z<m_iDimensions[2];z++)
			t(z) = (*this)(i,j,z,k);
		return t;
	}

	CDTensor1 Reduce134(int i, int j, int k) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (k < 0) || (i >= m_iDimensions[0]) || (j >= m_iDimensions[2]) || (k >= m_iDimensions[3])) {
				eprintf("CDTensor4 Reduce134() boundary error (%d/%d, %d/%d, %d/%d).\n",i,m_iDimensions[0],j,m_iDimensions[2],k,m_iDimensions[3]);
				abort();
			}
		#endif
		CDTensor1 t(m_iDimensions[1]);
		for (int z=0;z<m_iDimensions[1];z++)
			t(z) = (*this)(i,z,j,k);
		return t;
	}

	CDTensor1 Reduce234(int i, int j, int k) const {
		#ifdef TENSOR_BOUND_CHECK
			if ((i < 0) || (j < 0) || (k < 0) || (i >= m_iDimensions[1]) || (j >= m_iDimensions[2]) || (k >= m_iDimensions[3])) {
				eprintf("CDTensor4 Reduce234() boundary error (%d/%d, %d/%d, %d/%d).\n",i,m_iDimensions[1],j,m_iDimensions[2],k,m_iDimensions[3]);
				abort();
			}
		#endif
		CDTensor1 t(m_iDimensions[0]);
		for (int z=0;z<m_iDimensions[0];z++)
			t(z) = (*this)(z,i,j,k);
		return t;
	}

	void setZero() {
		std::fill(m_pData.begin(),m_pData.end(),0);
	}

	int getDimension(int i) const {
		return m_iDimensions[i];
	}

	bool DimensionsMatch(const CDTensor4 &t) const {
		if ((m_iDimensions[0] != t.m_iDimensions[0]) ||
			(m_iDimensions[1] != t.m_iDimensions[1]) ||
			(m_iDimensions[2] != t.m_iDimensions[2]) ||
			(m_iDimensions[3] != t.m_iDimensions[3])) {
			eprintf("Dimension mismatch ( %d, %d, %d, %d ) <--> ( %d, %d, %d, %d ).\n",m_iDimensions[0],m_iDimensions[1],m_iDimensions[2],m_iDimensions[3],t.m_iDimensions[0],t.m_iDimensions[1],t.m_iDimensions[2],t.m_iDimensions[3]);
			return false;
		}
		return true;
	}

private:
	std::vector<double> m_pData;
	int m_iDimensions[4];
};


inline CDTensor4 operator+( CDTensor4 lhs, CDTensor4 const & rhs ) {
   lhs += rhs;
   return lhs;
}


inline CDTensor4 operator-( CDTensor4 lhs, CDTensor4 const & rhs ) {
   lhs -= rhs;
   return lhs;
}


inline CDTensor4 operator*( CDTensor4 lhs, double f ) {
   lhs *= f;
   return lhs;
}


inline CDTensor4 operator/( CDTensor4 lhs, double f ) {
   lhs /= f;
   return lhs;
}


inline CDTensor4 operator*( double f, CDTensor4 lhs ) {
   lhs *= f;
   return lhs;
}


inline CDTensor4 operator/( double f, CDTensor4 lhs ) {
   lhs /= f;
   return lhs;
}



#endif




