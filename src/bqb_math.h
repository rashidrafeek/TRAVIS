/***********************************************************************************

    LibBQB - File Format and Compression Algorithms for Trajectories of
             Volumetric Data and Atom Positions

    https://brehm-research.de/bqb

    Free software, licensed under GNU LGPL v3

    Copyright (c) Martin Brehm and Martin Thomas,
                  Martin Luther University Halle-Wittenberg, Germany,
                  2016 - 2021.

    Please cite:  M. Brehm, M. Thomas: "An Efficient Lossless Compression Algorithm
                  for Trajectories of Atom Positions and Volumetric Data",
                  J. Chem. Inf. Model. 2018, 58 (10), pp 2092-2107.

    --------------------------------------------------------------------------------

    LibBQB is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

***********************************************************************************/



#ifndef BQB_MATH_H
#define BQB_MATH_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"



class CBQBInterface;



#ifdef BQB_INSIDE_TRAVIS


#include "xdvector3.h"
#include "xdvectorn.h"
#include "xdmatrix3.h"
#include "xdmatrixmn.h"

#define bqbpow2 pow2
#define bqbpow3 pow3
#define bqbpow4 pow4
#define bqbpow5 pow5

#define CBQBDVector3  CxDVector3
#define CBQBDVectorN  CxDVectorN
#define CBQBDMatrix3  CxDMatrix3
#define CBQBDMatrixMN CxDMatrixMN


#else


class CBQBDMatrix3;
class CBQBDVector3;
class CBQBDMatrixMN;
class CBQBDVectorN;


inline double bqbpow2(double x) { return x*x; }
inline double bqbpow3(double x) { return x*x*x; }
inline double bqbpow4(double x) { return (x*x)*(x*x); }
inline double bqbpow5(double x) { return (x*x*x)*(x*x); }


#endif



class CBQBMath { // Encapsulates everything that may print to log
public:

	explicit CBQBMath(CBQBInterface &i) : m_IF(i) {
	}

	void DumpMatrix(const CBQBDMatrix3 &m) const;

	void DumpMatrixSmall(const CBQBDMatrix3 &m) const;

	void DumpMatrix(const CBQBDMatrixMN &m) const;

	void DumpMatrixSmall(const CBQBDMatrixMN &m) const;

	void InvertMatrix(CBQBDMatrix3 &m);

	double Angle(const CBQBDVector3 &vec1, const CBQBDVector3 &vec2) const;

	double Angle_Deg(const CBQBDVector3 &vec1, const CBQBDVector3 &vec2) const;

	double DotP(const CBQBDVectorN &vec1, const CBQBDVectorN &vec2) const;

	double Angle(const CBQBDVectorN &vec1, const CBQBDVectorN &vec2) const;

	double Angle_Deg(const CBQBDVectorN &vec1, const CBQBDVectorN &vec2) const;

private:
	CBQBInterface &m_IF;
};



#ifndef BQB_INSIDE_TRAVIS


//double Dihedral(const CBQBDVector3 &vec1, const CBQBDVector3 &vec2, const CBQBDVector3 &norm, bool absolute);



class CBQBDVector3
{
public:
	CBQBDVector3() {
	}


	~CBQBDVector3() {
	}
	

	CBQBDVector3(const CBQBDVector3 &v) {
		m_pData[0] = v.m_pData[0];
		m_pData[1] = v.m_pData[1];
		m_pData[2] = v.m_pData[2];
	}


	explicit CBQBDVector3(double f) {
		m_pData[0] = f;
		m_pData[1] = f;
		m_pData[2] = f;
	}


	CBQBDVector3(double x, double y, double z) {
		m_pData[0] = x;
		m_pData[1] = y;
		m_pData[2] = z;
	}


	CBQBDVector3 & operator = (const CBQBDVector3 &v) {
		m_pData[0] = v.m_pData[0];
		m_pData[1] = v.m_pData[1];
		m_pData[2] = v.m_pData[2];
		return *this;
	}


	double &GetAt(int i) {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}


	double &operator [] (int i) {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	double GetAt(int i) const {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::GetAt(int): %d...",i);
		#endif
		#ifdef DEBUG_ARRAYS
			if (i > 2) {
				printf("CBQBDVector3::GetAt(int): Boundary Error (%d/3)...",i);
				abort();
			}
		#endif
		return m_pData[i];
	}


	double operator [] (int i) const {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	CBQBDVector3 operator + (const CBQBDVector3 &v) const {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::operator + (CBQBDVector3&)\n");
		#endif
		return CBQBDVector3(m_pData[0]+v.m_pData[0],m_pData[1]+v.m_pData[1],m_pData[2]+v.m_pData[2]);
	}


	CBQBDVector3 operator - (const CBQBDVector3 &v) const {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::operator - (CBQBDVector3&)\n");
		#endif
		return CBQBDVector3(m_pData[0]-v.m_pData[0],m_pData[1]-v.m_pData[1],m_pData[2]-v.m_pData[2]);
	}


	CBQBDVector3 operator * (double f) const {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::operator * (double)\n");
		#endif
		return CBQBDVector3(m_pData[0]*f,m_pData[1]*f,m_pData[2]*f);
	}


	CBQBDVector3 operator / (double f) const {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::operator / (double)\n");
		#endif
		return CBQBDVector3(m_pData[0]/f,m_pData[1]/f,m_pData[2]/f);
	}


	void operator += (const CBQBDVector3 &v) {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::operator += (CBQBDVector3&)\n");
		#endif
		m_pData[0] += v.m_pData[0];
		m_pData[1] += v.m_pData[1];
		m_pData[2] += v.m_pData[2];
	}


	void operator -= (const CBQBDVector3 &v) {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::operator -= (CBQBDVector3&)\n");
		#endif
		m_pData[0] -= v.m_pData[0];
		m_pData[1] -= v.m_pData[1];
		m_pData[2] -= v.m_pData[2];
	}


	void operator *= (double f) {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::operator *= (double)\n");
		#endif
		m_pData[0] *= f;
		m_pData[1] *= f;
		m_pData[2] *= f;
	}


	void operator /= (double f) {
		#ifdef DEBUG_CDVECTOR3
			printf("@ CBQBDVector3::operator /= (double)\n");
		#endif
		m_pData[0] /= f;
		m_pData[1] /= f;
		m_pData[2] /= f;
	}


	double GetLength() const {
		return sqrt(m_pData[0]*m_pData[0]+m_pData[1]*m_pData[1]+m_pData[2]*m_pData[2]);
	}


	double GetLengthSqr() const {
		return m_pData[0]*m_pData[0]+m_pData[1]*m_pData[1]+m_pData[2]*m_pData[2];
	}


	void Normalize() {
		double l;
		l = GetLength();
		m_pData[0] /= l;
		m_pData[1] /= l;
		m_pData[2] /= l;
	}


	void Chop(double d) {
		if ((m_pData[0] != 0) && (fabs(m_pData[0]) < d))
			m_pData[0] = 0;
		if ((m_pData[1] != 0) && (fabs(m_pData[1]) < d))
			m_pData[1] = 0;
		if ((m_pData[2] != 0) && (fabs(m_pData[2]) < d))
			m_pData[2] = 0;
	}


	void PointRoot(const CBQBDVector3 &vec1, const CBQBDVector3 &vec2, const CBQBDVector3 &point);


private:
	double m_pData[3];
};



inline void Swap(CBQBDVector3 &vec1, CBQBDVector3 &vec2) {
	CBQBDVector3 t;
	t = vec1;
	vec1 = vec2;
	vec2 = t;
}


inline CBQBDVector3 operator * (double f, const CBQBDVector3 &v) {
	return v*f;
}


inline CBQBDVector3 CrossP(const CBQBDVector3 &vec1, const CBQBDVector3 &vec2) {
	return CBQBDVector3(vec1[1]*vec2[2] - vec1[2]*vec2[1],
		vec1[2]*vec2[0] - vec1[0]*vec2[2],
		vec1[0]*vec2[1] - vec1[1]*vec2[0]);
}


inline double DotP(const CBQBDVector3 &vec1, const CBQBDVector3 &vec2) {
	return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
}


inline double VecDist(const CBQBDVector3 &vec1, const CBQBDVector3 &vec2) {
	return sqrt((vec1[0]-vec2[0])*(vec1[0]-vec2[0]) + (vec1[1]-vec2[1])*(vec1[1]-vec2[1]) + (vec1[2]-vec2[2])*(vec1[2]-vec2[2]));
}


inline CBQBDVector3 Normalize(const CBQBDVector3 &vec) {
	double tf;
	tf = vec.GetLength();
	return CBQBDVector3(vec[0]/tf,vec[1]/tf,vec[2]/tf);
}



class CBQBDMatrix3 {
public:

	CBQBDMatrix3() {
	}

	~CBQBDMatrix3() {
	}


	CBQBDMatrix3(const CBQBDMatrix3 &m) {
		m_pData[0] = m.m_pData[0];
		m_pData[1] = m.m_pData[1];
		m_pData[2] = m.m_pData[2];
		m_pData[3] = m.m_pData[3];
		m_pData[4] = m.m_pData[4];
		m_pData[5] = m.m_pData[5];
		m_pData[6] = m.m_pData[6];
		m_pData[7] = m.m_pData[7];
		m_pData[8] = m.m_pData[8];
	}


	explicit CBQBDMatrix3(double f) {
		m_pData[0] = f;
		m_pData[1] = f;
		m_pData[2] = f;
		m_pData[3] = f;
		m_pData[4] = f;
		m_pData[5] = f;
		m_pData[6] = f;
		m_pData[7] = f;
		m_pData[8] = f;
	}


	CBQBDMatrix3(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
		m_pData[0] = a;
		m_pData[1] = b;
		m_pData[2] = c;
		m_pData[3] = d;
		m_pData[4] = e;
		m_pData[5] = f;
		m_pData[6] = g;
		m_pData[7] = h;
		m_pData[8] = i;
	}


	CBQBDMatrix3 & operator = (const CBQBDMatrix3 &v) {
		m_pData[0] = v.m_pData[0];
		m_pData[1] = v.m_pData[1];
		m_pData[2] = v.m_pData[2];
		m_pData[3] = v.m_pData[3];
		m_pData[4] = v.m_pData[4];
		m_pData[5] = v.m_pData[5];
		m_pData[6] = v.m_pData[6];
		m_pData[7] = v.m_pData[7];
		m_pData[8] = v.m_pData[8];
		return *this;
	}


	double &GetAt(int i) {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}


	double &GetAt(int i, int j) {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::GetAt(int, int): %d, %d...",i,j);
		#endif
		return m_pData[i*3+j];
	}


	double &operator [] (int i) {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	double GetAt(int i) const {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}


	double GetAt(int i, int j) const {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::GetAt(int, int): %d, %d...",i,j);
		#endif
		return m_pData[i*3+j];
	}


	double operator [] (int i) const {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	double &operator () (int i, int j) {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::operator () (int, int): %d, %d\n",i,j);
		#endif
		return GetAt(i,j);
	}


	double operator () (int i, int j) const {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::operator () (int, int): %d, %d\n",i,j);
		#endif
		return GetAt(i,j);
	}


	void operator *= (double f) {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::operator *= (double)\n");
		#endif
		m_pData[0] *= f;
		m_pData[1] *= f;
		m_pData[2] *= f;
		m_pData[3] *= f;
		m_pData[4] *= f;
		m_pData[5] *= f;
		m_pData[6] *= f;
		m_pData[7] *= f;
		m_pData[8] *= f;
	}


	CBQBDMatrix3 operator * (const CBQBDMatrix3 &m) const {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::operator * (CBQBDMatrix3)\n");
		#endif
		return CBQBDMatrix3(
			GetAt(0,0)*m.GetAt(0,0) + GetAt(1,0)*m.GetAt(0,1) + GetAt(2,0)*m.GetAt(0,2),
			GetAt(0,1)*m.GetAt(0,0) + GetAt(1,1)*m.GetAt(0,1) + GetAt(2,1)*m.GetAt(0,2),
			GetAt(0,2)*m.GetAt(0,0) + GetAt(1,2)*m.GetAt(0,1) + GetAt(2,2)*m.GetAt(0,2),
			GetAt(0,0)*m.GetAt(1,0) + GetAt(1,0)*m.GetAt(1,1) + GetAt(2,0)*m.GetAt(1,2),
			GetAt(0,1)*m.GetAt(1,0) + GetAt(1,1)*m.GetAt(1,1) + GetAt(2,1)*m.GetAt(1,2),
			GetAt(0,2)*m.GetAt(1,0) + GetAt(1,2)*m.GetAt(1,1) + GetAt(2,2)*m.GetAt(1,2),
			GetAt(0,0)*m.GetAt(2,0) + GetAt(1,0)*m.GetAt(2,1) + GetAt(2,0)*m.GetAt(2,2),
			GetAt(0,1)*m.GetAt(2,0) + GetAt(1,1)*m.GetAt(2,1) + GetAt(2,1)*m.GetAt(2,2),
			GetAt(0,2)*m.GetAt(2,0) + GetAt(1,2)*m.GetAt(2,1) + GetAt(2,2)*m.GetAt(2,2)
			);
	}


	CBQBDVector3 operator * (const CBQBDVector3 &v) const {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::operator * (CBQBDVector3)\n");
		#endif
		return CBQBDVector3(
			GetAt(0,0)*v[0] + GetAt(1,0)*v[1] + GetAt(2,0)*v[2],
			GetAt(0,1)*v[0] + GetAt(1,1)*v[1] + GetAt(2,1)*v[2],
			GetAt(0,2)*v[0] + GetAt(1,2)*v[1] + GetAt(2,2)*v[2]);
	}


	CBQBDMatrix3 operator * (double f) const {
		#ifdef DEBUG_CMATRIX3
			printf("@ CBQBDMatrix3::operator * (double)\n");
		#endif
		return CBQBDMatrix3(m_pData[0]*f,m_pData[1]*f,m_pData[2]*f,m_pData[3]*f,m_pData[4]*f,m_pData[5]*f,m_pData[6]*f,m_pData[7]*f,m_pData[8]*f);
	}


	void RotMat(const CBQBDVector3 &vec, double a) {
		double ca, sa;
		ca = cos(a);
		sa = sin(a);
		GetAt(0,0) = (vec[0]*vec[0]*(1.0-ca) + ca);
		GetAt(1,0) = (vec[0]*vec[1]*(1.0-ca) - vec[2]*sa);
		GetAt(2,0) = (vec[0]*vec[2]*(1.0-ca) + vec[1]*sa);
		GetAt(0,1) = (vec[1]*vec[0]*(1.0-ca) + vec[2]*sa);
		GetAt(1,1) = (vec[1]*vec[1]*(1.0-ca) + ca);
		GetAt(2,1) = (vec[1]*vec[2]*(1.0-ca) - vec[0]*sa);
		GetAt(0,2) = (vec[2]*vec[0]*(1.0-ca) - vec[1]*sa);
		GetAt(1,2) = (vec[2]*vec[1]*(1.0-ca) + vec[0]*sa);
		GetAt(2,2) = (vec[2]*vec[2]*(1.0-ca) + ca);
	}


	void Unity() { 
		GetAt(0,0) = 1;
		GetAt(1,0) = 0;
		GetAt(2,0) = 0;
		GetAt(0,1) = 0;
		GetAt(1,1) = 1;
		GetAt(2,1) = 0;
		GetAt(0,2) = 0;
		GetAt(1,2) = 0;
		GetAt(2,2) = 1;
	}


	CBQBDMatrix3 Transpose() {
		return CBQBDMatrix3(
			GetAt(0,0), GetAt(1,0), GetAt(2,0),
			GetAt(0,1), GetAt(1,1), GetAt(2,1),
			GetAt(0,2), GetAt(1,2), GetAt(2,2));
	}


	double Det() {
		#define _deta11 GetAt(0,0)
		#define _deta12 GetAt(0,1)
		#define _deta13 GetAt(0,2)
		#define _deta21 GetAt(1,0)
		#define _deta22 GetAt(1,1)
		#define _deta23 GetAt(1,2)
		#define _deta31 GetAt(2,0)
		#define _deta32 GetAt(2,1)
		#define _deta33 GetAt(2,2)

		return  _deta11 * (_deta33*_deta22 - _deta32*_deta23)
			  - _deta21 * (_deta33*_deta12 - _deta32*_deta13)
			  + _deta31 * (_deta23*_deta12 - _deta22*_deta13);
	}


private:
	double m_pData[9];
};



class CBQBDVectorN {
public:

	CBQBDVectorN() {
		m_iDim = 0;
		m_pData = NULL;
	}


	~CBQBDVectorN() {
		if (m_pData != NULL) {
			delete[] m_pData;
			m_pData = NULL;
		}
	}


	CBQBDVectorN(const CBQBDVectorN &v) {
		m_iDim = v.m_iDim;
		m_pData = new double[m_iDim];
		memcpy(m_pData,v.m_pData,sizeof(double)*m_iDim);
	}


	CBQBDVectorN & operator = (const CBQBDVectorN &v) {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::operator = (const CBQBDVectorN &)\n");
		#endif
		if (m_pData != NULL)
			delete[] m_pData;
		m_iDim = v.m_iDim;
		m_pData = new double[m_iDim];
		memcpy(m_pData,v.m_pData,sizeof(double)*m_iDim);
		return *this;
	}


	CBQBDVectorN(int i, double f) {
		int z;
		m_iDim = i;
		m_pData = new double[m_iDim];
		for (z=0;z<m_iDim;z++)
			m_pData[z] = f;
	}


	explicit CBQBDVectorN(int i) {
		m_iDim = i;
		m_pData = new double[m_iDim];
	}


	void ZeroVector() {
		memset(m_pData,0,sizeof(double)*m_iDim);
	}


	int GetDim() const {
		return m_iDim;
	}


	double &GetAt(int i) {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::GetAt(int): %d...",i);
		#endif
		#ifdef DEBUG_ARRAYS
			if (i >= m_iDim) {
				printf("& CBQBDVectorN::GetAt(int): Boundary Error (%d/%d).\n",i,m_iDim);
				abort();
			}
		#endif
		return m_pData[i];
	}


	double &operator [] (int i) {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::operator [] (int): %d\n",i);
		#endif
		#ifdef DEBUG_ARRAYS
			if ((i < 0) || (i >= m_iDim)) {
				printf("& CBQBDVectorN::operator [] (int): Boundary Error (%d/%d).\n",i,m_iDim);
				abort();
			}
		#endif
		return m_pData[i];
	}


	double GetAt(int i) const {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::GetAt(int): %d...",i);
		#endif
		#ifdef DEBUG_ARRAYS
			if (i >= m_iDim) {
				printf("CBQBDVectorN::GetAt(int): Boundary Error (%d/%d).\n",i,m_iDim);
				abort();
			}
		#endif
		return m_pData[i];
	}


	double operator [] (int i) const {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::operator [] (int): %d\n",i);
		#endif
		#ifdef DEBUG_ARRAYS
			if ((i < 0) || (i >= m_iDim)) {
				printf("CBQBDVectorN::operator [] (int): Boundary Error (%d/%d).\n",i,m_iDim);
				abort();
			}
		#endif
		return m_pData[i];
	}


	CBQBDVectorN operator + (const CBQBDVectorN &v) const {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::operator + (CBQBDVectorN&)\n");
		#endif
		if (m_iDim != v.m_iDim) {
			printf("CBQBDVectorN::operator + (CBQBDVectorN&): Dimension mismatch (%d vs %d).\n",m_iDim,v.m_iDim);
			abort();
		}
		CBQBDVectorN r(m_iDim);
		int z;
		for (z=0;z<m_iDim;z++)
			r.m_pData[z] = m_pData[z] + v.m_pData[z];
		return r;
	}


	CBQBDVectorN operator - (const CBQBDVectorN &v) const {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::operator - (CBQBDVectorN&)\n");
		#endif
		if (m_iDim != v.m_iDim) {
			printf("CBQBDVectorN::operator - (CBQBDVectorN&): Dimension mismatch (%d vs %d).\n",m_iDim,v.m_iDim);
			abort();
		}
		CBQBDVectorN r(m_iDim);
		int z;
		for (z=0;z<m_iDim;z++)
			r.m_pData[z] = m_pData[z] - v.m_pData[z];
		return r;
	}


	CBQBDVectorN operator * (double f) const {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::operator * (double)\n");
		#endif
		CBQBDVectorN r(m_iDim);
		int z;
		for (z=0;z<m_iDim;z++)
			r.m_pData[z] = m_pData[z] * f;
		return r;
	}


	CBQBDVectorN operator / (double f) const {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::operator / (double)\n");
		#endif
		CBQBDVectorN r(m_iDim);
		int z;
		for (z=0;z<m_iDim;z++)
			r.m_pData[z] = m_pData[z] / f;
		return r;
	}


	void operator += (const CBQBDVectorN &v) {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::operator += (CBQBDVectorN&)\n");
		#endif
		if (m_iDim != v.m_iDim) {
			printf("CBQBDVectorN::operator += (CBQBDVectorN&): Dimension mismatch (%d vs %d).\n",m_iDim,v.m_iDim);
			abort();
		}
		int z;
		for (z=0;z<m_iDim;z++)
			m_pData[z] += v.m_pData[z];
	}


	void operator -= (const CBQBDVectorN &v) {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::operator -= (CBQBDVectorN&)\n");
		#endif
		if (m_iDim != v.m_iDim) {
			printf("CBQBDVectorN::operator -= (CBQBDVectorN&): Dimension mismatch (%d vs %d).\n",m_iDim,v.m_iDim);
			abort();
		}
		int z;
		for (z=0;z<m_iDim;z++)
			m_pData[z] -= v.m_pData[z];
	}


	void operator *= (double f) {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::operator *= (double)\n");
		#endif
		int z;
		for (z=0;z<m_iDim;z++)
			m_pData[z] *= f;
	}


	void operator /= (double f) {
		#ifdef DEBUG_CDVECTORN
			printf("@ CBQBDVectorN::operator /= (double)\n");
		#endif
		int z;
		for (z=0;z<m_iDim;z++)
			m_pData[z] /= f;
	}


	double GetLength() const {
		int z;
		double r=0;
		for (z=0;z<m_iDim;z++)
			r += m_pData[z]*m_pData[z];
		return sqrt(r);
	}


	double GetLengthSqr() const {
		int z;
		double r=0;
		for (z=0;z<m_iDim;z++)
			r += m_pData[z]*m_pData[z];
		return r;
	}


	void Normalize() {
		int z;
		double l;
		l = GetLength();
		for (z=0;z<m_iDim;z++)
			m_pData[z] /= l;
	}


private:
	double *m_pData;
	int m_iDim;
};


inline CBQBDVectorN operator * (double f, const CBQBDVectorN &v) {
	return v*f;
}


inline double VecDist(const CBQBDVectorN &vec1, const CBQBDVectorN &vec2) {
	return (vec1-vec2).GetLength();
}



class CBQBDMatrixMN {
public:

	CBQBDMatrixMN() {
		m_iRows = 0;
		m_iCols = 0;
		m_pData = NULL;
	}


	~CBQBDMatrixMN() {
		if (m_pData != NULL) {
			delete[] m_pData;
			m_pData = NULL;
		}
	}


	CBQBDMatrixMN(const CBQBDMatrixMN &m) {
		m_iCols = m.m_iCols;
		m_iRows = m.m_iRows;
		m_pData = new double[m_iRows*m_iCols];
		memcpy(m_pData,m.m_pData,sizeof(double)*m_iRows*m_iCols);
	}


	CBQBDMatrixMN & operator = (const CBQBDMatrixMN &m) {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::operator = (const CBQBDMatrixMN &)\n");
		#endif
		if (m_pData != NULL)
			delete[] m_pData;
		m_iRows = m.m_iRows;
		m_iCols = m.m_iCols;
		m_pData = new double[m_iRows*m_iCols];
		memcpy(m_pData,m.m_pData,sizeof(double)*m_iRows*m_iCols);
		return *this;
	}


	CBQBDMatrixMN(int m, int n, double f) {
		int z;
		m_iRows = m;
		m_iCols = n;
		m_pData = new double[m_iRows*m_iCols];
		for (z=0;z<m_iRows*m_iCols;z++)
			m_pData[z] = f;
	}


	CBQBDMatrixMN(int m, int n) {
		m_iRows = m;
		m_iCols = n;
		m_pData = new double[m_iRows*m_iCols];
	}


	void ZeroMatrix() {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::ZeroMatrix().\n");
		#endif
		memset(m_pData,0,sizeof(double)*m_iCols*m_iRows);
	}


	int GetRows() const {
		return m_iRows;
	}


	int GetCols() const {
		return m_iCols;
	}


	double &GetAt(int i) {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::GetAt(int): %d.\n",i);
		#endif
		#ifdef BOUND_MATVEC
			if ((i < 0) || (i >= m_iRows*m_iCols)) {
				printf("& CBQBDMatrixMN::GetAt(int): Boundary Error (%d/%d).\n",i,m_iRows*m_iCols);
				abort();
			}
		#endif
		return m_pData[i];
	}


	double &GetAt(int i, int j) {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::GetAt(int, int): %d, %d.\n",i,j);
		#endif
		#ifdef BOUND_MATVEC
			if ((i < 0) || (i >= m_iRows) || (j < 0) || (j >= m_iCols)) {
				printf("& CBQBDMatrixMN::GetAt(int,int): Boundary Error (%d|%d / %d|%d).\n",i,j,m_iRows,m_iCols);
				abort();
			}
		#endif
		return m_pData[i*m_iCols+j];
	}


	double &operator [] (int i) {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::operator [] (int): %d.\n",i);
		#endif
		return GetAt(i);
	}


	double &operator () (int i, int j) {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::operator () (int,int): %d, %d.\n",i,j);
		#endif
		#ifdef BOUND_MATVEC
			if ((i < 0) || (i >= m_iRows) || (j < 0) || (j >= m_iCols)) {
				printf("& CBQBDMatrixMN::operator () (int,int): Boundary Error (%d|%d / %d|%d).\n",i,j,m_iRows,m_iCols);
				abort();
			}
		#endif
		return m_pData[i*m_iCols+j];
	}


	double GetAt(int i) const {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::GetAt(int): %d.\n",i);
		#endif
		#ifdef BOUND_MATVEC
			if ((i < 0) || (i >= m_iRows*m_iCols)) {
				printf("CBQBDMatrixMN::GetAt(int): Boundary Error (%d/%d).\n",i,m_iRows*m_iCols);
				abort();
			}
		#endif
		return m_pData[i];
	}


	double GetAt(int i, int j) const {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::GetAt(int,int): %d, %d.\n",i,j);
		#endif
		#ifdef BOUND_MATVEC
			if ((i < 0) || (i >= m_iRows) || (j < 0) || (j >= m_iCols)) {
				printf("CBQBDMatrixMN::GetAt(int,int): Boundary Error (%d|%d / %d|%d).\n",i,j,m_iRows,m_iCols);
				abort();
			}
		#endif
		return m_pData[i*m_iCols+j];
	}


	double operator [] (int i) const {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::operator [] (int): %d.\n",i);
		#endif
		#ifdef BOUND_MATVEC
			if ((i < 0) || (i >= m_iRows*m_iCols)) {
				printf("CBQBDMatrixMN::operator [] (int): Boundary Error (%d/%d).\n",i,m_iRows*m_iCols);
				abort();
			}
		#endif
		return m_pData[i];
	}


	double operator () (int i, int j) const {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::operator () (int,int): %d, %d.\n",i,j);
		#endif
		#ifdef BOUND_MATVEC
			if ((i < 0) || (i >= m_iRows) || (j < 0) || (j >= m_iCols)) {
				printf("CBQBDMatrixMN::operator () (int,int): Boundary Error (%d|%d / %d|%d).\n",i,j,m_iRows,m_iCols);
				abort();
			}
		#endif
		return m_pData[i*m_iCols+j];
	}


	void operator *= (double f) {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::operator *= (double)\n");
		#endif
		int z;
		for (z=0;z<m_iRows*m_iCols;z++)
			m_pData[z] *= f;
	}


	void operator += (const CBQBDMatrixMN &m) {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::operator += (const CBQBDMatrixMN)\n");
		#endif
		if ((m_iCols != m.m_iCols) || (m_iRows != m.m_iRows)) {
			printf("@ CBQBDMatrixMN::operator += (CBQBDMatrixMN): Dimension mismatch (%d x %d vs %d x %d).\n",m_iCols,m_iRows,m.m_iCols,m.m_iRows);
			abort();
		}
		int z;
		for (z=0;z<m_iRows*m_iCols;z++)
			m_pData[z] += m.m_pData[z];
	}


	CBQBDMatrixMN operator * (const CBQBDMatrixMN &m) const {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::operator * (CBQBDMatrixMN)\n");
		#endif
		if (m_iCols != m.m_iRows) {
			printf("@ CBQBDMatrixMN::operator * (CBQBDMatrixMN): Dimension mismatch (%d vs %d).\n",m_iCols,m.m_iRows);
			abort();
		}
		int z1, z2, z3;
		CBQBDMatrixMN r = CBQBDMatrixMN(m_iRows,m.m_iCols);
		r.ZeroMatrix();
		for (z1=0;z1<m_iRows;z1++)
			for (z2=0;z2<m.m_iCols;z2++)
				for (z3=0;z3<m_iCols;z3++)
					r.GetAt(z1,z2) += GetAt(z1,z3) * m.GetAt(z3,z2);
		return r;
	}


	CBQBDVectorN operator * (const CBQBDVectorN &v) const {
		#ifdef DEBUG_CDMATRIXMN
			printf("@ CBQBDMatrixMN::operator * (CBQBDVector3)\n");
		#endif
		if (m_iCols != v.GetDim()) {
			printf("@ CBQBDMatrixMN::operator * (CBQBDVectorN): Dimension mismatch (%d vs %d).\n",m_iCols,v.GetDim());
			abort();
		}
		int z1, z3;
		CBQBDVectorN r = CBQBDVectorN(v.GetDim());
		r.ZeroVector();
		for (z1=0;z1<m_iRows;z1++)
			for (z3=0;z3<m_iCols;z3++)
				r.GetAt(z1) += GetAt(z1,z3) * v.GetAt(z3);
		return r;
	}


	void Unity() { 
		if (m_iCols != m_iRows) {
			printf("CBQBDMatrixMN::Unity(): Only defined for square matrices (have %dx%d).\n",m_iRows,m_iCols);
			abort();
		}
		int z;
		ZeroMatrix();
		for (z=0;z<m_iCols;z++)
			GetAt(z,z) = 1.0;
	}


	CBQBDMatrixMN Transpose() const {
		int z1, z2;
		CBQBDMatrixMN r=CBQBDMatrixMN(m_iCols,m_iRows);
		for (z1=0;z1<m_iRows;z1++)
			for (z2=0;z2<m_iCols;z2++)
				r.GetAt(z2,z1) = GetAt(z1,z2);
		return r;
	}


	void Dump(FILE *a) const;

	void DumpSmall(FILE *a) const;

private:
	double *m_pData;
	int m_iCols;
	int m_iRows;
};



#endif // Not BQB_INSIDE_TRAVIS


#endif



