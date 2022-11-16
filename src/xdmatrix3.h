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


#ifndef CXDMATRIX3_DEFINED
#define CXDMATRIX3_DEFINED


// This must always be the first include directive
#include "config.h"

#include "tools.h"
#include "xobject.h"
#include "xdvector3.h"
#include "backtrace.h"


class CxDMatrix3 : public CxObject {

public:

	CxDMatrix3() { }

	~CxDMatrix3() { }


/*	CxDMatrix3(CxDMatrix3 m)
	{
		m_pData[0] = m.m_pData[0];
		m_pData[1] = m.m_pData[1];
		m_pData[2] = m.m_pData[2];
		m_pData[3] = m.m_pData[3];
		m_pData[4] = m.m_pData[4];
		m_pData[5] = m.m_pData[5];
		m_pData[6] = m.m_pData[6];
		m_pData[7] = m.m_pData[7];
		m_pData[8] = m.m_pData[8];
	}*/
	

	CxDMatrix3(const CxDMatrix3 &m) : CxObject()
	{
		BXIN;
		m_pData[0] = m.m_pData[0];
		m_pData[1] = m.m_pData[1];
		m_pData[2] = m.m_pData[2];
		m_pData[3] = m.m_pData[3];
		m_pData[4] = m.m_pData[4];
		m_pData[5] = m.m_pData[5];
		m_pData[6] = m.m_pData[6];
		m_pData[7] = m.m_pData[7];
		m_pData[8] = m.m_pData[8];
		BXOUT;
	}


	explicit CxDMatrix3(double f)
	{
		BXIN;
		m_pData[0] = f;
		m_pData[1] = f;
		m_pData[2] = f;
		m_pData[3] = f;
		m_pData[4] = f;
		m_pData[5] = f;
		m_pData[6] = f;
		m_pData[7] = f;
		m_pData[8] = f;
		BXOUT;
	}


	CxDMatrix3(double a, double b, double c, double d, double e, double f, double g, double h, double i)
	{
		BXIN;
		m_pData[0] = a;
		m_pData[1] = b;
		m_pData[2] = c;
		m_pData[3] = d;
		m_pData[4] = e;
		m_pData[5] = f;
		m_pData[6] = g;
		m_pData[7] = h;
		m_pData[8] = i;
		BXOUT;
	}


	CxDMatrix3 & operator = (const CxDMatrix3 &v) {
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


	double &GetAt(int i)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}


	double &GetAt(int i, int j)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::GetAt(int, int): %d, %d...",i,j);
		#endif
		return m_pData[i*3+j];
	}


	double &operator [] (int i)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	double GetAt(int i) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}


	double GetAt(int i, int j) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::GetAt(int, int): %d, %d...",i,j);
		#endif
		return m_pData[i*3+j];
	}


	double operator [] (int i) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	double &operator () (int i, int j)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator () (int, int): %d, %d\n",i,j);
		#endif
		return GetAt(i,j);
	}


	double operator () (int i, int j) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator () (int, int): %d, %d\n",i,j);
		#endif
		return GetAt(i,j);
	}


	void operator *= (double f)
	{
		BXIN;
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator *= (double)\n");
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
		BXOUT;
	}


	CxDMatrix3 operator * (const CxDMatrix3 &m) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator * (CxDMatrix3)\n");
		#endif
		return CxDMatrix3(
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


	void operator += (const CxDMatrix3 &m)
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator += (CxDMatrix3)\n");
		#endif
		m_pData[0] += m.m_pData[0];
		m_pData[1] += m.m_pData[1];
		m_pData[2] += m.m_pData[2];
		m_pData[3] += m.m_pData[3];
		m_pData[4] += m.m_pData[4];
		m_pData[5] += m.m_pData[5];
		m_pData[6] += m.m_pData[6];
		m_pData[7] += m.m_pData[7];
		m_pData[8] += m.m_pData[8];
	}


	CxDVector3 operator * (const CxDVector3 &v) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator * (CxDVector3)\n");
		#endif
		return CxDVector3(
			GetAt(0,0)*v[0] + GetAt(1,0)*v[1] + GetAt(2,0)*v[2],
			GetAt(0,1)*v[0] + GetAt(1,1)*v[1] + GetAt(2,1)*v[2],
			GetAt(0,2)*v[0] + GetAt(1,2)*v[1] + GetAt(2,2)*v[2]);
	}


	CxDMatrix3 operator * (double f) const
	{
		#ifdef DEBUG_CMATRIX3
		mprintf("@ CxDMatrix3::operator * (double)\n");
		#endif
		return CxDMatrix3(m_pData[0]*f,m_pData[1]*f,m_pData[2]*f,m_pData[3]*f,m_pData[4]*f,m_pData[5]*f,m_pData[6]*f,m_pData[7]*f,m_pData[8]*f);
	}


	void RotMat(const CxDVector3 &vec, double a)
	{
		BXIN;
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
		BXOUT;
	}


	void Unity()
	{ 
		BXIN;
		GetAt(0,0) = 1;
		GetAt(1,0) = 0;
		GetAt(2,0) = 0;
		GetAt(0,1) = 0;
		GetAt(1,1) = 1;
		GetAt(2,1) = 0;
		GetAt(0,2) = 0;
		GetAt(1,2) = 0;
		GetAt(2,2) = 1;
		BXOUT;
	}


	CxDMatrix3 Transpose()
	{
		return CxDMatrix3(
			GetAt(0,0), GetAt(1,0), GetAt(2,0),
			GetAt(0,1), GetAt(1,1), GetAt(2,1),
			GetAt(0,2), GetAt(1,2), GetAt(2,2));
	}


	double Det() {
		#define a11 GetAt(0,0)
		#define a12 GetAt(0,1)
		#define a13 GetAt(0,2)
		#define a21 GetAt(1,0)
		#define a22 GetAt(1,1)
		#define a23 GetAt(1,2)
		#define a31 GetAt(2,0)
		#define a32 GetAt(2,1)
		#define a33 GetAt(2,2)

		return  a11 * (a33*a22 - a32*a23)
			  - a21 * (a33*a12 - a32*a13)
			  + a31 * (a23*a12 - a22*a13);
	}


	void Dump() const;

	void DumpSmall() const;

	bool Invert();

	void MatUltra(const CxDVector3 &vec1, const CxDVector3 &vec2);

	void FindRotZ(const CxDVector3 &v);

private:
	double m_pData[9];
};

double Abs_Dihedral(const CxDVector3 &v1, const CxDVector3 &v2, const CxDVector3 &v3);

#endif


