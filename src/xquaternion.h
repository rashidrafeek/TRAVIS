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


#ifndef CXQUATERNION_DEFINED
#define CXQUATERNION_DEFINED


// This must always be the first include directive
#include "config.h"

#include "tools.h"
#include "xobject.h"
#include "backtrace.h"


class CxQuaternion : public CxObject
{
public:

	CxQuaternion() { }


	~CxQuaternion() { }

	
	CxQuaternion(const CxQuaternion &v) : CxObject()
	{
		BXIN;
		m_pData[0] = v.m_pData[0];
		m_pData[1] = v.m_pData[1];
		m_pData[2] = v.m_pData[2];
		m_pData[3] = v.m_pData[3];
		BXOUT;
	}


	explicit CxQuaternion(double f)
	{
		BXIN;
		m_pData[0] = f;
		m_pData[1] = f;
		m_pData[2] = f;
		m_pData[3] = f;
		BXOUT;
	}


	CxQuaternion(double a, double b, double c, double d)
	{
		BXIN;
		m_pData[0] = a;
		m_pData[1] = b;
		m_pData[2] = c;
		m_pData[3] = d;
		BXOUT;
	}


	CxQuaternion & operator = (const CxQuaternion &v) {
		m_pData[0] = v.m_pData[0];
		m_pData[1] = v.m_pData[1];
		m_pData[2] = v.m_pData[2];
		m_pData[3] = v.m_pData[3];
		return *this;
	}


	double &GetAt(int i)
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}


	double &operator [] (int i)
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	double GetAt(int i) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::GetAt(int): %d...",i);
		#endif
		#ifdef DEBUG_ARRAYS
		if (i > 3)
		{
			mprintf("CxQuaternion::GetAt(int): Boundary Error (%d/3)...",i);
			abort();
		}
		#endif
		return m_pData[i];
	}


	double operator [] (int i) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	CxQuaternion operator + (const CxQuaternion &v) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator + (CxQuaternion&)\n");
		#endif
		return CxQuaternion(m_pData[0]+v.m_pData[0],m_pData[1]+v.m_pData[1],m_pData[2]+v.m_pData[2],m_pData[3]+v.m_pData[3]);
	}


	CxQuaternion operator - (const CxQuaternion &v) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator - (CxQuaternion&)\n");
		#endif
		return CxQuaternion(m_pData[0]-v.m_pData[0],m_pData[1]-v.m_pData[1],m_pData[2]-v.m_pData[2],m_pData[3]-v.m_pData[3]);
	}


	CxQuaternion operator * (const CxQuaternion &v) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator + (CxQuaternion&)\n");
		#endif
		return CxQuaternion(
			m_pData[0]*v.m_pData[0] - m_pData[1]*v.m_pData[1] - m_pData[2]*v.m_pData[2] - m_pData[3]*v.m_pData[3],
			m_pData[0]*v.m_pData[1] + m_pData[1]*v.m_pData[0] + m_pData[2]*v.m_pData[3] - m_pData[3]*v.m_pData[2],
			m_pData[0]*v.m_pData[2] - m_pData[1]*v.m_pData[3] + m_pData[2]*v.m_pData[0] + m_pData[3]*v.m_pData[1],
			m_pData[0]*v.m_pData[3] + m_pData[1]*v.m_pData[2] - m_pData[2]*v.m_pData[1] + m_pData[3]*v.m_pData[0] );
	}


	CxQuaternion operator * (double f) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator * (double)\n");
		#endif
		return CxQuaternion(m_pData[0]*f,m_pData[1]*f,m_pData[2]*f,m_pData[3]*f);
	}


	CxQuaternion operator / (double f) const
	{
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator / (double)\n");
		#endif
		return CxQuaternion(m_pData[0]/f,m_pData[1]/f,m_pData[2]/f,m_pData[3]/f);
	}


	void operator += (const CxQuaternion &v)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator += (CxQuaternion&)\n");
		#endif
		m_pData[0] += v.m_pData[0];
		m_pData[1] += v.m_pData[1];
		m_pData[2] += v.m_pData[2];
		m_pData[3] += v.m_pData[3];
		BXOUT;
	}


	void operator -= (const CxQuaternion &v)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator -= (CxQuaternion&)\n");
		#endif
		m_pData[0] -= v.m_pData[0];
		m_pData[1] -= v.m_pData[1];
		m_pData[2] -= v.m_pData[2];
		m_pData[3] -= v.m_pData[3];
		BXOUT;
	}


	void operator *= (double f)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator *= (double)\n");
		#endif
		m_pData[0] *= f;
		m_pData[1] *= f;
		m_pData[2] *= f;
		m_pData[3] *= f;
		BXOUT;
	}


	void operator /= (double f)
	{
		BXIN;
		#ifdef DEBUG_CVECTOR3
		mprintf("@ CxQuaternion::operator /= (double)\n");
		#endif
		m_pData[0] /= f;
		m_pData[1] /= f;
		m_pData[2] /= f;
		m_pData[3] /= f;
		BXOUT;
	}


	double GetLength() const
	{
		return sqrt(m_pData[0]*m_pData[0]+m_pData[1]*m_pData[1]+m_pData[2]*m_pData[2]+m_pData[3]*m_pData[3]);
	}


	double GetLengthSqr() const
	{
		return m_pData[0]*m_pData[0]+m_pData[1]*m_pData[1]+m_pData[2]*m_pData[2]+m_pData[3]*m_pData[3];
	}


	CxQuaternion Conjugate() const
	{
		return CxQuaternion(m_pData[0],-m_pData[1],-m_pData[2],-m_pData[3]);
	}


	void BuildRotation(CxDVector3 axis, double angle)
	{
		CxDVector3 tmp;
		tmp = axis;
		tmp.Normalize();
		m_pData[0] = cos(angle/2.0);
		m_pData[1] = sin(angle/2.0) * axis[0];
		m_pData[2] = sin(angle/2.0) * axis[1];
		m_pData[3] = sin(angle/2.0) * axis[2];
	}


/*	inline CxVector3 Transform(const CxDVector3 &v) const
	{
		CxQuaternion out;
		out = this->Conjugate() * CxQuaternion(0,v[0],v[1],v[2]) * (*this);
		return CxDVector3(out[1],out[2],out[3]);
	}*/


	CxDVector3 Transform(const CxDVector3 &v) const
	{
		CxQuaternion out;
		out = this->Conjugate() * CxQuaternion(0,v[0],v[1],v[2]) * (*this);
		return CxDVector3(out[1],out[2],out[3]);
	}


	void Normalize()
	{
		BXIN;
		double l;
		l = GetLength();
		m_pData[0] /= l;
		m_pData[1] /= l;
		m_pData[2] /= l;
		m_pData[3] /= l;
		BXOUT;
	}


	void Unity()
	{
		BXIN;
		m_pData[0] = 1.0;
		m_pData[1] = 0;
		m_pData[2] = 0;
		m_pData[3] = 0;
		BXOUT;
	}


private:
	double m_pData[4];
};


#endif
