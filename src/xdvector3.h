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


#ifndef CXDVECTOR3_DEFINED
#define CXDVECTOR3_DEFINED


// This must always be the first include directive
#include "config.h"

#include "tools.h"
#include "xobject.h"
#include "backtrace.h"



class CxDVector3 : public CxObject
{
public:
	CxDVector3()
	{
	}


	~CxDVector3()
	{
	}
	

	CxDVector3(const CxDVector3 &v) : CxObject()
	{
		BXIN;
		m_pData[0] = v.m_pData[0];
		m_pData[1] = v.m_pData[1];
		m_pData[2] = v.m_pData[2];
		BXOUT;
	}


//	CxDVector3(const CxVector3 &v);


	explicit CxDVector3(double f)
	{
		BXIN;
		m_pData[0] = f;
		m_pData[1] = f;
		m_pData[2] = f;
		BXOUT;
	}


	CxDVector3(double x, double y, double z)
	{
		BXIN;
		m_pData[0] = x;
		m_pData[1] = y;
		m_pData[2] = z;
		BXOUT;
	}


	CxDVector3 & operator = (const CxDVector3 &v) {
		m_pData[0] = v.m_pData[0];
		m_pData[1] = v.m_pData[1];
		m_pData[2] = v.m_pData[2];
		return *this;
	}


	double &GetAt(int i)
	{
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::GetAt(int): %d...",i);
		#endif
		return m_pData[i];
	}


	double &operator [] (int i)
	{
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	double GetAt(int i) const
	{
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::GetAt(int): %d...",i);
		#endif
		#ifdef DEBUG_ARRAYS
		if (i > 2)
		{
			mprintf("CxDVector3::GetAt(int): Boundary Error (%d/3)...",i);
			abort();
		}
		#endif
		return m_pData[i];
	}


	double operator [] (int i) const
	{
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::operator [] (int): %d\n",i);
		#endif
		return GetAt(i);
	}


	CxDVector3 operator + (const CxDVector3 &v) const
	{
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::operator + (CxDVector3&)\n");
		#endif
		return CxDVector3(m_pData[0]+v.m_pData[0],m_pData[1]+v.m_pData[1],m_pData[2]+v.m_pData[2]);
	}


	CxDVector3 operator - (const CxDVector3 &v) const
	{
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::operator - (CxDVector3&)\n");
		#endif
		return CxDVector3(m_pData[0]-v.m_pData[0],m_pData[1]-v.m_pData[1],m_pData[2]-v.m_pData[2]);
	}


	CxDVector3 operator * (double f) const
	{
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::operator * (double)\n");
		#endif
		return CxDVector3(m_pData[0]*f,m_pData[1]*f,m_pData[2]*f);
	}


	CxDVector3 operator / (double f) const
	{
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::operator / (double)\n");
		#endif
		return CxDVector3(m_pData[0]/f,m_pData[1]/f,m_pData[2]/f);
	}


	void operator += (const CxDVector3 &v)
	{
		BXIN;
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::operator += (CxDVector3&)\n");
		#endif
		m_pData[0] += v.m_pData[0];
		m_pData[1] += v.m_pData[1];
		m_pData[2] += v.m_pData[2];
		BXOUT;
	}


	void operator -= (const CxDVector3 &v)
	{
		BXIN;
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::operator -= (CxDVector3&)\n");
		#endif
		m_pData[0] -= v.m_pData[0];
		m_pData[1] -= v.m_pData[1];
		m_pData[2] -= v.m_pData[2];
		BXOUT;
	}


	void operator *= (double f)
	{
		BXIN;
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::operator *= (double)\n");
		#endif
		m_pData[0] *= f;
		m_pData[1] *= f;
		m_pData[2] *= f;
		BXOUT;
	}


	void operator /= (double f)
	{
		BXIN;
		#ifdef DEBUG_CDVECTOR3
		mprintf("@ CxDVector3::operator /= (double)\n");
		#endif
		m_pData[0] /= f;
		m_pData[1] /= f;
		m_pData[2] /= f;
		BXOUT;
	}


	double GetLength() const
	{
		return sqrt(m_pData[0]*m_pData[0]+m_pData[1]*m_pData[1]+m_pData[2]*m_pData[2]);
	}


	double GetLengthSqr() const
	{
		return m_pData[0]*m_pData[0]+m_pData[1]*m_pData[1]+m_pData[2]*m_pData[2];
	}


	void Normalize()
	{
		BXIN;
		double l;
		l = GetLength();
		m_pData[0] /= l;
		m_pData[1] /= l;
		m_pData[2] /= l;
		BXOUT;
	}


	void Chop(double d)
	{
		if ((m_pData[0] != 0) && (fabs(m_pData[0]) < d))
			m_pData[0] = 0;
		if ((m_pData[1] != 0) && (fabs(m_pData[1]) < d))
			m_pData[1] = 0;
		if ((m_pData[2] != 0) && (fabs(m_pData[2]) < d))
			m_pData[2] = 0;
	}


	void PointRoot(const CxDVector3 &vec1, const CxDVector3 &vec2, const CxDVector3 &point);


	void Dump() const;


private:
	double m_pData[3];
};


double Dihedral(const CxDVector3 &vec1, const CxDVector3 &vec2, const CxDVector3 &norm, bool absolute);


inline void Swap(CxDVector3 &vec1, CxDVector3 &vec2)
{
	BXIN;
	CxDVector3 t;
	t = vec1;
	vec1 = vec2;
	vec2 = t;
	BXOUT;
}


inline CxDVector3 operator * (double f, const CxDVector3 &v)
{
	return v*f;
}


inline CxDVector3 CrossP(const CxDVector3 &vec1, const CxDVector3 &vec2)
{
	return CxDVector3(vec1[1]*vec2[2] - vec1[2]*vec2[1],
		vec1[2]*vec2[0] - vec1[0]*vec2[2],
		vec1[0]*vec2[1] - vec1[1]*vec2[0]);
}


inline double DotP(const CxDVector3 &vec1, const CxDVector3 &vec2)
{
	return vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2];
}


inline double Angle(const CxDVector3 &vec1, const CxDVector3 &vec2)
{
	BXIN;
	double t;

	if ((vec1.GetLength() == 0) || (vec2.GetLength() == 0))
	{
		mprintf("\nAngle(): Warning: Indeterminate angle between ( %f | %f | %f ) and ( %f | %f | %f ).  ",vec1[0],vec1[1],vec1[2],vec2[0],vec2[1],vec2[2]);
		return 0;
	}

	t = DotP(vec1,vec2) / vec1.GetLength() / vec2.GetLength();
	if (t > 1.0)
		t = 1.0;
	if (t < -1.0)
		t = -1.0;
	BXOUT;
	return acos(t);
}


inline double Angle_Deg(const CxDVector3 &vec1, const CxDVector3 &vec2)
{
	BXIN;
	double t;

	if ((vec1.GetLength() == 0) || (vec2.GetLength() == 0))
	{
		mprintf("\nAngle_Deg(): Warning: Indeterminate angle between ( %f | %f | %f ) and ( %f | %f | %f ).  ",vec1[0],vec1[1],vec1[2],vec2[0],vec2[1],vec2[2]);
		return 0;
	}

	t = DotP(vec1,vec2) / vec1.GetLength() / vec2.GetLength();
	if (t > 1.0)
		t = 1.0;
	if (t < -1.0)
		t = -1.0;
	t = acos(t);
	BXOUT;
	return fabs(t*180.0 / Pi);
}


inline double VecDist(const CxDVector3 &vec1, const CxDVector3 &vec2)
{
	return sqrt((vec1[0]-vec2[0])*(vec1[0]-vec2[0]) + (vec1[1]-vec2[1])*(vec1[1]-vec2[1]) + (vec1[2]-vec2[2])*(vec1[2]-vec2[2]));
}


inline CxDVector3 Normalize(const CxDVector3 &vec)
{
	double tf;

	tf = vec.GetLength();

	return CxDVector3(vec[0]/tf,vec[1]/tf,vec[2]/tf);
}


CxDVector3 PointFromRAD(CxDVector3 r1, CxDVector3 r2, CxDVector3 r3, double r, double a, double d);

void PointRootCoefficients(const CxDVector3 &vec1, const CxDVector3 &vec2, const CxDVector3 &point, double &a, double &b);



#endif


