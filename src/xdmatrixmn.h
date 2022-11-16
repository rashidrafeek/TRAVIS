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


#ifndef XDMATRIXMN_H
#define XDMATRIXMN_H


// This must always be the first include directive
#include "config.h"

#include "tools.h"
#include "xdvectorn.h"


#define BOUND_MATVEC


class CxDMatrixMN : public CxObject
{
public:

	CxDMatrixMN()
	{
		m_iRows = 0;
		m_iCols = 0;
		m_pData = NULL;
	}


	~CxDMatrixMN()
	{
		if (m_pData != NULL)
		{
			delete[] m_pData;
			m_pData = NULL;
		}
	}


	CxDMatrixMN(const CxDMatrixMN &m)
	{
		BXIN;
		m_iCols = m.m_iCols;
		m_iRows = m.m_iRows;
		m_pData = new double[m_iRows*m_iCols];
		memcpy(m_pData,m.m_pData,sizeof(double)*m_iRows*m_iCols);
		BXOUT;
	}


	CxDMatrixMN & operator = (const CxDMatrixMN &m)
	{
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::operator = (const CxDMatrixMN &)\n");
		#endif
		BXIN;
		if (m_pData != NULL)
			delete[] m_pData;
		m_iRows = m.m_iRows;
		m_iCols = m.m_iCols;
		m_pData = new double[m_iRows*m_iCols];
		memcpy(m_pData,m.m_pData,sizeof(double)*m_iRows*m_iCols);
		BXOUT;
		return *this;
	}


	CxDMatrixMN(int m, int n, double f)
	{
		BXIN;
		int z;
		m_iRows = m;
		m_iCols = n;
		m_pData = new double[m_iRows*m_iCols];
		for (z=0;z<m_iRows*m_iCols;z++)
			m_pData[z] = f;
		BXOUT;
	}


	CxDMatrixMN(int m, int n)
	{
		BXIN;
		m_iRows = m;
		m_iCols = n;
		m_pData = new double[m_iRows*m_iCols];
		BXOUT;
	}


	void ZeroMatrix()
	{
		BXIN;
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::ZeroMatrix().\n");
		#endif
		memset(m_pData,0,sizeof(double)*m_iCols*m_iRows);
		BXOUT;
	}


	int GetRows() const
	{
		return m_iRows;
	}


	int GetCols() const
	{
		return m_iCols;
	}


	double &GetAt(int i)
	{
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::GetAt(int): %d.\n",i);
		#endif
		#ifdef BOUND_MATVEC
		if ((i < 0) || (i >= m_iRows*m_iCols))
		{
			eprintf("& CxDMatrixMN::GetAt(int): Boundary Error (%d/%d).\n",i,m_iRows*m_iCols);
			abort();
		}
		#endif
		return m_pData[i];
	}


	double &GetAt(int i, int j)
	{
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::GetAt(int, int): %d, %d.\n",i,j);
		#endif
		#ifdef BOUND_MATVEC
		if ((i < 0) || (i >= m_iRows) || (j < 0) || (j >= m_iCols))
		{
			eprintf("& CxDMatrixMN::GetAt(int,int): Boundary Error (%d|%d / %d|%d).\n",i,j,m_iRows,m_iCols);
			abort();
		}
		#endif
		return m_pData[i*m_iCols+j];
	}


	double &operator [] (int i)
	{
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::operator [] (int): %d.\n",i);
		#endif
		return GetAt(i);
	}


	double &operator () (int i, int j)
	{
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::operator () (int,int): %d, %d.\n",i,j);
		#endif
		#ifdef BOUND_MATVEC
		if ((i < 0) || (i >= m_iRows) || (j < 0) || (j >= m_iCols))
		{
			eprintf("& CxDMatrixMN::operator () (int,int): Boundary Error (%d|%d / %d|%d).\n",i,j,m_iRows,m_iCols);
			abort();
		}
		#endif
		return m_pData[i*m_iCols+j];
	}


	double GetAt(int i) const
	{
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::GetAt(int): %d.\n",i);
		#endif
		#ifdef BOUND_MATVEC
		if ((i < 0) || (i >= m_iRows*m_iCols))
		{
			eprintf("CxDMatrixMN::GetAt(int): Boundary Error (%d/%d).\n",i,m_iRows*m_iCols);
			abort();
		}
		#endif
		return m_pData[i];
	}


	double GetAt(int i, int j) const
	{
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::GetAt(int,int): %d, %d.\n",i,j);
		#endif
		#ifdef BOUND_MATVEC
		if ((i < 0) || (i >= m_iRows) || (j < 0) || (j >= m_iCols))
		{
			eprintf("CxDMatrixMN::GetAt(int,int): Boundary Error (%d|%d / %d|%d).\n",i,j,m_iRows,m_iCols);
			abort();
		}
		#endif
		return m_pData[i*m_iCols+j];
	}


	double operator [] (int i) const
	{
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::operator [] (int): %d.\n",i);
		#endif
		#ifdef BOUND_MATVEC
		if ((i < 0) || (i >= m_iRows*m_iCols))
		{
			eprintf("CxDMatrixMN::operator [] (int): Boundary Error (%d/%d).\n",i,m_iRows*m_iCols);
			abort();
		}
		#endif
		return m_pData[i];
	}


	double operator () (int i, int j) const
	{
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::operator () (int,int): %d, %d.\n",i,j);
		#endif
		#ifdef BOUND_MATVEC
		if ((i < 0) || (i >= m_iRows) || (j < 0) || (j >= m_iCols))
		{
			eprintf("CxDMatrixMN::operator () (int,int): Boundary Error (%d|%d / %d|%d).\n",i,j,m_iRows,m_iCols);
			abort();
		}
		#endif
		return m_pData[i*m_iCols+j];
	}


	void operator *= (double f)
	{
		BXIN;
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::operator *= (double)\n");
		#endif
		int z;
		for (z=0;z<m_iRows*m_iCols;z++)
			m_pData[z] *= f;
		BXOUT;
	}


	void operator += (const CxDMatrixMN &m)
	{
		BXIN;
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::operator += (const CxDMatrixMN)\n");
		#endif
		if ((m_iCols != m.m_iCols) || (m_iRows != m.m_iRows))
		{
			eprintf("@ CxDMatrixMN::operator += (CxDMatrixMN): Dimension mismatch (%d x %d vs %d x %d).\n",m_iCols,m_iRows,m.m_iCols,m.m_iRows);
			abort();
		}
		int z;
		for (z=0;z<m_iRows*m_iCols;z++)
			m_pData[z] += m.m_pData[z];
		BXOUT;
	}


	CxDMatrixMN operator * (const CxDMatrixMN &m) const
	{
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::operator * (CxDMatrixMN)\n");
		#endif
		if (m_iCols != m.m_iRows)
		{
			eprintf("@ CxDMatrixMN::operator * (CxDMatrixMN): Dimension mismatch (%d vs %d).\n",m_iCols,m.m_iRows);
			abort();
		}
		int z1, z2, z3;
		CxDMatrixMN r = CxDMatrixMN(m_iRows,m.m_iCols);
		r.ZeroMatrix();
		for (z1=0;z1<m_iRows;z1++)
			for (z2=0;z2<m.m_iCols;z2++)
				for (z3=0;z3<m_iCols;z3++)
					r.GetAt(z1,z2) += GetAt(z1,z3) * m.GetAt(z3,z2);
		return r;
	}


	CxDVectorN operator * (const CxDVectorN &v) const
	{
		#ifdef DEBUG_CDMATRIXMN
		mprintf("@ CxDMatrixMN::operator * (CxDVector3)\n");
		#endif
		if (m_iCols != v.GetDim())
		{
			eprintf("@ CxDMatrixMN::operator * (CxDVectorN): Dimension mismatch (%d vs %d).\n",m_iCols,v.GetDim());
			abort();
		}
		int z1, z3;
		CxDVectorN r = CxDVectorN(v.GetDim());
		r.ZeroVector();
		for (z1=0;z1<m_iRows;z1++)
		{
			for (z3=0;z3<m_iCols;z3++)
				r.GetAt(z1) += GetAt(z1,z3) * v.GetAt(z3);
		}
		return r;
	}


	void Unity()
	{ 
		BXIN;
		if (m_iCols != m_iRows)
		{
			eprintf("CxDMatrixMN::Unity(): Only defined for square matrices (have %dx%d).\n",m_iRows,m_iCols);
			abort();
		}
		int z;
		ZeroMatrix();
		for (z=0;z<m_iCols;z++)
			GetAt(z,z) = 1.0;
		BXOUT;
	}


	CxDMatrixMN Transpose() const
	{
		int z1, z2;
		CxDMatrixMN r=CxDMatrixMN(m_iCols,m_iRows);
		for (z1=0;z1<m_iRows;z1++)
			for (z2=0;z2<m_iCols;z2++)
				r.GetAt(z2,z1) = GetAt(z1,z2);
		return r;
	}


	void Dump() const;

	void Dump(FILE *a) const;

	void DumpSmall() const;

	void DumpSmall(FILE *a) const;

private:
	double *m_pData;
	int m_iCols;
	int m_iRows;
};


#endif
