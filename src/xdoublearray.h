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


#ifndef XDOUBLEARRAY_H
#define XDOUBLEARRAY_H


// This must always be the first include directive
#include "config.h"

#include "tools.h"
#include "xobject.h"
#include "backtrace.h"


class CxDoubleArray : public CxObject
{
public:
	
	CxDoubleArray();
	explicit CxDoubleArray(const char *name);
	void SetName(const char *name);
	~CxDoubleArray();
	CxDoubleArray(const CxDoubleArray &o);
	void Add(double f);
	void SetAt(unsigned long pos, double f);
	void SetSize(unsigned long i);
	void SetMaxSize(unsigned long i);
	void SetGrow(unsigned long i);
	void RemoveAll();
	void RemoveAll_KeepSize();
	void RemoveAt(unsigned long pos, unsigned long count);
	void InsertAt(double f, unsigned long pos);
	void CopyFrom(CxDoubleArray *d);


	double &GetAt(unsigned long i)
	{
		BXIN;
	#ifdef DEBUG_CDOUBLEARRAY
		mprintf("@ CxDoubleArray::GetAt(int): %d...",i);
	#endif
#ifdef DEBUG_ARRAYS
		if (i >= m_iSize)
		{
			if (m_sName != NULL)
				eprintf("CxDoubleArray \"%s\" Boundary Error (%lu/%lu).\n",m_sName,i,m_iSize);
					else eprintf("CxDoubleArray Boundary Error (%lu/%lu).\n",i,m_iSize);
			abort();
		}
#endif
	#ifdef DEBUG_CDOUBLEARRAY
		mprintf("done.\n");
	#endif
		BXOUT;
		return m_pData[i];
	}
		

	double & operator [] (unsigned long i)
	{
	#ifdef DEBUG_CDOUBLEARRAY
		mprintf("@ CxDoubleArray::operator [] (int): %d\n",i);
	#endif
		return GetAt(i);
	}


	double GetAt(unsigned long i) const
	{
		BXIN;
	#ifdef DEBUG_CDOUBLEARRAY
		mprintf("@ CxDoubleArray::GetAt(int): %d...",i);
	#endif
#ifdef DEBUG_ARRAYS
		if (i >= m_iSize)
		{
			if (m_sName != NULL)
				eprintf("CxDoubleArray \"%s\" Boundary Error (%lu/%lu).\n",m_sName,i,m_iSize);
					else eprintf("CxDoubleArray Boundary Error (%lu/%lu).\n",i,m_iSize);
			abort();
		}
#endif
	#ifdef DEBUG_CDOUBLEARRAY
		mprintf("done.\n");
	#endif
		BXOUT;
		return m_pData[i];
	}
		

	double operator [] (unsigned long i) const
	{
	#ifdef DEBUG_CDOUBLEARRAY
		mprintf("@ CxDoubleArray::operator [] (int): %d\n",i);
	#endif
		return GetAt(i);
	}


	int GetSize()
	{
	#ifdef DEBUG_CDOUBLEARRAY
		mprintf("@ CxDoubleArray::GetSize(): %d\n",m_iSize);
	#endif
		return m_iSize;
	}	
	

private:	
	double *m_pData;
	unsigned long m_iSize;
	unsigned long m_iMaxSize;
	unsigned long m_iGrow;
#ifdef DEBUG_ARRAYS
	char *m_sName;
#endif
};

#endif
