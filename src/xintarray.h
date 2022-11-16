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


#ifndef XINTARRAY_H
#define XINTARRAY_H


// This must always be the first include directive
#include "config.h"

#include "tools.h"
#include "xobject.h"
#include "backtrace.h"


class CxIntArray : public CxObject
{
public:
	CxIntArray();
	~CxIntArray();
	explicit CxIntArray(const char *name);
	void SetName(const char *name);
	CxIntArray(const CxIntArray &o);
	void CopyFrom(CxIntArray *o);
	void Add(int f);
	void Append(CxIntArray *o);
	void SetAt(unsigned int pos, int f);
	void SetSize(unsigned int i);
	void SetMaxSize(unsigned int i);
	void SetGrow(unsigned int i);
	void RemoveAll();
	void RemoveAll_KeepSize();
	void RemoveAt(unsigned int pos, unsigned int count);
	void RemoveAt_KeepSize(unsigned int pos, unsigned int count);
	void InsertAt(int f, unsigned int pos);
	int Pop_KeepSize();


	bool Contains(int i)
	{
		BXIN;
		int z;
		for (z=0;z<(int)m_iSize;z++)
			if (m_pData[z] == i)
				return true;
		BXOUT;
		return false;
	}


	int GetPosition(int i)
	{
		BXIN;
		int z;
		for (z=0;z<(int)m_iSize;z++)
			if (m_pData[z] == i)
				return z;
		BXOUT;
		return -1;
	}


	int &GetAt(unsigned int i)
	{
		BXIN;
#ifdef DEBUG_CINTARRAY
		mprintf("@ CxIntArray::GetAt(int): %d...",i);
#endif
#ifdef DEBUG_ARRAYS
		if (i >= m_iSize)
		{
			if (m_sName != NULL)
				eprintf("CxIntArray \"%s\" Boundary Error (%d/%d).\n",m_sName,i,m_iSize);
					else eprintf("CxIntArray Boundary Error (%d/%d).\n",i,m_iSize);
			abort();
		}
#endif
#ifdef DEBUG_CINTARRAY
		mprintf("done.\n");
#endif
		BXOUT;
		return m_pData[i];
	}


	int &operator [] (unsigned int i)
	{
#ifdef DEBUG_CINTARRAY
		mprintf("@ CxIntArray::operator [] (int): %d\n",i);
#endif
		return GetAt(i);
	}


	int operator [] (unsigned int i) const
	{
		BXIN;
#ifdef DEBUG_CINTARRAY
		mprintf("@ const CxIntArray::operator [] (int): %d...",i);
#endif
#ifdef DEBUG_ARRAYS
		if (i >= m_iSize)
		{
			if (m_sName != NULL)
				eprintf("CxIntArray \"%s\" Boundary Error (%d/%d).\n",m_sName,i,m_iSize);
					else eprintf("CxIntArray Boundary Error (%d/%d).\n",i,m_iSize);
			abort();
		}
#endif
#ifdef DEBUG_CINTARRAY
		mprintf("done.\n");
#endif
		BXOUT;
		return m_pData[i];
	}


	int GetAt(unsigned int i) const
	{
		BXIN;
#ifdef DEBUG_CINTARRAY
		mprintf("@ const CxIntArray::GetAt(int): %d...",i);
#endif
#ifdef DEBUG_ARRAYS
		if (i >= m_iSize)
		{
			if (m_sName != NULL)
				eprintf("CxIntArray \"%s\" Boundary Error (%d/%d).\n",m_sName,i,m_iSize);
					else eprintf("CxIntArray Boundary Error (%d/%d).\n",i,m_iSize);
			abort();
		}
#endif
#ifdef DEBUG_CINTARRAY
		mprintf("done.\n");
#endif
		BXOUT;
		return m_pData[i];
	}


	int GetSize() const
	{
	#ifdef DEBUG_CINTARRAY
		mprintf("@ CxIntArray::GetSize(): %d\n",m_iSize);
	#endif
		return m_iSize;
	}	

	
private:	
	int *m_pData;
	unsigned int m_iSize;
	unsigned int m_iMaxSize;
	unsigned int m_iGrow;
#ifdef DEBUG_ARRAYS
	char *m_sName;
#endif
};

#endif
