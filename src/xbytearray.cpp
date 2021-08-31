/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

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


// This must always be the first include directive
#include "config.h"

#include "xbytearray.h"


const char *GetRevisionInfo_xbytearray(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_xbytearray() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


#ifndef DEBUG_ARRAYS
#define m_sName (NULL)
#endif


CxByteArray::CxByteArray()
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::CxByteArray()\n");
#endif
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
	m_iGrow = 16;
#ifdef DEBUG_ARRAYS
	m_sName = NULL;
#endif
	BXOUT;
}


CxByteArray::CxByteArray(const char *name)
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::CxByteArray(const char *s)\n");
#endif
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
	m_iGrow = 16;
#ifdef DEBUG_ARRAYS
	m_sName = NULL;
#endif
	SetName(name);
	BXOUT;
}

	
CxByteArray::~CxByteArray()
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::~CxByteArray()\n");
#endif
	RemoveAll();
#ifdef DEBUG_ARRAYS
	if (m_sName != NULL)
	{
		delete[] m_sName;
		m_sName = NULL;
	}
#endif
	BXOUT;
}

	
CxByteArray::CxByteArray(const CxByteArray &o) : CxObject()
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::CxByteArray(CxByteArray &)...");
#endif
	unsigned long z;
	m_iGrow = o.m_iGrow;

	if ((o.m_iMaxSize != 0) && (o.m_pData != NULL)) {

		m_iSize = o.m_iSize;
		m_iMaxSize = o.m_iMaxSize;

		try { m_pData = new unsigned char[m_iMaxSize]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(unsigned char),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
		
		for (z=0;z<m_iSize;z++)
			m_pData[z] = o.m_pData[z];

	} else {
		m_iSize = 0;
		m_iMaxSize = 0;
		m_pData = NULL;
	}

#ifdef DEBUG_ARRAYS
	if (o.m_sName != NULL)
		SetName(o.m_sName);
	else
		m_sName = NULL;
#endif

#ifdef DEBUG_CBYTEARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxByteArray::CopyFrom(const CxByteArray *o)
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::CopyFrom(CxByteArray*)...");
#endif
	if (m_iMaxSize != o->m_iMaxSize)
	{
		m_iMaxSize = o->m_iMaxSize;
		if (m_pData != NULL)
			delete[] m_pData;

		try { m_pData = new unsigned char[m_iMaxSize]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(unsigned char),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	}
	m_iSize = o->m_iSize;
	m_iGrow = o->m_iGrow;
	if (o->m_pData != NULL)
		memcpy(m_pData,o->m_pData,m_iSize*sizeof(unsigned char));
#ifdef DEBUG_CBYTEARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}	


void CxByteArray::Append(CxByteArray *o)
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::Append(CxByteArray*)...");
#endif
	if (m_iMaxSize != o->m_iMaxSize)
	{
		m_iMaxSize = o->m_iMaxSize;
		if (m_pData != NULL)
			delete[] m_pData;

		try { m_pData = new unsigned char[m_iMaxSize]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(unsigned char),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	}
	m_iGrow = o->m_iGrow;
	m_iSize = o->m_iSize;
	if (o->m_pData != NULL)
		memcpy(m_pData,o->m_pData,m_iSize*sizeof(unsigned char));
#ifdef DEBUG_CBYTEARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}	


void CxByteArray::SetAt(unsigned long pos, unsigned char f)
{
	BXIN;
	unsigned long z;
#ifdef DEBUG_CBYTEARRAY
	bool s = false;
	mprintf("@ CxByteArray::Add(unsigned short)...");
#endif
	if (pos >= m_iMaxSize)
	{
#ifdef DEBUG_CBYTEARRAY
		s = true;
		mprintf("\n");
#endif
		SetMaxSize(pos+1);
	}
	m_pData[pos] = f;
	if (pos >= m_iSize)
	{
		for (z=m_iSize;z<pos;z++)
			m_pData[z] = 0;
		m_iSize = pos+1;
	}
#ifdef DEBUG_CBYTEARRAY
	if (s)
		mprintf("@ done.\n");
	else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxByteArray::Add(unsigned char f)
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	bool s = false;
	mprintf("@ CxByteArray::Add(unsigned short)...");
#endif
	if (m_iSize+1 > m_iMaxSize)
	{
#ifdef DEBUG_CBYTEARRAY
		s = true;
		mprintf("\n");
#endif
//		SetMaxSize(m_iMaxSize + m_iGrow);
		if (m_iMaxSize == 0)
			SetMaxSize(m_iMaxSize + m_iGrow);
				else SetMaxSize(m_iMaxSize*2);
	}
	m_pData[m_iSize] = f;
	m_iSize++;
#ifdef DEBUG_CBYTEARRAY
	if (s)
		mprintf("@ done.\n");
			else mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxByteArray::SetSize(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::SetSize(int): %d...",i);
#endif
	if (m_iSize == i)
		return;
	unsigned char *temp;

	try { temp = new unsigned char[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(unsigned char),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		if (m_iSize > i)
			memcpy(temp,m_pData,i*sizeof(unsigned char));
		else
			memcpy(temp,m_pData,m_iSize*sizeof(unsigned char));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
	m_iSize = i;
#ifdef DEBUG_CBYTEARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxByteArray::SetMaxSize(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::SetMaxSize(int): %d...",i);
#endif
	unsigned char *temp;

	try { temp = new unsigned char[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(unsigned char),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,m_iSize*sizeof(unsigned char));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
#ifdef DEBUG_CBYTEARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxByteArray::SetGrow(unsigned long i)
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::SetGrow(int): %d\n",i);
#endif
	m_iGrow = i;
	BXOUT;
}
		
	
void CxByteArray::RemoveAll()
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::RemoveAll():...");
#endif
	if (m_pData != NULL)
		delete[] m_pData;
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
#ifdef DEBUG_CBYTEARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxByteArray::RemoveAll_KeepSize()
{
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::RemoveAll_KeepSize():...");
#endif
	m_iSize = 0;
#ifdef DEBUG_CBYTEARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxByteArray::RemoveAt(unsigned long pos, unsigned long count)
{
	BXIN;
	unsigned char *temp;
		
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::RemoveAt(int, int): %d, %d...",pos,count);
#endif

	try { temp = new unsigned char[m_iSize-count]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize-count)*sizeof(unsigned char),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(unsigned char));
		memcpy(&temp[pos],&m_pData[pos+count],(m_iSize-pos-count)*sizeof(unsigned char));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize-=count;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CBYTEARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}

	
void CxByteArray::InsertAt(unsigned char f, unsigned long pos)
{
	BXIN;
	unsigned char *temp;
		
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::InsertAt(unsigned char, unsigned long): %d...",f);
#endif

	try { temp = new unsigned char[m_iSize+1]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize+1)*sizeof(unsigned char),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(unsigned char));
		memcpy(&temp[pos+1],&m_pData[pos],(m_iSize-pos)*sizeof(unsigned char));
		delete[] m_pData;
	}
	temp[pos] = f;
	m_pData = temp;
	m_iSize++;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_CBYTEARRAY
	mprintf("done.\n");
#endif
	BXOUT;
}


void CxByteArray::SetName(const char *name)
{
	UNUSED(name);
#ifdef DEBUG_ARRAYS
	BXIN;
#ifdef DEBUG_CBYTEARRAY
	mprintf("@ CxByteArray::SetName(const char *): \"%s\"...",s);
#endif
	if (m_sName != NULL)
		delete[] m_sName;
	try { m_sName = new char[strlen(name)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sName,name);
#ifdef DEBUG_CBYTEARRAY
	mprintf("done.\n");
#endif
	BXOUT;
#endif
}
