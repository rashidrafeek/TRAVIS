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


// This must always be the first include directive
#include "config.h"

#include "xobarray.h"


const char *GetRevisionInfo_xobarray(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_xobarray() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


#ifndef DEBUG_ARRAYS
#define m_sName (NULL)
#endif


CxObArray::CxObArray()
{
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::CxObArray()\n");
#endif
	m_pData = NULL;
#ifdef DEBUG_ARRAYS
	m_sName = NULL;
#endif
	m_iSize = 0;
	m_iMaxSize = 0;
	m_iGrow = 16;
}


CxObArray::CxObArray(const char *name)
{
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::CxObArray(char *name)\n");
#endif
	m_pData = NULL;
#ifdef DEBUG_ARRAYS
	m_sName = NULL;
#endif
	m_iSize = 0;
	m_iMaxSize = 0;
	m_iGrow = 16;
	SetName(name);
}

	
CxObArray::~CxObArray()
{
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::~CxObArray()\n");
#endif
	RemoveAll();
#ifdef DEBUG_ARRAYS
	if (m_sName != NULL)
	{
		delete[] m_sName;
		m_sName = NULL;
	}
#endif
}


CxObArray::CxObArray(const CxObArray &o) : CxObject()
{
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::CxObArray(CxObArray &)...");
#endif
	unsigned long z;
	m_iGrow = o.m_iGrow;

	if ((o.m_iMaxSize != 0) && (o.m_pData != NULL)) {

		m_iSize = o.m_iSize;
		m_iMaxSize = o.m_iMaxSize;

		try { m_pData = new CxObject*[m_iMaxSize]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)m_iMaxSize*sizeof(CxObject*),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
		
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

#ifdef DEBUG_COBARRAY
	mprintf("done.\n");
#endif
}

	
void CxObArray::SetName(const char *name)
{
	UNUSED(name);
#ifdef DEBUG_ARRAYS
	BXIN;
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::SetName(const char *): \"%s\"...",name);
#endif
	if (m_sName != NULL)
		delete[] m_sName;
	try { m_sName = new char[strlen(name)+1]; } catch(...) { m_sName = NULL; }
	if (m_sName == NULL) NewException((double)(strlen(name)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	strcpy(m_sName,name);
#ifdef DEBUG_COBARRAY
	mprintf("done.\n");
#endif
	BXOUT;
#endif
}

	
void CxObArray::SetAt(unsigned long pos, CxObject *o)
{
	unsigned long z;
#ifdef DEBUG_COBARRAY
	bool s = false;
	mprintf("@ CxObArray::Add(CObject *)...");
#endif
	if (pos >= m_iMaxSize)
	{
#ifdef DEBUG_COBARRAY
		s = true;
		mprintf("\n");
#endif
		SetMaxSize(pos+1);
	}
	m_pData[pos] = o;
	if (pos >= m_iSize)
	{
		for (z=m_iSize;z<pos;z++)
			m_pData[z] = NULL;
		m_iSize = pos+1;
	}
#ifdef DEBUG_COBARRAY
	if (s)
		mprintf("@ done.\n");
	else mprintf("done.\n");
#endif
}

	
void CxObArray::Add(CxObject *o)
{
#ifdef DEBUG_COBARRAY
	bool s = false;
	mprintf("@ CxObArray::Add(CObject *)...");
#endif
	if (m_iSize+1 > m_iMaxSize)
	{
#ifdef DEBUG_COBARRAY
		s = true;
		mprintf("\n");
#endif
		if (m_iGrow < 1)
			m_iGrow = 1;
//		SetMaxSize(m_iMaxSize + m_iGrow);
		if (m_iMaxSize == 0)
			SetMaxSize(m_iMaxSize + m_iGrow);
				else SetMaxSize(m_iMaxSize*2);
	}
	m_pData[m_iSize] = o;
	m_iSize++;
#ifdef DEBUG_COBARRAY
	if (s)
		mprintf("@ done.\n");
			else mprintf("done.\n");
#endif
}

	
void CxObArray::SetSize(unsigned long i)
{
	unsigned long z;
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::SetSize(int): %d...",i);
#endif
	if (i == m_iSize)
		return;
	CxObject **temp;

	try { temp = new CxObject*[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(CxObject*),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		if (m_iSize > i)
			memcpy(temp,m_pData,i*sizeof(CxObject*));
		else
			memcpy(temp,m_pData,m_iSize*sizeof(CxObject*));
		delete[] m_pData;
	}
// Neu und heikel (Ballmer-Peak ^^)
	for (z=m_iSize;z<i;z++)
		temp[z] = NULL;
// Ende neu
	m_pData = temp;
	m_iMaxSize = i;
	m_iSize = i;
#ifdef DEBUG_COBARRAY
	mprintf("done.\n");
#endif
}


void CxObArray::SetSize_NoShrink(unsigned long i)
{
	unsigned long z;
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::SetSize_NoShrink(int): %d...",i);
#endif
	if (i == m_iSize)
		return;

	if (m_iSize < i)
	{
		CxObject **temp;

		try { temp = new CxObject*[i]; } catch(...) { temp = NULL; }
		if (temp == NULL) NewException((double)i*sizeof(CxObject*),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
		
		if (m_pData != NULL) {
			memcpy(temp,m_pData,m_iSize*sizeof(CxObject*));

			delete[] m_pData;
		}

		for (z=m_iSize;z<i;z++)
			temp[z] = NULL;

		m_pData = temp;
		m_iMaxSize = i;
	} else
	{
		m_iSize = i;
	}
#ifdef DEBUG_COBARRAY
	mprintf("done.\n");
#endif
}


void CxObArray::SetMaxSize(unsigned long i)
{
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::SetMaxSize(int): %d...",i);
#endif
	if (i == m_iMaxSize)
		return;
	CxObject **temp;

	try { temp = new CxObject*[i]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)i*sizeof(CxObject*),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,m_iSize*sizeof(CxObject*));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iMaxSize = i;
#ifdef DEBUG_COBARRAY
	mprintf("done.\n");
#endif
}

	
void CxObArray::SetGrow(unsigned long i)
{
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::SetGrow(int): %d\n",i);
#endif
	m_iGrow = i;
}
		
	
void CxObArray::RemoveAll()
{
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::RemoveAll():...");
#endif
	
	if (m_pData != NULL)
		delete[] m_pData;
	m_pData = NULL;
	m_iSize = 0;
	m_iMaxSize = 0;
#ifdef DEBUG_COBARRAY
	mprintf("done.\n");
#endif
}


void CxObArray::RemoveAll_KeepSize()
{
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::RemoveAll_KeepSize():...");
#endif
	m_iSize = 0;
#ifdef DEBUG_COBARRAY
	mprintf("done.\n");
#endif
}


void CxObArray::RemoveAt(unsigned long pos, unsigned long count)
{
	CxObject **temp;
		
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::RemoveAt(int, int): %d, %d...",pos,count);
#endif

	try { temp = new CxObject*[m_iSize-count]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize-count)*sizeof(CxObject*),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
//	for (z=pos;z<pos+count;z++)
//		delete m_pData[z];
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(CxObject*));
		memcpy(&temp[pos],&m_pData[pos+count],(m_iSize-pos-count)*sizeof(CxObject*));
		delete[] m_pData;
	}
	m_pData = temp;
	m_iSize-=count;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_COBARRAY
	mprintf("done.\n");
#endif
}


void CxObArray::RemoveAt_NoShrink(unsigned long pos, unsigned long count)
{
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::RemoveAt_NoShrink(int, int): %d, %d...",pos,count);
#endif
	memmove(&m_pData[pos],&m_pData[pos+count],(m_iSize-pos-count)*sizeof(CxObject*));
	m_iSize-=count;
#ifdef DEBUG_COBARRAY
	mprintf("done.\n");
#endif
}	


void CxObArray::InsertAt(CxObject *o, unsigned long pos)
{
	CxObject **temp;
		
#ifdef DEBUG_COBARRAY
	mprintf("@ CxObArray::InsertAt(CObject *, int): %d...");
#endif

	try { temp = new CxObject*[m_iSize+1]; } catch(...) { temp = NULL; }
	if (temp == NULL) NewException((double)(m_iSize+1)*sizeof(CxObject*),__FILE__,__LINE__,__PRETTY_FUNCTION__,m_sName);
	
	if (m_pData != NULL) {
		memcpy(temp,m_pData,pos*sizeof(CxObject*));
		memcpy(&temp[pos+1],&m_pData[pos],(m_iSize-pos)*sizeof(CxObject*));
		delete[] m_pData;
	}
	temp[pos] = o;
	m_pData = temp;
	m_iSize++;
	m_iMaxSize = m_iSize;
#ifdef DEBUG_COBARRAY
	mprintf("done.\n");
#endif
}

