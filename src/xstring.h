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


#ifndef XSTRING_H
#define XSTRING_H


// This must always be the first include directive
#include "config.h"

#include "tools.h"
#include "xobject.h"
#include "backtrace.h"


#define XVSPRINTF_LINUX(obj, format, pre)                   \
{                                                           \
	va_list XXXparams;                                      \
	int XXXi;                                               \
                                                            \
	va_start(XXXparams, pre);                               \
	XXXi = obj.Format_Internal(format,0,XXXparams);         \
	va_end(XXXparams);                                      \
                                                            \
	if (XXXi != 0)                                          \
	{                                                       \
		va_start(XXXparams, pre);                           \
		obj.Format_Internal(format,XXXi,XXXparams);         \
		va_end(XXXparams);                                  \
	}                                                       \
}


#define XVSPRINTF_WINDOWS(obj, format, pre)                 \
{                                                           \
	va_list XXXparams;                                      \
	int XXXi;                                               \
                                                            \
	XXXi = 0;                                               \
	do {                                                    \
		va_start(XXXparams, pre);                           \
		XXXi = obj.Format_Internal(format,XXXi,XXXparams);  \
		va_end(XXXparams);                                  \
	} while (XXXi != 0);                                    \
}



class CxString : public CxObject
{
public:
	//CxString();
	//~CxString();
	//CxString(const char *s);
	//CxString(const CxString &s);
	CxString(const CxString &s1, const CxString &s2);
	//operator const char*() const;
	//CxString & operator = (const CxString &s);
	CxString operator + (const CxString &s) const;
	bool operator == (const CxString &s) const;
	bool operator != (const CxString &s) const;
	bool operator > (const CxString &s) const;
	bool operator < (const CxString &s) const;
	bool operator >= (const CxString &s) const;
	bool operator <= (const CxString &s) const;
	void operator += (const CxString &s);
//	char& operator [] (int i);
	//char& operator () (int i);
	//int GetLength() const;
	int FindFirst(char c) const;
	int FindNext(int i, char c) const;
	int FindLast(char c) const;
	CxString Mid(int pos, int count) const;
	CxString Mid(int pos) const;
	CxString Left(int count) const;

#ifdef __GNUG__ // Variadic Argument Type Checking of GCC

	void Format(const char *s, ...)   __attribute__ ((format (printf, 2, 3)));
	void sprintf(const char *s, ...)  __attribute__ ((format (printf, 2, 3)));

#else

	void Format(const char *s, ...);
	void sprintf(const char *s, ...);

#endif

//	void vsprintf(const char *s, va_list params);
	void strcat(const char *s);
	void strcpy(const char *s);
	void memcpy(const char *s, int size);
	void Dump();
	void SetBufSize(int i);
	const char* fgets(int len, FILE *a);
	void ToLowercase();
//	void ToUppercase();
	char* GetWritePointer();

	int Format_Internal(const char *s, int length, va_list params);


	CxString()
	{
		BXIN;
	#ifdef DEBUG_CxString
		mprintf("@ CxString::CxString()\n");
	#endif

		try { m_pData = new char[1]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)1*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		m_pData[0] = 0;
		m_iBufLen = 1;

		BXOUT;
	}

		
	~CxString()
	{
		BXIN;
	#ifdef DEBUG_CxString
		mprintf("@ CxString::~CxString(): \"%s\"...",m_pData);
	#endif
		if (m_pData != NULL)
		{
			delete[] m_pData;
			m_pData = NULL;
		}
		m_iBufLen = 0;
	#ifdef DEBUG_CxString
		mprintf("done.\n");
	#endif
		BXOUT;
	}

		
	CxString(const char *s)
	{
		BXIN;
	#ifdef DEBUG_CxString
		mprintf("@ CxString::CxString(const char *): \"%s\"...",s);
	#endif

		int i;

/*		if (s == NULL)
		{
			eprintf("CxString::CxString(const char *): NULL pointer passed.\n");
			abort();
		}

		i = (int)strlen(s);*/

		if (s != NULL) {

			i = (int)strlen(s);

			try { m_pData = new char[i+1]; } catch(...) { m_pData = NULL; }
			if (m_pData == NULL) NewException((double)(i+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			::strcpy(m_pData,s);

		} else {

			i = 0;

			try { m_pData = new char[i+1]; } catch(...) { m_pData = NULL; }
			if (m_pData == NULL) NewException((double)(i+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			m_pData[0] = 0;
		}

		m_iBufLen = i+1;

	#ifdef DEBUG_CxString
		mprintf("done.\n");
	#endif
		BXOUT;
	}

		
	CxString(const CxString &s) : CxObject()
	{
		BXIN;
	#ifdef DEBUG_CxString
		mprintf("@ CxString::CxString(const CxString &): \"%s\"...",(const char*)s);
	#endif

		int i;

		i = s.GetLength();

		try { m_pData = new char[i+1]; } catch(...) { m_pData = NULL; }
		if (m_pData == NULL) NewException((double)(i+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		::strcpy(m_pData,s.m_pData);

		m_iBufLen = i+1;

	#ifdef DEBUG_CxString
		mprintf("done.\n");
	#endif
		BXOUT;
	}


	CxString & operator = (const CxString &s)
	{
		int i;

		i = s.GetLength();

		if (i+1 > m_iBufLen)
		{
			if (m_pData != NULL)
				delete[] m_pData;

			try { m_pData = new char[i+1]; } catch(...) { m_pData = NULL; }
			if (m_pData == NULL) NewException((double)(i+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			::strcpy(m_pData,s.m_pData);

			m_iBufLen = i+1;
		} else
			::strcpy(m_pData,s.m_pData);

		return *this;
	}


	operator const char*() const
	{
	#ifdef DEBUG_CxString
		mprintf("@ CxString::operator const char* () const: \"%s\"\n",m_pData);
	#endif
		return m_pData;
	}

		
	char& operator () (int i)
	{
	#ifdef DEBUG_CxString
		mprintf("@ CxString::operator () (int i): \"%s\", %d\n",m_pData,i);
	#endif
		if ((i < 0) || (i >= m_iBufLen))
		{
			eprintf("Error: CxString::operator (): Boundary Error (%d/%d).\n",i,m_iBufLen);
			abort();
		}
		return m_pData[i];
	}

  
	char operator () (int i) const
	{
	#ifdef DEBUG_CxString
		mprintf("@ CxString::operator () (int i): \"%s\", %d\n",m_pData,i);
	#endif
		if ((i < 0) || (i >= m_iBufLen))
		{
			eprintf("Error: CxString::operator (): Boundary Error (%d/%d).\n",i,m_iBufLen);
			abort();
		}
		return m_pData[i];
	}

  
	int GetLength() const
	{
	#ifdef DEBUG_CxString
		mprintf("@ CxString::GetLength() const: \"%s\", %d\n",m_pData,strlen(m_pData));
	#endif
		if (m_pData != NULL)
			return (int)strlen(m_pData);
		else
			return 0;
	}

	

	
private:

	char *m_pData;
	int m_iBufLen;
};


CxString operator + (const char *s1, const CxString& s2);


#endif
