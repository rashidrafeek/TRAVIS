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


#ifndef MEMFILE_H
#define MEMFILE_H


// This must always be the first include directive
#include "config.h"

#include <stdio.h>
#include <string.h>
#include "xobject.h"


class CxMemFile : public CxObject
{
public:
	CxMemFile();
	~CxMemFile();
	explicit CxMemFile(CxMemFile *p);

	int fgets(char *buf, int len);

	bool ReadFileSuccessive(const char *s, int lines, bool verbose);

	bool SetSize(int i);

	void ReleaseBuffer();

	bool ReadFile(const char *s, bool text);

	int scanf(const char *s, void* data);

	void printf(const char *s, ...);

	void WriteFile(const char *s, bool text);

	void SkipEmpty();
	void SkipWord();
	void ReadWord(char *buf, int len);

	void SkipEmpty(const char *sep);
	void ReadWord(char *buf, int len, const char *sep);

	void ReverseSkipEmpty();

	bool Match(const char*s);

	void Seek(int pos);

	bool Eof();

	void Create(int length);

	unsigned char *m_pBuffer;
	unsigned char *m_pPointer;
	long m_iBufSize;
};


#endif
