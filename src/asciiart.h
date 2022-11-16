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


#ifndef ASCIIART_H
#define ASCIIART_H


// This must always be the first include directive
#include "config.h"

#include "xobject.h"
#include "xobarray.h"
#include "backtrace.h"


class CAsciiPiece : public CxObject
{
public:
	char GetAt(int x, int y);


	CAsciiPiece() {
		m_pBuf = NULL;
	}


	~CAsciiPiece() {
		if (m_pBuf != NULL) {
			delete[] m_pBuf;
			m_pBuf = NULL;
		}
	}


	unsigned short m_iWidth;
	unsigned short m_iHeight;
	char *m_pBuf;
};


class CAsciiArt : public CxObject
{
public:
	void NewLine();
	char GetChar();
	int GetWidth();
	int GetHeight();
	int m_iPosX;
	int m_iPosY;
	void SelectOne();
	CAsciiArt();
	~CAsciiArt();
	void AddPiece(const char *s);
	int m_iSelectedPiece;
	CxObArray m_oaAsciiPieces;
};

#endif
