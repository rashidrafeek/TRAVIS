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

#include "xdmatrixmn.h"
#include "tools.h"


const char *GetRevisionInfo_xdmatrixmn(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_xdmatrixmn() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


void CxDMatrixMN::Dump() const
{
	BXIN;
	int z, z2;
	for (z=0;z<m_iRows;z++)
	{
		for (z2=0;z2<m_iCols;z2++)
			if (fabs(GetAt(z,z2)) < 0.000005)
				mprintf(" %11.5f",0.0);
			else
				mprintf(" %11.5f",GetAt(z,z2));
		mprintf("\n");
	}
	BXOUT;
}


void CxDMatrixMN::Dump(FILE *a) const
{
	BXIN;
	int z, z2;
	for (z=0;z<m_iRows;z++)
	{
		for (z2=0;z2<m_iCols;z2++)
			if (fabs(GetAt(z,z2)) < 0.000005)
				mfprintf(a," %11.5f",0.0);
			else
				mfprintf(a," %11.5f",GetAt(z,z2));
		mfprintf(a,"\n");
	}
	BXOUT;
}


void CxDMatrixMN::DumpSmall() const
{
	BXIN;
	int z, z2;
	for (z=0;z<m_iRows;z++)
	{
		for (z2=0;z2<m_iCols;z2++)
			if (fabs(GetAt(z,z2)) < 0.00000005)
				mprintf(" %10.7f",0.0);
			else
				mprintf(" %10.7f",GetAt(z,z2));
		mprintf("\n");
	}
	BXOUT;
}


void CxDMatrixMN::DumpSmall(FILE *a) const
{
	BXIN;
	int z, z2;
	for (z=0;z<m_iRows;z++)
	{
		for (z2=0;z2<m_iCols;z2++)
			if (fabs(GetAt(z,z2)) < 0.00000005)
				mfprintf(a," %10.7f",0.0);
			else
				mfprintf(a," %10.7f",GetAt(z,z2));
		mfprintf(a,"\n");
	}
	BXOUT;
}


