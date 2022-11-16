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

#include "spectrum.h"
#include <math.h>
#include "travis.h"


const char *GetRevisionInfo_spectrum(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_spectrum() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


CSpectrum::CSpectrum()
{
	BTIN;
	m_pData = NULL;
	m_pComplex = NULL;
	BTOUT;
}

CSpectrum::~CSpectrum()
{
	BTIN;
	if (m_pData != NULL)
	{
		delete[] m_pData;
		m_pData = NULL;
	}
	if (m_pComplex != NULL)
	{
		delete[] m_pComplex;
		m_pComplex = NULL;
	}
	BTOUT;
}


void CSpectrum::Create(int i)
{
	m_iSize = i;

	try { m_pData = new double[i]; } catch(...) { m_pData = NULL; }
	if (m_pData == NULL) NewException((double)i*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { m_pComplex = new double[i*2]; } catch(...) { m_pComplex = NULL; }
	if (m_pComplex == NULL) NewException((double)i*2*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
}


void CSpectrum::FromComplex(double *f)
{
	int z;

	if (m_pComplex != NULL)
		memcpy(m_pComplex,f,m_iSize*2*sizeof(double));

	for (z=0;z<m_iSize;z++)
		m_pData[z] = f[z*2];
//		m_pData[z] = (double)sqrt(f[z*2]*f[z*2]+f[z*2+1]*f[z*2+1]);
}


void CSpectrum::Write(const char *pre, const char *s, const char *post)
{
	BTIN;
	FILE *a;
	int z;
//	char buf[256];
	CxString buf;

//	sprintf(buf,"%s%s%s",pre,s,post);
	buf.sprintf("%s%s%s",pre,s,post);

	a = OpenFileWrite(buf,true);

/*	mfprintf(a,"# Wave Number / cm^-1; Abs; Abs^2");
	if (m_pComplex != NULL)
		mfprintf(a,"; Realteil; Imaginaerteil");
	mfprintf(a,"\n");*/
//	mfprintf(a,"# Size=%d, MaxRWL=%f, WN=%f\n",m_iSize,m_fMaxRWL,m_fWaveNumber);

	mfprintf(a,"# wave number / cm^-1;  intensity;  intensity^2\n");

	for (z=0;z<m_iSize;z++)
	{
		if ((z*m_fMaxRWL/m_iSize > m_fWaveNumber) && (m_fWaveNumber > 0))
			break;
/*		mfprintf(a,"%.2f;  %.5f; %.5f",z*m_fMaxRWL/m_iSize,m_pData[z],m_pData[z]*m_pData[z]);
		if (m_pComplex != NULL)
			mfprintf(a,";  %.5f;  %.5f",m_pComplex[z*2],m_pComplex[z*2+1]);
		mfprintf(a,"\n");*/
		mfprintf(a,"%G;  %G;  %G\n",z*m_fMaxRWL/m_iSize,m_pData[z],m_pData[z]*m_pData[z]);
	}

	fclose(a);
	BTOUT;
}


void CSpectrum::SetMaxRWL(double f)
{
	m_fMaxRWL = f;
}


void CSpectrum::MakeDB()
{
	double f;
	int z;

	f = 0;
	for (z=0;z<m_iSize;z++)
		f += m_pData[z]*m_pData[z];

	for (z=0;z<m_iSize;z++)
		m_pData[z] = (10.0 * log10(m_pData[z]*m_pData[z] / f));
}


void CSpectrum::Multiply(double f)
{
	int z;

	for (z=0;z<m_iSize;z++)
		m_pData[z] *= f;
}


void CSpectrum::SetMaxValue(double f)
{
	int z;
	double m;

	m = 0;
	for (z=0;z<m_iSize;z++)
		if (m_pData[z] > m)
			m = m_pData[z];
	m = f/m;
	for (z=0;z<m_iSize;z++)
		m_pData[z] *= m;
}


void CSpectrum::SetIntegral(double f)
{
	int z;
	double m;

	m = 0;
	for (z=0;z<m_iSize;z++)
		m += m_pData[z];

	m *= m_fMaxRWL / m_iSize / f;

	for (z=0;z<m_iSize;z++)
		m_pData[z] /= m;
}
