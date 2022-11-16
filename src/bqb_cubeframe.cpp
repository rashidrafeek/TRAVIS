/***********************************************************************************

    LibBQB - File Format and Compression Algorithms for Trajectories of
             Volumetric Data and Atom Positions

    https://brehm-research.de/bqb

    Free software, licensed under GNU LGPL v3

    Copyright (c) Martin Brehm and Martin Thomas,
                  Martin Luther University Halle-Wittenberg, Germany,
                  2016 - 2022.

    Please cite:  M. Brehm, M. Thomas: "An Efficient Lossless Compression Algorithm
                  for Trajectories of Atom Positions and Volumetric Data",
                  J. Chem. Inf. Model. 2018, 58 (10), pp 2092-2107.

    --------------------------------------------------------------------------------

    LibBQB is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

***********************************************************************************/



// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"
#include <math.h>
#include "bqb_cubeframe.h"
#include "bqb_fastatof.h"
#include "bqb_interface.h"


const char *GetRevisionInfo_bqb_cubeframe(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_cubeframe() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



bool CBQBCubeFrame::ReadFrame(FILE *a, int eps, int csigni, int asigni, bool verbose) {

	int ac, z, z2, i, trunc;
	unsigned long ul;
	char buf[256], buf2[256], obuf[256], *p, *q;
	bool eol;
	double tf, mi, ma;
	CBQBAtom *ta;


	m_iEps = eps;
	m_iSigni = csigni;

	if (m_pAtoms != NULL)
		delete m_pAtoms;
	m_pAtoms = new CBQBAtomSet(m_IF);

	m_pAtoms->m_iSigni = asigni;
	m_pAtoms->m_bOrd = true;

	if (feof(a))
		return false;

	(void)!fgets(buf,256,a);

	if (feof(a))
		return false;

	buf[strlen(buf)-1] = 0;

	(void)!fgets(buf2,256,a);

	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Unexpected end of file (0).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}

	buf2[strlen(buf2)-1] = 0;
	m_pAtoms->m_sComment = new char[strlen(buf)+strlen(buf2)+2];

	strcpy(m_pAtoms->m_sComment,buf);
	strcat(m_pAtoms->m_sComment,"\n");
	strcat(m_pAtoms->m_sComment,buf2);

	(void)!fgets(buf,256,a);
	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Unexpected end of file (A).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	buf[strlen(buf)-1] = 0;
	strcpy(obuf,buf);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (A):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	ac = atoi(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (B):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	m_fCenter[0] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (C):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	m_fCenter[1] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	*q = 0;
	m_fCenter[2] = atof(p);

	(void)!fgets(buf,256,a);
	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Unexpected end of file (B).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	buf[strlen(buf)-1] = 0;
	strcpy(obuf,buf);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (D):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	m_iRes[0] = atoi(p);
	if (m_iRes[0] == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error reading X resolution:\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (E):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	m_fStrideA[0] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (F):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	m_fStrideA[1] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	*q = 0;
	m_fStrideA[2] = atof(p);

	(void)!fgets(buf,256,a);
	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Unexpected end of file (C).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	buf[strlen(buf)-1] = 0;
	strcpy(obuf,buf);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (G):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	m_iRes[1] = atoi(p);
	if (m_iRes[1] == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error reading Y resolution:\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (H):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	m_fStrideB[0] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (I):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	m_fStrideB[1] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	*q = 0;
	m_fStrideB[2] = atof(p);

	(void)!fgets(buf,256,a);
	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Unexpected end of file (D).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	buf[strlen(buf)-1] = 0;
	strcpy(obuf,buf);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (J):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	m_iRes[2] = atoi(p);
	if (m_iRes[2] == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error reading Z resolution:\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (K):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	m_fStrideC[0] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	if (*q == 0) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (L):\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}
	*q = 0;
	m_fStrideC[1] = atof(p);
	p = q+1;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	*q = 0;
	m_fStrideC[2] = atof(p);

	for (z=0;z<3;z++) {
		if (fabs(m_fStrideA[z]) > 2000.0) {
			m_IF.printf("CBQBCubeFrame::ReadFrame(): Error: m_fStrideA[%d] > 2000 a.u. (is %f a.u.).\n",z,m_fStrideA[z]);
			return false;
		}
		if (fabs(m_fStrideB[z]) > 2000.0) {
			m_IF.printf("CBQBCubeFrame::ReadFrame(): Error: m_fStrideB[%d] > 2000 a.u. (is %f a.u.).\n",z,m_fStrideB[z]);
			return false;
		}
		if (fabs(m_fStrideC[z]) > 2000.0) {
			m_IF.printf("CBQBCubeFrame::ReadFrame(): Error: m_fStrideC[%d] > 2000 a.u. (is %f a.u.).\n",z,m_fStrideC[z]);
			return false;
		}
		if (fabs(m_fCenter[z]) > 2000.0) {
			m_IF.printf("CBQBCubeFrame::ReadFrame(): Error: m_fCenter[%d] > 2000 a.u. (is %f a.u.).\n",z,m_fCenter[z]);
			return false;
		}
		m_iStrideA[z] = FloatToFixed(m_fStrideA[z],6);
		m_iStrideB[z] = FloatToFixed(m_fStrideB[z],6);
		m_iStrideC[z] = FloatToFixed(m_fStrideC[z],6);
		m_iCenter[z] = FloatToFixed(m_fCenter[z],6);
	}

	for (z=0;z<ac;z++) {
		ta = new CBQBAtom();
		m_pAtoms->m_oaAtoms.push_back(ta);
		(void)!fgets(buf,256,a);
		if (feof(a)) {
			m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Unexpected end of file (Atom %d/%d).\n",z+1,ac);
			m_IF.eprintf("The CUBE file is corrupted.\n");
			return false;
		}
		buf[strlen(buf)-1] = 0;
		strcpy(obuf,buf);
		p = buf;
		while (*p == ' ')
			p++;
		q = p;
		while ((*q != ' ') && (*q != 0))
			q++;
		if (*q == 0) {
			m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (Atom %d/%d A):\n",z+1,ac);
			m_IF.eprintf("    \"%s\"\n",obuf);
			m_IF.eprintf("The CUBE file is corrupted.\n");
			return false;
		}
		*q = 0;
		ta->m_iOrd = atoi(p);
		p = q+1;
		while (*p == ' ')
			p++;
		q = p;
		while ((*q != ' ') && (*q != 0))
			q++;
		if (*q == 0) {
			m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (Atom %d/%d B):\n",z+1,ac);
			m_IF.eprintf("    \"%s\"\n",obuf);
			m_IF.eprintf("The CUBE file is corrupted.\n");
			return false;
		}
		p = q+1;
		while (*p == ' ')
			p++;
		q = p;
		while ((*q != ' ') && (*q != 0))
			q++;
		if (*q == 0) {
			m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (Atom %d/%d C):\n",z+1,ac);
			m_IF.eprintf("    \"%s\"\n",obuf);
			m_IF.eprintf("The CUBE file is corrupted.\n");
			return false;
		}
		*q = 0;
		ta->m_fCoord[0] = atof(p);
		p = q+1;
		while (*p == ' ')
			p++;
		q = p;
		while ((*q != ' ') && (*q != 0))
			q++;
		if (*q == 0) {
			m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Incomplete line (Atom %d/%d D):\n",z+1,ac);
			m_IF.eprintf("    \"%s\"\n",obuf);
			m_IF.eprintf("The CUBE file is corrupted.\n");
			return false;
		}
		*q = 0;
		ta->m_fCoord[1] = atof(p);
		p = q+1;
		while (*p == ' ')
			p++;
		q = p;
		while ((*q != ' ') && (*q != 0))
			q++;
		*q = 0;
		ta->m_fCoord[2] = atof(p);

		for (z2=0;z2<3;z2++) {
			if (fabs(ta->m_fCoord[z2]) > 2000.0) {
				m_IF.printf("CBQBCubeFrame::ReadFrame(): Error: m_oaAtoms[%d][%d] > 2000 a.u. (is %f a.u.).\n",z,z2,ta->m_fCoord[z2]);
				return false;
			}
			m_pAtoms->m_faCOM[z2] += ta->m_fCoord[z2];
			ta->m_iCoord[z2] = FloatToFixed(ta->m_fCoord[z2],asigni);
	//		ta->m_fCoord[z2] = FixedToFloat(ta->m_iCoord[z2],asigni);
	//		ta->m_fCoord[z2] *= 0.52917721;
		}
//		m_IF.printf("%d  %f  %f  %f\n",ta->m_iOrd,ta->m_fCoord[0],ta->m_fCoord[1],ta->m_fCoord[2]);
	}

	for (z2=0;z2<3;z2++) {
		m_pAtoms->m_faCOM[z2] /= m_pAtoms->m_oaAtoms.size();
		m_pAtoms->m_iaCOM[z2] = FloatToFixed(m_pAtoms->m_faCOM[z2],asigni);
		m_pAtoms->m_faCOM[z2] = FixedToFloat(m_pAtoms->m_iaCOM[z2],asigni);
	}

	for (z=0;z<(int)m_pAtoms->m_oaAtoms.size();z++) {
		for (z2=0;z2<3;z2++) {
			//m_pAtoms->m_oaAtoms[z]->m_fRelCoord[z2] = m_pAtoms->m_oaAtoms[z]->m_fCoord[z2] - m_pAtoms->m_faCOM[z2];
			//m_pAtoms->m_oaAtoms[z]->m_iRelCoord[z2] = FloatToFixed(m_pAtoms->m_oaAtoms[z]->m_fRelCoord[z2],asigni);
			m_pAtoms->m_oaAtoms[z]->m_iRelCoord[z2] = m_pAtoms->m_oaAtoms[z]->m_iCoord[z2] - m_pAtoms->m_iaCOM[z2];
			m_pAtoms->m_oaAtoms[z]->m_fRelCoord[z2] = FixedToFloat(m_pAtoms->m_oaAtoms[z]->m_iRelCoord[z2],asigni);
		}
	}

/*	if (feof(a)) {
		m_IF.printf("Error: Unexpected end of cube file.\n");
		return false;
	}*/

/*	for (z=0;z<3;z++) {
		m_fCenter[z] *= 0.52917721;
		m_fStride[z] *= 0.52917721;
	}*/

	m_iResXY = m_iRes[0] * m_iRes[1];
	m_iResYZ = m_iRes[1] * m_iRes[2];
	m_iResXYZ = m_iRes[0] * m_iRes[1] * m_iRes[2];

	m_fMinVal[0] = m_fCenter[0];
	m_fMaxVal[0] = m_fCenter[0] + m_fStrideA[0] * m_iRes[0];
	m_fMinVal[1] = m_fCenter[1];
	m_fMaxVal[1] = m_fCenter[1] + m_fStrideB[1] * m_iRes[1];
	m_fMinVal[2] = m_fCenter[2];
	m_fMaxVal[2] = m_fCenter[2] + m_fStrideC[2] * m_iRes[2];

	if (verbose) {
		m_IF.printf("Res: %d %d %d\n",m_iRes[0],m_iRes[1],m_iRes[2]);
		m_IF.printf("Center: %f %f %f a.u.\n",m_fCenter[0],m_fCenter[1],m_fCenter[2]);
		m_IF.printf("Stride A: %f %f %f a.u.\n",m_fStrideA[0],m_fStrideA[1],m_fStrideA[2]);
		m_IF.printf("Stride B: %f %f %f a.u.\n",m_fStrideB[0],m_fStrideB[1],m_fStrideB[2]);
		m_IF.printf("Stride C: %f %f %f a.u.\n",m_fStrideC[0],m_fStrideC[1],m_fStrideC[2]);
		m_IF.printf("Range: X %f - %f, Y %f - %f, Z %f - %f a.u.\n",m_fMinVal[0],m_fMaxVal[0],m_fMinVal[1],m_fMaxVal[1],m_fMinVal[2],m_fMaxVal[2]);
		m_IF.printf("\n");
		m_IF.printf("    Reading cube file (resolution %d x %d x %d, %d atoms)...\n",m_iRes[0],m_iRes[1],m_iRes[2],ac);
		m_IF.printf("      [");
	}

	m_faBin.resize(m_iResXYZ);
	m_iaExpo.resize(m_iResXYZ);
	m_iaMantis.resize(m_iResXYZ);

	i = 0;
	mi = 1.0e30f;
	ma = -1.0e30f;
	trunc = 0;
	while (!feof(a))
	{
_read:
		(void)!fgets(buf,256,a);
		if (feof(a))
			break;
		buf[strlen(buf)-1] = 0;
		if ((ul = (unsigned long)strspn(buf,"0123456789eE.+- \r\n")) != strlen(buf)) {
			m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Data line contains invalid character at grid index %d/%d: 0x%02X (\"%c",i+1,m_iResXYZ,buf[ul],buf[ul]);
			m_IF.eprintf("\").\n");
			m_IF.eprintf("  Line: \"%s\"\n",buf);
			m_IF.eprintf("The CUBE file is corrupted.\n");
			return false;
		}
		strcpy(obuf,buf);
		p = buf;
		eol = false;

_next:
		while (*p == ' ')
			p++;
		if (*p == 0)
			goto _read;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
	//	if ((p-q) < 8)
	//		goto _read;

		if (*p == 0)
			eol = true;
		else
			*p = 0;

		#ifdef BQB_FAST_ATOF
			tf = bqb_fast_atof(q);
		#else
			tf = atof(q);
		#endif

/*		if (tf == 0) {
			if (strspn(q,"0123456789eE.+-") != strlen(q)) {
				m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error reading floating point number from \"%s\" at index %d/%d.\n",q,i+1,m_iResXYZ);
				m_IF.eprintf("  Line: \"%s\"\n",obuf);
				m_IF.eprintf("The CUBE file is corrupted.\n");
				return false;
			}
		}*/

		if (tf < mi)
			mi = tf;
		if (tf > ma)
			ma = tf;

		m_faBin[i] = tf;

		m_iaExpo[i] = (char)((int)floor(log10(fabs(tf)))-csigni+1);
		if (m_iaExpo[i] < -eps) {
			m_iaExpo[i] = (char)-eps;
			trunc++;
//			m_IF.printf("%s --> %15.10G --> %6d E %2d\n",q,tf,(int)(tf*pow(10,-m_pExpo[i])+0.5),m_pExpo[i]);
		}
		m_iaMantis[i] = (int)floor(tf*pow10(-m_iaExpo[i])+0.5);

		// MB Hack 26.01.2018
		m_faBin[i] = (double)m_iaMantis[i] * pow10(m_iaExpo[i]);

		i++;

		if (verbose) {
			if (fmod(i,m_iResXYZ/60.0) < 1.0) {
				m_IF.printf("#");
				fflush(stdout);
			}
		}

		if (i == m_iResXYZ)
			break;

		p++;

		if (eol)
			goto _read;
		else
			goto _next;
	}

	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::ReadFrame(): Error: Unexpected end of file (E).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}

	if (verbose) {
		m_IF.printf("]\n");
		m_IF.printf("      %d values truncated due to EPS.\n",trunc);
		m_IF.printf("      Value range: %10.6G ... %10.6G\n",mi,ma);
		m_IF.printf("\n");
	}

	return true;
}



void CBQBCubeFrame::FloatToInt() {

	unsigned long i;

	for (i=0;i<m_faBin.size();i++) {

		m_iaExpo[i] = (char)((int)floor(log10(fabs(m_faBin[i])))-m_iSigni+1);

		if (m_iaExpo[i] < -m_iEps)
			m_iaExpo[i] = (char)-m_iEps;

		m_iaMantis[i] = (int)floor(m_faBin[i]*pow10(-m_iaExpo[i])+0.5);

		m_faBin[i] = (double)m_iaMantis[i] * pow10(m_iaExpo[i]);
	}
}



bool CBQBCubeFrame::SkipFrame(FILE *a, bool verbose) const {

	UNUSED(verbose);

	int ac, z, i, resx, resy, resz;
	char buf[256], obuf[256], *p, *q;
	bool eol;


	if (feof(a))
		return false;

	(void)!fgets(buf,256,a);

	if (feof(a))
		return false;

	(void)!fgets(buf,256,a);

	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::SkipFrame(): Error: Unexpected end of file (A).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}

	(void)!fgets(buf,256,a);

	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::SkipFrame(): Error: Unexpected end of file (B).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}

	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	*q = 0;
	ac = atoi(p);

	(void)!fgets(buf,256,a);

	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::SkipFrame(): Error: Unexpected end of file (C).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}

	buf[strlen(buf)-1] = 0;
	strcpy(obuf,buf);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	*q = 0;
	resx = atoi(p);
	if (resx == 0) {
		m_IF.eprintf("CBQBCubeFrame::SkipFrame(): Error reading X resolution:\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}

	(void)!fgets(buf,256,a);

	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::SkipFrame(): Error: Unexpected end of file (D).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}

	buf[strlen(buf)-1] = 0;
	strcpy(obuf,buf);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	*q = 0;
	resy = atoi(p);
	if (resy == 0) {
		m_IF.eprintf("CBQBCubeFrame::SkipFrame(): Error reading Y resolution:\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}

	(void)!fgets(buf,256,a);

	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::SkipFrame(): Error: Unexpected end of file (E).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}

	buf[strlen(buf)-1] = 0;
	strcpy(obuf,buf);
	p = buf;
	while (*p == ' ')
		p++;
	q = p;
	while ((*q != ' ') && (*q != 0))
		q++;
	*q = 0;
	resz = atoi(p);
	if (resz == 0) {
		m_IF.eprintf("CBQBCubeFrame::SkipFrame(): Error reading Z resolution:\n");
		m_IF.eprintf("    \"%s\"\n",obuf);
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}

	for (z=0;z<ac;z++) {

		(void)!fgets(buf,256,a);

		if (feof(a)) {
			m_IF.eprintf("CBQBCubeFrame::SkipFrame(): Error: Unexpected end of file (Atom %d/%d).\n",z+1,ac);
			m_IF.eprintf("The CUBE file is corrupted.\n");
			return false;
		}
	}


	i = 0;
	while (!feof(a)) {
_read:
		(void)!fgets(buf,256,a);
		if (feof(a))
			break;
		p = buf;
		eol = false;

_next:
		while (*p == ' ')
			p++;
		q = p;
		while ((*p != ' ') && (*p != 0))
			p++;
		if (*p == 0)
			eol = true;
		if ((p-q) < 8)
			goto _read;

		i++;

		if (i == resx*resy*resz)
			break;

		p++;
		if (eol)
			goto _read;
		else
			goto _next;
	}

	if (feof(a)) {
		m_IF.eprintf("CBQBCubeFrame::SkipFrame(): Error: Unexpected end of file (F).\n");
		m_IF.eprintf("The CUBE file is corrupted.\n");
		return false;
	}

	return true;
}


void CBQBCubeFrame::WriteFrame(FILE *a, bool verbose) {

	UNUSED(verbose);

	int z, ix, iy, iz, i;

	if (m_pAtoms->m_sComment != NULL)
		fprintf(a,"%s\n",m_pAtoms->m_sComment);
	else
		fprintf(a,"(Comment 1)\n(Comment 2)\n");

	fprintf(a,"%5lu %11.6f %11.6f %11.6f\n",(unsigned long)m_pAtoms->m_oaAtoms.size(),m_fCenter[0],m_fCenter[1],m_fCenter[2]);
	fprintf(a,"%5d %11.6f %11.6f %11.6f\n",m_iRes[0],m_fStrideA[0],m_fStrideA[1],m_fStrideA[2]);
	fprintf(a,"%5d %11.6f %11.6f %11.6f\n",m_iRes[1],m_fStrideB[0],m_fStrideB[1],m_fStrideB[2]);
	fprintf(a,"%5d %11.6f %11.6f %11.6f\n",m_iRes[2],m_fStrideC[0],m_fStrideC[1],m_fStrideC[2]);

	for (z=0;z<(int)m_pAtoms->m_oaAtoms.size();z++)
		fprintf(a,"%5d %11.6f %11.6f %11.6f %11.6f\n",m_pAtoms->m_oaAtoms[z]->m_iOrd,0.0,m_pAtoms->m_oaAtoms[z]->m_fCoord[0],m_pAtoms->m_oaAtoms[z]->m_fCoord[1],m_pAtoms->m_oaAtoms[z]->m_fCoord[2]);

	i = 0;
	for (ix=0;ix<m_iRes[0];ix++) {
		for (iy=0;iy<m_iRes[1];iy++) {
			for (iz=0;iz<m_iRes[2];iz++) {
//				fm_IF.printf(a,"  %11.6E",m_faBin[ix*m_iResYZ+iy*m_iRes[2]+iz]);
				fprintf_expo(a,m_iaMantis[ix*m_iResYZ+iy*m_iRes[2]+iz],m_iSigni,m_iaExpo[ix*m_iResYZ+iy*m_iRes[2]+iz]);
				i++;
				if (i == 6) {
					i = 0;
					fprintf(a,"\n");
				}
			}
			if (i != 0) {
				fprintf(a,"\n");
				i = 0;
			}
		}
	}
}


void CBQBCubeFrame::WriteFrame_Double(FILE *a, bool verbose) {

	UNUSED(verbose);

	int z, ix, iy, iz, i;

	if (m_pAtoms->m_sComment != NULL)
		fprintf(a,"%s\n",m_pAtoms->m_sComment);
	else
		fprintf(a,"(Comment 1)\n(Comment 2)\n");

	fprintf(a,"%5lu %11.6f %11.6f %11.6f\n",(unsigned long)m_pAtoms->m_oaAtoms.size(),m_fCenter[0],m_fCenter[1],m_fCenter[2]);
	fprintf(a,"%5d %11.6f %11.6f %11.6f\n",m_iRes[0],m_fStrideA[0],m_fStrideA[1],m_fStrideA[2]);
	fprintf(a,"%5d %11.6f %11.6f %11.6f\n",m_iRes[1],m_fStrideB[0],m_fStrideB[1],m_fStrideB[2]);
	fprintf(a,"%5d %11.6f %11.6f %11.6f\n",m_iRes[2],m_fStrideC[0],m_fStrideC[1],m_fStrideC[2]);

	for (z=0;z<(int)m_pAtoms->m_oaAtoms.size();z++)
		fprintf(a,"%5d %11.6f %11.6f %11.6f %11.6f\n",m_pAtoms->m_oaAtoms[z]->m_iOrd,0.0,m_pAtoms->m_oaAtoms[z]->m_fCoord[0],m_pAtoms->m_oaAtoms[z]->m_fCoord[1],m_pAtoms->m_oaAtoms[z]->m_fCoord[2]);

	i = 0;
	for (ix=0;ix<m_iRes[0];ix++) {
		for (iy=0;iy<m_iRes[1];iy++) {
			for (iz=0;iz<m_iRes[2];iz++) {
				fprintf(a,"  %11.6E",m_faBin[ix*m_iResYZ+iy*m_iRes[2]+iz]);
//				fprintf_expo(a,m_iaMantis[ix*m_iResYZ+iy*m_iRes[2]+iz],m_iSigni,m_iaExpo[ix*m_iResYZ+iy*m_iRes[2]+iz]);
				i++;
				if (i == 6) {
					i = 0;
					fprintf(a,"\n");
				}
			}
			if (i != 0) {
				fprintf(a,"\n");
				i = 0;
			}
		}
	}
}


void CBQBAtomSet::WriteXYZ(FILE *a, int signi) {

	int z, l;

	fprintf(a,"  %lu\n",(unsigned long)m_oaAtoms.size());
	if (m_sComment != NULL)
		fprintf(a,"%s\n",m_sComment);
	else
		fprintf(a,"\n");

	l = signi+4;

	for (z=0;z<(int)m_oaAtoms.size();z++)
		fprintf(a,"%2s %*.*f %*.*f %*.*f\n",m_oaAtoms[z]->m_sLabel.c_str(),l,signi,m_oaAtoms[z]->m_fCoord[0],l,signi,m_oaAtoms[z]->m_fCoord[1],l,signi,m_oaAtoms[z]->m_fCoord[2]);
}


bool CBQBAtomSet::ReadXYZ(FILE *a, int signi, FILE *ref) {

	char buf[256], obuf[256], *p, *q;
	int i, z, z2;
	CBQBAtom *at;

	m_iSigni = signi;

	m_bLabels = true;

	for (z=0;z<(int)m_oaAtoms.size();z++)
		delete m_oaAtoms[z];
	m_oaAtoms.clear();

	(void)!fgets(buf,256,a);
	buf[strlen(buf)-1] = 0;
	i = atoi(buf);

	if (feof(a))
		return false;

	(void)!fgets(buf,256,a);

	if (feof(a)) {
		m_IF.eprintf("CBQBAtomSet::ReadXYZ(): Error: Unexpected end of file (A).\n");
		m_IF.eprintf("The XYZ file is corrupted.\n");
		return false;
	}

	buf[strlen(buf)-1] = 0;
	m_sComment = new char[strlen(buf)+1];
	strcpy(m_sComment,buf);

	for (z=0;z<3;z++)
		m_faCOM[z] = 0;

	for (z=0;z<i;z++) {

		at = new CBQBAtom();
		m_oaAtoms.push_back(at);

		(void)!fgets(buf,256,a);

		if (feof(a)) {
			m_IF.eprintf("CBQBAtomSet::ReadXYZ(): Error: Unexpected end of file (Atom %d/%d).\n",z+1,i);
			m_IF.eprintf("The XYZ file is corrupted.\n");
			return false;
		}

		buf[strlen(buf)-1] = 0;
		strcpy(obuf,buf);

		p = buf;

		while (*p == ' ')
			p++;
		q = p;
		while ((*q != ' ') && (*q != 0))
			q++;
		if (*q == 0) {
			m_IF.eprintf("CBQBAtomSet::ReadXYZ(): Error: Incomplete line (Atom %d/%d A):\n",z+1,i);
			m_IF.eprintf("    \"%s\"\n",obuf);
			m_IF.eprintf("The XYZ file is corrupted.\n");
			return false;
		}
		*q = 0;

		at->m_sLabel = p;
		at->m_iOrd = GetAtomOrd(p);

		p = q+1;
		while (*p == ' ')
			p++;
		q = p;
		while ((*q != ' ') && (*q != 0))
			q++;
		if (*q == 0) {
			m_IF.eprintf("CBQBAtomSet::ReadXYZ(): Error: Incomplete line (Atom %d/%d B):\n",z+1,i);
			m_IF.eprintf("    \"%s\"\n",obuf);
			m_IF.eprintf("The XYZ file is corrupted.\n");
			return false;
		}
		*q = 0;

		//at->m_fCoord[0] = atof(p);
		#ifdef BQB_FAST_ATOF
			at->m_fCoord[0] = bqb_fast_atof(p);
		#else
			at->m_fCoord[0] = atof(p);
		#endif

		p = q+1;
		while (*p == ' ')
			p++;
		q = p;
		while ((*q != ' ') && (*q != 0))
			q++;
		if (*q == 0) {
			m_IF.eprintf("CBQBAtomSet::ReadXYZ(): Error: Incomplete line (Atom %d/%d C):\n",z+1,i);
			m_IF.eprintf("    \"%s\"\n",obuf);
			m_IF.eprintf("The XYZ file is corrupted.\n");
			return false;
		}
		*q = 0;

		//at->m_fCoord[1] = atof(p);
		#ifdef BQB_FAST_ATOF
			at->m_fCoord[1] = bqb_fast_atof(p);
		#else
			at->m_fCoord[1] = atof(p);
		#endif

		p = q+1;
		while (*p == ' ')
			p++;

		//at->m_fCoord[2] = atof(p);
		#ifdef BQB_FAST_ATOF
			at->m_fCoord[2] = bqb_fast_atof(p);
		#else
			at->m_fCoord[2] = atof(p);
		#endif

		for (z2=0;z2<3;z2++) {
			at->m_iCoord[z2] = FloatToFixed(at->m_fCoord[z2],signi);
			at->m_fCoord[z2] = FixedToFloat(at->m_iCoord[z2],signi);
			m_faCOM[z2] += at->m_fCoord[z2];
		}
	}

	for (z=0;z<3;z++) {
		m_faCOM[z] /= m_oaAtoms.size();
		m_iaCOM[z] = FloatToFixed(m_faCOM[z],signi);
		m_faCOM[z] = FixedToFloat(m_iaCOM[z],signi);
	}

	for (z=0;z<(int)m_oaAtoms.size();z++) {
		for (z2=0;z2<3;z2++) {
			//m_oaAtoms[z]->m_fRelCoord[z2] = m_oaAtoms[z]->m_fCoord[z2] - m_faCOM[z2];
			//m_oaAtoms[z]->m_iRelCoord[z2] = FloatToFixed(m_oaAtoms[z]->m_fRelCoord[z2],signi);
			m_oaAtoms[z]->m_iRelCoord[z2] = m_oaAtoms[z]->m_iCoord[z2] - m_iaCOM[z2];
			m_oaAtoms[z]->m_fRelCoord[z2] = FixedToFloat(m_oaAtoms[z]->m_iRelCoord[z2],signi);
		}
	}

	if (ref != NULL)
		WriteXYZ(ref,signi);

	return true;
}



void CBQBAtomSet::FloatToInt() {

	unsigned long z, z2;
	CBQBAtom *at;


	for (z=0;z<3;z++)
		m_faCOM[z] = 0;

	for (z=0;z<m_oaAtoms.size();z++) {

		at = m_oaAtoms[z];

		for (z2=0;z2<3;z2++) {
			at->m_iCoord[z2] = FloatToFixed(at->m_fCoord[z2],m_iSigni);
			at->m_fCoord[z2] = FixedToFloat(at->m_iCoord[z2],m_iSigni);
			m_faCOM[z2] += at->m_fCoord[z2];
		}

	}

	for (z=0;z<3;z++) {
		m_faCOM[z] /= m_oaAtoms.size();
		m_iaCOM[z] = FloatToFixed(m_faCOM[z],m_iSigni);
		m_faCOM[z] = FixedToFloat(m_iaCOM[z],m_iSigni);
	}

	for (z=0;z<m_oaAtoms.size();z++) {
		for (z2=0;z2<3;z2++) {
			m_oaAtoms[z]->m_iRelCoord[z2] = m_oaAtoms[z]->m_iCoord[z2] - m_iaCOM[z2];
			m_oaAtoms[z]->m_fRelCoord[z2] = FixedToFloat(m_oaAtoms[z]->m_iRelCoord[z2],m_iSigni);
		}
	}
}



bool CBQBAtomSet::SkipXYZ(FILE *a) {

	char buf[256];
	int i, z;

	(void)!fgets(buf,256,a);
	buf[strlen(buf)-1] = 0;
	i = atoi(buf);

	(void)!fgets(buf,256,a);

	if (feof(a)) {
		m_IF.eprintf("CBQBAtomSet::SkipXYZ(): Error: Unexpected end of file (A).\n");
		m_IF.eprintf("The XYZ file is corrupted.\n");
		return false;
	}

	for (z=0;z<i;z++) {
		(void)!fgets(buf,256,a);
		if (feof(a)) {
			m_IF.eprintf("CBQBAtomSet::SkipXYZ(): Error: Unexpected end of file (Atom %d/%d).\n",z+1,i);
			m_IF.eprintf("The XYZ file is corrupted.\n");
			return false;
		}
	}

	return true;
}

