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



#ifndef BQB_CUBEFRAME_H
#define BQB_CUBEFRAME_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"
#include <vector>
#include <string>
#include <string.h>



class CBQBInterface;



class CBQBAtom {
public:

	CBQBAtom()
		: m_iOrd(0) { }


	int m_iOrd;
	std::string m_sLabel;
	double m_fCoord[3];
	double m_fRelCoord[3];
	long m_iCoord[3];
	long m_iRelCoord[3];
};



class CBQBAtomSet {
public:

	explicit CBQBAtomSet(CBQBInterface &i)
		: m_bOrd(false), m_bLabels(false), m_sComment(NULL), m_IF(i) {
		for (int z=0;z<3;z++) {
			m_faCOM[z] = 0;
			m_iaCOM[z] = 0;
		}
	}

	~CBQBAtomSet() {
		for (int z=0;z<(int)m_oaAtoms.size();z++)
			delete m_oaAtoms[z];
		if (m_sComment != NULL) {
			delete[] m_sComment;
			m_sComment = NULL;
		}
	}


	bool ReadXYZ(FILE *a, int signi, FILE *ref);
	void WriteXYZ(FILE *a, int signi);
	bool SkipXYZ(FILE *a);

	void FloatToInt();

	bool m_bOrd;
	bool m_bLabels;
	int m_iSigni;
	double m_faCOM[3];
	int m_iaCOM[3];
	std::vector<CBQBAtom*> m_oaAtoms;
	char *m_sComment;

private:
	CBQBInterface &m_IF;
};



class CBQBCubeFrame {
public:

	explicit CBQBCubeFrame(CBQBInterface &i)
		: m_iResXY(0), m_iResYZ(0), m_iResXYZ(0), m_pAtoms(NULL), m_iID(-1), m_IF(i) {
		for (int z=0;z<3;z++) {
			m_iRes[z] = 0;
			m_fMinVal[z] = 0;
			m_fMaxVal[z] = 0;
			m_fCenter[z] = 0;
			m_fStrideA[z] = 0;
			m_fStrideB[z] = 0;
			m_fStrideC[z] = 0;
			m_iCenter[z] = 0;
			m_iStrideA[z] = 0;
			m_iStrideB[z] = 0;
			m_iStrideC[z] = 0;
		}
	}


	~CBQBCubeFrame() {
//		if (m_pAtoms != NULL)
//			delete m_pAtoms;
	}


	void CopyHeader(const CBQBCubeFrame *p) {

		m_iRes[0] = p->m_iRes[0];
		m_iRes[1] = p->m_iRes[1];
		m_iRes[2] = p->m_iRes[2];
		m_fStrideA[0] = p->m_fStrideA[0];
		m_fStrideA[1] = p->m_fStrideA[1];
		m_fStrideA[2] = p->m_fStrideA[2];
		m_iStrideA[0] = p->m_iStrideA[0];
		m_iStrideA[1] = p->m_iStrideA[1];
		m_iStrideA[2] = p->m_iStrideA[2];
		m_fStrideB[0] = p->m_fStrideB[0];
		m_fStrideB[1] = p->m_fStrideB[1];
		m_fStrideB[2] = p->m_fStrideB[2];
		m_iStrideB[0] = p->m_iStrideB[0];
		m_iStrideB[1] = p->m_iStrideB[1];
		m_iStrideB[2] = p->m_iStrideB[2];
		m_fStrideC[0] = p->m_fStrideC[0];
		m_fStrideC[1] = p->m_fStrideC[1];
		m_fStrideC[2] = p->m_fStrideC[2];
		m_iStrideC[0] = p->m_iStrideC[0];
		m_iStrideC[1] = p->m_iStrideC[1];
		m_iStrideC[2] = p->m_iStrideC[2];
		m_fCenter[0] = p->m_fCenter[0];
		m_fCenter[1] = p->m_fCenter[1];
		m_fCenter[2] = p->m_fCenter[2];
		m_iCenter[0] = p->m_iCenter[0];
		m_iCenter[1] = p->m_iCenter[1];
		m_iCenter[2] = p->m_iCenter[2];

		m_iResXY = m_iRes[0] * m_iRes[1];
		m_iResYZ = m_iRes[1] * m_iRes[2];
		m_iResXYZ = m_iRes[0] * m_iRes[1] * m_iRes[2];

		m_faBin.resize(m_iResXYZ);
		m_iaExpo.resize(m_iResXYZ);
		m_iaMantis.resize(m_iResXYZ);
	}

	bool ReadFrame(FILE *a, int eps, int csigni, int asigni, bool verbose=false);
	bool SkipFrame(FILE *a, bool verbose=false) const;

	void WriteFrame(FILE *a, bool verbose=false);
	void WriteFrame_Double(FILE *a, bool verbose=false);

	void FloatToInt();

	int m_iEps;
	int m_iSigni;
	std::vector<double> m_faBin;
	std::vector<int> m_iaMantis;
	std::vector<char> m_iaExpo;
	int m_iRes[3];
	int m_iResXY;
	int m_iResYZ;
	int m_iResXYZ;
	double m_fMinVal[3];
	double m_fMaxVal[3];
	double m_fCenter[3];
	double m_fStrideA[3];
	double m_fStrideB[3];
	double m_fStrideC[3];
	long m_iCenter[3];
	long m_iStrideA[3];
	long m_iStrideB[3];
	long m_iStrideC[3];
	CBQBAtomSet *m_pAtoms;
	int m_iID;

private:
	CBQBInterface &m_IF;
};



#endif

