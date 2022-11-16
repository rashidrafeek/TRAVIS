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



#ifndef BQB_EXTRAPOLATOR_H
#define BQB_EXTRAPOLATOR_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_cubeframe.h"
#include "bqb_bitset.h"
#include <vector>
#include "bqb_parmset.h"



class CBQBInterface;


//#define CHECK_EXTRAPOLATOR



void TestExtrapolator();



class CBQBExtraPattern {
public:
	std::vector<int> m_iaSpatialIndex;
	std::vector<int> m_iaTemporalIndex;
	std::vector<double> m_faCoeff;
};




class CBQBExtrapolator {
public:

	explicit CBQBExtrapolator(CBQBInterface &i) : m_IF(i) {
	}


	~CBQBExtrapolator() {
	}


	void Reset();


	void Initialize(
		int resx,
		int resy,
		int resz,
		int srangex,
		int srangey,
		int srangez,
		int trange,
		int sorder,
		int torder,
		int offsetx,
		int offsety,
		int offsetz,
		bool crosss,
		bool crosst,
		bool wrap,
		bool crossranges,
		bool crossranget,
		double distexpo,
		double timeexpo,
		bool silent
	);


	void PushCubeFrame(const CBQBCubeFrame *frame);

	void PushAtomFrame(const CBQBAtomSet *frame);


	void InitializeXYZ(int trange, int torder, double timeexpo, bool silent) {

		Initialize(
			0,
			0,
			0,
			1,
			1,
			1,
			trange,
			0,
			torder,
			0,
			0,
			0,
			false,
			false,
			false,
			false,
			false,
			0.0,
			timeexpo,
			//verbose,
			silent
		);
	}


	void InitializeXYZ(CBQBParameterSet_Position *p, bool silent);


	const CBQBCubeFrame* GetCubeFrame(int depth) const;


	const CBQBAtomSet* GetAtomFrame(int depth) const;


	double Extrapolate(int x, int y, int z) const {

		return ExtrapolateKnown( FindPattern(m_iTRange-m_iCubeFrameVisibleCount,x,y,z), x*m_iRes12 + y*m_iRes[2] + z );
	}


	double Extrapolate(int x, int y, int z, int index) const {

		return ExtrapolateKnown( FindPattern(m_iTRange-m_iCubeFrameVisibleCount,x,y,z), index );
	}


	double ExtrapolatePred(int index) const {

		return ExtrapolateKnownPred( m_iTRange-m_iCubeFrameVisibleCount, index );
	}


	double ExtrapolateCorr(const std::vector<double> &fa, int x, int y, int z, int index) const {

		return ExtrapolateKnownCorr( fa, FindPatternCorr(x,y,z), index );
	}


	double ExtrapolateXYZ(int index, int index2) const {

		return ExtrapolateKnownXYZ( m_iTRange-m_iCubeFrameVisibleCount, index, index2 );
	}


	double ExtrapolateKnown(int pattern, int index) const;


	double ExtrapolateKnownCorr(const std::vector<double> &fa, int pattern, int index) const;


	double ExtrapolateKnownPred(int pattern, int index) const {

		double tf;
		unsigned int z;

		const CBQBExtraPattern &p = m_oaPatterns[pattern];
		tf = 0;
		for (z=0;z<p.m_iaSpatialIndex.size();z++)
			tf += p.m_faCoeff[z] * GetCubeFrame(p.m_iaTemporalIndex[z])->m_faBin[index];
		return tf;
	}


	double ExtrapolateKnownXYZ(int pattern, int index, int index2) const {

		double tf;
		unsigned int z;

		const CBQBExtraPattern &p = m_oaPatterns[pattern];
		tf = 0;
		for (z=0;z<p.m_iaSpatialIndex.size();z++) {
			const CBQBAtomSet *afr = GetAtomFrame(p.m_iaTemporalIndex[z]);
			if (index < (int)afr->m_oaAtoms.size())
				tf += p.m_faCoeff[z] * afr->m_oaAtoms[index]->m_fRelCoord[index2];
			else
				tf += p.m_faCoeff[z] * afr->m_faCOM[index2];
		}
		return tf;
	}


	int FindPattern(int t, int x, int y, int z) const {

		int ix, iy, iz;

		if (x < m_iSOffset[0])
			ix = m_iSOffset[0] - x;
		else if (x > m_iTempVal[0]/*m_iRes[0]-m_iSRange[0]+m_iSOffset[0]*/)
			ix = m_iSOffset[0] + x - m_iTempVal[0]/*m_iRes[0]+m_iSRange[0]-m_iSOffset[0]*/;
		else
			ix = 0;

		if (y < m_iSOffset[1])
			iy = m_iSOffset[1] - y;
		else if (y > m_iTempVal[1]/*m_iRes[1]-m_iSRange[1]+m_iSOffset[1]*/)
			iy = m_iSOffset[1] + y - m_iTempVal[1]/*m_iRes[1]+m_iSRange[1]-m_iSOffset[1]*/;
		else
			iy = 0;

		if (z < m_iSOffset[2])
			iz = m_iSOffset[2] - z;
		else if (z > m_iTempVal[2]/*m_iRes[2]-m_iSRange[2]+m_iSOffset[2]*/)
			iz = m_iSOffset[2] + z - m_iTempVal[2]/*m_iRes[2]+m_iSRange[2]-m_iSOffset[2]*/;
		else
			iz = 0;

		return t * m_iSRange012/*m_iSRange[0]*m_iSRange[1]*m_iSRange[2]*/
			+ iz * m_iSRange01/*m_iSRange[0]*m_iSRange[1]*/ + iy * m_iSRange[0] + ix;
	}


	int FindPatternCorr(int x, int y, int z) const {

		int ix, iy, iz;

		if (x < m_iSOffset[0])
			ix = m_iSOffset[0] - x;
		else if (x > m_iTempVal[0]/*m_iRes[0]-m_iSRange[0]+m_iSOffset[0]*/)
			ix = m_iSOffset[0] + x - m_iTempVal[0]/*m_iRes[0]+m_iSRange[0]-m_iSOffset[0]*/;
		else
			ix = 0;

		if (y < m_iSOffset[1])
			iy = m_iSOffset[1] - y;
		else if (y > m_iTempVal[1]/*m_iRes[1]-m_iSRange[1]+m_iSOffset[1]*/)
			iy = m_iSOffset[1] + y - m_iTempVal[1]/*m_iRes[1]+m_iSRange[1]-m_iSOffset[1]*/;
		else
			iy = 0;

		if (z < m_iSOffset[2])
			iz = m_iSOffset[2] - z;
		else if (z > m_iTempVal[2]/*m_iRes[2]-m_iSRange[2]+m_iSOffset[2]*/)
			iz = m_iSOffset[2] + z - m_iTempVal[2]/*m_iRes[2]+m_iSRange[2]-m_iSOffset[2]*/;
		else
			iz = 0;

		return iz * m_iSRange01/*m_iSRange[0]*m_iSRange[1]*/ + iy * m_iSRange[0] + ix;
	}



	int m_iRes[3];
	int m_iSOffset[3];
	int m_iTempVal[3];


	int m_iSRange[3];
	int m_iSOrder;
	int m_iTRange;
	int m_iTOrder;
	bool m_bCrossS;
	bool m_bCrossRangeS;
	bool m_bWrap;
	bool m_bCrossT;
	bool m_bCrossRangeT;
	double m_fDistExpo;
	double m_fTimeExpo;

	double m_fDistExpoAdd;
	double m_fTimeExpoAdd;

	int m_iCubeFramePos;
	int m_iCubeFrameCount;
	int m_iCubeFrameVisibleCount;


private:

	int m_iSRange01;
	int m_iSRange012;
	int m_iRes012;
	int m_iRes12;

	std::vector<const CBQBCubeFrame*> m_oaCubeFrames;
	std::vector<const CBQBAtomSet*> m_oaAtomFrames;
	std::vector<CBQBExtraPattern> m_oaPatterns;

	CBQBInterface &m_IF;
};



#endif



