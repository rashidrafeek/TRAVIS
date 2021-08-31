/***********************************************************************************

    LibBQB - File Format and Compression Algorithms for Trajectories of
             Volumetric Data and Atom Positions

    https://brehm-research.de/bqb

    Free software, licensed under GNU LGPL v3

    Copyright (c) Martin Brehm and Martin Thomas,
                  Martin Luther University Halle-Wittenberg, Germany,
                  2016 - 2021.

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



#ifndef BQB_ARITHMETICCODING_MODEL_H
#define BQB_ARITHMETICCODING_MODEL_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_bitset.h"
#include <vector>



#ifdef BQI_DEVELOPMENT



#ifdef TARGET_WINDOWS

	#define BQB_AC_CODE_VALUE  unsigned int

	#define BQB_AC_CODE_VALUE_BITS    (17)  // (std::numeric_limits<CODE_VALUE_>::digits + 3) / 2

	#define BQB_AC_MAX_CODE       (131071)  // = (CODE_VALUE(1) << CODE_VALUE_BITS) - 1;
	#define BQB_AC_MAX_FREQ        (32767)  // = (CODE_VALUE(1) << FREQUENCY_BITS) - 1;
	#define BQB_AC_ONE_FOURTH      (32768)  // = CODE_VALUE(1) << (CODE_VALUE_BITS - 2);
	#define BQB_AC_ONE_HALF        (65536)  // = 2 * ONE_FOURTH;
	#define BQB_AC_THREE_FOURTHS   (98304)  // = 3 * ONE_FOURTH;

#endif



#ifdef TARGET_LINUX

	#define BQB_AC_CODE_VALUE  unsigned long

	#define BQB_AC_CODE_VALUE_BITS          (34)  // (std::numeric_limits<CODE_VALUE_>::digits + 3) / 2

	#define BQB_AC_MAX_CODE        (17179869183)  // = (CODE_VALUE(1) << CODE_VALUE_BITS) - 1;
	#define BQB_AC_MAX_FREQ         (1073741823)  // = (CODE_VALUE(1) << FREQUENCY_BITS) - 1;
	#define BQB_AC_ONE_FOURTH       (4294967296)  // = CODE_VALUE(1) << (CODE_VALUE_BITS - 2);
	#define BQB_AC_ONE_HALF         (8589934592)  // = 2 * ONE_FOURTH;
	#define BQB_AC_THREE_FOURTHS   (12884901888)  // = 3 * ONE_FOURTH;

#endif



//#define BQB_AC_FREQUENCY_BITS     (15)  // std::numeric_limits<CODE_VALUE_>::digits - CODE_VALUE_BITS_>
//#define BQB_AC_PRECISION          (32)  // std::numeric_limits<CODE_VALUE>::digits


//#define BQB_AC_MAX_FREQ         (4095)  // = (CODE_VALUE(1) << FREQUENCY_BITS) - 1;
//#define BQB_AC_MAX_FREQ         (8191)  // = (CODE_VALUE(1) << FREQUENCY_BITS) - 1;
//#define BQB_AC_MAX_FREQ        (16383)  // = (CODE_VALUE(1) << FREQUENCY_BITS) - 1;
//#define BQB_AC_MAX_CODE       (1048575)  // = (CODE_VALUE(1) << CODE_VALUE_BITS) - 1;
//#define BQB_AC_ONE_FOURTH      (262144)  // = CODE_VALUE(1) << (CODE_VALUE_BITS - 2);
//#define BQB_AC_ONE_HALF        (524288)  // = 2 * ONE_FOURTH;
//#define BQB_AC_THREE_FOURTHS   (786432)  // = 3 * ONE_FOURTH;



class CBQBACProbability {
public:

	CBQBACProbability() : m_iLow(0), m_iHigh(0), m_iCount(0) {
	}


	CBQBACProbability( BQB_AC_CODE_VALUE low, BQB_AC_CODE_VALUE high, BQB_AC_CODE_VALUE count ) :
		m_iLow(low), m_iHigh(high), m_iCount(count) {
	}


	BQB_AC_CODE_VALUE m_iLow;
	BQB_AC_CODE_VALUE m_iHigh;
	BQB_AC_CODE_VALUE m_iCount;
};



class CBQBIntRingBuffer {
public:

	CBQBIntRingBuffer() : m_iSymbolRange(0), m_iPos(-1), m_iDepth(-1) {
	}


	void Clear() {
		m_iPos = -1;
		m_iaEntries.clear();
		m_iaEntries.reserve(m_iDepth);
	}


	void AddInterval(int i) {
		m_iaIntervals.push_back(i);
	}


	void Initialize(int depth, int symrange) {
		unsigned int z, z2;
		m_iDepth = depth;
		m_iSymbolRange = symrange;
		m_iaEntries.reserve(m_iDepth);
		m_iaaCount.resize(m_iaIntervals.size());
		for (z=0;z<m_iaIntervals.size();z++) {
			if (m_iaIntervals[z] >= m_iDepth) {
				printf("CBQBIntRingBuffer::Initialize(): Error: Depth needs to be larger than any interval size.\n");
				abort();
			}
			m_iaaCount[z].resize(m_iSymbolRange);
			for (z2=0;(int)z2<m_iSymbolRange;z2++)
				m_iaaCount[z][z2] = 0;
		}
	}


	void Push(int p) {
		unsigned int z;
		m_iPos++;
		if (m_iPos >= m_iDepth)
			m_iPos = 0;
		if ((int)m_iaEntries.size() < m_iDepth)
			m_iaEntries.push_back(p);
		else
			m_iaEntries[m_iPos] = p;
		for (z=0;z<m_iaIntervals.size();z++) {
			m_iaaCount[z][p]++;
			if (GetEntries() >= m_iaIntervals[z]+1)
				m_iaaCount[z][Get(m_iaIntervals[z])]--;
		}
	}


	int GetEntries() const {
		if ((int)m_iaEntries.size() < m_iDepth)
			return m_iaEntries.size();
		else
			return m_iDepth;
	}


	int Get(int i) const {
		int t;
		t = m_iPos - i;
		if (t < 0)
			t += m_iDepth;
		return m_iaEntries[t];
	}


	int m_iSymbolRange;
	int m_iPos;
	int m_iDepth;
	std::vector<int> m_iaEntries;
	std::vector<int> m_iaIntervals;
	std::vector<std::vector<int> > m_iaaCount;
};



class CBQBACModelBase {
public:

	explicit CBQBACModelBase( unsigned int symrange ) {
		m_iSymbolRange = symrange;
	}


	virtual ~CBQBACModelBase() {
	}


	virtual CBQBACProbability GetProbability( unsigned int sym ) const {
		return CBQBACProbability(
			m_iaCumulativeCount[sym],
			m_iaCumulativeCount[sym+1],
			m_iSum
		);
	}


	BQB_AC_CODE_VALUE GetSum() const {
		return m_iSum;
	}


	virtual CBQBACProbability GetSymbol( BQB_AC_CODE_VALUE scaled_value, unsigned int &sym ) = 0;

	virtual void Update( unsigned int sym ) = 0;



	unsigned int m_iSymbolRange;

	std::vector<BQB_AC_CODE_VALUE> m_iaCumulativeCount;
	BQB_AC_CODE_VALUE m_iSum;
};



class CBQBACModel1 : public CBQBACModelBase {
public:

	explicit CBQBACModel1( unsigned int sc ) : CBQBACModelBase( sc ) {
	}


	void InitUniform();

	void InitFromSample( std::vector<int> &ia );

	CBQBACProbability GetSymbol( BQB_AC_CODE_VALUE scaled_value, unsigned int &sym );

	void Update( unsigned int sym );



	std::vector<BQB_AC_CODE_VALUE> m_iaCount;
};



class CBQBACModel2 : public CBQBACModelBase {
public:

	explicit CBQBACModel2( unsigned int sc ) : CBQBACModelBase( sc ) {
	}


	void InitUniform();

	void InitFromSample( std::vector<int> &ia );

	CBQBACProbability GetSymbol( BQB_AC_CODE_VALUE scaled_value, unsigned int &sym );

	void Update( unsigned int sym );



	std::vector<BQB_AC_CODE_VALUE> m_iaCount;
	std::vector<BQB_AC_CODE_VALUE> m_iaCount2;

	CBQBIntRingBuffer *m_pRingBuffer;
};



class CBQBACModel3 : public CBQBACModelBase {
public:

	explicit CBQBACModel3( unsigned int sc ) : CBQBACModelBase( sc ) {
	}


	void InitUniform();

	void InitFromSample( std::vector<int> &ia );

	CBQBACProbability GetSymbol( BQB_AC_CODE_VALUE scaled_value, unsigned int &sym );

	void Update( unsigned int sym );



	std::vector<std::vector<BQB_AC_CODE_VALUE> > m_iaaCount;
	std::vector<BQB_AC_CODE_VALUE> m_iaSum;

	BQB_AC_CODE_VALUE m_iLastSym;
};



class CBQBACModelB : public CBQBACModelBase {
public:

	explicit CBQBACModelB( unsigned int sc ) : CBQBACModelBase( sc ) {
	}


	void InitUniform();

	void InitFromSample( std::vector<int> &ia );

	CBQBACProbability GetProbability( unsigned int sym ) const;

	CBQBACProbability GetSymbol( BQB_AC_CODE_VALUE scaled_value, unsigned int &sym );

	void Update( unsigned int sym );

	void ChangeState( int sym );


	unsigned int m_iZeroCounter;
	unsigned int m_iZeroSum;

	unsigned int m_iSignCounter;
	unsigned int m_iSignSum;

	std::vector<unsigned int> m_iaExpoCounter;
	std::vector<unsigned int> m_iaExpoSum;

	std::vector<unsigned int> m_iaMantisCounter;
	std::vector<unsigned int> m_iaMantisSum;

	int m_iState;
	int m_iSubState;
};



#endif



#endif


