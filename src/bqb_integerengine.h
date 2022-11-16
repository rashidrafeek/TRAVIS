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



#ifndef BQB_INTEGERENGINE_H
#define BQB_INTEGERENGINE_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"
#include <vector>
#include "bqb_bitset.h"
#include "bqb_hufftree.h"
#include "bqb_alphabet.h"
#include "bqb_parmset.h"



class CBQBInterface;



class CBQBBWPair {
public:
	CBQBBWPair(int symbol, int index) : m_iSymbol(symbol), m_iIndex(index) { }
	int m_iSymbol;
	int m_iIndex;
};



class CBQBIntegerEngine {
public:

	explicit CBQBIntegerEngine(CBQBInterface &i) : m_IF(i) {
	}


	bool Compress(
		std::vector<int> &inp,
		CBQBBitSet *outp,
		bool bw,
		bool mtf,
		bool coderun,
		int blocklength,
		int tables,
		bool opttables,
		bool chr,
		int maxiter,
		int maxchunk,
		CBQBStatistics *stat
	);


	bool Compress(
		std::vector<int> &inp,
		CBQBBitSet *outp,
		const CBQBParameterSet_Compressor *parm,
		bool chr,
		CBQBStatistics *stat
	);


	bool Decompress(
		CBQBBitSet *inp,
		std::vector<int> &outp,
		CBQBParameterSet_Compressor *parm
	);


	bool CompressSingle(
		std::vector<int> &inp,
		CBQBBitSet *outp,
		bool bw,
		bool mtf,
		bool coderun,
		int blocklength,
		int tables,
		bool opttables,
		bool chr,
		int maxiter,
		CBQBStatistics *stat
	);


	bool CompressSingle(
		std::vector<int> &inp,
		CBQBBitSet *outp,
		const CBQBParameterSet_Compressor *parm,
		bool chr,
		CBQBStatistics *stat
	);


	bool DecompressSingle(
		CBQBBitSet *inp,
		std::vector<int> &outp,
		CBQBParameterSet_Compressor *parm
	);


	std::vector<int> *m_iaBW;

	std::vector<int> m_iaBWRunTable;


	// These structs are a hack to enable the comparison function for std::sort
	//   to access class members of CBQBIntegerEngine.
	// Can't use global variables here, because multiple CBQBIntegerEngines
	//   might be running in different threads...

	struct SSortBW {

		explicit SSortBW( const CBQBIntegerEngine& ie ) : m_IE(ie) {
		}

		const CBQBIntegerEngine &m_IE;

		bool operator()( const int & i1, const int & i2  ) {

			int k1, k2, i;

			if (i1 == i2)
				return false;

			k1 = i1;
			k2 = i2;

			i = 0;
			while ((*m_IE.m_iaBW)[k1] == (*m_IE.m_iaBW)[k2]) {
				k1++;
				k2++;
				i++;
				if (k1 >= (int)m_IE.m_iaBW->size())
					k1 = 0;
				if (k2 >= (int)m_IE.m_iaBW->size())
					k2 = 0;
				if (i >= (int)m_IE.m_iaBW->size())
					return false;
			}

			return ((*m_IE.m_iaBW)[k1] < (*m_IE.m_iaBW)[k2]);
		}
	};


	struct SSortBWRuntable {

		explicit SSortBWRuntable( const CBQBIntegerEngine& ie ) : m_IE(ie) {
		}

		const CBQBIntegerEngine &m_IE;

		bool operator()( const int & i1, const int & i2  ) {

			int k1, k2, i, j;

			if (i1 == i2)
				return false;

			k1 = i1;
			k2 = i2;

			i = 0;
			while ((*m_IE.m_iaBW)[k1] == (*m_IE.m_iaBW)[k2]) {
				j = MIN( m_IE.m_iaBWRunTable[k1] , m_IE.m_iaBWRunTable[k2] );
				i += j;
				k1 += j;
				k2 += j;
				if (k1 >= (int)m_IE.m_iaBW->size())
					k1 -= (int)m_IE.m_iaBW->size();
				if (k2 >= (int)m_IE.m_iaBW->size())
					k2 -= (int)m_IE.m_iaBW->size();
				if (i >= (int)m_IE.m_iaBW->size())
					return false;
			}

			return ((*m_IE.m_iaBW)[k1] < (*m_IE.m_iaBW)[k2]);
		}
	};


private:

	CBQBInterface &m_IF;


	//void MultiHuffmanOptimize(int tables, std::vector<CBQBHuffmanTree*> &hta, CBQBAlphabet *alp, std::vector<int> &tia, std::vector<int> &asi, std::vector<int> &iasi, int blocklength, bool verbose);

};


#endif


