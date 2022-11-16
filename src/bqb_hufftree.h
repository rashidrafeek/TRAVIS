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



#ifndef BQB_HUFFTREE_H
#define BQB_HUFFTREE_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"
#include <vector>
#include "bqb_bitset.h"



class CBQBInterface;



class CBQBHuffmanSymbol {
public:

	CBQBHuffmanSymbol() 
		: m_iSymbol(0), m_pBitString(NULL) { }

	~CBQBHuffmanSymbol() {
		for (int z=0;z<(int)m_oaChildren.size();z++)
			delete m_oaChildren[z];
		if (m_pBitString != NULL)
			delete m_pBitString;
	}

	void REC_BuildBitstrings(int depth);

	int m_iSymbol;
	int m_iDepth;
	int m_iFrequency;
	bool m_bVirtual;
	std::vector<CBQBHuffmanSymbol*> m_oaChildren;
	CBQBBitSet *m_pBitString;
};



class CBQBHuffmanEstimator {

public:
	void Init(int alphasize);
	unsigned int EstimateBitLength(const std::vector<int> &ia);
	unsigned int EstimateBitLength(const std::vector<int> &ia, int i2);
	unsigned int EstimateBitLength(const std::vector<int> &ia, const std::vector<int> &ia2);
	unsigned int EstimateBitLength(const std::vector<int> &ia, const std::vector<int> &ia2, const std::vector<int> &ia3);
	void BuildBitLengthTable(const std::vector<int> &ia);
	std::vector<int> m_iaLength;

private:
	std::vector<unsigned int> m_iaTempFrequencies;
	std::vector<int> m_iaTempHeap;
	std::vector<int> m_iaTempParent;
	int m_iAlphaSize;
};



class CBQBHuffmanTree {
public:

	explicit CBQBHuffmanTree(CBQBInterface &i)
		: m_pTree(NULL), m_iMaxBitLength(0), m_IF(i) { }


	~CBQBHuffmanTree() {
		if (m_pTree != NULL)
			delete m_pTree;
	}


	int DecodeSymbol(CBQBBitSet *bs) const {
		CBQBHuffmanSymbol *hs = m_pTree;
//		printf("********\n");
		while (hs->m_oaChildren.size() == 2)
			if (bs->ReadBit()) {
//				printf("A\n");
				hs = hs->m_oaChildren[1];
			} else {
//				printf("B\n");
				hs = hs->m_oaChildren[0];
			}
		return hs->m_iSymbol;
	}


	void Init(int alphasize);
	void BuildTree(bool canonical=false, bool showsymbols=false);
	void BuildPrelimTree();

	void REC_BuildBitstringLengths(int depth, CBQBHuffmanSymbol *sym);

	int ExportTree(CBQBBitSet *bs, bool chr);
	int REC_ExportTree(CBQBBitSet *bs, CBQBHuffmanSymbol *sym, int bits);

	void ImportTree(CBQBBitSet *bs, bool chr);
	int REC_ImportTree(CBQBBitSet *bs, CBQBHuffmanSymbol *sym, int bits);

	void REC_PushCanonicalSymbol(CBQBHuffmanSymbol *sym, CBQBHuffmanSymbol *ts, unsigned long b, int depth);

	bool m_bCanonical;
	CBQBHuffmanSymbol *m_pTree;
	int m_iMaxBitLength;
	std::vector<int> m_iaFrequencies;
	std::vector<int> m_iaLengths;
	std::vector<CBQBBitSet*> m_oaBitStrings;
	std::vector<CBQBHuffmanSymbol*> m_oaSymbols;
	std::vector<CBQBHuffmanSymbol*> m_oaTempSymbols;

private:
	CBQBInterface &m_IF;
};



#endif


