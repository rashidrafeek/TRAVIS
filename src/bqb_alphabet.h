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



#ifndef BQB_ALPHABET_H
#define BQB_ALPHABET_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"
#include <vector>
#include "bqb_bitset.h"



class CBQBInterface;



class CBQBAlphabetEntry {
public:

	CBQBAlphabetEntry()
		: m_iSymbol(0), m_iFrequency(0), m_iIndex(0) { }

	int m_iSymbol;
	int m_iFrequency;
	int m_iIndex;
};



class CBQBAlphabet {
public:

	explicit CBQBAlphabet(CBQBInterface &i) : m_IF(i) {
	}

	~CBQBAlphabet() {
		for (int z=0;z<(int)m_oaAlphabet.size();z++)
			delete m_oaAlphabet[z];
	}

	void BuildAlphabet(const std::vector<int> &inp, bool reserved=true);

	void Export(CBQBBitSet *bs, bool huffman, bool chr) const;

	void Import(CBQBBitSet *bs, bool chr);

	int FindIndex(int symbol) const;

	void RecalcFrequencies(const std::vector<int> &ia);

	std::vector<CBQBAlphabetEntry*> m_oaAlphabet;
	std::vector<CBQBAlphabetEntry*> m_oaAlphabetSortFreq;
	std::vector<int> m_iaIndices;

private:
	CBQBInterface &m_IF;
};



#endif

