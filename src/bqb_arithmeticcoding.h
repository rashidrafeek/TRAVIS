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



#ifndef BQB_ARITHMETICCODING_H
#define BQB_ARITHMETICCODING_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_bitset.h"
#include <vector>
#include "bqb_arithmeticcoding_model.h"



#ifdef BQI_DEVELOPMENT



class CBQBInterface;



inline void put_bit_plus_pending( unsigned char bit, unsigned int &pending_bits, CBQBBitSet *bs ) {

	unsigned int i;

	bs->WriteBit( bit );
	for ( i=0; i<pending_bits; i++ )
		bs->WriteBit( (bit==0)?1:0 ); // Flipped bit
	pending_bits = 0;
}



class CBQBACEngine {
public:

	explicit CBQBACEngine(CBQBInterface &i) : m_pModel(NULL), m_iSymbolRange(0), m_IF(i) {
	}



	~CBQBACEngine() {
	}



	void InitCompress() {

		m_iPendingBits = 0;
		m_iLow = 0;
		m_iHigh = BQB_AC_MAX_CODE;
	}



	void Compress_PushSymbol( CBQBBitSet *bs, BQB_AC_CODE_VALUE plow, BQB_AC_CODE_VALUE phigh, BQB_AC_CODE_VALUE pcount ) {

		BQB_AC_CODE_VALUE range;


		if (pcount == 0) {
			printf("CBQBACEngine::Compress(): Error: Model returned zero count.\n");
			abort();
		}

		if (plow == phigh) {
			printf("CBQBACEngine::Compress(): Error: Model returned low == high.\n");
			abort();
		}

		range = m_iHigh - m_iLow + 1;
		m_iHigh = m_iLow + (range * phigh / pcount) - 1;
		m_iLow = m_iLow + (range * plow / pcount);

		while (true) {

			if (m_iHigh < BQB_AC_ONE_HALF)
				put_bit_plus_pending( 0, m_iPendingBits, bs );
			else if (m_iLow >= BQB_AC_ONE_HALF)
				put_bit_plus_pending( 1, m_iPendingBits, bs );
			else if ((m_iLow >= BQB_AC_ONE_FOURTH) && (m_iHigh < BQB_AC_THREE_FOURTHS)) {
				m_iPendingBits++;
				m_iLow -= BQB_AC_ONE_FOURTH;  
				m_iHigh -= BQB_AC_ONE_FOURTH;  
			} else
				break;

			m_iHigh <<= 1;
			m_iHigh++;
			m_iLow <<= 1;
			m_iHigh &= BQB_AC_MAX_CODE;
			m_iLow &= BQB_AC_MAX_CODE;
		}
	}



	void FinishCompress( CBQBBitSet *bs ) {

		m_iPendingBits++;

		if (m_iLow < BQB_AC_ONE_FOURTH)
			put_bit_plus_pending( 0, m_iPendingBits, bs );
		else
			put_bit_plus_pending( 1, m_iPendingBits, bs );
	}



	bool Compress(
		std::vector<int> *iain,
		CBQBBitSet *bsout,
		unsigned int symrange
	);



	bool Decompress(
		CBQBBitSet *bsin,
		std::vector<int> *iaout
	);



	CBQBACModelBase *m_pModel;
	unsigned int m_iSymbolRange;

	unsigned int m_iPendingBits;
	BQB_AC_CODE_VALUE m_iHigh;
	BQB_AC_CODE_VALUE m_iLow;

private:
	CBQBInterface &m_IF;
};



#endif



#endif


