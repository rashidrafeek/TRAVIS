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



#ifndef BQB_CRC32_H
#define BQB_CRC32_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"
#include <vector>
#include <stdlib.h>



/* This code is adapted from  http://www.rajivchakravorty.com/source-code/uncertainty/multimedia-sim/html/crc8_8c-source.html  */


class CBQBCRC8 {
public:

	CBQBCRC8();


	void Reset() {
		m_iCRC8 = 0xFF;
	}


	void Finalize() {
		m_iCRC8 ^= 0xFF;
	}


	void PushByte(unsigned char m) {
		// For a byte array whose accumulated crc value is stored in *crc, computes
		// resultant crc obtained by appending m to the byte array
		m_iCRC8 = m_iTable[m_iCRC8 ^ m];
		m_iCRC8 &= 0xFF;
	}


	unsigned char m_iCRC8;
	unsigned char m_iTable[256];
};



/* This code is adapted from   http://www.networkdls.com/Software/View/CRC32   */


class CBQBCRC32 {
public:

	CBQBCRC32();

	unsigned int ComputeCRC32(const std::vector<unsigned char> &data, int from=0, int to=-1) const;

	unsigned int ComputeCRC32_Begin(const std::vector<unsigned char> &data, int from=0, int to=-1) const;
	unsigned int ComputeCRC32_Continue(const std::vector<unsigned char> &data, unsigned int last, int from=0, int to=-1) const;
	unsigned int ComputeCRC32_Finish(unsigned int last) const;

private:
	unsigned int Reflect(unsigned int iReflect, const char cChar) const;
	unsigned int iTable[256]; // CRC lookup table array.
};


#endif

