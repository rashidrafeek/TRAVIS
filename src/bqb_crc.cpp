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



// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_crc.h"


const char *GetRevisionInfo_bqb_crc(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_crc() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



/* This code is adapted from  http://www.rajivchakravorty.com/source-code/uncertainty/multimedia-sim/html/crc8_8c-source.html  */


CBQBCRC8::CBQBCRC8() {

	#define CRC8_DI  0x2F    /* x^8 + x^5 + x^3 + x^2 + x + 1 */

	int i, j;
	unsigned char crc;

	for (i=0; i<256; i++) {
		crc = (unsigned char)i;
		for (j=0; j<8; j++)
			crc = (crc << 1) ^ ((crc & 0x80) ? CRC8_DI : 0);
		m_iTable[i] = crc & 0xFF;
		/* printf("table[%d] = %d (0x%X)\n", i, crc, crc); */
	}

	Reset();
}




/* This code is adapted from   http://www.networkdls.com/Software/View/CRC32   */


CBQBCRC32::CBQBCRC32() {

	//0x04C11DB7 is the official polynomial used by PKZip, WinZip and Ethernet.
	unsigned int iPolynomial = 0x04C11DB7;

	memset(&this->iTable, 0, sizeof(this->iTable));

	// 256 values representing ASCII character codes.
	for(int iCodes = 0; iCodes <= 0xFF; iCodes++) {
		this->iTable[iCodes] = this->Reflect(iCodes, 8) << 24;

		for(int iPos = 0; iPos < 8; iPos++) {
			this->iTable[iCodes] = (this->iTable[iCodes] << 1)
				^ ((this->iTable[iCodes] & (1 << 31)) ? iPolynomial : 0);
		}

		this->iTable[iCodes] = this->Reflect(this->iTable[iCodes], 32);
	}
}


unsigned int CBQBCRC32::Reflect(unsigned int iReflect, const char cChar) const {

	unsigned int iValue = 0;

	// Swap bit 0 for bit 7, bit 1 for bit 6, etc....
	for (int iPos = 1; iPos < (cChar + 1); iPos++) {
		if (iReflect & 1)
			iValue |= (1 << (cChar - iPos));
		iReflect >>= 1;
	}

	return iValue;
}


unsigned int CBQBCRC32::ComputeCRC32(const std::vector<unsigned char> &data, int from, int to) const {

	unsigned int crc;
	int z;

	crc = 0xffffffff; // Initialize the CRC.

	if (to == -1)
		to = (int)data.size()-1;

	for (z=from;z<=to;z++)
		crc = (crc >> 8) ^ this->iTable[(crc & 0xFF) ^ data[z]];

	crc ^= 0xffffffff; // Finalize the CRC.

	return crc;
}


unsigned int CBQBCRC32::ComputeCRC32_Begin(const std::vector<unsigned char> &data, int from, int to) const {

	unsigned int crc;
	int z;

	crc = 0xffffffff; // Initialize the CRC.

	if (to == -1)
		to = (int)data.size()-1;

	for (z=from;z<=to;z++)
		crc = (crc >> 8) ^ this->iTable[(crc & 0xFF) ^ data[z]];

	return crc;
}


unsigned int CBQBCRC32::ComputeCRC32_Continue(const std::vector<unsigned char> &data, unsigned int last, int from, int to) const {

	unsigned int crc;
	int z;

	crc = last; // Initialize the CRC.

	if (to == -1)
		to = (int)data.size()-1;

	for (z=from;z<=to;z++)
		crc = (crc >> 8) ^ this->iTable[(crc & 0xFF) ^ data[z]];

	return crc;
}


unsigned int CBQBCRC32::ComputeCRC32_Finish(unsigned int last) const {

	return last ^ 0xffffffff; // Finalize the CRC.
}


