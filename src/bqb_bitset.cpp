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

#include "bqb_bitset.h"
#include "bqb_crc.h"
#include "bqb_interface.h"


const char *GetRevisionInfo_bqb_bitset(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_bitset() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



const char *g_sBQBKeyAlpha = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+-";


std::string CBQBBitSet::ExportKey() {

	std::string key;
	CBQBCRC8 crc;
	unsigned char c1, c2, c3;
	int z, len;


	len = GetLength();
	//printf("Export: len=%d\n",len);
	if (len > 1023) {
		m_IF.eprintf("CBQBBitSet::ExportKey(): Error: Maximum key length is 1023 bits (Bitset is %d bits).\n",len);
		abort();
	}

	c1 = len%64;
	c2 = (len/64)%16;

	for (z=0;z<GetByteLength();z++)
		crc.PushByte(m_iaData[z]);

	crc.Finalize();

	c2 += (crc.m_iCRC8 % 4) * 16;
	c3 = (crc.m_iCRC8 / 4);

	key = "@";
	key += std::string(1,g_sBQBKeyAlpha[c1]);
	key += std::string(1,g_sBQBKeyAlpha[c2]);
	key += std::string(1,g_sBQBKeyAlpha[c3]);

	ResetReadPos();

	for (z=0;z<len/6;z++) {
		c1 = (unsigned char)ReadBitsInteger(6);
		key += std::string(1,g_sBQBKeyAlpha[c1]);
	}
	if ((len % 6) != 0) {
		c1 = (unsigned char)ReadBitsInteger(len % 6);
		key += std::string(1,g_sBQBKeyAlpha[c1]);
	}

	return key;
}


bool CBQBBitSet::ImportKey(std::string s) {

	CBQBCRC8 crc;
	int z, len;
	const char *p;
	unsigned char c;
	std::vector<unsigned char> tia;
	

	Clear();

	if (s.length() < 4) {
		m_IF.eprintf("CBQBBitSet::ImportKey(): Error: Key requires at least four characters (has only %lu).\n",(unsigned long)s.length());
		return false;
	}
	if (s[0] != '@') {
		m_IF.eprintf("CBQBBitSet::ImportKey(): Error: Key required to start with \"@\" (found \"%c\").\n",s[0]);
		return false;
	}

	tia.resize(s.length()-1);
	for (z=1;z<(int)s.length();z++) {
		p = strchr(g_sBQBKeyAlpha,s[z]);
		if (p == NULL) {
			m_IF.eprintf("CBQBBitSet::ImportKey(): Error: Found invalid character \"%c\" in key.\n",s[z]);
			return false;
		}
		tia[z-1] = (unsigned char)(p - g_sBQBKeyAlpha);
	}

	len = tia[0] + (tia[1]%16)*64;

	//printf("Import: len=%d\n",len);

	if ((len+5)/6 < (int)s.length()-4)
		m_IF.printf("CBQBBitSet::ImportKey(): Warning: Key contains %d extra characters at the end.\n",(int)s.length()-4-(len+5)/6);

	if ((len+5)/6 > (int)s.length()-4) {
		m_IF.eprintf("CBQBBitSet::ImportKey(): Error: Key is too short by %d characters.\n",(len+5)/6-(int)s.length()+4);
		return false;
	}

	c = (tia[1]/16) + tia[2]*4;

	for (z=0;z<len/6;z++)
		WriteBits(tia[z+3],6);

	if ((len % 6) != 0) {
		if (tia.back() >= pow2i(len%6)) {
			m_IF.eprintf("CBQBBitSet::ImportKey(): Error: Last byte too large (%d, %d bits).\n",tia.back(),len%6);
			return false;
		}
		WriteBits(tia.back(),len % 6);
	}

	for (z=0;z<GetByteLength();z++)
		crc.PushByte(m_iaData[z]);

	crc.Finalize();

	if (c != crc.m_iCRC8) {
		m_IF.eprintf("CBQBBitSet::ImportKey(): Error: CRC check failed (%02X vs. %02X).\n",c,crc.m_iCRC8);
		return false;
	}

	return true;
}


void CBQBBitSet::WriteBits(unsigned long i, int bits) {
	if ((int)mylog2(i+1) > bits) {
		m_IF.eprintf("CBitSet::WriteBits(): Error: %lu does not fit into %d bits.\n",i,bits);
		abort();
	}
	for (int z=0;z<bits;z++)
		WriteBit((i&(((unsigned long)1)<<z))!=0?1:0);
}


void CBQBBitSet::WriteSignedBits(long i, int bits) {
	unsigned long ul;
	ul = bqbabs(i);
	if ((int)mylog2(ul+1) > bits-1) {
		m_IF.eprintf("CBitSet::WriteSignedBits(): Error: %ld does not fit into %d bits.\n",i,bits);
		abort();
	}
	for (int z=0;z<bits-1;z++)
		WriteBit((ul&(((unsigned long)1)<<z))!=0?1:0);
	if (i >= 0)
		WriteBit(0);
	else
		WriteBit(1);
}


bool CBQBBitSet::ReadBit() {
	bool b;
	if (m_iReadPosBytes >= (int)m_iaData.size()) {
		m_IF.eprintf("CBQBBitSet::ReadBit(): Error: End of BitSet reached.\n");
		abort();
	}
	b = (m_iaData[m_iReadPosBytes] & (1<<m_iReadPosExcess)) != 0;
	m_iReadPosExcess++;
	if (m_iReadPosExcess == 8) {
		m_iReadPosExcess = 0;
		m_iReadPosBytes++;
	}
	return b;
}

	
bool CBQBBitSet::ReadBitOrZero() {
	bool b;
	if (m_iReadPosBytes >= (int)m_iaData.size())
		return false;
	b = (m_iaData[m_iReadPosBytes] & (1<<m_iReadPosExcess)) != 0;
	m_iReadPosExcess++;
	if (m_iReadPosExcess == 8) {
		m_iReadPosExcess = 0;
		m_iReadPosBytes++;
	}
	return b;
}

	
void CBQBBitSet::Dump() const {
	
	int z, z2, ti;

	m_IF.printf("\n");
	m_IF.printf("BitSet Dump:\n");
	if (m_iExcessBits == 0)
		ti = (int)m_iaData.size();
	else
		ti = (int)m_iaData.size()-1;
	for (z=0;z<ti;z++) {
		m_IF.printf("%6d: ",z*8);
		for (z2=0;z2<8;z2++)
			m_IF.printf("%d ",(m_iaData[z]&(1<<(7-z2)))!=0?1:0);
		m_IF.printf("   (%3d)\n",m_iaData[z]);
	}
	if (m_iExcessBits != 0) {
		m_IF.printf("%6d: ",z*8);
		for (z2=0;z2<8-m_iExcessBits;z2++)
			m_IF.printf("  ");
		for (z2=0;z2<m_iExcessBits;z2++)
			m_IF.printf("%d ",(m_iaData[z]&(1<<(m_iExcessBits-z2-1)))!=0?1:0);
		m_IF.printf("   (%3d)\n",m_iaData[z]);
	}
	m_IF.printf("\n");
}


void CBQBBitSet::DumpPlain() const {
	
	int z, z2, ti;

	if (m_iExcessBits == 0)
		ti = (int)m_iaData.size();
	else
		ti = (int)m_iaData.size()-1;
	for (z=0;z<ti;z++)
		for (z2=7;z2>=0;z2--)
			m_IF.printf("%d",(m_iaData[z]&(1<<(7-z2)))!=0?1:0);
	for (z2=m_iExcessBits-1;z2>=0;z2--)
		m_IF.printf("%d",(m_iaData[z]&(1<<(m_iExcessBits-z2-1)))!=0?1:0);
}


void CBQBBitSet::ExportToFile(FILE *a) const {

	int z;
//	unsigned char c;

	for (z=0;z<(int)m_iaData.size()/4096;z++)
		(void)!fwrite(&m_iaData[z*4096],4096,1,a);
	if ((m_iaData.size()%4096) != 0)
		(void)!fwrite(&m_iaData[z*4096],m_iaData.size()%4096,1,a);

/*	for (z=0;z<(int)m_iaData.size();z++) {
		c = m_iaData[z];
		(void)!fwrite(&c,1,1,a);
	}*/
}


bool CBQBBitSet::ImportFromFile(FILE *a, int bytes) {

	int z;
//	unsigned char c;

	m_iaData.resize(bytes);
//	m_iaData.clear();
	m_iExcessBits = 0;
	m_iReadPosBytes = 0;
	m_iReadPosExcess = 0;

	for (z=0;z<bytes/4096;z++) {

		(void)!fread(&m_iaData[z*4096],4096,1,a);

		if (feof(a))
			return false;
	}

	if ((bytes%4096) != 0)
		(void)!fread(&m_iaData[z*4096],bytes%4096,1,a); 

/*	for (z=0;z<bytes;z++) {
		(void)!fread(&c,1,1,a);
		m_iaData.push_back(c);
	}*/

	if (feof(a))
		return false;

	return true;
}


