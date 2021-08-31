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



#ifndef BQB_BITSET_H
#define BQB_BITSET_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"
#include <vector>



class CBQBInterface;



class CBQBBitSet {
public:

	explicit CBQBBitSet(CBQBInterface &i)
		: m_iExcessBits(0), m_iReadPosBytes(0), m_iReadPosExcess(0), m_IF(i) {
	}


	CBQBBitSet(const CBQBBitSet &set)
		: m_iExcessBits(set.m_iExcessBits), m_iReadPosBytes(set.m_iReadPosBytes),
			m_iReadPosExcess(set.m_iReadPosExcess), m_IF(set.m_IF) {
		m_iaData.assign(set.m_iaData.begin(),set.m_iaData.end());
	}


	explicit CBQBBitSet(CBQBBitSet *set)
		: m_iExcessBits(set->m_iExcessBits), m_IF(set->m_IF) {
		m_iaData.assign(set->m_iaData.begin(),set->m_iaData.end());
	}


	CBQBBitSet& operator=(const CBQBBitSet &rhs) {
		m_iExcessBits =  rhs.m_iExcessBits;
		m_iReadPosBytes = rhs.m_iReadPosBytes;
		m_iReadPosExcess = rhs.m_iReadPosExcess;
		m_iaData.assign(rhs.m_iaData.begin(),rhs.m_iaData.end());
		return *this;
	}


	int GetLength() const {
		if (m_iExcessBits == 0)
			return (int)(8*m_iaData.size());
		else
			return (int)(8*(m_iaData.size()-1) + m_iExcessBits);
	}


	int GetByteLength() const {
		return (int)m_iaData.size();
	}


	void ResetReadPos() {
		m_iReadPosBytes = 0;
		m_iReadPosExcess = 0;
	}


	void Dump() const;

	void DumpPlain() const;

	void ExportToFile(FILE *a) const;

	bool ImportFromFile(FILE *a, int bytes);

	std::string ExportKey();

	bool ImportKey(std::string s);


	void Clear() {
		m_iExcessBits = 0;
		m_iReadPosBytes = 0;
		m_iReadPosExcess = 0;
		m_iaData.clear();
	}


	void WriteBit(unsigned char i) {
		if (m_iExcessBits == 0)
			m_iaData.push_back(i);
		else
			m_iaData[m_iaData.size()-1] |= i<<m_iExcessBits;

		m_iExcessBits++;

		if (m_iExcessBits == 8)
			m_iExcessBits = 0;
	}


	void WriteBits(unsigned long i, int bits);


	void WriteSignedBits(long i, int bits);


	bool ReadBit();


	bool ReadBitOrZero();


	void WriteBitsFloat(float f) {
		unsigned char *uc;
		int z;
		uc = reinterpret_cast<unsigned char*>(&f);
		for (z=0;z<4;z++)
			WriteBits(uc[z],8);
	}


	void WriteBitsDouble(double f) {
		unsigned char *uc;
		int z;
		uc = reinterpret_cast<unsigned char*>(&f);
		for (z=0;z<8;z++)
			WriteBits(uc[z],8);
	}


	float ReadBitsFloat() {
		unsigned char uc[4];
		char *p;

		int z;
		for (z=0;z<4;z++)
			uc[z] = (unsigned char)ReadBitsInteger(8);
		p = (char*)uc; // This is to comply with strict-aliasing rules (char* may alias to anything)
		return *reinterpret_cast<float*>(p);
	}


	double ReadBitsDouble() {
		unsigned char uc[8];
		int z;
		char *p;

		for (z=0;z<8;z++)
			uc[z] = (unsigned char)ReadBitsInteger(8);
		p = (char*)uc; // This is to comply with strict-aliasing rules (char* may alias to anything)
		return *reinterpret_cast<double*>(p);
	}


	void WriteBits(CBQBBitSet *set) {
		int z, z2;
		if (set->m_iExcessBits == 0) {
			for (z=0;z<(int)set->m_iaData.size();z++)
				for (z2=0;z2<8;z2++)
					WriteBit((set->m_iaData[z]&(((unsigned long)1)<<z2))!=0?1:0);
		} else {
			for (z=0;z<(int)set->m_iaData.size()-1;z++)
				for (z2=0;z2<8;z2++)
					WriteBit((set->m_iaData[z]&(((unsigned long)1)<<z2))!=0?1:0);
			for (z2=0;z2<set->m_iExcessBits;z2++)
				WriteBit((set->m_iaData[set->m_iaData.size()-1]&(((unsigned long)1)<<z2))!=0?1:0);
		}
	}


	void WriteProgressiveUnsignedChar(unsigned char i) {
		if (i < 4) {
			WriteBit(0);
			WriteBits(i,2);
		} else if (i < 16) {
			WriteBit(1);
			WriteBit(0);
			WriteBits(i,4);
		} else if (i < 64) {
			WriteBit(1);
			WriteBit(1);
			WriteBit(0);
			WriteBits(i,6);
		} else {
			WriteBit(1);
			WriteBit(1);
			WriteBit(1);
			WriteBit(0);
			WriteBits(i,8);
		}
	}


	unsigned char ReadProgressiveUnsignedChar() {
		int z = 1;
		while (ReadBit())
			z++;
		return (unsigned char)ReadBitsInteger(z*2);
	}


	void WriteProgressiveUnsignedInt(unsigned int i) {
		if (i < 256) {
			WriteBit(0);
			WriteBits(i,8);
		} else if (i < 65536) {
			WriteBit(1);
			WriteBit(0);
			WriteBits(i,16);
		} else if (i < 16777216) {
			WriteBit(1);
			WriteBit(1);
			WriteBit(0);
			WriteBits(i,24);
		} else {
			WriteBit(1);
			WriteBit(1);
			WriteBit(1);
			WriteBit(0);
			WriteBits(i,32);
		}
	}


	unsigned int ReadProgressiveUnsignedInt() {
		int z = 1;
		while (ReadBit())
			z++;
		return ReadBitsInteger(z*8);
	}


	unsigned long ReadBitsInteger(int bits) {
		unsigned long u = 0;
		for (int z=0;z<bits;z++)
			if (ReadBit())
				u |= ((unsigned long)1)<<z;
		return u;
	}


	long ReadBitsSignedInteger(int bits) {
		unsigned long u = 0;
		for (int z=0;z<bits-1;z++)
			if (ReadBit())
				u |= ((unsigned long)1)<<z;
		if (ReadBit())
			return -((long)u);
		else
			return (long)u;
	}


	int GetReadPos() {
		return m_iReadPosBytes*8 + m_iReadPosExcess;
	}


	std::vector<unsigned char> m_iaData;
	int m_iExcessBits;
	int m_iReadPosBytes;
	int m_iReadPosExcess;

private:
	CBQBInterface &m_IF;
};



#endif


