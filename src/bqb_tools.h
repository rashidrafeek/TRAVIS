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



#ifndef BQB_TOOLS_H
#define BQB_TOOLS_H


// This must always be the first include directive
#include "bqb_config.h"

#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <cstdlib>
#include "bqb_largeinteger.h"



class CBQBInterface;



#ifndef M_LOG2E
  #define M_LOG2E (1.44269504088896340736) // log2(e)
#endif


#define C_RUNA          (0x7FFFFFFF)
#define C_RUNB          (0x7FFFFFFE)
#define C_FULLNUMBER    (0x7FFFFFFD)
#define C_SPLIT         (0x7FFFFFFC)
#define C_OVERFLOW      (0x7FFFFFFB)
#define C_UNDERFLOW     (0x7FFFFFFA)

#define C_MIN_RESERVED  (0x7FFFFFFA)
#define C_MAX_RESERVED  (0x7FFFFFFF)


int GetAtomOrd(const char *s);

const char* GetAtomOrdLabel( int ord );


#define UNUSED(a) ((void)a)


#define Pi (3.1415926535897932385)

#ifndef M_PI
	#define M_PI (3.1415926535897932385)
#endif

#define MIN(a,b) (((a)<=(b)) ? (a) : (b))
#define MAX(a,b) (((a)>=(b)) ? (a) : (b))

#define MIN3(a,b,c) MIN(MIN(a,b),c)
#define MAX3(a,b,c) MAX(MAX(a,b),c)

#define MIN4(a,b,c,d) MIN(MIN(a,b),MIN(c,d))
#define MAX4(a,b,c,d) MAX(MAX(a,b),MAX(c,d))


#define SQR(x)   ((x)*(x))


#ifdef _MSC_VER // Visual Studio (at least old versions) do not have std::abs
	#define bqbabs(a) (abs(a))
#else
	#define bqbabs(a) (std::abs(a))
#endif


typedef struct {
	unsigned long m_iN;
	unsigned long m_iH;
} THilbertPair;



int CompareHilbertPair( const void *arg1, const void *arg2 );



class CBQBStatisticsSample {
public:

	CBQBStatisticsSample() : m_lSize(0), m_lOverhead(0), m_lAlphabet(0), m_lHuffmanTables(0),
		m_lTableSwitch(0), m_lHuffmanData(0) {
	}

	void reset() {
		m_lSize = 0;
		m_lOverhead = 0;
		m_lAlphabet = 0;
		m_lHuffmanTables = 0;
		m_lTableSwitch = 0;
		m_lHuffmanData = 0;
	}

	BQBLargeInteger m_lSize;
	BQBLargeInteger m_lOverhead;
	BQBLargeInteger m_lAlphabet;
	BQBLargeInteger m_lHuffmanTables;
	BQBLargeInteger m_lTableSwitch;
	BQBLargeInteger m_lHuffmanData;
};



class CBQBStatistics {
public:

	explicit CBQBStatistics(CBQBInterface &i) : m_IF(i) {
	}

	void PushStatistics();
	void PopStatistics();
	void PopDiffStatistics();
	void PopIgnoreStatistics();
	void ResetStatistics();

	CBQBStatisticsSample m_oStat;
	std::vector<CBQBStatisticsSample> m_oaStatStack;

private:
	CBQBInterface &m_IF;
};



// This helper class encapsulates all functions
// which need to print to the screen (via CBQBInterface)
class CBQBTools {
public:

	explicit CBQBTools(CBQBInterface &i) : m_IF(i) {
	}

	FILE* BQBOpenFileWrite(const char *buf, bool text);

	long FindIndexPosition(FILE *a);

private:
	CBQBInterface &m_IF;
};



int bqb_strcmp_nocase(const char *s1, const char *s2);

bool BQBFileExist(const char *s);

void BQBFormatTime(unsigned long eta, char *buf);


inline int MoveToFrontValue(std::vector<int> &ia, int val) {

	int rtmp, rll_i, rtmp2;
	int* ryy_j;

	rtmp  = ia[1];
	ia[1] = ia[0];
	ryy_j = &(ia[1]);
	rll_i = val;
	while (rll_i != rtmp) {
		ryy_j++;
		rtmp2  = rtmp;
		rtmp   = *ryy_j;
		*ryy_j = rtmp2;
	}
	ia[0] = rtmp;
	return (int)(ryy_j - &(ia[0]));
}


inline void MoveToFrontIndex(std::vector<int> &ia, int idx) {
/*	int ti = ia[idx];
	ia.insert(ia.begin(),ti);
	ia.erase(ia.begin()+idx+1);*/
	int z, t;

	t = ia[idx];

	for (z=idx-1;z>=0;z--)
		ia[z+1] = ia[z];
	ia[0] = t;
}


inline bool IsAllIdentical(const std::vector<int> &ia) {

	unsigned int z;


	if (ia.size() < 2)
		return true;

	for (z=1;z<ia.size();z++)
		if (ia[z] != ia[0])
			return false;

	return true;
}


inline void BuildRunTable(const std::vector<int> &ia, std::vector<int> &runtable) {
	
	int i, v, j, z;
	bool b;

	runtable.resize(ia.size());

	i = 0;
	j = 0;
	v = ia[0];
	b = false;
	while (true) {
		i++;
		if (i >= (int)ia.size()) {
			i = 0;
			b = true;
		}
		if (ia[i] != v) {
			if (b) {
				for (z=j;z<(int)ia.size();z++)
					runtable[z] = (int)(i-z+ia.size());
				for (z=0;z<i;z++)
					runtable[z] = i-z;
				break;
			} else {
				for (z=j;z<i;z++)
					runtable[z] = i-z;
				j = i;
				v = ia[i];
			}
		}
	}
}


inline double mylog2(const double x) {
	return log(x) * M_LOG2E;
}


inline int positive_modulo(int i, unsigned int n) {
    return (i % n + n) % n;
}


/*
inline int positive_modulo(int a, int b) {
	const int result = a % b;
	return (result >= 0) ? result : (result + b);
}
*/


inline unsigned long reverse_bit_order(unsigned long u, int l) {

	unsigned long r;
	int z;

	r = 0;
	for (z=0;z<l;z++)
		if ((u&(((unsigned long)1)<<z)) != 0)
			r |= ((unsigned long)1)<<(l-z-1);

	return r;
}


inline double pow10(int i) {

	static double table[29] = {
		1.0E-14,
		1.0E-13,
		1.0E-12,
		1.0E-11,
		1.0E-10,
		1.0E-9,
		1.0E-8,
		1.0E-7,
		1.0E-6,
		1.0E-5,
		1.0E-4,
		1.0E-3,
		1.0E-2,
		1.0E-1,
		1.0,
		1.0E1,
		1.0E2,
		1.0E3,
		1.0E4,
		1.0E5,
		1.0E6,
		1.0E7,
		1.0E8,
		1.0E9,
		1.0E10,
		1.0E11,
		1.0E12,
		1.0E13,
		1.0E14
	};
	return table[i+14];
}


inline long pow10i(int i) {

	static int table[10] = {
		1,
		10,
		100,
		1000,
		10000,
		100000,
		1000000,
		10000000,
		100000000,
		1000000000
	};
	return table[i];
}


inline unsigned int pow2i(int i) {

	static unsigned int table[32] = {
		1,
		2,
		4,
		8,
		16,
		32,
		64,
		128,
		256,
		512,
		1024,
		2048,
		4096,
		8192,
		16384,
		32768,
		65536,
		131072,
		262144,
		524288,
		1048576,
		2097152,
		4194304,
		8388608,
		16777216,
		33554432,
		67108864,
		134217728,
		268435456,
		536870912,
		1073741824,
		2147483648
	};
	return table[i];
}



inline int log10i(long i) {

	if (i == 0)
		return 0;
	else if (bqbabs(i) <= 10)
		return 1;
	else if (bqbabs(i) <= 100)
		return 2;
	else if (bqbabs(i) <= 1000)
		return 3;
	else if (bqbabs(i) <= 10000)
		return 4;
	else if (bqbabs(i) <= 100000)
		return 5;
	else if (bqbabs(i) <= 1000000)
		return 6;
	else if (bqbabs(i) <= 10000000)
		return 7;
	else if (bqbabs(i) <= 100000000)
		return 8;
	else if (bqbabs(i) <= 1000000000)
		return 9;
	else
		return 10;
}


inline int log2i(long i) {

	if (i == 0)
		return 0;
	else if (bqbabs(i) <= 2)
		return 1;
	else if (bqbabs(i) <= 4)
		return 2;
	else if (bqbabs(i) <= 8)
		return 3;
	else if (bqbabs(i) <= 16)
		return 4;
	else if (bqbabs(i) <= 32)
		return 5;
	else if (bqbabs(i) <= 64)
		return 6;
	else if (bqbabs(i) <= 128)
		return 7;
	else if (bqbabs(i) <= 256)
		return 8;
	else if (bqbabs(i) <= 512)
		return 9;
	else if (bqbabs(i) <= 1024)
		return 10;
	else if (bqbabs(i) <= 2048)
		return 11;
	else if (bqbabs(i) <= 4096)
		return 12;
	else if (bqbabs(i) <= 8192)
		return 13;
	else if (bqbabs(i) <= 16384)
		return 14;
	else if (bqbabs(i) <= 32768)
		return 15;
	else if (bqbabs(i) <= 65536)
		return 16;
	else if (bqbabs(i) <= 131072)
		return 17;
	else if (bqbabs(i) <= 262144)
		return 18;
	else if (bqbabs(i) <= 524288)
		return 19;
	else if (bqbabs(i) <= 1048576)
		return 20;
	else if (bqbabs(i) <= 2097152)
		return 21;
	else if (bqbabs(i) <= 4194304)
		return 22;
	else if (bqbabs(i) <= 8388608)
		return 23;
	else if (bqbabs(i) <= 16777216)
		return 24;
	else if (bqbabs(i) <= 33554432)
		return 25;
	else if (bqbabs(i) <= 67108864)
		return 26;
	else if (bqbabs(i) <= 134217728)
		return 27;
	else if (bqbabs(i) <= 268435456)
		return 28;
	else if (bqbabs(i) <= 536870912)
		return 29;
	else if (bqbabs(i) <= 1073741824)
		return 30;
	else
		return 31;
}


inline bool ExpMantisEqual(int exp1, int mantis1, int exp2, int mantis2) {

	if (exp1 > exp2)
		mantis1 *= pow10i(exp1-exp2);
	else if (exp2 > exp1)
		mantis2 *= pow10i(exp2-exp1);

//	printf("@ %d vs %d\n",mantis1,mantis2);

	return (mantis1 == mantis2);
}


inline void fprintf_expo(FILE *a, long mantis, int signi, int expo) {

	fprintf(
		a,
		" %c0.%0*luE%c%02d",
		(mantis>=0)?' ':'-',
		signi,
		(unsigned long)bqbabs(mantis),
		((expo+signi)>=0)?'+':'-',
		bqbabs(expo+signi)
	);
}


inline long FloatToFixed(double f, int digits) {

	return (long)floor(f*pow10i(digits)+0.5);
}


inline double FixedToFloat(long l, int digits) {

	return (double)l*pow10(-digits);
}


/*
inline void DumpFixed(long f, int digits) {

	m_IF.printf("%ld",f/pow10i(digits));
	m_IF.printf(".");
	m_IF.printf("%0*ld",digits,f%pow10i(digits));
}
*/


inline int RoundCenter(double f) {

	if (f >= 0)
		return (int)floor(f);
	else
		return (int)ceil(f);
}


void PushNumeration(std::vector<int> &ia, int n, int s1, int s2);


bool bqbisinteger(const char *s);



#endif


