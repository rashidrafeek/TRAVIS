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

#include "bqb_tools.h"
#include "bqb_hufftree.h"
#include <algorithm>
#include <math.h>
#include "bqb_integerengine.h"
#include "bqb_interface.h"


const char *GetRevisionInfo_bqb_hufftree(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_hufftree() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



void CBQBHuffmanTree::Init(int alphasize) {

	int z;

	m_iaLengths.resize(alphasize);
	m_oaBitStrings.resize(alphasize);
	m_iaFrequencies.resize(alphasize);
	for (z=0;z<(int)m_iaFrequencies.size();z++)
		m_iaFrequencies[z] = 0;
}


bool SORT_SymbolsAsc(CBQBHuffmanSymbol *s1, CBQBHuffmanSymbol *s2) {

	return (s1->m_iFrequency < s2->m_iFrequency);
}


bool SORT_SymbolsDesc(CBQBHuffmanSymbol *s1, CBQBHuffmanSymbol *s2) {

	return (s1->m_iFrequency > s2->m_iFrequency);
}


bool SORT_SymbolsAlpha(CBQBHuffmanSymbol *s1, CBQBHuffmanSymbol *s2) {

	return (s1->m_iSymbol < s2->m_iSymbol);
}


bool SORT_SymbolsCanonical(CBQBHuffmanSymbol *s1, CBQBHuffmanSymbol *s2) {

	if (s1->m_pBitString->GetLength() < s2->m_pBitString->GetLength())
		return true;
	else if (s1->m_pBitString->GetLength() > s2->m_pBitString->GetLength())
		return false;
	else
		return (s1->m_iSymbol < s2->m_iSymbol);
}


void CBQBHuffmanEstimator::Init(int alphasize) {

	m_iAlphaSize = alphasize;

	m_iaLength.resize(alphasize);
	m_iaTempHeap.resize(alphasize+2);
	m_iaTempFrequencies.resize(alphasize*2);
	m_iaTempParent.resize(alphasize*2);
}


// These macros are adapted from the BZIP2 source code

#define WEIGHTOF(zz0)  ((zz0) & 0xffffff00)
#define DEPTHOF(zz1)   ((zz1) & 0x000000ff)
#define MYMAX(zz2,zz3) ((zz2) > (zz3) ? (zz2) : (zz3))


#define ADDWEIGHTS(zw1,zw2)               \
   (WEIGHTOF(zw1)+WEIGHTOF(zw2)) |        \
   (1 + MYMAX(DEPTHOF(zw1),DEPTHOF(zw2)))


#define UPHEAP(z)                                                                  \
{                                                                                  \
   unsigned int zz, tmp;                                                           \
   zz = z; tmp = m_iaTempHeap[zz];                                                 \
   while (m_iaTempFrequencies[tmp] < m_iaTempFrequencies[m_iaTempHeap[zz >> 1]]) { \
      m_iaTempHeap[zz] = m_iaTempHeap[zz >> 1];                                    \
      zz >>= 1;                                                                    \
   }                                                                               \
   m_iaTempHeap[zz] = tmp;                                                         \
}


#define DOWNHEAP(z)                                                                         \
{                                                                                           \
   unsigned int zz, yy, tmp;                                                                \
   zz = z; tmp = m_iaTempHeap[zz];                                                          \
   while (true) {                                                                           \
      yy = zz << 1;                                                                         \
      if (yy > nHeap) break;                                                                \
      if (yy < nHeap &&                                                                     \
          m_iaTempFrequencies[m_iaTempHeap[yy+1]] < m_iaTempFrequencies[m_iaTempHeap[yy]])  \
         yy++;                                                                              \
      if (m_iaTempFrequencies[tmp] < m_iaTempFrequencies[m_iaTempHeap[yy]]) break;          \
      m_iaTempHeap[zz] = m_iaTempHeap[yy];                                                  \
      zz = yy;                                                                              \
   }                                                                                        \
   m_iaTempHeap[zz] = tmp;                                                                  \
}


void CBQBHuffmanEstimator::BuildBitLengthTable(const std::vector<int> &ia) {

	unsigned int nNodes, nHeap, n1, n2, i, j, k;

	for (i = 0; i < (unsigned int)m_iAlphaSize; i++)
		m_iaTempFrequencies[i+1] = ia[i] << 8;
//		m_iaTempFrequencies[i+1] = (ia[i] == 0 ? 1 : ia[i]) << 8;

	nNodes = m_iAlphaSize;
	nHeap = 0;

	m_iaTempHeap[0] = 0;
	m_iaTempFrequencies[0] = 0;
	m_iaTempParent[0] = -2;

	for (i = 1; i <= (unsigned int)m_iAlphaSize; i++) {
		m_iaTempParent[i] = -1;
		nHeap++;
		m_iaTempHeap[nHeap] = i;
		UPHEAP(nHeap);
	}

	while (nHeap > 1) {
		n1 = m_iaTempHeap[1]; m_iaTempHeap[1] = m_iaTempHeap[nHeap]; nHeap--; DOWNHEAP(1);
		n2 = m_iaTempHeap[1]; m_iaTempHeap[1] = m_iaTempHeap[nHeap]; nHeap--; DOWNHEAP(1);
		nNodes++;
		m_iaTempParent[n1] = m_iaTempParent[n2] = nNodes;
		m_iaTempFrequencies[nNodes] = ADDWEIGHTS(m_iaTempFrequencies[n1], m_iaTempFrequencies[n2]);
		m_iaTempParent[nNodes] = -1;
		nHeap++;
		m_iaTempHeap[nHeap] = nNodes;
		UPHEAP(nHeap);
	}

	for (i = 1; i <= (unsigned int)m_iAlphaSize; i++) {
		j = 0;
		k = i;
		while (m_iaTempParent[k] >= 0) { 
			k = m_iaTempParent[k];
			j++;
		}
		m_iaLength[i-1] = j;
	}
}


unsigned int CBQBHuffmanEstimator::EstimateBitLength(const std::vector<int> &ia) {

	unsigned int nNodes, nHeap, n1, n2, i, j, k, ret;

	for (i = 0; i < (unsigned int)m_iAlphaSize; i++)
		m_iaTempFrequencies[i+1] = ia[i] << 8;
//		m_iaTempFrequencies[i+1] = (ia[i] == 0 ? 1 : ia[i]) << 8;

	nNodes = m_iAlphaSize;
	nHeap = 0;

	m_iaTempHeap[0] = 0;
	m_iaTempFrequencies[0] = 0;
	m_iaTempParent[0] = -2;

	for (i = 1; i <= (unsigned int)m_iAlphaSize; i++) {
		m_iaTempParent[i] = -1;
		nHeap++;
		m_iaTempHeap[nHeap] = i;
		UPHEAP(nHeap);
	}

	while (nHeap > 1) {
		n1 = m_iaTempHeap[1]; m_iaTempHeap[1] = m_iaTempHeap[nHeap]; nHeap--; DOWNHEAP(1);
		n2 = m_iaTempHeap[1]; m_iaTempHeap[1] = m_iaTempHeap[nHeap]; nHeap--; DOWNHEAP(1);
		nNodes++;
		m_iaTempParent[n1] = m_iaTempParent[n2] = nNodes;
		m_iaTempFrequencies[nNodes] = ADDWEIGHTS( m_iaTempFrequencies[n1], m_iaTempFrequencies[n2] );
		m_iaTempParent[nNodes] = -1;
		nHeap++;
		m_iaTempHeap[nHeap] = nNodes;
		UPHEAP(nHeap);
	}

	ret = 0;
	for (i = 0; i < (unsigned int)m_iAlphaSize; i++) {
		if (ia[i] == 0)
			continue;
		j = 0;
		k = i+1;
		while (m_iaTempParent[k] >= 0) { 
			k = m_iaTempParent[k];
			j++;
		}
		ret += j * ia[i];
	}

	return ret;
}


unsigned int CBQBHuffmanEstimator::EstimateBitLength(const std::vector<int> &ia, int i2) {

	std::vector<int> tia;

	tia.assign(ia.begin(),ia.end());

	tia[i2]++;

	return EstimateBitLength(tia);
}


unsigned int CBQBHuffmanEstimator::EstimateBitLength(const std::vector<int> &ia, const std::vector<int> &ia2) {

	std::vector<int> tia;
	int z;

	tia.assign(ia.begin(),ia.end());

	for (z=0;z<(int)tia.size();z++)
		tia[z] += ia2[z];

	return EstimateBitLength(tia);
}


unsigned int CBQBHuffmanEstimator::EstimateBitLength(const std::vector<int> &ia, const std::vector<int> &ia2, const std::vector<int> &ia3) {

	std::vector<int> tia;
	int z;

	tia.assign(ia.begin(),ia.end());

	for (z=0;z<(int)tia.size();z++)
		tia[z] += ia2[z] + ia3[z];

	return EstimateBitLength(tia);
}


void CBQBHuffmanTree::BuildTree(bool canonical, bool showsymbols) {

	int z, ti, c;
	CBQBHuffmanSymbol *sym;
	unsigned long l;

	m_bCanonical = canonical;

	m_oaSymbols.clear();
	if (m_pTree != NULL)
		delete m_pTree;
	m_pTree = NULL;

	for (z=0;z<(int)m_iaFrequencies.size();z++) {
		if (m_iaFrequencies[z] != 0) {
			sym = new CBQBHuffmanSymbol();
			sym->m_iSymbol = z;
			sym->m_iFrequency = m_iaFrequencies[z];
//			sym->m_iFrequency = (m_iaFrequencies[z]==0)?1:m_iaFrequencies[z];
			sym->m_bVirtual = false;
			sym->m_iDepth = 0;
			m_oaSymbols.push_back(sym);
		}
	}

	if (m_oaSymbols.size() == 0)
		return;

	std::sort(m_oaSymbols.begin(),m_oaSymbols.end(),SORT_SymbolsAsc);

	m_oaTempSymbols.assign(m_oaSymbols.begin(),m_oaSymbols.end());

	while (m_oaTempSymbols.size() > 1) {
		sym = new CBQBHuffmanSymbol();
		sym->m_iSymbol = 0;
		sym->m_bVirtual = true;
		sym->m_iFrequency = m_oaTempSymbols[0]->m_iFrequency + m_oaTempSymbols[1]->m_iFrequency;
		sym->m_oaChildren.push_back(m_oaTempSymbols[0]);
		sym->m_oaChildren.push_back(m_oaTempSymbols[1]);
		sym->m_iDepth = MAX( m_oaTempSymbols[0]->m_iDepth, m_oaTempSymbols[1]->m_iDepth ) + 1;
		m_oaTempSymbols.erase(m_oaTempSymbols.begin(),m_oaTempSymbols.begin()+2);
		for (z=0;z<(int)m_oaTempSymbols.size();z++) {
			if (m_oaTempSymbols[z]->m_iFrequency > sym->m_iFrequency)
				break;
			if ((m_oaTempSymbols[z]->m_iFrequency == sym->m_iFrequency) && (m_oaTempSymbols[z]->m_iDepth > sym->m_iDepth))
				break;
		}
		m_oaTempSymbols.insert(m_oaTempSymbols.begin()+z,sym);
	}

	m_pTree = m_oaTempSymbols[0];
	m_oaTempSymbols.clear();

	m_pTree->m_pBitString = new CBQBBitSet(m_IF);

	m_pTree->REC_BuildBitstrings(0);

	if (m_bCanonical) {

		if (m_pTree->m_oaChildren.size() == 0) // Hack 02.04.2017
			m_pTree->m_pBitString->WriteBits(0,1);

		std::sort(m_oaSymbols.begin(),m_oaSymbols.end(),SORT_SymbolsCanonical);
		ti = 0;
		l = 0;
		for (z=0;z<(int)m_oaSymbols.size();z++) {
			if (m_oaSymbols[z]->m_pBitString->GetLength() > ti) {
				l <<= m_oaSymbols[z]->m_pBitString->GetLength() - ti;
				ti = m_oaSymbols[z]->m_pBitString->GetLength();
			}
			m_oaSymbols[z]->m_pBitString->Clear();
			m_oaSymbols[z]->m_pBitString->WriteBits(reverse_bit_order(l,ti),ti);
			l++;
		}
		std::sort(m_oaSymbols.begin(),m_oaSymbols.end(),SORT_SymbolsAlpha);
	}

	for (z=0;z<(int)m_iaLengths.size();z++)
		m_iaLengths[z] = 0;

	m_iMaxBitLength = 0;
	for (z=0;z<(int)m_oaSymbols.size();z++) {
		m_iaLengths[m_oaSymbols[z]->m_iSymbol] = m_oaSymbols[z]->m_pBitString->GetLength();
		if (m_iMaxBitLength < m_iaLengths[m_oaSymbols[z]->m_iSymbol])
			m_iMaxBitLength = m_iaLengths[m_oaSymbols[z]->m_iSymbol];
		m_oaBitStrings[m_oaSymbols[z]->m_iSymbol] = m_oaSymbols[z]->m_pBitString;
	}

	if (showsymbols) {
		std::sort(m_oaSymbols.begin(),m_oaSymbols.end(),SORT_SymbolsDesc);
		m_IF.printf("*** Output of Huffman Table ***\n");
		ti = 0;
		c = 0;
		for (z=0;z<(int)m_oaSymbols.size();z++) {
			m_IF.printf("    %6d  (%9d):  ",m_oaSymbols[z]->m_iSymbol,m_oaSymbols[z]->m_iFrequency);
			m_oaSymbols[z]->m_pBitString->DumpPlain();
			m_IF.printf(" (%d)\n",m_oaSymbols[z]->m_pBitString->GetLength());
			ti += m_oaSymbols[z]->m_iFrequency * m_oaSymbols[z]->m_pBitString->GetLength();
			c += m_oaSymbols[z]->m_iFrequency;
		}
		m_IF.printf("      %lu symbol types, %d symbols, %.4f bits per entry, %.4f MiB per frame.\n",(unsigned long)m_oaSymbols.size(),c,((double)ti)/m_pTree->m_iFrequency,ti/8.0/1024.0/1024.0);
		m_IF.printf("\n");
		m_IF.printf("\n");
	}
}


void CBQBHuffmanTree::BuildPrelimTree() {

	int z;
	CBQBHuffmanSymbol *sym;

	m_oaSymbols.clear();
	if (m_pTree != NULL)
		delete m_pTree;
	m_pTree = NULL;

	for (z=0;z<(int)m_iaFrequencies.size();z++) {
		sym = new CBQBHuffmanSymbol();
		sym->m_iSymbol = z;
		sym->m_iFrequency = m_iaFrequencies[z];
//		sym->m_iFrequency = (m_iaFrequencies[z]==0)?1:m_iaFrequencies[z];
		sym->m_bVirtual = false;
		sym->m_iDepth = 0;
		m_oaSymbols.push_back(sym);
	}

	std::sort(m_oaSymbols.begin(),m_oaSymbols.end(),SORT_SymbolsAsc);

	m_oaTempSymbols.assign(m_oaSymbols.begin(),m_oaSymbols.end());

	while (m_oaTempSymbols.size() > 1) {
		sym = new CBQBHuffmanSymbol();
		sym->m_iSymbol = 0;
		sym->m_bVirtual = true;
		sym->m_iFrequency = m_oaTempSymbols[0]->m_iFrequency + m_oaTempSymbols[1]->m_iFrequency;
		sym->m_iDepth = MAX( m_oaTempSymbols[0]->m_iDepth, m_oaTempSymbols[1]->m_iDepth ) + 1;
		sym->m_oaChildren.push_back(m_oaTempSymbols[0]);
		sym->m_oaChildren.push_back(m_oaTempSymbols[1]);
		m_oaTempSymbols.erase(m_oaTempSymbols.begin(),m_oaTempSymbols.begin()+2);
		for (z=0;z<(int)m_oaTempSymbols.size();z++) {
			if (m_oaTempSymbols[z]->m_iFrequency > sym->m_iFrequency)
				break;
			if ((m_oaTempSymbols[z]->m_iFrequency == sym->m_iFrequency) && (m_oaTempSymbols[z]->m_iDepth > sym->m_iDepth))
				break;
		}
		m_oaTempSymbols.insert(m_oaTempSymbols.begin()+z,sym);
	}

	m_pTree = m_oaTempSymbols[0];
	m_oaTempSymbols.clear();

	REC_BuildBitstringLengths(0,m_pTree);
}


void CBQBHuffmanSymbol::REC_BuildBitstrings(int depth) {

	if (m_oaChildren.size() == 2) {
		m_oaChildren[0]->m_pBitString = new CBQBBitSet(m_pBitString);
		m_oaChildren[0]->m_pBitString->WriteBit(0);
		m_oaChildren[0]->REC_BuildBitstrings(depth+1);
		m_oaChildren[1]->m_pBitString = new CBQBBitSet(m_pBitString);
		m_oaChildren[1]->m_pBitString->WriteBit(1);
		m_oaChildren[1]->REC_BuildBitstrings(depth+1);
	}
}


void CBQBHuffmanTree::REC_BuildBitstringLengths(int depth, CBQBHuffmanSymbol *sym) {

	if (sym->m_oaChildren.size() == 2) {
		REC_BuildBitstringLengths(depth+1,sym->m_oaChildren[0]);
		REC_BuildBitstringLengths(depth+1,sym->m_oaChildren[1]);
	} else
//		if (sym->m_iFrequency == 0)
//			m_iaLengths[sym->m_iSymbol] = depth+5;
//		else
			m_iaLengths[sym->m_iSymbol] = depth;
}


int CBQBHuffmanTree::ExportTree(CBQBBitSet *bs, bool chr) {

	int b2, b2a, b2s, i, z, z2, z0, nz, orig, alt, alt2, alt20;
	int alt3, method, bborder, iborder, si, tc;
	std::vector<int> mtf, ia, ia2, border;
	CBQBHuffmanTree ht(m_IF);
	//CBQBBitSet bs2(m_IF);


	si = 0;

	if (m_bCanonical) { // Canonical

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf(">>> ExportTree Canonical Version >>>\n");

		si++;
		if (bs != NULL)
			bs->WriteBit(1);

		i = 0;
		nz = 0;
		for (z=0;z<(int)m_iaFrequencies.size();z++) {
			if (m_iaFrequencies[z] != 0) {
				i = z;
				nz++;
			}
		}
		i++;
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Symbol type count is %d (%d non-zero symbols).\n",i,nz);

		b2 = (int)ceil(mylog2(m_iMaxBitLength+1));
		b2a = (int)ceil(mylog2(m_iMaxBitLength+2));
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Max bit length is %d, using %d bits per length token.\n",m_iMaxBitLength,b2);

		ia.resize(i);
		for (z=0;z<i;z++) {
			if (m_iaFrequencies[z] != 0) {
				for (z2=0;z2<(int)m_oaSymbols.size();z2++)
					if (m_oaSymbols[z2]->m_iSymbol == z)
						break;
				if (z2 == (int)m_oaSymbols.size()) {
					m_IF.eprintf("CBQBHuffmanTree::ExportTree(): Internal error.\n");
					abort();
				}
				ia[z] = m_oaSymbols[z2]->m_pBitString->GetLength();
			} else
				ia[z] = 0;
		}

		mtf.resize(m_iMaxBitLength+1);
		for (z=0;z<=m_iMaxBitLength;z++)
			mtf[z] = z;

		ia2.clear();
		tc = 0;
		for (z=0;z<i;z++) {

			for (z2=0;z2<(int)mtf.size();z2++)
				if (mtf[z2] == ia[z])
					break;

			if (z2 != 0) {
				if (tc != 0)
					PushNumeration(ia2,tc,0,1);
				tc = 0;
				ia2.push_back(z2+1);

				MoveToFrontIndex(mtf,z2);
				//std::iter_swap(mtf.begin(),mtf.begin()+z2);

	//			mtf.insert(mtf.begin(),ia[z]);
	//			mtf.erase(mtf.begin()+z2+1);
			} else
				tc++;

	/*		if (z2 != 0) {
				ia2.push_back(z2);
				mtf.insert(mtf.begin(),ia[z]);
				mtf.erase(mtf.begin()+z2+1);
			} else
				ia2.push_back(0);*/

		}
		if (tc != 0)
			PushNumeration(ia2,tc,0,1);

		alt20 = 0;
		for (z=0;z<(int)ia2.size();z++) {
			if (ia2[z] == 0)
				alt20++;
			else
				alt20 += b2a+1;
		}

		border.resize(m_iMaxBitLength+1);
		bborder = 1000000000;
		iborder = -1;
		for (z0=0;z0<=m_iMaxBitLength;z0++) {
			alt2 = 0;
			for (z=0;z<(int)ia2.size();z++) {
				if (ia2[z] == 0)
					alt2++;
				else if ((int)ceil(mylog2(ia2[z]+1)) <= z0)
					alt2 += z0+2;
				else
					alt2 += b2+2;
			}
			if (alt2 < bborder) {
				bborder = alt2;
				iborder = z0;
			}
		}

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Method 2: Best border storage %d bits (border %d), no border %d bits.\n",bborder,iborder,alt20);

		alt2 = bborder;
		if (alt20 < alt2) {
			alt2 = alt20;
			iborder = 0;
		}

		if (i > 10) {
			alt3 = 0;
			ht.Init(m_iMaxBitLength+2);
			for (z=0;z<(int)ia2.size();z++)
				ht.m_iaFrequencies[ia2[z]]++;
			ht.BuildTree(true);
			alt3 = ht.ExportTree(NULL,chr);
			for (z=0;z<(int)ia2.size();z++)
				alt3 += ht.m_iaLengths[ia2[z]];
		} else {
			alt3 = 1000000000;
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Low symbol type count, not trying alt3.\n");
		}

		b2s = (int)ceil(mylog2(i+1));

		if (chr) {
			if (i < 64) {
				orig = b2*i + 13;
				alt = b2*nz + i + 13;
				alt3 += 7;
				if (iborder == 0)
					alt2 += 14;
				else
					alt2 += 17;
			} else {
				orig = b2*i + 16;
				alt = b2*nz + i + 16;
				alt3 += 10;
				if (iborder == 0)
					alt2 += 17;
				else
					alt2 += 20;
			}
		} else {
			if (i < 64) {
				orig = b2*i + 13;
				alt = b2*nz + i + 13;
				alt3 += 7;
				if (iborder == 0)
					alt2 += 14;
				else
					alt2 += 19;
			} else {
				orig = b2*i + b2s + 13;
				alt = b2*nz + i + b2s + 13;
				alt3 += 7+b2s;
				if (iborder == 0)
					alt2 += b2s + 14;
				else
					alt2 += b2s + 19;
			}
		}

		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			m_IF.printf("Original method (0):    %5d bits.\n",orig);
			m_IF.printf("Alternative method (1): %5d bits.\n",alt);
			m_IF.printf("Alternative method (2): %5d bits (border %d).\n",alt2,iborder);
			if (alt3 < 1e9)
				m_IF.printf("Alternative method (3): %5d bits.\n",alt3);
		}
		if ((orig <= alt) && (orig <= alt2) && (orig <= alt3))
			method = 0;
		else if ((alt <= orig) && (alt <= alt2) && (alt <= alt3))
			method = 1;
		else if ((alt2 <= orig) && (alt2 <= alt) && (alt2 <= alt3))
			method = 2;
		else 
			method = 3;

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("--> Using method %d.\n",method);
		si += 2;
		if (bs != NULL)
			bs->WriteBits(method,2);

		nz = si;

		if ((method == 0) || (method == 1)) {

			if (i < 64) {
				if (m_IF.IsPL(BQB_PL_DEBUG))
					m_IF.printf("Symbol type count is %d (%d non-zero symbols), using short storage.\n",i,nz);
				si += 7;
				if (bs != NULL) {
					bs->WriteBit(0);
					bs->WriteBits(i,6);
				}
			} else {
				if (chr) {
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("Symbol type count is %d (%d non-zero symbols), using 9 bits (char).\n",i,nz);
					si += 10;
					if (bs != NULL) {
						bs->WriteBit(1);
						bs->WriteBits(i,9);
					}
				} else {
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("Symbol type count is %d (%d non-zero symbols), using %d bits.\n",i,nz,b2s);
					si += 7 + b2s;
					if (bs != NULL) {
						bs->WriteBit(1);
						bs->WriteBits(b2s,6);
						bs->WriteBits(i,b2s);
					}
				}
			}

	//		if (!chr) {
				si += 6;
				if (bs != NULL)
					bs->WriteBits(b2,6);
	//		} else {
	//			si += 4;
	//			if (bs != NULL)
	//				bs->WriteBits(b2,4);
	//		}

			for (z=0;z<i;z++) {

				if (method == 1) {
					if (m_iaFrequencies[z] != 0) {
						si += b2+1;
						if (bs != NULL) {
							bs->WriteBit(1);
							bs->WriteBits(ia[z],b2);
						}
					} else {
						si++;
						if (bs != NULL)
							bs->WriteBit(0);
					}
				} else {
					si += b2;
					if (bs != NULL)
						bs->WriteBits(ia[z],b2);
				}

	/*			if (m_iaFrequencies[z] != 0) {
					if (method == 1) {
						si++;
						if (bs != NULL)
							bs->WriteBit(1);
					}
					si += b2;
					if (bs != NULL)
						bs->WriteBits(ia[z],b2);
				} else {
					if (method == 1) {
						si ++;
						if (bs != NULL)
							bs->WriteBit(0);
					} else {
						si += b2;
						if (bs != NULL)
							bs->WriteBits(0,b2);
					}
				}*/

			}

		} else if (method == 2) {

			b2s = (int)ceil(mylog2((double)ia2.size()+1));

			if (ia2.size() < 64) {
				if (m_IF.IsPL(BQB_PL_DEBUG))
					m_IF.printf("Symbol type count is %lu, using short storage.\n",(unsigned long)ia2.size());
				si += 7;
				if (bs != NULL) {
					bs->WriteBit(0);
					bs->WriteBits((unsigned long)ia2.size(),6);
				}
			} else {
				if (!chr) {
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("Transmission symbol type count is %lu, using %d bits.\n",(unsigned long)ia2.size(),b2s);
					si += 7 + b2s;
					if (bs != NULL) {
						bs->WriteBit(1);
						bs->WriteBits(b2s,6);
						bs->WriteBits((unsigned long)ia2.size(),b2s);
					}
				} else {
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("Transmission symbol type count is %lu, using 9 bits (char).\n",(unsigned long)ia2.size());
					si += 10;
					if (bs != NULL) {
						bs->WriteBit(1);
						bs->WriteBits((unsigned long)ia2.size(),9);
					}
				}
			}

	//		if (!chr) {
				si += 6;
				if (bs != NULL)
					bs->WriteBits(b2a,6);
	//		} else {
	//			si += 4;
	//			if (bs != NULL)
	//				bs->WriteBits(b2a,4);
	//		}

			if (iborder != 0) {
				si++;
				if (bs != NULL)
					bs->WriteBit(1);
				if (chr) {
					if (bs != NULL)
						bs->WriteBits(iborder,3);
					si += 3;
				} else {
					if (bs != NULL)
						bs->WriteBits(iborder,5);
					si += 5;
				}
			} else {
				si++;
				if (bs != NULL)
					bs->WriteBit(0);
			}

			for (z=0;z<(int)ia2.size();z++) {
				if (ia2[z] != 0) {
					si++;
					if (bs != NULL)
						bs->WriteBit(1);
					if ((int)ceil(mylog2(ia2[z]+1)) <= iborder) {
						if (iborder != 0) {
							si++;
							if (bs != NULL)
								bs->WriteBit(0);
						}
						si += iborder;
						if (bs != NULL)
							bs->WriteBits(ia2[z],iborder);
					} else {
						if (iborder != 0) {
							si++;
							if (bs != NULL)
								bs->WriteBit(1);
						}
						si += b2a;
						if (bs != NULL)
							bs->WriteBits(ia2[z],b2a);
					}
				} else {
					si++;
					if (bs != NULL)
						bs->WriteBit(0);
				}
			}

		} else { // Method 3
		
			b2s = (int)ceil(mylog2((double)ia2.size()+1));

			if (ia2.size() < 64) {
				if (m_IF.IsPL(BQB_PL_DEBUG))
					m_IF.printf("Symbol type count is %lu, using short storage.\n",(unsigned long)ia2.size());
			//	si += 7;
				if (bs != NULL) {
					bs->WriteBit(0);
					bs->WriteBits((unsigned long)ia2.size(),6);
				}
			} else {
				if (!chr) {
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("Transmission symbol type count is %lu, using %d bits.\n",(unsigned long)ia2.size(),b2s);
			//		si += 7 + b2s;
					if (bs != NULL) {
						bs->WriteBit(1);
						bs->WriteBits(b2s,6);
						bs->WriteBits((unsigned long)ia2.size(),b2s);
					}
				} else {
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("Transmission symbol type count is %lu, using 9 bits (char).\n",(unsigned long)ia2.size());
			//		si += 9;
					if (bs != NULL) {
						bs->WriteBit(1);
						bs->WriteBits((unsigned long)ia2.size(),9);
					}
				}
			}
			si += alt3;

			if (bs != NULL)
				ht.ExportTree(bs,chr);
			for (z=0;z<(int)ia2.size();z++) {
		//		si += ht.m_oaBitStrings[ia2[z]]->GetLength();
				if (bs != NULL)
					bs->WriteBits(ht.m_oaBitStrings[ia2[z]]);
			}
		}

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("%d Bits written.\n",si-nz);

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("<<< ExportTree Canonical Version <<<\n");

	} else { // Classical

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("ExportTree Normal Version.\n");
		si++;
		if (bs != NULL)
			bs->WriteBit(0);
		b2 = (int)ceil(mylog2((double)m_iaFrequencies.size()+1));
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Using %d bits.\n",b2);
		si+=6;
		if (bs != NULL)
			bs->WriteBits(b2,6);
		if (m_pTree == NULL) {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("m_pTree==NULL.\n");
			si+=b2+1;
			if (bs != NULL) {
				bs->WriteBit(1);
				bs->WriteBits(0,b2);
			}
		} else {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Starting REC_ExportTree() (si=%d).\n",si);
			si += REC_ExportTree(bs,m_pTree,b2);
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Returned from REC_ExportTree() (si=%d).\n",si);
		}
	}

	return si;
}


int CBQBHuffmanTree::REC_ExportTree(CBQBBitSet *bs, CBQBHuffmanSymbol *sym, int bits) {

	int t;

	if (sym->m_oaChildren.size() == 2) {
		if (bs != NULL)
			bs->WriteBit(0);
		t = REC_ExportTree(bs,sym->m_oaChildren[0],bits);
		t += REC_ExportTree(bs,sym->m_oaChildren[1],bits);
	} else {
		t = bits+1;
		if (bs != NULL) {
			bs->WriteBit(1);
			bs->WriteBits(sym->m_iSymbol,bits);
		}
	}

	return t;
}


void CBQBHuffmanTree::ImportTree(CBQBBitSet *bs, bool chr) {

	int z, z2, b2, i, ti, /*ti2,*/ method, border, tia, tib, ts;
	unsigned long l;
	CBQBHuffmanSymbol *hs;
	CBQBHuffmanTree ht(m_IF);
	std::vector<int> mtf;


	if (bs->ReadBit()) { // Canonical

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf(">>> ImportTree Canonical Version >>>\n");

		method = bs->ReadBitsInteger(2);

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Using method %d.\n",method);

		if (bs->ReadBit()) {
			if (chr) {
				if (m_IF.IsPL(BQB_PL_DEBUG))
					m_IF.printf("Reading 9 bits of symbol type count (char, long storage)...\n");
				i = bs->ReadBitsInteger(9);
			} else {
				b2 = bs->ReadBitsInteger(6);
				if (m_IF.IsPL(BQB_PL_DEBUG))
					m_IF.printf("Reading %d bits of symbol type count (long storage)...\n",b2);
				i = bs->ReadBitsInteger(b2);
			}
		} else {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Reading symbol type count (short storage)...\n");
			i = bs->ReadBitsInteger(6);
		}

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Symbol type count is %d.\n",i);

		m_oaSymbols.clear();

		if ((method == 0) || (method == 1)) {

	//		if (chr)
	//			b2 = bs->ReadBitsInteger(4);
	//		else
				b2 = bs->ReadBitsInteger(6);

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Reading %d bits per length token...\n",b2);

			for (z=0;z<i;z++) {
				if (method == 1)
					if (!bs->ReadBit())
						continue;
				ti = bs->ReadBitsInteger(b2);
				if (ti == 0)
					continue;
				hs = new CBQBHuffmanSymbol();
				hs->m_iSymbol = z;
				hs->m_pBitString = new CBQBBitSet(m_IF);
				hs->m_pBitString->WriteBits(0,ti);
				m_oaSymbols.push_back(hs);
			}

		} else if (method == 2) {

	//		if (chr)
	//			b2 = bs->ReadBitsInteger(4);
	//		else
				b2 = bs->ReadBitsInteger(6);

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Reading %d bits per length token...\n",b2);

			if (bs->ReadBit()) {
				if (chr)
					border = bs->ReadBitsInteger(3);
				else
					border = bs->ReadBitsInteger(5);
			} else
				border = 0;

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Border = %d\n",border);

			mtf.resize(pow2i(b2));
			for (z=0;z<(int)mtf.size();z++)
				mtf[z] = z;

			tia = 1;
			tib = 0;
			ts = 0;
			for (z=0;z<i;z++) {

				if (bs->ReadBit()) {
					if (border != 0) {
						if (bs->ReadBit())
							ti = bs->ReadBitsInteger(b2);
						else
							ti = bs->ReadBitsInteger(border);
					} else
						ti = bs->ReadBitsInteger(b2);
				} else
					ti = 0;

		//		m_IF.printf("K %3d: %3d\n",z,ti);

				if (ti == 0) {
					tib += tia;
					tia *= 2;
				} else if (ti == 1) {
					tib += 2*tia;
					tia *= 2;
				} else {
					ti--;
					if (mtf[0] != 0) {
						for (z2=0;z2<tib;z2++) {
							hs = new CBQBHuffmanSymbol();
		//					m_IF.printf("%3lu: %6d (%3d)\n",m_oaSymbols.size(),ts,mtf[0]);
							hs->m_iSymbol = ts++;
							hs->m_pBitString = new CBQBBitSet(m_IF);
							hs->m_pBitString->WriteBits(0,mtf[0]);
							m_oaSymbols.push_back(hs);
						}
					} else
						ts += tib;

					tia = 1;
					tib = 0;

			//		ti2 = mtf[ti];

					MoveToFrontIndex(mtf,ti);
					//std::iter_swap(mtf.begin(),mtf.begin()+ti);

			//		mtf.insert(mtf.begin(),ti2);
			//		mtf.erase(mtf.begin()+ti+1);

					if (mtf[0] != 0) {
						hs = new CBQBHuffmanSymbol();
		//				m_IF.printf("%3lu: %6d (%3d)\n",m_oaSymbols.size(),ts,mtf[0]);
						hs->m_iSymbol = ts++;
						hs->m_pBitString = new CBQBBitSet(m_IF);
						hs->m_pBitString->WriteBits(0,mtf[0]);
						m_oaSymbols.push_back(hs);
					} else
						ts++;
				}

		/*		ti2 = mtf[ti];
				mtf.insert(mtf.begin(),ti2);
				mtf.erase(mtf.begin()+ti+1);

				hs = new CBQBHuffmanSymbol();
				hs->m_iSymbol = ts++;
				hs->m_pBitString = new CBitSet();
				hs->m_pBitString->WriteBits(0,mtf[0]);
				m_oaSymbols.push_back(hs);*/

			}

			if (mtf[0] != 0) {
				for (z2=0;z2<tib;z2++) {
					hs = new CBQBHuffmanSymbol();
		//			m_IF.printf("%3lu: %6d (%3d)\n",m_oaSymbols.size(),ts,mtf[0]);
					hs->m_iSymbol = ts++;
					hs->m_pBitString = new CBQBBitSet(m_IF);
					hs->m_pBitString->WriteBits(0,mtf[0]);
					m_oaSymbols.push_back(hs);
				}
			}

		} else { // Method 3

			ht.ImportTree(bs,chr);

			z2 = 0;
			for (z=0;z<(int)ht.m_oaSymbols.size();z++) {
		//		m_IF.printf("I %3d: %6d\n",z,ht.m_oaSymbols[z]->m_iSymbol);
				if (ht.m_oaSymbols[z]->m_iSymbol > z2)
					z2 = ht.m_oaSymbols[z]->m_iSymbol;
			}

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Found largest symbol %d in tree.\n",z2);

			mtf.resize(z2+1);
			for (z=0;z<(int)mtf.size();z++)
				mtf[z] = z;

			tia = 1;
			tib = 0;
			ts = 0;
			for (z=0;z<i;z++) {

				ti = ht.DecodeSymbol(bs);

				if (ti == 0) {
					tib += tia;
					tia *= 2;
				} else if (ti == 1) {
					tib += 2*tia;
					tia *= 2;
				} else {

					ti--;
					if (mtf[0] != 0) {
						for (z2=0;z2<tib;z2++) {
							hs = new CBQBHuffmanSymbol();
							hs->m_iSymbol = ts++;
							hs->m_pBitString = new CBQBBitSet(m_IF);
							hs->m_pBitString->WriteBits(0,mtf[0]);
							m_oaSymbols.push_back(hs);
						}
					} else
						ts += tib;

					tia = 1;
					tib = 0;

//					ti2 = mtf[ti];

					MoveToFrontIndex(mtf,ti);
					//std::iter_swap(mtf.begin(),mtf.begin()+ti);

			//		mtf.insert(mtf.begin(),ti2);
			//		mtf.erase(mtf.begin()+ti+1);

					if (mtf[0] != 0) {
						hs = new CBQBHuffmanSymbol();
						hs->m_iSymbol = ts++;
						hs->m_pBitString = new CBQBBitSet(m_IF);
						hs->m_pBitString->WriteBits(0,mtf[0]);
						m_oaSymbols.push_back(hs);
					} else
						ts++;
				}

			/*	ti2 = mtf[ti];
				mtf.insert(mtf.begin(),ti2);
				mtf.erase(mtf.begin()+ti+1);


				hs = new CBQBHuffmanSymbol();
				hs->m_iSymbol = ts++;
				hs->m_pBitString = new CBitSet();
				hs->m_pBitString->WriteBits(0,mtf[0]);
				m_oaSymbols.push_back(hs);*/
			}

			if (mtf[0] != 0) {
				for (z2=0;z2<tib;z2++) {
					hs = new CBQBHuffmanSymbol();
					hs->m_iSymbol = ts++;
					hs->m_pBitString = new CBQBBitSet(m_IF);
					hs->m_pBitString->WriteBits(0,mtf[0]);
					m_oaSymbols.push_back(hs);
				}
			}
		}

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Sorting symbols in canonical order...\n");

		std::sort(m_oaSymbols.begin(),m_oaSymbols.end(),SORT_SymbolsCanonical);

	//	m_IF.printf("Have %lu symbols:\n",(unsigned long)m_oaSymbols.size());
	//	for (z=0;z<(int)m_oaSymbols.size();z++)
	//		m_IF.printf("@%d: %d (Len %d)\n",z+1,m_oaSymbols[z]->m_iSymbol,m_oaSymbols[z]->m_pBitString->GetLength());

//		if (m_oaSymbols.size() < 50)
//			for (z=0;z<(int)m_oaSymbols.size();z++)
//				m_IF.printf("A %3d: %6d (%3d)\n",z,m_oaSymbols[z]->m_iSymbol,m_oaSymbols[z]->m_pBitString->GetLength());

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Re-creating bit strings and tree...\n");

		ti = 0;
		l = 0;
		if (m_pTree != NULL)
			delete m_pTree;
		m_pTree = new CBQBHuffmanSymbol();
		for (z=0;z<(int)m_oaSymbols.size();z++) {
			if (m_oaSymbols[z]->m_pBitString->GetLength() > ti) {
				l <<= m_oaSymbols[z]->m_pBitString->GetLength() - ti;
				ti = m_oaSymbols[z]->m_pBitString->GetLength();
			}
	//		if (m_oaSymbols.size() < 50)
	//			for (z2=0;z2<(int)m_oaSymbols.size();z2++)
	//				m_IF.printf("#%d# %3d: %6d\n",z,z2,m_oaSymbols[z2]->m_iSymbol);
			m_oaSymbols[z]->m_pBitString->Clear();
	//		m_IF.printf("%3d (%2d) >>> ",z,ti);
			REC_PushCanonicalSymbol(m_pTree,m_oaSymbols[z],reverse_bit_order(l,ti),ti);
	//		m_IF.printf("\n");
			l++;
		}

//		if (m_oaSymbols.size() < 50)
//			for (z=0;z<(int)m_oaSymbols.size();z++)
//				m_IF.printf("B %3d: %6d\n",z,m_oaSymbols[z]->m_iSymbol);

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("<<< ImportTree Canonical Version <<<\n");

	} else { // Classical

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("ImportTree Normal Version.\n");
		b2 = bs->ReadBitsInteger(6);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Using %d bits.\n",b2);
		m_pTree = new CBQBHuffmanSymbol();
		z = REC_ImportTree(bs,m_pTree,b2);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("REC_ImportTree(): Imported %d symbols.\n",z);
	}
}


void CBQBHuffmanTree::REC_PushCanonicalSymbol(CBQBHuffmanSymbol *sym, CBQBHuffmanSymbol *ts, unsigned long b, int depth) {

//	m_IF.printf("%d",b&1);
	if (depth > 1) {
		if (sym->m_oaChildren.size() == 0) {
			sym->m_oaChildren.push_back(new CBQBHuffmanSymbol());
			sym->m_oaChildren.push_back(new CBQBHuffmanSymbol());
		}
		REC_PushCanonicalSymbol(sym->m_oaChildren[b&1],ts,b>>1,depth-1);
	} else {
		if (sym->m_oaChildren.size() == 0) {
			if ((b&1) == 0) {
				sym->m_oaChildren.push_back(ts);
				sym->m_oaChildren.push_back(new CBQBHuffmanSymbol());
			} else {
				sym->m_oaChildren.push_back(new CBQBHuffmanSymbol());
				sym->m_oaChildren.push_back(ts);
			}
		} else {
			delete sym->m_oaChildren[b&1];
			sym->m_oaChildren[b&1] = ts;
		}
	}
}


int CBQBHuffmanTree::REC_ImportTree(CBQBBitSet *bs, CBQBHuffmanSymbol *sym, int bits) {

	int i;

	i = 0;
	if (!bs->ReadBit()) {
		sym->m_oaChildren.push_back(new CBQBHuffmanSymbol());
		sym->m_oaChildren.push_back(new CBQBHuffmanSymbol());
		i += REC_ImportTree(bs,sym->m_oaChildren[0],bits);
		i += REC_ImportTree(bs,sym->m_oaChildren[1],bits);
	} else {
		sym->m_iSymbol = bs->ReadBitsInteger(bits);
		i++;
	}
	return i;
}

