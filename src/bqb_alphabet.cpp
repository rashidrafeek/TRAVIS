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
#include <algorithm>
#include <set>
#include <math.h>
#include "bqb_alphabet.h"
#include "bqb_hufftree.h"
#include "bqb_interface.h"


const char *GetRevisionInfo_bqb_alphabet(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_alphabet() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



bool SORT_Alphabet_Freq(const CBQBAlphabetEntry *e1, const CBQBAlphabetEntry *e2) {

	return (e1->m_iFrequency > e2->m_iFrequency);
}



bool SORT_Alphabet_Symbol(const CBQBAlphabetEntry *e1, const CBQBAlphabetEntry *e2) {

	return (e1->m_iSymbol < e2->m_iSymbol);
}



void CBQBAlphabet::BuildAlphabet(const std::vector<int> &inp, bool reserved) {

	int z, z2;
	CBQBAlphabetEntry *ae;
	std::vector<std::vector<CBQBAlphabetEntry*> > ta;
	int ti2;
	const unsigned int hash_size = 1024;


	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("Building alphabet (%lu symbols)...\n",(unsigned long)inp.size());
		m_IF.printf("    [");
	}

	for (z=0;z<(int)m_oaAlphabet.size();z++)
		delete m_oaAlphabet[z];
	m_oaAlphabet.clear();
	m_oaAlphabetSortFreq.clear();

	ta.resize(hash_size);

	if (reserved) {
		for (z=C_MIN_RESERVED;true;z++) {
			ti2 = positive_modulo(z, hash_size);
			ae = new CBQBAlphabetEntry();
			ae->m_iSymbol = z;
			ae->m_iFrequency = 0;
			m_oaAlphabet.push_back(ae);
			ta[ti2].push_back(ae);
			if (z == C_MAX_RESERVED)
				break;
		}
	}

	for (z=0;z<(int)inp.size();z++) {
		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			if (fmod(z,inp.size()/60.0) < 1.0) {
				m_IF.printf("#");
				m_IF.FlushLog();
			}
		}
		ti2 = positive_modulo(inp[z], hash_size);
		for (z2=0;z2<(int)ta[ti2].size();z2++)
			if (ta[ti2][z2]->m_iSymbol == inp[z])
				goto _found2;
		z2 = (int)ta[ti2].size();
		ae = new CBQBAlphabetEntry();
		ae->m_iSymbol = inp[z];
		ae->m_iFrequency = 0;
		m_oaAlphabet.push_back(ae);
		ta[ti2].push_back(ae);
_found2:
		ta[ti2][z2]->m_iFrequency++;
	}
	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("]\n");
		m_IF.printf("Found %lu symbol types.\n",(unsigned long)m_oaAlphabet.size());
		m_IF.printf("Sorting alphabet by symbols...\n");
	}

	m_oaAlphabetSortFreq.assign(m_oaAlphabet.begin(),m_oaAlphabet.end());

	std::sort(m_oaAlphabet.begin(),m_oaAlphabet.end(),SORT_Alphabet_Symbol);

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Sorting alphabet by frequencies...\n");

	std::sort(m_oaAlphabetSortFreq.begin(),m_oaAlphabetSortFreq.end(),SORT_Alphabet_Freq);

	for (z=0;z<(int)m_oaAlphabet.size();z++)
		m_oaAlphabet[z]->m_iIndex = z;

	if (m_IF.IsPL(BQB_PL_DEBUG))  {
		m_IF.printf("Creating alphabet indices...\n");
		m_IF.printf("    [");
	}

//	m_IF.printf("Assigning output alphabet...\n");
	m_iaIndices.resize(inp.size());
	for (z=0;z<(int)inp.size();z++) {
		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			if (fmod(z,inp.size()/60.0) < 1.0) {
				m_IF.printf("#");
				m_IF.FlushLog();
			}
		}

		//ti2 = inp[z] % hash_size;
		ti2 = positive_modulo(inp[z], hash_size);
		for (z2=0;z2<(int)ta[ti2].size();z2++)
			if (ta[ti2][z2]->m_iSymbol == inp[z])
				break;
		m_iaIndices[z] = ta[ti2][z2]->m_iIndex;
	}

	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("]\n");
		m_IF.printf("Finished building alphabet.\n");
	}
}


void CBQBAlphabet::Export(CBQBBitSet *bs, bool huffman, bool chr) const {

	int z, bits, ti, border, bborder, bsize, bold;
	long mis, mas, s, m, cht;
	std::vector<int> ia, ia2;
	CBQBHuffmanTree *ht, *ht2;


	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		if (chr)
			m_IF.printf("Exporting alphabet (char flag)...\n");
		else
			m_IF.printf("Exporting alphabet (no char flag)...\n");
	}

	mis = 0x7FFFFFF0;
	mas = -0x7FFFFFF0;
	for (z=0;z<(int)m_oaAlphabet.size();z++) {
		if (m_oaAlphabet[z]->m_iSymbol >= C_MIN_RESERVED)
			break;
		if (m_oaAlphabet[z]->m_iSymbol > mas)
			mas = m_oaAlphabet[z]->m_iSymbol;
		if (m_oaAlphabet[z]->m_iSymbol < mis)
			mis = m_oaAlphabet[z]->m_iSymbol;
	}

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Will write %d symbol types.\n",z);

	if (!chr) {
		bits = (int)ceil(mylog2(z+1));
		if (bits <= 8) {
			bs->WriteBit(0);
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Using short format (8 bits) for symbol type count.\n");
			bs->WriteBits(z,8);
		} else {
			bs->WriteBit(1);
			bs->WriteBits(bits,6);
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Using %d bits for symbol type count.\n",bits);
			bs->WriteBits(z,bits);
		}
	} else
		bs->WriteBits(z,8);

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Symbols range from %ld to %ld.\n",mis,mas);

	if (huffman) {

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf(">>> Huffman coding of alphabet >>>\n");

		bs->WriteBit(1);

		if (!chr) {
			bits = (int)ceil(mylog2(abs(m_oaAlphabet[0]->m_iSymbol)+1));

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("First alphabet entry is %d, using %d bits.\n",m_oaAlphabet[0]->m_iSymbol,bits);

			bs->WriteBits(bits,6);

			if (m_oaAlphabet[0]->m_iSymbol < 0)
				bs->WriteBit(1);
			else
				bs->WriteBit(0);

			bs->WriteBits(abs(m_oaAlphabet[0]->m_iSymbol),bits);
		} else
			bs->WriteBits(m_oaAlphabet[0]->m_iSymbol,8);

		m = 0;
		for (z=0;z<(int)m_oaAlphabet.size()-1;z++) {
			if (m_oaAlphabet[z+1]->m_iSymbol >= C_MIN_RESERVED)
				break;
			if (m_oaAlphabet[z+1]->m_iSymbol-m_oaAlphabet[z]->m_iSymbol > m)
				m = m_oaAlphabet[z+1]->m_iSymbol-m_oaAlphabet[z]->m_iSymbol;
			ia.push_back(m_oaAlphabet[z+1]->m_iSymbol-m_oaAlphabet[z]->m_iSymbol);
		}

		bits = (int)ceil(mylog2(m+1));

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Largest difference is %ld, using %d bits.\n",m,bits);

		ia2.resize(bits+1);

		for (z=0;z<(int)ia2.size();z++)
			ia2[z] = 0;

		for (z=0;z<(int)ia.size();z++)
			ia2[(int)ceil(mylog2(ia[z]+1))]++;

		if (m_IF.IsPL(BQB_PL_DEBUG))
			for (z=0;z<(int)ia2.size();z++)
				m_IF.printf("      %2d bits: %6d\n",z,ia2[z]);

		bborder = -1;
		bsize = 1000000000;
		for (border=0;border<bits;border++) {
			if (chr)
				ti = 0;
			else
				ti = 10;
			for (z=0;z<(int)ia.size();z++)
				if ((int)ceil(mylog2(ia[z]+1)) <= border)
					ti += border+1;
				else
					ti += bits+1;
			if (ti < bsize) {
				bsize = ti;
				bborder = border;
			}
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("      Border=%2d: %.2f Bytes.\n",border,ti/8.0);
		}
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Best border storage would require %.2f bytes (border=%d).\n",bsize/8.0,bborder);

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Creating Normal Huffman tree...\n");

		ht = new CBQBHuffmanTree(m_IF);
		ht->Init(m+1);
		for (z=0;z<(int)ia.size();z++)
			ht->m_iaFrequencies[ia[z]]++;
		ht->BuildTree(false);
		ti = ht->ExportTree(NULL,chr);
		for (z=0;z<(int)ia.size();z++)
			ti += ht->m_oaBitStrings[ia[z]]->GetLength();

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Creating Canonical Huffman tree...\n");

		ht2 = new CBQBHuffmanTree(m_IF);
		ht2->Init(m+1);
		for (z=0;z<(int)ia.size();z++)
			ht2->m_iaFrequencies[ia[z]]++;
		ht2->BuildTree(true);
		cht = ht2->ExportTree(NULL,chr);
		for (z=0;z<(int)ia.size();z++)
			cht += ht2->m_oaBitStrings[ia[z]]->GetLength();

		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			m_IF.printf("Best border storage would require %8.2f bytes.\n",bsize/8.0);
			m_IF.printf("Huffman alphabet storage requires %8.2f Bytes.\n",ti/8.0);
			m_IF.printf("Canonical Huffman would require   %8.2f Bytes.\n",cht/8.0);
		}

		if ((bsize <= ti) && (bsize <= cht)) {

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Falling back to more efficient border algorithm.\n");

			bs->WriteBit(0);

			if (chr) {
				bs->WriteBits(bits,4);
				bs->WriteBits(bborder,3);
			} else {
				bs->WriteBits(bits,6);
				bs->WriteBits(bborder,5);
			}

			for (z=0;z<(int)ia.size();z++) {
				if ((int)ceil(mylog2(ia[z]+1)) <= bborder) {
					bs->WriteBit(0);
					bs->WriteBits(ia[z],bborder);
				} else {
					bs->WriteBit(1);
					bs->WriteBits(ia[z],bits);
				}
			}

		} else if (ti <= cht) {

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Exporting Normal Huffman tree...\n");

			bs->WriteBit(1);

			bold = bs->GetLength();

			ht->ExportTree(bs,chr);

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("%d bits written to Huffman tree.\n",bs->GetLength()-bold);

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Exporting Normal Huffman data (count %lu)...\n",(unsigned long)ia.size());

			bold = bs->GetLength();

			for (z=0;z<(int)ia.size();z++)
				bs->WriteBits(ht->m_oaBitStrings[ia[z]]);

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("%d bits written to Huffman data.\n",bs->GetLength()-bold);

		} else {

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Exporting Canonical Huffman tree...\n");

			bs->WriteBit(1);

			ht2->ExportTree(bs,chr);

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Exporting Canonical Huffman data...\n");

			for (z=0;z<(int)ia.size();z++)
				bs->WriteBits(ht2->m_oaBitStrings[ia[z]]);
		}

		delete ht;
		delete ht2;

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Done.\n");

	} else { // Not Huffman

		bs->WriteBit(0);

		mis = bqbabs(mis);
		mas = bqbabs(mas);
		if (mis > mas)
			mas = mis;
		bits = (int)ceil(mylog2(mas+1))+1;
		s = ((unsigned long)1)<<(bits-1);

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Requires %d bits per symbol, s=%ld.\n",bits,s);

		bs->WriteBits(bits,6);

		for (z=0;z<(int)m_oaAlphabet.size();z++) {
			if (m_oaAlphabet[z]->m_iSymbol >= C_MIN_RESERVED)
				break;
			bs->WriteBits(m_oaAlphabet[z]->m_iSymbol+s,bits);
		}
	}

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Exporting alphabet done.\n");
}


void CBQBAlphabet::Import(CBQBBitSet *bs, bool chr) {

	int z, ti, bits, border, bold;
	unsigned long i;
	long l;
	long s;
	CBQBAlphabetEntry *ae;
	CBQBHuffmanTree *ht;


	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		if (chr)
			m_IF.printf("Importing alphabet (char flag)...\n");
		else
			m_IF.printf("Importing alphabet (no char flag)...\n");
	}

	if (!chr) {
		if (bs->ReadBit()) {
			bits = bs->ReadBitsInteger(6);
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Symbol type count is stored with %d bits.\n",bits);
			i = bs->ReadBitsInteger(bits);
		} else {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Symbol type count short format (8 bits).\n");
			i = bs->ReadBitsInteger(8);
		}
	} else
		i = bs->ReadBitsInteger(8);

	if (chr && (i == 0)) {
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Warning: Symbol type count was 0, assuming 256.\n");
		i = 256;
	}

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Will read %lu symbol types.\n",i);

	m_oaAlphabet.reserve(i);

	if (bs->ReadBit()) { // Huffman

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Reading Huffman-coded alphabet.\n");

		if (!chr) {
			bits = bs->ReadBitsInteger(6);

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("First alphabet symbol has %d bits.\n",bits);

			if (bs->ReadBit())
				l = -((long)bs->ReadBitsInteger(bits));
			else
				l = bs->ReadBitsInteger(bits);
		} else
			l = bs->ReadBitsInteger(8);

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("First alphabet symbol is %ld.\n",l);

		ae = new CBQBAlphabetEntry();
		ae->m_iFrequency = 0;
		ae->m_iIndex = 0;
		ae->m_iSymbol = l;
		m_oaAlphabet.push_back(ae);

		if (bs->ReadBit()) { // Huffman

			ht = new CBQBHuffmanTree(m_IF);

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Reading Huffman table...\n");

			bold = bs->GetReadPos();

			ht->ImportTree(bs,chr);

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Read %d bits from Huffman table.\n",bs->GetReadPos()-bold);

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Tree root element has %lu children.\n",(unsigned long)ht->m_pTree->m_oaChildren.size());

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Reading Huffman data (%ld symbols)...\n",((long)i)-1);

			bold = bs->GetReadPos();

			for (z=0;z<((long)i)-1;z++) {
				ae = new CBQBAlphabetEntry();
				ae->m_iFrequency = 0;
				ae->m_iIndex = (int)m_oaAlphabet.size();
				ae->m_iSymbol = m_oaAlphabet[m_oaAlphabet.size()-1]->m_iSymbol + ht->DecodeSymbol(bs);
				m_oaAlphabet.push_back(ae);
			}

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Read %d bits from Huffman data.\n",bs->GetReadPos()-bold);

			delete ht;

		} else { // Fallback Border Storage

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Fallback border storage.\n");

			if (chr) {
				bits = bs->ReadBitsInteger(4);
				border = bs->ReadBitsInteger(3);
			} else {
				bits = bs->ReadBitsInteger(6);
				border = bs->ReadBitsInteger(5);
			}

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("%d bits maximum, %d bits border.\n",bits,border);

			for (z=0;z<((long)i)-1;z++) {

				if (bs->ReadBit())
					ti = bs->ReadBitsInteger(bits);
				else
					ti = bs->ReadBitsInteger(border);

				ae = new CBQBAlphabetEntry();
				ae->m_iFrequency = 0;
				ae->m_iIndex = (int)m_oaAlphabet.size();
				ae->m_iSymbol = m_oaAlphabet[m_oaAlphabet.size()-1]->m_iSymbol + ti;
				m_oaAlphabet.push_back(ae);
			}

		}

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Done.\n");

	} else { // Not Huffman

		bits = bs->ReadBitsInteger(6);
		s = ((unsigned long)1)<<(bits-1);

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("%lu symbol types, %d bits per symbol type, s=%ld.\n",i,bits,s);

		for (z=0;z<(long)i;z++) {
			ae = new CBQBAlphabetEntry();
			ae->m_iIndex = z;
			m_oaAlphabet.push_back(ae);
			ae->m_iSymbol = bs->ReadBitsInteger(bits) - s;
		}
	}

	for (z=C_MIN_RESERVED;true;z++) {
		ae = new CBQBAlphabetEntry();
		ae->m_iSymbol = z;
		ae->m_iFrequency = 0;
		m_oaAlphabet.push_back(ae);
		if (z == C_MAX_RESERVED)
			break;
	}

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Importing alphabet done.\n");
}


int CBQBAlphabet::FindIndex(int symbol) const {

	int z;

	for (z=0;z<(int)m_oaAlphabet.size();z++)
		if (m_oaAlphabet[z]->m_iSymbol == symbol)
			return z;

	return -1;
}


void CBQBAlphabet::RecalcFrequencies(const std::vector<int> &ia) {

	int z;

	for (z=0;z<(int)m_oaAlphabet.size();z++)
		m_oaAlphabet[z]->m_iFrequency = 0;

	for (z=0;z<(int)ia.size();z++)
		m_oaAlphabet[ia[z]]->m_iFrequency++;

	std::sort(m_oaAlphabetSortFreq.begin(),m_oaAlphabetSortFreq.end(),SORT_Alphabet_Freq);
}


