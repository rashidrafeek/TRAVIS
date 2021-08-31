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
#include "bqb_integerengine.h"
#include "bqb_bitset.h"
#include "bqb_parmset.h"
#include "bqb_hufftree.h"
#include "bqb_alphabet.h"
#include "bqb_interface.h"


const char *GetRevisionInfo_bqb_integerengine(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_integerengine() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}




bool SORT_InverseBW(const CBQBBWPair &p1, const CBQBBWPair &p2) {

	return p1.m_iSymbol < p2.m_iSymbol;
}


bool CBQBIntegerEngine::Compress(
		std::vector<int> &inp,
		CBQBBitSet *outp,
		bool bw,
		bool mtf,
		bool coderun,
		int blocklength,
		int tables,
		bool opttables,
		bool chr,
		int maxiter,
		int maxchunk,
		CBQBStatistics *stat
	) {

	CBQBParameterSet_Compressor parm(m_IF);

	parm.SetBW(bw);
	parm.SetMTF(mtf);
	parm.SetRLE(coderun);
	parm.SetBlockLength(blocklength);
	parm.SetTableCount(tables);
	parm.SetOptTables(opttables);
	parm.SetMaxIter(maxiter);
	parm.SetMaxChunk(maxchunk);

	return Compress(
		inp,
		outp,
		&parm,
		chr,
		stat
	);
}


bool CBQBIntegerEngine::Compress(
		std::vector<int> &inp,
		CBQBBitSet *outp,
		const CBQBParameterSet_Compressor *parm,
		bool chr,
		CBQBStatistics *stat
	) {

	int i, z;
	std::vector<int> tia;
	bool b;
	unsigned long tpos;
	int maxchunk;


	if (parm == NULL) {
		m_IF.eprintf("CBQBIntegerEngine::Compress(): Error: parm == NULL.\n");
		abort();
	}

	maxchunk = parm->GetMaxChunk();

	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("################################################\n");
		m_IF.printf("###   Entering CBQBIntegerEngine::Compress   ###\n");
		m_IF.printf("################################################\n");
		m_IF.printf("Data size %lu, chunk size %d.\n",(unsigned long)inp.size(),maxchunk);
	}

	if (inp.size() == 0) {
		m_IF.eprintf("CBQBIntegerEngine::Compress(): Error: Input array is empty.\n");
		abort();
	}

	tpos = outp->GetLength();

	if (inp.size() > 1000) // Reserves 2 bits per integer in the beginning
		outp->m_iaData.reserve(outp->m_iaData.size()+inp.size()/4);

	// Magic number :-)
	outp->WriteBits(42,8);

	b = true;

	if ((maxchunk > 0) && ((int)inp.size() > maxchunk)) {

		i = (int)ceil((double)inp.size()/maxchunk);
		if (i > 256) {
			m_IF.eprintf("CBQBIntegerEngine::Compress(): More than 256 chunks not supported. Increase chunk size. (%d * %d = %lu)\n",i,maxchunk,(unsigned long)inp.size());
			return false;
		}
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Will write %d chunks.\n",i);
		outp->WriteBit(1);
		outp->WriteBits(i-1,8);

		i = 0;
		z = 0;
		while (i < (int)inp.size()) {

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("### Chunk %d ###\n",z+1);

			if (i+maxchunk >= (int)inp.size())
				tia.assign(inp.begin()+i,inp.end());
			else
				tia.assign(inp.begin()+i,inp.begin()+i+maxchunk);

			if (stat != NULL)
				stat->m_oStat.m_lOverhead += outp->GetLength()-tpos;
			tpos = outp->GetLength();

			if (!CompressSingle(
					tia,
					outp,
					parm,
					chr,
					stat
			)) {
				m_IF.eprintf("CBQBIntegerEngine::Compress(): CompressSingle returned an error in chunk %d.\n",z+1);
				b = false;
				goto _end;
			}

			i += maxchunk;
			z++;
		}

	} else {

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Will write one chunk.\n");

		outp->WriteBit(0); // Only one chunk

		if (stat != NULL)
			stat->m_oStat.m_lOverhead += outp->GetLength()-tpos;
		//tpos = outp->GetLength();

		b = CompressSingle(
				inp,
				outp,
				parm,
				chr,
				stat
		);

	}

_end:
	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("###############################################\n");
		m_IF.printf("###   Leaving CBQBIntegerEngine::Compress   ###\n");
		m_IF.printf("###############################################\n");
	}
	return b;
}


bool CBQBIntegerEngine::Decompress(
		CBQBBitSet *inp,
		std::vector<int> &outp,
		CBQBParameterSet_Compressor *parm
	) {

	unsigned char uc;
	int i, z;
	bool b;


	b = true;

	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("##################################################\n");
		m_IF.printf("###   Entering CBQBIntegerEngine::Decompress   ###\n");
		m_IF.printf("##################################################\n");
	}

	uc = (unsigned char)inp->ReadBitsInteger(8);

	if (uc != 42) {
		m_IF.eprintf("CBQBIntegerEngine::Decompress(): Error in block begin marker: %d != 42.\n",uc);
		b = false;
		goto _end;
	}

	if (inp->ReadBit()) {

		i = inp->ReadBitsInteger(8)+1;

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Will read %d chunks.\n",i);

		for (z=0;z<i;z++) {

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("### Chunk %d ###\n",z+1);

			if (!DecompressSingle(
					inp,
					outp,
					parm
				)) {
				m_IF.eprintf("CBQBIntegerEngine::Decompress(): DecompressSingle returned an error in chunk %d.\n",z+1);
				b = false;
				goto _end;
			}
		}

	} else {

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Will read one chunk.\n");

		b = DecompressSingle(
				inp,
				outp,
				parm
			);
	}

_end:
	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("#################################################\n");
		m_IF.printf("###   Leaving CBQBIntegerEngine::Decompress   ###\n");
		m_IF.printf("#################################################\n");
	}
	return b;
}


bool CBQBIntegerEngine::CompressSingle(
		std::vector<int> &inp,
		CBQBBitSet *outp,
		bool bw,
		bool mtf,
		bool coderun,
		int blocklength,
		int tables,
		bool opttables,
		bool chr,
		int maxiter,
		CBQBStatistics *stat
	) {

	CBQBParameterSet_Compressor parm(m_IF);

	parm.SetBW(bw);
	parm.SetMTF(mtf);
	parm.SetRLE(coderun);
	parm.SetBlockLength(blocklength);
	parm.SetTableCount(tables);
	parm.SetOptTables(opttables);
	parm.SetMaxIter(maxiter);

	return CompressSingle( inp, outp, &parm, chr, stat );
}


bool CBQBIntegerEngine::CompressSingle(
		std::vector<int> &inp,
		CBQBBitSet *outp,
		const CBQBParameterSet_Compressor *parm,
		bool chr,
		CBQBStatistics *stat
	) {

	CBQBAlphabet *alp;
	std::vector<int> iasub, gu, asi, iasi, galp, tia, tia2, tiatrans, tfra;
	std::vector<int> /*bestgu,*/ bestasi, bestiasi;
	std::vector<std::vector<int> > /*lastfreq,*/ bestfreq;
	int ibestfreq;
	int z, z2, i, k, ti, ti2, ti3, tio, nr, iter, beg, end, sw;
	int runtc, itables, stables, ttables, utables, ztables;
	int iruna, irunb, bits, bwindex;
	unsigned long tpos;
	double tf, tf2, tfx, fbestfreq;
	bool b, change;
	CBQBHuffmanTree *ht;
	std::vector<CBQBHuffmanTree*> hta;


	if (parm == NULL) {
		m_IF.eprintf("CBQBIntegerEngine::CompressSingle(): Error: parm == NULL.\n");
		abort();
	}

	bool bw = parm->GetBW();
	bool mtf = parm->GetMTF();
	bool coderun = parm->GetRLE();
	int blocklength = parm->GetBlockLength();
	int tables = parm->GetTableCount();
	bool opttables = parm->GetOptTables();
	int maxiter = parm->GetMaxIter();

	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("\n");
		m_IF.printf(">>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		m_IF.printf(">>>   CompressSingle   >>>\n");
		m_IF.printf(">>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		m_IF.printf("\n");
	}

	if (inp.size() == 0) {
		m_IF.eprintf("CBQBIntegerEngine::CompressSingle(): Error: Input array is empty.\n");
		abort();
	}

	tpos = outp->GetLength();

	// Magic number :-)
	outp->WriteBits(23,6);

	if (chr) {
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Char flag was set.\n");
		outp->WriteBit(1);
	} else {
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("No char flag set.\n");
		outp->WriteBit(0);
	}

/*	if (prerle) {
		if (verbose)
			m_IF.printf("Performing Pre-RLE transformation...\n");
		outp->WriteBit(1);
		ti = 0;
		ti2 = 0;
		tiatrans.clear();
		for (z=0;z<(int)inp.size();z++) {
			if ((inp[z] != ti) || (ti2 >= 255)) {
				if (ti2 >= 4) {
					for (z2=0;z2<4;z2++)
						tiatrans.push_back(ti);
					tiatrans.push_back(ti2-4);
				} else {
					for (z2=0;z2<ti2;z2++)
						tiatrans.push_back(ti);
				}
				ti = inp[z];
				ti2 = 1;
			} else
				ti2++;
		}
		if (ti2 >= 4) {
			for (z2=0;z2<4;z2++)
				tiatrans.push_back(ti);
			tiatrans.push_back(ti2-4);
		} else {
			for (z2=0;z2<ti2;z2++)
				tiatrans.push_back(ti);
		}
		if (verbose)
			m_IF.printf("Done. Reduced %lu Bytes to %lu Bytes.\n",(unsigned long)inp.size(),(unsigned long)tiatrans.size());
	} else {
		if (verbose)
			m_IF.printf("Pre-RLE disabled.\n");
		outp->WriteBit(0);
		tiatrans.assign(inp.begin(),inp.end());
	}*/

	tiatrans.assign(inp.begin(),inp.end());

	alp = new CBQBAlphabet(m_IF);
	alp->BuildAlphabet(tiatrans);

	if (bw) {
		if (IsAllIdentical(alp->m_iaIndices)) {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("!!! Array is \"all identical\" - skipping BW step !!!\n");
			bw = false;
		}
	}

	if (bw) {
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("*** Burrows-Wheeler Transformation ***\n");
		outp->WriteBit(1);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Preparing index buffer...\n");
		//g_iaBW = &alp->m_iaIndices;
		m_iaBW = &alp->m_iaIndices;
		tia.resize(alp->m_iaIndices.size());
		for (z=0;z<(int)tia.size();z++)
			tia[z] = z;
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Building run table...\n");
		//BuildRunTable(*g_iaBW,g_iaBWRunTable);
		BuildRunTable(*m_iaBW,m_iaBWRunTable);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Sorting...\n");

		std::sort(
			tia.begin(),
			tia.end(),
			SSortBWRuntable(*this)
		);

		//std::sort(tia.begin(),tia.end(),SORT_BW_Runtable);
//		std::sort(tia.begin(),tia.end(),SORT_BW);
//		std::stable_sort(tia.begin(),tia.end(),SORT_BW);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Sorting done.\n");
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("BW Step...\n");
		tiatrans.resize(tia.size());
		bwindex = -1;
		for (z=0;z<(int)tia.size();z++) {
			if (tia[z] == 0) {
				bwindex = z;
				tiatrans[z] = alp->m_iaIndices[tia.size()-1];
			} else
				tiatrans[z] = alp->m_iaIndices[tia[z]-1];
		}
/*		FILE *tfi;
		unsigned char tuc;
		tfi = fopen("E:\\bwout.txt","wb");
		for (z=0;z<(int)tia.size();z++) {
			tuc = alp->m_oaAlphabet[tiatrans[z]]->m_iSymbol;
			(void)!fwrite(&tuc,1,1,tfi);
		}
		fclose(tfi);*/
		if (bwindex < 256) {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("BW Index is %d, using short storage.\n",bwindex);
			outp->WriteBit(0);
			outp->WriteBits(bwindex,8);
		} else {
			bits = (int)ceil(mylog2(bwindex+1));
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("BW Index is %d, storing with %d bits.\n",bwindex,bits);
			outp->WriteBit(1);
			outp->WriteBits(bits,6);
			outp->WriteBits(bwindex,bits);
		}
		tia.clear();
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("*** Burrows-Wheeler Finished ***\n");
	} else {
		outp->WriteBit(0);
		tiatrans.assign(alp->m_iaIndices.begin(),alp->m_iaIndices.end());
	}

	if (mtf) {
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Move-To-Front Transformation...\n");
		outp->WriteBit(1);
		galp.resize(alp->m_oaAlphabet.size());
		tia.resize(tiatrans.size());
		for (z=0;z<(int)galp.size();z++)
			galp[z] = z;
		for (z=0;z<(int)tiatrans.size();z++) {
			for (z2=0;z2<(int)galp.size();z2++)
				if (tiatrans[z] == galp[z2])
					break;
	//		if ((z >= 4040534-10) && (z <= 4040534+10))
	//			m_IF.printf("%9d: %4d --> %9d\n",z,tiatrans[z],z2);
			tia[z] = z2;
			if (z2 != 0)
				MoveToFrontIndex(galp,z2);
/*			if (tiatrans[z] == galp[0])
				tia[z] = 0;
			else
				tia[z] = MoveToFront(galp,tiatrans[z]);*/
/*			for (z2=0;z2<(int)galp.size();z2++)
				if (tiatrans[z] == galp[z2])
					break;
			tia[z] = z2;
			if (z2 != 0)
				std::iter_swap(galp.begin(),galp.begin()+z2);*/
		}
		tiatrans.assign(tia.begin(),tia.end());
		tia.clear();
	} else
		outp->WriteBit(0);

	if (coderun) {
		outp->WriteBit(1);
		iruna = alp->FindIndex(C_RUNA);
		irunb = alp->FindIndex(C_RUNB);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Encoding code runs (RUNA=%d, RUNB=%d)...\n",iruna,irunb);
		runtc = 0;
		ti = (int)alp->m_oaAlphabet.size()+100;
		for (i=0;i<(int)tiatrans.size();i++) {

			if (i != 0) {
				if (tiatrans[i] != ti) {
			//		if (runtc != 0) {
						//if (/*(runtc+1 >= 5) && */(ti == 0))
						if (ti == 0)
							PushNumeration(tia,runtc+1,iruna,irunb);
						else
							for (z=0;z<runtc+1;z++)
								tia.push_back(ti);
						runtc = 0;
			//		} else
			//			tia.push_back(ti);
				} else
					runtc++;
			}
			ti = tiatrans[i];
		}

	//	if (runtc != 0) {
			//if (/*(runtc+1 >= 5) && */(ti == 0))
			if (ti == 0)
				PushNumeration(tia,runtc+1,iruna,irunb);
			else
				for (z=0;z<runtc+1;z++)
					tia.push_back(ti);
	//	} else
	//		tia.push_back(ti);

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Recalculating alphabet frequencies...\n");

		alp->RecalcFrequencies(tia);

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Done. Compressed %lu to %lu symbols.\n",(unsigned long)tiatrans.size(),(unsigned long)tia.size());

	//	for (z=0;z<(int)tia.size();z++)
	//		m_IF.printf("@ %2d: %d\n",z+1,tia[z]);

		tiatrans.assign(tia.begin(),tia.end());
		tia.clear();
	} else 
		outp->WriteBit(0);

	if (stat != NULL)
		stat->m_oStat.m_lOverhead += outp->GetLength()-tpos;
	tpos = outp->GetLength();

	tia.assign(tiatrans.begin(),tiatrans.end());

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Writing alphabet...\n");
	i = outp->GetLength();
	alp->Export(outp,true,chr);
	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("%d Bytes written.\n",(outp->GetLength()-i)/8);

	if (stat != NULL)
		stat->m_oStat.m_lAlphabet += outp->GetLength()-tpos;
	tpos = outp->GetLength();

	if (tables == 1) {

_onetab:
		outp->WriteBit(0);
		outp->WriteBits(0,3);

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Creating Huffman tree...\n");

		ht = new CBQBHuffmanTree(m_IF);
		ht->Init((int)alp->m_oaAlphabet.size());

		for (z=0;z<(int)tia.size();z++)
			ht->m_iaFrequencies[tia[z]]++;

		ht->BuildTree(true);

		if (stat != NULL)
			stat->m_oStat.m_lOverhead += outp->GetLength()-tpos;
		tpos = outp->GetLength();

		i = outp->GetLength();

		ht->ExportTree(outp,chr);

		if (stat != NULL)
			stat->m_oStat.m_lHuffmanTables += outp->GetLength()-tpos;
		tpos = outp->GetLength();

		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			m_IF.printf("%d Bytes (%d Bits) written.\n",(outp->GetLength()-i)/8,outp->GetLength()-i);
			m_IF.printf("Writing Huffman-encoded data...\n");
		}

		bits = (int)ceil(mylog2((double)tia.size()+1));

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("%lu symbols to be written, using %d bits for symbol count.\n",(unsigned long)tia.size(),bits);

		outp->WriteBits(bits,6);
		outp->WriteBits((unsigned long)tia.size(),bits);

		if (stat != NULL)
			stat->m_oStat.m_lOverhead += outp->GetLength()-tpos;
		tpos = outp->GetLength();

		k = outp->GetLength();

		for (z=0;z<(int)tia.size();z++)
			outp->WriteBits(ht->m_oaBitStrings[tia[z]]);

		if (m_IF.IsPL(BQB_PL_DEBUG)) 
			m_IF.printf("%d Bytes (%d Bits) written.\n",(outp->GetLength()-k)/8,outp->GetLength()-k);

		if (stat != NULL)
			stat->m_oStat.m_lHuffmanData += outp->GetLength()-tpos;
		//tpos = outp->GetLength();

		delete ht;

	} else {

		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/*		if (preopt) { // Alternative multi-table optimization
//		if (true) { // Alternative multi-table optimization

			MultiHuffmanOptimize(tables,hta,alp,tia,asi,iasi,blocklength,verbose);

		//	itables = tables;
			itables = (int)hta.size();

			for (z=0;z<itables;z++)
				hta[z]->BuildPrelimTree();

		} else { // Classical multi-table optimization*/

			if (opttables) {

				stables = 1000000000;
				itables = 1;
				for (ztables=1;ztables<=tables;ztables++) {

					ttables = ztables;
					utables = 0;

					if (ttables == 1) {

						utables += 4;

						if (m_IF.IsPL(BQB_PL_DEBUG))
							m_IF.printf("$1 Creating Huffman tree...\n");

						ht = new CBQBHuffmanTree(m_IF);
						ht->Init((int)alp->m_oaAlphabet.size());

						for (z=0;z<(int)tia.size();z++)
							ht->m_iaFrequencies[tia[z]]++;

						ht->BuildTree(true);

						utables += ht->ExportTree(NULL,chr);

						if (m_IF.IsPL(BQB_PL_DEBUG))
							m_IF.printf("$1 Writing Huffman-encoded data...\n");

						bits = (int)ceil(mylog2((double)tia.size()+1));

						if (m_IF.IsPL(BQB_PL_DEBUG))
							m_IF.printf("$1 %lu symbols to be written, using %d bits for symbol count.\n",(unsigned long)tia.size(),bits);

						utables += 5+bits;

						for (z=0;z<(int)tia.size();z++)
							utables += ht->m_oaBitStrings[tia[z]]->GetLength();

						delete ht;

					} else {

						ti = 0;
						tio = 0;
						nr = (int)tia.size();
						hta.clear();
						if (m_IF.IsPL(BQB_PL_DEBUG))
							m_IF.printf("$%d Initializing table populations...\n",ttables);
						for (i=0;i<ttables;i++) {
							ht = new CBQBHuffmanTree(m_IF);
							ht->Init((int)alp->m_oaAlphabetSortFreq.size());
							hta.push_back(ht);

							ti2 = 0;
							while ((ti2 < nr/(ttables-i)) && (ti < (int)alp->m_oaAlphabetSortFreq.size()))
								ti2 += alp->m_oaAlphabetSortFreq[ti++]->m_iFrequency;
							if ((i != 0) && (i != ttables) && ((i%2) == 1) && (ti < (int)alp->m_oaAlphabetSortFreq.size()) && (ti > tio+1)) {
								ti--;
								ti2 -= alp->m_oaAlphabetSortFreq[ti]->m_iFrequency;
							}
							if (m_IF.IsPL(BQB_PL_DEBUG))
								m_IF.printf("$%d     Group %2d: Symbols %5d .. %5d / %lu, Count %8d (%7.3f%%)\n",ttables,i+1,tio+1,ti,(unsigned long)alp->m_oaAlphabetSortFreq.size(),ti2,(double)ti2/alp->m_iaIndices.size()*100.0);

							for (z=0;z<(int)alp->m_oaAlphabetSortFreq.size();z++)
								if ((z >= tio) && (z < ti))
									ht->m_iaLengths[alp->m_oaAlphabetSortFreq[z]->m_iIndex] = 0;
								else
									ht->m_iaLengths[alp->m_oaAlphabetSortFreq[z]->m_iIndex] = 20;

							tio = ti;
							nr -= ti2;
						}

						if (m_IF.IsPL(BQB_PL_DEBUG))
							m_IF.printf("$%d Starting optimization...\n",ttables);
						iasub.resize(blocklength);
						asi.resize(tia.size()/blocklength+1);
						iasi.resize(tia.size()/blocklength+1);
						gu.resize(ttables);
						galp.resize(ttables);

						tf2 = -3456.0;
						iter = 0;
						while (true) {
							i = 0;
							for (z=0;z<ttables;z++)
								gu[z] = 0;
							for (k=0;k<ttables;k++)
								for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
									hta[k]->m_iaFrequencies[z] = 0;
							tf = 0;
							tfx = 0;
							sw = 0;
							beg = 0;
							for (z=0;z<ttables;z++)
								galp[z] = z;
							while (true) {
								end = beg+blocklength-1;
								if (end >= (int)tia.size())
									end = ((int)tia.size())-1;

								ti = 5000000;
								ti2 = -1;
								ti3 = -1;
								for (k=0;k<ttables;k++) {
									tio = 0;
									for (z=beg;z<=end;z++)
										tio += hta[k]->m_iaLengths[tia[z]];
									for (z=0;z<ttables;z++)
										if (galp[z] == k)
											break;
							//		tio += z+1;
									if (ttables > 1)
										tio += z+1;
									if (tio < ti) {
										ti = tio;
										ti2 = k;
										ti3 = z;
									}
								}
								asi[i] = ti2;
								iasi[i] = ti3;
								tfx += ti3+1;
								gu[ti2]++;
								if (ti3 != 0) {
									MoveToFrontIndex(galp,ti3);
									//std::iter_swap(galp.begin(),galp.begin()+ti3);
									sw++;
								}
								for (z=beg;z<=end;z++)
									hta[ti2]->m_iaFrequencies[tia[z]]++;
								tf += ti/8.0;

								beg = end+1;
								i++;
								if (beg == (int)tia.size())
									break;
							}
							if (m_IF.IsPL(BQB_PL_DEBUG))
								m_IF.printf("$%d    Iteration %2d, Total size: %8.3f KiB\n",ttables,iter+1,tf/1024.0);
							for (z=0;z<ttables;z++)
								hta[z]->BuildPrelimTree();
							if (tf2 == tf)
								break;
							tf2 = tf;
							iter++;
							if (iter > 40) {
								if (m_IF.IsPL(BQB_PL_DEBUG))
									m_IF.printf("$%d    Warning: Multi-Huffman optimization did not converge after 40 iterations.\n",ttables);
								break;
							}
						}
						if (m_IF.IsPL(BQB_PL_DEBUG))
							m_IF.printf("$%d Group statistics:\n",ttables);
						tf2 = tfx/8.0;
						b = false;
						for (z=0;z<ttables;z++) {
							hta[z]->BuildTree(true);
							if (m_IF.IsPL(BQB_PL_DEBUG))
								m_IF.printf("$%d    %d: %5.2f%%, %6d blocks, %5lu/%5lu symbols, maxlen %d bits\n",ttables,z+1,(double)gu[z]/tia.size()*blocklength*100.0,gu[z],(unsigned long)hta[z]->m_oaSymbols.size(),(unsigned long)alp->m_oaAlphabet.size(),hta[z]->m_iMaxBitLength);
							if (hta[z]->m_pTree == NULL) {
								if (m_IF.IsPL(BQB_PL_DEBUG))
									m_IF.printf("$%d    Warning: Empty Huffman tree (%d/%d).\n",ttables,z+1,ttables);
								b = true;
							}
							tf2 += tf;
						}

						if (m_IF.IsPL(BQB_PL_DEBUG))
							m_IF.printf("$%d Exporting Huffman trees...\n",ttables);

						if (b) {
							if (m_IF.IsPL(BQB_PL_DEBUG))
								m_IF.printf("$%d Removing empty trees and reorganizing structures...\n",ttables);
							tia2.clear();
							z2 = 0;
							for (z=0;z<ttables;z++) {
								tia2.push_back(z2);
								if (hta[z]->m_pTree != NULL)
									z2++;
							}

							ttables = z2;

							if (m_IF.IsPL(BQB_PL_DEBUG))
								m_IF.printf("$%d Now only %d tables left.\n",ttables,ttables);

							for (z=0;z<(int)hta.size();z++) {
								if (hta[z]->m_pTree == NULL) {
									delete hta[z];
									hta.erase(hta.begin()+z);
									z--;
								}
							}

							galp.resize(ttables);

							for (z=0;z<ttables;z++)
								galp[z] = z;

							for (z=0;z<(int)asi.size();z++) {

								asi[z] = tia2[asi[z]];

								for (z2=0;z2<ttables;z2++)
									if (galp[z2] == asi[z])
										break;

								iasi[z] = z2;

								if (z2 != 0)
									MoveToFrontIndex(galp,z2);
									//std::iter_swap(galp.begin(),galp.begin()+z2);
							}
						}

						utables += 12;

						for (z=0;z<ttables;z++) {
							if (m_IF.IsPL(BQB_PL_DEBUG))
								m_IF.printf("$%d *** Tree %d ***\n",ttables,z+1);
							utables += hta[z]->ExportTree(NULL,chr);
						}

						bits = (int)ceil(mylog2((double)tia.size()+1));
						if (m_IF.IsPL(BQB_PL_DEBUG)) {
							m_IF.printf("$%d Writing Huffman-encoded data...\n",ttables);
							m_IF.printf("$%d %lu symbols to be written, using %d bits for symbol count.\n",ttables,(unsigned long)tia.size(),bits);
						}

						utables += 5+bits;

						beg = 0;
						i = 0;
						while (true) {

							end = beg+blocklength-1;
							if (end >= (int)tia.size())
								end = ((int)tia.size())-1;

							utables += iasi[i]+1;

							for (z=beg;z<=end;z++)
								utables += hta[asi[i]]->m_oaBitStrings[tia[z]]->GetLength();

							beg = end+1;
							i++;
							if (beg == (int)tia.size())
								break;
						}

						for (i=0;i<ttables;i++)
							delete hta[i];
					}

					if (utables < stables) {
						stables = utables;
						itables = ztables;
					}
				}

				if (m_IF.IsPL(BQB_PL_DEBUG)) {
					m_IF.printf("\n");
					m_IF.printf("### Optimization result: Using %d tables.\n",itables);
					m_IF.printf("\n");
				}
				else if (m_IF.IsPL(BQB_PL_STANDARD))
					m_IF.printf("      Using %d tables.\n",itables);

				if (itables == 1)
					goto _onetab;

			} else
				itables = tables;

			ti = 0;
			tio = 0;
			nr = (int)tia.size();
			hta.clear();
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Initializing table populations...\n");
			for (i=0;i<itables;i++) {
				ht = new CBQBHuffmanTree(m_IF);
				ht->Init((int)alp->m_oaAlphabetSortFreq.size());
				hta.push_back(ht);

				ti2 = 0;
				while ((ti2 < nr/(itables-i)) && (ti < (int)alp->m_oaAlphabetSortFreq.size()))
					ti2 += alp->m_oaAlphabetSortFreq[ti++]->m_iFrequency;
				if ((i != 0) && (i != itables) && ((i%2) == 1) && (ti < (int)alp->m_oaAlphabetSortFreq.size()) && (ti > tio+1)) {
					ti--;
					ti2 -= alp->m_oaAlphabetSortFreq[ti]->m_iFrequency;
				}
				if (m_IF.IsPL(BQB_PL_DEBUG))
					m_IF.printf("    Group %2d: Symbols %5d .. %5d / %lu, Count %8d (%7.3f%%)\n",i+1,tio+1,ti,(unsigned long)alp->m_oaAlphabetSortFreq.size(),ti2,(double)ti2/tia.size()*100.0);

				for (z=0;z<(int)alp->m_oaAlphabetSortFreq.size();z++)
					if ((z >= tio) && (z < ti))
						ht->m_iaLengths[alp->m_oaAlphabetSortFreq[z]->m_iIndex] = 0;
					else
						ht->m_iaLengths[alp->m_oaAlphabetSortFreq[z]->m_iIndex] = 20;

				tio = ti;
				nr -= ti2;
			}

			for (k=0;k<itables;k++)
				for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
					hta[k]->m_iaFrequencies[z] = 0;

	//	} // End if preopt

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Starting optimization...\n");

			iasub.resize(blocklength);
			asi.resize(tia.size()/blocklength+1);
			iasi.resize(tia.size()/blocklength+1);
			gu.resize(itables);
			galp.resize(itables);

		//	lastfreq.resize(itables);
			bestfreq.resize(itables);

/*			CHuffmanEstimator he;
			std::vector<std::vector<unsigned char> > tmatrix;

			m_IF.printf("Building tmatrix...\n");

			he.Init(alp->m_oaAlphabet.size());

			tmatrix.resize(tia.size()/blocklength+1);

			for (z=0;z<tia.size()/blocklength+1;z++) {
				tmatrix[z].resize(alp->m_oaAlphabet.size());
				for (z2=0;z2<alp->m_oaAlphabet.size();z2++)
					tmatrix[z][z2] = 0;
			}

			for (z=0;z<(int)tia.size();z++)
				tmatrix[z/blocklength][tia[z]]++;

			m_IF.printf("Done.\n");*/

			tf2 = -3456.0;
			iter = 0;
			fbestfreq = 1.0e20;
			ibestfreq = -10;
			while (true) {
				i = 0;
				for (z=0;z<itables;z++)
					gu[z] = 0;

				for (k=0;k<itables;k++)
					for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
						hta[k]->m_iaFrequencies[z] /= 2;
				//		hta[k]->m_iaFrequencies[z] = 0;

				for (z=0;z<itables;z++)
					galp[z] = z;

				tf = 0;
				tfx = 0;
				sw = 0;
				beg = 0;
				change = false;
				while (true) {

					end = beg+blocklength-1;
					if (end >= (int)tia.size())
						end = ((int)tia.size())-1;

					ti = 5000000;
					ti2 = -1;
					ti3 = -1;
					for (k=0;k<itables;k++) {
						tfra.assign(hta[k]->m_iaFrequencies.begin(),hta[k]->m_iaFrequencies.end());
						tio = 0;
						for (z=beg;z<=end;z++) {
							tio += hta[k]->m_iaLengths[tia[z]];
							if (tfra[tia[z]] == 0) {
								tfra[tia[z]]++;
								tio += 1;
					//			tio += 1+(gu[k]*5)/(i+1);
							}
						}
						for (z=0;z<itables;z++)
							if (galp[z] == k)
								break;
						if (itables > 1)
							tio += z+1;
						if (tio < ti) {
							ti = tio;
							ti2 = k;
							ti3 = z;
						}
					}
					if (asi[i] != ti2)
						change = true;
					asi[i] = ti2;
					iasi[i] = ti3;
					tfx += ti3+1;
					gu[ti2]++;
					if (ti3 != 0) {
						MoveToFrontIndex(galp,ti3);
						//std::iter_swap(galp.begin(),galp.begin()+ti3);
						sw++;
					}
					for (z=beg;z<=end;z++)
						hta[ti2]->m_iaFrequencies[tia[z]]++;
					tf += ti/8.0;

					beg = end+1;
					i++;
					if (beg == (int)tia.size())
						break;
				}

/*				for (i=0;i<tia.size()/blocklength+1;i++) {
			//		printf("B%d\n",i);
					ti = 5000000;
					ti2 = -1;
					ti3 = -1;
					for (k=0;k<itables;k++) {
			//			tfra.assign(hta[k]->m_iaFrequencies.begin(),hta[k]->m_iaFrequencies.end());

						tio = he.EstimateBitLength(hta[k]->m_iaFrequencies,tmatrix[i]);

			//			printf("C%d\n",k);

						for (z=0;z<itables;z++)
							if (galp[z] == k)
								break;
						if (itables > 1)
							tio += z+1;
						if (tio < ti) {
							ti = tio;
							ti2 = k;
							ti3 = z;
						}
					}
					asi[i] = ti2;
					iasi[i] = ti3;
					tfx += ti3+1;
					gu[ti2]++;
					if (ti3 != 0) {
						std::iter_swap(galp.begin(),galp.begin()+ti3);
						sw++;
					}
			//		for (z=beg;z<=end;z++)
			//			hta[ti2]->m_iaFrequencies[tia[z]]++;
					for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
						hta[ti2]->m_iaFrequencies[z] += tmatrix[i][z];
					tf += ti/8.0;
				}*/


				if (m_IF.IsPL(BQB_PL_DEBUG))
					m_IF.printf("    Iteration %2d, Total size: %8.3f KiB\n",iter+1,tf/1024.0);

//				m_IF.printf("$ %d\n",hta[1]->m_iaFrequencies[1492]);

				if (tf < fbestfreq) {
					fbestfreq = tf;
					for (z=0;z<itables;z++)
						bestfreq[z].assign(hta[z]->m_iaFrequencies.begin(),hta[z]->m_iaFrequencies.end());
		//				bestfreq[z].assign(lastfreq[z].begin(),lastfreq[z].end());
		//			bestgu.assign(gu.begin(),gu.end());
					bestasi.assign(asi.begin(),asi.end());
					bestiasi.assign(iasi.begin(),iasi.end());
					ibestfreq = iter;
				}
				for (z=0;z<itables;z++)
					hta[z]->BuildPrelimTree();
				if (!change)
					break;
				tf2 = tf;
				iter++;
		//		if (alp->m_oaAlphabet.size() > 500)
		//			break; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				if (iter >= maxiter) {
					if (m_IF.IsPL(BQB_PL_DEBUG))
						m_IF.printf("Multi-Huffman optimization did not converge.\n");
					break;
				}
			}
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Best result in iteration %d.\n",ibestfreq+1);
			if (iter != 1) {
	//			gu.assign(bestgu.begin(),bestgu.end());
				asi.assign(bestasi.begin(),bestasi.end());
				iasi.assign(bestiasi.begin(),bestiasi.end());
				for (z=0;z<itables;z++) {
					hta[z]->m_iaFrequencies.assign(bestfreq[z].begin(),bestfreq[z].end());
					hta[z]->BuildPrelimTree();
				}
			}
//		}


		for (k=0;k<itables;k++)
			for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
				hta[k]->m_iaFrequencies[z] = 0;

		beg = 0;
		i = 0;
		while (true) {
			end = beg+blocklength-1;
			if (end >= (int)tia.size())
				end = ((int)tia.size())-1;
			for (z=beg;z<=end;z++)
				hta[asi[i]]->m_iaFrequencies[tia[z]]++;
			beg = end+1;
			if (beg == (int)tia.size())
				break;
			i++;
		}


//		m_IF.printf("$$ %d\n",hta[1]->m_iaFrequencies[1492]);

		gu.resize(itables);
		for (z=0;z<itables;z++)
			gu[z] = 0;
		for (z=0;z<(int)asi.size();z++)
			gu[asi[z]]++;

//		verbose = true;

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Group statistics:\n");
		tf2 = tfx/8.0;
		b = false;
		for (z=0;z<itables;z++) {
			hta[z]->BuildTree(true);
			tf = 0;
			for (z2=0;z2<(int)alp->m_oaAlphabet.size();z2++)
				tf += hta[z]->m_iaLengths[z2] * hta[z]->m_iaFrequencies[z2];
			if (m_IF.IsPL(BQB_PL_DEBUG) && (hta[z]->m_pTree != NULL))
				m_IF.printf("    %2d: %5.2f%%, %6d blocks, %5lu/%5lu symbols, %6.3f bits/symbol, maxlen %d bits\n",z+1,(double)gu[z]/tia.size()*blocklength*100.0,gu[z],(unsigned long)hta[z]->m_oaSymbols.size(),(unsigned long)alp->m_oaAlphabet.size(),tf/gu[z]/blocklength,hta[z]->m_iMaxBitLength);
			if (hta[z]->m_pTree == NULL) {
		//		if (verbose)
		//			m_IF.printf("Warning: Empty Huffman tree (%d/%d).\n",z+1,itables);
				b = true;
			}
			tf2 += tf;
		}

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Exporting Huffman trees...\n");

		if (b) {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Removing empty trees and reorganizing structures...\n");
			tia2.clear();
			z2 = 0;
			for (z=0;z<itables;z++) {
				tia2.push_back(z2);
				if (hta[z]->m_pTree != NULL)
					z2++;
			}

			itables = z2;

			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Now only %d tables left.\n",itables);

			for (z=0;z<(int)hta.size();z++) {
				if (hta[z]->m_pTree == NULL) {
					delete hta[z];
					hta.erase(hta.begin()+z);
					z--;
				}
			}

			galp.resize(itables);

			for (z=0;z<itables;z++)
				galp[z] = z;

			for (z=0;z<(int)asi.size();z++) {

				asi[z] = tia2[asi[z]];

				for (z2=0;z2<itables;z2++)
					if (galp[z2] == asi[z])
						break;

				iasi[z] = z2;

				if (z2 != 0)
					MoveToFrontIndex(galp,z2);
					//std::iter_swap(galp.begin(),galp.begin()+z2);
			}
		}

		if (itables == 1) {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Only 1 table left, going to dedicated one-table code.\n");
			goto _onetab;
		}

		if (itables <= 8) {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Have %d <= 8 Huffman trees, using short storage.\n",itables);
			outp->WriteBit(0);
			outp->WriteBits(itables-1,3);
		} else {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Have %d > 8 Huffman trees, using long storage.\n",itables);
			if (itables >= 1024) {
				m_IF.eprintf("CBQBIntegerEngine::CompressSingle(): Error: More than 1023 Huffman trees not supported (have %d).\n",itables);
				abort();
			}
			outp->WriteBit(1);
			outp->WriteBits(itables-1,10);
		}
		outp->WriteBits(blocklength,8);

		if (stat != NULL)
			stat->m_oStat.m_lOverhead += outp->GetLength()-tpos;
		tpos = outp->GetLength();

		i = outp->GetLength();
		for (z=0;z<itables;z++) {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("*** Tree %d ***\n",z+1);
			hta[z]->ExportTree(outp,chr);
		}

		if (stat != NULL)
			stat->m_oStat.m_lHuffmanTables += outp->GetLength()-tpos;
		tpos = outp->GetLength();

		bits = (int)ceil(mylog2((double)tia.size()+1));
		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			m_IF.printf("%d Bytes written.\n",(outp->GetLength()-i)/8);
			m_IF.printf("Writing Huffman-encoded data...\n");
			m_IF.printf("%lu symbols to be written, using %d bits for symbol count.\n",(unsigned long)tia.size(),bits);
		}

		outp->WriteBits(bits,6);
		outp->WriteBits((unsigned long)tia.size(),bits);

/*		std::vector<int> blu, blu2;

		blu.resize(alp->m_oaAlphabet.size());
		blu2.resize(alp->m_oaAlphabet.size());

		for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
			blu[z] = 0;*/

		if (stat != NULL)
			stat->m_oStat.m_lOverhead += outp->GetLength()-tpos;
		tpos = outp->GetLength();

		k = outp->GetLength();
		beg = 0;
		i = 0;
		ti = 0;
		while (true) {

//			m_IF.printf("*** Block %d ***\n",i+1);

			end = beg+blocklength-1;
			if (end >= (int)tia.size())
				end = ((int)tia.size())-1;

/*			for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
				blu2[z] = 0;
			for (z=beg;z<=end;z++)
				blu2[tia[z]]++;
			ti = 0;
			for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
				if (blu2[z] != 0)
					ti++;
			blu[ti]++; */

			for (z=0;z<iasi[i];z++)
				outp->WriteBit(1);
			outp->WriteBit(0);
			ti += iasi[i]+1;

			if (stat != NULL)
				stat->m_oStat.m_lTableSwitch += outp->GetLength()-tpos;
			tpos = outp->GetLength();

//			m_IF.printf("A %d\n",asi[i]);

			for (z=beg;z<=end;z++) {
/*				if (i == 5804) {
					m_IF.printf("    %3d / %3lu\n",tia[z],hta[asi[i]]->m_oaBitStrings.size());
					m_IF.printf("    @ %08X %d\n",hta[asi[i]]->m_oaBitStrings[tia[z]],hta[asi[i]]->m_iaFrequencies[tia[z]]);
					m_IF.printf("    @ %d\n",hta[asi[i]]->m_oaBitStrings[tia[z]]->GetLength());
				}*/
				outp->WriteBits(hta[asi[i]]->m_oaBitStrings[tia[z]]);
			}

			if (stat != NULL)
				stat->m_oStat.m_lHuffmanData += outp->GetLength()-tpos;
			tpos = outp->GetLength();

			beg = end+1;
			i++;
			if (beg == (int)tia.size())
				break;
		}

		if (m_IF.IsPL(BQB_PL_DEBUG)) 
			m_IF.printf("%d Bytes written (including %d Bytes Switching penalty).\n",(outp->GetLength()-k)/8,ti/8);

		ti = i;

		for (i=0;i<itables;i++)
			delete hta[i];

/*		if (alp->m_oaAlphabet.size() < 10000) {
			m_IF.printf("Block statistics (%d blocks in total):\n",ti);
			for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
				if (blu[z] != 0)
					m_IF.printf("  %6d/%6lu symbols: %6d (%5.2f%%)\n",z,alp->m_oaAlphabet.size(),blu[z],blu[z]*100.0/ti);
		}*/

/*		if (alp->m_oaAlphabet.size() == 130) {
			FILE *tta;
			tta = fopen("matrix.nb","wt");
			fprintf(tta,"{ \n");

			beg = 0;
			while (true) {

				end = beg+blocklength-1;
				if (end >= (int)tia.size())
					end = ((int)tia.size())-1;

				for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
					blu2[z] = 0;
				for (z=beg;z<=end;z++)
					blu2[tia[z]]++;

				fprintf(tta,"{ ");

				for (z=0;z<(int)alp->m_oaAlphabet.size();z++) {
					fprintf(tta,"%d",blu2[z]);
					if (z+1 < alp->m_oaAlphabet.size())
						fprintf(tta,",");
				}
				fprintf(tta,"}");

				beg = end+1;
				if (beg == (int)tia.size())
					break;
				fprintf(tta,",\n");
			}
			fprintf(tta," }\n");

			fclose(tta);
			abort();
		}*/
	}

	delete alp;

	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("\n");
		m_IF.printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
		m_IF.printf("<<<   CompressSingle Done   <<<\n");
		m_IF.printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
		m_IF.printf("\n");
	}

	return true;
}


bool CBQBIntegerEngine::DecompressSingle(
		CBQBBitSet *inp,
		std::vector<int> &outp,
		CBQBParameterSet_Compressor *parm
	) {

	CBQBAlphabet *alp;
	int i, k, z, z2, c, dc, ti, ti2, bits;
	int blocklength, tables, iruna, irunb, bwindex;
	std::vector<CBQBHuffmanTree*> hta;
	CBQBHuffmanTree *ht;
	std::vector<int> galp, tia, tiatrans;
	std::vector<CBQBBWPair> tbwpairs;
	bool bw, mtf, coderun, chr;


	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("\n");
		m_IF.printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		m_IF.printf(">>>   DecompressSingle   >>>\n");
		m_IF.printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		m_IF.printf("\n");
	}

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Checking block begin magic number (6 bits)...\n");

	i = inp->ReadBitsInteger(6);

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Found %d.\n",i);

	if (i != 23) {
		m_IF.eprintf("CBQBIntegerEngine::Decompress(): Error in block begin: %d != 23.\n",i);
		return false;
	}

	chr = inp->ReadBit();
	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		if (chr)
			m_IF.printf("Found char flag.\n");
		else
			m_IF.printf("No char flag found.\n");
	}

	bwindex = -1;
	if (inp->ReadBit()) {
		bw = true;
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Burrows-Wheeler transform was applied.\n");
		if (inp->ReadBit()) {
			bits = inp->ReadBitsInteger(6);
			bwindex = inp->ReadBitsInteger(bits);
		} else
			bwindex = inp->ReadBitsInteger(8);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Burrows-Wheeler index is %d.\n",bwindex);
	} else {
		bw = false;
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("No Burrows-Wheeler transform.\n");
	}
	if (parm != NULL)
		parm->SetBW(bw);

	if (inp->ReadBit()) {
		mtf = true;
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Move-to-Front was applied.\n");
	} else {
		mtf = false;
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("No Move-to-Front.\n");
	}
	if (parm != NULL)
		parm->SetMTF(mtf);

	if (inp->ReadBit()) {
		coderun = true;
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("CodeRun was applied.\n");
	} else {
		coderun = false;
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("No CodeRun.\n");
	}
	if (parm != NULL)
		parm->SetRLE(coderun);

	alp = new CBQBAlphabet(m_IF);

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Reading alphabet...\n");
	i = inp->GetReadPos();
	alp->Import(inp,chr);
	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("%d Bytes read, %lu symbol types.\n",(inp->GetReadPos()-i)/8,(unsigned long)alp->m_oaAlphabet.size());

	if (inp->ReadBit()) {
		tables = inp->ReadBitsInteger(10)+1;
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Using %d Huffman tables (long storage).\n",tables);
	} else {
		tables = inp->ReadBitsInteger(3)+1;
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Using %d Huffman tables (short storage).\n",tables);
	}
	if (parm != NULL)
		parm->SetTableCount(tables);

	if (tables > 1) {

		blocklength = inp->ReadBitsInteger(8);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Using a block length of %d.\n",blocklength);
		if (parm != NULL)
			parm->SetBlockLength(blocklength);

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Importing Huffman trees...\n");
		i = inp->GetReadPos();
		hta.resize(tables);
		for (z=0;z<tables;z++) {
			if (m_IF.IsPL(BQB_PL_DEBUG))
				m_IF.printf("Tree %d...\n",z+1);
			hta[z] = new CBQBHuffmanTree(m_IF);
			hta[z]->ImportTree(inp,chr);
		}
		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			m_IF.printf("%d Bytes read.\n",(inp->GetReadPos()-i)/8);
			m_IF.printf("Importing and decoding data stream...\n");
		}
		bits = inp->ReadBitsInteger(6);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Symbol count stored with %d bits.\n",bits);
		dc = inp->ReadBitsInteger(bits);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Expecting %d symbols.\n",dc);
		k = inp->GetReadPos();
		galp.resize(tables);
		for (z=0;z<tables;z++)
			galp[z] = z;
		i = 0;
		tia.clear();
		while (true) {
			z = 0;
			while (inp->ReadBit())
				z++;
			ti = galp[z];
			if (z != 0)
				MoveToFrontIndex(galp,z);
				//std::iter_swap(galp.begin(),galp.begin()+z);
			if ((i+1)*blocklength < dc)
				c = blocklength;
			else
				c = dc-(i*blocklength);
			for (z=0;z<c;z++)
				tia.push_back(hta[ti]->DecodeSymbol(inp));
			i++;
			if (i*blocklength >= dc)
				break;
		}

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("%d Bytes read.\n",(inp->GetReadPos()-k)/8);

		for (z=0;z<tables;z++)
			delete hta[z];

	} else {

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Importing Huffman tree...\n");
		i = inp->GetReadPos();
		ht = new CBQBHuffmanTree(m_IF);
		ht->ImportTree(inp,chr);
		if (m_IF.IsPL(BQB_PL_DEBUG)) {
			m_IF.printf("%d Bytes read.\n",(inp->GetReadPos()-i)/8);
			m_IF.printf("Importing and decoding data stream...\n");
		}
		bits = inp->ReadBitsInteger(6);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Symbol count stored with %d bits.\n",bits);
		dc = inp->ReadBitsInteger(bits);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Expecting %d symbols.\n",dc);
		k = inp->GetReadPos();
		tia.clear();
		for (z=0;z<dc;z++) {
			tia.push_back(ht->DecodeSymbol(inp));
	//		if (dc < 5)
	//			m_IF.printf("@@@ %d\n",(int)tia[tia.size()-1]);
		}

		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("%d Bytes (%d bits) read.\n",(inp->GetReadPos()-k)/8,inp->GetReadPos()-k);

		delete ht;
	}

	tiatrans.assign(tia.begin(),tia.end());

	if (coderun) {
		tia.clear();
		iruna = alp->FindIndex(C_RUNA);
		irunb = alp->FindIndex(C_RUNB);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Decoding code runs (RUNA=%d, RUNB=%d)...\n",iruna,irunb);
		ti = 1;
		ti2 = 0;
		for (z=0;z<(int)tiatrans.size();z++) {
			if (tiatrans[z] == iruna) {
				ti2 += ti;
				ti *= 2;
			} else if (tiatrans[z] == irunb) {
				ti2 += 2*ti;
				ti *= 2;
			} else {
				ti = 1;
				for (z2=0;z2<ti2;z2++)
					tia.push_back(0);
				ti2 = 0;
		//		if (tia.size() == 4040534)
		//			m_IF.printf("@@@ Stems from rle position %d.\n",z);
				tia.push_back(tiatrans[z]);
			}
		}
		for (z2=0;z2<ti2;z2++)
			tia.push_back(0);
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Done. Expanded %lu symbols to %lu symbols.\n",(unsigned long)tiatrans.size(),(unsigned long)tia.size());
		tiatrans.assign(tia.begin(),tia.end());
	}

	if (mtf) {
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Performing inverse Move-to-Front transform...\n");
		tia.assign(tiatrans.begin(),tiatrans.end());
		galp.resize(alp->m_oaAlphabet.size());
		for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
			galp[z] = z;
		for (z=0;z<(int)tia.size();z++) {
//			if ((z >= 4040534-10) && (z <= 4040534+10))
//				m_IF.printf("%9d: %4d --> %9d\n",z,tia[z],galp[tia[z]]);
			if (tia[z] != 0)
				MoveToFrontIndex(galp,tia[z]);
				//std::iter_swap(galp.begin(),galp.begin()+tia[z]);
			tiatrans[z] = galp[0];
		}
	}

	if (bw) {
		if (m_IF.IsPL(BQB_PL_DEBUG))
			m_IF.printf("Performing inverse Burrows-Wheeler transform...\n");
		tia.assign(tiatrans.begin(),tiatrans.end());
		tbwpairs.clear();
		for (z=0;z<(int)tia.size();z++)
			tbwpairs.push_back(CBQBBWPair(tia[z],z));
		std::stable_sort(tbwpairs.begin(),tbwpairs.end(),SORT_InverseBW);
		ti = bwindex;
		for (z=0;z<(int)tia.size();z++) {
			tiatrans[z] = tbwpairs[ti].m_iSymbol;
			ti = tbwpairs[ti].m_iIndex;
		}
	}

	if (m_IF.IsPL(BQB_PL_DEBUG))
		m_IF.printf("Converting indices to original alphabet symbols...\n");

	outp.reserve(outp.size()+tiatrans.size());

	for (z=0;z<(int)tiatrans.size();z++)
		outp.push_back(alp->m_oaAlphabet[tiatrans[z]]->m_iSymbol);

	delete alp;

	if (m_IF.IsPL(BQB_PL_DEBUG)) {
		m_IF.printf("\n");
		m_IF.printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
		m_IF.printf("<<<   DecompressSingle Done   <<<\n");
		m_IF.printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
		m_IF.printf("\n");
	}

	return true;
}


/*
void CBQBIntegerEngine::MultiHuffmanOptimize(int tables, std::vector<CBQBHuffmanTree*> &hta, CBQBAlphabet *alp, std::vector<int> &tia, std::vector<int> &asi, std::vector<int> &iasi, int blocklength, bool verbose) {

	UNUSED(tables);
	UNUSED(asi);
	UNUSED(iasi);

	//std::vector<std::vector<int> > tmatrix;
	std::vector<int> galp;
	int z, z2, z0, ti, ti2, tii, zi;
	int sc, payload, trees, blockendsum;
	CBQBHuffmanTree ht(m_IF), *pht;
	CBQBHuffmanEstimator he;
	double nsizeold, nsizenew, nsize, minval;
	std::vector<int> tfra, tsize, tpop, imintfra;
	std::vector<std::vector<int> > told;
	int lasttree, treemalus, symbmalus, ii, switchmalussum, iminval, lookahead, lastpos, imintii, iminsc, iskip;
	int blc, iminti, iminti2, ti4;
	//int z3, za, tibest, ti3;
	//int mic, mac, mam, last;
	//int blockcount, tblockcount;
	//bool b;


	m_IF.eprintf("CBQBIntegerEngine::MultiHuffmanOptimize(): Internal error: This code is not functional.\n\n");
	abort();

	bool switchrl = true;

	iminti = 0; // Avoid uninitialized variable compiler warning

//	blocklength = 50;
	treemalus = 60;
	symbmalus = 3;
	lookahead = 100000;

	if (verbose) {
		m_IF.printf("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
		m_IF.printf(">>>   MultiHuffmanOptimize   >>>\n");
		m_IF.printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n");
	}

	if (blocklength != 0) {

		hta.clear();

		he.Init((int)alp->m_oaAlphabet.size());

		tfra.resize(alp->m_oaAlphabet.size());

		for (zi=0;zi<1;zi++) {

			if (verbose)
				m_IF.printf("*** Iteration %d ***\n",zi+1);

			lastpos = 0;
			blc = 0;
			lasttree = -1;
			switchmalussum = 0;
			nsizeold = 0;

			for (z=0;z<(int)hta.size();z++) {
				for (z2=0;z2<(int)alp->m_oaAlphabet.size();z2++)
					hta[z]->m_iaFrequencies[z2] = 0;
				tsize[z] = he.EstimateBitLength(told[z]);
				tpop[z] = 0;
				galp[z] = z;
			}

			while (true) {

				if (lastpos+blocklength < (int)tia.size())
					sc = blocklength;
				else
					sc = (int)tia.size()-lastpos;

				for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
					tfra[z] = 0;

				for (z0=lastpos;z0<lastpos+sc;z0++)
					tfra[tia[z0]]++;

				ti = 1000000000;
				tii = -1;
				for (z=0;z<(int)galp.size();z++) {

					ii = galp[z];

					if ((int)told.size() > ii)
						ti2 = he.EstimateBitLength(hta[ii]->m_iaFrequencies,tfra,told[ii]);
					else
						ti2 = he.EstimateBitLength(hta[ii]->m_iaFrequencies,tfra);

					if (ti2 - tsize[ii] < ti) {
						ti = ti2 - tsize[ii];
						tii = z;
					}
				}

				if (galp.size() != 0) {
			//		if (switchrl) {
						nsizeold = (double)(ti + tii + 1) / sc;
				}

				ti2 = he.EstimateBitLength(tfra);

			//	if (switchrl)
					nsizenew = (double)(ti2 + treemalus + hta.size() + 1) / sc;
			//	else
			//		nsizenew = (double)(ti2 + treemalus + 8) / sc;

		//		m_IF.printf("Block %5d: New %.2f, old %.2f (%d)\n",blc,nsizenew,nsizeold,tii);

				if ((nsizeold < nsizenew) && (galp.size() != 0)) {
					nsize = nsizeold;
				} else {
					nsize = nsizenew;
					tii = -1;
				}

				blc++;

				if (verbose)
					if ((blc % 100) == 0)
						m_IF.printf("%7.3f%% (%lu tables)...\n",(double)lastpos/tia.size()*100.0,(unsigned long)hta.size());

				if (tii == -1) {
					pht = new CBQBHuffmanTree(m_IF);
					pht->Init((int)alp->m_oaAlphabet.size());
					pht->m_iaFrequencies.assign(tfra.begin(),tfra.end());
					lasttree = (int)hta.size();
					hta.push_back(pht);
					tsize.push_back(ti2);
					tpop.push_back(sc);
					galp.insert(galp.begin(),lasttree);
			//		if (switchrl) {
						switchmalussum += lasttree+1;
				} else {
					for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
						hta[galp[tii]]->m_iaFrequencies[z] += tfra[z];
					tsize[galp[tii]] += ti;
					tpop[galp[tii]] += sc;
					lasttree = tii;
					if (lasttree != 0) {
						MoveToFrontIndex(galp,lasttree);
				//		if (switchrl)
							switchmalussum += lasttree+1;
				//		else
				//			switchmalussum += 8;
					} else
						switchmalussum++;
				}

				lastpos += sc;

				if (lastpos >= (int)tia.size())
					break;
			}

			told.resize(hta.size());
			for (z=0;z<(int)hta.size();z++) {
				told[z].resize(alp->m_oaAlphabet.size());
				for (z2=0;z2<(int)alp->m_oaAlphabet.size();z2++)
					told[z][z2] = hta[z]->m_iaFrequencies[z2];
			}

			payload = 0;
			trees = 0;
			blockendsum = 0;

			for (z=0;z<(int)hta.size();z++)
				hta[z]->BuildTree(true);

			if (verbose) {
//				ti3 = 0;
				for (z=0;z<(int)hta.size();z++) {
					ti = 0;
					ti2 = 0;
					for (z2=0;z2<(int)alp->m_oaAlphabet.size();z2++) {
						ti += hta[z]->m_iaFrequencies[z2];
						ti2 += hta[z]->m_iaFrequencies[z2] * hta[z]->m_iaLengths[z2];
					}
					payload += ti2;
					ti4 = hta[z]->ExportTree(NULL,false,false);
					trees += ti4;
					m_IF.printf("  * Tree %3d: %4lu/%lu symbols, %7d payload = %8.3f KiB, %4d ends = %7.3f KiB, tree = %6.3f KiB\n",z,(unsigned long)hta[z]->m_oaSymbols.size(),(unsigned long)alp->m_oaAlphabet.size()+1,ti,ti2/8.0/1024.0,0,0.0,ti4/8.0/1024.0);
				}
				m_IF.printf("\n");
				m_IF.printf("Total switching malus:  %7.0f Bytes (%7.3f%%).\n",switchmalussum/8.0,(double)switchmalussum/(payload+trees+switchmalussum+blockendsum)*100.0);
				m_IF.printf("Total block end label:  %7.0f Bytes (%7.3f%%).\n",blockendsum/8.0,(double)blockendsum/(payload+trees+switchmalussum+blockendsum)*100.0);
				m_IF.printf("Total Payload:          %7.0f Bytes (%7.3f%%).\n",payload/8.0,(double)payload/(payload+trees+switchmalussum+blockendsum)*100.0);
				m_IF.printf("Total Trees:            %7.0f Bytes (%7.3f%%).\n",trees/8.0,(double)trees/(payload+trees+switchmalussum+blockendsum)*100.0);
				m_IF.printf("\n");
				m_IF.printf("Total size:             %7.0f Bytes.\n",(payload+trees+switchmalussum+blockendsum)/8.0);
			}
		}

	} else {

		hta.clear();

		he.Init((int)alp->m_oaAlphabet.size()+1);
		tfra.resize(alp->m_oaAlphabet.size()+1);
		for (z=0;z<(int)alp->m_oaAlphabet.size()+1;z++)
			tfra[z] = 0;
		tfra[0] = 1;

		sc = 0;
		minval = 1.0e20;
		lasttree = -1;
		switchmalussum = 0;
		iskip = 0;
		blc = 0;
		lastpos = 0;
		nsizeold = 0;
		iminval = 0;
		imintii = -1;
		iminsc = 0;

		for (z0=0;z0<(int)tia.size();z0++) {

			tfra[tia[z0]+1]++;
	//		tfra[tia[z0]]++;
			sc++;

			ti = 1000000000;
			tii = -1;
			for (z=0;z<(int)galp.size();z++) {

				ii = galp[z];

				ti2 = he.EstimateBitLength(hta[ii]->m_iaFrequencies,tfra);

				for (z2=0;z2<(int)alp->m_oaAlphabet.size()+1;z2++)
	//			for (z2=0;z2<(int)alp->m_oaAlphabet.size();z2++)
	//				if ((hta[ii]->m_iaFrequencies[z2] != 0) || (tfra[z2] != 0))
					if ((hta[ii]->m_iaFrequencies[z2] == 0) && (tfra[z2] != 0))
						ti2 += symbmalus;

				if (ti2 - tsize[ii] < ti) {
					ti = ti2 - tsize[ii];
					tii = z;
				}
			}

			if (galp.size() != 0) {
				if (switchrl) {
					nsizeold = (double)(ti + tii + 1) / sc;
				} else {
					if (tii == 0)
						nsizeold = (double)(ti + 1) / sc;
					else
						nsizeold = (double)(ti + 8) / sc;
				}
			}

			ti2 = he.EstimateBitLength(tfra);

			for (z=0;z<(int)alp->m_oaAlphabet.size()+1;z++)
	//		for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
				if (tfra[z] != 0)
					ti2 += symbmalus;

			if (switchrl)
				nsizenew = (double)(ti2 + treemalus + hta.size() + 1) / sc;
			else
				nsizenew = (double)(ti2 + treemalus + 8) / sc;

			if ((nsizeold < nsizenew) && (galp.size() != 0)) {
				nsize = nsizeold;
			} else {
				nsize = nsizenew;
				tii = -1;
			}

			if (nsize < minval) {
				iskip = 0;
				minval = nsize;
				iminval = z0;
				imintii = tii;
				iminsc = sc;
				iminti = ti;
				iminti2 = ti2;
				imintfra.assign(tfra.begin(),tfra.end());
			} else
				iskip++;

			if ((z0-lastpos == lookahead) || (iskip == 500) || (z0+1 == (int)tia.size())) {

				blc++;

				if (imintii == -1) {
					pht = new CBQBHuffmanTree(m_IF);
					pht->Init((int)alp->m_oaAlphabet.size()+1);
	//				pht->Init(alp->m_oaAlphabet.size());
					pht->m_iaFrequencies.assign(imintfra.begin(),imintfra.end());
					lasttree = (int)hta.size();
					hta.push_back(pht);
					tsize.push_back(iminti2);
					tpop.push_back(iminsc);
					galp.insert(galp.begin(),lasttree);
					if (switchrl) {
						switchmalussum += lasttree+1;
					} else {
						if (lasttree != 0)
							switchmalussum += 8;
						else
							switchmalussum++;
					}
				} else {
					for (z=0;z<(int)alp->m_oaAlphabet.size()+1;z++)
	//				for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
						hta[galp[imintii]]->m_iaFrequencies[z] += imintfra[z];
					tsize[galp[imintii]] += iminti ;
					tpop[galp[imintii]] += iminsc;
					lasttree = imintii;
					if (lasttree != 0) {
						MoveToFrontIndex(galp,lasttree);
						//std::iter_swap(galp.begin(),galp.begin()+lasttree);
						if (switchrl)
							switchmalussum += 8;
						else
							switchmalussum += lasttree+1;
					} else
						switchmalussum++;
				}

				if ((blc % 10) == 0) {
					if (imintii == -1)
						m_IF.printf("*** Block %d (pos %.3f%%, len %d) new tree ***\n",blc,(double)iminval/tia.size()*100.0,iminval-lastpos+1);
					else
						m_IF.printf("*** Block %d (pos %.3f%%, len %d) using tree %d (%d) ***\n",blc,(double)iminval/tia.size()*100.0,iminval-lastpos+1,imintii,galp[0]);
					for (z=0;z<(int)hta.size();z++)
						m_IF.printf("  * Tree %3d: %7d symbols, %8d bits, %6.3f b/s.\n",z,tpop[z],tsize[z],(double)tsize[z]/tpop[z]);
					m_IF.printf("    Total switching malus: %d bits.\n\n",switchmalussum);
				}

		//		if (iminval-lastpos+1 == 1)
		//			abort();

				for (z=1;z<(int)alp->m_oaAlphabet.size()+1;z++)
	//			for (z=0;z<(int)alp->m_oaAlphabet.size();z++)
					tfra[z] = 0;
				tfra[0] = 1;
				sc = 0;

				z0 = iminval;
				lastpos = iminval+1;
				minval = 1.0e20;
			}
		}
	}

	for (z=0;z<(int)hta.size();z++)
		hta[z]->BuildTree(true);

	if (verbose) {
		m_IF.printf("\n*** Loop finished ***\n\n");
		payload = 0;
		trees = 0;
		blockendsum = 0;


		m_IF.printf("\n<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
		m_IF.printf("<<<   MultiHuffmanOptimize   <<<\n");
		m_IF.printf("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n");
	}

//	abort();
}

*/






