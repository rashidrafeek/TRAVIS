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

#include <sstream>
#include <iomanip>
#include "bqb_parmset.h"
#include "bqb_tools.h"
#include "bqb_interface.h"


const char *GetRevisionInfo_bqb_parmset(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_parmset() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}




void BQBChkRngInt(CBQBInterface &in, int &i, int mi, int ma, const char *s) {

	if ((i < mi) || (i > ma)) {
		in.eprintf("%s: Warning: Allowed range is %d .. %d, found %d. Clipping.\n",s,mi,ma,i);
		if (i < mi)
			i = mi;
		if (i > ma)
			i = ma;
	}
}


void BQBChkRngFlt(CBQBInterface &in, double &f, double mi, double ma, const char *s) {

	if ((f < mi) || (f > ma)) {
		in.eprintf("%s: Warning: Allowed range is %.3f .. %.3f, found %.3f. Clipping.\n",s,mi,ma,f);
		if (f < mi)
			f = mi;
		if (f > ma)
			f = ma;
	}
}


void CBQBParameterSet_Compressor::ExportBits(CBQBBitSet *bs) const {

	bs->WriteBit(m_bOptTables);
	bs->WriteBit(m_bRLE);
	bs->WriteBit(m_bBW);
	bs->WriteBit(m_bMTF);
	bs->WriteBits(m_iTableCount,6);
	if (m_iTableCount != 0) {
		bs->WriteBits(m_iBlock,7);
		bs->WriteBits(m_iMaxIter,6);
	}
	bs->WriteBits(m_iMaxChunk,4);
}


bool CBQBParameterSet_Compressor::ImportBits(CBQBBitSet *bs, int version) {

	if (version > 0) {
		m_IF.eprintf("CBQBParameterSet_Compressor::ImportBits(): Invalid parameter set version: %d (allowed: 0).\n",version);
		m_IF.eprintf("Probably the LibBQB version is too old for this parameter set.\n");
		return false;
	}

	m_bOptTables = bs->ReadBit();
	m_bRLE = bs->ReadBit();
	m_bBW = bs->ReadBit();
	m_bMTF = bs->ReadBit();
	m_iTableCount = bs->ReadBitsInteger(6);
	if (m_iTableCount != 0) {
		m_iBlock = bs->ReadBitsInteger(7);
		m_iMaxIter = bs->ReadBitsInteger(6);
	} else {
		m_iBlock = 1;
		m_iMaxIter = 0;
	}
	m_iMaxChunk = bs->ReadBitsInteger(4);

	return true;
}


std::string CBQBParameterSet_Compressor::ToString(int indent) const {

	std::string i;
	int z;
	std::ostringstream o;

	for (z=0;z<indent;z++)
		i += " ";

	o << i << "Maximum Chunk Length:      -maxchunk " << GetMaxChunk() << "\n";
	o << i << "Huffman Table Count:       -tables " << GetTableCount() << "\n";
	if (GetTableCount() > 1) {
		o << i << "Huffman Block Length:      -block " << GetBlockLength() << "\n";
		o << i << "Table Opt. Iterations:     -iter " << GetMaxIter() << "\n";
	}
	o << i << "Optimize Table Count:      -opttables " << (GetOptTables()?"yes":"no") << "\n";
	o << i << "Burrows-Wheeler:           -bw " << (GetBW()?"yes":"no") << "\n";
	o << i << "Move-to-Front:             -mtf " << (GetMTF()?"yes":"no") << "\n";
	o << i << "Zero Run-length Encoding:  -rle " << (GetRLE()?"yes":"no") << "\n";

	return o.str();
}


std::string CBQBParameterSet_Compressor::ToKey() const {

	CBQBBitSet bs(m_IF);

	bs.WriteBits(0,4); // Version 0
	bs.WriteBit(0); // Position data
	bs.WriteBit(0); // Volumetric data

	ExportBits(&bs);

	return bs.ExportKey();
}


bool CBQBParameterSet_Compressor::FromKey(std::string key) {

	CBQBBitSet bs(m_IF);
	bool bp, bv;
	int ver;

	if (!bs.ImportKey(key)) {
		m_IF.eprintf("CBQBParameterSet_Compressor::FromKey(): Error importing key.\n");
		return false;
	}

	ver = bs.ReadBitsInteger(4);
	if (ver > 0) {
		m_IF.eprintf("CBQBParameterSet_Compressor::FromKey(): Invalid parameter set version: %d (allowed: 0).\n",ver);
		m_IF.eprintf("Probably the LibBQB version is too old for this key.\n");
		return false;
	}

	bp = bs.ReadBit();
	bv = bs.ReadBit();
	if (bp || bv) {
		m_IF.eprintf("CBQBParameterSet_Compressor::FromKey(): Error: This key is not for compressor data.\n");
		return false;
	}

	if (!ImportBits(&bs,ver)) {
		m_IF.eprintf("CBQBParameterSet_Compressor::FromKey(): Error importing compressor parameters from key.\n");
		return false;
	}

	return true;
}



/*******************************************************************************************************************/



void CBQBParameterSet_Position::ExportBits(CBQBBitSet *bs) const {

	m_oCompressor.ExportBits(bs);

	bs->WriteBits(m_iPrec,3);
	bs->WriteBits(m_iSplit,5);
	bs->WriteBit(m_bSortAtom);

	bs->WriteBit(m_bUseExtra);

	if (m_bUseExtra) {
		bs->WriteBits(m_iExtraTRange,4);
		bs->WriteBits(m_iExtraTOrder,4);
		bs->WriteBits(m_iExtraTimeExpo,14);
	} else {
		bs->WriteBit(m_bOptOrder);
		bs->WriteBits(m_iOrder,4);
	}
}


bool CBQBParameterSet_Position::ImportBits(CBQBBitSet *bs, int version) {

	if (version > 0) {
		m_IF.eprintf("CBQBParameterSet_Position::ImportBits(): Invalid parameter set version: %d (allowed: 0).\n",version);
		m_IF.eprintf("Probably the LibBQB version is too old for this parameter set.\n");
		return false;
	}

	if (!m_oCompressor.ImportBits(bs,version)) {
		m_IF.eprintf("CBQBParameterSet_Position::ImportBits(): Error reading compressor parameters.\n");
		return false;
	}

	m_iPrec = bs->ReadBitsInteger(3);
	m_iSplit = bs->ReadBitsInteger(5);
	m_bSortAtom = bs->ReadBit();

	m_bUseExtra = bs->ReadBit();

	if (m_bUseExtra) {
		m_iExtraTRange = bs->ReadBitsInteger(4);
		m_iExtraTOrder = bs->ReadBitsInteger(4);
		m_iExtraTimeExpo = bs->ReadBitsInteger(14);
	} else {
		m_bOptOrder = bs->ReadBit();
		m_iOrder = bs->ReadBitsInteger(4);
	}

	return true;
}


std::string CBQBParameterSet_Position::ToString(int indent) const {

	std::string i;
	int z;
	std::ostringstream o;

	for (z=0;z<indent;z++)
		i += " ";

	o << i << "Concerning the position prediction:" << "\n";
	o << i << "  Absolute Precision:        -prec " << GetPrecision() << "\n";
	o << i << "  Extrapolation Order:       -order " << GetOrder() << "\n";
	o << i << "  Optimize Order:            -optorder " << (GetOptOrder()?"yes":"no") << "\n";
	o << i << "  Symbol Modulo:             -split " << GetSplit() << "\n";
	o << i << "  Sort Atoms by Mass:        -sortatoms " << (GetSortAtom()?"yes":"no") << "\n";
	o << i << "  Use Extrapolator:          -extra " << (GetUseExtra()?"yes":"no") << "\n";
	o << i << "  Temporal Range:            -extrange " << GetExtraTRange() << "\n";
	o << i << "  Max. Temporal Order:       -extorder " << GetExtraTOrder() << "\n";
	o << i << "  Time Weighting Exponent:   -extimeexpo " << std::fixed << std::setprecision(3) << GetExtraTimeExpo() << "\n";
	
	o << i << "Concerning the position compression:" << "\n";
	o << i << "  Maximum Chunk Length:      -maxchunk " << GetMaxChunk() << "\n";
	o << i << "  Huffman Table Count:       -tables " << GetTableCount() << "\n";
	o << i << "  Huffman Block Length:      -block " << GetBlockLength() << "\n";
	o << i << "  Table Opt. Iterations:     -iter " << GetMaxIter() << "\n";
	o << i << "  Optimize Table Count:      -opttables " << (GetOptTables()?"yes":"no") << "\n";
	o << i << "  Burrows-Wheeler:           -bw " << (GetBW()?"yes":"no") << "\n";
	o << i << "  Move-to-Front:             -mtf " << (GetMTF()?"yes":"no") << "\n";
	o << i << "  Zero Run-length Encoding:  -rle " << (GetRLE()?"yes":"no") << "\n";

	return o.str();
}


std::string CBQBParameterSet_Position::ToKey() const {

	CBQBBitSet bs(m_IF);

	bs.WriteBits(0,4); // Version 0
	bs.WriteBit(1); // Position data
	bs.WriteBit(0); // Volumetric data

	ExportBits(&bs);

	return bs.ExportKey();
}


bool CBQBParameterSet_Position::FromKey(std::string key) {

	CBQBBitSet bs(m_IF);
	bool bp, bv;
	int ver;

	if (!bs.ImportKey(key)) {
		m_IF.eprintf("CBQBParameterSet_Position::FromKey(): Error importing key.\n");
		return false;
	}

	ver = bs.ReadBitsInteger(4);
	if (ver > 0) {
		m_IF.eprintf("CBQBParameterSet_Position::FromKey(): Invalid parameter set version: %d (allowed: 0).\n",ver);
		m_IF.eprintf("Probably the LibBQB version is too old for this key.\n");
		return false;
	}

	bp = bs.ReadBit();
	bv = bs.ReadBit();
	if (!bp || bv) {
		m_IF.eprintf("CBQBParameterSet_Position::FromKey(): Error: This key is not for position data.\n");
		return false;
	}

	if (!ImportBits(&bs,ver)) {
		m_IF.eprintf("CBQBParameterSet_Position::FromKey(): Error importing position parameters from key.\n");
		return false;
	}

	return true;
}



/*******************************************************************************************************************/



void CBQBParameterSet_Volumetric::ExportBits(CBQBBitSet *bs) const {

	m_oCompressor.ExportBits(bs);

	bs->WriteBits(m_iSigni,4);
	bs->WriteBits(m_iEps,6);
	bs->WriteBits(m_iSplit,5);
	bs->WriteBit(m_bHilbert);

	bs->WriteBit(m_bUseExtra);

	if (m_bUseExtra) {

		if ((m_iExtraSRangeX == m_iExtraSRangeY) && (m_iExtraSRangeX == m_iExtraSRangeZ)) {
			bs->WriteBit(0);
			bs->WriteBits(m_iExtraSRangeX,4);
		} else {
			bs->WriteBit(1);
			bs->WriteBits(m_iExtraSRangeX,4);
			bs->WriteBits(m_iExtraSRangeY,4);
			bs->WriteBits(m_iExtraSRangeZ,4);
		}

		bs->WriteBits(m_iExtraTRange,4);
		bs->WriteBits(m_iExtraSOrder,4);
		bs->WriteBits(m_iExtraTOrder,4);

		if ((m_iExtraOffsetX == m_iExtraOffsetY) && (m_iExtraOffsetX == m_iExtraOffsetZ)) {
			bs->WriteBit(0);
			bs->WriteBits(m_iExtraOffsetX,4);
		} else {
			bs->WriteBit(1);
			bs->WriteBits(m_iExtraOffsetX,4);
			bs->WriteBits(m_iExtraOffsetY,4);
			bs->WriteBits(m_iExtraOffsetZ,4);
		}

		bs->WriteBit(m_bExtraCrossS);
		bs->WriteBit(m_bExtraCrossT);
		bs->WriteBit(m_bExtraWrap);
		bs->WriteBit(m_bExtraCrossRangeS);
		bs->WriteBit(m_bExtraCrossRangeT);
		bs->WriteBits(m_iExtraDistExpo,14);
		bs->WriteBits(m_iExtraTimeExpo,14);
		bs->WriteBit(m_bExtraPredCorr);

	} else {

		bs->WriteBits(m_iOrder,4);
		bs->WriteBit(m_bOptOrder);
		bs->WriteBits(m_iNbhFac,12);
	}
}

	
bool CBQBParameterSet_Volumetric::ImportBits(CBQBBitSet *bs, int version) {

	if (version > 0) {
		m_IF.eprintf("CBQBParameterSet_Volumetric::ImportBits(): Invalid parameter set version: %d (allowed: 0).\n",version);
		m_IF.eprintf("Probably the LibBQB version is too old for this parameter set.\n");
		return false;
	}

	if (!m_oCompressor.ImportBits(bs,version)) {
		m_IF.eprintf("CBQBParameterSet_Volumetric::ImportBits(): Error reading compressor parameters.\n");
		return false;
	}

	m_iSigni = bs->ReadBitsInteger(4);
	m_iEps = bs->ReadBitsInteger(6);
	m_iSplit = bs->ReadBitsInteger(5);
	m_bHilbert = bs->ReadBit();

	m_bUseExtra = bs->ReadBit();

	if (m_bUseExtra) {

		if (bs->ReadBit()) {
			m_iExtraSRangeX = bs->ReadBitsInteger(4);
			m_iExtraSRangeY = bs->ReadBitsInteger(4);
			m_iExtraSRangeZ = bs->ReadBitsInteger(4);
		} else {
			m_iExtraSRangeX = bs->ReadBitsInteger(4);
			m_iExtraSRangeY = m_iExtraSRangeX;
			m_iExtraSRangeZ = m_iExtraSRangeX;
		}

		m_iExtraTRange = bs->ReadBitsInteger(4);
		m_iExtraSOrder = bs->ReadBitsInteger(4);
		m_iExtraTOrder = bs->ReadBitsInteger(4);

		if (bs->ReadBit()) {
			m_iExtraOffsetX = bs->ReadBitsInteger(4);
			m_iExtraOffsetY = bs->ReadBitsInteger(4);
			m_iExtraOffsetZ = bs->ReadBitsInteger(4);
		} else {
			m_iExtraOffsetX = bs->ReadBitsInteger(4);
			m_iExtraOffsetY = m_iExtraOffsetX;
			m_iExtraOffsetZ = m_iExtraOffsetX;
		}

		m_bExtraCrossS = bs->ReadBit();
		m_bExtraCrossT = bs->ReadBit();
		m_bExtraWrap = bs->ReadBit();
		m_bExtraCrossRangeS = bs->ReadBit();
		m_bExtraCrossRangeT = bs->ReadBit();
		m_iExtraDistExpo = bs->ReadBitsInteger(14);
		m_iExtraTimeExpo = bs->ReadBitsInteger(14);
		m_bExtraPredCorr = bs->ReadBit();

	} else {

		m_iOrder = bs->ReadBitsInteger(4);
		m_bOptOrder = bs->ReadBit();
		m_iNbhFac = bs->ReadBitsInteger(12);
	}

	return true;
}


std::string CBQBParameterSet_Volumetric::ToString(int indent) const {

	std::string i;
	int z;
	std::ostringstream o;


	for (z=0;z<indent;z++)
		i += " ";

	o << i << "Concerning the volumetric data prediction:" << "\n";
	o << i << "  Significant Digits:        -vsigni " << GetSigni() << "\n";
	o << i << "  Smallest Exponent:         -veps " << GetEps() << "\n";
	o << i << "  Extrapolation Order:       -vorder " << GetOrder() << "\n";
	o << i << "  Optimize Order:            -voptorder " << (GetOptOrder()?"yes":"no") << "\n";
	o << i << "  Hilbert Curve:             -vhilbert " << (GetHilbert()?"yes":"no") << "\n";
	o << i << "  Neighborhood Factor:       -vnbhfac " << GetNbhFac() << "\n";
	o << i << "  Symbol Modulo:             -vsplit " << GetSplit() << "\n";
	o << i << "  Use Extrapolator:          -vextra " << (GetUseExtra()?"yes":"no") << "\n";
	if (!IsExtraSRangeXYZEqual()) {
		o << i << "  Spatial X Range:           -vexsrangex " << GetExtraSRangeX() << "\n";
		o << i << "  Spatial Y Range:           -vexsrangey " << GetExtraSRangeY() << "\n";
		o << i << "  Spatial Z Range:           -vexsrangez " << GetExtraSRangeZ() << "\n";
	} else
		o << i << "  Spatial Range:             -vexsrangez " << GetExtraSRangeX() << "\n";
	o << i << "  Temporal Range:            -vextrange " << GetExtraTRange() << "\n";
	o << i << "  Max. Spatial Order:        -vexsorder " << GetExtraSOrder() << "\n";
	o << i << "  Max. Temporal Order:       -vextorder " << GetExtraTOrder() << "\n";
	if (!IsExtraOffsetXYZEqual()) {
		o << i << "  Spatial X Offset:          -vexoffsetx " << GetExtraOffsetX() << "\n";
		o << i << "  Spatial Y Offset:          -vexoffsety " << GetExtraOffsetY() << "\n";
		o << i << "  Spatial Z Offset:          -vexoffsetz " << GetExtraOffsetZ() << "\n";
	} else
		o << i << "  Spatial Offset:            -vexoffsetz " << GetExtraOffsetX() << "\n";
	o << i << "  Spatial Cross-Terms:       -vexscross " << (GetExtraCrossS()?"yes":"no") << "\n";
	o << i << "  S/T Cross-Terms:           -vextcross " << (GetExtraCrossT()?"yes":"no") << "\n";
	o << i << "  Spatial Wrapping:          -vexwrap " << (GetExtraWrap()?"yes":"no") << "\n";
	o << i << "  Spatial Diagonal Data:     -vexscrossrange " << (GetExtraCrossRangeS()?"yes":"no") << "\n";
	o << i << "  S/T Diagonal Data:         -vextcrossrange " << (GetExtraCrossRangeT()?"yes":"no") << "\n";
	o << i << "  Dist. Weighting Exponent:  -vexdistexpo " << std::fixed << std::setprecision(3) << GetExtraDistExpo() << "\n";
	o << i << "  Time Weighting Exponent:   -vextimeexpo " << std::fixed << std::setprecision(3) << GetExtraTimeExpo() << "\n";
	o << i << "  Predictor-Corrector:       -vexpredcorr " << (GetExtraPredCorr()?"yes":"no") << "\n";

	o << i << "Concerning the volumetric data compression:" << "\n";
	o << i << "  Maximum Chunk Length:      -vmaxchunk " << GetMaxChunk() << "\n";
	o << i << "  Huffman Table Count:       -vtables " << GetTableCount() << "\n";
	o << i << "  Huffman Block Length:      -vblock " << GetBlockLength() << "\n";
	o << i << "  Table Opt. Iterations:     -viter " << GetMaxIter() << "\n";
	o << i << "  Optimize Table Count:      -vopttables " << (GetOptTables()?"yes":"no") << "\n";
	o << i << "  Burrows-Wheeler:           -vbw " << (GetBW()?"yes":"no") << "\n";
	o << i << "  Move-to-Front:             -vmtf " << (GetMTF()?"yes":"no") << "\n";
	o << i << "  Zero Run-length Encoding:  -vrle " << (GetRLE()?"yes":"no") << "\n";

	return o.str();
}


std::string CBQBParameterSet_Volumetric::ToKey() const {

	CBQBBitSet bs(m_IF);


	bs.WriteBits(0,4); // Version 0
	bs.WriteBit(0); // Position data
	bs.WriteBit(1); // Volumetric data

	ExportBits(&bs);

	return bs.ExportKey();
}


bool CBQBParameterSet_Volumetric::FromKey(std::string key) {

	CBQBBitSet bs(m_IF);
	bool bp, bv;
	int ver;


	if (!bs.ImportKey(key)) {
		m_IF.eprintf("CBQBParameterSet_Volumetric::FromKey(): Error importing key.\n");
		return false;
	}

	ver = bs.ReadBitsInteger(4);
	if (ver > 0) {
		m_IF.eprintf("CBQBParameterSet_Volumetric::FromKey(): Invalid parameter set version: %d (allowed: 0).\n",ver);
		m_IF.eprintf("Probably the LibBQB version is too old for this key.\n");
		return false;
	}

	bp = bs.ReadBit();
	bv = bs.ReadBit();
	if (bp || !bv) {
		m_IF.eprintf("CBQBParameterSet_Volumetric::FromKey(): Error: This key is not for volumetric data.\n");
		return false;
	}

	if (!ImportBits(&bs,ver)) {
		m_IF.eprintf("CBQBParameterSet_Volumetric::FromKey(): Error importing parameters from key.\n");
		return false;
	}

	return true;
}



/*******************************************************************************************************************/



void CBQBParameterSet_PosAndVol::ExportBits(CBQBBitSet *bs) const {

	m_oPosition.ExportBits(bs);

	m_oVolumetric.ExportBits(bs);
}

	
bool CBQBParameterSet_PosAndVol::ImportBits(CBQBBitSet *bs, int version) {

	if (version > 0) {
		m_IF.eprintf("CBQBParameterSet_PosAndVol::ImportBits(): Invalid parameter set version: %d (allowed: 0).\n",version);
		m_IF.eprintf("Probably the LibBQB version is too old for this parameter set.\n");
		return false;
	}

	if (!m_oPosition.ImportBits(bs,version)) {
		m_IF.eprintf("CBQBParameterSet_PosAndVol::ImportBits(): Error reading position parameters.\n");
		return false;
	}

	if (!m_oVolumetric.ImportBits(bs,version)) {
		m_IF.eprintf("CBQBParameterSet_PosAndVol::ImportBits(): Error reading volumetric parameters.\n");
		return false;
	}

	return true;
}


std::string CBQBParameterSet_PosAndVol::ToString(int indent) const {

	std::string i;
	int z;
	std::ostringstream o;

	for (z=0;z<indent;z++)
		i += " ";

	o << i << "Concerning the volumetric data prediction:" << "\n";
	o << i << "  Significant Digits:        -vsigni " << GetVolSigni() << "\n";
	o << i << "  Smallest Exponent:         -veps " << GetVolEps() << "\n";
	o << i << "  Extrapolation Order:       -vorder " << GetVolOrder() << "\n";
	o << i << "  Optimize Order:            -voptorder " << (GetVolOptOrder()?"yes":"no") << "\n";
	o << i << "  Hilbert Curve:             -vhilbert " << (GetVolHilbert()?"yes":"no") << "\n";
	o << i << "  Neighborhood Factor:       -vnbhfac " << GetVolNbhFac() << "\n";
	o << i << "  Symbol Modulo:             -vsplit " << GetVolSplit() << "\n";
	o << i << "  Use Extrapolator:          -vextra " << (GetVolUseExtra()?"yes":"no") << "\n";
	if (!IsVolExtraSRangeXYZEqual()) {
		o << i << "  Spatial X Range:           -vexsrangex " << GetVolExtraSRangeX() << "\n";
		o << i << "  Spatial Y Range:           -vexsrangey " << GetVolExtraSRangeY() << "\n";
		o << i << "  Spatial Z Range:           -vexsrangez " << GetVolExtraSRangeZ() << "\n";
	} else
		o << i << "  Spatial Range:             -vexsrangez " << GetVolExtraSRangeX() << "\n";
	o << i << "  Temporal Range:            -vextrange " << GetVolExtraTRange() << "\n";
	o << i << "  Max. Spatial Order:        -vexsorder " << GetVolExtraSOrder() << "\n";
	o << i << "  Max. Temporal Order:       -vextorder " << GetVolExtraTOrder() << "\n";
	if (!IsVolExtraOffsetXYZEqual()) {
		o << i << "  Spatial X Offset:          -vexoffsetx " << GetVolExtraOffsetX() << "\n";
		o << i << "  Spatial Y Offset:          -vexoffsety " << GetVolExtraOffsetY() << "\n";
		o << i << "  Spatial Z Offset:          -vexoffsetz " << GetVolExtraOffsetZ() << "\n";
	} else
		o << i << "  Spatial Offset:            -vexoffsetz " << GetVolExtraOffsetX() << "\n";
	o << i << "  Spatial Cross-Terms:       -vexscross " << (GetVolExtraCrossS()?"yes":"no") << "\n";
	o << i << "  S/T Cross-Terms:           -vextcross " << (GetVolExtraCrossT()?"yes":"no") << "\n";
	o << i << "  Spatial Wrapping:          -vexwrap " << (GetVolExtraWrap()?"yes":"no") << "\n";
	o << i << "  Spatial Diagonal Data:     -vexscrossrange " << (GetVolExtraCrossRangeS()?"yes":"no") << "\n";
	o << i << "  S/T Diagonal Data:         -vextcrossrange " << (GetVolExtraCrossRangeT()?"yes":"no") << "\n";
	o << i << "  Dist. Weighting Exponent:  -vexdistexpo " << std::fixed << std::setprecision(3) << GetVolExtraDistExpo() << "\n";
	o << i << "  Time Weighting Exponent:   -vextimeexpo " << std::fixed << std::setprecision(3) << GetVolExtraTimeExpo() << "\n";
	o << i << "  Predictor-Corrector:       -vexpredcorr " << (GetVolExtraPredCorr()?"yes":"no") << "\n";

	o << i << "Concerning the volumetric data compression:" << "\n";
	o << i << "  Maximum Chunk Length:      -vmaxchunk " << GetVolMaxChunk() << "\n";
	o << i << "  Huffman Table Count:       -vtables " << GetVolTableCount() << "\n";
	o << i << "  Huffman Block Length:      -vblock " << GetVolBlockLength() << "\n";
	o << i << "  Table Opt. Iterations:     -viter " << GetVolMaxIter() << "\n";
	o << i << "  Optimize Table Count:      -vopttables " << (GetVolOptTables()?"yes":"no") << "\n";
	o << i << "  Burrows-Wheeler:           -vbw " << (GetVolBW()?"yes":"no") << "\n";
	o << i << "  Move-to-Front:             -vmtf " << (GetVolMTF()?"yes":"no") << "\n";
	o << i << "  Zero Run-length Encoding:  -vrle " << (GetVolRLE()?"yes":"no") << "\n";

	o << i << "Concerning the position prediction:" << "\n";
	o << i << "  Absolute Precision:        -pprec " << GetPosPrecision() << "\n";
	o << i << "  Extrapolation Order:       -porder " << GetPosOrder() << "\n";
	o << i << "  Optimize Order:            -poptorder " << (GetPosOptOrder()?"yes":"no") << "\n";
	o << i << "  Symbol Modulo:             -psplit " << GetPosSplit() << "\n";
	o << i << "  Sort Atoms by Mass:        -psortatoms " << (GetPosSortAtom()?"yes":"no") << "\n";
	o << i << "  Use Extrapolator:          -pextra " << (GetPosUseExtra()?"yes":"no") << "\n";
	o << i << "  Temporal Range:            -pextrange " << GetPosExtraTRange() << "\n";
	o << i << "  Max. Temporal Order:       -pextorder " << GetPosExtraTOrder() << "\n";
	o << i << "  Time Weighting Exponent:   -pextimeexpo " << std::fixed << std::setprecision(3) << GetPosExtraTimeExpo() << "\n";
	
	o << i << "Concerning the position compression:" << "\n";
	o << i << "  Maximum Chunk Length:      -pmaxchunk " << GetPosMaxChunk() << "\n";
	o << i << "  Huffman Table Count:       -ptables " << GetPosTableCount() << "\n";
	o << i << "  Huffman Block Length:      -pblock " << GetPosBlockLength() << "\n";
	o << i << "  Table Opt. Iterations:     -piter " << GetPosMaxIter() << "\n";
	o << i << "  Optimize Table Count:      -popttables " << (GetPosOptTables()?"yes":"no") << "\n";
	o << i << "  Burrows-Wheeler:           -pbw " << (GetPosBW()?"yes":"no") << "\n";
	o << i << "  Move-to-Front:             -pmtf " << (GetPosMTF()?"yes":"no") << "\n";
	o << i << "  Zero Run-length Encoding:  -prle " << (GetPosRLE()?"yes":"no") << "\n";

	return o.str();
}


std::string CBQBParameterSet_PosAndVol::ToKey() const {

	CBQBBitSet bs(m_IF);

	bs.WriteBits(0,4); // Version 0
	bs.WriteBit(1); // Position data
	bs.WriteBit(1); // Volumetric data

	ExportBits(&bs);

	return bs.ExportKey();
}


bool CBQBParameterSet_PosAndVol::FromKey(std::string key) {

	CBQBBitSet bs(m_IF);
	bool bp, bv;
	int ver;

	if (!bs.ImportKey(key)) {
		m_IF.eprintf("CBQBParameterSet_PosAndVol::FromKey(): Error importing key.\n");
		return false;
	}

	ver = bs.ReadBitsInteger(4);
	if (ver > 0) {
		m_IF.eprintf("CBQBParameterSet_PosAndVol::FromKey(): Invalid parameter set version: %d (allowed: 0).\n",ver);
		m_IF.eprintf("Probably the LibBQB version is too old for this key.\n");
		return false;
	}

	bp = bs.ReadBit();
	bv = bs.ReadBit();
	if (!bp || !bv) {
		m_IF.eprintf("CBQBParameterSet_PosAndVol::FromKey(): Error: This key is not for position & volumetric data.\n");
		return false;
	}

	if (!ImportBits(&bs,ver)) {
		m_IF.eprintf("CBQBParameterSet_PosAndVol::FromKey(): Error importing parameters from key.\n");
		return false;
	}

	return true;
}


