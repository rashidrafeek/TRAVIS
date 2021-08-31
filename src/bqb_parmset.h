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



#ifndef BQB_PARMSET_H
#define BQB_PARMSET_H


// This must always be the first include directive
#include "bqb_config.h"

#include <string>
#include "bqb_bitset.h"



class CBQBInterface;



void BQBChkRngInt(CBQBInterface &in, int &i, int mi, int ma, const char *s);


void BQBChkRngFlt(CBQBInterface &in, double &f, double mi, double ma, const char *s);



// Current Version: 0
class CBQBParameterSet_Compressor {
friend class CBQBParameterSet_Position;
friend class CBQBParameterSet_Volumetric;
friend class CBQBParameterSet_PosAndVol;
public:

	// No pointers as members --> standard copy constructor and assignment operator will work

	explicit CBQBParameterSet_Compressor(CBQBInterface &i) : m_IF(i) {
		// Default values for file compression
		SetOptTables(false);
		SetRLE(true);
		SetBW(true);
		SetMTF(true);
		SetTableCount(6);
		SetBlockLength(50);
		SetMaxIter(10);
		SetMaxChunk(4194304);
		m_iVersion = 0;
	}

	void ExportBits(CBQBBitSet *bs) const;

	bool ImportBits(CBQBBitSet *bs, int version);

	std::string ToString(int indent) const;

	std::string ToKey() const;

	bool FromKey(std::string key);

	void SetOptTables(bool b) { m_bOptTables = b; }
	bool GetOptTables() const { return m_bOptTables; }

	void SetRLE(bool b) { m_bRLE = b; }
	bool GetRLE() const { return m_bRLE; }

	void SetBW(bool b) { m_bBW = b; }
	bool GetBW() const { return m_bBW; }

	void SetMTF(bool b) { m_bMTF = b; }
	bool GetMTF() const { return m_bMTF; }

	void SetTableCount(int i) { BQBChkRngInt(m_IF,i,1,64,"CBQBParameterSet_Compressor::SetTableCount()"); m_iTableCount = i-1; }
	int GetTableCount() const { return m_iTableCount+1; }

	void SetBlockLength(int i) { BQBChkRngInt(m_IF,i,1,128,"CBQBParameterSet_Compressor::SetBlockLength()"); m_iBlock = i-1; }
	int GetBlockLength() const { return m_iBlock+1; }

	void SetMaxIter(int i) { BQBChkRngInt(m_IF,i,0,126,"CBQBParameterSet_Compressor::SetMaxIter()"); m_iMaxIter = i/2; }
	int GetMaxIter() const { return m_iMaxIter*2; }

	void SetMaxChunk(int i) { BQBChkRngInt(m_IF,i,0,16777216,"CBQBParameterSet_Compressor::SetMaxChunk()"); if (i == 0) m_iMaxChunk = 0; else m_iMaxChunk = log2i(i/1024)+1; }
	int GetMaxChunk() const { if (m_iMaxChunk == 0) return 0; else return pow2i(m_iMaxChunk-1)*1024; }

	int GetVersion() const { return m_iVersion; }

private:
	CBQBInterface &m_IF;
	bool m_bOptTables;  //  1 Bit
	bool m_bRLE;        //  1 Bit
	bool m_bBW;         //  1 Bit
	bool m_bMTF;        //  1 Bit
	int m_iTableCount;  //  6 Bits, 1 .. 64, +1
	int m_iBlock;       //  7 Bits, 1 .. 128, + 1
	int m_iMaxIter;     //  6 Bits, 0 .. 126, * 2
	int m_iMaxChunk;    //  4 Bits, 1 kiSymbol .. 16384 kiSymbol, 2^(i-1), disable: i=0

	int m_iVersion;     //  4 Bits, 0 .. 15
};



// Current Version: 0
class CBQBParameterSet_Position {
friend class CBQBParameterSet_PosAndVol;
public:

	// No pointers as members --> standard copy constructor and assignment operator will work

	explicit CBQBParameterSet_Position(CBQBInterface &i) : m_IF(i), m_oCompressor(i) {

		SetPrecision(5);
		SetOrder(8);
		SetOptOrder(true);
		SetSplit(14);
		SetTableCount(1);
		SetOptTables(false);
		SetBlockLength(40);
		SetBW(false);
		SetMTF(false);
		SetRLE(true);
		SetMaxIter(10);
		SetSortAtom(true);
		SetMaxChunk(0);
		SetUseExtra(true);
		SetExtraTRange(9);
		SetExtraTOrder(6);
		SetExtraTimeExpo(4.0);
	}

	void ExportBits(CBQBBitSet *bs) const;
	bool ImportBits(CBQBBitSet *bs, int version);

	std::string ToString(int indent) const;

	std::string ToKey() const;

	bool FromKey(std::string key);

	void SetOptTables(bool b) { m_oCompressor.m_bOptTables = b; }
	bool GetOptTables() const { return m_oCompressor.m_bOptTables; }

	void SetRLE(bool b) { m_oCompressor.m_bRLE = b; }
	bool GetRLE() const { return m_oCompressor.m_bRLE; }

	void SetBW(bool b) { m_oCompressor.m_bBW = b; }
	bool GetBW() const { return m_oCompressor.m_bBW; }

	void SetMTF(bool b) { m_oCompressor.m_bMTF = b; }
	bool GetMTF() const { return m_oCompressor.m_bMTF; }

	void SetTableCount(int i) { BQBChkRngInt(m_IF,i,1,64,"CBQBParameterSet_Position::SetTableCount()"); m_oCompressor.m_iTableCount = i-1; }
	int GetTableCount() const { return m_oCompressor.m_iTableCount+1; }

	void SetBlockLength(int i) { BQBChkRngInt(m_IF,i,1,128,"CBQBParameterSet_Position::SetBlockLength()"); m_oCompressor.m_iBlock = i-1; }
	int GetBlockLength() const { return m_oCompressor.m_iBlock+1; }

	void SetMaxIter(int i) { BQBChkRngInt(m_IF,i,0,126,"CBQBParameterSet_Position::SetMaxIter()"); m_oCompressor.m_iMaxIter = i/2; }
	int GetMaxIter() const { return m_oCompressor.m_iMaxIter*2; }

	void SetMaxChunk(int i) { BQBChkRngInt(m_IF,i,0,16777216,"CBQBParameterSet_Position::SetMaxChunk()"); if (i == 0) m_oCompressor.m_iMaxChunk = 0; else m_oCompressor.m_iMaxChunk = log2i(i/1024)+1; }
	int GetMaxChunk() const { if (m_oCompressor.m_iMaxChunk == 0) return 0; else return pow2i(m_oCompressor.m_iMaxChunk-1)*1024; }



	void SetOptOrder(bool b) { m_bOptOrder = b; }
	bool GetOptOrder() const { return m_bOptOrder; }

	void SetSortAtom(bool b) { m_bSortAtom = b; }
	bool GetSortAtom() const { return m_bSortAtom; }

	void SetUseExtra(bool b) { m_bUseExtra = b; }
	bool GetUseExtra() const { return m_bUseExtra; }

	void SetPrecision(int i) { BQBChkRngInt(m_IF,i,0,7,"CBQBParameterSet_Position::SetPrecision()"); m_iPrec = i; }
	int GetPrecision() const { return m_iPrec; }

	void SetOrder(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Position::SetOrder()"); m_iOrder = i; }
	int GetOrder() const { return m_iOrder; }

	void SetSplit(int i) { BQBChkRngInt(m_IF,i,1,31,"CBQBParameterSet_Position::SetSplit()"); m_iSplit = i; }
	int GetSplit() const { return m_iSplit; }

	void SetExtraTRange(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Position::SetExtraTRange()"); m_iExtraTRange = i; }
	int GetExtraTRange() const { return m_iExtraTRange; }

	void SetExtraTOrder(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Position::SetExtraTOrder()"); m_iExtraTOrder = i; }
	int GetExtraTOrder() const { return m_iExtraTOrder; }

	void SetExtraTimeExpo(double f) { BQBChkRngFlt(m_IF,f,0,16.0,"CBQBParameterSet_Position::SetExtraTimeExpo()"); m_iExtraTimeExpo = (int)(f*1000.0+0.5); }
	double GetExtraTimeExpo() const { return m_iExtraTimeExpo/1000.0; }


	int GetVersion() const { return m_iVersion; }

	CBQBParameterSet_Compressor* GetCompressorParameterSet() { return &m_oCompressor; }

private:
	CBQBInterface &m_IF;
	CBQBParameterSet_Compressor m_oCompressor;

	int m_iPrec;           //  3 Bits, 0 .. 7
	int m_iOrder;          //  4 Bits, 0 .. 15
	bool m_bOptOrder;      //  1 Bit
	int m_iSplit;          //  5 Bits, 1 .. 31
	bool m_bSortAtom;      //  1 Bit
	bool m_bUseExtra;      //  1 Bit
	int m_iExtraTRange;    //  4 Bits, 0 .. 15
	int m_iExtraTOrder;    //  4 Bits, 0 .. 15
	int m_iExtraTimeExpo;  // 14 Bits, 0.000 .. 16.000, / 1000.0

	int m_iVersion;        //  4 Bits, 0 .. 15
};



// Current Version: 0
class CBQBParameterSet_Volumetric {
friend class CBQBParameterSet_PosAndVol;
public:

	// No pointers as members --> standard copy constructor and assignment operator will work

	explicit CBQBParameterSet_Volumetric(CBQBInterface &i) : m_IF(i), m_oCompressor(i) {

		SetSigni(5);
		SetEps(12);
		SetOrder(8);
		SetOptOrder(true);
		SetHilbert(true);
		SetNbhFac(1.075);
		SetSplit(10);
		SetTableCount(6);
		SetOptTables(false);
		SetBlockLength(20);
		SetBW(false);
		SetMTF(false);
		SetMaxIter(10);
		SetRLE(true);
		SetMaxChunk(0);
		SetUseExtra(true);
		SetExtraSRange(7);
		SetExtraTRange(6);
		SetExtraSOrder(3);
		SetExtraTOrder(2);
		SetExtraOffset(3);
		SetExtraCrossS(true);
		SetExtraCrossT(false);
		SetExtraWrap(true);
		SetExtraCrossRangeS(true);
		SetExtraCrossRangeT(false);
		SetExtraDistExpo(3.0);
		SetExtraTimeExpo(1.0);
		SetExtraPredCorr(true);
	}

	void ExportBits(CBQBBitSet *bs) const;
	bool ImportBits(CBQBBitSet *bs, int version);

	std::string ToString(int indent) const;

	std::string ToKey() const;

	bool FromKey(std::string key);

	void SetOptTables(bool b) { m_oCompressor.m_bOptTables = b; }
	bool GetOptTables() const { return m_oCompressor.m_bOptTables; }

	void SetRLE(bool b) { m_oCompressor.m_bRLE = b; }
	bool GetRLE() const { return m_oCompressor.m_bRLE; }

	void SetBW(bool b) { m_oCompressor.m_bBW = b; }
	bool GetBW() const { return m_oCompressor.m_bBW; }

	void SetMTF(bool b) { m_oCompressor.m_bMTF = b; }
	bool GetMTF() const { return m_oCompressor.m_bMTF; }

	void SetTableCount(int i) { BQBChkRngInt(m_IF,i,1,64,"CBQBParameterSet_Volumetric::SetTableCount()"); m_oCompressor.m_iTableCount = i-1; }
	int GetTableCount() const { return m_oCompressor.m_iTableCount+1; }

	void SetBlockLength(int i) { BQBChkRngInt(m_IF,i,1,128,"CBQBParameterSet_Volumetric::SetBlockLength()"); m_oCompressor.m_iBlock = i-1; }
	int GetBlockLength() const { return m_oCompressor.m_iBlock+1; }

	void SetMaxIter(int i) { BQBChkRngInt(m_IF,i,0,126,"CBQBParameterSet_Volumetric::SetMaxIter()"); m_oCompressor.m_iMaxIter = i/2; }
	int GetMaxIter() const { return m_oCompressor.m_iMaxIter*2; }

	void SetMaxChunk(int i) { BQBChkRngInt(m_IF,i,0,16777216,"CBQBParameterSet_Volumetric::SetMaxChunk()"); if (i == 0) m_oCompressor.m_iMaxChunk = 0; else m_oCompressor.m_iMaxChunk = log2i(i/1024)+1; }
	int GetMaxChunk() const { if (m_oCompressor.m_iMaxChunk == 0) return 0; else return pow2i(m_oCompressor.m_iMaxChunk-1)*1024; }



	void SetOptOrder(bool b) { m_bOptOrder = b; }
	bool GetOptOrder() const { return m_bOptOrder; }

	void SetHilbert(bool b) { m_bHilbert = b; }
	bool GetHilbert() const { return m_bHilbert; }

	void SetUseExtra(bool b) { m_bUseExtra = b; }
	bool GetUseExtra() const { return m_bUseExtra; }

	void SetExtraCrossS(bool b) { m_bExtraCrossS = b; }
	bool GetExtraCrossS() const { return m_bExtraCrossS; }

	void SetExtraCrossT(bool b) { m_bExtraCrossT = b; }
	bool GetExtraCrossT() const { return m_bExtraCrossT; }

	void SetExtraWrap(bool b) { m_bExtraWrap = b; }
	bool GetExtraWrap() const { return m_bExtraWrap; }

	void SetExtraCrossRangeS(bool b) { m_bExtraCrossRangeS = b; }
	bool GetExtraCrossRangeS() const { return m_bExtraCrossRangeS; }

	void SetExtraCrossRangeT(bool b) { m_bExtraCrossRangeT = b; }
	bool GetExtraCrossRangeT() const { return m_bExtraCrossRangeT; }

	void SetExtraPredCorr(bool b) { m_bExtraPredCorr = b; }
	bool GetExtraPredCorr() const { return m_bExtraPredCorr; }



	void SetSigni(int i) { BQBChkRngInt(m_IF,i,1,9,"CBQBParameterSet_Volumetric::SetSigni()"); m_iSigni = i; }
	int GetSigni() const { return m_iSigni; }

	void SetOrder(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetOrder()"); m_iOrder = i; }
	int GetOrder() const { return m_iOrder; }

	void SetEps(int i) { BQBChkRngInt(m_IF,i,0,63,"CBQBParameterSet_Volumetric::SetEps()"); m_iEps = i; }
	int GetEps() const { return m_iEps; }

	void SetSplit(int i) { BQBChkRngInt(m_IF,i,1,31,"CBQBParameterSet_Volumetric::SetSplit()"); m_iSplit = i; }
	int GetSplit() const { return m_iSplit; }

	void SetExtraSRangeX(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetExtraSRangeX()"); m_iExtraSRangeX = i; }
	int GetExtraSRangeX() const { return m_iExtraSRangeX; }

	void SetExtraSRangeY(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetExtraSRangeY()"); m_iExtraSRangeY = i; }
	int GetExtraSRangeY() const { return m_iExtraSRangeY; }

	void SetExtraSRangeZ(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetExtraSRangeZ()"); m_iExtraSRangeZ = i; }
	int GetExtraSRangeZ() const { return m_iExtraSRangeZ; }

	void SetExtraSRange(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetExtraSRange()"); m_iExtraSRangeX = i; m_iExtraSRangeY = i; m_iExtraSRangeZ = i; }
	bool IsExtraSRangeXYZEqual() const { if ((m_iExtraSRangeX == m_iExtraSRangeY) && (m_iExtraSRangeX == m_iExtraSRangeZ)) return true; else return false; }

	void SetExtraTRange(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetExtraTRange()"); m_iExtraTRange = i; }
	int GetExtraTRange() const { return m_iExtraTRange; }

	void SetExtraSOrder(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetExtraSOrder()"); m_iExtraSOrder = i; }
	int GetExtraSOrder() const { return m_iExtraSOrder; }

	void SetExtraTOrder(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetExtraTOrder()"); m_iExtraTOrder = i; }
	int GetExtraTOrder() const { return m_iExtraTOrder; }

	void SetExtraOffsetX(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetExtraOffsetX()"); m_iExtraOffsetX = i; }
	int GetExtraOffsetX() const { return m_iExtraOffsetX; }

	void SetExtraOffsetY(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetExtraOffsetY()"); m_iExtraOffsetY = i; }
	int GetExtraOffsetY() const { return m_iExtraOffsetY; }

	void SetExtraOffsetZ(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetExtraOffsetZ()"); m_iExtraOffsetZ = i; }
	int GetExtraOffsetZ() const { return m_iExtraOffsetZ; }

	void SetExtraOffset(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_Volumetric::SetExtraOffset()"); m_iExtraOffsetX = i; m_iExtraOffsetY = i; m_iExtraOffsetZ = i; }
	bool IsExtraOffsetXYZEqual() const { if ((m_iExtraOffsetX == m_iExtraOffsetY) && (m_iExtraOffsetX == m_iExtraOffsetZ)) return true; else return false; }

	void SetNbhFac(double f) { BQBChkRngFlt(m_IF,f,0,4.0,"CBQBParameterSet_Volumetric::SetNbhFac()"); m_iNbhFac = (int)(f*1000.0+0.5); }
	double GetNbhFac() const { return m_iNbhFac/1000.0; }

	void SetExtraDistExpo(double f) { BQBChkRngFlt(m_IF,f,0,16.0,"CBQBParameterSet_Volumetric::SetExtraDistExpo()"); m_iExtraDistExpo = (int)(f*1000.0+0.5); }
	double GetExtraDistExpo() const { return m_iExtraDistExpo/1000.0; }

	void SetExtraTimeExpo(double f) { BQBChkRngFlt(m_IF,f,0,16.0,"CBQBParameterSet_Volumetric::SetExtraTimeExpo()"); m_iExtraTimeExpo = (int)(f*1000.0+0.5); }
	double GetExtraTimeExpo() const { return m_iExtraTimeExpo/1000.0; }


	int GetVersion() const { return m_iVersion; }

	CBQBParameterSet_Compressor* GetCompressorParameterSet() { return &m_oCompressor; }

private:
	CBQBInterface &m_IF;
	CBQBParameterSet_Compressor m_oCompressor;

	int m_iSigni;              //  4 Bits, 1 .. 9
	int m_iOrder;              //  4 Bits, 0 .. 15
	bool m_bOptOrder;          //  1 Bit
	int m_iEps;                //  6 Bits, 0 .. 63
	int m_iSplit;              //  5 Bits, 1 .. 31
	bool m_bHilbert;           //  1 Bit
	int m_iNbhFac;             // 12 Bits, 0.000 .. 4.000, / 1000.0
	bool m_bUseExtra;          //  1 Bit
	int m_iExtraSRangeX;       //  4 Bits, 0 .. 15
	int m_iExtraSRangeY;       //  4 Bits, 0 .. 15
	int m_iExtraSRangeZ;       //  4 Bits, 0 .. 15
	int m_iExtraTRange;        //  4 Bits, 0 .. 15
	int m_iExtraSOrder;        //  4 Bits, 0 .. 15
	int m_iExtraTOrder;        //  4 Bits, 0 .. 15
	int m_iExtraOffsetX;       //  4 Bits, 0 .. 15
	int m_iExtraOffsetY;       //  4 Bits, 0 .. 15
	int m_iExtraOffsetZ;       //  4 Bits, 0 .. 15
	bool m_bExtraCrossS;       //  1 Bit
	bool m_bExtraCrossT;       //  1 Bit
	bool m_bExtraWrap;         //  1 Bit
	bool m_bExtraCrossRangeS;  //  1 Bit
	bool m_bExtraCrossRangeT;  //  1 Bit
	int m_iExtraDistExpo;      // 14 Bits, 0.000 .. 16.000, / 1000.0
	int m_iExtraTimeExpo;      // 14 Bits, 0.000 .. 16.000, / 1000.0
	bool m_bExtraPredCorr;     //  1 Bit

	int m_iVersion;            //  4 Bits, 0 .. 15
};



// Current Version: 0
class CBQBParameterSet_PosAndVol {
public:

	// No pointers as members --> standard copy constructor and assignment operator will work

	explicit CBQBParameterSet_PosAndVol(CBQBInterface &i) : m_IF(i), m_oPosition(i), m_oVolumetric(i), m_iVersion(0) {
	}

	void ExportBits(CBQBBitSet *bs) const;
	bool ImportBits(CBQBBitSet *bs, int version);

	std::string ToString(int indent) const;

	std::string ToKey() const;

	bool FromKey(std::string key);

	void SetVolOptTables(bool b) { m_oVolumetric.m_oCompressor.m_bOptTables = b; }
	bool GetVolOptTables() const { return m_oVolumetric.m_oCompressor.m_bOptTables; }

	void SetVolRLE(bool b) { m_oVolumetric.m_oCompressor.m_bRLE = b; }
	bool GetVolRLE() const { return m_oVolumetric.m_oCompressor.m_bRLE; }

	void SetVolBW(bool b) { m_oVolumetric.m_oCompressor.m_bBW = b; }
	bool GetVolBW() const { return m_oVolumetric.m_oCompressor.m_bBW; }

	void SetVolMTF(bool b) { m_oVolumetric.m_oCompressor.m_bMTF = b; }
	bool GetVolMTF() const { return m_oVolumetric.m_oCompressor.m_bMTF; }

	void SetVolTableCount(int i) { BQBChkRngInt(m_IF,i,1,64,"CBQBParameterSet_PosAndVol::SetVolTableCount()"); m_oVolumetric.m_oCompressor.m_iTableCount = i-1; }
	int GetVolTableCount() const { return m_oVolumetric.m_oCompressor.m_iTableCount+1; }

	void SetVolBlockLength(int i) { BQBChkRngInt(m_IF,i,1,128,"CBQBParameterSet_PosAndVol::SetVolBlockLength()"); m_oVolumetric.m_oCompressor.m_iBlock = i-1; }
	int GetVolBlockLength() const { return m_oVolumetric.m_oCompressor.m_iBlock+1; }

	void SetVolMaxIter(int i) { BQBChkRngInt(m_IF,i,0,126,"CBQBParameterSet_PosAndVol::SetVolMaxIter()"); m_oVolumetric.m_oCompressor.m_iMaxIter = i/2; }
	int GetVolMaxIter() const { return m_oVolumetric.m_oCompressor.m_iMaxIter*2; }

	void SetVolMaxChunk(int i) { BQBChkRngInt(m_IF,i,0,16777216,"CBQBParameterSet_PosAndVol::SetVolMaxChunk()"); if (i == 0) m_oVolumetric.m_oCompressor.m_iMaxChunk = 0; else m_oVolumetric.m_oCompressor.m_iMaxChunk = log2i(i/1024)+1; }
	int GetVolMaxChunk() const { if (m_oVolumetric.m_oCompressor.m_iMaxChunk == 0) return 0; else return pow2i(m_oVolumetric.m_oCompressor.m_iMaxChunk-1)*1024; }



	void SetVolOptOrder(bool b) { m_oVolumetric.m_bOptOrder = b; }
	bool GetVolOptOrder() const { return m_oVolumetric.m_bOptOrder; }

	void SetVolHilbert(bool b) { m_oVolumetric.m_bHilbert = b; }
	bool GetVolHilbert() const { return m_oVolumetric.m_bHilbert; }

	void SetVolUseExtra(bool b) { m_oVolumetric.m_bUseExtra = b; }
	bool GetVolUseExtra() const { return m_oVolumetric.m_bUseExtra; }

	void SetVolExtraCrossS(bool b) { m_oVolumetric.m_bExtraCrossS = b; }
	bool GetVolExtraCrossS() const { return m_oVolumetric.m_bExtraCrossS; }

	void SetVolExtraCrossT(bool b) { m_oVolumetric.m_bExtraCrossT = b; }
	bool GetVolExtraCrossT() const { return m_oVolumetric.m_bExtraCrossT; }

	void SetVolExtraWrap(bool b) { m_oVolumetric.m_bExtraWrap = b; }
	bool GetVolExtraWrap() const { return m_oVolumetric.m_bExtraWrap; }

	void SetVolExtraCrossRangeS(bool b) { m_oVolumetric.m_bExtraCrossRangeS = b; }
	bool GetVolExtraCrossRangeS() const { return m_oVolumetric.m_bExtraCrossRangeS; }

	void SetVolExtraCrossRangeT(bool b) { m_oVolumetric.m_bExtraCrossRangeT = b; }
	bool GetVolExtraCrossRangeT() const { return m_oVolumetric.m_bExtraCrossRangeT; }

	void SetVolExtraPredCorr(bool b) { m_oVolumetric.m_bExtraPredCorr = b; }
	bool GetVolExtraPredCorr() const { return m_oVolumetric.m_bExtraPredCorr; }



	void SetVolSigni(int i) { BQBChkRngInt(m_IF,i,1,9,"CBQBParameterSet_PosAndVol::SetVolSigni()"); m_oVolumetric.m_iSigni = i; }
	int GetVolSigni() const { return m_oVolumetric.m_iSigni; }

	void SetVolOrder(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolOrder()"); m_oVolumetric.m_iOrder = i; }
	int GetVolOrder() const { return m_oVolumetric.m_iOrder; }

	void SetVolEps(int i) { BQBChkRngInt(m_IF,i,0,63,"CBQBParameterSet_PosAndVol::SetVolEps()"); m_oVolumetric.m_iEps = i; }
	int GetVolEps() const { return m_oVolumetric.m_iEps; }

	void SetVolSplit(int i) { BQBChkRngInt(m_IF,i,1,31,"CBQBParameterSet_PosAndVol::SetVolSplit()"); m_oVolumetric.m_iSplit = i; }
	int GetVolSplit() const { return m_oVolumetric.m_iSplit; }

	void SetVolExtraSRangeX(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolExtraSRangeX()"); m_oVolumetric.m_iExtraSRangeX = i; }
	int GetVolExtraSRangeX() const { return m_oVolumetric.m_iExtraSRangeX; }

	void SetVolExtraSRangeY(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolExtraSRangeY()"); m_oVolumetric.m_iExtraSRangeY = i; }
	int GetVolExtraSRangeY() const { return m_oVolumetric.m_iExtraSRangeY; }

	void SetVolExtraSRangeZ(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolExtraSRangeZ()"); m_oVolumetric.m_iExtraSRangeZ = i; }
	int GetVolExtraSRangeZ() const { return m_oVolumetric.m_iExtraSRangeZ; }

	void SetVolExtraSRange(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolExtraSRange()"); m_oVolumetric.m_iExtraSRangeX = i; m_oVolumetric.m_iExtraSRangeY = i; m_oVolumetric.m_iExtraSRangeZ = i; }
	bool IsVolExtraSRangeXYZEqual() const { if ((m_oVolumetric.m_iExtraSRangeX == m_oVolumetric.m_iExtraSRangeY) && (m_oVolumetric.m_iExtraSRangeX == m_oVolumetric.m_iExtraSRangeZ)) return true; else return false; }

	void SetVolExtraTRange(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolExtraTRange()"); m_oVolumetric.m_iExtraTRange = i; }
	int GetVolExtraTRange() const { return m_oVolumetric.m_iExtraTRange; }

	void SetVolExtraSOrder(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolExtraSOrder()"); m_oVolumetric.m_iExtraSOrder = i; }
	int GetVolExtraSOrder() const { return m_oVolumetric.m_iExtraSOrder; }

	void SetVolExtraTOrder(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolExtraTOrder()"); m_oVolumetric.m_iExtraTOrder = i; }
	int GetVolExtraTOrder() const { return m_oVolumetric.m_iExtraTOrder; }

	void SetVolExtraOffsetX(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolExtraOffsetX()"); m_oVolumetric.m_iExtraOffsetX = i; }
	int GetVolExtraOffsetX() const { return m_oVolumetric.m_iExtraOffsetX; }

	void SetVolExtraOffsetY(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolExtraOffsetY()"); m_oVolumetric.m_iExtraOffsetY = i; }
	int GetVolExtraOffsetY() const { return m_oVolumetric.m_iExtraOffsetY; }

	void SetVolExtraOffsetZ(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolExtraOffsetZ()"); m_oVolumetric.m_iExtraOffsetZ = i; }
	int GetVolExtraOffsetZ() const { return m_oVolumetric.m_iExtraOffsetZ; }

	void SetVolExtraOffset(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetVolExtraOffset()"); m_oVolumetric.m_iExtraOffsetX = i; m_oVolumetric.m_iExtraOffsetY = i; m_oVolumetric.m_iExtraOffsetZ = i; }
	bool IsVolExtraOffsetXYZEqual() const { if ((m_oVolumetric.m_iExtraOffsetX == m_oVolumetric.m_iExtraOffsetY) && (m_oVolumetric.m_iExtraOffsetX == m_oVolumetric.m_iExtraOffsetZ)) return true; else return false; }

	void SetVolNbhFac(double f) { BQBChkRngFlt(m_IF,f,0,4.0,"CBQBParameterSet_PosAndVol::SetVolNbhFac()"); m_oVolumetric.m_iNbhFac = (int)(f*1000.0+0.5); }
	double GetVolNbhFac() const { return m_oVolumetric.m_iNbhFac/1000.0; }

	void SetVolExtraDistExpo(double f) { BQBChkRngFlt(m_IF,f,0,16.0,"CBQBParameterSet_PosAndVol::SetVolExtraDistExpo()"); m_oVolumetric.m_iExtraDistExpo = (int)(f*1000.0+0.5); }
	double GetVolExtraDistExpo() const { return m_oVolumetric.m_iExtraDistExpo/1000.0; }

	void SetVolExtraTimeExpo(double f) { BQBChkRngFlt(m_IF,f,0,16.0,"CBQBParameterSet_PosAndVol::SetVolExtraTimeExpo()"); m_oVolumetric.m_iExtraTimeExpo = (int)(f*1000.0+0.5); }
	double GetVolExtraTimeExpo() const { return m_oVolumetric.m_iExtraTimeExpo/1000.0; }





	void SetPosOptTables(bool b) { m_oPosition.m_oCompressor.m_bOptTables = b; }
	bool GetPosOptTables() const { return m_oPosition.m_oCompressor.m_bOptTables; }

	void SetPosRLE(bool b) { m_oPosition.m_oCompressor.m_bRLE = b; }
	bool GetPosRLE() const { return m_oPosition.m_oCompressor.m_bRLE; }

	void SetPosBW(bool b) { m_oPosition.m_oCompressor.m_bBW = b; }
	bool GetPosBW() const { return m_oPosition.m_oCompressor.m_bBW; }

	void SetPosMTF(bool b) { m_oPosition.m_oCompressor.m_bMTF = b; }
	bool GetPosMTF() const { return m_oPosition.m_oCompressor.m_bMTF; }

	void SetPosTableCount(int i) { BQBChkRngInt(m_IF,i,1,64,"CBQBParameterSet_PosAndVol::SetPosTableCount()"); m_oPosition.m_oCompressor.m_iTableCount = i-1; }
	int GetPosTableCount() const { return m_oPosition.m_oCompressor.m_iTableCount+1; }

	void SetPosBlockLength(int i) { BQBChkRngInt(m_IF,i,1,128,"CBQBParameterSet_PosAndVol::SetPosBlockLength()"); m_oPosition.m_oCompressor.m_iBlock = i-1; }
	int GetPosBlockLength() const { return m_oPosition.m_oCompressor.m_iBlock+1; }

	void SetPosMaxIter(int i) { BQBChkRngInt(m_IF,i,0,126,"CBQBParameterSet_PosAndVol::SetPosMaxIter()"); m_oPosition.m_oCompressor.m_iMaxIter = i/2; }
	int GetPosMaxIter() const { return m_oPosition.m_oCompressor.m_iMaxIter*2; }

	void SetPosMaxChunk(int i) { BQBChkRngInt(m_IF,i,0,16777216,"CBQBParameterSet_PosAndVol::SetPosMaxChunk()"); if (i == 0) m_oPosition.m_oCompressor.m_iMaxChunk = 0; else m_oPosition.m_oCompressor.m_iMaxChunk = log2i(i/1024)+1; }
	int GetPosMaxChunk() const { if (m_oPosition.m_oCompressor.m_iMaxChunk == 0) return 0; else return pow2i(m_oPosition.m_oCompressor.m_iMaxChunk-1)*1024; }



	void SetPosOptOrder(bool b) { m_oPosition.m_bOptOrder = b; }
	bool GetPosOptOrder() const { return m_oPosition.m_bOptOrder; }

	void SetPosSortAtom(bool b) { m_oPosition.m_bSortAtom = b; }
	bool GetPosSortAtom() const { return m_oPosition.m_bSortAtom; }

	void SetPosUseExtra(bool b) { m_oPosition.m_bUseExtra = b; }
	bool GetPosUseExtra() const { return m_oPosition.m_bUseExtra; }

	void SetPosPrecision(int i) { BQBChkRngInt(m_IF,i,0,7,"CBQBParameterSet_PosAndVol::SetPosPrecision()"); m_oPosition.m_iPrec = i; }
	int GetPosPrecision() const { return m_oPosition.m_iPrec; }

	void SetPosOrder(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetPosOrder()"); m_oPosition.m_iOrder = i; }
	int GetPosOrder() const { return m_oPosition.m_iOrder; }

	void SetPosSplit(int i) { BQBChkRngInt(m_IF,i,1,31,"CBQBParameterSet_PosAndVol::SetPosSplit()"); m_oPosition.m_iSplit = i; }
	int GetPosSplit() const { return m_oPosition.m_iSplit; }

	void SetPosExtraTRange(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetPosExtraTRange()"); m_oPosition.m_iExtraTRange = i; }
	int GetPosExtraTRange() const { return m_oPosition.m_iExtraTRange; }

	void SetPosExtraTOrder(int i) { BQBChkRngInt(m_IF,i,0,15,"CBQBParameterSet_PosAndVol::SetPosExtraTOrder()"); m_oPosition.m_iExtraTOrder = i; }
	int GetPosExtraTOrder() const { return m_oPosition.m_iExtraTOrder; }

	void SetPosExtraTimeExpo(double f) { BQBChkRngFlt(m_IF,f,0,16.0,"CBQBParameterSet_PosAndVol::SetPosExtraTimeExpo()"); m_oPosition.m_iExtraTimeExpo = (int)(f*1000.0+0.5); }
	double GetPosExtraTimeExpo() const { return m_oPosition.m_iExtraTimeExpo/1000.0; }


	int GetVersion() const { return m_iVersion; }

	CBQBParameterSet_Position* GetPositionParameterSet() { return &m_oPosition; }

	CBQBParameterSet_Volumetric* GetVolumetricParameterSet() { return &m_oVolumetric; }

	CBQBParameterSet_Compressor* GetVolumetricCompressorParameterSet() { return &m_oVolumetric.m_oCompressor; }

	CBQBParameterSet_Compressor* GetPositionCompressorParameterSet() { return &m_oPosition.m_oCompressor; }

private:
	CBQBInterface &m_IF;
	CBQBParameterSet_Position m_oPosition;
	CBQBParameterSet_Volumetric m_oVolumetric;
	int m_iVersion;            //  4 Bits, 0 .. 15
};


#endif


