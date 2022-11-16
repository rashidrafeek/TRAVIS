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



#ifndef BQB_FORMAT_H
#define BQB_FORMAT_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_math.h"
#include "bqb_parmset.h"

#include <string>
#include <vector>



class CBQBInterface;


/* Frame Types:

	 0 - Reserved
	 1 - General Container Frame
	 2 - Index Frame
	 3 - Compressed Index Frame
	 4 - Trajectory Frame
	 5 - Compressed Trajectory Frame
	 6 - Volumetric Data Frame
	 7 - Compressed Volumetric Data Frame

*/

#define BQB_FRAMETYPE_GENERAL         (1)
#define BQB_FRAMETYPE_IDX             (2)
#define BQB_FRAMETYPE_COMPIDX         (3)
#define BQB_FRAMETYPE_TRAJ            (4)
#define BQB_FRAMETYPE_COMPTRAJSTART   (5)
#define BQB_FRAMETYPE_COMPTRAJ        (6)
#define BQB_FRAMETYPE_CUBE            (7)
#define BQB_FRAMETYPE_COMPCUBESTART   (8)
#define BQB_FRAMETYPE_COMPCUBE        (9)
#define BQB_FRAMETYPE_FILE           (10)
#define BQB_FRAMETYPE_COMPFILE       (11)
#define BQB_FRAMETYPE_COMPTRAJKEY    (12)
#define BQB_FRAMETYPE_COMPCUBEKEY    (13)


#define BQB_TYPE_STRING   (1)
#define BQB_TYPE_FLOAT    (2)
#define BQB_TYPE_DOUBLE   (3)
#define BQB_TYPE_INT8     (4)
#define BQB_TYPE_INT16    (5)
#define BQB_TYPE_INT32    (6)
#define BQB_TYPE_UINT8    (7)
#define BQB_TYPE_UINT16   (8)
#define BQB_TYPE_UINT32   (9)


#define BQB_FRAMETYPE_COMPCUBE_VERSION  (2)   // Currently Version 2
#define BQB_FRAMETYPE_COMPTRAJ_VERSION  (2)   // Currently Version 2
#define BQB_FRAMETYPE_COMPFILE_VERSION  (0)   // Currently Version 0
#define BQB_FRAMETYPE_COMPIDX_VERSION   (1)   // Currently Version 1



//extern bool g_bBQBVerbose;


class CBQBEngine;


//void SetBQBVerbose(bool v);

bool WriteArrayToFile(FILE *a, const std::vector<unsigned char> &ia);

bool ReadArrayFromFile(FILE *a, std::vector<unsigned char> &ia, int length);

void ReadArrayFromFile(FILE *a, std::vector<unsigned char> &ia);


class CBQBFile;
class CBQBListEntry;

class CBQBBitSet;



const char *GetFrameTypeString(int type);




class CBQBIndex {
public:

	explicit CBQBIndex(CBQBInterface &i) : m_IF(i) {
	}

	CBQBIndex(const CBQBIndex &bi)
		: m_iaFrameLengths(bi.m_iaFrameLengths), m_iaFrameTypes(bi.m_iaFrameTypes), m_IF(bi.m_IF) {
	}

	~CBQBIndex() {
	}

	bool ImportFromArray(bool compressed, const std::vector<unsigned char> &ia, int ver);
	bool ExportToArray(bool compressed, std::vector<unsigned char> &ia, CBQBStatistics *stat);

	void Dump();

	std::vector<int> m_iaFrameLengths;
	std::vector<int> m_iaFrameTypes;
	std::vector<int> m_iaFrameIDs;

private:
	CBQBInterface &m_IF;
};



class CBQBShortFrame {
public:

	CBQBShortFrame()
		: m_iFrameType(0), m_iFrameTypeVersion(0), m_iID(0), m_iCRC32(0) {
	}

	~CBQBShortFrame() {
	}

	int m_iFrameType;
	int m_iFrameTypeVersion;
	int m_iID;
	unsigned long m_iCRC32;

	std::vector<unsigned char> m_iaPayload;

};



class CBQBListEntry {
public:

	CBQBListEntry()
		: m_iFrameCount(-1), m_iFullFrameCount(-1), m_pIndex(NULL), m_pFile(NULL), m_iFrameStart(-1), m_iFrameEnd(-1) {
	}

	~CBQBListEntry();

	std::string m_sFileName;
	int m_iFrameCount;
	int m_iFullFrameCount;
	CBQBIndex *m_pIndex;
	CBQBFile *m_pFile;
	int m_iFrameStart;
	int m_iFrameEnd;
};



class CBQBFile {
public:

	explicit CBQBFile(CBQBInterface &i)
		: m_oParmSetCube(i), m_oParmSetPos(i), m_pEngine(NULL), m_IF(i), m_bEOF(false),
		m_bShortFrame(false), m_bListFile(false), m_pFile(NULL), m_bOpenRead(false),
		m_bOpenWrite(false), m_pShortFrame(NULL), m_pIndex(NULL), m_iTotalFrameCount(-1) {

		#ifndef _MSC_VER
			m_bIndexWritten = false;
			m_iIndexOffset = 0;
		#endif
	}

	~CBQBFile();

	bool OpenRead(std::string s);
	bool OpenWriteAppend(std::string s);
	bool OpenWriteReplace(std::string s);
	bool Close();
	bool Rewind();

	bool SeekFrame(int i);
	bool SkipFrames(int i);

	bool IsEOF();

	bool OpenListFile(FILE *a);

	// General Frame Functions
	bool CreateGeneralFrame();

	// Short Frame Functions
	bool CreateShortFrame(int type, int version, int id);
	bool PushPayload(const std::vector<unsigned char> &ia);

	bool FinalizeFrame(CBQBStatistics *stat);

	bool WriteIndexFrame(bool compressed, CBQBStatistics *stat);

	int GetFrameID() const;
	int GetFrameType() const;
	int GetFrameTypeVersion() const;
	const std::vector<unsigned char>* GetFramePayload() const;

	void GetFileOffset(unsigned long &high, unsigned long &low);

	bool ReadFrame();

	void DumpIndex();

	bool CheckIntegrity(std::string s, bool verbose);

	bool CompareCoords(std::string infile, std::string reffile, bool verbose);

	int GetTotalFrameCount() const {
		return m_iTotalFrameCount;
	}

	CBQBEngine* GetEngine();



	CBQBParameterSet_PosAndVol m_oParmSetCube;
	CBQBParameterSet_Position m_oParmSetPos;

private:
	CBQBEngine *m_pEngine;
	CBQBInterface &m_IF;
	int m_iFrameCounter;
	bool m_bEOF;
	int m_iListIndex;
	std::vector<CBQBListEntry*> m_oaBQBList;
	bool m_bShortFrame;
	bool m_bListFile;
	std::string m_sFileName;
	std::string m_sDirectory;
	FILE *m_pFile;
	bool m_bOpenRead;
	bool m_bOpenWrite;
	CBQBShortFrame *m_pShortFrame;
	CBQBIndex *m_pIndex;
	int m_iTotalFrameCount;

	#ifndef _MSC_VER
		bool m_bIndexWritten;
		off_t m_iIndexOffset;
	#endif
};



class CBQBTrajectoryFrameColumn {
public:

	explicit CBQBTrajectoryFrameColumn(CBQBInterface &i) : m_iType(255), m_IF(i) {
	}

	CBQBTrajectoryFrameColumn(CBQBInterface &i, unsigned char type, std::string label) : 
		m_iType(type), m_sLabel(label), m_IF(i) {
	}

	~CBQBTrajectoryFrameColumn() {
	}

	bool ReadColumn(int ac, CBQBBitSet *bs);
	void WriteColumn(int ac, CBQBBitSet *bs);

	unsigned char m_iType;
	std::string m_sLabel;

	std::vector<std::string> m_aString;
	std::vector<double> m_aReal;
	std::vector<long> m_aSignedInt;
	std::vector<unsigned long> m_aUnsignedInt;

private:
	CBQBInterface &m_IF;
};



class CBQBTrajectoryFrame {
public:

	explicit CBQBTrajectoryFrame(CBQBInterface &i) : m_iAtomCount(0), m_pCellMatrix(NULL), m_IF(i) {
	}

	~CBQBTrajectoryFrame() {
		int z;
		for (z=0;z<(int)m_oaColumns.size();z++)
			delete m_oaColumns[z];
		if (m_pCellMatrix != NULL) {
			delete m_pCellMatrix;
			m_pCellMatrix = NULL;
		}
	}

	bool ReadFrame(const std::vector<unsigned char> *data);
	void WriteFrame(std::vector<unsigned char> *data);

	CBQBTrajectoryFrameColumn* GetColumn(std::string label);
	CBQBTrajectoryFrameColumn* AddColumn(int type, std::string label);


	unsigned long m_iAtomCount;
	std::vector<CBQBTrajectoryFrameColumn*> m_oaColumns;
	CBQBDMatrix3 *m_pCellMatrix;

private:
	CBQBInterface &m_IF;
};



#endif

