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

#include "bqb_format.h"
#include "bqb_crc.h"
#include "bqb_integerengine.h"
#include "bqb_engine.h"


#ifdef _WIN32
	#include <io.h>
#endif


const char *GetRevisionInfo_bqb_format(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqb_format() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}




CBQBListEntry::~CBQBListEntry() {

	if (m_pIndex != NULL) {
		delete m_pIndex;
		m_pIndex = NULL;
	}
	if (m_pFile != NULL) {
		delete m_pFile;
		m_pFile = NULL;
	}
}


CBQBFile::~CBQBFile() {

	int z;

	for (z=0;z<(int)m_oaBQBList.size();z++)
		delete m_oaBQBList[z];

	if (m_pShortFrame != NULL) {
		if (!m_bListFile) // Otherwise it is deleted elsewhere
			delete m_pShortFrame;
		m_pShortFrame = NULL;
	}
	if (m_pIndex != NULL) {
		delete m_pIndex;
		m_pIndex = NULL;
	}
	if (m_pEngine != NULL) {
		delete m_pEngine;
		m_pEngine = NULL;
	}
}


bool WriteArrayToFile(FILE *a, const std::vector<unsigned char> &ia) {

	unsigned int i;

	i = 0;
	while (i < ia.size())
		i += (unsigned int)fwrite(&ia[i],1,(i+4096<=ia.size())?4096:(ia.size()-i),a);

	return true;
}


bool ReadArrayFromFile(FILE *a, std::vector<unsigned char> &ia, int length) {

	unsigned int i, k;
	int z;
	static unsigned char buf[4096];

	i = (unsigned int)ia.size();
	ia.resize(ia.size()+length);
	while (i < ia.size()) {
		k = (unsigned int)fread(buf,1,(i+4096<=ia.size())?4096:(ia.size()-i),a);
		for (z=0;z<(int)k;z++)
			ia[i+z] = buf[z];
		i += k;
		if (feof(a))
			return false;
	}

	return true;
}


void ReadArrayFromFile(FILE *a, std::vector<unsigned char> &ia) {

	unsigned int k;
	int z;
	static unsigned char buf[4096];

	while (!feof(a)) {
		k = (unsigned int)fread(buf,1,4096,a);
		for (z=0;z<(int)k;z++)
			ia.push_back(buf[z]);
	}
}


bool CBQBIndex::ImportFromArray(bool compressed, const std::vector<unsigned char> &ia, int ver) {

	CBQBIntegerEngine ie(m_IF);
	CBQBBitSet bs(m_IF);
	std::vector<int> diff;
	int z;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBIndex::ImportFromArray >>>\n");

	if (compressed) {

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Importing compressed index frame...\n");

		m_iaFrameTypes.clear();
		m_iaFrameLengths.clear();

		bs.m_iaData.assign( ia.begin(), ia.end() );

		ie.DecompressSingle(
			&bs,
			m_iaFrameTypes,
			NULL
		);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Decompressed %lu frame types.\n",(unsigned long)m_iaFrameTypes.size());

	//	for (z=0;z<(int)m_iaFrameTypes.size();z++)
	//		m_IF.printf("@@@ %d\n",m_iaFrameTypes[z]);

		ie.DecompressSingle(
			&bs,
			m_iaFrameLengths,
			NULL
		);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Decompressed %lu frame lengths.\n",(unsigned long)m_iaFrameLengths.size());

		if (ver >= 1) {

			ie.DecompressSingle(
				&bs,
				diff,
				NULL
			);

			// Reconstruct values from differences
			m_iaFrameIDs.resize(diff.size());
			m_iaFrameIDs[0] = diff[0];
			for (z=1;z<(int)diff.size();z++)
				m_iaFrameIDs[z] = m_iaFrameIDs[z-1] + diff[z];

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("    Decompressed %lu frame IDs.\n",(unsigned long)m_iaFrameIDs.size());

		} else {

			m_iaFrameIDs.resize(m_iaFrameLengths.size());
			for (z=0;z<(int)m_iaFrameIDs.size();z++)
				m_iaFrameIDs[z] = 0;
		}

		for (z=0;z<(int)m_iaFrameTypes.size();z++)
			if ((m_iaFrameTypes[z] == 2) || (m_iaFrameTypes[z] == 3))
				m_IF.eprintf("CBQBIndex::ImportFromArray(): Warning: Index frame is on the index (%d/%lu).\n",z+1,(unsigned long)m_iaFrameTypes.size());

	} else {

		m_IF.eprintf("CBQBIndex::ImportFromArray(): Uncompressed index not yet implemented.\n");
		abort();

	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBIndex::ImportFromArray <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


bool CBQBIndex::ExportToArray(bool compressed, std::vector<unsigned char> &ia, CBQBStatistics *stat) {

	CBQBIntegerEngine ie(m_IF);
	CBQBBitSet bs(m_IF);
	std::vector<int> diff;
	int i, z;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBIndex::ExportToArray >>>\n");

	if (stat != NULL)
		stat->PushStatistics();

	if (compressed) {

	//	for (i=0;i<(int)m_iaFrameTypes.size();i++)
	//		m_IF.printf("@@@ %d\n",m_iaFrameTypes[i]);

		ie.CompressSingle(
			m_iaFrameTypes,
			&bs,
			false,
			true,
			true,
			50,
			1,
			false,
			false,
			10,
			stat
		);

		i = bs.GetByteLength();
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Compressed %lu frame types to %d bytes.\n",(unsigned long)m_iaFrameTypes.size(),i);

		ie.CompressSingle(
			m_iaFrameLengths,
			&bs,
			false,
			false,
			true,
			50,
			1,
			false,
			false,
			10,
			stat
		);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Compressed %lu frame lengths to %d bytes.\n",(unsigned long)m_iaFrameLengths.size(),bs.GetByteLength()-i);
		i = bs.GetByteLength();

		// Store differences instead of values
		diff.resize(m_iaFrameIDs.size());
		diff[0] = m_iaFrameIDs[0];
		for (z=1;z<(int)m_iaFrameIDs.size();z++)
			diff[z] = m_iaFrameIDs[z] - m_iaFrameIDs[z-1];

		ie.CompressSingle(
			diff,
			&bs,
			false,
			true,
			true,
			50,
			1,
			false,
			false,
			10,
			stat
		);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Compressed %lu frame IDs to %d bytes.\n",(unsigned long)m_iaFrameIDs.size(),bs.GetByteLength()-i);

		ia.insert(ia.end(),bs.m_iaData.begin(),bs.m_iaData.end());

	} else {

		m_IF.eprintf("CBQBIndex::ExportToArray(): Uncompressed index not yet implemented.\n");
		abort();
	}

	if (stat != NULL)
		stat->PopStatistics();

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBIndex::ExportToArray <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


void CBQBIndex::Dump() {

	int z;
	unsigned long u;

	m_IF.printf("\n");
	m_IF.printf("    Offset    Length    Type    Ver          ID\n");
	m_IF.printf("---------------------------------------------------\n");

	u = 0;
	for (z=0;z<(int)m_iaFrameTypes.size();z++) {

		m_IF.printf("%10lu  %8d     %3d    %3d  %10d\n",
			u,
			m_iaFrameLengths[z],
			m_iaFrameTypes[z]>>3,
			m_iaFrameTypes[z]%7,
			m_iaFrameIDs[z]
		);

		u += m_iaFrameLengths[z];
	}

	m_IF.printf("\n");
}


bool CBQBFile::OpenRead(std::string s) {

	long ip;
	unsigned char uc[5];
	CBQBTools bqbtools(m_IF);
	CBQBBitSet bs(m_IF);
	unsigned long fohi, folo, rohi, rolo;
	std::vector<unsigned char> tuca;
	std::string::size_type sp1, sp2;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::OpenRead >>>\n");

	if (m_bOpenRead || m_bOpenWrite) {
		m_IF.eprintf("CBQBFile::OpenRead(): Error: A file is already open.\n");
		return false;
	}

	m_pFile = fopen(s.c_str(),"rb");
	if (m_pFile == NULL) {
		m_IF.eprintf("CBQBFile::OpenRead(): Error: Could not open \"%s\" for reading.\n",s.c_str());
		return false;
	}
	m_sFileName = s;
	m_bOpenRead = true;
	m_bEOF = false;
	m_iFrameCounter = 0;

	// Extract directory part
	sp1 = m_sFileName.rfind('\\');
	sp2 = m_sFileName.rfind('/');
	if ((sp1 != std::string::npos) || (sp2 != std::string::npos)) {
		if (sp2 == std::string::npos)
			m_sDirectory = m_sFileName.substr(0,sp1+1);
		else if (sp1 == std::string::npos)
			m_sDirectory = m_sFileName.substr(0,sp2+1);
		else if (sp1 > sp2)
			m_sDirectory = m_sFileName.substr(0,sp1+1);
		else
			m_sDirectory = m_sFileName.substr(0,sp2+1);
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Directory is \"%s\".\n",m_sDirectory.c_str());
	} else
		m_sDirectory = "";

	memset(uc,0,5);
	(void)!fread(uc,5,1,m_pFile);
	if (memcmp(uc,"BQ",2) == 0) {
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    File header is \"BQ\" --> Single File.\n");
		m_bListFile = false;
	} else if (memcmp(uc,"BLIST",5) == 0) {
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    File header is \"BLIST\" --> List File.\n");
		m_bListFile = true;
		if (!OpenListFile(m_pFile)) {
			m_IF.eprintf("CBQBFile::OpenRead(): Error: Failed to open list file.\n");
			return false;
		}
		goto _end;
	} else {
		m_IF.eprintf("CBQBFile::OpenRead(): Error: File does not start with \"BQ\" or \"BLIST\".\n");
		return false;
	}

	if (m_pEngine == NULL) {
		m_pEngine = m_IF.CreateEngine(0);
		m_pEngine->PrepareDecompressCube();
	} else
		m_pEngine->Reset();

	ip = bqbtools.FindIndexPosition(m_pFile);

	if (ip != 0) {

		m_pIndex = new CBQBIndex(m_IF);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Found index frame at offset %ld, reading...\n",ip);

		fseek(m_pFile,ip,SEEK_SET);

		GetFileOffset(fohi,folo);

		if (!ReadFrame()) {
			m_IF.eprintf("CBQBFile::OpenRead(): Error: Could not read index frame.\n");
			return false;
		}

		if (m_pShortFrame->m_iFrameType == BQB_FRAMETYPE_COMPIDX) {

			if (m_pShortFrame->m_iFrameTypeVersion > 1) {
				m_IF.eprintf("CBQBFile::OpenRead(): Error: Unexpected frame version %d of compressed index frame.\n",
					m_pShortFrame->m_iFrameTypeVersion);
				return false;
			}

			if (m_pShortFrame->m_iFrameTypeVersion >= 1) {

				bs.m_iaData.assign( m_pShortFrame->m_iaPayload.begin(), m_pShortFrame->m_iaPayload.begin()+8 );
				bs.m_iExcessBits = 0;
				bs.m_iReadPosBytes = 0;
				bs.m_iReadPosExcess = 0;

				rohi = bs.ReadBitsInteger(32);
				rolo = bs.ReadBitsInteger(32);

				if ((rohi != fohi) || (rolo != folo))  {
					m_IF.eprintf("CBQBFile::OpenRead(): Warning: Index frame found at wrong offset; ignoring (found at %02lX:%08lX, expected at %02lX:%08lX).\n",
						fohi, folo, rohi, rolo );
					delete m_pIndex;
					m_pIndex = NULL;
					goto _noindex;
				}

				tuca.assign( m_pShortFrame->m_iaPayload.begin()+8, m_pShortFrame->m_iaPayload.end() );

			} else
				tuca.assign( m_pShortFrame->m_iaPayload.begin(), m_pShortFrame->m_iaPayload.end() );

			m_pIndex->ImportFromArray(
				true,
				tuca,
				m_pShortFrame->m_iFrameTypeVersion
			);

		} else if (m_pShortFrame->m_iFrameType == BQB_FRAMETYPE_IDX) {

			if (m_pShortFrame->m_iFrameTypeVersion > 1) {
				m_IF.eprintf("CBQBFile::OpenRead(): Error: Unexpected frame version %d of uncompressed index frame.\n",
					m_pShortFrame->m_iFrameTypeVersion);
				return false;
			}

			m_pIndex->ImportFromArray(
				false,
				m_pShortFrame->m_iaPayload,
				m_pShortFrame->m_iFrameTypeVersion
			);

		} else {

			m_IF.eprintf("CBQBFile::OpenRead(): Error: Unexpected frame type %d of index frame.\n",
				m_pShortFrame->m_iFrameType);
			return false;

		}

		m_iTotalFrameCount = (int)m_pIndex->m_iaFrameLengths.size();

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Positioning file pointer to start of file.\n");

	} else if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("    Found no index frame.\n");

_noindex:

	fseek(m_pFile,0,SEEK_SET);

	if (m_pShortFrame != NULL) {
		delete m_pShortFrame;
		m_pShortFrame = NULL;
	}

_end:

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::OpenRead <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


bool CBQBFile::ReadFrame() {

	unsigned char uc;
	unsigned char uci[2];
	unsigned int id, len, crc;
	int ty, ve;
	CBQBCRC32 ecrc;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::ReadFrame >>>\n");

	if (!m_bOpenRead) {
		m_IF.eprintf("\n");
		m_IF.eprintf("CBQBFile::ReadFrame(): Error: File not open for reading.\n");
		abort();
	}

	if (m_bListFile) {

		if (m_oaBQBList.size() == 0) {
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("    LIST: Error: No entries in list.\n");
			m_bEOF = true;
			return false;
		}

		if (m_iListIndex == -1) {
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("    LIST: Opening first entry from list...\n");
			m_iListIndex = 0;
			if (!m_oaBQBList[m_iListIndex]->m_pFile->OpenRead(m_oaBQBList[m_iListIndex]->m_sFileName)) {
				if (m_IF.IsPL(BQB_PL_VERBOSE))
					m_IF.printf("    LIST: Failed.\n");
				return false;
			}
			if (m_oaBQBList[m_iListIndex]->m_iFrameStart > 0) {
				if (m_IF.IsPL(BQB_PL_VERBOSE))
					m_IF.printf("    LIST: Skipping first %d entries...\n",m_oaBQBList[m_iListIndex]->m_iFrameStart);
				if (!m_oaBQBList[m_iListIndex]->m_pFile->SeekFrame(m_oaBQBList[m_iListIndex]->m_iFrameStart))
					return false;
			}
		}

_again:
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    LIST: Reading frame from list entry %d...\n",m_iListIndex+1);
		if (!m_oaBQBList[m_iListIndex]->m_pFile->ReadFrame()) {
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("    LIST: Could not read further frame from list entry %d. Closing entry.\n",m_iListIndex+1);
			if (!m_oaBQBList[m_iListIndex]->m_pFile->Close()) {
				if (m_IF.IsPL(BQB_PL_VERBOSE))
					m_IF.printf("    LIST: Failed to close entry.\n");
				return false;
			}
			if (m_iListIndex+1 == (int)m_oaBQBList.size()) {
				if (m_IF.IsPL(BQB_PL_VERBOSE))
					m_IF.printf("    LIST: No more list entries; leaving.\n");
				m_bEOF = true;
				m_iListIndex = -1;
				return false;
			}
			m_iListIndex++;
			if (!m_oaBQBList[m_iListIndex]->m_pFile->OpenRead(m_oaBQBList[m_iListIndex]->m_sFileName)) {
				if (m_IF.IsPL(BQB_PL_VERBOSE))
					m_IF.printf("    LIST: Error: Could not open list entry %d; leaving.\n",m_iListIndex+1);
				return false;
			}
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("    LIST: Opened list entry %d. Trying to read frame...\n",m_iListIndex+1);
			goto _again;
		}

		if ((m_oaBQBList[m_iListIndex]->m_pFile->m_pShortFrame->m_iFrameType == 2) ||
			(m_oaBQBList[m_iListIndex]->m_pFile->m_pShortFrame->m_iFrameType == 3)) {
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("    LIST: Is an index frame, skipping.\n");
			goto _again;
		}

		m_pShortFrame = m_oaBQBList[m_iListIndex]->m_pFile->m_pShortFrame;

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    LIST: Frame successfully read.\n");

	} else {

		if (m_pFile == NULL) {
			m_IF.eprintf("\n");
			m_IF.eprintf("CBQBFile::ReadFrame(): Error: File pointer is NULL.\n");
			abort();
		}

		(void)!fread(uci,2,1,m_pFile);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Offset is now %ld.\n",ftell(m_pFile));

		if (feof(m_pFile)) {
			m_bEOF = true;
			if (m_IF.IsPL(BQB_PL_VERBOSE)) {
				m_IF.printf("    Found EOF while reading frame header.\n");
				m_IF.printf("<<< CBQBFile::ReadFrame <<<\n");
				m_IF.printf("\n");
			}
			return false;
		}

		if (memcmp(uci,"BQ",2) != 0) {
			m_IF.eprintf("\n");
			m_IF.eprintf("CBQBFile::ReadFrame(): Error: Frame does not start with \"BQ\" (\"%s\", Frame %d).\n",
				m_sFileName.c_str(), m_iFrameCounter+1 );
			abort();
		}

		(void)!fread(&uc,1,1,m_pFile);

		ty = uc >> 3;
		ve = uc & 7;

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Found frame type %d, version %d.\n",ty,ve);

		if ((ty < 2) || (ty > 13)) {
			m_IF.eprintf("\n");
			m_IF.eprintf("CBQBFile::ReadFrame(): Error: Frame type %d is not a valid short frame (\"%s\", Frame %d).\n",
				ty, m_sFileName.c_str(), m_iFrameCounter+1 );
			abort();
		}

		if (m_pShortFrame != NULL)
			delete m_pShortFrame;
		m_pShortFrame = new CBQBShortFrame();

		m_pShortFrame->m_iFrameType = ty;
		m_pShortFrame->m_iFrameTypeVersion = ve;

		(void)!fread(&id,4,1,m_pFile);
		(void)!fread(&len,4,1,m_pFile);
		(void)!fread(&crc,4,1,m_pFile);

		m_pShortFrame->m_iID = id;

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Found id %u, payload length %u, crc32 %u.\n",id,len,crc);

		if (!ReadArrayFromFile(m_pFile,m_pShortFrame->m_iaPayload,len)) {
			m_IF.eprintf("\n");
			m_IF.eprintf("CBQBFile::ReadFrame(): Error: Unexpected end of file while reading frame payload (%u bytes; \"%s\", Frame %d).\n",
				len, m_sFileName.c_str(), m_iFrameCounter+1 );
			abort();
		}

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    %u bytes of payload read.\n",len);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Finished reading payload, CRC check...\n");

		m_pShortFrame->m_iCRC32 = ecrc.ComputeCRC32(m_pShortFrame->m_iaPayload);

		if (crc != m_pShortFrame->m_iCRC32) {
			m_IF.eprintf("\n");
			m_IF.eprintf("CBQBFile::ReadFrame(): Error: CRC check failed for frame payload (computed %08lX, read %08X; \"%s\", Frame %d).\n",
				m_pShortFrame->m_iCRC32, crc, m_sFileName.c_str(), m_iFrameCounter+1 );
			abort();
		}

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    CRC matches. All done.\n");
	}

	m_iFrameCounter++;

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::ReadFrame <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


bool CBQBFile::OpenWriteAppend(std::string s) {

	long ip;
	CBQBTools bqbtools(m_IF);
	CBQBBitSet bs(m_IF);
	unsigned long fohi, folo, rohi, rolo;
	std::vector<unsigned char> tuca;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::OpenWriteAppend >>>\n");

	if (m_bOpenRead || m_bOpenWrite) {
		m_IF.eprintf("CBQBFile::OpenWriteAppend(): Error: A file is already open.\n");
		abort();
	}

	m_pIndex = new CBQBIndex(m_IF);

	m_pFile = fopen(s.c_str(),"r+b");

	if (m_pFile != NULL) { // File exists

		if (m_IF.IsPL(BQB_PL_STANDARD)) {
	//		m_IF.bprintf("\n");
			m_IF.bprintf("Output BQB file \"%s\" already exists, appending.\n",s.c_str());
	//		m_IF.bprintf("\n");
		}

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    File already exists, looking for index frame...\n");

		ip = bqbtools.FindIndexPosition(m_pFile);

		if (ip != 0) {

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("    Found index frame at offset %ld, reading...\n",ip);

			m_bOpenRead = true;

			fseek(m_pFile,ip,SEEK_SET);

			GetFileOffset(fohi,folo);

			ReadFrame();

			m_bOpenRead = false;

			if (m_pShortFrame->m_iFrameType == BQB_FRAMETYPE_COMPIDX) {

				if (m_pShortFrame->m_iFrameTypeVersion > 1) {
					m_IF.eprintf("CBQBFile::OpenWriteAppend(): Error: Unexpected frame version %d of compressed index frame.\n",
						m_pShortFrame->m_iFrameTypeVersion);
					abort();
				}

				if (m_pShortFrame->m_iFrameTypeVersion >= 1) {

					bs.m_iaData.assign( m_pShortFrame->m_iaPayload.begin(), m_pShortFrame->m_iaPayload.begin()+8 );
					bs.m_iExcessBits = 0;
					bs.m_iReadPosBytes = 0;
					bs.m_iReadPosExcess = 0;

					rohi = bs.ReadBitsInteger(32);
					rolo = bs.ReadBitsInteger(32);

					if ((rohi != fohi) || (rolo != folo))  {
						m_IF.eprintf("CBQBFile::OpenWriteAppend(): Warning: Index frame found at wrong offset; ignoring (found at %02lX:%08lX, expected at %02lX:%08lX).\n",
							fohi, folo, rohi, rolo );
						goto _noindex;
					}

					tuca.assign( m_pShortFrame->m_iaPayload.begin()+8, m_pShortFrame->m_iaPayload.end() );

				} else
					tuca.assign( m_pShortFrame->m_iaPayload.begin(), m_pShortFrame->m_iaPayload.end() );

				m_pIndex->ImportFromArray(
					true,
					tuca,
					m_pShortFrame->m_iFrameTypeVersion
				);

			} else if (m_pShortFrame->m_iFrameType == BQB_FRAMETYPE_IDX) {

				if (m_pShortFrame->m_iFrameTypeVersion > 1) {
					m_IF.eprintf("CBQBFile::OpenWriteAppend(): Error: Unexpected frame version %d of uncompressed index frame.\n",
						m_pShortFrame->m_iFrameTypeVersion);
					abort();
				}

				m_pIndex->ImportFromArray(
					false,
					m_pShortFrame->m_iaPayload,
					m_pShortFrame->m_iFrameTypeVersion
				);

			} else {

				m_IF.eprintf("CBQBFile::OpenWriteAppend(): Error: Unexpected frame type %d of index frame.\n",
					m_pShortFrame->m_iFrameType);
				abort();

			}

			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("    Positioning file pointer to start of old index frame.\n");

			fseek(m_pFile,ip,SEEK_SET);

		} else if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Found no index frame.\n");

_noindex:;

	} else { // File does not exist

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    File does not exist, creating...\n");

		m_pFile = fopen(s.c_str(),"wb");
		if (m_pFile == NULL) {
			m_IF.eprintf("CBQBFile::OpenWriteAppend(): Error: Could not create \"%s\" for writing.\n",
				s.c_str());
			abort();
		}

	}

	m_sFileName = s;
	m_bOpenWrite = true;

	if (m_pShortFrame != NULL) {
		delete m_pShortFrame;
		m_pShortFrame = NULL;
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::OpenWriteAppend <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


bool CBQBFile::OpenWriteReplace(std::string s) {

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::OpenWriteReplace >>>\n");

	if (m_bOpenRead || m_bOpenWrite) {
		m_IF.eprintf("CBQBFile::OpenWriteReplace(): Error: A file is already open.\n");
		abort();
	}

	m_pIndex = new CBQBIndex(m_IF);

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("    Creating file...\n");

	m_pFile = fopen(s.c_str(),"wb");
	if (m_pFile == NULL) {
		m_IF.eprintf("CBQBFile::OpenWriteReplace(): Error: Could not create \"%s\" for writing.\n",s.c_str());
		abort();
	}

	m_sFileName = s;
	m_bOpenWrite = true;

	if (m_pShortFrame != NULL) {
		delete m_pShortFrame;
		m_pShortFrame = NULL;
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::OpenWriteReplace <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


bool CBQBFile::Rewind() {

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::Rewind() >>>\n");

	if (!m_bOpenRead && !m_bOpenWrite) {
		m_IF.eprintf("CBQBFile::Rewind(): Error: No file is open.\n");
		abort();
	}

	if (m_bListFile) {

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    LIST: Rewinding list file...\n");
		if (m_iListIndex != -1)
			m_oaBQBList[m_iListIndex]->m_pFile->Close();
		m_iListIndex = -1;

	} else {

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Rewinding bqb file...\n");

		if (m_pFile == NULL) {
			m_IF.eprintf("CBQBFile::Rewind(): Error: File pointer is NULL.\n");
			abort();
		}
		
		rewind(m_pFile);

		if (m_pEngine == NULL) {
			//m_pEngine = new CBQBEngine();
			m_pEngine = m_IF.CreateEngine(0);
			m_pEngine->PrepareDecompressCube();
		} else
			m_pEngine->Reset();
	}

	m_bEOF = false;
	m_iFrameCounter = 0;

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::Rewind <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


void CBQBFile::GetFileOffset(unsigned long &high, unsigned long &low) {

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::GetFileOffset >>>\n");

	if (!m_bOpenRead && !m_bOpenWrite) {
		m_IF.eprintf("CBQBFile::GetFileOffset(): Error: No file is open.\n");
		abort();
	}

	// It is not so easy to get the 64 bit file offset
	// in a platform-independent portable way...
	#ifdef __CYGWIN__
		low = ftell(m_pFile) & 0xFFFFFFFF;
		high = ftell(m_pFile) >> 32;
	#elif defined (_WIN32)
		low = _telli64(_fileno(m_pFile)) & 0xFFFFFFFF;
		high = _telli64(_fileno(m_pFile)) >> 32;
	#else
		low = ftello(m_pFile) & 0xFFFFFFFF;
		high = ftello(m_pFile) >> 32;
	#endif

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("<<< CBQBFile::GetFileOffset <<<\n");
}


CBQBEngine* CBQBFile::GetEngine() {

	if (!m_bOpenRead && !m_bOpenWrite) {
		m_IF.eprintf("CBQBFile::GetEngine(): Error: No file is open.\n");
		abort();
	}

	if (m_bListFile) {

		if (m_oaBQBList.size() == 0) {
			m_IF.eprintf("CBQBFile::GetEngine(): Error: No entries in list.\n");
			abort();
		}

		if (m_iListIndex == -1) {
			m_IF.eprintf("CBQBFile::GetEngine(): Error: No entry from list file is open.\n");
			abort();
		}

		if (m_oaBQBList[m_iListIndex]->m_pFile == NULL) {
			m_IF.eprintf("CBQBFile::GetEngine(): Error: m_pFile == NULL in current list entry.\n");
			abort();
		}

		if (m_oaBQBList[m_iListIndex]->m_pFile->m_pEngine == NULL) {
			m_IF.eprintf("CBQBFile::GetEngine(): Error: m_pFile->m_pEngine == NULL in current list entry.\n");
			abort();
		}

		return m_oaBQBList[m_iListIndex]->m_pFile->m_pEngine;
	}

	if (m_pEngine == NULL) {
		m_IF.eprintf("CBQBFile::GetEngine(): Error: m_pEngine == NULL.\n");
		abort();
	}

	return m_pEngine;
}


bool CBQBFile::Close() {

	int z;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::Close >>>\n");

	if (!m_bOpenRead && !m_bOpenWrite) {
		m_IF.eprintf("CBQBFile::Close(): Error: No file is open.\n");
		abort();
	}

	if (m_bListFile) {

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    LIST: Closing list file...\n");
		if (m_iListIndex != -1)
			m_oaBQBList[m_iListIndex]->m_pFile->Close();
		m_iListIndex = -1;
		for (z=0;z<(int)m_oaBQBList.size();z++)
			delete m_oaBQBList[z];
		m_oaBQBList.clear();
		m_pShortFrame = NULL; // Managed/deleted by subsequent instances

	} else {

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Closing BQB file...\n");

		if (m_pFile == NULL) {
			m_IF.eprintf("CBQBFile::Close(): Error: File pointer is NULL.\n");
			abort();
		}
		
		fclose(m_pFile);
		m_pFile = NULL;
	}

	m_sFileName = "";

	if (m_pIndex != NULL) {
		delete m_pIndex;
		m_pIndex = NULL;
	}

	if (m_pShortFrame != NULL) {
		delete m_pShortFrame;
		m_pShortFrame = NULL;
	}

	m_bOpenRead = false;
	m_bOpenWrite = false;
	m_iTotalFrameCount = -1;
	m_bListFile = false;

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::Close <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


bool CBQBFile::CreateGeneralFrame() {

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::CreateGeneralFrame >>>\n");

	m_IF.eprintf("CBQBFile::CreateGeneralFrame(): Not yet implemented.\n");
	abort();

	return true;
}


bool CBQBFile::CreateShortFrame(int type, int version, int id) {

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::CreateShortFrame >>>\n");

	if (!m_bOpenWrite) {
		m_IF.eprintf("CBQBFile::CreateShortFrame(): Error: No file is open for writing.\n");
		abort();
	}

	if ((type < 2) || (type > 13)) {
		m_IF.eprintf("CBQBFile::CreateShortFrame(): Error: Short frame type must be in range 2 .. 9 (specified: %d).\n",type);
		abort();
	}

	if ((version < 0) || (version > 7)) {
		m_IF.eprintf("CBQBFile::CreateShortFrame(): Error: Frame type version must be in range 0 .. 7 (specified: %d).\n",version);
		abort();
	}

	m_bShortFrame = true;

	if (m_pShortFrame != NULL)
		delete m_pShortFrame;
	m_pShortFrame = new CBQBShortFrame();
	m_pShortFrame->m_iFrameType = type;
	m_pShortFrame->m_iFrameTypeVersion = version;
	m_pShortFrame->m_iID = id;

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::CreateShortFrame <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


bool CBQBFile::PushPayload(const std::vector<unsigned char> &ia) {

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::PushPayload >>>\n");

	if (!m_bOpenWrite) {
		m_IF.eprintf("CBQBFile::PushPayload(): Error: No file is open for writing.\n");
		abort();
	}

	if (!m_bShortFrame) {
		m_IF.eprintf("CBQBFile::PushPayload(): Error: Only applicable to short frames.\n");
		abort();
	}

	if (m_pShortFrame == NULL) {
		m_IF.eprintf("CBQBFile::PushPayload(): Error: No frame has been created before.\n");
		abort();
	}

	m_pShortFrame->m_iaPayload.insert( m_pShortFrame->m_iaPayload.end(), ia.begin(), ia.end() );

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::PushPayload <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


bool CBQBFile::FinalizeFrame(CBQBStatistics *stat) {

	unsigned char uc, tv;
	unsigned int ui;
	CBQBCRC32 crc;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::FinalizeFrame >>>\n");

	if (!m_bOpenWrite) {
		m_IF.eprintf("CBQBFile::FinalizeFrame(): Error: No file is open for writing.\n");
		abort();
	}

	if (m_bShortFrame) {

		if (m_pShortFrame == NULL) {
			m_IF.eprintf("CBQBFile::FinalizeFrame(): Error: No frame has been created before.\n");
			abort();
		}

		if (m_pShortFrame->m_iaPayload.size() == 0) {
			m_IF.eprintf("CBQBFile::FinalizeFrame(): Error: Payload is empty.\n");
			abort();
		}

		#ifndef _MSC_VER
			if (m_bIndexWritten) {
	//			m_IF.printf("@ Restoring old offset (%lu) to overwrite index frame.\n",m_iIndexOffset);
				m_bIndexWritten = false;
				fseeko( m_pFile, m_iIndexOffset, SEEK_SET );
				m_pIndex->m_iaFrameLengths.pop_back();
				m_pIndex->m_iaFrameTypes.pop_back();
				m_pIndex->m_iaFrameIDs.pop_back();
			}
		#endif

		uc = 'B';
		(void)!fwrite(&uc,1,1,m_pFile);
		uc = 'Q';
		(void)!fwrite(&uc,1,1,m_pFile);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Offset %ld, writing frame type %d.\n",ftell(m_pFile),m_pShortFrame->m_iFrameType);

		tv = (unsigned char)((m_pShortFrame->m_iFrameType << 3) + m_pShortFrame->m_iFrameTypeVersion);
		(void)!fwrite(&tv,1,1,m_pFile);

		ui = m_pShortFrame->m_iID;
		(void)!fwrite(&ui,4,1,m_pFile);

		ui = (unsigned int)m_pShortFrame->m_iaPayload.size();
		(void)!fwrite(&ui,4,1,m_pFile);

		m_pShortFrame->m_iCRC32 = crc.ComputeCRC32(m_pShortFrame->m_iaPayload);

		ui = m_pShortFrame->m_iCRC32;
		(void)!fwrite(&ui,4,1,m_pFile);

		if (stat != NULL)
			stat->m_oStat.m_lOverhead += 15 * 8;

		WriteArrayToFile(m_pFile,m_pShortFrame->m_iaPayload);

		fflush(m_pFile);

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    %lu bytes of payload written, CRC32 is %08lX.\n",(unsigned long)m_pShortFrame->m_iaPayload.size(),m_pShortFrame->m_iCRC32);

		m_pIndex->m_iaFrameLengths.push_back((int)m_pShortFrame->m_iaPayload.size()+15);
		m_pIndex->m_iaFrameTypes.push_back(tv);
		m_pIndex->m_iaFrameIDs.push_back(m_pShortFrame->m_iID);

		delete m_pShortFrame;
		m_pShortFrame = NULL;

	} else {

		m_IF.eprintf("CBQBFile::FinalizeFrame(): Error: Not yet implemented for general frames.\n");
		abort();

	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::FinalizeFrame <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


bool CBQBFile::WriteIndexFrame(bool compressed, CBQBStatistics *stat) {

	std::vector<unsigned char> uca;
	unsigned long ul, fohi, folo;
	CBQBBitSet bs(m_IF);


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::WriteIndexFrame >>>\n");

	if (!m_bOpenWrite) {
		m_IF.eprintf("CBQBFile::WriteIndexFrame(): Error: No file is open for writing.\n");
		abort();
	}

	if (m_pIndex->m_iaFrameTypes.size() == 0) {
		m_IF.eprintf("CBQBFile::WriteIndexFrame(): Warning: Will not write index for empty file.\n");
		goto _ende;
	}

	#ifndef _MSC_VER
		m_iIndexOffset = ftello(m_pFile);
//		m_IF.printf("@ Saving offset (%lu) before writing index frame.\n",m_iIndexOffset);
	#endif

	if (compressed) {

		GetFileOffset(fohi,folo);

		// Compressed Index Frame Version 1
		CreateShortFrame(BQB_FRAMETYPE_COMPIDX,1,0);

		bs.WriteBits(fohi,32);
		bs.WriteBits(folo,32);

		PushPayload(bs.m_iaData);

		m_pIndex->ExportToArray(compressed,uca,stat);

		PushPayload(uca);

		if (stat != NULL)
			stat->m_oStat.m_lOverhead += (int)((8+uca.size())*8);

		ul = (unsigned long)(8 + uca.size() + 22);

		uca.clear();

		uca.push_back((ul>>24)&0xFF);
		uca.push_back((ul>>16)&0xFF);
		uca.push_back((ul>>8)&0xFF);
		uca.push_back(ul&0xFF);
		uca.push_back('I');
		uca.push_back('D');
		uca.push_back('X');

		PushPayload(uca);

		if (stat != NULL)
			stat->m_oStat.m_lOverhead += (int)(uca.size()*8);

		FinalizeFrame(stat);

	} else {

		m_IF.eprintf("CBQBFile::WriteIndexFrame(): Non-compressed index frames not yet implemented.\n");
		abort();
	}

	#ifndef _MSC_VER
		m_bIndexWritten = true;
	#endif

_ende:
	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::WriteIndexFrame <<<\n");
		m_IF.printf("\n");
	}

	return true;
}


void CBQBFile::DumpIndex() {

	m_pIndex->Dump();
}


bool CBQBFile::CheckIntegrity(std::string s, bool verbose) {

	UNUSED(verbose);

	int i;
	unsigned long crc;
	CBQBCRC32 ecrc;
	bool b;


	m_IF.printf("####################################\n");
	m_IF.printf("###   BQB File Integrity Check   ###\n");
	m_IF.printf("####################################\n");
	m_IF.printf("\n");

	m_IF.printf("Trying to open %s ...\n",s.c_str());

	if (!OpenRead(s)) {
		m_IF.eprintf("CBQBFile::CheckIntegrity(): Error: Failed to open %s\n",s.c_str());
		return false;
	}

	if (m_pIndex->m_iaFrameLengths.size() != 0) {

		m_IF.printf("\n");
		m_IF.printf("File contains an index with %lu/%lu entries.\n",(unsigned long)m_pIndex->m_iaFrameTypes.size(),(unsigned long)m_pIndex->m_iaFrameLengths.size());
		m_IF.printf("\n");

		i = 0;
		b = false;
		while (ReadFrame()) {
			m_IF.printf("    Frame %6d,  type %2d v%d (%s),  ID %6d,  payload size %8lu bytes",
				i+1,
				m_pShortFrame->m_iFrameType,
				m_pShortFrame->m_iFrameTypeVersion,
				GetFrameTypeString(m_pShortFrame->m_iFrameType),
				m_pShortFrame->m_iID,
				(unsigned long)m_pShortFrame->m_iaPayload.size()
			);
			if (b) {
				m_IF.eprintf("\n");
				m_IF.eprintf("CBQBFile::CheckIntegrity(): Error: Frame %d not on the index, but was not the last frame.\n",i);
				Close();
				return false;
			}
			crc = ecrc.ComputeCRC32(m_pShortFrame->m_iaPayload);
			if (crc != m_pShortFrame->m_iCRC32) {
				m_IF.eprintf(",  CRC fail.\n");
				m_IF.eprintf("CBQBFile::CheckIntegrity(): Error in frame %d: CRC mismatch (computed %08lX, read %08lX).\n",i+1,crc,m_pShortFrame->m_iCRC32);
				Close();
				return false;
			}
			m_IF.printf(",  CRC matches.\n");
			if ((int)m_pIndex->m_iaFrameTypes.size() > i) {
				if ((m_pShortFrame->m_iFrameType != (int)(m_pIndex->m_iaFrameTypes[i]>>3)) || (m_pShortFrame->m_iFrameTypeVersion != (int)(m_pIndex->m_iaFrameTypes[i]&7)) || ((int)m_pShortFrame->m_iaPayload.size()+15 != m_pIndex->m_iaFrameLengths[i])) {
					m_IF.eprintf("CBQBFile::CheckIntegrity(): Error in frame %d: Mismatch between frame header and index: Type %d<-->%d, Version %d<-->%d, Length %lu<-->%d.\n",i+1,m_pShortFrame->m_iFrameType,m_pIndex->m_iaFrameTypes[i]>>3,m_pShortFrame->m_iFrameTypeVersion,m_pIndex->m_iaFrameTypes[i]&7,(unsigned long)m_pShortFrame->m_iaPayload.size()+15,m_pIndex->m_iaFrameLengths[i]);
					Close();
					return false;
				}
			} else {
				b = true;
				m_IF.printf("      This frame is not on the index...\n");
			}
			i++;
		}

	} else {

		m_IF.printf("File does not contain an index.\n");

		i = 0;
		while (ReadFrame()) {
			m_IF.printf("    Frame %6d,  type %2d v%d (%s),  ID %6d,  payload size %8lu bytes",
				i+1,
				m_pShortFrame->m_iFrameType,
				m_pShortFrame->m_iFrameTypeVersion,
				GetFrameTypeString(m_pShortFrame->m_iFrameType),
				m_pShortFrame->m_iID,
				(unsigned long)m_pShortFrame->m_iaPayload.size()
			);
			crc = ecrc.ComputeCRC32(m_pShortFrame->m_iaPayload);
			if (crc != m_pShortFrame->m_iCRC32) {
				m_IF.eprintf(",  CRC fail.\n");
				m_IF.eprintf("CBQBFile::CheckIntegrity(): Error in frame %d: CRC mismatch (computed %08lX, read %08lX).\n",i+1,crc,m_pShortFrame->m_iCRC32);
				Close();
				return false;
			}
			m_IF.printf(",  CRC matches.\n");
			i++;
		}

	}

	Close();

	m_IF.printf("\n");
	m_IF.printf("Integrity check passed.\n");
	m_IF.printf("\n");

	return true;
}


int CBQBFile::GetFrameID() const {

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::GetFrameID >>>\n");

	if (!m_bOpenRead) {
		m_IF.eprintf("CBQBFile::GetFrameID(): Error: File not open for reading.\n");
		abort();
	}

	if (m_pShortFrame == NULL) {
		m_IF.eprintf("CBQBFile::GetFrameID(): Error: No current frame.\n");
		abort();
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::GetFrameID <<<\n");
		m_IF.printf("\n");
	}

	return m_pShortFrame->m_iID;
}


int CBQBFile::GetFrameType() const {

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> GetFrameType::GetFrameIndex >>>\n");

	if (!m_bOpenRead) {
		m_IF.eprintf("CBQBFile::GetFrameType(): Error: File not open for reading.\n");
		abort();
	}

	if (m_pShortFrame == NULL) {
		m_IF.eprintf("CBQBFile::GetFrameType(): Error: No current frame.\n");
		abort();
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::GetFrameType <<<\n");
		m_IF.printf("\n");
	}

	return m_pShortFrame->m_iFrameType;
}


int CBQBFile::GetFrameTypeVersion() const {

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::GetFrameTypeVersion >>>\n");

	if (!m_bOpenRead) {
		m_IF.eprintf("CBQBFile::GetFrameTypeVersion(): Error: File not open for reading.\n");
		abort();
	}

	if (m_pShortFrame == NULL) {
		m_IF.eprintf("CBQBFile::GetFrameTypeVersion(): Error: No current frame.\n");
		abort();
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::GetFrameTypeVersion <<<\n");
		m_IF.printf("\n");
	}

	return m_pShortFrame->m_iFrameTypeVersion;
}


const std::vector<unsigned char>* CBQBFile::GetFramePayload() const {

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::GetFramePayload >>>\n");

	if (!m_bOpenRead) {
		m_IF.eprintf("CBQBFile::GetFramePayload(): Error: File not open for reading.\n");
		abort();
	}

	if (m_pShortFrame == NULL) {
		m_IF.eprintf("CBQBFile::GetFramePayload(): Error: No current frame.\n");
		abort();
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		m_IF.printf("<<< CBQBFile::GetFramePayload <<<\n");
		m_IF.printf("\n");
	}

	return &m_pShortFrame->m_iaPayload;
}


const char *GetFrameTypeString(int type) {

	const char* names[] = {
		"RESERVED     ",
		"GENERAL      ",
		"IDX          ",
		"COMPIDX      ",
		"TRAJ         ",
		"COMPTRAJSTART",
		"COMPTRAJ     ",
		"CUBE         ",
		"COMPCUBESTART",
		"COMPCUBE     ",
		"FILE         ",
		"COMPFILE     ",
		"COMPTRAJKEY  ",
		"COMPCUBEKEY  "
	};

	if ((type >= 0) && (type < 14))
		return names[type];
	else
		return "UNKNOWN      ";
}


bool CBQBFile::OpenListFile(FILE *a) {

	char buf[1024], buf2[32], bstart[16], bend[16], *p, *q, *r;
	CBQBListEntry *le;
	bool b;
	int z, ti;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::OpenListFile >>>\n");

	fseek(a,0,SEEK_SET);

	(void)!fgets(buf,1024,a);
	while ((buf[strlen(buf)-1] == '\r') || (buf[strlen(buf)-1] == '\n'))
		buf[strlen(buf)-1] = 0;
	if (strcmp(buf,"BLIST") != 0) {
		m_IF.eprintf("CBQBFile::OpenListFile(): Error: File header mismatch (\"%s\" != \"BLIST\").\n",
			buf);
		return false;
	}

	b = false;
	m_iTotalFrameCount = 0;
	m_iListIndex = -1;
	while (!feof(a)) {
		if (fgets(buf,1024,a) == NULL) {
			if (!feof(a)) {
				m_IF.eprintf("CBQBFile::OpenListFile(): Error while reading line %lu from BQB list file.\n",
					(unsigned long)m_oaBQBList.size()+1);
				return false;
			} else
				break;
		}
		while ((buf[strlen(buf)-1] == '\r') || (buf[strlen(buf)-1] == '\n'))
			buf[strlen(buf)-1] = 0;
		if (strlen(buf) == 0)
			continue;
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    LIST: Checking list entry %3d: \"%s\"...\n",
				(int)m_oaBQBList.size()+1,buf);
		le = new CBQBListEntry();

		p = buf;
		while (*p == ' ')
			p++;
		if (*p == '[') {
			q = p+1;
			while ((*q != ']') && (*q != 0))
				q++;
			if (*q == 0) {
				m_IF.eprintf("CBQBFile::OpenListFile(): Missing \"]\" in line %lu: \"%s\".\n",
					(unsigned long)m_oaBQBList.size()+1,buf);
				return false;
			}
			r = p+1;
			while ((r < q) && (*r != '-'))
				r++;
			if (r == q) {
				m_IF.eprintf("CBQBFile::OpenListFile(): Missing \"-\" in line %lu: \"%s\".\n",
					(unsigned long)m_oaBQBList.size()+1,buf);
				return false;
			}
			if (p+1 < r)
				memcpy(buf2,p+1,r-p-1);
			buf2[r-p-1] = 0;
			le->m_iFrameStart = atoi(buf2)-1;
			if (r+1 < q)
				memcpy(buf2,r+1,q-r-1);
			buf2[q-r-1] = 0;
			le->m_iFrameEnd = atoi(buf2)-1;
			le->m_sFileName = q+1;
		} else
			le->m_sFileName = buf;

		if (le->m_sFileName[0] == '/')
			goto _abs; // Absolute path (Unix)

		if (le->m_sFileName.length() > 2)
			if (le->m_sFileName[1] == ':')
				goto _abs; // Absolute path (Windows)

		// Prepend directory of list file to relative path of entry
		if (m_sDirectory.length() > 0) {
			le->m_sFileName = m_sDirectory + le->m_sFileName;
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("      LIST: Reconstructed path: \"%s\"\n", le->m_sFileName.c_str() );
		}
_abs:

		le->m_pFile = new CBQBFile(m_IF);
		if (!le->m_pFile->OpenRead(le->m_sFileName)) {
			m_IF.eprintf("CBQBFile::OpenListFile(): Error while opening entry %d: \"%s\".\n",
				(int)m_oaBQBList.size()+1,buf);
			return false;
		}

		if (le->m_pFile->m_pIndex != NULL) {
			le->m_pIndex = new CBQBIndex(*le->m_pFile->m_pIndex);
			ti = 0;
			for (z=0;z<(int)le->m_pIndex->m_iaFrameLengths.size();z++)
				if ((le->m_pIndex->m_iaFrameTypes[z] != 2) && (le->m_pIndex->m_iaFrameTypes[z] != 3))
					ti++;
			le->m_iFullFrameCount = ti/*le->m_pIndex->m_iaFrameLengths.size()*/;
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("      LIST: Found index frame: %d frames, %d payload frames.\n",
					(int)le->m_pIndex->m_iaFrameLengths.size(),le->m_iFullFrameCount);
		} else {
			b = true;
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("      LIST: Found no index.\n");
		}

		if ((le->m_iFrameStart != -1) && (le->m_iFrameEnd != -1)) {
			le->m_iFrameCount = le->m_iFrameEnd - le->m_iFrameStart + 1;
			if ((le->m_iFullFrameCount != -1) && (le->m_iFrameCount > le->m_iFullFrameCount)) {
				m_IF.eprintf("CBQBFile::OpenListFile(): Error in entry %d: Frame range %d - %d is more than total frame count (%d).\n",
					(int)m_oaBQBList.size()+1,le->m_iFrameStart+1,le->m_iFrameEnd+1,le->m_iFullFrameCount);
				return false;
			}
		} else if ((le->m_iFrameStart != -1) && (le->m_iFullFrameCount != -1)) {
			le->m_iFrameCount = le->m_iFullFrameCount - le->m_iFrameStart;
			if (le->m_iFrameStart+1 > le->m_iFullFrameCount) {
				m_IF.eprintf("CBQBFile::OpenListFile(): Error in entry %d: Starting frame number %d is larger than total frame count (%d).\n",
					(int)m_oaBQBList.size()+1,le->m_iFrameStart+1,le->m_iFullFrameCount);
				return false;
			}
		} else if ((le->m_iFrameEnd != -1) && (le->m_iFullFrameCount != -1)) {
			le->m_iFrameCount = le->m_iFrameEnd + 1;
			if (le->m_iFrameEnd+1 > le->m_iFullFrameCount) {
				m_IF.eprintf("CBQBFile::OpenListFile(): Error in entry %d: Ending frame number %d is larger than total frame count (%d).\n",
					(int)m_oaBQBList.size()+1,le->m_iFrameEnd+1,le->m_iFullFrameCount);
				return false;
			}
		} else
			le->m_iFrameCount = le->m_iFullFrameCount;

		if (le->m_iFrameCount != -1)
			m_iTotalFrameCount += le->m_iFrameCount;

		le->m_pFile->Close();
		m_oaBQBList.push_back(le);
	}

	if (!b) {
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("LIST: Building overall index...\n");
		m_pIndex = new CBQBIndex(m_IF);
		for (z=0;z<(int)m_oaBQBList.size();z++) {
			if (m_oaBQBList[z]->m_pIndex == NULL) {
				m_IF.eprintf("CBQBFile::OpenListFile(): Internal Error: Index of list element %d is NULL.\n",z+1);
				return false;
			}
			m_pIndex->m_iaFrameLengths.insert(
				m_pIndex->m_iaFrameLengths.end(),
				m_oaBQBList[z]->m_pIndex->m_iaFrameLengths.begin(),
				m_oaBQBList[z]->m_pIndex->m_iaFrameLengths.end()
			);
			m_pIndex->m_iaFrameTypes.insert(
				m_pIndex->m_iaFrameTypes.end(),
				m_oaBQBList[z]->m_pIndex->m_iaFrameTypes.begin(),
				m_oaBQBList[z]->m_pIndex->m_iaFrameTypes.end()
			);
		}
	} else {
		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("LIST: Some list entries have no index; not building overall index.\n");
		m_iTotalFrameCount = -1;
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("LIST: Added %lu entries to list.\n",(unsigned long)m_oaBQBList.size());

	if (m_IF.IsPL(BQB_PL_VERBOSE)) {
		if (m_iTotalFrameCount != -1)
			m_IF.printf("LIST: Found %d frames in total.\n",m_iTotalFrameCount);
		else
			m_IF.printf("LIST: Total frame count is unknown.\n");
	}

	fclose(a);

	m_IF.printf("\n");
	m_IF.printf("    Contents of BQB list file %s:\n",m_sFileName.c_str());
	for (z=0;z<(int)m_oaBQBList.size();z++) {

		if (m_oaBQBList[z]->m_iFrameStart != -1)
			sprintf(bstart,"%5d",m_oaBQBList[z]->m_iFrameStart+1);
		else
			sprintf(bstart,"start");

		if (m_oaBQBList[z]->m_iFrameEnd != -1)
			sprintf(bend,"%5d",m_oaBQBList[z]->m_iFrameEnd+1);
		else
			sprintf(bend,"  end");

		if (m_oaBQBList[z]->m_pIndex != NULL)
			m_IF.printf("    %3d.) %5d payload frames, using %5d (%5s - %5s), %s\n",
				z+1,m_oaBQBList[z]->m_iFullFrameCount,m_oaBQBList[z]->m_iFrameCount,bstart,
				bend,m_oaBQBList[z]->m_sFileName.c_str());
		else if (m_oaBQBList[z]->m_iFrameCount != -1)
			m_IF.printf("    %3d.) (no index),           using %5d (%5s - %5s), %s\n",
				z+1,m_oaBQBList[z]->m_iFrameCount,bstart,bend,m_oaBQBList[z]->m_sFileName.c_str());
		else
			m_IF.printf("    %3d.) (no index),                       (%5s - %5s), %s\n",
				z+1,bstart,bend,m_oaBQBList[z]->m_sFileName.c_str());
	}
	m_IF.printf("\n");

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("<<< CBQBFile::OpenListFile <<<\n");

	return true;
}


bool CBQBFile::IsEOF() {

	return m_bEOF;
}


bool CBQBFile::CompareCoords(std::string infile, std::string reffile, bool verbose) {

	UNUSED(infile);
	UNUSED(reffile);
	UNUSED(verbose);

/*	FILE *a;
	int i, z, z2, ty, ver;
	CTimeStep ts;
	CxDVec3Array da;
	double cv[3], m, t1, t2, tv[3];
	bool wrap;
	CCubeFrame *cfr;
	CAtomSet *afr;
	bool fail;

	m_IF.printf("\n");
	m_IF.printf(WHITE,"    *********************************\n");
	m_IF.printf(WHITE,"    ***   Coordinate Comparison   ***\n");
	m_IF.printf(WHITE,"    *********************************\n\n");

	m_IF.printf("    Will compare the coordinates in BQB file %s\n",infile.c_str());
	m_IF.printf("    to the reference coordinates in XYZ file %s.\n",reffile.c_str());
	m_IF.printf("\n");

	SetBarbecubeVerbose(verbose);

	m_IF.printf("    Opening BQB file...\n");
	if (!OpenRead(infile))
		return false;

	m_IF.printf("\n");
	m_IF.printf("    Opening XYZ file...\n");
	a = fopen(reffile.c_str(),"rb");
	if (a == NULL) {
		m_IF.eprintf("Error: Could not open file for reading.\n");
		return false;
	}

	g_pCCEngine = new CCCEngine();
	g_pCCEngine->PrepareDecompressCube(false);

	m_IF.printf("\n");
	if (AskYesNo("    Allow wrapping (y) or compare absolute coordinates (n)? [yes] ",true)) {

		wrap = true;

		m_IF.printf("\n");
		m_IF.printf("    Reading first frame...\n");
		if (!ReadFrame()) {
			fclose(a);
			return false;
		}

		ty = GetFrameType();
		ver = GetFrameTypeVersion();

		m_IF.printf("    First frame is of type %d v%d.\n",ty,ver);

		if ((ty != 8) || (ver > 1)) {
			m_IF.eprintf("Error: Frame type not yet supported.\n");
			fclose(a);
			return false;
		}

		if (!g_pCCEngine->DecompressCubeStep(GetFramePayload(),ver,false)) {
			fclose(a);
			return false;
		}

		cfr = g_pCCEngine->GetOutputCubeFrame(0);

		m_IF.printf("\n    *** Input of cell vector (in pm) ***\n\n");
		cv[0] = AskFloat("    Enter X component: [%11.4f] ",cfr->m_fStrideA[0]*cfr->m_iRes[0]*LEN_AU2PM,cfr->m_fStrideA[0]*cfr->m_iRes[0]*LEN_AU2PM);
		cv[1] = AskFloat("    Enter Y component: [%11.4f] ",cfr->m_fStrideB[1]*cfr->m_iRes[1]*LEN_AU2PM,cfr->m_fStrideB[1]*cfr->m_iRes[1]*LEN_AU2PM);
		cv[2] = AskFloat("    Enter Z component: [%11.4f] ",cfr->m_fStrideC[2]*cfr->m_iRes[2]*LEN_AU2PM,cfr->m_fStrideC[2]*cfr->m_iRes[2]*LEN_AU2PM);

		m_IF.printf("\n");
		m_IF.printf("    Rewinding BQB file...\n");

		if (!Rewind()) {
			fclose(a);
			return false;
		}

	} else
		wrap = false;

	i = 0;
	fail = false;
	m_IF.printf("\n");
	m_IF.printf("Starting comparison process.\n");
	while (true) {

		if (!ReadFrame())
			break;

		m_IF.printf("  Frame %4d ...  ",i+1);

		m_IF.printf("Type %2dv%1d,  ",GetFrameType(),GetFrameTypeVersion());

		if ((GetFrameType() == 2) || (GetFrameType() == 3)) {
			m_IF.printf("Skipping index frame.\n");
			continue;
		} else if (((GetFrameType() != 8) && (GetFrameType() != 9)) || (GetFrameTypeVersion() != 0)) {
			m_IF.eprintf("\nError: Frame type %d v%d not yet supported.\n",GetFrameType(),GetFrameTypeVersion());
			fclose(a);
			return false;
		}

		if (!g_pCCEngine->DecompressCubeStep(GetFramePayload(),verbose,false))
			return false;

		afr = g_pCCEngine->GetOutputAtomFrame(0);

		if (!ts.ReadXYZ(a,true,&da)) {
			m_IF.eprintf("\nError: Could not read frame %d from XYZ trajectory.\n",i+1);
			fclose(a);
			return false;
		}
		m = 0;
		for (z=0;z<(int)afr->m_oaAtoms.size();z++) {
//			m_IF.printf("@  %10.4f  %10.4f  %10.4f  vs  %10.4f  %10.4f  %10.4f\n",afr->m_oaAtoms[z]->m_fCoord[0]*LEN_AU2PM,afr->m_oaAtoms[z]->m_fCoord[1]*LEN_AU2PM,afr->m_oaAtoms[z]->m_fCoord[2]*LEN_AU2PM,da[z][0],da[z][1],da[z][2]);
			for (z2=0;z2<3;z2++) {
				t1 = da[z][z2];
				t2 = afr->m_oaAtoms[z]->m_fCoord[z2]*LEN_AU2PM;
				if (wrap) {
					while (t1 < 0)
						t1 += cv[z2];
					while (t1 >= cv[z2])
						t1 -= cv[z2];
					while (t2 < 0)
						t2 += cv[z2];
					while (t2 >= cv[z2])
						t2 -= cv[z2];
				}
				tv[z2] = t2 - t1;
			}
			t1 = sqrt(tv[0]*tv[0] + tv[1]*tv[1] + tv[2]*tv[2]);
			if (t1 > m)
				m = t1;
		}

		m_IF.printf("Max. dev. %11.5f pm  -->  ",m);

		if (m >= 0.001) {
			fail = true;
			m_IF.printf(RED,"fail\n");
		} else
			m_IF.printf(GREEN,"match\n");

		i++;
	}

	m_IF.printf("Compared %d frames.\n",i);
	m_IF.printf("Result: ");
	if (fail)
		m_IF.printf(RED,"Fail.\n");
	else
		m_IF.printf(GREEN,"Success.\n");

	fclose(a);
	Close();

	delete g_pCCEngine;

	return (!fail);*/

	return false;
}


bool CBQBTrajectoryFrameColumn::ReadColumn(int ac, CBQBBitSet *bs) {

	unsigned char uc;
	int z, z2;


	m_iType = (unsigned char)bs->ReadBitsInteger(8);

	uc = (unsigned char)bs->ReadBitsInteger(8);
	m_sLabel.resize(uc);
	for (z=0;z<(int)uc;z++)
		m_sLabel[z] = (char)bs->ReadBitsInteger(8);

	switch(m_iType) {
		case BQB_TYPE_STRING:
			m_aString.resize(ac);
			for (z=0;z<ac;z++) {
				uc = (unsigned char)bs->ReadBitsInteger(8);
				m_aString[z].resize(uc+1);
				//m_aString[z].SetBufSize(uc+1);
				for (z2=0;z2<(int)uc;z2++)
					m_aString[z][z2] = (char)bs->ReadBitsInteger(8);
					//m_aString[z](z2) = (char)bs->ReadBitsInteger(8);
				m_aString[z][uc] = 0;
				//m_aString[z](uc) = 0;
			}
			break;
		case BQB_TYPE_FLOAT:
			m_aReal.resize(ac);
			for (z=0;z<ac;z++)
				m_aReal[z] = bs->ReadBitsFloat();
			break;
		case BQB_TYPE_DOUBLE:
			m_aReal.resize(ac);
			for (z=0;z<ac;z++)
				m_aReal[z] = bs->ReadBitsDouble();
			break;
		case BQB_TYPE_UINT8:
			m_aUnsignedInt.resize(ac);
			for (z=0;z<ac;z++)
				m_aUnsignedInt[z] = bs->ReadBitsInteger(8);
			break;
		case BQB_TYPE_UINT16:
			for (z=0;z<ac;z++)
				m_aUnsignedInt[z] = bs->ReadBitsInteger(16);
			break;
		case BQB_TYPE_UINT32:
			for (z=0;z<ac;z++)
				m_aUnsignedInt[z] = bs->ReadBitsInteger(32);
			break;
		default:
			m_IF.eprintf("CBQBTrajectoryFrameColumn::ReadColumn(): Error: Type %u not yet implemented.\n",m_iType);
			abort();
	}

	return true;
}


void CBQBTrajectoryFrameColumn::WriteColumn(int ac, CBQBBitSet *bs) {

	unsigned char uc;
	int z, z2;


	bs->WriteBits(m_iType,8);

	uc = (unsigned char)m_sLabel.length();
	bs->WriteBits(uc,8);
	for (z=0;z<(int)m_sLabel.length();z++) {
		uc = m_sLabel[z];
		bs->WriteBits(uc,8);
	}

	switch(m_iType) {
		case BQB_TYPE_STRING:
			for (z=0;z<ac;z++) {
				uc = (unsigned char)m_aString[z].length();
				//uc = (unsigned char)m_aString[z].GetLength();
				bs->WriteBits(uc,8);
				for (z2=0;z2<(int)m_aString[z].length();z2++) {
				//for (z2=0;z2<(int)m_aString[z].GetLength();z2++) {
					uc = m_aString[z][z2];
					bs->WriteBits(uc,8);
				}
			}
			break;
		case BQB_TYPE_FLOAT:
			for (z=0;z<ac;z++)
				bs->WriteBitsFloat((float)m_aReal[z]);
			break;
		case BQB_TYPE_DOUBLE:
			for (z=0;z<ac;z++)
				bs->WriteBitsDouble(m_aReal[z]);
			break;
		case BQB_TYPE_UINT8:
			for (z=0;z<ac;z++)
				bs->WriteBits(m_aUnsignedInt[z],8);
			break;
		case BQB_TYPE_UINT16:
			for (z=0;z<ac;z++)
				bs->WriteBits(m_aUnsignedInt[z],16);
			break;
		case BQB_TYPE_UINT32:
			for (z=0;z<ac;z++)
				bs->WriteBits(m_aUnsignedInt[z],32);
			break;
		default:
			m_IF.eprintf("CBQBTrajectoryFrameColumn::WriteColumn(): Error: Type %u not yet implemented.\n",m_iType);
			abort();
	}
}


bool CBQBTrajectoryFrame::ReadFrame(const std::vector<unsigned char> *data) {

	CBQBBitSet bs(m_IF);
	int i, z;
	CBQBTrajectoryFrameColumn *col;


	bs.m_iaData.assign(data->begin(),data->end());

	m_iAtomCount = bs.ReadBitsInteger(32);

	if (bs.ReadBitsInteger(8) != 0) {
		if (m_pCellMatrix == NULL)
			m_pCellMatrix = new CBQBDMatrix3();
		for (z=0;z<9;z++)
			m_pCellMatrix->GetAt(z) = bs.ReadBitsDouble();
	} else {
		if (m_pCellMatrix != NULL) {
			delete m_pCellMatrix;
			m_pCellMatrix = NULL;
		}
	}

	i = bs.ReadBitsInteger(32);

	for (z=0;z<(int)m_oaColumns.size();z++)
		if (m_oaColumns[z] != NULL)
			delete m_oaColumns[z];

	m_oaColumns.resize(i);
	for (z=0;z<i;z++) {
		col = new CBQBTrajectoryFrameColumn(m_IF);
		m_oaColumns[z] = col;
		col->ReadColumn(m_iAtomCount,&bs);
	}

	return true;
}


void CBQBTrajectoryFrame::WriteFrame(std::vector<unsigned char> *data) {

	CBQBBitSet bs(m_IF);
	unsigned long ul;
	int z;


	bs.WriteBits(m_iAtomCount,32);

	if (m_pCellMatrix != NULL) {
		bs.WriteBits(1,8);
		for (z=0;z<9;z++)
			bs.WriteBitsDouble(m_pCellMatrix->GetAt(z));
	} else
		bs.WriteBits(0,8);

	ul = (unsigned long)m_oaColumns.size();

	bs.WriteBits(ul,32);

	for (z=0;z<(int)m_oaColumns.size();z++)
		m_oaColumns[z]->WriteColumn(m_iAtomCount,&bs);

	data->insert(data->end(),bs.m_iaData.begin(),bs.m_iaData.end());
}


CBQBTrajectoryFrameColumn* CBQBTrajectoryFrame::GetColumn(std::string label) {

	int z;

	for (z=0;z<(int)m_oaColumns.size();z++)
		if (bqb_strcmp_nocase(label.c_str(),m_oaColumns[z]->m_sLabel.c_str()) == 0)
			return m_oaColumns[z];

	return NULL;
}


CBQBTrajectoryFrameColumn* CBQBTrajectoryFrame::AddColumn(int type, std::string label) {

	CBQBTrajectoryFrameColumn *col;

	col = new CBQBTrajectoryFrameColumn(m_IF,(unsigned char)type,label);
	m_oaColumns.push_back(col);

	return col;
}


bool CBQBFile::SeekFrame(int i) {

	int z, k;
	unsigned long ul;
	CBQBListEntry *le;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::SeekFrame() i=%d >>>\n",i);

	if (!m_bOpenRead) {
		m_IF.eprintf("CBQBFile::SeekFrame(): Error: File not open for reading.\n");
		abort();
	}

	if (m_bListFile) {

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Is a list file.\n");

		if (m_oaBQBList.size() == 0) {
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.eprintf("CBQBFile::SeekFrame(): Error: No entries in list.\n");
			return false;
		}

		if (m_iListIndex != -1) {
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("    List entry %d is currently open; closing.\n",m_iListIndex+1);
			if (!m_oaBQBList[m_iListIndex]->m_pFile->Close()) {
				m_IF.eprintf("CBQBFile::SeekFrame(): Error: Failed to close list entry %d.\n",
					m_iListIndex+1);
				return false;
			}
		} else if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    No list entry is currently open.\n");

		m_pShortFrame = NULL; // Deleted inside of the corresponding list entry

		k = 0;
		for (z=0;z<(int)m_oaBQBList.size();z++) {
			if (m_IF.IsPL(BQB_PL_VERBOSE))
				m_IF.printf("    Processing list entry %d...\n",z+1);
			le = m_oaBQBList[z];
			if (le->m_pIndex == NULL) {
				m_IF.eprintf("CBQBFile::SeekFrame(): Error: List entry %d does not contain an index.\n",
					z+1);
				return false;
			}
			if (i >= (int)le->m_pIndex->m_iaFrameLengths.size()+k) {
				k += (int)le->m_pIndex->m_iaFrameLengths.size();
				if (m_IF.IsPL(BQB_PL_VERBOSE))
					m_IF.printf("      Entry has %lu frames, so end of entry is frame %d, which is < %d.\n",
						(unsigned long)le->m_pIndex->m_iaFrameLengths.size(),k,i);
				continue;
			} else {
				if (m_IF.IsPL(BQB_PL_VERBOSE))
					m_IF.printf("      Entry has %lu frames, so end of entry is frame %d, which is >= %d --> Found.\n",
						(unsigned long)le->m_pIndex->m_iaFrameLengths.size(),(int)le->m_pIndex->m_iaFrameLengths.size()+k,i);
				m_iListIndex = z;
				if (!m_oaBQBList[m_iListIndex]->m_pFile->OpenRead(m_oaBQBList[m_iListIndex]->m_sFileName)) {
					m_IF.eprintf("CBQBFile::SeekFrame(): Error: Could not open list entry %d for reading.\n",
						m_iListIndex+1);
					return false;
				}
				if (i-k != 0) {
					if (m_IF.IsPL(BQB_PL_VERBOSE))
						m_IF.printf("      Seeking to frame %d of list entry %d...\n",
							i-k+1,m_iListIndex+1);
					if (!m_oaBQBList[m_iListIndex]->m_pFile->SeekFrame(i-k)) {
						m_IF.eprintf("CBQBFile::SeekFrame(): Error: Could not seek to frame %d within list entry %d.\n",
							i-k,m_iListIndex+1);
						return false;
					}
				} else if (m_IF.IsPL(BQB_PL_VERBOSE))
					m_IF.printf("      Required frame is first frame of list entry %d, not seeking.\n",
						m_iListIndex+1);

				break;
			}
		}

	} else {

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Is a simple file.\n");

		if (m_pFile == NULL) {
			m_IF.eprintf("CBQBFile::SeekFrame(): Error: File pointer is NULL.\n");
			abort();
		}

		if (m_pIndex == NULL) {
			m_IF.eprintf("CBQBFile::SeekFrame(): Error: File does not contain an index.\n");
			abort();
		}

		if (i >= (int)m_pIndex->m_iaFrameLengths.size()) {
			m_IF.eprintf("CBQBFile::SeekFrame(): Error: Tried to seek frame beyond index (%d/%lu).\n",
				i,(unsigned long)m_pIndex->m_iaFrameLengths.size());
			abort();
		}

		ul = 0;
		for (z=0;z<i;z++)
			ul += (unsigned long)m_pIndex->m_iaFrameLengths[z];

		if (m_IF.IsPL(BQB_PL_VERBOSE))
			m_IF.printf("    Offset of frame %d is %lu.\n",i,ul);

		fseek(m_pFile,ul,SEEK_SET);
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("<<< CBQBFile::SeekFrame() i=%d <<<\n",i);

	return true;
}


bool CBQBFile::SkipFrames(int i) {

	int z;


	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf(">>> CBQBFile::SkipFrames() i=%d >>>\n",i);

	if (!m_bOpenRead) {
		m_IF.eprintf("CBQBFile::SkipFrames(): Error: File not open for reading.\n");
		abort();
	}

	for (z=0;z<i;z++) {
		if (!ReadFrame()) {
			m_IF.eprintf("CBQBFile::SkipFrames(): Error while skipping frame %d/%d.\n",
				z+1,i);
			return false;
		}
	}

	if (m_IF.IsPL(BQB_PL_VERBOSE))
		m_IF.printf("<<< CBQBFile::SkipFrames() i=%d <<<\n",i);

	return true;
}


