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



#ifndef BQB_INTERFACE_H
#define BQB_INTERFACE_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_math.h"
#include "bqb_parmset.h"
#include "bqb_format.h"
#include "bqb_cubeframe.h"

#include <vector>
#include <set>



#define BQB_PL_SILENT     0  // Do not print anything on the screen
#define BQB_PL_QUIET      1
#define BQB_PL_STANDARD   2  // As in TRAVIS / bqbtool
#define BQB_PL_VERBOSE    3  // Print all information from Front-Ends
#define BQB_PL_DEBUG      4  // Also print all information from Compressor


#define BQB_OK                       0
#define BQB_ERROR_FILE             100
#define BQB_ERROR_COMPRESS_FAILED  200
#define BQB_ERROR_PARAMETERS       300


#define BQB_UNIT_BOHR     1


#define BQB_ORDER_ZYX     1
#define BQB_ORDER_XYZ     2


#define BQB_CONTENTS_ATOMPOS      1
#define BQB_CONTENTS_VOLUMETRIC   2
#define BQB_CONTENTS_CELLINFO     4


#define BQB_OVERWRITE   1
#define BQB_APPEND      2



class CBQBDriver;
class CBQBInterface;
class CBQBEngine;
class CBQBFileWriter;
class CBQBFileReader;
class CBQBVolHeader;
class CBQBAtomHeader;
class CBQBCellInfo;
class CBQBVolFrame;
class CBQBAtomPosFrame;



CBQBInterface* BQBCreateInterface(unsigned long flags);


bool BQBReleaseInterface(CBQBInterface *i);


unsigned long BQBGetLastError();



class CBQBInterface {

	friend CBQBInterface* BQBCreateInterface(unsigned long flags);
	friend bool BQBReleaseInterface(CBQBInterface *i);
	friend class CBQBFileWriter;
	friend class CBQBFileReader;
	friend class CBQBVolHeader;
	friend class CBQBAtomHeader;
	friend class CBQBCellInfo;
	friend class CBQBVolFrame;
	friend class CBQBAtomPosFrame;

public:
	CBQBDriver* CreateDriver(unsigned long flags);

	bool DestroyDriver(CBQBDriver *d);

	CBQBEngine* CreateEngine(unsigned long flags);

	bool DestroyEngine(CBQBEngine *e);

	void SetPrintCallback( void (*fp)(const char *) );
	void SetPrintCallback( void (*fp)(const char *, ...) );
	void SetPrintFile( FILE *f );
	void ResetPrint();

	void SetEPrintCallback( void (*fp)(const char *) );
	void SetEPrintCallback( void (*fp)(const char *, ...) );
	void SetEPrintFile( FILE *f );
	void ResetEPrint();

	void SetBPrintCallback( void (*fp)(const char *) );
	void SetBPrintCallback( void (*fp)(const char *, ...) );
	void SetBPrintFile( FILE *f );
	void ResetBPrint();

	void FlushLog();


	bool WriteCubeFrame(
		FILE *outfile,
		const CBQBVolHeader *volhead,
		const CBQBVolFrame *volframe,
		const CBQBAtomHeader *atomhead,
		const CBQBAtomPosFrame *atomframe,
		const CBQBCellInfo *cellinfo,
		unsigned long flags=0
	);


	CBQBFileWriter* CreateBQBFileWriter( const char *name, int history, unsigned long contents, unsigned long flags=0 );

	CBQBFileReader* CreateBQBFileReader( const char *name, unsigned long flags=0 );

	CBQBVolHeader* CreateVolumetricHeader();

	CBQBAtomHeader* CreateAtomHeader();

	CBQBCellInfo* CreateCellInfo();

	CBQBVolFrame* CreateVolumetricFrame( const CBQBVolHeader *header, int id=-1 );

	CBQBAtomPosFrame* CreateAtomPosFrame( const CBQBAtomHeader *header );

	void RecycleVolumetricFrame( CBQBVolFrame *frame ) {
		m_oaRecycledVolFrames.push_back( frame );
	}

	void RecycleAtomPositionFrame( CBQBAtomPosFrame *frame ) {
		m_oaRecycledAtomPosFrames.push_back( frame );
	}


	#ifdef __GNUG__  // Variadic Argument Type Checking of GCC

		// Note the implicit first "this" argument!

		void printf(const char *s, ...)  const __attribute__ ((format (printf, 2, 3)));

		void bprintf(const char *s, ...) const __attribute__ ((format (printf, 2, 3)));

		void eprintf(const char *s, ...) const __attribute__ ((format (printf, 2, 3)));

	#else

		void printf(const char *s, ...)  const;

		void bprintf(const char *s, ...) const;

		void eprintf(const char *s, ...) const;

	#endif

	void SetPrintLevel(int i);

	int GetPrintLevel() const;

	bool IsPL(int i) const { return (m_iPrintLevel >= i); }

	unsigned long GetLastError() const { return m_iLastError; }


private:

	CBQBInterface();

	~CBQBInterface();

	void SetLastError( unsigned long e );

	int m_iPrintLevel;

	unsigned long m_iLastError;

	void (*m_pPrintCallback)(const char *);
	void (*m_pEPrintCallback)(const char *);
	void (*m_pBPrintCallback)(const char *);

	void (*m_pPrintCallbackVar)(const char *, ...);
	void (*m_pEPrintCallbackVar)(const char *, ...);
	void (*m_pBPrintCallbackVar)(const char *, ...);

	FILE *m_pPrintFile;
	FILE *m_pBPrintFile;
	FILE *m_pEPrintFile;

	std::set<CBQBFileWriter*> m_oaListFileWriter;
	std::set<CBQBFileReader*> m_oaListFileReader;
	std::set<CBQBVolHeader*> m_oaListVolHeader;
	std::set<CBQBAtomHeader*> m_oaListAtomHeader;
	std::set<CBQBCellInfo*> m_oaListCellInfo;
	std::set<CBQBVolFrame*> m_oaListVolFrame;
	std::set<CBQBAtomPosFrame*> m_oaListAtomPosFrame;

	std::vector<CBQBVolFrame*> m_oaRecycledVolFrames;
	std::vector<CBQBAtomPosFrame*> m_oaRecycledAtomPosFrames;

	bool m_bActive;
};



class CBQBFileWriter {

	friend class CBQBInterface;
	friend CBQBInterface* BQBCreateInterface( unsigned long flags );
	friend bool BQBReleaseInterface( CBQBInterface *i );

public:

	void CloseFile();

	void SetVolumetricHeader( CBQBVolHeader *header );

	void SetAtomHeader( CBQBAtomHeader *header );

	void SetCellInfo( CBQBCellInfo *info );

	const CBQBVolHeader* GetVolumetricHeader() const;

	const CBQBAtomHeader* GetAtomHeader() const;

	const CBQBCellInfo* GetCellInfo() const;

	void PushVolumetricFrame( CBQBVolFrame *frame );

	void PushAtomPosFrame( CBQBAtomPosFrame *frame );

	void PushCommentLine( const char *line1, const char *line2=NULL );

	void Optimize( int mode );

	bool WritePendingFrames( bool check );

	bool WriteIndex();

	bool HavePendingVolumetricFrames() const;

	int GetPendingVolumetricFrameCount() const;

	bool HavePendingAtomPosFrames() const;

	int GetPendingAtomPosFrameCount() const;

	bool HavePendingCommentLines() const;

	int GetPendingCommentLineCount() const;


	CBQBParameterSet_PosAndVol *m_pParmPosVol;


private:
	explicit CBQBFileWriter( CBQBInterface &i ) : m_pParmPosVol(NULL), m_bContainsAtomPosition(false), m_bContainsVolumetricData(false),
		m_bContainsCellInfo(false), m_pVolHeader(NULL), m_pAtomHeader(NULL), m_pCellInfo(NULL), m_pEngine(NULL), m_pBQBFile(NULL), m_iHistoryDepth(0),
		m_iReadPos(0), m_iVolumetricWritePos(0), m_iAtomPosWritePos(0), m_iCommentLineWritePos(0), m_bStartLock(false), m_bAlreadyPushed(false),
		m_iFramesWritten(0), m_iVolumetricFramesPending(0), m_iAtomPosFramesPending(0), m_iCommentLinesPending(0), m_bActive(true), m_IF(i) {
	}

	~CBQBFileWriter();

	bool OpenFileWrite( const char *s, int history, unsigned long contents, unsigned long flags );

	bool m_bContainsAtomPosition;
	bool m_bContainsVolumetricData;
	bool m_bContainsCellInfo;

	CBQBVolHeader *m_pVolHeader;
	CBQBAtomHeader *m_pAtomHeader;
	CBQBCellInfo *m_pCellInfo;

	CBQBEngine *m_pEngine;
	CBQBFile *m_pBQBFile;

	int m_iHistoryDepth;

	int m_iReadPos;

	int m_iVolumetricWritePos;
	std::vector<CBQBVolFrame*> m_oaVolumetricFrameHistory;

	int m_iAtomPosWritePos;
	std::vector<CBQBAtomPosFrame*> m_oaAtomPosFrameHistory;

	int m_iCommentLineWritePos;
	std::vector<std::string> m_saCommentLine1;
	std::vector<std::string> m_saCommentLine2;

	bool m_bStartLock;

	bool m_bAlreadyPushed;

	int m_iFramesWritten;
	
	int m_iVolumetricFramesPending;
	int m_iAtomPosFramesPending;
	int m_iCommentLinesPending;

	int m_iVolMaxOrder;
	int m_iAtomPosMaxOrder;
	int m_iUseHistory;

	bool m_bActive;
	CBQBInterface &m_IF;
};



class CBQBFileReader {

	friend class CBQBInterface;
	friend CBQBInterface* BQBCreateInterface( unsigned long flags );
	friend bool BQBReleaseInterface( CBQBInterface *i );

public:

	void CloseFile();


	bool ContainsAtomPositions() const {
		return m_bContainsAtomPosition;
	}


	bool ContainsVolumetricData() const {
		return m_bContainsVolumetricData;
	}


	bool ContainsCellInfo() const {
		return m_bContainsCellInfo;
	}


	bool HasIndex() const {
		return false;
	}


	unsigned long GetTotalFrameCount() const {
		return 0;
	}


	bool ReadFrame( unsigned long flags=0 );


	bool IsEOF() const {
		return false;
	}


	const CBQBVolHeader* GetVolumetricHeader() const;

	const CBQBAtomHeader* GetAtomHeader() const;

	const CBQBCellInfo* GetCellInfo() const;

	const CBQBVolFrame* GetCurrentVolumetricFrame() const;

	const CBQBAtomPosFrame* GetCurrentAtomPosFrame() const;


private:
	explicit CBQBFileReader( CBQBInterface &i ) : m_bContainsAtomPosition(false), m_bContainsVolumetricData(false), m_bContainsCellInfo(false),
		m_pVolHeader(NULL), m_pAtomHeader(NULL), m_pCellInfo(NULL), m_bActive(true), m_IF(i) {
	}

	~CBQBFileReader();

	bool OpenFileRead( const char *s, unsigned long flags );

	bool m_bContainsAtomPosition;
	bool m_bContainsVolumetricData;
	bool m_bContainsCellInfo;

	const CBQBVolHeader *m_pVolHeader;
	const CBQBAtomHeader *m_pAtomHeader;
	const CBQBCellInfo *m_pCellInfo;

	bool m_bActive;
	CBQBInterface &m_IF;
};



class CBQBVolHeader {

	friend class CBQBInterface;
	friend class CBQBFileWriter;

public:

	void SetResolution( unsigned int x, unsigned int y, unsigned int z ) {
		if (m_bAssigned) {
			m_IF.eprintf("CBQBVolHeader::SetResolution(): Error: Cannot modify cell info after it has been used. Create new CBQBVolHeader instance instead.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_iRes[0] = x;
		m_iRes[1] = y;
		m_iRes[2] = z;
		m_iResXY = x*y;
		m_iResXYZ = (unsigned long)x * y * z;
	}


	unsigned int GetResolutionX() const {
		return m_iRes[0];
	}


	unsigned int GetResolutionY() const {
		return m_iRes[1];
	}


	unsigned int GetResolutionZ() const {
		return m_iRes[2];
	}


	unsigned int GetResolutionXY() const {
		return m_iResXY;
	}


	unsigned long GetResolutionXYZ() const {
		return m_iResXYZ;
	}


	unsigned long GetSigni() const {
		return m_iSigni;
	}


	unsigned long GetEps() const {
		return m_iEps;
	}


	void SetSigni( int signi ) {
		if (m_bAssigned) {
			m_IF.eprintf("CBQBVolHeader::SetSigni(): Error: Cannot modify cell info after it has been used. Create new CBQBVolHeader instance instead.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_iSigni = signi;
	}


	void SetEps( int eps ) {
		if (m_bAssigned) {
			m_IF.eprintf("CBQBVolHeader::SetEps(): Error: Cannot modify cell info after it has been used. Create new CBQBVolHeader instance instead.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_iEps = eps;
	}


private:
	explicit CBQBVolHeader( CBQBInterface &i ) : m_iResXY(0), m_iResXYZ(0), m_iSigni(0), m_iEps(0), m_bAssigned(false), m_bActive(true), m_IF(i) {
		m_iRes[0] = 0;
		m_iRes[1] = 0;
		m_iRes[2] = 0;
		m_faStrideA[0] = 0;
		m_faStrideA[1] = 0;
		m_faStrideA[2] = 0;
		m_faStrideB[0] = 0;
		m_faStrideB[1] = 0;
		m_faStrideB[2] = 0;
		m_faStrideC[0] = 0;
		m_faStrideC[1] = 0;
		m_faStrideC[2] = 0;
		m_faCenter[0] = 0;
		m_faCenter[1] = 0;
		m_faCenter[2] = 0;
		m_iaStrideA[0] = 0;
		m_iaStrideA[1] = 0;
		m_iaStrideA[2] = 0;
		m_iaStrideB[0] = 0;
		m_iaStrideB[1] = 0;
		m_iaStrideB[2] = 0;
		m_iaStrideC[0] = 0;
		m_iaStrideC[1] = 0;
		m_iaStrideC[2] = 0;
		m_iaCenter[0] = 0;
		m_iaCenter[1] = 0;
		m_iaCenter[2] = 0;
	}

	~CBQBVolHeader();

	unsigned int m_iRes[3];
	unsigned int m_iResXY;
	unsigned long m_iResXYZ;

	double m_faStrideA[3];
	double m_faStrideB[3];
	double m_faStrideC[3];
	double m_faCenter[3];

	long m_iaStrideA[3];
	long m_iaStrideB[3];
	long m_iaStrideC[3];
	long m_iaCenter[3];

	int m_iSigni;
	int m_iEps;

	bool m_bAssigned;

	bool m_bActive;
	CBQBInterface &m_IF;
};



class CBQBAtomHeader {

	friend class CBQBInterface;
	friend class CBQBFileWriter;

public:

	void SetAtomCount( unsigned int i ) {
		if (m_bAssigned) {
			m_IF.eprintf("CBQBAtomHeader::SetAtomCount(): Error: Cannot modify cell info after it has been used. Create new CBQBAtomHeader instance instead.\n");
			m_IF.eprintf("\n");
			abort();
		}
		unsigned int z;
		m_iAtomCount = i;
		m_iaOrd.resize( i );
		m_saLabel.resize( i );
		for (z=0;z<i;z++)
			m_iaOrd[ z ] = 0;
	}


	unsigned int GetAtomCount() const {
		return m_iAtomCount;
	}


	void SetAtomOrd( unsigned int i, unsigned char ord ) {
		if (m_bAssigned) {
			m_IF.eprintf("CBQBAtomHeader::SetAtomOrd(): Error: Cannot modify cell info after it has been used. Create new CBQBAtomHeader instance instead.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_iaOrd[ i ] = ord;
		m_saLabel[ i ] = GetAtomOrdLabel( ord );
	}


	void SetAtomLabel( unsigned int i, const char *s ) {
		if (m_bAssigned) {
			m_IF.eprintf("CBQBAtomHeader::SetAtomLabel(): Error: Cannot modify cell info after it has been used. Create new CBQBAtomHeader instance instead.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_iaOrd[ i ] = (unsigned char)GetAtomOrd( s );
		m_saLabel[ i ] = s;
	}


	unsigned long GetPrecision() const {
		return m_iPrecision;
	}


	void SetPrecision( int prec ) {
		if (m_bAssigned) {
			m_IF.eprintf("CBQBAtomHeader::SetPrecision(): Error: Cannot modify cell info after it has been used. Create new CBQBAtomHeader instance instead.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_iPrecision = prec;
	}


private:
	explicit CBQBAtomHeader( CBQBInterface &i ) : m_iAtomCount(0), m_iPrecision(0), m_bAssigned(false), m_bActive(true), m_IF(i) {
	}

	~CBQBAtomHeader();

	unsigned int m_iAtomCount;
	std::vector<unsigned char> m_iaOrd;
	std::vector<std::string> m_saLabel;

	int m_iPrecision;

	bool m_bAssigned;

	bool m_bActive;
	CBQBInterface &m_IF;
};



class CBQBCellInfo {

	friend class CBQBInterface;
	friend class CBQBFileWriter;

public:

	void SetCubic( double val, int unit ) {
		if (m_bAssigned) {
			m_IF.eprintf("CBQBCellInfo::SetCubic(): Error: Cannot modify cell info after it has been used. Create new CBQBCellInfo instance instead.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_faValues.resize( 1 );
		switch( unit ) {
			case BQB_UNIT_BOHR:
				m_faValues[ 0 ] = val;
				break;
			default:
				m_IF.eprintf( "CBQBCellInfo::SetCubic(): Error: Unknown unit %d.\n", unit );
				abort();
		}
		m_bOrthorhombic = true;
		m_bCubic = true;
	}


	void SetOrthorhombic( double x, double y, double z, int unit ) {
		if (m_bAssigned) {
			m_IF.eprintf("CBQBCellInfo::SetOrthorhombic(): Error: Cannot modify cell info after it has been used. Create new CBQBCellInfo instance instead.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_faValues.resize( 3 );
		switch( unit ) {
			case BQB_UNIT_BOHR:
				m_faValues[ 0 ] = x;
				m_faValues[ 1 ] = y;
				m_faValues[ 2 ] = z;
				break;
			default:
				m_IF.eprintf( "CBQBCellInfo::SetOrthorhombic(): Error: Unknown unit %d.\n", unit );
				abort();
		}
		m_bOrthorhombic = true;
		m_bCubic = false;
	}


	void SetGeneric( double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz, int unit ) {
		if (m_bAssigned) {
			m_IF.eprintf("CBQBCellInfo::SetGeneric(): Error: Cannot modify cell info after it has been used. Create new CBQBCellInfo instance instead.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_faValues.resize( 9 );
		switch( unit ) {
			case BQB_UNIT_BOHR:
				m_faValues[ 0 ] = ax;
				m_faValues[ 1 ] = ay;
				m_faValues[ 2 ] = az;
				m_faValues[ 3 ] = bx;
				m_faValues[ 4 ] = by;
				m_faValues[ 5 ] = bz;
				m_faValues[ 6 ] = cx;
				m_faValues[ 7 ] = cy;
				m_faValues[ 8 ] = cz;
				break;
			default:
				m_IF.eprintf( "CBQBCellInfo::SetGeneric(): Error: Unknown unit %d.\n", unit );
				abort();
		}
		m_bOrthorhombic = false;
		m_bCubic = false;
	}


	bool IsCubic() const {
		return m_bCubic;
	}


	bool IsOrthorhombic() const {
		return m_bOrthorhombic;
	}


	double GetCellAX() const {
		if (m_bCubic)
			return m_faValues[0];
		else if (m_bOrthorhombic)
			return m_faValues[0];
		else
			return m_faValues[0];
	}


	double GetCellAY() const {
		if (m_bCubic)
			return 0;
		else if (m_bOrthorhombic)
			return 0;
		else
			return m_faValues[1];
	}


	double GetCellAZ() const {
		if (m_bCubic)
			return 0;
		else if (m_bOrthorhombic)
			return 0;
		else
			return m_faValues[2];
	}


	double GetCellBX() const {
		if (m_bCubic)
			return 0;
		else if (m_bOrthorhombic)
			return 0;
		else
			return m_faValues[3];
	}


	double GetCellBY() const {
		if (m_bCubic)
			return m_faValues[0];
		else if (m_bOrthorhombic)
			return m_faValues[1];
		else
			return m_faValues[4];
	}


	double GetCellBZ() const {
		if (m_bCubic)
			return 0;
		else if (m_bOrthorhombic)
			return 0;
		else
			return m_faValues[5];
	}


	double GetCellCX() const {
		if (m_bCubic)
			return 0;
		else if (m_bOrthorhombic)
			return 0;
		else
			return m_faValues[6];
	}


	double GetCellCY() const {
		if (m_bCubic)
			return 0;
		else if (m_bOrthorhombic)
			return 0;
		else
			return m_faValues[7];
	}


	double GetCellCZ() const {
		if (m_bCubic)
			return m_faValues[0];
		else if (m_bOrthorhombic)
			return m_faValues[2];
		else
			return m_faValues[8];
	}


private:
	explicit CBQBCellInfo( CBQBInterface &i ) : m_bCubic(false), m_bOrthorhombic(false), m_bAssigned(false), m_bActive(true), m_IF(i) {
	}

	~CBQBCellInfo();

	std::vector<double> m_faValues;
	bool m_bCubic;
	bool m_bOrthorhombic;

	bool m_bAssigned;

	bool m_bActive;
	CBQBInterface &m_IF;
};



class CBQBVolFrame {

	friend class CBQBInterface;
	friend class CBQBFileWriter;

public:

	unsigned int GetResolutionX() const {
		return m_oCubeFrame.m_iRes[0];
	}


	unsigned int GetResolutionY() const {
		return m_oCubeFrame.m_iRes[1];
	}


	unsigned int GetResolutionZ() const {
		return m_oCubeFrame.m_iRes[2];
	}


	double GetVoxel( unsigned int ix, unsigned int iy, unsigned int iz ) const {
		return m_oCubeFrame.m_faBin[ iz*m_oCubeFrame.m_iResXY + iy*m_oCubeFrame.m_iRes[0] + ix ];
	}


	void SetVoxel( unsigned int ix, unsigned int iy, unsigned int iz, double val) {
		m_oCubeFrame.m_faBin[ iz*m_oCubeFrame.m_iResXY + iy*m_oCubeFrame.m_iRes[0] + ix ] = val;
	}


	void LoadFromArray( const double *dp, int mode );


	void StoreToArray( double *dp, int mode ) const;


private:
	explicit CBQBVolFrame( CBQBInterface &i ) : m_pVolHeader(NULL), m_oCubeFrame(i), m_bActive(true), m_IF(i) {
	}

	~CBQBVolFrame();

	const CBQBVolHeader *m_pVolHeader;

	CBQBCubeFrame m_oCubeFrame;

//	std::vector<double> m_faData;

	bool m_bActive;
	CBQBInterface &m_IF;
};



class CBQBAtomPosFrame {

	friend class CBQBInterface;
	friend class CBQBFileWriter;

public:

	void SetAtomPosition( unsigned int i, double px, double py, double pz, int unit ) {
		switch( unit ) {
			case BQB_UNIT_BOHR:
				m_oAtomSet.m_oaAtoms[ i ]->m_fCoord[ 0 ] = px;
				m_oAtomSet.m_oaAtoms[ i ]->m_fCoord[ 1 ] = py;
				m_oAtomSet.m_oaAtoms[ i ]->m_fCoord[ 2 ] = pz;
				break;
			default:
				m_IF.eprintf( "CBQBAtomPosFrame::SetAtomPosition(): Error: Unknown unit %d.\n", unit );
				abort();
		}
	}


	unsigned int GetAtomCount() const {
		return m_pAtomHeader->GetAtomCount();
	}


private:
	explicit CBQBAtomPosFrame( CBQBInterface &i ) : m_pAtomHeader(NULL), m_oAtomSet(i), m_bActive(true), m_IF(i) {
	}

	~CBQBAtomPosFrame();

	const CBQBAtomHeader *m_pAtomHeader;

	//std::vector<CBQBDVector3> m_vaCoord;

	CBQBAtomSet m_oAtomSet;

	bool m_bActive;
	CBQBInterface &m_IF;
};




#endif


