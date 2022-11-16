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



#ifndef BQB_ENGINE_H
#define BQB_ENGINE_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_tools.h"
#include <stdio.h>
#include "bqb_cubeframe.h"
#include "bqb_hufftree.h"
#include "bqb_extrapolator.h"
#include "bqb_parmset.h"
#include "bqb_interface.h"



class CBQBInterface;



class CBQBOldCellInfo {
public:

	CBQBOldCellInfo() : m_mCell(0) {
	}


	void SetCellInfo(double x, double y, double z) {
		m_mCell.Unity();
		m_mCell(0,0) = x;
		m_mCell(1,1) = y;
		m_mCell(2,2) = z;
	}


	void SetCellInfo(double ax, double ay, double az, double bx, double by, double bz, double cx, double cy, double cz) {
		m_mCell(0,0) = ax;
		m_mCell(0,1) = ay;
		m_mCell(0,2) = az;
		m_mCell(1,0) = bx;
		m_mCell(1,1) = by;
		m_mCell(1,2) = bz;
		m_mCell(2,0) = cx;
		m_mCell(2,1) = cy;
		m_mCell(2,2) = cz;
	}


	bool IsCubic() const {
		if ((m_mCell(0,0) == m_mCell(1,1)) && (m_mCell(0,0) == m_mCell(2,2)) &&
			(m_mCell(0,1) == 0) && (m_mCell(0,2) == 0) && (m_mCell(1,0) == 0) &&
			(m_mCell(1,2) == 0) && (m_mCell(2,0) == 0) && (m_mCell(2,1) == 0))
			return true;
		else
			return false;
	}


	bool IsOrthorhombic() const {
		if ((m_mCell(0,1) == 0) && (m_mCell(0,2) == 0) && (m_mCell(1,0) == 0) &&
			(m_mCell(1,2) == 0) && (m_mCell(2,0) == 0) && (m_mCell(2,1) == 0))
			return true;
		else
			return false;
	}


	CBQBDMatrix3 m_mCell;
};



/********************************************************************************************************************

template<class FT> class CBQBCache {
public:

	CBQBCache(CBQBInterface &i) : m_IF(i), m_iHistoryDepth(0), m_iCachePos(0), m_iReadPos(0),
		m_bStartLock(false) {
	}

	virtual ~CBQBCache() {
		m_oaFrames.clear();
	}

	void SetStartLock() {
		if (m_bStartLock) {
			m_IF.eprintf("CBQBCache::SetStartLock(): Error: Start lock already set.\n");
			abort();
		}
		m_bStartLock = true;
	}

	void LiftStartLock() {
		if (!m_bStartLock) {
			m_IF.eprintf("CBQBCache::LiftStartLock(): Error: Start lock not in place.\n");
			abort();
		}
		m_bStartLock = false;
	}

	void SetHistoryDepth(int i) {
		m_oaFrames.resize(i,NULL);
	}

	int GetCachedSteps() const {
		if (m_iCachePos >= m_iReadPos)
			return m_iCachePos - m_iReadPos;
		else
			return m_iCachePos - m_iReadPos + (int)m_oaFrames.size();
	}

	const FT* GetNewFrameHistory(int i) const {
		if (i > m_iHistoryDepth) {
			m_IF.eprintf("CBQBCache::GetFrameHistory(): Error: Index out of range (%d/%d).\n",i,m_iHistoryDepth);
			abort();
		}
		int t = m_iReadPos - i;
		if (t < 0)
			t += (int)m_oaFrames.size();
		return m_oaFrames[t];
	}

	void RewindReadPos() {
		if (!m_bStartLock) {
			m_IF.eprintf("CBQBCache::RewindReadPos(): Error: Requires start lock.\n");
			abort();
		}
		m_iReadPos = -1;
	}

	void Reset() {
		m_iCachePos = -1;
		m_iReadPos = -1;
	}

	virtual void PushFrame( FT *frame ) = 0;

protected:
	CBQBInterface &m_IF;

	int m_iHistoryDepth;
	int m_iCachePos; // Frame which has been cached last
	int m_iReadPos; // Frame which was read last
	std::vector<FT*> m_oaFrames;
	bool m_bStartLock;
};



class CBQBCubeCache : public CBQBCache<CBQBVolFrame> {
public:

	CBQBCubeCache(CBQBInterface &i)
		: CBQBCache<CBQBVolFrame>(i), m_iEps(-1), m_iCSigni(-1), m_iASigni(-1) {
	}

	~CBQBCubeCache() {
	}

	void SetReadParameters(int eps, int csigni, int asigni) {
		m_iEps = eps;
		m_iCSigni = csigni;
		m_iASigni = asigni;
	}

	void PushFrame( CBQBVolFrame *frame ) {
		if ((m_iCachePos < 0) || (m_iCachePos >= (int)m_oaFrames.size())) {
			m_IF.eprintf("CBQBCubeCache::CacheOneStep(): Error: m_iCachePos out of range (%d/%lu).\n",m_iCachePos,(unsigned long)m_oaFrames.size());
			abort();
		}
		if (m_oaFrames[m_iCachePos] != NULL)
			m_IF.RecycleVolumetricFrame( m_oaFrames[m_iCachePos] );
		m_oaFrames[m_iCachePos] = frame;
	}

private:

	int m_iEps;
	int m_iCSigni;
	int m_iASigni;
};



template<class FT> class CBQBReadCache : public CBQBCache<FT> {
public:

	CBQBReadCache(CBQBInterface &i) : CBQBCache<FT>(i), m_fRef(NULL), m_fInput(NULL) {
	}

	virtual ~CBQBReadCache() {
		if (m_fInput != NULL)
			CloseFile();
	}

	bool CacheSteps(int i, bool progress) {
		if (progress) {
			m_IF.printf("    Caching %d frames...\n",i);
			m_IF.printf("      [");
		}
		double tfs = i / 60.0;
		int k=0;
		bool b=false;
		for (int z=0;z<i;z++) {
			m_iCachePos++;
			if (m_iCachePos >= (int)m_oaFrames.size()) {
				if (m_bStartLock) {
					m_IF.eprintf("CBQBReadCache::CacheSteps(): Error: Violation of start lock.\n");
					abort();
				}
				m_iCachePos = 0;
				b = true;
			}
			if (!CacheOneStep()) {
				m_iCachePos--;
				if (b && (m_iCachePos < 0))
					m_iCachePos = (int)m_oaFrames.size()-1;
				m_IF.printf("] Only %d frames read.\n",k);
				return false;
			}
			k++;
			if (progress && (fmod(z,tfs) < 1.0))
				m_IF.printf("#");
		}
		if (progress)
			m_IF.printf("] Done.\n");
		return true;
	}

	bool FileOpenRead(const char *inp, const char *ref = NULL) {
		m_iCachePos = -1;
		m_iReadPos = -1;
		m_fInput = fopen(inp,"rb");
		if (m_fInput == NULL) {
			m_IF.eprintf("CBQBReadCache::FileOpenRead(): Error: Could not open input file for reading.\n");
			return false;
		}
		if (ref != NULL) {
			m_fRef = fopen(inp,"wb");
			if (m_fRef == NULL) {
				m_IF.eprintf("CBQBReadCache::FileOpenRead(): Error: Could not open reference file for writing.\n");
				fclose(m_fInput);
				m_fInput = NULL;
				return false;
			}
		}
		return true;
	}

	void CloseFile() {
		fclose(m_fInput);
		m_fInput = NULL;
		if (m_fRef != NULL) {
			fclose(m_fRef);
			m_fRef = NULL;
		}
	}

	long ftell() {
		if (m_fInput == NULL)
			return 0;
		return ::ftell(m_fInput);
	}

	void PushFrame( FT *frame ) {
		if ((m_iCachePos < 0) || (m_iCachePos >= (int)m_oaFrames.size())) {
			m_IF.eprintf("CBQBReadCache::PushFrame(): Error: m_iCachePos out of range (%d/%lu).\n",m_iCachePos,(unsigned long)m_oaFrames.size());
			abort();
		}
		if (m_oaFrames[m_iCachePos] != NULL)
			delete m_oaFrames[m_iCachePos];
		m_oaFrames[m_iCachePos] = frame;
	}

	const FT* GetNextFrame() {
		if (m_iReadPos == m_iCachePos) {
			m_iCachePos++;
			if (m_iCachePos >= (int)m_oaFrames.size()) {
				if (m_bStartLock) {
					m_IF.eprintf("CBQBCache::GetNextFrame(): Error: Violation of start lock.\n");
					abort();
				}
				m_iCachePos = 0;
			}
			if (!CacheOneStep())
				return NULL;
		}
		m_iReadPos++;
		if (m_iReadPos >= (int)m_oaFrames.size())
			m_iReadPos = 0;
		if (m_iHistoryDepth+1 < (int)m_oaFrames.size())
			m_iHistoryDepth++;
		return m_oaFrames[m_iReadPos];
	}

	virtual bool SkipOneStep() = 0;


protected:

	virtual bool CacheOneStep() = 0;

	FILE *m_fRef;
	FILE *m_fInput;
};



class CBQBReadXYZCache : public CBQBReadCache<CBQBAtomSet> {
public:

	CBQBReadXYZCache(CBQBInterface &i) : CBQBReadCache<CBQBAtomSet>(i), m_iSigni(-1) {
	}

	~CBQBReadXYZCache() {
	}

	void SetReadParameters(int signi) {
		m_iSigni = signi;
	}

	bool SkipOneStep() {
		return CBQBAtomSet(m_IF).SkipXYZ(m_fInput);
	}

private:
	bool CacheOneStep() {
		if ((m_iCachePos < 0) || (m_iCachePos >= (int)m_oaFrames.size())) {
			m_IF.eprintf("CBQBReadXYZCache::CacheOneStep(): Error: m_iCachePos out of range (%d/%lu).\n",m_iCachePos,(unsigned long)m_oaFrames.size());
			abort();
		}
		if (m_iSigni < 0) {
			m_IF.eprintf("CBQBReadXYZCache::CacheOneStep(): Error: m_iSigni out of range (%d).\n",m_iSigni);
			abort();
		}
		if (m_oaFrames[m_iCachePos] != NULL)
			delete m_oaFrames[m_iCachePos];
		m_oaFrames[m_iCachePos] = new CBQBAtomSet(m_IF);
		if (m_fInput == NULL) {
			m_IF.eprintf("CBQBReadXYZCache::CacheOneStep(): Error: m_fInput == NULL.\n");
			abort();
		}
		if (!m_oaFrames[m_iCachePos]->ReadXYZ(m_fInput,m_iSigni,m_fRef))
			return false;
		return true;
	}

	int m_iSigni;
};



class CBQBReadCubeCache : public CBQBReadCache<CBQBCubeFrame> {
public:

	CBQBReadCubeCache(CBQBInterface &i)
		: CBQBReadCache<CBQBCubeFrame>(i), m_iEps(-1), m_iCSigni(-1), m_iASigni(-1) {
	}

	~CBQBReadCubeCache() {
	}

	void SetReadParameters(int eps, int csigni, int asigni) {
		m_iEps = eps;
		m_iCSigni = csigni;
		m_iASigni = asigni;
	}

	bool SkipOneStep() {
		return CBQBCubeFrame(m_IF).SkipFrame(m_fInput);
	}

private:
	bool CacheOneStep() {
		if ((m_iCachePos < 0) || (m_iCachePos >= (int)m_oaFrames.size())) {
			m_IF.eprintf("CBQBReadCubeCache::CacheOneStep(): Error: m_iCachePos out of range (%d/%lu).\n",m_iCachePos,(unsigned long)m_oaFrames.size());
			abort();
		}
		if (m_oaFrames[m_iCachePos] != NULL)
			delete m_oaFrames[m_iCachePos];
		m_oaFrames[m_iCachePos] = new CBQBCubeFrame(m_IF);
		if (m_fInput == NULL) {
			m_IF.eprintf("CBQBReadCubeCache::CacheOneStep(): Error: m_fInput == NULL.\n");
			abort();
		}
		if (!m_oaFrames[m_iCachePos]->ReadFrame(m_fInput,m_iEps,m_iCSigni,m_iASigni,false))
			return false;
		if (m_fRef != NULL)
			m_oaFrames[m_iCachePos]->WriteFrame(m_fRef,false);
		return true;
	}

	int m_iEps;
	int m_iCSigni;
	int m_iASigni;
};

***************************************************************************************************************************/




template<class FT> class CBQBReadCache {
public:

	CBQBReadCache(CBQBInterface &i) : m_IF(i), m_iHistoryDepth(0), m_iCachePos(0), m_iReadPos(0),
		m_fRef(NULL), m_fInput(NULL), m_bStartLock(false), m_bInterfaceMode(false) {
	}

	virtual ~CBQBReadCache() {
		if (m_fInput != NULL)
			CloseFile();
		if (!m_bInterfaceMode)
			for (int z=0;z<(int)m_oaFrames.size();z++)
				if (m_oaFrames[z] != NULL)
					delete m_oaFrames[z];
		m_oaFrames.clear();
	}

	void SetStartLock() {
		if (m_bStartLock) {
			m_IF.eprintf("CBQBReadCache::SetStartLock(): Error: Start lock already set.\n");
			abort();
		}
		m_bStartLock = true;
	}

	void LiftStartLock() {
		if (!m_bStartLock) {
			m_IF.eprintf("CBQBReadCache::LiftStartLock(): Error: Start lock not in place.\n");
			abort();
		}
		m_bStartLock = false;
	}

	void SetHistoryDepth(int i) {
		m_oaFrames.resize(i,NULL);
	}

	int GetCachedSteps() const {
		if (m_iCachePos >= m_iReadPos)
			return m_iCachePos - m_iReadPos;
		else
			return m_iCachePos - m_iReadPos + (int)m_oaFrames.size();
	}

	const FT* GetFrameHistory(int i) const {
		if (i > m_iHistoryDepth) {
			m_IF.eprintf("CBQBReadCache::GetFrameHistory(): Error: Index out of range (%d/%d).\n",i,m_iHistoryDepth);
			abort();
		}
		int t = m_iReadPos - i;
		if (t < 0)
			t += (int)m_oaFrames.size();
		if (m_IF.IsPL( BQB_PL_VERBOSE ))
			m_IF.printf("GetFrameHistory(%d): readpos=%d, t=%d, size=%lu\n",i,m_iReadPos,t,(unsigned long)m_oaFrames.size());
		return m_oaFrames[t];
	}

	const FT* GetNextFrame() {
		if (m_iReadPos == m_iCachePos) {
			if (m_bInterfaceMode && (m_oaFrames.size() > 1)) {
				m_IF.eprintf("CBQBReadCache::GetNextFrame(): Error: ReadPos == CachePos.\n");
				abort();
			}
			m_iCachePos++;
			if (m_iCachePos >= (int)m_oaFrames.size()) {
				if (m_bStartLock) {
					m_IF.eprintf("CBQBReadCache::GetNextFrame(): Error: Violation of start lock.\n");
					abort();
				}
				m_iCachePos = 0;
			}
			if (!m_bInterfaceMode)
				if (!CacheOneStep())
					return NULL;
		}
		m_iReadPos++;
		if (m_iReadPos >= (int)m_oaFrames.size())
			m_iReadPos = 0;
		if (m_iHistoryDepth+1 < (int)m_oaFrames.size())
			m_iHistoryDepth++;
		return m_oaFrames[m_iReadPos];
	}

	void PushFrame( FT *frame ) {
		m_iCachePos++;
		if (m_iCachePos >= (int)m_oaFrames.size()) {
			if (m_bStartLock) {
				m_IF.eprintf("CBQBReadCache::GetNextFrame(): Error: Violation of start lock.\n");
				abort();
			}
			m_iCachePos = 0;
		}
		m_oaFrames[m_iCachePos] = frame;
		if (m_iHistoryDepth+1 < (int)m_oaFrames.size())
			m_iHistoryDepth++;
	}

	bool CacheSteps(int i, bool progress) {
		if (m_bInterfaceMode) {
			m_IF.eprintf("CBQBReadCache::CacheSteps(): Error: Can't call this function in interface mode.\n");
			m_IF.eprintf("\n");
			abort();
		}
		if (progress) {
			m_IF.printf("    Caching %d frames...\n",i);
			m_IF.printf("      [");
		}
		double tfs = i / 60.0;
		int k=0;
		bool b=false;
		for (int z=0;z<i;z++) {
			m_iCachePos++;
			if (m_iCachePos >= (int)m_oaFrames.size()) {
				if (m_bStartLock) {
					m_IF.eprintf("CBQBReadCache::CacheSteps(): Error: Violation of start lock.\n");
					abort();
				}
				m_iCachePos = 0;
				b = true;
			}
			if (!CacheOneStep()) {
				m_iCachePos--;
				if (b && (m_iCachePos < 0))
					m_iCachePos = (int)m_oaFrames.size()-1;
				m_IF.printf("] Only %d frames read.\n",k);
				return false;
			}
			k++;
			if (progress && (fmod(z,tfs) < 1.0))
				m_IF.printf("#");
		}
		if (progress)
			m_IF.printf("] Done.\n");
		return true;
	}

	void RewindReadPos() {
		if (!m_bStartLock) {
			m_IF.eprintf("CBQBReadCache::RewindReadPos(): Error: Requires start lock.\n");
			abort();
		}
		m_iReadPos = -1;
	}

	bool FileOpenRead(const char *inp, const char *ref = NULL) {
		if (m_bInterfaceMode) {
			m_IF.eprintf("CBQBReadCache::FileOpenRead(): Error: Can't call this function in interface mode.\n");
			m_IF.eprintf("\n");
			abort();
		}
		m_iCachePos = -1;
		m_iReadPos = -1;
		m_fInput = fopen(inp,"rb");
		if (m_fInput == NULL) {
			m_IF.eprintf("CBQBReadCache::FileOpenRead(): Error: Could not open input file for reading.\n");
			return false;
		}
		if (ref != NULL) {
			m_fRef = fopen(inp,"wb");
			if (m_fRef == NULL) {
				m_IF.eprintf("CBQBReadCache::FileOpenRead(): Error: Could not open reference file for writing.\n");
				fclose(m_fInput);
				m_fInput = NULL;
				return false;
			}
		}
		return true;
	}

	void Reset() {
		m_iCachePos = -1;
		m_iReadPos = -1;
	}

	void CloseFile() {
		if (m_bInterfaceMode) {
			m_IF.eprintf("CBQBReadCache::CloseFile(): Error: Can't call this function in interface mode.\n");
			m_IF.eprintf("\n");
			abort();
		}
		fclose(m_fInput);
		m_fInput = NULL;
		if (m_fRef != NULL) {
			fclose(m_fRef);
			m_fRef = NULL;
		}
	}

	long ftell() {
		if (m_bInterfaceMode) {
			m_IF.eprintf("CBQBReadCache::ftell(): Error: Can't call this function in interface mode.\n");
			m_IF.eprintf("\n");
			abort();
		}
		if (m_fInput == NULL)
			return 0;
		return ::ftell(m_fInput);
	}

	void ActivateInterfaceMode() {
		m_bInterfaceMode = true;
	}

	virtual bool SkipOneStep() = 0;


protected:
	CBQBInterface &m_IF;

	virtual bool CacheOneStep() = 0;

	int m_iHistoryDepth;
	int m_iCachePos; // Frame which has been cached last
	int m_iReadPos; // Frame which was read last
	std::vector<FT*> m_oaFrames;
	FILE *m_fRef;
	FILE *m_fInput;
	bool m_bStartLock;
	bool m_bInterfaceMode;
};



class CBQBReadXYZCache : public CBQBReadCache<CBQBAtomSet> {
	public:

	explicit CBQBReadXYZCache(CBQBInterface &i) : CBQBReadCache<CBQBAtomSet>(i), m_iSigni(-1) {
	}

	~CBQBReadXYZCache() {
	}

	void SetReadParameters(int signi) {
		m_iSigni = signi;
	}

	bool SkipOneStep() {
		return CBQBAtomSet(m_IF).SkipXYZ(m_fInput);
	}

private:
	bool CacheOneStep() {
		if ((m_iCachePos < 0) || (m_iCachePos >= (int)m_oaFrames.size())) {
			m_IF.eprintf("CBQBReadXYZCache::CacheOneStep(): Error: m_iCachePos out of range (%d/%lu).\n",m_iCachePos,(unsigned long)m_oaFrames.size());
			abort();
		}
		if (m_iSigni < 0) {
			m_IF.eprintf("CBQBReadXYZCache::CacheOneStep(): Error: m_iSigni out of range (%d).\n",m_iSigni);
			abort();
		}
		if (m_oaFrames[m_iCachePos] != NULL)
			delete m_oaFrames[m_iCachePos];
		m_oaFrames[m_iCachePos] = new CBQBAtomSet(m_IF);
		if (m_fInput == NULL) {
			m_IF.eprintf("CBQBReadXYZCache::CacheOneStep(): Error: m_fInput == NULL.\n");
			abort();
		}
		if (!m_oaFrames[m_iCachePos]->ReadXYZ(m_fInput,m_iSigni,m_fRef))
			return false;
		return true;
	}

	int m_iSigni;
};



class CBQBReadCubeCache : public CBQBReadCache<CBQBCubeFrame> {
public:

	explicit CBQBReadCubeCache(CBQBInterface &i)
		: CBQBReadCache<CBQBCubeFrame>(i), m_iEps(-1), m_iCSigni(-1), m_iASigni(-1) {
	}

	~CBQBReadCubeCache() {
	}

	void SetReadParameters(int eps, int csigni, int asigni) {
		m_iEps = eps;
		m_iCSigni = csigni;
		m_iASigni = asigni;
	}

	bool SkipOneStep() {
		if (m_bInterfaceMode) {
			m_IF.eprintf("CBQBReadCubeCache::SkipOneStep(): Error: Can't call this function in interface mode.\n");
			m_IF.eprintf("\n");
			abort();
		}
		return CBQBCubeFrame(m_IF).SkipFrame(m_fInput);
	}

private:
	bool CacheOneStep() {
		if (m_bInterfaceMode) {
			m_IF.eprintf("CBQBReadCubeCache::CacheOneStep(): Error: Can't call this function in interface mode.\n");
			m_IF.eprintf("\n");
			abort();
		}
		if ((m_iCachePos < 0) || (m_iCachePos >= (int)m_oaFrames.size())) {
			m_IF.eprintf("CBQBReadCubeCache::CacheOneStep(): Error: m_iCachePos out of range (%d/%lu).\n",m_iCachePos,(unsigned long)m_oaFrames.size());
			abort();
		}
		if (m_oaFrames[m_iCachePos] != NULL)
			delete m_oaFrames[m_iCachePos];
		m_oaFrames[m_iCachePos] = new CBQBCubeFrame(m_IF);
		if (m_fInput == NULL) {
			m_IF.eprintf("CBQBReadCubeCache::CacheOneStep(): Error: m_fInput == NULL.\n");
			abort();
		}
		if (!m_oaFrames[m_iCachePos]->ReadFrame(m_fInput,m_iEps,m_iCSigni,m_iASigni,false))
			return false;
		if (m_fRef != NULL)
			m_oaFrames[m_iCachePos]->WriteFrame(m_fRef,false);
		return true;
	}

	int m_iEps;
	int m_iCSigni;
	int m_iASigni;
};





class CBQBEngine {
public:

	explicit CBQBEngine(CBQBInterface &i);

	~CBQBEngine();

	void Reset();

	void DisableDeleteOutFrames() {
		m_bDeleteOutFrames = false;
	}


	void CompressString(
		const char *s,
		CBQBBitSet *bs,
		CBQBStatistics *stat
	);


	void DecompressString(
		CBQBBitSet *bs,
		char **s
	);


	void ExportCellInfo(CBQBOldCellInfo *ci, CBQBBitSet *bs);

	bool ImportCellInfo(CBQBOldCellInfo *ci, CBQBBitSet *bs);



	bool CubeToIntegerArray(
		std::vector<int> &outp,
		int order,
		CBQBParameterSet_Volumetric *parm,
		int &histused,
		bool skipcomp=false,
		double *resi=NULL
	);


	bool IntegerArrayToCube(
		std::vector<int> &inp,
		int order,
		CBQBParameterSet_Volumetric *parm,
		int ver,
		bool second
	);


	bool CompressCubeFrame(
		CBQBBitSet *bs,
		int corder,
		CBQBParameterSet_Volumetric *parm,
		int &histused,
		bool staticinfo,
		CBQBStatistics *stat,
		bool skipcomp=false,
		double *resi=NULL,
		std::vector<std::vector<int> > *tciaa=NULL
	);


	bool DecompressCubeFrame(
		CBQBBitSet *bs,
		int ver,
		CBQBParameterSet_Volumetric *parm,
		int *histused,
		bool second,
		bool check
	);


	bool CompressAtomFrame(
		CBQBBitSet *bs,
		int order,
		CBQBParameterSet_Position *parm,
		int &histused,
		bool staticinfo,
		bool comment,
		CBQBStatistics *stat,
		bool skipcomp=false,
		double *resi=NULL,
		std::vector<std::vector<int> > *tciaa=NULL
	);


	bool DecompressAtomStep(
		const std::vector<unsigned char> *inp,
		int ver,
		CBQBParameterSet_Position *parm,
		int *histused,
		bool second,
		bool check
	);


	bool DecompressAtomFrame(
		CBQBBitSet *bs,
		int ver,
		CBQBParameterSet_Position *parm,
		int *histused,
		bool second,
		bool check
	);


	bool PrepareDecompressCube();


	bool DecompressCubeStep(
		const std::vector<unsigned char> *inp,
		int ver,
		CBQBParameterSet_PosAndVol *parm,
		int *ahistused,
		int *chistused,
		bool skipvol=false
	);


	void ExportCubeHeader(
		CBQBBitSet *bs,
		int order
	);

	void ImportCubeHeader(
		CBQBBitSet *bs,
		bool second
	);

	void ExportAtomLabels(
		const CBQBAtomSet *as,
		CBQBBitSet *bs,
		CBQBStatistics *stat
	);

	void ImportAtomLabels(
		CBQBAtomSet *as,
		CBQBBitSet *bs
	);


	void AtomsToIntegerArray(
		std::vector<int> &outp,
		int order,
		CBQBParameterSet_Position *parm,
		int &histused,
		bool skipcomp=false,
		double *resi=NULL
	);


	void IntegerArrayToAtoms(
		const std::vector<int> &inp,
		int order,
		CBQBParameterSet_Position *parm,
		int ver,
		bool second
	);


	void ExportPOVRayScene(int degree, const char *s);


	void ExportXYZScene(int degree, const char *s);



	// This struct is a hack to enable the comparison function for std::sort
	//   to access class members of CBQBEngine.
	// Can't use global variables here, because multiple CBQBEngines
	//   might be running in different threads...

	struct SSortAtomOrd {

		explicit SSortAtomOrd(const CBQBEngine &e) : m_E(e) {
		}

		const CBQBEngine &m_E;

		bool operator()( const int & i1, const int & i2  ) {
			return m_E.m_iaAtomOrd[i1] > m_E.m_iaAtomOrd[i2];
		}
	};


	void ExportExtrapolatorSettings(
		CBQBBitSet *bs,
		bool volumetric
	);


	void ImportExtrapolatorSettings(
		CBQBBitSet *bs,
		bool volumetric
	);


	void PushOutputCubeFrame(CBQBCubeFrame *cfr);
	CBQBCubeFrame* GetOutputCubeFrame(int depth);

	void PushOutput2CubeFrame(CBQBCubeFrame *cfr);
	CBQBCubeFrame* GetOutput2CubeFrame(int depth);

	void PushOutputAtomFrame(CBQBAtomSet *cfr);
	CBQBAtomSet* GetOutputAtomFrame(int depth);

	void PushOutput2AtomFrame(CBQBAtomSet *cfr);
	CBQBAtomSet* GetOutput2AtomFrame(int depth);

/*	void PushInputAtomFrame_NoDelete(CBQBAtomSet *cfr);
	void PushInputCubeFrame_NoDelete(CBQBCubeFrame *cfr);*/

	void BuildHilbertIdx(
		int resx,
		int resy,
		int resz
	);

	void BuildAtomSort(const CBQBAtomSet *as);

	void BuildIdentityAtomSort(const CBQBAtomSet *as);

	int m_iOutputCubeBufPos;
	std::vector<CBQBCubeFrame*> m_oaOutputCubeBuf;

	int m_iOutput2CubeBufPos;
	std::vector<CBQBCubeFrame*> m_oaOutput2CubeBuf;

	int m_iOutputAtomBufPos;
	std::vector<CBQBAtomSet*> m_oaOutputAtomBuf;

	int m_iOutput2AtomBufPos;
	std::vector<CBQBAtomSet*> m_oaOutput2AtomBuf;

	std::vector<int> m_iaHilbertIdx;

	std::vector<double> m_faTempPred;
	std::vector<double> m_faTempPred2;

	std::vector<int> m_iaAtomSort;

	std::vector<int> m_iaTempCompressCube;
	std::vector<int> m_iaTempCompressXYZ;

	CBQBExtrapolator *m_pExtrapolator;
	CBQBExtrapolator *m_pExtrapolatorCorr;
	CBQBExtrapolator *m_pExtrapolatorXYZ;

	CBQBReadXYZCache *m_pReadCacheXYZ;
	CBQBReadCubeCache *m_pReadCacheCube;

	std::vector<int> m_iaAtomOrd;

private:

	bool m_bDeleteOutFrames;

	CBQBInterface &m_IF;

};


#endif


