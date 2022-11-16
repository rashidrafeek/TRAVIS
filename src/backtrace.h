/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Martin Brehm.

    ---------------------------------------------------------------------------

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*****************************************************************************/


#ifndef BACKTRACE_H
#define BACKTRACE_H


// This must always be the first include directive
#include "config.h"

#include <stdlib.h>
#include <signal.h>
#include "tools.h"
#include "xobject.h"
#include "xobarray.h"


class CxTracePoint : public CxObject
{
public:
	CxTracePoint() { }
	~CxTracePoint() { }
	char m_sFunction[256];
	char m_sFile[64];
	unsigned long m_iLine;
	char m_sMessage[256];
};

class CxObArray;

extern CxObArray g_oaBackTrace;
extern char *g_sExeName;

#ifdef DEBUG_BACKTRACE

#define BTDUMP DumpBacktrace()

#ifdef TARGET_WINDOWS

#define BTIN { if (g_oaBackTrace.GetSize() > 500) { eprintf("*** Trace to deep ***\n"); DumpBacktrace(); abort(); } CxTracePoint *ctp = new CxTracePoint(); strcpy(ctp->m_sFile,__FILE__); ctp->m_sFunction[0] = 0; ctp->m_iLine = __LINE__; ctp->m_sMessage[0] = 0; g_oaBackTrace.Add(ctp); }
#define BTOUT { if (g_oaBackTrace.GetSize() == 0) { eprintf("*** BTOUT Error ***\n"); abort(); } delete (CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1]; g_oaBackTrace.RemoveAt(g_oaBackTrace.GetSize()-1,1); }
#define BTL ((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_iLine = __LINE__;
#define BTMSG(msg) { ((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_iLine = __LINE__; strcpy(((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_sMessage,msg); }

#else // not TARGET_WINDOWS

#define BTIN { if (g_oaBackTrace.GetSize() > 500) { eprintf("*** Trace to deep ***\n"); DumpBacktrace(); abort(); } CxTracePoint *ctp = new CxTracePoint(); strcpy(ctp->m_sFile,__FILE__); strcpy(ctp->m_sFunction,__PRETTY_FUNCTION__); ctp->m_iLine = __LINE__; ctp->m_sMessage[0] = 0; g_oaBackTrace.Add(ctp); }
#define BTOUT { if (g_oaBackTrace.GetSize() == 0) { eprintf("*** BTOUT Error ***\n"); abort(); } delete (CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1]; g_oaBackTrace.RemoveAt(g_oaBackTrace.GetSize()-1,1); }
#define BTL ((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_iLine = __LINE__;
#define BTMSG(msg) { ((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_iLine = __LINE__; strcpy(((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_sMessage,msg); }

#endif // end TARGET_WINDOWS

#else // not DEBUG_BACKTRACE

#define BTDUMP
#define BTIN
#define BTOUT
#define BTL
#define BTMSG(msg)

#endif // end DEBUG_BACKTRACE

#ifdef DEBUG_EXTENDED_BACKTRACE

#ifdef TARGET_WINDOWS

#define BXIN { if (g_oaBackTrace.GetSize() > 500) { eprintf("*** Trace to deep ***\n"); DumpBacktrace(); abort(); } CxTracePoint *ctp = new CxTracePoint(); strcpy(ctp->m_sFile,__FILE__); ctp->m_sFunction[0] = 0; ctp->m_iLine = __LINE__; ctp->m_sMessage[0] = 0; g_oaBackTrace.Add(ctp); }
#define BXOUT { if (g_oaBackTrace.GetSize() == 0) { eprintf("*** BXOUT Error ***\n"); abort(); } delete (CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1]; g_oaBackTrace.RemoveAt(g_oaBackTrace.GetSize()-1,1); }
#define BXL ((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_iLine = __LINE__;
#define BXMSG(msg) { ((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_iLine = __LINE__; strcpy(((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_sMessage,msg); }

#else // not TARGET_WINDOWS

#define BXIN { if (g_oaBackTrace.GetSize() > 500) { eprintf("*** Trace to deep ***\n"); DumpBacktrace(); abort(); } CxTracePoint *ctp = new CxTracePoint(); strcpy(ctp->m_sFile,__FILE__); strcpy(ctp->m_sFunction,__PRETTY_FUNCTION__); ctp->m_iLine = __LINE__; ctp->m_sMessage[0] = 0; g_oaBackTrace.Add(ctp); }
#define BXOUT { if (g_oaBackTrace.GetSize() == 0) { eprintf("*** BXOUT Error ***\n"); abort(); } delete (CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1]; g_oaBackTrace.RemoveAt(g_oaBackTrace.GetSize()-1,1); }
#define BXL ((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_iLine = __LINE__;
#define BXMSG(msg) { ((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_iLine = __LINE__; strcpy(((CxTracePoint*)g_oaBackTrace[g_oaBackTrace.GetSize()-1])->m_sMessage,msg); }

#endif // end TARGET_WINDOWS

#else // not DEBUG_EXTENDED_BACKTRACE

#define BXIN
#define BXOUT
#define BXL
#define BXMSG(msg)

#endif // end DEBUG_EXTENDED_BACKTRACE


void SIGNAL_SEGV(int param);
void InstallSignalHandler();
void UninstallSignalHandler();
void DumpBacktrace();


void NewException(double d, const char *filename, int line, const char *function)
#if __cplusplus > 199711L
__attribute__((noreturn))
#endif
;

void NewException(double d, const char *filename, int line, const char *function, const char *info)
#if __cplusplus > 199711L
__attribute__((noreturn))
#endif
;

void BoundsException(int i, int j, const char *filename, int line, const char *function)
#if __cplusplus > 199711L
__attribute__((noreturn))
#endif
;

void BoundsException(int i, int j, const char *filename, int line, const char *function, const char *info)
#if __cplusplus > 199711L
__attribute__((noreturn))
#endif
;

#endif

