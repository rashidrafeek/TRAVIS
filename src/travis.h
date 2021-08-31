/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

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


#ifndef TRAVIS_H
#define TRAVIS_H

// This must always be the first include directive
#include "config.h"

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include "bintools.h"
#include "moltools.h"
#include "xobject.h"
#include "xobarray.h"
#include "asciiart.h"
#include "backtrace.h"
#include "grace.h"
#include "acf.h"
#include "fft.h"
#include "spectrum.h"
#include "analysis.h"
#include "analysisgroup.h"
#include "timestep.h"
#include "database.h"
#include "xintarray.h"
#include "globalvar.h"


void mprintf(const char *s, ...);
void strtolower(char *s);
/*double AtomMass(char *s);
double AtomRadius(char *s);
int AtomOrd(char *s);*/
//double AtomVDWRadius(char *s);
bool GatherInfos();
bool ParseAtom(const char *s, int refmol, int &ty, int &rty, int &atom);
void DumpAnalyses();
bool ParseFunctions(const char *s);
bool ParseRefSystem(int refmol, const char *s, int points);
bool ParsePeriodic(const char *s);
double GuessBoxSize();
void AddElement(const char *s, int ord, double mass, double radius, double vdw);
void AddElement(const char *s, int ord, double mass, double radius, double vdw, double ncs);
void AddElementData();
void SetElementColor(const char *s, unsigned char r, unsigned char g, unsigned char b, unsigned char br, unsigned char bb, unsigned char bg);




#endif
