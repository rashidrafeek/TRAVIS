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


#ifndef RANDOM_H
#define RANDOM_H


// This must always be the first include directive
#include "config.h"

#include "xobject.h"
#include "ziggurat.h"
#include <stdlib.h>


class CRandom : public CxObject
{
public:
	CRandom();
	~CRandom();

	double RandomUniform()
	{
//		return r4_uni(&seed);
		return ((double)(rand()%RAND_MAX)/RAND_MAX/RAND_MAX) + ((double)(rand()%RAND_MAX)/RAND_MAX);
	}

	double RandomNormal()
	{
		return r4_nor(&seed,kn,fn,wn);
	}

	double RandomExp()
	{
		return r4_exp(&seed,ke,fe,we);
	}

	double fn[128];
	int kn[128];
	double wn[128];
	double fe[256];
	int ke[256];
	double we[256];
	unsigned long int seed;
};

#endif
