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



#ifndef BQB_LINALG_H
#define BQB_LINALG_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_math.h"



class CBQBInterface;



class CBQBLinAlg {
public:

	explicit CBQBLinAlg(CBQBInterface &i) : m_IF(i) {
	}

	int BQB_ComputeSVD_Flat(double *a, int m, int n, double *w, double *v);

	void BQB_TestSVD();

	void BQB_ComputePseudoInverse(int m, double *m_in, double *m_out);

	CBQBDMatrixMN BQB_ComputePseudoInverse(CBQBDMatrixMN input);

	void BQB_TestPseudoInverse();

private:
	CBQBInterface &m_IF;
};



#endif




