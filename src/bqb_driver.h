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



#ifndef BQB_DRIVER_H
#define BQB_DRIVER_H


// This must always be the first include directive
#include "bqb_config.h"

#include "bqb_engine.h"
#include "bqb_interface.h"
#include <vector>



class CBQBDriver {
public:

	explicit CBQBDriver(CBQBInterface &i);

	~CBQBDriver();

	bool MergeBQB(const char *outfile, const std::vector<const char*> &infile);

//	bool SplitBQB(const char *infile, const char *outbase, int splitlength, int steps);

	bool CompressFile(
		const char *inp,
		const char *outp,
		const CBQBParameterSet_Compressor *parm
	);
	
	bool DecompressFile(const char *inp, const char *outp, const char *ref);

	bool CompressCube(
		const char *inp,
		const char *outp,
		const char *ref,
		int start,
		int steps,
		int stride,
		int keyfreq,
		CBQBParameterSet_PosAndVol *parm,
		int optimize,
		int optsteps,
		bool onlyopt,
		bool comment,
		bool compare,
		bool dummyread
	);


	double CompressCubeBenchmark(
		CBQBParameterSet_Volumetric *parm,
		bool realsize,
		double *presi,
		std::vector<std::vector<int> > *tciaa
	);


	double OptimizeCube_ObjectiveTExpo(double texpo);

	double OptimizeCube_ObjectiveSExpo(double sexpo);


	bool OptimizeCubeParameters(
		int     level,
		CBQBParameterSet_Volumetric *parm
	);

	
	bool DecompressCube(
		const char *inp,
		const char *outp,
		const char *readref,
		int steps,
		int stride,
		CBQBParameterSet_PosAndVol *parm
	);


	bool CompressXYZ(
		const char *inp, 
		const char *outp, 
		const char *ref,
		int start,
		int steps, 
		int stride, 
		int keyfreq,
		CBQBParameterSet_Position *parm,
		bool comment, 
		bool compare,
		int optimize,
		int optsteps,
		bool onlyopt
	);


	double CompressXYZBenchmark(
		CBQBParameterSet_Position *parm,
		bool comment, 
		bool realsize,
		double *presi,
		std::vector<std::vector<int> > *tciaa
	);


	double GoldenSectionSearch(double maxl, double eps, double (CBQBDriver::*fobj)(double));


	double REC_GoldenSectionSearch(
		double a,
		double b,
		double c,
		double fa,
		double fb,
		double fc,
		double (CBQBDriver::*fobj)(double),
		double eps,
		int depth
	);

		
	double OptimizeXYZ_ObjectiveTExpo(double texpo);


	bool OptimizeXYZParameters(
		int     level,
		CBQBParameterSet_Position *parm,
		bool    comment, 
		bool    cubepart
	);


	bool DecompressXYZ(
		const char *inp,
		const char *outp,
		const char *readref,
		int steps,
		int stride,
		CBQBParameterSet_Position *parm
	);



	CBQBEngine *m_pEngine;


	CBQBParameterSet_Volumetric *m_pVolParm;
	CBQBParameterSet_Position *m_pPosParm;


	CBQBParameterSet_PosAndVol m_ParmPosAndVol;
	CBQBParameterSet_Position m_ParmPos;


	std::vector<int> m_iaAtomOrd;
	int m_iHistogramCounter;


	bool m_bOptCubeRealSize;

	bool m_bOptXYZComment;
	bool m_bOptXYZRealSize;


	int m_iOptIncludeFirst;

private:
	CBQBInterface &m_IF;
};




#endif


