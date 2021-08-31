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



// This must always be the first include directive
#include "bqb_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include "bqb.h"


const char *GetRevisionInfo_bqbtool(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_bqbtool() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



// Headers required for Linux stack trace
#ifndef BQB_INSIDE_TRAVIS
	#if !defined(_WIN32) && (defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__)))
		#include <execinfo.h>
		#include <signal.h>
	#endif
#endif


// Headers required for obtaining host name and working directory
#if !defined(_WIN32) && (defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__)))
	#include <unistd.h>
#elif defined(_WIN32)
	#include <windows.h>
	#include <direct.h>
#endif



char *gc_sExeName = NULL;
char *gc_sWorkingDir = NULL;
char *gc_sHostName = NULL;
FILE *gc_fLogFile = NULL;

CBQBInterface *gc_pBQBInterface = NULL;

bool gc_bCompress;
bool gc_bDecompress;
bool gc_bCheck;
bool gc_bCompare;
bool gc_bMerge;
bool gc_bSplit;

bool gc_bFile;
bool gc_bCube;
bool gc_bXYZ;

const char *gc_sInFile;
const char *gc_sOutFile;
const char *gc_sRefFile;
std::vector<const char*> gc_sInFileList;


CBQBParameterSet_Compressor *gc_pParmFile = NULL;
int gc_iFileVerbose;

CBQBParameterSet_PosAndVol *gc_pParmVol = NULL;

int gc_iCubeStart;
int gc_iCubeSteps;
int gc_iCubeStride;
bool gc_bCubeKeepComment;
bool gc_bCubeCompare;
const char *gc_sCubeRefOut;
int gc_iCubeVerbose;
bool gc_bCubeDummyRead;
int gc_iCubeKey;


int gc_iCubeOptimize;
int gc_iCubeOptSteps;
bool gc_bCubeOnlyOpt;
int gc_iCubeOptIncludeFirst;


CBQBParameterSet_Position *gc_pParmPos = NULL;

int gc_iAtomStart;
int gc_iAtomSteps;
int gc_iAtomStride;
bool gc_bAtomKeepComment;
bool gc_bAtomCompare;
const char *gc_sAtomRefOut;
int gc_iAtomVerbose;
int gc_iAtomKey;


int gc_iAtomOptimize;
int gc_iAtomOptSteps;
bool gc_bAtomOnlyOpt;
int gc_iAtomOptIncludeFirst;

int gc_iCheckVerbose;
int gc_iCompareVerbose;
int gc_iMergeVerbose;
int gc_iSplitVerbose;

int gc_iSplitLength;
int gc_iSplitSteps;

bool gc_bOptTouched;
bool gc_bDryRun;





void BQBWriteRevisionInfo() {

	std::vector<const char*> sa;
	unsigned int z, len;

	gc_pBQBInterface->printf("Revision information of the source code files:\n\n");

	len = 0;

_begin:
	sa.push_back( GetRevisionInfo_bqb_alphabet( len ) );
	sa.push_back( GetRevisionInfo_bqb_bitset( len ) );
	sa.push_back( GetRevisionInfo_bqb_crc( len ) );
	sa.push_back( GetRevisionInfo_bqb_cubeframe( len ) );
	sa.push_back( GetRevisionInfo_bqb_driver( len ) );
	sa.push_back( GetRevisionInfo_bqb_engine( len ) );
	sa.push_back( GetRevisionInfo_bqb_extrapolator( len ) );
	sa.push_back( GetRevisionInfo_bqb_format( len ) );
	sa.push_back( GetRevisionInfo_bqb_hilbert( len ) );
	sa.push_back( GetRevisionInfo_bqb_hufftree( len ) );
	sa.push_back( GetRevisionInfo_bqb_integerengine( len ) );
	sa.push_back( GetRevisionInfo_bqb_interface( len ) );
	sa.push_back( GetRevisionInfo_bqb_largeinteger( len ) );
	sa.push_back( GetRevisionInfo_bqb_linalg( len ) );
	sa.push_back( GetRevisionInfo_bqb_math( len ) );
	sa.push_back( GetRevisionInfo_bqb_parmset( len ) );
	sa.push_back( GetRevisionInfo_bqb_tools( len ) );
	sa.push_back( GetRevisionInfo_bqbtool( len ) );

	if (len == 0) {
		for (z=0;z<sa.size();z++)
			if (strlen(sa[z]) > len)
				len = (unsigned int)strlen(sa[z]);
		sa.clear();
		goto _begin;
	}

	for (z=0;z<sa.size();z++)
		gc_pBQBInterface->printf( "  %s\n", sa[z] );

	gc_pBQBInterface->printf( "\n");
}


void BQBCheckSourceVersion() {

	const char *p;
	bool b;


	b = false;

	p = GetSourceVersion_bqb_alphabet();

	if (strcmp( p, GetSourceVersion_bqb_bitset() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_crc() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_cubeframe() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_driver() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_engine() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_extrapolator() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_format() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_hilbert() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_hufftree() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_integerengine() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_interface() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_largeinteger() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_linalg() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_math() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_parmset() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqb_tools() ) != 0) b = true;
	if (strcmp( p, GetSourceVersion_bqbtool() ) != 0) b = true;

	if (b) {

		gc_pBQBInterface->printf("\n");
		gc_pBQBInterface->eprintf("!!!!    Warning    !!!!\n");
		gc_pBQBInterface->printf("The object files do not belong to the same source code version.\n");
		gc_pBQBInterface->printf("Something went wrong with compiling. Expect problems and crashes.\n");
		gc_pBQBInterface->printf("\n");

		BQBWriteRevisionInfo();
	}
}


void CC_InitGlobalVars() {

	gc_bOptTouched = false;

	gc_sInFile = NULL;
	gc_sOutFile = NULL;
	gc_sRefFile = NULL;

	gc_bCompress = false;
	gc_bDecompress = false;
	gc_bCheck = false;
	gc_bCompare = false;
	gc_bMerge = false;
	gc_bSplit = false;
	gc_bFile = false;
	gc_bCube = false;
	gc_bXYZ = false;

	gc_pParmFile = new CBQBParameterSet_Compressor(*gc_pBQBInterface);

	gc_pParmFile->SetTableCount( 6 );
	gc_pParmFile->SetOptTables( false );
	gc_pParmFile->SetBlockLength( 50 );
	gc_pParmFile->SetBW( true );
	gc_pParmFile->SetMTF( true );
	gc_pParmFile->SetRLE( true );
	gc_pParmFile->SetMaxIter( 10 );
	gc_pParmFile->SetMaxChunk( 4194304 );

	gc_iFileVerbose = 0;

	gc_iCubeStart = 0;
	gc_iCubeSteps = -1;
	gc_iCubeStride = 1;
	gc_bCubeKeepComment = true;
	gc_bCubeCompare = true;
	gc_sCubeRefOut = NULL;
	gc_bCubeDummyRead = false;
	gc_iCubeVerbose = 0;
	gc_iCubeKey = 0;

	gc_pParmVol = new CBQBParameterSet_PosAndVol(*gc_pBQBInterface);

	gc_pParmVol->SetVolSigni(5);
	gc_pParmVol->SetVolEps(12);
	gc_pParmVol->SetVolOrder(8);
	gc_pParmVol->SetVolOptOrder(true);
	gc_pParmVol->SetVolHilbert(true);
	gc_pParmVol->SetVolNbhFac(1.075);
	gc_pParmVol->SetVolSplit(10);
	gc_pParmVol->SetVolTableCount(6);
	gc_pParmVol->SetVolOptTables(false);
	gc_pParmVol->SetVolBlockLength(20);
	gc_pParmVol->SetVolBW(false);
	gc_pParmVol->SetVolMTF(false);
	gc_pParmVol->SetVolMaxIter(10);
	gc_pParmVol->SetVolRLE(true);
	gc_pParmVol->SetVolMaxChunk(0);
	gc_pParmVol->SetPosPrecision(6);
	gc_pParmVol->SetPosOrder(8);
	gc_pParmVol->SetPosOptOrder(true);
	gc_pParmVol->SetPosSplit(14);
	gc_pParmVol->SetPosTableCount(1);
	gc_pParmVol->SetPosOptTables(false);
	gc_pParmVol->SetPosBlockLength(40);
	gc_pParmVol->SetPosBW(false);
	gc_pParmVol->SetPosMTF(false);
	gc_pParmVol->SetPosRLE(true);
	gc_pParmVol->SetPosMaxIter(10);
	gc_pParmVol->SetPosSortAtom(true);
	gc_pParmVol->SetPosMaxChunk(0);

	gc_pParmVol->SetVolUseExtra(true);
	gc_pParmVol->SetVolExtraSRange(7);
	gc_pParmVol->SetVolExtraTRange(6);
	gc_pParmVol->SetVolExtraSOrder(3);
	gc_pParmVol->SetVolExtraTOrder(2);
	gc_pParmVol->SetVolExtraOffset(3);
	gc_pParmVol->SetVolExtraCrossS(true);
	gc_pParmVol->SetVolExtraCrossT(false);
	gc_pParmVol->SetVolExtraWrap(true);
	gc_pParmVol->SetVolExtraCrossRangeS(true);
	gc_pParmVol->SetVolExtraCrossRangeT(false);
	gc_pParmVol->SetVolExtraDistExpo(3.0);
	gc_pParmVol->SetVolExtraTimeExpo(1.0);
	gc_pParmVol->SetVolExtraPredCorr(true);

	gc_pParmVol->SetPosUseExtra(true);
	gc_pParmVol->SetPosExtraTRange(9);
	gc_pParmVol->SetPosExtraTOrder(6);
	gc_pParmVol->SetPosExtraTimeExpo(4.0);

	gc_iCubeOptimize = 1;
	gc_iCubeOptSteps = -1;
	gc_bCubeOnlyOpt = false;
	gc_iCubeOptIncludeFirst = -1;

	gc_iAtomStart = 0;
	gc_iAtomSteps = -1;
	gc_iAtomStride = 1;
	gc_bAtomKeepComment = true;
	gc_bAtomCompare = true;
	gc_sAtomRefOut = NULL;
	gc_iAtomVerbose = 0;
	gc_iAtomKey = 0;

	gc_bDryRun = false;

	gc_pParmPos = new CBQBParameterSet_Position(*gc_pBQBInterface);

	gc_pParmPos->SetPrecision(5);
	gc_pParmPos->SetOrder(8);
	gc_pParmPos->SetOptOrder(true);
	gc_pParmPos->SetSplit(14);
	gc_pParmPos->SetTableCount(1);
	gc_pParmPos->SetOptTables(false);
	gc_pParmPos->SetBlockLength(40);
	gc_pParmPos->SetBW(false);
	gc_pParmPos->SetMTF(false);
	gc_pParmPos->SetRLE(true);
	gc_pParmPos->SetMaxIter(10);
	gc_pParmPos->SetSortAtom(true);
	gc_pParmPos->SetMaxChunk(0);
	gc_pParmPos->SetUseExtra(true);
	gc_pParmPos->SetExtraTRange(9);
	gc_pParmPos->SetExtraTOrder(6);
	gc_pParmPos->SetExtraTimeExpo(4.0);

	gc_iAtomOptimize = 2;
	gc_iAtomOptSteps = -1;
	gc_bAtomOnlyOpt = false;
	gc_iAtomOptIncludeFirst = -1;

	gc_iCheckVerbose = 0;
	gc_iCompareVerbose = 0;
	gc_iMergeVerbose = 0;
	gc_iSplitVerbose = 0;

	gc_iSplitLength = -1;
	gc_iSplitSteps = -1;
}


void CC_CommandLineHelp() {

	gc_pBQBInterface->printf("\n");
	gc_pBQBInterface->bprintf("Command line usage:\n\n");
	gc_pBQBInterface->printf("    compress   {file|voltraj|postraj} [options] <infile> <outfile.bqb>\n");
	gc_pBQBInterface->printf("    decompress {file|voltraj|postraj} [options] <infile.bqb> <outfile>\n");
	gc_pBQBInterface->printf("    check      [options] <infile.bqb>\n");
/*	gc_pBQBInterface->printf("    compare    [options] <infile.bqb> <xyz file>\n");
	gc_pBQBInterface->printf("    merge      [options] <outfile.bqb> <infile1.bqb> [infile2.bqb ...]\n");
	gc_pBQBInterface->printf("    split      [options] <infile.bqb> <outbasename> <steps/piece>\n\n");*/
	gc_pBQBInterface->printf("\n");

	gc_pBQBInterface->printf("More operating modes (check, compare, split, merge) will be added soon.\n");
	gc_pBQBInterface->printf("\n");

	gc_pBQBInterface->printf("For volumetric trajectories, currently only the Gaussian Cube file format\n");
	gc_pBQBInterface->printf("is supported.\n");
	gc_pBQBInterface->printf("\n");

	gc_pBQBInterface->printf("For position trajectories, currently only the XMol XYZ file format\n");
	gc_pBQBInterface->printf("is supported. More formats will be added soon.\n");
	gc_pBQBInterface->printf("\n");

	gc_pBQBInterface->printf("In case of \"compress\" or \"decompress\", <outfile> may be left out,\n");
	gc_pBQBInterface->printf("then no output is written (e.g. to check compression ratio or timings).\n\n");

	gc_pBQBInterface->printf("To see a list of possible options and their default values,\n");
	gc_pBQBInterface->printf("enter the corresponding command without input and output file names.\n");

	gc_pBQBInterface->printf("\n");
}


bool CC_ParseArguments(int argc, const char *argv[]) {

	int z, i, ti;
	double tf;
	bool b;

	gc_pBQBInterface->printf("Parsing command line arguments...\n");

	if (argc < 3) {
		gc_pBQBInterface->eprintf("Error: Not enough command line arguments.\n");
		return false;
	}

	if (strcmp(argv[1],"compress") == 0)
		gc_bCompress = true;
	else if (strcmp(argv[1],"decompress") == 0)
		gc_bDecompress = true;
	else if (strcmp(argv[1],"check") == 0)
		gc_bCheck = true;
/*	else if (strcmp(argv[1],"compare") == 0)
		gc_bCompare = true;
	else if (strcmp(argv[1],"merge") == 0)
		gc_bMerge = true;
	else if (strcmp(argv[1],"split") == 0)
		gc_bSplit = true;*/
	else {
		gc_pBQBInterface->eprintf("Error: Unknown first argument \"%s\".\n",argv[1]);
		return false;
	}

	if (gc_bCheck) {

		b = true;
		i = 0;
		for (z=2;z<argc;z++) {

			if ((argv[z][0] == '-') && b) {

				if (strcmp(argv[z],"-verbose") == 0) {
					if (z+1 == argc) {
						gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
						return false;
					}
					if (strcmp(argv[z+1],"2") == 0) {
						gc_iCheckVerbose = 2;
						z++;
					} else if (strcmp(argv[z+1],"1") == 0) {
						gc_iCheckVerbose = 1;
						z++;
					} else if (strcmp(argv[z+1],"0") == 0) {
						gc_iCheckVerbose = 0;
						z++;
					} else
						gc_iCheckVerbose = 1;
					gc_pBQBInterface->printf("  -verbose %d\n",gc_iCheckVerbose);
					continue;
				}

				gc_pBQBInterface->eprintf("Error: Unknown command line argument \"%s\".\n",argv[z]);
				return false;

			} else {

				b = false;
				if (i == 0) {
					gc_sInFile = argv[z];
					gc_pBQBInterface->printf("  infile %s\n",gc_sInFile);
				} else {
					gc_pBQBInterface->eprintf("Error: Unrecognized trailing argument: \"%s\".\n",argv[z]);
					return false;
				}
				i++;
			}
		}

	} else if (gc_bSplit) {

		b = true;
		i = 0;
		for (z=2;z<argc;z++) {

			if ((argv[z][0] == '-') && b) {

				if (strcmp(argv[z],"-steps") == 0) {
					if (z+1 == argc) {
						gc_pBQBInterface->eprintf("Error: Missing value after \"-steps\".\n");
						return false;
					}
					z++;
					gc_iSplitSteps = atoi(argv[z]);
					if ((gc_iSplitSteps < -1) || (gc_iSplitSteps == 0)) {
						gc_pBQBInterface->eprintf("Error: Invalid value for \"-steps\": \"%s\".\n",argv[z]);
						gc_pBQBInterface->printf("       Allowed values are -1 (=all) or > 0.\n");
						return false;
					}
					gc_pBQBInterface->printf("  -steps %d\n",gc_iSplitSteps);
					continue;
				}

				if (strcmp(argv[z],"-verbose") == 0) {
					if (z+1 == argc) {
						gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
						return false;
					}
					if (strcmp(argv[z+1],"2") == 0) {
						gc_iSplitVerbose = 2;
						z++;
					} else if (strcmp(argv[z+1],"1") == 0) {
						gc_iSplitVerbose = 1;
						z++;
					} else if (strcmp(argv[z+1],"0") == 0) {
						gc_iSplitVerbose = 0;
						z++;
					} else
						gc_iSplitVerbose = 1;
					gc_pBQBInterface->printf("  -verbose %d\n",gc_iSplitVerbose);
					continue;
				}

				gc_pBQBInterface->eprintf("Error: Unknown command line argument \"%s\".\n",argv[z]);
				return false;

			} else {

				b = false;
				if (i == 0) {
					gc_sInFile = argv[z];
					gc_pBQBInterface->printf("  infile %s\n",gc_sInFile);
				} else if (i == 1) {
					gc_sOutFile = argv[z];
					gc_pBQBInterface->printf("  outbasename %s\n",gc_sOutFile);
				} else if (i == 2) {
					gc_iSplitLength = atoi(argv[z]);
					gc_pBQBInterface->printf("  splitlength %d\n",gc_iSplitLength);
				} else {
					gc_pBQBInterface->eprintf("Error: Unrecognized trailing argument: \"%s\".\n",argv[z]);
					return false;
				}
				i++;
			}
		}

	} else if (gc_bCompare) {

		b = true;
		i = 0;
		for (z=2;z<argc;z++) {

			if ((argv[z][0] == '-') && b) {

				if (strcmp(argv[z],"-verbose") == 0) {
					if (z+1 == argc) {
						gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
						return false;
					}
					if (strcmp(argv[z+1],"2") == 0) {
						gc_iCompareVerbose = 2;
						z++;
					} else if (strcmp(argv[z+1],"1") == 0) {
						gc_iCompareVerbose = 1;
						z++;
					} else if (strcmp(argv[z+1],"0") == 0) {
						gc_iCompareVerbose = 0;
						z++;
					} else
						gc_iCompareVerbose = 1;
					gc_pBQBInterface->printf("  -verbose %d\n",gc_iCompareVerbose);
					continue;
				}

				gc_pBQBInterface->eprintf("Error: Unknown command line argument \"%s\".\n",argv[z]);
				return false;

			} else {

				b = false;
				if (i == 0) {
					gc_sInFile = argv[z];
					gc_pBQBInterface->printf("  infile %s\n",gc_sInFile);
				} else if (i == 1) {
					gc_sRefFile = argv[z];
					gc_pBQBInterface->printf("  reffile %s\n",gc_sRefFile);
				} else {
					gc_pBQBInterface->eprintf("Error: Unrecognized trailing argument: \"%s\".\n",argv[z]);
					return false;
				}
				i++;
			}
		}

	} else if (gc_bMerge) {

		b = true;
		i = 0;
		for (z=2;z<argc;z++) {

			if ((argv[z][0] == '-') && b) {

				if (strcmp(argv[z],"-verbose") == 0) {
					if (z+1 == argc) {
						gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
						return false;
					}
					if (strcmp(argv[z+1],"2") == 0) {
						gc_iMergeVerbose = 2;
						z++;
					} else if (strcmp(argv[z+1],"1") == 0) {
						gc_iMergeVerbose = 1;
						z++;
					} else if (strcmp(argv[z+1],"0") == 0) {
						gc_iMergeVerbose = 0;
						z++;
					} else
						gc_iMergeVerbose = 1;
					gc_pBQBInterface->printf("  -verbose %d\n",gc_iMergeVerbose);
					continue;
				}

				gc_pBQBInterface->eprintf("Error: Unknown command line argument \"%s\".\n",argv[z]);
				return false;

			} else {

				b = false;
				if (i == 0) {
					gc_sOutFile = argv[z];
					gc_pBQBInterface->printf("  outfile %s\n",gc_sInFile);
				} else {
					gc_sInFileList.push_back(argv[z]);
					gc_pBQBInterface->printf("  infile %s\n",argv[z]);
				}
				i++;
			}
		}

	} else {

		if (strcmp(argv[2],"file") == 0)
			gc_bFile = true;
		else if (strcmp(argv[2],"voltraj") == 0)
			gc_bCube = true;
		else if (strcmp(argv[2],"postraj") == 0)
			gc_bXYZ = true;
		else {
			gc_pBQBInterface->eprintf("Error: Unknown second argument \"%s\".\n",argv[2]);
			return false;
		}

		b = true;
		i = 0;
		for (z=3;z<argc;z++) {

			if (gc_bCompress) {

				if (gc_bFile) {

					if ((argv[z][0] == '-') && b) {

						if (strcmp(argv[z],"-dryrun") == 0) {
							if (z+1 != argc) {
								if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
									gc_bDryRun = true;
									z++;
								} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
									gc_bDryRun = false;
									z++;
								} else
									gc_bDryRun = true;
							} else
								gc_bDryRun = true;
							gc_pBQBInterface->printf("  -dryrun %s\n",gc_bDryRun?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-pkey") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pkey\".\n");
								return false;
							}
							z++;
							if (!gc_pParmFile->FromKey(argv[z])) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-pkey\": \"%s\".\n",argv[z]);
								return false;
							}
							gc_pBQBInterface->printf("  -pkey %s\n",argv[z]);
							continue;
						}

						if (strcmp(argv[z],"-tables") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-tables\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 31)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-tables\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 31.\n");
								return false;
							}
							gc_pParmFile->SetTableCount(ti);
							gc_pBQBInterface->printf("  -tables %d\n",gc_pParmFile->GetTableCount());
							continue;
						}

						if (strcmp(argv[z],"-opttables") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-opttables\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmFile->SetOptTables(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmFile->SetOptTables(false);
								z++;
							} else
								gc_pParmFile->SetOptTables(true);
							gc_pBQBInterface->printf("  -opttables %s\n",gc_pParmFile->GetOptTables()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-block") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-block\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 31)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-block\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 255.\n");
								return false;
							}
							gc_pParmFile->SetBlockLength(ti);
							gc_pBQBInterface->printf("  -block %d\n",gc_pParmFile->GetBlockLength());
							continue;
						}

						if (strcmp(argv[z],"-bw") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-bw\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmFile->SetBW(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmFile->SetBW(false);
								z++;
							} else
								gc_pParmFile->SetBW(true);
							gc_pBQBInterface->printf("  -bw %s\n",gc_pParmFile->GetBW()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-mtf") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-mtf\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmFile->SetMTF(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmFile->SetMTF(false);
								z++;
							} else
								gc_pParmFile->SetMTF(true);
							gc_pBQBInterface->printf("  -mtf %s\n",gc_pParmFile->GetMTF()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-rle") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-rle\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmFile->SetRLE(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmFile->SetRLE(false);
								z++;
							} else
								gc_pParmFile->SetRLE(true);
							gc_pBQBInterface->printf("  -rle %s\n",gc_pParmFile->GetRLE()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-verbose") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
								return false;
							}
							if (strcmp(argv[z+1],"2") == 0) {
								gc_iFileVerbose = 2;
								z++;
							} else if (strcmp(argv[z+1],"1") == 0) {
								gc_iFileVerbose = 1;
								z++;
							} else if (strcmp(argv[z+1],"0") == 0) {
								gc_iFileVerbose = 0;
								z++;
							} else
								gc_iFileVerbose = 1;
							gc_pBQBInterface->printf("  -verbose %d\n",gc_iFileVerbose);
							continue;
						}

						if (strcmp(argv[z],"-iter") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-iter\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 126)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-iter\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 126.\n");
								return false;
							}
							gc_pParmFile->SetMaxIter(ti);
							gc_pBQBInterface->printf("  -iter %d\n",gc_pParmFile->GetMaxIter());
							continue;
						}

						if (strcmp(argv[z],"-maxchunk") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-maxchunk\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 16777216)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-maxchunk\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 16777216.\n");
								return false;
							}
							gc_pParmFile->SetMaxChunk(ti);
							gc_pBQBInterface->printf("  -maxchunk %d\n",gc_pParmFile->GetMaxChunk());
							continue;
						}

						gc_pBQBInterface->eprintf("Error: Unknown command line argument \"%s\".\n",argv[z]);
						return false;

					} else {

						b = false;
						if (i == 0) {
							gc_sInFile = argv[z];
							gc_pBQBInterface->printf("  infile %s\n",gc_sInFile);
						} else if (i == 1) {
							gc_sOutFile = argv[z];
							gc_pBQBInterface->printf("  outfile %s\n",gc_sOutFile);
						} else {
							gc_pBQBInterface->eprintf("Error: Unrecognized trailing argument: \"%s\".\n",argv[z]);
							return false;
						}
						i++;
					}
				} // END IF FILE

				if (gc_bCube) {

					if ((argv[z][0] == '-') && b) {

						if (strcmp(argv[z],"-dryrun") == 0) {
							if (z+1 != argc) {
								if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
									gc_bDryRun = true;
									z++;
								} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
									gc_bDryRun = false;
									z++;
								} else
									gc_bDryRun = true;
							} else
								gc_bDryRun = true;
							gc_pBQBInterface->printf("  -dryrun %s\n",gc_bDryRun?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-pkey") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pkey\".\n");
								return false;
							}
							z++;
							if (!gc_pParmVol->FromKey(argv[z])) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-pkey\": \"%s\".\n",argv[z]);
								return false;
							}
							gc_pBQBInterface->printf("  -pkey %s\n",argv[z]);
							gc_iCubeOptimize = 0;
							continue;
						}

						if (strcmp(argv[z],"-start") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-start\".\n");
								return false;
							}
							z++;
							gc_iCubeStart = atoi(argv[z]);
							if (gc_iCubeStart < 1) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-start\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are >= 1.\n");
								return false;
							}
							gc_iCubeStart -= 1;
							gc_pBQBInterface->printf("  -start %d\n",gc_iCubeStart+1);
							continue;
						}

						if (strcmp(argv[z],"-steps") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-steps\".\n");
								return false;
							}
							z++;
							gc_iCubeSteps = atoi(argv[z]);
							if ((gc_iCubeSteps < -1) || (gc_iCubeSteps == 0)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-steps\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are -1 (=all) or > 0.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -steps %d\n",gc_iCubeSteps);
							continue;
						}

						if (strcmp(argv[z],"-stride") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-stride\".\n");
								return false;
							}
							z++;
							gc_iCubeStride = atoi(argv[z]);
							if (gc_iCubeStride < 1) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-stride\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are > 0.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -stride %d\n",gc_iCubeStride);
							continue;
						}

						if (strcmp(argv[z],"-keyframe") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-keyframe\".\n");
								return false;
							}
							z++;
							gc_iCubeKey = atoi(argv[z]);
							if (gc_iCubeKey < 1) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-keyframe\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are >= 1.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -keyframe %d\n",gc_iCubeKey);
							continue;
						}

						if (strcmp(argv[z],"-comment") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-comment\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_bCubeKeepComment = true;
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_bCubeKeepComment = false;
								z++;
							} else
								gc_bCubeKeepComment = true;
							gc_pBQBInterface->printf("  -comment %s\n",gc_bCubeKeepComment?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-check") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-check\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_bCubeCompare = true;
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_bCubeCompare = false;
								z++;
							} else
								gc_bCubeCompare = true;
							gc_pBQBInterface->printf("  -check %s\n",gc_bCubeCompare?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-dummyread") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-dummyread\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_bCubeDummyRead = true;
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_bCubeDummyRead = false;
								z++;
							} else
								gc_bCubeDummyRead = true;
							gc_pBQBInterface->printf("  -dummyread %s\n",gc_bCubeDummyRead?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-verbose") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
								return false;
							}
							if (strcmp(argv[z+1],"2") == 0) {
								gc_iCubeVerbose = 2;
								z++;
							} else if (strcmp(argv[z+1],"1") == 0) {
								gc_iCubeVerbose = 1;
								z++;
							} else if (strcmp(argv[z+1],"0") == 0) {
								gc_iCubeVerbose = 0;
								z++;
							} else
								gc_iCubeVerbose = 1;
							gc_pBQBInterface->printf("  -verbose %d\n",gc_iCubeVerbose);
							continue;
						}

						if (strcmp(argv[z],"-refout") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
								return false;
							}
							z++;
							gc_sCubeRefOut = argv[z];
							gc_pBQBInterface->printf("  -refout %s\n",gc_sCubeRefOut);
							continue;
						}

						if (strcmp(argv[z],"-vsigni") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vsigni\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 9)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vsigni\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 9.\n");
								return false;
							}
							gc_pParmVol->SetVolSigni(ti);
							gc_pBQBInterface->printf("  -vsigni %d\n",gc_pParmVol->GetVolSigni());
							continue;
						}

						if (strcmp(argv[z],"-veps") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-veps\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 63)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-veps\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 63.\n");
								return false;
							}
							gc_pParmVol->SetVolEps(ti);
							gc_pBQBInterface->printf("  -veps %d\n",gc_pParmVol->GetVolEps());
							continue;
						}

						if (strcmp(argv[z],"-vorder") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vorder\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vorder\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetVolOrder(ti);
							gc_pBQBInterface->printf("  -vorder %d\n",gc_pParmVol->GetVolOrder());
							continue;
						}

						if (strcmp(argv[z],"-voptorder") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-voptorder\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolOptOrder(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolOptOrder(false);
								z++;
							} else
								gc_pParmVol->SetVolOptOrder(true);
							gc_pBQBInterface->printf("  -voptorder %s\n",gc_pParmVol->GetVolOptOrder()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vhilbert") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vhilbert\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolHilbert(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolHilbert(false);
								z++;
							} else
								gc_pParmVol->SetVolHilbert(true);
							gc_pBQBInterface->printf("  -vhilbert %s\n",gc_pParmVol->GetVolHilbert()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vnbhfac") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vnbhfac\".\n");
								return false;
							}
							z++;
							tf = atof(argv[z]);
							if ((tf < 0.0) || (tf > 4.0)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vnbhfac\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0.0 .. 4.0\n");
								return false;
							}
							gc_pParmVol->SetVolNbhFac(tf);
							gc_pBQBInterface->printf("  -vnbhfac %.3f\n",gc_pParmVol->GetVolNbhFac());
							continue;
						}

						if (strcmp(argv[z],"-vsplit") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vsplit\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 31)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vsplit\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 31.\n");
								return false;
							}
							gc_pParmVol->SetVolSplit(ti);
							gc_pBQBInterface->printf("  -vsplit %d\n",gc_pParmVol->GetVolSplit());
							continue;
						}

						if (strcmp(argv[z],"-vtables") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vtables\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 64)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vtables\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 64.\n");
								return false;
							}
							gc_pParmVol->SetVolTableCount(ti);
							gc_pBQBInterface->printf("  -vtables %d\n",gc_pParmVol->GetVolTableCount());
							continue;
						}

						if (strcmp(argv[z],"-vopttables") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vopttables\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolOptTables(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolOptTables(false);
								z++;
							} else
								gc_pParmVol->SetVolOptTables(true);
							gc_pBQBInterface->printf("  -vopttables %s\n",gc_pParmVol->GetVolOptTables()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vblock") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vblock\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 128)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vblock\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 128.\n");
								return false;
							}
							gc_pParmVol->SetVolBlockLength(ti);
							gc_pBQBInterface->printf("  -vblock %d\n",gc_pParmVol->GetVolBlockLength());
							continue;
						}

						if (strcmp(argv[z],"-vbw") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vbw\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolBW(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolBW(false);
								z++;
							} else
								gc_pParmVol->SetVolBW(true);
							gc_pBQBInterface->printf("  -vbw %s\n",gc_pParmVol->GetVolBW()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vmtf") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vmtf\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolMTF(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolMTF(false);
								z++;
							} else
								gc_pParmVol->SetVolMTF(true);
							gc_pBQBInterface->printf("  -vmtf %s\n",gc_pParmVol->GetVolMTF()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vrle") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vrle\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolRLE(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolRLE(false);
								z++;
							} else
								gc_pParmVol->SetVolRLE(true);
							gc_pBQBInterface->printf("  -vrle %s\n",gc_pParmVol->GetVolRLE()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-viter") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-viter\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 126)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-viter\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 126.\n");
								return false;
							}
							gc_pParmVol->SetVolMaxIter(ti);
							gc_pBQBInterface->printf("  -viter %d\n",gc_pParmVol->GetVolMaxIter());
							continue;
						}

						if (strcmp(argv[z],"-vmaxchunk") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vmaxchunk\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 16777216)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vmaxchunk\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 16777216.\n");
								return false;
							}
							gc_pParmVol->SetVolMaxChunk(ti);
							gc_pBQBInterface->printf("  -vmaxchunk %d\n",gc_pParmVol->GetVolMaxChunk());
							continue;
						}

						if (strcmp(argv[z],"-optimize") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-optimize\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"0") == 0) || (strcmp(argv[z+1],"1") == 0) || (strcmp(argv[z+1],"2") == 0) ||
								(strcmp(argv[z+1],"3") == 0) || (strcmp(argv[z+1],"4") == 0)) {
								gc_iCubeOptimize = atoi(argv[z+1]);
								z++;
							} else if (bqbisinteger(argv[z+1])) {
								gc_pBQBInterface->eprintf("Error: Invalid value \"%s\" after \"-optimize\". Allowed range is 0 to 4.\n",
									argv[z+1]);
								return false;
							} else
								gc_iCubeOptimize = 2;
							if ((gc_iCubeOptimize < 0) || (gc_iCubeOptimize > 4)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-optimize\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 4.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -optimize %d\n",gc_iCubeOptimize);
							gc_bOptTouched = true;
							continue;
						}

						if (strcmp(argv[z],"-optsteps") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-optsteps\".\n");
								return false;
							}
							z++;
							gc_iCubeOptSteps = atoi(argv[z]);
							if (gc_iCubeOptSteps < 1) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-optsteps\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are > 0.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -optsteps %d\n",gc_iCubeOptSteps);
							continue;
						}

						if (strcmp(argv[z],"-onlyopt") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-onlyopt\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_bCubeOnlyOpt = true;
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_bCubeOnlyOpt = false;
								z++;
							} else
								gc_bCubeOnlyOpt = true;
							gc_pBQBInterface->printf("  -onlyopt %s\n",gc_bCubeOnlyOpt?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-optstart") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-optstart\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_iCubeOptIncludeFirst = 1;
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_iCubeOptIncludeFirst = 0;
								z++;
							} else
								gc_iCubeOptIncludeFirst = 1;
							gc_pBQBInterface->printf("  -optstart %s\n",(gc_iCubeOptIncludeFirst>0)?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vextra") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vextra\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolUseExtra(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolUseExtra(false);
								z++;
							} else
								gc_pParmVol->SetVolUseExtra(true);
							gc_pBQBInterface->printf("  -vextra %s\n",gc_pParmVol->GetVolUseExtra()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vexsrange") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexsrange\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vexsrange\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetVolExtraSRange(ti);
							gc_pBQBInterface->printf("  -vexsrange %d\n",gc_pParmVol->GetVolExtraSRangeX());
							continue;
						}

						if (strcmp(argv[z],"-vexsrangex") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexsrangex\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vexsrangex\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetVolExtraSRangeX(ti);
							gc_pBQBInterface->printf("  -vexsrangex %d\n",gc_pParmVol->GetVolExtraSRangeX());
							continue;
						}

						if (strcmp(argv[z],"-vexsrangey") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexsrangey\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vexsrangey\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetVolExtraSRangeY(ti);
							gc_pBQBInterface->printf("  -vexsrangey %d\n",gc_pParmVol->GetVolExtraSRangeY());
							continue;
						}

						if (strcmp(argv[z],"-vexsrangez") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexsrangez\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vexsrangez\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetVolExtraSRangeZ(ti);
							gc_pBQBInterface->printf("  -vexsrangez %d\n",gc_pParmVol->GetVolExtraSRangeZ());
							continue;
						}

						if (strcmp(argv[z],"-vextrange") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vextrange\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vextrange\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetVolExtraTRange(ti);
							gc_pBQBInterface->printf("  -vextrange %d\n",gc_pParmVol->GetVolExtraTRange());
							continue;
						}

						if (strcmp(argv[z],"-vexsorder") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexsorder\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vexsorder\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetVolExtraSOrder(ti);
							gc_pBQBInterface->printf("  -vexsorder %d\n",gc_pParmVol->GetVolExtraSOrder());
							continue;
						}

						if (strcmp(argv[z],"-vextorder") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vextorder\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vextorder\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetVolExtraTOrder(ti);
							gc_pBQBInterface->printf("  -vextorder %d\n",gc_pParmVol->GetVolExtraTOrder());
							continue;
						}

						if (strcmp(argv[z],"-vexoffset") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexoffset\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vexoffset\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetVolExtraOffset(ti);
							gc_pBQBInterface->printf("  -vexoffset %d\n",gc_pParmVol->GetVolExtraOffsetX());
							continue;
						}

						if (strcmp(argv[z],"-vexoffsetx") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexoffsetx\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vexoffsetx\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 15?.\n");
								return false;
							}
							gc_pParmVol->SetVolExtraOffsetX(ti);
							gc_pBQBInterface->printf("  -vexoffsetx %d\n",gc_pParmVol->GetVolExtraOffsetX());
							continue;
						}

						if (strcmp(argv[z],"-vexoffsety") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexoffsety\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vexoffsety\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetVolExtraOffsetY(ti);
							gc_pBQBInterface->printf("  -vexoffsety %d\n",gc_pParmVol->GetVolExtraOffsetY());
							continue;
						}

						if (strcmp(argv[z],"-vexoffsetz") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexoffsetz\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vexoffsetz\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetVolExtraOffsetZ(ti);
							gc_pBQBInterface->printf("  -vexoffsetz %d\n",gc_pParmVol->GetVolExtraOffsetZ());
							continue;
						}

						if (strcmp(argv[z],"-vexscross") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexscross\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolExtraCrossS(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolExtraCrossS(false);
								z++;
							} else
								gc_pParmVol->SetVolExtraCrossS(true);
							gc_pBQBInterface->printf("  -vexscross %s\n",gc_pParmVol->GetVolExtraCrossS()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vextcross") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vextcross\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolExtraCrossT(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolExtraCrossT(false);
								z++;
							} else
								gc_pParmVol->SetVolExtraCrossT(true);
							gc_pBQBInterface->printf("  -vextcross %s\n",gc_pParmVol->GetVolExtraCrossT()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vexwrap") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexwrap\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolExtraWrap(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolExtraWrap(false);
								z++;
							} else
								gc_pParmVol->SetVolExtraWrap(true);
							gc_pBQBInterface->printf("  -vexwrap %s\n",gc_pParmVol->GetVolExtraWrap()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vexscrossrange") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexscrossrange\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolExtraCrossRangeS(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolExtraCrossRangeS(false);
								z++;
							} else
								gc_pParmVol->SetVolExtraCrossRangeS(true);
							gc_pBQBInterface->printf("  -vexscrossrange %s\n",gc_pParmVol->GetVolExtraCrossRangeS()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vextcrossrange") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vextcrossrange\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolExtraCrossRangeT(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolExtraCrossRangeT(false);
								z++;
							} else
								gc_pParmVol->SetVolExtraCrossRangeT(true);
							gc_pBQBInterface->printf("  -vextcrossrange %s\n",gc_pParmVol->GetVolExtraCrossRangeT()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-vexdistexpo") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexdistexpo\".\n");
								return false;
							}
							z++;
							tf = atof(argv[z]);
							if ((tf < 0.0) || (tf > 16.0)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vexdistexpo\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0.0 .. 16.0\n");
								return false;
							}
							gc_pParmVol->SetVolExtraDistExpo(tf);
							gc_pBQBInterface->printf("  -vexdistexpo %.3f\n",gc_pParmVol->GetVolExtraDistExpo());
							continue;
						}

						if (strcmp(argv[z],"-vextimeexpo") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vextimeexpo\".\n");
								return false;
							}
							z++;
							tf = atof(argv[z]);
							if ((tf < 0.0) || (tf > 16.0)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-vextimeexpo\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0.0 .. 16.0\n");
								return false;
							}
							gc_pParmVol->SetVolExtraTimeExpo(tf);
							gc_pBQBInterface->printf("  -vextimeexpo %.3f\n",gc_pParmVol->GetVolExtraTimeExpo());
							continue;
						}

						if (strcmp(argv[z],"-vexpredcorr") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-vexpredcorr\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetVolExtraPredCorr(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetVolExtraPredCorr(false);
								z++;
							} else
								gc_pParmVol->SetVolExtraPredCorr(true);
							gc_pBQBInterface->printf("  -vexpredcorr %s\n",gc_pParmVol->GetVolExtraPredCorr()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-pextra") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pextra\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetPosUseExtra(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetPosUseExtra(false);
								z++;
							} else
								gc_pParmVol->SetPosUseExtra(true);
							gc_pBQBInterface->printf("  -pextra %s\n",gc_pParmVol->GetPosUseExtra()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-pextrange") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pextrange\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-pextrange\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetPosExtraTRange(ti);
							gc_pBQBInterface->printf("  -pextrange %d\n",gc_pParmVol->GetPosExtraTRange());
							continue;
						}

						if (strcmp(argv[z],"-pextorder") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pextorder\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-pextorder\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetPosExtraTOrder(ti);
							gc_pBQBInterface->printf("  -pextorder %d\n",gc_pParmVol->GetPosExtraTOrder());
							continue;
						}

						if (strcmp(argv[z],"-pextimeexpo") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pextimeexpo\".\n");
								return false;
							}
							z++;
							tf = atof(argv[z]);
							if ((tf < 0.0) || (tf > 16.0)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-pextimeexpo\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0.0 .. 16.0\n");
								return false;
							}
							gc_pParmVol->SetPosExtraTimeExpo(tf);
							gc_pBQBInterface->printf("  -pextimeexpo %.3f\n",gc_pParmVol->GetPosExtraTimeExpo());
							continue;
						}

						if (strcmp(argv[z],"-pprec") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pprec\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 7)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-pprec\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 7.\n");
								return false;
							}
							gc_pParmVol->SetPosPrecision(ti);
							gc_pBQBInterface->printf("  -pprec %d\n",gc_pParmVol->GetPosPrecision());
							continue;
						}

						if (strcmp(argv[z],"-porder") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-porder\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-porder\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 15.\n");
								return false;
							}
							gc_pParmVol->SetPosOrder(ti);
							gc_pBQBInterface->printf("  -porder %d\n",gc_pParmVol->GetPosOrder());
							continue;
						}

						if (strcmp(argv[z],"-poptorder") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-poptorder\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetPosOptOrder(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetPosOptOrder(false);
								z++;
							} else
								gc_pParmVol->SetPosOptOrder(true);
							gc_pBQBInterface->printf("  -poptorder %s\n",gc_pParmVol->GetPosOptOrder()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-psplit") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-psplit\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 31)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-psplit\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 31.\n");
								return false;
							}
							gc_pParmVol->SetPosSplit(ti);
							gc_pBQBInterface->printf("  -psplit %d\n",gc_pParmVol->GetPosSplit());
							continue;
						}

						if (strcmp(argv[z],"-ptables") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-ptables\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 64)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-ptables\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 64.\n");
								return false;
							}
							gc_pParmVol->SetPosTableCount(ti);
							gc_pBQBInterface->printf("  -ptables %d\n",gc_pParmVol->GetPosTableCount());
							continue;
						}

						if (strcmp(argv[z],"-popttables") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-popttables\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetPosOptTables(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetPosOptTables(false);
								z++;
							} else
								gc_pParmVol->SetPosOptTables(true);
							gc_pBQBInterface->printf("  -popttables %s\n",gc_pParmVol->GetPosOptTables()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-pblock") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pblock\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 128)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-pblock\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 128.\n");
								return false;
							}
							gc_pParmVol->SetPosBlockLength(ti);
							gc_pBQBInterface->printf("  -pblock %d\n",gc_pParmVol->GetPosBlockLength());
							continue;
						}

						if (strcmp(argv[z],"-pbw") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pbw\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetPosBW(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetPosBW(false);
								z++;
							} else
								gc_pParmVol->SetPosBW(true);
							gc_pBQBInterface->printf("  -pbw %s\n",gc_pParmVol->GetPosBW()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-pmtf") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pmtf\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetPosMTF(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetPosMTF(false);
								z++;
							} else
								gc_pParmVol->SetPosMTF(true);
							gc_pBQBInterface->printf("  -pmtf %s\n",gc_pParmVol->GetPosMTF()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-prle") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-prle\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetPosRLE(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetPosRLE(false);
								z++;
							} else
								gc_pParmVol->SetPosRLE(true);
							gc_pBQBInterface->printf("  -prle %s\n",gc_pParmVol->GetPosRLE()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-piter") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-piter\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 126)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-piter\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 126.\n");
								return false;
							}
							gc_pParmVol->SetPosMaxIter(ti);
							gc_pBQBInterface->printf("  -piter %d\n",gc_pParmVol->GetPosMaxIter());
							continue;
						}

						if (strcmp(argv[z],"-psortatom") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-psortatom\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmVol->SetPosSortAtom(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmVol->SetPosSortAtom(false);
								z++;
							} else
								gc_pParmVol->SetPosSortAtom(true);
							gc_pBQBInterface->printf("  -psortatom %s\n",gc_pParmVol->GetPosSortAtom()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-pmaxchunk") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pmaxchunk\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 16777216)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-pmaxchunk\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 16777216.\n");
								return false;
							}
							gc_pParmVol->SetPosMaxChunk(ti);
							gc_pBQBInterface->printf("  -pmaxchunk %d\n",gc_pParmVol->GetPosMaxChunk());
							continue;
						}

						gc_pBQBInterface->eprintf("Error: Unknown command line argument \"%s\".\n",argv[z]);
						return false;

					} else {

						b = false;
						if (i == 0) {
							gc_sInFile = argv[z];
							gc_pBQBInterface->printf("  infile %s\n",gc_sInFile);
						} else if (i == 1) {
							gc_sOutFile = argv[z];
							gc_pBQBInterface->printf("  outfile %s\n",gc_sOutFile);
						} else {
							gc_pBQBInterface->eprintf("Error: Unrecognized trailing argument: \"%s\".\n",argv[z]);
							return false;
						}
						i++;
					}
				} // END IF CUBE

				if (gc_bXYZ) {

					if ((argv[z][0] == '-') && b) {

						if (strcmp(argv[z],"-dryrun") == 0) {
							if (z+1 != argc) {
								if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
									gc_bDryRun = true;
									z++;
								} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
									gc_bDryRun = false;
									z++;
								} else
									gc_bDryRun = true;
							} else
								gc_bDryRun = true;
							gc_pBQBInterface->printf("  -dryrun %s\n",gc_bDryRun?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-pkey") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-pkey\".\n");
								return false;
							}
							z++;
							if (!gc_pParmPos->FromKey(argv[z])) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-pkey\": \"%s\".\n",argv[z]);
								return false;
							}
							gc_pBQBInterface->printf("  -pkey %s\n",argv[z]);
							gc_iAtomOptimize = 0;
							continue;
						}

						if (strcmp(argv[z],"-start") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-start\".\n");
								return false;
							}
							z++;
							gc_iAtomStart = atoi(argv[z]);
							if (gc_iAtomStart < 1) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-start\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are > 1.\n");
								return false;
							}
							gc_iAtomStart -= 1;
							gc_pBQBInterface->printf("  -start %d\n",gc_iAtomStart+1);
							continue;
						}

						if (strcmp(argv[z],"-steps") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-steps\".\n");
								return false;
							}
							z++;
							gc_iAtomSteps = atoi(argv[z]);
							if ((gc_iAtomSteps < -1) || (gc_iAtomSteps == 0)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-steps\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are -1 (=all) or > 0.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -steps %d\n",gc_iAtomSteps);
							continue;
						}

						if (strcmp(argv[z],"-stride") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-stride\".\n");
								return false;
							}
							z++;
							gc_iAtomStride = atoi(argv[z]);
							if (gc_iAtomStride < 1) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-stride\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are > 0.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -stride %d\n",gc_iAtomStride);
							continue;
						}

						if (strcmp(argv[z],"-keyframe") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-keyframe\".\n");
								return false;
							}
							z++;
							gc_iAtomKey = atoi(argv[z]);
							if (gc_iAtomKey < 1) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-keyframe\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are >= 1.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -keyframe %d\n",gc_iAtomKey);
							continue;
						}

						if (strcmp(argv[z],"-comment") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-comment\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_bAtomKeepComment = true;
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_bAtomKeepComment = false;
								z++;
							} else
								gc_bAtomKeepComment = true;
							gc_pBQBInterface->printf("  -comment %s\n",gc_bAtomKeepComment?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-check") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-check\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_bAtomCompare = true;
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_bAtomCompare = false;
								z++;
							} else
								gc_bAtomCompare = true;
							gc_pBQBInterface->printf("  -check %s\n",gc_bAtomCompare?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-verbose") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
								return false;
							}
							if (strcmp(argv[z+1],"2") == 0) {
								gc_iAtomVerbose = 2;
								z++;
							} else if (strcmp(argv[z+1],"1") == 0) {
								gc_iAtomVerbose = 1;
								z++;
							} else if (strcmp(argv[z+1],"0") == 0) {
								gc_iAtomVerbose = 0;
								z++;
							} else
								gc_iAtomVerbose = 1;
							gc_pBQBInterface->printf("  -verbose %d\n",gc_iAtomVerbose);
							continue;
						}

						if (strcmp(argv[z],"-refout") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
								return false;
							}
							z++;
							gc_sAtomRefOut = argv[z];
							gc_pBQBInterface->printf("  -refout %s\n",gc_sAtomRefOut);
							continue;
						}

						if (strcmp(argv[z],"-prec") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-prec\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 7)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-prec\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 7.\n");
								return false;
							}
							gc_pParmPos->SetPrecision(ti);
							gc_pBQBInterface->printf("  -prec %d\n",gc_pParmPos->GetPrecision());
							continue;
						}

						if (strcmp(argv[z],"-order") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-order\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-order\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 15.\n");
								return false;
							}
							gc_pParmPos->SetOrder(ti);
							gc_pBQBInterface->printf("  -order %d\n",gc_pParmPos->GetOrder());
							continue;
						}

						if (strcmp(argv[z],"-optorder") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-optorder\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmPos->SetOptOrder(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmPos->SetOptOrder(false);
								z++;
							} else
								gc_pParmPos->SetOptOrder(true);
							gc_pBQBInterface->printf("  -optorder %s\n",gc_pParmPos->GetOptOrder()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-split") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-split\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 31)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-split\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 31.\n");
								return false;
							}
							gc_pParmPos->SetSplit(ti);
							gc_pBQBInterface->printf("  -split %d\n",gc_pParmPos->GetSplit());
							continue;
						}

						if (strcmp(argv[z],"-tables") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-tables\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 64)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-tables\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 64.\n");
								return false;
							}
							gc_pParmPos->SetTableCount(ti);
							gc_pBQBInterface->printf("  -tables %d\n",gc_pParmPos->GetTableCount());
							continue;
						}

						if (strcmp(argv[z],"-opttables") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-opttables\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmPos->SetOptTables(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmPos->SetOptTables(false);
								z++;
							} else
								gc_pParmPos->SetOptTables(true);
							gc_pBQBInterface->printf("  -opttables %s\n",gc_pParmPos->GetOptTables()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-block") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-block\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 128)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-block\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 128.\n");
								return false;
							}
							gc_pParmPos->SetBlockLength(ti);
							gc_pBQBInterface->printf("  -block %d\n",gc_pParmPos->GetBlockLength());
							continue;
						}

						if (strcmp(argv[z],"-bw") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-bw\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmPos->SetBW(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmPos->SetBW(false);
								z++;
							} else
								gc_pParmPos->SetBW(true);
							gc_pBQBInterface->printf("  -bw %s\n",gc_pParmPos->GetBW()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-mtf") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-mtf\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmPos->SetMTF(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmPos->SetMTF(false);
								z++;
							} else
								gc_pParmPos->SetMTF(true);
							gc_pBQBInterface->printf("  -mtf %s\n",gc_pParmPos->GetMTF()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-rle") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-rle\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmPos->SetRLE(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmPos->SetRLE(false);
								z++;
							} else
								gc_pParmPos->SetRLE(true);
							gc_pBQBInterface->printf("  -rle %s\n",gc_pParmPos->GetRLE()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-iter") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-iter\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 126)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-iter\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 126.\n");
								return false;
							}
							gc_pParmPos->SetMaxIter(ti);
							gc_pBQBInterface->printf("  -iter %d\n",gc_pParmPos->GetMaxIter());
							continue;
						}

						if (strcmp(argv[z],"-sortatom") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-sortatom\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmPos->SetSortAtom(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmPos->SetSortAtom(false);
								z++;
							} else
								gc_pParmPos->SetSortAtom(true);
							gc_pBQBInterface->printf("  -sortatom %s\n",gc_pParmPos->GetSortAtom()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-maxchunk") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-maxchunk\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 16777216)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-maxchunk\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 16777216.\n");
								return false;
							}
							gc_pParmPos->SetMaxChunk(ti);
							gc_pBQBInterface->printf("  -maxchunk %d\n",gc_pParmPos->GetMaxChunk());
							continue;
						}

						if (strcmp(argv[z],"-extra") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-extra\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_pParmPos->SetUseExtra(true);
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_pParmPos->SetUseExtra(false);
								z++;
							} else
								gc_pParmPos->SetUseExtra(true);
							gc_pBQBInterface->printf("  -extra %s\n",gc_pParmPos->GetUseExtra()?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-extrange") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-extrange\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 1) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-extrange\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 1 .. 15.\n");
								return false;
							}
							gc_pParmPos->SetExtraTRange(ti);
							gc_pBQBInterface->printf("  -extrange %d\n",gc_pParmPos->GetExtraTRange());
							continue;
						}

						if (strcmp(argv[z],"-extorder") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-extorder\".\n");
								return false;
							}
							z++;
							ti = atoi(argv[z]);
							if ((ti < 0) || (ti > 15)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-extorder\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 15.\n");
								return false;
							}
							gc_pParmPos->SetExtraTOrder(ti);
							gc_pBQBInterface->printf("  -extorder %d\n",gc_pParmPos->GetExtraTOrder());
							continue;
						}

						if (strcmp(argv[z],"-extimeexpo") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-extimeexpo\".\n");
								return false;
							}
							z++;
							tf = atof(argv[z]);
							if ((tf < 0) || (tf > 16.0)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-extimeexpo\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0.0 .. 16.0\n");
								return false;
							}
							gc_pParmPos->SetExtraTimeExpo(tf);
							gc_pBQBInterface->printf("  -extimeexpo %.3f\n",gc_pParmPos->GetExtraTimeExpo());
							continue;
						}

						if (strcmp(argv[z],"-optimize") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-optimize\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"0") == 0) || (strcmp(argv[z+1],"1") == 0) || (strcmp(argv[z+1],"2") == 0) ||
								(strcmp(argv[z+1],"3") == 0) || (strcmp(argv[z+1],"4") == 0)) {
								gc_iAtomOptimize = atoi(argv[z+1]);
								z++;
							} else if (bqbisinteger(argv[z+1])) {
								gc_pBQBInterface->eprintf("Error: Invalid value \"%s\" after \"-optimize\". Allowed range is 0 to 3.\n",
									argv[z+1]);
								return false;
							} else
								gc_iAtomOptimize = 3;
							if ((gc_iAtomOptimize < 0) || (gc_iAtomOptimize > 4)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-optimize\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are 0 .. 4.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -optimize %d\n",gc_iAtomOptimize);
							gc_bOptTouched = true;
							continue;
						}

						if (strcmp(argv[z],"-optsteps") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-optsteps\".\n");
								return false;
							}
							z++;
							gc_iAtomOptSteps = atoi(argv[z]);
							if (gc_iAtomOptSteps < 1) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-optsteps\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are > 0.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -optsteps %d\n",gc_iAtomOptSteps);
							continue;
						}

						if (strcmp(argv[z],"-onlyopt") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-onlyopt\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_bAtomOnlyOpt = true;
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_bAtomOnlyOpt = false;
								z++;
							} else
								gc_bAtomOnlyOpt = true;
							gc_pBQBInterface->printf("  -onlyopt %s\n",gc_bAtomOnlyOpt?"yes":"no");
							continue;
						}

						if (strcmp(argv[z],"-optstart") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-optstart\".\n");
								return false;
							}
							if ((strcmp(argv[z+1],"yes") == 0) || (strcmp(argv[z+1],"y") == 0)) {
								gc_iAtomOptIncludeFirst = 1;
								z++;
							} else if ((strcmp(argv[z+1],"no") == 0) || (strcmp(argv[z+1],"n") == 0)) {
								gc_iAtomOptIncludeFirst = 0;
								z++;
							} else
								gc_iAtomOptIncludeFirst = 1;
							gc_pBQBInterface->printf("  -optstart %s\n",(gc_iAtomOptIncludeFirst>0)?"yes":"no");
							continue;
						}

						gc_pBQBInterface->eprintf("Error: Unknown command line argument \"%s\".\n",argv[z]);
						return false;

					} else {

						b = false;
						if (i == 0) {
							gc_sInFile = argv[z];
							gc_pBQBInterface->printf("  infile %s\n",gc_sInFile);
						} else if (i == 1) {
							gc_sOutFile = argv[z];
							gc_pBQBInterface->printf("  outfile %s\n",gc_sOutFile);
						} else {
							gc_pBQBInterface->eprintf("Error: Unrecognized trailing argument: \"%s\".\n",argv[z]);
							return false;
						}
						i++;
					}
				} // END IF XYZ

			} // END IF COMPRESS


			if (gc_bDecompress) {

				if (gc_bFile) {

					if ((argv[z][0] == '-') && b) {

						if (strcmp(argv[z],"-verbose") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
								return false;
							}
							if (strcmp(argv[z+1],"2") == 0) {
								gc_iFileVerbose = 2;
								z++;
							} else if (strcmp(argv[z+1],"1") == 0) {
								gc_iFileVerbose = 1;
								z++;
							} else if (strcmp(argv[z+1],"0") == 0) {
								gc_iFileVerbose = 0;
								z++;
							} else
								gc_iFileVerbose = 1;
							gc_pBQBInterface->printf("  -verbose %d\n",gc_iFileVerbose);
							continue;
						}

						if (strcmp(argv[z],"-refcompare") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-refcompare\".\n");
								return false;
							}
							z++;
							gc_sRefFile = argv[z];
							gc_pBQBInterface->printf("  -refcompare %s\n",gc_sRefFile);
							continue;
						}

						gc_pBQBInterface->eprintf("Error: Unknown command line argument \"%s\".\n",argv[z]);
						return false;

					} else {

						b = false;
						if (i == 0) {
							gc_sInFile = argv[z];
							gc_pBQBInterface->printf("  infile %s\n",gc_sInFile);
						} else if (i == 1) {
							gc_sOutFile = argv[z];
							gc_pBQBInterface->printf("  outfile %s\n",gc_sOutFile);
						} else {
							gc_pBQBInterface->eprintf("Error: Unrecognized trailing argument: \"%s\".\n",argv[z]);
							return false;
						}
						i++;
					}
				} // END IF FILE

				if (gc_bCube) {

					if ((argv[z][0] == '-') && b) {

						if (strcmp(argv[z],"-steps") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-steps\".\n");
								return false;
							}
							z++;
							gc_iCubeSteps = atoi(argv[z]);
							if ((gc_iCubeSteps < -1) || (gc_iCubeSteps == 0)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-steps\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are -1 (=all) or > 0.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -steps %d\n",gc_iCubeSteps);
							continue;
						}

						if (strcmp(argv[z],"-stride") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-stride\".\n");
								return false;
							}
							z++;
							gc_iCubeStride = atoi(argv[z]);
							if (gc_iCubeStride < 1) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-stride\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are > 0.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -stride %d\n",gc_iCubeStride);
							continue;
						}

						if (strcmp(argv[z],"-verbose") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
								return false;
							}
							if (strcmp(argv[z+1],"2") == 0) {
								gc_iCubeVerbose = 2;
								z++;
							} else if (strcmp(argv[z+1],"1") == 0) {
								gc_iCubeVerbose = 1;
								z++;
							} else if (strcmp(argv[z+1],"0") == 0) {
								gc_iCubeVerbose = 0;
								z++;
							} else
								gc_iCubeVerbose = 1;
							gc_pBQBInterface->printf("  -verbose %d\n",gc_iCubeVerbose);
							continue;
						}

						if (strcmp(argv[z],"-refcompare") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-refcompare\".\n");
								return false;
							}
							z++;
							gc_sRefFile = argv[z];
							gc_pBQBInterface->printf("  -refcompare %s\n",gc_sRefFile);
							continue;
						}

						gc_pBQBInterface->eprintf("Error: Unknown command line argument \"%s\".\n",argv[z]);
						return false;

					} else {

						b = false;
						if (i == 0) {
							gc_sInFile = argv[z];
							gc_pBQBInterface->printf("  infile %s\n",gc_sInFile);
						} else if (i == 1) {
							gc_sOutFile = argv[z];
							gc_pBQBInterface->printf("  outfile %s\n",gc_sOutFile);
						} else {
							gc_pBQBInterface->eprintf("Error: Unrecognized trailing argument: \"%s\".\n",argv[z]);
							return false;
						}
						i++;
					}
				} // END IF CUBE

				if (gc_bXYZ) {

					if ((argv[z][0] == '-') && b) {

						if (strcmp(argv[z],"-steps") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-steps\".\n");
								return false;
							}
							z++;
							gc_iAtomSteps = atoi(argv[z]);
							if ((gc_iAtomSteps < -1) || (gc_iAtomSteps == 0)) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-steps\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are -1 (=all) or > 0.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -steps %d\n",gc_iAtomSteps);
							continue;
						}

						if (strcmp(argv[z],"-stride") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-stride\".\n");
								return false;
							}
							z++;
							gc_iAtomStride = atoi(argv[z]);
							if (gc_iAtomStride < 1) {
								gc_pBQBInterface->eprintf("Error: Invalid value for \"-stride\": \"%s\".\n",argv[z]);
								gc_pBQBInterface->printf("       Allowed values are > 0.\n");
								return false;
							}
							gc_pBQBInterface->printf("  -stride %d\n",gc_iAtomStride);
							continue;
						}

						if (strcmp(argv[z],"-verbose") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-verbose\".\n");
								return false;
							}
							if (strcmp(argv[z+1],"2") == 0) {
								gc_iAtomVerbose = 2;
								z++;
							} else if (strcmp(argv[z+1],"1") == 0) {
								gc_iAtomVerbose = 1;
								z++;
							} else if (strcmp(argv[z+1],"0") == 0) {
								gc_iAtomVerbose = 0;
								z++;
							} else
								gc_iAtomVerbose = 1;
							gc_pBQBInterface->printf("  -verbose %d\n",gc_iAtomVerbose);
							continue;
						}

						if (strcmp(argv[z],"-refcompare") == 0) {
							if (z+1 == argc) {
								gc_pBQBInterface->eprintf("Error: Missing value after \"-refcompare\".\n");
								return false;
							}
							z++;
							gc_sRefFile = argv[z];
							gc_pBQBInterface->printf("  -refcompare %s\n",gc_sRefFile);
							continue;
						}

						gc_pBQBInterface->eprintf("Error: Unknown command line argument \"%s\".\n",argv[z]);
						return false;

					} else {

						b = false;
						if (i == 0) {
							gc_sInFile = argv[z];
							gc_pBQBInterface->printf("  infile %s\n",gc_sInFile);
						} else if (i == 1) {
							gc_sOutFile = argv[z];
							gc_pBQBInterface->printf("  outfile %s\n",gc_sOutFile);
						} else {
							gc_pBQBInterface->eprintf("Error: Unrecognized trailing argument: \"%s\".\n",argv[z]);
							return false;
						}
						i++;
					}
				} // END IF XYZ

			} // END IF DECOMPRESS

		}
	}

	gc_pBQBInterface->printf("Finished.\n\n");

	return true;
}


static const char  BQB_BASE64_table[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";


bool bqb_isbase64(char c) {
   return c && strchr(BQB_BASE64_table, c) != NULL;
}


inline char BQB_BASE64_value(char c) {
   const char *p = strchr(BQB_BASE64_table, c);
   if(p) {
      return (char)(p-BQB_BASE64_table);
   } else {
      return 0;
   }
}


int BQB_UnBase64(unsigned char *dest, const unsigned char *src, int srclen) {

   *dest = 0;
   if(*src == 0) 
   {
      return 0;
   }
   unsigned char *p = dest;
   do
   {

      char a = BQB_BASE64_value(src[0]);
      char b = BQB_BASE64_value(src[1]);
      char c = BQB_BASE64_value(src[2]);
      char d = BQB_BASE64_value(src[3]);
      *p++ = (a << 2) | (b >> 4);
      *p++ = (b << 4) | (c >> 2);
      *p++ = (c << 6) | d;
      if(!bqb_isbase64(src[1])) 
      {
         p -= 2;
         break;
      } 
      else if(!bqb_isbase64(src[2])) 
      {
         p -= 2;
         break;
      } 
      else if(!bqb_isbase64(src[3])) 
      {
         p--;
         break;
      }
      src += 4;
      while(*src && (*src == 13 || *src == 10)) src++;
   }
   while(srclen-= 4);
   *p = 0;
   return (int)(p-dest);
}


void BQB_PrintBMode() {

	const char *stext[15];
	char buf[256];
	int z;

	stext[0]   = "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIF8sIC4u";
	stext[1]   = "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICdfKF8wbylvKCku";
	stext[2]   = "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAsbyhfXyksXyksKF9fKW8=";
	stext[3]   = "ICAgICAgICAgXy4oLSkuXyAgICAgICAgICAgICAgICAgICAgbyhfLC1vKF8gKSgoXylvTylf";
	stext[4]   = "ICAgICAgLicgICAgICAgICAnLiAgICAgICAgICAgICAgICAgLk8oX18pbyxfXylvKF8pT29fKQ==";
	stext[5]   = "ICAgICAvICAgICAgICAgICAgIFwgICAgICAgICAgICAuLS0tLXwgICB8ICAgfCAgIHxfKTA=";
	stext[6]   = "ICAgICB8Jy0uLi5fX18uLi4tJ3wgICAgICAgICAgIC8gIC4tLXwgICB8ICAgfCAgIHwsXyk=";
	stext[7]   = "ICAgICAgXCAgICAnPScgICAgLyAgICAgICAgICAgfCAgLyAgIHwgICB8ICAgfCAgIHxvKF8p";
	stext[8]   = "ICAgICAgIGAnLl9fX19fLidgICAgICAgICAgICAgfCAgfCAgIHwgICB8ICAgfCAgIHxfL2Ap";
	stext[9]   = "ICAgICAgICAvICAgfCAgIFwgICAgICAgICAgICAgfCAgfCAgIHwgICB8ICAgfCAgIHxPXyk=";
	stext[10]  = "ICAgICAgIC8uLS0nfCctLS5cICAgICAgICAgICAgfCAgXCAgIHwgICB8ICAgfCAgIHw=";
	stext[11]  = "ICAgIFtdLyctLl9ffF9fLi0nXFtdICAgICAgICAgIFwgICctLXwgICB8ICAgfCAgIHw=";
	stext[12]  = "ICAgICAgICAgICAgfCAgICAgICAgICAgICAgICAgICAnLS0tLXwgICB8ICAgfCAgIHw=";
	stext[13]  = "ICAgICAgICAgICBbXSAgICAgICAgICAgICAgICAgICAgICAgIFwgICBcICAgLyAgIC8=";
	stext[14]  = "ICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICBgIiIiIiIiIiIiYA==";

	//gc_pBQBInterface->printf("\n");
	for (z=0;z<15;z++) {
		BQB_UnBase64((unsigned char*)buf,(const unsigned char*)stext[z],(int)strlen(stext[z]));
		gc_pBQBInterface->bprintf("    %s\n",buf);
	}
	gc_pBQBInterface->printf("\n");
}




// Functions for Linux signal handling / stack trace
#ifndef BQB_INSIDE_TRAVIS

	void bqb_printf(const char *s, ...) {
		va_list args;
		va_start(args,s);
		vprintf(s,args);
		va_end(args);
		fflush(stdout);
		if (gc_fLogFile != NULL) {
			va_start(args,s);
			vfprintf(gc_fLogFile,s,args);
			va_end(args);
			fflush(gc_fLogFile);
		}
	}

	#if !defined(_WIN32) && (defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__)))

		void UninstallSignalHandler() {

			signal(SIGSEGV,SIG_DFL);
			signal(SIGFPE,SIG_DFL);
			signal(SIGILL,SIG_DFL);
			signal(SIGINT,SIG_DFL);
			signal(SIGABRT,SIG_DFL);
			signal(SIGTERM,SIG_DFL);
			signal(SIGHUP,SIG_DFL);
		}


		void DumpLinuxBacktrace() {

			void *buffer[4096];
			char **strings, buf[256], *p;
			int n, z, i, j;
			bool showall;
			FILE *a;
			

			n = backtrace(buffer,4096);
			if (n <= 0) {
				bqb_printf("DumpLinuxBacktrace(): Error: Could not retrieve stack trace.\n");
				return;
			}
			strings = backtrace_symbols(buffer,n);
			if (strings == NULL)
				bqb_printf("Could not retrieve debug symbol strings.\n");
			bqb_printf("Decoding line numbers... ");
			bqb_printf("(\"addr2line <addr> -f -i -C -e %s\")\n",gc_sExeName);
			a = fopen(gc_sExeName,"rb");
			if (a == NULL) {
				bqb_printf("Could not locate \"%s\" :-(\n",gc_sExeName);
			} else {
				fclose(a);
				j = 0;
				for (z=0;z<n;z++)
					if ((buffer[z] < (void*)0xF00000) && (buffer[z] > (void*)0x400000))
						j++;
				if (j < 3)
					showall = true;
				else
					showall = false;
				for (z=0;z<n;z++) {
					if (showall || ((buffer[z] < (void*)0xF00000) && (buffer[z] > (void*)0x400000))) {
						bqb_printf(" # %2d ",n-z);
						sprintf(buf,"addr2line %p -f -i -C -e %s > tmp.txt",buffer[z],gc_sExeName);
						(void)!system(buf);
						a = fopen("tmp.txt","rt");
						if (a == NULL) {
							bqb_printf("Error. Make sure you have the addr2line tool installed (it is part of \"binutils\")\nand writing permission to this directory.\n");
							continue;
						}
						bqb_printf("- [%10p] ",buffer[z]);
						i = 0;
						while (true) {
							(void)!fgets(buf,256,a);
							if (feof(a))
								break;
							if (strlen(buf) > 2) {
								if (i == 2) {
									i = 0;
									bqb_printf("\n                     ");
								}
								buf[strlen(buf)-1] = 0;
								p = strrchr(buf,'/');
								if (p == NULL)
									p = buf;
										else p++;
								bqb_printf("- %s ",p);
								i++;
							}
						}
						fclose(a);
						bqb_printf("\n");
					}
				}
				(void)!system("rm tmp.txt");
			}
		}


		void Crash() {

			UninstallSignalHandler();
			bqb_printf("bqbtool apparently has crashed :-( Sorry.\n");
			bqb_printf("Please send the screen output to the developers,\n");
			bqb_printf("making them able to analyze and fix this error.\n");
			DumpLinuxBacktrace();
			bqb_printf("Delivering control to operating system.\n");
			abort();
		}


		void CrashAbort() {

			UninstallSignalHandler();
			bqb_printf("bqbtool had to abort execution.\n");
			bqb_printf("Hopefully, the reason was printed above this message.\n");
			DumpLinuxBacktrace();
			bqb_printf("Delivering control to operating system.\n");
			abort();
		}


		void CrashTerm() {

			UninstallSignalHandler();
			bqb_printf("bqbtool received a termination request (SIGTERM) from the operating system\n");
			bqb_printf("and therefore has to stop execution.\n");
			DumpLinuxBacktrace();
			bqb_printf("Delivering control to operating system.\n");
			abort();
		}


		void CrashInt() {

			UninstallSignalHandler();
			DumpLinuxBacktrace();
			bqb_printf("Delivering control to operating system.\n");
			abort();
		}


		void SIGNAL_SEGV(int param) {
			UNUSED(param);  // Suppress "unused parameter" warning
			bqb_printf("\n\n*** SIGSEGV caught: Segmentation fault ***\n");
			Crash();
		}


		void SIGNAL_FPE(int param) {
			UNUSED(param);  // Suppress "unused parameter" warning
			bqb_printf("\n\n*** SIGFPE caught: Floating point exception ***\n");
			Crash();
		}


		void SIGNAL_ILL(int param) {
			UNUSED(param);  // Suppress "unused parameter" warning
			bqb_printf("\n\n*** SIGILL caught: Illegal instruction ***\n");
			Crash();
		}


		void SIGNAL_INT(int param) {
			UNUSED(param);  // Suppress "unused parameter" warning
			bqb_printf("\n\n*** SIGINT caught: Hard interrupt by user ***\n");
			bqb_printf("Stopping execution.\n");
			CrashInt();
			exit(0);
		}


		void SIGNAL_ABRT(int param) {
			UNUSED(param);  // Suppress "unused parameter" warning
			bqb_printf("\n\n*** SIGABRT caught: Abnormal program termination ***\n");
			CrashAbort();
		}


		void SIGNAL_TERM(int param) {
			UNUSED(param);  // Suppress "unused parameter" warning
			bqb_printf("\n\n*** SIGTERM caught: Terminating ***\n");
			CrashTerm();
		}


		void SIGNAL_HANGUP(int param) {
			UNUSED(param);  // Suppress "unused parameter" warning
			bqb_printf("\n\n*** SIGHUP caught: Trying to ignore hangup ^^ ***\n");
			signal(SIGHUP,SIGNAL_HANGUP);
		}


		void InstallSignalHandler() {
			signal(SIGSEGV,SIGNAL_SEGV);
			signal(SIGFPE,SIGNAL_FPE);
			signal(SIGILL,SIGNAL_ILL);
			signal(SIGINT,SIGNAL_INT);
			signal(SIGABRT,SIGNAL_ABRT);
			signal(SIGTERM,SIGNAL_TERM);
			signal(SIGHUP,SIGNAL_HANGUP);
		}

	#endif

#endif





#ifdef BQB_INSIDE_TRAVIS

	int bqbtool_main(int argc, const char *argv[]) {

#else

	int main(int argc, const char *argv[]) {

#endif

	CBQBDriver *driver = NULL;
	int z;
	bool err, b;
	struct tm *today;
	time_t ltime;
	char buf[256];


	gc_sExeName = new char[strlen(argv[0])+1];
	strcpy(gc_sExeName,argv[0]);


	// Install Linux signal handlers
	#ifndef BQB_INSIDE_TRAVIS
		#if !defined(_WIN32) && (defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__)))
			InstallSignalHandler();
		#endif
	#endif


	// Before this command, no printing should occur
	gc_pBQBInterface = BQBCreateInterface( 0 );

	if (gc_pBQBInterface == NULL) {
		printf("Error: Could not create BQB interface.\n");
		printf("Aborting execution.\n\n");
		return 1;
	}

	
	// Use TRAVIS functions for screen output
	#ifdef BQB_INSIDE_TRAVIS
		gc_pBQBInterface->SetPrintCallback(  &mprintf );
		gc_pBQBInterface->SetBPrintCallback( &bprintf );
		gc_pBQBInterface->SetEPrintCallback( &eprintf );
	#else
		gc_fLogFile = fopen("bqbtool.log","wt");
		if (gc_fLogFile == NULL) {
			printf("Error: Could not open \"bqbtool.log\" for writing.\n");
			printf("Aborting execution.\n\n");
			return 1;
		}
		gc_pBQBInterface->SetPrintCallback(  &bqb_printf );
		gc_pBQBInterface->SetBPrintCallback( &bqb_printf );
		gc_pBQBInterface->SetEPrintCallback( &bqb_printf );
	#endif


	// Get host name and working directory
	#if !defined(_WIN32) && (defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__)))
		if (gethostname(buf,256)==0) {
			gc_sHostName = new char[strlen(buf)+1];
			strcpy(gc_sHostName,buf);
		} else
			gc_sHostName = NULL;
		gc_sWorkingDir = getcwd(NULL,1024);
	#elif defined(_WIN32)
		unsigned long l;
		l = 256;
		if (GetComputerNameA(buf,&l)) { 
			gc_sHostName = new char[strlen(buf)+1];
			strcpy(gc_sHostName,buf);
		} else
			gc_sHostName = NULL;
		gc_sWorkingDir = _getcwd(NULL,1024);
	#else
		gc_sHostName = NULL;
		gc_sWorkingDir = NULL;
	#endif


	gc_pBQBInterface->bprintf("\n");
	gc_pBQBInterface->bprintf("     __                    __           __                             __\n");
	gc_pBQBInterface->bprintf("    /  |                  /  |         /  |                           /  |\n");
	gc_pBQBInterface->bprintf("    ## |____     ______   ## |____    _## |_      ______     ______   ## |\n");
	gc_pBQBInterface->bprintf("    ##      \\   /      \\  ##      \\  / ##   |    /      \\   /      \\  ## |\n");
	gc_pBQBInterface->bprintf("    #######  | /######  | #######  | ######/    /######  | /######  | ## |\n");
	gc_pBQBInterface->bprintf("    ## |  ## | ## |  ## | ## |  ## |   ## | __  ## |  ## | ## |  ## | ## |\n");
	gc_pBQBInterface->bprintf("    ## |__## | ## \\__## | ## |__## |   ## |/  | ## \\__## | ## \\__## | ## |\n");
	gc_pBQBInterface->bprintf("    ##    ##/  ##    ## | ##    ##/    ##  ##/  ##    ##/  ##    ##/  ## \\\n");
	gc_pBQBInterface->bprintf("    #######/    ####### | #######/      ####/    ######/    ######/    ##/\n");
	gc_pBQBInterface->bprintf("                     ## |\n");
	gc_pBQBInterface->bprintf("                     ## |\n");
	gc_pBQBInterface->bprintf("                     ##/\n");
	gc_pBQBInterface->printf("\n");
	gc_pBQBInterface->printf("    LibBQB - File Format and Compression Algorithms for Trajectories of\n");
	gc_pBQBInterface->printf("             Volumetric Data and Atom Positions\n");
	gc_pBQBInterface->printf("\n");
	gc_pBQBInterface->bprintf("    https://brehm-research.de/bqb\n");
	gc_pBQBInterface->printf("\n");
	gc_pBQBInterface->printf("    Free software, licensed under GNU LGPL v3\n");
	gc_pBQBInterface->printf("\n");
	gc_pBQBInterface->printf("    Copyright (c) Martin Brehm and Martin Thomas,\n");
	gc_pBQBInterface->printf("                  Martin Luther University Halle-Wittenberg (Germany), 2016 - 2021.\n");
	gc_pBQBInterface->printf("\n");


	time(&ltime);
	today = localtime(&ltime);
	strcpy(buf,asctime(today));
	buf[strlen(buf)-1] = 0;
	if (gc_sHostName != NULL)
		gc_pBQBInterface->printf("    # Running on %s at %s",gc_sHostName,buf);
	else
		gc_pBQBInterface->printf("    # Running at %s",buf);

	#if !defined(_WIN32) && (defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__)))
		gc_pBQBInterface->printf(" (PID %d)\n",getpid());
	#elif defined(_WIN32)
		gc_pBQBInterface->printf(" (PID %d)\n",GetCurrentProcessId());
	#else
		gc_pBQBInterface->printf("\n");
	#endif

	if (gc_sWorkingDir != NULL)
		gc_pBQBInterface->printf("    # Running in %s\n",gc_sWorkingDir);

	gc_pBQBInterface->printf("    # Source code version: ");
	gc_pBQBInterface->printf("%s\n",BQB_SOURCE_VERSION);
	gc_pBQBInterface->printf("    # Compiled at ");
	gc_pBQBInterface->printf("%s, %s\n",__DATE__,__TIME__);

	#if defined(__VERSION) || defined(__GNUC__) || defined(_MSC_VER) || defined(__clang__) || defined(__llvm__)
		gc_pBQBInterface->printf("    # Compiler");
		b = false;
		#ifdef __VERSION__
			gc_pBQBInterface->printf(" \"%s\"",__VERSION__);
			b = true;
		#endif
		#ifdef __llvm__
			if (b)
				gc_pBQBInterface->printf(",");
			gc_pBQBInterface->printf(" LLVM");
			b = true;
		#endif
		#ifdef __clang__
			if (b)
				gc_pBQBInterface->printf(",");
			#if defined(__clang_major__) && defined(__clang_minor__) && defined(__clang_patchlevel__)
				gc_pBQBInterface->printf(" CLANG %d.%d.%d",__clang_major__,__clang_minor__,__clang_patchlevel__);
			#elif defined(__clang_version__)
				gc_pBQBInterface->printf(" CLANG \"%s\"",__clang_version__);
			#else
				gc_pBQBInterface->printf(" CLANG");
			#endif
			b = true;
		#endif
		#ifdef __GNUC__
			if (b)
				gc_pBQBInterface->printf(",");
			#ifdef __MINGW32__
				gc_pBQBInterface->printf(" MinGW GCC %d.%d.%d",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__);
			#else
				gc_pBQBInterface->printf(" GCC %d.%d.%d",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__);
			#endif
			b = true;
		#endif
		#ifdef _MSC_VER
			if (b)
				gc_pBQBInterface->printf(",");
			gc_pBQBInterface->printf(" MSVC++ ");
			// See  https://docs.microsoft.com/de-de/cpp/preprocessor/predefined-macros?view=vs-2019
			switch(_MSC_VER) {
				case 1928: gc_pBQBInterface->printf("16.8, VS 2019"); break;
				case 1927: gc_pBQBInterface->printf("16.7, VS 2019"); break;
				case 1926: gc_pBQBInterface->printf("16.6, VS 2019"); break;
				case 1925: gc_pBQBInterface->printf("16.5, VS 2019"); break;
				case 1924: gc_pBQBInterface->printf("16.4, VS 2019"); break;
				case 1923: gc_pBQBInterface->printf("16.3, VS 2019"); break;
				case 1922: gc_pBQBInterface->printf("16.2, VS 2019"); break;
				case 1921: gc_pBQBInterface->printf("16.1, VS 2019"); break;
				case 1920: gc_pBQBInterface->printf("16.0, VS 2019"); break;
				case 1916: gc_pBQBInterface->printf("15.9, VS 2017"); break;
				case 1915: gc_pBQBInterface->printf("15.8, VS 2017"); break;
				case 1914: gc_pBQBInterface->printf("15.7, VS 2017"); break;
				case 1913: gc_pBQBInterface->printf("15.6, VS 2017"); break;
				case 1912: gc_pBQBInterface->printf("15.5, VS 2017"); break;
				case 1911: gc_pBQBInterface->printf("15.3, VS 2017"); break;
				case 1910: gc_pBQBInterface->printf("15.1, VS 2017"); break;
				case 1900: gc_pBQBInterface->printf("14.0, VS 2015"); break;
				case 1800: gc_pBQBInterface->printf("12.0, VS 2013"); break;
				case 1700: gc_pBQBInterface->printf("11.0, VS 2012"); break;
				case 1600: gc_pBQBInterface->printf("10.0, VS 2010"); break;
				case 1500: gc_pBQBInterface->printf("9.0, VS 2008"); break;
				case 1400: gc_pBQBInterface->printf("8.0, VS 2005"); break;
				case 1310: gc_pBQBInterface->printf("7.1, VS 2003"); break;
				case 1300: gc_pBQBInterface->printf("7.0"); break;
				case 1200: gc_pBQBInterface->printf("6.0"); break;
				case 1100: gc_pBQBInterface->printf("5.0"); break;
				default: gc_pBQBInterface->printf("%d",_MSC_VER); break;
			}
			b = true;
		#endif
		gc_pBQBInterface->printf("\n");
	#endif

	gc_pBQBInterface->printf("\n");


	gc_pBQBInterface->bprintf("    Please cite:");
	gc_pBQBInterface->printf("  M. Brehm, M. Thomas: \"An Efficient Lossless Compression Algorithm\n");
	gc_pBQBInterface->printf("                  for Trajectories of Atom Positions and Volumetric Data\",\n");
	gc_pBQBInterface->printf("                  J. Chem. Inf. Model. 2018, 58 (10), pp 2092-2107.\n");
	gc_pBQBInterface->printf("\n");

	err = false;

	gc_pBQBInterface->printf("Command line:\n");
	gc_pBQBInterface->printf("  \"");
	for (z=1;z<argc;z++) {
		gc_pBQBInterface->printf("%s",argv[z]);
		if (z+1 < argc)
			gc_pBQBInterface->printf(" ");
	}
	gc_pBQBInterface->printf("\"\n\n");

	CC_InitGlobalVars();

	BQBCheckSourceVersion();

	if (argc == 2) {
		if (strcmp(argv[1],"-revision") == 0) {
			BQBWriteRevisionInfo();
			goto _end;
		}
	}

	if (!CC_ParseArguments(argc,argv)) {
		CC_CommandLineHelp();
		return 1;
	}

	driver = gc_pBQBInterface->CreateDriver( 0 );


	if (gc_bCompress) {

		if (gc_bFile) {

			gc_pBQBInterface->printf("Will compress file.\n\n");

			gc_pBQBInterface->bprintf("*** Parameters ***\n");
			gc_pBQBInterface->printf("    Input file:                %s\n",
				(gc_sInFile!=NULL)?gc_sInFile:"(NULL)");
			gc_pBQBInterface->printf("    Output file:               %s\n",
				(gc_sOutFile!=NULL)?gc_sOutFile:"(NULL)");
			gc_pBQBInterface->printf("    Verbose:                   -verbose %d\n",
				gc_iFileVerbose);
			gc_pBQBInterface->printf("    Dry run:                   -dryrun %s\n",
				gc_bDryRun?"yes":"no");
			gc_pBQBInterface->printf("%s",gc_pParmFile->ToString(4).c_str());
			gc_pBQBInterface->printf("\n");

			gc_pBQBInterface->bprintf("Parameter Key:");
			gc_pBQBInterface->printf(" %s\n\n",gc_pParmFile->ToKey().c_str());

			if (gc_bDryRun) {
				gc_pBQBInterface->printf("This is a dry run. Not processing file. Leaving now.\n\n");
				goto _end;
			}

			if (gc_sInFile == NULL) {
				gc_pBQBInterface->eprintf("Error: Input file name is required.\n");
				err = true;
				goto _end;
			}

			if (gc_sOutFile != NULL) {
				if (strlen(gc_sOutFile) > 4) {
					if (bqb_strcmp_nocase(&gc_sOutFile[strlen(gc_sOutFile)-4],".bqb") == 0)
						goto _okcfile;
					if (bqb_strcmp_nocase(&gc_sOutFile[strlen(gc_sOutFile)-4],".bbq") == 0) {
						BQB_PrintBMode();
						goto _okcfile;
					}
				}
				gc_pBQBInterface->eprintf("Error: Output file extension needs to be \".bqb\".\n");
				err = true;
				goto _end;
_okcfile:;
			}

			if (gc_iFileVerbose == 2)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_DEBUG);
			else if (gc_iFileVerbose == 1)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_VERBOSE);
			else
				gc_pBQBInterface->SetPrintLevel(BQB_PL_STANDARD);

			if (!driver->CompressFile(
				gc_sInFile,  // Input File
				gc_sOutFile,  // Output File
				gc_pParmFile
			)) {
				gc_pBQBInterface->eprintf("CompressFile returned an error.\n");
				err = true;
				goto _end;
			} else
				gc_pBQBInterface->printf("CompressFile successfully returned.\n");
		}

		if (gc_bCube) {

			gc_pBQBInterface->printf("Will compress cube file.\n\n");

			b = false;

			if (gc_iCubeOptimize != 0) {

				if (gc_iCubeOptSteps == -1) {
					gc_pBQBInterface->printf("    -optsteps not specified, defaulting to 20.\n");
					gc_iCubeOptSteps = 20;
					b = true;
				}

				if ((gc_iCubeOptSteps > gc_iCubeSteps) && (gc_iCubeSteps != -1)) {
					gc_pBQBInterface->printf("    -optsteps larger than -steps, reducing to %d.\n",gc_iCubeSteps);
					gc_iCubeOptSteps = gc_iCubeSteps;
					b = true;
				}
			}

			if (b)
				gc_pBQBInterface->printf("\n");

			gc_pBQBInterface->bprintf("*** Parameters ***\n");
			gc_pBQBInterface->printf("  General:\n");
			gc_pBQBInterface->printf("    Infile:                    %s\n",
				(gc_sInFile!=NULL)?gc_sInFile:"(NULL)");
			gc_pBQBInterface->printf("    Outfile:                   %s\n",
				(gc_sOutFile!=NULL)?gc_sOutFile:"(NULL)");
			gc_pBQBInterface->printf("    Output Ref. Traj.:         -refout %s\n",
				(gc_sCubeRefOut!=NULL)?gc_sCubeRefOut:"(NULL)");
			gc_pBQBInterface->printf("    Start:                     -start %d\n",
				gc_iCubeStart+1);
			gc_pBQBInterface->printf("    Steps:                     -steps %d\n",
				gc_iCubeSteps);
			gc_pBQBInterface->printf("    Stride:                    -stride %d\n",
				gc_iCubeStride);
			gc_pBQBInterface->printf("    Keyframe frequency:        -keyframe %d\n",
				gc_iCubeKey);
			gc_pBQBInterface->printf("    Keep Comment Line:         -comment %s\n",
				gc_bCubeKeepComment?"yes":"no");
			gc_pBQBInterface->printf("    Check Compressed Data:     -check %s\n",
				gc_bCubeCompare?"yes":"no");
			gc_pBQBInterface->printf("    Dummy Read (no comp.):     -dummyread %s\n",
				gc_bCubeDummyRead?"yes":"no");
			gc_pBQBInterface->printf("    Optimize Parameters:       -optimize %d (0=off, 1=quick, 2=patient, 3=precise, 4=exhaustive)\n",
				gc_iCubeOptimize);
			gc_pBQBInterface->printf("    Optimization Frames:       -optsteps %d\n",
				gc_iCubeOptSteps);
			gc_pBQBInterface->printf("    Quit after Optim.:         -onlyopt %s\n",
				gc_bCubeOnlyOpt?"yes":"no");
			gc_pBQBInterface->printf("    Opt. from Start on:        -optstart %s\n",
				(gc_iCubeOptIncludeFirst==1)?"yes":((gc_iCubeOptIncludeFirst==0)?"no":"(undef)"));
			gc_pBQBInterface->printf("    Verbose:                   -verbose %d\n",
				gc_iCubeVerbose);
			gc_pBQBInterface->printf("    Dry run:                   -dryrun %s\n",
				gc_bDryRun?"yes":"no");
			gc_pBQBInterface->printf("%s",gc_pParmVol->ToString(2).c_str());
			gc_pBQBInterface->printf("\n");

			gc_pBQBInterface->bprintf("Parameter Key:");
			gc_pBQBInterface->printf(" %s\n\n",gc_pParmVol->ToKey().c_str());

			if (gc_bDryRun) {
				gc_pBQBInterface->printf("This is a dry run. Not processing file. Leaving now.\n\n");
				goto _end;
			}

			if (gc_sInFile == NULL) {
				gc_pBQBInterface->eprintf("Error: Input file name is required.\n");
				err = true;
				goto _end;
			}

			if (strlen(gc_sInFile) < 6) {
				gc_pBQBInterface->eprintf("Error: Input file name is too short (needs to be *.cube).\n");
				err = true;
				goto _end;
			}

			if (bqb_strcmp_nocase( &gc_sInFile[strlen(gc_sInFile)-5], ".cube" ) != 0) {
				gc_pBQBInterface->eprintf("Error: Only Gaussian Cube trajectories are currently supported as input.\n");
				err = true;
				goto _end;
			}

			if (!gc_bOptTouched)
				gc_pBQBInterface->printf("To disable the automatic parameter optimization, add \"-optimize 0\" to the command line.\n\n");

			if (gc_iCubeOptimize == 3)
				gc_pBQBInterface->printf("To achieve highest compression ratio (at the cost of compression time), add \"-vsplit 14 -viter 40\".\n\n");

			if (gc_sOutFile != NULL) {
				if (strlen(gc_sOutFile) > 4) {
					if (bqb_strcmp_nocase(&gc_sOutFile[strlen(gc_sOutFile)-4],".bqb") == 0)
						goto _okcvol;
					if (bqb_strcmp_nocase(&gc_sOutFile[strlen(gc_sOutFile)-4],".bbq") == 0) {
						BQB_PrintBMode();
						goto _okcvol;
					}
				}
				gc_pBQBInterface->eprintf("Error: Output file extension needs to be \".bqb\".\n");
				err = true;
				goto _end;
_okcvol:;
			}

			driver->m_iOptIncludeFirst = gc_iCubeOptIncludeFirst;

			if (gc_iCubeVerbose == 2)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_DEBUG);
			else if (gc_iCubeVerbose == 1)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_VERBOSE);
			else
				gc_pBQBInterface->SetPrintLevel(BQB_PL_STANDARD);

			if (!driver->CompressCube(
				gc_sInFile,  // Cube Trajectory
				gc_sOutFile,
				gc_sCubeRefOut,     // Output reference file
				gc_iCubeStart,      // Start
				gc_iCubeSteps,       // Steps
				gc_iCubeStride,        // Stride
				gc_iCubeKey,
				gc_pParmVol,
				gc_iCubeOptimize,
				gc_iCubeOptSteps,
				gc_bCubeOnlyOpt,
				gc_bCubeKeepComment,    // Keep Comment Line
				gc_bCubeCompare,     // Compare decompressed results
				gc_bCubeDummyRead   // Dummy Read (no compression or output)
			)) {
				gc_pBQBInterface->eprintf("CompressCube returned an error.\n");
				err = true;
				goto _end;
			} else
				gc_pBQBInterface->printf("CompressCube successfully returned.\n");
		}

		if (gc_bXYZ) {

			gc_pBQBInterface->printf("Will compress XYZ file.\n\n");

			b = false;

			if (gc_iAtomOptimize != 0) {

				if (gc_iAtomOptSteps == -1) {
					gc_pBQBInterface->printf("    -optsteps not specified, defaulting to 50.\n");
					gc_iAtomOptSteps = 50;
					b = true;
				}

				if ((gc_iAtomOptSteps > gc_iAtomSteps) && (gc_iAtomSteps != -1)) {
					gc_pBQBInterface->printf("    -optsteps larger than -steps, reducing to %d.\n",gc_iAtomSteps);
					gc_iAtomOptSteps = gc_iAtomSteps;
					b = true;
				}
			}

			if (b)
				gc_pBQBInterface->printf("\n");

			gc_pBQBInterface->bprintf("*** Parameters ***\n");
			gc_pBQBInterface->printf("  General:\n");
			gc_pBQBInterface->printf("    Infile:                    %s\n",
				(gc_sInFile!=NULL)?gc_sInFile:"(NULL)");
			gc_pBQBInterface->printf("    Outfile:                   %s\n",
				(gc_sOutFile!=NULL)?gc_sOutFile:"(NULL)");
			gc_pBQBInterface->printf("    Output Ref. Traj.:         -refout %s\n",
				(gc_sAtomRefOut!=NULL)?gc_sAtomRefOut:"(NULL)");
			gc_pBQBInterface->printf("    Start:                     -start %d\n",
				gc_iAtomStart);
			gc_pBQBInterface->printf("    Steps:                     -steps %d\n",
				gc_iAtomSteps);
			gc_pBQBInterface->printf("    Stride:                    -stride %d\n",
				gc_iAtomStride);
			gc_pBQBInterface->printf("    Keyframe frequency:        -keyframe %d\n",
				gc_iAtomKey);
			gc_pBQBInterface->printf("    Keep Comment Line:         -comment %s\n",
				gc_bAtomKeepComment?"yes":"no");
			gc_pBQBInterface->printf("    Check Compressed Data:     -check %s\n",
				gc_bAtomCompare?"yes":"no");
			gc_pBQBInterface->printf("    Optimize Parameters:       -optimize %d (0=off, 1=quick, 2=patient, 3=precise, 4=exhaustive)\n",
				gc_iAtomOptimize);
			gc_pBQBInterface->printf("    Optimization Frames:       -optsteps %d\n",
				gc_iAtomOptSteps);
			gc_pBQBInterface->printf("    Quit after Optim.:         -onlyopt %s\n",
				gc_bAtomOnlyOpt?"yes":"no");
			gc_pBQBInterface->printf("    Opt. from Start on:        -optstart %s\n",
				(gc_iAtomOptIncludeFirst==1)?"yes":((gc_iAtomOptIncludeFirst==0)?"no":"(undef)"));
			gc_pBQBInterface->printf("    Verbose:                   -verbose %d\n",
				gc_iAtomVerbose);
			gc_pBQBInterface->printf("    Dry run:                   -dryrun %s\n",
				gc_bDryRun?"yes":"no");
			gc_pBQBInterface->printf("%s",gc_pParmPos->ToString(2).c_str());
			gc_pBQBInterface->printf("\n");

			gc_pBQBInterface->bprintf("Parameter Key:");
			gc_pBQBInterface->printf(" %s\n\n",gc_pParmPos->ToKey().c_str());

			if (gc_bDryRun) {
				gc_pBQBInterface->printf("This is a dry run. Not processing file. Leaving now.\n\n");
				goto _end;
			}

			if (gc_sInFile == NULL) {
				gc_pBQBInterface->eprintf("Error: Input file name is required.\n");
				err = true;
				goto _end;
			}

			if (strlen(gc_sInFile) < 5) {
				gc_pBQBInterface->eprintf("Error: Input file name is too short (needs to be *.xyz).\n");
				gc_pBQBInterface->printf("       More formats will be added in future releases.\n");
				err = true;
				goto _end;
			}

			if (bqb_strcmp_nocase( &gc_sInFile[strlen(gc_sInFile)-4], ".xyz" ) != 0) {
				gc_pBQBInterface->eprintf("Error: Only XMol XYZ trajectories are currently supported as input.\n");
				gc_pBQBInterface->printf("       More formats will be added in future releases.\n");
				err = true;
				goto _end;
			}

			if (!gc_bOptTouched)
				gc_pBQBInterface->printf("To disable the automatic parameter optimization, add \"-optimize 0\" to the command line.\n\n");

			if (gc_sOutFile != NULL) {
				if (strlen(gc_sOutFile) > 4) {
					if (bqb_strcmp_nocase(&gc_sOutFile[strlen(gc_sOutFile)-4],".bqb") == 0)
						goto _okcpos;
					if (bqb_strcmp_nocase(&gc_sOutFile[strlen(gc_sOutFile)-4],".btr") == 0)
						goto _okcpos;
					if (bqb_strcmp_nocase(&gc_sOutFile[strlen(gc_sOutFile)-4],".bbq") == 0) {
						BQB_PrintBMode();
						goto _okcpos;
					}
				}
				gc_pBQBInterface->eprintf("Error: Output file extension needs to be \".bqb\" or \".btr\".\n");
				err = true;
				goto _end;
_okcpos:;
			}

			driver->m_iOptIncludeFirst = gc_iAtomOptIncludeFirst;

			if (gc_iAtomVerbose == 2)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_DEBUG);
			else if (gc_iAtomVerbose == 1)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_VERBOSE);
			else
				gc_pBQBInterface->SetPrintLevel(BQB_PL_STANDARD);

			if (!driver->CompressXYZ(
				gc_sInFile,  // XYZ file name
				gc_sOutFile,  //
				gc_sAtomRefOut,  // Output reference trajectory
				gc_iAtomStart,     // Start
				gc_iAtomSteps,     // Steps
				gc_iAtomStride,        // Stride
				gc_iAtomKey,
				gc_pParmPos,
				gc_bAtomKeepComment,    // Keep Comment Line
				gc_bAtomCompare,     // Compare decompressed results
				gc_iAtomOptimize,
				gc_iAtomOptSteps,
				gc_bAtomOnlyOpt
			)) {
				gc_pBQBInterface->eprintf("CompressXYZ returned an error.\n");
				err = true;
				goto _end;
			} else
				gc_pBQBInterface->printf("CompressXYZ successfully returned.\n");
		}
	}


	if (gc_bDecompress) {

		if (gc_bFile) {

			gc_pBQBInterface->printf("Will decompress file.\n\n");

			gc_pBQBInterface->bprintf("*** Parameters ***\n");
			gc_pBQBInterface->printf("    Infile:                %s\n",
				(gc_sInFile!=NULL)?gc_sInFile:"(NULL)");
			gc_pBQBInterface->printf("    Outfile:               %s\n",
				(gc_sOutFile!=NULL)?gc_sOutFile:"(NULL)");
			gc_pBQBInterface->printf("    Reference File:        -refcompare %s\n",
				(gc_sRefFile!=NULL)?gc_sRefFile:"(NULL)");
			gc_pBQBInterface->printf("    Verbose:               -verbose %d\n",
				gc_iFileVerbose);
			gc_pBQBInterface->printf("\n");

			if (gc_sInFile == NULL) {
				gc_pBQBInterface->eprintf("Error: Input file name is required.\n");
				err = true;
				goto _end;
			}

			if (strlen(gc_sInFile) > 4) {
				if (bqb_strcmp_nocase(&gc_sInFile[strlen(gc_sInFile)-4],".bqb") == 0)
					goto _okdfile;
				if (bqb_strcmp_nocase(&gc_sInFile[strlen(gc_sInFile)-4],".bbq") == 0) {
					BQB_PrintBMode();
					goto _okdfile;
				}
			}
			gc_pBQBInterface->eprintf("Error: Input file extension needs to be \".bqb\".\n");
			err = true;
			goto _end;
_okdfile:

			if (gc_iFileVerbose == 2)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_DEBUG);
			else if (gc_iFileVerbose == 1)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_VERBOSE);
			else
				gc_pBQBInterface->SetPrintLevel(BQB_PL_STANDARD);

			if (!driver->DecompressFile(
				gc_sInFile,  // Input File
				gc_sOutFile,  // Output File
				gc_sRefFile // Ref File
			)) {
				gc_pBQBInterface->eprintf("DecompressFile returned an error.\n");
				err = true;
				goto _end;
			} else
				gc_pBQBInterface->printf("DecompressFile successfully returned.\n");
		}

		if (gc_bCube) {

			gc_pBQBInterface->printf("Will decompress Cube file.\n\n");

			gc_pBQBInterface->bprintf("*** Parameters ***\n");
			gc_pBQBInterface->printf("    Infile:                %s\n",
				(gc_sInFile!=NULL)?gc_sInFile:"(NULL)");
			gc_pBQBInterface->printf("    Outfile:               %s\n",
				(gc_sOutFile!=NULL)?gc_sOutFile:"(NULL)");
			gc_pBQBInterface->printf("    Reference File:        -refcompare %s\n",
				(gc_sRefFile!=NULL)?gc_sRefFile:"(NULL)");
			gc_pBQBInterface->printf("    Steps:                 -steps %d\n",
				gc_iCubeSteps);
			gc_pBQBInterface->printf("    Stride:                -stride %d\n",
				gc_iCubeStride);
			gc_pBQBInterface->printf("    Verbose:               -verbose %d\n",
				gc_iCubeVerbose);
			gc_pBQBInterface->printf("\n");

			if (gc_sInFile == NULL) {
				gc_pBQBInterface->eprintf("Error: Input file name is required.\n");
				err = true;
				goto _end;
			}

			if (gc_sOutFile != NULL) {

				if (strlen(gc_sOutFile) < 6) {
					gc_pBQBInterface->eprintf("Error: Output file name is too short (needs to be *.cube).\n");
					err = true;
					goto _end;
				}

				if (bqb_strcmp_nocase( &gc_sOutFile[strlen(gc_sOutFile)-5], ".cube" ) != 0) {
					gc_pBQBInterface->eprintf("Error: Only Gaussian Cube trajectories are currently supported as output.\n");
					err = true;
					goto _end;
				}
			}

			if (strlen(gc_sInFile) > 4) {
				if (bqb_strcmp_nocase(&gc_sInFile[strlen(gc_sInFile)-4],".bqb") == 0)
					goto _okdcube;
				if (strlen(gc_sInFile) > 6)
					if (bqb_strcmp_nocase(&gc_sInFile[strlen(gc_sInFile)-6],".blist") == 0)
						goto _okdcube;
				if (bqb_strcmp_nocase(&gc_sInFile[strlen(gc_sInFile)-4],".bbq") == 0) {
					BQB_PrintBMode();
					goto _okdcube;
				}
			}
			gc_pBQBInterface->eprintf("Error: Input file extension needs to be \".bqb\" or \".blist\".\n");
			err = true;
			goto _end;
_okdcube:

			if (gc_iCubeVerbose == 2)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_DEBUG);
			else if (gc_iCubeVerbose == 1)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_VERBOSE);
			else
				gc_pBQBInterface->SetPrintLevel(BQB_PL_STANDARD);

			if (!driver->DecompressCube(
				gc_sInFile,
				gc_sOutFile,
				gc_sRefFile,           // Input reference file
				gc_iCubeSteps,             // Steps
				gc_iCubeStride,              // Stride
				NULL
			)) {
				gc_pBQBInterface->eprintf("DecompressCube returned an error.\n");
				err = true;
				goto _end;
			} else
				gc_pBQBInterface->printf("DecompressCube successfully returned.\n");
		}

		if (gc_bXYZ) {

			gc_pBQBInterface->printf("Will decompress XYZ file.\n\n");

			gc_pBQBInterface->bprintf("*** Parameters ***\n");
			gc_pBQBInterface->printf("    Infile:                %s\n",
				(gc_sInFile!=NULL)?gc_sInFile:"(NULL)");
			gc_pBQBInterface->printf("    Outfile:               %s\n",
				(gc_sOutFile!=NULL)?gc_sOutFile:"(NULL)");
			gc_pBQBInterface->printf("    Reference File:        -refcompare %s\n",
				(gc_sRefFile!=NULL)?gc_sRefFile:"(NULL)");
			gc_pBQBInterface->printf("    Steps:                 -steps %d\n",
				gc_iAtomSteps);
			gc_pBQBInterface->printf("    Stride:                -stride %d\n",
				gc_iAtomStride);
			gc_pBQBInterface->printf("    Verbose:               -verbose %d\n",
				gc_iAtomVerbose);
			gc_pBQBInterface->printf("\n");

			if (gc_sInFile == NULL) {
				gc_pBQBInterface->eprintf("Error: Input file name is required.\n");
				err = true;
				goto _end;
			}

			if (gc_sOutFile != NULL) {

				if (strlen(gc_sOutFile) < 5) {
					gc_pBQBInterface->eprintf("Error: Output file name is too short (needs to be *.xyz).\n");
					gc_pBQBInterface->printf("       More formats will be added in future releases.\n");
					err = true;
					goto _end;
				}

				if (bqb_strcmp_nocase( &gc_sOutFile[strlen(gc_sOutFile)-4], ".xyz" ) != 0) {
					gc_pBQBInterface->eprintf("Error: Only XMol XYZ trajectories are currently supported as output.\n");
					gc_pBQBInterface->printf("       More formats will be added in future releases.\n");
					err = true;
					goto _end;
				}
			}

			if (strlen(gc_sInFile) > 4) {
				if (bqb_strcmp_nocase(&gc_sInFile[strlen(gc_sInFile)-4],".bqb") == 0)
					goto _okdpos;
				if (bqb_strcmp_nocase(&gc_sInFile[strlen(gc_sInFile)-4],".btr") == 0)
					goto _okdpos;
				if (strlen(gc_sInFile) > 6)
					if (bqb_strcmp_nocase(&gc_sInFile[strlen(gc_sInFile)-6],".blist") == 0)
						goto _okdpos;
				if (bqb_strcmp_nocase(&gc_sInFile[strlen(gc_sInFile)-4],".bbq") == 0) {
					BQB_PrintBMode();
					goto _okdpos;
				}
			}
			gc_pBQBInterface->eprintf("Error: Input file extension needs to be \".bqb\", \".btr\", or \".blist\".\n");
			err = true;
			goto _end;
_okdpos:

			if (gc_iAtomVerbose == 2)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_DEBUG);
			else if (gc_iAtomVerbose == 1)
				gc_pBQBInterface->SetPrintLevel(BQB_PL_VERBOSE);
			else
				gc_pBQBInterface->SetPrintLevel(BQB_PL_STANDARD);

			if (!driver->DecompressXYZ(
				gc_sInFile,  //
				gc_sOutFile,  // XYZ file name
				gc_sRefFile,  //
				gc_iAtomSteps,       // Steps
				gc_iAtomStride,        // Stride
				NULL
			)) {
				gc_pBQBInterface->eprintf("DecompressXYZ returned an error.\n");
				err = true;
				goto _end;
			} else
				gc_pBQBInterface->printf("DecompressXYZ successfully returned.\n");
		}
	}


	if (gc_bCheck) {

	/*	gc_pBQBInterface->eprintf("\"check\": Not yet implemented.\n");
		err = true;
		goto _end;*/


/*		gc_pBQBInterface->printf("\n###### Seek Test ######\n\n");

		FILE *a;

		gc_pBQBInterface->printf("\n");
		gc_pBQBInterface->printf("##########################\n");
		gc_pBQBInterface->printf("###   1.) Seek to 17   ###\n");
		gc_pBQBInterface->printf("##########################\n");
		bf.OpenRead(gc_sInFile);
		SetBQBVerbose(true);
		bf.SeekFrame(17);
		bf.ReadFrame();
		SetBQBVerbose(false);
		a = OpenFileWrite("seek_17.bqb",false);
		WriteArrayToFile(a,*bf.GetFramePayload());
		fclose(a);
		bf.Close();

		gc_pBQBInterface->printf("\n");
		gc_pBQBInterface->printf("##########################\n");
		gc_pBQBInterface->printf("###   2.) Skip to 17   ###\n");
		gc_pBQBInterface->printf("##########################\n");
		bf.OpenRead(gc_sInFile);
		bf.SkipFrames(17);
		bf.ReadFrame();
		a = OpenFileWrite("skip_17.bqb",false);
		WriteArrayToFile(a,*bf.GetFramePayload());
		fclose(a);
		bf.Close();

		gc_pBQBInterface->printf("\n");
		gc_pBQBInterface->printf("#############################\n");
		gc_pBQBInterface->printf("###   3.) Seek to 19999   ###\n");
		gc_pBQBInterface->printf("#############################\n");
		bf.OpenRead(gc_sInFile);
		SetBQBVerbose(true);
		bf.SeekFrame(19999);
		bf.ReadFrame();
		SetBQBVerbose(false);
		a = OpenFileWrite("seek_19999.bqb",false);
		WriteArrayToFile(a,*bf.GetFramePayload());
		fclose(a);
		bf.Close();

		gc_pBQBInterface->printf("\n");
		gc_pBQBInterface->printf("#############################\n");
		gc_pBQBInterface->printf("###   4.) Skip to 19999   ###\n");
		gc_pBQBInterface->printf("#############################\n");
		bf.OpenRead(gc_sInFile);
		bf.SkipFrames(19999);
		bf.ReadFrame();
		a = OpenFileWrite("skip_19999.bqb",false);
		WriteArrayToFile(a,*bf.GetFramePayload());
		fclose(a);
		bf.Close();

		gc_pBQBInterface->printf("\n");
		gc_pBQBInterface->printf("#############################\n");
		gc_pBQBInterface->printf("###   5.) Seek to 20000   ###\n");
		gc_pBQBInterface->printf("#############################\n");
		bf.OpenRead(gc_sInFile);
		SetBQBVerbose(true);
		bf.SeekFrame(20000);
		bf.ReadFrame();
		SetBQBVerbose(false);
		a = OpenFileWrite("seek_20000.bqb",false);
		WriteArrayToFile(a,*bf.GetFramePayload());
		fclose(a);
		bf.Close();

		gc_pBQBInterface->printf("\n");
		gc_pBQBInterface->printf("#############################\n");
		gc_pBQBInterface->printf("###   6.) Skip to 20000   ###\n");
		gc_pBQBInterface->printf("#############################\n");
		bf.OpenRead(gc_sInFile);
		bf.SkipFrames(20000);
		bf.ReadFrame();
		a = OpenFileWrite("skip_20000.bqb",false);
		WriteArrayToFile(a,*bf.GetFramePayload());
		fclose(a);
		bf.Close();

		gc_pBQBInterface->printf("\n");
		gc_pBQBInterface->printf("#############################\n");
		gc_pBQBInterface->printf("###   7.) Seek to 20017   ###\n");
		gc_pBQBInterface->printf("#############################\n");
		bf.OpenRead(gc_sInFile);
		SetBQBVerbose(true);
		bf.SeekFrame(20017);
		bf.ReadFrame();
		SetBQBVerbose(false);
		a = OpenFileWrite("seek_20017.bqb",false);
		WriteArrayToFile(a,*bf.GetFramePayload());
		fclose(a);
		bf.Close();

		gc_pBQBInterface->printf("\n");
		gc_pBQBInterface->printf("#############################\n");
		gc_pBQBInterface->printf("###   8.) Skip to 20017   ###\n");
		gc_pBQBInterface->printf("#############################\n");
		bf.OpenRead(gc_sInFile);
		bf.SkipFrames(20017);
		bf.ReadFrame();
		a = OpenFileWrite("skip_20017.bqb",false);
		WriteArrayToFile(a,*bf.GetFramePayload());
		fclose(a);
		bf.Close();

		gc_pBQBInterface->printf("\n###### Seek Test done ######\n\n");

		return 0;*/

		gc_pBQBInterface->printf("Will check BQB file.\n\n");

		gc_pBQBInterface->bprintf("*** Parameters ***\n");
		gc_pBQBInterface->printf("    Infile:                %s\n",(gc_sInFile!=NULL)?gc_sInFile:"(NULL)");
		gc_pBQBInterface->printf("    Verbose:               -verbose %d\n",gc_iCheckVerbose);
		gc_pBQBInterface->printf("\n");

		if (gc_sInFile == NULL) {
			gc_pBQBInterface->eprintf("Error: Input file name is required.\n");
			err = true;
			goto _end;
		}

		gc_pBQBInterface->printf("Checking integrity of compressed file \"%s\".\n\n",gc_sInFile);

		if (gc_iCheckVerbose == 2)
			gc_pBQBInterface->SetPrintLevel(BQB_PL_DEBUG);
		else if (gc_iCheckVerbose == 1)
			gc_pBQBInterface->SetPrintLevel(BQB_PL_VERBOSE);
		else
			gc_pBQBInterface->SetPrintLevel(BQB_PL_STANDARD);

		CBQBFile bf( *gc_pBQBInterface );

		if (!bf.CheckIntegrity(gc_sInFile,gc_iCheckVerbose!=0)) {
			gc_pBQBInterface->eprintf("Integrity check failed.\n");
			err = true;
		}
	}


	if (gc_bCompare) {

		gc_pBQBInterface->eprintf("\"compare\": Not yet implemented.\n");
		goto _end;

		gc_pBQBInterface->printf("Will compare BQB coordinates to XYZ file.\n\n");

		gc_pBQBInterface->bprintf("*** Parameters ***\n");
		gc_pBQBInterface->printf("    Infile:                %s\n",(gc_sInFile!=NULL)?gc_sInFile:"(NULL)");
		gc_pBQBInterface->printf("    Reffile:               %s\n",(gc_sRefFile!=NULL)?gc_sRefFile:"(NULL)");
		gc_pBQBInterface->printf("    Verbose:               -verbose %d\n",gc_iCompareVerbose);
		gc_pBQBInterface->printf("\n");

		if (gc_sInFile == NULL) {
			gc_pBQBInterface->eprintf("Error: Input file name is required.\n");
			err = true;
			goto _end;
		}
		if (gc_sRefFile == NULL) {
			gc_pBQBInterface->eprintf("Error: Reference file name is required.\n");
			err = true;
			goto _end;
		}

		gc_pBQBInterface->printf("Comparing atom coordinates between %s and %s.\n\n",gc_sInFile,gc_sRefFile);

		if (gc_iCompareVerbose == 2)
			gc_pBQBInterface->SetPrintLevel(BQB_PL_DEBUG);
		else if (gc_iCompareVerbose == 1)
			gc_pBQBInterface->SetPrintLevel(BQB_PL_VERBOSE);
		else
			gc_pBQBInterface->SetPrintLevel(BQB_PL_STANDARD);

/*		if (!bf.CompareCoords(gc_sInFile,gc_sRefFile,gc_bCompareVerbose)) {
			gc_pBQBInterface->eprintf("Comparison failed.\n");
			err = true;
		}*/
	}


	if (gc_bMerge) {

		gc_pBQBInterface->eprintf("\"merge\": Not yet implemented.\n");
		goto _end;

		gc_pBQBInterface->printf("Will merge BQB files.\n\n");

		gc_pBQBInterface->bprintf("*** Parameters ***\n");
		gc_pBQBInterface->printf("    Outfile:               %s\n",(gc_sOutFile!=NULL)?gc_sOutFile:"(NULL)");
		if (gc_sInFileList.size() != 0) {
			gc_pBQBInterface->printf("    Infile:                %s\n",(gc_sInFileList[0]!=NULL)?gc_sInFileList[0]:"(NULL)");
			for (z=1;z<(int)gc_sInFileList.size();z++)
				gc_pBQBInterface->printf("                           %s\n",(gc_sInFileList[z]!=NULL)?gc_sInFileList[z]:"(NULL)");
		} else
			gc_pBQBInterface->printf("    Infile:                (NULL)\n");
		gc_pBQBInterface->printf("    Verbose:               -verbose %d\n",gc_iMergeVerbose);
		gc_pBQBInterface->printf("\n");

		if (gc_sOutFile == NULL) {
			gc_pBQBInterface->eprintf("Error: Output file name is required.\n");
			err = true;
			goto _end;
		}
		if (gc_sInFileList.size() == 0) {
			gc_pBQBInterface->eprintf("Error: At least one input file name is required.\n");
			err = true;
			goto _end;
		}
		gc_pBQBInterface->printf("Merging %lu BQB files into \"%s\".\n\n",(unsigned long)gc_sInFileList.size(),gc_sOutFile);

		if (gc_iMergeVerbose == 2)
			gc_pBQBInterface->SetPrintLevel(BQB_PL_DEBUG);
		else if (gc_iMergeVerbose == 1)
			gc_pBQBInterface->SetPrintLevel(BQB_PL_VERBOSE);
		else
			gc_pBQBInterface->SetPrintLevel(BQB_PL_STANDARD);

		if (!driver->MergeBQB(
				gc_sOutFile,
				gc_sInFileList
				//gc_bMergeVerbose
			)) {
			gc_pBQBInterface->eprintf("Merging failed.\n");
			err = true;
			goto _end;
		}
	}

	if (gc_bSplit) {

		gc_pBQBInterface->eprintf("\"split\": Not yet implemented.\n");
		goto _end;

		gc_pBQBInterface->printf("Will split BQB file.\n\n");

		gc_pBQBInterface->bprintf("*** Parameters ***\n");
		gc_pBQBInterface->printf("    Infile:                %s\n",(gc_sInFile!=NULL)?gc_sInFile:"(NULL)");
		gc_pBQBInterface->printf("    Outfile:               %s\n",(gc_sOutFile!=NULL)?gc_sOutFile:"(NULL)");
		gc_pBQBInterface->printf("    Split length:          %d\n",gc_iSplitLength);
		gc_pBQBInterface->printf("    Steps:                 -steps %d\n",gc_iSplitSteps);
		gc_pBQBInterface->printf("    Verbose:               -verbose %d\n",gc_iSplitVerbose);
		gc_pBQBInterface->printf("\n");

		if (gc_sInFile == NULL) {
			gc_pBQBInterface->eprintf("Error: Input file name is required.\n");
			err = true;
			goto _end;
		}
		if (gc_sOutFile == NULL) {
			gc_pBQBInterface->eprintf("Error: Output base name is required.\n");
			err = true;
			goto _end;
		}
		if (gc_iSplitLength == -1) {
			gc_pBQBInterface->eprintf("Error: Split length is required.\n");
			err = true;
			goto _end;
		}
		gc_pBQBInterface->printf("Splitting BQB file \"%s\" into pieces of %d frames with base name \"%s\".\n\n",gc_sInFile,gc_iSplitLength,gc_sOutFile);

		if (gc_iSplitVerbose == 2)
			gc_pBQBInterface->SetPrintLevel(BQB_PL_DEBUG);
		else if (gc_iSplitVerbose == 1)
			gc_pBQBInterface->SetPrintLevel(BQB_PL_VERBOSE);
		else
			gc_pBQBInterface->SetPrintLevel(BQB_PL_STANDARD);

/*		if (!engine->SplitBQB(gc_sInFile,gc_sOutFile,gc_iSplitLength,gc_iSplitSteps,gc_bSplitVerbose)) {
			gc_pBQBInterface->eprintf("Splitting failed.\n");
			err = true;
		}*/
	}

_end:

	gc_pBQBInterface->DestroyDriver( driver );

	gc_pBQBInterface->printf("\n");
	gc_pBQBInterface->bprintf("    Please cite:");
	gc_pBQBInterface->printf("  M. Brehm, M. Thomas: \"An Efficient Lossless Compression Algorithm\n");
	gc_pBQBInterface->printf("                  for Trajectories of Atom Positions and Volumetric Data\",\n");
	gc_pBQBInterface->printf("                  J. Chem. Inf. Model. 2018, 58 (10), pp 2092-2107.\n");

	gc_pBQBInterface->bprintf("\n    *** bqbtool leaving ***\n\n");

	// After this command, no more printing should occur
	BQBReleaseInterface( gc_pBQBInterface );

	#ifndef BQB_INSIDE_TRAVIS
		if (gc_fLogFile != NULL)
			fclose(gc_fLogFile);
	#endif


	if (err)
		return 1;
	else
		return 0;
}



