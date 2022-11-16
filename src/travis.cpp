/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2022 Martin Brehm
                  2012-2022 Martin Thomas
                  2016-2022 Sascha Gehrke

    Please cite:  J. Chem. Phys. 2020, 152 (16), 164105.         (DOI 10.1063/5.0005078 )
                  J. Chem. Inf. Model. 2011, 51 (8), 2007-2023.  (DOI 10.1021/ci200217w )

    This file was written by Martin Brehm and Martin Thomas.

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


// This must always be the first include directive
#include "config.h"

#include "travis.h"
#include "tools.h"
#include "database.h"
#include "statistics.h"
#include "maintools.h"
#include <float.h>
#include "dacf.h"
#include "lmwrapper.h"
#include "conversion.h"
#include "bqb_config.h"
#include "cluster.h"
#include "xdvector3.h"
#include "xdmatrix3.h"
#include "xquaternion.h"
#include "interface.h"
#include "random.h"
#include "plproj.h"
#include "sdfmap.h"
#include "linalg.h"
#include "fixplproj.h"
#include "matrixplot.h"
#include "sankey.h"

#ifdef TARGET_WINDOWS
	#include <windows.h>
	#include <direct.h>
#endif

#ifdef TARGET_LINUX
	#include <unistd.h>
#endif

#include <complex>
#include "largeinteger.h"

#ifdef USE_NETCDF
#include <netcdf.h>
#endif


const char *GetRevisionInfo_travis(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_travis() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


/***************************************************************************************************************
 ********************** Beginn Quellcode ***********************************************************************
 ***************************************************************************************************************/


int main(int argc, const char *argv[])
{
	FILE *a, *tfi;
	CxString buf, buf2, buf3, procsplitbuf;
	char *p, *q;
	const char *cp, *pn, *upn;
	int ti=-1, ti2=-1, ti3, cc, cc2, ticomb, tic_r=-1, tic_o=-1;
	int z0, z, z2, z3, z4, z5, z6, z7, z8, zr, zs, zi, z2b, z3b;
	int z1a, z1t, z2a, z2t;
	int procsplitc, procspliti;
	CxIntArray *a1, *a2;
	CAtomGroup *g1, *g2;
	int tia[3];
	CxDVector3 vec0, vec1, vec2, vec3, vec4, vec5, vecv, vecc;
	CxDMatrix3 mat;
	CxDMatrix3 dmat;
	double tf, tf2, tf3, tf4, tfs, tfx, tfy, tfz;
	double c0, c1, r;
	double *pd;
	double *pf;
	CMolecule *m, *m2;
	CSingleMolecule *sm, *smfix, *sm2;
	int ix, iy;
	char tc;
	unsigned long t0, t1, eta;
	int showinterval;
	CObservation *o;
	CTimeStep rot_ts;
	C3DF<double> *tempSDF;
	C3DF<double> *temp3DF;
	CDF *tdf;
	int sic;
	CxIntArray templa, *tla;
	CxDoubleArray **apfa; // "A"rray von "P"ointern auf "F"loat-"A"rrays
	double *tda;          // "T"emporaeres "D"ouble-"Array"
	CxByteArray **apba;
	CxDVec3Array tempvel;
	CAtomGroup *atgr, *ag;
	struct tm *today;
	time_t ltime;
	bool secondmolrun, tb, tbs;
	int multicounter;
	char multibuf[64];
	long fpos=0, fpos2=0;
	CMolBondGroup *bg;
	CMolBond *bond;
	CxDoubleArray tempfa, tempfa2, *ptfa, *ptfab;
	CxDoubleArray *ptfa2, *ptfa3, *ptfa2b;
	CGrace *gc, *gc2;
	CFFT *fft;
	C2DF *temp2df, *tempc2df, **temp2dfa;
	CConditionGroup *cg;
	int *tpi;
	CAutoCorrelation *ac;
	CCrossCorrelation *ccr;
	CDACFSub *dacfsub;
	CxDVector3 vecsnap1, vecsnap2, vecsnap3;
	CxDMatrix3 matsnap;
	CxString bufcom;
	std::vector<int> tiasnap;
	unsigned int *pSwapMatrix;
	CxDMatrix3 tdmat, tdmat2;
	CxDVector3 dvec0, dvec1, dvec2, dvec3, dvec4;
	CxQuaternion tq, tq2;
	CReorDyn *trdyn;
	char obuf[64];

	#ifdef TARGET_WINDOWS
		unsigned long tul;
	#endif






	InitGlobalVars();

	/**** Initialize Local Variables ***/
	procsplitc = 0;
	procspliti = 0;
	tdmat2.Unity();
	tq2.Unity();
	multicounter = 0;
	pSwapMatrix = NULL;
	apfa = NULL;
	tda = NULL;
	/***********************************/


	AddElementData();

	SortElementsLabel();

	try { g_sExeName = new char[strlen(argv[0])+1]; } catch(...) { g_sExeName = NULL; }
	if (g_sExeName == NULL) NewException((double)(strlen(argv[0])+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	strcpy(g_sExeName,argv[0]);

	GetTravisPath();

	InstallSignalHandler();

	InitColor();

	g_pLogFile = OpenFileWrite("travis.log",true);

	fprintf(g_pLogFile,"Command line:\n\"");
	for (z=0;z<argc;z++)
	{
		fprintf(g_pLogFile,"%s",argv[z]);
		if (z < argc-1)
			fprintf(g_pLogFile," ");
	}
	fprintf(g_pLogFile,"\"\n\n");
	fflush(g_pLogFile);

	showinterval = 1;

	srand((unsigned int)time(NULL));

	g_pBQBInterface = BQBCreateInterface(0);

	g_pBQBInterface->SetPrintCallback(  &mprintf );
	g_pBQBInterface->SetBPrintCallback( &bprintf );
	g_pBQBInterface->SetEPrintCallback( &eprintf );

#ifdef TARGET_WINDOWS
	tul = 256;
	g_sHomeDir = getenv("APPDATA");
	buf.SetBufSize(256);
	if (GetComputerNameA(buf.GetWritePointer(),&tul))
	{ 
		try { g_sHostName = new char[strlen(buf)+1]; } catch(...) { g_sHostName = NULL; }
		if (g_sHostName == NULL) NewException((double)(strlen(buf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		strcpy(g_sHostName,buf);
	} else g_sHostName = NULL;

	g_sWorkingDir = _getcwd(NULL,1024);
#elif defined(TARGET_LINUX)
	g_sHomeDir = getenv("HOME");
//	buf[0] = 0;
	buf.SetBufSize(256);
	if (gethostname(buf.GetWritePointer(),256)==0)
	{
		try { g_sHostName = new char[strlen(buf)+1]; } catch(...) { g_sHostName = NULL; }
		if (g_sHostName == NULL) NewException((double)(strlen(buf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		strcpy(g_sHostName,buf);
	} else g_sHostName = NULL;

	g_sWorkingDir = getcwd(NULL,1024);
#else
	g_sHomeDir = NULL;
	g_sHostName = NULL;
	g_sWorkingDir = NULL;
#endif

	ParsePassiveArgs(argc,argv);

	InitAnalyses();
	
	WriteHeader();

	InitDatabase();

	if (argc == 2) {
		if (strcmp(argv[1],"-revision") == 0) {
			WriteRevisionInfo();
			goto _ende;
		}
	}

	if (argc > 1) {

		// Call the bqbtool from TRAVIS
		if (
			(strcmp(argv[1],"compress") == 0) || 
			(strcmp(argv[1],"decompress") == 0) || 
			(strcmp(argv[1],"check") == 0) 
		//	(strcmp(argv[1],"compare") == 0) || 
		//	(strcmp(argv[1],"merge") == 0) || 
		//	(strcmp(argv[1],"split") == 0)
			) {

			g_bUseBQB = true;

			bqbtool_main( argc, argv );

			goto _ende;
		}

		// Call the CubeTool from TRAVIS
		if (strcmp(argv[1],"cubetool") == 0) {

			cubetool_main( argc-1, &argv[1] );

			goto _ende;
		}
	}





	CheckSourceVersion();

	if (g_bControlRun) {

		mprintf("\n");
		mprintf(WHITE,"    #####################################\n");
		mprintf(WHITE,"    ####    This is a Control Run    ####\n");
		mprintf(WHITE,"    #####################################\n\n\n");

		mprintf("    Control runs are a way for direct interfacing between TRAVIS and other program packages.\n");
		mprintf("    This is an experimental feature; there is no support and no documentation on control runs.\n");
		mprintf("    If you are a software developer and would like to have a direct interface to TRAVIS, please\n");
		mprintf("    contact the developers.\n\n\n");

		g_bInputRedirected = true;

		mprintf("    Trying to read database file from \"%s\"...\n",g_sControlFile);

		try { g_pControlDB = new CDatabase(); } catch(...) { g_pControlDB = NULL; }
		if (g_pControlDB == NULL) NewException((double)sizeof(CDatabase),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		g_pControlDB->ParseInputFile(g_sControlFile);

		mprintf("\n    Done. Dump of database tree:\n\n");

		g_pControlDB->DumpTree();

		mprintf("    Trying to determine runtype...\n\n");

		if (!g_pControlDB->ExistNode("/CONTROL")) {
			eprintf("Error: Node \"/CONTROL\" is missing in control file.\n");
			goto _ende;
		}

		cp = g_pControlDB->GetString("/CONTROL/RUNTYPE");

		if (cp == NULL) {
			eprintf("Error: Value \"/CONTROL/RUNTYPE\" is missing or not of string type (%d).\n",g_pControlDB->GetElementType("/CONTROL/RUNTYPE"));
			goto _ende;
		}

		mprintf("    The runtype is \"%s\".\n\n",cp);

		if (mystricmp(cp,"ANALYZE") == 0)
			g_iControlRT = RUNTYPE_ANALYZE;
		else if (mystricmp(cp,"SCANMOL") == 0)
			g_iControlRT = RUNTYPE_SCANMOL;
		else {
			eprintf("Error: Invalid runtype in control file.\n");
			goto _ende;
		}

		mprintf("\n");
		mprintf(WHITE,"    ###################################################\n");
		mprintf(WHITE,"    ####    Control Run Initialization finished    ####\n");
		mprintf(WHITE,"    ###################################################\n\n\n");
	}

	if (!g_bControlRun) {

		g_bInputRedirected = false;
		if (g_sInputFile != NULL)
		{
			g_fInputFile = fopen(g_sInputFile,"rt");
			if (g_fInputFile == NULL)
			{
				eprintf("Could not open input file \"%s\".\n",g_sInputFile);
				goto _ende;
			}
			g_bInputRedirected = true;
		}

		if (IsTTY(stdin) && (!g_bInputRedirected))
		{
			g_fInput = OpenFileWrite("input.txt",true);
			g_bInputRedirected = false;
			inpprintf("! TRAVIS input file\n");
			inpprintf("! Created with TRAVIS version compiled at %s %s\n",__DATE__,__TIME__);
			inpprintf("! Source code version: %s\n",SOURCE_VERSION);

			time(&ltime);
			today = localtime(&ltime);
	//		strcpy(buf,asctime(today));
			buf.strcpy(asctime(today));
			buf((int)(strlen(buf)-1)) = 0;
			inpprintf("! Input file written at %s.\n",(const char*)buf);
		} else
			g_bInputRedirected = true;
	}




	if (g_bDipolGrimme)
	{
		if (argc < 3)
			DipolGrimme(NULL);
		else
			DipolGrimme(argv[2]);
		return 0;
	}

	if (g_bSankey) {
		CSankeyDiagramEngine ske;
		if (strlen(g_sSankeyFile) != 0)
			ske.BuildSankeyDiagram( g_sSankeyFile );
		else
			ske.BuildSankeyDiagram();
		goto _ende;
	}

	
	if (g_bRamanFromPolarizability) {

		for (z=1;z<argc;z++)
			if ((mystricmp(argv[z],"-ramanfrompola") == 0) || (mystricmp(argv[z],"-ramanfrompolarizability") == 0))
				break;

		if (z == argc) {
			eprintf("Internal error: Did not find \"-ramanfrompola\" argument.\n\n");
			abort();
		}

		RamanFromPolarizability( argc-z, &argv[z] );

		// For the literature reference
		g_bRaman = true;

		goto _ende;
	}

		if (!ParseArgs(argc,argv))
			goto _ende;

	SortElementsMass();

	if (g_sInputTraj == NULL)
	{
		eprintf("    No trajectory file specified.\n");
		mprintf("    Please use the \"-p\" flag to specify an input trajectory in the command line.\n\n");

		if ((argc == 1) && !g_bControlRun)
		{
			mprintf("    Enter the file name of the trajectory file to open: [Quit] ");
			myget(&buf);
			mprintf("\n");
			if (strlen(buf)==0)
			{
				CommandLineHelp();
				goto _ende;
			}

			try { g_sInputTraj = new char[strlen(buf)+1]; } catch(...) { g_sInputTraj = NULL; }
			if (g_sInputTraj == NULL) NewException((double)(strlen(buf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

			strcpy(g_sInputTraj,buf);
		} else
		{
			CommandLineHelp();
			goto _ende;
		}
	}

	if (!DetermineTrajFormat())
		goto _ende;

	if (g_iTrajFormat == 6) // Amber
	{
//		sprintf(buf,"%sprmtop",g_sInputTraj);
#ifdef USE_NETCDF
		buf.sprintf("%sparm7",g_sInputTraj);
#else
		buf.sprintf("%sprmtop",g_sInputTraj);
#endif
		if (!FileExist(buf))
		{
			eprintf("\n    Error: File \"%s\" does not exist.\n",(const char*)buf);
#ifdef USE_NETCDF
			printf("    When opening Amber trajectories, both .parm7 and .nc files need to have the same base file name.\n\n");
#else
			printf("    When opening Amber trajectories, both .prmtop and .mdcrd files need to have the same base file name.\n\n");
#endif
			goto _ende;
		}
//		strcpy(g_sAmberParmFile,buf);
		g_sAmberParmFile.strcpy(buf);

//		sprintf(buf,"%smdcrd",g_sInputTraj);
#ifdef USE_NETCDF
		buf.sprintf("%snc",g_sInputTraj);
#else
		buf.sprintf("%smdcrd",g_sInputTraj);
#endif
		if (!FileExist(buf))
		{
			eprintf("\n    Error: File \"%s\" does not exist.\n",(const char*)buf);
#ifdef USE_NETCDF
			printf("    When opening Amber trajectories, both .parm7 and .nc files need to have the same base file name.\n\n");
#else
			printf("    When opening Amber trajectories, both .prmtop and .mdcrd files need to have the same base file name.\n\n");
#endif
			goto _ende;
		}

		try { g_sInputTraj = new char[strlen(buf)+1]; } catch(...) { g_sInputTraj = NULL; }
		if (g_sInputTraj == NULL) NewException((double)(strlen(buf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);

		strcpy(g_sInputTraj,buf);
	}

	if (g_sInputTraj[0] != 0)
	{
		mprintf("    Opening position trajectory %s ...",g_sInputTraj);
		if (!FileExist(g_sInputTraj))
		{
			eprintf("\nError. File does not exist or cannot be read.\n");
			goto _ende;
		}
		mprintf("\n");
	} else {
		eprintf("Error: No trajectory file specified.\n");
		goto _ende;
	}

	if (g_bReadVelocity)
	{
		mprintf("    Opening velocity trajectory %s ...",(const char*)g_sInputVel);
		if (!FileExist(g_sInputVel))
		{
			eprintf("\nError. File does not exist or cannot be read.\n");
			goto _ende;
		}
		if (mystricmp(GetFileExtension((const char*)g_sInputVel),"xyz") != 0) {
			eprintf("\nError. File is not in .xyz format. Only xyz files are currently supported for velocities.\n");
			goto _ende;
		}
		mprintf("\n");

		mprintf("\n");
		mprintf("    You have to enter the unit of the velocities in the data file.\n\n");
		mprintf("    Note: CP2k prints velocities in a.u. = 2.187691*10^6 m/s\n");
		mprintf("          LAMMPS (\"units real\") prints velocities in Angstrom/fs = 1.0*10^5 m/s\n\n");

		g_fVelocityConversion = AskFloat("    One velocity unit in the file corresponds to how many m/s? [2.187691e6] ",2.187691E6);
	}
	if (g_sInputForce[0] != 0)
	{
		mprintf("    Opening force trajectory %s ...",(const char*)g_sInputForce);
		if (!FileExist(g_sInputForce))
		{
			eprintf("\nError. File does not exist or cannot be read.\n");
			goto _ende;
		}
		mprintf("\n");
	}

	if (!OpenInputTrajectory())
		goto _ende;
/*	g_fPos = fopen(g_sInputTraj,"rt");
	if (g_fPos == NULL)
	{
		eprintf("Error. Could not open \"%s\".\n",g_sInputTraj);
		goto _ende;
	}*/

	if (g_iTrajFormat == 10) { // DCD
		g_bDCDFirstFrame = true;
	}

	if (!g_TimeStep.ReadTimestep(g_fPos,true))
	{
		eprintf("\nError reading first time step from trajectory. Leaving.\n");
		goto _ende;
	}

	if (g_iTrajFormat == 10) { // DCD
		if (g_fBoxX == 0)
			mprintf("\n    DCD file does not contain cell vector information (at least one cell vector was zero).\n");
		mprintf("\n    DCD file contains %lu comment lines:\n",(unsigned long)g_saDCDRemarkLines.size());
		for (z=0;z<(int)g_saDCDRemarkLines.size();z++)
			mprintf("      \"%s\"\n",g_saDCDRemarkLines[z].c_str());
	}

	if (g_iTrajFormat == 9) { // .voronoi
		mprintf("\n    This is a CP2k .voronoi trajectory.\n");
		g_bElMagProperties = true;
	}

	if (g_iTrajFormat == 7) {

		g_bUseBQB = true;

		if (g_pBQBFile->GetFrameType() == 4) {

			mprintf("\n    This is an uncompressed BQB trajectory (%dv%d) and contains the following fields:\n",g_pBQBFile->GetFrameType(),g_pBQBFile->GetFrameTypeVersion());
			CBQBTrajectoryFrame btf(*g_pBQBInterface);
			if (!btf.ReadFrame(g_pBQBFile->GetFramePayload())) {
				eprintf("Error: Could not read trajectory frame from BQB payload.\n");
				goto _ende;
			}

			tb = false;
			mprintf("      ");
			if (btf.GetColumn("Label") != NULL) { if (tb) mprintf(", "); mprintf("Labels"); tb=true; }
			if (btf.GetColumn("PosX") != NULL) { if (tb) mprintf(", "); mprintf("Positions"); tb=true; }
			if (btf.GetColumn("VelX") != NULL) { if (tb) mprintf(", "); mprintf("Velocities"); tb=true; }
			if (btf.GetColumn("EChg") != NULL) { if (tb) mprintf(", "); mprintf("El.Charges"); tb=true; }
			if (btf.GetColumn("EDipX") != NULL) { if (tb) mprintf(", "); mprintf("El.Dipoles"); tb=true; g_bElMagProperties = true; }
			if (btf.GetColumn("EQXX") != NULL) { if (tb) mprintf(", "); mprintf("El.Quadrupoles"); tb=true; g_bElMagProperties = true; }
			if (btf.GetColumn("ECurX") != NULL) { if (tb) mprintf(", "); mprintf("El.Currents"); tb=true; g_bElMagProperties = true; }
			if (btf.GetColumn("MDipX") != NULL) { if (tb) mprintf(", "); mprintf("Mag.Dipoles"); g_bElMagProperties = true; }
			mprintf("\n");

		} else if ((g_pBQBFile->GetFrameType() == 5) || (g_pBQBFile->GetFrameType() == 6)) {

			mprintf("\n    This is a compressed BQB trajectory (%dv%d).\n",g_pBQBFile->GetFrameType(),g_pBQBFile->GetFrameTypeVersion());

		} else if (g_pBQBFile->GetFrameType() == 7) {

			mprintf("\n    This is an uncompressed volumetric BQB trajectory (%dv%d).\n",g_pBQBFile->GetFrameType(),g_pBQBFile->GetFrameTypeVersion());

			g_bVolumetricData = true;

		} else if ((g_pBQBFile->GetFrameType() == 8) || (g_pBQBFile->GetFrameType() == 9)) {

			mprintf("\n    This is a compressed volumetric BQB trajectory (%dv%d).\n",g_pBQBFile->GetFrameType(),g_pBQBFile->GetFrameTypeVersion());

			g_bVolumetricData = true;

		} else {

			mprintf("\n");
			eprintf("    Warning: ");
			mprintf("This is a BQB trajectory of unknown format (%dv%d).\n",g_pBQBFile->GetFrameType(),g_pBQBFile->GetFrameTypeVersion());

		}
	}

	if (g_bXYZ4thCol)
		mprintf("\n    Found 4th column of numbers in XYZ file.\n");

	if (g_bLAMMPSCharge)
		mprintf("\n    Found atomic partial charges in LAMMPS file.\n");

	if (g_bLAMMPSForce)
		mprintf("\n    Found atomic forces in LAMMPS trajectory.\n");

	if (g_bXYZComment3Numbers)
		mprintf("\n    Found 3 numeric values in comment line of XYZ file.\n");

	if (g_bXYZComment6Numbers)
		mprintf("\n    Found 6 numeric values in comment line of XYZ file.\n");

	if (g_bXYZComment9Numbers)
		mprintf("\n    Found 9 numeric values in comment line of XYZ file.\n");

	if (g_bFoundNonOrtho)
		mprintf("\n    Found non-orthorhombic cell geometry in trajectory file.\n");

	g_TimeStep.CalcMinMax();

	if (g_iTrajFormat == 7) {

		g_iTrajSteps = g_pBQBFile->GetTotalFrameCount();
		if (g_iTrajSteps == -1)
			mprintf("\n    Could not determine trajectory step count (no BQB index).\n\n");
		else
			mprintf("\n    Trajectory contains %d frames (indexed).\n\n",g_iTrajSteps);

#ifdef USE_NETCDF
	} else if (g_iTrajFormat == 6) {
		int ncid;
		nc_open(g_sInputTraj, NC_NOWRITE, &ncid); 

        	int frame_id;
        	size_t nof;
        	nc_inq_dimid(ncid, "frame", &frame_id);
        	nc_inq_dimlen(ncid, frame_id, &nof);
		g_iTrajSteps = (int)nof;
		mprintf("\n    Trajectory contains %d frames (indexed).\n\n",g_iTrajSteps);
#endif
	} else {

		if (!g_bStreamInput)
		{
		//	fgetpos(g_fPos,&fpos);
			fpos = ftell(g_fPos);
			fseek(g_fPos,0,SEEK_END);
			fpos2 = ftell(g_fPos);
		//	fgetpos(g_fPos,&fpos2);
			fclose(g_fPos);
		} else fpos2 = -1;

		if (fpos2 < 0)
		{
			mprintf("\n    Could not determine trajectory file size.\n\n");
			g_iTrajSteps = -1;
		} else
		{
			if (g_iTrajFormat == 10) // DCD
				g_iTrajSteps = fpos2/(fpos-280);
			else
				g_iTrajSteps = fpos2/fpos;
			mprintf("\n    Trajectory contains approximately %d frames (just a *guess* from the file size of %s).\n\n",g_iTrajSteps,FormatBytes(fpos2));
		}
	}

	g_iGesAtomCount = g_TimeStep.m_iGesAtomCount;

	if (g_iGesAtomCount == 0)
	{
		eprintf("\n\nNo atoms found. This is probably not what you want. Leaving.\n");
		goto _ende;
	}

	g_iGesVirtAtomCount = g_iGesAtomCount;
	g_TimeStep.AddAtoms();
	if (g_bUnknownElements)
		mprintf("\n");

	SortAtoms();
	mprintf(WHITE,"    %d atoms in the system: ",g_iGesAtomCount);
	for (z=0;z<g_oaAtoms.GetSize();z++)
	{
		mprintf("%dx %s",((CAtom*)g_oaAtoms[z])->m_iCount,(const char*)((CAtom*)g_oaAtoms[z])->m_sName);
		if (z < (int)g_oaAtoms.GetSize()-1)
			mprintf(", ");
	}
	mprintf("\n\n");
	mprintf("    System extent:  X = { %+6.0f .. %+6.0f pm }, dX = %.0f pm\n",g_TimeStep.m_vMin[0],g_TimeStep.m_vMax[0],g_TimeStep.m_vMax[0]-g_TimeStep.m_vMin[0]);
	mprintf("      (in step 1)   Y = { %+6.0f .. %+6.0f pm }, dY = %.0f pm\n",g_TimeStep.m_vMin[1],g_TimeStep.m_vMax[1],g_TimeStep.m_vMax[1]-g_TimeStep.m_vMin[1]);
	mprintf("                    Z = { %+6.0f .. %+6.0f pm }, dZ = %.0f pm\n",g_TimeStep.m_vMin[2],g_TimeStep.m_vMax[2],g_TimeStep.m_vMax[2]-g_TimeStep.m_vMin[2]);

	g_iVirtAtomType = (unsigned char)g_oaAtoms.GetSize();

	xAddAtom("#",true);

/*	if (g_iTrajFormat == 3)
	{
		mprintf(YELLOW,"\n*** Important note: ");
		mprintf(" Trajectory file is in LAMMPS format.\n\n");
		mprintf("    TRAVIS requires that the atom ordering (sorting) is identical in each time step of the trajectory.\n");
		mprintf("    However, LAMMPS' defaults may order the atoms randomly in each step when dumping a trajectory.\n");
		mprintf("    To read a LAMMPS trajectory in TRAVIS, you *need* to make sure the atom ordering stays constant.\n");
		mprintf("    This can be achieved by specifying \"dump modify ... sort id\" to your trajectory dump command.\n\n");
		while (!AskYesNo("    I have acknowledged this note (y/n). [yes] ",true)) { }
	}*/

	if ((g_iTrajFormat == 7) && (mystricmp(&g_sInputTraj[strlen(g_sInputTraj)-4],".bbq") == 0))
		PrintBMode();

	/****************************/

	if (!GatherInfos())
		goto _ende;

	if (g_bMinimizeChargeVar2)
		goto _ende2;

	/************** Beginn Analyse **************/

	if (g_bStreamInput)
	{
		if (g_bScanVelocities)
		{
			eprintf("\nError: Velocity pre-analysis not compatible with stream input.\n\n");
			goto _ende;
		}

		if ((g_bSaveRefEnv) && (g_iNbhMode == 3))
		{
			eprintf("\nError: Neighborhood pre-analysis not compatible with stream input.\n\n");
			goto _ende;
		}
	}

	g_iStartTime = time(NULL);

	if (g_iMaxStep != -1)
	{
		g_iTrajSteps = g_iMaxStep;
	} else
	{
		if (g_iTrajSteps != -1)
		{
			g_iTrajSteps -= g_iBeginStep;
			if (g_iTrajSteps < 1)
				g_iTrajSteps = 1;
		}
	}


	if (g_bReRa)
		goto _ende2;

	if (g_bROA) { // This *replaces* the main loop
		g_pROAEngine->MainLoop();
		goto _ende2;
	}

	if (((g_iTrajFormat == 5) || (g_iTrajFormat == 7)) && !g_bCubeTimeDev)
		g_iStepHistory = 1;
	else
		g_iStepHistory = 3;

	multicounter = 0;

_multiintervalstart:
	if (g_bMultiInterval)
	{
		g_iBeginStep = g_laMultiIntervalStart[multicounter];
		g_iScanNbhStart = g_iBeginStep;
		g_iMaxStep = g_laMultiIntervalEnd[multicounter] - g_laMultiIntervalStart[multicounter] + 1;
		g_iScanNbhSteps = g_iMaxStep;
		mprintf(YELLOW,"***********************************************\n");
		mprintf(YELLOW,"*** Interval %2d (%6d - %6d) starting. ***\n",multicounter+1,g_laMultiIntervalStart[multicounter]+1,g_laMultiIntervalEnd[multicounter]+1);
		mprintf(YELLOW,"***********************************************\n\n");
		sprintf(multibuf,"_I%d_%d-%d",multicounter+1,g_laMultiIntervalStart[multicounter],g_laMultiIntervalEnd[multicounter]);
	} else
		multibuf[0] = 0;

	mprintf(WHITE,">>> Initialization >>>\n\n");

	g_bAbortAnalysis = false;

	if (!g_bStreamInput)
	{
		if (g_iTrajFormat == 10)
			g_bDCDSkipHeader = true;
		if (!OpenInputTrajectory())
			goto _ende;
/*		g_fPos = fopen(g_sInputTraj,"rb"); // Eingabedatei erneut oeffnen
		if (g_fPos == NULL)
		{
			eprintf("\nError. Input Trajectory suddenly vanished ^^\n");
			goto _ende;
		}*/
		if ((g_bNPT) && (g_sNPTFile[0] != 0))
		{
			g_fNPTFile = fopen(g_sNPTFile,"rt");
			if (g_fNPTFile == NULL)
			{
				eprintf("\nCould not open cell vector file.\n");
				goto _ende;
			}
		}
	}

	if (g_iBeginStep != 0)
	{
		mprintf("Fast-forwarding to step %d...\n",g_iBeginStep+1);

/*		if (g_iTrajFormat == 7) {

			if (!SeekInputTrajectory(g_iBeginStep))
				goto _ende;

		}*/

		if ((g_iBeginStep != g_iScanMolStep) || (g_iTrajFormat == 7))
		{
			if (g_iTrajFormat == 7)
				mprintf("Warning: Index-based seeking for bqb files not yet implemented.\n");
			mprintf(WHITE,"  [");
			for (z=0;z<g_iBeginStep;z++)
			{
				if (fmod(z,g_iBeginStep/50.0) < 1.0)
					mprintf(WHITE,"#");
				if (!g_TimeStep.SkipTimestep(g_fPos))
				{
					eprintf("\nError: Unexpected end of position trajectory.\n");
					goto _ende;
				}
				if ((g_bNPT) && (g_sNPTFile[0] != 0))
//					(void)!fgets(buf,256,g_fNPTFile);
					(void)!buf.fgets(256,g_fNPTFile);
			}
			mprintf(WHITE,"]\n");
			if (!g_bStreamInput && (g_iTrajFormat != 7))
				g_iFastForwardPos = ftell(g_fPos);
					else g_iFastForwardPos = 0;
		} else
		{
			mprintf("  Position already known from molecule recognition, directly seeking.\n");
			fseek(g_fPos,g_iFastForwardPos,SEEK_SET);
		}
		if (!g_bStreamInput && (g_iTrajFormat != 7))
			mprintf("    Step %d begins at offset %lu (%.1f MB).\n\n",g_iBeginStep+1,g_iFastForwardPos,g_iFastForwardPos/1024.0/1024.0);
	}
//	fclose(g_fPos);
//	g_fPos = fopen(g_sInputTraj,"rb");
//	mprintf("Seek: %d.\n",g_iFastForwardPos);
//	fseek(g_fPos,g_iFastForwardPos,SEEK_SET);

	if (!g_bStreamInput)
	{
		g_TimeStep.ReadTimestep(g_fPos,true);
		if ((g_bNPT) && (g_sNPTFile[0] != 0))
		{
			if (!g_TimeStep.ReadCellVector(g_fNPTFile))
				abort();
			fclose(g_fNPTFile);
		}
	//	fclose(g_fPos);
		if (!CloseInputTrajectory())
			goto _ende;
	}

	g_TimeStep.CalcCenters();

	if (g_bPeriodic)
		g_TimeStep.UniteMolecules(false);

	if (g_bScanVelocities && (!g_bStreamInput))
	{
		mprintf(WHITE,"\n>>> Pre-Analysis for velocity distribution >>>\n\n");

/*		g_fPos = fopen(g_sInputTraj,"rt"); // Eingabedatei erneut Oeffnen
		if (g_fPos == NULL)
		{
			eprintf("Error. Input Trajectory suddenly vanished ^^\n");
			return 0;
		}*/
		if (!OpenInputTrajectory())
			goto _ende;
		g_iSteps = 0; // Der Zaehler der Zeitschritte
		if (g_iScanVelStart != 0)
		{
			mprintf("Fast-forwarding to step %d...\n",g_iScanVelStart+1);
			if (g_iTrajFormat == 7)
				mprintf("Warning: Index-based seeking for bqb files not yet implemented.\n");
			mprintf(WHITE,"  [");
			for (z=0;z<g_iScanVelStart;z++)
			{
				if (fmod(z,g_iScanVelStart/60.0) < 1.0)
					mprintf(WHITE,"#");
				if (!g_TimeStep.SkipTimestep(g_fPos))
					break;
			}
			mprintf(WHITE,"]\n");
		}

//		while (!feof(g_fPos)) // Zeitschritt fuer Zeitschritt die Trajektorie durchgehen
		while (!InputTrajectoryEOF()) // Zeitschritt fuer Zeitschritt die Trajektorie durchgehen
		{
			if (!g_TimeStep.ReadTimestep(g_fPos,false))
				goto _endvel;
			if (g_TimeStep.m_iGesAtomCount == 0)
				goto _endvel;

			g_TimeStep.CalcCenters();
	
			if ((g_iSteps % 4) == 0)
			{
				if ((g_iSteps % 200) == 0) 
					mprintf("\nStep %6lu...",g_iSteps);
						else mprintf(".");
			}
	
			g_iSteps++;
	
			if ((g_iScanVelSteps > 0) && ((int)g_iSteps >= g_iScanVelSteps))
				break;
		}
		_endvel:

//		fclose(g_fPos);
		if (!CloseInputTrajectory())
			goto _ende;

		mprintf(WHITE,"\n<<< End of Pre-Analysis for velocity distribution <<<\n\n");
	} // END IF g_bScanVelocities


	if (g_bSaveRefEnv && (g_iNbhMode == 1))
	{
		mprintf("Creating statical neighborhood...\n");
		g_pNbSet->Scan((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[g_iSaveRefMol]],&g_TimeStep);
		if (g_bSaveRefWithEnv)
		{
			mprintf("Adding reference molecule to neighborhood...\n");
			g_pNbSet->AddMolecule(g_iFixMol,g_iSaveRefMol);
		}
		mprintf("\n");
		g_pNbSet->Dump();
		mprintf("\n");
	}

	if (g_bVACF && g_bGlobalVACF)
	{
		mprintf("Initializing Global Velocity Autocorrelation Function...\n");
		g_pGlobalVACF->Create();

		if (!g_bVACFCacheMode)
		{
			if (g_pGlobalVACF->m_iSize > g_iStepHistory)
				g_iStepHistory = g_pGlobalVACF->m_iSize;
		} else
		{
			mprintf("    VACF Cache: Trying to allocate %s of memory...\n",FormatBytes((double)g_iGesAtomCount*g_iTrajSteps*3.1*sizeof(double)));
			for (z=0;z<g_iGesAtomCount;z++)
			{
//				mprintf("%d: %s: %f\n",z,((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_sLabel,((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_fRadius);
				if (g_pGlobalVACF->m_bExcludeR0 && (g_waAtomRealElement[z] > 1000))
					continue;
//				mprintf("  --> Ok\n");

				try { ptfa = new CxDoubleArray("main():ptfa"); } catch(...) { ptfa = NULL; }
				if (ptfa == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				if (g_iTrajSteps != -1)
				{
					ptfa->SetMaxSize((long)(g_iTrajSteps*3.1));
					ptfa->SetGrow((long)(g_iTrajSteps*0.1));
				} else ptfa->SetGrow(1000);
				g_pGlobalVACF->m_oaCache.Add(ptfa);
			}
		}
	}

	if (g_bIRSpec && g_bGlobalIR)
	{
		mprintf("  Creating global IR spectrum...\n");
		g_pGlobalIR->m_pRDyn->m_fMinVal = 0;
		g_pGlobalIR->m_pRDyn->m_fMaxVal = g_pGlobalIR->m_iDepth * g_fTimestepLength / 1000.0;
		g_pGlobalIR->m_pRDyn->m_iResolution = g_pGlobalIR->m_iDepth/g_pGlobalIR->m_iStride;
		g_pGlobalIR->m_pRDyn->SetLabelX("Tau / ps");
		g_pGlobalIR->m_pRDyn->SetLabelY("Dipole autocorrelation");
		g_pGlobalIR->m_pRDyn->Create();
		if (g_bRDynCacheMode)
		{
			if (g_iTrajSteps != -1)
				mprintf("    RDyn Cache: Trying to allocate %s of memory...\n",FormatBytes((double)g_oaSingleMolecules.GetSize()*g_iTrajSteps/g_iStride*3.1*sizeof(double)));
					else mprintf("    RDyn Cache: Trying to allocate %s of memory...\n",FormatBytes((double)g_oaSingleMolecules.GetSize()*10000/g_iStride*3.1*sizeof(double)));
			for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++)
			{
				try { ptfa = new CxDoubleArray("main():ptfa"); } catch(...) { ptfa = NULL; }
				if (ptfa == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				if (g_iTrajSteps != -1)
				{
					ptfa->SetMaxSize((long)(g_iTrajSteps/g_iStride*3.3));
					ptfa->SetGrow((long)(g_iTrajSteps/g_iStride*0.3));
				} else ptfa->SetGrow(10000);
				g_pGlobalIR->m_oaCache.Add(ptfa);
			}
		} else
		{
			try { g_pGlobalIR->m_pCount = new double[g_pGlobalIR->m_pRDyn->m_iResolution]; } catch(...) { g_pGlobalIR->m_pCount = NULL; }
			if (g_pGlobalIR->m_pCount == NULL) NewException((double)g_pGlobalIR->m_pRDyn->m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			for (z2=0;z2<g_pGlobalIR->m_pRDyn->m_iResolution;z2++)
				g_pGlobalIR->m_pCount[z2] = 0;
			if (!g_bRDynCacheMode)
			{
				if (g_pGlobalIR->m_iDepth > g_iStepHistory)
					g_iStepHistory = g_pGlobalIR->m_iDepth;
			}
		}
	}

/*	if (g_bSFac)
	{
		mprintf("  Creating Structure Factor Analysis...\n");
		g_pSFac->Create();
	}*/

	if (g_bSaveCondSnapshot)
		g_fSaveCondFile = OpenFileWrite("savecondition.xyz",true);

	for (z=0;z<g_oaObserv.GetSize();z++)
	{
		mprintf("Initializing observation %d...\n",z+1);
		o = (CObservation*)g_oaObserv[z];

		if (o->m_bPercTimeDev) {
			buf.sprintf("fraction_timedev_%d_%s.csv",z+1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
			mprintf("  Opening Fraction Temporal Development file \"%s\" ...\n",(const char*)buf);
			o->m_pPercTimeDevFile = OpenFileWrite((const char*)buf,true);
			mfprintf(
				o->m_pPercTimeDevFile,
				"# Step;  Count;  Percentage\n"
			);
		}

		if (g_bAggregation)
		{
			mprintf("  Creating aggregation functions (%d value sets)...\n",o->m_pDACF->m_oaSubDACFs.GetSize());
			for (z2=0;z2<o->m_pDACF->m_oaSubDACFs.GetSize();z2++)
			{
				dacfsub = (CDACFSub*)o->m_pDACF->m_oaSubDACFs[z2];
				if (g_bDLDisp)
				{
					try { dacfsub->m_pDLDisp = new C2DF(); } catch(...) { dacfsub->m_pDLDisp = NULL; }
					if (dacfsub->m_pDLDisp == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					dacfsub->m_pDLDisp->m_iRes[0] = o->m_pDACF->m_iLifetimeRes;
					dacfsub->m_pDLDisp->m_fMinVal[0] = 0.0;
					dacfsub->m_pDLDisp->m_fMaxVal[0] = o->m_pDACF->m_fLargestLifetime;
					dacfsub->m_pDLDisp->m_iRes[1] = o->m_pDACF->m_iDisplacementRes;
					dacfsub->m_pDLDisp->m_fMinVal[1] = 0.0;
					dacfsub->m_pDLDisp->m_fMaxVal[1] = o->m_pDACF->m_fLargestDisplacement;
					dacfsub->m_pDLDisp->Create();
					dacfsub->m_pDLDisp->SetLabelX("Lifetime / ps");
					dacfsub->m_pDLDisp->SetLabelY("Displacement / pm");
				}

				if (g_bPairMSD)
				{
					try { dacfsub->m_pPairMSD = new CAF(); } catch(...) { dacfsub->m_pPairMSD = NULL; }
					if (dacfsub->m_pPairMSD == NULL) NewException((double)sizeof(CAF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					dacfsub->m_pPairMSD->m_iResolution = o->m_pDACF->m_iDACFRes;
					dacfsub->m_pPairMSD->m_fMinVal = 0.0;
					dacfsub->m_pPairMSD->m_fMaxVal = o->m_pDACF->m_fLargestLifetime;
					dacfsub->m_pPairMSD->Create();
				}

				if (g_bDDisp)
				{
					try { dacfsub->m_pDDisp = new CDF(); } catch(...) { dacfsub->m_pDDisp = NULL; }
					if (dacfsub->m_pDDisp == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					dacfsub->m_pDDisp->m_iResolution = o->m_pDACF->m_iDisplacementRes;
					dacfsub->m_pDDisp->m_fMinVal = 0.0;
					dacfsub->m_pDDisp->m_fMaxVal = o->m_pDACF->m_fLargestDisplacement;
					dacfsub->m_pDDisp->SetLabelX("Dimer displacement / pm");
					dacfsub->m_pDDisp->SetLabelY("Occurrence");
					dacfsub->m_pDDisp->Create();
				}

				if (g_bDLDF)
				{
					try { dacfsub->m_pDLDF = new CDF(); } catch(...) { dacfsub->m_pDLDF = NULL; }
					if (dacfsub->m_pDLDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					dacfsub->m_pDLDF->m_bLeft = true;
					dacfsub->m_pDLDF->m_iResolution = o->m_pDACF->m_iLifetimeRes;
					dacfsub->m_pDLDF->m_fMinVal = 0.0;
					dacfsub->m_pDLDF->m_fMaxVal = o->m_pDACF->m_fLargestLifetime;
					dacfsub->m_pDLDF->SetLabelX("Tau / ps");
					dacfsub->m_pDLDF->SetLabelY("Occurrence");
					dacfsub->m_pDLDF->Create();
				}

				if (g_bDACF)
				{
					try { dacfsub->m_pDACF = new CDF(); } catch(...) { dacfsub->m_pDACF = NULL; }
					if (dacfsub->m_pDACF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					dacfsub->m_pDACF->m_bLeft = true;
					dacfsub->m_pDACF->m_iResolution = o->m_pDACF->m_iDACFRes;
					dacfsub->m_pDACF->m_fMinVal = 0;
					dacfsub->m_pDACF->m_fMaxVal = o->m_pDACF->m_iDACFRes * g_fTimestepLength / 1000.0;
					dacfsub->m_pDACF->SetLabelX("Tau / ps");
					dacfsub->m_pDACF->SetLabelY("Occurrence");
					dacfsub->m_pDACF->Create();

					if (dacfsub->m_bNewMode)
					{
				//		mprintf("Create: %d*%d=%d\n",((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize(),((CMolecule*)g_oaMolecules[dacfsub->m_iShowMol])->m_laSingleMolIndex.GetSize(),((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize()*((CMolecule*)g_oaMolecules[dacfsub->m_iShowMol])->m_laSingleMolIndex.GetSize());
						try { dacfsub->m_piaIntervals = new CxIntArray[((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize()*((CMolecule*)g_oaMolecules[dacfsub->m_iShowMol])->m_laSingleMolIndex.GetSize()]; } catch(...) { dacfsub->m_piaIntervals = NULL; }
						if (dacfsub->m_piaIntervals == NULL) NewException((double)((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize()*((CMolecule*)g_oaMolecules[dacfsub->m_iShowMol])->m_laSingleMolIndex.GetSize()*sizeof(CxIntArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				}

				try { dacfsub->m_pNDF = new CDF(); } catch(...) { dacfsub->m_pNDF = NULL; }
				if (dacfsub->m_pNDF == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				dacfsub->m_pNDF->m_iResolution = ((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+1;
				dacfsub->m_pNDF->m_fMinVal = 0.0;
				dacfsub->m_pNDF->m_fMaxVal = ((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+1;
				dacfsub->m_pNDF->Create();
			}
		}

		if (g_bSDF)
		{
			mprintf("  Creating SDF...\n");
			o->m_pSDF->m_pSDF->m_fMinVal[0] = -o->m_pSDF->m_fRadius;
			o->m_pSDF->m_pSDF->m_fMaxVal[0] = o->m_pSDF->m_fRadius;
			o->m_pSDF->m_pSDF->m_fMinVal[1] = -o->m_pSDF->m_fRadius;
			o->m_pSDF->m_pSDF->m_fMaxVal[1] = o->m_pSDF->m_fRadius;
			o->m_pSDF->m_pSDF->m_fMinVal[2] = -o->m_pSDF->m_fRadius;
			o->m_pSDF->m_pSDF->m_fMaxVal[2] = o->m_pSDF->m_fRadius;
			o->m_pSDF->m_pSDF->m_iRes[0] = o->m_pSDF->m_iResolution;
			o->m_pSDF->m_pSDF->m_iRes[1] = o->m_pSDF->m_iResolution;
			o->m_pSDF->m_pSDF->m_iRes[2] = o->m_pSDF->m_iResolution;
			o->m_pSDF->m_pSDF->m_iHistogramRes = o->m_pSDF->m_iHistogramRes;
			o->m_pSDF->m_pSDF->Create();
//			mprintf("Observation %d: Creating %d VecArrays.\n",z+1,o->m_iShowMolCount);

			try { o->m_pSDF->m_vaData = new CxDVec3Array[o->m_iShowMolCount]; } catch(...) { o->m_pSDF->m_vaData = NULL; }
			if (o->m_pSDF->m_vaData == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDVec3Array),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			if (o->m_pSDF->m_bVdWSpheres)
			{
				try { o->m_pSDF->m_faRadius = new CxDoubleArray[o->m_iShowMolCount]; } catch(...) { o->m_pSDF->m_faRadius = NULL; }
				if (o->m_pSDF->m_faRadius == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}

			if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
			{
				try { o->m_pSDF->m_baDataEnabled = new CxByteArray[o->m_iShowMolCount]; } catch(...) { o->m_pSDF->m_baDataEnabled = NULL; }
				if (o->m_pSDF->m_baDataEnabled == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}
		}

		if (g_bPlProj)
		{
			mprintf("  Creating plane projection DF...\n");
			o->m_pPlProj->m_p2DF->m_fMinVal[0] = o->m_pPlProj->m_fMinVal[0];
			o->m_pPlProj->m_p2DF->m_fMaxVal[0] = o->m_pPlProj->m_fMaxVal[0];
			o->m_pPlProj->m_p2DF->m_fMinVal[1] = o->m_pPlProj->m_fMinVal[1];
			o->m_pPlProj->m_p2DF->m_fMaxVal[1] = o->m_pPlProj->m_fMaxVal[1];
			o->m_pPlProj->m_p2DF->m_iRes[0] = o->m_pPlProj->m_iResolution[0];
			o->m_pPlProj->m_p2DF->m_iRes[1] = o->m_pPlProj->m_iResolution[1];
			o->m_pPlProj->m_p2DF->m_iHistogramRes = o->m_pPlProj->m_iHistogramRes;
			o->m_pPlProj->m_p2DF->SetLabelX("X / pm");
			o->m_pPlProj->m_p2DF->SetLabelY("Y / pm");
			o->m_pPlProj->m_p2DF->Create();
//			mprintf("Observation %d: Creating %d VecArrays.\n",z+1,o->m_iShowMolCount);

			try { o->m_pPlProj->m_vaData = new CxDVec3Array[o->m_iShowMolCount]; } catch(...) { o->m_pPlProj->m_vaData = NULL; }
			if (o->m_pPlProj->m_vaData == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDVec3Array),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
			{
				try { o->m_pPlProj->m_baDataEnabled = new CxByteArray[o->m_iShowMolCount]; } catch(...) { o->m_pPlProj->m_baDataEnabled = NULL; }
				if (o->m_pPlProj->m_baDataEnabled == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}

			if (o->m_pPlProj->m_bVector) {
				o->m_pPlProj->m_p2DF->m_bContourLines = false;
				o->m_pPlProj->m_faVectorData.resize(o->m_pPlProj->m_iVecRes[0]*o->m_pPlProj->m_iVecRes[1]*2);
				o->m_pPlProj->m_faVectorCount.resize(o->m_pPlProj->m_iVecRes[0]*o->m_pPlProj->m_iVecRes[1]);
				for (z2=0;z2<o->m_pPlProj->m_iVecRes[0]*o->m_pPlProj->m_iVecRes[1]*2;z2++)
					o->m_pPlProj->m_faVectorData[z2] = 0;
				for (z2=0;z2<o->m_pPlProj->m_iVecRes[0]*o->m_pPlProj->m_iVecRes[1];z2++)
					o->m_pPlProj->m_faVectorCount[z2] = 0;
				try { o->m_pPlProj->m_vaVectorData = new CxDVec3Array[o->m_iShowMolCount]; } catch(...) { o->m_pPlProj->m_vaVectorData = NULL; }
				if (o->m_pPlProj->m_vaVectorData == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDVec3Array),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}
		}

		if (g_bRDyn)
		{
			mprintf("  Creating reorientation dynamics...\n");
			o->m_pRDyn->m_pRDyn->m_fMinVal = 0;
			o->m_pRDyn->m_pRDyn->m_fMaxVal = o->m_pRDyn->m_iDepth * g_fTimestepLength / 1000.0;
			o->m_pRDyn->m_pRDyn->m_iResolution = o->m_pRDyn->m_iDepth/o->m_pRDyn->m_iStride;
			o->m_pRDyn->m_pRDyn->SetLabelX("Tau / ps");
			o->m_pRDyn->m_pRDyn->SetLabelY("Vector autocorrelation");
			o->m_pRDyn->m_pRDyn->Create();
			if (g_bRDynCacheMode)
			{
				if (g_iTrajSteps != -1)
					mprintf("    RDyn Cache: Trying to allocate %s of memory...\n",FormatBytes((double)((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pRDyn->m_iCombinations*g_iTrajSteps/g_iStride*3.1*sizeof(double)));
						else mprintf("    RDyn Cache: Trying to allocate %s of memory...\n",FormatBytes((double)((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pRDyn->m_iCombinations*10000/g_iStride*3.1*sizeof(double)));
				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pRDyn->m_iCombinations;z2++)
				{
					try { ptfa = new CxDoubleArray("main():ptfa"); } catch(...) { ptfa = NULL; }
					if (ptfa == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (g_iTrajSteps != -1)
					{
						ptfa->SetMaxSize((long)(g_iTrajSteps/g_iStride*3.3));
						ptfa->SetGrow((long)(g_iTrajSteps/g_iStride*0.3));
					} else ptfa->SetGrow(10000);
					o->m_pRDyn->m_oaCache.Add(ptfa);
				}
			} else
			{
				try { o->m_pRDyn->m_pCount = new double[o->m_pRDyn->m_pRDyn->m_iResolution]; } catch(...) { o->m_pRDyn->m_pCount = NULL; }
				if (o->m_pRDyn->m_pCount == NULL) NewException((double)o->m_pRDyn->m_pRDyn->m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				for (z2=0;z2<o->m_pRDyn->m_pRDyn->m_iResolution;z2++)
					o->m_pRDyn->m_pCount[z2] = 0;
				if (!g_bRDynCacheMode)
				{
					if (o->m_pRDyn->m_iDepth > g_iStepHistory)
						g_iStepHistory = o->m_pRDyn->m_iDepth;
				}
			}
		}

		if (g_bIRSpec)
		{
			mprintf("  Creating IR spectrum...\n");
			o->m_pIRSpec->m_pRDyn->m_fMinVal = 0;
			o->m_pIRSpec->m_pRDyn->m_fMaxVal = o->m_pIRSpec->m_iDepth * g_fTimestepLength / 1000.0;
			o->m_pIRSpec->m_pRDyn->m_iResolution = o->m_pIRSpec->m_iDepth/o->m_pIRSpec->m_iStride;
			o->m_pIRSpec->m_pRDyn->SetLabelX("Tau / ps");
			o->m_pIRSpec->m_pRDyn->SetLabelY("Dipole autocorrelation");
			o->m_pIRSpec->m_pRDyn->Create();
			if (g_bRDynCacheMode)
			{
				if (g_iTrajSteps != -1)
					mprintf("    RDyn Cache: Trying to allocate %s of memory...\n",FormatBytes((double)((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pIRSpec->m_iCombinations*g_iTrajSteps/g_iStride*3.1*sizeof(double)));
						else mprintf("    RDyn Cache: Trying to allocate %s of memory...\n",FormatBytes((double)((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pIRSpec->m_iCombinations*10000/g_iStride*3.1*sizeof(double)));
				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pIRSpec->m_iCombinations;z2++)
				{
					try { ptfa = new CxDoubleArray("main():ptfa"); } catch(...) { ptfa = NULL; }
					if (ptfa == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (g_iTrajSteps != -1)
					{
						ptfa->SetMaxSize((long)(g_iTrajSteps/g_iStride*3.3));
						ptfa->SetGrow((long)(g_iTrajSteps/g_iStride*0.3));
					} else ptfa->SetGrow(10000);
					o->m_pIRSpec->m_oaCache.Add(ptfa);
				}
			} else
			{
				try { o->m_pIRSpec->m_pCount = new double[o->m_pIRSpec->m_pRDyn->m_iResolution]; } catch(...) { o->m_pIRSpec->m_pCount = NULL; }
				if (o->m_pIRSpec->m_pCount == NULL) NewException((double)o->m_pIRSpec->m_pRDyn->m_iResolution*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				for (z2=0;z2<o->m_pIRSpec->m_pRDyn->m_iResolution;z2++)
					o->m_pIRSpec->m_pCount[z2] = 0;
				if (!g_bRDynCacheMode)
				{
					if (o->m_pIRSpec->m_iDepth > g_iStepHistory)
						g_iStepHistory = o->m_pIRSpec->m_iDepth;
				}
			}
		}

		if (g_bDens)
		{
			mprintf("  Creating density distribution function...\n");
			o->m_pDensityDF->m_pDensDF->m_fMinVal = o->m_pDensityDF->m_fMinDist;
			o->m_pDensityDF->m_pDensDF->m_fMaxVal = o->m_pDensityDF->m_fMaxDist;
			o->m_pDensityDF->m_pDensDF->m_iResolution = o->m_pDensityDF->m_iResolution;
			o->m_pDensityDF->m_pDensDF->m_iHistogramRes = o->m_pDensityDF->m_iHistogramRes;
			o->m_pDensityDF->m_pDensDF->SetLabelX("Distance / pm");
			if (o->m_pDensityDF->m_bDensityMass)
				o->m_pDensityDF->m_pDensDF->SetLabelY("Mass Density / g cm^-3");
					else o->m_pDensityDF->m_pDensDF->SetLabelY("Particle Density  / nm^-3");
			o->m_pDensityDF->m_pDensDF->Create();
		}

		if (g_bVHDF)
		{
			mprintf("  Creating VHCF...\n");
			o->m_pVHDF->m_pVHDF->m_fMinVal[1] = o->m_pVHDF->m_fMinDist;
			o->m_pVHDF->m_pVHDF->m_fMaxVal[1] = o->m_pVHDF->m_fMaxDist;
			o->m_pVHDF->m_pVHDF->m_iRes[1] = o->m_pVHDF->m_iResolution;
			o->m_pVHDF->m_pVHDF->m_fMinVal[0] = 0.0;
			o->m_pVHDF->m_pVHDF->m_fMaxVal[0] = o->m_pVHDF->m_iDepth * g_fTimestepLength / 1000.0;
			o->m_pVHDF->m_pVHDF->m_iRes[0] = o->m_pVHDF->m_iDepth / o->m_pVHDF->m_iStride;
			o->m_pVHDF->m_pVHDF->Create();

			try { o->m_pVHDF->m_pCount = new double[o->m_pVHDF->m_pVHDF->m_iRes[0]]; } catch(...) { o->m_pVHDF->m_pCount = NULL; }
			if (o->m_pVHDF->m_pCount == NULL) NewException((double)o->m_pVHDF->m_pVHDF->m_iRes[0]*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			for (z0=0;z0<o->m_pVHDF->m_pVHDF->m_iRes[0];z0++)
				o->m_pVHDF->m_pCount[z0] = 0;
			o->m_pVHDF->m_pVHDF->SetLabelY("Distance / pm");
			o->m_pVHDF->m_pVHDF->SetLabelX("Tau / ps");
//			o->m_pVHDF->m_pVHDF->m_fPlotExp = 1.0;
			if (o->m_pVHDF->m_iDepth > g_iStepHistory)
				g_iStepHistory = o->m_pVHDF->m_iDepth;
			mprintf("    Setting trajectory ring buffer to %s...\n",FormatBytes((double)(g_iStepHistory*g_iGesVirtAtomCount*sizeof(double)*3)));
		} 

		if (g_bNbAnalysis)
		{
			mprintf("  Creating neighborhood analysis...\n");

			try { o->m_pNbAnalysis->m_pNPFCount = new CDF(); } catch(...) { o->m_pNbAnalysis->m_pNPFCount = NULL; }
			if (o->m_pNbAnalysis->m_pNPFCount == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			o->m_pNbAnalysis->m_pNPFCount->m_fMinVal = o->m_pNbAnalysis->m_fMinDist;
			o->m_pNbAnalysis->m_pNPFCount->m_fMaxVal = o->m_pNbAnalysis->m_fMaxDist;
			o->m_pNbAnalysis->m_pNPFCount->m_iResolution = o->m_pNbAnalysis->m_iResolution;

			o->m_pNbAnalysis->m_pNPFCount->Create();

			try { tdf = new CDF(); } catch(...) { tdf = NULL; }
			if (tdf == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			tdf->m_fMinVal = o->m_pNbAnalysis->m_fMinDist;
			tdf->m_fMaxVal = o->m_pNbAnalysis->m_fMaxDist;
			tdf->m_iResolution = o->m_pNbAnalysis->m_iResolution;
			tdf->SetLabelX("Distance / pm");
			tdf->SetLabelY("Occurrence");
			tdf->Create();
			o->m_pNbAnalysis->m_oaNPF.Add(tdf);

			for (z2=0;z2<=o->m_pNbAnalysis->m_iNbCount;z2++)
			{
				try { tdf = new CDF(); } catch(...) { tdf = NULL; }
				if (tdf == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				tdf->m_fMinVal = o->m_pNbAnalysis->m_fMinDist;
				tdf->m_fMaxVal = o->m_pNbAnalysis->m_fMaxDist;
				tdf->m_iResolution = o->m_pNbAnalysis->m_iResolution;
				tdf->SetLabelX("Distance / pm");
				tdf->SetLabelY("Occurrence");
				tdf->Create();
				o->m_pNbAnalysis->m_oaDF.Add(tdf);

				try { tdf = new CDF(); } catch(...) { tdf = NULL; }
				if (tdf == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				tdf->m_fMinVal = o->m_pNbAnalysis->m_fMinDist;
				tdf->m_fMaxVal = o->m_pNbAnalysis->m_fMaxDist;
				tdf->m_iResolution = o->m_pNbAnalysis->m_iResolution;
				tdf->SetLabelX("Distance / pm");
				tdf->SetLabelY("Occurrence");
				tdf->Create();
				o->m_pNbAnalysis->m_oaNPF.Add(tdf);
			}
		} 

		for (z0=0;z0<g_iCDFChannels;z0++)
		{
			if (g_bRDF && (o->m_pRDF[z0] != NULL))
			{
				mprintf("  Creating RDF...\n");
				o->m_pRDF[z0]->m_pRDF->m_fMinVal = o->m_pRDF[z0]->m_fMinDist;
				o->m_pRDF[z0]->m_pRDF->m_fMaxVal = o->m_pRDF[z0]->m_fMaxDist;
				o->m_pRDF[z0]->m_pRDF->m_iResolution = o->m_pRDF[z0]->m_iResolution;
				o->m_pRDF[z0]->m_pRDF->m_iHistogramRes = o->m_pRDF[z0]->m_iHistogramRes;
				o->m_pRDF[z0]->m_pRDF->SetLabelX("Distance / pm");
				if (o->m_pRDF[z0]->m_bProbDens)
					o->m_pRDF[z0]->m_pRDF->SetLabelY("g(r) / nm^-3");
						else o->m_pRDF[z0]->m_pRDF->SetLabelY("g(r)");

				if (o->m_bObsCertain && o->m_bDecompDist)
				{
					o->m_pRDF[z0]->m_pRDF->CreateMulti(o->m_waObsRefList.GetSize()*o->m_waObsShowList.GetSize());

			/*		try { o->m_pRDF[z0]->m_pRDF->m_sLabelMulti = new char*[o->m_waObsRefList.GetSize()*o->m_waObsShowList.GetSize()]; } catch(...) { o->m_pRDF[z0]->m_pRDF->m_sLabelMulti = NULL; }
					if (o->m_pRDF[z0]->m_pRDF->m_sLabelMulti == NULL) NewException((double)o->m_waObsRefList.GetSize()*o->m_waObsShowList.GetSize()*sizeof(char*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					for (z=0;z<o->m_waObsRefList.GetSize()*o->m_waObsShowList.GetSize();z++)
					{
						try { o->m_pRDF[z0]->m_pRDF->m_sLabelMulti[z] = new char[128]; } catch(...) { o->m_pRDF[z0]->m_pRDF->m_sLabelMulti[z] = NULL; }
						if (o->m_pRDF[z0]->m_pRDF->m_sLabelMulti[z] == NULL) NewException((double)128*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
			*/
					for (z3=0;z3<o->m_waObsRefList.GetSize();z3++)
					{
						for (z2=0;z2<o->m_waObsShowList.GetSize();z2++)
						{
							if ((o->m_waObsRefList.GetSize() > 1) && (o->m_waObsShowList.GetSize() > 1))
//								sprintf(buf,"RM%d-OM%d",o->m_waObsRefList[z3]+1,o->m_waObsShowList[z2]+1);
								buf.sprintf("RM%d-OM%d",o->m_waObsRefList[z3]+1,o->m_waObsShowList[z2]+1);
							else if (o->m_waObsRefList.GetSize() > 1)
//								sprintf(buf,"RM%d",o->m_waObsRefList[z3]+1);
								buf.sprintf("RM%d",o->m_waObsRefList[z3]+1);
							else 
//								sprintf(buf,"OM%d",o->m_waObsShowList[z2]+1);
								buf.sprintf("OM%d",o->m_waObsShowList[z2]+1);

							o->m_pRDF[z0]->m_pRDF->SetLabelMulti(z3*o->m_waObsShowList.GetSize()+z2,buf);
						}
					}
				} else if (o->m_bDecompType)
				{
					m = (CMolecule*)g_oaMolecules[g_iFixMol];
					o->m_waDecompTypeRefOffs.SetSize(g_iGesVirtAtomCount);
					for (z3=0;z3<o->m_pRDF[z0]->m_oaVectors.GetSize()/2;z3++)
					{
						ag = (CAtomGroup*)o->m_pRDF[z0]->m_oaVectors[z3*2];
						for (z2=0;z2<ag->m_baRealAtomType.GetSize();z2++)
						{
							for (z3=0;z3<o->m_waDecompTypeRefList.GetSize();z3++)
								if (o->m_waDecompTypeRefList[z3] == ag->m_baRealAtomType[z2])
									goto _decomptype1;
							z3 = o->m_waDecompTypeRefList.GetSize();
							o->m_waDecompTypeRefList.Add(ag->m_baRealAtomType[z2]);
_decomptype1:
							for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
							{
								sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
								for (z5=0;z5<((CxIntArray*)ag->m_oaAtoms[z2])->GetSize();z5++)
									o->m_waDecompTypeRefOffs[((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z5))] = (unsigned short)z3;
							}
						}
					}

					if (o->m_bOthers)
						m = (CMolecule*)g_oaMolecules[o->m_iShowMol];
							else m = (CMolecule*)g_oaMolecules[g_iFixMol];
					o->m_waDecompTypeObsOffs.SetSize(g_iGesVirtAtomCount);
					for (z6=0;z6<o->m_pRDF[z0]->m_oaVectors.GetSize()/2;z6++)
					{
						ag = (CAtomGroup*)o->m_pRDF[z0]->m_oaVectors[z6*2+1];
						for (z2=0;z2<ag->m_baRealAtomType.GetSize();z2++)
						{
							for (z3=0;z3<o->m_waDecompTypeObsList.GetSize();z3++)
								if (o->m_waDecompTypeObsList[z3] == ag->m_baRealAtomType[z2])
									goto _decomptype2;
							z3 = o->m_waDecompTypeObsList.GetSize();
							o->m_waDecompTypeObsList.Add(ag->m_baRealAtomType[z2]);
_decomptype2:
							for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++)
							{
								sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
								for (z5=0;z5<((CxIntArray*)ag->m_oaAtoms[z2])->GetSize();z5++)
									o->m_waDecompTypeObsOffs[((CxIntArray*)sm->m_oaAtomOffset[ag->m_baAtomType[z2]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z2])->GetAt(z5))] = (unsigned short)z3;
							}
						}
					}

					o->m_pRDF[z0]->m_pRDF->CreateMulti(o->m_waDecompTypeRefList.GetSize()*o->m_waDecompTypeObsList.GetSize());

	/*				try { o->m_pRDF[z0]->m_pRDF->m_sLabelMulti = new char*[o->m_waDecompTypeRefList.GetSize()*o->m_waDecompTypeObsList.GetSize()]; } catch(...) { o->m_pRDF[z0]->m_pRDF->m_sLabelMulti = NULL; }
					if (o->m_pRDF[z0]->m_pRDF->m_sLabelMulti == NULL) NewException((double)o->m_waDecompTypeRefList.GetSize()*o->m_waDecompTypeObsList.GetSize()*sizeof(char*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					for (z=0;z<o->m_waDecompTypeRefList.GetSize()*o->m_waDecompTypeObsList.GetSize();z++)
					{
						try { o->m_pRDF[z0]->m_pRDF->m_sLabelMulti[z] = new char[128]; } catch(...) { o->m_pRDF[z0]->m_pRDF->m_sLabelMulti[z] = NULL; }
						if (o->m_pRDF[z0]->m_pRDF->m_sLabelMulti[z] == NULL) NewException((double)128*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}*/

					for (z3=0;z3<o->m_waDecompTypeRefList.GetSize();z3++)
					{
						for (z2=0;z2<o->m_waDecompTypeObsList.GetSize();z2++)
						{
//							sprintf(buf,"%s-%s",((CAtom*)g_oaAtoms[o->m_waDecompTypeRefList[z3]])->m_sName,((CAtom*)g_oaAtoms[o->m_waDecompTypeObsList[z2]])->m_sName);
							buf.sprintf("%s-%s",(const char*)((CAtom*)g_oaAtoms[o->m_waDecompTypeRefList[z3]])->m_sName,(const char*)((CAtom*)g_oaAtoms[o->m_waDecompTypeObsList[z2]])->m_sName);
							o->m_pRDF[z0]->m_pRDF->SetLabelMulti(z*o->m_waDecompTypeObsList.GetSize()+z2,buf);
						}
					}
				} else o->m_pRDF[z0]->m_pRDF->Create();

				if (o->m_bSecondShowMol && (z0 == 1))
				{
					try { o->m_pRDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMol2Count]; } catch(...) { o->m_pRDF[z0]->m_faData = NULL; }
					if (o->m_pRDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pRDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMol2Count]; } catch(...) { o->m_pRDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pRDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				} else
				{
					try { o->m_pRDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMolCount]; } catch(...) { o->m_pRDF[z0]->m_faData = NULL; }
					if (o->m_pRDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pRDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMolCount]; } catch(...) { o->m_pRDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pRDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				}
				if (o->m_bTimeDev)
				{
					if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
					{
						try { o->m_pRDF[z0]->m_fDist = new FILE*[1]; } catch(...) { o->m_pRDF[z0]->m_fDist = NULL; }
						if (o->m_pRDF[z0]->m_fDist == NULL) NewException((double)sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
//						sprintf(buf,"rdf_timedev_%s%s.csv",o->m_pRDF[z0]->m_sName,multibuf);
						buf.sprintf("rdf_timedev_%s%s.csv",o->m_pRDF[z0]->m_sName,multibuf);
						o->m_pRDF[z0]->m_fDist[0] = OpenFileWrite(buf,true);
					} else
					{
						try { o->m_pRDF[z0]->m_fDist = new FILE*[o->m_waSaveRefList.GetSize()]; } catch(...) { o->m_pRDF[z0]->m_fDist = NULL; }
						if (o->m_pRDF[z0]->m_fDist == NULL) NewException((double)o->m_waSaveRefList.GetSize()*sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
//							sprintf(buf,"rdf_timedev_%s_ref%d%s.csv",o->m_pRDF[z0]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							buf.sprintf("rdf_timedev_%s_ref%d%s.csv",o->m_pRDF[z0]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							o->m_pRDF[z0]->m_fDist[z2] = OpenFileWrite(buf,true);
						}
					}
					if (o->m_bCombinedPlot)
					{
//						mprintf("Combined: Reflist %d, Showlist %d, Combinations %d, Steps %d.\n",o->m_waSaveRefList.GetSize(),o->m_waSaveShowList.GetSize(),o->m_pRDF[z0]->m_iCombinations,g_iTrajSteps);

						try { o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot = new CGrace(); } catch(...) { o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot = NULL; }
						if (o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->SetTitle("Combined distance time development/histogram");
						o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->SetSubTitle(o->m_pRDF[z0]->m_sShortName);
				//		smfix = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[0]];
				//		sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex[0]];
				//		o->m_pRDF[z0]->BuildAtomList(smfix,sm,sm,&templa);
						mprintf("    Trying to reserve %s of memory for combined plot...\n",FormatBytes((double)sizeof(double)*o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pRDF[z0]->m_iCombinations*g_iTrajSteps/g_iStride));
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
							for (z3=0;z3<o->m_waSaveShowList.GetSize();z3++)
							{
			/*					for (z4=0;z4<templa.GetSize()/2;z4++)
								{
									o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->AddDataset();
									sprintf(buf,"%s[%d] %s%d - %s[%d] %s%d",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,o->m_waSaveRefList[z2]+1,((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[g_waAtomElement[templa[z4*2]]]])->m_sName,g_waAtomMolNumber[templa[z4*2]]+1,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName,o->m_waSaveShowList[z3]+1,((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_baAtomIndex[g_waAtomElement[templa[z4*2+1]]]])->m_sName,g_waAtomMolNumber[templa[z4*2+1]]+1);
									o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->SetDatasetName(buf);
									if (g_iTrajSteps != -1)
										o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->LastDataset()->m_faValues.SetMaxSize(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize() *g_iTrajSteps);
									o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->LastDataset()->m_faValues.SetGrow(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize() *100);
									if (o->m_bCombinedGreyMode)
									{
										ti = o->m_iCombinedGreyMin + ((z4+z3*o->m_pRDF[z0]->m_iCombinations+z2*o->m_pRDF[z0]->m_iCombinations*o->m_waSaveShowList.GetSize())%o->m_iCombinedGreyShades)*(o->m_iCombinedGreyMax-o->m_iCombinedGreyMin)/o->m_iCombinedGreyShades;
										o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->LastDataset()->m_iLineColor = ti*0x10000 + ti*0x100 + ti;
									}
								}*/

								ti2 = 0;
								for (z4=0;z4<o->m_pRDF[z0]->m_oaVectors.GetSize()/2;z4++)
								{
									g1 = (CAtomGroup*)o->m_pRDF[z0]->m_oaVectors[z4*2];
									for (z1t=0;z1t<g1->m_baAtomType.GetSize();z1t++)
									{
										a1 = (CxIntArray*)g1->m_oaAtoms[z1t];
										for (z1a=0;z1a<a1->GetSize();z1a++)
										{
											g2 = (CAtomGroup*)o->m_pRDF[z0]->m_oaVectors[z4*2+1];
											for (z2t=0;z2t<g2->m_baAtomType.GetSize();z2t++)
											{
												a2 = (CxIntArray*)g2->m_oaAtoms[z2t];
												for (z2a=0;z2a<a2->GetSize();z2a++)
												{
													o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->AddDataset();
													if (o->m_bOthers)
//														sprintf(buf,"%s[%d] %s%d - %s[%d] %s%d",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,o->m_waSaveRefList[z2]+1,((CAtom*)g_oaAtoms[g1->m_baRealAtomType[z1t]])->m_sName,a1->GetAt(z1a)+1,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName,o->m_waSaveShowList[z3]+1,((CAtom*)g_oaAtoms[g2->m_baRealAtomType[z2t]])->m_sName,a2->GetAt(z2a)+1);
														buf.sprintf("%s[%d] %s%d - %s[%d] %s%d",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,o->m_waSaveRefList[z2]+1,(const char*)((CAtom*)g_oaAtoms[g1->m_baRealAtomType[z1t]])->m_sName,a1->GetAt(z1a)+1,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName,o->m_waSaveShowList[z3]+1,(const char*)((CAtom*)g_oaAtoms[g2->m_baRealAtomType[z2t]])->m_sName,a2->GetAt(z2a)+1);
													else
//														sprintf(buf,"%s[%d] %s%d - %s%d",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,o->m_waSaveRefList[z2]+1,((CAtom*)g_oaAtoms[g1->m_baRealAtomType[z1t]])->m_sName,a1->GetAt(z1a)+1,((CAtom*)g_oaAtoms[g2->m_baRealAtomType[z2t]])->m_sName,a2->GetAt(z2a)+1);
														buf.sprintf("%s[%d] %s%d - %s%d",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,o->m_waSaveRefList[z2]+1,(const char*)((CAtom*)g_oaAtoms[g1->m_baRealAtomType[z1t]])->m_sName,a1->GetAt(z1a)+1,(const char*)((CAtom*)g_oaAtoms[g2->m_baRealAtomType[z2t]])->m_sName,a2->GetAt(z2a)+1);

													o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->SetDatasetName(buf);
													if (g_iTrajSteps != -1)
														o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->LastDataset()->m_faValues.SetMaxSize(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()/**o->m_pRDF[z0]->m_iCombinations*/ *g_iTrajSteps/g_iStride);
													o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->LastDataset()->m_faValues.SetGrow(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()/**o->m_pRDF[z0]->m_iCombinations*/ *100);
													if (o->m_bCombinedGreyMode)
													{
														ti = o->m_iCombinedGreyMin + ((ti2+z3*o->m_pRDF[z0]->m_iCombinations+z2*o->m_pRDF[z0]->m_iCombinations*o->m_waSaveShowList.GetSize())%o->m_iCombinedGreyShades)*(o->m_iCombinedGreyMax-o->m_iCombinedGreyMin)/o->m_iCombinedGreyShades;
														o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->SetSetLineColor((unsigned char)ti,(unsigned char)ti,(unsigned char)ti);
													}
													ti2++;
/*													if (o->m_pRDF[z0]->m_iRefOrSec[0])
														vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
															else vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g1->m_baAtomType[z1t]])->GetAt(a1->GetAt(z1a)));
													if (o->m_pRDF[z0]->m_iRefOrSec[1])
														vec->Add(((CxIntArray*)obs->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
															else vec->Add(((CxIntArray*)ref->m_oaAtomOffset[g2->m_baAtomType[z2t]])->GetAt(a2->GetAt(z2a)));
*/										//			mprintf("Vector z=%d, z1t=%d, z1a=%d, z2t=%d, z2a=%d.\n",z,z1t,z1a,z2t,z2a);
												}
											}
										}
									}
								}

							}
						}
						o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->SetLabelX("Time / ps; g(r)");
						o->m_pRDF[z0]->m_pRDF->m_pCombinedPlot->SetLabelY("Distance / pm");
					}
				} // END IF TIMEDEV

				if (o->m_bTimeDiff)
					o->CreateTimeDiff(o->m_pRDF[z0]->m_pRDF,o->m_pRDF[z0]->m_iCombinations);

				if (g_bDeriv)
					o->m_pRDF[z0]->InitDeriv();

			} // END IF RDF

			if (g_bADF && (o->m_pADF[z0] != NULL))
			{
				mprintf("  Creating ADF...\n");
				o->m_pADF[z0]->m_pADF->m_fMinVal = o->m_pADF[z0]->m_fMinAngle;
				o->m_pADF[z0]->m_pADF->m_fMaxVal = o->m_pADF[z0]->m_fMaxAngle;
				o->m_pADF[z0]->m_pADF->m_iResolution = o->m_pADF[z0]->m_iResolution;
				o->m_pADF[z0]->m_pADF->m_iHistogramRes = o->m_pADF[z0]->m_iHistogramRes;
				if (o->m_pADF[z0]->m_bCosine)
					o->m_pADF[z0]->m_pADF->SetLabelX("Cos(angle)");
						else o->m_pADF[z0]->m_pADF->SetLabelX("Angle (degree)");
				o->m_pADF[z0]->m_pADF->SetLabelY("Occurrence");
				if (o->m_bObsCertain && o->m_bDecompDist)
					o->m_pADF[z0]->m_pADF->CreateMulti(o->m_waObsRefList.GetSize()*o->m_waObsShowList.GetSize());
						else o->m_pADF[z0]->m_pADF->Create();
				if (o->m_bSecondShowMol && (z0 == 1))
				{
					try { o->m_pADF[z0]->m_faData = new CxDoubleArray[o->m_iShowMol2Count]; } catch(...) { o->m_pADF[z0]->m_faData = NULL; }
					if (o->m_pADF[z0]->m_faData == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pADF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMol2Count]; } catch(...) { o->m_pADF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pADF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				} else
				{
					try { o->m_pADF[z0]->m_faData = new CxDoubleArray[o->m_iShowMolCount]; } catch(...) { o->m_pADF[z0]->m_faData = NULL; }
					if (o->m_pADF[z0]->m_faData == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pADF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMolCount]; } catch(...) { o->m_pADF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pADF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				}
				if (o->m_bTimeDev)
				{
					if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
					{
						try { o->m_pADF[z0]->m_fAngle = new FILE*[1]; } catch(...) { o->m_pADF[z0]->m_fAngle = NULL; }
						if (o->m_pADF[z0]->m_fAngle == NULL) NewException((double)sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
//						sprintf(buf,"adf_timedev_%s%s.csv",o->m_pADF[z0]->m_sName,multibuf);
						buf.sprintf("adf_timedev_%s%s.csv",o->m_pADF[z0]->m_sName,multibuf);
						o->m_pADF[z0]->m_fAngle[0] = OpenFileWrite(buf,true);
					} else
					{
						try { o->m_pADF[z0]->m_fAngle = new FILE*[o->m_waSaveRefList.GetSize()]; } catch(...) { o->m_pADF[z0]->m_fAngle = NULL; }
						if (o->m_pADF[z0]->m_fAngle == NULL) NewException((double)o->m_waSaveRefList.GetSize()*sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
//							sprintf(buf,"adf_timedev_%s_ref%d%s.csv",o->m_pADF[z0]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							buf.sprintf("adf_timedev_%s_ref%d%s.csv",o->m_pADF[z0]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							o->m_pADF[z0]->m_fAngle[z2] = OpenFileWrite(buf,true);
						}
					}
					if (o->m_bCombinedPlot)
					{
						try { o->m_pADF[z0]->m_pADF->m_pCombinedPlot = new CGrace(); } catch(...) { o->m_pADF[z0]->m_pADF->m_pCombinedPlot = NULL; }
						if (o->m_pADF[z0]->m_pADF->m_pCombinedPlot == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						o->m_pADF[z0]->m_pADF->m_pCombinedPlot->SetTitle("Combined angle time development/histogram");
						o->m_pADF[z0]->m_pADF->m_pCombinedPlot->SetSubTitle(o->m_pADF[z0]->m_sShortName);
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
							for (z3=0;z3<o->m_waSaveShowList.GetSize();z3++)
							{
								for (z4=0;z4<o->m_pADF[z0]->m_iCombinations;z4++)
								{
									o->m_pADF[z0]->m_pADF->m_pCombinedPlot->AddDataset();
									if (g_iTrajSteps != -1)
										o->m_pADF[z0]->m_pADF->m_pCombinedPlot->LastDataset()->m_faValues.SetMaxSize(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pADF[z0]->m_iCombinations*g_iTrajSteps);
									o->m_pADF[z0]->m_pADF->m_pCombinedPlot->LastDataset()->m_faValues.SetGrow(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pADF[z0]->m_iCombinations*100);
									if (o->m_bCombinedGreyMode)
									{
										ti = o->m_iCombinedGreyMin + ((z4+z3*o->m_pADF[z0]->m_iCombinations+z2*o->m_pADF[z0]->m_iCombinations*o->m_waSaveShowList.GetSize())%o->m_iCombinedGreyShades)*(o->m_iCombinedGreyMax-o->m_iCombinedGreyMin)/o->m_iCombinedGreyShades;
										o->m_pADF[z0]->m_pADF->m_pCombinedPlot->SetSetLineColor((unsigned char)ti,(unsigned char)ti,(unsigned char)ti);
									}
					//				o->m_pADF[z0]->m_pADF->m_pCombinedPlot->SetSetLineWidth((z2*o->m_waSaveRefList.GetSize()+z3)*o->m_pADF[z0]->m_iCombinations+z4,2.0f);
								}
							}
						}
						o->m_pADF[z0]->m_pADF->m_pCombinedPlot->SetLabelX("Time / ps ; ADF(r)");
						o->m_pADF[z0]->m_pADF->m_pCombinedPlot->SetLabelY("Angle / Degree");
					}
				} // END IF TIMEDEV

				if (o->m_bTimeDiff)
					o->CreateTimeDiff(o->m_pADF[z0]->m_pADF,o->m_pADF[z0]->m_iCombinations);

				if (g_bDeriv)
					o->m_pADF[z0]->InitDeriv();
			} // END IF ADF

			if (g_bDipDF && (o->m_pDipDF[z0] != NULL))
			{
				mprintf("  Creating DipDF...\n");
				o->m_pDipDF[z0]->m_pDipoleDF->m_fMinVal = o->m_pDipDF[z0]->m_fDipoleMin;
				o->m_pDipDF[z0]->m_pDipoleDF->m_fMaxVal = o->m_pDipDF[z0]->m_fDipoleMax;
				o->m_pDipDF[z0]->m_pDipoleDF->m_iResolution = o->m_pDipDF[z0]->m_iResolution;
				o->m_pDipDF[z0]->m_pDipoleDF->m_iHistogramRes = o->m_pDipDF[z0]->m_iHistogramRes;
				o->m_pDipDF[z0]->m_pDipoleDF->SetLabelX("Dipole moment (Debye)");
				o->m_pDipDF[z0]->m_pDipoleDF->SetLabelY("Occurrence");
				if (o->m_bObsCertain && o->m_bDecompDist)
					o->m_pDipDF[z0]->m_pDipoleDF->CreateMulti(o->m_waObsRefList.GetSize()*o->m_waObsShowList.GetSize());
						else o->m_pDipDF[z0]->m_pDipoleDF->Create();
				if (o->m_bSecondShowMol && (z0 == 1))
				{
					try { o->m_pDipDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMol2Count]; } catch(...) { o->m_pDipDF[z0]->m_faData = NULL; }
					if (o->m_pDipDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pDipDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMol2Count]; } catch(...) { o->m_pDipDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pDipDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				} else
				{
					try { o->m_pDipDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMolCount]; } catch(...) { o->m_pDipDF[z0]->m_faData = NULL; }
					if (o->m_pDipDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pDipDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMolCount]; } catch(...) { o->m_pDipDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pDipDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				}
				if (o->m_bTimeDev)
				{
					if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
					{
						try { o->m_pDipDF[z0]->m_fDipole = new FILE*[1]; } catch(...) { o->m_pDipDF[z0]->m_fDipole = NULL; }
						if (o->m_pDipDF[z0]->m_fDipole == NULL) NewException((double)sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
//						sprintf(buf,"dipole_timedev_%s%s.csv",o->m_pDipDF[z0]->m_sName,multibuf);
						buf.sprintf("dipole_timedev_%s%s.csv",o->m_pDipDF[z0]->m_sName,multibuf);
						o->m_pDipDF[z0]->m_fDipole[0] = OpenFileWrite(buf,true);
					} else
					{
						try { o->m_pDipDF[z0]->m_fDipole = new FILE*[o->m_waSaveRefList.GetSize()]; } catch(...) { o->m_pDipDF[z0]->m_fDipole = NULL; }
						if (o->m_pDipDF[z0]->m_fDipole == NULL) NewException((double)o->m_waSaveRefList.GetSize()*sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
//							sprintf(buf,"dipole_timedev_%s_ref%d%s.csv",o->m_pDipDF[z0]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							buf.sprintf("dipole_timedev_%s_ref%d%s.csv",o->m_pDipDF[z0]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							o->m_pDipDF[z0]->m_fDipole[z2] = OpenFileWrite(buf,true);
						}
					}
					if (o->m_bCombinedPlot)
					{
						try { o->m_pDipDF[z0]->m_pDipoleDF->m_pCombinedPlot = new CGrace(); } catch(...) { o->m_pDipDF[z0]->m_pDipoleDF->m_pCombinedPlot = NULL; }
						if (o->m_pDipDF[z0]->m_pDipoleDF->m_pCombinedPlot == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						o->m_pDipDF[z0]->m_pDipoleDF->m_pCombinedPlot->SetTitle("Combined dipole moment time development/histogram");
						o->m_pDipDF[z0]->m_pDipoleDF->m_pCombinedPlot->SetSubTitle(o->m_pDipDF[z0]->m_sShortName);
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
							for (z3=0;z3<o->m_waSaveShowList.GetSize();z3++)
							{
								o->m_pDipDF[z0]->m_pDipoleDF->m_pCombinedPlot->AddDataset();
								if (g_iTrajSteps != -1)
									o->m_pDipDF[z0]->m_pDipoleDF->m_pCombinedPlot->LastDataset()->m_faValues.SetMaxSize(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*g_iTrajSteps);
								o->m_pDipDF[z0]->m_pDipoleDF->m_pCombinedPlot->LastDataset()->m_faValues.SetGrow(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*100);
								if (o->m_bCombinedGreyMode)
								{
									ti = o->m_iCombinedGreyMin + ((z3+z2*o->m_waSaveShowList.GetSize())%o->m_iCombinedGreyShades)*(o->m_iCombinedGreyMax-o->m_iCombinedGreyMin)/o->m_iCombinedGreyShades;
									o->m_pDipDF[z0]->m_pDipoleDF->m_pCombinedPlot->SetSetLineColor((unsigned char)ti,(unsigned char)ti,(unsigned char)ti);
								}
							}
						o->m_pDipDF[z0]->m_pDipoleDF->m_pCombinedPlot->SetLabelX("Time / ps; DipDF(r)");
						o->m_pDipDF[z0]->m_pDipoleDF->m_pCombinedPlot->SetLabelY("Dipole moment / Debye");
					}
				} // END IF TIMEDEV

				if (o->m_bTimeDiff)
					o->CreateTimeDiff(o->m_pDipDF[z0]->m_pDipoleDF,1);

				if (g_bDeriv)
					o->m_pDipDF[z0]->InitDeriv();
			} // END IF DIPOLE

			if (g_bVDF && (o->m_pVDF[z0] != NULL))
			{
				mprintf("  Creating VDF...\n");
				o->m_pVDF[z0]->m_pVDF->m_fMinVal = o->m_pVDF[z0]->m_fMinSpeed;
				o->m_pVDF[z0]->m_pVDF->m_fMaxVal = o->m_pVDF[z0]->m_fMaxSpeed;
				o->m_pVDF[z0]->m_pVDF->m_iResolution = o->m_pVDF[z0]->m_iResolution;
				o->m_pVDF[z0]->m_pVDF->m_iHistogramRes = o->m_pVDF[z0]->m_iHistogramRes;
				o->m_pVDF[z0]->m_pVDF->SetLabelX("Velocity / pm ps^-1");
				o->m_pVDF[z0]->m_pVDF->SetLabelY("Occurrence");
				if (o->m_bObsCertain && o->m_bDecompDist)
					o->m_pVDF[z0]->m_pVDF->CreateMulti(o->m_waObsRefList.GetSize()*o->m_waObsShowList.GetSize());
						else o->m_pVDF[z0]->m_pVDF->Create();
				if (o->m_pVDF[z0]->m_bSplitCart) {
					for (z2=0;z2<6;z2++) {
						o->m_pVDF[z0]->m_pVDFSplit[z2] = new CDF();
						o->m_pVDF[z0]->m_pVDFSplit[z2]->m_fMinVal = o->m_pVDF[z0]->m_fMinSpeed;
						o->m_pVDF[z0]->m_pVDFSplit[z2]->m_fMaxVal = o->m_pVDF[z0]->m_fMaxSpeed;
						o->m_pVDF[z0]->m_pVDFSplit[z2]->m_iResolution = o->m_pVDF[z0]->m_iResolution;
						o->m_pVDF[z0]->m_pVDFSplit[z2]->m_iHistogramRes = o->m_pVDF[z0]->m_iHistogramRes;
						switch(z2) {
							case 0: o->m_pVDF[z0]->m_pVDFSplit[z2]->SetLabelX("Velocity X Projection / pm ps^-1"); break;
							case 1: o->m_pVDF[z0]->m_pVDFSplit[z2]->SetLabelX("Velocity Y Projection / pm ps^-1"); break;
							case 2: o->m_pVDF[z0]->m_pVDFSplit[z2]->SetLabelX("Velocity Z Projection / pm ps^-1"); break;
							case 3: o->m_pVDF[z0]->m_pVDFSplit[z2]->SetLabelX("Velocity XY Projection / pm ps^-1"); break;
							case 4: o->m_pVDF[z0]->m_pVDFSplit[z2]->SetLabelX("Velocity XZ Projection / pm ps^-1"); break;
							case 5: o->m_pVDF[z0]->m_pVDFSplit[z2]->SetLabelX("Velocity YZ Projection / pm ps^-1"); break;
						}
						o->m_pVDF[z0]->m_pVDFSplit[z2]->SetLabelY("Occurrence");
						o->m_pVDF[z0]->m_pVDFSplit[z2]->Create();
					}
				}
				if (o->m_bSecondShowMol && (z0 == 1))
				{
					try { o->m_pVDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMol2Count]; } catch(...) { o->m_pVDF[z0]->m_faData = NULL; }
					if (o->m_pVDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pVDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMol2Count]; } catch(...) { o->m_pVDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pVDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				} else
				{
					try { o->m_pVDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMolCount]; } catch(...) { o->m_pVDF[z0]->m_faData = NULL; }
					if (o->m_pVDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pVDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMolCount]; } catch(...) { o->m_pVDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pVDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				}
				if (o->m_bTimeDev)
				{
					if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
					{
						try { o->m_pVDF[z0]->m_fSpeed = new FILE*[1]; } catch(...) { o->m_pVDF[z0]->m_fSpeed = NULL; }
						if (o->m_pVDF[z0]->m_fSpeed == NULL) NewException((double)sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
//						sprintf(buf,"vdf_timedev_%s%s.csv",o->m_pVDF[z0]->m_sName,multibuf);
						buf.sprintf("vdf_timedev_%s%s.csv",o->m_pVDF[z0]->m_sName,multibuf);
						o->m_pVDF[z0]->m_fSpeed[0] = OpenFileWrite(buf,true);
					} else
					{
						try { o->m_pVDF[z0]->m_fSpeed = new FILE*[o->m_waSaveRefList.GetSize()]; } catch(...) { o->m_pVDF[z0]->m_fSpeed = NULL; }
						if (o->m_pVDF[z0]->m_fSpeed == NULL) NewException((double)o->m_waSaveRefList.GetSize()*sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
//							sprintf(buf,"vdf_timedev_%s_ref%d%s.csv",o->m_pVDF[z0]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							buf.sprintf("vdf_timedev_%s_ref%d%s.csv",o->m_pVDF[z0]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							o->m_pVDF[z0]->m_fSpeed[z2] = OpenFileWrite(buf,true);
						}
					}
					if (o->m_bCombinedPlot)
					{
						try { o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot = new CGrace(); } catch(...) { o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot = NULL; }
						if (o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot->SetTitle("Combined velocity time development/histogram");
						o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot->SetSubTitle(o->m_pVDF[z0]->m_sShortName);
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
							for (z3=0;z3<o->m_waSaveShowList.GetSize();z3++)
								for (z4=0;z4<o->m_pVDF[z0]->m_iCombinations;z4++)
								{
									o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot->AddDataset();
									if (g_iTrajSteps != -1)
										o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot->LastDataset()->m_faValues.SetMaxSize(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pVDF[z0]->m_iCombinations*g_iTrajSteps);
									o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot->LastDataset()->m_faValues.SetGrow(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pVDF[z0]->m_iCombinations*100);
									if (o->m_bCombinedGreyMode)
									{
										ti = o->m_iCombinedGreyMin + ((z4+z3*o->m_pVDF[z0]->m_iCombinations+z2*o->m_pVDF[z0]->m_iCombinations*o->m_waSaveShowList.GetSize())%o->m_iCombinedGreyShades)*(o->m_iCombinedGreyMax-o->m_iCombinedGreyMin)/o->m_iCombinedGreyShades;
										o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot->SetSetLineColor((unsigned char)ti,(unsigned char)ti,(unsigned char)ti);
									}
//									o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot->SetSetLineWidth(z2*o->m_waSaveRefList.GetSize()+z3,3.0f);
								}
						o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot->SetLabelX("Time / ps ; VDF(r)");
						o->m_pVDF[z0]->m_pVDF->m_pCombinedPlot->SetLabelY("Velocity / pm ps^-1");
					}
				} // END IF TIMEDEV

				if (o->m_bTimeDiff)
					o->CreateTimeDiff(o->m_pVDF[z0]->m_pVDF,o->m_pVDF[z0]->m_iCombinations);

				if (g_bDeriv)
					o->m_pVDF[z0]->InitDeriv();
			} // END IF VDF

			if (g_bDDF && (o->m_pDDF[z0] != NULL))
			{
				mprintf("  Creating DDF...\n");
				o->m_pDDF[z0]->m_pDDF->m_fMinVal = o->m_pDDF[z0]->m_fMinAngle;
				o->m_pDDF[z0]->m_pDDF->m_fMaxVal = o->m_pDDF[z0]->m_fMaxAngle;
				o->m_pDDF[z0]->m_pDDF->m_iResolution = o->m_pDDF[z0]->m_iResolution;
				o->m_pDDF[z0]->m_pDDF->m_iHistogramRes = o->m_pDDF[z0]->m_iHistogramRes;
				if (o->m_pDDF[z0]->m_bCosine)
					o->m_pDDF[z0]->m_pDDF->SetLabelX("Cos(Dihedral Angle)");
						else o->m_pDDF[z0]->m_pDDF->SetLabelX("Dihedral Angle (Degree)");
				o->m_pDDF[z0]->m_pDDF->SetLabelY("Occurrence");
				if (o->m_bObsCertain && o->m_bDecompDist)
					o->m_pDDF[z0]->m_pDDF->CreateMulti(o->m_waObsRefList.GetSize()*o->m_waObsShowList.GetSize());
						else o->m_pDDF[z0]->m_pDDF->Create();
				if (o->m_bSecondShowMol && (z0 == 1))
				{
					try { o->m_pDDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMol2Count]; } catch(...) { o->m_pDDF[z0]->m_faData = NULL; }
					if (o->m_pDDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pDDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMol2Count]; } catch(...) { o->m_pDDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pDDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				} else
				{
					try { o->m_pDDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMolCount]; } catch(...) { o->m_pDDF[z0]->m_faData = NULL; }
					if (o->m_pDDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pDDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMolCount]; } catch(...) { o->m_pDDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pDDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
						
					if (o->m_pDDF[z0]->m_bRotate)
					{
						for (z2=0;z2<o->m_iShowMolCount * ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() * o->m_pDDF[z0]->m_iCombinations;z2++)
						{
							o->m_pDDF[z0]->m_faLastData.Add(0);
							o->m_pDDF[z0]->m_laRotation.Add(0);
						}
					}
				}
				if (o->m_bTimeDev)
				{
					if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
					{
						try { o->m_pDDF[z0]->m_fAngle = new FILE*[1]; } catch(...) { o->m_pDDF[z0]->m_fAngle = NULL; }
						if (o->m_pDDF[z0]->m_fAngle == NULL) NewException((double)sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
//						sprintf(buf,"ddf_timedev_%s%s.csv",o->m_pDDF[z0]->m_sName,multibuf);
						buf.sprintf("ddf_timedev_%s%s.csv",o->m_pDDF[z0]->m_sName,multibuf);
						o->m_pDDF[z0]->m_fAngle[0] = OpenFileWrite(buf,true);
					} else
					{
						try { o->m_pDDF[z0]->m_fAngle = new FILE*[o->m_waSaveRefList.GetSize()]; } catch(...) { o->m_pDDF[z0]->m_fAngle = NULL; }
						if (o->m_pDDF[z0]->m_fAngle == NULL) NewException((double)o->m_waSaveRefList.GetSize()*sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
//							sprintf(buf,"ddf_timedev_%s_ref%d%s.csv",o->m_pDDF[z0]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							buf.sprintf("ddf_timedev_%s_ref%d%s.csv",o->m_pDDF[z0]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							o->m_pDDF[z0]->m_fAngle[z2] = OpenFileWrite(buf,true);
						}
					}
					if (o->m_bCombinedPlot)
					{
						try { o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot = new CGrace(); } catch(...) { o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot = NULL; }
						if (o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot->SetTitle("Combined dihedral time development/histogram");
						o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot->SetSubTitle(o->m_pDDF[z0]->m_sShortName);
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
							for (z3=0;z3<o->m_waSaveShowList.GetSize();z3++)
								for (z4=0;z4<o->m_pDDF[z0]->m_iCombinations;z4++)
								{
									o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot->AddDataset();
									if (g_iTrajSteps != -1)
										o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot->LastDataset()->m_faValues.SetMaxSize(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pDDF[z0]->m_iCombinations*g_iTrajSteps);
									o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot->LastDataset()->m_faValues.SetGrow(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pDDF[z0]->m_iCombinations*100);
									if (o->m_bCombinedGreyMode)
									{
										ti = o->m_iCombinedGreyMin + ((z4+z3*o->m_pDDF[z0]->m_iCombinations+z2*o->m_pDDF[z0]->m_iCombinations*o->m_waSaveShowList.GetSize())%o->m_iCombinedGreyShades)*(o->m_iCombinedGreyMax-o->m_iCombinedGreyMin)/o->m_iCombinedGreyShades;
										o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot->SetSetLineColor((unsigned char)ti,(unsigned char)ti,(unsigned char)ti);
									}
				//					o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot->SetSetLineWidth((z2*o->m_waSaveRefList.GetSize()+z3)*o->m_pDDF[z0]->m_iCombinations+z4,2.0f);
								}
						o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot->SetLabelX("Time / ps ; DDF(r)");
						o->m_pDDF[z0]->m_pDDF->m_pCombinedPlot->SetLabelY("Dihedral / Degree");
					}
				} // END IF TIMEDEV

				if (o->m_bTimeDiff)
					o->CreateTimeDiff(o->m_pDDF[z0]->m_pDDF,o->m_pDDF[z0]->m_iCombinations);

				if (g_bDeriv)
					o->m_pDDF[z0]->InitDeriv();
			} // END IF DDF

			if (g_bPlDF && (o->m_pPlDF[z0] != NULL))
			{
				mprintf("  Creating PlDF...\n");
				o->m_pPlDF[z0]->m_pPlDF->m_fMinVal = o->m_pPlDF[z0]->m_fMinDist;
				o->m_pPlDF[z0]->m_pPlDF->m_fMaxVal = o->m_pPlDF[z0]->m_fMaxDist;
				o->m_pPlDF[z0]->m_pPlDF->m_iResolution = o->m_pPlDF[z0]->m_iResolution;
				o->m_pPlDF[z0]->m_pPlDF->m_iHistogramRes = o->m_pPlDF[z0]->m_iHistogramRes;
				o->m_pPlDF[z0]->m_pPlDF->SetLabelX("Distance from plane / pm");
				o->m_pPlDF[z0]->m_pPlDF->SetLabelY("Occurrence");
				o->m_pPlDF[z0]->m_pPlDF->Create();
				if (o->m_bSecondShowMol && (z0 == 1))
				{
					try { o->m_pPlDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMol2Count]; } catch(...) { o->m_pPlDF[z0]->m_faData = NULL; }
					if (o->m_pPlDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pPlDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMol2Count]; } catch(...) { o->m_pPlDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pPlDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				} else
				{
					try { o->m_pPlDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMolCount]; } catch(...) { o->m_pPlDF[z0]->m_faData = NULL; }
					if (o->m_pPlDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pPlDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMolCount]; } catch(...) { o->m_pPlDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pPlDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				}
			} // END IF PlDF

			if (g_bLiDF && (o->m_pLiDF[z0] != NULL))
			{
				mprintf("  Creating LiDF...\n");
				o->m_pLiDF[z0]->m_pLiDF->m_fMinVal = o->m_pLiDF[z0]->m_fMinDist;
				o->m_pLiDF[z0]->m_pLiDF->m_fMaxVal = o->m_pLiDF[z0]->m_fMaxDist;
				o->m_pLiDF[z0]->m_pLiDF->m_iResolution = o->m_pLiDF[z0]->m_iResolution;
				o->m_pLiDF[z0]->m_pLiDF->m_iHistogramRes = o->m_pLiDF[z0]->m_iHistogramRes;
				o->m_pLiDF[z0]->m_pLiDF->SetLabelX("Distance from line / pm");
				o->m_pLiDF[z0]->m_pLiDF->SetLabelY("Occurrence");
				o->m_pLiDF[z0]->m_pLiDF->Create();
				if (o->m_bSecondShowMol && (z0 == 1))
				{
					try { o->m_pLiDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMol2Count]; } catch(...) { o->m_pLiDF[z0]->m_faData = NULL; }
					if (o->m_pLiDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pLiDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMol2Count]; } catch(...) { o->m_pLiDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pLiDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMol2Count*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				} else
				{
					try { o->m_pLiDF[z0]->m_faData = new CxDoubleArray[o->m_iShowMolCount]; } catch(...) { o->m_pLiDF[z0]->m_faData = NULL; }
					if (o->m_pLiDF[z0]->m_faData == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
					{
						try { o->m_pLiDF[z0]->m_baDataEnabled = new CxByteArray[o->m_iShowMolCount]; } catch(...) { o->m_pLiDF[z0]->m_baDataEnabled = NULL; }
						if (o->m_pLiDF[z0]->m_baDataEnabled == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxByteArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					}
				}
			} // END IF LiDF

		} // END FOR z0

		if (g_bCDF)
		{
			mprintf("  Creating CDF...\n");
			if (g_iCDFChannels == 2)
			{
				try { o->m_pCDF->m_p2DF = new C2DF(); } catch(...) { o->m_pCDF->m_p2DF = NULL; }
				if (o->m_pCDF->m_p2DF == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
//				buf[0] = 0;
				buf = "";
				for (z2=0;z2<2;z2++)
				{
					o->m_pCDF->m_p2DF->m_iRes[z2] = o->m_pCDF->m_iResolution[z2];
					switch(g_iObsChannel[z2])
					{
						case 0:
//							strcat(buf,"_rdf");
//							strcat(buf,o->m_pRDF[z2]->m_sShortName);
							buf.strcat("_rdf");
							buf.strcat(o->m_pRDF[z2]->m_sShortName);
							o->m_pCDF->m_p2DF->m_fMinVal[z2] = o->m_pRDF[z2]->m_fMinDist;
							o->m_pCDF->m_p2DF->m_fMaxVal[z2] = o->m_pRDF[z2]->m_fMaxDist;
							o->m_pCDF->m_p2DF->m_pChannels[z2] = o->m_pRDF[z2]->m_pRDF;
							break;

						case 1:
//							strcat(buf,"_adf");
//							strcat(buf,o->m_pADF[z2]->m_sShortName);
							buf.strcat("_adf");
							buf.strcat(o->m_pADF[z2]->m_sShortName);
							o->m_pCDF->m_p2DF->m_fMinVal[z2] = o->m_pADF[z2]->m_fMinAngle;
							o->m_pCDF->m_p2DF->m_fMaxVal[z2] = o->m_pADF[z2]->m_fMaxAngle;
							o->m_pCDF->m_p2DF->m_pChannels[z2] = o->m_pADF[z2]->m_pADF;
							break;

						case 2:
//							strcat(buf,"_ddf");
//							strcat(buf,o->m_pDDF[z2]->m_sShortName);
							buf.strcat("_ddf");
							buf.strcat(o->m_pDDF[z2]->m_sShortName);
							o->m_pCDF->m_p2DF->m_fMinVal[z2] = o->m_pDDF[z2]->m_fMinAngle;
							o->m_pCDF->m_p2DF->m_fMaxVal[z2] = o->m_pDDF[z2]->m_fMaxAngle;
							o->m_pCDF->m_p2DF->m_pChannels[z2] = o->m_pDDF[z2]->m_pDDF;
							break;

						case 3:
//							strcat(buf,"_dipole");
//							strcat(buf,o->m_pDipDF[z2]->m_sShortName);
							buf.strcat("_dipole");
							buf.strcat(o->m_pDipDF[z2]->m_sShortName);
							o->m_pCDF->m_p2DF->m_fMinVal[z2] = o->m_pDipDF[z2]->m_fDipoleMin;
							o->m_pCDF->m_p2DF->m_fMaxVal[z2] = o->m_pDipDF[z2]->m_fDipoleMax;
							o->m_pCDF->m_p2DF->m_pChannels[z2] = o->m_pDipDF[z2]->m_pDipoleDF;
							break;

						case 4:
//							strcat(buf,"_vdf");
//							strcat(buf,o->m_pVDF[z2]->m_sShortName);
							buf.strcat("_vdf");
							buf.strcat(o->m_pVDF[z2]->m_sShortName);
							o->m_pCDF->m_p2DF->m_fMinVal[z2] = o->m_pVDF[z2]->m_fMinSpeed;
							o->m_pCDF->m_p2DF->m_fMaxVal[z2] = o->m_pVDF[z2]->m_fMaxSpeed;
							o->m_pCDF->m_p2DF->m_pChannels[z2] = o->m_pVDF[z2]->m_pVDF;
							break;

						case 5:
//							strcat(buf,"_pldf");
//							strcat(buf,o->m_pPlDF[z2]->m_sShortName);
							buf.strcat("_pldf");
							buf.strcat(o->m_pPlDF[z2]->m_sShortName);
							o->m_pCDF->m_p2DF->m_fMinVal[z2] = o->m_pPlDF[z2]->m_fMinDist;
							o->m_pCDF->m_p2DF->m_fMaxVal[z2] = o->m_pPlDF[z2]->m_fMaxDist;
							o->m_pCDF->m_p2DF->m_pChannels[z2] = o->m_pPlDF[z2]->m_pPlDF;
							break;

						case 6:
//							strcat(buf,"_lidf");
//							strcat(buf,o->m_pLiDF[z2]->m_sShortName);
							buf.strcat("_lidf");
							buf.strcat(o->m_pLiDF[z2]->m_sShortName);
							o->m_pCDF->m_p2DF->m_fMinVal[z2] = o->m_pLiDF[z2]->m_fMinDist;
							o->m_pCDF->m_p2DF->m_fMaxVal[z2] = o->m_pLiDF[z2]->m_fMaxDist;
							o->m_pCDF->m_p2DF->m_pChannels[z2] = o->m_pLiDF[z2]->m_pLiDF;
							break;
					}
				}
				o->m_pCDF->m_p2DF->m_iHistogramRes = o->m_pCDF->m_iHistogramRes;
				o->m_pCDF->m_p2DF->Create();
				switch(g_iObsChannel[0])
				{
					case 0:
//						sprintf(buf2,"%s Distance [pm]",o->m_pRDF[0]->m_sShortName);
						buf2.sprintf("%s Distance / pm",o->m_pRDF[0]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelX(buf2);
						break;
					case 1:
//						sprintf(buf2,"%s Angle [Degree]",o->m_pADF[0]->m_sShortName);
						buf2.sprintf("%s Angle / Degree",o->m_pADF[0]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelX(buf2);
						break;
					case 2:
//						sprintf(buf2,"%s Dihedral [Degree]",o->m_pDDF[0]->m_sShortName);
						buf2.sprintf("%s Dihedral / Degree",o->m_pDDF[0]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelX(buf2);
						break;
					case 3:
//						sprintf(buf2,"%s Dipole moment [Debye]",o->m_pDipDF[0]->m_sShortName);
						buf2.sprintf("%s Dipole moment / Debye",o->m_pDipDF[0]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelX(buf2);
						break;
					case 4:
//						sprintf(buf2,"%s Velocity [pm/ps]",o->m_pVDF[0]->m_sShortName);
						buf2.sprintf("%s Velocity / pm ps^-1",o->m_pVDF[0]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelX(buf2);
						break;
					case 5:
//						sprintf(buf2,"%s Distance from Plane [pm]",o->m_pPlDF[0]->m_sShortName);
						buf2.sprintf("%s Distance from Plane / pm",o->m_pPlDF[0]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelX(buf2);
						break;
					case 6:
//						sprintf(buf2,"%s Distance from Line [pm]",o->m_pLiDF[0]->m_sShortName);
						buf2.sprintf("%s Distance from Line / pm",o->m_pLiDF[0]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelX(buf2);
						break;
				}
				switch(g_iObsChannel[1])
				{
					case 0:
//						sprintf(buf2,"%s Distance [pm]",o->m_pRDF[1]->m_sShortName);
						buf2.sprintf("%s Distance / pm",o->m_pRDF[1]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelY(buf2);
						break;
					case 1:
//						sprintf(buf2,"%s Angle [Degree]",o->m_pADF[1]->m_sShortName);
						buf2.sprintf("%s Angle / Degree",o->m_pADF[1]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelY(buf2);
						break;
					case 2:
//						sprintf(buf2,"%s Dihedral [Degree]",o->m_pDDF[1]->m_sShortName);
						buf2.sprintf("%s Dihedral / Degree",o->m_pDDF[1]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelY(buf2);
						break;
					case 3:
//						sprintf(buf2,"%s Dipole moment [Debye]",o->m_pDipDF[1]->m_sShortName);
						buf2.sprintf("%s Dipole moment / Debye",o->m_pDipDF[1]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelY(buf2);
						break;
					case 4:
//						sprintf(buf2,"%s Velocity [pm/ps]",o->m_pVDF[1]->m_sShortName);
						buf2.sprintf("%s Velocity / pm ps^-1",o->m_pVDF[1]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelY(buf2);
						break;
					case 5:
//						sprintf(buf2,"%s Distance from Plane [pm]",o->m_pPlDF[1]->m_sShortName);
						buf2.sprintf("%s Distance from Plane / pm",o->m_pPlDF[1]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelY(buf2);
						break;
					case 6:
//						sprintf(buf2,"%s Distance from Line [pm]",o->m_pLiDF[1]->m_sShortName);
						buf2.sprintf("%s Distance from Line / pm",o->m_pLiDF[1]->m_sLabelName);
						o->m_pCDF->m_p2DF->SetLabelY(buf2);
						break;
				}
			} // END IF CHANNELS == 2

			if (g_iCDFChannels == 3)
			{
				try { o->m_pCDF->m_p3DF = new C3DF<double>(); } catch(...) { o->m_pCDF->m_p3DF = NULL; }
				if (o->m_pCDF->m_p3DF == NULL) NewException((double)sizeof(C3DF<double>),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
//				buf[0] = 0;
				buf = "";

				for (z2=0;z2<g_iCDFChannels;z2++)
				{
					o->m_pCDF->m_p3DF->m_iRes[z2] = o->m_pCDF->m_iResolution[z2];
					switch(g_iObsChannel[z2])
					{
						case 0:
//							strcat(buf,"_rdf");
//							strcat(buf,o->m_pRDF[z2]->m_sShortName);
							buf.strcat("_rdf");
							buf.strcat(o->m_pRDF[z2]->m_sShortName);
							o->m_pCDF->m_p3DF->m_fMinVal[z2] = o->m_pRDF[z2]->m_fMinDist;
							o->m_pCDF->m_p3DF->m_fMaxVal[z2] = o->m_pRDF[z2]->m_fMaxDist;
							o->m_pCDF->m_p3DF->m_pChannels[z2] = o->m_pRDF[z2]->m_pRDF;
							break;
						case 1:
//							strcat(buf,"_adf");
//							strcat(buf,o->m_pADF[z2]->m_sShortName);
							buf.strcat("_adf");
							buf.strcat(o->m_pADF[z2]->m_sShortName);
							o->m_pCDF->m_p3DF->m_fMinVal[z2] = o->m_pADF[z2]->m_fMinAngle;
							o->m_pCDF->m_p3DF->m_fMaxVal[z2] = o->m_pADF[z2]->m_fMaxAngle;
							o->m_pCDF->m_p3DF->m_pChannels[z2] = o->m_pADF[z2]->m_pADF;
							break;
						case 2:
//							strcat(buf,"_ddf");
//							strcat(buf,o->m_pDDF[z2]->m_sShortName);
							buf.strcat("_ddf");
							buf.strcat(o->m_pDDF[z2]->m_sShortName);
							o->m_pCDF->m_p3DF->m_fMinVal[z2] = o->m_pDDF[z2]->m_fMinAngle;
							o->m_pCDF->m_p3DF->m_fMaxVal[z2] = o->m_pDDF[z2]->m_fMaxAngle;
							o->m_pCDF->m_p3DF->m_pChannels[z2] = o->m_pDDF[z2]->m_pDDF;
							break;
						case 3:
//							strcat(buf,"_dipole");
//							strcat(buf,o->m_pDipDF[z2]->m_sShortName);
							buf.strcat("_dipole");
							buf.strcat(o->m_pDipDF[z2]->m_sShortName);
							o->m_pCDF->m_p3DF->m_fMinVal[z2] = o->m_pDipDF[z2]->m_fDipoleMin;
							o->m_pCDF->m_p3DF->m_fMaxVal[z2] = o->m_pDipDF[z2]->m_fDipoleMax;
							o->m_pCDF->m_p3DF->m_pChannels[z2] = o->m_pDipDF[z2]->m_pDipoleDF;
							break;
						case 4:
//							strcat(buf,"_vdf");
//							strcat(buf,o->m_pVDF[z2]->m_sShortName);
							buf.strcat("_vdf");
							buf.strcat(o->m_pVDF[z2]->m_sShortName);
							o->m_pCDF->m_p3DF->m_fMinVal[z2] = o->m_pVDF[z2]->m_fMinSpeed;
							o->m_pCDF->m_p3DF->m_fMaxVal[z2] = o->m_pVDF[z2]->m_fMaxSpeed;
							o->m_pCDF->m_p3DF->m_pChannels[z2] = o->m_pVDF[z2]->m_pVDF;
							break;
						case 5:
//							strcat(buf,"_pldf");
//							strcat(buf,o->m_pPlDF[z2]->m_sShortName);
							buf.strcat("_pldf");
							buf.strcat(o->m_pPlDF[z2]->m_sShortName);
							o->m_pCDF->m_p3DF->m_fMinVal[z2] = o->m_pPlDF[z2]->m_fMinDist;
							o->m_pCDF->m_p3DF->m_fMaxVal[z2] = o->m_pPlDF[z2]->m_fMaxDist;
							o->m_pCDF->m_p3DF->m_pChannels[z2] = o->m_pPlDF[z2]->m_pPlDF;
							break;
						case 6:
//							strcat(buf,"_lidf");
//							strcat(buf,o->m_pLiDF[z2]->m_sShortName);
							buf.strcat("_lidf");
							buf.strcat(o->m_pLiDF[z2]->m_sShortName);
							o->m_pCDF->m_p3DF->m_fMinVal[z2] = o->m_pLiDF[z2]->m_fMinDist;
							o->m_pCDF->m_p3DF->m_fMaxVal[z2] = o->m_pLiDF[z2]->m_fMaxDist;
							o->m_pCDF->m_p3DF->m_pChannels[z2] = o->m_pLiDF[z2]->m_pLiDF;
							break;
					}
				}
				o->m_pCDF->m_p3DF->m_iHistogramRes = o->m_pCDF->m_iHistogramRes;
				o->m_pCDF->m_p3DF->Create();
				switch(g_iObsChannel[0])
				{
					case 0:
//						sprintf(buf2,"%s Distance [pm]",o->m_pRDF[0]->m_sShortName);
						buf2.sprintf("%s Distance / pm",o->m_pRDF[0]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelX(buf2);
						break;
					case 1:
//						sprintf(buf2,"%s Angle [Degree]",o->m_pADF[0]->m_sShortName);
						buf2.sprintf("%s Angle / Degree",o->m_pADF[0]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelX(buf2);
						break;
					case 2:
//						sprintf(buf2,"%s Dihedral [Degree]",o->m_pDDF[0]->m_sShortName);
						buf2.sprintf("%s Dihedral / Degree",o->m_pDDF[0]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelX(buf2);
						break;
					case 3:
//						sprintf(buf2,"%s Dipole moment [Debye]",o->m_pDipDF[0]->m_sShortName);
						buf2.sprintf("%s Dipole moment / Debye",o->m_pDipDF[0]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelX(buf2);
						break;
					case 4:
//						sprintf(buf2,"%s Velocity [pm/ps]",o->m_pVDF[0]->m_sShortName);
						buf2.sprintf("%s Velocity / pm ps^-1",o->m_pVDF[0]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelX(buf2);
						break;
					case 5:
//						sprintf(buf2,"%s Distance from Plane [pm]",o->m_pPlDF[0]->m_sShortName);
						buf2.sprintf("%s Distance from Plane / pm",o->m_pPlDF[0]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelX(buf2);
						break;
					case 6:
//						sprintf(buf2,"%s Distance from Line [pm]",o->m_pLiDF[0]->m_sShortName);
						buf2.sprintf("%s Distance from Line / pm",o->m_pLiDF[0]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelX(buf2);
						break;
				}
				switch(g_iObsChannel[1])
				{
					case 0:
//						sprintf(buf2,"%s Distance [pm]",o->m_pRDF[1]->m_sShortName);
						buf2.sprintf("%s Distance / pm",o->m_pRDF[1]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelY(buf2);
						break;
					case 1:
//						sprintf(buf2,"%s Angle [Degree]",o->m_pADF[1]->m_sShortName);
						buf2.sprintf("%s Angle / Degree",o->m_pADF[1]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelY(buf2);
						break;
					case 2:
//						sprintf(buf2,"%s Dihedral [Degree]",o->m_pDDF[1]->m_sShortName);
						buf2.sprintf("%s Dihedral / Degree",o->m_pDDF[1]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelY(buf2);
						break;
					case 3:
//						sprintf(buf2,"%s Dipole moment [Debye]",o->m_pDipDF[1]->m_sShortName);
						buf2.sprintf("%s Dipole moment / Debye",o->m_pDipDF[1]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelY(buf2);
						break;
					case 4:
//						sprintf(buf2,"%s Velocity [pm/ps]",o->m_pVDF[1]->m_sShortName);
						buf2.sprintf("%s Velocity / pm ps^-1",o->m_pVDF[1]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelY(buf2);
						break;
					case 5:
//						sprintf(buf2,"%s Distance from Plane [pm]",o->m_pPlDF[1]->m_sShortName);
						buf2.sprintf("%s Distance from Plane / pm",o->m_pPlDF[1]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelY(buf2);
						break;
					case 6:
//						sprintf(buf2,"%s Distance from Line [pm]",o->m_pLiDF[1]->m_sShortName);
						buf2.sprintf("%s Distance from Line / pm",o->m_pLiDF[1]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelY(buf2);
						break;
				}
				switch(g_iObsChannel[2])
				{
					case 0:
//						sprintf(buf2,"%s Distance [pm]",o->m_pRDF[2]->m_sShortName);
						buf2.sprintf("%s Distance / pm",o->m_pRDF[2]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelZ(buf2);
						break;
					case 1:
//						sprintf(buf2,"%s Angle [Degree]",o->m_pADF[2]->m_sShortName);
						buf2.sprintf("%s Angle / Degree",o->m_pADF[2]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelZ(buf2);
						break;
					case 2:
//						sprintf(buf2,"%s Dihedral [Degree]",o->m_pDDF[2]->m_sShortName);
						buf2.sprintf("%s Dihedral / Degree",o->m_pDDF[2]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelZ(buf2);
						break;
					case 3:
//						sprintf(buf2,"%s Dipole moment [Debye]",o->m_pDipDF[2]->m_sShortName);
						buf2.sprintf("%s Dipole moment / Debye",o->m_pDipDF[2]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelZ(buf2);
						break;
					case 4:
//						sprintf(buf2,"%s Velocity [pm/ps]",o->m_pVDF[2]->m_sShortName);
						buf2.sprintf("%s Velocity / pm ps^-1",o->m_pVDF[2]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelZ(buf2);
						break;
					case 5:
//						sprintf(buf2,"%s Distance from Plane [pm]",o->m_pPlDF[2]->m_sShortName);
						buf2.sprintf("%s Distance from Plane / pm",o->m_pPlDF[2]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelZ(buf2);
						break;
					case 6:
//						sprintf(buf2,"%s Distance from Line [pm]",o->m_pLiDF[2]->m_sShortName);
						buf2.sprintf("%s Distance from Line / pm",o->m_pLiDF[2]->m_sLabelName);
						o->m_pCDF->m_p3DF->SetLabelZ(buf2);
						break;
				}

				// Fuer jede C3DF noch die 3 C2DFs erzeugen
				for (z3=0;z3<3;z3++)
				{
					try { o->m_pCDF->m_p3DF->m_p2DF[z3] = new C2DF(); } catch(...) { o->m_pCDF->m_p3DF->m_p2DF[z3] = NULL; }
					if (o->m_pCDF->m_p3DF->m_p2DF[z3] == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
//					buf3[0] = 0;
					buf3 = "";
					switch(z3)
					{
						case 0:
							tia[0] = 0;
							tia[1] = 1;
							break;
						case 1:
							tia[0] = 0;
							tia[1] = 2;
							break;
						case 2:
							tia[0] = 1;
							tia[1] = 2;
							break;
					}
					for (z2=0;z2<2;z2++)
					{
						o->m_pCDF->m_p3DF->m_p2DF[z3]->m_iRes[z2] = o->m_pCDF->m_iResolution[tia[z2]];
						switch(g_iObsChannel[tia[z2]])
						{
							case 0:
//								strcat(buf3,"_rdf");
//								strcat(buf3,o->m_pRDF[tia[z2]]->m_sShortName);
								buf3.strcat("_rdf");
								buf3.strcat(o->m_pRDF[tia[z2]]->m_sShortName);
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMinVal[z2] = o->m_pRDF[tia[z2]]->m_fMinDist;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMaxVal[z2] = o->m_pRDF[tia[z2]]->m_fMaxDist;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_pChannels[z2] = o->m_pRDF[tia[z2]]->m_pRDF;
								break;
							case 1:
//								strcat(buf3,"_adf");
//								strcat(buf3,o->m_pADF[tia[z2]]->m_sShortName);
								buf3.strcat("_adf");
								buf3.strcat(o->m_pADF[tia[z2]]->m_sShortName);
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMinVal[z2] = o->m_pADF[tia[z2]]->m_fMinAngle;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMaxVal[z2] = o->m_pADF[tia[z2]]->m_fMaxAngle;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_pChannels[z2] = o->m_pADF[tia[z2]]->m_pADF;
								break;
							case 2:
//								strcat(buf3,"_ddf");
//								strcat(buf3,o->m_pDDF[tia[z2]]->m_sShortName);
								buf3.strcat("_ddf");
								buf3.strcat(o->m_pDDF[tia[z2]]->m_sShortName);
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMinVal[z2] = o->m_pDDF[tia[z2]]->m_fMinAngle;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMaxVal[z2] = o->m_pDDF[tia[z2]]->m_fMaxAngle;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_pChannels[z2] = o->m_pDDF[tia[z2]]->m_pDDF;
								break;
							case 3:
//								strcat(buf3,"_dipole");
//								strcat(buf3,o->m_pDipDF[tia[z2]]->m_sShortName);
								buf3.strcat("_dipole");
								buf3.strcat(o->m_pDipDF[tia[z2]]->m_sShortName);
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMinVal[z2] = o->m_pDipDF[tia[z2]]->m_fDipoleMin;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMaxVal[z2] = o->m_pDipDF[tia[z2]]->m_fDipoleMax;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_pChannels[z2] = o->m_pDipDF[tia[z2]]->m_pDipoleDF;
								break;
							case 4:
//								strcat(buf3,"_vdf");
//								strcat(buf3,o->m_pVDF[tia[z2]]->m_sShortName);
								buf3.strcat("_vdf");
								buf3.strcat(o->m_pVDF[tia[z2]]->m_sShortName);
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMinVal[z2] = o->m_pVDF[tia[z2]]->m_fMinSpeed;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMaxVal[z2] = o->m_pVDF[tia[z2]]->m_fMaxSpeed;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_pChannels[z2] = o->m_pVDF[tia[z2]]->m_pVDF;
								break;
							case 5:
//								strcat(buf3,"_pldf");
//								strcat(buf3,o->m_pPlDF[tia[z2]]->m_sShortName);
								buf3.strcat("_pldf");
								buf3.strcat(o->m_pPlDF[tia[z2]]->m_sShortName);
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMinVal[z2] = o->m_pPlDF[tia[z2]]->m_fMinDist;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMaxVal[z2] = o->m_pPlDF[tia[z2]]->m_fMaxDist;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_pChannels[z2] = o->m_pPlDF[tia[z2]]->m_pPlDF;
								break;
							case 6:
//								strcat(buf3,"_lidf");
//								strcat(buf3,o->m_pLiDF[tia[z2]]->m_sShortName);
								buf3.strcat("_lidf");
								buf3.strcat(o->m_pLiDF[tia[z2]]->m_sShortName);
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMinVal[z2] = o->m_pLiDF[tia[z2]]->m_fMinDist;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fMaxVal[z2] = o->m_pLiDF[tia[z2]]->m_fMaxDist;
								o->m_pCDF->m_p3DF->m_p2DF[z3]->m_pChannels[z2] = o->m_pLiDF[tia[z2]]->m_pLiDF;
								break;
						}
					}
					o->m_pCDF->m_p3DF->m_p2DF[z3]->m_iHistogramRes = o->m_pCDF->m_iHistogramRes;
					o->m_pCDF->m_p3DF->m_p2DF[z3]->Create();
					switch(g_iObsChannel[tia[0]])
					{
						case 0:
//							sprintf(buf2,"%s Distance [pm]",o->m_pRDF[tia[0]]->m_sShortName);
							buf2.sprintf("%s Distance / pm",o->m_pRDF[tia[0]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelX(buf2);
							break;
						case 1:
//							sprintf(buf2,"%s Angle [Degree]",o->m_pADF[tia[0]]->m_sShortName);
							buf2.sprintf("%s Angle / Degree",o->m_pADF[tia[0]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelX(buf2);
							break;
						case 2:
//							sprintf(buf2,"%s Dihedral [Degree]",o->m_pDDF[tia[0]]->m_sShortName);
							buf2.sprintf("%s Dihedral / Degree",o->m_pDDF[tia[0]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelX(buf2);
							break;
						case 3:
//							sprintf(buf2,"%s Dipole moment [Debye]",o->m_pDipDF[tia[0]]->m_sShortName);
							buf2.sprintf("%s Dipole moment / Debye",o->m_pDipDF[tia[0]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelX(buf2);
							break;
						case 4:
//							sprintf(buf2,"%s Velocity [pm/ps]",o->m_pVDF[tia[0]]->m_sShortName);
							buf2.sprintf("%s Velocity / pm ps^-1",o->m_pVDF[tia[0]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelX(buf2);
							break;
						case 5:
//							sprintf(buf2,"%s Distance from Plane [pm]",o->m_pPlDF[tia[0]]->m_sShortName);
							buf2.sprintf("%s Distance from Plane / pm",o->m_pPlDF[tia[0]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelX(buf2);
							break;
						case 6:
//							sprintf(buf2,"%s Distance from Line [pm]",o->m_pLiDF[tia[0]]->m_sShortName);
							buf2.sprintf("%s Distance from Line / pm",o->m_pLiDF[tia[0]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelX(buf2);
							break;
					}
					switch(g_iObsChannel[tia[1]])
					{
						case 0:
//							sprintf(buf2,"%s Distance [pm]",o->m_pRDF[tia[1]]->m_sShortName);
							buf2.sprintf("%s Distance / pm",o->m_pRDF[tia[1]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelY(buf2);
							break;
						case 1:
//							sprintf(buf2,"%s Angle [Degree]",o->m_pADF[tia[1]]->m_sShortName);
							buf2.sprintf("%s Angle / Degree",o->m_pADF[tia[1]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelY(buf2);
							break;
						case 2:
//							sprintf(buf2,"%s Dihedral [Degree]",o->m_pDDF[tia[1]]->m_sShortName);
							buf2.sprintf("%s Dihedral / Degree",o->m_pDDF[tia[1]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelY(buf2);
							break;
						case 3:
//							sprintf(buf2,"%s Dipole moment [Debye]",o->m_pDipDF[tia[1]]->m_sShortName);
							buf2.sprintf("%s Dipole moment / Debye",o->m_pDipDF[tia[1]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelY(buf2);
							break;
						case 4:
//							sprintf(buf2,"%s Velocity [pm/ps]",o->m_pVDF[tia[1]]->m_sShortName);
							buf2.sprintf("%s Velocity / pm ps^-1",o->m_pVDF[tia[1]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelY(buf2);
							break;
						case 5:
//							sprintf(buf2,"%s Distance from Plane [pm]",o->m_pPlDF[tia[1]]->m_sShortName);
							buf2.sprintf("%s Distance from Plane / pm",o->m_pPlDF[tia[1]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelY(buf2);
							break;
						case 6:
//							sprintf(buf2,"%s Distance from Line [pm]",o->m_pLiDF[tia[1]]->m_sShortName);
							buf2.sprintf("%s Distance from Line / pm",o->m_pLiDF[tia[1]]->m_sLabelName);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->SetLabelY(buf2);
							break;
					}
					try { o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sShortName = new char[strlen(buf3)+1]; } catch(...) { o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sShortName = NULL; }
					if (o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sShortName == NULL) NewException((double)(strlen(buf3)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					strcpy(o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sShortName,buf3);

//					sprintf(buf3,"cdf_2");
//					strcat(buf3,o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sShortName);
					buf3.sprintf("cdf_2");
					buf3.strcat(o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sShortName);

					try { o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName = new char[strlen(buf3)+1]; } catch(...) { o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName = NULL; }
					if (o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName == NULL) NewException((double)(strlen(buf3)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					strcpy(o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName,buf3);
				} // END FOR z3
			} // END IF CHANNELS == 3

			try { o->m_pCDF->m_sShortName = new char[strlen(buf)+1]; } catch(...) { o->m_pCDF->m_sShortName = NULL; }
			if (o->m_pCDF->m_sShortName == NULL) NewException((double)(strlen(buf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(o->m_pCDF->m_sShortName,buf);

//			sprintf(buf,"cdf_%d",g_iCDFChannels);
//			strcat(buf,o->m_pCDF->m_sShortName);
			buf.sprintf("cdf_%d",g_iCDFChannels);
			buf.strcat(o->m_pCDF->m_sShortName);

			try { o->m_pCDF->m_sName = new char[strlen(buf)+1]; } catch(...) { o->m_pCDF->m_sName = NULL; }
			if (o->m_pCDF->m_sName == NULL) NewException((double)(strlen(buf)+1)*sizeof(char),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			strcpy(o->m_pCDF->m_sName,buf);

			if (o->m_pCDF->m_bWriteSnapshots) {

				if (!o->m_pCDF->m_bWriteSnapshotsNoWrite) {
					o->m_pCDF->m_sWriteSnapshotsFileName.sprintf( "cdf_%d%s_snapshots.xyz", z+1, o->m_pCDF->m_sShortName );
					mprintf("    Opening output file \"%s\" for writing snapshots...\n",(const char*)o->m_pCDF->m_sWriteSnapshotsFileName);
					o->m_pCDF->m_pWriteSnapshotsFile = OpenFileWrite( (const char*)o->m_pCDF->m_sWriteSnapshotsFileName, true );
				} else
					mprintf("    Not creating snapshot file, only counting configurations...\n");
			}

			if (o->m_pCDF->m_bDumpDat)
			{
//				sprintf(buf,"cdfdump_%dd%s%s.csv",g_iCDFChannels,o->m_pCDF->m_sShortName,multibuf);
				buf.sprintf("cdfdump_%dd%s%s.csv",g_iCDFChannels,o->m_pCDF->m_sShortName,multibuf);
				o->m_pCDF->m_fDump = OpenFileWrite(buf,true);
				mfprintf(o->m_pCDF->m_fDump,"#  step;  RM;  OM1;  OM2;  channels\n");
			}
			if (o->m_bTimeDev)
			{
				if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
				{
					try { o->m_pCDF->m_fTimeDev = new FILE*[1]; } catch(...) { o->m_pCDF->m_fTimeDev = NULL; }
					if (o->m_pCDF->m_fTimeDev == NULL) NewException((double)sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
//					sprintf(buf,"cdf_timedev_%dd%s%s.csv",g_iCDFChannels,o->m_pCDF->m_sShortName,multibuf);
					buf.sprintf("cdf_timedev_%dd%s%s.csv",g_iCDFChannels,o->m_pCDF->m_sShortName,multibuf);
					o->m_pCDF->m_fTimeDev[0] = OpenFileWrite(buf,true);
				} else
				{
					try { o->m_pCDF->m_fTimeDev = new FILE*[o->m_waSaveRefList.GetSize()]; } catch(...) { o->m_pCDF->m_fTimeDev = NULL; }
					if (o->m_pCDF->m_fTimeDev == NULL) NewException((double)o->m_waSaveRefList.GetSize()*sizeof(FILE*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
					{
//						sprintf(buf,"cdf_timedev_%dd%s_ref%d%s.csv",g_iCDFChannels,o->m_pCDF->m_sShortName,z2+1,multibuf);
						buf.sprintf("cdf_timedev_%dd%s_ref%d%s.csv",g_iCDFChannels,o->m_pCDF->m_sShortName,z2+1,multibuf);
						o->m_pCDF->m_fTimeDev[z2] = OpenFileWrite(buf,true);
					}
				}
				if (o->m_pCDF->m_bTDAnimation)
				{
					try { o->m_pCDF->m_pTDAPlot = new CGrace(); } catch(...) { o->m_pCDF->m_pTDAPlot = NULL; }
					if (o->m_pCDF->m_pTDAPlot == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					o->m_pCDF->m_pTDAPlot->SetTitle("CDF Time Development");
					o->m_pCDF->m_pTDAPlot->SetSubTitle(&o->m_pCDF->m_sShortName[1]);
					for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
					{
						for (z3=0;z3<o->m_waSaveShowList.GetSize();z3++)
						{
							for (z4=0;z4<o->m_pCDF->m_iCombinationsEnabled;z4++)
							{
								o->m_pCDF->m_pTDAPlot->AddDataset();
								if (g_iTrajSteps != -1)
									o->m_pCDF->m_pTDAPlot->LastDataset()->m_faValues.SetMaxSize(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pCDF->m_iCombinationsEnabled*g_iTrajSteps);
								o->m_pCDF->m_pTDAPlot->LastDataset()->m_faValues.SetGrow(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pCDF->m_iCombinationsEnabled*100);
							}
						}
					}
					o->m_pCDF->m_pTDAPlot->SetRangeX(o->m_pCDF->m_p2DF->m_fMinVal[0],o->m_pCDF->m_p2DF->m_fMaxVal[0]);
					o->m_pCDF->m_pTDAPlot->SetRangeY(o->m_pCDF->m_p2DF->m_fMinVal[1],o->m_pCDF->m_p2DF->m_fMaxVal[1]);
					o->m_pCDF->m_pTDAPlot->MakeTicks();
					o->m_pCDF->m_pTDAPlot->SetLabelX(o->m_pCDF->m_p2DF->m_sLabelX);
					o->m_pCDF->m_pTDAPlot->SetLabelY(o->m_pCDF->m_p2DF->m_sLabelY);
				}
			}
		} // END IF CDF

		if (g_bVACF)
		{
			mprintf("  Creating VACF...\n");
			o->m_pVACF->Create();
			if (g_bVACFCacheMode)
			{
				mprintf("    VACF Cache: Trying to allocate %s of memory...\n",FormatBytes((double)((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pVACF->m_iShowAtomGes*g_iTrajSteps*3.1*sizeof(double)));
				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pVACF->m_iShowAtomGes;z2++)
				{
					try { ptfa = new CxDoubleArray("main():ptfa"); } catch(...) { ptfa = NULL; }
					if (ptfa == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (g_iTrajSteps != -1)
					{
						ptfa->SetMaxSize((long)(g_iTrajSteps*3.1));
						ptfa->SetGrow((long)(g_iTrajSteps*0.1));
					} else ptfa->SetGrow(1000);
					o->m_pVACF->m_oaCache.Add(ptfa);
				}
			} else
			{
				if (o->m_pVACF->m_iSize > g_iStepHistory)
					g_iStepHistory = o->m_pVACF->m_iSize;
			}
		}

/*		if (g_bDipACF)
		{
			mprintf("  Creating DipACF...\n");
			o->m_pDipACF->Create();
		}*/

		if (g_bMSD)
		{
			mprintf("  Creating MSD...\n");
			o->m_pMSD->m_pMSD->m_fMinVal = 0.0;
			o->m_pMSD->m_pMSD->m_fMaxVal = o->m_pMSD->m_iResolution*g_fTimestepLength*g_iStride/1000.0;
			o->m_pMSD->m_pMSD->m_iResolution = o->m_pMSD->m_iResolution/o->m_pMSD->m_iStride;
			o->m_pMSD->m_pMSD->Create();
			if (g_bMSDCacheMode)
			{
				mprintf("    MSD Cache: Trying to reserve %s of memory...\n",FormatBytes((double)((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pMSD->m_iShowAtoms*g_iTrajSteps/g_iStride*3.1*sizeof(double)));
				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pMSD->m_iShowAtoms;z2++)
				{
					try { ptfa = new CxDoubleArray("main():ptfa"); } catch(...) { ptfa = NULL; }
					if (ptfa == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					if (g_iTrajSteps != -1)
					{
						ptfa->SetMaxSize((long)(g_iTrajSteps/g_iStride*3.1));
						ptfa->SetGrow((long)(g_iTrajSteps/g_iStride*0.1));
					} else ptfa->SetGrow(1000);
					o->m_pMSD->m_oaCache.Add(ptfa);
				}
				if (o->m_pMSD->m_bSplit)
				{
					mprintf("    MSD Split: Trying to reserve %s of memory...\n",FormatBytes((double)((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pMSD->m_iShowAtoms*o->m_pMSD->m_iResolution/o->m_pMSD->m_iStride*sizeof(double)));
					try { o->m_pMSD->m_pSplitMSD = new CAF*[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pMSD->m_iShowAtoms]; } catch(...) { o->m_pMSD->m_pSplitMSD = NULL; }
					if (o->m_pMSD->m_pSplitMSD == NULL) NewException((double)((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pMSD->m_iShowAtoms*sizeof(CAF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pMSD->m_iShowAtoms;z2++)
					{
						try { o->m_pMSD->m_pSplitMSD[z2] = new CAF(); } catch(...) { o->m_pMSD->m_pSplitMSD[z2] = NULL; }
						if (o->m_pMSD->m_pSplitMSD[z2] == NULL) NewException((double)sizeof(CAF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						o->m_pMSD->m_pSplitMSD[z2]->m_fMinVal = 0.0;
						o->m_pMSD->m_pSplitMSD[z2]->m_fMaxVal = o->m_pMSD->m_iResolution*g_fTimestepLength/1000.0;
						o->m_pMSD->m_pSplitMSD[z2]->m_iResolution = o->m_pMSD->m_iResolution/o->m_pMSD->m_iStride;
						o->m_pMSD->m_pSplitMSD[z2]->Create();
					}
				}
			} else // !g_bMSDCacheMode
			{
				if (o->m_pMSD->m_iResolution > g_iStepHistory)
					g_iStepHistory = o->m_pMSD->m_iResolution;
			}
		}

		if (g_bRevSDF)
		{
			mprintf("  Creating Pseudo SDF...\n");
			o->m_pRevSDF->m_p2DF->m_fMinVal[0] = -o->m_pRevSDF->m_fRadius;
			o->m_pRevSDF->m_p2DF->m_fMaxVal[0] = o->m_pRevSDF->m_fRadius;
			o->m_pRevSDF->m_p2DF->m_fMinVal[1] = -o->m_pRevSDF->m_fRadius;
			o->m_pRevSDF->m_p2DF->m_fMaxVal[1] = o->m_pRevSDF->m_fRadius;
			o->m_pRevSDF->m_p2DF->m_iRes[0] = o->m_pRevSDF->m_iResolution;
			o->m_pRevSDF->m_p2DF->m_iRes[1] = o->m_pRevSDF->m_iResolution;
			o->m_pRevSDF->m_p2DF->SetLabelX("X / pm");
			o->m_pRevSDF->m_p2DF->SetLabelY("Y / pm");
			o->m_pRevSDF->m_p2DF->m_iHistogramRes = o->m_pRevSDF->m_iHistogramRes;
			o->m_pRevSDF->m_p2DF->Create();

			try { o->m_pRevSDF->m_vaData = new CxDVec3Array[o->m_iShowMolCount]; } catch(...) { o->m_pRevSDF->m_vaData = NULL; }
			if (o->m_pRevSDF->m_vaData == NULL) NewException((double)o->m_iShowMolCount*sizeof(CxDVec3Array),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}

		if (g_bCond && o->m_bCondDevelopment) {
			buf.sprintf("cond_%d_%s_%s.csv",z+1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName);
			mprintf("  Creating Condition Temporal Development File \"%s\"...\n",(const char*)buf);
			o->m_pCondDevelopment = OpenFileWrite((const char*)buf,true);
			mfprintf(o->m_pCondDevelopment,"Step;  Passed RM Count;  Passed RM Percentage;  Avg. Passed OM Count");
			for (z2=1;z2<=o->m_iCondDevelopmentMax;z2++)
				mfprintf(o->m_pCondDevelopment,";  %d OM Count;  %d OM Percentage",z2,z2);
			mfprintf(o->m_pCondDevelopment,"\n");
			o->m_iaCondDevelopmentCounter.resize(o->m_iCondDevelopmentMax+1);
		}
 	}
	
	if (g_bSDF && g_bSDFMap)
	{
		mprintf("  Creating SDF surface maps...\n");
		for (z=0;z<g_oaSDFMaps.GetSize();z++)
			((CSDFMap*)g_oaSDFMaps[z])->Create();
	}

	if (g_bBondACF)
	{
		mprintf("  Creating BondACF...\n");
		for (z=0;z<g_oaSingleMolecules.GetSize();z++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[z];
			for (z2=0;z2<sm->m_oaBondGroups.GetSize();z2++)
			{
				bg = (CMolBondGroup*)sm->m_oaBondGroups[z2];
				for (z3=0;z3<bg->m_oaBonds.GetSize();z3++)
				{
					if (g_iTrajSteps != -1)
					{
						((CMolBond*)bg->m_oaBonds[z3])->m_faData.SetMaxSize(g_iTrajSteps);
						((CMolBond*)bg->m_oaBonds[z3])->m_faData.SetGrow(g_iTrajSteps/10);
					} else ((CMolBond*)bg->m_oaBonds[z3])->m_faData.SetGrow(10000);
				}
			}
		}
	}

	if (g_bClusterAnalysis)
	{
		mprintf("  Creating Cluster Analysis...\n");
		g_pClusterAnalysis->Create();
	}


	if (g_iRefSystemDim == 3) 
		g_pRefMol.SetSize(((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes);
	if ((g_iRefSystemDim == 3) && !g_bMiddleAvg) // Einfach das erstbeste Molekuel als Referenz nehmen
	{
		mprintf("Creating reference molecule from first time step...  ");
		g_TimeStep.CalcCenters();
		smfix = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[g_iRefMolNum]];
		vec1 = g_TimeStep.m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[0]])->GetAt(g_iFixAtom[0])];
		g_TimeStep.CenterPos(vec1);
		vec2 = g_TimeStep.m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[1]])->GetAt(g_iFixAtom[1])];
		vec3 = g_TimeStep.m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[2]])->GetAt(g_iFixAtom[2])];
//		vec2 -= vec1;
//		vec3 -= vec1;
		mat.MatUltra(vec2,vec3);
		cc = 0;
		// Jeden Atomtyp des Zielmolekuels durchgehen
		for (z3=0;z3<((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex.GetSize();z3++)
		{
			for (z4=0;z4<((CMolecule*)g_oaMolecules[g_iFixMol])->m_waAtomCount[z3];z4++)
			{
				vec2 = g_TimeStep.m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[z3])->GetAt(z4)];
			//	vec2 -= vec1;
				g_pRefMol[cc] = mat * vec2;
				cc++;
			}
		}
		if (g_bPlProj)
		{
			g_TimeStep.FoldAtoms();
			g_TimeStep.Transform(mat);
			for (z2=0;z2<g_oaObserv.GetSize();z2++)
			{
				o = (CObservation*)g_oaObserv[z2];
				if (!o->m_pPlProj->m_bAverageAtomPos)
				{
					ti = 0;
					for (z3=0;z3<o->m_pPlProj->m_oDrawAtoms.m_oaAtoms.GetSize();z3++)
					{
						for (z4=0;z4<((CxIntArray*)o->m_pPlProj->m_oDrawAtoms.m_oaAtoms[z3])->GetSize();z4++)
						{
							o->m_pPlProj->m_vaAtomPos[ti] = g_TimeStep.m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[o->m_pPlProj->m_oDrawAtoms.m_baAtomType[z3]])->GetAt(((CxIntArray*)o->m_pPlProj->m_oDrawAtoms.m_oaAtoms[z3])->GetAt(z4))];
							ti++;
						}
					}
				}
			}
		}
		mprintf("Done.\n");
	} // Ende Referenzbestimmung

	if (g_iSwapAtoms)
	{
		mprintf("Creating Reference Molecule Swap Matrix...\n");

		try { pSwapMatrix = new unsigned int[((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes*((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes]; } catch(...) { pSwapMatrix = NULL; }
		if (pSwapMatrix == NULL) NewException((double)((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes*((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes*sizeof(unsigned int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		for (z=0;z<((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes*((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes;z++)
			pSwapMatrix[z] = 0;
	}

	g_iCurrentTimeStep = -1;
	g_iNextTimeStep = -1;
	g_iLastTimeStep = -1;

	try { g_pTempTimestep = new CTimeStep(); } catch(...) { g_pTempTimestep = NULL; }
	if (g_pTempTimestep == NULL) NewException((double)sizeof(CTimeStep),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if ((g_bSaveRefEnv) && (g_iNbhMode == 3) && (!g_bStreamInput))
	{
		mprintf(WHITE,"\n>>> Pre-analysis for neighborhood search >>>\n");

/*		g_fPos = fopen(g_sInputTraj,"rt"); // Eingabedatei erneut Oeffnen
		if (g_fPos == NULL)
		{
			eprintf("\nError. Input Trajectory suddenly vanished ^^\n");
			return 0;
		}*/
		if (!OpenInputTrajectory())
			goto _ende;
		g_iSteps = 0; // Der Zaehler der Zeitschritte
		if (g_iScanNbhStart != 0)
		{
			mprintf("\nFast-forwarding to step %d...\n",g_iScanNbhStart+1);
			if (g_iTrajFormat == 7)
				mprintf("Warning: Index-based seeking for bqb files not yet implemented.\n");
			mprintf(WHITE,"  [");
			for (z=0;z<g_iScanNbhStart;z++)
			{
				if (fmod(z,g_iScanNbhStart/60.0) < 1.0)
					mprintf(WHITE,"#");
				if (!g_TimeStep.SkipTimestep(g_fPos))
					break;
			}
			mprintf(WHITE,"]\n");
		}
//		while (!feof(g_fPos)) // Zeitschritt fuer Zeitschritt die Trajektorie durchgehen
		while (!InputTrajectoryEOF()) // Zeitschritt fuer Zeitschritt die Trajektorie durchgehen
		{
			for (z=0;z<(g_iScanNbhStride-1);z++)
				if (!g_TimeStep.SkipTimestep(g_fPos))
					goto _endnbs;
			if (!g_TimeStep.ReadTimestep(g_fPos,false))
				goto _endnbs;
			if (g_TimeStep.m_iGesAtomCount == 0)
				goto _endnbs;
//			g_TimeStep.m_vaCoords.SetSize(g_iGesVirtAtomCount);
//			g_TimeStep.UniteMolecules();
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			g_TimeStep.CalcCenters();
	
			if ((g_iSteps % (4*g_iScanNbhStride)) == 0)
			{
				if ((g_iSteps % (200*g_iScanNbhStride)) == 0) 
					mprintf("\nStep %6lu...",g_iSteps);
						else mprintf(".");
			}
	
			g_iSteps+=g_iScanNbhStride;
	
//			vec1 = g_TimeStep.m_vaCoords[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[g_iSaveRefMol]])->m_oaAtomOffset[g_iFixAtomType[0]])->GetAt(g_iFixAtom[0])];
	
//			g_TimeStep.CenterPos(vec1);
	
//			if (g_bFold)
//				g_TimeStep.FoldMolecules();

//			for (z=0;z<g_oaNbSearches.GetSize();z++)
			g_pNbSet->Scan((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[g_iSaveRefMol]],&g_TimeStep);

			if ((g_iScanNbhSteps > 0) && ((int)g_iSteps >= g_iScanNbhSteps))
				break;
		}
		_endnbs:

//		fclose(g_fPos);
		if (!CloseInputTrajectory())
			goto _ende;

		mprintf(WHITE,"\n\n<<< Neighborhood search done <<<\n");
//		for (z0=0;z0<g_oaNbSearches.GetSize();z0++)
		{
			mprintf(YELLOW,"\n*** Choose Neighbors\n");
//			nbset = (CNbSet*)g_oaNbSearches[z0];
			// Sort Neighbors after time they have been neighbors
			for (z=0;z<g_pNbSet->m_oaConditionGroups.GetSize();z++)
			{
				if (g_pNbSet->m_oaConditionGroups[z] == NULL)
					continue;
				cg = (CConditionGroup*)g_pNbSet->m_oaConditionGroups[z];

				try { tpi = new int[((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize()]; } catch(...) { tpi = NULL; }
				if (tpi == NULL) NewException((double)((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize()*sizeof(int),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
				{
					cg->m_bAlwaysTrue[z2] = false;
					tpi[z2] = -1;
				}
				ti3 = 0;
				for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z2++)
				{
					ti = 0;
					for (z3=z2;z3<((CMolecule*)g_oaMolecules[z])->m_laSingleMolIndex.GetSize();z3++)
					{
						for (z4=0;z4<z2;z4++)
							if (tpi[z4] == z3)
								goto _nbhave;
						if (cg->m_iPassCounter[z3] > ti)
						{
							ti = cg->m_iPassCounter[z3];
							ti2 = z3;
						}
_nbhave:;
					}
					if (ti == 0)
						break;
					ti3++;
					tpi[z2] = ti2;
				}
				mprintf(WHITE,"\n  Molecule type %s. %d neighbors found in total:\n",((CMolecule*)g_oaMolecules[z])->m_sName,ti3);
				for (z2=0;z2<ti3;z2++)
					mprintf("    %2d.) Molecule %3d (%8.4f percent of the time, %ld hits)\n",z2+1,tpi[z2]+1,((double)cg->m_iPassCounter[tpi[z2]])*g_iScanNbhStride/g_iSteps*100.0,cg->m_iPassCounter[tpi[z2]]);

				if (ti3 != 0)
				{
					z3 = AskUnsignedInteger("\nUse how many of the frequentiest neighbors for molecule %s? [%d] ",ti3,((CMolecule*)g_oaMolecules[z])->m_sName,ti3);
					for (z2=0;z2<z3;z2++)
						cg->m_bAlwaysTrue[tpi[z2]] = true;
				}

				delete[] tpi;
				cg->m_bInactive = true;
			}
		}
		mprintf("\n");
//		for (z0=0;z0<g_oaNbSearches.GetSize();z0++)
		if (g_bSaveRefWithEnv)
		{
			mprintf("Adding reference molecule to neighborhood...\n");
			g_pNbSet->AddMolecule(g_iFixMol,g_iSaveRefMol);
		}

		g_pNbSet->Reset();
		g_pNbSet->Dump();
	}

	if (g_bVFDF)
	{
		g_iVFCorrCount = 0;
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			for (z2=0;z2<((CMolecule*)g_oaMolecules[z])->m_baAtomIndex.GetSize();z2++)
			{
//				sprintf(buf,"vfcorr_%s_%s%s.dat",((CMolecule*)g_oaMolecules[z])->m_sName,((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[z])->m_baAtomIndex[z2]])->m_sName,multibuf);
				buf.sprintf("vfcorr_%s_%s%s.dat",((CMolecule*)g_oaMolecules[z])->m_sName,(const char*)((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[z])->m_baAtomIndex[z2]])->m_sName,multibuf);
//				FreeFileName(buf);
				g_fVFCorr[g_iVFCorrCount] = OpenFileWrite(buf,true);
				g_iVFCorrCount++;
			}
		}
	}

	if (g_bAggrTopo)
		g_pAggrTopoEngine->Initialize();

	if (g_bContactMatrix)
		g_pContactMatrix->Initialize();

	if (g_bGeoDens)
		g_pGeoDensEngine->Initialize();

	if (g_bTDDF)
		g_pTDDFEngine->Initialize();


	if (g_bSaveRefEnv)
	{
		mprintf("\n");
//		sprintf(g_sRefEnv,"refenv_%s.%d.",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1);
		g_sRefEnv.sprintf("refenv_%s.%d.",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1);
/*		for (z=0;z<g_oaMolecules.GetSize();z++)
			if (g_pNbAll->m_waMolCount[z] != 0)
			{
				sprintf(buf2,"%dx%s.",g_pNbAll->m_waMolCount[z],((CMolecule*)g_oaMolecules[z])->m_sName);
				strcat(g_sRefEnv,buf2);
			}*/
/*		if (!g_bDynamicNeighbor)
			strcat(g_sRefEnv,"static.");
		if (g_bRefEnvCenter)
			strcat(g_sRefEnv,"center.");
		if (g_bRefEnvFix)
			strcat(g_sRefEnv,"fix.");
		if (g_bScanNeighbors)
			strcat(g_sRefEnv,"scannb.");*/

//		strcat(g_sRefEnv,multibuf);
//		strcat(g_sRefEnv,"xyz");

		g_sRefEnv.strcat(multibuf);
		g_sRefEnv.strcat("xyz");

		mprintf(">>> Saving reference environment as %s\n",(const char*)g_sRefEnv);
//		FreeFileName(g_sRefEnv);
		g_fRefEnv = OpenFileWrite(g_sRefEnv,true);
		mprintf("\n");
	}

	if (g_bCutCluster)
	{
		mprintf("\n");
//		sprintf(g_sRefEnv,"cluster_%s.",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
		g_sRefEnv.sprintf("cluster_%s.",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
/*		for (z=0;z<g_oaMolecules.GetSize();z++)
			if (g_pNbAll->m_waMolCount[z] != 0)
			{
				sprintf(buf2,"%dx%s.",g_pNbAll->m_waMolCount[z],((CMolecule*)g_oaMolecules[z])->m_sName);
				strcat(g_sRefEnv,buf2);
			}*/
/*		if (g_bRefEnvCenter)
			strcat(g_sRefEnv,"center.");
		if (g_bRefEnvFix)
			strcat(g_sRefEnv,"fix.");*/

//		strcat(g_sRefEnv,multibuf);
//		strcat(g_sRefEnv,"xyz");

		g_sRefEnv.strcat(multibuf);
		g_sRefEnv.strcat("xyz");

		mprintf(">>> Saving cluster list as %s\n",(const char*)g_sRefEnv);
//		FreeFileName(g_sRefEnv);
		g_fRefEnv = OpenFileWrite(g_sRefEnv,true);

		mprintf("\n");
	}

	if (g_bSaveJustTraj)
	{
//		strcpy(buf,g_sInputTraj);
		buf.strcpy(g_sInputTraj);
		p = strrchr(buf.GetWritePointer(),'/');
		q = strrchr(buf.GetWritePointer(),'\\');
		if (q > p)
			p = q;
		if (p == NULL)
			p = buf.GetWritePointer();
		else
			p++;
//		strcpy(buf2,p);
		buf2.strcpy(p);
		p = strrchr(buf2.GetWritePointer(),'.');
		if (p != NULL)
			*p = 0;

//		strcat(buf2,multibuf);
//		strcat(buf2,"_out.xyz");

		buf2.strcat(multibuf);

		if (g_bProcSplit) {
			procsplitbuf = buf2;
			buf2.strcat("_out.001.xyz");
			procspliti = 0;
			procsplitc = 0;
		} else
			buf2.strcat("_out.xyz");

//		FreeFileName(buf);
//		sprintf(buf,"traj_out.xyz");
		mprintf("Saving processed trajectory as %s ...\n",(const char*)buf2);
		g_fSaveJustTraj = OpenFileWrite(buf2,true);

		if (g_bProcDumpAwkScript) {

			buf = "";
			for (z=0;z<g_oaAtoms.GetSize();z++) {
				if (z == g_iVirtAtomType)
					continue;
				if (((CAtom*)g_oaAtoms[z])->m_iCount == 0)
					continue;
				buf2.sprintf("%s%d",(const char*)((CAtom*)g_oaAtoms[z])->m_sName,((CAtom*)g_oaAtoms[z])->m_iCount);
				buf.strcat(buf2);
			}
			buf2.sprintf("%s",(const char*)buf);
			buf.sprintf("resort_%s.sh",(const char*)buf2);
			mprintf("Saving resorting shell script as %s ...\n",(const char*)buf);
			a = OpenFileWrite(buf,true);

			mfprintf(a,"#!/usr/bin/gawk -f\n");
			mfprintf(a,"\n");
			mfprintf(a,"BEGIN {\n");

			mfprintf(a,"  order=\"");
			if (g_bWriteInputOrder) { // This implies that all non-virtual atoms are to be written
				for (z7=0;z7<g_iGesAtomCount;z7++)
					mfprintf(a,"%d ",z7+1);
			} else {
				if (g_bWriteAtomwise) {
					for (z7=0;z7<g_oaAtoms.GetSize();z7++) {
						for (z3=0;z3<g_oaSaveMolecules.GetSize();z3++) {
							atgr = (CAtomGroup*)g_oaSaveMolecules[z3];
							for (z6=0;z6<atgr->m_pMolecule->m_laSingleMolIndex.GetSize();z6++) {
								for (z4=0;z4<atgr->m_baAtomType.GetSize();z4++) {
									if (atgr->m_baRealAtomType[z4] != z7)
										continue;
									tla = (CxIntArray*)atgr->m_oaAtoms[z4];
									for (z5=0;z5<tla->GetSize();z5++) {
										if (atgr->m_baRealAtomType[z4] == g_iVirtAtomType)
											continue;
										mfprintf(a,"%d ",((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))+1);
									}
								}
							}
						}
					}
				} else {
					for (z3=0;z3<g_oaSaveMolecules.GetSize();z3++) {
						atgr = (CAtomGroup*)g_oaSaveMolecules[z3];
						for (z6=0;z6<atgr->m_pMolecule->m_laSingleMolIndex.GetSize();z6++) {
							for (z4=0;z4<atgr->m_baAtomType.GetSize();z4++) {
								tla = (CxIntArray*)atgr->m_oaAtoms[z4];
								for (z5=0;z5<tla->GetSize();z5++) {
									if (atgr->m_baRealAtomType[z4] == g_iVirtAtomType)
										continue;
									mfprintf(a,"%d ",((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))+1);
								}
							}
						}
					}
				}
			}
			mfprintf(a,"\";\n");


			mfprintf(a,"  labels=\"");
			if (g_bWriteInputOrder) { // This implies that all non-virtual atoms are to be written
				for (z7=0;z7<g_iGesAtomCount;z7++)
				{
					cp = ((CAtom*)g_oaAtoms[g_waAtomRealElement[z7]])->m_sName;
					if (g_bProcAlternativeLabels)
						if (((CxString*)g_oaProcAlternativeLabels[z7])->GetLength() != 0)
							cp = (const char*)*((CxString*)g_oaProcAlternativeLabels[z7]);
					mfprintf(a,"%s ",cp);
				}
			} else {
				if (g_bWriteAtomwise) {
					for (z7=0;z7<g_oaAtoms.GetSize();z7++) {
						for (z3=0;z3<g_oaSaveMolecules.GetSize();z3++) {
							atgr = (CAtomGroup*)g_oaSaveMolecules[z3];
							for (z6=0;z6<atgr->m_pMolecule->m_laSingleMolIndex.GetSize();z6++) {
								for (z4=0;z4<atgr->m_baAtomType.GetSize();z4++) {
									if (atgr->m_baRealAtomType[z4] != z7)
										continue;
									tla = (CxIntArray*)atgr->m_oaAtoms[z4];
									for (z5=0;z5<tla->GetSize();z5++) {
										if (atgr->m_baRealAtomType[z4] == g_iVirtAtomType)
											continue;
										cp = ((CAtom*)g_oaAtoms[atgr->m_baRealAtomType[z4]])->m_sName;
										if (g_bProcAlternativeLabels)
											if (((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))])->GetLength() != 0)
												cp = (const char*)*((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))]);
										mfprintf(a,"%s ",cp);
									}
								}
							}
						}
					}
				} else {
					for (z3=0;z3<g_oaSaveMolecules.GetSize();z3++) {
						atgr = (CAtomGroup*)g_oaSaveMolecules[z3];
						for (z6=0;z6<atgr->m_pMolecule->m_laSingleMolIndex.GetSize();z6++) {
							for (z4=0;z4<atgr->m_baAtomType.GetSize();z4++) {
								tla = (CxIntArray*)atgr->m_oaAtoms[z4];
								for (z5=0;z5<tla->GetSize();z5++) {
									if (atgr->m_baRealAtomType[z4] == g_iVirtAtomType)
										continue;
									cp = ((CAtom*)g_oaAtoms[atgr->m_baRealAtomType[z4]])->m_sName;
									if (g_bProcAlternativeLabels)
										if (((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))])->GetLength() != 0)
											cp = (const char*)*((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))]);
									mfprintf(a,"%s ",cp);
								}
							}
						}
					}
				}
			}
			mfprintf(a,"\";\n");


			char cbuf[256];

			time(&ltime);
			today = localtime(&ltime);
			strcpy(cbuf,asctime(today));
			cbuf[strlen(cbuf)-1] = 0;

			mfprintf(a,"  no=split(order,aorder);\n");
			mfprintf(a,"  nl=split(labels,alabels);\n");
			mfprintf(a,"  qn = 0;\n");
			mfprintf(a,"  i = 0;\n");
			mfprintf(a,"  s = 0;\n");
			mfprintf(a,"  print \"************************************************************\" > \"/dev/stderr\"\n");
			mfprintf(a,"  print \"*** Automatically generated trajectory reordering script ***\" > \"/dev/stderr\"\n");
			mfprintf(a,"  print \"*** Written by TRAVIS at %-31s ***\" > \"/dev/stderr\"\n",cbuf);
			mfprintf(a,"  print \"*** Licensed under GNU GPL v3                            ***\" > \"/dev/stderr\"\n");
			mfprintf(a,"  print \"*** (c) Martin Brehm, 2017                               ***\" > \"/dev/stderr\"\n");
			mfprintf(a,"  print \"*** See http://www.travis-analyzer.de/                   ***\" > \"/dev/stderr\"\n");
			mfprintf(a,"  print \"************************************************************\" > \"/dev/stderr\"\n");
			mfprintf(a,"  print \"\" > \"/dev/stderr\"\n");

			time(&ltime);
			today = localtime(&ltime);
			strcpy(obuf,asctime(today));
			obuf[strlen(obuf)-1] = 0;

			mfprintf(a,"  print \"Script created by TRAVIS at %s\" > \"/dev/stderr\"\n",(const char*)obuf);
			mfprintf(a,"  print \"Suitable for system %s\" > \"/dev/stderr\"\n",(const char*)buf2);
			mfprintf(a,"  print \"\" > \"/dev/stderr\"\n");
			mfprintf(a,"  print \"Usage:  cat inputfile.xyz | ./%s > outputfile.xyz\" > \"/dev/stderr\"\n",(const char*)buf);
			mfprintf(a,"  print \"\" > \"/dev/stderr\"\n");
			mfprintf(a,"  if (no != nl) {\n");
			mfprintf(a,"    print \"*** Internal error:\",no,\"= no != nl =\",nl,\".\" > \"/dev/stderr\"\n");
			mfprintf(a,"    print \"Leaving.\" > \"/dev/stderr\"\n");
			mfprintf(a,"    exit 1\n");
			mfprintf(a,"  }\n");
			mfprintf(a,"}\n");
			mfprintf(a,"\n");
			mfprintf(a,"{\n");
			mfprintf(a,"  if (qn == 0) {\n");
			mfprintf(a,"    n = $1;\n");
			mfprintf(a,"    qn = 1;\n");
			mfprintf(a,"    print \"Trajectory contains\",n,\"atoms.\" > \"/dev/stderr\"\n");
			mfprintf(a,"    if (n != %d) {\n",g_iGesAtomCount);
			mfprintf(a,"      print \"*** Error ***\" > \"/dev/stderr\"\n");
			mfprintf(a,"      print \"This script was created for system %s with %d atoms.\" > \"/dev/stderr\"\n",(const char*)buf2,g_iGesAtomCount);
			mfprintf(a,"      print \"It can only be used for this system.\" > \"/dev/stderr\"\n");
			mfprintf(a,"      print \"Leaving.\" > \"/dev/stderr\"\n");
			mfprintf(a,"      exit 1\n");
			mfprintf(a,"    }\n");
			mfprintf(a,"    print \"Will output\",no,\"atoms.\" > \"/dev/stderr\"\n");
			mfprintf(a,"  }\n");
			mfprintf(a,"  if (i == 0) {\n");
			mfprintf(a,"  } else if (i == 1)\n");
			mfprintf(a,"    comm = $0;\n");
			mfprintf(a,"  else {\n");
			mfprintf(a,"    f1[i-1] = $1;\n");
			mfprintf(a,"    f2[i-1] = $2;\n");
			mfprintf(a,"    f3[i-1] = $3;\n");
			mfprintf(a,"    f4[i-1] = $4;\n");
			mfprintf(a,"    if (NF > 4)\n");
			mfprintf(a,"      fr[i-1] = \"  \" substr($0, index($0,$5));\n");
			mfprintf(a,"    else\n");
			mfprintf(a,"      fr[i-1] = \"  \";\n");
			mfprintf(a,"  }\n");
			mfprintf(a,"  i++;\n");
			mfprintf(a,"  if (i >= n+2) {\n");
			mfprintf(a,"    printf(\"%%d\\n\",no);\n");
			mfprintf(a,"    printf(\"%%s\\n\",comm);\n");
			mfprintf(a,"    for (z=1;z<=no;z++)\n");
			mfprintf(a,"       printf(\"%%-4s %%12.8f  %%12.8f  %%12.8f%%s\\n\",alabels[z],f2[aorder[z]],f3[aorder[z]],f4[aorder[z]],fr[aorder[z]]);\n");
			mfprintf(a,"    i = 0;\n");
			mfprintf(a,"    s++;\n");
			mfprintf(a,"    if ((s %% 100) == 0)\n");
			mfprintf(a,"      print \"Processing step\",s,\"...\" > \"/dev/stderr\"\n");
			mfprintf(a,"  }\n");
			mfprintf(a,"}\n");
			mfprintf(a,"\n");
			mfprintf(a,"END {\n");
			mfprintf(a,"  if (i != 0) {\n");
			mfprintf(a,"    print \"\"\n");
			mfprintf(a,"    print \"Warning: Last frame of input trajectory only contained\",i-2,\"instead of\",n,\"atoms.\" > \"/dev/stderr\"\n");
			mfprintf(a,"    print \"         The incomplete frame has been skipped.\" > \"/dev/stderr\"\n");
			mfprintf(a,"    print \"\"\n");
			mfprintf(a,"  }\n");
			mfprintf(a,"  print \"Finished. Processed\",s,\"steps.\" > \"/dev/stderr\"\n");
			mfprintf(a,"}\n");

			fclose(a);

#ifdef TARGET_LINUX
			buf2.sprintf("chmod +x %s",(const char*)buf);
			mprintf("Executing \"%s\"...\n",(const char*)buf2);
			(void)!system(buf2);
#endif

			mprintf("Finished writing script.\n\n");
			mprintf("Usage:  cat inputfile.xyz | ./%s > outputfile.xyz\n\n",(const char*)buf);
		}
	}


	/******* Interface *******************/
	if (!Interface_Initialization())
		goto _ende;


	mprintf(WHITE,"\n<<< End of Initialization <<<\n\n");

		mprintf(WHITE,"\n##########   Starting Main Pass   ##########\n");
		if (g_bVHDF)
			mprintf("\n    Please note: The VHCF analysis will become slower while it proceeds.\n");

#ifdef TARGET_LINUX
	if (!g_bMultiInterval || (multicounter == 0))
	{
		mprintf(WHITE,"\nHint: ");
		mprintf("Press CTRL+C once to softly interrupt analysis and still write the results.\n");
		mprintf("      Creating an empty file named EXIT (\"touch EXIT\") has the same effect.\n");
	}
#endif


	if (!g_bStreamInput)
	{
/*		g_fPos = fopen(g_sInputTraj,"rb"); 
		if (g_fPos == NULL)
		{
			eprintf("\nCould not open position trajectory.\n");
			goto _ende;
		}*/
		if (g_iTrajFormat == 10)
			g_bDCDSkipHeader = true;

		if (!OpenInputTrajectory())
			goto _ende;
	}

	if ((g_bNPT) && (g_sNPTFile[0] != 0))
	{
		g_fNPTFile = fopen(g_sNPTFile,"rt");
		if (g_fNPTFile == NULL)
		{
			eprintf("\nCould not open cell vector file.\n");
			goto _ende;
		}
	}

	if (g_bUseVelocities && g_bReadVelocity)
	{
		g_fVel = fopen((const char*)g_sInputVel,"rt"); 
		if (g_fVel == NULL)
		{
			eprintf("\nCould not open velocity trajectory.\n");
			goto _ende;
		}
	}
	if (g_bUseForces && (g_sInputForce[0] != 0))
	{
		g_fForce = fopen(g_sInputForce,"rt"); 
		if (g_fForce == NULL)
		{
			eprintf("\nCould not open force trajectory.\n");
			goto _ende;
		}
	}

//	fff = fopen("dipole.txt","wt");

	if (g_iBeginStep != 0)
	{
		mprintf("\nFast-forwarding to step %d...\n",g_iBeginStep+1);

//		mprintf("Seek: %d.\n",g_iFastForwardPos);
//		fseek(g_fPos,g_iFastForwardPos,SEEK_SET);

		if (g_iTrajFormat == 7) {

			mprintf("Warning: Index-based seeking for bqb files not yet implemented.\n");

			mprintf(WHITE,"  [");
			for (z=0;z<g_iBeginStep;z++) {
				if (fmod(z,g_iBeginStep/60.0) < 1.0)
					mprintf(WHITE,"#");
				if (!g_TimeStep.SkipTimestep(g_fPos)) {
					eprintf("\nError: Unexpected end of position trajectory.\n");
					goto _ende;
				}
			}
			mprintf(WHITE,"]\n");

		} else
			SeekInputTrajectory(g_iBeginStep);

		if ((g_fVel != NULL) || (g_fForce != NULL) || ((g_bNPT) && (g_sNPTFile[0] != 0)))
		{
			mprintf(WHITE,"  [");
			for (z=0;z<g_iBeginStep;z++)
			{
/*				if (!g_TimeStep.SkipTimestep(g_fPos))
				{
					eprintf("Error. Unexpected end of position trajectory.\n");
					goto _endmainloop;
				}*/
				if (fmod(z,g_iBeginStep/60.0) < 1.0)
					mprintf(WHITE,"#");
				if (g_fVel != NULL)
				{
					if (!g_TimeStep.SkipTimestep(g_fVel))
					{
						eprintf("Error. Unexpected end of velocity trajectory.\n");
						goto _endmainloop;
					}
				}
				if (g_fForce != NULL)
				{
					if (!g_TimeStep.SkipTimestep(g_fForce))
					{
						eprintf("Error. Unexpected end of force trajectory.\n");
						goto _endmainloop;
					}
				}
				if (((g_bNPT) && (g_sNPTFile[0] != 0)))
//					(void)!fgets(buf,256,g_fNPTFile);
					(void)!buf.fgets(256,g_fNPTFile);
			}
			mprintf(WHITE,"]\n");
		}
	}

	g_oaTimeSteps.SetSize(g_iStepHistory);
	for (z=0;z<g_iStepHistory;z++)
		g_oaTimeSteps[z] = NULL;

	g_iSteps = 0; // Der Zaehler der Zeitschritte
	g_iCurrentTimeStep = 0;
	g_iClusterPos = 0;
	sic = 0;
	g_iSaveCondCount = 0;

	try { apfa = new CxDoubleArray*[g_iCDFChannels]; } catch(...) { apfa = NULL; }
	if (apfa == NULL) NewException((double)g_iCDFChannels*sizeof(CxDoubleArray*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
	
	try { apba = new CxByteArray*[g_iCDFChannels]; } catch(...) { apba = NULL; }
	if (apba == NULL) NewException((double)g_iCDFChannels*sizeof(CxByteArray*),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	try { tda = new double[g_iCDFChannels]; } catch(...) { tda = NULL; }
	if (tda == NULL) NewException((double)g_iCDFChannels*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);

	if (g_bDeriv)
	{
		g_iDerivNext = 2;
		g_iDerivCurr = 1;
		g_iDerivLast = 0;
	}

	if (g_bDipole && g_bDumpDipoleVector)
	{
		g_fDumpDipole = OpenFileWrite("dipole_vectors.csv",true);
		fprintf(g_fDumpDipole,"#Step");
		for (z=0;z<g_oaMolecules.GetSize();z++)
		{
			if (g_oaDumpDipoleVector[z] == NULL)
				continue;
			for (z2=0;z2<((CxIntArray*)g_oaDumpDipoleVector[z])->GetSize();z2++)
			{
				fprintf(g_fDumpDipole,";  %s[%d]_X;  Y;  Z",((CMolecule*)g_oaMolecules[z])->m_sName,((CxIntArray*)g_oaDumpDipoleVector[z])->GetAt(z2)+1);
				if (g_bDumpDipoleAbs)
					fprintf(g_fDumpDipole,";  Abs");
			}
		}
		fprintf(g_fDumpDipole,"\n");

		if (g_bDumpDipoleXYZ)
		{
			g_fDumpDipoleXYZ = OpenFileWrite("dipole_vectors.xyz",true);
			g_fDumpDipolePURE = OpenFileWrite("dipoles.xyz",true);

			fprintf(g_fDumpDipoleXYZ,"%d\n\n",g_iDumpDipoleXYZAtoms);

			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				if (g_oaDumpDipoleVector[z] == NULL)
					continue;
				m = (CMolecule*)g_oaMolecules[z];
				for (z2=0;z2<((CxIntArray*)g_oaDumpDipoleVector[z])->GetSize();z2++)
				{
					sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[((CxIntArray*)g_oaDumpDipoleVector[z])->GetAt(z2)]];

					for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
					{
						if (m->m_baAtomIndex[z3] == g_iVirtAtomType)
							continue;
						for (z4=0;z4<((CxIntArray*)sm->m_oaAtomOffset[z3])->GetSize();z4++)
						{
							ti = ((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(z4);
							fprintf(g_fDumpDipoleXYZ,"%s  %12f  %12f  %12f\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z3]])->m_sName,g_TimeStep.m_vaCoords[ti][0]/100.0,g_TimeStep.m_vaCoords[ti][1]/100.0,g_TimeStep.m_vaCoords[ti][2]/100.0);
						}
					}
				}
			}

			ti2 = 0;
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				if (g_oaDumpDipoleVector[z] == NULL)
					continue;
				m = (CMolecule*)g_oaMolecules[z];
				for (z2=0;z2<((CxIntArray*)g_oaDumpDipoleVector[z])->GetSize();z2++)
				{
					fprintf(g_fDumpDipoleXYZ,"B  %12f  0    0\n",100.0+ti2*3.0);
					fprintf(g_fDumpDipoleXYZ,"B  %12f  0.5  0\n",100.0+ti2*3.0);
					ti2++;
				}
			}
		}
	}


/*****************************************************************************
*************** Beginn Hauptschleife *****************************************
*****************************************************************************/

	g_iLastTimeStep = -1;
	g_iDotCounter = 0;
	g_bStepSkipped = false;
	g_iFirstStepSkipped = -1;
	g_iSaveTrajPatternPos = 0;
	tbs = false;

	t0 = (unsigned long)time(NULL);

	while (true) // Zeitschritt fuer Zeitschritt die Trajektorie durchgehen
	{
		if (g_bCubeStream) {
			char filename[1024];
			if (g_iSteps > 0) {
				printCubeMemFileName(filename, 1024, false);
				if (remove(filename) != 0) {
					eprintf("Could not remove \"%s\": %s\n", filename, strerror(errno));
				}
			}
			if (g_iSteps >= (unsigned int)g_iCubeMemFileSteps) {
				mprintf("\n\n%lu cube files read.\n", g_iSteps);
				break;
			}
			printCubeMemFileName(filename, 1024, true);
			g_fCubeMemFile->ReadFileSuccessive(filename, g_iCubeMemFileLines, false);
		}
		
//		if (feof(g_fPos))
		if (InputTrajectoryEOF())
		{
			mprintf("\nEnd of trajectory file reached.\n");
			break;
		}
		if (g_bAbortAnalysis)
		{
			mprintf("\nAnalysis aborted by user.\n");
			break;
		}
		g_iCurrentTimeStep++;
		if (g_iCurrentTimeStep >= g_iStepHistory)
			g_iCurrentTimeStep = 0;

		if (g_bDeriv)
		{
			g_iDerivNext++;
			if (g_iDerivNext > 2)
				g_iDerivNext = 0;
			g_iDerivCurr++;
			if (g_iDerivCurr > 2)
				g_iDerivCurr = 0;
			g_iDerivLast++;
			if (g_iDerivLast > 2)
				g_iDerivLast = 0;
		}

		for (z=0;z<(g_iStride-1);z++)
		{
			if (!g_TimeStep.SkipTimestep(g_fPos))
			{
//				if (!feof(g_fPos))
				if (!InputTrajectoryEOF())
					eprintf("\nError while skipping time step %lu.\n",g_iSteps+z+1);
				else
					mprintf("\nEnd of trajectory file reached while skipping steps.\n");
				goto _endmainloop;
			}
			if (((g_bNPT) && (g_sNPTFile[0] != 0)))
//				(void)!fgets(buf,256,g_fNPTFile);
				(void)!buf.fgets(256,g_fNPTFile);

			if (g_fVel != NULL)
			{
				if (!g_TimeStep.SkipTimestep(g_fVel))
				{
					eprintf("\nError while skipping velocity time step %lu.\n",g_iSteps+z+1);
					goto _endmainloop;
				}
			}
			if (g_fForce != NULL)
			{
				if (!g_TimeStep.SkipTimestep(g_fForce))
				{
					eprintf("\nError while skipping force time step %lu.\n",g_iSteps+z+1);
					goto _endmainloop;
				}
			}
		}

		if (g_bSaveTrajPattern) {

			while (g_sSaveTrajPattern[ g_iSaveTrajPatternPos % g_sSaveTrajPattern.GetLength() ] == '0') {

				if (!g_TimeStep.SkipTimestep(g_fPos)) {
					if (!InputTrajectoryEOF())
						eprintf("\nError while skipping time step %lu.\n",g_iSteps+z+1);
					else
						mprintf("\nEnd of trajectory file reached while skipping steps.\n");
					goto _endmainloop;
				}
				if (((g_bNPT) && (g_sNPTFile[0] != 0)))
					(void)!buf.fgets(256,g_fNPTFile);
				if (g_fVel != NULL) {
					if (!g_TimeStep.SkipTimestep(g_fVel)) {
						eprintf("\nError while skipping velocity time step %lu.\n",g_iSteps+z+1);
						goto _endmainloop;
					}
				}
				if (g_fForce != NULL) {
					if (!g_TimeStep.SkipTimestep(g_fForce)) {
						eprintf("\nError while skipping force time step %lu.\n",g_iSteps+z+1);
						goto _endmainloop;
					}
				}
				g_iSteps += g_iStride;
				g_iSaveTrajPatternPos++;
			}

			g_iSaveTrajPatternPos++;
		}

_readagain:
		if (g_bUseVelocities || g_bUseForces)
		{
			if (GetTimeStep(-1) == NULL)
			{
				try { *GetTimeStepAddress(-1) = new CTimeStep(); } catch(...) { *GetTimeStepAddress(-1) = NULL; }
				if (*GetTimeStepAddress(-1) == NULL) NewException((double)sizeof(CTimeStep),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}
			bool ok = false;
			if (g_bCubeStream)
				ok = GetTimeStep(-1)->ReadTimestep(g_fCubeMemFile);
			else
				ok = GetTimeStep(-1)->ReadTimestep(g_fPos,false);
			if (!ok)
			{
//				if (feof(g_fPos))
				if (InputTrajectoryEOF())
				{
					mprintf("\nEnd of Trajectory File reached.\n");
					break;
				}
				eprintf("\nError while reading time step %lu.\n",g_iSteps+1);
				goto _endmainloop;
			}
			if (((g_bNPT) && (g_sNPTFile[0] != 0)))
				if (!GetTimeStep(-1)->ReadCellVector(g_fNPTFile))
					break;
			if (g_bSkipDoubleSteps)
			{
				if (GetTimeStep(-1)->ExtractNumber(g_iNumberPos) <= g_iLastTimeStep)
				{
					if (!g_bStepSkipped)
					{
						g_bStepSkipped = true;
						g_iFirstStepSkipped = GetTimeStep(-1)->ExtractNumber(g_iNumberPos);
						mprintf("\nSkipping:");
					}
					mprintf("*");
					goto _readagain;
				} else if (g_iFirstStepSkipped != -1)
				{
					if (g_iFirstStepSkipped == GetTimeStep(-1)->ExtractNumber(g_iNumberPos))
						mprintf("\nRepeated step %d skipped.",g_iFirstStepSkipped);
							else mprintf("\nRepeated steps %d-%ld skipped.",g_iFirstStepSkipped,GetTimeStep(-1)->ExtractNumber(g_iNumberPos));
					g_iDotCounter = 0;
					g_iFirstStepSkipped = -1;
				}

				if (g_iLastTimeStep != -1)
				{
					if (g_iStrideDetect == -1)
					{
						g_iStrideDetect = GetTimeStep(-1)->ExtractNumber(g_iNumberPos) - g_iLastTimeStep;
					} else
					{
						if (g_iStrideDetect != (GetTimeStep(-1)->ExtractNumber(g_iNumberPos) - g_iLastTimeStep))
						{
							if (!tbs)
								eprintf("\n");
							eprintf("Warning: Two successive steps %lu (%d) and %lu (%ld) have different distance than seen before: %d.\n",g_iSteps-1,g_iLastTimeStep,g_iSteps,GetTimeStep(-1)->ExtractNumber(g_iNumberPos),g_iStrideDetect);
							tbs = true;
						} else tbs = false;
					}
				}

				if (GetTimeStep(-1)->ExtractNumber(g_iNumberPos) > g_iLastTimeStep)
					g_iLastTimeStep = GetTimeStep(-1)->ExtractNumber(g_iNumberPos);
			} 
			if (GetTimeStep(-1)->m_iGesAtomCount == 0)
			{
				eprintf("\nError: Atom count = 0 at time step %lu.\n",g_iSteps+1);
				goto _endmainloop;
			}
			GetTimeStep(-1)->UniteMolecules(false);
			if (g_bRemoveCOM)
				GetTimeStep(-1)->CenterCOM();
			GetTimeStep(-1)->CalcCenters();
			if (g_bDipole)
			{
				if (g_bWannier)
					GetTimeStep(-1)->ScanWannier(false);
// 				GetTimeStep(-1)->CalcDipoles(false);
// 				if (g_bDumpDipoleVector)
// 					GetTimeStep(-1)->DumpDipoles();
			}
		} else
		{
			if (GetTimeStep(0) == NULL)
			{
				try { *GetTimeStepAddress(0) = new CTimeStep(); } catch(...) { *GetTimeStepAddress(0) = NULL; }
				if (*GetTimeStepAddress(0) == NULL) NewException((double)sizeof(CTimeStep),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			}
			bool ok = false;
			if (g_bCubeStream)
				ok = GetTimeStep(0)->ReadTimestep(g_fCubeMemFile);
			else
				ok = GetTimeStep(0)->ReadTimestep(g_fPos,false);
			if (!ok)
			{
//				if (feof(g_fPos))
				if (InputTrajectoryEOF())
				{
					mprintf("\nEnd of trajectory file reached.\n");
					
					// This fix allows to use GetTimeStep(0) after the end of the main loop, like the Raman module does
					g_iCurrentTimeStep--;
					if (g_iCurrentTimeStep < 0)
						g_iCurrentTimeStep += g_iStepHistory;
					// End of Fix

					break;
				}
				eprintf("\nError while reading time step %lu.\n",g_iSteps+1);
				goto _endmainloop;
			}
			if (((g_bNPT) && (g_sNPTFile[0] != 0)))
				if (!GetTimeStep(0)->ReadCellVector(g_fNPTFile))
					break;
			if (g_bSkipDoubleSteps)
			{
				if (GetTimeStep(0)->ExtractNumber(g_iNumberPos) <= g_iLastTimeStep)
				{
//					mprintf("\n[Skip %d/%d]",GetTimeStep(0)->ExtractNumber(),g_iLastTimeStep);
					if (!g_bStepSkipped)
					{
						g_bStepSkipped = true;
						g_iFirstStepSkipped = GetTimeStep(0)->ExtractNumber(g_iNumberPos);
						mprintf("\nSkipping:");
					}
					mprintf("*");
					goto _readagain;
				} else if (g_iFirstStepSkipped != -1)
				{
					if (g_iFirstStepSkipped == GetTimeStep(0)->ExtractNumber(g_iNumberPos))
						mprintf("\nRepeated step %d skipped.",g_iFirstStepSkipped);
							else mprintf("\nRepeated steps %d-%ld skipped.",g_iFirstStepSkipped,GetTimeStep(0)->ExtractNumber(g_iNumberPos));
					g_iDotCounter = 0;
					g_iFirstStepSkipped = -1;
				}
//				mprintf("\nNumber %d.",GetTimeStep(0)->ExtractNumber());

				if (g_iLastTimeStep != -1)
				{
					if (g_iStrideDetect == -1)
					{
						g_iStrideDetect = GetTimeStep(0)->ExtractNumber(g_iNumberPos) - g_iLastTimeStep;
					} else
					{
						if (g_iStrideDetect != (GetTimeStep(0)->ExtractNumber(g_iNumberPos) - g_iLastTimeStep))
						{
							if (!tbs)
								eprintf("\n");
							eprintf("Warning: Successive steps %lu (%d) and %lu (%ld) have different distance than seen before: %d.\n",g_iSteps-1,g_iLastTimeStep,g_iSteps,GetTimeStep(0)->ExtractNumber(g_iNumberPos),g_iStrideDetect);
							tbs = true;
						} else tbs = false;
					}
				}

				if (GetTimeStep(0)->ExtractNumber(g_iNumberPos) > g_iLastTimeStep)
					g_iLastTimeStep = GetTimeStep(0)->ExtractNumber(g_iNumberPos);
			}
			if (GetTimeStep(0)->m_iGesAtomCount == 0)
			{
				eprintf("\nError: Atom count = 0 at time step %lu.\n",g_iSteps+1);
				goto _endmainloop;
			}
			if (!g_bSaveCoordsUnchanged)
			{
				GetTimeStep(0)->UniteMolecules(false);
				if (g_bRemoveCOM)
					GetTimeStep(0)->CenterCOM();
			}
			GetTimeStep(0)->CalcCenters();
			if (g_bDipole)
			{
				if (g_bWannier)
					GetTimeStep(0)->ScanWannier(false);
// 				GetTimeStep(0)->CalcDipoles(false);
// 				if (g_bDumpDipoleVector)
// 					GetTimeStep(0)->DumpDipoles();
			}
		}
		g_bWarnUnsteady = false;

		if (((g_iSteps-(sic*g_iStride*50)) % (showinterval*g_iStride)) == 0)
		{
			if ((g_iSteps == 0) || g_bStepSkipped)
			{
				g_bStepSkipped = false;
				if (!g_bSilentProgress)
					mprintf("\nStep %6lu ",g_iSteps);
			}
			if (g_bAsciiArt && (!g_bSilentProgress))
			{
				tc = g_oAsciiArt.GetChar();
				if (tc == 0)
				{
					if ((g_iSteps != 0) && (g_iTrajSteps != -1))
					{
						if ((time(NULL) - t0) > 5)
						{
							eta = (unsigned long)(((double)time(NULL) - t0) / g_iSteps * ((double)MAX(long(0),(((g_iMaxStep > 0)?g_iMaxStep:g_iTrajSteps) - ((long)g_iSteps)))));
							FormatTime(eta,&buf);
							mprintf(" ETA %s",(const char*)buf);
						}
					}
					mprintf("\nStep %6lu ",g_iSteps);
					if (g_iSteps != 0)
						g_oAsciiArt.NewLine();
				} else mprintf("%c",tc);
			} else
			{
				if (g_iDotCounter >= 50)
				{
					g_iDotCounter = 0;
					if (!g_bSilentProgress)
					{
						if ((g_iSteps != 0) && (g_iTrajSteps != -1))
						{
							if ((time(NULL) - t0) > 5)
							{
								eta = (unsigned long)(((double)time(NULL) - t0) / g_iSteps * ((double)MAX(long(0),g_iTrajSteps - ((long)g_iSteps))));
								FormatTime(eta,&buf);
								mprintf(" ETA %s",(const char*)buf);
							}
						}
						mprintf("\nStep %6lu ",g_iSteps);
					}
				}
				g_iDotCounter++;
				if (!g_bSilentProgress)
					mprintf(".");
			}
			if (FileExist("EXIT"))
			{
				mprintf("\n\n*** File \"EXIT\" detected. Aborting analysis. ***\n\n");
				remove("EXIT");
				break;
			}
		}

		if ((int)g_iSteps == /*showinterval**/g_iStride*50)
		{
			t1 = (unsigned long)time(NULL);
			if (t1 == t0)
				showinterval = 20;
					else showinterval = (int)(20.0/(t1-t0));
			if (showinterval == 0)
				showinterval = 1;
			sic = 1;
		}

		g_iSteps += g_iStride;

		if (GetTimeStep(0)==NULL)
			continue;

		if (((g_bUseVelocities && !g_bReadVelocity) || g_bUseForces) && (GetTimeStep(1)==NULL))
			continue;

		if (g_bUseVelocities)
		{
			if (g_fVel == NULL)
			{
				CalcVelocities();
			} else if (!GetTimeStep(-1)->ReadTimestepVel(g_fVel))
			{
				eprintf("\nError reading velocity time step %lu.\n",g_iSteps+1);
				goto _endmainloop;
			}
		}
		if (g_bUseForces)
		{
			if (g_fForce == NULL)
			{
				CalcForces();
			} else if (!GetTimeStep(-1)->ReadTimestepForce(g_fForce))
			{
				eprintf("\nError reading force time step %lu.\n",g_iSteps+1);
				goto _endmainloop;
			}
		}
		
		if (g_bCubeTimeDev) {
			CalcVolumetricDataTimeDev();
			CalcCurrentDensity();
		}
		
		if (g_bTegri) {
			g_pTetraPak->ProcessStep(GetTimeStep(0), 0);
		}
		
		if (g_bDipole) {
			GetTimeStep(0)->CalcDipoles(false);
			if (g_bDumpDipoleVector)
				GetTimeStep(0)->DumpDipoles();
		}
		
		if (g_bMagneticDipoleRestart || g_bVCD) {
			GetTimeStep(0)->CalcMagneticDipoles();
		}

		if (g_bTegri)
			if (g_pTetraPak->m_bSaveAtomIntegrals)
				DumpMolecularProps();
		
		if (g_bCutCluster)
		{
			if (((int)g_iSteps/g_iStride) >= g_iClusterSteps)
				break;

			if (g_iClusterPos >= g_iClusterCount)
				break;

			if (g_iaClusterSteps[g_iClusterPos] != ((int)g_iSteps/g_iStride))
				continue;
		}
	//	mprintf("\nStep %d for clusters.",((int)g_iSteps/g_iStride));

		if (g_bUnwrap && ((int)g_iSteps > g_iStride))
 // Nicht im ersten Schritt
		{
			for (z=0;z<g_oaMolecules.GetSize();z++)
			{
				m = (CMolecule*)g_oaMolecules[z];
				for (z3=0;z3<m->m_baAtomIndex.GetSize();z3++)
				{
					if (m->m_baAtomIndex[z3] != g_iVirtAtomType)
						continue;
					for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
					{
						sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
						if (GetTimeStep(0)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][0] - GetTimeStep(1)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][0] > g_fBoxX/2.0)
							g_vaUnwrapArray[m->m_laSingleMolIndex[z2]][0] -= g_fBoxX;
						if (GetTimeStep(0)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][0] - GetTimeStep(1)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][0] < -g_fBoxX/2.0)
							g_vaUnwrapArray[m->m_laSingleMolIndex[z2]][0] += g_fBoxX;
						if (GetTimeStep(0)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][1] - GetTimeStep(1)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][1] > g_fBoxY/2.0)
							g_vaUnwrapArray[m->m_laSingleMolIndex[z2]][1] -= g_fBoxY;
						if (GetTimeStep(0)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][1] - GetTimeStep(1)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][1] < -g_fBoxY/2.0)
							g_vaUnwrapArray[m->m_laSingleMolIndex[z2]][1] += g_fBoxY;
						if (GetTimeStep(0)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][2] - GetTimeStep(1)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][2] > g_fBoxZ/2.0)
							g_vaUnwrapArray[m->m_laSingleMolIndex[z2]][2] -= g_fBoxZ;
						if (GetTimeStep(0)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][2] - GetTimeStep(1)->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z3])->GetAt(0)][2] < -g_fBoxZ/2.0)
							g_vaUnwrapArray[m->m_laSingleMolIndex[z2]][2] += g_fBoxZ;
	//					mprintf("\n%2d.%2d: %f, %f, %f",z,z2,g_vaUnwrapArray[m->m_laSingleMolIndex[z2]][0],g_vaUnwrapArray[m->m_laSingleMolIndex[z2]][1],g_vaUnwrapArray[m->m_laSingleMolIndex[z2]][2]);
					}
				}
			}
		}

		if (g_bKeepUnfoldedCoords)
			GetTimeStep(0)->m_vaCoords_Unfolded.CopyFrom(&GetTimeStep(0)->m_vaCoords);

		if (g_bFixedPlProj)
			g_pFixedPlProj->Process(GetTimeStep(0));

		/******* Interface *********/
		Interface_ProcessStep(GetTimeStep(0));

		if (g_bMSD)
		{
			g_pT2Timestep = GetTimeStep(0);
			for (z0=0;z0<g_oaObserv.GetSize();z0++)
			{
				ti2 = 0;
				o = (CObservation*)g_oaObserv[z0];
				if (g_bMSDCacheMode)
				{
					for (z2=0;z2<o->m_pMSD->m_pAtomGroup->m_oaAtoms.GetSize();z2++)
					{
						for (z3=0;z3<((CxIntArray*)o->m_pMSD->m_pAtomGroup->m_oaAtoms[z2])->GetSize();z3++)
						{
							for (z4=0;z4<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();z4++)
							{
								ti = ((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex[z4]])->m_oaAtomOffset[o->m_pMSD->m_pAtomGroup->m_baAtomType[z2]])->GetAt(((CxIntArray*)o->m_pMSD->m_pAtomGroup->m_oaAtoms[z2])->GetAt(z3));
								ptfa = (CxDoubleArray*)o->m_pMSD->m_oaCache[ti2];

									vec2 = g_pT2Timestep->m_vaCoords[ti];

								if (g_bPeriodic && (ptfa->GetSize() > 0))
								{
									ti3 = ptfa->GetSize()-3;
									vec1[0] = ptfa->GetAt(ti3++);
									vec1[1] = ptfa->GetAt(ti3++);
									vec1[2] = ptfa->GetAt(ti3);
									tf = (vec2-vec1).GetLength();
										if (tf > g_fMinPeriodic/2.0)
											eprintf("\nDiscontinuity in step %lu: %s[%d] %s%d moved %.2f pm. Can't compute MSD from wrapped trajectory - please use unwrapped trajectory.",g_iSteps,((CMolecule*)g_oaMolecules[g_waAtomMolIndex[ti]])->m_sName,g_laAtomSMLocalIndex[ti]+1,(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[ti]])->m_sName,g_waAtomMolNumber[ti]+1,tf);
								}

								ptfa->Add( vec2[0] );
								ptfa->Add( vec2[1] );
								ptfa->Add( vec2[2] );

								ti2++;
							}
						}
					}
				} else
				{
					for (z=0;z<o->m_pMSD->m_iResolution/o->m_pMSD->m_iStride;z++)
					{
						g_pTempTimestep = GetTimeStep(z*o->m_pMSD->m_iStride);
						if (g_pTempTimestep == NULL)
							continue;
						for (z2=0;z2<o->m_pMSD->m_pAtomGroup->m_oaAtoms.GetSize();z2++)
						{
							for (z3=0;z3<((CxIntArray*)o->m_pMSD->m_pAtomGroup->m_oaAtoms[z2])->GetSize();z3++)
							{
								for (z4=0;z4<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();z4++)
								{
									ti = ((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex[z4]])->m_oaAtomOffset[o->m_pMSD->m_pAtomGroup->m_baAtomType[z2]])->GetAt(((CxIntArray*)o->m_pMSD->m_pAtomGroup->m_oaAtoms[z2])->GetAt(z3));
									o->m_pMSD->m_pMSD->AddToBin_Index(z,(double)(pow2(g_pTempTimestep->m_vaCoords[ti][0]-g_pT2Timestep->m_vaCoords[ti][0])+pow2(g_pTempTimestep->m_vaCoords[ti][1]-g_pT2Timestep->m_vaCoords[ti][1])+pow2(g_pTempTimestep->m_vaCoords[ti][2]-g_pT2Timestep->m_vaCoords[ti][2])));
								}
							}
						}
					}
				}
			}
		} // END IF MSD 

		if (g_bSaveJustTraj)
		{
			g_pTempTimestep->CopyFrom(GetTimeStep(0));

			if (g_bSaveJustCenter && (!g_bSaveCoordsUnchanged))
			{
				vec0 = g_pTempTimestep->m_vaCoords[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iSaveJustMol])->m_laSingleMolIndex[g_iSaveJustSM]])->m_oaAtomOffset[g_iSaveJustAtomType])->GetAt(g_iSaveJustAtom)];
				g_pTempTimestep->CenterPos(vec0);
			}

			if (g_bSaveTrajNoRot)
			{
				g_pTempTimestep->CenterCOM();

//				mfprintf(g_fSaveJustTraj,"  %d\n\n",(int)g_iSaveGesAtoms*3);

//				for (z=0;z<g_iGesAtomCount;z++)
//					mfprintf(g_fSaveJustTraj,"%2s  %8.5f  %8.5f  %8.5f\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,g_pTempTimestep->m_vaCoords[z][0]/100.0,g_pTempTimestep->m_vaCoords[z][1]/100.0,g_pTempTimestep->m_vaCoords[z][2]/100.0);

		//		g_pTempTimestep->Transform(tdmat2);
				g_pTempTimestep->Transform(tq2);

//				for (z=0;z<g_iGesAtomCount;z++)
//					mfprintf(g_fSaveJustTraj,"%2s  %8.5f  %8.5f  %8.5f\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,g_pTempTimestep->m_vaCoords[z][0]/100.0+10.0,g_pTempTimestep->m_vaCoords[z][1]/100.0,g_pTempTimestep->m_vaCoords[z][2]/100.0);

				if ((int)g_iSteps > g_iStride)
				{
					dvec1 = CxDVector3(0);
					for (z=0;z<g_iGesAtomCount;z++)
					{
						dvec2 = (g_pTempTimestep->m_vaCoords[z] - rot_ts.m_vaCoords[z]);
						dvec4 = (g_pTempTimestep->m_vaCoords[z] + rot_ts.m_vaCoords[z]) / 2.0;

		//				tf3 = dvec4.GetLength();
		//				tf2 = dvec2.GetLength();
		//				if (4*tf3*tf3 > 0.25*tf2*tf2)
		//				{
		//					dvec4 *= (1.0 + 1.0 - 0.5 / tf3 * sqrt(4*tf3*tf3 - 0.25*tf2*tf2) );
		//					mprintf("\n%.10g",1.0 + 1.0 - 0.5 / tf3 * sqrt(4*tf3*tf3 - 0.25*tf2*tf2));
		//				}

						dvec3 = CrossP(dvec4,dvec2);
						dvec1 += dvec3 * ((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_fMass;
//						mprintf("  %2d: ( %g, %g, %g ) x ( %g, %g, %g ) = ( %g, %g, %g ), m = %.2f\n",z,dvec2[0],dvec2[1],dvec2[2],rot_ts.m_vaCoords[z][0],rot_ts.m_vaCoords[z][1],rot_ts.m_vaCoords[z][2],dvec3[0],dvec3[1],dvec3[2],((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_fMass);
					}
//					mprintf("Drehimpuls: %g ( %g, %g, %g )\n",dvec1.GetLength(),dvec1[0],dvec1[1],dvec1[2]);
					dvec3 = dvec1;
					dvec3.Normalize();
	//				mprintf("      ( %g, %g, %g )\n",vec3[0],vec3[1],vec3[2]);

					tf = 0;
					for (z=0;z<g_iGesAtomCount;z++)
					{
						dvec4 = (g_pTempTimestep->m_vaCoords[z] + rot_ts.m_vaCoords[z]) / 2.0;

		//				tf3 = dvec4.GetLength();
		//				tf2 = dvec2.GetLength();
		//				dvec4 *= (1.0 + 1.0 - 0.5 / tf3 * sqrt(4*tf3*tf3 - 0.25*tf2*tf2) );

						dvec0 = CrossP(dvec4,dvec3);
	//					mprintf("#  ( %g, %g, %g ) x ( %g, %g, %g ) = %g\n",g_pTempTimestep->m_vaCoords[z][0],g_pTempTimestep->m_vaCoords[z][1],g_pTempTimestep->m_vaCoords[z][2],vec3[0],vec3[1],vec3[2],vec0.GetLength());
						tf += ((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_fMass * dvec0.GetLengthSqr();
		//				mprintf("\n  %2d: ( %g, %g, %g ) x ( %g, %g, %g ) = ( %g, %g, %g ), m = %.2f",z,vec2[0],vec2[1],vec2[2],rot_ts.m_vaCoords[z][0],rot_ts.m_vaCoords[z][1],rot_ts.m_vaCoords[z][2],vec3[0],vec3[1],vec3[2],((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_pElement->m_fMass);
					}
//					mprintf("  Traegheitsmoment: %g\n",tf);
					dvec1 /= tf;
//					mprintf("Angular momentum vector: ( %g, %g, %g ) (l = %g)\n",dvec1[0],dvec1[1],dvec1[2],dvec1.GetLength());
//					mprintf("  Winkelgeschw.: %g\n",dvec1.GetLength());
//					dvec3 = dvec1;
//					dvec3.Normalize();

//					mprintf("  Baue Matrix: ( %g | %g | %g ), %g\n",dvec3[0],dvec3[1],dvec3[2],-dvec1.GetLength());

//					tdmat.RotMat(dvec3,-dvec1.GetLength());

					tq.BuildRotation(dvec3,dvec1.GetLength());

//					mprintf("  Die Matrix:\n");
//					tdmat.Dump();

//					tdmat2 = tdmat2 * tdmat;
					tq2 = tq2 * tq;

					if (fabs(tq2.GetLength()-1.0) > 0.000001)
					{
						eprintf("\nRenormalizing rotation quaternion.\n");
						tq2.Normalize();
					}

//					g_pTempTimestep->Transform(tdmat);
					g_pTempTimestep->Transform(tq);
				} 

				rot_ts.CopyFrom(g_pTempTimestep);

//				for (z=0;z<g_iGesAtomCount;z++)
//					mfprintf(g_fSaveJustTraj,"%2s  %8.5f  %8.5f  %8.5f\n",((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,g_pTempTimestep->m_vaCoords[z][0]/100.0+20.0,g_pTempTimestep->m_vaCoords[z][1]/100.0,g_pTempTimestep->m_vaCoords[z][2]/100.0);

			} else
			{
				if (g_bFold && (!g_bSaveCoordsUnchanged))
				{
					if (g_bFoldAtomwise)
						g_pTempTimestep->FoldAtoms();
							else g_pTempTimestep->FoldMolecules();
				}
			}

			if (!g_bCenterZero && (!g_bSaveCoordsUnchanged) && g_bPeriodic)
			{
				if (g_bBoxNonOrtho)
				{
					vec1 = g_mBoxFromOrtho * CxDVector3(-0.5,-0.5,-0.5);
					g_pTempTimestep->CenterPos(vec1);
				} else
				{
					g_pTempTimestep->CenterPos(CxDVector3(-g_fBoxX/2,-g_fBoxY/2,-g_fBoxZ/2));
				}
			}

			ti = g_iSaveGesAtoms;

			if (g_bUnwrapWannier)
			{
				g_pTempTimestep->ScanWannier(false);

				for (z3=0;z3<g_oaSaveMolecules.GetSize();z3++)
				{
					atgr = (CAtomGroup*)g_oaSaveMolecules[z3];
					for (z6=0;z6<atgr->m_pMolecule->m_laSingleMolIndex.GetSize();z6++)
						ti += ((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_laWannier.GetSize();
				}
			}

			if (g_bProcAddMesh) {
				if (g_bProcAddMeshX)
					ti += (g_iProcAddMeshGrid-1)*(g_iProcAddMeshGrid-1);
				if (g_bProcAddMeshY)
					ti += (g_iProcAddMeshGrid-1)*(g_iProcAddMeshGrid-1);
				if (g_bProcAddMeshZ)
					ti += (g_iProcAddMeshGrid-1)*(g_iProcAddMeshGrid-1);
				if (g_bProcAddMeshX || g_bProcAddMeshY)
					ti += (g_iProcAddMeshGrid-1);
				if (g_bProcAddMeshX || g_bProcAddMeshZ)
					ti += (g_iProcAddMeshGrid-1);
				if (g_bProcAddMeshY || g_bProcAddMeshZ)
					ti += (g_iProcAddMeshGrid-1);
				if (g_bProcAddMeshX || g_bProcAddMeshY || g_bProcAddMeshZ)
					ti++;
			}

			mfprintf(g_fSaveJustTraj,"  %d\n",ti);
			if (g_bProcCellComment) {
				if (g_bProcCellCommentAngles)
					mfprintf( g_fSaveJustTraj, "    %.6f  %.6f  %.6f  %.6f  %.6f  %.6f\n",
						g_fBoxX/100.0, g_fBoxY/100.0, g_fBoxZ/100.0, g_fBoxAngleA, g_fBoxAngleB, g_fBoxAngleC );
				else
					mfprintf( g_fSaveJustTraj, "    %.6f  %.6f  %.6f\n",
						g_fBoxX/100.0, g_fBoxY/100.0, g_fBoxZ/100.0 );
			} else if (g_pTempTimestep->m_pComment != NULL)
				mfprintf(g_fSaveJustTraj,"%s\n",g_pTempTimestep->m_pComment);
			else
				mfprintf(g_fSaveJustTraj,"No Comment\n");

			if (g_bBoxNonOrtho && g_bWriteOrtho)
			{
				for (z7=0;z7<g_iGesAtomCount;z7++)
				{
					z3 = g_waAtomMolIndex[z7]; // Molecule
					z6 = g_laAtomSMLocalIndex[z7]; // SingleMolecule
					z4 = g_waAtomElement[z7]; // AtomType
					z5 = g_waAtomMolNumber[z7]; // Atom

					m = (CMolecule*)g_oaMolecules[z3];

					vec0 = g_pTempTimestep->m_vaCoords[z7];
					vec1 = g_mBoxToOrtho * vec0;

					while (vec1[0] < 0)
						vec1[0] += 1.0;
					while (vec1[0] >= 1.0)
						vec1[0] -= 1.0;
					while (vec1[1] < 0)
						vec1[1] += 1.0;
					while (vec1[1] >= 1.0)
						vec1[1] -= 1.0;
					while (vec1[2] < 0)
						vec1[2] += 1.0;
					while (vec1[2] >= 1.0)
						vec1[2] -= 1.0;

					vec1 *= g_fWriteOrthoFac; 

					cp = ((CAtom*)g_oaAtoms[g_waAtomRealElement[z7]])->m_sName;
					if (g_bProcAlternativeLabels)
						if (((CxString*)g_oaProcAlternativeLabels[z7])->GetLength() != 0)
							cp = (const char*)*((CxString*)g_oaProcAlternativeLabels[z7]);

					mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f\n",cp,vec1[0]/100.0,vec1[1]/100.0,vec1[2]/100.0);
				}
			} else
			{
				if (g_bWriteInputOrder) // This implies that all non-virtual atoms are to be written
				{
					for (z7=0;z7<g_iGesAtomCount;z7++)
					{
						z3 = g_waAtomMolIndex[z7]; // Molecule
						z6 = g_laAtomSMLocalIndex[z7]; // SingleMolecule
						z4 = g_waAtomElement[z7]; // AtomType
						z5 = g_waAtomMolNumber[z7]; // Atom

						m = (CMolecule*)g_oaMolecules[z3];

						vec0 = g_pTempTimestep->m_vaCoords[z7];
						if (g_bUnwrap)
						{
							vec0 += g_vaUnwrapArray[m->m_laSingleMolIndex[z6]];
				//			mprintf("\nUnwrap: (%G|%G|%G) <-- (%G|%G|%G)",vec0[0],vec0[1],vec0[2],g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]][0],g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]][1],g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]][2]);
						}

						cp = ((CAtom*)g_oaAtoms[g_waAtomRealElement[z7]])->m_sName;
						if (g_bProcAlternativeLabels)
							if (((CxString*)g_oaProcAlternativeLabels[z7])->GetLength() != 0)
								cp = (const char*)*((CxString*)g_oaProcAlternativeLabels[z7]);

						mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f",cp,vec0[0]/100.0,vec0[1]/100.0,vec0[2]/100.0);
						if (g_bProcWriteComments)
							mfprintf(g_fSaveJustTraj,"  # %s[%d] %s%d",m->m_sName,z6+1,(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[z7]])->m_sName,z5+1);
						mfprintf(g_fSaveJustTraj,"\n");
					}
				} else
				{
					if (g_bWriteAtomwise)
					{
						for (z7=0;z7<g_oaAtoms.GetSize();z7++)
						{
							for (z3=0;z3<g_oaSaveMolecules.GetSize();z3++)
							{
								atgr = (CAtomGroup*)g_oaSaveMolecules[z3];
								for (z6=0;z6<atgr->m_pMolecule->m_laSingleMolIndex.GetSize();z6++)
								{
									for (z4=0;z4<atgr->m_baAtomType.GetSize();z4++)
									{
										if (atgr->m_baRealAtomType[z4] != z7)
											continue;
										tla = (CxIntArray*)atgr->m_oaAtoms[z4];
										for (z5=0;z5<tla->GetSize();z5++)
										{
											vec0 = g_pTempTimestep->m_vaCoords[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))];
											if (g_bUnwrap)
											{
												vec0 += g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]];
									//			mprintf("\nUnwrap: (%G|%G|%G) <-- (%G|%G|%G)",vec0[0],vec0[1],vec0[2],g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]][0],g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]][1],g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]][2]);
											}
											if (atgr->m_baRealAtomType[z4] == g_iVirtAtomType)
											{
												if (g_bSaveVirtAtoms) {
													cp = ((CVirtualAtom*)g_oaVirtualAtoms[atgr->m_pMolecule->m_laVirtualAtoms[tla->GetAt(z5)]])->m_sLabel;
													if (g_bProcAlternativeLabels)
														if (((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))])->GetLength() != 0)
															cp = (const char*)*((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))]);
													mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f",cp,vec0[0]/100.0,vec0[1]/100.0,vec0[2]/100.0);
												}
											} else {
												cp = ((CAtom*)g_oaAtoms[atgr->m_baRealAtomType[z4]])->m_sName;
												if (g_bProcAlternativeLabels)
													if (((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))])->GetLength() != 0)
														cp = (const char*)*((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))]);

												mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f",cp,vec0[0]/100.0,vec0[1]/100.0,vec0[2]/100.0);
											}
											if (g_bProcWriteComments)
												mfprintf(g_fSaveJustTraj,"  # %s[%d] %s%d",atgr->m_pMolecule->m_sName,z6+1,(const char*)((CAtom*)g_oaAtoms[atgr->m_baRealAtomType[z4]])->m_sName,tla->GetAt(z5)+1);
											mfprintf(g_fSaveJustTraj,"\n");
										}
									}
								}
							}
						}
					} else
					{
						for (z3=0;z3<g_oaSaveMolecules.GetSize();z3++)
						{
							atgr = (CAtomGroup*)g_oaSaveMolecules[z3];
							for (z6=0;z6<atgr->m_pMolecule->m_laSingleMolIndex.GetSize();z6++)
							{
								for (z4=0;z4<atgr->m_baAtomType.GetSize();z4++)
								{
									tla = (CxIntArray*)atgr->m_oaAtoms[z4];
									for (z5=0;z5<tla->GetSize();z5++)
									{
										vec0 = g_pTempTimestep->m_vaCoords[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))];
										if (g_bUnwrap)
										{
											vec0 += g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]];
								//			mprintf("\nUnwrap: (%G|%G|%G) <-- (%G|%G|%G)",vec0[0],vec0[1],vec0[2],g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]][0],g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]][1],g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]][2]);
										}
										if (atgr->m_baRealAtomType[z4] == g_iVirtAtomType)
										{
											if (g_bSaveVirtAtoms) {
												cp = ((CVirtualAtom*)g_oaVirtualAtoms[atgr->m_pMolecule->m_laVirtualAtoms[tla->GetAt(z5)]])->m_sLabel;
												if (g_bProcAlternativeLabels)
													if (((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))])->GetLength() != 0)
														cp = (const char*)*((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))]);
												mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f",cp,vec0[0]/100.0,vec0[1]/100.0,vec0[2]/100.0);
											}
										} else {
											cp = ((CAtom*)g_oaAtoms[atgr->m_baRealAtomType[z4]])->m_sName;
											if (g_bProcAlternativeLabels)
												if (((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))])->GetLength() != 0)
													cp = (const char*)*((CxString*)g_oaProcAlternativeLabels[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_oaAtomOffset[atgr->m_baAtomType[z4]])->GetAt(tla->GetAt(z5))]);

											mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f",cp,vec0[0]/100.0,vec0[1]/100.0,vec0[2]/100.0);
										}
										if (g_bProcWriteComments)
											mfprintf(g_fSaveJustTraj,"  # %s[%d] %s%d",atgr->m_pMolecule->m_sName,z6+1,(const char*)((CAtom*)g_oaAtoms[atgr->m_baRealAtomType[z4]])->m_sName,tla->GetAt(z5)+1);
										mfprintf(g_fSaveJustTraj,"\n");
									}
								}
							}
						}
					}
				}
			}

			if (g_bUnwrapWannier)
			{
				for (z3=0;z3<g_oaSaveMolecules.GetSize();z3++)
				{
					atgr = (CAtomGroup*)g_oaSaveMolecules[z3];
					for (z6=0;z6<atgr->m_pMolecule->m_laSingleMolIndex.GetSize();z6++)
					{
						for (z4=0;z4<((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_laWannier.GetSize();z4++)
						{
							vec0 = g_pTempTimestep->m_vaCoords[((CSingleMolecule*)g_oaSingleMolecules[atgr->m_pMolecule->m_laSingleMolIndex[z6]])->m_laWannier[z4]];
							if (g_bUnwrap)
								vec0 += g_vaUnwrapArray[atgr->m_pMolecule->m_laSingleMolIndex[z6]];
							mfprintf(g_fSaveJustTraj,"X    %12.8f  %12.8f  %12.8f\n",vec0[0]/100.0,vec0[1]/100.0,vec0[2]/100.0);
						}
					}
				}	
			}

			if (g_bProcAddMesh) {

				if (g_bCenterZero) {
					tfx = -g_fBoxX/200.0;
					tfy = -g_fBoxY/200.0;
					tfz = -g_fBoxZ/200.0;
				} else {
					tfx = 0;
					tfy = 0;
					tfz = 0;
				}

				if (g_bProcAddMeshX || g_bProcAddMeshY || g_bProcAddMeshZ)
					mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f\n",(const char*)g_sProcAddMeshLabel,tfx+MeshRandom(),tfy+MeshRandom(),tfz+MeshRandom());

				if (g_bProcAddMeshX || g_bProcAddMeshY)
					for (z3=1;z3<g_iProcAddMeshGrid;z3++)
						mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f\n",(const char*)g_sProcAddMeshLabel,tfx+MeshRandom(),tfy+MeshRandom(),tfz+MeshRandom()+(double)z3/g_iProcAddMeshGrid*g_fBoxZ/100.0);

				if (g_bProcAddMeshX || g_bProcAddMeshZ)
					for (z3=1;z3<g_iProcAddMeshGrid;z3++)
						mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f\n",(const char*)g_sProcAddMeshLabel,tfx+MeshRandom(),tfy+MeshRandom()+(double)z3/g_iProcAddMeshGrid*g_fBoxY/100.0,tfz+MeshRandom());

				if (g_bProcAddMeshY || g_bProcAddMeshZ)
					for (z3=1;z3<g_iProcAddMeshGrid;z3++)
						mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f\n",(const char*)g_sProcAddMeshLabel,tfx+MeshRandom()+(double)z3/g_iProcAddMeshGrid*g_fBoxX/100.0,tfy+MeshRandom(),tfz+MeshRandom());

				if (g_bProcAddMeshX)
					for (z3=1;z3<g_iProcAddMeshGrid;z3++)
						for (z4=1;z4<g_iProcAddMeshGrid;z4++)
							mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f\n",(const char*)g_sProcAddMeshLabel,tfx+MeshRandom(),tfy+MeshRandom()+(double)z3/g_iProcAddMeshGrid*g_fBoxY/100.0,tfz+MeshRandom()+(double)z4/g_iProcAddMeshGrid*g_fBoxZ/100.0);

				if (g_bProcAddMeshY)
					for (z3=1;z3<g_iProcAddMeshGrid;z3++)
						for (z4=1;z4<g_iProcAddMeshGrid;z4++)
							mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f\n",(const char*)g_sProcAddMeshLabel,tfx+MeshRandom()+(double)z3/g_iProcAddMeshGrid*g_fBoxX/100.0,tfy+MeshRandom(),tfz+MeshRandom()+(double)z4/g_iProcAddMeshGrid*g_fBoxZ/100.0);

				if (g_bProcAddMeshZ)
					for (z3=1;z3<g_iProcAddMeshGrid;z3++)
						for (z4=1;z4<g_iProcAddMeshGrid;z4++)
							mfprintf(g_fSaveJustTraj,"%-4s %12.8f  %12.8f  %12.8f\n",(const char*)g_sProcAddMeshLabel,tfx+MeshRandom()+(double)z3/g_iProcAddMeshGrid*g_fBoxX/100.0,tfy+MeshRandom()+(double)z4/g_iProcAddMeshGrid*g_fBoxY/100.0,tfz+MeshRandom());
			}

			if (g_bProcSplit) {
				procsplitc++;
				if (procsplitc >= g_iProcSplitLength) {

					procsplitc = 0;
					procspliti++;

					fclose(g_fSaveJustTraj);

					buf2.sprintf("%s_out.%03d.xyz",(const char*)procsplitbuf,procspliti+1);

					g_fSaveJustTraj = OpenFileWrite(buf2,true);
				}
			}
		} // END IF g_bSaveJustTraj

		if (g_bBondACF)
		{
			g_pTempTimestep->CopyFrom(GetTimeStep(0));
			for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[z2];
				for (z3=0;z3<sm->m_oaBonds.GetSize();z3++)
				{
					bond = (CMolBond*)sm->m_oaBonds[z3];
					bond->m_faData.Add(VecDist(g_pTempTimestep->m_vaCoords[bond->m_iAtomOffset[0]],g_pTempTimestep->m_vaCoords[bond->m_iAtomOffset[1]]));
				}
			}
		}
		
		if (g_bAggregation || g_bNbExchange)
		{
			g_pTempTimestep->CopyFrom(GetTimeStep(0));

			for (z0=0;z0<g_oaObserv.GetSize();z0++)
			{
				o = (CObservation*)g_oaObserv[z0];
				for (zs=0;zs<o->m_pDACF->m_oaSubDACFs.GetSize();zs++)
				{
					dacfsub = (CDACFSub*)o->m_pDACF->m_oaSubDACFs[zs];
					for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize();z2++)
						dacfsub->m_pCondition->m_iPassCounter[z2] = 0;
				}
			}
			for (z0=0;z0<g_oaObserv.GetSize();z0++)
			{
				o = (CObservation*)g_oaObserv[z0];

				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_pDACF->m_iFirstMol])->m_laSingleMolIndex.GetSize();z2++)
				{
					smfix = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[o->m_pDACF->m_iFirstMol])->m_laSingleMolIndex[z2]];

					o->m_pDACF->m_pCondition->PreScanNeighborhoodAllOM(g_pTempTimestep,smfix);

					for (zs=0;zs<o->m_pDACF->m_oaSubDACFs.GetSize();zs++)
					{
						dacfsub = (CDACFSub*)o->m_pDACF->m_oaSubDACFs[zs];

						dacfsub->m_pCondition->CopyResults(o->m_pDACF->m_pCondition);
						dacfsub->m_pCondition->m_bAnyPassed = false;
						dacfsub->m_pCondition->m_iTempPassed = 0;
						dacfsub->m_pCondition->ReScan(smfix);
						if (dacfsub->m_pCondition->m_bAnyPassed)
							dacfsub->m_pCondition->m_iRMPassCounter[z2]++;
	
						if (g_bDACF)
							o->m_pDACF->UpdateDACFSub(z2,g_pTempTimestep,dacfsub);

						if (g_bNbExchange)
							o->m_pDACF->UpdateNbEx(z2,dacfsub);
					} // END FOR ZS
				} // END FOR Z2
			} // END FOR Z0
			for (z0=0;z0<g_oaObserv.GetSize();z0++)
			{
				o = (CObservation*)g_oaObserv[z0];
				for (zs=0;zs<o->m_pDACF->m_oaSubDACFs.GetSize();zs++)
				{
					dacfsub = (CDACFSub*)o->m_pDACF->m_oaSubDACFs[zs];
					for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize();z2++)
						if (dacfsub->m_pCondition->m_iPassCounter[z2] != 0)
							dacfsub->m_pCondition->m_iOMPassCounter[z2]++;
				}
			}
		} // END IF AGGREGATION OR NBEX

		if (g_bDens)
		{
			g_pTempTimestep->CopyFrom(GetTimeStep(0));
			for (z0=0;z0<g_oaObserv.GetSize();z0++)
			{
				o = (CObservation*)g_oaObserv[z0];
				m = (CMolecule*)g_oaMolecules[o->m_iShowMol];
				for (z3=0;z3<m->m_laSingleMolIndex.GetSize();z3++)
				{
					if (g_bRegionAnalysis)
						if ((!o->m_iaRMRegions.Contains(0)) && (!o->m_iaRMRegions.Contains(g_iaSMRegion[m->m_laSingleMolIndex[z3]])))
							continue;

					sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z3]];
					ti = ((CxIntArray*)sm->m_oaAtomOffset[o->m_pDensityDF->m_iCenterAtomType])->GetAt(o->m_pDensityDF->m_iCenterAtom);
					vec0 = g_pTempTimestep->m_vaCoords[ti];

					g_pTempTimestep->CenterPos(vec0);
					g_pTempTimestep->FoldAtoms();

					for (z4=0;z4<g_oaMolecules.GetSize();z4++)
					{
						if (!o->m_pDensityDF->m_pDensityMolSelect[z4])
							continue;
						m2 = (CMolecule*)g_oaMolecules[z4];
						ag = o->m_pDensityDF->m_pDensityMolAG[z4];
						for (z5=0;z5<m2->m_laSingleMolIndex.GetSize();z5++)
						{
							if (g_bRegionAnalysis)
								if ((!o->m_iaOM1Regions.Contains(0)) && (!o->m_iaOM1Regions.Contains(g_iaSMRegion[m2->m_laSingleMolIndex[z5]])))
									continue;

							sm2 = (CSingleMolecule*)g_oaSingleMolecules[m2->m_laSingleMolIndex[z5]];
							for (z6=0;z6<ag->m_baAtomType.GetSize();z6++)
							{
								if (ag->m_baRealAtomType[z6] == g_iVirtAtomType)
									continue;
								if (o->m_pDensityDF->m_bDensityMass)
								{
									tf2 = ((CAtom*)g_oaAtoms[ag->m_baRealAtomType[z6]])->m_pElement->m_fMass;
									for (z7=0;z7<((CxIntArray*)ag->m_oaAtoms[z6])->GetSize();z7++)
									{
										tf = g_pTempTimestep->m_vaCoords[((CxIntArray*)sm2->m_oaAtomOffset[ag->m_baAtomType[z6]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z6])->GetAt(z7))].GetLength();
										o->m_pDensityDF->m_pDensDF->AddToBin(tf,tf2);
									}
								} else
								{
									for (z7=0;z7<((CxIntArray*)ag->m_oaAtoms[z6])->GetSize();z7++)
									{
										tf = g_pTempTimestep->m_vaCoords[((CxIntArray*)sm2->m_oaAtomOffset[ag->m_baAtomType[z6]])->GetAt(((CxIntArray*)ag->m_oaAtoms[z6])->GetAt(z7))].GetLength();
										o->m_pDensityDF->m_pDensDF->AddToBin(tf);
									}
								}
							}
						}
					}
				}
			}
		} // END IF g_bDens

		if (g_bRDyn || g_bIRSpec)
		{
			g_pTempTimestep->CopyFrom(GetTimeStep(0));
			for (z6=0;z6<g_oaObserv.GetSize();z6++)
			{
				o = (CObservation*)g_oaObserv[z6];
				if (g_bIRSpec)
					trdyn = o->m_pIRSpec;
						else trdyn = o->m_pRDyn;
				ti = 0;
				for (z3=0;z3<o->m_iShowMolCount;z3++) // Alle anderen Molekuele durchgehen
				{
					sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex[z3]];

					if (trdyn->m_iVecType == 0)
					{
						trdyn->BuildAtomList(sm,&templa);
						for (z4=0;z4<templa.GetSize();z4+=3)
						{
							if (trdyn->m_bOrtho)
							{
								vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
								vec2 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
								vec3 = g_pTempTimestep->m_vaCoords[templa[z4+2]];
								vec1 = CrossP(vec2-vec0,vec3-vec0);
							} else
							{
								vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
								vec2 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
								vec1 = vec2-vec0;
							}
							vec1.Normalize();
							if (g_bRDynCacheMode)
							{
								((CxDoubleArray*)trdyn->m_oaCache[ti])->Add(vec1[0]);
								((CxDoubleArray*)trdyn->m_oaCache[ti])->Add(vec1[1]);
								((CxDoubleArray*)trdyn->m_oaCache[ti])->Add(vec1[2]);
								ti++;
							} else
							{
								for (z=0;z<trdyn->m_iDepth/trdyn->m_iStride;z++)
								{
									g_pT2Timestep = GetTimeStep(z*trdyn->m_iStride);
									if (g_pT2Timestep == NULL)
										continue;
									if (trdyn->m_bOrtho)
									{
										vec4 = g_pT2Timestep->m_vaCoords[templa[z4]];
										vec3 = g_pT2Timestep->m_vaCoords[templa[z4+1]];
										vec5 = g_pT2Timestep->m_vaCoords[templa[z4+2]];
										vec2 = CrossP(vec3-vec4,vec5-vec4);
									} else
									{
										vec4 = g_pT2Timestep->m_vaCoords[templa[z4]];
										vec3 = g_pT2Timestep->m_vaCoords[templa[z4+1]];
										vec2 = vec3-vec4;
									}
									vec2.Normalize();
									trdyn->m_pRDyn->AddToBin_Int(z,DotP(vec1,vec2));
									trdyn->m_pCount[z]++;
								}
							}
						}
					} else if (trdyn->m_iVecType == 1)
					{
						vec1 = sm->m_vDipole;
						((CxDoubleArray*)trdyn->m_oaCache[ti])->Add(vec1[0]);
						((CxDoubleArray*)trdyn->m_oaCache[ti])->Add(vec1[1]);
						((CxDoubleArray*)trdyn->m_oaCache[ti])->Add(vec1[2]);

	//					if (z3 == 0)
	//						mfprintf(fff,"%f, %f, %f\n",vec1[0],vec1[1],vec1[2]);

						ti++;
					} // END VECTYPE
				}
			}
		} // END IF RDYN

		if (g_bVACF)
		{
			g_pT2Timestep = GetTimeStep(0);
//			mprintf("\nStep=%d, adding.",g_iSteps);
			if (g_bGlobalVACF)
			{
				if (g_bVACFCacheMode)
				{
					ti = 0;
					for (z2=0;z2<g_iGesAtomCount;z2++)
					{
						if (g_pGlobalVACF->m_bExcludeR0 && (g_waAtomRealElement[z2] > 1000))
							continue;
		//				mprintf("Global: Step %d, Atom %d: %.4ff | %.4ff | %.4ff\n",g_iSteps,z2+1,g_pT2Timestep->m_vaVelocities[z2][0],g_pT2Timestep->m_vaVelocities[z2][1],g_pT2Timestep->m_vaVelocities[z2][2]);
						((CxDoubleArray*)g_pGlobalVACF->m_oaCache[ti])->Add(g_pT2Timestep->m_vaVelocities[z2][0]);
						((CxDoubleArray*)g_pGlobalVACF->m_oaCache[ti])->Add(g_pT2Timestep->m_vaVelocities[z2][1]);
						((CxDoubleArray*)g_pGlobalVACF->m_oaCache[ti])->Add(g_pT2Timestep->m_vaVelocities[z2][2]);
						ti++;
					}
				} else
				{
					for (z=0;z<g_pGlobalVACF->m_iSize;z++)
					{
						g_pTempTimestep = GetTimeStep(z);
						if (g_pTempTimestep == NULL)
							continue;
						if (g_pTempTimestep->m_vaVelocities.GetSize() == 0)
							continue;
						pd = &g_pGlobalVACF->m_pData[z];
						g_pGlobalVACF->m_pCounter[z] += g_iGesAtomCount;
						for (z2=0;z2<g_iGesAtomCount;z2++)
							*pd += DotP(g_pTempTimestep->m_vaVelocities[z2],g_pT2Timestep->m_vaVelocities[z2]);
					}
				}
			}
			for (z6=0;z6<g_oaObserv.GetSize();z6++)
			{
				o = (CObservation*)g_oaObserv[z6];
				for (z3=0;z3<o->m_iShowMolCount;z3++) // Alle anderen Molekuele durchgehen
				{
					sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex[z3]];
					o->m_pVACF->BuildAtomList(sm,&templa);
					for (z4=0;z4<templa.GetSize();z4++)
					{
						if (g_bVACFCacheMode)
						{
							ti = templa[z4];
							((CxDoubleArray*)o->m_pVACF->m_oaCache[z3*o->m_pVACF->m_iShowAtomGes+z4])->Add(g_pT2Timestep->m_vaVelocities[ti][0]);
							((CxDoubleArray*)o->m_pVACF->m_oaCache[z3*o->m_pVACF->m_iShowAtomGes+z4])->Add(g_pT2Timestep->m_vaVelocities[ti][1]);
							((CxDoubleArray*)o->m_pVACF->m_oaCache[z3*o->m_pVACF->m_iShowAtomGes+z4])->Add(g_pT2Timestep->m_vaVelocities[ti][2]);
				//			mprintf("%d*%d+%d=%d\n",z3,o->m_pVACF->m_iShowAtomGes,z4,z3*o->m_pVACF->m_iShowAtomGes+z4);
				//			mprintf("Lokal: Step %d, Atom %d: %.4ff | %.4ff | %.4ff\n",g_iSteps,ti+1,g_pT2Timestep->m_vaVelocities[ti][0],g_pT2Timestep->m_vaVelocities[ti][1],g_pT2Timestep->m_vaVelocities[ti][2]);
				//			mprintf("Step %d: Obs %d, Mol %d, Atom %d: Index %d.\n",g_iStepHistory,z6+1,z3+1,z4+1,ti+1);
						} else
						{
							for (z5=0;z5<o->m_pVACF->m_iSize;z5++)
							{
								o->m_pVACF->m_pCounter[z5]++;
								g_pTempTimestep = GetTimeStep(z5);
								if (g_pTempTimestep == NULL)
									continue;
								if (g_pTempTimestep->m_vaVelocities.GetSize() == 0)
									continue;
								o->m_pVACF->m_pData[z5] += DotP(g_pTempTimestep->m_vaVelocities[templa[z4]],g_pT2Timestep->m_vaVelocities[templa[z4]]);
							}
						}
					}
				}
			}
		} // END IF VACF

		if (g_bIRSpec && g_bGlobalIR)
		{
			for (z3=0;z3<g_oaSingleMolecules.GetSize();z3++) // Alle anderen Molekuele durchgehen
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[z3];
				vec1 = sm->m_vDipole;
				((CxDoubleArray*)g_pGlobalIR->m_oaCache[z3])->Add(vec1[0]);
				((CxDoubleArray*)g_pGlobalIR->m_oaCache[z3])->Add(vec1[1]);
				((CxDoubleArray*)g_pGlobalIR->m_oaCache[z3])->Add(vec1[2]);
			}
		} // END IF g_bGlobalIR


		if (g_bNbAnalysis)
		{
			g_pTempTimestep->CopyFrom(GetTimeStep(0));
			for (z2=0;z2<((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize();z2++)
			{
				smfix = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[z2]];
				for (z6=0;z6<g_oaObserv.GetSize();z6++)
				{
					o = (CObservation*)g_oaObserv[z6];
					o->m_pNbAnalysis->m_pNbSearch->ScanAllOM(smfix,g_pTempTimestep);
					o->m_pNbAnalysis->AnalyzeStep();
				}
			}
		}

		if (g_bNbExchange)
		{
		}

		if (g_bOrder)
			g_pOrderEngine->ProcessStep(GetTimeStep(0));

		if (g_bAggrTopo)
			g_pAggrTopoEngine->ProcessStep(GetTimeStep(0));

		if (g_bContactMatrix)
			g_pContactMatrix->ProcessStep(GetTimeStep(0));

		if (g_bTDDF)
			g_pTDDFEngine->ProcessStep(GetTimeStep(0));

		if (g_bGeoDens)
			g_pGeoDensEngine->ProcessStep(GetTimeStep(0));

		if (g_bClusterAnalysis)
		{
			g_pTempTimestep->CopyFrom(GetTimeStep(0));
			g_pClusterAnalysis->Process(g_pTempTimestep);
		}


		if (g_bVoro)
		{
			g_pTempTimestep->CopyFrom(GetTimeStep(0));

			if (g_pVoroWrapper->m_bVoroStat || g_pVoroWrapper->m_bVoroNbMatrix)
			{
				g_pTempTimestep->FoldAtomsPositive();
				g_pVoroWrapper->Build(g_pTempTimestep);
			}


			if (g_pVoroWrapper->m_bSurfCover)
				g_pVoroWrapper->ProcessSurfCover(g_pTempTimestep);

		}

		if (g_bDomA)
		{
			g_pTempTimestep->CopyFrom(GetTimeStep(0));
			g_pTempTimestep->FoldAtomsPositive();
			g_pDomainEngine->ProcessStep(g_pTempTimestep);
		}

		if (g_bVDF && (g_iFixMol == -1))
		{
			g_pTempTimestep->CopyFrom(GetTimeStep(0));
			for (z6=0;z6<g_oaObserv.GetSize();z6++)
			{
				o = (CObservation*)g_oaObserv[z6];
				for (z3=0;z3<o->m_iShowMolCount;z3++) // Alle anderen Molekuele durchgehen
				{
					sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex[z3]];
					o->m_pVDF[0]->BuildAtomList(NULL,sm,&templa);

					for (z4=0;z4<templa.GetSize();z4++)
						o->m_pVDF[0]->m_pVDF->AddToBin(g_pTempTimestep->m_vaVelocities[templa[z4]].GetLength());

					if (o->m_pVDF[0]->m_bSplitCart) {
						for (z4=0;z4<templa.GetSize();z4++) {
							o->m_pVDF[0]->m_pVDFSplit[0]->AddToBin( fabs(g_pTempTimestep->m_vaVelocities[templa[z4]][0]) );
							o->m_pVDF[0]->m_pVDFSplit[1]->AddToBin( fabs(g_pTempTimestep->m_vaVelocities[templa[z4]][1]) );
							o->m_pVDF[0]->m_pVDFSplit[2]->AddToBin( fabs(g_pTempTimestep->m_vaVelocities[templa[z4]][2]) );
							o->m_pVDF[0]->m_pVDFSplit[3]->AddToBin( sqrt(pow2(g_pTempTimestep->m_vaVelocities[templa[z4]][0])+pow2(g_pTempTimestep->m_vaVelocities[templa[z4]][1])) );
							o->m_pVDF[0]->m_pVDFSplit[4]->AddToBin( sqrt(pow2(g_pTempTimestep->m_vaVelocities[templa[z4]][0])+pow2(g_pTempTimestep->m_vaVelocities[templa[z4]][2])) );
							o->m_pVDF[0]->m_pVDFSplit[5]->AddToBin( sqrt(pow2(g_pTempTimestep->m_vaVelocities[templa[z4]][1])+pow2(g_pTempTimestep->m_vaVelocities[templa[z4]][2])) );
						}
					}
				}
			}
		}

		if ((!g_bCDF) && (!g_bRDF) && (!g_bSDF) && (!g_bPlProj) && (!g_bPlDF) && (!g_bLiDF) &&
           (!g_bADF) && (!g_bDDF) && (!g_bVDF) && (!g_bDipDF) && (!g_bCutCluster) && (!g_bSaveRefEnv) && (!g_bVHDF) && (!g_bRevSDF) && (!g_bCond))
			goto _endstep;


//		tfr = 9999.0;

		if (g_iFixMol == -1)
			goto _norefmol;

		for (z6=0;z6<g_oaObserv.GetSize();z6++)
		{
			o = (CObservation*)g_oaObserv[z6];

			if (o->m_pConditions != NULL)
				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();z2++)
					o->m_pConditions->m_iPassCounter[z2] = 0;

			if (o->m_pConditionsOM2 != NULL)
				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();z2++)
					o->m_pConditionsOM2->m_iPassCounter[z2] = 0;

			if (g_bCond) {
				if (o->m_bCondDevelopment) {
					for (z2=0;z2<=o->m_iCondDevelopmentMax;z2++)
						o->m_iaCondDevelopmentCounter[z2] = 0;
					o->m_iCondDevelopmentAvg = 0;
				}
			}

			if (o->m_bPercTimeDev) {
				o->m_iPercTimeDevTempCounter = 0;
				o->m_iPercTimeDevTempTotal = 0;
			}
		}


		if (g_bCDF) {
			for (z6=0;z6<g_oaObserv.GetSize();z6++) {
				o = (CObservation*)g_oaObserv[z6];
				if (o->m_pCDF != NULL)
					if (o->m_pCDF->m_bWriteSnapshots)
						o->m_pCDF->m_bWriteSnapshotDone = false;
			}
		}


		// Jedes Festhalte-Molekuel mal im Ursprung liegen lassen
		for (z2=0;z2<((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize();z2++)
		{
			if (g_bDoubleBox && (z2 >= ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()/g_iDoubleBoxFactor))
				continue;

			smfix = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[z2]];

			if (g_bCutCluster)
			{
				if (g_iClusterPos >= g_iClusterCount)
					break;

				if (((int)g_iSteps/g_iStride) >= g_iClusterSteps)
					break;

				if (g_iaClusterMol[g_iClusterPos] != z2)
					continue;

				if (g_iaClusterSteps[g_iClusterPos] != ((int)g_iSteps/g_iStride))
					break;

				g_iClusterPos++;
			}

			g_pTempTimestep->CopyFrom(GetTimeStep(0));

			if (g_bCDF) {
				for (z6=0;z6<g_oaObserv.GetSize();z6++) {
					o = (CObservation*)g_oaObserv[z6];
					if (o->m_pCDF != NULL)
						if (o->m_pCDF->m_bWriteSnapshots)
							o->m_pCDF->m_bWriteSnapshotsTempRM = false;
				}
			}

			// Zentrieren/Falten der ganzen Box nur noch noetig bei den folgenden Analysen
			if (g_bPlProj || g_bSDF || g_bRevSDF || g_bCutCluster || g_bSaveRefEnv || g_bSaveJustTraj ||
               g_bMiddleAvg)
			{
				vecc = g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[0]])->GetAt(g_iFixAtom[0])];

				if ((g_bCutCluster && !g_bRefEnvCenter) || (g_bSaveRefEnv && (!g_bRefEnvCenter) && (z2==g_iSaveRefMol)))
				{
					if (g_bFold)
					{
						g_pTempTimestep->CenterPos(vecc);

						if (g_bFoldAtomwise)
							g_pTempTimestep->FoldAtoms();
								else g_pTempTimestep->FoldMolecules();

						vec2 = -1.0 * vecc;

						g_pTempTimestep->CenterPos(vec2);
					}

					if (g_bCutCluster || (g_iNbhMode == 2))
					{
						g_pNbSet->Reset();

						if (g_bCutCluster)
							g_pNbSet->ResetAlwaysTrue();

						g_pNbSet->Scan(smfix,g_pTempTimestep);

						if (g_bSaveRefWithEnv)
							g_pNbSet->AddMolecule(g_iFixMol,z2);
					}

					g_pTempTimestep->WriteTimestepNb(g_fRefEnv,g_pNbSet);

					if (g_bTDO)
					{
						if (g_laTDOSteps.Contains(g_iSteps-1))
						{
//							sprintf(buf,"tdo_%s_%d_%06lu%s.xyz",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,g_iSteps-1,multibuf);
							buf.sprintf("tdo_%s_%d_%06lu%s.xyz",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,g_iSteps-1,multibuf);
							mprintf("\nSaving TDO step as %s.\n",(const char*)buf);
							tfi = OpenFileWrite(buf,true);
							g_pTempTimestep->WriteTimestepNb(tfi,g_pNbSet);
							fclose(tfi);
						}
					}
				}

				// Legt vecc genau in den Ursprung des Koordinatensystems
				g_pTempTimestep->CenterPos(vecc);


				if ((g_bCutCluster && g_bRefEnvCenter && (!g_bRefEnvFix)) || (g_bSaveRefEnv && g_bRefEnvCenter && (!g_bRefEnvFix) && (z2==g_iSaveRefMol)))
				{
					if (g_bFold)
					{
						if (g_bFoldAtomwise)
							g_pTempTimestep->FoldAtoms();
								else g_pTempTimestep->FoldMolecules();
					}

					if (g_bCutCluster || (g_iNbhMode == 2))
					{
						g_pNbSet->Reset();

						if (g_bCutCluster)
							g_pNbSet->ResetAlwaysTrue();

						g_pNbSet->Scan(smfix,g_pTempTimestep);

						if (g_bSaveRefWithEnv)
							g_pNbSet->AddMolecule(g_iFixMol,z2);
					}

					if (!g_bCenterZero && g_bPeriodic)
						g_pTempTimestep->CenterPos(CxDVector3(-g_fBoxX/2,-g_fBoxY/2,-g_fBoxZ/2));

					g_pTempTimestep->WriteTimestepNb(g_fRefEnv,g_pNbSet,((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[z2]);

					if (g_bTDO)
					{
						if (g_laTDOSteps.Contains(g_iSteps-1))
						{
//							sprintf(buf,"tdo_%s_%d_%06lu%s.xyz",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,g_iSteps-1,multibuf);
							buf.sprintf("tdo_%s_%d_%06lu%s.xyz",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,g_iSteps-1,multibuf);
							mprintf("\nSaving TDO step as %s.\n",(const char*)buf);
							tfi = OpenFileWrite(buf,true);
							g_pTempTimestep->WriteTimestepNb(tfi,g_pNbSet);
							fclose(tfi);
						}
					}
				}

				/* Fixed very bad SDF bug: Moved this block before the application of the transformation matrix */
				if (g_bFold)
				{
					if (g_bFoldAtomwise)
						g_pTempTimestep->FoldAtoms();
							else g_pTempTimestep->FoldMolecules();
				}
				/* End "this block" */

				if (g_bPlProj || g_bSDF ||
                   g_bAvg || ((g_bSaveRefEnv || g_bCutCluster) && g_bRefEnvCenter && g_bRefEnvFix)) // Wir brauchen ein neues Koordinatensystem und eine Drehmatrix
				{
					vec2 = g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[1]])->GetAt(g_iFixAtom[1])];
					vec3 = g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[2]])->GetAt(g_iFixAtom[2])];
					mat.MatUltra(vec2,vec3); // Dies erstellt uns die Drehmatrix

					if (g_bSDFVoro)
						dmat.MatUltra((CxDVector3)vec2,(CxDVector3)vec3); // Nochmal in Double Precision

					g_pTempTimestep->Transform(mat);

					if (g_bPlProj)
					{
						for (z6=0;z6<g_oaObserv.GetSize();z6++)
						{
							o = (CObservation*)g_oaObserv[z6];

							if (o->m_pPlProj->m_bAverageAtomPos)
							{
								o->m_pPlProj->m_iAverageCounter++;
								ti = 0;
								for (z3=0;z3<o->m_pPlProj->m_oDrawAtoms.m_oaAtoms.GetSize();z3++)
								{
									for (z4=0;z4<((CxIntArray*)o->m_pPlProj->m_oDrawAtoms.m_oaAtoms[z3])->GetSize();z4++)
									{
										o->m_pPlProj->m_vaAtomPos[ti] += g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[o->m_pPlProj->m_oDrawAtoms.m_baAtomType[z3]])->GetAt(((CxIntArray*)o->m_pPlProj->m_oDrawAtoms.m_oaAtoms[z3])->GetAt(z4))];
										ti++;
									}
								}
							}
						}
					}
				}

				if ((g_bCutCluster && g_bRefEnvCenter && g_bRefEnvFix) || (g_bSaveRefEnv && g_bRefEnvCenter && g_bRefEnvFix && (z2==g_iSaveRefMol)))
				{
					if (g_bCutCluster || (g_iNbhMode == 2))
					{
						g_pNbSet->Reset();
						if (g_bCutCluster)
							g_pNbSet->ResetAlwaysTrue();
						g_pNbSet->Scan(smfix,g_pTempTimestep);
						if (g_bSaveRefWithEnv)
							g_pNbSet->AddMolecule(g_iFixMol,z2);
					}

					if (!g_bCenterZero && g_bPeriodic)
						g_pTempTimestep->CenterPos(CxDVector3(-g_fBoxX/2,-g_fBoxY/2,-g_fBoxZ/2));

					g_pTempTimestep->WriteTimestepNb(g_fRefEnv,g_pNbSet,((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[z2]);

					if (g_bTDO)
					{
						if (g_laTDOSteps.Contains(g_iSteps-1))
						{
//							sprintf(buf,"tdo_%s_%d_%06lu%s.xyz",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,g_iSteps-1,multibuf);
							buf.sprintf("tdo_%s_%d_%06lu%s.xyz",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,g_iSteps-1,multibuf);
							mprintf("\nSaving TDO step as %s.\n",(const char*)buf);
							tfi = OpenFileWrite(buf,true);
							g_pTempTimestep->WriteTimestepNb(tfi,g_pNbSet);
							fclose(tfi);
						}
					}
				}
			} // END IF SDF, RevEnv, ...

/***************************************************************************************************************************/

			if (g_bSDF && g_bSDFMap)
			{
				for (z3=0;z3<g_oaSDFMaps.GetSize();z3++)
				{
					((CSDFMap*)g_oaSDFMaps[z3])->Process(z2,g_pTempTimestep);
				}
			}


			if (g_bPlProj || g_bDipDF || g_bRevSDF || g_bVHDF || g_bRDF || g_bPlDF || g_bLiDF || g_bSDF || g_bADF || g_bDDF || g_bVDF || g_bCond) 
			{
				for (z6=0;z6<g_oaObserv.GetSize();z6++)
				{
					o = (CObservation*)g_oaObserv[z6];

					if (o->m_bNormalizeCondition)
						o->m_iNormalizeConditionLocalCount = 0;

				//	mprintf("\nStep %d, z2=%d, Obs=%d, Reg=%d, ObsReg=%d.",g_iSteps,z2,z6,g_iaSMRegion[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[z2]],o->m_iaRMRegions[0]);

					if (g_bRegionAnalysis)
						if ((!o->m_iaRMRegions.Contains(0)) && (!o->m_iaRMRegions.Contains(g_iaSMRegion[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[z2]])))
							continue;

				//	mprintf("\n   ...RM ist in Region.");

					ti2 = o->m_waSaveRefList.GetPosition((unsigned short)z2);

					if (o->m_bObsCertain)
					{
				//		mprintf("\n    RefList: ");
				//		for (z7=0;z7<o->m_waObsRefList.GetSize();z7++)
				//			mprintf("%d, ",o->m_waObsRefList[z7]);

						if (!o->m_waObsRefList.Contains((unsigned short)z2))
							continue;
				//		mprintf("\n   ...RM ist auch in RefList.");
						tic_r = o->m_waObsRefList.GetPosition((unsigned short)z2);
					}

					if (g_bUseVelocities)
					{
						tempvel.CopyFrom(&g_pTempTimestep->m_vaVelocities);
						if (o->m_bVelocityRelToRef)
						{
							vecv = g_pTempTimestep->m_vaVelocities[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[0]])->GetAt(g_iFixAtom[0])];
							for (z3=0;z3<g_iGesVirtAtomCount;z3++)
								tempvel[z3] -= vecv;
						}
					}

					if ((g_bCDF) && (o->m_bTimeDev) && (ti2 != -1))
					{
						if (o->m_bSaveSeparateFiles)
							mfprintf(o->m_pCDF->m_fTimeDev[ti2],"%d",(int)g_iSteps);
								else if (ti2 == 0)
									mfprintf(o->m_pCDF->m_fTimeDev[0],"%d",(int)g_iSteps);
					}

					for (zr=0;zr<g_iCDFChannels;zr++)
					{
						if (g_bRDF && (o->m_pRDF[zr] != NULL))
						{
							for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
								o->m_pRDF[zr]->m_faData[z3].RemoveAll_KeepSize();
							if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
								for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
									o->m_pRDF[zr]->m_baDataEnabled[z3].RemoveAll_KeepSize();
							if ((o->m_bTimeDev) && (ti2 != -1) && (!g_bDeriv || (g_iSteps > 2)))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pRDF[zr]->m_fDist[ti2],"%d",(int)g_iSteps);
										else if (ti2 == 0)
											mfprintf(o->m_pRDF[zr]->m_fDist[0],"%d",(int)g_iSteps);
							}
						}
							
						if (g_bADF && (o->m_pADF[zr] != NULL))
						{
							for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
								o->m_pADF[zr]->m_faData[z3].RemoveAll_KeepSize();
							if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
								for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
									o->m_pADF[zr]->m_baDataEnabled[z3].RemoveAll_KeepSize();
							if ((o->m_bTimeDev) && (ti2 != -1) && (!g_bDeriv || (g_iSteps > 2)))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pADF[zr]->m_fAngle[ti2],"%d",(int)g_iSteps);
										else if (ti2 == 0)
											mfprintf(o->m_pADF[zr]->m_fAngle[0],"%d",(int)g_iSteps);
							}
						}

						if (g_bDDF && (o->m_pDDF[zr] != NULL))
						{
							for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
								o->m_pDDF[zr]->m_faData[z3].RemoveAll_KeepSize();
							if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
								for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
									o->m_pDDF[zr]->m_baDataEnabled[z3].RemoveAll_KeepSize();
							if ((o->m_bTimeDev) && (ti2 != -1) && (!g_bDeriv || (g_iSteps > 2)))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pDDF[zr]->m_fAngle[ti2],"%d",(int)g_iSteps);
										else if (ti2 == 0)
											mfprintf(o->m_pDDF[zr]->m_fAngle[0],"%d",(int)g_iSteps);
							}
						}

						if (g_bPlDF && (o->m_pPlDF[zr] != NULL))
						{
							for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
								o->m_pPlDF[zr]->m_faData[z3].RemoveAll_KeepSize();
							if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
								for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
									o->m_pPlDF[zr]->m_baDataEnabled[z3].RemoveAll_KeepSize();
/*							if ((o->m_bTimeDev) && (ti2 != -1) && (!g_bDeriv || (g_iSteps > 2)))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pPlDF[zr]->m_fAngle[ti2],"%d",(int)g_iSteps);
										else if (ti2 == 0)
											mfprintf(o->m_pDDF[zr]->m_fAngle[0],"%d",(int)g_iSteps);
							}*/
						}

						if (g_bLiDF && (o->m_pLiDF[zr] != NULL))
						{
							for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
								o->m_pLiDF[zr]->m_faData[z3].RemoveAll_KeepSize();
							if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
								for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
									o->m_pLiDF[zr]->m_baDataEnabled[z3].RemoveAll_KeepSize();
/*							if ((o->m_bTimeDev) && (ti2 != -1) && (!g_bDeriv || (g_iSteps > 2)))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pPlDF[zr]->m_fAngle[ti2],"%d",(int)g_iSteps);
										else if (ti2 == 0)
											mfprintf(o->m_pDDF[zr]->m_fAngle[0],"%d",(int)g_iSteps);
							}*/
						}

						if (g_bDipDF && (o->m_pDipDF[zr] != NULL))
						{
							for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
								o->m_pDipDF[zr]->m_faData[z3].RemoveAll_KeepSize();
							if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
								for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
									o->m_pDipDF[zr]->m_baDataEnabled[z3].RemoveAll_KeepSize();
							if ((o->m_bTimeDev) && (ti2 != -1) && (!g_bDeriv || (g_iSteps > 2)))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pDipDF[zr]->m_fDipole[ti2],"%d",(int)g_iSteps);
										else if (ti2 == 0)
											mfprintf(o->m_pDipDF[zr]->m_fDipole[0],"%d",(int)g_iSteps);
							}
						}

						if (g_bVDF && (o->m_pVDF[zr] != NULL))
						{
							for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
								o->m_pVDF[zr]->m_faData[z3].RemoveAll_KeepSize();
							if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
								for (z3=0;z3<((o->m_bSecondShowMol && (zr == 1))?o->m_iShowMol2Count:o->m_iShowMolCount);z3++)
									o->m_pVDF[zr]->m_baDataEnabled[z3].RemoveAll_KeepSize();
							if ((o->m_bTimeDev) && (ti2 != -1) && (!g_bDeriv || (g_iSteps > 2)))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pVDF[zr]->m_fSpeed[ti2],"%d",(int)g_iSteps);
										else if (ti2 == 0)
											mfprintf(o->m_pVDF[zr]->m_fSpeed[0],"%d",(int)g_iSteps);
							}
						}
					} // END FOR zr

					if (g_bSDF)
					{
						for (z3=0;z3<o->m_iShowMolCount;z3++)
						{
							o->m_pSDF->m_vaData[z3].RemoveAll_KeepSize();
							if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
								o->m_pSDF->m_baDataEnabled[z3].RemoveAll_KeepSize();
						}
					}

					if (g_bPlProj)
					{
						for (z3=0;z3<o->m_iShowMolCount;z3++)
						{
							o->m_pPlProj->m_vaData[z3].RemoveAll_KeepSize();
							if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
								o->m_pPlProj->m_baDataEnabled[z3].RemoveAll_KeepSize();
							if (o->m_pPlProj->m_bVector)
								o->m_pPlProj->m_vaVectorData[z3].RemoveAll_KeepSize();
						}
					}

					if (g_bRevSDF)
					{
						for (z3=0;z3<o->m_iShowMolCount;z3++)
						{
							o->m_pRevSDF->m_vaData[z3].RemoveAll_KeepSize();
							if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
								o->m_pRevSDF->m_baDataEnabled[z3].RemoveAll_KeepSize();
						}
					}

					if (o->m_bSelf)
					{
						if (g_bADF)
						{
							for (zr=0;zr<g_iCDFChannels;zr++)
							{
								if (o->m_pADF[zr] == NULL)
									continue;
								o->m_pADF[zr]->BuildAtomList(smfix,NULL,&templa);
								for (z4=0;z4<templa.GetSize();z4+=6)
								{
									if (o->m_pADF[zr]->m_bOrtho[0])
									{
										vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
										vec2 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
										vec3 = g_pTempTimestep->m_vaCoords[templa[z4+2]];
										vec1 = CrossP(FoldVector(vec2-vec0),FoldVector(vec3-vec0));
									} else
									{
										vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
										vec2 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
										vec1 = FoldVector(vec2-vec0);
									}
									if (o->m_pADF[zr]->m_bOrtho[1])
									{
										vec4 = g_pTempTimestep->m_vaCoords[templa[z4+3]];
										vec3 = g_pTempTimestep->m_vaCoords[templa[z4+4]];
										vec5 = g_pTempTimestep->m_vaCoords[templa[z4+5]];
										vec2 = CrossP(FoldVector(vec3-vec4),FoldVector(vec5-vec4));
									} else
									{
										vec4 = g_pTempTimestep->m_vaCoords[templa[z4+3]];
										vec3 = g_pTempTimestep->m_vaCoords[templa[z4+4]];
										vec2 = FoldVector(vec3-vec4);
									}
									tf = Angle_Deg(vec1,vec2);
									if ((tf > 90.0) && (o->m_pADF[zr]->m_bFoldAngle))
										tf = 180.0-tf;
									if (o->m_pADF[zr]->m_bCosine)
										tf = cos(tf/180.0*Pi);

									if (g_bDeriv)
									{
										o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pADF[zr]->m_iCombinations+z4/6) = tf;
										switch(o->m_pADF[zr]->m_iDeriv)
										{
											case 0:
												tf = o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(z2*o->m_pADF[zr]->m_iCombinations+z4/6);
												break;
											case 1:
												tf = (o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pADF[zr]->m_iCombinations+z4/6) - o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(z2*o->m_pADF[zr]->m_iCombinations+z4/6)) / (2*g_fTimestepLength);
												break;
											case 2:
												tf = (o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(z2*o->m_pADF[zr]->m_iCombinations+z4/6) + o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pADF[zr]->m_iCombinations+z4/6) - 2*o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(z2*o->m_pADF[zr]->m_iCombinations+z4/6)) / (g_fTimestepLength * g_fTimestepLength);
												break;
										}
									}

									if (!g_bDeriv || (g_iSteps > 2))
									{
										if (o->m_bObsCertain && o->m_bDecompDist)
											o->m_pADF[zr]->m_pADF->AddToBin_Multi(tic_r,tf);

										o->m_pADF[zr]->m_faData[0].Add(tf);

										if (o->m_pADF[zr]->m_bACF)
											o->m_pADF[zr]->m_pfaACFBuffer[z2*o->m_pADF[zr]->m_iCombinations+z4/6]->Add((double)tf);

					/*					if (o->m_pADF[zr]->m_bMirror)
										{
											if (o->m_bObsCertain && o->m_bDecompDist)
												o->m_pADF[zr]->m_pADF->AddToBin_Multi(tic_r,180.0-tf);
											o->m_pADF[zr]->m_faData[0].Add(180.0-tf);
										}*/

										if (o->m_bTimeDev && (ti2 != -1))
										{
											if (o->m_bSaveSeparateFiles)
												mfprintf(o->m_pADF[zr]->m_fAngle[ti2],"; %8.3f",tf);
													else mfprintf(o->m_pADF[zr]->m_fAngle[0],"; %8.3f",tf);
											if (o->m_bCombinedPlot)
												o->m_pADF[zr]->m_pADF->m_pCombinedPlot->AddXYTupel(ti2*o->m_pADF[zr]->m_iCombinations+z4/6,g_iSteps*g_fTimestepLength/1000.0,tf);
										}
										if (o->m_bTimeDiff)
											((CxDoubleArray*)o->m_pADF[zr]->m_pADF->m_oaTimeDiffBuf[z2*o->m_pADF[zr]->m_iCombinations+z4/6])->Add((double)tf);

										if (o->m_bPercTimeDev) {
											o->m_iPercTimeDevTempTotal++;
											if ((tf >= o->m_fPercTimeDevMin) && (tf <= o->m_fPercTimeDevMax))
												o->m_iPercTimeDevTempCounter++;
										}
									}
								} // END FOR z4
							} // END FOR zr
						} // Ende IF ADF

						if (g_bDDF)
						{
							for (zr=0;zr<g_iCDFChannels;zr++)
							{
								if (o->m_pDDF[zr] == NULL)
									continue;
								o->m_pDDF[zr]->BuildAtomList(smfix,smfix,&templa);
								for (z4=0;z4<templa.GetSize();z4+=9)
								{
									if (o->m_pDDF[zr]->m_bOrtho[0])
									{
										vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
										vec2 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
										vec3 = g_pTempTimestep->m_vaCoords[templa[z4+2]];
										vec1 = CrossP(FoldVector(vec2-vec0),FoldVector(vec3-vec0));
									} else
									{
										vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
										vec2 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
										vec1 = FoldVector(vec2-vec0);
									}
									if (o->m_pDDF[zr]->m_bOrtho[1])
									{
										vec4 = g_pTempTimestep->m_vaCoords[templa[z4+3]];
										vec3 = g_pTempTimestep->m_vaCoords[templa[z4+4]];
										vec5 = g_pTempTimestep->m_vaCoords[templa[z4+5]];
										vec2 = CrossP(FoldVector(vec3-vec4),FoldVector(vec5-vec4));
									} else
									{
										vec4 = g_pTempTimestep->m_vaCoords[templa[z4+3]];
										vec3 = g_pTempTimestep->m_vaCoords[templa[z4+4]];
										vec2 = FoldVector(vec3-vec4);
									}
									if (o->m_pDDF[zr]->m_bOrtho[2])
									{
										vec4 = g_pTempTimestep->m_vaCoords[templa[z4+6]];
										vec3 = g_pTempTimestep->m_vaCoords[templa[z4+7]];
										vec5 = g_pTempTimestep->m_vaCoords[templa[z4+8]];
										vec3 = CrossP(FoldVector(vec3-vec4),FoldVector(vec5-vec4));
									} else
									{
										vec4 = g_pTempTimestep->m_vaCoords[templa[z4+6]];
										vec3 = g_pTempTimestep->m_vaCoords[templa[z4+7]];
										vec3 = FoldVector(vec3-vec4);
									}
									tf = Dihedral(vec1,vec2,vec3,o->m_pDDF[zr]->m_bAbs);
									if (o->m_pDDF[zr]->m_bCosine)
										tf = cos(tf/180.0*Pi);
											else if (o->m_pDDF[zr]->m_bPositive)
												if (tf < 0)
													tf += 360.0;
									if (g_bDeriv)
									{
										o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pDDF[zr]->m_iCombinations+z4/9) = tf;
										switch(o->m_pDDF[zr]->m_iDeriv)
										{
											case 0:
												tf = o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(z2*o->m_pDDF[zr]->m_iCombinations+z4/9);
												break;
											case 1:
												tf = (o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pDDF[zr]->m_iCombinations+z4/9) - o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(z2*o->m_pDDF[zr]->m_iCombinations+z4/9)) / (2*g_fTimestepLength);
												break;
											case 2:
												tf = (o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(z2*o->m_pDDF[zr]->m_iCombinations+z4/9) + o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pDDF[zr]->m_iCombinations+z4/9) - 2*o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(z2*o->m_pDDF[zr]->m_iCombinations+z4/9)) / (g_fTimestepLength * g_fTimestepLength);
												break;
										}
									}
									if (!g_bDeriv || (g_iSteps > 2))
									{
										if (o->m_bObsCertain && o->m_bDecompDist)
											o->m_pDDF[zr]->m_pDDF->AddToBin_Multi(tic_r,tf);

										o->m_pDDF[zr]->m_faData[0].Add(tf);

										if (o->m_pDDF[zr]->m_bACF)
											o->m_pDDF[zr]->m_pfaACFBuffer[z2*o->m_pDDF[zr]->m_iCombinations+z4/9]->Add((double)tf);

										if (/*(!o->m_pDDF[zr]->m_bAbs) && */o->m_pDDF[zr]->m_bSymm)
										{
											o->m_pDDF[zr]->m_faData[0].Add(-tf);
											if (o->m_bObsCertain && o->m_bDecompDist)
												o->m_pDDF[zr]->m_pDDF->AddToBin_Multi(tic_r,-tf);
										}

										if (o->m_bTimeDiff)
											((CxDoubleArray*)o->m_pDDF[zr]->m_pDDF->m_oaTimeDiffBuf[z2*o->m_pDDF[zr]->m_iCombinations+z4/9])->Add((double)tf);

										if (o->m_bTimeDev && (ti2 != -1))
										{
											if (o->m_pDDF[zr]->m_bRotate)
											{
												if (g_iSteps > 1)
												{
													if ((tf - o->m_pDDF[zr]->m_faLastData[z2*o->m_pDDF[zr]->m_iCombinations+z4/9]) > 180.0)
														o->m_pDDF[zr]->m_laRotation[z2*o->m_pDDF[zr]->m_iCombinations+z4/9]--;
													if ((tf - o->m_pDDF[zr]->m_faLastData[z2*o->m_pDDF[zr]->m_iCombinations+z4/9]) < -180.0)
														o->m_pDDF[zr]->m_laRotation[z2*o->m_pDDF[zr]->m_iCombinations+z4/9]++;
												}
												o->m_pDDF[zr]->m_faLastData[z2*o->m_pDDF[zr]->m_iCombinations+z4/9] = tf;
												tf2 = tf + o->m_pDDF[zr]->m_laRotation[z2*o->m_pDDF[zr]->m_iCombinations+z4/9] * 360.0;
											} else tf2 = tf;
											if (o->m_bSaveSeparateFiles)
												mfprintf(o->m_pDDF[zr]->m_fAngle[ti2],"; %8.3f",tf2);
													else mfprintf(o->m_pDDF[zr]->m_fAngle[0],"; %8.3f",tf2);
											if (o->m_bCombinedPlot)
												o->m_pDDF[zr]->m_pDDF->m_pCombinedPlot->AddXYTupel(ti2*o->m_pDDF[zr]->m_iCombinations+z4/9,g_iSteps*g_fTimestepLength/1000.0,tf);
										}

										if (o->m_bPercTimeDev) {
											o->m_iPercTimeDevTempTotal++;
									//		mprintf("@ %f",tf);
											if ((tf >= o->m_fPercTimeDevMin) && (tf <= o->m_fPercTimeDevMax)) {
												o->m_iPercTimeDevTempCounter++;
									//			mprintf(" Ok");
											}
									//		mprintf("\n");
										}
									}
								}
							}
						} // END IF DDF

						if (g_bVHDF)
						{
							o->m_pVHDF->BuildAtomList(smfix,smfix,&templa);
							for (z4=0;z4<templa.GetSize();z4+=2)
							{
								vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];// + vecc;
								for (z0=0;z0<o->m_pVHDF->m_iDepth/o->m_pVHDF->m_iStride;z0++)
								{
									g_pT2Timestep = GetTimeStep(z0*o->m_pVHDF->m_iStride);
									if (g_pT2Timestep == NULL)
										continue;
									vec1 = g_pT2Timestep->m_vaCoords[templa[z4+1]];
									tf = FoldedLength(vec0-vec1);
									o->m_pVHDF->m_pVHDF->AddToBin_IntX_fast(z0,tf);
									o->m_pVHDF->m_pCount[z0]++;
								}
							} 
						} // Ende IF VHDF

						if (g_bRDF)
						{
							for (zr=0;zr<g_iCDFChannels;zr++)
							{
								if (o->m_pRDF[zr] == NULL)
									continue;
								o->m_pRDF[zr]->BuildAtomList(smfix,smfix,&templa);
								for (z4=0;z4<templa.GetSize();z4+=2)
								{
									vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
									vec1 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
									tf = FoldedLength(vec0-vec1);
									if (g_bDeriv)
									{
										o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pRDF[zr]->m_iCombinations+z4/2) = tf;
										switch(o->m_pRDF[zr]->m_iDeriv)
										{
											case 0:
												tf = o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(z2*o->m_pRDF[zr]->m_iCombinations+z4/2);
												break;
											case 1:
												tf = (o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pRDF[zr]->m_iCombinations+z4/2) - o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(z2*o->m_pRDF[zr]->m_iCombinations+z4/2)) / (2*g_fTimestepLength);
												break;
											case 2:
												tf = (o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(z2*o->m_pRDF[zr]->m_iCombinations+z4/2) + o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pRDF[zr]->m_iCombinations+z4/2) - 2*o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(z2*o->m_pRDF[zr]->m_iCombinations+z4/2)) / (g_fTimestepLength * g_fTimestepLength);
												break;
										}
									}
									if (!g_bDeriv || (g_iSteps > 2))
									{
										if (o->m_bObsCertain && o->m_bDecompDist)
											o->m_pRDF[zr]->m_pRDF->AddToBin_Multi(tic_r,tf);

										if (o->m_bDecompType)
											o->m_pRDF[zr]->m_pRDF->AddToBin_Multi(o->m_waDecompTypeRefOffs[templa[z4]]*o->m_waDecompTypeObsOffs.GetSize()+o->m_waDecompTypeObsOffs[templa[z4+1]],tf);

										o->m_pRDF[zr]->m_faData[0].Add(tf);

										if (o->m_pRDF[zr]->m_bACF)
											o->m_pRDF[zr]->m_pfaACFBuffer[z2*o->m_pRDF[zr]->m_iCombinations+z4/2]->Add((double)tf);

										if (o->m_bTimeDev && (ti2 != -1))
										{
											if (o->m_bSaveSeparateFiles)
												mfprintf(o->m_pRDF[zr]->m_fDist[ti2],"; %10.3f",tf);
													else mfprintf(o->m_pRDF[zr]->m_fDist[0],"; %10.3f",tf);
											if (o->m_bCombinedPlot)
												o->m_pRDF[zr]->m_pRDF->m_pCombinedPlot->AddXYTupel(ti2*o->m_pRDF[zr]->m_iCombinations+z4/2,g_iSteps*g_fTimestepLength/1000.0,tf);
										}
										if (o->m_bTimeDiff)
											((CxDoubleArray*)o->m_pRDF[zr]->m_pRDF->m_oaTimeDiffBuf[z2*o->m_pRDF[zr]->m_iCombinations+z4/2])->Add((double)tf);

										if (o->m_bPercTimeDev) {
											o->m_iPercTimeDevTempTotal++;
											if ((tf >= o->m_fPercTimeDevMin) && (tf <= o->m_fPercTimeDevMax))
												o->m_iPercTimeDevTempCounter++;
										}
									}
								} 
							}
						} // Ende IF RDF

						if (g_bPlDF)
						{
							for (zr=0;zr<g_iCDFChannels;zr++)
							{
								if (o->m_pPlDF[zr] == NULL)
									continue;
								o->m_pPlDF[zr]->BuildAtomList(smfix,smfix,&templa);
								for (z4=0;z4<templa.GetSize();z4+=4)
								{
									if (o->m_pPlDF[zr]->m_bNormal)
									{
										vec0 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+1]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
										vec0.Normalize();
										vec1 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+3]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
										tf = DotP(vec0,vec1);
									} else
									{
										vec0 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+1]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
										vec1 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+2]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
										vec3 = CrossP(vec0,vec1);
										vec3.Normalize();
										vec2 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+3]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
										tf = DotP(vec2,vec3);
									}

									if (!g_bDeriv || (g_iSteps > 2))
									{
										if (o->m_bObsCertain && o->m_bDecompDist)
											o->m_pPlDF[zr]->m_pPlDF->AddToBin_Multi(tic_r,tf);

										o->m_pPlDF[zr]->m_faData[0].Add(tf);
									}
								}
							}
						} // END IF PlDF

						if (g_bLiDF)
						{
							for (zr=0;zr<g_iCDFChannels;zr++)
							{
								if (o->m_pLiDF[zr] == NULL)
									continue;
								o->m_pLiDF[zr]->BuildAtomList(smfix,smfix,&templa);
								for (z4=0;z4<templa.GetSize();z4+=4)
								{
									if (o->m_pLiDF[zr]->m_bNormal)
									{
										vec2 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+1]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
										vec3 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+2]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
										vec0 = CrossP(vec2,vec3);
										vec0.Normalize();
										vec1 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+3]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
										vec0 *= DotP(vec0,vec1);
										tf = (vec1 - vec0).GetLength();
									} else
									{
										vec0 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+1]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
										vec0.Normalize();
										vec1 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+3]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
										vec0 *= DotP(vec0,vec1);
										tf = (vec1 - vec0).GetLength();
									}

									if (!g_bDeriv || (g_iSteps > 2))
									{
										if (o->m_bObsCertain && o->m_bDecompDist)
											o->m_pLiDF[zr]->m_pLiDF->AddToBin_Multi(tic_r,tf);

										o->m_pLiDF[zr]->m_faData[0].Add(tf);
									}
								}
							}
						} // END IF LiDF

						if (g_bDipDF)
						{
							for (zr=0;zr<g_iCDFChannels;zr++)
							{
								if (o->m_pDipDF[zr] == NULL)
									continue;
								tf = smfix->m_vDipole.GetLength();
								if (g_bDeriv)
								{
									o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2) = tf;
									switch(o->m_pDipDF[zr]->m_iDeriv)
									{
										case 0:
											tf = o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(z2);
											break;
										case 1:
											tf = (o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2) - o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(z2)) / (2*g_fTimestepLength);
											break;
										case 2:
											tf = (o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(z2) + o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2) - 2*o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(z2)) / (g_fTimestepLength * g_fTimestepLength);
											break;
									}
								}
								if (!g_bDeriv || (g_iSteps > 2))
								{
									if (o->m_bObsCertain && o->m_bDecompDist)
										o->m_pDipDF[zr]->m_pDipoleDF->AddToBin_Multi(tic_r,tf);

									o->m_pDipDF[zr]->m_faData[0].Add(tf);

									if (o->m_pDipDF[zr]->m_bACF)
										o->m_pDipDF[zr]->m_pfaACFBuffer[z2]->Add((double)tf);

									if (o->m_bTimeDev && (ti2 != -1))
									{
										if (o->m_bSaveSeparateFiles)
											mfprintf(o->m_pDipDF[zr]->m_fDipole[ti2],"; %8.3f",tf);
												else mfprintf(o->m_pDipDF[zr]->m_fDipole[0],"; %8.3f",tf);
										if (o->m_bCombinedPlot)
											o->m_pDipDF[zr]->m_pDipoleDF->m_pCombinedPlot->AddXYTupel(ti2,g_iSteps*g_fTimestepLength/1000.0,tf);
									}
									if (o->m_bTimeDiff)
										((CxDoubleArray*)o->m_pDipDF[zr]->m_pDipoleDF->m_oaTimeDiffBuf[z2])->Add((double)tf);
								}
							}
						} // Ende IF DIPDF

						if (g_bVDF)
						{
							for (zr=0;zr<g_iCDFChannels;zr++)
							{
								if (o->m_pVDF[zr] == NULL)
									continue;
								o->m_pVDF[zr]->BuildAtomList(smfix,smfix,&templa);
								for (z4=0;z4<templa.GetSize();z4++)
								{
									tf = tempvel[templa[z4]].GetLength();
									if (g_bDeriv)
									{
										o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pVDF[zr]->m_iCombinations+z4/2) = tf;
										switch(o->m_pVDF[zr]->m_iDeriv)
										{
											case 0:
												tf = o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(z2*o->m_pVDF[zr]->m_iCombinations+z4/2);
												break;
											case 1:
												tf = (o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pVDF[zr]->m_iCombinations+z4/2) - o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(z2*o->m_pVDF[zr]->m_iCombinations+z4/2)) / (2*g_fTimestepLength);
												break;
											case 2:
												tf = (o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(z2*o->m_pVDF[zr]->m_iCombinations+z4/2) + o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(z2*o->m_pVDF[zr]->m_iCombinations+z4/2) - 2*o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(z2*o->m_pVDF[zr]->m_iCombinations+z4/2)) / (g_fTimestepLength * g_fTimestepLength);
												break;
										}
									}
									if (!g_bDeriv || (g_iSteps > 2))
									{
										if (o->m_bObsCertain && o->m_bDecompDist)
											o->m_pVDF[zr]->m_pVDF->AddToBin_Multi(tic_r,tf);

										o->m_pVDF[zr]->m_faData[0].Add(tf);

										if (o->m_pVDF[zr]->m_bACF)
											o->m_pVDF[zr]->m_pfaACFBuffer[z2*o->m_pVDF[zr]->m_iCombinations+z4]->Add((double)tf);

										if (o->m_bTimeDev && (ti2 != -1))
										{
											if (o->m_bSaveSeparateFiles)
												mfprintf(o->m_pVDF[zr]->m_fSpeed[ti2],"; %8.3f",tf);
													else mfprintf(o->m_pVDF[zr]->m_fSpeed[0],"; %8.3f",tf);
											if (o->m_bCombinedPlot)
												o->m_pVDF[zr]->m_pVDF->m_pCombinedPlot->AddXYTupel(ti2*o->m_pVDF[zr]->m_iCombinations+z4,g_iSteps*g_fTimestepLength/1000.0,tf);
										}
										if (o->m_bTimeDiff)
											((CxDoubleArray*)o->m_pVDF[zr]->m_pVDF->m_oaTimeDiffBuf[z2*o->m_pVDF[zr]->m_iCombinations+z4])->Add((double)tf);
									}
								}
							}
						} // Ende IF VDF

						if (g_bSDF)
						{
							o->m_pSDF->BuildAtomList(smfix,smfix,&templa);
							for (z4=0;z4<templa.GetSize();z4++)
							{
								{
									o->m_pSDF->m_vaData[0].Add(g_pTempTimestep->m_vaCoords[templa[z4]]);

									if (o->m_pSDF->m_bVdWSpheres)
										o->m_pSDF->m_faRadius[0].Add(g_faVdWRadius[templa[z4]]);
								}
							} 
						} // Ende IF SDF

						if (g_bPlProj)
						{
							o->m_pPlProj->BuildAtomList(smfix,smfix,&templa);
							for (z4=0;z4<templa.GetSize();z4++)
								o->m_pPlProj->m_vaData[0].Add(g_pTempTimestep->m_vaCoords[templa[z4]]);
						} // Ende IF PlProj

						if (g_bRevSDF)
						{
							o->m_pRevSDF->BuildAtomList(smfix,smfix,&templa);
							for (z4=0;z4<templa.GetSize();z4++)
								o->m_pRevSDF->m_vaData[0].Add(g_pTempTimestep->m_vaCoords[templa[z4]]);
						} // Ende IF RevSDF

					} // Ende IF m_bSelf  

					if (o->m_bOthers)
					{
						secondmolrun = false;
_secondmolstart:
						if ((!secondmolrun) && (o->m_pConditions != NULL))
						{
							o->m_pConditions->m_bAnyPassed = false;
							o->m_pConditions->m_iTempPassed = 0;
							o->m_pConditions->ScanNeighborhoodAllOM(g_pTempTimestep,smfix);
							if (o->m_pConditions->m_bAnyPassed)
								o->m_pConditions->m_iRMPassCounter[z2]++;
						}

						if ((secondmolrun) && (o->m_pConditionsOM2 != NULL))
						{
							o->m_pConditionsOM2->m_bAnyPassed = false;
							o->m_pConditionsOM2->m_iTempPassed = 0;
							o->m_pConditionsOM2->ScanNeighborhoodAllOM(g_pTempTimestep,smfix);
							if (o->m_pConditionsOM2->m_bAnyPassed)
								o->m_pConditionsOM2->m_iRMPassCounter[z2]++;
						}

						for (z3=0;z3<(secondmolrun?o->m_iShowMol2Count:o->m_iShowMolCount);z3++) // Alle anderen Molekuele durchgehen
						{
							if (secondmolrun)
							{
								if ((g_iFixMol == o->m_iShowMol2) && (z3 == z2)) // Wir wollen nicht das Referenzmolekuel mitzaehlen
										continue;
								if (o->m_bObsCertain)
								{
									if (!o->m_waObsShow2List.Contains((unsigned short)z3))
										continue;
									tic_o = o->m_waObsShow2List.GetPosition((unsigned short)z3);
								}
							} else
							{
								if ((g_iFixMol == o->m_iShowMol) && (z3 == z2)) // Wir wollen nicht das Referenzmolekuel mitzaehlen
										continue;
								if (o->m_bObsCertain)
								{
									if (!o->m_waObsShowList.Contains((unsigned short)z3))
										continue;
									tic_o = o->m_waObsShowList.GetPosition((unsigned short)z3);
								}
							}

							if (g_bRegionAnalysis)
							{
								if (secondmolrun)
								{
									if ((!o->m_iaOM2Regions.Contains(0)) && (!o->m_iaOM2Regions.Contains(g_iaSMRegion[((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_laSingleMolIndex[z3]])))
										continue;
								} else
								{
									if ((!o->m_iaOM1Regions.Contains(0)) && (!o->m_iaOM1Regions.Contains(g_iaSMRegion[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex[z3]])))
										continue;
								}
							}
						//	mprintf("\n   ... OM %d ist drin!",z3);

							if (secondmolrun)
							{
								ti = -1;
								sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_laSingleMolIndex[z3]];
							} else
							{
								ti = o->m_waSaveShowList.GetPosition((unsigned short)z3);
								sm = (CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex[z3]];
							}

							if (o->m_bTimeDev && (o->m_pConditions != NULL))
							{
								if (!o->m_pConditions->Contains(z3) && (!g_bDeriv || (g_iSteps > 2)))
								{
									for (zr=0;zr<g_iCDFChannels;zr++)
									{
										if (g_bADF && (o->m_pADF[zr] != NULL))
										{
											if ((ti2 != -1) && (ti != -1))
											{
												if (o->m_bSaveSeparateFiles)
													mfprintf(o->m_pADF[zr]->m_fAngle[ti2],";     -   ");
														else mfprintf(o->m_pADF[zr]->m_fAngle[0],";     -   ");
											}
										}

										if (g_bDDF && (o->m_pDDF[zr] != NULL))
										{
											if ((ti2 != -1) && (ti != -1))
											{
												if (o->m_bSaveSeparateFiles)
													mfprintf(o->m_pADF[zr]->m_fAngle[ti2],";     -   ");
														else mfprintf(o->m_pDDF[zr]->m_fAngle[0],";     -   ");
											}
										}

										if (g_bDipDF && (o->m_pDipDF[zr] != NULL))
										{
											if ((ti2 != -1) && (ti != -1))
											{
												if (o->m_bSaveSeparateFiles)
													mfprintf(o->m_pDipDF[zr]->m_fDipole[ti2],";     -   ");
														else mfprintf(o->m_pDipDF[zr]->m_fDipole[0],";     -   ");
											}
										}

										if (g_bVDF && (o->m_pVDF[zr] != NULL))
										{
											if ((ti2 != -1) && (ti != -1))
											{
												if (o->m_bSaveSeparateFiles)
													mfprintf(o->m_pVDF[zr]->m_fSpeed[ti2],";     -   ");
														else mfprintf(o->m_pVDF[zr]->m_fSpeed[0],";     -   ");
											}
										}

										if (g_bRDF && (o->m_pRDF[zr] != NULL))
										{
											if ((ti2 != -1) && (ti != -1))
											{
												if (o->m_bSaveSeparateFiles)
													mfprintf(o->m_pRDF[zr]->m_fDist[ti2],";     -   ");
														else mfprintf(o->m_pRDF[zr]->m_fDist[0],";     -   ");
											}
										}
									}
									continue;
								}
							}

							if ((!secondmolrun) && (o->m_pConditions != NULL))
							{
								if (!o->m_pConditions->Contains(z3))
									continue;

								if (o->m_bBinOnlyPassedAtoms)
									o->m_pConditions->MarkPassedAtoms(z3,true);

								if (o->m_bBinOnlyNotPassedAtoms)
									o->m_pConditions->MarkPassedAtoms(z3,false);

								if (g_bSaveCondSnapshot)
								{
									g_iSaveCondCount++;
									if (g_bSaveCondWholeBox)
									{
										g_pTempTimestep->WriteTimestep(g_fSaveCondFile);
									} else
									{
										mfprintf(g_fSaveCondFile,"%d\nTimestep %lu, RM %d, OM %d\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes - ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laVirtualAtoms.GetSize() + ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_iAtomGes - ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laVirtualAtoms.GetSize(),g_iSteps,z2+1,z3+1);
										for (z4=0;z4<((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex.GetSize();z4++)
										{
											if ((!g_bSaveVirtAtoms) && (((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[z4] == g_iVirtAtomType))
												continue;
											for (z5=0;z5<((CxIntArray*)smfix->m_oaAtomOffset[z4])->GetSize();z5++)
												mfprintf(g_fSaveCondFile,"%s  %f  %f  %f\n",(const char*)((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[z4]])->m_sName,g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[z4])->GetAt(z5)][0]/100.0,g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[z4])->GetAt(z5)][1]/100.0,g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[z4])->GetAt(z5)][2]/100.0);
										}
										for (z4=0;z4<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_baAtomIndex.GetSize();z4++)
										{
											if ((!g_bSaveVirtAtoms) && (((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_baAtomIndex[z4] == g_iVirtAtomType))
												continue;
											for (z5=0;z5<((CxIntArray*)sm->m_oaAtomOffset[z4])->GetSize();z5++)
												mfprintf(g_fSaveCondFile,"%s  %f  %f  %f\n",(const char*)((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_baAtomIndex[z4]])->m_sName,g_pTempTimestep->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z4])->GetAt(z5)][0]/100.0,g_pTempTimestep->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z4])->GetAt(z5)][1]/100.0,g_pTempTimestep->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z4])->GetAt(z5)][2]/100.0);
										}
									}
								}
							}

							if (secondmolrun && (o->m_pConditionsOM2 != NULL))
							{
								if (!o->m_pConditionsOM2->Contains(z3))
									continue;

					//			mprintf("\nStep %d: z2=%d, z3=%d passed.",g_iSteps,z2,z3);

								if (o->m_bBinOnlyPassedAtoms)
									o->m_pConditionsOM2->MarkPassedAtoms(z3,true);

								if (o->m_bBinOnlyNotPassedAtoms)
									o->m_pConditionsOM2->MarkPassedAtoms(z3,false);

						/*		if (g_bSaveCondSnapshot)
								{
									g_iSaveCondCount++;
									if (g_bSaveCondWholeBox)
									{
										g_pTempTimestep->WriteTimestep(g_fSaveCondFile);
									} else
									{
										mfprintf(g_fSaveCondFile,"%d\nTimestep %lu, RM %d, OM %d\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes - ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laVirtualAtoms.GetSize() + ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_iAtomGes - ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laVirtualAtoms.GetSize(),g_iSteps,z2+1,z3+1);
										for (z4=0;z4<((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex.GetSize();z4++)
										{
											if ((!g_bSaveVirtAtoms) && (((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[z4] == g_iVirtAtomType))
												continue;
											for (z5=0;z5<((CxIntArray*)smfix->m_oaAtomOffset[z4])->GetSize();z5++)
												mfprintf(g_fSaveCondFile,"%s  %f  %f  %f\n",((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[z4]])->m_sName,g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[z4])->GetAt(z5)][0]/100.0,g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[z4])->GetAt(z5)][1]/100.0,g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[z4])->GetAt(z5)][2]/100.0);
										}
										for (z4=0;z4<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_baAtomIndex.GetSize();z4++)
										{
											if ((!g_bSaveVirtAtoms) && (((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_baAtomIndex[z4] == g_iVirtAtomType))
												continue;
											for (z5=0;z5<((CxIntArray*)sm->m_oaAtomOffset[z4])->GetSize();z5++)
												mfprintf(g_fSaveCondFile,"%s  %f  %f  %f\n",((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_baAtomIndex[z4]])->m_sName,g_pTempTimestep->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z4])->GetAt(z5)][0]/100.0,g_pTempTimestep->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z4])->GetAt(z5)][1]/100.0,g_pTempTimestep->m_vaCoords[((CxIntArray*)sm->m_oaAtomOffset[z4])->GetAt(z5)][2]/100.0);
										}
									}
								}*/
							}

							if (o->m_bNormalizeCondition)
								o->m_iNormalizeConditionLocalCount++;

							if (g_bVHDF)
							{
								o->m_pVHDF->BuildAtomList(smfix,sm,&templa);
								for (z4=0;z4<templa.GetSize();z4+=2)
								{
									vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];// + vecc;
									for (z0=0;z0<o->m_pVHDF->m_iDepth/o->m_pVHDF->m_iStride;z0++)
									{
										g_pT2Timestep = GetTimeStep(z0*o->m_pVHDF->m_iStride);
										if (g_pT2Timestep == NULL)
											continue;
										vec1 = g_pT2Timestep->m_vaCoords[templa[z4+1]];
										tf = FoldedLength(vec0-vec1);
										o->m_pVHDF->m_pVHDF->AddToBin_IntX_fast(z0,tf);
										o->m_pVHDF->m_pCount[z0]++;
									}
								} 
							} // Ende IF VHDF

							if (g_bADF)
							{
								for (zr=0;zr<g_iCDFChannels;zr++)
								{
									if (o->m_pADF[zr] == NULL)
										continue;
									if (o->m_bSecondShowMol)
									{
										if (secondmolrun && (zr == 0))
											continue;
										if (!secondmolrun && (zr == 1))
											continue;
									}
									o->m_pADF[zr]->BuildAtomList(smfix,sm,&templa);
									for (z4=0;z4<templa.GetSize();z4+=6)
									{
										if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
										{
											tb = true;
											if (g_baAtomPassedCondition[templa[z4]] == 0)
												tb = false;
											if (g_baAtomPassedCondition[templa[z4+1]] == 0)
												tb = false;
											if (o->m_pADF[zr]->m_bOrtho[0])
												if (g_baAtomPassedCondition[templa[z4+2]] == 0)
													tb = false;
											if (g_baAtomPassedCondition[templa[z4+3]] == 0)
												tb = false;
											if (g_baAtomPassedCondition[templa[z4+4]] == 0)
												tb = false;
											if (o->m_pADF[zr]->m_bOrtho[1])
												if (g_baAtomPassedCondition[templa[z4+5]] == 0)
													tb = false;
											if (tb)
											{
												o->m_pADF[zr]->m_baDataEnabled[z3].Add(1);
											} else
											{
									//			mprintf("Nicht bestanden; z2=%d, z3=%d; (%d,%d,%d) (%d,%d,%d).\n",z2,z3,g_baAtomPassedCondition[tempwa[z4]],g_baAtomPassedCondition[tempwa[z4+1]],g_baAtomPassedCondition[tempwa[z4+2]],g_baAtomPassedCondition[tempwa[z4+3]],g_baAtomPassedCondition[tempwa[z4+4]],g_baAtomPassedCondition[tempwa[z4+5]]);
												o->m_pADF[zr]->m_baDataEnabled[z3].Add(0);
											}
										}
										if (o->m_pADF[zr]->m_bOrtho[0])
										{
											vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
											vec2 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
											vec3 = g_pTempTimestep->m_vaCoords[templa[z4+2]];
											vec1 = CrossP(FoldVector(vec2-vec0),FoldVector(vec3-vec0));
										} else
										{
											vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
											vec2 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
											vec1 = FoldVector(vec2-vec0);
										}
										if (o->m_pADF[zr]->m_bOrtho[1])
										{
											vec4 = g_pTempTimestep->m_vaCoords[templa[z4+3]];
											vec3 = g_pTempTimestep->m_vaCoords[templa[z4+4]];
											vec5 = g_pTempTimestep->m_vaCoords[templa[z4+5]];
											vec2 = CrossP(FoldVector(vec3-vec4),FoldVector(vec5-vec4));
										} else
										{
											vec4 = g_pTempTimestep->m_vaCoords[templa[z4+3]];
											vec3 = g_pTempTimestep->m_vaCoords[templa[z4+4]];
											vec2 = FoldVector(vec3-vec4);
										}

										tf = Angle_Deg(vec1,vec2);
										if ((tf > 90.0) && (o->m_pADF[zr]->m_bFoldAngle))
											tf = 180.0-tf;
										if (o->m_pADF[zr]->m_bCosine)
											tf = cos(tf/180.0*Pi);

								//		mprintf("z4=%d, templa=%d,%d,%d,%d,%d,%d. vec1=%f|%f|%f, vec2=%f|%f|%f. tf=%f\n",z4,templa[z4],templa[z4+1],templa[z4+2],templa[z4+3],templa[z4+4],templa[z4+5],vec1[0],vec1[1],vec1[2],vec2[0],vec2[1],vec2[2],tf);

										if (g_bDeriv)
										{
											zi = (z2*o->m_iShowMolCount+z3)*o->m_pADF[zr]->m_iCombinations+z4/6;
											o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) = tf;
											switch(o->m_pADF[zr]->m_iDeriv)
											{
												case 0:
													tf = o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(zi);
													break;
												case 1:
													tf = (o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) - o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(zi)) / (2*g_fTimestepLength);
													break;
												case 2:
													tf = (o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(zi) + o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) - 2*o->m_pADF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(zi)) / (g_fTimestepLength * g_fTimestepLength);
													break;
											}
										}
										if (!g_bDeriv || (g_iSteps > 2))
										{
											if (o->m_bObsCertain && o->m_bDecompDist)
												o->m_pADF[zr]->m_pADF->AddToBin_Multi(tic_r*o->m_waObsShowList.GetSize()+tic_o,tf);

											o->m_pADF[zr]->m_faData[z3].Add(tf);

											if (o->m_pADF[zr]->m_bACF)
												o->m_pADF[zr]->m_pfaACFBuffer[(z2*o->m_iShowMolCount+z3)*o->m_pADF[zr]->m_iCombinations+z4/6]->Add((double)tf);

							/*				if (o->m_pADF[zr]->m_bMirror)
											{
												if (o->m_bObsCertain && o->m_bDecompDist)
													o->m_pADF[zr]->m_pADF->AddToBin_Multi(tic_r*o->m_waObsShowList.GetSize()+tic_o,180.0-tf);
												o->m_pADF[zr]->m_faData[z3].Add(180.0-tf);
											}*/

											if (o->m_bTimeDev && (ti2 != -1) && (ti != -1))
											{
										//		o->m_pAngleStat->AddValue(ti2,ti,tf);
												if (o->m_bSaveSeparateFiles)
													mfprintf(o->m_pADF[zr]->m_fAngle[ti2],"; %8.3f",tf);
														else mfprintf(o->m_pADF[zr]->m_fAngle[0],"; %8.3f",tf);
												if (o->m_bCombinedPlot)
													o->m_pADF[zr]->m_pADF->m_pCombinedPlot->AddXYTupel((ti2*o->m_waSaveShowList.GetSize()+ti)*o->m_pADF[zr]->m_iCombinations+z4/6,g_iSteps*g_fTimestepLength/1000.0,tf);
											}
											if (o->m_bTimeDiff)
												((CxDoubleArray*)o->m_pADF[zr]->m_pADF->m_oaTimeDiffBuf[(z2*o->m_iShowMolCount+z3)*o->m_pADF[zr]->m_iCombinations+z4/6])->Add((double)tf);

											if (o->m_bPercTimeDev) {
												o->m_iPercTimeDevTempTotal++;
												if ((tf >= o->m_fPercTimeDevMin) && (tf <= o->m_fPercTimeDevMax))
													o->m_iPercTimeDevTempCounter++;
											}
										}
									}
								}
							}  // Ende IF ADF

							if (g_bDDF)
							{
								for (zr=0;zr<g_iCDFChannels;zr++)
								{
									if (o->m_pDDF[zr] == NULL)
										continue;
									if (o->m_bSecondShowMol)
									{
										if (secondmolrun && (zr == 0))
											continue;
										if (!secondmolrun && (zr == 1))
											continue;
									}
									o->m_pDDF[zr]->BuildAtomList(smfix,sm,&templa);
									for (z4=0;z4<templa.GetSize();z4+=9)
									{
										if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
										{
											tb = true;
											for (z5=z4;z5<z4+9;z5++)
												if (g_baAtomPassedCondition[templa[z5]] == 0)
													tb = false;
											if (tb)
												o->m_pDDF[zr]->m_baDataEnabled[z3].Add(1);
													else o->m_pDDF[zr]->m_baDataEnabled[z3].Add(0);
										}
										if (o->m_pDDF[zr]->m_bOrtho[0])
										{
											vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
											vec2 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
											vec3 = g_pTempTimestep->m_vaCoords[templa[z4+2]];
											vec1 = CrossP(FoldVector(vec2-vec0),FoldVector(vec3-vec0));
										} else
										{
											vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
											vec2 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
											vec1 = FoldVector(vec2-vec0);
										}
										if (o->m_pDDF[zr]->m_bOrtho[1])
										{
											vec4 = g_pTempTimestep->m_vaCoords[templa[z4+3]];
											vec3 = g_pTempTimestep->m_vaCoords[templa[z4+4]];
											vec5 = g_pTempTimestep->m_vaCoords[templa[z4+5]];
											vec2 = CrossP(FoldVector(vec3-vec4),FoldVector(vec5-vec4));
										} else
										{
											vec4 = g_pTempTimestep->m_vaCoords[templa[z4+3]];
											vec3 = g_pTempTimestep->m_vaCoords[templa[z4+4]];
											vec2 = FoldVector(vec3-vec4);
										}
										if (o->m_pDDF[zr]->m_bOrtho[2])
										{
											vec4 = g_pTempTimestep->m_vaCoords[templa[z4+6]];
											vec3 = g_pTempTimestep->m_vaCoords[templa[z4+7]];
											vec5 = g_pTempTimestep->m_vaCoords[templa[z4+8]];
											vec3 = CrossP(FoldVector(vec3-vec4),FoldVector(vec5-vec4));
										} else
										{
											vec4 = g_pTempTimestep->m_vaCoords[templa[z4+6]];
											vec3 = g_pTempTimestep->m_vaCoords[templa[z4+7]];
											vec3 = FoldVector(vec3-vec4);
										}
										tf = Dihedral(vec1,vec2,vec3,o->m_pDDF[zr]->m_bAbs);
										if (o->m_pDDF[zr]->m_bCosine)
											tf = cos(tf/180.0*Pi);
												else if (o->m_pDDF[zr]->m_bPositive)
													if (tf < 0)
														tf += 360.0;
										if (g_bDeriv)
										{
											zi = (z2*o->m_iShowMolCount+z3)*o->m_pDDF[zr]->m_iCombinations+z4/9;
											o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) = tf;
											switch(o->m_pDDF[zr]->m_iDeriv)
											{
												case 0:
													tf = o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(zi);
													break;
												case 1:
													tf = (o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) - o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(zi)) / (2*g_fTimestepLength);
													break;
												case 2:
													tf = (o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(zi) + o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) - 2*o->m_pDDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(zi)) / (g_fTimestepLength * g_fTimestepLength);
													break;
											}
										}
										if (!g_bDeriv || (g_iSteps > 2))
										{
											if (o->m_bObsCertain && o->m_bDecompDist)
												o->m_pDDF[zr]->m_pDDF->AddToBin_Multi(tic_r*o->m_waObsShowList.GetSize()+tic_o,tf);

											o->m_pDDF[zr]->m_faData[z3].Add(tf);

											if (o->m_pDDF[zr]->m_bACF)
												o->m_pDDF[zr]->m_pfaACFBuffer[(z2*o->m_iShowMolCount+z3)*o->m_pDDF[zr]->m_iCombinations+z4/9]->Add((double)tf);

											if (/*(!o->m_pDDF[zr]->m_bAbs) && */o->m_pDDF[zr]->m_bSymm)
											{
												o->m_pDDF[zr]->m_faData[z3].Add(-tf);
												if (o->m_bObsCertain && o->m_bDecompDist)
													o->m_pDDF[zr]->m_pDDF->AddToBin_Multi(tic_r*o->m_waObsShowList.GetSize()+tic_o,-tf);
											}

											if (o->m_bTimeDiff)
												((CxDoubleArray*)o->m_pDDF[zr]->m_pDDF->m_oaTimeDiffBuf[(z2*o->m_iShowMolCount+z3)*o->m_pDDF[zr]->m_iCombinations+z4/9])->Add((double)tf);

											if (o->m_bTimeDev && (ti2 != -1) && (ti != -1))
											{
												if (o->m_pDDF[zr]->m_bRotate)
												{
													if (g_iSteps > 1)
													{
														if ((tf - o->m_pDDF[zr]->m_faLastData[(z2*o->m_iShowMolCount+z3)*o->m_pDDF[zr]->m_iCombinations+z4/9]) > 180.0)
															o->m_pDDF[zr]->m_laRotation[(z2*o->m_iShowMolCount+z3)*o->m_pDDF[zr]->m_iCombinations+z4/9]--;
														if ((tf - o->m_pDDF[zr]->m_faLastData[(z2*o->m_iShowMolCount+z3)*o->m_pDDF[zr]->m_iCombinations+z4/9]) < -180.0)
															o->m_pDDF[zr]->m_laRotation[(z2*o->m_iShowMolCount+z3)*o->m_pDDF[zr]->m_iCombinations+z4/9]++;
													}
													o->m_pDDF[zr]->m_faLastData[(z2*o->m_iShowMolCount+z3)*o->m_pDDF[zr]->m_iCombinations+z4/9] = tf;
													tf2 = tf + o->m_pDDF[zr]->m_laRotation[(z2*o->m_iShowMolCount+z3)*o->m_pDDF[zr]->m_iCombinations+z4/9] * 360.0;
												} else tf2 = tf;
												if (o->m_bSaveSeparateFiles)
													mfprintf(o->m_pDDF[zr]->m_fAngle[ti2],"; %8.3f",tf2);
														else mfprintf(o->m_pDDF[zr]->m_fAngle[0],"; %8.3f",tf2);
												if (o->m_bCombinedPlot)
													o->m_pDDF[zr]->m_pDDF->m_pCombinedPlot->AddXYTupel((ti2*o->m_waSaveShowList.GetSize()+ti)*o->m_pDDF[zr]->m_iCombinations+z4/9,g_iSteps*g_fTimestepLength/1000.0,tf);
											}

											if (o->m_bPercTimeDev) {
												o->m_iPercTimeDevTempTotal++;
												if ((tf >= o->m_fPercTimeDevMin) && (tf <= o->m_fPercTimeDevMax))
													o->m_iPercTimeDevTempCounter++;
											}
										}
									}
								}
							}  // Ende IF DDF

							if (g_bPlDF)
							{
								for (zr=0;zr<g_iCDFChannels;zr++)
								{
									if (o->m_pPlDF[zr] == NULL)
										continue;
									if (o->m_bSecondShowMol)
									{
										if (secondmolrun && (zr == 0))
											continue;
										if (!secondmolrun && (zr == 1))
											continue;
									}
									o->m_pPlDF[zr]->BuildAtomList(smfix,sm,&templa);
									for (z4=0;z4<templa.GetSize();z4+=4)
									{
										if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
										{
											tb = true;
											for (z5=z4;z5<z4+4;z5++)
												if (g_baAtomPassedCondition[templa[z5]] == 0)
													tb = false;
											if (tb)
												o->m_pPlDF[zr]->m_baDataEnabled[z3].Add(1);
													else o->m_pPlDF[zr]->m_baDataEnabled[z3].Add(0);
										}
										if (o->m_pPlDF[zr]->m_bNormal)
										{
											vec0 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+1]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
											vec0.Normalize();
											vec1 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+3]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
											tf = DotP(vec0,vec1);
										} else
										{
											vec0 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+1]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
											vec1 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+2]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
											vec3 = CrossP(vec0,vec1);
											vec3.Normalize();
											vec2 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+3]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
											tf = DotP(vec2,vec3);
										}

										if (!g_bDeriv || (g_iSteps > 2))
										{
											if (o->m_bObsCertain && o->m_bDecompDist)
												o->m_pPlDF[zr]->m_pPlDF->AddToBin_Multi(tic_r*o->m_waObsShowList.GetSize()+tic_o,tf);

											o->m_pPlDF[zr]->m_faData[z3].Add(tf);
										}
									}
								}
							}  // Ende IF PlDF

							if (g_bLiDF)
							{
								for (zr=0;zr<g_iCDFChannels;zr++)
								{
									if (o->m_pLiDF[zr] == NULL)
										continue;
									if (o->m_bSecondShowMol)
									{
										if (secondmolrun && (zr == 0))
											continue;
										if (!secondmolrun && (zr == 1))
											continue;
									}
									o->m_pLiDF[zr]->BuildAtomList(smfix,sm,&templa);
									for (z4=0;z4<templa.GetSize();z4+=4)
									{
										if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
										{
											tb = true;
											for (z5=z4;z5<z4+4;z5++)
												if (g_baAtomPassedCondition[templa[z5]] == 0)
													tb = false;
											if (tb)
												o->m_pLiDF[zr]->m_baDataEnabled[z3].Add(1);
													else o->m_pLiDF[zr]->m_baDataEnabled[z3].Add(0);
										}

										if (o->m_pLiDF[zr]->m_bNormal)
										{
											vec2 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+1]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
											vec3 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+2]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
											vec0 = CrossP(vec2,vec3);
											vec0.Normalize();
											vec1 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+3]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
											vec0 *= DotP(vec0,vec1);
											tf = (vec1 - vec0).GetLength();
										} else
										{
											vec0 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+1]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
											vec0.Normalize();
											vec1 = FoldVector(g_pTempTimestep->m_vaCoords[templa[z4+3]] - g_pTempTimestep->m_vaCoords[templa[z4]]);
											vec0 *= DotP(vec0,vec1);
											tf = (vec1 - vec0).GetLength();
										}

										if (!g_bDeriv || (g_iSteps > 2))
										{
											if (o->m_bObsCertain && o->m_bDecompDist)
												o->m_pLiDF[zr]->m_pLiDF->AddToBin_Multi(tic_r*o->m_waObsShowList.GetSize()+tic_o,tf);

											o->m_pLiDF[zr]->m_faData[z3].Add(tf);
										}
									}
								}
							}  // Ende IF LiDF

							if (g_bDipDF)
							{
								for (zr=0;zr<g_iCDFChannels;zr++)
								{
									if (o->m_pDipDF[zr] == NULL)
										continue;
									if (o->m_bSecondShowMol)
									{
										if (secondmolrun && (zr == 0))
											continue;
										if (!secondmolrun && (zr == 1))
											continue;
									}
									if (o->m_pDipDF[zr]->m_iRefOrSec == 0)
										tf = smfix->m_vDipole.GetLength();
											else tf = sm->m_vDipole.GetLength();
									if (g_bDeriv)
									{
										zi = z2*o->m_iShowMolCount+z3;
										o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) = tf;
										switch(o->m_pDipDF[zr]->m_iDeriv)
										{
											case 0:
												tf = o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(zi);
												break;
											case 1:
												tf = (o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) - o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(zi)) / (2*g_fTimestepLength);
												break;
											case 2:
												tf = (o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(zi) + o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) - 2*o->m_pDipDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(zi)) / (g_fTimestepLength * g_fTimestepLength);
												break;
										}
									}
									if (!g_bDeriv || (g_iSteps > 2))
									{
										if (o->m_bObsCertain && o->m_bDecompDist)
											o->m_pDipDF[zr]->m_pDipoleDF->AddToBin_Multi(tic_r*o->m_waObsShowList.GetSize()+tic_o,tf);

										o->m_pDipDF[zr]->m_faData[z3].Add(tf);

										if (o->m_pDipDF[zr]->m_bACF)
											o->m_pDipDF[zr]->m_pfaACFBuffer[z2*o->m_iShowMolCount+z3]->Add((double)tf);

										if (o->m_bTimeDev && (ti2 != -1) && (ti != -1))
										{
											if (o->m_bSaveSeparateFiles)
												mfprintf(o->m_pDipDF[zr]->m_fDipole[ti2],"; %8.3f",tf);
													else mfprintf(o->m_pDipDF[zr]->m_fDipole[0],"; %8.3f",tf);
											if (o->m_bCombinedPlot)
												o->m_pDipDF[zr]->m_pDipoleDF->m_pCombinedPlot->AddXYTupel(ti2*o->m_waSaveShowList.GetSize()+ti,g_iSteps*g_fTimestepLength/1000.0,tf);
										}
										if (o->m_bTimeDiff)
											((CxDoubleArray*)o->m_pDipDF[zr]->m_pDipoleDF->m_oaTimeDiffBuf[z2*o->m_iShowMolCount+z3])->Add((double)tf);
									}
								}
							} // Ende IF DIPOLE

							if (g_bVDF)
							{
								for (zr=0;zr<g_iCDFChannels;zr++)
								{
									if (o->m_pVDF[zr] == NULL)
										continue;
									if (o->m_bSecondShowMol)
									{
										if (secondmolrun && (zr == 0))
											continue;
										if (!secondmolrun && (zr == 1))
											continue;
									}
									o->m_pVDF[zr]->BuildAtomList(smfix,sm,&templa);
									for (z4=0;z4<templa.GetSize();z4++)
									{
										if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
										{
											tb = true;
											for (z5=z4;z5<z4+2;z5++)
												if (g_baAtomPassedCondition[templa[z5]] == 0)
													tb = false;
											if (tb)
												o->m_pVDF[zr]->m_baDataEnabled[z3].Add(1);
													else o->m_pVDF[zr]->m_baDataEnabled[z3].Add(0);
										}
										tf = tempvel[templa[z4]].GetLength();
										if (g_bDeriv)
										{
											zi = (z2*o->m_iShowMolCount+z3)*o->m_pVDF[zr]->m_iCombinations+z4/2;
											o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) = tf;
											switch(o->m_pVDF[zr]->m_iDeriv)
											{
												case 0:
													tf = o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(zi);
													break;
												case 1:
													tf = (o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) - o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(zi)) / (2*g_fTimestepLength);
													break;
												case 2:
													tf = (o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(zi) + o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) - 2*o->m_pVDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(zi)) / (g_fTimestepLength * g_fTimestepLength);
													break;
											}
										}
										if (!g_bDeriv || (g_iSteps > 2))
										{
											if (o->m_bObsCertain && o->m_bDecompDist)
												o->m_pVDF[zr]->m_pVDF->AddToBin_Multi(tic_r*o->m_waObsShowList.GetSize()+tic_o,tf);

											o->m_pVDF[zr]->m_faData[z3].Add(tf);

											if (o->m_pVDF[zr]->m_bACF)
												o->m_pVDF[zr]->m_pfaACFBuffer[(z2*o->m_iShowMolCount+z3)*o->m_pVDF[zr]->m_iCombinations+z4/2]->Add((double)tf);

											if (o->m_bTimeDev && (ti2 != -1) && (ti != -1))
											{
												if (o->m_bSaveSeparateFiles)
													mfprintf(o->m_pVDF[zr]->m_fSpeed[ti2],"; %8.3f",tf);
														else mfprintf(o->m_pVDF[zr]->m_fSpeed[0],"; %8.3f",tf);
												if (o->m_bCombinedPlot)
													o->m_pVDF[zr]->m_pVDF->m_pCombinedPlot->AddXYTupel((ti2*o->m_waSaveShowList.GetSize()+ti)*o->m_pVDF[zr]->m_iCombinations+z4/2,g_iSteps*g_fTimestepLength/1000.0,tf);
											}
											if (o->m_bTimeDiff)
												((CxDoubleArray*)o->m_pVDF[zr]->m_pVDF->m_oaTimeDiffBuf[(z2*o->m_iShowMolCount+z3)*o->m_pVDF[zr]->m_iCombinations+z4/2])->Add((double)tf);
										}
									}
								}
							} // Ende IF VDF

							if (g_bRDF)
							{
								for (zr=0;zr<g_iCDFChannels;zr++)
								{
									if (o->m_pRDF[zr] == NULL)
										continue;
									if (o->m_bSecondShowMol)
									{
										if (secondmolrun && (zr == 0))
											continue;
										if (!secondmolrun && (zr == 1))
											continue;
									}
									o->m_pRDF[zr]->BuildAtomList(smfix,sm,&templa);
									for (z4=0;z4<templa.GetSize();z4+=2)
									{
										if (o->m_bBinOnlyPassedAtoms || o->m_bBinOnlyNotPassedAtoms)
										{
											tb = true;
											for (z5=z4;z5<z4+2;z5++)
												if (g_baAtomPassedCondition[templa[z5]] == 0)
													tb = false;
											if (tb)
												o->m_pRDF[zr]->m_baDataEnabled[z3].Add(1);
													else o->m_pRDF[zr]->m_baDataEnabled[z3].Add(0);
										}
										vec0 = g_pTempTimestep->m_vaCoords[templa[z4]];
										vec1 = g_pTempTimestep->m_vaCoords[templa[z4+1]];
										tf = FoldedLength(vec0-vec1);
								/*		if (tf > 2000)
										{
											mprintf("\n@ z2=%d, z3=%d, z4=%d, templa[z4]=%d, templa[z4+1]=%d, dist=%f  ",z2,z3,z4,templa[z4],templa[z4+1],tf);
											mprintf("\n  vec0: "); vec0.Dump(); mprintf(", vec1: "); vec1.Dump(); mprintf(", vec2: "); vec2.Dump();
										}
						*/		//		mprintf("\n@ z2=%d, z3=%d, z4=%d, templa[z4]=%d, templa[z4+1]=%d, dist=%f  ",z2,z3,z4,templa[z4],templa[z4+1],tf);
								//		mprintf("\n  vec0: "); vec0.Dump(); mprintf(", vec1: "); vec1.Dump(); mprintf(", vec2: "); vec2.Dump();
										if (g_bDeriv)
										{
											zi = (z2*o->m_iShowMolCount+z3)*o->m_pRDF[zr]->m_iCombinations+z4/2;
											o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) = tf;
											switch(o->m_pRDF[zr]->m_iDeriv)
											{
												case 0:
													tf = o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(zi);
													break;
												case 1:
													tf = (o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) - o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(zi)) / (2*g_fTimestepLength);
													break;
												case 2:
													tf = (o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivLast]->GetAt(zi) + o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivNext]->GetAt(zi) - 2*o->m_pRDF[zr]->m_pfaDerivBuffer[g_iDerivCurr]->GetAt(zi)) / (g_fTimestepLength * g_fTimestepLength);
													break;
											}
										}
										if (!g_bDeriv || (g_iSteps > 2))
										{
											if (o->m_bObsCertain && o->m_bDecompDist)
												o->m_pRDF[zr]->m_pRDF->AddToBin_Multi(tic_r*o->m_waObsShowList.GetSize()+tic_o,tf);

											if (o->m_bDecompType)
												o->m_pRDF[zr]->m_pRDF->AddToBin_Multi(o->m_waDecompTypeRefOffs[templa[z4]]*o->m_waDecompTypeObsOffs.GetSize()+o->m_waDecompTypeObsOffs[templa[z4+1]],tf);

											o->m_pRDF[zr]->m_faData[z3].Add(tf);

											if (o->m_pRDF[zr]->m_bACF)
												o->m_pRDF[zr]->m_pfaACFBuffer[(z2*o->m_iShowMolCount+z3)*o->m_pRDF[zr]->m_iCombinations+z4/2]->Add((double)tf);

											if (o->m_bTimeDev && (ti2 != -1) && (ti != -1))
											{
												if (o->m_bSaveSeparateFiles)
													mfprintf(o->m_pRDF[zr]->m_fDist[ti2],"; %10.3f",tf);
														else mfprintf(o->m_pRDF[zr]->m_fDist[0],"; %10.3f",tf);
												if (o->m_bCombinedPlot)
													o->m_pRDF[zr]->m_pRDF->m_pCombinedPlot->AddXYTupel((ti2*o->m_waSaveShowList.GetSize()+ti)*o->m_pRDF[zr]->m_iCombinations+z4/2,g_iSteps*g_fTimestepLength/1000.0,tf);
											}
											if (o->m_bTimeDiff)
												((CxDoubleArray*)o->m_pRDF[zr]->m_pRDF->m_oaTimeDiffBuf[(z2*o->m_iShowMolCount+z3)*o->m_pRDF[zr]->m_iCombinations+z4/2])->Add((double)tf);

											if (o->m_bPercTimeDev) {
												o->m_iPercTimeDevTempTotal++;
												if ((tf >= o->m_fPercTimeDevMin) && (tf <= o->m_fPercTimeDevMax))
													o->m_iPercTimeDevTempCounter++;
											}
										}
									} 
								}
							} // Ende IF RDF

							if (g_bSDF)
							{
								o->m_pSDF->BuildAtomList(smfix,sm,&templa);
								for (z4=0;z4<templa.GetSize();z4++)
								{
									{
										o->m_pSDF->m_vaData[z3].Add(g_pTempTimestep->m_vaCoords[templa[z4]]);

										if (o->m_pSDF->m_bVdWSpheres)
											o->m_pSDF->m_faRadius[z3].Add(g_faVdWRadius[templa[z4]]);
									}
								}
							} // Ende IF SDF

							if (g_bPlProj)
							{
								o->m_pPlProj->BuildAtomList(smfix,sm,&templa);
								for (z4=0;z4<templa.GetSize();z4++)
									o->m_pPlProj->m_vaData[z3].Add(g_pTempTimestep->m_vaCoords[templa[z4]]);

								if (o->m_pPlProj->m_bVector) {
									o->m_pPlProj->BuildVectorAtomList(smfix,sm,&templa);
									for (z4=0;z4<templa.GetSize()/3;z4++) {
										o->m_pPlProj->m_vaVectorData[z3].Add(
											g_pTempTimestep->m_vaCoords[templa[z4*3+1]] - g_pTempTimestep->m_vaCoords[templa[z4*3]] );
										o->m_pPlProj->m_vaVectorData[z3].Add(
											g_pTempTimestep->m_vaCoords[templa[z4*3+2]] );
									}
								}
							} // Ende IF PlProj

							if (g_bRevSDF)
							{
								o->m_pRevSDF->BuildAtomList(smfix,sm,&templa);
								for (z4=0;z4<templa.GetSize();z4++)
									o->m_pRevSDF->m_vaData[z3].Add(g_pTempTimestep->m_vaCoords[templa[z4]]);
							} // Ende IF RevSDF

						} // Ende FOR Alle anderen Molekuele durchgehen
						if (o->m_bSecondShowMol && (!secondmolrun))
						{
							secondmolrun = true;
							goto _secondmolstart;
						}
					} // Ende IF m_bOthers

					if (g_bSDF && !g_bSDFVoro)
					{
						if (o->m_pSDF->m_bSDFMirrorBisect)
						{
							vec2 = g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[1]])->GetAt(g_iFixAtom[1])];
							vec3 = g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[2]])->GetAt(g_iFixAtom[2])];
							vec1 = CrossP(vec2,vec3);
							vec2 += vec3;
						}

						if (o->m_pSDF->m_bCutPlane)
						{
							o->m_pSDF->m_fAtom2PosX += g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[1]])->GetAt(g_iFixAtom[1])][0];
							o->m_pSDF->m_fAtom3PosX += g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[2]])->GetAt(g_iFixAtom[2])][0];
							o->m_pSDF->m_fAtom3PosY += g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[2]])->GetAt(g_iFixAtom[2])][1];
							o->m_pSDF->m_fPosCounter++;
						}

						if (o->m_bSelf)
						{
							for (z5=0;z5<o->m_pSDF->m_vaData[0].GetSize();z5++)
							{
								vec0 = o->m_pSDF->m_vaData[0][z5];

								if (o->m_pSDF->m_fMinRadiusSqr != 0)
									if (vec0.GetLengthSqr() < o->m_pSDF->m_fMinRadiusSqr)
										continue;

								if (o->m_pSDF->m_bRadiusCut)
									if (vec0.GetLengthSqr() > o->m_pSDF->m_fRadiusSqr)
										continue;

								if (o->m_pSDF->m_bVdWSpheres)
									o->m_pSDF->m_pSDF->AddToBin_Sphere(vec0,o->m_pSDF->m_faRadius[0][z5]);
								else
									o->m_pSDF->m_pSDF->AddToBin(vec0);

								if (o->m_pSDF->m_bSDFMirrorXY)
								{
									vec0[2] = -vec0[2];
									o->m_pSDF->m_pSDF->AddToBin(vec0);
								}

								if (o->m_pSDF->m_bSDFMirrorBisect)
								{
									vec3.PointRoot(vec1,vec2,vec0);
									vec3 -= vec0;
									vec3 *= 2.0;
									vec0 += vec3;
									o->m_pSDF->m_pSDF->AddToBin(vec0);
									if (o->m_pSDF->m_bSDFMirrorXY)
									{
										vec0[2] = -vec0[2];
										o->m_pSDF->m_pSDF->AddToBin(vec0);
									}
								}
							}
						} else
						{
							for (z4=0;z4<o->m_iShowMolCount;z4++)
							{
								for (z5=0;z5<o->m_pSDF->m_vaData[z4].GetSize();z5++)
								{
									vec0 = o->m_pSDF->m_vaData[z4][z5];

									if (o->m_pSDF->m_fMinRadiusSqr != 0)
										if (vec0.GetLengthSqr() < o->m_pSDF->m_fMinRadiusSqr)
											continue;

									if (o->m_pSDF->m_bRadiusCut)
										if (vec0.GetLengthSqr() > o->m_pSDF->m_fRadiusSqr)
											continue;

									if (o->m_pSDF->m_bVdWSpheres)
										o->m_pSDF->m_pSDF->AddToBin_Sphere(vec0,o->m_pSDF->m_faRadius[z4][z5]);
									else
										o->m_pSDF->m_pSDF->AddToBin(vec0);

									if (o->m_pSDF->m_bSDFMirrorXY)
									{
										vec0[2] = -vec0[2];
										o->m_pSDF->m_pSDF->AddToBin(vec0);
									}

									if (o->m_pSDF->m_bSDFMirrorBisect)
									{
										vec3.PointRoot(vec1,vec2,vec0);
										vec3 -= vec0;
										vec3 *= 2.0;
										vec0 += vec3;
										o->m_pSDF->m_pSDF->AddToBin(vec0);
										if (o->m_pSDF->m_bSDFMirrorXY)
										{
											vec0[2] = -vec0[2];
											o->m_pSDF->m_pSDF->AddToBin(vec0);
										}
									}
								}
							}
						}
					} // END IF SDF

					if (g_bPlProj)
					{
						if (o->m_bSelf)
						{
							for (z5=0;z5<o->m_pPlProj->m_vaData[0].GetSize();z5++)
							{
								vec0 = o->m_pPlProj->m_vaData[0][z5];

								if ((vec0[2] >= o->m_pPlProj->m_fSliceBorder[0]) && (vec0[2] <= o->m_pPlProj->m_fSliceBorder[1]))
									o->m_pPlProj->m_p2DF->AddToBin(vec0[0],vec0[1]);
							}
						} else
						{
							for (z4=0;z4<o->m_iShowMolCount;z4++)
							{
								for (z5=0;z5<o->m_pPlProj->m_vaData[z4].GetSize();z5++)
								{
									vec0 = o->m_pPlProj->m_vaData[z4][z5];

									if ((vec0[2] >= o->m_pPlProj->m_fSliceBorder[0]) && (vec0[2] <= o->m_pPlProj->m_fSliceBorder[1]))
										o->m_pPlProj->m_p2DF->AddToBin(vec0[0],vec0[1]);
								}

								if (o->m_pPlProj->m_bVector) {
									for (z5=0;z5<o->m_pPlProj->m_vaVectorData[z4].GetSize()/2;z5++) {

										vec0 = o->m_pPlProj->m_vaVectorData[z4][z5*2];
										vec1 = o->m_pPlProj->m_vaVectorData[z4][z5*2+1];

										if ((vec1[0] >= o->m_pPlProj->m_p2DF->m_fMinVal[0]) && (vec1[0] <= o->m_pPlProj->m_p2DF->m_fMaxVal[0]) &&
											(vec1[1] >= o->m_pPlProj->m_p2DF->m_fMinVal[1]) && (vec1[1] <= o->m_pPlProj->m_p2DF->m_fMaxVal[1]) &&
											(vec1[2] >= o->m_pPlProj->m_fSliceBorder[0]) && (vec1[2] <= o->m_pPlProj->m_fSliceBorder[1])) {

											ix = (int)floor((vec1[0]-o->m_pPlProj->m_p2DF->m_fMinVal[0])/(o->m_pPlProj->m_p2DF->m_fMaxVal[0]-o->m_pPlProj->m_p2DF->m_fMinVal[0])*o->m_pPlProj->m_iVecRes[0]);
											iy = (int)floor((vec1[1]-o->m_pPlProj->m_p2DF->m_fMinVal[1])/(o->m_pPlProj->m_p2DF->m_fMaxVal[1]-o->m_pPlProj->m_p2DF->m_fMinVal[1])*o->m_pPlProj->m_iVecRes[1]);

											vec2 = vec0;
											vec2[2] = 0;
											vec2.Normalize();

											o->m_pPlProj->m_faVectorData[(iy*o->m_pPlProj->m_iVecRes[0]+ix)*2]   += vec2[0];
											o->m_pPlProj->m_faVectorData[(iy*o->m_pPlProj->m_iVecRes[0]+ix)*2+1] += vec2[1];
											o->m_pPlProj->m_faVectorCount[iy*o->m_pPlProj->m_iVecRes[0]+ix]      += 1.0;
										}
									}
								}
							}
						}
					} // END IF PlProj

					if (g_bRevSDF)
					{
						vec2 = g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[1]])->GetAt(g_iFixAtom[1])];
						tf4 = g_pTempTimestep->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[g_iFixAtomType[1]])->GetAt(g_iFixAtom[1])].GetLength();
						o->m_pRevSDF->m_fSecondAtomPosX += tf4;
						o->m_pRevSDF->m_fSecondAtomCount++;
						if (o->m_bSelf)
						{
							for (z5=0;z5<o->m_pRevSDF->m_vaData[0].GetSize();z5++)
							{
								vec0 = o->m_pRevSDF->m_vaData[0][z5];
								tf = fabs(Angle(vec2,vec0));
								tf2 = vec0.GetLength();
								tf3 = 1.0;
								if (o->m_pRevSDF->m_bCorrectAngle)
								{
									if ((tf < 0.001) || (tf > Pi-0.001))
										tf3 /= 0.001;
											else tf3 /= sin(tf);
								}
								if (o->m_pRevSDF->m_bCorrectRadial)
									tf3 /= tf2; // Erste Potenz, da 2D-Grid (SDF hat "nullte Potenz" ^^)
								o->m_pRevSDF->m_p2DF->AddToBin(sin(tf)*tf2,cos(tf)*tf2,tf3);
								o->m_pRevSDF->m_p2DF->AddToBin(-sin(tf)*tf2,cos(tf)*tf2,tf3);
								if (o->m_pRevSDF->m_bMirrorY)
								{
									if (o->m_pRevSDF->m_bMirrorBond)
									{
										o->m_pRevSDF->m_p2DF->AddToBin(sin(tf)*tf2,tf4-cos(tf)*tf2,tf3);
										o->m_pRevSDF->m_p2DF->AddToBin(-sin(tf)*tf2,tf4-cos(tf)*tf2,tf3);
									} else
									{
										o->m_pRevSDF->m_p2DF->AddToBin(sin(tf)*tf2,-cos(tf)*tf2,tf3);
										o->m_pRevSDF->m_p2DF->AddToBin(-sin(tf)*tf2,-cos(tf)*tf2,tf3);
									}
								}
							}
						} else
						{
							for (z4=0;z4<o->m_iShowMolCount;z4++)
							{
								for (z5=0;z5<o->m_pRevSDF->m_vaData[z4].GetSize();z5++)
								{
									vec0 = o->m_pRevSDF->m_vaData[z4][z5];
									tf = fabs(Angle(vec2,vec0));
									tf2 = vec0.GetLength();
									tf3 = 1.0;
									if (o->m_pRevSDF->m_bCorrectAngle)
									{
										if ((tf < 0.001) || (tf > Pi-0.001))
											tf3 /= 0.001;
												else tf3 /= sin(tf);
									}
									if (o->m_pRevSDF->m_bCorrectRadial)
										tf3 /= tf2; // Erste Potenz, da 2D-Grid (SDF hat "nullte Potenz" ^^)
									o->m_pRevSDF->m_p2DF->AddToBin(sin(tf)*tf2,cos(tf)*tf2,tf3);
									o->m_pRevSDF->m_p2DF->AddToBin(-sin(tf)*tf2,cos(tf)*tf2,tf3);
									if (o->m_pRevSDF->m_bMirrorY)
									{
										if (o->m_pRevSDF->m_bMirrorBond)
										{
											o->m_pRevSDF->m_p2DF->AddToBin(sin(tf)*tf2,tf4-cos(tf)*tf2,tf3);
											o->m_pRevSDF->m_p2DF->AddToBin(-sin(tf)*tf2,tf4-cos(tf)*tf2,tf3);
										} else
										{
											o->m_pRevSDF->m_p2DF->AddToBin(sin(tf)*tf2,-cos(tf)*tf2,tf3);
											o->m_pRevSDF->m_p2DF->AddToBin(-sin(tf)*tf2,-cos(tf)*tf2,tf3);
										}
									}
								}
							}
						}
					} // END IF REVSDF

					if (g_bCDF && (!g_bDeriv || (g_iSteps > 2))) {

						if (o->m_pCDF->m_bWriteSnapshots) {

							if (o->m_pCDF->m_bWriteSnapshotsOnlyRMOM && !o->m_pCDF->m_bWriteSnapshotsNoWrite) {

								if (g_pTempTimestepSnap == NULL)
									g_pTempTimestepSnap = new CTimeStep();
								g_pTempTimestepSnap->CopyFrom( g_pTempTimestep );

								if (o->m_pCDF->m_bWriteSnapshotsCenterRM) {

									vecsnap1 = g_pTempTimestepSnap->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[ o->m_pCDF->m_iWriteSnapshotFixType[0] ])->GetAt( o->m_pCDF->m_iWriteSnapshotFixAtom[0] )];
									g_pTempTimestepSnap->CenterPos(vecsnap1);

									g_pTempTimestepSnap->FoldMolecules();

									if (o->m_pCDF->m_bWriteSnapshotsFixRM) {

										vecsnap2 = g_pTempTimestepSnap->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[ o->m_pCDF->m_iWriteSnapshotFixType[1] ])->GetAt( o->m_pCDF->m_iWriteSnapshotFixAtom[1] )];
										vecsnap3 = g_pTempTimestepSnap->m_vaCoords[((CxIntArray*)smfix->m_oaAtomOffset[ o->m_pCDF->m_iWriteSnapshotFixType[2] ])->GetAt( o->m_pCDF->m_iWriteSnapshotFixAtom[2] )];

										matsnap.MatUltra( vecsnap2, vecsnap3 );

										g_pTempTimestepSnap->Transform( matsnap );
									}
								}
							}
						}

						if (o->m_bSecondShowMol) {

							for (zr=0;zr<g_iCDFChannels;zr++) {

								switch(g_iObsChannel[zr]) {
									case 0: 
										apfa[zr] = o->m_pRDF[zr]->m_faData; 
										apba[zr] = o->m_pRDF[zr]->m_baDataEnabled; 
										break;
									case 1: 
										apfa[zr] = o->m_pADF[zr]->m_faData; 
										apba[zr] = o->m_pADF[zr]->m_baDataEnabled; 
										break;
									case 2: 
										apfa[zr] = o->m_pDDF[zr]->m_faData; 
										apba[zr] = o->m_pDDF[zr]->m_baDataEnabled; 
										break;
									case 3: 
										apfa[zr] = o->m_pDipDF[zr]->m_faData; 
										apba[zr] = o->m_pDipDF[zr]->m_baDataEnabled; 
										break;
									case 4: 
										apfa[zr] = o->m_pVDF[zr]->m_faData; 
										apba[zr] = o->m_pVDF[zr]->m_baDataEnabled; 
										break;
									case 5: 
										apfa[zr] = o->m_pPlDF[zr]->m_faData; 
										apba[zr] = o->m_pPlDF[zr]->m_baDataEnabled; 
										break;
									case 6: 
										apfa[zr] = o->m_pLiDF[zr]->m_faData; 
										apba[zr] = o->m_pLiDF[zr]->m_baDataEnabled; 
										break;
								}
							}

							if (g_iCDFChannels == 2) {

								for (z4=0;z4<o->m_iShowMolCount;z4++) {

									for (z5=0;z5<apfa[0][z4].GetSize();z5++) {

										for (z8=0;z8<o->m_iShowMol2Count;z8++) {

											if (o->m_bExclude1eq2 && (z4 == z8))
												continue;

											for (z7=0;z7<apfa[1][z8].GetSize();z7++) {

												if (o->m_bBinOnlyPassedAtoms)
													if ((apba[0][z4][z5] == 0) || (apba[1][z8][z7] == 0))
														goto _nocdfbin2;

												if (o->m_bBinOnlyNotPassedAtoms)
													if ((apba[0][z4][z5] != 0) && (apba[1][z8][z7] != 0))
														goto _nocdfbin2;

												o->m_pCDF->m_p2DF->AddToBin( apfa[0][z4][z5], apfa[1][z8][z7] );

												if (o->m_pCDF->m_bWriteSnapshots) {

													if ((apfa[0][z4][z5] >= o->m_pCDF->m_faWriteSnapshotsRange[0]) &&
														(apfa[0][z4][z5] <= o->m_pCDF->m_faWriteSnapshotsRange[1]) &&
														(apfa[1][z8][z7] >= o->m_pCDF->m_faWriteSnapshotsRange[2]) &&
														(apfa[1][z8][z7] <= o->m_pCDF->m_faWriteSnapshotsRange[3])) {

														if (!o->m_pCDF->m_bWriteSnapshotsTempRM) {
															o->m_pCDF->m_bWriteSnapshotsTempRM = true;
															o->m_pCDF->m_iWriteSnapshotsTotalPerRM++;
														}

														o->m_pCDF->m_iWriteSnapshotsTotalCounter++;

														if (o->m_pCDF->m_bWriteSnapshotsNoWrite)
															goto _nocdfwritesnapshot;

														if (o->m_pCDF->m_iWriteSnapshotStride > 1) {
															o->m_pCDF->m_iWriteSnapshotCounter++;
															if (o->m_pCDF->m_iWriteSnapshotCounter % o->m_pCDF->m_iWriteSnapshotStride != 0)
																goto _nocdfwritesnapshot;
														}

														bufcom.sprintf(
															"Step %lu, RM=%d, OM1=%d, OM2=%d, Val1=%f, Val2=%f",
															g_iSteps,
															smfix->m_iMolSMIndex+1,
															z4+1,
															z8+1,
															apfa[0][z4][z5],
															apfa[1][z8][z7]
														);

														if (o->m_pCDF->m_bWriteSnapshotsOnlyRMOM) {

															tiasnap.clear();
															tiasnap.push_back( smfix->m_iIndex );
															tiasnap.push_back( ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex[z4] );
															tiasnap.push_back( ((CMolecule*)g_oaMolecules[o->m_iShowMol2])->m_laSingleMolIndex[z8] );

															g_pTempTimestepSnap->ExportXYZConfiguration(
																o->m_pCDF->m_pWriteSnapshotsFile,
																(const char*)bufcom,
																o->m_pCDF->m_bWriteSnapshotsOriginalCoords,
																&tiasnap
															);

															o->m_pCDF->m_iWriteSnapshotsTotalWritten++;

														} else {

															if (!o->m_pCDF->m_bWriteSnapshotDone) {

																// Only write a snapshot once per trajectory frame, as all molecules are included
																o->m_pCDF->m_bWriteSnapshotDone = true;

																g_pTempTimestep->ExportXYZConfiguration(
																	o->m_pCDF->m_pWriteSnapshotsFile,
																	(const char*)bufcom,
																	o->m_pCDF->m_bWriteSnapshotsOriginalCoords,
																	NULL
																);

																o->m_pCDF->m_iWriteSnapshotsTotalWritten++;
															}
														}
													}
												}
_nocdfwritesnapshot:
_nocdfbin2:
												if (o->m_pCDF->m_bDumpDat)
													mfprintf(o->m_pCDF->m_fDump,"%d;  %d;  %d;  %d;  %10.3f;  %10.3f\n",(int)g_iSteps,z2+1,z4+1,z5+1,apfa[0][z4][z5],apfa[1][z8][z7]);
											}
										}
									}
								}
							}
						} else // IF NOT SECOND_OM
						{
							for (zr=0;zr<g_iCDFChannels;zr++)
							{
								switch(g_iObsChannel[zr])
								{
									case 0: 
										apfa[zr] = o->m_pRDF[zr]->m_faData; 
										apba[zr] = o->m_pRDF[zr]->m_baDataEnabled; 
										break;
									case 1: 
										apfa[zr] = o->m_pADF[zr]->m_faData; 
										apba[zr] = o->m_pADF[zr]->m_baDataEnabled; 
										break;
									case 2: 
										apfa[zr] = o->m_pDDF[zr]->m_faData; 
										apba[zr] = o->m_pDDF[zr]->m_baDataEnabled; 
										break;
									case 3: 
										apfa[zr] = o->m_pDipDF[zr]->m_faData; 
										apba[zr] = o->m_pDipDF[zr]->m_baDataEnabled; 
										break;
									case 4: 
										apfa[zr] = o->m_pVDF[zr]->m_faData; 
										apba[zr] = o->m_pVDF[zr]->m_baDataEnabled; 
										break;
									case 5: 
										apfa[zr] = o->m_pPlDF[zr]->m_faData; 
										apba[zr] = o->m_pPlDF[zr]->m_baDataEnabled; 
										break;
									case 6: 
										apfa[zr] = o->m_pLiDF[zr]->m_faData; 
										apba[zr] = o->m_pLiDF[zr]->m_baDataEnabled; 
										break;
								}
							}

							if (g_iCDFChannels == 2)
							{
								for (z4=0;z4<o->m_iShowMolCount;z4++)
								{
									ticomb = 0;
									for (z5=0;z5<apfa[0][z4].GetSize();z5++)
									{
										for (z7=0;z7<apfa[1][z4].GetSize();z7++)
										{

											/* Workaround for symmetrical DDFs (double data amount) */
											if (g_iObsChannel[0] == 2)
											{
												if (o->m_pDDF[0]->m_bSymm)
												{
													if (g_iObsChannel[1] == 2)
													{
														if (o->m_pDDF[1]->m_bSymm)
														{
															if (o->m_pCDF->m_pCombineList[(z5/2)*o->m_pCDF->m_iCombinations[1]+(z7/2)] == 0)
																continue;
															goto _ddfsymmdone;
														}
													}
												}
												if (o->m_pCDF->m_pCombineList[(z5/2)*o->m_pCDF->m_iCombinations[1]+z7] == 0)
													continue;
											} else
											{
												if (g_iObsChannel[1] == 2)
												{
													if (o->m_pDDF[1]->m_bSymm)
													{
														if (o->m_pCDF->m_pCombineList[z5*o->m_pCDF->m_iCombinations[1]+(z7/2)] == 0)
															continue;
														goto _ddfsymmdone;
													}
												}
												if (o->m_pCDF->m_pCombineList[z5*o->m_pCDF->m_iCombinations[1]+z7] == 0)
													continue;
											}
_ddfsymmdone:

											if (o->m_bBinOnlyPassedAtoms)
											{
												if ((apba[0][z4][z5] == 0) || (apba[1][z4][z7] == 0))
													goto _nocdfbin;
											}
											if (o->m_bBinOnlyNotPassedAtoms)
											{
												if ((apba[0][z4][z5] != 0) && (apba[1][z4][z7] != 0))
													goto _nocdfbin;
											}

											o->m_pCDF->m_p2DF->AddToBin(apfa[0][z4][z5],apfa[1][z4][z7]);

											if (o->m_pCDF->m_bWriteSnapshots) {

												if ((apfa[0][z4][z5] >= o->m_pCDF->m_faWriteSnapshotsRange[0]) &&
													(apfa[0][z4][z5] <= o->m_pCDF->m_faWriteSnapshotsRange[1]) &&
													(apfa[1][z4][z7] >= o->m_pCDF->m_faWriteSnapshotsRange[2]) &&
													(apfa[1][z4][z7] <= o->m_pCDF->m_faWriteSnapshotsRange[3])) {

													if (!o->m_pCDF->m_bWriteSnapshotsTempRM) {
														o->m_pCDF->m_bWriteSnapshotsTempRM = true;
														o->m_pCDF->m_iWriteSnapshotsTotalPerRM++;
													}

													o->m_pCDF->m_iWriteSnapshotsTotalCounter++;

													if (o->m_pCDF->m_bWriteSnapshotsNoWrite)
														goto _nocdfwritesnapshot2;

													if (o->m_pCDF->m_iWriteSnapshotStride > 1) {
														o->m_pCDF->m_iWriteSnapshotCounter++;
														if (o->m_pCDF->m_iWriteSnapshotCounter % o->m_pCDF->m_iWriteSnapshotStride != 0)
															goto _nocdfwritesnapshot2;
													}

													bufcom.sprintf(
														"Step %lu, RM=%d, OM=%d, Val1=%f, Val2=%f",
														g_iSteps,
														smfix->m_iMolSMIndex+1,
														z4+1,
														apfa[0][z4][z5],
														apfa[1][z4][z7]
													);

													if (o->m_pCDF->m_bWriteSnapshotsOnlyRMOM) {

														tiasnap.clear();
														tiasnap.push_back( smfix->m_iIndex );
														tiasnap.push_back( ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex[z4] );

														g_pTempTimestepSnap->ExportXYZConfiguration(
															o->m_pCDF->m_pWriteSnapshotsFile,
															(const char*)bufcom,
															o->m_pCDF->m_bWriteSnapshotsOriginalCoords,
															&tiasnap
														);

														o->m_pCDF->m_iWriteSnapshotsTotalWritten++;

													} else {

														if (!o->m_pCDF->m_bWriteSnapshotDone) {

															// Only write a snapshot once per trajectory frame, as all molecules are included
															o->m_pCDF->m_bWriteSnapshotDone = true;

															g_pTempTimestep->ExportXYZConfiguration(
																o->m_pCDF->m_pWriteSnapshotsFile,
																(const char*)bufcom,
																o->m_pCDF->m_bWriteSnapshotsOriginalCoords,
																NULL
															);

															o->m_pCDF->m_iWriteSnapshotsTotalWritten++;
														}
													}
												}
											}
_nocdfwritesnapshot2:
_nocdfbin:
											if (o->m_pCDF->m_bDumpDat)
												mfprintf(o->m_pCDF->m_fDump,"%d;  %d;  %d;  ;  %10.3f;  %10.3f\n",(int)g_iSteps,z2+1,z4+1,apfa[0][z4][z5],apfa[1][z4][z7]);

											if (o->m_bTimeDev)
											{
												ti = o->m_waSaveShowList.GetPosition((unsigned short)z4);
												if ((ti2 != -1) && (ti != -1))
												{
													if (o->m_bSaveSeparateFiles)
														mfprintf(o->m_pCDF->m_fTimeDev[ti2],"; %10.3f; %10.3f",apfa[0][z4][z5],apfa[1][z4][z7]);
															else mfprintf(o->m_pCDF->m_fTimeDev[0],"; %10.3f; %10.3f",apfa[0][z4][z5],apfa[1][z4][z7]);
													if (o->m_pCDF->m_bTDAnimation)
														o->m_pCDF->m_pTDAPlot->AddXYTupel((ti2*o->m_waSaveShowList.GetSize()+ti)*o->m_pCDF->m_iCombinationsEnabled+ticomb,apfa[0][z4][z5],apfa[1][z4][z7]);
												}
												ticomb++;
											}
										}
									}
								}
								if (o->m_bTimeDev && (ti2 != -1))
								{
									if (o->m_bSaveSeparateFiles)
										mfprintf(o->m_pCDF->m_fTimeDev[ti2],"\n");
											else if (ti2 == (int)o->m_waSaveRefList.GetSize()-1)
												mfprintf(o->m_pCDF->m_fTimeDev[0],"\n");
								}
							} // END IF CHANNELS == 2

							if (g_iCDFChannels == 3)
							{
								for (z4=0;z4<o->m_iShowMolCount;z4++)
								{
									ticomb = 0;
									for (z5=0;z5<apfa[0][z4].GetSize();z5++)
									{
										for (z8=0;z8<apfa[1][z4].GetSize();z8++)
										{
											for (z7=0;z7<apfa[2][z4].GetSize();z7++)
											{
												if (o->m_pCDF->m_pCombineList[z5*o->m_pCDF->m_iCombinations[1]*o->m_pCDF->m_iCombinations[2]+z8*o->m_pCDF->m_iCombinations[1]+z7] == 0)
													continue;
												if (o->m_bBinOnlyPassedAtoms)
												{
													if ((apba[0][z4][z5] == 0) || (apba[1][z4][z8] == 0) || (apba[2][z4][z7] == 0))
														goto _3nocdfbin;
												}
												if (o->m_bBinOnlyNotPassedAtoms)
												{
													if ((apba[0][z4][z5] != 0) && (apba[1][z4][z8] != 0) && (apba[2][z4][z7] != 0))
														goto _3nocdfbin;
												}

												o->m_pCDF->m_p3DF->AddToBin(apfa[0][z4][z5],apfa[1][z4][z8],apfa[2][z4][z7]);

												o->m_pCDF->m_p3DF->m_p2DF[0]->AddToBin(apfa[0][z4][z5],apfa[1][z4][z8]);
												o->m_pCDF->m_p3DF->m_p2DF[1]->AddToBin(apfa[0][z4][z5],apfa[2][z4][z7]);
												o->m_pCDF->m_p3DF->m_p2DF[2]->AddToBin(apfa[1][z4][z8],apfa[2][z4][z7]);

												if (o->m_pCDF->m_bWriteSnapshots) {

													if ((apfa[0][z4][z5] >= o->m_pCDF->m_faWriteSnapshotsRange[0]) &&
														(apfa[0][z4][z5] <= o->m_pCDF->m_faWriteSnapshotsRange[1]) &&
														(apfa[1][z4][z8] >= o->m_pCDF->m_faWriteSnapshotsRange[2]) &&
														(apfa[1][z4][z8] <= o->m_pCDF->m_faWriteSnapshotsRange[3]) &&
														(apfa[2][z4][z7] >= o->m_pCDF->m_faWriteSnapshotsRange[4]) &&
														(apfa[2][z4][z7] <= o->m_pCDF->m_faWriteSnapshotsRange[5])) {

														if (!o->m_pCDF->m_bWriteSnapshotsTempRM) {
															o->m_pCDF->m_bWriteSnapshotsTempRM = true;
															o->m_pCDF->m_iWriteSnapshotsTotalPerRM++;
														}

														o->m_pCDF->m_iWriteSnapshotsTotalCounter++;

														if (o->m_pCDF->m_bWriteSnapshotsNoWrite)
															goto _nocdfwritesnapshot3;

														if (o->m_pCDF->m_iWriteSnapshotStride > 1) {
															o->m_pCDF->m_iWriteSnapshotCounter++;
															if (o->m_pCDF->m_iWriteSnapshotCounter % o->m_pCDF->m_iWriteSnapshotStride != 0)
																goto _nocdfwritesnapshot3;
														}

														bufcom.sprintf(
															"Step %lu, RM=%d, OM=%d, Val1=%f, Val2=%f, Val3=%f",
															g_iSteps,
															smfix->m_iMolSMIndex+1,
															z4+1,
															apfa[0][z4][z5],
															apfa[1][z4][z8],
															apfa[2][z4][z7]
														);

														if (o->m_pCDF->m_bWriteSnapshotsOnlyRMOM) {

															tiasnap.clear();
															tiasnap.push_back( smfix->m_iIndex );
															tiasnap.push_back( ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex[z4] );

															g_pTempTimestepSnap->ExportXYZConfiguration(
																o->m_pCDF->m_pWriteSnapshotsFile,
																(const char*)bufcom,
																o->m_pCDF->m_bWriteSnapshotsOriginalCoords,
																&tiasnap
															);

															o->m_pCDF->m_iWriteSnapshotsTotalWritten++;

														} else {

															if (!o->m_pCDF->m_bWriteSnapshotDone) {

																// Only write a snapshot once per trajectory frame, as all molecules are included
																o->m_pCDF->m_bWriteSnapshotDone = true;

																g_pTempTimestep->ExportXYZConfiguration(
																	o->m_pCDF->m_pWriteSnapshotsFile,
																	(const char*)bufcom,
																	o->m_pCDF->m_bWriteSnapshotsOriginalCoords,
																	NULL
																);

																o->m_pCDF->m_iWriteSnapshotsTotalWritten++;
															}
														}
													}
												}
_nocdfwritesnapshot3:

_3nocdfbin:
												if (o->m_pCDF->m_bDumpDat)
													mfprintf(o->m_pCDF->m_fDump,"%d;  %d;  %d;  ;  %10.3f;  %10.3f;  %10.3f\n",(int)g_iSteps,z2+1,z4+1,apfa[0][z4][z5],apfa[1][z4][z8],apfa[2][z4][z7]);
											}
										}
									}
								}
							} // END IF CHANNELS == 3
						} // IF NOT SECOND_OM
					} // IF CDF

					for (zr=0;zr<g_iCDFChannels;zr++)
					{
						if (g_bRDF && !(o->m_bObsCertain && o->m_bDecompDist) && (!o->m_bDecompType) && (o->m_pRDF[zr] != NULL) && (!g_bDeriv || (g_iSteps > 2)))
						{
							if (o->m_bOthers)
							{
								if (o->m_bSecondShowMol && (zr == 1))
								{
									for (z4=0;z4<o->m_iShowMol2Count;z4++)
									{
										for (z5=0;z5<o->m_pRDF[zr]->m_faData[z4].GetSize();z5++)
										{
											if (o->m_bBinOnlyPassedAtoms)
												if (o->m_pRDF[zr]->m_baDataEnabled[z4][z5] == 0)
													continue;
											if (o->m_bBinOnlyNotPassedAtoms)
												if (o->m_pRDF[zr]->m_baDataEnabled[z4][z5] != 0)
													continue;
											o->m_pRDF[zr]->m_pRDF->AddToBin(o->m_pRDF[zr]->m_faData[z4][z5]);
										}
									}
								} else
								{
									for (z4=0;z4<o->m_iShowMolCount;z4++)
									{
										for (z5=0;z5<o->m_pRDF[zr]->m_faData[z4].GetSize();z5++)
										{
											if (o->m_bBinOnlyPassedAtoms)
												if (o->m_pRDF[zr]->m_baDataEnabled[z4][z5] == 0)
													continue;
											if (o->m_bBinOnlyNotPassedAtoms)
												if (o->m_pRDF[zr]->m_baDataEnabled[z4][z5] != 0)
													continue;
											o->m_pRDF[zr]->m_pRDF->AddToBin(o->m_pRDF[zr]->m_faData[z4][z5]);
										}
									}
								}
							} else
							{
								for (z5=0;z5<o->m_pRDF[zr]->m_faData[0].GetSize();z5++)
									o->m_pRDF[zr]->m_pRDF->AddToBin(o->m_pRDF[zr]->m_faData[0][z5]);
							}
							if (o->m_bTimeDev && (ti2 != -1))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pRDF[zr]->m_fDist[ti2],"\n");
										else if (ti2 == (int)o->m_waSaveRefList.GetSize()-1)
											mfprintf(o->m_pRDF[zr]->m_fDist[0],"\n");
							}
						}

						if (g_bVDF && !(o->m_bObsCertain && o->m_bDecompDist) && (o->m_pVDF[zr] != NULL) && (!g_bDeriv || (g_iSteps > 2)))
						{
							if (o->m_bOthers)
							{
								if (o->m_bSecondShowMol && (zr == 1))
								{
									for (z4=0;z4<o->m_iShowMol2Count;z4++)
										for (z5=0;z5<o->m_pVDF[zr]->m_faData[z4].GetSize();z5++)
											o->m_pVDF[zr]->m_pVDF->AddToBin(o->m_pVDF[zr]->m_faData[z4][z5]);
								} else
								{
									for (z4=0;z4<o->m_iShowMolCount;z4++)
									{
										for (z5=0;z5<o->m_pVDF[zr]->m_faData[z4].GetSize();z5++)
										{
											if (o->m_bBinOnlyPassedAtoms)
												if (o->m_pVDF[zr]->m_baDataEnabled[z4][z5] == 0)
													continue;
											if (o->m_bBinOnlyNotPassedAtoms)
												if (o->m_pVDF[zr]->m_baDataEnabled[z4][z5] != 0)
													continue;
											o->m_pVDF[zr]->m_pVDF->AddToBin(o->m_pVDF[zr]->m_faData[z4][z5]);
										}
									}
								}
							} else
							{
								for (z5=0;z5<o->m_pVDF[zr]->m_faData[0].GetSize();z5++)
									o->m_pVDF[zr]->m_pVDF->AddToBin(o->m_pVDF[zr]->m_faData[0][z5]);
							}
							if (o->m_bTimeDev && (ti2 != -1))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pVDF[zr]->m_fSpeed[ti2],"\n");
										else if (ti2 == (int)o->m_waSaveRefList.GetSize()-1)
											mfprintf(o->m_pVDF[zr]->m_fSpeed[0],"\n");
							}
						}

						if (g_bDipDF && !(o->m_bObsCertain && o->m_bDecompDist) && (o->m_pDipDF[zr] != NULL) && (!g_bDeriv || (g_iSteps > 2)))
						{
							if (o->m_bOthers)
							{
								if (o->m_bSecondShowMol && (zr == 1))
								{
									for (z4=0;z4<o->m_iShowMol2Count;z4++)
										for (z5=0;z5<o->m_pDipDF[zr]->m_faData[z4].GetSize();z5++)
											o->m_pDipDF[zr]->m_pDipoleDF->AddToBin(o->m_pDipDF[zr]->m_faData[z4][z5]);
								} else
								{
									for (z4=0;z4<o->m_iShowMolCount;z4++)
									{
										for (z5=0;z5<o->m_pDipDF[zr]->m_faData[z4].GetSize();z5++)
										{
											if (o->m_bBinOnlyPassedAtoms)
												if (o->m_pDipDF[zr]->m_baDataEnabled[z4][z5] == 0)
													continue;
											if (o->m_bBinOnlyNotPassedAtoms)
												if (o->m_pDipDF[zr]->m_baDataEnabled[z4][z5] != 0)
													continue;
											o->m_pDipDF[zr]->m_pDipoleDF->AddToBin(o->m_pDipDF[zr]->m_faData[z4][z5]);
										}
									}
								}
							} else
							{
								for (z5=0;z5<o->m_pDipDF[zr]->m_faData[0].GetSize();z5++)
									o->m_pDipDF[zr]->m_pDipoleDF->AddToBin(o->m_pDipDF[zr]->m_faData[0][z5]);
							}
							if (o->m_bTimeDev && (ti2 != -1))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pDipDF[zr]->m_fDipole[ti2],"\n");
										else if (ti2 == (int)o->m_waSaveRefList.GetSize()-1)
											mfprintf(o->m_pDipDF[zr]->m_fDipole[0],"\n");
							}
						}

						if (g_bADF && !(o->m_bObsCertain && o->m_bDecompDist) && (o->m_pADF[zr] != NULL) && (!g_bDeriv || (g_iSteps > 2)))
						{
							if (o->m_bOthers)
							{
								if (o->m_bSecondShowMol && (zr == 1))
								{
									for (z4=0;z4<o->m_iShowMol2Count;z4++)
										for (z5=0;z5<o->m_pADF[zr]->m_faData[z4].GetSize();z5++)
											o->m_pADF[zr]->m_pADF->AddToBin(o->m_pADF[zr]->m_faData[z4][z5]);
								} else
								{
									for (z4=0;z4<o->m_iShowMolCount;z4++)
									{
										for (z5=0;z5<o->m_pADF[zr]->m_faData[z4].GetSize();z5++)
										{
											if (o->m_bBinOnlyPassedAtoms)
												if (o->m_pADF[zr]->m_baDataEnabled[z4][z5] == 0)
													continue;
											if (o->m_bBinOnlyNotPassedAtoms)
												if (o->m_pADF[zr]->m_baDataEnabled[z4][z5] != 0)
													continue;
											o->m_pADF[zr]->m_pADF->AddToBin(o->m_pADF[zr]->m_faData[z4][z5]);
										}
									}
								}
							} else
							{
								for (z5=0;z5<o->m_pADF[zr]->m_faData[0].GetSize();z5++)
									o->m_pADF[zr]->m_pADF->AddToBin(o->m_pADF[zr]->m_faData[0][z5]);
							}
							if (o->m_bTimeDev && (ti2 != -1))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pADF[zr]->m_fAngle[ti2],"\n");
										else if (ti2 == (int)o->m_waSaveRefList.GetSize()-1)
											mfprintf(o->m_pADF[zr]->m_fAngle[0],"\n");
							}
						}

						if (g_bDDF && !(o->m_bObsCertain && o->m_bDecompDist) && (o->m_pDDF[zr] != NULL) && (!g_bDeriv || (g_iSteps > 2)))
						{
							if (o->m_bOthers)
							{
								if (o->m_bSecondShowMol && (zr == 1))
								{
									for (z4=0;z4<o->m_iShowMol2Count;z4++)
										for (z5=0;z5<o->m_pDDF[zr]->m_faData[z4].GetSize();z5++)
											o->m_pDDF[zr]->m_pDDF->AddToBin(o->m_pDDF[zr]->m_faData[z4][z5]);
								} else
								{
									for (z4=0;z4<o->m_iShowMolCount;z4++)
									{
										for (z5=0;z5<o->m_pDDF[zr]->m_faData[z4].GetSize();z5++)
										{
											if (o->m_bBinOnlyPassedAtoms)
												if (o->m_pDDF[zr]->m_baDataEnabled[z4][z5] == 0)
													continue;
											if (o->m_bBinOnlyNotPassedAtoms)
												if (o->m_pDDF[zr]->m_baDataEnabled[z4][z5] != 0)
													continue;
											o->m_pDDF[zr]->m_pDDF->AddToBin(o->m_pDDF[zr]->m_faData[z4][z5]);
										}
									}
								}
							} else
							{
								for (z5=0;z5<o->m_pDDF[zr]->m_faData[0].GetSize();z5++)
									o->m_pDDF[zr]->m_pDDF->AddToBin(o->m_pDDF[zr]->m_faData[0][z5]);
							}
							if (o->m_bTimeDev && (ti2 != -1))
							{
								if (o->m_bSaveSeparateFiles)
									mfprintf(o->m_pDDF[zr]->m_fAngle[ti2],"\n");
										else if (ti2 == (int)o->m_waSaveRefList.GetSize()-1)
											mfprintf(o->m_pDDF[zr]->m_fAngle[0],"\n");
							}
						} // End IF DDF

						if (g_bPlDF && !(o->m_bObsCertain && o->m_bDecompDist) && (o->m_pPlDF[zr] != NULL) && (!g_bDeriv || (g_iSteps > 2)))
						{
							if (o->m_bOthers)
							{
								if (o->m_bSecondShowMol && (zr == 1))
								{
									for (z4=0;z4<o->m_iShowMol2Count;z4++)
										for (z5=0;z5<o->m_pPlDF[zr]->m_faData[z4].GetSize();z5++)
											o->m_pPlDF[zr]->m_pPlDF->AddToBin(o->m_pPlDF[zr]->m_faData[z4][z5]);
								} else
								{
									for (z4=0;z4<o->m_iShowMolCount;z4++)
									{
										for (z5=0;z5<o->m_pPlDF[zr]->m_faData[z4].GetSize();z5++)
										{
											if (o->m_bBinOnlyPassedAtoms)
												if (o->m_pPlDF[zr]->m_baDataEnabled[z4][z5] == 0)
													continue;
											if (o->m_bBinOnlyNotPassedAtoms)
												if (o->m_pPlDF[zr]->m_baDataEnabled[z4][z5] != 0)
													continue;
											o->m_pPlDF[zr]->m_pPlDF->AddToBin(o->m_pPlDF[zr]->m_faData[z4][z5]);
										}
									}
								}
							} else
							{
								for (z5=0;z5<o->m_pPlDF[zr]->m_faData[0].GetSize();z5++)
									o->m_pPlDF[zr]->m_pPlDF->AddToBin(o->m_pPlDF[zr]->m_faData[0][z5]);
							}
						} // End IF PlDF

						if (g_bLiDF && !(o->m_bObsCertain && o->m_bDecompDist) && (o->m_pLiDF[zr] != NULL) && (!g_bDeriv || (g_iSteps > 2)))
						{
							if (o->m_bOthers)
							{
								if (o->m_bSecondShowMol && (zr == 1))
								{
									for (z4=0;z4<o->m_iShowMol2Count;z4++)
										for (z5=0;z5<o->m_pLiDF[zr]->m_faData[z4].GetSize();z5++)
											o->m_pLiDF[zr]->m_pLiDF->AddToBin(o->m_pLiDF[zr]->m_faData[z4][z5]);
								} else
								{
									for (z4=0;z4<o->m_iShowMolCount;z4++)
									{
										for (z5=0;z5<o->m_pLiDF[zr]->m_faData[z4].GetSize();z5++)
										{
											if (o->m_bBinOnlyPassedAtoms)
												if (o->m_pLiDF[zr]->m_baDataEnabled[z4][z5] == 0)
													continue;
											if (o->m_bBinOnlyNotPassedAtoms)
												if (o->m_pLiDF[zr]->m_baDataEnabled[z4][z5] != 0)
													continue;
											o->m_pLiDF[zr]->m_pLiDF->AddToBin(o->m_pLiDF[zr]->m_faData[z4][z5]);
										}
									}
								}
							} else
							{
								for (z5=0;z5<o->m_pLiDF[zr]->m_faData[0].GetSize();z5++)
									o->m_pLiDF[zr]->m_pLiDF->AddToBin(o->m_pLiDF[zr]->m_faData[0][z5]);
							}
						} // End IF LiDF

					} // FOR CDFChannels

					if (g_bCond && o->m_bCondDevelopment) {

						if (o->m_pConditions->m_bAnyPassed)
							o->m_iaCondDevelopmentCounter[0]++;

						for (z4=1;z4<=o->m_iCondDevelopmentMax;z4++)
							if (o->m_pConditions->m_iTempPassed == z4)
								o->m_iaCondDevelopmentCounter[z4]++;

						o->m_iCondDevelopmentAvg += o->m_pConditions->m_iTempPassed;
					}

					if (o->m_bNormalizeCondition)
						if (o->m_iNormalizeConditionLocalCount != 0)
							o->m_iNormalizeConditionCount++;

				} // Ende FOR g_oaObserv.GetSize

			} // Ende IF [blabla]	
	
			if (g_bMiddleAvg) // Referenzmolekuel mitteln
			{
				cc = 0;
				ti2 = (g_iSteps-2)*((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()+z2;
				// Jeden Atomtyp des Zielmolekuels durchgehen
				for (z3=0;z3<((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex.GetSize();z3++)
				{
					for (z4=0;z4<((CMolecule*)g_oaMolecules[g_iFixMol])->m_waAtomCount[z3];z4++)
					{
						// Atome vertauschen, aber nicht im allerersten Durchgang hier, und nur wenn mehr als 1 Atom dieser Sorte
						if ((g_iSwapAtoms) && (ti2 != 0) && (((CMolecule*)g_oaMolecules[g_iFixMol])->m_waAtomCount[z3]>1))
						{
							vec4 = g_pRefMol[cc];
							vec4 /= (double)ti2; // der bisherige Mittelwert des Referenzmolekuels
							tf = VecDist(g_pTempTimestep->m_vaCoords[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[z2]])->m_oaAtomOffset[z3])->GetAt(z4)],vec4); // der Abstand dieses Atoms von seinem Aequivalent im Referenzmolekuel
							ti3 = -1;
							for (z5=z4+1;z5<((CMolecule*)g_oaMolecules[g_iFixMol])->m_waAtomCount[z3];z5++) // Alle folgenden Atome dieser Sorte durchgehen
							{
								tf2 = VecDist(g_pTempTimestep->m_vaCoords[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[z2]])->m_oaAtomOffset[z3])->GetAt(z5)],vec4);
								if (tf2 < tf) // Das andere Atom ist naeher am Platz im Referenzmolekuel als das eigentlich vorgesehene
								{
									ti3 = z5;
									tf = tf2;
								}
							}
							if (ti3 != -1) // Ein Anderes ist naeher dran als unseres: Vertausche diese beiden
							{
								Swap(g_pTempTimestep->m_vaCoords[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[z2]])->m_oaAtomOffset[z3])->GetAt(z4)],g_pTempTimestep->m_vaCoords[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[z2]])->m_oaAtomOffset[z3])->GetAt(ti3)]);
								pSwapMatrix[cc*((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes+(cc-z4+ti3)] += 1;
								pSwapMatrix[(cc-z4+ti3)*((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes+cc] += 1;
							}
						} // Ende IF SwapAtoms
						g_pRefMol[cc] += g_pTempTimestep->m_vaCoords[((CxIntArray*)((CSingleMolecule*)g_oaSingleMolecules[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex[z2]])->m_oaAtomOffset[z3])->GetAt(z4)];
						cc++;
					} // Ende FOR Atom des Referenzmolekuels
				} // Ende FOR Atomtyp des Referenzmolekuels
			} // Ende IF Referenzmolekuel mitteln
		} // Ende FOR RefMol (z2)

		if (g_bCond) {

			for (z6=0;z6<g_oaObserv.GetSize();z6++) {

				o = (CObservation*)g_oaObserv[z6];

				if (o->m_bCondDevelopment) {

					mfprintf(o->m_pCondDevelopment,"%lu;  %d;  %.6f;  %.6f",
						g_iSteps,
						o->m_iaCondDevelopmentCounter[0],
						(double)o->m_iaCondDevelopmentCounter[0]/((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*100.0,
						(double)o->m_iCondDevelopmentAvg/((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()
						);

					for (z2=1;z2<=o->m_iCondDevelopmentMax;z2++)
						mfprintf(o->m_pCondDevelopment,";  %d;  %.6f",
							o->m_iaCondDevelopmentCounter[z2],
							(double)o->m_iaCondDevelopmentCounter[z2]/((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*100.0
							);

					mfprintf(o->m_pCondDevelopment,"\n");
				}
			}
		}


		for (z6=0;z6<g_oaObserv.GetSize();z6++)
		{
			o = (CObservation*)g_oaObserv[z6];

			if (o->m_pConditions != NULL)
				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();z2++)
					if (o->m_pConditions->m_iPassCounter[z2] != 0)
						o->m_pConditions->m_iOMPassCounter[z2]++;

			if (o->m_pConditionsOM2 != NULL)
				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();z2++)
					if (o->m_pConditionsOM2->m_iPassCounter[z2] != 0)
						o->m_pConditionsOM2->m_iOMPassCounter[z2]++;

			if (o->m_bPercTimeDev) {
				if (o->m_iPercTimeDevTempTotal != 0)
					tf = (double)o->m_iPercTimeDevTempCounter/o->m_iPercTimeDevTempTotal*100.0;
				else
					tf = 0;
				mfprintf(
					o->m_pPercTimeDevFile,
					"%lu;  %d;  %.10f\n",
					g_iSteps,
					o->m_iPercTimeDevTempCounter,
					tf
				);
			}

		}

//_end:;
_endstep:
_norefmol:
		if ((g_iMaxStep > 0) && ((int)g_iSteps >= g_iMaxStep))
		{
			mprintf("\n\nMaximum step count of %d reached, stopping.",g_iMaxStep);
			break;
		}
	} // Ende while

/*************************************************************************
*************** Ende Hauptschleife ***************************************
*************************************************************************/

//	fclose(fff);

_endmainloop:
	if (g_bUseVelocities || g_bUseForces)
		g_iSteps -= 2;

	mprintf(WHITE,"\n\n##########   Main pass done   ##########\n");

	if (g_iStride != 1) {
		if (g_iStride == 2)
			mprintf("\n%lu time steps processed (every %dnd from %lu).\n\n",g_iSteps/g_iStride,g_iStride,g_iSteps);
		else if (g_iStride == 3)
			mprintf("\n%lu time steps processed (every %drd from %lu).\n\n",g_iSteps/g_iStride,g_iStride,g_iSteps);
		else
			mprintf("\n%lu time steps processed (every %dth from %lu).\n\n",g_iSteps/g_iStride,g_iStride,g_iSteps);
	} else
		mprintf("\n%lu time steps processed.\n\n",g_iSteps);

//	fclose(g_fPos);
	CloseInputTrajectory();

	if (g_fVel != NULL)
		fclose(g_fVel);
	if (g_fForce != NULL)
		fclose(g_fForce);
	if ((g_bNPT) && (g_sNPTFile[0] != 0))
		fclose(g_fNPTFile);

	if (g_bDipole && g_bDumpDipoleVector)
	{
		fclose(g_fDumpDipole);
		mprintf("    Dipole vectors written to file \"dipole_vectors.csv\".\n");
		if (g_bDumpDipoleXYZ)
		{
			fclose(g_fDumpDipoleXYZ);
			mprintf("    Dipole XYZ trajectory written to file \"dipole_vectors.xyz\".\n");
		}
		mprintf("\n");
	}

	mprintf(WHITE,"\n\n    Finishing analyses:\n");
	mprintf(WHITE,"----------------------------------------------------------------------------------------------------------\n\n");

	/************ Interface **************/
	Interface_AfterAnalysis();

	if (g_bFixedPlProj)
		g_pFixedPlProj->Finish();

	if (g_bDomA)
		g_pDomainEngine->Finish();

	if (g_bVoro)
		g_pVoroWrapper->Finish();


	g_bAbortAnalysis = true;

	if (g_bSaveCondSnapshot)
	{
		fclose(g_fSaveCondFile);
		mprintf("%d condition snapshots saved as savecondition.xyz\n\n",g_iSaveCondCount);
	}

	if (g_bOrder)
		g_pOrderEngine->Finish();

	if (g_bAggrTopo)
		g_pAggrTopoEngine->Finish();

	if (g_bContactMatrix)
		g_pContactMatrix->Finish();

	if (g_bTDDF)
		g_pTDDFEngine->Finish();

	if (g_bGeoDens)
		g_pGeoDensEngine->Finish();

	if (g_bClusterAnalysis)
	{
		mprintf(WHITE,"*** Cluster Analysis\n\n");
		g_pClusterAnalysis->BuildClusterDistribution();
		g_pClusterAnalysis->WriteOutput(multibuf);
	}


	if (g_bBondACF)
	{
		mprintf(WHITE,"*** Bond autocorrelation function\n");
		if (g_bBondACFDebug)
		{
			mprintf("    Writing \"bondacf.agr\"...\n");

			try { gc = new CGrace(); } catch(...) { gc = NULL; }
			if (gc == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			gc->SetRangeX(0,g_iSteps*g_fTimestepLength/1000.0);
			gc->SetRangeY(0,250);
			gc->MakeTicks();
			gc->SetLabelX("Time / ps");
			gc->SetLabelY("Bond length / pm");
			for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[z2];
				for (z3=0;z3<sm->m_oaBonds.GetSize();z3++)
				{
					bond = (CMolBond*)sm->m_oaBonds[z3];
					gc->AddDataset();
					for (z4=0;z4<bond->m_faData.GetSize();z4++)
						gc->AddXYTupel(z4*g_fTimestepLength/1000.0,bond->m_faData[z4]);
				}
			}
			gc->WriteAgr("bondacf.agr",false);
			delete gc;
		}
		mprintf("    Differentiating bond length developments...\n");
		for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[z2];
			for (z3=0;z3<sm->m_oaBonds.GetSize();z3++)
			{
				bond = (CMolBond*)sm->m_oaBonds[z3];
				for (z4=0;z4<bond->m_faData.GetSize()-1;z4++)
					bond->m_faData[z4] = bond->m_faData[z4+1] - bond->m_faData[z4];
			}
		}
		if (g_bBondACFDebug)
		{
			mprintf("    Writing \"bondacf_diff.agr\"...\n");

			try { gc = new CGrace(); } catch(...) { gc = NULL; }
			if (gc == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			gc->SetRangeX(0,g_iSteps*g_fTimestepLength/1000.0);
			gc->SetRangeY(-15,15);
			gc->MakeTicks();
			gc->SetLabelX("Time / ps");
			gc->SetLabelY("Bond length change rate");
			for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[z2];
				for (z3=0;z3<sm->m_oaBonds.GetSize();z3++)
				{
					bond = (CMolBond*)sm->m_oaBonds[z3];
					gc->AddDataset();
					for (z4=0;z4<bond->m_faData.GetSize();z4++)
						gc->AddXYTupel(z4*g_fTimestepLength/1000.0,bond->m_faData[z4]);
				}
			}
			gc->WriteAgr("bondacf_diff.agr",false);
			delete gc;
		}
		mprintf("    Autocorrelating bond length developments...\n");
		mprintf(WHITE,"      [");
		tfs = 0;
		for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[z2];
			tfs += sm->m_oaBonds.GetSize();
		}
		tfs /= 60.0;
		ti2 = 0;
/*		fft = new CFFT();
		fft->PrepareFFT_C2C(g_iSteps*2);
		fft2 = new CFFT();
		fft2->PrepareInverseFFT_C2C(g_iSteps*2);*/

		try { ac = new CAutoCorrelation(); } catch(...) { ac = NULL; }
		if (ac == NULL) NewException((double)sizeof(CAutoCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		ac->Init(g_iSteps,g_iBondACFDepth,g_bACFFFT);
		for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++)
		{
			sm = (CSingleMolecule*)g_oaSingleMolecules[z2];
			for (z3=0;z3<sm->m_oaBonds.GetSize();z3++)
			{
				if (fmod(ti2,tfs) < 1.0)
					mprintf(WHITE,"#");
				bond = (CMolBond*)sm->m_oaBonds[z3];
				ac->AutoCorrelate(&bond->m_faData,&tempfa);
				bond->m_faData.CopyFrom(&tempfa);
				ti2++;
			}
		}
		delete ac;
		mprintf(WHITE,"]\n");
		if (g_bBondACFDebug)
		{
			mprintf("    Writing \"bondacf_autocorr.agr\"...\n");

			try { gc = new CGrace(); } catch(...) { gc = NULL; }
			if (gc == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			gc->SetRangeX(0,g_iBondACFDepth*g_fTimestepLength/1000.0);
			gc->SetRangeY(-25,25);
			gc->MakeTicks();
			gc->SetLabelX("Time / ps");
			gc->SetLabelY("ACF(Bond length change rate)");
			for (z2=0;z2<g_oaSingleMolecules.GetSize();z2++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[z2];
				for (z3=0;z3<sm->m_oaBonds.GetSize();z3++)
				{
					bond = (CMolBond*)sm->m_oaBonds[z3];
					gc->AddDataset();
					for (z4=0;z4<g_iBondACFDepth;z4++)
						gc->AddXYTupel(z4*g_fTimestepLength/1000.0,bond->m_faData[z4]);
				}
			}
			gc->WriteAgr("bondacf_autocorr.agr",false);
			delete gc;
		}
		mprintf("    Merging equivalent Bonds...\n");
		for (z0=0;z0<g_oaMolecules.GetSize();z0++)
		{
			m = (CMolecule*)g_oaMolecules[z0];
			sm2 = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
			for (z2=0;z2<m->m_laSingleMolIndex.GetSize();z2++)
			{
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z2]];
				for (z3=0;z3<sm->m_oaBondGroups.GetSize();z3++)
				{
					bg = (CMolBondGroup*)sm->m_oaBondGroups[z3];
					for (z4=1;z4<bg->m_oaBonds.GetSize();z4++)
					{
						bond = (CMolBond*)bg->m_oaBonds[z4];
						for (z5=0;z5<g_iBondACFDepth;z5++)
							((CMolBond*)((CMolBondGroup*)sm2->m_oaBondGroups[z3])->m_oaBonds[0])->m_faData[z5] += bond->m_faData[z5];
					}
				}
			}
		}
		if (g_bBondACFSymmetrize)
		{
			mprintf("    Symmetrizing bond ACFs...\n");
			for (z0=0;z0<g_oaMolecules.GetSize();z0++)
			{
				m = (CMolecule*)g_oaMolecules[z0];
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
				for (z3=0;z3<sm->m_oaBondGroups.GetSize();z3++)
				{
					bond = (CMolBond*)((CMolBondGroup*)sm->m_oaBondGroups[z3])->m_oaBonds[0];
					bond->m_faData.SetSize(2*g_iBondACFDepth);
					for (z4=0;z4<g_iBondACFDepth;z4++)
						bond->m_faData[z4+g_iBondACFDepth] = bond->m_faData[z4];
					for (z4=0;z4<g_iBondACFDepth;z4++)
						bond->m_faData[z4] = bond->m_faData[2*g_iBondACFDepth-z4-1];
				}
			}
			g_iBondACFDepth *= 2;
		}
		if (g_bBondACFWindow)
		{
			mprintf("    Applying Window Function...\n");
			for (z0=0;z0<g_oaMolecules.GetSize();z0++)
			{
				m = (CMolecule*)g_oaMolecules[z0];
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
				for (z3=0;z3<sm->m_oaBondGroups.GetSize();z3++)
				{
					bond = (CMolBond*)((CMolBondGroup*)sm->m_oaBondGroups[z3])->m_oaBonds[0];
					for (z4=0;z4<g_iBondACFDepth;z4++)
						bond->m_faData[z4] *= (double)pow2(sin(z4*Pi/(bond->m_faData.GetSize()-1)));
				}
			}
		}
		if (g_bBondACFNormalize)
		{
			mprintf("    Normalizing bond ACFs...\n");
			for (z0=0;z0<g_oaMolecules.GetSize();z0++)
			{
				m = (CMolecule*)g_oaMolecules[z0];
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
				for (z3=0;z3<sm->m_oaBondGroups.GetSize();z3++)
				{
					bond = (CMolBond*)((CMolBondGroup*)sm->m_oaBondGroups[z3])->m_oaBonds[0];
					tf = 0;
					for (z4=0;z4<g_iBondACFDepth;z4++)
						if (tf < bond->m_faData[z4])
							tf = bond->m_faData[z4];
					for (z4=0;z4<g_iBondACFDepth;z4++)
						bond->m_faData[z4] /= (double)tf;
				}
			}
		}
		mprintf("    Writing \"bondacf_ac_merged.agr\"...\n");

		try { gc = new CGrace(); } catch(...) { gc = NULL; }
		if (gc == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		gc->SetRangeX(0,g_iBondACFDepth*g_fTimestepLength/1000.0);
		gc->SetRangeY(-200,200);
		gc->MakeTicks();
		gc->SetLabelX("Time / ps");
		gc->SetLabelY("Sum ACF(Bond length change rate)");
		gc->CurrentGraph()->m_bLegend = true;
		for (z0=0;z0<g_oaMolecules.GetSize();z0++)
		{
			m = (CMolecule*)g_oaMolecules[z0];
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
			for (z3=0;z3<sm->m_oaBondGroups.GetSize();z3++)
			{
				bond = (CMolBond*)((CMolBondGroup*)sm->m_oaBondGroups[z3])->m_oaBonds[0];
				gc->AddDataset();
//				sprintf(buf,"%s %s%d - %s%d",m->m_sName,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[0]]])->m_sName,bond->m_iAtom[0]+1,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[1]]])->m_sName,bond->m_iAtom[1]+1);
				buf.sprintf("%s %s%d - %s%d",m->m_sName,(const char*)((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[0]]])->m_sName,bond->m_iAtom[0]+1,(const char*)((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[1]]])->m_sName,bond->m_iAtom[1]+1);
				gc->SetDatasetName(buf);
				((CGraceDataset*)gc->CurrentGraph()->m_oaDatasets[gc->CurrentGraph()->m_oaDatasets.GetSize()-1])->m_faValues.SetMaxSize(g_iBondACFDepth*2);
				for (z4=0;z4<g_iBondACFDepth;z4++)
					gc->AddXYTupel(z4*g_fTimestepLength/1000.0,bond->m_faData[z4]);
			}
		}
		gc->WriteAgr("bondacf_ac_merged.agr",false);
		delete gc;
		mprintf("    Applying Fourier Transformation...\n");

		try { fft = new CFFT(); } catch(...) { fft = NULL; }
		if (fft == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		fft->PrepareFFT_C2C(g_iBondACFDepth);
		tf3 = 0;
		for (z0=0;z0<g_oaMolecules.GetSize();z0++)
		{
			m = (CMolecule*)g_oaMolecules[z0];
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
			for (z3=0;z3<sm->m_oaBondGroups.GetSize();z3++)
			{
				bond = (CMolBond*)((CMolBondGroup*)sm->m_oaBondGroups[z3])->m_oaBonds[0];
				for (z4=0;z4<g_iBondACFDepth;z4++)
				{
					fft->m_pInput[z4*2] = bond->m_faData[z4];
					fft->m_pInput[z4*2+1] = 0;
				}
				fft->DoFFT();
				for (z4=0;z4<g_iBondACFDepth/2;z4++)
				{
					bond->m_faData[z4] = (double)(pow2(fft->m_pOutput[z4*2]) + pow2(fft->m_pOutput[z4*2+1])) / g_iBondACFDepth;
					if (bond->m_faData[z4] > tf3)
						tf3 = bond->m_faData[z4];
				}
			}
		}
		delete fft;
		mprintf("    Writing \"bond_spectrum.agr\"...\n");

		try { gc = new CGrace(); } catch(...) { gc = NULL; }
		if (gc == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
		
		gc->SetRangeX(0,4000.0);
		gc->SetRangeY(0,tf3*1.1);
		gc->MakeTicks();
		gc->SetLabelX("Wave number [1/cm]");
		gc->SetLabelY("Intensity");
		gc->CurrentGraph()->m_bInvertXAxis = true;
		gc->CurrentGraph()->m_bLegend = true;
		for (z0=0;z0<g_oaMolecules.GetSize();z0++)
		{
			m = (CMolecule*)g_oaMolecules[z0];
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
			for (z3=0;z3<sm->m_oaBondGroups.GetSize();z3++)
			{
				bond = (CMolBond*)((CMolBondGroup*)sm->m_oaBondGroups[z3])->m_oaBonds[0];
				gc->AddDataset();
//				sprintf(buf,"%s %s%d - %s%d",m->m_sName,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[0]]])->m_sName,bond->m_iAtom[0]+1,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[1]]])->m_sName,bond->m_iAtom[1]+1);
				buf.sprintf("%s %s%d - %s%d",m->m_sName,(const char*)((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[0]]])->m_sName,bond->m_iAtom[0]+1,(const char*)((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[1]]])->m_sName,bond->m_iAtom[1]+1);
				gc->SetDatasetName(buf);
				((CGraceDataset*)gc->CurrentGraph()->m_oaDatasets[gc->CurrentGraph()->m_oaDatasets.GetSize()-1])->m_faValues.SetMaxSize(g_iBondACFDepth*2);
				for (z4=0;z4<g_iBondACFDepth/2;z4++)
					gc->AddXYTupel(z4*2.0/g_iBondACFDepth*1E7/299.792/2.0/g_fTimestepLength/g_iStride,bond->m_faData[z4]);
			}
		}
		gc->WriteAgr("bond_spectrum.agr",false);
		delete gc;
		for (z5=1;z5<5;z5++)
		{
			mprintf("    Smoothing spectrum, degree %d...\n",z5);
			mprintf("    Writing \"bond_spectrum_s%d.agr\"...\n",z5);

			try { gc = new CGrace(); } catch(...) { gc = NULL; }
			if (gc == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			gc->SetRangeX(0,4000.0);
			gc->SetRangeY(0,tf3*1.1);
			gc->MakeTicks();
			gc->SetLabelX("Wave number [1/cm]");
			gc->SetLabelY("Intensity");
			gc->CurrentGraph()->m_bInvertXAxis = true;
			gc->CurrentGraph()->m_bLegend = true;

			try { gc2 = new CGrace(); } catch(...) { gc2 = NULL; }
			if (gc2 == NULL) NewException((double)sizeof(CGrace),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			gc2->SetRangeX(0,4000.0);
			gc2->SetRangeY(0,tf3*1.1);
			gc2->MakeTicks();
			gc2->SetLabelX("Wave number [1/cm]");
			gc2->SetLabelY("Intensity");
			gc2->CurrentGraph()->m_bInvertXAxis = true;
			gc2->CurrentGraph()->m_bLegend = true;
			
			tempfa.SetSize(g_iBondACFDepth/2);
			tempfa2.SetSize(g_iBondACFDepth/2);

			for (z4=0;z4<g_iBondACFDepth/2;z4++)
				tempfa2[z4] = 0;

			for (z0=0;z0<g_oaMolecules.GetSize();z0++)
			{
				m = (CMolecule*)g_oaMolecules[z0];
				sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
				for (z3=0;z3<sm->m_oaBondGroups.GetSize();z3++)
				{
					bond = (CMolBond*)((CMolBondGroup*)sm->m_oaBondGroups[z3])->m_oaBonds[0];

					for (z4=0;z4<g_iBondACFDepth/2;z4++)
					{
						tf = 0;
						tf2 = 0;
						for (z6=-z5;z6<=z5;z6++)
						{
							if ((z4+z6 < 0) || (z4+z6 >= g_iBondACFDepth/2))
								continue;
							tf += bond->m_faData[z4+z6] / (pow2((double)z6)+1.0);
							tf2 += 1.0 / (pow2((double)z6)+1.0);
						}
						tempfa[z4] = (double)(tf / tf2);
						tempfa2[z4] += tempfa[z4];
					}

					gc->AddDataset();
					gc2->AddDataset();
//					sprintf(buf,"%s %s%d - %s%d",m->m_sName,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[0]]])->m_sName,bond->m_iAtom[0]+1,((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[1]]])->m_sName,bond->m_iAtom[1]+1);
					buf.sprintf("%s %s%d - %s%d",m->m_sName,(const char*)((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[0]]])->m_sName,bond->m_iAtom[0]+1,(const char*)((CAtom*)g_oaAtoms[sm->m_baAtomIndex[bond->m_iAtomType[1]]])->m_sName,bond->m_iAtom[1]+1);
					gc->SetDatasetName(buf);
					gc2->SetDatasetName(buf);
					((CGraceDataset*)gc->CurrentGraph()->m_oaDatasets[gc->CurrentGraph()->m_oaDatasets.GetSize()-1])->m_faValues.SetMaxSize(g_iBondACFDepth);
					((CGraceDataset*)gc2->CurrentGraph()->m_oaDatasets[gc2->CurrentGraph()->m_oaDatasets.GetSize()-1])->m_faValues.SetMaxSize(g_iBondACFDepth);

					for (z4=0;z4<g_iBondACFDepth/2;z4++)
						gc->AddXYTupel(z4*2.0/g_iBondACFDepth*1E7/299.792/2.0/g_fTimestepLength/g_iStride,tempfa[z4]);

					for (z4=0;z4<g_iBondACFDepth/2;z4++)
						gc2->AddXYTupel(z4*2.0/g_iBondACFDepth*1E7/299.792/2.0/g_fTimestepLength/g_iStride,tempfa2[z4]);
				}
			}
	
//			sprintf(buf,"bond_spectrum_s%d.agr",z5);
			buf.sprintf("bond_spectrum_s%d.agr",z5);
			gc->WriteAgr(buf,false);
			delete gc;

			mprintf("    Writing \"bond_spectrum_cumulative_s%d.agr\"...\n",z5);
//			sprintf(buf,"bond_spectrum_cumulative_s%d.agr",z5);
			buf.sprintf("bond_spectrum_cumulative_s%d.agr",z5);
			gc2->WriteAgr(buf,false);
			delete gc2;
		}
	}

	if (g_bSaveRefEnv)
	{
		fclose(g_fRefEnv);
		mprintf("\nTrajectory of reference molecule was saved as %s.\n",(const char*)g_sRefEnv);
		if (g_bTDO)
		{
#ifdef TARGET_WINDOWS
//			sprintf(buf,"tdo_pymol_%s_%d%s.bat",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,multibuf);
			buf.sprintf("tdo_pymol_%s_%d%s.bat",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,multibuf);
#else
//			sprintf(buf,"tdo_pymol_%s_%d%s.sh",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,multibuf);
			buf.sprintf("tdo_pymol_%s_%d%s.sh",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,multibuf);
#endif

			mprintf("Saving TDO PyMol script as %s.\n",(const char*)buf);
			tfi = OpenFileWrite(buf,true);
#ifdef TARGET_WINDOWS
			mfprintf(tfi,"pymolwin ");
#else
			mfprintf(tfi,"pymol ");
#endif
			for (z=0;z<g_laTDOSteps.GetSize();z++)
				mfprintf(tfi,"tdo_%s_%d_%06d%s.xyz ",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,g_laTDOSteps[z],multibuf);
			mfprintf(tfi,"-d \"set sphere_scale,0.25\" ");
			mfprintf(tfi,"-d \"set stick_radius,0.17\" ");
			mfprintf(tfi,"-d \"set ray_trace_mode,1\" ");
			mfprintf(tfi,"-d \"set fog,0\" ");
			mfprintf(tfi,"-d \"set ray_trace_fog,0\" ");
			mfprintf(tfi,"-d \"set bg_rgb,(1,1,1)\" ");
			mfprintf(tfi,"-d \"set ray_shadow,0\" ");
			mfprintf(tfi,"-d \"set ray_shadows,0\" ");
			mfprintf(tfi,"-d \"set ray_interior_shadows,0\" ");
			mfprintf(tfi,"-d \"show sticks; show spheres\" ");
			
			for (z=0;z<g_laTDOSteps.GetSize();z++)
			{
				for (z2=0;z2<g_oaAtoms.GetSize();z2++)
				{
					if ((z2 == g_iVirtAtomType) && (!g_bSaveVirtAtoms))
						continue;
					mfprintf(tfi,"-d \"set_color mol_%s_%d, [%d,%d,%d]\" ",(const char*)((CAtom*)g_oaAtoms[z2])->m_sName,z+1,((CAtom*)g_oaAtoms[z2])->m_pElement->ColorR((double)z/(g_laTDOSteps.GetSize()-1)*g_fTDOBleaching),((CAtom*)g_oaAtoms[z2])->m_pElement->ColorG((double)z/(g_laTDOSteps.GetSize()-1)*g_fTDOBleaching),((CAtom*)g_oaAtoms[z2])->m_pElement->ColorB((double)z/(g_laTDOSteps.GetSize()-1)*g_fTDOBleaching));
				}
			}
			for (z=0;z<g_laTDOSteps.GetSize();z++)
			{
				for (z2=0;z2<g_oaAtoms.GetSize();z2++)
				{
					if ((z2 == g_iVirtAtomType) && (!g_bSaveVirtAtoms))
						continue;
					mfprintf(tfi,"-d \"color mol_%s_%d, tdo_%s_%d_%06d and name %s\" ",(const char*)((CAtom*)g_oaAtoms[z2])->m_sName,z+1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,g_laTDOSteps[z],(const char*)((CAtom*)g_oaAtoms[z2])->m_sName);
				}
			}
			mfprintf(tfi,"\n");
			fclose(tfi);
#ifdef TARGET_LINUX
			buf.sprintf("chmod 755 tdo_pymol_%s_%d%s.sh",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,g_iSaveRefMol+1,multibuf);
			mprintf("Executing \"%s\"...\n",(const char*)buf);
			(void)!system(buf);
#endif
		}
	}
	if (g_bCutCluster)
	{
		fclose(g_fRefEnv);
		mprintf("\nCluster trajectory was saved as %s.\n",(const char*)g_sRefEnv);
	}

	if (g_bSaveJustTraj)
	{
		fclose(g_fSaveJustTraj);
/*		strcpy(buf,g_sInputTraj);
		p = strrchr(buf,'.');
		*p = 0;
		strcat(buf,multibuf);
		strcat(buf,".out.xyz");*/
//		mprintf("\nProcessed output trajectory was saved as traj_out.xyz\n",buf);
	}

	if (g_bIRSpec && g_bGlobalIR)
	{
		mprintf(WHITE,"\n*** Global IR Spectrum\n");
		g_pGlobalIR->Finish(multibuf);
	}

	if (g_bVACF && g_bGlobalVACF)
	{
		mprintf(WHITE,"\n*** Global velocity autocorrelation function\n");
		if (g_pGlobalVACF->m_bDerivative)
		{
			mprintf("    Deriving velocities...\n");
			mprintf(WHITE,"      [");
			tfs = g_iGesAtomCount/60.0;
			ti = 0;
			for (z2=0;z2<g_iGesAtomCount;z2++)
			{
				if (fmod(z2,tfs) < 1)
					mprintf(WHITE,"#");

				if (g_pGlobalVACF->m_bExcludeR0 && (g_waAtomRealElement[z2] > 1000))
					continue;

				ptfa = (CxDoubleArray*)g_pGlobalVACF->m_oaCache[ti];

				for (z3=0;z3<(int)g_iSteps*3-3;z3++)
					(*ptfa)[z3] = (*ptfa)[z3+3] - (*ptfa)[z3];

				ti++;
			}
			mprintf(WHITE,"]\n");
		}
		if (g_bVACFCacheMode)
		{
			tfs = g_iGesAtomCount/60.0;
/*			if (g_bACFFFT)
			{*/
				mprintf("    Autocorrelating cached vectors...\n");
				mprintf(WHITE,"      [");
/*				fft = new CFFT();
				fft->PrepareFFT_C2C(2*g_iSteps);
				fft2 = new CFFT();
				fft2->PrepareInverseFFT_C2C(2*g_iSteps);*/

				try { ptfa2 = new CxDoubleArray("main():ptfa2"); } catch(...) { ptfa2 = NULL; }
				if (ptfa2 == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				ptfa2->SetSize(g_iSteps);

				try { ptfa3 = new CxDoubleArray("main():ptfa3"); } catch(...) { ptfa3 = NULL; }
				if (ptfa3 == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				ptfa3->SetSize(g_iSteps);

				try { ac = new CAutoCorrelation(); } catch(...) { ac = NULL; }
				if (ac == NULL) NewException((double)sizeof(CAutoCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				ac->Init(g_iSteps,g_pGlobalVACF->m_iSize,g_bACFFFT);
				ti = 0;
				for (z2=0;z2<g_iGesAtomCount;z2++)
				{
					if (fmod(z2,tfs) < 1.0)
						mprintf(WHITE,"#");

					if (g_pGlobalVACF->m_bExcludeR0 && (g_waAtomRealElement[z2] > 1000))
						continue;

					ptfa = (CxDoubleArray*)g_pGlobalVACF->m_oaCache[ti];

					if (g_pGlobalVACF->m_bMassWeight)
						tf = ((CAtom*)g_oaAtoms[g_waAtomRealElement[z2]])->m_pElement->m_fMass;
							else tf = 1.0;

					/* X */
					for (z3=0;z3<(int)g_iSteps;z3++)
						(*ptfa2)[z3] = (*ptfa)[z3*3];
					ac->AutoCorrelate(ptfa2,ptfa3);
		//			mprintf("Global Atom %d X: %f + %f = %f.\n",z2+1,g_pGlobalVACF->m_pData[0],(*ptfa3)[0]*tf,g_pGlobalVACF->m_pData[0]+(*ptfa3)[0]*tf);
					for (z3=0;z3<(int)g_pGlobalVACF->m_iSize;z3++)
						g_pGlobalVACF->m_pData[z3] += (*ptfa3)[z3] * tf;

					/* Y */
					for (z3=0;z3<(int)g_iSteps;z3++)
						(*ptfa2)[z3] = (*ptfa)[z3*3+1];
					ac->AutoCorrelate(ptfa2,ptfa3);
		//			mprintf("Global Atom %d Y: %f + %f = %f.\n",z2+1,g_pGlobalVACF->m_pData[0],(*ptfa3)[0]*tf,g_pGlobalVACF->m_pData[0]+(*ptfa3)[0]*tf);
					for (z3=0;z3<(int)g_pGlobalVACF->m_iSize;z3++)
						g_pGlobalVACF->m_pData[z3] += (*ptfa3)[z3] * tf;

					/* Z */
					for (z3=0;z3<(int)g_iSteps;z3++)
						(*ptfa2)[z3] = (*ptfa)[z3*3+2];
					ac->AutoCorrelate(ptfa2,ptfa3);
		//			mprintf("Global Atom %d Z: %f + %f = %f.\n",z2+1,g_pGlobalVACF->m_pData[0],(*ptfa3)[0]*tf,g_pGlobalVACF->m_pData[0]+(*ptfa3)[0]*tf);
					for (z3=0;z3<(int)g_pGlobalVACF->m_iSize;z3++)
						g_pGlobalVACF->m_pData[z3] += (*ptfa3)[z3] * tf;

					ti++;
				}
				delete ac;
				delete ptfa2;
				delete ptfa3;
	/*			delete fft2;
				delete fft;*/
/*			} else
			{
				mprintf("    Autocorrelating cached vectors...\n");
				mprintf(WHITE,"      [");
				for (z2=0;z2<g_iGesAtomCount;z2++)
				{
					if (fmod(z2,tfs) < 1.0)
						mprintf(WHITE,"#");
					ptfa = (CxFloatArray*)g_pGlobalVACF->m_oaCache[z2];
					for (z3=0;z3<g_pGlobalVACF->m_iSize;z3+=g_pGlobalVACF->m_iStride) // Das ist das Tau
					{
						tf = 0;
						for (z4=0;z4<(int)g_iSteps-z3;z4++) // Das ist der Startpunkt
							tf += (*ptfa)[z4*3]*(*ptfa)[(z4+z3)*3] + (*ptfa)[z4*3+1]*(*ptfa)[(z4+z3)*3+1] + (*ptfa)[z4*3+2]*(*ptfa)[(z4+z3)*3+2];
						g_pGlobalVACF->m_pData[z3/g_pGlobalVACF->m_iStride] += tf/(g_iSteps-z3);
					}
				}
			}*/
			mprintf(WHITE,"]\n");
			mprintf("      %d atoms, %lu time steps and %d correlation depths processed.\n",g_pGlobalVACF->m_oaCache.GetSize(),g_iSteps,g_pGlobalVACF->m_iSize);
			g_pGlobalVACF->MultiplyCached(1.0/g_iGesAtomCount);

			if (g_pGlobalVACF->m_bDecomposeModes)                                                                                                                                                                                                                               
			{                                                                                                                                                                                                                                                                   
				tfs = g_iGesAtomCount*g_iGesAtomCount*3/60.0;                                                                                                                                                                                                               
                                                                                                                                                                                                                                                                                           
				mprintf("    Allocating cross-correlation matrix (%s)...\n",FormatBytes((double)sizeof(double)*g_pGlobalVACF->m_iSize*g_pGlobalVACF->m_iParticles*3*g_pGlobalVACF->m_iParticles*3));                                                                         
                                                                                                                                                                                                                                                                                           
				g_pGlobalVACF->m_oaCCRMatrix.SetMaxSize(g_pGlobalVACF->m_iParticles*3*g_pGlobalVACF->m_iParticles*3);                                                                                                                                                       
				for (z2=0;z2<g_pGlobalVACF->m_iParticles*3;z2++)                                                                                                                                                                                                            
				{                                                                                                                                                                                                                                                           
					for (z3=0;z3<g_pGlobalVACF->m_iParticles*3;z3++)                                                                                                                                                                                                    
					{                                                                                                                                                                                                                                                   
			//			mprintf("X: %X\n",g_pGlobalVACF);                                                                                                                                                                                                           
			//			mprintf("%d - %d -->\n",z2,z3);                                                                                                                                                                                                             
						try { ptfa = new CxDoubleArray("main():ptfa"); } catch(...) { ptfa = NULL; }                                                                                                                                                                              
						if (ptfa == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);                                                                                                                                         
                                                                                                                                                                                                                                                                                           
			//			mprintf("A\n");                                                                                                                                                                                                                             
						ptfa->SetSize(g_pGlobalVACF->m_iSize);                                                                                                                                                                                                      
			//			mprintf("B\n");                                                                                                                                                                                                                             
			//			mprintf("X2: %X\n",g_pGlobalVACF);                                                                                                                                                                                                          
						g_pGlobalVACF->m_oaCCRMatrix.Add(ptfa);                                                                                                                                                                                                     
			//			mprintf("%d - %d <--\n",z2,z3);                                                                                                                                                                                                             
					}                                                                                                                                                                                                                                                   
				}                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                           
				mprintf("    Computing %d cross-correlations...\n",g_pGlobalVACF->m_iParticles*3*g_pGlobalVACF->m_iParticles*3);                                                                                                                                            
				mprintf(WHITE,"      [");                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                           
				try { ptfa2 = new CxDoubleArray("main():ptfa2"); } catch(...) { ptfa2 = NULL; }                                                                                                                                                                                            
				if (ptfa2 == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);                                                                                                                                                        
				                                                                                                                                                                                                                                                            
				ptfa2->SetSize(g_iSteps);                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                           
				try { ptfa2b = new CxDoubleArray("main():ptfa2b"); } catch(...) { ptfa2b = NULL; }                                                                                                                                                                                          
				if (ptfa2b == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);                                                                                                                                                       
				                                                                                                                                                                                                                                                            
				ptfa2b->SetSize(g_iSteps);                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                           
				try { ccr = new CCrossCorrelation(); } catch(...) { ccr = NULL; }                                                                                                                                                                                           
				if (ccr == NULL) NewException((double)sizeof(CCrossCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);                                                                                                                                                     
				                                                                                                                                                                                                                                                            
				ccr->Init(g_iSteps,g_pGlobalVACF->m_iSize,g_bACFFFT);                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                           
				ti = 0;                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                           
				for (z2=0;z2<g_iGesAtomCount;z2++)                                                                                                                                                                                                                          
				{                                                                                                                                                                                                                                                           
					if (g_pGlobalVACF->m_bExcludeR0 && (g_waAtomRealElement[z2] > 1000))                                                                                                                                                                                
						continue;                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                           
					ptfa = (CxDoubleArray*)g_pGlobalVACF->m_oaCache[ti];                                                                                                                                                                                                 
                                                                                                                                                                                                                                                                                           
					for (z2b=0;z2b<3;z2b++)                                                                                                                                                                                                                             
					{                                                                                                                                                                                                                                                   
						ti2 = 0;                                                                                                                                                                                                                                    
						for (z3=0;z3<g_iGesAtomCount;z3++)                                                                                                                                                                                                          
						{                                                                                                                                                                                                                                           
							if (fmod(z2*g_iGesAtomCount*3+z2b*g_iGesAtomCount+z3,tfs) < 1.0)                                                                                                                                                                    
								mprintf(WHITE,"#");                                                                                                                                                                                                         
                                                                                                                                                                                                                                                                                           
							if (g_pGlobalVACF->m_bExcludeR0 && (g_waAtomRealElement[z3] > 1000))                                                                                                                                                                
								continue;                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                           
							ptfab = (CxDoubleArray*)g_pGlobalVACF->m_oaCache[ti2];                                                                                                                                                                               
                                                                                                                                                                                                                                                                                           
							if (g_pGlobalVACF->m_bMassWeight)                                                                                                                                                                                                   
								tf = ((CAtom*)g_oaAtoms[g_waAtomRealElement[z2]])->m_pElement->m_fMass;                                                                                                                                                     
									else tf = 1.0;                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                           
							for (z3b=0;z3b<3;z3b++)                                                                                                                                                                                                             
							{                                                                                                                                                                                                                                   
								for (z4=0;z4<(int)g_iSteps;z4++)                                                                                                                                                                                            
								{                                                                                                                                                                                                                           
									(*ptfa2)[z4] = (*ptfa)[z4*3+z2b];                                                                                                                                                                                   
									(*ptfa2b)[z4] = (*ptfab)[z4*3+z3b];                                                                                                                                                                                 
								}                                                                                                                                                                                                                           
								ccr->CrossCorrelate(ptfa2,ptfa2b,(CxDoubleArray*)g_pGlobalVACF->m_oaCCRMatrix[ti*g_pGlobalVACF->m_iParticles*9 + z2b*g_pGlobalVACF->m_iParticles*3 + ti2*3 + z3b]);                                                          
					//			mprintf("z2=%d, z2b=%d, z3=%d, z3b=%d, I=%d\n",z2,z2b,z3,z3b,ti*g_pGlobalVACF->m_iParticles*9 + z2b*g_pGlobalVACF->m_iParticles*3 + ti2*3 + z3b);                                                                           
							}                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                           
							ti2++;                                                                                                                                                                                                                              
						}                                                                                                                                                                                                                                           
					}                                                                                                                                                                                                                                                   
					ti++;                                                                                                                                                                                                                                               
				}                                                                                                                                                                                                                                                           
				delete ccr;                                                                                                                                                                                                                                                 
				delete ptfa2;                                                                                                                                                                                                                                               
				delete ptfa2b;                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                           
				mprintf(WHITE,"]\n");                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                           
		/*		for (z2=0;z2<g_pGlobalVACF->m_iParticles*3;z2++)                                                                                                                                                                                                            
				{                                                                                                                                                                                                                                                           
					for (z3=0;z3<g_pGlobalVACF->m_iParticles*3;z3++)                                                                                                                                                                                                    
						mprintf("%6G  ",((CxFloatArray*)g_pGlobalVACF->m_oaCCRMatrix[z2*g_pGlobalVACF->m_iParticles*3+z3])->GetAt(0));                                                                                                                              
					mprintf("\n");                                                                                                                                                                                                                                      
				}*/                                                                                                                                                                                                                                                         
                                                                                                                                                                                                                                                                                           
				Interface_DecomposeModes(g_pGlobalVACF->m_iParticles*3, &g_pGlobalVACF->m_oaCCRMatrix);                                                                                                                                                                     
                                                                                                                                                                                                                                                                                           
			} // END IF m_bDecomposeModes                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                           
		} else g_pGlobalVACF->Multiply(1.0/g_iGesAtomCount);

		mprintf("    Saving global VACF as acf_global%s.csv ...\n",multibuf);
		g_pGlobalVACF->WriteACF("acf_global",multibuf,".csv");

		if (g_pGlobalVACF->m_iMirror != 0)
		{
			mprintf("    Mirroring global VACF...\n");
			g_pGlobalVACF->Mirror(g_pGlobalVACF->m_iMirror);
			mprintf("    Saving mirrored global VACF as acf_global.m%s.csv ...\n",multibuf);
			g_pGlobalVACF->WriteACF("acf_global.m",multibuf,".csv");
		}

		if (g_pGlobalVACF->m_bWindowFunction)
		{
			mprintf("    Applying window function...\n");
			g_pGlobalVACF->Window();
			mprintf("    Saving windowed global VACF as acf_global.w%s.csv ...\n",multibuf);
			g_pGlobalVACF->WriteACF("acf_global.w",multibuf,".csv");
		}

		if (g_pGlobalVACF->m_bSpectrum)
		{
			mprintf("    Performing Fourier transformation...\n");

			try { g_pFFT = new CFFT(); } catch(...) { g_pFFT = NULL; }
			if (g_pFFT == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			g_pFFT->PrepareFFT_C2C(g_pGlobalVACF->m_iSize+g_pGlobalVACF->m_iZeroPadding);
			g_pGlobalVACF->Transform(g_pFFT);
			delete g_pFFT;
			g_pGlobalVACF->m_pSpectrum->SetMaxRWL(1E7/299.792/g_fTimestepLength/g_iStride);
			if (g_pGlobalVACF->m_bACF_DB)
			{
				mprintf("    Normalising spectrum to decibel...\n");
				g_pGlobalVACF->m_pSpectrum->MakeDB();
			}/* else
			{
				mprintf("    Normalising integral of spectrum...\n");
				g_pGlobalVACF->m_pSpectrum->SetIntegral(1000000.0f);
			}*/
			if (g_pGlobalVACF->m_bWindowFunction)
			{
				mprintf("    Saving spectrum as power_global_w%s.csv ...\n",multibuf);
				g_pGlobalVACF->m_pSpectrum->Write("power_global_w",multibuf,".csv");
			} else
			{
				mprintf("    Saving spectrum as power_global%s.csv ...\n",multibuf);
				g_pGlobalVACF->m_pSpectrum->Write("power_global",multibuf,".csv");
			}
		}
		mprintf("\n");
	}

	if (g_bVFDF)
	{
		for (z=0;z<g_iVFCorrCount;z++)
			fclose(g_fVFCorr[z]);
	}


//	fclose(ff);

	for (z=0;z<g_oaObserv.GetSize();z++)
	{
		o = (CObservation*)g_oaObserv[z];
		mprintf(YELLOW,"\n>>> Observation %d >>>\n",z+1);

		if (o->m_pConditions != NULL)
		{
			if (o->m_bSecondShowMol)
				mprintf(WHITE,"\n#### Condition between RM and 1st OM ####\n");

			o->m_pConditions->PrintData();
			if (o->m_pConditions->m_oaConditionSubGroups.GetSize()==2)
				o->m_pConditions->PrintTable();
		}

		if ((o->m_bSecondShowMol) && (o->m_pConditionsOM2 != NULL))
		{
			mprintf(WHITE,"\n\n#### Condition between RM and 2nd OM ####\n");
			o->m_pConditionsOM2->PrintData();
			if (o->m_pConditionsOM2->m_oaConditionSubGroups.GetSize()==2)
				o->m_pConditionsOM2->PrintTable();
		}

		if (g_bCond && o->m_bCondDevelopment) {
			buf.sprintf("cond_%d_%s_%s.csv",z+1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_sName);
			mprintf("\n    Closing Condition Temporal Development File \"%s\"...\n",(const char*)buf);
			fclose(o->m_pCondDevelopment);
		}

		if (o->m_bPercTimeDev) {
			buf.sprintf("fraction_timedev_%d_%s.csv",z+1,((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
			mprintf("\n    Closing Fraction Temporal Development File \"%s\"...\n",(const char*)buf);
			fclose(o->m_pPercTimeDevFile);
			o->m_pPercTimeDevFile = NULL;
		}

		if (g_bAggregation)
		{
			mprintf(WHITE,"\n* Aggregation Functions\n");
			if ((g_iMaxStep > 0) && (((int)g_iSteps*g_iStride) >= g_iMaxStep))
				g_pTempTimestep->CopyFrom(GetTimeStep(0)); // Max. Schrittzahl
					else g_pTempTimestep->CopyFrom(GetTimeStep(1)); // End Of File

			if (g_bDACF) {
				if (o->m_pDACF->m_iDACFRes >= (int)g_iSteps) {
					mprintf("\n");
					mprintf(RED,"    Warning: ");
					mprintf("Requested ACF depth (%d) is larger than number of processed time steps (%lu).\n",o->m_pDACF->m_iDACFRes,g_iSteps);
					mprintf("             Reducing ACF depth to %lu.\n\n",g_iSteps-1);
					o->m_pDACF->m_iDACFRes = (int)g_iSteps-1;
					for (z2=0;z2<o->m_pDACF->m_oaSubDACFs.GetSize();z2++) {
						((CDACFSub*)o->m_pDACF->m_oaSubDACFs[z2])->m_pDACF->m_iResolution = o->m_pDACF->m_iDACFRes;
						((CDACFSub*)o->m_pDACF->m_oaSubDACFs[z2])->m_pDACF->m_fMaxVal = o->m_pDACF->m_iDACFRes * g_fTimestepLength / 1000.0;
					}
				}
			}

			for (zs=0;zs<o->m_pDACF->m_oaSubDACFs.GetSize();zs++)
			{
				dacfsub = (CDACFSub*)o->m_pDACF->m_oaSubDACFs[zs];
				mprintf("\n###########################################################################################################\n\n");
				mprintf(WHITE,"  > Value Set %d: %s\n",zs+1,dacfsub->m_sName);

				o->m_pDACF->FinishDACFSub(g_pTempTimestep,dacfsub);

//				sprintf(buf,"cond_%s.txt",dacfsub->m_sName);
				buf.sprintf("cond_%s.txt",dacfsub->m_sName);
				dacfsub->m_pCondition->PrintData(buf);

				mprintf(WHITE,"Neighbour Count Distribution\n");
				mprintf("    %.0f Bin entries have been made.\n",dacfsub->m_pNDF->m_fBinEntries);
//				sprintf(buf,"ncd_%s%s.csv",dacfsub->m_sName,multibuf);
				buf.sprintf("ncd_%s%s.csv",dacfsub->m_sName,multibuf);
				mprintf("    Writing Neighbor Count Distribution File \"%s\"...\n",(const char*)buf);
				dacfsub->m_pNDF->Write_Int("",buf,"");

				if (g_bDDisp)
				{
					mprintf(WHITE,"\nDimer Displacement Distribution Function\n");
					mprintf("    %.0f Bin entries have been made (%.0f skipped, %.0f total).\n",dacfsub->m_pDDisp->m_fBinEntries,dacfsub->m_pDDisp->m_fSkipEntries,dacfsub->m_pDDisp->m_fBinEntries+dacfsub->m_pDDisp->m_fSkipEntries);
//					sprintf(buf,"ddisp_%s%s.csv",dacfsub->m_sName,multibuf);
					buf.sprintf("ddisp_%s%s.csv",dacfsub->m_sName,multibuf);
					mprintf("    Saving DDisp File as \"%s\"...\n",(const char*)buf);
					dacfsub->m_pDDisp->Write("",buf,"",false);
//					sprintf(buf,"ddisp_%s%s.agr",dacfsub->m_sName,multibuf);
					buf.sprintf("ddisp_%s%s.agr",dacfsub->m_sName,multibuf);
					mprintf("    Saving DDisp Agr File as \"%s\"...\n",(const char*)buf);
					dacfsub->m_pDDisp->WriteAgr("",buf,"","",false);
				}

				if (g_bDACF)
				{
					mprintf(WHITE,"\nDimer Existence Autocorrelation Function\n");
		//			mprintf("    %.0f Bin entries have been made (%.0f skipped, %.0f total).\n",dacfsub->m_pDACF->m_fBinEntries,dacfsub->m_pDACF->m_fSkipEntries,dacfsub->m_pDACF->m_fBinEntries+dacfsub->m_pDACF->m_fSkipEntries);

		/*			if (dacfsub->m_bBorderMode)
					{
						mprintf("    *** New Border Mode ***\n");
						for (z2=0;z2<o->m_pDACF->m_iDACFRes;z2++)
							dacfsub->m_pDACF->m_pBin[z2] /= g_iSteps-z2;
					} else
					{
						for (z2=0;z2<o->m_pDACF->m_iDACFRes;z2++)
							dacfsub->m_pDACF->m_pBin[z2] /= g_iSteps;
					}*/


					if (dacfsub->m_bIntermittend)
					{
						if (dacfsub->m_bIntNew) {
							mprintf("    Building intermittent function (new ACF-less method)...\n");
							// "Cheating": Expand bin to trajectory length
							delete[] dacfsub->m_pDACF->m_pBin;
							dacfsub->m_pDACF->m_pBin = new double[g_iSteps];
							mprintf(WHITE,"      [");
							tfs = ((CMolecule*)g_oaMolecules[o->m_pDACF->m_iFirstMol])->m_laSingleMolIndex.GetSize() * ((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize() / 60.0;
							for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_pDACF->m_iFirstMol])->m_laSingleMolIndex.GetSize();z2++)
							{
								for (z3=0;z3<((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize();z3++)
								{
									if (fmod(z2*((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+z3,tfs) < 1.0)
										mprintf(WHITE,"#");

									if (dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+z3].GetSize() == 0)
										continue;

									dacfsub->ProcessPair(dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+z3]);
								}
							}
							mprintf(WHITE,"]\n");
							for (z2=0;z2<o->m_pDACF->m_iDACFRes;z2++)
								dacfsub->m_pDACF->m_pBin[z2] /= (double)(g_iSteps-z2);
						} else {

							mprintf("    Autocorrelating reconstructed time series...\n");
							if (dacfsub->m_bDebugOut)
								mprintf("      Writing Debug output to \"dacf_%s_debugout_*.csv...\n",dacfsub->m_sName);
							mprintf(WHITE,"      [");
							tfs = ((CMolecule*)g_oaMolecules[o->m_pDACF->m_iFirstMol])->m_laSingleMolIndex.GetSize() * ((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize() / 60.0;
							ptfa = new CxDoubleArray();
							ptfa2 = new CxDoubleArray();
							ac = new CAutoCorrelation();
							ac->Init(g_iSteps,o->m_pDACF->m_iDACFRes,true);

							ptfa->SetSize(g_iSteps);
							ptfa2->SetSize(o->m_pDACF->m_iDACFRes);

							for (z2=0;z2<o->m_pDACF->m_iDACFRes;z2++)
								dacfsub->m_pDACF->m_pBin[z2] = 0;

							ti = 0;
							ti2 = 0;
							for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_pDACF->m_iFirstMol])->m_laSingleMolIndex.GetSize();z2++)
							{
								for (z3=0;z3<((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize();z3++)
								{
									if (fmod(z2*((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+z3,tfs) < 1.0)
										mprintf(WHITE,"#");

									if (dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+z3].GetSize() == 0)
										continue;

							//		mprintf("*** %d <-> %d ***\n",z2,z3);

									for (z5=0;z5<(int)g_iSteps;z5++)
										(*ptfa)[z5] = 0;

									ti++;
									ti3 = 0;
									for (z4=0;z4<dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+z3].GetSize();z4+=2) {
										ti2++;
										ti3++;
										for (z5=dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+z3].GetAt(z4);z5<=dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+z3].GetAt(z4+1);z5++)
											(*ptfa)[z5-1] = 1;
									}

							//		if (dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+z3].GetSize() != 0)
							//			mprintf("\n%d-%d: %d ",z2,z3,dacfsub->m_piaIntervals[z2*((CMolecule*)g_oaMolecules[o->m_pDACF->m_iSecondMol])->m_laSingleMolIndex.GetSize()+z3].GetSize()/2);

									ac->AutoCorrelate(ptfa,ptfa2);

									if (dacfsub->m_bDebugOut) {
										if (ti3 != 0) {
											FILE *tfo;
											CxString tds;
											tds.sprintf("dacf_%s_debugout_%d_%d.csv",dacfsub->m_sName,z2,z3);
											tfo = OpenFileWrite((const char*)tds,true);
											mfprintf(tfo,"#Step;  TSeries;  ACF\n");
											for (z4=0;z4<(int)g_iSteps;z4++) {
												mfprintf(tfo,"%d;  %f",z4,(*ptfa)[z4]);
												if (z4 < o->m_pDACF->m_iDACFRes)
													mfprintf(tfo,";  %f",(*ptfa2)[z4]);
												else
													mfprintf(tfo,";  0");
												mfprintf(tfo,"\n");
											}
											fclose(tfo);
										}
									}

							/*		if ((z2==0) && (z3==28))
									{
										FILE *tfia;

										sprintf(buf,"condfunc_%d_%d.csv",z2,z3);
										tfia = fopen(buf,"wt");
										for (z4=0;z4<(int)g_iSteps;z4++)
											fprintf(tfia,"%d;  %f\n",z4,(*ptfa)[z4]);
										fclose(tfia);

										sprintf(buf,"condacf_%d_%d.csv",z2,z3);
										tfia = fopen(buf,"wt");
										for (z4=0;z4<o->m_pDACF->m_iDACFRes;z4++)
											fprintf(tfia,"%d;  %f\n",z4,(*ptfa2)[z4]);
										fclose(tfia);
									}*/
				
									for (z4=0;z4<o->m_pDACF->m_iDACFRes;z4++)
										dacfsub->m_pDACF->m_pBin[z4] += (*ptfa2)[z4];
								}
							}
							mprintf(WHITE,"]\n");
							mprintf("    %d pairs and %d intervals processed.\n",ti,ti2);
						}
					}

					if (dacfsub->m_bCorrectEq)
					{
						if (dacfsub->m_iRefMol == dacfsub->m_iShowMol)
						{
							tf = 100.0*dacfsub->m_fEqCounter/g_iSteps/(((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize()-1)/((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize();
							mprintf("    Finite-size ensemble average is %.6f %c.\n",tf,'%');

							tf = 1.0/(((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize()-1)/((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize();
							dacfsub->m_pDACF->MultiplyBin(tf);

							tf = dacfsub->m_fEqCounter/g_iSteps/(((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize()-1)/((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize();
							tf = tf*tf;
							mprintf("    Subtracting squared average value (%G - %G)...\n",dacfsub->m_pDACF->m_pBin[0],tf);
							dacfsub->m_pDACF->SubtractBin(tf);
						} else
						{
							tf = 100.0*dacfsub->m_fEqCounter/g_iSteps/((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize()/((CMolecule*)g_oaMolecules[dacfsub->m_iShowMol])->m_laSingleMolIndex.GetSize();
							mprintf("    Finite-size ensemble average is %.6f %c.\n",tf,'%');

							tf = 1.0/((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize()/((CMolecule*)g_oaMolecules[dacfsub->m_iShowMol])->m_laSingleMolIndex.GetSize();
							dacfsub->m_pDACF->MultiplyBin(tf);

							tf = dacfsub->m_fEqCounter/g_iSteps/((CMolecule*)g_oaMolecules[dacfsub->m_iRefMol])->m_laSingleMolIndex.GetSize()/((CMolecule*)g_oaMolecules[dacfsub->m_iShowMol])->m_laSingleMolIndex.GetSize();
							tf = tf*tf;
							mprintf("    Subtracting squared average value (%G - %G)...\n",dacfsub->m_pDACF->m_pBin[0],tf);
							dacfsub->m_pDACF->SubtractBin(tf);
						}
					}

					if (dacfsub->m_pDACF->m_pBin[0] != 0)
					{
						mprintf("    Normalizing ACF with factor %G ...\n",1.0/dacfsub->m_pDACF->m_pBin[0]);
						dacfsub->m_pDACF->MultiplyBin(1.0/dacfsub->m_pDACF->m_pBin[0]);
					} else
						mprintf("    Warning: Not normalizing this DACF, first entry is 0.\n");

					if (o->m_pDACF->m_bFitDACF)
					{
						mprintf("\n");

						try { dacfsub->m_pDACF->m_pParameters = new double*[o->m_pDACF->m_iFitDegreeMax+1]; } catch(...) { dacfsub->m_pDACF->m_pParameters = NULL; }
						if (dacfsub->m_pDACF->m_pParameters == NULL) NewException((double)(o->m_pDACF->m_iFitDegreeMax+1)*sizeof(double*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						try { dacfsub->m_pDACF->m_pFitIntegral = new double[o->m_pDACF->m_iFitDegreeMax+1]; } catch(...) { dacfsub->m_pDACF->m_pFitIntegral = NULL; }
						if (dacfsub->m_pDACF->m_pFitIntegral == NULL) NewException((double)(o->m_pDACF->m_iFitDegreeMax+1)*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						try { dacfsub->m_pDACF->m_pCorrCoeff = new double[o->m_pDACF->m_iFitDegreeMax+1]; } catch(...) { dacfsub->m_pDACF->m_pCorrCoeff = NULL; }
						if (dacfsub->m_pDACF->m_pCorrCoeff == NULL) NewException((double)(o->m_pDACF->m_iFitDegreeMax+1)*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						try { dacfsub->m_pDACF->m_pAdditionalSets = new double*[o->m_pDACF->m_iFitDegreeMax+1]; } catch(...) { dacfsub->m_pDACF->m_pAdditionalSets = NULL; }
						if (dacfsub->m_pDACF->m_pAdditionalSets == NULL) NewException((double)(o->m_pDACF->m_iFitDegreeMax+1)*sizeof(double*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						try { dacfsub->m_pDACF->m_pAdditionalSetLabels = new char*[o->m_pDACF->m_iFitDegreeMax+1]; } catch(...) { dacfsub->m_pDACF->m_pAdditionalSetLabels = NULL; }
						if (dacfsub->m_pDACF->m_pAdditionalSetLabels == NULL) NewException((double)(o->m_pDACF->m_iFitDegreeMax+1)*sizeof(char*),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						dacfsub->m_pDACF->m_iAdditionalSets = o->m_pDACF->m_iFitDegreeMax+1;

						for (z2=0;z2<=o->m_pDACF->m_iFitDegreeMax;z2++)
						{
							dacfsub->m_pDACF->m_pAdditionalSets[z2] = NULL;
							dacfsub->m_pDACF->m_pAdditionalSetLabels[z2] = NULL;
						}

						mprintf(YELLOW,"    Please note: ");
						mprintf("The lifetime is two times the integral value!\n\n");

						for (z2=o->m_pDACF->m_iFitDegreeMin;z2<=o->m_pDACF->m_iFitDegreeMax;z2++)
							dacfsub->m_pDACF->Fit_PolyExp(z2,100);

						if (o->m_pDACF->m_iFitDegreeMin != o->m_pDACF->m_iFitDegreeMax) {
							mprintf(YELLOW,"    DACF fit integral overview:\n");
							mprintf("      (use the row with the highest R value)\n\n");
							mprintf(WHITE,"      Degree  R             Integral / ps     Lifetime / ps\n");
							ti = -1;
							tf = 0;
							for (z2=o->m_pDACF->m_iFitDegreeMin;z2<=o->m_pDACF->m_iFitDegreeMax;z2++) {
								if (dacfsub->m_pDACF->m_pFitIntegral[z2] >= 0) {
									mprintf("      %d       %10.8f    %-10G        %-10G\n",z2,dacfsub->m_pDACF->m_pCorrCoeff[z2],dacfsub->m_pDACF->m_pFitIntegral[z2],dacfsub->m_pDACF->m_pFitIntegral[z2]*2.0);
									if (dacfsub->m_pDACF->m_pCorrCoeff[z2] > tf) {
										tf = dacfsub->m_pDACF->m_pCorrCoeff[z2];
										ti = z2;
									}
								} else
									mprintf("      %d       %10.8f             -             -\n",z2,dacfsub->m_pDACF->m_pCorrCoeff[z2]);
							}
							mprintf("\n");
						} else
							ti = o->m_pDACF->m_iFitDegreeMin;

						mprintf(WHITE,"    Resulting aggregate lifetime:");
						mprintf(" %-10G ps\n\n",dacfsub->m_pDACF->m_pFitIntegral[ti]*2.0);

						if (dacfsub->m_pDACF->m_pFitIntegral[ti]*2.0 > g_iTrajSteps*g_fTimestepLength/1000.0) {
							mprintf(RED,"    Warning:");
							mprintf(" The lifetime is larger than the total trajectory length! The value\n");
							mprintf("             is a pure guess! Use a (much) longer trajectory to obtain reliable lifetimes.\n\n\n");
						} else if (dacfsub->m_pDACF->m_pFitIntegral[ti]*2.0 > 0.1*g_iTrajSteps*g_fTimestepLength/1000.0) {
							mprintf(RED,"    Warning:");
							mprintf(" The lifetime is larger than 10%% of the total trajectory length. The value\n");
							mprintf("             is not reliable. Use a longer trajectory to obtain reliable lifetimes.\n\n\n");
						}
					}
//					sprintf(buf,"dacf_%s%s.csv",dacfsub->m_sName,multibuf);
					buf.sprintf("dacf_%s%s.csv",dacfsub->m_sName,multibuf);
					mprintf("    Saving DACF File as \"%s\"...\n",(const char*)buf);
					dacfsub->m_pDACF->Write("",buf,"",false);
//					sprintf(buf,"dacf_%s%s.agr",dacfsub->m_sName,multibuf);
					buf.sprintf("dacf_%s%s.agr",dacfsub->m_sName,multibuf);
					mprintf("    Saving DACF Agr File as \"%s\"...\n",(const char*)buf);
					dacfsub->m_pDACF->WriteAgr("",buf,"","",false);

				}

				if (g_bDLDF)
				{
					mprintf(WHITE,"\nDimer Lifetime Distribution Function\n");
					mprintf("    %.0f Bin entries have been made (%.0f skipped, %.0f total).\n",dacfsub->m_pDLDF->m_fBinEntries,dacfsub->m_pDLDF->m_fSkipEntries,dacfsub->m_pDLDF->m_fBinEntries+dacfsub->m_pDLDF->m_fSkipEntries);
//					sprintf(buf,"dldf_%s%s.csv",dacfsub->m_sName,multibuf);
					buf.sprintf("dldf_%s%s.csv",dacfsub->m_sName,multibuf);
					mprintf("    Saving DLDF File as \"%s\"...\n",(const char*)buf);
					dacfsub->m_pDLDF->Write("",buf,"",false);
//					sprintf(buf,"dldf_%s%s.agr",dacfsub->m_sName,multibuf);
					buf.sprintf("dldf_%s%s.agr",dacfsub->m_sName,multibuf);
					mprintf("    Saving DLDF Agr File as \"%s\"...\n",(const char*)buf);
					dacfsub->m_pDLDF->WriteAgr("",buf,"","",false);
				}

				if (g_bDLDisp)
				{
					mprintf(WHITE,"\nDimer Lifetime Displacement Distribution Function\n");
					mprintf("    %.0f Bin entries have been made (%.0f skipped, %.0f total).\n",dacfsub->m_pDLDisp->m_fBinEntries,dacfsub->m_pDLDisp->m_fSkipEntries,dacfsub->m_pDLDisp->m_fBinEntries+dacfsub->m_pDLDisp->m_fSkipEntries);
					mprintf("    Normalizing bin integral to 1000000...\n");
					dacfsub->m_pDLDisp->NormalizeBinIntegral(1000000.0);
//					sprintf(buf,"dldisp_%s%s",dacfsub->m_sName,multibuf);
					buf.sprintf("dldisp_%s%s",dacfsub->m_sName,multibuf);
					mprintf("    Saving DLDisp triples as \"%s_triples.csv\"...\n",(const char*)buf);
					dacfsub->m_pDLDisp->Write("",buf,"_triples.csv");
					mprintf("    Saving DLDisp matrix as \"%s_matrix.csv\"...\n",(const char*)buf);
					dacfsub->m_pDLDisp->WriteCSV("",buf,"_matrix.csv");
					mprintf("    Saving DLDisp Mathematica Notebook \"%s.nb\"...\n",(const char*)buf);
					dacfsub->m_pDLDisp->WriteMathematicaNb("",buf,".nb",false);
					mprintf("    Saving DLDisp Gnuplot Input \"%s.gp\"...\n",(const char*)buf);
					dacfsub->m_pDLDisp->WriteGnuplotInput("",buf,"",false);
				}

				if (g_bPairMSD)
				{
					mprintf(WHITE,"\nPair Mean Square Displacement\n");
					mprintf("    %.0f Bin entries have been made.\n",dacfsub->m_pPairMSD->m_fBinEntries);
					mprintf("    Calculating average values...\n");
					dacfsub->m_pPairMSD->BuildAverage();
//					sprintf(buf,"pairmsd_%s%s",dacfsub->m_sName,multibuf);
					buf.sprintf("pairmsd_%s%s",dacfsub->m_sName,multibuf);
					mprintf("    Saving Pair MSD as \"%s.csv\"...\n",(const char*)buf);
					dacfsub->m_pPairMSD->Write("",buf,".csv");
				}
			}
			if (o->m_pDACF->m_bDACFGrid && o->m_pDACF->m_bFitDACF)
			{
				mprintf(WHITE,"\n*** Condition Grid fitting overview:\n\n");
				mprintf("      Degree    R(min)        R(avg)        R(max)\n");
				o->m_pDACF->CalcGridFitParms();
				for (z2=o->m_pDACF->m_iFitDegreeMin;z2<=o->m_pDACF->m_iFitDegreeMax;z2++)
					mprintf("      %d         %10.8f    %10.8f    %10.8f\n",z2,o->m_pDACF->m_pFitRMin[z2],o->m_pDACF->m_pFitRAvg[z2],o->m_pDACF->m_pFitRMax[z2]);
				mprintf("\n");

				for (z2=o->m_pDACF->m_iFitDegreeMin;z2<=o->m_pDACF->m_iFitDegreeMax;z2++)
				{
					if ((o->m_pDACF->m_iGridMode == 3) || (o->m_pDACF->m_iGridMode == 5))
					{
						if (o->m_pDACF->m_bGridCon)
						{
							try { temp2df = new C2DF(); } catch(...) { temp2df = NULL; }
							if (temp2df == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
							
							o->m_pDACF->CreateGridFit2DF(temp2df,z2,false);
//							sprintf(buf,"dacf_gridint_%s%s_%dexp",o->m_pDACF->m_sName,multibuf,z2);
							buf.sprintf("dacf_gridint_%s%s_%dexp",o->m_pDACF->m_sName,multibuf,z2);
							mprintf("    Saving Grid Integral triples as \"%s_triples.csv\"...\n",(const char*)buf);
							temp2df->Write("",buf,"_triples.csv");
							mprintf("    Saving Grid Integral matrix as \"%s_matrix.csv\"...\n",(const char*)buf);
							temp2df->WriteCSV("",buf,"_matrix.csv");
							mprintf("    Saving Grid Integral Mathematica Notebook \"%s.nb\"...\n",(const char*)buf);
							temp2df->WriteMathematicaNb("",buf,".nb",false);
							mprintf("    Saving Grid Integral Gnuplot Input \"%s.gp\"...\n",(const char*)buf);
							temp2df->WriteGnuplotInput("",buf,"",false);
							delete temp2df;
						}
						if (o->m_pDACF->m_bGridInt)
						{
							try { temp2df = new C2DF(); } catch(...) { temp2df = NULL; }
							if (temp2df == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
							
							o->m_pDACF->CreateGridFit2DF(temp2df,z2,true);
//							sprintf(buf,"dacf_gridint_%s%s_%dexp_int%.2f",o->m_pDACF->m_sName,multibuf,z2,o->m_pDACF->m_fGridIntGap);
							buf.sprintf("dacf_gridint_%s%s_%dexp_int%.2f",o->m_pDACF->m_sName,multibuf,z2,o->m_pDACF->m_fGridIntGap);
							mprintf("    Saving Grid Integral triples as \"%s_triples.csv\"...\n",(const char*)buf);
							temp2df->Write("",buf,"_triples.csv");
							mprintf("    Saving Grid Integral matrix as \"%s_matrix.csv\"...\n",(const char*)buf);
							temp2df->WriteCSV("",buf,"_matrix.csv");
							mprintf("    Saving Grid Integral Mathematica Notebook \"%s.nb\"...\n",(const char*)buf);
							temp2df->WriteMathematicaNb("",buf,".nb",false);
							mprintf("    Saving Grid Integral Gnuplot Input \"%s.gp\"...\n",(const char*)buf);
							temp2df->WriteGnuplotInput("",buf,"",false);
							delete temp2df;
						}
					} else
					{
						if (o->m_pDACF->m_bGridCon)
						{
							try { tdf = new CDF(); } catch(...) { tdf = NULL; }
							if (tdf == NULL) NewException((double)sizeof(CDF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
							
//							o->m_pDACF->CreateGridFitDF(tdf,z2,false);
							eprintf("Internal abort call.\n");
							delete tdf;
							abort();
						}
					}
				}
			}
		} // END IF AGGREGATION

		if (g_bDens)
		{
			mprintf(WHITE,"\n* Density Distribution Function\n");
			mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pDensityDF->m_pDensDF->m_fBinEntries,o->m_pDensityDF->m_pDensDF->m_fSkipEntries,ZeroDivide(o->m_pDensityDF->m_pDensDF->m_fSkipEntries,o->m_pDensityDF->m_pDensDF->m_fBinEntries+o->m_pDensityDF->m_pDensDF->m_fSkipEntries)*100.0);

			mprintf("    Correcting radial distribution...\n");
			o->m_pDensityDF->m_pDensDF->CorrectRadialDist();

			mprintf("    Scaling values to match final density...\n");

			if (o->m_pDensityDF->m_bDensityMass)
			{
				o->m_pDensityDF->m_pDensDF->Integrate(true,1.0 / g_iSteps / (double)((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize());
				o->m_pDensityDF->m_pDensDF->MultiplyBin(1.66054e6 / (4.0/3.0*Pi) / g_iSteps / ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize());
			} else
			{
				o->m_pDensityDF->m_pDensDF->Integrate(true,1.0 / g_iSteps / (double)((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize());
				o->m_pDensityDF->m_pDensDF->MultiplyBin(1.0e9 / (4.0/3.0*Pi) / g_iSteps / ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize());
			}

	/*		if (g_bDoubleBox)
			{
				o->m_pDensityDF->m_pDensDF->MultiplyBin(g_iDoubleBoxFactor);
				o->m_pDensityDF->m_pDensDF->MultiplyIntegral(g_iDoubleBoxFactor);
			}*/

//			sprintf(buf,"density_df_%s%s.csv",o->m_pDensityDF->m_sName,multibuf);
			buf.sprintf("density_df_%s%s.csv",o->m_pDensityDF->m_sName,multibuf);
			mprintf("    Saving Density DF as \"%s\"...\n",(const char*)buf);
			o->m_pDensityDF->m_pDensDF->Write("",buf,"",true);
//			sprintf(buf,"density_df_%s%s.agr",o->m_pDensityDF->m_sName,multibuf);
			buf.sprintf("density_df_%s%s.agr",o->m_pDensityDF->m_sName,multibuf);
			mprintf("    Saving Density DF AGR file as \"%s\"...\n",(const char*)buf);
			o->m_pDensityDF->m_pDensDF->WriteAgr("",buf,"",o->m_pDensityDF->m_sName,true);
			if (o->m_pDensityDF->m_iHistogramRes != 0)
			{
				mprintf("    Calculating Histogram...\n");
				o->m_pDensityDF->m_pDensDF->CalcHistogram();
//				sprintf(buf,"his_density_df_%s%s.csv",o->m_pDensityDF->m_sName,multibuf);
				buf.sprintf("his_density_df_%s%s.csv",o->m_pDensityDF->m_sName,multibuf);
				mprintf("    Saving Histogram as \"%s\"...\n",(const char*)buf);
				o->m_pDensityDF->m_pDensDF->WriteHistogram("",buf,"");
			}
		} // END IF g_bDens

		if (g_bRDyn)
		{
			mprintf(WHITE,"\n* Vector Reorientation Dynamics\n");
			o->m_pRDyn->Finish(multibuf);
		}

		if (g_bIRSpec)
		{
			mprintf(WHITE,"\n* IR Spectrum\n");
			o->m_pIRSpec->Finish(multibuf);
		}

		if (g_bVACF)
		{
			mprintf(WHITE,"\n* Velocity Autocorrelation Function\n");
			if (o->m_pVACF->m_bDerivative)
			{
				mprintf("    Deriving velocities...\n");
				mprintf(WHITE,"      [");
				tfs = (((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pVACF->m_iShowAtomGes)/60.0;
				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pVACF->m_iShowAtomGes;z2++)
				{
					if (fmod(z2,tfs) < 1)
						mprintf(WHITE,"#");
					ptfa = (CxDoubleArray*)o->m_pVACF->m_oaCache[z2];

					for (z3=0;z3<(int)g_iSteps*3-3;z3++)
						(*ptfa)[z3] = (*ptfa)[z3+3] - (*ptfa)[z3];
				}
				mprintf(WHITE,"]\n");
			}
			if (g_bVACFCacheMode)
			{
				mprintf("    Autocorrelating cached vectors...\n");
				mprintf(WHITE,"      [");
				tfs = (((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pVACF->m_iShowAtomGes)/60.0;

				try { ptfa2 = new CxDoubleArray("main():ptfa2"); } catch(...) { ptfa2 = NULL; }
				if (ptfa2 == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				ptfa2->SetSize(g_iSteps);

				try { ptfa3 = new CxDoubleArray("main():ptfa3"); } catch(...) { ptfa3 = NULL; }
				if (ptfa3 == NULL) NewException((double)sizeof(CxDoubleArray),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				ptfa3->SetSize(g_iSteps);

				try { ac = new CAutoCorrelation(); } catch(...) { ac = NULL; }
				if (ac == NULL) NewException((double)sizeof(CAutoCorrelation),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				ac->Init(g_iSteps,o->m_pVACF->m_iSize,g_bACFFFT);
				ti = 0;
				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();z2++)
				{
					for (z3=0;z3<o->m_pVACF->m_oAtoms.m_baAtomType.GetSize();z3++)
					{
						if (o->m_pVACF->m_bMassWeight)
							tf = ((CAtom*)g_oaAtoms[o->m_pVACF->m_oAtoms.m_baRealAtomType[z3]])->m_pElement->m_fMass;
								else tf = 1.0;
						for (z4=0;z4<((CxIntArray*)o->m_pVACF->m_oAtoms.m_oaAtoms[z3])->GetSize();z4++)
						{
							if (fmod(ti,tfs) < 1)
								mprintf(WHITE,"#");
							ptfa = (CxDoubleArray*)o->m_pVACF->m_oaCache[ti];

							/* X */
							for (z5=0;z5<(int)g_iSteps;z5++)
								(*ptfa2)[z5] = (*ptfa)[z5*3];
							ac->AutoCorrelate(ptfa2,ptfa3);
			//				mprintf("Lokal Atom %d X: %f + %f = %f.\n",ti+1,o->m_pVACF->m_pData[0],(*ptfa3)[0]*tf,o->m_pVACF->m_pData[0]+(*ptfa3)[0]*tf);
							for (z5=0;z5<(int)o->m_pVACF->m_iSize;z5++)
								o->m_pVACF->m_pData[z5] += (*ptfa3)[z5]*tf;

							/* Y */
							for (z5=0;z5<(int)g_iSteps;z5++)
								(*ptfa2)[z5] = (*ptfa)[z5*3+1];
							ac->AutoCorrelate(ptfa2,ptfa3);
			//				mprintf("Lokal Atom %d Y: %f + %f = %f.\n",ti+1,o->m_pVACF->m_pData[0],(*ptfa3)[0]*tf,o->m_pVACF->m_pData[0]+(*ptfa3)[0]*tf);
							for (z5=0;z5<(int)o->m_pVACF->m_iSize;z5++)
								o->m_pVACF->m_pData[z5] += (*ptfa3)[z5]*tf;

							/* Z */
							for (z5=0;z5<(int)g_iSteps;z5++)
								(*ptfa2)[z5] = (*ptfa)[z5*3+2];
							ac->AutoCorrelate(ptfa2,ptfa3);
			//				mprintf("Lokal Atom %d Z: %f + %f = %f.\n",ti+1,o->m_pVACF->m_pData[0],(*ptfa3)[0]*tf,o->m_pVACF->m_pData[0]+(*ptfa3)[0]*tf);
							for (z5=0;z5<(int)o->m_pVACF->m_iSize;z5++)
								o->m_pVACF->m_pData[z5] += (*ptfa3)[z5]*tf;

							ti++;
						}
					}
				}
				delete ac;
				delete ptfa2;
				delete ptfa3;

			/*	} else
				{
					for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pVACF->m_iShowAtomGes;z2++)
					{
						if (fmod(z2,tfs) < 1)
							mprintf(WHITE,"#");
						ptfa = (CxFloatArray*)o->m_pVACF->m_oaCache[z2];
						for (z3=0;z3<o->m_pVACF->m_iSize;z3+=o->m_pVACF->m_iStride) // Das ist das Tau
						{
							tf = 0;
							for (z4=0;z4<(int)g_iSteps-z3;z4++) // Das ist der Startpunkt
								tf += (*ptfa)[z4*3]*(*ptfa)[(z4+z3)*3] + (*ptfa)[z4*3+1]*(*ptfa)[(z4+z3)*3+1] + (*ptfa)[z4*3+2]*(*ptfa)[(z4+z3)*3+2];
							o->m_pVACF->m_pData[z3/o->m_pVACF->m_iStride] += tf/(g_iSteps-z3);
						}
					}
				}*/
				mprintf(WHITE,"]\n");
//				mprintf("      %d atoms, %d time steps and %d correlation depths processed.\n",g_pGlobalVACF->m_oaCache.GetSize(),g_iSteps,g_pGlobalVACF->m_iSize);
				mprintf("      %d atoms, %lu time steps and %d correlation depths processed.\n",ti,g_iSteps,o->m_pVACF->m_iSize);
				o->m_pVACF->MultiplyCached(1.0/g_iGesAtomCount);
			} else o->m_pVACF->Multiply(1.0/g_iGesAtomCount);
/*			sprintf(buf,"acf_%s%s.csv",o->m_pVACF->m_sName,multibuf);
			mprintf("    Saving ACF as %s ...\n",buf);
			o->m_pVACF->WriteACF("",buf,"");*/

//			sprintf(buf,"acf_%s%s.csv",o->m_pVACF->m_sName,multibuf);
			buf.sprintf("acf_%s%s.csv",o->m_pVACF->m_sName,multibuf);
			mprintf("    Saving ACF as %s ...\n",(const char*)buf);
			o->m_pVACF->WriteACF("",buf,"");

			if (o->m_pVACF->m_iMirror != 0)
			{
				mprintf("    Mirroring ACF...\n");
				o->m_pVACF->Mirror(o->m_pVACF->m_iMirror);
//				sprintf(buf,"acf_%s%s.m.csv",o->m_pVACF->m_sName,multibuf);
				buf.sprintf("acf_%s%s.m.csv",o->m_pVACF->m_sName,multibuf);
				mprintf("    Saving mirrored ACF as %s ...\n",(const char*)buf);
				o->m_pVACF->WriteACF("",buf,"");
			}

			if (o->m_pVACF->m_bWindowFunction)
			{
				mprintf("    Applying window function to ACF...\n");
				o->m_pVACF->Window();
//				sprintf(buf,"acf_%s%s.w.csv",o->m_pVACF->m_sName,multibuf);
				buf.sprintf("acf_%s%s.w.csv",o->m_pVACF->m_sName,multibuf);
				mprintf("    Saving windowed ACF as %s ...\n",(const char*)buf);
				o->m_pVACF->WriteACF("",buf,"");
			}

			if (o->m_pVACF->m_bSpectrum)
			{
				mprintf("    Performing fourier transformation...\n");

				try { g_pFFT = new CFFT(); } catch(...) { g_pFFT = NULL; }
				if (g_pFFT == NULL) NewException((double)sizeof(CFFT),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
/*				if (o->m_pVACF->m_bMirror)
					g_pFFT->PrepareFFT_C2C((o->m_pVACF->m_iSize+o->m_pVACF->m_iZeroPadding)*2);
						else */g_pFFT->PrepareFFT_C2C(o->m_pVACF->m_iSize+o->m_pVACF->m_iZeroPadding);
				o->m_pVACF->Transform(g_pFFT);
				delete g_pFFT;
				o->m_pVACF->m_pSpectrum->SetMaxRWL(1E7/299.792/g_fTimestepLength/g_iStride);
				if (o->m_pVACF->m_bACF_DB)
				{
					mprintf("    Normalising spectrum to decibel...\n");
					o->m_pVACF->m_pSpectrum->MakeDB();
				}/* else
				{
					mprintf("    Normalising integral of spectrum...\n");
					o->m_pVACF->m_pSpectrum->SetIntegral(1000000.0f);
				}*/
/*				if (o->m_pVACF->m_bWindowFunction)
				{
					sprintf(buf,"power_%s_w%s.csv",o->m_pVACF->m_sName,multibuf);
					mprintf("    Saving spectrum as %s ...\n",buf);
					o->m_pVACF->m_pSpectrum->Write("",buf,"");
				} else
				{*/
//				sprintf(buf,"power_%s%s.csv",o->m_pVACF->m_sName,multibuf);
				buf.sprintf("power_%s%s.csv",o->m_pVACF->m_sName,multibuf);
				mprintf("    Saving spectrum as %s ...\n",(const char*)buf);
				o->m_pVACF->m_pSpectrum->Write("",buf,"");
//				}
			}

		} // End IF VACF

		if (g_bMSD)
		{
			mprintf(WHITE,"\n* Mean Square Displacement\n");
			if ((int)g_iSteps/g_iStride <= o->m_pMSD->m_iResolution) {
				eprintf("\n    Error: Requested MSD resolution of %d, but only %lu trajectory frames processed.\n",o->m_pMSD->m_iResolution,g_iSteps/g_iStride);
				eprintf("           Expect erroneous results. Reduce the MSD resolution or use a longer trajectory.\n\n");
			}
			if (g_bMSDCacheMode)
			{
				mprintf("    Computing square displacement of cached vectors...\n");
//				mprintf(WHITE,"      [");
//				tfs = (((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pMSD->m_iShowAtoms)/60.0;
				for (z2=0;z2<((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pMSD->m_iShowAtoms;z2++)
				{
//					if (fmod(z2,tfs) < 1.0)
//						mprintf(WHITE,"#");
					mprintf("      %4d/%d:  [",z2+1,((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pMSD->m_iShowAtoms);
					tfs = o->m_pMSD->m_iResolution/30.0;
					ptfa = (CxDoubleArray*)o->m_pMSD->m_oaCache[z2];
					for (z3=0;z3<o->m_pMSD->m_iResolution;z3+=o->m_pMSD->m_iStride) // Das ist das Tau
					{
						if (fmod(z3,tfs) < o->m_pMSD->m_iStride)
							mprintf(WHITE,"#");
						tf = 0;

						if (o->m_pMSD->m_bTakeX && o->m_pMSD->m_bTakeY && o->m_pMSD->m_bTakeZ)
						{
							for (z4=0;z4<(int)g_iSteps/g_iStride-z3-1;z4+=o->m_pMSD->m_iStride2) // Das ist der Startpunkt
								tf += pow2((*ptfa)[(z3+z4)*3]-(*ptfa)[z4*3]) + pow2((*ptfa)[(z3+z4)*3+1]-(*ptfa)[z4*3+1]) + pow2((*ptfa)[(z3+z4)*3+2]-(*ptfa)[z4*3+2]);
						} else
						{
							for (z4=0;z4<(int)g_iSteps/g_iStride-z3-1;z4+=o->m_pMSD->m_iStride2) // Das ist der Startpunkt
							{
								if (o->m_pMSD->m_bTakeX)
									tf += pow2((*ptfa)[(z3+z4)*3]-(*ptfa)[z4*3]);
								if (o->m_pMSD->m_bTakeY)
									tf += pow2((*ptfa)[(z3+z4)*3+1]-(*ptfa)[z4*3+1]);
								if (o->m_pMSD->m_bTakeZ)
									tf += pow2((*ptfa)[(z3+z4)*3+2]-(*ptfa)[z4*3+2]);
							}
						}

						o->m_pMSD->m_pMSD->AddToBin_Index(z3/o->m_pMSD->m_iStride,tf/(g_iSteps/g_iStride-z3-1)*o->m_pMSD->m_iStride2);
						o->m_pMSD->m_pMSD->m_fBinEntries += (g_iSteps/g_iStride-z3-1)/o->m_pMSD->m_iStride2;
						if (o->m_pMSD->m_bSplit)
						{
							o->m_pMSD->m_pSplitMSD[z2]->AddToBin_Index(z3/o->m_pMSD->m_iStride,tf/(g_iSteps/g_iStride-z3-1)*o->m_pMSD->m_iStride2);
							o->m_pMSD->m_pSplitMSD[z2]->m_fBinEntries += (g_iSteps/g_iStride-z3-1)/o->m_pMSD->m_iStride2;
						}
					}
					mprintf("]\n");
				}
				for (z3=0;z3<o->m_pMSD->m_iResolution;z3+=o->m_pMSD->m_iStride) // Das ist das Tau
					o->m_pMSD->m_pMSD->m_pBin[z3/o->m_pMSD->m_iStride] /= ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()*o->m_pMSD->m_iShowAtoms;
//				mprintf(WHITE,"]\n");
			} else
			{
				o->m_pMSD->m_pMSD->BuildAverage();
			}
			mprintf("    %.0f bin entries.\n",o->m_pMSD->m_pMSD->m_fBinEntries);
			o->m_pMSD->m_pMSD->CalcDeriv(1.0);
//			sprintf(buf,"msd_%s%s.csv",o->m_pMSD->m_sName,multibuf);
			buf.sprintf("msd_%s%s.csv",o->m_pMSD->m_sName,multibuf);
			mprintf("\n");
			mprintf("    Saving result as %s ...\n",(const char*)buf);
			if (o->m_pMSD->m_bSplit)
			{
				o->m_pMSD->WriteSplit(buf);
			} else o->m_pMSD->m_pMSD->Write("",buf,"");
			mprintf(WHITE,"\n    Performing linear regression on interval %.2f - %.2f ps...\n\n",o->m_pMSD->m_pMSD->m_fMinVal+0.5*(o->m_pMSD->m_pMSD->m_fMaxVal-o->m_pMSD->m_pMSD->m_fMinVal),o->m_pMSD->m_pMSD->m_fMaxVal);
			o->m_pMSD->m_pMSD->LinReg(o->m_pMSD->m_pMSD->m_iResolution/2,o->m_pMSD->m_pMSD->m_iResolution-1,&c0,&c1,&r);
			mprintf("      MSD(t) = %.6f + %.6f * t   (units: [MSD] = pm^2, [t] = ps).\n",c0,c1);
			mprintf("           R = %.6f   (correlation coefficient)\n",r);
			ti = 0;
			if (o->m_pMSD->m_bTakeX)
				ti++;
			if (o->m_pMSD->m_bTakeY)
				ti++;
			if (o->m_pMSD->m_bTakeZ)
				ti++;
			if (ti != 3)
				mprintf("     Please note: Dimensionality of the system is %d here.\n",ti);
			mprintf("\n      Diffusion coefficient D = %.6f pm^2/ps  =  %G m^2/s.\n",c1/2.0/ti,c1/2.0e12/ti);
			mprintf("        (assuming <x^2> = %d * D * t)\n\n",2*ti);

			buf.sprintf("msd_%s%s_fit.csv",o->m_pMSD->m_sName,multibuf);
			mprintf("    Saving linear regression curve as %s ...\n",(const char*)buf);

			tfi = OpenFileWrite( buf, true );

			mfprintf( tfi, "# Tau / ps;  Linear Fit to MSD / pm^2\n" );

			for (z2=0;z2<200;z2++) {
				tf = o->m_pMSD->m_pMSD->m_fMinVal+(0.5+(double)z2/399.0)*(o->m_pMSD->m_pMSD->m_fMaxVal-o->m_pMSD->m_pMSD->m_fMinVal);
				mfprintf( tfi, "%.4f;  %.4f\n", tf, c0 + c1*tf );
			}

			fclose( tfi );

		} // END IF MSD

		if (g_bNbAnalysis)
		{
			mprintf(WHITE,"\n* Neighborhood Analysis\n");
			mprintf("    %.0f bin entries.\n",o->m_pNbAnalysis->m_fBinEntries);
/*			for (z2=0;z2<o->m_pNbAnalysis->m_iMaxNbCount-o->m_pNbAnalysis->m_iMinNbCount;z2++)
				mprintf("    - %d: %.0f bin entries, %.0f skipped.\n",z2+1,((CDF*)o->m_pNbAnalysis->m_oaDF[z2])->m_fBinEntries,((CDF*)o->m_pNbAnalysis->m_oaDF[z2])->m_fSkipEntries);
*/			

			try { pf = new double[o->m_pNbAnalysis->m_iNbCount]; } catch(...) { pf = NULL; }
			if (pf == NULL) NewException((double)o->m_pNbAnalysis->m_iNbCount*sizeof(double),__FILE__,__LINE__,__PRETTY_FUNCTION__);
			
			for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount;z3++)
				pf[z3] = (double)((CDF*)o->m_pNbAnalysis->m_oaDF[z3])->NormBinIntegral();

//			sprintf(buf,"nbh_ncf_%s%s.csv",o->m_pNbAnalysis->m_sName,multibuf);
			buf.sprintf("nbh_ncf_%s%s.csv",o->m_pNbAnalysis->m_sName,multibuf);
			mprintf("    Saving neighbor count function as %s ...\n",(const char*)buf);
			a = OpenFileWrite(buf,true);
			mfprintf(a,"# distance / pm; ");
			for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount+1;z3++)
			{
				mfprintf(a,"%d neighbors",z3);
				if (z3 < o->m_pNbAnalysis->m_iNbCount)
					mfprintf(a,";  ");
			}
			mfprintf(a,"\n");
			for (z2=0;z2<o->m_pNbAnalysis->m_iResolution;z2++)
			{
				mfprintf(a,"%.4f;  ",o->m_pNbAnalysis->m_fMinDist+(z2+0.5)*(o->m_pNbAnalysis->m_fMaxDist-o->m_pNbAnalysis->m_fMinDist)/o->m_pNbAnalysis->m_iResolution);
				for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount+1;z3++)
				{
					mfprintf(a,"%.6f",((CDF*)o->m_pNbAnalysis->m_oaNPF[z3])->m_pBin[z2]/o->m_pNbAnalysis->m_pNPFCount->m_pBin[z2]*100.0);
					if (z3 < o->m_pNbAnalysis->m_iNbCount)
						mfprintf(a,";  ");
				}
				mfprintf(a,"\n");
			}
			fclose(a);

//			sprintf(buf,"nbh_dist_%s%s.csv",o->m_pNbAnalysis->m_sName,multibuf);
			buf.sprintf("nbh_dist_%s%s.csv",o->m_pNbAnalysis->m_sName,multibuf);
			mprintf("    Calculating mean values and standard deviation...\n");
			for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount;z3++)
				((CDF*)o->m_pNbAnalysis->m_oaDF[z3])->CalcMeanSD();
			mprintf("    Saving neighborhood distribution as %s ...\n",(const char*)buf);
			a = OpenFileWrite(buf,true);
			mfprintf(a,"# distance / pm; ");
			for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount;z3++)
			{
				mfprintf(a,"%d. neighbor",z3+1);
				if (z3 < o->m_pNbAnalysis->m_iNbCount-1)
					mfprintf(a,";  ");
			}
			mfprintf(a,"\n");
			for (z2=0;z2<o->m_pNbAnalysis->m_iResolution;z2++)
			{
				mfprintf(a,"%.4f;  ",o->m_pNbAnalysis->m_fMinDist+(z2+0.5)*(o->m_pNbAnalysis->m_fMaxDist-o->m_pNbAnalysis->m_fMinDist)/o->m_pNbAnalysis->m_iResolution);
				for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount;z3++)
				{
					mfprintf(a,"%.6f",((CDF*)o->m_pNbAnalysis->m_oaDF[z3])->m_pBin[z2]);
					if (z3 < o->m_pNbAnalysis->m_iNbCount-1)
						mfprintf(a,";  ");
				}
				mfprintf(a,"\n");
			}
			fclose(a);

//			sprintf(buf,"nbh_minmaxavgsd_%s%s.csv",o->m_pNbAnalysis->m_sName,multibuf);
			buf.sprintf("nbh_minmaxavgsd_%s%s.csv",o->m_pNbAnalysis->m_sName,multibuf);
			mprintf("    Saving neighborhood min/max/avg/sd as %s ...\n",(const char*)buf);
			a = OpenFileWrite(buf,true);
			mfprintf(a,"# n-th neighbor;  min. dist / pm; max. dist / pm; avg. dist / pm; standard deviation / pm\n");
			for (z2=0;z2<o->m_pNbAnalysis->m_iNbCount;z2++)
			{
				mfprintf(a,"%d;  %.4f;  %.4f;  %.4f;  %.4f",z2+1,o->m_pNbAnalysis->m_pDistMin[z2],o->m_pNbAnalysis->m_pDistMax[z2],((CDF*)o->m_pNbAnalysis->m_oaDF[z2])->m_fMean,((CDF*)o->m_pNbAnalysis->m_oaDF[z3])->m_fSD);
				mfprintf(a,"\n");
			}
			fclose(a);

//			sprintf(buf,"nbh_rdf_decomp_%s%s.csv",o->m_pNbAnalysis->m_sName,multibuf);
			buf.sprintf("nbh_rdf_decomp_%s%s.csv",o->m_pNbAnalysis->m_sName,multibuf);
			mprintf("    Saving neighborhood RDF decomposition as %s ...\n",(const char*)buf);
			a = OpenFileWrite(buf,true);
			mfprintf(a,"# distance / pm; ");
			for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount;z3++)
			{
				mfprintf(a,"%d. neighbor",z3+1);
				if (z3 < o->m_pNbAnalysis->m_iNbCount-1)
					mfprintf(a,";  ");
			}
			mfprintf(a,"\n");
			for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount;z3++)
			{
				((CDF*)o->m_pNbAnalysis->m_oaDF[z3])->MultiplyBin(pf[z3]);
				((CDF*)o->m_pNbAnalysis->m_oaDF[z3])->CorrectRadialDist();
			}
			for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount;z3++)
			{
				((CDF*)o->m_pNbAnalysis->m_oaDF[z3])->MultiplyBin(g_fBoxX*g_fBoxY*g_fBoxZ / (4.0/3.0*Pi) / g_iSteps / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize());
				if (g_bDoubleBox)
					((CDF*)o->m_pNbAnalysis->m_oaDF[z3])->MultiplyBin(g_iDoubleBoxFactor);
			}
			for (z2=0;z2<o->m_pNbAnalysis->m_iResolution;z2++)
			{
				mfprintf(a,"%.4f; ",o->m_pNbAnalysis->m_fMinDist+(z2+0.5)*(o->m_pNbAnalysis->m_fMaxDist-o->m_pNbAnalysis->m_fMinDist)/o->m_pNbAnalysis->m_iResolution);
				for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount;z3++)
				{
					mfprintf(a,"%.6f",((CDF*)o->m_pNbAnalysis->m_oaDF[z3])->m_pBin[z2]);
					if (z3 < o->m_pNbAnalysis->m_iNbCount-1)
						mfprintf(a,";  ");
				}
				mfprintf(a,"\n");
			}
			fclose(a);

//			sprintf(buf,"nbh_rdf_decomp_cumulative_%s%s.csv",o->m_pNbAnalysis->m_sName,multibuf);
			buf.sprintf("nbh_rdf_decomp_cumulative_%s%s.csv",o->m_pNbAnalysis->m_sName,multibuf);
			mprintf("    Saving cumulative neighborhood RDF decomposition as %s ...\n",(const char*)buf);
			a = OpenFileWrite(buf,true);
			mfprintf(a,"# distance / pm; ");
			for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount;z3++)
			{
				mfprintf(a,"%d. neighbor",z3+1);
				if (z3 < o->m_pNbAnalysis->m_iNbCount-1)
					mfprintf(a,";  ");
			}
			mfprintf(a,"\n");
			for (z2=0;z2<o->m_pNbAnalysis->m_iResolution;z2++)
				for (z3=1;z3<o->m_pNbAnalysis->m_iNbCount;z3++)
					((CDF*)o->m_pNbAnalysis->m_oaDF[z3])->m_pBin[z2] += ((CDF*)o->m_pNbAnalysis->m_oaDF[z3-1])->m_pBin[z2];
			for (z2=0;z2<o->m_pNbAnalysis->m_iResolution;z2++)
			{
				mfprintf(a,"%.4f; ",o->m_pNbAnalysis->m_fMinDist+(z2+0.5)*(o->m_pNbAnalysis->m_fMaxDist-o->m_pNbAnalysis->m_fMinDist)/o->m_pNbAnalysis->m_iResolution);
				for (z3=0;z3<o->m_pNbAnalysis->m_iNbCount;z3++)
				{
					mfprintf(a,"%.6f",((CDF*)o->m_pNbAnalysis->m_oaDF[z3])->m_pBin[z2]);
					if (z3 < o->m_pNbAnalysis->m_iNbCount-1)
						mfprintf(a,";  ");
				}
				mfprintf(a,"\n");
			}
			fclose(a);
			delete[] pf;
		}

		if (g_bDipDF)
		{
			for (zr=0;zr<g_iCDFChannels;zr++)
			{
				if (o->m_pDipDF[zr] == NULL)
					continue;
				mprintf(WHITE,"\n* Dipole Distribution Function\n");
				mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pDipDF[zr]->m_pDipoleDF->m_fBinEntries,o->m_pDipDF[zr]->m_pDipoleDF->m_fSkipEntries,ZeroDivide(o->m_pDipDF[zr]->m_pDipoleDF->m_fSkipEntries,o->m_pDipDF[zr]->m_pDipoleDF->m_fSkipEntries+o->m_pDipDF[zr]->m_pDipoleDF->m_fBinEntries)*100.0);
				o->m_pDipDF[zr]->m_pDipoleDF->CalcMinMax();
				mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pDipDF[zr]->m_pDipoleDF->m_fMaxEntry,o->m_pDipDF[zr]->m_pDipoleDF->m_fEps);
				if (o->m_pDipDF[zr]->m_pDipoleDF->m_fEps > 1.0E-4) {
					eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
					mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
				}
				o->m_pDipDF[zr]->m_pDipoleDF->CalcMeanSD();
				mprintf("    Mean value: %10G Debye    Standard deviation: %10G Debye\n",o->m_pDipDF[zr]->m_pDipoleDF->m_fMean,o->m_pDipDF[zr]->m_pDipoleDF->m_fSD);
				mprintf("    Min. value: %10G Debye    Max. value:         %10G Debye\n",o->m_pDipDF[zr]->m_pDipoleDF->m_fMinInput,o->m_pDipDF[zr]->m_pDipoleDF->m_fMaxInput);
				o->m_pDipDF[zr]->m_pDipoleDF->NormBinIntegral();
//				sprintf(buf,"dipole_%s%s.csv",o->m_pDipDF[zr]->m_sName,multibuf);
				buf.sprintf("dipole_%s%s.csv",o->m_pDipDF[zr]->m_sName,multibuf);
				mprintf("    Saving dipole distribution as %s ...\n",(const char*)buf);
				o->m_pDipDF[zr]->m_pDipoleDF->Write("",buf,"",false);
//				sprintf(buf,"dipole_%s%s.agr",o->m_pDipDF[zr]->m_sName,multibuf);
				buf.sprintf("dipole_%s%s.agr",o->m_pDipDF[zr]->m_sName,multibuf);
				mprintf("    Saving dipole distribution AGR file as \"%s\"...\n",(const char*)buf);
				o->m_pDipDF[zr]->m_pDipoleDF->WriteAgr("",buf,"",o->m_pDipDF[zr]->m_sName,false);
				if (o->m_pDipDF[zr]->m_iHistogramRes != 0)
				{
					mprintf("    Calculating Histogram...\n");
					o->m_pDipDF[zr]->m_pDipoleDF->CalcHistogram();
//					sprintf(buf,"his_dipdf_%s%s.csv",o->m_pDipDF[zr]->m_sName,multibuf);
					buf.sprintf("his_dipdf_%s%s.csv",o->m_pDipDF[zr]->m_sName,multibuf);
					mprintf("    Saving Histogram as \"%s\"...\n",(const char*)buf);
					o->m_pDipDF[zr]->m_pDipoleDF->WriteHistogram("",buf,"");
				}
				if (o->m_bTimeDev)
				{
					if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
					{
						fclose(o->m_pDipDF[zr]->m_fDipole[0]);
						mprintf("    Saving temporal development as dipole_timedev_%s%s.csv\n",o->m_pDipDF[zr]->m_sName,multibuf);
					} else
					{
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
							fclose(o->m_pDipDF[zr]->m_fDipole[z2]);
							mprintf("    Saving temporal development as dipole_timedev_%s_ref%d%s.csv\n",o->m_pDipDF[zr]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
						}
					}
					if (o->m_bCombinedPlot)
					{
//						sprintf(buf,"combined_%s%s.agr",o->m_pDipDF[zr]->m_sName,multibuf);
						buf.sprintf("combined_%s%s.agr",o->m_pDipDF[zr]->m_sName,multibuf);
						mprintf("    Saving combined plot as \"%s\"...\n",(const char*)buf);
						o->m_pDipDF[zr]->m_pDipoleDF->CreateCombinedPlot(false);
						o->m_pDipDF[zr]->m_pDipoleDF->m_pCombinedPlot->WriteAgr(buf,false);
					}
				} // END IF TIMEDEV

				if (o->m_bObsCertain && o->m_bDecompDist)
				{
//					sprintf(buf,"dipdf_decomp_%s%s.csv",o->m_pDipDF[zr]->m_sName,multibuf);
					buf.sprintf("dipdf_decomp_%s%s.csv",o->m_pDipDF[zr]->m_sName,multibuf);
					mprintf("    Saving DipDF decomposition as \"%s\"...\n",(const char*)buf);
					o->m_pDipDF[zr]->m_pDipoleDF->WriteMulti("",buf,"");
//					sprintf(buf,"dipdf_decomp_%s%s_cumulative.csv",o->m_pDipDF[zr]->m_sName,multibuf);
					buf.sprintf("dipdf_decomp_%s%s_cumulative.csv",o->m_pDipDF[zr]->m_sName,multibuf);
					mprintf("    Saving cumulative DipoleDF decomposition as \"%s\"...\n",(const char*)buf);
					o->m_pDipDF[zr]->m_pDipoleDF->WriteMulti_Cumulative("",buf,"");
				}

				if (o->m_bTimeDiff)
					o->WriteTimeDiff(o->m_pDipDF[zr]->m_pDipoleDF,"DipDF","dipdf",o->m_pDipDF[zr]->m_sName,multibuf,false);
			}
		} // END IF DIPOLE

		if (g_bVDF)
		{
			for (zr=0;zr<g_iCDFChannels;zr++)
			{
				if (o->m_pVDF[zr] == NULL)
					continue;
				mprintf(WHITE,"\n* Velocity Distribution Function\n");
				mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pVDF[zr]->m_pVDF->m_fBinEntries,o->m_pVDF[zr]->m_pVDF->m_fSkipEntries,ZeroDivide(o->m_pVDF[zr]->m_pVDF->m_fSkipEntries,o->m_pVDF[zr]->m_pVDF->m_fBinEntries+o->m_pVDF[zr]->m_pVDF->m_fSkipEntries)*100.0);
				o->m_pVDF[zr]->m_pVDF->CalcMinMax();
				mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pVDF[zr]->m_pVDF->m_fMaxEntry,o->m_pVDF[zr]->m_pVDF->m_fEps);
				if (o->m_pVDF[zr]->m_pVDF->m_fEps > 1.0E-4) {
					eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
					mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
				}
				o->m_pVDF[zr]->m_pVDF->CalcMeanSD();
				mprintf("    Mean value: %10G pm/ps    Standard deviation: %10G pm/ps\n",o->m_pVDF[zr]->m_pVDF->m_fMean,o->m_pVDF[zr]->m_pVDF->m_fSD);
				mprintf("    Min. value: %10G pm/ps    Max. value:         %10G pm/ps\n",o->m_pVDF[zr]->m_pVDF->m_fMinInput,o->m_pVDF[zr]->m_pVDF->m_fMaxInput);
				o->m_pVDF[zr]->m_pVDF->NormBinIntegral();
//				sprintf(buf,"vdf_%s%s.csv",o->m_pVDF[zr]->m_sName,multibuf);
				buf.sprintf("vdf_%s%s.csv",o->m_pVDF[zr]->m_sName,multibuf);
				mprintf("    Saving velocity distribution as %s ...\n",(const char*)buf);
				o->m_pVDF[zr]->m_pVDF->Write("",buf,"",false);
//				sprintf(buf,"vdf_%s%s.agr",o->m_pVDF[zr]->m_sName,multibuf);
				buf.sprintf("vdf_%s%s.agr",o->m_pVDF[zr]->m_sName,multibuf);
				mprintf("    Saving velocity distribution AGR file as \"%s\"...\n",(const char*)buf);
				o->m_pVDF[zr]->m_pVDF->WriteAgr("",buf,"",o->m_pVDF[zr]->m_sName,false);
				if (o->m_pVDF[zr]->m_iHistogramRes != 0) {
					mprintf("    Calculating Histogram...\n");
					o->m_pVDF[zr]->m_pVDF->CalcHistogram();
//					sprintf(buf,"his_vdf_%s%s.agr",o->m_pVDF[zr]->m_sName,multibuf);
					buf.sprintf("his_vdf_%s%s.csv",o->m_pVDF[zr]->m_sName,multibuf);
					mprintf("    Saving Histogram as \"%s\"...\n",(const char*)buf);
					o->m_pVDF[zr]->m_pVDF->WriteHistogram("",buf,"");
				}
				if (o->m_bTimeDev)
				{
					if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
					{
						fclose(o->m_pVDF[zr]->m_fSpeed[0]);
						mprintf("    Saving temporal development as vdf_timedev_%s%s.csv\n",o->m_pVDF[zr]->m_sName,multibuf);
					} else
					{
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
							fclose(o->m_pVDF[zr]->m_fSpeed[z2]);
							mprintf("    Saving temporal development as vdf_timedev_%s_ref%d%s.csv\n",o->m_pVDF[zr]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
						}
					}
					if (o->m_bCombinedPlot)
					{
//						sprintf(buf,"combined_%s%s.csv",o->m_pVDF[zr]->m_sName,multibuf);
						buf.sprintf("combined_%s%s.csv",o->m_pVDF[zr]->m_sName,multibuf);
						mprintf("    Saving combined plot as \"%s\"...\n",(const char*)buf);
						o->m_pVDF[zr]->m_pVDF->CreateCombinedPlot(false);
						o->m_pVDF[zr]->m_pVDF->m_pCombinedPlot->WriteAgr(buf,false);
					}
				} // END IF TIMEDEV

				if (o->m_bObsCertain && o->m_bDecompDist)
				{
//					sprintf(buf,"vdf_decomp_%s%s.csv",o->m_pVDF[zr]->m_sName,multibuf);
					buf.sprintf("vdf_decomp_%s%s.csv",o->m_pVDF[zr]->m_sName,multibuf);
					mprintf("    Saving VDF decomposition as \"%s\"...\n",(const char*)buf);
					o->m_pVDF[zr]->m_pVDF->WriteMulti("",buf,"");
//					sprintf(buf,"vdf_decomp_%s%s_cumulative.csv",o->m_pVDF[zr]->m_sName,multibuf);
					buf.sprintf("vdf_decomp_%s%s_cumulative.csv",o->m_pVDF[zr]->m_sName,multibuf);
					mprintf("    Saving cumulative VDF decomposition as \"%s\"...\n",(const char*)buf);
					o->m_pVDF[zr]->m_pVDF->WriteMulti_Cumulative("",buf,"");
				}

				if (o->m_bTimeDiff)
					o->WriteTimeDiff(o->m_pVDF[zr]->m_pVDF,"VDF","vdf",o->m_pVDF[zr]->m_sName,multibuf,false);

				if (o->m_pVDF[zr]->m_bSplitCart) {
					upn = "?";
					pn = "?";
					for (z2=0;z2<6;z2++) {
						switch(z2) {
							case 0: pn = "x"; upn = "X"; break;
							case 1: pn = "y"; upn = "Y"; break;
							case 2: pn = "z"; upn = "Z"; break;
							case 3: pn = "xy"; upn = "XY"; break;
							case 4: pn = "xz"; upn = "XZ"; break;
							case 5: pn = "yz"; upn = "YZ"; break;
						}
						mprintf(WHITE,"    * Velocity Distribution %s Projection\n",upn);
						mprintf("        %.0f bin entries, %.0f out of bin range (%.2f percent).\n",
							o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fBinEntries,
							o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fSkipEntries,
							ZeroDivide(
								o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fSkipEntries,
								o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fBinEntries+o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fSkipEntries
							)*100.0
						);
						o->m_pVDF[zr]->m_pVDFSplit[z2]->CalcMinMax();
						mprintf("        Max. bin entry: %.6E  -->  EPS: %.6E\n",
							o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fMaxEntry,
							o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fEps
						);
						if (o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fEps > 1.0E-4) {
							eprintf("\n        Warning: Very large bin entries - probably loss of accuracy occurred.\n");
							mprintf("        Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
						}
						o->m_pVDF[zr]->m_pVDFSplit[z2]->CalcMeanSD();
						mprintf("        Mean value: %10G pm/ps    Standard deviation: %10G pm/ps\n",
							o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fMean,
							o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fSD
						);
						mprintf("        Min. value: %10G pm/ps    Max. value:         %10G pm/ps\n",
							o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fMinInput,
							o->m_pVDF[zr]->m_pVDFSplit[z2]->m_fMaxInput
						);
						o->m_pVDF[zr]->m_pVDFSplit[z2]->NormBinIntegral();
						buf.sprintf("vdf_%s%s_%s.csv",o->m_pVDF[zr]->m_sName,multibuf,pn);
						mprintf("        Saving velocity distribution as \"%s\"...\n",(const char*)buf);
						o->m_pVDF[zr]->m_pVDFSplit[z2]->Write("",buf,"",false);
						buf.sprintf("vdf_%s%s_%s.agr",o->m_pVDF[zr]->m_sName,multibuf,pn);
						mprintf("        Saving velocity distribution AGR file as \"%s\"...\n",(const char*)buf);
						o->m_pVDF[zr]->m_pVDFSplit[z2]->WriteAgr("",buf,"",o->m_pVDF[zr]->m_sName,false);
						if (o->m_pVDF[zr]->m_iHistogramRes != 0) {
							mprintf("        Calculating Histogram...\n");
							o->m_pVDF[zr]->m_pVDFSplit[z2]->CalcHistogram();
							buf.sprintf("his_vdf_%s%s_%s.csv",o->m_pVDF[zr]->m_sName,multibuf,pn);
							mprintf("        Saving Histogram as \"%s\"...\n",(const char*)buf);
							o->m_pVDF[zr]->m_pVDFSplit[z2]->WriteHistogram("",buf,"");
						}
					}
				}
			}
		} // END IF VDF

		if (g_bPlDF)
		{
			for (zr=0;zr<g_iCDFChannels;zr++)
			{
				if (o->m_pPlDF[zr] == NULL)
					continue;
				mprintf(WHITE,"\n* Plane Distance Distribution Function\n");
				mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pPlDF[zr]->m_pPlDF->m_fBinEntries,o->m_pPlDF[zr]->m_pPlDF->m_fSkipEntries,ZeroDivide(o->m_pPlDF[zr]->m_pPlDF->m_fSkipEntries,o->m_pPlDF[zr]->m_pPlDF->m_fBinEntries+o->m_pPlDF[zr]->m_pPlDF->m_fSkipEntries)*100.0);
				o->m_pPlDF[zr]->m_pPlDF->CalcMinMax();
				mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pPlDF[zr]->m_pPlDF->m_fMaxEntry,o->m_pPlDF[zr]->m_pPlDF->m_fEps);
				if (o->m_pPlDF[zr]->m_pPlDF->m_fEps > 1.0E-4) {
					eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
					mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
				}
				o->m_pPlDF[zr]->m_pPlDF->CalcMeanSD();
				mprintf("    Mean value: %10G pm    Standard deviation: %10G pm\n",o->m_pPlDF[zr]->m_pPlDF->m_fMean,o->m_pPlDF[zr]->m_pPlDF->m_fSD);
				mprintf("    Min. value: %10G pm    Max. value:         %10G pm\n",o->m_pPlDF[zr]->m_pPlDF->m_fMinInput,o->m_pPlDF[zr]->m_pPlDF->m_fMaxInput);
				o->m_pPlDF[zr]->m_pPlDF->NormBinIntegral();
//				sprintf(buf,"pldf_%s%s.csv",o->m_pPlDF[zr]->m_sName,multibuf);
				buf.sprintf("pldf_%s%s.csv",o->m_pPlDF[zr]->m_sName,multibuf);
				mprintf("    Saving plane distance distribution as %s ...\n",(const char*)buf);
				o->m_pPlDF[zr]->m_pPlDF->Write("",buf,"",false);
//				sprintf(buf,"pldf_%s%s.agr",o->m_pPlDF[zr]->m_sName,multibuf);
				buf.sprintf("pldf_%s%s.agr",o->m_pPlDF[zr]->m_sName,multibuf);
				mprintf("    Saving plane distance distribution AGR file as \"%s\"...\n",(const char*)buf);
				o->m_pPlDF[zr]->m_pPlDF->WriteAgr("",buf,"",o->m_pPlDF[zr]->m_sName,false);
				if (o->m_pPlDF[zr]->m_iHistogramRes != 0)
				{
					mprintf("    Calculating Histogram...\n");
					o->m_pPlDF[zr]->m_pPlDF->CalcHistogram();
//					sprintf(buf,"his_pldf_%s%s.agr",o->m_pPlDF[zr]->m_sName,multibuf);
					buf.sprintf("his_pldf_%s%s.csv",o->m_pPlDF[zr]->m_sName,multibuf);
					mprintf("    Saving Histogram as \"%s\"...\n",(const char*)buf);
					o->m_pPlDF[zr]->m_pPlDF->WriteHistogram("",buf,"");
				}
			}
		} // END IF PlDF

		if (g_bLiDF)
		{
			for (zr=0;zr<g_iCDFChannels;zr++)
			{
				if (o->m_pLiDF[zr] == NULL)
					continue;
				mprintf(WHITE,"\n* Line Distance Distribution Function\n");
				mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pLiDF[zr]->m_pLiDF->m_fBinEntries,o->m_pLiDF[zr]->m_pLiDF->m_fSkipEntries,ZeroDivide(o->m_pLiDF[zr]->m_pLiDF->m_fSkipEntries,o->m_pLiDF[zr]->m_pLiDF->m_fBinEntries+o->m_pLiDF[zr]->m_pLiDF->m_fSkipEntries)*100.0);
				o->m_pLiDF[zr]->m_pLiDF->CalcMinMax();
				mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pLiDF[zr]->m_pLiDF->m_fMaxEntry,o->m_pLiDF[zr]->m_pLiDF->m_fEps);
				if (o->m_pLiDF[zr]->m_pLiDF->m_fEps > 1.0E-4) {
					eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
					mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
				}
				o->m_pLiDF[zr]->m_pLiDF->CalcMeanSD();
				mprintf("    Mean value: %10G pm    Standard deviation: %10G pm\n",o->m_pLiDF[zr]->m_pLiDF->m_fMean,o->m_pLiDF[zr]->m_pLiDF->m_fSD);
				mprintf("    Min. value: %10G pm    Max. value:         %10G pm\n",o->m_pLiDF[zr]->m_pLiDF->m_fMinInput,o->m_pLiDF[zr]->m_pLiDF->m_fMaxInput);
	
				if (o->m_pLiDF[zr]->m_bRadialCorrect)
				{
					mprintf("    Correcting radial distribution...\n");
					o->m_pLiDF[zr]->m_pLiDF->CorrectLiRadialDist();
				}
				if (o->m_bOthers && o->m_pLiDF[zr]->m_bRadialCorrect)
				{
					mprintf("    Scaling LiDF to uniform density...\n");
					if (g_bPeriodicX && g_bPeriodicY && g_bPeriodicZ)
					{
						if (o->m_bObsCertain)
						{
							o->m_pLiDF[zr]->m_pLiDF->Integrate(true,1.0 / g_iSteps / o->m_waObsRefList.GetSize());
							o->m_pLiDF[zr]->m_pLiDF->MultiplyBin(g_fBoxX*g_fBoxY*g_fBoxZ / (4.0/3.0*Pi) / g_iSteps / o->m_waObsRefList.GetSize() / o->m_pLiDF[zr]->m_iShowAtomGes / o->m_waObsShowList.GetSize() / o->m_pLiDF[zr]->m_iRefAtomGes);
						} else
						{
							if (g_bBoxNonOrtho) {
								o->m_pLiDF[zr]->m_pLiDF->Integrate(true,1.0 / g_iSteps / (double)((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize());
								o->m_pLiDF[zr]->m_pLiDF->MultiplyBin(g_fBoxVolume*1000000.0 / (4.0/3.0*Pi) / g_iSteps / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pLiDF[zr]->m_iShowAtomGes / ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize() / o->m_pLiDF[zr]->m_iRefAtomGes);
							} else {
								o->m_pLiDF[zr]->m_pLiDF->Integrate(true,1.0 / g_iSteps / (double)((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize());
								o->m_pLiDF[zr]->m_pLiDF->MultiplyBin(g_fBoxX*g_fBoxY*g_fBoxZ / (4.0/3.0*Pi) / g_iSteps / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pLiDF[zr]->m_iShowAtomGes / ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize() / o->m_pLiDF[zr]->m_iRefAtomGes);
							}
						}
						if (g_bDoubleBox)
						{
							o->m_pLiDF[zr]->m_pLiDF->MultiplyBin(g_iDoubleBoxFactor);
							o->m_pLiDF[zr]->m_pLiDF->MultiplyIntegral(g_iDoubleBoxFactor);
						}
					} else
					{
						eprintf("    Uniform density not defined if box is not XYZ-periodic!\n");
						goto _lidfint1;
					}
				} else
				{
_lidfint1:
					mprintf("    Scaling LiDF to integral value 1000 ...\n");
					o->m_pLiDF[zr]->m_pLiDF->NormBinIntegral(1000.0);
					o->m_pLiDF[zr]->m_pLiDF->Integrate(false,1000.0);
				}
				
				
		//		o->m_pLiDF[zr]->m_pLiDF->NormBinIntegral();
//				sprintf(buf,"lidf_%s%s.csv",o->m_pLiDF[zr]->m_sName,multibuf);
				buf.sprintf("lidf_%s%s.csv",o->m_pLiDF[zr]->m_sName,multibuf);
				mprintf("    Saving line distance distribution as %s ...\n",(const char*)buf);
				o->m_pLiDF[zr]->m_pLiDF->Write("",buf,"",false);
//				sprintf(buf,"lidf_%s%s.agr",o->m_pLiDF[zr]->m_sName,multibuf);
				buf.sprintf("lidf_%s%s.agr",o->m_pLiDF[zr]->m_sName,multibuf);
				mprintf("    Saving line distance distribution AGR file as \"%s\"...\n",(const char*)buf);
				o->m_pLiDF[zr]->m_pLiDF->WriteAgr("",buf,"",o->m_pLiDF[zr]->m_sName,false);
				if (o->m_pLiDF[zr]->m_iHistogramRes != 0)
				{
					mprintf("    Calculating Histogram...\n");
					o->m_pLiDF[zr]->m_pLiDF->CalcHistogram();
//					sprintf(buf,"his_lidf_%s%s.agr",o->m_pLiDF[zr]->m_sName,multibuf);
					buf.sprintf("his_lidf_%s%s.csv",o->m_pLiDF[zr]->m_sName,multibuf);
					mprintf("    Saving Histogram as \"%s\"...\n",(const char*)buf);
					o->m_pLiDF[zr]->m_pLiDF->WriteHistogram("",buf,"");
				}
			}
		} // END IF LiDF

		if (g_bADF)
		{
			for (zr=0;zr<g_iCDFChannels;zr++)
			{
				if (o->m_pADF[zr] == NULL)
					continue;
				mprintf(WHITE,"\n* Angular Distribution Function\n");
				mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pADF[zr]->m_pADF->m_fBinEntries,o->m_pADF[zr]->m_pADF->m_fSkipEntries,ZeroDivide(o->m_pADF[zr]->m_pADF->m_fSkipEntries,o->m_pADF[zr]->m_pADF->m_fBinEntries+o->m_pADF[zr]->m_pADF->m_fSkipEntries)*100.0);
				o->m_pADF[zr]->m_pADF->CalcMinMax();
				mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pADF[zr]->m_pADF->m_fMaxEntry,o->m_pADF[zr]->m_pADF->m_fEps);
				if (o->m_pADF[zr]->m_pADF->m_fEps > 1.0E-4) {
					eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
					mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
				}
				o->m_pADF[zr]->m_pADF->CalcMeanSD();
				if (o->m_pADF[zr]->m_bCosine)
				{
					mprintf("    Mean value: %10G    Standard deviation: %10G\n",o->m_pADF[zr]->m_pADF->m_fMean,o->m_pADF[zr]->m_pADF->m_fSD);
					mprintf("    Min. value: %10G    Max. value:         %10G\n",o->m_pADF[zr]->m_pADF->m_fMinInput,o->m_pADF[zr]->m_pADF->m_fMaxInput);
				} else
				{
					mprintf("    Mean value: %10G degree    Standard deviation: %10G degree\n",o->m_pADF[zr]->m_pADF->m_fMean,o->m_pADF[zr]->m_pADF->m_fSD);
					mprintf("    Min. value: %10G degree    Max. value:         %10G degree\n",o->m_pADF[zr]->m_pADF->m_fMinInput,o->m_pADF[zr]->m_pADF->m_fMaxInput);
				}
				if (o->m_pADF[zr]->m_bStat)
				{
					mprintf("    Applying cone correction...\n");
					o->m_pADF[zr]->m_pADF->AngleCorrect();
				}
				if (o->m_pADF[zr]->m_bMirror)
				{
					mprintf("    Making ADF mirror-symmetric...\n");
					if (o->m_pADF[zr]->m_bCosine)
						o->m_pADF[zr]->m_pADF->Mirror(0.0);
							else o->m_pADF[zr]->m_pADF->Mirror(90.0);
				}
				o->m_pADF[zr]->m_pADF->NormBinIntegral();
				if (o->m_bTimeDev)
				{
					if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
					{
						fclose(o->m_pADF[zr]->m_fAngle[0]);
						mprintf("    Saving temporal development as adf_timedev_%s%s.csv\n",o->m_pADF[zr]->m_sName,multibuf);
					} else
					{
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
							fclose(o->m_pADF[zr]->m_fAngle[z2]);
//							sprintf(buf,"adf_timedev_%s_ref%d%s.csv",o->m_pADF[zr]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							buf.sprintf("adf_timedev_%s_ref%d%s.csv",o->m_pADF[zr]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							mprintf("    Saving temporal development as %s\n",(const char*)buf);
						}
					}
	//				mprintf("        Speichere Winkel-Statistik als %s_angle_stat.txt\n",o->m_pADF->m_sName);
	//				FreeFileName("",o->m_pADF->m_sName,"_angle_stat.txt");
	//				o->m_pADF->m_pAngleStat->Evaluate();
	//				o->m_pADF->m_pAngleStat->Write("",o->m_pADF->m_sName,"_angle_stat.txt");
	//				delete o->m_pAngleStat;
					delete[] o->m_pADF[zr]->m_fAngle;
					o->m_pADF[zr]->m_fAngle = NULL;
					if (o->m_bCombinedPlot)
					{
//						sprintf(buf,"combined_%s%s.agr",o->m_pADF[zr]->m_sName,multibuf);
						buf.sprintf("combined_%s%s.agr",o->m_pADF[zr]->m_sName,multibuf);
						mprintf("    Saving combined plot as \"%s\"...\n",(const char*)buf);
						o->m_pADF[zr]->m_pADF->CreateCombinedPlot(false);
						o->m_pADF[zr]->m_pADF->m_pCombinedPlot->WriteAgr(buf,false);
					}
				} // END IF TIMEDEV
//				sprintf(buf,"adf_%s%s.csv",o->m_pADF[zr]->m_sName,multibuf);
				buf.sprintf("adf_%s%s.csv",o->m_pADF[zr]->m_sName,multibuf);
				mprintf("    Saving ADF as %s ...\n",(const char*)buf);
				o->m_pADF[zr]->m_pADF->NormBinIntegral();
				o->m_pADF[zr]->m_pADF->Write("",buf,"",false);
//				sprintf(buf,"adf_%s%s.agr",o->m_pADF[zr]->m_sName,multibuf);
				buf.sprintf("adf_%s%s.agr",o->m_pADF[zr]->m_sName,multibuf);
				mprintf("    Saving ADF AGR file as \"%s\"...\n",(const char*)buf);
				if (o->m_pADF[zr]->m_bCosine)
					o->m_pADF[zr]->m_pADF->WriteAgr("",buf,"",o->m_pADF[zr]->m_sName,false);
						else o->m_pADF[zr]->m_pADF->WriteAgr("",buf,"",o->m_pADF[zr]->m_sName,false);
				if (o->m_pADF[zr]->m_iHistogramRes != 0)
				{
					mprintf("    Calculating Histogram...\n");
					o->m_pADF[zr]->m_pADF->CalcHistogram();
//					sprintf(buf,"his_adf_%s%s.csv",o->m_pADF[zr]->m_sName,multibuf);
					buf.sprintf("his_adf_%s%s.csv",o->m_pADF[zr]->m_sName,multibuf);
					mprintf("    Saving Histogram as \"%s\"...\n",(const char*)buf);
					o->m_pADF[zr]->m_pADF->WriteHistogram("",buf,"");
				}

				if (o->m_bObsCertain && o->m_bDecompDist)
				{
//					sprintf(buf,"adf_decomp_%s%s.csv",o->m_pADF[zr]->m_sName,multibuf);
					buf.sprintf("adf_decomp_%s%s.csv",o->m_pADF[zr]->m_sName,multibuf);
					mprintf("    Saving ADF decomposition as \"%s\"...\n",(const char*)buf);
					o->m_pADF[zr]->m_pADF->WriteMulti("",buf,"");
//					sprintf(buf,"adf_decomp_%s%s_cumulative.csv",o->m_pADF[zr]->m_sName,multibuf);
					buf.sprintf("adf_decomp_%s%s_cumulative.csv",o->m_pADF[zr]->m_sName,multibuf);
					mprintf("    Saving cumulative ADF decomposition as \"%s\"...\n",(const char*)buf);
					o->m_pADF[zr]->m_pADF->WriteMulti_Cumulative("",buf,"");
				}

				if (o->m_bTimeDiff)
					o->WriteTimeDiff(o->m_pADF[zr]->m_pADF,"ADF","adf",o->m_pADF[zr]->m_sName,multibuf,false);
			}
		} // END IF ADF

		if (g_bDDF)
		{
			for (zr=0;zr<g_iCDFChannels;zr++)
			{
				if (o->m_pDDF[zr] == NULL)
					continue;
				mprintf(WHITE,"\n* Dihedral Distribution Function\n");
				mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pDDF[zr]->m_pDDF->m_fBinEntries,o->m_pDDF[zr]->m_pDDF->m_fSkipEntries,ZeroDivide(o->m_pDDF[zr]->m_pDDF->m_fSkipEntries,o->m_pDDF[zr]->m_pDDF->m_fBinEntries+o->m_pDDF[zr]->m_pDDF->m_fSkipEntries)*100.0);
				o->m_pDDF[zr]->m_pDDF->CalcMinMax();
				mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pDDF[zr]->m_pDDF->m_fMaxEntry,o->m_pDDF[zr]->m_pDDF->m_fEps);
				if (o->m_pDDF[zr]->m_pDDF->m_fEps > 1.0E-4) {
					eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
					mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
				}
				o->m_pDDF[zr]->m_pDDF->CalcMeanSD();
				if (o->m_pDDF[zr]->m_bCosine)
				{
					mprintf("    Mean value: %10G    Standard deviation: %10G\n",o->m_pDDF[zr]->m_pDDF->m_fMean,o->m_pDDF[zr]->m_pDDF->m_fSD);
					mprintf("    Min. value: %10G    Max.value:          %10G\n",o->m_pDDF[zr]->m_pDDF->m_fMinInput,o->m_pDDF[zr]->m_pDDF->m_fMaxInput);
				} else
				{
					mprintf("    Mean value: %10G degree    Standard deviation: %10G degree\n",o->m_pDDF[zr]->m_pDDF->m_fMean,o->m_pDDF[zr]->m_pDDF->m_fSD);
					mprintf("    Min. value: %10G degree    Max.value:          %10G degree\n",o->m_pDDF[zr]->m_pDDF->m_fMinInput,o->m_pDDF[zr]->m_pDDF->m_fMaxInput);
				}
				if (o->m_bTimeDev)
				{
					if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
					{
						fclose(o->m_pDDF[zr]->m_fAngle[0]);
						mprintf("    Saving temporal development as ddf_timedev_%s%s.csv\n",o->m_pDDF[zr]->m_sName,multibuf);
					} else
					{
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
							fclose(o->m_pDDF[zr]->m_fAngle[z2]);
//							sprintf(buf,"ddf_timedev_%s_ref%d%s.csv",o->m_pDDF[zr]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							buf.sprintf("ddf_timedev_%s_ref%d%s.csv",o->m_pDDF[zr]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							mprintf("    Saving temporal development as %s\n",(const char*)buf);
						}
					}
	/*				mprintf("        Speichere Winkel-Statistik als %s_angle_stat.txt\n",o->m_sName);
					FreeFileName("",o->m_sName,"_angle_stat.txt");
					o->m_pAngleStat->Evaluate();
					o->m_pAngleStat->Write("",o->m_sName,"_angle_stat.txt");
	//				delete o->m_pAngleStat;*/
					delete[] o->m_pDDF[zr]->m_fAngle;
					o->m_pDDF[zr]->m_fAngle = NULL;
					if (o->m_bCombinedPlot)
					{
//						sprintf(buf,"combined_%s%s.agr",o->m_pDDF[zr]->m_sName,multibuf);
						buf.sprintf("combined_%s%s.agr",o->m_pDDF[zr]->m_sName,multibuf);
						mprintf("    Saving combined plot as \"%s\"...\n",(const char*)buf);
						o->m_pDDF[zr]->m_pDDF->CreateCombinedPlot(false);
						o->m_pDDF[zr]->m_pDDF->m_pCombinedPlot->WriteAgr(buf,false);
					}
				} // END IF TIMEDEV
	/*			if (o->m_pDDF->m_iStat != 0)
				{
					mprintf("        Korrigiere statistische Verteilung...\n");
					o->m_pDDF->m_pDDF->AngleCorrect();
				}*/
//				sprintf(buf,"ddf_%s%s.csv",o->m_pDDF[zr]->m_sName,multibuf);
				buf.sprintf("ddf_%s%s.csv",o->m_pDDF[zr]->m_sName,multibuf);
				mprintf("    Saving DDF as %s ...\n",(const char*)buf);
				o->m_pDDF[zr]->m_pDDF->NormBinIntegral();
				o->m_pDDF[zr]->m_pDDF->Write("",buf,"",false);
//				sprintf(buf,"ddf_%s%s.agr",o->m_pDDF[zr]->m_sName,multibuf);
				buf.sprintf("ddf_%s%s.agr",o->m_pDDF[zr]->m_sName,multibuf);
				mprintf("    Saving DDF AGR file as \"%s\"...\n",(const char*)buf);
				if (o->m_pDDF[zr]->m_bCosine)
					o->m_pDDF[zr]->m_pDDF->WriteAgr("",buf,"",o->m_pDDF[zr]->m_sName,false);
						else o->m_pDDF[zr]->m_pDDF->WriteAgr("",buf,"",o->m_pDDF[zr]->m_sName,false);
				if (o->m_pDDF[zr]->m_iHistogramRes != 0)
				{
					mprintf("    Calculating Histogram...\n");
					o->m_pDDF[zr]->m_pDDF->CalcHistogram();
//					sprintf(buf,"his_ddf_%s%s.csv",o->m_pDDF[zr]->m_sName,multibuf);
					buf.sprintf("his_ddf_%s%s.csv",o->m_pDDF[zr]->m_sName,multibuf);
					mprintf("    Saving Histogram as \"%s\"...\n",(const char*)buf);
					o->m_pDDF[zr]->m_pDDF->WriteHistogram("",buf,"");
				}

				if (o->m_bObsCertain && o->m_bDecompDist)
				{
//					sprintf(buf,"ddf_decomp_%s%s.csv",o->m_pDDF[zr]->m_sName,multibuf);
					buf.sprintf("ddf_decomp_%s%s.csv",o->m_pDDF[zr]->m_sName,multibuf);
					mprintf("    Saving DDF decomposition as \"%s\"...\n",(const char*)buf);
					o->m_pDDF[zr]->m_pDDF->WriteMulti("",buf,"");
//					sprintf(buf,"ddf_decomp_%s%s_cumulative.csv",o->m_pDDF[zr]->m_sName,multibuf);
					buf.sprintf("ddf_decomp_%s%s_cumulative.csv",o->m_pDDF[zr]->m_sName,multibuf);
					mprintf("    Saving cumulative DDF decomposition as \"%s\"...\n",(const char*)buf);
					o->m_pDDF[zr]->m_pDDF->WriteMulti_Cumulative("",buf,"");
				}

				if (o->m_bTimeDiff)
					o->WriteTimeDiff(o->m_pDDF[zr]->m_pDDF,"DDF","ddf",o->m_pDDF[zr]->m_sName,multibuf,true);
			}
		} // END IF DDF

		if (g_bSDF)
		{	
			mprintf(WHITE,"\n* Spatial Distribution Function\n");
			mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pSDF->m_pSDF->m_fBinEntries,o->m_pSDF->m_pSDF->m_fSkipEntries,ZeroDivide(o->m_pSDF->m_pSDF->m_fSkipEntries,o->m_pSDF->m_pSDF->m_fBinEntries+o->m_pSDF->m_pSDF->m_fSkipEntries)*100.0);
			o->m_pSDF->m_pSDF->CalcMaxEntry();
			mprintf("    Raw data range from %.0f to %.0f hits.\n",o->m_pSDF->m_pSDF->m_fMinEntry,o->m_pSDF->m_pSDF->m_fMaxEntry);
			mprintf("    The volume of one bin is %.3f pm^3.\n",pow3(o->m_pSDF->m_fRadius*2.0/o->m_pSDF->m_iResolution));
			o->m_pSDF->m_pSDF->CalcMaxEntry();
			mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pSDF->m_pSDF->m_fMaxEntry,o->m_pSDF->m_pSDF->m_fEps);
			if (o->m_pSDF->m_pSDF->m_fEps > 1.0E-4) {
				eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
				mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
			}
			if (g_bPeriodicX && g_bPeriodicY && g_bPeriodicZ)
				mprintf("    The uniform particle density of this observation is %.6f nm^-3.\n",o->m_pSDF->m_fParticleDensity);
					else mprintf("    The uniform particle density of non-periodic boxes is not defined.\n");
/*			switch(g_iSDFScale)
			{
				case 0:
					mprintf("    Scaling values to ppm...\n");
					o->m_pSDF->m_pSDF->MultiplyBin(1000000.0 / (double)g_iSteps / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pSDF->m_iShowAtomGes * (g_bDoubleBox?g_iDoubleBoxFactor:1));
					break;
				case 1:
					mprintf("    Scaling values to pm^-3...\n");
					o->m_pSDF->m_pSDF->MultiplyBin(pow(o->m_pSDF->m_iResolution/o->m_pSDF->m_fRadius,3) / (double)g_iSteps / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pSDF->m_iShowAtomGes * (g_bDoubleBox?g_iDoubleBoxFactor:1));
					break;
				case 2:
					mprintf("    Scaling values to nm^-3...\n");*/
/*					break;
				case 3:
					mprintf("    Scaling values relative to average particle density...\n");
					o->m_pSDF->m_pSDF->MultiplyBin(pow(o->m_pSDF->m_iResolution/o->m_pSDF->m_fRadius*1000.0,3) / (double)g_iSteps / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pSDF->m_iShowAtomGes / o->m_pSDF->m_fParticleDensity * (g_bDoubleBox?g_iDoubleBoxFactor:1));
					break;
			}*/

			o->m_pSDF->m_pSDF->MultiplyBin(pow3(o->m_pSDF->m_iResolution/o->m_pSDF->m_fRadius/2.0*1000.0) / (double)g_iSteps * g_iStride / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pSDF->m_iShowAtomGes * (g_bDoubleBox?g_iDoubleBoxFactor:1));

			if (g_bSDFUniform)
			{
				if ((!g_bPeriodicX) || (!g_bPeriodicY) || (!g_bPeriodicZ))
					goto _sdf_nm;
				o->m_pSDF->m_pSDF->MultiplyBin(1.0/o->m_pSDF->m_fParticleDensity);
				mprintf(WHITE,"    Bin values are given relative to this uniform particle density.\n");
			} else
			{
_sdf_nm:
				mprintf(WHITE,"    Bin values are given in nm^-3 (particle density).\n");
			}
			o->m_pSDF->m_pSDF->CalcMaxEntry();
			mprintf("    Data range from %.6f to %.6f%s.\n",o->m_pSDF->m_pSDF->m_fMinEntry,o->m_pSDF->m_pSDF->m_fMaxEntry,g_bSDFUniform?"":" nm^-3");

			if (o->m_pSDF->m_bInvert)
			{
				mprintf("    Inverting SDF...\n");
				o->m_pSDF->m_pSDF->Invert();
			}
			if (o->m_pSDF->m_bCutPlane)
			{
				mprintf("    Creating Cut Plane...\n");
				o->m_pSDF->CreateCutPlane();
//				sprintf(buf,"sdf_cut_%s%s",o->m_pSDF->m_sName,multibuf);
				buf.sprintf("sdf_cut_%s%s",o->m_pSDF->m_sName,multibuf);
				mprintf("    Saving SDF Cut Plane triples as \"%s_triples.csv\"...\n",(const char*)buf);
				o->m_pSDF->m_pCutPlane->Write("",buf,"_triples.csv");
				mprintf("    Saving SDF Cut Plane matrix as \"%s_matrix.csv\"...\n",(const char*)buf);
				o->m_pSDF->m_pCutPlane->WriteCSV("",buf,"_matrix.csv");
				mprintf("    Saving CDF Mathematica Notebook as \"%s.nb\"...\n",(const char*)buf);
				o->m_pSDF->m_pCutPlane->WriteMathematicaNb("",buf,".nb",false);
				mprintf("    Saving CDF Gnuplot Input as \"%s.gp\"...\n",(const char*)buf);
				o->m_pSDF->m_pCutPlane->WriteGnuplotInput("",buf,"",false);
			}
			for (z2=0;z2<=g_iSDFSmoothGrade;z2++)
			{
				try { tempSDF = new C3DF<double>(); } catch(...) { tempSDF = NULL; }
				if (tempSDF == NULL) NewException((double)sizeof(C3DF<double>),__FILE__,__LINE__,__PRETTY_FUNCTION__);
				
				tempSDF->CopyFrom(o->m_pSDF->m_pSDF);

				if (z2 != 0)
					tempSDF->Smooth(z2);

				tempSDF->CalcMaxEntry();
				mprintf("        Data range from %.6f to %.6f%s.\n",tempSDF->m_fMinEntry,tempSDF->m_fMaxEntry,g_bSDFUniform?"":" nm^-3");
				if (o->m_pSDF->m_bClipPlane)
				{
					mprintf("        Creating Clip Plane in %c direction with value %.3f...\n",(o->m_pSDF->m_iClipDirection==0)?'X':((o->m_pSDF->m_iClipDirection==1)?'Y':'Z'),o->m_pSDF->m_fClipValue);
					o->m_pSDF->m_pSDF->ClipPlane(o->m_pSDF->m_iClipDirection,o->m_pSDF->m_fClipValue);
				}

				if (g_pDatabase->GetBool("/PLOT3D/FORMATS/WRITE_PLT"))
				{
					if (z2 != 0)
//						sprintf(buf,".s%d%s.plt",z2,multibuf);
						buf.sprintf(".s%d%s.plt",z2,multibuf);
					else
//						sprintf(buf,"%s.plt",multibuf);
						buf.sprintf("%s.plt",multibuf);
					mprintf("        Saving SDF as \"sdf_%s%s\"...\n",o->m_pSDF->m_sName,(const char*)buf);
					tempSDF->WritePLT("sdf_",o->m_pSDF->m_sName,buf,true);
				}

				if (g_pDatabase->GetBool("/PLOT3D/FORMATS/WRITE_CUBE"))
				{
					if (z2 != 0)
//						sprintf(buf,".s%d%s.cube",z2,multibuf);
						buf.sprintf(".s%d%s.cube",z2,multibuf);
					else
//						sprintf(buf,"%s.cube",multibuf);
						buf.sprintf("%s.cube",multibuf);
					mprintf("        Saving SDF as \"sdf_%s%s\"...\n",o->m_pSDF->m_sName,(const char*)buf);
					tempSDF->WriteCube("sdf_",o->m_pSDF->m_sName,buf,true);
				}

				if (o->m_pSDF->m_iHistogramRes != 0)
				{
					if (z2 != 0)
//						sprintf(buf,".s%d%s.csv",z2,multibuf);
						buf.sprintf(".s%d%s.csv",z2,multibuf);
					else
//						sprintf(buf,"%s.csv",multibuf);
						buf.sprintf("%s.csv",multibuf);
					mprintf("        Saving SDF Histogram as \"his_sdf_%s%s\"...\n",o->m_pSDF->m_sName,(const char*)buf);
					tempSDF->CalcHistogram();
					tempSDF->WriteHistogram("his_sdf_",o->m_pSDF->m_sName,buf);
				}
			}
			mprintf(YELLOW,"\n    Important: ");
			mprintf(WHITE,"The PLT/CUBE files only contain the volumetric data.\n");
			mprintf(WHITE,"    The atoms of the reference molecule are saved in the ref_*.xyz file!\n");
			mprintf(WHITE,"    You need to load *both files* into the visualization program (e.g. VMD).\n");
		} // END IF SDF

		if (g_bPlProj)
		{	
			mprintf(WHITE,"\n* Plane Projection Distribution Function\n");
			mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pPlProj->m_p2DF->m_fBinEntries,o->m_pPlProj->m_p2DF->m_fSkipEntries,ZeroDivide(o->m_pPlProj->m_p2DF->m_fSkipEntries,o->m_pPlProj->m_p2DF->m_fBinEntries+o->m_pPlProj->m_p2DF->m_fSkipEntries)*100.0);
			if (o->m_pPlProj->m_p2DF->m_fBinEntries == 0)
			{
				eprintf("    There were no bin entries. Check your function definition. Skipping this Plane Projection DF.\n\n");
				goto _skipplproj;
			}

			o->m_pPlProj->m_p2DF->CalcMaxEntry();
			mprintf("    Raw data range from %.0f to %.0f hits.\n",o->m_pPlProj->m_p2DF->m_fMinEntry,o->m_pPlProj->m_p2DF->m_fMaxEntry);
			mprintf("    The volume of one bin is %.3f pm^3.\n",(o->m_pPlProj->m_fMaxVal[0]-o->m_pPlProj->m_fMinVal[0])*(o->m_pPlProj->m_fMaxVal[1]-o->m_pPlProj->m_fMinVal[1])*(o->m_pPlProj->m_fSliceBorder[1]-o->m_pPlProj->m_fSliceBorder[0])/o->m_pPlProj->m_iResolution[0]/o->m_pPlProj->m_iResolution[1]);
			o->m_pPlProj->m_p2DF->CalcMaxEntry();
			mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pPlProj->m_p2DF->m_fMaxEntry,o->m_pPlProj->m_p2DF->m_fEps);
			if (o->m_pPlProj->m_p2DF->m_fEps > 1.0E-4) {
				eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
				mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
			}

			if (o->m_bNormalizeCondition) {

				mprintf("    Adapting normalization to condition (%.3f RMs per frame on average passed)...\n", (double)o->m_iNormalizeConditionCount / ((double)g_iSteps / g_iStride) );

				// Do the actual correction (steps * RMcount)
				o->m_pPlProj->m_p2DF->MultiplyBin( 1.0 / o->m_iNormalizeConditionCount );

				// Undo the correction that will be done below
				o->m_pPlProj->m_p2DF->MultiplyBin( (double)g_iSteps / g_iStride * (double)((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() );
			}

			if ((g_bPeriodicX) && (g_bPeriodicY) && (g_bPeriodicZ))
			{
				mprintf("    System is periodic, normalizing values to uniform density.\n");
				mprintf("    The uniform particle density of this observation is %.6f nm^-3.\n",o->m_pPlProj->m_fParticleDensity);
			//	mprintf(" (%f)\n",(double)o->m_pPlProj->m_iResolution[0]*o->m_pPlProj->m_iResolution[1]/(o->m_pPlProj->m_fMaxVal[0]-o->m_pPlProj->m_fMinVal[0])/(o->m_pPlProj->m_fMaxVal[1]-o->m_pPlProj->m_fMinVal[1])/(o->m_pPlProj->m_fSliceBorder[1]-o->m_pPlProj->m_fSliceBorder[0])*1.0e9 / (double)g_iSteps / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pPlProj->m_iShowAtomGes * (g_bDoubleBox?g_iDoubleBoxFactor:1));
			//	mprintf(" o->m_pPlProj->m_iShowAtomGes = %d\n",o->m_pPlProj->m_iShowAtomGes);
			//	mprintf(" 1/Vol = %f\n",(double)o->m_pPlProj->m_iResolution[0]*o->m_pPlProj->m_iResolution[1]/(o->m_pPlProj->m_fMaxVal[0]-o->m_pPlProj->m_fMinVal[0])/(o->m_pPlProj->m_fMaxVal[1]-o->m_pPlProj->m_fMinVal[1])/(o->m_pPlProj->m_fSliceBorder[1]-o->m_pPlProj->m_fSliceBorder[0])*1.0e9);
				o->m_pPlProj->m_p2DF->MultiplyBin((double)o->m_pPlProj->m_iResolution[0]*o->m_pPlProj->m_iResolution[1]/(o->m_pPlProj->m_fMaxVal[0]-o->m_pPlProj->m_fMinVal[0])/(o->m_pPlProj->m_fMaxVal[1]-o->m_pPlProj->m_fMinVal[1])/(o->m_pPlProj->m_fSliceBorder[1]-o->m_pPlProj->m_fSliceBorder[0])*1.0e9 / (double)g_iSteps * g_iStride / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() /* / o->m_pPlProj->m_iShowAtomGes*/ * (g_bDoubleBox?g_iDoubleBoxFactor:1));
			//	mprintf(" (%f)\n",1.0/o->m_pPlProj->m_fParticleDensity);
				o->m_pPlProj->m_p2DF->MultiplyBin(1.0/o->m_pPlProj->m_fParticleDensity);
			} else
			{
				mprintf("    System is non-periodic, normalizing integral value to %.2f.\n",1000000.0);
				o->m_pPlProj->m_p2DF->NormalizeBinIntegral(1000000.0);
			}

			o->m_pPlProj->m_p2DF->CalcMaxEntry();
			mprintf("    Resulting data range from %.6f to %.6f.\n",o->m_pPlProj->m_p2DF->m_fMinEntry,o->m_pPlProj->m_p2DF->m_fMaxEntry);


			if (o->m_pPlProj->m_bDrawAtoms) {

				std::vector<CAtomSort*> tasa;
				CAtomSort *tatomsort;

				mprintf("    Inserting reference atoms into the Plane Projection DF plot...\n");
				ti = 0;
				for (z3=0;z3<o->m_pPlProj->m_oDrawAtoms.m_oaAtoms.GetSize();z3++)
				{
					ti2 = o->m_pPlProj->m_oDrawAtoms.m_baRealAtomType[z3];
					for (z4=0;z4<((CxIntArray*)o->m_pPlProj->m_oDrawAtoms.m_oaAtoms[z3])->GetSize();z4++)
					{
						if (o->m_pPlProj->m_bAverageAtomPos)
							o->m_pPlProj->m_vaAtomPos[ti] /= o->m_pPlProj->m_iAverageCounter;

						tatomsort = new CAtomSort();
						tasa.push_back( tatomsort );

						tatomsort->m_fPos[0] = o->m_pPlProj->m_vaAtomPos[ti][0];
						tatomsort->m_fPos[1] = o->m_pPlProj->m_vaAtomPos[ti][1];
						tatomsort->m_fPos[2] = o->m_pPlProj->m_vaAtomPos[ti][2];
						if (((CAtom*)g_oaAtoms[ti2])->m_pElement->m_fRadius == 0) {
							tatomsort->m_fRadius = 25.0;
							tatomsort->m_fColor[0] = 0.4;
							tatomsort->m_fColor[1] = 0.4;
							tatomsort->m_fColor[2] = 0.4;
						} else {
							tatomsort->m_fRadius = ((CAtom*)g_oaAtoms[ti2])->m_pElement->m_fRadius*0.6;
							tatomsort->m_fColor[0] = ((CAtom*)g_oaAtoms[ti2])->m_pElement->m_iColorR/255.0;
							tatomsort->m_fColor[1] = ((CAtom*)g_oaAtoms[ti2])->m_pElement->m_iColorG/255.0;
							tatomsort->m_fColor[2] = ((CAtom*)g_oaAtoms[ti2])->m_pElement->m_iColorB/255.0;
						}

						//if (((CAtom*)g_oaAtoms[ti2])->m_pElement->m_fRadius == 0)
						//	o->m_pPlProj->m_p2DF->AddCircle(o->m_pPlProj->m_vaAtomPos[ti][0],o->m_pPlProj->m_vaAtomPos[ti][1],25.0,0.4,0.4,0.4);
						//		else o->m_pPlProj->m_p2DF->AddCircle(o->m_pPlProj->m_vaAtomPos[ti][0],o->m_pPlProj->m_vaAtomPos[ti][1],((CAtom*)g_oaAtoms[ti2])->m_pElement->m_fRadius*0.6,((CAtom*)g_oaAtoms[ti2])->m_pElement->m_iColorR/255.0,((CAtom*)g_oaAtoms[ti2])->m_pElement->m_iColorG/255.0,((CAtom*)g_oaAtoms[ti2])->m_pElement->m_iColorB/255.0);

						ti++;
					}
				}

				std::stable_sort( tasa.begin(), tasa.end(), SORT_AtomSort_Z );

				for (z3=0;z3<(int)tasa.size();z3++) {

					// Circle would exceed the plot border
					if ((tasa[z3]->m_fPos[0] - tasa[z3]->m_fRadius) < o->m_pPlProj->m_p2DF->m_fMinVal[0])
						continue;

					if ((tasa[z3]->m_fPos[0] + tasa[z3]->m_fRadius) > o->m_pPlProj->m_p2DF->m_fMaxVal[0])
						continue;

					if ((tasa[z3]->m_fPos[1] - tasa[z3]->m_fRadius) < o->m_pPlProj->m_p2DF->m_fMinVal[1])
						continue;

					if ((tasa[z3]->m_fPos[1] + tasa[z3]->m_fRadius) > o->m_pPlProj->m_p2DF->m_fMaxVal[1])
						continue;

					o->m_pPlProj->m_p2DF->AddCircle(
						tasa[z3]->m_fPos[0],
						tasa[z3]->m_fPos[1],
						tasa[z3]->m_fRadius,
						tasa[z3]->m_fColor[0],
						tasa[z3]->m_fColor[1],
						tasa[z3]->m_fColor[2]
					);
				}
			}

			if (o->m_pPlProj->m_bVector) {
				mprintf("    Preparing vector field data...\n");
				o->m_pPlProj->m_p2DF->m_bVectorField = true;
				o->m_pPlProj->m_p2DF->m_iVectorRes[0] = o->m_pPlProj->m_iVecRes[0];
				o->m_pPlProj->m_p2DF->m_iVectorRes[1] = o->m_pPlProj->m_iVecRes[1];
				o->m_pPlProj->m_p2DF->m_faVectorData.resize(o->m_pPlProj->m_iVecRes[0]*o->m_pPlProj->m_iVecRes[1]*2);

				if (o->m_bNormalizeCondition)
					ti = (int)(o->m_pPlProj->m_fVecThresh * o->m_iNormalizeConditionCount);
				else
					ti = (int)(o->m_pPlProj->m_fVecThresh * g_iSteps * ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize());

				mprintf("    Using a noise rejection count of %d.\n",ti);

				for (z3=0;z3<o->m_pPlProj->m_iVecRes[1];z3++) {
					for (z4=0;z4<o->m_pPlProj->m_iVecRes[0];z4++) {
						if (o->m_pPlProj->m_faVectorCount[z3*o->m_pPlProj->m_iVecRes[0]+z4] > ti) {
							o->m_pPlProj->m_p2DF->m_faVectorData[(z3*o->m_pPlProj->m_iVecRes[0]+z4)*2] =
								o->m_pPlProj->m_faVectorData[(z3*o->m_pPlProj->m_iVecRes[0]+z4)*2] /
								o->m_pPlProj->m_faVectorCount[z3*o->m_pPlProj->m_iVecRes[0]+z4];
							o->m_pPlProj->m_p2DF->m_faVectorData[(z3*o->m_pPlProj->m_iVecRes[0]+z4)*2+1] =
								o->m_pPlProj->m_faVectorData[(z3*o->m_pPlProj->m_iVecRes[0]+z4)*2+1] /
								o->m_pPlProj->m_faVectorCount[z3*o->m_pPlProj->m_iVecRes[0]+z4];
						} else {
							o->m_pPlProj->m_p2DF->m_faVectorData[(z3*o->m_pPlProj->m_iVecRes[0]+z4)*2] = 0;
							o->m_pPlProj->m_p2DF->m_faVectorData[(z3*o->m_pPlProj->m_iVecRes[0]+z4)*2+1] = 0;
						}
					}
				}
				mprintf("    Note that the vector field plot currently works *only* with Gnuplot!\n");
			}

			mprintf("    Saving Plane Projection DF triples as \"%s%s_triples.csv\"...\n",o->m_pPlProj->m_sName,multibuf);
			o->m_pPlProj->m_p2DF->Write(o->m_pPlProj->m_sName,multibuf,"_triples.csv");
			mprintf("    Saving Plane Projection DF matrix as \"%s%s_matrix.csv\"...\n",o->m_pPlProj->m_sName,multibuf);
			o->m_pPlProj->m_p2DF->WriteCSV(o->m_pPlProj->m_sName,multibuf,"_matrix.csv");
			mprintf("    Saving Plane Projection DF Mathematica Notebook as \"%s%s.nb\"...\n",o->m_pPlProj->m_sName,multibuf);
			o->m_pPlProj->m_p2DF->WriteMathematicaNb(o->m_pPlProj->m_sName,multibuf,".nb",false);
			mprintf("    Saving Plane Projection DF Gnuplot Input as \"%s%s.gp\"...\n",o->m_pPlProj->m_sName,multibuf);
			o->m_pPlProj->m_p2DF->WriteGnuplotInput(o->m_pPlProj->m_sName,multibuf,"",false);
			if (o->m_pPlProj->m_iHistogramRes != 0)
			{
				mprintf("    Calculating Histogram...\n");
				o->m_pPlProj->m_p2DF->CalcHistogram();
//				sprintf(buf,"his_%s",o->m_pPlProj->m_sName);
				buf.sprintf("his_%s",o->m_pPlProj->m_sName);
				mprintf("    Saving Histogram as \"%s%s.csv\"...\n",(const char*)buf,multibuf);
				o->m_pPlProj->m_p2DF->WriteHistogram(buf,multibuf,".csv");
			}
			mprintf("\n");
_skipplproj:;
		} // END IF PlProj

		if (g_bRDF)
		{
			for (zr=0;zr<g_iCDFChannels;zr++)
			{
				if (o->m_pRDF[zr] == NULL)
					continue;
				mprintf(WHITE,"\n* Radial Distribution Function\n");
				mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pRDF[zr]->m_pRDF->m_fBinEntries,o->m_pRDF[zr]->m_pRDF->m_fSkipEntries,ZeroDivide(o->m_pRDF[zr]->m_pRDF->m_fSkipEntries,o->m_pRDF[zr]->m_pRDF->m_fBinEntries+o->m_pRDF[zr]->m_pRDF->m_fSkipEntries)*100.0);
				o->m_pRDF[zr]->m_pRDF->CalcMinMax();
				mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pRDF[zr]->m_pRDF->m_fMaxEntry,o->m_pRDF[zr]->m_pRDF->m_fEps);
				if (o->m_pRDF[zr]->m_pRDF->m_fEps > 1.0E-4) {
					eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
					mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
				}
				o->m_pRDF[zr]->m_pRDF->CalcMeanSD();
				mprintf("    Mean value: %10G pm    Standard deviation: %10G pm\n",o->m_pRDF[zr]->m_pRDF->m_fMean,o->m_pRDF[zr]->m_pRDF->m_fSD);
				mprintf("    Min. value: %10G pm    Max.value:          %10G pm\n",o->m_pRDF[zr]->m_pRDF->m_fMinInput,o->m_pRDF[zr]->m_pRDF->m_fMaxInput);

				o->m_pRDF[zr]->m_pRDF->MultiplyBin(1.0 / g_iSteps * g_iStride);


				if (o->m_pRDF[zr]->m_bAdaptive)
				{
//					o->m_pRDF[zr]->m_pRDF->WriteHenry("rdf_",o->m_pRDF[zr]->m_sName,".csv");
					o->m_pRDF[zr]->m_pRDF->PrepareAdapt();
					o->m_pRDF[zr]->m_pRDF->BinTree_RadialDist();
					if (o->m_iShowMol == g_iFixMol)
						o->m_pRDF[zr]->m_pRDF->BinTree_MultiplyBin(3.0/4.0/Pi/g_iSteps*g_fBoxX*g_fBoxY*g_fBoxZ/(((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()-1.0)/o->m_pRDF[zr]->m_iShowAtomGes/((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()/o->m_pRDF[zr]->m_iRefAtomGes);
							else o->m_pRDF[zr]->m_pRDF->BinTree_MultiplyBin(3.0/4.0/Pi/g_iSteps*g_fBoxX*g_fBoxY*g_fBoxZ/((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()      /o->m_pRDF[zr]->m_iShowAtomGes/((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()/o->m_pRDF[zr]->m_iRefAtomGes);

					mprintf("    Saving RDF as \"rdf_%s.csv\"...\n",o->m_pRDF[zr]->m_sName);
					o->m_pRDF[zr]->m_pRDF->WriteAdapted("rdf_",o->m_pRDF[zr]->m_sName,".001.csv",5,10,0.001,true);
					o->m_pRDF[zr]->m_pRDF->WriteAdapted("rdf_",o->m_pRDF[zr]->m_sName,".002.csv",5,10,0.002,true);
					o->m_pRDF[zr]->m_pRDF->WriteAdapted("rdf_",o->m_pRDF[zr]->m_sName,".005.csv",5,10,0.005,true);
					o->m_pRDF[zr]->m_pRDF->WriteAdapted("rdf_",o->m_pRDF[zr]->m_sName,".010.csv",5,10,0.01,true);
					o->m_pRDF[zr]->m_pRDF->WriteAdapted("rdf_",o->m_pRDF[zr]->m_sName,".020.csv",5,10,0.02,true);
					o->m_pRDF[zr]->m_pRDF->WriteAdapted("rdf_",o->m_pRDF[zr]->m_sName,".050.csv",5,10,0.05,true);
					o->m_pRDF[zr]->m_pRDF->WriteAdapted("rdf_",o->m_pRDF[zr]->m_sName,".100.csv",5,10,0.1,true);

/*					o->m_pRDF[zr]->m_pRDF->WriteAdapted("rdf_",o->m_pRDF[zr]->m_sName,".50.csv",16,0,50,true);
					o->m_pRDF[zr]->m_pRDF->WriteAdapted("rdf_",o->m_pRDF[zr]->m_sName,".100.csv",16,0,100,true);
					o->m_pRDF[zr]->m_pRDF->WriteAdapted("rdf_",o->m_pRDF[zr]->m_sName,".500.csv",16,0,500,true);
					o->m_pRDF[zr]->m_pRDF->WriteAdapted("rdf_",o->m_pRDF[zr]->m_sName,".1000.csv",16,0,1000,true);*/
				} else
				{
					if (o->m_bOthers && o->m_pRDF[zr]->m_bRadialCorrect)
					{

						if (o->m_bNormalizeCondition) {

							mprintf("    Adapting normalization to condition (%.3f RMs per frame on average passed)...\n", (double)o->m_iNormalizeConditionCount / ((double)g_iSteps / g_iStride) );

							// Do the actual correction (steps * RMcount)
							o->m_pRDF[zr]->m_pRDF->MultiplyBin( 1.0 / o->m_iNormalizeConditionCount );

							// Undo the correction from above
							o->m_pRDF[zr]->m_pRDF->MultiplyBin( (double)g_iSteps / g_iStride);

							// Undo the RM count correction that will be done below ^^
							if (o->m_bObsCertain)
								o->m_pRDF[zr]->m_pRDF->MultiplyBin( (double)o->m_waObsRefList.GetSize() );
							else
								o->m_pRDF[zr]->m_pRDF->MultiplyBin( (double)((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() );
						}

						if (o->m_pRDF[zr]->m_bProbDens)
						{
							mprintf("    Scaling RDF to nm^(-3) ...\n");
							if (o->m_bObsCertain)
							{
								o->m_pRDF[zr]->m_pRDF->Integrate( /*true*/ false, 1.0 / o->m_waObsRefList.GetSize());
								o->m_pRDF[zr]->m_pRDF->MultiplyBin(1e9 / (4.0/3.0*Pi) / o->m_waObsRefList.GetSize() / o->m_pRDF[zr]->m_iRefAtomGes);
							} else
							{
								o->m_pRDF[zr]->m_pRDF->Integrate( /*true*/ false, 1.0 / (double)((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pRDF[zr]->m_iRefAtomGes);
								o->m_pRDF[zr]->m_pRDF->MultiplyBin(1e9 / (4.0/3.0*Pi) / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pRDF[zr]->m_iRefAtomGes);
							}
							if (g_bDoubleBox)
							{
								o->m_pRDF[zr]->m_pRDF->MultiplyBin(g_iDoubleBoxFactor);
								o->m_pRDF[zr]->m_pRDF->MultiplyIntegral(g_iDoubleBoxFactor);
							}
						} else
						{
							mprintf("    Scaling RDF to uniform density ...\n");
							if (g_bPeriodicX && g_bPeriodicY && g_bPeriodicZ)
							{
								// Bugfix MB 15.04.2018 --> reverted on 21.12.2018
								//if (g_iFixMol == o->m_iShowMol)
								//	ti = ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize() - 1;
								//else
								//	ti = ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();
								
								ti = ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize();

								if (o->m_bObsCertain)
								{
									o->m_pRDF[zr]->m_pRDF->Integrate( /*true*/ false, 1.0 / o->m_waObsRefList.GetSize());
									o->m_pRDF[zr]->m_pRDF->MultiplyBin(g_fBoxX*g_fBoxY*g_fBoxZ / (4.0/3.0*Pi) / o->m_waObsRefList.GetSize() / o->m_pRDF[zr]->m_iShowAtomGes / o->m_waObsShowList.GetSize() / o->m_pRDF[zr]->m_iRefAtomGes);
								} else
								{
									if (g_bRegionAnalysis)
									{
										o->m_pRDF[zr]->m_pRDF->Integrate( /*true*/ false, 1.0 * g_iSteps / (o->m_pRDF[zr]->m_pRDF->m_fBinEntries + o->m_pRDF[zr]->m_pRDF-> m_fSkipEntries) * o->m_pRDF[zr]->m_iShowAtomGes * ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize() * o->m_pRDF[zr]->m_iRefAtomGes);
										o->m_pRDF[zr]->m_pRDF->MultiplyBin(g_fBoxX*g_fBoxY*g_fBoxZ * g_iSteps / (4.0/3.0*Pi) / (o->m_pRDF[zr]->m_pRDF->m_fBinEntries + o->m_pRDF[zr]->m_pRDF-> m_fSkipEntries) );
									} else
									{
										if (g_bBoxNonOrtho)
										{
											o->m_pRDF[zr]->m_pRDF->Integrate( /*true*/ false, 1.0 / (double)((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pRDF[zr]->m_iRefAtomGes);
											o->m_pRDF[zr]->m_pRDF->MultiplyBin(g_fBoxVolume*1000000.0 / (4.0/3.0*Pi) / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pRDF[zr]->m_iShowAtomGes / ti / o->m_pRDF[zr]->m_iRefAtomGes);
										} else
										{
											o->m_pRDF[zr]->m_pRDF->Integrate( /*true*/ false, 1.0 / (double)((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pRDF[zr]->m_iRefAtomGes);
											o->m_pRDF[zr]->m_pRDF->MultiplyBin(g_fBoxX*g_fBoxY*g_fBoxZ / (4.0/3.0*Pi) / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / o->m_pRDF[zr]->m_iShowAtomGes / ti / o->m_pRDF[zr]->m_iRefAtomGes);
										}
									}
								}
								if (g_bDoubleBox)
								{
									o->m_pRDF[zr]->m_pRDF->MultiplyBin(g_iDoubleBoxFactor);
									o->m_pRDF[zr]->m_pRDF->MultiplyIntegral(g_iDoubleBoxFactor);
								}
							} else
							{
								eprintf("    Uniform density not defined if box is not XYZ-periodic!\n");
								goto _rdfint1;
							}
						}
					} else
					{
_rdfint1:
						mprintf("    Scaling RDF to integral value 1000 ...\n");
						o->m_pRDF[zr]->m_pRDF->NormBinIntegral(1000.0);
						o->m_pRDF[zr]->m_pRDF->Integrate(false,1000.0);
					}

					if (o->m_pRDF[zr]->m_bRadialCorrect)
					{

						if (o->m_pRDF[zr]->m_bKuehneLong) {
							mprintf("    Correcting radial distribution (Roehrig & Kuehne long-range mode)...\n");
							o->m_pRDF[zr]->m_pRDF->CorrectRadialDistKuehne();
						} else if (o->m_pRDF[zr]->m_bLongMode) {
							mprintf("    Correcting radial distribution (long-range mode)...\n");
							o->m_pRDF[zr]->m_pRDF->CorrectRadialDistLong();
						} else {
							mprintf("    Correcting radial distribution...\n");
							o->m_pRDF[zr]->m_pRDF->CorrectRadialDist();
						}
					}

					mprintf("    Performing automatic peak detection...\n");
					o->m_pRDF[zr]->m_pRDF->FindFirstMinMax(false);
					if (o->m_pRDF[zr]->m_pRDF->m_faMaxPos.size() != 0) {
						mprintf(WHITE,"    >>> Peak Data >>>\n");
						for (z2=0;z2<(int)o->m_pRDF[zr]->m_pRDF->m_faMaxPos.size();z2++) {
							mprintf("      %d. maximum at %7.2f pm, height %7.4f\n",z2+1,o->m_pRDF[zr]->m_pRDF->m_faMaxPos[z2],o->m_pRDF[zr]->m_pRDF->m_faMaxVal[z2]);
							if ((int)o->m_pRDF[zr]->m_pRDF->m_faMinPos.size() > z2) {
								mprintf("      %d. minimum at %7.2f pm, height %7.4f",z2+1,o->m_pRDF[zr]->m_pRDF->m_faMinPos[z2],o->m_pRDF[zr]->m_pRDF->m_faMinVal[z2]);
								if ((int)o->m_pRDF[zr]->m_pRDF->m_faShellInt.size() > z2)
									mprintf(", %7.4f particles in %d. solvation shell.",o->m_pRDF[zr]->m_pRDF->m_faShellInt[z2],z2+1);
								mprintf("\n");
							}
						}
						mprintf(WHITE,"    <<< Peak Data <<<\n");
					} else
						mprintf("    No peaks detected.\n");
//					sprintf(buf,"rdf_%s%s.csv",o->m_pRDF[zr]->m_sName,multibuf);
					buf.sprintf("rdf_%s%s.csv",o->m_pRDF[zr]->m_sName,multibuf);
					mprintf("    Saving RDF as \"%s\"...\n",(const char*)buf);
					o->m_pRDF[zr]->m_pRDF->Write("",buf,"",true);
					mprintf(YELLOW,"      Note:");
					mprintf(" The third column of the text file contains the number integral / coordination number.\n");
//					sprintf(buf,"rdf_%s%s.agr",o->m_pRDF[zr]->m_sName,multibuf);
					buf.sprintf("rdf_%s%s.agr",o->m_pRDF[zr]->m_sName,multibuf);
					mprintf("    Saving RDF AGR file as \"%s\"...\n",(const char*)buf);
					o->m_pRDF[zr]->m_pRDF->WriteAgr("",buf,"",o->m_pRDF[zr]->m_sName,o->m_pRDF[zr]->m_bLine);
					if (o->m_pRDF[zr]->m_iHistogramRes != 0)
					{
						mprintf("    Calculating Histogram...\n");
						o->m_pRDF[zr]->m_pRDF->CalcHistogram();
//						sprintf(buf,"his_rdf_%s%s.csv",o->m_pRDF[zr]->m_sName,multibuf);
						buf.sprintf("his_rdf_%s%s.csv",o->m_pRDF[zr]->m_sName,multibuf);
						mprintf("    Saving Histogram as \"%s\"...\n",(const char*)buf);
						o->m_pRDF[zr]->m_pRDF->WriteHistogram("",buf,"");
					}
					if (o->m_bObsCertain && o->m_bDecompDist)
					{
//						sprintf(buf,"rdf_decomp_mol_%s%s",o->m_pRDF[zr]->m_sName,multibuf);
						buf.sprintf("rdf_decomp_mol_%s%s",o->m_pRDF[zr]->m_sName,multibuf);
						mprintf("    Saving RDF molecule decomposition as \"%s.csv\"...\n",(const char*)buf);
						o->m_pRDF[zr]->m_pRDF->WriteMulti("",buf,".csv");
						mprintf("    Saving RDF molecule decomposition as \"%s.agr\"...\n",(const char*)buf);
						o->m_pRDF[zr]->m_pRDF->WriteMultiAgr("",buf,".agr",buf,true);

//						sprintf(buf,"rdf_decomp_nb_%s%s_cumulative",o->m_pRDF[zr]->m_sName,multibuf);
						buf.sprintf("rdf_decomp_nb_%s%s_cumulative",o->m_pRDF[zr]->m_sName,multibuf);
						mprintf("    Saving cumulative RDF molecule decomposition as \"%s.csv\"...\n",(const char*)buf);
						o->m_pRDF[zr]->m_pRDF->WriteMulti_Cumulative("",buf,".csv");
						mprintf("    Saving cumulative RDF molecule decomposition as \"%s.agr\"...\n",(const char*)buf);
						o->m_pRDF[zr]->m_pRDF->WriteMultiAgr_Cumulative("",buf,".agr",buf,true);
					}
					if (o->m_bDecompType)
					{
//						sprintf(buf,"rdf_decomp_type_%s%s",o->m_pRDF[zr]->m_sName,multibuf);
						buf.sprintf("rdf_decomp_type_%s%s",o->m_pRDF[zr]->m_sName,multibuf);
						mprintf("    Saving RDF type decomposition as \"%s.csv\"...\n",(const char*)buf);
						o->m_pRDF[zr]->m_pRDF->WriteMulti("",buf,".csv");
						mprintf("    Saving RDF type decomposition as \"%s.agr\"...\n",(const char*)buf);
						o->m_pRDF[zr]->m_pRDF->WriteMultiAgr("",buf,".agr",buf,true);

//						sprintf(buf,"rdf_decomp_type_%s%s_cumulative",o->m_pRDF[zr]->m_sName,multibuf);
						buf.sprintf("rdf_decomp_type_%s%s_cumulative",o->m_pRDF[zr]->m_sName,multibuf);
						mprintf("    Saving cumulative RDF type decomposition as \"%s.csv\"...\n",(const char*)buf);
						o->m_pRDF[zr]->m_pRDF->WriteMulti_Cumulative("",buf,".csv");
						mprintf("    Saving cumulative RDF type decomposition as \"%s.agr\"...\n",(const char*)buf);
						o->m_pRDF[zr]->m_pRDF->WriteMultiAgr_Cumulative("",buf,".agr",buf,true);
					}

					if (o->m_bTimeDiff)
						o->WriteTimeDiff(o->m_pRDF[zr]->m_pRDF,"RDF","rdf",o->m_pRDF[zr]->m_sName,multibuf,false);

					if (o->m_pRDF[zr]->m_bACF)
					{
						o->m_pRDF[zr]->Autocorrelate();
					}
				}
				if (o->m_bTimeDev)
				{
					if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
					{
						fclose(o->m_pRDF[zr]->m_fDist[0]);
						mprintf("    Saving temporal development as \"rdf_timedev_%s%s.csv\"...\n",o->m_pRDF[zr]->m_sName,multibuf);
					} else
					{
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
							fclose(o->m_pRDF[zr]->m_fDist[z2]);
//							sprintf(buf,"rdf_timedev_%s_ref%d%s.csv",o->m_pRDF[zr]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							buf.sprintf("rdf_timedev_%s_ref%d%s.csv",o->m_pRDF[zr]->m_sName,o->m_waSaveRefList[z2]+1,multibuf);
							mprintf("    Saving temporal development as \"%s\"...\n",(const char*)buf);
						}
					}
					delete[] o->m_pRDF[zr]->m_fDist;
					o->m_pRDF[zr]->m_fDist = NULL;
					if (o->m_bCombinedPlot)
					{
//						sprintf(buf,"combined_%s%s.agr",o->m_pRDF[zr]->m_sName,multibuf);
						buf.sprintf("combined_%s%s.agr",o->m_pRDF[zr]->m_sName,multibuf);
						mprintf("    Saving combined plot as \"%s\"...\n",(const char*)buf);
						o->m_pRDF[zr]->m_pRDF->CreateCombinedPlot(true);
						o->m_pRDF[zr]->m_pRDF->m_pCombinedPlot->WriteAgr(buf,false);
					}
				}
			}
		} // END IF RDF

		if (g_bVHDF)
		{
			mprintf(WHITE,"\n* Van Hove Correlation Function\n");
			mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pVHDF->m_pVHDF->m_fBinEntries,o->m_pVHDF->m_pVHDF->m_fSkipEntries,ZeroDivide(o->m_pVHDF->m_pVHDF->m_fSkipEntries,o->m_pVHDF->m_pVHDF->m_fBinEntries+o->m_pVHDF->m_pVHDF->m_fSkipEntries)*100.0);
			o->m_pVHDF->CorrectCount();
			if (o->m_pVHDF->m_bRadialCorrect)
			{
				mprintf("    Correcting radial distribution...\n");
				o->m_pVHDF->m_pVHDF->CorrectRadialDist(1);
				if (g_bPeriodicX && g_bPeriodicY && g_bPeriodicZ)
				{
					mprintf("    Normalizing bin to uniform density...\n");
					o->m_pVHDF->m_pVHDF->MultiplyBin(3.0/4.0/Pi*g_fBoxX*g_fBoxY*g_fBoxZ);
				}
				if (g_bDoubleBox)
					o->m_pVHDF->m_pVHDF->MultiplyBin(g_iDoubleBoxFactor);
			} else
			{
				mprintf("    Normalizing bin integral to 1000000...\n");
				o->m_pVHDF->m_pVHDF->NormalizeBinIntegral(1000000.0);
			}

//			sprintf(buf,"vhcf_%s%s",o->m_pVHDF->m_sName,multibuf);
			buf.sprintf("vhcf_%s%s",o->m_pVHDF->m_sName,multibuf);

			if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_TRIPLES"))
			{
				mprintf("    Saving VHCF triples as \"%s_triples.csv\"...\n",(const char*)buf);
				o->m_pVHDF->m_pVHDF->Write("",buf,"_triples.csv");
			}

			if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATRIX"))
			{
				mprintf("    Saving VHCF matrix as \"%s_matrix.csv\"...\n",(const char*)buf);
				o->m_pVHDF->m_pVHDF->WriteCSV("",buf,"_matrix.csv");
			}

			if (o->m_pVHDF->m_iGraceBunchDist > 0)
			{
//				sprintf(buf,"vhcf_dist_%s%s.agr",o->m_pVHDF->m_sName,multibuf);
				buf.sprintf("vhcf_dist_%s%s.agr",o->m_pVHDF->m_sName,multibuf);
				mprintf("    Saving VHCF distance Grace Stack as \"%s\"...\n",(const char*)buf);
				o->m_pVHDF->m_pVHDF->WriteGraceBunch(0,o->m_pVHDF->m_iGraceBunchDist,1.0,"",buf,"");
			}

			if (o->m_pVHDF->m_iGraceBunchTime > 0)
			{
//				sprintf(buf,"vhcf_time_%s%s.agr",o->m_pVHDF->m_sName,multibuf);
				buf.sprintf("vhcf_time_%s%s.agr",o->m_pVHDF->m_sName,multibuf);
				mprintf("    Saving VHCF time Grace Stack as \"%s\"...\n",(const char*)buf);
				o->m_pVHDF->m_pVHDF->WriteGraceBunch(1,o->m_pVHDF->m_iGraceBunchTime,1.0,"",buf,"");
			}

			if (o->m_pVHDF->m_bSwapAxes) {
				mprintf("    Transposing 2D histogram...\n");
				o->m_pVHDF->m_pVHDF->SwapAxes();
			}

			if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATHEMATICA"))
			{
//				sprintf(buf,"vhcf_%s%s.nb",o->m_pVHDF->m_sName,multibuf);
				buf.sprintf("vhcf_%s%s.nb",o->m_pVHDF->m_sName,multibuf);
				mprintf("    Saving VHCF Mathematica Notebook as \"%s\"...\n",(const char*)buf);
				o->m_pVHDF->m_pVHDF->WriteMathematicaNb("",buf,"",false);
			}

			if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_GNUPLOT"))
			{
//				sprintf(buf,"vhcf_%s%s",o->m_pVHDF->m_sName,multibuf);
				buf.sprintf("vhcf_%s%s",o->m_pVHDF->m_sName,multibuf);
				mprintf("    Saving VHCF Gnuplot Input as \"%s.gp\"...\n",(const char*)buf);
				o->m_pVHDF->m_pVHDF->WriteGnuplotInput("",buf,"",false);
			}
		} // END IF VHDF

		if (g_bCDF)
		{
			if (g_iCDFChannels == 2)
			{
				mprintf(WHITE,"\n* Combined Distribution Function (2D)\n");
			/*	if (o->m_pConditions != NULL)
				{
					mprintf("    Conditions: %.2f percent of the molecules passed.\n",o->m_pConditions->m_fPassed/(o->m_pConditions->m_fPassed+o->m_pConditions->m_fFailed)*100.0);
					mprintf("    After conditions: ");
				} else mprintf("    ");*/
				mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pCDF->m_p2DF->m_fBinEntries,o->m_pCDF->m_p2DF->m_fSkipEntries,ZeroDivide(o->m_pCDF->m_p2DF->m_fSkipEntries,o->m_pCDF->m_p2DF->m_fSkipEntries+o->m_pCDF->m_p2DF->m_fBinEntries)*100.0);
				if (o->m_pCDF->m_p2DF->m_fBinEntries == 0)
				{
					eprintf("    There were no bin entries. Check your CDF definition. Skipping this CDF.\n\n");
					goto _skipsdf;
				}
				o->m_pCDF->m_p2DF->CalcMaxEntry();
				mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pCDF->m_p2DF->m_fMaxEntry,o->m_pCDF->m_p2DF->m_fEps);
				if (o->m_pCDF->m_p2DF->m_fEps > 1.0E-4) {
					eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
					mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
				}
				for (z2=0;z2<2;z2++)
				{
					if (g_iObsChannel[z2] == 0)
					{
						if ((o->m_pRDF[z2]->m_bRadialCorrect) && (o->m_pCDF->m_iNormalize != 3))
						{
							mprintf("    Correcting radial distribution for CDF channel %d...\n",z2+1);
							o->m_pCDF->m_p2DF->CorrectRadialDist(z2);
						}
					}
					if (g_iObsChannel[z2] == 1)
					{
						if (o->m_pADF[z2]->m_bStat && (!o->m_pADF[z2]->m_bCosine) && (o->m_pCDF->m_iNormalize != 3))
						{
							mprintf("    Correcting angular distribution for CDF channel %d...\n",z2+1);
							o->m_pCDF->m_p2DF->CorrectAngle(z2);
						}
						if (o->m_pADF[z2]->m_bMirror)
						{
							mprintf("    Making channel %d mirror-symmetric...\n",z2+1);
							if (o->m_pADF[z2]->m_bCosine)
								o->m_pCDF->m_p2DF->Mirror(0.0,z2);
							else
								o->m_pCDF->m_p2DF->Mirror(90.0,z2);
						}
					}
					if (g_iObsChannel[z2] == 6)
					{
						if (o->m_pLiDF[z2]->m_bRadialCorrect)
						{
							mprintf("    Correcting radial distribution for CDF channel %d...\n",z2+1);
							o->m_pCDF->m_p2DF->CorrectLiRadialDist(z2);
						}
					}
				}
				o->m_pCDF->m_fFactor = 1.0;
				switch(o->m_pCDF->m_iNormalize)
				{
					case 3:
						mprintf("    Normalizing CDF to uniform density.\n");
						o->m_pCDF->m_p2DF->NormalizeUniform(g_fBoxX*g_fBoxY*g_fBoxZ/g_iSteps  / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize() / ((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize() / o->m_pCDF->m_iCombinationsEnabled);
						break;
					case 2:
						mprintf("    Normalizing maximum value to %.2f.\n",o->m_pCDF->m_fNormValue);
						o->m_pCDF->m_p2DF->NormalizeBin(0.0,o->m_pCDF->m_fNormValue);
						break;
					case 1:
						mprintf("    Normalizing integral value to %.2f.\n",o->m_pCDF->m_fNormValue);
						o->m_pCDF->m_fFactor = o->m_pCDF->m_p2DF->NormalizeBinIntegral(o->m_pCDF->m_fNormValue);
						break;
					case 0:
						mprintf("    Not normalizing this CDF.\n");
						break;
				}

				if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_TRIPLES"))
				{
					mprintf("    Saving CDF triples as \"%s%s_triples.csv\"...\n",o->m_pCDF->m_sName,multibuf);
					o->m_pCDF->m_p2DF->Write(o->m_pCDF->m_sName,multibuf,"_triples.csv");
				}

				if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATRIX"))
				{
					mprintf("    Saving CDF matrix as \"%s%s_matrix.csv\"...\n",o->m_pCDF->m_sName,multibuf);
					o->m_pCDF->m_p2DF->WriteCSV(o->m_pCDF->m_sName,multibuf,"_matrix.csv");
				}

				if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_MATHEMATICA"))
				{
					mprintf("    Saving CDF Mathematica Notebook as \"%s%s.nb\"...\n",o->m_pCDF->m_sName,multibuf);
					o->m_pCDF->m_p2DF->WriteMathematicaNb(o->m_pCDF->m_sName,multibuf,".nb",false);
				}

				if (g_pDatabase->GetBool("/PLOT2D/FORMATS/WRITE_GNUPLOT"))
				{
					mprintf("    Saving CDF Gnuplot Input as \"%s%s.gp\"...\n",o->m_pCDF->m_sName,multibuf);
					o->m_pCDF->m_p2DF->WriteGnuplotInput(o->m_pCDF->m_sName,multibuf,"",false);
				}

				mprintf("    Saving combined plot AGR file as \"%s%s_combined.agr\"...\n",o->m_pCDF->m_sName,multibuf);
				o->m_pCDF->m_p2DF->WriteCombinedPlot(o->m_pCDF->m_sName,multibuf,"_combined.agr");
				if (o->m_pCDF->m_iHistogramRes != 0)
				{
					mprintf("    Calculating Histogram...\n");
					o->m_pCDF->m_p2DF->CalcHistogram();
//					sprintf(buf,"his_%s",o->m_pCDF->m_sName);
					buf.sprintf("his_%s",o->m_pCDF->m_sName);
					mprintf("    Saving Histogram as \"%s%s.csv\"...\n",(const char*)buf,multibuf);
					o->m_pCDF->m_p2DF->WriteHistogram(buf,multibuf,".csv");
				}

				if (o->m_pCDF->m_bGraceBunch)
				{
					if (o->m_pCDF->m_iGraceBunchC1 > 0)
					{
						if ((g_iObsChannel[1] == 0) && (o->m_bOthers))
							tf = 1.0 / o->m_pCDF->m_fFactor * 3.0/4.0/Pi/g_iSteps*g_fBoxX*g_fBoxY*g_fBoxZ/((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()      /o->m_pRDF[1]->m_iShowAtomGes/((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()/o->m_pRDF[1]->m_iRefAtomGes;
								else tf = 1;
//						sprintf(buf,"cdf_c1_%s%s.agr",o->m_pCDF->m_sName,multibuf);
						buf.sprintf("cdf_c1_%s%s.agr",o->m_pCDF->m_sName,multibuf);
						mprintf("    Saving CDF Channel 1 Grace Stack as \"%s\"...\n",(const char*)buf);
						o->m_pCDF->m_p2DF->WriteGraceBunch(0,o->m_pCDF->m_iGraceBunchC1,(double)tf,"",buf,"");
					}
					if (o->m_pCDF->m_iGraceBunchC2 > 0)
					{
						if ((g_iObsChannel[0] == 0) && (o->m_bOthers))
							tf = 1.0 / o->m_pCDF->m_fFactor * 3.0/4.0/Pi/g_iSteps*g_fBoxX*g_fBoxY*g_fBoxZ/((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()      /o->m_pRDF[0]->m_iShowAtomGes/((CMolecule*)g_oaMolecules[o->m_iShowMol])->m_laSingleMolIndex.GetSize()/o->m_pRDF[0]->m_iRefAtomGes;
								else tf = 1;
//						sprintf(buf,"cdf_c2_%s%s.agr",o->m_pCDF->m_sName,multibuf);
						buf.sprintf("cdf_c2_%s%s.agr",o->m_pCDF->m_sName,multibuf);
						mprintf("    Saving CDF Channel 2 Grace Stack as \"%s\"...\n",(const char*)buf);
						o->m_pCDF->m_p2DF->WriteGraceBunch(1,o->m_pCDF->m_iGraceBunchC2,(double)tf,"",buf,"");
					}
				}

				if (o->m_pCDF->m_bAxisDivide)
				{
					mprintf(WHITE,"    Computing correlation plot...\n");
					mprintf("    Correlation coefficient: %f\n",o->m_pCDF->m_p2DF->CalcCorrelationFactor());
					mprintf("      (values > 0 indicate positive correlation, < 0 negative correlation of the two quantities in the CDF)\n");
	
					if (o->m_pCDF->m_bAxisDivideAll)
					{
						mprintf("    Saving CDF X projection as \"%s%s_pX.csv\"...\n",o->m_pCDF->m_sName,multibuf);
						o->m_pCDF->m_p2DF->WriteXProjection(o->m_pCDF->m_sName,multibuf,"_pX.csv");

						mprintf("    Saving CDF Y projection as \"%s%s_pY.csv\"...\n",o->m_pCDF->m_sName,multibuf);
						o->m_pCDF->m_p2DF->WriteYProjection(o->m_pCDF->m_sName,multibuf,"_pY.csv");

						mprintf("    Saving X-normalized CDF Mathematica Notebook as \"%s%s_nX.nb\"...\n",o->m_pCDF->m_sName,multibuf);

						try { tempc2df = new C2DF(); } catch(...) { tempc2df = NULL; }
						if (tempc2df == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						tempc2df->CopyFrom(o->m_pCDF->m_p2DF);
						if (g_iObsChannel[0] == 0)
							if (o->m_pRDF[0]->m_bRadialCorrect)
								tempc2df->UnCorrectRadialDist(0);
						if (g_iObsChannel[0] == 1)
							if (o->m_pADF[0]->m_bStat && (!o->m_pADF[0]->m_bCosine))
								tempc2df->UnCorrectAngle(0);
						tempc2df->NormalizeXCount();
						switch(o->m_pCDF->m_iNormalize)
						{
							case 2:
								tempc2df->NormalizeBin(0.0,o->m_pCDF->m_fNormValue);
								break;
							case 1:
								tempc2df->NormalizeBinIntegral(o->m_pCDF->m_fNormValue);
								break;
							case 0:
								break;
						}
						tempc2df->WriteMathematicaNb(o->m_pCDF->m_sName,multibuf,"_nX.nb",false);
						tempc2df->WriteGnuplotInput(o->m_pCDF->m_sName,multibuf,"_nX",false);
						delete tempc2df;

						mprintf("    Saving Y-normalized CDF Mathematica Notebook as \"%s%s_nY.nb\"...\n",o->m_pCDF->m_sName,multibuf);

						try { tempc2df = new C2DF(); } catch(...) { tempc2df = NULL; }
						if (tempc2df == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						tempc2df->CopyFrom(o->m_pCDF->m_p2DF);
						if (g_iObsChannel[1] == 0)
							if (o->m_pRDF[1]->m_bRadialCorrect)
								tempc2df->UnCorrectRadialDist(1);
						if (g_iObsChannel[1] == 1)
							if (o->m_pADF[1]->m_bStat && (!o->m_pADF[1]->m_bCosine))
								tempc2df->UnCorrectAngle(1);
						tempc2df->NormalizeYCount();
						switch(o->m_pCDF->m_iNormalize)
						{
							case 2:
								tempc2df->NormalizeBin(0.0,o->m_pCDF->m_fNormValue);
								break;
							case 1:
								tempc2df->NormalizeBinIntegral(o->m_pCDF->m_fNormValue);
								break;
							case 0:
								break;
						}
						tempc2df->WriteMathematicaNb(o->m_pCDF->m_sName,multibuf,"_nY.nb",false);
						tempc2df->WriteGnuplotInput(o->m_pCDF->m_sName,multibuf,"_nY",false);
						delete tempc2df;

						mprintf("    Saving XY-normalized CDF Mathematica Notebook as \"%s%s_nXY.nb\"...\n",o->m_pCDF->m_sName,multibuf);

						try { tempc2df = new C2DF(); } catch(...) { tempc2df = NULL; }
						if (tempc2df == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
						
						tempc2df->CopyFrom(o->m_pCDF->m_p2DF);
						for (z2=0;z2<2;z2++)
						{
							if (g_iObsChannel[z2] == 0)
								if (o->m_pRDF[z2]->m_bRadialCorrect)
									tempc2df->UnCorrectRadialDist(z2);
							if (g_iObsChannel[z2] == 1)
								if (o->m_pADF[z2]->m_bStat && (!o->m_pADF[z2]->m_bCosine))
									tempc2df->UnCorrectAngle(z2);
						}
						tempc2df->NormalizeXCount();
						tempc2df->NormalizeYCount();
						switch(o->m_pCDF->m_iNormalize)
						{
							case 2:
								tempc2df->NormalizeBin(0.0,o->m_pCDF->m_fNormValue);
								break;
							case 1:
								tempc2df->NormalizeBinIntegral(o->m_pCDF->m_fNormValue);
								break;
							case 0:
								break;
						}
						tempc2df->WriteMathematicaNb(o->m_pCDF->m_sName,multibuf,"_nXY.nb",false);
						tempc2df->WriteGnuplotInput(o->m_pCDF->m_sName,multibuf,"_nXY",false);
						delete tempc2df;
					}

					mprintf("    Calculating tensor product of CDF projections...\n");

					try { o->m_pCDF->m_pTensorProduct = new C2DF(); } catch(...) { o->m_pCDF->m_pTensorProduct = NULL; }
					if (o->m_pCDF->m_pTensorProduct == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					o->m_pCDF->m_pTensorProduct->CopyFrom(o->m_pCDF->m_p2DF);
					o->m_pCDF->m_pTensorProduct->MakeTensorProduct(o->m_pCDF->m_p2DF);
					switch(o->m_pCDF->m_iNormalize)
					{
						case 2:
							o->m_pCDF->m_pTensorProduct->NormalizeBin(0.0,o->m_pCDF->m_fNormValue);
							break;
						case 1:
							o->m_pCDF->m_pTensorProduct->NormalizeBinIntegral(o->m_pCDF->m_fNormValue);
							break;
						case 0:
							break;
					}
					if (o->m_pCDF->m_bAxisDivideAll)
					{
						mprintf("    Saving tensor product as \"%s%s_tensor.dat\"...\n",o->m_pCDF->m_sName,multibuf);
						o->m_pCDF->m_pTensorProduct->Write(o->m_pCDF->m_sName,multibuf,"_tensor.dat");
						mprintf("    Saving tensor product Mathematica Notebook as \"%s%s_tensor.nb\"...\n",o->m_pCDF->m_sName,multibuf);
						o->m_pCDF->m_pTensorProduct->WriteMathematicaNb(o->m_pCDF->m_sName,multibuf,"_tensor.nb",false);
						mprintf("    Saving tensor product Gnuplot Input as \"%s%s_tensor.gp\"...\n",o->m_pCDF->m_sName,multibuf);
						o->m_pCDF->m_pTensorProduct->WriteGnuplotInput(o->m_pCDF->m_sName,multibuf,"_tensor",false);
						mprintf("    Calculating difference between CDF and tensor product...\n");
					}
					o->m_pCDF->m_p2DF->Subtract(o->m_pCDF->m_pTensorProduct);
					o->m_pCDF->m_p2DF->CalcMaxEntry();
					o->m_pCDF->m_p2DF->m_iColorScale = 4; // PlusMinus Scale
					o->m_pCDF->m_p2DF->m_fPlotExp = 1.0;
	//				o->m_pCDF->m_p2DF->m_fMathematicaColorOffset = max(fabs(o->m_pCDF->m_p2DF->m_fMinEntry),fabs(o->m_pCDF->m_p2DF->m_fMaxEntry)) / 3.0;
	//				o->m_pCDF->m_p2DF->m_fMathematicaColorScale = o->m_pCDF->m_p2DF->m_fMathematicaColorOffset  * 2.0;
					mprintf("    Saving +/- correlation plot as \"%s%s_correlation.dat\"...\n",o->m_pCDF->m_sName,multibuf);
					o->m_pCDF->m_p2DF->Write(o->m_pCDF->m_sName,multibuf,"_correlation.dat");
					mprintf("    Saving +/- correlation plot Mathematica Notebook as \"%s%s_correlation.nb\"...\n",o->m_pCDF->m_sName,multibuf);
					o->m_pCDF->m_p2DF->WriteMathematicaNb(o->m_pCDF->m_sName,multibuf,"_correlation.nb",false);
					mprintf("    Saving +/- correlation plot Gnuplot Input as \"%s%s_correlation.gp\"...\n",o->m_pCDF->m_sName,multibuf);
					o->m_pCDF->m_p2DF->WriteGnuplotInput(o->m_pCDF->m_sName,multibuf,"_correlation",false);
				} // END IF AXISDIVIDE
				if (o->m_pCDF->m_bDumpDat)
				{
					fclose(o->m_pCDF->m_fDump);
					mprintf("    CDF raw data saved as cdfdump_%dd_%s%s.csv.\n",g_iCDFChannels,o->m_pCDF->m_sShortName,multibuf);
				}
				if (o->m_bTimeDev)
				{
					if ((o->m_waSaveRefList.GetSize() == 1) || (!o->m_bSaveSeparateFiles))
					{
						fclose(o->m_pCDF->m_fTimeDev[0]);
//						sprintf(buf,"cdf_timedev_%dd_%s_ref%d%s.csv",g_iCDFChannels,o->m_pCDF->m_sShortName,z2+1,multibuf);
						mprintf("    Temporal development saved as cdf_timedev_%dd%s%s.csv\n",g_iCDFChannels,o->m_pCDF->m_sShortName,multibuf);
					} else
					{
						for (z2=0;z2<o->m_waSaveRefList.GetSize();z2++)
						{
							fclose(o->m_pCDF->m_fTimeDev[z2]);
							mprintf("    Temporal development saved as cdf_timedev_%dd%s_ref%d%s.csv\n",g_iCDFChannels,o->m_pCDF->m_sShortName,z2+1,multibuf);
						}
					}
					if (o->m_pCDF->m_bTDAnimation)
					{
						mprintf(WHITE,"  * Time dependent animation *\n");
//						sprintf(buf,"animation_complete_%s%s.agr",o->m_pCDF->m_sName,multibuf);
						buf.sprintf("animation_complete_%s%s.agr",o->m_pCDF->m_sName,multibuf);
						mprintf("    Saving complete Plot as \"%s\"...\n",(const char*)buf);
						for (z2=0;z2<o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pCDF->m_iCombinationsEnabled;z2++)
						{
							o->m_pCDF->m_pTDAPlot->SetSetLineWidth(z2,2.0);
							o->m_pCDF->m_pTDAPlot->SetSetLineColorLong(z2,GraceColor(z2,0.0));
						}
						o->m_pCDF->m_pTDAPlot->WriteAgr(buf,false);
//						sprintf(buf,"render_%s%s",o->m_pCDF->m_sName,multibuf);
						buf.sprintf("render_%s%s",o->m_pCDF->m_sName,multibuf);
						mprintf("    Saving xmgrace render script for animation as \"%s\"...\n",(const char*)buf);
						a = OpenFileWrite("gracebatch",true);
						mfprintf(a,"PRINT TO \"output.png\"\n");
						mfprintf(a,"HARDCOPY DEVICE \"PNG\"\n");
						mfprintf(a,"PAGE SIZE %d, %d\n",o->m_pCDF->m_iTDAResX,o->m_pCDF->m_iTDAResY);
						mfprintf(a,"DEVICE \"PNG\" FONT ANTIALIASING on\n");
						mfprintf(a,"DEVICE \"PNG\" OP \"compression:9\"\n");
						mfprintf(a,"PRINT\n");
						fclose(a);
						a = OpenFileWrite(buf,true);
						if (o->m_pCDF->m_bTDATrace)
						{
							for (z2=0;z2<o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pCDF->m_iCombinationsEnabled;z2++)
							{
								o->m_pCDF->m_pTDAPlot->DuplicateSet(z2);
								o->m_pCDF->m_pTDAPlot->SetSetLineWidth(z2,2.0);
								o->m_pCDF->m_pTDAPlot->SetSetLineColorLong(z2,GraceColor(z2,0.7));
								o->m_pCDF->m_pTDAPlot->SetSetLineWidth(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pCDF->m_iCombinationsEnabled+z2,3.0);
								o->m_pCDF->m_pTDAPlot->SetSetLineColorLong(o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pCDF->m_iCombinationsEnabled+z2,GraceColor(z2,0.0));
							}
						}
						for (z2=0;z2<o->m_pCDF->m_iTDASteps;z2++)
						{
							if (o->m_pCDF->m_bTDATrace)
							{
								for (z3=0;z3<o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pCDF->m_iCombinationsEnabled;z3++)
								{
									o->m_pCDF->m_pTDAPlot->SetSetRange(z3+o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pCDF->m_iCombinationsEnabled,z2*o->m_pCDF->m_iTDAStride,z2*o->m_pCDF->m_iTDAStride+o->m_pCDF->m_iTDATail);
									o->m_pCDF->m_pTDAPlot->SetSetRange(z3,0,z2*o->m_pCDF->m_iTDAStride);
								}
							} else
							{
								for (z3=0;z3<o->m_waSaveRefList.GetSize()*o->m_waSaveShowList.GetSize()*o->m_pCDF->m_iCombinationsEnabled;z3++)
									o->m_pCDF->m_pTDAPlot->SetSetRange(z3,z2*o->m_pCDF->m_iTDAStride,z2*o->m_pCDF->m_iTDAStride+o->m_pCDF->m_iTDATail);
							}

//							sprintf(buf,"animation_%05d_%s%s.agr",z2+1,o->m_pCDF->m_sName,multibuf);
							buf.sprintf("animation_%05d_%s%s.agr",z2+1,o->m_pCDF->m_sName,multibuf);
							mprintf("    Saving frame \"%s\"...\n",(const char*)buf);
							o->m_pCDF->m_pTDAPlot->WriteAgr(buf,false);
							mfprintf(a,"echo 'Printing frame %d of %d'\n",z2+1,o->m_pCDF->m_iTDASteps);
							mfprintf(a,"xmgrace %s -batch gracebatch -nosafe -hardcopy\n",(const char*)buf);
							mfprintf(a,"mv output.png frame%04d.png\n",z2+1);
						}
						fclose(a);
					}
				}
			} // END IF CHANNELS == 2

			if (g_iCDFChannels == 3)
			{
				mprintf(WHITE,"\n* Combined Distribution Function (3D)\n");
				mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pCDF->m_p3DF->m_fBinEntries,o->m_pCDF->m_p3DF->m_fSkipEntries,ZeroDivide(o->m_pCDF->m_p3DF->m_fSkipEntries,o->m_pCDF->m_p3DF->m_fSkipEntries+o->m_pCDF->m_p3DF->m_fBinEntries)*100.0);
				if (o->m_pCDF->m_p3DF->m_fBinEntries == 0)
				{
					eprintf("    There were no bin entries. Check your CDF definition. Skipping this CDF.\n\n");
					goto _skipsdf;
				}
				o->m_pCDF->m_p3DF->CalcMaxEntry();
				mprintf("    Data range from %G to %G.\n",o->m_pCDF->m_p3DF->m_fMinEntry,o->m_pCDF->m_p3DF->m_fMaxEntry);
				mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pCDF->m_p3DF->m_fMaxEntry,o->m_pCDF->m_p3DF->m_fEps);
				if (o->m_pCDF->m_p3DF->m_fEps > 1.0E-4) {
					eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
					mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
				}
				for (z2=0;z2<g_iCDFChannels;z2++)
				{
					if (g_iObsChannel[z2] == 0)
					{
						if (o->m_pRDF[z2]->m_bRadialCorrect)
						{
							mprintf("    Correcting radial distribution for CDF channel %d...\n",z2+1);
							o->m_pCDF->m_p3DF->CorrectRadialDist(z2);
							o->m_pCDF->m_p3DF->CalcMaxEntry();
							mprintf("    Data range from %G to %G.\n",o->m_pCDF->m_p3DF->m_fMinEntry,o->m_pCDF->m_p3DF->m_fMaxEntry);
						}
					}
					if (g_iObsChannel[z2] == 1)
					{
						if (o->m_pADF[z2]->m_bStat && (!o->m_pADF[z2]->m_bCosine))
						{
							mprintf("    Correcting angular distribution for CDF channel %d...\n",z2+1);
							o->m_pCDF->m_p3DF->CorrectAngle(z2);
							o->m_pCDF->m_p3DF->CalcMaxEntry();
							mprintf("    Data range from %G to %G.\n",o->m_pCDF->m_p3DF->m_fMinEntry,o->m_pCDF->m_p3DF->m_fMaxEntry);
						}
						if (o->m_pADF[z2]->m_bMirror)
						{
							mprintf("    Making channel %d mirror-symmetric...\n",z2+1);
							if (o->m_pADF[z2]->m_bCosine)
								o->m_pCDF->m_p2DF->Mirror(0.0,z2);
									else o->m_pCDF->m_p2DF->Mirror(90.0,z2);
							o->m_pCDF->m_p3DF->CalcMaxEntry();
							mprintf("    Data range from %G to %G.\n",o->m_pCDF->m_p3DF->m_fMinEntry,o->m_pCDF->m_p3DF->m_fMaxEntry);
						}
					}
				}
//				o->m_pCDF->m_p3DF->CalcMaxEntry();
//				mprintf("    Data range now from %G to %G.\n",o->m_pCDF->m_p3DF->m_fMinEntry,o->m_pCDF->m_p3DF->m_fMaxEntry);
				o->m_pCDF->m_fFactor = 1.0;
				switch(o->m_pCDF->m_iNormalize)
				{
					case 2:
						mprintf("    Normalizing maximum value to %.2f.\n",o->m_pCDF->m_fNormValue);
						o->m_pCDF->m_p3DF->NormalizeBin(0.0,o->m_pCDF->m_fNormValue);
						o->m_pCDF->m_p3DF->CalcMaxEntry();
						mprintf("    Data range now from %G to %G.\n",o->m_pCDF->m_p3DF->m_fMinEntry,o->m_pCDF->m_p3DF->m_fMaxEntry);
						break;
					case 1:
						mprintf("    Normalizing integral value to %.2f.\n",o->m_pCDF->m_fNormValue);
						o->m_pCDF->m_fFactor = o->m_pCDF->m_p3DF->NormalizeBinIntegral(o->m_pCDF->m_fNormValue);
						mprintf("    Factor was %G.\n",o->m_pCDF->m_fFactor);
						o->m_pCDF->m_p3DF->CalcMaxEntry();
						mprintf("    Data range now from %G to %G.\n",o->m_pCDF->m_p3DF->m_fMinEntry,o->m_pCDF->m_p3DF->m_fMaxEntry);
						break;
					case 0:
						mprintf("    Not normalizing this CDF.\n");
						break;
				}

				for (z2=0;z2<=g_iSDFSmoothGrade;z2++)
				{
					try { temp3DF = new C3DF<double>(); } catch(...) { temp3DF = NULL; }
					if (temp3DF == NULL) NewException((double)sizeof(C3DF<double>),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					temp3DF->CopyFrom(o->m_pCDF->m_p3DF);
					if (z2 != 0)
					{
						temp3DF->Smooth(z2);
//						sprintf(buf,".s%d%s.plt",z2,multibuf);
						buf.sprintf(".s%d%s.plt",z2,multibuf);
					} else
//						sprintf(buf,"%s.plt",multibuf);
						buf.sprintf("%s.plt",multibuf);

					mprintf("    Saving 3D CDF as \"%s%s\"...\n",o->m_pCDF->m_sName,(const char*)buf);
					temp3DF->WritePLT("",o->m_pCDF->m_sName,buf,false);

					if (z2 != 0)
//						sprintf(buf,".s%d%s.cube",z2,multibuf);
						buf.sprintf(".s%d%s.cube",z2,multibuf);
					else
//						sprintf(buf,"%s.cube",multibuf);
						buf.sprintf("%s.cube",multibuf);

					mprintf("    Saving 3D CDF as \"%s%s\"...\n",o->m_pCDF->m_sName,(const char*)buf);
					temp3DF->WriteCube("",o->m_pCDF->m_sName,buf,false);
					if (o->m_pCDF->m_iHistogramRes != 0)
					{
						mprintf("    Calculating Histogram...\n");
						temp3DF->CalcHistogram();
						if (z2 != 0)
//							sprintf(buf,"his_%s.s%d",o->m_pCDF->m_sName,z2);
							buf.sprintf("his_%s.s%d",o->m_pCDF->m_sName,z2);
						else
//							sprintf(buf,"his_%s",o->m_pCDF->m_sName);
							buf.sprintf("his_%s",o->m_pCDF->m_sName);

						mprintf("    Saving Histogram as \"%s%s.csv\"...\n",(const char*)buf,multibuf);
						temp3DF->WriteHistogram(buf,multibuf,".csv");
					}

					if (o->m_pCDF->m_b3DSlices)
					{
						for (z3=0;z3<3;z3++)
						{
							if (o->m_pCDF->m_i3DSliceIntervals[z3] != 0)
							{
								if (z2 != 0)
									mprintf(WHITE,"    * Creating channel %d smooth degree %d slices...\n",z3+1,z2);
								else
									mprintf(WHITE,"    * Creating channel %d slices...\n",z3+1);

								try { temp2dfa = new C2DF*[o->m_pCDF->m_i3DSliceIntervals[z3]]; } catch(...) { temp2dfa = NULL; }
								if (temp2dfa == NULL) NewException((double)o->m_pCDF->m_i3DSliceIntervals[z3]*sizeof(C2DF*),__FILE__,__LINE__,__PRETTY_FUNCTION__);

								tf = 0;
								for (z4=0;z4<o->m_pCDF->m_i3DSliceIntervals[z3];z4++)
								{
									try { temp2dfa[z4] = new C2DF(); } catch(...) { temp2dfa[z4] = NULL; }
									if (temp2dfa[z4] == NULL) NewException((double)sizeof(C2DF),__FILE__,__LINE__,__PRETTY_FUNCTION__);
									
									temp3DF->CreateSlice(z3,(int)((double)z4/o->m_pCDF->m_i3DSliceIntervals[z3]*o->m_pCDF->m_p3DF->m_iRes[z3]),(int)((z4+1.0)/o->m_pCDF->m_i3DSliceIntervals[z3]*o->m_pCDF->m_p3DF->m_iRes[z3]-1),temp2dfa[z4]);
									mprintf("        - %2d: %9.4f - %9.4f (%3d - %3d), Intensity range %.3f ... %.3f\n",z4+1,o->m_pCDF->m_p3DF->m_fMinVal[z3]+(double)z4/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),o->m_pCDF->m_p3DF->m_fMinVal[z3]+(z4+1.0)/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),int((double)z4/o->m_pCDF->m_i3DSliceIntervals[z3]*o->m_pCDF->m_p3DF->m_iRes[z3]),int((z4+1.0)/o->m_pCDF->m_i3DSliceIntervals[z3]*o->m_pCDF->m_p3DF->m_iRes[z3]-1),temp2dfa[z4]->m_fMinEntry,temp2dfa[z4]->m_fMaxEntry);
									if (temp2dfa[z4]->m_fMaxEntry > tf)
										tf = temp2dfa[z4]->m_fMaxEntry;
								}

								mprintf("      Setting common intensity range to 0 ... %.3f.\n",tf);

								mprintf("      Writing out slice files (%d/3)...\n",z3+1);

								for (z4=0;z4<o->m_pCDF->m_i3DSliceIntervals[z3];z4++)
								{
									temp2dfa[z4]->m_fMinEntry = 0;
									temp2dfa[z4]->m_fMaxEntry = tf;

									if (z2 != 0)
//										sprintf(buf,"%s_slice_ch%d_%d_%.3f-%.3f.s%d.nb",o->m_pCDF->m_sName,z3+1,z4+1,o->m_pCDF->m_p3DF->m_fMinVal[z3]+(double)z4/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),o->m_pCDF->m_p3DF->m_fMinVal[z3]+(z4+1.0)/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),z2);
										buf.sprintf("%s_slice_ch%d_%d_%.3f-%.3f.s%d.nb",o->m_pCDF->m_sName,z3+1,z4+1,o->m_pCDF->m_p3DF->m_fMinVal[z3]+(double)z4/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),o->m_pCDF->m_p3DF->m_fMinVal[z3]+(z4+1.0)/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),z2);
									else
//										sprintf(buf,"%s_slice_ch%d_%d_%.3f-%.3f.nb",o->m_pCDF->m_sName,z3+1,z4+1,o->m_pCDF->m_p3DF->m_fMinVal[z3]+(double)z4/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),o->m_pCDF->m_p3DF->m_fMinVal[z3]+(z4+1.0)/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]));
										buf.sprintf("%s_slice_ch%d_%d_%.3f-%.3f.nb",o->m_pCDF->m_sName,z3+1,z4+1,o->m_pCDF->m_p3DF->m_fMinVal[z3]+(double)z4/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),o->m_pCDF->m_p3DF->m_fMinVal[z3]+(z4+1.0)/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]));

									mprintf("        %s ...\n",(const char*)buf);
									temp2dfa[z4]->WriteMathematicaNb("",buf,"",true);

									if (z2 != 0)
//										sprintf(buf,"%s_slice_ch%d_%d_%.3f-%.3f.s%d",o->m_pCDF->m_sName,z3+1,z4+1,o->m_pCDF->m_p3DF->m_fMinVal[z3]+(double)z4/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),o->m_pCDF->m_p3DF->m_fMinVal[z3]+(z4+1.0)/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),z2);
										buf.sprintf("%s_slice_ch%d_%d_%.3f-%.3f.s%d",o->m_pCDF->m_sName,z3+1,z4+1,o->m_pCDF->m_p3DF->m_fMinVal[z3]+(double)z4/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),o->m_pCDF->m_p3DF->m_fMinVal[z3]+(z4+1.0)/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),z2);
									else
//										sprintf(buf,"%s_slice_ch%d_%d_%.3f-%.3f",o->m_pCDF->m_sName,z3+1,z4+1,o->m_pCDF->m_p3DF->m_fMinVal[z3]+(double)z4/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),o->m_pCDF->m_p3DF->m_fMinVal[z3]+(z4+1.0)/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]));
										buf.sprintf("%s_slice_ch%d_%d_%.3f-%.3f",o->m_pCDF->m_sName,z3+1,z4+1,o->m_pCDF->m_p3DF->m_fMinVal[z3]+(double)z4/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]),o->m_pCDF->m_p3DF->m_fMinVal[z3]+(z4+1.0)/o->m_pCDF->m_i3DSliceIntervals[z3]*(o->m_pCDF->m_p3DF->m_fMaxVal[z3]-o->m_pCDF->m_p3DF->m_fMinVal[z3]));

									mprintf("        %s.gp ...\n",(const char*)buf);
									temp2dfa[z4]->WriteGnuplotInput("",buf,"",true);
								}

								for (z4=0;z4<o->m_pCDF->m_i3DSliceIntervals[z3];z4++)
									delete temp2dfa[z4];

								delete[] temp2dfa;
							}
						}
					}
				}

				mprintf("    Writing out 2D projections...\n");
				for (z3=0;z3<3;z3++)
				{
					mprintf(WHITE,"    * CDF 2D projection on channel %d\n",z3+1);
					mprintf("        %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fBinEntries,o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fSkipEntries,ZeroDivide(o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fSkipEntries,o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fSkipEntries+o->m_pCDF->m_p3DF->m_p2DF[z3]->m_fBinEntries)*100.0);
					switch(z3)
					{
						case 0:
							tia[0] = 0;
							tia[1] = 1;
							break;
						case 1:
							tia[0] = 0;
							tia[1] = 2;
							break;
						case 2:
							tia[0] = 1;
							tia[1] = 2;
							break;
					}
					for (z2=0;z2<2;z2++)
					{
						if (g_iObsChannel[tia[z2]] == 0)
						{
							if (o->m_pRDF[tia[z2]]->m_bRadialCorrect)
							{
								mprintf("        Correcting radial distribution for CDF channel %d...\n",z2+1);
								o->m_pCDF->m_p3DF->m_p2DF[z3]->CorrectRadialDist(z2);
							}
						}
						if (g_iObsChannel[tia[z2]] == 1)
						{
							if (o->m_pADF[tia[z2]]->m_bStat && (!o->m_pADF[tia[z2]]->m_bCosine))
							{
								mprintf("        Correcting angular distribution for CDF channel %d...\n",z2+1);
								o->m_pCDF->m_p3DF->m_p2DF[z3]->CorrectAngle(z2);
							}
		/*					if (o->m_pADF[z2]->m_bMirror)
							{
								mprintf("        Making channel %d mirror-symmetric...\n",z2+1);
								if (o->m_pADF[z2]->m_bCosine)
									o->m_pCDF->m_p2DF->Mirror(0.0f,z2);
										else o->m_pCDF->m_p2DF->Mirror(90.0f,z2);
							}*/
						}
					}
					o->m_pCDF->m_fFactor = 1.0;
					switch(o->m_pCDF->m_iNormalize)
					{
						case 2:
							mprintf("        Normalizing maximum value to %.2f.\n",o->m_pCDF->m_fNormValue);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->NormalizeBin(0.0,o->m_pCDF->m_fNormValue);
							break;
						case 1:
							mprintf("        Normalizing integral value to %.2f.\n",o->m_pCDF->m_fNormValue);
							o->m_pCDF->m_p3DF->m_p2DF[z3]->NormalizeBinIntegral(o->m_pCDF->m_fNormValue);
							break;
						case 0:
							mprintf("        Not normalizing this CDF.\n");
							break;
					}
					mprintf("        Saving CDF triples as \"%s%s_triples.csv\"...\n",o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName,multibuf);
					o->m_pCDF->m_p3DF->m_p2DF[z3]->Write(o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName,multibuf,"_triples.csv");
					mprintf("        Saving CDF matrix as \"%s%s_matrix.csv\"...\n",o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName,multibuf);
					o->m_pCDF->m_p3DF->m_p2DF[z3]->WriteCSV(o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName,multibuf,"_matrix.csv");
					mprintf("        Saving CDF Mathematica Notebook as \"%s%s.nb\"...\n",o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName,multibuf);
					o->m_pCDF->m_p3DF->m_p2DF[z3]->WriteMathematicaNb(o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName,multibuf,".nb",false);
					mprintf("        Saving CDF Gnuplot Input as \"%s%s.gp\"...\n",o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName,multibuf);
					o->m_pCDF->m_p3DF->m_p2DF[z3]->WriteGnuplotInput(o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName,multibuf,"",false);
					mprintf("        Saving combined plot AGR file as \"%s%s_combined.agr\"...\n",o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName,multibuf);
					o->m_pCDF->m_p3DF->m_p2DF[z3]->WriteCombinedPlot(o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName,multibuf,"_combined.agr");
					if (o->m_pCDF->m_iHistogramRes != 0)
					{
						mprintf("    Calculating Histogram...\n");
						o->m_pCDF->m_p3DF->m_p2DF[z3]->CalcHistogram();
//						sprintf(buf,"his_%s",o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName);
						buf.sprintf("his_%s",o->m_pCDF->m_p3DF->m_p2DF[z3]->m_sName);
						mprintf("    Saving Histogram as \"%s%s.csv\"...\n",(const char*)buf,multibuf);
						o->m_pCDF->m_p3DF->m_p2DF[z3]->WriteHistogram(buf,multibuf,".csv");
					}
				}
			} // END IF CHANNELS == 3
			if (o->m_pCDF->m_bWriteSnapshots) {
				mprintf(WHITE,"    * Snapshot Analysis\n");
				mprintf("    %ld configurations in total fulfilled the criteria.\n",o->m_pCDF->m_iWriteSnapshotsTotalCounter);
				mprintf("      (on average %.8f per RM and time step - multiple configurations per RM are possible)\n",
					(double)o->m_pCDF->m_iWriteSnapshotsTotalCounter/g_iSteps/((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()
				);
				mprintf("    The probability that any RM has at least one matching configuration is %.6f %%.\n",
					(double)o->m_pCDF->m_iWriteSnapshotsTotalPerRM/g_iSteps/((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize()*100.0
				);
				if (!o->m_pCDF->m_bWriteSnapshotsNoWrite) {
					if (o->m_pCDF->m_pWriteSnapshotsFile != NULL) {
						mprintf("    Closing CDF snapshot file \"%s\" ...\n",(const char*)o->m_pCDF->m_sWriteSnapshotsFileName);
						fclose( o->m_pCDF->m_pWriteSnapshotsFile );
						o->m_pCDF->m_pWriteSnapshotsFile = NULL;
						mprintf("    %ld snapshots have been written to file.\n",o->m_pCDF->m_iWriteSnapshotsTotalWritten);
					} else
						eprintf("    Warning: CDF snapshot file was not open.\n");
				} else
					mprintf("    Only counting, no snapshot file was written.\n");
			}
_skipsdf:;
		} // END IF CDF

		if (g_bRevSDF)
		{
			mprintf(WHITE,"\n* Pseudo SDF\n");
			mprintf("    %.0f bin entries, %.0f out of bin range (%.2f percent).\n",o->m_pRevSDF->m_p2DF->m_fBinEntries,o->m_pRevSDF->m_p2DF->m_fSkipEntries,ZeroDivide(o->m_pRevSDF->m_p2DF->m_fSkipEntries,o->m_pRevSDF->m_p2DF->m_fSkipEntries+o->m_pRevSDF->m_p2DF->m_fBinEntries)*100.0);
			if (o->m_pRevSDF->m_p2DF->m_fBinEntries == 0)
			{
				eprintf("    There were no bin entries. Check your function definition. Skipping this PseudoSDF.\n\n");
				goto _skiprevsdf;
			}
			o->m_pRevSDF->m_p2DF->CalcMaxEntry();
			mprintf("    Max. bin entry: %.6E  -->  EPS: %.6E\n",o->m_pRevSDF->m_p2DF->m_fMaxEntry,o->m_pRevSDF->m_p2DF->m_fEps);
			if (o->m_pRevSDF->m_p2DF->m_fEps > 1.0E-4) {
				eprintf("\n    Warning: Very large bin entries - probably loss of accuracy occurred.\n");
				mprintf("    Please reduce the bin counts (e.g. by analyzing only every 10th step).\n\n");
			}
			mprintf("    Normalizing integral value to %.2f.\n",1000000.0);
			o->m_pRevSDF->m_p2DF->NormalizeBinIntegral(1000000.0);
			if (o->m_pRevSDF->m_bDrawAtoms)
			{
				mprintf("    Inserting reference atoms into the Pseudo SDF plot...\n");

				if (((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_pElement->m_fRadius == 0)
					o->m_pRevSDF->m_p2DF->AddCircle(0,0,50.0,0.4,0.4,0.4);
				else
					o->m_pRevSDF->m_p2DF->AddCircle(0,0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_pElement->m_fRadius*0.75,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_pElement->m_iColorR/255.0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_pElement->m_iColorG/255.0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_pElement->m_iColorB/255.0);

				if (((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_pElement->m_fRadius == 0)
					o->m_pRevSDF->m_p2DF->AddCircle(0,o->m_pRevSDF->m_fSecondAtomPosX/o->m_pRevSDF->m_fSecondAtomCount,50.0,0.4,0.4,0.4);
				else
					o->m_pRevSDF->m_p2DF->AddCircle(0,o->m_pRevSDF->m_fSecondAtomPosX/o->m_pRevSDF->m_fSecondAtomCount,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_pElement->m_fRadius*0.75,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_pElement->m_iColorR/255.0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_pElement->m_iColorG/255.0,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_pElement->m_iColorB/255.0);
			}
			mprintf("    Saving PseudoSDF triples as \"%s%s_triples.csv\"...\n",o->m_pRevSDF->m_sName,multibuf);
			o->m_pRevSDF->m_p2DF->Write(o->m_pRevSDF->m_sName,multibuf,"_triples.csv");
			mprintf("    Saving PseudoSDF matrix as \"%s%s_matrix.csv\"...\n",o->m_pRevSDF->m_sName,multibuf);
			o->m_pRevSDF->m_p2DF->WriteCSV(o->m_pRevSDF->m_sName,multibuf,"_matrix.csv");
			mprintf("    Saving PseudoSDF Mathematica Notebook as \"%s%s.nb\"...\n",o->m_pRevSDF->m_sName,multibuf);
			o->m_pRevSDF->m_p2DF->WriteMathematicaNb(o->m_pRevSDF->m_sName,multibuf,".nb",false);
			mprintf("    Saving PseudoSDF Gnuplot Input as \"%s%s.gp\"...\n",o->m_pRevSDF->m_sName,multibuf);
			o->m_pRevSDF->m_p2DF->WriteGnuplotInput(o->m_pRevSDF->m_sName,multibuf,"",false);
			if (o->m_pRevSDF->m_iHistogramRes != 0)
			{
				mprintf("    Calculating Histogram...\n");
				o->m_pRevSDF->m_p2DF->CalcHistogram();
//				sprintf(buf,"his_%s",o->m_pRevSDF->m_sName);
				buf.sprintf("his_%s",o->m_pRevSDF->m_sName);
				mprintf("    Saving Histogram as \"%s%s.csv\"...\n",(const char*)buf,multibuf);
				o->m_pRevSDF->m_p2DF->WriteHistogram(buf,multibuf,".csv");
			}
			if (o->m_pRevSDF->m_bCreateRevSDF)
			{
				mprintf("    Creating volumetric Revolution SDF...\n");
				o->m_pRevSDF->CreateRevSDF();
//				sprintf(buf,"revsdf_ref_%s.xyz",o->m_pRevSDF->m_sName);
				buf.sprintf("revsdf_ref_%s.xyz",o->m_pRevSDF->m_sName);
				mprintf("    Saving Revolution SDF reference structure as \"%s\"...\n",(const char*)buf);
				a = OpenFileWrite(buf,true);
				mfprintf(a,"2\n\n%s   0.0  0.0  0.0\n%s  0.0  %f  0.0\n",(const char*)((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,(const char*)((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,o->m_pRevSDF->m_fSecondAtomPosX/o->m_pRevSDF->m_fSecondAtomCount/100.0);
				fclose(a);
				for (z2=0;z2<=g_iSDFSmoothGrade;z2++)
				{
					try { tempSDF = new C3DF<double>(); } catch(...) { tempSDF = NULL; }
					if (tempSDF == NULL) NewException((double)sizeof(C3DF<double>),__FILE__,__LINE__,__PRETTY_FUNCTION__);
					
					tempSDF->CopyFrom(o->m_pRevSDF->m_pRevSDF);
					if (z2 != 0)
					{
						tempSDF->Smooth(z2);
//						sprintf(buf,".s%d%s.plt",z2,multibuf);
						buf.sprintf(".s%d%s.plt",z2,multibuf);
					} else
//						sprintf(buf,"%s.plt",multibuf);
						buf.sprintf("%s.plt",multibuf);

					mprintf("    Saving Revolution SDF as \"revsdf_%s%s\"...\n",o->m_pRevSDF->m_sName,(const char*)buf);
					tempSDF->WritePLT("revsdf_",o->m_pRevSDF->m_sName,buf,true);

					if (z2 != 0)
//						sprintf(buf,".s%d%s.cube",z2,multibuf);
						buf.sprintf(".s%d%s.cube",z2,multibuf);
					else
//						sprintf(buf,"%s.cube",multibuf);
						buf.sprintf("%s.cube",multibuf);

					mprintf("    Saving Revolution SDF as \"revsdf_%s%s\"...\n",o->m_pRevSDF->m_sName,(const char*)buf);
					tempSDF->WriteCube("revsdf_",o->m_pRevSDF->m_sName,buf,true);

					if (z2 != 0)
//						sprintf(buf,".s%d%s.his",z2,multibuf);
						buf.sprintf(".s%d%s.his",z2,multibuf);
					else
//						sprintf(buf,"%s.his",multibuf);
						buf.sprintf("%s.his",multibuf);

					mprintf("    Saving Revolution SDF Histogram as \"revsdf_%s%s\"...\n",o->m_pRevSDF->m_sName,(const char*)buf);
					tempSDF->CalcHistogram();
					tempSDF->WriteHistogram("revsdf_",o->m_pRevSDF->m_sName,buf);
					if (o->m_pRevSDF->m_iHistogramRes != 0)
					{
						mprintf("    Calculating Histogram...\n");
						tempSDF->CalcHistogram();
//						sprintf(buf2,"his_revsdf_%s",o->m_pRevSDF->m_sName);
						buf2.sprintf("his_revsdf_%s",o->m_pRevSDF->m_sName);
						mprintf("    Saving Histogram as \"%s%s.csv\"...\n",(const char*)buf2,(const char*)buf);
						tempSDF->WriteHistogram(buf2,buf,".csv");
					}
				}
			}
_skiprevsdf:;
		} // END IF REVSDF
	}


	if (g_bSDF && g_bSDFMap)
	{
		mprintf(WHITE,"    ### Finishing SDF surface maps ###\n\n");
		for (z=0;z<g_oaSDFMaps.GetSize();z++)
			((CSDFMap*)g_oaSDFMaps[z])->Finish();
		mprintf("\n");
	}

	if (g_bSDF || g_bAvg
      )
	{
//		sprintf(buf,"ref_%s_%s%d_%s%d_%s%d%s.xyz",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,multibuf);
		buf.sprintf("ref_%s_%s%d_%s%d_%s%d%s.xyz",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName,(const char*)((CAtom*)g_oaAtoms[g_iFixRealAtomType[0]])->m_sName,g_iFixAtom[0]+1,(const char*)((CAtom*)g_oaAtoms[g_iFixRealAtomType[1]])->m_sName,g_iFixAtom[1]+1,(const char*)((CAtom*)g_oaAtoms[g_iFixRealAtomType[2]])->m_sName,g_iFixAtom[2]+1,multibuf);
//		strtolower(buf);
		buf.ToLowercase();
		mprintf(WHITE,"\n*** Reference molecule ***\n\n");
		FreeFileName(buf);
		mprintf("    Saving reference molecule as \"%s\" ...  ",(const char*)buf);
		a = OpenFileWrite(buf,true);
		if (a == NULL)
		{
			eprintf("\nError: Could not open %s for writing.\n",(const char*)buf);
			return 0;
		}

		if (!g_bSaveVirtAtoms)
			mfprintf(a,"  %d \n\n",(int)(((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes-((CMolecule*)g_oaMolecules[g_iFixMol])->m_laVirtualAtoms.GetSize()));
		else
			mfprintf(a,"  %d \n\n",((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes);

		cc = 0;
		for (z=0;z<((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex.GetSize();z++)
		{
			if ((!g_bSaveVirtAtoms) && (((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[z] == g_iVirtAtomType))
				continue;
			for (z2=0;z2<((CMolecule*)g_oaMolecules[g_iFixMol])->m_waAtomCount[z];z2++)
			{
				if (g_bMiddleAvg)
				{
					g_pRefMol[cc][0] /= ((double)g_iSteps * ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize());
					g_pRefMol[cc][1] /= ((double)g_iSteps * ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize());
					g_pRefMol[cc][2] /= ((double)g_iSteps * ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize());
				}

				if (((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[z] == g_iVirtAtomType)
					mfprintf(a,"  %s  %9.6f  %9.6f  %9.6f\n",((CVirtualAtom*)g_oaVirtualAtoms[((CMolecule*)g_oaMolecules[g_iFixMol])->m_laVirtualAtoms[z2]])->m_sLabel,g_pRefMol[cc][0]/100.0,g_pRefMol[cc][1]/100.0,g_pRefMol[cc][2]/100.0);
				else
					mfprintf(a,"  %s  %9.6f  %9.6f  %9.6f\n",(const char*)((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[z]])->m_sName,g_pRefMol[cc][0]/100.0,g_pRefMol[cc][1]/100.0,g_pRefMol[cc][2]/100.0);

				cc++;
			}
		}
		fclose(a);
		mprintf("Done.\n");
		if (g_iSwapAtoms)
		{
			mprintf("\n*** This is the swap matrix of reference molecule %s:\n  (It shows how often the atoms\n  have been swapped with each other)\n\n     ",((CMolecule*)g_oaMolecules[g_iFixMol])->m_sName);
			for (z2=0;z2<((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex.GetSize();z2++)
			{
				for (z3=0;z3<((CMolecule*)g_oaMolecules[g_iFixMol])->m_waAtomCount[z2];z3++)
				{
					if (((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[z2]])->m_sName[1] == 0)
						mprintf(" ");
					mprintf("%s%d ",(const char*)((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[z2]])->m_sName,z3+1);	
					if (z3 < 10)
						mprintf(" ");
				}
			}
			mprintf("\n");
			cc = 0;
			for (z2=0;z2<((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex.GetSize();z2++)
			{
				for (z3=0;z3<((CMolecule*)g_oaMolecules[g_iFixMol])->m_waAtomCount[z2];z3++)
				{
					mprintf("%s%d ",(const char*)((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[z2]])->m_sName,z3+1);
					if (((CAtom*)g_oaAtoms[((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex[z2]])->m_sName[1] == 0)
						mprintf(" ");
					if (z3 < 9)
						mprintf(" ");
					cc2 = 0;
					for (z4=0;z4<((CMolecule*)g_oaMolecules[g_iFixMol])->m_baAtomIndex.GetSize();z4++)
					{
						for (z5=0;z5<((CMolecule*)g_oaMolecules[g_iFixMol])->m_waAtomCount[z4];z5++)
						{
							if ((z2 == z4) && (z3 == z5))
							{
								mprintf(" *** ");
								goto _swm;
							}
							if (pSwapMatrix[cc*((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes+cc2] != 0)
								mprintf("%3.0f%c ",(double)pSwapMatrix[cc*((CMolecule*)g_oaMolecules[g_iFixMol])->m_iAtomGes+cc2] * 100.0 / g_iSteps / ((CMolecule*)g_oaMolecules[g_iFixMol])->m_laSingleMolIndex.GetSize(),'%');
									else mprintf("  -  ");
_swm:
							cc2++;
						}
					}	
					mprintf("\n");
					cc++;
				}
			}
		}
	}

_ende2:
	if (g_bMultiInterval)
	{
		mprintf("\n");
		mprintf(YELLOW,"***********************************************\n");
		mprintf(YELLOW,"*** Interval %2d (%6d - %6d) finished. ***\n",multicounter+1,g_laMultiIntervalStart[multicounter]+1,g_laMultiIntervalEnd[multicounter]+1);
		mprintf(YELLOW,"***********************************************\n\n");
		multicounter++;
		if (multicounter < g_laMultiIntervalStart.GetSize())
			goto _multiintervalstart;
	}

	mprintf(WHITE,"\n\n----------------------------------------------------------------------------------------------------------\n");
	mprintf(WHITE,"    Finished all analyses.\n\n");

	g_iEndTime = time(NULL);
	if (((g_iEndTime-g_iStartTime) > 5) && (g_iStartTime != 0))
	{
		mprintf(WHITE,"\n    Analysis was running for %luh %lum %lus.\n",(unsigned long)((g_iEndTime-g_iStartTime)/3600),(unsigned long)(((g_iEndTime-g_iStartTime)/60)%60),(unsigned long)((g_iEndTime-g_iStartTime)%60));
		mprintf(WHITE,"      (%.3f steps per second)\n\n",g_iSteps/g_iStride/((double)g_iEndTime-g_iStartTime));
	} //else mprintf("\n");

_ende:
	WriteCredits();

	if (g_bSMode)
	{
		PrintSMode();
	} else
	{
		mprintf(WHITE,"*** The End ***\n\n");
	}
	UninstallSignalHandler();

	if (g_bROA) {
		delete g_pROAEngine;
		g_pROAEngine = NULL;
	}


	if (g_fInput != NULL) {
		fclose(g_fInput);
		g_fInput = NULL;
	}

	if (g_pDatabase != NULL)
	{
		delete g_pDatabase;
		g_pDatabase = NULL;
	}

	if (g_pTempTimestep != NULL)
	{
		delete g_pTempTimestep;
		g_pTempTimestep = NULL;
	}


	if (apfa != NULL)
	{
		delete[] apfa;
		apfa = NULL;
	}

	if (tda != NULL)
	{
		delete[] tda;
		tda = NULL;
	}

	if (g_pBQBFile != NULL) {
		delete g_pBQBFile;
		g_pBQBFile = NULL;
	}

	if (g_sExeName != NULL) {
		delete[] g_sExeName;
		g_sExeName = NULL;
	}

	if (g_sHostName != NULL) {
		delete[] g_sHostName;
		g_sHostName = NULL;
	}

	if (g_pBQBInterface != NULL) {
		BQBReleaseInterface( g_pBQBInterface );
		g_pBQBInterface = NULL;
	}

	if (g_sWorkingDir != NULL) {
		free( g_sWorkingDir );
		g_sWorkingDir = NULL;
	}

	RemoveAllAnalyses();
	RemoveAllElements();
	RemoveAllAtoms();
	RemoveAllMolecules();
	RemoveAllObservations();

	if (g_pLogFile != NULL) {
		fclose( g_pLogFile );
		g_pLogFile = NULL;
	}

	BTOUT;
	return 0;
}




