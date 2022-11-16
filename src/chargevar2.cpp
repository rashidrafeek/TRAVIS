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


// This must always be the first include directive
#include "config.h"

#include "chargevar2.h"
#include "tools.h"
#include "globalvar.h"
#include "maintools.h"
#include "element.h"
#include "conversion.h"


const char *GetRevisionInfo_chargevar2(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_chargevar2() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}



#ifdef NEW_CHARGEVAR



CChargeVar2::CChargeVar2() {

	m_pBQBFile = NULL;
	m_iaRefining = { 1, 1, 1 };
	m_bVerboseTiming = false;
}



CChargeVar2::~CChargeVar2() {
	int z;
	for (z=0;z<(int)m_oaSnapshots.size();z++) {
		if (m_oaSnapshots[z] != NULL) {
			delete m_oaSnapshots[z];
			m_oaSnapshots[z] = NULL;
		}
	}
	if (m_pBQBFile != NULL) {
		delete m_pBQBFile;
		m_pBQBFile = NULL;
	}
}



bool CChargeVar2::loadData( const char *filename, int steps, int stride ) {

	const char *p;
	int i, n, z, z2, z3;
	CxDVector3 veca, vecb, vecc;
	double tfe, tfn, fac, tfac=0;
	bool isbqb;
	CSingleMolecule *sm;
	CChargeVar2_Snapshot *snap;
	CTimeStep ts;
	FILE *a=NULL;


	p = strrchr( filename, '.' );
	if (p == NULL) {
		eprintf("CChargeVar2::loadData(): Error: Could not determine file type of trajectory \"%s\".\n",filename);
		return false;
	}

	if (mystricmp(p,".cube") == 0) {

		isbqb = false;

	} else if ((mystricmp(p,".bqb") == 0) || (mystricmp(p,".bbq") == 0) || (mystricmp(p,".blist") == 0)) {

		isbqb = true;

	} else {

		eprintf("CChargeVar2::loadData(): Error: Could not determine file type of trajectory \"%s\".\n",filename);
		return false;
	}

	i = 0;

	mprintf(WHITE,"        [");

	if (isbqb) {

		if (m_pBQBFile != NULL)
			delete m_pBQBFile;
		m_pBQBFile = new CBQBFile(*g_pBQBInterface);

		if (!m_pBQBFile->OpenRead( filename )) {
			eprintf("Error: Not a valid BQB file.\n\n");
			return false;
		}

	} else {

		a = fopen( filename, "rb" );
		if (a == NULL) {
			eprintf("Error: Could not open \"%s\" for reading.\n",filename);
			return false;
		}
	}

	fac = 1.0e9 / (LEN_AU2PM * LEN_AU2PM * LEN_AU2PM);

	while (true) {

		if (i != 0) {
			for (z=0;z<(stride-1);z++) {
				if (isbqb) {
					if (!ts.SkipBQB( m_pBQBFile ))
						goto _readdone;
				} else {
					if (!ts.SkipCube( a ))
						goto _readdone;
				}
				mprintf(WHITE,"S");
			}
		}

		if (isbqb) {
			if (!ts.ReadBQB( m_pBQBFile, true ))
				goto _readdone;
		} else {
			if (!ts.ReadCube( a, true ))
				goto _readdone;
		}

		mprintf(WHITE,"#");

		if (m_bWrapCoords) {

			for (z=0;z<g_iGesAtomCount;z++) {

				CxDVector3 coord = g_mBoxToOrtho * ts.m_vaCoords[z];

				n = PeriodicImage1D( coord[0], 0, 1.0 );
				if (n != 0)
					coord[0] -= n;

				n = PeriodicImage1D( coord[1], 0, 1.0 );
				if (n != 0)
					coord[1] -= n;

				n = PeriodicImage1D( coord[2], 0, 1.0 );
				if (n != 0)
					coord[2] -= n;

				ts.m_vaCoords[z] = g_mBoxFromOrtho * coord;
			}
		}

		snap = new CChargeVar2_Snapshot();

		snap->m_vaAtomPositions.resize( g_iGesAtomCount );

		for (z=0;z<g_iGesAtomCount;z++)
			snap->m_vaAtomPositions[z] = Libtravis::Travis::Vector3( ts.m_vaCoords[z][0]/1000.0, ts.m_vaCoords[z][1]/1000.0, ts.m_vaCoords[z][2]/1000.0 );

		snap->m_pCellInfo = new Libtravis::Travis::CellInfo(
			Libtravis::Travis::Vector3( g_mCubeCell(0, 0)/1000.0, g_mCubeCell(0, 1)/1000.0, g_mCubeCell(0, 2)/1000.0 ),
			Libtravis::Travis::Vector3( g_mCubeCell(1, 0)/1000.0, g_mCubeCell(1, 1)/1000.0, g_mCubeCell(1, 2)/1000.0 ),
			Libtravis::Travis::Vector3( g_mCubeCell(2, 0)/1000.0, g_mCubeCell(2, 1)/1000.0, g_mCubeCell(2, 2)/1000.0 )
		);

		snap->m_pElectronDensity = new Libtravis::Travis::DataGrid3D<double>(
			ts.m_pVolumetricData->m_iRes[0],
			ts.m_pVolumetricData->m_iRes[1],
			ts.m_pVolumetricData->m_iRes[2]
		);

		if (i == 0) {
			veca = CxDVector3( g_fCubeXVector[0], g_fCubeXVector[1], g_fCubeXVector[2] );
			vecb = CxDVector3( g_fCubeYVector[0], g_fCubeYVector[1], g_fCubeYVector[2] );
			vecc = CxDVector3( g_fCubeZVector[0], g_fCubeZVector[1], g_fCubeZVector[2] );
			tfac = DotP( veca, CrossP( vecb, vecc ) ); // Bohr^3
			tfac *= 1.0e-9 * LEN_AU2PM * LEN_AU2PM * LEN_AU2PM;
		}

//		for (z=0;z<ts.m_pVolumetricData->m_iRes[0]*ts.m_pVolumetricData->m_iRes[1]*ts.m_pVolumetricData->m_iRes[2];z++)
//			(*snap->m_pElectronDensity)[z] = fac * ts.m_pVolumetricData->m_pBin[z];

		// 3D Transpose
		for (z=0;z<ts.m_pVolumetricData->m_iRes[0];z++)
			for (z2=0;z2<ts.m_pVolumetricData->m_iRes[1];z2++)
				for (z3=0;z3<ts.m_pVolumetricData->m_iRes[2];z3++)
					(*snap->m_pElectronDensity)[z*ts.m_pVolumetricData->m_iRes[1]*ts.m_pVolumetricData->m_iRes[2]+z2*ts.m_pVolumetricData->m_iRes[2]+z3]
						= fac * ts.m_pVolumetricData->m_pBin[z3*ts.m_pVolumetricData->m_iRes[0]*ts.m_pVolumetricData->m_iRes[1]+z2*ts.m_pVolumetricData->m_iRes[0]+z];

		m_oaSnapshots.push_back( snap );

		i++;
		if (i >= steps)
			break;
	}
_readdone:

	if (isbqb) {
		m_pBQBFile->Close();
		delete m_pBQBFile;
		m_pBQBFile = NULL;
	} else
		fclose(a);

	mprintf(WHITE,"] Done.\n");

	mprintf("\n");

	mprintf("    Analyzing first frame...\n");

	tfe = 0;
	for (z=0;z<ts.m_pVolumetricData->m_iRes[0]*ts.m_pVolumetricData->m_iRes[1]*ts.m_pVolumetricData->m_iRes[2];z++)
		tfe += (*(m_oaSnapshots[0]->m_pElectronDensity))[z];

	tfn = 0;
	for (z=0;z<g_iGesAtomCount;z++)
		tfn += ((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge;

	mprintf("      Electron density integral:  %20.10f\n",-tfe*tfac);
	mprintf("      Nuclear charge sum:         %20.10f\n",tfn);
	mprintf("      Total charge:               %20.10f\n",tfn-tfe*tfac);
	
	mprintf("\n");

	m_iaaMaskCellIds.resize( g_oaSingleMolecules.GetSize() );

	for (z=0;z<g_oaSingleMolecules.GetSize();z++) {
		sm = (CSingleMolecule*)g_oaSingleMolecules[z];
		for (z2=0;z2<(int)sm->m_baAtomIndex.GetSize();z2++) {
			if (sm->m_baAtomIndex[z2] == g_iVirtAtomType)
				continue;
			for (z3=0;z3<(int)((CxIntArray*)sm->m_oaAtomOffset[z2])->GetSize();z3++)
				m_iaaMaskCellIds[z].push_back( ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3) );
		}
	}


/*	TrajectoryType traj;
	int i;

	auto& f = traj.addDataSource<FileReaderCube>(filename);
	f.setCoordinatesConversionFactor(1e9 * Constants::bohrRadius);
	f.setAtomPositionsTargetAccessor(traj.getAccessor<AtomPositions>());
	f.setAtomicNumbersTargetAccessor(traj.getAccessor<AtomicNumbers>());
	f.setCellInfoTargetAccessor(traj.getAccessor<CellInfo>());
	f.setDataGridTargetAccessor(traj.getAccessor<ElectronDensity>());

	i = 0;

	mprintf("      Caching...\n");
	mprintf(WHITE,"        [");

	for (auto&& s = traj.begin(); s != traj.end(); ++s) {
		m_oaSnapshotHistory.push_back(*s);
		mprintf(WHITE,"#");
		i++;
		if (i >= steps)
			break;
	}

	mprintf(WHITE,"] Done.\n");

	const TrajectoryType::SnapshotType& s0 = m_oaSnapshotHistory[0];
	RecognizeMolecules molRecog;
	molRecog.setRadiusFactor(1e-3);
	m_oMolSet = molRecog.recognizeMolecules(s0.getData<AtomicNumbers>(), s0.getData<AtomPositions>(), s0.getData<CellInfo>());

	m_iaaMaskCellIds.reserve(m_oMolSet.getNumMolecules());
	for (auto&& m : m_oMolSet.getMolecules())
		m_iaaMaskCellIds.push_back(m.getAtomIds());*/

	return true;
}



std::vector<double> CChargeVar2::calcCharges(const std::vector<double>& radii1, const std::vector<double>& radii2, std::size_t frame, int threadid) {

	FILE *a;
	int z, z2, z3, i, j;
	struct tm *today;
	time_t ltime;
	char buf[256], cbuf[256];
	CSingleMolecule *sm;
	CChargeVar2_Snapshot *snapshot;
	double tfv1v, tfv1i, tfv2v, tfv2i;
	std::chrono::duration<double> duration;


	if ((frame == 0) && (threadid <= 0) && m_bVerboseTiming) {
		if (m_bSecPass)
			mprintf("    Second pass Voronoi integration...\n");
		else
			mprintf("    First pass Voronoi integration...\n");
	}

	const auto before = std::chrono::high_resolution_clock::now();

	snapshot = m_oaSnapshots[frame];

/*	mprintf("@ frame = %d\n",frame);
	mprintf("@ snapshot=%X\n",snapshot);
	mprintf("@ atompos[0] = %f %f %f\n",snapshot->m_vaAtomPositions[0][0],snapshot->m_vaAtomPositions[0][1],snapshot->m_vaAtomPositions[0][2]);
	mprintf("@ cellX = %f %f %f\n",(*(snapshot->m_pCellInfo)).getVectorA()[0],(*(snapshot->m_pCellInfo)).getVectorA()[1],(*(snapshot->m_pCellInfo)).getVectorA()[2]);
	mprintf("@ vol = %f %f %f ...\n",(*(snapshot->m_pElectronDensity))[0],(*(snapshot->m_pElectronDensity))[1],(*(snapshot->m_pElectronDensity))[2]);
*/
	std::vector<double> charges( snapshot->m_vaAtomPositions.size() );

	std::vector<Libtravis::Travis::Vector3> tempatompos;

	tempatompos.assign( snapshot->m_vaAtomPositions.begin(), snapshot->m_vaAtomPositions.end() );

	time(&ltime);
	today = localtime(&ltime);
	strcpy(cbuf,asctime(today));
	cbuf[strlen(cbuf)-1] = 0;
	if (threadid >= 0)
		sprintf(buf,"lastvoronoi_t%02d.txt",threadid+1);
	else
		sprintf(buf,"lastvoronoi.txt");
	a = fopen(buf,"a+t");
	fprintf(a,"%d\nVoronoi Tessellation Pass 1 - Frame %lu, written at %s\n",g_iGesAtomCount,(unsigned long)frame+1,cbuf);
	for (z=0;z<g_iGesAtomCount;z++)
		fprintf(a,"%20.14f  %20.14f  %20.14f  %20.14f\n",tempatompos[z].x(),tempatompos[z].y(),tempatompos[z].z(),radii1[z]);
	fclose(a);

	const auto tiv1a = std::chrono::high_resolution_clock::now();

	Libtravis::Travis::VoronoiTessellation vt = Libtravis::Travis::CalcVoronoi().calcVoronoi( tempatompos, radii1, *(snapshot->m_pCellInfo) );

	const auto tiv1b = std::chrono::high_resolution_clock::now();

	Libtravis::Travis::CalcVoronoiIntegrals::IntegrationConfig config;
	std::vector<Libtravis::Travis::CalcVoronoiIntegrals::IntegrationMask> maskSet;
	config.setRefining( m_iaRefining );
	if (m_bSecPass)
		config.setCalcMaskSet( maskSet, m_iaaMaskCellIds );
	else
		config.setCalcElectricMoments( *snapshot->m_pElectronDensity, true, false, false );

	Libtravis::Travis::CalcVoronoiIntegrals().calcIntegrals( vt, *(snapshot->m_pCellInfo), snapshot->m_pElectronDensity->getSizes(), config );

	duration = std::chrono::high_resolution_clock::now() - tiv1b;
	tfv1i = duration.count();
	duration = tiv1b - tiv1a;
	tfv1v = duration.count();

	//printf("@@ %f %f %f %f %f\n",vt.getCell(0).getIntegrals().charge,vt.getCell(1).getIntegrals().charge,vt.getCell(2).getIntegrals().charge,vt.getCell(3).getIntegrals().charge,vt.getCell(4).getIntegrals().charge);

	if (m_bSecPass) {

		tfv2v = 0;
		tfv2i = 0;

		time(&ltime);
		today = localtime(&ltime);
		strcpy(cbuf,asctime(today));
		cbuf[strlen(cbuf)-1] = 0;
		if (threadid >= 0)
			sprintf(buf,"lastvoronoi_t%02d.txt",threadid+1);
		else
			sprintf(buf,"lastvoronoi.txt");
		a = fopen(buf,"a+t");
		fprintf(a,"%d\nVoronoi Tessellation Pass 2 - Frame %lu, written at %s\n",g_iGesAtomCount,(unsigned long)frame+1,cbuf);
		for (z=0;z<g_iGesAtomCount;z++)
				fprintf(a,"%20.14f  %20.14f  %20.14f  %20.14f\n",tempatompos[z].x(),tempatompos[z].y(),tempatompos[z].z(),radii2[z]);
		fclose(a);

		for (i=0;i<g_oaSingleMolecules.GetSize();i++) {

			sm = (CSingleMolecule*)g_oaSingleMolecules[i];

			std::vector<Libtravis::Travis::Vector3> vt2Positions;
			std::vector<double> vt2Radii;

			for (z2=0;z2<(int)sm->m_baAtomIndex.GetSize();z2++) {
				if (sm->m_baAtomIndex[z2] == g_iVirtAtomType)
					continue;
				for (z3=0;z3<(int)((CxIntArray*)sm->m_oaAtomOffset[z2])->GetSize();z3++) {
					vt2Positions.push_back( tempatompos[ ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3) ] );
					vt2Radii.push_back( radii2[ ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3) ] );
				}
			}

			const auto tiv2a = std::chrono::high_resolution_clock::now();

			Libtravis::Travis::VoronoiTessellation vt2 = Libtravis::Travis::CalcVoronoi().calcVoronoi( vt2Positions, vt2Radii, *snapshot->m_pCellInfo );
	
			const auto tiv2b = std::chrono::high_resolution_clock::now();

			Libtravis::Travis::CalcVoronoiIntegrals::IntegrationConfig config2;
			config2.setRefining( m_iaRefining );
			config2.setIntegrationMask( maskSet.at(i) );
			config2.setCalcElectricMoments( *snapshot->m_pElectronDensity, true, false, false );

			Libtravis::Travis::CalcVoronoiIntegrals().calcIntegrals( vt2, *snapshot->m_pCellInfo, snapshot->m_pElectronDensity->getSizes(), config2 );

			duration = std::chrono::high_resolution_clock::now() - tiv2b;
			tfv2i += duration.count();
			duration = tiv2b - tiv2a;
			tfv2v += duration.count();
	
			j = 0;
			for (z2=0;z2<(int)sm->m_baAtomIndex.GetSize();z2++) {
				if (sm->m_baAtomIndex[z2] == g_iVirtAtomType)
					continue;
				for (z3=0;z3<(int)((CxIntArray*)sm->m_oaAtomOffset[z2])->GetSize();z3++) {
					charges[ ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3) ] = vt2.getCell(j).getIntegrals().charge;
					j++;
				}
			}
		}

		time(&ltime);
		today = localtime(&ltime);
		strcpy(cbuf,asctime(today));
		cbuf[strlen(cbuf)-1] = 0;
		if (threadid >= 0)
			sprintf(buf,"lastcharges_t%02d.txt",threadid+1);
		else
			sprintf(buf,"lastcharges.txt");
		a = fopen(buf,"a+t");
		duration = std::chrono::high_resolution_clock::now() - before;
		fprintf(a,"%d\nVoronoi Tessellation Pass 2 - Frame %lu, Duration %.3f s, Written at %s\n",g_iGesAtomCount,(unsigned long)frame+1,duration.count(),cbuf);
		for (z=0;z<g_iGesAtomCount;z++)
			fprintf(a,"%20.14f  %20.14f\n",vt.getCell(z).getIntegrals().charge,charges[z]);
		fclose(a);

		if ((frame == 0) && (threadid <= 0) && m_bVerboseTiming) {
			mprintf("\n");
			mprintf(WHITE,"    Timing (one frame):\n");
			mprintf("        Molecular Voronoi tessellation: %10.6f s\n",tfv1v);
			mprintf("        Molecular Voronoi integration:  %10.6f s\n",tfv1i);
			mprintf("        Atomic Voronoi tessellation:    %10.6f s  (%d x %9.6f s)\n",tfv2v,g_oaSingleMolecules.GetSize(),tfv2v/g_oaSingleMolecules.GetSize());
			mprintf("        Atomic Voronoi integration:     %10.6f s  (%d x %9.6f s)\n",tfv2i,g_oaSingleMolecules.GetSize(),tfv2i/g_oaSingleMolecules.GetSize());
			mprintf("        Overhead:                       %10.6f s\n",duration.count()-tfv1v-tfv1i-tfv2v-tfv2i);
			mprintf("        Total per frame:                %10.6f s\n",duration.count());
		}

	} else {

		for (j=0;j<g_iGesAtomCount;j++)
			charges.at(j) = vt.getCell(j).getIntegrals().charge;

		time(&ltime);
		today = localtime(&ltime);
		strcpy(cbuf,asctime(today));
		cbuf[strlen(cbuf)-1] = 0;
		if (threadid >= 0)
			sprintf(buf,"lastcharges_t%02d.txt",threadid+1);
		else
			sprintf(buf,"lastcharges.txt");
		a = fopen(buf,"a+t");
		duration = std::chrono::high_resolution_clock::now() - before;
		fprintf(a,"%d\nVoronoi Tessellation Pass 1 - Frame %lu, Duration %.3f s, Written at %s\n",g_iGesAtomCount,(unsigned long)frame+1,duration.count(),cbuf);
		for (z=0;z<g_iGesAtomCount;z++)
			fprintf(a,"%20.14f\n",charges[z]);
		fclose(a);

		if ((frame == 0) && (threadid <= 0) && m_bVerboseTiming) {
			mprintf("\n");
			mprintf(WHITE,"    Timing (one frame):\n");
			mprintf("        Voronoi tessellation: %10.6f s\n",tfv1v);
			mprintf("        Voronoi integration:  %10.6f s\n",tfv1i);
			mprintf("        Overhead:             %10.6f s\n",duration.count()-tfv1v-tfv1i);
			mprintf("        Total per frame:      %10.6f s\n",duration.count());
		}
	}

	return charges;
}



void CChargeVar2::ThreadMain( int id, int nfirst, int nlast ) {

	int zs, z;
	std::vector<double> charges;


	//mprintf("[T%02d] Online, Steps %d ... %d.\n",id+1,nfirst+1,nlast+1);

	for (zs=nfirst;zs<=nlast;zs++) {

		if (g_bAbortAnalysis)
			return;

		if (FileExist("EXIT")) {
			g_bAbortAnalysis = true;
			return;
		}

		charges = calcCharges( m_faRadii1, m_faRadii2, zs, id );

		//mprintf("[T%02d] Finished frame %d.\n",id+1,zs+1);

		for (z=0;z<g_iGesAtomCount;z++) {
			(*(CxDoubleArray*)m_oaAtomCharges[z])[zs] = ((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge - charges[z];
			if (fabs(((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge - charges[z]) > 10.0) {
				eprintf("\nWarning: ");
				mprintf("Extreme charge of %.3f - %.3f = %.3f found on atom %d/%d.  \n",((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge,charges[z],((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge - charges[z],z+1,g_iGesAtomCount);
			}
		}
	}

	//mprintf("[T%02d] Finished, leaving.\n",id+1);
}



void CChargeVar2::Parse() {

	int ti;

	mprintf(WHITE,"\n>>> New \"Double\" Charge Variance Minimization >>>\n\n");

	if (!g_bVolumetricData) {
		eprintf("    Error: This requires volumetric electron density data in each step of the trajectory (.cube or .bqb).\n\n");
		abort();
	}

	m_iLoadFrames = AskUnsignedInteger("    How many frames to use for the charge calculation? [10] ",10);
	m_iStride = AskUnsignedInteger("    Use every n-th frame from the volumetric data trajectory: [1] ",1);

	mprintf("\n");

	m_bUseThreads = AskYesNo("    Use multi-threading (y/n)? [no] ",false);

	if (m_bUseThreads)
		m_iThreads = AskUnsignedInteger("    How many threads to run simultaneously? [8] ",8);

	mprintf("\n");

	ti = AskUnsignedInteger("    Enter refinement factor for interpolation: [1] ",1);
	m_iaRefining[0] = ti;
	m_iaRefining[1] = ti;
	m_iaRefining[2] = ti;

	m_bWrapCoords = AskYesNo("    Wrap atom coordinates before sending them to Voro++ (should not be required, y/n)? [no] ",false);

	mprintf("\n");

	m_bOptimize = AskYesNo("    Run optimization (y), or compute charges from fixed atom radii (n)? [yes] ",true);
	m_bTwoStep = AskYesNo("    Use two-step Voronoi integration for molecular and atomic charges (y/n)? [yes] ",true);

	mprintf("\n");

	if (m_bOptimize) {

		m_bElement = (AskYesNo("    Vary only radii of elements (n) or of all non-equivalent atom types (y)? [yes] ",true) == false);

		mprintf("\n");

		if (m_bTwoStep) {
			m_iInitialRadii1 = AskRangeInteger("    Use unit radii (0), VdW radii (1), covalent radii (2), user input (3), or restart data (4) as initial radii for molecular charges? [1] ", 0, 4, 1 );
			m_iInitialRadii2 = AskRangeInteger("    Use unit radii (0), VdW radii (1), covalent radii (2), user input (3), or restart data (4) as initial radii for atom charges? [2] ", 0, 4, 2 );
		} else {
			m_iInitialRadii1 = AskRangeInteger("    Use unit radii (0), VdW radii (1), covalent radii (2), user input (3), or restart data (4) as initial radii? [1] ", 0, 4, 1 );
			mprintf("\n");
			m_bOneStepTargetAtomic = AskYesNo("    Minimize atomic (y) or molecular (n) charge variance in optimization? [yes] ",true);
		}

		mprintf("\n");

		m_bConjugateGradient = AskYesNo("    Use conjugate gradient (y) or steepest descent / line search method (n)? [yes] ",true);
		
		m_iGradientSamples = AskUnsignedInteger("    How many samples to take per direction for gradient smoothing? [2] ",2);
		m_bGradientOutlier = AskYesNo("    Delete sample that fits worst (\"outlier\") (y/n)? [yes] ",true);
		m_fGradientEps = (double)AskFloat("    Which finite difference to use for numerical gradient calculation (in pm)? [0.05] ",0.05);
		m_fLSEps = (double)AskFloat("    Perform line search until which smallest interval width (in pm)? [0.0001] ",0.0001);
		m_fVarEps = (double)AskFloat("    Which variance change to use as termination condition? [0.0000000001] ",0.0000000001);
		m_fInitialStepLength = (double)AskFloat("    Which initial maximal step length to use in line search (in pm)? [40.0] ",40.0);
		m_bWeightGradient = AskYesNo("    Weight gradient components by R^2 of fit (y/n)? [no] ",false);
		mprintf("\n");
		m_bOutputEveryStep = AskYesNo("    Print radii and charges after every optimization step (y/n)? [no] ",false);

		if (m_bTwoStep) {
			m_bPerformFirstPass = AskYesNo("    Perform first optimization pass (y/n)? [yes] ",true);
			m_bPerformSecondPass = AskYesNo("    Perform second optimization pass (y/n)? [yes] ",true);
		} else {
			m_bPerformFirstPass = true;
			m_bPerformSecondPass = false;
		}

	} else {

		m_bElement = (AskYesNo("    Set only radii of elements (n) or of all non-equivalent atom types (y)? [yes] ",true) == false);

		mprintf("\n");

		if (m_bTwoStep) {
			m_iInitialRadii1 = AskRangeInteger("    Use unit radii (0), VdW radii (1), covalent radii (2), user input (3), or restart data (4) for molecular charges? [1] ", 0, 4, 1 );
			m_iInitialRadii2 = AskRangeInteger("    Use unit radii (0), VdW radii (1), covalent radii (2), user input (3), or restart data (4) for atom charges? [2] ", 0, 4, 2 );
			m_bPerformFirstPass = false;
			m_bPerformSecondPass = true;
		} else {
			m_iInitialRadii1 = AskRangeInteger("    Use unit radii (0), VdW radii (1), covalent radii (2), user input (3), or restart data (4)? [1] ", 0, 4, 1 );
			m_bPerformFirstPass = true;
			m_bPerformSecondPass = false;
		}
	}


	parseCoreCharges();

	mprintf(WHITE,"\n<<< New \"Double\" Charge Variance Minimization <<<\n\n");
}



void CChargeVar2::WriteRadiiRestart() {

	int z;
	FILE *a;
	struct tm *today;
	time_t ltime;
	char cbuf[256];

	a = fopen("voro_radii.restart","wt");
	if (a == NULL) {
		eprintf("CChargeVar2::WriteRadiiRestart(): Error: Could not open \"voro_radii.restart\" for writing.\n");
		return;
	}
	fprintf(a,"%d\n",m_iDimensions);
	time(&ltime);
	today = localtime(&ltime);
	strcpy(cbuf,asctime(today));
	cbuf[strlen(cbuf)-1] = 0;
	fprintf(a,"# Radii restart file, written by TRAVIS at %s\n",cbuf);
	for (z=0;z<m_iDimensions;z++) {
		fprintf(a,"%-20s",(const char*)m_saComponentNames[z]);
		if (m_bTwoStep)
			fprintf(a,"  %16.10f  %16.10f\n",m_vPositionI1[z],m_vPositionI2[z]);
		else
			fprintf(a,"  %16.10f  %16.10f\n",m_vPositionI1[z],0.0);
	}
	fclose(a);
}



bool CChargeVar2::ReadRestartRadii( bool mol ) {

	FILE *a;
	char buf[256], *p, *q;
	int i, z;


	a = fopen( "voro_radii.restart", "rt" );
	if (a == NULL) {
		eprintf("CChargeVar2::ReadRestartRadii(): Error: Could not open \"voro_radii.restart\" for reading.\n" );
		return false;
	}

	(void)!fgets(buf,255,a);
	if (feof(a)) {
		eprintf("CChargeVar2::ReadRestartRadii(): Error: Unexpected end of restart file (1).\n" );
		fclose(a);
		return false;
	}

	i = atoi(buf);

	if (i != m_iDimensions) {
		eprintf("CChargeVar2::ReadRestartRadii(): Error: Dimensionality (%d) does not match restart file (%d).\n", m_iDimensions, i );
		fclose(a);
		return false;
	}

	(void)!fgets(buf,255,a);

	for (z=0;z<i;z++) {

		(void)!fgets(buf,255,a);
		if (feof(a)) {
			eprintf("CChargeVar2::ReadRestartRadii(): Error: Unexpected end of restart file (%d/%d).\n", z+1, i );
			fclose(a);
			return false;
		}

		p = &buf[21];
		q = p;
		while (*q == ' ')
			q++;
		while ((*q != ' ') && (*q != 0))
			q++;
		if (mol) {
			*q = 0;
			m_vPositionI1[z] = atof(p);
		} else {
			p = q;
			while (*q == ' ')
				q++;
			while ((*q != ' ') && (*q != 0))
				q++;
			*q = 0;
			m_vPositionI2[z] = atof(p);
		}
	}

	fclose(a);

	return true;
}



void CChargeVar2::PerformAnalysis() {

	int z, z2, z3, z4, z5, z6, z0, i, cyc, maxmollen, maxradlen;
	double sw, osw, oosw, f1, f2=0, beta;
	double tfmean, tfvar, tfmin, tfmax, tfsum;
	char buf[256], *p;
	CxDoubleArray *fa;
	CxDVectorN vtemp, grad, lastgrad, sdir, normsdir;
	CMolecule *m;
	CSingleMolecule *sm, *sm2;
	CMolAtom *ma, *ma0;
	FILE *a=NULL;
	bool downhill, forcedownhill;


	mprintf(WHITE,"\n>>> New \"Double\" Charge Variance Minimization starting >>>\n\n");

	// Rename existing files as backup
	if (m_bUseThreads) {
		for (z=0;z<m_iThreads;z++) {
			sprintf(buf,"lastvoronoi_t%02d.txt",z+1);
			a = OpenFileWrite(buf,true);
			fclose(a);
			sprintf(buf,"lastcharges_t%02d.txt",z+1);
			a = OpenFileWrite(buf,true);
			fclose(a);
		}
	} else {
		a = OpenFileWrite("lastvoronoi.txt",true);
		fclose(a);
		a = OpenFileWrite("lastcharges.txt",true);
		fclose(a);
	}

	g_bVoronoiMoments = false;
	g_bVoroSilent = true;
	downhill = false;
	forcedownhill = false;

	g_oaTimeSteps.SetSize(g_iStepHistory);
	for (z=0;z<g_iStepHistory;z++)
		g_oaTimeSteps[z] = NULL;

	mprintf("    Trying to cache %d frames with a stride of %d from input trajectory \"%s\"\n",m_iLoadFrames,m_iStride,g_sInputTraj);
	mprintf("      (resolution %d x %d x %d, requires %s of memory).\n",g_TimeStep.m_pVolumetricData->m_iRes[0],g_TimeStep.m_pVolumetricData->m_iRes[1],g_TimeStep.m_pVolumetricData->m_iRes[2],FormatBytes(((double)m_iLoadFrames)*g_TimeStep.m_pVolumetricData->m_iRes[0]*g_TimeStep.m_pVolumetricData->m_iRes[1]*g_TimeStep.m_pVolumetricData->m_iRes[2]*8));
	loadData(g_sInputTraj,m_iLoadFrames,m_iStride);
	mprintf("    %lu/%d frames cached.\n\n",(unsigned long)m_oaSnapshots.size(),m_iLoadFrames);

	maxradlen = 0;
	m_iaAtomVectorAssignment.SetSize(g_iGesAtomCount);
	if (m_bElement) {
		m_iDimensions = g_oaAtoms.GetSize() - 1;
		m_saComponentNames.SetSize(m_iDimensions);
		for (z=0;z<m_iDimensions;z++) {
			p = new char[strlen(((CAtom*)g_oaAtoms[z])->m_sName)+1];
			m_saComponentNames[z] = p;
			strcpy((char*)m_saComponentNames[z],((CAtom*)g_oaAtoms[z])->m_sName);
			if ((int)strlen(((CAtom*)g_oaAtoms[z])->m_sName) > maxradlen)
				maxradlen = (int)strlen(((CAtom*)g_oaAtoms[z])->m_sName);
		}
		for (z=0;z<g_iGesAtomCount;z++)
			m_iaAtomVectorAssignment[z] = g_waAtomRealElement[z];
	} else {
		for (z=0;z<g_oaMolecules.GetSize();z++) {
			m = (CMolecule*)g_oaMolecules[z];
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
		//	mprintf("Molecule %s ...\n",m->m_sName);
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++) {
				if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
					continue;
				ma = NULL;
				z4 = -1;
				//mprintf("@@@@ Element %s ...\n",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName);
				for (z3=0;z3<sm->m_oaMolAtomsSorted.GetSize();z3++) {
					//mprintf("@MolAtom %d/%d: Element=%2s, Number=%d\n",z3+1,sm->m_oaMolAtomsSorted.GetSize(),(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[((CMolAtom*)sm->m_oaMolAtomsSorted[z3])->m_iType]])->m_sName,((CMolAtom*)sm->m_oaMolAtomsSorted[z3])->m_iNumber+1);
					if (((CMolAtom*)sm->m_oaMolAtomsSorted[z3])->m_iType != z2)
						continue;
					ma0 = ma;
					ma = (CMolAtom*)sm->m_oaMolAtomsSorted[z3];
					if (ma0 == NULL) {
						z4 = ma->m_iNumber;
						continue;
					}
					if (ma->m_liAtomCode != ma0->m_liAtomCode) {
						if (z4 == ma0->m_iNumber) {
				//			mprintf("  [%s] %s%d\n",m->m_sName,((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z4+1);
							sprintf(buf,"[%s] %s%d",m->m_sName,(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z4+1);
						} else {
				//			mprintf("  [%s] %s%d-%d\n",m->m_sName,((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z4+1,ma0->m_iNumber+1);
							sprintf(buf,"[%s] %s%d-%d",m->m_sName,(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z4+1,ma0->m_iNumber+1);
						}
						for (z5=z4;z5<=ma0->m_iNumber;z5++) {
							for (z6=0;z6<m->m_laSingleMolIndex.GetSize();z6++) {
								sm2 = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z6]];
								i = ((CxIntArray*)sm2->m_oaAtomOffset[z2])->GetAt(z5);
								m_iaAtomVectorAssignment[i] = m_saComponentNames.GetSize();
							}
						}
						p = new char[strlen(buf)+1];
						if ((int)strlen(buf) > maxradlen)
							maxradlen = (int)strlen(buf);
						strcpy(p,buf);
						m_saComponentNames.Add(p);
						m_iaComponentElements.Add(m->m_baAtomIndex[z2]);
						z4 = ma->m_iNumber;
					}
				}
				if (ma != NULL) {
					if (z4 == ma->m_iNumber) {
				//		mprintf("* [%s] %s%d\n",m->m_sName,((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z4+1);
						sprintf(buf,"[%s] %s%d",m->m_sName,(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z4+1);
					} else {
				//		mprintf("* [%s] %s%d-%d\n",m->m_sName,((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z4+1,ma->m_iNumber+1);
						sprintf(buf,"[%s] %s%d-%d",m->m_sName,(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z4+1,ma->m_iNumber+1);
					}
					for (z5=z4;z5<=ma->m_iNumber;z5++) {
						for (z6=0;z6<m->m_laSingleMolIndex.GetSize();z6++) {
							sm2 = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z6]];
							i = ((CxIntArray*)sm2->m_oaAtomOffset[z2])->GetAt(z5);
							m_iaAtomVectorAssignment[i] = m_saComponentNames.GetSize();
						}
					}
					p = new char[strlen(buf)+1];
					if ((int)strlen(buf) > maxradlen)
						maxradlen = (int)strlen(buf);
					strcpy(p,buf);
					m_saComponentNames.Add(p);
					m_iaComponentElements.Add(m->m_baAtomIndex[z2]);
				}
			}
		}
		m_iDimensions = m_saComponentNames.GetSize();
/*		mprintf("\n*** Assignment Check:\n");
		for (z=0;z<g_iGesAtomCount;z++)
			mprintf("  %3d: %s, vector %d\n",z,((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,m_iaAtomVectorAssignment[z]);
*/	}

	mprintf("    This is a %d-dimensional optimization problem.\n\n",m_iDimensions);

	m_vPosition = CxDVectorN(m_iDimensions);
	m_vPositionI1 = CxDVectorN(m_iDimensions);
	m_vPositionI2 = CxDVectorN(m_iDimensions);
	grad = CxDVectorN(m_iDimensions);

	m_faRadii1.resize(g_iGesAtomCount);
	m_faRadii2.resize(g_iGesAtomCount);

	switch(m_iInitialRadii1) {

		case 0: // Unity radii
			for (z=0;z<m_iDimensions;z++)
				m_vPositionI1[z] = 100.0;
			break;

		case 1: // VdW radii
			for (z=0;z<m_iDimensions;z++)
				m_vPositionI1[z] = ((CAtom*)g_oaAtoms[m_iaComponentElements[z]])->m_pElement->m_fVdWRadius;
			break;

		case 2: // Covalent radii
			for (z=0;z<m_iDimensions;z++)
				m_vPositionI1[z] = ((CAtom*)g_oaAtoms[m_iaComponentElements[z]])->m_pElement->m_fRadius;
			break;

		case 3: // Custom radii
			for (z=0;z<m_iDimensions;z++)
				m_vPositionI1[z] = AskFloat_ND("    Please enter initial \"molecular\" atom radius for %s (in pm): ",(const char*)m_saComponentNames[z]);
			break;

		case 4: // Restart data
			if (!ReadRestartRadii( true )) {
				eprintf("Error reading restart radii.\n");
				abort();
			}
			break;

		default:
			eprintf("CChargeVar2::PerformAnalysis(): Internal error.\n");
			abort();
	}

	mprintf("\n");

	if (m_bTwoStep) {
		switch(m_iInitialRadii2) {

			case 0: // Unity radii
				for (z=0;z<m_iDimensions;z++)
					m_vPositionI2[z] = 100.0;
				break;

			case 1: // VdW radii
				for (z=0;z<m_iDimensions;z++)
					m_vPositionI2[z] = ((CAtom*)g_oaAtoms[m_iaComponentElements[z]])->m_pElement->m_fVdWRadius;
				break;

			case 2: // Covalent radii
				for (z=0;z<m_iDimensions;z++)
					m_vPositionI2[z] = ((CAtom*)g_oaAtoms[m_iaComponentElements[z]])->m_pElement->m_fRadius;
				break;

			case 3: // Custom radii
				for (z=0;z<m_iDimensions;z++)
					m_vPositionI2[z] = AskFloat_ND("    Please enter initial \"atomic\" atom radius for %s (in pm): ",(const char*)m_saComponentNames[z]);
				break;

			case 4: // Restart data
				if (!ReadRestartRadii( false )) {
					eprintf("Error reading restart radii.\n");
					abort();
				}
				break;

			default:
				eprintf("CChargeVar2::PerformAnalysis(): Internal error.\n");
				abort();
		}
	}


	/************************************************************************************/

	if (!m_bPerformFirstPass)
		m_bSecPass = true;
	else
		m_bSecPass = false;

_secpass:

	if (m_bOptimize) {
		if (m_bTwoStep) {
			if (m_bSecPass) {
				mprintf("\n");
				mprintf(GREEN,"***********************************************************\n");
				mprintf(GREEN,"****    Pass 2: Atomic Charge Variance Minimization    ****\n");
				mprintf(GREEN,"***********************************************************\n");
				mprintf("\n");
			} else {
				mprintf("\n");
				mprintf(GREEN,"**************************************************************\n");
				mprintf(GREEN,"****    Pass 1: Molecular Charge Variance Minimization    ****\n");
				mprintf(GREEN,"**************************************************************\n");
				mprintf("\n");
			}
		} else { // Not m_bTwoStep
			mprintf("\n");
			mprintf(GREEN,"**************************************************************\n");
			if (m_bOneStepTargetAtomic)
				mprintf(GREEN,"****     One-step Atomic Charge Variance Minimization     ****\n");
			else
				mprintf(GREEN,"****    One-step Molecular Charge Variance Minimization   ****\n");
			mprintf(GREEN,"**************************************************************\n");
			mprintf("\n");
		}
	} else { // Not optimize
		if (m_bTwoStep) {
			mprintf("\n");
			mprintf(GREEN,"***********************************************************\n");
			mprintf(GREEN,"****        Two-Step Atomic Charge Calculation         ****\n");
			mprintf(GREEN,"***********************************************************\n");
			mprintf("\n");
		} else {
			mprintf("\n");
			mprintf(GREEN,"***********************************************************\n");
			mprintf(GREEN,"****        One-Step Atomic Charge Calculation         ****\n");
			mprintf(GREEN,"***********************************************************\n");
			mprintf("\n");
		}
	}

	m_iGradCounter = 0;
	m_iLSCounter = 0;

	mprintf(WHITE,"    Initial radii:\n\n");
	for (z=0;z<m_iDimensions;z++) {
		mprintf("      r(%s",(const char*)m_saComponentNames[z]);
		for (z2=0;z2<maxradlen-(int)strlen((char*)m_saComponentNames[z]);z2++)
			mprintf(" ");
		if (m_bTwoStep)
			mprintf(") = %10.6f pm  /  %10.6f pm\n",m_vPositionI1[z],m_vPositionI2[z]);
		else
			mprintf(") = %10.6f pm\n",m_vPositionI1[z]);
	}
	mprintf("\n");

	WriteRadiiRestart();

	if (m_bSecPass) {
		for (z=0;z<m_iDimensions;z++)
			m_vPosition[z] = m_vPositionI2[z];
	} else {
		for (z=0;z<m_iDimensions;z++)
			m_vPosition[z] = m_vPositionI1[z];
	}

	m_oaAtomCharges.SetSize(g_iGesAtomCount);
	for (z=0;z<g_iGesAtomCount;z++) {
		fa = new CxDoubleArray();
		fa->SetSize((int)m_oaSnapshots.size());
		m_oaAtomCharges[z] = fa;
	}

	m_faMolChargeMax.SetSize(g_oaMolecules.GetSize());
	m_faMolChargeMin.SetSize(g_oaMolecules.GetSize());
	m_faMolChargeMean.SetSize(g_oaMolecules.GetSize());
	m_faMolChargeVar.SetSize(g_oaMolecules.GetSize());

	i = 0;
	for (z=0;z<g_oaMolecules.GetSize();z++) {
		m = (CMolecule*)g_oaMolecules[z];
		for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++) {
			if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
				continue;
			i += m->m_waAtomCount[z2];
		}
	}
	m_faAtomChargeMin.SetSize(i);
	m_faAtomChargeMax.SetSize(i);
	m_faAtomChargeMean.SetSize(i);
	m_faAtomChargeVar.SetSize(i);

	m_bVerboseTiming = true;

	f1 = CalculateObjective( &m_vPosition );

	m_bVerboseTiming = false;

	// Initial max. step length
	osw = m_fInitialStepLength;

	if (m_bSecPass)
		a = fopen("varopt_atomic.csv","wt");
	else
		a = fopen("varopt_molecular.csv","wt");

	fprintf(a,"Step; Variance; ");

	for (z=0;z<m_iDimensions;z++) {
		fprintf(a,"r(%s)",(char*)m_saComponentNames[z]);
		if (z < g_oaAtoms.GetSize()-2)
			fprintf(a,"; ");
	}
	fprintf(a,"\n");

	fprintf(a,"0; %.12f; ",f1);

	for (z=0;z<m_iDimensions;z++) {
		fprintf(a,"%.10f",m_vPosition[z]);
		if (z < m_iDimensions-1)
			fprintf(a,"; ");
	}
	fprintf(a,"\n");

	fflush(a);

	mprintf("\n    Initial variance:    %12.9f\n",f1);

	if (!m_bOptimize)
		goto _ende;

	cyc = 0;
	while (true) {

		cyc++;

		mprintf(YELLOW,"\n*** Optimization cycle %d ***\n\n",cyc);
		mprintf("    Computing gradient: ");

		lastgrad = grad;

		ComputeGradient(&m_vPosition,&grad);

		if (g_bAbortAnalysis)
			goto _ende;

		mprintf("\n    Gradient norm: %12.10f\n",grad.GetLength());

		grad *= -1.0; // Downhill direction

_newdir:
		if (m_bConjugateGradient) {

			if ((cyc == 1) || forcedownhill) { // First step: Steepest descent

				sdir = grad;
				mprintf("    Conjugate Gradient: Doing steepest descent step.\n");
				forcedownhill = false;
				downhill = true;
			} else {
				beta = DotP(grad, grad-lastgrad) / DotP(lastgrad, lastgrad); // Polak-Ribiere Formula
				if (beta < 0) {
					mprintf("    Conjugate Gradient: Beta = %8.5f, resetting direction (downhill step).\n",beta);
					beta = 0;
					downhill = true;
				} else {
					mprintf("    Conjugate Gradient: Beta = %8.5f.\n",beta);
					downhill = false;
				}
				sdir = grad + beta * sdir;
			}
		} else
			sdir = grad; // Just go downhill

		normsdir = sdir;
		normsdir.Normalize();

		mprintf("    Step direction: ");
		normsdir.Dump();
		mprintf("\n\n");

		oosw = osw;
_again:
		mprintf("    Performing line search (maximal step length = %.5f)...\n",osw);
		LineSearch(&m_vPosition,&normsdir,osw,&sw);

		if (g_bAbortAnalysis)
			goto _ende;

		mprintf("\n    Obtained %.5f as optimal step length.\n",sw);

		vtemp = m_vPosition + normsdir * sw;

		f2 = CalculateObjective(&vtemp);

		mprintf("    New variance:    %12.9f\n",f2);
		mprintf("    Variance change: %12.9f\n",f2-f1);

		if (g_bAbortAnalysis)
			goto _ende;

		if (f2 > f1) { // Reject step

			mprintf(RED,"    Error:");
			mprintf(" Variance increased, step rejected.\n");
			if (osw <= 0.01) {
				if (m_bConjugateGradient && !downhill) {
					mprintf(YELLOW,"\n    Cannot reduce step length any further. Trying downhill step.\n");
					osw = oosw;
					forcedownhill = true;
					goto _newdir;
				} else {
					mprintf(YELLOW,"\n*** Max. step length %.5f is already <= 0.01. Aborting optimization.\n",osw);
					break;
				}
			}
			osw /= 10.0;
			mprintf(YELLOW,"    Trying %.5f as new maximal step length.\n\n",osw);
			goto _again;

		} else { // Accept step

			m_vPosition = vtemp;
			if (sw < 0.2)
				osw = 1.0;
			else if (sw < 10.0)
				osw = 5.0*fabs(sw);
			else if (sw < 50.0)
				osw = 2.0*fabs(sw);
			else if (sw < 100.0)
				osw = sw;
			else osw = 100.0;
		}

		for (z=0;z<m_iDimensions;z++) {
			if (m_vPosition[z] < 0) {
				mprintf(YELLOW,"    Note:");
				mprintf("Swapped sign of negative radius for %s.\n",(const char*)m_saComponentNames[z]);
				m_vPosition[z] = -m_vPosition[z];	
			}
		}

		fprintf(a,"%d;  %.12f;  ",cyc,f2);

		for (z=0;z<m_iDimensions;z++) {
			fprintf(a,"%.10f",m_vPosition[z]);
			if (z < m_iDimensions-1)
				fprintf(a,";  ");
		}
		fprintf(a,"\n");
		fflush(a);

		mprintf("    New position:  ");
		m_vPosition.Dump();
		mprintf("\n");

		if (m_bSecPass) {
			for (z=0;z<m_iDimensions;z++)
				m_vPositionI2[z] = m_vPosition[z];
		} else {
			for (z=0;z<m_iDimensions;z++)
				m_vPositionI1[z] = m_vPosition[z];
		}

		WriteRadiiRestart();

		if (f1-f2 < m_fVarEps) {
			if (m_bConjugateGradient && !downhill) {
				mprintf(YELLOW,"\n    Variance decrease below EPS. Trying downhill step.\n");
				forcedownhill = true;
			} else {
				mprintf(YELLOW,"\n*** Optimization converged.\n");
				break;
			}
		}

		if (m_bOutputEveryStep) {
			mprintf(WHITE,"\n    *** Current radii:\n");
			for (z=0;z<m_iDimensions;z++) {
				mprintf("      r(%s",(const char*)m_saComponentNames[z]);
				for (z2=0;z2<maxradlen-(int)strlen((char*)m_saComponentNames[z]);z2++)
					mprintf(" ");
				mprintf(") = %10.6f pm\n",m_vPosition[z]);
			}
			mprintf("\n");

			mprintf(WHITE,"    *** Current molecular charges:\n");
			maxmollen = 0;
			for (z=0;z<g_oaMolecules.GetSize();z++)
				if ((int)strlen(((CMolecule*)g_oaMolecules[z])->m_sName) > maxmollen)
					maxmollen = (int)strlen(((CMolecule*)g_oaMolecules[z])->m_sName);
			for (z=0;z<g_oaMolecules.GetSize();z++) {
				mprintf("      %s:",((CMolecule*)g_oaMolecules[z])->m_sName);
				for (z2=0;z2<maxmollen-(int)strlen(((CMolecule*)g_oaMolecules[z])->m_sName);z2++)
					mprintf(" ");
				mprintf("   Mean = %12.8f,  Std.Dev. = %12.8f,  Var = %12.8f,  Min = %12.8f,  Max = %12.8f\n",m_faMolChargeMean[z],sqrt(m_faMolChargeVar[z]),m_faMolChargeVar[z],m_faMolChargeMin[z],m_faMolChargeMax[z]);
			}
			mprintf("\n");

			mprintf(WHITE,"    *** Current atomic partial charges (per individual atom):\n");
			if (maxmollen < (int)strlen("Molecule"))
				maxmollen = (int)strlen("Molecule");
			mprintf(YELLOW,"      Molecule / Atom  ");
			for (z4=0;z4<maxmollen-(int)strlen("Molecule");z4++)
				mprintf(" ");
			mprintf(YELLOW," Mean         Std.Dev.     Var          Min          Max\n");
			i = 0;
			for (z=0;z<g_oaMolecules.GetSize();z++) {
				m = (CMolecule*)g_oaMolecules[z];
				for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++) {
					if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
						continue;
					for (z3=0;z3<m->m_waAtomCount[z2];z3++) {
						mprintf("      [%s] ",m->m_sName);
						for (z4=0;z4<maxmollen-(int)strlen(m->m_sName);z4++)
							mprintf(" ");
						mprintf("%2s%-2d  ",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1);
						mprintf("%9.6f    %9.6f    %9.6f    %9.6f    %9.6f\n",m_faAtomChargeMean[i],sqrt(m_faAtomChargeVar[i]),m_faAtomChargeVar[i],m_faAtomChargeMin[i],m_faAtomChargeMax[i]);
						i++;
					}
				}
			}
			mprintf("\n");

			mprintf(WHITE,"    *** Current atomic partial charges (per equivalent atom type):\n");
			if (maxradlen < (int)strlen("Molecule"))
				maxradlen = (int)strlen("Molecule");
			mprintf(YELLOW,"      Molecule / Atom  ");
			for (z4=0;z4<maxradlen-(int)strlen("Molecule");z4++)
				mprintf(" ");
			mprintf(YELLOW," Mean         Std.Dev.     Var          Min          Max\n");
			for (z0=0;z0<m_iDimensions;z0++) {
				i = 0;
				tfmean = 0;
				tfvar = 0;
				tfmin = 1.0e30;
				tfmax = -1.0e30;
				tfsum = 0;
				for (z=0;z<g_oaMolecules.GetSize();z++) {
					m = (CMolecule*)g_oaMolecules[z];
					sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
					for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++) {
						if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
							continue;
						for (z3=0;z3<m->m_waAtomCount[z2];z3++) {
							if (m_iaAtomVectorAssignment[ ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3) ] == z0) {
								tfmean += m_faAtomChargeMean[i];
								tfvar += m_faAtomChargeVar[i];
								if (m_faAtomChargeMin[i] < tfmin)
									tfmin = m_faAtomChargeMin[i];
								if (m_faAtomChargeMax[i] > tfmax)
									tfmax = m_faAtomChargeMax[i];
								tfsum++;
							}
							i++;
						}
					}
				}
				mprintf("      %s ",(const char*)m_saComponentNames[z0]);
				for (z4=0;z4<maxradlen-(int)strlen((const char*)m_saComponentNames[z0]);z4++)
					mprintf(" ");
				mprintf("        %9.6f    %9.6f    %9.6f    %9.6f    %9.6f\n",
					tfmean/tfsum,
					sqrt(tfvar/tfsum),
					tfvar/tfsum,
					tfmin,
					tfmax
				);
			}
			mprintf("\n");

			mprintf(WHITE,"    *** Current overall charge variance:");
			mprintf("  %12.9f\n",f2);
		}

		f1 = f2;
	}

_ende:
	fclose(a);

	mprintf("\n");
	if (m_bOptimize) {
		mprintf(WHITE,"\n*** Final radii:\n\n");
		for (z=0;z<m_iDimensions;z++) {
			mprintf("      r(%s",(const char*)m_saComponentNames[z]);
			for (z2=0;z2<maxradlen-(int)strlen((char*)m_saComponentNames[z]);z2++)
				mprintf(" ");
			mprintf(") = %10.6f pm\n",m_vPosition[z]);
		}
		mprintf("\n");
		f2 = CalculateObjective(&m_vPosition);
	}

	mprintf(WHITE,"\n*** Final molecular charges:\n\n");
	maxmollen = 0;
	for (z=0;z<g_oaMolecules.GetSize();z++)
		if ((int)strlen(((CMolecule*)g_oaMolecules[z])->m_sName) > maxmollen)
			maxmollen = (int)strlen(((CMolecule*)g_oaMolecules[z])->m_sName);
	for (z=0;z<g_oaMolecules.GetSize();z++) {
		mprintf("      %s:",((CMolecule*)g_oaMolecules[z])->m_sName);
		for (z2=0;z2<maxmollen-(int)strlen(((CMolecule*)g_oaMolecules[z])->m_sName);z2++)
			mprintf(" ");
		mprintf("   Mean = %12.8f,  Std.Dev. = %12.8f,  Var = %12.8f,  Min = %12.8f,  Max = %12.8f\n",m_faMolChargeMean[z],sqrt(m_faMolChargeVar[z]),m_faMolChargeVar[z],m_faMolChargeMin[z],m_faMolChargeMax[z]);
	}
	mprintf("\n");

	mprintf(WHITE,"\n*** Final atomic partial charges (per individual atom):\n\n");
	if (maxmollen < (int)strlen("Molecule"))
		maxmollen = (int)strlen("Molecule");
	mprintf(YELLOW,"    Molecule / Atom  ");
	for (z4=0;z4<maxmollen-(int)strlen("Molecule");z4++)
		mprintf(" ");
	mprintf(YELLOW," Mean         Std.Dev.     Var          Min          Max\n");
	i = 0;
	for (z=0;z<g_oaMolecules.GetSize();z++) {
		m = (CMolecule*)g_oaMolecules[z];
		for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++) {
			if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
				continue;
			for (z3=0;z3<m->m_waAtomCount[z2];z3++) {
				mprintf("    [%s] ",m->m_sName);
				for (z4=0;z4<maxmollen-(int)strlen(m->m_sName);z4++)
					mprintf(" ");
				mprintf("%2s%-2d  ",(const char*)((CAtom*)g_oaAtoms[m->m_baAtomIndex[z2]])->m_sName,z3+1);
				mprintf("%9.6f    %9.6f    %9.6f    %9.6f    %9.6f\n",
					m_faAtomChargeMean[i],
					sqrt(m_faAtomChargeVar[i]),
					m_faAtomChargeVar[i],
					m_faAtomChargeMin[i],
					m_faAtomChargeMax[i]
				);
				i++;
			}
		}
	}
	mprintf("\n");

	mprintf(WHITE,"\n*** Final atomic partial charges (per equivalent atom type):\n\n");
	if (maxradlen < (int)strlen("Molecule"))
		maxradlen = (int)strlen("Molecule");
	mprintf(YELLOW,"    Molecule / Atom  ");
	for (z4=0;z4<maxradlen-(int)strlen("Molecule");z4++)
		mprintf(" ");
	mprintf(YELLOW," Mean         Std.Dev.     Var          Min          Max\n");
	for (z0=0;z0<m_iDimensions;z0++) {
		i = 0;
		tfmean = 0;
		tfvar = 0;
		tfmin = 1.0e30;
		tfmax = -1.0e30;
		tfsum = 0;
		for (z=0;z<g_oaMolecules.GetSize();z++) {
			m = (CMolecule*)g_oaMolecules[z];
			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[0]];
			for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++) {
				if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
					continue;
				for (z3=0;z3<m->m_waAtomCount[z2];z3++) {
					if (m_iaAtomVectorAssignment[ ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3) ] == z0) {
						tfmean += m_faAtomChargeMean[i];
						tfvar += m_faAtomChargeVar[i];
						if (m_faAtomChargeMin[i] < tfmin)
							tfmin = m_faAtomChargeMin[i];
						if (m_faAtomChargeMax[i] > tfmax)
							tfmax = m_faAtomChargeMax[i];
						tfsum++;
					}
					i++;
				}
			}
		}
		mprintf("    %s ",(const char*)m_saComponentNames[z0]);
		for (z4=0;z4<maxradlen-(int)strlen((const char*)m_saComponentNames[z0]);z4++)
			mprintf(" ");
		mprintf("        %9.6f    %9.6f    %9.6f    %9.6f    %9.6f\n",
			tfmean/tfsum,
			sqrt(tfvar/tfsum),
			tfvar/tfsum,
			tfmin,
			tfmax
		);
	}
	mprintf("\n\n");

	if (m_bSecPass) {

		if (m_bOptimize) {
			mprintf(WHITE,"*** Final overall atomic charge variance:");
			mprintf("  %12.9f\n",f2);
		}

		//for (z=0;z<m_iDimensions;z++)
		//	m_vPositionI2[z] = m_vPosition[z];

		if (m_bOptimize) {
			mprintf("\n");
			mprintf(GREEN,"*******************************\n");
			mprintf(GREEN,"****    Pass 2 Finished    ****\n");
			mprintf(GREEN,"*******************************\n");
			mprintf("\n");
		} else {
			mprintf("\n");
			mprintf(GREEN,"************************************\n");
			mprintf(GREEN,"****    Calculation Finished    ****\n");
			mprintf(GREEN,"************************************\n");
			mprintf("\n");
		}

		if (m_bOptimize) {
			mprintf(WHITE,"\n*** Final radii:\n\n");
			for (z=0;z<m_iDimensions;z++) {
				mprintf("      r(%s",(const char*)m_saComponentNames[z]);
				for (z2=0;z2<maxradlen-(int)strlen((char*)m_saComponentNames[z]);z2++)
					mprintf(" ");
				mprintf(") = %10.6f pm  /  %10.6f pm\n",m_vPositionI1[z],m_vPositionI2[z]);
			}
		}

	} else {

		if (m_bOptimize) {
			mprintf(WHITE,"*** Final overall molecular charge variance:");
			mprintf("  %12.9f\n",f2);
		}

		//for (z=0;z<m_iDimensions;z++)
		//	m_vPositionI1[z] = m_vPosition[z];

		if (m_bTwoStep && m_bPerformSecondPass && !g_bAbortAnalysis) {
			m_bSecPass = true;
			goto _secpass;
		}
	}

	mprintf(WHITE,"\n<<< New \"Double\" Charge Variance Minimization done <<<\n\n");
}



void CChargeVar2::ComputeGradient(CxDVectorN *pos, CxDVectorN *grad) {

	double *px, *py, *rsq, sxx, syy, sxy, mx, my, t0, t1, fa, fb, tf;
	int z, z2, i, skip, n;
	CxDVectorN temp;
	CxString buf;
	FILE *a;
	double **filebuf;
	unsigned long ti0;


	ti0 = (unsigned long)time(NULL);
	n = 0;

	px = new double[2*m_iGradientSamples+1];
	py = new double[2*m_iGradientSamples+1];
	rsq = new double[m_iDimensions];

	filebuf = new double*[m_iDimensions];
	for (z=0;z<m_iDimensions;z++)
		filebuf[z] = new double[2*m_iGradientSamples+1];

	mprintf(WHITE,"[");
	t0 = CalculateObjective(pos);
	n++;
	mprintf(WHITE,"#");
	for (z=0;z<m_iDimensions;z++) {
		skip = -1;
		for (z2=-m_iGradientSamples;z2<=m_iGradientSamples;z2++) {
			px[z2+m_iGradientSamples] = z2*m_fGradientEps;
			if (z2 == 0) {
				py[z2+m_iGradientSamples] = 0;
				filebuf[z][z2+m_iGradientSamples] = 0;
				continue;
			}
			temp = *pos;
			temp[z] += z2*m_fGradientEps;
			t1 = CalculateObjective(&temp);
			if (g_bAbortAnalysis)
				return;
			n++;
			py[z2+m_iGradientSamples] = t1 - t0;
			filebuf[z][z2+m_iGradientSamples] = t1 - t0;
			mprintf(WHITE,"#");
		}

_skipped:
		mx = 0;
		my = 0;
		i = 0;
		for (z2=0;z2<2*m_iGradientSamples+1;z2++) {
			if (z2 == skip)
				continue;
			mx += px[z2];
			my += py[z2];
			i++;
		}
		my /= i;
		mx /= i;

		sxx = 0;
		syy = 0;
		sxy = 0;

		for (z2=0;z2<2*m_iGradientSamples+1;z2++) {
			if (z2 == skip)
				continue;
			sxx += px[z2] * px[z2];
			sxy += px[z2] * (py[z2]-my);
			syy += (py[z2]-my) * (py[z2]-my);
		}

		fa = my;
		fb = sxy / sxx;

		if ((skip == -1) && m_bGradientOutlier) {
			tf = 0;
			for (z2=0;z2<2*m_iGradientSamples+1;z2++) {
				if (fabs(py[z2] - (fa + fb * px[z2])) > tf) {
					tf = fabs(py[z2] - (fa + fb * px[z2]));
					skip = z2;
				}
			}
			filebuf[z][skip] = 1.0e30;
			goto _skipped;
		}

		rsq[z] = pow2(sxy / sqrt(sxx * syy));
		if (m_bWeightGradient)
			(*grad)[z] = sxy / sxx * (0.1 + 0.9*rsq[z]);
		else
			(*grad)[z] = sxy / sxx;
	}
	delete[] px;
	delete[] py;

	mprintf(WHITE,"]\n");
	mprintf("      (R^2:  ");
	for (z=0;z<m_iDimensions;z++) {
		mprintf("%5.3f",rsq[z]);
		if (z < m_iDimensions-1)
			mprintf("  ");
	}
	mprintf(", %d samples, time %lus)",n,time(NULL)-ti0);
	delete[] rsq;

	if (m_bSecPass)
		buf.sprintf("gradient_atomic_%04d.csv",m_iGradCounter++);
	else
		buf.sprintf("gradient_molecular_%04d.csv",m_iGradCounter++);
	a = fopen(buf,"wt");

	if (a != NULL) {
		for (z=0;z<2*m_iGradientSamples+1;z++) {
			fprintf(a,"%.8f;  ",(z-m_iGradientSamples)*m_fGradientEps);
			for (z2=0;z2<m_iDimensions;z2++) {
				if (filebuf[z2][z] < 1.0e20)
					fprintf(a,"%.14f",filebuf[z2][z]);
				else
					fprintf(a,"0");

				if (z2 < m_iDimensions-1)
					fprintf(a,";  ");
			}
			fprintf(a,"\n");
		}
		fclose(a);
	}
}



double phi2 = (1 + sqrt(5)) / 2;
double resphi2 = 2 - phi2;



bool CChargeVar2::LineSearch(CxDVectorN *pos, CxDVectorN *dir, double maxlen, double *result) {

	double a, c, fa, fc;
	CxDVectorN temp;
	int z;
	CxString buf;
	FILE *f;

	m_faLSX.RemoveAll_KeepSize();
	m_faLSY.RemoveAll_KeepSize();

	a = -maxlen/5.0;
	temp = *pos + a * (*dir);
	fa = CalculateObjective(&temp);
	m_faLSX.Add(a);
	m_faLSY.Add(fa);

	c = maxlen;
	temp = *pos + c * (*dir);
	fc = CalculateObjective(&temp);
	m_faLSX.Add(c);
	m_faLSY.Add(fc);

	*result = REC_GoldenSectionSearch(a,c,c,fa,fc,fc,pos,dir,0);

	if (m_bSecPass)
		buf.sprintf("linesearch_atomic_%04d.csv",m_iLSCounter++);
	else
		buf.sprintf("linesearch_molecular_%04d.csv",m_iLSCounter++);
	f = fopen(buf,"wt");

	if (f != NULL) {
		for (z=0;z<m_faLSX.GetSize();z++)
			fprintf(f,"%.12f;  %.12f\n",m_faLSX[z],m_faLSY[z]);
		fclose(f);
	}

	return true;
}



double CChargeVar2::REC_GoldenSectionSearch(double a, double b, double c, double fa, double fb, double fc, CxDVectorN *pos, CxDVectorN *dir, int depth) {

	double x, fx;
	CxDVectorN temp;

	if ((c - b) > (b - a))
		x = b + resphi2 * (c - b);
	else
		x = b - resphi2 * (b - a);

	if (fabs(c - a) < m_fLSEps) {
		mprintf("      CChargeVar2::REC_GoldenSectionSearch(): EPS reached, returning.\n");
		return (c + a) / 2.0; 
	}

	temp = *pos + x * (*dir);
	fx = CalculateObjective( &temp );

	if (g_bAbortAnalysis)
		return (c + a) / 2.0; 

	m_faLSX.Add(x);
	m_faLSY.Add(fx);

	mprintf("      GoldenSection %2d,  Interval %8.5f .. %8.5f,  Len %8.5f,  Var %11.9f\n",depth,a,c,x,fx);

	if (fx == fb) {
		mprintf("      CChargeVar2::REC_GoldenSectionSearch(): fx == fb, returning.\n");
		return (c + a) / 2.0; 
	}

	if (fx < fb) {
		if ((c - b) > (b - a))
			return REC_GoldenSectionSearch(b, x, c, fb, fx, fc, pos, dir, depth+1);
		else
			return REC_GoldenSectionSearch(a, x, b, fa, fx, fb, pos, dir, depth+1);
	} else {
		if ((c - b) > (b - a))
			return REC_GoldenSectionSearch(a, b, x, fa, fb, fx, pos, dir, depth+1);
		else
			return REC_GoldenSectionSearch(x, b, c, fx, fb, fc, pos, dir, depth+1);
	}
}



double CChargeVar2::CalculateObjective( CxDVectorN *vec ) {

	double t;
	int z;

	if (m_bSecPass) {
		for (z=0;z<g_iGesAtomCount;z++) {
			m_faRadii1[z] = m_vPositionI1[m_iaAtomVectorAssignment[z]] / 1000.0;
			m_faRadii2[z] = (*vec)[m_iaAtomVectorAssignment[z]] / 1000.0;
		}
	} else {
		for (z=0;z<g_iGesAtomCount;z++) {
			m_faRadii1[z] = (*vec)[m_iaAtomVectorAssignment[z]] / 1000.0;
			m_faRadii2[z] = m_vPositionI2[m_iaAtomVectorAssignment[z]] / 1000.0;
		}
	}

	CalculateCharges();
	CalculateMinMaxVar();

	t = 0;

	if (m_bTwoStep) {
		if (m_bSecPass) {
			for (z=0;z<m_faAtomChargeVar.GetSize();z++)
				t += m_faAtomChargeVar[z] / ((double)m_faAtomChargeVar.GetSize());
		} else {
			for (z=0;z<g_oaMolecules.GetSize();z++)
				t += m_faMolChargeVar[z] / ((double)g_oaMolecules.GetSize());
		}
	} else {
		if (m_bOneStepTargetAtomic) {
			for (z=0;z<m_faAtomChargeVar.GetSize();z++)
				t += m_faAtomChargeVar[z] / ((double)m_faAtomChargeVar.GetSize());
		} else {
			for (z=0;z<g_oaMolecules.GetSize();z++)
				t += m_faMolChargeVar[z] / ((double)g_oaMolecules.GetSize());
		}
	}

	return t;
}



void CChargeVar2::CalculateMinMaxVar() {

	int z, z2, z3, z4, z5, i, k, o;
	double c;
	CMolecule *m;
	CSingleMolecule *sm;
	CxDoubleArray *fa;

	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];

		m_faMolChargeMean[z] = 0;
		m_faMolChargeMin[z] = 1.0e30;
		m_faMolChargeMax[z] = -1.0e30;
		m_faMolChargeVar[z] = 0;

		i = 0;
		for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++) {

			sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];

			for (z5=0;z5<(int)m_oaSnapshots.size();z5++) {

				c = 0;
				for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++) {

					if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
						continue;

					for (z3=0;z3<m->m_waAtomCount[z2];z3++) {
						o = ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3);
						fa = (CxDoubleArray*)m_oaAtomCharges[o];
						c += (*fa)[z5];
					}
				}

				if (c > m_faMolChargeMax[z])
					m_faMolChargeMax[z] = c;
				if (c < m_faMolChargeMin[z])
					m_faMolChargeMin[z] = c;

				m_faMolChargeMean[z] += c;
				m_faMolChargeVar[z] += c * c;
			}
			i += (int)m_oaSnapshots.size();
		}

		m_faMolChargeMean[z] /= i;
		m_faMolChargeVar[z] /= i;
		m_faMolChargeVar[z] -= m_faMolChargeMean[z] * m_faMolChargeMean[z];
	}

	k = 0;
	for (z=0;z<g_oaMolecules.GetSize();z++) {

		m = (CMolecule*)g_oaMolecules[z];
		for (z2=0;z2<m->m_baAtomIndex.GetSize();z2++) {

			if (m->m_baAtomIndex[z2] == g_iVirtAtomType)
				continue;

			for (z3=0;z3<m->m_waAtomCount[z2];z3++) {
				m_faAtomChargeMin[k] = 1.0e30;
				m_faAtomChargeMax[k] = -1.0e30;
				m_faAtomChargeMean[k] = 0;
				m_faAtomChargeVar[k] = 0;
				i = 0;
				for (z4=0;z4<m->m_laSingleMolIndex.GetSize();z4++) {
					sm = (CSingleMolecule*)g_oaSingleMolecules[m->m_laSingleMolIndex[z4]];
					o = ((CxIntArray*)sm->m_oaAtomOffset[z2])->GetAt(z3);
					fa = (CxDoubleArray*)m_oaAtomCharges[o];
					for (z5=0;z5<fa->GetSize();z5++) {
						if ((*fa)[z5] > m_faAtomChargeMax[k])
							m_faAtomChargeMax[k] = (*fa)[z5];
						if ((*fa)[z5] < m_faAtomChargeMin[k])
							m_faAtomChargeMin[k] = (*fa)[z5];
						m_faAtomChargeMean[k] += (*fa)[z5];
						m_faAtomChargeVar[k] += (*fa)[z5] * (*fa)[z5];
					}
					i += fa->GetSize();
				}
				m_faAtomChargeMean[k] /= i;
				m_faAtomChargeVar[k] /= i;
				m_faAtomChargeVar[k] -= m_faAtomChargeMean[k] * m_faAtomChargeMean[k];
				k++;
			}
		}
	}
}




bool CChargeVar2::CalculateCharges() {

	int zs, z, n;
	std::vector<double> charges;
	double tf=0;


	if (m_bUseThreads) {

		if (g_bAbortAnalysis) {
			mprintf("\nAnalysis aborted by user.\n");
			return false;
		}
	
		if (FileExist("EXIT")) {
			mprintf("\n\n*** File \"EXIT\" detected. Aborting analysis. ***\n\n");
			remove("EXIT");
			g_bAbortAnalysis = true;
			return false;
		}
	
		m_oaThreads.resize( m_iThreads-1 );

		n = (int)ceil( (double)m_oaSnapshots.size() / m_iThreads );

		for (z=0;z<m_iThreads-1;z++)
			m_oaThreads[z] = new std::thread( &CChargeVar2::ThreadMain, this, z+1, (z+1)*n, MIN( ((z+2)*n)-1, (int)m_oaSnapshots.size()-1 ) );

		// Thread 0 runs in the main thread
		ThreadMain( 0, 0, MIN( n-1, (int)m_oaSnapshots.size()-1 ) );

		for (z=0;z<m_iThreads-1;z++) {
			m_oaThreads[z]->join();
			delete m_oaThreads[z];
			m_oaThreads[z] = NULL;
		}

	} else {

		for (zs=0;zs<(int)m_oaSnapshots.size();zs++) {
	
			if (g_bAbortAnalysis) {
				mprintf("\nAnalysis aborted by user.\n");
				return false;
			}
	
			if (FileExist("EXIT")) {
				mprintf("\n\n*** File \"EXIT\" detected. Aborting analysis. ***\n\n");
				remove("EXIT");
				g_bAbortAnalysis = true;
				return false;
			}
	
			charges = calcCharges( m_faRadii1, m_faRadii2, zs, -1 );
	
			for (z=0;z<g_iGesAtomCount;z++) {
				(*(CxDoubleArray*)m_oaAtomCharges[z])[zs] = ((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge - charges[z];
				if (fabs(((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge - charges[z]) > 10.0) {
					eprintf("\nWarning: ");
					mprintf("Extreme charge of %.3f - %.3f = %.3f found on atom %d/%d.  \n",((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge,charges[z],((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge - charges[z],z+1,g_iGesAtomCount);
				}
				tf += charges[z];
//				if (z < 10)
//					mprintf("@ %s: %f - %f = %f\n",(const char*)((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_sName,((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge,charges[z],((CAtom*)g_oaAtoms[g_waAtomRealElement[z]])->m_fCharge - charges[z]);
			}
		}
	}

//	mprintf("@ Sum of charges: %.10E\n",tf);

	return true;
}



#endif





