/*****************************************************************************

    TRAVIS - Trajectory Analyzer and Visualizer

    http://www.travis-analyzer.de/

    Copyright (c) 2009-2021 Martin Brehm
                  2012-2021 Martin Thomas
                  2016-2021 Sascha Gehrke

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
#include "vorowrapper.h"
#include "globalvar.h"
#include "reactive.h"

#ifdef NEW_CHARGEVAR
	#include "chargevar2.h"
#endif




const char *GetRevisionInfo_globalvar(unsigned int len) {
	static char buf[256];
	GET_REVISION_INFO( buf, len );
	return buf;
}


const char *GetSourceVersion_globalvar() {
	static char buf[256];
	GET_SOURCE_VERSION( buf );
	return buf;
}


CxObArray g_oaAtoms;  // Die Atomsorten der Simulation
CxObArray g_oaMolecules; // Die Molekuelsorten der Simulation
CxObArray g_oaSingleMolecules;
CxObArray g_oaObserv;
CxObArray g_oaVirtualAtoms;
CTimeStep g_TimeStep;
char g_iFixMol; // Welches Molekuel fixieren wir?
int g_iFixAtomType[3]; // Welches Element fixieren wir? (3 mal wegen den 3 Fixpunkten)
int g_iFixRealAtomType[3]; // Wie Eben, nur Index auf g_pAtoms
int g_iFixAtom[3]; // Welche Atome in in dem Molekuel fixieren?
int g_iCurrentTimeStep, g_iLastTimeStep, g_iNextTimeStep;
unsigned long g_iSteps; // Zaehler fuer die schon verarbeiteten Zeitschritte
unsigned char g_iBinning;
bool g_bFold;
bool g_bWarnUnsteady;
double g_fUnsteadyLimit = 10000000.0; // Unit is pm/ps
FILE *g_fSaveJustTraj;
unsigned char g_iVirtAtomType;
double g_fVelLimit, g_fForceLimit;
double g_fMSVDFLevel, g_fMSFDFLevel;
bool g_bCalcVel, g_bCalcForces;
bool g_bManSVDFLevel, g_bManSFDFLevel;
bool g_bManSMeanVDFLevel, g_bManSMeanFDFLevel;
bool g_bManSMaxVDFLevel, g_bManSMaxFDFLevel;
double g_fSVDFLevel[16];
double g_fSMeanVDFLevel[16];
double g_fSMaxVDFLevel[16];
double g_fSFDFLevel[16];
double g_fSMeanFDFLevel[16];
double g_fSMaxFDFLevel[16];
//char g_sRefEnv[256];
CxString g_sRefEnv;
bool g_bUseVelocities;
bool g_bUseForces;
double g_fBoxX, g_fBoxY, g_fBoxZ;
CxDVec3Array g_pRefMol;  // Die Koordinaten des Referenzmolekueles (zum akkumulieren)
bool g_bAvg; // Gemitteltes Referenzmolekuel ausgeben
bool g_bMiddleAvg; // Soll das Referenzmolekuel wirklich gemittelt werden (true), oder soll einfach das erstbeste Molekuel genommen werden (false)?
int g_iGesAtomCount; // Gesamtzahl der Atome pro Zeitschritt - also in der Simulation
int g_iGesVirtAtomCount;
int g_iMaxStep; // Bis zu welchem Zeitschritt geht die Analyse? -1 fuer alle
unsigned char g_iSwapAtoms; // Sollen Atome, die sich frei drehen oder sich vertauschen, zurueckgetauscht werden?
unsigned char g_iNormSDF; // 0 - jede SDF einzeln auf 100, 1 - gar nicht normieren
unsigned char g_iSVDFLevel, g_iSFDFLevel;
//char g_pAtomMassLabels[][3] = { "H", "He", "Li", "B", "C", "N", "O", "F", "Na", "Al", "Si", "P", "S", "Cl", "Ru", "Pt", "Br", "Fe", "U" };
//unsigned char g_iAtomMassCount = 19;
//double g_fAtomMass[] = { 1.f, 4.f, 7.f, 11.f, 12.f, 14.f, 16.f, 19.f, 23.f, 27.f, 28.f, 31.f, 32.f, 35.45f, 101.07f, 195.1f, 79.9f, 55.8f, 238.0f };
//double g_fAtomRadius[] = { 0.37f, 0.32f, 1.34f, 0.9f, 0.82f, 0.77f, 0.75f, 0.73f, 0.71f, 1.54f, 1.18f, 1.11f, 1.06f, 1.02f, 0.99f, 1.3f, 1.f, 1.52f, 1.0f };
//double g_fAtomVDWRadius[] = { 1.20f, 1.22f, 0.f, 2.08f, 1.85f, 1.54f, 1.40f, 1.35f, 2.31f, 2.05f, 2.f, 1.9f, 1.85f, 1.81f, 1.46f, 10.f, 2.0f, 1.52f, 1.42f };
bool g_bMegaMat;
bool g_bMatOnlyBind;
bool g_bVFHisto;
unsigned short g_iHistogramRes = 256;
FILE *g_fRefTrajec, *g_fRefEnv, *g_fVRDF[32];
CTimeStep *g_pTempTimestep;
bool g_bFoldAtomwise;
double g_fBondFactor;
bool g_bSaveJustTraj;
//bool g_bRefEnvVirt;
bool g_bSaveVirtAtoms;
FILE *g_fVFCorr[64];
unsigned short g_iVFCorrCount;
//bool g_bDynamicNeighbor;
bool g_bSaveRefWithEnv;
bool g_bRefEnvCenter;
//bool g_bScanNeighbors;
bool g_bPeriodicX, g_bPeriodicY, g_bPeriodicZ;
bool g_bPeriodic;
int g_iStride;
CAsciiArt g_oAsciiArt;
CTimeStep *g_pT2Timestep;
//CNbSearch *g_pNbAll;
int g_iClusterPos;
CxByteArray g_baAtomIndex;
//bool g_bRefEnvAtomwise;
bool g_bScanMolecules;
CACF *g_pGlobalVACF;
CACF *g_pGlobalDipACF;
int g_iRefSystemDim;
//bool g_bSaveTrajVirt;
FILE *g_fPos, *g_fVel, *g_fForce;
double g_fTimestepLength;
CFFT *g_pFFT;
long g_iHighestStepNumber;
int g_iDotCounter;
bool g_bStepSkipped;
bool g_bAbortAnalysis;
//int g_iMSDDepth;
CxIntArray g_laBondBlackList;
CxIntArray g_laBondAddList;
bool g_bMolRecNoDomain;
int g_iBondBlackListUsed;
int g_iFirstStepSkipped;
bool g_bGlobalPsycho;
int g_iWannierAtomType;
double g_fWannierCharge;

bool g_bCombinedDF;
bool g_bSDF, g_bRDF; 
bool g_bVDF, g_bFDF, g_bADF;
bool g_bDipole;

CxObArray g_oaAnalysisGroups;
CxObArray g_oaAnalyses;
int g_iStepHistory;
CxObArray g_oaTimeSteps;
double g_fMaxVel, g_fMaxForce;
double g_fLMaxVel, g_fLMaxForce;
double g_fLMidVel, g_fLMidForce;
bool g_bAsciiArt;
//char g_sInputVel[256];
//char g_sInputForce[256];
//char g_sInputCtrl[256];
CxString g_sInputVel;
CxString g_sInputForce;
CxString g_sInputCtrl;
CDatabase *g_pDatabase;

CxIntArray g_laAtomSMIndex; // Enthaelt fuer jedes Atom, in welchem SingleMolecule es sich befindet, oder -1 wenn gar nicht
CxIntArray g_laAtomSMLocalIndex;
//CxWordArray g_waAtomMolElem; // Enthaelt fuer jedes Atom den Index im m_oaMolAtoms-Array des SingleMolecules
CxWordArray g_waAtomElement;
CxWordArray g_waAtomMolNumber;
CxWordArray g_waAtomMolIndex;
CxWordArray g_waAtomMolUID;
CxWordArray g_waAtomRealElement;
//CxDoubleArray g_faAtomCode;
std::vector<LargeInteger> g_liAtomCode;
CxDoubleArray g_faVdWRadius;
CxObArray g_oaMolAtoms;
std::vector<std::vector<int> > g_iaaMoleculeTypeSplit;

CxByteArray g_baAtomPassedCondition;

int g_iCDFChannels;
bool g_bWannier;

int g_iScanMolStep;
char *g_sInputTraj;
bool g_bCDF;
int *g_iObsChannel; // 1 - RDF, 2 - ADF, 3 - DDF, 4 - DipDF, 5 - VDF, 6 - FDF
bool g_bDDF;
bool g_bRevSDF;
bool g_bMSD;
bool g_bCutCluster;
bool g_bSaveRefEnv;
unsigned char g_iSaveRefMol;
bool g_bRefEnvFix;
//int g_iNbStride;
//bool g_bKeepNbCount;
int g_iClusterSteps, g_iClusterCount;
CxIntArray g_iaClusterSteps, g_iaClusterMol;
bool g_bNbAnalysis;
bool g_bVHDF;
//CxObArray g_oaNbSearches;
CNbSet *g_pNbSet;
bool g_bVACF;
bool g_bVFDF;
bool g_bGlobalVACF;
bool g_bGlobalDipACF;
int g_iSDFSmoothGrade;
bool g_bACFWindowFunction;
//int g_iMSDStride;
//int g_iVHDFDepth;
//int g_iVHDFStride;
//int g_iVHDFMinDepth;
bool g_bSaveVelForce;
double g_fVelPercentage, g_fForcePercentage;
//bool g_bRefNoVirt;
unsigned long g_iSaveGesAtoms;
bool g_bUnwrap;
CxObArray g_oaSaveMolecules;
//bool g_bTrajAtomwise;
CxDVec3Array g_vaUnwrapArray;
int g_iBeginStep;
bool g_bSkipDoubleSteps;
//bool g_bUseMassCenters;

int g_iNumberPos; // Die wievielte Nummer in der XYZ-Kommentarzeile ist die Schrittzahl?
FILE *g_fInput;

bool g_bDCDFirstFrame;
CTimeStep *g_pDCDLabelsTS;
bool g_bDCDSkipHeader;

bool g_bDipACF;
//int g_iMaxACFDepth;
bool g_bACF;
bool g_bDipDF;

bool g_bInputRedirected;

//int g_iSDFScale; // 0 = ppm, 1 = pm^-3, 2 = nm^-3, 3 = Rel. to Uniform Density
bool g_bSDFUniform;

bool g_bDDisp; 
bool g_bDLDF; 
bool g_bDLDisp; 


bool g_bDACF; 
bool g_bAggregation;
//CAggregation *g_pAggregation;

bool g_bKeepUnfoldedCoords;

CxObArray g_oaElements;

bool g_bTDO;
bool g_bTDOEqui;
int g_iTDOCount;
int g_iTDOStride;
int g_iTDOStart;
CxIntArray g_laTDOSteps;
double g_fTDOBleaching;


bool g_bMultiInterval;
int g_iMultiIntervalBegin;
int g_iMultiIntervalStride;
int g_iMultiIntervalLength;
CxIntArray g_laMultiIntervalStart;
CxIntArray g_laMultiIntervalEnd;

bool g_bDoubleBox;
int g_iDoubleBoxX;
int g_iDoubleBoxY;
int g_iDoubleBoxZ;
int g_iDoubleBoxFactor;

int g_iTrajSteps;

bool g_bRDyn;

bool g_bBondACF;
int g_iBondACFDepth;
bool g_bBondACFDebug;
bool g_bBondACFNormalize;
bool g_bBondACFSymmetrize;
bool g_bBondACFWindow;

bool g_bACFFFT = true;

bool g_bMSDCacheMode;
bool g_bRDynCacheMode;
bool g_bVACFCacheMode;

char *g_sInputFile;
FILE *g_fInputFile;

char *g_sControlFile;
bool g_bControlRun;
CDatabase *g_pControlDB;
int g_iControlRT;
bool g_bReactive;
CReactiveEngine *g_pReactEngine;

bool g_bSaveCondSnapshot;
bool g_bSaveCondWholeBox;
FILE *g_fSaveCondFile;
int g_iSaveCondCount;

bool g_bCond;
//CConditionGroup *g_pCondition;

bool g_bCombined;

unsigned long g_iFastForwardPos;

int g_iNbhMode;

bool g_bWriteAtomwise;


time_t g_iStartTime, g_iEndTime;

int g_iTrajFormat;

bool g_bLAMMPSColumnsAnalyzed = false;
bool g_bLAMMPSWrappedCoords = false;
bool g_bLAMMPSVelocity = false;
int g_iLAMMPSColumnCount;
char *g_sLAMMPSColumnString;
std::vector<int> g_iaLAMMPSOffsetToColumn;
std::vector<int> g_iaLAMMPSColumnToOffset;

bool g_bScanVelocities;

bool g_bSaveTrajPattern;
CxString g_sSaveTrajPattern;
int g_iSaveTrajPatternPos;

int g_iScanNbhStart;
int g_iScanNbhSteps;
int g_iScanNbhStride;

int g_iScanVelStart;
int g_iScanVelSteps;
int g_iScanVelStride;


bool g_bSilentProgress;

char *g_sHomeDir;
char *g_sSettingsFile;
char *g_sHostName;
char *g_sWorkingDir;

bool g_bSMode;

bool g_bCreateRevSDF;
bool g_bNoColor;

bool g_bAmberCellAsked;
bool g_bAmberCellInfo;

int g_iColorIntensity;
bool g_bNPT;
//char g_sNPTFile[256];
CxString g_sNPTFile;
FILE *g_fNPTFile;

bool g_bNbExchange;

bool g_bVerbose;
bool *g_pUniteTemp;

bool g_bCenterZero;

bool g_bSaveJustCenter;
int g_iSaveJustMol;
int g_iSaveJustSM;
int g_iSaveJustAtomType;
int g_iSaveJustRealAtomType;
int g_iSaveJustAtom;

bool g_bTimeDiff;
bool g_bDeriv;
int g_iDerivLast;
int g_iDerivCurr;
int g_iDerivNext;

int g_iCloseAtomCounter;

bool g_bRemoveCOM;

bool g_bVoro;
int g_iVoroMemory;
CVoroWrapper *g_pVoroWrapper;

bool g_bSaxonize;

bool g_bUnknownElements;


bool g_bNeedMoleculeWrap;


bool g_bAdvanced1;
bool g_bAdvanced2;

bool g_bVoid;


bool g_bDipoleDefined;

bool g_bDipolGrimme;

bool g_bVoidSDF;

bool g_bRamanFromPolarizability;


bool g_bSaveCoordsUnchanged;

bool g_bDens;

bool g_bSaveTrajNoRot;

bool g_bCheckWrite;

bool g_bRaman;
bool g_bIRSpec;
bool g_bPowerSpec;

bool g_bKeepOriginalCoords;

bool g_bShowConf;
bool g_bWriteConf;
char *g_sConfFile;

bool g_bPDF;

int g_iStrideDetect;

double g_fMinPeriodic;


bool g_bAggrTopo;
CAggrTopoEngine *g_pAggrTopoEngine;

bool g_bSankey;
char g_sSankeyFile[512];

bool g_bRegionAnalysis;
CxIntArray g_iaSMRegion;

bool g_bDumpDipoleVector;
bool g_bDumpDipoleAbs;
bool g_bDumpDipoleXYZ;
CxObArray g_oaDumpDipoleVector; // Contains CxIntAtrrays
FILE *g_fDumpDipole;
FILE *g_fDumpDipoleXYZ;
FILE *g_fDumpDipolePURE;
int g_iDumpDipoleSMCount;
int g_iDumpDipoleXYZAtoms;
double g_fDumpDipoleScale;

bool g_bPlDF;
bool g_bLiDF;

bool g_bStreamInput;

bool g_bGlobalIR;
CReorDyn *g_pGlobalIR;

bool g_bLMFitSilent;

int g_iFitDegree;
double *g_pExpSpecExpo;
bool g_bLMFitSmooth;
int g_iLMMaxIter;

bool g_bUnwrapWannier;
bool g_bDipoleRefFixed;

double *g_fLSpecEvolveBuf;

bool g_bXYZ4thCol;
bool g_bXYZComment3Numbers;
bool g_bXYZComment6Numbers;
bool g_bXYZComment9Numbers;
bool g_bReadChargesFrom4thXYZ;

bool g_bProcSplit;
int g_iProcSplitLength;

bool g_bPDFFixAtom;
int g_iPDFFixAtom;

int g_iRemoveCOMFixAtom;

bool g_bPairMSD;

bool g_bSFac;
//CStructureFactor *g_pSFac;

bool g_bEnvWriteDetailedInfo;
bool g_bEnvSortNb;
bool g_bEnvDisableSortNb;

bool g_bShowCredits;

bool g_bWriteInputOrder;

bool g_bPlProj;

bool g_bNormalCoordinate = false;
bool g_bChiral = false;
bool g_bSortWannier = false;
bool g_bVCD = false;
bool g_bEckartTransform = false;
bool g_bPower = false;
bool g_bIR = false;

CTetraPak *g_pTetraPak;

bool g_bTegri;
bool g_bDomA;
CDomainEngine *g_pDomainEngine;

CxDoubleArray g_faVoronoiRadii;
bool g_bBoxNonOrtho;
bool g_bFoundNonOrtho;
double g_fBoxAngleA;
double g_fBoxAngleB;
double g_fBoxAngleC;
CxDMatrix3 g_mBoxToOrtho;
CxDMatrix3 g_mBoxFromOrtho;
double g_fBoxMinDiamA;
double g_fBoxMinDiamB;
double g_fBoxMinDiamC;
double g_fBoxMinDiam;
double g_fBoxVolume;
bool g_bWriteOrtho;
double g_fWriteOrthoFac;

bool g_bSDFVoro;

bool g_bSDFMap;
CxObArray g_oaSDFMaps;
int g_iSDFMapSmoothGrade;

bool g_bCHDF = false;

int g_iCubeRes[3] = { -1, -1, -1 };
double g_fCubeXStep = 0.0;
double g_fCubeYStep = 0.0;
double g_fCubeZStep = 0.0;
double g_fCubeXVector[3] = { 0.0, 0.0, 0.0 };
double g_fCubeYVector[3] = { 0.0, 0.0, 0.0 };
double g_fCubeZVector[3] = { 0.0, 0.0, 0.0 };
CxDMatrix3 g_mCubeCell = CxDMatrix3(0.0);
int g_iCubeXStride = 1;
int g_iCubeYStride = 1;
int g_iCubeZStride = 1;
int g_iCubeXMismatch = 0;
int g_iCubeYMismatch = 0;
int g_iCubeZMismatch = 0;
bool g_bCubeTimeDev = false;
FILE *g_fPDESolverInfoFile = NULL;
double g_fBackgroundDensity = 1.0e-6;
double g_fPDEConvThresh = 0.05;
int g_iPDEMaxIter = 30;

bool g_bVoroSilent;
int g_iVoroPrintLevel;

bool g_bMinimizeChargeVar2;
#ifdef NEW_CHARGEVAR
CChargeVar2 *g_pChargeVar2;
#endif


bool g_bVoronoiMoments;
bool g_bVoroIntEquitable;
bool g_bVoroIntegrateCharge = false;
bool g_bVoroIntegrateDipoleMoment = false;
bool g_bVoroIntegrateQuadrupoleMoment = false;
bool g_bQuadrupoleKeepTrace = false;
bool g_bVoroIntegrateTotalCurrent = false;
bool g_bVoroIntegrateMagneticMoment = false;


//char g_sAmberParmFile[1024];
CxString g_sAmberParmFile;



bool g_bContactMatrix;
CContactMatrix *g_pContactMatrix;

bool g_bCubeStream = false;
CxMemFile *g_fCubeMemFile = NULL;
int g_iCubeMemFileLines = 0;
int g_iCubeMemFileSteps = 0;

bool g_bDipoleRestart = false;
bool g_bLoadDipoleRestart = false;
FILE *g_fDipoleRestartFile = NULL;
bool g_bMagneticDipoleRestart = false;
bool g_bLoadMagneticDipoleRestart = false;
FILE *g_fMagneticDipoleRestartFile = NULL;

double g_fIntegratorTimestep = 0;

bool g_bSetUpPolarizabilityCalc = false;

bool g_bPolarizabilityDefined = false;
int g_iPolarizabilityMode = 0;
int g_iPolarizabilityConf[3] = { 0, 0, 0 };
FILE *g_fPolarizabilityFile[6] = { NULL, NULL, NULL, NULL, NULL, NULL };
double g_fPolarizabilityFieldStrength = 0.0;

bool g_bLAMMPSCharge;
bool g_bReadLAMMPSCharges;
bool g_bLAMMPSForce;

bool g_bHBond = false;
bool g_bIonPair = false;
bool g_bHole = false;

int g_iRefMolNum;

bool g_bProcAlternativeLabels;
CxObArray g_oaProcAlternativeLabels;
bool g_bProcDumpAwkScript;
bool g_bProcWriteComments;
bool g_bProcCellComment;
bool g_bProcCellCommentAngles;
bool g_bInterWarning;

CBQBFile *g_pBQBFile = NULL;
//CBQBEngine *g_pBQBEngine = NULL;

bool g_bFlipCellVectorX;
bool g_bFlipCellVectorY;
bool g_bFlipCellVectorZ;

bool g_bOrder;
COrderEngine *g_pOrderEngine;


bool g_bTDDF;
CTDDFEngine *g_pTDDFEngine;

bool g_bGeoDens;
CGeoDensEngine *g_pGeoDensEngine;

bool g_bROA;
CROAEngine *g_pROAEngine;

FILE *g_fMolIntegralFile;

bool g_bVolumetricData;
bool g_bElMagProperties;

bool g_bReadVelocity;

bool g_bVACFNew;

double g_fVelocityConversion;

bool g_bStrideParsed;

bool g_bFixedPlProj;
CFixedPlProj *g_pFixedPlProj;

bool g_bProcAddMesh;
bool g_bProcAddMeshZ;
double g_fProcAddMeshZPos;
bool g_bProcAddMeshY;
double g_fProcAddMeshYPos;
bool g_bProcAddMeshX;
double g_fProcAddMeshXPos;
int g_iProcAddMeshGrid;
CxString g_sProcAddMeshLabel;
bool g_bProcAddMeshJitter;
double g_fProcAddMeshJitter;

CBQBInterface *g_pBQBInterface;
bool g_bUseBQB;

bool g_bKuehneLong;

int g_iCellVectorFileColumns;
std::vector<std::string> g_saDCDRemarkLines;










